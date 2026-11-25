// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file composition.hpp
 * @brief SISO composition: cascade and switched controller / experiment paths.
 *
 * - Cascade — outer.control(r, y_outer) → inner.control(outer_u, y_inner)
 * - SwitchedController — Normal | Experiment | Backup ownership of u
 *
 * SwitchedController stays on the SisoController concept so it nests and drops
 * into simulate_sampled / multi-rate harnesses. Cascade takes a distinct
 * measurement for each loop (control(r, y_outer, y_inner)) and is not a
 * SisoController — one shared y would regulate the inner loop to the outer
 * measurement. Experiments (relay, FRF drivers)
 * that are not SisoControllers are not stored inside SwitchedController;
 * call SwitchedController::set_mode and feed the experiment yourself, or
 * wrap a thin SisoController adapter that holds the experiment output.
 *
 * @see concepts.hpp SisoController
 */

#include "damp/backend.hpp"
#include "damp/concepts.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix_traits.hpp" // eps

namespace damp {

/**
 * @brief Series cascade of two SISO controllers: outer → inner reference.
 *
 * Tick:
 *
 *     u_outer = outer.control(r, y_outer)
 *     u       = inner.control(u_outer, y_inner)
 *
 * While the inner law's command is clamped, at its own u_lim or by a clamp reported with
 * back_calculate, the inner loop cannot follow its reference, and the outer law is told its
 * command was held at the last one the inner loop could follow, the outer command of the last
 * tick that went out uncut: back_calculate of it. Its memory holds as against a clamp of its own
 * (an outer PI's integral gives back what each tick added past it), where it ran on for as long
 * as the inner loop could not follow and overshot when it could. An inner law without a u_lim
 * tells only through back_calculate.
 *
 * Typical use: outer position / voltage loop, inner velocity / current.
 * Alias CascadePPI for PI+PI when both loops are PIController.
 *
 * @see "Advanced PID Control" (Åström & Hägglund, 2006), §3.5: anti-windup by conditional
 *      integration and by tracking
 *
 * @tparam Outer Outer controller type (SisoController)
 * @tparam Inner Inner controller type (SisoController)
 * @tparam T     Scalar type
 */
template<typename Outer, typename Inner, typename T = float>
    requires SisoController<Outer, T> && SisoController<Inner, T>
class Cascade {
public:
    constexpr Cascade() = default;

    constexpr Cascade(const Outer& outer, const Inner& inner)
        : outer_(outer), inner_(inner) {}

    /**
     * @brief One cascade tick.
     * @param r       Outer reference
     * @param y_outer Outer measurement (e.g. position, voltage)
     * @param y_inner Inner measurement (e.g. velocity, current)
     */
    [[nodiscard]] constexpr T control(T r, T y_outer, T y_inner) {
        if (!cut_) {
            r_held_ = r_inner_; // the last tick went out uncut: the inner loop could follow its reference
        }
        r_inner_ = outer_.control(r, y_outer);
        u_cmd_ = inner_.control(r_inner_, y_inner);
        cut_ = false;
        if constexpr (requires { inner_.u_lim.lower(0); }) {
            if (!(inner_.u_lim.lower(0) < u_cmd_ && u_cmd_ < inner_.u_lim.upper(0))) {
                hold_outer(); // the inner command is on its limit
            }
        }
        return u_cmd_;
    }

    /**
     * @brief Record the input the plant received, @p u_sat, when a stage after the cascade
     *        clamped its command, and the feedforward @p u_ff added to it before that clamp
     *
     * The inner law takes it; when the clamp cut the command off, the outer law holds as
     * against the inner law's own limit.
     */
    constexpr void back_calculate(T u_sat, T u_ff = T{0}) {
        inner_.back_calculate(u_sat, u_ff);
        const T asked = u_cmd_ + u_ff;
        if (!cut_ && damp::abs(u_sat - asked) > T{4} * eps<T>() * (damp::abs(u_sat) + damp::abs(asked))) {
            hold_outer();
        }
    }

    constexpr void reset() {
        outer_.reset();
        inner_.reset();
        r_held_ = T{0};
        r_inner_ = T{0};
        u_cmd_ = T{0};
        cut_ = false;
    }

    /** @brief The command the last control() returned. */
    [[nodiscard]] constexpr T last_command() const { return u_cmd_; }

    [[nodiscard]] constexpr Outer&       outer() { return outer_; }
    [[nodiscard]] constexpr const Outer& outer() const { return outer_; }
    [[nodiscard]] constexpr Inner&       inner() { return inner_; }
    [[nodiscard]] constexpr const Inner& inner() const { return inner_; }

private:
    /** @brief The inner command was cut off this tick: the outer command holds at the last one the inner loop could follow. */
    constexpr void hold_outer() {
        cut_ = true;
        outer_.back_calculate(r_held_, T{0});
    }

    Outer outer_{};
    Inner inner_{};
    T     r_held_{T{0}};  ///< The outer command of the last tick that went out uncut: where a clamp of the inner command holds it
    T     r_inner_{T{0}}; ///< The outer command of the last tick: the inner law's reference
    T     u_cmd_{T{0}};   ///< The command the last control() returned
    bool  cut_{false};    ///< The last tick's inner command was cut off
};

/**
 * @brief Two PI (or any Siso) loops in cascade — roadmap name CascadePPI.
 */
template<typename Outer, typename Inner, typename T = float>
using CascadePPI = Cascade<Outer, Inner, T>;

/**
 * @brief Which path owns the plant command.
 */
enum class SwitchMode : uint8_t {
    Normal,     ///< Primary closed-loop controller
    Experiment, ///< Commissioning / autotune experiment owns u
    Backup,     ///< Safe fallback controller
};

/**
 * @brief Switch between normal, experiment, and backup SISO laws.
 *
 * All three slots are SisoControllers. On mode change the newly active law is @c reset(), so
 * integrators do not inherit foreign state. A law with a tracking mode (`track(u, y)` and
 * `enable()`, as PIDController) is then preloaded at the next tick from the command the plant
 * last received, at that tick's measurement, so it resumes from it: the command does not jump
 * but by its proportional term on the error. A law without one starts from its reset state.
 * back_calculate goes to the active law. The experiment path is still a controller shape: wrap
 * relay/FRF drivers with a small adapter if needed, or set mode to Experiment and call
 * @ref inject when u comes from outside.
 *
 * @see "Advanced PID Control" (Åström & Hägglund, 2006), §3.5: bumpless transfer by tracking
 *
 * @tparam Normal     Primary controller
 * @tparam Experiment Experiment-period controller (or dummy passthrough)
 * @tparam Backup     Safe fallback
 * @tparam T          Scalar type
 */
template<typename Normal, typename Experiment, typename Backup, typename T = float>
    requires SisoController<Normal, T> && SisoController<Experiment, T> && SisoController<Backup, T>
class SwitchedController {
public:
    constexpr SwitchedController() = default;

    constexpr SwitchedController(const Normal& n, const Experiment& e, const Backup& b)
        : normal_(n), experiment_(e), backup_(b) {}

    [[nodiscard]] constexpr SwitchMode mode() const { return mode_; }

    /**
     * @brief Select active path; resets the newly selected controller, and preloads one with a
     *        tracking mode at the next tick.
     */
    constexpr void set_mode(SwitchMode m) {
        if (m == mode_) {
            return;
        }
        mode_ = m;
        on_active([](auto& law) { law.reset(); });
        handover_ = true;
    }

    /**
     * @brief One tick of the active controller.
     */
    [[nodiscard]] constexpr T control(T r, T y) {
        on_active([&](auto& law) {
            if constexpr (requires { law.track(applied_, y); law.enable(); }) {
                if (handover_) {
                    law.track(applied_, y); // the integrator takes what the plant last received
                    law.enable();
                }
            }
            last_u_ = law.control(r, y);
        });
        handover_ = false;
        applied_ = last_u_;
        return last_u_;
    }

    /**
     * @brief Record the input the plant received, @p u_sat, after a clamp downstream, and the
     *        feedforward @p u_ff added before it: the active law takes it, and a law switched in
     *        next starts from u_sat − u_ff.
     */
    constexpr void back_calculate(T u_sat, T u_ff = T{0}) {
        on_active([&](auto& law) { law.back_calculate(u_sat, u_ff); });
        applied_ = u_sat - u_ff;
    }

    /**
     * @brief Override the command when the experiment is an external command source.
     *
     * Only applied in Experiment mode; stores @p u as last command.
     */
    constexpr void inject(T u) {
        if (mode_ == SwitchMode::Experiment) {
            last_u_ = u;
            applied_ = u;
        }
    }

    /** @brief The command the last control() returned, or the one inject() set, as Cascade's and each law's. */
    [[nodiscard]] constexpr T last_command() const { return last_u_; }

    constexpr void reset() {
        normal_.reset();
        experiment_.reset();
        backup_.reset();
        last_u_ = T{0};
        applied_ = T{0};
        handover_ = false;
    }

    [[nodiscard]] constexpr Normal&           normal() { return normal_; }
    [[nodiscard]] constexpr const Normal&     normal() const { return normal_; }
    [[nodiscard]] constexpr Experiment&       experiment() { return experiment_; }
    [[nodiscard]] constexpr const Experiment& experiment() const { return experiment_; }
    [[nodiscard]] constexpr Backup&           backup() { return backup_; }
    [[nodiscard]] constexpr const Backup&     backup() const { return backup_; }

private:
    /** @brief Call @p f on the law of the active mode. */
    template<typename F>
    constexpr void on_active(F&& f) {
        switch (mode_) {
            case SwitchMode::Experiment:
                f(experiment_);
                break;
            case SwitchMode::Backup:
                f(backup_);
                break;
            case SwitchMode::Normal:
            default:
                f(normal_);
                break;
        }
    }

    Normal     normal_{};
    Experiment experiment_{};
    Backup     backup_{};
    SwitchMode mode_{SwitchMode::Normal};
    T          last_u_{T{0}};    ///< The command the last control() returned, or inject's
    T          applied_{T{0}};   ///< What the plant last received from the switch: a law switched in starts from it
    bool       handover_{false}; ///< A law was switched in and has not ticked yet
};

} // namespace damp
