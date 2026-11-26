// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file composition.hpp
 * @brief SISO composition: cascade and switched controller / experiment paths.
 *
 * - Cascade — outer.step(r, y_outer) → inner.step(outer_u, y_inner)
 * - SwitchedController — Normal | Experiment | Backup ownership of u
 *
 * SwitchedController stays on the SisoController concept so it nests and drops
 * into simulate_sampled / multi-rate harnesses, and is a SisoLaw itself. Cascade takes a distinct measurement for each loop
 * (step(r, y_outer, y_inner)) and is not a SisoController — one shared y would regulate the inner
 * loop to the outer measurement; its handover is track(u, r, y_outer, y_inner). Experiments
 * (relay, FRF drivers) that are not SisoControllers are not stored inside SwitchedController;
 * call SwitchedController::set_mode and feed the experiment yourself, or wrap a thin
 * SisoController adapter that holds the experiment output.
 *
 * @see concepts.hpp SisoController, SisoLaw
 */

#include "damp/backend.hpp"
#include "damp/concepts.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix_traits.hpp" // eps
#include "damp/toolbox/bounds.hpp"

namespace damp {

/**
 * @brief Series cascade of two SISO controllers: outer → inner reference.
 *
 * Tick:
 *
 *     u_outer = outer.step(r, y_outer)
 *     u       = inner.step(u_outer, y_inner)
 *
 * While the inner law's command is clamped, at its own output limit or by a clamp reported with
 * back_calculate, the inner loop cannot follow its reference, and the outer law is told its
 * command was held at the last one the inner loop could follow, the outer command of the last
 * tick that went out uncut: back_calculate of it. Its memory holds as against a clamp of its own
 * (an outer PI's integral gives back what each tick added past it), where it ran on for as long
 * as the inner loop could not follow and overshot when it could. The inner law's output limit is
 * its limits().input; an inner law without limits() tells only through back_calculate.
 *
 * | Member | What it is |
 * | --- | --- |
 * | step(r, y_outer, y_inner) | one tick: the outer law, then the inner law on its command |
 * | track(u, r, y_outer, y_inner) | in place of step while another source drives the plant |
 * | back_calculate(u_sat, u_ff) | a clamp after the cascade: the inner law takes it, the outer holds |
 * | limits(), set_limits(lim) | the inner law's: the cascade's command is the inner law's |
 * | last_command(), valid(), reset() | the last command; whether both laws run; reset both |
 * | outer(), inner() | the two laws: the outer's own limits are outer().limits() |
 *
 * The cascade's limits are the inner law's (limits() and set_limits() forward to it), since its
 * command is the inner law's; the outer law's limits, which bound the inner reference (a rate
 * command, a current command), are its own, `outer().set_limits(…)`.
 *
 * While another source drives the plant, track(u, r, y_outer, y_inner) takes the place of
 * step: as in manual mode on a process cascade, the inner reference follows the inner
 * measurement, so the outer law tracks y_inner at (r, y_outer) and the inner law tracks u at
 * (y_inner, y_inner). The next step at the same readings returns u when the outer law lands on
 * y_inner (an outer law with an integral to preload, inside its limits) and the inner law on u
 * (inside its limits); an outer law with nothing to preload starts from its own command.
 *
 * Typical use: outer position / voltage loop, inner velocity / current.
 * Alias CascadePPI for PI+PI when both loops are PIController.
 *
 * @see "Advanced PID Control" (Åström & Hägglund, 2006), §3.5: anti-windup by conditional
 *      integration and by tracking, and bumpless transfer
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
    [[nodiscard]] constexpr T step(T r, T y_outer, T y_inner) {
        if (!cut_) {
            r_held_ = r_inner_; // the last tick went out uncut: the inner loop could follow its reference
        }
        r_inner_ = outer_.step(r, y_outer);
        u_cmd_ = inner_.step(r_inner_, y_inner);
        cut_ = false;
        if constexpr (requires { inner_.limits().input.lower(0); }) {
            const auto lim = inner_.limits().input;
            if (!(lim.lower(0) < u_cmd_ && u_cmd_ < lim.upper(0))) {
                hold_outer(); // the inner command is on its limit
            }
        }
        return u_cmd_;
    }

    /**
     * @brief In place of step(@p r, @p y_outer, @p y_inner) while another source drives the plant, @p u the input applied
     *
     * The inner reference follows the inner measurement: the outer law tracks @p y_inner at
     * (@p r, @p y_outer), the inner law tracks @p u at (@p y_inner, @p y_inner), so the next step
     * at these readings starts from u (see the class). @p u is recorded as the command applied.
     */
    constexpr void track(T u, T r, T y_outer, T y_inner)
        requires SisoTracking<Outer, T> && SisoTracking<Inner, T>
    {
        outer_.track(y_inner, r, y_outer);
        inner_.track(u, y_inner, y_inner);
        r_held_ = y_inner;
        r_inner_ = y_inner;
        u_cmd_ = u;
        cut_ = false;
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

    /** @brief The inner law's limits: the cascade's command is the inner law's. */
    [[nodiscard]] constexpr auto limits() const
        requires HoldsLimits<Inner>
    {
        return inner_.limits();
    }

    /** @brief Set the inner law's limits at run time; refused (false, nothing changed) as the inner law refuses them. The outer law's are outer().set_limits(…). */
    template<typename Limits>
        requires HoldsLimits<Inner> && requires(Inner& in, const Limits& l) { in.set_limits(l); }
    [[nodiscard]] constexpr bool set_limits(const Limits& lim) {
        return inner_.set_limits(lim);
    }

    /** @brief Whether both laws run: each one's valid(), a law without one counted as running. */
    [[nodiscard]] constexpr bool valid() const { return runs(outer_) && runs(inner_); }

    constexpr void reset() {
        outer_.reset();
        inner_.reset();
        r_held_ = T{0};
        r_inner_ = T{0};
        u_cmd_ = T{0};
        cut_ = false;
    }

    /** @brief The command the last step() returned, or the input the last track() recorded. */
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

    /** @brief Whether @p law runs: its valid(), or true for a law without one. */
    template<typename L>
    [[nodiscard]] static constexpr bool runs(const L& law) {
        if constexpr (ReportsValidity<L>) {
            return law.valid();
        } else {
            static_cast<void>(law);
            return true;
        }
    }

    Outer outer_{};
    Inner inner_{};
    T     r_held_{T{0}};  ///< The outer command of the last tick that went out uncut: where a clamp of the inner command holds it
    T     r_inner_{T{0}}; ///< The outer command of the last tick: the inner law's reference
    T     u_cmd_{T{0}};   ///< The command the last step() returned, or the input track() recorded
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
 * @brief The limits of a SwitchedController: on the command, which every law it holds drives
 *
 * SwitchedController::limits reads the active law's; set_limits sets them on every law.
 */
template<typename T = float>
struct SwitchedLimits {
    Bounds<1, T> input{}; ///< The command u [lower, upper]

    /** @brief The same limits in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr SwitchedLimits<U> as() const {
        return SwitchedLimits<U>{input.template as<U>()};
    }
};

/**
 * @brief Switch between normal, experiment, and backup SISO laws.
 *
 * All three slots are SisoControllers. On mode change the newly active law is @c reset(), so
 * integrators do not inherit foreign state. A law with `track(u, r, y)` (every SisoLaw) then
 * tracks, at the next tick's reference and measurement, the command the plant last received, and
 * that tick's step starts from it: the first command is the last applied input wherever the law
 * has the memory to land there (a PID's integral, a super-twisting law's v). A law without track
 * starts from its reset state. back_calculate goes to the active law. The experiment path is
 * still a controller shape: wrap relay/FRF drivers with a small adapter if needed, or set mode to
 * Experiment and call @ref inject when u comes from outside.
 *
 * | Member | What it is |
 * | --- | --- |
 * | set_mode(m), mode() | the active path; a law switched in tracks the last applied input |
 * | step(r, y) | one tick of the active law |
 * | track(u, r, y) | in place of step while another source drives the plant: the active law tracks u |
 * | back_calculate(u_sat, u_ff) | a clamp after the switch: to the active law |
 * | limits(), set_limits(lim) | the command's limits (SwitchedLimits): the active law's; set on every law |
 * | inject(u), last_command(), valid(), reset() | an external command in Experiment mode; the last command; whether every law runs; reset all three |
 *
 * The three laws drive one actuator, so the switch's limit is the one on its command:
 * limits().input reads the active law's limits().input, and set_limits(lim) sets lim.input on
 * every law that holds limits, each law's other limits (a PID's integral) its own. A law's own
 * limits are still set on it directly, normal().set_limits(…).
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
     * @brief Select active path; resets the newly selected controller, and preloads one with
     *        track(u, r, y) at the next tick.
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
    [[nodiscard]] constexpr T step(T r, T y) {
        on_active([&](auto& law) {
            if constexpr (requires { law.track(applied_, r, y); }) {
                if (handover_) {
                    law.track(applied_, r, y); // the step below at this (r, y) returns what the plant last received
                }
            }
            last_u_ = law.step(r, y);
        });
        handover_ = false;
        applied_ = last_u_;
        return last_u_;
    }

    /**
     * @brief In place of step(@p r, @p y) while another source drives the plant, @p u the input applied
     *
     * The active law tracks @p u at (@p r, @p y), so its next step starts from u; a law without
     * track keeps its state. @p u is recorded as the command applied, and a law switched in next
     * starts from it.
     */
    constexpr void track(T u, T r, T y) {
        on_active([&](auto& law) {
            if constexpr (requires { law.track(u, r, y); }) {
                law.track(u, r, y);
            }
        });
        handover_ = false;
        last_u_ = u;
        applied_ = u;
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

    /** @brief The command the last step() returned, the one inject() set, or the input track() recorded. */
    [[nodiscard]] constexpr T last_command() const { return last_u_; }

    /** @brief The command's limits: the active law's limits().input (unbounded for a law without limits()). */
    [[nodiscard]] constexpr SwitchedLimits<T> limits() const {
        SwitchedLimits<T> out{};
        on_active([&](const auto& law) {
            if constexpr (requires { law.limits().input; }) {
                out.input = law.limits().input;
            }
        });
        return out;
    }

    /**
     * @brief Set the command's limits on every law that holds limits, each law's other limits kept
     *
     * Refused (false, nothing changed) when an end is NaN or ±Inf, or the lower end is above the
     * upper one.
     */
    [[nodiscard]] constexpr bool set_limits(const SwitchedLimits<T>& lim) {
        if (!lim.input.valid()) {
            return false;
        }
        const auto apply = [&](auto& law) {
            if constexpr (requires { law.limits().input; law.set_limits(law.limits()); }) {
                auto own = law.limits();
                own.input = lim.input;
                static_cast<void>(law.set_limits(own)); // its other limits are its own, already boxes
            }
        };
        apply(normal_);
        apply(experiment_);
        apply(backup_);
        return true;
    }

    /** @brief Whether every law runs: each one's valid(), a law without one counted as running. */
    [[nodiscard]] constexpr bool valid() const { return runs(normal_) && runs(experiment_) && runs(backup_); }

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

    /** @brief Call @p f on the law of the active mode, read only. */
    template<typename F>
    constexpr void on_active(F&& f) const {
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

    /** @brief Whether @p law runs: its valid(), or true for a law without one. */
    template<typename L>
    [[nodiscard]] static constexpr bool runs(const L& law) {
        if constexpr (ReportsValidity<L>) {
            return law.valid();
        } else {
            static_cast<void>(law);
            return true;
        }
    }

    Normal     normal_{};
    Experiment experiment_{};
    Backup     backup_{};
    SwitchMode mode_{SwitchMode::Normal};
    T          last_u_{T{0}};    ///< The command the last step() returned, inject's, or the input track() recorded
    T          applied_{T{0}};   ///< What the plant last received from the switch: a law switched in starts from it
    bool       handover_{false}; ///< A law was switched in and has not ticked yet
};

} // namespace damp
