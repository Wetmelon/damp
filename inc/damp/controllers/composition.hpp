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
 * Both stay on the SisoController concept so they nest and drop into
 * simulate_sampled / multi-rate harnesses. Experiments (relay, FRF drivers)
 * that are not SisoControllers are not stored inside SwitchedController;
 * call SwitchedController::set_mode and feed the experiment yourself, or
 * wrap a thin SisoController adapter that holds the experiment output.
 *
 * @see concepts.hpp SisoController
 * @see estimation/experiment_safety.hpp for arm/protect/rollback around a switch
 */

#include <cstddef>
#include <cstdint>

#include "damp/concepts.hpp"

namespace damp {

/**
 * @brief Series cascade of two SISO controllers: outer → inner reference.
 *
 * Tick:
 *
 *     u_outer = outer.control(r, y_outer)
 *     u       = inner.control(u_outer, y_inner)
 *
 * Typical use: outer position / voltage loop, inner velocity / current.
 * Alias CascadePPI for PI+PI when both are PIDController PI modes.
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
        const T r_inner = outer_.control(r, y_outer);
        return inner_.control(r_inner, y_inner);
    }

    /// SisoController-shaped tick when outer and inner share the same y (rare).
    [[nodiscard]] constexpr T control(T r, T y) { return control(r, y, y); }

    constexpr void reset() {
        outer_.reset();
        inner_.reset();
    }

    [[nodiscard]] constexpr Outer&       outer() { return outer_; }
    [[nodiscard]] constexpr const Outer& outer() const { return outer_; }
    [[nodiscard]] constexpr Inner&       inner() { return inner_; }
    [[nodiscard]] constexpr const Inner& inner() const { return inner_; }

private:
    Outer outer_{};
    Inner inner_{};
};

/**
 * @brief Two PI (or any Siso) loops in cascade — roadmap name CascadePPI.
 */
template<typename Outer, typename Inner, typename T = float>
using CascadePPI = Cascade<Outer, Inner, T>;

/**
 * @brief Which path owns the plant command.
 */
enum class SwitchMode : std::uint8_t {
    Normal,     ///< Primary closed-loop controller
    Experiment, ///< Commissioning / autotune experiment owns u
    Backup,     ///< Safe fallback controller
};

/**
 * @brief Bumpless-ish switch between normal, experiment, and backup SISO laws.
 *
 * All three slots are SisoControllers. On mode change, the newly active
 * controller is @c reset() so integrators do not inherit foreign state
 * (Hanus-style full bumpless needs more; this is the library's lightweight
 * default). Experiment path is still a controller shape: wrap relay/FRF
 * drivers with a small adapter if needed, or set mode to Experiment and
 * call @ref inject when u comes from outside.
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
     * @brief Select active path; resets the newly selected controller.
     */
    constexpr void set_mode(SwitchMode m) {
        if (m == mode_) {
            return;
        }
        mode_ = m;
        active_reset();
    }

    /**
     * @brief One tick of the active controller.
     */
    [[nodiscard]] constexpr T control(T r, T y) {
        switch (mode_) {
            case SwitchMode::Experiment:
                last_u_ = experiment_.control(r, y);
                break;
            case SwitchMode::Backup:
                last_u_ = backup_.control(r, y);
                break;
            case SwitchMode::Normal:
            default:
                last_u_ = normal_.control(r, y);
                break;
        }
        return last_u_;
    }

    /**
     * @brief Override the command when the experiment is external (not a SisoController).
     *
     * Only applied in Experiment mode; stores @p u as last command.
     */
    constexpr void inject(T u) {
        if (mode_ == SwitchMode::Experiment) {
            last_u_ = u;
        }
    }

    [[nodiscard]] constexpr T last_u() const { return last_u_; }

    constexpr void reset() {
        normal_.reset();
        experiment_.reset();
        backup_.reset();
        last_u_ = T{0};
    }

    [[nodiscard]] constexpr Normal&           normal() { return normal_; }
    [[nodiscard]] constexpr const Normal&     normal() const { return normal_; }
    [[nodiscard]] constexpr Experiment&       experiment() { return experiment_; }
    [[nodiscard]] constexpr const Experiment& experiment() const { return experiment_; }
    [[nodiscard]] constexpr Backup&           backup() { return backup_; }
    [[nodiscard]] constexpr const Backup&     backup() const { return backup_; }

private:
    constexpr void active_reset() {
        switch (mode_) {
            case SwitchMode::Experiment:
                experiment_.reset();
                break;
            case SwitchMode::Backup:
                backup_.reset();
                break;
            case SwitchMode::Normal:
            default:
                normal_.reset();
                break;
        }
    }

    Normal     normal_{};
    Experiment experiment_{};
    Backup     backup_{};
    SwitchMode mode_{SwitchMode::Normal};
    T          last_u_{T{0}};
};

} // namespace damp
