// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file experiment_safety.hpp
 * @brief Thin experiment policy: snapshot → protect commands → commit or rollback.
 *
 * Does not reinvent limits. Composes existing primitives:
 * - Bounds for command saturation
 * - optional slew rate on the protected command
 * - Timeout for experiment duration
 *
 * Use around sysid / autotune: arm with a snapshot of controller state (gains,
 * compensator bank, …), run the experiment under protect, then
 * commit or rollback.
 *
 * @see toolbox/bounds.hpp, toolbox/timing.hpp
 */

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/math/math.hpp"
#include "damp/toolbox/bounds.hpp"
#include "damp/toolbox/timing.hpp"

namespace damp {

/**
 * @brief Configuration for @ref ExperimentSafety.
 *
 * @tparam T Scalar type (command path).
 */
template<typename T = float>
struct ExperimentSafetyConfig {
    Bounds<1, T> command_bounds{T{-1}, T{1}}; ///< Hard clamp on protected u
    T            max_slew_rate{T{0}};         ///< |du/dt| cap; 0 = no slew limit
    T            timeout_s{T{0}};             ///< Experiment timeout; 0 = no timeout

    [[nodiscard]] constexpr bool valid() const {
        if (command_bounds.lower[0] > command_bounds.upper[0]) {
            return false;
        }
        if (!(max_slew_rate >= T{0}) || !(timeout_s >= T{0})) {
            return false;
        }
        return true;
    }
};

/**
 * @brief Experiment lifecycle + command protection for on-target commissioning.
 *
 * @tparam NSnapshot Fixed snapshot size (gains, flags, …) as scalar pack
 * @tparam T         Scalar type
 */
template<std::size_t NSnapshot, typename T = float>
class ExperimentSafety {
public:
    constexpr ExperimentSafety() = default;

    constexpr explicit ExperimentSafety(const ExperimentSafetyConfig<T>& cfg)
        : cfg_(cfg), valid_(cfg.valid()) {
        if (cfg_.timeout_s > T{0}) {
            timeout_ = Timeout<T>{cfg_.timeout_s};
        }
    }

    /**
     * @brief Begin an experiment: copy @p state into the rollback snapshot.
     */
    constexpr void arm(const damp::array<T, NSnapshot>& state) {
        snapshot_ = state;
        working_ = state;
        armed_ = valid_;
        committed_ = false;
        u_prev_ = T{0};
        have_u_ = false;
        timeout_.reset();
    }

    /// True after arm and before commit/rollback clears the run.
    [[nodiscard]] constexpr bool armed() const { return armed_; }

    /**
     * @brief Saturate and optional slew-limit a command sample.
     * @param u_cmd Raw experiment / controller command
     * @param dt    Sample period [s] (needed for slew)
     */
    [[nodiscard]] constexpr T protect(T u_cmd, T dt) {
        if (!valid_) {
            return T{0};
        }
        T u = cfg_.command_bounds.saturate(u_cmd);
        // Slew from last protected output (0 until the first call after arm).
        if (cfg_.max_slew_rate > T{0} && dt > T{0}) {
            const T max_du = cfg_.max_slew_rate * dt;
            const T du = u - u_prev_;
            if (du > max_du) {
                u = u_prev_ + max_du;
            } else if (du < -max_du) {
                u = u_prev_ - max_du;
            }
            u = cfg_.command_bounds.saturate(u);
        }
        u_prev_ = u;
        have_u_ = true;
        return u;
    }

    /**
     * @brief Advance the experiment timeout.
     * @return true if the timeout has expired (caller should rollback)
     */
    [[nodiscard]] constexpr bool tick_timeout(T dt) {
        if (!valid_ || !(cfg_.timeout_s > T{0}) || !armed_) {
            return false;
        }
        return timeout_.tick(dt);
    }

    [[nodiscard]] constexpr bool timed_out() const {
        return (cfg_.timeout_s > T{0}) && timeout_.expired();
    }

    /// Accept working state as final (clears armed).
    constexpr void commit() {
        armed_ = false;
        committed_ = true;
    }

    /**
     * @brief Restore snapshot into @p state and clear armed.
     */
    constexpr void rollback(damp::array<T, NSnapshot>& state) {
        state = snapshot_;
        working_ = snapshot_;
        armed_ = false;
        committed_ = false;
    }

    [[nodiscard]] constexpr const damp::array<T, NSnapshot>& snapshot() const { return snapshot_; }
    [[nodiscard]] constexpr damp::array<T, NSnapshot>&       working() { return working_; }
    [[nodiscard]] constexpr const damp::array<T, NSnapshot>& working() const { return working_; }
    [[nodiscard]] constexpr bool                             committed() const { return committed_; }

    [[nodiscard]] constexpr const ExperimentSafetyConfig<T>& config() const { return cfg_; }

private:
    ExperimentSafetyConfig<T> cfg_{};
    damp::array<T, NSnapshot> snapshot_{};
    damp::array<T, NSnapshot> working_{};
    Timeout<T>                timeout_{};
    T                         u_prev_{T{0}};
    bool                      have_u_{false};
    bool                      armed_{false};
    bool                      committed_{false};
    bool                      valid_{false};
};

/// Alias matching the roadmap name (same type).
template<std::size_t NSnapshot, typename T = float>
using SafetyEnvelope = ExperimentSafety<NSnapshot, T>;

} // namespace damp
