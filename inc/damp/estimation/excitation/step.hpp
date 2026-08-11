// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file step.hpp
 * @brief Step excitation — design payload and runtime generator
 */

#include <cstddef>
#include <cstdint>

#include "damp/backend.hpp"
#include "damp/math/math.hpp"

namespace damp {
namespace design {

/**
 * @struct StepConfig
 * @brief Configuration for a single step with optional settle detection.
 *
 * Holds `amplitude` from t = 0. Completion:
 * - If `settle_rate_tol <= 0`: done after `min_hold_s` (pure timed step).
 * - If `settle_rate_tol > 0`: call runtime `observe(y)` each tick; done when
 *   |dy/dt| stays below the tolerance for `settle_confirm_s` after
 *   `min_hold_s`, or when `max_hold_s` is reached (timeout).
 *
 * @tparam T Scalar type.
 */
template<typename T = double>
struct StepConfig {
    T amplitude{T{1}};                        ///< Step level (finite)
    T min_hold_s{static_cast<T>(0.1)};        ///< Minimum time before settle can fire (> 0)
    T max_hold_s{T{10}};                      ///< Timeout (>= min_hold_s)
    T settle_rate_tol{T{0}};                  ///< |dy/dt| band for settle; 0 = time-only
    T settle_confirm_s{static_cast<T>(0.05)}; ///< Time |dy/dt| must stay in band (>= 0)

    [[nodiscard]] constexpr bool valid() const {
        if (!damp::isfinite(amplitude)) {
            return false;
        }
        if (!damp::finite_positive(min_hold_s)) {
            return false;
        }
        if (!damp::isfinite(max_hold_s) || max_hold_s < min_hold_s) {
            return false;
        }
        if (!damp::finite_non_negative(settle_rate_tol)) {
            return false;
        }
        if (!damp::finite_non_negative(settle_confirm_s)) {
            return false;
        }
        return true;
    }
};

/**
 * @struct StepResult
 * @brief Single-step design payload.
 */
template<typename T = double>
struct StepResult {
    StepConfig<T> config{};
    bool          success{false};

    template<typename U>
    [[nodiscard]] constexpr StepResult<U> as() const {
        return StepResult<U>{
            StepConfig<U>{
                static_cast<U>(config.amplitude),
                static_cast<U>(config.min_hold_s),
                static_cast<U>(config.max_hold_s),
                static_cast<U>(config.settle_rate_tol),
                static_cast<U>(config.settle_confirm_s),
            },
            success,
        };
    }
};

/**
 * @brief Build a single-step design payload.
 */
template<typename T = double>
[[nodiscard]] constexpr StepResult<T>
step(const StepConfig<T>& config) {
    return StepResult<T>{config, config.valid()};
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Single step with optional |dy/dt| settle detection.
 *
 * step always returns the configured amplitude until @ref done.
 * For settle mode (`settle_rate_tol > 0`), call @ref observe with the plant
 * output each tick so the generator can estimate |dy/dt|.
 *
 * @tparam T Scalar type.
 */
template<typename T = float>
class Step {
public:
    constexpr Step() = default;

    constexpr explicit Step(const design::StepResult<T>& design, T Ts = T{0})
        : config_(design.config), Ts_(Ts), valid_(design.success) {}

    /**
     * @brief Return the step command and advance internal time by `Ts`.
     */
    [[nodiscard]] constexpr T step() {
        if (!valid_ || finished_) {
            return T{0};
        }
        const T out = config_.amplitude;
        if (Ts_ > T{0}) {
            t_ += Ts_;
            update_done_from_time();
        }
        return out;
    }

    /**
     * @brief Feed a plant measurement for settle detection.
     *
     * Ignored when `settle_rate_tol <= 0` or before the first sample.
     * Estimates |dy/dt| ≈ |y − y_prev| / Ts.
     */
    constexpr void observe(T y) {
        if (!valid_ || finished_ || !(config_.settle_rate_tol > T{0}) || !(Ts_ > T{0})) {
            y_prev_ = y;
            have_y_prev_ = true;
            return;
        }
        if (!have_y_prev_) {
            y_prev_ = y;
            have_y_prev_ = true;
            return;
        }
        const T rate = damp::abs(y - y_prev_) / Ts_;
        y_prev_ = y;
        if (t_ < config_.min_hold_s) {
            settle_hold_t_ = T{0};
            return;
        }
        if (rate <= config_.settle_rate_tol) {
            settle_hold_t_ += Ts_;
            if (settle_hold_t_ >= config_.settle_confirm_s) {
                finished_ = true;
                settled_ok_ = true;
            }
        } else {
            settle_hold_t_ = T{0};
        }
        update_done_from_time();
    }

    [[nodiscard]] constexpr bool done() const {
        return !valid_ || finished_;
    }

    /// true if completion was via settle band.
    [[nodiscard]] constexpr bool settled() const { return settled_ok_; }

    constexpr void reset() {
        t_ = T{0};
        settle_hold_t_ = T{0};
        finished_ = false;
        settled_ok_ = false;
        have_y_prev_ = false;
        y_prev_ = T{0};
    }

    [[nodiscard]] constexpr T time() const { return t_; }

private:
    constexpr void update_done_from_time() {
        if (finished_) {
            return;
        }
        if (config_.settle_rate_tol > T{0}) {
            if (t_ >= config_.max_hold_s) {
                finished_ = true;
            }
        } else if (t_ >= config_.min_hold_s) {
            finished_ = true;
            settled_ok_ = true; // timed step: treat as nominal completion
        }
    }

    design::StepConfig<T> config_{};
    T                     Ts_{T{0}};
    T                     t_{T{0}};
    T                     settle_hold_t_{T{0}};
    T                     y_prev_{T{0}};
    bool                  have_y_prev_{false};
    bool                  valid_{false};
    bool                  finished_{false};
    bool                  settled_ok_{false};
};

} // namespace damp
