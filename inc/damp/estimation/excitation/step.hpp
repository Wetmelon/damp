// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once


/**
 * @file step.hpp
 * @brief Step excitation — design payload and runtime generator
 */


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
struct [[nodiscard]] StepResult {
    StepConfig<T> config{};
    bool          success{false};
    T             Ts{0}; ///< [s] Sample period of the internal step() clock (> 0)

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

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
            static_cast<U>(Ts),
        };
    }
};

/**
 * @brief Build a single-step design payload.
 *
 * @param config Step configuration.
 * @param Ts     [s] Sample period of the runtime's step() clock (> 0: Step has no absolute-time step(t))
 */
template<typename T = double>
[[nodiscard]] constexpr StepResult<T>
step(const StepConfig<T>& config, T Ts) {
    return StepResult<T>{config, config.valid() && (Ts > T{0}) && damp::isfinite(Ts), Ts};
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

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr Step(const design::StepResult<U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : Step(result.template as<T>()) {}

    constexpr Step(const design::StepResult<T>& design) // NOLINT(google-explicit-constructor): assign-from-design
        : config_(design.config), Ts_(design.Ts), valid_(design.success) {
        if (!design.success && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    /**
     * @brief Return the step command and advance internal time by `Ts`.
     */
    [[nodiscard]] constexpr T step() {
        if (!valid_ || finished_) {
            return T{0};
        }
        const T out = config_.amplitude;
        ++k_;
        update_done_from_time();
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
        if (time() < config_.min_hold_s) {
            settle_k_ = 0;
            return;
        }
        if (rate <= config_.settle_rate_tol) {
            ++settle_k_;
            if (static_cast<T>(settle_k_) * Ts_ >= config_.settle_confirm_s) {
                finished_ = true;
                settled_ok_ = true;
            }
        } else {
            settle_k_ = 0;
        }
        update_done_from_time();
    }

    [[nodiscard]] constexpr bool done() const {
        return !valid_ || finished_;
    }

    /** @brief true if completion was via settle band. */
    [[nodiscard]] constexpr bool settled() const { return settled_ok_; }

    constexpr void reset() {
        k_ = 0;
        settle_k_ = 0;
        finished_ = false;
        settled_ok_ = false;
        have_y_prev_ = false;
        y_prev_ = T{0};
    }

    /** @brief Internal time k·Ts [s] after k samples (a count, not a running sum of Ts). */
    [[nodiscard]] constexpr T time() const { return static_cast<T>(k_) * Ts_; }

private:
    constexpr void update_done_from_time() {
        if (finished_) {
            return;
        }
        if (config_.settle_rate_tol > T{0}) {
            if (time() >= config_.max_hold_s) {
                finished_ = true;
            }
        } else if (time() >= config_.min_hold_s) {
            finished_ = true;
            settled_ok_ = true; // timed step: treat as nominal completion
        }
    }

    design::StepConfig<T> config_{};
    T                     Ts_{T{0}};
    size_t                k_{0};        ///< Samples emitted by step()
    size_t                settle_k_{0}; ///< Consecutive quiet samples observed
    T                     y_prev_{T{0}};
    bool                  have_y_prev_{false};
    bool                  valid_{false};
    bool                  finished_{false};
    bool                  settled_ok_{false};
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<typename T>
Step(const design::StepResult<T>&) -> Step<T>;

} // namespace damp
