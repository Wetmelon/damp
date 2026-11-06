// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file step_train.hpp
 * @brief Step-train excitation — design payload and runtime generator
 */

#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/math/math.hpp"

namespace damp {
namespace design {

/**
 * @struct StepTrainConfig
 * @brief Configuration for alternating +/- step excitation.
 *
 * One cycle is two plateaus: `+amplitude` for `hold_s`, then `-amplitude`
 * for `hold_s`.
 *
 * @tparam T Scalar type.
 */
template<typename T = double>
struct StepTrainConfig {
    T           amplitude{T{1}};             ///< Step magnitude (> 0)
    T           hold_s{static_cast<T>(0.1)}; ///< Hold time per plateau in seconds (> 0)
    std::size_t cycles{1};                   ///< Number of +/- pairs

    /**
     * @brief Validate step-train configuration.
     * @return true if amplitude, hold, and cycle count are valid.
     */
    [[nodiscard]] constexpr bool valid() const {
        if (!damp::finite_positive(amplitude)) {
            return false;
        }
        if (!damp::finite_positive(hold_s)) {
            return false;
        }
        if (cycles == 0) {
            return false;
        }
        return true;
    }
};

/**
 * @struct StepTrainResult
 * @brief Step-train design payload.
 *
 * Carries validated step-train settings into the runtime generator.
 *
 * @tparam T Scalar type.
 */
template<typename T = double>
struct StepTrainResult {
    StepTrainConfig<T> config{};       ///< Validated step-train configuration
    bool               success{false}; ///< true if `config.valid()`
    T                  Ts{0};          ///< [s] Sample period of the internal step() clock (0 = absolute-time step(t) only)

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /**
     * @brief Convert design payload to another scalar type.
     * @tparam U Target scalar type.
     * @return StepTrainResult\<U\> converted element-wise.
     */
    template<typename U>
    [[nodiscard]] constexpr StepTrainResult<U> as() const {
        return StepTrainResult<U>{
            StepTrainConfig<U>{
                static_cast<U>(config.amplitude),
                static_cast<U>(config.hold_s),
                config.cycles,
            },
            success,
            static_cast<U>(Ts),
        };
    }
};

/**
 * @brief Build a step-train design payload from a configuration.
 *
 * @param config Step-train configuration.
 * @param Ts     [s] Sample period of the runtime's step() clock (0 = absolute-time step(t) only; < 0 refused)
 * @return StepTrainResult with `success` = `config.valid()` and a valid Ts.
 */
template<typename T = double>
[[nodiscard]] constexpr StepTrainResult<T>
step_train(const StepTrainConfig<T>& config, T Ts = T{0}) {
    return StepTrainResult<T>{config, config.valid() && (Ts >= T{0}) && damp::isfinite(Ts), Ts};
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Alternating +/- step train runtime generator.
 *
 * Emits piecewise-constant plateaus of `+amplitude` and `-amplitude` with
 * hold time `hold_s`, for `cycles` complete +/- pairs.
 *
 * @tparam T Scalar type.
 */
template<typename T = float>
class StepTrain {
public:
    constexpr StepTrain() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!std::is_same_v<U, T>)
    constexpr StepTrain(const design::StepTrainResult<U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : StepTrain(result.template as<T>()) {}

    /**
     * @brief Construct from a step-train design payload.
     * @param design Validated design payload.
     */
    constexpr StepTrain(const design::StepTrainResult<T>& design) // NOLINT(google-explicit-constructor): assign-from-design
        : config_(design.config), Ts_(design.Ts), valid_(design.success) {
        if (!design.success && std::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    /**
     * @brief Evaluate step train at absolute time.
     * @param t Time in seconds.
     * @return +/- amplitude while active, or 0 after completion.
     */
    [[nodiscard]] constexpr T step(T t) const {
        if (!valid_) {
            return T{0};
        }

        const T tc = damp::max(t, T{0});
        if (done(tc)) {
            return T{0};
        }

        const auto segment = static_cast<std::size_t>(damp::floor(tc / config_.hold_s));
        return (segment % 2u == 0u) ? config_.amplitude : -config_.amplitude;
    }

    /**
     * @brief Evaluate step train at internal time and advance by `Ts`.
     * @return Excitation value at current internal time.
     */
    [[nodiscard]] constexpr T step() {
        const T out = step(elapsed());
        if (Ts_ > T{0} && !done()) {
            ++k_;
        }
        return out;
    }

    /**
     * @brief Query whether configured cycles are complete at absolute time.
     * @param t Time in seconds.
     * @return true if invalid or `t >= 2*cycles*hold_s`.
     */
    [[nodiscard]] constexpr bool done(T t) const {
        if (!valid_) {
            return true;
        }
        const T total = static_cast<T>(2u * config_.cycles) * config_.hold_s;
        return t >= total;
    }

    /**
     * @brief Query completion at internal time.
     * @return true if internal sequence is complete.
     */
    [[nodiscard]] constexpr bool done() const {
        return done(elapsed());
    }

    /**
     * @brief Reset internal time to zero.
     */
    constexpr void reset() {
        k_ = 0;
    }

private:
    /** @brief Internal time k·Ts (a sample count, not a running sum of Ts). */
    [[nodiscard]] constexpr T elapsed() const { return static_cast<T>(k_) * Ts_; }

    design::StepTrainConfig<T> config_{};
    T                          Ts_{T{0}};
    std::size_t                k_{0}; ///< Samples emitted by step()
    bool                       valid_{false};
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<typename T>
StepTrain(const design::StepTrainResult<T>&) -> StepTrain<T>;

} // namespace damp
