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
#include <cstdint>

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
        };
    }
};

/**
 * @brief Build a step-train design payload from a configuration.
 *
 * @param config Step-train configuration.
 * @return StepTrainResult with `success = config.valid()`.
 */
template<typename T = double>
[[nodiscard]] constexpr StepTrainResult<T>
step_train(const StepTrainConfig<T>& config) {
    return StepTrainResult<T>{config, config.valid()};
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

    /**
     * @brief Construct from a step-train design payload.
     * @param design Validated design payload.
     * @param Ts     Optional sample period for internal `step()` mode.
     */
    constexpr explicit StepTrain(const design::StepTrainResult<T>& design, T Ts = T{0})
        : config_(design.config), Ts_(Ts), valid_(design.success) {}

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
        const T out = step(t_);
        if (Ts_ > T{0} && !done()) {
            t_ += Ts_;
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
        return done(t_);
    }

    /**
     * @brief Reset internal time to zero.
     */
    constexpr void reset() {
        t_ = T{0};
    }

private:
    design::StepTrainConfig<T> config_{};
    T                          Ts_{T{0}};
    T                          t_{T{0}};
    bool                       valid_{false};
};

} // namespace damp
