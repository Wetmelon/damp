// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once


/**
 * @file ramp.hpp
 * @brief Ramp excitation — design payload and runtime generator
 */


#include "damp/backend.hpp"
#include "damp/math/math.hpp"

namespace damp {
namespace design {

/**
 * @struct RampConfig
 * @brief Configuration for a slew-rate-limited ramp excitation.
 *
 * Ramps from 0 toward `target` at fixed slope `rate`, then optionally holds
 * the final value for `hold_at_end_s`.
 *
 * @tparam T Scalar type.
 */
template<typename T = double>
struct RampConfig {
    T target{T{1}};        ///< Final target value (finite)
    T rate{T{1}};          ///< Ramp rate in units/second (> 0)
    T hold_at_end_s{T{0}}; ///< Post-ramp hold duration in seconds (>= 0)

    /**
     * @brief Validate ramp configuration.
     * @return true if target/rate/hold are finite and valid.
     */
    [[nodiscard]] constexpr bool valid() const {
        if (!damp::isfinite(target)) {
            return false;
        }
        if (!damp::finite_positive(rate)) {
            return false;
        }
        if (!damp::finite_non_negative(hold_at_end_s)) {
            return false;
        }
        return true;
    }
};

/**
 * @struct RampResult
 * @brief Ramp design payload.
 *
 * Carries validated ramp settings into the runtime generator.
 *
 * @tparam T Scalar type.
 */
template<typename T = double>
struct RampResult {
    RampConfig<T> config{};       ///< Validated ramp configuration
    bool          success{false}; ///< true if `config.valid()`
    T             Ts{0};          ///< [s] Sample period of the internal step() clock (0 = absolute-time step(t) only)

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /**
     * @brief Convert design payload to another scalar type.
     * @tparam U Target scalar type.
     * @return RampResult\<U\> converted element-wise.
     */
    template<typename U>
    [[nodiscard]] constexpr RampResult<U> as() const {
        return RampResult<U>{
            RampConfig<U>{
                static_cast<U>(config.target),
                static_cast<U>(config.rate),
                static_cast<U>(config.hold_at_end_s),
            },
            success,
            static_cast<U>(Ts),
        };
    }
};

/**
 * @brief Build a ramp design payload from a configuration.
 *
 * @param config Ramp configuration.
 * @param Ts     [s] Sample period of the runtime's step() clock (0 = absolute-time step(t) only; < 0 refused)
 * @return RampResult with `success` = `config.valid()` and a valid Ts.
 */
template<typename T = double>
[[nodiscard]] constexpr RampResult<T>
ramp(const RampConfig<T>& config, T Ts = T{0}) {
    return RampResult<T>{config, config.valid() && (Ts >= T{0}) && damp::isfinite(Ts), Ts};
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Rate-limited ramp runtime generator.
 *
 * Ramps from 0 to `target` at magnitude `rate` (units/s), then holds the
 * target for `hold_at_end_s` before reporting done.
 *
 * @tparam T Scalar type.
 */
template<typename T = float>
class Ramp {
public:
    constexpr Ramp() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr Ramp(const design::RampResult<U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : Ramp(result.template as<T>()) {}

    /**
     * @brief Construct from a ramp design payload.
     * @param design Validated design payload.
     */
    constexpr Ramp(const design::RampResult<T>& design) // NOLINT(google-explicit-constructor): assign-from-design
        : config_(design.config), Ts_(design.Ts), valid_(design.success) {
        if (!design.success && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    /**
     * @brief Evaluate ramp at absolute time.
     * @param t Time in seconds.
     * @return Slew-limited value, saturated at target after ramp duration.
     */
    [[nodiscard]] constexpr T step(T t) const {
        if (!valid_) {
            return T{0};
        }

        const T tc = damp::max(t, T{0});
        const T tr = ramp_duration();
        if (tc >= tr) {
            return config_.target;
        }

        const T slope = damp::copysign(config_.rate, config_.target);
        return slope * tc;
    }

    /**
     * @brief Evaluate ramp at internal time and advance by `Ts`.
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
     * @brief Query whether ramp and hold interval are complete at absolute time.
     * @param t Time in seconds.
     * @return true if invalid or `t >= ramp_duration() + hold_at_end_s`.
     */
    [[nodiscard]] constexpr bool done(T t) const {
        if (!valid_) {
            return true;
        }
        return t >= (ramp_duration() + config_.hold_at_end_s);
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

    /**
     * @brief Return the ramp-only duration.
     * @return `abs(target)/rate` for valid configuration.
     */
    [[nodiscard]] constexpr T ramp_duration() const {
        if (!valid_) {
            return T{0};
        }
        return damp::abs(config_.target) / config_.rate;
    }

private:
    /** @brief Internal time k·Ts (a sample count, not a running sum of Ts). */
    [[nodiscard]] constexpr T elapsed() const { return static_cast<T>(k_) * Ts_; }

    design::RampConfig<T> config_{};
    T                     Ts_{T{0}};
    size_t                k_{0}; ///< Samples emitted by step()
    bool                  valid_{false};
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<typename T>
Ramp(const design::RampResult<T>&) -> Ramp<T>;

} // namespace damp
