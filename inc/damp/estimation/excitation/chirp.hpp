// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file chirp.hpp
 * @brief Chirp excitation — design payload and runtime generator
 */

#include <cstddef>
#include <cstdint>
#include <limits>

#include "damp/backend.hpp"
#include "damp/math/math.hpp"

namespace damp {
namespace design {

/**
 * @brief Chirp sweep law.
 */
enum class ChirpMode : std::uint8_t {
    Linear, ///< Frequency increases linearly in time.
    Log,    ///< Frequency increases exponentially (constant ratio per unit time).
};

/**
 * @struct ChirpConfig
 * @brief Configuration for a sine chirp excitation.
 *
 * Defines `u(t) = amplitude * sin(phi(t))` with linear or logarithmic sweep.
 *
 * @tparam T Scalar type.
 */
template<typename T = double>
struct ChirpConfig {
    T         amplitude{T{1}};         ///< Signal amplitude (> 0)
    T         f_start_hz{T{1}};        ///< Start frequency in Hz (>= 0; > 0 for log mode)
    T         f_end_hz{T{10}};         ///< End frequency in Hz (>= 0; > 0 for log mode)
    T         duration_s{T{1}};        ///< Sweep duration in seconds (> 0)
    ChirpMode mode{ChirpMode::Linear}; ///< Sweep law (linear or logarithmic)

    /**
     * @brief Validate chirp configuration.
     * @return true if all parameters are finite and physically valid.
     */
    [[nodiscard]] constexpr bool valid() const {
        if (!damp::finite_positive(amplitude)) {
            return false;
        }
        if (!damp::finite_positive(duration_s)) {
            return false;
        }
        if (!damp::finite_non_negative(f_start_hz)) {
            return false;
        }
        if (!damp::finite_non_negative(f_end_hz)) {
            return false;
        }
        if (mode == ChirpMode::Log) {
            if (f_start_hz <= T{0} || f_end_hz <= T{0}) {
                return false;
            }
        }
        return true;
    }
};

/**
 * @struct ChirpResult
 * @brief Chirp design payload.
 *
 * Carries a validated chirp configuration into the runtime generator. Use
 * `.as<float>()` for embedded deployment after host-side design in double.
 *
 * @tparam T Scalar type.
 */
template<typename T = double>
struct ChirpResult {
    ChirpConfig<T> config{};       ///< Validated chirp configuration
    bool           success{false}; ///< true if `config.valid()`

    /**
     * @brief Convert design payload to another scalar type.
     * @tparam U Target scalar type.
     * @return ChirpResult\<U\> converted element-wise.
     */
    template<typename U>
    [[nodiscard]] constexpr ChirpResult<U> as() const {
        return ChirpResult<U>{
            ChirpConfig<U>{
                static_cast<U>(config.amplitude),
                static_cast<U>(config.f_start_hz),
                static_cast<U>(config.f_end_hz),
                static_cast<U>(config.duration_s),
                config.mode,
            },
            success,
        };
    }
};

/**
 * @brief Build a chirp design payload from a configuration.
 *
 * @param config Chirp configuration.
 * @return ChirpResult with `success = config.valid()`.
 */
template<typename T = double>
[[nodiscard]] constexpr ChirpResult<T>
chirp(const ChirpConfig<T>& config) {
    return ChirpResult<T>{config, config.valid()};
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Linear or logarithmic chirp runtime generator.
 *
 * Generates a bounded-time sinusoid with either linear or logarithmic
 * instantaneous frequency sweep. Supports absolute-time evaluation (`step(t)`)
 * and internal-sample-clock evaluation (`step()`).
 *
 * @tparam T Scalar type.
 */
template<typename T = float>
class Chirp {
public:
    constexpr Chirp() = default;

    /**
     * @brief Construct from a chirp design payload.
     * @param design Validated design payload.
     * @param Ts     Optional sample period for `step()` mode. If `Ts <= 0`,
     *               `step()` evaluates at a fixed internal time.
     */
    constexpr explicit Chirp(const design::ChirpResult<T>& design, T Ts = T{0})
        : config_(design.config), Ts_(Ts), valid_(design.success) {}

    /**
     * @brief Evaluate chirp at absolute time.
     * @param t Time in seconds.
     * @return Excitation value at clamped time `t in [0, duration]`.
     */
    [[nodiscard]] constexpr T step(T t) const {
        if (!valid_) {
            return T{0};
        }

        const T tc = clamp_time(t);
        const T phase = (config_.mode == design::ChirpMode::Linear)
                          ? linear_phase(tc)
                          : log_phase(tc);
        return config_.amplitude * damp::sin(phase);
    }

    /**
     * @brief Evaluate chirp at internal time and advance by `Ts`.
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
     * @brief Query whether the chirp duration has elapsed.
     * @param t Time in seconds.
     * @return true if generator is invalid or `t >= duration_s`.
     */
    [[nodiscard]] constexpr bool done(T t) const {
        if (!valid_) {
            return true;
        }
        return t >= config_.duration_s;
    }

    /**
     * @brief Query whether internal time has reached completion.
     * @return true if internal time has elapsed.
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

    /**
     * @brief Read current internal time.
     * @return Internal time in seconds.
     */
    [[nodiscard]] constexpr T time() const {
        return t_;
    }

private:
    [[nodiscard]] constexpr T clamp_time(T t) const {
        const T non_negative = damp::max(t, T{0});
        if (non_negative > config_.duration_s) {
            return config_.duration_s;
        }
        return non_negative;
    }

    [[nodiscard]] constexpr T linear_phase(T t) const {
        const T k = (config_.f_end_hz - config_.f_start_hz) / config_.duration_s;
        const T phase_cycles = (config_.f_start_hz * t) + (static_cast<T>(0.5) * k * t * t);
        return T{2} * damp::numbers::pi_v<T> * phase_cycles;
    }

    [[nodiscard]] constexpr T log_phase(T t) const {
        const T ratio = config_.f_end_hz / config_.f_start_hz;
        const T beta = damp::log(ratio) / config_.duration_s;
        if (damp::abs(beta) <= std::numeric_limits<T>::epsilon()) {
            return T{2} * damp::numbers::pi_v<T> * config_.f_start_hz * t;
        }
        const T phase_cycles = (config_.f_start_hz / beta) * (damp::exp(beta * t) - T{1});
        return T{2} * damp::numbers::pi_v<T> * phase_cycles;
    }

    design::ChirpConfig<T> config_{};

    T    Ts_{T{0}};
    T    t_{T{0}};
    bool valid_{false};
};

} // namespace damp
