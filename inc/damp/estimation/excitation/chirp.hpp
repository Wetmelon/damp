// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file chirp.hpp
 * @brief Chirp excitation — design payload and runtime generator
 */


#include "damp/backend.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix_traits.hpp" // eps

namespace damp {
namespace design {

/**
 * @brief Chirp sweep law.
 */
enum class ChirpMode : uint8_t {
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
    T              Ts{0};          ///< [s] Sample period of the internal step() clock (0 = absolute-time step(t) only)

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

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
            static_cast<U>(Ts),
        };
    }
};

/**
 * @brief Build a chirp design payload from a configuration.
 *
 * @param config Chirp configuration.
 * @param Ts     [s] Sample period of the runtime's step() clock (0 = absolute-time step(t) only; < 0 refused)
 * @return ChirpResult with `success` = `config.valid()`, a valid Ts, and (for Ts > 0) both
 *         sweep ends below the Nyquist frequency 1/(2·Ts).
 */
template<typename T = double>
[[nodiscard]] constexpr ChirpResult<T>
chirp(const ChirpConfig<T>& config, T Ts = T{0}) {
    const bool below_nyquist = (T{2} * config.f_start_hz * Ts < T{1}) && (T{2} * config.f_end_hz * Ts < T{1});
    return ChirpResult<T>{config, config.valid() && (Ts >= T{0}) && damp::isfinite(Ts) && below_nyquist, Ts};
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

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr Chirp(const design::ChirpResult<U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : Chirp(result.template as<T>()) {}

    /**
     * @brief Construct from a chirp design payload.
     * @param design Validated design payload.
     */
    constexpr Chirp(const design::ChirpResult<T>& design) // NOLINT(google-explicit-constructor): assign-from-design
        : config_(design.config), Ts_(design.Ts), valid_(design.success) {
        if (!design.success && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    /**
     * @brief Evaluate chirp at absolute time.
     * @param t Time in seconds.
     * @return Excitation value at `t` (a negative t evaluates at 0), or 0 once `t >= duration_s`.
     */
    [[nodiscard]] constexpr T step(T t) const {
        if (!valid_ || done(t)) {
            return T{0};
        }

        const T tc = damp::max(t, T{0});
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
        const T out = step(time());
        if (Ts_ > T{0} && !done()) {
            ++k_;
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
        return done(time());
    }

    /**
     * @brief Reset internal time to zero.
     */
    constexpr void reset() {
        k_ = 0;
    }

    /**
     * @brief Read current internal time.
     * @return Internal time in seconds, k·Ts after k samples (a count, not a running sum of Ts).
     */
    [[nodiscard]] constexpr T time() const {
        return static_cast<T>(k_) * Ts_;
    }

private:
    [[nodiscard]] constexpr T linear_phase(T t) const {
        const T k = (config_.f_end_hz - config_.f_start_hz) / config_.duration_s;
        const T phase_cycles = (config_.f_start_hz * t) + (static_cast<T>(0.5) * k * t * t);
        return T{2} * damp::numbers::pi_v<T> * phase_cycles;
    }

    [[nodiscard]] constexpr T log_phase(T t) const {
        // f(t) = f0·e^{βt}, β = ln(f1/f0)/T  ⇒  φ(t) = 2π·f0·(e^{βt} − 1)/β.
        // Test the dimensionless ln(f1/f0) (not β, which carries 1/s): when it is
        // within round-off of 0 the sweep is a constant tone.
        const T log_ratio = damp::log(config_.f_end_hz / config_.f_start_hz);
        if (damp::abs(log_ratio) <= eps<T>()) {
            return T{2} * damp::numbers::pi_v<T> * config_.f_start_hz * t;
        }
        const T beta = log_ratio / config_.duration_s;
        const T phase_cycles = (config_.f_start_hz / beta) * (damp::exp(beta * t) - T{1});
        return T{2} * damp::numbers::pi_v<T> * phase_cycles;
    }

    design::ChirpConfig<T> config_{};

    T      Ts_{T{0}};
    size_t k_{0}; ///< Samples emitted by step()
    bool   valid_{false};
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<typename T>
Chirp(const design::ChirpResult<T>&) -> Chirp<T>;

} // namespace damp
