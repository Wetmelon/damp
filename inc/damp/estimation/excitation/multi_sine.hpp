// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file multi_sine.hpp
 * @brief Multi-sine excitation — design payload and runtime generator
 */

#include <cstddef>
#include <cstdint>

#include "damp/backend.hpp"
#include "damp/math/math.hpp"

namespace damp {
namespace design {

/**
 * @struct Tone
 * @brief One sinusoidal component in a multi-sine excitation.
 *
 * @tparam T Scalar type.
 */
template<typename T = double>
struct Tone {
    T amplitude{T{0}}; ///< Tone amplitude (>= 0)
    T freq_hz{T{0}};   ///< Tone frequency in Hz (>= 0)
    T phase_rad{T{0}}; ///< Tone phase offset in radians (finite)

    /**
     * @brief Convert tone to another scalar type.
     * @tparam U Target scalar type.
     * @return Tone\<U\> converted element-wise.
     */
    template<typename U>
    [[nodiscard]] constexpr Tone<U> as() const {
        return Tone<U>{
            static_cast<U>(amplitude),
            static_cast<U>(freq_hz),
            static_cast<U>(phase_rad),
        };
    }
};

/**
 * @struct MultiSineConfig
 * @brief Configuration for fixed-component multi-sine excitation.
 *
 * Defines:
 *     u(t) = sum_i Ai * sin(2*pi*fi*t + phii)
 *
 * At least one tone must have positive amplitude.
 *
 * @tparam NTones Number of tones.
 * @tparam T      Scalar type.
 */
template<std::size_t NTones, typename T = double>
struct MultiSineConfig {
    damp::array<Tone<T>, NTones> tones{}; ///< Fixed tone table

    /**
     * @brief Validate multi-sine configuration.
     * @return true if all tones are finite/non-negative and at least one tone has energy.
     */
    [[nodiscard]] constexpr bool valid() const {
        if constexpr (NTones == 0) {
            return false;
        }

        bool has_energy = false;
        for (std::size_t i = 0; i < NTones; ++i) {
            const auto& tone = tones[i];
            if (!damp::finite_non_negative(tone.amplitude)) {
                return false;
            }
            if (!damp::finite_non_negative(tone.freq_hz)) {
                return false;
            }
            if (!damp::isfinite(tone.phase_rad)) {
                return false;
            }
            if (tone.amplitude > T{0}) {
                has_energy = true;
            }
        }
        return has_energy;
    }
};

/**
 * @struct MultiSineResult
 * @brief Multi-sine design payload.
 *
 * Carries validated tone definitions into the runtime generator.
 *
 * @tparam NTones Number of tones.
 * @tparam T      Scalar type.
 */
template<std::size_t NTones, typename T = double>
struct MultiSineResult {
    MultiSineConfig<NTones, T> config{};       ///< Validated multi-sine configuration
    bool                       success{false}; ///< true if `config.valid()`

    /**
     * @brief Convert design payload to another scalar type.
     * @tparam U Target scalar type.
     * @return MultiSineResult<NTones, U> converted element-wise.
     */
    template<typename U>
    [[nodiscard]] constexpr MultiSineResult<NTones, U> as() const {
        MultiSineConfig<NTones, U> cfg_u{};
        for (std::size_t i = 0; i < NTones; ++i) {
            cfg_u.tones[i] = config.tones[i].template as<U>();
        }
        return MultiSineResult<NTones, U>{cfg_u, success};
    }
};

/**
 * @brief Build a multi-sine design payload from a configuration.
 *
 * @tparam NTones Number of tones.
 * @tparam T      Scalar type.
 * @param config  Multi-sine configuration.
 * @return MultiSineResult with `success = config.valid()`.
 */
template<std::size_t NTones, typename T = double>
[[nodiscard]] constexpr MultiSineResult<NTones, T>
multi_sine(const MultiSineConfig<NTones, T>& config) {
    return MultiSineResult<NTones, T>{config, config.valid()};
}

/**
 * @brief Assign Schroeder low-crest-factor phases to a multi-sine tone table.
 *
 * For tone index k = 0 .. N-1 (in table order), with 1-based harmonic index
 * n = k+1:
 *
 *     φ_n = −π · n · (n − 1) / N_tones
 *         = −π · (k+1) · k / N_tones
 *
 * (Schroeder 1970; first tone has phase 0). Amplitudes and frequencies are
 * left unchanged. Equal-amplitude harmonic sets get the classic low crest
 * factor; arbitrary frequency tables still benefit as a cheap phase heuristic.
 *
 * @tparam NTones Number of tones.
 * @tparam T      Scalar type.
 * @param config  Multi-sine configuration (phases overwritten).
 * @return Config with Schroeder phases written into each tone.
 *
 * @see Schroeder, "Synthesis of low-peak-factor signals and binary sequences
 *      with low autocorrelation" (IEEE Trans. Inf. Theory, 1970)
 */
template<std::size_t NTones, typename T = double>
[[nodiscard]] constexpr MultiSineConfig<NTones, T>
with_schroeder_phases(MultiSineConfig<NTones, T> config) {
    if constexpr (NTones == 0) {
        return config;
    }
    const T inv_n = T{1} / static_cast<T>(NTones);
    for (std::size_t k = 0; k < NTones; ++k) {
        const T kk = static_cast<T>(k);
        // n = k+1: φ = −π n (n−1) / N = −π (k+1) k / N
        config.tones[k].phase_rad = -damp::numbers::pi_v<T> * (kk + T{1}) * kk * inv_n;
    }
    return config;
}

/**
 * @brief Build a multi-sine design payload with Schroeder phases applied.
 *
 * Equivalent to `multi_sine(with_schroeder_phases(config))`.
 */
template<std::size_t NTones, typename T = double>
[[nodiscard]] constexpr MultiSineResult<NTones, T>
schroeder_multi_sine(const MultiSineConfig<NTones, T>& config) {
    return multi_sine(with_schroeder_phases(config));
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Sum-of-tones multi-sine runtime generator.
 *
 * Emits:
 *     u(t) = sum_i Ai * sin(2*pi*fi*t + phii)
 *
 * Unlike finite-duration signals, MultiSine is open-ended (`done()` is false
 * for valid configuration).
 *
 * @tparam NTones Number of tones.
 * @tparam T      Scalar type.
 */
template<std::size_t NTones, typename T = float>
class MultiSine {
public:
    constexpr MultiSine() = default;

    /**
     * @brief Construct from a multi-sine design payload.
     * @param design Validated design payload.
     * @param Ts     Optional sample period for internal `step()` mode.
     */
    constexpr explicit MultiSine(const design::MultiSineResult<NTones, T>& design, T Ts = T{0})
        : config_(design.config), Ts_(Ts), valid_(design.success) {}

    /**
     * @brief Evaluate multi-sine at absolute time.
     * @param t Time in seconds.
     * @return Sum of all configured tones at time `t`.
     */
    [[nodiscard]] constexpr T step(T t) const {
        if (!valid_) {
            return T{0};
        }

        const T tc = damp::max(t, T{0});
        T       y = T{0};
        for (std::size_t i = 0; i < NTones; ++i) {
            const auto& tone = config_.tones[i];
            const T     omega_t = (T{2} * damp::numbers::pi_v<T> * tone.freq_hz * tc) + tone.phase_rad;
            y += tone.amplitude * damp::sin(omega_t);
        }
        return y;
    }

    /**
     * @brief Evaluate multi-sine at internal time and advance by `Ts`.
     * @return Excitation value at current internal time.
     */
    [[nodiscard]] constexpr T step() {
        const T out = step(t_);
        if (Ts_ > T{0}) {
            t_ += Ts_;
        }
        return out;
    }

    /**
     * @brief Query completion at absolute time.
     * @param t Time in seconds (ignored).
     * @return true only when configuration is invalid.
     */
    [[nodiscard]] constexpr bool done(T t) const {
        (void)t;
        return !valid_;
    }

    /**
     * @brief Query completion for internal-time mode.
     * @return true only when configuration is invalid.
     */
    [[nodiscard]] constexpr bool done() const {
        return !valid_;
    }

    /**
     * @brief Reset internal time to zero.
     */
    constexpr void reset() {
        t_ = T{0};
    }

private:
    design::MultiSineConfig<NTones, T> config_{};
    T                                  Ts_{T{0}};
    T                                  t_{T{0}};
    bool                               valid_{false};
};

} // namespace damp
