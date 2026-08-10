// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file stepped_sine.hpp
 * @brief Stepped-sine excitation — design payload and runtime generator
 */

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/math/math.hpp"

namespace damp {
namespace design {

/**
 * @struct SteppedSineConfig
 * @brief Configuration for single-tone stepped-sine excitation.
 *
 * Visits each frequency in @c frequencies_hz for @c cycles_per_freq periods
 * of a sinusoid at amplitude @c amplitude. Used as the stimulus for sequential
 * FRF estimation (FrequencyResponseEstimator).
 *
 * @tparam NFreq Number of frequency steps.
 * @tparam T     Scalar type.
 */
template<std::size_t NFreq, typename T = double>
struct SteppedSineConfig {
    damp::array<T, NFreq> frequencies_hz{};    ///< Frequency table [Hz], each > 0
    T                     amplitude{T{1}};     ///< Tone amplitude (>= 0)
    std::size_t           cycles_per_freq{10}; ///< Periods spent at each frequency (>= 1)
    T                     phase_rad{T{0}};     ///< Phase offset [rad] (finite)

    /**
     * @brief Validate stepped-sine configuration.
     * @return true if amplitude/freqs/cycles are usable.
     */
    [[nodiscard]] constexpr bool valid() const {
        if constexpr (NFreq == 0) {
            return false;
        }
        if (!damp::finite_non_negative(amplitude)) {
            return false;
        }
        if (cycles_per_freq == 0) {
            return false;
        }
        if (!damp::isfinite(phase_rad)) {
            return false;
        }
        for (std::size_t i = 0; i < NFreq; ++i) {
            if (!damp::finite_positive(frequencies_hz[i])) {
                return false;
            }
        }
        return true;
    }
};

/**
 * @struct SteppedSineResult
 * @brief Stepped-sine design payload.
 *
 * @tparam NFreq Number of frequency steps.
 * @tparam T     Scalar type.
 */
template<std::size_t NFreq, typename T = double>
struct SteppedSineResult {
    SteppedSineConfig<NFreq, T> config{};       ///< Validated configuration
    bool                        success{false}; ///< true if `config.valid()`

    /**
     * @brief Convert design payload to another scalar type.
     * @tparam U Target scalar type.
     */
    template<typename U>
    [[nodiscard]] constexpr SteppedSineResult<NFreq, U> as() const {
        SteppedSineConfig<NFreq, U> cfg_u{};
        for (std::size_t i = 0; i < NFreq; ++i) {
            cfg_u.frequencies_hz[i] = static_cast<U>(config.frequencies_hz[i]);
        }
        cfg_u.amplitude = static_cast<U>(config.amplitude);
        cfg_u.cycles_per_freq = config.cycles_per_freq;
        cfg_u.phase_rad = static_cast<U>(config.phase_rad);
        return SteppedSineResult<NFreq, U>{cfg_u, success};
    }
};

/**
 * @brief Build a stepped-sine design payload from a configuration.
 *
 * @tparam NFreq Number of frequency steps.
 * @tparam T     Scalar type.
 * @param config Stepped-sine configuration.
 * @return SteppedSineResult with `success = config.valid()`.
 */
template<std::size_t NFreq, typename T = double>
[[nodiscard]] constexpr SteppedSineResult<NFreq, T>
stepped_sine(const SteppedSineConfig<NFreq, T>& config) {
    return SteppedSineResult<NFreq, T>{config, config.valid()};
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Stepped-sine excitation — one pure tone at a time across a frequency table.
 *
 * Emits:
 *     u(t) = A * sin(2*pi*f_i*t + phi)
 *
 * for `cycles_per_freq` periods at each f_i, then advances to the next table
 * entry. Finite duration: @ref done is true after the last frequency finishes.
 * Pair with FrequencyResponseEstimator for on-target FRF measurement.
 *
 * @tparam NFreq Number of frequency steps.
 * @tparam T     Scalar type.
 */
template<std::size_t NFreq, typename T = float>
class SteppedSine {
public:
    constexpr SteppedSine() = default;

    /**
     * @brief Construct from a stepped-sine design payload.
     * @param design Validated design payload.
     * @param Ts     Sample period [s] for internal `step()` mode (must be > 0 to advance).
     */
    constexpr explicit SteppedSine(const design::SteppedSineResult<NFreq, T>& design, T Ts = T{0})
        : config_(design.config), Ts_(Ts), valid_(design.success) {
        if (valid_) {
            start_frequency(0);
        }
    }

    /**
     * @brief Emit the current tone and advance internal time by `Ts`.
     * @return Excitation sample at the current frequency.
     */
    [[nodiscard]] constexpr T step() {
        if (!valid_ || done()) {
            return T{0};
        }

        const T f = config_.frequencies_hz[index_];
        const T out = config_.amplitude * damp::sin((T{2} * damp::numbers::pi_v<T> * f * t_local_) + config_.phase_rad);

        if (Ts_ > T{0}) {
            t_local_ += Ts_;
            phase_accum_ += T{2} * damp::numbers::pi_v<T> * f * Ts_;
            const T target = T{2} * damp::numbers::pi_v<T> * static_cast<T>(config_.cycles_per_freq);
            if (phase_accum_ >= target) {
                if (index_ + 1 < NFreq) {
                    start_frequency(index_ + 1);
                } else {
                    finished_ = true;
                }
            }
        }
        return out;
    }

    /**
     * @brief Query whether the full frequency table has been visited.
     */
    [[nodiscard]] constexpr bool done() const {
        return !valid_ || finished_;
    }

    /**
     * @brief Reset to the first frequency.
     */
    constexpr void reset() {
        finished_ = false;
        if (valid_) {
            start_frequency(0);
        }
    }

    /// Index of the active frequency in the table.
    [[nodiscard]] constexpr std::size_t frequency_index() const { return index_; }

    /// Active frequency [Hz] (0 if invalid / finished).
    [[nodiscard]] constexpr T frequency_hz() const {
        if (!valid_ || finished_ || index_ >= NFreq) {
            return T{0};
        }
        return config_.frequencies_hz[index_];
    }

    /// Amplitude of the active tone.
    [[nodiscard]] constexpr T amplitude() const { return config_.amplitude; }

private:
    constexpr void start_frequency(std::size_t i) {
        index_ = i;
        t_local_ = T{0};
        phase_accum_ = T{0};
        finished_ = false;
    }

    design::SteppedSineConfig<NFreq, T> config_{};
    T                                   Ts_{T{0}};
    T                                   t_local_{T{0}};
    T                                   phase_accum_{T{0}};
    std::size_t                         index_{0};
    bool                                valid_{false};
    bool                                finished_{false};
};

} // namespace damp
