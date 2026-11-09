// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file stepped_sine.hpp
 * @brief Stepped-sine excitation — design payload and runtime generator
 */


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
template<size_t NFreq, typename T = double>
struct SteppedSineConfig {
    damp::array<T, NFreq> frequencies_hz{};    ///< Frequency table [Hz], each > 0
    T                     amplitude{T{1}};     ///< Tone amplitude (>= 0)
    size_t                cycles_per_freq{10}; ///< Periods spent at each frequency (>= 1)
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
        for (size_t i = 0; i < NFreq; ++i) {
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
template<size_t NFreq, typename T = double>
struct [[nodiscard]] SteppedSineResult {
    SteppedSineConfig<NFreq, T> config{};       ///< Validated configuration
    bool                        success{false}; ///< true if `config.valid()`
    T                           Ts{0};          ///< [s] Sample period of the internal step() clock (> 0)

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /**
     * @brief Convert design payload to another scalar type.
     * @tparam U Target scalar type.
     */
    template<typename U>
    [[nodiscard]] constexpr SteppedSineResult<NFreq, U> as() const {
        SteppedSineConfig<NFreq, U> cfg_u{};
        for (size_t i = 0; i < NFreq; ++i) {
            cfg_u.frequencies_hz[i] = static_cast<U>(config.frequencies_hz[i]);
        }
        cfg_u.amplitude = static_cast<U>(config.amplitude);
        cfg_u.cycles_per_freq = config.cycles_per_freq;
        cfg_u.phase_rad = static_cast<U>(config.phase_rad);
        return SteppedSineResult<NFreq, U>{cfg_u, success, static_cast<U>(Ts)};
    }
};

/**
 * @brief Build a stepped-sine design payload from a configuration.
 *
 * @tparam NFreq Number of frequency steps.
 * @tparam T     Scalar type.
 * @param config Stepped-sine configuration.
 * @param Ts     [s] Sample period of the runtime's step() clock (> 0: SteppedSine has no absolute-time step(t))
 * @return SteppedSineResult with `success` = `config.valid()`, a valid Ts, and every frequency
 *         below the Nyquist frequency 1/(2·Ts).
 */
template<size_t NFreq, typename T = double>
[[nodiscard]] constexpr SteppedSineResult<NFreq, T>
stepped_sine(const SteppedSineConfig<NFreq, T>& config, T Ts) {
    bool below_nyquist = true;
    for (size_t i = 0; i < NFreq; ++i) {
        below_nyquist = below_nyquist && (T{2} * config.frequencies_hz[i] * Ts < T{1});
    }
    return SteppedSineResult<NFreq, T>{config, config.valid() && (Ts > T{0}) && damp::isfinite(Ts) && below_nyquist, Ts};
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
template<size_t NFreq, typename T = float>
class SteppedSine {
public:
    constexpr SteppedSine() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr SteppedSine(const design::SteppedSineResult<NFreq, U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : SteppedSine(result.template as<T>()) {}

    /**
     * @brief Construct from a stepped-sine design payload.
     * @param design Validated design payload.
     */
    constexpr SteppedSine(const design::SteppedSineResult<NFreq, T>& design) // NOLINT(google-explicit-constructor): assign-from-design
        : config_(design.config), Ts_(design.Ts), valid_(design.success) {
        if (!design.success && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
        if (valid_) {
            start_frequency(0);
        }
    }

    /**
     * @brief Emit the current tone and advance internal time by `Ts`.
     *
     * The tone is evaluated at t = k·Ts from a sample count k, not a running sum of Ts.
     *
     * @return Excitation sample at the current frequency.
     */
    [[nodiscard]] constexpr T step() {
        if (!valid_ || done()) {
            return T{0};
        }

        const T f = config_.frequencies_hz[index_];
        const T t = static_cast<T>(k_) * Ts_;
        const T out = config_.amplitude * damp::sin((T{2} * damp::numbers::pi_v<T> * f * t) + config_.phase_rad);

        ++k_;
        if (f * static_cast<T>(k_) * Ts_ >= static_cast<T>(config_.cycles_per_freq)) {
            if (index_ + 1 < NFreq) {
                start_frequency(index_ + 1);
            } else {
                finished_ = true;
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

    /** @brief Index of the active frequency in the table. */
    [[nodiscard]] constexpr size_t frequency_index() const { return index_; }

    /** @brief Active frequency [Hz] (0 if invalid / finished). */
    [[nodiscard]] constexpr T frequency_hz() const {
        if (!valid_ || finished_ || index_ >= NFreq) {
            return T{0};
        }
        return config_.frequencies_hz[index_];
    }

    /** @brief Amplitude of the active tone. */
    [[nodiscard]] constexpr T amplitude() const { return config_.amplitude; }

private:
    constexpr void start_frequency(size_t i) {
        index_ = i;
        k_ = 0;
        finished_ = false;
    }

    design::SteppedSineConfig<NFreq, T> config_{};
    T                                   Ts_{T{0}};
    size_t                              k_{0}; ///< Samples emitted at the active frequency
    size_t                              index_{0};
    bool                                valid_{false};
    bool                                finished_{false};
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<size_t NFreq, typename T>
SteppedSine(const design::SteppedSineResult<NFreq, T>&) -> SteppedSine<NFreq, T>;

} // namespace damp
