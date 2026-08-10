// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file successive_compensator.hpp
 * @brief Successive residual commissioning: place notch or band-pass biquads.
 *
 * Greedy on-target procedure for servo / PE FRF features:
 *
 *     peak   (resonance)      → design::notch      (band-reject)
 *     valley (anti-resonance) → design::bandpass   (or peaking boost)
 *
 * Same fixed BiquadCascade; each section is a generic compensator, not
 * "notches only." Step-ring path still places peak→notch when FRF is unavailable.
 *
 * @see design::extract_frf_features, design::notch, design::bandpass
 * @see Bristow-Johnson biquad cookbook
 */

#include <cstddef>
#include <cstdint>

#include "damp/backend.hpp"
#include "damp/estimation/frequency_response.hpp"
#include "damp/filters/filters.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix_traits.hpp"

namespace damp {

/**
 * @brief Which biquad family to load for a feature.
 */
enum class CompensatorKind : std::uint8_t {
    Notch,    ///< Band-reject — kill a resonance peak
    Bandpass, ///< Band-pass — emphasize / shape an anti-resonance band
    Peaking,  ///< Peaking EQ cut/boost (optional valley map)
};

/**
 * @brief Configuration for successive compensator commissioning.
 *
 * @tparam T Scalar type.
 */
template<typename T = float>
struct SuccessiveCompensatorConfig {
    T               Ts{static_cast<T>(0.001)};      ///< Sample time [s] (> 0)
    T               notch_Q{T{5}};                  ///< RBJ notch Q (higher = narrower)
    T               bandpass_Q{T{5}};               ///< Band-pass Q for valleys
    T               peaking_Q{T{5}};                ///< Peaking Q when valley_map = Peaking
    T               peaking_gain_db{T{6}};          ///< Boost at valley for Peaking map [dB]
    T               min_peak{static_cast<T>(0.01)}; ///< Min residual peak |e| for step fit
    T               min_frequency_hz{T{1}};         ///< Reject rigid-body / drift
    T               max_frequency_hz{T{0}};         ///< 0 → ~0.45/Ts
    T               min_zeta{static_cast<T>(0.001)};
    T               max_zeta{static_cast<T>(0.5)};
    CompensatorKind valley_map{CompensatorKind::Bandpass}; ///< How valleys become biquads

    [[nodiscard]] constexpr bool valid() const {
        if (!(Ts > T{0}) || !(notch_Q > T{0}) || !(bandpass_Q > T{0}) || !(peaking_Q > T{0})) {
            return false;
        }
        if (!(min_peak >= T{0}) || !(min_frequency_hz >= T{0})) {
            return false;
        }
        if (!(min_zeta >= T{0}) || !(max_zeta > min_zeta)) {
            return false;
        }
        return damp::isfinite(Ts);
    }
};

namespace design {

/**
 * @brief Fit one underdamped second-order mode from a step / ring-down capture.
 *
 * Log-decrement on residual about steady state. Peak-only (not valleys).
 */
template<std::size_t N, typename T = float>
[[nodiscard]] constexpr ResonantMode<T> fit_second_order_step(
    const damp::array<T, N>& y,
    T                        Ts,
    T                        min_peak = T{0}
) {
    ResonantMode<T> mode{};
    mode.kind = FrfFeatureKind::Peak;
    if constexpr (N < 8) {
        return mode;
    }
    if (!(Ts > T{0})) {
        return mode;
    }

    const std::size_t n_ss = (N / 8) < 2 ? 2 : (N / 8);
    T                 sum = T{0};
    for (std::size_t i = N - n_ss; i < N; ++i) {
        sum += y[i];
    }
    const T y_ss = sum / static_cast<T>(n_ss);

    constexpr std::size_t              MaxPeaks = 8;
    damp::array<T, MaxPeaks>           peak_a{};
    damp::array<std::size_t, MaxPeaks> peak_k{};
    std::size_t                        n_peaks = 0;

    const std::size_t k0 = N / 10;
    for (std::size_t k = k0 + 1; k + 1 < N && n_peaks < MaxPeaks; ++k) {
        const T e = y[k] - y_ss;
        const T el = y[k - 1] - y_ss;
        const T er = y[k + 1] - y_ss;
        if (e > el && e >= er && e > min_peak) {
            peak_a[n_peaks] = e;
            peak_k[n_peaks] = k;
            ++n_peaks;
        }
    }
    if (n_peaks < 2) {
        n_peaks = 0;
        for (std::size_t k = k0 + 1; k + 1 < N && n_peaks < MaxPeaks; ++k) {
            const T e = damp::abs(y[k] - y_ss);
            const T el = damp::abs(y[k - 1] - y_ss);
            const T er = damp::abs(y[k + 1] - y_ss);
            if (e > el && e >= er && e > min_peak) {
                peak_a[n_peaks] = e;
                peak_k[n_peaks] = k;
                ++n_peaks;
            }
        }
    }
    if (n_peaks < 2) {
        return mode;
    }

    const T A1 = peak_a[0];
    const T A2 = peak_a[1];
    if (!(A1 > min_peak) || !(A2 > T{0})) {
        return mode;
    }

    const T Td = Ts * static_cast<T>(peak_k[1] - peak_k[0]);
    if (!(Td > T{0})) {
        return mode;
    }

    const T wd = (T{2} * damp::numbers::pi_v<T>) / Td;
    T       delta = T{0};
    if (A2 > default_tol<T>() && A1 > A2) {
        delta = damp::log(A1 / A2);
    }
    if (delta < T{0}) {
        delta = T{0};
    }

    const T two_pi = T{2} * damp::numbers::pi_v<T>;
    const T den = damp::sqrt((two_pi * two_pi) + (delta * delta));
    T       zeta = (den > default_tol<T>()) ? (delta / den) : T{0};
    if (zeta >= T{1}) {
        zeta = static_cast<T>(0.99);
    }
    const T wn = (zeta < T{1}) ? (wd / damp::sqrt(T{1} - (zeta * zeta))) : wd;

    mode.omega_n = wn;
    mode.frequency_hz = wn / (T{2} * damp::numbers::pi_v<T>);
    mode.zeta = zeta;
    mode.peak_gain = A1;
    mode.valid = mode.frequency_hz > T{0};
    return mode;
}

/**
 * @brief Map an FRF feature to RBJ biquad coefficients.
 */
template<typename T = float>
[[nodiscard]] constexpr design::SecondOrderCoeffs<T> compensator_from_feature(
    const ResonantMode<T>&                feature,
    const SuccessiveCompensatorConfig<T>& cfg
) {
    if (!feature.valid || !(feature.frequency_hz > T{0})) {
        return design::SecondOrderCoeffs<T>{T{1}, T{0}, T{0}, T{0}, T{0}};
    }
    if (feature.kind == FrfFeatureKind::Peak) {
        return design::notch(feature.frequency_hz, cfg.notch_Q, cfg.Ts);
    }
    // Valley
    if (cfg.valley_map == CompensatorKind::Peaking) {
        return design::peaking(feature.frequency_hz, cfg.peaking_Q, cfg.peaking_gain_db, cfg.Ts);
    }
    if (cfg.valley_map == CompensatorKind::Notch) {
        return design::notch(feature.frequency_hz, cfg.notch_Q, cfg.Ts);
    }
    return design::bandpass(feature.frequency_hz, cfg.bandpass_Q, cfg.Ts);
}

} // namespace design

/**
 * @brief Successive biquad compensator bank (notch and/or band-pass).
 *
 * @tparam NMax Maximum sections.
 * @tparam T    Scalar type.
 */
template<std::size_t NMax, typename T = float>
class SuccessiveCompensatorCommissioner {
public:
    static_assert(NMax >= 1, "SuccessiveCompensatorCommissioner needs NMax >= 1");

    constexpr SuccessiveCompensatorCommissioner() = default;

    constexpr explicit SuccessiveCompensatorCommissioner(const SuccessiveCompensatorConfig<T>& cfg)
        : cfg_(cfg), valid_(cfg.valid()) {
        cascade_.reset_identity();
    }

    [[nodiscard]] constexpr T compensate(T u) {
        if (!valid_ || n_ == 0) {
            return u;
        }
        return cascade_(u, n_);
    }

    template<std::size_t N>
    [[nodiscard]] constexpr bool assign_from_step_capture(const damp::array<T, N>& y) {
        auto mode = design::fit_second_order_step<N, T>(y, cfg_.Ts, cfg_.min_peak);
        mode.kind = FrfFeatureKind::Peak;
        return assign_feature(mode);
    }

    /**
     * @brief Install compensators for every peak (notch) and valley (band-pass/peaking)
     *        returned by @ref design::extract_frf_features, until full.
     */
    template<std::size_t MaxPeaks, std::size_t MaxValleys>
    [[nodiscard]] constexpr std::size_t assign_from_frf_features(
        const FrfFeatureResult<MaxPeaks, MaxValleys, T>& features
    ) {
        std::size_t added = 0;
        for (std::size_t i = 0; i < features.peak_count && n_ < NMax; ++i) {
            if (assign_feature(features.peaks[i])) {
                ++added;
            }
        }
        for (std::size_t i = 0; i < features.valley_count && n_ < NMax; ++i) {
            if (assign_feature(features.valleys[i])) {
                ++added;
            }
        }
        return added;
    }

    [[nodiscard]] constexpr bool assign_feature(const ResonantMode<T>& feature) {
        if (!valid_ || n_ >= NMax || !feature.valid) {
            return false;
        }
        T f_max = cfg_.max_frequency_hz;
        if (!(f_max > T{0})) {
            f_max = static_cast<T>(0.45) / cfg_.Ts;
        }
        if (feature.frequency_hz < cfg_.min_frequency_hz || feature.frequency_hz > f_max) {
            return false;
        }
        if (feature.kind == FrfFeatureKind::Peak) {
            if (feature.zeta < cfg_.min_zeta) {
                // still allow
            }
            if (feature.zeta > cfg_.max_zeta && feature.zeta > T{0}) {
                // oscillatory residual too damped for notch interest — still allow FRF peaks
            }
        }
        for (std::size_t i = 0; i < n_; ++i) {
            const T f0 = features_[i].frequency_hz;
            if (f0 > T{0} && damp::abs(feature.frequency_hz - f0) < static_cast<T>(0.15) * f0) {
                return false;
            }
        }

        const auto coeffs = design::compensator_from_feature(feature, cfg_);
        cascade_.set_section(n_, coeffs);
        cascade_.section(n_).reset();
        features_[n_] = feature;
        kinds_[n_] = (feature.kind == FrfFeatureKind::Peak) ? CompensatorKind::Notch : cfg_.valley_map;
        ++n_;
        return true;
    }

    [[nodiscard]] constexpr bool assign_mode(const ResonantMode<T>& mode) {
        ResonantMode<T> m = mode;
        if (m.kind != FrfFeatureKind::Valley) {
            m.kind = FrfFeatureKind::Peak;
        }
        return assign_feature(m);
    }

    [[nodiscard]] constexpr std::size_t count() const { return n_; }
    [[nodiscard]] constexpr bool        full() const { return n_ >= NMax; }
    [[nodiscard]] constexpr bool        empty() const { return n_ == 0; }

    [[nodiscard]] constexpr const ResonantMode<T>& feature(std::size_t i) const { return features_[i]; }
    [[nodiscard]] constexpr const ResonantMode<T>& mode(std::size_t i) const { return features_[i]; }
    [[nodiscard]] constexpr CompensatorKind        kind(std::size_t i) const { return kinds_[i]; }

    [[nodiscard]] constexpr BiquadCascade<NMax, T>&               cascade() { return cascade_; }
    [[nodiscard]] constexpr const BiquadCascade<NMax, T>&         cascade() const { return cascade_; }
    [[nodiscard]] constexpr const SuccessiveCompensatorConfig<T>& config() const { return cfg_; }

    constexpr void reset() {
        n_ = 0;
        cascade_.reset_identity();
        for (std::size_t i = 0; i < NMax; ++i) {
            features_[i] = ResonantMode<T>{};
            kinds_[i] = CompensatorKind::Notch;
        }
    }

private:
    SuccessiveCompensatorConfig<T>     cfg_{};
    BiquadCascade<NMax, T>             cascade_{};
    damp::array<ResonantMode<T>, NMax> features_{};
    damp::array<CompensatorKind, NMax> kinds_{};
    std::size_t                        n_{0};
    bool                               valid_{false};
};

} // namespace damp
