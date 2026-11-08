// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file spectral.hpp
 * @brief Single-frequency spectral primitives for control: a generalized
 *        Goertzel single-bin DFT and a thin harmonic analyzer (THD / RMS).
 *
 * Deliberately *basic* — this is the controls layer, not a DSP library. The use
 * case is measuring a *known* frequency component (the grid fundamental and its
 * harmonics, a motor's electrical/mechanical ripple order) cheaply, in an ISR,
 * without a full FFT. For general spectral analysis use a dedicated DSP library.
 *
 * The Goertzel algorithm computes one DFT bin with a two-state recurrence — no
 * sample buffer, O(1) memory, O(N) per block. It is exact for *coherent*
 * sampling (an integer number of cycles of the target frequency in the block);
 * choose the block length N ≈ integer/(f·Ts) for the best results.
 *
 * All constexpr, allocation-free, float/double.
 */


#include "damp/backend.hpp"   // damp::array, damp::numbers
#include "damp/math/math.hpp" // damp::sin/cos/sqrt/atan2

namespace damp {

/**
 * @ingroup filters
 * @brief Generalized Goertzel single-bin DFT — amplitude/phase at one frequency.
 *
 * Feed samples one at a time with @ref push; when it returns true a block of N
 * samples is complete and @ref amplitude / @ref phase / @ref power are valid
 * until the next @ref push (which begins the next block). Detecting the same
 * frequency continuously is just a loop of `if (g.push(x)) { use g.amplitude() }`.
 *
 * The recurrence is Reinsch's form of s[n] = x[n] + 2cos ω·s[n−1] − s[n−2]: it carries
 * d[n] = s[n] ∓ s[n−1] and the coefficient ∓4 sin²(ω/2) / 4 cos²(ω/2), which keep the
 * bin frequency to full precision near ω = 0 and ω = π, where 2cos ω does not (in float
 * at f·Ts = 1e-4, 1 − cos ω is about three float spacings below 1).
 *
 * @see "Introduction to Numerical Analysis" (Stoer & Bulirsch, 2002), §2.3.3 (the algorithms of Goertzel and Reinsch)
 *
 * @tparam T Scalar type (float or double)
 */
template<typename T = float>
class Goertzel {
public:
    constexpr Goertzel() = default;

    /**
     * @brief Configure for @p freq_hz sampled every @p Ts over @p block_size samples.
     *
     * Requires @p Ts > 0 (finite), @p block_size ≥ 1, and 0 ≤ @p freq_hz < 1/(2·Ts) (Nyquist).
     * Invalid arguments leave the detector inert: @ref push never completes and
     * @ref amplitude returns 0 (no division by zero when N = 0).
     *
     * @param freq_hz    Target frequency [Hz]
     * @param Ts         Sample time [s]
     * @param block_size Samples per block, N (≥ 1)
     */
    constexpr Goertzel(T freq_hz, T Ts, size_t block_size)
        : n_(block_size) {
        const bool ok = (block_size >= 1) && (Ts > T{0}) && damp::isfinite(Ts) && (freq_hz >= T{0})
                     && (freq_hz * Ts < T{0.5});
        if (!ok) {
            n_ = 0;
            cos_ = T{1};
            sin_ = T{0};
            lambda_ = T{0};
            amp_scale_ = T{0};
            return;
        }
        const T omega = T{2} * damp::numbers::pi_v<T> * freq_hz * Ts;
        const auto [s, c] = damp::sincos(omega);
        cos_ = c;
        sin_ = s;
        near_dc_ = !(c < T{0});
        const auto [hs, hc] = damp::sincos(omega / T{2});
        lambda_ = near_dc_ ? (T{-4} * hs * hs) : (T{4} * hc * hc);
        // The recurrence ends on e^{jω(N−1)}·X(ω); undo that rotation once per block.
        const auto [rs, rc] = damp::sincos(omega * static_cast<T>(block_size - 1));
        rot_cos_ = rc;
        rot_sin_ = rs;
        // A DC bin carries the full level in |X|/N; any other bin splits it with −ω.
        amp_scale_ = ((freq_hz > T{0}) ? T{2} : T{1}) / static_cast<T>(block_size);
    }

    /**
     * @brief Feed one sample.
     * @return true when the block just completed (results valid until next push).
     */
    constexpr bool push(T x) {
        if (n_ == 0) {
            return false;
        }
        if (complete_) {
            // Auto-restart for the next block.
            s_ = T{0};
            d_ = T{0};
            count_ = 0;
            complete_ = false;
        }
        if (near_dc_) {
            // d[n] = s[n] − s[n−1]: d ← d + λ·s + x, s ← s + d, λ = −4 sin²(ω/2).
            d_ += (lambda_ * s_) + x;
            s_ += d_;
        } else {
            // d[n] = s[n] + s[n−1]: d ← λ·s − d + x, s ← d − s, λ = 4 cos²(ω/2).
            d_ = (lambda_ * s_) - d_ + x;
            s_ = d_ - s_;
        }
        if (++count_ >= n_) {
            // y = s[N−1] − e^{−jω}·s[N−2] = e^{jω(N−1)}·X(ω); X = y·e^{−jω(N−1)}.
            // Re y = s[N−1]·(1 ∓ cos ω) ± cos ω·d, without the s[N−1] − cos ω·s[N−2] cancellation.
            const T s2 = near_dc_ ? (s_ - d_) : (d_ - s_);
            const T yr = near_dc_ ? ((cos_ * d_) - (lambda_ * s_ / T{2})) : ((lambda_ * s_ / T{2}) - (cos_ * d_));
            const T yi = s2 * sin_;
            real_ = (yr * rot_cos_) + (yi * rot_sin_);
            imag_ = (yi * rot_cos_) - (yr * rot_sin_);
            complete_ = true;
            return true;
        }
        return false;
    }

    /** @brief True once a full block has been accumulated. */
    [[nodiscard]] constexpr bool complete() const { return complete_; }

    /** @brief Magnitude of the DFT bin, |X_k|. */
    [[nodiscard]] constexpr T magnitude() const { return damp::hypot(real_, imag_); }

    /** @brief Power, |X_k|². */
    [[nodiscard]] constexpr T power() const { return (real_ * real_) + (imag_ * imag_); }

    /**
     * @brief Peak amplitude of a sinusoid at the bin frequency (2·|X_k| / N; the DC level
     * |X₀| / N for a 0 Hz bin). Zero if N = 0.
     */
    [[nodiscard]] constexpr T amplitude() const { return amp_scale_ * magnitude(); }

    /**
     * @brief Phase of the bin [rad]: φ for x[n] = A·cos(ωn + φ) on a coherent block,
     * with n = 0 at the first sample of the block.
     */
    [[nodiscard]] constexpr T phase() const { return damp::atan2(imag_, real_); }

    /** @brief Real part of the DFT bin X(ω) = Σ x[n]·e^{−jωn} (valid after @ref complete). */
    [[nodiscard]] constexpr T real() const { return real_; }

    /** @brief Imaginary part of the DFT bin X(ω) (valid after @ref complete). */
    [[nodiscard]] constexpr T imag() const { return imag_; }

    /** @brief Discard the in-progress block and start over. */
    constexpr void reset() {
        s_ = T{0};
        d_ = T{0};
        count_ = 0;
        complete_ = false;
        real_ = T{0};
        imag_ = T{0};
    }

    [[nodiscard]] constexpr size_t block_size() const { return n_; }

private:
    T      lambda_{T{0}};  ///< −4 sin²(ω/2) near DC, 4 cos²(ω/2) near Nyquist
    bool   near_dc_{true}; ///< cos ω ≥ 0
    T      cos_{T{1}};
    T      sin_{T{0}};
    T      rot_cos_{T{1}};   ///< cos(ω(N−1)): end-of-block phase rotation
    T      rot_sin_{T{0}};   ///< sin(ω(N−1))
    T      amp_scale_{T{1}}; ///< 2/N (1/N for the DC bin); 0 when inert
    size_t n_{1};

    T      s_{T{0}}; ///< s[n]
    T      d_{T{0}}; ///< s[n] − s[n−1] near DC, s[n] + s[n−1] near Nyquist
    size_t count_{0};

    T    real_{T{0}};
    T    imag_{T{0}};
    bool complete_{false};
};

/**
 * @ingroup filters
 * @brief Harmonic analyzer — a Goertzel bank over a fundamental and K−1 harmonics.
 *
 * Measures the amplitude/phase of the fundamental f0 and its harmonics 2·f0 …
 * K·f0 in one pass, and reports total harmonic distortion and RMS. The headline
 * uses are grid voltage/current THD and motor torque-ripple order analysis. The
 * fundamental @ref rms is leakage-immune on a coherently-sampled block (unlike a
 * boxcar RMS, which sees the harmonics).
 *
 * @tparam K Number of harmonics tracked, including the fundamental (K ≥ 1)
 * @tparam T Scalar type (float or double)
 */
template<size_t K, typename T = float>
class HarmonicAnalyzer {
public:
    static_assert(K >= 1, "HarmonicAnalyzer needs at least the fundamental (K >= 1)");

    constexpr HarmonicAnalyzer() = default;

    /**
     * @brief Configure for fundamental @p f0_hz [Hz] sampled every @p Ts [s] over @p block_size samples.
     */
    constexpr HarmonicAnalyzer(T f0_hz, T Ts, size_t block_size) {
        for (size_t h = 0; h < K; ++h) {
            bins_[h] = Goertzel<T>(f0_hz * static_cast<T>(h + 1), Ts, block_size);
        }
    }

    /**
     * @brief Feed one sample to every harmonic bin.
     * @return true when the block just completed (results valid until next push).
     */
    constexpr bool push(T x) {
        // Every live bin shares the block length, so the fundamental's flag is the
        // block flag (a harmonic at/above Nyquist is inert and never completes).
        const bool done = bins_[0].push(x);
        for (size_t h = 1; h < K; ++h) {
            (void)bins_[h].push(x);
        }
        return done;
    }

    /** @brief Peak amplitude of harmonic @p harmonic (1 = fundamental, … K). Zero if out of range. */
    [[nodiscard]] constexpr T amplitude(size_t harmonic) const {
        if (harmonic < 1 || harmonic > K) {
            return T{0};
        }
        return bins_[harmonic - 1].amplitude();
    }

    /** @brief Phase of harmonic @p harmonic [rad] (1 = fundamental, … K). Zero if out of range. */
    [[nodiscard]] constexpr T phase(size_t harmonic) const {
        if (harmonic < 1 || harmonic > K) {
            return T{0};
        }
        return bins_[harmonic - 1].phase();
    }

    /** @brief RMS of the fundamental (A₁ / √2). */
    [[nodiscard]] constexpr T rms() const {
        return bins_[0].amplitude() / damp::numbers::sqrt2_v<T>;
    }

    /** @brief RMS of the harmonics above the fundamental, √(Σ_{h≥2} (Aₕ/√2)²). */
    [[nodiscard]] constexpr T total_harmonic_rms() const {
        T sum_sq = T{0};
        for (size_t h = 1; h < K; ++h) {
            const T a = bins_[h].amplitude();
            sum_sq += a * a;
        }
        return damp::sqrt(sum_sq) / damp::numbers::sqrt2_v<T>;
    }

    /** @brief Total harmonic distortion, √(Σ_{h≥2} Aₕ²) / A₁ (0 if the fundamental is 0). */
    [[nodiscard]] constexpr T thd() const {
        const T fundamental = bins_[0].amplitude();
        if (!(fundamental > T{0})) {
            return T{0};
        }
        return (total_harmonic_rms() * damp::numbers::sqrt2_v<T>) / fundamental;
    }

    [[nodiscard]] constexpr bool complete() const { return bins_[0].complete(); }

    constexpr void reset() {
        for (size_t h = 0; h < K; ++h) {
            bins_[h].reset();
        }
    }

private:
    damp::array<Goertzel<T>, K> bins_{};
};

} // namespace damp
