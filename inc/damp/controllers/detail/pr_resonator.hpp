// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
// IWYU pragma: private, include "damp/controllers/pr.hpp"

/**
 * @file pr_resonator.hpp
 * @brief The resonant section PRController runs, and its coefficients at a fixed sample period
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/anti_windup.hpp"
#include "damp/math/math.hpp"

namespace damp::detail {

/**
 * @brief Pre-warped resonant frequency (2/Ts)·tan(ω₀Ts/2) that PRController discretizes at
 *
 * The bilinear map s → (2/Ts)(z−1)/(z+1) of a resonator rebuilt at this frequency places the
 * digital peak exactly on ω₀. Returns ω₀ unchanged when Ts is not positive or when ω₀ is at or
 * above Nyquist (tan has no finite value there).
 */
template<typename T>
[[nodiscard]] constexpr T pr_prewarped_w0(T w0, T Ts) {
    const T half = w0 * Ts / T{2};
    if (!damp::finite_positive(Ts) || !(half < damp::numbers::pi_v<T> / T{2})) {
        return w0;
    }
    return (T{2} / Ts) * damp::tan(half);
}

/**
 * @brief β = tan(ω₀Ts/2): the pre-warped resonant frequency of pr_prewarped_w0 over the Tustin factor 2/Ts
 *
 * What PRController builds its coefficients from: the frequency of @ref pr_prewarped_w0 without
 * the factor 2/Ts, whose square would swamp it in single precision. At or above Nyquist it is
 * ω₀Ts/2 itself, the frequency left unwarped there.
 */
template<typename T>
[[nodiscard]] constexpr T pr_tan_half_angle(T w0, T Ts) {
    const T half = (w0 * Ts) / T{2};
    return (half < (damp::numbers::pi_v<T> / T{2})) ? damp::tan(half) : half;
}

/** @brief The coefficients of a PRResonator: b̄₀, c₂ = 1 − ā₂ and κ = 1 + ā₁ + ā₂. */
template<typename T>
struct PRCoefficients {
    T b0{T{0}};    ///< b̄₀, the gain on x[n] − x[n−2]
    T c2{T{1}};    ///< c₂ = 1 − ā₂ (≥ 0; 0 = undamped)
    T kappa{T{1}}; ///< κ = 1 + ā₁ + ā₂ (> 0)

    /** @brief Every coefficient finite. */
    [[nodiscard]] constexpr bool finite() const { return damp::isfinite(b0) && damp::isfinite(c2) && damp::isfinite(kappa); }
};

/**
 * @brief The coefficients of the resonant term R(s) of gain @p Ki at @p w0 with bandwidth @p wc,
 *        sampled at @p Ts, each a ratio of positive terms
 *
 * Tustin of R(s), s = k(z−1)/(z+1) with k = 2/Ts, at the pre-warped frequency
 * @f$\omega_0' = k\beta@f$, @f$\beta = \tan(\omega_0 T_s/2)@f$, so that the digital peak sits
 * exactly on ω₀ (design::PRResult::to_discrete_ss uses the same frequency). Divided through by
 * k², with @f$\alpha = \omega_c T_s/2@f$, the denominator is
 * @f$(z-1)^2 + 2\alpha(z^2-1) + \beta^2(z+1)^2@f$, and normalizing by
 * @f$\mathrm{den} = 1 + 2\alpha + \beta^2@f$ gives the coefficients of @ref PRResonator:
 * @f[
 *   \kappa = \frac{4\beta^2}{\mathrm{den}}, \qquad c_2 = \frac{4\alpha}{\mathrm{den}}, \qquad
 *   \bar b_0 = \frac{2 K_i \alpha}{\mathrm{den}} = \frac{K_i c_2}{2}
 *   \quad\Bigl(\text{ideal, } \alpha = 0:\ \bar b_0 = \frac{K_i T_s}{2\,\mathrm{den}}\Bigr).
 * @f]
 * No k² is formed and nothing cancels: in float the sum k² + 2ω_c·k + ω₀'² drops ω₀'² once
 * it is below the ulp of k² (39 against 128 at 1 Hz and 20 kHz), and −2 + κ + c₂ cannot be
 * stored either (see @ref PRResonator). κ and c₂ are within a few ulps of their values. The
 * poles have |z|² = 1 − c₂ and are inside the unit circle while κ + 2c₂ < 4, which holds by
 * 4/den. A den that is not positive and finite gives the inert section (b̄₀ = 0, c₂ = κ = 1),
 * which finite() does not catch: the caller checks the design first (pr_resonance_ok).
 */
template<typename T>
[[nodiscard]] constexpr PRCoefficients<T> pr_coefficients(T Ki, T w0, T wc, T Ts) {
    const T alpha = (wc * Ts) / T{2};
    const T beta = pr_tan_half_angle(w0, Ts);
    const T beta2 = beta * beta;
    const T den = T{1} + (T{2} * alpha) + beta2;
    if (!damp::finite_positive(den)) {
        return PRCoefficients<T>{};
    }
    const T c2 = T{4} * (alpha / den);
    const T kappa = T{4} * (beta2 / den);
    const T b0 = (wc > T{0}) ? (Ki * (c2 / T{2})) : (Ki * (Ts / (T{2} * den)));
    return PRCoefficients<T>{.b0 = b0, .c2 = c2, .kappa = kappa};
}

/** @brief Whether a resonance at @p w [rad/s] can be placed at @p Ts: finite, positive and below Nyquist π/Ts (the bilinear map folds one at or above it onto another frequency). */
template<typename T>
[[nodiscard]] constexpr bool pr_resonance_ok(T w, T Ts) {
    return damp::finite_positive(w) && damp::finite_positive(Ts) && (w * Ts < damp::numbers::pi_v<T>);
}

/**
 * @brief Resonant section of PRController, run on its output y and the output's increment d
 *
 * The Tustin map of the resonant term is a second-order section with zeros at z = ±1,
 * @f[
 *   \frac{Y(z)}{X(z)} = \frac{\bar b_0\,(1 - z^{-2})}{1 + \bar a_1 z^{-1} + \bar a_2 z^{-2}},
 *   \qquad \bar a_1 = \kappa + c_2 - 2, \quad \bar a_2 = 1 - c_2 .
 * @f]
 * It is held by κ = 1 + ā₁ + ā₂ and c₂ = 1 − ā₂ and run on the output and its increment
 * @f$d[n] = y[n] - y[n-1]@f$:
 * @f[
 *   d[n] = d[n-1] + \bigl(\bar b_0\,(x[n] - x[n-2]) - c_2\, d[n-1] - \kappa\, y[n-1]\bigr),
 *   \qquad y[n] = y[n-1] + d[n] .
 * @f]
 * Substituting d gives the direct form
 * @f$y[n] = -\bar a_1 y[n-1] - \bar a_2 y[n-2] + \bar b_0 (x[n] - x[n-2])@f$ term for term, so
 * the transfer function and the meaning of the state (the last output, and through d the one
 * before it) are those of a Direct Form I section; only the rounding differs. It is the
 * section written around z = 1, an accumulator in place of each delay.
 *
 * What it guards against. With θ = ω₀Ts and a narrow resonance, the poles sit next to z = 1:
 * ā₁ ≈ −2 + θ² + 2ω_cTs and ā₂ ≈ 1 − 2ω_cTs. Single precision stores ā₁ to 1.2e-7, and θ² is
 * 9.9e-8 at 1 Hz and 20 kHz, so the direct form cannot hold the resonance: a float Direct Form I
 * section at 1 Hz and 20 kHz with ω_c = 5 rad/s ran to Inf in 36 s of a 0.01 sine at ω₀, and at
 * 50 Hz and 100 kHz with ω_c = 1 rad/s its gain at ω₀ was 0.78 of the design's. Here κ ≈ θ² and
 * c₂ ≈ 2ω_cTs are stored with their own relative precision, and each product in the recursion
 * is small next to the state it corrects, so its rounding is too.
 *
 * Limits. The resonance is placed to a few ulps of ω₀ and its peak is ω_c/(1 + tan²(θ/2)) wide
 * (the bilinear map narrows it toward Nyquist), so the gain and phase at ω₀ are the design's to
 * a few ε·Q, Q = ω₀(1 + tan²(θ/2))/(2ω_c): a float resonance with Q above about 1e6 is not
 * resolved. Roundoff in the recursion is about ε/√c₂ of the output, and the ideal form
 * (c₂ = 0), which is marginally stable by design, accumulates it. Next to Nyquist the poles
 * approach z = −1, where this form has the direct form's cancellation: the stability margin
 * 4 − κ − 2c₂ is 4/(1 + ω_cTs + tan²(θ/2)) and falls below rounding within about fs·√ε/(2π)
 * of Nyquist (5.5e-5·fs in float).
 *
 * @see "New Recursive Digital Filter Structures Having Very Low Sensitivity and Roundoff Noise"
 *      (Agarwal & Burrus, 1975), IEEE Trans. Circuits and Systems 22(12): second-order sections
 *      for poles near z = 1, whose coefficient sensitivity and roundoff noise do not grow with
 *      the sampling rate
 * @see "Improved Finite Word Length Characteristics in Digital Control Using Delta Operators"
 *      (Middleton & Goodwin, 1986), IEEE TAC 31(11): the delta operator, the increment in place
 *      of the shift, for better coefficient representation and less roundoff noise
 *
 * @tparam T Scalar type
 */
template<typename T>
class PRResonator {
public:
    constexpr PRResonator() = default;

    /** @brief Convert across scalar precision, preserving coefficients and state. */
    template<typename U>
    constexpr explicit PRResonator(const PRResonator<U>& o)
        : b0_(static_cast<T>(o.b0_)),
          c2_(static_cast<T>(o.c2_)),
          kappa_(static_cast<T>(o.kappa_)),
          x1_(static_cast<T>(o.x1_)),
          x2_(static_cast<T>(o.x2_)),
          y_(static_cast<T>(o.y_)),
          d_(static_cast<T>(o.d_)) {}

    /** @brief Process one sample: d += b̄₀·(x[n] − x[n−2]) − c₂·d − κ·y, then y += d. */
    constexpr T operator()(T x) {
        // The three small terms are summed first, so d takes one rounding at its own scale. Behind
        // a barrier: -ffast-math regroups d − c₂·d into (1 − c₂)·d, the direct form's ā₂ again
        // (GCC 14 with FMA: 1.3e-3 off the gain at ω₀ at 5 Hz and 100 kHz, ω_c = 1 rad/s).
        d_ += assoc_barrier(((b0_ * (x - x2_)) - (c2_ * d_)) - (kappa_ * y_));
        y_ += d_;
        x2_ = x1_;
        x1_ = x;
        return y_;
    }

    /** @brief Clear the inputs, the output and its increment; the coefficients are kept. */
    constexpr void reset() {
        x1_ = x2_ = y_ = d_ = T{0};
    }

    /** @brief Replace the coefficients in place, keeping the state, so a retune keeps the output and its increment. */
    constexpr void set_coefficients(const PRCoefficients<T>& k) {
        b0_ = k.b0;
        c2_ = k.c2;
        kappa_ = k.kappa;
    }

    /** @brief Most recent output y[n−1]. */
    [[nodiscard]] constexpr T last_output() const { return y_; }

    /** @brief Add @p dy to the most recent output y[n−1] and leave y[n−2] in place: the increment moves with it. */
    constexpr void shift_last_output(T dy) {
        y_ += dy;
        d_ += dy;
    }

    /**
     * @brief Take back what the last sample added to the output, d[n−1], as far as a clamp cut the
     *        command off (@p excess = u_sat − u_asked), and no further (detail::hold_against_clamp)
     * @return How far the output moved
     */
    [[nodiscard]] constexpr T hold_against_clamp(T excess) {
        const T held = y_;
        detail::hold_against_clamp(y_, d_, excess);
        return y_ - held;
    }

private:
    template<typename>
    friend class PRResonator;

    T b0_{0};    ///< b̄₀; with c₂ = κ = 1 the section is inert (y = 0) until designed
    T c2_{1};    ///< c₂ = 1 − ā₂
    T kappa_{1}; ///< κ = 1 + ā₁ + ā₂
    T x1_{0};    ///< x[n−1]
    T x2_{0};    ///< x[n−2]
    T y_{0};     ///< y[n−1]
    T d_{0};     ///< d[n−1] = y[n−1] − y[n−2]
};

} // namespace damp::detail
