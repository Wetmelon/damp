// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
// IWYU pragma: private, include "damp/controllers/pr.hpp"

/**
 * @file pr_resonator.hpp
 * @brief The resonant section PRController runs, and its coefficients at a fixed sample period, as
 *        members or as the constants of a `constexpr` design's Law
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

    /** @brief Convert across scalar precision, preserving the state. */
    template<typename U>
    constexpr explicit PRResonator(const PRResonator<U>& o)
        : x1_(static_cast<T>(o.x1_)),
          x2_(static_cast<T>(o.x2_)),
          y_(static_cast<T>(o.y_)),
          d_(static_cast<T>(o.d_)) {}

    /**
     * @brief Process one sample on the coefficients @p k: d += b̄₀·(x[n] − x[n−2]) − c₂·d − κ·y, then y += d
     *
     * The coefficients are the law's (a PRController's members, or the constants of a `constexpr`
     * design's Law), so a retune keeps the output and its increment.
     */
    constexpr T operator()(T x, const PRCoefficients<T>& k) {
        // The three small terms are summed first, so d takes one rounding at its own scale. Behind
        // a barrier: -ffast-math regroups d − c₂·d into (1 − c₂)·d, the direct form's ā₂ again
        // (GCC 14 with FMA: 1.3e-3 off the gain at ω₀ at 5 Hz and 100 kHz, ω_c = 1 rad/s).
        d_ += assoc_barrier(((k.b0 * (x - x2_)) - (k.c2 * d_)) - (k.kappa * y_));
        y_ += d_;
        x2_ = x1_;
        x1_ = x;
        return y_;
    }

    /** @brief Clear the inputs, the output and its increment. */
    constexpr void reset() {
        x1_ = x2_ = y_ = d_ = T{0};
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

    T x1_{0}; ///< x[n−1]
    T x2_{0}; ///< x[n−2]
    T y_{0};  ///< y[n−1]
    T d_{0};  ///< d[n−1] = y[n−1] − y[n−2]
};

/**
 * @brief The coefficients of every resonator of a PR bank: the h-th harmonic's resonant term,
 *        Kᵢ/h at h·ω₀, for each order h of @p orders, at the period @p Ts
 * @return The coefficients, or nullopt when Kp or Kᵢ is not finite, ω_c or the tracking time is
 *         below 0 or not finite, ω₀ or Ts is not finite and positive, an order is 0, a resonance
 *         is at or above Nyquist, or a coefficient overflows
 */
template<typename T, size_t NH>
[[nodiscard]] constexpr damp::optional<damp::array<PRCoefficients<T>, NH>>
pr_bank_coefficients(T Kp, T Ki, T w0, T wc, T tracking_time, T Ts, const damp::array<size_t, NH>& orders) {
    if (!damp::isfinite(Kp) || !damp::isfinite(Ki) || !damp::finite_positive(w0) || !damp::finite_non_negative(wc)
        || !damp::finite_non_negative(tracking_time) || !damp::finite_positive(Ts)) {
        return damp::nullopt;
    }
    damp::array<PRCoefficients<T>, NH> out{};
    for (size_t i = 0; i < NH; ++i) {
        if (orders[i] == 0) {
            return damp::nullopt; // the 0th harmonic (DC) is not a resonant target
        }
        const T order = static_cast<T>(orders[i]);
        if (!pr_resonance_ok(order * w0, Ts)) {
            return damp::nullopt;
        }
        out[i] = pr_coefficients(Ki / order, order * w0, wc, Ts);
        if (!out[i].finite()) {
            return damp::nullopt;
        }
    }
    return out;
}

/**
 * @brief A PRController's gains, harmonics, fixed sample period and the coefficients at it: the
 *        members a PRController holds, or the constants of a `constexpr` design's Law
 *        (ConstantPRGains)
 */
template<typename T, size_t NH>
struct PRGainSet {
    static constexpr bool constant = false; ///< Members of the runtime, settable

    T                                  Kp{T{0}};            ///< Proportional gain
    T                                  Ki{T{0}};            ///< Resonant gain at ω₀
    T                                  w0{T{0}};            ///< [rad/s] Resonant frequency, the fundamental's in a bank
    T                                  wc{T{0}};            ///< [rad/s] Resonant-term bandwidth
    T                                  tracking_time{T{0}}; ///< [s] Back-calculation tracking time (0: hold against a clamp)
    T                                  Ts{T{0}};            ///< [s] The tick's period
    T                                  unwind{T{0}};        ///< min(T_s/T_t, 1), 0 for T_t = 0 (hold)
    damp::array<size_t, NH>            orders{};            ///< The harmonic order of each resonator
    damp::array<PRCoefficients<T>, NH> res{};               ///< Each resonator's coefficients at T_s (inert until designed)
};

/**
 * @brief The gain set of Kp, Kᵢ, ω₀, ω_c and the tracking time at the period @p Ts for the
 *        harmonics @p orders, or nullopt as pr_bank_coefficients refuses them
 */
template<typename T, size_t NH>
[[nodiscard]] constexpr damp::optional<PRGainSet<T, NH>> pr_gain_set(T Kp, T Ki, T w0, T wc, T tracking_time, T Ts, const damp::array<size_t, NH>& orders) {
    const auto k = pr_bank_coefficients(Kp, Ki, w0, wc, tracking_time, Ts, orders);
    if (!k) {
        return damp::nullopt;
    }
    return PRGainSet<T, NH>{
        .Kp = Kp,
        .Ki = Ki,
        .w0 = w0,
        .wc = wc,
        .tracking_time = tracking_time,
        .Ts = Ts,
        .unwind = (tracking_time > T{0}) ? damp::min(Ts / tracking_time, T{1}) : T{0},
        .orders = orders,
        .res = *k,
    };
}

/**
 * @brief A `constexpr` PR design's gains as constants of the type: what a PRController reads in
 *        place of its gain members, in the runtime's scalar
 *
 * Each member is static, so the Law stores none of them: Kp, the tracking time's reciprocal,
 * T_s and every resonator's b̄₀, c₂ and κ sit in flash and fold into the tick. A result whose
 * gains are refused has no set, which the Law's static_assert reports.
 *
 * @tparam Res design::pr's result, `constexpr`
 * @tparam T   The runtime's scalar
 */
template<const auto& Res, typename T>
struct ConstantPRGains {
    static constexpr bool   constant = true;                                       ///< The gains are the type's: set_gains is absent
    static constexpr size_t kNH = damp::remove_cvref_t<decltype(Res)>::resonators; ///< The resonators

    static constexpr auto design = Res.template as<T>();                                                                                                                  ///< The design, in the runtime's scalar
    static constexpr auto opt = pr_gain_set(design.gains.Kp, design.gains.Ki, design.gains.w0, design.gains.wc, design.gains.tracking_time, design.Ts, design.harmonics); ///< The gains at T_s, if they are in their domain
    static constexpr bool ok = opt.has_value();                                                                                                                           ///< The gains are in their domain at T_s
    static constexpr auto set = ok ? *opt : PRGainSet<T, kNH>{};                                                                                                          ///< The gains at T_s

    static constexpr T    Kp = set.Kp;                       ///< Proportional gain
    static constexpr T    Ki = set.Ki;                       ///< Resonant gain at ω₀
    static constexpr T    w0 = set.w0;                       ///< [rad/s] Resonant frequency
    static constexpr T    wc = set.wc;                       ///< [rad/s] Resonant-term bandwidth
    static constexpr T    tracking_time = set.tracking_time; ///< [s] Back-calculation tracking time
    static constexpr T    Ts = set.Ts;                       ///< [s] The tick's period
    static constexpr T    unwind = set.unwind;               ///< min(T_s/T_t, 1), 0 for T_t = 0
    static constexpr auto orders = set.orders;               ///< The harmonic order of each resonator
    static constexpr auto res = set.res;                     ///< Each resonator's coefficients
};

} // namespace damp::detail
