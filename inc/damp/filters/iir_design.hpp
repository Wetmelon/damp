// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file iir_design.hpp
 * @brief IIR coefficient design (analog TFs, Tustin/RBJ/Butterworth) and to_coeffs
 */


#include <cstddef>

#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"


namespace damp {
namespace design {

/**
 * @brief DSP coefficients for first-order IIR filter
 */
template<typename T = float>
struct FirstOrderCoeffs {
    T b0{0}, b1{0}, a1{0};

    /**
     * @brief Convert coefficients to different scalar type
     * @tparam U Target scalar type
     * @return FirstOrderCoeffs\<U\> with converted coefficients
     */
    template<typename U>
    [[nodiscard]] constexpr FirstOrderCoeffs<U> as() const {
        return {static_cast<U>(b0), static_cast<U>(b1), static_cast<U>(a1)};
    }

    /**
     * @brief Discrete state-space realization (controllable canonical form).
     *
     * Returns the z-domain biquad H(z) = (b0 + b1·z⁻¹)/(1 + a1·z⁻¹) as a discrete
     * StateSpace<1,1,1> (Ts > 0) so the design drops into the discrete analysis /
     * simulation tooling. Inverse of `to_coeffs`.
     */
    [[nodiscard]] constexpr StateSpace<1, 1, 1, T> to_state_space(T Ts) const {
        return StateSpace<1, 1, 1, T>{
            .A = Matrix<1, 1, T>{{-a1}},
            .B = Matrix<1, 1, T>{{T{1}}},
            .C = Matrix<1, 1, T>{{b1 - (a1 * b0)}},
            .D = Matrix<1, 1, T>{{b0}},
            .Ts = Ts,
        };
    }
};

/**
 * @brief DSP coefficients for second-order IIR filter
 */
template<typename T = float>
struct SecondOrderCoeffs {
    T b0{0}, b1{0}, b2{0}, a1{0}, a2{0};

    /**
     * @brief Convert coefficients to different scalar type
     * @tparam U Target scalar type
     * @return SecondOrderCoeffs\<U\> with converted coefficients
     */
    template<typename U>
    [[nodiscard]] constexpr SecondOrderCoeffs<U> as() const {
        return {static_cast<U>(b0), static_cast<U>(b1), static_cast<U>(b2), static_cast<U>(a1), static_cast<U>(a2)};
    }

    /**
     * @brief Discrete state-space realization (controllable canonical form).
     *
     * Returns the z-domain biquad H(z) = (b0 + b1·z⁻¹ + b2·z⁻²)/(1 + a1·z⁻¹ +
     * a2·z⁻²) as a discrete StateSpace<2,1,1> (Ts > 0) so any biquad design
     * (notch/bandpass/peaking/shelf/…) drops into the discrete analysis /
     * simulation tooling. Inverse of `to_coeffs`.
     */
    [[nodiscard]] constexpr StateSpace<2, 1, 1, T> to_state_space(T Ts) const {
        return StateSpace<2, 1, 1, T>{
            .A = Matrix<2, 2, T>{{-a1, -a2}, {T{1}, T{0}}},
            .B = Matrix<2, 1, T>{{T{1}}, {T{0}}},
            .C = Matrix<1, 2, T>{{b1 - (a1 * b0), b2 - (a2 * b0)}},
            .D = Matrix<1, 1, T>{{b0}},
            .Ts = Ts,
        };
    }
};

/**
 * @defgroup filters Filter Design
 * @brief IIR/FIR coefficient design and related discrete filter helpers
 *
 * Functions for designing common filters used in control systems.
 * Filters can be designed in continuous-time or discrete-time, with
 * compile-time and runtime support through constexpr evaluation. Runtime
 * filter objects live in the sibling filter headers in this folder.
 */

namespace detail {
/**
 * @brief True when @p f [Hz] is strictly inside (0, Nyquist) for sample time @p Ts.
 *
 * Nyquist = 1/(2·Ts). Invalid designs return zero coefficients (no division by
 * zero / NaN poles).
 */
template<typename T>
[[nodiscard]] constexpr bool valid_cutoff(T f, T Ts) {
    return (Ts > T{0}) && (f > T{0}) && ((T{2} * f * Ts) < T{1});
}

template<typename T>
[[nodiscard]] constexpr bool valid_Q(T Q) {
    return Q > T{0};
}
} // namespace detail

/**
 * @brief First-order low-pass filter design
 *
 * Designs a discrete-time first-order low-pass filter with given cutoff frequency.
 * Returns DSP coefficients ready for runtime use.
 *
 * @param fc Cutoff frequency [Hz] (must satisfy 0 < fc < 1/(2·Ts))
 * @param Ts Sample time [s] (must be > 0)
 *  * @return FirstOrderCoeffs\<T\> DSP coefficients; all zero if @p fc / @p Ts invalid
 */
template<typename T = float>
[[nodiscard]] constexpr FirstOrderCoeffs<T> lowpass_1st(T fc, T Ts) {
    if (!detail::valid_cutoff(fc, Ts)) {
        return FirstOrderCoeffs<T>{};
    }
    // Bilinear / Tustin: s ← k (1 − z⁻¹)/(1 + z⁻¹), k = 2/Ts.
    // Pole gives a1; numerator is forced to unit DC by construction so
    // (b0 + b1) = (1 + a1) holds under -ffast-math reassociation (same pattern
    // as lowpass_2nd).
    const T omega_c = T{2} * damp::numbers::pi_v<T> * fc;
    const T k = T{2} / Ts;
    const T denom = omega_c + k;

    const T a1 = (omega_c - k) / denom;
    const T dc_sum = T{1} + a1; // = b0 + b1 for H(1) = 1
    return FirstOrderCoeffs<T>{
        .b0 = dc_sum / T{2},
        .b1 = dc_sum / T{2},
        .a1 = a1,
    };
}

/**
 * @brief First-order high-pass filter design (Tustin / bilinear)
 *
 * Discrete form of H(s) = s / (s + ωc), ωc = 2π·fc. Pair with @ref HighPass
 * or a first-order difference equation. Invalid @p fc / @p Ts → zero coeffs.
 *
 * @param fc Cutoff frequency [Hz] (0 < fc < Nyquist)
 * @param Ts Sample time [s] (> 0)
 * @return FirstOrderCoeffs\<T\> DSP coefficients
 * @see lowpass_1st()
 */
template<typename T = float>
[[nodiscard]] constexpr FirstOrderCoeffs<T> highpass_1st(T fc, T Ts) {
    if (!detail::valid_cutoff(fc, Ts)) {
        return FirstOrderCoeffs<T>{};
    }
    // Bilinear of H(s) = s/(s+ωc). Pin a1 from the pole; force exact DC null
    // (b1 = −b0) and unit Nyquist gain H(−1) = 1 by construction:
    //   (b0 − b1)/(1 − a1) = 1 with b1 = −b0  ⇒  b0 = (1 − a1)/2.
    const T omega_c = T{2} * damp::numbers::pi_v<T> * fc;
    const T k = T{2} / Ts;
    const T denom = omega_c + k;
    const T a1 = (omega_c - k) / denom;
    const T b0 = (T{1} - a1) / T{2};
    return FirstOrderCoeffs<T>{
        .b0 = b0,
        .b1 = -b0,
        .a1 = a1,
    };
}

/**
 * @brief First-order low-pass filter design (continuous-time)
 *
 * H(s) = 1 / (τs + 1) with τ = 1/(2π fc) (unit DC gain, pole at −1/τ).
 * Coefficients are ascending powers of s (den[0] constant).
 * Invalid @p fc (≤ 0) returns the zero TF (no NaN pole).
 *
 * @param fc Cutoff frequency [Hz] (must be > 0)
 * @return TransferFunction<2, 2, T> continuous-time transfer function
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunction<2, 2, T> lowpass_1st(T fc) {
    if (!(fc > T{0})) {
        // Zero gain, monic pole at −1 (realizable companion; no NaN τ)
        return TransferFunction<2, 2, T>{{T{0}, T{0}}, {T{1}, T{1}}};
    }
    const T tau = T{1} / (T{2} * damp::numbers::pi_v<T> * fc);
    // num = 1, den = 1 + τs  →  H = 1/(τs+1)
    return TransferFunction<2, 2, T>{{T{1}, T{0}}, {T{1}, tau}};
}

/**
 * @brief Second-order low-pass filter design
 *
 * Designs a discrete-time second-order low-pass filter with given cutoff frequency and damping.
 * Returns DSP coefficients ready for runtime use.
 *
 * @param fc Cutoff frequency [Hz] (0 < fc < Nyquist)
 * @param Ts Sample time [s] (> 0)
 * @param zeta Damping ratio (0.707 for Butterworth); must be > 0
 *  * @return SecondOrderCoeffs\<T\> DSP coefficients; all zero if arguments invalid
 */
template<typename T = float>
[[nodiscard]] constexpr SecondOrderCoeffs<T> lowpass_2nd(T fc, T Ts, T zeta = static_cast<T>(0.707)) {
    if (!detail::valid_cutoff(fc, Ts) || !(zeta > T{0})) {
        return SecondOrderCoeffs<T>{};
    }
    // Discretize H(s) = ω₀² / (s² + 2ζω₀s + ω₀²) with the bilinear transform
    // s ← k·(1 − z⁻¹)/(1 + z⁻¹), k = 2/Ts. Clearing (1 + z⁻¹)² gives a numerator
    // proportional to (1 + z⁻¹)² (tap shape 1 : 2 : 1) over a denominator with
    //   a0 = k² + 2ζω₀k + ω₀².
    const T omega_0 = T{2} * damp::numbers::pi_v<T> * fc;
    const T k = T{2} / Ts; // bilinear gain (no extra pre-warp; fc already gated vs Nyquist)
    const T k_sq = k * k;
    const T omega_0_sq = omega_0 * omega_0;
    const T two_zeta_omega = T{2} * zeta * omega_0;

    const T a0 = k_sq + two_zeta_omega * k + omega_0_sq;

    const T a1 = (T{2} * omega_0_sq - T{2} * k_sq) / a0;
    const T a2 = (k_sq - two_zeta_omega * k + omega_0_sq) / a0;

    // The numerator is ω₀²·(1 + z⁻¹)², so the taps are in fixed ratio 1 : 2 : 1 and
    // their only free parameter is the overall scale. Pin that scale from the
    // unit-DC-gain identity H(1) = 1 ⇔ (b0 + b1 + b2) = (1 + a1 + a2): deriving the
    // taps from that sum makes unity DC gain hold *by construction*, so it survives
    // -ffast-math reassociation in downstream builds instead of depending on the raw
    // bilinear numerator/denominator terms cancelling exactly.
    const T dc_sum = T{1} + a1 + a2;
    return SecondOrderCoeffs<T>{
        .b0 = dc_sum / T{4},
        .b1 = dc_sum / T{2},
        .b2 = dc_sum / T{4},
        .a1 = a1,
        .a2 = a2,
    };
}

/**
 * @brief Second-order low-pass filter design (continuous-time)
 *
 * H(s) = ω₀² / (s² + 2ζω₀ s + ω₀²), ω₀ = 2π fc.
 * Coefficients are ascending powers of s (den[0] constant, den[2] monic s²).
 *
 * Named `*_continuous` so it does not overload-clash with the discrete
 * `lowpass_2nd(fc, Ts, ζ)` (both would otherwise match a two-argument call).
 * One-argument `lowpass_2nd(fc)` is the default-ζ continuous convenience.
 *
 * @param fc Cutoff frequency [Hz]
 * @param zeta Damping ratio (≈1/√2 for 2nd-order Butterworth)
 * @return TransferFunction<3, 3, T> continuous-time transfer function
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunction<3, 3, T>
lowpass_2nd_continuous(T fc, T zeta = static_cast<T>(0.7071067811865476)) {
    if (!(fc > T{0}) || !(zeta > T{0})) {
        // Degenerate: return monic s² den with zero num (no NaN poles).
        return TransferFunction<3, 3, T>{{T{0}, T{0}, T{0}}, {T{0}, T{0}, T{1}}};
    }
    const T omega_0 = T{2} * damp::numbers::pi_v<T> * fc;
    const T omega_0_sq = omega_0 * omega_0;
    const T two_zeta_omega = T{2} * zeta * omega_0;

    return TransferFunction<3, 3, T>{
        {omega_0_sq, T{0}, T{0}},          // num: ω₀²
        {omega_0_sq, two_zeta_omega, T{1}} // den: ω₀² + 2ζω₀ s + s²
    };
}

/**
 * @brief Second-order continuous low-pass at default Butterworth damping (ζ = 1/√2)
 *
 * @see lowpass_2nd_continuous for a custom damping ratio
 * @see lowpass_2nd(fc, Ts, zeta) for the discrete coefficient form
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunction<3, 3, T> lowpass_2nd(T fc) {
    return lowpass_2nd_continuous(fc);
}

/**
 * @brief First-order high-pass filter design (continuous-time)
 *
 * H(s) = s / (s + ωc) with ωc = 2π fc (DC gain 0, high-frequency gain 1).
 * Coefficients are ascending powers of s (den[0] constant).
 * Invalid @p fc (≤ 0) returns the zero TF (no NaN pole).
 *
 * @param fc Cutoff frequency [Hz] (must be > 0)
 * @return TransferFunction<2, 2, T> continuous-time transfer function
 * @note Compare with MATLAB®'s butter(1, 2*pi*fc, 'high', 's').
 * @see highpass_1st(fc, Ts) for the discrete Tustin coefficient form
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunction<2, 2, T> highpass_1st(T fc) {
    if (!(fc > T{0})) {
        // Zero gain, monic pole at −1 (realizable companion; no NaN ωc)
        return TransferFunction<2, 2, T>{{T{0}, T{0}}, {T{1}, T{1}}};
    }
    const T omega_c = T{2} * damp::numbers::pi_v<T> * fc;
    // num = s, den = ωc + s  →  H = s/(s+ωc)
    return TransferFunction<2, 2, T>{{T{0}, T{1}}, {omega_c, T{1}}};
}

/**
 * @brief Second-order high-pass filter design (continuous-time)
 *
 * H(s) = s² / (s² + (ω₀/Q) s + ω₀²), ω₀ = 2π fc.
 * Coefficients are ascending powers of s (den[0] constant, den[2] monic s²).
 *
 * Named `*_continuous` so it does not overload-clash with the discrete
 * `highpass_2nd(fc, Ts, Q)` (both would otherwise match a two-argument call).
 * One-argument `highpass_2nd(fc)` is the default-Q Butterworth convenience.
 *
 * @param fc Cutoff frequency [Hz]
 * @param Q  Quality factor (1/√2 for 2nd-order Butterworth)
 * @return TransferFunction<3, 3, T> continuous-time transfer function
 * @note Compare with MATLAB®'s butter(2, 2*pi*fc, 'high', 's') (Q = 1/√2).
 * @see highpass_2nd(fc, Ts, Q) for the discrete RBJ coefficient form
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunction<3, 3, T> highpass_2nd_continuous(T fc, T Q) {
    if (!(fc > T{0}) || !detail::valid_Q(Q)) {
        // Degenerate: return monic s² den with zero num (no NaN poles).
        return TransferFunction<3, 3, T>{{T{0}, T{0}, T{0}}, {T{0}, T{0}, T{1}}};
    }
    const T omega_0 = T{2} * damp::numbers::pi_v<T> * fc;
    const T omega_0_sq = omega_0 * omega_0;
    const T omega_over_Q = omega_0 / Q;

    return TransferFunction<3, 3, T>{
        {T{0}, T{0}, T{1}},              // num: s²
        {omega_0_sq, omega_over_Q, T{1}} // den: ω₀² + (ω₀/Q) s + s²
    };
}

/**
 * @brief Second-order continuous high-pass at default Butterworth Q (1/√2)
 *
 * @see highpass_2nd_continuous for a custom quality factor
 * @see highpass_2nd(fc, Ts, Q) for the discrete coefficient form
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunction<3, 3, T> highpass_2nd(T fc) {
    return highpass_2nd_continuous(fc, damp::numbers::inv_sqrt2_v<T>);
}

/**
 * @brief Second-order band-pass filter design (continuous-time)
 *
 * H(s) = (ω₀/Q) s / (s² + (ω₀/Q) s + ω₀²), ω₀ = 2π f0 (unit gain at ω₀).
 * Coefficients are ascending powers of s (den[0] constant, den[2] monic s²).
 * Named `*_continuous` so the discrete `bandpass(f0, Q, Ts)` stays three-arg.
 *
 * @param f0 Center frequency [Hz]
 * @param Q  Quality factor (higher = narrower band)
 * @return TransferFunction<3, 3, T> continuous-time transfer function
 * @note Compare with MATLAB®'s iirpeak analog prototype (constant 0 dB peak).
 * @see bandpass(f0, Q, Ts) for the discrete RBJ coefficient form
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunction<3, 3, T> bandpass_continuous(T f0, T Q) {
    if (!(f0 > T{0}) || !detail::valid_Q(Q)) {
        return TransferFunction<3, 3, T>{{T{0}, T{0}, T{0}}, {T{0}, T{0}, T{1}}};
    }
    const T omega_0 = T{2} * damp::numbers::pi_v<T> * f0;
    const T omega_0_sq = omega_0 * omega_0;
    const T omega_over_Q = omega_0 / Q;

    return TransferFunction<3, 3, T>{
        {T{0}, omega_over_Q, T{0}},      // num: (ω₀/Q) s
        {omega_0_sq, omega_over_Q, T{1}} // den: ω₀² + (ω₀/Q) s + s²
    };
}

/**
 * @brief Second-order band-reject (notch) filter design (continuous-time)
 *
 * H(s) = (s² + ω₀²) / (s² + (ω₀/Q) s + ω₀²), ω₀ = 2π f0 (null at ω₀, unit DC/HF).
 * Coefficients are ascending powers of s (den[0] constant, den[2] monic s²).
 * Named `*_continuous` so the discrete `notch(f0, Q, Ts)` stays three-arg.
 *
 * @param f0 Notch (center) frequency [Hz]
 * @param Q  Quality factor (higher = narrower notch)
 * @return TransferFunction<3, 3, T> continuous-time transfer function
 * @note Compare with MATLAB®'s iirnotch analog prototype.
 * @see notch(f0, Q, Ts) for the discrete RBJ coefficient form
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunction<3, 3, T> notch_continuous(T f0, T Q) {
    if (!(f0 > T{0}) || !detail::valid_Q(Q)) {
        return TransferFunction<3, 3, T>{{T{0}, T{0}, T{0}}, {T{0}, T{0}, T{1}}};
    }
    const T omega_0 = T{2} * damp::numbers::pi_v<T> * f0;
    const T omega_0_sq = omega_0 * omega_0;
    const T omega_over_Q = omega_0 / Q;

    return TransferFunction<3, 3, T>{
        {omega_0_sq, T{0}, T{1}},        // num: ω₀² + s²
        {omega_0_sq, omega_over_Q, T{1}} // den: ω₀² + (ω₀/Q) s + s²
    };
}

/**
 * @brief Second-order all-pass filter design (continuous-time)
 *
 * H(s) = (s² − (ω₀/Q) s + ω₀²) / (s² + (ω₀/Q) s + ω₀²), ω₀ = 2π f0.
 * |H(jω)| = 1 for all ω; phase lags through 180° around ω₀ (Q sets the slope).
 * Coefficients are ascending powers of s (den[0] constant, den[2] monic s²).
 *
 * @param f0 Center frequency [Hz]
 * @param Q  Quality factor (higher = steeper phase around ω₀)
 * @return TransferFunction<3, 3, T> continuous-time transfer function
 * @note Compare with MATLAB®'s tf([1, −ω₀/Q, ω₀²], [1, ω₀/Q, ω₀²]) (descending powers).
 * @see allpass(f0, Q, Ts) for the discrete Tustin coefficient form
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunction<3, 3, T> allpass_2nd_continuous(T f0, T Q) {
    if (!(f0 > T{0}) || !detail::valid_Q(Q)) {
        return TransferFunction<3, 3, T>{{T{0}, T{0}, T{0}}, {T{0}, T{0}, T{1}}};
    }
    const T omega_0 = T{2} * damp::numbers::pi_v<T> * f0;
    const T omega_0_sq = omega_0 * omega_0;
    const T omega_over_Q = omega_0 / Q;

    return TransferFunction<3, 3, T>{
        {omega_0_sq, -omega_over_Q, T{1}}, // num: ω₀² − (ω₀/Q) s + s²
        {omega_0_sq, omega_over_Q, T{1}}   // den: ω₀² + (ω₀/Q) s + s²
    };
}

/**
 * @brief Butterworth low-pass filter design
 *
 * Continuous Butterworth of order 1–4 at cutoff fc [Hz], built as cascade of
 * first- and second-order factors of H(s) = 1 / B_n(s/ωc) with the standard
 * quadratic damping ratios:
 * - n=1: single real pole
 * - n=2: ζ = sin(π/4) = 1/√2
 * - n=3: real pole × quadratic with ζ = sin(π/6) = 1/2
 * - n=4: quadratics with ζ = sin(π/8), sin(3π/8)
 *
 * @param fc Cutoff frequency [Hz]
 * @return StateSpace system representing the filter
 */
template<size_t Order, typename T = double>
    requires(Order >= 1 && Order <= 4)
[[nodiscard]] constexpr auto butterworth_lowpass(T fc) {
    // Companion realization succeeds for monic/proper low-pass TFs when fc > 0.
    // Invalid fc falls back through the continuous designers' zero-gain monic dens.
    if constexpr (Order == 1) {
        auto opt = lowpass_1st<T>(fc).to_state_space();
        return opt ? *opt : StateSpace<1, 1, 1, T>{};
    } else if constexpr (Order == 2) {
        // Butterworth ζ = 1/√2 = sin(π/4)
        auto opt = lowpass_2nd_continuous(fc, damp::numbers::inv_sqrt2_v<T>).to_state_space();
        return opt ? *opt : StateSpace<2, 1, 1, T>{};
    } else if constexpr (Order == 3) {
        // B3(s) ∝ (s+1)(s² + s + 1) → ζ = 1/2 on the quadratic at the same ωc
        auto tf1 = lowpass_1st<T>(fc);
        auto tf2 = lowpass_2nd_continuous(fc, static_cast<T>(0.5));
        auto opt = (tf1 * tf2).to_state_space();
        return opt ? *opt : StateSpace<3, 1, 1, T>{};
    } else { // Order == 4
        // B4 factors: ζ_k = sin((2k−1)π / 8), k = 1, 2
        const T pi = damp::numbers::pi_v<T>;
        const T z1 = damp::sin(pi / T{8});
        const T z2 = damp::sin(T{3} * pi / T{8});
        auto    a = lowpass_2nd_continuous(fc, z1);
        auto    b = lowpass_2nd_continuous(fc, z2);
        auto    opt = (a * b).to_state_space();
        return opt ? *opt : StateSpace<4, 1, 1, T>{};
    }
}

/**
 * @brief Butterworth low-pass filter design (discrete-time)
 *
 * Creates a discrete-time Butterworth low-pass filter of specified order.
 * Butterworth filters have maximally flat magnitude response in passband.
 *
 * @param fc Cutoff frequency [Hz]
 * @param Ts Sample time [s]
 *  *  * @return StateSpace system representing the discrete-time filter
 */
template<size_t Order, typename T = float>
    requires(Order >= 1 && Order <= 4)
[[nodiscard]] constexpr auto butterworth_lowpass(T fc, T Ts) {
    // Fail closed on invalid rates: return continuous companion with Ts = 0
    // rather than UB from *nullopt on a failed Tustin step.
    if (!(Ts > T{0}) || !(fc > T{0}) || !((T{2} * fc * Ts) < T{1})) {
        return butterworth_lowpass<Order, T>(fc);
    }
    auto sys_c = butterworth_lowpass<Order, T>(fc);
    auto sys_d = discretize(sys_c, Ts, DiscretizationMethod::Tustin);
    if (!sys_d) {
        return sys_c;
    }
    return *sys_d;
}

/**
 * @brief First-order Pade approximation of time delay
 *
 * Approximates e^(-sT) with a first-order rational transfer function:
 * H(s) = (1 - sT/2) / (1 + sT/2)
 *
 * @param T_delay Time delay [s]
 *  * @return TransferFunction<2, 2, T> representing the delay approximation
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunction<2, 2, T> pade_delay_1st(T T_delay) {
    const T half_T = T_delay / T{2};
    return TransferFunction<2, 2, T>{
        {T{1}, -half_T}, // num: 1 - sT/2
        {T{1}, half_T}   // den: 1 + sT/2
    };
}

/**
 * @brief First-order Pade approximation of time delay (discrete-time)
 *
 * Discrete-time version of the first-order Pade delay approximation.
 *
 * @param T_delay Time delay [s]
 * @param Ts Sample time [s]
 *  * @return TransferFunction<2, 2, T> discrete-time delay approximation
 */
template<typename T = float>
[[nodiscard]] constexpr TransferFunction<2, 2, T> pade_delay_1st(T T_delay, T Ts) {
    if (!(Ts > T{0}) || !(T_delay >= T{0})) {
        return TransferFunction<2, 2, T>{{T{0}, T{0}}, {T{1}, T{0}}};
    }
    auto tf_c = pade_delay_1st<T>(T_delay);
    auto sys_c_opt = tf_c.to_state_space();
    if (!sys_c_opt) {
        return TransferFunction<2, 2, T>{{T{0}, T{0}}, {T{1}, T{0}}};
    }
    auto sys_d_opt = discretize(*sys_c_opt, Ts, DiscretizationMethod::Tustin);
    if (!sys_d_opt) {
        return TransferFunction<2, 2, T>{{T{0}, T{0}}, {T{1}, T{0}}};
    }
    return TransferFunction<2, 2, T>::from_state_space(*sys_d_opt);
}

/**
 * @brief Second-order Pade approximation of time delay
 *
 * Approximates e^(-sT) with a second-order rational transfer function:
 * H(s) = (1 - sT/2 + (sT)²/12) / (1 + sT/2 + (sT)²/12)
 *
 * @param T_delay Time delay [s]
 *  * @return TransferFunction<3, 3, T> representing the delay approximation
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunction<3, 3, T> pade_delay_2nd(T T_delay) {
    const T half_T = T_delay / T{2};
    const T T_sq_12 = T_delay * T_delay / T{12};
    return TransferFunction<3, 3, T>{
        {T{1}, -half_T, T_sq_12}, // num: 1 - sT/2 + (sT)²/12
        {T{1}, half_T, T_sq_12}   // den: 1 + sT/2 + (sT)²/12
    };
}

/**
 * @brief Second-order Pade approximation of time delay (discrete-time)
 *
 * Discrete-time version of the second-order Pade delay approximation.
 *
 * @param T_delay Time delay [s]
 * @param Ts Sample time [s]
 *  * @return TransferFunction<3, 3, T> discrete-time delay approximation
 */
template<typename T = float>
[[nodiscard]] constexpr TransferFunction<3, 3, T> pade_delay_2nd(T T_delay, T Ts) {
    if (!(Ts > T{0}) || !(T_delay >= T{0})) {
        return TransferFunction<3, 3, T>{{T{0}, T{0}, T{0}}, {T{0}, T{0}, T{1}}};
    }
    auto tf_c = pade_delay_2nd<T>(T_delay);
    auto sys_c_opt = tf_c.to_state_space();
    if (!sys_c_opt) {
        return TransferFunction<3, 3, T>{{T{0}, T{0}, T{0}}, {T{0}, T{0}, T{1}}};
    }
    auto sys_d_opt = discretize(*sys_c_opt, Ts, DiscretizationMethod::Tustin);
    if (!sys_d_opt) {
        return TransferFunction<3, 3, T>{{T{0}, T{0}, T{0}}, {T{0}, T{0}, T{1}}};
    }
    return TransferFunction<3, 3, T>::from_state_space(*sys_d_opt);
}

/**
 * @brief Convert StateSpace system to first-order DSP coefficients
 *
 * Assumes the system is first-order (1 state). If continuous-time, discretizes first.
 *
 * @param sys State-space system
 * @param Ts Sample time (only used if sys is continuous)
 *  * @return FirstOrderCoeffs\<T\> DSP coefficients
 */
template<typename T = float>
[[nodiscard]] constexpr FirstOrderCoeffs<T> to_coeffs(const StateSpace<1, 1, 1, T>& sys, T Ts = T{0}) {
    StateSpace<1, 1, 1, T> sys_d = sys;
    if (sys.Ts == T{0}) {
        if (!(Ts > T{0})) {
            return FirstOrderCoeffs<T>{};
        }
        auto opt = discretize(sys, Ts, DiscretizationMethod::Tustin);
        if (!opt) {
            return FirstOrderCoeffs<T>{};
        }
        sys_d = *opt;
    }

    return FirstOrderCoeffs<T>{
        .b0 = sys_d.D(0, 0),
        .b1 = (sys_d.C(0, 0) * sys_d.B(0, 0)) - (sys_d.D(0, 0) * sys_d.A(0, 0)),
        .a1 = -sys_d.A(0, 0),
    };
}

/**
 * @brief Convert TransferFunction to first-order DSP coefficients
 *
 * Assumes the transfer function is first-order. Discretizes using bilinear transform.
 *
 * @param tf Transfer function
 * @param Ts Sample time
 *  * @return FirstOrderCoeffs\<T\> DSP coefficients
 */
template<typename T = float>
[[nodiscard]] constexpr FirstOrderCoeffs<T> to_coeffs(const TransferFunction<1, 2, T>& tf, T Ts) {
    if (!(Ts > T{0})) {
        return FirstOrderCoeffs<T>{};
    }
    const auto sys_c = tf.to_state_space();
    if (!sys_c) {
        return FirstOrderCoeffs<T>{};
    }
    const auto sys_d = discretize(*sys_c, Ts, DiscretizationMethod::Tustin);
    if (!sys_d) {
        return FirstOrderCoeffs<T>{};
    }
    return to_coeffs(*sys_d);
}

template<typename T = float>
[[nodiscard]] constexpr FirstOrderCoeffs<T> to_coeffs(const TransferFunction<2, 2, T>& tf, T Ts) {
    if (!(Ts > T{0})) {
        return FirstOrderCoeffs<T>{};
    }
    const auto sys_c = tf.to_state_space();
    if (!sys_c) {
        return FirstOrderCoeffs<T>{};
    }
    const auto sys_d = discretize(*sys_c, Ts, DiscretizationMethod::Tustin);
    if (!sys_d) {
        return FirstOrderCoeffs<T>{};
    }
    return to_coeffs(*sys_d);
}

/**
 * @brief Convert StateSpace system to second-order DSP coefficients
 *
 * Assumes the system is second-order (2 states). If continuous-time, discretizes first.
 *
 * @param sys State-space system
 * @param Ts Sample time (only used if sys is continuous)
 *  * @return SecondOrderCoeffs\<T\> DSP coefficients
 */
template<typename T = float>
[[nodiscard]] constexpr SecondOrderCoeffs<T> to_coeffs(const StateSpace<2, 1, 1, T>& sys, T Ts = T{0}) {
    StateSpace<2, 1, 1, T> sys_d = sys;
    if (sys.Ts == T{0}) {
        if (!(Ts > T{0})) {
            return SecondOrderCoeffs<T>{};
        }
        auto opt = discretize(sys, Ts, DiscretizationMethod::Tustin);
        if (!opt) {
            return SecondOrderCoeffs<T>{};
        }
        sys_d = *opt;
    }

    const T A00 = sys_d.A(0, 0);
    const T A01 = sys_d.A(0, 1);
    const T A10 = sys_d.A(1, 0);
    const T A11 = sys_d.A(1, 1);
    const T B0 = sys_d.B(0, 0);
    const T B1 = sys_d.B(1, 0);
    const T C0 = sys_d.C(0, 0);
    const T C1 = sys_d.C(0, 1);
    const T D0 = sys_d.D(0, 0);

    // H(z) = D + C·adj(zI−A)·B / det(zI−A). Denominator is the characteristic
    // polynomial (a1 = −trace, a2 = det); the numerator is the full expansion,
    // not just C·B — the D·trace / D·det and cross terms matter whenever D ≠ 0
    // (e.g. any Tustin-discretized system).
    const T a1 = -(A00 + A11);
    const T a2 = (A00 * A11) - (A01 * A10);
    return SecondOrderCoeffs<T>{
        .b0 = D0,
        .b1 = (C0 * B0) + (C1 * B1) + (D0 * a1),
        .b2 = (D0 * a2) + (C0 * ((A01 * B1) - (A11 * B0))) + (C1 * ((A10 * B0) - (A00 * B1))),
        .a1 = a1,
        .a2 = a2,
    };
}

/**
 * @brief Convert TransferFunction to second-order DSP coefficients
 *
 * Assumes the transfer function is second-order. Discretizes using bilinear transform.
 *
 * @param tf Transfer function
 * @param Ts Sample time
 *  * @return SecondOrderCoeffs\<T\> DSP coefficients
 */
template<typename T = float>
[[nodiscard]] constexpr SecondOrderCoeffs<T> to_coeffs(const TransferFunction<3, 3, T>& tf, T Ts) {
    if (!(Ts > T{0})) {
        return SecondOrderCoeffs<T>{};
    }
    const auto sys_c = tf.to_state_space();
    if (!sys_c) {
        return SecondOrderCoeffs<T>{};
    }
    const auto sys_d = discretize(*sys_c, Ts, DiscretizationMethod::Tustin);
    if (!sys_d) {
        return SecondOrderCoeffs<T>{};
    }
    return to_coeffs(*sys_d);
}

// ============================================================================
// Biquad (second-order IIR) designs — RBJ audio-EQ cookbook formulas
// ============================================================================
//
// Coefficients use the normalized difference equation
//   y[n] = b0·x[n] + b1·x[n-1] + b2·x[n-2] − a1·y[n-1] − a2·y[n-2]
// matching SecondOrderCoeffs and the Biquad runtime. Frequencies are designed
// directly in the digital domain: ω₀ = 2π·f·Ts. Quality factor Q controls
// bandwidth (BW in octaves ≈ asinh(1/2Q)·2/ln2 near ω₀); larger Q is narrower.

namespace detail {
template<typename T>
[[nodiscard]] constexpr SecondOrderCoeffs<T> normalize_biquad(T b0, T b1, T b2, T a0, T a1, T a2) {
    // Guard a0 ≈ 0 (degenerate Q/freq) — return zero section, not ±∞ coeffs.
    if (!(damp::abs(a0) > default_tol<T>())) {
        return SecondOrderCoeffs<T>{};
    }
    const T inv = T{1} / a0;
    return SecondOrderCoeffs<T>{b0 * inv, b1 * inv, b2 * inv, a1 * inv, a2 * inv};
}
} // namespace detail

/**
 * @brief Second-order band-reject (notch) filter.
 *
 * Rejects a narrow band around f0 (gain → 0 at f0) while passing the rest
 * (gain → 1). Transfer function:
 *
 *     H(z) = (1 − 2cosω₀ z⁻¹ + z⁻²) / (1 + α) / (… z⁻¹ …),  α = sinω₀ / (2Q)
 *
 * @note Compare with MATLAB®'s iirnotch(w0, bw).
 * @param f0 Notch (center) frequency [Hz]
 * @param Q  Quality factor (higher = narrower notch)
 * @param Ts Sample time [s]
 * @return SecondOrderCoeffs\<T\> normalized biquad coefficients
 * @see "Cookbook formulae for audio EQ biquad filter coefficients" (Bristow-Johnson)
 */
template<typename T = float>
[[nodiscard]] constexpr SecondOrderCoeffs<T> notch(T f0, T Q, T Ts) {
    if (!detail::valid_cutoff(f0, Ts) || !detail::valid_Q(Q)) {
        return SecondOrderCoeffs<T>{};
    }
    const T w0 = T{2} * damp::numbers::pi_v<T> * f0 * Ts;
    const auto [s, c] = damp::sincos(w0);
    const T alpha = s / (T{2} * Q);
    return detail::normalize_biquad<T>(T{1}, T{-2} * c, T{1}, T{1} + alpha, T{-2} * c, T{1} - alpha);
}

/**
 * @brief Second-order band-pass filter (constant 0 dB peak gain).
 *
 * Passes a band around f0 (gain → 1 at f0), rejecting DC and high frequencies.
 *
 * @param f0 Center frequency [Hz]
 * @param Q  Quality factor (higher = narrower band)
 * @param Ts Sample time [s]
 * @return SecondOrderCoeffs\<T\> normalized biquad coefficients
 * @see "Cookbook formulae for audio EQ biquad filter coefficients" (Bristow-Johnson)
 */
template<typename T = float>
[[nodiscard]] constexpr SecondOrderCoeffs<T> bandpass(T f0, T Q, T Ts) {
    if (!detail::valid_cutoff(f0, Ts) || !detail::valid_Q(Q)) {
        return SecondOrderCoeffs<T>{};
    }
    const T w0 = T{2} * damp::numbers::pi_v<T> * f0 * Ts;
    const auto [s, c] = damp::sincos(w0);
    const T alpha = s / (T{2} * Q);
    return detail::normalize_biquad<T>(alpha, T{0}, -alpha, T{1} + alpha, T{-2} * c, T{1} - alpha);
}

/**
 * @brief Second-order high-pass filter (RBJ).
 *
 * Counterpart to lowpass_2nd. Note this family is parameterized by Q
 * (= 1 / 2ζ); the default Q = 1/√2 is the maximally-flat (Butterworth) response.
 *
 * @param fc Cutoff frequency [Hz]
 * @param Ts Sample time [s]
 * @param Q  Quality factor (default 1/√2 = Butterworth)
 * @return SecondOrderCoeffs\<T\> normalized biquad coefficients
 * @see lowpass_2nd()
 */
template<typename T = float>
[[nodiscard]] constexpr SecondOrderCoeffs<T> highpass_2nd(T fc, T Ts, T Q = damp::numbers::inv_sqrt2_v<T>) {
    if (!detail::valid_cutoff(fc, Ts) || !detail::valid_Q(Q)) {
        return SecondOrderCoeffs<T>{};
    }
    const T w0 = T{2} * damp::numbers::pi_v<T> * fc * Ts;
    const auto [s, c] = damp::sincos(w0);
    const T alpha = s / (T{2} * Q);
    // RBJ highpass: num taps ∝ (1+cos) · {1, −2, 1} → exact DC null by shape.
    // Normalize a0; DC identity H(1) = 0 holds independent of -ffast-math.
    const T b0 = (T{1} + c) / T{2};
    return detail::normalize_biquad<T>(b0, -(T{1} + c), b0, T{1} + alpha, T{-2} * c, T{1} - alpha);
}

/**
 * @brief Second-order all-pass filter (Tustin of the analog prototype).
 *
 * Discrete bilinear of @ref allpass_2nd_continuous. |H| ≈ 1 across the band;
 * phase wraps around f0. Invalid @p f0 / @p Q / @p Ts → zero coeffs.
 *
 * @param f0 Center frequency [Hz] (0 < f0 < Nyquist)
 * @param Q  Quality factor (higher = steeper phase around f0)
 * @param Ts Sample time [s] (> 0)
 * @return SecondOrderCoeffs\<T\> DSP coefficients
 * @see allpass_2nd_continuous
 */
template<typename T = float>
[[nodiscard]] constexpr SecondOrderCoeffs<T> allpass(T f0, T Q, T Ts) {
    if (!detail::valid_cutoff(f0, Ts) || !detail::valid_Q(Q)) {
        return SecondOrderCoeffs<T>{};
    }
    return to_coeffs(allpass_2nd_continuous<T>(f0, Q), Ts);
}

/**
 * @brief Peaking (bell) EQ filter: boost or cut a band around f0.
 *
 * @param f0      Center frequency [Hz]
 * @param Q       Quality factor (higher = narrower bell)
 * @param gain_db Peak gain at f0 [dB] (positive = boost, negative = cut)
 * @param Ts      Sample time [s]
 * @return SecondOrderCoeffs\<T\> normalized biquad coefficients
 * @see "Cookbook formulae for audio EQ biquad filter coefficients" (Bristow-Johnson)
 */
template<typename T = float>
[[nodiscard]] constexpr SecondOrderCoeffs<T> peaking(T f0, T Q, T gain_db, T Ts) {
    if (!detail::valid_cutoff(f0, Ts) || !detail::valid_Q(Q)) {
        return SecondOrderCoeffs<T>{};
    }
    const T A = damp::pow(T{10}, gain_db / T{40});
    const T w0 = T{2} * damp::numbers::pi_v<T> * f0 * Ts;
    const auto [s, c] = damp::sincos(w0);
    const T alpha = s / (T{2} * Q);
    return detail::normalize_biquad<T>(
        T{1} + (alpha * A), T{-2} * c, T{1} - (alpha * A),
        T{1} + (alpha / A), T{-2} * c, T{1} - (alpha / A)
    );
}

/**
 * @brief Low-shelf EQ filter: boost or cut everything below fc.
 *
 * @param fc      Shelf corner frequency [Hz]
 * @param gain_db Shelf gain [dB]
 * @param Ts      Sample time [s]
 * @param Q       Shelf shape (default 1/√2)
 * @return SecondOrderCoeffs\<T\> normalized biquad coefficients
 * @see "Cookbook formulae for audio EQ biquad filter coefficients" (Bristow-Johnson)
 */
template<typename T = float>
[[nodiscard]] constexpr SecondOrderCoeffs<T> lowshelf(T fc, T gain_db, T Ts, T Q = damp::numbers::inv_sqrt2_v<T>) {
    if (!detail::valid_cutoff(fc, Ts) || !detail::valid_Q(Q)) {
        return SecondOrderCoeffs<T>{};
    }
    const T A = damp::pow(T{10}, gain_db / T{40});
    const T w0 = T{2} * damp::numbers::pi_v<T> * fc * Ts;
    const auto [s, c] = damp::sincos(w0);
    const T alpha = s / (T{2} * Q);
    const T tsa = T{2} * damp::sqrt(A) * alpha;
    const T Ap = A + T{1};
    const T Am = A - T{1};
    return detail::normalize_biquad<T>(
        A * (Ap - (Am * c) + tsa),
        T{2} * A * (Am - (Ap * c)),
        A * (Ap - (Am * c) - tsa),
        Ap + (Am * c) + tsa,
        T{-2} * (Am + (Ap * c)),
        Ap + (Am * c) - tsa
    );
}

/**
 * @brief High-shelf EQ filter: boost or cut everything above fc.
 *
 * @param fc      Shelf corner frequency [Hz]
 * @param gain_db Shelf gain [dB]
 * @param Ts      Sample time [s]
 * @param Q       Shelf shape (default 1/√2)
 * @return SecondOrderCoeffs\<T\> normalized biquad coefficients
 * @see "Cookbook formulae for audio EQ biquad filter coefficients" (Bristow-Johnson)
 */
template<typename T = float>
[[nodiscard]] constexpr SecondOrderCoeffs<T> highshelf(T fc, T gain_db, T Ts, T Q = damp::numbers::inv_sqrt2_v<T>) {
    if (!detail::valid_cutoff(fc, Ts) || !detail::valid_Q(Q)) {
        return SecondOrderCoeffs<T>{};
    }
    const T A = damp::pow(T{10}, gain_db / T{40});
    const T w0 = T{2} * damp::numbers::pi_v<T> * fc * Ts;
    const auto [s, c] = damp::sincos(w0);
    const T alpha = s / (T{2} * Q);
    const T tsa = T{2} * damp::sqrt(A) * alpha;
    const T Ap = A + T{1};
    const T Am = A - T{1};
    return detail::normalize_biquad<T>(
        A * (Ap + (Am * c) + tsa),
        T{-2} * A * (Am + (Ap * c)),
        A * (Ap + (Am * c) - tsa),
        Ap - (Am * c) + tsa,
        T{2} * (Am - (Ap * c)),
        Ap - (Am * c) - tsa
    );
}
} // namespace design

} // namespace damp
