// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file iir_design.hpp
 * @brief IIR coefficient design (analog TFs, Tustin/RBJ/Butterworth) and to_coeffs
 */


#include "damp/backend.hpp"
#include "damp/design/model_result.hpp"
#include "damp/design/status.hpp"
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
struct [[nodiscard]] FirstOrderCoeffs {
    T b0{0}, b1{0}, a1{0};

    bool         success{true};            ///< status == DesignStatus::ok (a hand-built set is a valid design)
    DesignStatus status{DesignStatus::ok}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /**
     * @brief Convert coefficients to different scalar type
     * @tparam U Target scalar type
     * @return FirstOrderCoeffs\<U\> with converted coefficients
     */
    template<typename U>
    [[nodiscard]] constexpr FirstOrderCoeffs<U> as() const {
        return {static_cast<U>(b0), static_cast<U>(b1), static_cast<U>(a1), success, status};
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
struct [[nodiscard]] SecondOrderCoeffs {
    T b0{0}, b1{0}, b2{0}, a1{0}, a2{0};

    bool         success{true};            ///< status == DesignStatus::ok (a hand-built set is a valid design)
    DesignStatus status{DesignStatus::ok}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /**
     * @brief Convert coefficients to different scalar type
     * @tparam U Target scalar type
     * @return SecondOrderCoeffs\<U\> with converted coefficients
     */
    template<typename U>
    [[nodiscard]] constexpr SecondOrderCoeffs<U> as() const {
        return {static_cast<U>(b0), static_cast<U>(b1), static_cast<U>(b2), static_cast<U>(a1), static_cast<U>(a2), success, status};
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

// ---- Coefficients from a model ---------------------------------------------

/**
 * @brief Convert StateSpace system to first-order DSP coefficients
 *
 * Assumes the system is first-order (1 state). If continuous-time, discretizes first.
 *
 * @param sys State-space system
 * @param Ts Sample time (only used if sys is continuous)
 * @return FirstOrderCoeffs\<T\> DSP coefficients
 */
template<typename T = float>
[[nodiscard]] constexpr FirstOrderCoeffs<T> to_coeffs(const StateSpace<1, 1, 1, T>& sys, T Ts = T{0}) {
    StateSpace<1, 1, 1, T> sys_d = sys;
    if (sys.Ts == T{0}) {
        if (!(Ts > T{0})) {
            return ::damp::detail::refused<FirstOrderCoeffs<T>>(DesignStatus::invalid_sample_time);
        }
        auto opt = discretize(sys, Ts, DiscretizationMethod::Tustin);
        if (!opt) {
            return ::damp::detail::refused<FirstOrderCoeffs<T>>(DesignStatus::singular);
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
 * @return FirstOrderCoeffs\<T\> DSP coefficients
 */
template<typename T = float>
[[nodiscard]] constexpr FirstOrderCoeffs<T> to_coeffs(const TransferFunction<1, 2, T>& tf, T Ts) {
    if (!(Ts > T{0})) {
        return ::damp::detail::refused<FirstOrderCoeffs<T>>(DesignStatus::invalid_sample_time);
    }
    const auto sys_c = tf.to_state_space();
    if (!sys_c) {
        return ::damp::detail::refused<FirstOrderCoeffs<T>>(DesignStatus::singular);
    }
    const auto sys_d = discretize(*sys_c, Ts, DiscretizationMethod::Tustin);
    if (!sys_d) {
        return ::damp::detail::refused<FirstOrderCoeffs<T>>(DesignStatus::singular);
    }
    return to_coeffs(*sys_d);
}

/**
 * @brief Convert a proper first-order TransferFunction to DSP coefficients (Tustin)
 *
 * @param tf Transfer function (highest power of s first)
 * @param Ts Sample time [s]
 * @return FirstOrderCoeffs\<T\> DSP coefficients; refused if @p Ts ≤ 0 (invalid_sample_time) or the realization fails (singular)
 */
template<typename T = float>
[[nodiscard]] constexpr FirstOrderCoeffs<T> to_coeffs(const TransferFunction<2, 2, T>& tf, T Ts) {
    if (!(Ts > T{0})) {
        return ::damp::detail::refused<FirstOrderCoeffs<T>>(DesignStatus::invalid_sample_time);
    }
    const auto sys_c = tf.to_state_space();
    if (!sys_c) {
        return ::damp::detail::refused<FirstOrderCoeffs<T>>(DesignStatus::singular);
    }
    const auto sys_d = discretize(*sys_c, Ts, DiscretizationMethod::Tustin);
    if (!sys_d) {
        return ::damp::detail::refused<FirstOrderCoeffs<T>>(DesignStatus::singular);
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
 * @return SecondOrderCoeffs\<T\> DSP coefficients
 */
template<typename T = float>
[[nodiscard]] constexpr SecondOrderCoeffs<T> to_coeffs(const StateSpace<2, 1, 1, T>& sys, T Ts = T{0}) {
    StateSpace<2, 1, 1, T> sys_d = sys;
    if (sys.Ts == T{0}) {
        if (!(Ts > T{0})) {
            return ::damp::detail::refused<SecondOrderCoeffs<T>>(DesignStatus::invalid_sample_time);
        }
        auto opt = discretize(sys, Ts, DiscretizationMethod::Tustin);
        if (!opt) {
            return ::damp::detail::refused<SecondOrderCoeffs<T>>(DesignStatus::singular);
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
 * @return SecondOrderCoeffs\<T\> DSP coefficients
 */
template<typename T = float>
[[nodiscard]] constexpr SecondOrderCoeffs<T> to_coeffs(const TransferFunction<3, 3, T>& tf, T Ts) {
    if (!(Ts > T{0})) {
        return ::damp::detail::refused<SecondOrderCoeffs<T>>(DesignStatus::invalid_sample_time);
    }
    const auto sys_c = tf.to_state_space();
    if (!sys_c) {
        return ::damp::detail::refused<SecondOrderCoeffs<T>>(DesignStatus::singular);
    }
    const auto sys_d = discretize(*sys_c, Ts, DiscretizationMethod::Tustin);
    if (!sys_d) {
        return ::damp::detail::refused<SecondOrderCoeffs<T>>(DesignStatus::singular);
    }
    return to_coeffs(*sys_d);
}

// ---- Shared helpers --------------------------------------------------------

namespace detail {
/**
 * @brief True when @p f [Hz] is strictly inside (0, Nyquist) for sample time @p Ts.
 *
 * Nyquist = 1/(2·Ts). Invalid designs are refused (no division by zero / NaN poles).
 */
template<typename T>
[[nodiscard]] constexpr bool valid_cutoff(T f, T Ts) {
    return (Ts > T{0}) && (f > T{0}) && ((T{2} * f * Ts) < T{1});
}

/** @brief Why a discrete design was refused: invalid_sample_time if @p Ts is the cause, else invalid_parameter. */
template<typename T>
[[nodiscard]] constexpr DesignStatus cutoff_status(T Ts) {
    return ((Ts > T{0}) && damp::isfinite(Ts)) ? DesignStatus::invalid_parameter : DesignStatus::invalid_sample_time;
}

/** @brief A model refused for an out-of-domain parameter (fc ≤ 0, Q ≤ 0, ζ ≤ 0, delay < 0). */
template<typename Result>
[[nodiscard]] constexpr Result invalid() {
    return ::damp::detail::refused<Result>(DesignStatus::invalid_parameter);
}

template<typename T>
[[nodiscard]] constexpr bool valid_Q(T Q) {
    return Q > T{0};
}

/**
 * @brief Prewarped bilinear gain k = ω / tan(ω·Ts/2), ω = 2π·f
 *
 * With s ← k·(1 − z⁻¹)/(1 + z⁻¹) the analog response at ω lands exactly on the
 * digital frequency ω (plain Tustin, k = 2/Ts, compresses it toward DC). Requires
 * valid_cutoff(f, Ts).
 */
template<typename T>
[[nodiscard]] constexpr T prewarped_bilinear_gain(T f, T Ts) {
    const T omega = T{2} * damp::numbers::pi_v<T> * f;
    return omega / damp::tan(omega * Ts / T{2});
}

/**
 * @brief Analog design frequency that Tustin (k = 2/Ts) maps onto digital @p f
 *
 * f_a = tan(π·f·Ts)/(π·Ts). Requires valid_cutoff(f, Ts).
 */
template<typename T>
[[nodiscard]] constexpr T prewarped_frequency(T f, T Ts) {
    const T pi_Ts = damp::numbers::pi_v<T> * Ts;
    return damp::tan(pi_Ts * f) / pi_Ts;
}

// Biquad (second-order IIR) designs from the RBJ audio-EQ cookbook formulas: the discrete
// notch, band-pass and second-order high-pass, and the EQ section at the end.
//
// Coefficients use the normalized difference equation
//   y[n] = b0·x[n] + b1·x[n-1] + b2·x[n-2] − a1·y[n-1] − a2·y[n-2]
// matching SecondOrderCoeffs and the Biquad runtime. Frequencies are designed
// directly in the digital domain: ω₀ = 2π·f·Ts. Quality factor Q controls
// bandwidth (BW in octaves ≈ asinh(1/2Q)·2/ln2 near ω₀); larger Q is narrower.
template<typename T>
[[nodiscard]] constexpr SecondOrderCoeffs<T> normalize_biquad(T b0, T b1, T b2, T a0, T a1, T a2) {
    // Every caller's a0 is a sum of positive terms (> 0 for valid f/Q/gain); guard
    // only the degenerate zero — return a zero section, not ±∞ coeffs.
    if (!(damp::abs(a0) > T{0})) {
        return ::damp::detail::refused<SecondOrderCoeffs<T>>(DesignStatus::singular);
    }
    const T inv = T{1} / a0;
    return SecondOrderCoeffs<T>{b0 * inv, b1 * inv, b2 * inv, a1 * inv, a2 * inv};
}
} // namespace detail

// ---- First-order low-pass --------------------------------------------------

/**
 * @brief First-order low-pass filter design (continuous-time)
 *
 * H(s) = 1 / (τs + 1) with τ = 1/(2π fc) (unit DC gain, pole at −1/τ).
 * Coefficients are highest power of s first.
 * A cutoff fc ≤ 0 is refused with DesignStatus::invalid_parameter.
 *
 * @param fc Cutoff frequency [Hz] (must be > 0)
 * @return TransferFunctionResult<2, 2, T> with the continuous-time transfer function; refused (invalid_parameter) for an out-of-domain input
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunctionResult<2, 2, T> lowpass_1st(T fc) {
    if (!(fc > T{0})) {
        return detail::invalid<TransferFunctionResult<2, 2, T>>();
    }
    const T tau = T{1} / (T{2} * damp::numbers::pi_v<T> * fc);
    // num = 1, den = τs + 1  →  H = 1/(τs+1)
    return detail::designed(TransferFunction<2, 2, T>{{T{0}, T{1}}, {tau, T{1}}});
}

/**
 * @brief First-order low-pass filter design
 *
 * Bilinear (Tustin) discretization of H(s) = ωc/(s + ωc), prewarped so the
 * −3 dB corner sits exactly at @p fc. Returns DSP coefficients ready for runtime use.
 *
 * @param fc Cutoff frequency [Hz] (must satisfy 0 < fc < 1/(2·Ts))
 * @param Ts Sample time [s] (must be > 0)
 * @return FirstOrderCoeffs\<T\> DSP coefficients; refused (invalid_parameter / invalid_sample_time) if @p fc / @p Ts invalid
 * @note Compare with MATLAB®'s butter(1, 2*fc*Ts).
 */
template<typename T = float>
[[nodiscard]] constexpr FirstOrderCoeffs<T> lowpass_1st(T fc, T Ts) {
    if (!detail::valid_cutoff(fc, Ts)) {
        return ::damp::detail::refused<FirstOrderCoeffs<T>>(detail::cutoff_status(Ts));
    }
    // Bilinear: s ← k (1 − z⁻¹)/(1 + z⁻¹), k prewarped to ωc.
    // Pole gives a1; numerator is forced to unit DC by construction so
    // (b0 + b1) = (1 + a1) holds under -ffast-math reassociation (same pattern
    // as lowpass_2nd).
    const T omega_c = T{2} * damp::numbers::pi_v<T> * fc;
    const T k = detail::prewarped_bilinear_gain(fc, Ts);
    const T denom = omega_c + k;

    const T a1 = (omega_c - k) / denom;
    const T dc_sum = T{1} + a1; // = b0 + b1 for H(1) = 1
    return FirstOrderCoeffs<T>{
        .b0 = dc_sum / T{2},
        .b1 = dc_sum / T{2},
        .a1 = a1,
    };
}

// ---- First-order high-pass -------------------------------------------------

/**
 * @brief First-order high-pass filter design (continuous-time)
 *
 * H(s) = s / (s + ωc) with ωc = 2π fc (DC gain 0, high-frequency gain 1).
 * Coefficients are highest power of s first.
 * A cutoff fc ≤ 0 is refused with DesignStatus::invalid_parameter.
 *
 * @param fc Cutoff frequency [Hz] (must be > 0)
 * @return TransferFunctionResult<2, 2, T> with the continuous-time transfer function; refused (invalid_parameter) for an out-of-domain input
 * @note Compare with MATLAB®'s butter(1, 2*pi*fc, 'high', 's').
 * @see highpass_1st(fc, Ts) for the discrete Tustin coefficient form
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunctionResult<2, 2, T> highpass_1st(T fc) {
    if (!(fc > T{0})) {
        return detail::invalid<TransferFunctionResult<2, 2, T>>();
    }
    const T omega_c = T{2} * damp::numbers::pi_v<T> * fc;
    // num = s, den = s + ωc  →  H = s/(s+ωc)
    return detail::designed(TransferFunction<2, 2, T>{{T{1}, T{0}}, {T{1}, omega_c}});
}

/**
 * @brief First-order high-pass filter design (Tustin / bilinear)
 *
 * Discrete form of H(s) = s / (s + ωc), ωc = 2π·fc, prewarped so the −3 dB
 * corner sits exactly at @p fc. Pair with @ref FirstOrderSection or a first-order
 * difference equation. Invalid @p fc / @p Ts are refused.
 *
 * @param fc Cutoff frequency [Hz] (0 < fc < Nyquist)
 * @param Ts Sample time [s] (> 0)
 * @return FirstOrderCoeffs\<T\> DSP coefficients
 * @note Compare with MATLAB®'s butter(1, 2*fc*Ts, 'high').
 * @see lowpass_1st()
 */
template<typename T = float>
[[nodiscard]] constexpr FirstOrderCoeffs<T> highpass_1st(T fc, T Ts) {
    if (!detail::valid_cutoff(fc, Ts)) {
        return ::damp::detail::refused<FirstOrderCoeffs<T>>(detail::cutoff_status(Ts));
    }
    // Prewarped bilinear of H(s) = s/(s+ωc). Pin a1 from the pole; force exact
    // DC null (b1 = −b0) and unit Nyquist gain H(−1) = 1 by construction:
    //   (b0 − b1)/(1 − a1) = 1 with b1 = −b0  ⇒  b0 = (1 − a1)/2.
    const T omega_c = T{2} * damp::numbers::pi_v<T> * fc;
    const T k = detail::prewarped_bilinear_gain(fc, Ts);
    const T denom = omega_c + k;
    const T a1 = (omega_c - k) / denom;
    const T b0 = (T{1} - a1) / T{2};
    return FirstOrderCoeffs<T>{
        .b0 = b0,
        .b1 = -b0,
        .a1 = a1,
    };
}

// ---- Second-order low-pass -------------------------------------------------

/**
 * @brief Second-order low-pass filter design (continuous-time)
 *
 * H(s) = ω₀² / (s² + (ω₀/Q) s + ω₀²), ω₀ = 2π fc (damping ratio ζ = 1/(2Q)).
 * Coefficients are highest power of s first (den[0] = 1: monic s²).
 *
 * A trailing Ts, `lowpass_2nd(fc, Q, Ts)`, is the discrete (Tustin) coefficient form.
 *
 * @param fc Cutoff frequency [Hz]
 * @param Q  Quality factor (default 1/√2, 2nd-order Butterworth)
 * @return TransferFunctionResult<3, 3, T> with the continuous-time transfer function; refused (invalid_parameter) for an out-of-domain input
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunctionResult<3, 3, T> lowpass_2nd(T fc, T Q = damp::numbers::inv_sqrt2_v<T>) {
    if (!(fc > T{0}) || !detail::valid_Q(Q)) {
        return detail::invalid<TransferFunctionResult<3, 3, T>>();
    }
    const T omega_0 = T{2} * damp::numbers::pi_v<T> * fc;
    const T omega_0_sq = omega_0 * omega_0;
    const T omega_over_Q = omega_0 / Q;

    return detail::designed(TransferFunction<3, 3, T>{
        {T{0}, T{0}, omega_0_sq},        // num: ω₀²
        {T{1}, omega_over_Q, omega_0_sq} // den: s² + (ω₀/Q) s + ω₀²
    });
}

/**
 * @brief Second-order low-pass filter design
 *
 * Bilinear discretization of H(s) = ω₀² / (s² + (ω₀/Q)s + ω₀²), prewarped at
 * @p fc so the digital response at fc equals the analog one (−3 dB for the
 * Butterworth Q = 1/√2). Equivalent to the RBJ cookbook low-pass.
 *
 * @param fc Cutoff frequency [Hz] (0 < fc < Nyquist)
 * @param Q  Quality factor (1/√2 = Butterworth); must be > 0
 * @param Ts Sample time [s] (> 0)
 * @return SecondOrderCoeffs\<T\> DSP coefficients; refused (invalid_parameter / invalid_sample_time) if arguments invalid
 */
template<typename T = float>
[[nodiscard]] constexpr SecondOrderCoeffs<T> lowpass_2nd(T fc, T Q, T Ts) {
    if (!detail::valid_cutoff(fc, Ts) || !detail::valid_Q(Q)) {
        return ::damp::detail::refused<SecondOrderCoeffs<T>>(detail::cutoff_status(Ts));
    }
    // Discretize H(s) = ω₀² / (s² + (ω₀/Q)s + ω₀²) with the bilinear transform
    // s ← k·(1 − z⁻¹)/(1 + z⁻¹), k prewarped to ω₀. Clearing (1 + z⁻¹)² gives a
    // numerator proportional to (1 + z⁻¹)² (tap shape 1 : 2 : 1) over a
    // denominator with a0 = k² + (ω₀/Q)k + ω₀².
    const T omega_0 = T{2} * damp::numbers::pi_v<T> * fc;
    const T k = detail::prewarped_bilinear_gain(fc, Ts);
    const T k_sq = k * k;
    const T omega_0_sq = omega_0 * omega_0;
    const T omega_over_Q = omega_0 / Q;

    const T a0 = k_sq + omega_over_Q * k + omega_0_sq;

    const T a1 = (T{2} * omega_0_sq - T{2} * k_sq) / a0;
    const T a2 = (k_sq - omega_over_Q * k + omega_0_sq) / a0;

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

// ---- Second-order high-pass ------------------------------------------------

/**
 * @brief Second-order high-pass filter design (continuous-time)
 *
 * H(s) = s² / (s² + (ω₀/Q) s + ω₀²), ω₀ = 2π fc.
 * Coefficients are highest power of s first (den[0] = 1: monic s²).
 *
 * A trailing Ts, `highpass_2nd(fc, Q, Ts)`, is the discrete (RBJ) coefficient form.
 *
 * @param fc Cutoff frequency [Hz]
 * @param Q  Quality factor (default 1/√2, 2nd-order Butterworth)
 * @return TransferFunctionResult<3, 3, T> with the continuous-time transfer function; refused (invalid_parameter) for an out-of-domain input
 * @note Compare with MATLAB®'s butter(2, 2*pi*fc, 'high', 's') (Q = 1/√2).
 * @see highpass_2nd(fc, Q, Ts) for the discrete RBJ coefficient form
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunctionResult<3, 3, T> highpass_2nd(T fc, T Q = damp::numbers::inv_sqrt2_v<T>) {
    if (!(fc > T{0}) || !detail::valid_Q(Q)) {
        return detail::invalid<TransferFunctionResult<3, 3, T>>();
    }
    const T omega_0 = T{2} * damp::numbers::pi_v<T> * fc;
    const T omega_0_sq = omega_0 * omega_0;
    const T omega_over_Q = omega_0 / Q;

    return detail::designed(TransferFunction<3, 3, T>{
        {T{1}, T{0}, T{0}},              // num: s²
        {T{1}, omega_over_Q, omega_0_sq} // den: s² + (ω₀/Q) s + ω₀²
    });
}

/**
 * @brief Second-order high-pass filter (RBJ).
 *
 * Counterpart to lowpass_2nd, parameterized by Q
 * (= 1 / 2ζ); the default Q = 1/√2 is the maximally-flat (Butterworth) response.
 *
 * @param fc Cutoff frequency [Hz]
 * @param Q  Quality factor (1/√2 = Butterworth)
 * @param Ts Sample time [s]
 * @return SecondOrderCoeffs\<T\> normalized biquad coefficients
 * @see lowpass_2nd()
 */
template<typename T = float>
[[nodiscard]] constexpr SecondOrderCoeffs<T> highpass_2nd(T fc, T Q, T Ts) {
    if (!detail::valid_cutoff(fc, Ts) || !detail::valid_Q(Q)) {
        return ::damp::detail::refused<SecondOrderCoeffs<T>>(detail::cutoff_status(Ts));
    }
    const T w0 = T{2} * damp::numbers::pi_v<T> * fc * Ts;
    const auto [s, c] = damp::sincos(w0);
    const T alpha = s / (T{2} * Q);
    // RBJ highpass: num taps ∝ (1+cos) · {1, −2, 1} → exact DC null by shape.
    // Normalize a0; DC identity H(1) = 0 holds independent of -ffast-math.
    const T b0 = (T{1} + c) / T{2};
    return detail::normalize_biquad<T>(b0, -(T{1} + c), b0, T{1} + alpha, T{-2} * c, T{1} - alpha);
}

// ---- Band-pass -------------------------------------------------------------

/**
 * @brief Second-order band-pass filter design (continuous-time)
 *
 * H(s) = (ω₀/Q) s / (s² + (ω₀/Q) s + ω₀²), ω₀ = 2π f0 (unit gain at ω₀).
 * Coefficients are highest power of s first (den[0] = 1: monic s²).
 * A trailing Ts, `bandpass(f0, Q, Ts)`, is the discrete (RBJ) coefficient form.
 *
 * @param f0 Center frequency [Hz]
 * @param Q  Quality factor (higher = narrower band)
 * @return TransferFunctionResult<3, 3, T> with the continuous-time transfer function; refused (invalid_parameter) for an out-of-domain input
 * @note Compare with MATLAB®'s iirpeak analog prototype (constant 0 dB peak).
 * @see bandpass(f0, Q, Ts) for the discrete RBJ coefficient form
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunctionResult<3, 3, T> bandpass(T f0, T Q) {
    if (!(f0 > T{0}) || !detail::valid_Q(Q)) {
        return detail::invalid<TransferFunctionResult<3, 3, T>>();
    }
    const T omega_0 = T{2} * damp::numbers::pi_v<T> * f0;
    const T omega_0_sq = omega_0 * omega_0;
    const T omega_over_Q = omega_0 / Q;

    return detail::designed(TransferFunction<3, 3, T>{
        {T{0}, omega_over_Q, T{0}},      // num: (ω₀/Q) s
        {T{1}, omega_over_Q, omega_0_sq} // den: s² + (ω₀/Q) s + ω₀²
    });
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
        return ::damp::detail::refused<SecondOrderCoeffs<T>>(detail::cutoff_status(Ts));
    }
    const T w0 = T{2} * damp::numbers::pi_v<T> * f0 * Ts;
    const auto [s, c] = damp::sincos(w0);
    const T alpha = s / (T{2} * Q);
    return detail::normalize_biquad<T>(alpha, T{0}, -alpha, T{1} + alpha, T{-2} * c, T{1} - alpha);
}

// ---- Notch (band-reject) ---------------------------------------------------

/**
 * @brief Second-order band-reject (notch) filter design (continuous-time)
 *
 * H(s) = (s² + ω₀²) / (s² + (ω₀/Q) s + ω₀²), ω₀ = 2π f0 (null at ω₀, unit DC/HF).
 * Coefficients are highest power of s first (den[0] = 1: monic s²).
 * A trailing Ts, `notch(f0, Q, Ts)`, is the discrete (RBJ) coefficient form.
 *
 * @param f0 Notch (center) frequency [Hz]
 * @param Q  Quality factor (higher = narrower notch)
 * @return TransferFunctionResult<3, 3, T> with the continuous-time transfer function; refused (invalid_parameter) for an out-of-domain input
 * @note Compare with MATLAB®'s iirnotch analog prototype.
 * @see notch(f0, Q, Ts) for the discrete RBJ coefficient form
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunctionResult<3, 3, T> notch(T f0, T Q) {
    if (!(f0 > T{0}) || !detail::valid_Q(Q)) {
        return detail::invalid<TransferFunctionResult<3, 3, T>>();
    }
    const T omega_0 = T{2} * damp::numbers::pi_v<T> * f0;
    const T omega_0_sq = omega_0 * omega_0;
    const T omega_over_Q = omega_0 / Q;

    return detail::designed(TransferFunction<3, 3, T>{
        {T{1}, T{0}, omega_0_sq},        // num: s² + ω₀²
        {T{1}, omega_over_Q, omega_0_sq} // den: s² + (ω₀/Q) s + ω₀²
    });
}

/**
 * @brief Second-order band-reject (notch) filter.
 *
 * Rejects a narrow band around f0 (gain → 0 at f0) while passing the rest
 * (gain → 1). Transfer function:
 *
 *     H(z) = (1 − 2cosω₀ z⁻¹ + z⁻²) / ((1 + α) − 2cosω₀ z⁻¹ + (1 − α) z⁻²),
 *     ω₀ = 2π·f0·Ts,  α = sinω₀ / (2Q)
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
        return ::damp::detail::refused<SecondOrderCoeffs<T>>(detail::cutoff_status(Ts));
    }
    const T w0 = T{2} * damp::numbers::pi_v<T> * f0 * Ts;
    const auto [s, c] = damp::sincos(w0);
    const T alpha = s / (T{2} * Q);
    return detail::normalize_biquad<T>(T{1}, T{-2} * c, T{1}, T{1} + alpha, T{-2} * c, T{1} - alpha);
}

// ---- All-pass --------------------------------------------------------------

/**
 * @brief Second-order all-pass filter design (continuous-time)
 *
 * H(s) = (s² − (ω₀/Q) s + ω₀²) / (s² + (ω₀/Q) s + ω₀²), ω₀ = 2π f0.
 * |H(jω)| = 1 for all ω; phase lags through 180° around ω₀ (Q sets the slope).
 * Coefficients are highest power of s first (den[0] = 1: monic s²).
 *
 * @param f0 Center frequency [Hz]
 * @param Q  Quality factor (higher = steeper phase around ω₀)
 * @return TransferFunctionResult<3, 3, T> with the continuous-time transfer function; refused (invalid_parameter) for an out-of-domain input
 * @note Compare with MATLAB®'s tf([1, −ω₀/Q, ω₀²], [1, ω₀/Q, ω₀²]).
 * @see allpass(f0, Q, Ts) for the discrete RBJ coefficient form
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunctionResult<3, 3, T> allpass(T f0, T Q) {
    if (!(f0 > T{0}) || !detail::valid_Q(Q)) {
        return detail::invalid<TransferFunctionResult<3, 3, T>>();
    }
    const T omega_0 = T{2} * damp::numbers::pi_v<T> * f0;
    const T omega_0_sq = omega_0 * omega_0;
    const T omega_over_Q = omega_0 / Q;

    return detail::designed(TransferFunction<3, 3, T>{
        {T{1}, -omega_over_Q, omega_0_sq}, // num: s² − (ω₀/Q) s + ω₀²
        {T{1}, omega_over_Q, omega_0_sq}   // den: s² + (ω₀/Q) s + ω₀²
    });
}

/**
 * @brief Second-order all-pass filter (RBJ cookbook), phase −180° at f0
 *
 * RBJ all-pass at ω₀ = 2π f0 Ts, α = sin ω₀ / (2Q): the bilinear map of allpass(f0, Q)
 * prewarped at f0, so the phase is −180° exactly at @p f0. The numerator is the normalized
 * denominator mirrored (b0 = a2, b1 = a1, b2 = 1), so |H| = 1 at every frequency holds for
 * the rounded coefficients too.
 *
 * @param f0 Center frequency [Hz] (0 < f0 < Nyquist)
 * @param Q  Quality factor (higher = steeper phase around f0)
 * @param Ts Sample time [s] (> 0)
 * @return SecondOrderCoeffs\<T\> DSP coefficients; refused (invalid_parameter / invalid_sample_time) if arguments invalid
 * @see allpass(f0, Q) for the continuous prototype
 */
template<typename T = float>
[[nodiscard]] constexpr SecondOrderCoeffs<T> allpass(T f0, T Q, T Ts) {
    if (!detail::valid_cutoff(f0, Ts) || !detail::valid_Q(Q)) {
        return ::damp::detail::refused<SecondOrderCoeffs<T>>(detail::cutoff_status(Ts));
    }
    const auto [s, c] = damp::sincos(T{2} * damp::numbers::pi_v<T> * f0 * Ts);
    const T alpha = s / (T{2} * Q);
    const T a1 = (T{-2} * c) / (T{1} + alpha);
    const T a2 = (T{1} - alpha) / (T{1} + alpha);
    return SecondOrderCoeffs<T>{.b0 = a2, .b1 = a1, .b2 = T{1}, .a1 = a1, .a2 = a2};
}

// ---- Butterworth low-pass --------------------------------------------------

/**
 * @brief Butterworth low-pass filter design
 *
 * Continuous Butterworth of order 1–4 at cutoff fc [Hz], built as cascade of
 * first- and second-order factors of H(s) = 1 / B_n(s/ωc) with the standard
 * quadratic quality factors Q = 1/(2ζ):
 * - n=1: single real pole
 * - n=2: Q = 1/(2 sin(π/4)) = 1/√2
 * - n=3: real pole × quadratic with Q = 1/(2 sin(π/6)) = 1
 * - n=4: quadratics with Q = 1/(2 sin(π/8)), 1/(2 sin(3π/8))
 *
 * @param fc Cutoff frequency [Hz] (> 0)
 * @return StateSpaceResult with the continuous filter; refused (invalid_parameter) if fc ≤ 0
 */
template<size_t Order, typename T = double>
    requires(Order >= 1 && Order <= 4)
[[nodiscard]] constexpr StateSpaceResult<Order, 1, 1, T> butterworth_lowpass(T fc) {
    using Result = StateSpaceResult<Order, 1, 1, T>;
    if (!(fc > T{0}) || !damp::isfinite(fc)) {
        return detail::invalid<Result>();
    }
    damp::optional<StateSpace<Order, 1, 1, T>> opt{};
    if constexpr (Order == 1) {
        opt = lowpass_1st<T>(fc).tf.to_state_space();
    } else if constexpr (Order == 2) {
        // Butterworth Q = 1/√2
        opt = lowpass_2nd(fc, damp::numbers::inv_sqrt2_v<T>).tf.to_state_space();
    } else if constexpr (Order == 3) {
        // B3(s) ∝ (s+1)(s² + s + 1) → Q = 1 on the quadratic at the same ωc
        opt = (lowpass_1st<T>(fc).tf * lowpass_2nd(fc, T{1}).tf).to_state_space();
    } else { // Order == 4
        // B4 factors: Q_k = 1/(2 sin((2k−1)π / 8)), k = 1, 2
        const T pi = damp::numbers::pi_v<T>;
        const T q1 = T{1} / (T{2} * damp::sin(pi / T{8}));
        const T q2 = T{1} / (T{2} * damp::sin(T{3} * pi / T{8}));
        opt = (lowpass_2nd(fc, q1).tf * lowpass_2nd(fc, q2).tf).to_state_space();
    }
    if (!opt) {
        return ::damp::detail::refused<Result>(DesignStatus::singular);
    }
    return detail::designed(*opt);
}

/**
 * @brief Butterworth low-pass filter design (discrete-time)
 *
 * Creates a discrete-time Butterworth low-pass filter of specified order.
 * Butterworth filters have maximally flat magnitude response in passband.
 * The analog prototype is prewarped so the −3 dB corner lands exactly on @p fc
 * after the Tustin step.
 *
 * @param fc Cutoff frequency [Hz] (0 < fc < Nyquist)
 * @param Ts Sample time [s]
 * @return StateSpaceResult with the discrete filter; refused with invalid_sample_time if Ts ≤ 0,
 *         invalid_parameter if fc is outside (0, Nyquist), singular if the Tustin step fails
 * @note Compare with MATLAB®'s butter(Order, 2*fc*Ts).
 */
template<size_t Order, typename T = float>
    requires(Order >= 1 && Order <= 4)
[[nodiscard]] constexpr StateSpaceResult<Order, 1, 1, T> butterworth_lowpass(T fc, T Ts) {
    using Result = StateSpaceResult<Order, 1, 1, T>;
    if (!detail::valid_cutoff(fc, Ts)) {
        return ::damp::detail::refused<Result>(detail::cutoff_status(Ts));
    }
    const auto proto = butterworth_lowpass<Order, T>(detail::prewarped_frequency(fc, Ts));
    if (!proto.success) {
        return proto;
    }
    const auto sys_d = discretize(proto.sys, Ts, DiscretizationMethod::Tustin);
    if (!sys_d) {
        return ::damp::detail::refused<Result>(DesignStatus::singular);
    }
    return detail::designed(*sys_d);
}

// ---- Audio EQ (RBJ cookbook): peaking, low shelf, high shelf ---------------

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
        return ::damp::detail::refused<SecondOrderCoeffs<T>>(detail::cutoff_status(Ts));
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
        return ::damp::detail::refused<SecondOrderCoeffs<T>>(detail::cutoff_status(Ts));
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
        return ::damp::detail::refused<SecondOrderCoeffs<T>>(detail::cutoff_status(Ts));
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
