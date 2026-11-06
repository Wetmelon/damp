// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file pid_design.hpp
 * @brief PID design: pidtune, Ziegler–Nichols, pole placement → PIDGains
 *
 * @defgroup pid_design PID Design
 * @brief PID tuning rules, plant-aware loop shaping (pidtune) and PID pole placement
 *
 * The tuning rules use experimentally determined parameters (ultimate gain and
 * period, or a first-order-plus-dead-time step response); pidtune and
 * pid_pole_placement use a plant model. All functions return continuous PIDGains;
 * deploy via `.discretize(Ts)` → PIDController, or ContinuousPID for variable-rate.
 *
 * Every function reports failure through PIDGains::success / status: non-finite or
 * non-physical inputs (K = 0, L ≤ 0, tau ≤ 0, Ku ≤ 0, Tu ≤ 0, ωc ≤ 0, …) and modes the
 * cited rule does not define are refused with DesignStatus::invalid_parameter, and the
 * refused gains are zero, so no Inf/NaN reaches a runtime.
 *
 * Includes:
 * - Ziegler-Nichols (ultimate gain method and step response method)
 * - Tyreus-Luyben (conservative Ziegler-Nichols variant)
 * - AMIGO PI from the ultimate point and the static gain
 * - Cohen-Coon (first-order plus dead-time models)
 * - SIMC (Skogestad IMC) rules for FOPDT models
 * - Lambda tuning (setpoint tracking, robustness-focused)
 * - pidtune: crossover and phase margin on a state-space plant, with a closed-loop stability certificate
 * - Pole placement on a first- or second-order plant (continuous ω, ζ, ωᵢ, or s-plane poles sampled at Ts)
 */

#include "damp/backend.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/design/status.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/systems/zpk.hpp"

namespace damp {

namespace design {

namespace detail {

/** @brief Gains refused because a tuning input is outside its domain (K = 0, L ≤ 0, Tu ≤ 0, …). */
template<typename T>
[[nodiscard]] constexpr PIDGains<T> pid_invalid() {
    return ::damp::detail::refused<PIDGains<T>>(DesignStatus::invalid_parameter);
}

/** @brief x finite and > 0. The finiteness test is bit-pattern based, so it holds under -ffast-math. */
template<typename T>
[[nodiscard]] constexpr bool positive_finite(T x) {
    return damp::isfinite(x) && (x > T{0});
}

/** @brief A static gain K that is finite and nonzero (either sign). */
template<typename T>
[[nodiscard]] constexpr bool valid_static_gain(T K) {
    return damp::isfinite(K) && (K != T{0});
}

/** @brief A first-order-plus-dead-time model with K finite and nonzero, L > 0 and tau > 0. */
template<typename T>
[[nodiscard]] constexpr bool valid_fopdt(T K, T L, T tau) {
    return valid_static_gain(K) && positive_finite(L) && positive_finite(tau);
}

/**
 * @brief Back-calculation seed for a placed PI / PID: Kbc = sgn(Kᵢ)·max(|Kp|, |Kᵢ|·T_min)
 *
 * The tracking time T_t = Kbc/Kᵢ = max(|T_i|, T_min) is positive, so the correction opposes
 * windup even for a placement with Kp < 0, and T_min > 0 keeps Kbc nonzero as Kp → 0.
 * Ki = 0 (no integrator, back-calculation inert) gives |Kp|.
 */
template<typename T>
[[nodiscard]] constexpr T placed_kbc(T Kp, T Ki, T T_min = T{0}) {
    const T mag = damp::max(damp::abs(Kp), damp::abs(Ki) * T_min);
    return (Ki < T{0}) ? -mag : mag;
}

} // namespace detail

// ============================================================================
// Ziegler-Nichols Ultimate Gain Method
// ============================================================================

/**
 * @brief Ziegler-Nichols tuning from ultimate gain and ultimate period
 *
 * The user experimentally determines:
 * - Ku: ultimate (critical) gain where system oscillates continuously
 * - Tu: ultimate period of those oscillations
 *
 * Classic formulas from Ziegler & Nichols (1942).
 *
 * @param Ku   Ultimate gain (gain at which system marginally oscillates)
 * @param Tu   Ultimate period of oscillation (seconds)
 * @param mode Control law (P, PI, PD, or PID)
 * @return PIDGains with tuned gains; refused (invalid_parameter) if Ku ≤ 0, Tu ≤ 0 (PI/PD/PID) or an input is non-finite
 */
template<typename T = double>
[[nodiscard]] constexpr PIDGains<T>
ziegler_nichols(T Ku, T Tu, PIDMode mode = PIDMode::PID) {
    // PI/PD/PID divide by Tu; P-only uses Ku alone.
    if (!detail::positive_finite(Ku) || (mode != PIDMode::P && !detail::positive_finite(Tu))) {
        return detail::pid_invalid<T>();
    }

    T Kp{};
    T Ki{};
    T Kd{};

    switch (mode) {
        case PIDMode::P:
            Kp = static_cast<T>(0.5) * Ku;
            break;
        case PIDMode::PI:
            Kp = static_cast<T>(0.45) * Ku;
            Ki = Kp / (Tu / static_cast<T>(1.2));
            break;
        case PIDMode::PD:
            Kp = static_cast<T>(0.8) * Ku;
            Kd = Kp * Tu / T{8};
            break;
        case PIDMode::PID:
            Kp = static_cast<T>(0.6) * Ku;
            Ki = Kp / (Tu / T{2});
            Kd = Kp * Tu / T{8};
            break;
    }
    return PIDGains<T>{Kp, Ki, Kd};
}

/**
 * @brief Ziegler-Nichols step response method (reaction curve)
 *
 * Uses first-order-plus-dead-time (FOPDT) model parameters obtained
 * from a step response:
 * - K: static gain (output change / input change)
 * - L: apparent dead time (seconds)
 * - tau: time constant (seconds)
 *
 * @param K    Static gain
 * @param L    Apparent dead time (delay, seconds)
 * @param tau  Time constant (seconds)
 * @param mode Control law (P, PI, or PID; the reaction-curve table has no PD row)
 * @return PIDGains with tuned gains; refused (invalid_parameter) if K = 0, L ≤ 0, tau ≤ 0, an input
 *         is non-finite, or mode is PD
 *
 * @see "Optimum settings for automatic controllers" (Ziegler & Nichols, 1942)
 */
template<typename T = double>
[[nodiscard]] constexpr PIDGains<T>
ziegler_nichols_step(T K, T L, T tau, PIDMode mode = PIDMode::PID) {
    // a = tau/(K·L); Ki/Kd formulas also divide by L. tau ≤ 0 would flip the gain sign.
    if (!detail::valid_fopdt(K, L, tau) || mode == PIDMode::PD) {
        return detail::pid_invalid<T>();
    }

    T Kp{};
    T Ki{};
    T Kd{};

    T a = tau / (K * L); // normalized gain
    switch (mode) {
        case PIDMode::P:
        case PIDMode::PD: // refused above
            Kp = a;
            break;
        case PIDMode::PI:
            Kp = static_cast<T>(0.9) * a;
            Ki = Kp / (L * static_cast<T>(3.33));
            break;
        case PIDMode::PID:
            Kp = static_cast<T>(1.2) * a;
            Ki = Kp / (T{2} * L);
            Kd = Kp * static_cast<T>(0.5) * L;
            break;
    }
    return PIDGains<T>{Kp, Ki, Kd};
}

// ============================================================================
// Tyreus-Luyben (conservative Ziegler-Nichols variant)
// ============================================================================

/**
 * @brief Tyreus-Luyben tuning from ultimate gain and ultimate period
 *
 * A more conservative variant of Ziegler-Nichols that produces less aggressive
 * controllers with better robustness. Recommended for processes where
 * Ziegler-Nichols is too aggressive.
 *
 * @param Ku   Ultimate gain
 * @param Tu   Ultimate period (seconds)
 * @param mode Control law (PI or PID; the rule defines no P or PD setting)
 * @return PIDGains with tuned gains; refused (invalid_parameter) if Ku ≤ 0, Tu ≤ 0, an input is
 *         non-finite, or mode is P or PD
 *
 * @see "Tuning PI controllers for integrator/dead time processes" (Tyreus & Luyben, 1992),
 *      Ind. Eng. Chem. Res. 31(11)
 */
template<typename T = double>
[[nodiscard]] constexpr PIDGains<T>
tyreus_luyben(T Ku, T Tu, PIDMode mode = PIDMode::PID) {
    if (!detail::positive_finite(Ku) || !detail::positive_finite(Tu) || (mode != PIDMode::PI && mode != PIDMode::PID)) {
        return detail::pid_invalid<T>();
    }
    const T Kp = Ku / ((mode == PIDMode::PI) ? static_cast<T>(3.2) : static_cast<T>(2.2));
    const T Ki = Kp / (static_cast<T>(2.2) * Tu);
    const T Kd = (mode == PIDMode::PID) ? (Kp * Tu / static_cast<T>(6.3)) : T{0};
    return PIDGains<T>{Kp, Ki, Kd};
}

// ============================================================================
// AMIGO (from ultimate point + static gain)
// ============================================================================

/**
 * @brief AMIGO PI from ultimate gain/period and static gain Kₛ.
 *
 * Recovers an approximate FOPDT model from (Kₛ, Kᵤ, Tᵤ), then applies the
 * AMIGO PI rules (Åström & Hägglund, robust Ms ≈ 1.4 design). Prefer this over
 * Ziegler–Nichols when a biased/asymmetric relay has measured Kₛ.
 *
 * FOPDT recovery at ωᵤ = 2π/Tᵤ:
 *
 *     T = √((Kₛ Kᵤ)² − 1) / ωᵤ
 *     L = (π − atan(ωᵤ T)) / ωᵤ
 *
 * AMIGO PI (Advanced PID Control, §7.6):
 *
 *     Kp = 0.15/Kₛ + (0.35 − L T/(L + T)²) · T/(Kₛ L)
 *     Ti = 0.35 L + 13 L T²/(T² + 12 L T + 7 L²)
 *     Ki = Kp / Ti
 *
 * @param Ks Static process gain (from biased relay or step)
 * @param Ku Ultimate gain
 * @param Tu Ultimate period [s]
 * @return PIDGains PI gains (Kd = 0); refused (invalid_parameter) if inputs are non-physical
 *
 * @see "Advanced PID Control" (Åström & Hägglund, 2006), §7.6 (AMIGO PI)
 */
template<typename T = double>
[[nodiscard]] constexpr PIDGains<T>
amigo_kappa_tau(T Ks, T Ku, T Tu) {
    if (!detail::positive_finite(Ks) || !detail::positive_finite(Ku) || !detail::positive_finite(Tu)) {
        return detail::pid_invalid<T>();
    }
    const T wu = (T{2} * damp::numbers::pi_v<T>) / Tu;
    const T prod = Ks * Ku;
    if (prod <= T{1}) {
        // No real FOPDT consistent with this ultimate point.
        return detail::pid_invalid<T>();
    }
    const T tau = damp::sqrt((prod * prod) - T{1}) / wu;
    const T L = (damp::numbers::pi_v<T> - damp::atan(wu * tau)) / wu;
    if (!(L > T{0}) || !(tau > T{0})) {
        return detail::pid_invalid<T>();
    }
    const T LpT = L + tau;
    const T Kp = (static_cast<T>(0.15) / Ks) + ((static_cast<T>(0.35) - (L * tau / (LpT * LpT))) * tau / (Ks * L));
    const T Ti = (static_cast<T>(0.35) * L)
               + (T{13} * L * tau * tau / ((tau * tau) + (T{12} * L * tau) + (T{7} * L * L)));
    if (!(Ti > T{0}) || !(Kp > T{0})) {
        return detail::pid_invalid<T>();
    }
    const T Ki = Kp / Ti;
    return PIDGains<T>{Kp, Ki, T{0}};
}

// ============================================================================
// Cohen-Coon (FOPDT model)
// ============================================================================

/**
 * @brief Cohen-Coon tuning from first-order-plus-dead-time model
 *
 * Better than Ziegler-Nichols for processes with large dead time relative
 * to the time constant (L/tau > 0.25). Uses FOPDT model parameters.
 *
 * @param K    Static gain
 * @param L    Dead time (seconds)
 * @param tau  Time constant (seconds)
 * @param mode Control law (P, PI, PD, or PID)
 * @return PIDGains with tuned gains; refused (invalid_parameter) if K = 0, L ≤ 0, tau ≤ 0, an input
 *         is non-finite, or mode is PD with L/tau ≥ 3 (where the PD table's Td is not positive)
 *
 * @see "Theoretical consideration of retarded control" (Cohen & Coon, 1953)
 */
template<typename T = double>
[[nodiscard]] constexpr PIDGains<T>
cohen_coon(T K, T L, T tau, PIDMode mode = PIDMode::PID) {
    // r = L/tau and a = tau/(K·L) both require physical FOPDT parameters.
    if (!detail::valid_fopdt(K, L, tau)) {
        return detail::pid_invalid<T>();
    }

    T Kp{};
    T Ki{};
    T Kd{};

    T r = L / tau; // dead-time ratio
    T a = tau / (K * L);
    if (mode == PIDMode::PD && !(r < T{3})) {
        return detail::pid_invalid<T>(); // Td = L(6 − 2r)/(22 + 3r) ≤ 0
    }

    switch (mode) {
        case PIDMode::P:
            Kp = a * (T{1} + (r / T{3}));
            break;
        case PIDMode::PI: {
            Kp = a * (static_cast<T>(0.9) + (r / T{12}));
            T Ti = L * (T{30} + (T{3} * r)) / (T{9} + (T{20} * r));
            Ki = Kp / Ti;
        } break;
        case PIDMode::PD: {
            Kp = a * (static_cast<T>(1.25) + (r / T{6}));
            T Td = L * (T{6} - (T{2} * r)) / (T{22} + (T{3} * r));
            Kd = Kp * Td;
        } break;
        case PIDMode::PID: {
            Kp = a * ((T{4} / T{3}) + (r / T{4}));
            T Ti = L * (T{32} + (T{6} * r)) / (T{13} + (T{8} * r));
            T Td = L * T{4} / (T{11} + (T{2} * r));
            Ki = Kp / Ti;
            Kd = Kp * Td;
        } break;
    }
    return PIDGains<T>{Kp, Ki, Kd};
}

// ============================================================================
// SIMC (Skogestad IMC) Tuning Rules
// ============================================================================

/**
 * @brief SIMC (Skogestad Internal Model Control) tuning for FOPDT models
 *
 * The SIMC PI rule from Skogestad (2003): Kp = tau/(K·(tau_c + L)), Ti = min(tau, 4(tau_c + L)).
 * The user specifies a single tuning parameter tau_c (desired closed-loop time constant);
 * the SIMC default tau_c = L gives tight control with good robustness, and a larger tau_c is
 * slower and more robust. P mode is the PI gain without the integrator. The derivative of the
 * PD and PID modes, Kd = Kp·L/2, is not part of SIMC (which adds τ_D = τ₂ only for a
 * second-order plant); it is a dead-time heuristic and gives Kd = 0 at L = 0.
 *
 * @param K     Static gain
 * @param L     Dead time (seconds, ≥ 0)
 * @param tau   Time constant (seconds)
 * @param tau_c Desired closed-loop time constant (seconds). Larger = more robust.
 * @param mode  Control law (P, PI, PD or PID)
 * @return PIDGains with tuned gains; refused (invalid_parameter) if K = 0, L < 0, tau ≤ 0,
 *         (tau_c + L) ≤ 0, or an input is non-finite
 *
 * @see "Simple analytic rules for model reduction and PID controller tuning" (Skogestad, 2003),
 *      J. Process Control 13(4)
 */
template<typename T = double>
[[nodiscard]] constexpr PIDGains<T>
simc(T K, T L, T tau, T tau_c, PIDMode mode = PIDMode::PI) {
    // Kp = tau / (K·(tau_c+L)); PI/PID also need Ti > 0 for Ki = Kp/Ti.
    const T den = tau_c + L;
    if (!detail::valid_static_gain(K) || !damp::isfinite(L) || !(L >= T{0}) || !detail::positive_finite(tau)
        || !damp::isfinite(tau_c) || !(den > T{0})) {
        return detail::pid_invalid<T>();
    }

    T Kp{};
    T Ki{};
    T Kd{};

    switch (mode) {
        case PIDMode::P:
            Kp = tau / (K * den);
            break;
        case PIDMode::PI: {
            Kp = tau / (K * den);
            T Ti = (tau < T{4} * den) ? tau : T{4} * den;
            if (!(Ti > T{0})) {
                return detail::pid_invalid<T>();
            }
            Ki = Kp / Ti;
        } break;
        case PIDMode::PD: {
            Kp = tau / (K * den);
            Kd = Kp * L / T{2};
        } break;
        case PIDMode::PID: {
            Kp = tau / (K * den);
            T Ti = (tau < T{4} * den) ? tau : T{4} * den;
            if (!(Ti > T{0})) {
                return detail::pid_invalid<T>();
            }
            Ki = Kp / Ti;
            Kd = Kp * L / T{2};
        } break;
    }
    return PIDGains<T>{Kp, Ki, Kd};
}

// ============================================================================
// Lambda Tuning (IMC-based)
// ============================================================================

/**
 * @brief Lambda tuning for FOPDT model
 *
 * IMC-based approach where the user directly specifies the desired
 * closed-loop time constant (lambda). Produces non-oscillatory
 * setpoint tracking. Very robust, commonly used in process control.
 *
 * @param K      Static gain
 * @param L      Dead time (seconds)
 * @param tau    Time constant (seconds)
 * @param lambda Desired closed-loop time constant (seconds, > 0); λ ≥ L is the usual robust choice
 * @return PIDGains (PI controller gains); refused (invalid_parameter) if K = 0, L < 0, tau ≤ 0,
 *         λ ≤ 0, or an input is non-finite
 */
template<typename T = double>
[[nodiscard]] constexpr PIDGains<T>
lambda_tuning(T K, T L, T tau, T lambda) {
    const T den = lambda + L;
    if (!detail::valid_static_gain(K) || !damp::isfinite(L) || !(L >= T{0}) || !detail::positive_finite(tau)
        || !detail::positive_finite(lambda)) {
        return detail::pid_invalid<T>();
    }
    T Kp = tau / (K * den);
    T Ki = Kp / tau;
    return PIDGains<T>{Kp, Ki, T{0}};
}

/**
 * @brief Knobs for plant-aware @ref pidtune
 *
 * Crossover, phase margin, law, derivative filter Tf, and 2-DOF weights
 * b / c. Defaults match a 1-DOF PID at 60°.
 * When b or c is not 1, the loop crossover is raised so T (r→y) −3 dB matches
 * the 1-DOF design at wc.
 */
template<typename T = double>
struct PIDTuneSpec {
    T       wc{};                     ///< 1-DOF loop crossover ωc [rad/s]; 2-DOF matches T r→y −3 dB to this 1-DOF design
    T       phase_margin_deg = T{60}; ///< Phase margin at the loop crossover [deg], in (0, 180)
    PIDMode mode = PIDMode::PID;      ///< P / PI / PD / PID
    T       Tf = T{0};                ///< Derivative filter Tf [s]; 0 = off
    T       b = T{1};                 ///< P setpoint weight (0 = I-PD, 1 = on error)
    T       c = T{1};                 ///< D setpoint weight (0 = PI-D, 1 = on error)
};

namespace pidtune_detail {

template<typename T>
[[nodiscard]] constexpr bool weights_move_tracking(const PIDGains<T>& r) {
    const T    tol = T{64} * eps<T>(); // b, c are dimensionless weights
    const bool p_w = (damp::abs(r.b - T{1}) > tol) && (damp::abs(r.Kp) > T{0});
    const bool d_w = (damp::abs(r.c - T{1}) > tol) && (damp::abs(r.Kd) > T{0});
    if (!p_w && !d_w) {
        return false;
    }
    // I-P / I-PD with no integrator cannot hold a setpoint.
    if (!(damp::abs(r.Ki) > T{0}) && p_w && !(damp::abs(r.b) > tol)) {
        return false;
    }
    return true;
}

/** @brief The spec's scalar domain: ωc > 0, 0 < PM < 180°, Tf ≥ 0, all finite. */
template<typename T>
[[nodiscard]] constexpr bool valid_spec(const PIDTuneSpec<T>& spec) {
    return detail::positive_finite(spec.wc) && detail::positive_finite(spec.phase_margin_deg)
        && (spec.phase_margin_deg < T{180}) && damp::isfinite(spec.Tf) && (spec.Tf >= T{0}) && damp::isfinite(spec.b)
        && damp::isfinite(spec.c);
}

template<size_t NX, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr damp::optional<T>
tracking_mag(const damp::StateSpace<NX, 1, 1, T, NW, NV>& sys, const PIDGains<T>& pid, T w) {
    using Cplx = damp::complex<T>;
    if (!(w > T{0})) {
        return damp::nullopt;
    }
    const auto Gopt = damp::eval_frf(sys, Cplx{T{0}, w});
    if (!Gopt) {
        return damp::nullopt;
    }
    const Cplx P = (*Gopt)(0, 0);
    const Cplx L = pid.eval_fb(w) * P;
    const Cplx den = Cplx{T{1}, T{0}} + L;
    if (!(den.norm() > T{0})) { // closed-loop pole exactly on jω
        return damp::nullopt;
    }
    return damp::abs((pid.eval_ff(w) * P) / den);
}

template<size_t NX, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr damp::optional<T>
tracking_bandwidth(const damp::StateSpace<NX, 1, 1, T, NW, NV>& sys, const PIDGains<T>& pid, T w_hint) {
    const T    w_ref = (w_hint > T{0}) ? w_hint : T{1};
    const auto mag_dc = tracking_mag(sys, pid, w_ref * static_cast<T>(1e-3));
    if (!mag_dc || !(*mag_dc > T{0})) {
        return damp::nullopt;
    }
    const T thresh = *mag_dc / damp::sqrt(T{2});
    T       lo = w_ref * static_cast<T>(1e-3);
    T       hi = w_ref * static_cast<T>(1e3);
    auto    mag_hi = tracking_mag(sys, pid, hi);
    for (int k = 0; k < 8 && mag_hi && (*mag_hi >= thresh); ++k) {
        hi *= T{4};
        mag_hi = tracking_mag(sys, pid, hi);
        if (hi > (w_ref * static_cast<T>(1e6))) {
            return damp::nullopt;
        }
    }
    if (!mag_hi || (*mag_hi >= thresh)) {
        return damp::nullopt;
    }
    for (int k = 0; k < 16; ++k) {
        const T    mid = damp::sqrt(lo * hi);
        const auto m = tracking_mag(sys, pid, mid);
        if (!m || (*m >= thresh)) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return damp::sqrt(lo * hi);
}

/**
 * @brief The requested law with C(jω) = e^{j(−π + PM)}/P(jω), from the plant response P at ω = spec.wc
 *
 * The filtered derivative Kd·g, g = jω/(1 + jωTf), enters the fit, so the delivered C_fb meets
 * the crossover and the phase margin for any Tf. With g = g_r + j g_i:
 * PI: Kp = C_r, Kᵢ = −C_i ω. PD: Kd = C_i/g_i, Kp = C_r − Kd g_r. PID (Tᵢ = 4T_d, Kᵢ = Kp²/(4Kd)):
 * (4ω g_i − g_r²)Kd² + (2C_r g_r − 4ω C_i)Kd − C_r² = 0, Kp = C_r − Kd g_r. P meets the crossover
 * alone: Kp = 1/|P|. Refused with singular if P is zero or non-finite, and with not_placeable if
 * the law cannot realize C with positive gains (a PI cannot lead, a PD cannot lag, no law gives
 * Kp ≤ 0).
 */
template<typename T>
[[nodiscard]] constexpr PIDGains<T> fit_to_response(damp::complex<T> P, const PIDTuneSpec<T>& spec) noexcept {
    constexpr T pi = damp::numbers::pi_v<T>;
    const T     mag_P = damp::abs(P);
    if (!detail::positive_finite(mag_P)) {
        return ::damp::detail::refused<PIDGains<T>>(DesignStatus::singular);
    }
    const T want = -pi + (spec.phase_margin_deg * pi / T{180}) - damp::arg(P);
    const auto [s, c] = damp::sincos(want);
    const T     Cr = c / mag_P;
    const T     Ci = s / mag_P;
    const T     w = spec.wc;
    const T     wTf = w * spec.Tf;
    const T     gr = w * wTf / (T{1} + (wTf * wTf)); // Re g
    const T     gi = w / (T{1} + (wTf * wTf));       // Im g
    PIDGains<T> res{};
    switch (spec.mode) {
        case PIDMode::P:
            res.Kp = T{1} / mag_P;
            break;
        case PIDMode::PI:
            res.Kp = Cr;
            res.Ki = -Ci * w;
            break;
        case PIDMode::PD:
            res.Kd = Ci / gi;
            res.Kp = Cr - (res.Kd * gr);
            break;
        case PIDMode::PID: {
            const T a = (T{4} * w * gi) - (gr * gr); // = ω²(4 + 3ω²Tf²)/(1 + ω²Tf²)² > 0
            const T bb = (T{2} * Cr * gr) - (T{4} * w * Ci);
            const T sq = damp::sqrt((bb * bb) + (T{4} * a * Cr * Cr));
            // Positive root of a Kd² + bb Kd − C_r² = 0, without cancellation for either sign of bb.
            res.Kd = (bb >= T{0}) ? ((T{2} * Cr * Cr) / (bb + sq)) : ((sq - bb) / (T{2} * a));
            res.Kp = Cr - (res.Kd * gr);
            res.Ki = (res.Kd > T{0}) ? ((res.Kp * res.Kp) / (T{4} * res.Kd)) : T{0};
        } break;
    }
    if (!(res.Kp > T{0}) || res.Ki < T{0} || res.Kd < T{0} || ((spec.mode == PIDMode::PID) && !(res.Kd > T{0}))) {
        return ::damp::detail::refused<PIDGains<T>>(DesignStatus::not_placeable);
    }
    res.Kbc = res.Ki / w; // tracking time T_t = Kbc/Ki = 1/ω_c
    res.Tf = spec.Tf;
    res.b = spec.b;
    res.c = spec.c;
    return res;
}

/**
 * @brief fit_to_response at spec.wc with every gain of sign σ = ±1
 *
 * σ = −1 fits σ·P and negates the gains: the loop C·P is the same, and Kbc keeps Kᵢ's sign.
 */
template<size_t NX, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr PIDGains<T>
fit_at_crossover(const damp::StateSpace<NX, 1, 1, T, NW, NV>& sys, const PIDTuneSpec<T>& spec, T sigma) noexcept {
    const auto Gopt = damp::eval_frf(sys, damp::complex<T>{T{0}, spec.wc});
    if (!Gopt) {
        return ::damp::detail::refused<PIDGains<T>>(DesignStatus::singular);
    }
    auto g = fit_to_response((*Gopt)(0, 0) * sigma, spec);
    if (g.success && (sigma < T{0})) {
        g.Kp = -g.Kp;
        g.Ki = -g.Ki;
        g.Kd = -g.Kd;
        g.Kbc = -g.Kbc;
    }
    return g;
}

/**
 * @brief Routh–Hurwitz: every root of c[0] sⁿ + … + c[n] in the open left half plane
 *
 * Leading exact zeros are skipped. The polynomial is first rescaled in frequency, s = ω₀ s'
 * with ω₀ = |c[n]/c[0]|^{1/n} (cₖ → cₖ/ω₀ᵏ), so its roots have unit geometric mean magnitude
 * and the verdict does not depend on the time unit. A first-column entry must then exceed
 * N·64·ε times the largest entry of the two rows that form it, so a root on (or within
 * roundoff of) the imaginary axis fails, and so does a root at s = 0 (c[n] = 0).
 */
template<size_t N, typename T>
[[nodiscard]] constexpr bool hurwitz(const damp::array<T, N>& coeffs) {
    size_t lead = 0;
    while (lead < N && coeffs[lead] == T{0}) {
        ++lead;
    }
    if (lead >= N) {
        return false;
    }
    const size_t n = N - 1 - lead; // degree
    const T      sgn = (coeffs[lead] > T{0}) ? T{1} : T{-1};
    const T      tol = static_cast<T>(N) * T{64} * eps<T>();
    const T      ratio = damp::abs(coeffs[N - 1] / coeffs[lead]);
    if (!(ratio > T{0}) || !damp::isfinite(ratio)) {
        return false; // a root at s = 0, or a coefficient out of range
    }
    const T           inv_w0 = (n > 0) ? (T{1} / damp::pow(ratio, T{1} / static_cast<T>(n))) : T{1};
    T                 wk = T{1}; // ω₀⁻ⁱ
    damp::array<T, N> r1{};      // sⁿ row: c[0], c[2], …
    damp::array<T, N> r2{};      // sⁿ⁻¹ row: c[1], c[3], …
    for (size_t i = 0; i <= n; ++i) {
        const T ci = sgn * coeffs[lead + i] * wk;
        wk *= inv_w0;
        if (!damp::isfinite(ci)) {
            return false;
        }
        if ((i % 2) == 0) {
            r1[i / 2] = ci;
        } else {
            r2[i / 2] = ci;
        }
    }
    for (size_t row = 0; row < n; ++row) {
        T scale = T{0};
        for (size_t k = 0; k < N; ++k) {
            scale = damp::max(scale, damp::max(damp::abs(r1[k]), damp::abs(r2[k])));
        }
        if (!(r1[0] > tol * scale) || !(r2[0] > tol * scale)) {
            return false;
        }
        damp::array<T, N> next{};
        for (size_t k = 0; (k + 1) < N; ++k) {
            next[k] = r1[k + 1] - ((r1[0] / r2[0]) * r2[k + 1]);
        }
        r1 = r2;
        r2 = next;
    }
    return true;
}

/**
 * @brief Closed-loop certificate: den_P·den_C + num_P·num_C is Hurwitz for C_fb = Kp + Kᵢ/s + Kd s/(1 + Tf s)
 *
 * The integrator and the filter factors enter only when Kᵢ ≠ 0 and Kd·Tf ≠ 0, so the
 * polynomial has no root that a pole-zero cancellation inside C would put on the axis.
 */
template<size_t NP, typename T>
[[nodiscard]] constexpr bool closed_loop_stable(const TransferFunction<NP, NP, T>& plant, const PIDGains<T>& g) {
    const bool        has_i = g.Ki != T{0};
    const bool        has_f = (g.Kd != T{0}) && (g.Tf > T{0});
    const T           Tf = has_f ? g.Tf : T{0};
    damp::array<T, 3> num{(g.Kp * Tf) + g.Kd, g.Kp + (g.Ki * Tf), g.Ki}; // (Kp Tf + Kd) s² + (Kp + Ki Tf) s + Ki
    damp::array<T, 3> den{Tf, T{1}, T{0}};                               // Tf s² + s
    if (!has_i) {                                                        // no integrator: divide both by s
        num = {T{0}, num[0], num[1]};
        den = {T{0}, Tf, T{1}};
    }
    damp::array<T, NP + 2> cl{};
    for (size_t i = 0; i < NP; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            cl[i + j] += (plant.den[i] * den[j]) + (plant.num[i] * num[j]);
        }
    }
    return hurwitz(cl);
}

/**
 * @brief Raise the loop crossover so T −3 dB matches the 1-DOF target
 *
 * Keeps the sample whose bandwidth is closest to @p target. Each step is
 * clamped to [½, 2] so one update cannot jump the ridge where the gains change
 * sign. Every sample is fitted with gains of sign @p sigma, the sign of @p fitted.
 * The search stops at a crossover the law cannot realize or whose closed loop
 * fails the certificate.
 */
template<size_t NX, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr PIDGains<T> raise_loop_for_tracking(
    const damp::StateSpace<NX, 1, 1, T, NW, NV>& sys,
    const TransferFunction<NX + 1, NX + 1, T>&   plant,
    const PIDTuneSpec<T>&                        spec,
    PIDGains<T>                                  fitted,
    T                                            target,
    T                                            sigma
) noexcept {
    T           w = spec.wc;
    PIDGains<T> best = fitted;
    bool        scored = false;
    T           best_rel = T{0};
    const auto  bw0 = tracking_bandwidth(sys, fitted, spec.wc);
    if (bw0 && (target > T{0})) {
        best_rel = damp::abs(*bw0 - target) / target;
        scored = true;
    }
    const T wmin = spec.wc / T{4};
    const T wmax = spec.wc * T{64};
    for (int i = 0; i < 8; ++i) {
        PIDTuneSpec<T> s = spec;
        s.wc = w;
        const auto r = fit_at_crossover(sys, s, sigma);
        if (!r.success || !closed_loop_stable(plant, r)) {
            break;
        }
        const auto bw = tracking_bandwidth(sys, r, spec.wc);
        if (!bw || !(target > T{0})) {
            break;
        }
        const T rel = damp::abs(*bw - target) / target;
        if (!scored || (rel < best_rel)) {
            best = r;
            best_rel = rel;
            scored = true;
        }
        if (rel < static_cast<T>(0.03)) {
            break;
        }
        T step = target / *bw;
        if (step > T{2}) {
            step = T{2};
        }
        if (step < T{0.5}) {
            step = T{0.5};
        }
        T next = w * step;
        if (next < wmin) {
            next = wmin;
        }
        if (next > wmax) {
            next = wmax;
        }
        if (!(damp::abs(next - w) > (w * static_cast<T>(1e-4)))) {
            break;
        }
        w = next;
    }
    return best;
}

} // namespace pidtune_detail

/**
 * @brief Plant-aware PID tune at a crossover
 *
 * Chooses C_fb(jω) so |C_fb P|=1 and arg(C_fb P)= −π + PM, then maps onto the
 * requested law (Ti = 4 Td for PID). The derivative filter Tf is part of the fit, so the
 * delivered C_fb meets both at any Tf. P mode meets the crossover only (Kp = 1/|P(jωc)|);
 * its phase margin is the plant's. For 1-DOF (b = c = 1) that frequency is
 * spec.wc. When b or c is not 1, the loop crossover is raised so the tracking
 * map T_yr = C_ff P / (1 + C_fb P) has the same −3 dB bandwidth as the 1-DOF
 * design at spec.wc. The search keeps the closest match and steps by at most 2x.
 * Back-calculation is seeded to Kbc = Kᵢ/ω, ω the loop crossover of the result (spec.wc
 * for 1-DOF): tracking time T_t = 1/ω.
 *
 * The gains share one sign. The sign of the plant's low-frequency gain num(0) is tried
 * first, so a reverse-acting plant gets Kp, Kᵢ, Kd ≤ 0 (pidtune(−P) is −pidtune(P)), and the
 * other sign is taken only when it alone gives a stable loop.
 *
 * The result is certified: every pole of the closed loop 1 + C_fb P lies in the open left
 * half plane (Routh–Hurwitz on the characteristic polynomial from ss2tf).
 *
 * Status: invalid_parameter if ωc ≤ 0, PM ∉ (0°, 180°), Tf < 0, or an input is non-finite;
 * invalid_sample_time on a discrete plant; singular if P(jωc) is missing or 0; not_placeable if
 * the law cannot realize the crossover with gains of one sign (a PI cannot add phase, a PD
 * cannot remove it); not_stabilizable if the fitted loop is unstable. When neither sign gives
 * a stable loop, the status is that of the plant's sign.
 *
 * @note Compare with MATLAB®'s pidtune(sys, wc) (that call is 60° 1-DOF PID).
 * @see PIDTuneSpec
 */
template<size_t NX, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr PIDGains<T>
pidtune(const damp::StateSpace<NX, 1, 1, T, NW, NV>& sys, const PIDTuneSpec<T>& spec) noexcept {
    if (!pidtune_detail::valid_spec(spec)) {
        return detail::pid_invalid<T>();
    }
    if (sys.Ts != T{0}) {
        return ::damp::detail::refused<PIDGains<T>>(DesignStatus::invalid_sample_time);
    }
    const auto              plant = ::damp::ss2tf(sys);
    const T                 first = (plant.num[NX] < T{0}) ? T{-1} : T{1}; // sign of num(0), den monic
    const damp::array<T, 2> signs{first, -first};
    PIDGains<T>             fitted{};
    T                       sigma = first;
    for (const T s : signs) {
        auto g = pidtune_detail::fit_at_crossover(sys, spec, s);
        if (g.success && !pidtune_detail::closed_loop_stable(plant, g)) {
            g = ::damp::detail::refused<PIDGains<T>>(DesignStatus::not_stabilizable);
        }
        if (g.success || (s == first)) { // on a double failure, report the plant sign's status
            fitted = g;
            sigma = s;
        }
        if (g.success) {
            break;
        }
    }
    if (!fitted.success || !pidtune_detail::weights_move_tracking(fitted)) {
        return fitted;
    }
    PIDTuneSpec<T> spec1 = spec;
    spec1.b = T{1};
    spec1.c = T{1};
    auto target = pidtune_detail::tracking_bandwidth(sys, pidtune_detail::fit_at_crossover(sys, spec1, sigma), spec.wc);
    if (!target) {
        target = spec.wc;
    }
    return pidtune_detail::raise_loop_for_tracking(sys, plant, spec, fitted, *target, sigma);
}

/**
 * @brief Map percent overshoot target to equivalent damping ratio
 *
 * Uses the standard second-order relation:
 *
 *   OS = exp(-zeta*pi/sqrt(1-zeta^2))
 *
 * with OS as a unit fraction (e.g. 10% -> 0.10).
 *
 * @param overshoot_percent Desired percent overshoot in [0, 100)
 * @return Estimated damping ratio zeta
 */
template<typename T = double>
[[nodiscard]] constexpr T damping_ratio_from_overshoot_percent(T overshoot_percent) {
    constexpr T pi = damp::numbers::pi_v<T>;

    if (overshoot_percent <= T{0}) {
        return T{1};
    }

    T os = overshoot_percent / T{100};
    if (os >= T{1}) {
        os = static_cast<T>(0.999);
    }

    const T log_os = damp::log(os);
    return -log_os / damp::sqrt((pi * pi) + (log_os * log_os));
}

/**
 * @brief Approximate phase margin from damping ratio
 *
 * Exact for the standard second-order loop L(s) = ωₙ²/(s(s + 2ζωₙ)):
 *
 *   PM = atan(2*zeta / sqrt(sqrt(1 + 4*zeta^4) - 2*zeta^2))
 *
 * valid for every ζ > 0 (PM → 90° as ζ → ∞). The inner difference is evaluated
 * as 1/(√(1+4ζ⁴) + 2ζ²) to avoid cancellation at large ζ.
 *
 * @param zeta Damping ratio
 * @return Phase margin in degrees; 30° fallback for ζ ≤ 0
 */
template<typename T = double>
[[nodiscard]] constexpr T phase_margin_from_damping_ratio(T zeta) {
    constexpr T pi = damp::numbers::pi_v<T>;

    if (!(zeta > T{0})) {
        return T{30};
    }

    const T z2 = zeta * zeta;
    // √(1+4ζ⁴) − 2ζ² rewritten without the subtraction.
    const T inner = T{1} / (damp::sqrt(T{1} + (T{4} * z2 * z2)) + (T{2} * z2));
    const T pm_rad = damp::atan2(T{2} * zeta, damp::sqrt(inner));
    return pm_rad * T{180} / pi;
}

// ============================================================================
// PID pole placement
// ============================================================================

/**
 * @brief Closed-loop targets for the continuous pid_pole_placement forms
 *
 * The pole pair s² + 2ζω s + ω² and, on a second-order plant, the integrator pole at −ωᵢ.
 */
template<typename T = double>
struct PIDPolePlacementSpec {
    T omega{};       ///< [rad/s] Frequency ω of the placed pole pair (> 0)
    T zeta{T{1}};    ///< Damping ratio ζ of the pair (> 0); 1 = critically damped
    T omega_i{T{0}}; ///< [rad/s] Integrator pole at −ωᵢ, second-order plant only; 0 = PD
    T b{T{1}};       ///< Proportional setpoint weight: 1 = P on error, 0 = P on measurement
};

namespace detail {

/** @brief The spec's scalar domain: ω > 0, ζ > 0, ωᵢ ≥ 0, all finite. */
template<typename T>
[[nodiscard]] constexpr bool valid_placement_spec(const PIDPolePlacementSpec<T>& spec) {
    return (spec.omega > T{0}) && (spec.zeta > T{0}) && (spec.omega_i >= T{0}) && damp::isfinite(spec.omega)
        && damp::isfinite(spec.zeta) && damp::isfinite(spec.omega_i) && damp::isfinite(spec.b);
}

} // namespace detail

/**
 * @brief PI gains that place the closed-loop poles of a first-order plant b₀/(a₁ s + a₀)
 *
 * A PI has two free gains, so it places exactly the two closed-loop poles of a
 * first-order plant. With the plant made monic by its gain, a₁′ = a₁/b₀, a₀′ = a₀/b₀,
 * matching @f$ s^2 + \frac{a_0'+K_p}{a_1'}s + \frac{K_i}{a_1'} @f$ to
 * @f$ s^2 + 2\zeta\omega s + \omega^2 @f$ gives
 * @f[
 *   K_p = 2\zeta\omega a_1' - a_0', \qquad K_i = a_1'\,\omega^2 .
 * @f]
 * This is the kernel every loop in a motion cascade shares: (a₁, a₀) = (L, R) for a dq
 * current axis, (J, b_visc) for a velocity loop. The back-calculation gain is seeded to
 * Kbc = |Kp| (tracking time constant T_t = Kp/Ki = T_i; Kp > 0 here, since a pair that
 * needs Kp ≤ 0 is refused); the setpoint weight b leaves the poles fixed and only moves
 * the reference→output zero (b = 1 PI on error, b = 0 I-P).
 *
 * Status: invalid_parameter if b₀ = 0, a₁/b₀ ≤ 0, ω ≤ 0, ζ ≤ 0 or ωᵢ ≠ 0 (a PI has no third
 * pole); not_placeable if the placement needs Kp ≤ 0 (2ζω a₁′ ≤ a₀′: the plant pole is
 * already faster than the requested pair).
 *
 * @note ω is bounded by the eventual sample rate: keep it roughly an order of magnitude
 *       below the control frequency for the continuous placement to hold.
 *
 * @param plant First-order plant, `{.num = {b₀}, .den = {a₁, a₀}}`
 * @param spec  ω, ζ and the setpoint weight b
 * @return PIDGains with Kp, Ki, Kbc = |Kp| and b
 *
 * @see "Advanced PID Control" (Åström & Hägglund, 2006), §4.4 (setpoint weighting)
 */
template<typename T = double>
[[nodiscard]] constexpr PIDGains<T> pid_pole_placement(const TransferFunction<1, 2, T>& plant, const PIDPolePlacementSpec<T>& spec) {
    const T b0 = plant.num[0];
    if (!detail::valid_static_gain(b0) || !detail::valid_placement_spec(spec) || spec.omega_i != T{0}) {
        return detail::pid_invalid<T>();
    }
    const T a1 = plant.den[0] / b0;
    const T a0 = plant.den[1] / b0;
    if (!detail::positive_finite(a1) || !damp::isfinite(a0)) {
        return detail::pid_invalid<T>();
    }
    const T Kp = (T{2} * spec.zeta * spec.omega * a1) - a0;
    if (!(Kp > T{0})) {
        return ::damp::detail::refused<PIDGains<T>>(DesignStatus::not_placeable);
    }
    // Parallel form: T_t = Kbc/Ki. Kbc = |Kp| ⇒ T_t = T_i (Åström default).
    const T Ki = a1 * spec.omega * spec.omega;
    return PIDGains<T>{
        .Kp = Kp,
        .Ki = Ki,
        .Kbc = detail::placed_kbc(Kp, Ki),
        .b = spec.b,
    };
}

/**
 * @brief PD or PID gains that place the closed-loop poles of a second-order plant b₀/(a₂ s² + a₁ s + a₀)
 *
 * The derivative is filtered, Kd s/(1 + Tf s), and the filter pole is placed with the others:
 * with the plant made monic by its gain (a₂′, a₁′, a₀′ = a₂/b₀, a₁/b₀, a₀/b₀), a PID closes to
 * a₂′Tf s⁴ + (a₂′ + a₁′Tf) s³ + (a₁′ + a₀′Tf + Kp Tf + Kd) s² + (a₀′ + Kp + Kᵢ Tf) s + Kᵢ.
 * Matching a₂′Tf (s² + 2ζω s + ω²)(s + ωᵢ)(s + ω_f) = a₂′Tf (s⁴ + d₃ s³ + d₂ s² + d₁ s + d₀),
 * with the filter pole at ω_f = 10ω, gives
 * @f[
 *   \frac{1}{T_f} = d_3 - \frac{a_1'}{a_2'},\quad
 *   K_i = a_2' T_f d_0,\quad
 *   K_p = a_2' T_f d_1 - a_0' - K_i T_f,\quad
 *   K_d = a_2' T_f d_2 - a_1' - (a_0' + K_p) T_f ,
 * @f]
 * and ωᵢ = 0 is the PD that places the pair and the filter pole alone. A rigid inertia
 * 1/(J s²) is `{.num = {1}, .den = {J, 0, 0}}`; a common load-holding choice is ωᵢ = ω/8.
 *
 * Derivative is on the measurement (c = 0, no D kick on a reference step). Back-calculation
 * is seeded to Kbc = |Kp| (T_t = T_i).
 *
 * Status: invalid_parameter if b₀ = 0, a₂/b₀ ≤ 0, ω ≤ 0, ζ ≤ 0, ωᵢ < 0 or a coefficient is
 * non-finite; not_placeable if the placement needs Kp ≤ 0, Kd < 0 or Tf ≤ 0 (the plant is
 * already faster or more damped).
 *
 * @note Compare with MATLAB®'s pidtune on tf(1, [J 0 0]) for the PD part.
 *
 * @param plant Second-order plant, `{.num = {b₀}, .den = {a₂, a₁, a₀}}`
 * @param spec  ω, ζ, ωᵢ and the setpoint weight b
 * @return PIDGains with Kp, Ki, Kd, Kbc, b, c = 0 and Tf (0 when Kd = 0)
 *
 * @see "Feedback Control of Dynamic Systems" (Franklin, Powell, Emami-Naeini), PD on 1/s²
 */
template<typename T = double>
[[nodiscard]] constexpr PIDGains<T> pid_pole_placement(const TransferFunction<1, 3, T>& plant, const PIDPolePlacementSpec<T>& spec) {
    const T b0 = plant.num[0];
    if (!detail::valid_static_gain(b0) || !detail::valid_placement_spec(spec)) {
        return detail::pid_invalid<T>();
    }
    const T a2 = plant.den[0] / b0;
    const T a1 = plant.den[1] / b0;
    const T a0 = plant.den[2] / b0;
    if (!detail::positive_finite(a2) || !damp::isfinite(a1) || !damp::isfinite(a0)) {
        return detail::pid_invalid<T>();
    }
    const T w = spec.omega;
    const T wi = spec.omega_i;
    const T wf = T{10} * w;
    const T two_zw = T{2} * spec.zeta * w;
    // (s² + 2ζω s + ω²)(s + ωᵢ) = s³ + c₁ s² + c₂ s + c₃, then × (s + ω_f).
    const T c1 = two_zw + wi;
    const T c2 = (w * w) + (two_zw * wi);
    const T c3 = w * w * wi;
    const T inv_Tf = (c1 + wf) - (a1 / a2);
    if (!(inv_Tf > T{0})) {
        return ::damp::detail::refused<PIDGains<T>>(DesignStatus::not_placeable);
    }
    const T Tf = T{1} / inv_Tf;
    const T Ki = a2 * Tf * (c3 * wf);
    const T Kp = (a2 * Tf * (c3 + (c2 * wf))) - a0 - (Ki * Tf);
    const T Kd = (a2 * Tf * (c2 + (c1 * wf))) - a1 - ((a0 + Kp) * Tf);
    if (!(Kp > T{0}) || !(Kd >= T{0})) {
        return ::damp::detail::refused<PIDGains<T>>(DesignStatus::not_placeable);
    }
    PIDGains<T> r{};
    r.Kp = Kp;
    r.Ki = Ki;
    r.Kd = Kd;
    r.Kbc = detail::placed_kbc(Kp, Ki);
    r.b = spec.b;
    r.c = T{0};                     // PI-D / PD on measurement
    r.Tf = (Kd > T{0}) ? Tf : T{0}; // Kd = 0: the filter factor cancels from the loop
    return r;
}

namespace detail {

/** @brief 1 − e^{−x} over x, without the cancellation of 1 − e^{−x} at small |x|. */
template<typename T>
[[nodiscard]] constexpr T one_minus_exp_neg_over_x(T x) {
    // Below |x| = 1/16 the series through x⁷/8! is accurate to about 3ε in double, and the
    // direct form loses at most ε/|x| = 16ε above it.
    if (damp::abs(x) < static_cast<T>(0.0625)) {
        // Σ (−x)ᵏ/(k + 1)!, k = 0 … 7.
        T term = T{1};
        T sum = T{1};
        for (int k = 1; k <= 7; ++k) {
            term *= -x / static_cast<T>(k + 1);
            sum += term;
        }
        return sum;
    }
    return (T{1} - damp::exp(-x)) / x;
}

} // namespace detail

/**
 * @brief Discrete PI or PID gains that place the closed-loop poles of a first-order plant sampled at Ts
 *
 * The poles are s-plane poles, mapped to z = e^{s·Ts} as place(sys, poles, Ts) maps them. The
 * plant b₀/(a₁ s + a₀) is ZOH-discretized to G(z) = b/(z − a), a = e^{−x}, x = a₀Ts/a₁,
 * b = (b₀Ts/a₁)(1 − e^{−x})/x (evaluated without cancellation at small x, and b = b₀Ts/a₁ for
 * an integrator). Two poles give a PI (the integrator at z = 1 makes the closed loop second
 * order); three give a PID whose backward-difference derivative adds the third pole. The gains
 * are in the continuous form that PIDGains::discretize(Ts) maps back exactly, so discretize
 * with the same Ts.
 *
 * A pole product above the sampled plant pole (z₁z₂ > a for the PI) needs Kp < 0, and
 * z₁z₂ = a gives Kp = 0, an I-only law. The PID's derivative gain is Kd = −z₁z₂z₃·Ts/b, which
 * is negative for every s-plane set on a plant with b > 0. Back-calculation is seeded to
 * Kbc = sgn(Kᵢ)·max(|Kp|, |Kᵢ|·Ts): tracking time T_t = Kbc/Kᵢ = max(|T_i|, Ts), positive
 * whatever the sign of Kp, and never shorter than one sample (the tick's correction
 * Kᵢ·Ts/Kbc stays at most 1, so it cannot overshoot the saturation).
 *
 * Status: invalid_sample_time if Ts ≤ 0; invalid_parameter if b₀ = 0, a₁ = 0, a coefficient or
 * pole is non-finite, a pole has Re(s) ≥ 0 (as a z-plane value in [0, 1] read as s does), or a
 * complex pole has no conjugate; singular if the sampled gain b is zero.
 *
 * @tparam NP Number of placed poles: 2 (PI) or 3 (PID)
 * @param plant  First-order plant, `{.num = {b₀}, .den = {a₁, a₀}}`
 * @param poles  Desired closed-loop poles (s-plane, Re(s) < 0), real or in conjugate pairs
 * @param Ts     [s] Sample period
 * @return PIDGains (Tf = 0) with the Kbc seed above
 */
template<size_t NP, typename T = double>
    requires(NP == 2 || NP == 3)
[[nodiscard]] constexpr PIDGains<T> pid_pole_placement(
    const TransferFunction<1, 2, T>&         plant,
    const damp::array<damp::complex<T>, NP>& poles,
    T                                        Ts
) {
    using Cplx = damp::complex<T>;
    if (!detail::positive_finite(Ts)) {
        return ::damp::detail::refused<PIDGains<T>>(DesignStatus::invalid_sample_time);
    }
    const T b0 = plant.num[0];
    const T a1 = plant.den[0];
    const T a0 = plant.den[1];
    if (!detail::valid_static_gain(b0) || !detail::valid_static_gain(a1) || !damp::isfinite(a0)) {
        return detail::pid_invalid<T>();
    }
    // z = e^{s·Ts}; the elementary symmetric functions σ = Σz, ρ = Σzᵢzⱼ, π = Πz are real for a
    // conjugate-closed set.
    damp::array<Cplx, NP> z{};
    T                     zmag = T{0};
    for (size_t i = 0; i < NP; ++i) {
        // Re(s) < 0: the loop is stable, and a z-plane value in [0, 1] read as s is refused.
        if (!damp::isfinite(poles[i].real()) || !damp::isfinite(poles[i].imag()) || !(poles[i].real() < T{0})) {
            return detail::pid_invalid<T>();
        }
        const T er = damp::exp(poles[i].real() * Ts);
        const auto [sn, cs] = damp::sincos(poles[i].imag() * Ts);
        z[i] = Cplx{er * cs, er * sn};
        zmag += er;
    }
    Cplx sigma{};
    Cplx rho{};
    Cplx prod{T{1}};
    for (size_t i = 0; i < NP; ++i) {
        sigma += z[i];
        prod *= z[i];
        for (size_t j = i + 1; j < NP; ++j) {
            rho += z[i] * z[j];
        }
    }
    const T scale = damp::max(T{1}, zmag);
    const T tol = T{64} * static_cast<T>(NP) * eps<T>() * scale * scale * scale;
    if (damp::abs(sigma.imag()) > tol || damp::abs(rho.imag()) > tol || damp::abs(prod.imag()) > tol) {
        return detail::pid_invalid<T>(); // a complex pole without its conjugate
    }
    const T x = a0 * Ts / a1;
    const T a = damp::exp(-x);
    const T b = (b0 * Ts / a1) * detail::one_minus_exp_neg_over_x(x);
    if (b == T{0}) {
        return ::damp::detail::refused<PIDGains<T>>(DesignStatus::singular);
    }
    if constexpr (NP == 2) {
        // Deploy is PIDGains::discretize + PIDController: I ← I + (Ki Ts) e[k],
        // so C(z) = Kp + Ki Ts · z/(z−1) (backward Euler). Match
        //   (z−a)(z−1) + b (Kp (z−1) + Ki Ts z) = z² − σ z + π:
        //   a − b Kp = π                 => Kp = (a − π) / b
        //   b Kp + b Ki Ts = 1+a − σ     => Ki from the leftover.
        const T Kp = (a - prod.real()) / b;
        const T Ki = ((T{1} + a) - sigma.real() - (b * Kp)) / (b * Ts);
        return PIDGains<T>{.Kp = Kp, .Ki = Ki, .Kbc = detail::placed_kbc(Kp, Ki, Ts)};
    } else {
        // Deploy (Tf = 0): C(z) = Kp + Ki Ts z/(z−1) + (Kd/Ts)(z−1)/z. Closed loop:
        //   z³ + (−(1+a) + b Kp + b Ki Ts + b Kd/Ts) z² + (a − b Kp − 2 b Kd/Ts) z + b Kd/Ts
        // against z³ − σ z² + ρ z − π.
        const T Kd = -(prod.real() * Ts) / b;                                           // b Kd/Ts = −π
        const T Kp = (a - rho.real() - (T{2} * b * Kd / Ts)) / b;                       // b Kp + 2 b Kd/Ts = a − ρ
        const T Ki = ((T{1} + a - sigma.real()) - (b * Kp) - (b * Kd / Ts)) / (b * Ts); // b Kp + b Ki Ts + b Kd/Ts = 1 + a − σ
        return PIDGains<T>{.Kp = Kp, .Ki = Ki, .Kd = Kd, .Kbc = detail::placed_kbc(Kp, Ki, Ts)};
    }
}

namespace detail {

/** @brief A SISO plant as the all-pole transfer function b₀/den(s), or why it is not one. */
template<size_t NX, typename T>
struct AllPoleTF {
    TransferFunction<1, NX + 1, T> tf{};                                    ///< b₀/den(s), den monic
    DesignStatus                   status{DesignStatus::invalid_parameter}; ///< ok, or why the plant is refused
};

/**
 * @brief The all-pole transfer function b₀/den(s) of a SISO state-space plant
 *
 * ss2tf gives the monic denominator. A numerator coefficient above the constant one is a finite
 * zero unless it is cancellation residue, below √ε times the size of the terms that formed it,
 * so the verdict does not depend on the units of the states. Noise channels (G, H) do not enter.
 * invalid_sample_time on a discrete plant; invalid_parameter if D ≠ 0, the plant has a finite
 * zero, or there is no gain from u to y.
 */
template<size_t NX, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr AllPoleTF<NX, T> all_pole_tf(const StateSpace<NX, 1, 1, T, NW, NV>& sys) {
    const StateSpace<NX, 1, 1, T> plant{.A = sys.A, .B = sys.B, .C = sys.C, .D = sys.D, .Ts = sys.Ts};
    if (plant.Ts != T{0}) {
        return {.status = DesignStatus::invalid_sample_time};
    }
    if (plant.D(0, 0) != T{0}) {
        return {.status = DesignStatus::invalid_parameter};
    }
    const auto tf = ::damp::ss2tf(plant);
    const auto scale = ::damp::zpk_detail::leverrier_num_scale(plant);
    const T    tol = damp::sqrt(eps<T>());
    for (size_t j = 1; j < NX; ++j) {
        if (damp::abs(tf.num[j]) > (tol * scale[j])) {
            return {.status = DesignStatus::invalid_parameter}; // a finite zero
        }
    }
    if (!(damp::abs(tf.num[NX]) > (tol * scale[NX]))) {
        return {.status = DesignStatus::invalid_parameter}; // no path from u to y
    }
    return {.tf = {.num = {tf.num[NX]}, .den = tf.den}, .status = DesignStatus::ok};
}

} // namespace detail

/**
 * @brief PI (first-order) or PD / PID (second-order) gains that place the closed-loop poles of a state-space plant
 *
 * Converts the plant with ss2tf and forwards to the transfer-function form of the same order,
 * so one model feeds pid_pole_placement, pidtune and place. The plant must be continuous and
 * all-pole: D = 0 and no finite zero. Noise channels (G, H) are ignored.
 *
 * Status: invalid_sample_time on a discrete plant; invalid_parameter if D ≠ 0, the plant has a
 * finite zero, or it has no gain from u to y; otherwise as the transfer-function form.
 *
 * @param plant Continuous SISO plant with one or two states
 * @param spec  ω, ζ, ωᵢ and the setpoint weight b
 * @return PIDGains as the transfer-function form returns them
 */
template<size_t NX, typename T = double, size_t NW = 0, size_t NV = 0>
    requires(NX == 1 || NX == 2)
[[nodiscard]] constexpr PIDGains<T> pid_pole_placement(const StateSpace<NX, 1, 1, T, NW, NV>& plant, const PIDPolePlacementSpec<T>& spec) {
    const auto all_pole = detail::all_pole_tf(plant);
    if (all_pole.status != DesignStatus::ok) {
        return ::damp::detail::refused<PIDGains<T>>(all_pole.status);
    }
    return pid_pole_placement(all_pole.tf, spec);
}

/**
 * @brief Discrete PI or PID gains that place the s-plane poles of a first-order state-space plant sampled at Ts
 *
 * Converts the continuous plant with ss2tf and forwards to the transfer-function form. Status:
 * invalid_sample_time on a discrete plant or Ts ≤ 0; invalid_parameter if D ≠ 0 or there is no
 * gain from u to y; otherwise as the transfer-function form.
 *
 * @tparam NP Number of placed poles: 2 (PI) or 3 (PID)
 * @param plant Continuous SISO plant with one state
 * @param poles Desired closed-loop poles (s-plane, Re(s) < 0), real or in conjugate pairs
 * @param Ts    [s] Sample period
 * @return PIDGains (Tf = 0)
 */
template<size_t NP, typename T = double, size_t NW = 0, size_t NV = 0>
    requires(NP == 2 || NP == 3)
[[nodiscard]] constexpr PIDGains<T> pid_pole_placement(
    const StateSpace<1, 1, 1, T, NW, NV>&    plant,
    const damp::array<damp::complex<T>, NP>& poles,
    T                                        Ts
) {
    const auto all_pole = detail::all_pole_tf(plant);
    if (all_pole.status != DesignStatus::ok) {
        return ::damp::detail::refused<PIDGains<T>>(all_pole.status);
    }
    return pid_pole_placement(all_pole.tf, poles, Ts);
}

/**
 * @brief PI gains for an R–L current loop 1/(L s + R) at any bandwidth ω
 *
 * At or above the plant pole R/L, the pid_pole_placement pair: Kp = 2ζωL − R, Ki = Lω².
 * Below it, or where the pair would need Kp ≤ 0 (ζ ≤ ½), the PI zero cancels the plant
 * pole and the loop closes on one pole at −ω:
 * @f[
 *   K_p = \omega L, \qquad K_i = \omega R ,
 * @f]
 * the pair's gains at ω = R/L, ζ = 1. With b = 1 the reference response is then ω/(s + ω);
 * ζ has no effect. At ζ ≥ 1 both rules are LQ-optimal for some weights (|1 + L(jν)| ≥ 1 at
 * every ν).
 *
 * Status: as pid_pole_placement, without not_placeable.
 *
 * @param L    [H] inductance (> 0)
 * @param R    [Ω] resistance
 * @param spec ω, ζ and the setpoint weight b
 * @return PIDGains with Kp, Ki, Kbc = |Kp| and b
 *
 * @see pid_pole_placement
 * @see "When Is a Linear Control System Optimal?" (Kalman, 1964), the return-difference condition
 */
template<typename T = double>
[[nodiscard]] constexpr PIDGains<T> current_loop_pi(T L, T R, const PIDPolePlacementSpec<T>& spec) {
    const auto pair = pid_pole_placement(TransferFunction<1, 2, T>{.num = {T{1}}, .den = {L, R}}, spec);
    const bool slow = (pair.status == DesignStatus::not_placeable) || (pair.success && (spec.omega * L < R));
    if (!slow) {
        return pair;
    }
    // The R/L pole is faster than ω: the PI zero cancels it, one closed-loop pole at −ω.
    const T Kp = spec.omega * L;
    const T Ki = spec.omega * R;
    return PIDGains<T>{
        .Kp = Kp,
        .Ki = Ki,
        .Kbc = detail::placed_kbc(Kp, Ki),
        .b = spec.b,
    };
}

} // namespace design
} // namespace damp
