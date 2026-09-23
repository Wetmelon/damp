// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file pid_design.hpp
 * @brief PID design rules
 *
 * @defgroup pid_design Modelless PID Design
 * @brief Classical PID tuning rules and spec-driven synthesizers (ZN, SIMC, bandwidth, …)
 *
 * These methods use experimentally-determined parameters (ultimate gain/period,
 * step response characteristics, or desired closed-loop bandwidth) to compute
 * PID gains. All functions return continuous PIDResult; deploy via
 * `.discretize(Ts)` → PIDController, or ContinuousPID for variable-rate.
 *
 * @note PIDResult carries no success flag. Non-physical inputs (K = 0, L ≤ 0,
 *       tau ≤ 0, Tu ≤ 0, Ts ≤ 0, etc.) return all-zero gains so callers never
 *       see Inf/NaN from a closed-form divide-by-zero.
 *
 * Includes:
 * - Ziegler-Nichols (ultimate gain method and step response method)
 * - Cohen-Coon (first-order plus dead-time models)
 * - Lambda tuning (setpoint tracking, robustness-focused)
 * - Bandwidth-based design (specify desired closed-loop bandwidth)
 * - Direct pole placement (specify desired closed-loop poles for discretized PID)
 * - Double-integrator / inertia PID (bandwidth + optional I pole)
 * - SIMC (Skogestad IMC) rules for FOPDT/SOPDT models
 * - Tyreus-Luyben (conservative Ziegler-Nichols variant)
 */

#include <limits>

#include "damp/backend.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/math/math.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

namespace damp {

namespace design {

/**
 * @brief PID controller type selection for tuning methods
 */
enum class PIDType {
    P,
    PI,
    PD,
    PID,
};

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
 * @param Ts   Sample period [s]; unused by the continuous formulas (API parity)
 * @param type Controller type (P, PI, PD, or PID)
 * @return PIDResult with tuned gains; zeros if Tu ≤ 0 (PI/PD/PID) or inputs non-physical
 */
template<typename T = double>
[[nodiscard]] constexpr PIDResult<T>
ziegler_nichols(T Ku, T Tu, [[maybe_unused]] T Ts, PIDType type = PIDType::PID) {
    // PI/PD/PID divide by Tu; P-only uses Ku alone.
    if (type != PIDType::P && !(Tu > T{0})) {
        return PIDResult<T>{};
    }

    T Kp{};
    T Ki{};
    T Kd{};

    switch (type) {
        case PIDType::P:
            Kp = static_cast<T>(0.5) * Ku;
            break;
        case PIDType::PI:
            Kp = static_cast<T>(0.45) * Ku;
            Ki = Kp / (Tu / static_cast<T>(1.2));
            break;
        case PIDType::PD:
            Kp = static_cast<T>(0.8) * Ku;
            Kd = Kp * Tu / T{8};
            break;
        case PIDType::PID:
            Kp = static_cast<T>(0.6) * Ku;
            Ki = Kp / (Tu / T{2});
            Kd = Kp * Tu / T{8};
            break;
    }
    return PIDResult<T>{Kp, Ki, Kd};
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
 * @param Ts   Sample period [s]; unused by the continuous formulas (API parity)
 * @param type Controller type (P, PI, or PID)
 * @return PIDResult with tuned gains; zeros if K = 0 or L ≤ 0
 */
template<typename T = double>
[[nodiscard]] constexpr PIDResult<T>
ziegler_nichols_step(T K, T L, T tau, [[maybe_unused]] T Ts, PIDType type = PIDType::PID) {
    // a = tau/(K·L); Ki/Kd formulas also divide by L.
    if (K == T{0} || !(L > T{0})) {
        return PIDResult<T>{};
    }

    T Kp{};
    T Ki{};
    T Kd{};

    T a = tau / (K * L); // normalized gain
    switch (type) {
        case PIDType::P:
            Kp = a;
            break;
        case PIDType::PI:
            Kp = static_cast<T>(0.9) * a;
            Ki = Kp / (L * static_cast<T>(3.33));
            break;
        case PIDType::PD:
            Kp = a;
            Kd = Kp * L * static_cast<T>(0.5);
            break;
        case PIDType::PID:
            Kp = static_cast<T>(1.2) * a;
            Ki = Kp / (T{2} * L);
            Kd = Kp * static_cast<T>(0.5) * L;
            break;
    }
    return PIDResult<T>{Kp, Ki, Kd};
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
 * @param Ts   Sample period [s]; forwarded to Ziegler–Nichols on P/PD fallback
 * @param type Controller type (PI or PID only; others fall back to ZN)
 * @return PIDResult with tuned gains; zeros if Tu ≤ 0 for PI/PID
 */
template<typename T = double>
[[nodiscard]] constexpr PIDResult<T>
tyreus_luyben(T Ku, T Tu, T Ts, PIDType type = PIDType::PID) {
    T Kp{};
    T Ki{};
    T Kd{};

    switch (type) {
        case PIDType::PI:
            if (!(Tu > T{0})) {
                return PIDResult<T>{};
            }
            Kp = Ku / static_cast<T>(3.2);
            Ki = Kp / (static_cast<T>(2.2) * Tu);
            break;
        case PIDType::PID:
            if (!(Tu > T{0})) {
                return PIDResult<T>{};
            }
            Kp = Ku / static_cast<T>(2.2);
            Ki = Kp / (static_cast<T>(2.2) * Tu);
            Kd = Kp * Tu / static_cast<T>(6.3);
            break;
        default:
            // Fall back to Ziegler-Nichols for P and PD
            return ziegler_nichols(Ku, Tu, Ts, type);
    }
    return PIDResult<T>{Kp, Ki, Kd};
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
 * AMIGO PI (Advanced PID Control / JPC 2004 family):
 *
 *     Kp = (0.25 + 0.45 T/L) / Kₛ
 *     Ti = (0.4 L + 0.8 T) / (L + 0.1 T) · L
 *     Ki = Kp / Ti
 *
 * @param Ks Static process gain (from biased relay or step)
 * @param Ku Ultimate gain
 * @param Tu Ultimate period [s]
 * @param Ts Sample period [s]; unused by the continuous formulas (API parity)
 * @return PIDResult PI gains (Kd = 0); zeros if inputs are non-physical
 *
 * @see "Revisiting the Ziegler-Nichols step response method" (Åström & Hägglund,
 *      J. Process Control 14, 2004)
 * @see "Advanced PID Control" (Åström & Hägglund, 2006), ch. 6 & 8
 */
template<typename T = double>
[[nodiscard]] constexpr PIDResult<T>
amigo_kappa_tau(T Ks, T Ku, T Tu, [[maybe_unused]] T Ts = T{0}) {
    if (!(Ks > T{0}) || !(Ku > T{0}) || !(Tu > T{0})) {
        return PIDResult<T>{};
    }
    const T wu = (T{2} * damp::numbers::pi_v<T>) / Tu;
    const T prod = Ks * Ku;
    if (prod <= T{1}) {
        // No real FOPDT consistent with this ultimate point.
        return PIDResult<T>{};
    }
    const T tau = damp::sqrt((prod * prod) - T{1}) / wu;
    const T L = (damp::numbers::pi_v<T> - damp::atan(wu * tau)) / wu;
    if (!(L > T{0}) || !(tau > T{0})) {
        return PIDResult<T>{};
    }
    const T Kp = (static_cast<T>(0.25) + (static_cast<T>(0.45) * tau / L)) / Ks;
    const T Ti = ((static_cast<T>(0.4) * L) + (static_cast<T>(0.8) * tau)) / (L + (static_cast<T>(0.1) * tau)) * L;
    if (!(Ti > T{0}) || !(Kp > T{0})) {
        return PIDResult<T>{};
    }
    const T Ki = Kp / Ti;
    return PIDResult<T>{Kp, Ki, T{0}};
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
 * @param Ts   Sample period [s]; unused by the continuous formulas (API parity)
 * @param type Controller type (P, PI, PD, or PID)
 * @return PIDResult with tuned gains; zeros if K = 0, L ≤ 0, or tau ≤ 0
 */
template<typename T = double>
[[nodiscard]] constexpr PIDResult<T>
cohen_coon(T K, T L, T tau, [[maybe_unused]] T Ts, PIDType type = PIDType::PID) {
    // r = L/tau and a = tau/(K·L) both require physical FOPDT parameters.
    if (K == T{0} || !(L > T{0}) || !(tau > T{0})) {
        return PIDResult<T>{};
    }

    T Kp{};
    T Ki{};
    T Kd{};

    T r = L / tau; // dead-time ratio
    T a = tau / (K * L);

    switch (type) {
        case PIDType::P:
            Kp = a * (T{1} + (r / T{3}));
            break;
        case PIDType::PI: {
            Kp = a * (static_cast<T>(0.9) + (r / T{12}));
            T Ti = L * (T{30} + (T{3} * r)) / (T{9} + (T{20} * r));
            Ki = Kp / Ti;
        } break;
        case PIDType::PD: {
            Kp = a * (static_cast<T>(1.24) + (r / T{5}));
            T Td = L * T{6} / (T{22} + (T{3} * r));
            Kd = Kp * Td;
        } break;
        case PIDType::PID: {
            Kp = a * ((T{4} / T{3}) + (r / T{4}));
            T Ti = L * (T{32} + (T{6} * r)) / (T{13} + (T{8} * r));
            T Td = L * T{4} / (T{11} + (T{2} * r));
            Ki = Kp / Ti;
            Kd = Kp * Td;
        } break;
    }
    return PIDResult<T>{Kp, Ki, Kd};
}

// ============================================================================
// SIMC (Skogestad IMC) Tuning Rules
// ============================================================================

/**
 * @brief SIMC (Skogestad Internal Model Control) tuning for FOPDT models
 *
 * Simple, robust tuning rules from Skogestad (2003). The user specifies
 * a single tuning parameter tau_c (desired closed-loop time constant).
 * Rule of thumb: tau_c = max(tau, 8*L) for good robustness.
 *
 * @param K     Static gain
 * @param L     Dead time (seconds)
 * @param tau   Time constant (seconds)
 * @param tau_c Desired closed-loop time constant (seconds). Larger = more robust.
 * @param type  Controller type (PI or PID)
 * @return PIDResult with tuned gains; zeros if K = 0 or (tau_c + L) ≤ 0
 */
template<typename T = double>
[[nodiscard]] constexpr PIDResult<T>
simc(T K, T L, T tau, T tau_c, PIDType type = PIDType::PI) {
    // Kp = tau / (K·(tau_c+L)); PI/PID also need Ti > 0 for Ki = Kp/Ti.
    const T den = tau_c + L;
    if (K == T{0} || !(den > T{0})) {
        return PIDResult<T>{};
    }

    T Kp{};
    T Ki{};
    T Kd{};

    switch (type) {
        case PIDType::P:
            Kp = tau / (K * den);
            break;
        case PIDType::PI: {
            Kp = tau / (K * den);
            T Ti = (tau < T{4} * den) ? tau : T{4} * den;
            if (!(Ti > T{0})) {
                return PIDResult<T>{};
            }
            Ki = Kp / Ti;
        } break;
        case PIDType::PD: {
            Kp = tau / (K * den);
            Kd = Kp * L / T{2};
        } break;
        case PIDType::PID: {
            Kp = tau / (K * den);
            T Ti = (tau < T{4} * den) ? tau : T{4} * den;
            if (!(Ti > T{0})) {
                return PIDResult<T>{};
            }
            Ki = Kp / Ti;
            Kd = Kp * L / T{2};
        } break;
    }
    return PIDResult<T>{Kp, Ki, Kd};
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
 * @param lambda Desired closed-loop time constant (seconds). Must be > L.
 * @param Ts     Sampling time (seconds)
 * @return PIDResult (PI controller gains); zeros if K = 0, tau ≤ 0, or (lambda + L) ≤ 0
 */
template<typename T = double>
[[nodiscard]] constexpr PIDResult<T>
lambda_tuning(T K, T L, T tau, T lambda, [[maybe_unused]] T Ts) {
    const T den = lambda + L;
    if (K == T{0} || !(tau > T{0}) || !(den > T{0})) {
        return PIDResult<T>{};
    }
    T Kp = tau / (K * den);
    T Ki = Kp / tau;
    return PIDResult<T>{Kp, Ki, T{0}};
}

// ============================================================================
// Bandwidth-Based PID Design (model-free)
// ============================================================================

/**
 * @brief Design PID from desired bandwidth and phase margin
 *
 * Given a desired closed-loop bandwidth ωbw and phase margin φm,
 * computes PID gains that achieve approximately those specifications.
 * This is a model-free method — the user specifies desired performance
 * without needing a plant model.
 *
 * The controller is C(s) = Kp + Ki/s + Kd*s with:
 * - Crossover at ωbw: |C(jωbw)| = 1 (assuming unit-gain plant at crossover)
 * - Phase margin φm at crossover
 *
 * @param wbw          Desired bandwidth (crossover frequency, rad/s)
 * @param phase_margin Desired phase margin (degrees, default: 60°)
 * @param Ts           Sampling time (seconds)
 * @param type         Controller type (PI or PID)
 * @return PIDResult with tuned gains
 */
template<typename T = double>
[[nodiscard]] constexpr PIDResult<T>
pid_from_bandwidth(T wbw, T phase_margin, [[maybe_unused]] T Ts, PIDType type = PIDType::PID) {
    constexpr T pi = damp::numbers::pi_v<T>;
    T           phi = phase_margin * pi / T{180}; // Convert to radians
    T           desired_phase = phi - pi;         // Phase of C(jω) at crossover

    T Kp{};
    T Ki{};
    T Kd{};

    switch (type) {
        case PIDType::P:
            Kp = T{1}; // Unit gain at crossover
            break;
        case PIDType::PI: {
            // PI controller: C(jω) = Kp - j·Ki/ω
            // Phase of PI at crossover = -atan(Ki/(Kp·ω))
            // We want total loop phase margin = phi_m
            // For unit-gain plant: PM = 90° + atan(Kp·ω/Ki) [since PI phase is in [-90°,0°]]
            // So atan(Kp·wbw/Ki) = phi_m - 90° = phi - π/2
            // With |C(jωbw)| = 1: Kp² + (Ki/wbw)² = 1
            T alpha = phi - (pi / T{2}); // angle from the geometric relationship
            // If phi < π/2 (phase margin < 90°), alpha < 0, use |alpha|
            // tan(alpha) = Kp·wbw / Ki, and Kp² + Ki²/wbw² = 1
            // Let r = Ki/(Kp·wbw), then Kp²(1 + r²) = 1
            // r = 1/tan(alpha) = cos(alpha)/sin(alpha)
            if (alpha > T{0}) {
                const auto [s, c] = damp::sincos(alpha);
                T r = c / s; // cot(alpha)
                Kp = T{1} / damp::sqrt(T{1} + (r * r));
                Ki = Kp * r * wbw;
            } else {
                // Phase margin > 90°: PI alone overshoots; use simple heuristic
                const auto [s, c] = damp::sincos(phi);
                Kp = c;
                Ki = s * wbw;
            }
        } break;
        case PIDType::PD: {
            // C(jω) = Kp + Kd*jω
            // Phase = atan2(Kd*ω, Kp), Magnitude = sqrt(Kp² + (Kd*ω)²) = 1
            T tan_phi_d = damp::tan(desired_phase + pi);
            Kp = damp::cos(desired_phase + pi);
            Kd = Kp * tan_phi_d / wbw;
            if (Kp < T{0}) {
                Kp = -Kp;
                Kd = -Kd;
            }
        } break;
        case PIDType::PID: {
            // For PID, place the zero pair symmetrically around crossover
            // Ti = 1/(ωbw * tan(φ/2)), Td = Ti/4
            T half_phi = (phi - (pi / T{2})) / T{2};
            if (half_phi <= T{0}) {
                half_phi = pi / T{12}; // minimum 15° half-angle
            }
            T Ti = T{1} / (wbw * damp::tan(half_phi));
            T Td = Ti / T{4};
            Kp = T{1}; // Normalized for unit plant gain at crossover
            Ki = Kp / Ti;
            Kd = Kp * Td;
        } break;
    }
    return PIDResult<T>{Kp, Ki, Kd};
}

/**
 * @brief Knobs for plant-aware @ref pidtune
 *
 * Crossover, phase margin, law, derivative filter Tf, and 2-DOF weights b / c.
 * Defaults match a 1-DOF PID at 60°.
 * When b or c is not 1, the loop crossover is raised so T (r→y) −3 dB matches
 * the 1-DOF design at wc.
 */
template<typename T = double>
struct PidTuneSpec {
    T       wc{};                     ///< 1-DOF loop crossover ωc [rad/s]; 2-DOF matches T r→y −3 dB to this 1-DOF design
    T       phase_margin_deg = T{60}; ///< Phase margin at the loop crossover [deg]
    PIDType type = PIDType::PID;      ///< P / PI / PD / PID
    T       Tf = T{0};                ///< Derivative filter Tf [s]; 0 = off
    T       b = T{1};                 ///< P setpoint weight (0 = I-PD, 1 = on error)
    T       c = T{1};                 ///< D setpoint weight (0 = PI-D, 1 = on error)
};

namespace pidtune_detail {

template<typename T>
[[nodiscard]] constexpr bool weights_move_tracking(const PIDResult<T>& r) {
    const T    eps = static_cast<T>(1e-6);
    const bool p_w = (damp::abs(r.b - T{1}) > eps) && (damp::abs(r.Kp) > T{0});
    const bool d_w = (damp::abs(r.c - T{1}) > eps) && (damp::abs(r.Kd) > T{0});
    if (!p_w && !d_w) {
        return false;
    }
    // I-P / I-PD with no integrator cannot hold a setpoint.
    if (!(damp::abs(r.Ki) > T{0}) && p_w && !(damp::abs(r.b) > eps)) {
        return false;
    }
    return true;
}

template<size_t NX, typename T>
[[nodiscard]] constexpr damp::optional<T>
tracking_mag(const damp::StateSpace<NX, 1, 1, T>& sys, const PIDResult<T>& pid, T w) {
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
    if (!(den.norm() > static_cast<T>(1e-30))) {
        return damp::nullopt;
    }
    return damp::abs((pid.eval_ff(w) * P) / den);
}

template<size_t NX, typename T>
[[nodiscard]] constexpr damp::optional<T>
tracking_bandwidth(const damp::StateSpace<NX, 1, 1, T>& sys, const PIDResult<T>& pid, T w_hint) {
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

template<size_t NX, typename T>
[[nodiscard]] constexpr damp::optional<PIDResult<T>>
fit_at_crossover(const damp::StateSpace<NX, 1, 1, T>& sys, const PidTuneSpec<T>& spec) noexcept {
    using Cplx = damp::complex<T>;
    constexpr T pi = damp::numbers::pi_v<T>;
    if (!(spec.wc > T{0})) {
        return damp::nullopt;
    }
    const auto Gopt = damp::eval_frf(sys, Cplx{T{0}, spec.wc});
    if (!Gopt) {
        return damp::nullopt;
    }
    const Cplx G = (*Gopt)(0, 0);
    const T    mag_G = damp::abs(G);
    if (!(mag_G > T{0})) {
        return damp::nullopt;
    }
    const T pm = spec.phase_margin_deg * pi / T{180};
    const T want = -pi + pm - damp::arg(G);
    const T mag_C = T{1} / mag_G;
    const auto [s, c] = damp::sincos(want);
    const T      Cr = mag_C * c;
    const T      Ci = mag_C * s;
    const T      wc = spec.wc;
    PIDResult<T> res{};
    switch (spec.type) {
        case PIDType::P:
            res.Kp = Cr;
            break;
        case PIDType::PI:
            res.Kp = Cr;
            res.Ki = -Ci * wc;
            break;
        case PIDType::PD:
            res.Kp = Cr;
            res.Kd = Ci / wc;
            break;
        case PIDType::PID: {
            res.Kp = Cr;
            const T a = T{4} * wc * wc;
            const T bb = -T{4} * wc * Ci;
            const T cc = -(Cr * Cr);
            const T disc = (bb * bb) - (T{4} * a * cc);
            if (!(disc >= T{0}) || !(a > T{0})) {
                res.Ki = -Ci * wc;
                break;
            }
            const T Kd = (-bb + damp::sqrt(disc)) / (T{2} * a);
            const T Ki = (Cr != T{0} && Kd > T{0}) ? (Cr * Cr) / (T{4} * Kd) : -Ci * wc;
            res.Ki = Ki > T{0} ? Ki : T{0};
            res.Kd = Kd > T{0} ? Kd : T{0};
        } break;
    }
    res.Kbc = res.Ki;
    res.Tf = spec.Tf;
    res.b = spec.b;
    res.c = spec.c;
    return res;
}

/**
 * @brief Raise the loop crossover so T −3 dB matches the 1-DOF target
 *
 * Keeps the sample whose bandwidth is closest to @p target. Each step is
 * clamped to [½, 2] so one update cannot jump the positive-gain ridge. A
 * sample that is negative, and more so than the direct fit, is rejected.
 */
template<size_t NX, typename T>
[[nodiscard]] constexpr PIDResult<T> raise_loop_for_tracking(
    const damp::StateSpace<NX, 1, 1, T>& sys,
    const PidTuneSpec<T>&                spec,
    PIDResult<T>                         fitted,
    T                                    target
) noexcept {
    T            w = spec.wc;
    PIDResult<T> best = fitted;
    bool         scored = false;
    T            best_rel = T{0};
    const auto   bw0 = tracking_bandwidth(sys, fitted, spec.wc);
    if (bw0 && (target > T{0})) {
        best_rel = damp::abs(*bw0 - target) / target;
        scored = true;
    }
    const T wmin = spec.wc / T{4};
    const T wmax = spec.wc * T{64};
    for (int i = 0; i < 8; ++i) {
        PidTuneSpec<T> s = spec;
        s.wc = w;
        const auto r = fit_at_crossover(sys, s);
        if (!r) {
            break;
        }
        // Stop at a gain that is negative and more so than the direct fit.
        // A positive fit must not become negative, and a small negative fit
        // must not run away to a huge one.
        const bool worse_neg = ((r->Kp < T{0}) && (r->Kp < fitted.Kp)) ||
                               ((r->Ki < T{0}) && (r->Ki < fitted.Ki)) ||
                               ((r->Kd < T{0}) && (r->Kd < fitted.Kd));
        if (worse_neg) {
            break;
        }
        const auto bw = tracking_bandwidth(sys, *r, spec.wc);
        if (!bw || !(target > T{0})) {
            break;
        }
        const T rel = damp::abs(*bw - target) / target;
        if (!scored || (rel < best_rel)) {
            best = *r;
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
 * requested law (Ti = 4 Td for PID). For 1-DOF (b = c = 1) that frequency is
 * spec.wc. When b or c is not 1, the loop crossover is raised so the tracking
 * map T_yr = C_ff P / (1 + C_fb P) has the same −3 dB bandwidth as the 1-DOF
 * design at spec.wc. The search keeps the closest match, steps by at most 2x,
 * and does not replace a non-negative fit with a negative-gain sample.
 * Tf is copied onto the result and used in C_ff / C_fb.
 * Returns nullopt if ωc ≤ 0 or P(jωc) is missing / 0.
 *
 * @note Compare with MATLAB®'s pidtune(sys, wc) (that call is 60° 1-DOF PID).
 * @see PidTuneSpec
 */
template<size_t NX, typename T = double>
[[nodiscard]] constexpr damp::optional<PIDResult<T>>
pidtune(const damp::StateSpace<NX, 1, 1, T>& sys, const PidTuneSpec<T>& spec) noexcept {
    const auto fitted = pidtune_detail::fit_at_crossover(sys, spec);
    if (!fitted) {
        return fitted;
    }
    if (!pidtune_detail::weights_move_tracking(*fitted)) {
        return fitted;
    }
    PidTuneSpec<T> spec1 = spec;
    spec1.b = T{1};
    spec1.c = T{1};
    const auto one = pidtune_detail::fit_at_crossover(sys, spec1);
    if (!one) {
        return fitted;
    }
    auto target = pidtune_detail::tracking_bandwidth(sys, *one, spec.wc);
    if (!target) {
        target = spec.wc;
    }
    return pidtune_detail::raise_loop_for_tracking(sys, spec, *fitted, *target);
}

/// @overload ωc only (60° PID).
template<size_t NX, typename T = double>
[[nodiscard]] constexpr damp::optional<PIDResult<T>>
pidtune(const damp::StateSpace<NX, 1, 1, T>& sys, T wc) noexcept {
    return pidtune(sys, PidTuneSpec<T>{.wc = wc});
}

/// @overload ωc and phase margin (PID law).
template<size_t NX, typename T = double>
[[nodiscard]] constexpr damp::optional<PIDResult<T>>
pidtune(const damp::StateSpace<NX, 1, 1, T>& sys, T wc, T phase_margin_deg) noexcept {
    return pidtune(sys, PidTuneSpec<T>{.wc = wc, .phase_margin_deg = phase_margin_deg});
}

/// @overload ωc, phase margin, and law.
template<size_t NX, typename T = double>
[[nodiscard]] constexpr damp::optional<PIDResult<T>>
pidtune(
    const damp::StateSpace<NX, 1, 1, T>& sys,
    T                                    wc,
    T                                    phase_margin_deg,
    PIDType                              type
) noexcept {
    return pidtune(
        sys,
        PidTuneSpec<T>{.wc = wc, .phase_margin_deg = phase_margin_deg, .type = type}
    );
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
    if (os <= T{0}) {
        os = std::numeric_limits<T>::min();
    }

    const T log_os = damp::log(os);
    return -log_os / damp::sqrt((pi * pi) + (log_os * log_os));
}

/**
 * @brief Approximate phase margin from damping ratio
 *
 * Uses the standard second-order approximation:
 *
 *   PM ≈ atan(2*zeta / sqrt(sqrt(1 + 4*zeta^4) - 2*zeta^2))
 *
 * @param zeta Damping ratio
 * @return Approximate phase margin in degrees
 */
template<typename T = double>
[[nodiscard]] constexpr T phase_margin_from_damping_ratio(T zeta) {
    constexpr T pi = damp::numbers::pi_v<T>;

    if (zeta <= T{0}) {
        return T{30};
    }
    if (zeta >= T{1}) {
        return T{85};
    }

    const T z2 = zeta * zeta;
    const T inner = damp::sqrt(T{1} + (T{4} * z2 * z2)) - (T{2} * z2);
    if (inner <= T{0}) {
        return T{60};
    }

    const T pm_rad = damp::atan2(T{2} * zeta, damp::sqrt(inner));
    return pm_rad * T{180} / pi;
}

/**
 * @brief Map settling-time and damping-ratio targets to a bandwidth estimate
 *
 * Uses the 2% settling-time approximation Ts ≈ 4/(zeta*wn) and maps
 * desired bandwidth to wn for a practical one-parameter conversion.
 *
 * @param settling_time Desired 2% settling time in seconds
 * @param zeta          Damping ratio
 * @return Estimated bandwidth in rad/s
 */
template<typename T = double>
[[nodiscard]] constexpr T bandwidth_from_settling_time(T settling_time, T zeta) {
    if (settling_time <= T{0}) {
        return T{1};
    }
    if (zeta <= T{0}) {
        zeta = static_cast<T>(0.7);
    }
    return T{4} / (zeta * settling_time);
}

/**
 * @brief Time-domain performance targets for quick PID synthesis
 */
template<typename T = double>
struct PIDPerformanceSpec {
    T       settling_time{};       ///< Desired 2% settling time [s]
    T       overshoot_percent{};   ///< Desired percent overshoot [0..100)
    T       Ts{};                  ///< Controller sampling time [s]
    PIDType type{PIDType::PID};    ///< PI/PID family selection
    T       bandwidth_scale{T{1}}; ///< Optional tuning scale on estimated bandwidth
};

/**
 * @brief Design PID directly from settling-time and overshoot targets
 *
 * Converts time-domain targets to equivalent damping ratio/phase margin and
 * bandwidth, then calls pid_from_bandwidth().
 *
 * @param spec Performance specification bundle
 * @return PIDResult with tuned gains
 */
template<typename T = double>
[[nodiscard]] constexpr PIDResult<T> pid_from_performance_spec(const PIDPerformanceSpec<T>& spec) {
    const T zeta = damping_ratio_from_overshoot_percent(spec.overshoot_percent);
    const T pm = phase_margin_from_damping_ratio(zeta);
    const T wbw = spec.bandwidth_scale * bandwidth_from_settling_time(spec.settling_time, zeta);
    return pid_from_bandwidth(wbw, pm, spec.Ts, spec.type);
}

// ============================================================================
// Direct Pole Placement PID
// ============================================================================

/**
 * @brief Direct PID pole placement for a first-order-plus-dead-time model
 *
 * Places the closed-loop poles at the specified locations by computing PID
 * gains for a discretized FOPDT plant model.
 *
 * The FOPDT model is: G(s) = K * exp(-Ls) / (tau*s + 1)
 * Discretized with ZOH: G(z) = K*(1 - a) / (z - a)  where a = exp(-Ts/tau)
 * Dead time approximated as d = round(L/Ts) pure delays: G(z) = K*(1-a)*z^{-d} / (z - a)
 *
 * For a PI controller with one integrator pole at z=1, the closed-loop
 * characteristic polynomial is degree 2+d. The user specifies the desired
 * closed-loop poles and the gains are computed.
 *
 * This function handles the d=0 (no dead-time) case, placing closed-loop
 * poles for the 2nd-order closed-loop system (plant pole + integrator).
 *
 * @param K     Static gain
 * @param tau   Time constant (seconds)
 * @param p1    Desired closed-loop pole 1 (z-domain, inside unit circle)
 * @param p2    Desired closed-loop pole 2 (z-domain, inside unit circle)
 * @param Ts    Sampling time (seconds)
 * @return PIDResult (PI gains via pole placement); zeros if tau ≤ 0, Ts ≤ 0, or b = 0
 */
template<typename T = double>
[[nodiscard]] constexpr PIDResult<T>
pid_pole_placement(T K, T tau, T p1, T p2, T Ts) {
    if (!(tau > T{0}) || !(Ts > T{0}) || K == T{0}) {
        return PIDResult<T>{};
    }

    // Discretize the first-order plant: G(z) = K*(1-a)/(z-a)
    T a = damp::exp(-Ts / tau);
    T b = K * (T{1} - a);
    if (b == T{0}) {
        return PIDResult<T>{};
    }

    // Deploy is PIDResult::discretize + PIDController: I ← I + (Ki Ts) e[k],
    // so C(z) = Kp + Ki Ts · z/(z−1) (backward Euler). Match
    //   (z−a)(z−1) + b (Kp (z−1) + Ki Ts z) = (z−p1)(z−p2).
    // z² + (b Kp + b Ki Ts − a − 1) z + (a − b Kp) = z² − (p1+p2) z + p1 p2
    //   a − b Kp = p1 p2                 => Kp = (a − p1 p2) / b
    //   b Kp + b Ki Ts = 1+a − (p1+p2)   => Ki from the leftover.

    const T Kp = (a - (p1 * p2)) / b;
    const T Ki = ((T{1} + a) - (p1 + p2) - (b * Kp)) / (b * Ts);

    return PIDResult<T>{Kp, Ki, T{0}};
}

/**
 * @brief PID pole placement for first-order plant with 3 desired poles
 *
 * For PID control of a first-order plant (no dead time), the closed-loop
 * has 3 poles (plant + integrator + derivative filter).
 *
 * @param K     Static gain
 * @param tau   Time constant (seconds)
 * @param p1    Desired closed-loop pole 1 (z-domain)
 * @param p2    Desired closed-loop pole 2 (z-domain)
 * @param p3    Desired closed-loop pole 3 (z-domain)
 * @param Ts    Sampling time (seconds)
 * @return PIDResult with PID gains; zeros if tau ≤ 0, Ts ≤ 0, or b = 0
 */
template<typename T = double>
[[nodiscard]] constexpr PIDResult<T>
pid_pole_placement(T K, T tau, T p1, T p2, T p3, T Ts) {
    if (!(tau > T{0}) || !(Ts > T{0}) || K == T{0}) {
        return PIDResult<T>{};
    }

    // Discretize: G(z) = b/(z-a), a = exp(-Ts/tau), b = K*(1-a)
    T a = damp::exp(-Ts / tau);
    T b = K * (T{1} - a);
    if (b == T{0}) {
        return PIDResult<T>{};
    }

    // Deploy (Tf = 0): C(z) = Kp + Ki Ts z/(z−1) + (Kd/Ts)(z−1)/z.
    // Closed-loop:
    //   z(z−1)(z−a) + b [Kp z(z−1) + Ki Ts z² + (Kd/Ts)(z−1)²]
    //     = z³ + (−(1+a) + b Kp + b Ki Ts + b Kd/Ts) z²
    //         + (a − b Kp − 2 b Kd/Ts) z + b Kd/Ts
    // Desired: (z−p1)(z−p2)(z−p3) = z³ − σ z² + ρ z − π
    const T sigma = p1 + p2 + p3;
    const T rho = (p1 * p2) + (p1 * p3) + (p2 * p3);
    const T prod = p1 * p2 * p3;

    // b Kd/Ts = −π
    const T Kd = -(prod * Ts) / b;
    // b Kp + 2 b Kd/Ts = a − ρ
    const T Kp = (a - rho - (T{2} * b * Kd / Ts)) / b;
    // b Kp + b Ki Ts + b Kd/Ts = 1 + a − σ
    const T Ki = ((T{1} + a - sigma) - (b * Kp) - (b * Kd / Ts)) / (b * Ts);

    return PIDResult<T>{Kp, Ki, Kd};
}

// ============================================================================
// Continuous PI Pole Placement (first-order plant)
// ============================================================================

/**
 * @brief PI gains that place the closed-loop poles of a first-order plant
 *
 * The SISO pole-placement special case in which the controller structure is
 * fixed to a PI. A PI has two free gains, so it can place exactly two closed-loop
 * poles — which is well-posed precisely when the plant is first-order
 * @f$ G(s) = 1/(a_1 s + a_0) @f$ (the closed loop is then 2nd order). Matching
 * @f$ s^2 + \frac{a_0+K_p}{a_1}s + \frac{K_i}{a_1} @f$ to
 * @f$ s^2 + 2\zeta\omega_n s + \omega_n^2 @f$ gives the closed form
 * @f[
 *   K_p = 2\zeta\omega_n a_1 - a_0, \qquad K_i = a_1\,\omega_n^2 ,
 * @f]
 * so @p omega_bw is the closed-loop bandwidth and @f$ \zeta = 1 @f$ (default)
 * places a critically damped double pole at @f$ -\omega_n @f$. This is the kernel
 * every loop in a motion cascade shares — only the plant changes:
 * @f$ (a_1,a_0)=(L,R) @f$ for a dq current axis,
 * @f$ (J,b_{visc}) @f$ for a velocity loop on the inertia plant.
 *
 * Equivalent to running Ackermann (@ref place) on the *(plant + integrator)*
 * augmented single-input system, but the first-order case is analytic. The
 * back-calculation gain is seeded to @f$ K_p @f$ (tracking time constant
 * @f$ T_t = K_p/K_i @f$); the setpoint weight @p b leaves the poles fixed and
 * only moves the reference→output zero (@p b = 1 PI on error, @p b = 0 I-P).
 *
 * @note @p omega_bw is bounded by the eventual sample rate — keep it roughly an
 *       order of magnitude below the control frequency for the continuous
 *       placement to hold.
 *
 * @see place — general MIMO pole placement (full-state feedback A−BK).
 * @see Åström & Hägglund, "Advanced PID Control", 2006, §4.4 — setpoint weighting.
 *
 * @param a1       Plant denominator first-order coefficient (e.g. L, or J)
 * @param a0       Plant denominator zeroth-order coefficient (e.g. R, or viscous b)
 * @param omega_bw [rad/s]  Desired closed-loop pole frequency (bandwidth)
 * @param zeta     [-]      Closed-loop damping ratio (default 1, critically damped)
 * @param b        [-]      Proportional setpoint weight: 1 = PI, 0 = I-P
 * @return PIDResult with Kp, Ki, Kbc (= Kp) and the setpoint weight b set
 */
template<typename T = double>
[[nodiscard]] constexpr PIDResult<T>
pi_pole_placement_first_order(T a1, T a0, T omega_bw, T zeta = T{1}, T b = T{1}) {
    // Parallel form: T_t = Kbc/Ki. Setting Kbc = Kp ⇒ T_t = Kp/Ki = T_i (Åström default).
    // b: 1 = PI (P on error), 0 = I-P (P on measurement)
    const T Kp = (T{2} * zeta * omega_bw * a1) - a0;
    const T Ki = a1 * omega_bw * omega_bw;
    return PIDResult<T>{
        .Kp = Kp,
        .Ki = Ki,
        .Kbc = Kp,
        .b = b,
    };
}

/**
 * @brief PI pole placement from a first-order plant transfer function
 *
 * Convenience overload for callers working in the systems layer: takes the plant
 * @f$ G(s) = b_0/(a_1 s + a_0) @f$ as a first-order TransferFunction
 * (`num = {b₀}`, `den = {a₀, a₁}` in ascending powers) and reduces it to the monic
 * form by dividing through by the numerator gain @f$ b_0 @f$, then delegates to the
 * scalar @ref pi_pole_placement_first_order. The static `<1, 2>` shape constrains
 * the plant to first order with no zero — the only case a PI can pole-place.
 *
 * @note A zero numerator gain @f$ b_0 = 0 @f$ is a degenerate plant (no input
 *       authority — a default-constructed TransferFunction is one), so this
 *       returns an inert all-zero PI rather than dividing through by zero. Unlike
 *       the rest of this header's unchecked formulas, the guard is here because
 *       @f$ b_0 @f$ is a caller-supplied plant property that is easy to leave unset.
 *
 * @param plant    First-order plant b₀/(a₁ s + a₀)
 * @param omega_bw [rad/s]  Desired closed-loop bandwidth
 * @param zeta     [-]      Closed-loop damping ratio (default 1)
 * @param b        [-]      Proportional setpoint weight (default 1 = PI)
 * @return PIDResult with the placed PI gains, or all-zero gains if b₀ = 0
 */
template<typename T = double>
[[nodiscard]] constexpr PIDResult<T>
pi_pole_placement_first_order(const TransferFunction<1, 2, T>& plant, T omega_bw, T zeta = T{1}, T b = T{1}) {
    const T b0 = plant.num[0];
    if (b0 == T{0}) {
        return PIDResult<T>{}; // degenerate plant (no input gain) -> inert PI
    }
    return pi_pole_placement_first_order(plant.den[1] / b0, plant.den[0] / b0, omega_bw, zeta, b);
}

/**
 * @brief PID (or PD) pole placement for a rigid inertia @f$ G(s) = 1/(J s^2) @f$
 *
 * Sibling of @ref pi_pole_placement_first_order for the other everyday motion
 * plant: a double integrator (position loop on a flywheel / mass, no rate
 * measurement). A PD has two free gains, so it can place the two closed-loop
 * poles of @f$ J\ddot\theta = u @f$. Matching
 * @f$ s^2 + (K_d/J)s + K_p/J @f$ to @f$ s^2 + 2\zeta\omega s + \omega^2 @f$
 * gives
 * @f[
 *   K_p = J\omega^2,\qquad K_d = 2\zeta\omega J .
 * @f]
 * Optional @p omega_i > 0 adds an integrator and places the third pole at
 * @f$ -\omega_i @f$:
 * @f[
 *   (s^2 + 2\zeta\omega s + \omega^2)(s + \omega_i).
 * @f]
 * A common load-holding choice is @p omega_i = @p omega_bw / 8 (I a decade
 * slower than the PD pair). @p omega_i = 0 is pure PD.
 *
 * Derivative is on the measurement (@f$ c = 0 @f$, no D kick on a reference
 * step). @f$ T_f = 1/(10\omega) @f$ when there is a D term. Back-calculation
 * is seeded to @f$ K_p @f$ (same @f$ T_t = T_i @f$ default as the first-order
 * PI helper).
 *
 * @note Compare with MATLAB®'s pidtune on tf(1,[J 0 0]) for the PD part.
 * @see pi_pole_placement_first_order — velocity / current (one integrator)
 * @see "Feedback Control of Dynamic Systems" (Franklin et al.), PD on 1/s²
 *
 * @param J        [kg·m²] inertia (> 0)
 * @param omega_bw [rad/s] desired PD-pair frequency (> 0)
 * @param zeta     [-]     PD-pair damping (default 1, critically damped)
 * @param omega_i  [rad/s] integrator pole; 0 = PD only (default)
 * @return PIDResult with Kp, Ki, Kd, Kbc, c, Tf; zeros if J, ω, or ζ is not
 *         strictly positive, or if ωᵢ < 0
 */
template<typename T = double>
[[nodiscard]] constexpr PIDResult<T> pid_pole_placement_double_integrator(
    T J, T omega_bw, T zeta = T{1}, T omega_i = T{0}
) {
    if (!(J > T{0}) || !(omega_bw > T{0}) || !(zeta > T{0}) || omega_i < T{0}) {
        return PIDResult<T>{};
    }

    const T w = omega_bw;
    const T two_zw = T{2} * zeta * w;
    T       Kp = J * (w * w);
    T       Kd = J * two_zw;
    T       Ki = T{0};
    if (omega_i > T{0}) {
        // (s² + 2ζω s + ω²)(s + ωᵢ) = s³ + (2ζω+ωᵢ)s² + (ω²+2ζω ωᵢ)s + ω² ωᵢ
        Kp = J * ((w * w) + (two_zw * omega_i));
        Kd = J * (two_zw + omega_i);
        Ki = J * (w * w) * omega_i;
    }

    PIDResult<T> r{};
    r.Kp = Kp;
    r.Ki = Ki;
    r.Kd = Kd;
    r.Kbc = Kp;
    r.c = T{0}; // PI-D / PD on measurement
    r.Tf = (Kd > T{0}) ? (T{1} / (T{10} * w)) : T{0};
    return r;
}

/**
 * @brief Fixed-rate inductor current-loop PI (topology-agnostic)
 *
 * Places a continuous PI on residual plant @f$ 1/(sL + R_L) @f$ after volt-second
 * feedforward, clamps the continuous voltage command limits, and discretizes at
 * @p Ts. Shared by every DC–DC / PFC average-current stage — not named for a
 * topology. Pair the discrete PID with @ref power::CurrentStage and a map
 * (@ref power::BoostMap, @ref power::BuckMap, …).
 *
 * @param L        [H] inductance (> 0)
 * @param omega_i  [rad/s] current-loop bandwidth (> 0)
 * @param Ts       [s] control / PWM period (> 0)
 * @param R_L      [ohm] series resistance (default 0)
 * @param v_limit  [V] |v*| PI clamp (0 → loose default L·ωᵢ·100)
 * @param zeta     damping of pole placement (default 1)
 * @return @c InductorCurrentPIResult with discrete PID and @c success
 *
 * @see pi_pole_placement_first_order
 * @see power::CurrentStage
 */
template<typename T = double>
struct InductorCurrentPIResult {
    DiscretePIDResult<T> pid{};          ///< Discrete PI ready for PIController / CurrentStage
    bool                 success{false}; ///< true if placement and discretize succeeded

    template<typename U>
    [[nodiscard]] constexpr InductorCurrentPIResult<U> as() const {
        return InductorCurrentPIResult<U>{pid.template as<U>(), success};
    }
};

template<typename T = double>
[[nodiscard]] constexpr InductorCurrentPIResult<T> inductor_current_pi(
    T L, T omega_i, T Ts, T R_L = T{0}, T v_limit = T{0}, T zeta = T{1}
) {
    InductorCurrentPIResult<T> r{};
    if (!(L > T{0}) || !(omega_i > T{0}) || !(Ts > T{0}) || R_L < T{0}) {
        return r;
    }

    auto    cont = pi_pole_placement_first_order(L, R_L, omega_i, zeta, T{1});
    const T lim = (v_limit > T{0}) ? v_limit : (L * omega_i * T{100});
    cont.u_min = -lim;
    cont.u_max = lim;
    if (!(cont.Kp > T{0}) || !(cont.Ki > T{0})) {
        return r;
    }

    r.pid = cont.discretize(Ts);
    r.success = true;
    return r;
}

} // namespace design
} // namespace damp
