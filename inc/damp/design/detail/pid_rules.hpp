// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
// IWYU pragma: private, include "damp/design/pid_design.hpp"

/**
 * @file pid_rules.hpp
 * @brief The math of design::pid's rules: SIMC, AMIGO, the loop-shaping fit and its certificate,
 *        and the pole placements, each from a model's numbers to continuous gains or the reason
 *        there are none
 *
 * design::pid (design/pid_design.hpp) reads the model, the rule and the options, calls one of
 * these, and deploys the gains at Ts with the options' filter, weights, tracking time and limits.
 * Every function here takes scalars and systems, not the public model and rule types, and
 * refuses what it cannot tune through PIDTuning's status, with zero gains.
 */

#include "damp/backend.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/design/status.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/systems/zpk.hpp"

// Arduino.h defines PI as a macro, which would rewrite PIDMode::PI: set aside for this header.
#pragma push_macro("PI")
#undef PI

namespace damp::design::detail {

/** @brief A rule's continuous gains, or why it has none (zero gains). */
template<typename T>
struct PIDTuning {
    PIDGains<T>  gains{};                  ///< Kp, Ki, Kd, and the Tf, b, c a rule sets itself
    bool         success{true};            ///< status == DesignStatus::ok
    DesignStatus status{DesignStatus::ok}; ///< ok, or why the rule refused the model
};

/** @brief A rule's gains, accepted. */
template<typename T>
[[nodiscard]] constexpr PIDTuning<T> tuned(const PIDGains<T>& gains) {
    return PIDTuning<T>{.gains = gains};
}

/** @brief Refused because a model or rule number is outside its domain (K = 0, L < 0, τ ≤ 0, …). */
template<typename T>
[[nodiscard]] constexpr PIDTuning<T> pid_invalid() {
    return refused<PIDTuning<T>>(DesignStatus::invalid_parameter);
}

/** @brief A static gain K that is finite and nonzero (either sign). */
template<typename T>
[[nodiscard]] constexpr bool valid_static_gain(T K) {
    return damp::isfinite(K) && (K != T{0});
}

// ============================================================================
// SIMC
// ============================================================================

/**
 * @brief SIMC gains on the FOPDT model K·e^{−Ls}/(τs + 1)
 *
 * The SIMC PI rule, from the direct synthesis of a first-order closed loop e^{−Ls}/(τ_c s + 1):
 * @f[
 *   K_p = \frac{\tau}{K(\tau_c + L)}, \qquad T_i = \min\big(\tau,\ 4(\tau_c + L)\big), \qquad K_i = K_p/T_i .
 * @f]
 * τ_c is the one tuning number: τ_c = L is the rule's own choice, tight control with good
 * robustness, and a larger τ_c is slower and more robust. The min keeps the integral fast enough
 * to reject an input disturbance on a lag-dominant plant (τ ≫ L). P mode is the PI's Kp without
 * the integrator. The derivative of the PD and PID modes, K_d = K_p·L/2, is not SIMC's (which adds
 * τ_D = τ₂ only for a second-order plant): it is a dead-time heuristic, and 0 at L = 0.
 *
 * Refused with invalid_parameter if K is 0 or not finite, L < 0, τ ≤ 0, τ_c + L ≤ 0, or a number
 * is not finite.
 *
 * @see "Simple analytic rules for model reduction and PID controller tuning" (Skogestad, 2003),
 *      J. Process Control 13(4)
 */
template<typename T>
[[nodiscard]] constexpr PIDTuning<T> simc_gains(T K, T L, T tau, T tau_c, PIDMode mode) {
    const T den = tau_c + L;
    if (!valid_static_gain(K) || !damp::finite_non_negative(L) || !damp::finite_positive(tau) || !damp::isfinite(tau_c)
        || !damp::finite_positive(den)) {
        return pid_invalid<T>();
    }
    PIDGains<T> g{};
    g.Kp = tau / (K * den);
    if (mode == PIDMode::PI || mode == PIDMode::PID) {
        const T Ti = (tau < T{4} * den) ? tau : T{4} * den; // both > 0
        g.Ki = g.Kp / Ti;
    }
    if (mode == PIDMode::PD || mode == PIDMode::PID) {
        g.Kd = g.Kp * L / T{2};
    }
    if (!damp::isfinite(g.Kp) || !damp::isfinite(g.Ki) || !damp::isfinite(g.Kd)) {
        return pid_invalid<T>();
    }
    return tuned(g);
}

// ============================================================================
// AMIGO
// ============================================================================

/**
 * @brief The AMIGO PI on the FOPDT model K·e^{−Ls}/(τs + 1)
 *
 * The PI rule Åström and Hägglund fitted to designs that maximize integral gain under a
 * robustness constraint (M_s ≈ 1.4) over a test batch of process models:
 * @f[
 *   K_p = \frac{0.15}{K} + \Big(0.35 - \frac{L\tau}{(L + \tau)^2}\Big)\frac{\tau}{K L}, \qquad
 *   T_i = 0.35 L + \frac{13 L \tau^2}{\tau^2 + 12 L\tau + 7 L^2}, \qquad K_i = K_p/T_i .
 * @f]
 * The bracket is at least 0.1 (Lτ/(L + τ)² ≤ ¼), so K_p and K_i take K's sign, and T_i > 0: a
 * reverse-acting plant (K < 0) gets the negated gains. The rule divides by L, so a plant with no
 * dead time has no AMIGO PI.
 *
 * Refused with invalid_parameter if K is 0 or not finite, L ≤ 0, τ ≤ 0, or a number is not finite.
 *
 * @see "Revisiting the Ziegler–Nichols step response method for PID control" (Åström & Hägglund,
 *      2004), J. Process Control 14(6)
 * @see "Advanced PID Control" (Åström & Hägglund, 2006), §7.6 (AMIGO PI)
 */
template<typename T>
[[nodiscard]] constexpr PIDTuning<T> amigo_pi(T K, T L, T tau) {
    if (!valid_static_gain(K) || !damp::finite_positive(L) || !damp::finite_positive(tau)) {
        return pid_invalid<T>();
    }
    const T LpT = L + tau;
    const T Kp = (static_cast<T>(0.15) / K) + ((static_cast<T>(0.35) - (L * tau / (LpT * LpT))) * tau / (K * L));
    const T Ti = (static_cast<T>(0.35) * L) + (T{13} * L * tau * tau / ((tau * tau) + (T{12} * L * tau) + (T{7} * L * L)));
    if (!damp::isfinite(Kp) || Kp == T{0} || !damp::finite_positive(Ti)) {
        return pid_invalid<T>();
    }
    const T Ki = Kp / Ti;
    if (!damp::isfinite(Ki)) {
        return pid_invalid<T>();
    }
    return tuned(PIDGains<T>{.Kp = Kp, .Ki = Ki});
}

/**
 * @brief The AMIGO PI from a relay test's ultimate point (K_u, T_u) and the static gain K
 *
 * Recovers the FOPDT model that has the same static gain and the same ultimate point, then applies
 * amigo_pi. With ω_u = 2π/T_u, |P(jω_u)| = K/√(1 + ω_u²τ²) = 1/K_u and arg P(jω_u) = −π give
 * @f[
 *   \tau = \frac{\sqrt{(K K_u)^2 - 1}}{\omega_u}, \qquad L = \frac{\pi - \arctan(\omega_u \tau)}{\omega_u} .
 * @f]
 * A relay test measures |P|, so the plant is taken as direct-acting (K > 0).
 *
 * Refused with invalid_parameter if K, K_u or T_u is not finite and positive, or K·K_u ≤ 1 (no
 * FOPDT has that ultimate point: K = 0, the static gain not known, is refused here).
 *
 * @see "Advanced PID Control" (Åström & Hägglund, 2006), §7.6 (AMIGO PI)
 */
template<typename T>
[[nodiscard]] constexpr PIDTuning<T> amigo_from_ultimate(T K, T Ku, T Tu) {
    if (!damp::finite_positive(K) || !damp::finite_positive(Ku) || !damp::finite_positive(Tu)) {
        return pid_invalid<T>();
    }
    const T wu = (T{2} * damp::numbers::pi_v<T>) / Tu;
    const T prod = K * Ku;
    if (!damp::isfinite(prod) || !(prod > T{1})) {
        return pid_invalid<T>(); // no FOPDT has this ultimate point
    }
    const T tau = damp::sqrt((prod * prod) - T{1}) / wu;
    const T L = (damp::numbers::pi_v<T> - damp::atan(wu * tau)) / wu;
    return amigo_pi(K, L, tau);
}

// ============================================================================
// Loop shaping: a crossover and a phase margin
// ============================================================================

/** @brief What a LoopShape asks of the loop, with the options it reads: the derivative filter and the setpoint weights. */
template<typename T>
struct LoopSpec {
    T       wc{};                    ///< [rad/s] Loop crossover ωc (the 1-DOF design's)
    T       phase_margin_deg{T{60}}; ///< [deg] Phase margin at the crossover, in (0, 180)
    PIDMode mode{PIDMode::PID};      ///< P / PI / PD / PID
    T       Tf{T{0}};                ///< [s] Derivative filter (0: none)
    T       b{T{1}};                 ///< Proportional setpoint weight
    T       c{T{1}};                 ///< Derivative setpoint weight
};

/** @brief The spec's scalar domain: 0 < PM < 180°, Tf ≥ 0, b and c finite (ωc is the caller's to check). */
template<typename T>
[[nodiscard]] constexpr bool valid_loop_spec(const LoopSpec<T>& spec) {
    return damp::finite_positive(spec.phase_margin_deg) && (spec.phase_margin_deg < T{180}) && damp::finite_non_negative(spec.Tf)
        && damp::isfinite(spec.b) && damp::isfinite(spec.c);
}

/** @brief Whether the setpoint weights change the tracking map T_yr, so the crossover is raised for it. */
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

/** @brief |T_yr(jω)|, T_yr = C_ff P/(1 + C_fb P), or nullopt where P is missing or a closed-loop pole sits on jω. */
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

/** @brief The −3 dB bandwidth of T_yr, by bisection on a log grid about @p w_hint, or nullopt if it has none there. */
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
 * @brief The law of @p spec with C(jω) = e^{j(−π + PM)}/P(jω), from the plant's response P at ω = spec.wc
 *
 * The loop C·P then crosses the unit circle at ω with the phase margin PM: the plant's point at ω
 * is moved onto the unit circle at −π + PM. With the ultimate point P(jω_u) = −1/K_u this is the
 * phase-margin design of Åström and Hägglund's relay autotuner. The filtered derivative K_d·g,
 * g = jω/(1 + jωTf), enters the fit, so the delivered C_fb meets the crossover and the phase margin
 * for any Tf. With C = C_r + jC_i and g = g_r + j g_i:
 * - PI: K_p = C_r, Kᵢ = −C_i ω.
 * - PD: K_d = C_i/g_i, K_p = C_r − K_d g_r.
 * - PID, with Tᵢ = 4T_d (Kᵢ = K_p²/(4K_d)): (4ω g_i − g_r²)K_d² + (2C_r g_r − 4ω C_i)K_d − C_r² = 0,
 *   its positive root, and K_p = C_r − K_d g_r.
 * - P meets the crossover alone, K_p = 1/|P|; its phase margin is the plant's.
 *
 * Refused with singular if P is zero or not finite, and with not_placeable if the law cannot
 * realize C with positive gains (a PI cannot lead, a PD cannot lag, no law gives K_p ≤ 0).
 *
 * @see "Automatic tuning of simple regulators with specifications on phase and amplitude margins"
 *      (Åström & Hägglund, 1984), Automatica 20(5)
 */
template<typename T>
[[nodiscard]] constexpr PIDTuning<T> fit_to_response(damp::complex<T> P, const LoopSpec<T>& spec) noexcept {
    constexpr T pi = damp::numbers::pi_v<T>;
    const T     mag_P = damp::abs(P);
    if (!damp::finite_positive(mag_P)) {
        return refused<PIDTuning<T>>(DesignStatus::singular);
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
    if (!damp::finite_positive(res.Kp) || !damp::finite_non_negative(res.Ki) || !damp::finite_non_negative(res.Kd)
        || ((spec.mode == PIDMode::PID) && !(res.Kd > T{0}))) {
        return refused<PIDTuning<T>>(DesignStatus::not_placeable);
    }
    res.Tf = spec.Tf;
    res.b = spec.b;
    res.c = spec.c;
    return tuned(res);
}

/**
 * @brief fit_to_response at spec.wc with every gain of sign σ = ±1
 *
 * σ = −1 fits σ·P and negates the gains: the loop C·P is the same.
 */
template<size_t NX, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr PIDTuning<T>
fit_at_crossover(const damp::StateSpace<NX, 1, 1, T, NW, NV>& sys, const LoopSpec<T>& spec, T sigma) noexcept {
    const auto Gopt = damp::eval_frf(sys, damp::complex<T>{T{0}, spec.wc});
    if (!Gopt) {
        return refused<PIDTuning<T>>(DesignStatus::singular);
    }
    auto g = fit_to_response((*Gopt)(0, 0) * sigma, spec);
    if (g.success && (sigma < T{0})) {
        g.gains.Kp = -g.gains.Kp;
        g.gains.Ki = -g.gains.Ki;
        g.gains.Kd = -g.gains.Kd;
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
    if (!damp::finite_positive(ratio)) {
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
[[nodiscard]] constexpr PIDTuning<T> raise_loop_for_tracking(
    const damp::StateSpace<NX, 1, 1, T, NW, NV>& sys,
    const TransferFunction<NX + 1, NX + 1, T>&   plant,
    const LoopSpec<T>&                           spec,
    PIDTuning<T>                                 fitted,
    T                                            target,
    T                                            sigma
) noexcept {
    T            w = spec.wc;
    PIDTuning<T> best = fitted;
    bool         scored = false;
    T            best_rel = T{0};
    const auto   bw0 = tracking_bandwidth(sys, fitted.gains, spec.wc);
    if (bw0 && (target > T{0})) {
        best_rel = damp::abs(*bw0 - target) / target;
        scored = true;
    }
    const T wmin = spec.wc / T{4};
    const T wmax = spec.wc * T{64};
    for (int i = 0; i < 8; ++i) {
        LoopSpec<T> s = spec;
        s.wc = w;
        const auto r = fit_at_crossover(sys, s, sigma);
        if (!r.success || !closed_loop_stable(plant, r.gains)) {
            break;
        }
        const auto bw = tracking_bandwidth(sys, r.gains, spec.wc);
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
        const T next = damp::clamp(w * damp::clamp(target / *bw, T{0.5}, T{2}), wmin, wmax);
        if (!(damp::abs(next - w) > (w * static_cast<T>(1e-4)))) {
            break;
        }
        w = next;
    }
    return best;
}

/**
 * @brief The loop-shaping law on a continuous SISO plant: the crossover and phase margin of @p spec, certified
 *
 * Chooses C_fb(jω) so |C_fb P| = 1 and arg(C_fb P) = −π + PM (fit_to_response), the PID with
 * Tᵢ = 4T_d. For 1-DOF (b = c = 1) that frequency is spec.wc. When b or c is not 1, the loop
 * crossover is raised so the tracking map T_yr = C_ff P/(1 + C_fb P) has the same −3 dB bandwidth
 * as the 1-DOF design at spec.wc; the search keeps the closest match and steps by at most 2×.
 *
 * The gains share one sign. The sign of the plant's low-frequency gain num(0) is tried first, so a
 * reverse-acting plant gets K_p, Kᵢ, K_d ≤ 0 (the fit on −P is the negated fit on P), and the
 * other sign is taken only when it alone gives a stable loop.
 *
 * The result is certified: every pole of the closed loop 1 + C_fb P lies in the open left half
 * plane (Routh–Hurwitz on the characteristic polynomial from ss2tf).
 *
 * Refused with invalid_parameter if ωc ≤ 0, PM ∉ (0°, 180°), Tf < 0 or a number is not finite;
 * invalid_sample_time on a discrete plant; singular if P(jωc) is missing or 0; not_placeable if the
 * law cannot realize the crossover with gains of one sign; not_stabilizable if the fitted loop is
 * unstable. When neither sign gives a stable loop, the status is that of the plant's sign.
 */
template<size_t NX, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr PIDTuning<T> loop_shape(const damp::StateSpace<NX, 1, 1, T, NW, NV>& sys, const LoopSpec<T>& spec) noexcept {
    if (!damp::finite_positive(spec.wc) || !valid_loop_spec(spec)) {
        return pid_invalid<T>();
    }
    if (!damp::isfinite(sys.Ts) || sys.Ts != T{0}) {
        return refused<PIDTuning<T>>(DesignStatus::invalid_sample_time);
    }
    const auto              plant = ::damp::ss2tf(sys);
    const T                 first = (plant.num[NX] < T{0}) ? T{-1} : T{1}; // sign of num(0), den monic
    const damp::array<T, 2> signs{first, -first};
    PIDTuning<T>            fitted{};
    T                       sigma = first;
    for (const T s : signs) {
        auto g = fit_at_crossover(sys, spec, s);
        if (g.success && !closed_loop_stable(plant, g.gains)) {
            g = refused<PIDTuning<T>>(DesignStatus::not_stabilizable);
        }
        if (g.success || (s == first)) { // on a double failure, report the plant sign's status
            fitted = g;
            sigma = s;
        }
        if (g.success) {
            break;
        }
    }
    if (!fitted.success || !weights_move_tracking(fitted.gains)) {
        return fitted;
    }
    LoopSpec<T> spec1 = spec;
    spec1.b = T{1};
    spec1.c = T{1};
    const auto one_dof = fit_at_crossover(sys, spec1, sigma);
    auto       target = tracking_bandwidth(sys, one_dof.gains, spec.wc);
    if (!target) {
        target = spec.wc;
    }
    return raise_loop_for_tracking(sys, plant, spec, fitted, *target, sigma);
}

/**
 * @brief The loop-shaping law at a relay test's ultimate point: P(jω_u) = −1/K_u at ω_u = 2π/T_u
 *
 * fit_to_response at the one point of the plant the test measured, the crossover there: the
 * phase-margin design of Åström and Hägglund's autotuner, which moves the ultimate point onto the
 * unit circle at −π + PM. The point is all of the plant the design has, so the closed loop is not
 * certified, the setpoint weights do not move the crossover (they are the law's as given), and the
 * plant is taken as direct-acting (a relay measures |P|). At −π a PI cannot add the lead
 * (not_placeable), and a P is K_u itself, the loop on the stability boundary (not_stabilizable):
 * the useful laws are a PD and a PID.
 *
 * Refused with invalid_parameter if K_u or T_u is not finite and positive, spec.wc is not 0 (the
 * crossover is the ultimate frequency), PM ∉ (0°, 180°), Tf < 0 or a number is not finite.
 *
 * @see "Automatic tuning of simple regulators with specifications on phase and amplitude margins"
 *      (Åström & Hägglund, 1984), Automatica 20(5)
 */
template<typename T>
[[nodiscard]] constexpr PIDTuning<T> loop_shape_ultimate(T Ku, T Tu, LoopSpec<T> spec) noexcept {
    if (!damp::finite_positive(Ku) || !damp::finite_positive(Tu) || spec.wc != T{0} || !valid_loop_spec(spec)) {
        return pid_invalid<T>();
    }
    if (spec.mode == PIDMode::P) {
        return refused<PIDTuning<T>>(DesignStatus::not_stabilizable); // K_u: the loop at the ultimate point
    }
    spec.wc = (T{2} * damp::numbers::pi_v<T>) / Tu;
    return fit_to_response(damp::complex<T>{-T{1} / Ku, T{0}}, spec);
}

// ============================================================================
// Pole placement
// ============================================================================

/**
 * @brief The PI that places the closed-loop pair s² + 2ζω s + ω² of a first-order plant b₀/(a₁ s + a₀)
 *
 * A PI has two free gains, so it places exactly the two closed-loop poles of a first-order plant.
 * With the plant made monic by its gain, a₁′ = a₁/b₀, a₀′ = a₀/b₀, matching
 * @f$ s^2 + \frac{a_0'+K_p}{a_1'}s + \frac{K_i}{a_1'} @f$ to @f$ s^2 + 2\zeta\omega s + \omega^2 @f$ gives
 * @f[
 *   K_p = 2\zeta\omega a_1' - a_0', \qquad K_i = a_1'\,\omega^2 .
 * @f]
 * This is the kernel every loop in a motion cascade shares: (a₁, a₀) = (L, R) for a dq current
 * axis, (J, b_visc) for a velocity loop. The setpoint weight b leaves the poles fixed and only moves
 * the reference→output zero (b = 1 PI on error, b = 0 I-P). ω is bounded by the sample rate: keep it
 * roughly an order of magnitude below it for the continuous placement to hold at Ts.
 *
 * Refused with invalid_parameter if b₀ = 0, a₁/b₀ ≤ 0, ω ≤ 0, ζ ≤ 0, ωᵢ ≠ 0 (a PI has no third
 * pole) or a number is not finite; not_placeable if the placement needs K_p ≤ 0 (2ζω a₁′ ≤ a₀′: the
 * plant pole is already faster than the requested pair).
 *
 * @see "Advanced PID Control" (Åström & Hägglund, 2006), §4.4 (setpoint weighting)
 */
template<typename T>
[[nodiscard]] constexpr PIDTuning<T> pole_pair_pi(const TransferFunction<1, 2, T>& plant, T omega, T zeta, T omega_i) {
    const T b0 = plant.num[0];
    if (!valid_static_gain(b0) || !damp::finite_positive(omega) || !damp::finite_positive(zeta) || !damp::isfinite(omega_i)
        || omega_i != T{0}) {
        return pid_invalid<T>();
    }
    const T a1 = plant.den[0] / b0;
    const T a0 = plant.den[1] / b0;
    if (!damp::finite_positive(a1) || !damp::isfinite(a0)) {
        return pid_invalid<T>();
    }
    const T Kp = (T{2} * zeta * omega * a1) - a0;
    if (!damp::finite_positive(Kp)) {
        return refused<PIDTuning<T>>(DesignStatus::not_placeable);
    }
    return tuned(PIDGains<T>{.Kp = Kp, .Ki = a1 * omega * omega});
}

/**
 * @brief The PD or PID that places the closed-loop poles of a second-order plant b₀/(a₂ s² + a₁ s + a₀)
 *
 * The derivative is filtered, K_d s/(1 + Tf s), and the filter pole is placed with the others:
 * with the plant made monic by its gain (a₂′, a₁′, a₀′ = a₂/b₀, a₁/b₀, a₀/b₀), a PID closes to
 * a₂′Tf s⁴ + (a₂′ + a₁′Tf) s³ + (a₁′ + a₀′Tf + K_p Tf + K_d) s² + (a₀′ + K_p + Kᵢ Tf) s + Kᵢ.
 * Matching a₂′Tf (s² + 2ζω s + ω²)(s + ωᵢ)(s + ω_f) = a₂′Tf (s⁴ + d₃ s³ + d₂ s² + d₁ s + d₀),
 * with the filter pole at ω_f = 10ω, gives
 * @f[
 *   \frac{1}{T_f} = d_3 - \frac{a_1'}{a_2'},\quad
 *   K_i = a_2' T_f d_0,\quad
 *   K_p = a_2' T_f d_1 - a_0' - K_i T_f,\quad
 *   K_d = a_2' T_f d_2 - a_1' - (a_0' + K_p) T_f ,
 * @f]
 * and ωᵢ = 0 is the PD that places the pair and the filter pole alone. A rigid inertia 1/(J s²) is
 * `{.num = {1}, .den = {J, 0, 0}}`; a common load-holding choice is ωᵢ = ω/8. The derivative is on
 * the measurement (c = 0, no D kick on a reference step).
 *
 * Refused with invalid_parameter if b₀ = 0, a₂/b₀ ≤ 0, ω ≤ 0, ζ ≤ 0, ωᵢ < 0 or a number is not
 * finite; not_placeable if the placement needs K_p ≤ 0, K_d < 0 or Tf ≤ 0 (the plant is already
 * faster or more damped).
 *
 * @see "Feedback Control of Dynamic Systems" (Franklin, Powell & Emami-Naeini), PD on 1/s²
 */
template<typename T>
[[nodiscard]] constexpr PIDTuning<T> pole_placement_pid(const TransferFunction<1, 3, T>& plant, T omega, T zeta, T omega_i) {
    const T b0 = plant.num[0];
    if (!valid_static_gain(b0) || !damp::finite_positive(omega) || !damp::finite_positive(zeta) || !damp::finite_non_negative(omega_i)) {
        return pid_invalid<T>();
    }
    const T a2 = plant.den[0] / b0;
    const T a1 = plant.den[1] / b0;
    const T a0 = plant.den[2] / b0;
    if (!damp::finite_positive(a2) || !damp::isfinite(a1) || !damp::isfinite(a0)) {
        return pid_invalid<T>();
    }
    const T w = omega;
    const T wi = omega_i;
    const T wf = T{10} * w;
    const T two_zw = T{2} * zeta * w;
    // (s² + 2ζω s + ω²)(s + ωᵢ) = s³ + c₁ s² + c₂ s + c₃, then × (s + ω_f).
    const T c1 = two_zw + wi;
    const T c2 = (w * w) + (two_zw * wi);
    const T c3 = w * w * wi;
    const T inv_Tf = (c1 + wf) - (a1 / a2);
    if (!damp::finite_positive(inv_Tf)) {
        return refused<PIDTuning<T>>(DesignStatus::not_placeable);
    }
    const T Tf = T{1} / inv_Tf;
    const T Ki = a2 * Tf * (c3 * wf);
    const T Kp = (a2 * Tf * (c3 + (c2 * wf))) - a0 - (Ki * Tf);
    const T Kd = (a2 * Tf * (c2 + (c1 * wf))) - a1 - ((a0 + Kp) * Tf);
    if (!damp::finite_positive(Kp) || !damp::finite_non_negative(Kd)) {
        return refused<PIDTuning<T>>(DesignStatus::not_placeable);
    }
    PIDGains<T> r{};
    r.Kp = Kp;
    r.Ki = Ki;
    r.Kd = Kd;
    r.c = T{0};                     // PI-D / PD on measurement
    r.Tf = (Kd > T{0}) ? Tf : T{0}; // Kd = 0: the filter factor cancels from the loop
    return tuned(r);
}

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

/**
 * @brief The PI (two poles) or PID (three) that places s-plane poles on a first-order plant sampled at Ts
 *
 * The poles are s-plane poles, mapped to z = e^{s·Ts} as place(sys, poles, Ts) maps them. The
 * plant b₀/(a₁ s + a₀) is ZOH-discretized to G(z) = b/(z − a), a = e^{−x}, x = a₀Ts/a₁,
 * b = (b₀Ts/a₁)(1 − e^{−x})/x (evaluated without cancellation at small x, and b = b₀Ts/a₁ for an
 * integrator). The gains are matched against the runtime's own sampled law, the backward-Euler
 * integral and (Tf = 0) the backward-difference derivative, so the placement is exact at that Ts:
 * - PI: C(z) = K_p + Kᵢ T_s z/(z − 1). (z − a)(z − 1) + b (K_p (z − 1) + Kᵢ T_s z) = z² − σ z + π
 *   gives K_p = (a − π)/b and Kᵢ from 1 + a − σ.
 * - PID: C(z) = K_p + Kᵢ T_s z/(z − 1) + (K_d/T_s)(z − 1)/z, closed loop of third order, matched
 *   to z³ − σ z² + ρ z − π.
 *
 * σ, ρ, π are the elementary symmetric functions of the z poles, real for a conjugate-closed set.
 * A pole product above the sampled plant pole (z₁z₂ > a for the PI) needs K_p < 0, and z₁z₂ = a
 * gives K_p = 0, an I-only law. The PID's K_d = −z₁z₂z₃·T_s/b is negative for every s-plane set on a
 * plant with b > 0.
 *
 * Refused with invalid_parameter if b₀ = 0, a₁ = 0, a coefficient or pole is not finite, a pole has
 * Re(s) ≥ 0 (so a z-plane value in [0, 1] read as s is refused), or a complex pole has no conjugate;
 * singular if the sampled gain b is zero. Ts is the caller's to check.
 */
template<size_t NP, typename T>
    requires(NP == 2 || NP == 3)
[[nodiscard]] constexpr PIDTuning<T> sampled_placement(
    const TransferFunction<1, 2, T>&         plant,
    const damp::array<damp::complex<T>, NP>& poles,
    T                                        Ts
) {
    using Cplx = damp::complex<T>;
    const T b0 = plant.num[0];
    const T a1 = plant.den[0];
    const T a0 = plant.den[1];
    if (!valid_static_gain(b0) || !valid_static_gain(a1) || !damp::isfinite(a0)) {
        return pid_invalid<T>();
    }
    damp::array<Cplx, NP> z{};
    T                     zmag = T{0};
    for (size_t i = 0; i < NP; ++i) {
        // Re(s) < 0: the loop is stable, and a z-plane value in [0, 1] read as s is refused.
        if (!damp::isfinite(poles[i].real()) || !damp::isfinite(poles[i].imag()) || !(poles[i].real() < T{0})) {
            return pid_invalid<T>();
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
        return pid_invalid<T>(); // a complex pole without its conjugate
    }
    const T x = a0 * Ts / a1;
    const T a = damp::exp(-x);
    const T b = (b0 * Ts / a1) * one_minus_exp_neg_over_x(x);
    if (!damp::isfinite(b) || b == T{0}) {
        return refused<PIDTuning<T>>(DesignStatus::singular);
    }
    if constexpr (NP == 2) {
        const T Kp = (a - prod.real()) / b;                             // a − b Kp = π
        const T Ki = ((T{1} + a) - sigma.real() - (b * Kp)) / (b * Ts); // b Kp + b Ki Ts = 1 + a − σ
        return tuned(PIDGains<T>{.Kp = Kp, .Ki = Ki});
    } else {
        const T Kd = -(prod.real() * Ts) / b;                                           // b Kd/Ts = −π
        const T Kp = (a - rho.real() - (T{2} * b * Kd / Ts)) / b;                       // b Kp + 2 b Kd/Ts = a − ρ
        const T Ki = ((T{1} + a - sigma.real()) - (b * Kp) - (b * Kd / Ts)) / (b * Ts); // b Kp + b Ki Ts + b Kd/Ts = 1 + a − σ
        return tuned(PIDGains<T>{.Kp = Kp, .Ki = Ki, .Kd = Kd});
    }
}

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
    const StateSpace plant{.A = sys.A, .B = sys.B, .C = sys.C, .D = sys.D, .Ts = sys.Ts};
    if (!damp::isfinite(plant.Ts) || plant.Ts != T{0}) {
        return {.status = DesignStatus::invalid_sample_time};
    }
    if (plant.D(0, 0) != T{0}) {
        return {.status = DesignStatus::invalid_parameter};
    }
    const auto tf = ::damp::ss2tf(plant);
    const auto scale = ::damp::detail::leverrier_num_scale(plant);
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

} // namespace damp::design::detail

#pragma pop_macro("PI")
