// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
// IWYU pragma: private, include "damp/controllers/esc.hpp"

/**
 * @file esc_gains.hpp
 * @brief What ExtremumSeekingController ticks on: its dither, adaptation gain and filter corners
 *        at the fixed sample period, with the rotation and the filter coefficients computed once
 *        when they are set, as members or as the constants of a `constexpr` design's Law
 */

#include "damp/backend.hpp"
#include "damp/math/math.hpp"

namespace damp::detail {

/** @brief First-order discrete low-pass coefficient ω_cT_s/(1 + ω_cT_s) for the corner @p wc [rad/s] at @p Ts; 1 (pass-through) for ω_c = 0. */
template<typename T>
[[nodiscard]] constexpr T esc_lpf_alpha(T wc, T Ts) {
    if (!(wc > T{0})) {
        return T{1};
    }
    return (wc * Ts) / (T{1} + (wc * Ts));
}

/**
 * @brief An ExtremumSeekingController's gains, its fixed sample period and what it ticks on at
 *        it: the members an ExtremumSeekingController holds, or the constants of a `constexpr`
 *        design's Law (ConstantESCGains)
 */
template<typename T>
struct ESCGainSet {
    static constexpr bool constant = false; ///< Members of the runtime, settable

    T dither_amplitude{T{0}}; ///< a, the dither's amplitude (input units)
    T dither_omega{T{0}};     ///< [rad/s] ω, the dither's frequency
    T gain{T{0}};             ///< k, the adaptation (integrator) gain
    T hp_omega{T{0}};         ///< [rad/s] The high-pass corner in force
    T lp_omega{T{0}};         ///< [rad/s] The gradient low-pass corner (0: none)
    T direction{T{1}};        ///< +1 maximize, −1 minimize
    T demod_phase{T{0}};      ///< [rad] φ, the demodulation phase offset
    T u_init{T{0}};           ///< The initial operating point û₀, which reset() returns to
    T Ts{T{0}};               ///< [s] The tick's period
    T hp_alpha{T{0}};         ///< The high-pass's coefficient at T_s
    T lp_alpha{T{1}};         ///< The low-pass's coefficient at T_s (1: none)
    T rot_c{T{1}};            ///< cos(ωT_s): the dither phasor's rotation per tick
    T rot_s{T{0}};            ///< sin(ωT_s)
    T demod_c{T{1}};          ///< cos φ
    T demod_s{T{0}};          ///< sin φ
    T k_ts{T{0}};             ///< ±k·T_s: the operating point's step per unit gradient, signed by the direction
};

/**
 * @brief The gain set of the gains @p g (dither_amplitude, dither_omega, gain, hp_omega, lp_omega)
 *        with the direction, demodulation phase and initial point of a design, at the period @p Ts
 *
 * A high-pass corner of 0 is ω/5, as design::esc takes it. Trigonometry runs here, once: the
 * tick rotates the dither's phasor by (cos ωT_s, sin ωT_s) and demodulates with (cos φ, sin φ).
 *
 * @return The set, or nullopt when a, ω or k is not finite and positive, a corner is below 0 or
 *         not finite, Ts is not finite and positive, the dither has fewer than 8 samples per
 *         period (ω·T_s > π/4), the direction, φ or û₀ is not finite, or a coefficient overflows
 */
template<typename T, typename Gains>
[[nodiscard]] constexpr damp::optional<ESCGainSet<T>> esc_gain_set(const Gains& g, T direction, T demod_phase, T u_init, T Ts) {
    const bool in_domain = damp::finite_positive(g.dither_amplitude) && damp::finite_positive(g.dither_omega) && damp::finite_positive(g.gain)
                        && damp::finite_non_negative(g.hp_omega) && damp::finite_non_negative(g.lp_omega) && damp::finite_positive(Ts)
                        && damp::isfinite(direction) && damp::isfinite(demod_phase) && damp::isfinite(u_init);
    if (!in_domain || !(g.dither_omega * Ts <= damp::numbers::pi_v<T> / T{4})) {
        return damp::nullopt;
    }
    ESCGainSet<T> set{
        .dither_amplitude = g.dither_amplitude,
        .dither_omega = g.dither_omega,
        .gain = g.gain,
        .hp_omega = (g.hp_omega > T{0}) ? g.hp_omega : (g.dither_omega / T{5}),
        .lp_omega = g.lp_omega,
        .direction = direction,
        .demod_phase = demod_phase,
        .u_init = u_init,
        .Ts = Ts,
    };
    set.hp_alpha = esc_lpf_alpha(set.hp_omega, Ts);
    set.lp_alpha = esc_lpf_alpha(set.lp_omega, Ts);
    set.rot_c = damp::cos(g.dither_omega * Ts);
    set.rot_s = damp::sin(g.dither_omega * Ts);
    set.demod_c = damp::cos(demod_phase);
    set.demod_s = damp::sin(demod_phase);
    set.k_ts = direction * g.gain * Ts;
    if (!damp::isfinite(set.hp_alpha) || !damp::isfinite(set.lp_alpha) || !damp::isfinite(set.k_ts)) {
        return damp::nullopt;
    }
    return set;
}

/**
 * @brief A `constexpr` extremum-seeking design's gains as constants of the type: what an
 *        ExtremumSeekingController reads in place of its gain members, in the runtime's scalar
 *
 * Each member is static, so the Law stores none of them: the dither's rotation, the filters'
 * coefficients and ±k·T_s sit in flash and fold into the tick. A result whose gains are refused
 * has no set, which the Law's static_assert reports.
 *
 * @tparam Res design::esc's or design::esc_mppt's result, `constexpr`
 * @tparam T   The runtime's scalar
 */
template<const auto& Res, typename T>
struct ConstantESCGains {
    static constexpr bool constant = true; ///< The gains are the type's: set_gains is absent

    static constexpr auto design = Res.template as<T>();                                                                    ///< The design, in the runtime's scalar
    static constexpr auto opt = esc_gain_set(design.gains, design.direction, design.demod_phase, design.u_init, design.Ts); ///< The gains at T_s, if they are in their domain
    static constexpr bool ok = opt.has_value();                                                                             ///< The gains are in their domain at T_s
    static constexpr auto set = ok ? *opt : ESCGainSet<T>{};                                                                ///< The gains at T_s

    static constexpr T dither_amplitude = set.dither_amplitude; ///< a
    static constexpr T dither_omega = set.dither_omega;         ///< [rad/s] ω
    static constexpr T gain = set.gain;                         ///< k
    static constexpr T hp_omega = set.hp_omega;                 ///< [rad/s] The high-pass corner
    static constexpr T lp_omega = set.lp_omega;                 ///< [rad/s] The low-pass corner (0: none)
    static constexpr T direction = set.direction;               ///< +1 maximize, −1 minimize
    static constexpr T demod_phase = set.demod_phase;           ///< [rad] φ
    static constexpr T u_init = set.u_init;                     ///< û₀
    static constexpr T Ts = set.Ts;                             ///< [s] The tick's period
    static constexpr T hp_alpha = set.hp_alpha;                 ///< The high-pass's coefficient
    static constexpr T lp_alpha = set.lp_alpha;                 ///< The low-pass's coefficient
    static constexpr T rot_c = set.rot_c;                       ///< cos(ωT_s)
    static constexpr T rot_s = set.rot_s;                       ///< sin(ωT_s)
    static constexpr T demod_c = set.demod_c;                   ///< cos φ
    static constexpr T demod_s = set.demod_s;                   ///< sin φ
    static constexpr T k_ts = set.k_ts;                         ///< ±k·T_s
};

} // namespace damp::detail
