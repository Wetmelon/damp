// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
// IWYU pragma: private

/**
 * @file sliding_mode.hpp
 * @brief What the sliding-mode runtimes (SMCController, STSMCController) tick on: their gains at
 *        the fixed sample period, with the coefficients computed once when the gains are set, as
 *        members or as the constants of a `constexpr` design's Law
 */

#include "damp/backend.hpp"
#include "damp/math/math.hpp"

namespace damp::detail {

/**
 * @brief An SMCController's gains, its fixed sample period and the coefficients at it: the
 *        members an SMCController holds, or the constants of a `constexpr` design's Law
 */
template<typename T>
struct SMCGainSet {
    static constexpr bool constant = false; ///< Members of the runtime, settable

    T lambda{T{0}};  ///< [1/s] Sliding-surface slope λ
    T k{T{0}};       ///< Switching gain
    T b0{T{0}};      ///< Input gain b₀
    T phi{T{0}};     ///< Boundary-layer thickness φ (0: the hard relay)
    T Ts{T{0}};      ///< [s] The tick's period
    T relay{T{0}};   ///< k/b₀: the relay's magnitude
    T inv_ts{T{0}};  ///< 1/T_s: the backward difference's scale
    T inv_phi{T{0}}; ///< 1/φ, 0 for φ = 0 (the hard relay)
};

/**
 * @brief The gain set of λ, k, b₀, φ at the period @p Ts, or nullopt when a gain is outside its
 *        domain (λ, k, b₀ finite and positive, φ finite and ≥ 0), Ts is not finite and positive,
 *        or a coefficient (k/b₀, 1/T_s, 1/φ) overflows
 */
template<typename T>
[[nodiscard]] constexpr damp::optional<SMCGainSet<T>> smc_gain_set(T lambda, T k, T b0, T phi, T Ts) {
    const bool in_domain = damp::finite_positive(lambda) && damp::finite_positive(k) && damp::finite_positive(b0)
                        && damp::finite_non_negative(phi) && damp::finite_positive(Ts);
    if (!in_domain) {
        return damp::nullopt;
    }
    SMCGainSet<T> set{.lambda = lambda, .k = k, .b0 = b0, .phi = phi, .Ts = Ts, .relay = k / b0, .inv_ts = T{1} / Ts};
    if (phi > T{0}) {
        set.inv_phi = T{1} / phi;
    }
    if (!damp::isfinite(set.relay) || !damp::isfinite(set.inv_ts) || !damp::isfinite(set.inv_phi)) {
        return damp::nullopt;
    }
    return set;
}

/**
 * @brief A `constexpr` SMC design's gains as constants of the type: what an SMCController reads
 *        in place of its gain members, in the runtime's scalar
 *
 * Each member is static, so the Law stores none of them: they sit in flash and fold into the
 * tick. A result whose gains are refused has no set, which the Law's static_assert reports.
 */
template<const auto& Res, typename T>
struct ConstantSMCGains {
    static constexpr bool constant = true; ///< The gains are the type's: set_gains is absent

    static constexpr auto design = Res.template as<T>();                                                                         ///< The design, in the runtime's scalar
    static constexpr auto opt = smc_gain_set(design.gains.lambda, design.gains.k, design.gains.b0, design.gains.phi, design.Ts); ///< The gains at the design's T_s, if they are in their domain
    static constexpr bool ok = opt.has_value();                                                                                  ///< The gains are in their domain at T_s
    static constexpr auto set = ok ? *opt : SMCGainSet<T>{};                                                                     ///< The gains at the design's T_s

    static constexpr T lambda = set.lambda;   ///< [1/s] Sliding-surface slope λ
    static constexpr T k = set.k;             ///< Switching gain
    static constexpr T b0 = set.b0;           ///< Input gain b₀
    static constexpr T phi = set.phi;         ///< Boundary-layer thickness φ
    static constexpr T Ts = set.Ts;           ///< [s] The tick's period
    static constexpr T relay = set.relay;     ///< k/b₀
    static constexpr T inv_ts = set.inv_ts;   ///< 1/T_s
    static constexpr T inv_phi = set.inv_phi; ///< 1/φ (0: the hard relay)
};

/**
 * @brief An STSMCController's gains, its fixed sample period and the coefficients at it: the
 *        members an STSMCController holds, or the constants of a `constexpr` design's Law
 */
template<typename T>
struct STSMCGainSet {
    static constexpr bool constant = false; ///< Members of the runtime, settable

    T k1{T{0}};     ///< Gain on |s|^½·sign(s)
    T k2{T{0}};     ///< Integral gain on sign(s)
    T lambda{T{0}}; ///< [1/s] Slope of the (r, y) form's surface
    T k_lin{T{0}};  ///< Generalized-STA linear gain (0: the classic algorithm)
    T phi{T{0}};    ///< Half-width φ of the smoothed sign (0: the true sign)
    T Ts{T{0}};     ///< [s] The tick's period
    T inv_ts{T{0}}; ///< 1/T_s: the (r, y) form's backward difference
};

/**
 * @brief The gain set of k₁, k₂, λ, k_lin, φ at the period @p Ts, or nullopt when a gain is
 *        outside its domain (k₁, k₂ finite and positive; λ, k_lin, φ finite and ≥ 0), Ts is not
 *        finite and positive, or 1/T_s overflows
 */
template<typename T>
[[nodiscard]] constexpr damp::optional<STSMCGainSet<T>> stsmc_gain_set(T k1, T k2, T lambda, T k_lin, T phi, T Ts) {
    const bool in_domain = damp::finite_positive(k1) && damp::finite_positive(k2) && damp::finite_non_negative(lambda)
                        && damp::finite_non_negative(k_lin) && damp::finite_non_negative(phi) && damp::finite_positive(Ts);
    if (!in_domain) {
        return damp::nullopt;
    }
    const STSMCGainSet<T> set{.k1 = k1, .k2 = k2, .lambda = lambda, .k_lin = k_lin, .phi = phi, .Ts = Ts, .inv_ts = T{1} / Ts};
    if (!damp::isfinite(set.inv_ts)) {
        return damp::nullopt;
    }
    return set;
}

/**
 * @brief A `constexpr` super-twisting design's gains as constants of the type: what an
 *        STSMCController reads in place of its gain members, in the runtime's scalar
 */
template<const auto& Res, typename T>
struct ConstantSTSMCGains {
    static constexpr bool constant = true; ///< The gains are the type's: set_gains is absent

    static constexpr auto design = Res.template as<T>();                                                                                                ///< The design, in the runtime's scalar
    static constexpr auto opt = stsmc_gain_set(design.gains.k1, design.gains.k2, design.gains.lambda, design.gains.k_lin, design.gains.phi, design.Ts); ///< The gains at the design's T_s, if they are in their domain
    static constexpr bool ok = opt.has_value();                                                                                                         ///< The gains are in their domain at T_s
    static constexpr auto set = ok ? *opt : STSMCGainSet<T>{};                                                                                          ///< The gains at the design's T_s

    static constexpr T k1 = set.k1;         ///< Gain on |s|^½·sign(s)
    static constexpr T k2 = set.k2;         ///< Integral gain on sign(s)
    static constexpr T lambda = set.lambda; ///< [1/s] Slope of the (r, y) form's surface
    static constexpr T k_lin = set.k_lin;   ///< Generalized-STA linear gain
    static constexpr T phi = set.phi;       ///< Half-width φ of the smoothed sign
    static constexpr T Ts = set.Ts;         ///< [s] The tick's period
    static constexpr T inv_ts = set.inv_ts; ///< 1/T_s
};

} // namespace damp::detail
