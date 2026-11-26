// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
// IWYU pragma: private, include "damp/controllers/adrc.hpp"

/**
 * @file adrc_gains.hpp
 * @brief What ADRCController ticks on: its bandwidths, b₀ and PD gains at the fixed sample period,
 *        and where its ESO's maps live, as members or as the constants of a `constexpr` design's Law
 */

#include "damp/backend.hpp"
#include "damp/estimation/eso.hpp"
#include "damp/math/math.hpp"

namespace damp::detail {

/**
 * @brief An ADRCController's gains, its fixed sample period and its PD gains: the members an
 *        ADRCController holds, its ESO's maps members too (HeldEso), or the constants of a
 *        `constexpr` design's Law (ConstantADRCGains)
 */
template<size_t NX, typename T>
struct ADRCGainSet {
    static constexpr bool constant = false; ///< Members of the runtime, settable
    using Eso = HeldEso<NX, T>;             ///< The ESO's maps: members

    T wc{T{0}}; ///< [rad/s] Controller bandwidth ω_c
    T wo{T{0}}; ///< [rad/s] Observer bandwidth ω_o
    T b0{T{0}}; ///< Plant gain b₀, which the tick divides by
    T Kp{T{0}}; ///< Proportional gain (ω_c or ω_c²)
    T Kd{T{0}}; ///< Derivative gain (0 or 2ω_c)
    T Ts{T{0}}; ///< [s] The tick's period, which set_gains designs at
};

/**
 * @brief A `constexpr` design::adrc's gains as constants of the type: what an ADRCController reads
 *        in place of its gain members, in the runtime's scalar, beside its ESO's maps on the same
 *        design (ConstantEso)
 *
 * Each member is static, so the Law stores none of them: K_p, K_d and b₀ sit in flash and fold
 * into the tick, a zero K_d (a first-order law) is no code, and the ESO's chain of integrators
 * costs nothing.
 *
 * @tparam Res design::adrc's result, `constexpr`
 * @tparam T   The runtime's scalar
 */
template<const auto& Res, typename T>
struct ConstantADRCGains {
    static constexpr bool constant = true; ///< The gains are the type's: set_gains is absent
    using Eso = ConstantEso<Res, T>;       ///< The ESO's maps: constants of the type

    static constexpr auto design = Res.template as<T>(); ///< The design, in the runtime's scalar
    static constexpr T    wc = design.gains.wc;          ///< [rad/s] Controller bandwidth ω_c
    static constexpr T    wo = design.gains.wo;          ///< [rad/s] Observer bandwidth ω_o
    static constexpr T    b0 = design.gains.b0;          ///< Plant gain b₀
    static constexpr T    Kp = design.Kp;                ///< Proportional gain
    static constexpr T    Kd = design.Kd;                ///< Derivative gain
    static constexpr T    Ts = design.Ts;                ///< [s] The tick's period

    static constexpr bool ok = damp::finite_positive(b0) && damp::isfinite(Kp) && damp::isfinite(Kd); ///< The tick can divide by b₀ in the runtime's scalar
};

} // namespace damp::detail
