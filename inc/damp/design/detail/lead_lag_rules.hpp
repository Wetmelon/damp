// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
// IWYU pragma: private, include "damp/controllers/lead_lag.hpp"

/**
 * @file lead_lag_rules.hpp
 * @brief Where design::lead and design::lag place a section's zero, pole and gain
 */

#include "damp/backend.hpp"
#include "damp/math/math.hpp"

namespace damp::design::detail {

/** @brief A section K(s + z)/(s + p), or none (ok false) when the spec is out of its domain. */
template<typename T>
struct ZeroPoleGain {
    T    K{};       ///< High-frequency gain
    T    z{};       ///< [rad/s] Zero at s = −z
    T    p{};       ///< [rad/s] Pole at s = −p
    bool ok{false}; ///< The spec was in its domain
};

/**
 * @brief The lead section with its maximum phase advance @p phi_max [rad] at @p wc [rad/s] and unit gain there
 *
 * Places zero and pole symmetrically (in log-frequency) around ωc, so the maximum phase advance
 * is exactly at ωc:
 *
 *     α = (1 − sin φ_max) / (1 + sin φ_max),   z = ωc √α,   p = ωc / √α
 *
 * with α = z/p < 1, and K = 1/√α so |C(jωc)| = K·√α = 1: phase boost without shifting the
 * gain-crossover frequency.
 *
 * @see "Feedback Control of Dynamic Systems" (Franklin, Powell & Emami-Naeini), lead compensation
 * @return The section; none for φ_max outside (0, π/2) or an ωc that is not finite and positive
 */
template<typename T>
[[nodiscard]] constexpr ZeroPoleGain<T> lead_section(T phi_max, T wc) {
    if (!damp::finite_positive(phi_max) || !(phi_max < damp::numbers::pi_v<T> / T{2}) || !damp::finite_positive(wc)) {
        return {};
    }
    const T s = damp::sin(phi_max);
    const T alpha = (T{1} - s) / (T{1} + s);
    const T sqrt_alpha = damp::sqrt(alpha);
    return ZeroPoleGain<T>{.K = T{1} / sqrt_alpha, .z = wc * sqrt_alpha, .p = wc / sqrt_alpha, .ok = true};
}

/**
 * @brief The lag section that raises the DC gain by @p dc_gain_boost, a factor @p margin_factor below @p wc
 *
 *     z = ωc / margin_factor,   p = z / dc_gain_boost,   K = 1
 *
 * so the DC gain K·z/p is dc_gain_boost (z > p), and the phase lag at ωc stays small: about
 * atan(1/margin_factor) − atan(1/(margin_factor·dc_gain_boost)).
 *
 * @see "Feedback Control of Dynamic Systems" (Franklin, Powell & Emami-Naeini), lag compensation
 * @return The section; none for a boost not above 1, or an ωc or margin that is not finite and positive
 */
template<typename T>
[[nodiscard]] constexpr ZeroPoleGain<T> lag_section(T dc_gain_boost, T wc, T margin_factor) {
    if (!damp::isfinite(dc_gain_boost) || !(dc_gain_boost > T{1}) || !damp::finite_positive(wc) || !damp::finite_positive(margin_factor)) {
        return {};
    }
    const T z = wc / margin_factor;
    return ZeroPoleGain<T>{.K = T{1}, .z = z, .p = z / dc_gain_boost, .ok = true};
}

} // namespace damp::design::detail
