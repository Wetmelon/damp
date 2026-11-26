// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
// IWYU pragma: private, include "damp/controllers/lead_lag.hpp"

/**
 * @file lead_lag_stage.hpp
 * @brief One first-order section of LeadLagController, K(s + z)/(s + p) by Tustin, and its coefficients
 */

#include "damp/backend.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix_traits.hpp" // eps

namespace damp::detail {

/** @brief The coefficients of a first-order section u[n] = b₀e[n] + b₁e[n−1] − a₁u[n−1]. */
template<typename T>
struct LeadLagCoefficients {
    T b0{T{0}}; ///< Gain on e[n]
    T b1{T{0}}; ///< Gain on e[n−1]
    T a1{T{0}}; ///< Feedback on u[n−1]

    /** @brief Every coefficient finite. */
    [[nodiscard]] constexpr bool finite() const { return damp::isfinite(b0) && damp::isfinite(b1) && damp::isfinite(a1); }
};

/**
 * @brief Tustin coefficients of C(s) = K(s + z)/(s + p) at the period @p Ts
 *
 * With s = k(1 − z⁻¹)/(1 + z⁻¹), k = 2/T_s, the numerator is K((k + z) + (z − k)z⁻¹) and the
 * denominator (k + p) + (p − k)z⁻¹, both over 1 + z⁻¹, normalized by k + p:
 * b₀ = K(k + z)/(k + p), b₁ = K(z − k)/(k + p), a₁ = (p − k)/(k + p).
 * The caller checks @p z, @p p and @p Ts finite and positive first; for those every coefficient
 * is finite unless K(k + z) overflows.
 */
template<typename T>
[[nodiscard]] constexpr LeadLagCoefficients<T> lead_lag_coefficients(T K, T z, T p, T Ts) {
    const T k = T{2} / Ts; // Tustin substitution factor
    const T a0_raw = k + p;
    const T a1_raw = p - k;
    const T b0_raw = K * (k + z);
    const T b1_raw = K * (z - k);
    return LeadLagCoefficients<T>{.b0 = b0_raw / a0_raw, .b1 = b1_raw / a0_raw, .a1 = a1_raw / a0_raw};
}

/** @brief Whether a section K(s + z)/(s + p) is in its domain: K finite, z and p finite and positive. */
template<typename T>
[[nodiscard]] constexpr bool lead_lag_section_ok(T K, T z, T p) {
    return damp::isfinite(K) && damp::finite_positive(z) && damp::finite_positive(p);
}

/**
 * @brief A first-order section of LeadLagController and its memory: u[n] = b₀e[n] + b₁e[n−1] − a₁u[n−1]
 *
 * The last section of a series runs on the output that went out (its u_prev is clamped, or set
 * by a clamp reported after it); the others on their own outputs.
 */
template<typename T>
struct LeadLagStage {
    LeadLagCoefficients<T> k{};          ///< The coefficients at the law's T_s
    T                      e_prev{T{0}}; ///< The previous input e[n−1]
    T                      u_prev{T{0}}; ///< The previous output u[n−1]

    constexpr LeadLagStage() = default;

    /** @brief Convert across scalar precision, preserving coefficients and memory. */
    template<typename U>
    constexpr explicit LeadLagStage(const LeadLagStage<U>& o)
        : k{static_cast<T>(o.k.b0), static_cast<T>(o.k.b1), static_cast<T>(o.k.a1)},
          e_prev(static_cast<T>(o.e_prev)),
          u_prev(static_cast<T>(o.u_prev)) {}

    /** @brief The output for the input @p e, not committed. */
    [[nodiscard]] constexpr T raw(T e) const { return (k.b0 * e) + (k.b1 * e_prev) - (k.a1 * u_prev); }

    /** @brief One tick on @p e, its own output fed back. */
    constexpr T operator()(T e) {
        const T u = raw(e);
        e_prev = e;
        u_prev = u;
        return u;
    }

    /**
     * @brief The memory after a tick on @p e that makes the next tick's output @p target when its input is @p e_next
     *
     * u[n−1] = (b₀e_next + b₁e − target)/a₁. A section whose pole is at z = 0 within rounding
     * (|a₁| ≤ 4ε, p = 2/T_s) has no memory of its output to preload: u[n−1] = target, and the next
     * output is its own.
     */
    constexpr void preload(T e, T e_next, T target) {
        e_prev = e;
        u_prev = (damp::abs(k.a1) > (T{4} * eps<T>())) ? (((k.b0 * e_next) + (k.b1 * e) - target) / k.a1) : target;
    }

    /** @brief Clear the memory; the coefficients are kept. */
    constexpr void reset() {
        e_prev = T{0};
        u_prev = T{0};
    }
};

} // namespace damp::detail
