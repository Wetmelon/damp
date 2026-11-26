// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
// IWYU pragma: private, include "damp/controllers/lead_lag.hpp"

/**
 * @file lead_lag_stage.hpp
 * @brief One first-order section of LeadLagController, K(s + z)/(s + p) by Tustin, its memory and its
 *        coefficients, as members or as the constants of a `constexpr` design's Law
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
 * @brief The memory of a first-order section of LeadLagController, run on its coefficients:
 *        u[n] = b₀e[n] + b₁e[n−1] − a₁u[n−1]
 *
 * The coefficients are the law's (a LeadLagController's members, or the constants of a
 * `constexpr` design's Law), so a retune keeps the memory. The last section of a series runs on
 * the output that went out (its u_prev is clamped, or set by a clamp reported after it); the
 * others on their own outputs.
 */
template<typename T>
struct LeadLagStage {
    T e_prev{T{0}}; ///< The previous input e[n−1]
    T u_prev{T{0}}; ///< The previous output u[n−1]

    constexpr LeadLagStage() = default;

    /** @brief Convert across scalar precision, preserving the memory. */
    template<typename U>
    constexpr explicit LeadLagStage(const LeadLagStage<U>& o) : e_prev(static_cast<T>(o.e_prev)),
                                                                u_prev(static_cast<T>(o.u_prev)) {}

    /** @brief The output for the input @p e on the coefficients @p k, not committed. */
    [[nodiscard]] constexpr T raw(T e, const LeadLagCoefficients<T>& k) const { return (k.b0 * e) + (k.b1 * e_prev) - (k.a1 * u_prev); }

    /** @brief One tick on @p e and the coefficients @p k, its own output fed back. */
    constexpr T operator()(T e, const LeadLagCoefficients<T>& k) {
        const T u = raw(e, k);
        e_prev = e;
        u_prev = u;
        return u;
    }

    /**
     * @brief The memory after a tick on @p e that makes the next tick's output @p target when its
     *        input is @p e_next, on the coefficients @p k
     *
     * u[n−1] = (b₀e_next + b₁e − target)/a₁. A section whose pole is at z = 0 within rounding
     * (|a₁| ≤ 4ε, p = 2/T_s) has no memory of its output to preload: u[n−1] = target, and the next
     * output is its own.
     */
    constexpr void preload(T e, T e_next, T target, const LeadLagCoefficients<T>& k) {
        e_prev = e;
        u_prev = (damp::abs(k.a1) > (T{4} * eps<T>())) ? (((k.b0 * e_next) + (k.b1 * e) - target) / k.a1) : target;
    }

    /** @brief Clear the memory. */
    constexpr void reset() {
        e_prev = T{0};
        u_prev = T{0};
    }
};

/** @brief One section of a lead-lag law, K(s + z)/(s + p), and its Tustin coefficients at the law's T_s. */
template<typename T>
struct LeadLagSection {
    T                      K{T{0}}; ///< High-frequency gain
    T                      z{T{0}}; ///< [rad/s] Zero location
    T                      p{T{0}}; ///< [rad/s] Pole location
    LeadLagCoefficients<T> k{};     ///< b₀, b₁, a₁ at T_s
};

/**
 * @brief A LeadLagController's sections, its fixed sample period and the coefficients at it: the
 *        members a LeadLagController holds, or the constants of a `constexpr` design's Law
 *        (ConstantLeadLagGains)
 */
template<typename T, size_t NS>
struct LeadLagGainSet {
    static constexpr bool constant = false; ///< Members of the runtime, settable

    damp::array<LeadLagSection<T>, NS> sections{}; ///< Each section, in series from the input
    T                                  Ts{T{0}};   ///< [s] The tick's period
};

/**
 * @brief The gain set of the sections @p g (each with a K, z and p) at the period @p Ts, or
 *        nullopt when Ts is not finite and positive, a K is not finite, a z or p is not finite and
 *        positive, or a coefficient overflows
 */
template<typename T, size_t NS, typename Section>
[[nodiscard]] constexpr damp::optional<LeadLagGainSet<T, NS>> lead_lag_gain_set(const damp::array<Section, NS>& g, T Ts) {
    if (!damp::finite_positive(Ts)) {
        return damp::nullopt;
    }
    LeadLagGainSet<T, NS> set{.Ts = Ts};
    for (size_t i = 0; i < NS; ++i) {
        if (!lead_lag_section_ok(g[i].K, g[i].z, g[i].p)) {
            return damp::nullopt;
        }
        set.sections[i] = LeadLagSection<T>{.K = g[i].K, .z = g[i].z, .p = g[i].p, .k = lead_lag_coefficients(g[i].K, g[i].z, g[i].p, Ts)};
        if (!set.sections[i].k.finite()) {
            return damp::nullopt;
        }
    }
    return set;
}

/**
 * @brief A `constexpr` lead-lag design's sections as constants of the type: what a
 *        LeadLagController reads in place of its gain members, in the runtime's scalar
 *
 * Each member is static, so the Law stores none of them: each section's b₀, b₁ and a₁ sit in
 * flash and fold into the tick. A result whose sections are refused has no set, which the Law's
 * static_assert reports.
 *
 * @tparam Res A lead-lag design's result, `constexpr`
 * @tparam T   The runtime's scalar
 */
template<const auto& Res, typename T>
struct ConstantLeadLagGains {
    static constexpr bool   constant = true;                                     ///< The gains are the type's: set_gains is absent
    static constexpr size_t kNS = damp::remove_cvref_t<decltype(Res)>::sections; ///< The sections in series

    static constexpr auto design = Res.template as<T>();                    ///< The design, in the runtime's scalar
    static constexpr auto opt = lead_lag_gain_set(design.gains, design.Ts); ///< The sections at T_s, if they are in their domain
    static constexpr bool ok = opt.has_value();                             ///< The sections are in their domain at T_s
    static constexpr auto set = ok ? *opt : LeadLagGainSet<T, kNS>{};       ///< The sections at T_s

    static constexpr auto sections = set.sections; ///< Each section and its coefficients
    static constexpr T    Ts = set.Ts;             ///< [s] The tick's period
};

} // namespace damp::detail
