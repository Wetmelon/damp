// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
// IWYU pragma: private

/**
 * @file smith_model.hpp
 * @brief What the Smith predictor runs on: its first-order model sampled with the input delay
 *        (design::smith_predictor), as members or as the constants of a `constexpr` design's Law,
 *        and the ring of the commands in flight (SmithPredictor)
 */

#include "damp/backend.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/systems/discretization.hpp"

namespace damp {

namespace detail {

/**
 * @brief The slot @p j commands before the newest, @p newest, in a ring of @p last + 1 slots
 *        (j ≤ last): a compare and a wrap, no division
 */
[[nodiscard]] constexpr size_t ring_back(size_t newest, size_t j, size_t last) {
    return (newest >= j) ? (newest - j) : (newest + last + 1 - j);
}

/**
 * @brief A SmithPredictor's model: its pole, its gains and the delay's whole samples, the members
 *        a SmithPredictor holds or the constants of a `constexpr` design's Law (ConstantSmithGains)
 */
template<typename T>
struct SmithGainSet {
    static constexpr bool constant = false; ///< Members of the runtime

    T      a{T{0}};          ///< e^{−T_s/τ}
    T      b{T{0}};          ///< The delay-free model's gain on u_k
    T      b_fresh{T{0}};    ///< The delayed model's gain on u_{k−d}
    T      b_stale{T{0}};    ///< Its gain on u_{k−d−1} (0 for a delay of whole samples)
    size_t delay_samples{0}; ///< d, the delay's whole periods
};

/**
 * @brief The model of pole @p a, gains @p b, @p b_fresh, @p b_stale and @p d whole periods, or
 *        nullopt when it does not run in N delay states: an entry not finite, b = 0, d > N, or
 *        d = N with a gain on the command before the oldest held
 */
template<size_t N, typename T>
[[nodiscard]] constexpr damp::optional<SmithGainSet<T>> smith_gain_set(T a, T b, T b_fresh, T b_stale, size_t d) {
    const bool finite = damp::isfinite(a) && damp::isfinite(b) && damp::isfinite(b_fresh) && damp::isfinite(b_stale);
    if (!finite || (b == T{0}) || (d > N) || (d == N && b_stale != T{0})) {
        return damp::nullopt;
    }
    return SmithGainSet<T>{.a = a, .b = b, .b_fresh = b_fresh, .b_stale = b_stale, .delay_samples = d};
}

/**
 * @brief A `constexpr` Smith predictor design's model as constants of the type: what a
 *        SmithPredictor reads in place of its model members, in the runtime's scalar
 *
 * Each member is static, so the Law stores none of them: the pole, the gains and the delay's
 * whole samples sit in flash and fold into the tick, which reads the ring at a constant offset.
 *
 * @tparam Res design::smith_predictor's result, `constexpr`
 * @tparam T   The runtime's scalar
 * @tparam N   The runtime's delay states
 */
template<const auto& Res, typename T, size_t N>
struct ConstantSmithGains {
    static constexpr bool constant = true; ///< The model is the type's

    static constexpr auto design = Res.template as<T>();                                                                     ///< The design, in the runtime's scalar
    static constexpr auto opt = smith_gain_set<N>(design.a, design.b, design.b_fresh, design.b_stale, design.delay_samples); ///< The model, if it runs in N delay states
    static constexpr bool ok = opt.has_value();                                                                              ///< The model runs
    static constexpr auto set = ok ? *opt : SmithGainSet<T>{};                                                               ///< The model

    static constexpr T      a = set.a;                         ///< e^{−T_s/τ}
    static constexpr T      b = set.b;                         ///< The delay-free model's gain on u_k
    static constexpr T      b_fresh = set.b_fresh;             ///< The delayed model's gain on u_{k−d}
    static constexpr T      b_stale = set.b_stale;             ///< Its gain on u_{k−d−1}
    static constexpr size_t delay_samples = set.delay_samples; ///< d, the delay's whole periods
};

} // namespace detail

namespace design::detail {

using ::damp::detail::channel_delays;
using ::damp::detail::channel_seconds;
using ::damp::detail::split_delay;

/**
 * @brief The first-order model K/(τs + 1) under ZOH at T_s, with an input delay θ = d·T_s + τ_f
 *        (0 ≤ τ_f < T_s): the delay-free model x⁺ = a x + b u, and the delayed one
 *        x⁺ = a x + b_fresh u_{k−d} + b_stale u_{k−d−1}
 *
 * Over a period the delayed input holds u_{k−d−1} until τ_f, then u_{k−d}, so
 * b_stale = e^{−(T_s−τ_f)/τ}·K(1 − e^{−τ_f/τ}) and b_fresh = K(1 − e^{−(T_s−τ_f)/τ}), and
 * b = b_fresh + b_stale = K(1 − a), so both models settle on the same output. A delay of whole
 * samples has b_stale = 0.
 */
template<typename T>
struct SmithModel {
    T            a{T{0}};                  ///< e^{−T_s/τ}
    T            b{T{0}};                  ///< The delay-free model's gain on u_k: b_fresh + b_stale
    T            b_fresh{T{0}};            ///< The delayed model's gain on u_{k−d}, held from the switch
    T            b_stale{T{0}};            ///< Its gain on u_{k−d−1}, held until the switch (0 for whole samples)
    size_t       whole{0};                 ///< d, the delay's whole periods
    DesignStatus status{DesignStatus::ok}; ///< ok, invalid_sample_time or invalid_parameter
};

/**
 * @brief The model of a plant of gain @p K and time constant @p tau sampled at @p Ts, its input
 *        @p seconds late, in at most @p capacity delay states
 * @return status invalid_sample_time for a Ts that is not finite and positive; invalid_parameter
 *         for K not finite or 0, τ not finite and positive, a delay that is negative, not finite
 *         or longer than capacity·Ts, or a model gain that is 0 or not finite; else ok
 */
template<typename T>
[[nodiscard]] constexpr SmithModel<T> smith_model(T K, T tau, T seconds, T Ts, size_t capacity) {
    SmithModel<T> m{};
    if (!damp::finite_positive(Ts)) {
        m.status = DesignStatus::invalid_sample_time;
        return m;
    }
    const auto split = split_delay(seconds, Ts, capacity);
    if (!damp::isfinite(K) || K == T{0} || !damp::finite_positive(tau) || !split) {
        m.status = DesignStatus::invalid_parameter;
        return m;
    }
    m.whole = split->whole;
    m.a = damp::exp(-Ts / tau);
    const T held = damp::exp(-(Ts - split->fraction) / tau); // over the part of the period after the switch
    m.b_fresh = K * (T{1} - held);
    m.b_stale = (split->fraction > T{0}) ? (held * K * (T{1} - damp::exp(-split->fraction / tau))) : T{0};
    m.b = m.b_fresh + m.b_stale;
    if (!damp::isfinite(m.b) || !damp::isfinite(m.b_fresh) || !damp::isfinite(m.b_stale) || m.b == T{0}) {
        m.status = DesignStatus::invalid_parameter;
    }
    return m;
}

} // namespace design::detail

} // namespace damp
