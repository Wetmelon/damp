// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lqi.hpp
 * @brief The linear-quadratic regulator with integral action: design::lqi, LQIResult, the LQI runtime
 *
 * design::lqi(plant, [Ts,] {.Q, .R, .N, .Q_integral}[, IntegralAction<NI>{{i, …}}]) → LQIResult → LQI.
 * Q_integral is the cost on the discrete integral state, not a block of an augmented Q on [x; ξ].
 * step(r, x) is predict, update(x), control(r). The low-latency tick is update(x), control(r),
 * back_calculate, predict. With feedforward from a planned trajectory:
 * control(r, x_ref, u_ref), step(r, x_ref, u_ref, x).
 */

#include <cstddef>
#include <type_traits>

#include "damp/controllers/state_feedback_core.hpp" // IWYU pragma: export
#include "damp/matrix/colvec.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace design {

/** @brief design::lqi's result: the gain on [x; xi] and the integral action. */
template<size_t NX, size_t NU, size_t NI, typename T = double, size_t NY = NX>
using LQIResult = StateFeedbackResult<NX, NU, T, NI, NY>;

/**
 * @brief The linear-quadratic regulator with integral action on every output (xi += Ts·(r − Ci·y)): runs as LQI.
 * @note Compare with MATLAB®'s lqi(sys, Q, R, N).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr auto lqi(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const LQRWeights<NX, NU, NY, T>& weights) {
    return detail::state_feedback(sys, weights, detail::StateFeedbackOptions{.offset_free = IntegralAction{}});
}

/** @brief The sampled-data linear-quadratic regulator with integral action on every output: runs as LQI. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr auto lqi(const StateSpace<NX, NU, NY, T, NW, NV>& sys, std::type_identity_t<T> Ts, const LQRWeights<NX, NU, NY, T>& weights) {
    return detail::state_feedback(sys, Ts, weights, detail::StateFeedbackOptions{.offset_free = IntegralAction{}});
}

/**
 * @brief The linear-quadratic regulator with integral action on the listed outputs
 *        (`IntegralAction<2>{{0, 3}}`, in that order: r follows it): runs as LQI.
 */
template<size_t NI, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
    requires(NI > 0)
[[nodiscard]] constexpr auto lqi(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const LQRWeights<NX, NU, NI, T>&         weights,
    const IntegralAction<NI>&                outputs
) {
    return detail::state_feedback(sys, weights, detail::StateFeedbackOptions{.offset_free = outputs});
}

/** @brief The sampled-data linear-quadratic regulator with integral action on the listed outputs: runs as LQI. */
template<size_t NI, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
    requires(NI > 0)
[[nodiscard]] constexpr auto lqi(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    std::type_identity_t<T>                  Ts,
    const LQRWeights<NX, NU, NI, T>&         weights,
    const IntegralAction<NI>&                outputs
) {
    return detail::state_feedback(sys, Ts, weights, detail::StateFeedbackOptions{.offset_free = outputs});
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Linear-quadratic regulator with integral action, on the measured state
 *
 * Runs design::lqi. xi integrates r − Ci·C·x; the law is u = u* − K·[x − x*; xi]. The tick is
 * update(x) then control(r), or step(r, x), with x* = 0 and u* = 0; feedforward from a planned
 * trajectory (x_ref, u_ref) is control(r, x_ref, u_ref), or step(r, x_ref, u_ref, x). Then
 * back_calculate(u_sat) and predict() hold the integrators against a clamp.
 *
 * @note Compare with MATLAB®'s lqi.
 */
template<size_t NX, size_t NU, size_t NI, typename T = float, size_t NY = NX, typename Store = detail::HeldGains>
struct LQI : private detail::StateFeedbackCore<NX, NU, T, NI, NY, NoEstimator, 0, Store> {
    using Core = detail::StateFeedbackCore<NX, NU, T, NI, NY, NoEstimator, 0, Store>; ///< The shared implementation

    constexpr LQI() = default;

    /** @brief From design::lqi (any scalar). */
    template<size_t NXR, size_t NUR, typename U, size_t NIR, size_t NYR, typename Kalman, size_t NZR>
        requires(!Store::constant)
    constexpr LQI(const design::StateFeedbackResult<NXR, NUR, U, NIR, NYR, Kalman, NZR>& result) { // NOLINT
        constexpr bool kind = NIR > 0 && std::is_same_v<Kalman, NoEstimator>;
        static_assert(kind, "LQI runs design::lqi: this design has no integral action (LQR) or a Kalman filter (LQG, LQGI, OFLQG).");
        if constexpr (kind) {
            static_cast<Core&>(*this) = Core(result);
        }
    }

    using Core::integral;
    using Core::K;

    /** @brief Hold the measured state @p x for the next control. */
    constexpr void update(const ColVec<NX, T>& x) { Core::update(x); }

    /** @brief u = −K·[x; xi] at the state of the last update: the basic call, x* = 0 and u* = 0; holds r − Ci·C·x for the next predict. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NI, T>& r) { return Core::control(r, ColVec<NX, T>{}, ColVec<NU, T>{}); }

    /** @brief u = u* − K·[x − x*; xi] at the state of the last update: feedforward from a planned trajectory (x*, u*). */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NI, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref) {
        return Core::control(r, x_ref, u_ref);
    }

    /** @brief One tick: predict(), update(@p x), then control(@p r). The basic call. */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NI, T>& r, const ColVec<NX, T>& x) { return Core::step(r, ColVec<NX, T>{}, ColVec<NU, T>{}, x); }

    /** @brief One tick: predict(), update(@p x), then control(@p r, @p x_ref, @p u_ref), the feedforward form. */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NI, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NX, T>& x) {
        return Core::step(r, x_ref, u_ref, x);
    }

    using Core::applied;
    using Core::back_calculate;
    using Core::getK;
    using Core::last_command;
    using Core::predict;
    using Core::reset;
};

/** @brief LQI law{res}: the sizes and scalar of design::lqi's result. */
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, typename Kalman, size_t NZ>
LQI(const design::StateFeedbackResult<NX, NU, T, NI, NY, Kalman, NZ>&) -> LQI<NX, NU, NI, T, NY>;


namespace detail {

/** @brief An LQ design with integrators and no filter runs as LQI. */
template<size_t NX, size_t NU, typename U, size_t NI, size_t NY, typename T>
    requires(NI > 0)
struct lq_runtime<design::StateFeedbackResult<NX, NU, U, NI, NY, NoEstimator, 0>, T> {
    using type = LQI<NX, NU, NI, T, NY>; ///< LQI
    template<typename, typename Store>
    using on_constants = LQI<NX, NU, NI, T, NY, Store>; ///< LQI on the design's constants
};

} // namespace detail

} // namespace damp
