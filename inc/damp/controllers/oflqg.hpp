// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file oflqg.hpp
 * @brief The offset-free linear-quadratic Gaussian law: design::oflqg, OFLQGResult, the OFLQG runtime
 *
 * design::oflqg(plant, [Ts,] {.Q, .R, .N}, KalmanNoise{.process, .measurement, .disturbance}, disturbance) → OFLQGResult
 * → OFLQG. step(r, y) is predict, update(y), control(r). The low-latency tick is update(y),
 * control(r), back_calculate, predict. With feedforward from a planned
 * trajectory: control(r, x_ref, u_ref), step(r, x_ref, u_ref, y).
 */

#include <cstddef>
#include <type_traits>

#include "damp/controllers/state_feedback_core.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace design {

/** @brief design::oflqg's result: the gain, the Kalman filter on [x; d] and the steady-state target. */
template<size_t NX, size_t NU, size_t NY, size_t ND, size_t NZ, typename T = double, size_t NW = NX, size_t NV = NY>
using OFLQGResult = StateFeedbackResult<NX, NU, T, 0, NY, KalmanResult<NX + ND, NU, NY, T, NW + ND, NV, ND>, NZ>;

/**
 * @brief The offset-free linear-quadratic Gaussian law: a Kalman filter on [x; d] and the steady-state
 *        target that holds the controlled outputs at r under d̂: runs as OFLQG.
 *
 * @p disturbance says where d enters: InputDisturbance (B_d = B), OutputDisturbance (a bias on each
 * reading), PlantDisturbance (the plant's own G and H) or DisturbanceModel (explicit, with the
 * outputs held offset-free). The filter's model is
 * @f[
 *   x_{k+1} = A x_k + B u_k + B_d d_k,\qquad d_{k+1} = d_k,\qquad y_k = C x_k + C_d d_k + D u_k ,
 * @f]
 * the random walk of d is @p noise .disturbance (positive definite), and the target solves
 * @f[
 *   \begin{bmatrix} A-I & B \\ C_z & D_z \end{bmatrix}
 *   \begin{bmatrix} x_t \\ u_t \end{bmatrix}
 *   =
 *   \begin{bmatrix} -B_d\hat d \\ r - C_{d,z}\hat d \end{bmatrix},
 * @f]
 * so the law is @f$ u = u_t - K(\hat x - x_t) @f$ with the plain gain. The augment must be
 * detectable (ND ≤ NY), else DesignStatus::not_detectable; a target needs a unique steady state
 * for NZ = NU, else DesignStatus::singular.
 *
 * @note Compare with MATLAB®'s mpc with its default output disturbance model, without a QP.
 * @see "Disturbance Models for Offset-Free Model-Predictive Control" (Pannocchia & Rawlings,
 *      2003), AIChE Journal 49(2)
 * @see "Offset-free reference tracking with model predictive control" (Maeder, Borrelli &
 *      Morari, 2009), Automatica 45(9)
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Disturbance, size_t ND = detail::disturbance_traits<Disturbance, NU, NY, NW>::states>
    requires(detail::is_disturbance_model<Disturbance>)
[[nodiscard]] constexpr auto oflqg(
    const StateSpace<NX, NU, NY, T, NW, NV>&                sys,
    const LQRWeights<NX, NU, 0, T>&                         weights,
    const std::type_identity_t<KalmanNoise<NW, NV, T, ND>>& noise,
    const Disturbance&                                      disturbance
) {
    return detail::state_feedback(sys, weights, detail::StateFeedbackOptions{.noise = noise, .offset_free = disturbance});
}

/** @brief The sampled-data offset-free linear-quadratic Gaussian law (the process and disturbance noise densities): runs as OFLQG. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Disturbance, size_t ND = detail::disturbance_traits<Disturbance, NU, NY, NW>::states>
    requires(detail::is_disturbance_model<Disturbance>)
[[nodiscard]] constexpr auto oflqg(
    const StateSpace<NX, NU, NY, T, NW, NV>&                sys,
    std::type_identity_t<T>                                 Ts,
    const LQRWeights<NX, NU, 0, T>&                         weights,
    const std::type_identity_t<KalmanNoise<NW, NV, T, ND>>& noise,
    const Disturbance&                                      disturbance
) {
    return detail::state_feedback(sys, Ts, weights, detail::StateFeedbackOptions{.noise = noise, .offset_free = disturbance});
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Offset-free linear-quadratic Gaussian law: a Kalman filter on [x; d] and the steady-state target
 *
 * Runs design::oflqg. The estimator is a Luenberger observer on [x; d] with design::kalman's fixed
 * gain, holding deviations from the target (`estimator`). The target holds the controlled outputs
 * at r under the disturbance estimate d̂, taken about the reference point (x*, u*). The tick is
 * update(y) then control(r), or step(r, y), with x* = 0 and u* = 0; feedforward from a planned
 * trajectory (x_ref, u_ref) is control(r, x_ref, u_ref), or step(r, x_ref, u_ref, y). Then
 * back_calculate(u_sat) and predict().
 *
 * @note Compare with MATLAB®'s lqg with an integrated-white-noise disturbance model (setoutdist / setindist).
 */
template<size_t NX, size_t NU, size_t NY, size_t NZ, typename T, typename Est, typename Store = detail::HeldGains>
struct OFLQG : private detail::StateFeedbackCore<NX, NU, T, 0, NY, Est, NZ, Store> {
    using Core = detail::StateFeedbackCore<NX, NU, T, 0, NY, Est, NZ, Store>; ///< The shared implementation

    constexpr OFLQG() = default;

    /** @brief From design::oflqg (any scalar). */
    template<size_t NXR, size_t NUR, typename U, size_t NIR, size_t NYR, typename Kalman, size_t NZR>
        requires(!Store::constant)
    constexpr OFLQG(const design::StateFeedbackResult<NXR, NUR, U, NIR, NYR, Kalman, NZR>& result) { // NOLINT
        constexpr bool kind = design::StateFeedbackResult<NXR, NUR, U, NIR, NYR, Kalman, NZR>::ND > 0;
        static_assert(kind, "OFLQG runs design::oflqg: this design has no disturbance model (LQR, LQI, LQG, LQGI).");
        if constexpr (kind) {
            static_cast<Core&>(*this) = Core(result);
        }
    }

    using Core::estimator;
    using Core::K;
    using Core::target;

    /** @brief The filter's measurement update with @p y. */
    constexpr void update(const ColVec<NY, T>& y) { Core::update(y); }

    /** @brief The target for @p r under d̂ about the origin, then the law about it: the basic call. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NZ, T>& r) { return Core::control(r, ColVec<NX, T>{}, ColVec<NU, T>{}); }

    /** @brief The target for @p r under d̂ about (@p x_ref, @p u_ref), then the law about it: feedforward from a planned trajectory. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NZ, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref) {
        return Core::control(r, x_ref, u_ref);
    }

    /** @brief One tick: predict(), update(@p y), then control(@p r). The basic call. */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NZ, T>& r, const ColVec<NY, T>& y) { return Core::step(r, ColVec<NX, T>{}, ColVec<NU, T>{}, y); }

    /** @brief One tick: predict(), update(@p y), then control(@p r, @p x_ref, @p u_ref), the feedforward form. */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NZ, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NY, T>& y) {
        return Core::step(r, x_ref, u_ref, y);
    }

    using Core::applied;
    using Core::back_calculate;
    using Core::disturbance_estimate;
    using Core::getK;
    using Core::last_command;
    using Core::predict;
    using Core::reset;
    using Core::state_estimate;
};

/** @brief OFLQG law{res}: the sizes, scalar, filter and target of design::oflqg's result. */
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, typename Kalman, size_t NZ>
OFLQG(const design::StateFeedbackResult<NX, NU, T, NI, NY, Kalman, NZ>&) -> OFLQG<NX, NU, NY, NZ, T, detail::runtime_estimator_t<Kalman, T>>;


namespace detail {

/** @brief An LQ design with a disturbance model runs as OFLQG. */
template<size_t NX, size_t NU, typename U, size_t NY, typename Kalman, size_t NZ, typename T>
    requires(design::StateFeedbackResult<NX, NU, U, 0, NY, Kalman, NZ>::ND > 0)
struct lq_runtime<design::StateFeedbackResult<NX, NU, U, 0, NY, Kalman, NZ>, T> {
    using type = OFLQG<NX, NU, NY, NZ, T, runtime_estimator_t<Kalman, T>>; ///< OFLQG
    template<typename Est, typename Store>
    using on_constants = OFLQG<NX, NU, NY, NZ, T, Est, Store>; ///< OFLQG on the design's constants
};

} // namespace detail

} // namespace damp
