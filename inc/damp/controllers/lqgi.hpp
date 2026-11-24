// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lqgi.hpp
 * @brief The linear-quadratic Gaussian law with integral action: design::state_feedback with a KalmanNoise and an IntegralAction, LQGIResult, the LQGI runtime
 *
 * design::state_feedback(plant, [Ts,] {.Q, .R, .N}, KalmanNoise{.process, .measurement}, IntegralAction{.Q[, .outputs = {i, …}]}[, InputDelay<N>{θ}])
 * (state_feedback.hpp) → LQGIResult → LQGI. With an InputDelay the filter carries the inputs still on their way. step(r, y) is predict, update(y), control(r). The low-latency tick is
 * update(y), control(r), back_calculate, predict. With feedforward from a
 * planned trajectory: control(r, x_ref, u_ref), step(r, x_ref, u_ref, y).
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/lq_law.hpp" // IWYU pragma: export
#include "damp/controllers/detail/state_feedback_core.hpp"
#include "damp/controllers/law.hpp"                  // IWYU pragma: export
#include "damp/controllers/state_feedback_types.hpp" // IWYU pragma: export
#include "damp/estimation/kalman.hpp"                // IWYU pragma: export
#include "damp/matrix/colvec.hpp"

namespace damp {

namespace design {

/** @brief an LQGI design's result: the gain on [x; xi], the Kalman filter and the integral action. */
template<size_t NX, size_t NU, size_t NY, size_t NI, typename T = double, size_t NW = NX, size_t NV = NY>
using LQGIResult = StateFeedbackResult<NX, NU, T, NI, NY, KalmanResult<NX, NU, NY, T, NW, NV, 0>, 0>;

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Linear-quadratic Gaussian law with integral action on the measured outputs
 *
 * Runs design::state_feedback with a KalmanNoise and an IntegralAction. The estimator is a FixedGainObserver with design::kalman's fixed gain
 * (`estimator`). xi integrates r − Ci·y on the last update's measurement; the law is
 * u = u* − K·[x̂ − x*; xi]. The tick is step(r, y), or update(y), control(r), back_calculate,
 * predict() (predict moves the filter and xi), with x* = 0 and u* = 0; feedforward from a planned trajectory (x_ref, u_ref) is control(r, x_ref, u_ref), or
 * step(r, x_ref, u_ref, y). The command is clamped into u_lim and xi into integral_lim (the
 * design's, IntegralAction::lim; unbounded by default): predict() moves the filter on the clamped command
 * and holds the integrators against the clamp. back_calculate(u_sat) before it reports a stage
 * that clamps after the law.
 *
 * A design with an InputDelay (NRU > 0) has its filter on [x; v(k−1); …; v(k−NRU)]: predict
 * shifts the applied input into the delay states. x_ref is the plant's; the delay states' reference is
 * u_ref.
 *
 * @note Compare with MATLAB®'s lqgtrack.
 *
 * A design with a MeasuredDisturbance (NM > 0) feeds d forward as LQI does, and its filter takes d
 * as a known input after u: control(r, d), step(r, d, y), control(r, x_ref, u_ref, d),
 * step(r, x_ref, u_ref, d, y).
 *
 * @tparam NRU Input delay states of the filter's model (0: none)
 * @tparam NM  Measured-disturbance inputs (0: none)
 */
template<size_t NX, size_t NU, size_t NY, size_t NI, typename T, typename Est, typename Store = detail::HeldGains, size_t NRU = 0, size_t NM = 0>
struct LQGI : private detail::StateFeedbackCore<NX, NU, T, NI, NY, Est, 0, Store, NRU, NM> {
    using Core = detail::StateFeedbackCore<NX, NU, T, NI, NY, Est, 0, Store, NRU, NM>; ///< The shared implementation

    constexpr LQGI() = default;

    /** @brief From design::state_feedback(any scalar, design::IntegralAction{}). */
    template<detail::state_feedback_result R>
        requires(!Store::constant)
    constexpr LQGI(const R& result) { // NOLINT
        constexpr bool kind = R::kIntegral && R::kEstimated;
        static_assert(kind, "LQGI runs design::state_feedback with a KalmanNoise and an IntegralAction: this design has no Kalman filter (LQR, LQI) or no integral action (LQG, OFLQG).");
        if constexpr (kind) {
            static_cast<Core&>(*this) = Core(result);
        }
    }

    using Core::demand_lim;
    using Core::estimator;
    using Core::integral;
    using Core::integral_lim;
    using Core::K;
    using Core::tracking_time;
    using Core::u_lim;

    /** @brief The filter's measurement update with @p y (the integrators read it at the next control). */
    constexpr void update(const ColVec<NY, T>& y) { Core::update(y); }

    /** @brief u = −K·[x̂; xi]: the basic call, x* = 0 and u* = 0; holds r − Ci·y for the next predict. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NI, T>& r)
        requires(NM == 0)
    {
        return Core::control(r, ColVec<NX, T>{}, ColVec<NU, T>{});
    }

    /** @brief u = u* − K·[x̂ − x*; xi]: feedforward from a planned trajectory (x*, u*); holds r − Ci·y for the next predict. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NI, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref)
        requires(NM == 0)
    {
        return Core::control(r, x_ref, u_ref);
    }

    /** @brief One tick: predict(), update(@p y), then control(@p r). The basic call. */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NI, T>& r, const ColVec<NY, T>& y)
        requires(NM == 0)
    {
        return Core::step(r, ColVec<NX, T>{}, ColVec<NU, T>{}, y);
    }

    /** @brief One tick: predict(), update(@p y), then control(@p r, @p x_ref, @p u_ref), the feedforward form. */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NI, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NY, T>& y)
        requires(NM == 0)
    {
        return Core::step(r, x_ref, u_ref, y);
    }

    /** @brief control(@p r) under the measured disturbance @p d: the target under it, after the reference. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NI, T>& r, const ColVec<NM, T>& d)
        requires(NM > 0)
    {
        return Core::control(r, ColVec<NX, T>{}, ColVec<NU, T>{}, d);
    }

    /** @brief control(@p r, @p d) about the planned trajectory (@p x_ref, @p u_ref): the reference and its feedforward, then the measured disturbance. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NI, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NM, T>& d)
        requires(NM > 0)
    {
        return Core::control(r, x_ref, u_ref, d);
    }

    /** @brief One tick: predict(), update(@p y), then control(@p r, @p d). */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NI, T>& r, const ColVec<NM, T>& d, const ColVec<NY, T>& y)
        requires(NM > 0)
    {
        return Core::step(r, ColVec<NX, T>{}, ColVec<NU, T>{}, d, y);
    }

    /** @brief One tick: predict(), update(@p y), then control(@p r, @p x_ref, @p u_ref, @p d). */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NI, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NM, T>& d, const ColVec<NY, T>& y)
        requires(NM > 0)
    {
        return Core::step(r, x_ref, u_ref, d, y);
    }

    using Core::applied;
    using Core::back_calculate;
    using Core::getK;
    using Core::last_command;
    using Core::predict;
    using Core::reset;
    using Core::state_estimate;
    using Core::valid;
};

/** @brief LQGI law{res}: the sizes, scalar and filter of an LQGI design's result. */
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, typename Kalman, size_t NZ, size_t NRU, size_t NM>
LQGI(const design::StateFeedbackResult<NX, NU, T, NI, NY, Kalman, NZ, NRU, NM>&) -> LQGI<NX, NU, NY, NI, T, detail::runtime_estimator_t<Kalman, T>, detail::HeldGains, NRU, NM>;

namespace detail {

/** @brief An LQ design with a filter and integrators runs as LQGI. */
template<size_t NX, size_t NU, typename U, size_t NI, size_t NY, typename Kalman, size_t NRU, size_t NM, typename T>
    requires(NI > 0 && !damp::is_same_v<Kalman, NoEstimator>)
struct lq_runtime<design::StateFeedbackResult<NX, NU, U, NI, NY, Kalman, 0, NRU, NM>, T> {
    using type = LQGI<NX, NU, NY, NI, T, runtime_estimator_t<Kalman, T>, HeldGains, NRU, NM>; ///< LQGI
    template<typename Est, typename Store>
    using on_constants = LQGI<NX, NU, NY, NI, T, Est, Store, NRU, NM>; ///< LQGI on the design's constants
};

} // namespace detail

} // namespace damp
