// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lqi.hpp
 * @brief The linear-quadratic regulator with integral action: design::state_feedback with an IntegralAction, LQIResult, the LQI runtime
 *
 * design::state_feedback(plant, [Ts,] {.Q, .R, .N}, IntegralAction{.Q[, .outputs = {i, …}]}[, InputDelay<N>{θ}])
 * (state_feedback.hpp) → LQIResult → LQI; an IntegralAction that lists no outputs integrates every one.
 * The integrators are on the plant's outputs y = C·x + D·u, as MATLAB®'s lqi: every one, so the
 * plant's C names what r holds, or the ones IntegralAction lists, in that order.
 * Q_integral is the cost on the discrete integral state, not a block of an augmented Q on [x; ξ].
 * With an InputDelay the law holds the inputs still on their way, and update takes the plant state.
 * step(r, x) is predict, update(x), control(r). The low-latency tick is update(x), control(r),
 * back_calculate, predict. With feedforward from a planned trajectory:
 * control(r, x_ref, u_ref), step(r, x_ref, u_ref, x).
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/lq_law.hpp" // IWYU pragma: export
#include "damp/controllers/detail/state_feedback_core.hpp"
#include "damp/controllers/law.hpp"                  // IWYU pragma: export
#include "damp/controllers/state_feedback_types.hpp" // IWYU pragma: export
#include "damp/matrix/colvec.hpp"

namespace damp {

namespace design {

/** @brief an LQI design's result: the gain on [x; xi] and the integral action; with an input delay, on [x; v(k−1); …; v(k−NRU); xi]. */
template<size_t NX, size_t NU, size_t NI, typename T = double, size_t NY = NX, size_t NRU = 0>
using LQIResult = StateFeedbackResult<NX, NU, T, NI, NY, NoEstimator, 0, NRU>;

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Linear-quadratic regulator with integral action, on the measured state
 *
 * Runs design::state_feedback with an IntegralAction. xi integrates r − Ci·C·x; the law is u = u* − K·[x − x*; xi]. The tick is
 * step(r, x), or update(x), control(r), back_calculate, predict() (predict advances xi), with
 * x* = 0 and u* = 0; feedforward from a planned
 * trajectory (x_ref, u_ref) is control(r, x_ref, u_ref), or step(r, x_ref, u_ref, x). The
 * command is clamped into u_lim and xi into integral_lim (the design's, IntegralAction::lim; unbounded by
 * default), and predict() holds the integrators against that clamp; back_calculate(u_sat) before
 * it reports a stage that clamps after the law.
 *
 * A design with an InputDelay (NRU > 0) runs on [x; v(k−1); …; v(k−NRU); xi]: update takes the
 * plant state x, and predict shifts the input back_calculate recorded (else the last command)
 * into the registers, so it runs every tick. x_ref is the plant's; the registers' reference is u_ref.
 *
 * @note Compare with MATLAB®'s lqi.
 *
 * @tparam NRU Input registers the law holds (0: none)
 */
template<size_t NX, size_t NU, size_t NI, typename T = float, size_t NY = NX, typename Store = detail::HeldGains, size_t NRU = 0>
struct LQI : private detail::StateFeedbackCore<NX, NU, T, NI, NY, NoEstimator, 0, Store, NRU> {
    using Core = detail::StateFeedbackCore<NX, NU, T, NI, NY, NoEstimator, 0, Store, NRU>; ///< The shared implementation

    constexpr LQI() = default;

    /** @brief From design::state_feedback(any scalar, design::IntegralAction{}). */
    template<detail::state_feedback_result R>
        requires(!Store::constant)
    constexpr LQI(const R& result) { // NOLINT
        constexpr bool kind = R::kIntegral && !R::kEstimated;
        static_assert(kind, "LQI runs design::state_feedback with an IntegralAction: this design has no integral action (LQR) or a Kalman filter (LQG, LQGI, OFLQG).");
        if constexpr (kind) {
            static_cast<Core&>(*this) = Core(result);
        }
    }

    using Core::demand_lim;
    using Core::integral;
    using Core::integral_lim;
    using Core::K;
    using Core::tracking_time;
    using Core::u_lim;

    /** @brief Hold the measured plant state @p x for the next control (the input registers are the law's). */
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
    using Core::valid;
};

/** @brief LQI law{res}: the sizes and scalar of an LQI design's result. */
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, typename Kalman, size_t NZ, size_t NRU>
LQI(const design::StateFeedbackResult<NX, NU, T, NI, NY, Kalman, NZ, NRU>&) -> LQI<NX, NU, NI, T, NY, detail::HeldGains, NRU>;

namespace detail {

/** @brief An LQ design with integrators and no filter runs as LQI. */
template<size_t NX, size_t NU, typename U, size_t NI, size_t NY, size_t NRU, typename T>
    requires(NI > 0)
struct lq_runtime<design::StateFeedbackResult<NX, NU, U, NI, NY, NoEstimator, 0, NRU>, T> {
    using type = LQI<NX, NU, NI, T, NY, HeldGains, NRU>; ///< LQI
    template<typename, typename Store>
    using on_constants = LQI<NX, NU, NI, T, NY, Store, NRU>; ///< LQI on the design's constants
};

} // namespace detail

} // namespace damp
