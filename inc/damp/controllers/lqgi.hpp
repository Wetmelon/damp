// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lqgi.hpp
 * @brief The linear-quadratic Gaussian law with integral action: LQGI, the StateFeedback of design::state_feedback with a KalmanNoise and an IntegralAction, and LQGIResult
 *
 * design::state_feedback(plant, [Ts,] {.Q, .R, .N}, KalmanNoise{.process, .measurement}, IntegralAction{.Q[, .outputs = {i, …}]}[, InputDelay<N>{θ}])
 * (state_feedback.hpp) → LQGIResult → LQGI. With an InputDelay the filter carries the inputs still on their way. step(r, y) is predict, update(y), control(r). The low-latency tick is
 * update(y), control(r), back_calculate, predict. With feedforward from a
 * planned trajectory: control(r, x_ref, u_ref), step(r, x_ref, u_ref, y). While another source
 * drives the plant, track(u, r, y) takes step's place and preloads the integrators.
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/lq_law.hpp" // IWYU pragma: export
#include "damp/controllers/detail/state_feedback_core.hpp"
#include "damp/controllers/law.hpp"                    // IWYU pragma: export
#include "damp/controllers/state_feedback_runtime.hpp" // IWYU pragma: export
#include "damp/controllers/state_feedback_types.hpp"   // IWYU pragma: export
#include "damp/estimation/kalman.hpp"                  // IWYU pragma: export

namespace damp {

namespace design {

/** @brief an LQGI design's result: the gain on [x; xi], the Kalman filter and the integral action. */
template<size_t NX, size_t NU, size_t NY, size_t NI, typename T = double, size_t NW = NX, size_t NV = NY>
using LQGIResult = StateFeedbackResult<NX, NU, T, NI, NY, KalmanResult<NX, NU, NY, T, NW, NV, 0>, 0>;

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Linear-quadratic Gaussian law with integral action on the measured outputs: the StateFeedback of a design with a filter and integrators
 *
 * Runs design::state_feedback with a KalmanNoise and an IntegralAction. The estimator is a FixedGainObserver with design::kalman's fixed gain
 * (`estimator`). xi integrates r − Ci·y on the last update's measurement; the law is
 * u = u* − K·[x̂ − x*; xi]. The tick is step(r, y), or update(y), control(r), back_calculate,
 * predict() (predict moves the filter and xi), with x* = 0 and u* = 0; feedforward from a planned trajectory (x_ref, u_ref) is control(r, x_ref, u_ref), or
 * step(r, x_ref, u_ref, y). The command is clamped into the design's InputLimits, its move into an
 * InputRateLimits, and xi into IntegralAction::limits (unbounded by default; limits() reads them
 * and set_limits() moves them at run time): predict() moves the filter on the clamped command
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
 * While another source drives the plant (an operator, another law), track(u, …) takes the place of
 * each step(…), its arguments step's with the input applied first: track(u, r, y),
 * track(u, r, x_ref, u_ref, y), track(u, r, d, y), track(u, r, x_ref, u_ref, d, y). The filter
 * predicts with u (and the measured d), and ξ is preloaded on its estimate as on LQI, so that
 * control here would return u (u clamped into its input limits): exactly with NI ≥ NU, the
 * least-squares ξ with fewer. A design::Handover fades what the preload leaves.
 *
 * The alias takes a filter (Est a FixedGainObserver) that estimates no disturbance, and
 * integrators (NI > 0), so `LQGI law{res}` refuses a design without either. StateFeedback lists
 * the members and the tick.
 *
 * @tparam NX    Plant states
 * @tparam NU    Inputs
 * @tparam NY    Measured outputs
 * @tparam NI    Integrators (at least one)
 * @tparam T     Scalar
 * @tparam Est   The filter: a FixedGainObserver on the law's model, with no disturbance states
 * @tparam Store Where the matrices live (detail::HeldGains: members)
 * @tparam NRU   Input delay states per input of the filter's model (0: none)
 * @tparam NM    Measured-disturbance inputs (0: none)
 * @tparam NDU   Inputs with a move limit (an InputRateLimits: NU; 0: none)
 * @tparam HO    Whether the law fades what the preload leaves after a track (a Handover; false: none)
 */
template<size_t NX, size_t NU, size_t NY, size_t NI, typename T, typename Est, typename Store = detail::HeldGains, size_t NRU = 0, size_t NM = 0, size_t NDU = 0, bool HO = false>
    requires(NI > 0 && detail::runs_filter<Est> && !detail::estimates_disturbance<Est>)
using LQGI = StateFeedback<NX, NU, T, NI, NY, Est, 0, Store, NRU, NM, NDU, HO>;

} // namespace damp
