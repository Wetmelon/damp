// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lqi.hpp
 * @brief The linear-quadratic regulator with integral action: LQI, the StateFeedback of design::state_feedback with an IntegralAction, and LQIResult
 *
 * design::state_feedback(plant, [Ts,] {.Q, .R, .N}, IntegralAction{.Q[, .outputs = {i, …}]}[, InputDelay<N>{θ}])
 * (state_feedback.hpp) → LQIResult → LQI; an IntegralAction that lists no outputs integrates every one.
 * The integrators are on the plant's outputs y = C·x + D·u, as MATLAB®'s lqi: every one, so the
 * plant's C names what r holds, or the ones IntegralAction lists, in that order.
 * Q_integral is the cost on the discrete integral state, not a block of an augmented Q on [x; ξ].
 * With an InputDelay the law holds the inputs still on their way, and update takes the plant state.
 * step(r, x) is predict, update(x), control(r). The low-latency tick is update(x), control(r),
 * back_calculate, predict. With feedforward from a planned trajectory:
 * control(r, x_ref, u_ref), step(r, x_ref, u_ref, x). While another source drives the plant,
 * track(u, r, x) takes step's place and preloads the integrators, so the first step starts from u.
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/lq_law.hpp" // IWYU pragma: export
#include "damp/controllers/detail/state_feedback_core.hpp"
#include "damp/controllers/law.hpp"                    // IWYU pragma: export
#include "damp/controllers/state_feedback_runtime.hpp" // IWYU pragma: export
#include "damp/controllers/state_feedback_types.hpp"   // IWYU pragma: export

namespace damp {

namespace design {

/** @brief an LQI design's result: the gain on [x; xi] and the integral action; with an input delay, on [x; v(k−1); …; v(k−NRU); xi]. */
template<size_t NX, size_t NU, size_t NI, typename T = double, size_t NY = NX, size_t NRU = 0>
using LQIResult = StateFeedbackResult<NX, NU, T, NI, NY, NoEstimator, 0, NRU>;

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Linear-quadratic regulator with integral action, on the measured state: the StateFeedback of a design with integrators and no filter
 *
 * Runs design::state_feedback with an IntegralAction. xi integrates r − Ci·C·x; the law is u = u* − K·[x − x*; xi]. The tick is
 * step(r, x), or update(x), control(r), back_calculate, predict() (predict advances xi), with
 * x* = 0 and u* = 0; feedforward from a planned
 * trajectory (x_ref, u_ref) is control(r, x_ref, u_ref), or step(r, x_ref, u_ref, x). The
 * command is clamped into the design's InputLimits, its move into an InputRateLimits, and xi into
 * IntegralAction::limits (unbounded by default; limits() reads them and set_limits() moves them at
 * run time), and predict() holds the integrators against those clamps; back_calculate(u_sat) before it reports a stage that clamps after
 * the law.
 *
 * A design with an InputDelay (NRU > 0) runs on [x; v(k−1); …; v(k−NRU); xi]: update takes the
 * plant state x, and predict shifts the input back_calculate recorded (else the last command)
 * into the delay states, so it runs every tick. x_ref is the plant's; the delay states' reference is u_ref.
 *
 * @note Compare with MATLAB®'s lqi.
 *
 * A design with a MeasuredDisturbance (NM > 0) feeds d forward through the steady state under it
 * that holds the integrated outputs at 0, u = u_t(d) − K·[x − x_t(d); xi], and the integrators
 * hold r and take up what the model leaves: control(r, d), step(r, d, x), control(r, x_ref,
 * u_ref, d), step(r, x_ref, u_ref, d, x).
 *
 * While another source drives the plant (an operator, another law), track(u, …) takes the place of
 * each step(…), its arguments step's with the input applied first: track(u, r, x),
 * track(u, r, x_ref, u_ref, x), track(u, r, d, x), track(u, r, x_ref, u_ref, d, x). It holds the
 * state, the delay states take u, and ξ is preloaded through (∂u/∂ξ)⁺ so that control here would
 * return u (u clamped into its input limits): exactly with NI ≥ NU, the least-squares ξ with fewer,
 * then clamped into the integrators' limits. The first step after it moves from u by what the
 * period changed, at most an InputRateLimits' rate. A design::Handover fades what the preload
 * leaves.
 *
 * The alias takes integrators (NI > 0) on the measured state, so `LQI law{res}` refuses a design
 * without them or with a filter. StateFeedback lists the members and the tick.
 *
 * @tparam NX    Plant states
 * @tparam NU    Inputs
 * @tparam NI    Integrators (at least one)
 * @tparam T     Scalar
 * @tparam NY    Measured outputs, the ones the integrators read through C
 * @tparam Store Where the matrices live (detail::HeldGains: members)
 * @tparam NRU   Input delay states per input the law holds (0: none)
 * @tparam NM    Measured-disturbance inputs (0: none)
 * @tparam NDU   Inputs with a move limit (an InputRateLimits: NU; 0: none)
 * @tparam HO    Whether the law fades what the preload leaves after a track (a Handover; false: none)
 */
template<size_t NX, size_t NU, size_t NI, typename T = float, size_t NY = NX, typename Store = detail::HeldGains, size_t NRU = 0, size_t NM = 0, size_t NDU = 0, bool HO = false>
    requires(NI > 0)
using LQI = StateFeedback<NX, NU, T, NI, NY, NoEstimator, 0, Store, NRU, NM, NDU, HO>;

} // namespace damp
