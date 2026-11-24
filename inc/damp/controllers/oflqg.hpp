// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file oflqg.hpp
 * @brief The offset-free linear-quadratic Gaussian law: OFLQG, the StateFeedback of design::state_feedback with a KalmanNoise and a disturbance spec
 *
 * design::state_feedback(plant, [Ts,] {.Q, .R, .N}, KalmanNoise{.process, .measurement, .disturbance}, disturbance[, InputDelay<N>{θ}])
 * (state_feedback.hpp) → StateFeedbackResult → OFLQG. With an InputDelay the filter carries the inputs still on their way. step(r, y) is predict, update(y), control(r). The low-latency tick is update(y),
 * control(r), back_calculate, predict. With feedforward from a planned
 * trajectory: control(r, x_ref, u_ref), step(r, x_ref, u_ref, y). While another source drives the
 * plant, track(u, r, y) takes step's place: d̂ learns the load under u.
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/lq_law.hpp" // IWYU pragma: export
#include "damp/controllers/detail/state_feedback_core.hpp"
#include "damp/controllers/law.hpp"                    // IWYU pragma: export
#include "damp/controllers/state_feedback_runtime.hpp" // IWYU pragma: export
#include "damp/controllers/state_feedback_types.hpp"   // IWYU pragma: export
#include "damp/design/disturbance.hpp"                 // IWYU pragma: export
#include "damp/estimation/kalman.hpp"                  // IWYU pragma: export

namespace damp {

/**
 * @ingroup discrete_controllers
 * @brief Offset-free linear-quadratic Gaussian law, a Kalman filter on [x; d] and the steady-state target: the StateFeedback of a design with a disturbance model
 *
 * Runs design::state_feedback with a KalmanNoise and a disturbance spec. The estimator is a FixedGainObserver on [x; d] with design::kalman's fixed
 * gain, holding deviations from the target (`estimator`). The target holds the controlled outputs
 * at r under the disturbance estimate d̂, taken about the reference point (x*, u*). The tick is
 * step(r, y), or update(y), control(r), back_calculate, predict() (predict moves the filter),
 * with x* = 0 and u* = 0; feedforward from a planned
 * trajectory (x_ref, u_ref) is control(r, x_ref, u_ref), or step(r, x_ref, u_ref, y). The
 * command is clamped into the design's InputLimits and its move into an InputRateLimits (unbounded
 * by default), and predict() moves the filter with it, so the clamp is not taken for a disturbance;
 * back_calculate(u_sat) before it reports a stage that clamps after the law. update() clamps d̂
 * into StateFeedbackLimits::disturbance (the runtime's alone; unbounded by default), in d's units:
 * the most of a load the law leans against, as the integrators' limits on LQI. limits() reads them and set_limits() moves them at run time.
 *
 * A design with an InputDelay (NRU > 0) has its filter on [x; v(k−1); …; v(k−NRU); d]: predict
 * shifts the applied input into the delay states. x_ref is the plant's; the delay states' reference is
 * u_ref.
 *
 * @note Compare with MATLAB®'s lqg with an integrated-white-noise disturbance model (setoutdist / setindist).
 *
 * A design with a MeasuredDisturbance (NM > 0) has its filter take d as a known input after u and
 * its target hold r under [d̂; d], so d̂ estimates only the part of a load d does not measure:
 * control(r, d), step(r, d, y), control(r, x_ref, u_ref, d), step(r, x_ref, u_ref, d, y).
 *
 * While another source drives the plant (an operator, another law), track(u, …) takes the place of
 * each step(…), its arguments step's with the input applied first: track(u, r, y),
 * track(u, r, x_ref, u_ref, y), track(u, r, d, y), track(u, r, x_ref, u_ref, d, y). The filter
 * predicts with u, so d̂ learns the load the other source holds against: at a steady state the
 * target's input is then u and the first step returns it. Away from one the law's command differs
 * by its transient, which a design::Handover fades by exp(−Ts/fade_time) each tick, within its
 * input limits.
 *
 * The alias takes a filter that estimates a disturbance (Est a FixedGainObserver with ND > 0), so
 * `OFLQG law{res}` refuses a design without a disturbance model. StateFeedback lists the members
 * and the tick.
 *
 * @tparam NX    Plant states
 * @tparam NU    Inputs
 * @tparam NY    Measured outputs
 * @tparam NZ    Controlled outputs: the size of r
 * @tparam T     Scalar
 * @tparam Est   The filter: a FixedGainObserver on [x; d]
 * @tparam Store Where the matrices live (detail::HeldGains: members)
 * @tparam NRU   Input delay states per input of the filter's model (0: none)
 * @tparam NM    Measured-disturbance inputs (0: none)
 * @tparam NDU   Inputs with a move limit (an InputRateLimits: NU; 0: none)
 * @tparam HO    Whether the law fades the mismatch after a track (a Handover; false: none)
 */
template<size_t NX, size_t NU, size_t NY, size_t NZ, typename T, typename Est, typename Store = detail::HeldGains, size_t NRU = 0, size_t NM = 0, size_t NDU = 0, bool HO = false>
    requires(detail::estimates_disturbance<Est>)
using OFLQG = StateFeedback<NX, NU, T, 0, NY, Est, NZ, Store, NRU, NM, NDU, HO>;

} // namespace damp
