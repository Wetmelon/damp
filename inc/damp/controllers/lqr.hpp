// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lqr.hpp
 * @brief The linear-quadratic regulator: LQR, the StateFeedback of design::state_feedback's law on the measured state
 *
 * design::state_feedback(plant, [Ts,] {.Q, .R, .N} or Poles{…}[, TrackedOutputs{{i, …}}][, MeasuredDisturbance{Bd}][, InputDelay<N>{θ}])
 * (state_feedback.hpp) → StateFeedbackResult → LQR. discrete_lqr and design::place (riccati.hpp,
 * pole_placement.hpp) are the bare gains: no runtime takes their results.
 * With an InputDelay the LQR holds the inputs still on their way. step(x) is predict, update(x), control
 * (predict is empty on LQR without an InputDelay). The low-latency tick is update(x), control(),
 * then back_calculate and predict() with an InputDelay. While another source drives the plant,
 * track(u, …) takes step's place, its arguments step's with the input applied first. With tracked
 * outputs: control(r), step(r, x). With feedforward from a planned trajectory: control([r,]
 * x_ref, u_ref), step([r,] x_ref, u_ref, x). With a measured disturbance d, after the reference
 * and its feedforward: control([r,] [x_ref, u_ref,] d), step([r,] [x_ref, u_ref,] d, x).
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/lq_law.hpp" // IWYU pragma: export
#include "damp/controllers/detail/state_feedback_core.hpp"
#include "damp/controllers/law.hpp"                    // IWYU pragma: export
#include "damp/controllers/state_feedback_runtime.hpp" // IWYU pragma: export
#include "damp/controllers/state_feedback_types.hpp"   // IWYU pragma: export
#include "damp/design/riccati.hpp"

namespace damp {

/**
 * @ingroup discrete_controllers
 * @brief Linear-quadratic regulator on the measured state, u = u* − K·(x − x*): the StateFeedback of a design with no integrators or filter
 *
 * Runs design::state_feedback, with LQRWeights or Poles, or a gain K of its own: `LQR law{res.as<float>()}`, or
 * `Law<kRes, float>` for a struct member. step(x) is predict, update(x), control; predict is empty
 * without an InputDelay, so that equals update then control. The regulator is x* = 0 and u* = 0. A design with
 * TrackedOutputs holds those outputs at r through its steady-state target: control(r), step(r, x).
 * Feedforward from a planned trajectory (x_ref, u_ref) adds to it: control([r,] x_ref, u_ref),
 * step([r,] x_ref, u_ref, x). The command is clamped into design::state_feedback's InputLimits, and
 * its move into an InputRateLimits (unbounded by default, and for a gain K of its own), and a
 * DemandLimit bounds what the outer terms ask of a state: limits() reads them and set_limits()
 * moves them at run time, as on an MPC. A failed or continuous design does not compile when
 * constant-evaluated.
 *
 * @code
 * static constinit LQR law{kRes.as<float>()};
 * law.update(x);
 * u = law.control(); // the regulator
 * @endcode
 *
 * @note Compare with MATLAB®'s lqr / lqrd.
 *
 * A design with an InputDelay (NRU > 0) runs on [x; v(k−1); …; v(k−NRU)]: update takes the plant
 * state, and predict shifts the applied input into the delay states, so it is no longer empty.
 *
 * A design with a MeasuredDisturbance (NM > 0) takes d each tick and holds the target under it,
 * u = u_t(r, d) − K·(x − x_t(r, d)): control(r, d), step(r, d, x), and with feedforward
 * control(r, x_ref, u_ref, d), step(r, x_ref, u_ref, d, x). The regulator (NZ = 0) holds every
 * output at 0 under d: control(d), step(d, x), control(x_ref, u_ref, d), step(x_ref, u_ref, d, x).
 *
 * While another source drives the plant (an operator, another law), track(u, …) takes the place of
 * each step(…), its arguments step's with the input applied first: track(u, x), track(u, r, x),
 * track(u, x_ref, u_ref, x), track(u, d, x), …. It holds the state, the delay states take u, and an
 * InputRateLimits' first move is from u. The regulator has nothing to preload, so its first step
 * after a track is its own command, unless the design took a design::Handover: then that command is
 * u, and the mismatch with the law's own decays by exp(−Ts/fade_time) each tick, within its input
 * limits.
 *
 * The alias's pattern is its constraint: no integrators and no filter, so `LQR law{res}` refuses
 * any other kind of design. StateFeedback lists the members and the tick.
 *
 * @tparam NX    Plant states
 * @tparam NU    Inputs
 * @tparam T     Scalar
 * @tparam NY    Measured outputs
 * @tparam NZ    Tracked outputs: the size of r (0: none, the regulator)
 * @tparam Store Where the matrices live (detail::HeldGains: members)
 * @tparam NRU   Input delay states per input the law holds (0: none)
 * @tparam NM    Measured-disturbance inputs (0: none)
 * @tparam NDU   Inputs with a move limit (an InputRateLimits: NU; 0: none)
 * @tparam HO    Whether the law fades the mismatch after a track (a Handover; false: none)
 */
template<size_t NX, size_t NU, typename T = float, size_t NY = NX, size_t NZ = 0, typename Store = detail::HeldGains, size_t NRU = 0, size_t NM = 0, size_t NDU = 0, bool HO = false>
using LQR = StateFeedback<NX, NU, T, 0, NY, NoEstimator, NZ, Store, NRU, NM, NDU, HO>;

} // namespace damp
