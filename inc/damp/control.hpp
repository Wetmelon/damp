// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file control.hpp
 * @brief Embeddable core — matrix, math, design-is-deploy control/estimation, CST design.
 *
 * Compile-time control design in the same tree as firmware
 * (variant gains as @c constexpr), with the same runtime objects usable in SIL.
 *
 * Start here on the MCU: @c StateSpace (or TF / ZPK), @c design::lqr / @c lqi / @c pidtune / @c kalman,
 * @c Law{res.as<float>()}, then @c step in the interrupt. Same objects in SIL.
 *
 * This umbrella is the embeddable surface:
 * - linear algebra + pluggable math backend (via matrix / math includes)
 * - LTI systems (SS / TF / ZPK) + @c design:: synthesis (@c discretize, @c InputDelay ≤ one sample)
 * - core design-is-deploy controllers and estimators (including MPC / MHE on the
 *   allocation-free QP, ESKF / attitude fusion and the strapdown INS)
 * - identification on target: excitation sources, frequency-response and relay
 *   experiments, grey-box parameter fits
 * - common filters, fixed-step ODE integrators, geometry / Clarke-Park, toolbox helpers
 *
 * Nothing reachable from this header allocates on the heap or pulls a
 * third-party dependency. Host tools — frequency-domain analysis, ODE solvers,
 * closed-loop simulation, MATLAB®-style aliases — live behind @ref workbench.hpp.
 *
 * @code
 * #include "damp/control.hpp"          // MCU / embeddable path
 * using namespace damp;
 *
 * constexpr auto res = design::lqr(sys_d, {.Q = Q, .R = R});
 * constinit LQR controller{res.as<float>()}; // a failed design does not compile
 * // ISR: auto u = controller.step(x);
 * @endcode
 *
 * @see workbench.hpp for host analysis / solvers / simulate / matlab::
 */

// --- Linear algebra + math core ---------------------------------------------
#include "damp/matrix/matrix.hpp" // IWYU pragma: export

// --- Runtime taxonomy (controller/estimator/source concepts) -----------------
#include "damp/concepts.hpp" // IWYU pragma: export

// --- LTI system types -------------------------------------------------------
#include "damp/systems/discretization.hpp"    // IWYU pragma: export
#include "damp/systems/realization.hpp"       // IWYU pragma: export
#include "damp/systems/state_space.hpp"       // IWYU pragma: export
#include "damp/systems/transfer_function.hpp" // IWYU pragma: export
#include "damp/systems/zeros.hpp"             // IWYU pragma: export
#include "damp/systems/zpk.hpp"               // IWYU pragma: export

// --- Design / CST synthesis (init / offline) --------------------------------
#include "damp/design/linearization.hpp"   // IWYU pragma: export
#include "damp/design/lyapunov.hpp"        // IWYU pragma: export
#include "damp/design/minreal.hpp"         // IWYU pragma: export
#include "damp/design/model_inverse.hpp"   // IWYU pragma: export
#include "damp/design/model_reduction.hpp" // IWYU pragma: export
#include "damp/design/options.hpp"         // IWYU pragma: export
#include "damp/design/pid_design.hpp"      // IWYU pragma: export
#include "damp/design/pole_placement.hpp"  // IWYU pragma: export
#include "damp/design/qp.hpp"              // IWYU pragma: export
#include "damp/design/riccati.hpp"         // IWYU pragma: export
#include "damp/design/stability.hpp"       // IWYU pragma: export
#include "damp/design/status.hpp"          // IWYU pragma: export

// --- Core DiD controllers ---------------------------------------------------
#include "damp/controllers/action_governor.hpp"      // IWYU pragma: export
#include "damp/controllers/adrc.hpp"                 // IWYU pragma: export
#include "damp/controllers/composition.hpp"          // IWYU pragma: export
#include "damp/controllers/esc.hpp"                  // IWYU pragma: export
#include "damp/controllers/harmonic_suppression.hpp" // IWYU pragma: export
#include "damp/controllers/law_detail.hpp"           // IWYU pragma: export
#include "damp/controllers/lead_lag.hpp"             // IWYU pragma: export
#include "damp/controllers/lqg.hpp"                  // IWYU pragma: export
#include "damp/controllers/lqgi.hpp"                 // IWYU pragma: export
#include "damp/controllers/lqi.hpp"                  // IWYU pragma: export
#include "damp/controllers/lqr.hpp"                  // IWYU pragma: export
#include "damp/controllers/mpc.hpp"                  // IWYU pragma: export
#include "damp/controllers/oflqg.hpp"                // IWYU pragma: export
#include "damp/controllers/ofmpc.hpp"                // IWYU pragma: export
#include "damp/controllers/pid.hpp"                  // IWYU pragma: export
#include "damp/controllers/pr.hpp"                   // IWYU pragma: export
#include "damp/controllers/repetitive.hpp"           // IWYU pragma: export
#include "damp/controllers/smc.hpp"                  // IWYU pragma: export
#include "damp/controllers/smith_predictor.hpp"      // IWYU pragma: export
#include "damp/controllers/stsmc.hpp"                // IWYU pragma: export

// --- Core design-is-deploy estimators ---------------------------------------
#include "damp/estimation/dob.hpp"                   // IWYU pragma: export
#include "damp/estimation/ekf.hpp"                   // IWYU pragma: export
#include "damp/estimation/eskf.hpp"                  // IWYU pragma: export
#include "damp/estimation/eso.hpp"                   // IWYU pragma: export
#include "damp/estimation/ins_eskf.hpp"              // IWYU pragma: export
#include "damp/estimation/kalman.hpp"                // IWYU pragma: export
#include "damp/estimation/luenberger.hpp"            // IWYU pragma: export
#include "damp/estimation/mhe.hpp"                   // IWYU pragma: export
#include "damp/estimation/rls.hpp"                   // IWYU pragma: export
#include "damp/estimation/sensor_fusion.hpp"         // IWYU pragma: export
#include "damp/estimation/state_variable_filter.hpp" // IWYU pragma: export
#include "damp/estimation/ukf.hpp"                   // IWYU pragma: export

// --- Identification on target (excitation, experiments, grey-box fits) ------
#include "damp/estimation/excitation.hpp"           // IWYU pragma: export
#include "damp/estimation/frequency_response.hpp"   // IWYU pragma: export
#include "damp/estimation/parameter_estimation.hpp" // IWYU pragma: export
#include "damp/estimation/relay_autotune.hpp"       // IWYU pragma: export

// --- Filters, phase-locked loops + fixed-step ODE integrators ---------------
#include "damp/filters/filters.hpp"       // IWYU pragma: export
#include "damp/filters/pll.hpp"           // IWYU pragma: export
#include "damp/simulation/integrator.hpp" // IWYU pragma: export

// --- Frame math (3D attitude + three-phase Clarke/Park) ---------------------
#include "damp/math/geometry.hpp"   // IWYU pragma: export
#include "damp/math/transforms.hpp" // IWYU pragma: export

// --- Encoders: counts → exact position, smooth speed ------------------------
#include "damp/encoder/axis_position.hpp" // IWYU pragma: export
#include "damp/encoder/gate.hpp"          // IWYU pragma: export
#include "damp/encoder/quadrature.hpp"    // IWYU pragma: export
#include "damp/encoder/sincos.hpp"        // IWYU pragma: export
#include "damp/encoder/tracker.hpp"       // IWYU pragma: export
#include "damp/encoder/unwrap.hpp"        // IWYU pragma: export

// --- Embed helpers (controls-adjacent; pair with ETL for generic plumbing) --
#include "damp/toolbox/actuator.hpp"     // IWYU pragma: export
#include "damp/toolbox/bounds.hpp"       // IWYU pragma: export
#include "damp/toolbox/conditioning.hpp" // IWYU pragma: export
#include "damp/toolbox/iec61131.hpp"     // IWYU pragma: export
#include "damp/toolbox/io.hpp"           // IWYU pragma: export
#include "damp/toolbox/logic.hpp"        // IWYU pragma: export
#include "damp/toolbox/lookup.hpp"       // IWYU pragma: export
#include "damp/toolbox/scaling.hpp"      // IWYU pragma: export
#include "damp/toolbox/thermistor.hpp"   // IWYU pragma: export
#include "damp/toolbox/timing.hpp"       // IWYU pragma: export
