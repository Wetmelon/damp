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
 * This umbrella is the embeddable surface:
 * - linear algebra + pluggable math backend (via matrix / math includes)
 * - LTI systems (SS / TF / ZPK) + @c design:: synthesis
 * - core design-is-deploy controllers and estimators (including ESKF / attitude fusion)
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
 * constexpr auto art = design::lqi_bundle(sys_d, Q_aug, R);
 * static_assert(art.success);
 * constinit LQI controller{art.runtime.controller};
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
#include "damp/systems/state_space.hpp"       // IWYU pragma: export
#include "damp/systems/transfer_function.hpp" // IWYU pragma: export
#include "damp/systems/zpk.hpp"               // IWYU pragma: export

// --- Design / CST synthesis (init / offline) --------------------------------
#include "damp/design/linearization.hpp"   // IWYU pragma: export
#include "damp/design/lyapunov.hpp"        // IWYU pragma: export
#include "damp/design/minreal.hpp"         // IWYU pragma: export
#include "damp/design/model_reduction.hpp" // IWYU pragma: export
#include "damp/design/pid_design.hpp"      // IWYU pragma: export
#include "damp/design/pole_placement.hpp"  // IWYU pragma: export
#include "damp/design/qp.hpp"              // IWYU pragma: export
#include "damp/design/riccati.hpp"         // IWYU pragma: export
#include "damp/design/stability.hpp"       // IWYU pragma: export
#include "damp/design/synthesis.hpp"       // IWYU pragma: export

// --- Core DiD controllers ---------------------------------------------------
#include "damp/controllers/action_governor.hpp" // IWYU pragma: export
#include "damp/controllers/adrc.hpp"            // IWYU pragma: export
#include "damp/controllers/composition.hpp"     // IWYU pragma: export
#include "damp/controllers/lead_lag.hpp"        // IWYU pragma: export
#include "damp/controllers/lqg.hpp"             // IWYU pragma: export
#include "damp/controllers/lqgi.hpp"            // IWYU pragma: export
#include "damp/controllers/lqi.hpp"             // IWYU pragma: export
#include "damp/controllers/lqr.hpp"             // IWYU pragma: export
#include "damp/controllers/pid.hpp"             // IWYU pragma: export
#include "damp/controllers/pr.hpp"              // IWYU pragma: export
#include "damp/controllers/smc.hpp"             // IWYU pragma: export
#include "damp/controllers/smith_predictor.hpp" // IWYU pragma: export
#include "damp/controllers/stsmc.hpp"           // IWYU pragma: export

// --- Core design-is-deploy estimators ---------------------------------------
#include "damp/estimation/dob.hpp"           // IWYU pragma: export
#include "damp/estimation/ekf.hpp"           // IWYU pragma: export
#include "damp/estimation/eskf.hpp"          // IWYU pragma: export
#include "damp/estimation/kalman.hpp"        // IWYU pragma: export
#include "damp/estimation/luenberger.hpp"    // IWYU pragma: export
#include "damp/estimation/rls.hpp"           // IWYU pragma: export
#include "damp/estimation/sensor_fusion.hpp" // IWYU pragma: export
#include "damp/estimation/ukf.hpp"           // IWYU pragma: export

// --- Filters + fixed-step ODE integrators -----------------------------------
#include "damp/filters/differentiator.hpp" // IWYU pragma: export
#include "damp/filters/filters.hpp"        // IWYU pragma: export
#include "damp/simulation/integrator.hpp"  // IWYU pragma: export

// --- Frame math (3D attitude + three-phase Clarke/Park) ---------------------
#include "damp/math/geometry.hpp"   // IWYU pragma: export
#include "damp/math/transforms.hpp" // IWYU pragma: export

// --- Embed helpers (controls-adjacent; pair with ETL for generic plumbing) --
#include "damp/toolbox/actuator.hpp"     // IWYU pragma: export
#include "damp/toolbox/bounds.hpp"       // IWYU pragma: export
#include "damp/toolbox/conditioning.hpp" // IWYU pragma: export
#include "damp/toolbox/encoder.hpp"      // IWYU pragma: export
#include "damp/toolbox/iec61131.hpp"     // IWYU pragma: export
#include "damp/toolbox/io.hpp"           // IWYU pragma: export
#include "damp/toolbox/logic.hpp"        // IWYU pragma: export
#include "damp/toolbox/lookup.hpp"       // IWYU pragma: export
#include "damp/toolbox/scaling.hpp"      // IWYU pragma: export
#include "damp/toolbox/thermal.hpp"      // IWYU pragma: export
#include "damp/toolbox/timing.hpp"       // IWYU pragma: export
