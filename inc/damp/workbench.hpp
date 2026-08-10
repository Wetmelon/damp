// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file workbench.hpp
 * @brief Host product surface — CST analysis, ODE solvers, closed-loop sim.
 *
 * Includes the slim embeddable @ref control.hpp plus workstation tools:
 * frequency-domain analysis (Bode, Nyquist, margins, sweeps), ODE solvers on
 * top of the integrators in @c simulation/integrator.hpp, closed-loop
 * @c simulate helpers, and MATLAB®-style aliases.
 *
 * Host-only — not for target firmware. These extras allocate on the heap
 * (`std::vector`). Deploy with @ref control.hpp; keep analysis / solvers /
 * SIL here on the host. @c make embedded-check enforces the control.hpp contract.
 *
 * @note The Plotly backend (`damp/simulation/plot_plotly.hpp`) is intentionally
 *       not pulled in here: it requires plotlypp and nlohmann-json on the
 *       include path and must never leak into the embeddable umbrella. Include
 *       it directly when you need HTML figures.
 *
 * @code
 * #include "damp/workbench.hpp"          // host: design + analyze + solve + simulate
 * using namespace damp;
 *
 * const auto omega = analysis::logspace(1.0, 1000.0, 200);
 * const auto metrics = analysis::loop_metrics(loop, omega);
 * auto traj = fixed_solve\<RK4\>(f, x0, {0.0, 1.0}, 1e-3);
 * @endcode
 *
 * @see control.hpp for the slim embeddable core
 */

#include "damp/control.hpp" // IWYU pragma: export

// --- Host CST analysis + ODE solvers + closed-loop sim (heap OK) ------------
#include "damp/analysis/analysis.hpp"           // IWYU pragma: export  (Bode/Nyquist/margins/sweeps)
#include "damp/controllers/mpc.hpp"             // IWYU pragma: export  (constrained MPC; allocation-free — include directly on-target)
#include "damp/controllers/offset_free_mpc.hpp" // IWYU pragma: export  (MPC + disturbance-augmented Kalman bundle)
#include "damp/estimation/mhe.hpp"              // IWYU pragma: export  (moving-horizon estimation; allocation-free — include on-target)
#include "damp/matlab.hpp"                      // IWYU pragma: export  (MATLAB®-style aliases)
#include "damp/simulation/cached_zoh.hpp"       // IWYU pragma: export  (cached Ad/Bd Exact ZOH steps)
#include "damp/simulation/hybrid.hpp"           // IWYU pragma: export  (piecewise-LTI Exact event sim)
#include "damp/simulation/multirate.hpp"        // IWYU pragma: export  (multi-rate cascade SIL harness)
#include "damp/simulation/npy_export.hpp"       // IWYU pragma: export  (dense trace → .npy)
#include "damp/simulation/simulate.hpp"         // IWYU pragma: export  (closed-loop simulation)
#include "damp/simulation/solver.hpp"           // IWYU pragma: export  (fixed/adaptive ODE solvers)
