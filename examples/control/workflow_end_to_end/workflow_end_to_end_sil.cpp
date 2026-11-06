// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file workflow_end_to_end_sil.cpp
 * @brief Workflow end-to-end host SIL (calls workflow_end_to_end_controller.hpp)
 *
 * =============================================================================
 * 1. Design → Deploy  (design results + float runtime; gate on success)
 * =============================================================================
 * 2. Host SIL         (analysis + nonlinear plant calling the same runtime.step)
 * =============================================================================
 */

#include <cmath>
#include <numbers>

#include "damp/analysis/analysis.hpp"
#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/controllers/pr.hpp"
#include "damp/design/linearization.hpp"
#include "damp/design/pid_design.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/simulation/integrator.hpp"
#include "damp/simulation/simulate.hpp"
#include "damp/simulation/solver.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "fmt/base.h"
#include "fmt/core.h"
#include "workflow_end_to_end_controller.hpp"

using namespace damp;
using namespace damp::sim;
using namespace damp::examples_workflow_e2e;

// Design + host SIL (plant / kTs / ops from workflow_end_to_end_controller.hpp)

int main() {
    const auto lin = design::linearize<2, 1, 1>(plant_nonlinear, plant_output, kXop, kUop);

    StateSpace<2, 1, 1, double, 2, 1> sys_c{
        .A = lin.A,
        .B = lin.B,
        .C = lin.C,
        .D = lin.D,
        .G = Matrix<2, 2>::identity(),
        .H = Matrix<1, 1>::identity(),
        .Ts = 0.0
    };

    const auto sys_d = *discretize(sys_c, kTs, DiscretizationMethod::ZOH);

    // PID seed (side path): 20 rad/s crossover, phase margin from a 10% overshoot target.
    const double zeta_os = design::damping_ratio_from_overshoot_percent(10.0);
    const auto   pid_seed = design::pidtune(sys_c, {.wc = 20.0, .phase_margin_deg = design::phase_margin_from_damping_ratio(zeta_os)});

    const auto         Q = Matrix<2, 2>::diagonal({20.0, 2.0});
    const double       Q_integral = 80.0;
    const double       R = 0.25;
    const Matrix<2, 2> Q_kf{{1e-3, 0.0}, {0.0, 1e-2}};
    const double       R_kf = 5e-3;

    const auto lqgi_design = design::lqgi(sys_d, {.Q = Q, .R = R, .Q_integral = Q_integral}, {.process = Q_kf, .measurement = R_kf});
    if (!lqgi_design.success) {
        fmt::print("LQGI synthesis failed — aborting (do not run runtime on a bad design)\n");
        return 1;
    }

    const auto pr_design = design::pr({.Kp = 0.0, .Ki = 10.0, .w0 = 2.0 * std::numbers::pi, .wc = 6.0}, kTs);
    if (!pr_design.success) {
        fmt::print("PR design failed — aborting\n");
        return 1;
    }

    // Deploy-shaped float runtime (same objects a target would constinit after design)
    LQGI                runtime{lqgi_design.as<float>()};
    PRController<float> pr_runtime(pr_design.as<float>());

    // Frequency-domain analysis on the discrete plant (host)
    const auto omega = analysis::geomspace(1.0, 1000.0, 150);
    const auto bode_open = analysis::bode(sys_d, omega);
    const auto nyq_open = analysis::nyquist(sys_d, omega);

    // Nonlinear closed-loop: SIL calls the same float LQGI + PR path
    RK4<2>          rk4;
    FixedStepSolver solver(rk4, kTs);

    const float reference = 0.5f;
    auto        controller = [&](const ColVec<1>& y) -> ColVec<1> {
        const auto  y_f = static_cast<float>(y(0, 0));
        const float u_lqgi = runtime.step(ColVec<1, float>{reference}, ColVec<1, float>{y_f})(0, 0);
        const float u_pr = pr_runtime.control(reference - y_f);
        return ColVec<1>{static_cast<double>(u_lqgi + u_pr)};
    };

    // The plant is the same Dynamics object linearize took; simulate's output map is y(x), and
    // plant_output has no feedthrough.
    const auto      sensor = [](const ColVec<2>& x) { return plant_output(x, ColVec<1>{}); };
    const ColVec<2> x0{0.25, 0.0};
    const auto      sim = simulate<2, 1, 1>(plant_nonlinear, sensor, controller, solver, x0, {0.0, 2.0});

    fmt::print("Workflow end-to-end summary\n");
    fmt::print("  Linearized A(1,0): {:.6f}\n", lin.A(1, 0));
    fmt::print("  PID seed (pidtune, PID): success={}, Kp={:.6f}, Ki={:.6f}, Kd={:.6f}\n", pid_seed.success, pid_seed.Kp, pid_seed.Ki, pid_seed.Kd);
    fmt::print("  LQGI synthesis success: true\n");
    fmt::print("  Bode points: {}, Nyquist points: {}\n", bode_open.points.size(), nyq_open.points.size());
    fmt::print("  Final output: {:.6f}\n", sim.y.back()(0, 0));
    fmt::print("  Final control: {:.6f}\n", sim.u.back()(0, 0));

    return 0;
}
