// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file pendulum_sil.cpp
 * @brief Upright pendulum LQR — host nonlinear SIL (calls pendulum_controller.hpp)
 *
 * Host only. Flashable path: pendulum_sketch.cpp + pendulum_controller.hpp.
 * Also prints gains and a host-only redesign about another linearization point
 * (folded from the former example_lqr_pendulum gain-print demo).
 */

#include "damp/backend.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/simulation/integrator.hpp"
#include "damp/simulation/plot_plotly.hpp"
#include "damp/simulation/simulate.hpp"
#include "damp/simulation/solver.hpp"
#include "fmt/base.h"
#include "fmt/core.h"
#include "pendulum_controller.hpp"

using namespace damp;
using namespace damp::examples_pendulum;
using namespace damp::sim;

// Same design as sketch; double for host plant numerics
static constinit LQR<2, 1, double> controller{lqr_d};

/**
 * Nonlinear pendulum: θ̈ = (g/L)sinθ − (b/mL²)θ̇ + (1/mL²)u
 * Jacobian at the origin matches pendulum_controller.hpp linearization.
 */
static ColVec<2, double> pendulum_f(double /*t*/, const ColVec<2, double>& x, const ColVec<1, double>& u) {
    const double theta = x(0);
    const double theta_dot = x(1);
    const double torque = u(0);

    const double theta_ddot = ((g / L) * damp::sin(theta)) - ((b_damp / (m * L * L)) * theta_dot)
                            + ((1.0 / (m * L * L)) * torque);

    return ColVec<2, double>{theta_dot, theta_ddot};
}

int main() {
    fmt::print("===== Pendulum LQR (host SIL) =====\n\n");
    fmt::print("g={:.2f}, L={:.1f} m, m={:.1f} kg, b={:.1f}, Ts={:.0f} Hz\n", g, L, m, b_damp, 1.0 / Ts);
    fmt::print("K = [{:.4f}, {:.4f}]\n\n", controller.K(0, 0), controller.K(0, 1));

    // One-tick smoke (same float deploy path as sketch / former lqr print)
    {
        ColVec<2, float> x_smoke{0.1f, 0.0f};
        LQR<2, 1, float> ctrl_f{lqr_d.as<float>()};
        const float      u_smoke = control_period(ctrl_f, x_smoke)(0);
        fmt::print("u at θ=0.1 rad (float): {:.4f}\n\n", u_smoke);
    }

    auto output = [](const ColVec<2, double>& x) -> ColVec<2, double> { return x; };
    auto ctrl = [](const ColVec<2, double>& x) -> ColVec<1, double> {
        return control_period(controller, x);
    };

    RK4<2>          rk4;
    FixedStepSolver solver(rk4, 0.001);

    ColVec<2, double> x0{0.5236, 0.0}; // ~30° tip

    auto sim = simulate_state_feedback<2, 1, 2>(pendulum_f, output, ctrl, solver, x0, {0.0, 5.0});

    fmt::print("Simulated {} time steps\n", sim.t.size());
    fmt::print(
        "Final state: theta = {:.4f} rad, theta_dot = {:.4f} rad/s\n",
        sim.x.back()(0),
        sim.x.back()(1)
    );

    auto fig = plot::plot_simulation(sim, "Pendulum LQR Simulation");
    plot::write_html(fig, "plots/control/pendulum_sim.html");
    fmt::print("\nPlot written to plots/control/pendulum_sim.html\n");

    // Runtime redesign (host only) — still gates on success
    fmt::print("\n===== Runtime LQR Redesign (Q = 50 I) =====\n");
    auto Q2 = Matrix<2, 2>::identity() * 50.0;
    auto lqr_d2 = design::discrete_lqr_from_continuous(sys.A, sys.B, Q2, R, Ts);
    if (!lqr_d2.success) {
        fmt::print("LQR redesign failed\n");
        return 1;
    }
    LQR<2, 1, double> controller2{lqr_d2};
    fmt::print("New LQR Gain K = [{:.4f}, {:.4f}]\n", lqr_d2.K(0, 0), lqr_d2.K(0, 1));

    auto ctrl2 = [&](const ColVec<2, double>& x) -> ColVec<1, double> {
        return controller2.control(x);
    };

    auto sim2 = simulate_state_feedback<2, 1, 2>(pendulum_f, output, ctrl2, solver, x0, {0.0, 5.0});
    fmt::print(
        "Final state: theta = {:.4f} rad, theta_dot = {:.4f} rad/s\n",
        sim2.x.back()(0),
        sim2.x.back()(1)
    );

    auto fig2 = plot::plot_simulation(sim2, "Pendulum LQR (Higher Q) Simulation");
    plot::write_html(fig2, "plots/control/pendulum_sim_high_q.html");
    fmt::print("Plot written to plots/control/pendulum_sim_high_q.html\n");

    {
        std::vector<double> th1, thd1, th2, thd2;
        th1.reserve(sim.x.size());
        thd1.reserve(sim.x.size());
        th2.reserve(sim2.x.size());
        thd2.reserve(sim2.x.size());
        for (const auto& x : sim.x) {
            th1.push_back(x(0));
            thd1.push_back(x(1));
        }
        for (const auto& x : sim2.x) {
            th2.push_back(x(0));
            thd2.push_back(x(1));
        }
        auto fig_phase = plot::plot_xy(
            {
                {.name = "Q = 10 I", .x = th1, .y = thd1, .lines = true, .markers = false, .color = "#1f77b4"},
                {.name = "Q = 50 I", .x = th2, .y = thd2, .lines = true, .markers = false, .color = "#ff7f0e"},
            },
            "Pendulum phase portrait (θ vs θ̇)",
            "θ [rad]",
            "θ̇ [rad/s]",
            /*equal_aspect=*/false
        );
        plot::write_html(fig_phase, "plots/control/pendulum_phase.html");
        fmt::print("Plot written to plots/control/pendulum_phase.html\n");
    }

    // Host-only: redesign about another angle (not the flash path)
    fmt::print("\nHost-only redesign at θ=0.5 rad:\n");
    const double theta_op = 0.5;
    const double a21 = (g / L) * damp::cos(theta_op);
    const double a22 = -b_damp / (m * L * L);
    Matrix<2, 2> A_op{{0.0, 1.0}, {a21, a22}};
    Matrix<2, 1> B_op{{0.0}, {1.0 / (m * L * L)}};
    const auto   res_op = design::lqrd(A_op, B_op, Q, R, Ts);
    if (res_op.success) {
        fmt::print("  K_op = [{:.4f}, {:.4f}]\n", res_op.K(0, 0), res_op.K(0, 1));
    }

    return 0;
}
