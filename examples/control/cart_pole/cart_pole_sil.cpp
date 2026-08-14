// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file cart_pole_sil.cpp
 * @brief Cart-pole LQR — host nonlinear SIL (calls cart_pole_controller.hpp)
 *
 * Host only. Flashable path: cart_pole_sketch.cpp + cart_pole_controller.hpp.
 */

#include "cart_pole_controller.hpp"
#include "damp/backend.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/simulation/integrator.hpp"
#include "damp/simulation/plot_plotly.hpp"
#include "damp/simulation/simulate.hpp"
#include "damp/simulation/solver.hpp"
#include "fmt/base.h"
#include "fmt/core.h"

using namespace damp;
using namespace damp::examples_cart_pole;
using namespace damp::sim;

// Same design as sketch; double for host plant numerics
static constinit LQR<4, 1, double> controller{lqr_d};

/**
 * Nonlinear cart-pole, θ = 0 upright.
 * Jacobian at the origin matches cart_pole_controller.hpp (a22,a23,a42,a43,b2,b4).
 */
static ColVec<4, double> cart_pole_f(double /*t*/, const ColVec<4, double>& x, const ColVec<1, double>& u) {
    const double x_dot = x(1);
    const double theta = x(2);
    const double theta_dot = x(3);
    const double force = u(0);

    const auto [s, c] = damp::sincos(theta);
    // denom = M + m sin²θ  (→ M at upright)
    const double denom = M + (m_pole * s * s);

    // ẍ = (u − b ẋ + m sin(θ) (L θ̇² − g cos(θ))) / denom
    const double x_ddot = (force - (b_fric * x_dot) + (m_pole * s * ((L * theta_dot * theta_dot) - (g * c)))) / denom;
    // θ̈ = (g sin(θ) − ẍ cos(θ)) / L
    const double theta_ddot = ((g * s) - (x_ddot * c)) / L;

    return ColVec<4, double>{x_dot, x_ddot, theta_dot, theta_ddot};
}

int main() {
    fmt::print("===== Cart-pole LQR (host SIL) =====\n\n");
    fmt::print("M={:.1f} kg, m={:.1f} kg, L={:.1f} m, Ts={:.0f} Hz\n", M, m_pole, L, 1.0 / Ts);
    fmt::print(
        "K = [{:.4f}, {:.4f}, {:.4f}, {:.4f}]\n\n",
        controller.K(0, 0),
        controller.K(0, 1),
        controller.K(0, 2),
        controller.K(0, 3)
    );

    auto output = [](const ColVec<4, double>& x) -> ColVec<4, double> { return x; };
    auto ctrl = [](const ColVec<4, double>& x) -> ColVec<1, double> {
        return control_period(controller, x);
    };

    RK4<4>          rk4;
    FixedStepSolver solver(rk4, 0.001); // 1 ms integration; control re-evaluated each step

    // Small tip: 5 cm cart offset + ~3° pole
    ColVec<4, double> x0{0.05, 0.0, 0.05, 0.0};

    auto sim = simulate_state_feedback<4, 1, 4>(cart_pole_f, output, ctrl, solver, x0, {0.0, 8.0});

    fmt::print("Simulated {} samples\n", sim.t.size());
    fmt::print(
        "Final: x={:.4f} m, θ={:.4f} rad ({:.2f}°)\n",
        sim.x.back()(0),
        sim.x.back()(2),
        sim.x.back()(2) * 180.0 / damp::numbers::pi_v<double>
    );

    const double xf = damp::abs(sim.x.back()(0));
    const double thf = damp::abs(sim.x.back()(2));
    if (xf > 0.05 || thf > 0.05) {
        fmt::print("WARNING: final state not near origin (xf={:.4f}, |θ|={:.4f})\n", xf, thf);
    } else {
        fmt::print("Regulated to origin band (|x|,|θ| < 0.05)\n");
    }

    auto fig = plot::plot_simulation(sim, "Cart-pole LQR (nonlinear SIL)");
    plot::write_html(fig, "plots/control/cart_pole_lqr.html");
    fmt::print("\nPlot written to plots/control/cart_pole_lqr.html\n");

    return 0;
}
