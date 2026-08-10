// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>
#include <vector>

#include "damp/backend.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/simulation/integrator.hpp"
#include "damp/simulation/simulate.hpp"
#include "damp/simulation/solver.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using namespace damp::sim;

// Pendulum parameters for simulation tests
constexpr double g_pend = 9.81;
constexpr double L_pend = 1.0;
constexpr double m_pend = 1.0;
constexpr double b_pend = 0.1;

TEST_CASE("simulate_state_feedback - pendulum converges under LQR") {
    // Linearize at upright equilibrium
    Matrix<2, 2> A{{0.0, 1.0}, {g_pend / L_pend, -b_pend / (m_pend * L_pend * L_pend)}};
    Matrix<2, 1> B{{0.0}, {1.0 / (m_pend * L_pend * L_pend)}};

    auto Q = Matrix<2, 2>::identity() * 10.0;
    auto R = Matrix<1, 1>{{1.0}};
    auto Ts = 0.01;

    auto lqr_result = design::discrete_lqr_from_continuous(A, B, Q, R, Ts);
    CHECK(lqr_result.success);
    LQR<2, 1, double> controller{lqr_result};

    // Nonlinear plant
    auto plant = [](double /*t*/, const ColVec<2>& x, const ColVec<1>& u) -> ColVec<2> {
        double theta = x(0, 0);
        double theta_dot = x(1, 0);
        double torque = u(0, 0);

        double theta_ddot = ((g_pend / L_pend) * damp::sin(theta))
                          - ((b_pend / (m_pend * L_pend * L_pend)) * theta_dot)
                          + ((1.0 / (m_pend * L_pend * L_pend)) * torque);

        return ColVec<2>{theta_dot, theta_ddot};
    };

    auto output = [](const ColVec<2>& x) -> ColVec<2> { return x; };
    auto ctrl = [&](const ColVec<2>& x) -> ColVec<1> { return controller.control(x); };

    RK4<2>          rk4;
    FixedStepSolver solver(rk4, 0.001);

    // Start at 15 degrees
    ColVec<2> x0{0.2618, 0.0};
    auto      sim = simulate_state_feedback<2, 1, 2>(plant, output, ctrl, solver, x0, {0.0, 5.0});

    // Should converge to near zero by 5 seconds
    CHECK(sim.x.back()(0, 0) == doctest::Approx(0.0).epsilon(0.01));
    CHECK(sim.x.back()(1, 0) == doctest::Approx(0.0).epsilon(0.01));

    // Check that we recorded everything
    CHECK(sim.t.size() == sim.x.size());
    CHECK(sim.t.size() == sim.y.size());
    CHECK(sim.t.size() == sim.u.size());
}

TEST_CASE("simulate_lti - matches step_response for zero-input") {
    // Simple 1st-order system: dx/dt = -x, y = x
    Matrix<1, 1> A{{-1.0}};
    Matrix<1, 1> B{{1.0}};
    Matrix<1, 1> C{{1.0}};
    StateSpace   sys{.A = A, .B = B, .C = C, .Ts = 0.0};

    // Zero control (free response)
    auto controller = [](const ColVec<1>& /*y*/) -> ColVec<1> {
        return ColVec<1>{0.0};
    };

    RK4<1>          rk4;
    FixedStepSolver solver(rk4, 0.01);

    ColVec<1> x0{1.0};
    auto      sim = simulate_lti(sys, controller, solver, x0, {0.0, 3.0});

    // x(t) = exp(-t), so x(3) ≈ 0.0498
    CHECK(sim.x.back()(0, 0) == doctest::Approx(std::exp(-3.0)).epsilon(1e-6));
}

TEST_CASE("SIL: PIDController tracks a step on a first-order plant (simulate_lti_siso)") {
    // Plant G(s) = 1/(s+1):  ẋ = -x + u, y = x.  PI designed for ~unit DC gain tracking.
    constexpr double Ts = 0.01;
    StateSpace       sys{
              .A = Matrix<1, 1>{{-1.0}},
              .B = Matrix<1, 1>{{1.0}},
              .C = Matrix<1, 1>{{1.0}},
              .Ts = 0.0,
    };

    // design::… → discretize → controller (same shape as on target; double for SIL).
    auto pid = PIController<double>{design::pid(4.0, 8.0, 0.0).discretize(Ts)};

    RK4<1>          rk4;
    FixedStepSolver fine{rk4, Ts / 10.0}; // 10× plant substeps per control tick

    ColVec<1>    x0{0.0};
    const double r = 1.0;
    auto         sil = simulate_lti_siso(sys, pid, r, fine, Ts, x0, {0.0, 3.0});

    CHECK(sil.t.size() > 10);
    CHECK(sil.y.back()[0] == doctest::Approx(r).epsilon(0.05)); // settled near the step
    CHECK(sil.u.back()[0] == doctest::Approx(r).epsilon(0.15)); // DC: u_ss = r for this plant
}

TEST_CASE("simulate_discrete - basic feedback") {
    // Discrete integrator: x_{k+1} = x_k + Ts * u_k
    double       Ts = 0.1;
    Matrix<1, 1> A_d{{1.0}};
    Matrix<1, 1> B_d{{Ts}};
    Matrix<1, 1> C_d{{1.0}};
    StateSpace   sys{.A = A_d, .B = B_d, .C = C_d, .Ts = Ts};

    // Simple proportional controller: u = -10 * y (drives x to 0)
    auto controller = [](const ColVec<1>& y) -> ColVec<1> {
        return ColVec<1>{-10.0 * y(0, 0)};
    };

    ColVec<1> x0{1.0};
    auto      sim = simulate_discrete(sys, controller, x0, 100);

    CHECK(sim.t.size() == 101);
    CHECK(sim.t.front() == 0.0);
    CHECK(sim.t.back() == doctest::Approx(10.0).epsilon(1e-10));

    // State should converge to zero
    CHECK(sim.x.back()(0, 0) == doctest::Approx(0.0).epsilon(1e-6));
}

TEST_CASE("simulate - output feedback with D matrix") {
    // Controller: u = 0 (open loop, step input via plant)
    auto plant = [](double /*t*/, const ColVec<1>& x, const ColVec<1>& u) -> ColVec<1> {
        return ColVec<1>{(-2.0 * x(0, 0)) + u(0, 0)};
    };
    auto output = [](const ColVec<1>& x) -> ColVec<1> { return x; };
    auto ctrl = [](const ColVec<1>& /*y*/) -> ColVec<1> { return ColVec<1>{1.0}; }; // Constant input

    RK4<1>          rk4;
    FixedStepSolver solver(rk4, 0.01);

    ColVec<1> x0{0.0};
    auto      sim = simulate<1, 1, 1>(plant, output, ctrl, solver, x0, {0.0, 5.0});

    // Steady state: dx/dt=0 => x = u/2 = 0.5
    CHECK(sim.x.back()(0, 0) == doctest::Approx(0.5).epsilon(0.01));
}

TEST_CASE("simulate_discrete_nonlinear - logistic-style nonlinear map") {
    // Natively-discrete nonlinear plant with no continuous-time analogue:
    //   x[k+1] = x[k] + Ts*(-x[k]^3 + u[k])
    // A cubic restoring term — not expressible as Ax+Bu, so this must use the
    // discrete nonlinear path, not simulate_discrete (linear) or the ODE solver.
    const double Ts = 0.01;

    auto plant = [&](std::size_t /*k*/, const ColVec<1>& x, const ColVec<1>& u) -> ColVec<1> {
        const double xv = x(0, 0);
        return ColVec<1>{xv + (Ts * (-(xv * xv * xv) + u(0, 0)))};
    };
    auto output = [](const ColVec<1>& x) -> ColVec<1> { return x; };
    // Constant input u = 8 → equilibrium where x^3 = 8 → x = 2.
    auto ctrl = [](const ColVec<1>& /*y*/) -> ColVec<1> { return ColVec<1>{8.0}; };

    auto sim = simulate_discrete_nonlinear<1, 1, 1>(plant, output, ctrl, ColVec<1>{0.0}, Ts, 5000);

    CHECK(sim.t.size() == 5001);
    CHECK(sim.t.front() == 0.0);
    CHECK(sim.t.back() == doctest::Approx(50.0).epsilon(1e-9));
    CHECK(sim.x.back()(0, 0) == doctest::Approx(2.0).epsilon(1e-3)); // settles at cube root of 8
}

TEST_CASE("simulate_discrete_nonlinear - matches simulate_discrete on a linear plant") {
    // When f(k,x,u) = Ax + Bu, the nonlinear path must reproduce the linear one.
    const double Ts = 0.1;
    Matrix<1, 1> A{{1.0}};
    Matrix<1, 1> B{{Ts}};
    Matrix<1, 1> C{{1.0}};
    StateSpace   sys{.A = A, .B = B, .C = C, .Ts = Ts};

    auto ctrl = [](const ColVec<1>& y) -> ColVec<1> { return ColVec<1>{-10.0 * y(0, 0)}; };

    auto lin = simulate_discrete(sys, ctrl, ColVec<1>{1.0}, 100);

    auto plant = [&](std::size_t, const ColVec<1>& x, const ColVec<1>& u) -> ColVec<1> {
        return A * x + B * u;
    };
    auto output = [&](const ColVec<1>& x) -> ColVec<1> { return C * x; };
    auto nl = simulate_discrete_nonlinear<1, 1, 1>(plant, output, ctrl, ColVec<1>{1.0}, Ts, 100);

    REQUIRE(lin.x.size() == nl.x.size());
    for (std::size_t k = 0; k < lin.x.size(); ++k) {
        CHECK(nl.x[k](0, 0) == doctest::Approx(lin.x[k](0, 0)));
        CHECK(nl.t[k] == doctest::Approx(lin.t[k]));
    }
}

TEST_CASE("SimulationResult has consistent sizes") {
    auto plant = [](double, const ColVec<1>& x, const ColVec<1>& u) -> ColVec<1> { return -x + u; };
    auto output = [](const ColVec<1>& x) -> ColVec<1> { return x; };
    auto ctrl = [](const ColVec<1>&) -> ColVec<1> { return ColVec<1>{0.0}; };

    RK4<1>          rk4;
    FixedStepSolver solver(rk4, 0.1);

    auto sim = simulate<1, 1, 1>(plant, output, ctrl, solver, ColVec<1>{1.0}, {0.0, 1.0});

    CHECK(sim.t.size() == sim.x.size());
    CHECK(sim.t.size() == sim.y.size());
    CHECK(sim.t.size() == sim.u.size());
    CHECK(sim.t.front() == 0.0);
    CHECK(sim.t.back() == doctest::Approx(1.0).epsilon(0.01));
}

TEST_CASE("simulate_discrete - biproper plant includes D*u in y") {
    // Pure feedthrough biproper: A=B=C=0, D=2 → y = 2u for every sample.
    const double Ts = 0.1;
    StateSpace   sys{
          .A = Matrix<1, 1>{{0.0}},
          .B = Matrix<1, 1>{{0.0}},
          .C = Matrix<1, 1>{{0.0}},
          .D = Matrix<1, 1>{{2.0}},
          .Ts = Ts,
    };
    auto ctrl = [](const ColVec<1>& /*y*/) -> ColVec<1> { return ColVec<1>{1.5}; };

    auto sim = simulate_discrete(sys, ctrl, ColVec<1>{0.0}, 5);
    REQUIRE(sim.y.size() == 6);
    for (const auto& y : sim.y) {
        CHECK(y(0, 0) == doctest::Approx(3.0)); // 2 * 1.5
    }
    CHECK(sim.u.back()(0, 0) == doctest::Approx(1.5));
}

TEST_CASE("simulate_lti - biproper plant includes D*u in y") {
    // Continuous biproper: ẋ = -x + u, y = x + 0.5·u (D ≠ 0).
    StateSpace sys{
        .A = Matrix<1, 1>{{-1.0}},
        .B = Matrix<1, 1>{{1.0}},
        .C = Matrix<1, 1>{{1.0}},
        .D = Matrix<1, 1>{{0.5}},
        .Ts = 0.0,
    };
    auto ctrl = [](const ColVec<1>& /*y*/) -> ColVec<1> { return ColVec<1>{2.0}; };

    RK4<1>          rk4;
    FixedStepSolver solver(rk4, 0.01);
    auto            sim = simulate_lti(sys, ctrl, solver, ColVec<1>{0.0}, {0.0, 5.0});

    // Steady state: x_ss = u = 2, y_ss = x + 0.5·u = 3
    CHECK(sim.x.back()(0, 0) == doctest::Approx(2.0).epsilon(0.01));
    CHECK(sim.y.back()(0, 0) == doctest::Approx(3.0).epsilon(0.01));
    // Initial sample after control update: x=0, u=2 → y = 0 + 0.5*2 = 1
    CHECK(sim.y.front()(0, 0) == doctest::Approx(1.0).epsilon(1e-12));
}
