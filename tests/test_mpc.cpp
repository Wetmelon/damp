// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>
#include <limits>

#include "damp/controllers/lqr.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/design/qp.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

constexpr double Ts = 0.1;

// Double integrator, ZOH-discretized at Ts = 0.1, position output.
constexpr StateSpace<2, 1, 1> double_integrator{
    .A = {{1.0, Ts}, {0.0, 1.0}},
    .B = {{0.5 * Ts * Ts}, {Ts}},
    .C = {{1.0, 0.0}},
    .Ts = Ts,
};

// Same plant with full-state output (C = I), for the LQR-equivalence test.
constexpr StateSpace<2, 1, 2> double_integrator_full{
    .A = {{1.0, Ts}, {0.0, 1.0}},
    .B = {{0.5 * Ts * Ts}, {Ts}},
    .C = {{1.0, 0.0}, {0.0, 1.0}},
    .Ts = Ts,
};

constexpr ColVec<2> plant_step(const ColVec<2>& x, double u) {
    return ColVec<2>{
        x(0) + (Ts * x(1)) + (0.5 * Ts * Ts * u),
        x(1) + (Ts * u),
    };
}

} // namespace

TEST_CASE("unconstrained MPC with a DARE terminal weight reproduces the LQR law") {
    // Finite-horizon LQ with terminal cost S (the DARE solution) gives
    // u0 = -K x for any horizon; the MPC must match move for move.
    const Matrix<2, 2> Q{{1.0, 0.0}, {0.0, 0.1}};
    const Matrix<1, 1> R{{0.5}};

    const auto lqr = design::discrete_lqr(double_integrator_full.A, double_integrator_full.B, Q, R);
    REQUIRE(lqr.success);

    design::MPCWeights<1, 2> weights{};
    weights.Qy = Q;
    weights.Qy_terminal = lqr.S;
    weights.Rdu = Matrix<1, 1>{}; // no move penalty: pure LQ cost
    weights.Ru = R;

    const auto art = design::state_mpc<8>(double_integrator_full, weights);
    REQUIRE(art.success);

    MPC             controller{art};
    ColVec<2>       x{1.0, -0.5};
    const ColVec<2> r{}; // regulation

    for (size_t k = 0; k < 30; ++k) {
        const auto   u_mpc = controller.control(r, x);
        const double u_lqr = -((lqr.K(0, 0) * x(0)) + (lqr.K(0, 1) * x(1)));
        CHECK(u_mpc(0) == doctest::Approx(u_lqr).epsilon(1e-6));
        CHECK(controller.last_status() == design::QPStatus::Success);
        CHECK(controller.last_slack() == doctest::Approx(0.0));
        x = plant_step(x, u_mpc(0));
    }
}

TEST_CASE("constrained MPC tracks a step within u and du limits") {
    design::MPCWeights<1, 1> weights{};
    weights.Rdu = Matrix<1, 1>{{0.05}};

    design::MPCConstraints<1, 1> limits{};
    limits.u_min = ColVec<1>{-0.5};
    limits.u_max = ColVec<1>{0.5};
    limits.du_min = ColVec<1>{-0.2};
    limits.du_max = ColVec<1>{0.2};

    const auto art = design::state_mpc<20, 8>(double_integrator, weights, limits);
    REQUIRE(art.success);

    MPC             controller{art};
    ColVec<2>       x{};
    const ColVec<1> r{1.0};

    double u_last = 0.0;
    double u_peak = 0.0;
    for (size_t k = 0; k < 200; ++k) {
        const auto u = controller.control(r, x);
        CHECK(controller.last_status() == design::QPStatus::Success);
        CHECK(controller.last_iterations() <= art.max_qp_iterations);

        // Hard limits hold at every tick.
        CHECK(u(0) <= 0.5 + 1e-9);
        CHECK(u(0) >= -0.5 - 1e-9);
        CHECK(u(0) - u_last <= 0.2 + 1e-9);
        CHECK(u(0) - u_last >= -0.2 - 1e-9);

        u_peak = (damp::abs(u(0)) > u_peak) ? damp::abs(u(0)) : u_peak;
        u_last = u(0);
        x = plant_step(x, u(0));
    }

    // The step was aggressive enough that the input limit actually bound...
    CHECK(u_peak == doctest::Approx(0.5));
    // ...and the loop still converged to the reference.
    CHECK(x(0) == doctest::Approx(1.0).epsilon(0.02));
    CHECK(x(1) == doctest::Approx(0.0).epsilon(0.02));
}

TEST_CASE("hard output constraint caps the predicted response") {
    design::MPCWeights<1, 1> weights{};
    weights.Rdu = Matrix<1, 1>{{0.01}};

    design::MPCConstraints<1, 1> limits{};
    limits.u_min = ColVec<1>{-1.0};
    limits.u_max = ColVec<1>{1.0};
    limits.y_max = ColVec<1>{1.02};    // overshoot limiter on the position output
    limits.y_max_ecr = ColVec<1>{0.0}; // hard: no slack allowed

    const auto art = design::state_mpc<20, 8>(double_integrator, weights, limits);
    REQUIRE(art.success);

    MPC             controller{art};
    ColVec<2>       x{};
    const ColVec<1> r{1.0};

    for (size_t k = 0; k < 300; ++k) {
        const auto u = controller.control(r, x);
        x = plant_step(x, u(0));
        CHECK(x(0) <= 1.02 + 1e-6); // plant == model, so the predicted cap holds exactly
    }
    CHECK(x(0) == doctest::Approx(1.0).epsilon(0.02));
}

TEST_CASE("soft output constraint recovers from an infeasible start") {
    // Start above the output ceiling with limited authority: no admissible
    // input sequence satisfies the hard rows, so a hard configuration reports
    // Infeasible while the default (soft) configuration keeps solving.
    design::MPCWeights<1, 1> weights{};
    weights.Rdu = Matrix<1, 1>{{0.01}};

    design::MPCConstraints<1, 1> limits{};
    limits.u_min = ColVec<1>{-1.0};
    limits.u_max = ColVec<1>{1.0};
    limits.y_max = ColVec<1>{1.0};

    const ColVec<2> x0{2.0, 0.0}; // y(0) = 2 > y_max
    const ColVec<1> r{0.5};

    // Hard variant: first solve is infeasible and the fallback holds the input.
    {
        design::MPCConstraints<1, 1> hard = limits;
        hard.y_max_ecr = ColVec<1>{0.0};
        const auto art = design::state_mpc<20, 8>(double_integrator, weights, hard);
        REQUIRE(art.success);
        MPC        controller{art};
        const auto u = controller.control(r, x0);
        CHECK(controller.last_status() == design::QPStatus::Infeasible);
        CHECK(u(0) == doctest::Approx(0.0)); // Δu falls back to zero
    }

    // Soft variant (default ECR = 1): solvable from the violated start, slack
    // active initially, and the loop still converges to the reference.
    {
        const auto art = design::state_mpc<20, 8>(double_integrator, weights, limits);
        REQUIRE(art.success);
        MPC controller{art};

        ColVec<2>  x = x0;
        const auto u0 = controller.control(r, x);
        CHECK(controller.last_status() == design::QPStatus::Success);
        CHECK(controller.last_slack() > 0.0);
        x = plant_step(x, u0(0));

        for (size_t k = 0; k < 400; ++k) {
            const auto u = controller.control(r, x);
            CHECK(controller.last_status() == design::QPStatus::Success);
            x = plant_step(x, u(0));
        }
        CHECK(x(0) == doctest::Approx(0.5).epsilon(0.02));
        CHECK(controller.last_slack() == doctest::Approx(0.0).epsilon(1e-6));
    }
}

TEST_CASE("measured-disturbance feedforward cancels a known load") {
    // The disturbance enters like the input (Bd = B). With the MD channel the
    // controller sees it in the prediction and cancels it exactly at steady
    // state; without it, the same load leaves a residual error.
    const Matrix<2, 1> Bd = double_integrator.B;

    design::MPCWeights<1, 1> weights{};
    weights.Rdu = Matrix<1, 1>{{0.05}};

    design::MPCConstraints<1, 1> limits{};
    limits.u_min = ColVec<1>{-2.0};
    limits.u_max = ColVec<1>{2.0};

    const ColVec<1> r{0.5};
    const double    d_load = 0.4;

    // With MD feedforward:
    const auto art_md = design::state_mpc<15, 5>(double_integrator, Bd, weights, limits);
    REQUIRE(art_md.success);
    MPC       md_controller{art_md};
    ColVec<2> x_md{};
    double    u_md = 0.0;
    for (size_t k = 0; k < 400; ++k) {
        const auto u = md_controller.control(r, x_md, ColVec<1>{d_load});
        u_md = u(0);
        x_md = plant_step(x_md, u_md + d_load);
    }

    // Baseline without MD, same load applied to the plant:
    const auto art_base = design::state_mpc<15, 5>(double_integrator, weights, limits);
    REQUIRE(art_base.success);
    MPC       base_controller{art_base};
    ColVec<2> x_base{};
    for (size_t k = 0; k < 400; ++k) {
        const auto u = base_controller.control(r, x_base);
        x_base = plant_step(x_base, u(0) + d_load);
    }

    // MD version: exact tracking, input settles at -d (net force zero).
    CHECK(x_md(0) == doctest::Approx(0.5).epsilon(1e-5));
    CHECK(u_md == doctest::Approx(-d_load).epsilon(1e-4));
    // And it beats the disturbance-blind baseline.
    CHECK(damp::abs(x_md(0) - 0.5) < damp::abs(x_base(0) - 0.5));
}

TEST_CASE("input target pulls the input when tracking is indifferent") {
    design::MPCWeights<1, 1> weights{};
    weights.Qy = Matrix<1, 1>{};          // no output cost:
    weights.Qy_terminal = Matrix<1, 1>{}; // pure input shaping
    weights.Ru = Matrix<1, 1>{{1.0}};
    weights.Rdu = Matrix<1, 1>{{0.1}};
    weights.u_target = ColVec<1>{0.3};

    const auto art = design::state_mpc<5>(double_integrator, weights);
    REQUIRE(art.success);

    MPC             controller{art};
    ColVec<2>       x{};
    const ColVec<1> r{};

    ColVec<1> u{};
    for (size_t k = 0; k < 100; ++k) {
        u = controller.control(r, x);
        x = plant_step(x, u(0));
    }
    CHECK(u(0) == doctest::Approx(0.3).epsilon(1e-6));

    // Runtime override takes effect without re-synthesis.
    controller.set_input_target(ColVec<1>{-0.2});
    for (size_t k = 0; k < 100; ++k) {
        u = controller.control(r, x);
        x = plant_step(x, u(0));
    }
    CHECK(u(0) == doctest::Approx(-0.2).epsilon(1e-6));
}

TEST_CASE("automatic scaling is consistent with the unscaled problem") {
    // Finite ranges give su = sy = 20; compensating the weights by the scale
    // factors must reproduce the sentinel-range (unit-scale) trajectory.
    design::MPCConstraints<1, 1> wide{};
    wide.u_min = ColVec<1>{-10.0};
    wide.u_max = ColVec<1>{10.0};
    wide.y_min = ColVec<1>{-10.0};
    wide.y_max = ColVec<1>{10.0};

    design::MPCWeights<1, 1> weights_scaled{};
    weights_scaled.Qy = Matrix<1, 1>{{400.0}}; // 1 · sy²
    weights_scaled.Qy_terminal = weights_scaled.Qy;
    weights_scaled.Rdu = Matrix<1, 1>{{0.05 * 400.0}}; // 0.05 · su²

    design::MPCWeights<1, 1> weights_unit{};
    weights_unit.Rdu = Matrix<1, 1>{{0.05}};

    const auto art_scaled = design::state_mpc<12, 4>(double_integrator, weights_scaled, wide);
    const auto art_unit = design::state_mpc<12, 4>(double_integrator, weights_unit);
    REQUIRE(art_scaled.success);
    REQUIRE(art_unit.success);
    CHECK(art_scaled.scale_u(0) == doctest::Approx(20.0));
    CHECK(art_scaled.scale_y(0) == doctest::Approx(20.0));

    MPC             scaled{art_scaled};
    MPC             unit{art_unit};
    ColVec<2>       xs{};
    ColVec<2>       xu{};
    const ColVec<1> r{1.0};

    for (size_t k = 0; k < 100; ++k) {
        const auto us = scaled.control(r, xs);
        const auto uu = unit.control(r, xu);
        CHECK(us(0) == doctest::Approx(uu(0)).epsilon(1e-8));
        xs = plant_step(xs, us(0));
        xu = plant_step(xu, uu(0));
    }
}

TEST_CASE("scaling handles badly-conditioned mixed units") {
    // Milli-scale actuator driving a kilo-scale output: DC gain 1e6.
    constexpr StateSpace<1, 1, 1> plant{
        .A = {{0.9}},
        .B = {{1000.0}},
        .C = {{100.0}},
        .Ts = 0.01,
    };

    design::MPCConstraints<1, 1> limits{};
    limits.u_min = ColVec<1>{-1.0e-3};
    limits.u_max = ColVec<1>{1.0e-3};
    limits.y_min = ColVec<1>{0.0};
    limits.y_max = ColVec<1>{1000.0};

    const auto art = design::state_mpc<10, 4>(plant, design::MPCWeights<1, 1>{}, limits);
    REQUIRE(art.success);

    MPC             controller{art};
    ColVec<1>       x{};
    const ColVec<1> r{500.0};

    for (size_t k = 0; k < 200; ++k) {
        const auto u = controller.control(r, x);
        CHECK(controller.last_status() == design::QPStatus::Success);
        CHECK(u(0) <= 1.0e-3 + 1e-12);
        CHECK(u(0) >= -1.0e-3 - 1e-12);
        x = ColVec<1>{0.9 * x(0) + 1000.0 * u(0)};
    }
    CHECK(100.0 * x(0) == doctest::Approx(500.0).epsilon(0.02));
}

TEST_CASE("mpc rejects invalid plants") {
    // Continuous-time (Ts = 0): must be discretized first.
    StateSpace<2, 1, 1> continuous = double_integrator;
    continuous.Ts = 0.0;
    CHECK(!design::state_mpc<5>(continuous).success);

    // Feedthrough (D != 0) is unsupported.
    StateSpace<2, 1, 1> feedthrough = double_integrator;
    feedthrough.D = Matrix<1, 1>{{1.0}};
    CHECK(!design::state_mpc<5>(feedthrough).success);

    // Non-positive slack penalty cannot keep the Hessian PD.
    design::MPCWeights<1, 1> weights{};
    weights.ecr_weight = 0.0;
    CHECK(!design::state_mpc<5>(double_integrator, weights).success);
}

TEST_CASE("mpc is constexpr") {
    constexpr auto art = design::state_mpc<3>(double_integrator);
    static_assert(art.success);
    static_assert(art.max_qp_iterations > 0);
}

TEST_CASE("artifact float conversion maps the unbounded sentinels safely") {
    constexpr auto art = design::state_mpc<3>(double_integrator);
    constexpr auto art_f = art.as<float>();
    static_assert(art_f.success);

    // Default constraints are unbounded: the double sentinel must arrive as
    // the float sentinel, not as an overflowed cast.
    static_assert(art_f.constraints.u_max(0) == std::numeric_limits<float>::max());
    static_assert(art_f.constraints.u_min(0) == std::numeric_limits<float>::lowest());
    static_assert(art_f.b_bound(0) == std::numeric_limits<float>::max());
}

TEST_CASE("MPC runtime runs in float") {
    design::MPCWeights<1, 1> weights{};
    weights.Rdu = Matrix<1, 1>{{0.05}};

    design::MPCConstraints<1, 1> limits{};
    limits.u_min = ColVec<1>{-0.5};
    limits.u_max = ColVec<1>{0.5};

    const auto art = design::state_mpc<10, 4>(double_integrator, weights, limits);
    REQUIRE(art.success);

    MPC                    controller{art.as<float>()};
    ColVec<2, float>       x{};
    const ColVec<1, float> r{1.0F};

    for (size_t k = 0; k < 150; ++k) {
        const auto u = controller.control(r, x);
        CHECK(u(0) <= 0.5F + 1e-6F);
        CHECK(u(0) >= -0.5F - 1e-6F);
        x = ColVec<2, float>{
            x(0) + (0.1F * x(1)) + (0.005F * u(0)),
            x(1) + (0.1F * u(0)),
        };
    }
    CHECK(x(0) == doctest::Approx(1.0).epsilon(0.05));
}

namespace {

// Solver policy that counts invocations and forwards to the default solver —
// exercises the pluggable-solver seam with a stateful policy.
struct CountingSolver {
    size_t calls = 0;

    template<size_t NV, size_t NI, typename T>
    [[nodiscard]] constexpr design::QPResult<NV, NI, T> operator()(
        const Matrix<NV, NV, T>& H,
        const ColVec<NV, T>&     f,
        const Matrix<NI, NV, T>& A,
        const ColVec<NI, T>&     b,
        size_t                   max_iterations
    ) {
        ++calls;
        return design::solve_qp(H, f, A, b, max_iterations);
    }
};

} // namespace

TEST_CASE("MPC accepts a custom QP solver policy") {
    design::MPCWeights<1, 1> weights{};
    weights.Rdu = Matrix<1, 1>{{0.05}};

    const auto art = design::state_mpc<8, 4>(double_integrator, weights);
    REQUIRE(art.success);

    MPC reference{art};
    MPC custom{art, CountingSolver{}};

    ColVec<2>       xr{0.5, 0.0};
    ColVec<2>       xc{0.5, 0.0};
    const ColVec<1> r{1.0};

    for (size_t k = 0; k < 20; ++k) {
        const auto ur = reference.control(r, xr);
        const auto uc = custom.control(r, xc);
        CHECK(uc(0) == doctest::Approx(ur(0)).epsilon(1e-12));
        xr = plant_step(xr, ur(0));
        xc = plant_step(xc, uc(0));
    }
    CHECK(custom.qp_solver().calls == 20);
}

TEST_CASE("horizon advisor reads the dominant time constant") {
    // First-order plant, pole at 0.9: tau = -Ts/ln(0.9) = 0.9491 s;
    // 4*tau = 3.796 s at Ts = 0.1 -> NP = 38, NC = 20% = 8.
    constexpr StateSpace<1, 1, 1> plant{
        .A = {{0.9}},
        .B = {{1.0}},
        .C = {{1.0}},
        .Ts = 0.1,
    };

    constexpr auto h = design::suggest_mpc_horizon(plant);
    static_assert(h.success);
    CHECK(h.time_constant == doctest::Approx(0.94912).epsilon(1e-4));
    CHECK(h.prediction_horizon == 38);
    CHECK(h.control_horizon == 8);
}

TEST_CASE("horizon advisor falls back to the settling-time overload on integrating plants") {
    // A pure integrator chain has no strictly stable pole to read tau from.
    const auto from_poles = design::suggest_mpc_horizon(double_integrator);
    CHECK(!from_poles.success);

    // The explicit target works anywhere: 2 s at Ts = 0.1 -> NP = 20, NC = 4.
    const auto h = design::suggest_mpc_horizon(2.0, 0.1);
    REQUIRE(h.success);
    CHECK(h.prediction_horizon == 20);
    CHECK(h.control_horizon == 4);
}

TEST_CASE("unconstrained analysis gain matches the LQR-equivalent design") {
    const Matrix<2, 2> Q{{1.0, 0.0}, {0.0, 0.1}};
    const Matrix<1, 1> R{{0.5}};
    const auto         lqr = design::discrete_lqr(double_integrator_full.A, double_integrator_full.B, Q, R);
    REQUIRE(lqr.success);

    design::MPCWeights<1, 2> weights{};
    weights.Qy = Q;
    weights.Qy_terminal = lqr.S;
    weights.Rdu = Matrix<1, 1>{};
    weights.Ru = R;

    const auto art = design::state_mpc<8>(double_integrator_full, weights);
    REQUIRE(art.success);

    const auto models = design::build_mpc_analysis_models(double_integrator_full, art);
    REQUIRE(models.success);

    // du = -K [x; u_prev]: the x block is the LQR gain, the u_prev block is I
    // (u = u_prev + du = -K_lqr x exactly).
    CHECK(models.K(0, 0) == doctest::Approx(lqr.K(0, 0)).epsilon(1e-8));
    CHECK(models.K(0, 1) == doctest::Approx(lqr.K(0, 1)).epsilon(1e-8));
    CHECK(models.K(0, 2) == doctest::Approx(1.0).epsilon(1e-8));
}

TEST_CASE("unconstrained closed-loop model is stable with unit DC gain") {
    design::MPCWeights<1, 1> weights{};
    weights.Rdu = Matrix<1, 1>{{0.05}};

    const auto art = design::state_mpc<15, 5>(double_integrator, weights);
    REQUIRE(art.success);

    const auto models = design::build_mpc_analysis_models(double_integrator, art);
    REQUIRE(models.success);

    // All closed-loop poles inside the unit circle...
    const auto eig = mat::compute_eigenvalues(models.closed_loop.A);
    for (size_t i = 0; i < 3; ++i) {
        CHECK(damp::abs(eig.values(i)) < 1.0);
    }
    // ...and the r -> y DC gain is exactly 1 (integral action).
    const auto dc = *eval_frf(models.closed_loop, damp::complex<double>{1.0, 0.0});
    CHECK(dc(0, 0).real() == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(dc(0, 0).imag() == doctest::Approx(0.0));
}

TEST_CASE("MPC reset and bumpless seeding") {
    const auto art = design::state_mpc<5>(double_integrator);
    REQUIRE(art.success);

    MPC controller{art};
    controller.set_previous_control(ColVec<1>{0.3});
    CHECK(controller.previous_control()(0) == doctest::Approx(0.3));

    controller.reset();
    CHECK(controller.previous_control()(0) == doctest::Approx(0.0));
    CHECK(controller.last_iterations() == 0);
    CHECK(controller.last_slack() == doctest::Approx(0.0));
}
