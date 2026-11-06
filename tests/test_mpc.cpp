// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>
#include <limits>

#include "damp/backend.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/design/qp.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/status.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

constexpr double Ts = 0.1;

// Double integrator, ZOH-discretized at Ts = 0.1, the state measured, the position the controlled output.
constexpr StateSpace<2, 1, 1> double_integrator{
    .A = {{1.0, Ts}, {0.0, 1.0}},
    .B = {{0.5 * Ts * Ts}, {Ts}},
    .C = {{1.0, 0.0}},
    .Ts = Ts,
};

// Same plant with every state a controlled output (C = I), for the LQR-equivalence tests.
constexpr StateSpace<2, 1, 2> double_integrator_full{
    .A = {{1.0, Ts}, {0.0, 1.0}},
    .B = {{0.5 * Ts * Ts}, {Ts}},
    .C = {{1.0, 0.0}, {0.0, 1.0}},
    .Ts = Ts,
};

/// The move weight @p rdu, the position weighted by the default Qy = 1.
design::MPCWeights<1, 1> move_weight(double rdu) {
    design::MPCWeights<1, 1> w{};
    w.Rdu = rdu;
    return w;
}

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

    const auto art = design::mpc<8, 8>(double_integrator_full, weights, {}); // a free move every step: equals LQR
    REQUIRE(art.success);

    MPC             controller{art};
    ColVec<2>       x{1.0, -0.5};
    const ColVec<2> r{}; // regulation

    for (size_t k = 0; k < 30; ++k) {
        const auto   u_mpc = controller.step(r, x);
        const double u_lqr = -((lqr.K(0, 0) * x(0)) + (lqr.K(0, 1) * x(1)));
        CHECK(u_mpc(0) == doctest::Approx(u_lqr).epsilon(1e-6));
        CHECK(controller.last_status() == design::QPStatus::Success);
        CHECK(controller.last_slack() == doctest::Approx(0.0));
        x = plant_step(x, u_mpc(0));
    }
}

TEST_CASE("constrained MPC tracks a step within u and du limits") {
    auto weights = move_weight(0.05);

    design::MPCConstraints<1, 1> limits{};
    limits.u_lim = {
        .lower = -0.5,
        .upper = 0.5,
    };
    limits.du_lim = {
        .lower = -0.2,
        .upper = 0.2,
    };

    const auto art = design::mpc<20, 8>(double_integrator, weights, limits);
    REQUIRE(art.success);

    MPC             controller{art};
    ColVec<2>       x{};
    const ColVec<1> r{1.0};

    double u_last = 0.0;
    double u_peak = 0.0;
    for (size_t k = 0; k < 200; ++k) {
        const auto u = controller.step(r, x);
        CHECK(controller.last_status() == design::QPStatus::Success);
        CHECK(controller.last_iterations() <= art.mpc.max_qp_iterations);

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

TEST_CASE("set_limits changes the input limit at run time; a bad limit is refused") {
    design::MPCConstraints<1, 1> limits{};
    limits.u_lim = {
        .lower = -0.5,
        .upper = 0.5,
    };
    const auto art = design::mpc<20, 8>(double_integrator, move_weight(0.05), limits);
    REQUIRE(art.success);

    MPC  controller{art};
    auto lim = controller.limits();
    CHECK(lim.u_lim.upper(0) == 0.5);
    CHECK(lim.du_lim.upper(0) == std::numeric_limits<double>::max()); // the design left Δu unbounded

    lim.u_lim.upper(0) = 0.2;
    REQUIRE(controller.set_limits(lim));
    CHECK(controller.limits().u_lim.upper(0) == 0.2);

    ColVec<2>       x{};
    const ColVec<1> r{1.0};
    double          u_peak = 0.0;
    for (size_t k = 0; k < 400; ++k) {
        const auto u = controller.step(r, x);
        REQUIRE(controller.last_status() == design::QPStatus::Success);
        CHECK(u(0) <= 0.2 + 1e-9); // the new limit, not the design's 0.5
        u_peak = (u(0) > u_peak) ? u(0) : u_peak;
        x = plant_step(x, u(0));
    }
    CHECK(u_peak == doctest::Approx(0.2));
    CHECK(x(0) == doctest::Approx(1.0).epsilon(0.02));

    auto bad = controller.limits();
    bad.u_lim.lower(0) = 0.3; // above u.upper
    CHECK_FALSE(controller.set_limits(bad));
    CHECK(controller.limits().u_lim.lower(0) == -0.5); // unchanged
}

TEST_CASE("hard output constraint caps the predicted response") {
    auto weights = move_weight(0.01);

    design::MPCConstraints<1, 1> limits{};
    limits.u_lim = {
        .lower = -1.0,
        .upper = 1.0,
    };
    limits.y_lim.upper(0) = 1.02; // overshoot limiter on the position
    limits.y_upper_ecr(0) = 0.0;  // hard: no slack allowed

    const auto art = design::mpc<20, 8>(double_integrator, weights, limits);
    REQUIRE(art.success);

    MPC             controller{art};
    ColVec<2>       x{};
    const ColVec<1> r{1.0};

    for (size_t k = 0; k < 300; ++k) {
        const auto u = controller.step(r, x);
        x = plant_step(x, u(0));
        CHECK(x(0) <= 1.02 + 1e-6); // plant == model, so the predicted cap holds exactly
    }
    CHECK(x(0) == doctest::Approx(1.0).epsilon(0.02));
}

TEST_CASE("soft output constraint recovers from an infeasible start") {
    // Start above the output ceiling with limited authority: no admissible
    // input sequence satisfies the hard rows, so a hard configuration reports
    // Infeasible while the default (soft) configuration keeps solving.
    auto weights = move_weight(0.01);

    design::MPCConstraints<1, 1> limits{};
    limits.u_lim = {
        .lower = -1.0,
        .upper = 1.0,
    };
    limits.y_lim.upper(0) = 1.0;

    const ColVec<2> x0{2.0, 0.0}; // x₀ = 2 > y.upper
    const ColVec<1> r{0.5};

    // Hard variant: first solve is infeasible and the fallback holds the input.
    {
        design::MPCConstraints<1, 1> hard = limits;
        hard.y_upper_ecr(0) = 0.0;
        const auto art = design::mpc<20, 8>(double_integrator, weights, hard);
        REQUIRE(art.success);
        MPC        controller{art};
        const auto u = controller.step(r, x0);
        CHECK(controller.last_status() == design::QPStatus::Infeasible);
        CHECK(u(0) == doctest::Approx(0.0)); // Δu falls back to zero
    }

    // Soft variant (default ECR = 1): solvable from the violated start, slack
    // active initially, and the loop still converges to the reference.
    {
        const auto art = design::mpc<20, 8>(double_integrator, weights, limits);
        REQUIRE(art.success);
        MPC controller{art};

        ColVec<2>  x = x0;
        const auto u0 = controller.step(r, x);
        CHECK(controller.last_status() == design::QPStatus::Success);
        CHECK(controller.last_slack() > 0.0);
        x = plant_step(x, u0(0));

        for (size_t k = 0; k < 400; ++k) {
            const auto u = controller.step(r, x);
            CHECK(controller.last_status() == design::QPStatus::Success);
            x = plant_step(x, u(0));
        }
        CHECK(x(0) == doctest::Approx(0.5).epsilon(0.02));
        CHECK(controller.last_slack() == doctest::Approx(0.0).epsilon(1e-6));
    }
}

TEST_CASE("an infeasible tick holds the last applied input inside the u box") {
    design::MPCConstraints<1, 1> hard{};
    hard.u_lim = {.lower = -1.0, .upper = 1.0};
    hard.y_lim.upper(0) = 1.0;
    hard.y_upper_ecr(0) = 0.0;
    const auto art = design::mpc<20, 8>(double_integrator, move_weight(0.01), hard);
    REQUIRE(art.success);
    MPC             controller{art};
    ColVec<2>       x{};
    const ColVec<1> r{0.5};
    double          prev = 0.0;
    for (size_t k = 0; k < 20; ++k) {
        const auto u = controller.step(r, x);
        CHECK(u(0) >= -1.0);
        CHECK(u(0) <= 1.0);
        REQUIRE(controller.last_status() == design::QPStatus::Success);
        prev = u(0);
        x = plant_step(x, u(0));
    }
    x = ColVec<2>{2.0, 0.0};
    const auto held = controller.step(r, x);
    CHECK(controller.last_status() == design::QPStatus::Infeasible);
    CHECK(held(0) == doctest::Approx(prev));
    CHECK(held(0) >= -1.0);
    CHECK(held(0) <= 1.0);
}

TEST_CASE("MaxIterations with hard output rows holds unless the iterate is primal feasible") {
    design::MPCConstraints<1, 1> limits{};
    limits.u_lim = {.lower = -1.0, .upper = 1.0};
    limits.y_lim.upper(0) = 1.02;
    limits.y_upper_ecr(0) = 0.0;
    limits.max_qp_iterations = 1;
    const auto art = design::mpc<30, 10>(double_integrator, move_weight(0.01), limits);
    REQUIRE(art.success);
    MPC             controller{art};
    ColVec<2>       x{};
    const ColVec<1> r{1.0};
    const auto      u0 = controller.step(r, x);
    CHECK(u0(0) >= -1.0);
    CHECK(u0(0) <= 1.0);
    const auto st0 = controller.last_status();
    if (st0 == design::QPStatus::MaxIterations || st0 == design::QPStatus::Degenerate) {
        CHECK(u0(0) == doctest::Approx(0.0)); // no previous move
    }
    x = plant_step(x, u0(0));
    const auto u1 = controller.step(r, x);
    CHECK(u1(0) >= -1.0);
    CHECK(u1(0) <= 1.0);
    const auto st1 = controller.last_status();
    if (st1 != design::QPStatus::Success) {
        CHECK(u1(0) == doctest::Approx(u0(0)));
    }
}

TEST_CASE("hard output rows with ill-scaled u and y boxes stay inside both") {
    design::MPCConstraints<1, 1> limits{};
    limits.u_lim = {.lower = -1.0e6, .upper = 1.0e6};
    limits.y_lim.upper(0) = 1.02;
    limits.y_upper_ecr(0) = 0.0;
    const auto art = design::mpc<20, 8>(double_integrator, move_weight(0.01), limits);
    REQUIRE(art.success);
    MPC             controller{art};
    ColVec<2>       x{};
    const ColVec<1> r{1.0};
    for (size_t k = 0; k < 80; ++k) {
        const auto u = controller.step(r, x);
        CHECK(u(0) >= -1.0e6);
        CHECK(u(0) <= 1.0e6);
        x = plant_step(x, u(0));
        CHECK(x(0) <= 1.02 + 1e-6);
    }
}

TEST_CASE("measured-disturbance feedforward cancels a known load") {
    // The disturbance enters like the input (Bd = B). With the MD channel the
    // controller sees it in the prediction and cancels it exactly at steady
    // state; without it, the same load leaves a residual error.
    const Matrix<2, 1> Bd = double_integrator.B;

    auto weights = move_weight(0.05);

    design::MPCConstraints<1, 1> limits{};
    limits.u_lim = {
        .lower = -2.0,
        .upper = 2.0,
    };

    const ColVec<1> r{0.5};
    const double    d_load = 0.4;

    // With MD feedforward:
    const auto art_md = design::mpc<15, 5>(double_integrator, weights, limits, design::MeasuredDisturbance{Bd});
    REQUIRE(art_md.success);
    MPC       md_controller{art_md};
    ColVec<2> x_md{};
    double    u_md = 0.0;
    for (size_t k = 0; k < 400; ++k) {
        const auto u = md_controller.step(r, ColVec<1>{d_load}, x_md);
        u_md = u(0);
        x_md = plant_step(x_md, u_md + d_load);
    }

    // Baseline without MD, same load applied to the plant:
    const auto art_base = design::mpc<15, 5>(double_integrator, weights, limits);
    REQUIRE(art_base.success);
    MPC       base_controller{art_base};
    ColVec<2> x_base{};
    for (size_t k = 0; k < 400; ++k) {
        const auto u = base_controller.step(r, x_base);
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
    weights.Ru = 1.0;
    weights.Rdu = 0.1;
    weights.u_target = ColVec<1>{0.3};

    const auto art = design::mpc<5, 5>(double_integrator, weights, {});
    REQUIRE(art.success);

    MPC             controller{art};
    ColVec<2>       x{};
    const ColVec<1> r{};

    ColVec<1> u{};
    for (size_t k = 0; k < 100; ++k) {
        u = controller.step(r, x);
        x = plant_step(x, u(0));
    }
    CHECK(u(0) == doctest::Approx(0.3).epsilon(1e-6));

    // Runtime override takes effect without re-synthesis.
    controller.set_input_target(ColVec<1>{-0.2});
    for (size_t k = 0; k < 100; ++k) {
        u = controller.step(r, x);
        x = plant_step(x, u(0));
    }
    CHECK(u(0) == doctest::Approx(-0.2).epsilon(1e-6));
}

TEST_CASE("automatic scaling is consistent with the unscaled problem") {
    // Finite ranges give su = sy = 20; compensating the weights by the scale
    // factors must reproduce the sentinel-range (unit-scale) trajectory.
    design::MPCConstraints<1, 1> wide{};
    wide.u_lim = {
        .lower = -10.0,
        .upper = 10.0,
    };
    wide.y_lim.lower(0) = -10.0;
    wide.y_lim.upper(0) = 10.0;

    auto weights_scaled = move_weight(0.05 * 400.0); // 0.05 · su²
    weights_scaled.Qy(0, 0) = 400.0;                 // 1 · sy²
    weights_scaled.Qy_terminal = weights_scaled.Qy;

    const auto weights_unit = move_weight(0.05);

    const auto art_scaled = design::mpc<12, 4>(double_integrator, weights_scaled, wide);
    const auto art_unit = design::mpc<12, 4>(double_integrator, weights_unit, {});
    REQUIRE(art_scaled.success);
    REQUIRE(art_unit.success);
    CHECK(art_scaled.mpc.scale_u(0) == doctest::Approx(20.0));
    CHECK(art_scaled.mpc.scale_y(0) == doctest::Approx(20.0));

    MPC             scaled{art_scaled};
    MPC             unit{art_unit};
    ColVec<2>       xs{};
    ColVec<2>       xu{};
    const ColVec<1> r{1.0};

    for (size_t k = 0; k < 100; ++k) {
        const auto us = scaled.step(r, xs);
        const auto uu = unit.step(r, xu);
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
    limits.u_lim = {
        .lower = -1.0e-3,
        .upper = 1.0e-3,
    };
    limits.y_lim = {
        .lower = 0.0,
        .upper = 1000.0,
    };

    const auto art = design::mpc<10, 4>(plant, design::MPCWeights<1, 1>{}, limits);
    REQUIRE(art.success);

    MPC             controller{art};
    ColVec<1>       x{};
    const ColVec<1> r{500.0};

    for (size_t k = 0; k < 200; ++k) {
        const auto u = controller.step(r, x);
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
    CHECK(design::mpc<5, 5>(continuous, {}, {}).status == design::DesignStatus::invalid_sample_time);

    // Feedthrough (D != 0) is unsupported.
    StateSpace<2, 1, 1> feedthrough = double_integrator;
    feedthrough.D(0, 0) = 1.0;
    CHECK(!design::mpc<5, 5>(feedthrough, {}, {}).success);

    // Non-positive slack penalty cannot keep the Hessian PD.
    design::MPCWeights<1, 1> weights{};
    weights.ecr_weight = 0.0;
    CHECK(!design::mpc<5, 5>(double_integrator, weights, {}).success);
}

TEST_CASE("mpc is constexpr") {
    constexpr auto art = design::mpc<3, 3>(double_integrator, {}, {});
    static_assert(art.success);
    static_assert(art.mpc.max_qp_iterations > 0);
}

TEST_CASE("Law names an MPC member's type by its design") {
    static constexpr auto art = design::mpc<3, 3>(double_integrator, {}, {});
    struct Loop {
        Law<art, double> mpc; // the design's law
    };
    Loop            loop{};
    MPC             plain{art};
    const ColVec<1> r{1.0};
    const ColVec<2> x{0.0, 0.0};
    CHECK(loop.mpc.step(r, x)(0) == plain.step(r, x)(0));
}

TEST_CASE("artifact float conversion maps the unbounded sentinels safely") {
    constexpr auto art = design::mpc<3, 3>(double_integrator, {}, {});
    constexpr auto art_f = art.as<float>();
    static_assert(art_f.success);

    // Default constraints are unbounded: the double sentinel must arrive as
    // the float sentinel, not as an overflowed cast.
    static_assert(art_f.mpc.constraints.u_lim.upper(0) == std::numeric_limits<float>::max());
    static_assert(art_f.mpc.constraints.u_lim.lower(0) == -std::numeric_limits<float>::max());
    static_assert(art_f.mpc.b_bound(0) == std::numeric_limits<float>::max());
}

TEST_CASE("MPC runtime runs in float") {
    auto weights = move_weight(0.05);

    design::MPCConstraints<1, 1> limits{};
    limits.u_lim = {
        .lower = -0.5,
        .upper = 0.5,
    };

    const auto art = design::mpc<10, 4>(double_integrator, weights, limits);
    REQUIRE(art.success);

    MPC                    controller{art.as<float>()};
    ColVec<2, float>       x{};
    const ColVec<1, float> r{1.0F};

    for (size_t k = 0; k < 150; ++k) {
        const auto u = controller.step(r, x);
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
    auto weights = move_weight(0.05);

    const auto art = design::mpc<8, 4>(double_integrator, weights, {});
    REQUIRE(art.success);

    MPC reference{art};
    MPC custom{art, CountingSolver{}};

    ColVec<2>       xr{0.5, 0.0};
    ColVec<2>       xc{0.5, 0.0};
    const ColVec<1> r{1.0};

    for (size_t k = 0; k < 20; ++k) {
        const auto ur = reference.step(r, xr);
        const auto uc = custom.step(r, xc);
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

TEST_CASE("to_ss of the LQR-equivalent design is the LQR gain on x and no memory of u_prev") {
    const Matrix<2, 2> Q{{1.0, 0.0}, {0.0, 0.1}};
    const Matrix<1, 1> R{{0.5}};
    const auto         lqr = design::discrete_lqr(double_integrator_full.A, double_integrator_full.B, Q, R);
    REQUIRE(lqr.success);

    design::MPCWeights<1, 2> weights{};
    weights.Qy = Q;
    weights.Qy_terminal = lqr.S;
    weights.Rdu = Matrix<1, 1>{};
    weights.Ru = R;

    const auto art = design::mpc<8, 8>(double_integrator_full, weights, {}); // a free move every step: equals LQR
    REQUIRE(art.success);

    // State u_prev, inputs [r; x]: Δu = −K_x x − K_u u_prev with K_x the LQR gain and K_u = 1,
    // so u = u_prev + Δu = −K_lqr x and the block keeps nothing of u_prev (A_c = C_c = 0).
    const auto ss = art.mpc.to_ss();
    CHECK(ss.Ts == doctest::Approx(Ts));
    CHECK(ss.A(0, 0) == doctest::Approx(0.0).scale(1.0).epsilon(1e-8));
    CHECK(ss.C(0, 0) == doctest::Approx(0.0).scale(1.0).epsilon(1e-8));
    CHECK(ss.D(0, 2) == doctest::Approx(-lqr.K(0, 0)).epsilon(1e-8));
    CHECK(ss.D(0, 3) == doctest::Approx(-lqr.K(0, 1)).epsilon(1e-8));
}

TEST_CASE("to_ss realizes the unconstrained runtime, and closed around the plant is stable with unit DC gain") {
    const auto art = design::mpc<15, 5>(double_integrator, move_weight(0.05), {});
    REQUIRE(art.success);
    const auto ss = art.mpc.to_ss(); // state u_prev, in [r; x], out u

    MPC       runtime{art};
    ColVec<1> xc{};
    ColVec<2> xp{0.3, -0.1};
    for (int k = 0; k < 12; ++k) {
        const ColVec<1> r{(k < 6) ? 1.0 : -0.5};
        const ColVec<3> in{r[0], xp[0], xp[1]};
        const ColVec<1> u = ss.C * xc + ss.D * in;
        CHECK(u[0] == doctest::Approx(runtime.step(r, xp)[0]).epsilon(1e-9));
        xc = ss.A * xc + ss.B * in;
        xp = plant_step(xp, u[0]);
    }

    // The plant reports [0; x], so the law's inputs [r; x] close by unity positive feedback.
    const StateSpace<2, 1, 3> plant_out{
        .A = double_integrator.A,
        .B = double_integrator.B,
        .C = Matrix<3, 2>{{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}},
        .Ts = Ts,
    };
    StateSpace<0, 3, 3> positive{.Ts = Ts};
    positive.D = -Matrix<3, 3>::identity();
    const auto loop = series(ss, plant_out);
    REQUIRE(loop.has_value());
    const auto cl = feedback(*loop, positive);
    REQUIRE(cl.has_value());

    // All closed-loop poles inside the unit circle...
    const auto eig = mat::compute_eigenvalues(cl->A);
    for (size_t i = 0; i < 3; ++i) {
        CHECK(damp::abs(eig.values(i)) < 1.0);
    }
    // ...and the r → y DC gain is exactly 1 (integral action).
    const auto dc = *eval_frf(*cl, damp::complex<double>{1.0, 0.0});
    CHECK(dc(1, 0).real() == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(dc(1, 0).imag() == doctest::Approx(0.0));
}

TEST_CASE("to_ss takes a measured disturbance between r and x, in the tick's order") {
    const auto art = design::mpc<15, 5>(double_integrator, move_weight(0.05), {}, design::MeasuredDisturbance{double_integrator.B});
    REQUIRE(art.success);
    const auto ss = art.mpc.to_ss(); // in [r; d; x]
    MPC        runtime{art};
    ColVec<1>  xc{};
    ColVec<2>  xp{};
    for (int k = 0; k < 12; ++k) {
        const ColVec<1> r{0.5};
        const ColVec<1> d{(k < 4) ? 0.0 : 0.4};
        const ColVec<4> in{r[0], d[0], xp[0], xp[1]};
        const ColVec<1> u = ss.C * xc + ss.D * in;
        CHECK(u[0] == doctest::Approx(runtime.step(r, d, xp)[0]).epsilon(1e-9));
        xc = ss.A * xc + ss.B * in;
        xp = plant_step(xp, u[0] + d[0]);
    }
}

TEST_CASE("MPC reset and bumpless seeding") {
    const auto art = design::mpc<5, 5>(double_integrator, {}, {});
    REQUIRE(art.success);

    MPC controller{art};
    controller.back_calculate(ColVec<1>{0.3}); // the input another controller applied
    controller.predict();
    CHECK(controller.applied()(0) == doctest::Approx(0.3));

    controller.reset();
    CHECK(controller.applied()(0) == doctest::Approx(0.0));
    CHECK(controller.last_iterations() == 0);
    CHECK(controller.last_slack() == doctest::Approx(0.0));
}

TEST_CASE("MPC step is predict, update, then control; the no-input step holds the last command") {
    auto                         weights = move_weight(0.05);
    design::MPCConstraints<1, 1> limits{};
    limits.u_lim = {
        .lower = -1.0,
        .upper = 1.0,
    };
    const auto art = design::mpc<15, 5>(double_integrator, weights, limits);
    REQUIRE(art.success);

    // A downstream clamp tighter than the MPC's: the applied input is the base of the next move.
    MPC             fused{art};
    MPC             split{art};
    ColVec<2>       x{};
    ColVec<1>       u_applied{};
    const ColVec<1> r{0.5};
    for (int k = 0; k < 40; ++k) {
        fused.back_calculate(u_applied);
        const ColVec<1> u_fused = fused.step(r, x);
        split.back_calculate(u_applied);
        split.predict();
        split.update(x);
        const ColVec<1> u_split = split.control(r);
        CHECK(u_fused(0) == doctest::Approx(u_split(0)).epsilon(1e-12));
        CHECK(fused.applied()(0) == u_applied(0));
        u_applied = ColVec<1>{damp::clamp(u_fused(0), -0.3, 0.3)};
        x = plant_step(x, u_applied(0));
    }

    // Without an applied input the step holds the command it returned last.
    MPC       implicit{art};
    MPC       explicit_u{art};
    ColVec<2> xi{};
    ColVec<2> xe{};
    for (int k = 0; k < 40; ++k) {
        const ColVec<1> ui = implicit.step(r, xi);
        explicit_u.back_calculate(explicit_u.last_command());
        const ColVec<1> ue = explicit_u.step(r, xe);
        CHECK(ui(0) == doctest::Approx(ue(0)).epsilon(1e-12));
        xi = plant_step(xi, ui(0));
        xe = plant_step(xe, ue(0));
    }
    CHECK(xi(0) == doctest::Approx(0.5).epsilon(1e-3));
}

namespace {

/// A damped mass, the position controlled, continuous.
constexpr StateSpace<2, 1, 1> kDampedMass{.A = {{0.0, 1.0}, {-2.0, -0.6}}, .B = {{0.0}, {1.5}}, .C = {{1.0, 0.0}}};

/// The position at every @p every-th tick of a 3 s step to 1 under @p art's runtime, on the plant sampled at @p h.
template<typename Result>
damp::array<double, 30> step_response(const Result& art, double h, size_t every) {
    const auto              sampled = discretize(kDampedMass, h).value();
    MPC                     law{art};
    ColVec<2>               x{};
    damp::array<double, 30> y{};
    for (size_t k = 0; k < 30 * every; ++k) {
        if (k % every == 0) {
            y[k / every] = x[0];
        }
        x = sampled.A * x + sampled.B * law.step(ColVec<1>{1.0}, x);
    }
    return y;
}

} // namespace

TEST_CASE("the continuous-plant design's QP is the Van Loan integral of its continuous cost") {
    // J = Σ_k ∫ q·(y − r)² dt over each period, + R_Δu·Δu²/Ts per move, + Ts·R_u·u² per move.
    constexpr double         h = 0.1;
    constexpr size_t         NP = 6;
    constexpr size_t         NC = 3;
    design::MPCWeights<1, 1> w{.Qy = 2.0, .Rdu = 0.3, .Ru = 0.5, .Qy_terminal = 2.0};
    const auto               art = design::mpc<NP, NC>(kDampedMass, h, w, {});
    REQUIRE(art.success);

    // Hand-built: the reference a constant state, the period's cost one Van Loan exponential.
    const Matrix<3, 3> A_aug{{0.0, 1.0, 0.0}, {-2.0, -0.6, 0.0}, {0.0, 0.0, 0.0}};
    const Matrix<3, 1> B_aug{{0.0}, {1.5}, {0.0}};
    const Matrix<1, 3> E{{1.0, 0.0, -1.0}}; // y − r on [x; r]
    const auto         period = design::discretize_lqr_cost(A_aug, B_aug, E.transpose() * (w.Qy * E), Matrix<1, 1>{}, h).value();
    const auto         sampled = discretize(kDampedMass, h).value();
    const ColVec<2>    x0{0.4, -0.2};
    const double       u_prev = 0.3;
    const double       r = 1.0;
    const auto         J = [&](const ColVec<NC>& moves) {
        double    cost = 0.0;
        ColVec<2> x = x0;
        double    u = u_prev;
        for (size_t k = 0; k < NP; ++k) {
            if (k < NC) {
                u += moves[k];
                cost += ((w.Rdu(0, 0) / h) * moves[k] * moves[k]) + (h * w.Ru(0, 0) * u * u);
            }
            const ColVec<3> p{x[0], x[1], r};
            cost += (p.transpose() * period.Q * p)(0, 0) + (2.0 * (p.transpose() * period.N)(0, 0) * u) + (period.R(0, 0) * u * u);
            x = sampled.A * x + sampled.B * ColVec<1>{u};
        }
        return cost;
    };

    // The QP at that point, unit scales and u_target = 0: ½ZᵀHZ + fᵀZ with f = F_x·x̃₀ − F_r·R′.
    const ColVec<3> xa{x0[0], x0[1], u_prev};
    ColVec<NP>      R{};
    for (size_t k = 0; k < NP; ++k) {
        R[k] = r;
    }
    const ColVec<NC> f = (art.mpc.F_x * xa) - (art.mpc.F_r * R);
    const double     scale = damp::abs(J(ColVec<NC>{}));
    for (size_t i = 0; i < NC; ++i) {
        ColVec<NC> ei{};
        ei[i] = 1.0;
        const double g = 0.5 * (J(ei) - J(ei * -1.0));
        CHECK(f[i] == doctest::Approx(g).epsilon(1e-9).scale(scale));
        for (size_t j = 0; j < NC; ++j) {
            ColVec<NC> ej{};
            ej[j] = 1.0;
            const double hij = J(ei + ej) - J(ei) - J(ej) + J(ColVec<NC>{});
            CHECK(art.mpc.H(i, j) == doctest::Approx(hij).epsilon(1e-9).scale(scale));
        }
    }
}

TEST_CASE("halving the period at fixed continuous weights keeps the closed loop") {
    // The same 2 s horizon and 0.4 s of moves at 0.1 s and 0.05 s. The weights are per second,
    // so the cost is the same functional of the response and the closed loop barely moves; the
    // same numbers taken per step (the discrete-plant form) weigh the output twice as often at
    // the finer period and give a different loop.
    const design::MPCWeights<1, 1> w{.Qy = 50.0, .Rdu = 0.5};
    const auto                     coarse = step_response(design::mpc<20, 4>(kDampedMass, 0.1, w, {}), 0.1, 1);
    const auto                     fine = step_response(design::mpc<40, 8>(kDampedMass, 0.05, w, {}), 0.05, 2);
    const auto                     coarse_d = step_response(design::mpc<20, 4>(discretize(kDampedMass, 0.1).value(), w, {}), 0.1, 1);
    const auto                     fine_d = step_response(design::mpc<40, 8>(discretize(kDampedMass, 0.05).value(), w, {}), 0.05, 2);
    double                         gap = 0.0;
    double                         gap_per_step = 0.0;
    for (size_t j = 0; j < 30; ++j) {
        gap = damp::max(gap, damp::abs(coarse[j] - fine[j]));
        gap_per_step = damp::max(gap_per_step, damp::abs(coarse_d[j] - fine_d[j]));
    }
    MESSAGE("largest position gap, 0.1 s against 0.05 s: continuous weights ", gap, ", per-step weights ", gap_per_step);
    CHECK(fine[29] == doctest::Approx(1.0).epsilon(1e-2)); // settling onto the step at 2.9 s
    CHECK(gap < 0.03);
    CHECK(gap < 0.5 * gap_per_step);
}

TEST_CASE("the continuous-plant terminal cost is the cost-to-go of the period cost: a free move every step gives one law at any horizon") {
    // With NC = NP and V_f the exact cost-to-go of the stage the QP sums, the first move is the
    // infinite-horizon law whatever NP is (Bellman), if and only if the two costs agree.
    const design::MPCWeights<1, 1>     w{.Qy = 5.0, .Rdu = 0.2, .Ru = 0.1};
    const design::MPCConstraints<1, 1> dare{.terminal = design::MPCTerminal::dare_cost};
    const auto                         short_h = design::mpc<3, 3>(kDampedMass, 0.1, w, dare);
    const auto                         long_h = design::mpc<7, 7>(kDampedMass, 0.1, w, dare);
    REQUIRE(short_h.success);
    REQUIRE(long_h.success);
    const auto a = short_h.mpc.to_ss(); // the feedback from x and u_prev is exact with a terminal cost
    const auto b = long_h.mpc.to_ss();
    CHECK(a.A(0, 0) == doctest::Approx(b.A(0, 0)).epsilon(1e-8));
    CHECK(a.D(0, 1) == doctest::Approx(b.D(0, 1)).epsilon(1e-8));
    CHECK(a.D(0, 2) == doctest::Approx(b.D(0, 2)).epsilon(1e-8));
}

TEST_CASE("mpc rejects any nonzero feedthrough, whatever its units") {
    StateSpace<2, 1, 1> sys{
        .A = {{1.0, 0.1}, {0.0, 1.0}},
        .B = {{0.005}, {0.1}},
        .C = {{1.0, 0.0}},
        .D = {{1e-13}}, // e.g. y in V, u in MW: real feedthrough, tiny number
        .Ts = 0.1,
    };
    CHECK_FALSE(design::mpc<5, 2>(sys, {}, {}).success);
    sys.D = Matrix<1, 1>{};
    CHECK(design::mpc<5, 2>(sys, {}, {}).success);
}

TEST_CASE("suggest_mpc_horizon does not read a defective integrator as a slow pole") {
    // Double integrator in a non-orthogonal basis with a large off-diagonal:
    // the computed eigenvalues sit O(√eps) off 1.
    const Matrix<2, 2> J{{1.0, 100.0}, {0.0, 1.0}};
    const Matrix<2, 2> Tm{{1.0, 0.3}, {0.7, 1.9}};
    const auto         Ti = Tm.inverse();
    REQUIRE(Ti.has_value());
    const StateSpace<2, 1, 1> sys{.A = Tm * J * *Ti, .B = {{0.0}, {1.0}}, .C = {{1.0, 0.0}}, .Ts = 0.1};
    CHECK_FALSE(design::suggest_mpc_horizon(sys).success);

    // A triple integrator's poles land O(eps^(1/3)) = 1e-6 off 1, outside a √eps band: it
    // suggested NP = 8686359.
    const Matrix<3, 3> J3{{1.0, 0.1, 0.005}, {0.0, 1.0, 0.1}, {0.0, 0.0, 1.0}};
    const Matrix<3, 3> T3{{1.0, 0.3, -0.2}, {0.7, 1.9, 0.4}, {0.1, -0.5, 1.3}};
    const auto         T3i = T3.inverse();
    REQUIRE(T3i.has_value());
    const StateSpace<3, 1, 1> triple{.A = T3 * J3 * *T3i, .B = {{0.0}, {0.0}, {1.0}}, .C = {{1.0, 0.0, 0.0}}, .Ts = 0.1};
    CHECK_FALSE(design::suggest_mpc_horizon(triple).success);
    const StateSpace<2, 1, 1> slow{.A = {{0.99, 0.0}, {0.0, 0.5}}, .B = {{1.0}, {1.0}}, .C = {{1.0, 1.0}}, .Ts = 0.1};
    CHECK(design::suggest_mpc_horizon(slow).time_constant == doctest::Approx(-0.1 / damp::log(0.99)));

    // A genuinely stable pole still yields a suggestion.
    const StateSpace<1, 1, 1> stable{.A = {{0.9}}, .B = {{1.0}}, .C = {{1.0}}, .Ts = 0.1};
    CHECK(design::suggest_mpc_horizon(stable).success);
    CHECK_FALSE(design::suggest_mpc_horizon(1e300, 1e-300).success); // NP not representable
}

TEST_CASE("a measured disturbance with a terminal cost settles on r: the target is the steady state under d") {
    //! The target ignored d, so V_f pulled the plan toward the d = 0 steady state: y = 0.969 for r = 1.
    const StateSpace<1, 1, 1>    plant{.A = {{0.9}}, .B = {{0.1}}, .C = {{1.0}}, .Ts = 0.1};
    design::MPCConstraints<1, 1> limits{};
    limits.terminal = design::MPCTerminal::dare_cost;
    const auto art = design::mpc<10, 10>(plant, move_weight(0.1), limits, design::MeasuredDisturbance{Matrix<1, 1>{{0.1}}});
    REQUIRE(art.success);
    MPC       law{art};
    double    x = 0.0;
    ColVec<1> u{};
    for (size_t k = 0; k < 400; ++k) {
        u = law.step(ColVec<1>{1.0}, ColVec<1>{2.0}, ColVec<1>{x});
        x = (0.9 * x) + (0.1 * u(0)) + (0.1 * 2.0);
    }
    CHECK(x == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(u(0) == doctest::Approx(-1.0).epsilon(1e-6));
}

TEST_CASE("as<float>() refuses a Hessian that does not factor in float") {
    //! It kept success: constinit compiled and every float tick returned NotPositiveDefinite.
    design::MPCConstraints<1, 1> limits{};
    limits.u_lim = {.lower = -1.0, .upper = 1.0};
    const auto art = design::mpc<40, 10>(double_integrator, move_weight(1e-3), limits);
    REQUIRE(art.success);
    const auto art_f = art.as<float>();
    CHECK_FALSE(art_f.success);
    CHECK(art_f.status == design::DesignStatus::singular);
    CHECK_FALSE(art_f.mpc.success);
}

TEST_CASE("one hard output side by designated initializer leaves the other side soft") {
    //! y_lim_ecr was a Bounds: {.upper = {0.0}} set the lower softness to −max(), a hard lower
    //! bound and ±inf under as<float>().
    const design::MPCConstraints<1, 1> limits{
        .u_lim = {.lower = -1.0, .upper = 1.0},
        .y_lim = {.lower = 0.2, .upper = 1.0},
        .y_upper_ecr = 0.0,
    };
    const auto art = design::mpc<10, 3>(double_integrator, move_weight(0.1), limits);
    REQUIRE(art.success);
    MPC        controller{art};
    const auto u = controller.step(ColVec<1>{0.5}, ColVec<2>{}); // y = 0 starts below the soft lower bound
    CHECK(controller.last_status() == design::QPStatus::Success);
    CHECK(u(0) > 0.0);
    constexpr size_t r_y_hi = 4 * 3;             // after the Δu and u rows of NC = 3 moves
    constexpr size_t r_y_lo = r_y_hi + (10 * 1); // after the NP = 10 upper output rows
    CHECK(art.mpc.A_con(r_y_hi, 3) == 0.0);      // the slack column: the upper side hard
    CHECK(art.mpc.A_con(r_y_lo, 3) == -1.0);     // the lower side soft
    CHECK(art.mpc.as<float>().A_con(r_y_lo, 3) == -1.0F);
}

TEST_CASE("design::mpc refuses NaN or inverted limits, negative softness and indefinite weights") {
    //! The design accepted what set_limits refuses (u_lim {1, −1}: std::clamp with lo > hi on
    //! every tick), and output weights design::lqr refuses.
    design::MPCConstraints<1, 1> inverted{};
    inverted.u_lim = {.lower = 1.0, .upper = -1.0};
    CHECK(design::mpc<6, 2>(double_integrator, move_weight(1.0), inverted).status == design::DesignStatus::invalid_parameter);
    design::MPCConstraints<1, 1> nan_du{};
    nan_du.du_lim.upper(0) = std::numeric_limits<double>::quiet_NaN();
    CHECK(design::mpc<6, 2>(double_integrator, move_weight(1.0), nan_du).status == design::DesignStatus::invalid_parameter);
    design::MPCConstraints<1, 1> negative_ecr{};
    negative_ecr.y_upper_ecr(0) = -1.0;
    CHECK(design::mpc<6, 2>(double_integrator, move_weight(1.0), negative_ecr).status == design::DesignStatus::bad_weights);

    design::MPCWeights<1, 1> negative{};
    negative.Qy = -0.05;
    negative.Qy_terminal = -0.05;
    CHECK(design::mpc<6, 2>(double_integrator, negative, {}).status == design::DesignStatus::bad_weights);
    design::MPCWeights<1, 2> asymmetric{};
    asymmetric.Qy = Matrix<2, 2>{{1.0, 5.0}, {0.0, 1.0}};
    CHECK(design::mpc<6, 2>(double_integrator_full, asymmetric, {}).status == design::DesignStatus::bad_weights);
    CHECK(design::mpc<6, 2>(double_integrator, move_weight(1.0), {}).success);
}

TEST_CASE("float MPC with soft output bounds, started outside them, converges") {
    //! The float active set reported Infeasible on most such ticks (κ at roundoff for soft rows
    //! of adjacent steps), the input was held and the double integrator ran off to y = −356.
    design::MPCConstraints<1, 1> limits{};
    limits.u_lim = {.lower = -1.0, .upper = 1.0};
    limits.du_lim = {.lower = -0.2, .upper = 0.2};
    limits.y_lim = {.lower = -0.5, .upper = 0.5};
    const auto art = design::mpc<30, 10>(double_integrator, move_weight(0.01), limits);
    REQUIRE(art.success);
    MPC       controller{art.as<float>()};
    ColVec<2> x{1.5, 0.0};
    size_t    infeasible = 0;
    for (size_t k = 0; k < 300; ++k) {
        const auto u = controller.step(ColVec<1, float>{0.0F}, x.as<float>());
        infeasible += (controller.last_status() == design::QPStatus::Infeasible) ? 1U : 0U;
        x = plant_step(x, static_cast<double>(u(0)));
    }
    CHECK(infeasible == 0);
    CHECK(x(0) == doctest::Approx(0.0).scale(1.0).epsilon(1e-3));
    CHECK(x(1) == doctest::Approx(0.0).scale(1.0).epsilon(1e-3));
}

TEST_CASE("back_calculate with a feedforward: the command is the law's share, and the loop settles on r") {
    //! The held input already carried u_ff, the user added it again: y = 1.092 for u_ff = 0.2.
    const StateSpace<1, 1, 1> plant{.A = {{0.9}}, .B = {{0.1}}, .C = {{1.0}}, .Ts = 0.1};
    const auto                art = design::mpc<10, 3>(plant, move_weight(0.1), {});
    REQUIRE(art.success);
    MPC    law{art};
    double x = 0.0;
    double total = 0.0;
    for (size_t k = 0; k < 400; ++k) {
        const auto command = law.step(ColVec<1>{1.0}, ColVec<1>{x});
        total = command(0) + 0.2;
        law.back_calculate(ColVec<1>{total}, ColVec<1>{0.2});
        x = (0.9 * x) + (0.1 * total);
    }
    CHECK(x == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(total == doctest::Approx(1.0).epsilon(1e-6));
}
