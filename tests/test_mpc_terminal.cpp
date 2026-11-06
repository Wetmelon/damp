// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/design/invariant_set.hpp"
#include "damp/design/qp.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

constexpr double Ts = 0.1;

/// The state measured; the position is the controlled output the tests track and bound.
constexpr StateSpace<2, 1, 1> double_integrator{
    .A = {{1.0, Ts}, {0.0, 1.0}},
    .B = {{0.5 * Ts * Ts}, {Ts}},
    .C = {{1.0, 0.0}},
    .Ts = Ts,
};

/// Every state a controlled output (C = I), for the LQR-equivalence test.
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

/// |u| ≤ 1 and the position within ±1 (@p hard: no slack).
design::MPCConstraints<1, 1> unit_box(bool hard) {
    design::MPCConstraints<1, 1> limits{};
    limits.u_lim = {
        .lower = -1.0,
        .upper = 1.0,
    };
    limits.y_lim.lower(0) = -1.0;
    limits.y_lim.upper(0) = 1.0;
    if (hard) {
        limits.y_lower_ecr(0) = 0.0;
        limits.y_upper_ecr(0) = 0.0;
    }
    return limits;
}

/// The four-state plant: two double integrators, the state measured; both positions tracked and bounded.
constexpr StateSpace<4, 2, 2> two_axes{
    .A = {{1.0, Ts, 0.0, 0.0}, {0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, Ts}, {0.0, 0.0, 0.0, 1.0}},
    .B = {{0.5 * Ts * Ts, 0.0}, {Ts, 0.0}, {0.0, 0.5 * Ts * Ts}, {0.0, Ts}},
    .C = {{1.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}},
    .Ts = Ts,
};

constexpr design::MPCWeights<2, 2> two_axes_weights() {
    design::MPCWeights<2, 2> weights{};
    weights.Rdu = Matrix<2, 2>::identity() * 0.1;
    return weights;
}

constexpr ColVec<2> plant_step(const ColVec<2>& x, double u) {
    return ColVec<2>{
        x(0) + (Ts * x(1)) + (0.5 * Ts * Ts * u),
        x(1) + (Ts * u),
    };
}

[[nodiscard]] bool has_facet(
    const Matrix<64, 2>& G,
    const ColVec<64>&    h,
    size_t               n,
    double               a0,
    double               a1,
    double               beta
) {
    const double na = damp::hypot(a0, a1);
    if (!(na > 0.0)) {
        return false;
    }
    const double tol = 64.0 * eps<double>() * (1.0 + na + damp::abs(beta));
    for (size_t i = 0; i < n; ++i) {
        const double g0 = G(i, 0);
        const double g1 = G(i, 1);
        const double ng = damp::hypot(g0, g1);
        if (!(ng > 0.0)) {
            continue;
        }
        const double s = ng / na;
        if ((g0 * a0) + (g1 * a1) <= 0.0) {
            continue;
        }
        if (damp::abs(g0 - (s * a0)) <= tol && damp::abs(g1 - (s * a1)) <= tol && damp::abs(h(i) - (s * beta)) <= tol) {
            return true;
        }
    }
    return false;
}

} // namespace

TEST_CASE("default mpc: no terminal rows, NI is the box formula") {
    const auto art = design::mpc<8, 4>(double_integrator, {}, {});
    REQUIRE(art.success);
    CHECK(art.status == design::DesignStatus::ok);
    CHECK(art.mpc.n_terminal == 0);
    CHECK(art.mpc.terminal == design::MPCTerminal::none);
    CHECK(std::remove_cvref_t<decltype(art.mpc)>::NI == (4 * 4) + (2 * 8));
    CHECK(art.mpc.P(0, 0) == doctest::Approx(0.0));
    CHECK(art.mpc.K_terminal(0, 0) == doctest::Approx(0.0));
}

TEST_CASE("unconstrained DARE terminal reproduces the LQR law") {
    const Matrix<2, 2> Q{{1.0, 0.0}, {0.0, 0.1}};
    const Matrix<1, 1> R{{0.5}};
    const auto         lqr = design::discrete_lqr(double_integrator_full.A, double_integrator_full.B, Q, R);
    REQUIRE(lqr.success);

    design::MPCWeights<1, 2> weights{};
    weights.Qy = Q;
    weights.Rdu = Matrix<1, 1>{};
    weights.Ru = R;

    design::MPCConstraints<1, 2> limits{};
    limits.terminal = design::MPCTerminal::dare_cost;

    const auto art = design::mpc<8, 8>(double_integrator_full, weights, limits);
    REQUIRE(art.success);
    CHECK(art.status == design::DesignStatus::ok);
    CHECK(art.mpc.n_terminal == 0);
    //! With no move penalty the applied input u = v + Δ is free each step, so the Δu-form
    //! cost-to-go on z = [e; v] is the LQR one on e alone: P = blkdiag(S, 0), Δ = −K e − v.
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            CHECK(art.mpc.P(i, j) == doctest::Approx(lqr.S(i, j)).epsilon(1e-8));
        }
        CHECK(art.mpc.P(i, 2) == doctest::Approx(0.0).scale(1.0).epsilon(1e-8));
        CHECK(art.mpc.K_terminal(0, i) == doctest::Approx(lqr.K(0, i)).epsilon(1e-8));
    }
    CHECK(art.mpc.P(2, 2) == doctest::Approx(0.0).scale(1.0).epsilon(1e-8));
    CHECK(art.mpc.K_terminal(0, 2) == doctest::Approx(1.0).epsilon(1e-8));

    MPC             controller{art};
    ColVec<2>       x{1.0, -0.5};
    const ColVec<2> r{};
    double          V_prev = -1.0;
    for (size_t k = 0; k < 30; ++k) {
        const auto   u_mpc = controller.step(r, x);
        const double u_lqr = -((lqr.K(0, 0) * x(0)) + (lqr.K(0, 1) * x(1)));
        CHECK(u_mpc(0) == doctest::Approx(u_lqr).epsilon(1e-6));
        CHECK(controller.last_status() == design::QPStatus::Success);
        CHECK(damp::isfinite(controller.last_cost()));
        const double V = (x(0) * ((lqr.S(0, 0) * x(0)) + (lqr.S(0, 1) * x(1)))) + (x(1) * ((lqr.S(1, 0) * x(0)) + (lqr.S(1, 1) * x(1))));
        if (k > 0) {
            CHECK(V <= V_prev + 1e-9);
        }
        V_prev = V;
        x = plant_step(x, u_mpc(0));
    }
}

TEST_CASE("short-horizon SIL is infeasible without the set and recursively feasible with it") {
    const auto weights = move_weight(0.01);
    auto       limits = unit_box(true);

    const ColVec<1> r{0.9}; // the target (x_t, u_t) = ([0.9, 0], 0)
    // Fast approach to a near-limit setpoint: 4-step reachable to the MPI set,
    // but a terminal-free short horizon over-commits on speed and later QP-fails.
    const ColVec<2> x0{0.5, 0.8};

    bool saw_infeasible = false;
    {
        const auto art = design::mpc<4, 2>(double_integrator, weights, limits);
        REQUIRE(art.success);
        MPC       controller{art};
        ColVec<2> x = x0;
        for (size_t k = 0; k < 40; ++k) {
            const auto u = controller.step(r, x);
            if (controller.last_status() == design::QPStatus::Infeasible) {
                saw_infeasible = true;
                break;
            }
            x = plant_step(x, u(0));
        }
    }
    CHECK(saw_infeasible);

    limits.terminal = design::MPCTerminal::invariant_set;
    const auto art = design::mpc<4, 4, 64>(double_integrator, weights, limits);
    REQUIRE(art.success);
    REQUIRE(art.mpc.n_terminal > 0);
    CHECK(art.mpc.terminal == design::MPCTerminal::invariant_set);

    MPC       controller{art};
    ColVec<2> x = x0;
    for (size_t k = 0; k < 80; ++k) {
        const auto u = controller.step(r, x);
        CHECK(controller.last_status() == design::QPStatus::Success);
        CHECK(u(0) <= 1.0 + 1e-9);
        CHECK(u(0) >= -1.0 - 1e-9);
        x = plant_step(x, u(0));
        CHECK(x(0) <= 1.0 + 1e-6);
        CHECK(x(0) >= -1.0 - 1e-6);
    }
}

TEST_CASE("maximal_invariant_set matches a hand-derived 2-state MPI polytope") {
    const Matrix<2, 2> A{{0.5, 0.0}, {0.0, 0.5}};
    const Matrix<4, 2> G0{{1.0, 0.0}, {-1.0, 0.0}, {0.0, 1.0}, {0.0, -1.0}};
    const ColVec<4>    h0{1.0, 1.0, 1.0, 1.0};
    const auto         set = design::maximal_invariant_set<64>(A, G0, h0);
    REQUIRE(set.success);
    CHECK(set.status == design::DesignStatus::ok);
    CHECK(set.n_rows >= 4);
    CHECK(has_facet(set.G, set.h, set.n_rows, 1.0, 0.0, 1.0));
    CHECK(has_facet(set.G, set.h, set.n_rows, -1.0, 0.0, 1.0));
    CHECK(has_facet(set.G, set.h, set.n_rows, 0.0, 1.0, 1.0));
    CHECK(has_facet(set.G, set.h, set.n_rows, 0.0, -1.0, 1.0));
}

TEST_CASE("maximal_invariant_set adds the first propagated cut of a shear map") {
    const Matrix<2, 2> A{{0.8, 0.5}, {0.0, 0.8}};
    const Matrix<4, 2> G0{{1.0, 0.0}, {-1.0, 0.0}, {0.0, 1.0}, {0.0, -1.0}};
    const ColVec<4>    h0{1.0, 1.0, 1.0, 1.0};
    const auto         set = design::maximal_invariant_set<64>(A, G0, h0);
    REQUIRE(set.success);
    CHECK(set.n_rows > 4);
    CHECK(has_facet(set.G, set.h, set.n_rows, 0.8, 0.5, 1.0));
    CHECK(has_facet(set.G, set.h, set.n_rows, -0.8, -0.5, 1.0));
}

TEST_CASE("maximal_invariant_set reports empty on inconsistent G0 h0") {
    const Matrix<1, 1> A{{0.5}};
    const Matrix<2, 1> G0{{1.0}, {-1.0}};
    const ColVec<2>    h0{-1.0, -1.0};
    const auto         set = design::maximal_invariant_set<8>(A, G0, h0);
    CHECK_FALSE(set.success);
    CHECK(set.status == design::DesignStatus::not_converged);
}

TEST_CASE("float conversion keeps terminal rows and the float tick stays in budget") {
    const auto weights = move_weight(0.01);
    auto       limits = unit_box(true);
    limits.terminal = design::MPCTerminal::invariant_set;

    const auto art = design::mpc<4, 4, 64>(double_integrator, weights, limits);
    REQUIRE(art.success);
    const auto art_f = art.as<float>();
    CHECK(art_f.success);
    CHECK(art_f.status == art.status);
    CHECK(art_f.mpc.n_terminal == art.mpc.n_terminal);
    CHECK(art_f.mpc.terminal == art.mpc.terminal);
    CHECK(art_f.mpc.P(0, 0) == doctest::Approx(art.mpc.P(0, 0)).epsilon(1e-5));
    CHECK(art_f.mpc.K_terminal(0, 0) == doctest::Approx(art.mpc.K_terminal(0, 0)).epsilon(1e-5));
    REQUIRE(art.mpc.n_terminal > 0);
    CHECK(art_f.mpc.G_terminal(0, 0) == doctest::Approx(art.mpc.G_terminal(0, 0)).epsilon(1e-5));
    CHECK(art_f.mpc.h_terminal(0) == doctest::Approx(art.mpc.h_terminal(0)).epsilon(1e-5));

    MPC                    controller{art_f};
    ColVec<2, float>       x{0.5F, 0.8F};
    const ColVec<1, float> r{0.9F};
    for (size_t k = 0; k < 40; ++k) {
        const auto u = controller.step(r, x);
        CHECK(controller.last_status() == design::QPStatus::Success);
        CHECK(controller.last_iterations() <= art_f.mpc.max_qp_iterations);
        x = ColVec<2, float>{
            x(0) + (0.1F * x(1)) + (0.005F * u(0)),
            x(1) + (0.1F * u(0)),
        };
    }
}

TEST_CASE("invariant_set with NT_max = 0 is refused as bad_weights") {
    design::MPCConstraints<1, 1> limits{};
    limits.terminal = design::MPCTerminal::invariant_set;
    const auto art = design::mpc<4, 2>(double_integrator, {}, limits);
    CHECK_FALSE(art.success);
    CHECK(art.status == design::DesignStatus::bad_weights);
}

TEST_CASE("terminal set moves with the target: every admissible target is inside, a near-bound overshoot is not") {
    const auto weights = move_weight(0.01);
    auto       limits = unit_box(false);
    limits.terminal = design::MPCTerminal::invariant_set;
    const auto art = design::mpc<4, 4, 64>(double_integrator, weights, limits);
    REQUIRE(art.success);
    REQUIRE(art.mpc.n_terminal > 0);

    //! Rows are G·z + G_target·[y_t; u_t] ≤ h, z = [e; v]: the largest row excess at (z, y_t).
    const auto excess = [&](const ColVec<3>& z, double y_t) {
        double w = -1.0;
        for (size_t r = 0; r < art.mpc.n_terminal; ++r) {
            double g = art.mpc.G_terminal_target(r, 0) * y_t;
            for (size_t j = 0; j < 3; ++j) {
                g += art.mpc.G_terminal(r, j) * z(j);
            }
            w = damp::max(w, g - art.mpc.h_terminal(r));
        }
        return w;
    };
    bool target_enters = false;
    for (size_t r = 0; r < art.mpc.n_terminal; ++r) {
        target_enters = target_enters || (art.mpc.G_terminal_target(r, 0) != 0.0);
    }
    CHECK(target_enters);
    //! The target itself (z = 0) is in its own set anywhere in the 5 %-tightened box.
    for (const double y_t : {-0.95, -0.5, 0.0, 0.5, 0.9, 0.95}) {
        CHECK(excess(ColVec<3>{}, y_t) <= 1e-12);
    }
    //! y = y_t + 0.11 = 1.01 is past the bound: excluded whatever the input.
    CHECK(excess(ColVec<3>{0.11, 0.0, 0.0}, 0.9) > 0.0);
}

TEST_CASE("4-state NP=10 dare_cost design is constexpr under default ops") {
    // The Δu-form terminal cost is a DARE on z = [e; v] (6 states here): cheap enough for constinit.
    constexpr auto art = [] {
        constexpr auto               weights = two_axes_weights();
        design::MPCConstraints<2, 2> limits{};
        limits.u_lim = {
            .lower = {-1.0, -1.0},
            .upper = {1.0, 1.0},
        };
        limits.terminal = design::MPCTerminal::dare_cost;
        return design::mpc<10, 4>(two_axes, weights, limits);
    }();
    static_assert(art.success);
    static_assert(art.mpc.P(4, 4) > 0.0); // the v block carries the move cost
}

TEST_CASE("4-state NP=10 invariant_set design") {
    // The target-augmented set lives on [z; q] (NX + NU + NY + NU = 10 here) and needs ~160
    // facets: it is designed at initialization with the capacity raised accordingly.
    const auto art = [] {
        constexpr auto               weights = two_axes_weights();
        design::MPCConstraints<2, 2> limits{};
        limits.u_lim = {
            .lower = {-1.0, -1.0},
            .upper = {1.0, 1.0},
        };
        limits.y_lim = {
            .lower = {-1.0, -1.0},
            .upper = {1.0, 1.0},
        };
        limits.terminal = design::MPCTerminal::invariant_set;
        return design::mpc<10, 10, 192>(two_axes, weights, limits);
    }();
    REQUIRE(art.success);
    CHECK(art.mpc.n_terminal > 64);
    CHECK(art.status == design::DesignStatus::ok);
}

TEST_CASE("maximal_invariant_set does not depend on the units of h") {
    //! The ridge LP stood in for max aᵀe with an absolute weight: past h ≈ 1e7 (1e3 in float) a
    //! direction where the set was still unbounded read as bounded, propagated facets were
    //! dropped and a non-invariant set came back ok (2 facets at 1e8 against 26).
    const Matrix<2, 2> A{{0.9, 0.5}, {0.0, 0.9}};
    const Matrix<2, 2> G0{{1.0, 0.0}, {-1.0, 0.0}};
    const auto         unit = design::maximal_invariant_set<64>(A, G0, ColVec<2>{1.0, 1.0});
    REQUIRE(unit.success);
    for (const double scale : {1e5, 1e8}) {
        const auto big = design::maximal_invariant_set<64>(A, G0, ColVec<2>{scale, scale});
        REQUIRE(big.success);
        CHECK(big.n_rows == unit.n_rows);
        for (size_t r = 0; r < unit.n_rows; ++r) {
            CHECK(big.h(r) == doctest::Approx(scale * unit.h(r)).epsilon(1e-9));
        }
    }
    const auto unit_f = design::maximal_invariant_set<64>(A.as<float>(), G0.as<float>(), ColVec<2, float>{1.0F, 1.0F});
    const auto big_f = design::maximal_invariant_set<64>(A.as<float>(), G0.as<float>(), ColVec<2, float>{1000.0F, 1000.0F});
    REQUIRE(unit_f.success);
    REQUIRE(big_f.success);
    CHECK(big_f.n_rows == unit_f.n_rows);
}

TEST_CASE("invariant_set with a control horizon shorter than the prediction horizon is refused") {
    //! The shifted plan holds the input over moves NC … NP − 1, so it is not admissible and the
    //! hard terminal rows went infeasible on the nominal plant (24 of 359 starts, NP = 4, NC = 2).
    auto limits = unit_box(true);
    limits.terminal = design::MPCTerminal::invariant_set;
    const auto short_nc = design::mpc<4, 2, 64>(double_integrator, move_weight(0.01), limits);
    CHECK_FALSE(short_nc.success);
    CHECK(short_nc.status == design::DesignStatus::invalid_parameter);
    CHECK(design::mpc<4, 4, 64>(double_integrator, move_weight(0.01), limits).success);
}

TEST_CASE("MaxIterations with an invariant_set holds unless the iterate is primal feasible") {
    auto limits = unit_box(true);
    limits.terminal = design::MPCTerminal::invariant_set;
    const ColVec<1> r{0.9};
    const ColVec<2> x0{0.5, 0.8};

    const auto art_ok = design::mpc<4, 4, 64>(double_integrator, move_weight(0.01), limits);
    REQUIRE(art_ok.success);
    REQUIRE(art_ok.mpc.n_terminal > 0);
    MPC        ok{art_ok};
    const auto u_ok = ok.step(r, x0);
    CHECK(ok.last_status() == design::QPStatus::Success);
    CHECK(u_ok(0) >= -1.0);
    CHECK(u_ok(0) <= 1.0);

    limits.max_qp_iterations = 1;
    const auto art = design::mpc<4, 4, 64>(double_integrator, move_weight(0.01), limits);
    REQUIRE(art.success);
    MPC        controller{art};
    const auto u0 = controller.step(r, x0);
    CHECK(u0(0) >= -1.0);
    CHECK(u0(0) <= 1.0);
    if (controller.last_status() != design::QPStatus::Success) {
        CHECK(u0(0) == doctest::Approx(0.0));
    }
    const auto u1 = controller.step(r, plant_step(x0, u0(0)));
    CHECK(u1(0) >= -1.0);
    CHECK(u1(0) <= 1.0);
    if (controller.last_status() != design::QPStatus::Success) {
        CHECK(u1(0) == doctest::Approx(u0(0)));
    }
}
