// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/design/qp.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

// MATLAB® quadprog documentation example: x* = [2/3; 4/3], fval = -74/9,
// lambda = [28/9, 4/9, 0]. Shared known answer across all solvers.
struct QuadprogExample {
    Matrix<2, 2> H{{1.0, -1.0}, {-1.0, 2.0}};
    ColVec<2>    f{-2.0, -6.0};
    Matrix<3, 2> A{{1.0, 1.0}, {-1.0, 2.0}, {2.0, 1.0}};
    ColVec<3>    b{2.0, 2.0, 3.0};
};

constexpr double Ts = 0.1;

constexpr StateSpace<2, 1, 1> double_integrator{
    .A = {{1.0, Ts}, {0.0, 1.0}},
    .B = {{0.5 * Ts * Ts}, {Ts}},
    .C = {{1.0, 0.0}},
    .Ts = Ts,
};

constexpr ColVec<2> plant_step(const ColVec<2>& x, double u) {
    return ColVec<2>{
        x(0) + (Ts * x(1)) + (0.5 * Ts * Ts * u),
        x(1) + (Ts * u),
    };
}

} // namespace

// --- Warm-started active set ---------------------------------------------------

TEST_CASE("warm-start policy matches the stateless solve exactly") {
    const QuadprogExample                  qp;
    design::WarmStartActiveSetSolver<2, 3> solver{};

    const auto cold_ref = design::solve_qp(qp.H, qp.f, qp.A, qp.b);
    const auto first = solver(qp.H, qp.f, qp.A, qp.b, 100);
    REQUIRE(first.success);
    CHECK(first.x(0) == doctest::Approx(cold_ref.x(0)).epsilon(1e-12));
    CHECK(first.x(1) == doctest::Approx(cold_ref.x(1)).epsilon(1e-12));
    CHECK(first.lambda(0) == doctest::Approx(cold_ref.lambda(0)).epsilon(1e-12));

    // Second solve runs off the harvested hint: same answer, no extra work.
    const auto second = solver(qp.H, qp.f, qp.A, qp.b, 100);
    REQUIRE(second.success);
    CHECK(second.x(0) == doctest::Approx(first.x(0)).epsilon(1e-12));
    CHECK(second.iterations <= first.iterations);
}

TEST_CASE("warm-start solver stays correct across a shifting constraint set") {
    // Sweep the bounds the way an MPC tick does (same H/A, moving b): every
    // solve must satisfy the KKT conditions even as the active set changes and
    // the hint goes stale.
    const QuadprogExample                  qp;
    design::WarmStartActiveSetSolver<2, 3> solver{};

    for (int step = 0; step <= 20; ++step) {
        const double    shift = 0.2 * step; // loosens until nothing binds
        const ColVec<3> b{2.0 + shift, 2.0 + shift, 3.0 + shift};

        const auto res = solver(qp.H, qp.f, qp.A, b, 100);
        REQUIRE(res.success);

        // Primal feasibility + complementary slackness.
        for (size_t i = 0; i < 3; ++i) {
            double s = -b(i);
            for (size_t c = 0; c < 2; ++c) {
                s += qp.A(i, c) * res.x(c);
            }
            CHECK(s <= 1e-9);
            CHECK(res.lambda(i) * s == doctest::Approx(0.0).epsilon(1e-8));
        }
        // Stationarity.
        for (size_t r = 0; r < 2; ++r) {
            double g = qp.f(r);
            for (size_t c = 0; c < 2; ++c) {
                g += qp.H(r, c) * res.x(c);
            }
            for (size_t i = 0; i < 3; ++i) {
                g += qp.A(i, r) * res.lambda(i);
            }
            CHECK(g == doctest::Approx(0.0).epsilon(1e-8));
        }
    }
}

TEST_CASE("warm-start solver is constexpr") {
    constexpr auto twice = [] {
        constexpr QuadprogExample              qp;
        design::WarmStartActiveSetSolver<2, 3> solver{};
        const auto                             a = solver(qp.H, qp.f, qp.A, qp.b, 100);
        const auto                             b2 = solver(qp.H, qp.f, qp.A, qp.b, 100);
        return damp::pair{a, b2};
    }();
    static_assert(twice.first.success);
    static_assert(twice.second.success);
    static_assert(twice.second.iterations <= twice.first.iterations);
}

TEST_CASE("warm-start default MPC matches the stateless policy tick for tick") {
    design::MPCWeights<1, 1> weights{};
    weights.Rdu = Matrix<1, 1>{{0.05}};
    design::MPCConstraints<1, 1> limits{};
    limits.u_min = ColVec<1>{-0.5};
    limits.u_max = ColVec<1>{0.5};
    limits.du_min = ColVec<1>{-0.2};
    limits.du_max = ColVec<1>{0.2};

    const auto art = design::state_mpc<15, 5>(double_integrator, weights, limits);
    REQUIRE(art.success);

    MPC                                                     warm{art}; // WarmStartActiveSetSolver is the default policy
    MPC<2, 1, 1, 15, 5, 0, double, design::ActiveSetSolver> cold{art};

    ColVec<2>       xw{};
    ColVec<2>       xc{};
    const ColVec<1> r{1.0};
    size_t          warm_settled = 0;
    size_t          cold_settled = 0;
    for (size_t k = 0; k < 200; ++k) {
        const auto uw = warm.control(r, xw);
        const auto uc = cold.control(r, xc);
        CHECK(uw(0) == doctest::Approx(uc(0)).epsilon(1e-9));
        // Transient pick orders may differ either way; the guaranteed win is at
        // steady state, where the hint equals the active set.
        if (k >= 150) {
            warm_settled += warm.last_iterations();
            cold_settled += cold.last_iterations();
        }
        xw = plant_step(xw, uw(0));
        xc = plant_step(xc, uc(0));
    }
    CHECK(warm_settled <= cold_settled);
    CHECK(xw(0) == doctest::Approx(1.0).epsilon(0.02));
}

// --- ADMM --------------------------------------------------------------------

TEST_CASE("ADMM reproduces the quadprog example") {
    const QuadprogExample qp;
    const auto            res = design::solve_qp_admm(qp.H, qp.f, qp.A, qp.b);
    REQUIRE(res.success);
    CHECK(res.x(0) == doctest::Approx(2.0 / 3.0).epsilon(1e-6));
    CHECK(res.x(1) == doctest::Approx(4.0 / 3.0).epsilon(1e-6));
    CHECK(res.objective == doctest::Approx(-74.0 / 9.0).epsilon(1e-6));
    CHECK(res.lambda(0) == doctest::Approx(28.0 / 9.0).epsilon(1e-4));
    CHECK(res.lambda(1) == doctest::Approx(4.0 / 9.0).epsilon(1e-4));
}

TEST_CASE("ADMM warm start cuts the iteration count on a repeated solve") {
    const QuadprogExample    qp;
    design::AdmmSolver<2, 3> solver{};

    const auto cold = solver(qp.H, qp.f, qp.A, qp.b, 4000);
    REQUIRE(cold.success);
    const auto warm = solver(qp.H, qp.f, qp.A, qp.b, 4000);
    REQUIRE(warm.success);
    CHECK(warm.iterations < cold.iterations);
    CHECK(warm.x(0) == doctest::Approx(cold.x(0)).epsilon(1e-8));
}

TEST_CASE("ADMM solver policy drives the MPC like the active-set default") {
    design::MPCWeights<1, 1> weights{};
    weights.Rdu = Matrix<1, 1>{{0.05}};
    design::MPCConstraints<1, 1> limits{};
    limits.u_min = ColVec<1>{-0.5};
    limits.u_max = ColVec<1>{0.5};

    const auto art = design::state_mpc<15, 5>(double_integrator, weights, limits, 20000);
    REQUIRE(art.success);

    using Ctl = MPC<2, 1, 1, 15, 5, 0, double>;
    MPC                                                                            exact{art};
    MPC<2, 1, 1, 15, 5, 0, double, design::AdmmSolver<Ctl::NZ1, Ctl::NI1, double>> admm{art};

    ColVec<2>       xe{};
    ColVec<2>       xa{};
    const ColVec<1> r{1.0};
    for (size_t k = 0; k < 120; ++k) {
        const auto ue = exact.control(r, xe);
        const auto ua = admm.control(r, xa);
        CHECK(ua(0) == doctest::Approx(ue(0)).epsilon(1e-4));
        xe = plant_step(xe, ue(0));
        xa = plant_step(xa, ua(0));
    }
    CHECK(xa(0) == doctest::Approx(1.0).epsilon(0.02));
}

// --- Interior point ------------------------------------------------------------

TEST_CASE("interior point reproduces the quadprog example") {
    const QuadprogExample qp;
    const auto            res = design::solve_qp_interior_point(qp.H, qp.f, qp.A, qp.b);
    REQUIRE(res.success);
    CHECK(res.x(0) == doctest::Approx(2.0 / 3.0).epsilon(1e-6));
    CHECK(res.x(1) == doctest::Approx(4.0 / 3.0).epsilon(1e-6));
    CHECK(res.lambda(0) == doctest::Approx(28.0 / 9.0).epsilon(1e-4));
    CHECK(res.lambda(1) == doctest::Approx(4.0 / 9.0).epsilon(1e-4));
    CHECK(res.lambda(2) == doctest::Approx(0.0).epsilon(1e-6));
}

TEST_CASE("interior point returns the unconstrained minimum on sentinel rows") {
    const QuadprogExample qp;
    ColVec<3>             b_open{};
    for (size_t i = 0; i < 3; ++i) {
        b_open(i) = design::unbounded_bound<double>();
    }
    const auto res = design::solve_qp_interior_point(qp.H, qp.f, qp.A, b_open);
    REQUIRE(res.success);
    CHECK(res.x(0) == doctest::Approx(10.0));
    CHECK(res.x(1) == doctest::Approx(8.0));
}

TEST_CASE("interior point flags an infeasible constraint set") {
    const Matrix<1, 1> H{{1.0}};
    const ColVec<1>    f{0.0};
    const Matrix<2, 1> A{{1.0}, {-1.0}};
    const ColVec<2>    b{-1.0, -1.0}; // x <= -1 and x >= 1

    const auto res = design::solve_qp_interior_point(H, f, A, b, 100);
    CHECK(!res.success);
    CHECK(res.status != design::QPStatus::Success);
}

TEST_CASE("interior-point solver policy drives the MPC like the active-set default") {
    design::MPCWeights<1, 1> weights{};
    weights.Rdu = Matrix<1, 1>{{0.05}};
    design::MPCConstraints<1, 1> limits{};
    limits.u_min = ColVec<1>{-0.5};
    limits.u_max = ColVec<1>{0.5};

    const auto art = design::state_mpc<15, 5>(double_integrator, weights, limits);
    REQUIRE(art.success);

    MPC                                                         exact{art};
    MPC<2, 1, 1, 15, 5, 0, double, design::InteriorPointSolver> ip{art};

    ColVec<2>       xe{};
    ColVec<2>       xi{};
    const ColVec<1> r{1.0};
    for (size_t k = 0; k < 120; ++k) {
        const auto ue = exact.control(r, xe);
        const auto ui = ip.control(r, xi);
        CHECK(ui(0) == doctest::Approx(ue(0)).epsilon(1e-5));
        xe = plant_step(xe, ue(0));
        xi = plant_step(xi, ui(0));
    }
    CHECK(xi(0) == doctest::Approx(1.0).epsilon(0.02));
}

// --- MIQP ---------------------------------------------------------------------

TEST_CASE("MIQP matches the continuous solve when the relaxation is integral") {
    // Unconstrained minimum at exactly (3, -2): already integral, zero branching.
    const auto      H = Matrix<2, 2>::identity();
    const ColVec<2> f{-3.0, 2.0};
    Matrix<1, 2>    A{{1.0, 0.0}};
    ColVec<1>       b{design::unbounded_bound<double>()};

    damp::array<bool, 2> all_int{};
    all_int[0] = true;
    all_int[1] = true;

    const auto mi = design::solve_miqp(H, f, A, b, all_int);
    REQUIRE(mi.success);
    CHECK(mi.status == design::QPStatus::Success);
    CHECK(mi.x(0) == doctest::Approx(3.0));
    CHECK(mi.x(1) == doctest::Approx(-2.0));
    CHECK(mi.iterations == 1); // root relaxation only
}

TEST_CASE("MIQP branches to the best lattice point under a coupling constraint") {
    // min ||x - (2.4, 1.6)||^2 with x integer and x1 + x2 <= 3.
    // Relaxed optimum (1.9, 1.1) is fractional; best integer point is (2, 1).
    const auto         H = Matrix<2, 2>::identity();
    const ColVec<2>    f{-2.4, -1.6};
    const Matrix<1, 2> A{{1.0, 1.0}};
    const ColVec<1>    b{3.0};

    damp::array<bool, 2> all_int{};
    all_int[0] = true;
    all_int[1] = true;

    const auto mi = design::solve_miqp(H, f, A, b, all_int);
    REQUIRE(mi.success);
    CHECK(mi.status == design::QPStatus::Success);
    CHECK(mi.x(0) == doctest::Approx(2.0));
    CHECK(mi.x(1) == doctest::Approx(1.0));
    CHECK(mi.iterations > 1); // branching actually happened
}

TEST_CASE("MIQP handles mixed integer and continuous variables") {
    // x0 integer, x1 continuous: min (x0-0.6)^2 + (x1-0.7)^2 -> (1, 0.7).
    const auto      H = Matrix<2, 2>::identity();
    const ColVec<2> f{-0.6, -0.7};
    Matrix<1, 2>    A{{1.0, 0.0}};
    ColVec<1>       b{design::unbounded_bound<double>()};

    damp::array<bool, 2> mask{};
    mask[0] = true;

    const auto mi = design::solve_miqp(H, f, A, b, mask);
    REQUIRE(mi.success);
    CHECK(mi.x(0) == doctest::Approx(1.0));
    CHECK(mi.x(1) == doctest::Approx(0.7));
}

TEST_CASE("MIQP proves infeasibility when no lattice point fits") {
    // 0.2 <= x <= 0.8 contains no integer.
    const Matrix<1, 1> H{{1.0}};
    const ColVec<1>    f{0.0};
    const Matrix<2, 1> A{{1.0}, {-1.0}};
    const ColVec<2>    b{0.8, -0.2};

    damp::array<bool, 1> mask{};
    mask[0] = true;

    const auto mi = design::solve_miqp(H, f, A, b, mask);
    CHECK(!mi.success);
    CHECK(mi.status == design::QPStatus::Infeasible);
}

// --- constexpr coverage ---------------------------------------------------------

TEST_CASE("alternative solvers are constexpr") {
    constexpr QuadprogExample qp;

    constexpr auto ip = design::solve_qp_interior_point(qp.H, qp.f, qp.A, qp.b);
    static_assert(ip.success);
    static_assert(ip.x(0) - (2.0 / 3.0) < 1e-6 && ip.x(0) - (2.0 / 3.0) > -1e-6);

    constexpr auto admm = design::solve_qp_admm(qp.H, qp.f, qp.A, qp.b);
    static_assert(admm.success);

    constexpr auto mi = [] {
        constexpr auto         H = Matrix<1, 1>{{1.0}};
        constexpr ColVec<1>    f{-1.4};
        constexpr Matrix<1, 1> A{{1.0}};
        constexpr ColVec<1>    b{design::unbounded_bound<double>()};
        damp::array<bool, 1>   mask{};
        mask[0] = true;
        return design::solve_miqp(H, f, A, b, mask);
    }();
    static_assert(mi.success);
    static_assert(mi.x(0) - 1.0 < 1e-9 && mi.x(0) - 1.0 > -1e-9);
}
