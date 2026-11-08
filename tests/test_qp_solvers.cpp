// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>
#include <limits>

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
    .C = {{1.0, 0.0}}, // the state measured, the position controlled
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
    weights.Rdu = 0.05;
    design::MPCConstraints<1, 1> limits{};
    limits.u_lim = {
        .lower = -0.5,
        .upper = 0.5,
    };
    limits.du_lim = {
        .lower = -0.2,
        .upper = 0.2,
    };

    const auto art = design::mpc<15, 5>(double_integrator, weights, limits);
    REQUIRE(art.success);

    MPC warm{art}; // WarmStartActiveSetSolver is the default policy
    MPC cold{art, design::ActiveSetSolver{}};

    ColVec<2>       xw{};
    ColVec<2>       xc{};
    const ColVec<1> r{1.0};
    size_t          warm_settled = 0;
    size_t          cold_settled = 0;
    for (size_t k = 0; k < 200; ++k) {
        const auto uw = warm.step(r, xw);
        const auto uc = cold.step(r, xc);
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
    weights.Rdu = 0.05;
    design::MPCConstraints<1, 1> limits{};
    limits.u_lim = {
        .lower = -0.5,
        .upper = 0.5,
    };
    limits.max_qp_iterations = 20000;

    const auto art = design::mpc<15, 5>(double_integrator, weights, limits);
    REQUIRE(art.success);

    using Ctl = decltype(MPC{art});
    MPC exact{art};
    MPC admm{art, design::AdmmSolver<Ctl::NZ1, Ctl::NI1, double>{}};

    ColVec<2>       xe{};
    ColVec<2>       xa{};
    const ColVec<1> r{1.0};
    for (size_t k = 0; k < 120; ++k) {
        const auto ue = exact.step(r, xe);
        const auto ua = admm.step(r, xa);
        CHECK(ua(0) == doctest::Approx(ue(0)).epsilon(1e-4));
        xe = plant_step(xe, ue(0));
        xa = plant_step(xa, ua(0));
    }
    CHECK(xa(0) == doctest::Approx(1.0).epsilon(0.02));
}

TEST_CASE("ADMM certifies a primal-infeasible QP from its dual iterate difference") {
    // The quadprog example plus x1 + x2 >= 3, against row 0's x1 + x2 <= 2.
    const QuadprogExample qp;
    const Matrix<4, 2>    A{{1.0, 1.0}, {-1.0, 2.0}, {2.0, 1.0}, {-1.0, -1.0}};
    const ColVec<4>       b{2.0, 2.0, 3.0, -3.0};

    design::AdmmSolver<2, 4> solver{};
    const auto               res = solver(qp.H, qp.f, A, b, 4000);
    CHECK(!res.success);
    CHECK(res.status == design::QPStatus::Infeasible);
    CHECK(res.iterations < 200); // certified long before the 4000-iteration budget
    CHECK(design::solve_qp(qp.H, qp.f, A, b).status == design::QPStatus::Infeasible);

    // The certificate clears the warm start: the next feasible solve (same A,
    // row 3 disabled) runs as on a fresh solver: the same iterations, and the same
    // iterate up to rounding (the two calls are inlined at different sites, which
    // -ffast-math and FMA contraction may round differently).
    const ColVec<4> b_ok{2.0, 2.0, 3.0, design::unbounded_bound<double>()};
    const auto      after = solver(qp.H, qp.f, A, b_ok, 4000);
    const auto      fresh = design::solve_qp_admm(qp.H, qp.f, A, b_ok);
    REQUIRE(after.success);
    CHECK(after.iterations == fresh.iterations);
    const double ulps = (after.x - fresh.x).norm() / (eps<double>() * fresh.x.norm());
    INFO("after vs fresh: ", ulps, " ulps of |x|"); // GCC 15: 3
    CHECK(ulps <= 16.0);
    CHECK(after.x(0) == doctest::Approx(2.0 / 3.0).epsilon(1e-6));
}

TEST_CASE("ADMM primal certificate: contradictory bounds, and one infeasible variable beside an active one") {
    const Matrix<1, 1> H1{{1.0}};
    const ColVec<1>    f1{0.0};
    const Matrix<2, 1> A1{{1.0}, {-1.0}};
    const ColVec<2>    b1{-1.0, -1.0}; // x <= -1 and x >= 1
    const auto         one = design::solve_qp_admm(H1, f1, A1, b1);
    CHECK(one.status == design::QPStatus::Infeasible);
    CHECK(one.iterations < 10);

    // x1 infeasible as above; x2 pushed onto its bound x2 <= 1, whose multiplier settles.
    const auto         H2 = Matrix<2, 2>::identity();
    const ColVec<2>    f2{0.0, -5.0};
    const Matrix<3, 2> A2{{1.0, 0.0}, {-1.0, 0.0}, {0.0, 1.0}};
    const ColVec<3>    b2{-1.0, -1.0, 1.0};
    const auto         two = design::solve_qp_admm(H2, f2, A2, b2);
    CHECK(two.status == design::QPStatus::Infeasible);
    CHECK(two.iterations < 400);

    // Scaling a row (aᵢ, bᵢ) or the whole problem leaves the certificate in place.
    const Matrix<2, 1> A1s{{1.0e3}, {-1.0}};
    const ColVec<2>    b1s{-1.0e3, -1.0};
    CHECK(design::solve_qp_admm(H1, f1, A1s, b1s).status == design::QPStatus::Infeasible);
    CHECK(design::solve_qp_admm(H2, f2, A2 * 1.0e4, b2 * 1.0e4).status == design::QPStatus::Infeasible);
}

TEST_CASE("ADMM certifies a dual-infeasible (unbounded) QP from its primal iterate difference") {
    // LP min −x1 with only x2 <= 1: x1 → +∞.
    const Matrix<2, 2> H0{};
    const ColVec<2>    f_lp{-1.0, 0.0};
    const Matrix<1, 2> A_lp{{0.0, 1.0}};
    const ColVec<1>    b_lp{1.0};
    const auto         lp = design::solve_qp_admm(H0, f_lp, A_lp, b_lp);
    CHECK(!lp.success);
    CHECK(lp.status == design::QPStatus::Unbounded);

    // Semidefinite H with a coupled null direction: ½(x1 − x2)² − x1 − 2x2 with
    // x1 − x2 <= 1 and x1 >= 0 decreases without bound along (1, 1).
    const Matrix<2, 2> H{{1.0, -1.0}, {-1.0, 1.0}};
    const ColVec<2>    f{-1.0, -2.0};
    const Matrix<2, 2> A{{1.0, -1.0}, {-1.0, 0.0}};
    const ColVec<2>    b{1.0, 0.0};
    const auto         res = design::solve_qp_admm(H, f, A, b);
    CHECK(res.status == design::QPStatus::Unbounded);
    CHECK(res.iterations < 200);
    // The factoring solvers cannot factor this H.
    CHECK(design::solve_qp(H, f, A, b).status == design::QPStatus::NotPositiveDefinite);
    CHECK(design::solve_qp_interior_point(H, f, A, b).status == design::QPStatus::NotPositiveDefinite);
}

TEST_CASE("ADMM feasible results do not depend on the certificate tolerance") {
    // A certificate is only tested after the convergence test fails, so on a
    // feasible problem it changes nothing: a loose tolerance gives the same bits.
    design::AdmmSettings<double> loose{};
    loose.eps_infeasible = 1.0e-2;

    const QuadprogExample qp;
    const auto            tight_qp = design::solve_qp_admm(qp.H, qp.f, qp.A, qp.b);
    const auto            loose_qp = design::solve_qp_admm(qp.H, qp.f, qp.A, qp.b, 4000, loose);
    REQUIRE(tight_qp.success);
    REQUIRE(loose_qp.success);
    CHECK(loose_qp.iterations == tight_qp.iterations);
    CHECK(loose_qp.x(0) == tight_qp.x(0));
    CHECK(loose_qp.x(1) == tight_qp.x(1));

    // A bounded LP (H = 0, which the factoring solvers refuse) still solves: min −x on [0, 2].
    const Matrix<1, 1> H0{};
    const ColVec<1>    f{-1.0};
    const Matrix<2, 1> A{{1.0}, {-1.0}};
    const ColVec<2>    b{2.0, 0.0};
    const auto         lp = design::solve_qp_admm(H0, f, A, b);
    const auto         lp_loose = design::solve_qp_admm(H0, f, A, b, 4000, loose);
    REQUIRE(lp.success);
    CHECK(lp.x(0) == doctest::Approx(2.0).epsilon(1e-6));
    CHECK(lp_loose.x(0) == lp.x(0));
    CHECK(lp_loose.iterations == lp.iterations);
}

TEST_CASE("ADMM reports an MPC's infeasible constraint set like the active-set default") {
    // The held input 2.0 lies outside |u| <= 0.5 and |Δu| <= 0.1 cannot bring it
    // back in one move: the hard boxes admit no move. Both solvers say so, the
    // command holds the input clamped to the box, and the next tick (held input
    // back in range) solves again.
    design::MPCWeights<1, 1> weights{};
    weights.Rdu = 0.05;
    design::MPCConstraints<1, 1> limits{};
    limits.u_lim = {
        .lower = -0.5,
        .upper = 0.5,
    };
    limits.du_lim = {
        .lower = -0.1,
        .upper = 0.1,
    };
    limits.max_qp_iterations = 20000;

    const auto art = design::mpc<15, 5>(double_integrator, weights, limits);
    REQUIRE(art.success);

    using Ctl = decltype(MPC{art});
    MPC exact{art};
    MPC admm{art, design::AdmmSolver<Ctl::NZ1, Ctl::NI1, double>{}};

    const ColVec<2> x{};
    const ColVec<1> r{1.0};
    exact.back_calculate(ColVec<1>{2.0});
    admm.back_calculate(ColVec<1>{2.0});
    const auto ue = exact.step(r, x);
    const auto ua = admm.step(r, x);
    CHECK(exact.last_status() == design::QPStatus::Infeasible);
    CHECK(admm.last_status() == design::QPStatus::Infeasible);
    CHECK(admm.last_iterations() < 1000); // was the whole 20000 budget, reported as MaxIterations
    CHECK(ue(0) == doctest::Approx(0.5));
    CHECK(ua(0) == doctest::Approx(0.5));

    exact.back_calculate(ColVec<1>{0.3});
    admm.back_calculate(ColVec<1>{0.3});
    const auto ue2 = exact.step(r, x);
    const auto ua2 = admm.step(r, x);
    CHECK(exact.last_status() == design::QPStatus::Success);
    CHECK(admm.last_status() == design::QPStatus::Success);
    CHECK(ua2(0) == doctest::Approx(ue2(0)).epsilon(1e-6));
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
    weights.Rdu = 0.05;
    design::MPCConstraints<1, 1> limits{};
    limits.u_lim = {
        .lower = -0.5,
        .upper = 0.5,
    };

    const auto art = design::mpc<15, 5>(double_integrator, weights, limits);
    REQUIRE(art.success);

    MPC exact{art};
    MPC ip{art, design::InteriorPointSolver{}};

    ColVec<2>       xe{};
    ColVec<2>       xi{};
    const ColVec<1> r{1.0};
    for (size_t k = 0; k < 120; ++k) {
        const auto ue = exact.step(r, xe);
        const auto ui = ip.step(r, xi);
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
    constexpr auto admm_infeasible = design::solve_qp_admm(Matrix<1, 1>{{1.0}}, ColVec<1>{0.0}, Matrix<2, 1>{{1.0}, {-1.0}}, ColVec<2>{-1.0, -1.0});
    static_assert(admm_infeasible.status == design::QPStatus::Infeasible);

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

TEST_CASE("interior-point QP is invariant to scaling the cost") {
    for (const double s : {1e-6, 1.0, 1e9}) {
        const Matrix<2, 2> H{{s, -s}, {-s, 2.0 * s}};
        const ColVec<2>    f{-2.0 * s, -6.0 * s};
        const Matrix<3, 2> A{{1.0, 1.0}, {-1.0, 2.0}, {2.0, 1.0}};
        const ColVec<3>    b{2.0, 2.0, 3.0};
        const auto         res = design::solve_qp_interior_point(H, f, A, b);
        CAPTURE(s);
        REQUIRE(res.success);
        CHECK(res.x(0) == doctest::Approx(2.0 / 3.0).epsilon(1e-6));
        CHECK(res.x(1) == doctest::Approx(4.0 / 3.0).epsilon(1e-6));
    }
}

TEST_CASE("ADMM reports a NaN input as NonFinite and its warm start recovers on finite data") {
    //! The residual maxima dropped the NaN: Success at iteration 1, and x = NaN on every later
    //! call until reset().
    const QuadprogExample        qp;
    design::AdmmSolver<2, 3>     solver;
    const ColVec<3>              b_nan{2.0, std::numeric_limits<double>::quiet_NaN(), 3.0};
    const design::QPResult<2, 3> bad = solver(qp.H, qp.f, qp.A, b_nan, 1000);
    CHECK(bad.status == design::QPStatus::NonFinite);
    CHECK_FALSE(bad.success);
    const design::QPResult<2, 3> good = solver(qp.H, qp.f, qp.A, qp.b, 1000);
    REQUIRE(good.success);
    CHECK(good.x(0) == doctest::Approx(2.0 / 3.0).epsilon(1e-6));
    CHECK(good.x(1) == doctest::Approx(4.0 / 3.0).epsilon(1e-6));
}

TEST_CASE("float ADMM reaches its default tolerances") {
    //! eps_abs = eps_rel = 1e-8 sat below float's ε: every float solve ran the budget out.
    const Matrix<2, 2, float>       H{{1.0F, -1.0F}, {-1.0F, 2.0F}};
    const ColVec<2, float>          f{-2.0F, -6.0F};
    const Matrix<3, 2, float>       A{{1.0F, 1.0F}, {-1.0F, 2.0F}, {2.0F, 1.0F}};
    const ColVec<3, float>          b{2.0F, 2.0F, 3.0F};
    design::AdmmSolver<2, 3, float> solver;
    const auto                      res = solver(H, f, A, b, 1000);
    REQUIRE(res.success);
    CHECK(res.iterations < 1000);
    CHECK(res.x(0) == doctest::Approx(2.0 / 3.0).epsilon(1e-3));
    CHECK(res.x(1) == doctest::Approx(4.0 / 3.0).epsilon(1e-3));
}

TEST_CASE("float MIQP solves at 1e4, where a bound sits inside the relaxation's slack tolerance") {
    //! It returned MaxIterations with no incumbent: x₁ ≥ 10000 at x₁ = 9999.6 lies within the
    //! active-set solver's (NV+1)·64·ε of 2e4, so the up child never bound and was re-branched.
    const float               t = 9999.6F;
    const Matrix<2, 2, float> H = Matrix<2, 2, float>::identity();
    const ColVec<2, float>    f{-t, -(t + 0.3F)};
    const Matrix<1, 2, float> A{{1.0F, 1.0F}};
    const ColVec<1, float>    b{(2.0F * t) + 0.1F};
    const auto                mi = design::solve_miqp(H, f, A, b, damp::array<bool, 2>{true, true});
    REQUIRE(mi.success);
    CHECK(mi.status == design::QPStatus::Success);
    CHECK(mi.x(0) == 9999.0F);
    CHECK(mi.x(1) == 10000.0F);
}

TEST_CASE("float MIQP does not take a half-integer at 1.5e5 as integral") {
    //! The relative integrality term 8·(NV+1)·ε·|x| exceeded 0.5 past 1.05e5 for NV = 4, so the
    //! root relaxation x₀ = 150000.5 returned as Success.
    const Matrix<4, 4, float> H = Matrix<4, 4, float>::identity();
    const ColVec<4, float>    f{-150000.5F, 0.0F, 0.0F, 0.0F};
    const Matrix<1, 4, float> A{};
    const ColVec<1, float>    b{design::unbounded_bound<float>()};
    const auto                mi = design::solve_miqp(H, f, A, b, damp::array<bool, 4>{true, true, true, true});
    REQUIRE(mi.success);
    CHECK((mi.x(0) == 150000.0F || mi.x(0) == 150001.0F));
}

TEST_CASE("MIQP refuses non-finite data") {
    //! A NaN relaxation reached floor_scalar's long long cast (undefined behavior).
    const Matrix<2, 2> H = Matrix<2, 2>::identity();
    const ColVec<2>    f{std::numeric_limits<double>::quiet_NaN(), -1.0};
    const Matrix<1, 2> A{{1.0, 1.0}};
    const auto         mi = design::solve_miqp(H, f, A, ColVec<1>{2.0}, damp::array<bool, 2>{true, true});
    CHECK(mi.status == design::QPStatus::NonFinite);
    CHECK_FALSE(mi.success);
}
