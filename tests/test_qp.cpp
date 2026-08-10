// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>

#include "damp/design/qp.hpp"
#include "damp/matlab.hpp"
#include "damp/matrix/matrix.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

// MATLAB® quadprog documentation example:
//   H = [1 -1; -1 2], f = [-2; -6], A = [1 1; -1 2; 2 1], b = [2; 2; 3]
//   -> x = [2/3; 4/3], fval = -74/9, constraints 1 and 2 active.
struct QuadprogExample {
    Matrix<2, 2> H{{1.0, -1.0}, {-1.0, 2.0}};
    ColVec<2>    f{-2.0, -6.0};
    Matrix<3, 2> A{{1.0, 1.0}, {-1.0, 2.0}, {2.0, 1.0}};
    ColVec<3>    b{2.0, 2.0, 3.0};
};

} // namespace

TEST_CASE("solve_qp reproduces the MATLAB® quadprog example") {
    const QuadprogExample qp;

    const auto res = design::solve_qp(qp.H, qp.f, qp.A, qp.b);
    REQUIRE(res.success);
    CHECK(res.status == design::QPStatus::Success);

    CHECK(res.x(0) == doctest::Approx(2.0 / 3.0).epsilon(1e-9));
    CHECK(res.x(1) == doctest::Approx(4.0 / 3.0).epsilon(1e-9));
    CHECK(res.objective == doctest::Approx(-74.0 / 9.0).epsilon(1e-9));

    // Multipliers (hand KKT solve): lambda = [28/9, 4/9, 0].
    CHECK(res.lambda(0) == doctest::Approx(28.0 / 9.0).epsilon(1e-9));
    CHECK(res.lambda(1) == doctest::Approx(4.0 / 9.0).epsilon(1e-9));
    CHECK(res.lambda(2) == doctest::Approx(0.0));
}

TEST_CASE("solve_qp satisfies the KKT conditions") {
    const QuadprogExample qp;
    const auto            res = design::solve_qp(qp.H, qp.f, qp.A, qp.b);
    REQUIRE(res.success);

    // Primal feasibility + dual feasibility + complementary slackness.
    for (size_t i = 0; i < 3; ++i) {
        double slack = -qp.b(i);
        for (size_t c = 0; c < 2; ++c) {
            slack += qp.A(i, c) * res.x(c);
        }
        CHECK(slack <= 1e-9);          // A x <= b
        CHECK(res.lambda(i) >= -1e-9); // lambda >= 0
        CHECK(res.lambda(i) * slack == doctest::Approx(0.0).epsilon(1e-9));
    }

    // Stationarity: H x + f + A' lambda = 0.
    for (size_t r = 0; r < 2; ++r) {
        double g = qp.f(r);
        for (size_t c = 0; c < 2; ++c) {
            g += qp.H(r, c) * res.x(c);
        }
        for (size_t i = 0; i < 3; ++i) {
            g += qp.A(i, r) * res.lambda(i);
        }
        CHECK(g == doctest::Approx(0.0).epsilon(1e-9));
    }
}

TEST_CASE("solve_qp is constexpr") {
    constexpr QuadprogExample qp;
    constexpr auto            res = design::solve_qp(qp.H, qp.f, qp.A, qp.b);
    static_assert(res.success);
    static_assert(res.x(0) - (2.0 / 3.0) < 1e-9 && res.x(0) - (2.0 / 3.0) > -1e-9);
    static_assert(res.x(1) - (4.0 / 3.0) < 1e-9 && res.x(1) - (4.0 / 3.0) > -1e-9);
}

TEST_CASE("solve_qp returns the unconstrained minimum when nothing binds") {
    const QuadprogExample qp;

    // All bounds at the unbounded sentinel: rows are skipped entirely.
    ColVec<3> b_open{};
    for (size_t i = 0; i < 3; ++i) {
        b_open(i) = design::unbounded_bound<double>();
    }
    const auto res = design::solve_qp(qp.H, qp.f, qp.A, b_open);
    REQUIRE(res.success);
    CHECK(res.iterations == 0);

    // x = -H^{-1} f = [10; 8]
    CHECK(res.x(0) == doctest::Approx(10.0));
    CHECK(res.x(1) == doctest::Approx(8.0));
    CHECK(res.lambda(0) == doctest::Approx(0.0));
}

TEST_CASE("solve_qp clamps against simple box constraints") {
    // min 1/2||x - (2,2)||^2 s.t. x <= (1,1): solution x = (1,1), lambda = (1,1).
    const auto         H = Matrix<2, 2>::identity();
    const ColVec<2>    f{-2.0, -2.0};
    const Matrix<2, 2> A{{1.0, 0.0}, {0.0, 1.0}};
    const ColVec<2>    b{1.0, 1.0};

    const auto res = design::solve_qp(H, f, A, b);
    REQUIRE(res.success);
    CHECK(res.x(0) == doctest::Approx(1.0));
    CHECK(res.x(1) == doctest::Approx(1.0));
    CHECK(res.lambda(0) == doctest::Approx(1.0));
    CHECK(res.lambda(1) == doctest::Approx(1.0));
    CHECK(res.objective == doctest::Approx(-3.0));
}

TEST_CASE("solve_qp detects infeasibility") {
    // x <= -1 and -x <= -1 (i.e. x >= 1) cannot both hold.
    const Matrix<1, 1> H{{1.0}};
    const ColVec<1>    f{0.0};
    const Matrix<2, 1> A{{1.0}, {-1.0}};
    const ColVec<2>    b{-1.0, -1.0};

    const auto res = design::solve_qp(H, f, A, b);
    CHECK(!res.success);
    CHECK(res.status == design::QPStatus::Infeasible);
}

TEST_CASE("solve_qp rejects an indefinite Hessian") {
    const Matrix<1, 1> H{{-1.0}};
    const ColVec<1>    f{1.0};
    const Matrix<1, 1> A{{1.0}};
    const ColVec<1>    b{1.0};

    const auto res = design::solve_qp(H, f, A, b);
    CHECK(!res.success);
    CHECK(res.status == design::QPStatus::NotPositiveDefinite);
}

TEST_CASE("solve_qp works in float") {
    const QuadprogExample qp;
    const auto            res = design::solve_qp(
        qp.H.as<float>(), ColVec<2, float>{qp.f.as<float>()},
        qp.A.as<float>(), ColVec<3, float>{qp.b.as<float>()}
    );
    REQUIRE(res.success);
    CHECK(res.x(0) == doctest::Approx(2.0F / 3.0F).epsilon(1e-5));
    CHECK(res.x(1) == doctest::Approx(4.0F / 3.0F).epsilon(1e-5));
}

TEST_CASE("equality-constrained solve_qp finds the minimum-norm allocation") {
    // min 1/2||x||^2 s.t. x1 + x2 + x3 = 1: x = (1/3, 1/3, 1/3), equality
    // multiplier -1/3 (exercises the negative-multiplier lower twin).
    const auto         H = Matrix<3, 3>::identity();
    const ColVec<3>    f{};
    const Matrix<1, 3> A{{1.0, 0.0, 0.0}};
    const ColVec<1>    b{design::unbounded_bound<double>()}; // no inequality binds
    const Matrix<1, 3> Aeq{{1.0, 1.0, 1.0}};
    const ColVec<1>    beq{1.0};

    const auto res = design::solve_qp(H, f, A, b, Aeq, beq);
    REQUIRE(res.success);
    for (size_t v = 0; v < 3; ++v) {
        CHECK(res.x(v) == doctest::Approx(1.0 / 3.0).epsilon(1e-9));
    }
    // Equality multiplier = upper twin - lower twin = -1/3.
    const double lambda_eq = res.lambda(1) - res.lambda(2);
    CHECK(lambda_eq == doctest::Approx(-1.0 / 3.0).epsilon(1e-9));
}

TEST_CASE("equality constraint composes with active inequalities") {
    // The quadprog example restricted to the line x1 = x2: along the line the
    // cost decreases until inequality rows 1 and 3 clamp it at x = (1, 1).
    const QuadprogExample qp;
    const Matrix<1, 2>    Aeq{{1.0, -1.0}};
    const ColVec<1>       beq{0.0};

    const auto res = design::solve_qp(qp.H, qp.f, qp.A, qp.b, Aeq, beq);
    REQUIRE(res.success);
    CHECK(res.x(0) == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(res.x(1) == doctest::Approx(1.0).epsilon(1e-9));
}

TEST_CASE("inconsistent equalities are reported infeasible") {
    const Matrix<2, 2> H{{1.0, 0.0}, {0.0, 1.0}};
    const ColVec<2>    f{};
    const Matrix<1, 2> A{{1.0, 0.0}};
    const ColVec<1>    b{design::unbounded_bound<double>()};
    const Matrix<2, 2> Aeq{{1.0, 1.0}, {1.0, 1.0}};
    const ColVec<2>    beq{1.0, 2.0}; // x1+x2 = 1 and = 2

    const auto res = design::solve_qp(H, f, A, b, Aeq, beq);
    CHECK(!res.success);
    CHECK(res.status == design::QPStatus::Infeasible);
}

TEST_CASE("equality-constrained solve_qp is constexpr") {
    constexpr auto res = [] {
        constexpr auto         H = Matrix<2, 2>::identity();
        constexpr ColVec<2>    f{-1.0, 0.0};
        constexpr Matrix<1, 2> A{{1.0, 0.0}};
        constexpr ColVec<1>    b{design::unbounded_bound<double>()};
        constexpr Matrix<1, 2> Aeq{{1.0, 1.0}};
        constexpr ColVec<1>    beq{1.0};
        return design::solve_qp(H, f, A, b, Aeq, beq);
    }();
    static_assert(res.success);
    // min 1/2||x||^2 - x1 s.t. x1+x2=1 -> x = (1, 0).
    static_assert(res.x(0) - 1.0 < 1e-9 && res.x(0) - 1.0 > -1e-9);
    static_assert(res.x(1) < 1e-9 && res.x(1) > -1e-9);
}

TEST_CASE("matlab::quadprog equality overload returns the same minimizer") {
    const auto         H = Matrix<3, 3>::identity();
    const ColVec<3>    f{};
    const Matrix<1, 3> A{{1.0, 0.0, 0.0}};
    const ColVec<1>    b{design::unbounded_bound<double>()};
    const Matrix<1, 3> Aeq{{1.0, 1.0, 1.0}};
    const ColVec<1>    beq{1.0};

    const auto x = matlab::quadprog(H, f, A, b, Aeq, beq);
    for (size_t v = 0; v < 3; ++v) {
        CHECK(x(v) == doctest::Approx(1.0 / 3.0).epsilon(1e-9));
    }
}

TEST_CASE("matlab::quadprog alias returns the same minimizer") {
    const QuadprogExample qp;
    const auto            x = matlab::quadprog(qp.H, qp.f, qp.A, qp.b);
    const auto            ref = design::solve_qp(qp.H, qp.f, qp.A, qp.b).x;
    CHECK(x(0) == doctest::Approx(ref(0)));
    CHECK(x(1) == doctest::Approx(ref(1)));
}
