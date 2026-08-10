// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>

#include "damp/backend.hpp"
#include "damp/design/lyapunov.hpp"
#include "damp/design/stability.hpp"
#include "damp/matlab.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

// Residual ‖A X + X Aᵀ + Q‖∞ for the continuous Lyapunov equation.
template<size_t N>
double cont_residual(const Matrix<N, N>& A, const Matrix<N, N>& X, const Matrix<N, N>& Q) {
    const Matrix<N, N> R = A * X + X * A.transpose() + Q;
    double             m = 0.0;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            m = damp::max(m, std::abs(R(i, j)));
        }
    }
    return m;
}

template<size_t N>
double disc_residual(const Matrix<N, N>& A, const Matrix<N, N>& X, const Matrix<N, N>& Q) {
    const Matrix<N, N> R = A * X * A.transpose() - X + Q;
    double             m = 0.0;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            m = damp::max(m, std::abs(R(i, j)));
        }
    }
    return m;
}

} // namespace

TEST_CASE("lyap solves the continuous Lyapunov equation") {
    const Matrix<2, 2> A{{-2.0, 1.0}, {0.0, -3.0}};
    const Matrix<2, 2> Q{{3.0, 1.0}, {1.0, 2.0}};

    const auto X = lyap(A, Q);
    REQUIRE(X.has_value());
    CHECK(cont_residual(A, X.value(), Q) < 1e-12);
}

TEST_CASE("dlyap solves the discrete Lyapunov equation") {
    const Matrix<2, 2> A{{0.5, 0.1}, {0.0, 0.3}};
    const Matrix<2, 2> Q{{3.0, 1.0}, {1.0, 2.0}};

    const auto X = dlyap(A, Q);
    REQUIRE(X.has_value());
    CHECK(disc_residual(A, X.value(), Q) < 1e-12);
}

TEST_CASE("lyap reports no unique solution when A and -A share an eigenvalue") {
    // Eigenvalues {+1, -1}: λ_i + λ_j = 0, so the operator is singular.
    const Matrix<2, 2> A{{1.0, 0.0}, {0.0, -1.0}};
    const Matrix<2, 2> Q = Matrix<2, 2>::identity();
    CHECK_FALSE(lyap(A, Q).has_value());
}

TEST_CASE("gram gives the controllability/observability Gramians") {
    // Stable continuous SISO plant.
    StateSpace<2, 1, 1> sys;
    sys.A = Matrix<2, 2>{{-2.0, 1.0}, {0.0, -3.0}};
    sys.B = ColVec<2>{1.0, 1.0};
    sys.C = RowVec<2>{1.0, 0.0};
    sys.Ts = 0.0;

    const auto Wc = matlab::gram(sys, 'c');
    const auto Wo = matlab::gram(sys, 'o');
    REQUIRE(Wc.has_value());
    REQUIRE(Wo.has_value());

    // Wc must satisfy A Wc + Wc Aᵀ + B Bᵀ = 0; symmetric positive definite.
    CHECK(cont_residual(sys.A, Wc.value(), sys.B * sys.B.transpose()) < 1e-12);
    CHECK(Wc.value()(0, 0) > 0.0);
    CHECK(Wo.value()(0, 0) > 0.0);
    // Controllable + observable ⇒ both Gramians full rank (positive determinant 2×2).
    const auto det = [](const Matrix<2, 2>& M) { return (M(0, 0) * M(1, 1)) - (M(0, 1) * M(1, 0)); };
    CHECK(det(Wc.value()) > 0.0);
    CHECK(det(Wo.value()) > 0.0);
}

TEST_CASE("H2 and Hinf norms of a first-order lag") {
    // G(s) = 1/(s+1): A=-1, B=1, C=1, D=0.
    // ‖G‖₂ = sqrt(1/2) ≈ 0.70711 ; ‖G‖∞ = |G(0)| = 1.
    StateSpace<1, 1, 1> sys;
    sys.A = Matrix<1, 1>{{-1.0}};
    sys.B = Matrix<1, 1>{{1.0}};
    sys.C = Matrix<1, 1>{{1.0}};
    sys.Ts = 0.0;

    const auto h2 = matlab::norm(sys);
    REQUIRE(h2.has_value());
    CHECK(h2.value() == doctest::Approx(0.7071067812).epsilon(1e-6));
    CHECK(matlab::hinfnorm(sys) == doctest::Approx(1.0).epsilon(1e-3));
}

TEST_CASE("Hinf norm captures a lightly-damped resonant peak") {
    // G(s) = wn² / (s² + 2ζwn s + wn²), wn=10, ζ=0.05.
    // DC gain 1; resonant peak ≈ 1/(2ζ√(1−ζ²)) ≈ 10.01 near ω≈wn.
    const double        wn = 10.0;
    const double        zeta = 0.05;
    StateSpace<2, 1, 1> sys; // controllable canonical form
    sys.A = Matrix<2, 2>{{0.0, 1.0}, {-wn * wn, -2.0 * zeta * wn}};
    sys.B = ColVec<2>{0.0, 1.0};
    sys.C = RowVec<2>{wn * wn, 0.0};
    sys.Ts = 0.0;

    const double peak = matlab::hinfnorm(sys);
    const double expected = 1.0 / (2.0 * zeta * std::sqrt(1.0 - zeta * zeta));
    CHECK(peak == doctest::Approx(expected).epsilon(2e-3));
}

TEST_CASE("H2 norm of a continuous system with feedthrough is infinite (nullopt)") {
    StateSpace<1, 1, 1> sys;
    sys.A = Matrix<1, 1>{{-1.0}};
    sys.B = Matrix<1, 1>{{1.0}};
    sys.C = Matrix<1, 1>{{1.0}};
    sys.D = Matrix<1, 1>{{0.5}}; // D ≠ 0 ⇒ ‖G‖₂ = ∞
    sys.Ts = 0.0;
    CHECK_FALSE(matlab::norm(sys).has_value());
}
