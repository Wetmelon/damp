// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>

#include "damp/backend.hpp"
#include "damp/design/lyapunov.hpp"
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

TEST_SUITE("Lyapunov on the Schur form") {
    template<size_t N>
    double maxabs(const Matrix<N, N>& M) {
        double m = 0.0;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                m = damp::max(m, std::abs(M(i, j)));
            }
        }
        return m;
    }

    // A random 6×6 with complex eigenvalues, shifted to be Hurwitz and scaled to be Schur
    // stable, a nonsymmetric Q and a two-column B.
    struct Case {
        Matrix<6, 6> Ac{};
        Matrix<6, 6> Ad{};
        Matrix<6, 6> Q{};
        Matrix<6, 2> B{};
    };

    Case make_case(unsigned& seed) {
        auto rnd = [&seed]() {
            seed = (seed * 1103515245U) + 12345U;
            return (static_cast<double>((seed >> 8) & 0xFFFFU) / 32768.0) - 1.0;
        };
        Case         c;
        Matrix<6, 6> A{};
        for (size_t i = 0; i < 6; ++i) {
            for (size_t j = 0; j < 6; ++j) {
                A(i, j) = rnd();
                c.Q(i, j) = rnd();
            }
            c.B(i, 0) = rnd();
            c.B(i, 1) = rnd();
        }
        const auto e = mat::compute_eigenvalues(A);
        double     shift = -1e300;
        double     rho = 0.0;
        for (size_t i = 0; i < 6; ++i) {
            shift = damp::max(shift, e.values[i].real());
            rho = damp::max(rho, e.values[i].abs());
        }
        c.Ac = A;
        for (size_t i = 0; i < 6; ++i) {
            c.Ac(i, i) -= shift + 0.1;
        }
        c.Ad = A * (0.95 / rho);
        return c;
    }

    TEST_CASE("Bartels–Stewart residuals are at rounding level") {
        unsigned seed = 3U;
        for (int trial = 0; trial < 100; ++trial) {
            const auto c = make_case(seed);
            const auto Xc = lyap(c.Ac, c.Q);
            const auto Xd = dlyap(c.Ad, c.Q);
            REQUIRE(Xc.has_value());
            REQUIRE(Xd.has_value());
            const double sc = (maxabs(c.Ac) * maxabs(*Xc)) + maxabs(c.Q);
            CHECK(cont_residual(c.Ac, *Xc, c.Q) <= 1e-13 * sc);
            const double sd = maxabs(*Xd) + maxabs(c.Q);
            CHECK(disc_residual(c.Ad, *Xd, c.Q) <= 1e-13 * sd);
        }
    }

    // The Kronecker form this replaced held an NX²×NX² system: 1.3 MB of stack at NX = 20.
    TEST_CASE("a 20-state Lyapunov equation solves in O(NX²) storage") {
        Matrix<20, 20> A{};
        Matrix<20, 20> Q{};
        for (size_t i = 0; i < 20; ++i) {
            A(i, i) = -1.0 - (0.1 * static_cast<double>(i));
            if (i + 1 < 20) {
                A(i, i + 1) = 1.0;
                A(i + 1, i) = -0.5;
            }
            Q(i, i) = 1.0;
        }
        const auto X = lyap(A, Q);
        REQUIRE(X.has_value());
        CHECK(cont_residual(A, *X, Q) <= 1e-13 * ((maxabs(A) * maxabs(*X)) + 1.0));
    }

    TEST_CASE("Hammarling's factor reproduces the Bartels–Stewart solution") {
        unsigned seed = 11U;
        for (int trial = 0; trial < 100; ++trial) {
            const auto         c = make_case(seed);
            const Matrix<6, 6> BBt = c.B * c.B.transpose();
            const auto         Rc = lyapunov_cholesky(c.Ac, c.B);
            const auto         Rd = discrete_lyapunov_cholesky(c.Ad, c.B);
            const auto         Xc = lyap(c.Ac, BBt);
            const auto         Xd = dlyap(c.Ad, BBt);
            REQUIRE(Rc.has_value());
            REQUIRE(Rd.has_value());
            REQUIRE(Xc.has_value());
            REQUIRE(Xd.has_value());
            CHECK(maxabs(Rc->transpose() * *Rc - *Xc) <= 1e-13 * maxabs(*Xc));
            CHECK(maxabs(Rd->transpose() * *Rd - *Xd) <= 1e-13 * maxabs(*Xd));
            for (size_t i = 0; i < 6; ++i) {
                CHECK((*Rc)(i, i) >= 0.0);
                CHECK((*Rd)(i, i) >= 0.0);
                for (size_t j = 0; j < i; ++j) {
                    CHECK((*Rc)(i, j) == 0.0);
                    CHECK((*Rd)(i, j) == 0.0);
                }
            }
        }
    }

    TEST_CASE("Hammarling takes a B wider than A and refuses an unstable A") {
        Matrix<3, 3> A{{-0.1, 1.0, 0.0}, {0.0, -1.0, 0.0}, {0.0, 0.0, -2.0}};
        Matrix<3, 5> B{};
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 5; ++j) {
                B(i, j) = static_cast<double>((i + 1) * (j + 2) % 7) - 3.0;
            }
        }
        const auto R = lyapunov_cholesky(A, B);
        const auto X = lyap(A, B * B.transpose());
        REQUIRE(R.has_value());
        REQUIRE(X.has_value());
        CHECK(maxabs(R->transpose() * *R - *X) <= 1e-13 * maxabs(*X));
        A(0, 0) = 0.1;
        CHECK_FALSE(lyapunov_cholesky(A, B).has_value());
        CHECK_FALSE(discrete_lyapunov_cholesky(A * 10.0, B).has_value());
        CHECK(matlab::lyapchol(Matrix<3, 3>{{-1.0, 0.0, 0.0}, {0.0, -2.0, 0.0}, {0.0, 0.0, -3.0}}, B).has_value());
    }

    // A = −diag(1, …, N), B = 1: the Gramian Wᵢⱼ = 1/(i + j) is a Cauchy matrix with Cholesky
    // pivots dₖ = 1/(2k)·Π_{j<k} ((k − j)/(k + j))², so R(k,k) = √dₖ is known exactly. At
    // N = 12 the factor of the computed Gramian has lost all but two digits of R's diagonal;
    // at N = 14 the computed Gramian no longer factors. Hammarling never forms it.
    template<size_t N>
    double hammarling_diagonal_error() {
        Matrix<N, N> A{};
        Matrix<N, 1> B{};
        for (size_t i = 0; i < N; ++i) {
            A(i, i) = -static_cast<double>(i + 1);
            B(i, 0) = 1.0;
        }
        const auto R = lyapunov_cholesky(A, B);
        REQUIRE(R.has_value());
        double worst = 0.0;
        for (size_t k = 0; k < N; ++k) {
            const double ak = static_cast<double>(k + 1);
            double       d = 1.0 / (2.0 * ak);
            for (size_t j = 0; j < k; ++j) {
                const double aj = static_cast<double>(j + 1);
                const double r = (ak - aj) / (ak + aj);
                d *= r * r;
            }
            worst = damp::max(worst, std::abs((*R)(k, k) - std::sqrt(d)) / std::sqrt(d));
        }
        return worst;
    }

    TEST_CASE("Hammarling resolves a Gramian that cholesky(lyap) cannot") {
        CHECK(hammarling_diagonal_error<12>() <= 1e-8);
        CHECK(hammarling_diagonal_error<14>() <= 1e-6);
        Matrix<14, 14> A{};
        Matrix<14, 1>  B{};
        for (size_t i = 0; i < 14; ++i) {
            A(i, i) = -static_cast<double>(i + 1);
            B(i, 0) = 1.0;
        }
        const auto X = lyap(A, B * B.transpose());
        REQUIRE(X.has_value());
        CHECK_FALSE(mat::cholesky(*X).has_value());
    }
}
