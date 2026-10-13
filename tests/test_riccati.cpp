// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/design/riccati.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

// The Riccati solvers (care/dare, in namespace damp) are
// foundational: LQR, LQI, LQG and LQGI all stand on top of dare(), yet they had
// no direct test of their own. These tests validate by residual (plug the
// returned P back into the defining equation and check the residual is ~0) so
// they don't depend on a hand-derived closed-form solution, plus a cross-check
// that design::discrete_lqr's Riccati solution S matches dare() on the same
// system.
//
// @see "Optimal Control" (Anderson & Moore, 1990), §3.3 (CARE) and Ch. 4 (DARE).

TEST_SUITE("Riccati Solvers") {
    // ---- Continuous-time algebraic Riccati equation: AᵀP + PA − PBR⁻¹BᵀP + Q = 0 ----

    TEST_CASE("CARE scalar has known solution P = 1") {
        // a = 0, b = 1, q = 1, r = 1  ⇒  −P² + 1 = 0  ⇒  P = 1.
        constexpr Matrix<1, 1> A{{0.0}};
        constexpr Matrix<1, 1> B{{1.0}};
        constexpr Matrix<1, 1> Q{{1.0}};
        constexpr Matrix<1, 1> R{{1.0}};

        const auto P = care(A, B, Q, R);
        REQUIRE(P.has_value());
        CHECK(P.value()(0, 0) == doctest::Approx(1.0));
    }

    TEST_CASE("CARE residual is ~0 for a double integrator") {
        // Double integrator: position/velocity, acceleration input.
        constexpr Matrix<2, 2> A{{0.0, 1.0}, {0.0, 0.0}};
        constexpr Matrix<2, 1> B{{0.0}, {1.0}};
        constexpr auto         Q = Matrix<2, 2>::identity();
        constexpr auto         R = Matrix<1, 1>::identity(); // R = I ⇒ R⁻¹ = I, no inverse needed

        const auto P_opt = care(A, B, Q, R);
        REQUIRE(P_opt.has_value());
        const Matrix<2, 2>& P = P_opt.value();

        // Residual: AᵀP + PA − PB(R⁻¹)BᵀP + Q  (R⁻¹ = I)
        const Matrix<2, 2> residual = A.t() * P + P * A - P * B * B.t() * P + Q;
        CHECK(residual.norm() == doctest::Approx(0.0).epsilon(1e-6));

        // P of a CARE solution is symmetric positive definite.
        CHECK(P(0, 1) == doctest::Approx(P(1, 0)));
        CHECK(P(0, 0) > 0.0);
    }

    TEST_CASE("CARE is usable at compile time (constexpr)") {
        constexpr Matrix<1, 1> A{{0.0}};
        constexpr Matrix<1, 1> B{{1.0}};
        constexpr Matrix<1, 1> Q{{1.0}};
        constexpr Matrix<1, 1> R{{1.0}};
        constexpr auto         P = care(A, B, Q, R);
        static_assert(P.has_value(), "scalar CARE must converge at compile time");
    }

    TEST_CASE("CARE scalar cross-term has the stabilizing root") {
        // (X + N)² = Q with N = 1/2, Q = R = 1 ⇒ X = 1/2 (closed loop pole −1).
        // X = −3/2 also has a zero residual and an unstable closed loop.
        constexpr Matrix<1, 1> A{{0.0}};
        constexpr Matrix<1, 1> B{{1.0}};
        constexpr Matrix<1, 1> Q{{1.0}};
        constexpr Matrix<1, 1> R{{1.0}};
        constexpr Matrix<1, 1> N{{0.5}};
        const auto             P = care(A, B, Q, R, N);
        REQUIRE(P.has_value());
        CHECK(P.value()(0, 0) == doctest::Approx(0.5));
    }

    TEST_CASE("CARE tiny input weight matches the analytic double integrator") {
        // q1 = q2 = 1, r = 1e−8. Stabilizing root:
        //   p12 = √(q1 r),  p22 = √(r (q2 + 2 p12)),  p11 = p12 p22 / r.
        constexpr Matrix<2, 2> A{{0.0, 1.0}, {0.0, 0.0}};
        constexpr Matrix<2, 1> B{{0.0}, {1.0}};
        constexpr auto         Q = Matrix<2, 2>::identity();
        constexpr double       r = 1e-8;
        constexpr Matrix<1, 1> R{{r}};
        const auto             P = care(A, B, Q, R);
        REQUIRE(P.has_value());
        const double p12 = damp::sqrt(r);
        const double p22 = damp::sqrt(r * (1.0 + 2.0 * p12));
        const double p11 = p12 * p22 / r;
        CHECK(P.value()(0, 0) == doctest::Approx(p11).epsilon(1e-8));
        CHECK(P.value()(0, 1) == doctest::Approx(p12).epsilon(1e-8));
        CHECK(P.value()(1, 0) == doctest::Approx(p12).epsilon(1e-8));
        CHECK(P.value()(1, 1) == doctest::Approx(p22).epsilon(1e-8));
    }

    // CAREX cases from Benner, Laub, Mehrmann, SPC 95_23. Reference solutions are
    // scipy.linalg.solve_continuous_are (1.17, balanced). A returned X must also
    // have a small relative residual; CAREX #12 is refused because even the
    // balanced residual sits near 3×10⁻⁴, above the certificate.

    TEST_CASE("CARE CAREX 7 matches the balanced reference") {
        constexpr Matrix<2, 2> A{{1.0, 0.0}, {0.0, -2.0}};
        constexpr Matrix<2, 1> B{{1e-6}, {0.0}};
        constexpr Matrix<2, 2> Q{{1.0, 1.0}, {1.0, 1.0}};
        constexpr Matrix<1, 1> R{{1.0}};
        const auto             P = care(A, B, Q, R);
        REQUIRE(P.has_value());
        CHECK(P.value()(0, 0) == doctest::Approx(2e12).epsilon(1e-8));
        CHECK(P.value()(0, 1) == doctest::Approx(1.0 / 3.0).epsilon(1e-6));
        CHECK(P.value()(1, 1) == doctest::Approx(0.25).epsilon(1e-8));
    }

    TEST_CASE("CARE CAREX 9 matches the balanced reference") {
        constexpr Matrix<2, 2> A{{0.0, 1e6}, {0.0, 0.0}};
        constexpr Matrix<2, 1> B{{0.0}, {1.0}};
        constexpr auto         Q = Matrix<2, 2>::identity();
        constexpr Matrix<1, 1> R{{1.0}};
        const auto             P = care(A, B, Q, R);
        REQUIRE(P.has_value());
        const double s = damp::sqrt(2.0);
        CHECK(P.value()(0, 0) == doctest::Approx(s * 1e-3).epsilon(1e-6));
        CHECK(P.value()(0, 1) == doctest::Approx(1.0).epsilon(1e-8));
        CHECK(P.value()(1, 1) == doctest::Approx(s * 1e3).epsilon(1e-6));
    }

    TEST_CASE("CARE CAREX 13 matches the balanced reference") {
        constexpr Matrix<4, 4> A{{0.0, 0.4, 0.0, 0.0}, {0.0, 0.0, 0.345, 0.0}, {0.0, -5.24e5, -4.65e5, 2.62e5}, {0.0, 0.0, 0.0, -1e6}};
        constexpr Matrix<4, 1> B{{0.0}, {0.0}, {0.0}, {1e6}};
        Matrix<4, 4>           Q = Matrix<4, 4>::zeros();
        Q(0, 0) = 1.0;
        Q(2, 2) = 1.0;
        constexpr Matrix<1, 1> R{{1.0}};
        const auto             P = care(A, B, Q, R);
        REQUIRE(P.has_value());
        constexpr Matrix<4, 4> Xref{{7.38402467, 5.90476206, 3.99308234e-6, 1.00000000e-6}, {5.90476206, 7.15160630, 3.79969886e-6, 8.61234718e-7}, {3.99308234e-6, 3.79969886e-6, 1.04029358e-6, 1.80359619e-7}, {1.00000000e-6, 8.61234718e-7, 1.80359619e-7, 4.61875742e-8}};
        CHECK((P.value() - Xref).norm() / Xref.norm() < 1e-7);
    }

    TEST_CASE("CARE CAREX 14 stays near the identity with a tiny residual") {
        // Eigenvalues sit on the imaginary axis. The residual is flat, so X
        // agrees with SciPy to about 3×10⁻⁴; a false subspace used to return
        // a residual of 0.46.
        constexpr Matrix<4, 4> A{{-1e-6, 1.0, 0.0, 0.0}, {-1.0, -1e-6, 0.0, 0.0}, {0.0, 0.0, 1e-6, 1.0}, {0.0, 0.0, -1.0, 1e-6}};
        constexpr Matrix<4, 1> B{{1.0}, {1.0}, {1.0}, {1.0}};
        constexpr Matrix<4, 4> Q{{1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}};
        constexpr Matrix<1, 1> R{{1.0}};
        const auto             P = care(A, B, Q, R);
        REQUIRE(P.has_value());
        const Matrix<4, 4> G = B * B.transpose();
        const auto         AtX = A.transpose() * P.value();
        const auto         XA = P.value() * A;
        const auto         XGX = P.value() * G * P.value();
        const auto         Res = AtX + XA - XGX + Q;
        const double       rel = Res.norm() / (AtX.norm() + XA.norm() + XGX.norm() + Q.norm());
        CHECK(rel < 1e-8);
        for (size_t i = 0; i < 4; ++i) {
            CHECK(P.value()(i, i) == doctest::Approx(1.0).epsilon(1e-3));
            for (size_t j = 0; j < 4; ++j) {
                if (i != j) {
                    CHECK(damp::abs(P.value()(i, j)) < 1e-2);
                }
            }
        }
    }

    TEST_CASE("float CARE succeeds on a double integrator") {
        constexpr Matrix<2, 2, float> A{{0.0f, 1.0f}, {0.0f, 0.0f}};
        constexpr Matrix<2, 1, float> B{{0.0f}, {1.0f}};
        constexpr auto                Q = Matrix<2, 2, float>::identity();
        constexpr auto                R = Matrix<1, 1, float>::identity();
        const auto                    P = care(A, B, Q, R);
        REQUIRE(P.has_value());
        CHECK(P.value()(0, 0) == doctest::Approx(damp::sqrt(3.0f)).epsilon(1e-4));
        CHECK(P.value()(0, 1) == doctest::Approx(1.0f).epsilon(1e-4));
    }

    // A = V·diag(1, 2, 3)·10⁶·V and Q = V·diag(10⁻⁶, 1, 10⁶)·V share the eigenvectors of the
    // reflection V = M/3, so with B = I and R = 10⁶·I the equation decouples into three scalar
    // ones: xᵢ = r(λᵢ + √(λᵢ² + qᵢ/r)). Until the Schur form deflated by Ahues–Tisseur this
    // example came back with a residual the certificate refused; SciPy's solve_continuous_are
    // is off by 5e-4 relative here.
    TEST_CASE("CARE CAREX 12 matches its closed form") {
        constexpr Matrix<3, 3> A{{7000000.0 / 3.0, 2000000.0 / 3.0, 0.0}, {2000000.0 / 3.0, 6000000.0 / 3.0, -2000000.0 / 3.0}, {0.0, -2000000.0 / 3.0, 5000000.0 / 3.0}};
        constexpr auto         B = Matrix<3, 3>::identity();
        constexpr Matrix<3, 3> M{{1.0, -2.0, -2.0}, {-2.0, 1.0, -2.0}, {-2.0, -2.0, 1.0}};
        Matrix<3, 3>           D = Matrix<3, 3>::zeros();
        D(0, 0) = 1e-6;
        D(1, 1) = 1.0;
        D(2, 2) = 1e6;
        const Matrix<3, 3> Q = (M * D * M) * (1.0 / 9.0);
        constexpr auto     R = Matrix<3, 3>::identity() * 1e6;
        const auto         X = care(A, B, Q, R);
        REQUIRE(X.has_value());
        const Matrix<3, 3> V = M * (1.0 / 3.0);
        Matrix<3, 3>       Xd = Matrix<3, 3>::zeros();
        const double       lam[3] = {1e6, 2e6, 3e6};
        const double       q[3] = {1e-6, 1.0, 1e6};
        for (size_t i = 0; i < 3; ++i) {
            Xd(i, i) = 1e6 * (lam[i] + std::sqrt((lam[i] * lam[i]) + (q[i] / 1e6)));
        }
        const Matrix<3, 3> exact = V * Xd * V;
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                CHECK(std::fabs((*X)(i, j) - exact(i, j)) <= 1e-5 * 6e12);
            }
        }
    }

    // ---- Discrete-time algebraic Riccati equation ----
    // AᵀPA − P − AᵀPB(R + BᵀPB)⁻¹BᵀPA + Q = 0

    TEST_CASE("DARE residual is ~0 for a scalar system") {
        constexpr Matrix<1, 1> A{{1.0}};
        constexpr Matrix<1, 1> B{{1.0}};
        constexpr Matrix<1, 1> Q{{1.0}};
        constexpr Matrix<1, 1> R{{1.0}};

        const auto P_opt = dare(A, B, Q, R);
        REQUIRE(P_opt.has_value());
        const Matrix<1, 1>& P = P_opt.value();

        const Matrix<1, 1> S = R + B.t() * P * B; // scalar, positive
        const Matrix<1, 1> Sinv{{1.0 / S(0, 0)}};
        const Matrix<1, 1> residual = A.t() * P * A - P - A.t() * P * B * Sinv * B.t() * P * A + Q;
        CHECK(residual.norm() == doctest::Approx(0.0).epsilon(1e-6));
        CHECK(P(0, 0) > 0.0);
    }

    TEST_CASE("DARE residual is ~0 for a discretized double integrator") {
        // x[k+1] = A x[k] + B u[k], Ts = 0.1 s.
        constexpr Matrix<2, 2> A{{1.0, 0.1}, {0.0, 1.0}};
        constexpr Matrix<2, 1> B{{0.005}, {0.1}};
        constexpr auto         Q = Matrix<2, 2>::identity();
        constexpr auto         R = Matrix<1, 1>::identity();

        const auto P_opt = dare(A, B, Q, R);
        REQUIRE(P_opt.has_value());
        const Matrix<2, 2>& P = P_opt.value();

        const Matrix<1, 1> S = R + B.t() * P * B;
        const Matrix<1, 1> Sinv{{1.0 / S(0, 0)}};
        const Matrix<2, 2> residual = A.t() * P * A - P - A.t() * P * B * Sinv * B.t() * P * A + Q;
        CHECK(residual.norm() == doctest::Approx(0.0).epsilon(1e-6));

        // Symmetric, positive definite.
        CHECK(P(0, 1) == doctest::Approx(P(1, 0)));
        CHECK(P(0, 0) > 0.0);
    }

    TEST_CASE("DARE solution matches discrete_lqr's Riccati solution S") {
        // discrete_lqr solves the same DARE internally; its S field must agree
        // with calling design::dare directly on the same (A, B, Q, R).
        constexpr Matrix<2, 2> A{{1.0, 0.1}, {0.0, 1.0}};
        constexpr Matrix<2, 1> B{{0.005}, {0.1}};
        constexpr auto         Q = Matrix<2, 2>::identity();
        constexpr auto         R = Matrix<1, 1>::identity();

        const auto P_opt = dare(A, B, Q, R);
        REQUIRE(P_opt.has_value());

        const auto lqr = design::discrete_lqr(A, B, Q, R);
        REQUIRE(lqr.success);

        const Matrix<2, 2> diff = P_opt.value() - lqr.S;
        CHECK(diff.norm() == doctest::Approx(0.0).epsilon(1e-6));
    }

    TEST_CASE("DARE is usable at compile time (constexpr)") {
        constexpr Matrix<1, 1> A{{1.0}};
        constexpr Matrix<1, 1> B{{1.0}};
        constexpr Matrix<1, 1> Q{{1.0}};
        constexpr Matrix<1, 1> R{{1.0}};
        constexpr auto         P = damp::dare(A, B, Q, R);
        static_assert(P.has_value(), "scalar DARE must converge at compile time");
    }

    TEST_CASE("float DARE succeeds on a simple discrete plant") {
        // Hard-coded 1e-12 tol fails float DARE; default_tol<float>() = 1e-6.
        constexpr Matrix<1, 1, float> A{{1.0f}};
        constexpr Matrix<1, 1, float> B{{1.0f}};
        constexpr Matrix<1, 1, float> Q{{1.0f}};
        constexpr Matrix<1, 1, float> R{{1.0f}};

        const auto P_opt = dare(A, B, Q, R);
        REQUIRE(P_opt.has_value());
        CHECK(P_opt.value()(0, 0) > 0.0f);
        CHECK(std::isfinite(P_opt.value()(0, 0)));
    }

    TEST_CASE("float DARE succeeds on discretized double integrator") {
        constexpr Matrix<2, 2, float> A{{1.0f, 0.1f}, {0.0f, 1.0f}};
        constexpr Matrix<2, 1, float> B{{0.005f}, {0.1f}};
        constexpr auto                Q = Matrix<2, 2, float>::identity();
        constexpr auto                R = Matrix<1, 1, float>::identity();

        const auto P_opt = dare(A, B, Q, R);
        REQUIRE(P_opt.has_value());
        CHECK(P_opt.value()(0, 0) > 0.0f);
        CHECK(std::isfinite(P_opt.value()(0, 0)));
        CHECK(std::isfinite(P_opt.value()(1, 1)));
    }
}

TEST_SUITE("DARE: Cross-Term N Support") {
    TEST_CASE("dare with zero cross-term N matches no-N variant") {
        // Verify that dare(A, B, Q, R, N={0}) produces same result as dare(A, B, Q, R)
        Matrix<2, 2> A{{0.95, 0.1}, {0.0, 0.9}};
        Matrix<2, 1> B{{0.1}, {0.1}};
        Matrix<2, 2> Q = Matrix<2, 2>::identity();
        Matrix<1, 1> R = {{1.0}};
        Matrix<2, 1> N{};

        auto P_no_N = dare(A, B, Q, R);
        auto P_with_zero_N = dare(A, B, Q, R, N);

        CHECK(P_no_N.has_value());
        CHECK(P_with_zero_N.has_value());

        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                CHECK(doctest::Approx(P_no_N.value()(i, j)).epsilon(1e-10) == P_with_zero_N.value()(i, j));
            }
        }
    }

    TEST_CASE("dare with non-zero cross-term N produces valid result") {
        Matrix<1, 1> A = {{1.0}};
        Matrix<1, 1> B = {{1.0}};
        Matrix<1, 1> Q = {{1.0}};
        Matrix<1, 1> R = {{1.0}};
        Matrix<1, 1> N = {{0.1}};

        auto P_opt = dare(A, B, Q, R, N);
        auto P0 = dare(A, B, Q, R);
        REQUIRE(P_opt.has_value());
        REQUIRE(P0.has_value());

        double P = P_opt.value()(0, 0);
        CHECK(P > 0.0);
        // Non-zero N must change the solution
        CHECK(P != doctest::Approx(P0.value()(0, 0)));
    }

    TEST_CASE("dare rejects a cross term that makes the joint cost indefinite") {
        // Q − N R⁻¹ Nᵀ = 1 − 4 < 0 although Q and R are each positive definite.
        const Matrix<1, 1> A{{0.9}}, B{{1.0}}, Q{{1.0}}, R{{1.0}}, N{{2.0}};
        CHECK_FALSE(dare(A, B, Q, R, N).has_value());
    }
}

TEST_SUITE("DARE: stabilizing certificate") {
    TEST_CASE("Q = 0 with unstable A fails closed instead of returning X = 0") {
        // SDA and RDE both sit at the non-stabilizing fixed point X = 0 (K = 0,
        // closed-loop pole 2). The stabilizing solution is X = 3.
        const Matrix<1, 1> A{{2.0}}, B{{1.0}}, Q{{0.0}}, R{{1.0}};
        CHECK_FALSE(dare(A, B, Q, R).has_value());
        CHECK_FALSE(dare(A, B, Q, R, Matrix<1, 1>{}, DareMethod::RDE).has_value());
        CHECK_FALSE(design::discrete_lqr(A, B, Q, R).success);
    }

    TEST_CASE("a uniformly tiny Q is judged like a unit one") {
        const Matrix<1, 1> A{{0.9}}, B{{1.0}}, R{{1.0}};
        CHECK(dare(A, B, Matrix<1, 1>{{1e-14}}, R).has_value());
        CHECK_FALSE(dare(A, B, Matrix<1, 1>{{-1e-14}}, R).has_value());
    }
}

TEST_SUITE("DARE: R=0 (Positive Semidefinite R)") {
    TEST_CASE("dare with R=0, scalar system (deadbeat)") {
        // A=1, B=1, Q=1, R=0 → X = A'XA + Q - A'XB(B'XB)⁻¹B'XA = X - X + 1 = 1
        Matrix<1, 1> A{{1.0}};
        Matrix<1, 1> B{{1.0}};
        Matrix<1, 1> Q{{1.0}};
        Matrix<1, 1> R{{0.0}};

        auto P = dare(A, B, Q, R);
        REQUIRE(P.has_value());
        CHECK(P.value()(0, 0) == doctest::Approx(1.0).epsilon(1e-10));
    }

    TEST_CASE("dare with R=0, 2x2 system") {
        // scipy: solve_discrete_are(A.T, B.T, Q, zeros(1,1))
        // A = [[0.9, 0.1],[0, 0.8]], B = [[0.1],[0.2]], Q = I, R = 0
        Matrix<2, 2> A{{0.9, 0.1}, {0.0, 0.8}};
        Matrix<2, 1> B{{0.1}, {0.2}};
        Matrix<2, 2> Q = Matrix<2, 2>::identity();
        Matrix<1, 1> R{{0.0}};

        auto P = dare(A, B, Q, R);
        REQUIRE(P.has_value());
        // Verify DARE residual: A'PA - P - A'PB(B'PB)⁻¹B'PA + Q ≈ 0
        const auto& Pv = P.value();
        auto        BtP = B.transpose() * Pv;
        auto        BtPB = BtP * B;
        auto        BtPA = BtP * A;
        auto        residual = A.transpose() * Pv * A - Pv
                      - A.transpose() * Pv * B * (BtPA * (1.0 / BtPB(0, 0))) + Q;
        CHECK(residual.norm() < 1e-8);
    }

    TEST_CASE("dare with diagonal PSD R (partially singular)") {
        Matrix<2, 2> A{{0.9, 0.1}, {0.0, 0.8}};
        Matrix<2, 2> B{{1.0, 0.0}, {0.0, 1.0}};
        Matrix<2, 2> Q = Matrix<2, 2>::identity();
        Matrix<2, 2> R{{1.0, 0.0}, {0.0, 0.0}};

        auto P = dare(A, B, Q, R);
        REQUIRE(P.has_value());
        // Verify symmetry and positive semidefiniteness
        auto Pv = P.value();
        CHECK(Pv(0, 1) == doctest::Approx(Pv(1, 0)).epsilon(1e-10));
        CHECK(Pv(0, 0) > 0.0);
        CHECK(Pv(1, 1) > 0.0);
    }

    TEST_CASE("DareMethod::SDA rejects R=0") {
        Matrix<2, 2> A{{0.9, 0.1}, {0.0, 0.8}};
        Matrix<2, 1> B{{0.1}, {0.2}};
        Matrix<2, 2> Q = Matrix<2, 2>::identity();
        Matrix<1, 1> R{{0.0}};

        auto P = dare(A, B, Q, R, Matrix<2, 1>{}, DareMethod::SDA);
        CHECK_FALSE(P.has_value());
    }

    TEST_CASE("DareMethod::RDE succeeds for R=0") {
        Matrix<2, 2> A{{0.9, 0.1}, {0.0, 0.8}};
        Matrix<2, 1> B{{0.1}, {0.2}};
        Matrix<2, 2> Q = Matrix<2, 2>::identity();
        Matrix<1, 1> R{{0.0}};

        auto P = dare(A, B, Q, R, Matrix<2, 1>{}, DareMethod::RDE);
        REQUIRE(P.has_value());
    }

    TEST_CASE("SDA and RDE agree for R>0") {
        Matrix<2, 2> A{{0.9, 0.1}, {0.0, 0.8}};
        Matrix<2, 1> B{{0.1}, {0.2}};
        Matrix<2, 2> Q = Matrix<2, 2>::identity();
        Matrix<1, 1> R{{1.0}};

        auto P_sda = dare(A, B, Q, R, Matrix<2, 1>{}, DareMethod::SDA);
        auto P_rde = dare(A, B, Q, R, Matrix<2, 1>{}, DareMethod::RDE);

        REQUIRE(P_sda.has_value());
        REQUIRE(P_rde.has_value());

        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                CHECK(doctest::Approx(P_sda.value()(i, j)).epsilon(1e-8) == P_rde.value()(i, j));
            }
        }
    }

    TEST_CASE("dare with R>0 regression (SDA path unchanged)") {
        constexpr Matrix<2, 2> A{{0.9, 0.1}, {0.0, 0.8}};
        constexpr Matrix<2, 1> B{{0.1}, {0.2}};
        constexpr Matrix<2, 2> Q{{1.0, 0.0}, {0.0, 1.0}};
        constexpr Matrix<1, 1> R{{1.0}};

        constexpr auto P = dare(A, B, Q, R);
        static_assert(P.has_value());
        CHECK(doctest::Approx(P.value()(0, 0)).epsilon(1e-6) == 4.201439333611463);
        CHECK(doctest::Approx(P.value()(0, 1)).epsilon(1e-6) == 0.5954889098912259);
        CHECK(doctest::Approx(P.value()(1, 1)).epsilon(1e-6) == 2.543807455993601);
    }
}

TEST_SUITE("is_stabilizable") {
    TEST_CASE("one large-gain input does not hide the others (column scaling)") {
        // Double integrator in each of two axes. Input 0 drives axis 0 only, with a gain 10⁸
        // times input 1's, which drives axis 1 only. Scaling a column of B cannot change
        // stabilizability, so this pair is stabilizable whatever the gains.
        Matrix<4, 4> A = Matrix<4, 4>::identity();
        A(0, 1) = 1.0;
        A(2, 3) = 1.0;
        Matrix<4, 2> B{};
        B(1, 0) = 1.0e8;
        B(3, 1) = 1.0;
        CHECK(is_stabilizable(A, B));
        Matrix<4, 2> B_axis0_only{};
        B_axis0_only(1, 0) = 1.0e8;
        CHECK_FALSE(is_stabilizable(A, B_axis0_only)); // axis 1 has no input at all
    }

    TEST_CASE("a badly scaled change of coordinates does not change the verdict (balancing)") {
        // T = diag(1, 1e6, 1, 1e-4) mixes scales across A's rows the way a delay register copying an
        // accelerometer row does. Stabilizability is invariant under it.
        Matrix<4, 4> A = Matrix<4, 4>::identity();
        A(0, 1) = 1.0;
        A(2, 3) = 1.0;
        Matrix<4, 2> B{};
        B(1, 0) = 1.0;
        B(3, 1) = 1.0;
        const double s[4] = {1.0, 1.0e6, 1.0, 1.0e-4};
        Matrix<4, 4> As{};
        Matrix<4, 2> Bs{};
        for (size_t r = 0; r < 4; ++r) {
            for (size_t c = 0; c < 4; ++c) {
                As(r, c) = s[r] * A(r, c) / s[c];
            }
            for (size_t c = 0; c < 2; ++c) {
                Bs(r, c) = s[r] * B(r, c);
            }
        }
        CHECK(is_stabilizable(As, Bs));
        Matrix<4, 2> Bs_axis0{};
        Bs_axis0(1, 0) = Bs(1, 0);
        CHECK_FALSE(is_stabilizable(As, Bs_axis0)); // axis 1 still has no input
    }

    TEST_CASE_TEMPLATE("an accelerometer register next to an integrator disturbance: the dual is detectable", T, float, double) {
        // Mass-spring at 1 kHz with k/m = 1e5. A register copies the accelerometer row −k/m·x₁ + d one
        // sample late (entries of 1e5 next to rows of 1e-3); d is an integrator disturbance. Outputs: the
        // position and the register. Row scaling alone called the dual not detectable; the pair itself is
        // not stabilizable (d has no input) in either precision.
        const T         dt = static_cast<T>(1e-3);
        const T         km = static_cast<T>(1e5);
        Matrix<4, 4, T> A{};
        A(0, 0) = T{1};
        A(0, 1) = dt;
        A(1, 0) = -dt * km;
        A(1, 1) = T{1};
        A(1, 3) = dt;
        A(2, 0) = -km;
        A(2, 3) = T{1};
        A(3, 3) = T{1};
        Matrix<4, 1, T> B{};
        B(1, 0) = dt;
        Matrix<2, 4, T> C{};
        C(0, 0) = T{1};
        C(1, 2) = T{1};
        CHECK(is_stabilizable(Matrix<4, 4, T>(A.t()), Matrix<4, 2, T>(C.t())));
        CHECK_FALSE(is_stabilizable(A, B));
        A(3, 3) = T{0}; // without the disturbance the plant alone is stabilizable
        A(1, 3) = T{0};
        A(2, 3) = T{0};
        CHECK(is_stabilizable(A, B));
    }

    TEST_CASE_TEMPLATE("an undamped mode with no input is never skipped as stable", T, float, double) {
        // |λ| = 1 exactly; the computed magnitude lands on either side of 1 by roundoff. Within √eps of the
        // unit circle counts as unstable, so every angle gives the same verdict.
        for (int k = 1; k <= 20; ++k) {
            const T         th = static_cast<T>(0.05 * k);
            Matrix<3, 3, T> A{};
            A(0, 0) = damp::cos(th);
            A(0, 1) = -damp::sin(th);
            A(1, 0) = damp::sin(th);
            A(1, 1) = damp::cos(th);
            A(2, 2) = static_cast<T>(0.5);
            Matrix<3, 1, T> B{};
            B(2, 0) = T{1};
            CHECK_FALSE(is_stabilizable(A, B));
        }
    }
}

TEST_SUITE("DARE stabilizability certificate") {
    TEST_CASE("an uncontrollable mode on the unit circle is refused after the solve") {
        // Axis 1 is a double integrator with no input; Q does not weight it, so a
        // Riccati fixed point exists with that mode left at z = 1 in A − BK. The
        // closed-loop PBH test must still refuse the pair.
        Matrix<4, 4> A = Matrix<4, 4>::identity();
        A(0, 1) = 0.1;
        A(2, 3) = 0.1;
        Matrix<4, 1> B{};
        B(1, 0) = 0.1;
        Matrix<4, 4> Q{};
        Q(0, 0) = 1.0;
        Q(1, 1) = 1.0;
        CHECK_FALSE(dare(A, B, Q, Matrix<1, 1>{{1.0}}).has_value());
        CHECK_FALSE(dare(A, B, Matrix<4, 4>::identity(), Matrix<1, 1>{{1.0}}).has_value());

        // An uncontrollable but stable mode is fine.
        Matrix<4, 4> A_stable = A;
        A_stable(2, 2) = 0.5;
        A_stable(3, 3) = 0.5;
        CHECK(dare(A_stable, B, Matrix<4, 4>::identity(), Matrix<1, 1>{{1.0}}).has_value());
    }
}

namespace {
// Chain of five unit masses on springs (10 states), forward-Euler at 10 ms; forces on
// the end masses, their positions measured. Lightly damped, so every open-loop pole
// sits near the unit circle — the costly case for the design-time checks.
constexpr StateSpace<10, 2, 2> mass_chain() {
    StateSpace<10, 2, 2> s{};
    const double         Ts = 0.01;
    const double         k = 50.0;
    for (size_t i = 0; i < 10; ++i) {
        s.A(i, i) = 1.0;
    }
    for (size_t m = 0; m < 5; ++m) {
        s.A(2 * m, (2 * m) + 1) = Ts;
        double kk = 0.0;
        if (m > 0) {
            s.A((2 * m) + 1, 2 * (m - 1)) += Ts * k;
            kk += k;
        }
        if (m + 1 < 5) {
            s.A((2 * m) + 1, 2 * (m + 1)) += Ts * k;
            kk += k;
        }
        s.A((2 * m) + 1, 2 * m) -= Ts * kk;
        s.A((2 * m) + 1, (2 * m) + 1) -= Ts * 0.2;
    }
    s.B(1, 0) = Ts;
    s.B(9, 1) = Ts;
    s.C(0, 0) = 1.0;
    s.C(1, 8) = 1.0;
    s.Ts = Ts;
    return s;
}
} // namespace

TEST_CASE("10-state LQI and steady-state Kalman designs fit the default constexpr budget") {
    // Tripwire for -fconstexpr-ops-limit (2^25): both must stay designable at compile time.
    constexpr auto sys = mass_chain();
    constexpr auto lqi = design::discrete_lqi(sys, Matrix<12, 12>::identity(), Matrix<2, 2>::identity());
    static_assert(lqi.success);
    constexpr StateSpace<10, 2, 2, double, 10, 2> noisy{.A = sys.A, .B = sys.B, .C = sys.C, .G = Matrix<10, 10>::identity(), .H = Matrix<2, 2>::identity(), .Ts = sys.Ts};
    constexpr auto                                kf = design::kalman(noisy, Matrix<10, 10>::identity() * 1e-4, Matrix<2, 2>::identity() * 1e-2);
    static_assert(kf.success);
    CHECK(lqi.success);
    CHECK(kf.success);
}

TEST_SUITE("DARE scale invariance") {
    // Scaling Q and R together scales X and leaves K unchanged. A stop test with an absolute
    // floor (tol · max(1, ‖X‖)) declared a covariance-sized X (‖X‖ ≪ 1) converged after one
    // iteration and returned a gain five orders of magnitude off.
    TEST_CASE_TEMPLATE("a solution far below unit size converges to the same gain", T, float, double) {
        Matrix<4, 4, T> A = Matrix<4, 4, T>::identity();
        A(0, 1) = static_cast<T>(1e-3);
        A(2, 3) = static_cast<T>(1e-3);
        Matrix<4, 2, T> B{};
        B(1, 0) = static_cast<T>(1e-3);
        B(3, 1) = static_cast<T>(1e-3);
        const Matrix<4, 4, T> Q0 = Matrix<4, 4, T>::identity();
        const Matrix<2, 2, T> R0 = Matrix<2, 2, T>::identity() * static_cast<T>(0.1);
        const auto            X0 = dare(A, B, Q0, R0);
        REQUIRE(X0.has_value());
        const auto K0 = lqr_gain(A, B, *X0, R0);
        REQUIRE(K0.has_value());
        const T    scale = std::is_same_v<T, float> ? static_cast<T>(1e-8) : static_cast<T>(1e-16);
        const auto X = dare(A, B, Matrix<4, 4, T>(Q0 * scale), Matrix<2, 2, T>(R0 * scale));
        REQUIRE(X.has_value());
        const auto K = lqr_gain(A, B, *X, Matrix<2, 2, T>(R0 * scale));
        REQUIRE(K.has_value());
        for (size_t r = 0; r < 2; ++r) {
            for (size_t c = 0; c < 4; ++c) {
                CHECK((*K)(r, c) == doctest::Approx((*K0)(r, c)).epsilon(static_cast<double>(T{100} * damp::sqrt(eps<T>()))));
            }
        }
        CHECK(X->norm() / scale == doctest::Approx(X0->norm()).epsilon(1e-6));
    }
}

TEST_SUITE("real Schur block swap") {
    // A standardized quasi-triangular window: P×P block (eigenvalue re or re ± i·im), Q×Q
    // block at re + gap, coupling of size beta.
    template<size_t P, size_t Q>
    Matrix<P + Q, P + Q> window(double re, double gap, double im, double beta, double skew, unsigned& seed) {
        auto rnd = [&seed]() {
            seed = (seed * 1103515245U) + 12345U;
            return (static_cast<double>((seed >> 8) & 0xFFFFU) / 32768.0) - 1.0;
        };
        Matrix<P + Q, P + Q> W{};
        const auto           put = [&W](size_t o, double r, double i, double k) {
            W(o, o) = r;
            W(o + 1, o + 1) = r;
            W(o, o + 1) = i * k;
            W(o + 1, o) = -i / k;
        };
        if constexpr (P == 2) {
            put(0, re, im, skew);
        } else {
            W(0, 0) = re;
        }
        if constexpr (Q == 2) {
            put(P, re + gap, im, 1.0 / skew);
        } else {
            W(P, P) = re + gap;
        }
        for (size_t r = 0; r < P; ++r) {
            for (size_t c = P; c < P + Q; ++c) {
                W(r, c) = beta * rnd();
            }
        }
        return W;
    }

    // An accepted swap is an orthogonal similarity that puts the Q block first, stores the
    // vanished block as exact zeros and leaves standardized 2×2 blocks; a refused one leaves
    // the form and the Schur vectors as they were.
    template<size_t P, size_t Q>
    void check_swaps() {
        constexpr size_t S = P + Q;
        unsigned         seed = 99U;
        int              accepted = 0;
        for (int trial = 0; trial < 400; ++trial) {
            const double gap = std::pow(10.0, -static_cast<double>(trial % 13));
            const double beta = std::pow(10.0, static_cast<double>(trial % 7));
            const double skew = std::pow(10.0, static_cast<double>((trial % 5) - 2));
            const auto   T0 = window<P, Q>(0.3, gap, 0.8, beta, skew, seed);
            auto         Tm = T0;
            auto         Z = Matrix<S, S>::identity();
            if (!damp::detail::swap_schur_blocks<P, Q>(Tm, Z, 0)) {
                for (size_t r = 0; r < S; ++r) {
                    for (size_t c = 0; c < S; ++c) {
                        CHECK(Tm(r, c) == T0(r, c));
                        CHECK(Z(r, c) == (r == c ? 1.0 : 0.0));
                    }
                }
                continue;
            }
            ++accepted;
            double tmax = 0.0;
            for (size_t r = 0; r < S; ++r) {
                for (size_t c = 0; c < S; ++c) {
                    tmax = std::max(tmax, std::fabs(T0(r, c)));
                }
            }
            const Matrix<S, S> back = Z * Tm * Z.transpose();
            const Matrix<S, S> ztz = Z.transpose() * Z;
            for (size_t r = 0; r < S; ++r) {
                for (size_t c = 0; c < S; ++c) {
                    CHECK(std::fabs(back(r, c) - T0(r, c)) <= 100.0 * damp::eps<double>() * tmax);
                    CHECK(std::fabs(ztz(r, c) - (r == c ? 1.0 : 0.0)) <= 100.0 * damp::eps<double>());
                }
            }
            for (size_t r = Q; r < S; ++r) {
                for (size_t c = 0; c < Q; ++c) {
                    CHECK(Tm(r, c) == 0.0);
                }
            }
            CHECK(Tm(0, 0) == doctest::Approx(0.3 + gap).epsilon(1e-9)); // the Q block now leads
            if constexpr (Q == 2) {
                CHECK(Tm(0, 0) == Tm(1, 1));
                CHECK(Tm(0, 1) * Tm(1, 0) < 0.0);
            }
            if constexpr (P == 2) {
                CHECK(Tm(Q, Q) == Tm(Q + 1, Q + 1));
                CHECK(Tm(Q, Q + 1) * Tm(Q + 1, Q) < 0.0);
            }
        }
        CHECK(accepted > 300);
    }

    TEST_CASE("1×1 with 1×1") {
        check_swaps<1, 1>();
    }
    TEST_CASE("1×1 with 2×2") {
        check_swaps<1, 2>();
    }
    TEST_CASE("2×2 with 1×1") {
        check_swaps<2, 1>();
    }
    TEST_CASE("2×2 with 2×2") {
        check_swaps<2, 2>();
    }

    // A window that is not block upper triangular has no swap: the similarity built from the
    // Sylvester solution leaves the old coupling in the block that should vanish. The weak
    // stability test refuses it instead of storing a zero over it.
    TEST_CASE("a swap that would not leave a Schur form is refused and changes nothing") {
        unsigned seed = 5U;
        auto     T0 = window<2, 2>(0.3, 0.5, 0.8, 1.0, 1.0, seed);
        T0(2, 0) = 1e-3; // not a Schur form
        auto Tm = T0;
        auto Z = Matrix<4, 4>::identity();
        CHECK_FALSE(damp::detail::swap_schur_blocks<2, 2>(Tm, Z, 0));
        for (size_t r = 0; r < 4; ++r) {
            for (size_t c = 0; c < 4; ++c) {
                CHECK(Tm(r, c) == T0(r, c));
            }
        }
    }

    // The block test once used max(scale, 1) as its scale, an absolute floor: a complex pair
    // with entries of 1e-20 read as two 1×1 blocks. It is now relative to the block.
    TEST_CASE("a 2×2 block is recognized at any scale") {
        for (const double s : {1e-20, 1.0, 1e20}) {
            Matrix<3, 3> T{};
            T(0, 0) = s * 0.5;
            T(0, 1) = s * 2.0;
            T(1, 0) = -s * 0.5;
            T(1, 1) = s * 0.5;
            T(2, 2) = -s;
            T(0, 2) = s;
            CHECK(damp::detail::schur_subdiagonal(T, 0));
            CHECK_FALSE(damp::detail::schur_subdiagonal(T, 1));
        }
    }
}
