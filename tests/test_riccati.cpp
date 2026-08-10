// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "damp/backend.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/design/riccati.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
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
