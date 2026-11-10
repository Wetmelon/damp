// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

/**
 * @brief Tests for matrix function implementations
 *
 * Verifies exponentials, trigonometric functions, and power operations
 * on matrices.
 */
TEST_SUITE("Matrix Functions") {
    TEST_CASE("Matrix exponential - identity") {
        Matrix<2, 2> zeros = Matrix<2, 2>::zeros();

        auto exp_zeros = mat::expm(zeros).value();

        // exp(0) = I
        CHECK(doctest::Approx(exp_zeros(0, 0)).epsilon(1e-10) == 1.0);
        CHECK(doctest::Approx(exp_zeros(0, 1)).epsilon(1e-10) == 0.0);
        CHECK(doctest::Approx(exp_zeros(1, 0)).epsilon(1e-10) == 0.0);
        CHECK(doctest::Approx(exp_zeros(1, 1)).epsilon(1e-10) == 1.0);
    }

    TEST_CASE("Matrix exponential - diagonal") {
        Matrix<2, 2> A{{-1.0, 0.0}, {0.0, -2.0}};
        auto         exp_A = mat::expm(A).value();

        // exp(diag(-1,-2)) = diag(exp(-1), exp(-2))
        // Relaxed tolerance for Padé approximation
        CHECK(doctest::Approx(exp_A(0, 0)).epsilon(1e-3) == std::exp(-1.0));
        CHECK(doctest::Approx(exp_A(0, 1)).epsilon(1e-10) == 0.0);
        CHECK(doctest::Approx(exp_A(1, 0)).epsilon(1e-10) == 0.0);
        CHECK(doctest::Approx(exp_A(1, 1)).epsilon(1e-3) == std::exp(-2.0));
    }

    TEST_CASE("Matrix exponential - rotation") {
        // Rotation matrix: exp([0 -θ; θ 0]) = [cos(θ) -sin(θ); sin(θ) cos(θ)]
        double       theta = 0.5;
        Matrix<2, 2> A{{0.0, -theta}, {theta, 0.0}};
        auto         exp_A = mat::expm(A).value();

        // Relaxed tolerance for Padé approximation
        CHECK(doctest::Approx(exp_A(0, 0)).epsilon(1e-3) == std::cos(theta));
        CHECK(doctest::Approx(exp_A(0, 1)).epsilon(1e-3) == -std::sin(theta));
        CHECK(doctest::Approx(exp_A(1, 0)).epsilon(1e-3) == std::sin(theta));
        CHECK(doctest::Approx(exp_A(1, 1)).epsilon(1e-3) == std::cos(theta));
    }

    TEST_CASE("Matrix exponential - wide nilpotent") {
        // exp([[0, 1e16], [0, 0]]) = I + A. Scaling must not drop the off-diagonal.
        Matrix<2, 2> A{{0.0, 1.0e16}, {0.0, 0.0}};
        auto         E = mat::expm(A).value();
        CHECK(E(0, 0) == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(E(0, 1) == doctest::Approx(1.0e16).epsilon(1e-12));
        CHECK(E(1, 0) == doctest::Approx(0.0).epsilon(1e-12));
        CHECK(E(1, 1) == doctest::Approx(1.0).epsilon(1e-12));
    }

    TEST_CASE("Matrix exponential - wide triangular entries") {
        // Closed form: t_12 * (e^{λ1} − e^{λ2}) / (λ1 − λ2). Infinity-norm
        // scaling of this matrix drops the diagonal; the triangular rewrite keeps it.
        const double e1 = std::exp(-1.0);
        const double e2 = std::exp(-2.0);
        const double off = 1.0e16 * (e1 - e2);

        Matrix<2, 2> upper{{-1.0, 1.0e16}, {0.0, -2.0}};
        auto         Eu = mat::expm(upper).value();
        CHECK(Eu(0, 0) == doctest::Approx(e1).epsilon(1e-8));
        CHECK(Eu(0, 1) == doctest::Approx(off).epsilon(1e-8));
        CHECK(Eu(1, 0) == doctest::Approx(0.0).epsilon(1e-12));
        CHECK(Eu(1, 1) == doctest::Approx(e2).epsilon(1e-8));

        Matrix<2, 2> lower{{-1.0, 0.0}, {1.0e16, -2.0}};
        auto         El = mat::expm(lower).value();
        CHECK(El(0, 0) == doctest::Approx(e1).epsilon(1e-8));
        CHECK(El(0, 1) == doctest::Approx(0.0).epsilon(1e-12));
        CHECK(El(1, 0) == doctest::Approx(off).epsilon(1e-8));
        CHECK(El(1, 1) == doctest::Approx(e2).epsilon(1e-8));
    }

    TEST_CASE("Matrix exponential - diagonal similarity with wide entries") {
        // A = D M D^{-1}, D = diag(1e16, 1). Digits from scipy.linalg.expm.
        constexpr double d = 1.0e16;
        Matrix<2, 2>     M{{-0.5, 0.2}, {0.1, -0.3}};
        Matrix<2, 2>     A{{M(0, 0), M(0, 1) * d}, {M(1, 0) / d, M(1, 1)}};
        auto             E = mat::expm(A).value();

        const double ref00 = 0.613032341;
        const double ref01 = 1.34735335e15;
        const double ref10 = 6.73676677e-18;
        const double ref11 = 0.747767677;
        CHECK(E(0, 0) == doctest::Approx(ref00).epsilon(1e-8));
        CHECK(E(0, 1) == doctest::Approx(ref01).epsilon(1e-8));
        CHECK(std::abs(E(1, 0) - ref10) <= 1e-6 * std::abs(ref10));
        CHECK(E(1, 1) == doctest::Approx(ref11).epsilon(1e-8));
    }

    TEST_CASE("Matrix exponential - wide triangular entries, float") {
        Matrix<2, 2, float> A{{-1.0f, 1.0e8f}, {0.0f, -2.0f}};
        auto                E = mat::expm(A).value();
        const float         e1 = std::exp(-1.0f);
        const float         e2 = std::exp(-2.0f);
        const float         off = 1.0e8f * (e1 - e2);
        CHECK(E(0, 0) == doctest::Approx(e1).epsilon(1e-5));
        CHECK(E(0, 1) == doctest::Approx(off).epsilon(1e-5));
        CHECK(E(1, 0) == doctest::Approx(0.0f).epsilon(1e-6));
        CHECK(E(1, 1) == doctest::Approx(e2).epsilon(1e-5));
    }

    TEST_CASE("Matrix exponential - complex diagonal") {
        using C = damp::complex<double>;
        Matrix<2, 2, C> A{{C{0.0, 1.0}, C{0.0, 0.0}}, {C{0.0, 0.0}, C{0.0, -1.0}}};
        auto            E = mat::expm(A).value();
        const double    c = std::cos(1.0);
        const double    s = std::sin(1.0);
        CHECK(E(0, 0).real() == doctest::Approx(c).epsilon(1e-12));
        CHECK(E(0, 0).imag() == doctest::Approx(s).epsilon(1e-12));
        CHECK(E(1, 1).real() == doctest::Approx(c).epsilon(1e-12));
        CHECK(E(1, 1).imag() == doctest::Approx(-s).epsilon(1e-12));
        CHECK(E(0, 1).real() == doctest::Approx(0.0).epsilon(1e-12));
        CHECK(E(1, 0).imag() == doctest::Approx(0.0).epsilon(1e-12));
    }

    TEST_CASE("Matrix trig identity: sin^2 + cos^2 = I") {
        Matrix<2, 2> A{{0.3, 0.1}, {-0.1, 0.4}};

        auto sinA = mat::sin(A);
        auto cosA = mat::cos(A);
        auto sum = sinA * sinA + cosA * cosA;
        auto I = Matrix<2, 2>::identity();

        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                CHECK(doctest::Approx(sum(i, j)).epsilon(1e-6) == I(i, j));
            }
        }
    }

    TEST_CASE("Matrix sin/cos - diagonal, large norm") {
        // sin/cos of a diagonal matrix = diagonal of sin/cos of eigenvalues
        // Norm ≈ 10, exercises the scaling-and-doubling path
        Matrix<2, 2> A{{5.0, 0.0}, {0.0, 10.0}};

        auto sinA = mat::sin(A);
        auto cosA = mat::cos(A);

        CHECK(doctest::Approx(sinA(0, 0)).epsilon(1e-6) == std::sin(5.0));
        CHECK(doctest::Approx(sinA(0, 1)).epsilon(1e-10) == 0.0);
        CHECK(doctest::Approx(sinA(1, 0)).epsilon(1e-10) == 0.0);
        CHECK(doctest::Approx(sinA(1, 1)).epsilon(1e-6) == std::sin(10.0));

        CHECK(doctest::Approx(cosA(0, 0)).epsilon(1e-6) == std::cos(5.0));
        CHECK(doctest::Approx(cosA(0, 1)).epsilon(1e-10) == 0.0);
        CHECK(doctest::Approx(cosA(1, 0)).epsilon(1e-10) == 0.0);
        CHECK(doctest::Approx(cosA(1, 1)).epsilon(1e-6) == std::cos(10.0));
    }

    TEST_CASE("Matrix trig identity: large norm sin^2 + cos^2 = I") {
        Matrix<2, 2> A{{3.0, 1.5}, {-1.0, 4.0}};

        auto sinA = mat::sin(A);
        auto cosA = mat::cos(A);
        auto sum = sinA * sinA + cosA * cosA;
        auto I = Matrix<2, 2>::identity();

        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                CHECK(doctest::Approx(sum(i, j)).epsilon(1e-6) == I(i, j));
            }
        }
    }

    TEST_CASE("Matrix sqrt - identity") {
        auto I = Matrix<2, 2>::identity();
        auto sqrt_I = mat::sqrt(I);

        REQUIRE(sqrt_I.has_value());
        CHECK(doctest::Approx(sqrt_I.value()(0, 0)).epsilon(1e-6) == 1.0);
        CHECK(doctest::Approx(sqrt_I.value()(1, 1)).epsilon(1e-6) == 1.0);
    }

    TEST_CASE("Matrix power - integer") {
        Matrix<2, 2> A{{1.0, 1.0}, {0.0, 1.0}};

        auto A2_opt = mat::pow(A, 2);
        REQUIRE(A2_opt.has_value());
        auto A2 = A2_opt.value();
        auto A_sq = A * A;

        CHECK(doctest::Approx(A2(0, 0)).epsilon(1e-10) == A_sq(0, 0));
        CHECK(doctest::Approx(A2(0, 1)).epsilon(1e-10) == A_sq(0, 1));
        CHECK(doctest::Approx(A2(1, 0)).epsilon(1e-10) == A_sq(1, 0));
        CHECK(doctest::Approx(A2(1, 1)).epsilon(1e-10) == A_sq(1, 1));

        // Negative power of a singular matrix must fail closed (no fake I).
        Matrix<2, 2> singular{{1.0, 2.0}, {2.0, 4.0}};
        CHECK_FALSE(mat::pow(singular, -1).has_value());
    }

    TEST_CASE("Matrix hyperbolic identity: cosh^2 - sinh^2 = I") {
        Matrix<2, 2> A{{0.2, 0.1}, {0.1, 0.3}};

        auto sinhA = mat::sinh(A);
        auto coshA = mat::cosh(A);
        auto diff = coshA * coshA - sinhA * sinhA;
        auto I = Matrix<2, 2>::identity();

        // Relaxed tolerance for numerical precision through exp
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                CHECK(doctest::Approx(diff(i, j)).epsilon(1e-3) == I(i, j));
            }
        }
    }

    TEST_CASE("Matrix norms") {
        Matrix<2, 2> A{{1.0, 2.0}, {3.0, 4.0}};

        // One norm: max column sum = max(|1|+|3|, |2|+|4|) = max(4, 6) = 6
        CHECK(doctest::Approx(mat::one_norm(A)).epsilon(1e-10) == 6.0);

        // Infinity norm: max row sum = max(|1|+|2|, |3|+|4|) = max(3, 7) = 7
        CHECK(doctest::Approx(mat::infinity_norm(A)).epsilon(1e-10) == 7.0);

        // Frobenius norm: sqrt(1^2 + 2^2 + 3^2 + 4^2) = sqrt(1+4+9+16) = sqrt(30)
        CHECK(doctest::Approx(A.norm()).epsilon(1e-10) == std::sqrt(30.0));

        // Two norm (spectral norm): largest singular value of A = [[1,2],[3,4]]
        // sigma_max = sqrt(lambda_max(A^T A)) ≈ 5.46499
        CHECK(doctest::Approx(mat::two_norm(A)).epsilon(1e-4) == 5.4649857);
    }

    TEST_CASE("Matrix norms - identity") {
        auto I = Matrix<2, 2>::identity();

        CHECK(doctest::Approx(mat::one_norm(I)).epsilon(1e-10) == 1.0);
        CHECK(doctest::Approx(mat::infinity_norm(I)).epsilon(1e-10) == 1.0);
        CHECK(doctest::Approx(I.norm()).epsilon(1e-10) == std::sqrt(2.0));
        // Spectral norm of identity is 1 (all singular values are 1)
        CHECK(doctest::Approx(mat::two_norm(I)).epsilon(1e-10) == 1.0);
    }

    TEST_CASE("Matrix norms - zero matrix") {
        Matrix<2, 2> Z = Matrix<2, 2>::zeros();

        CHECK(doctest::Approx(mat::one_norm(Z)).epsilon(1e-10) == 0.0);
        CHECK(doctest::Approx(mat::infinity_norm(Z)).epsilon(1e-10) == 0.0);
        CHECK(doctest::Approx(Z.norm()).epsilon(1e-10) == 0.0);
        CHECK(doctest::Approx(mat::two_norm(Z)).epsilon(1e-10) == 0.0);
    }

    TEST_CASE("Matrix norms - complex matrices return real scalars via A^H A") {
        using C = damp::complex<double>;
        // Hermitian diagonal: A = diag(i, 2i) → singular values 1 and 2
        Matrix<2, 2, C> A{
            {C{0.0, 1.0}, C{0.0, 0.0}},
            {C{0.0, 0.0}, C{0.0, 2.0}},
        };

        // |i| + 0 = 1, 0 + |2i| = 2 → one/inf norms both 2
        CHECK(doctest::Approx(mat::one_norm(A)).epsilon(1e-12) == 2.0);
        CHECK(doctest::Approx(mat::infinity_norm(A)).epsilon(1e-12) == 2.0);
        // Frobenius: sqrt(1^2 + 2^2) = sqrt(5)
        CHECK(doctest::Approx(A.norm()).epsilon(1e-12) == std::sqrt(5.0));
        // Spectral: max |diag| = 2 (uses conjugate_transpose, not plain transpose)
        CHECK(doctest::Approx(mat::two_norm(A)).epsilon(1e-6) == 2.0);

        // Static type check: norms are real, not complex
        static_assert(std::is_same_v<decltype(mat::two_norm(A)), double>);
        static_assert(std::is_same_v<decltype(A.norm()), double>);
    }

    TEST_CASE("Matrix determinant") {
        // 1x1 matrix
        Matrix<1, 1> A1{{5.0}};
        CHECK(doctest::Approx(mat::det(A1)).epsilon(1e-10) == 5.0);

        // 2x2 matrix: det([1 2; 3 4]) = 1*4 - 2*3 = -2
        Matrix<2, 2> A2{{1.0, 2.0}, {3.0, 4.0}};
        CHECK(doctest::Approx(mat::det(A2)).epsilon(1e-10) == -2.0);

        // 3x3 matrix: det([1 2 3; 4 5 6; 7 8 9]) = 0 (singular)
        Matrix<3, 3> A3{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
        CHECK(doctest::Approx(mat::det(A3)).epsilon(1e-10) == 0.0);

        // Identity matrix
        auto I2 = Matrix<2, 2>::identity();
        CHECK(doctest::Approx(mat::det(I2)).epsilon(1e-10) == 1.0);

        auto I3 = Matrix<3, 3>::identity();
        CHECK(doctest::Approx(mat::det(I3)).epsilon(1e-10) == 1.0);
    }

    TEST_CASE("Matrix rank") {
        // Full rank matrices
        auto I2 = Matrix<2, 2>::identity();
        CHECK(mat::rank(I2) == 2);

        Matrix<2, 2> A2{{1.0, 2.0}, {3.0, 4.0}};
        CHECK(mat::rank(A2) == 2);

        // Rank deficient matrices
        Matrix<2, 2> singular{{1.0, 2.0}, {2.0, 4.0}}; // Second row is 2*first row
        CHECK(mat::rank(singular) == 1);

        Matrix<3, 3> rank1{{1.0, 2.0, 3.0}, {2.0, 4.0, 6.0}, {3.0, 6.0, 9.0}};
        CHECK(mat::rank(rank1) == 1);

        // Zero matrix
        Matrix<2, 2> Z = Matrix<2, 2>::zeros();
        CHECK(mat::rank(Z) == 0);
    }

    // Largest |A - B| element, for comparing matrix-function results.
    static constexpr auto max_abs_diff = []<size_t N>(const Matrix<N, N>& A, const Matrix<N, N>& B) {
        double worst = 0.0;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                worst = damp::max(worst, std::abs(A(i, j) - B(i, j)));
            }
        }
        return worst;
    };

    TEST_CASE("Matrix logarithm is the inverse of expm") {
        // log(exp(A)) == A for a non-trivial, well-conditioned A.
        Matrix<3, 3> A = {
            {0.1, 0.2, 0.0},
            {-0.1, 0.3, 0.1},
            {0.0, -0.2, 0.2},
        };
        auto log_exp = mat::log(mat::expm(A).value());
        REQUIRE(log_exp.has_value());
        CHECK(max_abs_diff(log_exp.value(), A) < 1e-9);

        // exp(log(B)) == B for an SPD matrix B (log well-defined).
        Matrix<2, 2> B = {
            {4.0, 1.0},
            {1.0, 3.0},
        };
        auto logB = mat::log(B);
        REQUIRE(logB.has_value());
        Matrix<2, 2> back = mat::expm(logB.value()).value();
        CHECK(max_abs_diff(back, B) < 1e-9);

        // log(identity) == 0
        auto logI = mat::log(Matrix<3, 3>::identity());
        REQUIRE(logI.has_value());
        CHECK(max_abs_diff(logI.value(), Matrix<3, 3>::zeros()) < 1e-12);

        // log(diagonal) == diagonal of logs
        Matrix<2, 2> D = {{2.0, 0.0}, {0.0, 5.0}};
        auto         logD_opt = mat::log(D);
        REQUIRE(logD_opt.has_value());
        Matrix<2, 2> logD = logD_opt.value();
        // logm uses inverse scaling-and-squaring; accuracy degrades for
        // eigenvalues far from 1 (more squaring steps), ~1e-7 here.
        CHECK(logD(0, 0) == doctest::Approx(std::log(2.0)).epsilon(1e-6));
        CHECK(logD(1, 1) == doctest::Approx(std::log(5.0)).epsilon(1e-6));

        // Singular / non-principal-branch matrices fail closed.
        Matrix<2, 2> singular{{1.0, 2.0}, {2.0, 4.0}};
        CHECK_FALSE(mat::log(singular).has_value());
        Matrix<2, 2> neg_eig{{-1.0, 0.0}, {0.0, 2.0}};
        CHECK_FALSE(mat::log(neg_eig).has_value());

        // MATLAB-style aliases share the fail-closed contract.
        CHECK(mat::log(B).has_value());
        CHECK_FALSE(mat::log(singular).has_value());
    }

    TEST_CASE("Matrix square root: sqrt(A)^2 == A") {
        Matrix<2, 2> A = {
            {5.0, 4.0},
            {1.0, 2.0},
        };
        auto S = mat::sqrt(A);
        REQUIRE(S.has_value());
        Matrix<2, 2> squared = S.value() * S.value();
        CHECK(max_abs_diff(squared, A) < 1e-8);

        // SPD 3x3
        Matrix<3, 3> B = {
            {4.0, 1.0, 0.0},
            {1.0, 3.0, 1.0},
            {0.0, 1.0, 2.0},
        };
        auto SB = mat::sqrt(B);
        REQUIRE(SB.has_value());
        CHECK(max_abs_diff(SB.value() * SB.value(), B) < 1e-7);

        // sqrt(identity) == identity
        auto SI = mat::sqrt(Matrix<3, 3>::identity());
        REQUIRE(SI.has_value());
        CHECK(max_abs_diff(SI.value(), Matrix<3, 3>::identity()) < 1e-10);
    }

    TEST_CASE("Matrix real-exponent power") {
        Matrix<2, 2> A = {
            {5.0, 4.0},
            {1.0, 2.0},
        };
        // A^0.5 squared == A
        auto half_opt = mat::pow(A, 0.5);
        REQUIRE(half_opt.has_value());
        Matrix<2, 2> half = half_opt.value();
        CHECK(max_abs_diff(half * half, A) < 1e-7);

        // A^2.0 == A*A
        auto sq_opt = mat::pow(A, 2.0);
        REQUIRE(sq_opt.has_value());
        CHECK(max_abs_diff(sq_opt.value(), A * A) < 1e-7);

        // A^1.0 == A, A^0.0 == I
        REQUIRE(mat::pow(A, 1.0).has_value());
        REQUIRE(mat::pow(A, 0.0).has_value());
        CHECK(max_abs_diff(mat::pow(A, 1.0).value(), A) < 1e-9);
        CHECK(max_abs_diff(mat::pow(A, 0.0).value(), Matrix<2, 2>::identity()) < 1e-9);

        // Non-integer power of a matrix with a real negative eigenvalue fails.
        Matrix<2, 2> neg{{-1.0, 0.0}, {0.0, 4.0}};
        CHECK_FALSE(mat::pow(neg, 0.5).has_value());
    }

    TEST_CASE("Matrix real-exponent power: a double exponent on a float or complex matrix is not truncated to an int") {
        // Regression: pow(A, T p) deduced T from both arguments, so pow(Matrix<float>, 0.5)
        // took the integer overload as pow(A, 0) and returned the identity.
        const Matrix<2, 2, float> Af{{4.0F, 0.0F}, {0.0F, 9.0F}};
        const auto                hf = mat::pow(Af, 0.5);
        REQUIRE(hf.has_value());
        CHECK((*hf)(0, 0) == doctest::Approx(2.0F));
        CHECK((*hf)(1, 1) == doctest::Approx(3.0F));

        using C = damp::complex<double>;
        const Matrix<2, 2, C> Ac{{C{4.0}, C{0.0}}, {C{0.0}, C{9.0}}};
        const auto            hc = mat::pow(Ac, 0.5);
        REQUIRE(hc.has_value());
        CHECK(damp::abs((*hc)(0, 0) - C{2.0}) < 1e-12);
        CHECK(damp::abs((*hc)(1, 1) - C{3.0}) < 1e-12);

        const auto cube = mat::pow(Af, 3U); // any integral exponent takes the integer path
        REQUIRE(cube.has_value());
        CHECK((*cube)(1, 1) == doctest::Approx(729.0F));
    }

    TEST_CASE("sincos / sin / cos consistency") {
        Matrix<3, 3> A = {
            {0.3, 0.1, 0.0},
            {-0.2, 0.4, 0.1},
            {0.0, -0.1, 0.5},
        };
        auto [S, C] = mat::sincos(A);

        // sin(A) and cos(A) agree with the paired primitive ...
        CHECK(max_abs_diff(mat::sin(A), S) < 1e-12);
        // ... and the standalone cos path matches sincos's cos exactly.
        CHECK(max_abs_diff(mat::cos(A), C) < 1e-12);

        // Pythagorean identity sin^2 + cos^2 == I
        Matrix<3, 3> id = (S * S) + (C * C);
        CHECK(max_abs_diff(id, Matrix<3, 3>::identity()) < 1e-9);

        // Larger norm exercises the scaling/doubling path
        Matrix<2, 2> B = {{2.0, 1.0}, {0.5, 3.0}};
        auto [SB, CB] = mat::sincos(B);
        CHECK(max_abs_diff(mat::cos(B), CB) < 1e-10);
        CHECK(max_abs_diff(mat::sin(B), SB) < 1e-10);
        CHECK(max_abs_diff((SB * SB) + (CB * CB), Matrix<2, 2>::identity()) < 1e-8);
    }

    TEST_CASE("expm refuses a non-finite input instead of returning the identity") {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        CHECK_FALSE(mat::expm(Matrix<2, 2>{{nan, 0.0}, {1.0, -1.0}}).has_value());
        CHECK_FALSE(mat::expm(Matrix<2, 2>{{0.0, std::numeric_limits<double>::infinity()}, {0.0, 0.0}}).has_value());
        CHECK(mat::expm(Matrix<2, 2>{{0.0, 1.0}, {-1.0, 0.0}}).has_value());
    }

    TEST_CASE("det of a nearly singular 5×5 and 6×6 keeps its small value") {
        // lu_decomposition's cancellation floor refused these, and det returned exactly 0.
        Matrix<5, 5> A5 = Matrix<5, 5>::identity();
        A5(0, 1) = 1.0;
        A5(1, 0) = 1.0;
        A5(1, 1) = 1.0 + 1e-13;
        CHECK(mat::det(A5) == doctest::Approx(1e-13).epsilon(1e-3));
        Matrix<6, 6, float> A6 = Matrix<6, 6, float>::identity();
        A6(0, 1) = 1.0f;
        A6(1, 0) = 1.0f;
        A6(1, 1) = 1.0f + 1e-5f;
        CHECK(mat::det(A6) == doctest::Approx(1e-5).epsilon(0.1));
        // Sign of the row permutation: a cyclic shift of the 5×5 identity has det +1.
        Matrix<5, 5> P{};
        for (size_t i = 0; i < 5; ++i) {
            P(i, (i + 1) % 5) = 2.0;
        }
        CHECK(mat::det(P) == doctest::Approx(32.0));
        CHECK(mat::det(Matrix<5, 5>{}) == 0.0);
    }

    TEST_CASE("pow with a negative exponent powers A⁻¹, keeping κ(A) instead of κ(A)^|p|") {
        const double       c = std::cos(0.7);
        const double       sn = std::sin(0.7);
        const Matrix<2, 2> Q{{c, -sn}, {sn, c}};
        const double       lambda = 1e-4;
        const Matrix<2, 2> A = Q * Matrix<2, 2>{{1.0, 0.0}, {0.0, lambda}} * Q.t();
        // A is rounded: its small eigenvalue moves by ε‖A‖, a relative κ(A)·ε, so A^p is no closer
        // to ref than |p|·κ(A)·ε. Inverting A^|p| instead would cost κ(A)^|p|·ε: 2e-4 and 2 here.
        const double kappa = 1.0 / lambda;
        for (const int p : {-3, -4}) {
            const Matrix<2, 2> ref = Q * Matrix<2, 2>{{1.0, 0.0}, {0.0, std::pow(lambda, p)}} * Q.t();
            const auto         got = mat::pow(A, p);
            REQUIRE(got.has_value());
            CHECK((got.value() - ref).norm() <= 2.0 * static_cast<double>(-p) * kappa * damp::eps<double>() * ref.norm());
        }
    }

    TEST_CASE("pow rounds a near-integer exponent from either side to the integer path") {
        // diag(−1, 4) has no real logarithm, so only the integer path succeeds.
        const Matrix<2, 2> A{{-1.0, 0.0}, {0.0, 4.0}};
        for (const double p : {std::nextafter(3.0, 0.0), std::nextafter(3.0, 4.0)}) {
            const auto got = mat::pow(A, p);
            REQUIRE(got.has_value());
            CHECK(got.value()(0, 0) == doctest::Approx(-1.0));
            CHECK(got.value()(1, 1) == doctest::Approx(64.0));
        }
        const auto inv2 = mat::pow(A, std::nextafter(-2.0, 0.0));
        REQUIRE(inv2.has_value());
        CHECK(inv2.value()(1, 1) == doctest::Approx(1.0 / 16.0));
    }

    TEST_CASE_TEMPLATE("sinh keeps its relative accuracy for a small argument", T, float, double) {
        // (exp(A) − exp(−A))/2 cancelled to eps/‖A‖: 78 % error at ‖A‖ = 1e-6 in float.
        const T               s = std::is_same_v<T, float> ? T(1e-6) : T(1e-12);
        const Matrix<2, 2, T> A{{T(0), s}, {T(-2) * s, T(0.5) * s}};
        const Matrix<2, 2, T> ref = A + (A * A * A * (T(1) / T(6)));
        const Matrix<2, 2, T> got = mat::sinh(A);
        CHECK(Matrix<2, 2, T>(got - ref).norm() <= T(16) * eps<T>() * ref.norm());
        const Matrix<2, 2, T> D{{T(1e-10), T(0)}, {T(0), T(2e-10)}};
        CHECK(mat::sinh(D)(0, 0) == doctest::Approx(1e-10).epsilon(static_cast<double>(T(8) * eps<T>())));
        // cosh² − sinh² = I at a larger norm (the doubling path)
        const Matrix<2, 2, T> B{{T(1.5), T(0.4)}, {T(-0.3), T(2.0)}};
        const auto            sh = mat::sinh(B);
        const auto            ch = mat::cosh(B);
        CHECK(Matrix<2, 2, T>((ch * ch) - (sh * sh) - Matrix<2, 2, T>::identity()).norm() <= T(1e3) * eps<T>() * Matrix<2, 2, T>(ch * ch).norm());
    }

    TEST_CASE_TEMPLATE("sinh and cosh keep expm's accuracy above ‖A‖∞ = 1/2", T, float, double) {
        // The double-angle recovery from A/2ˢ loses a factor 2 per step: 2.2–2.9·λ·eps on the
        // symmetric A and 2350–10200 eps on the non-normal one, against 0.3–0.75·λ·eps and 30 eps
        // by (expm(A) ∓ expm(−A))/2. Under -ffast-math each translation unit may round expm's
        // squaring its own way, and the test runner links one of those copies: up to 1.5·λ·eps
        // in float. The bound sits between that and the recovery's 2.2.
        const auto tol = static_cast<double>(eps<T>());
        // Q = H/2 (H the 4×4 Hadamard matrix) is orthogonal with entries ±1/2, so A = Q·D·Qᵀ is exact.
        const Matrix<4, 4> Q = Matrix<4, 4>{{1.0, 1.0, 1.0, 1.0}, {1.0, -1.0, 1.0, -1.0}, {1.0, 1.0, -1.0, -1.0}, {1.0, -1.0, -1.0, 1.0}} * 0.5;
        for (const double m : {1.0, 4.0}) {
            const double lambda = 10.0 * m; // spectrum λ·(1, −0.8, 0.3, 0.1)
            const double ev[4] = {lambda, -8.0 * m, 3.0 * m, m};
            Matrix<4, 4> D{};
            Matrix<4, 4> Ds{};
            Matrix<4, 4> Dc{};
            for (size_t i = 0; i < 4; ++i) {
                D(i, i) = ev[i];
                Ds(i, i) = std::sinh(ev[i]);
                Dc(i, i) = std::cosh(ev[i]);
            }
            const Matrix<4, 4, T> A(Matrix<4, 4>(Q * D * Q.t()));
            const Matrix<4, 4>    rs = Q * Ds * Q.t();
            const Matrix<4, 4>    rc = Q * Dc * Q.t();
            CHECK(Matrix<4, 4>(Matrix<4, 4>(mat::sinh(A)) - rs).norm() <= 2.0 * lambda * tol * rs.norm());
            CHECK(Matrix<4, 4>(Matrix<4, 4>(mat::cosh(A)) - rc).norm() <= 2.0 * lambda * tol * rc.norm());
        }
        // Non-normal, ‖A‖_F = 43: A = S·U·S⁻¹ with S = [[1,0,0],[1,1,0],[−1,2,1]] and U upper
        // triangular; the reference is S·f(U)·S⁻¹, f(U) by Parlett's recurrence.
        const Matrix<3, 3> U{{1.5, 4.0, -2.0}, {0.0, -2.0, 3.0}, {0.0, 0.0, 0.5}};
        const Matrix<3, 3> S{{1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {-1.0, 2.0, 1.0}};
        const Matrix<3, 3> S_inv{{1.0, 0.0, 0.0}, {-1.0, 1.0, 0.0}, {3.0, -2.0, 1.0}};
        const auto         parlett = [&](double f0, double f1, double f2) {
            Matrix<3, 3> F{{f0, 0.0, 0.0}, {0.0, f1, 0.0}, {0.0, 0.0, f2}};
            F(0, 1) = U(0, 1) * (f0 - f1) / (U(0, 0) - U(1, 1));
            F(1, 2) = U(1, 2) * (f1 - f2) / (U(1, 1) - U(2, 2));
            F(0, 2) = ((U(0, 2) * (f0 - f2)) + (F(0, 1) * U(1, 2)) - (U(0, 1) * F(1, 2))) / (U(0, 0) - U(2, 2));
            return Matrix<3, 3>(S * F * S_inv);
        };
        const Matrix<3, 3, T> A(Matrix<3, 3>(S * U * S_inv));
        const Matrix<3, 3>    rs = parlett(std::sinh(1.5), std::sinh(-2.0), std::sinh(0.5));
        const Matrix<3, 3>    rc = parlett(std::cosh(1.5), std::cosh(-2.0), std::cosh(0.5));
        CHECK(Matrix<3, 3>(Matrix<3, 3>(mat::sinh(A)) - rs).norm() <= 200.0 * tol * rs.norm());
        CHECK(Matrix<3, 3>(Matrix<3, 3>(mat::cosh(A)) - rc).norm() <= 200.0 * tol * rc.norm());
    }

    TEST_CASE("sin and cos return for an infinite entry or an overflowing row sum") {
        // The scaling loop halved an infinite norm forever.
        const double       inf = std::numeric_limits<double>::infinity();
        const Matrix<2, 2> A{{inf, 0.0}, {0.0, 1.0}};
        static_cast<void>(mat::cos(A));
        static_cast<void>(mat::sin(A));
        const Matrix<2, 2> big{{1e308, 1e308}, {0.0, 0.0}};
        static_cast<void>(mat::sincos(big));
        static_cast<void>(mat::cosh(big));
        CHECK(mat::detail::trig_scale_exponent(A) == 0);
        CHECK(mat::detail::trig_scale_exponent(big) == 0);
        CHECK(mat::detail::trig_scale_exponent(Matrix<2, 2>{{1e300, 0.0}, {0.0, 0.0}}) == 998);
    }

    TEST_CASE_TEMPLATE("sqrtm on well-posed matrices with κ(A) up to 1e9: residual at roundoff, float included", T, float, double) {
        // Denman–Beavers inverted iterates of condition κ(A) and returned nullopt in float.
        std::mt19937                           gen(7);
        std::uniform_real_distribution<double> u(-1.0, 1.0);
        for (int trial = 0; trial < 50; ++trial) {
            Matrix<4, 4> S{};
            for (size_t i = 0; i < 4; ++i) {
                for (size_t j = 0; j < 4; ++j) {
                    S(i, j) = u(gen) + ((i == j) ? 2.0 : 0.0);
                }
            }
            Matrix<4, 4> D{};
            D(0, 0) = std::sqrt(std::is_same_v<T, float> ? 1e-4 : 1e-9);
            D(1, 1) = 1.0;
            D(2, 2) = 1.2;
            D(3, 3) = 1.2;
            if (trial % 2 == 1) {
                D(2, 3) = 0.7; // a complex pair
                D(3, 2) = -0.7;
            }
            const Matrix<4, 4>    X = S * D * S.inverse().value();
            const Matrix<4, 4, T> A(Matrix<4, 4>(X * X));
            const auto            R = mat::sqrt(A);
            REQUIRE(R.has_value());
            CHECK(Matrix<4, 4, T>((R.value() * R.value()) - A).norm() <= T(64) * eps<T>() * R.value().norm() * R.value().norm());
        }
        // SPD, eigenvalues 1e-3, 1, 1e3, rotated: κ = 1e6
        const Matrix<3, 3, T> U{{T(0.6), T(0.8), T(0)}, {T(-0.48), T(0.36), T(0.8)}, {T(0.64), T(-0.48), T(0.6)}};
        const Matrix<3, 3, T> A = U * Matrix<3, 3, T>{{T(1e3), T(0), T(0)}, {T(0), T(1), T(0)}, {T(0), T(0), T(1e-3)}} * U.t();
        const auto            R = mat::sqrt(A);
        REQUIRE(R.has_value());
        const Matrix<3, 3, T> ref = U * Matrix<3, 3, T>{{std::sqrt(T(1e3)), T(0), T(0)}, {T(0), T(1), T(0)}, {T(0), T(0), std::sqrt(T(1e-3))}} * U.t();
        CHECK(Matrix<3, 3, T>(R.value() - ref).norm() <= T(1e4) * eps<T>() * ref.norm());
        CHECK(mat::log(A).has_value());
        // A rotation by π/2 has the root rotation by π/4; a negative eigenvalue has no real principal root.
        const auto rot = mat::sqrt(Matrix<2, 2, T>{{T(0), T(-1)}, {T(1), T(0)}});
        REQUIRE(rot.has_value());
        CHECK(rot.value()(0, 0) == doctest::Approx(std::sqrt(0.5)).epsilon(1e-6));
        CHECK(rot.value()(1, 0) == doctest::Approx(std::sqrt(0.5)).epsilon(1e-6));
        CHECK_FALSE(mat::sqrt(Matrix<2, 2, T>{{T(-1), T(0)}, {T(0), T(2)}}).has_value());
    }

    TEST_CASE("sqrtm and logm take a complex matrix") {
        // The real Schur method has no complex form here; a complex A takes Denman–Beavers.
        using C = damp::complex<double>;
        const Matrix<2, 2, C> A{{C{4.0, 1.0}, C{1.0, 0.0}}, {C{0.0, 0.5}, C{3.0, -1.0}}};
        const auto            X = mat::sqrt(A);
        REQUIRE(X.has_value());
        CHECK(((X.value() * X.value()) - A).norm() <= 64.0 * eps<double>() * A.norm());
        const auto L = mat::log(A);
        REQUIRE(L.has_value());
        const auto E = mat::expm(L.value());
        REQUIRE(E.has_value());
        CHECK((E.value() - A).norm() <= 1e-12 * A.norm());
        const C zero{0.0, 0.0};
        CHECK_FALSE(mat::sqrt(Matrix<2, 2, C>{{C{-1.0, 0.0}, zero}, {zero, C{2.0, 0.0}}}).has_value());
    }

    TEST_CASE_TEMPLATE("logm of a graded matrix keeps the accuracy of the ungraded one", T, float, double) {
        // log(D·E·D⁻¹) = D·log(E)·D⁻¹ (D powers of two, exact). Unbalanced, the stop test
        // ‖A_s − I‖∞ ≤ 1/2 saw the grading and took 20 square roots: 3.6e6 eps in double and
        // a 10 % error in float, reported as success.
        const Matrix<3, 3, T> M{{T(-0.375), T(0.25), T(0.125)}, {T(0.25), T(-0.5), T(0.375)}, {T(-0.125), T(0.25), T(-0.25)}};
        const Matrix<3, 3, T> D{{T(1) / T(1024), T(0), T(0)}, {T(0), T(1), T(0)}, {T(0), T(0), T(1024)}};
        const Matrix<3, 3, T> D_inv{{T(1024), T(0), T(0)}, {T(0), T(1), T(0)}, {T(0), T(0), T(1) / T(1024)}};
        const auto            E = mat::expm(M);
        REQUIRE(E.has_value());
        // Balancing takes the grading out in full, up to a common shift, so the two logs run the
        // same arithmetic. Stopping short (a spread of 16 left) cost float 43–82 eps, by host.
        const auto e = mat::detail::balance_exponents(Matrix<3, 3, T>(D * E.value() * D_inv));
        CHECK(e[1] - e[0] == 10);
        CHECK(e[2] - e[1] == 10);
        const auto L = mat::log(Matrix<3, 3, T>(D * E.value() * D_inv));
        const auto L_ref = mat::log(E.value());
        REQUIRE(L.has_value());
        REQUIRE(L_ref.has_value());
        CHECK(Matrix<3, 3, T>((D_inv * L.value() * D) - L_ref.value()).norm() <= T(64) * eps<T>() * L_ref.value().norm());
    }

    TEST_CASE("logm takes one Schur form, so float keeps a few eps over several square roots") {
        // 2I + E, E uniform in [−1/2, 1/2]; the float logm against the double logm of the same
        // matrix. A Schur form per square root repeated its roundoff at every step: up to 46–82 eps.
        std::uint32_t state = 12345u;
        const auto    next = [&state] {
            state = (1664525u * state) + 1013904223u;
            return (static_cast<float>(state >> 8u) / 16777216.0f) - 0.5f;
        };
        float worst = 0.0f;
        for (int t = 0; t < 10; ++t) {
            Matrix<6, 6, float> A{};
            for (size_t i = 0; i < 6; ++i) {
                for (size_t j = 0; j < 6; ++j) {
                    A(i, j) = (i == j ? 2.0f : 0.0f) + next();
                }
            }
            const auto Lf = mat::log(A);
            const auto Ld = mat::log(A.as<double>());
            REQUIRE(Lf.has_value());
            REQUIRE(Ld.has_value());
            const double err = (Lf.value().as<double>() - Ld.value()).norm() / Ld.value().norm();
            worst = damp::max(worst, static_cast<float>(err));
        }
        CHECK(worst <= 24.0f * eps<float>());
    }

    TEST_CASE_TEMPLATE("logm refuses a non-finite entry", T, float, double) {
        // A NaN norm passes the stop test ‖A_s − I‖∞ > 1/2 (the comparison is false), which skipped
        // the square roots; |NaN + 0i| is 0, so a complex NaN entry was not caught by its magnitude.
        const T nan = std::numeric_limits<T>::quiet_NaN();
        CHECK_FALSE(mat::log(Matrix<2, 2, T>{{T(2), nan}, {T(0.1), T(3)}}).has_value());
        using C = damp::complex<T>;
        const Matrix<2, 2, C> Ac{{C{T(2), T(0)}, C{nan, T(0)}}, {C{T(0.1), T(0)}, C{T(3), T(0)}}};
        CHECK_FALSE(mat::log(Ac).has_value());
    }

    TEST_CASE_TEMPLATE("sinh and cosh carry an overflow into the result when expm refuses", T, float, double) {
        // The column sums overflow (expm refuses) while the row sums are finite; the series then
        // scales A by 2⁻¹²⁹ (float). A scalar 2⁻ˢ is subnormal and -ffast-math flushed it to zero,
        // so sinh came back 0 and cosh I.
        const T               big = static_cast<T>(0.6) * std::numeric_limits<T>::max();
        const Matrix<2, 2, T> A{{big, T(0)}, {big, T(0)}};
        REQUIRE_FALSE(mat::expm(A).has_value());
        CHECK_FALSE(mat::all_of(mat::sinh(A), [](T v) { return damp::isfinite(v); }));
        CHECK_FALSE(mat::all_of(mat::cosh(A), [](T v) { return damp::isfinite(v); }));
    }
}

TEST_SUITE("Matrix functions at any scale") {
    // √(s·A) = √s·√A and log(s·A) = log(s)·I + log A. The real Schur form behind both squared
    // the entries of its 2×2 blocks, so a complex pair past √max (1.8e19 in float) or below √min
    // came back wrong and the residual test refused the root.
    TEST_CASE_TEMPLATE("sqrtm and logm of s·A", T, float, double) {
        const Matrix<3, 3, T> A{{T(4), T(1), T(0.5)}, {T(-1), T(3), T(0.25)}, {T(0), T(0.5), T(2)}};
        const auto            R0 = mat::sqrt(A);
        const auto            L0 = mat::log(A);
        REQUIRE(R0.has_value());
        REQUIRE(L0.has_value());
        const double tol = 256.0 * static_cast<double>(eps<T>());
        // Even powers of two: √s is exact, so the reference √s·√A carries no rounding of its own.
        // Named arrays, not a ?: of two initializer lists: those lists' arrays die with the
        // expression, before the loop reads them (GCC 11 read garbage).
        const std::array<int, 4> float_exponents{-100, -64, 64, 100};
        const std::array<int, 4> double_exponents{-1000, -520, 520, 1000};
        for (const int e : (damp::is_same_v<T, float> ? float_exponents : double_exponents)) {
            CAPTURE(e);
            const T    s = mat::detail::scale_by_pow2(T{1}, e);
            const auto R = mat::sqrt(Matrix<3, 3, T>(A * s));
            REQUIRE(R.has_value());
            const T root_s = mat::detail::scale_by_pow2(T{1}, e / 2);
            CHECK(static_cast<double>(Matrix<3, 3, T>((R.value() / root_s) - R0.value()).norm()) <= tol * static_cast<double>(R0.value().norm()));
            const auto L = mat::log(Matrix<3, 3, T>(A * s));
            REQUIRE(L.has_value());
            const T log_s = static_cast<T>(static_cast<double>(e) * std::log(2.0));
            CHECK(static_cast<double>(Matrix<3, 3, T>(L.value() - (Matrix<3, 3, T>::identity() * log_s) - L0.value()).norm()) <= tol * static_cast<double>(damp::abs(log_s) + L0.value().norm()));
        }
    }
}
