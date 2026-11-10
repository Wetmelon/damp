// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>
#include <limits>
#include <random>

#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using namespace damp::mat;
using damp::complex;

TEST_SUITE("cholesky") {
    TEST_CASE("Cholesky decomposition - basic 2x2 SPD matrix") {
        // Test matrix: [[4, 2], [2, 3]] - known to be SPD
        Matrix<2, 2, double> A = {{4.0, 2.0}, {2.0, 3.0}};

        auto L_opt = cholesky(A);
        REQUIRE(L_opt.has_value());

        const auto& L = L_opt.value();

        // Check that L is lower triangular
        CHECK(L(0, 1) == 0.0);

        // Check that A = L * L^T
        auto L_LT = L * L.transpose();
        CHECK(L_LT(0, 0) == doctest::Approx(A(0, 0)));
        CHECK(L_LT(0, 1) == doctest::Approx(A(0, 1)));
        CHECK(L_LT(1, 0) == doctest::Approx(A(1, 0)));
        CHECK(L_LT(1, 1) == doctest::Approx(A(1, 1)));
    }

    TEST_CASE("Cholesky decomposition - 3x3 SPD matrix") {
        // Test matrix: [[4, 2, 1], [2, 5, 3], [1, 3, 6]]
        Matrix<3, 3, double> A = {
            {4.0, 2.0, 1.0},
            {2.0, 5.0, 3.0},
            {1.0, 3.0, 6.0},
        };

        auto L_opt = cholesky(A);
        REQUIRE(L_opt.has_value());

        const auto& L = L_opt.value();

        // Check that L is lower triangular
        CHECK(L(0, 1) == 0.0);
        CHECK(L(0, 2) == 0.0);
        CHECK(L(1, 2) == 0.0);

        // Check that A = L * L^T
        auto L_LT = L * L.transpose();
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                CHECK(L_LT(i, j) == doctest::Approx(A(i, j)));
            }
        }
    }

    TEST_CASE("Cholesky decomposition - identity matrix") {
        Matrix<3, 3, double> A = Matrix<3, 3, double>::identity();

        auto L_opt = cholesky(A);
        REQUIRE(L_opt.has_value());

        const auto& L = L_opt.value();

        // For identity, L should also be identity
        CHECK(L(0, 0) == doctest::Approx(1.0));
        CHECK(L(1, 1) == doctest::Approx(1.0));
        CHECK(L(2, 2) == doctest::Approx(1.0));
        CHECK(L(0, 1) == 0.0);
        CHECK(L(0, 2) == 0.0);
        CHECK(L(1, 2) == 0.0);
    }

    TEST_CASE("Cholesky decomposition - non-SPD matrix (negative eigenvalue)") {
        // Matrix with negative eigenvalue: [[1, 2], [2, 1]]
        Matrix<2, 2, double> A = {{1.0, 2.0}, {2.0, 1.0}};

        auto L_opt = cholesky(A);
        CHECK(!L_opt.has_value());
    }

    TEST_CASE("Cholesky decomposition - non-SPD matrix (not positive definite)") {
        // Matrix that's not SPD: [[0, 0], [0, 0]]
        Matrix<2, 2, double> A = {{0.0, 0.0}, {0.0, 0.0}};

        auto L_opt = cholesky(A);
        CHECK(!L_opt.has_value());
    }

    TEST_CASE("Cholesky decomposition - non-symmetric matrix") {
        // Non-symmetric matrix: [[4, 1], [2, 3]]
        // Cholesky doesn't check symmetry, so it may succeed or fail
        Matrix<2, 2, double> A = {{4.0, 1.0}, {2.0, 3.0}};

        auto L_opt = cholesky(A);
        // The algorithm doesn't enforce symmetry, so result is implementation-dependent
        // Just check that it either succeeds or fails gracefully
        (void)L_opt; // Suppress unused variable warning
    }

    TEST_CASE("Forward substitution - basic test") {
        // Lower triangular matrix L
        Matrix<3, 3, double> L = {{2.0, 0.0, 0.0}, {1.0, 3.0, 0.0}, {0.5, 0.5, 4.0}};

        // Right-hand side b
        ColVec<3, double> b = {6.0, 8.0, 9.0};

        // Expected solution: L * x = b
        // x(0) = 6/2 = 3
        // x(1) = (8 - 1*3)/3 = 5/3 ≈ 1.6667
        // x(2) = (9 - 0.5*3 - 0.5*5/3)/4 = (9 - 1.5 - 5/6)/4 ≈ 1.7917
        auto x = forward_substitute(L, b);

        CHECK(x(0) == doctest::Approx(3.0));
        CHECK(x(1) == doctest::Approx(5.0 / 3.0));
        CHECK(x(2) == doctest::Approx((9.0 - 1.5 - (5.0 / 6.0)) / 4.0));
    }

    TEST_CASE("Backward substitution - basic test") {
        // Upper triangular matrix U (L^T)
        Matrix<3, 3, double> U = {{2.0, 1.0, 0.5}, {0.0, 3.0, 0.5}, {0.0, 0.0, 4.0}};

        // Right-hand side b
        ColVec<3, double> b = {7.0, 5.5, 4.0};

        // Expected solution: U * x = b (treating U as L in backward_substitute)
        // x(2) = 4/4 = 1
        // x(1) = 5.5/3 ≈ 1.8333 (since U(2,1) = 0)
        // x(0) = 7/2 = 3.5 (since U(1,0) = 0, U(2,0) = 0)
        auto x = backward_substitute_transpose(U, b);

        CHECK(x(0) == doctest::Approx(3.5));
        CHECK(x(1) == doctest::Approx(5.5 / 3.0));
        CHECK(x(2) == doctest::Approx(1.0));
    }

    TEST_CASE("Cholesky solve - single RHS") {
        // SPD matrix A
        Matrix<2, 2, double> A = {{4.0, 2.0}, {2.0, 3.0}};

        // Right-hand side B (single column)
        Matrix<2, 1, double> B = {{8.0}, {7.0}};

        // Expected solution: A * x = B
        // x should be approximately [1.25, 1.5]
        auto X_opt = cholesky_solve(A, B);
        REQUIRE(X_opt.has_value());

        const auto& X = X_opt.value();
        CHECK(X(0, 0) == doctest::Approx(1.25));
        CHECK(X(1, 0) == doctest::Approx(1.5));
    }

    TEST_CASE("Cholesky solve - multiple RHS") {
        // SPD matrix A
        Matrix<3, 3, double> A = {{4.0, 2.0, 1.0}, {2.0, 5.0, 3.0}, {1.0, 3.0, 6.0}};

        // Right-hand side B (two columns)
        Matrix<3, 2, double> B = {{7.0, 10.0}, {9.0, 13.0}, {11.0, 15.0}};

        auto X_opt = cholesky_solve(A, B);
        REQUIRE(X_opt.has_value());

        const auto& X = X_opt.value();

        // Verify A * X = B for both columns
        auto AX = A * X;
        CHECK(AX(0, 0) == doctest::Approx(B(0, 0)));
        CHECK(AX(0, 1) == doctest::Approx(B(0, 1)));
        CHECK(AX(1, 0) == doctest::Approx(B(1, 0)));
        CHECK(AX(1, 1) == doctest::Approx(B(1, 1)));
        CHECK(AX(2, 0) == doctest::Approx(B(2, 0)));
        CHECK(AX(2, 1) == doctest::Approx(B(2, 1)));
    }

    TEST_CASE("Cholesky solve - non-SPD matrix") {
        // Non-SPD matrix
        Matrix<2, 2, double> A = {{1.0, 2.0}, {2.0, 1.0}};

        Matrix<2, 1, double> B = {{3.0}, {4.0}};

        auto X_opt = cholesky_solve(A, B);
        CHECK(!X_opt.has_value());
    }

    TEST_CASE("Cholesky decomposition - constexpr") {
        // Test that Cholesky works at compile time
        constexpr Matrix<2, 2, double> A = {{4.0, 2.0}, {2.0, 3.0}};

        constexpr auto L_opt = cholesky(A);
        static_assert(L_opt.has_value(), "Cholesky should succeed at compile time");

        // We can't access .value() in constexpr context for static_assert,
        // but we can verify the function compiles and returns a value
    }

    TEST_CASE("Forward substitution - constexpr") {
        constexpr Matrix<2, 2, double> L = {{2.0, 0.0}, {1.0, 1.0}};
        constexpr ColVec<2, double>    b = {4.0, 1.0};

        constexpr auto x = forward_substitute(L, b);

        static_assert(x(0) > 1.9 && x(0) < 2.1, "x(0) should be ~2.0");
        static_assert(x(1) > -1.1 && x(1) < -0.9, "x(1) should be ~-1.0");
    }

    TEST_CASE("Cholesky solve - float precision") {
        Matrix<2, 2, float> A = {{4.0f, 2.0f}, {2.0f, 3.0f}};
        Matrix<2, 1, float> B = {{8.0f}, {7.0f}};

        auto X_opt = cholesky_solve(A, B);
        REQUIRE(X_opt.has_value());

        const auto& X = X_opt.value();
        CHECK(X(0, 0) == doctest::Approx(1.25));
        CHECK(X(1, 0) == doctest::Approx(1.5));
    }

    TEST_CASE("LU decomposition - basic 2x2 matrix") {
        Matrix<2, 2, double> A = {{4.0, 2.0}, {2.0, 3.0}};

        auto lu_opt = lu_decomposition(A);
        REQUIRE(lu_opt.has_value());

        const auto& [L, U, piv] = lu_opt.value();

        // Check that L is lower triangular with unit diagonal
        CHECK(L(0, 0) == doctest::Approx(1.0));
        CHECK(L(1, 0) == doctest::Approx(0.5));
        CHECK(L(0, 1) == 0.0);
        CHECK(L(1, 1) == doctest::Approx(1.0));

        // Check that U is upper triangular
        CHECK(U(0, 0) == doctest::Approx(4.0));
        CHECK(U(0, 1) == doctest::Approx(2.0));
        CHECK(U(1, 0) == 0.0);
        CHECK(U(1, 1) == doctest::Approx(2.0));

        // Check that A = L * U
        auto LU = L * U;
        CHECK(LU(0, 0) == doctest::Approx(A(0, 0)));
        CHECK(LU(0, 1) == doctest::Approx(A(0, 1)));
        CHECK(LU(1, 0) == doctest::Approx(A(1, 0)));
        CHECK(LU(1, 1) == doctest::Approx(A(1, 1)));
    }

    TEST_CASE("LU decomposition - 3x3 matrix") {
        Matrix<3, 3, double> A = {
            {4.0, 1.0, 1.0},
            {2.0, 3.0, 3.0},
            {2.0, 1.0, 9.0}
        };

        auto lu_opt = lu_decomposition(A);
        REQUIRE(lu_opt.has_value());

        const auto& [L, U, piv] = lu_opt.value();

        // Check that L is lower triangular with unit diagonal
        CHECK(L(0, 0) == doctest::Approx(1.0));
        CHECK(L(1, 0) == doctest::Approx(0.5));
        CHECK(L(2, 0) == doctest::Approx(0.5));
        CHECK(L(1, 1) == doctest::Approx(1.0));
        CHECK(L(2, 1) == doctest::Approx(0.2));
        CHECK(L(2, 2) == doctest::Approx(1.0));
        CHECK(L(0, 1) == 0.0);
        CHECK(L(0, 2) == 0.0);
        CHECK(L(1, 2) == 0.0);

        // Check that U is upper triangular
        CHECK(U(0, 0) == doctest::Approx(4.0));
        CHECK(U(0, 1) == doctest::Approx(1.0));
        CHECK(U(0, 2) == doctest::Approx(1.0));
        CHECK(U(1, 1) == doctest::Approx(2.5));
        CHECK(U(1, 2) == doctest::Approx(2.5));
        CHECK(U(2, 2) == doctest::Approx(8.0));
        CHECK(U(1, 0) == 0.0);
        CHECK(U(2, 0) == 0.0);
        CHECK(U(2, 1) == 0.0);

        // Check that A = L * U
        auto LU = L * U;
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                CHECK(LU(i, j) == doctest::Approx(A(i, j)));
            }
        }
    }

    TEST_CASE("LU decomposition - singular matrix") {
        // Singular matrix: [[1, 2], [2, 4]] - second row is 2*first row
        Matrix<2, 2, double> A = {{1.0, 2.0}, {2.0, 4.0}};

        auto lu_opt = lu_decomposition(A);
        CHECK(!lu_opt.has_value());
    }

    TEST_CASE("LU decomposition - float precision") {
        Matrix<2, 2, float> A = {{4.0f, 2.0f}, {2.0f, 3.0f}};

        auto lu_opt = lu_decomposition(A);
        REQUIRE(lu_opt.has_value());

        const auto& [L, U, piv] = lu_opt.value();

        // Check that A = L * U
        auto LU = L * U;
        CHECK(LU(0, 0) == doctest::Approx(A(0, 0)));
        CHECK(LU(0, 1) == doctest::Approx(A(0, 1)));
        CHECK(LU(1, 0) == doctest::Approx(A(1, 0)));
        CHECK(LU(1, 1) == doctest::Approx(A(1, 1)));
    }

    TEST_CASE("LU decomposition - complex matrix") {
        using Cplx = complex<double>;
        Matrix<2, 2, Cplx> A = {
            {Cplx{4.0, 1.0}, Cplx{2.0, 0.5}},
            {Cplx{2.0, 0.5}, Cplx{3.0, 2.0}}
        };

        auto lu_opt = lu_decomposition(A);
        REQUIRE(lu_opt.has_value());

        const auto& [L, U, piv] = lu_opt.value();

        // Check that L is lower triangular with unit diagonal
        CHECK(L(0, 0).real() == doctest::Approx(1.0));
        CHECK(L(0, 0).imag() == doctest::Approx(0.0));
        CHECK(L(1, 1).real() == doctest::Approx(1.0));
        CHECK(L(1, 1).imag() == doctest::Approx(0.0));
        CHECK(L(0, 1).real() == doctest::Approx(0.0));
        CHECK(L(0, 1).imag() == doctest::Approx(0.0));

        // Check that U is upper triangular
        CHECK(U(1, 0).real() == doctest::Approx(0.0));
        CHECK(U(1, 0).imag() == doctest::Approx(0.0));

        // Check that A = L * U
        auto LU = L * U;
        CHECK(LU(0, 0).real() == doctest::Approx(A(0, 0).real()));
        CHECK(LU(0, 0).imag() == doctest::Approx(A(0, 0).imag()));
        CHECK(LU(0, 1).real() == doctest::Approx(A(0, 1).real()));
        CHECK(LU(0, 1).imag() == doctest::Approx(A(0, 1).imag()));
        CHECK(LU(1, 0).real() == doctest::Approx(A(1, 0).real()));
        CHECK(LU(1, 0).imag() == doctest::Approx(A(1, 0).imag()));
        CHECK(LU(1, 1).real() == doctest::Approx(A(1, 1).real()));
        CHECK(LU(1, 1).imag() == doctest::Approx(A(1, 1).imag()));
    }

    TEST_CASE("Solve function - SPD matrix (uses Cholesky)") {
        Matrix<2, 2, double> A = {{4.0, 2.0}, {2.0, 3.0}};
        Matrix<2, 1, double> B = {{8.0}, {7.0}};

        auto X_opt = solve(A, B);
        REQUIRE(X_opt.has_value());

        const auto& X = X_opt.value();
        CHECK(X(0, 0) == doctest::Approx(1.25));
        CHECK(X(1, 0) == doctest::Approx(1.5));

        // Verify A * X = B
        auto AX = A * X;
        CHECK(AX(0, 0) == doctest::Approx(B(0, 0)));
        CHECK(AX(1, 0) == doctest::Approx(B(1, 0)));
    }

    TEST_CASE("Solve function - non-SPD matrix (uses LU)") {
        Matrix<2, 2, double> A = {{1.0, 2.0}, {3.0, 4.0}}; // Not SPD
        Matrix<2, 1, double> B = {{5.0}, {11.0}};

        auto X_opt = solve(A, B);
        REQUIRE(X_opt.has_value());

        const auto& X = X_opt.value();
        CHECK(X(0, 0) == doctest::Approx(1.0));
        CHECK(X(1, 0) == doctest::Approx(2.0));

        // Verify A * X = B
        auto AX = A * X;
        CHECK(AX(0, 0) == doctest::Approx(B(0, 0)));
        CHECK(AX(1, 0) == doctest::Approx(B(1, 0)));
    }

    TEST_CASE("Solve function - singular matrix") {
        Matrix<2, 2, double> A = {{1.0, 2.0}, {2.0, 4.0}}; // Singular
        Matrix<2, 1, double> B = {{3.0}, {6.0}};

        auto X_opt = solve(A, B);
        CHECK(!X_opt.has_value());
    }

    TEST_CASE("Solve function - float precision") {
        Matrix<2, 2, float> A = {{4.0f, 2.0f}, {2.0f, 3.0f}};
        Matrix<2, 1, float> B = {{8.0f}, {7.0f}};

        auto X_opt = solve(A, B);
        REQUIRE(X_opt.has_value());

        const auto& X = X_opt.value();
        CHECK(X(0, 0) == doctest::Approx(1.25f));
        CHECK(X(1, 0) == doctest::Approx(1.5f));
    }

    TEST_CASE("Solve function - complex matrix") {
        using Cplx = complex<double>;
        Matrix<2, 2, Cplx> A = {
            {Cplx{4.0, 1.0}, Cplx{2.0, 0.5}},
            {Cplx{2.0, 0.5}, Cplx{3.0, 2.0}}
        };
        Matrix<2, 1, Cplx> B = {
            {Cplx{8.0, 2.0}},
            {Cplx{7.0, 1.5}}
        };

        auto X_opt = solve(A, B);
        REQUIRE(X_opt.has_value());

        const auto& X = X_opt.value();

        // Verify A * X = B
        auto AX = A * X;
        CHECK(AX(0, 0).real() == doctest::Approx(B(0, 0).real()));
        CHECK(AX(0, 0).imag() == doctest::Approx(B(0, 0).imag()));
        CHECK(AX(1, 0).real() == doctest::Approx(B(1, 0).real()));
        CHECK(AX(1, 0).imag() == doctest::Approx(B(1, 0).imag()));
    }

    TEST_CASE("LU decomposition - constexpr") {
        constexpr Matrix<2, 2, double> A = {{4.0, 2.0}, {2.0, 3.0}};

        constexpr auto lu_opt = lu_decomposition(A);
        static_assert(lu_opt.has_value(), "LU decomposition should succeed at compile time");
    }

    TEST_CASE("Solve function - constexpr") {
        constexpr Matrix<2, 2, double> A = {{4.0, 2.0}, {2.0, 3.0}};
        constexpr Matrix<2, 1, double> B = {{8.0}, {7.0}};

        constexpr auto X_opt = solve(A, B);
        static_assert(X_opt.has_value(), "Solve should succeed at compile time");

        // We can't access .value() in constexpr context for static_assert,
        // but we can verify the function compiles
    }

    TEST_CASE("Triangular solve - lower triangular 3x3") {
        // L * x = b  →  x = [1, 2, 3]
        // L = [[2,0,0],[1,3,0],[4,2,5]], b = L * [1,2,3]
        Matrix<3, 3, double> A = {{2.0, 0.0, 0.0}, {1.0, 3.0, 0.0}, {4.0, 2.0, 5.0}};
        Matrix<3, 1, double> B = {{2.0}, {7.0}, {23.0}};

        auto L = A.lower_triangle();
        auto X_opt = solve(L, B);
        REQUIRE(X_opt.has_value());

        const auto& X = X_opt.value();
        CHECK(X(0, 0) == doctest::Approx(1.0));
        CHECK(X(1, 0) == doctest::Approx(2.0));
        CHECK(X(2, 0) == doctest::Approx(3.0));

        // Verify L * X = B
        auto LX = A * X;
        CHECK(LX(0, 0) == doctest::Approx(B(0, 0)));
        CHECK(LX(1, 0) == doctest::Approx(B(1, 0)));
        CHECK(LX(2, 0) == doctest::Approx(B(2, 0)));
    }

    TEST_CASE("Triangular solve - upper triangular 3x3") {
        // U * x = b  →  x = [1, 2, 3]
        // U = [[2,1,4],[0,3,2],[0,0,5]], b = U * [1,2,3]
        Matrix<3, 3, double> A = {{2.0, 1.0, 4.0}, {0.0, 3.0, 2.0}, {0.0, 0.0, 5.0}};
        Matrix<3, 1, double> B = {{16.0}, {12.0}, {15.0}};

        auto U = A.upper_triangle();
        auto X_opt = solve(U, B);
        REQUIRE(X_opt.has_value());

        const auto& X = X_opt.value();
        CHECK(X(0, 0) == doctest::Approx(1.0));
        CHECK(X(1, 0) == doctest::Approx(2.0));
        CHECK(X(2, 0) == doctest::Approx(3.0));

        // Verify U * X = B
        auto UX = A * X;
        CHECK(UX(0, 0) == doctest::Approx(B(0, 0)));
        CHECK(UX(1, 0) == doctest::Approx(B(1, 0)));
        CHECK(UX(2, 0) == doctest::Approx(B(2, 0)));
    }

    TEST_CASE("Triangular solve - singular lower triangular") {
        Matrix<2, 2, double> A = {{0.0, 0.0}, {1.0, 2.0}}; // zero on diagonal
        Matrix<2, 1, double> B = {{1.0}, {3.0}};

        auto L = A.lower_triangle();
        auto X_opt = solve(L, B);
        CHECK(!X_opt.has_value());
    }

    TEST_CASE("Triangular solve - singular upper triangular") {
        Matrix<2, 2, double> A = {{1.0, 2.0}, {0.0, 0.0}}; // zero on diagonal
        Matrix<2, 1, double> B = {{3.0}, {0.0}};

        auto U = A.upper_triangle();
        auto X_opt = solve(U, B);
        CHECK(!X_opt.has_value());
    }

    TEST_CASE("Triangular solve - constexpr lower triangular") {
        static constexpr Matrix<2, 2, double> A = {{2.0, 0.0}, {1.0, 3.0}};
        static constexpr Matrix<2, 1, double> B = {{4.0}, {7.0}};

        constexpr auto L = A.lower_triangle();
        constexpr auto X_opt = solve(L, B);
        static_assert(X_opt.has_value(), "Lower triangular solve should succeed at compile time");
    }

    TEST_CASE("Triangular solve - constexpr upper triangular") {
        static constexpr Matrix<2, 2, double> A = {{2.0, 1.0}, {0.0, 3.0}};
        static constexpr Matrix<2, 1, double> B = {{5.0}, {6.0}};

        constexpr auto U = A.upper_triangle();
        constexpr auto X_opt = solve(U, B);
        static_assert(X_opt.has_value(), "Upper triangular solve should succeed at compile time");
    }

    TEST_CASE("Cholesky decomposition - complex Hermitian positive-definite") {
        using Cplx = complex<double>;
        // Hermitian PD matrix: A = [[2, 1-j], [1+j, 3]]
        // Diagonal is real, A(0,1) == conj(A(1,0))
        Matrix<2, 2, Cplx> A{
            {Cplx{2.0, 0.0}, Cplx{1.0, -1.0}},
            {Cplx{1.0, 1.0}, Cplx{3.0, 0.0}}
        };

        auto L_opt = cholesky(A);
        REQUIRE(L_opt.has_value());

        // Verify A = L * Lᴴ
        const auto& L = L_opt.value();
        auto        reconstructed = L * L.conjugate_transpose();
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                CHECK(damp::abs(reconstructed(i, j) - A(i, j)) < 1e-10);
            }
        }
    }

    TEST_CASE("Cholesky solve - complex Hermitian system") {
        using Cplx = complex<double>;
        Matrix<2, 2, Cplx> A{
            {Cplx{4.0, 0.0}, Cplx{1.0, -1.0}},
            {Cplx{1.0, 1.0}, Cplx{3.0, 0.0}}
        };
        Matrix<2, 1, Cplx> B{
            {Cplx{3.0, 1.0}},
            {Cplx{2.0, -1.0}}
        };

        auto X_opt = cholesky_solve(A, B);
        REQUIRE(X_opt.has_value());

        // Verify A * X ≈ B
        auto AX = A * X_opt.value();
        for (size_t i = 0; i < 2; ++i) {
            CHECK(damp::abs(AX(i, 0) - B(i, 0)) < 1e-10);
        }
    }

    TEST_CASE("Hermitian check - 1x1 complex with imaginary part rejected") {
        using Cplx = complex<double>;
        // A 1x1 complex matrix with nonzero imaginary part is NOT Hermitian
        Matrix<1, 1, Cplx> A{{Cplx{1.0, 0.5}}};
        CHECK_FALSE(is_symmetric_or_hermitian(A));
    }

    TEST_CASE("LU solve - complex non-Hermitian system") {
        using Cplx = complex<double>;
        // Non-Hermitian complex matrix — must use LU path
        Matrix<2, 2, Cplx> A{
            {Cplx{1.0, 1.0}, Cplx{2.0, 0.0}},
            {Cplx{0.0, 1.0}, Cplx{3.0, -1.0}}
        };
        Matrix<2, 1, Cplx> B{
            {Cplx{5.0, 1.0}},
            {Cplx{3.0, 2.0}}
        };

        auto X_opt = lu_solve(A, B);
        REQUIRE(X_opt.has_value());

        auto AX = A * X_opt.value();
        for (size_t i = 0; i < 2; ++i) {
            CHECK(damp::abs(AX(i, 0) - B(i, 0)) < 1e-10);
        }
    }
}

TEST_SUITE("full_qr") {
    // Largest |A - B| element.
    static constexpr auto max_abs_diff = []<size_t R, size_t C>(const Matrix<R, C>& A, const Matrix<R, C>& B) {
        double worst = 0.0;
        for (size_t i = 0; i < R; ++i) {
            for (size_t j = 0; j < C; ++j) {
                worst = damp::abs(A(i, j) - B(i, j)) > worst ? damp::abs(A(i, j) - B(i, j)) : worst;
            }
        }
        return worst;
    };

    TEST_CASE("full_qr reconstructs A with an orthogonal full Q") {
        SUBCASE("square 3x3") {
            Matrix<3, 3> A = {
                {12.0, -51.0, 4.0},
                {6.0, 167.0, -68.0},
                {-4.0, 24.0, -41.0},
            };
            auto qr = mat::full_qr(A);
            // Q·R == A
            CHECK(max_abs_diff(qr.Q * qr.R, A) < 1e-9);
            // QᵀQ == I (full orthogonality)
            CHECK(max_abs_diff(qr.Q.transpose() * qr.Q, Matrix<3, 3>::identity()) < 1e-9);
            // R upper triangular
            CHECK(damp::abs(qr.R(1, 0)) < 1e-9);
            CHECK(damp::abs(qr.R(2, 0)) < 1e-9);
            CHECK(damp::abs(qr.R(2, 1)) < 1e-9);
        }

        SUBCASE("tall 4x2 — trailing Q columns span the orthogonal complement") {
            Matrix<4, 2> B = {
                {1.0, 2.0},
                {3.0, 4.0},
                {5.0, 6.0},
                {7.0, 9.0},
            };
            auto qr = mat::full_qr(B);
            CHECK(max_abs_diff(qr.Q * qr.R, B) < 1e-9);
            CHECK(max_abs_diff(qr.Q.transpose() * qr.Q, Matrix<4, 4>::identity()) < 1e-9);

            // Columns 2..3 of Q (the complement U₁) are orthogonal to range(B):
            // U₁ᵀ·B == 0.
            for (size_t c = 2; c < 4; ++c) {
                for (size_t j = 0; j < 2; ++j) {
                    double d = 0.0;
                    for (size_t i = 0; i < 4; ++i) {
                        d += qr.Q(i, c) * B(i, j);
                    }
                    CHECK(damp::abs(d) < 1e-9);
                }
            }
        }
    }

    TEST_CASE("full_qr is constexpr-evaluable") {
        constexpr auto ortho_ok = []() consteval {
            Matrix<3, 2> B = {{1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
            auto         qr = mat::full_qr(B);
            // U₁ (column 2) ⟂ range(B)
            double d0 = 0.0;
            double d1 = 0.0;
            for (size_t i = 0; i < 3; ++i) {
                d0 += qr.Q(i, 2) * B(i, 0);
                d1 += qr.Q(i, 2) * B(i, 1);
            }
            return damp::abs(d0) < 1e-12 && damp::abs(d1) < 1e-12;
        };
        static_assert(ortho_ok(), "full_qr complement must be constexpr-correct");
        CHECK(ortho_ok());
    }

    TEST_CASE("full_qr factors complex matrices with a unitary Q") {
        using Cplx = complex<double>;

        // Largest |·| deviation of a complex matrix from a reference.
        auto worst = []<size_t R, size_t Cc>(const Matrix<R, Cc, Cplx>& X, const Matrix<R, Cc, Cplx>& Y) {
            double w = 0.0;
            for (size_t i = 0; i < R; ++i) {
                for (size_t j = 0; j < Cc; ++j) {
                    w = damp::abs(X(i, j) - Y(i, j)) > w ? damp::abs(X(i, j) - Y(i, j)) : w;
                }
            }
            return w;
        };

        SUBCASE("square 3x3") {
            const Matrix<3, 3, Cplx> A = {
                {{2, 1}, {0, -1}, {1, 0}},
                {{1, 1}, {3, 0}, {0, 2}},
                {{0, 1}, {1, -1}, {2, 1}},
            };
            const auto qr = mat::full_qr(A);
            CHECK(worst(qr.Q * qr.R, A) < 1e-9);                                                    // A = Q·R
            CHECK(worst(qr.Q.conjugate_transpose() * qr.Q, Matrix<3, 3, Cplx>::identity()) < 1e-9); // Qᴴ·Q = I
            CHECK(damp::abs(qr.R(1, 0)) < 1e-9);                                                    // R upper-triangular
            CHECK(damp::abs(qr.R(2, 0)) < 1e-9);
            CHECK(damp::abs(qr.R(2, 1)) < 1e-9);
        }

        SUBCASE("tall 5x3 — trailing Q columns span ker(Aᴴ)") {
            // This is the [A−λI | B]ᴴ shape that arbitrary pole placement uses.
            const Matrix<5, 3, Cplx> A = {
                {{1, 1}, {2, -1}, {0, 2}},
                {{0, 2}, {1, 0}, {3, -1}},
                {{3, 0}, {-1, 1}, {2, 2}},
                {{-1, 1}, {0, 3}, {1, 0}},
                {{2, -2}, {1, 1}, {0, 1}},
            };
            const auto qr = mat::full_qr(A);
            CHECK(worst(qr.Q * qr.R, A) < 1e-9);
            CHECK(worst(qr.Q.conjugate_transpose() * qr.Q, Matrix<5, 5, Cplx>::identity()) < 1e-9);
            // Columns 3..4 of Q are orthogonal to range(A): Qᴴ_⊥·A = 0.
            for (size_t c = 3; c < 5; ++c) {
                for (size_t j = 0; j < 3; ++j) {
                    Cplx d{0, 0};
                    for (size_t i = 0; i < 5; ++i) {
                        d += damp::conj(qr.Q(i, c)) * A(i, j);
                    }
                    CHECK(damp::abs(d) < 1e-9);
                }
            }
        }

        SUBCASE("leading-zero pivot exercises the phase guard") {
            const Matrix<3, 3, Cplx> A = {
                {{0, 0}, {1, 0}, {2, 1}},
                {{1, 1}, {0, 0}, {1, 0}},
                {{2, 0}, {1, 1}, {0, 0}},
            };
            const auto qr = mat::full_qr(A);
            CHECK(worst(qr.Q * qr.R, A) < 1e-9);
            CHECK(worst(qr.Q.conjugate_transpose() * qr.Q, Matrix<3, 3, Cplx>::identity()) < 1e-9);
        }
    }
}

TEST_SUITE("scale-invariant singularity and symmetry") {
    TEST_CASE("solve does not mistake a tiny non-symmetric matrix for symmetric") {
        // An absolute 1e-12 symmetry tolerance sent this to Cholesky, which reads
        // only the lower triangle and returned a residual of 0.5.
        const Matrix<2, 2> A{{2e-13, 1e-13}, {0.0, 2e-13}};
        const Matrix<2, 1> b{{1.0}, {1.0}};
        CHECK_FALSE(is_symmetric_or_hermitian(A));
        const auto x = solve(A, b);
        REQUIRE(x.has_value());
        const auto r = A * x.value() - b;
        CHECK(damp::abs(r(0, 0)) < 1e-12);
        CHECK(damp::abs(r(1, 0)) < 1e-12);
    }

    TEST_CASE("lu_solve accepts a well-posed system at any scale") {
        for (const double s : {1e-20, 1.0, 1e20}) {
            const Matrix<2, 2> A{{2.0 * s, 1.0 * s}, {1.0 * s, 3.0 * s}};
            const Matrix<2, 1> b{{s}, {2.0 * s}};
            const auto         x = lu_solve(A, b);
            REQUIRE(x.has_value());
            CHECK(x.value()(0, 0) == doctest::Approx(0.2));
            CHECK(x.value()(1, 0) == doctest::Approx(0.6));
        }
    }

    TEST_CASE("lu_decomposition rejects roundoff-singular input at any scale") {
        for (const double s : {1e-20, 1.0, 1e20}) {
            const Matrix<3, 3> A = Matrix<3, 3>{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}} * s;
            CHECK_FALSE(lu_decomposition(A).has_value());
        }
    }

    TEST_CASE("a float upper-triangular matrix with a 1e19–1e20 diagonal is not taken as symmetric") {
        // √(|a₀₀||a₁₁|) formed as a product overflowed float, so every pair passed and solve()
        // ran Cholesky on the lower triangle only.
        const Matrix<2, 2, float> A{{2e19f, 1e19f}, {0.0f, 3e19f}};
        const Matrix<2, 1, float> b{{1e19f}, {1e19f}};
        CHECK_FALSE(is_symmetric_or_hermitian(A));
        const Matrix<2, 2, float> big{{1e20f, 1e19f}, {-1e19f, 1e20f}};
        CHECK_FALSE(is_symmetric_or_hermitian(big));
        const auto x = solve(A, b);
        REQUIRE(x.has_value());
        CHECK(x.value()(0, 0) == doctest::Approx(1.0f / 3.0f).epsilon(1e-5));
        CHECK(x.value()(1, 0) == doctest::Approx(1.0f / 3.0f).epsilon(1e-5));
        const auto inv = A.inverse();
        REQUIRE(inv.has_value());
        CHECK(inv.value()(0, 1) == doctest::Approx(-1.0f / 6.0f * 1e-19f).epsilon(1e-5));
    }

    TEST_CASE("solve takes LU for a nearly symmetric matrix, so its residual stays at roundoff") {
        // An asymmetry of 1e-13 relative passes is_symmetric_or_hermitian (N·64·eps), but
        // Cholesky reads only the lower triangle and would solve a different matrix.
        Matrix<3, 3> A{{4.0, 1.0, 0.5}, {1.0, 3.0, 0.2}, {0.5, 0.2, 2.0}};
        A(0, 1) *= 1.0 + 1e-13;
        A(1, 2) *= 1.0 - 1e-13;
        REQUIRE(is_symmetric_or_hermitian(A));
        const Matrix<3, 1> b{{1.0}, {-2.0}, {0.5}};
        const auto         x = solve(A, b);
        REQUIRE(x.has_value());
        const Matrix<3, 1> r = (A * x.value()) - b;
        CHECK(r.norm() <= 16.0 * eps<double>() * A.norm() * x.value().norm());
    }

    TEST_CASE("lu_decomposition factors a row-graded well-conditioned matrix") {
        // diag(1, 1e-3, 1e-22)·B with κ₂(B) ≈ 6: the largest entry of column 1 after the first
        // step is cancellation roundoff in row 1, while row 2 holds the genuine pivot.
        const Matrix<3, 3> B{{3.0, 1.1, 0.0}, {1.0, 1.1 / 3.0, 1.0}, {0.0, 1.0, 1.0}};
        const Matrix<3, 3> D{{1.0, 0.0, 0.0}, {0.0, 1e-3, 0.0}, {0.0, 0.0, 1e-22}};
        const Matrix<3, 3> A = D * B;
        const auto         lu = lu_decomposition(A);
        REQUIRE(lu.has_value());
        const auto& [L, U, piv] = lu.value();
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < i; ++j) {
                CHECK(U(i, j) == 0.0);
                CHECK(damp::abs(L(i, j)) <= 1.0);
            }
        }
        const Matrix<3, 1> b = D * Matrix<3, 1>{{1.0}, {2.0}, {3.0}};
        const auto         x = solve(A, b);
        REQUIRE(x.has_value());
        const auto xb = solve(B, Matrix<3, 1>{{1.0}, {2.0}, {3.0}});
        REQUIRE(xb.has_value());
        CHECK((x.value() - xb.value()).norm() <= 1e-12 * xb.value().norm());
    }

    TEST_CASE("lu_decomposition refuses a row-graded matrix whose digit-carrying entry outranks a clean pivot") {
        // diag(1e30, s, 1)·B with det B ≈ −1. After the first step U(1, 1) ≈ s·100ε is 50ε of
        // its magnitude: above the 32ε roundoff bound, so it is a candidate and beats the clean 1
        // in row 2 by absolute size, but below the N·64·eps floor. Pivot choice by absolute size
        // is not row-scale invariant: with the row scaled to s = 1e10 it no longer outranks the 1.
        const double d = 100.0 * eps<double>();
        const auto   graded = [d](double s) {
            return Matrix<3, 3>{{1e30, 1e30, 1e30}, {s, s * (1.0 + d), 2.0 * s}, {0.0, 1.0, 3.0}};
        };
        CHECK_FALSE(lu_decomposition(graded(1e15)).has_value());
        CHECK_FALSE(lu_decomposition(graded(1e20)).has_value());
        CHECK(lu_decomposition(graded(1e10)).has_value());
    }

    TEST_CASE_TEMPLATE("lu_decomposition eliminates a lost entry smaller than the pivot instead of dropping it", T, float, double) {
        // Symmetric indefinite, κ ≈ 6. After the first step U(1, 1) = δ is exact but below the
        // N·64·eps floor; taking it as zero left ‖PA − LU‖ at 48·eps·‖A‖ and a forward error of 118 eps.
        const T               d = T(128) * eps<T>();
        const Matrix<3, 3, T> A{{T(1), T(1), T(0)}, {T(1), T(1) + d, T(1)}, {T(0), T(1), T(1)}};
        const auto            lu = lu_decomposition(A);
        REQUIRE(lu.has_value());
        const auto& [L, U, piv] = lu.value();
        Matrix<3, 3, T> PA{};
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                PA(i, j) = A(piv[i], j);
            }
        }
        CHECK(Matrix<3, 3, T>(PA - (L * U)).norm() <= T(4) * eps<T>() * A.norm());
        const Matrix<3, 1, T> x{{T(1)}, {T(2)}, {T(3)}};
        const auto            xs = lu_solve(A, Matrix<3, 1, T>(A * x));
        REQUIRE(xs.has_value());
        CHECK(Matrix<3, 1, T>(xs.value() - x).norm() <= T(8) * eps<T>() * x.norm());
    }

    TEST_CASE_TEMPLATE("a NaN or ±Inf entry is not symmetric, and solve, inverse and pow(A, −1) refuse it", T, float, double) {
        // The bound comparison is false for NaN and for ±Inf (whose scale is infinite), so a
        // non-finite strict-upper entry passed the symmetry test and Cholesky, which reads only
        // the lower triangle, returned inv([[2, 0.1], [0.1, 3]]) as a success.
        for (const T x : {std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::infinity(), -std::numeric_limits<T>::infinity()}) {
            const Matrix<2, 2, T> A{{T(2), x}, {T(0.1), T(3)}};
            CHECK_FALSE(is_symmetric_or_hermitian(A));
            CHECK_FALSE(A.inverse().has_value());
            CHECK_FALSE(solve(A, Matrix<2, 2, T>::identity()).has_value());
            CHECK_FALSE(mat::pow(A, -1).has_value());
            // An infinite pivot read as roundoff of its infinite magnitude, and LU zeroed it.
            const Matrix<2, 2, T> D{{x, T(0.1)}, {T(0.1), T(3)}};
            CHECK_FALSE(solve(D, Matrix<2, 2, T>::identity()).has_value());
        }
        // |NaN + 0i| is 0, so the complex test checks the parts.
        using C = complex<T>;
        const Matrix<2, 2, C> Ac{{C{T(2), T(0)}, C{std::numeric_limits<T>::quiet_NaN(), T(0)}}, {C{T(0.1), T(0)}, C{T(3), T(0)}}};
        CHECK_FALSE(is_symmetric_or_hermitian(Ac));
        CHECK_FALSE(Ac.inverse().has_value());
    }

    TEST_CASE("lu_decomposition refuses a column whose skipped entry still carries digits; solve falls back to Cholesky") {
        // Near-symmetric SPD (asymmetry 4.9·eps, κ∞ = 1.4·10¹⁴). At step 2 the largest entry, 203 eps
        // of its magnitude, sits below the N·64·eps floor but far above roundoff; zeroing it to pivot
        // on an entry 6 times smaller returned x with a forward error of 4, as a success.
        const Matrix<4, 4> A{{0x1.e63ce298bcaddp-4, -0x1.b1bcdbf0afa1ap-3, -0x1.d1f0c07dcf05cp-3, -0x1.0e385c0ac0bcfp-6}, {-0x1.b1bcdbf0afa22p-3, 0x1.9aafbdf5f078cp-2, 0x1.bfad96223139dp-2, 0x1.e28341504c77ap-6}, {-0x1.d1f0c07dcf065p-3, 0x1.bfad9622313a5p-2, 0x1.e9acf007d391ap-2, 0x1.0339f699ea5bbp-5}, {-0x1.0e385c0ac0bcap-6, 0x1.e28341504c771p-6, 0x1.0339f699ea5b6p-5, 0x1.2c58a3c703785p-9}};
        const Matrix<4, 1> b{{0x1.7990adc263206p+0}, {-0x1.75c8513ad6794p-1}, {-0x1.77b70205d24d0p-2}, {0x1.487948145a374p+1}};
        // The exact solution of this (A, b), by rational elimination.
        const Matrix<4, 1> x_exact{{0x1.3fa332e1260e7p+45}, {0x1.0e53d9ede8353p+45}, {-0x1.93e08520ee188p+44}, {0x1.93bea8c899cb6p+47}};
        CHECK_FALSE(lu_decomposition(A).has_value());
        CHECK_FALSE(lu_solve(A, b).has_value());
        // Symmetric to N·64·eps, so solve() takes Cholesky once LU refuses (κ∞·eps = 0.03).
        const auto x = solve(A, b);
        REQUIRE(x.has_value());
        CHECK((x.value() - x_exact).norm() <= 1e-2 * x_exact.norm());
    }

    TEST_CASE("backward_substitute_transpose solves Lᴴx = b for a complex Cholesky factor") {
        using C = damp::complex<double>;
        const Matrix<2, 2, C> A{{C{2.0, 0.0}, C{1.0, -1.0}}, {C{1.0, 1.0}, C{3.0, 0.0}}};
        const ColVec<2, C>    b{C{1.0, 0.0}, C{0.0, 1.0}};
        const auto            L = cholesky(A);
        REQUIRE(L.has_value());
        const auto x = backward_substitute_transpose(*L, forward_substitute(*L, b));
        const auto r = (A * x) - b;
        CHECK(r(0).abs() <= 1e-14);
        CHECK(r(1).abs() <= 1e-14);
    }
}

namespace {

/// No diagonal of R was stored as zero, the mark qr_decompose leaves on a column it drops as roundoff.
template<typename T, size_t N, size_t M>
bool no_column_dropped(const QRDecomposition<T, N, M>& qr) {
    for (size_t j = 0; j < M; ++j) {
        if (qr.R(j, j) == T{0}) {
            return false;
        }
    }
    return true;
}

} // namespace

TEST_SUITE("rank, two_norm and QR are scale-invariant") {
    TEST_CASE("rank uses relative SVD threshold at any scale") {
        for (const double s : {1e-20, 1.0, 1e20}) {
            const Matrix<3, 3> A = Matrix<3, 3>{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}} * s;
            CHECK(rank(A) == 2);
            CHECK(rank(Matrix<3, 3>::identity() * s) == 3);
        }
        CHECK(rank(Matrix<2, 3>{}) == 0);
    }

    TEST_CASE("two_norm when the all-ones vector is a null vector") {
        CHECK(two_norm(Matrix<2, 2>{{1.0, -1.0}, {0.0, 0.0}}) == doctest::Approx(std::sqrt(2.0)));
        CHECK(two_norm(Matrix<2, 2>::identity() * 1e-9) == doctest::Approx(1e-9));
    }

    TEST_CASE_TEMPLATE("float QR factors columns of 1e-20 and 1e20 without overflow or flush", T, float, double) {
        // A plain sum of squares overflows past 1.8e19 in float and flushes below 1e-19.
        for (const T s : {T(1e-20), T(1e20)}) {
            const Matrix<3, 2, T> A = Matrix<3, 2, T>{{T(1), T(2)}, {T(3), T(4)}, {T(5), T(6)}} * s;
            const auto            q = qr_decompose(A);
            CHECK(no_column_dropped(q));
            const Matrix<3, 2, T> QR = q.Q * q.R;
            CHECK(Matrix<3, 2, T>(QR - A).norm() <= T(64) * eps<T>() * A.norm());
            const auto f = full_qr(A);
            CHECK(f.R(1, 0) == T(0));
            CHECK(f.R(2, 0) == T(0));
            CHECK(f.R(2, 1) == T(0));
            CHECK(Matrix<3, 2, T>((f.Q * f.R) - A).norm() <= T(64) * eps<T>() * A.norm());
        }
        const Matrix<3, 2, T> mixed{{T(1), T(0)}, {T(0), T(1e-20)}, {T(0), T(1e-20)}};
        const auto            qm = qr_decompose(mixed);
        CHECK(no_column_dropped(qm));
        CHECK(Matrix<3, 2, T>((qm.Q * qm.R) - mixed).norm() <= T(64) * eps<T>() * T(1e-20));
    }

    TEST_CASE("qr_decompose drops a column as roundoff relative to its own norm") {
        const Matrix<3, 2> tiny{{1e-14, 2e-14}, {0.0, 1e-14}, {0.0, 0.0}};
        CHECK(no_column_dropped(qr_decompose(tiny)));
        const Matrix<3, 2> dep{{1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0}};
        CHECK(qr_decompose(dep).R(1, 1) == 0.0);
    }
}

TEST_SUITE("qr_decompose orthogonality") {
    /// max |QᴴQ − I| over the thin Q's columns.
    template<typename T, size_t N, size_t M>
    double orthogonality_loss(const Matrix<N, M, T>& Q) {
        double worst = 0.0;
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < M; ++j) {
                T dot{};
                for (size_t r = 0; r < N; ++r) {
                    dot += damp::conj(Q(r, i)) * Q(r, j);
                }
                worst = damp::max(worst, static_cast<double>(damp::abs(dot - ((i == j) ? T{1} : T{0}))));
            }
        }
        return worst;
    }

    /// max |QR − A| / max |A|, and R upper-triangular with a real non-negative diagonal.
    template<typename T, size_t N, size_t M>
    void check_factors(const Matrix<N, M, T>& A, const QRDecomposition<T, N, M>& qr) {
        const Matrix<N, M, T> back = qr.Q * qr.R;
        double                err = 0.0;
        double                amax = 0.0;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                err = damp::max(err, static_cast<double>(damp::abs(back(i, j) - A(i, j))));
                amax = damp::max(amax, static_cast<double>(damp::abs(A(i, j))));
            }
        }
        CHECK(err <= 16.0 * damp::eps<double>() * amax);
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < i; ++j) {
                CHECK(damp::abs(qr.R(i, j)) == 0.0);
            }
            CHECK(damp::abs(qr.R(i, i) - T{damp::abs(qr.R(i, i))}) == 0.0); // real, ≥ 0
        }
    }

    // Läuchli's matrix [1 1 1; δ 0 0; 0 δ 0; 0 0 δ] has κ ≈ √3/δ. Modified Gram–Schmidt left
    // its Q columns 1.6e-9 from orthogonal at δ = 1e-7 (the κ·ε loss); Householder keeps them
    // within a few ε for any δ.
    TEST_CASE("an ill-conditioned Läuchli matrix gets an orthonormal Q") {
        constexpr double   delta = 1e-7;
        const Matrix<4, 3> L{{1.0, 1.0, 1.0}, {delta, 0.0, 0.0}, {0.0, delta, 0.0}, {0.0, 0.0, delta}};
        const auto         qr = qr_decompose(L);
        REQUIRE(no_column_dropped(qr));
        CHECK(orthogonality_loss(qr.Q) <= 16.0 * damp::eps<double>());
        check_factors(L, qr);
    }

    TEST_CASE("complex Läuchli: orthonormal Q, real non-negative diagonal") {
        using C = damp::complex<double>;
        constexpr double      delta = 1e-7;
        const Matrix<4, 3, C> L{
            {C(1.0, 0.0), C(0.0, 1.0), C(1.0, 1.0)},
            {C(0.0, delta), C(0.0, 0.0), C(0.0, 0.0)},
            {C(0.0, 0.0), C(delta, 0.0), C(0.0, 0.0)},
            {C(0.0, 0.0), C(0.0, 0.0), C(0.0, -delta)},
        };
        const auto qr = qr_decompose(L);
        REQUIRE(no_column_dropped(qr));
        CHECK(orthogonality_loss(qr.Q) <= 16.0 * damp::eps<double>());
        check_factors(L, qr);
    }

    // On rank loss R(j, j) is stored as zero and Q keeps an orthonormal column, so QR still
    // reconstructs A and QᴴQ = I.
    TEST_CASE("rank-deficient input keeps Q orthonormal") {
        const Matrix<3, 2> dep{{1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0}};
        const auto         qr = qr_decompose(dep);
        CHECK(qr.R(1, 1) == 0.0);
        CHECK(orthogonality_loss(qr.Q) <= 16.0 * damp::eps<double>());
        check_factors(dep, qr);
    }

    // More columns than rows: the trailing R diagonal is zero and its Q column is zero.
    TEST_CASE("a wide matrix factors with zero trailing columns of Q") {
        const Matrix<2, 3> W{{3.0, 1.0, 2.0}, {4.0, -1.0, 5.0}};
        const auto         qr = qr_decompose(W);
        CHECK(qr.R(0, 0) == doctest::Approx(5.0));
        CHECK(qr.R(2, 2) == 0.0);
        CHECK(qr.Q(0, 2) == 0.0);
        CHECK(qr.Q(1, 2) == 0.0);
        const Matrix<2, 3> back = qr.Q * qr.R;
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                CHECK(damp::abs(back(i, j) - W(i, j)) <= 16.0 * damp::eps<double>() * 5.0);
            }
        }
    }
}

TEST_SUITE("svd on rank-deficient input") {
    TEST_CASE("Jacobi SVD converges when a column collapses to roundoff") {
        // Controllability matrix of a system with an uncontrollable mode: the zero
        // middle row leaves a ~1e-159 residue column whose Gram sums underflow.
        const Matrix<3, 3> Co{{1.0, -1.0 / 3.0, 1.0 / 9.0}, {0.0, 0.0, 0.0}, {1.0, -1.0, 1.0}};
        const auto         s = svd(Co);
        CHECK(s.success);
        CHECK(rank(Co) == 2);
        CHECK(s.singular_values[2] < 1e-15 * s.singular_values[0]);
    }

    TEST_CASE("constexpr complex SVD matches runtime (wide [λI − A, B])") {
        using C = damp::complex<double>;
        constexpr Matrix<2, 3, C> M{{C(1.0e-3, 0.0), C(-1.0e-3, 0.0), C(5.0e-7, 0.0)}, {C(-1.0e-3, 0.0), C(1.0e-3, 0.0), C(1.0e-3, 0.0)}};
        constexpr auto            ce = svd(M);
        static_assert(ce.success);
        const auto rt = svd(M);
        CHECK(ce.singular_values[0] == doctest::Approx(rt.singular_values[0]).epsilon(1e-12));
        CHECK(ce.singular_values[1] == doctest::Approx(rt.singular_values[1]).epsilon(1e-12));
    }
}

TEST_SUITE("svd of a single column") {
    TEST_CASE("one column is its own left vector; sigma is its norm") {
        using C = damp::complex<double>;
        const Matrix<2, 1, C> a{{C(3.0, 0.0)}, {C(0.0, 4.0)}};
        const auto            s = svd(a);
        CHECK(s.success);
        CHECK(s.singular_values[0] == doctest::Approx(5.0));
        CHECK(rank(a) == 1);
        CHECK(rank(Matrix<2, 1, C>{}) == 0);
    }
}

TEST_SUITE("cholesky rank deficiency") {
    // The Gram matrix of N−1 vectors in Rᴺ is singular. Its last pivot is roundoff of either
    // sign; a positivity-only test accepted the positive half (about 50 % of draws) and factored
    // it into entries of 1/√eps, so mat::solve (Cholesky first for symmetric input) returned x of
    // size 10¹⁹ instead of nullopt. No floor refuses them all: the last pivot also carries the
    // forward error of the leading block, and about 3 % of draws clear 16·(i+1)·eps of the
    // pivot's terms. Rank decisions belong to mat::rank; this guards the positivity-only half.
    TEST_CASE_TEMPLATE("a singular Gram matrix is refused", T, float, double) {
        std::mt19937_64                   gen(7);
        std::uniform_real_distribution<T> U(T{-1}, T{1});
        constexpr int                     trials = 1000;
        int                               accepted = 0;
        for (int trial = 0; trial < trials; ++trial) {
            Matrix<5, 4, T> V{};
            for (size_t i = 0; i < 5; ++i) {
                for (size_t j = 0; j < 4; ++j) {
                    V(i, j) = U(gen);
                }
            }
            const Matrix<5, 5, T> G = V * V.t();
            if (cholesky(G).has_value()) {
                ++accepted;
            }
        }
        CHECK(accepted <= trials / 20);
    }

    TEST_CASE_TEMPLATE("a well-conditioned SPD matrix with graded diagonal still factors", T, float, double) {
        // diag(1, 1e-4, 1e-8) plus a small coupling: the pivots are small relative to each other
        // but not relative to their own diagonal entries.
        Matrix<3, 3, T> A{};
        A(0, 0) = T{1};
        A(1, 1) = static_cast<T>(1e-4);
        A(2, 2) = static_cast<T>(1e-8);
        A(0, 1) = A(1, 0) = static_cast<T>(1e-3);
        A(1, 2) = A(2, 1) = static_cast<T>(1e-7);
        REQUIRE(cholesky(A).has_value());
        const auto            L = *cholesky(A);
        const Matrix<3, 3, T> back = L * L.t();
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                CHECK(back(i, j) == doctest::Approx(A(i, j)).epsilon(static_cast<double>(T{16} * eps<T>())).scale(static_cast<double>(A(i, i))));
            }
        }
    }
}

namespace {

template<typename T>
constexpr double eps_of() {
    return static_cast<double>(eps<T>());
}

/// 1/(1 − α), α = (1 + √17)/8: the rook bound on |lᵢⱼ|.
constexpr double rook_L_bound = 1.0 / (1.0 - 0.6403882032022076);

/// x through a volatile, so the optimizer cannot fold a NaN test on a known constant.
template<typename T>
T opaque_value(T x) {
    volatile T v = x;
    return v;
}

/// D of an ldlt() factorization as a full matrix, in double.
template<size_t N, typename T>
Matrix<N, N, damp::complex<double>> ldlt_D(const LDLT<N, T>& f) {
    using C = damp::complex<double>;
    Matrix<N, N, C> D{};
    for (size_t i = 0; i < N; ++i) {
        D(i, i) = C{static_cast<double>(f.d[i]), 0.0};
        if (i + 1 < N) {
            const C e{static_cast<double>(damp::real(complex<scalar_type_t<T>>{f.e[i]})), static_cast<double>(damp::imag(complex<scalar_type_t<T>>{f.e[i]}))};
            D(i + 1, i) = e;
            D(i, i + 1) = damp::conj(e);
        }
    }
    return D;
}

/// L in complex double.
template<size_t N, typename T>
Matrix<N, N, damp::complex<double>> ldlt_L(const LDLT<N, T>& f) {
    using C = damp::complex<double>;
    Matrix<N, N, C> L{};
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            const complex<scalar_type_t<T>> l{f.L(i, j)};
            L(i, j) = C{static_cast<double>(l.real()), static_cast<double>(l.imag())};
        }
    }
    return L;
}

/// Backward error and growth of an ldlt() factorization, in double: ‖Pᵀ·L·D·Lᴴ·P − A‖ and ‖|L||D||Lᴴ|‖.
struct LdltError {
    double backward = 0.0; ///< ‖Pᵀ·L·D·Lᴴ·P − A‖_F
    double growth = 0.0;   ///< ‖|L||D||Lᴴ|‖_F
    double max_l = 0.0;    ///< max |lᵢⱼ|
    size_t blocks = 0;     ///< 2×2 blocks
};

template<size_t N, typename T>
LdltError ldlt_error(const Matrix<N, N, T>& A, const LDLT<N, T>& f) {
    const auto L = ldlt_L(f);
    const auto D = ldlt_D(f);
    const auto PAP = L * D * L.conjugate_transpose();
    LdltError  out{};
    double     back2 = 0.0;
    double     grow2 = 0.0;
    for (size_t i = 0; i < N; ++i) {
        out.blocks += (i + 1 < N && f.e[i] != T{0}) ? 1U : 0U;
        for (size_t j = 0; j < N; ++j) {
            const complex<scalar_type_t<T>> a{A(f.perm[i], f.perm[j])};
            const damp::complex<double>     r = PAP(i, j) - damp::complex<double>{static_cast<double>(a.real()), static_cast<double>(a.imag())};
            back2 += damp::real(r * damp::conj(r));
            double g = 0.0;
            for (size_t p = 0; p < N; ++p) {
                for (size_t q = 0; q < N; ++q) {
                    g += damp::abs(L(i, p)) * damp::abs(D(p, q)) * damp::abs(L(j, q));
                }
            }
            grow2 += g * g;
            out.max_l = damp::max(out.max_l, damp::abs(L(i, j)));
        }
    }
    out.backward = std::sqrt(back2);
    out.growth = std::sqrt(grow2);
    return out;
}

/// Relative 2-norm distance of x from a double reference.
template<size_t N, typename T>
double rel_error(const Matrix<N, 1, T>& x, const Matrix<N, 1, double>& ref) {
    double e = 0.0;
    double n = 0.0;
    for (size_t i = 0; i < N; ++i) {
        const double d = static_cast<double>(x(i, 0)) - ref(i, 0);
        e += d * d;
        n += ref(i, 0) * ref(i, 0);
    }
    return std::sqrt(e / n);
}

/// [[δ, 1, 0.5, 0.2], [1, 1, 0.3, 0.4], [0.5, 0.3, −1, 0.1], [0.2, 0.4, 0.1, 2]]: κ₂ = 4.6 for small δ.
template<typename T>
Matrix<4, 4, T> graded_corner_4(T delta) {
    return Matrix<4, 4, T>{{delta, T(1), T(0.5), T(0.2)}, {T(1), T(1), T(0.3), T(0.4)}, {T(0.5), T(0.3), T(-1), T(0.1)}, {T(0.2), T(0.4), T(0.1), T(2)}};
}

} // namespace

TEST_SUITE("ldlt") {
    TEST_CASE("a KKT saddle-point matrix factors and solves where Cholesky cannot") {
        // [[H, Aᵀ], [A, 0]] with H = diag(2, 3), A = [1, 1]: indefinite, and no pivot is needed
        const Matrix<3, 3> K{{2.0, 0.0, 1.0}, {0.0, 3.0, 1.0}, {1.0, 1.0, 0.0}};
        CHECK_FALSE(cholesky(K).has_value());
        const auto f = ldlt(K);
        REQUIRE(f.has_value());
        for (size_t i = 0; i < 3; ++i) {
            CHECK(f->perm[i] == i);
            CHECK(f->e[i] == 0.0);
        }
        CHECK(f->d[0] == doctest::Approx(2.0));
        CHECK(f->d[1] == doctest::Approx(3.0));
        CHECK(f->d[2] == doctest::Approx(-(1.0 / 2.0 + 1.0 / 3.0))); // Schur complement −A H⁻¹ Aᵀ
        CHECK(ldlt_error(K, *f).backward <= 8.0 * 3.0 * eps<double>() * K.norm());
        const Matrix<3, 1> b{{1.0}, {2.0}, {0.5}};
        const auto         x = ldlt_solve(K, b);
        const auto         x_lu = lu_solve(K, b);
        REQUIRE(x.has_value());
        REQUIRE(x_lu.has_value());
        for (size_t i = 0; i < 3; ++i) {
            CHECK((*x)(i, 0) == doctest::Approx((*x_lu)(i, 0)).epsilon(1e-12));
        }
    }

    TEST_CASE_TEMPLATE("[[δ, 1], [1, 0]] (κ = 1) takes a 2×2 pivot and solves to working precision", T, float, double) {
        // Unpivoted, d₁ = −1/δ and l = 1/δ: the factors grew as 1/δ, and ldlt_solve returned x with
        // a relative error of 3.6·10⁻³ in float at δ = 10⁻⁵ and 2.5·10⁻⁷ in double at δ = 10⁻¹⁰.
        const T               delta = damp::is_same_v<T, float> ? T(1e-5) : T(1e-10);
        const Matrix<2, 2, T> A{{delta, T(1)}, {T(1), T(0)}};
        const auto            f = ldlt(A);
        REQUIRE(f.has_value());
        CHECK(f->e[0] == T(1)); // a single 2×2 block, D = A
        CHECK(f->d[0] == delta);
        CHECK(f->d[1] == T(0));
        CHECK(f->perm[0] == 0);
        CHECK(f->perm[1] == 1);
        const Matrix<2, 1, T> b{{T(0.3)}, {T(0.7)}};
        const auto            x = ldlt_solve(A, b);
        REQUIRE(x.has_value());
        // x = (b₂, b₁ − δb₂), in double from the T data
        const Matrix<2, 1> ref{{static_cast<double>(b(1, 0))}, {static_cast<double>(b(0, 0)) - (static_cast<double>(delta) * static_cast<double>(b(1, 0)))}};
        CHECK(rel_error(x.value(), ref) <= 2.0 * 2.0 * eps_of<T>());
    }

    TEST_CASE_TEMPLATE("a 4×4 with κ = 4.6 and a tiny a₀₀ solves as accurately as LU", T, float, double) {
        // Unpivoted, the first pivot δ grew the factors by 1/δ: ldlt_solve lost 6.4·10⁻³ in float
        // at δ = 10⁻⁵ and 2.1·10⁻⁷ in double at δ = 10⁻¹⁰, where lu_solve keeps 4·10⁻⁸ and 10⁻¹⁶.
        const T               delta = damp::is_same_v<T, float> ? T(1e-5) : T(1e-10);
        const auto            A = graded_corner_4(delta);
        const Matrix<4, 1, T> b{{T(0.3)}, {T(-0.7)}, {T(1.1)}, {T(0.2)}};
        const auto            f = ldlt(A);
        REQUIRE(f.has_value());
        const LdltError err = ldlt_error(A, *f);
        CHECK(err.max_l <= rook_L_bound);
        CHECK(err.growth <= 2.0 * Matrix<4, 4>(A).norm());
        CHECK(err.backward <= 8.0 * 4.0 * eps_of<T>() * err.growth);
        const auto x = ldlt_solve(A, b);
        const auto ref = lu_solve(Matrix<4, 4>(A), Matrix<4, 1>(b)); // double LU on the same data: κ·10⁻¹⁶
        REQUIRE(x.has_value());
        REQUIRE(ref.has_value());
        CHECK(rel_error(x.value(), ref.value()) <= 4.0 * 4.0 * 4.6 * eps_of<T>());
    }

    TEST_CASE_TEMPLATE("a 2×2 block with rows below: exact multipliers, then a 1×1 pivot", T, float, double) {
        // E = [[0, 1], [1, 0]] is the first block (no 1×1 pivot clears α); Lᵢ = [s₁ s₂]E⁻¹ = [s₂ s₁].
        const Matrix<4, 4, T> A{{T(0), T(1), T(0.5), T(0.25)}, {T(1), T(0), T(0.75), T(-0.5)}, {T(0.5), T(0.75), T(2), T(0.5)}, {T(0.25), T(-0.5), T(0.5), T(-1)}};
        const auto            f = ldlt(A);
        REQUIRE(f.has_value());
        CHECK(f->e[0] == T(1));
        CHECK(f->e[1] == T(0));
        CHECK(f->e[2] == T(0));
        for (size_t i = 0; i < 4; ++i) {
            CHECK(f->perm[i] == i);
        }
        CHECK(f->L(2, 0) == T(0.75));
        CHECK(f->L(2, 1) == T(0.5));
        CHECK(f->L(3, 0) == T(-0.5));
        CHECK(f->L(3, 1) == T(0.25));
        CHECK(f->d[2] == T(1.25)); // 2 − 2·0.5·0.75
        CHECK(ldlt_error(A, *f).backward <= 8.0 * 4.0 * eps_of<T>() * Matrix<4, 4>(A).norm());
        const Matrix<4, 1, T> b{{T(1)}, {T(-2)}, {T(0.5)}, {T(3)}};
        const auto            x = ldlt_solve(A, b);
        const auto            ref = lu_solve(Matrix<4, 4>(A), Matrix<4, 1>(b));
        REQUIRE(x.has_value());
        REQUIRE(ref.has_value());
        CHECK(rel_error(x.value(), ref.value()) <= 64.0 * eps_of<T>());
    }

    TEST_CASE("float: a 2×2 block at 10⁻²⁵ and 10²⁵ solves without forming |w|²") {
        // |w|² flushes to zero below 10⁻¹⁹ and overflows past 1.8·10¹⁹ in float; the block is solved
        // scaled by |w|, as LU (which never squares an entry) does.
        for (const float s : {1e-25f, 1e25f}) {
            const Matrix<2, 2, float> A{{s * 1e-3f, s}, {s, s * 2e-3f}};
            const Matrix<2, 1, float> b{{0.3f}, {0.7f}};
            const auto                x = ldlt_solve(A, b);
            const auto                ref = lu_solve(Matrix<2, 2>(A), Matrix<2, 1>(b));
            REQUIRE(x.has_value());
            REQUIRE(ref.has_value());
            CHECK(rel_error(x.value(), ref.value()) <= 8.0 * eps_of<float>());
        }
    }

    TEST_CASE_TEMPLATE("singular matrices are refused, semidefinite or indefinite", T, float, double) {
        const Matrix<2, 2, T> S2{{T(1), T(2)}, {T(2), T(4)}};
        const Matrix<3, 3, T> S3{{T(0), T(1), T(1)}, {T(1), T(0), T(0)}, {T(1), T(0), T(0)}}; // a 2×2 pivot, then a lost row
        const Matrix<3, 3, T> D3{{T(1), T(0), T(0)}, {T(0), T(-1), T(0)}, {T(0), T(0), T(0)}};
        // u·uᵀ − v·vᵀ: rank 2, indefinite, reached through pivots that do not cancel exactly
        const Matrix<4, 1, T> u{{T(1)}, {T(2)}, {T(0.5)}, {T(1)}};
        const Matrix<4, 1, T> v{{T(0.25)}, {T(1)}, {T(1)}, {T(2)}};
        const Matrix<4, 4, T> R4 = (u * u.t()) - (v * v.t());
        CHECK_FALSE(ldlt(S2).has_value());
        CHECK_FALSE(ldlt(S3).has_value());
        CHECK_FALSE(ldlt(D3).has_value());
        CHECK_FALSE(ldlt(R4).has_value());
        CHECK_FALSE(ldlt(Matrix<3, 3, T>::zeros()).has_value());
        CHECK_FALSE(ldlt_solve(S3, Matrix<3, 1, T>{{T(1)}, {T(1)}, {T(1)}}).has_value());
        CHECK_FALSE(ldlt_solve(R4, Matrix<4, 1, T>{{T(1)}, {T(1)}, {T(1)}, {T(1)}}).has_value());
    }

    TEST_CASE_TEMPLATE("NaN and ±Inf are refused in every position, and so is an overflowing Schur complement", T, float, double) {
        // Under -ffast-math GCC folds !(|d| > floor) to |d| ≤ floor, false for NaN: a NaN pivot
        // passed, and ldlt returned NaN factors as a success.
        for (const T bad : {std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::infinity(), -std::numeric_limits<T>::infinity()}) {
            for (size_t i = 0; i < 4; ++i) {
                for (size_t j = 0; j <= i; ++j) {
                    INFO("entry (", i, ", ", j, ")");
                    auto A = graded_corner_4(T(1e-3));
                    A(i, j) = opaque_value(bad);
                    A(j, i) = A(i, j);
                    CHECK_FALSE(ldlt(A).has_value());
                    CHECK_FALSE(ldlt_solve(A, Matrix<4, 1, T>{{T(1)}, {T(1)}, {T(1)}, {T(1)}}).has_value());
                }
            }
        }
        // A NaN in either part of a complex entry
        using C = complex<T>;
        for (const C bad : {C{std::numeric_limits<T>::quiet_NaN(), T(0)}, C{T(0), std::numeric_limits<T>::quiet_NaN()}}) {
            Matrix<2, 2, C> A{{C{T(2), T(0)}, C{T(0)}}, {C{T(0)}, C{T(3), T(0)}}};
            A(1, 0) = C{opaque_value(bad.real()), opaque_value(bad.imag())};
            A(0, 1) = damp::conj(A(1, 0));
            CHECK_FALSE(ldlt(A).has_value());
            A = Matrix<2, 2, C>{{C{T(2), T(0)}, C{T(0)}}, {C{T(0)}, C{T(3), T(0)}}};
            A(1, 1) = C{opaque_value(bad.real()), opaque_value(bad.imag())};
            CHECK_FALSE(ldlt(A).has_value());
        }
        // Finite entries whose Schur complement −0.8·max − 0.19²/0.125·max overflows to −Inf. (The
        // pivot stays below max/4: -ffast-math divides by a pivot through its reciprocal, which
        // is subnormal past max/4 and flushes to zero.)
        const T               big = opaque_value(std::numeric_limits<T>::max());
        const Matrix<2, 2, T> O{{T(0.125) * big, static_cast<T>(0.19) * big}, {static_cast<T>(0.19) * big, static_cast<T>(-0.8) * big}};
        CHECK_FALSE(ldlt(O).has_value());
    }

    TEST_CASE("constexpr: a permuted factorization, a 2×2 block, and a solve") {
        // [[0, 1, 0], [1, 2, 1], [0, 1, 4]]: a₀₀ = 0 fails α·ω₀, a₁₁ = 2 pivots first
        static constexpr Matrix<3, 3> A{{0.0, 1.0, 0.0}, {1.0, 2.0, 1.0}, {0.0, 1.0, 4.0}};
        static constexpr auto         f = ldlt(A);
        static_assert(f.has_value());
        static_assert(f->perm[0] == 1 && f->perm[1] == 0 && f->perm[2] == 2);
        static_assert(f->d[0] == 2.0 && f->d[1] == -0.5 && f->d[2] == 4.0);
        static_assert(f->L(1, 0) == 0.5 && f->L(2, 0) == 0.5 && f->L(2, 1) == 1.0);
        static_assert(f->e[0] == 0.0 && f->e[1] == 0.0);
        static constexpr auto x = ldlt_solve(A, Matrix<3, 1>{{2.0}, {8.0}, {14.0}});
        static_assert(x.has_value() && (*x)(0, 0) == 1.0 && (*x)(1, 0) == 2.0 && (*x)(2, 0) == 3.0);
        // [[4, 2, 2], [2, 1, 3], [2, 3, 1]]: a 1×1 pivot leaves [[0, 2], [2, 0]], a 2×2 block
        static constexpr auto g = ldlt(Matrix<3, 3>{{4.0, 2.0, 2.0}, {2.0, 1.0, 3.0}, {2.0, 3.0, 1.0}});
        static_assert(g.has_value());
        static_assert(g->d[0] == 4.0 && g->d[1] == 0.0 && g->d[2] == 0.0 && g->e[0] == 0.0 && g->e[1] == 2.0);
        static_assert(g->L(1, 0) == 0.5 && g->L(2, 0) == 0.5 && g->L(2, 1) == 0.0);
        static_assert(!ldlt(Matrix<2, 2>{{1.0, 2.0}, {2.0, 4.0}}).has_value());
        static_assert(ldlt(Matrix<1, 1>{{-3.0}})->d[0] == -3.0 && !ldlt(Matrix<1, 1>{{0.0}}).has_value());
    }

    TEST_CASE("Hermitian: complex off-diagonal, real pivots") {
        const Matrix<2, 2, complex<double>> A{
            {complex<double>{4.0, 0.0}, complex<double>{1.0, -2.0}},
            {complex<double>{1.0, 2.0}, complex<double>{3.0, 0.0}},
        };
        const auto f = ldlt(A);
        REQUIRE(f.has_value());
        CHECK(f->d[0] == doctest::Approx(4.0));
        CHECK(f->d[1] == doctest::Approx(3.0 - 5.0 / 4.0)); // 3 − |1+2i|²/4
        const Matrix<2, 1, complex<double>> b{{complex<double>{1.0, 0.0}}, {complex<double>{0.0, 1.0}}};
        const auto                          x = ldlt_solve(A, b);
        REQUIRE(x.has_value());
        const auto r = A * (*x);
        for (size_t i = 0; i < 2; ++i) {
            CHECK((r(i, 0) - b(i, 0)).abs() <= 1e-12);
        }
    }

    TEST_CASE_TEMPLATE("Hermitian: a complex 2×2 block with a row below", T, float, double) {
        using C = complex<T>;
        const Matrix<3, 3, C> A{
            {C{T(0), T(0)}, C{T(1), T(-2)}, C{T(0.5), T(0)}},
            {C{T(1), T(2)}, C{T(0.25), T(0)}, C{T(0), T(1)}},
            {C{T(0.5), T(0)}, C{T(0), T(-1)}, C{T(1), T(0)}},
        };
        const auto f = ldlt(A);
        REQUIRE(f.has_value());
        const LdltError err = ldlt_error(A, *f);
        CHECK(err.blocks == 1);
        CHECK(err.max_l <= rook_L_bound);
        CHECK(err.backward <= 8.0 * 3.0 * eps_of<T>() * err.growth);
        const Matrix<3, 1, C> b{{C{T(1), T(0)}}, {C{T(0), T(1)}}, {C{T(-1), T(0.5)}}};
        const auto            x = ldlt_solve(A, b);
        REQUIRE(x.has_value());
        const Matrix<3, 1, C> r = (A * (*x)) - b;
        for (size_t i = 0; i < 3; ++i) {
            CHECK(static_cast<double>(r(i, 0).abs()) <= 64.0 * eps_of<T>() * static_cast<double>(x->norm() + b.norm()));
        }
    }

    TEST_CASE_TEMPLATE("random indefinite matrices: |L| ≤ 1/(1 − α), backward stable, solved", T, float, double) {
        // Off-diagonal entries in [−1, 1], diagonal shrunk by up to 10⁶ so that 2×2 pivots are common.
        std::mt19937_64                        gen(11);
        std::uniform_real_distribution<T>      U(T(-1), T(1));
        std::uniform_real_distribution<double> E(-6.0, 0.0);
        constexpr size_t                       N = 7;
        size_t                                 blocks = 0;
        double                                 worst_back = 0.0;
        double                                 worst_res = 0.0;
        for (int trial = 0; trial < 300; ++trial) {
            Matrix<N, N, T> A{};
            for (size_t i = 0; i < N; ++i) {
                for (size_t j = 0; j < i; ++j) {
                    A(i, j) = A(j, i) = U(gen);
                }
                A(i, i) = U(gen) * static_cast<T>(std::pow(10.0, E(gen)));
            }
            const auto f = ldlt(A);
            REQUIRE(f.has_value());
            const LdltError err = ldlt_error(A, *f);
            blocks += err.blocks;
            CHECK(err.max_l <= rook_L_bound * (1.0 + (16.0 * eps_of<T>())));
            worst_back = damp::max(worst_back, err.backward / (static_cast<double>(N) * eps_of<T>() * err.growth));
            Matrix<N, 1, T> b{};
            for (size_t i = 0; i < N; ++i) {
                b(i, 0) = U(gen);
            }
            const auto x = ldlt_solve(A, b);
            REQUIRE(x.has_value());
            const Matrix<N, 1> r = (Matrix<N, N>(A) * Matrix<N, 1>(x.value())) - Matrix<N, 1>(b);
            worst_res = damp::max(worst_res, r.norm() / (static_cast<double>(N) * eps_of<T>() * ((err.growth * Matrix<N, 1>(x.value()).norm()) + Matrix<N, 1>(b).norm())));
        }
        INFO("2×2 blocks ", blocks, ", backward error ", worst_back, " and residual ", worst_res, " units of n·eps");
        CHECK(blocks >= 300);
        CHECK(worst_back <= 1.0); // 0.12 measured
        CHECK(worst_res <= 1.0);  // 0.05 measured
    }
}
