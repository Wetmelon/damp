// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>

#include "damp/matrix/matrix.hpp"
#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

TEST_SUITE("quadratic_form") {
    // quadratic_form(M, X) computes S = M X Mᵀ for symmetric X, evaluating only
    // the lower triangle and mirroring it. Because it accumulates in the same
    // order as the general product M*X*M.t(), the lower triangle is bit-identical
    // to that product; the upper triangle is an exact copy of the lower, so the
    // result is symmetric to the last bit.

    TEST_CASE("matches M*X*Mᵀ on the lower triangle (non-square M)") {
        // M is 3×2, X is 2×2 symmetric → S is 3×3.
        const Matrix<3, 2, double> M = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
        const Matrix<2, 2, double> X = {{2.0, 1.0}, {1.0, 3.0}};

        const auto S = quadratic_form(M, X);
        const auto ref = M * X * M.t();

        // Lower triangle + diagonal: identical accumulation order → bit-exact.
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j <= i; ++j) {
                CHECK(S(i, j) == ref(i, j));
            }
        }
    }

    TEST_CASE("result is exactly symmetric") {
        const Matrix<3, 3, double> M = {{1.0, -2.0, 0.5}, {0.0, 3.0, -1.0}, {4.0, 1.0, 2.0}};
        const Matrix<3, 3, double> X = {{5.0, 2.0, 1.0}, {2.0, 4.0, 3.0}, {1.0, 3.0, 6.0}};

        const auto S = quadratic_form(M, X);

        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                CHECK(S(i, j) == S(j, i)); // bit-exact, not Approx
            }
        }
    }

    TEST_CASE("identity congruence returns X") {
        const Matrix<3, 3, double> X = {{5.0, 2.0, 1.0}, {2.0, 4.0, 3.0}, {1.0, 3.0, 6.0}};
        const auto                 I = Matrix<3, 3, double>::identity();

        const auto S = quadratic_form(I, X);

        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                CHECK(S(i, j) == doctest::Approx(X(i, j)));
            }
        }
    }

    TEST_CASE("row-vector M yields the scalar form mᵀ X m") {
        // M is 1×3 → S is 1×1 equal to m X mᵀ with m = [2, -1, 3].
        const Matrix<1, 3, double> m = {{2.0, -1.0, 3.0}};
        const Matrix<3, 3, double> X = {{5.0, 2.0, 1.0}, {2.0, 4.0, 3.0}, {1.0, 3.0, 6.0}};

        const auto S = quadratic_form(m, X);

        // Hand expansion: sum_ij m_i X_ij m_j.
        double expected = 0.0;
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                expected += m(0, i) * X(i, j) * m(0, j);
            }
        }
        CHECK(S(0, 0) == doctest::Approx(expected));
    }

    TEST_CASE("zero-dimension inner (NW = 0) gives a zero result") {
        // Mirrors the Kalman default NW = 0: G is NX×0, Q is 0×0 → GQGᵀ = 0.
        const Matrix<3, 0, double> G{};
        const Matrix<0, 0, double> Q{};

        const auto S = quadratic_form(G, Q);

        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                CHECK(S(i, j) == 0.0);
            }
        }
    }

    TEST_CASE("works for float element type") {
        const Matrix<2, 2, float> M = {{1.0F, 2.0F}, {3.0F, 4.0F}};
        const Matrix<2, 2, float> X = {{2.0F, 1.0F}, {1.0F, 3.0F}};

        const auto S = quadratic_form(M, X);
        const auto ref = M * X * M.t();

        CHECK(S(0, 0) == ref(0, 0));
        CHECK(S(1, 0) == ref(1, 0));
        CHECK(S(1, 1) == ref(1, 1));
        CHECK(S(0, 1) == S(1, 0));
    }

    TEST_CASE("evaluates at compile time") {
        constexpr Matrix<2, 2, double> M = {{1.0, 2.0}, {0.0, 3.0}};
        constexpr Matrix<2, 2, double> X = {{4.0, 1.0}, {1.0, 5.0}};
        constexpr auto                 S = quadratic_form(M, X);

        // S(0,0) = [1 2] X [1 2]ᵀ = 4 + 2(1)(2)... compute: row0 M = [1,2]
        // (M X)_0 = [1*4+2*1, 1*1+2*5] = [6, 11]; S00 = 6*1 + 11*2 = 28.
        static_assert(S(0, 0) == 28.0);
        static_assert(S(0, 1) == S(1, 0));
    }
}
