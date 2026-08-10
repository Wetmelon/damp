// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

/**
 * @file test_robustness.cpp
 * @brief Adversarial edge-case coverage for damp/matrix: singular and
 *        ill-conditioned inverses, broadcasting, view aliasing, and the
 *        identities tying solve/inverse/det/rank together. Complements the
 *        happy-path coverage in test_matrix*.cpp.
 */

namespace {

// Largest |A·A⁻¹ − I| element, for verifying inverses on conditioned matrices.
template<size_t N, typename T>
T inverse_residual(const Matrix<N, N, T>& A) {
    auto inv = A.inverse();
    REQUIRE(inv.has_value());
    Matrix<N, N, T> prod = A * inv.value();
    Matrix<N, N, T> I = Matrix<N, N, T>::identity();
    T               worst = T{0};
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            worst = damp::max(worst, std::abs(prod(i, j) - I(i, j)));
        }
    }
    return worst;
}

} // namespace

TEST_SUITE("robustness") {

    TEST_CASE("inverse of structurally singular matrices returns nullopt") {
        SUBCASE("zero matrix 3x3") {
            Matrix<3, 3, double> Z{};
            CHECK_FALSE(Z.inverse().has_value());
        }
        SUBCASE("zero matrix 5x5") {
            Matrix<5, 5, double> Z{};
            CHECK_FALSE(Z.inverse().has_value());
        }
        SUBCASE("duplicate rows") {
            Matrix<3, 3, double> A = {
                {1.0, 2.0, 3.0},
                {1.0, 2.0, 3.0}, // identical to row 0
                {4.0, 5.0, 7.0},
            };
            CHECK_FALSE(A.inverse().has_value());
        }
        SUBCASE("duplicate columns") {
            Matrix<3, 3, double> A = {
                {1.0, 1.0, 3.0},
                {4.0, 4.0, 6.0},
                {7.0, 7.0, 9.0},
            };
            CHECK_FALSE(A.inverse().has_value());
        }
        SUBCASE("linearly dependent row (r2 = r0 + r1)") {
            Matrix<3, 3, double> A = {
                {1.0, 2.0, 3.0},
                {4.0, 5.0, 6.0},
                {5.0, 7.0, 9.0}, // = row0 + row1
            };
            CHECK_FALSE(A.inverse().has_value());
        }
        SUBCASE("zero row inside larger system") {
            Matrix<4, 4, double> A = Matrix<4, 4, double>::identity();
            A(2, 2) = 0.0; // diagonal zero -> singular
            CHECK_FALSE(A.inverse().has_value());
        }
    }

    TEST_CASE("det and rank agree with singularity") {
        Matrix<3, 3, double> singular = {
            {1.0, 2.0, 3.0},
            {4.0, 5.0, 6.0},
            {7.0, 8.0, 9.0},
        };
        CHECK(mat::det(singular) == doctest::Approx(0.0).epsilon(1e-10));
        CHECK(mat::rank(singular) == 2);
        CHECK_FALSE(singular.inverse().has_value());

        Matrix<4, 4, double> full = {
            {2.0, 0.0, 1.0, 0.0},
            {0.0, 3.0, 0.0, 1.0},
            {1.0, 0.0, 2.0, 0.0},
            {0.0, 1.0, 0.0, 3.0},
        };
        CHECK(mat::rank(full) == 4);
        CHECK(std::abs(mat::det(full)) > 1.0);
        CHECK(full.inverse().has_value());
    }

    TEST_CASE("ill-conditioned but invertible matrices still round-trip") {
        // Scaled diagonal spanning ~12 orders of magnitude (condition ~1e12).
        SUBCASE("widely scaled diagonal") {
            Matrix<3, 3, double> A = {
                {1e6, 0.0, 0.0},
                {0.0, 1.0, 0.0},
                {0.0, 0.0, 1e-6},
            };
            CHECK(inverse_residual(A) < 1e-9);
        }
        SUBCASE("4x4 Hilbert matrix (notoriously ill-conditioned)") {
            Matrix<4, 4, double> H;
            for (size_t i = 0; i < 4; ++i) {
                for (size_t j = 0; j < 4; ++j) {
                    H(i, j) = 1.0 / static_cast<double>(i + j + 1);
                }
            }
            // cond(H_4) ~ 1.5e4; inverse should still be accurate.
            CHECK(inverse_residual(H) < 1e-6);
        }
        SUBCASE("near-singular 2x2 (small but nonzero pivot)") {
            Matrix<2, 2, double> A = {
                {1.0, 1.0},
                {1.0, 1.0 + 1e-9},
            };
            CHECK(inverse_residual(A) < 1e-3); // loose: condition ~ 1e9
        }
    }

    TEST_CASE("inverse round-trips across sizes (well-conditioned)") {
        // Diagonally dominant => well-conditioned at every size.
        auto make_dd = []<size_t N>(std::integral_constant<size_t, N>) {
            Matrix<N, N, double> A;
            for (size_t i = 0; i < N; ++i) {
                for (size_t j = 0; j < N; ++j) {
                    A(i, j) = (i == j) ? static_cast<double>(N + 2) : 0.5;
                }
            }
            CHECK(inverse_residual(A) < 1e-10);
        };
        make_dd(std::integral_constant<size_t, 2>{});
        make_dd(std::integral_constant<size_t, 3>{});
        make_dd(std::integral_constant<size_t, 5>{});
        make_dd(std::integral_constant<size_t, 8>{});
    }

    TEST_CASE("solve(A,B) is consistent with inverse(A)*B") {
        Matrix<3, 3, double> A = {
            {4.0, 1.0, 2.0},
            {1.0, 5.0, 0.0},
            {2.0, 0.0, 6.0},
        };
        Matrix<3, 2, double> B = {
            {1.0, 2.0},
            {3.0, 4.0},
            {5.0, 6.0},
        };
        auto X = mat::solve(A, B);
        REQUIRE(X.has_value());
        Matrix<3, 2, double> viaInv = A.inverse().value() * B;
        // And the residual A*X - B should vanish.
        Matrix<3, 2, double> resid = (A * X.value()) - B;
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                CHECK(X.value()(i, j) == doctest::Approx(viaInv(i, j)).epsilon(1e-10));
                CHECK(resid(i, j) == doctest::Approx(0.0).epsilon(1e-10));
            }
        }
    }

    TEST_CASE("broadcasting operators") {
        Matrix<3, 3, double> A = Matrix<3, 3, double>::constant(1.0);
        SUBCASE("row vector broadcasts down the rows") {
            Matrix<1, 3, double> r = {{10.0, 20.0, 30.0}};
            Matrix<3, 3, double> S = A + r;
            for (size_t i = 0; i < 3; ++i) {
                CHECK(S(i, 0) == doctest::Approx(11.0));
                CHECK(S(i, 1) == doctest::Approx(21.0));
                CHECK(S(i, 2) == doctest::Approx(31.0));
            }
        }
        SUBCASE("column vector broadcasts across the columns") {
            Matrix<3, 1, double> c = {{1.0}, {2.0}, {3.0}};
            Matrix<3, 3, double> S = A + c;
            for (size_t j = 0; j < 3; ++j) {
                CHECK(S(0, j) == doctest::Approx(2.0));
                CHECK(S(1, j) == doctest::Approx(3.0));
                CHECK(S(2, j) == doctest::Approx(4.0));
            }
        }
        SUBCASE("subtraction broadcasts too") {
            Matrix<1, 3, double> r = {{1.0, 2.0, 3.0}};
            Matrix<3, 3, double> S = A - r;
            CHECK(S(0, 0) == doctest::Approx(0.0));
            CHECK(S(1, 2) == doctest::Approx(-2.0));
        }
    }

    TEST_CASE("view aliasing and consistency") {
        Matrix<3, 3, double> A = {
            {1.0, 2.0, 3.0},
            {4.0, 5.0, 6.0},
            {7.0, 8.0, 9.0},
        };

        SUBCASE("transpose view matches concrete transpose") {
            auto                 tv = A.t();
            Matrix<3, 3, double> tc = A.transpose();
            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = 0; j < 3; ++j) {
                    CHECK(tv(i, j) == tc(i, j));
                    CHECK(tv(i, j) == A(j, i));
                }
            }
        }

        SUBCASE("diagonal view sums to trace") {
            auto   d = A.diagonal();
            double s = 0.0;
            for (size_t i = 0; i < 3; ++i) {
                s += d(i);
            }
            CHECK(s == doctest::Approx(A.trace()));
        }

        SUBCASE("block write aliases the parent storage") {
            Matrix<4, 4, double> M = Matrix<4, 4, double>::identity();
            auto                 blk = M.block<2, 2>(1, 1);
            blk = Matrix<2, 2, double>::constant(7.0);
            CHECK(M(1, 1) == doctest::Approx(7.0));
            CHECK(M(2, 2) == doctest::Approx(7.0));
            CHECK(M(1, 2) == doctest::Approx(7.0));
            CHECK(M(0, 0) == doctest::Approx(1.0)); // outside the block untouched
            CHECK(M(3, 3) == doctest::Approx(1.0));
        }

        SUBCASE("head/tail blocks read the right rows") {
            auto h = A.head<1>();
            CHECK(h(0, 0) == doctest::Approx(1.0));
            CHECK(h(0, 2) == doctest::Approx(3.0));
            auto t = A.tail<1>();
            CHECK(t(0, 0) == doctest::Approx(7.0));
            CHECK(t(0, 2) == doctest::Approx(9.0));
        }
    }

    TEST_CASE("complex conjugate-transpose and norms") {
        using Cplx = damp::complex<double>;
        Matrix<2, 2, Cplx> A = {
            {Cplx(1.0, 1.0), Cplx(2.0, -1.0)},
            {Cplx(0.0, 3.0), Cplx(-1.0, 2.0)},
        };
        auto Ah = A.conjugate_transpose();
        // (Aᴴ)_{ij} = conj(A_{ji})
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                CHECK(Ah(i, j).real() == doctest::Approx(A(j, i).real()));
                CHECK(Ah(i, j).imag() == doctest::Approx(-A(j, i).imag()));
            }
        }
        // Frobenius norm is real and equals sqrt(sum |a|^2).
        double expect = std::sqrt(2.0 + 5.0 + 9.0 + 5.0);
        CHECK(A.norm() == doctest::Approx(expect));
    }

    TEST_CASE("math contracts at compile time (singular guards)") {
        // The constexpr math layer is documented to degrade gracefully rather
        // than produce NaN where NaN is unavailable at compile time.
        static_assert(damp::sqrt(-4.0) == 0.0, "constexpr sqrt of negative -> 0");
        static_assert(damp::fmod(3.0, 0.0) == 0.0, "constexpr fmod by zero -> 0");
        static_assert(damp::pow(-2.0, 0.5) == 0.0, "constexpr pow of negative base -> 0");
        static_assert(damp::pow(5.0, 0.0) == 1.0, "x^0 == 1");
        static_assert(damp::log(0.0) == 0.0, "constexpr log(0) -> 0");
        static_assert(damp::asin(2.0) > 1.5 && damp::asin(2.0) < 1.6, "asin clamps to pi/2");
    }
}
