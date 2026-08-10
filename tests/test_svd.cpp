// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>

#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/svd.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using namespace damp::mat;
using damp::complex;

namespace {

// Max-abs reconstruction error ‖A − U·Σ·Vᴴ‖_∞.
template<size_t M, size_t N, typename T>
double reconstruction_error(const Matrix<M, N, T>& A, const SVDResult<M, N, T>& s) {
    constexpr size_t K = (M < N) ? M : N;
    Matrix<M, N, T>  Sigma = Matrix<M, N, T>::zeros();
    for (size_t k = 0; k < K; ++k) {
        Sigma(k, k) = T{s.singular_values[k]};
    }
    const Matrix<M, N, T> R = s.singular_U * Sigma * s.singular_V.conjugate_transpose();
    double                e = 0.0;
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            e = std::max(e, static_cast<double>(damp::abs(R(i, j) - A(i, j))));
        }
    }
    return e;
}

// Max-abs deviation of QᴴQ from identity.
template<size_t N, typename T>
double unitary_error(const Matrix<N, N, T>& Q) {
    const Matrix<N, N, T> G = Q.conjugate_transpose() * Q;
    double                e = 0.0;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            const T t = G(i, j) - (i == j ? T{1} : T{0});
            e = std::max(e, static_cast<double>(damp::abs(t)));
        }
    }
    return e;
}

} // namespace

TEST_SUITE("svd") {
    TEST_CASE("Tall real matrix: known singular values and factor properties") {
        // [[1,2],[3,4],[5,6]] has σ = {9.5255180..., 0.5143006...} (MATLAB® svd).
        const Matrix<3, 2, double> A = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
        const auto                 s = svd(A);

        CHECK(s.converged);
        CHECK(s.singular_values[0] == doctest::Approx(9.52551809));
        CHECK(s.singular_values[1] == doctest::Approx(0.51430058));
        CHECK(s.singular_values[0] >= s.singular_values[1]); // descending
        CHECK(reconstruction_error(A, s) < 1e-12);
        CHECK(unitary_error(s.singular_U) < 1e-12);
        CHECK(unitary_error(s.singular_V) < 1e-12);
    }

    TEST_CASE("Wide real matrix decomposes via the transpose path") {
        const Matrix<2, 3, double> A = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
        const auto                 s = svd(A);

        CHECK(s.converged);
        CHECK(reconstruction_error(A, s) < 1e-12);
        CHECK(unitary_error(s.singular_U) < 1e-12); // 2×2
        CHECK(unitary_error(s.singular_V) < 1e-12); // 3×3
    }

    TEST_CASE("Diagonal matrix: singular values are |diagonal|, sorted") {
        const Matrix<3, 3, double> A = {{3.0, 0.0, 0.0}, {0.0, -5.0, 0.0}, {0.0, 0.0, 2.0}};
        const auto                 s = svd(A);

        CHECK(s.singular_values[0] == doctest::Approx(5.0));
        CHECK(s.singular_values[1] == doctest::Approx(3.0));
        CHECK(s.singular_values[2] == doctest::Approx(2.0));
        CHECK(reconstruction_error(A, s) < 1e-12);
    }

    TEST_CASE("Rank-deficient matrix: zero singular value and correct rank") {
        // Row 2 = 2·row 1 ⇒ rank 2.
        const Matrix<3, 3, double> A = {{1.0, 2.0, 3.0}, {2.0, 4.0, 6.0}, {1.0, 0.0, 1.0}};
        const auto                 s = svd(A);

        CHECK(reconstruction_error(A, s) < 1e-12);
        CHECK(unitary_error(s.singular_U) < 1e-12);
        CHECK(unitary_error(s.singular_V) < 1e-12);
        CHECK(s.singular_values[2] < 1e-10); // numerically zero
        CHECK(rank_from_svd(s) == 2);
    }

    TEST_CASE("Complex matrix reconstructs with unitary factors") {
        const Matrix<3, 2, complex<double>> A = {
            {{1.0, 1.0}, {2.0, -1.0}},
            {{0.0, 2.0}, {1.0, 0.0}},
            {{3.0, 0.0}, {-1.0, 1.0}},
        };
        const auto s = svd(A);

        CHECK(s.converged);
        CHECK(reconstruction_error(A, s) < 1e-12);
        CHECK(unitary_error(s.singular_U) < 1e-12);
        CHECK(unitary_error(s.singular_V) < 1e-12);
    }

    TEST_CASE("Zero matrix: all singular values zero, factors identity-like") {
        const Matrix<2, 2, double> A = Matrix<2, 2, double>::zeros();
        const auto                 s = svd(A);
        CHECK(s.singular_values[0] == doctest::Approx(0.0));
        CHECK(s.singular_values[1] == doctest::Approx(0.0));
        CHECK(rank_from_svd(s) == 0);
        CHECK(unitary_error(s.singular_U) < 1e-12);
    }

    TEST_CASE("constexpr evaluation") {
        constexpr Matrix<2, 2, double> A = {{2.0, 0.0}, {0.0, 3.0}};
        constexpr auto                 s = svd(A);
        static_assert(s.converged);
        // Largest singular value of diag(2,3) is 3.
        static_assert(s.singular_values[0] > 2.9 && s.singular_values[0] < 3.1);
    }
}

TEST_SUITE("pseudo_inverse") {
    TEST_CASE("Overdetermined: Moore–Penrose conditions hold") {
        const Matrix<4, 2, double> A = {{1.0, 0.0}, {1.0, 1.0}, {1.0, 2.0}, {1.0, 3.0}};
        const auto                 P = pseudo_inverse(A);

        // A·A⁺·A = A and A⁺·A·A⁺ = A⁺.
        const Matrix<4, 2, double> APA = A * P * A;
        const Matrix<2, 4, double> PAP = P * A * P;
        double                     e1 = 0.0;
        double                     e2 = 0.0;
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                e1 = std::max(e1, std::abs(APA(i, j) - A(i, j)));
                e2 = std::max(e2, std::abs(PAP(j, i) - P(j, i)));
            }
        }
        CHECK(e1 < 1e-12);
        CHECK(e2 < 1e-12);

        // A⁺·A is symmetric (real Hermitian) for a full-column-rank A it equals I₂.
        const Matrix<2, 2, double> PA = P * A;
        CHECK(PA(0, 0) == doctest::Approx(1.0));
        CHECK(PA(1, 1) == doctest::Approx(1.0));
        CHECK(std::abs(PA(0, 1)) < 1e-12);
    }

    TEST_CASE("Least-squares solution matches normal equations") {
        // Fit y = a + b·x to (0,1),(1,3),(2,2),(3,5): A⁺·b is the LS fit.
        const Matrix<4, 2, double> A = {{1.0, 0.0}, {1.0, 1.0}, {1.0, 2.0}, {1.0, 3.0}};
        const Matrix<4, 1, double> b = {{1.0}, {3.0}, {2.0}, {5.0}};
        const auto                 P = pseudo_inverse(A);
        const Matrix<2, 1, double> x = P * b;

        // Closed-form LS for this data (normal equations): intercept = slope = 1.1.
        CHECK(x(0, 0) == doctest::Approx(1.1));
        CHECK(x(1, 0) == doctest::Approx(1.1));
    }

    TEST_CASE("Square nonsingular: pseudoinverse equals inverse") {
        const Matrix<3, 3, double> A = {{2.0, 0.0, 1.0}, {1.0, 3.0, 2.0}, {1.0, 0.0, 4.0}};
        const auto                 P = pseudo_inverse(A);
        const auto                 Inv = A.inverse().value();
        double                     e = 0.0;
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                e = std::max(e, std::abs(P(i, j) - Inv(i, j)));
            }
        }
        CHECK(e < 1e-12);
    }

    TEST_CASE("Complex pseudoinverse satisfies A·A⁺·A = A") {
        const Matrix<3, 2, complex<double>> A = {
            {{1.0, 1.0}, {2.0, -1.0}},
            {{0.0, 2.0}, {1.0, 0.0}},
            {{3.0, 0.0}, {-1.0, 1.0}},
        };
        const auto                          P = pseudo_inverse(A);
        const Matrix<3, 2, complex<double>> APA = A * P * A;
        double                              e = 0.0;
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                e = std::max(e, damp::abs(APA(i, j) - A(i, j)));
            }
        }
        CHECK(e < 1e-12);
    }
}

TEST_SUITE("null_space") {
    TEST_CASE("Wide full-rank matrix: 1-D kernel orthogonal to the rows") {
        const Matrix<2, 3, double> A = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
        const auto                 ns = null_space(A);
        REQUIRE(ns.dim == 1);

        // Kernel basis is the trailing column of ns.vectors.
        Matrix<3, 1, double> n;
        for (size_t i = 0; i < 3; ++i) {
            n(i, 0) = ns.vectors(i, 3 - ns.dim);
        }
        const Matrix<2, 1, double> An = A * n;
        CHECK(An.norm() < 1e-12);
        CHECK(n.norm() == doctest::Approx(1.0)); // orthonormal
    }

    TEST_CASE("Full-rank square matrix has trivial kernel") {
        const Matrix<3, 3, double> A = {{2.0, 0.0, 1.0}, {1.0, 3.0, 2.0}, {1.0, 0.0, 4.0}};
        const auto                 ns = null_space(A);
        CHECK(ns.dim == 0);
    }

    TEST_CASE("Rank-deficient matrix: kernel dimension matches the defect") {
        // Column 3 = column 1 + column 2 ⇒ 1-D kernel [1, 1, −1].
        const Matrix<3, 3, double> A = {{1.0, 0.0, 1.0}, {0.0, 1.0, 1.0}, {2.0, 1.0, 3.0}};
        const auto                 ns = null_space(A);
        REQUIRE(ns.dim == 1);
        Matrix<3, 1, double> n;
        for (size_t i = 0; i < 3; ++i) {
            n(i, 0) = ns.vectors(i, 2);
        }
        const Matrix<3, 1, double> An = A * n;
        CHECK(An.norm() < 1e-10);
    }
}
