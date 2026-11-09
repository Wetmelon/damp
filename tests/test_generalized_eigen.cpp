// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Generalized eigenvalues (QZ): hand-built pencils with known answers. The scipy comparisons are
// in test_oracle_spectral.cpp.

#include <cmath>
#include <limits>

#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/generalized_eigen.hpp"
#include "damp/matrix/matrix.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

/// Every expected root has a distinct computed root within tol·max(1, |root|).
template<size_t N, size_t M>
bool contains_roots(const damp::array<damp::complex<double>, N>& got, size_t n_got, const damp::array<damp::complex<double>, M>& want, double tol) {
    if (n_got != M) {
        return false;
    }
    damp::array<bool, N> used{};
    for (size_t i = 0; i < M; ++i) {
        bool found = false;
        for (size_t j = 0; j < n_got && !found; ++j) {
            if (!used[j] && damp::abs(got[j] - want[i]) <= tol * damp::max(1.0, damp::abs(want[i]))) {
                used[j] = true;
                found = true;
            }
        }
        if (!found) {
            return false;
        }
    }
    return true;
}

} // namespace

TEST_SUITE("Generalized eigenvalues") {
    TEST_CASE("B = I gives the eigenvalues of A") {
        const Matrix<4, 4> A{{1.0, 2.0, 0.0, -1.0}, {0.5, -3.0, 1.0, 0.0}, {1.0, 0.0, 2.0, 1.5}, {0.0, 1.0, -1.0, 0.5}};
        const auto         ge = mat::generalized_eigenvalues(A, Matrix<4, 4>::identity());
        const auto         ev = mat::compute_eigenvalues(A);
        REQUIRE(ge.success);
        REQUIRE(ge.n_finite == 4);
        CHECK(contains_roots(ge.values, ge.n_finite, ev.values, 1e-12));
    }

    TEST_CASE("a singular B gives an infinite eigenvalue, constant-evaluated") {
        // det(A − λ·diag(1, 1, 0)) = 2(λ − 2)²
        constexpr Matrix<3, 3> A{{1.0, 2.0, 0.0}, {0.0, 3.0, 1.0}, {1.0, 0.0, 2.0}};
        constexpr Matrix<3, 3> B{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}};
        constexpr auto         ev = mat::generalized_eigenvalues(A, B);
        static_assert(ev.success && ev.n_finite == 2);
        CHECK(ev.values[0].real() == doctest::Approx(2.0).epsilon(1e-7)); // a double root: √ε
        CHECK(ev.values[1].real() == doctest::Approx(2.0).epsilon(1e-7));
    }

    TEST_CASE("a zero leading diagonal entry of B is chased down and deflated") {
        // det(A − λ·diag(0, 1, 1)) = λ² − 5λ + 8: 5/2 ± i√7/2
        const Matrix<3, 3> A{{1.0, 2.0, 0.0}, {0.0, 3.0, 1.0}, {1.0, 0.0, 2.0}};
        const Matrix<3, 3> B{{0.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
        const auto         ev = mat::generalized_eigenvalues(A, B);
        REQUIRE(ev.success);
        REQUIRE(ev.n_finite == 2);
        CHECK(ev.values[0].real() == doctest::Approx(2.5));
        CHECK(damp::abs(ev.values[0].imag()) == doctest::Approx(std::sqrt(7.0) / 2.0));
        CHECK(ev.values[1].imag() == doctest::Approx(-ev.values[0].imag())); // an exact conjugate pair
    }

    TEST_CASE("B = 0 has every eigenvalue at infinity") {
        const Matrix<3, 3> A{{1.0, 2.0, 0.0}, {0.0, 3.0, 1.0}, {1.0, 0.0, 2.0}};
        const auto         ev = mat::generalized_eigenvalues(A, Matrix<3, 3>{});
        CHECK(ev.success);
        CHECK(ev.n_finite == 0);
    }

    TEST_CASE("a non-finite entry is refused") {
        Matrix<2, 2> A = Matrix<2, 2>::identity();
        A(0, 1) = std::numeric_limits<double>::quiet_NaN();
        CHECK_FALSE(mat::generalized_eigenvalues(A, Matrix<2, 2>::identity()).success);
    }

    TEST_CASE("float") {
        const Matrix<3, 3, float> A{{1.0f, 2.0f, 0.0f}, {0.0f, 3.0f, 1.0f}, {1.0f, 0.0f, 2.0f}};
        const Matrix<3, 3, float> B{{0.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, {0.0f, 0.0f, 1.0f}};
        const auto                ev = mat::generalized_eigenvalues(A, B);
        REQUIRE(ev.success);
        REQUIRE(ev.n_finite == 2);
        CHECK(ev.values[0].real() == doctest::Approx(2.5f).epsilon(1e-5));
        const auto evd = ev.as<double>();
        CHECK(evd.n_finite == 2);
    }
}
