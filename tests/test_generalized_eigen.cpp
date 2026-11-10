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

TEST_SUITE("Generalized eigenvalues: close pairs and scaling") {
    // The 2×2 block's eigenvalues came from the unshifted quadratic det(A − λB) = 0, whose
    // discriminant q₁² − 4q₂q₀ cancels when the two eigenvalues are close: their imaginary part
    // lost ε·λ²/gap, and whether they were real followed the sign of the cancelled difference.
    // A rotation by θ has eigenvalues c ± i·s exactly (c, s the stored cos θ, sin θ); at θ = 1e-4
    // in float and 1e-8 in double the old code returned a real double eigenvalue.
    TEST_CASE_TEMPLATE("a rotation pencil keeps its small imaginary part", T, float, double) {
        const T               theta = damp::is_same_v<T, float> ? static_cast<T>(1e-4) : static_cast<T>(1e-8);
        const T               c = damp::cos(theta);
        const T               s = damp::sin(theta);
        const Matrix<2, 2, T> A{{c, -s}, {s, c}};
        const T               e = eps<T>();
        const auto            ev = mat::generalized_eigenvalues(A, Matrix<2, 2, T>::identity());
        REQUIRE(ev.success);
        REQUIRE(ev.n_finite == 2);
        CHECK(damp::abs(ev.values[0].real() - c) <= T{4} * e);
        CHECK(damp::abs(ev.values[0].imag() - s) <= T{8} * e * s); // to the digits of s, not of 1
        CHECK(ev.values[1].real() == ev.values[0].real());
        CHECK(ev.values[1].imag() == -ev.values[0].imag());

        // The same eigenvalues from (A·Bt, Bt): the rounding of A·Bt moves them by a few ε.
        const Matrix<2, 2, T> Bt{{T{2}, T{0.5}}, {T{0}, T{0.25}}};
        const auto            evt = mat::generalized_eigenvalues(Matrix<2, 2, T>(A * Bt), Bt);
        REQUIRE(evt.success);
        REQUIRE(evt.n_finite == 2);
        const damp::complex<T> want{c, s};
        const T                err0 = damp::min(damp::abs(evt.values[0] - want), damp::abs(evt.values[1] - want));
        CHECK(err0 <= T{32} * e);
        CHECK(damp::abs(evt.values[0].imag()) > s / T{2});
    }

    // A symmetric A with B = I has real eigenvalues a + g/2 ± √(g²/4 + d²). With d and g small
    // the old discriminant cancelled to a negative number and reported a complex pair (216 of
    // 40,000 such double pencils, more under x86-64-v3); the shifted one is (g/2)² + d².
    TEST_CASE_TEMPLATE("a symmetric pencil with close eigenvalues stays real", T, float, double) {
        const bool            f = damp::is_same_v<T, float>;
        const T               a = T{1} + T{0.0009765625};                                                          // 1 + 2⁻¹⁰, exact
        const T               d = f ? static_cast<T>(3.814697265625e-06) : static_cast<T>(1.4551915228366852e-11); // 2⁻¹⁸, 2⁻³⁶
        const T               g = f ? static_cast<T>(6.103515625e-05) : static_cast<T>(7.450580596923828125e-09);  // 2⁻¹⁴, 2⁻²⁷
        const Matrix<2, 2, T> A{{a, d}, {d, a + g}};
        const auto            ev = mat::generalized_eigenvalues(A, Matrix<2, 2, T>::identity());
        REQUIRE(ev.success);
        REQUIRE(ev.n_finite == 2);
        const double mid = static_cast<double>(a) + (static_cast<double>(g) / 2.0);
        const double half_gap = std::hypot(static_cast<double>(g) / 2.0, static_cast<double>(d));
        const double want_hi = mid + half_gap;
        const double want_lo = mid - half_gap;
        const double tol = 4.0 * static_cast<double>(eps<T>());
        for (size_t i = 0; i < 2; ++i) {
            CHECK(ev.values[i].imag() == T{0});
        }
        const double hi = static_cast<double>(damp::max(ev.values[0].real(), ev.values[1].real()));
        const double lo = static_cast<double>(damp::min(ev.values[0].real(), ev.values[1].real()));
        CHECK(damp::abs(hi - want_hi) <= tol);
        CHECK(damp::abs(lo - want_lo) <= tol);
    }
}
