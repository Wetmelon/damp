// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstdlib>
#include <numbers>

#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using namespace damp::mat;

TEST_CASE("Complex Matrix Support") {
    using Cplx = damp::complex<double>;

    SUBCASE("Complex matrix creation") {
        Matrix<2, 2, Cplx> A;
        A(0, 0) = Cplx{1.0, 2.0};
        A(0, 1) = Cplx{3.0, -1.0};
        A(1, 0) = Cplx{-2.0, 1.0};
        A(1, 1) = Cplx{0.0, 1.0};

        CHECK(A(0, 0).real() == 1.0);
        CHECK(A(0, 0).imag() == 2.0);
    }
}

TEST_CASE("Eigen Decomposition") {
    SUBCASE("2x2 symmetric matrix") {
        Mat2 A = {{4.0, 1.0}, {1.0, 3.0}};
        auto result = compute_eigenvalues(A);

        CHECK(result.converged);

        auto lambda1 = result.values[0];
        auto lambda2 = result.values[1];

        CHECK(std::abs(lambda1.imag()) < 1e-6);
        CHECK(std::abs(lambda2.imag()) < 1e-6);

        double ev1 = lambda1.real();
        double ev2 = lambda2.real();

        if (ev1 < ev2) {
            damp::swap(ev1, ev2);
        }

        CHECK(ev1 == doctest::Approx(4.618).epsilon(0.01));
        CHECK(ev2 == doctest::Approx(2.382).epsilon(0.01));
    }

    SUBCASE("3x3 symmetric: Schur vectors satisfy A v ≈ λ v") {
        // Non-diagonal SPD-ish matrix so the Schur basis is a real eigenbasis.
        Mat3 A = {
            {4.0, 1.0, 0.5},
            {1.0, 3.0, 0.2},
            {0.5, 0.2, 2.0},
        };
        auto result = compute_eigenvalues(A);
        REQUIRE(result.converged);

        for (size_t i = 0; i < 3; ++i) {
            ColVec<3, double> v{
                result.vectors(0, i).real(),
                result.vectors(1, i).real(),
                result.vectors(2, i).real(),
            };
            const double nrm = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            REQUIRE(nrm > 1e-12);
            v = v * (1.0 / nrm);
            const auto   Av = A * v;
            const double lam = result.values[i].real();
            CHECK(Av[0] == doctest::Approx(lam * v[0]).epsilon(1e-5));
            CHECK(Av[1] == doctest::Approx(lam * v[1]).epsilon(1e-5));
            CHECK(Av[2] == doctest::Approx(lam * v[2]).epsilon(1e-5));
        }
    }

    SUBCASE("Diagonal matrix") {
        Mat3<double> D = {{2.0, 0.0, 0.0}, {0.0, 5.0, 0.0}, {0.0, 0.0, -3.0}};

        auto result = compute_eigenvalues(D);
        CHECK(result.converged);

        std::array<double, 3> evs{};
        for (size_t i = 0; i < 3; ++i) {
            evs[i] = result.values[i].real();
            CHECK(std::abs(result.values[i].imag()) < 1e-6);
        }

        std::ranges::sort(evs);

        CHECK(evs[0] == doctest::Approx(-3.0));
        CHECK(evs[1] == doctest::Approx(2.0));
        CHECK(evs[2] == doctest::Approx(5.0));
    }

    SUBCASE("2x2 general matrix (Buck converter closed-loop)") {
        // This is approximately what the Buck converter A-B*K looks like
        // with large values typical of power electronics
        Mat2 A_cl = {{-3e4, -5e3}, {6.8e4, -6.8e3}};

        auto result = compute_eigenvalues(A_cl);
        CHECK(result.converged);

        // Both eigenvalues should be computed (non-zero)
        double ev1 = result.values[0].real();
        double ev2 = result.values[1].real();

        CHECK(ev1 != 0.0);
        CHECK(ev2 != 0.0);
        // trace = -3e4 - 6.8e3 = -36800, so sum of eigenvalues should be this
        CHECK(ev1 + ev2 == doctest::Approx(-36800.0).epsilon(0.01));
    }

    SUBCASE("2x2 consteval eigenvalue computation") {
        // Test that eigenvalues work in consteval context
        constexpr Mat2 A_cl = {{-3e4, -5e3}, {6.8e4, -6.8e3}};
        constexpr auto result = compute_eigenvalues(A_cl);

        static_assert(result.converged, "Eigenvalue computation should converge");

        constexpr double ev1 = result.values[0].real();
        constexpr double ev2 = result.values[1].real();

        // At compile time, verify eigenvalues are non-zero
        static_assert(ev1 != 0.0 || ev2 != 0.0, "At least one eigenvalue should be non-zero");

        CHECK(ev1 != 0.0);
        CHECK(ev2 != 0.0);
        CHECK(ev1 + ev2 == doctest::Approx(-36800.0).epsilon(0.01));
    }
}

TEST_CASE("compute_eigenvalues - 1x1 matrix") {
    Matrix<1, 1> A = {{5.0}};

    auto result = compute_eigenvalues(A);

    CHECK(result.converged);
    CHECK(result.values[0].real() == doctest::Approx(5.0));
}

TEST_CASE("compute_eigenvalues - 2x2 diagonal matrix") {
    Matrix<2, 2> A = {
        {3.0, 0.0},
        {0.0, 7.0}
    };

    auto result = compute_eigenvalues(A);

    CHECK(result.converged);
    // Eigenvalues should be 3 and 7 (on diagonal)
    CHECK((result.values[0].real() == doctest::Approx(3.0) || result.values[0].real() == doctest::Approx(7.0)));
    CHECK((result.values[1].real() == doctest::Approx(7.0) || result.values[1].real() == doctest::Approx(3.0)));
}

TEST_CASE("compute_eigenvalues - 2x2 symmetric matrix") {
    // Symmetric matrix with known eigenvalues
    Matrix<2, 2> A = {
        {4.0, 1.0},
        {1.0, 4.0}
    };

    auto result = compute_eigenvalues(A);

    CHECK(result.converged);
    // Eigenvalues should be 3 and 5
    double lambda1 = result.values[0].real();
    double lambda2 = result.values[1].real();

    // Check that we got both eigenvalues (order may vary)
    bool has_3 = (std::abs(lambda1 - 3.0) < 1e-6 || std::abs(lambda2 - 3.0) < 1e-6);
    bool has_5 = (std::abs(lambda1 - 5.0) < 1e-6 || std::abs(lambda2 - 5.0) < 1e-6);

    CHECK(has_3);
    CHECK(has_5);
}

TEST_CASE("compute_eigenvalues - 2x2 identity matrix") {
    Matrix<2, 2> A = Matrix<2, 2>::identity();

    auto result = compute_eigenvalues(A);

    CHECK(result.converged);
    CHECK(result.values[0].real() == doctest::Approx(1.0));
    CHECK(result.values[1].real() == doctest::Approx(1.0));
}

TEST_CASE("compute_eigenvalues - 3x3 symmetric matrix") {
    // Symmetric matrix
    Matrix<3, 3> A = {
        {6.0, 2.0, 1.0},
        {2.0, 3.0, 1.0},
        {1.0, 1.0, 1.0}
    };

    auto result = compute_eigenvalues(A);

    CHECK(result.converged);

    // Sum of eigenvalues should equal trace
    double trace = A(0, 0) + A(1, 1) + A(2, 2);
    double sum_eigenvalues = result.values[0].real() + result.values[1].real() + result.values[2].real();
    CHECK(sum_eigenvalues == doctest::Approx(trace).epsilon(1e-5));
}

TEST_CASE("compute_eigenvalues - positive definite matrix") {
    // Positive definite matrix (all eigenvalues should be positive)
    Matrix<2, 2> A = {
        {2.0, 0.5},
        {0.5, 2.0}
    };

    auto result = compute_eigenvalues(A);

    CHECK(result.converged);
    CHECK(result.values[0].real() > 0.0);
    CHECK(result.values[1].real() > 0.0);
}

TEST_CASE("compute_eigenvalues - negative eigenvalues") {
    // Matrix with negative eigenvalues
    Matrix<2, 2> A = {
        {-3.0, 0.0},
        {0.0, -5.0}
    };

    auto result = compute_eigenvalues(A);

    CHECK(result.converged);
    CHECK(result.values[0].real() == doctest::Approx(-3.0));
    CHECK(result.values[1].real() == doctest::Approx(-5.0));
}

TEST_CASE("compute_eigenvalues - check for positive definiteness") {
    SUBCASE("positive definite matrix") {
        Matrix<3, 3> Q = {
            {2.0, 0.1, 0.1},
            {0.1, 2.0, 0.1},
            {0.1, 0.1, 2.0}
        };

        auto result = compute_eigenvalues(Q);

        if (result.converged) {
            bool all_positive = true;
            for (size_t i = 0; i < 3; ++i) {
                if (result.values[i].real() <= 0.0) {
                    all_positive = false;
                    break;
                }
            }
            CHECK(all_positive);
        }
    }

    SUBCASE("non-positive definite matrix") {
        Matrix<2, 2> Q = {
            {1.0, 2.0},
            {2.0, 1.0}
        };

        auto result = compute_eigenvalues(Q);

        if (result.converged) {
            // Should have eigenvalues 3 and -1
            double lambda1 = result.values[0].real();
            double lambda2 = result.values[1].real();

            bool has_negative = (lambda1 < 0.0 || lambda2 < 0.0);
            CHECK(has_negative);
        }
    }
}

TEST_CASE("compute_eigenvalues - zero matrix") {
    Matrix<2, 2> A = Matrix<2, 2>::zeros();

    auto result = compute_eigenvalues(A);

    // Zero matrix should converge immediately (already diagonal)
    // or may not converge if QR decomposition fails
    // Either way, eigenvalues should be zero
    CHECK(result.values[0].real() == doctest::Approx(0.0));
    CHECK(result.values[1].real() == doctest::Approx(0.0));
}

namespace {

// Largest |A - B| element over an N x M pair.
template<size_t N, size_t M>
double max_abs_diff(const Matrix<N, M>& A, const Matrix<N, M>& B) {
    double worst = 0.0;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            worst = damp::max(worst, std::abs(A(i, j) - B(i, j)));
        }
    }
    return worst;
}

// Returns {real, imag} eigenvalue pairs sorted by real then imag part, reading
// the QR/Schur result's diagonal real part and corresponding imaginary part.
template<size_t N>
std::array<damp::pair<double, double>, N> eig_pairs(const EigenResult<N, double>& r) {
    std::array<damp::pair<double, double>, N> out;
    for (size_t i = 0; i < N; ++i) {
        out[i] = {r.values[i].real(), r.values[i].imag()};
    }
    std::sort(out.begin(), out.end());
    return out;
}

} // namespace

TEST_CASE("qr_decompose reconstructs A and Q is orthonormal") {
    SUBCASE("square 3x3") {
        Matrix<3, 3> A = {
            {12.0, -51.0, 4.0},
            {6.0, 167.0, -68.0},
            {-4.0, 24.0, -41.0},
        };
        auto qr = qr_decompose(A);
        REQUIRE(qr.is_valid());

        // Q * R == A
        CHECK(max_abs_diff(qr.Q * qr.R, A) < 1e-9);

        // Qᵀ Q == I (orthonormal columns)
        Matrix<3, 3> QtQ = qr.Q.transpose() * qr.Q;
        CHECK(max_abs_diff(QtQ, Matrix<3, 3>::identity()) < 1e-9);

        // R is upper triangular
        CHECK(std::abs(qr.R(1, 0)) < 1e-9);
        CHECK(std::abs(qr.R(2, 0)) < 1e-9);
        CHECK(std::abs(qr.R(2, 1)) < 1e-9);
    }

    SUBCASE("tall 4x2") {
        Matrix<4, 2> A = {
            {1.0, 2.0},
            {3.0, 4.0},
            {5.0, 6.0},
            {7.0, 9.0},
        };
        auto qr = qr_decompose(A);
        REQUIRE(qr.is_valid());
        CHECK(max_abs_diff(qr.Q * qr.R, A) < 1e-9);
        // Columns of Q orthonormal: QᵀQ == I_2
        Matrix<2, 2> QtQ = qr.Q.transpose() * qr.Q;
        CHECK(max_abs_diff(QtQ, Matrix<2, 2>::identity()) < 1e-9);
    }
}

TEST_CASE("compute_eigenvalues resolves complex conjugate eigenvalues") {
    SUBCASE("2x2 pure rotation has eigenvalues ±i") {
        Matrix<2, 2> A = {
            {0.0, -1.0},
            {1.0, 0.0},
        };
        auto r = compute_eigenvalues(A);
        REQUIRE(r.converged);
        CHECK(std::abs(r.values[0].real()) < 1e-9);
        CHECK(std::abs(r.values[1].real()) < 1e-9);
        CHECK(std::abs(std::abs(r.values[0].imag()) - 1.0) < 1e-9);
        CHECK(r.values[0].imag() == doctest::Approx(-r.values[1].imag())); // conjugate pair
    }

    SUBCASE("2x2 damped oscillator: λ = -1 ± 2i") {
        // [[-1, -2], [2, -1]] has eigenvalues -1 ± 2i.
        Matrix<2, 2> A = {
            {-1.0, -2.0},
            {2.0, -1.0},
        };
        auto r = compute_eigenvalues(A);
        REQUIRE(r.converged);
        CHECK(r.values[0].real() == doctest::Approx(-1.0).epsilon(1e-9));
        CHECK(r.values[1].real() == doctest::Approx(-1.0).epsilon(1e-9));
        CHECK(std::abs(std::abs(r.values[0].imag()) - 2.0) < 1e-9);
    }

    SUBCASE("compute_eigenvalues now resolves complex pairs (Francis double-shift)") {
        // Previously the QR path returned only real parts. With the Francis
        // double-shift it resolves the 2x2 block into a true conjugate pair.
        Matrix<2, 2> A = {
            {0.0, -1.0},
            {1.0, 0.0},
        };
        auto r = compute_eigenvalues(A);
        REQUIRE(r.converged);
        CHECK(std::abs(r.values[0].real()) < 1e-9);
        CHECK(std::abs(r.values[1].real()) < 1e-9);
        CHECK(std::abs(std::abs(r.values[0].imag()) - 1.0) < 1e-9);
        CHECK(r.values[0].imag() == doctest::Approx(-r.values[1].imag()));
    }

    SUBCASE("compute_eigenvalues - 6x6 mixed real and two complex pairs") {
        // Block diagonal: rotation(4±3i), real 2 and -1, rotation(±5i).
        Matrix<6, 6> A = {
            {4.0, -3.0, 0.0, 0.0, 0.0, 0.0},
            {3.0, 4.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 2.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, -1.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0, 0.0, -5.0},
            {0.0, 0.0, 0.0, 0.0, 5.0, 0.0},
        };
        auto r = compute_eigenvalues(A);
        REQUIRE(r.converged);
        auto                                      got = eig_pairs(r);
        std::array<damp::pair<double, double>, 6> want = {{
            {-1.0, 0.0},
            {2.0, 0.0},
            {4.0, -3.0},
            {4.0, 3.0},
            {0.0, -5.0},
            {0.0, 5.0},
        }};
        std::ranges::sort(want);
        for (size_t i = 0; i < 6; ++i) {
            CHECK(got[i].first == doctest::Approx(want[i].first).epsilon(1e-9));
            CHECK(got[i].second == doctest::Approx(want[i].second).epsilon(1e-9));
        }
    }

    SUBCASE("companion matrix recovers polynomial roots 1..5") {
        // Companion of (x-1)(x-2)(x-3)(x-4)(x-5): eigenvalues 1,2,3,4,5.
        Matrix<5, 5> A = {
            {15.0, -85.0, 225.0, -274.0, 120.0},
            {1.0, 0.0, 0.0, 0.0, 0.0},
            {0.0, 1.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 1.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 1.0, 0.0},
        };
        auto r = compute_eigenvalues(A);
        REQUIRE(r.converged);
        auto got = eig_pairs(r);
        for (size_t i = 0; i < 5; ++i) {
            CHECK(got[i].first == doctest::Approx(static_cast<double>(i + 1)).epsilon(1e-6));
            CHECK(std::abs(got[i].second) < 1e-6);
        }
    }
}

TEST_CASE("compute_eigenvalues is constexpr-evaluable") {
    // Locks that the Hessenberg + Francis path contains no construct that
    // breaks constant evaluation (design-time eigenvalue computation).
    constexpr auto spectrum_ok = []() consteval {
        Matrix<3, 3> A = {
            {2.0, -1.0, 0.0},
            {-1.0, 2.0, -1.0},
            {0.0, -1.0, 2.0},
        };
        auto r = compute_eigenvalues(A);
        if (!r.converged) {
            return false;
        }
        // trace is preserved: sum of eigenvalues == 6
        double s = r.values[0].real() + r.values[1].real() + r.values[2].real();
        return damp::abs(s - 6.0) < 1e-9;
    };
    static_assert(spectrum_ok(), "compute_eigenvalues must work at compile time");
    CHECK(spectrum_ok());
}

TEST_CASE("compute_eigenvalues - 4x4 symmetric matches known spectrum") {
    // Tridiagonal [2 on diag, -1 off-diag] has eigenvalues
    // 2 - 2cos(k·π/5), k=1..4.
    Matrix<4, 4> A = {
        {2.0, -1.0, 0.0, 0.0},
        {-1.0, 2.0, -1.0, 0.0},
        {0.0, -1.0, 2.0, -1.0},
        {0.0, 0.0, -1.0, 2.0},
    };
    auto r = compute_eigenvalues(A);
    REQUIRE(r.converged);

    std::array<double, 4> want{};
    for (size_t k = 1; k <= 4; ++k) {
        want[k - 1] = 2.0 - (2.0 * std::cos(static_cast<double>(k) * std::numbers::pi / 5.0));
    }

    auto got = eig_pairs(r);
    std::ranges::sort(want);
    for (size_t i = 0; i < 4; ++i) {
        CHECK(got[i].first == doctest::Approx(want[i]).epsilon(1e-6));
        CHECK(std::abs(got[i].second) < 1e-6); // symmetric => real spectrum
    }
}
