// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
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

TEST_SUITE("compute_eigenvalues balancing") {
    // T = diag(1, 1e6, 1e-6) applied as a similarity leaves the spectrum alone but puts entries
    // of 1e12 next to entries of 1e-12. Without balancing the QR iteration's backward error
    // eps·‖A‖ ≈ 2e-4 lands on eigenvalues of size 1; with it the eigenvalues come back to
    // relative precision, and the first Schur vector is still an eigenvector.
    TEST_CASE("a badly scaled similarity keeps the eigenvalues and the leading eigenvector") {
        const Matrix<3, 3> A0{{2.0, 1.0, 0.0}, {0.0, 3.0, 1.0}, {0.5, 0.0, 5.0}}; // real spectrum
        const double       s[3] = {1.0, 1.0e6, 1.0e-6};
        Matrix<3, 3>       A{};
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                A(i, j) = s[i] * A0(i, j) / s[j];
            }
        }
        const auto ref = compute_eigenvalues(A0);
        const auto got = compute_eigenvalues(A);
        REQUIRE(ref.converged);
        REQUIRE(got.converged);
        for (size_t i = 0; i < 3; ++i) {
            // match each eigenvalue of A to the nearest of A0
            double best = 1e300;
            for (size_t k = 0; k < 3; ++k) {
                best = std::min(best, (got.values[i] - ref.values[k]).abs());
            }
            CHECK(best <= 1e-10 * ref.values[i].abs());
        }
        // The first Schur vector is an eigenvector of A: ‖A v − λ v‖ small relative to ‖A v‖.
        ColVec<3, double> v{got.vectors(0, 0).real(), got.vectors(1, 0).real(), got.vectors(2, 0).real()};
        const auto        Av = A * v;
        const double      lam = got.values[0].real();
        double            res = 0.0;
        double            scale = 0.0;
        for (size_t i = 0; i < 3; ++i) {
            res = std::max(res, std::fabs(Av[i] - (lam * v[i])));
            scale = std::max(scale, std::fabs(Av[i]));
        }
        CHECK(res <= 1e-9 * scale);
    }
}

TEST_SUITE("Francis QR on the Hessenberg form") {
    // Run the Schur reduction without balancing, so the test sees the iteration itself.
    template<size_t N>
    struct Schur {
        Matrix<N, N>           T{};
        Matrix<N, N>           Z{};
        damp::array<double, N> wr{};
        damp::array<double, N> wi{};
        bool                   converged{false};
    };

    template<size_t N>
    Schur<N> schur(const Matrix<N, N>& A) {
        Schur<N> s;
        s.T = A;
        mat::detail::hessenberg_reduce(s.T, s.Z);
        s.converged = mat::detail::francis_qr(s.T, s.Z, s.wr, s.wi);
        return s;
    }

    // h₃₂ = 1e-16 passes the classic test |h₃₂| ≤ ε(|h₂₂| + |h₃₃|) and used to be zeroed,
    // returning h₃₃ = 1e-20 as an eigenvalue. The eigenvalue is 1e-20 − h₂₃h₃₂/(h₂₂ − h₃₃)
    // ≈ −9e-20: the product h₂₃h₃₂ = 1e-19 is what moves it, and Ahues–Tisseur compares that
    // product, not h₃₂ alone, against ε|h₃₃||h₂₂ − h₃₃|.
    TEST_CASE("a graded Hessenberg matrix keeps its tiny eigenvalue (Ahues–Tisseur)") {
        const Matrix<3, 3>     H{{2.0, 1.0, 1.0}, {0.0, 1.0, 1e-3}, {0.0, 1e-16, 1e-20}};
        auto                   T = H;
        auto                   Z = Matrix<3, 3>::identity();
        damp::array<double, 3> wr{};
        damp::array<double, 3> wi{};
        REQUIRE(mat::detail::francis_qr(T, Z, wr, wi));
        double smallest = wr[0];
        for (const double w : wr) {
            if (std::fabs(w) < std::fabs(smallest)) {
                smallest = w;
            }
        }
        const double exact = 1e-20 - (1e-19 / (1.0 - 1e-20 + 1e-19));
        CHECK(std::fabs(smallest - exact) <= 1e-12 * std::fabs(exact)); // relative: Approx would add an absolute 1e-12
    }

    // A cyclic permutation has all its eigenvalues on the unit circle and stalls the Francis
    // shift; only the exceptional shifts move it. Two cycles take the schedule once per
    // eigenvalue (it passed before the count restarted at each deflation, too).
    TEST_CASE("two permutation cycles converge through exceptional shifts") {
        Matrix<6, 6> P{};
        P(1, 0) = P(2, 1) = P(0, 2) = 1.0; // 3-cycle
        P(4, 3) = P(5, 4) = P(3, 5) = 1.0; // 3-cycle
        const auto s = schur(P);
        REQUIRE(s.converged);
        for (size_t i = 0; i < 6; ++i) {
            CHECK(std::hypot(s.wr[i], s.wi[i]) == doctest::Approx(1.0).epsilon(1e-12));
        }
    }

    // Every 2×2 block left by the iteration is in standard form: a complex pair with equal
    // diagonals and off-diagonals of opposite sign; real pairs are triangular. Zᵀ A Z = T with
    // Z orthogonal, and the returned eigenvalues are the ones T displays.
    TEST_CASE("converged 2×2 blocks are standardized") {
        unsigned seed = 12345U;
        auto     rnd = [&seed]() {
            seed = (seed * 1103515245U) + 12345U;
            return (static_cast<double>((seed >> 8) & 0xFFFFU) / 32768.0) - 1.0;
        };
        for (int trial = 0; trial < 200; ++trial) {
            Matrix<6, 6> A{};
            for (size_t i = 0; i < 6; ++i) {
                for (size_t j = 0; j < 6; ++j) {
                    A(i, j) = rnd();
                }
            }
            const auto s = schur(A);
            REQUIRE(s.converged);
            const Matrix<6, 6> back = s.Z.transpose() * A * s.Z;
            const Matrix<6, 6> ztz = s.Z.transpose() * s.Z;
            for (size_t i = 0; i < 6; ++i) {
                for (size_t j = 0; j < 6; ++j) {
                    CHECK(back(i, j) == doctest::Approx(s.T(i, j)).scale(1.0).epsilon(1e-12));
                    CHECK(ztz(i, j) == doctest::Approx(i == j ? 1.0 : 0.0).scale(1.0).epsilon(1e-13));
                    if (i > j + 1) {
                        CHECK(s.T(i, j) == 0.0);
                    }
                }
            }
            for (size_t i = 0; i + 1 < 6; ++i) {
                if (s.T(i + 1, i) == 0.0) {
                    CHECK(s.wi[i] == 0.0);
                    CHECK(s.wr[i] == s.T(i, i));
                    continue;
                }
                CHECK(s.T(i, i) == s.T(i + 1, i + 1));
                CHECK(s.T(i, i + 1) * s.T(i + 1, i) < 0.0);
                CHECK(s.wr[i] == s.T(i, i));
                CHECK(s.wi[i] == doctest::Approx(std::sqrt(-s.T(i, i + 1) * s.T(i + 1, i))));
                CHECK(s.wi[i + 1] == -s.wi[i]);
                if (i + 2 < 6) {
                    CHECK(s.T(i + 2, i + 1) == 0.0); // blocks never overlap
                }
                ++i;
            }
        }
    }
}

TEST_SUITE("compute_eigenvalues balancing in float") {
    // The similarity diag(1, 1e10, 1e-10) puts entries of 1e20 next to entries of 1e-20. Their
    // squares overflow and underflow float, so a balancing step built on plain sums of squares
    // saw infinite row norms and left the matrix unbalanced: the QR iteration then carried a
    // backward error of ε‖A‖ ≈ 1e14 onto eigenvalues of size 1. Scaled norms balance it.
    TEST_CASE("entries whose squares overflow float are still balanced") {
        const Matrix<3, 3, float> A0{{2.0f, 1.0f, 0.0f}, {0.0f, 3.0f, 1.0f}, {0.5f, 0.0f, 5.0f}};
        const float               s[3] = {1.0f, 1.0e10f, 1.0e-10f};
        Matrix<3, 3, float>       A{};
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                A(i, j) = s[i] * A0(i, j) / s[j];
            }
        }
        const auto ref = compute_eigenvalues(A0);
        const auto got = compute_eigenvalues(A);
        REQUIRE(ref.converged);
        REQUIRE(got.converged);
        for (size_t i = 0; i < 3; ++i) {
            float best = 1e30f;
            for (size_t k = 0; k < 3; ++k) {
                best = damp::min(best, (got.values[i] - ref.values[k]).abs());
            }
            CHECK(best <= 1e-4f * ref.values[i].abs());
        }
    }

    TEST_CASE("a non-finite entry leaves the matrix unbalanced") {
        Matrix<2, 2> A{{1.0, 2.0}, {3.0, 4.0}};
        A(0, 1) = std::numeric_limits<double>::quiet_NaN();
        damp::array<int, 2> e{};
        CHECK_FALSE(mat::detail::gebal_log2_scale(A, e));
        CHECK(e[0] == 0);
        CHECK(e[1] == 0);
    }
}
