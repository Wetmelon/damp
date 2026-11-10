// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>

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
constexpr double kEps = std::numeric_limits<double>::epsilon();

template<size_t R, size_t C>
double max_abs_diff(const Matrix<R, C, double>& a, const Matrix<R, C, double>& b) {
    double m = 0.0;
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            m = std::fmax(m, std::fabs(a(i, j) - b(i, j)));
        }
    }
    return m;
}

template<size_t M, size_t N>
Matrix<M, N, double> svd_reconstruct(const mat::SVDResult<M, N, double>& s) {
    Matrix<M, N, double> out{};
    constexpr size_t     K = (M < N) ? M : N;
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            double acc = 0.0;
            for (size_t k = 0; k < K; ++k) {
                acc += s.singular_U(i, k) * s.singular_values[k] * s.singular_V(j, k);
            }
            out(i, j) = acc;
        }
    }
    return out;
}

} // namespace

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

        CHECK(s.success);
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

        CHECK(s.success);
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

        CHECK(s.success);
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
        static_assert(s.success);
        // Largest singular value of diag(2,3) is 3.
        static_assert(s.singular_values[0] > 2.9 && s.singular_values[0] < 3.1);
    }

    TEST_CASE("SVD reconstructs to machine precision with tiny singular values") {
        // A = U diag(1, 1e-13, 1e-14) Vᵀ with dense orthogonal factors.
        const auto           Q = mat::svd(Matrix<3, 3, double>{{2.0, -1.0, 0.5}, {1.0, 3.0, -1.0}, {0.5, 1.0, 4.0}});
        Matrix<3, 3, double> S{};
        S(0, 0) = 1.0;
        S(1, 1) = 1e-13;
        S(2, 2) = 1e-14;
        const Matrix<3, 3, double> A = Q.singular_U * S * Q.singular_V.transpose();
        const auto                 s = mat::svd(A);
        CHECK(max_abs_diff<3, 3>(svd_reconstruct<3, 3>(s), A) <= 64 * kEps);
        CHECK(s.singular_values[1] == doctest::Approx(1e-13).epsilon(1e-6));
    }

    TEST_CASE("SVD: zero, rank-deficient, wide, badly scaled columns") {
        const auto z = mat::svd(Matrix<3, 2, double>{});
        CHECK(z.success);
        CHECK(z.singular_values[0] == 0.0);
        CHECK(mat::rank_from_svd(z) == 0);
        CHECK(max_abs_diff<3, 3>(z.singular_U * z.singular_U.transpose(), Matrix<3, 3, double>::identity()) <= 8 * kEps);

        // Rank 2 (third row = first + second), wide 3×4.
        const Matrix<3, 4, double> W{{1.0, 2.0, 3.0, 4.0}, {-1.0, 0.5, 2.0, 1.0}, {0.0, 2.5, 5.0, 5.0}};
        const auto                 sw = mat::svd(W);
        CHECK(sw.success);
        CHECK(mat::rank_from_svd(sw) == 2);
        CHECK(max_abs_diff<3, 4>(svd_reconstruct<3, 4>(sw), W) <= 64 * kEps * 10.0);
        CHECK(max_abs_diff<4, 4>(sw.singular_V.transpose() * sw.singular_V, Matrix<4, 4, double>::identity()) <= 64 * kEps);
        CHECK(max_abs_diff<3, 3>(sw.singular_U.transpose() * sw.singular_U, Matrix<3, 3, double>::identity()) <= 64 * kEps);

        // Columns graded over 30 decades: each singular value keeps relative accuracy.
        const Matrix<3, 3, double> G{{1.0, 1e-15, 1e-30}, {0.5, -1e-15, 2e-30}, {0.25, 3e-15, -1e-30}};
        const auto                 sg = mat::svd(G);
        CHECK(mat::rank_from_svd(sg, 0.0) == 3);
        for (size_t j = 0; j < 3; ++j) {
            double colmax = 0.0;
            for (size_t i = 0; i < 3; ++i) {
                double acc = 0.0;
                for (size_t k = 0; k < 3; ++k) {
                    acc += sg.singular_U(i, k) * sg.singular_values[k] * sg.singular_V(j, k);
                }
                colmax = std::fmax(colmax, std::fabs(acc - G(i, j)) / std::fabs(G(0, j)));
            }
            CHECK(colmax <= 64 * kEps);
        }
    }

    TEST_CASE("SVD of a complex matrix") {
        using Cd = damp::complex<double>;
        const Matrix<3, 2, Cd> A{{Cd{1.0, 1.0}, Cd{2.0, 0.0}}, {Cd{0.0, -1.0}, Cd{1.0, 3.0}}, {Cd{2.0, 0.5}, Cd{-1.0, 1.0}}};
        const auto             s = mat::svd(A);
        CHECK(s.success);
        double err = 0.0;
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                Cd acc{};
                for (size_t k = 0; k < 2; ++k) {
                    acc += s.singular_U(i, k) * s.singular_values[k] * damp::conj(s.singular_V(j, k));
                }
                err = std::fmax(err, (acc - A(i, j)).abs());
            }
        }
        CHECK(err <= 64 * kEps * 5.0);
        // σ_max of a scaled complex matrix scales exactly.
        const auto big = mat::svd(A * 1e200);
        CHECK(big.singular_values[0] == doctest::Approx(s.singular_values[0] * 1e200));
    }

    TEST_CASE("constexpr complex SVD matches runtime (small, nearly orthogonal columns)") {
        // [λI − A, B]-shaped 2×3 with millisecond-scale entries: |γ| ≪ √(αβ),
        // so the Jacobi angle needs √(ζ² + 1) with ζ far from 1.
        using Cd = damp::complex<double>;
        constexpr Matrix<2, 3, Cd> M{
            {Cd{1e-3, 1e-3}, Cd{-1e-3, 0.0}, Cd{0.0, 0.0}},
            {Cd{5e-4, 0.0}, Cd{1e-3, -2e-3}, Cd{1e-9, 0.0}},
        };
        constexpr auto cs = mat::svd(M);
        static_assert(cs.success);
        const auto rs = mat::svd(M);
        REQUIRE(rs.success);
        for (size_t k = 0; k < 2; ++k) {
            CHECK(std::fabs(cs.singular_values[k] - rs.singular_values[k]) <= 64 * kEps * rs.singular_values[0]);
        }
        CHECK(cs.singular_values[1] > 1e-4); // full rank, not the 1e-13 garbage
    }

    TEST_CASE("Column norms whose ratio underflows: the small column is still orthogonalized") {
        // The norm ratio ρ underflows to 0, and the pair was skipped and reported converged with
        // UᵀU − I at 0.707 and σ_min off by √2. The small column takes the limiting Gram–Schmidt step.
        constexpr Matrix<2, 2, double> A{{1e200, 1e-120}, {0.0, 1e-120}};
        static_assert(mat::svd(A).success);
        const auto s = mat::svd(A);
        REQUIRE(s.success);
        CHECK(std::fabs(s.singular_values[0] - 1e200) <= 1e-12 * 1e200);
        CHECK(std::fabs(s.singular_values[1] - 1e-120) <= 1e-12 * 1e-120);
        CHECK(unitary_error(s.singular_U) <= 8 * kEps);
        CHECK(unitary_error(s.singular_V) <= 8 * kEps);

        const Matrix<2, 2, float> Af{{1e30f, 1e-10f}, {0.0f, 1e-10f}};
        const auto                sf = mat::svd(Af);
        REQUIRE(sf.success);
        CHECK(std::fabs(sf.singular_values[1] - 1e-10f) <= 1e-5f * 1e-10f);
        CHECK(unitary_error(sf.singular_U) <= 8 * static_cast<double>(eps<float>()));
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

    TEST_CASE("Kernel of a wide matrix: orthogonal to the rows to a few ε·σ₁") {
        // The completion of U projected each candidate once: the kernel leaked up to ~30·ε·σ₁
        // into the row space at 3×20. A second pass on the chosen column restores ε.
        std::mt19937 gen(5); // NOLINT raw mt19937 output is fully specified, so the matrices are portable
        double       worst = 0.0;
        for (int trial = 0; trial < 20; ++trial) {
            Matrix<3, 20, double> A{};
            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = 0; j < 20; ++j) {
                    A(i, j) = (static_cast<double>(gen()) / 4294967296.0) - 0.5;
                }
            }
            const auto ns = null_space(A);
            REQUIRE(ns.dim == 17);
            const double s1 = svd(A).singular_values[0];
            for (size_t k = 20 - ns.dim; k < 20; ++k) {
                Matrix<20, 1, double> x{};
                for (size_t i = 0; i < 20; ++i) {
                    x(i, 0) = ns.vectors(i, k);
                }
                worst = std::max(worst, (A * x).norm() / (kEps * s1));
            }
        }
        CHECK(worst <= 8.0);
    }
}

TEST_SUITE("svd scaling") {
    // A = B·diag(1e19, 1, 0.1) in float: the two small columns sit 1e-19 and 1e-20 below the
    // largest entry. The sweep used to scale A by its largest entry and test Gram sums ‖a‖²;
    // for those columns the sums fell below min/eps and the columns were reported as zero, so
    // σ₂ and σ₃ came back as 0. One-sided Jacobi gets the singular values of a column-graded
    // matrix to relative precision (Demmel–Veselić); the double SVD of the same entries is the
    // reference.
    TEST_CASE("a float column-graded matrix keeps its small singular values") {
        const Matrix<3, 3, double> B{{1.0, 0.3, -0.2}, {0.4, 1.1, 0.5}, {-0.1, 0.2, 0.9}};
        const double               d[3] = {1e19, 1.0, 0.1};
        Matrix<3, 3, float>        Af{};
        Matrix<3, 3, double>       Ad{};
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                Af(i, j) = static_cast<float>(B(i, j) * d[j]);
                Ad(i, j) = static_cast<double>(Af(i, j));
            }
        }
        const auto sf = svd(Af);
        const auto sd = svd(Ad);
        REQUIRE(sd.singular_values[2] > 0.0);
        for (size_t k = 0; k < 3; ++k) {
            const double got = static_cast<double>(sf.singular_values[k]);
            const double ref = sd.singular_values[k];
            CHECK(damp::abs(got - ref) <= 1e-5 * ref);
        }
    }

    TEST_CASE("SVD is invariant to scaling the input") {
        const Matrix<3, 2, double> A{{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
        const auto                 ref = mat::svd(A);
        for (const double scale : {1e-160, 1e160}) {
            const auto s = mat::svd(A * scale);
            CHECK(s.success);
            CHECK(s.singular_values[0] == doctest::Approx(ref.singular_values[0] * scale));
            CHECK(s.singular_values[1] == doctest::Approx(ref.singular_values[1] * scale));
        }
    }

    // The prescale multiplied A by 1/max|a|. Past 1/min (8.5e37 in float, 4.5e307 in double) that
    // reciprocal is subnormal and -ffast-math flushed it to zero, so A became 0 and every σ came
    // back 0 as a success; for a subnormal max|a| it overflowed. The scale is now an exact power of
    // two, and only as large as the column norms need, so a graded matrix keeps its small σ.
    // volatile keeps the inputs out of constant folding, so this runs the suite's -ffast-math code.
    TEST_CASE_TEMPLATE("singular values near the top of the range are kept", T, float, double) {
        const bool      f = damp::is_same_v<T, float>;
        volatile T      top_in = f ? static_cast<T>(1e38) : static_cast<T>(1e308);
        volatile T      big_in = f ? static_cast<T>(3e38) : static_cast<T>(1.5e308);
        const T         top = top_in;
        const T         big = big_in;
        const T         tenth = static_cast<T>(0.1);
        const T         e = damp::eps<T>();
        const auto      near = [&](T got, T want) { return damp::abs(got - want) <= T{4} * e * want; };
        Matrix<2, 2, T> A{};
        A(0, 0) = top;
        A(1, 1) = top * tenth;
        const auto s = svd(A);
        REQUIRE(s.success);
        CHECK(near(s.singular_values[0], top));
        CHECK(near(s.singular_values[1], top * tenth));
        CHECK(rank_from_svd(s) == 2);

        Matrix<3, 3, T> G{}; // graded: σ₂ and σ₃ sit 38 (float) or 308 (double) decades below σ₁
        G(0, 0) = big;
        G(1, 1) = T{1};
        G(2, 2) = tenth;
        const auto sg = svd(G);
        REQUIRE(sg.success);
        CHECK(near(sg.singular_values[0], big));
        CHECK(near(sg.singular_values[1], T{1}));
        CHECK(near(sg.singular_values[2], tenth));

        Matrix<2, 2, complex<T>> C{};
        C(0, 0) = complex<T>{top, T{0}};
        C(1, 1) = complex<T>{T{0}, top * tenth};
        const auto sc = svd(C);
        REQUIRE(sc.success);
        CHECK(near(sc.singular_values[0], top));
        CHECK(near(sc.singular_values[1], top * tenth));
    }

    TEST_CASE_TEMPLATE("a singular value above the largest finite value fails the SVD", T, float, double) {
        // σ of [[2, 1], [1, 1]] are (3 ± √5)/2; σ₁ of [[3, 1], [2, 1]] is 3.86. With u = 1e38 in
        // float and 5e307 in double the first fits in T and the second does not: it came back
        // +Inf (or, before that, 0) with success = true.
        volatile T unit_in = damp::is_same_v<T, float> ? static_cast<T>(1e38) : static_cast<T>(5e307);
        const T    u = unit_in;
        const auto fits = svd(Matrix<2, 2, T>{{T{2} * u, u}, {u, u}});
        REQUIRE(fits.success);
        const long double root5 = std::sqrt(5.0L);
        const long double s0 = (3.0L + root5) / 2.0L * static_cast<long double>(u);
        const long double s1 = (3.0L - root5) / 2.0L * static_cast<long double>(u);
        const long double tol = 8.0L * static_cast<long double>(damp::eps<T>());
        CHECK(std::fabs(static_cast<long double>(fits.singular_values[0]) - s0) <= tol * s0);
        CHECK(std::fabs(static_cast<long double>(fits.singular_values[1]) - s1) <= tol * s0);
        const auto over = svd(Matrix<2, 2, T>{{T{3} * u, u}, {T{2} * u, u}});
        CHECK_FALSE(over.success);
        CHECK_FALSE(damp::isfinite(over.singular_values[0]));
        CHECK(rank_from_svd(over) == 0);
    }

    TEST_CASE("subnormal singular values come back exactly, constant-evaluated") {
        // Constant evaluation keeps subnormals (no flush-to-zero); 1/max|a| overflowed there.
        constexpr auto sd = svd(Matrix<2, 2, double>{{3e-310, 0.0}, {0.0, 1e-310}});
        static_assert(sd.success && sd.singular_values[0] == 3e-310 && sd.singular_values[1] == 1e-310);
        constexpr auto sf = svd(Matrix<2, 2, float>{{3e-40f, 0.0f}, {0.0f, 1e-40f}});
        static_assert(sf.success && sf.singular_values[0] == 3e-40f && sf.singular_values[1] == 1e-40f);
        CHECK(sd.success);
    }

    TEST_CASE("a NaN or ±Inf entry fails the SVD, and rank, pinv, null space and two-norm do not read it") {
        // Regression: the sweeps zeroed a NaN column and reported success with finite,
        // plausible singular values, so rank and pinv returned plain garbage.
        const double               nan = std::numeric_limits<double>::quiet_NaN();
        const double               inf = std::numeric_limits<double>::infinity();
        const Matrix<3, 3, double> A{{1.0, 2.0, 3.0}, {4.0, nan, 6.0}, {7.0, 8.0, 10.0}};
        const Matrix<2, 3, double> W{{1.0, inf, 0.0}, {0.0, 1.0, 0.0}};
        const auto                 s = mat::svd(A);
        CHECK_FALSE(s.success);
        CHECK_FALSE(static_cast<bool>(s));
        CHECK_FALSE(damp::isfinite(s.singular_values[0]));
        CHECK_FALSE(mat::svd(W).success);
        CHECK(mat::rank(A) == 0);
        CHECK(mat::rank(W) == 0);
        CHECK_FALSE(damp::isfinite(mat::pseudo_inverse(A)(0, 0)));
        CHECK_FALSE(damp::isfinite(mat::two_norm(A)));
        const auto ns = mat::null_space(A);
        CHECK(ns.dim == 3);
        CHECK_FALSE(damp::isfinite(ns.vectors(0, 0)));
        const Matrix<2, 2, complex<double>> C{{complex<double>{1.0, nan}, complex<double>{0.0}}, {complex<double>{0.0}, complex<double>{1.0}}};
        CHECK_FALSE(mat::svd(C).success);
        CHECK(mat::svd(Matrix<3, 3, double>{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 10.0}}).success);
        CHECK_FALSE(mat::SVDResult<2, 2, double>{}.success);
    }
}
