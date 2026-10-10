// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Regression tests for the geometry / scalar-math / matrix-core review (group g1).

#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>

#include "damp/math/complex.hpp"
#include "damp/math/geometry.hpp"
#include "damp/math/math.hpp"
#include "damp/math/trig.hpp"
#include "damp/matrix/matrix.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

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

template<EulerOrder Order>
void check_euler_round_trip(double a1, double a2, double a3) {
    const Euler<double, Order> e{a1, a2, a3};
    const DCM<double>          R = e.to_dcm();
    const DCM<double>          R_back = Euler<double, Order>::from_dcm(R).to_dcm();
    CHECK(max_abs_diff<3, 3>(R, R_back) <= 64 * kEps);
}

template<EulerOrder Order>
void check_euler_order() {
    constexpr double half_pi = std::numbers::pi / 2.0;
    // Generic, near gimbal lock (cos θ ≈ 1e-3 and 1e-6), and exactly at it.
    for (const double a2 : {0.4, -1.1, half_pi - 1e-3, -(half_pi - 1e-3), half_pi - 1e-6, half_pi, -half_pi}) {
        check_euler_round_trip<Order>(0.3, a2, 0.7);
        check_euler_round_trip<Order>(-2.5, a2, -0.2);
    }
}

template<size_t N>
double eigen_residual(const Matrix<N, N, double>& A, const mat::EigenResult<N, double>& r, size_t i) {
    double res = 0.0;
    for (size_t row = 0; row < N; ++row) {
        double av = 0.0;
        for (size_t k = 0; k < N; ++k) {
            av += A(row, k) * r.vectors(k, i).real();
        }
        res = std::fmax(res, std::fabs(av - (r.values[i].real() * r.vectors(row, i).real())));
    }
    return res;
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

TEST_SUITE("review g1: constexpr scalar math") {
    TEST_CASE("series sqrt and cbrt converge across the exponent range") {
        for (const double x : {1e-300, 1e-100, 1e-10, 0.3, 2.0, 1e10, 1e100, 1e300}) {
            CHECK(std::fabs(detail::sqrt(x) - std::sqrt(x)) <= 2 * kEps * std::sqrt(x));
            CHECK(std::fabs(detail::cbrt(x) - std::cbrt(x)) <= 2 * kEps * std::cbrt(x));
            CHECK(std::fabs(detail::cbrt(-x) + std::cbrt(x)) <= 2 * kEps * std::cbrt(x));
        }
        CHECK(std::fabs(detail::sqrt(1e30f) - 1e15f) <= 2 * std::numeric_limits<float>::epsilon() * 1e15f);
    }

    TEST_CASE("constant-evaluated roots are accurate far from 1, subnormals included") {
        constexpr double denorm = 4.9406564584124654e-324; // smallest subnormal
        static_assert(detail::abs(damp::sqrt(1e150) - 1e75) <= 2 * kEps * 1e75);
        static_assert(detail::abs(damp::sqrt(1e-150) - 1e-75) <= 2 * kEps * 1e-75);
        static_assert(detail::abs(damp::sqrt(1e300) - 1e150) <= 2 * kEps * 1e150);
        static_assert(detail::abs(damp::sqrt(denorm) - 2.2227587494850775e-162) <= 2 * kEps * 2.3e-162);
        static_assert(detail::abs(damp::cbrt(1e-300) - 1e-100) <= 2 * kEps * 1e-100);
        static_assert(detail::abs(damp::cbrt(-8e150) + 2e50) <= 2 * kEps * 2e50);
        static_assert(detail::abs(damp::sqrt(1e30f) - 1e15f) <= 2 * std::numeric_limits<float>::epsilon() * 1e15f);
        static_assert(detail::abs(damp::log(1e300) - 690.77552789821368) <= 4 * kEps * 690.8);
        static_assert(detail::abs(damp::log(denorm) + 744.44007192138122) <= 4 * kEps * 744.5);
        static_assert(detail::abs(damp::hypot(3e300, 4e300) - 5e300) <= 4 * kEps * 5e300);
    }

    TEST_CASE("series log/exp/pow/atan2 hold accuracy far from 1") {
        for (const double x : {1e-300, 1e-150, 1e-10, 0.7, 1.0, 1.3, 1e10, 1e150, 1e300}) {
            CHECK(std::fabs(detail::log(x) - std::log(x)) <= (4 * kEps * std::fabs(std::log(x))) + kEps);
        }
        for (const double x : {-700.0, -300.0, -1.0, 0.5, 1.0, 300.0, 700.0}) {
            CHECK(std::fabs(detail::exp(x) - std::exp(x)) <= 8 * kEps * std::exp(x));
        }
        // exp(y·ln b): the rounding of y·ln b is amplified by |y·ln b| (≈ 690 here).
        CHECK(std::fabs(detail::pow(10.0, 300.0) - 1e300) <= 1024 * kEps * 1e300);
        CHECK(std::fabs(detail::pow(10.0, -300.0) - 1e-300) <= 1024 * kEps * 1e-300);
        CHECK(std::fabs(detail::atan2(1e300, 1e-300) - (std::numbers::pi / 2.0)) <= 2 * kEps);
        CHECK(std::fabs(detail::atan2(1e-300, 1e300)) <= 1e-300);
        CHECK(std::fabs(detail::atan2(-1e-300, -1e-300) + (0.75 * std::numbers::pi)) <= 4 * kEps);
    }

    TEST_CASE("series sin/cos/sincos are accurate over several periods") {
        double worst = 0.0;
        for (int i = -2000; i <= 2000; ++i) {
            const double x = 0.01 * i;
            const auto   sc = detail::sincos(x);
            worst = std::fmax(worst, std::fabs(detail::sin(x) - std::sin(x)));
            worst = std::fmax(worst, std::fabs(detail::cos(x) - std::cos(x)));
            worst = std::fmax(worst, std::fabs(sc.sin - std::sin(x)));
            worst = std::fmax(worst, std::fabs(sc.cos - std::cos(x)));
        }
        CHECK(worst <= 4 * kEps);
        CHECK(std::fabs(detail::cos(std::numbers::pi) + 1.0) <= 2 * kEps);
        // Same kernels under constant evaluation (the consteval design path).
        static_assert(detail::abs(damp::cos(std::numbers::pi) + 1.0) <= 2 * kEps);
        static_assert(detail::abs(damp::sqrt(1e300) - 1e150) <= 2 * kEps * 1e150);
    }

    TEST_CASE("series tan keeps relative accuracy next to pi/2") {
        for (const double x : {1.5707963267, 1.57079632679, -1.5707963267, 1.0, 0.1, 4.0}) {
            CHECK(std::fabs(detail::tan(x) - std::tan(x)) <= 16 * kEps * std::fabs(std::tan(x)));
        }
    }

    TEST_CASE("integer pow handles INT_MIN without signed overflow") {
        static_assert(damp::pow(1.0, std::numeric_limits<int>::min()) == 1.0);
        static_assert(damp::pow(-1.0, std::numeric_limits<int>::min()) == 1.0);
        CHECK(damp::pow(2.0, std::numeric_limits<int>::min()) == 0.0);
    }

    TEST_CASE("fast atan2 is exact in ratio for tiny normal arguments") {
        constexpr float pi_4 = std::numbers::pi_v<float> / 4.0f;
        CHECK(std::fabs(detail::fast_atan2(1e-36f, 1e-36f) - pi_4) <= 1e-5f);
        CHECK(std::fabs(detail::fast_atan2(1e-306, 1e-306) - (std::numbers::pi / 4.0)) <= 4e-12); // poly accuracy
    }
}

TEST_SUITE("review g1: complex") {
    using C = damp::complex<double>;

    TEST_CASE("division does not overflow or underflow in the denominator") {
        const C big = C{1e200, 1e200} / C{1e200, 1e200};
        CHECK(big.real() == doctest::Approx(1.0));
        CHECK(std::fabs(big.imag()) <= kEps);
        const C small = C{1e-200, 2e-200} / C{1e-200, 0.0};
        CHECK(small.real() == doctest::Approx(1.0));
        CHECK(small.imag() == doctest::Approx(2.0));
        C q{3e200, 4e200};
        q /= C{0.0, 1e200};
        CHECK(q.real() == doctest::Approx(4.0));
        CHECK(q.imag() == doctest::Approx(-3.0));
        const C r = 1e-300 / C{1e-300, 1e-300};
        CHECK(r.real() == doctest::Approx(0.5));
        CHECK(r.imag() == doctest::Approx(-0.5));
    }

    TEST_CASE("abs does not overflow") {
        CHECK(C{3e200, 4e200}.abs() == doctest::Approx(5e200));
        CHECK(damp::abs(C{3e-200, 4e-200}) == doctest::Approx(5e-200));
    }

    TEST_CASE("sqrt keeps the small component (no cancellation)") {
        const C s = damp::sqrt(C{1.0, 1e-10});
        CHECK(s.real() == doctest::Approx(1.0));
        CHECK(s.imag() == doctest::Approx(5e-11));
        const C t = damp::sqrt(C{-1.0, -1e-10});
        CHECK(t.real() == doctest::Approx(5e-11));
        CHECK(t.imag() == doctest::Approx(-1.0));
    }
}

TEST_SUITE("review g1: geometry") {
    TEST_CASE("Euler extraction round-trips near and at gimbal lock, every order") {
        check_euler_order<EulerOrder::ZYX>();
        check_euler_order<EulerOrder::XYZ>();
        check_euler_order<EulerOrder::ZXY>();
        check_euler_order<EulerOrder::YXZ>();
        check_euler_order<EulerOrder::YZX>();
        check_euler_order<EulerOrder::XZY>();
    }

    TEST_CASE("quaternion normalization and inverse are scale-invariant") {
        const Quatd small{1e-5, 0.0, 0.0, 0.0};
        REQUIRE(small.normalized_safe().has_value());
        CHECK(small.normalized().w() == doctest::Approx(1.0));
        const auto inv = Quatd{0.0, 2e-5, 0.0, 0.0}.inverse();
        REQUIRE(inv.has_value());
        CHECK(inv->x() == doctest::Approx(-5e4));
        CHECK_FALSE(Quatd{0.0, 0.0, 0.0, 0.0}.normalized_safe().has_value());
        CHECK_FALSE(Quatd{0.0, 0.0, 0.0, 0.0}.inverse().has_value());
    }

    TEST_CASE("axis-angle accepts any nonzero axis and log keeps tiny rotations") {
        const Vec3d tiny_axis{0.0, 1e-12, 0.0};
        REQUIRE(Quatd::from_axis_angle(tiny_axis, 0.3).has_value());
        REQUIRE(DCMd::from_axis_angle(tiny_axis, 0.3).has_value());
        const auto q = Quatd::from_axis_angle(Vec3d{0.0, 0.0, 1.0}, 1e-9);
        REQUIRE(q.has_value());
        const Vec3d rv = q->log();
        CHECK(rv[2] == doctest::Approx(1e-9));
        const auto [axis, angle] = q->to_axis_angle();
        CHECK(axis[2] == doctest::Approx(1.0));
        CHECK(angle == doctest::Approx(1e-9));
    }
}

TEST_SUITE("review g1: matrix aliasing") {
    TEST_CASE("compound assignment from a view of itself") {
        Matrix<2, 2, double> A{{1.0, 2.0}, {3.0, 4.0}};
        A += A.t();
        CHECK(A == Matrix<2, 2, double>{{2.0, 5.0}, {5.0, 8.0}});
        A -= A.t();
        CHECK(A == Matrix<2, 2, double>{{0.0, 0.0}, {0.0, 0.0}});
    }

    TEST_CASE("block assigned from an overlapping block temporary") {
        Matrix<3, 1, double> v{{1.0}, {2.0}, {3.0}};
        v.block<2, 1>(1, 0) = v.block<2, 1>(0, 0);
        CHECK(v == Matrix<3, 1, double>{{1.0}, {1.0}, {2.0}});
        Matrix<3, 1, double> w{{1.0}, {2.0}, {3.0}};
        auto                 dst = w.block<2, 1>(1, 0);
        dst += w.block<2, 1>(0, 0);
        CHECK(w == Matrix<3, 1, double>{{1.0}, {3.0}, {5.0}});
    }

    TEST_CASE("row/column view assigned from an overlapping view") {
        Matrix<3, 3, double> A{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
        A.row(2) = A.col(0);
        CHECK(A.row_vector(2) == Matrix<1, 3, double>{{1.0, 4.0, 7.0}});
        Matrix<3, 3, double> B{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
        B.col(2) = B.row(0);
        CHECK(B.col_vector(2) == Matrix<3, 1, double>{{1.0}, {2.0}, {3.0}});
    }
}

TEST_SUITE("review g1: norms, eigen, svd") {
    TEST_CASE("norm neither overflows nor underflows") {
        CHECK(ColVec<2, double>{3e200, 4e200}.norm() == doctest::Approx(5e200));
        CHECK(ColVec<2, double>{3e-200, 4e-200}.norm() == doctest::Approx(5e-200));
        CHECK(Matrix<1, 2, double>{{3e200, 4e200}}.norm() == doctest::Approx(5e200));
        CHECK(RowVec<2, double>{3e-200, 4e-200}.normalized()[1] == doctest::Approx(0.8));
    }

    TEST_CASE("2x2 eigenvalues: small eigenvalue keeps relative accuracy") {
        const auto   r = mat::compute_eigenvalues(Matrix<2, 2, double>{{1e8, 1.0}, {0.0, 1e-8}});
        const double lo = std::fmin(r.values[0].real(), r.values[1].real());
        CHECK(lo == doctest::Approx(1e-8).epsilon(1e-12));
    }

    TEST_CASE("Schur vectors are eigenvectors of symmetric matrices") {
        const Matrix<2, 2, double> A2{{2.0, 1.0}, {1.0, 2.0}};
        const auto                 r2 = mat::compute_eigenvalues(A2);
        for (size_t i = 0; i < 2; ++i) {
            CHECK(eigen_residual<2>(A2, r2, i) <= 64 * kEps * 3.0);
        }
        const Matrix<4, 4, double> A4{
            {4.0, 1.0, -2.0, 2.0},
            {1.0, 2.0, 0.0, 1.0},
            {-2.0, 0.0, 3.0, -2.0},
            {2.0, 1.0, -2.0, -1.0},
        };
        const auto r4 = mat::compute_eigenvalues(A4);
        REQUIRE(r4.converged);
        for (size_t i = 0; i < 4; ++i) {
            CHECK(eigen_residual<4>(A4, r4, i) <= 64 * kEps * 8.0);
        }
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
        CHECK(z.converged);
        CHECK(z.singular_values[0] == 0.0);
        CHECK(mat::rank_from_svd(z) == 0);
        CHECK(max_abs_diff<3, 3>(z.singular_U * z.singular_U.transpose(), Matrix<3, 3, double>::identity()) <= 8 * kEps);

        // Rank 2 (third row = first + second), wide 3×4.
        const Matrix<3, 4, double> W{{1.0, 2.0, 3.0, 4.0}, {-1.0, 0.5, 2.0, 1.0}, {0.0, 2.5, 5.0, 5.0}};
        const auto                 sw = mat::svd(W);
        CHECK(sw.converged);
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
        CHECK(s.converged);
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
        static_assert(cs.converged);
        const auto rs = mat::svd(M);
        REQUIRE(rs.converged);
        for (size_t k = 0; k < 2; ++k) {
            CHECK(std::fabs(cs.singular_values[k] - rs.singular_values[k]) <= 64 * kEps * rs.singular_values[0]);
        }
        CHECK(cs.singular_values[1] > 1e-4); // full rank, not the 1e-13 garbage
    }

    TEST_CASE("SVD is invariant to scaling the input") {
        const Matrix<3, 2, double> A{{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
        const auto                 ref = mat::svd(A);
        for (const double scale : {1e-160, 1e160}) {
            const auto s = mat::svd(A * scale);
            CHECK(s.converged);
            CHECK(s.singular_values[0] == doctest::Approx(ref.singular_values[0] * scale));
            CHECK(s.singular_values[1] == doctest::Approx(ref.singular_values[1] * scale));
        }
    }
}
