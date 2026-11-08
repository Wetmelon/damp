// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// PatternMatrix: zeros and ones known at compile time, the other entries runtime values. Every product
// has to be the dense one to rounding, and the filter on it has to be the filter on the dense Jacobians.

#include <cmath>
#include <cstddef>
#include <type_traits>

#include "damp/estimation/eskf.hpp"
#include "damp/estimation/ins_eskf.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/pattern_matrix.hpp"
#include "doctest.h"

using namespace damp;

namespace {

/// Deterministic, nothing special about the values: no accidental zeros or ones.
double wiggle(size_t i, size_t j, double seed) {
    return std::sin((1.3 * static_cast<double>(i)) + (0.7 * static_cast<double>(j)) + seed);
}

/// A symmetric positive definite matrix.
template<size_t N>
Matrix<N, N> spd(double seed) {
    Matrix<N, N> a{};
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            a(i, j) = wiggle(i, j, seed);
        }
    }
    Matrix<N, N> s = a * a.transpose();
    for (size_t i = 0; i < N; ++i) {
        s(i, i) += 1.0;
    }
    return s;
}

/// Fill every runtime value of a sparse matrix.
template<const auto& P>
void fill(PatternMatrix<P, double>& m, double seed) {
    for (size_t k = 0; k < PatternMatrix<P, double>::NV; ++k) {
        m.values[k] = 0.1 + (0.05 * std::sin((2.1 * static_cast<double>(k)) + seed));
    }
}

template<size_t R, size_t C>
double gap(const Matrix<R, C>& a, const Matrix<R, C>& b) {
    double w = 0.0;
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            w = std::fmax(w, std::abs(a(i, j) - b(i, j)));
        }
    }
    return w;
}

// Shapes that exercise the sandwich's cases.
inline constexpr auto kIns = kInsTransitionPattern;
inline constexpr auto kAllUnit = SparsityPattern<5, 5>::identity();
inline constexpr auto kNoUnit = SparsityPattern<4, 4>{}.values(0, 0, 4, 4);
inline constexpr auto kApartUnits = SparsityPattern<6, 6>::identity().value_diagonal(1, 4, 1).values(3, 0, 1, 6).values_off_diagonal(4, 0, 2);
inline constexpr auto kZeroRow = SparsityPattern<4, 4>{}.unit_diagonal(0, 0, 2).values(3, 1, 1, 2); // row 2 is all zero: F S Fᵀ loses that state
inline constexpr auto kSelector = SparsityPattern<3, 7>{}.unit_diagonal(0, 2, 3);
inline constexpr auto kLever = SparsityPattern<3, 7>{}.unit_diagonal(0, 2, 3).values(0, 5, 3, 2);
inline constexpr auto kMissesState = SparsityPattern<2, 5>{}.values(0, 1, 1, 2).value_diagonal(1, 4, 1);

template<const auto& P>
void check_sandwich(double seed) {
    constexpr size_t         N = PatternMatrix<P, double>::NR;
    PatternMatrix<P, double> F{};
    fill(F, seed);
    const Matrix<N, N> S = spd<N>(seed);
    const Matrix<N, N> dense = F.to_matrix() * S * F.to_matrix().transpose();
    const Matrix<N, N> out = quadratic_form(F, S);
    CHECK(gap(out, dense) <= 1e-12 * (1.0 + dense.norm()));
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            CHECK(out(i, j) == out(j, i)); // exactly
        }
    }
    Matrix<N, N> in_place = S;
    quadratic_form_in_place(F, in_place);
    CHECK(gap(in_place, out) <= 1e-12 * (1.0 + dense.norm())); // the same arithmetic, which fast-math may order apart
}

} // namespace

TEST_SUITE("SparsityPattern") {
    TEST_CASE("it says what each entry is, at compile time") {
        static constexpr auto p = SparsityPattern<3, 4>{}.unit_diagonal(0, 1, 2).values(2, 0, 1, 2).value_diagonal(0, 3, 1);
        static_assert(p(0, 1) == PatternEntry::one && p(1, 2) == PatternEntry::one);
        static_assert(p(2, 0) == PatternEntry::value && p(2, 1) == PatternEntry::value && p(0, 3) == PatternEntry::value);
        static_assert(p(0, 0) == PatternEntry::zero && p(1, 3) == PatternEntry::zero);
        static_assert(p.value_count() == 3);
        static_assert(detail::nonzero_columns<p, 0>.n == 2 && detail::nonzero_columns<p, 0>.col[1] == 3);
        static_assert(detail::nonzero_columns<p, 1>.n == 1);
        static_assert(detail::pattern_read_columns<p>[1] && detail::pattern_read_columns<p>[3]);
        static constexpr auto sparse_reads = SparsityPattern<1, 3>{}.unit_diagonal(0, 2, 1);
        static_assert(!detail::pattern_read_columns<sparse_reads>[0] && detail::pattern_read_columns<sparse_reads>[2]);
    }

    TEST_CASE("the off-diagonal builder leaves a diagonal one alone; a unit row is exactly eᵣ") {
        static constexpr auto p = SparsityPattern<4, 4>::identity().values_off_diagonal(0, 0, 3);
        static_assert(p(0, 0) == PatternEntry::one && p(0, 1) == PatternEntry::value && p(3, 3) == PatternEntry::one);
        static_assert(p.value_count() == 6);
        static_assert(!p.unit_row(0) && p.unit_row(3));
    }

    TEST_CASE("the INS transition has the shape the filter relies on") {
        static_assert(kInsTransitionPattern.value_count() == 30);
        static_assert(kInsTransitionPattern.unit_row(ins_err::dbg) && kInsTransitionPattern.unit_row(ins_err::dba + 2));
        static_assert(!kInsTransitionPattern.unit_row(ins_err::dp));
        static_assert(sizeof(PatternMatrix<kInsTransitionPattern, float>) == 30 * sizeof(float));
    }
}

TEST_SUITE("PatternMatrix") {
    TEST_CASE("a known entry has no storage, and reads as what the pattern says") {
        static constexpr auto    p = SparsityPattern<2, 3>{}.unit_diagonal(0, 0, 1).values(1, 1, 1, 2);
        PatternMatrix<p, double> m{};
        m.at<1, 1>() = 2.5;
        m.at<1, 2>() = -3.0;
        CHECK(m(0, 0) == 1.0);
        CHECK(m(0, 1) == 0.0);
        CHECK(m(1, 1) == 2.5);
        CHECK(m(1, 2) == -3.0);
        CHECK(m.to_matrix()(1, 2) == -3.0);
        static_assert(std::is_same_v<decltype(m.values), damp::array<double, 2>>);
    }

    TEST_CASE("set_block writes only the runtime values") {
        static constexpr auto    p = SparsityPattern<3, 3>::identity().values_off_diagonal(0, 0, 3);
        PatternMatrix<p, double> m{};
        m.set_block<0, 0>(Matrix<3, 3>{{9.0, 1.0, 2.0}, {3.0, 9.0, 4.0}, {5.0, 6.0, 9.0}});
        CHECK(m(0, 0) == 1.0); // the pattern's one, not the block's 9
        CHECK(m(0, 1) == 1.0);
        CHECK(m(2, 1) == 6.0);
    }

    TEST_CASE("its products are the dense ones") {
        PatternMatrix<kLever, double> h{};
        fill(h, 0.3);
        const Matrix<3, 7> d = h.to_matrix();
        const ColVec<7>    x{1.0, -2.0, 3.0, 0.5, -0.25, 7.0, 0.1};
        CHECK(gap(Matrix<3, 1>(h * x), Matrix<3, 1>(d * x)) <= 1e-13);
        const Matrix<7, 7> s = spd<7>(0.4);
        CHECK(gap(h * s, d * s) <= 1e-13);
        CHECK(gap(times_transpose(s, h), s * d.transpose()) <= 1e-13);
        CHECK(gap(PatternMatrix<kMissesState, double>{}.to_matrix(), Matrix<2, 5>{}) == 0.0);
    }

    TEST_CASE("an all-zero row of the pattern gives zero, in a product and a sandwich") {
        PatternMatrix<kZeroRow, double> f{};
        fill(f, 0.2);
        const ColVec<4> y = f * ColVec<4>{1.0, 2.0, 3.0, 4.0};
        CHECK(y[2] == 0.0);
        const Matrix<4, 4> out = quadratic_form(f, spd<4>(0.1));
        for (size_t j = 0; j < 4; ++j) {
            CHECK(out(2, j) == 0.0);
            CHECK(out(j, 2) == 0.0);
        }
    }

    TEST_CASE("F S Fᵀ is the dense one, exactly symmetric, in place or not, on every shape") {
        check_sandwich<kIns>(0.1);
        check_sandwich<kAllUnit>(0.2);
        check_sandwich<kNoUnit>(0.3);
        check_sandwich<kApartUnits>(0.4);
        check_sandwich<kZeroRow>(0.5);
    }

    TEST_CASE("the identity leaves S alone: no write, no rounding") {
        const Matrix<5, 5> s = spd<5>(0.6);
        CHECK(gap(quadratic_form(PatternMatrix<kAllUnit, double>{}, s), s) == 0.0);
    }

    TEST_CASE("it works in constant evaluation") {
        static constexpr auto kP = SparsityPattern<2, 2>::identity().value_diagonal(0, 1, 1);
        static constexpr auto out = [] {
            PatternMatrix<kP, double> f{};
            f.at<0, 1>() = 2.0;
            return quadratic_form(f, Matrix<2, 2>{{1.0, 0.5}, {0.5, 3.0}});
        }();
        static_assert(out(0, 0) == 1.0 + 2.0 + 12.0); // 1 + 2·2·0.5 + 4·3
        static_assert(out(0, 1) == 0.5 + 6.0 && out(1, 0) == 0.5 + 6.0);
        static_assert(out(1, 1) == 3.0);
    }
}

TEST_SUITE("ESKF on sparse Jacobians") {
    constexpr size_t N = 15;

    /// The same run on the sparse and the dense paths, comparing the whole state after every step.
    template<const auto& HPattern, size_t NY>
    void check_run(const Matrix<N, N>& Q) {
        const ErrorStateKalmanFilterConfig<N, double> config{.Q = Q, .P0 = spd<N>(0.7)};
        ErrorStateKalmanFilter<N, double>             sparse{config};
        ErrorStateKalmanFilter<N, double>             dense{config};
        double                                        worst = 0.0;
        double                                        scale = 1.0;
        for (int k = 0; k < 60; ++k) {
            const double                                 seed = 0.1 * static_cast<double>(k);
            PatternMatrix<kInsTransitionPattern, double> F{};
            fill(F, seed);
            sparse.predict(F);
            const Matrix<N, N> F_dense = F.to_matrix();
            dense.predict([&](double) { return ErrorStateJacobian<double, N>{F_dense, Matrix<N, N>::identity()}; }, 0.01);
            worst = std::fmax(worst, gap(sparse.covariance(), dense.covariance()));

            PatternMatrix<HPattern, double> H{};
            fill(H, seed + 1.0);
            const Matrix<NY, NY> R = Matrix<NY, NY>::identity() * 0.5;
            ColVec<NY>           innov{};
            for (size_t i = 0; i < NY; ++i) {
                innov[i] = std::cos(seed + static_cast<double>(i));
            }
            const Matrix<NY, N> H_dense = H.to_matrix();
            REQUIRE(sparse.update_innovation(innov, H, R));
            REQUIRE(dense.update_innovation(innov, H_dense, R));
            worst = std::fmax(worst, gap(sparse.covariance(), dense.covariance()));
            worst = std::fmax(worst, gap(Matrix<N, 1>(sparse.error_state()), Matrix<N, 1>(dense.error_state())));

            PatternMatrix<kInsResetPattern, double> G{};
            fill(G, seed + 2.0);
            sparse.reset_error_state(G);
            dense.reset_error_state(G.to_matrix());
            worst = std::fmax(worst, gap(sparse.covariance(), dense.covariance()));
            scale = std::fmax(scale, dense.covariance().norm());

            const auto& P = sparse.covariance();
            for (size_t i = 0; i < N; ++i) {
                for (size_t j = 0; j < N; ++j) {
                    CHECK(P(i, j) == P(j, i)); // exactly
                }
            }
        }
        CAPTURE(worst);
        CAPTURE(scale);
        CHECK(worst <= 1e-9 * scale);
    }

    TEST_CASE("predict, update and reset are the dense filter's, on selector and runtime measurement shapes") {
        const Matrix<N, N> Q = Matrix<N, N>::identity() * 1e-3;
        check_run<kInsPositionPattern, 3>(Q);
        check_run<kInsVelocityPattern, 3>(Q);
        check_run<kInsOrientationPattern, 3>(Q);
        check_run<kInsHeadingPattern, 1>(Q);
    }

    TEST_CASE("a Q with off-diagonal entries is added whole, and symmetrized") {
        Matrix<N, N> Q = Matrix<N, N>::identity() * 1e-3;
        Q(3, 7) = 2e-4;
        Q(7, 3) = 2e-4;
        check_run<kInsPositionPattern, 3>(Q);

        Matrix<N, N> lopsided = Matrix<N, N>::identity() * 1e-3;
        lopsided(1, 2) = 4e-4; // stored asymmetry: the dense path scrubs it, and so does this one
        const ErrorStateKalmanFilterConfig<N, double> config{.Q = lopsided, .P0 = Matrix<N, N>::identity()};
        ErrorStateKalmanFilter<N, double>             f{config};
        f.predict(PatternMatrix<kInsResetPattern, double>{});
        CHECK(f.covariance()(1, 2) == doctest::Approx(2e-4));
        CHECK(f.covariance()(2, 1) == doctest::Approx(2e-4));
    }

    TEST_CASE("a measurement that cannot be solved leaves the filter as it was") {
        static constexpr auto                   kP = SparsityPattern<1, N>{}.unit_diagonal(0, 0, 1);
        ErrorStateKalmanFilterConfig<N, double> config{.Q = Matrix<N, N>{}, .P0 = Matrix<N, N>{}};
        ErrorStateKalmanFilter<N, double>       f{config};
        const ColVec<N>                         dx0{};
        const Matrix<N, N>                      before = f.covariance();
        CHECK_FALSE(f.update_innovation(ColVec<1>{1.0}, PatternMatrix<kP, double>{}, Matrix<1, 1>{})); // S = 0
        CHECK(gap(f.covariance(), before) == 0.0);
        CHECK(gap(Matrix<N, 1>(f.error_state()), Matrix<N, 1>(dx0)) == 0.0);
    }

    TEST_CASE("the δx not yet injected is carried through F and the update") {
        ErrorStateKalmanFilter<N, double> f{ErrorStateKalmanFilterConfig<N, double>{.Q = Matrix<N, N>{}, .P0 = Matrix<N, N>::identity()}};
        ColVec<N>                         dx{};
        dx[ins_err::dv] = 2.0;
        f.set_error_state(dx);
        PatternMatrix<kInsTransitionPattern, double> F{};
        F.at<ins_err::dp, ins_err::dv>() = 0.5;
        f.predict(F);
        CHECK(f.error_state()[ins_err::dp] == doctest::Approx(1.0));
        CHECK(f.error_state()[ins_err::dv] == doctest::Approx(2.0));
    }
}

TEST_SUITE("The INS Jacobians") {
    TEST_CASE("the transition's pattern holds every entry the textbook F has") {
        // The dense F, written out independently: anything the pattern dropped shows up here.
        const Quaternion<double> q = Quaternion<double>::from_axis_angle(Vec3<double>{1.0, 2.0, -1.0}.normalized(), 0.7).value();
        const Vec3<double>       omega{0.3, -0.2, 0.5};
        const Vec3<double>       a{1.1, -0.4, 9.7};
        constexpr double         dt = 0.01;

        const Matrix<3, 3> R = q.to_dcm().matrix();
        const Matrix<3, 3> W{{0.0, -omega[2], omega[1]}, {omega[2], 0.0, -omega[0]}, {-omega[1], omega[0], 0.0}};
        const Matrix<3, 3> A{{0.0, -a[2], a[1]}, {a[2], 0.0, -a[0]}, {-a[1], a[0], 0.0}};
        Matrix<15, 15>     dense = Matrix<15, 15>::identity();
        dense.block<3, 3>(ins_err::dtheta, ins_err::dtheta) = Matrix<3, 3>::identity() - (W * dt);
        dense.block<3, 3>(ins_err::dtheta, ins_err::dbg) = Matrix<3, 3>::identity() * (-dt);
        dense.block<3, 3>(ins_err::dv, ins_err::dtheta) = (R * A) * (-dt);
        dense.block<3, 3>(ins_err::dv, ins_err::dba) = R * (-dt);
        dense.block<3, 3>(ins_err::dp, ins_err::dv) = Matrix<3, 3>::identity() * dt;

        const auto F = ins_error_jacobian(q, omega, a, dt);
        CHECK(gap(F.to_matrix(), dense) <= 1e-14);
    }
}
