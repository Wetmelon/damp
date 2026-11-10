// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// A NaN or ±Inf that reaches a factorization or a solve is refused, under the runner's
// -ffast-math as in IEEE arithmetic. -ffinite-math-only lets GCC fold a negated test
// !(x > floor) into x ≤ floor, which a NaN fails, so a NaN pivot passed a test the IEEE build
// refuses it with. The values here arrive at run time (through a volatile), as a sensor's
// would: a NaN constant can be folded at compile time.

#include <bit>
#include <cstddef>
#include <cstdint>
#include <utility>

#include "damp/estimation/kalman.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

/// A quiet NaN made at run time.
template<typename T>
T runtime_nan() {
    volatile T zero = T{0};
    return zero / zero;
}

/// +Inf made at run time.
template<typename T>
T runtime_inf() {
    volatile T zero = T{0};
    return T{1} / zero;
}

/// The three non-finite values, made at run time; index 0 … 2.
template<typename T>
T non_finite(int k) {
    if (k == 0) {
        return runtime_nan<T>();
    }
    return (k == 1) ? runtime_inf<T>() : -runtime_inf<T>();
}

/// An entry whose double overflows: 3·10³⁸ in float, 1.5·10³⁰⁸ in double.
template<typename T>
T near_overflow() {
    if constexpr (sizeof(T) == 4) {
        return 3e38f;
    } else {
        return 1.5e308;
    }
}

/// The bits of a float or double.
template<typename T>
auto bits_of(T x) {
    if constexpr (sizeof(T) == 4) {
        return std::bit_cast<std::uint32_t>(x);
    } else {
        return std::bit_cast<std::uint64_t>(x);
    }
}

/// Every entry of @p a and @p b has the same bits.
template<size_t R, size_t C, typename T>
bool same_bits(const Matrix<R, C, T>& a, const Matrix<R, C, T>& b) {
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            if (bits_of(a(i, j)) != bits_of(b(i, j))) {
                return false;
            }
        }
    }
    return true;
}

/// A symmetric positive-definite 3 × 3 matrix: every pivot well above the floor.
template<typename T>
Matrix<3, 3, T> spd3() {
    return Matrix<3, 3, T>{{T(4), T(1), T(0.5)}, {T(1), T(3), T(0.25)}, {T(0.5), T(0.25), T(2)}};
}

} // namespace

TEST_SUITE("non-finite input to the factorizations") {
    TEST_CASE_TEMPLATE("cholesky and cholesky_solve refuse a NaN or ±Inf in each lower-triangle position", T, float, double) {
        // A NaN off the diagonal reaches pivot i through lᵢⱼ: under -ffast-math the NaN pivot
        // passed the floor test and cholesky returned a NaN factor, cholesky_solve a NaN x.
        const Matrix<3, 1, T> b{{T(1)}, {T(2)}, {T(3)}};
        REQUIRE(mat::cholesky(spd3<T>()).has_value());
        for (int k = 0; k < 3; ++k) {
            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = 0; j <= i; ++j) {
                    INFO("value ", k, " at (", i, ", ", j, ")");
                    Matrix<3, 3, T> lower = spd3<T>(); // the upper triangle is never read
                    lower(i, j) = non_finite<T>(k);
                    CHECK_FALSE(mat::cholesky(lower).has_value());
                    CHECK_FALSE(mat::cholesky_solve(lower, b).has_value());
                    Matrix<3, 3, T> sym = lower;
                    sym(j, i) = lower(i, j);
                    CHECK_FALSE(mat::cholesky(sym).has_value());
                    CHECK_FALSE(mat::cholesky_solve(sym, b).has_value());
                    CHECK_FALSE(mat::solve(sym, b).has_value());
                }
            }
        }
    }

    TEST_CASE_TEMPLATE("complex cholesky refuses a NaN or ±Inf in either part of a lower-triangle entry", T, float, double) {
        using C = damp::complex<T>;
        const Matrix<2, 2, C> A{{C(T(4), T(0)), C(T(1), T(-1))}, {C(T(1), T(1)), C(T(3), T(0))}};
        REQUIRE(mat::cholesky(A).has_value());
        for (int k = 0; k < 3; ++k) {
            for (size_t i = 0; i < 2; ++i) {
                for (size_t j = 0; j <= i; ++j) {
                    for (const bool imag : {false, true}) {
                        INFO("value ", k, " at (", i, ", ", j, "), imaginary part ", imag);
                        Matrix<2, 2, C> bad = A;
                        bad(i, j) = imag ? C(A(i, j).real(), non_finite<T>(k)) : C(non_finite<T>(k), A(i, j).imag());
                        CHECK_FALSE(mat::cholesky(bad).has_value());
                    }
                }
            }
        }
    }

    TEST_CASE_TEMPLATE("lu_decomposition refuses a pivot that a finite A overflows to NaN", T, float, double) {
        // Row 1 reaches −Inf at step 0 and Inf − Inf = NaN at step 1. At -O1 and -O2 the
        // roundoff and floor tests let the NaN pivot through and solve returned NaN.
        volatile T            big_v = near_overflow<T>();
        const T               big = big_v;
        const Matrix<3, 3, T> A{{T(1), T(0), big}, {T(1), T(1), -big}, {T(1), T(2), -big}};
        const Matrix<3, 1, T> b{{T(1)}, {T(1)}, {T(1)}};
        CHECK_FALSE(mat::lu_decomposition(A).has_value());
        CHECK_FALSE(mat::lu_solve(A, b).has_value());
        CHECK_FALSE(mat::solve(A, b).has_value());
        CHECK_FALSE(A.inverse().has_value());
    }

    TEST_CASE_TEMPLATE("triangular solves refuse a NaN or ±Inf on the diagonal", T, float, double) {
        const Matrix<2, 1, T> b{{T(1)}, {T(1)}};
        for (int k = 0; k < 3; ++k) {
            for (size_t i = 0; i < 2; ++i) {
                INFO("value ", k, " at (", i, ", ", i, ")");
                Matrix<2, 2, T> L{{T(2), T(0)}, {T(1), T(3)}};
                L(i, i) = non_finite<T>(k);
                CHECK_FALSE(mat::solve(L.lower_triangle(), b).has_value());
                Matrix<2, 2, T> U{{T(2), T(1)}, {T(0), T(3)}};
                U(i, i) = non_finite<T>(k);
                CHECK_FALSE(mat::solve(U.upper_triangle(), b).has_value());
            }
        }
    }
}

TEST_SUITE("KalmanFilter refuses a non-finite reading covariance") {
    TEST_CASE_TEMPLATE("update(y, C, R) with a NaN or Inf in R returns false and leaves x and P bit for bit", T, float, double) {
        // S = C P Cᵀ + R inherits the NaN. Cholesky factored the NaN pivot under -ffast-math, and
        // update returned true and wrote NaN into x and P.
        StateSpace<2, 1, 1, T, 2, 1> sys{};
        sys.A = Matrix<2, 2, T>::identity();
        sys.B = {{T(0)}, {T(1)}};
        sys.C = {{T(1), T(0)}};
        sys.G = Matrix<2, 2, T>::identity();
        sys.H = Matrix<1, 1, T>::identity();
        sys.Ts = T(1);
        for (int k = 0; k < 2; ++k) {
            for (const auto& at : {std::pair<size_t, size_t>{0, 0}, {1, 1}, {1, 0}}) {
                const size_t i = at.first;
                const size_t j = at.second;
                INFO("value ", k, " at (", i, ", ", j, ")");
                KalmanFilter kf{sys, {.Q = Matrix<2, 2, T>::identity() * T(0.25), .R = T(1)}};
                kf.predict(ColVec<1, T>{T(0.5)});
                REQUIRE(kf.update(ColVec<1, T>{T(0.75)})); // a finite x and a full P to keep
                const ColVec<2, T>    x0 = kf.state();
                const Matrix<2, 2, T> P0 = kf.covariance();
                Matrix<2, 2, T>       R = Matrix<2, 2, T>::diagonal({T(1), T(3)});
                R(i, j) = non_finite<T>(k);
                R(j, i) = R(i, j);
                CHECK_FALSE(kf.update(ColVec<2, T>{T(2), T(4)}, Matrix<2, 2, T>::identity(), R));
                CHECK(same_bits(kf.state(), x0));
                CHECK(same_bits(kf.covariance(), P0));
            }
        }
    }
}
