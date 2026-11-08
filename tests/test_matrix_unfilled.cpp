// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// A Matrix about to be overwritten is not zero-filled at run time. Every entry still has to be what
// the plain arithmetic says, with junk left on the stack where the fill used to be, and constant
// evaluation still starts from zeros.

#include <cmath>
#include <cstddef>

#include "damp/math/complex.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/rowvec.hpp"
#include "doctest.h"

using namespace damp;

namespace {

constexpr size_t kBig = 16; // past for_each_entry's 64 entries, so loops do the writing

/// Leave non-zero junk where the next call's frame goes, so a read of an entry never written shows.
[[gnu::noinline]] float dirty_stack() {
    volatile float junk[16384];
    for (size_t i = 0; i < 16384; ++i) {
        junk[i] = 1.25e30f;
    }
    return junk[16383];
}

/// @p op, a noinline lambda holding its results in its own frame, on a dirty stack.
template<typename Op>
[[gnu::noinline]] size_t on_dirty_stack(const Op& op) {
    static_cast<void>(dirty_stack());
    return op();
}

template<typename T>
T wiggle(size_t i, size_t j, double seed) {
    return static_cast<T>(std::sin((1.3 * static_cast<double>(i)) + (0.7 * static_cast<double>(j)) + seed));
}

template<size_t R, size_t C, typename T>
Matrix<R, C, T> wiggly(double seed) {
    Matrix<R, C, T> m;
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            m(i, j) = wiggle<T>(i, j, seed);
        }
    }
    return m;
}

/// How many entries of @p got differ from @p want(i, j) by more than rounding.
template<typename Got, typename Want>
size_t wrong_entries(const Got& got, const Want& want) {
    size_t wrong = 0;
    for (size_t i = 0; i < Got::rows(); ++i) {
        for (size_t j = 0; j < Got::cols(); ++j) {
            const auto expected = static_cast<double>(want(i, j));
            if (!(std::abs(static_cast<double>(got(i, j)) - expected) <= 1e-3 * (1.0 + std::abs(expected)))) {
                ++wrong;
            }
        }
    }
    return wrong;
}

template<size_t R, size_t K, size_t C, typename T, typename RhsAt>
auto product_of(const Matrix<R, K, T>& a, const RhsAt& rhs_at) {
    return [&a, rhs_at](size_t i, size_t j) {
        double sum = 0.0;
        for (size_t k = 0; k < K; ++k) {
            sum += static_cast<double>(a(i, k)) * static_cast<double>(rhs_at(k, j));
        }
        return sum;
    };
}

} // namespace

TEST_SUITE("Matrix without a zero-fill") {
    TEST_CASE_TEMPLATE("products are the plain arithmetic on a dirty stack", T, float, double) {
        const auto a = wiggly<kBig, kBig, T>(0.1);
        const auto b = wiggly<kBig, kBig, T>(0.9);
        const auto c = wiggly<kBig, 5, T>(1.7);
        const auto b_at = [&](size_t i, size_t j) { return b(i, j); };
        const auto bt_at = [&](size_t i, size_t j) { return b(j, i); };
        const auto c_at = [&](size_t i, size_t j) { return c(i, j); };
        const auto ab = product_of<kBig, kBig, kBig, T>(a, b_at);
        const auto abt = product_of<kBig, kBig, kBig, T>(a, bt_at);
        const auto ac = product_of<kBig, kBig, 5, T>(a, c_at);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(a * b, ab); }) == 0);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(a * b.t(), abt); }) == 0);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(a * c, ac); }) == 0);
    }

    TEST_CASE_TEMPLATE("entry-wise results are the plain arithmetic on a dirty stack", T, float, double) {
        const auto a = wiggly<kBig, kBig, T>(0.1);
        const auto b = wiggly<kBig, kBig, T>(0.9);
        const T    s = static_cast<T>(2.5);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(a + b, [&](size_t i, size_t j) { return a(i, j) + b(i, j); }); }) == 0);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(a - b, [&](size_t i, size_t j) { return a(i, j) - b(i, j); }); }) == 0);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(-a, [&](size_t i, size_t j) { return -a(i, j); }); }) == 0);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(a + s, [&](size_t i, size_t j) { return a(i, j) + s; }); }) == 0);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(a - s, [&](size_t i, size_t j) { return a(i, j) - s; }); }) == 0);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(s - a, [&](size_t i, size_t j) { return s - a(i, j); }); }) == 0);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(a * s, [&](size_t i, size_t j) { return a(i, j) * s; }); }) == 0);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(a / s, [&](size_t i, size_t j) { return a(i, j) / s; }); }) == 0);
    }

    TEST_CASE_TEMPLATE("copies and rearrangements are whole on a dirty stack", T, float, double) {
        const auto a = wiggly<kBig, kBig, T>(0.1);
        const auto a_at = [&](size_t i, size_t j) { return a(i, j); };
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) {
            const Matrix<kBig, kBig, T> copy{a};
            return wrong_entries(copy, a_at);
        }) == 0);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(a.transpose(), [&](size_t i, size_t j) { return a(j, i); }); }) == 0);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(a.conjugate_transpose(), [&](size_t i, size_t j) { return a(j, i); }); }) == 0);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(a.template reshape<8, 32>(), [&](size_t i, size_t j) { const size_t f = (i * 32) + j; return a(f / kBig, f % kBig); }); }) == 0);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(a.row_vector(3), [&](size_t, size_t j) { return a(3, j); }); }) == 0);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(a.col_vector(3), [&](size_t i, size_t) { return a(i, 3); }); }) == 0);
    }

    TEST_CASE_TEMPLATE("conversions and small results are whole on a dirty stack", T, float, double) {
        //! A conversion fills through elements_of, and a result of at most 64 entries through
        //! for_each_entry's constant-index fold rather than a loop.
        const auto big = wiggly<kBig, kBig, T>(0.1);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(big.template as<double>(), [&](size_t i, size_t j) { return big(i, j); }); }) == 0);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(Matrix<kBig, kBig, T>(big.t()), [&](size_t i, size_t j) { return big(j, i); }); }) == 0);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(Matrix<4, 4, T>(big.template block<4, 4>(2, 3)), [&](size_t i, size_t j) { return big(2 + i, 3 + j); }); }) == 0);
        const auto a3 = wiggly<3, 3, T>(0.2);
        const auto b3 = wiggly<3, 3, T>(0.8);
        const auto a4 = wiggly<4, 4, T>(0.2);
        const auto b4 = wiggly<4, 4, T>(0.8);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(a3 + b3, [&](size_t i, size_t j) { return a3(i, j) + b3(i, j); }); }) == 0);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(-a4, [&](size_t i, size_t j) { return -a4(i, j); }); }) == 0);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(a4 - b4, [&](size_t i, size_t j) { return a4(i, j) - b4(i, j); }); }) == 0);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(a4 * b4, product_of<4, 4, 4, T>(a4, [&](size_t i, size_t j) { return b4(i, j); })); }) == 0);
    }

    TEST_CASE("a product summing to -0 is +0, as a zero-filled sum gives") {
        const Matrix<1, 2> a{{-0.0, 0.0}};
        const Matrix<2, 1> b{{1.0}, {-1.0}};
        CHECK_FALSE(std::signbit((a * b)(0, 0)));
    }

    TEST_CASE_TEMPLATE("vector products and congruences are the plain arithmetic on a dirty stack", T, float, double) {
        const auto a = wiggly<kBig, kBig, T>(0.1);
        const auto m = wiggly<6, kBig, T>(1.3);
        const auto x = ColVec<kBig, T>(wiggly<kBig, 1, T>(0.4));
        const auto y = RowVec<kBig, T>(wiggly<1, kBig, T>(0.4));
        const auto x_at = [&](size_t k, size_t) { return x(k); };
        const auto ax = product_of<kBig, kBig, 1, T>(a, x_at);
        const auto ya = [&](size_t, size_t j) {
            double sum = 0.0;
            for (size_t k = 0; k < kBig; ++k) {
                sum += static_cast<double>(y(k)) * static_cast<double>(a(k, j));
            }
            return sum;
        };
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(a * x, ax); }) == 0);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(y * a, ya); }) == 0);

        Matrix<kBig, kBig, T> sym{};
        for (size_t i = 0; i < kBig; ++i) {
            for (size_t j = 0; j < kBig; ++j) {
                sym(i, j) = a(i, j) + a(j, i);
            }
        }
        const auto want = [&](size_t i, size_t j) {
            double sum = 0.0;
            for (size_t k = 0; k < kBig; ++k) {
                for (size_t l = 0; l < kBig; ++l) {
                    sum += static_cast<double>(m(i, k)) * static_cast<double>(sym(k, l)) * static_cast<double>(m(j, l));
                }
            }
            return sum;
        };
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) { return wrong_entries(quadratic_form(m, sym), want); }) == 0);
    }

    TEST_CASE_TEMPLATE("triangular solves satisfy their systems on a dirty stack", T, float, double) {
        Matrix<kBig, kBig, T> spd = wiggly<kBig, kBig, T>(0.1) * wiggly<kBig, kBig, T>(0.1).transpose();
        for (size_t i = 0; i < kBig; ++i) {
            spd(i, i) += static_cast<T>(kBig);
        }
        const auto b = wiggly<kBig, 5, T>(0.6);
        const auto l = mat::cholesky(spd);
        REQUIRE(l.has_value());
        const auto u = l->transpose();
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) {
            const auto x = mat::solve(l->lower_triangle(), b);
            return x.has_value() ? wrong_entries((*l) * (*x), [&](size_t i, size_t j) { return b(i, j); }) : size_t{1};
        }) == 0);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) {
            const auto x = mat::solve(u.upper_triangle(), b);
            return x.has_value() ? wrong_entries(u * (*x), [&](size_t i, size_t j) { return b(i, j); }) : size_t{1};
        }) == 0);
        CHECK(on_dirty_stack([&]() __attribute__((noinline)) {
            const auto x = mat::cholesky_solve(spd, b);
            return x.has_value() ? wrong_entries(spd * (*x), [&](size_t i, size_t j) { return b(i, j); }) : size_t{1};
        }) == 0);
    }

    TEST_CASE("a triangular solve refuses a zero anywhere on the diagonal") {
        Matrix<3, 3> factor = Matrix<3, 3>::identity();
        factor(1, 1) = 0.0;
        const Matrix<3, 2> rhs{{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
        CHECK_FALSE(mat::solve(factor.lower_triangle(), rhs).has_value());
        CHECK_FALSE(mat::solve(factor.upper_triangle(), rhs).has_value());
    }

    TEST_CASE("complex products take the same path") {
        using C = complex<double>;
        Matrix<2, 2, C> a{};
        a(0, 0) = C{1.0, 2.0};
        a(0, 1) = C{0.0, -1.0};
        a(1, 0) = C{3.0, 0.0};
        a(1, 1) = C{-1.0, 1.0};
        const Matrix<2, 2, C> p = a * a;
        CHECK(p(0, 0).real() == doctest::Approx(-3.0 + 0.0));
        CHECK(p(0, 0).imag() == doctest::Approx(4.0 - 3.0));
        CHECK(p(1, 1).real() == doctest::Approx(0.0 + 0.0));
        CHECK(p(1, 1).imag() == doctest::Approx(-3.0 - 2.0));
    }

    TEST_CASE("an empty sum is a zero matrix, on every entry") {
        const Matrix<3, 0> a{};
        const Matrix<0, 2> b{};
        const Matrix<3, 2> p = a * b;
        CHECK(p == Matrix<3, 2>{});
    }

    TEST_CASE("constant evaluation still starts from zeros") {
        static_assert(Matrix<2, 2>{detail::unfilled}(1, 1) == 0.0);
        static constexpr Matrix<2, 2> a{{1.0, 0.0}, {0.0, 0.0}};
        static constexpr Matrix<2, 2> b{{5.0, 6.0}, {7.0, 8.0}};
        static constexpr auto         p = a * b;
        static_assert(p(0, 0) == 5.0 && p(0, 1) == 6.0 && p(1, 0) == 0.0 && p(1, 1) == 0.0);
        static constexpr auto q = a * b.t();
        static_assert(q(0, 0) == 5.0 && q(0, 1) == 7.0 && q(1, 0) == 0.0 && q(1, 1) == 0.0);
        static constexpr auto copy = [] {
            Matrix<2, 2> c{b};
            return c;
        }();
        static_assert(copy == b);
        static constexpr auto solved = mat::cholesky_solve(Matrix<2, 2>{{4.0, 0.0}, {0.0, 9.0}}, Matrix<2, 1>{{8.0}, {27.0}});
        static_assert(solved.has_value() && (*solved)(0, 0) == 2.0 && (*solved)(1, 0) == 3.0);
        CHECK(p(0, 0) == 5.0);
    }
}
