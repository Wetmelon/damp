// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file test_quantize.cpp
 * @brief ADC quantizer for SIL readings (sim::quantize)
 *
 * Exact checks use a dyadic lsb (a power of two), so x/lsb and lsb·n are exact in
 * float and double; the non-dyadic case is checked against a tolerance scaled to lsb.
 */

#include "damp/backend.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "damp/simulation/quantize.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using namespace damp::sim;

// Constant-evaluated: a tie at 1.5 codes goes to the even code 2.
static_assert(quantize(0.375, 0.25, -1.0, 1.0) == 0.5);
static_assert(quantize(0.375F, 0.25F, -1.0F, 1.0F) == 0.5F);

TEST_SUITE("ADC quantizer") {

    TEST_CASE_TEMPLATE("rounds to the nearest multiple of lsb", T, float, double) {
        const T lsb = T{0.25};
        const T lo = T{-8};
        const T hi = T{8};
        CHECK(quantize(T{0}, lsb, lo, hi) == T{0});
        CHECK(quantize(static_cast<T>(0.1), lsb, lo, hi) == T{0});
        CHECK(quantize(static_cast<T>(0.13), lsb, lo, hi) == T{0.25});
        CHECK(quantize(static_cast<T>(0.9), lsb, lo, hi) == T{1});
        CHECK(quantize(static_cast<T>(1.12), lsb, lo, hi) == T{1});
        CHECK(quantize(T{3.25}, lsb, lo, hi) == T{3.25}); // on the grid: unchanged
    }

    TEST_CASE_TEMPLATE("negative values round symmetrically", T, float, double) {
        const T lsb = T{0.25};
        const T lo = T{-8};
        const T hi = T{8};
        CHECK(quantize(static_cast<T>(-0.1), lsb, lo, hi) == T{0});
        CHECK(quantize(static_cast<T>(-0.13), lsb, lo, hi) == T{-0.25});
        CHECK(quantize(static_cast<T>(-0.9), lsb, lo, hi) == T{-1});
        CHECK(quantize(static_cast<T>(-1.12), lsb, lo, hi) == T{-1});
        for (const T x : {static_cast<T>(0.13), static_cast<T>(0.62), static_cast<T>(2.4), T{0.375}}) {
            CHECK(quantize(-x, lsb, lo, hi) == -quantize(x, lsb, lo, hi));
        }
    }

    TEST_CASE_TEMPLATE("ties round to the even code", T, float, double) {
        const T lsb = T{0.25};
        const T lo = T{-8};
        const T hi = T{8};
        CHECK(quantize(T{0.125}, lsb, lo, hi) == T{0});     // 0.5 codes → 0
        CHECK(quantize(T{0.375}, lsb, lo, hi) == T{0.5});   // 1.5 codes → 2
        CHECK(quantize(T{0.625}, lsb, lo, hi) == T{0.5});   // 2.5 codes → 2
        CHECK(quantize(T{0.875}, lsb, lo, hi) == T{1});     // 3.5 codes → 4
        CHECK(quantize(T{-0.125}, lsb, lo, hi) == T{0});    // −0.5 codes → 0
        CHECK(quantize(T{-0.375}, lsb, lo, hi) == T{-0.5}); // −1.5 codes → −2
        CHECK(quantize(T{-0.625}, lsb, lo, hi) == T{-0.5}); // −2.5 codes → −2
    }

    TEST_CASE_TEMPLATE("clamps at both rails of a 12-bit bidirectional converter", T, float, double) {
        // ±8 full scale: lsb = 16 / 4096 = 2⁻⁸, codes −2048 … 2047.
        const T lsb = T{16} / T{4096};
        const T lo = T{-2048} * lsb;
        const T hi = T{2047} * lsb;
        CHECK(lo == T{-8});
        CHECK(hi == T{7.99609375});

        CHECK(quantize(T{100}, lsb, lo, hi) == hi);
        CHECK(quantize(T{-100}, lsb, lo, hi) == lo);
        CHECK(quantize(static_cast<T>(7.999), lsb, lo, hi) == hi); // rounds to 8, above the top code
        CHECK(quantize(static_cast<T>(-8.001), lsb, lo, hi) == lo);
        CHECK(quantize(hi, lsb, lo, hi) == hi);
        CHECK(quantize(lo, lsb, lo, hi) == lo);
        CHECK(quantize(hi - lsb, lsb, lo, hi) == hi - lsb);
    }

    TEST_CASE_TEMPLATE("non-dyadic lsb lands on the code grid", T, float, double) {
        const T lsb = static_cast<T>(0.1);
        const T lo = T{-10};
        const T hi = T{10};
        const T tol = T{8} * damp::eps<T>() * T{10}; // a few ulps of the largest reading
        CHECK(damp::abs(quantize(static_cast<T>(0.26), lsb, lo, hi) - static_cast<T>(0.3)) <= tol);
        CHECK(damp::abs(quantize(static_cast<T>(-4.04), lsb, lo, hi) - static_cast<T>(-4.0)) <= tol);
        for (const T x : {static_cast<T>(0.26), static_cast<T>(-4.04), static_cast<T>(7.77)}) {
            const T codes = quantize(x, lsb, lo, hi) / lsb;
            CHECK(damp::abs(codes - damp::nearbyint(codes)) <= T{8} * damp::eps<T>() * T{100});
            CHECK(damp::abs(quantize(x, lsb, lo, hi) - x) <= (lsb / T{2}) + tol);
        }
    }

    TEST_CASE_TEMPLATE("lsb ≤ 0 clamps without rounding", T, float, double) {
        CHECK(quantize(static_cast<T>(0.13), T{0}, T{-1}, T{1}) == static_cast<T>(0.13));
        CHECK(quantize(static_cast<T>(0.13), T{-1}, T{-1}, T{1}) == static_cast<T>(0.13));
        CHECK(quantize(T{3}, T{0}, T{-1}, T{1}) == T{1});
        CHECK(quantize(T{-3}, T{0}, T{-1}, T{1}) == T{-1});
    }
}
