// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "damp/backend.hpp"
#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

/**
 * @brief Tests for the backend-agnostic min/max/minmax value helpers.
 *
 * Guards the initializer_list overloads in particular: ETL's min/max/minmax are
 * binary-only and dangle (return references), so damp supplies by-value versions.
 * These must behave identically on both backends.
 */
TEST_SUITE("backend min/max/minmax") {
    TEST_CASE("minmax(a, b) orders the pair and is stable on ties") {
        CHECK(damp::minmax(3, 1) == pair{1, 3});
        CHECK(damp::minmax(1, 3) == pair{1, 3});
        CHECK(damp::minmax(2, 2) == pair{2, 2});
    }

    TEST_CASE("minmax(initializer_list) finds the extremes") {
        const auto [lo, hi] = damp::minmax({3.0, 1.0, 2.0, 5.0, -1.0});
        CHECK(lo == doctest::Approx(-1.0));
        CHECK(hi == doctest::Approx(5.0));
        // constexpr-usable
        static_assert(damp::minmax({3, 1, 2}).first == 1);
        static_assert(damp::minmax({3, 1, 2}).second == 3);
    }

    TEST_CASE("min/max(initializer_list) match the std/etl-parity contract") {
        CHECK(damp::min({4, 2, 7, 1, 9}) == 1);
        CHECK(damp::max({4, 2, 7, 1, 9}) == 9);
        static_assert(damp::min({4, 2, 7}) == 2);
        static_assert(damp::max({4, 2, 7}) == 7);
    }

    TEST_CASE("single-element lists return that element") {
        CHECK(damp::min({42}) == 42);
        CHECK(damp::max({42}) == 42);
        CHECK(damp::minmax({42}) == pair{42, 42});
    }
}
