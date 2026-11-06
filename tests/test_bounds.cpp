// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <limits>

#include "damp/matrix/matrix.hpp"
#include "damp/toolbox/bounds.hpp"
#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

TEST_CASE("Bounds: default is unbounded (not a zero box)") {
    constexpr Bounds<2> b{};
    // The trap this guards: a zero-initialised box would clamp everything to 0.
    CHECK(b.saturate(ColVec<2>{1e9, -1e9}) == ColVec<2>{1e9, -1e9});
    CHECK(b.contains(ColVec<2>{1e9, -1e9}));
}

TEST_CASE("Bounds: per-channel saturate and contains") {
    const Bounds<3> b{
        .lower = {-1.0, 0.0, -5.0},
        .upper = {1.0, 2.0, 5.0},
    };
    CHECK(b.saturate(ColVec<3>{-2.0, 3.0, 0.0}) == ColVec<3>{-1.0, 2.0, 0.0});
    CHECK(b.contains(ColVec<3>{0.0, 1.0, 4.9}));
    CHECK_FALSE(b.contains(ColVec<3>{0.0, 1.0, 5.1}));
}

TEST_CASE("Bounds::symmetric is the |x|<=mag convention") {
    const auto b = Bounds<2>::symmetric({3.0, 4.0});
    CHECK(b.lower == ColVec<2>{-3.0, -4.0});
    CHECK(b.upper == ColVec<2>{3.0, 4.0});
    CHECK(b.saturate(ColVec<2>{-9.0, 9.0}) == ColVec<2>{-3.0, 4.0});
}

TEST_CASE("Bounds<1>: SISO scalar conveniences") {
    constexpr Bounds<1> b{
        .lower = -2.0,
        .upper = 2.0,
    };
    CHECK(b.saturate(5.0) == 2.0);
    CHECK(b.saturate(-5.0) == -2.0);
    CHECK(b.saturate(1.0) == 1.0);
    CHECK(b.contains(1.5));
    CHECK_FALSE(b.contains(2.5));
}

TEST_CASE("Bounds::as maps the unbounded ends to the target's max and casts the finite ones") {
    constexpr Bounds<2> b{
        .lower = {-std::numeric_limits<double>::max(), -0.5},
        .upper = {0.25, std::numeric_limits<double>::max()},
    };
    constexpr Bounds<2, float> f = b.as<float>();
    static_assert(f.lower(0) == -std::numeric_limits<float>::max());
    static_assert(f.lower(1) == -0.5f);
    static_assert(f.upper(0) == 0.25f);
    static_assert(f.upper(1) == std::numeric_limits<float>::max());
    constexpr Bounds<2> back = f.as<double>();
    static_assert(back.lower(0) == -std::numeric_limits<double>::max());
    static_assert(back.upper(1) == std::numeric_limits<double>::max());
    CHECK(back.upper(0) == 0.25);
}
