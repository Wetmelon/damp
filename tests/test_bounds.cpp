// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "damp/backend.hpp"
#include "damp/toolbox/bounds.hpp"
#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

TEST_CASE("Bounds: default is unbounded (not a zero box)") {
    constexpr Bounds<2> b{};
    // The trap this guards: a zero-initialised box would clamp everything to 0.
    CHECK(b.saturate(damp::array<double, 2>{1e9, -1e9}) == damp::array<double, 2>{1e9, -1e9});
    CHECK(b.contains(damp::array<double, 2>{1e9, -1e9}));
}

TEST_CASE("Bounds: per-channel saturate and contains") {
    const Bounds<3> b{{-1.0, 0.0, -5.0}, {1.0, 2.0, 5.0}};
    CHECK(b.saturate(damp::array<double, 3>{-2.0, 3.0, 0.0}) == damp::array<double, 3>{-1.0, 2.0, 0.0});
    CHECK(b.contains(damp::array<double, 3>{0.0, 1.0, 4.9}));
    CHECK_FALSE(b.contains(damp::array<double, 3>{0.0, 1.0, 5.1}));
}

TEST_CASE("Bounds::symmetric is the |x|<=mag convention") {
    const auto b = Bounds<2>::symmetric({3.0, 4.0});
    CHECK(b.lower == damp::array<double, 2>{-3.0, -4.0});
    CHECK(b.upper == damp::array<double, 2>{3.0, 4.0});
    CHECK(b.saturate(damp::array<double, 2>{-9.0, 9.0}) == damp::array<double, 2>{-3.0, 4.0});
}

TEST_CASE("Bounds<1>: SISO scalar conveniences") {
    constexpr Bounds<1> b{-2.0, 2.0};
    CHECK(b.saturate(5.0) == 2.0);
    CHECK(b.saturate(-5.0) == -2.0);
    CHECK(b.saturate(1.0) == 1.0);
    CHECK(b.contains(1.5));
    CHECK_FALSE(b.contains(2.5));
}
