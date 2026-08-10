// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Exercises the DAMP_BACKEND_ETL profile in isolation — the only place the ETL
// container backend is actually compiled and run. The damp:: vocabulary aliases
// resolve to etl:: types here, so this guards that the embeddable core builds
// against ETL and that the backend-agnostic helpers (the initializer_list
// min/max/minmax overloads in backend.hpp) behave identically to the std path.
// DAMP_BACKEND_ETL is selected via this directory's damp_profile.hpp.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "damp/backend.hpp"
#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

// Confirm we really are on the ETL backend, not silently the std default.
#ifndef DAMP_BACKEND_ETL
#error "ETL backend suite compiled without DAMP_BACKEND_ETL — check the include order of damp_profile.hpp"
#endif

#include <etl/type_traits.h>

using namespace damp;

TEST_SUITE("etl_backend") {
    TEST_CASE("damp vocabulary aliases resolve to ETL types") {
        static_assert(etl::is_same<damp::array<int, 3>, etl::array<int, 3>>::value);
        static_assert(etl::is_same<damp::pair<int, float>, etl::pair<int, float>>::value);

        damp::array<int, 3> a{1, 2, 3};
        CHECK(a[0] == 1);
        CHECK(a.size() == 3);

        damp::optional<int> o = 7;
        REQUIRE(o.has_value());
        CHECK(*o == 7);
        o = damp::nullopt;
        CHECK_FALSE(o.has_value());

        CHECK(damp::clamp(5, 0, 3) == 3);
    }

    TEST_CASE("binary min/max/minmax on the ETL backend") {
        CHECK(damp::min(2, 9) == 2);
        CHECK(damp::max(2, 9) == 9);
        CHECK(damp::minmax(9, 2) == damp::pair<int, int>{2, 9});
        CHECK(damp::minmax(2, 9) == damp::pair<int, int>{2, 9});
    }

    // The reason this suite exists: ETL's min/max/minmax are binary-only and
    // dangle, so backend.hpp supplies by-value initializer_list overloads. They
    // must match the std backend's behaviour (leftmost min, rightmost max).
    TEST_CASE("initializer_list min/max/minmax overloads work under ETL") {
        CHECK(damp::min({4, 2, 7, 1, 9}) == 1);
        CHECK(damp::max({4, 2, 7, 1, 9}) == 9);

        const auto [lo, hi] = damp::minmax({3.0, 1.0, 2.0, 5.0, -1.0});
        CHECK(lo == doctest::Approx(-1.0));
        CHECK(hi == doctest::Approx(5.0));

        // constexpr-usable on this backend too
        static_assert(damp::min({4, 2, 7}) == 2);
        static_assert(damp::max({4, 2, 7}) == 7);
        static_assert(damp::minmax({3, 1, 2}).first == 1);
        static_assert(damp::minmax({3, 1, 2}).second == 3);

        // single-element lists
        CHECK(damp::min({42}) == 42);
        CHECK(damp::minmax({42}) == damp::pair<int, int>{42, 42});
    }
}
