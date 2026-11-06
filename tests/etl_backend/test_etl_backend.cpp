// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Exercises the DAMP_BACKEND_ETL + DAMP_MATH_BACKEND_FREESTANDING profile in isolation:
// the only place the ETL container backend and the series math backend are compiled
// and run. The damp:: vocabulary aliases resolve to etl:: types here, so this guards
// that the embeddable core builds against ETL, that the backend-agnostic helpers (the
// initializer_list min/max/minmax overloads in backend.hpp) behave identically to the
// std path, and that the series math gives at run time what it gives in constant
// evaluation. Both macros are selected via this directory's damp_profile.hpp.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <cstddef>

#include "damp/backend.hpp"
#include "damp/control.hpp"
#include "damp/math/math.hpp"
#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

// Confirm we really are on the ETL backend, not silently the std default.
#ifndef DAMP_BACKEND_ETL
#error "ETL backend suite compiled without DAMP_BACKEND_ETL — check the include order of damp_profile.hpp"
#endif
#ifndef DAMP_MATH_BACKEND_FREESTANDING
#error "ETL backend suite compiled without DAMP_MATH_BACKEND_FREESTANDING — check damp_profile.hpp"
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

namespace {

/** A value the optimizer cannot fold, so the call runs the run-time backend. */
template<typename T>
T opaque(T v) {
    volatile T x = v;
    return x;
}

constexpr damp::array<double, 7> kXs{0.3, -1.7, 2.0, 8.0, 12.5, -40.25, 1.0e4};

/** f applied to every input in constant evaluation. */
template<typename F>
constexpr damp::array<double, 7> constant_evaluated(F f) {
    damp::array<double, 7> out{};
    for (size_t i = 0; i < kXs.size(); ++i) {
        out[i] = f(kXs[i]);
    }
    return out;
}

} // namespace

TEST_SUITE("freestanding math") {
    TEST_CASE("the series backend gives at run time what it gives in constant evaluation") {
        constexpr auto sin_ct = constant_evaluated([](double x) { return damp::sin(x); });
        constexpr auto cos_ct = constant_evaluated([](double x) { return damp::cos(x); });
        constexpr auto exp_ct = constant_evaluated([](double x) { return damp::exp(x * 0.01); });
        constexpr auto log_ct = constant_evaluated([](double x) { return damp::log(damp::abs(x)); });
        constexpr auto sqrt_ct = constant_evaluated([](double x) { return damp::sqrt(damp::abs(x)); });
        constexpr auto atan2_ct = constant_evaluated([](double x) { return damp::atan2(x, 1.5); });
        constexpr auto fmod_ct = constant_evaluated([](double x) { return damp::fmod(x, 0.5); });
        constexpr auto pow_ct = constant_evaluated([](double x) { return damp::pow(x, 4.0); });
        for (size_t i = 0; i < kXs.size(); ++i) {
            const double x = opaque(kXs[i]);
            CHECK(damp::sin(x) == doctest::Approx(sin_ct[i]).epsilon(1e-12));
            CHECK(damp::cos(x) == doctest::Approx(cos_ct[i]).epsilon(1e-12));
            CHECK(damp::exp(x * 0.01) == doctest::Approx(exp_ct[i]).epsilon(1e-12));
            CHECK(damp::log(damp::abs(x)) == doctest::Approx(log_ct[i]).epsilon(1e-12));
            CHECK(damp::sqrt(damp::abs(x)) == doctest::Approx(sqrt_ct[i]).epsilon(1e-12));
            CHECK(damp::atan2(x, 1.5) == doctest::Approx(atan2_ct[i]).epsilon(1e-12));
            CHECK(damp::fmod(x, 0.5) == doctest::Approx(fmod_ct[i]).epsilon(1e-12));
            CHECK(damp::pow(x, 4.0) == doctest::Approx(pow_ct[i]).epsilon(1e-12));
        }
        // Exact cases: x/y a power of two, and an even power of a negative base.
        CHECK(damp::fmod(opaque(2.0), 1.0) == 0.0);
        CHECK(damp::fmod(opaque(8.0F), 2.0F) == 0.0F);
        CHECK(damp::pow(opaque(-2.0), 4.0) == doctest::Approx(16.0));
    }

    TEST_CASE("a constexpr design and its tick build and run on ETL + freestanding math") {
        constexpr auto sys = damp::StateSpace{
            .A = damp::Matrix<2, 2>{{0.0, 1.0}, {9.81, -0.1}},
            .B = damp::Matrix<2, 1>{{0.0}, {1.0}},
            .C = damp::Matrix<1, 2>{{1.0, 0.0}},
        };
        constexpr auto res = damp::design::lqr(sys, 0.01, {.Q = damp::Matrix<2, 2>::identity() * 10.0, .R = 1.0});
        static_assert(res.success);
        damp::LQR law{res.as<float>()};
        law.update(damp::ColVec<2, float>{opaque(0.1F), 0.0F});
        const float u = law.control()(0, 0);
        CHECK(damp::isfinite(u));
        CHECK(u < 0.0F); // the angle is pushed back toward upright
    }
}
