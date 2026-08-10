// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "damp/matrix/rowvec.hpp"
#include "damp/toolbox/scaling.hpp"
#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

// Scaling / interpolation / calibration helpers.

TEST_SUITE("Scaling") {
    TEST_CASE("lerp and inverse_lerp round-trip") {
        CHECK(lerp(10.0, 20.0, 0.0) == doctest::Approx(10.0));
        CHECK(lerp(10.0, 20.0, 1.0) == doctest::Approx(20.0));
        CHECK(lerp(10.0, 20.0, 0.25) == doctest::Approx(12.5));
        CHECK(lerp(0.0, 10.0, -0.5) == doctest::Approx(-5.0)); // extrapolates
        // inverse recovers the fraction.
        CHECK(inverse_lerp(10.0, 20.0, 12.5) == doctest::Approx(0.25));
        const double x = 7.3;
        CHECK(lerp(2.0, 9.0, inverse_lerp(2.0, 9.0, x)) == doctest::Approx(x));
    }

    TEST_CASE("inverse_lerp degenerate range returns 0 (no Inf/NaN)") {
        CHECK(inverse_lerp(5.0, 5.0, 5.0) == doctest::Approx(0.0));
        CHECK(inverse_lerp(5.0, 5.0, 99.0) == doctest::Approx(0.0));
        // rescale through a collapsed input span stays at the lower output endpoint.
        CHECK(rescale(3.0, 1.0, 1.0, 10.0, 20.0) == doctest::Approx(10.0));
    }

    TEST_CASE("rescale maps between ranges (ADC counts to volts)") {
        CHECK(rescale(0.0f, 0.0f, 4095.0f, 0.0f, 3.3f) == doctest::Approx(0.0f));
        CHECK(rescale(4095.0f, 0.0f, 4095.0f, 0.0f, 3.3f) == doctest::Approx(3.3f));
        CHECK(rescale(2047.5f, 0.0f, 4095.0f, 0.0f, 3.3f) == doctest::Approx(1.65f));
        // Inverted output range works too.
        CHECK(rescale(0.0f, 0.0f, 10.0f, 100.0f, 0.0f) == doctest::Approx(100.0f));
    }

    TEST_CASE("AffineCal apply/invert round-trip") {
        const AffineCal<double> cal{2.0, -3.0}; // y = 2x - 3
        CHECK(cal.apply(5.0) == doctest::Approx(7.0));
        CHECK(cal.invert(7.0) == doctest::Approx(5.0));
        CHECK(cal.invert(cal.apply(1.234)) == doctest::Approx(1.234));
        // .as<U>() re-casts.
        const auto f = cal.as<float>();
        CHECK(f.apply(5.0f) == doctest::Approx(7.0f));
    }

    TEST_CASE("two_point_cal reproduces both anchor points") {
        // 0.5 V -> 0 °C, 4.5 V -> 100 °C.
        const auto cal = two_point_cal(0.5, 0.0, 4.5, 100.0);
        CHECK(cal.apply(0.5) == doctest::Approx(0.0));
        CHECK(cal.apply(4.5) == doctest::Approx(100.0));
        CHECK(cal.apply(2.5) == doctest::Approx(50.0)); // midpoint
    }

    TEST_CASE("poly_horner evaluates ascending-order coefficients") {
        // 1 + 2x + 3x^2 at x = 2 -> 1 + 4 + 12 = 17.
        constexpr RowVec<3, double> c{1.0, 2.0, 3.0};
        CHECK(poly_horner(c, 2.0) == doctest::Approx(17.0));
        CHECK(poly_horner(c, 0.0) == doctest::Approx(1.0)); // constant term
        // Single coefficient is a constant.
        constexpr RowVec<1, double> k{42.0};
        CHECK(poly_horner(k, 999.0) == doctest::Approx(42.0));
    }
}
