// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <array>
#include <cmath>
#include <initializer_list>
#include <numbers>
#include <utility>

// Product path: public damp::sin/cos via MathBackend → detail::fast_* kernels.
#define DAMP_MATH_BACKEND_DAMP
#include "damp/math/math.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

// ~8 ULP accuracy on float32; epsilon(1e-5) gives comfortable margin
static constexpr float  kEps = 1e-5f;
static constexpr double kEpsD = 1e-12;

TEST_SUITE("trig") {

    TEST_CASE("sin - known values") {
        using std::numbers::pi_v;
        CHECK(damp::sin(0.0f) == doctest::Approx(std::sin(0.0f)).epsilon(kEps));
        CHECK(damp::sin(pi_v<float> / 6.0f) == doctest::Approx(0.5f).epsilon(kEps));
        CHECK(damp::sin(pi_v<float> / 4.0f) == doctest::Approx(std::sin(pi_v<float> / 4.0f)).epsilon(kEps));
        CHECK(damp::sin(pi_v<float> / 3.0f) == doctest::Approx(std::sin(pi_v<float> / 3.0f)).epsilon(kEps));
        CHECK(damp::sin(pi_v<float> / 2.0f) == doctest::Approx(1.0f).epsilon(kEps));
        CHECK(damp::sin(pi_v<float>) == doctest::Approx(std::sin(pi_v<float>)).epsilon(kEps));
        CHECK(damp::sin(3.0f * pi_v<float> / 2.0f) == doctest::Approx(-1.0f).epsilon(kEps));
        CHECK(damp::sin(2.0f * pi_v<float>) == doctest::Approx(std::sin(2.0f * pi_v<float>)).epsilon(kEps));
        CHECK(damp::sin(-pi_v<float> / 2.0f) == doctest::Approx(-1.0f).epsilon(kEps));
    }

    TEST_CASE("sin double - known values") {
        using std::numbers::pi_v;
        CHECK(damp::sin(0.0) == doctest::Approx(std::sin(0.0)).epsilon(kEpsD));
        CHECK(damp::sin(pi_v<double> / 6.0) == doctest::Approx(0.5).epsilon(kEpsD));
        CHECK(damp::sin(pi_v<double> / 2.0) == doctest::Approx(1.0).epsilon(kEpsD));
        CHECK(damp::sin(pi_v<double>) == doctest::Approx(std::sin(pi_v<double>)).epsilon(kEpsD));
        CHECK(damp::sin(-pi_v<double> / 2.0) == doctest::Approx(-1.0).epsilon(kEpsD));
        CHECK(damp::sin(100.0) == doctest::Approx(std::sin(100.0)).epsilon(kEpsD));
    }

    TEST_CASE("sincos double - matches sin and cos") {
        for (double x : {0.0, 0.5, 1.0, 1.5, 2.0, -0.7, 100.0}) {
            auto [s, c] = damp::sincos(x);
            CHECK(s == doctest::Approx(damp::sin(x)).epsilon(kEpsD));
            CHECK(c == doctest::Approx(damp::cos(x)).epsilon(kEpsD));
        }
    }

    TEST_CASE("asin/atan double - sample points") {
        CHECK(damp::asin(0.5) == doctest::Approx(std::asin(0.5)).epsilon(1e-10));
        CHECK(damp::acos(0.5) == doctest::Approx(std::acos(0.5)).epsilon(1e-10));
        CHECK(damp::atan(0.7) == doctest::Approx(std::atan(0.7)).epsilon(1e-10));
        CHECK(damp::atan2(0.7, -0.3) == doctest::Approx(std::atan2(0.7, -0.3)).epsilon(1e-10));
    }

    TEST_CASE("sin - odd symmetry") {
        for (float x : {0.1f, 0.5f, 1.0f, 1.5f, 2.0f, 3.0f}) {
            CHECK(damp::sin(-x) == doctest::Approx(-damp::sin(x)).epsilon(kEps));
        }
    }

    TEST_CASE("sin - matches std::sin") {
        for (float x : {0.1f, 0.3f, 0.7f, 1.0f, 1.3f, 1.9f, 2.4f, 3.0f, -0.5f, -1.2f, -2.7f}) {
            CHECK(damp::sin(x) == doctest::Approx(std::sin(x)).epsilon(kEps));
        }
    }

    TEST_CASE("sin - large arguments") {
        CHECK(damp::sin(10.0f * std::numbers::pi_v<float>) == doctest::Approx(std::sin(10.0f * std::numbers::pi_v<float>)).epsilon(kEps));
        CHECK(damp::sin(100.0f) == doctest::Approx(std::sin(100.0f)).epsilon(kEps));
        CHECK(damp::sin(-100.0f) == doctest::Approx(std::sin(-100.0f)).epsilon(kEps));
    }

    // -------------------------------------------------------------------------

    TEST_CASE("cos - known values") {
        using std::numbers::pi_v;
        CHECK(damp::cos(0.0f) == doctest::Approx(1.0f).epsilon(kEps));
        CHECK(damp::cos(pi_v<float> / 3.0f) == doctest::Approx(0.5f).epsilon(kEps));
        CHECK(damp::cos(pi_v<float> / 4.0f) == doctest::Approx(std::cos(pi_v<float> / 4.0f)).epsilon(kEps));
        CHECK(damp::cos(pi_v<float> / 2.0f) == doctest::Approx(std::cos(pi_v<float> / 2.0f)).epsilon(kEps));
        CHECK(damp::cos(pi_v<float>) == doctest::Approx(-1.0f).epsilon(kEps));
        CHECK(damp::cos(2.0f * pi_v<float>) == doctest::Approx(1.0f).epsilon(kEps));
    }

    TEST_CASE("cos - even symmetry") {
        for (float x : {0.1f, 0.5f, 1.0f, 1.5f, 2.0f, 3.0f}) {
            CHECK(damp::cos(-x) == doctest::Approx(damp::cos(x)).epsilon(kEps));
        }
    }

    TEST_CASE("cos - matches std::cos") {
        for (float x : {0.1f, 0.3f, 0.7f, 1.0f, 1.3f, 1.9f, 2.4f, 3.0f, -0.5f, -1.2f, -2.7f}) {
            CHECK(damp::cos(x) == doctest::Approx(std::cos(x)).epsilon(kEps));
        }
    }

    // -------------------------------------------------------------------------

    TEST_CASE("sincos - matches individual sin and cos") {
        for (float x : {0.0f, 0.5f, 1.0f, 1.5f, 2.0f, -0.7f, -2.3f}) {
            auto [s, c] = damp::sincos(x);
            CHECK(s == doctest::Approx(damp::sin(x)).epsilon(kEps));
            CHECK(c == doctest::Approx(damp::cos(x)).epsilon(kEps));
        }
    }

    TEST_CASE("sincos - Pythagorean identity") {
        for (float x : {0.0f, 0.5f, 1.0f, 1.5f, 2.0f, 2.5f, 3.0f, -1.0f, -2.5f}) {
            auto [s, c] = damp::sincos(x);
            CHECK((s * s) + (c * c) == doctest::Approx(1.0f).epsilon(kEps));
        }
    }

    // -------------------------------------------------------------------------

    TEST_CASE("asin - known values") {
        using std::numbers::pi_v;
        CHECK(damp::asin(0.0f) == doctest::Approx(0.0f).epsilon(kEps));
        CHECK(damp::asin(0.5f) == doctest::Approx(pi_v<float> / 6.0f).epsilon(kEps));
        CHECK(damp::asin(1.0f) == doctest::Approx(pi_v<float> / 2.0f).epsilon(kEps));
        CHECK(damp::asin(-0.5f) == doctest::Approx(-pi_v<float> / 6.0f).epsilon(kEps));
        CHECK(damp::asin(-1.0f) == doctest::Approx(-pi_v<float> / 2.0f).epsilon(kEps));
    }

    TEST_CASE("asin - clamped at boundaries") {
        using std::numbers::pi_v;
        CHECK(damp::asin(1.5f) == doctest::Approx(pi_v<float> / 2.0f).epsilon(kEps));
        CHECK(damp::asin(-1.5f) == doctest::Approx(-pi_v<float> / 2.0f).epsilon(kEps));
    }

    TEST_CASE("asin - matches std::asin") {
        for (float x : {0.0f, 0.1f, 0.3f, 0.5f, 0.7f, 0.9f, -0.2f, -0.6f, -0.95f}) {
            CHECK(damp::asin(x) == doctest::Approx(std::asin(x)).epsilon(kEps));
        }
    }

    // -------------------------------------------------------------------------

    TEST_CASE("acos - known values") {
        using std::numbers::pi_v;
        CHECK(damp::acos(1.0f) == doctest::Approx(0.0f).epsilon(kEps));
        CHECK(damp::acos(0.5f) == doctest::Approx(pi_v<float> / 3.0f).epsilon(kEps));
        CHECK(damp::acos(0.0f) == doctest::Approx(pi_v<float> / 2.0f).epsilon(kEps));
        CHECK(damp::acos(-0.5f) == doctest::Approx(2.0f * pi_v<float> / 3.0f).epsilon(kEps));
        CHECK(damp::acos(-1.0f) == doctest::Approx(pi_v<float>).epsilon(kEps));
    }

    TEST_CASE("acos - clamped at boundaries") {
        using std::numbers::pi_v;
        CHECK(damp::acos(1.5f) == doctest::Approx(0.0f).epsilon(kEps));
        CHECK(damp::acos(-1.5f) == doctest::Approx(pi_v<float>).epsilon(kEps));
    }

    TEST_CASE("acos - matches std::acos") {
        for (float x : {0.0f, 0.1f, 0.3f, 0.5f, 0.7f, 0.9f, -0.2f, -0.6f, -0.95f}) {
            CHECK(damp::acos(x) == doctest::Approx(std::acos(x)).epsilon(kEps));
        }
    }

    TEST_CASE("acos - asin + acos = pi/2") {
        for (float x : {0.0f, 0.2f, 0.5f, 0.8f, 1.0f}) {
            CHECK(damp::asin(x) + damp::acos(x) == doctest::Approx(std::numbers::pi_v<float> / 2.0f).epsilon(kEps));
        }
    }

    // -------------------------------------------------------------------------

    TEST_CASE("atan - known values") {
        using std::numbers::pi_v;
        CHECK(damp::atan(0.0f) == doctest::Approx(0.0f).epsilon(kEps));
        CHECK(damp::atan(1.0f) == doctest::Approx(pi_v<float> / 4.0f).epsilon(kEps));
        CHECK(damp::atan(-1.0f) == doctest::Approx(-pi_v<float> / 4.0f).epsilon(kEps));
    }

    TEST_CASE("atan - |x| > 1 complement path") {
        for (float x : {2.0f, 5.0f, 10.0f, 100.0f, -3.0f, -7.0f}) {
            CHECK(damp::atan(x) == doctest::Approx(std::atan(x)).epsilon(kEps));
        }
    }

    TEST_CASE("atan - matches std::atan") {
        for (float x : {0.0f, 0.1f, 0.3f, 0.5f, 0.7f, 0.9f, 1.0f, 1.5f, 2.0f, -0.4f, -1.0f, -2.5f}) {
            CHECK(damp::atan(x) == doctest::Approx(std::atan(x)).epsilon(kEps));
        }
    }

    TEST_CASE("atan - odd symmetry") {
        for (float x : {0.3f, 1.0f, 2.5f, 10.0f}) {
            CHECK(damp::atan(-x) == doctest::Approx(-damp::atan(x)).epsilon(kEps));
        }
    }

    // -------------------------------------------------------------------------

    TEST_CASE("atan2 - quadrant coverage") {
        using std::numbers::pi_v;
        // Q1
        CHECK(damp::atan2(1.0f, 1.0f) == doctest::Approx(pi_v<float> / 4.0f).epsilon(kEps));
        // Positive axes
        CHECK(damp::atan2(1.0f, 0.0f) == doctest::Approx(pi_v<float> / 2.0f).epsilon(kEps));
        CHECK(damp::atan2(0.0f, 1.0f) == doctest::Approx(0.0f).epsilon(kEps));
        // Q2
        CHECK(damp::atan2(1.0f, -1.0f) == doctest::Approx(3.0f * pi_v<float> / 4.0f).epsilon(kEps));
        // Negative x axis
        CHECK(damp::atan2(0.0f, -1.0f) == doctest::Approx(pi_v<float>).epsilon(kEps));
        // Q3 / negative y
        CHECK(damp::atan2(-1.0f, 0.0f) == doctest::Approx(-pi_v<float> / 2.0f).epsilon(kEps));
        CHECK(damp::atan2(-1.0f, -1.0f) == doctest::Approx(-3.0f * pi_v<float> / 4.0f).epsilon(kEps));
        // Q4
        CHECK(damp::atan2(-1.0f, 1.0f) == doctest::Approx(-pi_v<float> / 4.0f).epsilon(kEps));
    }

    TEST_CASE("atan2 - matches std::atan2") {
        const std::array pairs = std::to_array<damp::pair<float, float>>({
            {1.0f, 2.0f},
            {-1.0f, 2.0f},
            {1.0f, -2.0f},
            {-1.0f, -2.0f},
            {3.0f, 1.0f},
            {0.5f, 0.5f},
            {0.1f, 10.0f},
            {10.0f, 0.1f},
        });

        for (auto [y, x] : pairs) {
            CHECK(damp::atan2(y, x) == doctest::Approx(std::atan2(y, x)).epsilon(kEps));
        }
    }

    TEST_CASE("atan2 - consistent with atan for x > 0") {
        for (float t : {0.1f, 0.5f, 1.0f, 2.0f, 5.0f}) {
            CHECK(damp::atan2(t, 1.0f) == doctest::Approx(damp::atan(t)).epsilon(kEps));
        }
    }

    TEST_CASE("sqrt - matches std::sqrt") {
        for (float x : {0.0f, 1e-6f, 0.25f, 1.0f, 2.0f, 9.0f, 1234.5f, 1e6f}) {
            CHECK(damp::sqrt(x) == doctest::Approx(std::sqrt(x)).epsilon(kEps));
        }
    }

    TEST_CASE("nearbyint - rounds to nearest, ties to even") {
        // Exact integers and clear non-ties match std exactly (no rounding slack).
        for (float x : {-3.4f, -1.6f, -0.2f, 0.0f, 0.2f, 1.6f, 3.4f, 100.0f, -100.0f}) {
            CHECK(damp::nearbyint(x) == std::nearbyint(x));
        }
        // Half-integer ties round to even (FPSCR default / __builtin_nearbyintf).
        CHECK(damp::nearbyint(2.5f) == 2.0f);
        CHECK(damp::nearbyint(3.5f) == 4.0f);
        CHECK(damp::nearbyint(-2.5f) == -2.0f);
        CHECK(damp::nearbyint(0.5f) == 0.0f);
    }

} // TEST_SUITE("trig")
