// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Product path: public damp::sin/cos via MathBackend → detail::fast_* kernels. Built into
// test_damp_backend.exe with test_damp_backend.cpp: the DAMP MathBackend is ODR-incompatible
// with the std backend the main test_runner links against.
#define DAMP_MATH_BACKEND_DAMP

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <numbers>

#include "damp/backend.hpp"
#include "damp/math/math.hpp"
#include "damp/math/trig.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"
using namespace damp;

// Float kernels are good to about 1e-6 absolute; epsilon(1e-5) gives comfortable margin
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
        CHECK(damp::asin(0.5) == doctest::Approx(std::asin(0.5)).epsilon(kEpsD));
        CHECK(damp::acos(0.5) == doctest::Approx(std::acos(0.5)).epsilon(kEpsD));
        CHECK(damp::atan(0.7) == doctest::Approx(std::atan(0.7)).epsilon(kEpsD));
        CHECK(damp::atan2(0.7, -0.3) == doctest::Approx(std::atan2(0.7, -0.3)).epsilon(kEpsD));
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

    TEST_CASE("asin and acos are within 6e-7 absolute on [-1, 1]") {
        double worst_asin = 0.0;
        double worst_acos = 0.0;
        for (int i = 0; i <= 20000; ++i) {
            const float  x = static_cast<float>(-1.0 + (static_cast<double>(i) * 1e-4));
            const double xd = x;
            worst_asin = std::fmax(worst_asin, std::fabs(static_cast<double>(damp::asin(x)) - std::asin(xd)));
            worst_acos = std::fmax(worst_acos, std::fabs(static_cast<double>(damp::acos(x)) - std::acos(xd)));
        }
        CHECK(worst_asin <= 6e-7);
        CHECK(worst_acos <= 6e-7);
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


    TEST_CASE("fast atan2 is exact in ratio for tiny normal arguments") {
        constexpr float pi_4 = std::numbers::pi_v<float> / 4.0f;
        CHECK(std::fabs(detail::fast_atan2(1e-36f, 1e-36f) - pi_4) <= 1e-5f);
    }

    TEST_CASE("atan, asin and atan2 are odd: exactly 0 at 0, relative accuracy near 0") {
        CHECK(damp::atan(0.0f) == 0.0f);
        CHECK(damp::asin(0.0f) == 0.0f);
        CHECK(damp::atan2(0.0f, 1.0f) == 0.0f);
        for (const float x : {1e-7f, 1e-6f, 1e-5f, 1e-3f, 0.3f, 0.5f}) {
            const double rel = 4e-7;
            CHECK(std::fabs(static_cast<double>(damp::atan(x)) - std::atan(static_cast<double>(x))) <= rel * std::atan(static_cast<double>(x)));
            CHECK(std::fabs(static_cast<double>(damp::asin(x)) - std::asin(static_cast<double>(x))) <= rel * std::asin(static_cast<double>(x)));
            CHECK(std::fabs(static_cast<double>(damp::atan2(x, 1.0f)) - std::atan(static_cast<double>(x))) <= rel * std::atan(static_cast<double>(x)));
            CHECK(damp::atan(-x) == -damp::atan(x));
            CHECK(damp::asin(-x) == -damp::asin(x));
        }
    }

    TEST_CASE("double inverse trig is odd and exact at and near 0") {
        CHECK(damp::atan(0.0) == 0.0);
        CHECK(damp::asin(0.0) == 0.0);
        CHECK(damp::atan2(0.0, 1.0) == 0.0);
        CHECK(damp::acos(1.0) == 0.0);
        CHECK(damp::acos(-1.0) == std::numbers::pi_v<double>);
        CHECK(damp::asin(1.0) == std::numbers::pi_v<double> / 2.0);
        for (const double x : {1e-300, 1e-100, 1e-20, 1e-9}) { // x ∓ x³/3, x³/6 round to x
            CHECK(damp::atan(x) == x);
            CHECK(damp::asin(x) == x);
            CHECK(damp::atan2(x, 1.0) == x);
            CHECK(damp::atan(-x) == -x);
        }
        for (double x = 1e-6; x < 2e3; x *= 1.37) {
            CHECK(damp::atan(-x) == -damp::atan(x));
            CHECK(damp::atan2(-x, 1.0) == -damp::atan2(x, 1.0));
            CHECK(damp::atan2(-x, -1.0) == -damp::atan2(x, -1.0));
            if (x <= 1.0) {
                CHECK(damp::asin(-x) == -damp::asin(x));
            }
        }
    }

    TEST_CASE("double inverse trig is within a few ulp of a long double reference") {
        // Oracle in long double; asin and acos through (1 − x)(1 + x), which stays exact near ±1.
        auto ulps = [](double got, long double ref) {
            const double r = std::fabs(static_cast<double>(ref));
            const double ulp = std::nextafter(r, 2.0 * r + 1.0) - r;
            return static_cast<double>(std::fabs(static_cast<long double>(got) - ref)) / ulp;
        };
        auto ref_asin = [](double x) {
            const long double lx = x;
            return std::atan2(lx, std::sqrt((1.0L - lx) * (1.0L + lx)));
        };
        auto ref_acos = [](double x) {
            const long double lx = x;
            return std::atan2(std::sqrt((1.0L - lx) * (1.0L + lx)), lx);
        };
        double worst_asin = 0, worst_acos = 0, worst_atan = 0, worst_atan2 = 0;
        for (int i = -20000; i <= 20000; ++i) {
            const double x = i / 20000.0;
            worst_asin = std::fmax(worst_asin, ulps(damp::asin(x), ref_asin(x)));
            worst_acos = std::fmax(worst_acos, ulps(damp::acos(x), ref_acos(x)));
            worst_atan = std::fmax(worst_atan, ulps(damp::atan(x), std::atan(static_cast<long double>(x))));
        }
        for (double d = 1e-16; d < 1e-2; d *= 1.3) { // the sqrt branch near ±1
            for (const double x : {1.0 - d, -1.0 + d}) {
                worst_asin = std::fmax(worst_asin, ulps(damp::asin(x), ref_asin(x)));
                worst_acos = std::fmax(worst_acos, ulps(damp::acos(x), ref_acos(x)));
            }
        }
        for (double x = 1e-300; x < 1e300; x *= 1.07) { // all three atan reductions, both signs
            worst_atan = std::fmax(worst_atan, ulps(damp::atan(x), std::atan(static_cast<long double>(x))));
            worst_atan = std::fmax(worst_atan, ulps(damp::atan(-x), std::atan(static_cast<long double>(-x))));
        }
        for (int i = 0; i < 8192; ++i) { // every octant, |y| ≠ |x|
            const double th = std::numbers::pi * ((i + 0.5) / 4096.0 - 1.0);
            for (const double r : {1e-200, 1.0, 1e200}) {
                const double      y = r * std::sin(th), x = r * std::cos(th);
                const long double ref = std::atan2(static_cast<long double>(y), static_cast<long double>(x));
                worst_atan2 = std::fmax(worst_atan2, ulps(damp::atan2(y, x), ref));
            }
        }
        // Measured against libquadmath: asin 3.3, acos 3.1, atan 2.1, atan2 2.1 ulp.
        CHECK(worst_asin <= 4.0);
        CHECK(worst_acos <= 4.0);
        CHECK(worst_atan <= 3.0);
        CHECK(worst_atan2 <= 3.0);
    }

    TEST_CASE("double atan2 quadrants and axes") {
        using std::numbers::pi_v;
        CHECK(damp::atan2(0.0, -1.0) == pi_v<double>);
        CHECK(damp::atan2(-0.0, -1.0) == -pi_v<double>);
        CHECK(damp::atan2(1.0, 0.0) == pi_v<double> / 2.0);
        CHECK(damp::atan2(-1.0, 0.0) == -pi_v<double> / 2.0);
        CHECK(damp::atan2(0.0, 0.0) == 0.0);
        CHECK(damp::atan2(1.0, 1.0) == doctest::Approx(pi_v<double> / 4.0).epsilon(1e-15));
        CHECK(damp::atan2(1.0, -1.0) == doctest::Approx(3.0 * pi_v<double> / 4.0).epsilon(1e-15));
        CHECK(damp::atan2(-1.0, -1.0) == doctest::Approx(-3.0 * pi_v<double> / 4.0).epsilon(1e-15));
        CHECK(damp::atan(1e300) == pi_v<double> / 2.0);
        CHECK(damp::asin(1.5) == pi_v<double> / 2.0); // clamped, as the float kernel
        CHECK(damp::acos(-1.5) == pi_v<double>);
    }

    TEST_CASE("sin and cos propagate NaN and Inf as NaN through the reduction") {
        // Run-time NaN and Inf: clang -ffast-math folds arithmetic on a NaN or Inf constant.
        volatile float  zf = 0.0f;
        volatile double zd = 0.0;
        for (const float a : {zf / zf, 1.0f / zf, -1.0f / zf}) {
            CHECK_FALSE(damp::isfinite(damp::sin(a)));
            CHECK_FALSE(damp::isfinite(damp::cos(a)));
            const auto [s, c] = damp::sincos(a);
            CHECK_FALSE(damp::isfinite(s));
            CHECK_FALSE(damp::isfinite(c));
        }
        for (const double a : {zd / zd, 1.0 / zd, -1.0 / zd}) {
            CHECK_FALSE(damp::isfinite(damp::sin(a)));
            CHECK_FALSE(damp::isfinite(damp::cos(a)));
            const auto [s, c] = damp::sincos(a);
            CHECK_FALSE(damp::isfinite(s));
            CHECK_FALSE(damp::isfinite(c));
        }
    }

    TEST_CASE("float reduction: angles whose quotient rounds to 2^21 half-turns take the parity path") {
        constexpr float inv_pi = std::numbers::inv_pi_v<float>;
        const float     last_hot = std::bit_cast<float>(0x4ac90fd7U);
        const float     first_cold = std::bit_cast<float>(0x4ac90fd8U);
        CHECK(std::nearbyint(last_hot * inv_pi) < 0x1p21f);
        CHECK(std::nearbyint(first_cold * inv_pi) == 0x1p21f);
        for (const float sign : {1.0f, -1.0f}) {
            CHECK(std::fabs(damp::sin(sign * last_hot)) <= 1.0f);
            CHECK(std::fabs(damp::cos(sign * last_hot)) <= 1.0f);
            // Past the split the angle is a multiple of π with n's parity: sin ±0, cos ±1.
            for (const float x : {first_cold, 6588400.5f, 1e9f, 3e38f}) {
                const float n = std::nearbyint(x * inv_pi);
                const float even = std::fmod(n, 2.0f) == 0.0f ? 1.0f : -1.0f;
                const auto [s, c] = damp::sincos(sign * x);
                CHECK(s == 0.0f);
                CHECK(c == even);
                CHECK(damp::cos(sign * x) == even);
                CHECK(damp::sin(sign * x) == 0.0f);
            }
        }
    }

    TEST_CASE("double reduction: both sides of the 2^30·π hot/cold split agree with std::sin") {
        const double split = 0x1p30 * std::numbers::pi;
        for (const double x : {split * (1.0 - 1e-12), split * (1.0 + 1e-12), -split * (1.0 - 1e-12), -split * (1.0 + 1e-12), 1e10}) {
            const double tol = 4.0 * (std::nextafter(std::fabs(x), 2.0 * std::fabs(x)) - std::fabs(x)); // two-word π past 2²¹ half-turns
            CHECK(std::fabs(damp::sin(x) - std::sin(x)) <= tol);
            CHECK(std::fabs(damp::cos(x) - std::cos(x)) <= tol);
        }
        // Past 2⁵⁰ half-turns: a multiple of π with n's parity (cos_poly(0) is 1 − 7e-16).
        for (const double x : {0x1p60, -0x1p60, 1e300}) {
            const auto [s, c] = damp::sincos(x);
            CHECK(s == 0.0);
            CHECK(std::fabs(std::fabs(c) - 1.0) <= 1e-15);
        }
        // Double π split: the low word is π − HI (it was fl(π) − HI, worth n·1.2e-16 per half-turn).
        CHECK(std::fabs(damp::sin(1000.0) - std::sin(1000.0)) <= 1e-15);
        CHECK(std::fabs(damp::sin(1e6) - std::sin(1e6)) <= 1e-14);
    }

    TEST_CASE("float sin stays accurate at 2e4 rad under -ffast-math (reassociation barrier)") {
        for (const float x : {20000.3f, -19999.7f, 24344.2f}) {
            CHECK(std::fabs(static_cast<double>(damp::sin(x)) - std::sin(static_cast<double>(x))) <= 2e-6);
        }
    }

    TEST_CASE("float sin and cos keep their accuracy out to 2.5e4 rad (exact π split)") {
        double worst_sin = 0.0;
        double worst_cos = 0.0;
        for (int i = 0; i <= 200000; ++i) {
            const float  x = static_cast<float>(-2.5e4 + (static_cast<double>(i) * 0.25));
            const double xd = x;
            worst_sin = std::fmax(worst_sin, std::fabs(static_cast<double>(damp::sin(x)) - std::sin(xd)));
            worst_cos = std::fmax(worst_cos, std::fabs(static_cast<double>(damp::cos(x)) - std::cos(xd)));
        }
        CHECK(worst_sin <= 9e-7);
        CHECK(worst_cos <= 3.5e-7);
    }
} // TEST_SUITE("trig")
