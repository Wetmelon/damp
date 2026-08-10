// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <complex>
#include <initializer_list>
#include <limits>
#include <numbers>
#include <utility>

#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

TEST_SUITE("constexpr_math") {
    TEST_CASE("damp::sqrt matches std::sqrt") {
        // Test positive values
        CHECK(damp::sqrt(0.0) == doctest::Approx(std::sqrt(0.0)));
        CHECK(damp::sqrt(1.0) == doctest::Approx(std::sqrt(1.0)));
        CHECK(damp::sqrt(2.0) == doctest::Approx(std::sqrt(2.0)));
        CHECK(damp::sqrt(4.0) == doctest::Approx(std::sqrt(4.0)));
        CHECK(damp::sqrt(9.0) == doctest::Approx(std::sqrt(9.0)));
        CHECK(damp::sqrt(100.0) == doctest::Approx(std::sqrt(100.0)));
        CHECK(damp::sqrt(0.25) == doctest::Approx(std::sqrt(0.25)));
        CHECK(damp::sqrt(1e-10) == doctest::Approx(std::sqrt(1e-10)));
        CHECK(damp::sqrt(1e10) == doctest::Approx(std::sqrt(1e10)));

        // Test float
        CHECK(damp::sqrt(2.0f) == doctest::Approx(std::sqrt(2.0)));
        CHECK(damp::sqrt(0.5f) == doctest::Approx(std::sqrt(0.5)));
    }

    TEST_CASE("damp::abs matches std::abs") {
        CHECK(damp::abs(0.0) == std::abs(0.0));
        CHECK(damp::abs(1.0) == std::abs(1.0));
        CHECK(damp::abs(-1.0) == std::abs(-1.0));
        CHECK(damp::abs(3.14159) == std::abs(3.14159));
        CHECK(damp::abs(-3.14159) == std::abs(-3.14159));
        CHECK(damp::abs(1e-15) == std::abs(1e-15));
        CHECK(damp::abs(-1e-15) == std::abs(-1e-15));
        CHECK(damp::abs(1e15) == std::abs(1e15));
        CHECK(damp::abs(-1e15) == std::abs(-1e15));
    }

    TEST_CASE("damp::hypot matches std::hypot and avoids overflow") {
        CHECK(damp::hypot(3.0, 4.0) == doctest::Approx(5.0));
        CHECK(damp::hypot(0.0, 0.0) == doctest::Approx(0.0));
        CHECK(damp::hypot(-3.0, 4.0) == doctest::Approx(5.0));
        CHECK(damp::hypot(5.0, 0.0) == doctest::Approx(5.0)); // exact on the axis
        CHECK(damp::hypot(1.0, 1.0) == doctest::Approx(std::sqrt(2.0)));
        // Overflow guard: 1e200^2 would overflow double, hypot must stay finite.
        CHECK(damp::hypot(1e200, 1e200) == doctest::Approx(std::hypot(1e200, 1e200)));
        // constexpr usable (consteval design path).
        static_assert(damp::hypot(3.0, 4.0) == 5.0);
    }

    TEST_CASE("damp::cbrt matches std::cbrt") {
        // Positive values
        CHECK(damp::cbrt(0.0) == doctest::Approx(std::cbrt(0.0)));
        CHECK(damp::cbrt(1.0) == doctest::Approx(std::cbrt(1.0)));
        CHECK(damp::cbrt(8.0) == doctest::Approx(std::cbrt(8.0)));
        CHECK(damp::cbrt(27.0) == doctest::Approx(std::cbrt(27.0)));
        CHECK(damp::cbrt(64.0) == doctest::Approx(std::cbrt(64.0)));
        CHECK(damp::cbrt(1000.0) == doctest::Approx(std::cbrt(1000.0)));

        // Negative values
        CHECK(damp::cbrt(-1.0) == doctest::Approx(std::cbrt(-1.0)));
        CHECK(damp::cbrt(-8.0) == doctest::Approx(std::cbrt(-8.0)));
        CHECK(damp::cbrt(-27.0) == doctest::Approx(std::cbrt(-27.0)));

        // Non-perfect cubes
        CHECK(damp::cbrt(2.0) == doctest::Approx(std::cbrt(2.0)));
        CHECK(damp::cbrt(10.0) == doctest::Approx(std::cbrt(10.0)));
        CHECK(damp::cbrt(-10.0) == doctest::Approx(std::cbrt(-10.0)));

        // Small and large values
        CHECK(damp::cbrt(1e-9) == doctest::Approx(std::cbrt(1e-9)));
        CHECK(damp::cbrt(1e9) == doctest::Approx(std::cbrt(1e9)));
    }

    TEST_CASE("damp::sqrt matches std::sqrt for complex") {
        using Cplx = damp::complex<double>;
        using StdComplex = std::complex<double>;

        // Real positive
        CHECK(damp::sqrt(Cplx(4.0, 0.0)).real() == doctest::Approx(std::sqrt(StdComplex(4.0, 0.0)).real()));
        CHECK(damp::sqrt(Cplx(4.0, 0.0)).imag() == doctest::Approx(std::sqrt(StdComplex(4.0, 0.0)).imag()));

        // Real negative (pure imaginary result)
        CHECK(damp::sqrt(Cplx(-4.0, 0.0)).real() == doctest::Approx(std::sqrt(StdComplex(-4.0, 0.0)).real()));
        CHECK(damp::sqrt(Cplx(-4.0, 0.0)).imag() == doctest::Approx(std::sqrt(StdComplex(-4.0, 0.0)).imag()));

        // Pure imaginary
        CHECK(damp::sqrt(Cplx(0.0, 4.0)).real() == doctest::Approx(std::sqrt(StdComplex(0.0, 4.0)).real()));
        CHECK(damp::sqrt(Cplx(0.0, 4.0)).imag() == doctest::Approx(std::sqrt(StdComplex(0.0, 4.0)).imag()));

        // General complex
        CHECK(damp::sqrt(Cplx(3.0, 4.0)).real() == doctest::Approx(std::sqrt(StdComplex(3.0, 4.0)).real()));
        CHECK(damp::sqrt(Cplx(3.0, 4.0)).imag() == doctest::Approx(std::sqrt(StdComplex(3.0, 4.0)).imag()));

        // Negative imaginary
        CHECK(damp::sqrt(Cplx(3.0, -4.0)).real() == doctest::Approx(std::sqrt(StdComplex(3.0, -4.0)).real()));
        CHECK(damp::sqrt(Cplx(3.0, -4.0)).imag() == doctest::Approx(std::sqrt(StdComplex(3.0, -4.0)).imag()));

        // Zero
        CHECK(damp::sqrt(Cplx(0.0, 0.0)).real() == doctest::Approx(0.0));
        CHECK(damp::sqrt(Cplx(0.0, 0.0)).imag() == doctest::Approx(0.0));
    }

    TEST_CASE("damp::atan2 matches std::atan2") {
        // Standard quadrants
        CHECK(damp::atan2(0.0, 1.0) == doctest::Approx(std::atan2(0.0, 1.0)));     // 0
        CHECK(damp::atan2(1.0, 1.0) == doctest::Approx(std::atan2(1.0, 1.0)));     // π/4
        CHECK(damp::atan2(1.0, 0.0) == doctest::Approx(std::atan2(1.0, 0.0)));     // π/2
        CHECK(damp::atan2(1.0, -1.0) == doctest::Approx(std::atan2(1.0, -1.0)));   // 3π/4
        CHECK(damp::atan2(0.0, -1.0) == doctest::Approx(std::atan2(0.0, -1.0)));   // π
        CHECK(damp::atan2(-1.0, -1.0) == doctest::Approx(std::atan2(-1.0, -1.0))); // -3π/4
        CHECK(damp::atan2(-1.0, 0.0) == doctest::Approx(std::atan2(-1.0, 0.0)));   // -π/2
        CHECK(damp::atan2(-1.0, 1.0) == doctest::Approx(std::atan2(-1.0, 1.0)));   // -π/4

        // Various ratios
        CHECK(damp::atan2(3.0, 4.0) == doctest::Approx(std::atan2(3.0, 4.0)));
        CHECK(damp::atan2(4.0, 3.0) == doctest::Approx(std::atan2(4.0, 3.0)));
        CHECK(damp::atan2(0.5, 0.1) == doctest::Approx(std::atan2(0.5, 0.1)));
        CHECK(damp::atan2(0.1, 0.5) == doctest::Approx(std::atan2(0.1, 0.5)));
    }

    TEST_CASE("damp::cos matches std::cos") {
        constexpr double pi = 3.14159265358979323846;

        // Standard angles
        CHECK(damp::cos(0.0) == doctest::Approx(std::cos(0.0)));
        CHECK(damp::cos(pi / 6) == doctest::Approx(std::cos(pi / 6)));                        // 30°
        CHECK(damp::cos(pi / 4) == doctest::Approx(std::cos(pi / 4)));                        // 45°
        CHECK(damp::cos(pi / 3) == doctest::Approx(std::cos(pi / 3)));                        // 60°
        CHECK(damp::cos(pi / 2) == doctest::Approx(std::cos(pi / 2)).epsilon(1e-10));         // 90°
        CHECK(damp::cos(pi) == doctest::Approx(std::cos(pi)));                                // 180°
        CHECK(damp::cos(3 * pi / 2) == doctest::Approx(std::cos(3 * pi / 2)).epsilon(1e-10)); // 270°

        // Negative angles
        CHECK(damp::cos(-pi / 4) == doctest::Approx(std::cos(-pi / 4)));
        CHECK(damp::cos(-pi / 2) == doctest::Approx(std::cos(-pi / 2)).epsilon(1e-10));
        CHECK(damp::cos(-pi) == doctest::Approx(std::cos(-pi)));

        // Arbitrary angles
        CHECK(damp::cos(0.5) == doctest::Approx(std::cos(0.5)));
        CHECK(damp::cos(1.0) == doctest::Approx(std::cos(1.0)));
        CHECK(damp::cos(2.0) == doctest::Approx(std::cos(2.0)));
        CHECK(damp::cos(3.0) == doctest::Approx(std::cos(3.0)));
    }

    TEST_CASE("damp::sin matches std::sin") {
        constexpr double pi = 3.14159265358979323846;

        // Standard angles
        CHECK(damp::sin(0.0) == doctest::Approx(std::sin(0.0)));
        CHECK(damp::sin(pi / 6) == doctest::Approx(std::sin(pi / 6)));         // 30°
        CHECK(damp::sin(pi / 4) == doctest::Approx(std::sin(pi / 4)));         // 45°
        CHECK(damp::sin(pi / 3) == doctest::Approx(std::sin(pi / 3)));         // 60°
        CHECK(damp::sin(pi / 2) == doctest::Approx(std::sin(pi / 2)));         // 90°
        CHECK(damp::sin(pi) == doctest::Approx(std::sin(pi)).epsilon(1e-10));  // 180°
        CHECK(damp::sin(3 * pi / 2) == doctest::Approx(std::sin(3 * pi / 2))); // 270°

        // Negative angles
        CHECK(damp::sin(-pi / 4) == doctest::Approx(std::sin(-pi / 4)));
        CHECK(damp::sin(-pi / 2) == doctest::Approx(std::sin(-pi / 2)));
        CHECK(damp::sin(-pi) == doctest::Approx(std::sin(-pi)).epsilon(1e-10));

        // Arbitrary angles
        CHECK(damp::sin(0.5) == doctest::Approx(std::sin(0.5)));
        CHECK(damp::sin(1.0) == doctest::Approx(std::sin(1.0)));
        CHECK(damp::sin(2.0) == doctest::Approx(std::sin(2.0)));
        CHECK(damp::sin(3.0) == doctest::Approx(std::sin(3.0)));
    }

    TEST_CASE("damp::tan matches std::tan") {
        constexpr double pi = std::numbers::pi;

        // Small angles (direct continued fraction path)
        CHECK(damp::tan(0.0) == doctest::Approx(std::tan(0.0)));
        CHECK(damp::tan(0.1) == doctest::Approx(std::tan(0.1)));
        CHECK(damp::tan(0.5) == doctest::Approx(std::tan(0.5)));
        CHECK(damp::tan(1.0) == doctest::Approx(std::tan(1.0)));

        // Standard angles
        CHECK(damp::tan(pi / 6) == doctest::Approx(std::tan(pi / 6))); // 30°
        CHECK(damp::tan(pi / 4) == doctest::Approx(std::tan(pi / 4))); // 45°
        CHECK(damp::tan(pi / 3) == doctest::Approx(std::tan(pi / 3))); // 60°

        // Near π/2 (complementary angle path)
        CHECK(damp::tan(1.5) == doctest::Approx(std::tan(1.5)).epsilon(1e-8));
        CHECK(damp::tan(1.55) == doctest::Approx(std::tan(1.55)).epsilon(1e-6));
        CHECK(damp::tan(1.57) == doctest::Approx(std::tan(1.57)).epsilon(1e-4));

        // Negative angles
        CHECK(damp::tan(-pi / 4) == doctest::Approx(std::tan(-pi / 4)));
        CHECK(damp::tan(-pi / 3) == doctest::Approx(std::tan(-pi / 3)));
        CHECK(damp::tan(-1.5) == doctest::Approx(std::tan(-1.5)).epsilon(1e-8));

        // Beyond first period (range reduction)
        CHECK(damp::tan(pi) == doctest::Approx(std::tan(pi)).epsilon(1e-10));
        CHECK(damp::tan(2.0) == doctest::Approx(std::tan(2.0)));
        CHECK(damp::tan(3.0) == doctest::Approx(std::tan(3.0)));
        CHECK(damp::tan(5.0) == doctest::Approx(std::tan(5.0)));
        CHECK(damp::tan(-5.0) == doctest::Approx(std::tan(-5.0)));

        // Constexpr verification
        constexpr double tan_val = damp::tan(pi / 4);
        CHECK(tan_val == doctest::Approx(1.0));
    }

    TEST_CASE("damp::asin matches std::asin") {
        for (int i = 0; i <= 20; ++i) {
            const double x = -1.0 + (0.1 * i);
            CHECK(damp::asin(x) == doctest::Approx(std::asin(x)).epsilon(1e-9));
        }
        CHECK(damp::asin(2.0) == doctest::Approx(std::numbers::pi / 2)); // clamped
    }

    TEST_CASE("damp::acos matches std::acos") {
        for (int i = 0; i <= 20; ++i) {
            const double x = -1.0 + (0.1 * i);
            CHECK(damp::acos(x) == doctest::Approx(std::acos(x)).epsilon(1e-9));
        }
        CHECK(damp::acos(-2.0) == doctest::Approx(std::numbers::pi)); // clamped
    }

    TEST_CASE("damp::atan matches std::atan") {
        for (int i = 0; i <= 20; ++i) {
            const double x = -5.0 + (0.5 * i);
            CHECK(damp::atan(x) == doctest::Approx(std::atan(x)).epsilon(1e-9));
        }
    }

    TEST_CASE("damp::fmod matches std::fmod") {
        CHECK(damp::fmod(7.0, 3.0) == doctest::Approx(std::fmod(7.0, 3.0)));
        CHECK(damp::fmod(-7.0, 3.0) == doctest::Approx(std::fmod(-7.0, 3.0))); // sign of dividend
        CHECK(damp::fmod(7.0, -3.0) == doctest::Approx(std::fmod(7.0, -3.0)));
        CHECK(damp::fmod(5.5, 2.0) == doctest::Approx(std::fmod(5.5, 2.0)));
        CHECK(damp::fmod(1.0, 0.0) == doctest::Approx(0.0)); // guarded
    }

    TEST_CASE("floor/ceil/nearbyint/fmod full range (beyond long long)") {
        // Regression: these used static_cast<long long> with no range check — UB for
        // |x| ≳ 2⁶³. Every finite binary float/double that large is already integral.
        constexpr double big = 1.0e20; // > 2⁶³ ≈ 9e18
        constexpr double huge = 1.0e300;
        constexpr float  bigf = 1.0e30f;

        CHECK(damp::floor(big) == std::floor(big));
        CHECK(damp::ceil(big) == std::ceil(big));
        CHECK(damp::nearbyint(big) == std::nearbyint(big));
        CHECK(damp::floor(-big) == std::floor(-big));
        CHECK(damp::ceil(-big) == std::ceil(-big));
        CHECK(damp::nearbyint(-big) == std::nearbyint(-big));
        CHECK(damp::floor(huge) == std::floor(huge));
        CHECK(damp::ceil(-huge) == std::ceil(-huge));
        CHECK(damp::floor(bigf) == std::floor(bigf));
        CHECK(damp::ceil(bigf) == std::ceil(bigf));

        // Constexpr path (forces detail::, not the runtime backend).
        static_assert(damp::floor(big) == big);
        static_assert(damp::ceil(big) == big);
        static_assert(damp::nearbyint(big) == big);
        static_assert(damp::floor(-big) == -big);
        static_assert(damp::ceil(-big) == -big);
        static_assert(damp::floor(3.7) == 3.0);
        static_assert(damp::floor(-3.7) == -4.0);
        static_assert(damp::ceil(3.2) == 4.0);
        static_assert(damp::ceil(-3.2) == -3.0);
        static_assert(damp::floor(bigf) == bigf);

        // Quotient x/y also used to overflow the long-long cast inside fmod.
        CHECK(damp::fmod(big, 3.0) == doctest::Approx(std::fmod(big, 3.0)));
        CHECK(damp::fmod(-big, 3.0) == doctest::Approx(std::fmod(-big, 3.0)));
        constexpr double fmod_big = damp::fmod(1.0e20, 3.0);
        CHECK(fmod_big == doctest::Approx(std::fmod(1.0e20, 3.0)));
    }

    TEST_CASE("damp::copysign matches std::copysign") {
        CHECK(damp::copysign(3.0, -2.0) == std::copysign(3.0, -2.0));
        CHECK(damp::copysign(3.0, 2.0) == std::copysign(3.0, 2.0));
        CHECK(damp::copysign(-3.0, 2.0) == std::copysign(-3.0, 2.0));
    }

    TEST_CASE("damp::isfinite: bit-pattern safe under -ffast-math") {
        // Test runner builds with -ffast-math. Bit-pattern isfinite must still
        // reject Inf/NaN constants at runtime and compile time (std::isfinite
        // is allowed to be a no-op under -ffinite-math-only).
        CHECK(damp::isfinite(1.0));
        CHECK(damp::isfinite(-1e300));
        CHECK_FALSE(damp::isfinite(std::numeric_limits<double>::infinity()));
        CHECK_FALSE(damp::isfinite(-std::numeric_limits<double>::infinity()));
        CHECK_FALSE(damp::isfinite(std::numeric_limits<double>::quiet_NaN()));
        CHECK_FALSE(damp::isfinite(std::numeric_limits<float>::infinity()));
        CHECK_FALSE(damp::isfinite(std::numeric_limits<float>::quiet_NaN()));

        static_assert(!damp::isfinite(std::numeric_limits<double>::infinity()));
        static_assert(!damp::isfinite(-std::numeric_limits<double>::infinity()));
        static_assert(!damp::isfinite(std::numeric_limits<double>::quiet_NaN()));
    }

    // Forces evaluation of the constexpr (series / Newton) path — NOT the runtime
    // MathBackend — by capturing into a constexpr variable, then compares against
    // std at runtime. The plain `damp::f(x)` runtime sweeps elsewhere in this file
    // exercise the backend; these exercise the compile-time implementation.
#define CEXPR_APPROX(expr, ref, eps)                    \
    do {                                                \
        constexpr double v_ = (expr);                   \
        CHECK(v_ == doctest::Approx(ref).epsilon(eps)); \
    } while (0)

    TEST_CASE("constexpr exp matches std (argument reduction over full range)") {
        CEXPR_APPROX(damp::exp(0.0), std::exp(0.0), 1e-15);
        CEXPR_APPROX(damp::exp(1.0), std::exp(1.0), 1e-14);
        CEXPR_APPROX(damp::exp(-1.0), std::exp(-1.0), 1e-14);
        CEXPR_APPROX(damp::exp(0.5), std::exp(0.5), 1e-14);
        CEXPR_APPROX(damp::exp(5.0), std::exp(5.0), 1e-13);
        CEXPR_APPROX(damp::exp(-5.0), std::exp(-5.0), 1e-13);
        CEXPR_APPROX(damp::exp(20.0), std::exp(20.0), 1e-12);
        CEXPR_APPROX(damp::exp(-20.0), std::exp(-20.0), 1e-12);
        CEXPR_APPROX(damp::exp(50.0), std::exp(50.0), 1e-12);
        CEXPR_APPROX(damp::exp(700.0), std::exp(700.0), 1e-11);
        CEXPR_APPROX(damp::exp(-700.0), std::exp(-700.0), 1e-11);

        // Overflow / underflow guards: saturate to max() / 0 (no ±inf under
        // -ffinite-math-only).
        static_assert(damp::exp(800.0) == std::numeric_limits<double>::max());
        static_assert(damp::exp(-800.0) == 0.0);

        // float specialization stays finite/zero past its narrower range
        CEXPR_APPROX(damp::exp(10.0f), std::exp(10.0), 1e-5);
        // Type-scaled saturation (float overflows ~88.7, not at the double 709 gate)
        static_assert(damp::exp(100.0f) == std::numeric_limits<float>::max());
        static_assert(damp::exp(-120.0f) == 0.0f);
    }

    TEST_CASE("constexpr log / pow / log10 match std") {
        CEXPR_APPROX(damp::log(1.0), std::log(1.0), 1e-15);
        CEXPR_APPROX(damp::log(2.0), std::log(2.0), 1e-13);
        CEXPR_APPROX(damp::log(0.5), std::log(0.5), 1e-13);
        CEXPR_APPROX(damp::log(1000.0), std::log(1000.0), 1e-13);
        CEXPR_APPROX(damp::log(1e-9), std::log(1e-9), 1e-13);
        CEXPR_APPROX(damp::log10(1000.0), std::log10(1000.0), 1e-13);
        CEXPR_APPROX(damp::log10(1e-9), std::log10(1e-9), 1e-13);
        CEXPR_APPROX(damp::pow(2.0, 10.0), std::pow(2.0, 10.0), 1e-12);
        CEXPR_APPROX(damp::pow(2.0, 0.5), std::pow(2.0, 0.5), 1e-13);
        CEXPR_APPROX(damp::pow(10.0, -3.0), std::pow(10.0, -3.0), 1e-12);
        static_assert(damp::log(0.0) == 0.0);  // guarded
        static_assert(damp::log(-1.0) == 0.0); // guarded
        static_assert(damp::log10(0.0) == 0.0);
        static_assert(damp::log10(-1.0) == 0.0);
        static_assert(damp::sqrt(-1.0) == 0.0); // domain guard before dispatch
        static_assert(damp::sqrt(-4.0f) == 0.0f);
        static_assert(damp::pow(-2.0, 3.0) == 0.0); // base ≤ 0
        static_assert(damp::pow(2.0, 0.0) == 1.0);  // exponent 0
        // Runtime path must agree (guards are outside is_constant_evaluated)
        CHECK(damp::sqrt(-1.0) == 0.0);
        CHECK(damp::log(-1.0) == 0.0);
        CHECK(damp::log10(0.0) == 0.0);
        CHECK(damp::pow(-1.0, 2.5) == 0.0);
    }

    TEST_CASE("constexpr sin/cos/tan match std for large arguments") {
        // The previous subtract-2π-in-a-loop reduction lost precision (and looped
        // ~10^6 times) for large arguments; the Cody–Waite reduction is O(1).
        CEXPR_APPROX(damp::sin(100.0), std::sin(100.0), 1e-12);
        CEXPR_APPROX(damp::cos(100.0), std::cos(100.0), 1e-12);
        CEXPR_APPROX(damp::sin(1000.0), std::sin(1000.0), 1e-11);
        CEXPR_APPROX(damp::cos(1000.0), std::cos(1000.0), 1e-11);
        CEXPR_APPROX(damp::sin(-500.0), std::sin(-500.0), 1e-11);
        CEXPR_APPROX(damp::cos(1e6), std::cos(1e6), 1e-7);
        CEXPR_APPROX(damp::sin(1e6), std::sin(1e6), 1e-7);
        CEXPR_APPROX(damp::tan(1000.0), std::tan(1000.0), 1e-10);

        // Regression: the complementary-angle path (|reduced r| > 1.2) once had a
        // sign error — tan(r) = 1/tan(π/2−r), not −1/tan(...). These compile-time
        // points sit squarely in that branch and would flip sign if it returns.
        CEXPR_APPROX(damp::tan(1.5), std::tan(1.5), 1e-9);
        CEXPR_APPROX(damp::tan(-1.6), std::tan(-1.6), 1e-9);
        static_assert(damp::tan(1.5) > 0.0, "tan(1.5) must be positive");
        static_assert(damp::tan(-1.6) > 0.0, "tan(-1.6) must be positive");
    }

    // Property-based checks evaluated entirely at compile time across a sweep.
    // These need no std reference values and so run inside consteval.
    TEST_CASE("constexpr identities hold at compile time") {
        constexpr auto pythag_ok = []() consteval {
            for (int i = -200; i <= 200; ++i) {
                const double x = i * 0.37;
                const double s = damp::sin(x);
                const double c = damp::cos(x);
                if (damp::abs((s * s) + (c * c) - 1.0) > 1e-12) {
                    return false;
                }
            }
            return true;
        };
        static_assert(pythag_ok(), "sin^2 + cos^2 == 1 across sweep");

        constexpr auto exp_log_roundtrip = []() consteval {
            for (int i = 1; i <= 100; ++i) {
                const double x = i * 0.5;
                if (damp::abs(damp::log(damp::exp(x)) - x) > 1e-9 * x) {
                    return false;
                }
            }
            return true;
        };
        static_assert(exp_log_roundtrip(), "log(exp(x)) == x across sweep");

        constexpr auto exp_addition = []() consteval {
            for (int i = -20; i <= 20; ++i) {
                const double a = i * 0.3;
                const double b = i * 0.17;
                if (damp::abs(damp::exp(a + b) - (damp::exp(a) * damp::exp(b))) > 1e-10 * damp::exp(a + b)) {
                    return false;
                }
            }
            return true;
        };
        static_assert(exp_addition(), "exp(a+b) == exp(a)*exp(b) across sweep");

        constexpr auto tan_quotient = []() consteval {
            for (int i = -30; i <= 30; ++i) {
                const double x = i * 0.1; // avoid exact π/2
                const double c = damp::cos(x);
                if (damp::abs(c) < 1e-3) {
                    continue;
                }
                const double t = damp::tan(x);
                if (damp::abs(t - (damp::sin(x) / c)) > 1e-9 * (1.0 + damp::abs(t))) {
                    return false;
                }
            }
            return true;
        };
        static_assert(tan_quotient(), "tan == sin/cos across sweep");

        constexpr auto sin_periodic = []() consteval {
            constexpr double two_pi = 6.283185307179586476925286766559;
            for (int i = -50; i <= 50; ++i) {
                const double x = i * 0.13;
                if (damp::abs(damp::sin(x) - damp::sin(x + (10.0 * two_pi))) > 1e-9) {
                    return false;
                }
            }
            return true;
        };
        static_assert(sin_periodic(), "sin is 2π-periodic across sweep");
    }

    TEST_CASE("constexpr verification") {
        // Verify all functions can be used in constexpr context
        constexpr double sqrt_val = damp::sqrt(4.0);
        constexpr double abs_val = damp::abs(-5.0);
        constexpr double cbrt_val = damp::cbrt(8.0);
        constexpr double atan2_val = damp::atan2(1.0, 1.0);
        constexpr double cos_val = damp::cos(0.0);
        constexpr double sin_val = damp::sin(0.0);
        constexpr auto   csqrt_val = damp::sqrt(damp::complex<double>(4.0, 0.0));

        constexpr double asin_val = damp::asin(1.0);
        constexpr double acos_val = damp::acos(0.0);
        constexpr double atan_val = damp::atan(1.0);
        constexpr double fmod_val = damp::fmod(7.0, 3.0);
        constexpr double csign_val = damp::copysign(3.0, -1.0);
        constexpr bool   finite_val = damp::isfinite(1.0);

        CHECK(sqrt_val == doctest::Approx(2.0));
        CHECK(abs_val == doctest::Approx(5.0));
        CHECK(cbrt_val == doctest::Approx(2.0));
        CHECK(atan2_val == doctest::Approx(std::atan2(1.0, 1.0)));
        CHECK(cos_val == doctest::Approx(1.0));
        CHECK(sin_val == doctest::Approx(0.0));
        CHECK(csqrt_val.real() == doctest::Approx(2.0));
        CHECK(csqrt_val.imag() == doctest::Approx(0.0));
        CHECK(asin_val == doctest::Approx(std::numbers::pi / 2));
        CHECK(acos_val == doctest::Approx(std::numbers::pi / 2));
        CHECK(atan_val == doctest::Approx(std::numbers::pi / 4));
        CHECK(fmod_val == doctest::Approx(1.0));
        CHECK(csign_val == doctest::Approx(-3.0));
        CHECK(finite_val);
    }

    TEST_CASE("damp::nearbyint ties to even (compile-time and runtime)") {
        // Non-tie values agree on every path and with std::nearbyint.
        for (double x : {-3.4, -1.6, -0.2, 0.2, 1.6, 3.4, 100.7, -100.7}) {
            CHECK(damp::nearbyint(x) == std::nearbyint(x));
        }

        // Half-integer ties → even integer (IEEE default / FE_TONEAREST).
        CHECK(damp::nearbyint(0.5) == 0.0);
        CHECK(damp::nearbyint(1.5) == 2.0);
        CHECK(damp::nearbyint(2.5) == 2.0);
        CHECK(damp::nearbyint(3.5) == 4.0);
        CHECK(damp::nearbyint(-0.5) == 0.0);
        CHECK(damp::nearbyint(-1.5) == -2.0);
        CHECK(damp::nearbyint(-2.5) == -2.0);

        // Compile-time path matches the same policy (detail::nearbyint).
        static_assert(damp::nearbyint(0.5) == 0.0);
        static_assert(damp::nearbyint(1.5) == 2.0);
        static_assert(damp::nearbyint(2.5) == 2.0);
        static_assert(damp::nearbyint(3.5) == 4.0);
        static_assert(damp::nearbyint(-0.5) == 0.0);
        static_assert(damp::nearbyint(-1.5) == -2.0);
        static_assert(damp::nearbyint(-2.5) == -2.0);
        static_assert(damp::nearbyint(2.4) == 2.0);
        static_assert(damp::nearbyint(2.6) == 3.0);
    }

    TEST_CASE("damp::wrap folds into the half-open interval") {
        constexpr double pi = std::numbers::pi;

        // Values already inside [min, max) are returned unchanged.
        CHECK(damp::wrap(0.0, -pi, pi) == doctest::Approx(0.0));
        CHECK(damp::wrap(1.0, -pi, pi) == doctest::Approx(1.0));
        CHECK(damp::wrap(-1.0, -pi, pi) == doctest::Approx(-1.0));

        // Out-of-range angles fold back by whole periods.
        CHECK(damp::wrap(pi + 0.5, -pi, pi) == doctest::Approx(0.5 - pi));
        CHECK(damp::wrap(-pi - 0.5, -pi, pi) == doctest::Approx(pi - 0.5));
        CHECK(damp::wrap(7.0, -pi, pi) == doctest::Approx(7.0 - (2.0 * pi)));

        // Non-symmetric interval (degrees).
        CHECK(damp::wrap(370.0, 0.0, 360.0) == doctest::Approx(10.0));
        CHECK(damp::wrap(-10.0, 0.0, 360.0) == doctest::Approx(350.0));

        // Equivalent to the atan2(sin, cos) phase reduction the PLL relies on.
        for (int i = -50; i <= 50; ++i) {
            const double a = i * 0.37;
            CHECK(damp::wrap(a, -pi, pi) == doctest::Approx(std::atan2(std::sin(a), std::cos(a))).epsilon(1e-12));
        }

        // Result always lands in [min, max).
        for (int i = -100; i <= 100; ++i) {
            const double w = damp::wrap(i * 1.1, -pi, pi);
            CHECK(w >= -pi);
            CHECK(w < pi);
        }
    }

    TEST_CASE("damp::wrap_pi is wrap(θ, −π, π)") {
        constexpr double pi = std::numbers::pi;
        CHECK(damp::wrap_pi(0.0) == doctest::Approx(0.0));
        CHECK(damp::wrap_pi(pi + 0.5) == doctest::Approx(0.5 - pi));
        CHECK(damp::wrap_pi(-pi - 0.5) == doctest::Approx(pi - 0.5));
        CHECK(damp::wrap_pi(7.0) == doctest::Approx(damp::wrap(7.0, -pi, pi)));
        for (int i = -40; i <= 40; ++i) {
            const double a = i * 0.41;
            CHECK(damp::wrap_pi(a) == doctest::Approx(damp::wrap(a, -pi, pi)));
            CHECK(damp::wrap_pi(a) == doctest::Approx(std::atan2(std::sin(a), std::cos(a))).epsilon(1e-12));
        }
        static_assert(damp::wrap_pi(0.0) == 0.0);
    }

    TEST_CASE("damp::wrap_two_pi is wrap(θ, 0, 2π)") {
        constexpr double pi = std::numbers::pi;
        constexpr double two_pi = 2.0 * pi;
        CHECK(damp::wrap_two_pi(0.0) == doctest::Approx(0.0));
        CHECK(damp::wrap_two_pi(-0.1) == doctest::Approx(two_pi - 0.1));
        CHECK(damp::wrap_two_pi(two_pi + 0.25) == doctest::Approx(0.25));
        CHECK(damp::wrap_two_pi(pi) == doctest::Approx(pi)); // 180° interior, not a seam
        for (int i = -40; i <= 40; ++i) {
            const double a = i * 0.41;
            const double w = damp::wrap_two_pi(a);
            CHECK(w == doctest::Approx(damp::wrap(a, 0.0, two_pi)));
            CHECK(w >= 0.0);
            CHECK(w < two_pi);
        }
        static_assert(damp::wrap_two_pi(0.0) == 0.0);
    }
}

TEST_SUITE("Utility") {
    TEST_CASE("unit conversion helpers support float") {
        const float deg = damp::rad2deg(3.14159265358979323846f);
        const float rad = damp::deg2rad(180.0f);

        CHECK(deg == doctest::Approx(180.0f).epsilon(1e-5));
        CHECK(rad == doctest::Approx(3.14159265358979323846f).epsilon(1e-5));
    }

    TEST_CASE("mag and db conversions support float") {
        const float db = damp::mag2db(10.0f);
        const float mag = damp::db2mag(20.0f);

        CHECK(db == doctest::Approx(20.0f).epsilon(1e-5));
        CHECK(mag == doctest::Approx(10.0f).epsilon(1e-5));
    }
}