// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <bit>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <numbers>

#include "damp/math/complex.hpp"
#include "damp/math/constexpr_math.hpp"
#include "damp/math/math.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {
constexpr double kEps = std::numeric_limits<double>::epsilon();

using C = damp::complex<double>;

} // namespace

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

    TEST_CASE("damp::fmod is 0 when x/y is a power of two, and pow keeps the sign of an even exponent") {
        // Regression: the binary reduction stopped doubling at d = x/2 and returned y.
        static_assert(damp::fmod(2.0, 1.0) == 0.0);
        static_assert(damp::fmod(8.0, 2.0) == 0.0);
        static_assert(damp::fmod(-2.0, 0.25) == 0.0);
        static_assert(damp::fmod(4.0f, 1.0f) == 0.0f);
        static_assert(damp::fmod(2.0, -1.0) == 0.0);
        static_assert(damp::fmod(9.0, 2.0) == 1.0);
        static_assert(damp::abs(damp::pow(-2.0, 4.0) - 16.0) < 1e-12);
        static_assert(damp::pow(-3.0, 8.0) > 0.0);
        static_assert(damp::pow(-2.0, 3.0) < 0.0);
        volatile double two = 2.0;
        volatile double one = 1.0;
        CHECK(damp::fmod(two, one) == 0.0);
        CHECK(damp::pow(-two, 4.0) == doctest::Approx(16.0));
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
        static_assert(damp::pow(-2.0, 2.5) == 0.0); // base < 0, non-integral exponent
        static_assert(damp::pow(0.0, 3.0) == 0.0);  // base 0
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

    TEST_CASE("damp::wrap keeps [min, max) for large float arguments and an empty interval") {
        constexpr float pi = std::numbers::pi_v<float>;
        for (float x = 1e3F; x < 3e38F; x *= 1.9F) {
            for (const float a : {x, -x}) {
                const float w = damp::wrap_pi(a);
                CHECK(w >= -pi);
                CHECK(w < pi);
                const float w2 = damp::wrap_two_pi(a);
                CHECK(w2 >= 0.0F);
                CHECK(w2 < 2.0F * pi);
            }
        }
        CHECK(damp::wrap(1.0, 2.0, 2.0) == 2.0);
        static_assert(damp::wrap_pi(1e9F) >= -pi && damp::wrap_pi(1e9F) < pi);
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

    TEST_CASE("series sqrt and cbrt converge across the exponent range") {
        for (const double x : {1e-300, 1e-100, 1e-10, 0.3, 2.0, 1e10, 1e100, 1e300}) {
            CHECK(std::fabs(detail::sqrt(x) - std::sqrt(x)) <= 2 * kEps * std::sqrt(x));
            CHECK(std::fabs(detail::cbrt(x) - std::cbrt(x)) <= 2 * kEps * std::cbrt(x));
            CHECK(std::fabs(detail::cbrt(-x) + std::cbrt(x)) <= 2 * kEps * std::cbrt(x));
        }
        CHECK(std::fabs(detail::sqrt(1e30f) - 1e15f) <= 2 * std::numeric_limits<float>::epsilon() * 1e15f);
    }

    TEST_CASE("constant-evaluated roots are accurate far from 1, subnormals included") {
        constexpr double denorm = 4.9406564584124654e-324; // smallest subnormal
        static_assert(detail::abs(damp::sqrt(1e150) - 1e75) <= 2 * kEps * 1e75);
        static_assert(detail::abs(damp::sqrt(1e-150) - 1e-75) <= 2 * kEps * 1e-75);
        static_assert(detail::abs(damp::sqrt(1e300) - 1e150) <= 2 * kEps * 1e150);
        static_assert(detail::abs(damp::sqrt(denorm) - 2.2227587494850775e-162) <= 2 * kEps * 2.3e-162);
        static_assert(detail::abs(damp::cbrt(1e-300) - 1e-100) <= 2 * kEps * 1e-100);
        static_assert(detail::abs(damp::cbrt(-8e150) + 2e50) <= 2 * kEps * 2e50);
        static_assert(detail::abs(damp::sqrt(1e30f) - 1e15f) <= 2 * std::numeric_limits<float>::epsilon() * 1e15f);
        static_assert(detail::abs(damp::log(1e300) - 690.77552789821368) <= 4 * kEps * 690.8);
        static_assert(detail::abs(damp::log(denorm) + 744.44007192138122) <= 4 * kEps * 744.5);
        static_assert(detail::abs(damp::hypot(3e300, 4e300) - 5e300) <= 4 * kEps * 5e300);
    }

    TEST_CASE("series log/exp/pow/atan2 hold accuracy far from 1") {
        for (const double x : {1e-300, 1e-150, 1e-10, 0.7, 1.0, 1.3, 1e10, 1e150, 1e300}) {
            CHECK(std::fabs(detail::log(x) - std::log(x)) <= (4 * kEps * std::fabs(std::log(x))) + kEps);
        }
        for (const double x : {-700.0, -300.0, -1.0, 0.5, 1.0, 300.0, 700.0}) {
            CHECK(std::fabs(detail::exp(x) - std::exp(x)) <= 8 * kEps * std::exp(x));
        }
        // exp(y·ln b): the rounding of y·ln b is amplified by |y·ln b| (≈ 690 here).
        CHECK(std::fabs(detail::pow(10.0, 300.0) - 1e300) <= 1024 * kEps * 1e300);
        CHECK(std::fabs(detail::pow(10.0, -300.0) - 1e-300) <= 1024 * kEps * 1e-300);
        CHECK(std::fabs(detail::atan2(1e300, 1e-300) - (std::numbers::pi / 2.0)) <= 2 * kEps);
        CHECK(std::fabs(detail::atan2(1e-300, 1e300)) <= 1e-300);
        CHECK(std::fabs(detail::atan2(-1e-300, -1e-300) + (0.75 * std::numbers::pi)) <= 4 * kEps);
    }

    TEST_CASE("series sin/cos/sincos are accurate over several periods") {
        double worst = 0.0;
        for (int i = -2000; i <= 2000; ++i) {
            const double x = 0.01 * i;
            const auto   sc = detail::sincos(x);
            worst = std::fmax(worst, std::fabs(detail::sin(x) - std::sin(x)));
            worst = std::fmax(worst, std::fabs(detail::cos(x) - std::cos(x)));
            worst = std::fmax(worst, std::fabs(sc.sin - std::sin(x)));
            worst = std::fmax(worst, std::fabs(sc.cos - std::cos(x)));
        }
        CHECK(worst <= 4 * kEps);
        CHECK(std::fabs(detail::cos(std::numbers::pi) + 1.0) <= 2 * kEps);
        // Same kernels under constant evaluation (the consteval design path).
        static_assert(detail::abs(damp::cos(std::numbers::pi) + 1.0) <= 2 * kEps);
        static_assert(detail::abs(damp::sqrt(1e300) - 1e150) <= 2 * kEps * 1e150);
    }

    TEST_CASE("series tan keeps relative accuracy next to pi/2") {
        for (const double x : {1.5707963267, 1.57079632679, -1.5707963267, 1.0, 0.1, 4.0}) {
            CHECK(std::fabs(detail::tan(x) - std::tan(x)) <= 16 * kEps * std::fabs(std::tan(x)));
        }
    }

    TEST_CASE("real pow keeps the sign of a negative base with an integral exponent") {
        CEXPR_APPROX(damp::pow(-2.0, 3.0), -8.0, 1e-13);
        CEXPR_APPROX(damp::pow(-2.0, -3.0), -0.125, 1e-13);
        CEXPR_APPROX(damp::pow(-0.5, 2.0), 0.25, 1e-13);
        CHECK(std::fabs(detail::pow(-2.0, 3.0) + 8.0) <= 16 * kEps); // series backend body
        volatile double base = -2.0;                                 // runtime backend path
        CHECK(damp::pow(static_cast<double>(base), 3.0) == doctest::Approx(-8.0));
        CHECK(damp::pow(static_cast<double>(base), 2.0) == doctest::Approx(4.0));
        CHECK(damp::pow(-1e200, 3.0) == -std::numeric_limits<double>::max());
    }

    TEST_CASE("pow with a mixed float/double exponent keeps the real exponent") {
        volatile float xf = 2.0F;
        CHECK(damp::pow(static_cast<float>(xf), 0.5) == doctest::Approx(std::sqrt(2.0F)));
        CHECK(damp::pow(2.0, 2.5F) == doctest::Approx(std::pow(2.0, 2.5)));
        constexpr float cf = damp::pow(2.0F, 0.5);
        CHECK(cf == doctest::Approx(std::sqrt(2.0F)));
        static_assert(damp::pow(2.0F, 3) == 8.0F); // integer exponent stays exact
        static_assert(damp::pow(-2.0, 3) == -8.0);
        static_assert(damp::pow(2.0, std::size_t{10}) == 1024.0);
    }

    TEST_CASE("integer pow is a constant expression up to and past the finite range") {
        static_assert(damp::pow(1e100, 2) == 1e200); // the discarded last square used to overflow
        static_assert(damp::pow(1.0e5F, 4) == 1e20F);
        static_assert(damp::pow(1e200, 2) == std::numeric_limits<double>::max());
        static_assert(damp::pow(-1e200, 3) == -std::numeric_limits<double>::max());
        static_assert(damp::pow(10.0, -310) == 0.0);
        static_assert(damp::pow(1e-200, -2) == std::numeric_limits<double>::max());
        // |base|^|n| subnormal: 1/p would round past max (no FTZ in constant evaluation).
        static_assert(damp::pow(1e-160, -2) == std::numeric_limits<double>::max());
        static_assert(damp::pow(1e-20F, -2) == std::numeric_limits<float>::max());
        static_assert(damp::pow(1e-10F, -4) == std::numeric_limits<float>::max());
        static_assert(damp::pow(-1e-160, -3) == -std::numeric_limits<double>::max());
        static_assert(damp::pow(2.0, -3) == 0.125);
        CHECK(damp::pow(1e200, 2) == std::numeric_limits<double>::max());
        volatile double tiny = 1e-160;
        volatile float  tiny_f = 1e-20F;
        CHECK(damp::pow(static_cast<double>(tiny), -2) == std::numeric_limits<double>::max());
        CHECK(damp::pow(static_cast<float>(tiny_f), -2) == std::numeric_limits<float>::max());
    }

    TEST_CASE("float series reduction uses float splits: sin/cos/tan/exp sweeps") {
        const auto ulp = [](double v) {
            // Next float by bits, difference in double: under -ffast-math nextafter is unreliable
            // and a float difference below FLT_MIN flushes to zero.
            const float f = static_cast<float>(std::fabs(v));
            return static_cast<double>(std::bit_cast<float>(std::bit_cast<std::uint32_t>(f) + 1U)) - static_cast<double>(f);
        };
        double worst = 0.0;     // in ulp of the result, |x| ≤ 100
        double worst_abs = 0.0; // absolute, |x| ≤ 1e5
        for (int i = -200000; i <= 200000; ++i) {
            const float  x = static_cast<float>(i) * 5e-4F;
            const double s = std::sin(static_cast<double>(x));
            const double c = std::cos(static_cast<double>(x));
            worst = std::fmax(worst, std::fabs(static_cast<double>(detail::sin(x)) - s) / ulp(s));
            worst = std::fmax(worst, std::fabs(static_cast<double>(detail::cos(x)) - c) / ulp(c));
            const float y = static_cast<float>(i) * 0.5F;
            worst_abs = std::fmax(worst_abs, std::fabs(static_cast<double>(detail::sin(y)) - std::sin(static_cast<double>(y))));
        }
        CHECK(worst <= 4.0);
        CHECK(worst_abs <= 1.2e-7);
        double worst_tan = 0.0;
        double worst_exp = 0.0;
        for (int i = -100000; i <= 100000; ++i) {
            const float  t = static_cast<float>(i) * 1.55e-5F;
            const double rt = std::tan(static_cast<double>(t));
            worst_tan = std::fmax(worst_tan, std::fabs(static_cast<double>(detail::tan(t)) - rt) / ulp(rt));
            const float  x = (static_cast<float>(i) * 8e-4F) + 4.0F; // [-76, 84]: results stay normal (fast-math flushes subnormals)
            const double re = std::exp(static_cast<double>(x));
            worst_exp = std::fmax(worst_exp, std::fabs(static_cast<double>(detail::exp(x)) - re) / ulp(re));
        }
        CHECK(worst_tan <= 4.0);
        CHECK(worst_exp <= 4.0);
        // Relative accuracy next to the zeros and poles (the sign was wrong before).
        constexpr float c = damp::cos(1.57079637F);
        static_assert(c < 0.0F);
        CHECK(c == doctest::Approx(std::cos(static_cast<double>(1.57079637F))).epsilon(1e-6));
        constexpr float t = damp::tan(1.57079637F);
        static_assert(t < 0.0F);
        CHECK(t == doctest::Approx(std::tan(static_cast<double>(1.57079637F))).epsilon(1e-6));
        // The largest float with a finite exp is not saturated, and is a constant expression.
        constexpr float e = damp::exp(88.72283172607421875F);
        static_assert(e < std::numeric_limits<float>::max());
        CHECK(e == doctest::Approx(std::exp(88.72283172607421875)).epsilon(1e-6));
    }

    TEST_CASE("series sin/cos stay in [-1, 1] and tan stays finite for every finite argument") {
        for (float x = 1e5F; x < 3e38F; x *= 1.37F) {
            for (const float a : {x, -x}) {
                const auto sc = detail::sincos(a);
                CHECK(std::fabs(sc.sin) <= 1.0F);
                CHECK(std::fabs(sc.cos) <= 1.0F);
                CHECK(damp::isfinite(detail::tan(a)));
            }
        }
        for (double x = 1e5; x < 1e308; x *= 13.7) {
            const auto sc = detail::sincos(x);
            CHECK(std::fabs(sc.sin) <= 1.0);
            CHECK(std::fabs(sc.cos) <= 1.0);
        }
        constexpr double s = damp::sin(1e20); // was not a constant expression
        static_assert(s <= 1.0 && s >= -1.0);
        constexpr float sf = damp::sin(1e30F);
        static_assert(sf <= 1.0F && sf >= -1.0F);
    }

    TEST_CASE("integer pow handles INT_MIN without signed overflow") {
        static_assert(damp::pow(1.0, std::numeric_limits<int>::min()) == 1.0);
        static_assert(damp::pow(-1.0, std::numeric_limits<int>::min()) == 1.0);
        CHECK(damp::pow(2.0, std::numeric_limits<int>::min()) == 0.0);
    }

    TEST_CASE("division does not overflow or underflow in the denominator") {
        const C big = C{1e200, 1e200} / C{1e200, 1e200};
        CHECK(big.real() == doctest::Approx(1.0));
        CHECK(std::fabs(big.imag()) <= kEps);
        const C small = C{1e-200, 2e-200} / C{1e-200, 0.0};
        CHECK(small.real() == doctest::Approx(1.0));
        CHECK(small.imag() == doctest::Approx(2.0));
        C q{3e200, 4e200};
        q /= C{0.0, 1e200};
        CHECK(q.real() == doctest::Approx(4.0));
        CHECK(q.imag() == doctest::Approx(-3.0));
        const C r = 1e-300 / C{1e-300, 1e-300};
        CHECK(r.real() == doctest::Approx(0.5));
        CHECK(r.imag() == doctest::Approx(-0.5));
    }

    TEST_CASE("abs does not overflow") {
        CHECK(C{3e200, 4e200}.abs() == doctest::Approx(5e200));
        CHECK(damp::abs(C{3e-200, 4e-200}) == doctest::Approx(5e-200));
    }

    TEST_CASE("sqrt keeps the small component (no cancellation)") {
        const C s = damp::sqrt(C{1.0, 1e-10});
        CHECK(s.real() == doctest::Approx(1.0));
        CHECK(s.imag() == doctest::Approx(5e-11));
        const C t = damp::sqrt(C{-1.0, -1e-10});
        CHECK(t.real() == doctest::Approx(5e-11));
        CHECK(t.imag() == doctest::Approx(-1.0));
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

    TEST_CASE("mag2db of zero or a negative magnitude is the finite minus infinity") {
        static_assert(damp::mag2db(0.0) == -std::numeric_limits<double>::max());
        static_assert(damp::mag2db(-1.0F) == -std::numeric_limits<float>::max());
        CHECK(damp::mag2db(0.0) == -std::numeric_limits<double>::max());
        CHECK(damp::mag2db(1.0) == 0.0);
    }
}