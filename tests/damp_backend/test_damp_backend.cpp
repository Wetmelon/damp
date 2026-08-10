// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Exercises the DAMP_MATH_BACKEND_DAMP profile in isolation. The macro is set
// here rather than in tests/damp_profile.hpp on purpose: the Damp MathBackend<float>
// specialization is ODR-incompatible with the default std backend the rest of the
// suite links against, so it gets its own executable. This is the only place
// damp_backend.hpp (and the trig.hpp fast-float path it forwards to) is compiled
// and run — the guard that catches a future break like the truncated nearbyint.
#define DAMP_MATH_BACKEND_DAMP

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <cmath>

#include "damp/math/math.hpp"
#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

static constexpr float kEps = 1e-5f;

TEST_SUITE("damp_backend") {
    TEST_CASE("MathBackend<float> routes to the fast-float trig path") {
        CHECK(damp::MathBackend<float>::sin(1.0f) == doctest::Approx(std::sin(1.0f)).epsilon(kEps));
        CHECK(damp::MathBackend<float>::cos(1.0f) == doctest::Approx(std::cos(1.0f)).epsilon(kEps));
        CHECK(damp::MathBackend<float>::atan(0.7f) == doctest::Approx(std::atan(0.7f)).epsilon(kEps));
        CHECK(damp::MathBackend<float>::atan2(0.7f, -0.3f) == doctest::Approx(std::atan2(0.7f, -0.3f)).epsilon(kEps));
        CHECK(damp::MathBackend<float>::asin(0.4f) == doctest::Approx(std::asin(0.4f)).epsilon(kEps));
        CHECK(damp::MathBackend<float>::acos(0.4f) == doctest::Approx(std::acos(0.4f)).epsilon(kEps));
        CHECK(damp::MathBackend<float>::sqrt(2.0f) == doctest::Approx(std::sqrt(2.0f)).epsilon(kEps));

        auto [s, c] = damp::MathBackend<float>::sincos(0.9f);
        CHECK(s == doctest::Approx(std::sin(0.9f)).epsilon(kEps));
        CHECK(c == doctest::Approx(std::cos(0.9f)).epsilon(kEps));
    }

    TEST_CASE("MathBackend<float>::nearbyint rounds to even") {
        CHECK(damp::MathBackend<float>::nearbyint(2.5f) == 2.0f);
        CHECK(damp::MathBackend<float>::nearbyint(3.5f) == 4.0f);
        CHECK(damp::MathBackend<float>::nearbyint(-1.4f) == -1.0f);
        CHECK(damp::MathBackend<float>::nearbyint(100.0f) == 100.0f);
    }

    TEST_CASE("non-overridden float ops fall through to the std base") {
        CHECK(damp::MathBackend<float>::tan(0.5f) == doctest::Approx(std::tan(0.5f)).epsilon(kEps));
        CHECK(damp::MathBackend<float>::exp(1.0f) == doctest::Approx(std::exp(1.0f)).epsilon(kEps));
        CHECK(damp::MathBackend<float>::log(2.0f) == doctest::Approx(std::log(2.0f)).epsilon(kEps));
    }

    TEST_CASE("MathBackend<double> routes to the fast double trig path") {
        constexpr double kEpsD = 1e-12;
        CHECK(damp::MathBackend<double>::sin(1.0) == doctest::Approx(std::sin(1.0)).epsilon(kEpsD));
        CHECK(damp::MathBackend<double>::cos(1.0) == doctest::Approx(std::cos(1.0)).epsilon(kEpsD));
        CHECK(damp::MathBackend<double>::atan(0.7) == doctest::Approx(std::atan(0.7)).epsilon(1e-10));
        CHECK(damp::MathBackend<double>::asin(0.4) == doctest::Approx(std::asin(0.4)).epsilon(1e-10));
        CHECK(damp::sin(1.0) == doctest::Approx(std::sin(1.0)).epsilon(kEpsD));
        auto [s, c] = damp::MathBackend<double>::sincos(0.9);
        CHECK(s == doctest::Approx(std::sin(0.9)).epsilon(kEpsD));
        CHECK(c == doctest::Approx(std::cos(0.9)).epsilon(kEpsD));
        CHECK(damp::wrap(7.0, -3.14159265358979, 3.14159265358979) == doctest::Approx(std::atan2(std::sin(7.0), std::cos(7.0))));
    }
}
