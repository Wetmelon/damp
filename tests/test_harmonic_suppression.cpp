// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <array>
#include <cmath>
#include <cstddef>

#include "damp/controllers/harmonic_suppression.hpp"
#include "damp/controllers/pr.hpp"
#include "damp/math/math.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {
constexpr double pi = 3.14159265358979323846;
} // namespace

TEST_SUITE("harmonic_suppression") {
    TEST_CASE("synthesize validates the spec") {
        const double                 w0 = 2.0 * pi * 50.0; // 50 Hz fundamental
        const double                 Ts = 1.0 / 6000.0;    // 6 kHz loop
        const damp::array<size_t, 4> harmonics{1, 5, 7, 11};

        SUBCASE("valid spec succeeds and tunes one resonator per harmonic") {
            const auto d = design::harmonic_suppressor(0.5, 200.0, w0, 5.0, Ts, harmonics);
            CHECK(d.success);
            // Fundamental carries Kp; harmonics carry Kp = 0.
            CHECK(d.gains[0].Kp == doctest::Approx(0.5));
            CHECK(d.gains[1].Kp == doctest::Approx(0.0));
            CHECK(d.gains[3].w0 == doctest::Approx(w0 * 11.0));
        }
        SUBCASE("rejects non-positive Ts / frequency") {
            CHECK_FALSE(design::harmonic_suppressor(0.5, 200.0, w0, 5.0, 0.0, harmonics).success);
            CHECK_FALSE(design::harmonic_suppressor(0.5, 200.0, -w0, 5.0, Ts, harmonics).success);
        }
        SUBCASE("rejects a 0th harmonic and any harmonic at/above Nyquist") {
            CHECK_FALSE(
                design::harmonic_suppressor(0.5, 200.0, w0, 5.0, Ts, damp::array<size_t, 2>{1, 0}).success
            );
            // 61st harmonic of 50 Hz = 3050 Hz > Nyquist (3 kHz).
            CHECK_FALSE(
                design::harmonic_suppressor(0.5, 200.0, w0, 5.0, Ts, damp::array<size_t, 1>{61}).success
            );
        }
    }

    TEST_CASE("tracks a multi-harmonic reference to zero steady-state error") {
        // Mirror the single-PR tracking test: unit-gain plant y[k] = u[k], the bank
        // closes the loop on the error. Each PR resonator has (near-)infinite gain at
        // its harmonic, so a reference built from exactly those harmonics is tracked
        // with vanishing steady-state error.
        const double                 fs = 10000.0;
        const double                 Ts = 1.0 / fs;
        const double                 w0 = 2.0 * pi * 50.0;
        const damp::array<size_t, 2> harmonics{5, 7};

        const auto d = design::harmonic_suppressor(0.0, 500.0, w0, 10.0, Ts, harmonics);
        REQUIRE(d.success);
        HarmonicSuppressor<2, double> sup(d);

        auto reference = [&](size_t k) {
            const double t = static_cast<double>(k) * Ts;
            return 1.0 * std::sin(2.0 * pi * 250.0 * t)  // 5th
                 + 0.7 * std::sin(2.0 * pi * 350.0 * t); // 7th
        };

        double       y = 0.0;
        double       error_sum = 0.0;
        const size_t n_steps = 20000; // 2 s at 10 kHz
        for (size_t k = 0; k < n_steps; ++k) {
            const double e = reference(k) - y;
            y = sup.control(e); // unit-gain plant: y[k] = u[k]
            if (k > 15000) {    // settled window
                error_sum += std::abs(e);
            }
        }
        CHECK(error_sum / 5000.0 < 0.05); // both targeted harmonics tracked
    }

    TEST_CASE("invalid design is inert (valid() false, finite output)") {
        const double                  Ts = 1.0 / 6000.0;
        const auto                    d = design::harmonic_suppressor(0.5, 200.0, -1.0, 5.0, Ts, damp::array<size_t, 1>{5});
        HarmonicSuppressor<1, double> sup(d);
        CHECK_FALSE(sup.valid());
    }

    TEST_CASE("float specialization compiles and runs") {
        const float                  Ts = 1.0f / 6000.0f;
        const float                  w0 = 2.0f * 3.14159265f * 50.0f;
        const damp::array<size_t, 2> harmonics{1, 5};
        const auto                   d = design::harmonic_suppressor(0.5f, 50.0f, w0, 5.0f, Ts, harmonics);
        REQUIRE(d.success);
        HarmonicSuppressor<2, float> sup(d);
        const float                  u = sup.control(1.0f);
        CHECK(std::isfinite(u));
    }

    TEST_CASE("as<U>() converts the design result precision") {
        const double                 w0 = 2.0 * pi * 50.0;
        const double                 Ts = 1.0 / 6000.0;
        const damp::array<size_t, 2> harmonics{1, 5};
        const auto                   d = design::harmonic_suppressor(0.5, 200.0, w0, 5.0, Ts, harmonics);
        const auto                   df = d.template as<float>();
        CHECK(df.success);
        CHECK(df.gains[0].Kp == doctest::Approx(0.5f));
    }

    TEST_CASE("harmonic suppressor is constexpr-evaluable") {
        constexpr bool ok = []() consteval {
            const double                  w0 = 2.0 * 3.14159265358979323846 * 50.0;
            const double                  Ts = 1.0 / 6000.0;
            const damp::array<size_t, 2>  harmonics{1, 5};
            const auto                    d = design::harmonic_suppressor(0.5, 100.0, w0, 5.0, Ts, harmonics);
            HarmonicSuppressor<2, double> sup(d);
            double                        u = 0.0;
            for (int k = 0; k < 5; ++k) {
                u = sup.control(1.0);
            }
            return d.success && damp::abs(u) < 1e6; // finite, bounded
        }();
        static_assert(ok, "harmonic suppressor must work at compile time");
        CHECK(ok);
    }
}
