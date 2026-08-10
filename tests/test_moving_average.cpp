// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "damp/filters/filters.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {
constexpr double pi = 3.14159265358979323846;
} // namespace

TEST_SUITE("moving_average") {
    TEST_CASE("comb_notch_window computes N = round(fs / f_notch)") {
        CHECK(comb_notch_window(1200.0, 120.0) == 10); // HV bus: 2f ripple at 120 Hz
        CHECK(comb_notch_window(1000.0, 47.0) == 21);  // 21.28 -> 21
        CHECK(comb_notch_window(1000.0, 1000.0) == 0); // f >= fs invalid
        CHECK(comb_notch_window(1000.0, -50.0) == 0);  // invalid
    }

    TEST_CASE("passes DC unchanged") {
        MovingAverage<32, double> ma(8);
        double                    y = 0.0;
        for (int k = 0; k < 50; ++k) {
            y = ma(5.0);
        }
        CHECK(y == doctest::Approx(5.0)); // constant in -> constant out
    }

    TEST_CASE("comb notch: strips a ripple and its harmonics, keeps the DC (HV-bus case)") {
        // fs = 1200 Hz, ripple fundamental 120 Hz (2f for a 60 Hz line) + 240 Hz.
        // N = fs / 120 = 10 nulls every harmonic of 120 Hz, unity at DC.
        const double fs = 1200.0;
        const double f_ripple = 120.0;
        const size_t N = comb_notch_window(fs, f_ripple);
        REQUIRE(N == 10);

        const double dc = 100.0;
        auto         bus = [&](size_t k) {
            const double t = static_cast<double>(k) / fs;
            return dc + 5.0 * std::sin(2.0 * pi * 120.0 * t) // ripple fundamental
                 + 2.0 * std::sin(2.0 * pi * 240.0 * t)      // 2nd ripple harmonic
                 + 1.0 * std::cos(2.0 * pi * 360.0 * t);     // 3rd
        };

        MovingAverage<32, double> ma(N);
        double                    worst = 0.0;
        for (size_t k = 0; k < 200; ++k) {
            const double y = ma(bus(k));
            if (k >= N) { // past warm-up: ripple fully averaged out
                worst = std::max(worst, std::abs(y - dc));
            }
        }
        CHECK(worst < 1e-9); // DC preserved, every ripple harmonic nulled
    }

    TEST_CASE("acts as a boxcar smoother (mean of the window)") {
        MovingAverage<8, double> ma(4);
        // Push 1,2,3,4 -> average of last 4 = 2.5.
        ma(1.0);
        ma(2.0);
        ma(3.0);
        const double y = ma(4.0);
        CHECK(y == doctest::Approx((1.0 + 2.0 + 3.0 + 4.0) / 4.0));
    }

    TEST_CASE("set_window and reset") {
        MovingAverage<16, double> ma;
        ma.set_window(4);
        CHECK(ma.window() == 4);
        ma.set_window(100); // clamped to MaxN
        CHECK(ma.window() == 16);

        ma.set_window(4);
        ma(10.0);
        ma.reset();
        // After reset the window is empty: first sample averages against zeros.
        CHECK(ma(8.0) == doctest::Approx(8.0 / 4.0));
    }

    TEST_CASE("MovingAverage is constexpr-evaluable") {
        constexpr double y = []() consteval {
            MovingAverage<16, double> ma(4);
            double                    out = 0.0;
            for (int k = 0; k < 20; ++k) {
                out = ma(3.0);
            }
            return out;
        }();
        static_assert(y > 2.999 && y < 3.001, "MovingAverage must work at compile time");
        CHECK(y == doctest::Approx(3.0));
    }
}
