// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <initializer_list>

#include "damp/filters/filters.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

/**
 * @brief Tests for the simple runtime filters added to filters.hpp:
 *        HighPass (washout), MedianFilter (despiker), Complementary (scalar fusion).
 */

TEST_SUITE("Simple filters") {

    TEST_CASE("HighPass rejects DC exactly, passes a step transient") {
        HighPass<double> hp{2.0, 1.0 / 1000.0}; // 2 Hz corner @ 1 kHz (Tustin highpass_1st)

        // Feed a constant; after the first sample the washout must decay to ~0.
        double y = 0.0;
        for (int i = 0; i < 20000; ++i) {
            y = hp(5.0);
        }
        CHECK(y == doctest::Approx(0.0).epsilon(1e-6)); // DC fully removed

        // A fresh step: first sample is b0 (Tustin HF gain ≈ 1 for fc ≪ Nyquist).
        hp.reset();
        const auto   c = design::highpass_1st<double>(2.0, 1.0 / 1000.0);
        const double first = hp(1.0);
        CHECK(first == doctest::Approx(c.b0));
        CHECK(first == doctest::Approx(1.0).epsilon(0.02)); // near-unity HF gain
    }

    TEST_CASE("HighPass invalid Ts/fc yields zero output") {
        HighPass<double> bad_ts{10.0, 0.0};
        HighPass<double> bad_fc{0.0, 0.001};
        CHECK(bad_ts(1.0) == doctest::Approx(0.0));
        CHECK(bad_fc(1.0) == doctest::Approx(0.0));
    }

    TEST_CASE("MedianFilter rejects a single-sample spike") {
        MedianFilter<5, double> med{5};
        for (double v : {1.0, 1.0, 1.0, 1.0}) {
            (void)med(v);
        }
        // Inject a large spike; median of {1,1,1,1,100} = 1 (spike rejected).
        CHECK(med(100.0) == doctest::Approx(1.0));
        // A moving average would have jumped to ~20.8 here.
    }

    TEST_CASE("MedianFilter even window averages the two middle samples") {
        MedianFilter<4, double> med{4};
        (void)med(1.0);
        (void)med(2.0);
        (void)med(3.0);
        // window {1,2,3,4} -> mean(2,3) = 2.5
        CHECK(med(4.0) == doctest::Approx(2.5));
    }

    TEST_CASE("Complementary seeds to measurement then blends") {
        Complementary<double> cf{0.5}; // tau = 0.5 s
        const double          dt = 0.01;

        // First call seeds to the measurement exactly.
        CHECK(cf(10.0, 0.0, dt) == doctest::Approx(10.0));

        // With zero rate and a constant measurement it must converge to it.
        double y = 10.0;
        for (int i = 0; i < 100000; ++i) {
            y = cf(2.0, 0.0, dt);
        }
        CHECK(y == doctest::Approx(2.0).epsilon(1e-6));

        // Pure rate integration (measurement == current value) advances by rate*dt.
        cf.reset();
        (void)cf(0.0, 0.0, dt); // seed at 0
        const double after = cf(0.0, 1.0, dt);
        CHECK(after > 0.0); // gyro term moved it positive
    }
}
