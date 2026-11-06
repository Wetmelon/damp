// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstdlib>

#include "damp/toolbox/conditioning.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

/**
 * @brief Tests for nonlinear signal-conditioning blocks (conditioning.hpp):
 *        deadband, inverse_deadband, SlewLimiter, Hysteresis.
 */

TEST_SUITE("Signal conditioning") {

    TEST_CASE("deadband matches Simulink® Dead Zone (symmetric and [lower, upper])") {
        // Symmetric overload: [-0.5, 0.5].
        CHECK(deadband(0.3, 0.5) == doctest::Approx(0.0));       // inside
        CHECK(deadband(-0.5, 0.5) == doctest::Approx(0.0));      // on edge
        CHECK(deadband(0.6, 0.5) == doctest::Approx(0.1));       // x - upper
        CHECK(deadband(-0.6, 0.5) == doctest::Approx(-0.1));     // x - lower
        CHECK(deadband(0.5001, 0.5) == doctest::Approx(0.0001)); // continuous at edge

        // Asymmetric [lower, upper] like the Dead Zone block.
        CHECK(deadband(2.0, 1.5, 2.5) == doctest::Approx(0.0));  // inside
        CHECK(deadband(2.7, 1.5, 2.5) == doctest::Approx(0.2));  // x - upper
        CHECK(deadband(1.0, 1.5, 2.5) == doctest::Approx(-0.5)); // x - lower (below band)
    }

    TEST_CASE("inverse_deadband boosts past a physical dead zone (symmetric + asymmetric)") {
        // Symmetric overload.
        CHECK(inverse_deadband(0.0, 0.2) == doctest::Approx(0.0)); // exact zero held
        CHECK(inverse_deadband(1.0, 0.2) == doctest::Approx(1.2)); // +band
        CHECK(inverse_deadband(-1.0, 0.2) == doctest::Approx(-1.2));

        // Asymmetric [lower, upper] with independent offsets + threshold.
        CHECK(inverse_deadband(0.5, -0.1, 0.3) == doctest::Approx(0.8));       // x>0 -> +upper
        CHECK(inverse_deadband(-0.5, -0.1, 0.3) == doctest::Approx(-0.6));     // x<0 -> +lower
        CHECK(inverse_deadband(0.05, -0.1, 0.3, 0.1) == doctest::Approx(0.0)); // within threshold
    }

    TEST_CASE("soft_sign is odd, unit-bounded, and approaches sgn") {
        CHECK(soft_sign(0.0, 0.1) == doctest::Approx(0.0));
        CHECK(soft_sign(10.0, 0.1) == doctest::Approx(1.0).epsilon(1e-4));
        CHECK(soft_sign(-10.0, 0.1) == doctest::Approx(-1.0).epsilon(1e-4));
        CHECK(soft_sign(0.1, 0.1) == doctest::Approx(0.1 / std::sqrt(0.02)).epsilon(1e-12));
        CHECK(soft_sign(-0.3, 0.05) == doctest::Approx(-soft_sign(0.3, 0.05)));
    }

    TEST_CASE("inverse_deadband_soft approaches hard inverse for small eps") {
        CHECK(inverse_deadband_soft(0.0, 0.2, 1e-9) == doctest::Approx(0.0));
        CHECK(inverse_deadband_soft(1.0, 0.2, 1e-6) == doctest::Approx(1.2).epsilon(1e-5));
        CHECK(inverse_deadband_soft(-1.0, 0.2, 1e-6) == doctest::Approx(-1.2).epsilon(1e-5));
        // Softer near zero than hard step
        const double soft = inverse_deadband_soft(0.05, 0.2, 0.1);
        const double hard = inverse_deadband(0.05, 0.2);
        CHECK(std::abs(soft) < std::abs(hard));
        CHECK(inverse_deadband_soft(0.01, 0.2, 0.1, 0.05) == doctest::Approx(0.0)); // threshold
    }

    TEST_CASE("inverse_deadband_sine ramps to full band at width") {
        CHECK(inverse_deadband_sine(0.0, 0.3, 0.2) == doctest::Approx(0.0));
        // |x| >= width → full ±band offset
        CHECK(inverse_deadband_sine(0.5, 0.3, 0.2) == doctest::Approx(0.8).epsilon(1e-12));
        CHECK(inverse_deadband_sine(-0.5, 0.3, 0.2) == doctest::Approx(-0.8).epsilon(1e-12));
        // Mid ramp: sin(π/4)=√2/2 → offset = 0.3 * √2/2
        const double mid = inverse_deadband_sine(0.1, 0.3, 0.2); // s = 0.5
        CHECK(mid == doctest::Approx(0.1 + 0.3 * std::sin(0.5 * 3.141592653589793 * 0.5)).epsilon(1e-9));
        CHECK(inverse_deadband_sine(0.01, 0.3, 0.2, 0.05) == doctest::Approx(0.0));
    }

    TEST_CASE("scaled_deadband rescales surviving range to full span") {
        CHECK(scaled_deadband(0.05, 0.1) == doctest::Approx(0.0)); // inside dead zone
        CHECK(scaled_deadband(1.0, 0.1) == doctest::Approx(1.0));  // extreme stays at 1
        CHECK(scaled_deadband(-1.0, 0.1) == doctest::Approx(-1.0));
        CHECK(scaled_deadband(0.1, 0.1) == doctest::Approx(0.0)); // on edge
        // Just past the edge starts from 0: (0.55-0.1)/(1-0.1) = 0.5.
        CHECK(scaled_deadband(0.55, 0.1) == doctest::Approx(0.5));
    }

    TEST_CASE("expo softens center, preserves endpoints and sign") {
        CHECK(expo(0.0, 0.5) == doctest::Approx(0.0));
        CHECK(expo(1.0, 0.5) == doctest::Approx(1.0)); // endpoint preserved
        CHECK(expo(-1.0, 0.5) == doctest::Approx(-1.0));
        CHECK(expo(0.5, 0.0) == doctest::Approx(0.5));   // k=0 linear
        CHECK(expo(0.5, 1.0) == doctest::Approx(0.125)); // k=1 cubic: 0.5^3
        // Mid expo: (1-0.5)*0.5 + 0.5*0.125 = 0.3125, softer than linear.
        CHECK(expo(0.5, 0.5) == doctest::Approx(0.3125));
    }

    TEST_CASE("SlewLimiter caps the per-step change, asymmetric up/down") {
        SlewLimiter<double> sl{10.0, 100.0}; // up 10/s, down 100/s
        const double        dt = 0.1;        // max +1.0 / −10.0 per step

        CHECK(sl(5.0, dt) == doctest::Approx(5.0));      // first sample seeds to target
        CHECK(sl(100.0, dt) == doctest::Approx(6.0));    // limited rise: +1.0
        CHECK(sl(100.0, dt) == doctest::Approx(7.0));    // +1.0 again
        CHECK(sl(-100.0, dt) == doctest::Approx(-3.0));  // fast fall: −10.0 (7 → −3)
        CHECK(sl(-100.0, dt) == doctest::Approx(-13.0)); // −10.0 again
        // A small move within the step limit snaps straight to target.
        CHECK(sl(-13.5, dt) == doctest::Approx(-13.5));

        // Small moves within the limit pass through exactly.
        SlewLimiter<double> s2{1000.0};
        (void)s2(0.0, dt);
        CHECK(s2(0.5, dt) == doctest::Approx(0.5));
    }

    TEST_CASE("SlewLimiter holds output when dt <= 0") {
        SlewLimiter<double> sl{10.0};
        CHECK(sl(0.0, 0.1) == doctest::Approx(0.0));  // seed
        CHECK(sl(10.0, 0.1) == doctest::Approx(1.0)); // rise by rate*dt
        // Non-positive dt must not reverse or jump — hold last output.
        CHECK(sl(100.0, 0.0) == doctest::Approx(1.0));
        CHECK(sl(100.0, -0.05) == doctest::Approx(1.0));
        CHECK(sl.value() == doctest::Approx(1.0));
        // Positive dt resumes ramping.
        CHECK(sl(100.0, 0.1) == doctest::Approx(2.0));
    }

    TEST_CASE("classify_range / RangeMonitor bands per NE43 [0.25 (0.5 4.5) 4.75]") {
        // Pure kernel; a braced range takes its scalar from x.
        const SignalRange<double> ne43{.valid = {
                                           .lower = 0.5,
                                           .upper = 4.5,
                                       },
                                       .fault = {
                                           .lower = 0.25,
                                           .upper = 4.75,
                                       }};
        CHECK(classify_range(0.1, {.valid = {
                                       .lower = 0.5,
                                       .upper = 4.5,
                                   },
                                   .fault = {
                                       .lower = 0.25,
                                       .upper = 4.75,
                                   }})
              == SignalStatus::FaultLow);
        CHECK(classify_range(0.3, ne43) == SignalStatus::UnderRange);
        CHECK(classify_range(2.5, ne43) == SignalStatus::Valid);
        CHECK(classify_range(4.6, ne43) == SignalStatus::OverRange);
        CHECK(classify_range(4.9, ne43) == SignalStatus::FaultHigh);
        // Band edges: valid span is inclusive, fault floor belongs to UnderRange.
        CHECK(classify_range(0.5, ne43) == SignalStatus::Valid);
        CHECK(classify_range(0.25, ne43) == SignalStatus::UnderRange);

        RangeMonitor<double> ai{ne43};
        CHECK(ai.faulted(0.1) == true);  // broken wire
        CHECK(ai.faulted(0.3) == false); // saturated low is NOT a wire fault
        CHECK(ai.valid(2.5) == true);
        CHECK(is_fault(ai(4.9)) == true);
    }

    TEST_CASE("RangeMonitor hysteresis suppresses boundary chatter") {
        RangeMonitor<double> ai{{.valid = {
                                     .lower = 0.5,
                                     .upper = 4.5,
                                 },
                                 .fault = {
                                     .lower = 0.25,
                                     .upper = 4.75,
                                 }},
                                0.05};
        CHECK(ai(2.5) == SignalStatus::Valid);
        // Just below valid.lower: without hysteresis this would flip to UnderRange,
        // but it must drop a further 0.05 (below 0.45) to leave Valid.
        CHECK(ai(0.48) == SignalStatus::Valid);
        CHECK(ai(0.44) == SignalStatus::UnderRange); // crossed valid.lower - hyst
        // And to recover it must rise back above valid.lower + hyst = 0.55.
        CHECK(ai(0.52) == SignalStatus::UnderRange);
        CHECK(ai(0.56) == SignalStatus::Valid);
    }

    TEST_CASE("Hysteresis latches with separate trip/release thresholds") {
        Hysteresis<double> h{1.0, 2.0}; // release < 1, trip > 2

        CHECK(h(0.0) == false);
        CHECK(h(1.5) == false); // in the band, still low
        CHECK(h(2.5) == true);  // trips above high
        CHECK(h(1.5) == true);  // holds in the band
        CHECK(h(0.5) == false); // releases below low
        CHECK(h(1.5) == false); // holds low again
    }

    TEST_CASE("soft_sign does not collapse to 0 when x² overflows") {
        CHECK(soft_sign(1e20f, 1.0f) == doctest::Approx(1.0f));
        CHECK(soft_sign(-1e200, 1.0) == doctest::Approx(-1.0));
        CHECK(soft_sign(0.0, 0.0) == 0.0);
        CHECK(soft_sign(1e-30, 0.0) == doctest::Approx(1.0)); // eps = 0 is the hard sign
    }
}
