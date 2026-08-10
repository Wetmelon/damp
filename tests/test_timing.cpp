// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "damp/toolbox/timing.hpp"
#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

// Non-blocking software timers.

TEST_SUITE("Timing") {
    TEST_CASE("Stopwatch accumulates and resets") {
        Stopwatch<double> sw;
        CHECK(sw.elapsed() == doctest::Approx(0.0));
        sw.tick(0.1);
        sw.tick(0.2);
        CHECK(sw.elapsed() == doctest::Approx(0.3));
        sw.reset();
        CHECK(sw.elapsed() == doctest::Approx(0.0));
        sw.reset(5.0);
        CHECK(sw.elapsed() == doctest::Approx(5.0));
    }

    TEST_CASE("Timeout expires once and clamps its accumulator") {
        Timeout<double> to{0.25};
        CHECK_FALSE(to.tick(0.1)); // 0.1
        CHECK_FALSE(to.tick(0.1)); // 0.2
        CHECK(to.remaining() == doctest::Approx(0.05));
        CHECK(to.tick(0.1)); // 0.3 -> expired
        CHECK(to.expired());
        CHECK(to.remaining() == doctest::Approx(0.0));
        // Bounded: feeding more time doesn't grow the accumulator past duration.
        for (int i = 0; i < 1000; ++i) {
            to.tick(1.0);
        }
        CHECK(to.remaining() == doctest::Approx(0.0));
        to.reset();
        CHECK_FALSE(to.expired());
    }

    TEST_CASE("Periodic fires once per period without drift") {
        Periodic<double> p{0.30};
        const double     dt = 0.10;
        int              fires = 0;
        for (int i = 0; i < 30; ++i) { // 3.0 s total
            if (p(dt)) {
                ++fires;
            }
        }
        // 3.0 s / 0.3 s = 10 firings, no drift accumulation.
        CHECK(fires == 10);
    }

    TEST_CASE("Periodic catches up when dt spans multiple periods") {
        Periodic<double> p{0.01};
        // A single big step is one fire; remainder carries so the next small
        // ticks fire again.
        CHECK(p(0.025)); // phase 0.025 -> fire, phase 0.015
        CHECK(p(0.0));   // phase 0.015 -> fire, phase 0.005
        CHECK_FALSE(p(0.0));
    }

    TEST_CASE("Periodic with period <= 0 never fires") {
        Periodic<double>   zero{0.0};
        Periodic<double>   neg{-1.0};
        Periodic<uint32_t> zero_u{0u};
        for (int i = 0; i < 100; ++i) {
            CHECK_FALSE(zero(0.1));
            CHECK_FALSE(neg(0.1));
            CHECK_FALSE(zero_u(1u));
        }
        // Default-constructed period is 0 — also never fires.
        Periodic<> def;
        CHECK_FALSE(def(1u));
        CHECK_FALSE(def(100u));
    }
}
