// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "damp/toolbox/logic.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

/**
 * @brief Tests for the modern discrete logic/timing blocks (logic.hpp).
 */

TEST_SUITE("Logic & timing blocks") {

    TEST_CASE("Edge detectors fire on the transition only") {
        RisingEdge re;
        CHECK(re(false) == false);
        CHECK(re(true) == true);  // rising
        CHECK(re(true) == false); // held high, no edge
        CHECK(re(false) == false);

        FallingEdge fe;
        CHECK(fe(true) == false);
        CHECK(fe(false) == true); // falling
        CHECK(fe(false) == false);
    }

    TEST_CASE("Latch dominance") {
        Latch<> sd;                      // set-dominant
        CHECK(sd(true, true) == true);   // set wins
        CHECK(sd(false, true) == false); // reset clears
        CHECK(sd(true, false) == true);  // set holds

        ResetDominantLatch rd;
        CHECK(rd(true, false) == true); // set
        CHECK(rd(true, true) == false); // reset wins
    }

    TEST_CASE("OnDelayTimer waits then latches; drops immediately") {
        OnDelayTimer<double> t{0.5};
        const double         dt = 0.1;
        for (int i = 0; i < 4; ++i) {
            CHECK(t(true, dt) == false); // 0.1..0.4 s, not yet
        }
        CHECK(t(true, dt) == true);   // 0.5 s reached
        CHECK(t(true, dt) == true);   // stays
        CHECK(t(false, dt) == false); // input drops -> output drops at once
    }

    TEST_CASE("OnDelayTimer with integral T counts ticks exactly (dt = 1)") {
        OnDelayTimer<int> t{3};
        CHECK(t(true, 1) == false); // 1
        CHECK(t(true, 1) == false); // 2
        CHECK(t(true, 1) == true);  // 3 -> fire
    }

    TEST_CASE("OffDelayTimer holds output after input drops") {
        OffDelayTimer<double> t{0.3};
        const double          dt = 0.1;
        CHECK(t(true, dt) == true);   // on immediately
        CHECK(t(false, dt) == true);  // 0.1 s into off-delay, still on
        CHECK(t(false, dt) == true);  // 0.2 s
        CHECK(t(false, dt) == false); // 0.3 s -> off
    }

    TEST_CASE("PulseTimer emits a fixed-width, non-retriggerable pulse") {
        PulseTimer<double> t{0.25};
        const double       dt = 0.1;
        CHECK(t(true, dt) == true);  // rising edge starts pulse
        CHECK(t(true, dt) == true);  // 0.2 s, still high (retrigger ignored)
        CHECK(t(true, dt) == false); // 0.3 s >= width -> done
    }

    TEST_CASE("Debounce ignores brief glitches") {
        Debounce<double> db{0.3}; // need 0.3 s stable to flip
        const double     dt = 0.1;
        CHECK(db(true, dt) == false);  // 0.1 s of "true", not yet
        CHECK(db(false, dt) == false); // glitch back to false resets the timer
        CHECK(db(true, dt) == false);  // 0.1
        CHECK(db(true, dt) == false);  // 0.2
        CHECK(db(true, dt) == true);   // 0.3 -> committed
    }

    TEST_CASE("Toggle flips on each rising edge") {
        Toggle tg;
        CHECK(tg(true) == true); // flip
        CHECK(tg(true) == true); // held, no flip
        CHECK(tg(false) == true);
        CHECK(tg(true) == false); // flip back
    }

    TEST_CASE("Counter counts up/down edges") {
        Counter<int> c;
        CHECK(c(true, false) == 1); // up edge
        CHECK(c(true, false) == 1); // held, no edge
        CHECK(c(false, false) == 1);
        CHECK(c(true, true) == 1); // simultaneous up+down edges net zero
        CHECK(c(false, false) == 1);
        CHECK(c(false, true) == 0); // down edge
    }
}
