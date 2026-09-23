// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>

#include "damp/backend.hpp"
#include "damp/motor/estimators/encoder_tracker.hpp"
#include "damp/toolbox/encoder.hpp"
#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using namespace damp::motor;

namespace {
// Quadrature cycles in (A, B), each starting and ending at state 00 so repeated
// cycles chain cleanly from the decoder's initial 00 state.
// Forward: 00 -> 01 -> 11 -> 10 -> 00.
constexpr std::array<damp::pair<bool, bool>, 4> kForward{{{false, true}, {true, true}, {true, false}, {false, false}}};
// Reverse: 00 -> 10 -> 11 -> 01 -> 00.
constexpr std::array<damp::pair<bool, bool>, 4> kReverse{{{true, false}, {true, true}, {false, true}, {false, false}}};

template<typename Dec>
constexpr void drive(Dec& d, const std::array<damp::pair<bool, bool>, 4>& seq, int cycles) {
    for (int c = 0; c < cycles; ++c) {
        for (const auto& s : seq) {
            d.update(s.first, s.second);
        }
    }
}
} // namespace

TEST_SUITE("Quadrature encoder & tachometer") {
    TEST_CASE("X4 counts every edge, both directions") {
        QuadratureDecoder fwd{QuadMode::X4};
        drive(fwd, kForward, 3);
        CHECK(fwd.position() == 12); // 4 edges * 3 cycles

        QuadratureDecoder rev{QuadMode::X4};
        drive(rev, kReverse, 3);
        CHECK(rev.position() == -12);
    }

    TEST_CASE("X1 counts one detent per cycle") {
        QuadratureDecoder x1{QuadMode::X1};
        drive(x1, kForward, 5);
        CHECK(x1.position() == 5);
    }

    TEST_CASE("X2 counts two edges per cycle") {
        QuadratureDecoder x2{QuadMode::X2};
        drive(x2, kForward, 5);
        CHECK(x2.position() == 10);
    }

    TEST_CASE("illegal (double-bit) transition is ignored") {
        QuadratureDecoder d{QuadMode::X4};
        d.update(false, false); // state 00
        d.update(true, true);   // 00 -> 11 : illegal, no count
        CHECK(d.position() == 0);
    }

    TEST_CASE("index resets position to zero") {
        QuadratureDecoder d{QuadMode::X4};
        drive(d, kForward, 2);
        CHECK(d.position() == 8);
        d.index(true);
        CHECK(d.position() == 0);
    }

    TEST_CASE("wrapped_delta handles counter rollover") {
        CHECK(wrapped_delta<uint16_t>(0xFFFE, 0x0002) == 4);  // wrap up
        CHECK(wrapped_delta<uint16_t>(0x0002, 0xFFFE) == -4); // wrap down
        CHECK(wrapped_delta<uint16_t>(100, 150) == 50);
    }

    TEST_CASE("Tachometer frequency and period methods agree") {
        Tachometer<double> tach{1000.0}; // 1000 counts/rev
        // 1000 counts in 1 s = 1 rev/s.
        CHECK(tach.rev_per_s_from_counts(1000, 1.0) == doctest::Approx(1.0));
        // At 1 rev/s with 1000 cpr, edges are 1 ms apart.
        CHECK(tach.rev_per_s_from_period(0.001) == doctest::Approx(1.0));
        // Conversions.
        CHECK(Tachometer<double>::to_rpm(1.0) == doctest::Approx(60.0));
        CHECK(Tachometer<double>::to_rad_per_s(1.0) == doctest::Approx(2.0 * std::numbers::pi_v<double>));
        // Guard against divide-by-zero.
        CHECK(tach.rev_per_s_from_counts(5, 0.0) == doctest::Approx(0.0));
    }

    TEST_CASE("AngleUnwrapper accumulates continuously across the wrap seam") {
        AngleUnwrapper<double> unwrap{1.0};                // turns
        CHECK(unwrap.update(0.9) == doctest::Approx(0.9)); // first sample seeds, no step
        CHECK(unwrap.update(0.1) == doctest::Approx(1.1)); // +0.2 forward across 1.0 seam, not −0.8
        CHECK(unwrap.update(0.9) == doctest::Approx(0.9)); // −0.2 reverse across the seam
        CHECK(unwrap.update(0.8) == doctest::Approx(0.8)); // plain −0.1, no wrap
    }

    TEST_CASE("AngleUnwrapper works in radians and counts down multiple turns") {
        AngleUnwrapper<double> unwrap{2.0 * std::numbers::pi_v<double>}; // radians
        const double           step = 2.0;                               // rad/sample, under half-period (π)
        double                 pos = unwrap.update(0.0);
        for (int i = 1; i <= 10; ++i) {
            pos = unwrap.update(std::fmod(static_cast<double>(i) * step, 2.0 * std::numbers::pi_v<double>));
        }
        CHECK(pos == doctest::Approx(10.0 * step)); // 20 rad continuous, several wraps
    }
}

TEST_SUITE("EncoderTracker") {
    TEST_CASE("ctor seeds omega from x0 like reset") {
        const ColVec<2, double> x0{{0.25}, {12.5}}; // [turns, turns/s]
        EncoderTracker<double>  obs{2000.0, x0};
        CHECK(obs.theta() == doctest::Approx(0.25));
        CHECK(obs.omega() == doctest::Approx(12.5));
        CHECK(obs.state()[0] == doctest::Approx(obs.theta()));
        CHECK(obs.state()[1] == doctest::Approx(obs.omega()));

        // Default x0 leaves omega at 0
        EncoderTracker<double> zeroed{1000.0};
        CHECK(zeroed.omega() == doctest::Approx(0.0));

        // Config path and reset() agree with ctor seed
        EncoderTrackerConfig<double> cfg{.bandwidth = 1500.0, .x0 = ColVec<2, double>{{-0.1}, {-3.0}}};
        EncoderTracker<double>       from_cfg{cfg};
        CHECK(from_cfg.theta() == doctest::Approx(-0.1));
        CHECK(from_cfg.omega() == doctest::Approx(-3.0));

        zeroed.reset(x0);
        CHECK(zeroed.theta() == doctest::Approx(0.25));
        CHECK(zeroed.omega() == doctest::Approx(12.5));
    }

    TEST_CASE("predict advances theta with seeded omega") {
        EncoderTracker<double> obs{2000.0, ColVec<2, double>{{0.0}, {10.0}}}; // 10 turns/s
        const auto             x = obs.predict(0.01);                         // 0.1 turns
        CHECK(x[0] == doctest::Approx(0.1));
        CHECK(x[1] == doctest::Approx(10.0));
        CHECK(obs.theta() == doctest::Approx(0.1));
    }

    TEST_CASE("update tracks constant speed and multi-turn position") {
        // Critically-damped tracker: predict (free run) then update (encoder sample).
        EncoderTracker<double> obs{100.0};
        const double           dt = 0.001;
        const double           omega = 0.5; // turns/s
        double                 true_pos = 0.0;

        for (int k = 0; k < 3000; ++k) { // 3 s → 1.5 turns
            true_pos += omega * dt;
            const double meas = true_pos - std::floor(true_pos + 0.5);
            (void)obs.predict(dt);
            (void)obs.update(meas, dt);
        }
        CHECK(obs.omega() == doctest::Approx(omega).epsilon(0.05));
        CHECK(obs.position() == doctest::Approx(true_pos).epsilon(0.02));
        // Crossed at least one ±0.5 seam → multi-turn position > 1 turn
        CHECK(obs.position() > 1.0);
    }
}
