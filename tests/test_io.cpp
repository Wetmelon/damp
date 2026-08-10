// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "damp/toolbox/conditioning.hpp"
#include "damp/toolbox/io.hpp"
#include "damp/toolbox/lookup.hpp"
#include "damp/toolbox/scaling.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

/**
 * @brief Tests for the tier-2 I/O appliances (io.hpp).
 */

TEST_SUITE("Tier-2 IO appliances") {

    TEST_CASE("AnalogInput validates raw then scales to engineering units") {
        // 0.5–4.5 V -> 0–100 %, with NE43 fault bands [0.25 (0.5 4.5) 4.75].
        io::AnalogInput<double> ai{
            two_point_cal(0.5, 0.0, 4.5, 100.0),
            RangeMonitor<double>{0.25, 0.5, 4.5, 4.75}
        };

        CHECK(ai.update(2.5) == doctest::Approx(50.0)); // mid-scale
        CHECK(ai.valid() == true);

        ai.update(0.1); // below fault floor
        CHECK(ai.faulted() == true);
        CHECK(ai.status() == SignalStatus::FaultLow);

        ai.update(0.3); // saturated low, not a wire fault
        CHECK(ai.faulted() == false);
        CHECK(ai.status() == SignalStatus::UnderRange);
    }

    TEST_CASE("cross_check_analog same-slope agrees and averages") {
        const auto r = cross_check_analog(
            100.0, SignalStatus::Valid, 102.0, SignalStatus::Valid, AnalogCrossMode::SameSlope, 5.0
        );
        CHECK(r.agree());
        CHECK(r.has_value);
        CHECK(r.value == doctest::Approx(101.0));
    }

    TEST_CASE("cross_check_analog same-slope disagree yields no value") {
        const auto r = cross_check_analog(
            100.0, SignalStatus::Valid, 120.0, SignalStatus::Valid, AnalogCrossMode::SameSlope, 5.0
        );
        CHECK(r.status == AnalogCrossStatus::Disagree);
        CHECK_FALSE(r.has_value);
    }

    TEST_CASE("cross_check_analog opposite-slope sums to span") {
        // Complementary 0–100%: a=30, b=70 → sum 100.
        const auto r = cross_check_analog(
            30.0, SignalStatus::Valid, 70.0, SignalStatus::Valid, AnalogCrossMode::OppositeSlope, 1.0, 100.0
        );
        CHECK(r.agree());
        CHECK(r.value == doctest::Approx(30.0)); // (30 + (100-70))/2
    }

    TEST_CASE("cross_check_analog opposite-slope detects mismatch") {
        const auto r = cross_check_analog(
            30.0, SignalStatus::Valid, 50.0, SignalStatus::Valid, AnalogCrossMode::OppositeSlope, 1.0, 100.0
        );
        CHECK(r.status == AnalogCrossStatus::Disagree);
        CHECK_FALSE(r.has_value);
    }

    TEST_CASE("cross_check_analog falls back to the healthy channel") {
        const auto a_only = cross_check_analog(
            55.0, SignalStatus::Valid, 0.0, SignalStatus::FaultLow, AnalogCrossMode::SameSlope, 2.0
        );
        CHECK(a_only.status == AnalogCrossStatus::ChannelAOnly);
        CHECK(a_only.value == doctest::Approx(55.0));

        const auto b_only = cross_check_analog(
            0.0, SignalStatus::FaultHigh, 40.0, SignalStatus::Valid, AnalogCrossMode::SameSlope, 2.0
        );
        CHECK(b_only.status == AnalogCrossStatus::ChannelBOnly);
        CHECK(b_only.value == doctest::Approx(40.0));

        const auto none = cross_check_analog(
            0.0, SignalStatus::FaultLow, 0.0, SignalStatus::FaultHigh, AnalogCrossMode::SameSlope, 2.0
        );
        CHECK(none.status == AnalogCrossStatus::BothUnusable);
        CHECK_FALSE(none.has_value);
    }

    TEST_CASE("RedundantAnalogPair wires two AnalogInputs through cross-check") {
        const auto                      mon = RangeMonitor<double>{0.25, 0.5, 4.5, 4.75};
        const auto                      cal = two_point_cal(0.5, 0.0, 4.5, 100.0);
        io::RedundantAnalogPair<double> pair{
            io::AnalogInput<double>{cal, mon},
            io::AnalogInput<double>{cal, mon},
            AnalogCrossMode::SameSlope,
            2.0
        };

        // Both mid-scale → 50%.
        auto r = pair.update(2.5, 2.5);
        CHECK(r.agree());
        CHECK(r.value == doctest::Approx(50.0));

        // Channel B open-wire → A-only degraded.
        r = pair.update(2.5, 0.1);
        CHECK(r.status == AnalogCrossStatus::ChannelAOnly);
        CHECK(r.value == doctest::Approx(50.0));
    }

    TEST_CASE("ProportionalCurrentMap maps [0,1] command to mA window") {
        io::ProportionalCurrentMap<double> valve{200.0, 800.0};
        CHECK(valve(0.0) == doctest::Approx(200.0));
        CHECK(valve(1.0) == doctest::Approx(800.0));
        CHECK(valve(0.5) == doctest::Approx(500.0));
        CHECK(valve(-1.0) == doctest::Approx(200.0)); // clamp low
        CHECK(valve(2.0) == doctest::Approx(800.0));  // clamp high
        CHECK(valve.invert(500.0) == doctest::Approx(0.5));
    }

    TEST_CASE("AxisInput chains cal, dead zone, expo, and scale") {
        // Identity cal (already normalized), 10% dead zone, k=1 expo, unit scale.
        io::AxisInput<double> axis{AffineCal<double>{1.0, 0.0}, 0.1, 1.0, 1.0};

        CHECK(axis.update(0.05) == doctest::Approx(0.0));  // inside dead zone
        CHECK(axis.update(1.0) == doctest::Approx(1.0));   // full deflection -> full output
        CHECK(axis.update(-1.0) == doctest::Approx(-1.0)); // symmetric

        // No dead zone / no expo / scaled & inverted: pure -2x.
        io::AxisInput<double> inv{AffineCal<double>{1.0, 0.0}, 0.0, 0.0, 2.0, true};
        CHECK(inv.update(0.5) == doctest::Approx(-1.0));
    }

    TEST_CASE("Button reports edges, level, and hold") {
        io::Button<double> btn{0.02, 0.5}; // 20 ms debounce, 500 ms hold
        const double       dt = 0.01;

        btn.update(true, dt); // 10 ms, not debounced yet
        CHECK(btn.down() == false);
        btn.update(true, dt); // 20 ms -> commits down
        CHECK(btn.down() == true);
        CHECK(btn.pressed() == true); // rising edge this tick
        btn.update(true, dt);
        CHECK(btn.pressed() == false); // edge is one-shot
        CHECK(btn.held() == false);    // not held long enough

        // Hold it down past 500 ms total.
        for (int i = 0; i < 60; ++i) {
            btn.update(true, dt);
        }
        CHECK(btn.held() == true);

        // Release: falling edge, then down() clears after debounce.
        btn.update(false, dt);
        btn.update(false, dt);
        CHECK(btn.down() == false);
        CHECK(btn.released() == true);
    }

    TEST_CASE("differential_drive mixes and preserves ratio on saturation") {
        // Straight line: equal wheels, no turn component.
        auto straight = io::differential_drive(0.5, 0.0);
        CHECK(straight.left == doctest::Approx(0.5));
        CHECK(straight.right == doctest::Approx(0.5));

        // Zero input -> zero output.
        auto zero = io::differential_drive(0.0, 0.0);
        CHECK(zero.left == doctest::Approx(0.0));
        CHECK(zero.right == doctest::Approx(0.0));

        // In-range combined command passes through unscaled.
        auto ok = io::differential_drive(0.6, 0.3);
        CHECK(ok.left == doctest::Approx(0.9));
        CHECK(ok.right == doctest::Approx(0.3));

        // Saturating command: peak = |0.8+0.5| = 1.3, both scaled by 1/1.3.
        // The turn/throttle *ratio* (left/right) must survive normalization.
        auto sat = io::differential_drive(0.8, 0.5);
        CHECK(sat.left == doctest::Approx(1.3 / 1.3)); // 1.0
        CHECK(sat.right == doctest::Approx(0.3 / 1.3));
        CHECK(sat.left / sat.right == doctest::Approx(1.3 / 0.3)); // ratio preserved
    }

    TEST_CASE("arcade_drive is differential_drive with single-stick args") {
        auto a = io::arcade_drive(0.5, 0.8); // x=steer, y=throttle
        auto d = io::differential_drive(0.8, 0.5);
        CHECK(a.left == doctest::Approx(d.left));
        CHECK(a.right == doctest::Approx(d.right));
    }

    TEST_CASE("h_pattern passes each lever to its own track, clamped") {
        auto straight = io::h_pattern(0.5, 0.5);
        CHECK(straight.left == doctest::Approx(0.5));
        CHECK(straight.right == doctest::Approx(0.5));

        auto spin = io::h_pattern(1.0, -1.0); // opposed levers -> spin in place
        CHECK(spin.left == doctest::Approx(1.0));
        CHECK(spin.right == doctest::Approx(-1.0));

        auto clamped = io::h_pattern(1.5, -2.0); // input overshoot guarded
        CHECK(clamped.left == doctest::Approx(1.0));
        CHECK(clamped.right == doctest::Approx(-1.0));
    }

    TEST_CASE("iso_s_drive is the speed-summed (differential) travel stick") {
        auto s = io::iso_s_drive(0.8, 0.5);
        auto d = io::differential_drive(0.8, 0.5);
        CHECK(s.left == doctest::Approx(d.left));
        CHECK(s.right == doctest::Approx(d.right));

        // Rotation direction is independent of travel direction: right steer
        // gives CW (left track leads) both forward and in reverse.
        auto fwd_r = io::iso_s_drive(0.5, 0.5);
        auto rev_r = io::iso_s_drive(-0.5, 0.5);
        CHECK(fwd_r.left > fwd_r.right); // CW forward
        CHECK(rev_r.left > rev_r.right); // still CW in reverse
    }

    TEST_CASE("iso_c_drive scales turn with throttle (no zero-speed spin)") {
        // Straight: equal wheels.
        auto straight = io::iso_c_drive(0.8, 0.0);
        CHECK(straight.left == doctest::Approx(0.8));
        CHECK(straight.right == doctest::Approx(0.8));

        // Full turn at zero throttle -> no motion (the ISO-C defining property).
        auto pivot0 = io::iso_c_drive(0.0, 1.0);
        CHECK(pivot0.left == doctest::Approx(0.0));
        CHECK(pivot0.right == doctest::Approx(0.0));

        // Partial: left = 0.5*(1+0.5)=0.75, right = 0.5*(1-0.5)=0.25, no clip.
        auto part = io::iso_c_drive(0.5, 0.5);
        CHECK(part.left == doctest::Approx(0.75));
        CHECK(part.right == doctest::Approx(0.25));

        // Full throttle + full turn: pre-clip {2,0}, peak 2 -> pivot about inner track.
        auto pivot = io::iso_c_drive(1.0, 1.0);
        CHECK(pivot.left == doctest::Approx(1.0));
        CHECK(pivot.right == doctest::Approx(0.0));

        // Car-like: rotation direction flips with travel direction. Right steer
        // is CW forward but CCW in reverse.
        auto fwd_r = io::iso_c_drive(0.5, 0.5);
        auto rev_r = io::iso_c_drive(-0.5, 0.5);
        CHECK(fwd_r.left > fwd_r.right); // CW forward
        CHECK(rev_r.left < rev_r.right); // CCW in reverse
    }

    TEST_CASE("steer_map reads inner-track geometry from a 2D curve") {
        // 3x3 inner-track LUT over (|throttle|, |turn|): the calibration where a
        // full-turn pivot happens at half throttle, yet the inner track holds
        // 0.5 at full throttle. z(row = |throttle|, col = |turn|).
        constexpr Lut2D<3, 3, double> inner{
            {0.0, 0.5, 1.0},   // |throttle| breakpoints
            {0.0, 0.5, 1.0},   // |turn| breakpoints
            {{0.0, 0.00, 0.0}, // throttle 0.0
             {0.5, 0.25, 0.0}, // throttle 0.5 -> pivot at full turn
             {1.0, 0.75, 0.5}} // throttle 1.0 -> inner holds 0.5 at full turn
        };

        // Straight: both tracks = throttle.
        auto straight = io::steer_map(inner, 0.8, 0.0);
        CHECK(straight.left == doctest::Approx(0.8));
        CHECK(straight.right == doctest::Approx(0.8));

        // Company calibration points; right turn -> right is the inner track.
        auto pivot = io::steer_map(inner, 0.5, 1.0);
        CHECK(pivot.left == doctest::Approx(0.5));  // outer runs at throttle
        CHECK(pivot.right == doctest::Approx(0.0)); // inner pivots
        auto fast = io::steer_map(inner, 1.0, 1.0);
        CHECK(fast.left == doctest::Approx(1.0));  // outer
        CHECK(fast.right == doctest::Approx(0.5)); // inner holds 0.5

        // Left turn mirrors: left becomes the inner track.
        auto left = io::steer_map(inner, 1.0, -1.0);
        CHECK(left.left == doctest::Approx(0.5));
        CHECK(left.right == doctest::Approx(1.0));

        // Reverse-right is car-like: inner still the right track, travel reversed.
        auto rev = io::steer_map(inner, -0.5, 1.0);
        CHECK(rev.left == doctest::Approx(-0.5));
        CHECK(rev.right == doctest::Approx(0.0));
    }

    TEST_CASE("steer_map interpolation is swappable via the surface callable") {
        // Any (|throttle|,|turn|) -> inner callable works; here a closed-form
        // expo-shaped lambda instead of a LUT.
        auto expo_curve = [](double thr, double turn) {
            return thr * (1.0 - (turn * turn * turn)); // soft near-center turn falloff
        };
        auto mid = io::steer_map(expo_curve, 1.0, 0.5);
        CHECK(mid.left == doctest::Approx(1.0));          // outer
        CHECK(mid.right == doctest::Approx(1.0 - 0.125)); // inner = 1*(1-0.5^3)

        // A SplineSurface is a drop-in surface too (smooth interpolation).
        constexpr SplineSurface<3, 3, double> inner{
            {0.0, 0.5, 1.0},
            {0.0, 0.5, 1.0},
            {{0.0, 0.0, 0.0}, {0.5, 0.25, 0.0}, {1.0, 0.75, 0.5}}
        };
        auto sp = io::steer_map(inner, 1.0, 1.0); // grid node -> exact
        CHECK(sp.left == doctest::Approx(1.0));
        CHECK(sp.right == doctest::Approx(0.5));
    }

    TEST_CASE("mecanum_drive recovers cardinal patterns and normalizes") {
        // Pure forward: all four wheels equal and forward.
        auto fwd = io::mecanum_drive(0.0, 1.0, 0.0);
        CHECK(fwd.front_left == doctest::Approx(1.0));
        CHECK(fwd.front_right == doctest::Approx(1.0));
        CHECK(fwd.rear_left == doctest::Approx(1.0));
        CHECK(fwd.rear_right == doctest::Approx(1.0));

        // Pure strafe right: FL/RR forward, FR/RL reverse (mecanum X pattern).
        auto strafe = io::mecanum_drive(1.0, 0.0, 0.0);
        CHECK(strafe.front_left == doctest::Approx(1.0));
        CHECK(strafe.front_right == doctest::Approx(-1.0));
        CHECK(strafe.rear_left == doctest::Approx(-1.0));
        CHECK(strafe.rear_right == doctest::Approx(1.0));

        // Pure yaw: left side one way, right side the other.
        auto yaw = io::mecanum_drive(0.0, 0.0, 1.0);
        CHECK(yaw.front_left == doctest::Approx(1.0));
        CHECK(yaw.front_right == doctest::Approx(-1.0));
        CHECK(yaw.rear_left == doctest::Approx(1.0));
        CHECK(yaw.rear_right == doctest::Approx(-1.0));

        // Zero input -> zero output.
        auto zero = io::mecanum_drive(0.0, 0.0, 0.0);
        CHECK(zero.front_left == doctest::Approx(0.0));

        // Saturating blend: peak = |1 + 0.5 + 0.5| = 2, all scaled by 1/2.
        auto sat = io::mecanum_drive(0.5, 1.0, 0.5);
        CHECK(sat.front_left == doctest::Approx(2.0 / 2.0)); // 1.0
        CHECK(sat.front_right == doctest::Approx(0.0 / 2.0));
        CHECK(sat.rear_left == doctest::Approx(1.0 / 2.0));
        CHECK(sat.rear_right == doctest::Approx(1.0 / 2.0));

        // mecanum_drive / omniwheel_drive are aliases for holonomic_drive.
        auto h = io::holonomic_drive(0.5, 1.0, 0.5);
        auto o = io::omniwheel_drive(0.5, 1.0, 0.5);
        CHECK(sat.front_left == doctest::Approx(h.front_left));
        CHECK(o.front_left == doctest::Approx(h.front_left));
        CHECK(o.rear_right == doctest::Approx(h.rear_right));
    }

    TEST_CASE("Switch debounces and flags changes") {
        io::Switch<double> sw{0.02}; // 20 ms debounce
        const double       dt = 0.01;

        CHECK(sw.update(false, dt) == false);
        sw.update(true, dt);                // 10 ms, not yet
        CHECK(sw.update(true, dt) == true); // 20 ms -> on
        CHECK(sw.changed() == true);        // flipped this tick
        CHECK(sw.update(true, dt) == true);
        CHECK(sw.changed() == false); // steady

        // A glitch shorter than the debounce window is ignored.
        sw.update(false, dt);
        CHECK(sw.on() == true); // single false sample didn't flip it
        CHECK(sw.changed() == false);
    }
}
