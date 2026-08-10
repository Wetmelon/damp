// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <array>
#include <numbers>

#include "damp/toolbox/actuator.hpp"
#include "damp/trajectory/trajectory_types.hpp"
#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

// Bridge from motion-profile output (SI joint units) to drive-native servo commands.

TEST_SUITE("Actuator") {
    constexpr double two_pi = 2.0 * std::numbers::pi;

    TEST_CASE("rotary_gearbox maps rad to motor turns") {
        // 10:1 reduction: one joint revolution (2π rad) → 10 motor turns.
        const auto axis = rotary_gearbox<double>(10.0);
        CHECK(axis.turns_per_unit == doctest::Approx(10.0 / two_pi));
        CHECK(axis.position_command(two_pi) == doctest::Approx(10.0));
        CHECK(axis.velocity_command(two_pi) == doctest::Approx(10.0)); // rad/s → turns/s
        // Homing offset is added to position only, not velocity.
        const auto homed = rotary_gearbox<double>(10.0, 1.0, 2.5);
        CHECK(homed.position_command(0.0) == doctest::Approx(2.5));
        CHECK(homed.velocity_command(0.0) == doctest::Approx(0.0));
    }

    TEST_CASE("direction flip reflects into position, velocity, and torque sign") {
        const auto fwd = rotary_gearbox<double>(5.0, +1.0);
        const auto rev = rotary_gearbox<double>(5.0, -1.0);
        CHECK(rev.position_command(1.0) == doctest::Approx(-fwd.position_command(1.0)));
        CHECK(rev.velocity_command(1.0) == doctest::Approx(-fwd.velocity_command(1.0)));
        // A positive joint torque must reflect to opposite motor-torque signs.
        CHECK(rev.reflect_torque(3.0) == doctest::Approx(-fwd.reflect_torque(3.0)));
    }

    TEST_CASE("torque reflection is power-consistent (gearbox divides torque by n)") {
        // τ_motor = τ_joint / (2π·r) = τ_joint / n for a rotary n:1 gearbox.
        const auto axis = rotary_gearbox<double>(20.0);
        CHECK(axis.reflect_torque(40.0) == doctest::Approx(40.0 / 20.0));
        // Power balance: τ_motor·ω_motor == τ_joint·q̇ for any rate.
        const double qdot = 3.7;
        const double tau_joint = 8.0;
        const double omega_motor = axis.velocity_command(qdot) * two_pi; // turns/s → rad/s
        CHECK(axis.reflect_torque(tau_joint) * omega_motor == doctest::Approx(tau_joint * qdot));
    }

    TEST_CASE("linear_screw: force reflects to τ = F·L/2π") {
        const double lead = 0.005; // 5 mm/rev
        const auto   axis = linear_screw<double>(lead);
        CHECK(axis.position_command(lead) == doctest::Approx(1.0)); // one lead → one turn
        CHECK(axis.reflect_torque(100.0) == doctest::Approx(100.0 * lead / two_pi));
    }

    TEST_CASE("ConstantInertiaFeedforward: J·a + b·v + τc·sign(v) + g") {
        ConstantInertiaFeedforward<2, double> ff{};
        ff.inertia = {0.5, 2.0};
        ff.viscous = {0.1, 0.0};
        ff.coulomb = {0.3, 0.0};
        ff.gravity = {0.0, 9.81};

        const TrajectoryState<double> s0{0.0, 2.0, 4.0, 0.0}; // v=2, a=4
        CHECK(ff(0, s0) == doctest::Approx((0.5 * 4.0) + (0.1 * 2.0) + (0.3 * 1.0)));
        // Axis 1 is gravity-dominated; Coulomb vanishes at zero velocity.
        const TrajectoryState<double> s1{0.0, 0.0, 1.0, 0.0};
        CHECK(ff(1, s1) == doctest::Approx((2.0 * 1.0) + 9.81));
        // Coulomb term flips with velocity sign.
        const TrajectoryState<double> sneg{0.0, -2.0, 0.0, 0.0};
        CHECK(ff(0, sneg) == doctest::Approx((0.1 * -2.0) + (0.3 * -1.0)));
    }

    TEST_CASE("ServoBank maps a multi-axis state array to commands") {
        ServoBank<2, double> servos{{rotary_gearbox<double>(10.0), rotary_gearbox<double>(4.0)}};

        ConstantInertiaFeedforward<2, double> ff{};
        ff.inertia = {0.5, 0.25};

        damp::array<TrajectoryState<double>, 2> states{
            TrajectoryState<double>{two_pi, two_pi, 1.0, 0.0}, // axis 0
            TrajectoryState<double>{0.0, 0.0, 2.0, 0.0},       // axis 1
        };

        const auto cmds = servos.command(states, ff);
        CHECK(cmds[0].position == doctest::Approx(10.0));
        CHECK(cmds[0].velocity == doctest::Approx(10.0));
        CHECK(cmds[0].torque == doctest::Approx((0.5 * 1.0) / 10.0)); // joint Nm → motor Nm

        // No-policy overload leaves torque at zero.
        const auto cmds_nt = servos.command(states);
        CHECK(cmds_nt[0].torque == doctest::Approx(0.0));
        CHECK(cmds_nt[1].position == doctest::Approx(0.0));
    }

    TEST_CASE("constexpr and precision rebind") {
        constexpr auto   axis = rotary_gearbox<double>(8.0);
        constexpr double p = axis.position_command(1.0);
        static_assert(p > 0.0);

        const auto faxis = axis.as<float>();
        CHECK(static_cast<double>(faxis.turns_per_unit) == doctest::Approx(axis.turns_per_unit));

        ServoBank<1, double> servos{{axis}};
        const auto           fservos = servos.rebind<float>();
        CHECK(static_cast<double>(fservos.axis(0).turns_per_unit) == doctest::Approx(axis.turns_per_unit));
    }
}
