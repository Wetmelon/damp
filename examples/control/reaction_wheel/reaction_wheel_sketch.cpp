// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file reaction_wheel_sketch.cpp
 * @brief Reaction-wheel cascade PI — flashable sketch (Design Is Deploy)
 *
 * Host: examples toolchain. Target: copy + reaction_wheel_controller.hpp;
 * replace mock sensors / wheel driver.
 *
 * SIL (gains + all laws one-tick): reaction_wheel_sil.cpp
 */

#include "damp/controllers/pid.hpp"
#include "reaction_wheel_controller.hpp"

using namespace damp;
using namespace damp::examples_reaction_wheel;

#if defined(ARDUINO)
#include <Arduino.h>
#else
static void pinMode(int, int) {}
static void delay(int) {}

static constexpr int OUTPUT = 1;
#define setup setup_
#define loop loop_
#endif

static constinit PIController<float> rate_pi = rate_pi_design.discretize(Ts).as<float>();
static constinit PIController<float> angle_pi = angle_pi_design.discretize(Ts).as<float>();

static float theta_ref = 0.0f;

static float read_theta() {
    return 0.0f;
}

static float read_omega() {
    return 0.0f;
}

static void write_wheel_torque(float /*u_nm*/) {}

void setup() {
    pinMode(5, OUTPUT);
    rate_pi.reset();
    angle_pi.reset();
}

void loop() {
    const float theta = read_theta();
    const float omega = read_omega();
    const float u = cascade_period(angle_pi, rate_pi, theta_ref, theta, omega);
    write_wheel_torque(u);
    delay(static_cast<int>(Ts * 1000.0));
}

#if !defined(ARDUINO)
int main() {
    setup();
    for (int i = 0; i < 10; ++i) {
        loop();
    }
    return 0;
}
#endif
