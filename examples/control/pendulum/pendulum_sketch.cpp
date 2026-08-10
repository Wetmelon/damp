// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file pendulum_sketch.cpp
 * @brief Upright pendulum LQR — flashable sketch (Design Is Deploy)
 *
 * Host: examples toolchain. Target: copy this file + pendulum_controller.hpp;
 * replace mock sensors/actuators.
 *
 * Closed-loop SIL: pendulum_sil.cpp (same controller header).
 */

#include "damp/controllers/lqr.hpp"
#include "damp/matrix/colvec.hpp"
#include "pendulum_controller.hpp"

using namespace damp;
using namespace damp::examples_pendulum;

#if defined(ARDUINO)
#include <Arduino.h>
#else
static void pinMode(int, int) {}
static void delay(int) {}

static constexpr int OUTPUT = 1;
#define setup setup_
#define loop loop_
#endif

// Deploy runtime: float on target
static constinit LQR<2, 1, float> controller{lqr_d.as<float>()};

// Mock full state (replace with sensors on hardware)
static ColVec<2, float> read_state() {
    return ColVec<2, float>{0.0f, 0.0f};
}

static void write_torque(float /*u_nm*/) {
    // Map torque to motor/PWM on target
}

void setup() {
    pinMode(5, OUTPUT);
}

void loop() {
    const ColVec<2, float> x = read_state();
    const float            u = control_period(controller, x)(0);
    write_torque(u);
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
