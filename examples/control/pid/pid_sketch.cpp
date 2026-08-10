// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file pid_sketch.cpp
 * @brief PI — flashable sketch (Design Is Deploy)
 *
 * Host: examples toolchain (mocks for IO).
 * Target: copy this file + pid_controller.hpp; replace mock ADC/PWM.
 *
 * SIL: pid_sil.cpp (same controller header).
 */

#include "damp/controllers/pid.hpp"
#include "pid_controller.hpp"

using namespace damp;
using namespace damp::examples_pid;

#if defined(ARDUINO)
#include <Arduino.h>
#else
static void pinMode(int, int) {}
static void delay(int) {}

static constexpr int OUTPUT = 1;
#define setup setup_
#define loop loop_
#endif

static constinit PIController<float> controller = pi_design.discretize(static_cast<double>(Ts)).as<float>();

static float reference = 1.0f;

static float read_process_variable() {
    return 0.0f; // mock: plant at rest
}

static void write_actuator(float /*u*/) {}

void setup() {
    pinMode(5, OUTPUT);
    controller.reset();
}

void loop() {
    const float y = read_process_variable();
    const float u = control_period(controller, reference, y);
    write_actuator(u);
    delay(static_cast<int>(Ts * 1000.0f));
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
