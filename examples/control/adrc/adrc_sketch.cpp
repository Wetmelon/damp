// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file adrc_sketch.cpp
 * @brief Position ADRC — flashable sketch (Design Is Deploy)
 *
 * Host: examples toolchain. Target: copy this file + adrc_controller.hpp;
 * replace mock sensors/actuators.
 *
 * Closed-loop SIL: adrc_sil.cpp (same ADRC object; PI-D is SIL-only compare).
 */

#include "adrc_controller.hpp"
#include "damp/controllers/adrc.hpp"

using namespace damp;
using namespace damp::examples_adrc;

#if defined(ARDUINO)
#include <Arduino.h>
#else
static void pinMode(int, int) {}
static void delay(int) {}

static constexpr int OUTPUT = 1;
#define setup setup_
#define loop loop_
#endif

static constinit ADRCController<2, float> controller{adrc_d.as<float>(), static_cast<float>(Ts)};

static float reference = 1.0f;

static float read_angle() {
    return 0.0f; // mock encoder
}

static void write_torque(float /*u_nm*/) {}

void setup() {
    pinMode(5, OUTPUT);
    controller.reset();
}

void loop() {
    const float y = read_angle();
    const float u = control_period_adrc(controller, reference, y);
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
