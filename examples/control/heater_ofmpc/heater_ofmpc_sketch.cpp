// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file heater_ofmpc_sketch.cpp
 * @brief Heater with dead time — offset-free MPC, flashable sketch (Design Is Deploy)
 *
 * Target: copy this file + heater_ofmpc_controller.hpp; replace the thermistor read and the
 * heater PWM write. Closed-loop SIL: heater_ofmpc_sil.cpp (same tick).
 */

#include "damp/controllers/offset_free_mpc.hpp"
#include "heater_ofmpc_controller.hpp"

using namespace damp;
using namespace damp::examples_heater_ofmpc;

#if defined(ARDUINO)
#include <Arduino.h>
#else
static void pinMode(int, int) {}
static void analogWrite(int, int) {}
static void delay(int) {}

static constexpr int OUTPUT = 1;
#define setup setup_
#define loop loop_
#endif

static constexpr int kHeaterPin = 5;

// Designed at compile time; float on target.
static constinit OffsetFreeMPC controller{kLaw.as<float>()};
static constexpr float         kSetpointC = 210.0f;

// Mock thermistor (replace with the ADC read and the thermistor curve).
static float read_temperature_c() {
    return 25.0f;
}

void setup() {
    pinMode(kHeaterPin, OUTPUT);
}

void loop() {
    const float power_w = heater_period(controller, kSetpointC, read_temperature_c());
    analogWrite(kHeaterPin, static_cast<int>((power_w / static_cast<float>(kPmax)) * 255.0f)); // duty = P / Pmax
    delay(static_cast<int>(kTs * 1000.0));
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
