// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file lqi_sketch.cpp
 * @brief Position LQI with a clamp — flashable sketch
 *
 * Target: copy this file + lqi_controller.hpp; replace read_state / write_force.
 * SIL: lqi_sil.cpp (same tick).
 */

#include "damp/controllers/lqi.hpp"
#include "damp/matrix/colvec.hpp"
#include "lqi_controller.hpp"

using namespace damp;
using namespace damp::examples_lqi;

#ifdef ARDUINO
#include <Arduino.h>
#else

// Stand-ins for Arduino functions
static void pinMode(int, int) {}
static void delay(int) {}

static constexpr int OUTPUT = 1;

#define setup setup_
#define loop loop_
#endif

// Create an instance of an LQI runtime from the LQI "Law"
static constinit Law<lqi_result, float> controller; // the design's LQI in float, its gains constants of the type

static ColVec<2, float> read_state() {
    return ColVec<2, float>{0.0f, 0.0f};
}

static void write_force(float /*u*/) {}

void setup() {
    pinMode(5, OUTPUT);
}

void loop() {
    const auto  x = read_state();
    const float u = control_period(controller, x, 1.0f);
    write_force(u);
    delay(static_cast<int>(Ts * 1000.0));
}

#ifndef ARDUINO
int main() {
    setup();
    for (int i = 0; i < 10; ++i) {
        loop();
    }
    return 0;
}
#endif
