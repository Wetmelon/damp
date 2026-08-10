// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file workflow_end_to_end_sketch.cpp
 * @brief Workflow e2e — flashable float runtime smoke (after host design)
 *
 * Full linearize + LQGI synthesis is host-side (workflow_end_to_end_sil.cpp).
 * Sketch shows the per-tick shape once gains are baked.
 */

#include "damp/controllers/pr.hpp"
#include "workflow_end_to_end_controller.hpp"

using namespace damp;
using namespace damp::examples_workflow_e2e;

#if defined(ARDUINO)
#include <Arduino.h>
#else
static void pinMode(int, int) {}
static void delay(int) {}

static constexpr int OUTPUT = 1;
#define setup setup_
#define loop loop_
#endif

// On target: constinit from design artifacts.as<float>() after offline/host design.
// Smoke: PR only with a fixed design (LQGI runtime needs host synthesis).
static constexpr auto kPr = design::pr(0.0, 10.0, 2.0 * damp::numbers::pi_v<double>, 6.0, kTs);
static_assert(kPr.success);
static constinit PRController<float> pr{kPr.as<float>()};
static float                         y = 0.25f;
static float                         reference = 0.5f;

void setup() {
    pinMode(5, OUTPUT);
}

void loop() {
    const float u = pr.control(reference - y);
    y += 0.01f * u; // mock
    delay(static_cast<int>(kTs * 1000.0));
}

#if !defined(ARDUINO)
int main() {
    setup();
    for (int i = 0; i < 20; ++i) {
        loop();
    }
    return 0;
}
#endif
