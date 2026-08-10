// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file encoder_velocity_sketch.cpp
 * @brief Encoder PLL velocity — thin deploy-shaped smoke
 *
 * Host comparison plots: encoder_velocity_sil.cpp.
 */

#include "damp/filters/differentiator.hpp"
#include "encoder_velocity_estimator.hpp"

using namespace damp;
using namespace damp::examples_encoder_velocity;

#if defined(ARDUINO)
#include <Arduino.h>
#else
static void delayMicroseconds(int) {}
#define setup setup_
#define loop loop_
#endif

// 14-bit absolute encoder scenario defaults
static constexpr double kCpr = 16384.0;
static constexpr double kBw = 1000.0;

static PllObserver                       pll{kBw};
static RobustExactDifferentiator<double> levant{2.0 * 0.25 * (2.0 * 3.141592653589793 * 1.0) * (2.0 * 3.141592653589793 * 1.0), dt};

// Mock quantized position [turns] — replace with encoder counts / cpr
static double read_encoder_turns() {
    return 0.0;
}

void setup() {}

void loop() {
    const double th = read_encoder_turns();
    const double v_pll = estimate_period_pll(pll, th);
    const double v_lev = estimate_period_levant(levant, th);
    (void)v_pll;
    (void)v_lev;
    delayMicroseconds(static_cast<int>(dt * 1.0e6));
}

#if !defined(ARDUINO)
int main() {
    setup();
    for (int i = 0; i < 100; ++i) {
        loop();
    }
    return 0;
}
#endif
