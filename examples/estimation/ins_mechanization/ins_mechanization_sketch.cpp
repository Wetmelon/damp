// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file ins_mechanization_sketch.cpp
 * @brief Strapdown mechanization — thin rest-frame smoke
 *
 * Host animation: ins_mechanization_sil.cpp.
 */

#include "damp/estimation/ins_mechanization.hpp"
#include "damp/math/geometry.hpp"
#include "ins_mechanization_estimator.hpp"

using namespace damp;
using namespace damp::examples_ins_mechanization;

#if defined(ARDUINO)
#include <Arduino.h>
#else
static void delay(int) {}
#define setup setup_
#define loop loop_
#endif

static InsState<float> x{};

void setup() {
    x.q = Quaternion<float>::identity();
}

void loop() {
    const ImuSample<float> imu{
        .gyro = {0.0f, 0.0f, 0.0f},
        .accel = specific_force_at_rest(x.q, kFrame),
    };
    x = estimate_period(x, imu, static_cast<float>(kTsSketch), kFrame);
    (void)x;
    delay(static_cast<int>(kTsSketch * 1000.0));
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
