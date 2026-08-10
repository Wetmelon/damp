// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file imu_pose_sketch.cpp
 * @brief Open-loop gyro attitude — thin smoke (teaching)
 *
 * Host animation: imu_pose_sil.cpp.
 */

#include "damp/math/geometry.hpp"
#include "damp/matrix/matrix.hpp"
#include "imu_pose_estimator.hpp"

using namespace damp;
using namespace damp::examples_imu_pose;

#if defined(ARDUINO)
#include <Arduino.h>
#else
static void delay(int) {}
#define setup setup_
#define loop loop_
#endif

static Quaternion<float> q = Quaternion<float>::identity();

void setup() {}

void loop() {
    // Mock: zero rate on target until a real gyro is wired
    q = estimate_period(q, Vec3<float>{0.0f, 0.0f, 0.0f}, static_cast<float>(kTsSketch));
    (void)q;
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
