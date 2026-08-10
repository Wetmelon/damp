// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file eskf_sketch.cpp
 * @brief MARG attitude ESKF — flashable sketch (Design Is Deploy)
 *
 * Host: examples toolchain. Target: copy this file + eskf_estimator.hpp;
 * replace mock sensors. Finite-tick SIL: eskf_sil.cpp (same estimator header).
 */

#include "damp/estimation/sensor_fusion.hpp"
#include "damp/math/geometry.hpp"
#include "damp/matrix/matrix.hpp"
#include "eskf_estimator.hpp"

using namespace damp;
using namespace damp::examples_eskf;

#if defined(ARDUINO)
#include <Arduino.h>
#else
static void delay(int) {}
#define setup setup_
#define loop loop_
#endif

// Deploy runtime: float on target
static constinit ESKFOrientationFilter<float, 6> filt{kDesign.as<float>()};

// Mock sensors at rest, identity attitude, ENU — replace with drivers
static Vec3<float> read_accelerometer() {
    return {0.0f, 0.0f, 9.81f}; // specific force ≈ +g on body z
}
static Vec3<float> read_gyroscope() {
    return {0.0f, 0.0f, 0.0f};
}
static Vec3<float> read_magnetometer() {
    return {0.0f, 1.0f, 0.0f};
}

static Euler<float> euler_angles;

void setup() {}

void loop() {
    estimate_period(filt, read_accelerometer(), read_gyroscope(), read_magnetometer());
    euler_angles = filt.orientation().to_euler<EulerOrder::ZYX>();
    (void)euler_angles;
    delay(static_cast<int>(dt * 1000.0f));
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
