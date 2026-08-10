// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file ins_navigator_sketch.cpp
 * @brief 15-state INS navigator — flashable sketch (Design Is Deploy)
 *
 * Firmware-shaped. Copy onto an MCU with an IMU + absolute aids (GPS position
 * and optional dual-antenna heading). Host build runs a finite number of ticks.
 *
 * Host 3D SIL of the same densities: estimation/ins_eskf/
 */

#include "damp/estimation/ins_eskf.hpp"
#include "damp/estimation/ins_mechanization.hpp"
#include "damp/math/geometry.hpp"
#include "damp/matrix/matrix.hpp"
#include "ins_navigator_estimator.hpp"

using namespace damp;
using namespace damp::examples_ins_navigator;

#if defined(ARDUINO)
#include <Arduino.h>
#else
static void delay(int) {}
#define setup setup_
#define loop loop_
#endif

static constinit InsNavigator<float> nav{kInsDesign.as<float>(), kFrame};

// Mock sensors — replace with driver reads on target
static Vec3<float> read_gyro() {
    return {0.0f, 0.0f, 0.0f};
}
static Vec3<float> read_accel() {
    // ENU at rest, q = I → specific force ≈ (0,0,-g) in nav; body via state
    return specific_force_at_rest(nav.state().q, kFrame);
}
static bool read_gps_position(Vec3<float>& p_out) {
    p_out = nav.state().p; // mock: “perfect” hold — replace with GPS ECEF→ENU
    return true;
}
static bool read_dual_antenna_heading(float& psi_out) {
    psi_out = nav.heading(kBaselineBody);
    return true;
}

void setup() {
    // Optional: seed pose from first GPS fix
}

void loop() {
    const ImuSample<float> imu{.gyro = read_gyro(), .accel = read_accel()};
    Vec3<float>            p_gps{};
    const bool             have_gps = read_gps_position(p_gps);
    float                  psi = 0.0f;
    const bool             have_psi = read_dual_antenna_heading(psi);

    estimate_period(nav, imu, have_gps, p_gps, have_psi, psi, kBaselineBody);

    // Use nav.state().p / .v / .q / .b_g / .b_a for guidance, logging, etc.
    (void)nav.state();
    delay(static_cast<int>(kTs * 1000.0f));
}

#if !defined(ARDUINO)
int main() {
    setup();
    for (int i = 0; i < 200; ++i) {
        loop();
    }
    return 0;
}
#endif
