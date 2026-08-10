// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file ins_eskf_sketch.cpp
 * @brief InsNavigator thin smoke — same design as ins_eskf_sil (float deploy)
 *
 * Full 3D SIL: ins_eskf_sil.cpp. Firmware product sketch also lives at
 * estimation/ins_navigator/ (same densities).
 */

#include "damp/estimation/ins_eskf.hpp"
#include "damp/estimation/ins_mechanization.hpp"
#include "damp/math/geometry.hpp"
#include "damp/matrix/matrix.hpp"
#include "ins_eskf_estimator.hpp"

using namespace damp;
using namespace damp::examples_ins_eskf;

#if defined(ARDUINO)
#include <Arduino.h>
#else
static void delay(int) {}
#define setup setup_
#define loop loop_
#endif

static constinit InsNavigator<float> nav{kInsDesign.as<float>(), kFrame};

static Vec3<float> read_gyro() {
    return {0.0f, 0.0f, 0.0f};
}
static Vec3<float> read_accel() {
    return specific_force_at_rest(nav.state().q, kFrame);
}

void setup() {}

void loop() {
    const ImuSample<float> imu{.gyro = read_gyro(), .accel = read_accel()};
    nav.predict(imu, static_cast<float>(kTs));
    // Sparse absolute aids on target (GPS / dual-antenna) — mock hold
    (void)nav.update_position(nav.state().p);
    const Matrix<1, 1, float> R_psi{{static_cast<float>(kHeadingStd * kHeadingStd)}};
    (void)nav.update_heading(nav.heading(), R_psi);
    (void)nav.state();
    delay(static_cast<int>(kTs * 1000.0));
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
