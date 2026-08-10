// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file fo_plant_inertia_sketch.cpp
 * @brief FO plant inertia ID — thin deploy-shaped smoke
 *
 * Host SIL with true plant + print: fo_plant_inertia_sil.cpp.
 */

#include "damp/estimation/parameter_estimation.hpp"
#include "fo_plant_inertia_estimator.hpp"

using namespace damp;
using namespace damp::examples_fo_plant_inertia;

#if defined(ARDUINO)
#include <Arduino.h>
#else
static void delay(int) {}
#define setup setup_
#define loop loop_
#endif

static constexpr design::FirstOrderPlantEstimatorConfig<float> kCfgF{
    .Ts = static_cast<float>(Ts),
    .forgetting = 0.999f,
    .initial_covariance = 1e4f,
};
static FirstOrderPlantEstimator est{kCfgF};
static float                    omega = 0.0f;
static std::size_t              k = 0;

// Mock measured omega / commanded torque on target
static float read_omega() {
    return omega;
}
static float command_torque() {
    return static_cast<float>(0.5 * excitation(k, 10));
}

void setup() {}

void loop() {
    const float torque = command_torque();
    // On target the plant is the motor; here we only update the estimator
    (void)est.update(torque, read_omega());
    ++k;
    delay(static_cast<int>(Ts * 1000.0));
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
