// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file heater_lqg_controller.hpp
 * @brief Heater with dead time — offset-free LQG, the dead time in one sample: the fewest lines that design and step it
 *
 * The same hotend as heater_ofmpc/ (40 W cartridge, thermistor 5 s away, 5 s period) under the
 * unconstrained law the MPC reduces to: LQR on the block temperature, the disturbance-augmented
 * Kalman filter, a steady-state target from the estimated heat loss. The law knows nothing of
 * the 0 … Pmax limit, so the tick clamps its request at the actuator and tells the estimator
 * what was applied. The dead time is one control period, one register state in the model, as in
 * heater_ofmpc/.
 */

#pragma once

#include "damp/backend.hpp"
#include "damp/controllers/oflqg.hpp"
#include "damp/controllers/state_feedback_core.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/core.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

namespace damp::examples_heater_lqg {

// Nameplate
inline constexpr double kGain = 6.0;     // [°C/W] steady rise per watt
inline constexpr double kTau = 30.0;     // [s] block time constant
inline constexpr double kDeadTime = 5.0; // [s] heater → thermistor
inline constexpr double kPmax = 40.0;    // [W]
inline constexpr double kAmbient = 25.0; // [°C] assumed; the disturbance estimate absorbs the error
inline constexpr double kTs = 5.0;       // [s] control period: the dead time, so it is one sample

// Plant: first order in degrees above ambient. The dead time is one period: the power decided now
// reaches the block over the next period, so the model sampled at Ts carries the last period's
// power as one register state.
inline constexpr StateSpace<1, 1, 1, double, 1, 1> kPlant{.A = {{-1.0 / kTau}}, .B = {{kGain / kTau}}, .C = {{1.0}}}; // the model error on the block, the noise on the reading
inline constexpr auto                              kModel = *discretize(kPlant, kTs, InputDelay{kDeadTime});          // [T; P(k−1)], per period

// Design on the model, per period: weight on the temperature and on the power itself (the
// register unweighted); model error on the block, thermistor noise, and how fast the unknown heat
// loss may drift (0.01 and 0.05 per second over the 5 s period).
inline constexpr auto kLaw = design::oflqg(
    kModel,
    {
        .Q = {{1.0, 0.0}, {0.0, 0.0}},
        .R = 0.02,
    },
    {
        .process = 0.05,
        .measurement = 0.01,
        .disturbance = 0.25,
    },
    design::InputDisturbance{}
); // a failed design does not compile into a constinit controller

// Runtime: the design's law, float on the board; its sizes come from the result.
using Heater = Law<kLaw, float>;

/// One control period: thermistor reading in, heater power out [W]. The request is clamped
/// to 0 … Pmax and the estimator is propagated with what was applied (update, control, then
/// back_calculate and predict once the power is out).
[[nodiscard]] inline float heater_period(Heater& law, float setpoint_c, float temperature_c) {
    const auto ambient = static_cast<float>(kAmbient);
    law.update(ColVec<1, float>{temperature_c - ambient});
    const float requested = law.control(ColVec<1, float>{setpoint_c - ambient})(0);
    const float applied = damp::clamp(requested, 0.0f, static_cast<float>(kPmax));
    law.back_calculate(ColVec<1, float>{applied});
    law.predict();
    return applied;
}

} // namespace damp::examples_heater_lqg
