// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file heater_oflqg_controller.hpp
 * @brief Heater with dead time — offset-free LQG, predictor form: the fewest lines that design and step it
 *
 * The same hotend as heater_ofmpc/ (40 W cartridge, thermistor 5 s away, 2 s period) under the
 * unconstrained law the MPC reduces to: LQR on the block temperature, the disturbance-augmented
 * Kalman filter, a steady-state target from the estimated heat loss. The law knows nothing of
 * the 0 … Pmax limit, so the tick clamps its request at the actuator and tells the estimator
 * what was applied. The dead time is handled by prediction, as in heater_ofmpc/.
 */

#pragma once

#include "damp/backend.hpp"
#include "damp/controllers/offset_free_lqg.hpp"
#include "damp/systems/discretization.hpp"

namespace damp::examples_heater_oflqg {

// Nameplate
inline constexpr double kGain = 6.0;     // [°C/W] steady rise per watt
inline constexpr double kTau = 30.0;     // [s] block time constant
inline constexpr double kDeadTime = 5.0; // [s] heater → thermistor
inline constexpr double kPmax = 40.0;    // [W]
inline constexpr double kAmbient = 25.0; // [°C] assumed; the disturbance estimate absorbs the error
inline constexpr double kTs = 2.0;       // [s] control period

// Plant: first order in degrees above ambient. The dead time is 2.5 periods: two predicted at
// run time from the powers already sent, half a period left in the model as one register state.
inline constexpr StateSpace kPlant{.A = {{-1.0 / kTau}}, .B = {{kGain / kTau}}, .C = {{1.0}}};

// Design: weight on the temperature and on the power itself; model error on the block,
// thermistor noise, and how fast the unknown heat loss may drift.
inline constexpr auto kLaw = design::discrete_oflqg<kDeadTime / kTs>(
    kPlant, kTs,
    Matrix<1, 1>{{1.0}},
    Matrix<1, 1>{{0.02}},
    {.process = {{0.01}}, .measurement = {{0.01}}, .disturbance = {{0.05}}}
); // a failed design does not compile into a constinit controller

// Runtime: float on the board; sizes and the delay come from the design result.
using Heater = decltype(OFLQG{kLaw.as<float>()});

/// One control period: thermistor reading in, heater power out [W]. The request is clamped
/// to 0 … Pmax and the estimator is propagated with what was applied (feedback / commit).
[[nodiscard]] inline float heater_period(Heater& law, float setpoint_c, float temperature_c) {
    const auto  ambient = static_cast<float>(kAmbient);
    const float requested = law.feedback(setpoint_c - ambient, temperature_c - ambient);
    const float applied = damp::clamp(requested, 0.0f, static_cast<float>(kPmax));
    law.commit(applied);
    return applied;
}

} // namespace damp::examples_heater_oflqg
