// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file heater_ofmpc_controller.hpp
 * @brief Heater with dead time — offset-free MPC, the dead time in one sample (flashable)
 *
 * Folder layout:
 *
 *   heater_ofmpc_controller.hpp  — this file (nameplate, plant, design, the flashable tick)
 *   heater_ofmpc_sketch.cpp      — Arduino setup/loop
 *   heater_ofmpc_sil.cpp         — host plant: the MPC vs the same MPC blind to the delay,
 *                                  and the setpoint known ahead vs held
 *   heater_ofmpc_derivation.md   — FOPDT, the delay state, offset-free MPC, preview
 *
 * A 3D-printer hotend: a 40 W cartridge in an aluminium block, read by a thermistor. Heat
 * takes 5 s to reach the thermistor, one control period. The model carries the power already sent
 * as a delay state, so the MPC plans from where the block will be when the power it decides
 * now arrives, and stops heating before the thermistor has seen the heat it already put in.
 */

#pragma once

#include <cstddef>

#include "damp/control.hpp"

namespace damp::examples_heater_ofmpc {

// Nameplate
inline constexpr double kGain = 6.0;     // [°C/W] steady rise per watt
inline constexpr double kTau = 30.0;     // [s] block time constant
inline constexpr double kDeadTime = 5.0; // [s] heater → thermistor
inline constexpr double kPmax = 40.0;    // [W]
inline constexpr double kAmbient = 25.0; // [°C] assumed; the disturbance estimate absorbs the error
inline constexpr double kTs = 5.0;       // [s] control period: the dead time, so it is one sample

// Plant: first order in degrees above ambient, the model error on the block temperature (G = 1)
// and the thermistor noise on the reading (H = 1). The dead time is one period, one delay state
// in the model.
inline constexpr StateSpace kPlant{.A = {{-1.0 / kTau}}, .B = {{kGain / kTau}}, .C = {{1.0}}, .G = Matrix{{1.0}}, .H = Matrix{{1.0}}};
inline constexpr auto       kModel = *discretize(kPlant, kTs, InputDelay{kDeadTime}); // [T; P(k−1)], per period

inline constexpr std::size_t kHorizon = 12; // [periods] 60 s, two time constants

/// Per period, on the discrete model: the temperature error at Qy = 1, the power's change at
/// Rdu = 1.6, which moves the power gently (see the derivation).
inline constexpr design::MPCWeights kWeights{.Rdu = {1.6}};
/// The cartridge: 0 … Pmax.
inline constexpr design::InputLimits kLimits{.lower = {0.0}, .upper = {kPmax}};
/// Model error on the block, thermistor noise (one reading's variance), and how fast the unknown
/// heat loss may drift, per period (0.01 and 0.05 per second over 5 s).
inline constexpr design::KalmanNoise kNoise{
    .process = 0.05,
    .measurement = 0.01,
    .disturbance = 0.25,
};

/// 12 periods ahead, 4 moves, offset-free under a heat loss at the heater's input, on the model
/// with the dead time. Designed at compile time; a failed design does not compile into a
/// constinit controller.
inline constexpr auto kLaw = design::mpc<kHorizon, 4>(kModel, kWeights, kLimits, kNoise, design::InputDisturbance{});

/// One control period: thermistor reading in, heater power out [W].
template<typename Mpc>
[[nodiscard]] inline float heater_period(Mpc& mpc, float setpoint_c, float temperature_c) {
    const auto ambient = static_cast<float>(kAmbient);
    return mpc.step(setpoint_c - ambient, temperature_c - ambient)[0];
}

/// The setpoint over the next kHorizon periods (°C), when the print job knows it ahead.
using SetpointAhead = damp::array<float, kHorizon>;

/// The same period with the setpoint known over the horizon: the MPC plans a change in the
/// setpoint before it happens (it stops heating ahead of a step down) instead of after.
template<typename Mpc>
[[nodiscard]] inline float heater_period(Mpc& mpc, const SetpointAhead& setpoint_c, float temperature_c) {
    const auto            ambient = static_cast<float>(kAmbient);
    typename Mpc::Preview r{};
    for (std::size_t k = 0; k < kHorizon; ++k) {
        r[k] = setpoint_c[k] - ambient;
    }
    return mpc.step(r, temperature_c - ambient)[0];
}

} // namespace damp::examples_heater_ofmpc
