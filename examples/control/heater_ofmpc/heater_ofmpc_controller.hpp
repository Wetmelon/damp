// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file heater_ofmpc_controller.hpp
 * @brief Heater with dead time — offset-free MPC, the delay predicted (flashable)
 *
 * Folder layout:
 *
 *   heater_ofmpc_controller.hpp  — this file (nameplate, plant, design, the flashable tick)
 *   heater_ofmpc_sketch.cpp      — Arduino setup/loop
 *   heater_ofmpc_sil.cpp         — host plant: the MPC vs the same MPC blind to the delay,
 *                                  and the setpoint known ahead vs held
 *   heater_ofmpc_derivation.md   — FOPDT, the predicted delay, offset-free MPC, preview
 *
 * A 3D-printer hotend: a 40 W cartridge in an aluminium block, read by a thermistor. Heat
 * takes 5 s to reach the thermistor, two and a half control periods. The MPC predicts where
 * the block will be when the power it decides now arrives, so it stops heating before the
 * thermistor has seen the heat it already put in.
 */

#pragma once

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/controllers/offset_free_mpc.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/systems/state_space.hpp"

namespace damp::examples_heater_ofmpc {

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

inline constexpr std::size_t kHorizon = 30; // [periods] 60 s, two time constants

/// 30 periods ahead, 10 moves, the dead time in periods. Move the power gently (Rdu = 10, see
/// the derivation), 0 … Pmax; model error on the block, thermistor noise, and how fast the
/// unknown heat loss may drift.
template<double DelayPeriods>
[[nodiscard]] constexpr auto design_mpc() {
    return design::mpc<kHorizon, 10, DelayPeriods>(
        kPlant, kTs,
        {.Rdu = {{10.0}}},
        {.u_min = {0.0}, .u_max = {kPmax}},
        {.process = {{0.01}}, .measurement = {{0.01}}, .disturbance = {{0.05}}}
    );
}

inline constexpr auto kLaw = design_mpc<kDeadTime / kTs>(); // designed at compile time; a failed design does not compile into a constinit controller

/// One control period: thermistor reading in, heater power out [W].
template<typename Mpc>
[[nodiscard]] inline float heater_period(Mpc& mpc, float setpoint_c, float temperature_c) {
    const auto ambient = static_cast<float>(kAmbient);
    return mpc.control(setpoint_c - ambient, temperature_c - ambient);
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
        r[k] = ColVec<1, float>{setpoint_c[k] - ambient};
    }
    return mpc.control(r, ColVec<1, float>{temperature_c - ambient})[0];
}

} // namespace damp::examples_heater_ofmpc
