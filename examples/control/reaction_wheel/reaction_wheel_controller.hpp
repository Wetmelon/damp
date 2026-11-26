// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file reaction_wheel_controller.hpp
 * @brief Single-axis reaction-wheel stabilization — controller synthesis gallery
 *
 * Folder layout:
 *
 *   reaction_wheel_controller.hpp  — this file (all designs + cascade_period)
 *   reaction_wheel_sketch.cpp      — one-tick cascade smoke (float deploy)
 *   reaction_wheel_sil.cpp         — host: print gains + one tick each law
 *   reaction_wheel_derivation.md   — plant + design map
 *
 * Teaching path (one plant, several designs, same deploy shape):
 *
 *   1. Cascade PI (θ → ω* → u)     — industry default for a rate-capable axis
 *   2. ADRC / SMC / STSMC          — light-model SISO alternatives on θ
 *   3. placed / LQR / LQI          — full-state linear design
 *   4. LQG / LQGI                  — state feedback + Kalman from y
 *   5. Offset-free MPC             — constrained predictive (synthesis only here)
 *
 * Do not put plant ODE, std::vector, or plotly here.
 */

#pragma once

#include <cstddef>

#include "damp/control.hpp"

namespace damp::examples_reaction_wheel {

inline constexpr double Ts = 0.001; // [s]  1 kHz control loop
inline constexpr double w0 = 1.0;   // [rad/s] open-loop unstable natural frequency
inline constexpr double b_u = 1.0;  // [rad/s² per N·m]  ω̇ = … + b_u·u

inline constexpr size_t NX = 2;
inline constexpr size_t NU = 1;
inline constexpr size_t NY = 1;

// Continuous plant + identity noise channels so Kalman-based designs have a
// well-posed process/measurement model (G, H = I).
inline constexpr StateSpace plant{
    .A = {{0.0, 1.0}, {w0 * w0, 0.0}},
    .B = {{0.0}, {b_u}},
    .C = {{1.0, 0.0}},
    .G = Matrix<NX, NX>::identity(),
    .H = Matrix<NY, NY>::identity(),
};

inline constexpr auto sysd = *discretize(plant, Ts, DiscretizationMethod::ZOH);

// Shared LQ* weights
inline constexpr design::LQRWeights lqr_weights = {
    .Q = Matrix<NX, NX>::diagonal({40.0, 4.0}),
    .R = {1.0},
};

inline constexpr design::IntegralAction integral = {.Q = 80.0};

inline constexpr design::KalmanNoise estimator_noise = {
    .process = Matrix<NX, NX>::diagonal({1e-4, 1e-3}),
    .measurement = 1e-3,
};

// Cascade PI
inline constexpr double w_rate = 25.0; // [rad/s] inner rate-loop bandwidth
inline constexpr double w_angle = 5.0; // [rad/s] outer angle-loop bandwidth
inline constexpr float  u_max = 12.0f;
inline constexpr float  w_cmd_max = 8.0f;

inline constexpr auto rate_pi_design = design::pid_pole_placement(TransferFunction{.num = {1.0}, .den = {1.0 / b_u, 0.0}}, {.omega = w_rate});
inline constexpr auto angle_pi_design = design::pid_pole_placement(TransferFunction{.num = {1.0}, .den = {1.0, 0.0}}, {.omega = w_angle});

/// The loops' output limits, which discretize gives each PI: the motor's ±u_max, the rate command's ±w_cmd_max.
inline constexpr design::PIDLimits<double> rate_pi_limits{.input = {.lower = -static_cast<double>(u_max), .upper = static_cast<double>(u_max)}};
inline constexpr design::PIDLimits<double> angle_pi_limits{.input = {.lower = -static_cast<double>(w_cmd_max), .upper = static_cast<double>(w_cmd_max)}};

// Light-model SISO on θ
inline constexpr double wc_adrc = 8.0;
inline constexpr double wo_adrc = 40.0;
inline constexpr auto   adrc_res = design::adrc<2>(wc_adrc, wo_adrc, b_u, Ts);
static_assert(adrc_res.success);

// Sliding Mode Controllers
inline constexpr auto smc_res = design::smc(10.0, 8.0, b_u, Ts, 0.05); // λ, k, b0, Ts, boundary layer φ
inline constexpr auto stsmc_res = design::stsmc(4.0, Ts, 10.0);

// Full-state linear design: the same design call, placed at s-plane poles or weighted
inline constexpr design::Poles place_poles_s{-8.0, -12.0};

inline constexpr auto place_res = design::state_feedback(plant, Ts, place_poles_s); // s-plane poles, sampled at Ts: an LQR placed
static_assert(place_res.success);

// Implement full state feedback controllers
inline constexpr auto lqr_res = design::state_feedback(sysd, lqr_weights);
inline constexpr auto lqi_res = design::state_feedback(sysd, lqr_weights, integral);
inline constexpr auto lqg_res = design::state_feedback(sysd, lqr_weights, estimator_noise);
inline constexpr auto lqgi_res = design::state_feedback(sysd, lqr_weights, estimator_noise, integral);

// Offset-free MPC — synthesis only
inline constexpr auto mpc_art = design::mpc<10, 4>(
    sysd,
    design::MPCWeights{
        .Qy = {{20.0}},
        .Rdu = {{0.1}},
    },
    design::InputLimits{
        .lower = {-static_cast<double>(u_max)},
        .upper = {static_cast<double>(u_max)},
    },
    design::KalmanNoise{
        .process = estimator_noise.process,
        .measurement = estimator_noise.measurement,
        .disturbance = 1e-2,
    },
    design::InputDisturbance{}
);

/**
 * One control period: outer θ → ω*, inner ω → u.
 * On target: θ from fused roll estimate, ω from gyro; write u to the wheel driver.
 */
[[nodiscard]] inline float cascade_period(
    PIController<float>& angle_pi,
    PIController<float>& rate_pi,
    float                theta_ref,
    float                theta,
    float                omega
) {
    const float omega_cmd = angle_pi.step(theta_ref, theta);
    return rate_pi.step(omega_cmd, omega);
}

} // namespace damp::examples_reaction_wheel
