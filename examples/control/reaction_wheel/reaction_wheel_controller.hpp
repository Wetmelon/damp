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
 *   3. place / LQR / LQI           — full-state linear design
 *   4. LQG / LQGI                  — state feedback + Kalman from y
 *   5. Offset-free MPC             — constrained predictive (synthesis only here)
 *
 * Do not put plant ODE, std::vector, or plotly here.
 */

#pragma once

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/controllers/adrc.hpp"
#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/controllers/ofmpc.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/controllers/smc.hpp"
#include "damp/controllers/stsmc.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/pid_design.hpp"
#include "damp/design/pole_placement.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/math/complex.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

namespace damp::examples_reaction_wheel {

inline constexpr double Ts = 0.001; // [s]  1 kHz control loop
inline constexpr double w0 = 1.0;   // [rad/s] open-loop unstable natural frequency
inline constexpr double b_u = 1.0;  // [rad/s² per N·m]  ω̇ = … + b_u·u

inline constexpr size_t NX = 2;
inline constexpr size_t NU = 1;
inline constexpr size_t NY = 1;

// Continuous plant + identity noise channels so Kalman-based designs have a
// well-posed process/measurement model (G, H = I).
inline constexpr StateSpace<NX, NU, NY, double, NX, NY> plant{
    .A = {{0.0, 1.0}, {w0 * w0, 0.0}},
    .B = {{0.0}, {b_u}},
    .C = {{1.0, 0.0}},
    .D = {},
    .G = Matrix<NX, NX>::identity(),
    .H = Matrix<NY, NY>::identity(),
    .Ts = 0.0,
};

inline constexpr auto sysd = *discretize(plant, Ts, DiscretizationMethod::ZOH);

// Shared LQ* weights
inline constexpr auto Q = Matrix<NX, NX>::diagonal({40.0, 4.0});
inline constexpr auto R = Matrix<NU, NU>{{1.0}};
inline constexpr auto Q_integral = Matrix<NY, NY>{{80.0}};
inline constexpr auto Q_kf = Matrix<NX, NX>::diagonal({1e-4, 1e-3});
inline constexpr auto R_kf = Matrix<NY, NY>{{1e-3}};

// Cascade PI
inline constexpr double w_rate = 25.0; // [rad/s] inner rate-loop bandwidth
inline constexpr double w_angle = 5.0; // [rad/s] outer angle-loop bandwidth
inline constexpr float  u_max = 12.0f;
inline constexpr float  w_cmd_max = 8.0f;

[[nodiscard]] constexpr design::PIDGains<double> design_rate_pi() {
    auto r = design::pid_pole_placement(TransferFunction<1, 2>{.num = {1.0}, .den = {1.0 / b_u, 0.0}}, {.omega = w_rate});
    r.u_lim = {
        .lower = -static_cast<double>(u_max),
        .upper = static_cast<double>(u_max),
    };
    return r;
}

[[nodiscard]] constexpr design::PIDGains<double> design_angle_pi() {
    auto r = design::pid_pole_placement(TransferFunction<1, 2>{.num = {1.0}, .den = {1.0, 0.0}}, {.omega = w_angle});
    r.u_lim = {
        .lower = -static_cast<double>(w_cmd_max),
        .upper = static_cast<double>(w_cmd_max),
    };
    return r;
}

inline constexpr auto rate_pi_design = design_rate_pi();
inline constexpr auto angle_pi_design = design_angle_pi();

// Light-model SISO on θ
inline constexpr double wc_adrc = 8.0;
inline constexpr double wo_adrc = 40.0;
inline constexpr auto   adrc_res = design::adrc<2>(wc_adrc, wo_adrc, b_u, Ts);
static_assert(adrc_res.success);

// Sliding Mode Controllers
inline constexpr auto smc_res = design::smc(10.0, 8.0, b_u, Ts, 0.05); // λ, k, b0, Ts, boundary layer φ
inline constexpr auto stsmc_res = design::stsmc(4.0, Ts, 10.0);

// Full-state linear design
inline constexpr damp::array<damp::complex<double>, NX> place_poles_s{-8.0, -12.0};

inline constexpr auto place_res = design::place(plant, place_poles_s, Ts); // s-plane poles, sampled at Ts
static_assert(place_res.success);

// Implement full state feedback controllers
inline constexpr auto lqr_res = design::lqr(sysd, {.Q = Q, .R = R});
inline constexpr auto lqi_res = design::lqi(sysd, {.Q = Q, .R = R, .Q_integral = Q_integral});

// Estimator noise in <units> (per tick?  Van Loan?)
inline constexpr auto estimator_noise = design::KalmanNoise{.process = Q_kf, .measurement = R_kf};
inline constexpr auto lqg_res = design::lqg(sysd, {.Q = Q, .R = R}, estimator_noise);
inline constexpr auto lqgi_res = design::lqgi(sysd, {.Q = Q, .R = R, .Q_integral = Q_integral}, estimator_noise);

// Offset-free MPC — synthesis only
inline constexpr design::MPCWeights<NU, NY> mpc_weights{
    .Qy = Matrix<NY, NY>{{20.0}},
    .Rdu = Matrix<NU, NU>{{0.1}},
};

inline constexpr design::MPCConstraints<NU, NY> mpc_limits{
    .u_lim = {
        .lower = -static_cast<double>(u_max),
        .upper = static_cast<double>(u_max),
    },
};

inline constexpr auto mpc_art = design::ofmpc<10, 4>(
    sysd, mpc_weights, mpc_limits, {.process = Q_kf, .measurement = R_kf, .disturbance = 1e-2}, design::InputDisturbance{}
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
    const float omega_cmd = angle_pi.control(theta_ref, theta);
    return rate_pi.control(omega_cmd, omega);
}

} // namespace damp::examples_reaction_wheel
