// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file pendulum_controller.hpp
 * @brief Upright 2-state pendulum LQR — flashable design + tick
 *
 * Folder layout (Arduino-style example):
 *
 *   pendulum_controller.hpp  — this file (nameplate, design::, control_period)
 *   pendulum_sketch.cpp      — thin setup/loop smoke (float deploy)
 *   pendulum_sil.cpp         — host nonlinear plant + plots
 *   pendulum_derivation.md   — plant derivation (SS / TF)
 *
 * Copy this header + sketch onto target; replace read_state / write_torque.
 * Do not put plant, std::vector, or plotly here.
 */

#pragma once

#include "damp/controllers/lqr.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

namespace damp::examples_pendulum {

// Nameplate (linearization about upright equilibrium θ = 0)
inline constexpr double g = 9.81;     // gravity [m/s²]
inline constexpr double L = 1.0;      // pivot → bob length [m]
inline constexpr double m = 1.0;      // point mass [kg]
inline constexpr double b_damp = 0.1; // pivot viscous damping [N·m·s/rad]
inline constexpr double Ts = 0.01;    // 100 Hz control

/**
 * @brief Linearize pendulum at upright equilibrium (θ = 0)
 *
 * Continuous A, B for design; states [θ, θ̇], input torque u.
 */
[[nodiscard]] constexpr auto linearize_upright() {
    const double a21 = g / L;
    const double a22 = -b_damp / (m * L * L);

    Matrix<2, 2> A{{0.0, 1.0}, {a21, a22}};
    Matrix<2, 1> B{{0.0}, {1.0 / (m * L * L)}};
    Matrix<1, 2> C{{1.0, 0.0}};

    return StateSpace{.A = A, .B = B, .C = C};
}

inline constexpr auto sys = linearize_upright();
inline constexpr auto Q = Matrix<2, 2>::identity() * 10.0;
inline constexpr auto R = Matrix<1, 1>{{1.0}};

inline constexpr auto lqr_d = design::discrete_lqr_from_continuous(sys.A, sys.B, Q, R, Ts);
static_assert(lqr_d.success);

/**
 * @brief One control period: full-state feedback u = −Kx
 *
 * On target, replace x with encoder/IMU fusion. SIL may use double; sketch uses float.
 */
template<typename T>
[[nodiscard]] ColVec<1, T> control_period(LQR<2, 1, T>& controller, const ColVec<2, T>& x) {
    return controller.control(x);
}

} // namespace damp::examples_pendulum
