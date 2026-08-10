// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file cart_pole_controller.hpp
 * @brief Cart-pole LQR — flashable design + tick (single source of truth)
 *
 * Folder layout (Arduino-style example):
 *
 *   cart_pole_controller.hpp  — this file (nameplate, design::, control_period)
 *   cart_pole_sketch.cpp      — thin setup/loop smoke (float deploy)
 *   cart_pole_sil.cpp         — host nonlinear plant + plots
 *   cart_pole_derivation.md   — plant derivation (SS / TF)
 *
 * Copy this header + sketch onto target; replace read_state / write_force.
 * Do not put plant, std::vector, or plotly here.
 */

#pragma once

#include "damp/controllers/lqr.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

namespace damp::examples_cart_pole {

// Nameplate
inline constexpr double M = 1.0;      // cart mass [kg]
inline constexpr double m_pole = 0.1; // pole mass [kg]
inline constexpr double L = 0.5;      // pole half-length / CoM distance [m]
inline constexpr double g = 9.81;     // gravity [m/s²]
inline constexpr double b_fric = 0.1; // cart viscous friction [N·s/m]
inline constexpr double Ts = 0.01;    // 100 Hz control

/**
 * @brief Linearize cart-pole at upright equilibrium (θ = 0)
 *
 * Continuous A, B for design; matches the nonlinear plant Jacobian at rest.
 */
[[nodiscard]] constexpr auto linearize_upright() {
    // Jacobian of cart_pole_sil nonlinear ODE at (x,u)=0 (θ from upright).
    // ẍ gets gravity coupling −(m g / M) θ; θ̈ uses ẍ cosθ / ℓ so a43 = g/ℓ.
    const double denom = M;
    const double a22 = -b_fric / denom;
    const double a23 = -(m_pole * g) / denom;
    const double a42 = a22 / L; // = −b/(M ℓ)
    const double a43 = g / L;   // = (M+m)g/(M ℓ) + a23/ℓ
    const double b2 = 1.0 / denom;
    const double b4 = b2 / L; // = 1/(M ℓ)

    Matrix<4, 4> A{
        {0.0, 1.0, 0.0, 0.0},
        {0.0, a22, a23, 0.0},
        {0.0, 0.0, 0.0, 1.0},
        {0.0, a42, a43, 0.0}
    };
    ColVec<4>    B{0.0, b2, 0.0, b4};
    Matrix<2, 4> C{
        {1.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0}
    };
    return StateSpace{.A = A, .B = B, .C = C};
}

inline constexpr auto sys = linearize_upright();

inline constexpr auto Q = Matrix<4, 4>{
    {10.0, 0.0, 0.0, 0.0},
    {0.0, 0.1, 0.0, 0.0},
    {0.0, 0.0, 100.0, 0.0},
    {0.0, 0.0, 0.0, 10.0}
};

inline constexpr auto R = Matrix<1, 1>{{1.0}};

// 4×4 DARE at compile time — may need a higher -fconstexpr-ops-limit on some hosts
inline constexpr auto lqr_d = design::discrete_lqr_from_continuous(sys.A, sys.B, Q, R, Ts);
static_assert(lqr_d.success);

/**
 * @brief One control period: full-state feedback u = −Kx
 *
 * On target, replace x with encoder/IMU fusion. SIL may use double; sketch uses float.
 */
template<typename T>
[[nodiscard]] ColVec<1, T> control_period(LQR<4, 1, T>& controller, const ColVec<4, T>& x) {
    return controller.control(x);
}

} // namespace damp::examples_cart_pole
