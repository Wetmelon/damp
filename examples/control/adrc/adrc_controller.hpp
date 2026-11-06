// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file adrc_controller.hpp
 * @brief Position ADRC — flashable design + tick (compare SIL uses the same objects)
 *
 * Folder layout:
 *
 *   adrc_controller.hpp  — this file (nameplate, ADRC + PI-D designs, ticks)
 *   adrc_sketch.cpp      — thin setup/loop smoke (float ADRC)
 *   adrc_sil.cpp         — host plant: ADRC vs well-tuned PI-D + plot
 *   adrc_derivation.md   — plant, why PI-D loses, gain map
 *
 * Product path is ADRC. The PI-D is a fair nominal-plant tune in the same
 * header so the SIL does not invent a second law.
 */

#pragma once

#include "damp/controllers/adrc.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/design/pid_design.hpp"
#include "damp/systems/transfer_function.hpp"

namespace damp::examples_adrc {

inline constexpr double Ts = 0.001;     // 1 kHz
inline constexpr double J_nom = 0.05;   // [kg·m²] nameplate inertia
inline constexpr double B_visc = 0.03;  // [N·m·s/rad]
inline constexpr double tau_c = 0.12;   // [N·m] Coulomb (smooth tanh)
inline constexpr double omega_c = 0.04; // [rad/s] tanh scale
inline constexpr double u_max = 2.5;    // [N·m]
inline constexpr double wc = 8.0;       // [rad/s] shared position bandwidth
inline constexpr double wo = 40.0;      // [rad/s] ESO (5× wc)

/// High-frequency gain the ADRC is told about: θ̈ ≈ b0·u + f.
inline constexpr double b0 = 1.0 / J_nom;

[[nodiscard]] constexpr design::ADRCResult<2, double> make_adrc() {
    return design::adrc<2>(wc, wo, b0, Ts);
}

inline constexpr auto adrc_d = make_adrc();

/// Critically-damped PD at wc plus an I pole at wc/8 (load hold).
[[nodiscard]] constexpr design::PIDGains<double> make_pid() {
    auto r = design::pid_pole_placement(TransferFunction<1, 3>{.num = {1.0}, .den = {J_nom, 0.0, 0.0}}, {.omega = wc, .zeta = 1.0, .omega_i = wc / 8.0});
    r.u_lim = {
        .lower = -u_max,
        .upper = u_max,
    };
    r.integral_lim = {
        .lower = -u_max,
        .upper = u_max,
    };
    return r;
}

inline constexpr auto pid_d = make_pid();

template<typename T>
[[nodiscard]] T control_period_adrc(ADRCController<2, T>& ctrl, T r, T y) {
    return ctrl.control(r, y);
}

template<typename T>
[[nodiscard]] T control_period_pid(PIDController<T>& ctrl, T r, T y) {
    return ctrl.control(r, y);
}

} // namespace damp::examples_adrc
