// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file foc_controller.hpp
 * @brief PMSM FOC current-loop setup (host teaching; prefer foc_switching for DiD)
 *
 * Shared nameplate + FOController factory used by foc_sil (PI vs I-P plots)
 * and foc_sketch (thin smoke). Not the product Design-Is-Deploy path —
 * see motor/foc_switching/.
 */

#pragma once

#include "damp/backend.hpp"
#include "damp/math/transforms.hpp"
#include "damp/motor/foc.hpp"
#include "damp/motor/spm.hpp"

namespace damp::examples_foc {

using namespace damp::motor;
using T = double;

// ---- Motor nameplate (Turnigy D5065 270KV surface-PMSM, Ld = Lq) -----------
inline constexpr T Rs = 0.039;      // [ohm] phase resistance (phase-neutral)
inline constexpr T Ls = 16e-6;      // [H]   phase inductance (Ld = Lq, phase-neutral)
inline constexpr T Kt_spec = 0.031; // [Nm/A] torque constant (amplitude / peak per-phase)
inline constexpr T pole_pairs = 7.0;
inline constexpr T Vdc = 24.0;      // [V]  DC bus
inline constexpr T omega_e = 600.0; // [erad/s] fixed electrical speed

// ---- Loop timing -----------------------------------------------------------
inline constexpr T   fsw = 20000.0;  // [Hz] PWM / current-loop rate
inline constexpr T   Ts = 1.0 / fsw; // [s]  control period
inline constexpr int substeps = 100; // plant Euler sub-steps per control tick
inline constexpr T   t_end = 6.0e-3; // [s]

// omega_bw ~ fsw/13 keeps the continuous pole placement valid (see
// damp::design::pi_pole_placement_first_order's sampling note).
inline constexpr T bw = 2.0 * numbers::pi_v<T> * fsw / 13.0;

// ---- Scenario --------------------------------------------------------------
inline constexpr T Iq_ref = 8.0;    // [A]   torque-current step target
inline constexpr T t_step = 0.5e-3; // [s]   reference step time
inline constexpr T t_dist = 3.0e-3; // [s]   disturbance onset
inline constexpr T vq_dist = -1.0;  // [V]   unmodeled q-axis voltage disturbance

// Plant flux linkage recovered from the datasheet Kt.
inline constexpr T lambda_pm = damp::motor::spm::flux_from_torque_constant(pole_pairs, Kt_spec);
inline constexpr T Vmax = damp::motor::voltage_circle_radius(Vdc);

// Inductance in DQ frame
inline constexpr DirectQuadrature<T> Ldq{Ls, Ls};

} // namespace damp::examples_foc
