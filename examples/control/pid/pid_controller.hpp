// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file pid_controller.hpp
 * @brief PI deploy — flashable design (Design Is Deploy)
 *
 * Folder layout:
 *
 *   pid_controller.hpp  — this file (nameplate + design result)
 *   pid_sketch.cpp      — setup/loop mock I/O
 *   pid_sil.cpp         — thin host loop (same control API)
 *   pid_derivation.md   — plant notes
 *
 * Lean includes only — do not pull damp/control.hpp on a small MCU for this demo.
 */

#pragma once

#include "damp/control.hpp"

namespace damp::examples_pid {

inline constexpr float Ts = 0.01f; // 100 Hz

// Prefer: design → discretize → as<float> → assign into the runtime type.
// discretize in double so static constinit is a constant expression under GCC.
// For setpoint weights / filter / tracking time, name them:
// design::PIDGains{.Kp = 4.0, .Ki = 8.0, .b = 0.5}; limits go to discretize,
// .discretize(Ts, {.input = {.lower = -12.0, .upper = 12.0}}), and the law takes new gains and
// limits at run time with set_gains / set_limits.
inline constexpr auto pi_design = design::pid(4.0, 8.0, 0.0);

/**
 * @brief One control period: u = PI(r, y)
 */
template<typename T>
[[nodiscard]] T control_period(PIController<T>& controller, T r, T y) {
    return controller.step(r, y);
}

} // namespace damp::examples_pid
