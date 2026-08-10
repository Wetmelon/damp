// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file eskf_estimator.hpp
 * @brief MARG attitude ESKF — flashable design + tick (single source of truth)
 *
 * Folder layout (Arduino-style example):
 *
 *   eskf_estimator.hpp  — this file (nameplate, design::eskf_marg, estimate_period)
 *   eskf_sketch.cpp     — thin setup/loop smoke (float deploy)
 *   eskf_sil.cpp        — finite-tick host smoke
 *   eskf_derivation.md  — filter / measurement model notes
 *
 * Copy this header + sketch onto target; replace mock IMU/mag reads.
 * Do not put plant, std::vector, or plotly here.
 */

#pragma once

#include "damp/estimation/eskf.hpp"
#include "damp/estimation/sensor_fusion.hpp"
#include "damp/matrix/matrix.hpp"

namespace damp::examples_eskf {

// Nameplate / rates
inline constexpr float dt = 0.01f;       // 100 Hz
inline constexpr float kGyroNd = 0.003f; // rad/s/√Hz
inline constexpr float kAccelNd = 0.03f; // m/s²/√Hz
inline constexpr float kMagNd = 0.3f;    // mag noise density
inline constexpr float kBgRw = 0.0001f;  // rad/s²/√Hz-class gyro bias RW

// design::eskf_marg → constinit design payload (same DiD shape as LQR etc.)
inline constexpr auto kDesign = design::eskf_marg(kGyroNd, kAccelNd, kMagNd, kBgRw, dt);

/**
 * @brief One estimation period: MARG update → orientation
 *
 * On target, pass driver samples. SIL/sketch use mock sensors at rest.
 */
template<typename T>
void estimate_period(
    ESKFOrientationFilter<T, 6>& filt,
    const Vec3<T>&               accel,
    const Vec3<T>&               gyro,
    const Vec3<T>&               mag
) {
    filt.update(accel, gyro, mag, static_cast<T>(dt));
}

} // namespace damp::examples_eskf
