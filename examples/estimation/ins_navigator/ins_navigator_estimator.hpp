// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file ins_navigator_estimator.hpp
 * @brief 15-state INS navigator — flashable design + tick (single source of truth)
 *
 * Folder layout (Arduino-style product demo):
 *
 *   ins_navigator_estimator.hpp  — this file (nameplate, design::ins_eskf_design, estimate_period)
 *   ins_navigator_sketch.cpp     — thin setup/loop smoke (float deploy)
 *   ins_navigator_sil.cpp        — finite-tick host smoke
 *   ins_navigator_derivation.md  — design / deploy model
 *
 * Same sensor densities as estimation/ins_eskf/ (DiD parity). Copy header + sketch
 * onto target; replace mock IMU/GPS/heading. Do not put plant or plotly here.
 *
 * @see estimation/ins_eskf/ for the 3D SIL of this design
 */

#pragma once

#include "damp/estimation/ins_eskf.hpp"
#include "damp/estimation/ins_mechanization.hpp"
#include "damp/math/geometry.hpp"
#include "damp/matrix/matrix.hpp"

namespace damp::examples_ins_navigator {

// Nameplate / rates (shared with ins_eskf DiD SIL)
inline constexpr float    kTs = 0.01f; // 100 Hz IMU
inline constexpr NavFrame kFrame = NavFrame::ENU;
inline constexpr float    kGyroNd = 0.003f; // rad/s/√Hz
inline constexpr float    kAccelNd = 0.03f; // m/s²/√Hz
inline constexpr float    kBgRw = 0.0001f;  // rad/s^{3/2} (≡ rad/s²/√Hz)
inline constexpr float    kBaRw = 0.001f;   // m/s^{5/2} class (≡ m/s³/√Hz)
inline constexpr float    kPosStd = 1.0f;   // m (position aid)

// design::ins_eskf_design → float design payload (Q/R/P0 baked at compile time)
inline constexpr auto kInsDesign = design::ins_eskf_design(kGyroNd, kAccelNd, kBgRw, kBaRw, kTs, kPosStd);

// Optional dual-antenna baseline along body +x [m]
inline constexpr Vec3<float> kBaselineBody{1.0f, 0.0f, 0.0f};

// Heading aid 1-σ (~1°); R_ψ = kHeadingStd² in estimate_period
inline constexpr float kHeadingStd = 0.02f;

/**
 * @brief One estimation period: IMU predict + optional absolute aids
 *
 * Call from a 100 Hz ISR / RTOS task. No heap. Pass nullptr-style flags via bools.
 */
template<typename T>
void estimate_period(
    InsNavigator<T>&    nav,
    const ImuSample<T>& imu,
    bool                have_gps,
    const Vec3<T>&      p_gps,
    bool                have_heading,
    T                   psi,
    const Vec3<T>&      baseline_body
) {
    nav.predict(imu, static_cast<T>(kTs));
    if (have_gps) {
        (void)nav.update_position(p_gps);
    }
    if (have_heading) {
        const Matrix<1, 1, T> R_psi{{static_cast<T>(kHeadingStd) * static_cast<T>(kHeadingStd)}};
        (void)nav.update_heading(psi, R_psi, baseline_body);
    }
}

} // namespace damp::examples_ins_navigator
