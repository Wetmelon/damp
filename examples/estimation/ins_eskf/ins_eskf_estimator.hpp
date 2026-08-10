// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file ins_eskf_estimator.hpp
 * @brief InsNavigator design nameplate — shared by 3D SIL and thin sketch
 *
 * Folder layout:
 *
 *   ins_eskf_estimator.hpp  — this file (sensor densities, design::ins_eskf_design)
 *   ins_eskf_sketch.cpp     — thin navigator smoke (same design, float deploy)
 *   ins_eskf_sil.cpp        — primary host 3D animation (truth / free-run / aided)
 *   ins_eskf_derivation.md  — plant / filter notes
 *
 * Densities match estimation/ins_navigator/ (DiD parity). Flashable product path
 * for deploy is also documented there; this folder’s sketch is a smoke twin.
 */

#pragma once

#include "damp/estimation/ins_eskf.hpp"
#include "damp/estimation/ins_mechanization.hpp"
#include "damp/math/geometry.hpp"
#include "damp/matrix/matrix.hpp"

namespace damp::examples_ins_eskf {

// Same IMU period / densities as ins_navigator (DiD contract)
inline constexpr double   kTs = 0.01; // 100 Hz
inline constexpr NavFrame kFrame = NavFrame::ENU;
inline constexpr double   kGyroNd = 0.003;
inline constexpr double   kAccelNd = 0.03;
inline constexpr double   kBgRw = 0.0001;
inline constexpr double   kBaRw = 0.001;
inline constexpr double   kPosStd = 1.0;
inline constexpr double   kHeadingStd = 0.02; // 1-σ ~1°; R_ψ = kHeadingStd² at call site

inline constexpr auto kInsDesign = design::ins_eskf_design(kGyroNd, kAccelNd, kBgRw, kBaRw, kTs, kPosStd);

inline constexpr Vec3<double> kBaselineBody{1.0, 0.0, 0.0};

/**
 * @brief One period: predict + optional position/heading aids (deploy-shaped)
 */
template<typename T>
void estimate_period(
    InsNavigator<T>&       nav,
    const ImuSample<T>&    imu,
    bool                   have_pose_heading,
    const Vec3<T>&         p_meas,
    T                      psi_meas,
    const Matrix<3, 3, T>& R_pos,
    const Matrix<1, 1, T>& R_psi
) {
    nav.predict(imu, static_cast<T>(kTs));
    if (have_pose_heading) {
        (void)nav.update_pose_heading(p_meas, psi_meas, R_pos, R_psi);
    }
}

} // namespace damp::examples_ins_eskf
