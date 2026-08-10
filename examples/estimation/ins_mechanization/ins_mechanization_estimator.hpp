// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file ins_mechanization_estimator.hpp
 * @brief Strapdown mechanization nameplate + open-loop tick
 *
 * Folder layout:
 *
 *   ins_mechanization_estimator.hpp  — this file (frame, rates, mechanize tick)
 *   ins_mechanization_sketch.cpp     — thin smoke at rest
 *   ins_mechanization_sil.cpp        — primary host 3D animation
 *   ins_mechanization_derivation.md  — open-loop model notes
 *
 * No filter — IMU → p,v,q by dead-reckoning. Bias peel-away motivates
 * estimation/ins_eskf/.
 */

#pragma once

#include "damp/estimation/ins_mechanization.hpp"
#include "damp/math/geometry.hpp"
#include "damp/matrix/matrix.hpp"

namespace damp::examples_ins_mechanization {

inline constexpr NavFrame kFrame = NavFrame::ENU;
// Sketch uses a convenient IMU period; SIL animation uses 1/60 for smooth frames
inline constexpr double kTsSketch = 0.01;
inline constexpr double kTsSil = 1.0 / 60.0;

/** Ideal body specific force for known nav accel and orientation. */
template<typename T>
[[nodiscard]] ImuSample<T>
synthetic_imu(const Quaternion<T>& q, const Vec3<T>& a_nav, const Vec3<T>& omega_body, NavFrame frame) {
    const Vec3<T> g_n = gravity_nav(frame);
    const Vec3<T> a_b = q.conjugate().rotate(a_nav - g_n);
    return ImuSample<T>{.gyro = omega_body, .accel = a_b};
}

/**
 * @brief One open-loop mechanization period
 */
template<typename T>
[[nodiscard]] InsState<T> estimate_period(
    const InsState<T>&  x,
    const ImuSample<T>& imu,
    T                   dt,
    NavFrame            frame = kFrame
) {
    return mechanize_step(x, imu, dt, frame);
}

} // namespace damp::examples_ins_mechanization
