// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file imu_pose_estimator.hpp
 * @brief Open-loop gyro attitude nameplate (teaching — no filter)
 *
 * Folder layout:
 *
 *   imu_pose_estimator.hpp  — this file (rates, bias, integrate tick)
 *   imu_pose_sketch.cpp     — thin smoke
 *   imu_pose_sil.cpp        — primary host animation
 *   imu_pose_derivation.md  — kinematics notes
 *
 * Drift under constant gyro bias is why MARG ESKF / dual-antenna heading exist
 * (estimation/eskf/, estimation/ins_eskf/).
 */

#pragma once

#include "damp/math/geometry.hpp"
#include "damp/matrix/matrix.hpp"

namespace damp::examples_imu_pose {

inline constexpr double kTsSketch = 0.01;
inline constexpr double kTsSil = 1.0 / 60.0;

// Constant yaw bias on the raw gyro [rad/s]
inline constexpr double kBgZ = 0.05;

/**
 * @brief Integrate body rates one step (open-loop attitude)
 */
template<typename T>
[[nodiscard]] Quaternion<T> estimate_period(const Quaternion<T>& q, const Vec3<T>& omega, T dt) {
    return q.integrate_body_rates(omega, dt);
}

/** True body rate schedule used by the SIL [rad/s]. */
template<typename T>
[[nodiscard]] constexpr Vec3<T> omega_true(T t) {
    if (t < T{3}) {
        return {T{0}, T{0}, static_cast<T>(0.6)};
    }
    if (t < T{6}) {
        return {T{0}, static_cast<T>(0.5), T{0}};
    }
    if (t < T{9}) {
        return {static_cast<T>(0.4), T{0}, T{0}};
    }
    return {T{0}, T{0}, T{0}};
}

template<typename T>
[[nodiscard]] constexpr Vec3<T> gyro_bias() {
    return {T{0}, T{0}, static_cast<T>(kBgZ)};
}

} // namespace damp::examples_imu_pose
