// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file imu_pose_sil.cpp
 * @brief Animated body-frame triad from gyro integration (IMU attitude only)
 *
 * Host-only. Orientation dead-reckoning — no filter, no position.
 * Uses imu_pose_estimator.hpp for rates / bias / integrate tick.
 *
 * Open plots/estimation/imu_pose_3d.html → Play.
 */

#include <array>
#include <cstddef>
#include <vector>

#include "animate_3d.hpp"
#include "damp/backend.hpp"
#include "damp/math/geometry.hpp"
#include "damp/math/math.hpp"
#include "fmt/base.h"
#include "fmt/core.h"
#include "imu_pose_estimator.hpp"

using namespace damp;
using namespace damp::examples_imu_pose;
using namespace damp::examples_plot3d;

int main() {
    fmt::print("===== Animated IMU pose (ideal body triad vs biased open-loop) =====\n\n");

    constexpr double t_end = 10.0;
    constexpr double dt = kTsSil;
    const int        n = static_cast<int>(t_end / dt) + 1;
    constexpr double axis_len = 0.45;
    constexpr double world_len = 0.55;

    const Vec3<double> b_g = gyro_bias<double>();

    Quaternion<double> q_ideal = Quaternion<double>::identity();
    Quaternion<double> q_open_loop = Quaternion<double>::identity();

    std::vector<std::vector<Polyline3>> frames(static_cast<size_t>(n));
    std::vector<double>                 times(static_cast<size_t>(n));

    const auto w_ax = world_axes(world_len);

    for (int k = 0; k < n; ++k) {
        const double t = static_cast<double>(k) * dt;
        times[static_cast<size_t>(k)] = t;

        const auto b_ideal = body_axes_at(q_ideal, Vec3<double>{}, axis_len);
        const auto b_ol = body_axes_at(q_open_loop, Vec3<double>{}, axis_len * 0.92);

        frames[static_cast<size_t>(k)] = {
            w_ax[0],
            w_ax[1],
            w_ax[2],
            b_ideal[0],
            b_ideal[1],
            b_ideal[2],
            b_ol[0],
            b_ol[1],
            b_ol[2],
        };

        const Vec3<double> w = omega_true(t);
        q_ideal = estimate_period(q_ideal, w, dt);
        q_open_loop = estimate_period(q_open_loop, w + b_g, dt);
    }

    const std::vector<TraceStyle> styles = {
        {.name = "world x (fixed)", .color = "#aa4444", .mode = "lines", .line_width = 3, .dash = "longdash", .static_geometry = true},
        {.name = "world y (fixed)", .color = "#44aa44", .mode = "lines", .line_width = 3, .dash = "longdash", .static_geometry = true},
        {.name = "world z (fixed)", .color = "#4444aa", .mode = "lines", .line_width = 3, .dash = "longdash", .static_geometry = true},
        {.name = "ideal body x", .color = "#ff2222", .mode = "lines", .line_width = 8, .dash = "solid"},
        {.name = "ideal body y", .color = "#22cc22", .mode = "lines", .line_width = 8, .dash = "solid"},
        {.name = "ideal body z", .color = "#4488ff", .mode = "lines", .line_width = 8, .dash = "solid"},
        {.name = "open-loop body x (ω+bias)", .color = "#ff2222", .mode = "lines", .line_width = 5, .dash = "dot"},
        {.name = "open-loop body y (ω+bias)", .color = "#22cc22", .mode = "lines", .line_width = 5, .dash = "dot"},
        {.name = "open-loop body z (ω+bias)", .color = "#4488ff", .mode = "lines", .line_width = 5, .dash = "dot"},
    };

    SceneBox box;
    box.x = {-0.7, 0.7};
    box.y = {-0.7, 0.7};
    box.z = {-0.7, 0.7};
    box.eye_x = 1.6;
    box.eye_y = 1.4;
    box.eye_z = 1.2;

    write_animated_scatter3d(
        "plots/estimation/imu_pose_3d.html",
        styles,
        frames,
        times,
        dt,
        "IMU attitude — ideal body (solid) vs open-loop gyro+bias (dotted); world fixed",
        box
    );

    fmt::print("Frames: {}  → plots/estimation/imu_pose_3d.html\n", n);
    fmt::print("World = fixed reference frame. Ideal = integrate true ω. Dotted = integrate(ω+bias), not a filter.\n");
    return 0;
}
