// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file ins_mechanization_sil.cpp
 * @brief Animated strapdown INS path: IMU samples → p, v, q (mechanization)
 *
 * Host-only. Synthesizes a simple vehicle motion in ENU, builds the matching
 * ideal IMU stream, and runs mechanize_step via ins_mechanization_estimator.hpp.
 *
 * Open plots/estimation/ins_mechanization_3d.html → Play.
 */

#include <array>
#include <cstddef>
#include <vector>

#include "animate_3d.hpp"
#include "damp/backend.hpp"
#include "damp/estimation/ins_mechanization.hpp"
#include "damp/math/geometry.hpp"
#include "damp/math/math.hpp"
#include "fmt/base.h"
#include "fmt/core.h"
#include "ins_mechanization_estimator.hpp"

using namespace damp;
using namespace damp::examples_ins_mechanization;
using namespace damp::examples_plot3d;

namespace {

[[nodiscard]] Polyline3 ground_grid() {
    return Polyline3{
        {-2.0, 12.0, 12.0, -2.0, -2.0},
        {-4.0, -4.0, 4.0, 4.0, -4.0},
        {0.0, 0.0, 0.0, 0.0, 0.0},
    };
}

} // namespace

int main() {
    fmt::print("===== Animated INS mechanization (path + body axes + bias drift) =====\n\n");

    constexpr double t_end = 12.0;
    constexpr double dt = kTsSil;
    const int        n = static_cast<int>(t_end / dt) + 1;
    constexpr double axis_len = 0.8;

    auto a_nav_true = [](double t) -> Vec3<double> {
        if (t < 3.0) {
            return {1.2, 0.0, 0.0};
        }
        return {0.0, 0.0, 0.0};
    };
    auto omega_true = [](double t) -> Vec3<double> {
        if (t >= 4.0 && t < 7.0) {
            return {0.0, 0.0, 0.35};
        }
        return {0.0, 0.0, 0.0};
    };

    InsState<double> truth{};
    truth.q = Quaternion<double>::identity();

    InsState<double> mech{};
    mech.q = Quaternion<double>::identity();

    InsState<double> mech_biased{};
    mech_biased.q = Quaternion<double>::identity();
    const Vec3<double> sensor_ba{0.05, 0.0, 0.0};

    std::vector<std::vector<Polyline3>> frames(static_cast<size_t>(n));
    std::vector<double>                 times(static_cast<size_t>(n));

    Polyline3  path_truth;
    Polyline3  path_mech;
    Polyline3  path_biased;
    const auto ground = ground_grid();

    for (int k = 0; k < n; ++k) {
        const double t = static_cast<double>(k) * dt;
        times[static_cast<size_t>(k)] = t;

        const Vec3<double> a_n = a_nav_true(t);
        const Vec3<double> w_b = omega_true(t);

        {
            const ImuSample<double> imu_ideal = synthetic_imu(truth.q, a_n, w_b, kFrame);
            truth = estimate_period(truth, imu_ideal, dt, kFrame);
        }

        {
            const ImuSample<double> imu = synthetic_imu(mech.q, a_n, w_b, kFrame);
            mech = estimate_period(mech, imu, dt, kFrame);
        }

        {
            ImuSample<double> imu = synthetic_imu(mech_biased.q, a_n, w_b, kFrame);
            imu.accel = imu.accel + sensor_ba;
            mech_biased = estimate_period(mech_biased, imu, dt, kFrame);
        }

        path_truth.x.push_back(truth.p[0]);
        path_truth.y.push_back(truth.p[1]);
        path_truth.z.push_back(truth.p[2]);
        path_mech.x.push_back(mech.p[0]);
        path_mech.y.push_back(mech.p[1]);
        path_mech.z.push_back(mech.p[2]);
        path_biased.x.push_back(mech_biased.p[0]);
        path_biased.y.push_back(mech_biased.p[1]);
        path_biased.z.push_back(mech_biased.p[2]);

        const auto ax = body_axes_at(truth.q, truth.p, axis_len);

        frames[static_cast<size_t>(k)] = {
            ground,
            path_truth,
            path_biased,
            ax[0],
            ax[1],
            ax[2],
            Polyline3{{truth.p[0]}, {truth.p[1]}, {truth.p[2]}},
        };
    }

    const std::vector<TraceStyle> styles = {
        {.name = "ground", .color = "#555555", .mode = "lines", .line_width = 2, .dash = "solid", .static_geometry = true},
        {.name = "path (ideal mechanization)", .color = "#2ca02c", .mode = "lines", .line_width = 5, .dash = "solid"},
        {.name = "path (open-loop, accel bias)", .color = "#ff7f0e", .mode = "lines", .line_width = 5, .dash = "dot"},
        {.name = "ideal body x", .color = "#d62728", .mode = "lines", .line_width = 6, .dash = "solid"},
        {.name = "ideal body y", .color = "#2ca02c", .mode = "lines", .line_width = 6, .dash = "solid"},
        {.name = "ideal body z", .color = "#1f77b4", .mode = "lines", .line_width = 6, .dash = "solid"},
        {.name = "vehicle", .color = "#ffffff", .mode = "markers", .line_width = 2, .marker_size = 6, .dash = "solid"},
    };

    SceneBox box;
    box.x = {-3.0, 14.0};
    box.y = {-5.0, 6.0};
    box.z = {-1.0, 4.0};
    box.eye_x = 1.8;
    box.eye_y = -2.2;
    box.eye_z = 1.4;
    box.center_x = 4.0;
    box.center_y = 0.0;
    box.center_z = 0.3;

    write_animated_scatter3d(
        "plots/estimation/ins_mechanization_3d.html",
        styles,
        frames,
        times,
        dt,
        "INS mechanization — ideal path (solid green) vs open-loop accel bias (dotted orange)",
        box
    );

    fmt::print("Frames: {}  → plots/estimation/ins_mechanization_3d.html\n", n);
    fmt::print("Solid green = ideal mechanization; dotted orange = open-loop with accel bias (not a filter).\n");
    fmt::print("RGB triad = ideal body axes at the vehicle (body → ENU).\n");
    return 0;
}
