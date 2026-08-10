// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file ins_eskf_sil.cpp
 * @brief Animated INS: free-run bias vs InsNavigator (position + dual-antenna heading)
 *
 * Host-only SIL of the design in ins_eskf_estimator.hpp (same densities as
 * estimation/ins_navigator/). IMU period matches deploy (100 Hz). Biased IMU;
 * absolute position + dual-antenna heading at 2 Hz on the aided path.
 *
 *   green  — truth path
 *   orange — open-loop mechanization (accel bias, no filter)
 *   cyan   — InsNavigator (position + heading aided)
 *   RGB    — body axes at the true vehicle pose
 *
 * Open plots/estimation/ins_eskf_3d.html → Play.
 */

#include <algorithm>
#include <cstddef>
#include <vector>

#include "animate_3d.hpp"
#include "damp/backend.hpp"
#include "damp/estimation/ins_eskf.hpp"
#include "damp/estimation/ins_mechanization.hpp"
#include "damp/math/geometry.hpp"
#include "damp/math/math.hpp"
#include "fmt/core.h"
#include "ins_eskf_estimator.hpp"

using namespace damp;
using namespace damp::examples_ins_eskf;
using namespace damp::examples_plot3d;

namespace {

[[nodiscard]] ImuSample<double>
synthetic_imu(const Quaternion<double>& q, const Vec3<double>& a_nav, const Vec3<double>& omega_body, NavFrame frame) {
    const Vec3<double> g_n = gravity_nav(frame);
    const Vec3<double> a_b = q.conjugate().rotate(a_nav - g_n);
    return ImuSample<double>{.gyro = omega_body, .accel = a_b};
}

void expand_bounds(
    double&          xmin,
    double&          xmax,
    double&          ymin,
    double&          ymax,
    double&          zmin,
    double&          zmax,
    const Polyline3& pl
) {
    for (std::size_t i = 0; i < pl.x.size(); ++i) {
        xmin = std::min(xmin, pl.x[i]);
        xmax = std::max(xmax, pl.x[i]);
        ymin = std::min(ymin, pl.y[i]);
        ymax = std::max(ymax, pl.y[i]);
        zmin = std::min(zmin, pl.z[i]);
        zmax = std::max(zmax, pl.z[i]);
    }
}

[[nodiscard]] Polyline3 path_prefix(const Polyline3& full, int k) {
    const auto n = static_cast<std::size_t>(k) + 1;
    return Polyline3{
        std::vector<double>(full.x.begin(), full.x.begin() + static_cast<std::ptrdiff_t>(n)),
        std::vector<double>(full.y.begin(), full.y.begin() + static_cast<std::ptrdiff_t>(n)),
        std::vector<double>(full.z.begin(), full.z.begin() + static_cast<std::ptrdiff_t>(n)),
    };
}

} // namespace

int main() {
    fmt::print("===== Animated InsNavigator (truth vs free-run vs position+heading) =====\n\n");

    constexpr double   t_end = 18.0;
    constexpr double   dt = kTs;
    const int          n = static_cast<int>(t_end / dt) + 1;
    constexpr double   axis_len = 0.8;
    constexpr int      aid_every = 50; // 2 Hz at 100 Hz IMU
    const Vec3<double> sensor_ba{0.05, 0.0, 0.0};

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

    InsState<double> free_run{};
    free_run.q = Quaternion<double>::identity();

    // Design from header (DiD parity with ins_navigator / sketch)
    InsNavigator<double> nav(kInsDesign, kFrame);

    const Matrix<1, 1, double> R_psi{{kHeadingStd * kHeadingStd}};

    Polyline3                       path_truth;
    Polyline3                       path_free;
    Polyline3                       path_aid;
    std::vector<Quaternion<double>> truth_q(static_cast<std::size_t>(n));

    double max_free = 0.0;
    double max_aid = 0.0;
    double max_aid_late = 0.0;

    for (int k = 0; k < n; ++k) {
        const double t = static_cast<double>(k) * dt;

        const Vec3<double> a_n = a_nav_true(t);
        const Vec3<double> w_b = omega_true(t);

        {
            const ImuSample<double> imu_ideal = synthetic_imu(truth.q, a_n, w_b, kFrame);
            truth = mechanize_step(truth, imu_ideal, dt, kFrame);
        }

        ImuSample<double> imu = synthetic_imu(truth.q, a_n, w_b, kFrame);
        imu.accel = imu.accel + sensor_ba;

        free_run = mechanize_step(free_run, imu, dt, kFrame);

        nav.predict(imu, dt);
        if ((k % aid_every) == 0) {
            (void)nav.update_pose_heading(truth.p, ins_heading(truth, kFrame), kInsDesign.R, R_psi);
        }

        truth_q[static_cast<std::size_t>(k)] = truth.q;
        path_truth.x.push_back(truth.p[0]);
        path_truth.y.push_back(truth.p[1]);
        path_truth.z.push_back(truth.p[2]);
        path_free.x.push_back(free_run.p[0]);
        path_free.y.push_back(free_run.p[1]);
        path_free.z.push_back(free_run.p[2]);
        path_aid.x.push_back(nav.state().p[0]);
        path_aid.y.push_back(nav.state().p[1]);
        path_aid.z.push_back(nav.state().p[2]);

        const auto err = [](const Vec3<double>& a, const Vec3<double>& b) {
            const double dx = a[0] - b[0];
            const double dy = a[1] - b[1];
            const double dz = a[2] - b[2];
            return damp::sqrt(dx * dx + dy * dy + dz * dz);
        };
        const double e_free = err(free_run.p, truth.p);
        const double e_aid = err(nav.state().p, truth.p);
        max_free = damp::max(max_free, e_free);
        max_aid = damp::max(max_aid, e_aid);
        if (t > 10.0) {
            max_aid_late = damp::max(max_aid_late, e_aid);
        }
    }

    double xmin = path_truth.x.front();
    double xmax = xmin;
    double ymin = path_truth.y.front();
    double ymax = ymin;
    double zmin = path_truth.z.front();
    double zmax = zmin;
    expand_bounds(xmin, xmax, ymin, ymax, zmin, zmax, path_truth);
    expand_bounds(xmin, xmax, ymin, ymax, zmin, zmax, path_aid);
    zmin = damp::min(zmin, -0.5);
    zmax = damp::max(zmax, 1.5);

    const Polyline3 ground = ground_rect(xmin, xmax, ymin, ymax);
    const SceneBox  box = fit_scene_box(xmin, xmax, ymin, ymax, zmin, zmax);

    std::vector<std::vector<Polyline3>> frames(static_cast<std::size_t>(n));
    std::vector<double>                 times(static_cast<std::size_t>(n));
    for (int k = 0; k < n; ++k) {
        times[static_cast<std::size_t>(k)] = static_cast<double>(k) * dt;
        const Vec3<double> p_t{
            path_truth.x[static_cast<std::size_t>(k)],
            path_truth.y[static_cast<std::size_t>(k)],
            path_truth.z[static_cast<std::size_t>(k)],
        };
        const Vec3<double> p_a{
            path_aid.x[static_cast<std::size_t>(k)],
            path_aid.y[static_cast<std::size_t>(k)],
            path_aid.z[static_cast<std::size_t>(k)],
        };
        const auto ax = body_axes_at(truth_q[static_cast<std::size_t>(k)], p_t, axis_len);
        frames[static_cast<std::size_t>(k)] = {
            ground,
            path_prefix(path_truth, k),
            path_prefix(path_free, k),
            path_prefix(path_aid, k),
            ax[0],
            ax[1],
            ax[2],
            Polyline3{{p_t[0]}, {p_t[1]}, {p_t[2]}},
            Polyline3{{p_a[0]}, {p_a[1]}, {p_a[2]}},
        };
    }

    const std::vector<TraceStyle> styles = {
        {.name = "ground", .color = "#555555", .mode = "lines", .line_width = 2, .dash = "solid", .static_geometry = true},
        {.name = "path (truth)", .color = "#2ca02c", .mode = "lines", .line_width = 5, .dash = "solid"},
        {.name = "path (free-run, accel bias)", .color = "#ff7f0e", .mode = "lines", .line_width = 5, .dash = "dot"},
        {.name = "path (InsNavigator p+ψ)", .color = "#17becf", .mode = "lines", .line_width = 5, .dash = "solid"},
        {.name = "body x", .color = "#d62728", .mode = "lines", .line_width = 6, .dash = "solid"},
        {.name = "body y", .color = "#2ca02c", .mode = "lines", .line_width = 6, .dash = "solid"},
        {.name = "body z", .color = "#1f77b4", .mode = "lines", .line_width = 6, .dash = "solid"},
        {.name = "vehicle (truth)", .color = "#ffffff", .mode = "markers", .line_width = 2, .marker_size = 7, .dash = "solid"},
        {.name = "vehicle (aided)", .color = "#17becf", .mode = "markers", .line_width = 2, .marker_size = 6, .dash = "solid"},
    };

    write_animated_scatter3d(
        "plots/estimation/ins_eskf_3d.html",
        styles,
        frames,
        times,
        dt,
        "InsNavigator — truth (green) vs free-run bias (orange) vs position+heading (cyan)",
        box
    );

    fmt::print("Frames: {}  → plots/estimation/ins_eskf_3d.html\n", n);
    fmt::print(
        "Scene (truth+aided) XY [{:.1f},{:.1f}]×[{:.1f},{:.1f}]  camera target ({:.1f},{:.1f},{:.1f})\n",
        box.x[0],
        box.x[1],
        box.y[0],
        box.y[1],
        box.center_x,
        box.center_y,
        box.center_z
    );
    fmt::print(
        "Max |p| error free-run: {:.2f} m   aided peak: {:.2f} m   aided after t>10s: {:.3f} m\n",
        max_free,
        max_aid,
        max_aid_late
    );
    fmt::print(
        "Heading truth: {:.3f} rad   nav: {:.3f} rad\n"
        "b_a estimate: ({:.4f}, {:.4f}, {:.4f})  true body-x bias: {:.2f}\n",
        ins_heading(truth, kFrame),
        nav.heading(),
        nav.state().b_a[0],
        nav.state().b_a[1],
        nav.state().b_a[2],
        sensor_ba[0]
    );
    return 0;
}
