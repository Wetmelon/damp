// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file ins_navigator_sil.cpp
 * @brief 15-state INS navigator — finite-tick host smoke (calls ins_navigator_estimator.hpp)
 *
 * Host only. Flashable path: ins_navigator_sketch.cpp. Full 3D animation SIL:
 * estimation/ins_eskf/ins_eskf_sil.cpp (same sensor densities).
 */

#include "damp/estimation/ins_eskf.hpp"
#include "damp/estimation/ins_mechanization.hpp"
#include "damp/math/geometry.hpp"
#include "damp/matrix/matrix.hpp"
#include "fmt/core.h"
#include "ins_navigator_estimator.hpp"

using namespace damp;
using namespace damp::examples_ins_navigator;

int main() {
    fmt::print("===== InsNavigator (host finite-tick SIL) =====\n\n");
    fmt::print(
        "Ts={:.0f} Hz, frame=ENU, gyro_nd={:.4f}, accel_nd={:.3f}, pos_std={:.1f} m\n\n",
        1.0f / kTs,
        kGyroNd,
        kAccelNd,
        kPosStd
    );

    InsNavigator<double> nav{kInsDesign.as<double>(), kFrame};
    const Vec3<double>   baseline{1.0, 0.0, 0.0};

    constexpr int n = 200;
    for (int i = 0; i < n; ++i) {
        const ImuSample<double> imu{
            .gyro = {0.0, 0.0, 0.0},
            .accel = specific_force_at_rest(nav.state().q, kFrame),
        };
        const Vec3<double> p_gps = nav.state().p;
        const double       psi = nav.heading(baseline);
        estimate_period(nav, imu, true, p_gps, true, psi, baseline);
    }

    const auto& s = nav.state();
    fmt::print("After {} ticks at rest (mock GPS hold):\n", n);
    fmt::print("  p = ({:.4f}, {:.4f}, {:.4f}) m\n", s.p[0], s.p[1], s.p[2]);
    fmt::print("  v = ({:.4f}, {:.4f}, {:.4f}) m/s\n", s.v[0], s.v[1], s.v[2]);
    fmt::print("  heading = {:.4f} rad\n", nav.heading(baseline));
    return 0;
}
