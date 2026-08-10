// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file eskf_sil.cpp
 * @brief MARG attitude ESKF — finite-tick host smoke (calls eskf_estimator.hpp)
 *
 * Host only. Flashable path: eskf_sketch.cpp + eskf_estimator.hpp.
 * Rest-frame mock IMU/mag; prints Euler after a short run.
 */

#include "damp/estimation/sensor_fusion.hpp"
#include "damp/math/geometry.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "eskf_estimator.hpp"
#include "fmt/core.h"

using namespace damp;
using namespace damp::examples_eskf;

int main() {
    fmt::print("===== MARG ESKF (host finite-tick SIL) =====\n\n");
    fmt::print("dt={:.0f} Hz, gyro_nd={:.4f}, accel_nd={:.3f}, mag_nd={:.2f}\n\n", 1.0f / dt, kGyroNd, kAccelNd, kMagNd);

    // Same design as sketch; double for host numerics
    ESKFOrientationFilter<double, 6> filt{kDesign.as<double>()};

    const Vec3<double> accel{0.0, 0.0, 9.81};
    const Vec3<double> gyro{0.0, 0.0, 0.0};
    const Vec3<double> mag{0.0, 1.0, 0.0};

    constexpr int n = 200;
    for (int i = 0; i < n; ++i) {
        estimate_period(filt, accel, gyro, mag);
    }

    const auto e = filt.orientation().to_euler<EulerOrder::ZYX>();
    fmt::print("After {} ticks at rest (mock IMU):\n", n);
    fmt::print(
        "  yaw={:.4f} rad, pitch={:.4f} rad, roll={:.4f} rad\n",
        e.yaw(),
        e.pitch(),
        e.roll()
    );
    fmt::print(
        "  |yaw|,|pitch|,|roll| expected near 0 (identity attitude)\n"
    );
    return 0;
}
