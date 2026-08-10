// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file fo_plant_inertia_sil.cpp
 * @brief Online J,b from FirstOrderPlantEstimator — host identification demo
 *
 * Calls fo_plant_inertia_estimator.hpp for config, discrete plant, tick, and map.
 */

#include "damp/estimation/parameter_estimation.hpp"
#include "damp/math/math.hpp"
#include "fmt/core.h"
#include "fo_plant_inertia_estimator.hpp"

using namespace damp;
using namespace damp::examples_fo_plant_inertia;

int main() {
    const double a = discrete_a();
    const double b_disc = discrete_b();

    FirstOrderPlantEstimator est{kCfg};

    double omega = 0.0;
    double torque_prev = 0.0;
    for (std::size_t k = 0; k < 1200; ++k) {
        const double torque = 0.5 * excitation(k, 10);
        omega = estimate_period(est, omega, torque_prev, a, b_disc, torque);
        torque_prev = torque;
    }

    const auto mech = map_mechanical(est);
    if (!mech.valid) {
        fmt::print("FO plant estimate not valid yet (need excitation).\n");
        return 1;
    }

    fmt::print(
        "First-order plant: K={:.4f}, tau={:.4f} s (conf={:.3f})\n",
        est.gain(),
        est.time_constant(),
        est.confidence()
    );
    fmt::print(
        "Mechanical map:    J={:.5f} kg m^2 (true {:.5f}), b={:.5f} N m s/rad (true {:.5f})\n",
        mech.J,
        J_true,
        mech.b,
        b_true
    );
    return 0;
}
