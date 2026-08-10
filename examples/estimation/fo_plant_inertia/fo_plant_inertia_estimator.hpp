// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file fo_plant_inertia_estimator.hpp
 * @brief Online J,b from FirstOrderPlantEstimator — nameplate + tick + map
 *
 * Mechanical axis $`J\dot\omega + b\,\omega = \tau`$ is a first-order plant
 * with $`K = 1/b`$ and $`\tau_m = J/b`$. Use the shared grey-box RLS estimator,
 * then map physical units in application code:
 *
 * $`b = 1/K`$, $`J = \tau_m / K`$ ($`K > 0`$).
 *
 * Folder layout:
 *
 *   fo_plant_inertia_estimator.hpp  — this file
 *   fo_plant_inertia_sketch.cpp     — thin smoke
 *   fo_plant_inertia_sil.cpp        — host identification demo
 *   fo_plant_inertia_derivation.md  — plant / estimator notes
 */

#pragma once

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/estimation/parameter_estimation.hpp"
#include "damp/math/math.hpp"

namespace damp::examples_fo_plant_inertia {

inline constexpr double Ts = 0.01;
inline constexpr double J_true = 0.01; // kg·m²
inline constexpr double b_true = 0.05; // N·m·s/rad

inline constexpr design::FirstOrderPlantEstimatorConfig<double> kCfg{
    .Ts = Ts,
    .forgetting = 0.999,
    .initial_covariance = 1e4,
};

// Short PRBS-like levels (hold long enough to see the mechanical pole).
inline constexpr damp::array<double, 8> kPattern{1, -1, 1, 1, -1, -1, 1, -1};

[[nodiscard]] constexpr double excitation(std::size_t k, std::size_t hold) {
    return kPattern[(k / hold) % kPattern.size()];
}

/**
 * @brief One identification period: plant step + estimator update
 *
 * @return updated omega (discrete plant state)
 */
template<typename Est>
[[nodiscard]] double estimate_period(
    Est&   est,
    double omega,
    double torque_prev,
    double a,
    double b_disc,
    double torque_cmd
) {
    omega = (a * omega) + (b_disc * torque_prev);
    (void)est.update(torque_cmd, omega);
    return omega;
}

/** Application map — not a dedicated library estimator. */
struct MechanicalParams {
    double J{};
    double b{};
    bool   valid{false};
};

template<typename Est>
[[nodiscard]] MechanicalParams map_mechanical(const Est& est) {
    MechanicalParams out{};
    if (!est.valid() || !(est.gain() > 0.0)) {
        return out;
    }
    out.b = 1.0 / est.gain();
    out.J = est.time_constant() / est.gain();
    out.valid = true;
    return out;
}

[[nodiscard]] constexpr double discrete_a() {
    return damp::exp(-b_true * Ts / J_true);
}

[[nodiscard]] constexpr double discrete_b() {
    return (1.0 / b_true) * (1.0 - discrete_a());
}

} // namespace damp::examples_fo_plant_inertia
