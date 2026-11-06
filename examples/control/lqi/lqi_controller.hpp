// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file lqi_controller.hpp
 * @brief Position LQI with a clamp — flashable design + tick
 *
 * Canonical path: TF → to_state_space → discretize → design::lqi → LQI → saturate →
 * back_calculate. SIL calls this same tick.
 *
 *   lqi_controller.hpp  — this file (nameplate, design::, control_period)
 *   lqi_sketch.cpp      — thin setup/loop
 *   lqi_sil.cpp         — host plant, clamp-aware vs naive
 *   lqi_derivation.md
 */

#pragma once

#include "damp/controllers/lqi.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/transfer_function.hpp"

namespace damp::examples_lqi {

inline constexpr double Ts = 0.01;    // 100 Hz
inline constexpr double u_min = -1.2; // just above k·r = 1 at the setpoint
inline constexpr double u_max = 1.2;

namespace {
constexpr double m = 1.0; // kg
constexpr double k = 1.0; // N/m
constexpr double c = 0.5; // N·s/m

// G(s) = 1 / (m s² + c s + k). Companion form is [position, speed] when m = 1.
constexpr auto G = TransferFunction<1, 3>{.num = {1.0}, .den = {m, c, k}};
constexpr auto sys = *G.to_state_space();

} // namespace

inline constexpr auto sysd = discretize(sys, Ts);
static_assert(sysd);

// Design an LQI "law" on the discretized plant, and check for success
inline constexpr auto lqi_result = design::lqi(
    *sysd,
    {
        .Q = Matrix<2, 2>{
            {10.0, 0.0},
            {0.0, 1.0},
        },
        .R = 1.0,
        .Q_integral = 80.0,
    }
);
static_assert(lqi_result.success);

/**
 * @brief One control period: update → control → clamp → back_calculate → predict.
 *
 * Pass the voltage (or force) that actually left the driver, not the unsaturated request.
 */
template<typename T>
[[nodiscard]] T control_period(LQI<2, 1, 1, T, 1>& law, const ColVec<2, T>& x, T r) {
    law.update(x);
    const T u = law.control(ColVec<1, T>{{r}})[0];
    const T u_sat = damp::clamp(u, static_cast<T>(u_min), static_cast<T>(u_max));
    law.back_calculate(ColVec<1, T>{{u_sat}});
    law.predict();
    return u_sat;
}

} // namespace damp::examples_lqi
