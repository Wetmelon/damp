// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file encoder_velocity_estimator.hpp
 * @brief Encoder velocity nameplate + critically-damped PLL observer tick
 *
 * Folder layout:
 *
 *   encoder_velocity_estimator.hpp  — this file (dt, PLL observer)
 *   encoder_velocity_sketch.cpp     — thin smoke
 *   encoder_velocity_sil.cpp        — host comparison + plots
 *   encoder_velocity_derivation.md  — plant / estimator notes
 *
 * SIL also runs LPF, raw FD, and Levant differentiator for teaching comparison.
 */

#pragma once

#include "damp/filters/differentiator.hpp"

namespace damp::examples_encoder_velocity {

inline constexpr double dt = 1.0 / 8000.0; // 8 kHz current-loop rate [s]

/**
 * @brief 2nd-order critically-damped PLL position/velocity tracking observer
 *
 * Gains: $`k_p = 2\,\mathrm{bw}`$, $`k_i = 0.25\,k_p^2`$ (repeated pole at $`-\mathrm{bw}`$).
 */
struct PllObserver {
    double kp;
    double ki;
    double pos{0.0};
    double vel{0.0};

    explicit PllObserver(double bw) : kp(2.0 * bw), ki(0.25 * (2.0 * bw) * (2.0 * bw)) {}

    double update(double meas) {
        pos += dt * vel;
        const double e = meas - pos;
        pos += dt * kp * e;
        vel += dt * ki * e;
        return vel;
    }
};

/**
 * @brief One period: PLL velocity from quantized position [turns]
 */
inline double estimate_period_pll(PllObserver& pll, double theta_meas) {
    return pll.update(theta_meas);
}

/**
 * @brief One period: Levant robust exact differentiator velocity
 */
template<typename T>
T estimate_period_levant(RobustExactDifferentiator<T>& red, T theta_meas) {
    return red.update(theta_meas);
}

} // namespace damp::examples_encoder_velocity
