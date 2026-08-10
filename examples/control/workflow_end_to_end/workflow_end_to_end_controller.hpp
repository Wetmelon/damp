// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file workflow_end_to_end_controller.hpp
 * @brief Linearize → LQGI + PR workflow — nameplate plant + rates
 *
 * Heavy design runs on host (sil); float runtime tick is the deploy surface.
 */

#pragma once

#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"

namespace damp::examples_workflow_e2e {

inline constexpr double kTs = 0.001; // 1 kHz

/** Nonlinear plant used for linearization and SIL. */
[[nodiscard]] inline ColVec<2, double> plant_nonlinear(
    double /*t*/, const ColVec<2, double>& x, const ColVec<1, double>& u
) {
    const double x1 = x(0, 0);
    const double x2 = x(1, 0);
    const double uu = u(0, 0);
    return ColVec<2, double>{x2, (-0.8 * x2) - (2.0 * damp::sin(x1)) + (1.5 * uu)};
}

[[nodiscard]] inline ColVec<1, double> plant_output(const ColVec<2, double>& x) {
    return ColVec<1, double>{x(0, 0)};
}

inline constexpr ColVec<2, double> kXop{0.0, 0.0};
inline constexpr ColVec<1, double> kUop{0.0};

} // namespace damp::examples_workflow_e2e
