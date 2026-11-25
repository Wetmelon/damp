// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file state_limits.hpp
 * @brief StateLimits, the limits on the states a QP design lists: a trailing option of
 *        design::mpc and design::mhe, soft by default in both
 */

#include "damp/backend.hpp"
#include "damp/design/options.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/core.hpp"

namespace damp::design {

/**
 * @struct StateLimits
 * @brief Limits on the states a QP design lists: among design::mpc's and design::mhe's options
 *
 * lower and upper bound the NXC states of the design's model listed in states (a delayed model's
 * delay states come after the plant's states). `StateLimits{.states = {2}, .lower = {-0.26},
 * .upper = {0.26}}` bounds state 2; states is an array, so the limits take NXC from its list, and
 * an end left out is unbounded. The QP carries rows for the listed states alone, so a design
 * without a StateLimits has none.
 *
 * | Design | What the limits are |
 * | --- | --- |
 * | design::mpc (MPC, OFMPC) | rows on the predicted states over the whole prediction horizon: a lean, a speed or a temperature the law must keep within without tracking it |
 * | design::mhe (MHE) | rows on every state of the estimation window: a level, a pressure or a state of charge the estimate must respect |
 *
 * Soft by default, the same way in both (MATLAB®'s equal concern for relaxation): one slack
 * @f$ \varepsilon \ge 0 @f$, shared by every soft row of the QP, relaxes each row to
 * @f$ x_i / s_i \le \bar x_i / s_i + V_i\,\varepsilon @f$ at the cost
 * @f$ \rho_\varepsilon \varepsilon^2 @f$, @f$ s_i @f$ the state's scale (its range, upper − lower,
 * where both ends are finite, else 1) and @f$ V_i @f$ the row's ECR coefficient: lower_ecr(i) for
 * the lower(i) rows and upper_ecr(i) for the upper(i) rows. A coefficient of 0 makes that side
 * hard, and a larger one relaxes it more readily relative to the others. @f$ \rho_\varepsilon @f$
 * is MPCWeights::ecr_weight in the MPC, and 10⁵ in the MHE, that weight's default, against a
 * window cost whose residuals are counted in standard deviations of their noise.
 *
 * A listed state left unbounded has its rows all the same, for set_limits to move at run time
 * (MPCLimits::state, MHELimits::state, in the order of states). Either design refuses an index
 * past the model's states or listed twice, an end NaN or ±Inf or a lower end above its upper one
 * (invalid_parameter), and a negative or non-finite ECR (bad_weights). An MPC with an
 * InvariantSet, whose rows do not cover the states, does not compile with one.
 */
template<size_t NXC, typename T = double>
struct StateLimits {
    using option_kind = kind::StateLimits; ///< Its kind of design option
    static constexpr size_t listed = NXC;  ///< The states it bounds
    static_assert(NXC > 0, "StateLimits lists the states it bounds: .states = {i, ...}, at least one.");

    size_t         states[NXC]{};                                                        ///< Indices of the model's states bounded, each once (an array, so bare braces give NXC)
    ColVec<NXC, T> lower = Matrix<NXC, 1, T>::constant(-damp::numeric_limits<T>::max()); ///< Their lowest values (default unbounded; soft through lower_ecr)
    ColVec<NXC, T> upper = Matrix<NXC, 1, T>::constant(damp::numeric_limits<T>::max());  ///< Their highest values (default unbounded; soft through upper_ecr)
    ColVec<NXC, T> lower_ecr = Matrix<NXC, 1, T>::constant(T{1});                        ///< Softness of the lower rows, ≥ 0 (0 = hard)
    ColVec<NXC, T> upper_ecr = Matrix<NXC, 1, T>::constant(T{1});                        ///< Softness of the upper rows, ≥ 0 (0 = hard)
};

} // namespace damp::design
