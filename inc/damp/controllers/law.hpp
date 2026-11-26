// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file law.hpp
 * @brief Law<kRes, T>: the runtime a design result deploys as, named by the result
 *
 * Each runtime's header specializes it for the results it runs, and exports this header; it exports
 * the parts those specializations share (detail/law_detail.hpp).
 */

#include "damp/controllers/detail/law_detail.hpp" // IWYU pragma: export

namespace damp {

/**
 * @brief The runtime a design result deploys as, named by the result: `Law<kRes, float> law;`
 *
 * A struct member's type cannot be deduced, so it is named by its design instead of its sizes.
 * kRes is a named result with static storage (namespace scope, or a `static` member or local).
 *
 * On a `constexpr` kRes the design's matrices are constants of the type: they take no RAM, and a
 * coefficient that is exactly 0 is no code, which removes the sparse work of a delay state, a
 * plant with an integrator or a measurement of a single state. The ticks and the numbers are the
 * runtime's own. Nothing can replace them, so a member a redesign assigns is
 * `Law<kRes, float>::Runtime law{kRes.as<float>()};`. The coefficients move from RAM to flash, as
 * constants in the code: where flash has wait states and no cache, keep the tick in RAM or tightly
 * coupled memory.
 *
 * On a kRes initialized at start-up, from values measured or identified at init, the matrices are
 * members, and assigning a runtime of the same shape replaces them.
 *
 * @code
 * inline constexpr auto kRes = design::state_feedback(plant, {.Q = Q, .R = R}, design::KalmanNoise{.process = Qn, .measurement = Rn});
 *
 * struct Axis {
 *     Law<kRes, float> law; // the design's law, in float
 *     float            offset{0.0f};
 * };
 * @endcode
 *
 * The runtime's header defines it for the results that runtime is built from: StateFeedback, the
 * LQ runtime (controllers/state_feedback_runtime.hpp, which lqr.hpp, lqi.hpp, lqg.hpp, lqgi.hpp and
 * oflqg.hpp include), and the observers
 * (estimation/fixed_gain_observer.hpp: FixedGainObserver, ReducedOrderObserver), design::eso (estimation/eso.hpp),
 * design::adrc (controllers/adrc.hpp) and the model inverses (design/model_inverse.hpp:
 * StateFeedforward, OutputInverse), which all run on constants when kRes is `constexpr`, MPC
 * (controllers/mpc.hpp), whose QP data stays in members and whose filter, when it has one, is on
 * constants, and PIDController (controllers/pid.hpp: design::pid), whose gains
 * stay members, settable at run time, in the narrowest mode a `constexpr` design's gains need.
 * SMCController and STSMCController (controllers/smc.hpp, controllers/stsmc.hpp: design::smc,
 * design::stsmc) run on their gains as constants of the type when kRes is `constexpr`, and
 * ExtremumSeekingController and SmithPredictor (controllers/esc.hpp, controllers/smith_predictor.hpp:
 * design::esc, design::smith_predictor) on members either way.
 * design::kalman's Law also names its recursive filter, `Law<kRes, float>::Recursive`: the KalmanFilter
 * of the same model and noise, on the same constants when kRes is `constexpr`.
 *
 * @tparam Res The design result
 * @tparam T   The runtime scalar
 */
template<const auto& Res, typename T>
struct Law {
    static_assert(sizeof(T) == 0, "Law<res, T>: res is not a design result a law is built from; include its runtime's header (controllers/state_feedback_runtime.hpp or lqr.hpp, lqi.hpp, lqg.hpp, lqgi.hpp, oflqg.hpp, mpc.hpp, adrc.hpp, pid.hpp, smc.hpp, stsmc.hpp, esc.hpp, smith_predictor.hpp, estimation/fixed_gain_observer.hpp, estimation/eso.hpp, design/model_inverse.hpp)");
};

} // namespace damp
