// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lqg.hpp
 * @brief The linear-quadratic Gaussian law: LQG, the StateFeedback of design::state_feedback with a KalmanNoise, LQGResult and lqg_from_parts
 *
 * design::state_feedback(plant, [Ts,] {.Q, .R, .N}, KalmanNoise{.process, .measurement}[, TrackedOutputs{{i, …}}][, InputDelay<N>{θ}])
 * (state_feedback.hpp) → LQGResult → LQG. With an InputDelay the filter carries the inputs still on their way.
 * step(y) is predict, update(y), control. The low-latency tick is update(y), control(),
 * back_calculate, predict. With tracked outputs: control(r), step(r, y). With
 * feedforward from a planned trajectory: control([r,] x_ref, u_ref), step([r,] x_ref, u_ref, y).
 * While another source drives the plant, track(u, …) takes step's place: the filter runs on u.
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/lq_law.hpp" // IWYU pragma: export
#include "damp/controllers/detail/state_feedback_core.hpp"
#include "damp/controllers/law.hpp"                    // IWYU pragma: export
#include "damp/controllers/state_feedback_runtime.hpp" // IWYU pragma: export
#include "damp/controllers/state_feedback_types.hpp"   // IWYU pragma: export
#include "damp/design/riccati.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/kalman.hpp" // IWYU pragma: export

namespace damp {

namespace design {

/** @brief an LQG design's result: the gain and the steady-state Kalman filter. */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = NX, size_t NV = NY>
using LQGResult = StateFeedbackResult<NX, NU, T, 0, NY, KalmanResult<NX, NU, NY, T, NW, NV, 0>, 0>;

/**
 * @brief Assemble an LQG design from separately computed Kalman and LQR designs
 *
 * The regulator must be discrete: a continuous one (Ts = 0) is refused with
 * DesignStatus::invalid_sample_time. A design on the plant (design::state_feedback) carries its sample time,
 * which must match kest.sys.Ts; a gain on bare matrices (discrete_lqr, matlab::dlqr) is discrete
 * with no sample time (Ts = 0), so the law takes the filter's. A failure already
 * reported by either half is kept (the regulator's first).
 *
 * @note Compare with MATLAB®'s lqgreg(kest, k) — exposed as the lqgreg() alias in matlab.hpp.
 *
 * @param kest  Steady-state Kalman design (design::kalman, no options)
 * @param reg   The regulator on the measured state: the LQRResult of design::state_feedback, discrete_lqr or
 *              matlab::dlqr
 * @return The LQG design: deploy with `LQG law{res}`
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Reg>
    requires damp::is_same_v<Reg, StateFeedbackResult<NX, NU, T>> || damp::is_same_v<Reg, LQRResult<NX, NU, T>>
[[nodiscard]] constexpr StateFeedbackResult<NX, NU, T, 0, NY, KalmanResult<NX, NU, NY, T, NW, NV>> lqg_from_parts(
    const KalmanResult<NX, NU, NY, T, NW, NV>& kest,
    const Reg&                                 reg
) {

    auto out = StateFeedbackResult<NX, NU, T, 0, NY, KalmanResult<NX, NU, NY, T, NW, NV>>{reg.K, reg.S, reg.e, false, DesignStatus::not_converged, kest.sys.Ts, kest};
    auto status = reg.success ? kest.status : reg.status;

    if (status == DesignStatus::ok && (!::damp::detail::discrete_gain(reg) || (reg.Ts != T{0} && reg.Ts != kest.sys.Ts))) {
        status = DesignStatus::invalid_sample_time; // a continuous regulator, or the halves disagree on Ts
    }
    if (status == DesignStatus::ok && (!reg.limits.input.valid() || !reg.limits.demand.template valid<NX>())) {
        status = DesignStatus::invalid_parameter; // a regulator built by hand with a limit that is not a box
    }
    if (status == DesignStatus::ok) {
        out.limits.input = reg.limits.input;
        out.limits.demand = reg.limits.demand;
    }
    return detail::set_status(out, status);
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Linear-quadratic Gaussian law, LQR on the Kalman filter's estimate: the StateFeedback of a design with a filter and no integrators or disturbance model
 *
 * Runs design::state_feedback with a KalmanNoise. The estimator is a FixedGainObserver with design::kalman's fixed gain
 * (`estimator`). The tick is step(y), or update(y), control(), back_calculate, predict()
 * (predict moves the filter): the regulator. A design with
 * TrackedOutputs holds those outputs at r through its steady-state target: control(r),
 * step(r, y). Feedforward from a planned trajectory (x_ref, u_ref) adds to it:
 * control([r,] x_ref, u_ref), step([r,] x_ref, u_ref, y). The command is clamped into the
 * design's InputLimits, and its move into an InputRateLimits (unbounded by default; limits() reads
 * them and set_limits() moves them at run time), and predict() moves the observer with it,
 * or with the input back_calculate(u_sat) recorded for a stage that clamps after the law.
 *
 * A design with an InputDelay (NRU > 0) has its filter on [x; v(k−1); …; v(k−NRU)]: predict
 * shifts the applied input into the delay states. x_ref is the plant's; the delay states' reference is
 * u_ref.
 *
 * @note Compare with MATLAB®'s lqg.
 *
 * A design with a MeasuredDisturbance (NM > 0) takes d each tick, after the reference and its
 * feedforward: the target under it, as on LQR (the regulator holding every output at 0), and the
 * filter's known input after u, so the estimate does not read the load as a model error:
 * control([r,] d), step([r,] d, y), control([r,] x_ref, u_ref, d), step([r,] x_ref, u_ref, d, y).
 *
 * While another source drives the plant (an operator, another law), track(u, …) takes the place of
 * each step(…), its arguments step's with the input applied first: track(u, y), track(u, r, y),
 * track(u, x_ref, u_ref, y), track(u, d, y), …. The filter predicts with u (and the measured d), so
 * its estimate is current at the handover. The law has nothing to preload, so its first step is its
 * own command, unless the design took a design::Handover: then that command is u, and the mismatch
 * with the law's own decays by exp(−Ts/fade_time) each tick, within its input limits.
 *
 * The alias takes a filter (Est a FixedGainObserver) that estimates no disturbance, and no
 * integrators, so `LQG law{res}` refuses an LQR's, an LQI's, an LQGI's or an OFLQG's design.
 * StateFeedback lists the members and the tick.
 *
 * @tparam NX    Plant states
 * @tparam NU    Inputs
 * @tparam NY    Measured outputs
 * @tparam T     Scalar
 * @tparam Est   The filter: a FixedGainObserver on the law's model, with no disturbance states
 * @tparam NZ    Tracked outputs: the size of r (0: none, the regulator)
 * @tparam Store Where the matrices live (detail::HeldGains: members)
 * @tparam NRU   Input delay states per input of the filter's model (0: none)
 * @tparam NM    Measured-disturbance inputs (0: none)
 * @tparam NDU   Inputs with a move limit (an InputRateLimits: NU; 0: none)
 * @tparam HO    Whether the law fades the mismatch after a track (a Handover; false: none)
 */
template<size_t NX, size_t NU, size_t NY, typename T, typename Est, size_t NZ = 0, typename Store = detail::HeldGains, size_t NRU = 0, size_t NM = 0, size_t NDU = 0, bool HO = false>
    requires(detail::runs_filter<Est> && !detail::estimates_disturbance<Est>)
using LQG = StateFeedback<NX, NU, T, 0, NY, Est, NZ, Store, NRU, NM, NDU, HO>;

} // namespace damp
