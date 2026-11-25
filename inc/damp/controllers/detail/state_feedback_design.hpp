// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
// IWYU pragma: private: a public header that exports it provides its names

/**
 * @file state_feedback_design.hpp
 * @brief The design behind design::state_feedback: detail::designed, the gain on the plant or on its
 *        model with input delay states, then the filter and the target the options ask for
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/law_detail.hpp"
#include "damp/controllers/detail/lq_cost.hpp"
#include "damp/controllers/detail/placed_design.hpp"
#include "damp/controllers/state_feedback_types.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/input_limits.hpp"
#include "damp/design/options.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/status.hpp"
#include "damp/design/steady_state_target.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/decomposition.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace design {
namespace detail {

/**
 * @brief (∂u/∂ξ)⁺: the pseudo-inverse of −K on the integrators, the last NI columns of @p K. A
 *        clamp's excess Δu moves ξ by its least-squares share, ∂u/∂ξ·Δξ ≈ Δu; empty when the SVD fails.
 */
template<size_t NI, size_t NU, size_t NK, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NI, NU, T>> integrator_unwind(const Matrix<NU, NK, T>& K) {
    Matrix<NU, NI, T> du_dxi{};
    for (size_t i = 0; i < NU; ++i) {
        for (size_t j = 0; j < NI; ++j) {
            du_dxi(i, j) = -K(i, NK - NI + j);
        }
    }
    return mat::pseudo_inverse(du_dxi);
}

/** @brief The gain on the law's state. */
template<typename Result, size_t NK, size_t NU, typename T>
[[nodiscard]] constexpr Result with_gain(Result out, const LQRResult<NK, NU, T>& lq) {
    out.K = lq.K;
    out.S = lq.S;
    out.e = lq.e;
    return set_status(out, lq.status);
}

/**
 * @brief The law on a measured state: the gain on [x; xi] for the cost of @p w, and the
 *        integrators @p integ (NI > 0), in the design's @p Result. A continuous plant with no Ts
 *        and no integrators gets the continuous law (CARE).
 */
template<typename Result, Domain D, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NI>
[[nodiscard]] constexpr Result law_design(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T Ts, const LQRWeights<NX, NU, T>& w, const Integrators<NI, NY, T>& integ) {
    const auto& Ci = integ.Ci;
    if (!mat::cholesky(w.R).has_value()) {
        return refused<Result>(DesignStatus::bad_weights); // R not positive definite, on every route
    }
    if constexpr ((D == Domain::plant) && NI == 0) {
        if (sys.is_continuous()) {
            const auto lq = continuous_lqr(sys.A, sys.B, w.Q, w.R, w.N);
            return with_gain(Result{}, lq);
        }
    }
    const auto m = sampled_model(sys, Ts, D);
    if (m.status != DesignStatus::ok) {
        return refused<Result>(m.status);
    }
    const T h = (D == Domain::plant) ? sys.Ts : Ts;
    if constexpr (NI == 0) {
        LQRCost<NX, NU, T> cd{w.Q, w.R, w.N};
        if constexpr (D == Domain::sampled) {
            const auto c = discretize_lqr_cost(sys.A, sys.B, w.Q, w.R, Ts, w.N);
            if (!c) {
                return refused<Result>(DesignStatus::bad_weights); // the plant was sampled: a non-finite weight
            }
            cd = *c;
        }
        return with_gain(Result{.Ts = h}, discrete_lqr(m.law.A, m.law.B, cd));
    } else {
        if (!mat::cholesky(integ.Q).has_value()) {
            return refused<Result>(DesignStatus::bad_weights); // an unweighted integrator: a mode at z = 1 the cost does not see
        }
        // Forward-Euler integrator xi += h·e (MATLAB®'s discrete lqi): xi ≈ ∫e dt, so its cost
        // means the same on either route. The runtime integrates e = r − Ci·(C·x + D·u), u the
        // input applied over the period, and so does this model on both routes.
        const auto              im = integral_model(m.law.A, m.law.B, m.law.C, m.law.D, Ci, h);
        LQRCost<NX + NI, NU, T> cd = stacked_weights(w, integ.Q);
        if constexpr (D == Domain::sampled) {
            const auto c = sampled_integral_cost(sys.A, sys.B, sys.C, sys.D, Ci, cd, Ts);
            if (!c) {
                return refused<Result>(DesignStatus::bad_weights); // the plant was sampled: a non-finite weight
            }
            cd = *c;
        }
        Result out{.Ts = h};
        out.integral = IntegralDesign<NI, NX, NU, NY, T>{Ci, Ci * sys.C, Ci * sys.D, h};
        return with_gain(out, discrete_lqr(im.A, im.B, cd));
    }
}

/**
 * @brief The law on a measured state of a continuous plant whose input arrives @p delay late:
 *        the gain on [x; v(k−1); …; v(k−N); ξ] for the continuous cost of @p w, exact (delayed_cost),
 *        and the integrators @p integ (NI > 0), in the design's @p Result.
 */
template<typename Result, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t N, typename S, size_t NI>
[[nodiscard]] constexpr Result delayed_law_design(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        Ts,
    const InputDelay<N, S>&                  delay,
    const LQRWeights<NX, NU, T>&             w,
    const Integrators<NI, NY, T>&            integ
) {
    const auto& Ci = integ.Ci;
    if (!mat::cholesky(w.R).has_value()) {
        return refused<Result>(DesignStatus::bad_weights); // R not positive definite
    }
    if constexpr (NI > 0) {
        if (!mat::cholesky(integ.Q).has_value()) {
            return refused<Result>(DesignStatus::bad_weights); // an unweighted integrator
        }
    }
    const auto m = sampled_model(sys, Ts, Domain::sampled, delay); // a continuous plant, a usable Ts, finite entries, a delay discretize takes
    if (m.status != DesignStatus::ok) {
        return refused<Result>(m.status);
    }
    const auto& model = m.law;
    const auto  inputs = ::damp::detail::channel_delay_splits<NU>(delay.seconds, Ts, N); // as discretize accepted them
    const auto  cost = inputs ? delayed_cost(sys, model, *inputs, Ci, w, integ.Q, Ts) : damp::nullopt;
    if (!cost) {
        return refused<Result>(DesignStatus::bad_weights); // a non-finite weight on the sampled plant
    }
    Result out{.Ts = Ts};
    if constexpr (NI == 0) {
        return with_gain(out, discrete_lqr(model.A, model.B, *cost));
    } else {
        const auto im = integral_model(model.A, model.B, model.C, model.D, Ci, Ts);
        out.integral = IntegralDesign<NI, NX + (N * NU), NU, NY, T>{Ci, Ci * model.C, Ci * model.D, Ts};
        return with_gain(out, discrete_lqr(im.A, im.B, *cost));
    }
}

/**
 * @brief The gain of a design in its @p Result: the LQ gain of @p w (@p gain LQGain), on the plant
 *        (law_design) or on the model with the input delay states of @p delay
 *        (delayed_law_design), or the gain placed at @p gain's poles on the same model
 *        (placed_law_design); with the integrators @p integ (NI > 0).
 */
template<typename Result, Domain D, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NI, typename Delay, typename Gain>
[[nodiscard]] constexpr Result gain_design(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        Ts,
    const LQRWeights<NX, NU, T>&             w,
    const Integrators<NI, NY, T>&            integ,
    const Delay&                             delay,
    const Gain&                              gain
) {
    using In = damp::remove_cvref_t<decltype(input_delay_of(delay))>;
    if constexpr (!damp::is_same_v<Gain, LQGain>) {
        static_cast<void>(w); // placed at its poles: the weights are the LQ gain's
        return placed_law_design<Result, D, NI>(sys, Ts, gain, integ.Ci, input_delay_of(delay));
    } else if constexpr (damp::is_same_v<In, NoDelay>) {
        static_cast<void>(gain);
        static_cast<void>(delay);
        return law_design<Result, D>(sys, Ts, w, integ); // a late reading leaves the gain the plant's: only the filter carries it
    } else {
        static_cast<void>(gain);
        return delayed_law_design<Result>(sys, Ts, input_delay_of(delay), w, integ); // a delay only on the sampled form: check_delays
    }
}

/**
 * @brief The target of an LQ law with a measured disturbance, on the model of @p sys with @p delay:
 *        under [d̂; d] holding at r the outputs its spec or TrackedOutputs lists (@p Held), else
 *        the steady state under d that holds at 0 the integrated outputs @p Ci (an LQI's, NI > 0,
 *        whose integrators then hold r) or every output (a regulator), with no map of r. Either
 *        holds at most NU outputs (check_held_outputs): a regulator on a plant with more outputs
 *        than inputs does not compile under a measured d, and lists the ones it holds at r = 0.
 *
 * The law is u = u_t − K·(x − x_t) about it: for a d held constant, as the MPC holds it over its
 * horizon, that is the LQ law in the deviation from the steady state under d, so no term beyond
 * the target is optimal.
 */
template<bool Held, Domain D, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Spec, size_t NI, size_t NM, typename Delay>
[[nodiscard]] constexpr auto measured_target(
    const StateSpace<NX, NU, NY, T, NW, NV>&                               sys,
    T                                                                      Ts,
    const Spec&                                                            spec,
    const Matrix<NI, NY, T>&                                               Ci,
    const damp::optional<Matrix<NX + delay_states<Delay, NU, NY>, NM, T>>& Bm,
    const Delay&                                                           delay
) {
    constexpr size_t NXE = NX + delay_states<Delay, NU, NY>;
    if constexpr (Held) {
        static_cast<void>(Ci);
        using Result = SteadyStateTargetResult<NXE, NU, disturbance_states<Spec, NU, NY, NW> + NM, target_outputs<Spec, NY>, T>;
        if (!Bm) {
            return refused<Result>(DesignStatus::bad_weights); // a non-finite B_d
        }
        return target_design<D>(sys, Ts, spec, delay, *Bm);
    } else {
        static_cast<void>(spec);
        if (!Bm) {
            return refused<SteadyStateTargetResult<NXE, NU, NM, 0, T>>(DesignStatus::bad_weights); // a non-finite B_d
        }
        if constexpr (NI > 0) {
            return held_at_zero<D>(sys, Ts, Ci, delay, *Bm);
        } else {
            return held_at_zero<D>(sys, Ts, Matrix<NY, NY, T>::identity(), delay, *Bm);
        }
    }
}

/** @brief Ci of the integrators on the listed outputs, one row each; nullopt for an index out of range or listed twice. */
template<size_t NI, size_t NY, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NI, NY, T>> output_selector(const damp::array<size_t, NI>& outputs) {
    const auto rows = selected_outputs<NY>(outputs);
    if (!rows) {
        return damp::nullopt;
    }
    Matrix<NI, NY, T> Ci{};
    for (size_t i = 0; i < NI; ++i) {
        Ci(i, (*rows)[i]) = T{1};
    }
    return Ci;
}

/**
 * @brief @p law with the @p limits it runs within (the command's, its move's, the integrators' and
 *        the demand limit), the integrators' @p tracking_time and the fade time of @p handover (a
 *        Handover; Unset for none). A refused design ships none.
 * @return Refused with DesignStatus::invalid_parameter for a limit that is not a box (an end NaN
 *         or ±Inf, or a lower end above its upper one), a negative tracking time or a fade time
 *         that is not finite and positive; not_converged when back-calculation is asked for and
 *         the gain has no SVD; else @p law, a failed one as it is
 */
template<typename Result, size_t NX, size_t NU, typename T, size_t NI, size_t NDU, typename Fade>
[[nodiscard]] constexpr Result with_limits(Result law, const StateFeedbackLimits<NX, NU, T, NI, 0, NDU>& limits, T tracking_time, const Fade& handover) {
    if (!law.success) {
        return law;
    }
    bool usable = limits.input.valid() && limits.input_rate.valid() && limits.integral.valid() && limits.demand.template valid<NX>() && damp::finite_non_negative(tracking_time);
    if constexpr (!is_unset<Fade>) {
        law.fade_time = static_cast<T>(handover.fade_time);
        usable = usable && damp::finite_positive(handover.fade_time) && damp::finite_positive(law.fade_time); // and still > 0 in the plant's scalar
    } else {
        static_cast<void>(handover);
    }
    if (!usable) {
        return refused<Result>(DesignStatus::invalid_parameter);
    }
    law.limits.input = limits.input;
    law.limits.input_rate = limits.input_rate;
    law.limits.integral = limits.integral;
    law.limits.demand = limits.demand;
    if constexpr (NI > 0) {
        law.tracking_time = tracking_time;
        if (const auto unwind = integrator_unwind<NI>(law.K)) {
            law.integral.unwind = *unwind; // formed for any tracking time: the runtime's may be set later
        } else if (tracking_time > T{0}) {
            return refused<Result>(DesignStatus::not_converged); // back-calculation asked for, and no SVD of the gain
        }
    }
    return law;
}

/**
 * @brief The state-feedback design: the law (the LQ gain of @p w, or placed at @p gain's poles,
 *        LQGain or PlacedPoles) with the integrators @p integ (NI > 0) and its @p limits
 *        (with_limits), then the filter of @p filter (a KalmanNoise, the PlacedPoles of
 *        ObserverPoles, or Unset for none) and the target its options ask for: under the
 *        disturbance model @p offset_free, else on the outputs @p track lists (TrackedOutputs),
 *        under the measured disturbance @p measured (NM = 0: none); on the plant in its own domain
 *        (@p D plant) or sampled at @p Ts, with the delay states of @p delay. NoDisturbance is the
 *        absent offset-free spec and track. The result's type is worked out first and filled once.
 * @return The result; status the first failure: the gain's, the limits', the filter's, the target's
 */
template<
    Domain D, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NI, size_t NDU, typename Filter, typename OffsetFree, typename Track, size_t NM,
    typename Delay, typename Gain, typename Fade>
[[nodiscard]] constexpr auto designed(
    const StateSpace<NX, NU, NY, T, NW, NV>&          sys,
    T                                                 Ts,
    const LQRWeights<NX, NU, T>&                      w,
    const Integrators<NI, NY, T>&                     integ,
    const StateFeedbackLimits<NX, NU, T, NI, 0, NDU>& limits,
    const Filter&                                     filter,
    const OffsetFree&                                 offset_free,
    const Track&                                      track,
    const MeasuredDisturbance<NX, NM, T>&             measured,
    const Delay&                                      delay,
    const Gain&                                       gain,
    const Fade&                                       handover
) {
    constexpr bool kEstimated = !is_unset<Filter>;
    constexpr bool kHeld = DisturbanceSpec<OffsetFree> || OptionOf<Track, kind::TrackedOutputs>; // r is the target's
    constexpr bool kTarget = kHeld || (NM > 0);                                                  // a measured d moves every law's target
    const auto     spec = target_spec_of(offset_free, track);
    const auto     Bm = measured_channel<D>(sys, Ts, measured, delay); // d on the model, none on the input's delay states
    const auto     estimator = [&] {
        if constexpr (kEstimated) {
            return with_measured_inputs(filter_design<D>(sys, Ts, filter, spec, delay), Bm); // d a known input after u
        } else {
            static_cast<void>(filter);
            return NoEstimator{};
        }
    }();
    constexpr size_t NZ = kHeld ? target_outputs<damp::remove_cvref_t<decltype(spec)>, NY> : 0;
    constexpr size_t NYL = (kEstimated || NI > 0) ? NY : NX; // what the law reads: y with a filter or integrators, else x
    using Result = StateFeedbackResult<NX, NU, T, NI, NYL, damp::remove_cvref_t<decltype(estimator)>, NZ, input_delay_states<Delay>, NM, NDU, !is_unset<Fade>>;

    Result       out = with_limits(gain_design<Result, D>(sys, Ts, w, integ, delay, gain), limits, integ.tracking_time, handover);
    DesignStatus status = out.status;
    out.kalman = estimator;
    if constexpr (kEstimated) {
        status = (status == DesignStatus::ok) ? estimator.status : status;
    }
    if constexpr (kTarget && NM == 0) {
        out.target = target_design<D>(sys, Ts, spec, delay);
    } else if constexpr (kTarget) {
        out.target = measured_target<kHeld, D>(sys, Ts, spec, integ.Ci, Bm, delay);
    } else {
        static_cast<void>(Bm);
    }
    if constexpr (kTarget) {
        status = (status == DesignStatus::ok) ? out.target.status : status;
    }
    return set_status(out, status);
}

} // namespace detail
} // namespace design
} // namespace damp
