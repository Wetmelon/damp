// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file state_feedback.hpp
 * @brief design::state_feedback: the linear-quadratic state-feedback designs, the law picked by its options
 *
 * design::state_feedback(plant, [Ts,] weights, opts...) → StateFeedbackResult, which runs as
 * StateFeedback (`Law<kRes, T>`, or `StateFeedback law{res}`) of the kind its options name, each
 * kind's name an alias of it that takes only that kind's design (`LQI law{res}`). The weights are LQRWeights
 * {.Q, .R, .N}, the LQ gain, or Poles{…}, the gain placed there; a KalmanNoise among the options
 * gives the law a Kalman filter, ObserverPoles{…} a placed one, any of them with any weights:
 *
 * | KalmanNoise (or ObserverPoles) among the options | Other options | Kind |
 * | --- | --- | --- |
 * | — | none | LQR on the measured state (MATLAB®'s lqr, lqrd) |
 * | — | IntegralAction{.Q[, .outputs = {i, …}]}: its own cost, every output unless it lists them | LQI (lqi) |
 * | — | TrackedOutputs{{i, …}} | LQR holding the outputs at r through the steady-state target |
 * | KalmanNoise | none | LQG: LQR on a Kalman filter's estimate (lqg) |
 * | KalmanNoise | IntegralAction | LQGI (lqgtrack) |
 * | KalmanNoise | TrackedOutputs | LQG holding the outputs at r |
 * | KalmanNoise | a disturbance spec | OFLQG: a filter on [x; d] and the target under d̂ |
 * | either | MeasuredDisturbance{Bd}, beside any of these | the same law, its target under the d it reads each tick |
 *
 * With Ts, InputDelay<N>{θ} among the options designs any of these on the plant whose input
 * arrives θ late, and OutputDelay<M>{θ} the ones with a KalmanNoise on readings that arrive θ
 * late. The options trail in any order, at most one of each kind (design/options.hpp).
 * This header includes StateFeedback and the five names, so the law a design names deploys from it. The matlab::
 * aliases (lqr, dlqr, lqrd, lqi, lqg, lqgtrack) take MATLAB®'s arguments and call through.
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/state_feedback_design.hpp"
#include "damp/controllers/law.hpp"
#include "damp/controllers/lqg.hpp"                    // IWYU pragma: export
#include "damp/controllers/lqgi.hpp"                   // IWYU pragma: export
#include "damp/controllers/lqi.hpp"                    // IWYU pragma: export
#include "damp/controllers/lqr.hpp"                    // IWYU pragma: export
#include "damp/controllers/oflqg.hpp"                  // IWYU pragma: export
#include "damp/controllers/state_feedback_runtime.hpp" // IWYU pragma: export
#include "damp/controllers/state_feedback_types.hpp"
#include "damp/design/detail/config.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/input_limits.hpp" // IWYU pragma: export
#include "damp/design/options.hpp"
#include "damp/design/steady_state_target.hpp" // IWYU pragma: export
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/toolbox/bounds.hpp" // IWYU pragma: export

namespace damp {
namespace design {
namespace detail {

/**
 * @brief The integrators of an IntegralAction as the design reads them, on a plant of NY outputs in
 *        its scalar @p T: Ci of the outputs it lists (every one, in order, when it lists none), its
 *        cost and tracking time; nullopt for an output out of range or listed twice, or left out
 *        with fewer integrators than outputs.
 */
template<size_t NY, typename T, typename TQ, typename TI>
[[nodiscard]] constexpr auto integrators_of(const ::damp::design::IntegralAction<TQ, TI>& action) {
    constexpr size_t NI = ::damp::design::IntegralAction<TQ, TI>::NI;
    static_assert(
        NI > 0,
        "IntegralAction::Q is the cost on the integrators: IntegralAction{.Q = ...}, a scalar for one, an NI x NI Matrix for more. "
        "Given one inside a template? There clang deduces nothing from a .Q set to a literal or a constant (llvm/llvm-project#114870): "
        "name the type, IntegralAction<double>{.Q = ...}."
    );
    static_assert(square_value<TQ, NI>, "IntegralAction::Q is NI x NI, one integrator per integrated output: a Matrix, or a scalar for one.");
    static_assert(NI <= NY, "IntegralAction::Q has one integrator per integrated output: at most one per measured output.");
    damp::array<size_t, NI> outputs = action.outputs;
    if constexpr (NI == NY) {
        if (outputs[0] == kEveryOutput) {
            for (size_t i = 0; i < NI; ++i) {
                outputs[i] = i; // left out: every output, in order
            }
        }
    }
    const auto Ci = output_selector<NI, NY, T>(outputs);
    if (!Ci) {
        return damp::optional<Integrators<NI, NY, T>>{};
    }
    return damp::optional<Integrators<NI, NY, T>>{Integrators<NI, NY, T>{
        .Ci = *Ci,
        .Q = square_matrix<NI, T>(action.Q, Matrix<NI, NI, T>{}),
        .tracking_time = static_cast<T>(action.tracking_time),
    }};
}

/**
 * @brief The options of a state-feedback design: a filter spec (a KalmanNoise or ObserverPoles)
 *        makes it estimate the state. A weights spec among them is the weights given twice.
 */
using state_feedback_options = ::damp::detail::option_set<
    is_input_delay,
    is_output_delay,
    is_integral_option,
    is_tracked_option,
    is_disturbance,
    is_measured_disturbance,
    is_filter_spec,
    is_input_limits,
    is_input_rate_limits,
    is_demand_limit,
    is_handover,
    is_weight_spec>;

/** @brief The integrators of a state-feedback design's options: an IntegralAction's NI, else none. */
template<typename... Opts>
consteval size_t integrator_count() {
    using Action = ::damp::detail::option_t<is_integral_option, Unset, Opts...>;
    if constexpr (is_unset<Action>) {
        return 0;
    } else {
        return Action::NI;
    }
}

/** @brief The poles a spec lists: none for anything but Poles and ObserverPoles. */
template<typename S>
inline constexpr size_t listed_poles = 0;

/** @brief Poles: N. */
template<size_t N, typename U>
inline constexpr size_t listed_poles<Poles<N, U>> = N;

/** @brief ObserverPoles: N. */
template<size_t N, typename U>
inline constexpr size_t listed_poles<ObserverPoles<N, U>> = N;

/** @brief The gain a design's weights spec asks for: the LQ gain of LQRWeights. */
template<typename T>
[[nodiscard]] constexpr LQGain gain_spec_of(const LQGain& lq) {
    return lq;
}

/** @brief Poles: the gain placed there, in the plant's scalar. */
template<typename T, size_t N, typename U>
[[nodiscard]] constexpr PlacedPoles<N, T> gain_spec_of(const Poles<N, U>& poles) {
    return {pole_list<T>(poles.poles)};
}

/** @brief ObserverPoles: the filter's gain placed there, in the plant's scalar. */
template<typename T, size_t N, typename U>
[[nodiscard]] constexpr PlacedPoles<N, T> observer_spec_of(const ObserverPoles<N, U>& poles) {
    return {pole_list<T>(poles.poles)};
}

/** @brief The filter spec among @p opts: the KalmanNoise sized to the plant, or ObserverPoles (observer_spec_of). */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename... Opts>
[[nodiscard]] constexpr auto filter_spec_of(const Opts&... opts) {
    if constexpr (::damp::detail::option_count<is_observer_poles, Opts...> > 0) {
        return observer_spec_of<T>(::damp::detail::option<is_observer_poles>(Unset{}, opts...));
    } else {
        return noise_of<NX, NU, NY, T, NW, NV>(opts...);
    }
}

/**
 * @brief state_feedback_design_of once the options pass their checks: the options read by kind
 *        and the design on them (designed), its gain the LQ one or placed (@p gain).
 */
template<bool Sampled, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Gain, typename... Opts>
[[nodiscard]] constexpr auto state_feedback_counted(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        Ts,
    const LQRWeights<NX, NU, T>&             weights,
    const Gain&                              gain,
    const Opts&... opts
) {
    constexpr bool   kFilter = ::damp::detail::option_count<is_filter_spec, Opts...> > 0;
    constexpr size_t NI = integrator_count<Opts...>();
    constexpr size_t NDU = rate_limited_inputs<NU, Opts...>;
    const auto       filter = [&] {
        if constexpr (kFilter) {
            return filter_spec_of<NX, NU, NY, T, NW, NV>(opts...);
        } else {
            return Unset{};
        }
    }();
    const auto                                 offset_free = ::damp::detail::option<is_disturbance>(NoDisturbance{}, opts...);
    const auto                                 track = ::damp::detail::option<is_tracked_option>(NoDisturbance{}, opts...);
    const auto                                 measured = measured_of<NX, T>(opts...);
    const auto                                 delay = delays_of(opts...);
    const auto                                 handover = ::damp::detail::option<is_handover>(Unset{}, opts...);
    StateFeedbackLimits<NX, NU, T, NI, 0, NDU> limits{
        .input = input_limits_of<NU, T>(opts...),
        .demand = demand_terms_of<::damp::detail::DemandTerms<NX + NI, T>>(opts...), // its outer terms may be the integrators
    };
    if constexpr (NDU > 0) {
        limits.input_rate = input_rate_limits_of<NU, T>(opts...);
    }
    const auto design = [&](const Integrators<NI, NY, T>& integ) {
        if constexpr (kFilter) {
            return designed<!Sampled>(noise_model(sys), Ts, weights, integ, limits, filter, offset_free, track, measured, delay, gain, handover);
        } else {
            return designed<!Sampled>(sys, Ts, weights, integ, limits, filter, offset_free, track, measured, delay, gain, handover);
        }
    };
    if constexpr (::damp::detail::option_count<is_integral_option, Opts...> > 0) { // there without a .Q too: integrators_of refuses it
        const auto action = ::damp::detail::option<is_integral_option>(Unset{}, opts...);
        limits.integral = action.limits.template as<T>();
        const auto integ = integrators_of<NY, T>(action);
        auto       out = design(integ ? *integ : Integrators<NI, NY, T>{});
        if (!integ) {
            out = ::damp::detail::refused<decltype(out)>(DesignStatus::bad_weights); // an integrated output out of range or listed twice
        }
        return out;
    } else {
        return design(Integrators<0, NY, T>{});
    }
}

/** @brief Whether a disturbance spec fits the plant: a DisturbanceModel on its states and outputs, in its scalar. */
template<typename Spec, size_t NX, size_t NY, typename T>
inline constexpr bool disturbance_fits = true;

/** @brief A DisturbanceModel fits a plant of its sizes and scalar (NY 0: C_d left out). */
template<size_t NXD, size_t NYD, size_t ND, size_t NZ, typename U, size_t NX, size_t NY, typename T>
inline constexpr bool disturbance_fits<DisturbanceModel<NXD, NYD, ND, NZ, U>, NX, NY, T> = (NXD == NX) && (NYD == NY || NYD == 0) && damp::is_same_v<U, T>;

/**
 * @brief design::state_feedback in the plant's own domain (@p Sampled false) or at @p Ts: the LQ
 *        gain of @p lq (@p weights LQRWeights) or placed at Poles, on the measured state, or on a
 *        filter's estimate with a KalmanNoise or ObserverPoles among @p opts.
 */
template<bool Sampled, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Gain, typename... Opts>
[[nodiscard]] constexpr auto state_feedback_design_of(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        Ts,
    const LQRWeights<NX, NU, T>&             lq,
    const Gain&                              weights,
    const Opts&... opts
) {
    check_held_outputs<NU, NY, integrator_count<Opts...>(), Opts...>(); // the target's: one TrackedOutputs, none beside a disturbance spec, at most NU held
    state_feedback_options::check<Opts...>();
    static_assert(
        ::damp::detail::option_count<is_weight_spec, Opts...> == 0,
        "The weights given twice: LQRWeights{...} or Poles{...} is the one argument after the plant (and Ts), state_feedback(plant, [Ts,] weights, opts...)."
    );
    constexpr bool kFilter = ::damp::detail::option_count<is_filter_spec, Opts...> > 0;
    constexpr bool kIntegral = ::damp::detail::option_count<is_integral_option, Opts...> > 0;
    constexpr bool kDisturbance = ::damp::detail::option_count<is_disturbance, Opts...> > 0;
    constexpr bool kDelay = ::damp::detail::option_count<is_input_delay, Opts...> + ::damp::detail::option_count<is_output_delay, Opts...> > 0;
    constexpr bool kLate = ::damp::detail::option_count<is_output_delay, Opts...> > 0;
    constexpr bool kFits = disturbance_fits<::damp::detail::option_t<is_disturbance, NoDisturbance, Opts...>, NX, NY, T>;
    constexpr bool kTracked = ::damp::detail::option_count<is_tracked_option, Opts...> > 0;
    static_assert(kFilter || !kDisturbance, "A disturbance spec is estimated by the law's filter: give a KalmanNoise or ObserverPoles among the options.");
    static_assert(!(kIntegral && kDisturbance), "Integral action and a disturbance model are two offset-free choices: give an IntegralAction or a disturbance spec, not both.");
    static_assert(!(kIntegral && kTracked), "Integral action already holds r: TrackedOutputs is for a law without it (LQR, LQG).");
    static_assert(kFits, "The DisturbanceModel is on the plant's states and outputs, in the plant's scalar.");
    static_assert(Sampled || !kDelay, "A delay is sampled with the plant: give the continuous plant and Ts, fn(plant, Ts, config, ..., InputDelay<N>{theta}, OutputDelay<M>{theta}).");
    static_assert(!kLate || kFilter, "An OutputDelay is a late reading: it needs a filter to carry it, a KalmanNoise or ObserverPoles.");
    constexpr bool   kOptionsFit = (kFilter || !kDisturbance) && !(kIntegral && (kDisturbance || kTracked)) && kFits && (Sampled || !kDelay) && (!kLate || kFilter);
    constexpr size_t NK = NX + integrator_count<Opts...>();                     // the poles Poles lists
    constexpr size_t NKF = NX + disturbance_states_in<NX, NU, NY, NW, Opts...>; // and ObserverPoles
    constexpr size_t NF = listed_poles<::damp::detail::option_t<is_observer_poles, Unset, Opts...>>;
    constexpr bool   kGainCount = damp::is_same_v<Gain, LQGain> || (listed_poles<Gain> == NK);
    constexpr bool   kFilterCount = (::damp::detail::option_count<is_observer_poles, Opts...> == 0) || (NF == NKF);
    if constexpr (!kGainCount || !kFilterCount || !kOptionsFit) {
        if constexpr (!kGainCount) {
            static_cast<void>(controller_pole_count<NK, listed_poles<Gain>>());
        } else if constexpr (!kFilterCount) {
            static_cast<void>(observer_pole_count<NKF, NF>());
        }
        static_cast<void>(weights);
        static_cast<void>(sys);
        static_cast<void>(Ts);
        static_cast<void>(lq);
        (static_cast<void>(opts), ...);
        return StateFeedbackResult<NX, NU, T>{}; // the message says which: nothing further
    } else {
        return state_feedback_counted<Sampled>(sys, Ts, lq, gain_spec_of<T>(weights), opts...);
    }
}

} // namespace detail

/**
 * @brief The linear-quadratic state-feedback law (DARE on a discrete plant, CARE on a continuous one): runs as StateFeedback, an LQR, LQI, LQG, LQGI or OFLQG as its options say.
 *
 * The cost is @p weights in the plant's own time domain: LQRWeights {.Q, .R, .N}, the same with
 * or without integral action, which brings its own cost. Poles{…} in their place places the gain
 * instead, with the same options (state_feedback(sys, poles, opts…)). @p opts trail in any order,
 * at most one of each kind:
 * - InputLimits{.lower, .upper}: the limits the law's command runs within, clamped into them
 *   (InputLimits).
 * - InputRateLimits{.lower, .upper}: the limits on the command's move per period: each tick the
 *   command is clamped into them about the input applied last tick, and into the InputLimits
 *   (InputRateLimits). The law's memory holds against both clamps.
 * - DemandLimit{.state, .outer, .lower, .upper}: a limit on what the outer terms, listed by index (the plant
 *   states, then the integrators), ask of one state (DemandLimit). The gain is the same without
 *   either.
 * - IntegralAction{.Q, .outputs, .limits, .tracking_time}: integrators on the listed
 *   outputs y = C·x + D·u, in that order (r follows it), every output when it lists none, with
 *   their cost .Q (NI × NI, a scalar for one), their limits and tracking time: runs as LQI,
 *   xi += Ts·(r − Ci·y).
 * - TrackedOutputs{{i, …}}: the listed outputs held at r through the model's steady-state target
 *   (x* = Nₓ·r, u* = Nᵤ·r): step(r, x). Refused with DesignStatus::singular when the steady state is
 *   not unique (NZ = NU) or none holds every r (NZ < NU); more tracked outputs than inputs does not
 *   compile, as no steady state holds them all. A model error or a constant disturbance leaves an
 *   offset: an IntegralAction removes it.
 * - KalmanNoise{.process, .measurement[, .disturbance]}: the law runs on a steady-state Kalman
 *   filter's estimate of the state (design::kalman's): LQG, LQGI with an IntegralAction. The noise is
 *   sized to the plant: NX × NX process noise (the columns of G when it names them), NY × NY readings.
 * - ObserverPoles{…}, in place of the KalmanNoise: the same filter on the same model, its gain
 *   placed. The poles are those of its error dynamics, (I − M·C)·A in the current form it runs
 *   (the same as A − L·C's, L = A·M), one per plant state and disturbance state; the delay
 *   states take theirs at 0. Mixed with either weights: LQRWeights with ObserverPoles, Poles with
 *   a KalmanNoise.
 * - A disturbance spec, with a KalmanNoise: runs as OFLQG, the filter on [x; d] and the target under
 *   d̂. It says where d enters: InputDisturbance (B_d = B), OutputDisturbance (a bias on each
 *   reading), PlantDisturbance (the plant's own G and H) or DisturbanceModel (explicit, with the
 *   outputs held offset-free); the noise's .disturbance is the random walk of d (ND × ND). The filter's model is
 *   @f[
 *     x_{k+1} = A x_k + B u_k + B_d d_k,\qquad d_{k+1} = d_k,\qquad y_k = C x_k + C_d d_k + D u_k ,
 *   @f]
 *   and the target solves
 *   @f[
 *     \begin{bmatrix} A-I & B \\ C_z & D_z \end{bmatrix}
 *     \begin{bmatrix} x_t \\ u_t \end{bmatrix}
 *     =
 *     \begin{bmatrix} -B_d\hat d \\ r - C_{d,z}\hat d \end{bmatrix},
 *   @f]
 *   so the law is @f$ u = u_t - K(\hat x - x_t) @f$ with the plain gain. The augment must be
 *   detectable (ND ≤ NY), else DesignStatus::not_detectable; a target needs a unique steady state
 *   for NZ = NU, else DesignStatus::singular, and more held outputs than inputs (a spec that lists
 *   none, on a plant with more outputs than inputs) does not compile. An IntegralAction and a
 *   disturbance spec are two offset-free choices, refused together.
 * - MeasuredDisturbance{.Bd}: a disturbance the law reads each tick, entering the plant's state
 *   through B_d (continuous with Ts, ZOH'd with the plant: d held over the period). The gain is the
 *   same; the target is the steady state under d, so the law is u = u_t(r, d) − K·(x − x_t(r, d)),
 *   the target the MPC poses its QP about. For a d held constant, as the MPC holds it over its
 *   horizon, that is the LQ-optimal law in the deviation from that steady state: no term beyond the
 *   target. A law whose r is not the target's takes the steady state under d that holds at 0 every
 *   output (a regulator, control(d)) or its integrated outputs (LQI, LQGI: d is fed forward, and the
 *   integrators hold r and take up what the model leaves). That is at most NU outputs, as any
 *   target's: a regulator on a plant with more outputs than inputs does not compile under a
 *   MeasuredDisturbance, and lists the ones it holds, TrackedOutputs{{i, …}} with r = 0. A filter
 *   takes d as a known input after u, its model's B [B, B_d] and its noise unchanged, and with a
 *   disturbance spec the target is under [d̂; d], so d̂ estimates only what d does not account for.
 *   The runtime takes d after the reference and its feedforward: control(r, d), step(r, d, m),
 *   control(r, x_ref, u_ref, d), step(r, x_ref, u_ref, d, m). With an InputDelay d enters the
 *   plant, not the delay states.
 * - Handover{.fade_time}: the first command after the runtime's track(u, …) is u, the mismatch
 *   with the law's own decaying by exp(−Ts/fade_time) each tick (Handover), for the laws with
 *   nothing to preload: LQR, LQG and an OFLQG's transient (an LQI's or LQGI's preload is exact with
 *   NI ≥ NU; with it they fade what the preload leaves). The gain is the same; a fade_time that is not
 *   finite and positive is refused with DesignStatus::invalid_parameter.
 * - InputDelay<N>{θ}, on the sampled form state_feedback(plant, Ts, weights, …): the plant's input
 *   arrives θ late. The weights and the noise are the plant's, as without the delay; the cost is
 *   discretized exactly on the delayed model; a law on the measured state holds the inputs still on
 *   their way, a filter carries them. An InputDisturbance adds to u, so it reaches the plant through
 *   the delay states; a load on the plant itself, which acts at once, is DisturbanceModel{.Bd =
 *   sys.B, .Cd = sys.D}.
 * - OutputDelay<M>{θ}, on the sampled form with a KalmanNoise (LQG, LQGI, OFLQG): the readings
 *   arrive θ late. The filter is design::kalman(plant, Ts, noise, …, OutputDelay<M>{θ})'s, which
 *   carries the delay states of the readings already taken, so the estimate is the current
 *   state; the gain is the plant's (the input-delayed plant's with an InputDelay too) and reads
 *   none of those delay states, and the target holds each at the steady output. update(y) takes
 *   the reading as it arrives. An LQGI's integrators integrate that late reading, so they agree
 *   with the gain's model at a steady state, not in a transient. A law on the measured state has
 *   no filter to carry the late readings and is refused at compile time.
 *
 * @code
 * constexpr auto lqi = design::state_feedback(plant, 0.01, {.Q = Q, .R = {R}}, IntegralAction{.Q = 80.0}, InputDelay<3>{0.022});
 * constexpr auto lqg = design::state_feedback(plant, 0.01, {.Q = Q, .R = {R}}, KalmanNoise{.process = Qn, .measurement = 0.01}, OutputDelay<2>{0.015});
 * constexpr auto oflqg = design::state_feedback(plant, 0.01, {.Q = Q, .R = {R}}, KalmanNoise{.process = Qn, .measurement = 0.01, .disturbance = 0.25}, InputDisturbance{});
 * static constinit Law<lqi, float> law;
 * // each tick: law.update(x); u = law.control(r); write(sat(u)); law.back_calculate(sat(u)); law.predict();
 * @endcode
 * @note Compare with MATLAB®'s lqr(sys, Q, R, N), lqi, lqg(sys, QXU, QWV) and lqgtrack (the matlab::
 *       aliases); the offset-free law has none.
 * @see "Disturbance Models for Offset-Free Model-Predictive Control" (Pannocchia & Rawlings,
 *      2003), AIChE Journal 49(2)
 * @see "Offset-free reference tracking with model predictive control" (Maeder, Borrelli &
 *      Morari, 2009), Automatica 45(9)
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename... Opts>
[[nodiscard]] constexpr auto state_feedback(
    const StateSpace<NX, NU, NY, T, NW, NV>&            sys,
    const damp::type_identity_t<LQRWeights<NX, NU, T>>& weights,
    const Opts&... opts
) {
    return detail::state_feedback_design_of<false>(sys, T{0}, weights, detail::LQGain{}, opts...);
}

/**
 * @brief The sampled-data linear-quadratic state-feedback law: the continuous plant, cost and noise densities, the law at @p Ts (ZOH, Van Loan), with the options of state_feedback(sys, weights, opts…).
 * @note Compare with MATLAB®'s lqrd(sys, Q, R, Ts) (matlab::lqrd).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename... Opts>
[[nodiscard]] constexpr auto state_feedback(
    const StateSpace<NX, NU, NY, T, NW, NV>&            sys,
    damp::type_identity_t<T>                            Ts,
    const damp::type_identity_t<LQRWeights<NX, NU, T>>& weights,
    const Opts&... opts
) {
    return detail::state_feedback_design_of<true>(sys, Ts, weights, detail::LQGain{}, opts...);
}

/**
 * @brief The state-feedback law placed at @p poles: the gain with eig(A − B·K) where they say, and
 *        every option of state_feedback(sys, weights, opts…), in the plant's own time domain.
 *
 * In place of the weights, Poles{…} (design::place's machinery: KNV, the Jordan construction,
 * Ackermann's gain for one input). The poles are the plant's states', then the integrators' of an
 * IntegralAction, z-plane on a discrete plant (s-plane, a continuous gain for analysis, on a
 * continuous one with no option that samples it); a wrong count does not compile. The options
 * are the LQ design's and mean the same: the target, the limits, the integrators (their .Q, a cost,
 * is not read: it sizes them), a filter (a KalmanNoise, or ObserverPoles for a placed one), a
 * disturbance spec, a MeasuredDisturbance. The result is the LQ design's type, S zero, and runs
 * on the same runtime.
 *
 * Refused with DesignStatus::invalid_parameter for a pole that is not finite, not_placeable for a
 * mode the input does not reach or a pole repeated more times than placement takes (design::place),
 * singular for a rank-deficient B; the other statuses as the LQ design's.
 *
 * @code
 * constexpr auto placed = design::state_feedback(discrete_plant, design::Poles{0.9, 0.95});
 * static constinit Law<placed, float> law; // the LQR of the same shape
 * @endcode
 * @note Compare with MATLAB®'s place(A, B, p) and acker, then reg or lqgtrack's law around it.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t N, typename U, typename... Opts>
[[nodiscard]] constexpr auto state_feedback(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const Poles<N, U>& poles, const Opts&... opts) {
    return detail::state_feedback_design_of<false>(sys, T{0}, LQRWeights<NX, NU, T>{}, poles, opts...);
}

/**
 * @brief The sampled-data state-feedback law placed at @p poles: s-plane poles in rad/s, mapped
 *        to z = e^{s·Ts} on the ZOH model, with the options of state_feedback(sys, Ts, weights, opts…).
 *
 * The gain is design::place(sys, poles, Ts)'s on the plant, and on the law's model with its
 * integrators and an InputDelay's delay states otherwise: the delay states, which only hold the
 * commands still on their way, take their poles at 0 (deadbeat), on top of the list. With one
 * input the gain is unique; with several, a pole is placed at most NU times, so an InputDelay of
 * N ≥ 2 delay states per input is refused with not_placeable.
 *
 * @code
 * constexpr auto lqi = design::state_feedback(plant, 0.01, design::Poles{-8.0, -12.0, -4.0}, IntegralAction{.Q = 80.0});
 * constexpr auto adrc_like = design::state_feedback(chain, 0.001, design::Poles{-8.0, -8.0}, design::ObserverPoles{-40.0, -40.0, -40.0}, design::InputDisturbance{});
 * @endcode
 * @note Compare with MATLAB®'s place(A, B, exp(p·Ts)) on c2d(sys, Ts).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t N, typename U, typename... Opts>
[[nodiscard]] constexpr auto state_feedback(const StateSpace<NX, NU, NY, T, NW, NV>& sys, damp::type_identity_t<T> Ts, const Poles<N, U>& poles, const Opts&... opts) {
    return detail::state_feedback_design_of<true>(sys, Ts, LQRWeights<NX, NU, T>{}, poles, opts...);
}

} // namespace design
} // namespace damp
