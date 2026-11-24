// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file state_feedback.hpp
 * @brief design::state_feedback: the linear-quadratic state-feedback designs, the law picked by its options
 *
 * design::state_feedback(plant, [Ts,] config, opts...) → StateFeedbackResult, which runs
 * as the law its options name (`Law<kRes, T>`, or `LQI law{res}`):
 *
 * | KalmanNoise among the options | Other options | Law |
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
 * This header includes the five runtimes, so the law a design names deploys from it. The matlab::
 * aliases (lqr, dlqr, lqrd, lqi, lqg, lqgtrack) take MATLAB®'s arguments and call through.
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/state_feedback_design.hpp"
#include "damp/controllers/detail/state_feedback_options.hpp"
#include "damp/controllers/law.hpp"
#include "damp/controllers/lqg.hpp"   // IWYU pragma: export
#include "damp/controllers/lqgi.hpp"  // IWYU pragma: export
#include "damp/controllers/lqi.hpp"   // IWYU pragma: export
#include "damp/controllers/lqr.hpp"   // IWYU pragma: export
#include "damp/controllers/oflqg.hpp" // IWYU pragma: export
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

/** @brief The cost and the integrated outputs of a design with an IntegralAction: what its runtime is built from. */
template<size_t NX, size_t NU, size_t NI, typename T>
struct IntegralCost {
    LQIConfig<NX, NU, NI, T> cost{};     ///< the full cost with the integrators' weight and limits
    IntegralAction<NI>       integral{}; ///< The integrated outputs, in order
};

/** @brief The full cost and an IntegralAction as the design reads them: one cost on [x; ξ], and the outputs listed (every one left out). */
template<size_t NY, size_t NX, size_t NU, typename T, typename TQ, typename TI>
[[nodiscard]] constexpr auto integral_cost(const LQRConfig<NX, NU, T>& config, const ::damp::design::IntegralAction<TQ, TI>& action) {
    constexpr size_t NI = ::damp::design::IntegralAction<TQ, TI>::NI;
    static_assert(
        NI > 0,
        "IntegralAction::Q is the cost on the integrators: IntegralAction{.Q = ...}, a scalar for one, an NI x NI Matrix for more. "
        "Given one inside a template? There clang deduces nothing from a .Q set to a literal or a constant (llvm/llvm-project#114870): "
        "name the type, IntegralAction<double>{.Q = ...}."
    );
    static_assert(square_value<TQ, NI>, "IntegralAction::Q is NI x NI, one integrator per integrated output: a Matrix, or a scalar for one.");
    static_assert(NI <= NY, "IntegralAction::Q has one integrator per integrated output: at most one per measured output.");
    IntegralCost<NX, NU, NI, T> out{};
    out.integral.outputs = action.outputs;
    if constexpr (NI == NY) {
        if (out.integral.outputs[0] == kEveryOutput) {
            for (size_t i = 0; i < NI; ++i) {
                out.integral.outputs[i] = i; // left out: every output, in order
            }
        }
    }
    out.cost.Q = config.Q;
    out.cost.R = config.R;
    out.cost.N = config.N;
    out.cost.Q_integral = square_matrix<NI, T>(action.Q, Matrix<NI, NI, T>{});
    out.cost.u_lim = config.u_lim;
    out.cost.integral_lim = action.limits.template as<T>();
    out.cost.tracking_time = static_cast<T>(action.tracking_time);
    return out;
}

/**
 * @brief @p law with a Handover: the same design, its fade_time the option's. Refused with
 *        DesignStatus::invalid_parameter for a fade_time that is not finite and positive (a design
 *        that failed keeps its own status).
 */
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, typename Kalman, size_t NZ, size_t NRU, size_t NM, size_t NDU>
[[nodiscard]] constexpr StateFeedbackResult<NX, NU, T, NI, NY, Kalman, NZ, NRU, NM, NDU, true> with_handover(
    const StateFeedbackResult<NX, NU, T, NI, NY, Kalman, NZ, NRU, NM, NDU>& law,
    const Handover&                                                         handover
) {
    StateFeedbackResult<NX, NU, T, NI, NY, Kalman, NZ, NRU, NM, NDU, true> out{law.K, law.S, law.e, law.success, law.status, law.Ts, law.kalman};
    if constexpr (NI > 0) {
        out.integral = law.integral;
        out.tracking_time = law.tracking_time;
    }
    if constexpr (NZ > 0 || NM > 0) {
        out.target = law.target;
    }
    out.limits = law.limits;
    out.fade_time = static_cast<T>(handover.fade_time);
    const bool usable = damp::finite_positive(handover.fade_time) && damp::finite_positive(out.fade_time); // and still > 0 in the plant's scalar
    return ::damp::detail::set_status(out, (law.status == DesignStatus::ok && !usable) ? DesignStatus::invalid_parameter : law.status);
}

/** @brief The options of a state-feedback design: a KalmanNoise makes it estimate the state. */
using state_feedback_options = ::damp::detail::
    option_set<is_input_delay, is_output_delay, is_integral_option, is_tracked_option, is_disturbance, is_measured_disturbance, is_kalman_noise, is_input_limits, is_input_rate_limits, is_demand_limit, is_handover>;

/**
 * @brief A design's result with the InputRateLimits among @p opts (NDU = NU), or itself without one:
 *        refused with DesignStatus::invalid_parameter when the box is not one, as the InputLimits are.
 */
template<size_t NU, typename T, typename Result, typename... Opts>
[[nodiscard]] constexpr auto with_rate_limits(const Result& res, const Opts&... opts) {
    if constexpr (rate_limited_inputs<NU, Opts...> == 0) {
        (static_cast<void>(opts), ...);
        return res;
    } else {
        const Bounds<NU, T> rate = input_rate_limits_of<NU, T>(opts...);
        const auto          out = res.with_input_rate(rate);
        if (res.success && !rate.valid()) {
            return ::damp::detail::refused<damp::remove_cvref_t<decltype(out)>>(DesignStatus::invalid_parameter); // an end NaN or ±Inf, or the lower above the upper
        }
        return out;
    }
}

/** @brief The integrators of a state-feedback design's options: an IntegralAction's NI, else none. */
template<typename... Opts>
consteval size_t integrators_of() {
    using Action = ::damp::detail::option_t<is_integral_option, Unset, Opts...>;
    if constexpr (is_unset<Action>) {
        return 0;
    } else {
        return Action::NI;
    }
}

/** @brief design::state_feedback in the plant's domain (@p Sampled false) or at @p Ts: on the measured state, or on a Kalman filter's estimate with a KalmanNoise among @p opts. */
template<bool Sampled, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename... Opts>
[[nodiscard]] constexpr auto state_feedback_design_of(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        Ts,
    const LQRConfig<NX, NU, T>&              config,
    const Opts&... opts
) {
    check_held_outputs<NU, NY, integrators_of<Opts...>(), Opts...>(); // the target's: one TrackedOutputs, none beside a disturbance spec, at most NU held
    state_feedback_options::check<Opts...>();
    constexpr bool kFilter = ::damp::detail::option_count<is_kalman_noise, Opts...> > 0;
    constexpr bool kIntegral = ::damp::detail::option_count<is_integral_option, Opts...> > 0;
    constexpr bool kDisturbance = ::damp::detail::option_count<is_disturbance, Opts...> > 0;
    static_assert(kFilter || !kDisturbance, "A disturbance spec is estimated by a Kalman filter: give its KalmanNoise among the options.");
    static_assert(!(kIntegral && kDisturbance), "Integral action and a disturbance model are two offset-free choices: give an IntegralAction or a disturbance spec, not both.");
    const auto track = ::damp::detail::option<is_tracked_option>(Unset{}, opts...);
    static_cast<void>(measured_of<NX, T>(opts...)); // its checks: one, on the plant's states
    const auto measured = ::damp::detail::option<is_measured_disturbance>(Unset{}, opts...);
    const auto run = [&](const auto& cost, const auto& offset_free) {
        if constexpr (kFilter) {
            const StateFeedbackOptions options{.noise = noise_of<NX, NU, NY, T, NW, NV>(opts...), .offset_free = offset_free, .track = track, .measured = measured};
            return state_feedback_with<Sampled>(noise_model(sys), Ts, cost, options, opts...);
        } else {
            const StateFeedbackOptions options{.offset_free = offset_free, .track = track, .measured = measured};
            return state_feedback_with<Sampled>(sys, Ts, cost, options, opts...);
        }
    };
    const auto law = [&] {
        if constexpr (kIntegral) {
            auto integral = integral_cost<NY>(config, ::damp::detail::option<is_integral_option>(Unset{}, opts...));
            integral.cost.demand_lim = demand_terms_of<decltype(integral.cost.demand_lim)>(opts...); // its outer terms may be the integrators
            return with_rate_limits<NU, T>(run(integral.cost, integral.integral), opts...);
        } else {
            return with_rate_limits<NU, T>(run(config, ::damp::detail::option<is_disturbance>(Unset{}, opts...)), opts...);
        }
    }();
    if constexpr (::damp::detail::option_count<is_handover, Opts...> == 0) {
        return law;
    } else {
        return with_handover(law, ::damp::detail::option<is_handover>(Handover{}, opts...));
    }
}

} // namespace detail

/**
 * @brief The linear-quadratic state-feedback law (DARE on a discrete plant, CARE on a continuous one): runs as LQR, LQI, LQG, LQGI or OFLQG, as its options say.
 *
 * The cost is @p weights in the plant's own time domain: LQRWeights {.Q, .R, .N}, the same with
 * or without integral action, which brings its own cost. @p opts trail in any order, at most one
 * of each kind:
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
    return detail::state_feedback_design_of<false>(sys, T{0}, detail::lq_cost_of(weights, opts...), opts...);
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
    return detail::state_feedback_design_of<true>(sys, Ts, detail::lq_cost_of(weights, opts...), opts...);
}

} // namespace design
} // namespace damp
