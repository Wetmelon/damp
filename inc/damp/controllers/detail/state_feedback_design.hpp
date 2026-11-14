// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file state_feedback_design.hpp
 * @brief The design behind lqr, lqi, lqg, lqgi and oflqg: detail::state_feedback, the gain on the
 *        plant or on its model with input registers, then the filter and the target the options ask for
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/lq_cost.hpp"
#include "damp/controllers/detail/state_feedback_options.hpp"
#include "damp/controllers/law_detail.hpp"
#include "damp/controllers/state_feedback_types.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/status.hpp"
#include "damp/design/steady_state_target.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/decomposition.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace design {
namespace detail {

/** @brief The gain on the law's state. */
template<typename Result, size_t NK, size_t NU, typename T>
[[nodiscard]] constexpr Result with_gain(Result out, const LQRResult<NK, NU, T>& lq) {
    out.K = lq.K;
    out.S = lq.S;
    out.e = lq.e;
    return ::damp::detail::set_status(out, lq.status);
}

/**
 * @brief The law on a measured state: the gain on [x; xi] for the cost of @p w, and the
 *        integrators of @p Ci (NI > 0). A continuous plant with no Ts and no integrators gets the
 *        continuous law (CARE).
 */
template<bool Discrete, size_t NI, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr auto law_design(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T Ts, const LQRConfig<NX, NU, NI, T>& w, const Matrix<NI, NY, T>& Ci) {
    using Result = StateFeedbackResult<NX, NU, T, NI, ((NI > 0) ? NY : NX)>;
    if (!mat::cholesky(w.R).has_value()) {
        return ::damp::detail::refused<Result>(DesignStatus::bad_weights); // R not positive definite, on every route
    }
    if constexpr (Discrete && NI == 0) {
        if (sys.Ts == T{0}) {
            const auto lq = continuous_lqr(sys.A, sys.B, w.Q, w.R, w.N);
            return with_gain(Result{}, lq);
        }
    }
    const auto m = sampled_model(sys, Ts, Discrete);
    if (m.status != DesignStatus::ok) {
        return ::damp::detail::refused<Result>(m.status);
    }
    const T h = Discrete ? sys.Ts : Ts;
    if constexpr (NI == 0) {
        LQRCost<NX, NU, T> cd{w.Q, w.R, w.N};
        if constexpr (!Discrete) {
            const auto c = discretize_lqr_cost(sys.A, sys.B, w.Q, w.R, Ts, w.N);
            if (!c) {
                return ::damp::detail::refused<Result>(DesignStatus::bad_weights); // the plant was sampled: a non-finite weight
            }
            cd = *c;
        }
        return with_gain(Result{.Ts = h}, discrete_lqr(m.law.A, m.law.B, cd));
    } else {
        if (!mat::cholesky(w.Q_integral).has_value()) {
            return ::damp::detail::refused<Result>(DesignStatus::bad_weights); // an unweighted integrator: a mode at z = 1 the cost does not see
        }
        // Forward-Euler integrator xi += h·e (MATLAB®'s discrete lqi): xi ≈ ∫e dt, so Q_integral
        // means the same on either route. The runtime integrates e = r − Ci·(C·x + D·u), u the
        // input applied over the period, and so does this model on both routes.
        const auto              im = integral_model(m.law.A, m.law.B, m.law.C, m.law.D, Ci, h);
        LQRCost<NX + NI, NU, T> cd = stacked_weights(w);
        if constexpr (!Discrete) {
            const auto c = sampled_integral_cost(sys.A, sys.B, sys.C, sys.D, Ci, cd, Ts);
            if (!c) {
                return ::damp::detail::refused<Result>(DesignStatus::bad_weights); // the plant was sampled: a non-finite weight
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
 *        and the integrators of @p Ci (NI > 0).
 */
template<size_t NI, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t N, typename S>
[[nodiscard]] constexpr auto delayed_law_design(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        Ts,
    const InputDelay<N, S>&                  delay,
    const LQRConfig<NX, NU, NI, T>&          w,
    const Matrix<NI, NY, T>&                 Ci
) {
    using Result = StateFeedbackResult<NX, NU, T, NI, ((NI > 0) ? NY : NX), NoEstimator, 0, N>;
    if (!mat::cholesky(w.R).has_value()) {
        return ::damp::detail::refused<Result>(DesignStatus::bad_weights); // R not positive definite
    }
    if constexpr (NI > 0) {
        if (!mat::cholesky(w.Q_integral).has_value()) {
            return ::damp::detail::refused<Result>(DesignStatus::bad_weights); // an unweighted integrator
        }
    }
    const auto m = sampled_model(sys, Ts, false, delay); // a continuous plant, a usable Ts, finite entries, a delay discretize takes
    if (m.status != DesignStatus::ok) {
        return ::damp::detail::refused<Result>(m.status);
    }
    const auto& model = m.law;
    const auto  inputs = ::damp::detail::channel_register_delays<NU>(delay.seconds, Ts, N); // as discretize accepted them
    const auto  cost = inputs ? delayed_cost(sys, model, *inputs, Ci, w, Ts) : damp::nullopt;
    if (!cost) {
        return ::damp::detail::refused<Result>(DesignStatus::bad_weights); // a non-finite weight on the sampled plant
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

/** @brief The gain of a design: on the plant (law_design), or on the model with the input registers of @p delay (delayed_law_design). */
template<bool Discrete, size_t NI, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Delay>
[[nodiscard]] constexpr auto gain_design(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T Ts, const LQRConfig<NX, NU, NI, T>& w, const Matrix<NI, NY, T>& Ci, const Delay& delay) {
    if constexpr (damp::is_same_v<Delay, NoDelay>) {
        static_cast<void>(delay);
        return law_design<Discrete, NI>(sys, Ts, w, Ci);
    } else {
        static_assert(!Discrete, "A delay is sampled with the plant: a continuous plant and Ts.");
        return delayed_law_design<NI>(sys, Ts, delay, w, Ci);
    }
}

/**
 * @brief The design on checked options: the law, then the filter and the target they ask for, on
 *        the plant or on its model with the input registers of @p delay; status the first failure.
 */
template<bool Discrete, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Options, typename Delay = NoDelay>
[[nodiscard]] constexpr auto designed(
    const StateSpace<NX, NU, NY, T, NW, NV>&  sys,
    T                                         Ts,
    const config_for<NX, NU, NY, T, Options>& w,
    const Options&                            opts,
    const Delay&                              delay = Delay{}
) {
    using O = options_traits<Options>;
    constexpr size_t NI = integrator_count<NY, Options>;
    constexpr size_t NRU = delay_registers<Delay>;
    constexpr bool   kEstimated = !is_unset<typename O::noise>;
    constexpr bool   kTarget = !is_unset<typename O::disturbance> || !is_unset<typename O::track>;
    const auto       spec = disturbance_of(opts);

    const auto Ci = integrated_selector<NI, NY, T>(integral_action_of(opts));
    auto       law = gain_design<Discrete>(sys, Ts, w, Ci ? *Ci : Matrix<NI, NY, T>{}, delay);
    if (!Ci) {
        law = ::damp::detail::refused<decltype(law)>(DesignStatus::bad_weights); // an integrated output out of range or listed twice
    } else if (law.success && !(w.u_lim.valid() && w.integral_lim.valid())) {
        law = ::damp::detail::refused<decltype(law)>(DesignStatus::invalid_parameter); // a limit that is NaN, or a lower one above its upper one
    }
    law.u_lim = w.u_lim; // the runtime's: the gain does not depend on them
    if constexpr (NI > 0) {
        law.integral.lim = w.integral_lim;
    }
    if constexpr (!kEstimated && !kTarget) {
        return law;
    } else {
        using SpecT = damp::remove_cvref_t<decltype(spec)>;
        constexpr size_t NDS = disturbance_traits<SpecT, NU, NY, NW>::listed;
        constexpr size_t NZS = (NDS == 0) ? NY : NDS; // the spec's outputs, or every output
        constexpr size_t NZ = kTarget ? NZS : 0;
        const auto       filter = [&] {
            if constexpr (kEstimated) {
                return kalman_design<Discrete>(sys, Ts, opts.noise, spec, delay);
            } else {
                return NoEstimator{};
            }
        }();
        using Result = StateFeedbackResult<NX, NU, T, NI, (kEstimated ? NY : NX), damp::remove_cvref_t<decltype(filter)>, NZ, NRU>;
        Result out{law.K, law.S, law.e, law.success, law.status, law.Ts, filter};
        if constexpr (NI > 0) {
            out.integral = law.integral;
        }
        out.u_lim = law.u_lim;
        DesignStatus status = law.status;
        if constexpr (kEstimated) {
            status = (status == DesignStatus::ok) ? filter.status : status;
        }
        if constexpr (kTarget) {
            out.target = target_design<Discrete>(sys, Ts, spec, target_weight<SpecT, NU, NY, NW, T>::identity(), delay);
            status = (status == DesignStatus::ok) ? out.target.status : status;
        }
        return ::damp::detail::set_status(out, status);
    }
}

/**
 * @brief The state-feedback design on @p sys in its own time domain (@p Discrete) or sampled at
 *        @p Ts, with its input @p delay late when that is an InputDelay, once the options pass
 *        their checks: a mistake stops here with its message alone.
 */
template<bool Discrete, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Options, typename Delay = NoDelay>
[[nodiscard]] constexpr auto state_feedback_design(
    const StateSpace<NX, NU, NY, T, NW, NV>&  sys,
    T                                         Ts,
    const config_for<NX, NU, NY, T, Options>& w,
    const Options&                            opts,
    const Delay&                              delay = Delay{}
) {
    if constexpr (!options_traits<Options>::named) {
        static_assert(options_traits<Options>::named, "The options are design::StateFeedbackOptions{.noise = ..., .offset_free = ...}: one named field per feature.");
        return StateFeedbackResult<NX, NU, T>{};
    } else if constexpr (!fields_valid<NX, NY, T, Options>()) {
        return StateFeedbackResult<NX, NU, T>{}; // fields_valid says which field
    } else if constexpr (!features_valid<NU, NY, T, NW, NV, Options>()) {
        return StateFeedbackResult<NX, NU, T>{}; // features_valid says why
    } else {
        return designed<Discrete>(sys, Ts, w, opts, delay);
    }
}

/**
 * @brief The linear-quadratic state-feedback law for a plant, in the plant's own time domain
 *
 * The law u = −K·x minimizes the weighted cost of @p config; @p options adds the features of the
 * law, one named field each (StateFeedbackOptions). A discrete plant (@c sys.Ts > 0) gets the
 * discrete law (DARE) at its own Ts; a continuous plant with no options gets the continuous law
 * (CARE, Ts = 0), for analysis: only a discrete law deploys. For a digital law from a continuous
 * plant, pass the sample time: state_feedback(sys, Ts, config, options).
 *
 * @f[
 *   J = \sum_k \bigl( x_k^\top Q x_k + u_k^\top R u_k + 2 x_k^\top N u_k + \xi_k^\top Q_\xi \xi_k \bigr)
 * @f]
 *
 * With `.offset_free = IntegralAction{}`, the integrators ξ advance by xi += Ts·(r − Ci·y)
 * (MATLAB®'s discrete lqi, so ξ ≈ ∫e dt and Q_integral means what it does on the sampled route).
 * An input delay on a discrete plant is part of its model: design on the one
 * discretize(plant, Ts, InputDelay<N>{θ}) returns, its registers states like any other, or from
 * the continuous plant with state_feedback(sys, Ts, InputDelay<N>{θ}, config, options).
 * The other options are state_feedback(sys, Ts, config, options)'s.
 *
 * Refused with DesignStatus::invalid_sample_time for a negative or non-finite sys.Ts, or a
 * continuous plant with any option; with DesignStatus::bad_weights for a listed output out
 * of range or listed twice.
 *
 * @note Compare with MATLAB®'s lqr(sys, Q, R, N), lqi, lqg and lqgtrack, which also follow the
 *       plant's domain.
 *
 * @param sys      Plant (A, B, C, D, G, H and Ts)
 * @param config   The cost {.Q, .R, .N, .Q_integral}, in the plant's time domain, and the limits {.u_lim, .integral_lim}
 * @param options  design::StateFeedbackOptions{…}: the law's features (default none)
 * @return StateFeedbackResult: the runtime its design names runs it
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Options>
[[nodiscard]] constexpr auto state_feedback(
    const StateSpace<NX, NU, NY, T, NW, NV>&  sys,
    const config_for<NX, NU, NY, T, Options>& config,
    const Options&                            options
) {
    return state_feedback_design<true>(sys, T{0}, config, options);
}

/**
 * @brief The sampled-data state-feedback law: a continuous plant and its cost, the digital law at
 *        @p Ts
 *
 * Discretizes the dynamics (ZOH) and the cost integral exactly (Van Loan, discretize_lqr_cost),
 * then solves the discrete problem, so the weights mean the same at any rate. @p options adds
 * the law's features, one named field each:
 *
 * | Option | The law |
 * | --- | --- |
 * | `.noise = KalmanNoise{.process, .measurement}` | estimates x with a Kalman filter from y (design::kalman(sys, Ts, noise)): the process noise a density, the measurement a variance per reading |
 * | `.offset_free = IntegralAction{}` or `IntegralAction<NI>{{i, …}}` | integrates r − Ci·y: ξ̇ = r − Ci·y, so ξ = ∫e dt and Q_integral holds at any rate |
 * | `.offset_free = InputDisturbance{}`, … | estimates a disturbance d and holds r through the steady-state target under it (OFLQG, below; needs `.noise`) |
 *
 * @par Integral action
 * The integrator advances by xi += Ts·(r − Ci·y) on each sampled error, and the cost integral
 * is taken over each sample with the plant running on the held input and ξ ramping on the held
 * error (Van Loan), so the closed-loop poles on the result are those of the deployed loop.
 * Constant references and disturbances leave no steady-state error on the integrated outputs.
 *
 * @par Input delay
 * A delay is a model of the plant: state_feedback(sys, Ts, InputDelay<N>{θ}, config, options)
 * designs on discretize(sys, Ts, InputDelay<N>{θ}) with the cost discretized exactly, every
 * option as here. A measured-state law holds the registers v(k−1) … v(k−N); a filter carries
 * them as states of its model, with no process noise on them.
 *
 * @par Offset-free control by a disturbance model
 * A disturbance model in `.offset_free` (InputDisturbance, OutputDisturbance, PlantDisturbance,
 * DisturbanceModel) augments the filter's model with integrating disturbances,
 * @f[
 *   x_{k+1} = A x_k + B u_k + B_d d_k,\qquad d_{k+1} = d_k,\qquad y_k = C x_k + C_d d_k + D u_k ,
 * @f]
 * the filter estimates @f$ [\hat x;\, \hat d] @f$ (the random walk of d is noise.disturbance, positive definite),
 * and the steady-state target holds the controlled outputs z (the rows the spec lists, all
 * without a list) at r under d̂:
 * @f[
 *   \begin{bmatrix} A-I & B \\ C_z & D_z \end{bmatrix}
 *   \begin{bmatrix} x_t \\ u_t \end{bmatrix}
 *   =
 *   \begin{bmatrix} -B_d\hat d \\ r - C_{d,z}\hat d \end{bmatrix},
 * @f]
 * so the law is @f$ u = u_t - K(\hat x - x_t) @f$ with the plain gain. Where d enters: an input
 * disturbance adds to u (B_d = B), an output disturbance is a bias on each reading (C_d = I),
 * PlantDisturbance is the plant's own G and H, DisturbanceModel is explicit. The augment must be
 * detectable (ND ≤ NY), else DesignStatus::not_detectable; a target needs a unique steady state
 * for NZ = NU, else DesignStatus::singular.
 *
 * Refused with DesignStatus::invalid_sample_time when @p Ts is not finite and positive, or when
 * @p sys is already discrete, even at the same Ts (a second discretization is the bug this
 * signature rules out); with DesignStatus::bad_weights for a listed output out of range or listed twice. The status is
 * the first failure: the gain's, the filter's, then the target's.
 *
 * @code
 * const auto res = design::state_feedback(plant, Ts, {.Q = Q, .R = R},
 *     design::detail::StateFeedbackOptions{.noise = design::KalmanNoise{.process = Qn, .measurement = Rn, .disturbance = Qd},
 *                                          .offset_free = design::InputDisturbance{}});
 * static constinit OFLQG law{res.as<float>()}; // the public form: design::oflqg
 * // each tick: auto u = law.step(r, y);
 * @endcode
 *
 * @note Compare with MATLAB®'s lqrd(sys, Q, R, Ts); lqi, lqg and lqgtrack on a continuous plant,
 *       followed by c2d; and for the offset-free law, its mpc with a disturbance model, without
 *       a QP.
 * @see "Optimal Control" (Anderson & Moore, 1990), §4 and §8
 * @see "Computing Integrals Involving the Matrix Exponential" (Van Loan, 1978)
 * @see "Disturbance Models for Offset-Free Model-Predictive Control" (Pannocchia & Rawlings,
 *      2003), AIChE Journal 49(2)
 * @see "Offset-free reference tracking with model predictive control" (Maeder, Borrelli &
 *      Morari, 2009), Automatica 45(9)
 *
 * @param sys      Continuous plant (sys.Ts = 0; A, B, C, D, G, H)
 * @param Ts       Sampling time of the law [s]
 * @param config   The continuous cost {.Q, .R, .N, .Q_integral}, a cost integral, and the limits {.u_lim, .integral_lim}
 * @param options  design::StateFeedbackOptions{…}: the law's features (default none)
 * @return StateFeedbackResult on the discretized plant: the runtime its design names runs it
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Options>
[[nodiscard]] constexpr auto state_feedback(
    const StateSpace<NX, NU, NY, T, NW, NV>&  sys,
    damp::type_identity_t<T>                  Ts,
    const config_for<NX, NU, NY, T, Options>& config,
    const Options&                            options
) {
    return state_feedback_design<false>(sys, Ts, config, options);
}

/**
 * @brief The sampled-data state-feedback law of a continuous plant whose input arrives @p delay
 *        late: the weights and the noise on the plant, the law on [x; v(k−1); …; v(k−N)]
 *
 * The model is discretize(sys, Ts, delay), and the cost integral is discretized on it exactly:
 * the plant runs on the input it holds within each period, which switches at the delay's
 * fraction of a sample (delayed_cost). The registers take no weight of their own, so the weights
 * of @p config are the plant's, as for the design without a delay, and mean the same at any Ts
 * and any delay.
 *
 * The options are state_feedback(sys, Ts, config, options)'s, each on the delayed model. A
 * measured-state law (LQR, LQI) holds the registers: update takes the plant state, and predict
 * shifts in the input back_calculate recorded. A filter (`.noise`) is design::kalman(sys, Ts,
 * delay, noise): the process noise acts on the plant, none on the registers, which it carries
 * exactly. A target holds each register at u_t. x_ref is the plant's and the registers'
 * reference is u_ref, on every law.
 *
 * Refused as state_feedback(sys, Ts, config, options), and with discretize's status for a delay
 * it refuses (negative, not finite or longer than N samples).
 *
 * @param sys      Continuous plant (sys.Ts = 0)
 * @param Ts       Sampling time of the law [s]
 * @param delay    InputDelay<N>{θ}: θ ≤ N·Ts, one for every input or a ColVec of one per input
 * @param config   The continuous cost {.Q, .R, .N, .Q_integral} on the plant, a cost integral, and the limits {.u_lim, .integral_lim}
 * @param options  design::StateFeedbackOptions{…}: the law's features (default none)
 * @return StateFeedbackResult with NRU = N: the runtime its design names runs it
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t N, typename S, typename Options>
[[nodiscard]] constexpr auto state_feedback(
    const StateSpace<NX, NU, NY, T, NW, NV>&  sys,
    damp::type_identity_t<T>                  Ts,
    const InputDelay<N, S>&                   delay,
    const config_for<NX, NU, NY, T, Options>& config,
    const Options&                            options
) {
    return state_feedback_design<false>(sys, static_cast<T>(Ts), config, options, delay);
}

} // namespace detail
} // namespace design
} // namespace damp
