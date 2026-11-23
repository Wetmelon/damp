// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file config.hpp
 * @brief A helper over the trailing design options: one config struct per design, each option a
 *        named field
 *
 * A design takes its options as a trailing pack recognized by type (options.hpp), which is the
 * spelling of Damp's docs and examples:
 *
 * @code
 * design::state_feedback(plant, Ts, {.Q = Q, .R = 0.1}, design::IntegralAction{.Q = 80.0}, InputDelay<3>{0.022});
 * @endcode
 *
 * This header adds, for a caller who would name each option or keep a design's inputs in one
 * object, `fn(plant, [Ts,] config)`: the config a struct whose fields hold those same objects.
 * It passes the fields that were given on as the options, so it is the same design with the same
 * result:
 *
 * @code
 * design::state_feedback(plant, Ts, design::LQConfig{
 *     .weights = design::LQRWeights{.Q = Q, .R = {0.1}},
 *     .integral = design::IntegralAction{.Q = 80.0},
 *     .delay = InputDelay<3>{0.022},
 * });
 * @endcode
 *
 * | Design | Config | Fields, in this order |
 * | --- | --- | --- |
 * | state_feedback | LQConfig | weights (required), limits, demand, integral, noise, disturbance, tracked, delay |
 * | kalman | KalmanConfig | noise (required), disturbance, delay |
 * | mpc | MPCConfig | weights, constraints, noise, disturbance, measured |
 * | steady_state_target | TargetConfig | disturbance, delay |
 *
 * A field left out is off (MPC's weights and constraints: their defaults). A field's type is a
 * template parameter, so a config names no sizes, and the plant cannot size bare braces in it:
 * every field holds a named object, the weights too. `LQRWeights{.Q = Q, .R = {0.1}}` takes its states
 * and its scalar from Q and its inputs from R, one input when R is braced; write the sizes when
 * the fields do not give them (`LQRWeights<1, 1>{…}`). The config itself is named too
 * (`LQConfig{…}`): bare braces in its place are the pack's weights, as without this header.
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/state_feedback_options.hpp"
#include "damp/controllers/law_detail.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/controllers/state_feedback.hpp"
#include "damp/controllers/state_feedback_types.hpp"
#include "damp/design/detail/config.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/options.hpp"
#include "damp/design/steady_state_target.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace design {

/**
 * @brief design::state_feedback's weights and options as named fields
 *
 * | Field | Holds | Makes the law |
 * | --- | --- | --- |
 * | weights | LQRWeights{.Q, .R, .N}, in the plant's time domain (required) | LQR |
 * | limits | LQRLimits{.u_lim}: the command clamped | the same law, its command bounded |
 * | demand | DemandLimit{.state, .outer, .lim}: what the outer terms ask of one state, bounded | the same law, read as a cascade |
 * | integral | IntegralAction{.Q, .outputs, .lim, .tracking_time} | LQI |
 * | noise | KalmanNoise{.process, .measurement[, .disturbance]} | LQG (LQGI with integral) |
 * | disturbance | InputDisturbance{}, OutputDisturbance{}, PlantDisturbance{} or DisturbanceModel{...}, with a noise | OFLQG |
 * | tracked | TrackedOutputs{} or TrackedOutputs{{i, …}} | LQR or LQG holding the outputs at r |
 * | delay | InputDelay<N>{θ}, with Ts | the same law on the plant whose input arrives θ late |
 */
template<
    typename TWeights = detail::Unset,
    typename TLimits = detail::Unset,
    typename TDemand = detail::Unset,
    typename TIntegral = detail::Unset,
    typename TNoise = detail::Unset,
    typename TDisturbance = detail::Unset,
    typename TTracked = detail::Unset,
    typename TDelay = detail::Unset>
struct LQConfig {
    TWeights     weights{};     ///< LQRWeights{.Q, .R, .N}: the cost (required)
    TLimits      limits{};      ///< LQRLimits{.u_lim}: the limit the command runs within
    TDemand      demand{};      ///< DemandLimit{.state, .outer, .lim}: a limit on what the outer terms ask of one state
    TIntegral    integral{};    ///< IntegralAction{.Q, …}: integrators on r − y (LQI)
    TNoise       noise{};       ///< KalmanNoise{.process, .measurement}: the law on a Kalman filter's estimate (LQG)
    TDisturbance disturbance{}; ///< Where an integrating disturbance enters, with a noise (OFLQG)
    TTracked     tracked{};     ///< TrackedOutputs{…}: the outputs held at r through the steady-state target
    TDelay       delay{};       ///< InputDelay<N>{θ}: the plant's input arrives θ late (with Ts)
};

/**
 * @brief design::kalman's options as named fields
 *
 * | Field | Holds |
 * | --- | --- |
 * | noise | KalmanNoise{.process, .measurement[, .disturbance]} (required) |
 * | disturbance | InputDisturbance{}, OutputDisturbance{}, PlantDisturbance{} or DisturbanceModel{...}: the filter on [x; d] |
 * | delay | InputDelay<N>{θ}, with Ts: the plant's input arrives θ late |
 */
template<typename TNoise = detail::Unset, typename TDisturbance = detail::Unset, typename TDelay = detail::Unset>
struct KalmanConfig {
    TNoise       noise{};       ///< KalmanNoise{.process, .measurement}: the process and measurement noise (required)
    TDisturbance disturbance{}; ///< Where an integrating disturbance enters, estimated on [x; d]
    TDelay       delay{};       ///< InputDelay<N>{θ}: the plant's input arrives θ late (with Ts)
};

/**
 * @brief design::mpc's weights, constraints and options as named fields
 *
 * | Field | Holds | Left out |
 * | --- | --- | --- |
 * | weights | MPCWeights{.Qy, .Rdu, …}, per prediction step | Qy = I, RΔu = I |
 * | constraints | MPCConstraints{.u_lim, .du_lim, .y_lim, …} | unbounded |
 * | noise | KalmanNoise{.process, .measurement[, .disturbance]}: the MPC on a Kalman filter's estimate | the state measured |
 * | disturbance | a disturbance spec, with a noise: OFMPC, offset-free | no disturbance estimated |
 * | measured | MeasuredDisturbance{Bd}: a disturbance input measured each tick | none |
 */
template<
    typename TWeights = detail::Unset,
    typename TConstraints = detail::Unset,
    typename TNoise = detail::Unset,
    typename TDisturbance = detail::Unset,
    typename TMeasured = detail::Unset>
struct MPCConfig {
    TWeights     weights{};     ///< MPCWeights{…}: the cost on the controlled outputs and the moves (left out: Qy = I, RΔu = I)
    TConstraints constraints{}; ///< MPCConstraints{…}: the limits, the terminal ingredients and the QP budget (left out: unbounded)
    TNoise       noise{};       ///< KalmanNoise{…}: the MPC on a Kalman filter's estimate (left out: the state measured)
    TDisturbance disturbance{}; ///< Where an integrating disturbance enters, with a noise: OFMPC
    TMeasured    measured{};    ///< MeasuredDisturbance{Bd}: a disturbance input measured each tick
};

/** @brief design::steady_state_target's options as named fields: the disturbance the target holds the outputs under, and the plant's input delay. */
template<typename TDisturbance = detail::Unset, typename TDelay = detail::Unset>
struct TargetConfig {
    TDisturbance disturbance{}; ///< Where d enters: InputDisturbance{}, OutputDisturbance{}, PlantDisturbance{} or DisturbanceModel{...}
    TDelay       delay{};       ///< InputDelay<N>{θ}: the plant's input arrives θ late (with Ts)
};

namespace detail {

/** @brief Call @p f with no fields left. */
template<typename F>
[[nodiscard]] constexpr auto forward_set(const F& f) {
    return f();
}

/** @brief Call @p f with the fields of a config that were given, in order: an Unset one is dropped. */
template<typename F, typename A, typename... Rest>
[[nodiscard]] constexpr auto forward_set(const F& f, const A& first, const Rest&... rest) {
    if constexpr (is_unset<A>) {
        static_cast<void>(first);
        return forward_set(f, rest...);
    } else {
        return forward_set([&](const auto&... later) { return f(first, later...); }, rest...);
    }
}

/** @brief design::state_feedback(sys, [Ts,] weights, opts…) from an LQConfig: each field checked for its kind, the given ones passed on. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename TW, typename TL, typename TDm, typename TI, typename TN, typename TD, typename TT, typename TDl, typename... Period>
[[nodiscard]] constexpr auto state_feedback_from(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const LQConfig<TW, TL, TDm, TI, TN, TD, TT, TDl>& config, Period... Ts) {
    static_assert(!is_unset<TW>, "LQConfig::weights is the design's cost: .weights = LQRWeights{.Q = ..., .R = {...}}.");
    static_assert(damp::is_same_v<TW, LQRWeights<NX, NU, T>>, "LQConfig::weights is the plant's LQRWeights<NX, NU, T>: sized to its states and inputs, in its scalar.");
    static_assert(is_unset<TL> || is_lqr_limits<TL>::value, "LQConfig::limits is LQRLimits{.u_lim = {...}}.");
    static_assert(is_unset<TDm> || is_demand_limit<TDm>::value, "LQConfig::demand is DemandLimit{.state = s, .outer = {i, ...}, .lim = {...}}.");
    static_assert(is_unset<TI> || is_integral_option<TI>::value, "LQConfig::integral is IntegralAction{.Q = ...}.");
    static_assert(is_unset<TN> || is_kalman_noise<TN>::value, "LQConfig::noise is KalmanNoise{.process = ..., .measurement = ...}.");
    static_assert(is_unset<TD> || is_disturbance_model<TD>, "LQConfig::disturbance is InputDisturbance{}, OutputDisturbance{}, PlantDisturbance{} or DisturbanceModel{.Bd, .Cd}.");
    static_assert(is_unset<TT> || is_tracked_outputs<TT>, "LQConfig::tracked is TrackedOutputs{} or TrackedOutputs{{i, ...}}.");
    static_assert(is_unset<TDl> || is_input_delay<TDl>::value, "LQConfig::delay is InputDelay<N>{theta}.");
    return forward_set(
        [&](const auto&... options) { return design::state_feedback(sys, Ts..., config.weights, options...); },
        config.limits, config.demand, config.integral, config.noise, config.disturbance, config.tracked, config.delay
    );
}

/** @brief design::kalman(sys, [Ts,] opts…) from a KalmanConfig: each field checked for its kind, the given ones passed on. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename TN, typename TD, typename TDl, typename... Period>
[[nodiscard]] constexpr auto kalman_from(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const KalmanConfig<TN, TD, TDl>& config, Period... Ts) {
    static_assert(is_kalman_noise<TN>::value, "KalmanConfig::noise is the filter's noise: .noise = KalmanNoise{.process = ..., .measurement = ...}.");
    static_assert(is_unset<TD> || is_disturbance_model<TD>, "KalmanConfig::disturbance is InputDisturbance{}, OutputDisturbance{}, PlantDisturbance{} or DisturbanceModel{.Bd, .Cd}.");
    static_assert(is_unset<TDl> || is_input_delay<TDl>::value, "KalmanConfig::delay is InputDelay<N>{theta}.");
    return forward_set([&](const auto&... options) { return design::kalman(sys, Ts..., options...); }, config.noise, config.disturbance, config.delay);
}

/** @brief design::mpc(sys, [Ts,] weights, constraints, opts…) from an MPCConfig: the weights and the limits sized to the plant (their defaults left out), the other fields passed on. */
template<size_t NP, size_t NC, size_t NT_max, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename TW, typename TC, typename TN, typename TD, typename TM, typename... Period>
[[nodiscard]] constexpr auto mpc_from(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const MPCConfig<TW, TC, TN, TD, TM>& config, Period... Ts) {
    constexpr size_t NZ = mpc_outputs_in<NX, NU, NY, NW, TD>;
    static_assert(is_unset<TW> || damp::is_same_v<TW, MPCWeights<NU, NZ, T>>, "MPCConfig::weights is the plant's MPCWeights<NU, NZ, T>: one weight per input and controlled output, in its scalar.");
    static_assert(is_unset<TC> || is_mpc_constraints_on<TC, NU, NZ, T>, "MPCConfig::constraints is the plant's MPCConstraints<NU, NZ, T, NXC>: limits on its inputs and controlled outputs, in its scalar.");
    static_assert(is_unset<TN> || is_kalman_noise<TN>::value, "MPCConfig::noise is KalmanNoise{.process = ..., .measurement = ...}.");
    static_assert(is_unset<TD> || is_disturbance_model<TD>, "MPCConfig::disturbance is InputDisturbance{}, OutputDisturbance{}, PlantDisturbance{} or DisturbanceModel{.Bd, .Cd}.");
    static_assert(is_unset<TM> || is_measured_disturbance<TM>::value, "MPCConfig::measured is MeasuredDisturbance{Bd}.");
    const auto weights = [&] {
        if constexpr (is_unset<TW>) {
            return MPCWeights<NU, NZ, T>{};
        } else {
            return config.weights;
        }
    }();
    const auto constraints = [&] {
        if constexpr (is_unset<TC>) {
            return MPCConstraints<NU, NZ, T>{};
        } else {
            return config.constraints;
        }
    }();
    return forward_set(
        [&](const auto&... options) { return design::mpc<NP, NC, NT_max>(sys, Ts..., weights, constraints, options...); },
        config.noise, config.disturbance, config.measured
    );
}

/** @brief design::steady_state_target(sys, [Ts,] opts…) from a TargetConfig: each field checked for its kind, the given ones passed on. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename TD, typename TDl, typename... Period>
[[nodiscard]] constexpr auto target_from(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const TargetConfig<TD, TDl>& config, Period... Ts) {
    static_assert(is_unset<TD> || is_disturbance_model<TD>, "TargetConfig::disturbance is InputDisturbance{}, OutputDisturbance{}, PlantDisturbance{} or DisturbanceModel{.Bd, .Cd}.");
    static_assert(is_unset<TDl> || is_input_delay<TDl>::value, "TargetConfig::delay is InputDelay<N>{theta}.");
    return forward_set([&](const auto&... options) { return design::steady_state_target(sys, Ts..., options...); }, config.disturbance, config.delay);
}

} // namespace detail

// A config is a named object. Its first field's type has no default here, so braces in a config's
// place deduce none and stay the pack's: mpc<NP, NC>(sys, {}, {}) is the default weights and
// limits, not a period and a config.

/** @brief design::state_feedback(sys, weights, opts…) with the weights and the options the fields of @p config: the same design. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename C0, typename... C>
[[nodiscard]] constexpr auto state_feedback(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const LQConfig<C0, C...>& config) {
    return detail::state_feedback_from(sys, config);
}

/** @brief design::state_feedback(sys, Ts, weights, opts…), the sampled-data design, with the weights and the options the fields of @p config. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename C0, typename... C>
[[nodiscard]] constexpr auto state_feedback(const StateSpace<NX, NU, NY, T, NW, NV>& sys, damp::type_identity_t<T> Ts, const LQConfig<C0, C...>& config) {
    return detail::state_feedback_from(sys, config, Ts);
}

/** @brief design::kalman(sys, opts…) with the noise and the options the fields of @p config: the same design. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename C0, typename... C>
[[nodiscard]] constexpr auto kalman(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const KalmanConfig<C0, C...>& config) {
    return detail::kalman_from(sys, config);
}

/** @brief design::kalman(sys, Ts, opts…), the filter of a continuous plant sampled at @p Ts, with the noise and the options the fields of @p config. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename C0, typename... C>
[[nodiscard]] constexpr auto kalman(const StateSpace<NX, NU, NY, T, NW, NV>& sys, damp::type_identity_t<T> Ts, const KalmanConfig<C0, C...>& config) {
    return detail::kalman_from(sys, config, Ts);
}

/** @brief design::mpc<NP, NC>(sys, weights, constraints, opts…) with them the fields of @p config: the same design. */
template<size_t NP, size_t NC, size_t NT_max = 0, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename C0, typename... C>
[[nodiscard]] constexpr auto mpc(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const MPCConfig<C0, C...>& config) {
    return detail::mpc_from<NP, NC, NT_max>(sys, config);
}

/** @brief design::mpc<NP, NC>(plant, Ts, weights, constraints, opts…), on a continuous plant sampled at @p Ts, with them the fields of @p config. */
template<size_t NP, size_t NC, size_t NT_max = 0, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename C0, typename... C>
[[nodiscard]] constexpr auto mpc(const StateSpace<NX, NU, NY, T, NW, NV>& plant, damp::type_identity_t<T> Ts, const MPCConfig<C0, C...>& config) {
    return detail::mpc_from<NP, NC, NT_max>(plant, config, Ts);
}

/** @brief design::steady_state_target(sys, opts…) with the options the fields of @p config: the same design. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename C0, typename... C>
[[nodiscard]] constexpr auto steady_state_target(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const TargetConfig<C0, C...>& config) {
    return detail::target_from(sys, config);
}

/** @brief design::steady_state_target(sys, Ts, opts…), of a continuous plant sampled at @p Ts, with the options the fields of @p config. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename C0, typename... C>
[[nodiscard]] constexpr auto steady_state_target(const StateSpace<NX, NU, NY, T, NW, NV>& sys, damp::type_identity_t<T> Ts, const TargetConfig<C0, C...>& config) {
    return detail::target_from(sys, config, Ts);
}

} // namespace design
} // namespace damp
