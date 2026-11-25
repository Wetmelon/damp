// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file observer.hpp
 * @brief The placed filter (design::observer), design::kalman's counterpart: the same
 *        FixedGainObserver on the same model, its gain placed at a list of error poles
 *        (ObserverPoles), with the options of the filter designs
 *
 * design::place_observer places the bare gain on the plant it is given. design::observer builds the
 * model design::kalman would with the same options, the delay states of an InputDelay and an
 * OutputDelay and a disturbance spec's states included, and places the filter's error poles on it:
 *
 * @code
 * using namespace damp;
 * // continuous plant: s-plane poles in rad/s, one per plant state and per disturbance state
 * constexpr auto res = design::observer(plant, 0.01, {-30.0, -35.0, -40.0}, design::InputDisturbance{});
 * // discrete plant: z-plane poles
 * constexpr auto res_d = design::observer(plant_d, {0.6, 0.65});
 * static constinit Law<res, float> obs;   // the FixedGainObserver on the design's constants
 * // each tick: auto xhat = obs.step(y, u);
 * @endcode
 */

#include "damp/backend.hpp"
#include "damp/controllers/law.hpp"    // IWYU pragma: export
#include "damp/design/disturbance.hpp" // IWYU pragma: export
#include "damp/design/options.hpp"
#include "damp/estimation/detail/placed_filter.hpp"
#include "damp/estimation/fixed_gain_observer.hpp" // IWYU pragma: export
#include "damp/estimation/kalman.hpp"              // IWYU pragma: export
#include "damp/math/complex.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace design {

/**
 * @struct ObserverPoles
 * @brief The error poles of a placed filter: design::observer's poles spelled out, and a trailing
 *        option of design::state_feedback and design::mpc in place of KalmanNoise
 *
 * `design::observer(plant, Ts, {-30.0, -35.0, -40.0}, InputDisturbance{})` places a filter's error
 * poles; `design::state_feedback(plant, Ts, weights, ObserverPoles{-40.0, -40.0, -40.0},
 * InputDisturbance{})` gives the law's filter a placed gain where a KalmanNoise gives it the
 * steady-state Kalman gain. Either is the same FixedGainObserver in the current form,
 * x̂[k|k] = x̂[k|k−1] + M·(y − C·x̂[k|k−1] − D·u), on the same model, so the result and its runtime
 * are the Kalman design's types. The poles are the eigenvalues of its error dynamics,
 * e[k|k] = (I − M·C)·A·e[k−1|k−1], the same as those of the predictor error (A − L·C), L = A·M.
 * Any weights take it, LQRWeights or Poles.
 *
 * The count is the plant's states plus the disturbance spec's (none without one), refused at
 * compile time otherwise. The delay states take theirs at 0: an InputDelay's hold inputs already
 * sent, which the readings do not correct, and a late reading's (OutputDelay) are placed at 0
 * with the rest.
 *
 * The list is in the plant's own time domain: rad/s (s-plane) on the sampled form, mapped to
 * z = e^{s·Ts}; z-plane on a discrete plant. Complex poles come as adjacent conjugate pairs,
 * `ObserverPoles{damp::complex{-30.0, 20.0}, damp::complex{-30.0, -20.0}}`, and a real one
 * converts. The count and the scalar are deduced from the entries, with bare braces (a deduction
 * guide over the C array): the first entry's scalar, real or complex, is the list's.
 *
 * @tparam N The poles listed
 * @tparam T Their scalar (the plant's, or converted to it)
 */
template<size_t N, typename T = double>
struct ObserverPoles {
    using option_kind = kind::ObserverPoles; ///< Its kind of design option
    static constexpr size_t listed = N;      ///< The poles it lists
    static_assert(damp::is_floating_point_v<T>, "A pole is a floating-point number or a damp::complex of one: ObserverPoles{-40.0, -40.0}.");

    damp::complex<T> poles[N]{}; ///< The filter's error poles, adjacent conjugate pairs
};

/** @brief ObserverPoles{p, …}: the count from the list, the scalar from its first entry. */
template<typename P, typename... Rest>
ObserverPoles(P, Rest...) -> ObserverPoles<1 + sizeof...(Rest), typename detail::pole_scalar<P>::type>;

namespace detail {

/** @brief ObserverPoles: the filter's gain placed there, in the plant's scalar. */
template<typename T, size_t N, typename U>
[[nodiscard]] constexpr PlacedPoles<N, T> observer_spec_of(const ObserverPoles<N, U>& poles) {
    return {pole_list<T>(poles.poles)};
}

/**
 * @brief design::observer's poles as bare braces, sized from the plant and the options: the N it
 *        places, from the entries `{p, …}`, each real or complex, in any scalar.
 *
 * The parameter's type is never deduced from the argument, so the braces need no type and the
 * count is the design's: a list of another count does not compile (observer_pole_count), where an
 * aggregate of N entries would take a shorter list and fill the rest with poles at 0. The
 * constructor is implicit, as the braces' conversion.
 */
template<size_t N, typename T>
struct PoleList {
    damp::array<damp::complex<T>, N> poles{}; ///< The poles listed, in the plant's scalar

    /** @brief From the entries of bare braces, each a real or complex number. */
    template<typename... P>
        requires(damp::is_constructible_v<damp::complex<T>, const P&> && ...)
    constexpr PoleList(const P&... listed) { // NOLINT(google-explicit-constructor): the braces of a pole list
        if constexpr (observer_pole_count<N, sizeof...(P)>()) {
            size_t i = 0;
            ((poles[i++] = damp::complex<T>(listed)), ...);
        }
    }
};

/** @brief The poles design::observer takes as braces on this plant with the options @p Opts: the plant's states and the disturbance's, never deduced. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, typename... Opts>
using poles_for = damp::type_identity_t<PoleList<NX + disturbance_states_in<NX, NU, NY, NW, Opts...>, T>>;

/** @brief A brace list's poles as the placement takes them. */
template<size_t N, typename T>
[[nodiscard]] constexpr PlacedPoles<N, T> placed_of(const PoleList<N, T>& listed) {
    return {listed.poles};
}

/** @brief The options design::observer takes: a disturbance spec, a MeasuredDisturbance, an InputDelay, an OutputDelay (a filter spec only for its message). */
using observer_options = option_set<kind::InputDelay, kind::OutputDelay, kind::Disturbance, kind::MeasuredDisturbance, kind::KalmanNoise, kind::ObserverPoles>;

/**
 * @brief design::observer with its trailing options, in the plant's own domain (@p D plant) or at
 *        @p Ts: the filter placed at @p placed on the Kalman design's model (placed_filter_design),
 *        then a measured disturbance's known inputs after u (measured_filter).
 */
template<Domain D, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NP, typename... Opts>
[[nodiscard]] constexpr auto observer_with(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T Ts, const PlacedPoles<NP, T>& placed, const Opts&... opts) {
    observer_options::check<Opts...>();
    constexpr bool kNoNoise = option_count<kind::KalmanNoise, Opts...> == 0;
    constexpr bool kNoPoles = option_count<kind::ObserverPoles, Opts...> == 0;
    static_assert(kNoNoise, "design::observer places its filter at the poles it is given: the Kalman filter of a KalmanNoise is design::kalman(plant, [Ts,] noise, opts...).");
    static_assert(kNoPoles, "The poles are design::observer's argument after the plant (and Ts), observer(plant, [Ts,] {p, ...}, opts...), not an option.");
    constexpr bool   kModelFits = check_delays<D, Opts...>() && check_disturbance<NX, NY, T, Opts...>();
    constexpr bool   kOptionsFit = observer_options::known<Opts...> && observer_options::unique<Opts...> && kNoNoise && kNoPoles && kModelFits;
    constexpr size_t NKF = NX + disturbance_states_in<NX, NU, NY, NW, Opts...>; // the poles listed: the plant's and the disturbance's
    const auto       delay = delays_of(opts...);
    const auto       measured = measured_of<NX, T>(opts...);
    if constexpr (!kOptionsFit || NP != NKF) {
        if constexpr (kOptionsFit) {
            static_cast<void>(observer_pole_count<NKF, NP>()); // ObserverPoles{…} of another count (bare braces are sized: poles_for)
        }
        static_cast<void>(sys);
        static_cast<void>(Ts);
        static_cast<void>(placed);
        static_cast<void>(delay);
        static_cast<void>(measured);
        return KalmanResult<NX, NU, NY, T>{}; // the message says why: nothing further
    } else {
        const auto kf = placed_filter_design<D>(noise_model(sys), Ts, placed, option<kind::Disturbance>(NoDisturbance{}, opts...), delay);
        return measured_filter<D>(kf, sys, Ts, measured, delay);
    }
}

} // namespace detail

/**
 * @brief The placed filter on a discrete plant: the FixedGainObserver design::kalman would run on
 *        the same model, its error poles eig((I − M·C)·A) at @p poles (z-plane)
 *
 * design::kalman's counterpart, with its options and its result: where the Kalman design solves
 * the filter DARE for the gain, this one places it. The filter is the current estimator,
 * x̂[k|k] = x̂[k|k−1] + M·(y − C·x̂[k|k−1] − D·u), its error after an update
 * e[k|k] = (I − M·C)·A·e[k−1|k−1], so the poles are those of (I − M·C)·A, the same as those of the
 * predictor error A − L·C with L = A·M. Placement is design::place_observer's: place on
 * (Aᵀ, Cᵀ) by Kautsky–Nichols–Van Dooren, then A M = L, or M placed directly on (A, C·A) on a
 * numerically singular A. Without options it is design::place_observer(sys, poles)'s M and L in
 * the Kalman design's type.
 *
 * @p poles lists one pole per plant state and per state of the disturbance spec, its type sized
 * from the plant and the options (as design::kalman's noise is), so bare braces take it,
 * `observer(sys, {0.6, 0.65})`, and so does ObserverPoles{…} in any scalar. Real entries convert,
 * complex ones come as adjacent conjugate pairs; another count does not compile.
 *
 * @p opts trail in any order, at most one of each kind:
 * - A disturbance spec (InputDisturbance, OutputDisturbance, PlantDisturbance, DisturbanceModel):
 *   the filter is on [x; d], d after the plant's states, as design::kalman's with the same spec;
 *   the poles listed are the plant's then d's. ND ≤ NY, or it does not compile; an augment the
 *   readings do not observe is refused with not_placeable.
 * - MeasuredDisturbance{.Bd}: a disturbance read each tick, a known input after u: the model's B is
 *   [B, B_d] and D is [D, 0], the gain the same without it. Predict with [u; d].
 *
 * The result is a KalmanResult on the filter's model, M and L = A·M placed, Q, R and P zero (no
 * noise), so `FixedGainObserver obs{res}` and `Law<kRes, T>` deploy it as they deploy
 * design::kalman's. A delay among this form's options is refused at compile time: a delay is
 * sampled with the continuous plant (observer(plant, Ts, poles, opts…)).
 *
 * Status: invalid_sample_time for a plant that is not discrete, bad_weights for a non-finite
 * model entry or measured channel, invalid_parameter for a pole that is not finite, then
 * design::place_observer's: not_placeable for an unobservable mode or a spectrum that cannot be
 * placed (a pole repeated more than NY times, complex poles not in adjacent conjugate pairs),
 * singular for a rank-deficient C.
 *
 * @note Compare with MATLAB®'s L = place(A', C', p)' and the current estimator's Mx = A \ L, on the
 *       model augmented by hand.
 * @see "Digital Control of Dynamic Systems" (Franklin, Powell & Workman, 1998), §8.2 (current estimator)
 * @see "Disturbance Models for Offset-Free Model-Predictive Control" (Pannocchia & Rawlings, 2003), AIChE Journal 49(2)
 * @see design::place_observer for the bare gain on the plant; design::kalman for the Kalman gain on the same model
 *
 * @param sys   Discrete plant (Ts > 0)
 * @param poles The error poles, z-plane: the plant's, then the disturbance's
 * @param opts  A disturbance spec and a MeasuredDisturbance, each or none
 * @return KalmanResult on the filter's model, with M, L and status
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename... Opts>
[[nodiscard]] constexpr auto observer(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const detail::poles_for<NX, NU, NY, T, NW, Opts...>& poles, const Opts&... opts) {
    return detail::observer_with<detail::Domain::plant>(sys, T{0}, detail::placed_of(poles), opts...);
}

/** @brief The placed filter on a discrete plant, its poles spelled ObserverPoles{…} in any scalar: as observer(sys, {p, …}, opts…). */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NP, typename U, typename... Opts>
[[nodiscard]] constexpr auto observer(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const ObserverPoles<NP, U>& poles, const Opts&... opts) {
    return detail::observer_with<detail::Domain::plant>(sys, T{0}, detail::observer_spec_of<T>(poles), opts...);
}

/**
 * @brief The placed filter for a continuous plant sampled at @p Ts: s-plane @p poles in rad/s,
 *        mapped to z = e^{s·Ts} on the model design::kalman(plant, Ts, …) would be on
 *
 * The plant is discretized by ZOH at @p Ts, with its delays when the options hold them, and the
 * filter placed on that model as observer(sys, poles, opts…) places it. The same poles in rad/s
 * give the same error dynamics at any rate. A discrete @p sys, or a Ts that is not positive and
 * finite, is refused with invalid_sample_time.
 *
 * @p opts trail in any order, at most one of each kind:
 * - A disturbance spec, as on the discrete form: InputDisturbance adds to u, so with an InputDelay
 *   it reaches the plant through the delay states; another spec's d is a state of the plant held
 *   over the period (DisturbanceModel's B_d continuous, ZOH'd with the plant).
 * - MeasuredDisturbance{.Bd}: B_d continuous, ZOH'd with the plant, none on the input's delay
 *   states, and on a late reading's what the plant moved under it until the reading was taken.
 * - InputDelay<N>{θ}: the filter is on discretize(plant, Ts, InputDelay<N>{θ}), the plant then
 *   the delay states v(k−1) … v(k−N) of the inputs still on their way. M has no rows on them
 *   (predict shifts the applied input in, the readings do not correct it), and their error poles
 *   are at 0; with an InputDisturbance they carry d's correction.
 * - OutputDelay<M>{θ}: each reading arrives θ late (one θ, or a ColVec of one per output); the
 *   filter carries the readings already taken in M delay states per output, their error poles at 0
 *   with the rest, and update(y) takes the reading as it arrives. With several readings a pole is
 *   placed at most NY times, so two or more delay states per reading are refused with
 *   not_placeable.
 *
 * @code
 * constexpr auto obs = design::observer(plant, 0.01, {-30.0, -35.0}, InputDelay<3>{0.022});
 * constexpr auto late = design::observer(plant, 0.01, {-30.0, -35.0, -40.0}, design::OutputDisturbance{}, OutputDelay{0.004});
 * @endcode
 *
 * @see "Digital Control of Dynamic Systems" (Franklin, Powell & Workman, 1998), §8.2 (current estimator)
 *
 * @param sys   Continuous plant (sys.Ts = 0)
 * @param Ts    [s] Sampling time of the filter
 * @param poles The error poles, s-plane in rad/s: the plant's, then the disturbance's
 * @param opts  A disturbance spec, a MeasuredDisturbance, an InputDelay and an OutputDelay, each or none
 * @return KalmanResult on the sampled model, its delay states and d, with M, L and status
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename... Opts>
[[nodiscard]] constexpr auto observer(
    const StateSpace<NX, NU, NY, T, NW, NV>&             sys,
    damp::type_identity_t<T>                             Ts,
    const detail::poles_for<NX, NU, NY, T, NW, Opts...>& poles,
    const Opts&... opts
) {
    return detail::observer_with<detail::Domain::sampled>(sys, Ts, detail::placed_of(poles), opts...);
}

/** @brief The placed filter for a continuous plant sampled at @p Ts, its poles spelled ObserverPoles{…} in any scalar: as observer(sys, Ts, {p, …}, opts…). */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NP, typename U, typename... Opts>
[[nodiscard]] constexpr auto observer(const StateSpace<NX, NU, NY, T, NW, NV>& sys, damp::type_identity_t<T> Ts, const ObserverPoles<NP, U>& poles, const Opts&... opts) {
    return detail::observer_with<detail::Domain::sampled>(sys, Ts, detail::observer_spec_of<T>(poles), opts...);
}

} // namespace design
} // namespace damp
