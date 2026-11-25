// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file state_feedback_types.hpp
 * @brief The types of an LQ design: what it is given (LQRWeights, IntegralAction; InputLimits and
 *        InputRateLimits in design/input_limits.hpp) and what it hands its runtime
 *        (StateFeedbackResult, NoEstimator; the limits it takes at run time, StateFeedbackLimits, in
 *        design/input_limits.hpp)
 *
 * The design is design::state_feedback (state_feedback.hpp); the runtime is StateFeedback
 * (state_feedback_runtime.hpp), which LQR, LQI, LQG, LQGI and OFLQG (lqr.hpp, lqi.hpp, lqg.hpp,
 * lqgi.hpp, oflqg.hpp) alias, over detail/state_feedback_core.hpp.
 */

#include "damp/backend.hpp"
#include "damp/controllers/law.hpp"
#include "damp/design/detail/config.hpp"
#include "damp/design/input_limits.hpp" // IWYU pragma: export
#include "damp/design/options.hpp"
#include "damp/design/status.hpp"
#include "damp/design/steady_state_target.hpp"
#include "damp/estimation/fixed_gain_observer.hpp"
#include "damp/math/complex.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/toolbox/bounds.hpp"

namespace damp {
namespace design {

/** @brief The filter half of an estimated design (defined in kalman.hpp). */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t ND>
struct KalmanResult;

} // namespace design

/**
 * @brief The estimator of a law that measures its state: none, an empty type, so a
 *        `[[no_unique_address]]` member of it takes no storage.
 */
struct NoEstimator {
    /** @brief The estimator in another scalar: still none. */
    template<typename U>
    using as_t = NoEstimator;

    /** @brief Any scalar: still nothing. */
    template<typename U>
    [[nodiscard]] constexpr NoEstimator as() const {
        return {};
    }
};

namespace detail {

/** @brief What an estimator type carries: none (NoEstimator). */
template<typename E>
struct estimator_traits {
    static constexpr bool   present = false; ///< A filter runs
    static constexpr size_t NXE = 0;         ///< Filter states before the disturbance
    static constexpr size_t ND = 0;          ///< Disturbance states
    static constexpr size_t NUK = 0;         ///< Inputs it takes: u, then a measured disturbance's
};

/** @brief A FixedGainObserver's sizes. */
template<size_t NXK, size_t NU, size_t NY, typename T, size_t ND_, typename Store>
struct estimator_traits<FixedGainObserver<NXK, NU, NY, T, ND_, Store>> {
    static constexpr bool   present = true;  ///< A filter runs
    static constexpr size_t NXE = NXK - ND_; ///< Filter states before the disturbance
    static constexpr size_t ND = ND_;        ///< Disturbance states
    static constexpr size_t NUK = NU;        ///< Inputs it takes: u, then a measured disturbance's
};

/** @brief The runtime of a design's filter half: none for a measured state. */
template<typename Kalman, typename T>
struct runtime_estimator {
    using type = NoEstimator; ///< NoEstimator
};

/** @brief A KalmanResult runs as the FixedGainObserver of its sizes (its fixed gain M), in the runtime's scalar @p T. */
template<size_t NXK, size_t NU, size_t NY, typename U, size_t NW, size_t NV, size_t ND, typename T>
struct runtime_estimator<design::KalmanResult<NXK, NU, NY, U, NW, NV, ND>, T> {
    using type = FixedGainObserver<NXK, NU, NY, T, ND>; ///< The observer and the disturbance
};

/** @brief The runtime of a design's filter half (@p Kalman: KalmanResult or NoEstimator), in the scalar @p T. */
template<typename Kalman, typename T>
using runtime_estimator_t = runtime_estimator<Kalman, T>::type;

} // namespace detail

namespace design {

namespace detail {

/**
 * @brief Integral action as a design reads it: integrators on r − Ci·y, Ci an NI × NY matrix, with
 *        their cost and tracking time. design::IntegralAction{…} makes one from the outputs it
 *        lists; matlab::lqi(sys, Q, R, Ci) takes any Ci. NI = 0: none.
 */
template<size_t NI, size_t NY, typename T>
struct Integrators {
    Matrix<NI, NY, T> Ci{};            ///< The integrated outputs: e = r − Ci·y
    Matrix<NI, NI, T> Q{};             ///< The cost on ξ ≈ ∫e dt (positive definite)
    T                 tracking_time{}; ///< [s] Back-calculation tracking time (0: conditional integration)
};

} // namespace detail

/**
 * @struct LQRWeights
 * @brief The cost of an LQ design (design::state_feedback): what to penalize and how much
 *
 * The weights Q, R and N, in the plant's own time domain: a cost integral on a continuous plant
 * (with a Ts argument it is discretized exactly), per step on a discrete one. Braces size it from
 * the call: `{.Q = Q, .R = R}`. Named, it takes its sizes from its fields:
 * `LQRWeights{.Q = diag(1.0, 0.1), .R = {0.5}}` is LQRWeights<2, 1>. NU is 1 unless a field says
 * otherwise, so one input's weight is braced, `{0.5}`. Integral action brings its own cost:
 * IntegralAction{.Q, …} among the options. The limits the law runs within are not the cost's:
 * InputLimits among the options.
 *
 * | Field | Penalizes | Size |
 * | --- | --- | --- |
 * | Q | the plant states x | NX × NX, positive semidefinite |
 * | R | the inputs u | NU × NU, positive definite |
 * | N | the cross term 2xᵀNu (default zero) | NX × NU |
 *
 * An R that is not positive definite is refused with DesignStatus::bad_weights.
 */
template<size_t NX, size_t NU = 1, typename T = double>
struct LQRWeights {
    using option_kind = kind::Weights; ///< Its kind of design option
    Matrix<NX, NX, T> Q{};             ///< On the plant states x (positive semidefinite)
    Matrix<NU, NU, T> R{};             ///< On the inputs u (positive definite)
    Matrix<NX, NU, T> N{};             ///< Cross term between x and u (default zero)
};

namespace detail {

/** @brief The scalar of a pole as written: a floating-point number is its own. */
template<typename P>
struct pole_scalar {
    using type = P; ///< The number's type
};

/** @brief A damp::complex: its parts'. */
template<typename U>
struct pole_scalar<damp::complex<U>> {
    using type = U; ///< The parts' type
};

} // namespace detail

/**
 * @struct Poles
 * @brief The closed-loop poles of a placed design (design::state_feedback): the weights' slot,
 *        in place of LQRWeights
 *
 * `design::state_feedback(plant, Ts, Poles{-8.0, -12.0}, opts...)` places eig(A − B·K) where the
 * list says, by design::place's machinery (KNV, Ackermann for one input), and builds the same law
 * the LQ weights would, with every option: the result and its runtime are the LQ design's types.
 * The count is the plant's states plus the integrators of an IntegralAction, refused at compile
 * time otherwise; an InputDelay's delay states, which only hold inputs already sent, take theirs
 * at 0 (deadbeat), and a disturbance spec's states none, since the target cancels them.
 *
 * The list is in the plant's own time domain, as design::place takes it: rad/s (s-plane) on the
 * sampled form, mapped to z = e^{s·Ts}; z-plane on a discrete plant. Complex poles come as
 * adjacent conjugate pairs, `Poles{damp::complex{-1.0, 2.0}, damp::complex{-1.0, -2.0}}`, and a
 * real one converts, `Poles{-5.0, -5.0}`. The count and the scalar are deduced from the entries,
 * with bare braces (a deduction guide over the C array): the first entry's scalar, real or
 * complex, is the list's.
 *
 * @tparam N The poles listed
 * @tparam T Their scalar (the plant's, or converted to it)
 */
template<size_t N, typename T = double>
struct Poles {
    using option_kind = kind::Weights; ///< Its kind of design option
    static_assert(damp::is_floating_point_v<T>, "A pole is a floating-point number or a damp::complex of one: Poles{-5.0, -5.0}, Poles{damp::complex{-1.0, 2.0}, damp::complex{-1.0, -2.0}}.");

    damp::complex<T> poles[N]{}; ///< The closed-loop poles, adjacent conjugate pairs
};

/** @brief Poles{p, …}: the count from the list, the scalar from its first entry. */
template<typename P, typename... Rest>
Poles(P, Rest...) -> Poles<1 + sizeof...(Rest), typename detail::pole_scalar<P>::type>;

/**
 * @struct ObserverPoles
 * @brief The error poles of a placed filter (design::state_feedback): a trailing option, in place
 *        of KalmanNoise
 *
 * `design::state_feedback(plant, Ts, weights, ObserverPoles{-40.0, -40.0, -40.0}, InputDisturbance{})`
 * gives the law's filter a placed gain where a KalmanNoise gives it the steady-state Kalman gain:
 * the same FixedGainObserver in the current form, x̂[k|k] = x̂[k|k−1] + M·(y − C·x̂[k|k−1] − D·u), on
 * the same model, so the result and its runtime are the Kalman design's types. The poles are the
 * eigenvalues of its error dynamics, e[k|k] = (I − M·C)·A·e[k−1|k−1], the same as those of the
 * predictor error (A − L·C), L = A·M. Any weights take it, LQRWeights or Poles.
 *
 * The count is the plant's states plus the disturbance spec's (none without one), refused at
 * compile time otherwise. The delay states take theirs at 0: an InputDelay's hold inputs already
 * sent, which the readings do not correct, and a late reading's (OutputDelay) are placed at 0
 * with the rest. Units, pairs and deduction are Poles'.
 *
 * @tparam N The poles listed
 * @tparam T Their scalar (the plant's, or converted to it)
 */
template<size_t N, typename T = double>
struct ObserverPoles {
    using option_kind = kind::ObserverPoles; ///< Its kind of design option
    static_assert(damp::is_floating_point_v<T>, "A pole is a floating-point number or a damp::complex of one: ObserverPoles{-40.0, -40.0}.");

    damp::complex<T> poles[N]{}; ///< The filter's error poles, adjacent conjugate pairs
};

/** @brief ObserverPoles{p, …}: the count from the list, the scalar from its first entry. */
template<typename P, typename... Rest>
ObserverPoles(P, Rest...) -> ObserverPoles<1 + sizeof...(Rest), typename detail::pole_scalar<P>::type>;

/** @brief The spec of a law's gain: LQRWeights for the LQ gain, or Poles for a placed one. The weights are an argument, never an option. */
template<typename W>
concept GainSpec = OptionOf<W, kind::Weights>;

/** @brief The spec of a law's filter: a KalmanNoise for the Kalman gain, or ObserverPoles for a placed one. */
template<typename F>
concept FilterSpec = OptionOf<F, kind::KalmanNoise> || OptionOf<F, kind::ObserverPoles>;

namespace detail {

/** @brief A list of poles in the scalar @p T, as the placement takes them. */
template<typename T, size_t N, typename U>
[[nodiscard]] constexpr damp::array<damp::complex<T>, N> pole_list(const damp::complex<U> (&poles)[N]) {
    damp::array<damp::complex<T>, N> out{};
    for (size_t i = 0; i < N; ++i) {
        out[i] = damp::complex<T>(poles[i]);
    }
    return out;
}

/** @brief The poles of a placed gain or filter, in the plant's scalar: the gain's or the filter's, and the design adds 0 for the delay states. */
template<size_t N, typename T>
struct PlacedPoles {
    damp::array<damp::complex<T>, N> poles{}; ///< The poles listed, in the time domain of the plant
};

/** @brief Not a placed spec: the LQ gain, or the Kalman filter. */
template<typename G>
inline constexpr bool is_placed = false;

/** @brief Placed. */
template<size_t N, typename T>
inline constexpr bool is_placed<PlacedPoles<N, T>> = true;

/** @brief The gain of a design is the LQ one, from its weights. */
struct LQGain {};

/** @brief The demand limit among @p opts on the gain's terms (@p Terms, a DemandTerms), or off without one. */
template<typename Terms, typename... Opts>
[[nodiscard]] constexpr Terms demand_terms_of(const Opts&... opts) {
    if constexpr (option_count<kind::DemandLimit, Opts...> == 0) {
        return Terms{};
    } else {
        return Terms(the<kind::DemandLimit>(opts...));
    }
}

} // namespace detail

/**
 * @brief Integral action: integrators on r − y, with the cost and the limits that go with them.
 *        Among a design's trailing options it makes the law an LQI (an LQGI with a KalmanNoise).
 *
 * The cost on the integrators is the action's own, so the design's LQRWeights are the same with
 * or without it: `state_feedback(plant, Ts, {.Q = Q, .R = {1.0}}, IntegralAction{.Q = 80.0})`.
 * `.Q`'s type is its value's: a floating-point scalar for one integrator, an NI × NI Matrix for
 * NI; the fields it sizes take braces.
 *
 * Inside a template, clang (18 to 23 at least) deduces nothing from a `.Q` set to a literal or a
 * constant (llvm/llvm-project#114870): `IntegralAction{.Q = 80.0}` there is `IntegralAction<>`, and
 * the design does not compile. Name the type there, `IntegralAction<double>{.Q = 80.0}` or
 * `IntegralAction<Matrix<2, 2>>{…}`. The other options deduced from their fields share it
 * (docs/known_limitations.md).
 *
 * | Field | What it is | Size | Left out |
 * | --- | --- | --- | --- |
 * | Q | the cost on the discrete integral state ξ ≈ ∫(r − y) dt | NI × NI, positive definite | required |
 * | outputs | the integrated outputs y = C·x + D·u, in that order (r follows it) | NI | every output (NI = NY) |
 * | limits | limits on ξ, per integrator, a Bounds | NI | unbounded |
 * | tracking_time | the integrators' back-calculation tracking time T_t [s] | ≥ 0 | 0: conditional integration |
 *
 * Q is the cost on the Euler-integrated error ξ ≈ ∫e dt, not a block of a generic augmented Q.
 * LQRWeights' Q, R and N are the plant cost (Van Loan on a continuous plant); Q keeps that meaning
 * on both the continuous and discrete routes. The stacked DARE cost is block-diagonal in the two:
 * x–ξ cross terms are not a public weight. The closed loop still couples x and ξ, so raising Q
 * speeds the integral and moves the regulator poles too. A Q that is not positive definite is
 * refused with DesignStatus::bad_weights, and so are outputs out of range, listed twice, or left
 * out with fewer integrators than outputs.
 *
 * The law holds its integrators against the clamp of the InputLimits (and of an InputRateLimits),
 * with no back_calculate. limits is on ξ itself, in the units of ∫e dt (the error's times seconds):
 * through the gain, ξ moves input j by K(j, NX + i)·ξᵢ. The runtime moves it at run time with
 * set_limits (StateFeedbackLimits::integral), and refuses it as it refuses the InputLimits.
 *
 * tracking_time picks how the integrators meet a clamp, as PID's Kbc does. At 0, the default, an
 * integrator holds where integrating its error would push a clamped input further (conditional
 * integration): it stops, and keeps what it had. Above 0 it integrates as ever and is pulled toward
 * the value that would have asked for the input applied (back-calculation),
 * @f[
 *   \xi \leftarrow \xi + T_s\, e + \min\!\left(\tfrac{T_s}{T_t}, 1\right)
 *   \left(\tfrac{\partial u}{\partial \xi}\right)^{+} \left(u_{\text{applied}} - u_{\text{asked}}\right),
 * @f]
 * the pseudo-inverse spreading the excess over the integrators by least squares, T_t = Kbc/Ki on a
 * PID. A tracking time up to Ts takes the whole excess each tick: the integrators track the clamp.
 * A longer one takes a share, never more than the whole: past that it would overshoot the
 * realizable command, and past twice it diverges.
 *
 * The two meet a clamp differently. Back-calculation recovers from a long one without the
 * overshoot of holding: on a cart-pole with its force limited, a step that overshot by 55% while
 * the integrators held settles without overshoot at a tracking time of Ts to 0.2 s. It also
 * unwinds what the integrators carry against a constant load, integrated back once the clamp lets
 * go, and it takes back each excursion of noise that reaches the clamp, from that side only: an LQI
 * at 96% of its limit with noise on its reading settles with an offset about 13 times holding's,
 * and a slow integrator gives back for seconds what a short clamp took from it. Holding keeps what
 * the integrators carry, which is why it is the default, as Kbc = 0 is on a PID. It is a member of
 * the result and of the runtime, `law.tracking_time`, set at run time; one that is negative or not
 * finite is refused.
 *
 * @see "Advanced PID Control" (Åström & Hägglund, 2006), §3.5: back-calculation and conditional
 *      integration
 */
template<typename TQ = detail::Unset, typename T = detail::scalar_of_t<TQ>>
struct IntegralAction {
    using option_kind = kind::IntegralAction;               ///< Its kind of design option
    static constexpr size_t NI = detail::weight_size<TQ>(); ///< Integrators

    TQ                                             Q{};                                  ///< On ξ ≈ ∫(r − y) dt (NI × NI, positive definite)
    damp::type_identity_t<damp::array<size_t, NI>> outputs = detail::every_output<NI>(); ///< The integrated outputs, in order (left out: every one)
    damp::type_identity_t<Bounds<NI, T>>           limits{};                             ///< Limits on ξ, per integrator (default unbounded)
    damp::type_identity_t<T>                       tracking_time{};                      ///< [s] Back-calculation tracking time (default 0: conditional integration)
};

/**
 * @brief Bumpless handover for a law with nothing to preload: the first command after track() is
 *        the input tracked, and the law's own takes over with the time constant fade_time. Among
 *        design::state_feedback's trailing options.
 *
 * While another source drives the plant (an operator, another law), the caller calls the law's
 * track(u, …) in place of step(…), u the input applied. An LQI or LQGI preloads its integrators so
 * that its first command is u, and an OFLQG at a steady state has learned the load into d̂; an LQR
 * or LQG has nothing to preload, and an OFLQG not at a steady state has a transient. With this
 * option the first step after a track stores the mismatch e = u − u_law and adds it to the command,
 * e decaying by exp(−Ts/fade_time) each tick. The sum is clamped into the input limits, and the
 * integrators take e as a feedforward, not as saturation. An LQI or LQGI takes it too, for what the
 * preload leaves: the least-squares remainder with fewer integrators than inputs, or a ξ held at
 * its limit. Without it the first step is the law's own command, and the law pays nothing for it:
 * no storage, no per-tick work. With it the law holds NU values, the factor and a flag, and adds e
 * to each command. design::mpc refuses it: an offset outside the QP would break its constraints,
 * and the MPC's handover is its held input, the base of its first move.
 *
 * | Field | What it is | Left out |
 * | --- | --- | --- |
 * | fade_time | [s] the time constant e decays with, finite and > 0 | required |
 *
 * A fade_time that is not finite and positive is refused with DesignStatus::invalid_parameter.
 */
struct Handover {
    using option_kind = kind::Handover; ///< Its kind of design option
    double fade_time{0};                ///< [s] The time constant the mismatch decays with: exp(−Ts/fade_time) per tick
};

/** @brief The integral action of a design: the integrated outputs, their maps of the state and the input, the per-sample step, and what back-calculation moves the integrators by. */
template<size_t NI, size_t NX, size_t NU, size_t NY, typename T = double>
struct IntegralDesign {
    Matrix<NI, NY, T> Ci{};     ///< The integrated outputs: e = r − Ci·y, y = C·x + D·u
    Matrix<NI, NX, T> Ci_x{};   ///< Ci·C: the integrated outputs as a map of the plant state (what a law that measures its state integrates)
    Matrix<NI, NU, T> Ci_u{};   ///< Ci·D: the integrated outputs' feedthrough of the input applied over the period
    T                 step{1};  ///< Per-sample integrator gain: the design's Ts (xi += Ts·e ≈ ∫e dt)
    Matrix<NI, NU, T> unwind{}; ///< (∂u/∂ξ)⁺, the pseudo-inverse of −K on the integrators: the move of ξ that back-calculation makes per unit of a clamp's excess

    /** @brief The same integral action in the scalar U. */
    template<typename U>
    [[nodiscard]] constexpr IntegralDesign<NI, NX, NU, NY, U> as() const {
        return {Ci.template as<U>(), Ci_x.template as<U>(), Ci_u.template as<U>(), static_cast<U>(step), unwind.template as<U>()};
    }
};

/**
 * @struct StateFeedbackResult
 * @brief A state-feedback design: the gain, and the parts the law's options add
 *
 * Its template arguments are the runtime's: StateFeedback deduces them from it, and so does the name of
 * its kind, LQR, LQI, LQG, LQGI or OFLQG, which refuses another kind's. The
 * gain is on the law's state [x; xi]: the plant, then the integrators (NI > 0).
 *
 * | Member | What it is | Present |
 * | --- | --- | --- |
 * | K, S, e | the gain, the Riccati solution (zero for a gain placed at Poles), the closed-loop poles | always |
 * | state_gain(), delay_gain(), integral_gain() | K's columns on x, on the delay states (NRU > 0) and on ξ (NI > 0) | always, and as noted |
 * | success, status, Ts | ok or why not; the law's sample time | always |
 * | limits | the limits the law starts with, its limits(): input (InputLimits), input_rate (InputRateLimits, NDU > 0), integral (IntegralAction::limits, NI > 0), demand (DemandLimit), disturbance (ND > 0, unbounded: the runtime's) | always, default unbounded, the demand limit off |
 * | kalman | the Kalman filter (design::kalman's result), or the filter placed at ObserverPoles in its shape | lqg, lqgi, oflqg |
 * | integral | the integrated outputs Ci, Ci·C, Ci·D and the step | NI > 0 |
 * | tracking_time | the integrators' back-calculation tracking time (IntegralAction::tracking_time) | NI > 0, default 0 |
 * | target | the steady-state target (design::steady_state_target's result) | oflqg; lqr and lqg with TrackedOutputs (NZ > 0); any law with a MeasuredDisturbance (NM > 0) |
 * | fade_time | the handover's fade time (Handover::fade_time) | HO, a design with a Handover |
 *
 * With a MeasuredDisturbance (NM > 0) the target is the steady state under it: target_d is [x_t; u_t]
 * per unit [d̂; d], the filter's disturbance then the measured one, and on a law whose r is not the
 * target's (a regulator, NZ = 0, holding every output at 0; an LQI or LQGI, its integrated outputs)
 * it has no target_r. The filter takes d as a known input after u: its model's B is [B, B_d].
 *
 * A design on a plant whose input arrives late (an InputDelay) is on the law's model
 * [x; v(k−1); …; v(k−NRU)]: K, S, e, the integrated outputs, the filter's model and the target
 * cover the delay states. A measured-state runtime holds them, a filter carries them.
 *
 * @tparam NX     Plant states
 * @tparam NI     Integrators (0: none)
 * @tparam NY     Measured outputs (NX for a law without integrators or estimator)
 * @tparam Kalman The filter's design: KalmanResult, or NoEstimator for a measured state
 * @tparam NZ     Controlled outputs of the target (0: none)
 * @tparam NRU    Delay states per input, the commands still on their way: NRU·NU states after the plant's (0: none)
 * @tparam NM     Measured-disturbance inputs (MeasuredDisturbance; 0: none)
 * @tparam NDU    Inputs with a move limit: NU with an InputRateLimits (0: none)
 * @tparam HO     Whether the law fades the mismatch after a track (a Handover among the options; false: none)
 *
 * @see "Optimal Control" (Anderson & Moore, 1990), §4 and §8
 */
template<size_t NX, size_t NU, typename T = double, size_t NI = 0, size_t NY = NX, typename Kalman = NoEstimator, size_t NZ = 0, size_t NRU = 0, size_t NM = 0, size_t NDU = 0, bool HO = false>
struct [[nodiscard]] StateFeedbackResult {
    using Filter = ::damp::detail::estimator_traits<::damp::detail::runtime_estimator_t<Kalman, T>>;
    static constexpr bool   kEstimated = Filter::present;         ///< A Kalman filter estimates the state
    static constexpr bool   kIntegral = NI > 0;                   ///< The law has integral action
    static constexpr bool   kTarget = (NZ > 0) || (NM > 0);       ///< The law has a steady-state target: held outputs, or a measured disturbance
    static constexpr bool   kHandover = HO;                       ///< The law fades the mismatch after a track (Handover)
    static constexpr size_t NXL = NX + (NRU * NU);                ///< The law's model: the plant, then its delay states
    static constexpr size_t NK = NXL + NI;                        ///< The gain's state: the law's model and the integrators
    static constexpr size_t ND = Filter::ND;                      ///< Disturbance states of the filter
    static constexpr size_t NXF = kEstimated ? Filter::NXE : NXL; ///< The filter's model before d: the law's, then the delay states of late readings
    using Integral = ::damp::detail::maybe<(NI > 0), IntegralDesign<NI, NXL, NU, NY, T>, ::damp::detail::OffSlot::design_integral>;
    using Target = ::damp::detail::maybe<kTarget, SteadyStateTargetResult<NXF, NU, ND + NM, NZ, T>, ::damp::detail::OffSlot::design_target>;
    using TrackingTime = ::damp::detail::maybe<(NI > 0), T, ::damp::detail::OffSlot::design_tracking_time>;
    using FadeTime = ::damp::detail::maybe<HO, T, ::damp::detail::OffSlot::design_fade_time>;

    Matrix<NU, NK, T>                           K{};                                 ///< Optimal gain: u = −K·[x; xi]
    Matrix<NK, NK, T>                           S{};                                 ///< Riccati solution (positive semidefinite; zero for a gain placed at Poles)
    damp::array<damp::complex<T>, NK>           e{};                                 ///< Closed-loop poles of the law's model
    bool                                        success{false};                      ///< status == DesignStatus::ok
    DesignStatus                                status{DesignStatus::not_converged}; ///< ok, or the first failure: the gain's, the limits', the filter's, the target's
    T                                           Ts{0};                               ///< Sample time of the law: 0 continuous (e in the s-plane), > 0 discrete
    [[no_unique_address]] Kalman                kalman{};                            ///< The Kalman filter (state_feedback with a KalmanNoise: LQG, LQGI, OFLQG)
    [[no_unique_address]] Integral              integral{};                          ///< The integral action (NI > 0)
    [[no_unique_address]] Target                target{};                            ///< The steady-state target (NZ > 0, or NM > 0: under the measured disturbance)
    StateFeedbackLimits<NX, NU, T, NI, ND, NDU> limits{};                            ///< The limits the law starts with, its limits() (default unbounded, the demand limit off)
    [[no_unique_address]] TrackingTime          tracking_time{};                     ///< [s] The integrators' back-calculation tracking time (IntegralAction::tracking_time; NI > 0, default 0)
    [[no_unique_address]] FadeTime              fade_time{};                         ///< [s] The handover's fade time (Handover::fade_time; HO)

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief K's columns on the plant's states x (NU × NX): u = −state_gain()·x − … */
    [[nodiscard]] constexpr Matrix<NU, NX, T> state_gain() const { return K.template block<NU, NX>(0, 0).to_matrix(); }

    /**
     * @brief K's columns on the delay states v(k−1) … v(k−NRU) of a plant whose input arrives late
     *        (an InputDelay): the commands sent that have not reached it yet (NU × NRU·NU; NRU > 0).
     */
    [[nodiscard]] constexpr Matrix<NU, NRU * NU, T> delay_gain() const
        requires(NRU > 0)
    {
        return K.template block<NU, NRU * NU>(0, NX).to_matrix();
    }

    /** @brief K's columns on the integrators ξ (NU × NI; NI > 0). */
    [[nodiscard]] constexpr Matrix<NU, NI, T> integral_gain() const
        requires(NI > 0)
    {
        return K.template block<NU, NI>(0, NXL).to_matrix();
    }

    /** @brief The same design in the scalar @p U: `res.as<float>()` for the embedded runtime. */
    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        using Out = StateFeedbackResult<NX, NU, U, NI, NY, typename Kalman::template as_t<U>, NZ, NRU, NM, NDU, HO>;
        Out out{K.template as<U>(), S.template as<U>(), array_as<damp::complex<U>>(e), success, status, static_cast<U>(Ts), kalman.template as<U>()};
        if constexpr (NI > 0) {
            out.integral = integral.template as<U>();
            out.tracking_time = static_cast<U>(tracking_time);
        }
        if constexpr (kTarget) {
            out.target = target.template as<U>();
        }
        if constexpr (HO) {
            out.fade_time = static_cast<U>(fade_time);
        }
        out.limits = limits.template as<U>();
        return out;
    }

    /**
     * @brief Check if the closed loop is stable: every pole strictly inside the unit circle
     *        (Ts > 0), or strictly in the left half plane (Ts = 0). False for a failed design.
     */
    [[nodiscard]] constexpr bool is_stable() const {
        if (!success) {
            return false; // a failed design has no closed loop (e is zero)
        }
        for (size_t i = 0; i < NK; ++i) {
            if (Ts > T{0} ? (e[i].abs() >= T{1}) : (e[i].real() >= T{0})) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief The law as a discrete state-space block, for closing the loop in analysis
     *
     * The realization of the runtime's tick, by the law's features:
     *
     * | Law | State | Inputs | Output |
     * | --- | --- | --- | --- |
     * | measured state | none | x | u = −Kx |
     * | measured state, integrators | ξ | [r; x] | u |
     * | Kalman filter | x̂(k−1\|k−1) | y | u |
     * | Kalman filter, disturbance model (OFLQG) | [x̂; d̂](k−1\|k−1) | y | u |
     * | Kalman filter, integrators | [x̂(k\|k−1); ξ; u_prev] | [r; y] | u |
     *
     * With input delay states, x is the law's model [x; v(k−1); …], the state of the delayed plant the
     * loop closes on.
     *
     * With integrators, @f$ u = -K_x x - K_i \xi @f$ reads the current ξ, then
     * @f$ \xi \leftarrow \xi + h(r - C_i (C x + D u)) @f$ (h = integral.step), the measurement
     * less the feedthrough of the input it was taken under standing in for C x with a filter. With a filter it is the
     * current estimator: predict with the last command, the measurement update with the
     * innovation gain M (KalmanResult::M), then the law,
     * @f[
     *   A_c = (I-MC)(A-BK) + MDK,\quad B_c = M,\quad C_c = -K A_c,\quad D_c = -KM.
     * @f]
     * With a disturbance model the target folds into the gain: at r = 0 the target for d̂ is
     * x_t = T_x d̂, u_t = T_u d̂ (SteadyStateTargetResult::target_d), so the law
     * u = u_t − K(x̂ − x_t) is the same current estimator on [x; d] with the gain
     * [K, −(T_u + K T_x)] on [x̂; d̂]. The steady-state target is a feedforward on r outside this
     * block, with TrackedOutputs and with a disturbance model alike: the block is the regulator
     * (r = 0). It closes on any plant its y is read from, one with another delay among them:
     * u = C_c s + D_c y, its sign included. The block is the
     * law within its limits: the limits on the command, its move and the integrators, and the demand limit, are not part of it. A
     * measured disturbance (NM > 0) is a feedforward outside it as r is: the block is the law at
     * d = 0, its filter on the plant's inputs alone.
     */
    [[nodiscard]] constexpr auto to_ss() const {
        if constexpr (!kEstimated && NI == 0) {
            return StateSpace{.A = Matrix<0, 0, T>{}, .B = Matrix<0, NXL, T>{}, .C = Matrix<NU, 0, T>{}, .D = -K, .Ts = Ts};
        } else if constexpr (!kEstimated) {
            return integral_ss();
        } else if constexpr (NI == 0) {
            return estimator_ss();
        } else {
            return integral_estimator_ss();
        }
    }

private:
    /**
     * @brief The law with integrators on a measured state: u = C_c ξ + D_c [r; x], C_c = −K_i,
     *        D_c = [0, −K_x], and ξ ← ξ + h·(r − Ci·C·x − Ci·D·u): A_c = I − h·Ci·D·C_c,
     *        B_c = [hI, −h·Ci·C] − h·Ci·D·D_c.
     */
    [[nodiscard]] constexpr StateSpace<NI, NI + NXL, NU, T> integral_ss() const {
        const T                 h = integral.step;
        const Matrix<NU, NI, T> Cc = -K.template block<NU, NI>(0, NXL);
        Matrix<NU, NI + NXL, T> Dc{};
        Dc.template block<NU, NXL>(0, NI) = -K.template block<NU, NXL>(0, 0);
        Matrix<NI, NI + NXL, T> Bc{};
        Bc.template block<NI, NI>(0, 0) = Matrix<NI, NI, T>::identity() * h;
        Bc.template block<NI, NXL>(0, NI) = -integral.Ci_x * h;
        const Matrix<NI, NU, T> hCiD = integral.Ci_u * h;
        const Matrix<NI, NI, T> Ac = Matrix<NI, NI, T>::identity() - (hCiD * Cc);
        return {.A = Ac, .B = Matrix<NI, NI + NXL, T>(Bc - (hCiD * Dc)), .C = Cc, .D = Dc, .Ts = Ts};
    }

    /** @brief The gain on the filter's state: K, and with a disturbance model −(T_u + K T_x) on d̂, the target's at r = 0 folded in. */
    [[nodiscard]] constexpr Matrix<NU, NXL + ND, T> filter_gain() const {
        if constexpr (ND == 0) {
            return K;
        } else {
            Matrix<NU, NXL + ND, T> Ks{};
            Ks.template block<NU, NXL>(0, 0) = K;
            const Matrix<NXL, ND, T> Tx = target.target_d.template block<NXL, ND>(0, 0);
            const Matrix<NU, ND, T>  Tu = target.target_d.template block<NU, ND>(NXL, 0);
            Ks.template block<NU, ND>(0, NXL) = -(Tu + (K * Tx));
            return Ks;
        }
    }

    /** @brief The filter's model on the plant's inputs: itself, or without the measured disturbance's columns of B and D after them (NM > 0). */
    [[nodiscard]] constexpr decltype(auto) commanded_model() const {
        if constexpr (NM == 0) {
            return (kalman.sys);
        } else {
            constexpr size_t NF = NXL + ND;
            return StateSpace<NF, NU, NY, T>{
                .A = kalman.sys.A,
                .B = kalman.sys.B.template block<NF, NU>(0, 0).to_matrix(),
                .C = kalman.sys.C,
                .D = kalman.sys.D.template block<NY, NU>(0, 0).to_matrix(),
                .Ts = kalman.sys.Ts,
            };
        }
    }

    /** @brief The current-estimator law: the filter's estimate x̂(k−1|k−1) is the state ([x̂; d̂] with a disturbance model). */
    [[nodiscard]] constexpr StateSpace<NXL + ND, NY, NU, T> estimator_ss() const {
        constexpr size_t NF = NXL + ND;
        const auto&      sys = commanded_model();
        const auto       Kf = filter_gain();
        const auto       ImMC = Matrix<NF, NF, T>::identity() - (kalman.M * sys.C);
        const auto       Ac = (ImMC * (sys.A - (sys.B * Kf))) + ((kalman.M * sys.D) * Kf);
        return {.A = Ac, .B = kalman.M, .C = -(Kf * Ac), .D = -(Kf * kalman.M), .Ts = Ts};
    }

    /** @brief The current-estimator law with integrators: the state is [x̂(k|k−1); ξ; u_prev], the inputs [r; y] (to_ss()). */
    [[nodiscard]] constexpr StateSpace<NXL + NI + NU, NI + NY, NU, T> integral_estimator_ss() const {
        constexpr size_t         NS = NXL + NI + NU;
        const auto&              sys = commanded_model();
        const auto&              M = kalman.M;
        const Matrix<NU, NXL, T> Kx = K.template block<NU, NXL>(0, 0);
        const Matrix<NU, NI, T>  Ki = K.template block<NU, NI>(0, NXL);
        const auto               ImMC = Matrix<NXL, NXL, T>::identity() - (M * sys.C);
        const auto               MD = M * sys.D;
        const auto               KxImMC = Kx * ImMC;
        const auto               AmBKx = sys.A - (sys.B * Kx);
        const auto               KxM = Kx * M;
        const T                  h = integral.step;

        Matrix<NS, NS, T> Ac{};
        Ac.template block<NXL, NXL>(0, 0) = AmBKx * ImMC;
        Ac.template block<NXL, NI>(0, NXL) = -(sys.B * Ki);
        Ac.template block<NXL, NU>(0, NXL + NI) = -(AmBKx * MD);
        Ac.template block<NI, NI>(NXL, NXL) = Matrix<NI, NI, T>::identity();
        Ac.template block<NU, NXL>(NXL + NI, 0) = -KxImMC;
        Ac.template block<NU, NI>(NXL + NI, NXL) = -Ki;
        Ac.template block<NU, NU>(NXL + NI, NXL + NI) = Kx * MD;
        Matrix<NS, NI + NY, T> Bc{};
        Bc.template block<NI, NI>(NXL, 0) = Matrix<NI, NI, T>::identity() * h;
        Bc.template block<NXL, NY>(0, NI) = AmBKx * M;
        Bc.template block<NI, NY>(NXL, NI) = -integral.Ci * h;
        Bc.template block<NU, NY>(NXL + NI, NI) = -KxM;
        Matrix<NU, NS, T> Cc{};
        Cc.template block<NU, NXL>(0, 0) = -KxImMC;
        Cc.template block<NU, NI>(0, NXL) = -Ki;
        Cc.template block<NU, NU>(0, NXL + NI) = Kx * MD;
        Matrix<NU, NI + NY, T> Dc{};
        Dc.template block<NU, NY>(0, NI) = -KxM;
        //! ξ ← ξ + h·(r − Ci·y + Ci·D·u_prev − Ci·D·u): the measurement less the feedthrough of the
        //! input it was taken under, then the feedthrough of the input the law applies, u = C_c s + D_c [r; y].
        const Matrix<NI, NU, T> hCiD = integral.Ci_u * h;
        Ac.template block<NI, NU>(NXL, NXL + NI) = hCiD;
        const Matrix<NI, NS, T>      fed_a = hCiD * Cc;
        const Matrix<NI, NI + NY, T> fed_b = hCiD * Dc;
        for (size_t i = 0; i < NI; ++i) {
            for (size_t j = 0; j < NS; ++j) {
                Ac(NXL + i, j) -= fed_a(i, j);
            }
            for (size_t j = 0; j < NI + NY; ++j) {
                Bc(NXL + i, j) -= fed_b(i, j);
            }
        }
        return {.A = Ac, .B = Bc, .C = Cc, .D = Dc, .Ts = Ts};
    }
};

} // namespace design

namespace detail {

/** @brief Not a StateFeedbackResult. */
template<typename R>
struct is_state_feedback_result : damp::false_type {};

/** @brief A StateFeedbackResult. */
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, typename Kalman, size_t NZ, size_t NRU, size_t NM, size_t NDU, bool HO>
struct is_state_feedback_result<design::StateFeedbackResult<NX, NU, T, NI, NY, Kalman, NZ, NRU, NM, NDU, HO>> : damp::true_type {};

/**
 * @brief A StateFeedbackResult of any sizes: what StateFeedback's constructor takes, so a parameter
 *        added to the result does not reach the constructor's signature.
 */
template<typename R>
concept state_feedback_result = is_state_feedback_result<R>::value;

} // namespace detail

} // namespace damp
