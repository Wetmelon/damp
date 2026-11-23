// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file state_feedback_types.hpp
 * @brief The types of an LQ design: what it is given (LQRWeights, LQRLimits, IntegralAction) and what it hands
 *        its runtime (StateFeedbackResult, NoEstimator)
 *
 * The design is design::state_feedback (state_feedback.hpp); the runtimes are in lqr.hpp, lqi.hpp,
 * lqg.hpp, lqgi.hpp and oflqg.hpp, over the one they share in detail/state_feedback_core.hpp.
 */

#include "damp/backend.hpp"
#include "damp/controllers/law_detail.hpp"
#include "damp/design/detail/config.hpp"
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
};

/** @brief A FixedGainObserver's sizes. */
template<size_t NXK, size_t NU, size_t NY, typename T, size_t ND_, typename Store>
struct estimator_traits<FixedGainObserver<NXK, NU, NY, T, ND_, Store>> {
    static constexpr bool   present = true;  ///< A filter runs
    static constexpr size_t NXE = NXK - ND_; ///< Filter states before the disturbance
    static constexpr size_t ND = ND_;        ///< Disturbance states
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
 * @brief Integral action as a design reads it: integrators on r − y, every measured output
 *        (`IntegralAction{}`) or the listed ones (`IntegralAction<2>{{0, 3}}`, in that order:
 *        the references r follow it). design::IntegralAction{.Q, …} makes it.
 * @tparam NI Outputs integrated, listed in @ref outputs (0: every measured output)
 */
template<size_t NI = 0>
struct IntegralAction {
    damp::array<size_t, NI> outputs{}; ///< Indices of the integrated outputs, each once (none listed: all)
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
 * LQRLimits among the options.
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
    Matrix<NX, NX, T> Q{}; ///< On the plant states x (positive semidefinite)
    Matrix<NU, NU, T> R{}; ///< On the inputs u (positive definite)
    Matrix<NX, NU, T> N{}; ///< Cross term between x and u (default zero)
};

/**
 * @struct LQRLimits
 * @brief The limit an LQ law's command runs within: among design::state_feedback's options
 *
 * The runtime's, not the cost's: the gain is the same with or without it, so it is an option and
 * not part of LQRWeights. Named, it takes its sizes from its field, one input unless u_lim says
 * otherwise (`.u_lim = Bounds<2>{…}`, or `LQRLimits<2>{…}`).
 *
 * `LQRLimits{.u_lim = {.lower = -12.0, .upper = 12.0}}` clamps the command control returns to ±12,
 * per input (default unbounded); back_calculate stays for a stage that clamps after the law. u_lim
 * is a member of the runtime, `law.u_lim`, so a rail that moves is set at run time. A limit that is
 * NaN, or a lower limit above its upper one, is refused with DesignStatus::invalid_parameter.
 *
 * The other limits are options of their own: a DemandLimit on what the outer terms ask of one state,
 * and the IntegralAction's on the integrators (.lim, .tracking_time).
 */
template<size_t NU = 1, typename T = double>
struct LQRLimits {
    Bounds<NU, T> u_lim{}; ///< Limits on the command u, per input (default unbounded)
};

namespace detail {

/** @brief An LQRLimits: a kind of state_feedback option. */
template<typename O>
struct is_lqr_limits : damp::false_type {};

/** @brief LQRLimits of any size. */
template<size_t NU, typename T>
struct is_lqr_limits<LQRLimits<NU, T>> : damp::true_type {};

/** @brief A DemandLimit: a kind of state_feedback option. */
template<typename O>
struct is_demand_limit : damp::false_type {};

/** @brief DemandLimit of any length. */
template<size_t M, typename T>
struct is_demand_limit<DemandLimit<M, T>> : damp::true_type {};

/** @brief The inputs an LQRLimits' u_lim is on. */
template<typename L>
inline constexpr size_t limits_inputs = 0;

/** @brief LQRLimits<NU, …>: NU. */
template<size_t NU, typename T>
inline constexpr size_t limits_inputs<LQRLimits<NU, T>> = NU;

/** @brief The demand limit among @p opts on the gain's terms (@p Terms, a DemandTerms), or off without one. */
template<typename Terms, typename... Opts>
[[nodiscard]] constexpr Terms demand_terms_of(const Opts&... opts) {
    if constexpr (::damp::detail::option_count<is_demand_limit, Opts...> == 0) {
        return Terms{};
    } else {
        return Terms(::damp::detail::option<is_demand_limit>(Unset{}, opts...));
    }
}

/**
 * @brief The full cost a design reads: LQRWeights and LQRLimits as one, sized to the plant
 *
 * Its u_lim comes from the LQRLimits among the options and demand_lim from the DemandLimit
 * (unbounded and off without them); the integrators add theirs in LQIConfig.
 */
template<size_t NX, size_t NU = 1, typename T = double>
struct LQRConfig {
    Matrix<NX, NX, T>                  Q{};          ///< On the plant states x (positive semidefinite)
    Matrix<NU, NU, T>                  R{};          ///< On the inputs u (positive definite)
    Matrix<NX, NU, T>                  N{};          ///< Cross term between x and u (default zero)
    Bounds<NU, T>                      u_lim{};      ///< Limits on the command u, per input (default unbounded)
    ::damp::detail::DemandTerms<NX, T> demand_lim{}; ///< A limit on what the outer states ask of one state (default off)
};

/** @brief The weights, the LQRLimits and the DemandLimit among @p opts (none: unbounded, off) as the design's full cost, the demand on the plant states. */
template<size_t NX, size_t NU, typename T, typename... Opts>
[[nodiscard]] constexpr LQRConfig<NX, NU, T> lq_cost_of(const LQRWeights<NX, NU, T>& weights, const Opts&... opts) {
    const auto limits = ::damp::detail::option<is_lqr_limits>(LQRLimits<NU, T>{}, opts...);
    static_assert(limits_inputs<damp::remove_cvref_t<decltype(limits)>> == NU, "LQRLimits::u_lim is one Bounds per input of the plant: LQRLimits<NU>{...}, or .u_lim = Bounds<NU>{...}.");
    return {.Q = weights.Q, .R = weights.R, .N = weights.N, .u_lim = limits.u_lim.template as<T>(), .demand_lim = demand_terms_of<::damp::detail::DemandTerms<NX, T>>(opts...)};
}

} // namespace detail

namespace detail {

/** @brief The full cost and an IntegralAction sized to the plant, NI integrators: what a design reads them as. */
template<size_t NX, size_t NU = 1, size_t NI = 1, typename T = double>
struct LQIConfig {
    Matrix<NX, NX, T>                       Q{};             ///< On the plant states x (positive semidefinite)
    Matrix<NU, NU, T>                       R{};             ///< On the inputs u (positive definite)
    Matrix<NX, NU, T>                       N{};             ///< Cross term between x and u (default zero)
    Matrix<NI, NI, T>                       Q_integral{};    ///< On the discrete integral state ξ ≈ ∫e dt (positive definite)
    Bounds<NU, T>                           u_lim{};         ///< Limits on the command u, per input (default unbounded)
    Bounds<NI, T>                           integral_lim{};  ///< Limits on the integral state ξ, per integrator (default unbounded)
    ::damp::detail::DemandTerms<NX + NI, T> demand_lim{};    ///< A limit on what the outer states and integrators ask of one state (default off)
    T                                       tracking_time{}; ///< [s] The integrators' back-calculation tracking time (default 0: conditional integration)
};

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
 * | Field | What it is | Size | Left out |
 * | --- | --- | --- | --- |
 * | Q | the cost on the discrete integral state ξ ≈ ∫(r − y) dt | NI × NI, positive definite | required |
 * | outputs | the integrated outputs y = C·x + D·u, in that order (r follows it) | NI | every output (NI = NY) |
 * | lim | limits on ξ, per integrator | NI | unbounded |
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
 * The law holds its integrators against the clamp of LQRLimits::u_lim, with no back_calculate. lim is on ξ
 * itself, in the units of ∫e dt (the error's times seconds): through the gain, ξ moves input j by
 * K(j, NX + i)·ξᵢ. It is a member of the result and of the runtime, `law.integral_lim`, set at run
 * time like `law.u_lim`, and refused as u_lim is.
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
    static constexpr size_t NI = detail::size_of<TQ>; ///< Integrators

    TQ                                             Q{};                                  ///< On ξ ≈ ∫(r − y) dt (NI × NI, positive definite)
    damp::type_identity_t<damp::array<size_t, NI>> outputs = detail::every_output<NI>(); ///< The integrated outputs, in order (left out: every one)
    damp::type_identity_t<Bounds<NI, T>>           lim{};                                ///< Limits on ξ, per integrator (default unbounded)
    damp::type_identity_t<T>                       tracking_time{};                      ///< [s] Back-calculation tracking time (default 0: conditional integration)
};

namespace detail {

/** @brief The config of a law with @p NI integrators: the full cost (LQRConfig) for none, LQIConfig with them. */
template<size_t NX, size_t NU, size_t NI, typename T>
using lq_config = damp::conditional_t<(NI == 0), LQRConfig<NX, NU, T>, LQIConfig<NX, NU, NI, T>>;

} // namespace detail

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
 * Its template arguments are the runtime's: LQR, LQI, LQG, LQGI or OFLQG, the one its design names, deduces them from it. The
 * gain is on the law's state [x; xi]: the plant, then the integrators (NI > 0).
 *
 * | Member | What it is | Present |
 * | --- | --- | --- |
 * | K, S, e | the gain, the Riccati solution, the closed-loop poles | always |
 * | state_gain(), delay_gain(), integral_gain() | K's columns on x, on the delay states (NRU > 0) and on ξ (NI > 0) | always, and as noted |
 * | success, status, Ts | ok or why not; the law's sample time | always |
 * | u_lim | the limits on the command (LQRLimits::u_lim) | always, default unbounded |
 * | demand_lim | the limit on what the outer terms ask of one state (the design's DemandLimit) | always, default off |
 * | kalman | the Kalman filter (design::kalman's result) | lqg, lqgi, oflqg |
 * | integral | the integrated outputs Ci, Ci·C, Ci·D and the step | NI > 0 |
 * | integral_lim | the limits on ξ (IntegralAction::lim) | NI > 0, default unbounded |
 * | tracking_time | the integrators' back-calculation tracking time (IntegralAction::tracking_time) | NI > 0, default 0 |
 * | target | the steady-state target (design::steady_state_target's result) | oflqg; lqr and lqg with TrackedOutputs (NZ > 0) |
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
 *
 * @see "Optimal Control" (Anderson & Moore, 1990), §4 and §8
 */
template<size_t NX, size_t NU, typename T = double, size_t NI = 0, size_t NY = NX, typename Kalman = NoEstimator, size_t NZ = 0, size_t NRU = 0>
struct [[nodiscard]] StateFeedbackResult {
    using Filter = ::damp::detail::estimator_traits<::damp::detail::runtime_estimator_t<Kalman, T>>;
    static constexpr bool   kEstimated = Filter::present; ///< A Kalman filter estimates the state
    static constexpr bool   kIntegral = NI > 0;           ///< The law has integral action
    static constexpr size_t NXL = NX + (NRU * NU);        ///< The law's model: the plant, then its delay states
    static constexpr size_t NK = NXL + NI;                ///< The gain's state: the law's model and the integrators
    static constexpr size_t ND = Filter::ND;              ///< Disturbance states of the filter
    using Integral = ::damp::detail::maybe<(NI > 0), IntegralDesign<NI, NXL, NU, NY, T>, ::damp::detail::OffSlot::design_integral>;
    using Target = ::damp::detail::maybe<(NZ > 0), SteadyStateTargetResult<NXL, NU, ND, NZ, T>, ::damp::detail::OffSlot::design_target>;
    using IntegralLim = ::damp::detail::maybe<(NI > 0), Bounds<NI, T>, ::damp::detail::OffSlot::design_integral_lim>;
    using TrackingTime = ::damp::detail::maybe<(NI > 0), T, ::damp::detail::OffSlot::design_tracking_time>;

    Matrix<NU, NK, T>                       K{};                                 ///< Optimal gain: u = −K·[x; xi]
    Matrix<NK, NK, T>                       S{};                                 ///< Riccati solution (positive semidefinite)
    damp::array<damp::complex<T>, NK>       e{};                                 ///< Closed-loop poles of the law's model
    bool                                    success{false};                      ///< status == DesignStatus::ok
    DesignStatus                            status{DesignStatus::not_converged}; ///< ok, or the first failure: the gain's, the filter's, the target's
    T                                       Ts{0};                               ///< Sample time of the law: 0 continuous (e in the s-plane), > 0 discrete
    [[no_unique_address]] Kalman            kalman{};                            ///< The Kalman filter (state_feedback with a KalmanNoise: LQG, LQGI, OFLQG)
    [[no_unique_address]] Integral          integral{};                          ///< The integral action (NI > 0)
    [[no_unique_address]] Target            target{};                            ///< The steady-state target (NZ > 0)
    Bounds<NU, T>                           u_lim{};                             ///< Limits on the command the runtime returns (default unbounded)
    ::damp::detail::DemandTerms<NX + NI, T> demand_lim{};                        ///< A limit on what the outer terms ask of one state (the design's DemandLimit; default off)
    [[no_unique_address]] IntegralLim       integral_lim{};                      ///< Limits on the integral state ξ (IntegralAction::lim; NI > 0, default unbounded)
    [[no_unique_address]] TrackingTime      tracking_time{};                     ///< [s] The integrators' back-calculation tracking time (IntegralAction::tracking_time; NI > 0, default 0)

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
        using Out = StateFeedbackResult<NX, NU, U, NI, NY, damp::remove_cvref_t<decltype(kalman.template as<U>())>, NZ, NRU>;
        Out out{K.template as<U>(), S.template as<U>(), array_as<damp::complex<U>>(e), success, status, static_cast<U>(Ts), kalman.template as<U>()};
        if constexpr (NI > 0) {
            out.integral = integral.template as<U>();
            out.integral_lim = integral_lim.template as<U>();
            out.tracking_time = static_cast<U>(tracking_time);
        }
        if constexpr (NZ > 0) {
            out.target = target.template as<U>();
        }
        out.u_lim = u_lim.template as<U>();
        out.demand_lim = demand_lim.template as<U>();
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
     * With input registers, x is the law's model [x; v(k−1); …], the state of the delayed plant the
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
     * law within its limits: u_lim, the integrators' limits and demand_lim are not part of it.
     */
    [[nodiscard]] constexpr auto to_ss() const {
        if constexpr (!kEstimated && NI == 0) {
            return StateSpace<0, NXL, NU, T>{.A = Matrix<0, 0, T>{}, .B = Matrix<0, NXL, T>{}, .C = Matrix<NU, 0, T>{}, .D = -K, .Ts = Ts};
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

    /** @brief The current-estimator law: the filter's estimate x̂(k−1|k−1) is the state ([x̂; d̂] with a disturbance model). */
    [[nodiscard]] constexpr StateSpace<NXL + ND, NY, NU, T> estimator_ss() const {
        constexpr size_t NF = NXL + ND;
        const auto&      sys = kalman.sys;
        const auto       Kf = filter_gain();
        const auto       ImMC = Matrix<NF, NF, T>::identity() - (kalman.M * sys.C);
        const auto       Ac = (ImMC * (sys.A - (sys.B * Kf))) + ((kalman.M * sys.D) * Kf);
        return {.A = Ac, .B = kalman.M, .C = -(Kf * Ac), .D = -(Kf * kalman.M), .Ts = Ts};
    }

    /** @brief The current-estimator law with integrators: the state is [x̂(k|k−1); ξ; u_prev], the inputs [r; y] (to_ss()). */
    [[nodiscard]] constexpr StateSpace<NXL + NI + NU, NI + NY, NU, T> integral_estimator_ss() const {
        constexpr size_t         NS = NXL + NI + NU;
        const auto&              sys = kalman.sys;
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
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, typename Kalman, size_t NZ, size_t NRU>
struct is_state_feedback_result<design::StateFeedbackResult<NX, NU, T, NI, NY, Kalman, NZ, NRU>> : damp::true_type {};

/**
 * @brief A StateFeedbackResult of any sizes. A runtime's constructor takes one and reads its kind
 *        from the result's own constants (kIntegral, kEstimated, ND), so a parameter added to the
 *        result does not reach the runtimes' constructors.
 */
template<typename R>
concept state_feedback_result = is_state_feedback_result<R>::value;

} // namespace detail

} // namespace damp
