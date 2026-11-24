// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
// IWYU pragma: private: a public header that exports it provides its names

/**
 * @file state_feedback_core.hpp
 * @brief The runtime the LQ laws share: detail::StateFeedbackCore, one tick with every feature
 *        (integral action, a filter, a steady-state target, delay states) optional
 *
 * The runtimes built on it are in lqr.hpp, lqi.hpp, lqg.hpp, lqgi.hpp and oflqg.hpp, the design that
 * feeds them is design::state_feedback (state_feedback.hpp), and the types it hands them are in
 * state_feedback_types.hpp. Law<kRes, T> is in lq_law.hpp, beside this header.
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/law_detail.hpp"
#include "damp/controllers/state_feedback_types.hpp"
#include "damp/design/status.hpp"
#include "damp/design/steady_state_target.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "damp/toolbox/bounds.hpp"

namespace damp {
namespace detail {

/**
 * @brief The integral action of a law: the integrators, the integrated outputs Ci (on the
 *        measurement) and their feedthrough Ci·D, the per-sample step, the error and the unclamped
 *        command held for the next predict, and with an estimator the measurement held from the
 *        last update (NYH = NY; 0: none). Their limits are the law's (StateFeedbackLimits::integral).
 */
template<size_t NI, size_t NC, size_t NU, size_t NYH, typename T, typename CiMatrix = Matrix<NI, NC, T>, typename CiuMatrix = Matrix<NI, NU, T>>
struct IntegralAction {
    using HeldY = maybe<(NYH > 0), ColVec<NYH, T>, OffSlot::held_y>;

    ColVec<NI, T>                   xi{};     ///< Integral of the tracking error
    [[no_unique_address]] CiMatrix  Ci{};     ///< Integrated outputs: e = r − Ci·(measurement)
    [[no_unique_address]] CiuMatrix Ci_u{};   ///< Ci·D: the integrated outputs' feedthrough of the input
    T                               step{1};  ///< Per-sample integrator gain: the design's Ts (xi += Ts·e ≈ ∫e dt)
    ColVec<NI, T>                   e{};      ///< Error held from control for the next predict, before the input's feedthrough
    [[no_unique_address]] HeldY     y{};      ///< Measurement from the last update (with an estimator)
    ColVec<NU, T>                   asked{};  ///< The law's last command before its limits: what the applied input is compared with
    Matrix<NI, NU, T>               unwind{}; ///< (∂u/∂ξ)⁺: what back-calculation moves ξ by per unit of a clamp's excess
};

/**
 * @brief The fade of a law designed with a Handover: the mismatch its command carries after a
 *        track, and the factor it decays by each tick, held (@p Held) or a constant of the type on
 *        a `constexpr` design's Law.
 */
template<size_t NU, typename T, bool Held>
struct HandoverFade {
    using Decay = maybe<Held, T, OffSlot::fade_decay>;

    ColVec<NU, T>               offset{};       ///< While pending, the input tracked; then e, the mismatch the command carries
    [[no_unique_address]] Decay decay{};        ///< exp(−Ts/fade_time): e's factor per tick
    bool                        pending{false}; ///< A track since the last command: the next takes e = u_tracked − u_law
};

/** @brief A gain design result other than StateFeedbackResult: K, success and Ts (PlaceResult, LQRResult). */
template<typename R>
concept gain_result = requires(const R& r) {
    r.success;
    r.Ts;
    r.K;
    requires MatrixLike<damp::remove_cvref_t<decltype(r.K)>>;
} && !is_state_feedback_result<R>::value;

/** @brief A gain design is discrete: its @c discrete flag (LQRResult), else a sample time Ts > 0 (PlaceResult, StateFeedbackResult). */
template<typename R>
[[nodiscard]] constexpr bool discrete_gain(const R& r) {
    if constexpr (requires { r.discrete; }) {
        return r.discrete;
    } else {
        return r.Ts != decltype(r.Ts){0};
    }
}

/** @brief Where a state-feedback law keeps its matrices: in the object, any value at run time. */
struct HeldGains {
    static constexpr bool constant = false; ///< The matrices are members

    template<size_t R, size_t C, typename T>
    using K = Matrix<R, C, T>; ///< The gain
    template<size_t R, size_t C, typename T>
    using Ci = Matrix<R, C, T>; ///< The integrated outputs
    template<size_t R, size_t C, typename T>
    using Ci_u = Matrix<R, C, T>; ///< Their feedthrough
    template<size_t R, size_t C, typename T>
    using TargetD = Matrix<R, C, T>; ///< The target per unit d̂
    template<size_t R, size_t C, typename T>
    using TargetR = Matrix<R, C, T>; ///< The target per unit r
    template<size_t N, typename T>
    using Demand = ::damp::detail::DemandTerms<N, T>; ///< The limit on what the outer terms ask of a state
};

/**
 * @brief The implementation the state-feedback runtimes share: u = u* − K·(x̂ − x*), every feature
 *        optional and free when off.
 *
 * LQR, LQI, LQG, LQGI and OFLQG wrap it, each exposing one update, one control and one step.
 *
 * | Feature | Parameter | Off |
 * | --- | --- | --- |
 * | integral action on r − Ci·y | NI, NY | 0 |
 * | estimator (FixedGainObserver with design::kalman's gain, disturbance model) | Est | NoEstimator: the state is measured |
 * | steady-state target for r, under the disturbance model if there is one (SteadyStateTarget) | NZ | 0 |
 * | input delay states: the inputs applied over the last NRU periods, held by the law (a measured state) | NRU | 0 |
 * | a measured disturbance d: the target under it, and the filter's known input after u | NMD | 0 |
 * | limits on the command and on the integrators | @ref limits and @ref set_limits: input, integral | unbounded |
 * | a limit on the command's move per period (InputRateLimits) | NDU, then @ref limits and @ref set_limits: input_rate | 0: none |
 * | back-calculation for the integrators against a clamp | @ref tracking_time | 0: conditional integration |
 * | a limit on the disturbance estimate d̂ (ND > 0) | @ref limits and @ref set_limits: disturbance | unbounded |
 * | a limit on what the outer terms ask of one state (DemandLimit) | @ref limits and @ref set_limits: demand | off |
 * | a fade of the mismatch after a @ref track (design::Handover) | HO | false: the first step is the law's own command |
 *
 * One tick is @ref step(): @ref predict (the time update with the input @ref back_calculate
 * recorded, else the last command), then @ref update with the measurement at sample k (the
 * state itself without an estimator), then @ref control, which returns u_k. Discrete-time
 * index: predict uses the input held over [k−1, k); update consumes m_k; control emits u_k.
 * The low-latency ISR is update(m_k) → control → write → back_calculate(u) → predict(): the
 * trailing predict is the next sample's time update. For a measured-state LQR, predict is
 * empty, so step(x) equals update then control.
 *
 * While another source drives the plant (an operator, another law), @ref track(u, …) takes step's
 * place, its arguments step's with the input applied first: the same predict and update, then,
 * instead of a command, the memory that makes control at this state, with these arguments, return
 * u. u is the applied input the next predict takes (as back_calculate(u) records it): the filter
 * predicts with it, and with the measured d, so d̂ keeps learning the load, the input's delay
 * states hold it, and a move limit (NDU > 0) is about it, so the first command after a track moves
 * at most input_rate from u. The integrators are preloaded through (∂u/∂ξ)⁺, the map
 * back-calculation moves them by, so that the command before the clamps is u: exactly with at
 * least as many integrators as inputs, the least-squares ξ (the command nearest u, Euclidean in
 * u's units) with fewer, then clamped into the integrators' limits, which leaves the command short
 * of u by what that clamp took. They compare the applied u with u itself at the next predict, so
 * they neither hold nor unwind on the gap between the law's command and u. A u outside the input
 * limits is one the law cannot emit: the preload is toward u clamped into them, the input the first
 * command hands over from, while the filter, the delay states and the move limit take the u the
 * plant received. There is no mode: the law tracks for as long as the caller calls track, and the
 * first step after it hands over. track adds no storage and no work to step. It runs the law once
 * (twice with integrators, the second at the preloaded ξ), so last_command() is the law's own
 * command at this state, u's mismatch with it the bump a handover without a preload would make. A
 * law ticked update → control → write → back_calculate → predict calls predict() once after the
 * last track, as step does at its start.
 *
 * A law with nothing to preload (LQR, LQG, an OFLQG away from a steady state) hands over with the
 * bump of its own command, within a move limit's box about u, unless its design took a
 * design::Handover (HO): then its first command after a track stores e = u_tracked − u_law,
 * u_tracked clamped into the input limits, and adds it, e decaying by exp(−Ts/fade_time) each tick;
 * the sum is clamped into the limits, and the integrators take e as a feedforward, not as
 * saturation (an LQI's, for what its preload leaves). Without HO there is no storage and no
 * per-tick work for it.
 *
 * With input delay states (NRU > 0, a design on a plant whose input arrives late) the law's model is
 * [x; v(k−1); …; v(k−NRU)]. On a measured state update takes the plant state x, and the law holds
 * the inputs it applied: @ref predict shifts the input back_calculate recorded (else the last
 * command) into the first delay state, so it runs every tick, as step does. A filter carries them as
 * states of its model, shifted by its own predict. The reference point of a delay state is u*: an
 * input held there is the one the law would send at the reference, so x_ref is the plant's.
 *
 * The law:
 * - about a reference point (x*, u*): u = u* − K·[x̂ − x*; xi]; x* = 0, u* = 0 is the regulator;
 * - with a target and no disturbance model (TrackedOutputs): the target's steady state for r,
 *   x_t = Nₓ·r and u_t = Nᵤ·r, adds to (x*, u*);
 * - with a disturbance model (ND > 0): the target is the steady state that holds the controlled
 *   outputs at r under d̂, taken about (x*, u*); the estimator runs in deviations from it
 *   (SteadyStateTarget::anchor), and u = u_a − K·(x̂ − x_a) at the deviation the command meets;
 * - with integral action: xi integrates r − Ci·y, xi += Ts·(r − Ci·y) per tick (xi ≈ ∫e dt),
 *   y = C·x + D·u the output over the period, u the input applied over it (the design's model):
 *   control holds the part it knows, r − Ci·C·x (from the measured state, or from the
 *   measurement less the feedthrough D·u of the input it was taken under), and @ref predict
 *   subtracts Ci·D·u for the input applied, then advances xi unless the applied input shows the
 *   actuator clamping against it (conditional integration), or with a @ref tracking_time advances
 *   it and pulls it toward the value that asks for the input applied (back-calculation), and
 *   clamps xi into the integrators' limits;
 * - with a measured disturbance (NMD > 0): the target is the steady state under it too, x_t and
 *   u_t moving with d (under [d̂; d] with a disturbance model, which then estimates only what d
 *   does not account for); a law whose r is not the target's (a regulator, an LQI) takes the
 *   steady state under d that holds every output, or the integrated outputs, at 0, its r still
 *   the integrators'. The control and step that take d record it, and a filter's next
 *   @ref predict holds it over the period with the input, as the plant does;
 * - within limits (design::StateFeedbackLimits): @ref control returns the command clamped into
 *   input, and that is the input the next @ref predict takes as applied unless @ref back_calculate
 *   recorded another. With a move limit (NDU > 0) the command is first clamped into
 *   u_prev + input_rate, u_prev the law's share of the input applied last tick (back_calculate's,
 *   else the last command; 0 after construction and @ref reset, as an MPC's held input), and where
 *   the two boxes do not meet the input box wins. The integrators compare the applied input with
 *   the command before either clamp, so they hold against the law's own limits as against a clamp
 *   downstream, and a filter and the delay states move on the clamped command. @ref limits
 *   reads the limits and @ref set_limits moves them at run time, for a rail that moves, and
 *   refuses one that is not a box. demand bounds what the outer terms ask of one state
 *   (DemandLimit), and the integrators hold against it the same way. With a disturbance model,
 *   @ref update clamps d̂ into disturbance, so a load past it (a wall the plant is held
 *   against) does not wind d̂, and with it the target, up the way an integrator would.
 *
 * A failed or continuous design leaves the law invalid at run time: valid() is false and control
 * returns 0, as every law does (when constant-evaluated it does not compile).
 *
 * Readings, references and reported inputs are the caller's to check: a NaN or ±Inf that reaches
 * the law stays in its memory (the integrators, a filter's estimate and its anchor) until that is
 * reset. @ref reset clears the integrators and the input memory; a filter restarts with
 * estimator.reset(x0, d0) (absolute, as state_estimate() reports it), or the whole law with
 * `law = decltype(law){}`. damp::isfinite tests on the
 * bits, so it holds under -ffast-math.
 *
 * @see "Optimal Control" (Anderson & Moore, 1990), §4.1 and §8
 * @see "Advanced PID Control" (Åström & Hägglund, 2006), §3.5: back-calculation anti-windup
 *
 * @tparam NX  Plant states
 * @tparam NU  Inputs
 * @tparam T   Scalar
 * @tparam NI  Integrators (0: none)
 * @tparam NY  Measured outputs (with no estimator, the outputs the integrators read through C)
 * @tparam Est Estimator: FixedGainObserver (the fixed Kalman gain of design::kalman), or NoEstimator for a measured state
 * @tparam NZ  Controlled outputs of the steady-state target (0: none)
 * @tparam Store Where the matrices live: HeldGains (members), or constants of the type (Law of a `constexpr` design)
 * @tparam NRU Input delay states the law holds, NRU·NU states after the plant's (0: none)
 * @tparam NMD  Measured-disturbance inputs (0: none): a filter takes them after u
 * @tparam NDU  Inputs with a move limit (an InputRateLimits: NU; 0: none)
 * @tparam HO   Whether the command fades the mismatch after a track (design::Handover; false: none)
 */
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, typename Est, size_t NZ, typename Store = HeldGains, size_t NRU = 0, size_t NMD = 0, size_t NDU = 0, bool HO = false>
struct StateFeedbackCore {
    using Traits = detail::estimator_traits<Est>;
    static constexpr bool   kEstimated = Traits::present;         ///< A filter estimates the state
    static constexpr size_t NXL = NX + (NRU * NU);                ///< The law's model: the plant, then its input delay states
    static constexpr size_t ND = Traits::ND;                      ///< Disturbance states of the estimator
    static constexpr size_t NXE = kEstimated ? Traits::NXE : NXL; ///< The filter's model before the disturbance (the law's without one)
    static constexpr size_t NIn = kEstimated ? NY : NX;           ///< What update takes: y, or the plant state
    static constexpr size_t NM = kEstimated ? NY : NXL;           ///< What the integrators read: y, or the law's state
    static constexpr size_t NR = (NI > 0) ? NI : NZ;              ///< The reference r: integrated or tracked outputs (0: none)
    static_assert(!kEstimated || Traits::NXE >= NXL, "The estimator is on the law's model, the delay states of late readings after it.");
    static_assert(NI == 0 || NZ == 0, "Integral action and a disturbance model are two offset-free choices: an IntegralAction or a disturbance spec, not both.");
    static_assert(ND == 0 || NZ > 0, "A disturbance model tracks through its steady-state target.");
    static_assert(!kEstimated || Traits::NUK == NU + NMD, "The filter takes the inputs, then the measured disturbance: design::state_feedback's with a MeasuredDisturbance.");
    static_assert(NDU == 0 || NDU == NU, "A move limit bounds every input: an InputRateLimits on the plant's NU.");
    static constexpr bool kTarget = (NZ > 0) || (NMD > 0);               ///< A steady-state target: for r, or under a measured d
    static constexpr bool kRecord = (NI > 0) || kEstimated || (NRU > 0); ///< The law keeps its command and the applied input
    using Record = damp::conditional_t<kRecord, detail::AppliedRecord<NU, T>, detail::NoAppliedRecord>;

    using Target = detail::maybe<
        kTarget,
        SteadyStateTarget<NXE, NU, ND + NMD, NZ, T, typename Store::template TargetD<NXE + NU, ND + NMD, T>, typename Store::template TargetR<NXE + NU, NZ, T>>,
        detail::OffSlot::target>;
    using Integral = detail::maybe<(NI > 0), detail::IntegralAction<NI, NM, NU, (kEstimated ? NY : 0), T, typename Store::template Ci<NI, NM, T>, typename Store::template Ci_u<NI, NU, T>>, detail::OffSlot::integral>;
    using HeldX = detail::maybe<!kEstimated, ColVec<NXL, T>, detail::OffSlot::held_x>;
    using DisturbanceLim = detail::maybe<(ND > 0), Bounds<ND, T>, detail::OffSlot::disturbance_lim>;
    using IntegralLim = detail::maybe<(NI > 0), Bounds<NI, T>, detail::OffSlot::integral_lim>;
    using TrackingTime = detail::maybe<(NI > 0), T, detail::OffSlot::tracking_time>;
    using HeldMeasured = detail::maybe<(kEstimated && NMD > 0), ColVec<NMD, T>, detail::OffSlot::held_measured>;
    using RateLim = detail::maybe<(NDU > 0), Bounds<NDU, T>, detail::OffSlot::input_rate>;
    using HeldCommand = detail::maybe<(NDU > 0 && !kRecord), ColVec<NU, T>, detail::OffSlot::held_command>; ///< The base of the move limit where no record keeps it
    using Demand = Store::template Demand<NX + NI, T>;                                                      ///< The demand limit as the store holds it: a member, or constants of the type
    using Fade = detail::maybe<HO, HandoverFade<NU, T, !Store::constant>, detail::OffSlot::handover>;
    using Limits = design::StateFeedbackLimits<NX, NU, T, NI, ND, NDU>; ///< What limits() returns and set_limits() takes

    [[no_unique_address]] Store::template K<NU, NXL + NI, T> K{}; ///< Gain on [x; xi]

    [[no_unique_address]] Est      estimator{}; ///< The filter (none for a measured state)
    [[no_unique_address]] Target   target{};    ///< The steady-state target (NZ > 0, or NMD > 0: under the measured d)
    [[no_unique_address]] Integral integral{};  ///< The integrators (NI > 0)

    [[no_unique_address]] TrackingTime tracking_time{}; ///< [s] The integrators' back-calculation tracking time (NI > 0; 0, the default: conditional integration; up to Ts: the whole excess each tick)

    constexpr StateFeedbackCore() = default;

    /** @brief From a gain on the law's state (no estimator): K = @p K_. */
    constexpr explicit StateFeedbackCore(const Matrix<NU, NXL + NI, T>& K_)
        requires(!kEstimated && !Store::constant)
        : K(K_) {}

    /**
     * @brief From a gain design (design::place's PlaceResult, LQRResult; any scalar): converts K.
     *
     * A failed or continuous (not discrete) result does not compile when constant-evaluated; at run
     * time either leaves the law invalid.
     */
    template<typename R>
        requires(detail::gain_result<R> && !kEstimated && NI == 0 && NZ == 0 && NRU == 0 && NMD == 0 && !Store::constant)
    constexpr StateFeedbackCore(const R& result) // NOLINT
        : valid_(usable(result, detail::discrete_gain(result))) {
        if (valid_) {
            K = result.K.template as<T>();
            if constexpr (requires { result.limits; }) {
                u_lim_ = result.limits.input.template as<T>();
                demand_lim_ = result.limits.demand.template as<T>();
            }
        }
    }

    /**
     * @brief From an LQ design (any scalar): the gain, and the filter, integrators and target the
     *        result carries.
     *
     * A failed or continuous result (Ts = 0) does not compile when constant-evaluated; at run
     * time either leaves the law invalid.
     */
    template<size_t NXR, size_t NUR, typename U, size_t NIR, size_t NYR, typename Kalman, size_t NZR, size_t NRR, size_t NMR, size_t NDR, bool HOR>
        requires(!Store::constant)
    constexpr StateFeedbackCore(const design::StateFeedbackResult<NXR, NUR, U, NIR, NYR, Kalman, NZR, NRR, NMR, NDR, HOR>& result) { // NOLINT
        constexpr bool same = NXR == NX && NUR == NU && NIR == NI && NYR == NY && NZR == NZ && NRR == NRU && NMR == NMD && NDR == NDU && HOR == HO
                           && damp::is_same_v<Est, detail::runtime_estimator_t<Kalman, T>>;
        static_assert(same, "This law's sizes are not the design's: let them deduce, LQG law{res}, or Law<res, T> for a member.");
        if constexpr (same) {
            valid_ = usable(result, damp::finite_positive(result.Ts));
            if (valid_) {
                K = result.K.template as<T>();
            }
            if constexpr (kEstimated) {
                estimator = Est(result.kalman);
            }
            if constexpr (NI > 0 && kEstimated) {
                integral.Ci = result.integral.Ci.template as<T>(); // on the measurement y
            } else if constexpr (NI > 0) {
                integral.Ci = result.integral.Ci_x.template as<T>(); // on the measured state: Ci·C
            }
            if constexpr (NI > 0) {
                integral.Ci_u = result.integral.Ci_u.template as<T>();
                integral.step = static_cast<T>(result.integral.step);
                integral.unwind = result.integral.unwind.template as<T>();
                integral_lim_ = result.limits.integral.template as<T>();
                tracking_time = static_cast<T>(result.tracking_time);
            }
            if constexpr (NDU > 0) {
                rate_lim_ = result.limits.input_rate.template as<T>();
            }
            u_lim_ = result.limits.input.template as<T>();
            demand_lim_ = result.limits.demand.template as<T>();
            if constexpr (kTarget) {
                target = Target(result.target);
            }
            if constexpr (HO) {
                if (valid_) {
                    handover_.decay = static_cast<T>(damp::exp(-result.Ts / result.fade_time)); // a design that succeeded: fade_time > 0
                }
            }
        }
    }

    // ---- the tick ------------------------------------------------------------------------

    /**
     * @brief Record the input the plant received, @p u_sat (after the actuator's clamp), and the
     *        feedforward @p u_ff added outside the law, for the next @ref predict.
     *
     * The clamp's excess is u_sat − (u_cmd + u_ff), u_cmd the law's last command before its own
     * limits, so a feedforward is not taken for saturation. The estimator moves on @p u_sat, and a
     * move limit's next box is about u_sat − u_ff, the law's share of it. A law that keeps no
     * memory (no integrator, estimator or move limit) ignores it. Without it the applied input is
     * the command control returned, clamped into its limits.
     */
    constexpr void back_calculate(const ColVec<NU, T>& u_sat, const ColVec<NU, T>& u_ff = ColVec<NU, T>{}) {
        if constexpr (kRecord) {
            applied_.record(u_sat, u_ff);
        } else if constexpr (NDU > 0) {
            last_ = ColVec<NU, T>(u_sat - u_ff); // the base of the next move: the law's share of what went out
        } else {
            static_cast<void>(u_sat);
            static_cast<void>(u_ff);
        }
    }

    /**
     * @brief Time update with the input the plant received over the last period (back_calculate's,
     *        else the last command): the integrators advance unless that would wind them up
     *        against a clamp, the estimator predicts.
     */
    constexpr void predict() {
        if constexpr (kRecord) {
            applied_.take();
        }
        if constexpr (NI > 0) {
            integrate();
        }
        if constexpr (kEstimated && NMD > 0) {
            estimator.predict(detail::stacked(applied_.applied, measured_)); // the measured d the last control took, held over the period
        } else if constexpr (kEstimated) {
            estimator.predict(applied_.applied);
        }
        if constexpr (NRU > 0 && !kEstimated) {
            for (size_t i = NXL; i-- > NX + NU;) { // oldest first: each delay state takes the one before it
                x_(i) = x_(i - NU);
            }
            for (size_t j = 0; j < NU; ++j) {
                x_(NX + j) = applied_.applied(j);
            }
        }
    }

    /**
     * @brief Measurement update: the filter's with the outputs @p m (the integrators read them at
     *        the next control), or, with no filter, the plant state @p m held for the next control
     *        (the input delay states after it are the law's own).
     */
    constexpr void update(const ColVec<NIn, T>& m) {
        if constexpr (kEstimated) {
            if constexpr (NI > 0) {
                integral.y = m;
            }
            estimator.update(m);
            if constexpr (ND > 0) {
                estimator.limit_disturbance(disturbance_lim_);
            }
        } else if constexpr (NRU == 0) {
            x_ = m;
        } else {
            for (size_t i = 0; i < NX; ++i) {
                x_(i) = m(i);
            }
        }
    }

    /**
     * @brief The law at the last update's measurement, about the reference point (@p x_ref, @p u_ref),
     *        clamped into its limits: integral action on @p r (NI > 0), or the steady-state target for
     *        @p r (NZ > 0, under d̂ when ND > 0).
     */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NR, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref)
        requires(NMD == 0)
    {
        if (!valid_) {
            return ColVec<NU, T>{}; // a refused design: no command
        }
        const ColVec<NXL, T> ref = detail::on_law_model<NRU>(x_ref, u_ref);
        if constexpr (ND > 0) {
            return anchored_law(r, ref, u_ref);
        } else if constexpr (NZ > 0) {
            const auto p = target.equilibrium(r); // on the law's model: each delay state at u_t
            return reference_law(law_state(), ColVec<NXL, T>(ref + detail::leading<NXL>(p.x)), ColVec<NU, T>(u_ref + p.u));
        } else if constexpr (NI > 0 && kEstimated) {
            // y was measured under the input the last predict applied: take its feedthrough out.
            return integral_law(r, law_state(), ref, u_ref, integral.y, ColVec<NI, T>(integral.Ci_u * applied_.applied));
        } else if constexpr (NI > 0) {
            return integral_law(r, x_, ref, u_ref, x_, ColVec<NI, T>{});
        } else {
            static_cast<void>(r);
            return reference_law(law_state(), ref, u_ref);
        }
    }

    /**
     * @brief control(@p r, @p x_ref, @p u_ref) under the measured disturbance @p d (NMD > 0): the target
     *        for @p r (none on a regulator; an LQI's integrators keep r) under d, and under d̂ too
     *        with a disturbance model, about (x*, u*). @p d is held for a filter's next predict.
     */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NR, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NMD, T>& d)
        requires(NMD > 0)
    {
        if constexpr (kEstimated) {
            measured_ = d; // the plant holds it over the period, as it holds the input
        }
        if (!valid_) {
            return ColVec<NU, T>{}; // a refused design: no command
        }
        const ColVec<NXL, T> ref = detail::on_law_model<NRU>(x_ref, u_ref);
        if constexpr (ND > 0) {
            target.anchor(estimator, r, detail::on_filter_model<NXE>(ref), u_ref, d); // under [d̂; d]
            const ColVec<NXL, T> z = detail::leading<NXL>(estimator.at_command());
            const ColVec<NU, T>  u = detail::leading<NU>(estimator.anchor_input()) - (K * z);
            return record(u, limit_demand(u, z));
        } else {
            const auto           p = target.equilibrium(target_reference(r), d); // on the law's model: each delay state at u_t
            const ColVec<NXL, T> x_a = ref + detail::leading<NXL>(p.x);
            const ColVec<NU, T>  u_a = u_ref + p.u;
            if constexpr (NI > 0 && kEstimated) {
                return integral_law(r, law_state(), x_a, u_a, integral.y, ColVec<NI, T>(integral.Ci_u * applied_.applied));
            } else if constexpr (NI > 0) {
                return integral_law(r, x_, x_a, u_a, x_, ColVec<NI, T>{});
            } else {
                return reference_law(law_state(), x_a, u_a);
            }
        }
    }

    /** @brief One tick: predict(), update(@p m), then control(@p r, @p x_ref, @p u_ref). */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NR, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NIn, T>& m)
        requires(NMD == 0)
    {
        predict();
        update(m);
        return control(r, x_ref, u_ref);
    }

    /** @brief One tick: predict(), update(@p m), then control(@p r, @p x_ref, @p u_ref, @p d): the reference and its feedforward, then what was read. */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NR, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NMD, T>& d, const ColVec<NIn, T>& m)
        requires(NMD > 0)
    {
        predict();
        update(m);
        return control(r, x_ref, u_ref, d);
    }

    /**
     * @brief In place of step(@p r, @p x_ref, @p u_ref, @p m) while another source drives the plant:
     *        predict(), update(@p m), then the memory that makes the next command start from @p u,
     *        the input applied now (the class comment).
     */
    constexpr void track(const ColVec<NU, T>& u, const ColVec<NR, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NIn, T>& m)
        requires(NMD == 0)
    {
        predict();
        update(m);
        follow(u, [&] { return control(r, x_ref, u_ref); });
    }

    /** @brief track(@p u, @p r, @p x_ref, @p u_ref, @p m) under the measured disturbance @p d, held for the filter's next predict. */
    constexpr void track(const ColVec<NU, T>& u, const ColVec<NR, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NMD, T>& d, const ColVec<NIn, T>& m)
        requires(NMD > 0)
    {
        predict();
        update(m);
        follow(u, [&] { return control(r, x_ref, u_ref, d); });
    }

    // ---- limits --------------------------------------------------------------------------

    /**
     * @brief Change the limits on u, on its move, on the integrators ξ and on the disturbance
     *        estimate d̂, and the limit on what the outer terms ask of a state
     *        (design::StateFeedbackLimits)
     *
     * The gain, the filter and the target keep their design; the next control clamps into the new
     * input, input_rate and demand, the next predict ξ into integral and the next update d̂ into
     * disturbance. Refused (false, limits unchanged) when an end is NaN or ±Inf, a lower end
     * is above its upper one, or the demand limit is not one this law can hold (a state that is not
     * the plant's or that is among its own outer terms, an index past the terms, one listed twice),
     * and on a `constexpr` design's Law, whose demand limit is a constant of the type, when it is
     * not the design's.
     */
    [[nodiscard]] constexpr bool set_limits(const Limits& lim) {
        if (!lim.input.valid() || !lim.input_rate.valid() || !lim.integral.valid() || !lim.disturbance.valid() || !lim.demand.template valid<NX>()) {
            return false;
        }
        if constexpr (Store::constant) {
            if (!lim.demand.same(demand_lim_.terms())) {
                return false; // a constant of the type
            }
        } else {
            demand_lim_ = Demand(lim.demand);
        }
        // Each limit is copied, then moved in: a copy assignment's self-test compares addresses,
        // which g++ cannot evaluate on a `constinit` Law under construction (its constructor sets
        // the design's limits here).
        u_lim_ = Bounds<NU, T>(lim.input);
        if constexpr (NDU > 0) {
            rate_lim_ = Bounds<NDU, T>(lim.input_rate);
        }
        if constexpr (NI > 0) {
            integral_lim_ = Bounds<NI, T>(lim.integral);
        }
        if constexpr (ND > 0) {
            disturbance_lim_ = Bounds<ND, T>(lim.disturbance);
        }
        return true;
    }

    /** @brief The limits in force (the design's, or the last set_limits). */
    [[nodiscard]] constexpr Limits limits() const {
        return Limits{u_lim_, held_bounds<NDU>(rate_lim_), held_bounds<NI>(integral_lim_), held_demand(), held_bounds<ND>(disturbance_lim_)};
    }

    // ---- state ---------------------------------------------------------------------------

    /** @brief The command the last control returned: the input a predict without back_calculate uses. */
    [[nodiscard]] constexpr const ColVec<NU, T>& last_command() const
        requires kRecord
    {
        return applied_.command;
    }

    /** @brief The input the last predict used: back_calculate's, else the command before it. */
    [[nodiscard]] constexpr const ColVec<NU, T>& applied() const
        requires kRecord
    {
        return applied_.applied;
    }

    /** @brief The estimated plant state x̂ (absolute). */
    [[nodiscard]] constexpr ColVec<NX, T> state_estimate() const
        requires kEstimated
    {
        return plant_estimate();
    }

    /** @brief The estimated disturbance d̂. */
    [[nodiscard]] constexpr ColVec<ND, T> disturbance_estimate() const
        requires(ND > 0)
    {
        return estimator.disturbance_estimate();
    }

    /** @brief The gain K on [x; xi]. */
    [[nodiscard]] constexpr const Matrix<NU, NXL + NI, T>& getK() const { return K; }

    /** @brief False when built from a design it refuses (failed, or continuous): control returns 0. */
    [[nodiscard]] constexpr bool valid() const { return valid_; }

    /** @brief Clear the integrators and the input memory (an estimate and its anchor stay). */
    constexpr void reset() {
        if constexpr (kRecord) {
            applied_ = Record{};
        }
        if constexpr (kEstimated) {
            estimator.reset();
        }
        if constexpr (NI > 0) {
            integral.xi = ColVec<NI, T>{};
            integral.e = ColVec<NI, T>{};
            integral.asked = ColVec<NU, T>{};
            if constexpr (kEstimated) {
                integral.y = ColVec<NY, T>{};
            }
        }
        if constexpr (NRU > 0 && !kEstimated) {
            x_ = ColVec<NXL, T>{};
        }
        if constexpr (kEstimated && NMD > 0) {
            measured_ = ColVec<NMD, T>{};
        }
        if constexpr (NDU > 0 && !kRecord) {
            last_ = ColVec<NU, T>{}; // the next move from 0, as an MPC's
        }
        if constexpr (HO) {
            handover_.offset = ColVec<NU, T>{}; // the factor is the design's
            handover_.pending = false;
        }
    }

private:
    Bounds<NU, T>                        u_lim_{};           ///< Limits on the command control returns (default unbounded)
    [[no_unique_address]] RateLim        rate_lim_{};        ///< Limits on the command's move per period (NDU > 0; default unbounded)
    [[no_unique_address]] HeldCommand    last_{};            ///< The base of the next move where no record keeps it: the last command, or back_calculate's share (NDU > 0)
    [[no_unique_address]] IntegralLim    integral_lim_{};    ///< Limits on the integrators ξ (NI > 0; default unbounded)
    [[no_unique_address]] DisturbanceLim disturbance_lim_{}; ///< Limits on the disturbance estimate d̂, in d's units (ND > 0; default unbounded)
    [[no_unique_address]] Demand         demand_lim_{};      ///< A limit on what the outer terms ask of a state (default off; a constant of a `constexpr` design's Law)
    [[no_unique_address]] Record         applied_{};         ///< The last command and the input back_calculate recorded
    [[no_unique_address]] HeldX          x_{};               ///< The measured state from the last update (no estimator)
    [[no_unique_address]] HeldMeasured   measured_{};        ///< The measured disturbance the last control took, for the filter's next predict (NMD > 0)
    [[no_unique_address]] Fade           handover_{};        ///< The mismatch the command fades after a track (HO)
    bool                                 valid_{true};       ///< False when built from a design it refuses

    /** @brief A limit the law holds (N > 0), or the empty box of one it does not have: what limits() reports. */
    template<size_t N, typename Held>
    [[nodiscard]] static constexpr Bounds<N, T> held_bounds(const Held& held) {
        if constexpr (N > 0) {
            return held;
        } else {
            static_cast<void>(held);
            return Bounds<0, T>{};
        }
    }

    /** @brief The demand limit as limits() reports it: the member, or a `constexpr` design's constant. */
    [[nodiscard]] constexpr ::damp::detail::DemandTerms<NX + NI, T> held_demand() const {
        if constexpr (Store::constant) {
            return demand_lim_.terms();
        } else {
            return demand_lim_;
        }
    }

    /** @brief The target's reference: @p r when the target holds it, none when the law's r is the integrators' (NI > 0) or it has none. */
    [[nodiscard]] static constexpr ColVec<NZ, T> target_reference(const ColVec<NR, T>& r) {
        if constexpr (NI > 0) {
            static_cast<void>(r);
            return ColVec<NZ, T>{};
        } else {
            return r;
        }
    }

    /** @brief Whether a design's gain runs on this tick: it succeeded and is discrete. Otherwise it does not compile when constant-evaluated (the error names its status) and the law is invalid. */
    template<typename R>
    [[nodiscard]] static constexpr bool usable(const R& result, bool discrete) {
        const bool ok = result.success && discrete;
        if (!ok && damp::is_constant_evaluated()) {
            if constexpr (requires { result.status; }) {
                detail::refuse_failed_design(result.success ? design::DesignStatus::invalid_sample_time : result.status);
            } else {
                detail::refuse_failed_design(); // a gain result with no status
            }
        }
        return ok;
    }

    [[nodiscard]] constexpr ColVec<NX, T> plant_estimate() const
        requires kEstimated
    {
        return detail::leading<NX>(estimator.state_estimate());
    }

    [[nodiscard]] constexpr ColVec<NXL, T> law_state() const {
        if constexpr (kEstimated) {
            return detail::leading<NXL>(estimator.state_estimate());
        } else {
            return x_;
        }
    }

    [[nodiscard]] constexpr auto integrators() const {
        if constexpr (NI > 0) {
            return integral.xi;
        } else {
            return ColVec<0, T>{};
        }
    }

    /**
     * @brief The law u = u* − K·z with the outer terms' demand on the demand limit's state clamped into
     *        its box, per input that acts on it; u itself where the demand is inside.
     */
    [[nodiscard]] constexpr ColVec<NU, T> limit_demand(const ColVec<NU, T>& u, const ColVec<NXL + NI, T>& z) const {
        const size_t s = demand_lim_.state;
        if (s >= NX) {
            return u; // off (none), or not a plant state
        }
        ColVec<NU, T> out = u;
        for (size_t i = 0; i < NU; ++i) {
            const T k_s = K(i, s);
            T       outer{0}; // −Σ K_j·z_j over the outer terms: the command they make
            T       row{0};   // the row's largest gain, for the test on k_s
            for (size_t j = 0; j < NX + NI; ++j) {
                const size_t c = (j < NX) ? j : (NXL + (j - NX)); // the integrators sit after the delay states
                row = damp::max(row, damp::abs(K(i, c)));
                if (demand_lim_.outer(j)) {
                    outer -= K(i, c) * z(c);
                }
            }
            if (damp::abs(k_s) <= T{4} * eps<T>() * row) {
                continue; // this input does not act on s: its gain there is rounding
            }
            const T demand = outer / k_s;
            const T allowed = demand_lim_.lim.saturate(demand);
            if (allowed != demand) {
                out(i) += (k_s * allowed) - outer;
            }
        }
        return out;
    }

    /**
     * @brief The command: @p u clamped into its limits, with HO the handover's mismatch added first. The
     *        law keeps it, and @p asked (u before any clamp, the mismatch with it) for the integrators.
     */
    [[nodiscard]] constexpr ColVec<NU, T> record(const ColVec<NU, T>& asked, const ColVec<NU, T>& u) {
        if constexpr (HO) {
            const ColVec<NU, T> e = fade(u);
            return emit(ColVec<NU, T>(asked + e), ColVec<NU, T>(u + e)); // a feedforward to the integrators, not saturation
        } else {
            return emit(asked, u);
        }
    }

    /** @brief The handover's mismatch this command carries: e = u_tracked − @p u on the first after a track, then e decayed once per command. */
    constexpr ColVec<NU, T> fade(const ColVec<NU, T>& u)
        requires HO
    {
        if (handover_.pending) {
            handover_.offset = ColVec<NU, T>(handover_.offset - u);
            handover_.pending = false;
        }
        const ColVec<NU, T> e = handover_.offset;
        if constexpr (Store::constant) {
            handover_.offset = ColVec<NU, T>(e * Store::kFadeDecay);
        } else {
            handover_.offset = ColVec<NU, T>(e * handover_.decay);
        }
        return e;
    }

    /**
     * @brief The memory a track leaves: the law's command at this state (@p law runs control), the
     *        integrators preloaded toward @p u clamped into the input limits and the command again
     *        at them, and @p u recorded as the input applied: the base of the next move, and the
     *        command the integrators compare it with. With HO the next command takes the mismatch
     *        with that clamped u.
     */
    template<typename Control>
    constexpr void follow(const ColVec<NU, T>& u, const Control& law) {
        const ColVec<NU, T> from = u_lim_.saturate(u); // the law cannot emit more: what it hands over from
        if constexpr (HO) {
            handover_.offset = ColVec<NU, T>{}; // the law's own command, with no fade
            handover_.pending = false;
        }
        if constexpr (kRecord) {
            static_cast<void>(law());
            if constexpr (NI > 0) {
                if (valid_) {
                    const ColVec<NI, T> xi = integral.xi + (integral.unwind * ColVec<NU, T>(from - integral.asked));
                    integral.xi = integral_lim_.saturate(xi); // least squares with NI < NU, short of u where the limit binds
                    static_cast<void>(law());                 // the error, the command before the clamps and the command at that ξ
                }
                integral.asked = u; // what the applied u is compared with: no excess to hold or unwind
            }
            applied_.record(u, ColVec<NU, T>{}); // all of it the law's from here: the next move is from u
        } else if constexpr (NDU > 0) {
            static_cast<void>(law); // no memory but the move's base
            last_ = u;
        } else {
            static_cast<void>(law); // no memory: the state is all update holds
        }
        if constexpr (HO) {
            handover_.offset = from;
            handover_.pending = true;
        }
    }

    /**
     * @brief The command: @p u clamped into its move's box about u_prev (NDU > 0), then into the
     *        input box, which wins where the two do not meet. The law keeps it, and @p asked (u
     *        before any clamp) for the integrators.
     */
    [[nodiscard]] constexpr ColVec<NU, T> emit(const ColVec<NU, T>& asked, const ColVec<NU, T>& u) {
        const ColVec<NU, T> out = u_lim_.saturate(moved(u));
        if constexpr (NI > 0) {
            integral.asked = asked;
        } else {
            static_cast<void>(asked);
        }
        if constexpr (kRecord) {
            applied_.command = out;
        } else if constexpr (NDU > 0) {
            last_ = out;
        }
        return out;
    }

    /**
     * @brief @p u with each move from u_prev clamped into the move limit (NDU > 0), @p u itself
     *        without one. u_prev is the law's share of the input applied last tick: the record's
     *        (back_calculate's, else the last command, as the last predict took it), or on a law
     *        that keeps no record the last command, or what back_calculate reported since.
     */
    [[nodiscard]] constexpr decltype(auto) moved(const ColVec<NU, T>& u) const {
        if constexpr (NDU == 0) {
            return (u);
        } else {
            ColVec<NU, T> out = u;
            for (size_t j = 0; j < NU; ++j) {
                T base{};
                if constexpr (kRecord) {
                    base = applied_.applied(j) - applied_.feedforward(j);
                } else {
                    base = last_(j);
                }
                const T move = u(j) - base; // an unbounded end is never passed, so no sum with ±max() is formed
                if (move > rate_lim_.upper(j)) {
                    out(j) = base + rate_lim_.upper(j);
                } else if (move < rate_lim_.lower(j)) {
                    out(j) = base + rate_lim_.lower(j);
                }
            }
            return out;
        }
    }

    /** @brief u = u* − K·z, z = [x − x*; xi], x and x* on the law's model. */
    [[nodiscard]] constexpr ColVec<NU, T> reference_law(const ColVec<NXL, T>& x, const ColVec<NXL, T>& x_ref, const ColVec<NU, T>& u_ref) {
        const ColVec<NXL + NI, T> z = detail::stacked(ColVec<NXL, T>(x - x_ref), integrators());
        const ColVec<NU, T>       u = u_ref - (K * z);
        return record(u, limit_demand(u, z));
    }

    /** @brief The law, then the integrators' error from @p r and the measurement, with the feedthrough @p fed (Ci·D·u) it carries added back. */
    template<size_t NMeas>
    [[nodiscard]] constexpr ColVec<NU, T> integral_law(
        const ColVec<NI, T>&    r,
        const ColVec<NXL, T>&   x,
        const ColVec<NXL, T>&   x_ref,
        const ColVec<NU, T>&    u_ref,
        const ColVec<NMeas, T>& meas,
        const ColVec<NI, T>&    fed
    ) {
        const ColVec<NU, T> u = reference_law(x, x_ref, u_ref);
        integral.e = ColVec<NI, T>((r - (integral.Ci * meas)) + fed);
        return u;
    }

    /** @brief The estimator anchored on the target for @p r about (x*, u*), then u = u_a − K·(the deviation the command meets). */
    [[nodiscard]] constexpr ColVec<NU, T> anchored_law(const ColVec<NZ, T>& r, const ColVec<NXL, T>& x_ref, const ColVec<NU, T>& u_ref)
        requires(ND > 0)
    {
        target.anchor(estimator, r, detail::on_filter_model<NXE>(x_ref), u_ref);
        const ColVec<NXL, T> z = detail::leading<NXL>(estimator.at_command());
        const ColVec<NU, T>  u = estimator.anchor_input() - (K * z);
        return record(u, limit_demand(u, z));
    }

    /**
     * @brief Integrate the held error less the applied input's feedthrough Ci·D·u, into the integrators' limits,
     *        against a clamp (the input and move limits, the demand limit, or one back_calculate reported): integrator i holds
     *        when −K(j, i)·eᵢ points into a clamped input j (conditional integration), or with a
     *        tracking time T_t > 0 integrates and takes back min(Ts/T_t, 1) of the excess through
     *        (∂u/∂ξ)⁺ (back-calculation, IntegralAction::tracking_time).
     */
    constexpr void integrate()
        requires(NI > 0)
    {
        const ColVec<NI, T> e = integral.e - (integral.Ci_u * applied_.applied);
        ColVec<NU, T>       excess{};
        for (size_t j = 0; j < NU; ++j) {
            const T u_req = integral.asked(j) + applied_.feedforward(j);
            const T diff = applied_.applied(j) - u_req;
            const T tol = T{4} * eps<T>() * (damp::abs(applied_.applied(j)) + damp::abs(u_req));
            excess(j) = (damp::abs(diff) > tol) ? diff : T{0};
        }
        if (tracking_time > T{0}) {
            const T share = damp::min(integral.step / tracking_time, T{1}); // the whole excess at most: past it, an overshoot
            integral.xi = ColVec<NI, T>((integral.xi + (integral.step * e)) + (share * (integral.unwind * excess)));
        } else {
            for (size_t i = 0; i < NI; ++i) {
                bool winds_up = false;
                for (size_t j = 0; j < NU; ++j) {
                    winds_up = winds_up || (K(j, NXL + i) * e(i) * excess(j) > T{0});
                }
                if (!winds_up) {
                    integral.xi(i) += integral.step * e(i);
                }
            }
        }
        integral.xi = integral_lim_.saturate(integral.xi); // held or not: a limit narrowed at run time pulls xi in
        integral.e = ColVec<NI, T>{};
    }
};

} // namespace detail
} // namespace damp
