// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file state_feedback_runtime.hpp
 * @brief StateFeedback: the runtime of every LQ design, and the Law that names it
 *
 * design::state_feedback(plant, [Ts,] weights, opts...) (state_feedback.hpp) → StateFeedbackResult
 * → StateFeedback, whichever options the design took. LQR, LQI, LQG, LQGI and OFLQG (lqr.hpp,
 * lqi.hpp, lqg.hpp, lqgi.hpp, oflqg.hpp) name what the options produce: each is an alias of
 * StateFeedback that takes only its own kind of design.
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/lq_law.hpp" // IWYU pragma: export
#include "damp/controllers/detail/state_feedback_core.hpp"
#include "damp/controllers/law.hpp"                  // IWYU pragma: export
#include "damp/controllers/state_feedback_types.hpp" // IWYU pragma: export
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/core.hpp"

namespace damp {

/**
 * @ingroup discrete_controllers
 * @brief The runtime of every LQ design, u = u* − K·(x̂ − x*) on the measured state or a filter's estimate: LQR, LQI, LQG, LQGI and OFLQG name its kinds
 *
 * Runs design::state_feedback, whichever options it took, with LQRWeights or Poles, and a gain K of
 * its own: `StateFeedback law{res.as<float>()}` deduces the law from the result, and
 * `Law<kRes, float>` names it for a struct member. The options pick the kind, and each kind's name
 * is an alias of this type that takes only that kind's design, so `LQG law{res}` refuses an LQI's:
 *
 * | The design's options | Kind | update(m) takes | r |
 * | --- | --- | --- | --- |
 * | none | LQR | the plant state x | none: the regulator |
 * | TrackedOutputs | LQR | x | the tracked outputs, held through the steady-state target |
 * | IntegralAction | LQI | x | the integrated outputs |
 * | KalmanNoise or ObserverPoles | LQG | the reading y | none, or the tracked outputs with TrackedOutputs |
 * | a filter and an IntegralAction | LQGI | y | the integrated outputs |
 * | a filter and a disturbance spec | OFLQG | y | the controlled outputs, held through the target under d̂ |
 *
 * A MeasuredDisturbance, an InputDelay, an OutputDelay (with a filter), InputLimits,
 * InputRateLimits, a DemandLimit or a Handover adds to any kind and leaves it that kind. What the
 * law has follows from its design:
 *
 * | Member | What it is | Present |
 * | --- | --- | --- |
 * | update(m) | the measurement update: hold the state x, or the filter's with the reading y | always |
 * | control(), step(m) | the regulator at the last update; one tick | no r and no measured d |
 * | control(x_ref, u_ref), step(x_ref, u_ref, m) | the same about a planned trajectory (x*, u*) | no r and no measured d |
 * | control(r), step(r, m) | the law holding r; one tick | r: integrators (NI > 0) or a target (NZ > 0) |
 * | control(r, x_ref, u_ref), step(r, x_ref, u_ref, m) | the same about (x*, u*), which adds to the target's | r |
 * | control([r,] [x_ref, u_ref,] d), step([r,] [x_ref, u_ref,] d, m) | each form under the measured disturbance d, after the reference and its feedforward | a MeasuredDisturbance (NM > 0), in place of the forms without d |
 * | track(u, …) | in place of each step while another source drives the plant: step's arguments with the input applied first | beside each step |
 * | predict(), back_calculate(u_sat, u_ff) | the time update on the input applied; the input a clamp after the law let through | always |
 * | limits(), set_limits(lim) | the limits in force (design::StateFeedbackLimits); change them at run time | always |
 * | reset(), valid(), getK() | clear the integrators and the input memory; whether the design ran; K | always |
 * | K | the gain on [x; ξ] | always |
 * | last_command(), applied() | the last command; the input the last predict used | the law keeps them: integrators, a filter or delay states |
 * | integral, tracking_time | the integrators; their back-calculation tracking time | NI > 0: LQI, LQGI |
 * | estimator, state_estimate() | the FixedGainObserver with the design's fixed gain; x̂ | a filter: LQG, LQGI, OFLQG |
 * | target | the steady-state target | NZ > 0 (TrackedOutputs, OFLQG) or NM > 0 |
 * | disturbance_estimate() | d̂ | a disturbance model: OFLQG |
 *
 * A data member the design does not have is an empty type that takes no storage (`estimator` a
 * NoEstimator, the others a detail::Off), and a function it does not have does not compile.
 *
 * The tick is step, predict then update then control, or the low-latency order update(m),
 * control, write, back_calculate, predict(): the trailing predict is the next sample's time update.
 * On a measured state without delay states predict is empty, so step(x) equals update then
 * control. The regulator is x* = 0 and u* = 0.
 *
 * The command is clamped into the design's InputLimits, and its move into an InputRateLimits
 * (unbounded by default, and for a gain K of its own); a DemandLimit bounds what the outer terms
 * ask of a state; the integrators are clamped into IntegralAction::limits, and d̂ into
 * StateFeedbackLimits::disturbance (the runtime's alone). limits() reads them and set_limits()
 * moves them at run time, as on an MPC. predict() holds the integrators against those clamps
 * (conditional integration, or back-calculation with a tracking_time) and moves the filter and the
 * delay states on the clamped command, so a clamp is not taken for a disturbance;
 * back_calculate(u_sat) before it reports a stage that clamps after the law. A failed or
 * continuous design does not compile when constant-evaluated; at run time valid() is false and
 * control returns 0.
 *
 * With integral action, xi integrates r − Ci·y: on the measured state Ci·C·x, through a filter the
 * last update's reading; the law is u = u* − K·[x̂ − x*; xi]. With a disturbance model the filter is
 * on [x; d] and holds deviations from the target, which holds the controlled outputs at r under d̂,
 * taken about (x*, u*). With a measured disturbance the target is the steady state under the d
 * the tick takes too, and a filter takes d as a known input after u.
 *
 * A design with an InputDelay (NRU > 0) runs on [x; v(k−1); …; v(k−NRU)]. On the measured state
 * update takes the plant state, and the law holds the delay states: predict shifts the input
 * back_calculate recorded (else the last command) into them, so it runs every tick. A filter carries
 * them as states of its model, shifted by its own predict. x_ref is the plant's; the delay states'
 * reference is u_ref.
 *
 * While another source drives the plant (an operator, another law), track(u, …) takes the place of
 * each step(…): its predict and update, then the memory that starts the next command from u. An
 * LQI's or LQGI's integrators are preloaded, a filter predicts with u, and a law with nothing to
 * preload hands over with its own command unless its design took a design::Handover, which fades
 * the mismatch (detail::StateFeedbackCore says how).
 *
 * @code
 * static constinit StateFeedback law{kRes.as<float>()}; // LQI law{…} is the same type, and refuses another kind's design
 * law.update(x);
 * u = law.control(r);
 * @endcode
 *
 * Deducing through one of the five names (`LQG law{res}`) is C++20's deduction for alias
 * templates, which GCC has from 10 and Clang from 19.
 *
 * @note Compare with MATLAB®'s lqr / lqrd, lqi, lqg and lqgtrack.
 *
 * @tparam NX    Plant states
 * @tparam NU    Inputs
 * @tparam T     Scalar
 * @tparam NI    Integrators (0: none)
 * @tparam NY    Measured outputs (NX for a law without integrators or a filter)
 * @tparam Est   The filter: a FixedGainObserver, or NoEstimator for a measured state
 * @tparam NZ    Controlled outputs of the steady-state target: the size of r without integrators (0: none)
 * @tparam Store Where the matrices live: detail::HeldGains (members), or constants of the type (the Law of a `constexpr` design)
 * @tparam NRU   Input delay states per input of the law's model, NRU·NU after the plant's (0: none)
 * @tparam NM    Measured-disturbance inputs (0: none)
 * @tparam NDU   Inputs with a move limit (an InputRateLimits: NU; 0: none)
 * @tparam HO    Whether the law fades the mismatch after a track (a Handover; false: none)
 */
template<size_t NX, size_t NU, typename T = float, size_t NI = 0, size_t NY = NX, typename Est = NoEstimator, size_t NZ = 0, typename Store = detail::HeldGains, size_t NRU = 0, size_t NM = 0, size_t NDU = 0, bool HO = false>
struct StateFeedback : private detail::StateFeedbackCore<NX, NU, T, NI, NY, Est, NZ, Store, NRU, NM, NDU, HO> {
    using Core = detail::StateFeedbackCore<NX, NU, T, NI, NY, Est, NZ, Store, NRU, NM, NDU, HO>; ///< The shared implementation

    constexpr StateFeedback() = default;

    /** @brief From a gain K on the measured state. */
    constexpr explicit StateFeedback(const Matrix<NU, NX, T>& K_)
        requires(!Core::kEstimated && NI == 0 && NZ == 0 && NRU == 0 && NM == 0 && !HO && !Store::constant)
        : Core(K_) {}

    /** @brief From design::state_feedback, any options and any scalar: the gain, and the filter, integrators, target and limits the result carries. */
    template<detail::state_feedback_result R>
        requires(!Store::constant)
    constexpr StateFeedback(const R& result) : Core(result) {} // NOLINT

    using Core::estimator;
    using Core::integral;
    using Core::K;
    using Core::target;
    using Core::tracking_time;

    using Core::update;

    // ---- no r, no measured d: the regulator ------------------------------------------------

    /** @brief u = −K·x̂ at the last update: the regulator. */
    [[nodiscard]] constexpr ColVec<NU, T> control()
        requires(Core::NR == 0 && NM == 0)
    {
        return Core::control(ColVec<0, T>{}, ColVec<NX, T>{}, ColVec<NU, T>{});
    }

    /** @brief u = u* − K·(x̂ − x*) at the last update: feedforward from a planned trajectory (x*, u*). */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref)
        requires(Core::NR == 0 && NM == 0)
    {
        return Core::control(ColVec<0, T>{}, x_ref, u_ref);
    }

    /** @brief One tick: predict(), update(@p m), then control(): the regulator. */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<Core::NIn, T>& m)
        requires(Core::NR == 0 && NM == 0)
    {
        return Core::step(ColVec<0, T>{}, ColVec<NX, T>{}, ColVec<NU, T>{}, m);
    }

    /** @brief One tick: predict(), update(@p m), then control(@p x_ref, @p u_ref), the feedforward form. */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<Core::NIn, T>& m)
        requires(Core::NR == 0 && NM == 0)
    {
        return Core::step(ColVec<0, T>{}, x_ref, u_ref, m);
    }

    /** @brief In place of step(m) while another source drives the plant: its predict() and update(m), then the memory that starts the next command from @p u, the input applied. */
    constexpr void track(const ColVec<NU, T>& u, const ColVec<Core::NIn, T>& m)
        requires(Core::NR == 0 && NM == 0)
    {
        Core::track(u, ColVec<0, T>{}, ColVec<NX, T>{}, ColVec<NU, T>{}, m);
    }

    /** @brief In place of step(x_ref, u_ref, m) while another source drives the plant: its predict() and update(m), then the memory that starts the next command from @p u, the input applied. */
    constexpr void track(const ColVec<NU, T>& u, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<Core::NIn, T>& m)
        requires(Core::NR == 0 && NM == 0)
    {
        Core::track(u, ColVec<0, T>{}, x_ref, u_ref, m);
    }

    // ---- r, no measured d ------------------------------------------------------------------

    /** @brief The law holding @p r at the last update: on the integrated outputs, or the target's for the tracked or controlled outputs (under d̂ with a disturbance model). */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<Core::NR, T>& r)
        requires(Core::NR > 0 && NM == 0)
    {
        return Core::control(r, ColVec<NX, T>{}, ColVec<NU, T>{});
    }

    /** @brief control(@p r) about the planned trajectory (x*, u*), which adds to the target's. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<Core::NR, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref)
        requires(Core::NR > 0 && NM == 0)
    {
        return Core::control(r, x_ref, u_ref);
    }

    /** @brief One tick: predict(), update(@p m), then control(@p r). The basic call. */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<Core::NR, T>& r, const ColVec<Core::NIn, T>& m)
        requires(Core::NR > 0 && NM == 0)
    {
        return Core::step(r, ColVec<NX, T>{}, ColVec<NU, T>{}, m);
    }

    /** @brief One tick: predict(), update(@p m), then control(@p r, @p x_ref, @p u_ref), the feedforward form. */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<Core::NR, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<Core::NIn, T>& m)
        requires(Core::NR > 0 && NM == 0)
    {
        return Core::step(r, x_ref, u_ref, m);
    }

    /** @brief In place of step(r, m) while another source drives the plant: its predict() and update(m), then the memory that starts the next command from @p u, the input applied. */
    constexpr void track(const ColVec<NU, T>& u, const ColVec<Core::NR, T>& r, const ColVec<Core::NIn, T>& m)
        requires(Core::NR > 0 && NM == 0)
    {
        Core::track(u, r, ColVec<NX, T>{}, ColVec<NU, T>{}, m);
    }

    /** @brief In place of step(r, x_ref, u_ref, m) while another source drives the plant: its predict() and update(m), then the memory that starts the next command from @p u, the input applied. */
    constexpr void track(const ColVec<NU, T>& u, const ColVec<Core::NR, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<Core::NIn, T>& m)
        requires(Core::NR > 0 && NM == 0)
    {
        Core::track(u, r, x_ref, u_ref, m);
    }

    // ---- no r, a measured d: the regulator under it ------------------------------------------

    /** @brief u = u_t(d) − K·(x̂ − x_t(d)): the regulator under the measured disturbance @p d, every output held at 0. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NM, T>& d)
        requires(Core::NR == 0 && NM > 0)
    {
        return Core::control(ColVec<0, T>{}, ColVec<NX, T>{}, ColVec<NU, T>{}, d);
    }

    /** @brief control(@p d) about the planned trajectory (x*, u*). */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NM, T>& d)
        requires(Core::NR == 0 && NM > 0)
    {
        return Core::control(ColVec<0, T>{}, x_ref, u_ref, d);
    }

    /** @brief One tick: predict(), update(@p m), then control(@p d). */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NM, T>& d, const ColVec<Core::NIn, T>& m)
        requires(Core::NR == 0 && NM > 0)
    {
        return Core::step(ColVec<0, T>{}, ColVec<NX, T>{}, ColVec<NU, T>{}, d, m);
    }

    /** @brief One tick: predict(), update(@p m), then control(@p x_ref, @p u_ref, @p d). */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NM, T>& d, const ColVec<Core::NIn, T>& m)
        requires(Core::NR == 0 && NM > 0)
    {
        return Core::step(ColVec<0, T>{}, x_ref, u_ref, d, m);
    }

    /** @brief In place of step(d, m) while another source drives the plant: its predict() and update(m), then the memory that starts the next command from @p u, the input applied. */
    constexpr void track(const ColVec<NU, T>& u, const ColVec<NM, T>& d, const ColVec<Core::NIn, T>& m)
        requires(Core::NR == 0 && NM > 0)
    {
        Core::track(u, ColVec<0, T>{}, ColVec<NX, T>{}, ColVec<NU, T>{}, d, m);
    }

    /** @brief In place of step(x_ref, u_ref, d, m) while another source drives the plant: its predict() and update(m), then the memory that starts the next command from @p u, the input applied. */
    constexpr void track(const ColVec<NU, T>& u, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NM, T>& d, const ColVec<Core::NIn, T>& m)
        requires(Core::NR == 0 && NM > 0)
    {
        Core::track(u, ColVec<0, T>{}, x_ref, u_ref, d, m);
    }

    // ---- r and a measured d ------------------------------------------------------------------

    /** @brief The law holding @p r under the measured disturbance @p d: the target under d, after the reference. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<Core::NR, T>& r, const ColVec<NM, T>& d)
        requires(Core::NR > 0 && NM > 0)
    {
        return Core::control(r, ColVec<NX, T>{}, ColVec<NU, T>{}, d);
    }

    /** @brief control(@p r, @p d) about the planned trajectory (x*, u*): the reference and its feedforward, then the measured disturbance. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<Core::NR, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NM, T>& d)
        requires(Core::NR > 0 && NM > 0)
    {
        return Core::control(r, x_ref, u_ref, d);
    }

    /** @brief One tick: predict(), update(@p m), then control(@p r, @p d). */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<Core::NR, T>& r, const ColVec<NM, T>& d, const ColVec<Core::NIn, T>& m)
        requires(Core::NR > 0 && NM > 0)
    {
        return Core::step(r, ColVec<NX, T>{}, ColVec<NU, T>{}, d, m);
    }

    /** @brief One tick: predict(), update(@p m), then control(@p r, @p x_ref, @p u_ref, @p d). */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<Core::NR, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NM, T>& d, const ColVec<Core::NIn, T>& m)
        requires(Core::NR > 0 && NM > 0)
    {
        return Core::step(r, x_ref, u_ref, d, m);
    }

    /** @brief In place of step(r, d, m) while another source drives the plant: its predict() and update(m), then the memory that starts the next command from @p u, the input applied. */
    constexpr void track(const ColVec<NU, T>& u, const ColVec<Core::NR, T>& r, const ColVec<NM, T>& d, const ColVec<Core::NIn, T>& m)
        requires(Core::NR > 0 && NM > 0)
    {
        Core::track(u, r, ColVec<NX, T>{}, ColVec<NU, T>{}, d, m);
    }

    /** @brief In place of step(r, x_ref, u_ref, d, m) while another source drives the plant: its predict() and update(m), then the memory that starts the next command from @p u, the input applied. */
    constexpr void track(const ColVec<NU, T>& u, const ColVec<Core::NR, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NM, T>& d, const ColVec<Core::NIn, T>& m)
        requires(Core::NR > 0 && NM > 0)
    {
        Core::track(u, r, x_ref, u_ref, d, m);
    }

    using Core::applied;
    using Core::back_calculate;
    using Core::disturbance_estimate;
    using Core::getK;
    using Core::last_command;
    using Core::limits;
    using Core::predict;
    using Core::reset;
    using Core::set_limits;
    using Core::state_estimate;
    using Core::valid;
};

/** @brief StateFeedback law{res}: the sizes, scalar, filter and target of design::state_feedback's result. */
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, typename Kalman, size_t NZ, size_t NRU, size_t NM, size_t NDU, bool HO>
StateFeedback(const design::StateFeedbackResult<NX, NU, T, NI, NY, Kalman, NZ, NRU, NM, NDU, HO>&)
    -> StateFeedback<NX, NU, T, NI, NY, detail::runtime_estimator_t<Kalman, T>, NZ, detail::HeldGains, NRU, NM, NDU, HO>;

/** @brief StateFeedback law{K}: the sizes of a gain matrix, on the measured state. */
template<size_t NU, size_t NX, typename T>
StateFeedback(const Matrix<NU, NX, T>&) -> StateFeedback<NX, NU, T>;

namespace detail {

/** @brief An estimator that runs a filter (a FixedGainObserver): the kinds LQG, LQGI and OFLQG. */
template<typename Est>
concept runs_filter = estimator_traits<Est>::present;

/** @brief A filter that estimates a disturbance too (ND > 0): the kind OFLQG. */
template<typename Est>
concept estimates_disturbance = estimator_traits<Est>::ND > 0;

/** @brief Every LQ design runs as StateFeedback: with members, or on the design's constants. */
template<size_t NX, size_t NU, typename U, size_t NI, size_t NY, typename Kalman, size_t NZ, size_t NRU, size_t NM, size_t NDU, bool HO, typename T>
struct lq_runtime<design::StateFeedbackResult<NX, NU, U, NI, NY, Kalman, NZ, NRU, NM, NDU, HO>, T> {
    using type = StateFeedback<NX, NU, T, NI, NY, runtime_estimator_t<Kalman, T>, NZ, HeldGains, NRU, NM, NDU, HO>; ///< StateFeedback
    template<typename Est, typename Store>
    using on_constants = StateFeedback<NX, NU, T, NI, NY, Est, NZ, Store, NRU, NM, NDU, HO>; ///< StateFeedback on the design's constants
};

} // namespace detail

} // namespace damp
