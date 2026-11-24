// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file oflqg.hpp
 * @brief The offset-free linear-quadratic Gaussian law: design::state_feedback with a KalmanNoise and a disturbance spec, OFLQGResult, the OFLQG runtime
 *
 * design::state_feedback(plant, [Ts,] {.Q, .R, .N}, KalmanNoise{.process, .measurement, .disturbance}, disturbance[, InputDelay<N>{θ}])
 * (state_feedback.hpp) → OFLQGResult → OFLQG. With an InputDelay the filter carries the inputs still on their way. step(r, y) is predict, update(y), control(r). The low-latency tick is update(y),
 * control(r), back_calculate, predict. With feedforward from a planned
 * trajectory: control(r, x_ref, u_ref), step(r, x_ref, u_ref, y). While another source drives the
 * plant, track(u, r, y) takes step's place: d̂ learns the load under u.
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/lq_law.hpp" // IWYU pragma: export
#include "damp/controllers/detail/state_feedback_core.hpp"
#include "damp/controllers/law.hpp"                  // IWYU pragma: export
#include "damp/controllers/state_feedback_types.hpp" // IWYU pragma: export
#include "damp/design/disturbance.hpp"               // IWYU pragma: export
#include "damp/estimation/kalman.hpp"                // IWYU pragma: export
#include "damp/matrix/colvec.hpp"

namespace damp {

namespace design {

/** @brief an OFLQG design's result: the gain, the Kalman filter on [x; d] and the steady-state target. */
template<size_t NX, size_t NU, size_t NY, size_t ND, size_t NZ, typename T = double, size_t NW = NX, size_t NV = NY>
using OFLQGResult = StateFeedbackResult<NX, NU, T, 0, NY, KalmanResult<NX + ND, NU, NY, T, NW + ND, NV, ND>, NZ>;

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Offset-free linear-quadratic Gaussian law: a Kalman filter on [x; d] and the steady-state target
 *
 * Runs design::state_feedback with a KalmanNoise and a disturbance spec. The estimator is a FixedGainObserver on [x; d] with design::kalman's fixed
 * gain, holding deviations from the target (`estimator`). The target holds the controlled outputs
 * at r under the disturbance estimate d̂, taken about the reference point (x*, u*). The tick is
 * step(r, y), or update(y), control(r), back_calculate, predict() (predict moves the filter),
 * with x* = 0 and u* = 0; feedforward from a planned
 * trajectory (x_ref, u_ref) is control(r, x_ref, u_ref), or step(r, x_ref, u_ref, y). The
 * command is clamped into the design's InputLimits and its move into an InputRateLimits (unbounded
 * by default), and predict() moves the filter with it, so the clamp is not taken for a disturbance;
 * back_calculate(u_sat) before it reports a stage that clamps after the law. update() clamps d̂
 * into StateFeedbackLimits::disturbance (the runtime's alone; unbounded by default), in d's units:
 * the most of a load the law leans against, as the integrators' limits on LQI. limits() reads them and set_limits() moves them at run time.
 *
 * A design with an InputDelay (NRU > 0) has its filter on [x; v(k−1); …; v(k−NRU); d]: predict
 * shifts the applied input into the delay states. x_ref is the plant's; the delay states' reference is
 * u_ref.
 *
 * @note Compare with MATLAB®'s lqg with an integrated-white-noise disturbance model (setoutdist / setindist).
 *
 * A design with a MeasuredDisturbance (NM > 0) has its filter take d as a known input after u and
 * its target hold r under [d̂; d], so d̂ estimates only the part of a load d does not measure:
 * control(r, d), step(r, d, y), control(r, x_ref, u_ref, d), step(r, x_ref, u_ref, d, y).
 *
 * While another source drives the plant (an operator, another law), track(u, …) takes the place of
 * each step(…), its arguments step's with the input applied first: track(u, r, y),
 * track(u, r, x_ref, u_ref, y), track(u, r, d, y), track(u, r, x_ref, u_ref, d, y). The filter
 * predicts with u, so d̂ learns the load the other source holds against: at a steady state the
 * target's input is then u and the first step returns it. Away from one the law's command differs
 * by its transient, which a design::Handover fades by exp(−Ts/fade_time) each tick, within its
 * input limits.
 *
 * @tparam NRU Input delay states of the filter's model (0: none)
 * @tparam NM  Measured-disturbance inputs (0: none)
 * @tparam NDU Inputs with a move limit (an InputRateLimits: NU; 0: none)
 * @tparam HO  Whether the law fades the mismatch after a track (a Handover; false: none)
 */
template<size_t NX, size_t NU, size_t NY, size_t NZ, typename T, typename Est, typename Store = detail::HeldGains, size_t NRU = 0, size_t NM = 0, size_t NDU = 0, bool HO = false>
struct OFLQG : private detail::StateFeedbackCore<NX, NU, T, 0, NY, Est, NZ, Store, NRU, NM, NDU, HO> {
    using Core = detail::StateFeedbackCore<NX, NU, T, 0, NY, Est, NZ, Store, NRU, NM, NDU, HO>; ///< The shared implementation

    constexpr OFLQG() = default;

    /** @brief From design::state_feedback(any scalar). */
    template<detail::state_feedback_result R>
        requires(!Store::constant)
    constexpr OFLQG(const R& result) { // NOLINT
        constexpr bool kind = R::ND > 0;
        static_assert(kind, "OFLQG runs design::state_feedback with a KalmanNoise and a disturbance spec: this design has no disturbance model (LQR, LQI, LQG, LQGI).");
        if constexpr (kind) {
            static_cast<Core&>(*this) = Core(result);
        }
    }

    using Core::estimator;
    using Core::K;
    using Core::target;

    /** @brief The filter's measurement update with @p y. */
    constexpr void update(const ColVec<NY, T>& y) { Core::update(y); }

    /** @brief The target for @p r under d̂ about the origin, then the law about it: the basic call. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NZ, T>& r)
        requires(NM == 0)
    {
        return Core::control(r, ColVec<NX, T>{}, ColVec<NU, T>{});
    }

    /** @brief The target for @p r under d̂ about (@p x_ref, @p u_ref), then the law about it: feedforward from a planned trajectory. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NZ, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref)
        requires(NM == 0)
    {
        return Core::control(r, x_ref, u_ref);
    }

    /** @brief One tick: predict(), update(@p y), then control(@p r). The basic call. */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NZ, T>& r, const ColVec<NY, T>& y)
        requires(NM == 0)
    {
        return Core::step(r, ColVec<NX, T>{}, ColVec<NU, T>{}, y);
    }

    /** @brief In place of step(r, y) while another source drives the plant: its predict() and update(y), then the memory that starts the next command from @p u, the input applied. */
    constexpr void track(const ColVec<NU, T>& u, const ColVec<NZ, T>& r, const ColVec<NY, T>& y)
        requires(NM == 0)
    {
        Core::track(u, r, ColVec<NX, T>{}, ColVec<NU, T>{}, y);
    }

    /** @brief One tick: predict(), update(@p y), then control(@p r, @p x_ref, @p u_ref), the feedforward form. */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NZ, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NY, T>& y)
        requires(NM == 0)
    {
        return Core::step(r, x_ref, u_ref, y);
    }

    /** @brief In place of step(r, x_ref, u_ref, y) while another source drives the plant: its predict() and update(y), then the memory that starts the next command from @p u, the input applied. */
    constexpr void track(const ColVec<NU, T>& u, const ColVec<NZ, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NY, T>& y)
        requires(NM == 0)
    {
        Core::track(u, r, x_ref, u_ref, y);
    }

    /** @brief control(@p r) under the measured disturbance @p d: the target under it, after the reference. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NZ, T>& r, const ColVec<NM, T>& d)
        requires(NM > 0)
    {
        return Core::control(r, ColVec<NX, T>{}, ColVec<NU, T>{}, d);
    }

    /** @brief control(@p r, @p d) about the planned trajectory (@p x_ref, @p u_ref): the reference and its feedforward, then the measured disturbance. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NZ, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NM, T>& d)
        requires(NM > 0)
    {
        return Core::control(r, x_ref, u_ref, d);
    }

    /** @brief One tick: predict(), update(@p y), then control(@p r, @p d). */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NZ, T>& r, const ColVec<NM, T>& d, const ColVec<NY, T>& y)
        requires(NM > 0)
    {
        return Core::step(r, ColVec<NX, T>{}, ColVec<NU, T>{}, d, y);
    }

    /** @brief In place of step(r, d, y) while another source drives the plant: its predict() and update(y), then the memory that starts the next command from @p u, the input applied. */
    constexpr void track(const ColVec<NU, T>& u, const ColVec<NZ, T>& r, const ColVec<NM, T>& d, const ColVec<NY, T>& y)
        requires(NM > 0)
    {
        Core::track(u, r, ColVec<NX, T>{}, ColVec<NU, T>{}, d, y);
    }

    /** @brief One tick: predict(), update(@p y), then control(@p r, @p x_ref, @p u_ref, @p d). */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NZ, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NM, T>& d, const ColVec<NY, T>& y)
        requires(NM > 0)
    {
        return Core::step(r, x_ref, u_ref, d, y);
    }

    /** @brief In place of step(r, x_ref, u_ref, d, y) while another source drives the plant: its predict() and update(y), then the memory that starts the next command from @p u, the input applied. */
    constexpr void track(const ColVec<NU, T>& u, const ColVec<NZ, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NM, T>& d, const ColVec<NY, T>& y)
        requires(NM > 0)
    {
        Core::track(u, r, x_ref, u_ref, d, y);
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

/** @brief OFLQG law{res}: the sizes, scalar, filter and target of an OFLQG design's result. */
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, typename Kalman, size_t NZ, size_t NRU, size_t NM, size_t NDU, bool HO>
OFLQG(const design::StateFeedbackResult<NX, NU, T, NI, NY, Kalman, NZ, NRU, NM, NDU, HO>&) -> OFLQG<NX, NU, NY, NZ, T, detail::runtime_estimator_t<Kalman, T>, detail::HeldGains, NRU, NM, NDU, HO>;

namespace detail {

/** @brief An LQ design with a disturbance model runs as OFLQG. */
template<size_t NX, size_t NU, typename U, size_t NY, typename Kalman, size_t NZ, size_t NRU, size_t NM, size_t NDU, bool HO, typename T>
    requires(design::StateFeedbackResult<NX, NU, U, 0, NY, Kalman, NZ, NRU, NM, NDU, HO>::ND > 0)
struct lq_runtime<design::StateFeedbackResult<NX, NU, U, 0, NY, Kalman, NZ, NRU, NM, NDU, HO>, T> {
    using type = OFLQG<NX, NU, NY, NZ, T, runtime_estimator_t<Kalman, T>, HeldGains, NRU, NM, NDU, HO>; ///< OFLQG
    template<typename Est, typename Store>
    using on_constants = OFLQG<NX, NU, NY, NZ, T, Est, Store, NRU, NM, NDU, HO>; ///< OFLQG on the design's constants
};

} // namespace detail

} // namespace damp
