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
 * trajectory: control(r, x_ref, u_ref), step(r, x_ref, u_ref, y).
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
 * command is clamped into u_lim (the design's InputLimits; unbounded by default), and predict()
 * moves the filter with it, so the clamp is not taken for a disturbance; back_calculate(u_sat)
 * before it reports a stage that clamps after the law. update() clamps d̂ into disturbance_lim, a
 * member of the runtime (unbounded by default), in d's units: the most of a load the law leans
 * against, as integral_lim on LQI.
 *
 * A design with an InputDelay (NRU > 0) has its filter on [x; v(k−1); …; v(k−NRU); d]: predict
 * shifts the applied input into the registers. x_ref is the plant's; the registers' reference is
 * u_ref.
 *
 * @note Compare with MATLAB®'s lqg with an integrated-white-noise disturbance model (setoutdist / setindist).
 *
 * A design with a MeasuredDisturbance (NM > 0) has its filter take d as a known input after u and
 * its target hold r under [d̂; d], so d̂ estimates only the part of a load d does not measure:
 * control(r, d), step(r, d, y), control(r, x_ref, u_ref, d), step(r, x_ref, u_ref, d, y).
 *
 * @tparam NRU Input registers of the filter's model (0: none)
 * @tparam NM  Measured-disturbance inputs (0: none)
 */
template<size_t NX, size_t NU, size_t NY, size_t NZ, typename T, typename Est, typename Store = detail::HeldGains, size_t NRU = 0, size_t NM = 0>
struct OFLQG : private detail::StateFeedbackCore<NX, NU, T, 0, NY, Est, NZ, Store, NRU, NM> {
    using Core = detail::StateFeedbackCore<NX, NU, T, 0, NY, Est, NZ, Store, NRU, NM>; ///< The shared implementation

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

    using Core::demand_lim;
    using Core::disturbance_lim;
    using Core::estimator;
    using Core::K;
    using Core::target;
    using Core::u_lim;

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

    /** @brief One tick: predict(), update(@p y), then control(@p r, @p x_ref, @p u_ref), the feedforward form. */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NZ, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NY, T>& y)
        requires(NM == 0)
    {
        return Core::step(r, x_ref, u_ref, y);
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

    /** @brief One tick: predict(), update(@p y), then control(@p r, @p x_ref, @p u_ref, @p d). */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NZ, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NM, T>& d, const ColVec<NY, T>& y)
        requires(NM > 0)
    {
        return Core::step(r, x_ref, u_ref, d, y);
    }

    using Core::applied;
    using Core::back_calculate;
    using Core::disturbance_estimate;
    using Core::getK;
    using Core::last_command;
    using Core::predict;
    using Core::reset;
    using Core::state_estimate;
    using Core::valid;
};

/** @brief OFLQG law{res}: the sizes, scalar, filter and target of an OFLQG design's result. */
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, typename Kalman, size_t NZ, size_t NRU, size_t NM>
OFLQG(const design::StateFeedbackResult<NX, NU, T, NI, NY, Kalman, NZ, NRU, NM>&) -> OFLQG<NX, NU, NY, NZ, T, detail::runtime_estimator_t<Kalman, T>, detail::HeldGains, NRU, NM>;

namespace detail {

/** @brief An LQ design with a disturbance model runs as OFLQG. */
template<size_t NX, size_t NU, typename U, size_t NY, typename Kalman, size_t NZ, size_t NRU, size_t NM, typename T>
    requires(design::StateFeedbackResult<NX, NU, U, 0, NY, Kalman, NZ, NRU, NM>::ND > 0)
struct lq_runtime<design::StateFeedbackResult<NX, NU, U, 0, NY, Kalman, NZ, NRU, NM>, T> {
    using type = OFLQG<NX, NU, NY, NZ, T, runtime_estimator_t<Kalman, T>, HeldGains, NRU, NM>; ///< OFLQG
    template<typename Est, typename Store>
    using on_constants = OFLQG<NX, NU, NY, NZ, T, Est, Store, NRU, NM>; ///< OFLQG on the design's constants
};

} // namespace detail

} // namespace damp
