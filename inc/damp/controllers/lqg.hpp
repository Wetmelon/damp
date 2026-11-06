// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lqg.hpp
 * @brief The linear-quadratic Gaussian law: design::lqg, LQGResult, the LQG runtime
 *
 * design::lqg(plant, [Ts,] {.Q, .R, .N}, KalmanNoise{.process, .measurement}[, TrackedOutputs<NZ>{{i, …}}])
 * → LQGResult → LQG.
 * step(y) is predict, update(y), control. The low-latency tick is update(y), control(),
 * back_calculate, predict. With tracked outputs: control(r), step(r, y). With
 * feedforward from a planned trajectory: control([r,] x_ref, u_ref), step([r,] x_ref, u_ref, y).
 */

#include <cstddef>
#include <type_traits>

#include "damp/controllers/state_feedback_core.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace design {

/** @brief design::lqg's result: the gain and the steady-state Kalman filter. */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = NX, size_t NV = NY>
using LQGResult = StateFeedbackResult<NX, NU, T, 0, NY, KalmanResult<NX, NU, NY, T, NW, NV, 0>, 0>;

/**
 * @brief The linear-quadratic Gaussian law: LQR on a steady-state Kalman filter's estimate (design::kalman(sys, noise)): runs as LQG.
 * @note Compare with MATLAB®'s lqg(sys, QXU, QWV).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr auto lqg(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const LQRWeights<NX, NU, T>&             weights,
    const KalmanNoise<NW, NV, T>&            noise
) {
    return detail::state_feedback(sys, weights, detail::StateFeedbackOptions{.noise = noise});
}

/** @brief The sampled-data linear-quadratic Gaussian law (the process noise a density, the measurement a variance per reading): runs as LQG. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr auto lqg(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    std::type_identity_t<T>                  Ts,
    const LQRWeights<NX, NU, T>&             weights,
    const KalmanNoise<NW, NV, T>&            noise
) {
    return detail::state_feedback(sys, Ts, weights, detail::StateFeedbackOptions{.noise = noise});
}

/**
 * @brief The linear-quadratic Gaussian law holding the listed outputs at r through the model's
 *        steady-state target (x* = Nₓ·r, u* = Nᵤ·r, see lqr(sys, weights, outputs)): runs as LQG, step(r, y).
 *
 * As for lqr: refused with DesignStatus::singular when no unique steady state holds r (NZ = NU)
 * or none holds every r (NZ < NU); with more tracked outputs than inputs (NZ > NU, as
 * TrackedOutputs{} gives on a plant with NY > NU) the target is the least-squares compromise,
 * which holds none of them exactly.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NZ>
[[nodiscard]] constexpr auto lqg(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const LQRWeights<NX, NU, T>&             weights,
    const KalmanNoise<NW, NV, T>&            noise,
    const TrackedOutputs<NZ>&                outputs
) {
    return detail::state_feedback(sys, weights, detail::StateFeedbackOptions{.noise = noise, .track = outputs});
}

/** @brief The sampled-data linear-quadratic Gaussian law holding the listed outputs at r: runs as LQG, step(r, y). */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NZ>
[[nodiscard]] constexpr auto lqg(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    std::type_identity_t<T>                  Ts,
    const LQRWeights<NX, NU, T>&             weights,
    const KalmanNoise<NW, NV, T>&            noise,
    const TrackedOutputs<NZ>&                outputs
) {
    return detail::state_feedback(sys, Ts, weights, detail::StateFeedbackOptions{.noise = noise, .track = outputs});
}

/**
 * @brief Assemble an LQG design from separately computed Kalman and LQR designs
 *
 * The regulator must be discrete: a continuous one (Ts = 0) is refused with
 * DesignStatus::invalid_sample_time. A design on the plant (design::lqr) carries its sample time,
 * which must match kest.sys.Ts; a gain on bare matrices (discrete_lqr, matlab::dlqr) is discrete
 * with no sample time (Ts = 0), so the law takes the filter's. A failure already
 * reported by either half is kept (the regulator's first).
 *
 * @note Compare with MATLAB®'s lqgreg(kest, k) — exposed as the lqgreg() alias in matlab.hpp.
 *
 * @param kest  Steady-state Kalman design (design::kalman, no options)
 * @param reg   The regulator on the measured state: the LQRResult of design::lqr, discrete_lqr or
 *              matlab::dlqr
 * @return The LQG design: deploy with `LQG law{res}`
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Reg>
    requires std::is_same_v<Reg, StateFeedbackResult<NX, NU, T>> || std::is_same_v<Reg, LQRResult<NX, NU, T>>
[[nodiscard]] constexpr StateFeedbackResult<NX, NU, T, 0, NY, KalmanResult<NX, NU, NY, T, NW, NV>> lqg_from_parts(
    const KalmanResult<NX, NU, NY, T, NW, NV>& kest,
    const Reg&                                 reg
) {

    auto out = StateFeedbackResult<NX, NU, T, 0, NY, KalmanResult<NX, NU, NY, T, NW, NV>>{reg.K, reg.S, reg.e, false, DesignStatus::not_converged, kest.sys.Ts, kest};
    auto status = reg.success ? kest.status : reg.status;

    if (status == DesignStatus::ok && (!::damp::detail::discrete_gain(reg) || (reg.Ts != T{0} && reg.Ts != kest.sys.Ts))) {
        status = DesignStatus::invalid_sample_time; // a continuous regulator, or the halves disagree on Ts
    }
    return ::damp::detail::set_status(out, status);
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Linear-quadratic Gaussian law: LQR on the Kalman filter's estimate
 *
 * Runs design::lqg. The estimator is a Luenberger observer with design::kalman's fixed gain
 * (`estimator`). The tick is update(y) then control(), or step(y): the regulator. A design with
 * TrackedOutputs holds those outputs at r through its steady-state target: control(r),
 * step(r, y). Feedforward from a planned trajectory (x_ref, u_ref) adds to it:
 * control([r,] x_ref, u_ref), step([r,] x_ref, u_ref, y). Then back_calculate(u_sat) and
 * predict() move the observer with the input the plant received.
 *
 * @note Compare with MATLAB®'s lqg.
 *
 * @tparam NZ Tracked outputs: the size of r (0: none, the regulator)
 */
template<size_t NX, size_t NU, size_t NY, typename T, typename Est, size_t NZ = 0>
struct LQG : private detail::StateFeedbackCore<NX, NU, T, 0, NY, Est, NZ> {
    using Core = detail::StateFeedbackCore<NX, NU, T, 0, NY, Est, NZ>; ///< The shared implementation

    constexpr LQG() = default;

    /** @brief From design::lqg or lqg_from_parts (any scalar). */
    template<size_t NXR, size_t NUR, typename U, size_t NIR, size_t NYR, typename Kalman, size_t NZR>
    constexpr LQG(const design::StateFeedbackResult<NXR, NUR, U, NIR, NYR, Kalman, NZR>& result) { // NOLINT
        constexpr bool kind = NIR == 0 && design::StateFeedbackResult<NXR, NUR, U, NIR, NYR, Kalman, NZR>::ND == 0 && !std::is_same_v<Kalman, NoEstimator>;
        static_assert(kind, "LQG runs design::lqg: this design has no Kalman filter (LQR, LQI), integral action (LQGI) or a disturbance model (OFLQG).");
        if constexpr (kind) {
            static_cast<Core&>(*this) = Core(result);
        }
    }

    using Core::estimator;
    using Core::K;

    /** @brief The filter's measurement update with @p y. */
    constexpr void update(const ColVec<NY, T>& y) { Core::update(y); }

    /** @brief u = −K·x̂: the regulator. */
    [[nodiscard]] constexpr ColVec<NU, T> control()
        requires(NZ == 0)
    {
        return Core::control(ColVec<0, T>{}, ColVec<NX, T>{}, ColVec<NU, T>{});
    }

    /** @brief u = u* − K·(x̂ − x*): feedforward from a planned trajectory (x*, u*). */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref)
        requires(NZ == 0)
    {
        return Core::control(ColVec<0, T>{}, x_ref, u_ref);
    }

    /** @brief u = Nᵤ·r − K·(x̂ − Nₓ·r): the tracked outputs at @p r. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NZ, T>& r)
        requires(NZ > 0)
    {
        return Core::control(r, ColVec<NX, T>{}, ColVec<NU, T>{});
    }

    /** @brief control(@p r) about the planned trajectory (x*, u*), which adds to the target's. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NZ, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref)
        requires(NZ > 0)
    {
        return Core::control(r, x_ref, u_ref);
    }

    /** @brief One tick: predict(), update(@p y), then control(): the regulator. */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NY, T>& y)
        requires(NZ == 0)
    {
        return Core::step(ColVec<0, T>{}, ColVec<NX, T>{}, ColVec<NU, T>{}, y);
    }

    /** @brief One tick: predict(), update(@p y), then control(@p x_ref, @p u_ref), the feedforward form. */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NY, T>& y)
        requires(NZ == 0)
    {
        return Core::step(ColVec<0, T>{}, x_ref, u_ref, y);
    }

    /** @brief One tick: predict(), update(@p y), then control(@p r). */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NZ, T>& r, const ColVec<NY, T>& y)
        requires(NZ > 0)
    {
        return Core::step(r, ColVec<NX, T>{}, ColVec<NU, T>{}, y);
    }

    /** @brief One tick: predict(), update(@p y), then control(@p r, @p x_ref, @p u_ref), the feedforward form. */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NZ, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NY, T>& y)
        requires(NZ > 0)
    {
        return Core::step(r, x_ref, u_ref, y);
    }

    using Core::applied;
    using Core::back_calculate;
    using Core::getK;
    using Core::last_command;
    using Core::predict;
    using Core::reset;
    using Core::state_estimate;
};

/** @brief LQG law{res}: the sizes, scalar and filter of design::lqg's result. */
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, typename Kalman, size_t NZ>
LQG(const design::StateFeedbackResult<NX, NU, T, NI, NY, Kalman, NZ>&) -> LQG<NX, NU, NY, T, detail::runtime_estimator_t<Kalman, T>, NZ>;


namespace detail {

/** @brief An LQ design with a filter and no integrators or disturbance model runs as LQG (with its tracked outputs, NZ > 0). */
template<size_t NX, size_t NU, typename U, size_t NY, typename Kalman, size_t NZ, typename T>
    requires(!std::is_same_v<Kalman, NoEstimator> && design::StateFeedbackResult<NX, NU, U, 0, NY, Kalman, NZ>::ND == 0)
struct lq_runtime<design::StateFeedbackResult<NX, NU, U, 0, NY, Kalman, NZ>, T> {
    using type = LQG<NX, NU, NY, T, runtime_estimator_t<Kalman, T>, NZ>; ///< LQG
};

} // namespace detail

} // namespace damp
