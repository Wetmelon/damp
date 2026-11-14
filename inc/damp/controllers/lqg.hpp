// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lqg.hpp
 * @brief The linear-quadratic Gaussian law: design::lqg, LQGResult, the LQG runtime
 *
 * design::lqg(plant, [Ts, [InputDelay<N>{θ},]] {.Q, .R, .N}, KalmanNoise{.process, .measurement}[, TrackedOutputs<NZ>{{i, …}}])
 * → LQGResult → LQG. With an InputDelay the filter carries the inputs still on their way.
 * step(y) is predict, update(y), control. The low-latency tick is update(y), control(),
 * back_calculate, predict. With tracked outputs: control(r), step(r, y). With
 * feedforward from a planned trajectory: control([r,] x_ref, u_ref), step([r,] x_ref, u_ref, y).
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/lq_law.hpp"
#include "damp/controllers/detail/state_feedback_design.hpp"
#include "damp/controllers/state_feedback_core.hpp"
#include "damp/controllers/state_feedback_types.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/systems/discretization.hpp"
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
    const LQRConfig<NX, NU, 0, T>&           config,
    const KalmanNoise<NW, NV, T>&            noise
) {
    return detail::state_feedback(sys, config, detail::StateFeedbackOptions{.noise = noise});
}

/** @brief The sampled-data linear-quadratic Gaussian law (the process noise a density, the measurement a variance per reading): runs as LQG. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr auto lqg(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    damp::type_identity_t<T>                 Ts,
    const LQRConfig<NX, NU, 0, T>&           config,
    const KalmanNoise<NW, NV, T>&            noise
) {
    return detail::state_feedback(sys, Ts, config, detail::StateFeedbackOptions{.noise = noise});
}

/**
 * @brief The linear-quadratic Gaussian law holding the listed outputs at r through the model's
 *        steady-state target (x* = Nₓ·r, u* = Nᵤ·r, see lqr(sys, config, outputs)): runs as LQG, step(r, y).
 *
 * As for lqr: refused with DesignStatus::singular when no unique steady state holds r (NZ = NU)
 * or none holds every r (NZ < NU); with more tracked outputs than inputs (NZ > NU, as
 * TrackedOutputs{} gives on a plant with NY > NU) the target is the least-squares compromise,
 * which holds none of them exactly.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NZ>
[[nodiscard]] constexpr auto lqg(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const LQRConfig<NX, NU, 0, T>&           config,
    const KalmanNoise<NW, NV, T>&            noise,
    const TrackedOutputs<NZ>&                outputs
) {
    return detail::state_feedback(sys, config, detail::StateFeedbackOptions{.noise = noise, .track = outputs});
}

/** @brief The sampled-data linear-quadratic Gaussian law holding the listed outputs at r: runs as LQG, step(r, y). */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NZ>
[[nodiscard]] constexpr auto lqg(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    damp::type_identity_t<T>                 Ts,
    const LQRConfig<NX, NU, 0, T>&           config,
    const KalmanNoise<NW, NV, T>&            noise,
    const TrackedOutputs<NZ>&                outputs
) {
    return detail::state_feedback(sys, Ts, config, detail::StateFeedbackOptions{.noise = noise, .track = outputs});
}

/**
 * @brief The sampled-data linear-quadratic Gaussian law of a plant whose input arrives @p delay
 *        late: runs as LQG, its filter carrying the inputs still on their way.
 *
 * The weights and the noise are the plant's, as without the delay: the gain is lqr(sys, Ts,
 * delay, config)'s and the filter design::kalman(sys, Ts, delay, noise)'s, on
 * discretize(sys, Ts, delay) (detail::state_feedback with an InputDelay).
 *
 * @code
 * constexpr auto res = design::lqg(plant, 0.01, InputDelay<3>{0.022}, {.Q = Q, .R = R}, {.process = Qn, .measurement = Rn});
 * static constinit Law<res, float> law;
 * // each tick: law.update(y); u = law.control(); write(sat(u)); law.back_calculate(sat(u)); law.predict();
 * @endcode
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t N, typename S>
[[nodiscard]] constexpr auto lqg(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    damp::type_identity_t<T>                 Ts,
    const InputDelay<N, S>&                  delay,
    const LQRConfig<NX, NU, 0, T>&           config,
    const KalmanNoise<NW, NV, T>&            noise
) {
    return detail::state_feedback(sys, Ts, delay, config, detail::StateFeedbackOptions{.noise = noise});
}

/** @brief The sampled-data linear-quadratic Gaussian law of a plant whose input arrives @p delay late, holding the listed outputs at r: runs as LQG, step(r, y). */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t N, typename S, size_t NZ>
[[nodiscard]] constexpr auto lqg(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    damp::type_identity_t<T>                 Ts,
    const InputDelay<N, S>&                  delay,
    const LQRConfig<NX, NU, 0, T>&           config,
    const KalmanNoise<NW, NV, T>&            noise,
    const TrackedOutputs<NZ>&                outputs
) {
    return detail::state_feedback(sys, Ts, delay, config, detail::StateFeedbackOptions{.noise = noise, .track = outputs});
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
    requires damp::is_same_v<Reg, StateFeedbackResult<NX, NU, T>> || damp::is_same_v<Reg, LQRResult<NX, NU, T>>
[[nodiscard]] constexpr StateFeedbackResult<NX, NU, T, 0, NY, KalmanResult<NX, NU, NY, T, NW, NV>> lqg_from_parts(
    const KalmanResult<NX, NU, NY, T, NW, NV>& kest,
    const Reg&                                 reg
) {

    auto out = StateFeedbackResult<NX, NU, T, 0, NY, KalmanResult<NX, NU, NY, T, NW, NV>>{reg.K, reg.S, reg.e, false, DesignStatus::not_converged, kest.sys.Ts, kest};
    out.u_lim = reg.u_lim;
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
 * Runs design::lqg. The estimator is a FixedGainObserver with design::kalman's fixed gain
 * (`estimator`). The tick is update(y) then control(), or step(y): the regulator. A design with
 * TrackedOutputs holds those outputs at r through its steady-state target: control(r),
 * step(r, y). Feedforward from a planned trajectory (x_ref, u_ref) adds to it:
 * control([r,] x_ref, u_ref), step([r,] x_ref, u_ref, y). The command is clamped into u_lim
 * (the design's, LQRConfig; unbounded by default), and predict() moves the observer with it,
 * or with the input back_calculate(u_sat) recorded for a stage that clamps after the law.
 *
 * A design with an InputDelay (NRU > 0) has its filter on [x; v(k−1); …; v(k−NRU)]: predict
 * shifts the applied input into the registers. x_ref is the plant's; the registers' reference is
 * u_ref.
 *
 * @note Compare with MATLAB®'s lqg.
 *
 * @tparam NZ Tracked outputs: the size of r (0: none, the regulator)
 * @tparam NRU Input registers of the filter's model (0: none)
 */
template<size_t NX, size_t NU, size_t NY, typename T, typename Est, size_t NZ = 0, typename Store = detail::HeldGains, size_t NRU = 0>
struct LQG : private detail::StateFeedbackCore<NX, NU, T, 0, NY, Est, NZ, Store, NRU> {
    using Core = detail::StateFeedbackCore<NX, NU, T, 0, NY, Est, NZ, Store, NRU>; ///< The shared implementation

    constexpr LQG() = default;

    /** @brief From design::lqg or lqg_from_parts (any scalar). */
    template<size_t NXR, size_t NUR, typename U, size_t NIR, size_t NYR, typename Kalman, size_t NZR, size_t NRR>
        requires(!Store::constant)
    constexpr LQG(const design::StateFeedbackResult<NXR, NUR, U, NIR, NYR, Kalman, NZR, NRR>& result) { // NOLINT
        constexpr bool kind = NIR == 0 && design::StateFeedbackResult<NXR, NUR, U, NIR, NYR, Kalman, NZR, NRR>::ND == 0 && !damp::is_same_v<Kalman, NoEstimator>;
        static_assert(kind, "LQG runs design::lqg: this design has no Kalman filter (LQR, LQI), integral action (LQGI) or a disturbance model (OFLQG).");
        if constexpr (kind) {
            static_cast<Core&>(*this) = Core(result);
        }
    }

    using Core::estimator;
    using Core::K;
    using Core::u_lim;

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
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, typename Kalman, size_t NZ, size_t NRU>
LQG(const design::StateFeedbackResult<NX, NU, T, NI, NY, Kalman, NZ, NRU>&) -> LQG<NX, NU, NY, T, detail::runtime_estimator_t<Kalman, T>, NZ, detail::HeldGains, NRU>;

namespace detail {

/** @brief An LQ design with a filter and no integrators or disturbance model runs as LQG (with its tracked outputs, NZ > 0). */
template<size_t NX, size_t NU, typename U, size_t NY, typename Kalman, size_t NZ, size_t NRU, typename T>
    requires(!damp::is_same_v<Kalman, NoEstimator> && design::StateFeedbackResult<NX, NU, U, 0, NY, Kalman, NZ, NRU>::ND == 0)
struct lq_runtime<design::StateFeedbackResult<NX, NU, U, 0, NY, Kalman, NZ, NRU>, T> {
    using type = LQG<NX, NU, NY, T, runtime_estimator_t<Kalman, T>, NZ, HeldGains, NRU>; ///< LQG
    template<typename Est, typename Store>
    using on_constants = LQG<NX, NU, NY, T, Est, NZ, Store, NRU>; ///< LQG on the design's constants
};

} // namespace detail

} // namespace damp
