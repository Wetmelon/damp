// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lqgi.hpp
 * @brief The linear-quadratic Gaussian law with integral action: design::lqgi, LQGIResult, the LQGI runtime
 *
 * design::lqgi(plant, [Ts, [InputDelay<N>{θ},]] {.Q, .R, .N, .Q_integral}, KalmanNoise{.process, .measurement}[, IntegralAction<NI>{.outputs = {i, …}}])
 * → LQGIResult → LQGI. With an InputDelay the filter carries the inputs still on their way. step(r, y) is predict, update(y), control(r). The low-latency tick is
 * update(y), control(r), back_calculate, predict. With feedforward from a
 * planned trajectory: control(r, x_ref, u_ref), step(r, x_ref, u_ref, y).
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/lq_law.hpp"
#include "damp/controllers/detail/state_feedback_design.hpp"
#include "damp/controllers/state_feedback_core.hpp"
#include "damp/controllers/state_feedback_types.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace design {

/** @brief design::lqgi's result: the gain on [x; xi], the Kalman filter and the integral action. */
template<size_t NX, size_t NU, size_t NY, size_t NI, typename T = double, size_t NW = NX, size_t NV = NY>
using LQGIResult = StateFeedbackResult<NX, NU, T, NI, NY, KalmanResult<NX, NU, NY, T, NW, NV, 0>, 0>;

/**
 * @brief The linear-quadratic Gaussian law with integral action on every output: runs as LQGI.
 * @note Compare with MATLAB®'s lqgtrack.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr auto lqgi(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const LQRConfig<NX, NU, NY, T>&          config,
    const KalmanNoise<NW, NV, T>&            noise
) {
    return detail::state_feedback(sys, config, detail::StateFeedbackOptions{.noise = noise, .offset_free = IntegralAction{}});
}

/** @brief The sampled-data linear-quadratic Gaussian law with integral action on every output: runs as LQGI. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr auto lqgi(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    damp::type_identity_t<T>                 Ts,
    const LQRConfig<NX, NU, NY, T>&          config,
    const KalmanNoise<NW, NV, T>&            noise
) {
    return detail::state_feedback(sys, Ts, config, detail::StateFeedbackOptions{.noise = noise, .offset_free = IntegralAction{}});
}

/** @brief The linear-quadratic Gaussian law with integral action on the listed outputs: runs as LQGI. */
template<size_t NI, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
    requires(NI > 0)
[[nodiscard]] constexpr auto lqgi(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const LQRConfig<NX, NU, NI, T>&          config,
    const KalmanNoise<NW, NV, T>&            noise,
    const IntegralAction<NI>&                outputs
) {
    return detail::state_feedback(sys, config, detail::StateFeedbackOptions{.noise = noise, .offset_free = outputs});
}

/** @brief The sampled-data linear-quadratic Gaussian law with integral action on the listed outputs: runs as LQGI. */
template<size_t NI, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
    requires(NI > 0)
[[nodiscard]] constexpr auto lqgi(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    damp::type_identity_t<T>                 Ts,
    const LQRConfig<NX, NU, NI, T>&          config,
    const KalmanNoise<NW, NV, T>&            noise,
    const IntegralAction<NI>&                outputs
) {
    return detail::state_feedback(sys, Ts, config, detail::StateFeedbackOptions{.noise = noise, .offset_free = outputs});
}

/**
 * @brief The sampled-data linear-quadratic Gaussian law with integral action on every output, of
 *        a plant whose input arrives @p delay late: runs as LQGI, its filter carrying the inputs
 *        still on their way.
 *
 * The weights and the noise are the plant's, as without the delay: the gain is lqi(sys, Ts,
 * delay, config)'s and the filter design::kalman(sys, Ts, delay, noise)'s.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t N, typename S>
[[nodiscard]] constexpr auto lqgi(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    damp::type_identity_t<T>                 Ts,
    const InputDelay<N, S>&                  delay,
    const LQRConfig<NX, NU, NY, T>&          config,
    const KalmanNoise<NW, NV, T>&            noise
) {
    return detail::state_feedback(sys, Ts, delay, config, detail::StateFeedbackOptions{.noise = noise, .offset_free = IntegralAction{}});
}

/** @brief The sampled-data linear-quadratic Gaussian law with integral action on the listed outputs, of a plant whose input arrives @p delay late: runs as LQGI. */
template<size_t NI, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t N, typename S>
    requires(NI > 0)
[[nodiscard]] constexpr auto lqgi(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    damp::type_identity_t<T>                 Ts,
    const InputDelay<N, S>&                  delay,
    const LQRConfig<NX, NU, NI, T>&          config,
    const KalmanNoise<NW, NV, T>&            noise,
    const IntegralAction<NI>&                outputs
) {
    return detail::state_feedback(sys, Ts, delay, config, detail::StateFeedbackOptions{.noise = noise, .offset_free = outputs});
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Linear-quadratic Gaussian law with integral action on the measured outputs
 *
 * Runs design::lqgi. The estimator is a FixedGainObserver with design::kalman's fixed gain
 * (`estimator`). xi integrates r − Ci·y on the last update's measurement; the law is
 * u = u* − K·[x̂ − x*; xi]. The tick is update(y) then control(r), or step(r, y), with x* = 0 and
 * u* = 0; feedforward from a planned trajectory (x_ref, u_ref) is control(r, x_ref, u_ref), or
 * step(r, x_ref, u_ref, y). The command is clamped into u_lim and xi into integral.lim (the
 * design's, LQRConfig; unbounded by default): predict() moves the filter on the clamped command
 * and holds the integrators against the clamp. back_calculate(u_sat) before it reports a stage
 * that clamps after the law.
 *
 * A design with an InputDelay (NRU > 0) has its filter on [x; v(k−1); …; v(k−NRU)]: predict
 * shifts the applied input into the registers. x_ref is the plant's; the registers' reference is
 * u_ref.
 *
 * @note Compare with MATLAB®'s lqgtrack.
 *
 * @tparam NRU Input registers of the filter's model (0: none)
 */
template<size_t NX, size_t NU, size_t NY, size_t NI, typename T, typename Est, typename Store = detail::HeldGains, size_t NRU = 0>
struct LQGI : private detail::StateFeedbackCore<NX, NU, T, NI, NY, Est, 0, Store, NRU> {
    using Core = detail::StateFeedbackCore<NX, NU, T, NI, NY, Est, 0, Store, NRU>; ///< The shared implementation

    constexpr LQGI() = default;

    /** @brief From design::lqgi (any scalar). */
    template<size_t NXR, size_t NUR, typename U, size_t NIR, size_t NYR, typename Kalman, size_t NZR, size_t NRR>
        requires(!Store::constant)
    constexpr LQGI(const design::StateFeedbackResult<NXR, NUR, U, NIR, NYR, Kalman, NZR, NRR>& result) { // NOLINT
        constexpr bool kind = NIR > 0 && !damp::is_same_v<Kalman, NoEstimator>;
        static_assert(kind, "LQGI runs design::lqgi: this design has no Kalman filter (LQR, LQI) or no integral action (LQG, OFLQG).");
        if constexpr (kind) {
            static_cast<Core&>(*this) = Core(result);
        }
    }

    using Core::estimator;
    using Core::integral;
    using Core::K;
    using Core::u_lim;

    /** @brief The filter's measurement update with @p y (the integrators read it at the next control). */
    constexpr void update(const ColVec<NY, T>& y) { Core::update(y); }

    /** @brief u = −K·[x̂; xi]: the basic call, x* = 0 and u* = 0; holds r − Ci·y for the next predict. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NI, T>& r) { return Core::control(r, ColVec<NX, T>{}, ColVec<NU, T>{}); }

    /** @brief u = u* − K·[x̂ − x*; xi]: feedforward from a planned trajectory (x*, u*); holds r − Ci·y for the next predict. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NI, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref) {
        return Core::control(r, x_ref, u_ref);
    }

    /** @brief One tick: predict(), update(@p y), then control(@p r). The basic call. */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NI, T>& r, const ColVec<NY, T>& y) { return Core::step(r, ColVec<NX, T>{}, ColVec<NU, T>{}, y); }

    /** @brief One tick: predict(), update(@p y), then control(@p r, @p x_ref, @p u_ref), the feedforward form. */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NI, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NY, T>& y) {
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

/** @brief LQGI law{res}: the sizes, scalar and filter of design::lqgi's result. */
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, typename Kalman, size_t NZ, size_t NRU>
LQGI(const design::StateFeedbackResult<NX, NU, T, NI, NY, Kalman, NZ, NRU>&) -> LQGI<NX, NU, NY, NI, T, detail::runtime_estimator_t<Kalman, T>, detail::HeldGains, NRU>;

namespace detail {

/** @brief An LQ design with a filter and integrators runs as LQGI. */
template<size_t NX, size_t NU, typename U, size_t NI, size_t NY, typename Kalman, size_t NRU, typename T>
    requires(NI > 0 && !damp::is_same_v<Kalman, NoEstimator>)
struct lq_runtime<design::StateFeedbackResult<NX, NU, U, NI, NY, Kalman, 0, NRU>, T> {
    using type = LQGI<NX, NU, NY, NI, T, runtime_estimator_t<Kalman, T>, HeldGains, NRU>; ///< LQGI
    template<typename Est, typename Store>
    using on_constants = LQGI<NX, NU, NY, NI, T, Est, Store, NRU>; ///< LQGI on the design's constants
};

} // namespace detail

} // namespace damp
