// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lqr.hpp
 * @brief The linear-quadratic regulator: design::lqr, LQRResult, the LQR runtime
 *
 * design::lqr(plant, [Ts,] {.Q, .R, .N}) → LQRResult (also discrete_lqr's, riccati.hpp) → LQR;
 * with TrackedOutputs<NZ>{{i, …}}, a StateFeedbackResult. step(x) is predict, update(x), control
 * (predict is empty on LQR). The low-latency tick is update(x), control(). With tracked
 * outputs: control(r), step(r, x). With feedforward from a planned trajectory: control([r,]
 * x_ref, u_ref), step([r,] x_ref, u_ref, x).
 */

#include "damp/backend.hpp"
#include "damp/controllers/state_feedback_core.hpp" // IWYU pragma: export
#include "damp/design/riccati.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/core.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace design {
namespace detail {

/** @brief design::lqr's result: the gain, the Riccati solution, the closed-loop poles, the status and Ts. */
template<size_t NX, size_t NU, typename T, size_t NY>
[[nodiscard]] constexpr LQRResult<NX, NU, T> lqr_result(const StateFeedbackResult<NX, NU, T, 0, NY>& r) {
    return {r.K, r.S, r.e, r.success, r.status, r.Ts, r.Ts > T{0}};
}

} // namespace detail

/**
 * @brief The linear-quadratic regulator on the plant (DARE on a discrete plant, CARE on a continuous one): runs as LQR.
 * @note Compare with MATLAB®'s lqr(sys, Q, R, N).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr LQRResult<NX, NU, T> lqr(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const LQRWeights<NX, NU, 0, T>& weights) {
    return detail::lqr_result(detail::state_feedback(sys, weights, detail::StateFeedbackOptions{}));
}

/**
 * @brief The sampled-data linear-quadratic regulator: a continuous plant and its cost, the law at @p Ts (ZOH, Van Loan): runs as LQR.
 * @note Compare with MATLAB®'s lqrd(sys, Q, R, Ts).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr LQRResult<NX, NU, T> lqr(const StateSpace<NX, NU, NY, T, NW, NV>& sys, damp::type_identity_t<T> Ts, const LQRWeights<NX, NU, 0, T>& weights) {
    return detail::lqr_result(detail::state_feedback(sys, Ts, weights, detail::StateFeedbackOptions{}));
}

/**
 * @brief The linear-quadratic regulator holding the listed outputs at r through the model's
 *        steady-state target (x* = Nₓ·r, u* = Nᵤ·r): runs as LQR, step(r, x).
 *
 * The target is design::steady_state_target's on the listed outputs (TrackedOutputs{}: all of
 * them). A design is refused with DesignStatus::singular when the steady state is not unique
 * (NZ = NU) or none holds every r (NZ < NU). With more tracked outputs than inputs (NZ > NU,
 * as TrackedOutputs{} gives on a plant with NY > NU) the target is the least-squares
 * compromise, which holds none of them exactly. A model error or a constant disturbance
 * leaves an offset: lqi removes it.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NZ>
[[nodiscard]] constexpr auto lqr(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const LQRWeights<NX, NU, 0, T>& weights, const TrackedOutputs<NZ>& outputs) {
    return detail::state_feedback(sys, weights, detail::StateFeedbackOptions{.track = outputs});
}

/** @brief The sampled-data linear-quadratic regulator holding the listed outputs at r (see lqr(sys, weights, outputs)): runs as LQR, step(r, x). */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NZ>
[[nodiscard]] constexpr auto lqr(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    damp::type_identity_t<T>                 Ts,
    const LQRWeights<NX, NU, 0, T>&          weights,
    const TrackedOutputs<NZ>&                outputs
) {
    return detail::state_feedback(sys, Ts, weights, detail::StateFeedbackOptions{.track = outputs});
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Linear-quadratic regulator on the measured state: u = u* − K·(x − x*)
 *
 * Runs design::lqr, or a gain design (design::place, discrete_lqr): `LQR law{res.as<float>()}`, or
 * `Law<kRes, float>` for a struct member. step(x) is predict, update(x), control; predict is
 * empty, so that equals update then control. The regulator is x* = 0 and u* = 0. A design with
 * TrackedOutputs holds those outputs at r through its steady-state target: control(r),
 * step(r, x). Feedforward from a planned trajectory
 * (x_ref, u_ref) adds to it: control([r,] x_ref, u_ref), step([r,] x_ref, u_ref, x). A failed or
 * continuous design does not compile when constant-evaluated.
 *
 * @code
 * static constinit LQR law{kRes.as<float>()};
 * law.update(x);
 * u = law.control(); // the regulator
 * @endcode
 *
 * @note Compare with MATLAB®'s lqr / lqrd.
 *
 * @tparam NZ Tracked outputs: the size of r (0: none, the regulator)
 */
template<size_t NX, size_t NU, typename T = float, size_t NY = NX, size_t NZ = 0, typename Store = detail::HeldGains>
struct LQR : private detail::StateFeedbackCore<NX, NU, T, 0, NY, NoEstimator, NZ, Store> {
    using Core = detail::StateFeedbackCore<NX, NU, T, 0, NY, NoEstimator, NZ, Store>; ///< The shared implementation

    constexpr LQR() = default;

    /** @brief From a gain K. */
    constexpr explicit LQR(const Matrix<NU, NX, T>& K_)
        requires(NZ == 0 && !Store::constant)
        : Core(K_) {}

    /** @brief From a gain design (design::place's PlaceResult, LQRResult; any scalar). */
    template<typename R>
        requires(detail::gain_result<R> && NZ == 0 && !Store::constant)
    constexpr LQR(const R& result) : Core(result) {} // NOLINT

    /** @brief From design::lqr (any scalar). */
    template<size_t NXR, size_t NUR, typename U, size_t NIR, size_t NYR, typename Kalman, size_t NZR>
        requires(!Store::constant)
    constexpr LQR(const design::StateFeedbackResult<NXR, NUR, U, NIR, NYR, Kalman, NZR>& result) { // NOLINT
        constexpr bool kind = NIR == 0 && damp::is_same_v<Kalman, NoEstimator>;
        static_assert(kind, "LQR runs design::lqr: this design has integral action (LQI, LQGI) or a Kalman filter (LQG, OFLQG).");
        if constexpr (kind) {
            static_cast<Core&>(*this) = Core(result);
        }
    }

    using Core::K;

    /** @brief Hold the measured state @p x for the next control. */
    constexpr void update(const ColVec<NX, T>& x) { Core::update(x); }

    /** @brief u = −K·x at the state of the last update: the regulator. */
    [[nodiscard]] constexpr ColVec<NU, T> control()
        requires(NZ == 0)
    {
        return Core::control(ColVec<0, T>{}, ColVec<NX, T>{}, ColVec<NU, T>{});
    }

    /** @brief u = u* − K·(x − x*) at the state of the last update: feedforward from a planned trajectory (x*, u*). */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref)
        requires(NZ == 0)
    {
        return Core::control(ColVec<0, T>{}, x_ref, u_ref);
    }

    /** @brief u = Nᵤ·r − K·(x − Nₓ·r) at the state of the last update: the tracked outputs at @p r. */
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

    /** @brief One tick: predict(), update(@p x), then control(): the regulator. */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NX, T>& x)
        requires(NZ == 0)
    {
        return Core::step(ColVec<0, T>{}, ColVec<NX, T>{}, ColVec<NU, T>{}, x);
    }

    /** @brief One tick: predict(), update(@p x), then control(@p x_ref, @p u_ref), the feedforward form. */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NX, T>& x)
        requires(NZ == 0)
    {
        return Core::step(ColVec<0, T>{}, x_ref, u_ref, x);
    }

    /** @brief One tick: predict(), update(@p x), then control(@p r). */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NZ, T>& r, const ColVec<NX, T>& x)
        requires(NZ > 0)
    {
        return Core::step(r, ColVec<NX, T>{}, ColVec<NU, T>{}, x);
    }

    /** @brief One tick: predict(), update(@p x), then control(@p r, @p x_ref, @p u_ref), the feedforward form. */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NZ, T>& r, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref, const ColVec<NX, T>& x)
        requires(NZ > 0)
    {
        return Core::step(r, x_ref, u_ref, x);
    }

    using Core::back_calculate;
    using Core::getK;
    using Core::predict;
    using Core::reset;
};

/** @brief LQR law{res}: the sizes and scalar of design::lqr's result. */
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, typename Kalman, size_t NZ>
LQR(const design::StateFeedbackResult<NX, NU, T, NI, NY, Kalman, NZ>&) -> LQR<NX, NU, T, NY, NZ>;

/** @brief LQR law{K}: the sizes of a gain matrix. */
template<size_t NU, size_t NX, typename T>
LQR(const Matrix<NU, NX, T>&) -> LQR<NX, NU, T>;

/** @brief LQR law{res}: the sizes of a gain design (design::place, LQRResult). */
template<typename R>
    requires detail::gain_result<R>
LQR(const R& r) -> LQR<damp::remove_cvref_t<decltype(r.K)>::cols(), damp::remove_cvref_t<decltype(r.K)>::rows(), typename damp::remove_cvref_t<decltype(r.K)>::value_type>;

namespace detail {

/** @brief A gain design (LQRResult, design::place's PlaceResult) runs as LQR. */
template<typename R, typename T>
    requires gain_result<R>
struct lq_runtime<R, T> {
    using Gain = damp::remove_cvref_t<decltype(damp::declval<const R&>().K)>;
    using type = LQR<Gain::cols(), Gain::rows(), T>; ///< LQR
    template<typename, typename Store>
    using on_constants = LQR<Gain::cols(), Gain::rows(), T, Gain::cols(), 0, Store>; ///< LQR on the design's constants
};

/** @brief An LQ design with no integrators or filter runs as LQR (with its tracked outputs, NZ > 0). */
template<size_t NX, size_t NU, typename U, size_t NY, size_t NZ, typename T>
struct lq_runtime<design::StateFeedbackResult<NX, NU, U, 0, NY, NoEstimator, NZ>, T> {
    using type = LQR<NX, NU, T, NY, NZ>; ///< LQR
    template<typename, typename Store>
    using on_constants = LQR<NX, NU, T, NY, NZ, Store>; ///< LQR on the design's constants
};

} // namespace detail

} // namespace damp
