// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lqr.hpp
 * @brief The linear-quadratic regulator: design::state_feedback, LQRResult, the LQR runtime
 *
 * design::state_feedback(plant, [Ts,] {.Q, .R, .N}[, TrackedOutputs{{i, …}}][, InputDelay<N>{θ}])
 * (state_feedback.hpp) → StateFeedbackResult → LQR; discrete_lqr (riccati.hpp) → LQRResult → LQR.
 * With an InputDelay the LQR holds the inputs still on their way. step(x) is predict, update(x), control
 * (predict is empty on LQR without an InputDelay). The low-latency tick is update(x), control(),
 * then back_calculate and predict() with an InputDelay. With tracked
 * outputs: control(r), step(r, x). With feedforward from a planned trajectory: control([r,]
 * x_ref, u_ref), step([r,] x_ref, u_ref, x).
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/lq_law.hpp" // IWYU pragma: export
#include "damp/controllers/detail/state_feedback_core.hpp"
#include "damp/controllers/law.hpp"                  // IWYU pragma: export
#include "damp/controllers/state_feedback_types.hpp" // IWYU pragma: export
#include "damp/design/riccati.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/core.hpp"

namespace damp {
namespace design {
namespace detail {

/** @brief design::state_feedback's result: the gain, the Riccati solution, the closed-loop poles, the status and Ts. */
template<size_t NX, size_t NU, typename T, size_t NY>
[[nodiscard]] constexpr LQRResult<NX, NU, T> lqr_result(const StateFeedbackResult<NX, NU, T, 0, NY>& r) {
    return {r.K, r.S, r.e, r.success, r.status, r.Ts, r.Ts > T{0}, r.u_lim, r.demand_lim};
}

} // namespace detail

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Linear-quadratic regulator on the measured state: u = u* − K·(x − x*)
 *
 * Runs design::state_feedback, or a gain design (design::place, discrete_lqr): `LQR law{res.as<float>()}`, or
 * `Law<kRes, float>` for a struct member. step(x) is predict, update(x), control; predict is
 * empty without an InputDelay, so that equals update then control. The regulator is x* = 0 and u* = 0. A design with
 * TrackedOutputs holds those outputs at r through its steady-state target: control(r),
 * step(r, x). Feedforward from a planned trajectory
 * (x_ref, u_ref) adds to it: control([r,] x_ref, u_ref), step([r,] x_ref, u_ref, x). The command
 * is clamped into u_lim (design::state_feedback's LQRLimits; unbounded by default, and for a gain K or a
 * design::place result). A failed or continuous design does not compile when constant-evaluated.
 *
 * @code
 * static constinit LQR law{kRes.as<float>()};
 * law.update(x);
 * u = law.control(); // the regulator
 * @endcode
 *
 * @note Compare with MATLAB®'s lqr / lqrd.
 *
 * A design with an InputDelay (NRU > 0) runs on [x; v(k−1); …; v(k−NRU)]: update takes the plant
 * state, and predict shifts the applied input into the registers, so it is no longer empty.
 *
 * @tparam NZ Tracked outputs: the size of r (0: none, the regulator)
 * @tparam NRU Input registers the law holds (0: none)
 */
template<size_t NX, size_t NU, typename T = float, size_t NY = NX, size_t NZ = 0, typename Store = detail::HeldGains, size_t NRU = 0>
struct LQR : private detail::StateFeedbackCore<NX, NU, T, 0, NY, NoEstimator, NZ, Store, NRU> {
    using Core = detail::StateFeedbackCore<NX, NU, T, 0, NY, NoEstimator, NZ, Store, NRU>; ///< The shared implementation

    constexpr LQR() = default;

    /** @brief From a gain K. */
    constexpr explicit LQR(const Matrix<NU, NX, T>& K_)
        requires(NZ == 0 && NRU == 0 && !Store::constant)
        : Core(K_) {}

    /** @brief From a gain design (design::place's PlaceResult, LQRResult; any scalar). */
    template<typename R>
        requires(detail::gain_result<R> && NZ == 0 && NRU == 0 && !Store::constant)
    constexpr LQR(const R& result) : Core(result) {} // NOLINT

    /** @brief From design::state_feedback(any scalar). */
    template<detail::state_feedback_result R>
        requires(!Store::constant)
    constexpr LQR(const R& result) { // NOLINT
        constexpr bool kind = !R::kIntegral && !R::kEstimated;
        static_assert(kind, "LQR runs design::state_feedback: this design has integral action (LQI, LQGI) or a Kalman filter (LQG, OFLQG).");
        if constexpr (kind) {
            static_cast<Core&>(*this) = Core(result);
        }
    }

    using Core::demand_lim;
    using Core::K;
    using Core::u_lim;

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
    using Core::valid;
};

/** @brief LQR law{res}: the sizes and scalar of design::state_feedback's result. */
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, typename Kalman, size_t NZ, size_t NRU>
LQR(const design::StateFeedbackResult<NX, NU, T, NI, NY, Kalman, NZ, NRU>&) -> LQR<NX, NU, T, NY, NZ, detail::HeldGains, NRU>;

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
template<size_t NX, size_t NU, typename U, size_t NY, size_t NZ, size_t NRU, typename T>
struct lq_runtime<design::StateFeedbackResult<NX, NU, U, 0, NY, NoEstimator, NZ, NRU>, T> {
    using type = LQR<NX, NU, T, NY, NZ, HeldGains, NRU>; ///< LQR
    template<typename, typename Store>
    using on_constants = LQR<NX, NU, T, NY, NZ, Store, NRU>; ///< LQR on the design's constants
};

} // namespace detail

} // namespace damp
