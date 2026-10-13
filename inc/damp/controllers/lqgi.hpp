// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lqgi.hpp
 * @brief LQG with integral action design and runtime controller
 */

#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#include "lqi.hpp"

namespace damp {
namespace design {
/**
 * @struct LQGIResult
 * @brief LQGI design result
 *
 * LQI gain design plus steady-state Kalman design for output tracking with
 * integral action on estimated state. Use .as<float>() for deploy.
 *
 * @see discrete_lqgi()
 * @see "Optimal Control" (Anderson & Moore, 1990), §8
 * @see "Optimal State Estimation" (Simon, 2006), §5–7
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = NX, size_t NV = NY, size_t NI = NY>
struct LQGIResult {
    LQIResult<NX, NU, NY, T, NI>        lqi{};
    KalmanResult<NX, NU, NY, T, NW, NV> kalman{};
    bool                                success{false};

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return LQGIResult<NX, NU, NY, U, NW, NV, NI>{
            lqi.template as<U>(),
            kalman.template as<U>(),
            success
        };
    }

    /**
     * @brief Convert the LQGI compensator to a discrete state-space block
     *
     * Current-estimator realization of @ref LQGI::step. State is the runtime
     * registers at the start of the tick, @f$ [\hat x(k|k-1);\; \xi_i;\; u_{\mathrm{prev}}] @f$:
     * measurement update with the previous applied input, then
     * @f$ u = -[K_x\; K_i][\hat x(k|k);\; \xi_i] @f$, then @f$ \xi_i \leftarrow \xi_i + (r - C_i y) @f$
     * and predict. Inputs are @f$ [r; y] @f$ (r for the NI integrated outputs).
     *
     * @return StateSpace with NX+NI+NU states, NI+NY inputs ([r; y]), NU outputs (u)
     */
    [[nodiscard]] constexpr StateSpace<NX + NI + NU, NI + NY, NU, T> to_ss() const {
        constexpr size_t NS = NX + NI + NU;
        const auto&      A = kalman.sys.A;
        const auto&      B = kalman.sys.B;
        const auto&      C = kalman.sys.C;
        const auto&      D = kalman.sys.D;
        const auto&      L = kalman.L;

        const Matrix<NU, NX, T> Kx = lqi.K.template block<NU, NX>(0, 0);
        const Matrix<NU, NI, T> Ki = lqi.K.template block<NU, NI>(0, NX);

        const auto I = Matrix<NX, NX, T>::identity();
        const auto ImLC = I - (L * C);
        const auto LD = L * D;
        const auto KxImLC = Kx * ImLC;
        const auto AmBKx = A - (B * Kx);
        const auto KxL = Kx * L;

        Matrix<NS, NS, T> Ac{};
        Ac.template block<NX, NX>(0, 0) = AmBKx * ImLC;
        Ac.template block<NX, NI>(0, NX) = -(B * Ki);
        Ac.template block<NX, NU>(0, NX + NI) = -(AmBKx * LD);
        Ac.template block<NI, NI>(NX, NX) = Matrix<NI, NI, T>::identity();
        Ac.template block<NU, NX>(NX + NI, 0) = -KxImLC;
        Ac.template block<NU, NI>(NX + NI, NX) = -Ki;
        Ac.template block<NU, NU>(NX + NI, NX + NI) = Kx * LD;

        Matrix<NS, NI + NY, T> Bc{};
        Bc.template block<NI, NI>(NX, 0) = Matrix<NI, NI, T>::identity();
        Bc.template block<NX, NY>(0, NI) = AmBKx * L;
        Bc.template block<NI, NY>(NX, NI) = -lqi.Ci;
        Bc.template block<NU, NY>(NX + NI, NI) = -KxL;

        Matrix<NU, NS, T> Cc{};
        Cc.template block<NU, NX>(0, 0) = -KxImLC;
        Cc.template block<NU, NI>(0, NX) = -Ki;
        Cc.template block<NU, NU>(0, NX + NI) = Kx * LD;

        Matrix<NU, NI + NY, T> Dc{};
        Dc.template block<NU, NY>(0, NI) = -KxL;

        return StateSpace<NS, NI + NY, NU, T>{
            .A = Ac,
            .B = Bc,
            .C = Cc,
            .D = Dc,
            .Ts = kalman.sys.Ts,
        };
    }
};

/**
 * @brief Discrete LQG with integral action (LQI + Kalman) for output tracking
 *
 * Combines @ref discrete_lqi on the plant with @ref kalman for state estimation.
 * Runtime law: u = −K [x̂; xi] with xi integrating (r − Ci·y). Separation applies
 * between the LQI gain on the true augmented plant and the Kalman estimator.
 * Ci selects the integrated outputs exactly as in @ref discrete_lqi: pass a
 * selector matrix or a list of output indices as the last argument, or leave it
 * out to integrate every output. On a model with delay states (discretize with
 * InputDelay / OutputDelay) Q_aug may be sized to the plant,
 * (NX_plant+NI)×(NX_plant+NI), as in @ref discrete_lqi; the process noise is already
 * plant-sized through G.
 *
 * @note Compare with MATLAB®'s lqgtrack(...) — exposed as the lqgtrack() alias in matlab.hpp.
 *
 * @see discrete_lqi()
 * @see kalman()
 * @see discrete_lqg() for regulation without integral action
 * @see "Optimal Control" (Anderson & Moore, 1990), §8
 * @see "Optimal State Estimation" (Simon, 2006), §5.3–5.4
 *
 * @param sys      Discrete-time state-space plant
 * @param Q_aug    Augmented state cost ((NX+NI)×(NX+NI), or plant-sized on a delayed model; positive semidefinite)
 * @param R        Input cost (NU×NU, positive definite)
 * @param noise    Process (NW×NW, positive semidefinite) and measurement (NV×NV, positive
 *                 definite) noise covariances
 * @param Ci       Integrated-output selector (NI×NY)
 * @return LQGIResult with success = both LQI and Kalman designs succeeded
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NI, size_t NQA = NX + NI>
[[nodiscard]] constexpr LQGIResult<NX, NU, NY, T, NW, NV, NI> discrete_lqgi(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NQA, NQA, T>&               Q_aug,
    const Matrix<NU, NU, T>&                 R,
    const KalmanNoise<NW, NV, T>&            noise,
    const Matrix<NI, NY, T>&                 Ci
) {
    const auto lqi_result = discrete_lqi(sys, Q_aug, R, Ci);
    // design::kalman requires Ts > 0; hand-discretized plants often leave Ts = 0.
    auto sys_d = sys;
    if (!(sys_d.Ts > T{0})) {
        sys_d.Ts = T{1};
    }
    auto kalman_result = kalman(sys_d, noise.process, noise.measurement);
    kalman_result.sys.Ts = sys.Ts; // keep caller's sample-time tag on the result model
    return LQGIResult<NX, NU, NY, T, NW, NV, NI>{lqi_result, kalman_result, lqi_result.success && kalman_result.success};
}

/// Every output integrated (Ci = I); see the selector overload.
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0, size_t NQA = NX + NY>
[[nodiscard]] constexpr LQGIResult<NX, NU, NY, T, NW, NV> discrete_lqgi(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NQA, NQA, T>&               Q_aug,
    const Matrix<NU, NU, T>&                 R,
    const KalmanNoise<NW, NV, T>&            noise = {}
) {
    return discrete_lqgi(sys, Q_aug, R, noise, Matrix<NY, NY, T>::identity());
}

/// Integrate the listed outputs only, e.g. `discrete_lqgi(sys, Q_aug, R, {.process = Qn, .measurement = Rn}, {0, 3})`.
/// An index outside [0, NY) reports success = false.
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NI, size_t NQA = NX + NI>
[[nodiscard]] constexpr LQGIResult<NX, NU, NY, T, NW, NV, NI> discrete_lqgi(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NQA, NQA, T>&               Q_aug,
    const Matrix<NU, NU, T>&                 R,
    const KalmanNoise<NW, NV, T>&            noise,
    const size_t (&outputs)[NI]
) {
    const auto Ci = damp::detail::output_selector<NI, NY, T>(outputs);
    if (!Ci) {
        return LQGIResult<NX, NU, NY, T, NW, NV, NI>{};
    }
    return discrete_lqgi(sys, Q_aug, R, noise, Ci.value());
}
} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Linear-Quadratic-Gaussian-Integral (LQGI) controller
 *
 * Combines LQI output tracking with a steady-state Kalman estimator (fixed L
 * from design). Classical LQGI separation: u = −K [x̂; xi] with x̂ from L.
 *
 * @tparam NX Number of states
 * @tparam NU Number of control inputs
 * @tparam NY Number of outputs
 * @tparam NW Number of process noise inputs (default: NX)
 * @tparam NV Number of measurement noise inputs (default: NY)
 * @tparam T  Scalar type (default: float — embedded runtime)
 * @tparam NI Number of integrated outputs (default: NY)
 */
template<size_t NX, size_t NU, size_t NY, typename T = float, size_t NW = NX, size_t NV = NY, size_t NI = NY>
struct LQGI {
    LQI<NX, NU, NY, T, NI>                         lqi{};
    SteadyStateKalmanFilter<NX, NU, NY, T, NW, NV> kf{};         ///< Fixed-L estimator (designed L)
    ColVec<NU, T>                                  u_prev{};     ///< Last applied input, for the estimator's predict step
    ColVec<NX, T>                                  x_ref_prev{}; ///< Last full state reference, to extend a plant-sized one over delay states

    constexpr LQGI() = default;

    constexpr LQGI(
        const LQI<NX, NU, NY, T, NI>&                         lqi_,
        const SteadyStateKalmanFilter<NX, NU, NY, T, NW, NV>& kf_
    )
        : lqi(lqi_), kf(kf_) {}

    /// From design result (any scalar); converts via nested runtime ctors.
    template<typename U>
    constexpr LQGI(const design::LQGIResult<NX, NU, NY, U, NW, NV, NI>& result) // NOLINT
        : lqi(result.lqi), kf(result.kalman) {
        if (!result.success && std::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    template<typename U>
    constexpr LQGI(const LQGI<NX, NU, NY, U, NW, NV, NI>& other)
        : lqi(other.lqi), kf(other.kf), u_prev(other.u_prev.template as<T>()), x_ref_prev(other.x_ref_prev.template as<T>()) {} // NOLINT

    constexpr void predict(const ColVec<NU, T>& u = ColVec<NU, T>{}) { kf.predict(u); }
    constexpr bool update(const ColVec<NY, T>& y, const ColVec<NU, T>& u = ColVec<NU, T>{}) {
        return kf.update(y, u);
    }

    /**
     * @brief Compute control with integral action from the augmented state.
     *
     * @param x_aug Augmented state [x̂; xi] — estimated plant state stacked on
     *              the integral-of-tracking-error state (size NX + NI).
     * @return Control input u = −K·x_aug
     */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NX + NI, T>& x_aug) const {
        return lqi.control(x_aug);
    }

    /**
     * @brief Compute control using the estimate and the controller's own integrator
     *
     * Pulls the plant estimate from the Kalman filter and the integral state from
     * the embedded LQI, so the caller only supplies the reference and measurement:
     * u = -[Kx Ki]·[x̂; xi], then xi advances by (r − Ci·y).
     *
     * @warning This does NOT advance the estimator — call predict(u) and
     *          update(y) first each tick (caller-sequenced). It is not an
     *          OutputFeedbackController tick; @ref step is the self-contained
     *          counterpart, and @ref feedback / commit the saturation-aware one.
     *
     * @param r Reference for the integrated outputs
     * @param y Measured output
     * @return Control input u
     */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NI, T>& r, const ColVec<NY, T>& y) {
        return lqi.control(r, y, kf.state());
    }

    /// @copydoc control(const ColVec<NI,T>&, const ColVec<NY,T>&)
    /// @param x_ref State reference the feedback regulates toward (see LQI::control).
    [[nodiscard]] constexpr ColVec<NU, T>
    control(const ColVec<NI, T>& r, const ColVec<NY, T>& y, const ColVec<NX, T>& x_ref) {
        return lqi.control(r, y, kf.state(), x_ref);
    }

    /**
     * @brief Measurement update + control law: the raw (pre-saturation) command for this tick.
     *
     * Runs the steady-state Kalman measurement update with the previously applied input,
     * then the LQI control law at the estimate. The tick is not complete until commit
     * is called with the input the plant actually received (post-saturation, including any
     * feedforward): commit propagates the estimator with it and advances the integrators
     * except those that would wind up against the clamp (see LQI::commit).
     * @code
     * const auto u_req = ctrl.feedback(r, y);
     * const auto u     = clamp_each(u_req + u_ff);
     * ctrl.commit(u, u_ff);
     * @endcode
     *
     * @param r Reference for the integrated outputs
     * @param y Measured output
     * @return Raw control input u (pre-saturation)
     */
    [[nodiscard]] constexpr ColVec<NU, T> feedback(const ColVec<NI, T>& r, const ColVec<NY, T>& y) {
        kf.update(y, u_prev);
        return lqi.feedback(r, y, kf.state());
    }

    /// @copydoc feedback(const ColVec<NI,T>&, const ColVec<NY,T>&)
    /// @param x_ref State reference the feedback regulates toward (see LQI::control).
    [[nodiscard]] constexpr ColVec<NU, T>
    feedback(const ColVec<NI, T>& r, const ColVec<NY, T>& y, const ColVec<NX, T>& x_ref) {
        kf.update(y, u_prev);
        x_ref_prev = x_ref;
        return lqi.feedback(r, y, kf.state(), x_ref);
    }

    /**
     * @brief Measurement update + LQI law on a model with delay states (discretize with
     *        InputDelay / OutputDelay): @p x_ref covers the NR plant states only.
     *
     * The delay states that follow are filled in from the reference history with no
     * feedforward input (the integrator supplies it): past inputs regulate toward zero, past
     * outputs toward those the reference produced.
     */
    template<size_t NR>
        requires(NR < NX)
    [[nodiscard]] constexpr ColVec<NU, T> feedback(const ColVec<NI, T>& r, const ColVec<NY, T>& y, const ColVec<NR, T>& x_ref) {
        return feedback(r, y, detail::extend_reference<NX>(kf.model().A, kf.model().B, x_ref_prev, ColVec<NU, T>{}, x_ref));
    }

    /**
     * @brief Complete the tick: integrate (anti-windup) and propagate the estimator with the applied input.
     * @param u_applied The realized plant input (post-saturation, including feedforward).
     * @param u_ff      Feedforward added to the command before the clamp (default: none).
     */
    constexpr void commit(const ColVec<NU, T>& u_applied, const ColVec<NU, T>& u_ff = ColVec<NU, T>{}) {
        lqi.commit(u_applied, u_ff);
        kf.predict(u_applied);
        u_prev = u_applied;
    }

    /**
     * @brief One self-contained tick: @ref feedback then commit of the raw command.
     *
     * Convenience for plants without actuator limits. If the command can saturate, use the
     * split @ref feedback / commit pair so the estimator sees the applied input.
     */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NI, T>& r, const ColVec<NY, T>& y) {
        const auto u = feedback(r, y);
        commit(u);
        return u;
    }

    /**
     * @brief One self-contained tick with an explicit state reference.
     *
     * Same as the two-argument @c step, but passes @p x_ref into LQI feedback.
     *
     * @param r     Output reference
     * @param y     Measured output
     * @param x_ref State reference the feedback regulates toward (see LQI::control)
     */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NI, T>& r, const ColVec<NY, T>& y, const ColVec<NX, T>& x_ref) {
        const auto u = feedback(r, y, x_ref);
        commit(u);
        return u;
    }

    /// One self-contained tick, @p x_ref on the plant states only (see feedback).
    template<size_t NR>
        requires(NR < NX)
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NI, T>& r, const ColVec<NY, T>& y, const ColVec<NR, T>& x_ref) {
        const auto u = feedback(r, y, x_ref);
        commit(u);
        return u;
    }

    /// Clear the integral state and the applied-input and reference memory (the estimator state is left untouched).
    constexpr void reset() {
        lqi.reset();
        u_prev = ColVec<NU, T>{};
        x_ref_prev = ColVec<NX, T>{};
    }
};

} // namespace damp