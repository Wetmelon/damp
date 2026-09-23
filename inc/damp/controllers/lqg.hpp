// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lqg.hpp
 * @brief Linear-Quadratic-Gaussian design and runtime controller
 *
 * @code
 * using namespace damp;
 * constexpr auto res = design::discrete_lqg(
 *     sys, Q_lqr, R_lqr, Q_kf, R_kf);
 * static_assert(res.success);
 * LQG ctrl(res);
 * auto u = ctrl.step(y); // update + control + predict
 * @endcode
 */

#include <cstddef>

#include "damp/estimation/kalman.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#include "lqr.hpp"

namespace damp {

namespace design {

/**
 * @struct LQGResult
 * @brief LQG design result
 *
 * Holds the LQR gain design and steady-state Kalman design that, by the
 * separation principle, form the discrete LQG regulator
 * u = −K x̂ with x̂ from the Kalman filter. Use .as<float>() for deploy.
 *
 * @see discrete_lqg()
 * @see "Optimal Control" (Anderson & Moore, 1990), §8
 * @see "Optimal State Estimation" (Simon, 2006), §5–7
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = NX, size_t NV = NY>
struct LQGResult {
    LQRResult<NX, NU, T>                lqr{};          ///< LQR design result
    KalmanResult<NX, NU, NY, T, NW, NV> kalman{};       ///< Kalman filter result
    bool                                success{false}; ///< true if both designs succeeded

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return LQGResult<NX, NU, NY, U, NW, NV>{
            lqr.template as<U>(),
            kalman.template as<U>(),
            success
        };
    }

    /**
     * @brief Convert the LQG regulator to a discrete state-space block
     *
     * Current-estimator realization of @ref LQG::step: measurement update, then
     * @f$ u = -K \hat x(k|k) @f$, then predict. The compensator state is the
     * previous filtered estimate @f$ \hat x(k-1|k-1) @f$ (zero at start, matching
     * a zero predicted state and @f$ u_{\mathrm{prev}} = 0 @f$). Direct term
     * @f$ D_c = -KL @f$ is the current-estimator feedthrough.
     * @f[
     *   A_c = (I-LC)(A-BK) + LDK,\quad B_c = L,\quad C_c = -K A_c,\quad D_c = -KL.
     * @f]
     *
     * @return StateSpace with NX states, NY inputs (y), NU outputs (u)
     */
    [[nodiscard]] constexpr StateSpace<NX, NY, NU, T> to_ss() const {
        const auto& A = kalman.sys.A;
        const auto& B = kalman.sys.B;
        const auto& C = kalman.sys.C;
        const auto& D = kalman.sys.D;
        const auto& L = kalman.L;
        const auto& K = lqr.K;

        const auto I = Matrix<NX, NX, T>::identity();
        const auto ImLC = I - (L * C);
        const auto AmBK = A - (B * K);
        const auto Ac = (ImLC * AmBK) + ((L * D) * K);
        return StateSpace<NX, NY, NU, T>{
            .A = Ac,
            .B = L,
            .C = -(K * Ac),
            .D = -(K * L),
            .Ts = kalman.sys.Ts,
        };
    }
};

/**
 * @brief Discrete Linear-Quadratic-Gaussian regulator design
 *
 * Designs the optimal full-state feedback gain K (discrete LQR) and the
 * steady-state Kalman gain L for the discrete plant
 *
 *     x[k+1] = A x[k] + B u[k] + G w[k]
 *     y[k]   = C x[k] + D u[k] + H v[k]
 *
 * with w ~ N(0, Q_kf), v ~ N(0, R_kf). Separation applies: u = −K x̂ with x̂ from
 * the filter using L. Control cost is the discrete LQR criterion with weights
 * Q_lqr, R_lqr, and optional cross-term N.
 *
 * @note Compare with MATLAB®'s lqg(sys, ...) — exposed as the lqg() alias in matlab.hpp.
 *
 * @see discrete_lqr() for the control half
 * @see kalman() for the estimator half
 * @see "Optimal Control" (Anderson & Moore, 1990), §8
 * @see "Optimal State Estimation" (Simon, 2006), §5.3–5.4
 *
 * @param sys     Discrete-time state-space plant (Ts > 0)
 * @param Q_lqr   State cost (NX × NX, positive semidefinite)
 * @param R_lqr   Input cost (NU × NU, positive definite)
 * @param Q_kf    Process noise covariance (NW × NW, positive semidefinite)
 * @param R_kf    Measurement noise covariance (NV × NV, positive definite)
 * @param N       Cross-term cost (NX × NU, default zero)
 * @return LQGResult with lqr, kalman, and success = both designs succeeded
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr LQGResult<NX, NU, NY, T, NW, NV> discrete_lqg(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NX, NX, T>&                 Q_lqr,
    const Matrix<NU, NU, T>&                 R_lqr,
    const Matrix<NW, NW, T>&                 Q_kf,
    const Matrix<NV, NV, T>&                 R_kf,
    const Matrix<NX, NU, T>&                 N = Matrix<NX, NU, T>{}
) {
    const auto lqr_result = discrete_lqr(sys.A, sys.B, Q_lqr, R_lqr, N);
    // design::kalman requires Ts > 0; hand-discretized plants often leave Ts = 0.
    auto sys_d = sys;
    if (!(sys_d.Ts > T{0})) {
        sys_d.Ts = T{1};
    }
    auto kalman_result = kalman(sys_d, Q_kf, R_kf);
    kalman_result.sys.Ts = sys.Ts; // keep caller's sample-time tag on the result model
    return LQGResult<NX, NU, NY, T, NW, NV>{lqr_result, kalman_result, lqr_result.success && kalman_result.success};
}

/**
 * @brief Assemble an LQG design from separately computed Kalman and LQR results
 *
 * @note Compare with MATLAB®'s lqgreg(kest, k) — exposed as the lqgreg() alias in matlab.hpp.
 *
 * @param kest        Steady-state Kalman design (success required for overall success)
 * @param lqr_result  Discrete LQR design (success required for overall success)
 * @return LQGResult with success = kest.success && lqr_result.success
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = NX, size_t NV = NY>
[[nodiscard]] constexpr LQGResult<NX, NU, NY, T, NW, NV> lqg_from_parts(
    const KalmanResult<NX, NU, NY, T, NW, NV>& kest,
    const LQRResult<NX, NU, T>&                lqr_result
) {
    return LQGResult<NX, NU, NY, T, NW, NV>{lqr_result, kest, lqr_result.success && kest.success};
}
} // namespace design
/**
 * @ingroup discrete_controllers
 * @brief Linear-Quadratic-Gaussian (LQG) controller
 *
 * Combines LQR optimal control with a steady-state Kalman estimator (fixed L
 * from design). Separation: u = −K x̂ with x̂ from L (classical LQG). The full
 * time-varying KalmanFilter remains available for non-LQG use.
 *
 * Prefer step for a self-contained tick from the measurement alone (mirrors
 * @ref LQGI::step). The split @ref predict / @ref update / @ref control path remains
 * for multirate use and for feeding a post-saturation applied input into the filter.
 *
 * @tparam NX Number of states
 * @tparam NU Number of control inputs
 * @tparam NY Number of outputs
 * @tparam NW Number of process noise inputs (default: NX)
 * @tparam NV Number of measurement noise inputs (default: NY)
 * @tparam T  Scalar type (default: float — embedded runtime)
 */
template<size_t NX, size_t NU, size_t NY, typename T = float, size_t NW = NX, size_t NV = NY>
struct LQG {
    StateFeedback<NX, NU, T>                       lqr{};
    SteadyStateKalmanFilter<NX, NU, NY, T, NW, NV> kf{};     ///< Fixed-L estimator (designed L)
    ColVec<NU, T>                                  u_prev{}; ///< Last applied input, for the estimator's predict step

    constexpr LQG() = default;

    constexpr LQG(
        const StateFeedback<NX, NU, T>&                       lqr_,
        const SteadyStateKalmanFilter<NX, NU, NY, T, NW, NV>& kf_
    )
        : lqr(lqr_), kf(kf_) {}

    /// From design result (any scalar); converts via nested runtime ctors.
    template<typename U>
    constexpr LQG(const design::LQGResult<NX, NU, NY, U, NW, NV>& result) // NOLINT
        : lqr(result.lqr), kf(result.kalman) {}

    template<typename U>
    constexpr LQG(const LQG<NX, NU, NY, U, NW, NV>& other)
        : lqr(other.lqr), kf(other.kf), u_prev(other.u_prev.template as<T>()) {} // NOLINT

    constexpr void predict(const ColVec<NU, T>& u = ColVec<NU, T>{}) { kf.predict(u); }
    constexpr bool update(const ColVec<NY, T>& y, const ColVec<NU, T>& u = ColVec<NU, T>{}) {
        return kf.update(y, u);
    }

    /// Regulator control at the current estimate: u = -K*x̂.
    [[nodiscard]] constexpr ColVec<NU, T> control() const { return lqr.control(kf.state()); }

    /// Servo control at the current estimate: u = -K*(x̂ - x_ref). The feedback half of a
    /// 2-DOF law with no integral action, so it droops without the u_ff feedforward — see
    /// StateFeedback::control(const ColVec<NX,T>&, const ColVec<NX,T>&).
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NX, T>& x_ref) const {
        return lqr.control(x_ref, kf.state());
    }

    /**
     * @brief Measurement update + control law: the raw (pre-saturation) command for this tick.
     *
     * Runs the steady-state Kalman measurement update with the previously applied input,
     * then the regulator at the estimate. Complete the tick with commit of the input
     * the plant actually received (post-saturation). For plants without actuator limits,
     * prefer step.
     *
     * @param y Measured output
     * @return Raw control input u (pre-saturation)
     */
    [[nodiscard]] constexpr ColVec<NU, T> feedback(const ColVec<NY, T>& y) {
        kf.update(y, u_prev);
        return control();
    }

    /// @copydoc feedback(const ColVec<NY,T>&)
    /// @param x_ref State reference the feedback regulates toward.
    [[nodiscard]] constexpr ColVec<NU, T> feedback(const ColVec<NY, T>& y, const ColVec<NX, T>& x_ref) {
        kf.update(y, u_prev);
        return control(x_ref);
    }

    /**
     * @brief Complete the tick: propagate the estimator with the input the plant actually received.
     * @param u_applied The realized plant input (post-saturation, including feedforward).
     */
    constexpr void commit(const ColVec<NU, T>& u_applied) {
        kf.predict(u_applied);
        u_prev = u_applied;
    }

    /**
     * @brief One self-contained regulator tick from the measurement alone.
     *
     * @ref feedback then commit of the raw command. Use the split pair when the
     * command can saturate so the estimator sees the applied input.
     *
     * @param y Measured output
     * @return Control input u
     */
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NY, T>& y) {
        const auto u = feedback(y);
        commit(u);
        return u;
    }

    /// @copydoc step(const ColVec<NY,T>&)
    /// @param x_ref State reference the feedback regulates toward.
    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NY, T>& y, const ColVec<NX, T>& x_ref) {
        const auto u = feedback(y, x_ref);
        commit(u);
        return u;
    }

    /// Clear the applied-input memory (the estimator state is left untouched).
    constexpr void reset() { u_prev = ColVec<NU, T>{}; }
};

} // namespace damp