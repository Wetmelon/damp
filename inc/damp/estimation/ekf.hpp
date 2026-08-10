// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file ekf.hpp
 * @brief Extended Kalman Filter for nonlinear systems
 *
 * The user provides callable objects that compute both the nonlinear
 * prediction/measurement AND the corresponding Jacobians in one call.
 * This avoids redundant computation when the Jacobians share intermediate
 * terms with the function evaluation.
 *
 * @note The callable returns a struct containing the function value and its
 *       Jacobians — it is NOT the Jacobian function itself.
 *
 * @see "Optimal State Estimation" (Simon, 2006), Chapter 13
 * @see Solà et al., "Quaternion kinematics for the error-state Kalman filter" (2017)
 */

#include <concepts>
#include <cstddef>

#include "damp/matrix/matrix.hpp"

namespace damp {

/**
 * @brief State prediction result from the user's dynamics function
 *
 * The user's state function computes x[k+1] = f(x[k], u[k]) and returns
 * both the predicted state and the Jacobians needed for covariance propagation.
 *
 * @tparam T  Scalar type
 * @tparam NX Number of states
 */
template<typename T, size_t NX>
struct StateJacobian {
    ColVec<NX, T>     x_pred{};                         ///< Predicted state x[k+1] = f(x[k], u[k])
    Matrix<NX, NX, T> F{};                              ///< State Jacobian ∂f/∂x
    Matrix<NX, NX, T> G{Matrix<NX, NX, T>::identity()}; ///< Process noise Jacobian ∂f/∂w
};

/**
 * @brief Concept for EKF state functions
 *
 * A valid state function takes (x, u) and returns a StateJacobian containing
 * the predicted state and linearization Jacobians F, G.
 */
template<typename Fn, typename T, size_t NX, size_t NU>
concept EKFStateFn = requires(Fn&& fn, const ColVec<NX, T>& x, const ColVec<NU, T>& u) {
    { fn(x, u) } -> std::convertible_to<StateJacobian<T, NX>>;
};

/**
 * @brief Measurement prediction result from the user's observation function
 *
 * The user's measurement function computes y_pred = h(x, u) and returns
 * both the predicted measurement and the Jacobians needed for the Kalman update.
 *
 * @tparam T  Scalar type
 * @tparam NY Number of outputs
 * @tparam NX Number of states
 */
template<typename T, size_t NY, size_t NX>
struct MeasJacobian {
    ColVec<NY, T>     y_pred{};                         ///< Predicted measurement y = h(x, u)
    Matrix<NY, NX, T> H{};                              ///< Measurement Jacobian ∂h/∂x
    Matrix<NY, NY, T> M{Matrix<NY, NY, T>::identity()}; ///< Measurement noise Jacobian ∂h/∂v
};

/**
 * @brief Concept for EKF measurement functions
 *
 * A valid measurement function takes (x, u) and returns a MeasJacobian
 * containing the predicted measurement and linearization Jacobians H, M.
 */
template<typename Fn, typename T, size_t NX, size_t NU, size_t NY>
concept EKFMeasFn = requires(Fn&& fn, const ColVec<NX, T>& x, const ColVec<NU, T>& u) {
    { fn(x, u) } -> std::convertible_to<MeasJacobian<T, NY, NX>>;
};

/**
 * @brief Extended Kalman Filter for nonlinear discrete-time systems
 *
 * Implements the standard EKF predict/update cycle for systems of the form:
 *
 *     x[k+1] = f(x[k], u[k]) + G·w[k],   w ~ N(0, Q)
 *     y[k]   = h(x[k], u[k]) + M·v[k],   v ~ N(0, R)
 *
 * The user provides two callables:
 *   - state_fn(x, u) → StateJacobian{x_pred, F, G}
 *   - meas_fn(x, u)  → MeasJacobian{y_pred, H, M}
 *
 * Each callable evaluates the nonlinear function AND returns its Jacobians
 * in one call. This is more efficient than separate function/Jacobian calls
 * when they share intermediate computations.
 *
 * The covariance update uses the Joseph form for numerical stability.
 *
 * @see "Optimal State Estimation" (Simon, 2006), §13.3
 *
 * @tparam NX Number of states
 * @tparam NU Number of inputs
 * @tparam NY Number of outputs
 * @tparam T  Scalar type (default: float — embedded runtime)
 */
template<size_t NX, size_t NU, size_t NY, typename T = float>
struct ExtendedKalmanFilter {
    constexpr ExtendedKalmanFilter() = default;

    constexpr ExtendedKalmanFilter(
        const ColVec<NX, T>&     x0,
        const Matrix<NX, NX, T>& P0,
        const Matrix<NX, NX, T>& Q_
    ) : x(x0), P(P0), Q(Q_) {}

    // Type conversion constructor
    template<typename U>
    constexpr ExtendedKalmanFilter(const ExtendedKalmanFilter<NX, NU, NY, U>& other)
        : x(other.state()), P(other.covariance()), Q(other.process_noise_covariance()), innov(other.innovation()) {}

    /**
     * @brief Predict step: propagate state and covariance through dynamics
     *
     *     x[k+1|k] = f(x[k|k], u[k])
     *     P[k+1|k] = F·P[k|k]·Fᵀ + G·Q·Gᵀ
     *
     * @param state_fn Callable (x, u) → StateJacobian{x_pred, F, G}
     * @param u        Control input vector
     */
    template<typename StateFn>
        requires EKFStateFn<StateFn, T, NX, NU>
    constexpr void predict(const StateFn& state_fn, const ColVec<NU, T>& u = ColVec<NU, T>{}) {
        StateJacobian<T, NX> sj = state_fn(x, u);
        x = sj.x_pred;
        // quadratic_form forces exact P symmetry (mirrors Joseph path in KalmanFilter)
        P = quadratic_form(sj.F, P) + quadratic_form(sj.G, Q);
    }

    /**
     * @brief Update step: correct state estimate from measurement
     *
     *     S = H·P·Hᵀ + M·R·Mᵀ
     *     K = P·Hᵀ·S⁻¹  (Cholesky; LU fallback if S is only semi-definite)
     *     x[k|k] = x[k|k−1] + K·(y − ŷ)
     *     P[k|k] = (I − KH)·P·(I − KH)ᵀ + K·M·R·Mᵀ·Kᵀ  [Joseph form]
     *
     * @param meas_fn Callable (x, u) → MeasJacobian{y_pred, H, M}
     * @param y       Actual measurement vector
     * @param R       Measurement noise covariance (must match sensor model)
     * @param u       Control input vector
     * @return true if update succeeded, false if S is singular
     */
    template<typename MeasFn>
        requires EKFMeasFn<MeasFn, T, NX, NU, NY>
    constexpr bool update(const MeasFn& meas_fn, const ColVec<NY, T>& y, const Matrix<NY, NY, T>& R, const ColVec<NU, T>& u = ColVec<NU, T>{}) {
        MeasJacobian<T, NY, NX> mj = meas_fn(x, u);
        innov = y - mj.y_pred;
        // S = HPHᵀ + MRMᵀ with exact symmetry on each term
        const Matrix<NY, NY, T> S = quadratic_form(mj.H, P) + quadratic_form(mj.M, R);

        // K = PHᵀS⁻¹ → solve S Kᵀ = H P (fail closed if S is singular)
        auto K_opt = mat::cholesky_solve(S, mj.H * P);
        if (!K_opt) {
            K_opt = mat::lu_solve(S, mj.H * P);
        }
        if (!K_opt) {
            return false;
        }

        const Matrix<NX, NY, T> K = K_opt.value().transpose();
        x = x + K * innov;

        // Joseph form: P = (I − KH) P (I − KH)ᵀ + (KM) R (KM)ᵀ
        const auto I_KH = Matrix<NX, NX, T>::identity() - K * mj.H;
        const auto KM = K * mj.M;
        P = quadratic_form(I_KH, P) + quadratic_form(KM, R);
        return true;
    }

    // Accessors
    [[nodiscard]] constexpr const auto& innovation() const { return innov; }
    [[nodiscard]] constexpr const auto& state() const { return x; }
    [[nodiscard]] constexpr const auto& covariance() const { return P; }
    [[nodiscard]] constexpr const auto& process_noise_covariance() const { return Q; }

    // Mutators. These exist so the caller can intervene *between* filter steps —
    // most importantly to support sequential scalar updates with inter-measurement
    // state clamping: run a scalar update() (NY == 1), clamp/saturate the affected
    // state component, write it back, then run the next scalar update against the
    // clamped estimate. Clamping after a fused vector update can't enforce the
    // constraint each measurement sees; doing it between scalar updates can.
    constexpr void set_state(const ColVec<NX, T>& x_new) { x = x_new; }
    constexpr void set_state(size_t i, T value) { x[i] = value; }
    constexpr void set_covariance(const Matrix<NX, NX, T>& P_new) { P = P_new; }
    constexpr void set_process_noise_covariance(const Matrix<NX, NX, T>& Q_new) { Q = Q_new; }

    // Re-initialize the estimate and clear the innovation (Q is kept).
    constexpr void reset(
        const ColVec<NX, T>&     x0 = ColVec<NX, T>{},
        const Matrix<NX, NX, T>& P0 = Matrix<NX, NX, T>::identity()
    ) {
        x = x0;
        P = P0;
        innov = ColVec<NY, T>{};
    }

private:
    ColVec<NX, T>     x{};
    Matrix<NX, NX, T> P{};
    Matrix<NX, NX, T> Q{};
    ColVec<NY, T>     innov{};
};

} // namespace damp