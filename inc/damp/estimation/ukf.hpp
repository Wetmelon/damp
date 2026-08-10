// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file ukf.hpp
 * @brief Unscented (sigma-point) Kalman Filter for nonlinear systems
 *
 * The Unscented Kalman Filter (UKF) propagates a deterministic set of
 * sigma points through the exact nonlinear dynamics/measurement
 * functions and reconstructs the posterior mean and covariance from the
 * transformed points. Unlike the EKF (@ref ekf.hpp), it never forms a
 * Jacobian: the user supplies only the plain nonlinear maps
 * `f(x, u)` and `h(x, u)`. The unscented transform is accurate to at least
 * second order for any nonlinearity (third order for Gaussian priors),
 * which is why the UKF typically outperforms the EKF on strongly nonlinear
 * problems where first-order linearization is poor.
 *
 * The sigma points are the columns of \f$\pm\sqrt{(n+\lambda)P}\f$ offset
 * about the current mean, with the matrix square root taken via the lower
 * Cholesky factor. The scaling parameters follow the standard
 * \f$(\alpha,\beta,\kappa)\f$ parameterization of the scaled unscented
 * transform.
 *
 * MATLAB® equivalent: `unscentedKalmanFilter` (Control System / Sensor Fusion
 * toolboxes). The short alias UKF is provided.
 *
 * @see Julier & Uhlmann, "Unscented Filtering and Nonlinear Estimation,"
 *      Proc. IEEE 92(3), 2004. https://doi.org/10.1109/JPROC.2003.823141
 * @see Wan & van der Merwe, "The Unscented Kalman Filter for Nonlinear
 *      Estimation," IEEE AS-SPCC, 2000. https://doi.org/10.1109/ASSPCC.2000.882463
 * @see "Optimal State Estimation" (Simon, 2006), §14.3
 */

#include <concepts>
#include <cstddef>

#include "damp/backend.hpp"
#include "damp/matrix/matrix.hpp"

namespace damp {

/**
 * @brief Concept for UKF state (process) functions
 *
 * A valid state function maps the current state and input to the next state,
 * @f$x[k+1] = f(x[k], u[k])@f$. No Jacobian is required — the UKF propagates
 * sigma points through @p f directly.
 */
template<typename Fn, typename T, size_t NX, size_t NU>
concept UKFStateFn = requires(Fn&& fn, const ColVec<NX, T>& x, const ColVec<NU, T>& u) {
    { fn(x, u) } -> std::convertible_to<ColVec<NX, T>>;
};

/**
 * @brief Concept for UKF measurement functions
 *
 * A valid measurement function maps state and input to the predicted output,
 * @f$y = h(x, u)@f$. No Jacobian is required.
 */
template<typename Fn, typename T, size_t NX, size_t NU, size_t NY>
concept UKFMeasFn = requires(Fn&& fn, const ColVec<NX, T>& x, const ColVec<NU, T>& u) {
    { fn(x, u) } -> std::convertible_to<ColVec<NY, T>>;
};

/**
 * @brief Tuning parameters for the scaled unscented transform
 *
 * @f$\lambda = \alpha^2 (n + \kappa) - n@f$ sets the spread of the sigma
 * points about the mean.
 *
 * - @p alpha — spread of the sigma points (typically small, e.g. 1e-3).
 * - @p beta  — prior-distribution knowledge; @p beta = 2 is optimal for
 *   Gaussians.
 * - @p kappa — secondary scaling, usually 0 or @f$3 - n@f$.
 */
template<typename T = double>
struct UnscentedParams {
    T alpha{static_cast<T>(1e-3)};
    T beta{static_cast<T>(2)};
    T kappa{static_cast<T>(0)};
};

/**
 * @brief Unscented (sigma-point) Kalman Filter for nonlinear discrete-time systems
 *
 * Estimates the state of a system of the form
 *
 *     x[k+1] = f(x[k], u[k]) + w[k],   w ~ N(0, Q)
 *     y[k]   = h(x[k], u[k]) + v[k],   v ~ N(0, R)
 *
 * The user provides two plain callables (no Jacobians):
 *   - `state_fn(x, u) → ColVec\<NX\>`   the nonlinear dynamics f
 *   - `meas_fn(x, u)  → ColVec\<NY\>`   the nonlinear measurement h
 *
 * mirroring the @ref ExtendedKalmanFilter API. Each step draws
 * @f$2n+1@f$ sigma points, transforms them through the exact nonlinearity,
 * and recovers the mean/covariance by the weighted unscented transform.
 *
 * @par Example
 * @code
 * // Pendulum: state [theta, omega], measure sin(theta) (nonlinear).
 * auto f = [&](const ColVec<2>& x, const ColVec<1>&) {
 *     return ColVec<2>{{x[0] + dt * x[1]},
 *                      {x[1] - dt * (g / L) * damp::sin(x[0])}};
 * };
 * auto h = [&](const ColVec<2>& x, const ColVec<1>&) {
 *     return ColVec<1>{{damp::sin(x[0])}};
 * };
 * UnscentedKalmanFilter<2, 1, 1> ukf{x0, P0, Q};
 * ukf.predict(f);
 * ukf.update(h, y_meas, R);
 * @endcode
 *
 * @see Wan & van der Merwe (2000); Julier & Uhlmann (2004); Simon (2006) §14.3.
 *
 * @tparam NX Number of states
 * @tparam NU Number of inputs
 * @tparam NY Number of outputs
 * @tparam T  Scalar type (default: float — embedded runtime)
 */
template<size_t NX, size_t NU, size_t NY, typename T = float>
struct UnscentedKalmanFilter {
    static constexpr size_t NSIG = (2 * NX) + 1; ///< Number of sigma points

    constexpr UnscentedKalmanFilter() { compute_weights(); }

    constexpr UnscentedKalmanFilter(
        const ColVec<NX, T>&      x0,
        const Matrix<NX, NX, T>&  P0,
        const Matrix<NX, NX, T>&  Q_,
        const UnscentedParams<T>& params = {}
    ) : x(x0), P(P0), Q(Q_), p(params) {
        compute_weights();
    }

    // Type conversion constructor
    template<typename U>
    constexpr UnscentedKalmanFilter(const UnscentedKalmanFilter<NX, NU, NY, U>& other)
        : x(other.state()),
          P(other.covariance()),
          Q(other.process_noise_covariance()),
          innov(other.innovation()),
          p{
              static_cast<T>(other.params().alpha),
              static_cast<T>(other.params().beta),
              static_cast<T>(other.params().kappa),
          } {
        compute_weights();
    }

    /**
     * @brief Predict step: propagate the state and covariance through the dynamics
     *
     *     {X_i}     = sigma points of (x, P)
     *     X_i[k+1]  = f(X_i, u)
     *     x[k+1|k]  = Σ Wm_i X_i[k+1]
     *     P[k+1|k]  = Σ Wc_i (X_i − x)(X_i − x)ᵀ + Q
     *
     * @param state_fn Callable (x, u) → ColVec\<NX\>
     * @param u        Control input vector
     * @return true on success, false if the prior covariance is not
     *         positive-definite (sigma points cannot be drawn)
     */
    template<typename StateFn>
        requires UKFStateFn<StateFn, T, NX, NU>
    constexpr bool predict(const StateFn& state_fn, const ColVec<NU, T>& u = ColVec<NU, T>{}) {
        const auto sigma_opt = sigma_points();
        if (!sigma_opt) {
            return false;
        }
        const auto& sig = sigma_opt.value();

        damp::array<ColVec<NX, T>, NSIG> Xp{};
        for (size_t i = 0; i < NSIG; ++i) {
            Xp[i] = state_fn(sig[i], u);
        }

        x = weighted_mean<NX>(Xp);
        Matrix<NX, NX, T> P_new = Q;
        for (size_t i = 0; i < NSIG; ++i) {
            const ColVec<NX, T> dx = Xp[i] - x;
            P_new = P_new + (dx * dx.t()) * Wc[i];
        }
        // Re-symmetrize: outer-product accumulation can break P = Pᵀ under float
        P = static_cast<T>(0.5) * (P_new + P_new.transpose());
        return true;
    }

    /**
     * @brief Update step: correct the state estimate from a measurement
     *
     *     {X_i}  = sigma points of (x, P)
     *     Y_i    = h(X_i, u)
     *     ŷ      = Σ Wm_i Y_i
     *     Pyy    = Σ Wc_i (Y_i − ŷ)(Y_i − ŷ)ᵀ + R
     *     Pxy    = Σ Wc_i (X_i − x)(Y_i − ŷ)ᵀ
     *     K      = Pxy Pyy⁻¹
     *     x[k|k] = x + K (y − ŷ)
     *     P[k|k] = P − K Pyy Kᵀ
     *
     * @param meas_fn Callable (x, u) → ColVec\<NY\>
     * @param y       Actual measurement vector
     * @param R       Measurement noise covariance
     * @param u       Control input vector
     * @return true on success, false if the covariance is not PD or the
     *         innovation covariance Pyy is singular
     */
    template<typename MeasFn>
        requires UKFMeasFn<MeasFn, T, NX, NU, NY>
    constexpr bool update(const MeasFn& meas_fn, const ColVec<NY, T>& y, const Matrix<NY, NY, T>& R, const ColVec<NU, T>& u = ColVec<NU, T>{}) {
        const auto sigma_opt = sigma_points();
        if (!sigma_opt) {
            return false;
        }
        const auto& sig = sigma_opt.value();

        damp::array<ColVec<NY, T>, NSIG> Yp{};
        for (size_t i = 0; i < NSIG; ++i) {
            Yp[i] = meas_fn(sig[i], u);
        }

        const ColVec<NY, T> y_pred = weighted_mean<NY>(Yp);

        Matrix<NY, NY, T> Pyy = R;
        Matrix<NX, NY, T> Pxy = Matrix<NX, NY, T>::zeros();
        for (size_t i = 0; i < NSIG; ++i) {
            const ColVec<NY, T> dy = Yp[i] - y_pred;
            const ColVec<NX, T> dx = sig[i] - x;
            Pyy = Pyy + (dy * dy.t()) * Wc[i];
            Pxy = Pxy + (dx * dy.t()) * Wc[i];
        }

        // Re-symmetrize Pyy before the SPD solve (sigma-point sum can skew off-diagonals).
        Pyy = static_cast<T>(0.5) * (Pyy + Pyy.transpose());

        // K = Pxy Pyy⁻¹ → solve Pyy Kᵀ = Pxyᵀ (Cholesky; LU if Pyy only semi-definite).
        // Fail closed if Pyy is singular (return false; state/P unchanged by early return).
        auto Kt_opt = mat::cholesky_solve(Pyy, Pxy.transpose());
        if (!Kt_opt) {
            Kt_opt = mat::lu_solve(Pyy, Pxy.transpose());
        }
        if (!Kt_opt) {
            return false;
        }
        const Matrix<NX, NY, T> K = Kt_opt.value().transpose();

        innov = y - y_pred;
        x = x + K * innov;
        // Joseph form: P − K Pxyᵀ − Pxy Kᵀ + K Pyy Kᵀ
        // (equiv. to P − K Pyy Kᵀ when K = Pxy Pyy⁻¹; better PD under float)
        const Matrix<NX, NX, T> P_post = P - (K * Pxy.t()) - (Pxy * K.t()) + quadratic_form(K, Pyy);
        P = static_cast<T>(0.5) * (P_post + P_post.transpose());
        return true;
    }

    // Accessors
    [[nodiscard]] constexpr const auto& innovation() const { return innov; }
    [[nodiscard]] constexpr const auto& state() const { return x; }
    [[nodiscard]] constexpr const auto& covariance() const { return P; }
    [[nodiscard]] constexpr const auto& process_noise_covariance() const { return Q; }
    [[nodiscard]] constexpr const auto& params() const { return p; }

    // Mutators. As with the EKF these support inter-step intervention —
    // most importantly sequential scalar updates with inter-measurement state
    // clamping (run a scalar update with NY == 1, clamp the affected state,
    // write it back, then run the next scalar update against the clamped estimate).
    constexpr void set_state(const ColVec<NX, T>& x_new) { x = x_new; }
    constexpr void set_state(size_t i, T value) { x[i] = value; }
    constexpr void set_covariance(const Matrix<NX, NX, T>& P_new) { P = P_new; }
    constexpr void set_process_noise_covariance(const Matrix<NX, NX, T>& Q_new) { Q = Q_new; }

    // Re-initialize the estimate and clear the innovation (Q and the unscented
    // weights are kept).
    constexpr void reset(
        const ColVec<NX, T>&     x0 = ColVec<NX, T>{},
        const Matrix<NX, NX, T>& P0 = Matrix<NX, NX, T>::identity()
    ) {
        x = x0;
        P = P0;
        innov = ColVec<NY, T>{};
    }

private:
    /// Compute the scaled-unscented-transform mean/covariance weights.
    constexpr void compute_weights() {
        const T n = static_cast<T>(NX);
        lambda = (p.alpha * p.alpha * (n + p.kappa)) - n;
        const T denom = n + lambda;
        Wm[0] = lambda / denom;
        Wc[0] = (lambda / denom) + (static_cast<T>(1) - (p.alpha * p.alpha) + p.beta);
        const T w = static_cast<T>(1) / (static_cast<T>(2) * denom);
        for (size_t i = 1; i < NSIG; ++i) {
            Wm[i] = w;
            Wc[i] = w;
        }
    }

    /// Draw the 2n+1 sigma points from (x, P); nullopt if P is not PD.
    [[nodiscard]] constexpr damp::optional<damp::array<ColVec<NX, T>, NSIG>> sigma_points() const {
        const auto L_opt = mat::cholesky(P * (static_cast<T>(NX) + lambda));
        if (!L_opt) {
            return damp::nullopt;
        }
        const Matrix<NX, NX, T> L = L_opt.value();

        damp::array<ColVec<NX, T>, NSIG> sig{};
        sig[0] = x;
        for (size_t j = 0; j < NX; ++j) {
            const ColVec<NX, T> dcol = L.col(j).to_vector();
            sig[j + 1] = x + dcol;
            sig[j + 1 + NX] = x - dcol;
        }
        return sig;
    }

    /// Weighted (mean-weight) combination of a sigma-point set.
    template<size_t N>
    [[nodiscard]] constexpr ColVec<N, T> weighted_mean(const damp::array<ColVec<N, T>, NSIG>& pts) const {
        ColVec<N, T> m{};
        for (size_t i = 0; i < NSIG; ++i) {
            m = m + pts[i] * Wm[i];
        }
        return m;
    }

    ColVec<NX, T>      x{};
    Matrix<NX, NX, T>  P{};
    Matrix<NX, NX, T>  Q{};
    ColVec<NY, T>      innov{};
    UnscentedParams<T> p{};

    T                    lambda{};
    damp::array<T, NSIG> Wm{};
    damp::array<T, NSIG> Wc{};
};

/**
 * @brief Short alias for @ref UnscentedKalmanFilter (MATLAB®-style name)
 *
 * @see UnscentedKalmanFilter
 */
template<size_t NX, size_t NU, size_t NY, typename T = double>
using UKF = UnscentedKalmanFilter<NX, NU, NY, T>;

} // namespace damp
