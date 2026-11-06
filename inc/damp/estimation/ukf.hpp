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
 * Jacobian. The model is the continuous Dynamics ẋ = f(x, u) of concepts.hpp, a struct or a
 * captureless lambda (the same object design::linearize and the sim::simulate* plants take),
 * and each reading comes with its Measurement y = h(x, u) and R; each sigma point is carried over the sample time
 * by one RK4 step of f. The unscented transform is accurate to at least
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
 * The short alias UKF is provided.
 *
 * @note Compare with MATLAB®'s unscentedKalmanFilter.
 * @see "Unscented Filtering and Nonlinear Estimation" (Julier & Uhlmann, 2004)
 * @see "The Unscented Kalman Filter for Nonlinear Estimation" (Wan & van der Merwe, 2000)
 * @see "Optimal State Estimation" (Simon, 2006), §14.3
 */


#include "damp/backend.hpp"
#include "damp/concepts.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/dynamics.hpp"

namespace damp {

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
template<typename T = float>
struct UnscentedParams {
    /**
     * @brief Spread of the sigma points. 1e-3 is a double-precision recipe; in float
     * it makes n+λ ~ 1e-6 n and a center weight ~ 1/α².
     */
    T alpha{damp::is_same_v<T, float> ? T{0.5} : static_cast<T>(1e-3)};
    T beta{static_cast<T>(2)};  ///< Prior-distribution knowledge (2 for a Gaussian)
    T kappa{static_cast<T>(0)}; ///< Secondary scaling
};

/**
 * @brief Configuration of an UnscentedKalmanFilter: its sample time, the process noise, the
 *        initial estimate and the unscented transform's tuning
 *
 * An aggregate; name the fields you set. The model is the filter's first constructor argument:
 * `UnscentedKalmanFilter ukf{f, {.Ts = 0.01F, .Q = Q, .P0 = P0}}`.
 *
 * @tparam Model Dynamics ẋ = f(x, u), a struct or a lambda whose call operator gives NX and NU
 * @tparam T     Scalar type (default: the model's)
 */
template<typename Model, typename T = detail::model_shape<Model>::scalar>
    requires detail::shaped_model<Model>
struct UnscentedKalmanFilterConfig {
    static constexpr size_t NX = detail::model_shape<Model>::states; ///< States

    T                  Ts{};                              ///< Sample time [s] predict integrates over, finite and positive
    Matrix<NX, NX, T>  Q{Matrix<NX, NX, T>::identity()};  ///< Process-noise covariance per step, on the full state
    ColVec<NX, T>      x0{};                              ///< Initial state estimate
    Matrix<NX, NX, T>  P0{Matrix<NX, NX, T>::identity()}; ///< Initial estimate covariance
    UnscentedParams<T> unscented{};                       ///< Scaled unscented transform tuning (α, β, κ)
};

/**
 * @brief Unscented (sigma-point) Kalman Filter: continuous nonlinear dynamics, sampled measurements
 *
 * Estimates the state of a system of the form
 *
 *     ẋ      = f(x, u)                       (Dynamics)
 *     x[k+1] = φ(x[k], u[k]) + w[k],   w ~ N(0, Q)
 *     y[k]   = h(x[k], u[k]) + v[k],   v ~ N(0, R)   (Measurement, per reading)
 *
 * where φ is one classical RK4 step of f over the sample time Ts with u held (zero-order
 * hold), as in the @ref ExtendedKalmanFilter. The filter holds its model, Ts and Q; each tick
 * is predict(u) and update(h, y, R) per reading, each reading of any size. The input h sees is
 * the one the last predict held. Each step draws @f$2n+1@f$ sigma points, carries them through the exact
 * nonlinearity, and recovers the mean/covariance by the weighted unscented transform.
 *
 * Q is the covariance of the state increment over one step, on the full state (G = I), as in
 * the EKF; each reading's R the covariance of an additive measurement noise.
 *
 * Cost per tick, NX states: predict evaluates f 4·(2·NX + 1) times (one RK4 step per sigma
 * point), an update evaluates h 2·NX + 1 times and factors Pyy once (Cholesky). The counts are fixed, so the worst-case time is
 * bounded.
 *
 * @par Example
 * @code
 * // Pendulum: state [θ, ω], measure sin θ (nonlinear).
 * constexpr auto pendulum = [](const ColVec<2, float>& x, const ColVec<1, float>& u) {
 *     return ColVec<2, float>{x[1], (-9.81F * damp::sin(x[0])) + u[0]};
 * };
 * constexpr auto sin_theta = [](const ColVec<2, float>& x, const ColVec<1, float>&) { return ColVec<1, float>{damp::sin(x[0])}; };
 * static constinit UnscentedKalmanFilter ukf{pendulum, {.Ts = 0.01F, .Q = Q, .x0 = x0, .P0 = P0}};
 * (void)ukf.update(sin_theta, y_meas, 1e-4F);
 * (void)ukf.predict(u);
 * @endcode
 *
 * @see "The Unscented Kalman Filter for Nonlinear Estimation" (Wan & van der Merwe, 2000)
 *
 * @tparam Model Dynamics ẋ = f(x, u), a struct or a lambda whose call operator gives NX and NU
 * @tparam T     Scalar type (default: the model's)
 */
template<typename Model, typename T = detail::model_shape<Model>::scalar>
    requires detail::shaped_model<Model> && Dynamics<Model, detail::model_shape<Model>::states, detail::model_shape<Model>::inputs, T>
struct UnscentedKalmanFilter {
    static constexpr size_t NX = detail::model_shape<Model>::states; ///< States
    static constexpr size_t NU = detail::model_shape<Model>::inputs; ///< Inputs
    static constexpr size_t NSIG = (2 * NX) + 1;                     ///< Number of sigma points

    constexpr UnscentedKalmanFilter() { compute_weights(); }

    /**
     * @brief Filter from its process model and its configuration
     *
     * A sample time that is not finite and positive is refused under constant evaluation; at run
     * time it becomes 0 and predict holds the sigma points.
     */
    constexpr UnscentedKalmanFilter(const Model& model, const UnscentedKalmanFilterConfig<Model, T>& config)
        : model_(model), x(config.x0), P(config.P0), Q(config.Q), Ts(detail::filter_sample_time(config.Ts)), p(config.unscented) {
        compute_weights();
    }

    /**
     * @brief Predict step: propagate the state and covariance through the dynamics over Ts
     *
     *     {X_i}     = sigma points of (x, P)
     *     X_i[k+1]  = RK4 step of f from X_i over Ts, u held
     *     x[k+1|k]  = Σ Wm_i X_i[k+1]
     *     P[k+1|k]  = Σ Wc_i (X_i − x)(X_i − x)ᵀ + Q
     *
     * @param u Input held over the step; the next updates' h sees it too
     * @return false if the prior covariance is not positive semidefinite (sigma points cannot be
     *         drawn); the estimate and P are unchanged then
     */
    [[nodiscard]] constexpr bool predict(const ColVec<NU, T>& u) {
        u_ = u;
        const auto sigma_opt = sigma_points();
        if (!sigma_opt) {
            return false;
        }
        const auto& sig = sigma_opt.value();

        damp::array<ColVec<NX, T>, NSIG> Xp{};
        for (size_t i = 0; i < NSIG; ++i) {
            Xp[i] = detail::rk4_step<NX, NU, T>(model_, sig[i], u, Ts);
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
     * @brief Measurement update with one reading y = h(x, u) + v of any size, u the input of the
     *        last predict
     *
     *     {X_i}  = sigma points of (x, P)
     *     Y_i    = h(X_i, u)
     *     ŷ      = Σ Wm_i Y_i
     *     Pyy    = Σ Wc_i (Y_i − ŷ)(Y_i − ŷ)ᵀ + R   (Cholesky: Pyy is symmetric positive definite)
     *     Pxy    = Σ Wc_i (X_i − x)(Y_i − ŷ)ᵀ
     *     K      = Pxy Pyy⁻¹
     *     x[k|k] = x + K (y − ŷ)
     *     P[k|k] = P − K Pyy Kᵀ
     *
     * @param h Measurement y = h(x, u) (Measurement), a struct or a lambda
     * @param y Reading (NY rows)
     * @param R Covariance of the reading's noise
     * @return false if P is not positive semidefinite or Pyy not positive definite: the estimate
     *         and P are unchanged
     */
    template<typename H, size_t NY>
        requires Measurement<H, NX, NU, NY, T>
    [[nodiscard]] constexpr bool update(const H& h, const ColVec<NY, T>& y, const damp::type_identity_t<Matrix<NY, NY, T>>& R) {
        const auto m = measurement_moments<NY>(h, R);
        if (!m) {
            return false;
        }
        // K = Pxy Pyy⁻¹ → solve Pyy Kᵀ = Pxyᵀ by Cholesky (fail closed if Pyy is not positive definite)
        const auto Kt_opt = mat::cholesky_solve(m->Pyy, Matrix<NY, NX, T>(m->Pxy.transpose()));
        if (!Kt_opt) {
            return false;
        }
        const Matrix<NX, NY, T> K = Kt_opt.value().transpose();
        x = x + (K * ColVec<NY, T>(y - m->y_pred));
        // P − K Pxyᵀ − Pxy Kᵀ + K Pyy Kᵀ: equal to P − K Pyy Kᵀ when K = Pxy Pyy⁻¹,
        // but (like the Joseph form) only second-order sensitive to rounding in K.
        const Matrix<NX, NX, T> P_post = P - (K * m->Pxy.t()) - (m->Pxy * K.t()) + quadratic_form(K, m->Pyy);
        P = static_cast<T>(0.5) * (P_post + P_post.transpose());
        return true;
    }

    /** @brief The state estimate x̂. */
    [[nodiscard]] constexpr const ColVec<NX, T>& state() const { return x; }
    /** @brief The estimate covariance P. */
    [[nodiscard]] constexpr const Matrix<NX, NX, T>& covariance() const { return P; }
    /** @brief The process-noise covariance Q. */
    [[nodiscard]] constexpr const Matrix<NX, NX, T>& process_noise_covariance() const { return Q; }
    /** @brief The process model. */
    [[nodiscard]] constexpr const Model& model() const { return model_; }
    /** @brief The unscented transform's tuning. */
    [[nodiscard]] constexpr const UnscentedParams<T>& params() const { return p; }
    /** @brief Sample time [s] that predict integrates over. */
    [[nodiscard]] constexpr T sample_time() const { return Ts; }

    /** @brief Overwrite the estimate between updates (a clamp, an angle wrap); P and the held input are kept. */
    constexpr void set_state(const ColVec<NX, T>& x_new) { x = x_new; }

    /** @brief Re-initialize the estimate and its covariance and forget the held input. The model, Ts, Q and the weights stay. */
    constexpr void reset(const ColVec<NX, T>& x0, const Matrix<NX, NX, T>& P0) {
        x = x0;
        P = P0;
        u_ = ColVec<NU, T>{};
    }

private:
    /** @brief Predicted measurement ŷ, innovation covariance Pyy and cross covariance Pxy of one update. */
    template<size_t NY>
    struct MeasurementMoments {
        ColVec<NY, T>     y_pred{}; ///< ŷ = Σ Wm_i Y_i
        Matrix<NY, NY, T> Pyy{};    ///< Σ Wc_i (Y_i − ŷ)(Y_i − ŷ)ᵀ + R, symmetrized
        Matrix<NX, NY, T> Pxy{};    ///< Σ Wc_i (X_i − x)(Y_i − ŷ)ᵀ
    };

    /** @brief Draw the sigma points, carry them through h, and form ŷ, Pyy, Pxy; nullopt if P is not positive semidefinite. */
    template<size_t NY, typename H>
    [[nodiscard]] constexpr damp::optional<MeasurementMoments<NY>> measurement_moments(const H& h, const Matrix<NY, NY, T>& R) const {
        const auto sigma_opt = sigma_points();
        if (!sigma_opt) {
            return damp::nullopt;
        }
        const auto& sig = sigma_opt.value();

        damp::array<ColVec<NY, T>, NSIG> Yp{};
        for (size_t i = 0; i < NSIG; ++i) {
            Yp[i] = h(sig[i], u_);
        }

        MeasurementMoments<NY> m{};
        m.y_pred = weighted_mean<NY>(Yp);
        m.Pyy = R;
        m.Pxy = Matrix<NX, NY, T>::zeros();
        for (size_t i = 0; i < NSIG; ++i) {
            const ColVec<NY, T> dy = Yp[i] - m.y_pred;
            const ColVec<NX, T> dx = sig[i] - x;
            m.Pyy = m.Pyy + (dy * dy.t()) * Wc[i];
            m.Pxy = m.Pxy + (dx * dy.t()) * Wc[i];
        }

        // Re-symmetrize Pyy before the factorization (sigma-point sum can skew off-diagonals).
        m.Pyy = static_cast<T>(0.5) * (m.Pyy + m.Pyy.transpose());
        return m;
    }

    /** @brief Compute the scaled-unscented-transform mean/covariance weights. */
    constexpr void compute_weights() {
        const T n = static_cast<T>(NX);
        lambda = (p.alpha * p.alpha * (n + p.kappa)) - n;
        // n + λ = α²(n + κ) must be positive (α = 0 or κ ≤ −n divides by zero):
        // fall back to the plain unscented transform, λ = 0.
        if (!(n + lambda > T{0})) {
            lambda = T{0};
        }
        const T denom = n + lambda;
        Wm[0] = lambda / denom;
        Wc[0] = (lambda / denom) + (static_cast<T>(1) - (p.alpha * p.alpha) + p.beta);
        const T w = static_cast<T>(1) / (static_cast<T>(2) * denom);
        for (size_t i = 1; i < NSIG; ++i) {
            Wm[i] = w;
            Wc[i] = w;
        }
    }

    /**
     * @brief L Lᵀ = S for a symmetric positive-semidefinite S.
     *
     * The Cholesky recursion, with a pivot within roundoff of zero (mat::cholesky's floor, either
     * sign) zeroing its column: a direction P does not spread in gets no sigma-point spread. An
     * entry below a zeroed pivot must then be within the bound |Sᵢⱼ|² ≤ Sᵢᵢ Sⱼⱼ of a
     * semidefinite matrix. nullopt on a clearly negative pivot, a violated bound or a non-finite
     * entry.
     */
    [[nodiscard]] static constexpr damp::optional<Matrix<NX, NX, T>> semidefinite_factor(const Matrix<NX, NX, T>& S) {
        Matrix<NX, NX, T>  L{};
        damp::array<T, NX> floor{};
        for (size_t j = 0; j < NX; ++j) {
            T pivot = S(j, j);
            T mag = damp::abs(S(j, j));
            for (size_t k = 0; k < j; ++k) {
                pivot -= L(j, k) * L(j, k);
                mag += L(j, k) * L(j, k);
            }
            floor[j] = T{16} * static_cast<T>(j + 1) * eps<T>() * mag;
            if (!damp::isfinite(pivot) || pivot < -floor[j]) {
                return damp::nullopt;
            }
            L(j, j) = (pivot > floor[j]) ? damp::sqrt(pivot) : T{0};
            for (size_t i = j + 1; i < NX; ++i) {
                T sum = S(i, j);
                for (size_t k = 0; k < j; ++k) {
                    sum -= L(i, k) * L(j, k);
                }
                if (L(j, j) > T{0}) {
                    L(i, j) = sum / L(j, j);
                } else if (sum * sum > T{4} * floor[j] * damp::abs(S(i, i))) {
                    return damp::nullopt; // not semidefinite
                }
            }
        }
        return L;
    }

    /** @brief Draw the 2n+1 sigma points from (x, P); nullopt if P is not positive semidefinite. */
    [[nodiscard]] constexpr damp::optional<damp::array<ColVec<NX, T>, NSIG>> sigma_points() const {
        const auto L_opt = semidefinite_factor(P * (static_cast<T>(NX) + lambda));
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

    /** @brief Weighted (mean-weight) combination of a sigma-point set. */
    template<size_t N>
    [[nodiscard]] constexpr ColVec<N, T> weighted_mean(const damp::array<ColVec<N, T>, NSIG>& pts) const {
        ColVec<N, T> m{};
        for (size_t i = 0; i < NSIG; ++i) {
            m = m + pts[i] * Wm[i];
        }
        return m;
    }

    Model              model_{}; ///< Continuous dynamics ẋ = f(x, u)
    ColVec<NX, T>      x{};      ///< Estimate
    Matrix<NX, NX, T>  P{};      ///< Estimate covariance
    Matrix<NX, NX, T>  Q{};      ///< Process-noise covariance per step
    T                  Ts{};     ///< Sample time [s]; 0 (default-constructed) holds the sigma points in predict
    UnscentedParams<T> p{};      ///< Unscented transform tuning
    ColVec<NU, T>      u_{};     ///< The input the last predict held, for the readings' h

    T                    lambda{}; ///< λ = α²(n + κ) − n
    damp::array<T, NSIG> Wm{};     ///< Mean weights
    damp::array<T, NSIG> Wc{};     ///< Covariance weights
};

/** @brief Deduce the filter from its model, a struct or a captureless lambda: `UnscentedKalmanFilter ukf{f, {.Ts = Ts, .Q = Q}}`. */
template<typename Model>
UnscentedKalmanFilter(const Model&, const UnscentedKalmanFilterConfig<Model>&) -> UnscentedKalmanFilter<Model>;

/**
 * @brief Short alias for @ref UnscentedKalmanFilter (MATLAB®-style name)
 *
 * @see UnscentedKalmanFilter
 */
template<typename Model, typename T = detail::model_shape<Model>::scalar>
using UKF = UnscentedKalmanFilter<Model, T>;

} // namespace damp
