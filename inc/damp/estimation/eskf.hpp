// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file eskf.hpp
 * @brief Error-state Kalman filter (ESKF) for attitude — design + runtime core
 *
 * Indirect (error-state) filter: the estimator tracks a small correction δx
 * while the full nominal state (quaternion, biases) is integrated outside the
 * KF. That avoids quaternion singularities and keeps the linearization valid
 * for small errors.
 *
 * Typical 6-state attitude error (Solà):
 * @f[
 *   \delta x = \bigl[\delta\theta^\top,\; \delta b_g^\top\bigr]^\top
 *   \in \mathbb{R}^{6}
 * @f]
 * with right-multiplicative attitude error on the body.
 *
 * @see Solà et al., "Quaternion kinematics for the error-state Kalman filter" (2017)
 * @see "Optimal State Estimation" (Simon, 2006), §14.2
 * @see ESKFOrientationFilter for a turnkey IMU/MARG runtime
 * @see ins_eskf.hpp for the 15-state navigation ESKF
 */

#include <concepts>
#include <cstddef>

#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "ekf.hpp"

namespace damp {
namespace design {
/**
 * @struct ESKFResult
 * @brief Design payload for an ESKF: Q, R, P₀, and success
 *
 * Pure data for @ref ErrorStateKalmanFilter / @ref ESKFOrientationFilter.
 * Use @c .as\<float\>() before embedding on target.
 *
 * @tparam NDX Error-state size (6 for attitude + gyro bias)
 * @tparam NY  Measurement size (3 = IMU accel, 6 = MARG accel+mag)
 */
template<size_t NDX, size_t NY, typename T = double>
struct ESKFResult {
    Matrix<NDX, NDX, T> Q{};            ///< Process noise covariance Q
    Matrix<NY, NY, T>   R{};            ///< Measurement noise covariance R
    Matrix<NDX, NDX, T> P0{};           ///< Initial error covariance P₀
    bool                success{false}; ///< true if densities/dt were valid


    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return ESKFResult<NDX, NY, U>{
            Q.template as<U>(),
            R.template as<U>(),
            P0.template as<U>(),
            success
        };
    }
};

namespace detail {

/// True if every entry is finite and >= 0 (densities, stds, dt).
template<typename T>
[[nodiscard]] constexpr bool nonneg_finite(T x) {
    return damp::finite_non_negative(x);
}

/**
 * Discrete process Q for attitude + gyro-bias RW (6×6: δθ then δb_g).
 *
 * First-order continuous–discrete form (σ_g [rad/s/√Hz], σ_bg [rad/s^{3/2}]):
 * @f[
 *   Q_{\theta\theta}=\sigma_g^2\Delta t + \sigma_{bg}^2\Delta t^3/3,\;
 *   Q_{\theta b}=-\sigma_{bg}^2\Delta t^2/2,\;
 *   Q_{bb}=\sigma_{bg}^2\Delta t
 * @f]
 * Matches @f$ \dot{\delta\theta}=-\delta b_g+\ldots @f$ integrated over fixed Δt.
 */
template<typename T>
constexpr void fill_attitude_bias_q(
    Matrix<6, 6, T>& Q,
    T                gyro_noise_density,
    T                gyro_bias_rw,
    T                dt
) {
    const T sg2 = gyro_noise_density * gyro_noise_density;
    const T sbg2 = gyro_bias_rw * gyro_bias_rw;
    const T dt2 = dt * dt;
    const T dt3 = dt2 * dt;
    const T q_tt = sg2 * dt + sbg2 * (dt3 / T{3});
    const T q_tb = -sbg2 * (dt2 / T{2});
    const T q_bb = sbg2 * dt;
    for (size_t i = 0; i < 3; ++i) {
        Q(i, i) = q_tt;
        Q(i + 3, i + 3) = q_bb;
        Q(i, i + 3) = q_tb;
        Q(i + 3, i) = q_tb;
    }
}

template<size_t NY, typename T>
[[nodiscard]] constexpr ESKFResult<6, NY, T> eskf_fill_qp0(
    T gyro_noise_density,
    T accel_noise_density,
    T gyro_bias_rw,
    T dt,
    T initial_attitude_uncertainty,
    T initial_bias_uncertainty,
    T mag_noise_density = T{0}
) {
    ESKFResult<6, NY, T> result{};
    result.success = nonneg_finite(dt) && (dt > T{0}) && nonneg_finite(gyro_noise_density)
                  && nonneg_finite(accel_noise_density) && nonneg_finite(gyro_bias_rw)
                  && nonneg_finite(initial_attitude_uncertainty) && nonneg_finite(initial_bias_uncertainty)
                  && nonneg_finite(mag_noise_density);
    if (!result.success) {
        return result;
    }

    // Discrete measurement variance from continuous density: R ≈ σ² / Δt (Hz).
    const T accel_var = (accel_noise_density * accel_noise_density) / dt;
    const T att0 = initial_attitude_uncertainty * initial_attitude_uncertainty;
    const T b0 = initial_bias_uncertainty * initial_bias_uncertainty;

    result.Q = Matrix<6, 6, T>::zeros();
    result.P0 = Matrix<6, 6, T>::zeros();
    result.R = Matrix<NY, NY, T>::zeros();
    fill_attitude_bias_q(result.Q, gyro_noise_density, gyro_bias_rw, dt);
    for (size_t i = 0; i < 3; ++i) {
        result.P0(i, i) = att0;
        result.P0(i + 3, i + 3) = b0;
        result.R(i, i) = accel_var;
    }
    if constexpr (NY == 6) {
        const T mag_var = (mag_noise_density * mag_noise_density) / dt;
        for (size_t i = 0; i < 3; ++i) {
            result.R(i + 3, i + 3) = mag_var;
        }
    }
    return result;
}

} // namespace detail

/**
 * @brief Design Q, R, P₀ for a 6-state IMU attitude ESKF (gyro + accel)
 *
 * Error state is δx = [δθ; δb_g]. Builds discrete process noise Q, 3×3
 * specific-force measurement noise R, and prior P₀.
 *
 * Observability: pitch/roll from gravity (specific force); yaw is free without
 * a heading aid (e.g. steel structures / excavator bases where mag is unreliable).
 *
 * Discrete Q uses the continuous–discrete bias random-walk form
 * (@ref detail::fill_attitude_bias_q). Measurement R is
 * @f$ \sigma_a^2 / \Delta t @f$ from accel density [m/s²/√Hz] at sample period @p dt.
 *
 * @param gyro_noise_density   Gyroscope noise density [rad/s/√Hz]
 * @param accel_noise_density  Accelerometer noise density [m/s²/√Hz]
 * @param gyro_bias_rw         Gyro bias random walk [rad/s^{3/2}]
 * @param dt                   Sample period [s] (forms discrete Q and R)
 * @param initial_attitude_uncertainty  Initial attitude 1-σ [rad]
 * @param initial_bias_uncertainty      Initial bias 1-σ [rad/s]
 * @return @c ESKFResult<6, 3, T> for @ref ESKFOrientationFilter with NY=3
 *
 * @note Compare with Madgwick's 6-DOF IMU filter; Solà IMU kinematics + accel aid.
 * @see eskf_marg for accel + magnetometer (MARG)
 * @see ESKFOrientationFilter
 */
template<typename T = double>
[[nodiscard]] constexpr ESKFResult<6, 3, T> eskf_imu(
    T gyro_noise_density,
    T accel_noise_density,
    T gyro_bias_rw,
    T dt,
    T initial_attitude_uncertainty = static_cast<T>(0.1),
    T initial_bias_uncertainty = static_cast<T>(0.01)
) {
    return detail::eskf_fill_qp0<3>(
        gyro_noise_density,
        accel_noise_density,
        gyro_bias_rw,
        dt,
        initial_attitude_uncertainty,
        initial_bias_uncertainty
    );
}

/**
 * @brief Design Q, R, P₀ for a 6-state MARG attitude ESKF (gyro + accel + mag)
 *
 * Same error state as @ref eskf_imu; measurement covariance is 6×6
 * (specific force 3 + magnetometer 3). Mag observes yaw when the local field
 * is usable.
 *
 * @param gyro_noise_density   Gyroscope noise density [rad/s/√Hz]
 * @param accel_noise_density  Accelerometer noise density [m/s²/√Hz]
 * @param mag_noise_density    Magnetometer noise density [1/√Hz] (same units as mag meas)
 * @param gyro_bias_rw         Gyro bias random walk [rad/s^{3/2}]
 * @param dt                   Sample period [s]
 * @param initial_attitude_uncertainty  Initial attitude 1-σ [rad]
 * @param initial_bias_uncertainty      Initial bias 1-σ [rad/s]
 * @return @c ESKFResult<6, 6, T> for @ref ESKFOrientationFilter with NY=6
 *
 * @note Compare with Madgwick MARG / a 9-axis AHRS; Solà ESKF + vector aids.
 * @see eskf_imu for mag-free (IMU-only) mounts
 * @see ESKFOrientationFilter
 */
template<typename T = double>
[[nodiscard]] constexpr ESKFResult<6, 6, T> eskf_marg(
    T gyro_noise_density,
    T accel_noise_density,
    T mag_noise_density,
    T gyro_bias_rw,
    T dt,
    T initial_attitude_uncertainty = static_cast<T>(0.1),
    T initial_bias_uncertainty = static_cast<T>(0.01)
) {
    return detail::eskf_fill_qp0<6>(
        gyro_noise_density,
        accel_noise_density,
        gyro_bias_rw,
        dt,
        initial_attitude_uncertainty,
        initial_bias_uncertainty,
        mag_noise_density
    );
}

} // namespace design

/**
 * @brief ESKF prediction Jacobians F and G (nominal state updated outside)
 *
 * @tparam T   Scalar type
 * @tparam NDX Error state dimension
 */
template<typename T, size_t NDX>
struct ErrorStateJacobian {
    Matrix<NDX, NDX, T> F{}; ///< δx transition: ∂(δx[k+1])/∂(δx[k])
    Matrix<NDX, NDX, T> G{}; ///< Process-noise input (often I for additive Q)
};

/**
 * @brief Callable that returns @ref ErrorStateJacobian for one predict step
 *
 * Signature: (dt) → ErrorStateJacobian. Nominal integration (e.g. gyro on q)
 * is the caller's job; only F, G are returned here.
 */
template<typename Fn, typename T, size_t NDX>
concept ESKFPredictFn = requires(Fn&& fn, T dt) {
    { fn(dt) } -> std::convertible_to<ErrorStateJacobian<T, NDX>>;
};

/**
 * @brief Callable that returns a measurement linearization for one update
 *
 * Signature: () → MeasJacobian (captures nominal state). Yields predicted
 * measurement and H (and optional M).
 */
template<typename Fn, typename T, size_t NDX, size_t NY>
concept ESKFMeasFn = requires(Fn&& fn) {
    { fn() } -> std::convertible_to<MeasJacobian<T, NY, NDX>>;
};

/**
 * @brief Runtime error-state KF: tracks δx and P, not the full nominal state
 *
 * Estimates a small error δx; the user keeps the nominal quaternion/biases
 * and injects corrections after @ref update. Typical layout:
 * δx = [δθ (3), δb_g (3), …].
 *
 * Workflow:
 * 1. Integrate nominal state (e.g. q ← q ⊗ Δq from gyro).
 * 2. @ref predict — @f$ P \leftarrow F P F^\top + G Q G^\top @f$.
 * 3. @ref update — @f$ \delta x = K(y - h(x_{\mathrm{nom}})) @f$.
 * 4. Inject δx into the nominal state and @ref reset_error_state.
 *
 * @see Solà et al., "Quaternion kinematics for the error-state Kalman filter" (2017), §5–6
 *
 * @tparam NDX Error state dimension
 * @tparam NY  Default measurement dimension for @ref update
 * @tparam T   Scalar (default float for embedded runtime)
 */
template<size_t NDX, size_t NY, typename T = float>
struct ErrorStateKalmanFilter {
    // Error state dimension typically: 3 (attitude) + 3 (gyro bias) + 3 (accel bias) = 9
    // Nominal state (quaternion + biases) is tracked separately and not part of the KF

    constexpr ErrorStateKalmanFilter() = default;

    constexpr ErrorStateKalmanFilter(
        const Matrix<NDX, NDX, T>& P0,
        const Matrix<NDX, NDX, T>& Q_,
        const Matrix<NY, NY, T>&   R_
    ) : P(P0), Q(Q_), R(R_), delta_x(ColVec<NDX, T>{}) {}

    // Type conversion constructor
    template<typename U>
    constexpr ErrorStateKalmanFilter(const ErrorStateKalmanFilter<NDX, NY, U>& other)
        : P(other.covariance()),
          Q(other.process_noise_covariance()),
          R(other.measurement_noise_covariance()),
          delta_x(other.error_state()),
          innov(other.innovation()) {}

    /// From design result (any scalar); converts P0/Q/R via @c .as\<T\>().
    template<typename U>
    constexpr ErrorStateKalmanFilter(const design::ESKFResult<NDX, NY, U>& result)
        : P(result.P0.template as<T>()),
          Q(result.Q.template as<T>()),
          R(result.R.template as<T>()),
          delta_x(ColVec<NDX, T>{}) {}

    // Predict: Propagate error covariance (nominal state updated externally by user)
    // User provides: (dt) -> ErrorStateJacobian{F, G}
    //   F = error state transition Jacobian (∂δx_next/∂δx)
    //   G = process noise Jacobian (maps Q to error covariance)
    // P[k+1|k] = F * P[k|k] * F' + G * Q * G'
    template<typename PredictFn>
        requires ESKFPredictFn<PredictFn, T, NDX>
    constexpr void predict(const PredictFn& propagate_nominal, const T dt) {
        ErrorStateJacobian<T, NDX> ej = propagate_nominal(dt);
        // P ← F P Fᵀ + G Q Gᵀ with exact symmetry (quadratic_form).
        // Typical IMU/INS paths use G = I (additive discrete Q).
        P = quadratic_form(ej.F, P);
        bool g_is_I = true;
        for (size_t i = 0; i < NDX && g_is_I; ++i) {
            for (size_t j = 0; j < NDX; ++j) {
                const T want = (i == j) ? T{1} : T{0};
                if (ej.G(i, j) != want) {
                    g_is_I = false;
                    break;
                }
            }
        }
        if (g_is_I) {
            // Q is assumed SPD/symmetric by the caller; scrub any stored asymmetry.
            P = P + static_cast<T>(0.5) * (Q + Q.transpose());
        } else {
            P = P + quadratic_form(ej.G, Q);
        }
        // After inject, δx ≈ 0; multi-predict without inject leaves δx unpropagated by design.
    }

    /**
     * @brief Joseph-form update from a precomputed innovation (any NY2)
     *
     * Use when this update's size differs from the filter's default NY
     * (e.g. scalar heading while stored R is 3×3 for position).
     *
     * @tparam NY2 Measurement dimension for this update
     * @param innovation  y − h(x) (form and wrap angles before calling)
     * @param H           ∂h/∂δx (NY2 × NDX)
     * @param R_meas      measurement covariance for this update
     * @return false if the innovation covariance is not SPD / solve fails
     */
    template<size_t NY2>
    [[nodiscard]] constexpr bool update_innovation(
        const ColVec<NY2, T>&      innovation,
        const Matrix<NY2, NDX, T>& H,
        const Matrix<NY2, NY2, T>& R_meas
    ) {
        // S = HPHᵀ + R with exact symmetry (quadratic_form); fail closed if singular.
        const Matrix<NY2, NY2, T> S = quadratic_form(H, P) + R_meas;
        auto                      K_opt = mat::cholesky_solve(S, H * P);
        if (!K_opt) {
            K_opt = mat::lu_solve(S, H * P);
        }
        if (!K_opt) {
            return false;
        }
        const Matrix K = K_opt.value().transpose(); // NDX × NY2
        delta_x = ColVec(K * innovation);
        // Joseph form: P = (I − KH) P (I − KH)ᵀ + K R Kᵀ (exact P symmetry)
        const auto I_KH = Matrix<NDX, NDX, T>::identity() - K * H;
        P = quadratic_form(I_KH, P) + quadratic_form(K, R_meas);
        if constexpr (NY2 == NY) {
            innov = innovation;
        }
        return true;
    }

    // Measurement update: Correct error state from innovation
    // User provides: () -> MeasJacobian{y_pred, H, M}
    //   y_pred = h(nominal_state) - predicted measurement
    //   H = measurement Jacobian (∂h/∂δx)
    //   M = measurement noise Jacobian (∂h/∂v, usually identity)
    // P update: P = (I - K*H) * P * (I - K*H)' + K*M*R*M'*K'  (Joseph form)
    template<typename MeasFn>
        requires ESKFMeasFn<MeasFn, T, NDX, NY>
    constexpr bool update(const MeasFn& meas_fn, const ColVec<NY, T>& y) {
        MeasJacobian<T, NY, NDX> mj = meas_fn();

        innov = y - mj.y_pred; // Innovation
        const Matrix<NY, NY, T> R_eff = mj.M * R * mj.M.t();
        return update_innovation(innov, mj.H, R_eff);
    }

    // Reset error state after applying corrections to nominal state
    // Call this after you've updated the nominal quaternion/biases with delta_x
    // Optionally provide a G matrix to adjust covariance (e.g., for attitude reset)
    constexpr void reset_error_state(const Matrix<NDX, NDX, T>& G = Matrix<NDX, NDX, T>::identity()) {
        delta_x = ColVec<NDX, T>{};
        P = G * P * G.t();
    }

    // Accessors
    [[nodiscard]] constexpr const auto& error_state() const { return delta_x; }
    [[nodiscard]] constexpr const auto& covariance() const { return P; }
    [[nodiscard]] constexpr const auto& process_noise_covariance() const { return Q; }
    [[nodiscard]] constexpr const auto& measurement_noise_covariance() const { return R; }
    [[nodiscard]] constexpr const auto& innovation() const { return innov; }

    // Setters for runtime tuning
    constexpr void set_process_noise_covariance(const Matrix<NDX, NDX, T>& Q_new) { Q = Q_new; }
    constexpr void set_measurement_noise_covariance(const Matrix<NY, NY, T>& R_new) { R = R_new; }
    constexpr void set_covariance(const Matrix<NDX, NDX, T>& P_new) { P = P_new; }

    // Error-state mutators. The ESKF estimates the *correction* δx; the nominal
    // state (quaternion, biases) lives outside the filter. Use these to bound or
    // zero a correction component before injecting it into the nominal state —
    // e.g. cap an attitude-error step so a bad measurement can't slew the
    // quaternion, or freeze a bias correction that has saturated. Apply between
    // update() and your nominal-state injection + reset_error_state().
    constexpr void set_error_state(const ColVec<NDX, T>& dx_new) { delta_x = dx_new; }
    constexpr void set_error_state(size_t i, T value) { delta_x[i] = value; }

private:
    Matrix<NDX, NDX, T> P{};       // Error covariance
    Matrix<NDX, NDX, T> Q{};       // Process noise covariance
    Matrix<NY, NY, T>   R{};       // Measurement noise covariance
    ColVec<NDX, T>      delta_x{}; // Error state δx
    ColVec<NY, T>       innov{};   // Innovation
};

// CTAD deduction guide for ErrorStateKalmanFilter
template<typename T, size_t NDX, size_t NY>
ErrorStateKalmanFilter(Matrix<NDX, NDX, T>, Matrix<NDX, NDX, T>, Matrix<NY, NY, T>) -> ErrorStateKalmanFilter<NDX, NY, T>;

} // namespace damp