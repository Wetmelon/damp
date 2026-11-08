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
 * @see "Quaternion kinematics for the error-state Kalman filter" (Solà, 2017), §5–6
 * @see ESKFOrientationFilter for a turnkey IMU/MARG runtime
 * @see ins_eskf.hpp for the 15-state navigation ESKF
 */

#include <concepts>
#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/math/geometry.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/pattern_matrix.hpp"

namespace damp {

namespace detail {

/** @brief Skew-symmetric (cross-product) matrix [v]×, so that [v]× w = v × w. */
template<typename T>
[[nodiscard]] constexpr Matrix<3, 3, T> skew_sym(const Vec3<T>& v) {
    return Matrix<3, 3, T>{
        {T{0}, -v[2], v[1]},
        {v[2], T{0}, -v[0]},
        {-v[1], v[0], T{0}},
    };
}

/** @brief Right-multiplicative attitude perturbation exp(δθ) as a unit quaternion. */
template<typename T>
[[nodiscard]] constexpr Quaternion<T> exp_delta_theta(const Vec3<T>& dth) {
    const T n = dth.norm();
    if (!(n > T{0})) {
        return Quaternion<T>::identity();
    }
    return Quaternion<T>::from_axis_angle(dth, n, T{0}).value_or(Quaternion<T>::identity());
}

/**
 * Error-state reset Jacobian for a right-multiplicative attitude error after
 * injecting δθ: G_θθ = I − ½[δθ]× (Solà 2017, §6.3).
 */
template<typename T>
[[nodiscard]] constexpr Matrix<3, 3, T> attitude_reset_jacobian(const Vec3<T>& dth) {
    return Matrix<3, 3, T>::identity() - (skew_sym(dth) * T{0.5});
}

} // namespace detail

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
    Matrix<NY, NY, T>   R{};            ///< Measurement noise covariance R of the update (update_innovation(ν, H, R))
    Matrix<NDX, NDX, T> P0{};           ///< Initial error covariance P₀
    bool                success{false}; ///< true if densities/dt were valid

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

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

/** @brief True if every entry is finite and >= 0 (densities, stds, dt). */
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
    // R must be positive definite: the attitude Jacobians [a]× and [m]× have rank 2,
    // so a zero density leaves S = H P Hᵀ + R singular and every update fails.
    result.success = nonneg_finite(dt) && (dt > T{0}) && nonneg_finite(gyro_noise_density)
                  && nonneg_finite(accel_noise_density) && (accel_noise_density > T{0})
                  && nonneg_finite(gyro_bias_rw) && nonneg_finite(initial_attitude_uncertainty)
                  && nonneg_finite(initial_bias_uncertainty) && nonneg_finite(mag_noise_density)
                  && (NY != 6 || mag_noise_density > T{0});
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
 * (@c detail::fill_attitude_bias_q). Measurement R is
 * @f$ \sigma_a^2 / \Delta t @f$ from accel density [m/s²/√Hz] at sample period @p dt.
 *
 * @param gyro_noise_density   Gyroscope noise density [rad/s/√Hz]
 * @param accel_noise_density  Accelerometer noise density [m/s²/√Hz], > 0
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
 * @param accel_noise_density  Accelerometer noise density [m/s²/√Hz], > 0
 * @param mag_noise_density    Magnetometer noise density [1/√Hz] (same units as mag meas), > 0
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
 * @brief Callable that returns @ref damp::ErrorStateJacobian "ErrorStateJacobian" for one predict step
 *
 * Signature: (dt) → ErrorStateJacobian. Nominal integration (e.g. gyro on q)
 * is the caller's job; only F, G are returned here.
 */
template<typename Fn, typename T, size_t NDX>
concept ESKFPredictFn = requires(Fn&& fn, T dt) {
    { fn(dt) } -> std::convertible_to<ErrorStateJacobian<T, NDX>>;
};


/**
 * @brief Configuration of an ErrorStateKalmanFilter: the process noise and the initial covariance
 *
 * An aggregate; name the fields you set: `ErrorStateKalmanFilter<15, float> eskf{{.Q = Q, .P0 = P0}}`.
 *
 * @tparam NDX Error-state size
 */
template<size_t NDX, typename T = float>
struct ErrorStateKalmanFilterConfig {
    Matrix<NDX, NDX, T> Q{Matrix<NDX, NDX, T>::identity()};  ///< Process-noise covariance per step (through G)
    Matrix<NDX, NDX, T> P0{Matrix<NDX, NDX, T>::identity()}; ///< Initial error covariance
};

/**
 * @brief Runtime error-state KF: tracks δx and P, not the full nominal state
 *
 * Estimates a small error δx; the user keeps the nominal quaternion/biases
 * and injects corrections after @ref update_innovation. Typical layout:
 * δx = [δθ (3), δb_g (3), …].
 *
 * Workflow:
 * 1. Integrate nominal state (e.g. q ← q ⊗ Δq from gyro).
 * 2. @ref predict — @f$ P \leftarrow F P F^\top + G Q G^\top @f$.
 * 3. @ref update_innovation with each reading, linearized at the nominal state —
 *    @f$ \delta x \leftarrow \delta x + K(\nu - H\,\delta x) @f$.
 * 4. Inject δx into the nominal state and @ref reset_error_state.
 *
 * The measurement model lives outside the filter with the nominal state, so each reading
 * arrives as its innovation ν = y − h(x_nom), its Jacobian H = ∂h/∂δx and its covariance R, of
 * any size.
 *
 * @see "Quaternion kinematics for the error-state Kalman filter" (Solà, 2017), §5–6
 *
 * @tparam NDX Error state dimension
 * @tparam T   Scalar (default float for embedded runtime)
 */
template<size_t NDX, typename T = float>
struct ErrorStateKalmanFilter {
    // Error state dimension typically: 3 (attitude) + 3 (gyro bias) + 3 (accel bias) = 9
    // Nominal state (quaternion + biases) is tracked separately and not part of the KF

    constexpr ErrorStateKalmanFilter() = default;

    /** @brief Filter from its process noise and initial covariance. */
    constexpr explicit ErrorStateKalmanFilter(const ErrorStateKalmanFilterConfig<NDX, T>& config) : P(config.P0), Q(config.Q), q_is_diagonal(is_diagonal(config.Q)) {}

    /** @brief The same filter in another scalar. */
    template<typename U>
    constexpr ErrorStateKalmanFilter(const ErrorStateKalmanFilter<NDX, U>& other) // NOLINT(google-explicit-constructor): a scalar conversion
        : P(other.covariance().template as<T>()),
          Q(other.process_noise_covariance().template as<T>()),
          delta_x(other.error_state().template as<T>()),
          q_is_diagonal(is_diagonal(Q)) {}

    /**
     * @brief From a design result (any scalar): its P₀ and Q via @c .as\<T\>(). The design's R
     *        is the update's: update_innovation(ν, H, res.R).
     */
    template<typename U, size_t NY>
    constexpr ErrorStateKalmanFilter(const design::ESKFResult<NDX, NY, U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : P(result.P0.template as<T>()),
          Q(result.Q.template as<T>()),
          q_is_diagonal(is_diagonal(Q)) {
        if (!result.success && std::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    /**
     * @brief Propagate the error mean and covariance over one step (the nominal state is
     *        integrated by the caller)
     *
     *     δx[k+1|k] = F δx[k|k];  P[k+1|k] = F P[k|k] Fᵀ + G Q Gᵀ
     *
     * @param propagate_nominal (dt) → ErrorStateJacobian{F, G}
     * @param dt                Step [s]
     */
    template<typename PredictFn>
        requires ESKFPredictFn<PredictFn, T, NDX>
    constexpr void predict(const PredictFn& propagate_nominal, const T dt) {
        ErrorStateJacobian<T, NDX> ej = propagate_nominal(dt);
        // After an inject δx = 0 and this is a no-op; a correction that has not
        // been injected yet must still be carried through the error dynamics.
        delta_x = ej.F * delta_x;
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
    }

    /**
     * @brief Propagate over one step with F's zeros and ones known at compile time (G = I):
     *        δx ← F δx, P ← F P Fᵀ + Q in place, with no term for a zero and no multiply for a one
     *
     * The same numbers as predict() with F dense, to rounding, on the exactly symmetric P the filter
     * keeps (P0 has to be symmetric). A diagonal Q adds only its diagonal.
     *
     * @tparam Pattern F's sparsity pattern
     * @param F The transition of δx, its runtime values filled in
     */
    template<const auto& Pattern>
    constexpr void predict(const PatternMatrix<Pattern, T>& F) {
        delta_x = F * delta_x;
        quadratic_form_in_place(F, P);
        if (q_is_diagonal) {
            for (size_t i = 0; i < NDX; ++i) {
                P(i, i) += Q(i, i);
            }
        } else {
            for (size_t i = 0; i < NDX; ++i) {
                for (size_t j = 0; j <= i; ++j) {
                    const T q = T{0.5} * (Q(i, j) + Q(j, i));
                    P(i, j) += q;
                    if (i != j) {
                        P(j, i) += q;
                    }
                }
            }
        }
    }

    /**
     * @brief Joseph-form update with one reading of any size, from its innovation
     *
     * ν̄ = ν − H δx is the residual against the corrected estimate, which is ν itself after
     * reset_error_state (δx = 0). S = H P Hᵀ + R is symmetric positive definite, so
     * K = P Hᵀ S⁻¹ comes from S Kᵀ = H P by Cholesky. Then δx ← δx + K ν̄ and
     * P ← (I − K H) P (I − K H)ᵀ + K R Kᵀ.
     *
     * @param innovation ν = y − h(x_nom), with any wrap the reading needs (a heading residual on (−π, π])
     * @param H          ∂h/∂δx at the nominal state (NY × NDX)
     * @param R          Covariance of the reading's noise
     * @return false if S is not positive definite: δx and P are unchanged
     */
    template<size_t NY>
    [[nodiscard]] constexpr bool update_innovation(
        const ColVec<NY, T>&                           innovation,
        const Matrix<NY, NDX, T>&                      H,
        const std::type_identity_t<Matrix<NY, NY, T>>& R
    ) {
        // S = HPHᵀ + R with exact symmetry (quadratic_form); fail closed if not positive definite.
        const Matrix<NY, NY, T> S = quadratic_form(H, P) + R;
        const auto              K_opt = mat::cholesky_solve(S, Matrix<NY, NDX, T>(H * P));
        if (!K_opt) {
            return false;
        }
        const Matrix<NDX, NY, T> K = K_opt.value().transpose();
        // Stack on any correction not yet injected. A single update after
        // reset_error_state (δx = 0) matches δx = K (y − h).
        const ColVec<NY, T> innov_bar = innovation - (H * delta_x);
        delta_x = delta_x + (K * innov_bar);
        // Joseph form: P = (I − KH) P (I − KH)ᵀ + K R Kᵀ (exact P symmetry)
        const Matrix<NDX, NDX, T> I_KH = Matrix<NDX, NDX, T>::identity() - (K * H);
        P = quadratic_form(I_KH, P) + quadratic_form(K, R);
        return true;
    }

    /**
     * @brief Joseph-form update with H's zeros and ones known at compile time: the same factors as
     *        update_innovation with H dense, with no term for what H does not read
     *
     * With A = I − K H, P ← A P Aᵀ + K R Kᵀ stays that sum of two positive semidefinite terms; it is
     * not the shorter P − K H P, which can lose symmetry and positive definiteness to rounding. Only
     * the entries of A P that the result reads are computed: its lower triangle and the columns H
     * reads. The result is exactly symmetric.
     *
     * @tparam Pattern H's sparsity pattern
     * @param innovation ν = y − h(x_nom), with any wrap the reading needs
     * @param H          ∂h/∂δx at the nominal state, its runtime values filled in
     * @param R          Covariance of the reading's noise
     * @return false if S = H P Hᵀ + R is not positive definite: δx and P are unchanged
     */
    template<const auto& Pattern, size_t NY>
    [[nodiscard]] constexpr bool update_innovation(
        const ColVec<NY, T>&                           innovation,
        const PatternMatrix<Pattern, T>&               H,
        const std::type_identity_t<Matrix<NY, NY, T>>& R
    ) {
        static_assert(PatternMatrix<Pattern, T>::NR == NY && PatternMatrix<Pattern, T>::NC == NDX, "H is NY × NDX.");
        const Matrix<NY, NDX, T> HP = H * P;
        const Matrix<NY, NY, T>  S = times_transpose(HP, H) + R;
        const auto               K_opt = mat::cholesky_solve(S, HP);
        if (!K_opt) {
            return false;
        }
        const Matrix<NDX, NY, T> K = K_opt.value().transpose();
        const ColVec<NY, T>      innov_bar = innovation - (H * delta_x);
        delta_x = delta_x + (K * innov_bar);

        const Matrix<NDX, NY, T> KR = K * R;
        constexpr auto&          read_by_H = ::damp::detail::pattern_read_columns<Pattern>;
        // A P, row by row: its lower triangle, and past the diagonal only the columns H reads.
        T AP[NDX][NDX]; // written before it is read, only where the result reads it
        for (size_t i = 0; i < NDX; ++i) {
            const auto ap_entry = [&](size_t l) __attribute__((always_inline)) {
                T s = P(i, l);
                [&]<size_t... M>(damp::index_sequence<M...>) __attribute__((always_inline)) {
                    ((s -= K(i, M) * HP(M, l)), ...);
                }(damp::make_index_sequence<NY>{});
                AP[i][l] = s;
            };
            for (size_t l = 0; l <= i; ++l) {
                ap_entry(l);
            }
            [&]<size_t... C>(damp::index_sequence<C...>) __attribute__((always_inline)) {
                ((read_by_H[C] && C > i ? ap_entry(C) : void()), ...);
            }(damp::make_index_sequence<NDX>{});
        }
        // P ← A P Aᵀ + K R Kᵀ, row by row: the part of each row that does not depend on the column once.
        for (size_t i = 0; i < NDX; ++i) {
            T d[NY];
            [&]<size_t... M>(damp::index_sequence<M...>) __attribute__((always_inline)) {
                ((d[M] = KR(i, M) - H.template row_dot<M>([&](auto c) __attribute__((always_inline)) { return AP[i][decltype(c)::value]; })), ...);
            }(damp::make_index_sequence<NY>{});
            for (size_t j = 0; j <= i; ++j) {
                T s = AP[i][j];
                [&]<size_t... M>(damp::index_sequence<M...>) __attribute__((always_inline)) {
                    ((s += K(j, M) * d[M]), ...);
                }(damp::make_index_sequence<NY>{});
                P(i, j) = s;
                P(j, i) = s;
            }
        }
        return true;
    }

    /**
     * @brief Zero the error state after its correction is injected into the nominal state
     *
     * G is the reset Jacobian (P ← G P Gᵀ); for a right-multiplicative attitude error its
     * attitude block is detail::attitude_reset_jacobian(δθ).
     */
    constexpr void reset_error_state(const Matrix<NDX, NDX, T>& G = Matrix<NDX, NDX, T>::identity()) {
        delta_x = ColVec<NDX, T>{};
        P = quadratic_form(G, P);
    }

    /** @brief Zero the error state, with the reset Jacobian's zeros and ones known at compile time: P ← G P Gᵀ in place. */
    template<const auto& Pattern>
    constexpr void reset_error_state(const PatternMatrix<Pattern, T>& G) {
        delta_x = ColVec<NDX, T>{};
        quadratic_form_in_place(G, P);
    }

    /** @brief The error-state estimate δx (zero after reset_error_state). */
    [[nodiscard]] constexpr const ColVec<NDX, T>& error_state() const { return delta_x; }
    /** @brief The error covariance P. */
    [[nodiscard]] constexpr const Matrix<NDX, NDX, T>& covariance() const { return P; }
    /** @brief The process-noise covariance Q. */
    [[nodiscard]] constexpr const Matrix<NDX, NDX, T>& process_noise_covariance() const { return Q; }

    /**
     * @brief Overwrite the correction before injecting it: cap an attitude-error step so a bad
     *        measurement cannot slew the quaternion, or freeze a saturated bias correction
     */
    constexpr void set_error_state(const ColVec<NDX, T>& dx_new) { delta_x = dx_new; }

    /**
     * @brief Replace the error state and its covariance: a re-initialization, or the result of a
     *        propagation done outside the filter. Q stays.
     */
    constexpr void reset(const ColVec<NDX, T>& dx0, const Matrix<NDX, NDX, T>& P0) {
        delta_x = dx0;
        P = P0;
    }

private:
    Matrix<NDX, NDX, T> P{};                 ///< Error covariance
    Matrix<NDX, NDX, T> Q{};                 ///< Process noise covariance
    ColVec<NDX, T>      delta_x{};           ///< Error state δx
    bool                q_is_diagonal{true}; ///< Q has no off-diagonal entries: predict adds only its diagonal

    [[nodiscard]] static constexpr bool is_diagonal(const Matrix<NDX, NDX, T>& m) {
        for (size_t i = 0; i < NDX; ++i) {
            for (size_t j = 0; j < NDX; ++j) {
                if (i != j && m(i, j) != T{0}) {
                    return false;
                }
            }
        }
        return true;
    }
};

/** @brief Deduce the filter from its design, keeping the result's scalar. */
template<size_t NDX, size_t NY, typename T>
ErrorStateKalmanFilter(const design::ESKFResult<NDX, NY, T>&) -> ErrorStateKalmanFilter<NDX, T>;

} // namespace damp
