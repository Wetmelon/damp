// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file sensor_fusion.hpp
 * @brief Orientation estimation filters for IMU sensor fusion
 *
 * Provides several filters of increasing sophistication:
 *   - ComplementaryFilter: simple α-blend of gyro and accelerometer
 *   - MadgwickFilter: gradient-descent AHRS (accel + mag)
 *   - MahonyFilter: nonlinear complementary filter with PI correction
 *   - ESKFOrientationFilter: Tier-3 attitude ESKF (IMU NY=3 or MARG NY=6)
 *
 * Attitude ESKF design is design::eskf_imu / design::eskf_marg; the
 * runtime tick is ESKFOrientationFilter::step (not free functions that
 * take a bare @c ErrorStateKalmanFilter).
 *
 * @see "An efficient orientation filter for inertial and inertial/magnetic
 *      sensor arrays" (Madgwick, 2010)
 * @see "Nonlinear Complementary Filters on the Special Orthogonal Group"
 *      (Mahony et al., 2008)
 * @see "Quaternion kinematics for the error-state Kalman filter" (Solà, 2017)
 */

#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/math/geometry.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/core.hpp"
#include "eskf.hpp"
#include "ins_mechanization.hpp"

namespace damp {

/**
 * @brief Simple complementary filter for orientation estimation
 *
 * Blends gyroscope integration (high-frequency) with accelerometer tilt
 * (low-frequency), crossing over at 1/τ: each tick slerps toward the accelerometer
 * tilt by dt/(τ + dt), i.e. α = τ/(τ + dt) as in the scalar @ref Complementary, so
 * the crossover does not depend on the sample rate. Yaw is carried by the gyro
 * alone (the accelerometer does not observe it).
 *
 * Conventions (shared by @ref MadgwickFilter and @ref MahonyFilter):
 * @c orientation() is body → world (@c q.rotate(v_body) is in the world frame,
 * z up), the same as @ref ESKFOrientationFilter and @ref InsState; @p accel is
 * the specific force in the body frame, the accelerometer output, as for
 * @ref ESKFOrientationFilter and @ref ImuSample — at rest and level it reads
 * (0, 0, +g).
 *
 * @tparam T Scalar type (default: float)
 */
template<typename T = float>
class ComplementaryFilter {
    T             tau{static_cast<T>(0.5)};
    Quaternion<T> q{Quaternion<T>::identity()};

public:
    constexpr ComplementaryFilter() = default;

    /** @param tau_s Crossover time constant τ [s] (larger ⇒ trust the gyro longer; default 0.5 s). */
    constexpr explicit ComplementaryFilter(T tau_s) : tau(tau_s) {}

    /**
     * @brief One tick: integrate @p gyro over @p dt, then blend toward the accelerometer tilt
     *
     * A zero or non-finite @p accel skips the blend (gyro only), as the
     * @ref ESKFOrientationFilter accel gate does.
     */
    constexpr void step(const Vec3<T>& accel, const Vec3<T>& gyro, T dt) {
        q = q.integrate_body_rates(gyro, dt);

        const T a_mag = accel.norm();
        if (!(a_mag > T{0}) || !damp::isfinite(a_mag)) {
            return;
        }

        // Specific force in body for body→world ZYX attitude (ψ, θ, φ):
        // Rᵀ·[0,0,1] = [−sinθ, sinφ·cosθ, cosφ·cosθ]. Invert for roll/pitch;
        // the pitch denominator is the horizontal (y,z) magnitude.
        const Vec3<T> a_norm = accel / a_mag;
        const T       roll = damp::atan2(a_norm[1], a_norm[2]);
        const T       pitch = damp::atan2(-a_norm[0], damp::hypot(a_norm[1], a_norm[2]));

        // Keep the gyro-integrated yaw: the accelerometer carries no heading.
        const T             yaw = Euler<T, EulerOrder::ZYX>::from_quaternion(q).yaw();
        const Quaternion<T> q_accel = Quaternion<T>::from_euler(Euler<T, EulerOrder::ZYX>(yaw, pitch, roll));

        const T den = tau + dt;
        const T weight = (den > T{0}) ? (dt / den) : T{1}; // 1 − α; τ = dt = 0: follow the accelerometer
        q = Quaternion<T>::slerp(q, q_accel, weight).normalized();
    }

    [[nodiscard]] constexpr const Quaternion<T>& orientation() const { return q; } ///< body→world
};

namespace detail {

/**
 * Mahony-style field-alignment error (measured × estimated, body frame) for a
 * body→world @p q: up from @p a_norm (unit specific force) and heading from
 * @p m_norm, whose world-frame reference is flattened to north = +Y.
 */
template<typename T>
[[nodiscard]] constexpr Vec3<T> attitude_field_error(const Quaternion<T>& q, const Vec3<T>& a_norm, const Vec3<T>& m_norm) {
    const Vec3<T>& up_meas = a_norm;
    const Vec3<T>  up_est = q.conjugate().rotate(Vec3<T>{T{0}, T{0}, T{1}});
    const Vec3<T>  h = q.rotate(m_norm);
    const Vec3<T>  m_est = q.conjugate().rotate(Vec3<T>{T{0}, damp::hypot(h[0], h[1]), h[2]});
    return up_meas.cross(up_est) + m_norm.cross(m_est);
}

} // namespace detail

/**
 * @brief Madgwick gradient-descent AHRS filter
 *
 * Fuses gyroscope, accelerometer, and magnetometer using a gradient-descent
 * algorithm to minimize orientation error. Conventions as
 * @ref ComplementaryFilter (body → world, specific-force accel).
 *
 * @see "An efficient orientation filter for inertial and inertial/magnetic
 *      sensor arrays" (Madgwick, 2010)
 *
 * @tparam T Scalar type (default: float)
 */
template<typename T = float>
class MadgwickFilter {
    T             beta;
    Quaternion<T> q;

public:
    constexpr MadgwickFilter(T beta_gain = static_cast<T>(0.1)) : beta(beta_gain), q(Quaternion<T>::identity()) {}

    /** @brief One tick: gyro rate plus a gradient-descent correction toward the accel/mag fields, integrated over @p dt. */
    constexpr void step(const Vec3<T>& accel, const Vec3<T>& gyro, const Vec3<T>& mag, T dt) {
        const Vec3<T> e = detail::attitude_field_error(q, accel.normalized(), mag.normalized());

        // q̇_gyro = ½ q ⊗ ω
        const Quaternion<T> q_dot_gyro = static_cast<T>(0.5) * q * Quaternion<T>{T{0}, gyro[0], gyro[1], gyro[2]};

        // Gradient-descent correction: ½ q ⊗ [0, e], normalized so the step size
        // is fixed at β regardless of error magnitude (the Madgwick characteristic).
        Quaternion<T> grad = static_cast<T>(0.5) * q * Quaternion<T>{T{0}, e[0], e[1], e[2]};
        const T       grad_norm = grad.norm();
        if (grad_norm > T{0}) {
            grad = grad / grad_norm;
        }

        const Quaternion<T> q_dot = Quaternion<T>(q_dot_gyro + beta * grad);
        q = (q + q_dot * dt).normalized();
    }

    [[nodiscard]] constexpr const Quaternion<T>& orientation() const { return q; } ///< body→world
};

/**
 * @brief Mahony nonlinear complementary filter with PI correction
 *
 * Uses a proportional-integral controller on the orientation error
 * to correct gyroscope drift. Conventions as @ref ComplementaryFilter
 * (body → world, specific-force accel).
 *
 * @see "Nonlinear Complementary Filters on the Special Orthogonal Group"
 *      (Mahony et al., 2008)
 *
 * @tparam T Scalar type (default: float)
 */
template<typename T = float>
class MahonyFilter {
    T             Kp, Ki;
    Vec3<T>       integral_error{};
    Quaternion<T> q;

public:
    constexpr MahonyFilter(T Kp_gain = static_cast<T>(0.5), T Ki_gain = T{0})
        : Kp(Kp_gain), Ki(Ki_gain), q(Quaternion<T>::identity()) {}

    /** @brief One tick: gyro rate plus the PI correction toward the accel/mag fields, integrated over @p dt. */
    constexpr void step(const Vec3<T>& accel, const Vec3<T>& gyro, const Vec3<T>& mag, T dt) {
        // Without the magnetometer term yaw drift is uncorrected.
        const Vec3<T> error = detail::attitude_field_error(q, accel.normalized(), mag.normalized());
        integral_error += error * dt;
        const Vec3<T> omega = gyro + Kp * error + Ki * integral_error;
        q = q.integrate_body_rates(omega, dt);
    }

    [[nodiscard]] constexpr const Quaternion<T>& orientation() const { return q; } ///< body→world
};

namespace detail {

/** @brief The shape of the 6-state attitude transition: I, but for δθ's off-diagonal and −dt on δb_g. */
inline constexpr auto kAttitudeTransitionPattern = SparsityPattern<6, 6>::identity().values_off_diagonal(0, 0, 3).value_diagonal(0, 3, 3);

/** @brief The shape of its reset Jacobian: I, but for δθ's off-diagonal. */
inline constexpr auto kAttitudeResetPattern = SparsityPattern<6, 6>::identity().values_off_diagonal(0, 0, 3);

/** @brief The shape of an accelerometer or magnetometer reading: [a]× on δθ, zero on its diagonal. */
inline constexpr auto kAttitudeVectorPattern = SparsityPattern<3, 6>{}.values_off_diagonal(0, 0, 3);

/**
 * Gyro predict + error covariance for 6-state attitude ESKF.
 * F uses first-order body-rate kinematics: δθ⁺ ≈ (I − [ω]×dt) δθ − dt δb_g
 * (same structure as @ref ins_error_jacobian attitude block).
 */
template<typename T>
constexpr void eskf_attitude_predict(
    ErrorStateKalmanFilter<6, T>& eskf, Quaternion<T>& q_nom, Vec3<T>& b_g_nom, const Vec3<T>& gyro_meas, T dt
) {
    if (!(dt > T{0})) {
        return;
    }
    const Vec3<T> omega = gyro_meas - b_g_nom;
    q_nom = q_nom.integrate_body_rates(omega, dt);

    SparseMatrix<kAttitudeTransitionPattern, T> F{};
    F.template set_block<0, 0>(Matrix<3, 3, T>(skew_sym(omega) * (-dt)));
    F.template at<0, 3>() = -dt;
    F.template at<1, 4>() = -dt;
    F.template at<2, 5>() = -dt;
    eskf.predict(F);
}

template<typename T>
constexpr void eskf_attitude_inject(ErrorStateKalmanFilter<6, T>& eskf, Quaternion<T>& q_nom, Vec3<T>& b_g_nom) {
    const auto    dx = eskf.error_state();
    const Vec3<T> dth{dx[0], dx[1], dx[2]};
    q_nom = (q_nom * exp_delta_theta(dth)).normalized();
    b_g_nom = b_g_nom + Vec3<T>{dx[3], dx[4], dx[5]};
    SparseMatrix<kAttitudeResetPattern, T> G_reset{};
    G_reset.template set_block<0, 0>(attitude_reset_jacobian(dth));
    eskf.reset_error_state(G_reset);
}

template<typename T>
[[nodiscard]] constexpr bool eskf_accel_gate_ok(const Vec3<T>& accel, const Vec3<T>& g_nav, T gate) {
    const T g_norm = g_nav.norm();
    const T a_norm = accel.norm();
    if (!(a_norm > T{0}) || !damp::isfinite(a_norm) || !(g_norm > T{0})) {
        return false;
    }
    if (gate > T{0} && damp::abs(a_norm / g_norm - T{1}) > gate) {
        return false;
    }
    return true;
}

template<typename T>
constexpr void eskf_apply_accel(
    ErrorStateKalmanFilter<6, T>& eskf,
    Quaternion<T>&                q_nom,
    Vec3<T>&                      b_g_nom,
    const Vec3<T>&                accel_meas,
    const Vec3<T>&                g_nav,
    const Matrix<3, 3, T>&        R_a
) {
    const Vec3<T> a_pred = q_nom.conjugate().rotate(Vec3<T>{-g_nav[0], -g_nav[1], -g_nav[2]});
    ColVec<3, T>  innov{};
    innov[0] = accel_meas[0] - a_pred[0];
    innov[1] = accel_meas[1] - a_pred[1];
    innov[2] = accel_meas[2] - a_pred[2];
    SparseMatrix<kAttitudeVectorPattern, T> H{};
    H.template set_block<0, 0>(skew_sym(a_pred));
    if (!eskf.update_innovation(innov, H, R_a)) {
        return;
    }
    eskf_attitude_inject(eskf, q_nom, b_g_nom);
}

/** @brief Gyro predict + gated specific-force (IMU, NY=3). Used only by ESKFOrientationFilter. */
template<typename T>
constexpr void eskf_tick_imu(
    ErrorStateKalmanFilter<6, T>& eskf,
    Quaternion<T>&                q_nom,
    Vec3<T>&                      b_g_nom,
    const Vec3<T>&                accel_meas,
    const Vec3<T>&                gyro_meas,
    T                             dt,
    const Vec3<T>&                g_nav,
    T                             accel_gate,
    const Matrix<3, 3, T>&        R
) {
    eskf_attitude_predict(eskf, q_nom, b_g_nom, gyro_meas, dt);
    if (!eskf_accel_gate_ok(accel_meas, g_nav, accel_gate)) {
        return;
    }
    eskf_apply_accel(eskf, q_nom, b_g_nom, accel_meas, g_nav, R);
}

/** @brief Gyro predict + gated accel + independent mag (MARG, NY=6). Used only by ESKFOrientationFilter. */
template<typename T>
constexpr void eskf_tick_marg(
    ErrorStateKalmanFilter<6, T>& eskf,
    Quaternion<T>&                q_nom,
    Vec3<T>&                      b_g_nom,
    const Vec3<T>&                accel_meas,
    const Vec3<T>&                gyro_meas,
    const Vec3<T>&                mag_meas,
    T                             dt,
    const Vec3<T>&                g_nav,
    const Vec3<T>&                m_nav,
    T                             accel_gate,
    const Matrix<6, 6, T>&        R
) {
    eskf_attitude_predict(eskf, q_nom, b_g_nom, gyro_meas, dt);

    if (eskf_accel_gate_ok(accel_meas, g_nav, accel_gate)) {
        eskf_apply_accel(eskf, q_nom, b_g_nom, accel_meas, g_nav, Matrix<3, 3, T>(R.template block<3, 3>(0, 0)));
    }

    const Vec3<T> m_pred = q_nom.conjugate().rotate(m_nav);
    ColVec<3, T>  innov_m{};
    innov_m[0] = mag_meas[0] - m_pred[0];
    innov_m[1] = mag_meas[1] - m_pred[1];
    innov_m[2] = mag_meas[2] - m_pred[2];
    SparseMatrix<kAttitudeVectorPattern, T> Hm{};
    Hm.template set_block<0, 0>(skew_sym(m_pred));
    const Matrix<3, 3, T> R_m = R.template block<3, 3>(3, 3);
    if (eskf.update_innovation(innov_m, Hm, R_m)) {
        eskf_attitude_inject(eskf, q_nom, b_g_nom);
    }
}

} // namespace detail

/**
 * @brief Turnkey 6-state attitude ESKF: owns q, b_g, and the error filter
 *
 * Nominal state and error ESKF in one object:
 * @f[
 *   x = \bigl\{ q,\; b_g \bigr\},
 *   \qquad
 *   \delta x = \bigl[ \delta\theta;\; \delta b_g \bigr] \in \mathbb{R}^{6}
 * @f]
 * with body→world quaternion @f$ q @f$ and gyro bias @f$ b_g @f$. Each step():
 *
 * 1. Predict — integrate @f$ \omega = \omega_{\mathrm{meas}} - b_g @f$
 *    (right-multiplicative, @c Quaternion::integrate_body_rates); propagate
 *    @f$ P @f$ with
 *    @f$ \delta\theta^+ \approx \bigl(I - [\omega]\times\Delta t\bigr)\delta\theta
 *        - \Delta t\,\delta b_g @f$
 *    (same attitude block as @ref ins_error_jacobian).
 * 2. Accel (tilt) — specific-force residual
 *    @f$ a_b \approx -R(q)^{\top} g_n @f$ (INS convention; matches
 *    specific_force_at_rest). Jacobian @f$ H = [a_{\mathrm{pred}}]\times @f$.
 *    Skipped when the magnitude gate
 *    @f$ \bigl|\,\|a\|/\|g_n\| - 1\bigr| > \gamma @f$ fires
 *    (@ref set_accel_gate; default @f$ \gamma = 0 @f$ = always update).
 * 3. Mag (MARG only, NY=6) — body field residual
 *    @f$ m_b \approx R(q)^{\top} m_n @f$, independent of the accel gate so
 *    heading can still correct under motion.
 * 4. Inject — @f$ q \leftarrow q \otimes \exp(\delta\theta) @f$, add
 *    @f$ \delta b_g @f$ into @f$ b_g @f$, reset the error state (Solà).
 *
 * @tparam T  Scalar (default float on target)
 * @tparam NY 3 = IMU (pitch/roll; yaw free), 6 = MARG (full attitude with mag)
 *
 * @note Compare with MATLAB®-style AHRS / Madgwick IMU vs MARG; measurement model
 *       follows strapdown specific force.
 *
 * @see design::eskf_imu
 * @see design::eskf_marg
 * @see ErrorStateKalmanFilter
 * @see "Quaternion kinematics for the error-state Kalman filter" (Solà, 2017), §5–6
 *
 * Example (design-is-deploy):
 * @code
 * static constinit ESKFOrientationFilter<float, 6> filt =
 *     design::eskf_marg(0.003f, 0.03f, 0.3f, 0.0001f, 0.01f);
 * // In loop: filt.step(accel, gyro, mag, dt);
 * @endcode
 */
template<typename T = float, size_t NY = 6>
    requires(NY == 3 || NY == 6)
class ESKFOrientationFilter {
public:
    static constexpr size_t kNdx = 6;
    static constexpr size_t kNy = NY;

    /**
     * @brief Default MARG design (NY=6 only)
     *
     * Uses modest default noise densities and @c dt = 0.01. Prefer constructing
     * from design::eskf_marg when tuning for a real sensor.
     */
    constexpr ESKFOrientationFilter()
        requires(NY == 6)
        : ESKFOrientationFilter(design::eskf_marg(static_cast<T>(0.003), static_cast<T>(0.03), static_cast<T>(0.3), static_cast<T>(0.0001), static_cast<T>(0.01))) {}

    /**
     * @brief Default IMU design (NY=3 only)
     *
     * Prefer constructing from design::eskf_imu when tuning for a real sensor.
     */
    constexpr ESKFOrientationFilter()
        requires(NY == 3)
        : ESKFOrientationFilter(design::eskf_imu(static_cast<T>(0.003), static_cast<T>(0.03), static_cast<T>(0.0001), static_cast<T>(0.01))) {}

    /**
     * @brief Construct from design payload
     *
     * Non-explicit so design-is-deploy copy-init works:
     * @code
     * static constinit ESKFOrientationFilter<float, 6> filt = design::eskf_marg(...);
     * @endcode
     *
     * @param design  @ref design::ESKFResult from design::eskf_imu or design::eskf_marg
     */
    template<typename U>
    constexpr ESKFOrientationFilter(const design::ESKFResult<kNdx, kNy, U>& design) // NOLINT(google-explicit-constructor): assign-from-design
        : eskf_(design), R_(design.R.template as<T>()), q_nom_(Quaternion<T>::identity()), b_g_nom_{} {
        if (!design.success && std::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    [[nodiscard]] constexpr const Quaternion<T>&         orientation() const { return q_nom_; } ///< body→world
    [[nodiscard]] constexpr const Vec3<T>&               gyro_bias() const { return b_g_nom_; } ///< [rad/s]
    [[nodiscard]] constexpr const Matrix<kNdx, kNdx, T>& covariance() const {
        return eskf_.covariance();
    } ///< error-state P
    [[nodiscard]] constexpr const Matrix<kNdx, kNdx, T>& process_noise() const {
        return eskf_.process_noise_covariance();
    }
    [[nodiscard]] constexpr const Matrix<kNy, kNy, T>& measurement_noise() const {
        return R_;
    }
    [[nodiscard]] constexpr T accel_gate() const { return accel_gate_; }

    constexpr void set_orientation(const Quaternion<T>& q) { q_nom_ = q.normalized(); }
    constexpr void set_gyro_bias(const Vec3<T>& b) { b_g_nom_ = b; }

    /**
     * @brief Specific-force magnitude gate fraction
     *
     * Skip the accel (tilt) update when
     * @f$ \bigl|\,\|a_{\mathrm{meas}}\|/\|g_n\| - 1\bigr| > \gamma @f$.
     * Zero (default) disables the gate. Mag updates (NY=6) are never gated by this.
     *
     * @param gate  γ ≥ 0
     */
    constexpr void set_accel_gate(T gate) { accel_gate_ = gate; }

    /** @brief Replace the error covariance; the error state is kept. */
    constexpr void set_covariance(const Matrix<kNdx, kNdx, T>& P) { eskf_.reset(eskf_.error_state(), P); }
    /** @brief Replace the accel (and mag) measurement covariance. */
    constexpr void set_measurement_noise(const Matrix<kNy, kNy, T>& R) { R_ = R; }

    /**
     * @brief Re-initialize nominal orientation, gyro bias, and error covariance
     *
     * Clears the ESKF error state. Q/R and the accel gate are kept.
     */
    constexpr void reset(const Quaternion<T>& q0, const Vec3<T>& b_g0, const Matrix<kNdx, kNdx, T>& P0) {
        q_nom_ = q0.normalized();
        b_g_nom_ = b_g0;
        eskf_.reset(ColVec<kNdx, T>{}, P0);
    }

    /**
     * @brief One period — IMU (NY=3): gyro predict + gated specific-force tilt
     *
     * Pitch/roll observable; yaw free. Accel is specific force (body)
     * (at rest ENU: @f$ a \approx (0,0,+g) @f$).
     *
     * @param accel  specific force [m/s²] (body)
     * @param gyro   angular rate [rad/s] (body)
     * @param dt     sample period [s]
     * @param g_nav  gravity in nav frame (default ENU (0,0,-g))
     */
    constexpr void step(
        const Vec3<T>& accel,
        const Vec3<T>& gyro,
        T              dt,
        const Vec3<T>& g_nav = gravity_nav(NavFrame::ENU, kStandardGravity<T>)
    )
        requires(NY == 3)
    {
        detail::eskf_tick_imu(eskf_, q_nom_, b_g_nom_, accel, gyro, dt, g_nav, accel_gate_, R_);
    }

    /**
     * @brief One period — MARG (NY=6): gyro predict + gated accel + independent mag
     *
     * Accel residual is magnitude-gated (tilt). Magnetometer residual uses
     * @f$ m_b \approx R(q)^{\top} m_n @f$ and is applied even when the accel gate
     * skips, so heading remains observable under dynamics.
     *
     * @param accel  specific force [m/s²] (body)
     * @param gyro   angular rate [rad/s] (body)
     * @param mag    magnetic field in body (same units as @p m_nav)
     * @param dt     sample period [s]
     * @param g_nav  gravity in nav (default ENU (0,0,-g))
     * @param m_nav  magnetic reference in nav (default (0,1,0))
     */
    constexpr void step(
        const Vec3<T>& accel,
        const Vec3<T>& gyro,
        const Vec3<T>& mag,
        T              dt,
        const Vec3<T>& g_nav = gravity_nav(NavFrame::ENU, kStandardGravity<T>),
        const Vec3<T>& m_nav = Vec3<T>{T{0}, T{1}, T{0}}
    )
        requires(NY == 6)
    {
        detail::eskf_tick_marg(eskf_, q_nom_, b_g_nom_, accel, gyro, mag, dt, g_nav, m_nav, accel_gate_, R_);
    }

private:
    ErrorStateKalmanFilter<kNdx, T> eskf_{};
    Matrix<kNy, kNy, T>             R_{}; ///< Accel (and mag) measurement covariance

    Quaternion<T> q_nom_{Quaternion<T>::identity()};
    Vec3<T>       b_g_nom_{};
    T             accel_gate_{T{0}};
};

template<typename T>
ESKFOrientationFilter(const design::ESKFResult<6, 6, T>&) -> ESKFOrientationFilter<T, 6>;

template<typename T>
ESKFOrientationFilter(const design::ESKFResult<6, 3, T>&) -> ESKFOrientationFilter<T, 3>;

} // namespace damp
