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
 * runtime tick is ESKFOrientationFilter::update (not free functions that
 * take a bare @c ErrorStateKalmanFilter).
 *
 * @see Madgwick, "An efficient orientation filter" (2010)
 * @see Mahony et al., "Nonlinear Complementary Filters on SO(3)" (2008)
 * @see Solà et al., "Quaternion kinematics for the error-state Kalman filter" (2017)
 */

#include <cstddef>

#include "damp/math/geometry.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/core.hpp"
#include "eskf.hpp"

namespace damp {

namespace detail {

/// Negative skew-symmetric (cross-product) matrix: skew_neg(v)·u = -(v × u) = u × v.
/// Used for the small-angle attitude Jacobian -[v]× and Solà reset I − ½[δθ]×.
template<typename T>
[[nodiscard]] constexpr Matrix<3, 3, T> skew_neg(const Vec3<T>& v) {
    Matrix<3, 3, T> m = Matrix<3, 3, T>::zeros();
    m(0, 1) = v[2];
    m(0, 2) = -v[1];
    m(1, 0) = -v[2];
    m(1, 2) = v[0];
    m(2, 0) = v[1];
    m(2, 1) = -v[0];
    return m;
}

} // namespace detail

/**
 * @brief Simple complementary filter for orientation estimation
 *
 * Blends gyroscope integration (high-frequency) with accelerometer tilt
 * (low-frequency) using a tunable alpha parameter.
 *
 * @tparam T Scalar type (default: float)
 */
template<typename T = float>
class ComplementaryFilter {
    T             alpha;
    Quaternion<T> q;

public:
    constexpr ComplementaryFilter(T alpha_gain = static_cast<T>(0.98)) : alpha(alpha_gain), q(Quaternion<T>::identity()) {}

    constexpr void update(const Vec3<T>& accel, const Vec3<T>& gyro, T dt) {
        Vec3<T> a_norm = accel.normalized();

        // Body gravity for attitude q (world→body) is R(q)·[0,0,-g] = g·[…],
        // giving a_norm = [-sinθ, sinφ·cosθ, -cosφ·cosθ]. Invert for roll/pitch.
        // The pitch denominator is the horizontal (y,z) magnitude, not the full
        // vector norm (which is 1 after normalization and would collapse the tilt).
        T roll = damp::atan2(a_norm[1], -a_norm[2]);
        T pitch = damp::atan2(-a_norm[0], damp::hypot(a_norm[1], a_norm[2]));

        Euler<T, EulerOrder::ZYX> euler(T{0}, pitch, roll);

        Quaternion<T> q_accel = Quaternion<T>::from_euler(euler);

        q = q.integrate_body_rates(gyro, dt);
        q = Quaternion<T>::slerp(q, q_accel, T{1} - alpha);
        q = q.normalized();
    }

    [[nodiscard]] constexpr const Quaternion<T>& orientation() const { return q; }
};

/**
 * @brief Madgwick gradient-descent AHRS filter
 *
 * Fuses gyroscope, accelerometer, and magnetometer using a gradient-descent
 * algorithm to minimize orientation error.
 *
 * @see Madgwick, "An efficient orientation filter for inertial and
 *      inertial/magnetic sensor arrays" (2010)
 *
 * @tparam T Scalar type (default: float)
 */
template<typename T = float>
class MadgwickFilter {
    T             beta;
    Quaternion<T> q;

public:
    constexpr MadgwickFilter(T beta_gain = static_cast<T>(0.1)) : beta(beta_gain), q(Quaternion<T>::identity()) {}

    constexpr void update(const Vec3<T>& accel, const Vec3<T>& gyro, const Vec3<T>& mag, T dt) {
        Vec3<T> a_norm = accel.normalized();
        Vec3<T> m_norm = mag.normalized();

        // Field-alignment error in the body frame (gravity + magnetometer), the
        // gradient the filter descends. Gravity: measured field crossed with the
        // estimated world-up axis in body, R(q)·[0,0,1] (antiparallel at the fixed
        // point). Magnetometer: refer the field to the world frame, flatten to the
        // horizontal reference (north = +Y), rotate the estimate back to body, and
        // cross (estimate × meas) — opposite order to gravity since the mag
        // estimate is parallel to the measurement at the fixed point.
        const Vec3<T> v = q.rotate(Vec3<T>{T{0}, T{0}, T{1}});
        Vec3<T>       e = a_norm.cross(v);

        const Vec3<T> h = q.conjugate().rotate(m_norm);
        const T       b_horizontal = damp::hypot(h[0], h[1]);
        const Vec3<T> w = q.rotate(Vec3<T>{T{0}, b_horizontal, h[2]});
        e += w.cross(m_norm);

        // q̇_gyro = ½ q ⊗ ω
        const Quaternion<T> q_dot_gyro = static_cast<T>(0.5) * q * Quaternion<T>{T{0}, gyro[0], gyro[1], gyro[2]};

        // Gradient-descent correction: ½ q ⊗ [0, e], normalized so the step size
        // is fixed at β regardless of error magnitude (the Madgwick characteristic).
        // +e drives the error toward zero, so the term is added to the gyro rate.
        Quaternion<T> grad = static_cast<T>(0.5) * q * Quaternion<T>{T{0}, e[0], e[1], e[2]};
        const T       grad_norm = grad.norm();
        if (grad_norm > static_cast<T>(1e-9)) {
            grad = grad / grad_norm;
        }

        const Quaternion<T> q_dot = Quaternion<T>(q_dot_gyro + beta * grad);
        q = q + q_dot * dt;
        q = q.normalized();
    }

    [[nodiscard]] constexpr const Quaternion<T>& orientation() const { return q; }
};

/**
 * @brief Mahony nonlinear complementary filter with PI correction
 *
 * Uses a proportional-integral controller on the orientation error
 * to correct gyroscope drift.
 *
 * @see Mahony et al., "Nonlinear Complementary Filters on the Special
 *      Orthogonal Group" (2008)
 *
 * @tparam T Scalar type (default: float)
 */
template<typename T = float>
class MahonyFilter {
    T             Kp, Ki;
    Vec3<T>       integral_error;
    Quaternion<T> q;

public:
    constexpr MahonyFilter(T Kp_gain = static_cast<T>(0.5), T Ki_gain = T{0})
        : Kp(Kp_gain), Ki(Ki_gain), q(Quaternion<T>::identity()) {}

    constexpr void update(const Vec3<T>& accel, const Vec3<T>& gyro, const Vec3<T>& mag, T dt) {
        Vec3<T> a_norm = accel.normalized();
        Vec3<T> m_norm = mag.normalized();

        // Gravity (accelerometer) correction: measured field crossed with the
        // estimated world-up axis in body, R(q)·[0,0,1]. (At the fixed point
        // a_norm is antiparallel to this axis, which is the stable equilibrium.)
        Vec3<T> v = q.rotate(Vec3<T>{T{0}, T{0}, T{1}});
        Vec3<T> error = a_norm.cross(v);

        // Magnetometer correction (yaw): refer the measured field to the world
        // frame, flatten it to the horizontal reference direction (north = +Y,
        // matching this file's world frame), rotate that estimate back to body,
        // and cross it with the measurement. The estimate is parallel to the
        // measurement at the fixed point, so the cross order is (estimate × meas)
        // — opposite to gravity — to keep the same (stable) sign. Without this
        // term yaw drift is uncorrected.
        const Vec3<T> h = q.conjugate().rotate(m_norm);
        const T       b_horizontal = damp::hypot(h[0], h[1]);
        const Vec3<T> w = q.rotate(Vec3<T>{T{0}, b_horizontal, h[2]});
        error += w.cross(m_norm);

        integral_error += error * dt;
        Vec3<T> feedback = Kp * error + Ki * integral_error;

        Vec3<T> omega = gyro + feedback;
        q = q.integrate_body_rates(omega, dt);
        q = q.normalized();
    }

    [[nodiscard]] constexpr const Quaternion<T>& orientation() const { return q; }
};

namespace detail {

/**
 * Gyro predict + error covariance for 6-state attitude ESKF.
 * F uses first-order body-rate kinematics: δθ⁺ ≈ (I − [ω]×dt) δθ − dt δb_g
 * (same structure as @ref ins_error_jacobian attitude block).
 */
template<size_t NY, typename T>
constexpr void eskf_attitude_predict(
    ErrorStateKalmanFilter<6, NY, T>& eskf, Quaternion<T>& q_nom, Vec3<T>& b_g_nom, const Vec3<T>& gyro_meas, T dt
) {
    if (!(dt > T{0})) {
        return;
    }
    const Vec3<T> omega = gyro_meas - b_g_nom;
    q_nom = q_nom.integrate_body_rates(omega, dt);

    // [ω]× = −skew_neg(ω);  I − [ω]×dt = I + skew_neg(ω)·dt
    Matrix<6, 6, T> F = Matrix<6, 6, T>::identity();
    F.template block<3, 3>(0, 0) = Matrix<3, 3, T>::identity() + (skew_neg(omega) * dt);
    F(0, 3) = -dt;
    F(1, 4) = -dt;
    F(2, 5) = -dt;
    const Matrix<6, 6, T> G = Matrix<6, 6, T>::identity();
    eskf.predict(
        [F, G](T) -> ErrorStateJacobian<T, 6> {
        return {F, G};
    },
        dt
    );
}

template<size_t NY, typename T>
constexpr void eskf_attitude_inject(ErrorStateKalmanFilter<6, NY, T>& eskf, Quaternion<T>& q_nom, Vec3<T>& b_g_nom) {
    const auto    dx = eskf.error_state();
    const Vec3<T> dth{dx[0], dx[1], dx[2]};
    const T       n = dth.norm();
    if (n > static_cast<T>(1e-12)) {
        const auto dq = Quaternion<T>::from_axis_angle(dth, n).value_or(Quaternion<T>::identity());
        q_nom = (q_nom * dq).normalized();
    }
    b_g_nom = b_g_nom + Vec3<T>{dx[3], dx[4], dx[5]};
    // Solà right-multiplicative reset: G_θθ ≈ I − ½[δθ]×  (arXiv:1711.02508 §5.4)
    Matrix<6, 6, T> G_reset = Matrix<6, 6, T>::identity();
    if (n > static_cast<T>(1e-12)) {
        G_reset.template block<3, 3>(0, 0) = Matrix<3, 3, T>::identity() + (skew_neg(dth) * static_cast<T>(0.5));
    }
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

template<size_t NY, typename T>
constexpr void eskf_apply_accel(
    ErrorStateKalmanFilter<6, NY, T>& eskf,
    Quaternion<T>&                    q_nom,
    Vec3<T>&                          b_g_nom,
    const Vec3<T>&                    accel_meas,
    const Vec3<T>&                    g_nav
) {
    const Vec3<T> a_pred = q_nom.conjugate().rotate(Vec3<T>{-g_nav[0], -g_nav[1], -g_nav[2]});
    ColVec<3, T>  innov{};
    innov[0] = accel_meas[0] - a_pred[0];
    innov[1] = accel_meas[1] - a_pred[1];
    innov[2] = accel_meas[2] - a_pred[2];
    Matrix<3, 6, T> H = Matrix<3, 6, T>::zeros();
    H.template block<3, 3>(0, 0) = skew_neg(a_pred) * T{-1};
    Matrix<3, 3, T> R_a = Matrix<3, 3, T>::zeros();
    for (size_t i = 0; i < 3; ++i) {
        R_a(i, i) = eskf.measurement_noise_covariance()(i, i);
    }
    if (!eskf.update_innovation(innov, H, R_a)) {
        return;
    }
    eskf_attitude_inject(eskf, q_nom, b_g_nom);
}

/// Gyro predict + gated specific-force (IMU, NY=3). Used only by ESKFOrientationFilter.
template<typename T>
constexpr void eskf_tick_imu(
    ErrorStateKalmanFilter<6, 3, T>& eskf,
    Quaternion<T>&                   q_nom,
    Vec3<T>&                         b_g_nom,
    const Vec3<T>&                   accel_meas,
    const Vec3<T>&                   gyro_meas,
    T                                dt,
    const Vec3<T>&                   g_nav,
    T                                accel_gate
) {
    eskf_attitude_predict(eskf, q_nom, b_g_nom, gyro_meas, dt);
    if (!eskf_accel_gate_ok(accel_meas, g_nav, accel_gate)) {
        return;
    }
    eskf_apply_accel(eskf, q_nom, b_g_nom, accel_meas, g_nav);
}

/// Gyro predict + gated accel + independent mag (MARG, NY=6). Used only by ESKFOrientationFilter.
template<typename T>
constexpr void eskf_tick_marg(
    ErrorStateKalmanFilter<6, 6, T>& eskf,
    Quaternion<T>&                   q_nom,
    Vec3<T>&                         b_g_nom,
    const Vec3<T>&                   accel_meas,
    const Vec3<T>&                   gyro_meas,
    const Vec3<T>&                   mag_meas,
    T                                dt,
    const Vec3<T>&                   g_nav,
    const Vec3<T>&                   m_nav,
    T                                accel_gate
) {
    eskf_attitude_predict(eskf, q_nom, b_g_nom, gyro_meas, dt);

    if (eskf_accel_gate_ok(accel_meas, g_nav, accel_gate)) {
        eskf_apply_accel(eskf, q_nom, b_g_nom, accel_meas, g_nav);
    }

    const Vec3<T> m_pred = q_nom.conjugate().rotate(m_nav);
    ColVec<3, T>  innov_m{};
    innov_m[0] = mag_meas[0] - m_pred[0];
    innov_m[1] = mag_meas[1] - m_pred[1];
    innov_m[2] = mag_meas[2] - m_pred[2];
    Matrix<3, 6, T> Hm = Matrix<3, 6, T>::zeros();
    Hm.template block<3, 3>(0, 0) = skew_neg(m_pred) * T{-1};
    Matrix<3, 3, T> R_m = Matrix<3, 3, T>::zeros();
    for (size_t i = 0; i < 3; ++i) {
        R_m(i, i) = eskf.measurement_noise_covariance()(i + 3, i + 3);
    }
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
 * with body→world quaternion @f$ q @f$ and gyro bias @f$ b_g @f$. Each update:
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
 *       follows strapdown specific force, not “gravity direction = accel”.
 *
 * @see design::eskf_imu
 * @see design::eskf_marg
 * @see ErrorStateKalmanFilter
 * @see Solà et al., "Quaternion kinematics for the error-state Kalman filter" (2017), §5–6
 *
 * Example (design-is-deploy):
 * @code
 * static constinit ESKFOrientationFilter<float, 6> filt =
 *     design::eskf_marg(0.003f, 0.03f, 0.3f, 0.0001f, 0.01f);
 * // In loop: filt.update(accel, gyro, mag, dt);
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
        : eskf_(design::eskf_marg(static_cast<T>(0.003), static_cast<T>(0.03), static_cast<T>(0.3), static_cast<T>(0.0001), static_cast<T>(0.01))),
          q_nom_(Quaternion<T>::identity()),
          b_g_nom_{} {}

    /**
     * @brief Default IMU design (NY=3 only)
     *
     * Prefer constructing from design::eskf_imu when tuning for a real sensor.
     */
    constexpr ESKFOrientationFilter()
        requires(NY == 3)
        : eskf_(design::eskf_imu(static_cast<T>(0.003), static_cast<T>(0.03), static_cast<T>(0.0001), static_cast<T>(0.01))),
          q_nom_(Quaternion<T>::identity()),
          b_g_nom_{} {}

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
    constexpr ESKFOrientationFilter(const design::ESKFResult<kNdx, kNy, T>& design)
        : eskf_(design), q_nom_(Quaternion<T>::identity()), b_g_nom_{} {}

    [[nodiscard]] constexpr const Quaternion<T>&         orientation() const { return q_nom_; } ///< body→world
    [[nodiscard]] constexpr const Vec3<T>&               gyro_bias() const { return b_g_nom_; } ///< [rad/s]
    [[nodiscard]] constexpr const Matrix<kNdx, kNdx, T>& covariance() const {
        return eskf_.covariance();
    } ///< error-state P

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

    constexpr void set_covariance(const Matrix<kNdx, kNdx, T>& P) { eskf_.set_covariance(P); }
    constexpr void set_process_noise(const Matrix<kNdx, kNdx, T>& Q) {
        eskf_.set_process_noise_covariance(Q);
    }
    constexpr void set_measurement_noise(const Matrix<kNy, kNy, T>& R) {
        eskf_.set_measurement_noise_covariance(R);
    }

    /**
     * @brief Re-initialize nominal orientation, gyro bias, and error covariance
     *
     * Clears the ESKF error state. Q/R and the accel gate are kept.
     */
    constexpr void reset(
        const Quaternion<T>&         q0 = Quaternion<T>::identity(),
        const Vec3<T>&               b_g0 = Vec3<T>{},
        const Matrix<kNdx, kNdx, T>& P0 = Matrix<kNdx, kNdx, T>::identity()
    ) {
        q_nom_ = q0.normalized();
        b_g_nom_ = b_g0;
        eskf_.set_error_state(ColVec<kNdx, T>{});
        eskf_.set_covariance(P0);
    }

    /**
     * @brief One period — IMU (NY=3): gyro predict + gated specific-force tilt
     *
     * Pitch/roll observable; yaw free. Accel must be specific force (body),
     * not “gravity direction” (at rest ENU: @f$ a \approx (0,0,+g) @f$).
     *
     * @param accel  specific force [m/s²] (body)
     * @param gyro   angular rate [rad/s] (body)
     * @param dt     sample period [s]
     * @param g_nav  gravity in nav frame (default ENU (0,0,-g))
     */
    constexpr void update(
        const Vec3<T>& accel,
        const Vec3<T>& gyro,
        T              dt,
        const Vec3<T>& g_nav = Vec3<T>{T{0}, T{0}, static_cast<T>(-9.81)}
    )
        requires(NY == 3)
    {
        detail::eskf_tick_imu(eskf_, q_nom_, b_g_nom_, accel, gyro, dt, g_nav, accel_gate_);
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
    constexpr void update(
        const Vec3<T>& accel,
        const Vec3<T>& gyro,
        const Vec3<T>& mag,
        T              dt,
        const Vec3<T>& g_nav = Vec3<T>{T{0}, T{0}, static_cast<T>(-9.81)},
        const Vec3<T>& m_nav = Vec3<T>{T{0}, T{1}, T{0}}
    )
        requires(NY == 6)
    {
        detail::eskf_tick_marg(eskf_, q_nom_, b_g_nom_, accel, gyro, mag, dt, g_nav, m_nav, accel_gate_);
    }

private:
    ErrorStateKalmanFilter<kNdx, kNy, T> eskf_{};

    Quaternion<T> q_nom_{Quaternion<T>::identity()};
    Vec3<T>       b_g_nom_{};
    T             accel_gate_{T{0}};
};

template<typename T>
ESKFOrientationFilter(const design::ESKFResult<6, 6, T>&) -> ESKFOrientationFilter<T, 6>;

template<typename T>
ESKFOrientationFilter(const design::ESKFResult<6, 3, T>&) -> ESKFOrientationFilter<T, 3>;

} // namespace damp
