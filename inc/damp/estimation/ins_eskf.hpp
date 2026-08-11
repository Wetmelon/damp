// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file damp/estimation/ins_eskf.hpp
 * @brief 15-state navigation ESKF for loosely-coupled strapdown INS aiding
 *
 * Nominal nav state lives in @ref InsState and is integrated by
 * @c mechanize_step. This header adds error dynamics, injection, sparse
 * absolute updates, and the Tier-3 @ref InsNavigator runtime.
 *
 * ## Error state (15)
 *
 * @f[
 *   \delta x =
 *   [\delta\theta^\top,\; \delta v^\top,\; \delta p^\top,\;
 *    \delta b_g^\top,\; \delta b_a^\top]^\top
 * @f]
 *
 * Right-multiplicative attitude error (body frame):
 * @f$ q_{\mathrm{true}} = q \otimes \exp(\delta\theta) @f$.
 * Flat Earth (no transport rate / Earth rate).
 *
 * Continuous error model (Solà local ESKF;
 * @f$ \omega = \omega_m - b_g @f$, @f$ a_b = a_m - b_a @f$,
 * @f$ R = R(q) @f$ body→nav):
 *
 * @f[
 *   \dot{\delta\theta} = -[\omega]_\times \delta\theta - \delta b_g,\quad
 *   \dot{\delta v} = -R[a_b]_\times \delta\theta - R\,\delta b_a,\quad
 *   \dot{\delta p} = \delta v,\quad
 *   \dot{\delta b}_g = 0,\quad
 *   \dot{\delta b}_a = 0.
 * @f]
 *
 * Discretization is first-order: @f$ F_d = I + F_c\,\Delta t @f$.
 *
 * ## Aiding API (measurement-shaped)
 *
 * | Function | Measurement | Notes |
 * | -------- | ----------- | ----- |
 * | ins_update_position | p (3) | optional body lever-arm |
 * | ins_update_velocity / ins_update_zupt | v (3) | |
 * | ins_update_heading | yaw ψ (1) | dual-antenna / mag yaw |
 * | ins_update_orientation | q (3 as rotvec) | mocap attitude |
 * | ins_update_pose | p then q | sequential sparse updates |
 *
 * Prefer small selector H rows over a dense 15-row zero-padded Jacobian.
 *
 * @see Solà, "Quaternion kinematics for the error-state Kalman filter," arXiv:1711.02508
 * @see Groves, Principles of GNSS, Inertial, and Multisensor Integrated Navigation, 2nd ed.
 * @see eskf.hpp for the 6-state attitude ESKF core
 *
 * Example (design → deploy):
 * @code
 * constexpr auto d = design::ins_eskf_design(1e-3, 1e-2, 1e-4, 1e-3, 0.01).as<float>();
 * static constinit InsNavigator<float> nav(d, NavFrame::ENU);
 * // ISR: nav.predict(imu, dt); nav.update_position(p_gps); nav.update_heading(psi);
 * @endcode
 */

#include <cstddef>

#include "damp/math/geometry.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "eskf.hpp"
#include "ins_mechanization.hpp"

namespace damp {

/// Error-state dimension for full navigation ESKF (attitude, v, p, gyro bias, accel bias).
inline constexpr std::size_t kInsErrorDim = 15;

/// Sparse aiding measurement size (position or velocity triad).
inline constexpr std::size_t kInsMeasDim = 3;

/**
 * @brief Start indices of each 3-vector block in the 15-state error δx
 */
namespace ins_err {
inline constexpr std::size_t dtheta = 0; ///< attitude error δθ [rad] (3)
inline constexpr std::size_t dv = 3;     ///< velocity error δv [m/s] (3)
inline constexpr std::size_t dp = 6;     ///< position error δp [m] (3)
inline constexpr std::size_t dbg = 9;    ///< gyro bias error δb_g [rad/s] (3)
inline constexpr std::size_t dba = 12;   ///< accel bias error δb_a [m/s²] (3)
} // namespace ins_err

namespace detail {

/// Skew-symmetric [v]× so [v]× w = v × w.
template<typename T>
[[nodiscard]] constexpr Matrix<3, 3, T> skew_sym(const Vec3<T>& v) {
    Matrix<3, 3, T> m = Matrix<3, 3, T>::zeros();
    m(0, 1) = -v[2];
    m(0, 2) = v[1];
    m(1, 0) = v[2];
    m(1, 2) = -v[0];
    m(2, 0) = -v[1];
    m(2, 1) = v[0];
    return m;
}

/// Small-angle right-multiplicative quaternion @f$ \exp(\delta\theta) @f$.
template<typename T>
[[nodiscard]] constexpr Quaternion<T> exp_delta_theta(const Vec3<T>& dth, T eps = static_cast<T>(1e-12)) {
    const T n2 = dth[0] * dth[0] + dth[1] * dth[1] + dth[2] * dth[2];
    if (n2 <= eps * eps) {
        // First-order: [1, δθ/2]
        return Quaternion<T>{T{1}, static_cast<T>(0.5) * dth[0], static_cast<T>(0.5) * dth[1], static_cast<T>(0.5) * dth[2]}.normalized();
    }
    const T    n = damp::sqrt(n2);
    const auto q = Quaternion<T>::from_axis_angle(dth, n, eps);
    return q.value_or(Quaternion<T>::identity());
}

/// Fill a 3×15 measurement Jacobian with @f$ I_3 @f$ at column @p col0.
template<typename T>
[[nodiscard]] constexpr Matrix<kInsMeasDim, kInsErrorDim, T> selector_H(std::size_t col0) {
    Matrix<kInsMeasDim, kInsErrorDim, T> H = Matrix<kInsMeasDim, kInsErrorDim, T>::zeros();
    H.template block<3, 3>(0, col0) = Matrix<3, 3, T>::identity();
    return H;
}

/**
 * Horizontal heading from a nav-frame direction.
 * ENU (E,N,U): ψ = atan2(E, N). NED (N,E,D): ψ = atan2(E, N).
 */
template<typename T>
[[nodiscard]] constexpr T heading_from_forward_nav(const Vec3<T>& f, NavFrame frame) {
    if (frame == NavFrame::NED) {
        return damp::atan2(f[1], f[0]); // E, N
    }
    return damp::atan2(f[0], f[1]); // E, N with (E,N,U)
}

/// ∂ψ/∂f (1×3) for the heading convention of @p frame (horizontal only).
template<typename T>
[[nodiscard]] constexpr Matrix<1, 3, T> dheading_dforward(const Vec3<T>& f, NavFrame frame) {
    Matrix<1, 3, T> J = Matrix<1, 3, T>::zeros();
    if (frame == NavFrame::NED) {
        // ψ = atan2(E, N) = atan2(f1, f0)
        const T n = f[0];
        const T e = f[1];
        const T r2 = n * n + e * e;
        if (r2 > static_cast<T>(1e-24)) {
            const T inv = T{1} / r2;
            J(0, 0) = -e * inv;
            J(0, 1) = n * inv;
        }
    } else {
        // ψ = atan2(E, N) = atan2(f0, f1)
        const T e = f[0];
        const T n = f[1];
        const T r2 = e * e + n * n;
        if (r2 > static_cast<T>(1e-24)) {
            const T inv = T{1} / r2;
            J(0, 0) = n * inv;
            J(0, 1) = -e * inv;
        }
    }
    return J;
}

} // namespace detail

namespace design {

/**
 * @brief Design result for the 15-state INS error-state filter
 *
 * @tparam T Scalar type
 */
template<typename T = double>
struct InsEskfResult {
    Matrix<kInsErrorDim, kInsErrorDim, T> Q{};  ///< Process noise (discrete)
    Matrix<kInsMeasDim, kInsMeasDim, T>   R{};  ///< Default position measurement noise
    Matrix<kInsErrorDim, kInsErrorDim, T> P0{}; ///< Initial error covariance
    bool                                  success{false};

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return InsEskfResult<U>{
            Q.template as<U>(),
            R.template as<U>(),
            P0.template as<U>(),
            success,
        };
    }
};

/**
 * @brief Build Q / R / P0 for a 15-state INS ESKF from IMU noise densities
 *
 * Discrete process noise (continuous–discrete, first-order integrated RW):
 * attitude/gyro-bias as in @ref design::detail::fill_attitude_bias_q; velocity/
 * accel-bias analogous; position gets @f$ Q_{pp}=\sigma_a^2\Delta t^3/3 @f$,
 * @f$ Q_{pv}=\sigma_a^2\Delta t^2/2 @f$ from double integration of white accel.
 * Default @p R is isotropic position std (m), not a density.
 *
 * @param gyro_noise_density   [rad/s/√Hz]
 * @param accel_noise_density  [m/s²/√Hz]
 * @param gyro_bias_rw         [rad/s^{3/2}] gyro bias random walk
 * @param accel_bias_rw        [m/s²^{3/2}] accel bias random walk
 * @param dt                   [s] sample period used to discretize Q
 * @param pos_noise_std        [m] 1-σ position measurement (default R)
 * @param init_att_std         [rad] initial attitude 1-σ
 * @param init_vel_std         [m/s]
 * @param init_pos_std         [m]
 * @param init_bg_std          [rad/s]
 * @param init_ba_std          [m/s²]
 */
template<typename T = double>
[[nodiscard]] constexpr InsEskfResult<T> ins_eskf_design(
    T gyro_noise_density,
    T accel_noise_density,
    T gyro_bias_rw,
    T accel_bias_rw,
    T dt,
    T pos_noise_std = T{1},
    T init_att_std = static_cast<T>(0.1),
    T init_vel_std = T{1},
    T init_pos_std = T{10},
    T init_bg_std = static_cast<T>(0.01),
    T init_ba_std = static_cast<T>(0.1)
) {
    InsEskfResult<T> r{};
    r.Q = Matrix<kInsErrorDim, kInsErrorDim, T>::zeros();
    r.P0 = Matrix<kInsErrorDim, kInsErrorDim, T>::zeros();
    r.R = Matrix<kInsMeasDim, kInsMeasDim, T>::zeros();

    r.success = detail::nonneg_finite(dt) && (dt > T{0}) && detail::nonneg_finite(gyro_noise_density)
             && detail::nonneg_finite(accel_noise_density) && detail::nonneg_finite(gyro_bias_rw)
             && detail::nonneg_finite(accel_bias_rw) && detail::nonneg_finite(pos_noise_std)
             && (pos_noise_std > T{0}) && detail::nonneg_finite(init_att_std)
             && detail::nonneg_finite(init_vel_std) && detail::nonneg_finite(init_pos_std)
             && detail::nonneg_finite(init_bg_std) && detail::nonneg_finite(init_ba_std);
    if (!r.success) {
        return r;
    }

    const T sg2 = gyro_noise_density * gyro_noise_density;
    const T sa2 = accel_noise_density * accel_noise_density;
    const T sbg2 = gyro_bias_rw * gyro_bias_rw;
    const T sba2 = accel_bias_rw * accel_bias_rw;
    const T dt2 = dt * dt;
    const T dt3 = dt2 * dt;

    const T q_tt = sg2 * dt + sbg2 * (dt3 / T{3});
    const T q_tb = -sbg2 * (dt2 / T{2});
    const T q_bbg = sbg2 * dt;
    const T q_vv = sa2 * dt + sba2 * (dt3 / T{3});
    const T q_vba = -sba2 * (dt2 / T{2});
    const T q_bba = sba2 * dt;
    const T q_pp = sa2 * (dt3 / T{3});
    const T q_pv = sa2 * (dt2 / T{2});

    for (std::size_t i = 0; i < 3; ++i) {
        r.Q(ins_err::dtheta + i, ins_err::dtheta + i) = q_tt;
        r.Q(ins_err::dbg + i, ins_err::dbg + i) = q_bbg;
        r.Q(ins_err::dtheta + i, ins_err::dbg + i) = q_tb;
        r.Q(ins_err::dbg + i, ins_err::dtheta + i) = q_tb;

        r.Q(ins_err::dv + i, ins_err::dv + i) = q_vv;
        r.Q(ins_err::dba + i, ins_err::dba + i) = q_bba;
        r.Q(ins_err::dv + i, ins_err::dba + i) = q_vba;
        r.Q(ins_err::dba + i, ins_err::dv + i) = q_vba;

        r.Q(ins_err::dp + i, ins_err::dp + i) = q_pp;
        r.Q(ins_err::dp + i, ins_err::dv + i) = q_pv;
        r.Q(ins_err::dv + i, ins_err::dp + i) = q_pv;

        r.P0(ins_err::dtheta + i, ins_err::dtheta + i) = init_att_std * init_att_std;
        r.P0(ins_err::dv + i, ins_err::dv + i) = init_vel_std * init_vel_std;
        r.P0(ins_err::dp + i, ins_err::dp + i) = init_pos_std * init_pos_std;
        r.P0(ins_err::dbg + i, ins_err::dbg + i) = init_bg_std * init_bg_std;
        r.P0(ins_err::dba + i, ins_err::dba + i) = init_ba_std * init_ba_std;

        r.R(i, i) = pos_noise_std * pos_noise_std;
    }

    return r;
}

} // namespace design

/**
 * @brief Discrete error-state transition and noise input for one IMU step
 *
 * First-order @f$ F_d = I + F_c\,\Delta t @f$ (flat Earth). Valid when
 * @f$ \|\omega\|\Delta t \ll 1 @f$ (typical IMU rates). Attitude block is
 * @f$ I-[\omega]_\times\Delta t @f$; for large steps prefer smaller @p dt.
 * G is identity (process noise already discrete-additive in @ref design::ins_eskf_design).
 *
 * @param q     nominal body→nav quaternion (before or after step; consistent use is fine for small dt)
 * @param omega bias-corrected body rate [rad/s]
 * @param a_b   bias-corrected body specific force [m/s²]
 * @param dt    [s]
 * @return F (15×15) and G (identity 15×15 for additive Q)
 */
template<typename T = double>
[[nodiscard]] constexpr ErrorStateJacobian<T, kInsErrorDim> ins_error_jacobian(
    const Quaternion<T>& q,
    const Vec3<T>&       omega,
    const Vec3<T>&       a_b,
    T                    dt
) {
    ErrorStateJacobian<T, kInsErrorDim> ej{};
    ej.F = Matrix<kInsErrorDim, kInsErrorDim, T>::identity();
    ej.G = Matrix<kInsErrorDim, kInsErrorDim, T>::identity();

    if (!(dt > T{0})) {
        return ej;
    }

    const Mat3<T> R = q.to_dcm().matrix();
    const Mat3<T> W = detail::skew_sym(omega);
    const Mat3<T> A = detail::skew_sym(a_b);
    const Mat3<T> R_A = R * A;

    // δθ⁺ ≈ (I − [ω]× dt) δθ − dt δb_g
    ej.F.template block<3, 3>(ins_err::dtheta, ins_err::dtheta) = Matrix<3, 3, T>::identity() - (W * dt);
    ej.F.template block<3, 3>(ins_err::dtheta, ins_err::dbg) = Matrix<3, 3, T>::identity() * (-dt);

    // δv⁺ ≈ δv − R[a_b]× dt δθ − R dt δb_a
    ej.F.template block<3, 3>(ins_err::dv, ins_err::dtheta) = R_A * (-dt);
    ej.F.template block<3, 3>(ins_err::dv, ins_err::dba) = R * (-dt);

    // δp⁺ ≈ δp + dt δv  (+ ½ a dt² is only on nominal; keep F first-order)
    ej.F.template block<3, 3>(ins_err::dp, ins_err::dv) = Matrix<3, 3, T>::identity() * dt;

    return ej;
}

/**
 * @brief Sparse y = F x for the 15-state INS first-order structure (G = I path)
 *
 * @p x may be any 15×1 @c MatrixLike (e.g. @c ColVec, @c ColView) or 1×15
 * row-shaped view; elements are read in storage order of the 15-vector.
 */
template<typename T, MatrixLike X>
    requires(
        (X::rows() == kInsErrorDim && X::cols() == 1)
        || (X::rows() == 1 && X::cols() == kInsErrorDim)
    )
[[nodiscard]] constexpr ColVec<kInsErrorDim, T>
ins_F_mul(const ErrorStateJacobian<T, kInsErrorDim>& ej, const X& x) {
    const auto x_at = [&x](std::size_t i) -> T {
        if constexpr (X::cols() == 1) {
            return static_cast<T>(x(i, 0));
        } else {
            return static_cast<T>(x(0, i));
        }
    };

    ColVec<kInsErrorDim, T> y{};
    for (std::size_t i = 0; i < kInsErrorDim; ++i) {
        y[i] = x_at(i);
    }

    // θ' = Fθθ θ + Fθbg bg
    for (std::size_t i = 0; i < 3; ++i) {
        T s = T{0};
        for (std::size_t j = 0; j < 3; ++j) {
            s += ej.F(ins_err::dtheta + i, ins_err::dtheta + j) * x_at(ins_err::dtheta + j);
            s += ej.F(ins_err::dtheta + i, ins_err::dbg + j) * x_at(ins_err::dbg + j);
        }
        y[ins_err::dtheta + i] = s;
    }
    // v' = v + Fvθ θ + Fvba ba
    for (std::size_t i = 0; i < 3; ++i) {
        T s = x_at(ins_err::dv + i);
        for (std::size_t j = 0; j < 3; ++j) {
            s += ej.F(ins_err::dv + i, ins_err::dtheta + j) * x_at(ins_err::dtheta + j);
            s += ej.F(ins_err::dv + i, ins_err::dba + j) * x_at(ins_err::dba + j);
        }
        y[ins_err::dv + i] = s;
    }
    // p' = p + Fp v
    for (std::size_t i = 0; i < 3; ++i) {
        T s = x_at(ins_err::dp + i);
        for (std::size_t j = 0; j < 3; ++j) {
            s += ej.F(ins_err::dp + i, ins_err::dv + j) * x_at(ins_err::dv + j);
        }
        y[ins_err::dp + i] = s;
    }
    // bg, ba unchanged (I) — already copied into y
    return y;
}

/**
 * @brief P ← F P Fᵀ + Q using sparse F (exact for @ref ins_error_jacobian structure)
 *
 * Applies @ref ins_F_mul to each column of @p P and each row of the intermediate
 * via @c .col / @c .row views (no manual pack/unpack of dense temporaries for I/O).
 */
template<typename T = double>
constexpr void ins_propagate_covariance(
    Matrix<kInsErrorDim, kInsErrorDim, T>&       P,
    const ErrorStateJacobian<T, kInsErrorDim>&   ej,
    const Matrix<kInsErrorDim, kInsErrorDim, T>& Q
) {
    // W = F P  (each column: W.col(j) = F * P.col(j))
    Matrix<kInsErrorDim, kInsErrorDim, T> W{};
    for (std::size_t j = 0; j < kInsErrorDim; ++j) {
        W.col(j) = ins_F_mul(ej, P.col(j));
    }
    // P ← W Fᵀ + Q.  (W Fᵀ)_i,: = (F * W.row(i)ᵀ)ᵀ → assign F*row into row i
    for (std::size_t i = 0; i < kInsErrorDim; ++i) {
        P.row(i) = ins_F_mul(ej, W.row(i));
    }
    for (std::size_t i = 0; i < kInsErrorDim; ++i) {
        for (std::size_t j = 0; j < kInsErrorDim; ++j) {
            P(i, j) += Q(i, j);
        }
    }
}

/**
 * @brief Inject δx into the nominal INS state (right-multiplicative q)
 *
 * @param x   nominal state (modified)
 * @param dx  error state from the filter
 * @return reset Jacobian G for @ref ErrorStateKalmanFilter::reset_error_state
 */
template<typename T = double>
[[nodiscard]] constexpr Matrix<kInsErrorDim, kInsErrorDim, T>
ins_inject(InsState<T>& x, const ColVec<kInsErrorDim, T>& dx) {
    const Vec3<T> dth{dx[0], dx[1], dx[2]};
    const Vec3<T> dv{dx[3], dx[4], dx[5]};
    const Vec3<T> dp{dx[6], dx[7], dx[8]};
    const Vec3<T> dbg{dx[9], dx[10], dx[11]};
    const Vec3<T> dba{dx[12], dx[13], dx[14]};

    x.q = (x.q * detail::exp_delta_theta(dth)).normalized();
    x.v = x.v + dv;
    x.p = x.p + dp;
    x.b_g = x.b_g + dbg;
    x.b_a = x.b_a + dba;

    // Solà right-multiplicative reset: G_θθ ≈ I − ½[δθ]×  (arXiv:1711.02508 §5.4)
    Matrix<kInsErrorDim, kInsErrorDim, T> G = Matrix<kInsErrorDim, kInsErrorDim, T>::identity();
    const T                               n2 = dth[0] * dth[0] + dth[1] * dth[1] + dth[2] * dth[2];
    if (n2 > static_cast<T>(1e-24)) {
        G.template block<3, 3>(ins_err::dtheta, ins_err::dtheta) = Matrix<3, 3, T>::identity() - (detail::skew_sym(dth) * static_cast<T>(0.5));
    }
    return G;
}

/**
 * @brief Mechanize nominal state and propagate the error covariance
 *
 * @param eskf  15-state filter (Q must be set)
 * @param x     nominal INS state (updated by mechanization)
 * @param imu   body-frame IMU sample
 * @param dt    [s]
 * @param frame NED or ENU
 * @param g     gravity magnitude
 */
template<typename T = double>
constexpr void ins_predict(
    ErrorStateKalmanFilter<kInsErrorDim, kInsMeasDim, T>& eskf,
    InsState<T>&                                          x,
    const ImuSample<T>&                                   imu,
    T                                                     dt,
    NavFrame                                              frame,
    T                                                     g = kStandardGravity<T>
) {
    const Vec3<T> omega = imu.gyro - x.b_g;
    const Vec3<T> a_b = imu.accel - x.b_a;
    // Jacobian at pre-step orientation (standard first-order practice)
    const auto ej = ins_error_jacobian(x.q, omega, a_b, dt);
    x = mechanize_step_from_corrected(x, omega, a_b, dt, frame, g);
    // Sparse F P Fᵀ + Q (G = I); avoids full 15×15×15 in ErrorStateKalmanFilter::predict
    auto P = eskf.covariance();
    ins_propagate_covariance(P, ej, eskf.process_noise_covariance());
    eskf.set_covariance(P);
}

/**
 * @brief Inject filter correction into @p x and reset the error state
 */
template<typename T = double>
constexpr void ins_apply_correction(
    ErrorStateKalmanFilter<kInsErrorDim, kInsMeasDim, T>& eskf,
    InsState<T>&                                          x
) {
    const auto G = ins_inject(x, eskf.error_state());
    eskf.reset_error_state(G);
}

/**
 * @brief Predicted antenna / marker position with body lever-arm
 *
 * @f$ p_{\mathrm{aid}} = p + R(q)\, \ell_{\mathrm{body}} @f$
 */
template<typename T = double>
[[nodiscard]] constexpr Vec3<T> ins_aid_position(
    const InsState<T>& x,
    const Vec3<T>&     lever_body = Vec3<T>{}
) {
    return x.p + x.q.rotate(lever_body);
}

/**
 * @brief Horizontal heading [rad] of a body axis expressed in the nav frame
 *
 * Default body axis is +x (vehicle forward). Dual-antenna baseline should match
 * the body vector used here (see @ref heading_from_baseline_nav).
 *
 * @param x          nominal state
 * @param frame      NED or ENU (defines atan2 axes)
 * @param body_axis  unit (or any non-zero) direction in the body frame
 */
template<typename T = double>
[[nodiscard]] constexpr T ins_heading(
    const InsState<T>& x,
    NavFrame           frame,
    const Vec3<T>&     body_axis = Vec3<T>{T{1}, T{0}, T{0}}
) {
    return detail::heading_from_forward_nav(x.q.rotate(body_axis), frame);
}

/**
 * @brief Heading [rad] from a nav-frame baseline vector (dual-antenna difference)
 *
 * @param baseline_nav  p_secondary − p_primary in the nav frame (along the body baseline)
 * @param frame         NED or ENU
 */
template<typename T = double>
[[nodiscard]] constexpr T heading_from_baseline_nav(const Vec3<T>& baseline_nav, NavFrame frame) {
    return detail::heading_from_forward_nav(baseline_nav, frame);
}

/**
 * @brief Sparse position update: y = p + Rℓ
 *
 * @param eskf        filter
 * @param x           nominal state (corrected on success)
 * @param p_meas      measured position of the antenna / marker (nav frame)
 * @param R_pos       3×3 measurement covariance
 * @param lever_body  body-frame lever-arm from IMU origin to the antenna (default zero)
 * @return false if the innovation covariance solve fails
 */
template<typename T = double>
[[nodiscard]] constexpr bool ins_update_position(
    ErrorStateKalmanFilter<kInsErrorDim, kInsMeasDim, T>& eskf,
    InsState<T>&                                          x,
    const Vec3<T>&                                        p_meas,
    const Matrix<kInsMeasDim, kInsMeasDim, T>&            R_pos,
    const Vec3<T>&                                        lever_body = Vec3<T>{}
) {
    const Vec3<T>          p_pred = ins_aid_position(x, lever_body);
    ColVec<kInsMeasDim, T> innov{
        p_meas[0] - p_pred[0],
        p_meas[1] - p_pred[1],
        p_meas[2] - p_pred[2],
    };

    // H_p = I on δp; H_θ = −R[ℓ]× for right-multiplicative attitude error
    Matrix<kInsMeasDim, kInsErrorDim, T> H = detail::selector_H<T>(ins_err::dp);
    const T                              ell2 = lever_body[0] * lever_body[0] + lever_body[1] * lever_body[1] + lever_body[2] * lever_body[2];
    if (ell2 > static_cast<T>(1e-24)) {
        const Mat3<T> R = x.q.to_dcm().matrix();
        H.template block<3, 3>(0, ins_err::dtheta) = (R * detail::skew_sym(lever_body)) * (-T{1});
    }

    if (!eskf.update_innovation(innov, H, R_pos)) {
        return false;
    }
    ins_apply_correction(eskf, x);
    return true;
}

/// Position update using the filter's stored R (default position noise from design).
template<typename T = double>
[[nodiscard]] constexpr bool ins_update_position(
    ErrorStateKalmanFilter<kInsErrorDim, kInsMeasDim, T>& eskf,
    InsState<T>&                                          x,
    const Vec3<T>&                                        p_meas,
    const Vec3<T>&                                        lever_body = Vec3<T>{}
) {
    return ins_update_position(eskf, x, p_meas, eskf.measurement_noise_covariance(), lever_body);
}

/**
 * @brief Sparse velocity update: y = v at the IMU origin (no lever-arm)
 *
 * @f$ H @f$ selects @f$ \delta v @f$ only. Body-lever / @f$ \omega\times R\ell @f$ is not modeled.
 */
template<typename T = double>
[[nodiscard]] constexpr bool ins_update_velocity(
    ErrorStateKalmanFilter<kInsErrorDim, kInsMeasDim, T>& eskf,
    InsState<T>&                                          x,
    const Vec3<T>&                                        v_meas,
    const Matrix<kInsMeasDim, kInsMeasDim, T>&            R_vel
) {
    ColVec<kInsMeasDim, T> innov{v_meas[0] - x.v[0], v_meas[1] - x.v[1], v_meas[2] - x.v[2]};
    const auto             H = detail::selector_H<T>(ins_err::dv);
    if (!eskf.update_innovation(innov, H, R_vel)) {
        return false;
    }
    ins_apply_correction(eskf, x);
    return true;
}

/**
 * @brief Zero-velocity update (ZUPT): velocity measurement of zero
 */
template<typename T = double>
[[nodiscard]] constexpr bool ins_update_zupt(
    ErrorStateKalmanFilter<kInsErrorDim, kInsMeasDim, T>& eskf,
    InsState<T>&                                          x,
    const Matrix<kInsMeasDim, kInsMeasDim, T>&            R_zupt
) {
    return ins_update_velocity(eskf, x, Vec3<T>{T{0}, T{0}, T{0}}, R_zupt);
}

/**
 * @brief Sparse dual-antenna / yaw heading update (scalar)
 *
 * Innovation is the wrapped heading residual. Jacobian couples into @f$ \delta\theta @f$
 * via the body baseline direction (default body +x).
 *
 * @param eskf       filter
 * @param x          nominal state
 * @param psi_meas   measured heading [rad] (same convention as @ref ins_heading)
 * @param R_psi      1×1 heading variance [rad²]
 * @param frame      NED or ENU
 * @param body_axis  body direction of the antenna baseline
 */
template<typename T = double>
[[nodiscard]] constexpr bool ins_update_heading(
    ErrorStateKalmanFilter<kInsErrorDim, kInsMeasDim, T>& eskf,
    InsState<T>&                                          x,
    T                                                     psi_meas,
    const Matrix<1, 1, T>&                                R_psi,
    NavFrame                                              frame,
    const Vec3<T>&                                        body_axis = Vec3<T>{T{1}, T{0}, T{0}}
) {
    const Vec3<T> f = x.q.rotate(body_axis);
    const T       psi_pred = detail::heading_from_forward_nav(f, frame);
    ColVec<1, T>  innov{damp::wrap_pi(psi_meas - psi_pred)};

    // ∂f/∂δθ = −R [e]×  (right mult); H_θ = (∂ψ/∂f)(∂f/∂δθ)
    const Mat3<T>              R = x.q.to_dcm().matrix();
    const Mat3<T>              dfdth = (R * detail::skew_sym(body_axis)) * (-T{1});
    const auto                 dpsidf = detail::dheading_dforward(f, frame);
    Matrix<1, kInsErrorDim, T> H = Matrix<1, kInsErrorDim, T>::zeros();
    const Matrix<1, 3, T>      H_th = dpsidf * dfdth;
    H.template block<1, 3>(0, ins_err::dtheta) = H_th;

    if (!eskf.update_innovation(innov, H, R_psi)) {
        return false;
    }
    ins_apply_correction(eskf, x);
    return true;
}

/**
 * @brief Sparse orientation update from a measured body→nav quaternion
 *
 * Right-multiplicative residual: @f$ \delta\theta = \log(q^{-1} \otimes q_{\mathrm{meas}}) @f$.
 *
 * @param eskf   INS error-state filter (updated in place)
 * @param x      nominal INS state (corrected after a successful update)
 * @param q_meas measured body→nav quaternion
 * @param R_att  3×3 attitude measurement covariance [rad²]
 */
template<typename T = double>
[[nodiscard]] constexpr bool ins_update_orientation(
    ErrorStateKalmanFilter<kInsErrorDim, kInsMeasDim, T>& eskf,
    InsState<T>&                                          x,
    const Quaternion<T>&                                  q_meas,
    const Matrix<3, 3, T>&                                R_att
) {
    const Quaternion<T>        dq = (x.q.conjugate() * q_meas).normalized();
    const Vec3<T>              dth = dq.log();
    ColVec<3, T>               innov{dth[0], dth[1], dth[2]};
    Matrix<3, kInsErrorDim, T> H = Matrix<3, kInsErrorDim, T>::zeros();
    H.template block<3, 3>(0, ins_err::dtheta) = Matrix<3, 3, T>::identity();

    if (!eskf.update_innovation(innov, H, R_att)) {
        return false;
    }
    ins_apply_correction(eskf, x);
    return true;
}

/**
 * @brief Generic pose alias: position then orientation (two sparse updates)
 *
 * Camera / IR / mocap-style absolute pose. Lever-arm applies to the position half.
 */
template<typename T = double>
[[nodiscard]] constexpr bool ins_update_pose(
    ErrorStateKalmanFilter<kInsErrorDim, kInsMeasDim, T>& eskf,
    InsState<T>&                                          x,
    const Vec3<T>&                                        p_meas,
    const Quaternion<T>&                                  q_meas,
    const Matrix<kInsMeasDim, kInsMeasDim, T>&            R_pos,
    const Matrix<3, 3, T>&                                R_att,
    const Vec3<T>&                                        lever_body = Vec3<T>{}
) {
    if (!ins_update_position(eskf, x, p_meas, R_pos, lever_body)) {
        return false;
    }
    return ins_update_orientation(eskf, x, q_meas, R_att);
}

/**
 * @brief Pose alias: position + dual-antenna heading (two sparse updates)
 */
template<typename T = double>
[[nodiscard]] constexpr bool ins_update_pose_heading(
    ErrorStateKalmanFilter<kInsErrorDim, kInsMeasDim, T>& eskf,
    InsState<T>&                                          x,
    const Vec3<T>&                                        p_meas,
    T                                                     psi_meas,
    const Matrix<kInsMeasDim, kInsMeasDim, T>&            R_pos,
    const Matrix<1, 1, T>&                                R_psi,
    NavFrame                                              frame,
    const Vec3<T>&                                        lever_body = Vec3<T>{},
    const Vec3<T>&                                        body_axis = Vec3<T>{T{1}, T{0}, T{0}}
) {
    if (!ins_update_position(eskf, x, p_meas, R_pos, lever_body)) {
        return false;
    }
    return ins_update_heading(eskf, x, psi_meas, R_psi, frame, body_axis);
}

// =============================================================================
// Tier 3 — InsNavigator (Design → Deploy runtime)
// =============================================================================

/**
 * @brief Runtime strapdown navigator: nominal InsState + 15-state ESKF
 *
 * Lightweight ISR-friendly wrapper. Design Q/R/P0 offline with
 * @ref design::ins_eskf_design, convert with `.as<float>()`, then construct
 * once (`constinit` on target). Per tick: @ref predict then optional aiding.
 *
 * @tparam T Scalar type (default float for embedded)
 *
 * @see design::ins_eskf_design
 * @see ins_predict
 */
template<typename T = float>
class InsNavigator {
public:
    constexpr InsNavigator() = default;

    /**
     * @brief Construct from a design result
     * @param design  Q, R (default position), P0
     * @param frame   local-level frame for gravity / heading
     * @param g       gravity magnitude
     */
    constexpr explicit InsNavigator(
        const design::InsEskfResult<T>& design,
        NavFrame                        frame = NavFrame::NED,
        T                               g = kStandardGravity<T>
    )
        : eskf_(design.P0, design.Q, design.R),
          frame_(frame),
          g_(g) {
        x_.q = Quaternion<T>::identity();
    }

    constexpr InsNavigator(
        const Matrix<kInsErrorDim, kInsErrorDim, T>& P0,
        const Matrix<kInsErrorDim, kInsErrorDim, T>& Q,
        const Matrix<kInsMeasDim, kInsMeasDim, T>&   R,
        NavFrame                                     frame = NavFrame::NED,
        T                                            g = kStandardGravity<T>
    )
        : eskf_(P0, Q, R),
          frame_(frame),
          g_(g) {
        x_.q = Quaternion<T>::identity();
    }

    template<typename U>
    constexpr explicit InsNavigator(const InsNavigator<U>& other)
        : eskf_(other.filter()),
          x_{},
          frame_(other.frame()),
          g_(static_cast<T>(other.gravity())) {
        x_.q = Quaternion<T>{other.state().q};
        x_.v = Vec3<T>{other.state().v.template as<T>()};
        x_.p = Vec3<T>{other.state().p.template as<T>()};
        x_.b_g = Vec3<T>{other.state().b_g.template as<T>()};
        x_.b_a = Vec3<T>{other.state().b_a.template as<T>()};
    }

    /// One IMU step: mechanize nominal + propagate covariance.
    constexpr void predict(const ImuSample<T>& imu, T dt) {
        ins_predict(eskf_, x_, imu, dt, frame_, g_);
    }

    [[nodiscard]] constexpr bool update_position(
        const Vec3<T>& p_meas,
        const Vec3<T>& lever_body = Vec3<T>{}
    ) {
        return ins_update_position(eskf_, x_, p_meas, lever_body);
    }

    [[nodiscard]] constexpr bool update_position(
        const Vec3<T>&                             p_meas,
        const Matrix<kInsMeasDim, kInsMeasDim, T>& R_pos,
        const Vec3<T>&                             lever_body = Vec3<T>{}
    ) {
        return ins_update_position(eskf_, x_, p_meas, R_pos, lever_body);
    }

    [[nodiscard]] constexpr bool update_velocity(
        const Vec3<T>&                             v_meas,
        const Matrix<kInsMeasDim, kInsMeasDim, T>& R_vel
    ) {
        return ins_update_velocity(eskf_, x_, v_meas, R_vel);
    }

    [[nodiscard]] constexpr bool update_zupt(const Matrix<kInsMeasDim, kInsMeasDim, T>& R_zupt) {
        return ins_update_zupt(eskf_, x_, R_zupt);
    }

    [[nodiscard]] constexpr bool update_heading(
        T                      psi_meas,
        const Matrix<1, 1, T>& R_psi,
        const Vec3<T>&         body_axis = Vec3<T>{T{1}, T{0}, T{0}}
    ) {
        return ins_update_heading(eskf_, x_, psi_meas, R_psi, frame_, body_axis);
    }

    [[nodiscard]] constexpr bool update_orientation(
        const Quaternion<T>&   q_meas,
        const Matrix<3, 3, T>& R_att
    ) {
        return ins_update_orientation(eskf_, x_, q_meas, R_att);
    }

    [[nodiscard]] constexpr bool update_pose(
        const Vec3<T>&                             p_meas,
        const Quaternion<T>&                       q_meas,
        const Matrix<kInsMeasDim, kInsMeasDim, T>& R_pos,
        const Matrix<3, 3, T>&                     R_att,
        const Vec3<T>&                             lever_body = Vec3<T>{}
    ) {
        return ins_update_pose(eskf_, x_, p_meas, q_meas, R_pos, R_att, lever_body);
    }

    [[nodiscard]] constexpr bool update_pose_heading(
        const Vec3<T>&                             p_meas,
        T                                          psi_meas,
        const Matrix<kInsMeasDim, kInsMeasDim, T>& R_pos,
        const Matrix<1, 1, T>&                     R_psi,
        const Vec3<T>&                             lever_body = Vec3<T>{},
        const Vec3<T>&                             body_axis = Vec3<T>{T{1}, T{0}, T{0}}
    ) {
        return ins_update_pose_heading(
            eskf_, x_, p_meas, psi_meas, R_pos, R_psi, frame_, lever_body, body_axis
        );
    }

    [[nodiscard]] constexpr const InsState<T>& state() const { return x_; }
    [[nodiscard]] constexpr InsState<T>&       state() { return x_; }
    [[nodiscard]] constexpr const auto&        filter() const { return eskf_; }
    [[nodiscard]] constexpr auto&              filter() { return eskf_; }
    [[nodiscard]] constexpr NavFrame           frame() const { return frame_; }
    [[nodiscard]] constexpr T                  gravity() const { return g_; }

    [[nodiscard]] constexpr T heading(const Vec3<T>& body_axis = Vec3<T>{T{1}, T{0}, T{0}}) const {
        return ins_heading(x_, frame_, body_axis);
    }

    constexpr void set_state(const InsState<T>& x) { x_ = x; }
    constexpr void set_frame(NavFrame frame) { frame_ = frame; }

    /**
     * @brief Re-initialize nominal INS state and error-state covariance
     *
     * Clears the ESKF error state and sets P to @p P0 (default identity).
     * Does not change Q/R, frame, or gravity.
     */
    constexpr void reset(
        const InsState<T>&                           x0 = InsState<T>{},
        const Matrix<kInsErrorDim, kInsErrorDim, T>& P0 = Matrix<kInsErrorDim, kInsErrorDim, T>::identity()
    ) {
        x_ = x0;
        eskf_.set_error_state(ColVec<kInsErrorDim, T>{});
        eskf_.set_covariance(P0);
    }

    template<typename U>
    [[nodiscard]] constexpr InsNavigator<U> as() const {
        return InsNavigator<U>(*this);
    }

private:
    ErrorStateKalmanFilter<kInsErrorDim, kInsMeasDim, T> eskf_{};
    InsState<T>                                          x_{};
    NavFrame                                             frame_{NavFrame::NED};
    T                                                    g_{kStandardGravity<T>};
};

} // namespace damp
