// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file damp/estimation/ins_mechanization.hpp
 * @brief Strapdown INS step: IMU samples → attitude, velocity, position (NED or ENU)
 *
 * ## What “mechanization” means
 *
 * In inertial navigation, mechanization is the open-loop physics recipe that
 * turns a stream of IMU readings into a navigation state — not a Kalman filter
 * and not GPS. Each step:
 *
 * 1. Integrate the gyro → update orientation
 * 2. Rotate the accelerometer (specific force) into the map frame
 * 3. Add gravity
 * 4. Integrate once → velocity, twice → position
 *
 * The textbook name is “strapdown inertial mechanization” (sensors fixed to the
 * vehicle; the rotation is done in software). Plain language: dead-reckoning
 * from gyro + accel. mechanize_step is one discrete sample of that recipe.
 *
 * There is no aiding here. Biases and noise make position drift without bound;
 * absolute corrections (GPS, camera/mocap pose, ZUPT, dual-antenna heading, …) are
 * the ESKF half of roadmap #26 and live outside this header.
 *
 * ## Model
 *
 * Bias-corrected body rates and specific force, gravity in the chosen NavFrame,
 * flat-Earth (no Earth rate / transport rate / WGS84).
 *
 * ## Frames
 *
 * | NavFrame | Axes (x, y, z) | Gravity @f$ g_n @f$ |
 * | ------------- | -------------- | ------------------- |
 * | NED | North, East, Down | @f$ (0,0,+g) @f$ |
 * | ENU | East, North, Up | @f$ (0,0,-g) @f$ |
 *
 * Orientation @c q is body → nav: @c q.rotate(v_body) yields nav components.
 * IMU samples are body-frame (rad/s, m/s² specific force). With
 * @f$ a_n = R(q)\,a_b + g_n @f$, rest requires @f$ a_b = -R(q)^\top g_n @f$
 * (specific_force_at_rest). For @c q = I: ENU @f$ a_b \approx (0,0,+g) @f$
 * (Up), NED @f$ a_b \approx (0,0,-g) @f$ (Up in NED is −Down). Do not feed
 * @f$ -g_n @f$ as body accel without that rotation.
 *
 * @see Solà, "Quaternion kinematics for the error-state Kalman filter," 2017
 * @see Groves, Principles of GNSS, Inertial, and Multisensor Integrated Navigation, 2nd ed.
 *
 * Example:
 * @code
 * InsState<double> x{.q = Quaternion<double>::identity()};
 * ImuSample<double> imu{
 *     .gyro = {0, 0, 0},
 *     .accel = specific_force_at_rest(x.q, NavFrame::ENU), // ≈ (0,0,+g), not −g
 * };
 * x = mechanize_step(x, imu, 0.01, NavFrame::ENU);
 * @endcode
 */

#include <cstdint>

#include "damp/math/geometry.hpp"
#include "damp/matrix/colvec.hpp"

namespace damp {

/**
 * @brief Local-level navigation frame for strapdown mechanization
 *
 * Selects gravity direction and axis labeling only. One code path for NED and ENU
 * (MATLAB-style option; default documentation uses NED).
 */
enum class NavFrame : std::uint8_t {
    NED, ///< North-East-Down — @f$ g_n = (0,0,+g) @f$
    ENU, ///< East-North-Up — @f$ g_n = (0,0,-g) @f$
};

/// Standard gravity magnitude [m/s²] (CGPM conventional value).
template<typename T = double>
inline constexpr T kStandardGravity = static_cast<T>(9.80665);

/**
 * @brief Gravity vector in the local-level frame [m/s²]
 *
 * @param frame NED or ENU
 * @param g     gravity magnitude (default @ref kStandardGravity)
 */
template<typename T = double>
[[nodiscard]] constexpr Vec3<T> gravity_nav(NavFrame frame, T g = kStandardGravity<T>) {
    if (frame == NavFrame::NED) {
        return Vec3<T>{T{0}, T{0}, g};
    }
    return Vec3<T>{T{0}, T{0}, -g};
}

/**
 * @brief Map an ENU vector into NED (axis permute)
 *
 * ENU @f$ (E,N,U) @f$ → NED @f$ (N,E,D) = (N,E,-U) @f$.
 */
template<typename T = double>
[[nodiscard]] constexpr Vec3<T> ned_from_enu(const Vec3<T>& enu) {
    return Vec3<T>{enu[1], enu[0], -enu[2]};
}

/**
 * @brief Map a NED vector into ENU (axis permute)
 *
 * NED @f$ (N,E,D) @f$ → ENU @f$ (E,N,U) = (E,N,-D) @f$.
 */
template<typename T = double>
[[nodiscard]] constexpr Vec3<T> enu_from_ned(const Vec3<T>& ned) {
    return Vec3<T>{ned[1], ned[0], -ned[2]};
}

/**
 * @brief IMU sample in the body frame
 * @tparam T Scalar type
 */
template<typename T = double>
struct ImuSample {
    Vec3<T> gyro{};  ///< [rad/s] angular rate (body)
    Vec3<T> accel{}; ///< [m/s²] specific force (body); includes gravity reaction at rest
};

/**
 * @brief Nominal strapdown navigation state (external to the ESKF error state)
 *
 * @tparam T Scalar type
 *
 * @c q maps body → nav (@c q.rotate(v_body) is in nav). Position/velocity are in
 * the same NavFrame as the mechanization call.
 */
template<typename T = double>
struct InsState {
    Quaternion<T> q{Quaternion<T>::identity()}; ///< body → nav orientation
    Vec3<T>       v{};                          ///< [m/s] velocity in nav
    Vec3<T>       p{};                          ///< [m] position in nav
    Vec3<T>       b_g{};                        ///< [rad/s] gyro bias (body)
    Vec3<T>       b_a{};                        ///< [m/s²] accelerometer bias (body)
};

/**
 * @brief One dead-reckoning step from an IMU sample (strapdown mechanization)
 *
 * Open-loop only: updates @c q, @c v, @c p from gyro/accel. No GPS or other aid.
 *
 * @f[
 *   \omega = \omega_m - b_g,\quad
 *   a_b = a_m - b_a,\quad
 *   a_n = R(q)\,a_b + g_n,\quad
 *   v^+ = v + a_n\,\Delta t,\quad
 *   p^+ = p + v\,\Delta t + \tfrac12 a_n (\Delta t)^2,\quad
 *   q^+ = q \otimes \Delta q(\omega\,\Delta t)
 * @f]
 *
 * Specific force: at rest @f$ a_b = -R^\top g_n @f$ (specific_force_at_rest).
 * With @c q = I: ENU @f$ a_b \approx (0,0,+g) @f$; NED @f$ a_b \approx (0,0,-g) @f$.
 *
 * @param x     state at the start of the interval
 * @param omega body-frame angular rate [rad/s] (bias-corrected)
 * @param a_b   body-frame specific force [m/s²] (bias-corrected)
 * @param dt    [s] step (must be > 0 for a useful step; dt ≤ 0 returns @p x)
 * @param frame local-level frame for @f$ g_n @f$ and for interpreting @c v, @c p
 * @param g     gravity magnitude [m/s²]
 * @return state at the end of the interval
 */
template<typename T = double>
[[nodiscard]] constexpr InsState<T> mechanize_step_from_corrected(
    const InsState<T>& x,
    const Vec3<T>&     omega,
    const Vec3<T>&     a_b,
    T                  dt,
    NavFrame           frame,
    T                  g = kStandardGravity<T>
) {
    if (!(dt > T{0})) {
        return x;
    }
    const Vec3<T> a_n = x.q.rotate(a_b) + gravity_nav(frame, g);
    return InsState<T>{
        .q = x.q.integrate_body_rates(omega, dt),
        .v = x.v + (a_n * dt),
        .p = (x.p + (x.v * dt)) + (a_n * (static_cast<T>(0.5) * (dt * dt))),
        .b_g = x.b_g,
        .b_a = x.b_a,
    };
}

template<typename T = double>
[[nodiscard]] constexpr InsState<T> mechanize_step(
    const InsState<T>&  x,
    const ImuSample<T>& imu,
    T                   dt,
    NavFrame            frame,
    T                   g = kStandardGravity<T>
) {
    return mechanize_step_from_corrected(x, imu.gyro - x.b_g, imu.accel - x.b_a, dt, frame, g);
}

/**
 * @brief Specific force [m/s²] a stationary IMU measures for the given orientation
 *
 * At rest @f$ a_n = 0 @f$ ⇒ @f$ a_b = -R(q)^\top g_n @f$ (body reading counters gravity).
 * Useful for synthetic IMU generation in tests and SIL.
 */
template<typename T = double>
[[nodiscard]] constexpr Vec3<T> specific_force_at_rest(
    const Quaternion<T>& q,
    NavFrame             frame,
    T                    g = kStandardGravity<T>
) {
    // a_b such that q.rotate(a_b) + g_n = 0  ⇒  a_b = q.conjugate().rotate(-g_n)
    const Vec3<T> g_n = gravity_nav(frame, g);
    return q.conjugate().rotate(Vec3<T>{-g_n[0], -g_n[1], -g_n[2]});
}

} // namespace damp
