// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file damp/encoder/tracker.hpp
 * @brief Encoder position (exact) and speed (smoothed) from absolute counts
 */

#include "damp/math/split_real.hpp"
#include "damp/matrix/colvec.hpp"

namespace damp {

/**
 * @brief Configuration for EncoderTracker.
 * @tparam T Scalar type.
 */
template<typename T = float>
struct EncoderTrackerConfig {
    T            bandwidth{T{2000}};          ///< [rad/s] velocity-loop bandwidth (critically damped); larger = faster, noisier
    ColVec<2, T> x0{};                        ///< initial state [position (counts), speed (counts/s)]
    bool         interpolate_position{false}; ///< position() is the loop's sub-count estimate instead of the last reading

    template<typename U>
    [[nodiscard]] constexpr EncoderTrackerConfig<U> as() const {
        return {
            .bandwidth = static_cast<U>(bandwidth),
            .x0 = ColVec<2, U>{static_cast<U>(x0[0]), static_cast<U>(x0[1])},
            .interpolate_position = interpolate_position,
        };
    }
};

/**
 * @brief Encoder position (exact) and speed (smoothed) from absolute counts.
 *
 * Position is the last reading, extrapolated by the speed estimate between readings:
 * @f$ \theta = \theta_{meas} + \omega\,t_{since} @f$. Speed comes from a second-order
 * tracking loop (the encoder PLL) on the reading, critically damped at
 * EncoderTrackerConfig::bandwidth ω_bw: @f$ k_p = 2\omega_{bw} @f$, @f$ k_i = \omega_{bw}^2 @f$.
 * The loop is type 2 (no lag at constant speed); the estimate is the shaft speed through
 * @f$ \omega_{bw}^2/(s + \omega_{bw})^2 @f$, a delay of about 2/ω_bw.
 *
 * Setting ω_bw:
 * - a cascade ratio above the speed-loop bandwidth (motor::CascadeBandwidths::omega_tracker);
 * - at most motor::max_loop_bandwidth(Ts) = 0.5/Ts when run once per tick (the loop goes
 *   unstable at ω_bw·Ts ≈ 0.83).
 *
 * With EncoderTrackerConfig::interpolate_position, position() is instead the loop's own
 * estimate: it moves smoothly within a count rather than in one-count steps, and lags by
 * a/kᵢ while the shaft accelerates at a. Use it when a position loop's gain on a one-count
 * step is too much torque ripple.
 *
 * Readings are absolute counts in a SplitReal (sub-count resolution at any travel). Unwrap
 * single-turn readings (PhaseUnwrapper) and validate them (EncoderGate) upstream; convert to
 * turns outside (AxisPosition). A quadrature count n means the shaft is in [n, n + 1), so
 * either position carries a constant −½ count offset; a calibration taken through the same
 * tracker (encoder offset, homing) absorbs it.
 *
 * @see MechanicalEstimator — the torque-driven variant that also estimates load.
 * @see ODrive's encoder PLL (fw-v0.5.6): the same gains; it refuses Ts·kₚ ≥ 1.
 *
 * @tparam T Scalar type.
 */
template<typename T = float>
class EncoderTracker {
public:
    /**
     * @param bandwidth [rad/s] velocity-loop bandwidth
     * @param x0        initial [position (counts), speed (counts/s)]
     */
    constexpr explicit EncoderTracker(T bandwidth = EncoderTrackerConfig<T>{}.bandwidth, ColVec<2, T> x0 = {}) noexcept
        : EncoderTracker<T>(EncoderTrackerConfig<T>{.bandwidth = bandwidth, .x0 = x0}) {}

    constexpr explicit EncoderTracker(const EncoderTrackerConfig<T>& config) noexcept
        : kp_(T{2} * config.bandwidth),
          ki_(config.bandwidth * config.bandwidth),
          meas_(config.x0[0]),
          est_(meas_),
          omega_(config.x0[1]),
          interpolate_(config.interpolate_position) {}

    /** @brief Set the velocity-loop bandwidth [rad/s] (critically damped). The "smoother/faster" knob. */
    constexpr void set_bandwidth(T bandwidth) {
        kp_ = T{2} * bandwidth;
        ki_ = bandwidth * bandwidth;
    }

    /** @return [rad/s] velocity-loop bandwidth (kₚ / 2) */
    [[nodiscard]] constexpr T bandwidth() const { return kp_ / T{2}; }

    /** @brief Advance time by @p dt [s] without a reading (cheap; run every control tick). */
    constexpr void predict(T dt) {
        est_ += omega_ * dt;
        since_ += dt;
    }

    /** @brief Take an absolute reading [counts]; @p dt [s] since the previous reading. */
    constexpr void update(const SplitReal<T>& pos_meas, T dt) {
        const T delta = pos_meas - est_;
        est_ += kp_ * delta * dt;
        omega_ += ki_ * delta * dt;
        meas_ = pos_meas;
        since_ = T{0};
    }

    /** @return [counts] the last reading extrapolated by ω since it was taken, or the loop's estimate when interpolating */
    [[nodiscard]] constexpr SplitReal<T> position() const {
        if (interpolate_) {
            return est_;
        }
        return meas_ + (omega_ * since_);
    }

    /** @return [counts/s] smoothed speed */
    [[nodiscard]] constexpr T omega() const { return omega_; }

    /** @brief Restart from [position (counts), speed (counts/s)]. */
    constexpr void reset(const ColVec<2, T>& x0 = ColVec<2, T>{}) { reset(SplitReal<T>{x0[0]}, x0[1]); }

    /** @brief Restart from an absolute position [counts] (e.g. after homing) and speed [counts/s]. */
    constexpr void reset(const SplitReal<T>& pos, T omega) {
        meas_ = pos;
        est_ = pos;
        omega_ = omega;
        since_ = T{0};
    }

private:
    T            kp_{T{2} * T{2000}};    ///< velocity-loop proportional gain = 2·bandwidth
    T            ki_{T{2000} * T{2000}}; ///< velocity-loop integral gain = bandwidth²
    SplitReal<T> meas_{};                ///< [counts] last reading
    SplitReal<T> est_{};                 ///< [counts] velocity loop's internal position estimate
    T            omega_{T{0}};           ///< [counts/s] smoothed speed
    T            since_{T{0}};           ///< [s] time since the last reading
    bool         interpolate_{false};    ///< position() returns est_
};

} // namespace damp
