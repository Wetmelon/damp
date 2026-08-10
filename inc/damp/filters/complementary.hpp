// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file complementary.hpp
 * @brief Scalar complementary filter runtime
 */


namespace damp {
/**
 * @brief Scalar (1-D) complementary filter — fuse a fast rate with a slow absolute.
 *
 * The classic two-input sensor blend: high-pass an integrated *rate* signal
 * (low drift, but accumulates error) and low-pass an *absolute* measurement
 * (no drift, but noisy), crossing them over at `1/τ`:
 *
 *   y ← α·(y + rate·dt) + (1−α)·measurement,   α = τ / (τ + dt).
 *
 * The textbook example is tilt from a gyro (`rate`) and accelerometer
 * (`measurement`); also altitude from baro + vertical accel, etc. For full
 * 3-D orientation use @ref ComplementaryFilter / @ref MahonyFilter instead —
 * this is the cheap scalar case.
 *
 * The first sample seeds the state to `measurement` (no start-up ramp).
 */
template<typename T = float>
class Complementary {
public:
    constexpr Complementary() = default;

    /// @param tau Crossover time constant [s] (larger ⇒ trust the rate longer).
    constexpr explicit Complementary(T tau) : tau_(tau) {}

    /// @param measurement Absolute (slow/noisy) reading. @param rate Its
    /// derivative (fast/low-drift). @param dt Sample time [s].
    constexpr T operator()(T measurement, T rate, T dt) {
        if (!init_) {
            y_ = measurement;
            init_ = true;
            return y_;
        }
        const T alpha = tau_ / (tau_ + dt);
        y_ = (alpha * (y_ + (rate * dt))) + ((T{1} - alpha) * measurement);
        return y_;
    }

    [[nodiscard]] constexpr T value() const { return y_; }

    constexpr void reset() {
        y_ = T{0};
        init_ = false;
    }

private:
    T    tau_{T{1}};
    T    y_{T{0}};
    bool init_{false};
};
} // namespace damp
