// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file toolbox/encoder.hpp
 * @brief Quadrature decoding and speed (tachometer) measurement.
 *
 * Software A/B quadrature decode and pulse-based speed measurement — the
 * incremental-encoder front end for motion control. `constexpr`, allocation-free.
 * Counter wrap is handled deliberately (see wrapped_delta): position deltas
 * are computed in the unsigned domain and reinterpreted as signed so a rollover
 * never produces a spurious jump.
 *
 * @see motor_control.hpp for the Clarke/Park transforms these feed.
 */

#include <cstdint>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/math/math.hpp" // damp::wrap, numbers

namespace damp {

/**
 * @brief Signed difference between two unsigned counter readings, wrap-safe.
 *
 * Computes `curr − prev` in modular (unsigned) arithmetic and reinterprets the
 * result as signed, so a hardware counter rolling over its range yields the
 * correct small signed delta rather than a huge jump. Valid while the true
 * motion between reads is less than half the counter range.
 *
 * @tparam U Unsigned counter type (e.g. `uint16_t`, `uint32_t`).
 * @return The signed change, in the matching signed type.
 */
template<typename U>
[[nodiscard]] constexpr std::make_signed_t<U> wrapped_delta(U prev, U curr) {
    static_assert(std::is_unsigned_v<U>, "wrapped_delta requires an unsigned counter type");
    return static_cast<std::make_signed_t<U>>(static_cast<U>(curr - prev));
}

/// Quadrature decode resolution (edges counted per A/B cycle).
enum class QuadMode {
    X1, ///< Count one edge per cycle.
    X2, ///< Count A edges on both transitions.
    X4  ///< Count every A and B transition (full resolution).
};

/**
 * @brief Software A/B quadrature decoder with optional index.
 *
 * Feed the raw A and B channel levels each sample; the decoder tracks the
 * 2-bit Gray-coded state and accumulates a signed position. Direction follows
 * the standard quadrature convention (A leading B counts up). Invalid
 * transitions (both bits changing at once — a missed sample) are ignored rather
 * than corrupting the count. The position is a fixed-width signed integer and
 * wraps modulo its range by design; use wrapped_delta on it for rates.
 *
 * @code
 * QuadratureDecoder dec{QuadMode::X4};
 * // each sample: dec.update(a_pin, b_pin);
 * int32_t pos = dec.position();
 * @endcode
 */
class QuadratureDecoder {
public:
    constexpr QuadratureDecoder() = default;
    constexpr explicit QuadratureDecoder(QuadMode mode) : mode_(mode) {}

    /// Process one sample of the A and B channels; returns the new position.
    constexpr int32_t update(bool a, bool b) {
        const uint8_t state = static_cast<uint8_t>((static_cast<uint8_t>(a) << 1) | static_cast<uint8_t>(b));
        // Transition index: (previous << 2) | current, into a 16-entry table of
        // {-1, 0, +1} (0 for no-change and invalid double transitions).
        const uint8_t idx = static_cast<uint8_t>((prev_ << 2) | state);
        const int8_t  step = kTransition[idx];

        switch (mode_) {
            case QuadMode::X4:
                position_ += step;
                break;
            case QuadMode::X2:
                // Count only transitions of A (the high bit changed).
                if (((prev_ ^ state) & 0b10) != 0) {
                    position_ += step;
                }
                break;
            case QuadMode::X1:
                // Count one detent: a rising A edge.
                if ((prev_ & 0b10) == 0 && (state & 0b10) != 0) {
                    position_ += step;
                }
                break;
        }
        prev_ = state;
        return position_;
    }

    /// Latch the index/zero mark: reset position to 0 when @p z is true.
    constexpr void index(bool z) {
        if (z) {
            position_ = 0;
        }
    }

    [[nodiscard]] constexpr int32_t position() const { return position_; }

    constexpr void reset() {
        position_ = 0;
        prev_ = 0;
    }

private:
    // Quadrature transition table indexed by (prev<<2 | curr); +1 forward,
    // -1 reverse, 0 for no movement or an illegal (both-bits) transition.
    static constexpr damp::array<int8_t, 16> kTransition{
        0, +1, -1, 0,
        -1, 0, 0, +1,
        +1, 0, 0, -1,
        0, -1, +1, 0
    };

    QuadMode mode_{QuadMode::X4};
    uint8_t  prev_{0};
    int32_t  position_{0};
};

/**
 * @brief Pulse-based speed (tachometer) with frequency/period crossover.
 *
 * Recovers rotational speed from encoder pulses two ways and picks the better:
 * counting pulses over a fixed interval (accurate at high speed) and measuring
 * the period between edges (accurate at low speed, where few pulses land in an
 * interval). Configure the pulses per revolution; read out in your unit of
 * choice.
 *
 * @tparam T Scalar type (default: float)
 */
template<typename T = float>
class Tachometer {
public:
    T counts_per_rev{T{1}}; ///< Encoder counts (edges) per shaft revolution.

    constexpr Tachometer() = default;
    constexpr explicit Tachometer(T cpr) : counts_per_rev(cpr) {}

    /// Speed [rev/s] from a count delta over interval @p dt (frequency method).
    [[nodiscard]] constexpr T rev_per_s_from_counts(int32_t delta_counts, T dt) const {
        return (dt > T{0}) ? (static_cast<T>(delta_counts) / counts_per_rev) / dt : T{0};
    }

    /// Speed [rev/s] from the time @p edge_period between two successive counts.
    [[nodiscard]] constexpr T rev_per_s_from_period(T edge_period) const {
        return (edge_period > T{0}) ? T{1} / (counts_per_rev * edge_period) : T{0};
    }

    /// Convert a rev/s reading to RPM.
    [[nodiscard]] static constexpr T to_rpm(T rev_per_s) { return rev_per_s * T{60}; }

    /// Convert a rev/s reading to rad/s.
    [[nodiscard]] static constexpr T to_rad_per_s(T rev_per_s) {
        return rev_per_s * T{2} * damp::numbers::pi_v<T>;
    }
};

/**
 * @brief Wrap-safe accumulator: successive wrapped absolute readings → continuous position.
 *
 * An absolute angle sensor (magnetic encoder, resolver, sin/cos) reports a value that wraps
 * back to the start of its range every revolution. Naively differencing two readings across
 * that seam produces a full-scale jump; the classic bug. This turns the wrapped stream into a
 * continuous multi-turn position by folding each per-sample change onto the shortest arc,
 * @f[
 *   \Delta = \operatorname{wrap}\!\big(x_k - x_{k-1},\ -\tfrac{S}{2},\ \tfrac{S}{2}\big),
 *   \qquad p_k = p_{k-1} + \Delta,
 * @f]
 * for a wrap period @p full_scale @f$ S @f$ — 1.0 for turns, @f$ 2\pi @f$ for radians, or the
 * counts-per-revolution for a raw integer-scaled reading. Valid while the true motion between
 * samples is under half a period (the same Nyquist limit as wrapped_delta, its integer
 * sibling). No filtering or dynamics — pair with a velocity estimator or LPF if
 * you need a smoothed rate as well as continuous position.
 *
 * @tparam T Scalar type (default: float).
 */
template<typename T = float>
class AngleUnwrapper {
public:
    T full_scale{T{1}}; ///< wrap period S (1 = turns, 2π = radians, CPR = raw counts)

    constexpr AngleUnwrapper() = default;
    constexpr explicit AngleUnwrapper(T full_scale_in) : full_scale(full_scale_in) {}

    /**
     * @brief Feed one wrapped absolute reading; returns the continuous position.
     *
     * The first call after construction/reset seeds the reference without accumulating,
     * so an arbitrary starting phase does not inject a spurious first step.
     *
     * @param wrapped absolute reading in one period (any range; only differences matter)
     * @return continuous (multi-period) position in the same units
     */
    constexpr T update(T wrapped) {
        if (!started_) {
            prev_ = wrapped;
            position_ = wrapped;
            started_ = true;
            return position_;
        }
        const T half = full_scale / T{2};
        const T delta = wrap(wrapped - prev_, -half, half);
        position_ += delta;
        prev_ = wrapped;
        return position_;
    }

    [[nodiscard]] constexpr T position() const { return position_; } ///< continuous multi-period position

    constexpr void reset() {
        prev_ = T{0};
        position_ = T{0};
        started_ = false;
    }

private:
    T    prev_{T{0}};     ///< previous wrapped reading
    T    position_{T{0}}; ///< accumulated continuous position
    bool started_{false}; ///< false until the first sample seeds the reference
};

} // namespace damp
