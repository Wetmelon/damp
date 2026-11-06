// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file damp/encoder/axis_position.hpp
 * @brief The position of an axis, rotary or linear: an exact count carrying its own scale.
 */

#include <cstdint>

#include "damp/math/math.hpp"       // damp::floor
#include "damp/math/split_real.hpp" // SplitReal

namespace damp {

/**
 * @brief The position of an axis, rotary or linear: an exact count carrying its own scale.
 *
 * The count is exact (a @ref SplitReal "SplitReal": integer count plus fractional count). The scale is
 * a period: @ref counts_per_period counts make one period, and one period is
 * @ref period_length in the axis's units. The same type covers every axis:
 *
 * | Axis | Period | counts_per_period | period_length |
 * | --- | --- | --- | --- |
 * | motor shaft, in turns | a turn | 2²³ | 1 |
 * | the same shaft as belt travel | a turn | 2²³ | 0.040 (m, the pulley lead) |
 * | glass scale, 50 nm counts | a count | 1 | 50e-9 (m) |
 * | glass scale by grating line | a 20 µm line | 400 | 20e-6 (m) |
 * | linear motor | a pole pitch | counts per pitch | the pitch (m) |
 *
 * Every view is derived on demand: @ref phase (where in the period, e.g. the commutation
 * angle) and @ref distance_to (an exact displacement in the axis's units). Because the scale
 * travels with the value, a consumer needs no encoder configuration of its own.
 *
 * Why counts: a float holds 24 significant bits. A 40M count/rev encoder needs 25.3 bits per
 * turn, so a float position cannot resolve every count near half a turn, and a tracker adding
 * sub-count corrections there loses them to rounding. Counted in counts, the float part never
 * exceeds half a count. Differences stay exact; a single absolute number (@ref value) cannot.
 *
 * @tparam T Scalar type of the fractional count and the units.
 */
template<typename T = float>
struct AxisPosition {
    SplitReal<T> counts{};             ///< [counts] exact position
    std::int64_t counts_per_period{1}; ///< counts in one period (> 0)
    T            period_length{T{1}};  ///< one period in the axis's units (1 turn, 0.04 m of belt, …)

    /**
     * @brief Where in the period, as a fraction in [−0.5, 0.5) (the commutation angle on a rotary or
     * linear motor). Whole periods are removed exactly, in integers, before converting to float.
     */
    [[nodiscard]] constexpr T phase() const { return wrapped(counts.whole, counts.frac, counts_per_period); }

    /** @brief Displacement from here to @p target, in the axis's units (an exact count difference, then scaled). */
    [[nodiscard]] constexpr T distance_to(const SplitReal<T>& target) const {
        return (target - counts) * (period_length / static_cast<T>(counts_per_period));
    }

    /**
     * @brief Shortest displacement from here to @p target modulo @p wrap_counts (circular position mode),
     * in the axis's units, in [−½, ½)·wrap. Whole wraps are removed exactly, in integers.
     */
    [[nodiscard]] constexpr T distance_to(const SplitReal<T>& target, std::int64_t wrap_counts) const {
        const T wrap_length = static_cast<T>(wrap_counts) * (period_length / static_cast<T>(counts_per_period));
        return wrapped(target.whole - counts.whole, target.frac - counts.frac, wrap_counts) * wrap_length;
    }

    /**
     * @brief The position as one number in the axis's units. Lossy far from zero (about 7 significant
     * digits in float); prefer @ref distance_to.
     */
    [[nodiscard]] constexpr T value() const { return counts.value() * (period_length / static_cast<T>(counts_per_period)); }

private:
    /** @brief (whole + frac) / period, wrapped to [−0.5, 0.5); the integer part is reduced first. */
    [[nodiscard]] static constexpr T wrapped(std::int64_t whole, T frac, std::int64_t period) {
        if (period <= 0) {
            return T{0}; // invalid scale; never divide by zero
        }
        std::int64_t r = whole % period;
        if (r > period / 2) {
            r -= period;
        } else if (r < -(period / 2)) {
            r += period;
        }
        const T x = (static_cast<T>(r) + frac) / static_cast<T>(period);
        return x - damp::floor(x + static_cast<T>(0.5));
    }
};

} // namespace damp
