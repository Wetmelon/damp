// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file toolbox/scaling.hpp
 * @brief Scaling, interpolation, and calibration helpers.
 *
 * The everyday "turn a raw number into an engineering value" primitives:
 * linear interpolation, range-to-range rescaling, affine sensor calibration,
 * and polynomial cal-curve evaluation. All are `constexpr`, allocation-free,
 * and have no state — drop them straight into an ISR.
 *
 * @see toolbox/lookup.hpp for breakpoint-table (non-affine) sensor linearization.
 */

#include <cstddef>
#include <type_traits>

#include "damp/matrix/matrix.hpp" // RowVec

namespace damp {

/**
 * @brief Linear interpolation between @p a and @p b by fraction @p t.
 *
 * Returns `a + (b − a)·t`. `t = 0` gives `a`, `t = 1` gives `b`; values outside
 * [0, 1] extrapolate (use clamp on `t` first if that is unwanted).
 */
template<typename T>
    requires std::is_floating_point_v<T>
[[nodiscard]] constexpr T lerp(T a, T b, T t) {
    return a + ((b - a) * t);
}

/**
 * @brief Inverse of lerp: the fraction @p t such that `lerp(a, b, t) == x`.
 *
 * Returns `(x − a) / (b − a)`. When `a == b` (degenerate range) returns `0`
 * rather than Inf/NaN — the fraction is undefined, so the lower endpoint is
 * reported.
 */
template<typename T>
    requires std::is_floating_point_v<T>
[[nodiscard]] constexpr T inverse_lerp(T a, T b, T x) {
    const T den = b - a;
    if (den == T{0}) {
        return T{0};
    }
    return (x - a) / den;
}

/**
 * @brief Affine map of @p x from the input range to the output range.
 *
 * Composes @ref inverse_lerp and lerp: maps `[in_lo, in_hi]` onto
 * `[out_lo, out_hi]` linearly. The typed, range-explicit replacement for the
 * Arduino® `map()` idiom; the input range must be non-degenerate (`in_lo != in_hi`).
 *
 * @code
 * float volts = rescale(adc_counts, 0.0f, 4095.0f, 0.0f, 3.3f);
 * @endcode
 */
template<typename T>
    requires std::is_floating_point_v<T>
[[nodiscard]] constexpr T rescale(T x, T in_lo, T in_hi, T out_lo, T out_hi) {
    return lerp(out_lo, out_hi, inverse_lerp(in_lo, in_hi, x));
}

/**
 * @brief Affine sensor calibration `y = gain·x + offset`.
 *
 * The standard two-coefficient linearization between a raw reading and an
 * engineering value. Build one directly, or from two known points with
 * @ref two_point_cal.
 *
 * @tparam T Scalar type (default: float)
 */
template<typename T = float>
struct AffineCal {
    T gain{1};   ///< Slope applied to the raw value.
    T offset{0}; ///< Constant added after scaling.

    /// Raw → engineering: `gain·x + offset`.
    [[nodiscard]] constexpr T apply(T x) const { return (gain * x) + offset; }

    /// Engineering → raw: `(y − offset) / gain` (requires `gain != 0`).
    [[nodiscard]] constexpr T invert(T y) const { return (y - offset) / gain; }

    /// Re-cast the coefficients to a different scalar type.
    template<typename U>
    [[nodiscard]] constexpr AffineCal<U> as() const {
        return AffineCal<U>{static_cast<U>(gain), static_cast<U>(offset)};
    }
};

/**
 * @brief Fit an AffineCal through two `(raw, engineering)` points.
 *
 * @param raw0 Raw reading at the first point.
 * @param eng0 Engineering value at the first point.
 * @param raw1 Raw reading at the second point.
 * @param eng1 Engineering value at the second point.
 * @return AffineCal reproducing both points (requires `raw0 != raw1`).
 */
template<typename T>
    requires std::is_floating_point_v<T>
[[nodiscard]] constexpr AffineCal<T> two_point_cal(T raw0, T eng0, T raw1, T eng1) {
    const T gain = (eng1 - eng0) / (raw1 - raw0);
    return AffineCal<T>{gain, eng0 - (gain * raw0)};
}

/**
 * @brief Evaluate a polynomial at @p x by Horner's method.
 *
 * Coefficients are in ascending power order: `coeffs[0]` is the constant
 * term, so the result is `coeffs[0] + coeffs[1]·x + coeffs[2]·x² + …`. Horner's
 * scheme is both the fewest-multiply and the most numerically stable evaluation.
 *
 * @code
 * // Steinhart-style cal curve c0 + c1·r + c2·r²:
 * constexpr RowVec<3, float> c{0.1f, 2.0e-3f, -5.0e-7f};
 * float eng = poly_horner(c, raw);
 * @endcode
 *
 * @tparam N Number of coefficients (polynomial degree N−1; must be ≥ 1).
 * @tparam T Scalar type.
 */
template<size_t N, typename T>
    requires std::is_floating_point_v<T>
[[nodiscard]] constexpr T poly_horner(const RowVec<N, T>& coeffs, T x) {
    static_assert(N >= 1, "poly_horner needs at least one coefficient");
    T result = coeffs[N - 1];
    for (size_t i = N - 1; i > 0; --i) {
        result = (result * x) + coeffs[i - 1];
    }
    return result;
}

} // namespace damp
