// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file damp/math/split_real.hpp
 * @brief Real number with an exact integer part: constant resolution over any range.
 */

#include "damp/backend.hpp"
#include "damp/math/math.hpp" // damp::floor, damp::ceil, damp::abs, damp::isfinite

namespace damp {

/**
 * @brief A real number split into an exact integer part and a floating-point remainder.
 *
 * For quantities that accumulate small steps over a large range: a multi-turn shaft
 * angle (turns), a multi-period encoder or linear-scale position (periods), odometry
 * (metres), metered energy (joules over decades). A plain float keeps relative
 * precision, so it loses absolute resolution every time the value doubles — and once
 * the total dwarfs the step, adding the step changes nothing. Here the integer part
 * is exact and the remainder stays in [−0.5, 0.5), so the resolution is the
 * remainder's (≈6·10⁻⁸ units for float) whatever the magnitude, up to ±2⁶³ units
 * (the integer part wraps modulo 2⁶⁴ past that, like a hardware counter).
 *
 * It is an affine type, like a time point: the difference of two values is a plain
 * @c T, and a value plus a @c T is a value. Keep differences, not values, in vectors
 * and filter states — estimators hold their state relative to an anchor and
 * controllers work on errors — so the type never needs to be a matrix element. A
 * difference is exact while the two values are within 2²⁴ units (float) of each other.
 *
 * The remainder's range matches a wrapped angle in turns: for a rotary axis counted
 * in turns, @c frac is the single-turn angle.
 *
 * A plain @c T converts implicitly (no information is lost, as float → double), so a
 * reference or setpoint near zero can be written as a number. It is deliberately not an
 * aggregate: brace elision would otherwise let `{0, v, t}` fill @c whole and @c frac of
 * a SplitReal member from a struct's later fields.
 *
 * @tparam T Scalar type of the remainder (float on target).
 */
template<typename T = float>
struct SplitReal {
    int64_t whole{0}; ///< exact integer part
    T       frac{0};  ///< remainder in [−0.5, 0.5)

    constexpr SplitReal() noexcept = default;

    /** @brief From a plain value: split into integer part + remainder. Implicit by design (lossless). */
    constexpr SplitReal(T value) noexcept : frac(value) { normalize(); } // NOLINT(google-explicit-constructor)

    /** @brief From an integer part and a remainder (any range; normalized). */
    constexpr SplitReal(int64_t whole_in, T frac_in) noexcept : whole(whole_in), frac(frac_in) { normalize(); }

    /** @brief The value as one number. Lossy for large magnitudes; prefer differences. */
    [[nodiscard]] constexpr T value() const { return static_cast<T>(whole) + frac; }

    constexpr SplitReal& operator+=(T step) {
        frac += step;
        normalize();
        return *this;
    }
    constexpr SplitReal& operator-=(T step) { return *this += -step; }

    [[nodiscard]] friend constexpr SplitReal operator+(SplitReal r, T step) { return r += step; }
    [[nodiscard]] friend constexpr SplitReal operator-(SplitReal r, T step) { return r -= step; }

    /** @brief Difference @p a − @p b as a plain number. */
    [[nodiscard]] friend constexpr T operator-(const SplitReal& a, const SplitReal& b) {
        // Modular difference: defined for any two values (a signed subtraction could overflow).
        const auto dw = static_cast<int64_t>(static_cast<uint64_t>(a.whole) - static_cast<uint64_t>(b.whole));
        return static_cast<T>(dw) + (a.frac - b.frac);
    }

    [[nodiscard]] friend constexpr bool operator==(const SplitReal&, const SplitReal&) = default;

    template<typename U>
    [[nodiscard]] constexpr SplitReal<U> as() const {
        return {whole, static_cast<U>(frac)};
    }

private:
    /** @brief Move whole units out of the remainder so it lands in [−0.5, 0.5). */
    constexpr void normalize() {
        // A step this large is not an accumulation; leave it rather than cast an
        // out-of-range float to an integer (undefined behaviour). NaN and ±Inf are left too,
        // tested on the bit pattern: under -ffast-math a comparison may assume them away.
        constexpr T limit = static_cast<T>(int64_t{1} << 62);
        if (!damp::isfinite(frac) || !(damp::abs(frac) < limit)) {
            return;
        }
        // Exact split: frac − trunc(frac) rounds nothing for either sign, and the single fold
        // (r ∓ 1 with |r| in [0.5, 1)) is exact, so the value never moves and the remainder
        // lands in [−0.5, 0.5). (frac − floor(frac) would round a negative remainder onto
        // the grid of [0.5, 1) and drop small negative steps.)
        T n = (frac >= T{0}) ? damp::floor(frac) : damp::ceil(frac);
        T r = frac - n;
        if (r >= static_cast<T>(0.5)) {
            n += T{1};
            r -= T{1};
        } else if (r < static_cast<T>(-0.5)) {
            n -= T{1};
            r += T{1};
        }
        const auto step = static_cast<uint64_t>(static_cast<int64_t>(n));
        whole = static_cast<int64_t>(static_cast<uint64_t>(whole) + step); // modular: no signed overflow
        frac = r;
    }
};

} // namespace damp
