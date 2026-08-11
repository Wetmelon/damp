// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file math.hpp
 * @brief Public scalar math: dispatches between compile-time and runtime paths.
 *
 * Each function uses std::is_constant_evaluated() to choose:
 * - At compile time: the constexpr series / Newton-Raphson bodies in
 *   constexpr_math.hpp (necessary for consteval design functions).
 * - At runtime: MathBackend\<T\> (platform-specific or the std:: fallback).
 *
 * This keeps consteval design code working while runtime controllers get full
 * hardware-math performance. Shared domain guards (asin/acos clamping, fmod's
 * divide-by-zero) live in the dispatcher so both paths agree by construction.
 *
 * @see constexpr_math.hpp for the compile-time implementations
 * @see math_backend.hpp for backend selection via damp_profile.hpp
 */

#include <type_traits>

#include "constexpr_math.hpp"
#include "damp/backend.hpp"
#include "math_backend.hpp"

namespace damp {

/**
 * @brief Square root.
 *
 * Domain x < 0 is guarded before the compile/runtime dispatch so both paths
 * return 0 (NaN is unavailable under -ffinite-math-only / constant evaluation).
 * For x ≥ 0: compile-time Newton-Raphson, runtime MathBackend\<T\>::sqrt.
 *
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
constexpr T sqrt(T x) {
    if (x < T{0}) {
        return T{0};
    }
    if (std::is_constant_evaluated()) {
        return detail::sqrt(x);
    }
    return MathBackend<T>::sqrt(x);
}

/**
 * @brief Absolute value, |x|.
 * @tparam T Numeric type
 */
template<typename T>
constexpr T abs(T x) {
    if (std::is_constant_evaluated()) {
        return detail::abs(x);
    }
    return MathBackend<T>::abs(x);
}

/**
 * @brief Euclidean distance hypot(x, y) = √(x² + y²), without overflow.
 *
 * Scales by the larger magnitude before squaring, so it stays finite when x² or
 * y² would overflow (and is exact at the axes). Backend-independent — built on
 * damp::sqrt/abs, so the same body serves compile and run time.
 *
 * @note Compare with MATLAB®'s hypot(x, y).
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
constexpr T hypot(T x, T y) {
    const T ax = damp::abs(x);
    const T ay = damp::abs(y);
    const auto [lo, hi] = damp::minmax(ax, ay);
    if (hi == T{0}) {
        return T{0};
    }
    const T r = lo / hi;
    return hi * damp::sqrt<T>(T{1} + (r * r)); // <T>: dispatcher, not trig.hpp's non-constexpr sqrt(float)
}

/**
 * @brief Cube root (preserves sign for negative x).
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
constexpr T cbrt(T x) {
    if (std::is_constant_evaluated()) {
        return detail::cbrt(x);
    }
    return MathBackend<T>::cbrt(x);
}

/**
 * @brief Two-argument arctangent, atan2(y, x) ∈ [−π, π].
 * @note Compare with MATLAB®'s atan2(y, x).
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
constexpr T atan2(T y, T x) {
    if (std::is_constant_evaluated()) {
        return detail::atan2(y, x);
    }
    return MathBackend<T>::atan2(y, x);
}

/**
 * @brief Single-argument arctangent ∈ (−π/2, π/2).
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
constexpr T atan(T x) {
    if (std::is_constant_evaluated()) {
        return detail::atan(x);
    }
    return MathBackend<T>::atan(x);
}

/**
 * @brief Arcsine ∈ [−π/2, π/2]. Input is clamped to [−1, 1] in both paths so
 *        behavior matches at compile and run time (std::asin would return NaN
 *        for |x| > 1).
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
constexpr T asin(T x) {
    constexpr T half_pi = damp::numbers::pi_v<T> / T{2};
    if (x >= T{1}) {
        return half_pi;
    }
    if (x <= T{-1}) {
        return -half_pi;
    }
    if (std::is_constant_evaluated()) {
        return detail::asin(x);
    }
    return MathBackend<T>::asin(x);
}

/**
 * @brief Arccosine ∈ [0, π]. Input is clamped to [−1, 1] in both paths.
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
constexpr T acos(T x) {
    if (x >= T{1}) {
        return T{0};
    }
    if (x <= T{-1}) {
        return damp::numbers::pi_v<T>;
    }
    if (std::is_constant_evaluated()) {
        return detail::acos(x);
    }
    return MathBackend<T>::acos(x);
}

/**
 * @brief Cosine.
 * @tparam T Numeric type (floating-point)
 * @param x Angle in radians
 */
template<typename T>
constexpr T cos(T x) {
    if (std::is_constant_evaluated()) {
        return detail::cos(x);
    }
    return MathBackend<T>::cos(x);
}

/**
 * @brief Sine.
 * @tparam T Numeric type (floating-point)
 * @param x Angle in radians
 */
template<typename T>
constexpr T sin(T x) {
    if (std::is_constant_evaluated()) {
        return detail::sin(x);
    }
    return MathBackend<T>::sin(x);
}

/**
 * @brief Combined sine and cosine, {sin(x), cos(x)}.
 *
 * At runtime a platform MathBackend can compute both from a single range
 * reduction (~half the cost of separate calls), which is why FOC transforms
 * (Park) should prefer this over calling sin and cos.
 *
 * @tparam T Numeric type (floating-point)
 * @param x Angle in radians
 */
template<typename T>
constexpr damp::pair<T, T> sincos(T x) {
    if (std::is_constant_evaluated()) {
        const auto sc = detail::sincos(x);
        return {sc.sin, sc.cos};
    }
    return MathBackend<T>::sincos(x);
}

/**
 * @brief Tangent.
 * @note Compare with MATLAB®'s tan(x).
 * @tparam T Numeric type (floating-point)
 * @param x Angle in radians
 */
template<typename T>
constexpr T tan(T x) {
    if (std::is_constant_evaluated()) {
        return detail::tan(x);
    }
    return MathBackend<T>::tan(x);
}

/**
 * @brief Exponential function.
 *
 * Compile time: ln2 argument reduction + Taylor series, with over/underflow
 * saturated to max() / 0 (no ±inf — the library builds under -ffinite-math-only).
 * Runtime: MathBackend\<T\>::exp.
 *
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
constexpr T exp(T x) {
    if (std::is_constant_evaluated()) {
        return detail::exp(x);
    }
    return MathBackend<T>::exp(x);
}

/**
 * @brief Natural logarithm.
 *
 * Domain x ≤ 0 is guarded before the compile/runtime dispatch so both paths
 * return 0 (NaN/−inf unavailable under -ffinite-math-only / constant evaluation).
 * For x > 0: compile-time Newton-Raphson, runtime MathBackend\<T\>::log.
 *
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
constexpr T log(T x) {
    if (x <= T{0}) {
        return T{0};
    }
    if (std::is_constant_evaluated()) {
        return detail::log(x);
    }
    return MathBackend<T>::log(x);
}

/**
 * @brief Power function, base^exponent.
 *
 * Domain guards (before dispatch, both paths): exponent == 0 → 1; base ≤ 0 → 0.
 * For base > 0: compile-time exp(exponent · ln(base)), runtime MathBackend\<T\>::pow.
 *
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
constexpr T pow(T base, T exponent) {
    if (exponent == T{0}) {
        return T{1};
    }
    if (base <= T{0}) {
        return T{0};
    }
    if (std::is_constant_evaluated()) {
        return detail::pow(base, exponent);
    }
    return MathBackend<T>::pow(base, exponent);
}

/**
 * @brief Integer power, base^up, via binary exponentiation.
 *
 * Exact and backend-independent (the same implementation at compile and run
 * time), so it is not routed through MathBackend. base == 0 returns 0 for any
 * exponent (including negative — no ±inf under -ffinite-math-only).
 *
 * @param base Base value
 * @param up   Integer exponent
 */
template<typename T>
constexpr T pow(T base, int up) {
    if (up == 0) {
        return T{1};
    }
    if (base == T{0}) {
        return T{0};
    }
    T   result = T{1};
    T   b = base;
    int exponent = up >= 0 ? up : -up;
    while (exponent > 0) {
        if (exponent % 2 == 1) {
            result *= b;
        }
        b *= b;
        exponent /= 2;
    }
    return up >= 0 ? result : T{1} / result;
}

/**
 * @brief Floor — largest integer ≤ x.
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
constexpr T floor(T x) {
    if (std::is_constant_evaluated()) {
        return detail::floor(x);
    }
    return MathBackend<T>::floor(x);
}

/**
 * @brief Ceiling — smallest integer ≥ x.
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
constexpr T ceil(T x) {
    if (std::is_constant_evaluated()) {
        return detail::ceil(x);
    }
    return MathBackend<T>::ceil(x);
}

/**
 * @brief Round to nearest integer; ties to even (IEEE default / `FE_TONEAREST`).
 *
 * Compile-time path and freestanding series backend use the same policy as the
 * hosted runtime backends (`std::nearbyint`, `fast_nearbyint`).
 *
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
constexpr T nearbyint(T x) {
    if (std::is_constant_evaluated()) {
        return detail::nearbyint(x);
    }
    return MathBackend<T>::nearbyint(x);
}

/**
 * @brief Floating-point remainder, x − y·trunc(x/y) (sign of x), matching
 *        std::fmod's truncated-quotient convention.
 *
 * y == 0 is guarded in both paths and returns 0 rather than std::fmod's NaN.
 *
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
constexpr T fmod(T x, T y) {
    if (y == T{0}) {
        return T{0};
    }
    if (std::is_constant_evaluated()) {
        return detail::fmod(x, y);
    }
    return MathBackend<T>::fmod(x, y);
}

/**
 * @brief Base-10 logarithm, log10(x) = ln(x) / ln(10).
 *
 * Domain x ≤ 0 is guarded before the compile/runtime dispatch so both paths
 * return 0. For x > 0: compile-time ln(x)/ln(10), runtime MathBackend\<T\>::log10.
 *
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
constexpr T log10(T x) {
    if (x <= T{0}) {
        return T{0};
    }
    if (std::is_constant_evaluated()) {
        return detail::log10(x);
    }
    return MathBackend<T>::log10(x);
}

/**
 * @brief Sign function — −1 if val < 0, 1 if val > 0, 0 if val == 0.
 *
 * Backend-independent (pure comparisons), so it is not routed through
 * MathBackend.
 */
template<typename T>
constexpr int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

/**
 * @brief Copy sign — magnitude of @p mag with the sign of @p sgn_src.
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
constexpr T copysign(T mag, T sgn_src) {
    if (std::is_constant_evaluated()) {
        return detail::copysign(mag, sgn_src);
    }
    return MathBackend<T>::copysign(mag, sgn_src);
}

/**
 * @brief Finiteness test — false for NaN and ±∞.
 *
 * Always uses IEEE-754 exponent bit inspection (`detail::isfinite`), not
 * `std::isfinite` and not a pure FP comparison. Correct at compile time and at
 * runtime even under `-ffast-math` / `-ffinite-math-only` (where libm
 * `isfinite` is allowed to be a no-op). Divergence guards that only care about
 * large finite values should still use a magnitude bound.
 *
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
constexpr bool isfinite(T x) {
    return detail::isfinite(x);
}

/**
 * @brief True if @p x is finite and strictly greater than zero.
 *
 * Validation helper for design configs (cutoff frequencies, durations,
 * amplitudes, sample periods, …).
 *
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
[[nodiscard]] constexpr bool finite_positive(T x) {
    return damp::isfinite(x) && (x > T{0});
}

/**
 * @brief True if @p x is finite and non-negative (x ≥ 0).
 *
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
[[nodiscard]] constexpr bool finite_non_negative(T x) {
    return damp::isfinite(x) && (x >= T{0});
}

/**
 * @brief Magnitude to decibels, 20·log10(mag).
 * @note Compare with MATLAB®'s mag2db(mag).
 * @tparam T Numeric type (floating-point)
 */
template<typename T = double>
    requires std::is_floating_point_v<T>
constexpr T mag2db(T mag) {
    return T{20} * damp::log10(mag);
}

/**
 * @brief Decibels to magnitude, 10^(db/20).
 * @note Compare with MATLAB®'s db2mag(db).
 * @tparam T Numeric type (floating-point)
 */
template<typename T = double>
    requires std::is_floating_point_v<T>
constexpr T db2mag(T db) {
    return damp::pow(T{10}, db / T{20});
}

/**
 * @brief Radians to degrees, rad·180/π.
 * @note Compare with MATLAB®'s rad2deg(rad).
 * @tparam T Numeric type (floating-point)
 */
template<typename T = double>
    requires std::is_floating_point_v<T>
constexpr T rad2deg(T rad) {
    return rad * (T{180} / damp::numbers::pi_v<T>);
}

/**
 * @brief Degrees to radians, deg·π/180.
 * @note Compare with MATLAB®'s deg2rad(deg).
 * @tparam T Numeric type (floating-point)
 */
template<typename T = double>
    requires std::is_floating_point_v<T>
constexpr T deg2rad(T deg) {
    return deg * (damp::numbers::pi_v<T> / T{180});
}

/**
 * @brief Wrap @p x into the half-open interval [min, max) (period max − min).
 *
 * Round-to-nearest reduction about the interval midpoint: `x − range·nearbyint((x
 * − mid)/range)`. For the common angle case `wrap(θ, −π, π)` this is the cheap
 * phase-reduction kernel — when the bounds are compile-time constants (e.g. ±π)
 * the inlined form folds to one multiply by a constant reciprocal plus one
 * round-to-nearest, with no runtime divide.
 *
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
    requires std::is_floating_point_v<T>
constexpr T wrap(T x, T min, T max) {
    const T range = max - min;
    const T midpoint = min + (range / T{2});

    // Explicit \<T\>: force the constexpr dispatcher here. A bare damp::nearbyint(float) would bind
    // to the non-template fast-float overload in trig.hpp (when visible) and break constexpr eval;
    // the dispatcher still reaches that fast path at runtime via MathBackend<float>.
    T y = x - (range * damp::nearbyint<T>((x - midpoint) / range));

    // Rounding ties or floating error can land on max or just outside [min, max).
    // Fold once more so the result is half-open.
    if (y < min) {
        y += range;
    }
    if (y >= max) {
        y -= range;
    }
    return y;
}

/**
 * @brief Wrap an angle to [−π, π)
 *
 * Shortest-arc residual / signed phase error: heading innovation, PLL angle
 * error, joint nearest-branch, motor electrical angle when a bipolar wrap is
 * preferred. For a *stored* PLL phase that locks near 0, prefer @ref wrap_two_pi
 * so the discontinuity sits at a full turn rather than at ±π (180°).
 *
 * @note Compare with MATLAB®'s wrapToPi(θ).
 * @tparam T Numeric type (floating-point)
 * @param theta Angle [rad]
 * @return θ folded into [−π, π)
 */
template<typename T = double>
    requires std::is_floating_point_v<T>
[[nodiscard]] constexpr T wrap_pi(T theta) {
    return wrap(theta, -damp::numbers::pi_v<T>, damp::numbers::pi_v<T>);
}

/**
 * @brief Wrap an angle to [0, 2π)
 *
 * Natural range for a PLL / VCO phase integrator that locks about 0: noise near
 * lock does not flicker across a wrap seam, and split-phase 180° is interior to
 * the open interval of @ref wrap_pi.
 *
 * @note Compare with MATLAB®'s wrapTo2Pi(θ).
 * @tparam T Numeric type (floating-point)
 * @param theta Angle [rad]
 * @return θ folded into [0, 2π)
 */
template<typename T = double>
    requires std::is_floating_point_v<T>
[[nodiscard]] constexpr T wrap_two_pi(T theta) {
    return wrap(theta, T{0}, T{2} * damp::numbers::pi_v<T>);
}

} // namespace damp
