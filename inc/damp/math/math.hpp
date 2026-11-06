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
    if (damp::is_constant_evaluated()) {
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
    if (damp::is_constant_evaluated()) {
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
    return hi * damp::sqrt(T{1} + (r * r));
}

/**
 * @brief Cube root (preserves sign for negative x).
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
constexpr T cbrt(T x) {
    if (damp::is_constant_evaluated()) {
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
    if (damp::is_constant_evaluated()) {
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
    if (damp::is_constant_evaluated()) {
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
    if (damp::is_constant_evaluated()) {
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
    if (damp::is_constant_evaluated()) {
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
    if (damp::is_constant_evaluated()) {
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
    if (damp::is_constant_evaluated()) {
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
    if (damp::is_constant_evaluated()) {
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
    if (damp::is_constant_evaluated()) {
        return detail::tan(x);
    }
    return MathBackend<T>::tan(x);
}

/**
 * @brief Exponential function.
 *
 * Compile time: ln2 argument reduction + Taylor series, with over/underflow
 * saturated to max() / 0 (no ±inf — the library builds under -ffinite-math-only).
 * Runtime: MathBackend\<T\>::exp, with overflow saturated to max() the same way.
 *
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
constexpr T exp(T x) {
    if (damp::is_constant_evaluated()) {
        return detail::exp(x);
    }
    const T r = MathBackend<T>::exp(x);
    return detail::isfinite(r) ? r : damp::numeric_limits<T>::max();
}

/**
 * @brief Natural logarithm.
 *
 * Domain x ≤ 0 is guarded before the compile/runtime dispatch so both paths
 * return 0 (NaN/−inf unavailable under -ffinite-math-only / constant evaluation).
 * For x > 0: compile-time atanh series, runtime MathBackend\<T\>::log.
 *
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
constexpr T log(T x) {
    if (x <= T{0}) {
        return T{0};
    }
    if (damp::is_constant_evaluated()) {
        return detail::log(x);
    }
    return MathBackend<T>::log(x);
}

template<typename T>
constexpr T floor(T x);

/**
 * @brief Power function with a real exponent, base^exponent.
 *
 * Domain guards (before dispatch, both paths): exponent == 0 → 1; base == 0, or
 * base < 0 with a non-integral exponent, → 0. A negative base with an integral
 * exponent keeps its sign, as std::pow. Otherwise compile-time
 * exp(exponent · ln|base|), runtime MathBackend\<T\>::pow, with overflow saturated
 * to ±max().
 *
 * The exponent takes the type of @p base (it does not take part in deduction), so
 * a float base with a double exponent stays on this overload. Pass an integer
 * type for an exact integer power (the binary-exponentiation overload).
 *
 * @tparam T Floating-point type
 */
template<typename T>
    requires damp::is_floating_point_v<T>
constexpr T pow(T base, damp::type_identity_t<T> exponent) {
    if (exponent == T{0}) {
        return T{1};
    }
    const bool signed_integral = (base < T{0}) && (exponent == damp::floor(exponent));
    if ((base <= T{0}) && !signed_integral) {
        return T{0};
    }
    if (damp::is_constant_evaluated()) {
        return detail::pow(base, exponent);
    }
    const T r = MathBackend<T>::pow(base, exponent);
    if (detail::isfinite(r)) {
        return r;
    }
    return (r < T{0}) ? -damp::numeric_limits<T>::max() : damp::numeric_limits<T>::max();
}

/**
 * @brief Integer power, base^up, via binary exponentiation.
 *
 * Exact and backend-independent (the same implementation at compile and run
 * time), so it is not routed through MathBackend. base == 0 returns 0 for any
 * exponent (including negative — no ±inf under -ffinite-math-only). A result
 * past the finite range saturates to ±max(), and a negative power whose
 * magnitude would exceed max() returns 0; no intermediate overflows, so this
 * holds in constant evaluation too.
 *
 * @tparam T Floating-point type
 * @tparam I Integer type of the exponent
 * @param base Base value
 * @param up   Integer exponent
 */
template<typename T, damp::integral I>
    requires damp::is_floating_point_v<T>
constexpr T pow(T base, I up) {
    if (up == I{0}) {
        return T{1};
    }
    if (base == T{0}) {
        return T{0};
    }
    using U = damp::make_unsigned_t<I>;
    constexpr T max = damp::numeric_limits<T>::max();
    // Factors at or below lim cannot overflow a product (lim² ≤ max/4); only past it
    // is the exact max/a test needed.
    constexpr T lim = detail::sqrt(max) / T{2};
    // |up| in unsigned arithmetic: −min() overflows I.
    U          e = (up >= I{0}) ? static_cast<U>(up) : static_cast<U>(U{0} - static_cast<U>(up));
    const bool negative_result = (base < T{0}) && ((e & U{1}) != U{0});
    T          b = detail::abs(base);
    T          p = T{1};
    bool       saturated = false;
    while (true) {
        if ((e & U{1}) != U{0}) {
            if (((p > lim) || (b > lim)) && (p > max / b)) {
                saturated = true;
                break;
            }
            p *= b;
        }
        e >>= 1U;
        if (e == U{0}) {
            break;
        }
        if ((b > lim) && (b > max / b)) {
            saturated = true; // a later set bit multiplies by b² > max (and p ≥ 1 here)
            break;
        }
        b *= b;
    }
    T out{};
    if (up > I{0}) {
        out = saturated ? max : p;
    } else if (saturated) {
        out = T{0};
    } else {
        // 1/p would round past max: covers p == 0 and a subnormal p (kept without FTZ and
        // in constant evaluation). `<=`, since 1/fl(1/max) itself rounds past max.
        out = (p <= T{1} / max) ? max : T{1} / p;
    }
    return negative_result ? -out : out;
}

/**
 * @brief Floor — largest integer ≤ x.
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
constexpr T floor(T x) {
    if (damp::is_constant_evaluated()) {
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
    if (damp::is_constant_evaluated()) {
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
    if (damp::is_constant_evaluated()) {
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
    if (damp::is_constant_evaluated()) {
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
    if (damp::is_constant_evaluated()) {
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
 * @brief Smoothed sign: an odd quintic through ±@p width, exactly ±1 outside it.
 *
 * @f$ s(u) = \tfrac{u}{8}(15 - 10u^2 + 3u^4) @f$ on @f$ u = x/w \in [-1, 1] @f$: value ±1 and
 * zero first and second derivatives at the band edges (C²), slope 15/(8w) through zero.
 * Unlike tanh or x/(|x|+ε) it reaches ±1 exactly, so outside the band it is the true sign,
 * and it costs a few multiplies. @p width ≤ 0 gives the plain sign.
 */
template<typename T>
[[nodiscard]] constexpr T smooth_sign(T x, T width) {
    if (!(width > T{0})) {
        return static_cast<T>(sgn(x));
    }
    if (x >= width) {
        return T{1};
    }
    if (x <= -width) {
        return T{-1};
    }
    const T u = x / width;
    const T u2 = u * u;
    return u * (T{15} - (T{10} * u2) + (T{3} * u2 * u2)) / T{8};
}

/**
 * @brief Copy sign — magnitude of @p mag with the sign of @p sgn_src.
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
constexpr T copysign(T mag, T sgn_src) {
    if (damp::is_constant_evaluated()) {
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
 * @return 20·log10(mag); −max() (the finite stand-in for −∞) for mag ≤ 0
 */
template<typename T = double>
    requires damp::is_floating_point_v<T>
constexpr T mag2db(T mag) {
    if (!(mag > T{0})) {
        return -damp::numeric_limits<T>::max();
    }
    return T{20} * damp::log10(mag);
}

/**
 * @brief Decibels to magnitude, 10^(db/20).
 * @note Compare with MATLAB®'s db2mag(db).
 * @tparam T Numeric type (floating-point)
 */
template<typename T = double>
    requires damp::is_floating_point_v<T>
constexpr T db2mag(T db) {
    return damp::pow(T{10}, db / T{20});
}

/**
 * @brief Radians to degrees, rad·180/π.
 * @note Compare with MATLAB®'s rad2deg(rad).
 * @tparam T Numeric type (floating-point)
 */
template<typename T = double>
    requires damp::is_floating_point_v<T>
constexpr T rad2deg(T rad) {
    return rad * (T{180} / damp::numbers::pi_v<T>);
}

/**
 * @brief Degrees to radians, deg·π/180.
 * @note Compare with MATLAB®'s deg2rad(deg).
 * @tparam T Numeric type (floating-point)
 */
template<typename T = double>
    requires damp::is_floating_point_v<T>
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
 * Once ulp(x) approaches the period (|x| ≳ range/ε) the product above no longer
 * holds the remainder; the result then falls back to the exact remainder of x
 * and stays in [min, max), but x itself no longer resolves the angle. An empty
 * interval (max ≤ min) returns min.
 *
 * @tparam T Numeric type (floating-point)
 */
template<typename T>
    requires damp::is_floating_point_v<T>
constexpr T wrap(T x, T min, T max) noexcept {
    const T range = max - min;
    if (!(range > T{0})) {
        return min;
    }
    const T midpoint = min + (range / T{2});
    T       y = x - (range * damp::nearbyint((x - midpoint) / range));

    // Rounding ties or floating error can land on max or just outside [min, max).
    // Fold once more so the result is half-open.
    if (y < min) {
        y += range;
    }
    if (y >= max) {
        y -= range;
    }
    if ((y < min) || (y >= max)) { // |x| ≫ range: fall back to the exact remainder
        y = min + damp::fmod(x - min, range);
        if (y < min) {
            y += range;
        }
        if (y >= max) {
            y = min;
        }
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
    requires damp::is_floating_point_v<T>
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
    requires damp::is_floating_point_v<T>
[[nodiscard]] constexpr T wrap_two_pi(T theta) {
    return wrap(theta, T{0}, T{2} * damp::numbers::pi_v<T>);
}

} // namespace damp
