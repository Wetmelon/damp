// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file constexpr_math.hpp
 * @brief Compile-time scalar math for the damp:: dispatch layer.
 *
 * Bodies used when `std::is_constant_evaluated()` is true (and by the freestanding
 * `SeriesMathBackend`): Newton–Raphson (`sqrt`, `cbrt` after exact power-of-two
 * scaling), atanh series `log` after exact power-of-two scaling, Taylor kernels on |r| ≤ π/4 after π/2 quadrant reduction
 * (`sin`/`cos`/`sincos`/`tan`), Taylor `exp`, fdlibm-style multi-interval atan
 * (asin/acos via atan2), Cody–Waite two-part argument reduction,
 * cast-based `floor`/`ceil`/`fmod`, and a bit-pattern `isfinite`.
 *
 * No damp backend coupling: this header includes only freestanding standard
 * headers (`<bit>`, `<cstdint>`, `<limits>`, `<type_traits>`). It does not pull
 * `damp/backend.hpp` (std/ETL vocabulary) or `MathBackend` / libm. π and the
 * `sincos` result type are local; public wrappers map them to `damp::numbers` /
 * `damp::pair` in math.hpp / series_backend.hpp.
 *
 * Everything here lives in namespace damp::detail and is not part of the public
 * API; call the dispatching wrappers in math.hpp instead. Domain policies
 * (non-positive → 0, `|x|>1` clamp for inverse trig, `fmod` by zero → 0) are
 * self-contained here so freestanding callers match the public dispatcher.
 *
 * @see math.hpp for the public dispatch layer
 */

#include <bit>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace damp::detail {

/// Local π (same digits as `damp::numbers::pi_v`; kept here to avoid backend.hpp).
template<typename T>
inline constexpr T pi_v = static_cast<T>(3.141592653589793238462643383279502884L);

/// `{sin, cos}` result — plain aggregate so this header never needs `damp::pair`.
template<typename T>
struct sincos_result {
    T sin;
    T cos;
};

// ---------------------------------------------------------------------------
// Inverse-trig kernel — fdlibm / SunPro atan (public domain notice preserved
// in the coefficient block). Multi-interval reduction + odd/even split poly,
// ~1 ULP on double; evaluated entirely in constexpr without integer bit casts.
// ---------------------------------------------------------------------------

// Coefficients from FreeBSD msun s_atan.c (SunPro). Decimal values round-trip
// to the intended IEEE doubles.
inline constexpr double atan_aT[] = {
    +3.33333333333329318027e-01,
    -1.99999999998764832476e-01,
    +1.42857142725034663711e-01,
    -1.11111104054623557880e-01,
    +9.09088713343650656196e-02,
    -7.69187620504482999495e-02,
    +6.66107313738753120669e-02,
    -5.83357013379057348645e-02,
    +4.97687799461593236017e-02,
    -3.65315727442169155270e-02,
    +1.62858201153657823623e-02,
};

inline constexpr double atan_hi[] = {
    +4.63647609000806093515e-01, // atan(0.5)
    +7.85398163397448278999e-01, // atan(1)
    +9.82793723247329054082e-01, // atan(1.5)
    +1.57079632679489655800e+00, // atan(inf) = π/2
};

inline constexpr double atan_lo[] = {
    +2.26987774529616870924e-17,
    +3.06161699786838301793e-17,
    +1.39033110312309984516e-17,
    +6.12323399573676603587e-17,
};

/// True when @p x is strictly inside the open range of `long long` after cast to
/// `T` (so `static_cast<long long>(x)` is defined). Outside that range every
/// finite binary `float`/`double` is already an integer (ulp ≥ 1 long before
/// |x| reaches 2⁶³), so floor/ceil/trunc can return @p x unchanged.
template<typename T>
constexpr bool in_long_long_cast_range(T x) {
    // LLONG_MAX is not exactly representable as double (rounds to 2⁶³); LLONG_MIN
    // is exact (−2⁶³). Open bounds avoid a cast at the unrepresentable endpoint.
    constexpr T hi = static_cast<T>(std::numeric_limits<long long>::max());
    constexpr T lo = static_cast<T>(std::numeric_limits<long long>::min());
    return (x > lo) && (x < hi);
}

/**
 * @brief Round to nearest integer (ties away from zero), returned as long long.
 *
 * Used for argument reduction in the constexpr trig/exp paths. Values outside
 * the long-long range are clamped (trig reduction is already meaningless once
 * the quotient loses all fractional bits).
 */
template<typename T>
constexpr long long lround_away(T x) {
    const T y = x >= T{0} ? x + static_cast<T>(0.5) : x - static_cast<T>(0.5);
    if (!in_long_long_cast_range(y)) {
        return y >= T{0} ? std::numeric_limits<long long>::max()
                         : std::numeric_limits<long long>::min();
    }
    return static_cast<long long>(y);
}

/// Keep @p v as computed: stops `-ffast-math` from reassociating a Cody–Waite
/// reduction `(x − k·hi) − k·lo` into `x − k·(hi + lo)` when the series run at
/// runtime (constant evaluation is never reassociated).
template<typename T>
constexpr T assoc_barrier(T v) {
#if defined(__has_builtin)
#if __has_builtin(__builtin_assoc_barrier)
    return __builtin_assoc_barrier(v);
#endif
#endif
    return v;
}

/// |x|.
template<typename T>
constexpr T abs(T x) {
    return x >= T{0} ? x : -x;
}

template<typename T>
constexpr bool isfinite(T x); // bit-pattern test, defined at the end of this header

/// Argument and scale of @ref reduce_for_root.
template<typename T, int Root>
struct RootReduced {
    T x;     ///< Reduced argument in [1, 2^Root)
    T scale; ///< 2^k with root(input) = scale · root(x)
};

/**
 * @brief Exact power-of-two range reduction for the Newton root kernels.
 *
 * Rescales finite @p x > 0 by 2^(Root·k) into [1, 2^Root) and returns the
 * matching 2^k in @c scale, so root(x_in) = scale · root(x). Every multiply is
 * by a power of two (exact, subnormals included), and coarse 2^(±16·Root) steps
 * bound the loop count over the whole float/double exponent range.
 */
template<typename T, int Root>
constexpr RootReduced<T, Root> reduce_for_root(T x) {
    constexpr T big = static_cast<T>(1ULL << (Root * 16)); // 2^(16·Root)
    constexpr T big_root = T{65536};                       // 2^16
    constexpr T step = static_cast<T>(1ULL << Root);       // 2^Root
    T           scale = T{1};
    while (x >= big) {
        x /= big;
        scale *= big_root;
    }
    while (x < T{1} / big) {
        x *= big;
        scale /= big_root;
    }
    while (x >= step) {
        x /= step;
        scale *= T{2};
    }
    while (x < T{1}) {
        x *= step;
        scale /= T{2};
    }
    return {x, scale};
}

/// Square root via Newton–Raphson after exact power-of-four reduction into
/// [1, 4) (converges in a handful of steps for any finite @p x). Domain
/// @p x < 0 → 0 (library policy under `-ffinite-math-only` / constant
/// evaluation — no NaN).
template<typename T>
constexpr T sqrt(T x) {
    if (x <= T{0}) {
        return T{0};
    }
    if (!isfinite(x)) {
        return x; // +inf / NaN pass through (the reduction loops need a finite x)
    }
    const auto red = reduce_for_root<T, 2>(x);
    T          guess = (red.x + T{1}) / T{2}; // ≥ √x: Newton decreases monotonically
    for (int i = 0; i < 50; ++i) {
        const T next = (guess + (red.x / guess)) / T{2};
        if (next >= guess) {
            break;
        }
        guess = next;
    }
    return guess * red.scale;
}

/// Cube root via Newton–Raphson after exact power-of-eight reduction into
/// [1, 8); preserves sign for negative x.
template<typename T>
constexpr T cbrt(T x) {
    if (x == T{0} || !isfinite(x)) {
        return x; // ±0, ±inf, NaN pass through
    }
    const bool neg = x < T{0};
    const auto red = reduce_for_root<T, 3>(neg ? -x : x);
    T          guess = (red.x + T{2}) / T{3}; // ≥ ∛x: Newton decreases monotonically
    for (int i = 0; i < 50; ++i) {
        const T next = ((T{2} * guess) + (red.x / (guess * guess))) / T{3};
        if (next >= guess) {
            break;
        }
        guess = next;
    }
    const T r = guess * red.scale;
    return neg ? -r : r;
}

/**
 * @brief Arctangent on [0, ∞) — fdlibm range reduction + split odd/even poly.
 *
 * Breakpoints at 7/16, 11/16, 19/16, 39/16; kernel
 * @f$ x - x(s_1+s_2) @f$ or @f$ \mathrm{atanhi}_k - ((x(s_1+s_2)-\mathrm{atanlo}_k)-x) @f$.
 */
template<typename T>
constexpr T atan_nonneg(T x) {
    // |x| ≥ 2^66 → ±π/2 (argument is already non-negative here).
    if (x >= static_cast<T>(0x1p66)) {
        return static_cast<T>(atan_hi[3]) + static_cast<T>(atan_lo[3]);
    }

    int id = -1;
    if (x < static_cast<T>(0.4375)) { // [0, 7/16)
        if (x < static_cast<T>(0x1p-27)) {
            return x; // underflow / identity
        }
    } else if (x < static_cast<T>(1.1875)) { // [7/16, 19/16)
        if (x < static_cast<T>(0.6875)) {    // [7/16, 11/16)
            id = 0;
            x = ((T{2} * x) - T{1}) / (T{2} + x);
        } else { // [11/16, 19/16)
            id = 1;
            x = (x - T{1}) / (x + T{1});
        }
    } else if (x < static_cast<T>(2.4375)) { // [19/16, 39/16)
        id = 2;
        x = (x - static_cast<T>(1.5)) / (T{1} + (static_cast<T>(1.5) * x));
    } else { // [39/16, 2^66)
        id = 3;
        x = -T{1} / x;
    }

    const T z = x * x;
    const T w = z * z;
    // Split sum aT[i] z^{i+1} into odd/even parts in w = z² (fdlibm layout).
    T s1 = static_cast<T>(atan_aT[10]);
    s1 = static_cast<T>(atan_aT[8]) + (w * s1);
    s1 = static_cast<T>(atan_aT[6]) + (w * s1);
    s1 = static_cast<T>(atan_aT[4]) + (w * s1);
    s1 = static_cast<T>(atan_aT[2]) + (w * s1);
    s1 = static_cast<T>(atan_aT[0]) + (w * s1);
    s1 *= z;

    T s2 = static_cast<T>(atan_aT[9]);
    s2 = static_cast<T>(atan_aT[7]) + (w * s2);
    s2 = static_cast<T>(atan_aT[5]) + (w * s2);
    s2 = static_cast<T>(atan_aT[3]) + (w * s2);
    s2 = static_cast<T>(atan_aT[1]) + (w * s2);
    s2 *= w;

    if (id < 0) {
        return x - (x * (s1 + s2));
    }
    const T hi = static_cast<T>(atan_hi[id]);
    const T lo = static_cast<T>(atan_lo[id]);
    return hi - (((x * (s1 + s2)) - lo) - x);
}

/// Single-argument arctangent ∈ (−π/2, π/2). Odd; uses @ref atan_nonneg.
template<typename T>
constexpr T atan(T x) {
    const T ax = x >= T{0} ? x : -x;
    const T r = atan_nonneg(ax);
    return x < T{0} ? -r : r;
}

/**
 * @brief Two-argument arctangent ∈ [−π, π].
 *
 * Uses @ref atan_nonneg on the magnitude ratio (argument ≤ 1), then restores
 * quadrant from the signs of @p x and @p y.
 */
template<typename T>
constexpr T atan2(T y, T x) {
    constexpr T pi = pi_v<T>;
    constexpr T half_pi = pi / T{2};

    if (x == T{0}) {
        if (y > T{0}) {
            return half_pi;
        }
        if (y < T{0}) {
            return -half_pi;
        }
        return T{0};
    }

    const T ax = x >= T{0} ? x : -x;
    const T ay = y >= T{0} ? y : -y;
    T       a = (ay > ax) ? (half_pi - atan_nonneg(ax / ay)) : atan_nonneg(ay / ax);
    if (x < T{0}) {
        a = pi - a;
    }
    return y < T{0} ? -a : a;
}

/**
 * @brief Arcsine ∈ [−π/2, π/2]. Clamps |x| > 1.
 *
 * @f$ \arcsin x = \mathrm{atan2}\bigl(x,\sqrt{(1-x)(1+x)}\bigr) @f$ — product form
 * under the root avoids cancellation near |x| = 1.
 */
template<typename T>
constexpr T asin(T x) {
    constexpr T half_pi = pi_v<T> / T{2};
    if (x >= T{1}) {
        return half_pi;
    }
    if (x <= T{-1}) {
        return -half_pi;
    }
    return atan2(x, sqrt((T{1} - x) * (T{1} + x)));
}

/**
 * @brief Arccosine ∈ [0, π]. Clamps |x| > 1.
 *
 * @f$ \arccos x = \mathrm{atan2}\bigl(\sqrt{(1-x)(1+x)}, x\bigr) @f$.
 */
template<typename T>
constexpr T acos(T x) {
    if (x >= T{1}) {
        return T{0};
    }
    if (x <= T{-1}) {
        return pi_v<T>;
    }
    return atan2(sqrt((T{1} - x) * (T{1} + x)), x);
}

/**
 * @brief Quadrant reduction: x = k·π/2 + r with |r| ≲ π/4 (Cody–Waite, two-part π/2).
 *
 * The high word has its low mantissa bits cleared so k·hi is exact for
 * |k| < 2²⁰; accuracy degrades gracefully beyond |x| ~ 10⁶ rad.
 */
template<typename T>
struct QuadrantReduced {
    T   r;        ///< Remainder in about [−π/4, π/4]
    int quadrant; ///< k mod 4, in {0, 1, 2, 3}
};

template<typename T>
constexpr QuadrantReduced<T> reduce_half_pi(T x) {
    constexpr T     half_pi_hi = static_cast<T>(1.57079632673412561417);
    constexpr T     half_pi_lo = static_cast<T>(6.07710050650619224932e-11);
    constexpr T     two_over_pi = static_cast<T>(0.63661977236758134308);
    const long long k = lround_away(x * two_over_pi);
    const T         kreal = static_cast<T>(k);
    return {assoc_barrier(x - (kreal * half_pi_hi)) - (kreal * half_pi_lo), static_cast<int>(k & 3LL)};
}

/// Taylor sine on |r| ≤ π/4 (truncation < 10⁻²² there).
template<typename T>
constexpr T sin_kernel(T r) {
    const T r2 = r * r;
    T       result = r;
    T       term = r;
    for (int n = 1; n <= 10; ++n) {
        term *= -r2 / T((2 * n) * ((2 * n) + 1));
        result += term;
    }
    return result;
}

/// Taylor cosine on |r| ≤ π/4 (truncation < 10⁻²² there).
template<typename T>
constexpr T cos_kernel(T r) {
    const T r2 = r * r;
    T       result = T{1};
    T       term = T{1};
    for (int n = 1; n <= 10; ++n) {
        term *= -r2 / T((2 * n) * ((2 * n) - 1));
        result += term;
    }
    return result;
}

/// {sin(x), cos(x)} from one quadrant reduction and the |r| ≤ π/4 kernels.
template<typename T>
constexpr sincos_result<T> sincos(T x) {
    const auto red = reduce_half_pi(x);
    const T    s = sin_kernel(red.r);
    const T    c = cos_kernel(red.r);
    switch (red.quadrant) {
        case 0:
            return {s, c};
        case 1:
            return {c, -s};
        case 2:
            return {-s, -c};
        default:
            return {-c, s};
    }
}

/// Sine (quadrant reduction + Taylor kernel; ~1 ulp absolute).
template<typename T>
constexpr T sin(T x) {
    return sincos(x).sin;
}

/// Cosine (quadrant reduction + Taylor kernel; ~1 ulp absolute).
template<typename T>
constexpr T cos(T x) {
    return sincos(x).cos;
}

/// Tangent as a kernel ratio after quadrant reduction: sin r / cos r for even
/// k, −cos r / sin r for odd k. The reduced r keeps relative accuracy next to
/// the poles, so tan stays relatively accurate there.
template<typename T>
constexpr T tan(T x) {
    const auto red = reduce_half_pi(x);
    const T    s = sin_kernel(red.r);
    const T    c = cos_kernel(red.r);
    return ((red.quadrant & 1) == 0) ? (s / c) : (-c / s);
}

/// exp via ln2 argument reduction (@f$ e^x = 2^k e^r @f$, |r| ≤ ln2/2) and a
/// Taylor series on the remainder. Over/underflow saturates to `max()` / `0`
/// (no ±inf under `-ffinite-math-only`); thresholds are per-`T` (`float` vs
/// `double` range).
template<typename T>
constexpr T exp(T x) {
    // log(numeric_limits<T>::max()) and ~log(min subnormal), type-scaled.
    constexpr T overflow_x = std::is_same_v<T, float> ? static_cast<T>(88.722839f)
                                                      : static_cast<T>(709.782712893384);
    constexpr T underflow_x = std::is_same_v<T, float> ? static_cast<T>(-103.278929f)
                                                       : static_cast<T>(-745.133219101941);
    if (x > overflow_x) {
        return std::numeric_limits<T>::max();
    }
    if (x < underflow_x) {
        return T{0};
    }
    if (!isfinite(x)) {
        return x; // NaN (±inf were saturated above)
    }

    constexpr T     ln2_hi = static_cast<T>(6.93147180369123816490e-01);
    constexpr T     ln2_lo = static_cast<T>(1.90821492927058770002e-10);
    constexpr T     inv_ln2 = static_cast<T>(1.44269504088896340736);
    const long long k = lround_away(x * inv_ln2);
    const T         kreal = static_cast<T>(k);
    const T         r = assoc_barrier(x - (kreal * ln2_hi)) - (kreal * ln2_lo);

    T result = T{1};
    T term = T{1};
    for (int n = 1; n <= 14; ++n) {
        term *= r / T(n);
        result += term;
    }

    // Scale by 2^k via repeated doubling/halving on the O(1)-magnitude series
    // result so intermediates stay finite — forming 2^k directly would overflow
    // to ∞ on the final squaring for large k, which is not a constant expression
    // even when the final exp(x) is finite.
    if (k >= 0) {
        for (long long i = 0; i < k; ++i) {
            result *= T{2};
        }
    } else {
        for (long long i = 0; i < -k; ++i) {
            result *= static_cast<T>(0.5);
        }
    }
    return result;
}

/**
 * @brief Natural log: exact power-of-two reduction, then an atanh series.
 *
 * @f$ x = 2^m y @f$ with @f$ y \in [1/\sqrt2, \sqrt2) @f$ (exact scaling, so
 * subnormals and 1e±300 cost a few steps and lose nothing), then
 * @f$ \ln y = 2\,\mathrm{atanh}\,s @f$, @f$ s = (y-1)/(y+1) @f$, |s| < 0.172,
 * and @f$ \ln x = m \ln 2 + \ln y @f$ with a two-part ln 2.
 * Domain @p x ≤ 0 → 0; +inf / NaN pass through.
 */
template<typename T>
constexpr T log(T x) {
    if (x <= T{0}) {
        return T{0};
    }
    if (!isfinite(x)) {
        return x;
    }
    constexpr T big = static_cast<T>(1ULL << 32);
    constexpr T sqrt2 = static_cast<T>(1.41421356237309504880);
    T           m = T{0};
    while (x >= big) {
        x /= big;
        m += T{32};
    }
    while (x < T{1} / big) {
        x *= big;
        m -= T{32};
    }
    while (x >= sqrt2) {
        x /= T{2};
        m += T{1};
    }
    while (x < sqrt2 / T{2}) {
        x *= T{2};
        m -= T{1};
    }

    const T s = (x - T{1}) / (x + T{1});
    const T s2 = s * s;
    T       term = s;
    T       sum = s;
    for (int n = 1; n <= 12; ++n) { // s^25/25 < 1e-20
        term *= s2;
        sum += term / T((2 * n) + 1);
    }
    constexpr T ln2_hi = static_cast<T>(6.93147180369123816490e-01);
    constexpr T ln2_lo = static_cast<T>(1.90821492927058770002e-10);
    return (m * ln2_hi) + ((m * ln2_lo) + (T{2} * sum));
}

/// @f$ \mathrm{base}^{\mathrm{exponent}} = \exp(\mathrm{exponent}\cdot\ln\mathrm{base}) @f$.
/// `exponent == 0` → 1 (including 0⁰); `base ≤ 0` (and nonzero exponent) → 0.
template<typename T>
constexpr T pow(T base, T exponent) {
    if (exponent == T{0}) {
        return T{1};
    }
    if (base <= T{0}) {
        return T{0};
    }
    return exp(log(base) * exponent);
}

/// Largest integer ≤ @p x (full finite range; no long-long UB on huge |x|).
template<typename T>
constexpr T floor(T x) {
    if (!in_long_long_cast_range(x)) {
        // |x| ≥ 2⁶³ (or so): every finite binary float/double is already integral.
        return x;
    }
    T int_part = static_cast<T>(static_cast<long long>(x)); // toward zero
    if (x < T{0} && x != int_part) {
        int_part -= T{1};
    }
    return int_part;
}

/// Smallest integer ≥ @p x (full finite range).
template<typename T>
constexpr T ceil(T x) {
    return -floor(-x);
}

/**
 * @brief Round to nearest integer; ties to even (IEEE default / `std::nearbyint`
 *        under `FE_TONEAREST`).
 *
 * Matches the runtime backends (`fast_nearbyint`, libm). Full finite range: large
 * |x| are already integers so the value is returned unchanged.
 */
template<typename T>
constexpr T nearbyint(T x) {
    const T fl = floor(x);
    const T frac = x - fl;
    if (frac < static_cast<T>(0.5)) {
        return fl;
    }
    if (frac > static_cast<T>(0.5)) {
        return fl + T{1};
    }
    // Half-integer: choose the even integer of {fl, fl+1}.
    // |fl| is always well below 2⁶³ here (half-integers do not exist once ulp ≥ 1).
    const auto n = static_cast<long long>(fl);
    return ((n % 2) == 0) ? fl : (fl + T{1});
}

/// Non-negative remainder @p x mod @p y for @p x ≥ 0, @p y > 0, result in [0, y).
/// Binary (shift-and-subtract) reduction — exact for large |x/y| where
/// `x - y·trunc(x/y)` loses all precision in the product.
template<typename T>
constexpr T fmod_positive(T x, T y) {
    if (x < y) {
        return x;
    }
    T   d = y;
    int shifts = 0;
    // Cover the full binary64 exponent span with margin; each step doubles d.
    constexpr int max_shifts = 4096;
    while (shifts < max_shifts && d < (x * static_cast<T>(0.5)) && d < (std::numeric_limits<T>::max() * static_cast<T>(0.5))) {
        d *= T{2};
        ++shifts;
    }
    for (int i = 0; i <= shifts; ++i) {
        if (x >= d) {
            x -= d;
        }
        d *= static_cast<T>(0.5);
    }
    return x;
}

/// Floating-point remainder, @f$ x - y\cdot\mathrm{trunc}(x/y) @f$ (sign of x),
/// matching `std::fmod`'s truncated-quotient convention. @p y == 0 → 0.
/// Full finite range (no long-long UB; no catastrophic cancellation for huge x/y).
template<typename T>
constexpr T fmod(T x, T y) {
    if (y == T{0}) {
        return T{0};
    }
    const T ay = y >= T{0} ? y : -y;
    const T ax = x >= T{0} ? x : -x;
    const T r = fmod_positive(ax, ay);
    return x >= T{0} ? r : -r;
}

/// @f$ \log_{10}(x) = \ln(x)/\ln(10) @f$. Domain @p x ≤ 0 → 0.
template<typename T>
constexpr T log10(T x) {
    if (x <= T{0}) {
        return T{0};
    }
    return log(x) / log(T{10});
}

/// Magnitude of @p mag with the sign of @p sgn_src from a comparison
/// (`sgn_src < 0`), not IEEE `signbit` (so −0 is treated as non-negative).
template<typename T>
constexpr T copysign(T mag, T sgn_src) {
    const T m = mag >= T{0} ? mag : -mag;
    return sgn_src < T{0} ? -m : m;
}

/**
 * @brief Finiteness test — false for NaN and ±∞ (IEEE-754 bit patterns)
 *
 * Inspects the exponent field via `std::bit_cast` so the result stays correct
 * under `-ffast-math` / `-ffinite-math-only`. Those flags may make
 * `std::isfinite` and FP comparisons (`|x| ≤ max`) optimize to "always true"
 * because the compiler assumes no NaNs/Infs — integer bit tests do not.
 *
 * @note Prefer magnitude guards (`abs(x) > bound`) for *divergence* detection
 *       when values grow without necessarily becoming Inf.
 */
template<typename T>
constexpr bool isfinite(T x) {
    if constexpr (std::is_same_v<T, float>) {
        const auto bits = std::bit_cast<std::uint32_t>(x);
        return (bits & 0x7f800000u) != 0x7f800000u;
    } else if constexpr (std::is_same_v<T, double>) {
        const auto bits = std::bit_cast<std::uint64_t>(x);
        return (bits & 0x7ff0000000000000ull) != 0x7ff0000000000000ull;
    } else {
        // long double / extended: fall back to comparison (best-effort)
        return abs(x) <= std::numeric_limits<T>::max();
    }
}

} // namespace damp::detail
