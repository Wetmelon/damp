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
 * `SeriesMathBackend`): Newton–Raphson (`sqrt`, `cbrt`, `log`), Taylor series
 * (`sin`/`cos`/`exp`), fdlibm-style multi-interval atan (asin/acos via atan2),
 * continued fraction (`tan`), Cody–Waite two-part argument reduction,
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

/// |x|.
template<typename T>
constexpr T abs(T x) {
    return x >= T{0} ? x : -x;
}

/// Square root via Newton–Raphson. Domain @p x < 0 → 0 (library policy under
/// `-ffinite-math-only` / constant evaluation — no NaN).
template<typename T>
constexpr T sqrt(T x) {
    if (x <= T{0}) {
        return T{0};
    }
    T guess = x > T{1} ? x / T{2} : T{1};
    for (int i = 0; i < 50; ++i) {
        T next = (guess + (x / guess)) / T{2};
        if (next == guess) {
            break;
        }
        guess = next;
    }
    return guess;
}

/// Cube root via Newton-Raphson; preserves sign for negative x.
template<typename T>
constexpr T cbrt(T x) {
    if (x == T{0}) {
        return T{0};
    }
    bool neg = x < T{0};
    if (neg) {
        x = -x;
    }
    T guess = x > T{1} ? x / T{3} : T{1};
    for (int i = 0; i < 50; ++i) {
        T next = ((T{2} * guess) + (x / (guess * guess))) / T{3};
        if (next == guess) {
            break;
        }
        guess = next;
    }
    return neg ? -guess : guess;
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

/// Cody–Waite reduction of @p x mod 2π into [−π, π] (two-part 2π).
template<typename T>
constexpr T reduce_two_pi(T x) {
    constexpr T two_pi_hi = static_cast<T>(6.28318530693650245668);
    constexpr T two_pi_lo = static_cast<T>(2.43084020260247689728e-10);
    constexpr T inv_two_pi = static_cast<T>(0.15915494309189533577);
    const T     kreal = static_cast<T>(lround_away(x * inv_two_pi));
    return (x - (kreal * two_pi_hi)) - (kreal * two_pi_lo);
}

/// Cosine via Taylor series after Cody–Waite 2π reduction to [−π, π].
template<typename T>
constexpr T cos(T x) {
    x = reduce_two_pi(x);

    T x2 = x * x;
    T result = T{1};
    T term = T{1};
    for (int n = 1; n <= 12; ++n) {
        term *= -x2 / T(2 * n * ((2 * n) - 1));
        result += term;
    }
    return result;
}

/// Sine via Taylor series after Cody–Waite 2π reduction to [−π, π].
template<typename T>
constexpr T sin(T x) {
    x = reduce_two_pi(x);

    T x2 = x * x;
    T result = x;
    T term = x;
    for (int n = 1; n <= 12; ++n) {
        term *= -x2 / T((2 * n) * ((2 * n) + 1));
        result += term;
    }
    return result;
}

/// {sin(x), cos(x)} with a single 2π reduction (shared remainder).
template<typename T>
constexpr sincos_result<T> sincos(T x) {
    x = reduce_two_pi(x);

    T x2 = x * x;

    T cos_r = T{1};
    T cos_term = T{1};
    for (int n = 1; n <= 12; ++n) {
        cos_term *= -x2 / T(2 * n * ((2 * n) - 1));
        cos_r += cos_term;
    }

    T sin_r = x;
    T sin_term = x;
    for (int n = 1; n <= 12; ++n) {
        sin_term *= -x2 / T((2 * n) * ((2 * n) + 1));
        sin_r += sin_term;
    }
    return {sin_r, cos_r};
}

/// Tangent via continued fraction after π reduction to [−π/2, π/2], with a
/// complementary-angle identity when |r| is near π/2 (slow CF convergence).
/// @see Cuyt et al., "Handbook of Continued Fractions for Special Functions" §12.1
template<typename T>
constexpr T tan(T x) {
    constexpr T half_pi = pi_v<T> / T{2};

    constexpr T pi_hi = static_cast<T>(3.14159265346825122834);
    constexpr T pi_lo = static_cast<T>(1.21542010130123844986e-10);
    constexpr T inv_pi = static_cast<T>(0.31830988618379067154);
    const T     kreal = static_cast<T>(lround_away(x * inv_pi));
    T           r = (x - (kreal * pi_hi)) - (kreal * pi_lo);

    // Near ±π/2 the continued fraction converges slowly. Work with the
    // complementary angle via the cotangent identity, valid for r ∈ (0, π/2):
    //     tan(r) = cot(π/2 − r) = 1 / tan(π/2 − r)
    // The result's sign follows the sign of r (tan is odd), applied at the end.
    T abs_r = r >= T{0} ? r : -r;
    if (abs_r > static_cast<T>(1.2)) {
        T             comp = half_pi - abs_r;
        constexpr int N = 20;

        T x2 = comp * comp;
        T cf = T((2 * N) + 1);
        for (int i = N - 1; i >= 0; --i) {
            cf = T((2 * i) + 1) - (x2 / cf);
        }
        T tan_comp = comp / cf;
        T result = T{1} / tan_comp;
        return r >= T{0} ? result : -result;
    }

    constexpr int N = 20;

    T x2 = r * r;
    T cf = T((2 * N) + 1);
    for (int i = N - 1; i >= 0; --i) {
        cf = T((2 * i) + 1) - (x2 / cf);
    }

    return r / cf;
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

    constexpr T     ln2_hi = static_cast<T>(6.93147180369123816490e-01);
    constexpr T     ln2_lo = static_cast<T>(1.90821492927058770002e-10);
    constexpr T     inv_ln2 = static_cast<T>(1.44269504088896340736);
    const long long k = lround_away(x * inv_ln2);
    const T         kreal = static_cast<T>(k);
    const T         r = (x - (kreal * ln2_hi)) - (kreal * ln2_lo);

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

/// Natural log via Newton–Raphson on @f$ e^y = x @f$, after reducing @p x into
/// roughly [1/e, e] by multiplying/dividing by e. Domain @p x ≤ 0 → 0.
template<typename T>
constexpr T log(T x) {
    if (x <= T{0}) {
        return T{0};
    }

    constexpr T e_val = static_cast<T>(2.718281828459045235360287);
    constexpr T e_inv = static_cast<T>(0.367879441171442321595524);

    T k = T{0};
    T y = x;
    while (y > e_val) {
        y *= e_inv;
        k += T{1};
    }
    while (y < e_inv) {
        y *= e_val;
        k -= T{1};
    }

    T guess = y - T{1};
    for (int i = 0; i < 50; ++i) {
        T e_guess = exp(guess);
        T next = guess - ((e_guess - y) / e_guess);
        if (abs(next - guess) < static_cast<T>(1e-15)) {
            break;
        }
        guess = next;
    }
    return guess + k;
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
