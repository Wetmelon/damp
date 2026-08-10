// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file trig.hpp
 * @brief Fast float/double sine, cosine, and inverse trig with Cody–Waite reduction
 */

#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"

// Fast float/double implementations live in damp::detail as fast_* so they do
// not overload-shadow the constexpr public API in math.hpp (damp::sin/sqrt/...).
// MathBackend (damp_backend.hpp) routes runtime calls here under DAMP_MATH_BACKEND_DAMP.

namespace damp::detail {

// ---------------------------------------------------------------------------
// Coefficient tables (float: compact PE; double: ~1 ULP poly on the chart)
// ---------------------------------------------------------------------------

/// sin(pi*g) = g * P(u), u = g^2, g in [-0.5, 0.5]
constexpr float sin_coeffs_f[] = {
    +3.141582250595093e+0f,
    -5.167152404785156e+0f,
    +2.541962385177612e+0f,
    -5.547642707824707e-1f,
};

constexpr double sin_coeffs_d[] = {
    +3.1415926535897936e+0,
    -5.167712780049931e+0,
    +2.5501640398685796e+0,
    -5.99264528959652e-1,
    +8.214588017721217e-2,
    -7.370371648375263e-3,
    +4.6600698271522526e-4,
    -2.115533601999618e-5,
};

/// cos(pi*g) = Q(u), u = g^2
constexpr float cos_coeffs_f[] = {
    +1.000000000000000e+0f,
    -4.934792041778564e+0f,
    +4.058376789093018e+0f,
    -1.331601023674011e+0f,
    +2.190602868795395e-1f,
};

constexpr double cos_coeffs_d[] = {
    +0.9999999999999993e+0,
    -4.934802200543971e+0,
    +4.058712126348079e+0,
    -1.3352627664956782e+0,
    +2.3533059161430153e-1,
    -2.5806550193892084e-2,
    +1.9279206065477292e-3,
    -1.0047210059254057e-4,
};

/// asin: pi/2 - sqrt(1-x) * P(x)
constexpr float asin_coeffs_f[] = {
    +1.570795416831970e+00f,
    -2.145179212093353e-01f,
    +8.791970461606979e-02f,
    -4.508748650550842e-02f,
    +1.950908824801445e-02f,
    -4.407065920531750e-03f,
};

constexpr double asin_coeffs_d[] = {
    +1.5707963267930571e+0,
    -2.146018359616614e-1,
    +8.904858538557413e-2,
    -5.079196164551241e-2,
    +3.3671015677161804e-2,
    -2.429970001459448e-2,
    +1.8317342074084073e-2,
    -1.3721477009347843e-2,
    +9.480409700593866e-3,
    -5.48307619207655e-3,
    +2.368429826759506e-3,
    -6.56209332741492e-4,
    +8.571307314657812e-5,
};

/// atan on [0, 1]
constexpr float atan_coeffs_f[] = {
    +4.543576608284638e-07f,
    +9.999714493751526e-01f,
    +1.813957787817344e-04f,
    -3.311897814273834e-01f,
    -2.631079405546188e-02f,
    +3.089837431907654e-01f,
    -2.172826975584030e-01f,
    +5.104487016797066e-02f,
};

constexpr double atan_coeffs_d[] = {
    -1.095369548925466e-12,
    +1.0000000005219774e+0,
    -4.042843901583747e-8,
    -3.33332133030723e-1,
    -1.7816212181167113e-5,
    +2.001464697290472e-1,
    -6.515048882582934e-4,
    -1.4193701247823626e-1,
    +6.35043399470227e-3,
    +6.569399464769006e-2,
    +1.4928600708464673e-1,
    -3.9339790236438776e-1,
    +3.8906225836684943e-1,
    -2.096182097764598e-1,
    +6.164370400752765e-2,
    -7.83008577700179e-3,
};

/**
 * @brief Half-period chart: x = pi * (period_index + frac), frac in [-0.5, 0.5]
 */
template<typename T>
struct Reduced {
    T   frac;
    int period_index;
};

template<typename T>
static inline Reduced<T> reduce(T x);
template<typename T>
static inline T sin_poly(T frac);
template<typename T>
static inline T cos_poly(T frac);

template<typename T, size_t N>
static inline T horner_eval(T x, const T (&coeffs)[N]);
template<typename T, size_t N>
static inline T estrin_eval(T x, const T (&coeffs)[N]);

template<typename T>
static inline T reassoc_barrier(T v);

// --- Fast scalar math (bound by MathBackend when DAMP_MATH_BACKEND_DAMP) ---

/**
 * @brief Fast square root (float)
 *
 * Emits VSQRT on ARM FP targets; otherwise @c __builtin_sqrtf. Returns NaN for
 * x < 0 (IEEE FPU behaviour: no errno, no trap).
 */
inline float fast_sqrt(float x) {
#if defined(__ARM_FP)
    float result;
    asm("vsqrt.f32 %0, %1" : "=t"(result) : "t"(x));
    return result;
#else
    return __builtin_sqrtf(x);
#endif
}

/**
 * @brief Fast square root (double)
 *
 * Uses @c __builtin_sqrt (hardware when the target provides it).
 */
inline double fast_sqrt(double x) {
    return __builtin_sqrt(x);
}

/**
 * @brief Round-to-nearest, ties to even (float)
 *
 * ARM VFP: float→int32→float via vcvtr/vcvt. Valid for |x| in the int32 range
 * (covers range-reduction quotients used here). Otherwise @c __builtin_nearbyintf.
 */
inline float fast_nearbyint(float x) {
#if defined(__ARM_FP)
    float result;
    asm("vcvtr.s32.f32 %0, %1\n\t"
        "vcvt.f32.s32 %0, %0"
        : "=&t"(result)
        : "t"(x));
    return result;
#else
    return __builtin_nearbyintf(x);
#endif
}

/**
 * @brief Round-to-nearest, ties to even (double)
 *
 * Uses the compiler builtin (int32 VFP convert would truncate large quotients).
 */
inline double fast_nearbyint(double x) {
    return __builtin_nearbyint(x);
}

/**
 * @brief Sine (float), full-range Cody–Waite + poly; ~8 ULP
 * @see sincos
 */
inline float fast_sin(float angle_rad) {
    auto [frac, period_index] = reduce(angle_rad);
    float s = sin_poly(frac);
    return ((period_index & 1) != 0) ? -s : s;
}

/**
 * @brief Sine (double), full-range Cody–Waite + poly
 * @see sincos
 */
inline double fast_sin(double angle_rad) {
    auto [frac, period_index] = reduce(angle_rad);
    double s = sin_poly(frac);
    return ((period_index & 1) != 0) ? -s : s;
}

/**
 * @brief Cosine (float)
 * @see sincos
 */
inline float fast_cos(float angle_rad) {
    auto [frac, period_index] = reduce(angle_rad);
    float c = cos_poly(frac);
    return ((period_index & 1) != 0) ? -c : c;
}

/**
 * @brief Cosine (double)
 * @see sincos
 */
inline double fast_cos(double angle_rad) {
    auto [frac, period_index] = reduce(angle_rad);
    double c = cos_poly(frac);
    return ((period_index & 1) != 0) ? -c : c;
}

/**
 * @brief Paired sin/cos (float) from one reduction
 */
inline damp::pair<float, float> fast_sincos(float angle_rad) {
    auto [frac, period_index] = reduce(angle_rad);
    float s = sin_poly(frac);
    float c = cos_poly(frac);
    if ((period_index & 1) != 0) {
        s = -s;
        c = -c;
    }
    return {s, c};
}

/**
 * @brief Paired sin/cos (double) from one reduction
 */
inline damp::pair<double, double> fast_sincos(double angle_rad) {
    auto [frac, period_index] = reduce(angle_rad);
    double s = sin_poly(frac);
    double c = cos_poly(frac);
    if ((period_index & 1) != 0) {
        s = -s;
        c = -c;
    }
    return {s, c};
}

/**
 * @brief Arcsine, x in [-1, 1] → [-pi/2, pi/2] (float; clamped)
 */
inline float fast_asin(float x) {
    bool negate = false;
    if (x < 0.0f) {
        x = -x;
        negate = true;
    }
    x = damp::min(x, 1.0f);
    float p = estrin_eval(x, detail::asin_coeffs_f);
    float result = (damp::numbers::pi_v<float> / 2.0f) - (fast_sqrt(1.0f - x) * p);
    return negate ? -result : result;
}

/**
 * @brief Arcsine (double; clamped to [-1, 1])
 */
inline double fast_asin(double x) {
    bool negate = false;
    if (x < 0.0) {
        x = -x;
        negate = true;
    }
    x = damp::min(x, 1.0);
    double p = estrin_eval(x, detail::asin_coeffs_d);
    double result = (damp::numbers::pi_v<double> / 2.0) - (fast_sqrt(1.0 - x) * p);
    return negate ? -result : result;
}

/**
 * @brief Arccosine, x in [-1, 1] → [0, pi] (float; clamped)
 */
inline float fast_acos(float x) {
    bool negate = false;
    if (x < 0.0f) {
        x = -x;
        negate = true;
    }
    x = damp::min(x, 1.0f);
    float p = estrin_eval(x, detail::asin_coeffs_f);
    float result = fast_sqrt(1.0f - x) * p;
    return negate ? (damp::numbers::pi_v<float> - result) : result;
}

/**
 * @brief Arccosine (double; clamped to [-1, 1])
 */
inline double fast_acos(double x) {
    bool negate = false;
    if (x < 0.0) {
        x = -x;
        negate = true;
    }
    x = damp::min(x, 1.0);
    double p = estrin_eval(x, detail::asin_coeffs_d);
    double result = fast_sqrt(1.0 - x) * p;
    return negate ? (damp::numbers::pi_v<double> - result) : result;
}

/**
 * @brief Arctangent → [-pi/2, pi/2] (float)
 *
 * |x| > 1 reduced via atan(x) = pi/2 - atan(1/x).
 */
inline float fast_atan(float x) {
    bool negate = false;
    bool complement = false;
    if (x < 0.0f) {
        x = -x;
        negate = true;
    }
    if (__builtin_fabsf(x) > 1.0f) {
        x = 1.0f / x;
        complement = true;
    }
    float result = estrin_eval(x, detail::atan_coeffs_f);
    if (complement) {
        result = (damp::numbers::pi_v<float> / 2.0f) - result;
    }
    return negate ? -result : result;
}

/**
 * @brief Arctangent (double)
 */
inline double fast_atan(double x) {
    bool negate = false;
    bool complement = false;
    if (x < 0.0) {
        x = -x;
        negate = true;
    }
    if (__builtin_fabs(x) > 1.0) {
        x = 1.0 / x;
        complement = true;
    }
    double result = estrin_eval(x, detail::atan_coeffs_d);
    if (complement) {
        result = (damp::numbers::pi_v<double> / 2.0) - result;
    }
    return negate ? -result : result;
}

/**
 * @brief Two-argument arctangent → [-pi, pi] (float)
 */
inline float fast_atan2(float y, float x) {
    constexpr float min_normal = 1.17549435082228750797e-38f;
    float           ax = __builtin_fabsf(x);
    float           ay = __builtin_fabsf(y);
    float           lo = damp::min(ax, ay);
    float           hi = damp::max(ax, ay) + min_normal;
    float           t = estrin_eval(lo / hi, detail::atan_coeffs_f);
    float           r = (ay > ax) ? ((damp::numbers::pi_v<float> / 2.0f) - t) : t;
    r = (x >= 0.0f) ? r : (damp::numbers::pi_v<float> - r);
    return __builtin_copysignf(r, y);
}

/**
 * @brief Two-argument arctangent (double)
 */
inline double fast_atan2(double y, double x) {
    constexpr double min_normal = 2.2250738585072014e-308;
    double           ax = __builtin_fabs(x);
    double           ay = __builtin_fabs(y);
    double           lo = damp::min(ax, ay);
    double           hi = damp::max(ax, ay) + min_normal;
    double           t = estrin_eval(lo / hi, detail::atan_coeffs_d);
    double           r = (ay > ax) ? ((damp::numbers::pi_v<double> / 2.0) - t) : t;
    r = (x >= 0.0) ? r : (damp::numbers::pi_v<double> - r);
    return __builtin_copysign(r, y);
}

template<typename T>
static inline T reassoc_barrier(T v) {
#if defined(__has_builtin)
#if __has_builtin(__builtin_assoc_barrier)
    return __builtin_assoc_barrier(v);
#endif
#endif

#if (defined(__GNUC__) || defined(__clang__)) && defined(__ARM_FP)
    if constexpr (std::is_same_v<T, float>) {
        asm volatile("" : "+t"(v));
    } else {
        asm volatile("" : "+w"(v));
    }
#endif
    return v;
}

/**
 * @brief Cody–Waite reduction of x (radians) to {frac, period_index}
 *
 * Float: three π words, accurate to |x| ~ 2.5e4 rad.
 * Double: two π words (low part exact after HI truncation), large-arg path.
 */
template<typename T>
static inline Reduced<T> reduce(T x) {
    constexpr T inv_pi = damp::numbers::inv_pi_v<T>;

    T   n = fast_nearbyint(x * inv_pi);
    int period_index = static_cast<int>(n);

    T r;
    if constexpr (std::is_same_v<T, float>) {
        constexpr T PI_HI = 3.140625f;
        constexpr T PI_LO = 9.6765358467e-04f;
        constexpr T PI_LO2 = 5.1265658385e-12f;
        r = reassoc_barrier(x - (n * PI_HI));
        r = reassoc_barrier(r - (n * PI_LO));
        r = reassoc_barrier(r - (n * PI_LO2));
    } else {
        // HI: π with low 21 mantissa bits cleared (exact n*HI for moderate n)
        constexpr T PI_HI = 3.1415926534682512;
        constexpr T PI_LO = 1.2154188766544394e-10;
        r = reassoc_barrier(x - (n * PI_HI));
        r = reassoc_barrier(r - (n * PI_LO));
    }

    return {.frac = r * inv_pi, .period_index = period_index};
}

template<typename T>
static inline T sin_poly(T frac) {
    const T u = frac * frac;
    if constexpr (std::is_same_v<T, float>) {
        const T u2 = u * u;
        T       p = sin_coeffs_f[0];
        p += sin_coeffs_f[1] * u;
        p += sin_coeffs_f[2] * u2;
        p += sin_coeffs_f[3] * u2 * u;
        return frac * p;
    } else {
        // Horner in u (8 coeffs)
        T p = sin_coeffs_d[7];
        p = sin_coeffs_d[6] + u * p;
        p = sin_coeffs_d[5] + u * p;
        p = sin_coeffs_d[4] + u * p;
        p = sin_coeffs_d[3] + u * p;
        p = sin_coeffs_d[2] + u * p;
        p = sin_coeffs_d[1] + u * p;
        p = sin_coeffs_d[0] + u * p;
        return frac * p;
    }
}

template<typename T>
static inline T cos_poly(T frac) {
    const T u = frac * frac;
    if constexpr (std::is_same_v<T, float>) {
        const T u2 = u * u;
        T       p = cos_coeffs_f[0];
        p += cos_coeffs_f[1] * u;
        p += cos_coeffs_f[2] * u2;
        p += cos_coeffs_f[3] * u2 * u;
        p += cos_coeffs_f[4] * u2 * u2;
        return p;
    } else {
        T p = cos_coeffs_d[7];
        p = cos_coeffs_d[6] + u * p;
        p = cos_coeffs_d[5] + u * p;
        p = cos_coeffs_d[4] + u * p;
        p = cos_coeffs_d[3] + u * p;
        p = cos_coeffs_d[2] + u * p;
        p = cos_coeffs_d[1] + u * p;
        p = cos_coeffs_d[0] + u * p;
        return p;
    }
}

template<typename T, size_t N>
static inline T horner_eval(T x, const T (&coeffs)[N]) {
    static_assert(N >= 1, "horner_eval requires at least one coefficient");
    T p = coeffs[N - 1];
    for (int i = static_cast<int>(N) - 2; i >= 0; --i) {
        p = coeffs[i] + (x * p);
    }
    return p;
}

/**
 * @brief Estrin polynomial evaluation (float or double)
 */
template<typename T, size_t N>
static inline T estrin_eval(T x, const T (&coeffs)[N]) {
    static_assert(N >= 1, "estrin_eval requires at least one coefficient");
    T b[N];
    for (size_t i = 0; i < N; ++i) {
        b[i] = coeffs[i];
    }
    T power = x;
    for (size_t n = N; n > 1; n = (n + 1) / 2) {
        for (size_t i = 0; i < n / 2; ++i) {
            b[i] = b[2 * i] + (power * b[(2 * i) + 1]);
        }
        if (n % 2 != 0) {
            b[n / 2] = b[n - 1];
        }
        power = power * power;
    }
    return b[0];
}

} // namespace damp::detail
