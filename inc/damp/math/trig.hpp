// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file trig.hpp
 * @brief Fast float and double sine, cosine and inverse trig, with Cody–Waite reduction
 *
 * Measured absolute error, float: sin 8.9e-7 and cos 3.0e-7 for |x| ≤ 2.5e4 rad
 * (8.8e-7 and 2.7e-7 on [0, 2π]), asin 3.1e-7, acos 4.8e-7, atan 2.5e-7, atan2 4.2e-7
 * (atan and asin are odd with relative error ≤ 4e-7). Double: sin 1e-15, cos 3e-15
 * for |x| ≤ 6e6 rad; asin 3.3 ulp, acos 3.1 ulp, atan and atan2 2.1 ulp against a
 * quad-precision reference (relative error ≤ 6.8e-16; atan, asin and atan2 odd, exactly
 * 0 at 0). Float sin/cos: exhaustive over all floats; double: 3.2e7 samples per domain
 * (analysis/trig).
 */


#include "constexpr_math.hpp" // detail::isfinite (bit pattern, survives -ffast-math), detail::assoc_barrier
#include "damp/backend.hpp"

// Fast float/double implementations live in damp::detail as fast_* so they do
// not overload-shadow the constexpr public API in math.hpp (damp::sin/sqrt/...).
// MathBackend (damp_backend.hpp) routes runtime calls here under DAMP_MATH_BACKEND_DAMP.

namespace damp::detail {

// ---------------------------------------------------------------------------
// Coefficient tables (float: compact minimax; double sin 1e-15, cos 3e-15 on the chart)
// ---------------------------------------------------------------------------

/** @brief sin(pi*g) = g * P(u), u = g^2, g in [-0.5, 0.5] */
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

/** @brief cos(pi*g) = Q(u), u = g^2 */
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

/** @brief asin(x) = pi/2 - sqrt(1-x) * P(x) for x in (0.5, 1]; acos(x) = sqrt(1-x) * P(x) on [0, 1] */
constexpr float asin_coeffs_f[] = {
    +1.570796251296997e+00f,
    -2.145862877368927e-01f,
    +8.876860141754150e-02f,
    -4.885582253336906e-02f,
    +2.699527330696583e-02f,
    -1.122441701591015e-02f,
    +2.320223953574896e-03f,
};


/** @brief asin(x) = x * P(u), u = x^2, x in [0, 0.5] (relative minimax) */
constexpr float asin_small_coeffs_f[] = {
    +1.0000001192092896e+00f,
    +1.6664780676364900e-01f,
    +7.5550585985183716e-02f,
    +3.9067719131708145e-02f,
    +5.1992245018482208e-02f,
};

/** @brief atan(x) = x * P(u), u = x^2, x in [0, 1] (relative minimax; odd, so atan(0) = 0) */
constexpr float atan_coeffs_f[] = {
    +9.9999994039535522e-01f,
    -3.3331951498985291e-01f,
    +1.9967921078205109e-01f,
    -1.4007787406444550e-01f,
    +9.8823666572570801e-02f,
    -5.9047453105449677e-02f,
    +2.3951398208737373e-02f,
    -4.6112081035971642e-03f,
};

/**
 * @brief asin(x) = x + x·u·Q(u), u = x², |x| ≤ 0.5; above it asin|x| = π/2 − 2√z·(1 + z·Q(z)), z = (1 − |x|)/2
 *
 * Q of 1 + u·Q(u) ≈ asin(√u)/√u, relative minimax on u ∈ [0, 0.25]
 * (analysis/trig/tools/generate_coefficients.py --refit-double).
 */
constexpr double asin_q_coeffs_d[] = {
    0x1.5555555556cffp-3,
    0x1.33333330916e9p-4,
    0x1.6db6dd0de0b46p-5,
    0x1.f1c69c914b0e7p-6,
    0x1.6e97028548c24p-6,
    0x1.1bacceb5ccf35p-6,
    0x1.d54ed79064b4ep-7,
    0x1.336c86a4eee75p-7,
    0x1.2b086fb51a379p-6,
    -0x1.7f689d2a7f1ccp-7,
    0x1.0234a3e6f786cp-5,
};

/**
 * @brief atan(t) = t + t·u·Q(u), u = t², |t| ≤ tan(π/8)
 *
 * Q of 1 + u·Q(u) ≈ atan(√u)/√u, relative minimax on u ∈ [0, 0.1716]
 * (analysis/trig/tools/generate_coefficients.py --refit-double).
 */
constexpr double atan_q_coeffs_d[] = {
    -0x1.5555555555317p-2,
    0x1.99999999495c7p-3,
    -0x1.2492490574f9bp-3,
    0x1.c71c65ccef300p-4,
    -0x1.745bc2f167815p-4,
    0x1.3afc31e0a22cfp-4,
    -0x1.100510529a9c5p-4,
    0x1.d253b7ed1dd6ep-5,
    -0x1.651489d4fb63cp-5,
    0x1.5c4794e6d57d5p-6,
};

/** @brief k·π/4 for k = 0…4, rounded to double (hi) and the remainder k·π/4 − hi (lo) */
constexpr double quarter_pi_hi_d[] = {
    0x0.0p+0,
    0x1.921fb54442d18p-1,
    0x1.921fb54442d18p+0,
    0x1.2d97c7f3321d2p+1,
    0x1.921fb54442d18p+1,
};
constexpr double quarter_pi_lo_d[] = {
    0x0.0p+0,
    0x1.1a62633145c07p-55,
    0x1.1a62633145c07p-54,
    0x1.a79394c9e8a0ap-54,
    0x1.1a62633145c07p-53,
};

/**
 * @brief Half-period chart: x = pi * (period_index + frac), frac in [-0.5, 0.5]
 */
template<typename T>
struct Reduced {
    T   frac;
    int period_index;
};

inline Reduced<float>  reduce(float x);
inline Reduced<double> reduce(double x);
template<typename T>
inline T sin_poly(T frac);
template<typename T>
inline T cos_poly(T frac);

template<typename T, size_t N>
constexpr T estrin_eval(T x, const T (&coeffs)[N]);
template<typename T, size_t N>
inline T horner2_eval(T x, const T (&coeffs)[N]);

// --- Fast scalar math (bound by MathBackend when DAMP_MATH_BACKEND_DAMP) ---

/**
 * @brief Fast square root (float)
 *
 * Emits VSQRT on 32-bit ARM FP targets; otherwise @c __builtin_sqrtf. Returns NaN
 * for x < 0 (IEEE FPU behaviour: no errno, no trap).
 */
inline float fast_sqrt(float x) {
#if defined(__arm__) && defined(__ARM_FP) // AArch64 has no "t" constraint
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
 * 32-bit ARM VFP: float→int32→float via vcvtr/vcvt for |x| < 2²³; larger magnitudes
 * (already integers) and NaN are returned unchanged, so the int32 convert never
 * saturates. Otherwise @c __builtin_nearbyintf. Full float range.
 */
inline float fast_nearbyint(float x) {
#if defined(__arm__) && defined(__ARM_FP)
    // |x| ≥ 2²³, ±Inf and NaN, on the bit pattern: -ffast-math folds a float comparison
    // and lets NaN reach the convert (which returns 0).
    if ((damp::bit_cast<uint32_t>(x) & 0x7fffffffU) >= 0x4b000000U) {
        return x;
    }
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
 * @brief −v for an odd period index, else v (float)
 *
 * 32-bit ARM VFP: a VNEG predicated on the index's low bit, no branch (GCC -O2 for
 * Cortex-M7 branches around a plain VNEG). Thumb-2 needs the IT; in ARM state the VFP
 * instruction is conditional by encoding.
 */
inline float negate_if_odd(float v, int period_index) {
#if defined(__arm__) && defined(__ARM_FP)
    asm("lsls %1, %1, #31\n\t"
#if defined(__thumb__)
        "it ne\n\t"
#endif
        "vnegne.f32 %0, %0"
        : "+t"(v), "+r"(period_index)
        :
        : "cc");
    return v;
#else
    return ((period_index & 1) != 0) ? -v : v;
#endif
}

/**
 * @brief Sine (float), Cody–Waite + poly; argument range as @ref reduce
 * @see sincos
 */
inline float fast_sin(float angle_rad) {
    auto [frac, period_index] = reduce(angle_rad);
    return negate_if_odd(sin_poly(frac), period_index);
}

/**
 * @brief Sine (double), Cody–Waite + poly; argument range as @ref reduce
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
    return negate_if_odd(cos_poly(frac), period_index);
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
    // Branch kept: a two-value asm flip made GCC 14 add three VMOVs to the polynomials (+8 M7 cycles)
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
 *
 * Odd: |x| ≤ 0.5 uses x·P(x²) (relative error ≤ 4e-7, asin(0) = 0); above it
 * pi/2 − √(1−x)·P(x).
 *
 * @see "Approximations for Digital Computers" (Hastings, 1955)
 * @see "Handbook of Mathematical Functions" (Abramowitz & Stegun, 1964), §4.4.45
 */
inline float fast_asin(float x) {
    bool negate = false;
    if (x < 0.0f) {
        x = -x;
        negate = true;
    }
    float result = 0.0f;
    if (x <= 0.5f) {
        result = x * estrin_eval(x * x, detail::asin_small_coeffs_f);
    } else {
        x = damp::min(x, 1.0f);
        float p = estrin_eval(x, detail::asin_coeffs_f);
        result = (damp::numbers::pi_v<float> / 2.0f) - (fast_sqrt(1.0f - x) * p);
    }
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
 * @brief Arctangent → [-pi/2, pi/2] (float)
 *
 * Odd: x·P(x²) on [0, 1] (relative error ≤ 4e-7, atan(0) = 0); |x| > 1 reduced
 * via atan(x) = pi/2 - atan(1/x).
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
    float result = x * estrin_eval(x * x, detail::atan_coeffs_f);
    if (complement) {
        result = (damp::numbers::pi_v<float> / 2.0f) - result;
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
    float           hi = damp::max(damp::max(ax, ay), min_normal); // 0/0 guard, exact ratio otherwise
    float           q = lo / hi;
    float           t = q * estrin_eval(q * q, detail::atan_coeffs_f);
    float           r = (ay > ax) ? ((damp::numbers::pi_v<float> / 2.0f) - t) : t;
    r = (x >= 0.0f) ? r : (damp::numbers::pi_v<float> - r);
    return __builtin_copysignf(r, y);
}

/**
 * @brief k·π/4 + s·(1 + u·Q(u)) as (k·π/4 + s) + (lo + s·u·Q(u)), k = 0…4
 *
 * The leading sum is exact or rounded once at the result's scale; the low word of
 * k·π/4 joins the small term, so the polynomial's error is scaled by s·u.
 */
template<size_t N>
inline double quarter_pi_plus(int k, double s, double u, const double (&q_coeffs)[N]) {
    const double head = assoc_barrier(quarter_pi_hi_d[k] + s);
    const double tail = assoc_barrier(quarter_pi_lo_d[k] + ((s * u) * horner2_eval(u, q_coeffs)));
    return head + tail;
}

/**
 * @brief Arcsine, x in [-1, 1] → [-pi/2, pi/2] (double; clamped)
 *
 * Odd: |x| ≤ 0.5 uses x + x³·Q(x²), so asin(0) = 0 and tiny x return x; above it
 * asin|x| = π/2 − 2√z·(1 + z·Q(z)) with z = (1 − |x|)/2 ≤ 0.25, the same polynomial.
 */
inline double fast_asin(double x) {
    const double ax = damp::min(__builtin_fabs(x), 1.0);
    double       u = ax * ax;
    double       s = ax;
    int          k = 0; // asin|x| = k·π/4 + s·(1 + u·Q(u))
    if (ax > 0.5) {
        u = (1.0 - ax) * 0.5;
        s = -2.0 * fast_sqrt(u);
        k = 2;
    }
    return __builtin_copysign(quarter_pi_plus(k, s, u, asin_q_coeffs_d), x);
}

/**
 * @brief Arccosine, x in [-1, 1] → [0, pi] (double; clamped)
 *
 * π/2 − asin(x) for |x| ≤ 0.5; 2√z·(1 + z·Q(z)) above it and π minus that below
 * −0.5, z = (1 − |x|)/2, so acos(1) = 0. Same polynomial as @ref fast_asin.
 */
inline double fast_acos(double x) {
    const double ax = damp::min(__builtin_fabs(x), 1.0);
    double       u = ax * ax;
    double       s = -x;
    int          k = 2; // acos x = k·π/4 + s·(1 + u·Q(u))
    if (ax > 0.5) {
        u = (1.0 - ax) * 0.5;
        s = 2.0 * fast_sqrt(u);
        k = 0;
        if (x < 0.0) {
            s = -s;
            k = 4;
        }
    }
    return quarter_pi_plus(k, s, u, asin_q_coeffs_d);
}

/**
 * @brief Arctangent → [-pi/2, pi/2] (double)
 *
 * Odd: t + t³·Q(t²) for |t| ≤ tan(π/8), so atan(0) = 0 and tiny x return x. Larger
 * |x| use atan|x| = π/4 + atan((|x| − 1)/(|x| + 1)) up to tan(3π/8) and
 * π/2 + atan(−1/|x|) above it: one division at most.
 */
inline double fast_atan(double x) {
    constexpr double tan_pi_8 = 0x1.a827999fcef32p-2;  // √2 − 1
    constexpr double tan_3pi_8 = 0x1.3504f333f9de6p+1; // √2 + 1
    const double     ax = __builtin_fabs(x);
    double           t = ax;
    int              k = 0; // atan|x| = k·π/4 + atan(t)
    if (ax > tan_pi_8) {
        if (ax > tan_3pi_8) {
            t = -1.0 / ax;
            k = 2;
        } else {
            t = (ax - 1.0) / (ax + 1.0);
            k = 1;
        }
    }
    return __builtin_copysign(quarter_pi_plus(k, t, t * t, atan_q_coeffs_d), x);
}

/**
 * @brief Two-argument arctangent → [-pi, pi] (double)
 *
 * Same reduction as @ref fast_atan on min(|x|, |y|)/max(|x|, |y|), with one division:
 * (lo − hi)/(lo + hi) replaces lo/hi above tan(π/8). The quadrant fold adds a
 * multiple of π/4 once, in two words.
 */
inline double fast_atan2(double y, double x) {
    constexpr double tan_pi_8 = 0x1.a827999fcef32p-2;
    constexpr double min_normal = 0x1p-1022;
    const double     ax = __builtin_fabs(x);
    const double     ay = __builtin_fabs(y);
    const bool       swap = ay > ax;
    const double     lo = swap ? ax : ay;
    const double     hi = damp::max(swap ? ay : ax, min_normal); // 0/0 guard
    double           num = lo;
    double           den = hi;
    int              k = 0; // atan(lo/hi) = k·π/4 + atan(num/den)
    if (lo > tan_pi_8 * hi) {
        num = lo - hi;
        den = lo + hi;
        k = 1;
    }
    double t = num / den;
    if (swap) { // π/2 − atan(lo/hi)
        t = -t;
        k = 2 - k;
    }
    if (x < 0.0) { // π − …
        t = -t;
        k = 4 - k;
    }
    return __builtin_copysign(quarter_pi_plus(k, t, t * t, atan_q_coeffs_d), y);
}

/**
 * @brief Cody–Waite reduction of x (radians) to {frac, period_index} (float)
 *
 * Three π words. PI_HI (8 significant bits) and PI_MID (7) make n·PI_HI and n·PI_MID
 * exact for |n| < 2¹⁶, and both subtractions are exact: x − n·PI_HI − n·PI_MID is a
 * multiple of min(ulp(x), 2⁻¹⁷) smaller than 2 in magnitude. Only n·PI_LO
 * (|PI_LO| < 1.3e-6) is rounded, and no step needs an FMA, so targets with and without
 * one reduce alike. Accurate to |x| ~ 2.5e4 rad, where x·(1/π) starts to round n off by
 * one, then an error of about ulp(x). One unsigned compare on |x|'s bit pattern sends NaN,
 * ±Inf and every x whose quotient x/π rounds to 2²¹ or more (|x| ≳ 6.6e6 rad) to the
 * cold path, so the hot path has no float compare, rounds with the int32 convert (which
 * cannot saturate there) and takes the period index from that integer.
 *
 * Past 2²¹ half-turns ulp(x) is of order 1 rad and the angle has no fractional digits
 * left: the result is a multiple of π (frac 0) with n's parity, so sin and cos stay in
 * [−1, 1]. NaN and ±Inf return frac = x, so the polynomials give NaN.
 *
 * @see "Software Manual for the Elementary Functions" (Cody & Waite, 1980), Ch. 8
 */
inline Reduced<float> reduce(float x) {
    constexpr float    inv_pi = damp::numbers::inv_pi_v<float>;
    constexpr uint32_t cold_bits = 0x4ac90fd8U; // smallest |x| with fl(|x|·(1/π)) rounding to 2²¹
    const uint32_t     abs_bits = opaque_bits(damp::bit_cast<uint32_t>(x) & 0x7fffffffU);
    if (abs_bits >= cold_bits) [[unlikely]] {
        if (abs_bits >= 0x7f800000U) {
            return {.frac = x, .period_index = 0}; // NaN, ±Inf
        }
        const float n = fast_nearbyint(x * inv_pi);
        const float parity = n - (2.0f * fast_nearbyint(n * 0.5f));
        return {.frac = 0.0f, .period_index = (parity != 0.0f) ? 1 : 0};
    }
#if defined(__arm__) && defined(__ARM_FP)
    float q = x * inv_pi;
    asm("vcvtr.s32.f32 %0, %1" : "=t"(q) : "t"(q)); // round to nearest; |q| < 2²¹
    const int period_index = damp::bit_cast<int>(q);
    float     n;
    asm("vcvt.f32.s32 %0, %1" : "=t"(n) : "t"(q));
#else
    const float n = fast_nearbyint(x * inv_pi);
    const int   period_index = static_cast<int>(n);
#endif
    constexpr float PI_HI = 3.140625f;             // 201·2⁻⁶
    constexpr float PI_MID = 9.6893310546875e-04f; // 127·2⁻¹⁷
    constexpr float PI_LO = -1.2795156755e-06f;    // π − PI_HI − PI_MID
    const float     r_hi = assoc_barrier(x - assoc_barrier(n * PI_HI));
    const float     r_mid = assoc_barrier(r_hi - assoc_barrier(n * PI_MID));
    const float     r = assoc_barrier(r_mid - assoc_barrier(n * PI_LO));
    return {.frac = r * inv_pi, .period_index = period_index};
}

/**
 * @brief (x − n·π)/π in two π words (double); exact n·HI for |n| < 2²¹
 */
inline double pi_remainder(double x, double n) {
    // HI: π with low 21 mantissa bits cleared; LO = π − HI.
    constexpr double PI_HI = 3.1415926534682512;
    constexpr double PI_LO = 1.2154201013012384e-10;
    double           r = assoc_barrier(x - (n * PI_HI));
    r = assoc_barrier(r - (n * PI_LO));
    return r * damp::numbers::inv_pi_v<double>;
}

/**
 * @brief Cold path of reduce(double): |x| ≳ 2³⁰·π, ±Inf and NaN
 *
 * Out of line, so the hot path saves no registers for its library calls (targets
 * without a double round-to-integral instruction call nearbyint).
 */
[[gnu::noinline, gnu::cold]] inline Reduced<double> reduce_cold(double x) {
    if (!isfinite(x)) {
        return {.frac = x, .period_index = 0}; // NaN, ±Inf: the polynomials give NaN
    }
    const double n = fast_nearbyint(x * damp::numbers::inv_pi_v<double>);
    const double parity = n - (2.0 * fast_nearbyint(n * 0.5));
    const int    period_index = (parity != 0.0) ? 1 : 0;
    if (!(__builtin_fabs(n) < 0x1p50)) {
        return {.frac = 0.0, .period_index = period_index};
    }
    return {.frac = pi_remainder(x, n), .period_index = period_index};
}

/**
 * @brief Cody–Waite reduction of x (radians) to {frac, period_index} (double)
 *
 * Two π words, the first exact for |n| < 2²¹, then an error of about ulp(x). Below
 * about 2³⁰·π the high word of |x| alone selects the hot path, which rounds with the
 * int32 convert and takes the period index from that integer (no int64 conversion,
 * no library call). Above it, n is rounded in double; past 2⁵⁰ half-turns
 * (≳ 3.5e15 rad) the result is a multiple of π (frac 0) with n's parity. NaN and
 * ±Inf return frac = x.
 *
 * @see "Software Manual for the Elementary Functions" (Cody & Waite, 1980), Ch. 8
 */
inline Reduced<double> reduce(double x) {
    constexpr uint32_t cold_high = 0x41e921fbU; // high word of 2³⁰·π
    const uint32_t     abs_high = opaque_bits(static_cast<uint32_t>(damp::bit_cast<uint64_t>(x) >> 32U) & 0x7fffffffU);
    if (abs_high >= cold_high) [[unlikely]] {
        return reduce_cold(x);
    }
    const double q = x * damp::numbers::inv_pi_v<double>;
#if defined(__arm__) && defined(__ARM_FP) && ((__ARM_FP & 8) != 0)
    float qi;                                         // int32 in an S register
    asm("vcvtr.s32.f64 %0, %P1" : "=t"(qi) : "w"(q)); // round to nearest; |q| < 2³⁰
    const int period_index = damp::bit_cast<int>(qi);
    double    n;
    asm("vcvt.f64.s32 %P0, %1" : "=w"(n) : "t"(qi));
#else
    const double n = fast_nearbyint(q);
    const int    period_index = static_cast<int>(n);
#endif
    return {.frac = pi_remainder(x, n), .period_index = period_index};
}

template<typename T>
inline T sin_poly(T frac) {
    const T u = frac * frac;
    if constexpr (damp::is_same_v<T, float>) {
        const T u2 = u * u;
        T       p = sin_coeffs_f[0];
        p += sin_coeffs_f[1] * u;
        p += sin_coeffs_f[2] * u2;
        p += sin_coeffs_f[3] * u2 * u;
        return frac * p;
    } else {
        // Horner in u (8 coeffs)
        T p = sin_coeffs_d[7];
        p = sin_coeffs_d[6] + (u * p);
        p = sin_coeffs_d[5] + (u * p);
        p = sin_coeffs_d[4] + (u * p);
        p = sin_coeffs_d[3] + (u * p);
        p = sin_coeffs_d[2] + (u * p);
        p = sin_coeffs_d[1] + (u * p);
        p = sin_coeffs_d[0] + (u * p);
        return frac * p;
    }
}

template<typename T>
inline T cos_poly(T frac) {
    const T u = frac * frac;
    if constexpr (damp::is_same_v<T, float>) {
        const T u2 = u * u;
        T       p = cos_coeffs_f[0];
        p += cos_coeffs_f[1] * u;
        p += cos_coeffs_f[2] * u2;
        p += cos_coeffs_f[3] * u2 * u;
        p += cos_coeffs_f[4] * u2 * u2;
        return p;
    } else {
        T p = cos_coeffs_d[7];
        p = cos_coeffs_d[6] + (u * p);
        p = cos_coeffs_d[5] + (u * p);
        p = cos_coeffs_d[4] + (u * p);
        p = cos_coeffs_d[3] + (u * p);
        p = cos_coeffs_d[2] + (u * p);
        p = cos_coeffs_d[1] + (u * p);
        p = cos_coeffs_d[0] + (u * p);
        return p;
    }
}

/** @brief Terms left before level K of an N-term Estrin evaluation: ⌈N / 2ᴷ⌉ */
constexpr size_t estrin_terms(size_t n, size_t level) {
    return ((n - 1) >> level) + 1;
}

/**
 * @brief One Estrin level on b[0…Terms−1]: b[i] = b[2i] + p·b[2i + 1], an odd last term carried up;
 * then p = p² if another level follows
 *
 * Always inlined: a call per level would pass b through memory.
 */
template<size_t Terms, typename T, size_t N, size_t... I>
[[gnu::always_inline]] constexpr void estrin_level(T& p, T (&b)[N], damp::index_sequence<I...> /*pairs*/) {
    ((b[I] = assoc_barrier(b[2 * I] + (p * b[(2 * I) + 1]))), ...);
    if constexpr (Terms % 2 != 0) {
        b[Terms / 2] = b[Terms - 1];
    }
    if constexpr (Terms > 2) {
        p = assoc_barrier(p * p);
    }
}

/** @brief estrin_eval on a working copy b of the coefficients; level K multiplies by p = x^(2ᴷ) */
template<typename T, size_t N, size_t... I, size_t... K>
constexpr T estrin_levels(T x, const T (&coeffs)[N], damp::index_sequence<I...> /*terms*/, damp::index_sequence<K...> /*levels*/) {
    T                  b[N] = {coeffs[I]...};
    [[maybe_unused]] T p = x; // N = 1 has no level
    (estrin_level<estrin_terms(N, K)>(p, b, damp::make_index_sequence<estrin_terms(N, K) / 2>{}), ...);
    return b[0];
}

/**
 * @brief Estrin polynomial evaluation (float or double), unrolled at compile time
 *
 * (c₀ + x·c₁) + x²·(c₂ + x·c₃) + x⁴·(…) + …, level by level: b[i] = b[2i] + p·b[2i + 1] with
 * p = x, x², x⁴, …; N − 1 multiply-adds in ⌈log₂ N⌉ levels, every index a constant. Each sum
 * and power passes through assoc_barrier: -ffast-math would otherwise regroup the tree (into
 * Horner in x² over the pairs, for one) and the rounding would depend on the optimization level.
 */
template<typename T, size_t N>
constexpr T estrin_eval(T x, const T (&coeffs)[N]) {
    static_assert(N >= 1, "estrin_eval requires at least one coefficient");
    return estrin_levels(x, coeffs, damp::make_index_sequence<N>{}, damp::make_index_sequence<static_cast<size_t>(damp::bit_width(N - 1))>{});
}

/**
 * @brief Σ coeffs[First + Stride·j]·x^j by Horner, unrolled at compile time
 */
template<size_t First, size_t Stride, typename T, size_t N>
inline T horner_strided(T x, const T (&coeffs)[N]) {
    constexpr size_t count = (N - First + Stride - 1) / Stride;
    static_assert(count >= 1, "horner_strided requires at least one coefficient");
    return [&]<size_t... J>(damp::index_sequence<J...>) {
        T p = coeffs[First + (Stride * (count - 1))];
        ((p = coeffs[First + (Stride * (count - 2 - J))] + (x * p)), ...);
        return p;
    }(damp::make_index_sequence<count - 1>{});
}

/**
 * @brief Second-order Horner evaluation (float or double)
 *
 * Even and odd coefficients run as two independent Horner chains in x², joined by
 * one multiply-add: one operation more than Horner, half the dependency chain.
 */
template<typename T, size_t N>
inline T horner2_eval(T x, const T (&coeffs)[N]) {
    static_assert(N >= 2, "horner2_eval requires at least two coefficients");
    const T x2 = x * x;
    return horner_strided<0, 2>(x2, coeffs) + (x * horner_strided<1, 2>(x2, coeffs));
}

} // namespace damp::detail
