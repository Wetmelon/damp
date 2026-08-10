// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file constexpr_math.hpp
 * @brief Compile-time scalar math: series / Newton-Raphson implementations.
 *
 * These are the bodies that back the public damp:: math functions at compile
 * time (consteval design code, static_asserts). They are intentionally free of
 * any runtime dispatch or backend dependency — math.hpp routes to them via
 * std::is_constant_evaluated() and to MathBackend\<T\> at runtime.
 *
 * Everything here lives in namespace damp::detail and is not part of the public
 * API; call the dispatching wrappers in math.hpp instead.
 *
 * @see math.hpp for the public dispatch layer
 */

#include <bit>
#include <cstdint>
#include <limits>
#include <type_traits>

#include "damp/backend.hpp"

namespace damp::detail {

/**
 * @brief Round to nearest integer (ties away from zero), returned as long long.
 *
 * Used for argument reduction in the constexpr trig/exp paths. Avoids std::round
 * (not constexpr-usable for our purposes) and the slow, precision-losing
 * subtract-in-a-loop reduction. The caller is responsible for keeping the
 * quotient within the range of long long; for the elementary functions here the
 * inputs that matter for compile-time design code are well within that range.
 */
template<typename T>
constexpr long long lround_away(T x) {
    return static_cast<long long>(x >= T{0} ? x + static_cast<T>(0.5) : x - static_cast<T>(0.5));
}

/// |x|.
template<typename T>
constexpr T abs(T x) {
    return x >= T{0} ? x : -x;
}

/// Square root via Newton-Raphson. Returns 0 for x <= 0 (NaN is unavailable in
/// constant evaluation).
template<typename T>
constexpr T sqrt(T x) {
    if (x == T{0}) {
        return T{0};
    }
    if (x < T{0}) {
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

/// Two-argument arctangent ∈ [−π, π] via Taylor series with three-interval range
/// reduction. @see Cody & Waite, "Software Manual for the Elementary Functions".
template<typename T>
constexpr T atan2(T y, T x) {
    constexpr T pi = damp::numbers::pi_v<T>;

    if (x == T{0}) {
        if (y > T{0}) {
            return pi / T{2};
        }
        if (y < T{0}) {
            return -pi / T{2};
        }
        return T{0};
    }

    T ratio = y / x;
    T atan_val;

    T abs_ratio = ratio >= T{0} ? ratio : -ratio;

    if (abs_ratio <= static_cast<T>(0.4142135623730951)) {
        // |t| <= tan(π/8): Taylor series directly.
        T r2 = ratio * ratio;
        T term = ratio;
        atan_val = term;
        for (int n = 1; n <= 15; ++n) {
            term *= -r2;
            atan_val += term / T((2 * n) + 1);
        }
    } else if (abs_ratio <= static_cast<T>(2.4142135623730951)) {
        // tan(π/8) < |t| <= tan(3π/8): atan(t) = π/4 + atan((t−1)/(t+1)).
        T reduced = (abs_ratio - T{1}) / (abs_ratio + T{1});
        T r2 = reduced * reduced;
        T term = reduced;
        T atan_reduced = term;
        for (int n = 1; n <= 15; ++n) {
            term *= -r2;
            atan_reduced += term / T((2 * n) + 1);
        }
        atan_val = (pi / T{4}) + atan_reduced;
        if (ratio < T{0}) {
            atan_val = -atan_val;
        }
    } else {
        // |t| > tan(3π/8): atan(t) = π/2 − atan(1/t).
        T inv = T{1} / abs_ratio;
        T r2 = inv * inv;
        T term = inv;
        T atan_inv = term;
        for (int n = 1; n <= 15; ++n) {
            term *= -r2;
            atan_inv += term / T((2 * n) + 1);
        }
        atan_val = (pi / T{2}) - atan_inv;
        if (ratio < T{0}) {
            atan_val = -atan_val;
        }
    }

    // Adjust for quadrant.
    if (x < T{0}) {
        atan_val += (y >= T{0} ? pi : -pi);
    }

    return atan_val;
}

/// Single-argument arctangent ∈ (−π/2, π/2).
template<typename T>
constexpr T atan(T x) {
    return atan2(x, T{1});
}

/// Arcsine core, assumes |x| < 1 (domain clamping is done by the dispatcher).
/// sqrt((1−x)(1+x)) avoids catastrophic cancellation near |x| = 1.
template<typename T>
constexpr T asin(T x) {
    return atan2(x, sqrt((T{1} - x) * (T{1} + x)));
}

/// Arccosine core, assumes |x| < 1 (domain clamping is done by the dispatcher).
template<typename T>
constexpr T acos(T x) {
    return atan2(sqrt((T{1} - x) * (T{1} + x)), x);
}

/// Cosine via Taylor series, range-reduced to [−π, π] with a two-part 2π.
template<typename T>
constexpr T cos(T x) {
    constexpr T two_pi_hi = static_cast<T>(6.28318530693650245668);
    constexpr T two_pi_lo = static_cast<T>(2.43084020260247689728e-10);
    constexpr T inv_two_pi = static_cast<T>(0.15915494309189533577);
    const T     kreal = static_cast<T>(lround_away(x * inv_two_pi));
    x = (x - (kreal * two_pi_hi)) - (kreal * two_pi_lo);

    T x2 = x * x;
    T result = T{1};
    T term = T{1};
    for (int n = 1; n <= 12; ++n) {
        term *= -x2 / T(2 * n * ((2 * n) - 1));
        result += term;
    }
    return result;
}

/// Sine via Taylor series, range-reduced to [−π, π] with a two-part 2π.
template<typename T>
constexpr T sin(T x) {
    constexpr T two_pi_hi = static_cast<T>(6.28318530693650245668);
    constexpr T two_pi_lo = static_cast<T>(2.43084020260247689728e-10);
    constexpr T inv_two_pi = static_cast<T>(0.15915494309189533577);
    const T     kreal = static_cast<T>(lround_away(x * inv_two_pi));
    x = (x - (kreal * two_pi_hi)) - (kreal * two_pi_lo);

    T x2 = x * x;
    T result = x;
    T term = x;
    for (int n = 1; n <= 12; ++n) {
        term *= -x2 / T((2 * n) * ((2 * n) + 1));
        result += term;
    }
    return result;
}

/// {sin(x), cos(x)}.
template<typename T>
constexpr damp::pair<T, T> sincos(T x) {
    return {sin(x), cos(x)};
}

/// Tangent via continued fraction, range-reduced to [−π/2, π/2] with a two-part
/// π and a complementary-angle identity near ±π/2.
/// @see Cuyt et al., "Handbook of Continued Fractions for Special Functions" §12.1
template<typename T>
constexpr T tan(T x) {
    constexpr T pi = damp::numbers::pi_v<T>;
    constexpr T half_pi = pi / T{2};

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

/// exp via ln2 argument reduction (exp(x) = 2^k · exp(r), |r| ≤ ln2/2) and a
/// Taylor series on the remainder. Saturates over/underflow to max() / 0 (the
/// library builds under -ffinite-math-only, so ±inf must never be produced).
template<typename T>
constexpr T exp(T x) {
    if (x > static_cast<T>(709.782712893384)) {
        return std::numeric_limits<T>::max();
    }
    if (x < static_cast<T>(-745.133219101941)) {
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

/// Natural log via Newton-Raphson on eʸ = x, with argument reduction by powers
/// of e. Returns 0 for x <= 0.
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

/// base^exponent = exp(exponent · ln(base)). Returns 1 for exponent 0, 0 for
/// base <= 0.
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

/// Largest integer <= x.
template<typename T>
constexpr T floor(T x) {
    T int_part = static_cast<long long>(x);
    if (x < T{0} && x != int_part) {
        int_part -= T{1};
    }
    return int_part;
}

/// Smallest integer >= x.
template<typename T>
constexpr T ceil(T x) {
    T int_part = static_cast<long long>(x);
    if (x > T{0} && x != int_part) {
        int_part += T{1};
    }
    return int_part;
}

/// Round to nearest integer (ties away from zero). Differs from std::nearbyint
/// only in tie-breaking, which is immaterial for angle range reduction.
template<typename T>
constexpr T nearbyint(T x) {
    return x >= T{0} ? floor(x + static_cast<T>(0.5)) : ceil(x - static_cast<T>(0.5));
}

/// Floating-point remainder core (truncated-quotient convention), assumes
/// y != 0 (the zero guard is done by the dispatcher).
template<typename T>
constexpr T fmod(T x, T y) {
    const T q = x / y;
    const T truncated = static_cast<T>(static_cast<long long>(q)); // toward zero
    return x - (truncated * y);
}

/// log10(x) = ln(x) / ln(10). Returns 0 for x <= 0.
template<typename T>
constexpr T log10(T x) {
    if (x <= T{0}) {
        return T{0};
    }
    return log(x) / log(T{10});
}

/// Magnitude of @p mag with the sign of @p sgn_src.
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
