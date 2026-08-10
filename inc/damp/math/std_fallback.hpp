// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file std_fallback.hpp
 * @brief Composable std:: (\<cmath\>) base for the hosted math backends.
 *
 * Hosted-only: this is the one place the math layer pulls \<cmath\>. The std and
 * Damp (fast-float) backends inherit StdMathFallback and override only what their
 * platform provides; the freestanding/series backend does not use it, so a
 * freestanding build never reaches \<cmath\>.
 */

#include <cmath>

#include "constexpr_math.hpp" // detail::isfinite (bit-pattern, -ffast-math safe)
#include "damp/backend.hpp"   // damp::pair (sincos return type)

namespace damp {

/**
 * @brief Composable std:: base for math backends
 *
 * Backends inherit from this struct and override only the functions their
 * platform library provides. Any function not overridden falls through to the
 * corresponding \<cmath\> implementation.
 *
 * @tparam T Scalar type (float, double)
 */
template<typename T>
struct StdMathFallback {
    static T sin(T x) { return std::sin(x); }
    static T cos(T x) { return std::cos(x); }
    static T tan(T x) { return std::tan(x); }
    static T asin(T x) { return std::asin(x); }
    static T acos(T x) { return std::acos(x); }
    static T atan(T x) { return std::atan(x); }
    static T atan2(T y, T x) { return std::atan2(y, x); }

    /// Combined sin/cos. Returns {sin(x), cos(x)}. Platform backends should
    /// override this with a shared-range-reduction implementation.
    static damp::pair<T, T> sincos(T x) { return {std::sin(x), std::cos(x)}; }

    static T sqrt(T x) { return std::sqrt(x); }
    static T abs(T x) { return std::abs(x); }
    static T cbrt(T x) { return std::cbrt(x); }
    static T exp(T x) { return std::exp(x); }
    static T log(T x) { return std::log(x); }
    static T log10(T x) { return std::log10(x); }
    static T floor(T x) { return std::floor(x); }
    static T ceil(T x) { return std::ceil(x); }
    static T nearbyint(T x) { return std::nearbyint(x); }

    static T pow(T base, T exponent) { return std::pow(base, exponent); }
    static T fmod(T x, T y) { return std::fmod(x, y); }
    static T copysign(T mag, T sgn) { return std::copysign(mag, sgn); }
    /// Bit-pattern IEEE test (not `std::isfinite` — that can be a no-op under -ffast-math).
    static bool isfinite(T x) { return detail::isfinite(x); }
};

} // namespace damp
