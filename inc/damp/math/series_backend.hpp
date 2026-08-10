// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file series_backend.hpp
 * @brief Freestanding runtime math backend: routes scalar math to the constexpr
 *        series in constexpr_math.hpp, pulling no hosted headers (no \<cmath\>).
 *
 * Selected with DAMP_MATH_BACKEND_FREESTANDING. Runtime accuracy is the same
 * series used at compile time (machine-precision over the practical range). The
 * natural pairing for the ETL container backend on targets without libm.
 */

#include "constexpr_math.hpp" // damp::detail::* (series; freestanding)
#include "math_backend.hpp"   // MathBackend<T> primary

namespace damp {

/**
 * @brief Math backend that evaluates the constexpr series at runtime.
 * @tparam T Scalar type (float, double)
 */
template<typename T>
struct SeriesMathBackend {
    static constexpr T sin(T x) { return detail::sin(x); }
    static constexpr T cos(T x) { return detail::cos(x); }
    static constexpr T tan(T x) { return detail::tan(x); }
    static constexpr T asin(T x) { return detail::asin(x); }
    static constexpr T acos(T x) { return detail::acos(x); }
    static constexpr T atan(T x) { return detail::atan(x); }
    static constexpr T atan2(T y, T x) { return detail::atan2(y, x); }

    static constexpr damp::pair<T, T> sincos(T x) { return detail::sincos(x); }

    static constexpr T sqrt(T x) { return detail::sqrt(x); }
    static constexpr T abs(T x) { return detail::abs(x); }
    static constexpr T cbrt(T x) { return detail::cbrt(x); }
    static constexpr T exp(T x) { return detail::exp(x); }
    static constexpr T log(T x) { return detail::log(x); }
    static constexpr T log10(T x) { return detail::log10(x); }
    static constexpr T floor(T x) { return detail::floor(x); }
    static constexpr T ceil(T x) { return detail::ceil(x); }
    static constexpr T nearbyint(T x) { return detail::nearbyint(x); }

    static constexpr T    pow(T base, T exponent) { return detail::pow(base, exponent); }
    static constexpr T    fmod(T x, T y) { return (y == T{0}) ? T{0} : detail::fmod(x, y); }
    static constexpr T    copysign(T mag, T sgn) { return detail::copysign(mag, sgn); }
    static constexpr bool isfinite(T x) { return detail::isfinite(x); }
};

template<>
struct MathBackend<float> : SeriesMathBackend<float> {};

template<>
struct MathBackend<double> : SeriesMathBackend<double> {};

} // namespace damp
