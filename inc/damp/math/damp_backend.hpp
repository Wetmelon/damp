// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file damp_backend.hpp
 * @brief Fast float/double math backend (trig.hpp) over the std:: fallback.
 *
 * Routes MathBackend to damp::detail::fast_* kernels from trig.hpp. Those
 * kernels must not live as damp::sin/sqrt overloads — that would shadow the
 * constexpr public API in math.hpp and break compile-time design
 * (e.g. DARE → Cholesky → sqrt under DAMP_MATH_BACKEND_DAMP).
 */

// Not self-contained. Bound by math_backend.hpp after MathBackend\<T\> is declared.
#include "math_backend.hpp"
#include "std_fallback.hpp"
#include "trig.hpp"

namespace damp {

template<>
struct MathBackend<float> : StdMathFallback<float> {
    static float                    sin(float x) { return detail::fast_sin(x); }
    static float                    cos(float x) { return detail::fast_cos(x); }
    static float                    asin(float x) { return detail::fast_asin(x); }
    static float                    acos(float x) { return detail::fast_acos(x); }
    static float                    atan(float x) { return detail::fast_atan(x); }
    static float                    atan2(float y, float x) { return detail::fast_atan2(y, x); }
    static damp::pair<float, float> sincos(float x) { return detail::fast_sincos(x); }

    static float sqrt(float x) { return detail::fast_sqrt(x); }
    static float nearbyint(float x) { return detail::fast_nearbyint(x); }
};

template<>
struct MathBackend<double> : StdMathFallback<double> {
    static double                     sin(double x) { return detail::fast_sin(x); }
    static double                     cos(double x) { return detail::fast_cos(x); }
    static double                     asin(double x) { return detail::fast_asin(x); }
    static double                     acos(double x) { return detail::fast_acos(x); }
    static double                     atan(double x) { return detail::fast_atan(x); }
    static double                     atan2(double y, double x) { return detail::fast_atan2(y, x); }
    static damp::pair<double, double> sincos(double x) { return detail::fast_sincos(x); }

    static double sqrt(double x) { return detail::fast_sqrt(x); }
    static double nearbyint(double x) { return detail::fast_nearbyint(x); }
};

} // namespace damp
