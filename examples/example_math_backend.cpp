// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file example_math_backend.cpp
 * @brief Pluggable math backend example
 *
 * Demonstrates how to select and extend a math backend for embedded deployment.
 * The library's damp:: math functions (sin, cos, sqrt, etc.) dispatch through
 * MathBackend<T> at runtime. By providing a damp_profile.hpp anywhere in your
 * project's include path, you route math operations to hardware-accelerated
 * implementations without changing any library or application code.
 *
 * The compile-time path (consteval design functions, static_assert checks)
 * always uses the built-in constexpr series expansions regardless of backend.
 * Backends only affect runtime execution.
 *
 * Backend selection is done via a user-created, macro-only damp_profile.hpp (see
 * damp/config.hpp for the full surface):
 *
 * @code
 * // damp_profile.hpp  — place anywhere in your project's include path (macros only)
 *
 * // Option 1: std:: backend (desktop / host builds — the default, no macro needed)
 *
 * // Option 2: built-in fast float math (damp/math/trig.hpp)
 * #define DAMP_MATH_BACKEND_DAMP
 *
 * // Option 3: custom platform backend — point at a header that defines it:
 * #define DAMP_MATH_BACKEND_HEADER "my_platform_math.hpp"
 *
 * // ...where my_platform_math.hpp inherits StdMathFallback<T> and overrides only
 * // the functions the platform provides:
 * namespace damp {
 * template<>
 * struct MathBackend<float> : StdMathFallback<float> {
 *     static float sin(float x)  { return my_platform_sinf(x); }
 *     static float cos(float x)  { return my_platform_cosf(x); }
 *     static float sqrt(float x) { return my_platform_sqrtf(x); }
 *     // Unoverridden functions (cbrt, exp, log, etc.) fall through to std::
 * };
 * } // namespace damp
 * @endcode
 *
 * If no damp_profile.hpp is found, the library uses the std:: backend and emits a
 * one-time #warning so the choice is never silent.
 */

// This example uses the std:: backend via examples/damp_profile.hpp.
#include <cstdio>
#include <numbers>

#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"

using namespace damp;

int main() {
    // Runtime transcendentals go through MathBackend<float>. With
    // DAMP_MATH_BACKEND_DAMP (or a custom DAMP_MATH_BACKEND_HEADER), runtime math
    // uses that backend; this example's damp_profile keeps the default std:: path.
    constexpr float theta = std::numbers::pi_v<float> / 4.0f; // 45°

    const float s = damp::sin(theta);
    const float c = damp::cos(theta);
    const float r = damp::sqrt(s * s + c * c);

    // Small matrix product also exercises backend float arithmetic end-to-end.
    const ColVec<2, float> v{{s, c}};
    const float            n = v.norm();

    std::printf("\n=== MathBackend dispatch (damp::sin / cos / sqrt) ===\n");
    std::printf("theta = %.4f rad\n", double(theta));
    std::printf("sin   = %.4f\n", double(s));
    std::printf("cos   = %.4f\n", double(c));
    std::printf("hypot = %.4f  (expect 1)\n", double(r));
    std::printf("||v|| = %.4f  (expect 1)\n", double(n));

    return 0;
}
