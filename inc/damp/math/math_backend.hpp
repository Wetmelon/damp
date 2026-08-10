// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file math_backend.hpp
 * @brief Pluggable runtime math backend selection (freestanding-safe).
 *
 * Declares the MathBackend\<T\> customization point and binds the chosen backend
 * implementation from the profile macros (see config.hpp). This header is
 * intentionally free of \<cmath\> so a freestanding build never pulls it — the
 * std:: implementation lives in std_fallback.hpp.
 */

#include "damp/config.hpp" // pulls the profile's backend-selection macros

namespace damp {

/**
 * @brief Pluggable math backend for runtime scalar operations
 *
 * Primary template is intentionally undefined — a backend specialization must be
 * provided for each scalar type used at runtime. Selection is driven by the
 * profile macros read through damp/config.hpp (set in your damp_profile.hpp):
 *
 *   - default                       → std_fallback.hpp (the std:: backend)
 *   - DAMP_MATH_BACKEND_DAMP          → damp_backend.hpp (fast float math/trig.hpp)
 *   - DAMP_MATH_BACKEND_FREESTANDING → series_backend.hpp (constexpr series, no \<cmath\>)
 *   - DAMP_MATH_BACKEND_HEADER "h"   → include "h", which defines MathBackend\<T\>
 *
 * Backend authors: inherit from StdMathFallback\<T\> (std_fallback.hpp) and
 * override only the functions your platform library provides; the rest fall
 * through to \<cmath\>. A freestanding backend instead routes to the constexpr
 * series (see series_backend.hpp) and pulls no hosted headers.
 *
 * @see std_fallback.hpp for StdMathFallback
 * @see config.hpp for the unified profile / macro surface
 * @tparam T Scalar type (float, double)
 */
template<typename T>
struct MathBackend;

} // namespace damp

// Bind the runtime backend selected by the profile macros (see config.hpp),
// after the MathBackend primary above is declared. Default: std::.
#if defined(DAMP_MATH_BACKEND_HEADER)
#include DAMP_MATH_BACKEND_HEADER // IWYU pragma: keep
#elif defined(DAMP_MATH_BACKEND_FREESTANDING)
#include "series_backend.hpp" // IWYU pragma: keep
#elif defined(DAMP_MATH_BACKEND_DAMP)
#include "damp_backend.hpp" // IWYU pragma: keep
#else
// Default: the std:: (\<cmath\>) backend — StdMathFallback unmodified.
#include "std_fallback.hpp" // IWYU pragma: keep
namespace damp {

template<>
struct MathBackend<float> : StdMathFallback<float> {};

template<>
struct MathBackend<double> : StdMathFallback<double> {};

} // namespace damp
#endif
