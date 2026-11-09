// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file config.hpp
 * @brief Single configuration entry point for the Damp library.
 *
 * The library reads all of its compile-time configuration from one place. A
 * user-authored `damp_profile.hpp`, discovered on the include path, sets
 * per-facility selection macros; each facility (containers in backend.hpp,
 * scalar math in math/math_backend.hpp) keys off those macros.
 *
 * `damp_profile.hpp` must define macros only — it must not include a backend
 * implementation. That keeps this header safe to read both early (backend.hpp,
 * before any Damp types exist) and late (math_backend.hpp, after its types are
 * declared) without include-ordering hazards.
 *
 * Recognized macros (all optional — omit for the host defaults):
 *
 *   Containers (see backend.hpp):
 *     DAMP_BACKEND_ETL               damp::array/optional/... → ETL (else stdlib)
 *     (ETL_NO_STL goes in etl_profile.h, which ETL reads itself: see "Using ETL" in README.md)
 *
 *   Scalar math (see math/math_backend.hpp):
 *     DAMP_MATH_BACKEND_DAMP          runtime float math → damp/math/trig.hpp
 *     DAMP_MATH_BACKEND_FREESTANDING  constexpr series at run time, no \<cmath\> (series_backend.hpp)
 *     DAMP_MATH_BACKEND_HEADER "h"   include "h", which defines MathBackend\<T\>
 *     (none set → the std:: math backend; several set → HEADER, then FREESTANDING, then DAMP)
 *
 *   A target without libm sets DAMP_BACKEND_ETL and DAMP_MATH_BACKEND_FREESTANDING: the
 *   std:: and Damp math backends include \<cmath\> whatever the container backend, so a
 *   toolchain with no C++ standard headers takes FREESTANDING or HEADER too.
 *
 * If no `damp_profile.hpp` is found and none of these macros is set, the library
 * uses its host defaults — stdlib containers and the std:: math backend — without
 * a warning. A macro set with no profile (e.g. `-DDAMP_BACKEND_ETL`) draws a
 * warning: every file in the program must agree on the backend, and a profile on
 * the include path is the one place they all read it from.
 *
 * @code
 * // damp_profile.hpp  (user-created, anywhere on the include path) — macros only
 * #define DAMP_BACKEND_ETL                              // ETL containers
 * #define DAMP_MATH_BACKEND_HEADER "my_cmsis_math.hpp"  // custom math backend
 * @endcode
 */

#if __has_include("damp_profile.hpp")
#include "damp_profile.hpp" // IWYU pragma: keep
#elif defined(DAMP_BACKEND_ETL) || defined(DAMP_MATH_BACKEND_DAMP) || defined(DAMP_MATH_BACKEND_FREESTANDING) || defined(DAMP_MATH_BACKEND_HEADER)
#warning "A Damp backend macro is set but damp_profile.hpp is not on the include path. Every file in the program must agree on the backend: set it in a damp_profile.hpp rather than with -D."
#endif
