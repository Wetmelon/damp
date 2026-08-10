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
 *
 *   Scalar math (see math/math_backend.hpp):
 *     DAMP_MATH_BACKEND_DAMP          runtime float math → damp/math/trig.hpp
 *     DAMP_MATH_BACKEND_HEADER "h"   include "h", which defines MathBackend\<T\>
 *     (neither set → the std:: math backend)
 *
 * If no `damp_profile.hpp` is found, the library uses its host defaults — stdlib
 * containers and the std:: math backend — and emits a one-time warning so the
 * choice is never silent.
 *
 * @code
 * // damp_profile.hpp  (user-created, anywhere on the include path) — macros only
 * #define DAMP_BACKEND_ETL                              // ETL containers
 * #define DAMP_MATH_BACKEND_HEADER "my_cmsis_math.hpp"  // custom math backend
 * @endcode
 */

#if __has_include("damp_profile.hpp")
#include "damp_profile.hpp" // IWYU pragma: keep
#else
#warning "damp_profile.hpp not found in include path. Using host defaults (stdlib containers + std:: math). Create a damp_profile.hpp to select platform backends."
#endif
