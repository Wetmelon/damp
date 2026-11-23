// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file matrix.hpp
 * @brief Umbrella header for the fixed-size linear-algebra library. Pulls in the
 *        Matrix class (core.hpp) plus every extension that provides out-of-line
 *        members and free functions (block views, vectors, factorizations,
 *        solvers, eigen, matrix functions). Include this to get the full API.
 *        The leaves each include core.hpp directly, so include order here is
 *        irrelevant and the formatter may sort these freely.
 */

#include "damp/matrix/block.hpp"             // IWYU pragma: export
#include "damp/matrix/colvec.hpp"            // IWYU pragma: export
#include "damp/matrix/concat.hpp"            // IWYU pragma: export
#include "damp/matrix/const_matrix.hpp"      // IWYU pragma: export
#include "damp/matrix/core.hpp"              // IWYU pragma: export
#include "damp/matrix/decomposition.hpp"     // IWYU pragma: export
#include "damp/matrix/eigen.hpp"             // IWYU pragma: export
#include "damp/matrix/functions.hpp"         // IWYU pragma: export
#include "damp/matrix/generalized_eigen.hpp" // IWYU pragma: export
#include "damp/matrix/matrix_traits.hpp"     // IWYU pragma: export
#include "damp/matrix/pattern_matrix.hpp"    // IWYU pragma: export
#include "damp/matrix/rowvec.hpp"            // IWYU pragma: export
#include "damp/matrix/solve.hpp"             // IWYU pragma: export
#include "damp/matrix/svd.hpp"               // IWYU pragma: export
#include "damp/matrix/views.hpp"             // IWYU pragma: export
