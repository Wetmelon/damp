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

#include "block.hpp"         // IWYU pragma: export
#include "colvec.hpp"        // IWYU pragma: export
#include "core.hpp"          // IWYU pragma: export
#include "decomposition.hpp" // IWYU pragma: export
#include "eigen.hpp"         // IWYU pragma: export
#include "functions.hpp"     // IWYU pragma: export
#include "rowvec.hpp"        // IWYU pragma: export
#include "solve.hpp"         // IWYU pragma: export
#include "svd.hpp"           // IWYU pragma: export
#include "views.hpp"         // IWYU pragma: export
