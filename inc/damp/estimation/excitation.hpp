// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file excitation.hpp
 * @brief Allocation-free excitation generators — umbrella
 *
 * One generator per header (design + runtime together), matching controllers/:
 * - excitation/chirp.hpp
 * - excitation/prbs.hpp
 * - excitation/step_train.hpp
 * - excitation/ramp.hpp
 * - excitation/multi_sine.hpp
 * - excitation/impulse.hpp
 * - excitation/step.hpp
 * - excitation/stepped_sine.hpp
 */

#include "damp/estimation/excitation/chirp.hpp"        // IWYU pragma: export
#include "damp/estimation/excitation/impulse.hpp"      // IWYU pragma: export
#include "damp/estimation/excitation/multi_sine.hpp"   // IWYU pragma: export
#include "damp/estimation/excitation/prbs.hpp"         // IWYU pragma: export
#include "damp/estimation/excitation/ramp.hpp"         // IWYU pragma: export
#include "damp/estimation/excitation/step.hpp"         // IWYU pragma: export
#include "damp/estimation/excitation/step_train.hpp"   // IWYU pragma: export
#include "damp/estimation/excitation/stepped_sine.hpp" // IWYU pragma: export
