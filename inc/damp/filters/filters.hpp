// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file filters.hpp
 * @brief IIR/FIR filter design and runtime filters — umbrella
 *
 * - iir_design.hpp — design:: coefficient builders (shared First/SecondOrderCoeffs)
 * - lowpass.hpp, highpass.hpp, biquad.hpp — IIR runtimes
 * - fir.hpp — direct-form FIR + design::fir1 window method (#32)
 * - delay.hpp, moving_average.hpp, median.hpp, complementary.hpp — other blocks
 */

#include "damp/filters/biquad.hpp"         // IWYU pragma: export
#include "damp/filters/complementary.hpp"  // IWYU pragma: export
#include "damp/filters/delay.hpp"          // IWYU pragma: export
#include "damp/filters/fir.hpp"            // IWYU pragma: export  (#32 FIR + fir1)
#include "damp/filters/highpass.hpp"       // IWYU pragma: export
#include "damp/filters/iir_design.hpp"     // IWYU pragma: export
#include "damp/filters/lowpass.hpp"        // IWYU pragma: export
#include "damp/filters/median.hpp"         // IWYU pragma: export
#include "damp/filters/moving_average.hpp" // IWYU pragma: export
