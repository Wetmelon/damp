// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file filters.hpp
 * @brief Filters: the design:: functions and the runtimes that run them — umbrella
 *
 * Design at compile time or init, run the runtime each sample: `LowPass lpf{design::lowpass<2>(fc, Ts)}`,
 * then `y = lpf(x)`. Every runtime built from a failed design refuses it under constant evaluation.
 *
 * | Need | Design | Runtime | Header |
 * | ---- | ------ | ------- | ------ |
 * | Low-pass, order 1–4 | `design::lowpass<N>(fc, Ts)` | `LowPass` | lowpass.hpp |
 * | Low-pass, any order (Butterworth) | `design::lowpass_sections<N>(fc, Ts)` | `BiquadCascade` | biquad.hpp |
 * | High-pass / washout or low-pass, order 1 | `design::highpass_1st(fc, Ts)`, `design::lowpass_1st(fc, Ts)` | `FirstOrderSection` | first_order_section.hpp |
 * | High-pass, any order (Butterworth) | `design::highpass_sections<N>(fc, Ts)` | `BiquadCascade` | biquad.hpp |
 * | One second-order shape | `design::lowpass_2nd`, `highpass_2nd`, `bandpass`, `notch`, `allpass`, `peaking`, `lowshelf`, `highshelf` | `Biquad` | iir_design.hpp |
 * | Any transfer function in z, order ≥ 3 | `design::second_order_sections(tf)` | `BiquadCascade` | biquad.hpp |
 * | Fractional-sample delay (Thiran) | `fractional_delay<N>(τ, Ts)`, then `design::second_order_sections` | `BiquadCascade` | systems/discretization.hpp |
 * | Whole-sample delay | | `Delay` | delay.hpp |
 * | Rational delay model e^(−sT) (Padé), continuous or Tustin | `design::pade_delay_1st`, `design::pade_delay_2nd` | | pade.hpp |
 * | FIR, window method | `design::fir1`, `design::fir1_hz` | `FirFilter` | fir.hpp |
 * | Adaptive FIR | | `NlmsFilter` | fir.hpp |
 * | Moving average, comb notch | | `MovingAverage` | moving_average.hpp |
 * | Spike rejection | | `MedianFilter` | median.hpp |
 * | Blend a fast and a slow sensor | | `Complementary` | complementary.hpp |
 * | Derivative of a noisy signal | | `RobustExactDifferentiator` | differentiator.hpp |
 * | Quadrature signals, resonant band-pass | `design::sogi`, `design::mstogi` | `SOGI`, `MSTOGI`, `SogiFll` | sogi.hpp |
 * | One frequency's amplitude and phase | | `Goertzel`, `HarmonicAnalyzer` | spectral.hpp |
 *
 * Above second order, run a filter as sections rather than one direct form: rounding in a
 * direct form's coefficients moves poles near z = 1 far more than it moves them in a section.
 * The phase-locked loops (pll.hpp) close a PI loop, so they are exported by control.hpp.
 */

#include "damp/filters/biquad.hpp"              // IWYU pragma: export
#include "damp/filters/complementary.hpp"       // IWYU pragma: export
#include "damp/filters/delay.hpp"               // IWYU pragma: export
#include "damp/filters/differentiator.hpp"      // IWYU pragma: export
#include "damp/filters/fir.hpp"                 // IWYU pragma: export
#include "damp/filters/first_order_section.hpp" // IWYU pragma: export
#include "damp/filters/iir_design.hpp"          // IWYU pragma: export
#include "damp/filters/lowpass.hpp"             // IWYU pragma: export
#include "damp/filters/median.hpp"              // IWYU pragma: export
#include "damp/filters/moving_average.hpp"      // IWYU pragma: export
#include "damp/filters/pade.hpp"                // IWYU pragma: export
#include "damp/filters/sogi.hpp"                // IWYU pragma: export
#include "damp/filters/spectral.hpp"            // IWYU pragma: export
