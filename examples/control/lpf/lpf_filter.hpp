// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file lpf_filter.hpp
 * @brief Low-pass filter design — flashable coefficients + objects
 *
 * Folder layout:
 *
 *   lpf_filter.hpp      — this file (design + filter objects)
 *   lpf_sketch.cpp      — setup/loop smoke
 *   lpf_sil.cpp         — host step-response print
 *   lpf_derivation.md   — filter model notes
 */

#pragma once

#include "damp/filters/filters.hpp"

namespace damp::examples_lpf {

inline constexpr float fc = 10.0f;    // cutoff [Hz]
inline constexpr float Ts = 0.001f;   // 1 kHz
inline constexpr float zeta = 0.707f; // 2nd-order damping

// design::lowpass_* → static constinit: coefficients and filter objects must be constexpr.
inline constexpr auto coeffs1 = design::lowpass_1st(fc, Ts);
inline constexpr auto coeffs2 = design::lowpass_2nd(fc, Ts, zeta);

/**
 * @brief One sample of the first-order LPF
 */
template<typename T>
[[nodiscard]] T filter_period(LowPass<1, T>& lpf, T x) {
    return lpf(x);
}

} // namespace damp::examples_lpf
