// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file analysis.hpp
 * @brief Frequency-domain and time-domain LTI analysis (host-only) — umbrella
 *
 * @defgroup analysis System Analysis
 * @brief Host frequency- and time-domain LTI analysis (Bode, margins, poles, step/lsim)
 *
 * Host-only. Sweeps and time responses allocate @c std::vector; this is not
 * part of the embeddable @c control.hpp surface. Include via @c workbench.hpp
 * (or this header) on the workstation.
 *
 * Split for navigation; this header re-exports the full analysis surface:
 * - linspace.hpp — linspace / logspace / geomspace / arange
 * - frequency.hpp — Bode, Nyquist, Nichols, sigma, margins, loop metrics, impedance
 * - norms.hpp — dcgain, H2 / H∞
 * - poles.hpp — poles, damping, pzmap, rlocus
 * - time_response.hpp — step / impulse / lsim / stepinfo
 *
 * Controllability / observability rank tests live in damp/design/stability.hpp
 * (embeddable). FRF helpers omit contour poles honestly (@ref eval_frf
 * @c nullopt → sample skipped, no Inf/NaN fill).
 */

#include "damp/analysis/frequency.hpp"     // IWYU pragma: export
#include "damp/analysis/linspace.hpp"      // IWYU pragma: export
#include "damp/analysis/norms.hpp"         // IWYU pragma: export
#include "damp/analysis/poles.hpp"         // IWYU pragma: export
#include "damp/analysis/time_response.hpp" // IWYU pragma: export
