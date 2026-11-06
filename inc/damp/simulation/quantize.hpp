// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file quantize.hpp
 * @brief Uniform ADC quantizer for host SIL measurements
 *
 * Host-only SIL helper, reachable from @ref workbench.hpp (not @ref control.hpp).
 * A SIL reads the plant's continuous signals; a target reads ADC codes and scales
 * them back to engineering units. @ref damp::sim::quantize "quantize" gives the SIL
 * the second reading: the converter's resolution and rails.
 *
 * Sampling instant and hold belong to the harness: @c simulate_sampled and
 * @c simulate_multirate sample the plant at each tick and hold the command.
 */

#include "damp/backend.hpp"
#include "damp/math/math.hpp"

namespace damp::sim {

/**
 * @brief ADC reading of x: x rounded to the nearest multiple of lsb, then clamped to [lo, hi]
 *
 * @f[
 *   q(x) = \mathrm{clamp}\big(\mathrm{lsb}\cdot\mathrm{round}(x/\mathrm{lsb}),\ lo,\ hi\big)
 * @f]
 *
 * Ties round to even: an @f$x/\mathrm{lsb}@f$ exactly halfway between two integers goes
 * to the even one (IEEE round-half-to-even, @ref damp::nearbyint), so ties have no
 * bias. Rounding is odd-symmetric: @f$q(-x) = -q(x)@f$ between the rails.
 *
 * An N-bit bidirectional converter with full scale @f$\pm FS@f$ and zero at mid-scale has
 * @f$\mathrm{lsb} = 2\,FS/2^N@f$ and codes @f$-2^{N-1} \ldots 2^{N-1}-1@f$, so
 * @f$lo = -2^{N-1}\,\mathrm{lsb}@f$ and @f$hi = (2^{N-1}-1)\,\mathrm{lsb}@f$. Rails that are
 * multiples of lsb keep every reading on the code grid.
 *
 * Gain and offset errors are applied by the caller: @c quantize(g·x + o, lsb, lo, hi)
 * reads a channel with gain @c g and offset @c o, so the calibration the firmware
 * applies stays with the SIL.
 *
 * lsb ≤ 0 skips the rounding: the reading is x clamped to the rails.
 *
 * @code
 * // 12-bit phase-current sense, ±30 A full scale
 * constexpr double lsb = 60.0 / 4096.0; // 14.6 mA
 * const double     i_adc = sim::quantize(i_plant, lsb, -2048.0 * lsb, 2047.0 * lsb);
 * @endcode
 *
 * @tparam T  Floating-point scalar
 * @param x   Continuous value [plant units]
 * @param lsb One code [plant units]
 * @param lo  Bottom rail, the lowest reading
 * @param hi  Top rail, the highest reading (hi ≥ lo)
 * @return x quantized to lsb and clamped to [lo, hi]
 *
 * @note Compare with MATLAB®'s Simulink® Quantizer block (q·round(u/q), ties away from zero)
 *       followed by a Saturation block.
 * @see "Digital Control of Dynamic Systems" (Franklin, Powell & Workman, 1998), ch. 10
 */
template<typename T>
[[nodiscard]] constexpr T quantize(T x, T lsb, T lo, T hi) {
    if (!(lsb > T{0})) {
        return damp::clamp(x, lo, hi);
    }
    return damp::clamp(lsb * damp::nearbyint(x / lsb), lo, hi);
}

} // namespace damp::sim
