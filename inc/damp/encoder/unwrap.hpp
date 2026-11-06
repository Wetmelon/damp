// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file damp/encoder/unwrap.hpp
 * @brief Single-turn (or single-period) absolute readings → exact multi-period position.
 *
 * An absolute angle sensor (magnetic encoder, resolver, sin/cos) reports a value that
 * wraps every period. Differencing two readings across that seam gives a full-scale jump;
 * the classic bug. PhaseUnwrapper folds each per-sample change onto the shortest arc and
 * keeps the result as a SplitReal, so it stays exact however many periods accumulate.
 * Feed it the phase as a fraction of a period: divide radians by 2π, a raw single-turn
 * count by the counts per turn.
 *
 * `constexpr`, allocation-free, embeddable.
 *
 * @see wrapped_delta (encoder/quadrature.hpp), the integer sibling for counter rollover.
 */

#include "damp/math/math.hpp"       // damp::wrap
#include "damp/math/split_real.hpp" // SplitReal

namespace damp {

/**
 * @brief Counts whole periods from successive fine phases (no hardware counter).
 *
 * Each update takes the shortest-arc step from the previous phase, so the count is
 * right only while the axis moves less than half a period between samples. Use it for
 * a resolver (a period is a whole turn); use @ref combine_quadrature for multi-period
 * encoders and scales.
 *
 * @tparam T Scalar type.
 */
template<typename T = float>
class PhaseUnwrapper {
public:
    constexpr PhaseUnwrapper() noexcept = default;

    /** @brief Start from a known absolute position (e.g. after homing). */
    constexpr explicit PhaseUnwrapper(SplitReal<T> start) noexcept : pos_(start) {}

    /** @brief Advance to the period containing @p frac [periods, in [−0.5, 0.5)]; returns the position. */
    constexpr SplitReal<T> update(T frac) {
        constexpr T half = static_cast<T>(0.5);
        pos_ += damp::wrap(frac - pos_.frac, -half, half);
        return pos_;
    }

    [[nodiscard]] constexpr SplitReal<T> position() const { return pos_; }

    constexpr void reset(SplitReal<T> start = {}) { pos_ = start; }

private:
    SplitReal<T> pos_{};
};

} // namespace damp
