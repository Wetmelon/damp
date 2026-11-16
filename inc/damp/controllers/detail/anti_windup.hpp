// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file anti_windup.hpp
 * @brief Conditional integration against a clamp: the anti-windup the scalar laws with an
 *        integral share when no other is set
 */

#include "damp/backend.hpp"

namespace damp {
namespace detail {

/**
 * @brief Take back from @p state what its last step added, as far as a clamp cut the command
 *        off, and no further
 *
 * @p excess is u_sat − u_asked, negative when the command was clamped from above. A step that
 * pushed the command further past that limit is removed up to the excess, so the state stops
 * where the command meets the limit and never gives up what it held before the step: a jump of
 * another term through the limit does not reset it, as a short tracking time would. The state
 * enters the command with unit gain, so both are in the command's units.
 *
 * @param state  The integral
 * @param step   What the last tick added to it, reduced here by what is taken back: reports of
 *               a clamp take back no more than the step in all, but each takes up to its own
 *               excess, so a clamp is reported once
 * @param excess u_sat − u_asked
 *
 * @see "Advanced PID Control" (Åström & Hägglund, 2006), §3.5: conditional integration
 */
template<typename T>
constexpr void hold_against_clamp(T& state, T& step, T excess) {
    if (excess < T{0} && step > T{0}) {
        const T back = damp::min(step, -excess);
        state -= back;
        step -= back;
    } else if (excess > T{0} && step < T{0}) {
        const T back = damp::max(step, -excess);
        state -= back;
        step -= back;
    }
}

} // namespace detail
} // namespace damp
