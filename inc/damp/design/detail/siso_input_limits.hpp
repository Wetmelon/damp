// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
// IWYU pragma: private

/**
 * @file siso_input_limits.hpp
 * @brief The options a one-input law's design reads, an InputLimits on its one command: the check
 *        that refuses anything else at compile time, and the box in the design's scalar
 *
 * design::smc, design::stsmc, design::adrc and design::esc take `InputLimits{.lower = {lo},
 * .upper = {hi}}` among their trailing options, as design::pid does, and nothing else.
 */

#include "damp/backend.hpp"
#include "damp/design/input_limits.hpp"
#include "damp/design/options.hpp"
#include "damp/toolbox/bounds.hpp"

namespace damp::design::detail {

/**
 * @brief The compile-time check of a one-input law's options: an InputLimits at most, with one
 *        entry at each end given
 *
 * An option of another kind, a second InputLimits or one on two inputs does not compile, each
 * with its message.
 */
template<typename... Opts>
consteval void check_siso_input_limits() {
    option_set<kind::InputLimits>::check<Opts...>();
    if constexpr (option_count<kind::InputLimits, Opts...> > 0) {
        static_assert(option_t<kind::InputLimits, void, Opts...>::inputs == 1, "This law has one input: InputLimits{.lower = {lo}, .upper = {hi}}, one entry at each end.");
    }
}

/** @brief The InputLimits among @p opts on the law's one command, in its scalar @p T (none: unbounded). */
template<typename T, typename... Opts>
[[nodiscard]] constexpr Bounds<1, T> siso_input_limits_of(const Opts&... opts) {
    if constexpr (option_count<kind::InputLimits, Opts...> > 0) {
        if constexpr (option_t<kind::InputLimits, void, Opts...>::inputs == 1) { // more is refused by check_siso_input_limits
            return input_limits_of<1, T>(opts...);
        }
    }
    (static_cast<void>(opts), ...);
    return Bounds<1, T>{};
}

} // namespace damp::design::detail
