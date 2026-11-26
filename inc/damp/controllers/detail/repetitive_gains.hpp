// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once
// IWYU pragma: private, include "damp/controllers/repetitive.hpp"

/**
 * @file repetitive_gains.hpp
 * @brief Where a RepetitiveController's tuning lives: the period, gain, robustness-filter taps and
 *        lead as members, or as the constants of a `constexpr` design's Law
 */

namespace damp::detail {

/**
 * @brief A RepetitiveController's tuning as a member: the configuration it was built from
 *
 * @tparam Config design::RepetitiveConfig in the runtime's scalar
 */
template<typename Config>
struct RepetitiveGainSet {
    static constexpr bool constant = false; ///< A member of the runtime

    Config config{}; ///< The period N, the gain k_rc, the taps of Q and the lead m
};

/**
 * @brief A `constexpr` repetitive design's tuning as constants of the type: what a
 *        RepetitiveController reads in place of its configuration member, in the runtime's scalar
 *
 * The member is static, so the Law stores none of it: N, k_rc, the taps of Q and m sit in flash
 * and fold into the tick, and the ring's index arithmetic is on a constant period.
 *
 * @tparam Res design::repetitive's or repetitive_binomial's result, `constexpr`
 * @tparam T   The runtime's scalar
 */
template<const auto& Res, typename T>
struct ConstantRepetitiveGains {
    static constexpr bool constant = true; ///< The tuning is the type's

    static constexpr auto config = Res.config.template as<T>(); ///< The period N, the gain k_rc, the taps of Q and the lead m
};

} // namespace damp::detail
