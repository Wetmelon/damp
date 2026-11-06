// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file options.hpp
 * @brief Trailing design options recognized by their type: at most one of each kind, in any
 *        order, each with its default when absent.
 *
 * A design call takes its model features as a trailing pack, `fn(sys, …, opts...)`. Each kind
 * of option is a trait, `Is<T>::value`, true for the types of that kind (every KalmanNoise
 * for is_kalman_noise). option_set<Kinds...>::check<Opts...>() refuses at compile time a type of
 * no listed kind and a kind given twice; option<Is>(fallback, opts...) returns the option of
 * that kind, or the fallback, with its own type, so the sizes it carries (a count of listed
 * outputs) are template arguments the design reads.
 *
 * @code
 * template<typename... Opts>
 * constexpr auto fn(const Sys& sys, const Opts&... opts) {
 *     damp::detail::option_set<design::detail::is_kalman_noise>::check<Opts...>();
 *     const auto noise = damp::detail::option<design::detail::is_kalman_noise>(design::KalmanNoise<NX, NY>{}, opts...);
 *     // noise is the KalmanNoise given, or the fallback without one
 * }
 * @endcode
 */

#include <cstddef>
#include <type_traits>


namespace damp::detail {

/** @brief How many of @p Opts are of the kind @p Is. */
template<template<typename> class Is, typename... Opts>
inline constexpr size_t option_count = (size_t{0} + ... + (Is<std::remove_cvref_t<Opts>>::value ? size_t{1} : size_t{0}));

/** @brief The type of the first option of the kind @p Is among @p Opts, else @p Default. */
template<template<typename> class Is, typename Default, typename... Opts>
struct option_type {
    using type = Default; ///< No option of the kind: the default
};

/** @brief The first option, if it is of the kind, else the rest's. */
template<template<typename> class Is, typename Default, typename First, typename... Rest>
struct option_type<Is, Default, First, Rest...> {
    using type = std::conditional_t<Is<std::remove_cvref_t<First>>::value, std::remove_cvref_t<First>, typename option_type<Is, Default, Rest...>::type>; ///< The option's type
};

/** @brief The type of the option of the kind @p Is among @p Opts, else @p Default. */
template<template<typename> class Is, typename Default, typename... Opts>
using option_t = typename option_type<Is, Default, Opts...>::type;

/** @brief No option left: the fallback. */
template<template<typename> class Is, typename Default>
[[nodiscard]] constexpr Default option(const Default& fallback) {
    return fallback;
}

/**
 * @brief The option of the kind @p Is among @p first, @p rest…, or @p fallback when none is.
 * @return A copy of the option, with its own type (option_t), so its sizes are known
 */
template<template<typename> class Is, typename Default, typename First, typename... Rest>
[[nodiscard]] constexpr option_t<Is, Default, First, Rest...> option(const Default& fallback, const First& first, const Rest&... rest) {
    if constexpr (Is<std::remove_cvref_t<First>>::value) {
        static_cast<void>(fallback);
        (static_cast<void>(rest), ...);
        return first;
    } else {
        static_cast<void>(first);
        return option<Is>(fallback, rest...);
    }
}

/**
 * @brief The kinds of option a design call accepts: check<Opts...>() refuses, at compile time, an
 *        option of no listed kind and a kind given twice.
 */
template<template<typename> class... Kinds>
struct option_set {
    /** @brief How many listed kinds @p Opt is of: 1 for a known option. */
    template<typename Opt>
    static constexpr size_t kinds_of = (size_t{0} + ... + (Kinds<std::remove_cvref_t<Opt>>::value ? size_t{1} : size_t{0}));

    /** @brief Every option is of exactly one listed kind. */
    template<typename... Opts>
    static constexpr bool known = ((kinds_of<Opts> == 1) && ...);

    /** @brief No kind is given twice. */
    template<typename... Opts>
    static constexpr bool unique = ((option_count<Kinds, Opts...> <= 1) && ...);

    /** @brief The compile-time check of a call's options. */
    template<typename... Opts>
    static constexpr void check() {
        static_assert(known<Opts...>, "An option of a type this design does not take: the trailing options are recognized by type.");
        static_assert(unique<Opts...>, "An option given twice: each kind of option is given at most once.");
    }
};

} // namespace damp::detail

namespace damp::design {

/** @brief The noise of a Kalman design (kalman.hpp), declared here so option sets can recognize it. */
template<size_t NW, size_t NV, typename T, size_t ND>
struct KalmanNoise;

namespace detail {

/** @brief Not a KalmanNoise. */
template<typename>
struct is_kalman_noise : std::false_type {};

/** @brief KalmanNoise. */
template<size_t NW, size_t NV, typename T, size_t ND>
struct is_kalman_noise<KalmanNoise<NW, NV, T, ND>> : std::true_type {};

} // namespace detail

} // namespace damp::design
