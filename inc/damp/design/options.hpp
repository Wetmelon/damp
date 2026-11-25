// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file options.hpp
 * @brief Trailing design options recognized by their kind: at most one of each kind, in any
 *        order, each with its default when absent.
 *
 * A design call takes its model features as a trailing pack, `fn(sys, …, opts...)`. Each option
 * type names its kind, `using option_kind = design::kind::InputLimits;`, a tag of design::kind,
 * and OptionOf<O, Kind> tests it. option_set<Kinds...>::check<Opts...>() refuses at compile time
 * an option of no listed kind and a kind given twice, with that kind's message;
 * option<Kind>(fallback, opts...) returns the option of that kind, or the fallback, with its own
 * type, so the sizes it carries (a count of listed outputs) are template arguments the design
 * reads; the<Kind>(opts...) returns one that is there, inside an `if constexpr` that tested it.
 *
 * @code
 * template<typename... Opts>
 * constexpr auto fn(const Sys& sys, const Opts&... opts) {
 *     damp::detail::option_set<design::kind::KalmanNoise>::check<Opts...>();
 *     const auto noise = damp::detail::option<design::kind::KalmanNoise>(design::KalmanNoise<NX, NY>{}, opts...);
 *     // noise is the KalmanNoise given, or the fallback without one
 * }
 * @endcode
 */

#include "damp/backend.hpp"

namespace damp::design {

/** @brief The kinds of trailing design option: the tag each option type names as its option_kind. */
namespace kind {

struct InputLimits {};         ///< InputLimits: the command's limits
struct InputRateLimits {};     ///< InputRateLimits: the move's limits
struct DemandLimit {};         ///< DemandLimit: what an LQ gain's outer terms ask of a state
struct IntegralAction {};      ///< IntegralAction: integrators on r − y
struct TrackedOutputs {};      ///< TrackedOutputs: the outputs the target holds at r
struct Disturbance {};         ///< A disturbance spec: InputDisturbance, OutputDisturbance, PlantDisturbance, DisturbanceModel (or NoDisturbance)
struct MeasuredDisturbance {}; ///< MeasuredDisturbance: a disturbance read each tick
struct KalmanNoise {};         ///< KalmanNoise: a Kalman filter's noise
struct ObserverPoles {};       ///< ObserverPoles: a placed filter's poles
struct InputDelay {};          ///< InputDelay: the input arrives late
struct OutputDelay {};         ///< OutputDelay: the reading arrives late
struct Handover {};            ///< Handover: an LQ law's fade after track()
struct Weights {};             ///< LQRWeights or Poles: the weights, which are an argument, never an option
struct OutputLimits {};        ///< OutputLimits: an MPC's bounded outputs
struct StateLimits {};         ///< StateLimits: an MPC's bounded states
struct Terminal {};            ///< DareCost or InvariantSet: an MPC's terminal ingredient
struct QPIterations {};        ///< QPIterations: a QP's iteration budget

} // namespace kind

/** @brief An option of the kind @p Kind: a type that names it, `using option_kind = Kind;`. */
template<typename O, typename Kind>
concept OptionOf = requires { typename damp::remove_cvref_t<O>::option_kind; } && damp::is_same_v<typename damp::remove_cvref_t<O>::option_kind, Kind>;

} // namespace damp::design

namespace damp::detail {

/** @brief How many of @p Opts are of the kind @p Kind. */
template<typename Kind, typename... Opts>
inline constexpr size_t option_count = (size_t{0} + ... + (design::OptionOf<Opts, Kind> ? size_t{1} : size_t{0}));

/** @brief A type carried as a value, so that a fold over a call's options can pick one (option_t). */
template<typename Kind, typename T>
struct OptionPick {
    using type = T; ///< The type picked so far
};

/** @brief The left side when it is of the kind @p Kind, else the right: folded over the options onto a default, the first option of the kind. */
template<typename Kind, typename L, typename R>
[[nodiscard]] constexpr auto operator|(OptionPick<Kind, L> left, OptionPick<Kind, R> right) {
    if constexpr (design::OptionOf<L, Kind>) {
        return left;
    } else {
        return right;
    }
}

/** @brief The type of the first option of the kind @p Kind among @p Opts (a design takes at most one), else @p Default. */
template<typename Kind, typename Default, typename... Opts>
using option_t = typename decltype((OptionPick<Kind, damp::remove_cvref_t<Opts>>{} | ... | OptionPick<Kind, Default>{}))::type;

/**
 * @brief The option of the kind @p Kind among @p opts, which hold one: for code inside an
 *        `if constexpr` that tested it is there.
 * @return A copy of the first option of the kind, with its own type, so its sizes are known
 */
template<typename Kind, typename... Opts>
[[nodiscard]] constexpr option_t<Kind, void, Opts...> the(const Opts&... opts) {
    static_assert(option_count<Kind, Opts...> > 0, "the<Kind>(opts...) reads an option the call holds: test option_count<Kind, Opts...> first.");
    const option_t<Kind, void, Opts...>* found = nullptr;
    ([&] {
        if constexpr (design::OptionOf<Opts, Kind>) {
            found = (found == nullptr) ? &opts : found;
        }
    }(), ...);
    return *found;
}

/**
 * @brief The option of the kind @p Kind among @p opts, or @p fallback when none is.
 * @return A copy of the option, with its own type (option_t), so its sizes are known
 */
template<typename Kind, typename Default, typename... Opts>
[[nodiscard]] constexpr option_t<Kind, Default, Opts...> option(const Default& fallback, const Opts&... opts) {
    if constexpr (option_count<Kind, Opts...> == 0) {
        (static_cast<void>(opts), ...);
        return fallback;
    } else {
        static_cast<void>(fallback);
        return the<Kind>(opts...);
    }
}

/**
 * @brief Refuses, at compile time, @p N options of the kind @p Kind when that is more than one,
 *        with the kind's message: the spelling of the one that holds them all.
 */
template<typename Kind, size_t N>
consteval void at_most_one() {
    namespace kind = design::kind;
    if constexpr (damp::is_same_v<Kind, kind::InputLimits>) {
        static_assert(N <= 1, "One InputLimits bounds every input: InputLimits{.lower = {...}, .upper = {...}}.");
    } else if constexpr (damp::is_same_v<Kind, kind::InputRateLimits>) {
        static_assert(N <= 1, "One InputRateLimits bounds every input's move: InputRateLimits{.lower = {...}, .upper = {...}}.");
    } else if constexpr (damp::is_same_v<Kind, kind::OutputLimits>) {
        static_assert(N <= 1, "One OutputLimits lists every bounded output: OutputLimits{.outputs = {i, ...}, .lower = {...}, .upper = {...}}.");
    } else if constexpr (damp::is_same_v<Kind, kind::StateLimits>) {
        static_assert(N <= 1, "One StateLimits lists every bounded state: StateLimits{.states = {i, ...}, .lower = {...}, .upper = {...}}.");
    } else if constexpr (damp::is_same_v<Kind, kind::Terminal>) {
        static_assert(N <= 1, "One terminal ingredient: DareCost{} for the cost, or InvariantSet<N>{.target_margin = ...} for the cost and the set, not both.");
    } else if constexpr (damp::is_same_v<Kind, kind::QPIterations>) {
        static_assert(N <= 1, "One QP budget: QPIterations{n}.");
    } else if constexpr (damp::is_same_v<Kind, kind::TrackedOutputs>) {
        static_assert(N <= 1, "One TrackedOutputs lists every output the target holds at r: TrackedOutputs{{i, ...}}.");
    } else if constexpr (damp::is_same_v<Kind, kind::MeasuredDisturbance>) {
        static_assert(N <= 1, "One MeasuredDisturbance carries every measured input: MeasuredDisturbance{.Bd = Bd}, a column of Bd per input.");
    } else {
        static_assert(N <= 1, "An option given twice: each kind of option is given at most once.");
    }
}

/**
 * @brief The kinds of option a design call accepts: check<Opts...>() refuses, at compile time, an
 *        option of no listed kind and a kind given twice (at_most_one).
 */
template<typename... Kinds>
struct option_set {
    /** @brief @p Opt is of a listed kind. */
    template<typename Opt>
    static constexpr bool listed = (design::OptionOf<Opt, Kinds> || ...);

    /** @brief Every option is of a listed kind. */
    template<typename... Opts>
    static constexpr bool known = (listed<Opts> && ...);

    /** @brief No kind is given twice. */
    template<typename... Opts>
    static constexpr bool unique = ((option_count<Kinds, Opts...> <= 1) && ...);

    /** @brief The compile-time check of a call's options. */
    template<typename... Opts>
    static constexpr void check() {
        static_assert(known<Opts...>, "An option of a type this design does not take: the trailing options are recognized by type.");
        (at_most_one<Kinds, option_count<Kinds, Opts...>>(), ...);
    }
};

} // namespace damp::detail
