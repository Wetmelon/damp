// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file input_limits.hpp
 * @brief InputLimits and InputRateLimits, the limits on a law's command and on its move, trailing
 *        options of design::state_feedback and design::mpc; and StateFeedbackLimits, the limits an
 *        LQ law takes at run time
 *
 * Both design families read them the same way, input_limits_of<NU, T>(opts...) and
 * input_rate_limits_of<NU, T>(opts...): the box among the options on the plant's inputs in its
 * scalar, unbounded without one. What each does with it is the design's (InputLimits,
 * InputRateLimits).
 */

#include "damp/backend.hpp"
#include "damp/design/options.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/core.hpp"
#include "damp/toolbox/bounds.hpp"

namespace damp::design {

/**
 * @struct InputLimits
 * @brief The limits on a law's command u, per input: among design::state_feedback's and
 *        design::mpc's options
 *
 * Not part of the cost: the weights are the same with or without them, so they are an option,
 * and left out the inputs are unbounded. An end left out is unbounded too, as a Bounds' is. The
 * ends are vectors, written flat: `InputLimits{.lower = {-12.0}, .upper = {12.0}}` for one input,
 * `{.lower = {0.0, -1.0}, .upper = {40.0, 1.0}}` for two. Named, it takes its sizes from a field
 * whose type gives them (`.lower = ColVec<2>{…}`, or `InputLimits<2>{…}`), one input in double
 * otherwise: braces deduce nothing, so a braced end is one input's.
 *
 * | Design | What the limits are |
 * | --- | --- |
 * | design::state_feedback (LQR, LQI, LQG, LQGI, OFLQG) | the command control returns is clamped into them, per input, and the law's memory holds against that clamp |
 * | design::mpc (MPC, OFMPC) | hard rows on u over the control horizon |
 *
 * The runtime's set_limits moves them (StateFeedbackLimits::input, MPCLimits::input), and
 * back_calculate stays for a stage that clamps after the law. Either design refuses an end that
 * is NaN or ±Inf, or a lower end above its upper one, with DesignStatus::invalid_parameter.
 *
 * The other limits are options of their own: InputRateLimits on the move, a DemandLimit on what
 * an LQ gain's outer terms ask of one state, an IntegralAction's on its integrators, and an MPC's
 * OutputLimits and StateLimits.
 *
 * @note Mirrors the MATLAB® mpc object's ManipulatedVariables Min/Max.
 */
template<size_t NU = 1, typename T = double>
struct InputLimits {
    static constexpr size_t inputs = NU;                                               ///< The inputs it bounds
    using option_kind = kind::InputLimits;                                             ///< Its kind of design option
    ColVec<NU, T> lower = Matrix<NU, 1, T>::constant(-damp::numeric_limits<T>::max()); ///< The lowest command per input (default unbounded)
    ColVec<NU, T> upper = Matrix<NU, 1, T>::constant(damp::numeric_limits<T>::max());  ///< The highest command per input (default unbounded)
};

/**
 * @struct InputRateLimits
 * @brief The limits on a law's move Δu = u_k − u_{k−1} per period, per input: among
 *        design::state_feedback's and design::mpc's options
 *
 * In the input's units per period, as MATLAB®'s RateMin/RateMax are (a SlewLimiter's rates are
 * per second): `InputRateLimits{.lower = {-0.2}, .upper = {0.2}}` moves one input at most 0.2 a
 * tick either way. It deduces as InputLimits does, and an end left out is unbounded.
 *
 * | Design | What the limits are |
 * | --- | --- |
 * | design::state_feedback (LQR, LQI, LQG, LQGI, OFLQG) | a clamp on the command's move: each tick the command is clamped into u_prev + [lower, upper] and the InputLimits, u_prev the law's share of the input applied last tick (back_calculate's u_sat − u_ff, else the last command); the law's memory holds against that clamp as against the InputLimits' |
 * | design::mpc (MPC, OFMPC) | hard rows on Δu over the control horizon, so the plan is one a slew-limited actuator can make, and a clamp on the applied move |
 *
 * Only a design with one pays for it: an MPC's QP carries the 2·NC·NU rows on Δu, and an LQ law
 * the clamp and the record of u_prev, only with it. Left out, an MPC shapes its moves with
 * MPCWeights::Rdu alone. One left unbounded keeps the rows and the clamp, for set_limits to move
 * at run time (StateFeedbackLimits::input_rate, MPCLimits::input_rate). Either design refuses an
 * end that is NaN or ±Inf, or a lower end above its upper one, with
 * DesignStatus::invalid_parameter.
 *
 * Where u_prev sits outside the InputLimits (after a limit change), the two boxes do not meet and
 * the InputLimits win: the command is their end nearest the move's box, as the MPC's clamp gives. After construction and reset() u_prev is 0, as an MPC's held input is, until
 * back_calculate records another (a bumpless takeover).
 *
 * @note Mirrors the MATLAB® mpc object's ManipulatedVariables RateMin/RateMax, hard as they are by
 *       default.
 */
template<size_t NU = 1, typename T = double>
struct InputRateLimits {
    static constexpr size_t inputs = NU;                                               ///< The inputs whose move it bounds
    using option_kind = kind::InputRateLimits;                                         ///< Its kind of design option
    ColVec<NU, T> lower = Matrix<NU, 1, T>::constant(-damp::numeric_limits<T>::max()); ///< The lowest move per period, per input (default unbounded)
    ColVec<NU, T> upper = Matrix<NU, 1, T>::constant(damp::numeric_limits<T>::max());  ///< The highest move per period, per input (default unbounded)
};

/**
 * @brief The limits an LQ law takes at run time: u, its move, the integrators ξ, what the outer
 *        terms ask of a state and the disturbance estimate d̂, in physical units. A field left out
 *        is unbounded (the demand limit off), and one the law does not have is empty.
 *
 * The gain, the filter and the target stay as designed. `auto l = law.limits(); l.input.upper(0)
 * = …; law.set_limits(l);` changes one, as on an MPC (MPCLimits, whose fields start and end the
 * same way); `l.demand = DemandLimit{…}` sets the demand limit, and `l.demand = {}` turns it off.
 * A design's result holds the limits its law starts with in the same struct
 * (StateFeedbackResult::limits).
 *
 * @tparam NX  Plant states (the demand limit's terms: these, then the integrators)
 * @tparam NI  Integrators (0: none, integral empty)
 * @tparam ND  Disturbance states of the filter (0: none, disturbance empty)
 * @tparam NDU Inputs with a move limit: NU with an InputRateLimits among the design's options (0: none, input_rate empty)
 */
template<size_t NX, size_t NU, typename T = double, size_t NI = 0, size_t ND = 0, size_t NDU = 0>
struct StateFeedbackLimits {
    Bounds<NU, T>                           input{};       ///< Limits on the command (InputLimits)
    Bounds<NDU, T>                          input_rate{};  ///< Limits on the command's move per period (InputRateLimits; NDU > 0)
    Bounds<NI, T>                           integral{};    ///< Limits on the integrators ξ (IntegralAction::limits; NI > 0)
    ::damp::detail::DemandTerms<NX + NI, T> demand{};      ///< A limit on what the outer terms ask of one state (a DemandLimit; default off)
    Bounds<ND, T>                           disturbance{}; ///< Limits on the disturbance estimate d̂, in d's units (ND > 0)

    /** @brief The limits in another scalar (an unbounded end stays one). */
    template<typename U>
    [[nodiscard]] constexpr StateFeedbackLimits<NX, NU, U, NI, ND, NDU> as() const {
        return {input.template as<U>(), input_rate.template as<U>(), integral.template as<U>(), demand.template as<U>(), disturbance.template as<U>()};
    }
};

namespace detail {

/** @brief The box of a limit option's @p lower and @p upper ends, in the scalar @p T: Bounds::as, so an unbounded end stays one. */
template<typename T, size_t N, typename U>
[[nodiscard]] constexpr Bounds<N, T> box_of(const ColVec<N, U>& lower, const ColVec<N, U>& upper) {
    return Bounds<N, U>{.lower = lower, .upper = upper}.template as<T>();
}

/** @brief The box of the option of the kind @p Kind among @p opts, on the plant's NU inputs in its scalar @p T (none: unbounded). */
template<typename Kind, size_t NU, typename T, typename Fallback, typename... Opts>
[[nodiscard]] constexpr Bounds<NU, T> input_box_of(const Opts&... opts) {
    using Given = option_t<Kind, Fallback, Opts...>;
    if constexpr (Given::inputs != NU) {
        static_assert(
            Given::inputs == NU,
            "InputLimits and InputRateLimits bound each input of the plant: InputLimits<NU>{...}, or .lower = ColVec<NU>{...}. "
            "Given in a template? There clang deduces nothing from a field set to a constant declared outside it (llvm/llvm-project#114870): name the size."
        );
        return Bounds<NU, T>{};
    } else {
        const auto o = option<Kind>(Fallback{}, opts...);
        return box_of<T>(o.lower, o.upper);
    }
}

/** @brief The InputLimits among @p opts on the plant's NU inputs, in its scalar @p T (none: unbounded). */
template<size_t NU, typename T, typename... Opts>
[[nodiscard]] constexpr Bounds<NU, T> input_limits_of(const Opts&... opts) {
    return input_box_of<kind::InputLimits, NU, T, InputLimits<NU, T>>(opts...);
}

/** @brief Inputs with a move limit among @p opts: NU with an InputRateLimits, else 0. */
template<size_t NU, typename... Opts>
inline constexpr size_t rate_limited_inputs = (option_count<kind::InputRateLimits, Opts...> > 0) ? NU : 0;

/** @brief The InputRateLimits among @p opts on the plant's NU inputs, in its scalar @p T (none: unbounded). */
template<size_t NU, typename T, typename... Opts>
[[nodiscard]] constexpr Bounds<NU, T> input_rate_limits_of(const Opts&... opts) {
    return input_box_of<kind::InputRateLimits, NU, T, InputRateLimits<NU, T>>(opts...);
}

} // namespace detail

} // namespace damp::design
