// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file input_limits.hpp
 * @brief InputLimits: the limits on a law's command and on its moves, a trailing option of
 *        design::state_feedback and design::mpc
 *
 * Both design families read it the same way, input_limits_of<NU, T>(opts...): the one among the
 * options on the plant's inputs in its scalar, unbounded without one. What each does with it is
 * the design's (InputLimits).
 */

#include "damp/backend.hpp"
#include "damp/design/options.hpp"
#include "damp/math/math.hpp"
#include "damp/toolbox/bounds.hpp"

namespace damp::design {

/**
 * @struct InputLimits
 * @brief The limits on a law's command u and on its move Δu per period: among
 *        design::state_feedback's and design::mpc's options
 *
 * Not part of the cost: the weights are the same with or without them, so they are an option,
 * and left out the inputs are unbounded. Named, it takes its sizes from its fields, one
 * input unless a field says otherwise (`.u_lim = Bounds<2>{…}`, or `InputLimits<2>{…}`):
 * `InputLimits{.u_lim = {.lower = -12.0, .upper = 12.0}}`, or with `.du_lim = {.lower = -0.5,
 * .upper = 0.5}` too, is InputLimits<1>.
 *
 * | Field | design::state_feedback (LQR, LQI, LQG, LQGI, OFLQG) | design::mpc (MPC, OFMPC) |
 * | --- | --- | --- |
 * | u_lim | the command control returns is clamped into it, per input, and the law's memory holds against that clamp | hard rows on u over the control horizon |
 * | du_lim | refused unless unbounded (DesignStatus::invalid_parameter) | hard rows on Δu over the control horizon |
 *
 * An LQ law clamps its command; a move limit is the MPC's, whose QP plans the moves within it.
 * u_lim is a member of an LQ runtime, `law.u_lim`, so a rail that moves is set at run time, and
 * back_calculate stays for a stage that clamps after the law. An MPC's runtime has both as
 * members, `mpc.u_lim` and `mpc.du_lim`, which its next tick takes into the QP's rows. Either
 * design refuses a limit that is NaN, or a lower limit above its upper one, with
 * DesignStatus::invalid_parameter.
 *
 * The other limits are options of their own: a DemandLimit on what an LQ gain's outer terms ask
 * of one state, an IntegralAction's on its integrators, and an MPC's OutputLimits and StateLimits.
 *
 * @note Mirrors the MATLAB® mpc object's ManipulatedVariables Min/Max/RateMin/RateMax.
 */
template<size_t NU = 1, typename T = double>
struct InputLimits {
    Bounds<NU, T> u_lim{};  ///< Limits on the command u, per input (default unbounded)
    Bounds<NU, T> du_lim{}; ///< Limits on the move Δu per period, per input (default unbounded; an MPC's alone)

    /** @brief The limits in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr InputLimits<NU, U> as() const {
        return InputLimits<NU, U>{u_lim.template as<U>(), du_lim.template as<U>()};
    }
};

namespace detail {

/** @brief An InputLimits: a kind of design::state_feedback and design::mpc option. */
template<typename O>
struct is_input_limits : damp::false_type {};

/** @brief InputLimits of any size. */
template<size_t NU, typename T>
struct is_input_limits<InputLimits<NU, T>> : damp::true_type {};

/** @brief The inputs an InputLimits is on. */
template<typename L>
inline constexpr size_t input_limits_inputs = 0;

/** @brief InputLimits<NU, …>: NU. */
template<size_t NU, typename T>
inline constexpr size_t input_limits_inputs<InputLimits<NU, T>> = NU;

/**
 * @brief No limit at all: every end of @p b at or beyond ∓max()/4, the sentinels Bounds::as keeps
 *        unbounded. A NaN end is a limit, as Bounds::valid reads it: tested by its bits
 *        (damp::isfinite), since -ffast-math may fold a comparison with a NaN.
 */
template<size_t N, typename T>
[[nodiscard]] constexpr bool unlimited(const Bounds<N, T>& b) {
    constexpr T beyond = damp::numeric_limits<T>::max() / T{4};
    for (size_t i = 0; i < N; ++i) {
        if (!damp::isfinite(b.lower(i)) || !damp::isfinite(b.upper(i)) || b.lower(i) > -beyond || b.upper(i) < beyond) {
            return false;
        }
    }
    return true;
}

/** @brief The InputLimits among @p opts on the plant's NU inputs, in its scalar @p T (none: unbounded). */
template<size_t NU, typename T, typename... Opts>
[[nodiscard]] constexpr InputLimits<NU, T> input_limits_of(const Opts&... opts) {
    using Given = ::damp::detail::option_t<is_input_limits, InputLimits<NU, T>, Opts...>;
    if constexpr (input_limits_inputs<Given> != NU) {
        static_assert(input_limits_inputs<Given> == NU, "InputLimits bound each input of the plant: InputLimits<NU>{...}, or .u_lim = Bounds<NU>{...}.");
        return InputLimits<NU, T>{};
    } else {
        return ::damp::detail::option<is_input_limits>(InputLimits<NU, T>{}, opts...).template as<T>();
    }
}

} // namespace detail

} // namespace damp::design
