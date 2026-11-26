// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file poles.hpp
 * @brief Poles, the closed-loop poles a placed design takes in its weights' or rule's slot:
 *        design::state_feedback's gain, and design::pid's rule on a first-order plant
 */

#include "damp/backend.hpp"
#include "damp/design/options.hpp"
#include "damp/math/complex.hpp"

namespace damp::design {

namespace detail {

/** @brief The scalar of a pole as written: a floating-point number is its own. */
template<typename P>
struct pole_scalar {
    using type = P; ///< The number's type
};

/** @brief A damp::complex: its parts'. */
template<typename U>
struct pole_scalar<damp::complex<U>> {
    using type = U; ///< The parts' type
};

/** @brief A list of poles in the scalar @p T, as the placement takes them. */
template<typename T, size_t N, typename U>
[[nodiscard]] constexpr damp::array<damp::complex<T>, N> pole_list(const damp::complex<U> (&poles)[N]) {
    damp::array<damp::complex<T>, N> out{};
    for (size_t i = 0; i < N; ++i) {
        out[i] = damp::complex<T>(poles[i]);
    }
    return out;
}

} // namespace detail

/**
 * @struct Poles
 * @brief The closed-loop poles of a placed design: design::state_feedback's weights' slot, in
 *        place of LQRWeights, and design::pid's rule on a first-order plant
 *
 * `design::state_feedback(plant, Ts, Poles{-8.0, -12.0}, opts...)` places eig(A − B·K) where the
 * list says, by design::place's machinery (KNV, Ackermann for one input), and builds the same law
 * the LQ weights would, with every option: the result and its runtime are the LQ design's types.
 * The count is the plant's states plus the integrators of an IntegralAction, refused at compile
 * time otherwise; an InputDelay's delay states, which only hold inputs already sent, take theirs
 * at 0 (deadbeat), and a disturbance spec's states none, since the target cancels them.
 *
 * `design::pid(plant, Ts, Poles{-50.0, -50.0})` places the poles of a PI (two) or a PID (three) on
 * a first-order plant, exact at Ts.
 *
 * The list is in the plant's own time domain, as design::place takes it: rad/s (s-plane) on the
 * sampled form, mapped to z = e^{s·Ts}; z-plane on a discrete plant. Complex poles come as
 * adjacent conjugate pairs, `Poles{damp::complex{-1.0, 2.0}, damp::complex{-1.0, -2.0}}`, and a
 * real one converts, `Poles{-5.0, -5.0}`. The count and the scalar are deduced from the entries,
 * with bare braces (a deduction guide over the C array): the first entry's scalar, real or
 * complex, is the list's.
 *
 * @tparam N The poles listed
 * @tparam T Their scalar (the plant's, or converted to it)
 */
template<size_t N, typename T = double>
struct Poles {
    using option_kind = kind::Weights;  ///< Its kind of design option
    static constexpr size_t listed = N; ///< The poles it lists
    static_assert(damp::is_floating_point_v<T>, "A pole is a floating-point number or a damp::complex of one: Poles{-5.0, -5.0}, Poles{damp::complex{-1.0, 2.0}, damp::complex{-1.0, -2.0}}.");

    damp::complex<T> poles[N]{}; ///< The closed-loop poles, adjacent conjugate pairs
};

/** @brief Poles{p, …}: the count from the list, the scalar from its first entry. */
template<typename P, typename... Rest>
Poles(P, Rest...) -> Poles<1 + sizeof...(Rest), typename detail::pole_scalar<P>::type>;

} // namespace damp::design
