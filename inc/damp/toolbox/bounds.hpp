// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file bounds.hpp
 * @brief Per-channel interval limits — the one primitive for box constraints on
 *        a state, output, or input vector.
 *
 * Every "limit" in the library is the same shape: a closed interval
 * @f$[\,\text{lower}_i,\ \text{upper}_i\,]@f$ per channel. A state-limit set is a
 * `Bounds\<Nx\>`; an output/input-limit set is a `Bounds\<Nu\>` — same type, different
 * vector. The ends are column vectors, so they enter matrix expressions directly.
 * A scalar saturation pair is a `Bounds<1>`, written `{lo, hi}`: the PID output and
 * integrator limits, the ESC and relay-autotune clamps, the SOGI-FLL frequency range.
 * The vector form carries the MPC input, move and output limits and the MHE state
 * limits. The symmetric @f$|x_i|\le \text{max}_i@f$ convention is Bounds::symmetric.
 *
 * Scope is just the hard constraint. A soft/weighted constraint (MPC: a `Q`/`R`
 * cost coefficient or a slack penalty per channel) composes a `Bounds` with a
 * weight; it is not a member here.
 *
 * @see controllers/pid.hpp (Bounds<1> on output + integrator),
 *      controllers/mpc.hpp (MPCConstraints), estimation/mhe.hpp (MHEConstraints).
 */


#include "damp/backend.hpp" // clamp
#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/core.hpp"

namespace damp {

/**
 * @brief A per-channel closed-interval box constraint.
 *
 * An aggregate: `{.lower = -12.0, .upper = 12.0}` for one channel (a ColVec<1> takes its
 * scalar), `{.lower = {-1.0, -2.0}, .upper = {1.0, 2.0}}` for two. An end left out stays
 * unbounded. The default is unbounded: every end at the widest finite value, ∓max(). The
 * library builds under -ffinite-math-only, where producing ±inf is undefined, so
 * the ends are never ±inf; a zero-initialised box would clamp every channel to 0.
 *
 * @tparam N Number of channels
 * @tparam T Scalar type (floating point)
 */
template<size_t N, typename T = double>
struct Bounds {
    ColVec<N, T> lower = Matrix<N, 1, T>::constant(-damp::numeric_limits<T>::max()); ///< lower(i) ≤ x(i)
    ColVec<N, T> upper = Matrix<N, 1, T>::constant(damp::numeric_limits<T>::max());  ///< x(i) ≤ upper(i)

    /** @brief Symmetric box [−mag, +mag] per channel (the |xᵢ| ≤ magᵢ convention). */
    [[nodiscard]] static constexpr Bounds symmetric(const ColVec<N, T>& mag) {
        return Bounds{
            .lower = -mag,
            .upper = mag,
        };
    }

    /** @brief Per-channel clamp into the box. */
    [[nodiscard]] constexpr ColVec<N, T> saturate(const ColVec<N, T>& x) const {
        ColVec<N, T> out{};
        for (size_t i = 0; i < N; ++i) {
            out(i) = damp::clamp(x(i), lower(i), upper(i));
        }
        return out;
    }

    /** @brief True iff every channel is within its interval. */
    [[nodiscard]] constexpr bool contains(const ColVec<N, T>& x) const {
        for (size_t i = 0; i < N; ++i) {
            if (x(i) < lower(i) || x(i) > upper(i)) {
                return false;
            }
        }
        return true;
    }

    /** @brief True iff the box is usable: no end is NaN (∓max() for none) and no lower end is above its upper one. */
    [[nodiscard]] constexpr bool valid() const {
        for (size_t i = 0; i < N; ++i) {
            if (!damp::isfinite(lower(i)) || !damp::isfinite(upper(i)) || !(lower(i) <= upper(i))) {
                return false;
            }
        }
        return true;
    }

    /** @brief Clamp of a scalar into the SISO interval. */
    [[nodiscard]] constexpr T saturate(T x) const
        requires(N == 1)
    {
        return damp::clamp(x, lower(0), upper(0));
    }

    /** @brief True iff the scalar is within the SISO interval. */
    [[nodiscard]] constexpr bool contains(T x) const
        requires(N == 1)
    {
        return (x >= lower(0)) && (x <= upper(0));
    }

    /**
     * @brief The same box in scalar type U. An end at or beyond max()/4 in T is unbounded and
     * becomes ±max() of U; a plain cast of numeric_limits<double>::max() to float would overflow.
     */
    template<typename U>
    [[nodiscard]] constexpr Bounds<N, U> as() const {
        Bounds<N, U> out{};
        for (size_t i = 0; i < N; ++i) {
            out.lower(i) = convert_end<U>(lower(i));
            out.upper(i) = convert_end<U>(upper(i));
        }
        return out;
    }

private:
    template<typename U>
    [[nodiscard]] static constexpr U convert_end(T v) {
        constexpr T unbounded = damp::numeric_limits<T>::max() / T{4};
        if (v >= unbounded) {
            return damp::numeric_limits<U>::max();
        }
        if (v <= -unbounded) {
            return -damp::numeric_limits<U>::max();
        }
        return static_cast<U>(v);
    }
};

} // namespace damp
