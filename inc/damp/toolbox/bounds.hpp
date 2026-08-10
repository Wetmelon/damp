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
 * vector. This subsumes the scalar `u_min/u_max` saturation pairs (`Bounds<1>`)
 * and the symmetric @f$|x_i|\le \text{max}_i@f$ convention used by
 * `JointLimits` / `TrajectoryLimits` (see Bounds::symmetric).
 *
 * Scope is deliberately just the *hard* constraint. A soft/weighted constraint
 * (MPC: a `Q`/`R` cost coefficient or a slack penalty per channel) is a separate
 * concern that *composes* a `Bounds` with a weight — it is intentionally not a
 * member here, so the limit type stays usable in every controller today and the
 * MPC layer is purely additive later.
 *
 * @see controllers/pid.hpp (Bounds<1> on output + integrator),
 *      trajectory/cartesian_move.hpp (JointLimits), trajectory/trajectory_types.hpp.
 */

#include <cstddef>
#include <limits>

#include "damp/backend.hpp" // array, clamp

namespace damp {

/**
 * @brief A per-channel closed-interval box constraint.
 * @tparam N Number of channels
 * @tparam T Scalar type (floating point)
 */
template<size_t N, typename T = double>
struct Bounds {
    damp::array<T, N> lower; ///< lower[i] ≤ x[i]
    damp::array<T, N> upper; ///< x[i] ≤ upper[i]

    /// Default is *unbounded* — the widest finite interval. Uses min()/max()
    /// rather than ±inf because the library builds under -ffinite-math-only
    /// (-ffast-math), where producing ±inf is UB (see math/constexpr_math.hpp).
    /// A zero-initialised box would silently clamp every channel to 0.
    constexpr Bounds() {
        lower.fill(-std::numeric_limits<T>::max());
        upper.fill(std::numeric_limits<T>::max());
    }

    constexpr Bounds(const damp::array<T, N>& lo, const damp::array<T, N>& hi) : lower(lo), upper(hi) {}

    /// Scalar construction for the SISO case, e.g. a PID output limit.
    constexpr Bounds(T lo, T hi)
        requires(N == 1)
        : lower{lo}, upper{hi} {}

    /// Symmetric box [-mag, +mag] per channel (the |xᵢ|≤maxᵢ convention).
    [[nodiscard]] static constexpr Bounds symmetric(const damp::array<T, N>& mag) {
        Bounds b{};
        for (size_t i = 0; i < N; ++i) {
            b.lower[i] = -mag[i];
            b.upper[i] = mag[i];
        }
        return b;
    }

    /// Per-channel clamp into the box.
    [[nodiscard]] constexpr damp::array<T, N> saturate(const damp::array<T, N>& x) const {
        damp::array<T, N> out{};
        for (size_t i = 0; i < N; ++i) {
            out[i] = damp::clamp(x[i], lower[i], upper[i]);
        }
        return out;
    }

    /// True iff every channel is within its interval.
    [[nodiscard]] constexpr bool contains(const damp::array<T, N>& x) const {
        for (size_t i = 0; i < N; ++i) {
            if (x[i] < lower[i] || x[i] > upper[i]) {
                return false;
            }
        }
        return true;
    }

    // ---- SISO scalar conveniences (N == 1) ---------------------------------
    [[nodiscard]] constexpr T saturate(T x) const
        requires(N == 1)
    {
        return damp::clamp(x, lower[0], upper[0]);
    }
    [[nodiscard]] constexpr bool contains(T x) const
        requires(N == 1)
    {
        return (x >= lower[0]) && (x <= upper[0]);
    }
};

} // namespace damp
