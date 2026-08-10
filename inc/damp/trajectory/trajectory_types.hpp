// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file trajectory_types.hpp
 * @brief Shared value types and low-level helpers for the trajectory generators.
 *
 * Defines the data structures that appear across multiple profile families:
 * - TrajectoryLimits — kinematic bounds (Vmax, Amax, Dmax, Jmax)
 * - TrajectoryState — time-domain command (position…jerk) for servo / PTP / S-curve
 * - Jet — order-N derivative chain (scales past jerk for high-order polys / cams)
 * - TrajectoryBoundary — fixed-size endpoint BCs (p…jerk) for common poly orders
 * - design::detail::factorial / design::detail::falling_factorial — shared
 *   polynomial-coefficient helpers used by both @ref polynomial.hpp and @ref spline.hpp
 *
 * Include this header only when you need the shared types without pulling in a full
 * profile family. Normally just include the profile header you need (e.g.
 * `damp/trajectory/scurve.hpp`) — each one already pulls this in.
 */

#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"

namespace damp {

/**
 * @brief Asymmetric kinematic limits for a trapezoidal or S-curve motion profile.
 *
 * @p max_acceleration bounds the approach-to-cruise ramp and @p max_deceleration
 * the approach-to-target ramp (they may differ). All three must be > 0.
 *
 * Units are consistent with the caller's length unit LU: velocity LU/s,
 * acceleration LU/s², jerk LU/s³. Profile times (`Ta`, `Tf`, `duration`, and
 * runtime `step(dt)`) are in seconds. Use the same LU for position endpoints.
 *
 * @tparam T Scalar type
 */
template<typename T = double>
struct TrajectoryLimits {

    T max_velocity{T{0}};     ///< |v| ≤ max_velocity [LU/s] (> 0)
    T max_acceleration{T{0}}; ///< accel-region bound Amax [LU/s²] (> 0)
    T max_deceleration{T{0}}; ///< decel-region bound Dmax [LU/s²] (> 0)
    T max_jerk{T{0}};         ///< |j| ≤ max_jerk [LU/s³] (> 0; S-curve only, ignored by trapezoidal)

    /// Trapezoidal validity: positive v/a/d limits (jerk is not required).
    [[nodiscard]] constexpr bool valid() const {
        return (max_velocity > T{0}) && (max_acceleration > T{0}) && (max_deceleration > T{0});
    }

    /// S-curve validity: additionally requires a positive jerk limit.
    [[nodiscard]] constexpr bool valid_jerk_limited() const { return valid() && (max_jerk > T{0}); }
};

/**
 * @brief A point on a motion profile: commanded position, velocity, acceleration.
 *
 * Time-domain servo command (feedforward slot). Caps at jerk — enough for
 * trapezoid / S-curve / cascade FF. Higher-order polys (nonic, …) use Jet
 * so the type does not grow a snap/crackle field list.
 *
 * Units match TrajectoryLimits: position [LU], velocity [LU/s], acceleration
 * [LU/s²], jerk [LU/s³]. Time arguments that produce a TrajectoryState are
 * in seconds.
 *
 * @tparam T Scalar type
 * @see Jet for order-N derivative chains (design / high-order eval)
 */
template<typename T = double>
struct TrajectoryState {

    T position{T{0}};     ///< Commanded position [LU]
    T velocity{T{0}};     ///< Commanded velocity [LU/s]
    T acceleration{T{0}}; ///< Commanded acceleration [LU/s²]
    T jerk{T{0}};         ///< Commanded jerk [LU/s³] (meaningful for poly / S-curve)
};

/**
 * @brief Truncated jet of a scalar motion sample: @f$ d[k] = s^{(k)} @f$.
 *
 * @f$ d[0] @f$ is position, @f$ d[1] @f$ velocity, @f$ d[2] @f$ acceleration, …
 * Scales with @p N so septic/nonic (and cam parameter maps) can carry snap and
 * beyond without extending TrajectoryState.
 *
 * @tparam N Number of stored derivatives (including position as order 0)
 * @tparam T Scalar type
 */
template<size_t N, typename T = double>
struct Jet {
    static_assert(N >= 1, "Jet needs at least position (N >= 1)");

    static constexpr size_t size = N;

    damp::array<T, N> d{}; ///< d[k] = k-th derivative of the scalar sample

    [[nodiscard]] constexpr T&       operator[](size_t i) { return d[i]; }
    [[nodiscard]] constexpr const T& operator[](size_t i) const { return d[i]; }

    [[nodiscard]] constexpr T position() const { return d[0]; }
    [[nodiscard]] constexpr T velocity() const {
        if constexpr (N > 1) {
            return d[1];
        }
        return T{0};
    }

    [[nodiscard]] constexpr T acceleration() const {
        if constexpr (N > 2) {
            return d[2];
        }
        return T{0};
    }

    [[nodiscard]] constexpr T jerk() const {
        if constexpr (N > 3) {
            return d[3];
        }
        return T{0};
    }

    [[nodiscard]] constexpr T snap() const {
        if constexpr (N > 4) {
            return d[4];
        }
        return T{0};
    }

    /// Truncate or zero-pad into the time-domain command type (drops order ≥ 4).
    [[nodiscard]] constexpr TrajectoryState<T> as_trajectory_state() const {
        if constexpr (N > 3) {
            return TrajectoryState<T>{
                .position = d[0],
                .velocity = d[1],
                .acceleration = d[2],
                .jerk = d[3],
            };
        } else if constexpr (N > 2) {
            return TrajectoryState<T>{
                .position = d[0],
                .velocity = d[1],
                .acceleration = d[2],
            };
        } else if constexpr (N > 1) {
            return TrajectoryState<T>{.position = d[0], .velocity = d[1]};
        } else {
            return TrajectoryState<T>{.position = d[0]};
        }
    }

    template<typename U>
    [[nodiscard]] constexpr Jet<N, std::remove_const_t<U>> as() const {
        using O = std::remove_const_t<U>;
        Jet<N, O> out{};
        for (size_t i = 0; i < N; ++i) {
            out.d[i] = static_cast<O>(d[i]);
        }
        return out;
    }
};

/**
 * @brief Boundary conditions at one endpoint of a polynomial trajectory: a
 *        position and its time derivatives through jerk.
 *
 * Convenience for cubic…septic (design::poly_trajectory with Order ≤ 7).
 * For nonic and above, pass a Jet of length @f$ (Order+1)/2 @f$ instead.
 *
 * @tparam T Scalar type
 * @see Jet
 */
template<typename T = double>
struct TrajectoryBoundary {
    T position{T{0}};
    T velocity{T{0}};
    T acceleration{T{0}};
    T jerk{T{0}};

    /// Pack into a jet of length @p N (extra slots stay zero; N > 4 ignores snap+).
    template<size_t N>
    [[nodiscard]] constexpr Jet<N, T> as_jet() const {
        Jet<N, T> j{};
        j.d[0] = position;
        if constexpr (N > 1) {
            j.d[1] = velocity;
        }
        if constexpr (N > 2) {
            j.d[2] = acceleration;
        }
        if constexpr (N > 3) {
            j.d[3] = jerk;
        }
        return j;
    }
};

namespace design::detail {

/// k! (exact for modest k in float/double; used up through nonic BVP / jets).
template<typename T>
constexpr T factorial(size_t k) {
    T f{1};
    for (size_t i = 2; i <= k; ++i) {
        f *= static_cast<T>(i);
    }
    return f;
}

/// Falling factorial i·(i−1)···(i−k+1) = i! / (i−k)! — the k-th derivative
/// coefficient of tⁱ. Zero when i < k.
template<typename T>
constexpr T falling_factorial(size_t i, size_t k) {
    if (i < k) {
        return T{0};
    }
    T f{1};
    for (size_t m = 0; m < k; ++m) {
        f *= static_cast<T>(i - m);
    }
    return f;
}

} // namespace design::detail

} // namespace damp
