// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file pose.hpp
 * @brief Rigid-body pose — the shared interchange type for the kinematics solvers.
 *
 * A Pose is a translation-3-vector plus a unit quaternion, *not* a 4×4
 * homogeneous matrix: composing a kinematic chain as `(q, t)` pairs is roughly
 * half the FLOPs of stacking `Mat4×Mat4` per joint (≈31 vs 64 multiplies), carries
 * 7 scalars instead of 16, and the quaternion renormalizes cheaply where a 4×4
 * rotation block drifts from orthonormal along a long chain. Transform4
 * (in geometry.hpp) is retained only as an interop/export convenience — convert
 * with Pose::to_transform4 / Pose::from_transform4.
 *
 * Composition convention: `a * b` expresses pose `b` *in the frame of* `a` — i.e.
 * if `a` is the pose of frame B in world and `b` the pose of frame C in B, then
 * `a * b` is the pose of C in world. `transform_point(p)` maps a point from the
 * local frame into the parent.
 *
 * @see damp/math/geometry.hpp for Quaternion / DCM / Euler / Transform4.
 */

#include "damp/math/geometry.hpp" // Quaternion, Vec3, Transform4

namespace damp {

/**
 * @brief A 3-D translation — a thin @ref Vec3 with domain-named conveniences.
 * @tparam T Scalar type (floating point)
 */
template<typename T = double>
struct Translation3 : Vec3<T> {
    using Base = Vec3<T>;
    using Base::Base; // inherit the ColVec constructors

    constexpr Translation3() = default;
    constexpr Translation3(T x, T y, T z) {
        (*this)[0] = x;
        (*this)[1] = y;
        (*this)[2] = z;
    }
    constexpr explicit Translation3(const Base& v) : Base(v) {}

    /// Euclidean distance to another translation.
    [[nodiscard]] constexpr T distance(const Translation3& other) const {
        return (static_cast<const Base&>(*this) - static_cast<const Base&>(other)).norm();
    }
};

/**
 * @brief Rigid-body pose: a translation and an orientation (unit quaternion).
 * @tparam T Scalar type (floating point)
 */
template<typename T = double>
struct Pose {
    Quaternion<T>   orientation{Quaternion<T>::identity()};
    Translation3<T> translation{};

    /// The identity pose (origin, no rotation).
    [[nodiscard]] static constexpr Pose identity() { return Pose{}; }

    /// Compose: `*this * rhs` places `rhs` in this pose's frame (see file note).
    [[nodiscard]] constexpr Pose operator*(const Pose& rhs) const {
        const Vec3<T> rotated = orientation.rotate(static_cast<const Vec3<T>&>(rhs.translation));
        return Pose{
            .orientation = orientation * rhs.orientation,
            .translation = Translation3<T>(static_cast<const Vec3<T>&>(translation) + rotated),
        };
    }

    /// Map a point from this pose's local frame into the parent frame.
    [[nodiscard]] constexpr Vec3<T> transform_point(const Vec3<T>& p) const {
        return static_cast<const Vec3<T>&>(translation) + orientation.rotate(p);
    }

    /// The inverse pose (assumes a unit quaternion).
    [[nodiscard]] constexpr Pose inverse() const {
        const Quaternion<T> q_inv = orientation.conjugate();
        return Pose{
            .orientation = q_inv,
            .translation = Translation3<T>(-q_inv.rotate(static_cast<const Vec3<T>&>(translation))),
        };
    }

    /// Export to a 4×4 homogeneous transform (interop convenience).
    [[nodiscard]] constexpr Transform4<T> to_transform4() const {
        return Transform4<T>::from_quaternion_translation(orientation, static_cast<const Vec3<T>&>(translation));
    }

    /// Build from a 4×4 homogeneous transform.
    [[nodiscard]] static constexpr Pose from_transform4(const Transform4<T>& tf) {
        const auto& [q, t] = tf.to_quaternion_translation(); // {Quaternion, Vec3}
        return {q, t};
    }

    /// Rebind to another scalar type.
    /// Designated initializers follow declaration order: orientation, then translation.
    template<typename U>
    [[nodiscard]] constexpr Pose<U> as() const {
        return {
            .orientation = Quaternion<U>{
                static_cast<U>(orientation.w()),
                static_cast<U>(orientation.x()),
                static_cast<U>(orientation.y()),
                static_cast<U>(orientation.z()),
            },
            .translation = Translation3<U>{
                static_cast<U>(translation[0]),
                static_cast<U>(translation[1]),
                static_cast<U>(translation[2]),
            },
        };
    }
};

} // namespace damp
