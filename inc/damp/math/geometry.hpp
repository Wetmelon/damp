// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file geometry.hpp
 * @brief 3D rotations and homogeneous transforms: DCM, Quaternion, Euler, Transform4
 *
 * Working attitude types for estimation and robotics. Rigid-body pose
 * composition uses @c Quaternion + translation (@ref Pose in @c pose.hpp);
 * @ref Transform4 is retained for DH / homogeneous interop only.
 *
 * Convention: @c EulerZYX = aerospace yaw-pitch-roll; @c EulerXYZ = robotics
 * roll-pitch-yaw — do not conflate them.
 *
 * @see Solà et al., "Quaternion kinematics for the error-state Kalman filter" (2017)
 * @see Diebel, "Representing Attitude: Euler Angles, Unit Quaternions, and Rotation Vectors" (2006)
 */

#include <cstddef>
#include <initializer_list>
#include <span>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"

namespace damp {

/**
 * @brief Intrinsic Euler angle sequence (axis order of successive principal rotations)
 *
 * @c ZYX is aerospace yaw-pitch-roll; @c XYZ is robotics roll-pitch-yaw.
 */
enum class EulerOrder {
    XYZ, ///< Roll-pitch-yaw (robotics RPY)
    ZYX, ///< Yaw-pitch-roll (aerospace YPR)
    ZXY, ///< Z then X then Y
    YXZ, ///< Y then X then Z
    YZX, ///< Y then Z then X
    XZY, ///< X then Z then Y
};

// Forward declarations
template<typename T>
struct Quaternion;

template<typename T>
struct DCM;

template<typename T, EulerOrder Order>
struct Euler;

/**
 * @brief Direction cosine matrix — 3×3 rotation (SO(3) wrapper over @ref Mat3)
 * @tparam T Floating-point scalar
 */
template<typename T>
struct DCM : public Mat3<T> {

    static_assert(std::is_floating_point_v<T>, "DCM element type must be floating point");

    constexpr DCM() : Mat3<T>(Mat3<T>::identity()) {}
    constexpr DCM(const DCM&) = default;
    constexpr DCM& operator=(const DCM&) = default;
    constexpr DCM(DCM&&) = default;
    constexpr DCM& operator=(DCM&&) = default;
    constexpr ~DCM() = default;

    /// Construct from a raw @ref Mat3 (entries used as-is).
    constexpr explicit DCM(const Mat3<T>& m) : Mat3<T>(m) {}

    /// Row-nested initializer list (prvalue / RVO-friendly).
    constexpr DCM(std::initializer_list<std::initializer_list<T>> init) : Mat3<T>(init) {}

    /// Identity rotation.
    [[nodiscard]] static constexpr DCM identity() { return DCM{}; }

    /// Principal-axis rotation about X by @p angle [rad].
    [[nodiscard]] static constexpr DCM rotate_x(T angle) {
        const auto [s, c] = damp::sincos(angle);
        return DCM{
            {T{1}, T{0}, T{0}},
            {T{0}, c, -s},
            {T{0}, s, c},
        };
    }

    /// Principal-axis rotation about Y by @p angle [rad].
    [[nodiscard]] static constexpr DCM rotate_y(T angle) {
        const auto [s, c] = damp::sincos(angle);
        return DCM{
            {c, T{0}, s},
            {T{0}, T{1}, T{0}},
            {-s, T{0}, c},
        };
    }

    /// Principal-axis rotation about Z by @p angle [rad].
    [[nodiscard]] static constexpr DCM rotate_z(T angle) {
        const auto [s, c] = damp::sincos(angle);
        return DCM{
            {c, -s, T{0}},
            {s, c, T{0}},
            {T{0}, T{0}, T{1}},
        };
    }

    /// Compose rotations: R_this R_rhs.
    [[nodiscard]] constexpr DCM operator*(const DCM& rhs) const {
        return DCM(static_cast<const Mat3<T>&>(*this) * static_cast<const Mat3<T>&>(rhs));
    }

    /// In-place composition with @p rhs.
    constexpr DCM& operator*=(const DCM& rhs) {
        return *this = *this * rhs;
    }

    /// Rotate a 3-vector: R v.
    [[nodiscard]] constexpr Vec3<T> operator*(const Vec3<T>& v) const {
        return Vec3<T>(static_cast<const Mat3<T>&>(*this) * static_cast<const Matrix<3, 1, T>&>(v));
    }

    /// Matrix transpose (inverse when @c R is orthonormal).
    [[nodiscard]] constexpr DCM transpose() const {
        return DCM(Mat3<T>::transpose());
    }

    /// Inverse rotation (transpose for orthonormal @c R).
    [[nodiscard]] constexpr DCM inverse() const {
        return transpose();
    }

    /// Underlying @ref Mat3 view.
    [[nodiscard]] constexpr const Mat3<T>& matrix() const { return *this; }
    /// Underlying @ref Mat3 view (mutable).
    [[nodiscard]] constexpr Mat3<T>& matrix() { return *this; }

    /**
     * @brief Convert to a unit quaternion
     * @param eps Normalization guard forwarded to @ref Quaternion::from_dcm
     * @return Unit quaternion, or nullopt if the conversion is ill-conditioned
     */
    [[nodiscard]] constexpr damp::optional<Quaternion<T>> to_quaternion(T eps = static_cast<T>(1e-6)) const;

    /// Convert to Euler angles in sequence @p Order.
    template<EulerOrder Order>
    [[nodiscard]] constexpr Euler<T, Order> to_euler() const;

    /// Build from Euler angles in sequence @p Order.
    template<EulerOrder Order>
    [[nodiscard]] static constexpr DCM from_euler(const Euler<T, Order>& e);

    /// Build from a (unit) quaternion.
    [[nodiscard]] static constexpr DCM from_quaternion(const Quaternion<T>& q);

    /**
     * @brief Rodrigues formula from axis-angle
     *
     * @p eps is a linear tolerance on the axis length: the axis is normalized
     * internally, so any nonzero axis is valid. The guard rejects a near-zero
     * axis via @f$ \|a\|^2 \le \varepsilon^2 @f$,
     * so short axes from iterative solvers are kept (see @c kinematics/stewart.hpp).
     *
     * @param axis  Rotation axis (normalized internally)
     * @param angle Rotation angle [rad]
     * @param eps   Linear axis-length floor
     * @return DCM, or nullopt if the axis is (near) zero
     */
    [[nodiscard]] static constexpr damp::optional<DCM> from_axis_angle(const Vec3<T>& axis, T angle, T eps = static_cast<T>(1e-9)) {
        T axis_norm2 = (axis[0] * axis[0]) + (axis[1] * axis[1]) + (axis[2] * axis[2]);
        if (axis_norm2 <= (eps * eps)) {
            return damp::nullopt;
        }
        T inv_norm = T{1} / damp::sqrt(axis_norm2);
        T ux = axis[0] * inv_norm;
        T uy = axis[1] * inv_norm;
        T uz = axis[2] * inv_norm;

        const auto [s, c] = damp::sincos(angle);
        T t = T{1} - c;

        return DCM{
            {(t * ux * ux) + c, (t * ux * uy) - (s * uz), (t * ux * uz) + (s * uy)},
            {(t * ux * uy) + (s * uz), (t * uy * uy) + c, (t * uy * uz) - (s * ux)},
            {(t * ux * uz) - (s * uy), (t * uy * uz) + (s * ux), (t * uz * uz) + c},
        };
    }
};

/**
 * @brief Intrinsic Euler angles for sequence @p Order (default @c ZYX aerospace YPR)
 *
 * Storage is @c angle1, @c angle2, @c angle3 along the sequence axes. For
 * @c ZYX: yaw/pitch/roll; for @c XYZ: use @c roll_xyz / @c pitch_xyz / @c yaw_xyz.
 *
 * @tparam T     Floating-point scalar
 * @tparam Order Intrinsic rotation sequence
 */
template<typename T, EulerOrder Order = EulerOrder::ZYX>
struct Euler {

    static constexpr EulerOrder order = Order; ///< Compile-time sequence

    T angle1{}; ///< First rotation about the sequence's first axis [rad]
    T angle2{}; ///< Second rotation [rad]
    T angle3{}; ///< Third rotation [rad]

    constexpr Euler() = default;
    constexpr Euler(T a1, T a2, T a3) : angle1(a1), angle2(a2), angle3(a3) {}
    constexpr Euler(const Euler&) = default;
    constexpr Euler& operator=(const Euler&) = default;
    constexpr Euler(Euler&&) = default;
    constexpr Euler& operator=(Euler&&) = default;
    constexpr ~Euler() = default;

    // --- ZYX (aerospace yaw-pitch-roll): angle1=yaw, angle2=pitch, angle3=roll ---

    /// Yaw (ZYX only) [rad].
    [[nodiscard]] constexpr T yaw() const
        requires(Order == EulerOrder::ZYX)
    { return angle1; }
    /// Yaw (ZYX only) [rad], mutable.
    [[nodiscard]] constexpr T& yaw()
        requires(Order == EulerOrder::ZYX)
    { return angle1; }
    /// Pitch (ZYX only) [rad].
    [[nodiscard]] constexpr T pitch() const
        requires(Order == EulerOrder::ZYX)
    { return angle2; }
    /// Pitch (ZYX only) [rad], mutable.
    [[nodiscard]] constexpr T& pitch()
        requires(Order == EulerOrder::ZYX)
    { return angle2; }
    /// Roll (ZYX only) [rad].
    [[nodiscard]] constexpr T roll() const
        requires(Order == EulerOrder::ZYX)
    { return angle3; }
    /// Roll (ZYX only) [rad], mutable.
    [[nodiscard]] constexpr T& roll()
        requires(Order == EulerOrder::ZYX)
    { return angle3; }

    // --- XYZ (robotics roll-pitch-yaw): angle1=roll, angle2=pitch, angle3=yaw ---

    /// Roll (XYZ only) [rad].
    [[nodiscard]] constexpr T roll_xyz() const
        requires(Order == EulerOrder::XYZ)
    { return angle1; }
    /// Roll (XYZ only) [rad], mutable.
    [[nodiscard]] constexpr T& roll_xyz()
        requires(Order == EulerOrder::XYZ)
    { return angle1; }
    /// Pitch (XYZ only) [rad].
    [[nodiscard]] constexpr T pitch_xyz() const
        requires(Order == EulerOrder::XYZ)
    { return angle2; }
    /// Pitch (XYZ only) [rad], mutable.
    [[nodiscard]] constexpr T& pitch_xyz()
        requires(Order == EulerOrder::XYZ)
    { return angle2; }
    /// Yaw (XYZ only) [rad].
    [[nodiscard]] constexpr T yaw_xyz() const
        requires(Order == EulerOrder::XYZ)
    { return angle3; }
    /// Yaw (XYZ only) [rad], mutable.
    [[nodiscard]] constexpr T& yaw_xyz()
        requires(Order == EulerOrder::XYZ)
    { return angle3; }

    /// Compose principal rotations into a DCM for this sequence.
    [[nodiscard]] constexpr DCM<T> to_dcm() const {
        if constexpr (Order == EulerOrder::ZYX) {
            return DCM<T>::rotate_z(angle1) * DCM<T>::rotate_y(angle2) * DCM<T>::rotate_x(angle3);
        } else if constexpr (Order == EulerOrder::XYZ) {
            return DCM<T>::rotate_x(angle1) * DCM<T>::rotate_y(angle2) * DCM<T>::rotate_z(angle3);
        } else if constexpr (Order == EulerOrder::ZXY) {
            return DCM<T>::rotate_z(angle1) * DCM<T>::rotate_x(angle2) * DCM<T>::rotate_y(angle3);
        } else if constexpr (Order == EulerOrder::YXZ) {
            return DCM<T>::rotate_y(angle1) * DCM<T>::rotate_x(angle2) * DCM<T>::rotate_z(angle3);
        } else if constexpr (Order == EulerOrder::YZX) {
            return DCM<T>::rotate_y(angle1) * DCM<T>::rotate_z(angle2) * DCM<T>::rotate_x(angle3);
        } else if constexpr (Order == EulerOrder::XZY) {
            return DCM<T>::rotate_x(angle1) * DCM<T>::rotate_z(angle2) * DCM<T>::rotate_y(angle3);
        }
    }

    /// Convert to a unit quaternion (via DCM).
    [[nodiscard]] constexpr Quaternion<T> to_quaternion() const;

    /**
     * @brief Extract Euler angles from a DCM (sequence @p Order)
     *
     * Near gimbal lock the third angle is set to 0 and the first absorbs the
     * remaining freedom.
     */
    [[nodiscard]] static constexpr Euler from_dcm(const DCM<T>& R) {
        const auto clamp_unit = [](T sp) {
            if (sp > T{1}) {
                return T{1};
            }
            if (sp < T{-1}) {
                return T{-1};
            }
            return sp;
        };
        const T gimbal = T{1} - static_cast<T>(1e-6);

        if constexpr (Order == EulerOrder::ZYX) {
            // Yaw-Pitch-Roll
            const T sp = clamp_unit(-R(2, 0));
            const T pitch = damp::asin(sp);
            if (damp::abs(sp) < gimbal) {
                return Euler{damp::atan2(R(1, 0), R(0, 0)), pitch, damp::atan2(R(2, 1), R(2, 2))};
            }
            return Euler{damp::atan2(-R(0, 1), R(1, 1)), pitch, T{0}};
        } else if constexpr (Order == EulerOrder::XYZ) {
            // Roll-Pitch-Yaw
            const T sp = clamp_unit(R(0, 2));
            const T pitch = damp::asin(sp);
            if (damp::abs(sp) < gimbal) {
                return Euler{damp::atan2(-R(1, 2), R(2, 2)), pitch, damp::atan2(-R(0, 1), R(0, 0))};
            }
            return Euler{damp::atan2(R(2, 1), R(1, 1)), pitch, T{0}};
        } else if constexpr (Order == EulerOrder::ZXY) {
            const T sp = clamp_unit(R(2, 1));
            const T a2 = damp::asin(sp);
            if (damp::abs(sp) < gimbal) {
                return Euler{damp::atan2(-R(0, 1), R(1, 1)), a2, damp::atan2(-R(2, 0), R(2, 2))};
            }
            return Euler{damp::atan2(R(1, 0), R(0, 0)), a2, T{0}};
        } else if constexpr (Order == EulerOrder::YXZ) {
            const T sp = clamp_unit(-R(1, 2));
            const T a2 = damp::asin(sp);
            if (damp::abs(sp) < gimbal) {
                return Euler{damp::atan2(R(0, 2), R(2, 2)), a2, damp::atan2(R(1, 0), R(1, 1))};
            }
            return Euler{damp::atan2(-R(2, 0), R(0, 0)), a2, T{0}};
        } else if constexpr (Order == EulerOrder::YZX) {
            const T sp = clamp_unit(R(1, 0));
            const T a2 = damp::asin(sp);
            if (damp::abs(sp) < gimbal) {
                return Euler{damp::atan2(-R(2, 0), R(0, 0)), a2, damp::atan2(-R(1, 2), R(1, 1))};
            }
            return Euler{damp::atan2(R(0, 2), R(2, 2)), a2, T{0}};
        } else if constexpr (Order == EulerOrder::XZY) {
            const T sp = clamp_unit(-R(0, 1));
            const T a2 = damp::asin(sp);
            if (damp::abs(sp) < gimbal) {
                return Euler{damp::atan2(R(2, 1), R(1, 1)), a2, damp::atan2(R(0, 2), R(0, 0))};
            }
            return Euler{damp::atan2(-R(1, 2), R(2, 2)), a2, T{0}};
        }
    }

    /// Build from a unit quaternion (via DCM).
    [[nodiscard]] static constexpr Euler from_quaternion(const Quaternion<T>& q);
};

/// Aerospace yaw-pitch-roll (ψ, θ, φ); do not conflate with @ref EulerXYZ.
template<typename T>
using EulerZYX = Euler<T, EulerOrder::ZYX>;

/// Robotics roll-pitch-yaw; do not conflate with @ref EulerZYX.
template<typename T>
using EulerXYZ = Euler<T, EulerOrder::XYZ>;

/**
 * @brief Unit quaternion rotation (w, x, y, z) (Hamilton product)
 *
 * Stored as a 4×1 matrix with scalar part @c w first. Prefer @ref Pose for
 * rigid-body frames; use this type for pure orientation.
 *
 * @tparam T Floating-point scalar
 */
template<typename T>
struct Quaternion : public Matrix<4, 1, T> {

    static_assert(std::is_floating_point_v<T>, "Quaternion element type must be floating point");

    /// Scalar (real) part.
    constexpr T& w() { return this->data_[0]; }
    /// Scalar (real) part.
    constexpr const T& w() const { return this->data_[0]; }
    /// Vector part x.
    constexpr T& x() { return this->data_[1]; }
    /// Vector part x.
    constexpr const T& x() const { return this->data_[1]; }
    /// Vector part y.
    constexpr T& y() { return this->data_[2]; }
    /// Vector part y.
    constexpr const T& y() const { return this->data_[2]; }
    /// Vector part z.
    constexpr T& z() { return this->data_[3]; }
    /// Vector part z.
    constexpr const T& z() const { return this->data_[3]; }

    /// Identity quaternion (1, 0, 0, 0).
    constexpr Quaternion() : Matrix<4, 1, T>() { w() = T{1}; }
    constexpr Quaternion(const Quaternion&) = default;
    constexpr Quaternion& operator=(const Quaternion&) = default;
    constexpr Quaternion(Quaternion&&) = default;
    constexpr Quaternion& operator=(Quaternion&&) = default;
    constexpr ~Quaternion() = default;

    /// Construct from components (w, x, y, z); call @ref normalize if needed.
    constexpr Quaternion(T w_, T x_, T y_, T z_) : Matrix<4, 1, T>() {
        w() = w_;
        x() = x_;
        y() = y_;
        z() = z_;
    }

    /// Construct from up to four elements in (w, x, y, z) order.
    constexpr Quaternion(std::initializer_list<T> values) : Matrix<4, 1, T>() {
        size_t i = 0;
        for (const auto& val : values) {
            if (i < 4) {
                this->data_[i] = val;
            }
            ++i;
        }
    }

    /// Construct from a span of up to four elements in (w, x, y, z) order.
    template<typename SpanType>
        requires std::is_same_v<SpanType, std::span<const T>>
    constexpr explicit Quaternion(SpanType values) : Matrix<4, 1, T>() {
        size_t i = 0;
        for (const auto& val : values) {
            if (i < 4) {
                this->data_[i] = val;
            }
            ++i;
        }
    }

    /// Convert from another scalar precision.
    template<typename U>
    constexpr explicit Quaternion(const Quaternion<U>& other) : Matrix<4, 1, T>(other) {}

    /// Convert from a 4×1 matrix viewed as (w, x, y, z).
    template<typename U>
    constexpr explicit Quaternion(const Matrix<4, 1, U>& other) : Matrix<4, 1, T>(other) {}

    /// Assign from a 4×1 matrix viewed as (w, x, y, z).
    template<typename U>
    constexpr Quaternion& operator=(const Matrix<4, 1, U>& other) {
        Matrix<4, 1, T>::operator=(other);
        return *this;
    }

    /// Identity (no rotation).
    [[nodiscard]] static constexpr Quaternion identity() { return Quaternion{T{1}, T{0}, T{0}, T{0}}; }

    /// Squared Euclidean norm ‖q‖².
    [[nodiscard]] constexpr T norm_squared() const { return (w() * w()) + (x() * x()) + (y() * y()) + (z() * z()); }
    /// Euclidean norm ‖q‖.
    [[nodiscard]] constexpr T norm() const { return damp::sqrt(norm_squared()); }

    /**
     * @brief Unit quaternion, or nullopt if ‖q‖² ≤ ε
     * @param eps Squared-norm floor for the near-zero test
     */
    [[nodiscard]] constexpr damp::optional<Quaternion> normalized_safe(T eps = static_cast<T>(1e-9)) const {
        T n2 = norm_squared();
        if (n2 <= eps) {
            return damp::nullopt;
        }
        T inv_n = T{1} / damp::sqrt(n2);
        return Quaternion{w() * inv_n, x() * inv_n, y() * inv_n, z() * inv_n};
    }

    /**
     * @brief Normalize in place
     * @return false if the quaternion is near zero (left unchanged)
     */
    constexpr bool normalize_in_place(T eps = static_cast<T>(1e-9)) {
        T n2 = norm_squared();
        if (n2 <= eps) {
            return false;
        }
        T inv_n = T{1} / damp::sqrt(n2);
        w() *= inv_n;
        x() *= inv_n;
        y() *= inv_n;
        z() *= inv_n;
        return true;
    }

    /// Unit quaternion, or @c *this unchanged if near zero.
    [[nodiscard]] constexpr Quaternion normalized(T eps = static_cast<T>(1e-9)) const {
        auto n = normalized_safe(eps);
        return n.value_or(*this);
    }

    /// Conjugate (w, −x, −y, −z) (inverse for unit quaternions).
    [[nodiscard]] constexpr Quaternion conjugate() const { return Quaternion{w(), -x(), -y(), -z()}; }

    /**
     * @brief Multiplicative inverse q⁻¹ = q̄ / ‖q‖²
     * @return Inverse, or nullopt if near zero
     */
    [[nodiscard]] constexpr damp::optional<Quaternion> inverse(T eps = static_cast<T>(1e-9)) const {
        T n2 = norm_squared();
        if (n2 <= eps) {
            return damp::nullopt;
        }
        T inv_n2 = T{1} / n2;
        return Quaternion{w() * inv_n2, -x() * inv_n2, -y() * inv_n2, -z() * inv_n2};
    }

    /// Hamilton product.
    [[nodiscard]] constexpr Quaternion operator*(const Quaternion& rhs) const {
        return Quaternion{
            (w() * rhs.w()) - (x() * rhs.x()) - (y() * rhs.y()) - (z() * rhs.z()),
            (w() * rhs.x()) + (x() * rhs.w()) + (y() * rhs.z()) - (z() * rhs.y()),
            (w() * rhs.y()) - (x() * rhs.z()) + (y() * rhs.w()) + (z() * rhs.x()),
            (w() * rhs.z()) + (x() * rhs.y()) - (y() * rhs.x()) + (z() * rhs.w()),
        };
    }

    /// In-place Hamilton product.
    constexpr Quaternion& operator*=(const Quaternion& rhs) {
        return *this = (*this) * rhs;
    }

    /// Scale all components.
    [[nodiscard]] constexpr Quaternion operator*(T scalar) const {
        return Quaternion{w() * scalar, x() * scalar, y() * scalar, z() * scalar};
    }

    /// Divide all components by a scalar.
    [[nodiscard]] constexpr Quaternion operator/(T scalar) const {
        return Quaternion{w() / scalar, x() / scalar, y() / scalar, z() / scalar};
    }

    /// Rotate a 3-vector by this orientation (efficient sandwich form).
    [[nodiscard]] constexpr Vec3<T> rotate(const Vec3<T>& v) const {
        const Vec3<T> qv{x(), y(), z()};
        const T       t0 = T{2} * ((qv[1] * v[2]) - (qv[2] * v[1]));
        const T       t1 = T{2} * ((qv[2] * v[0]) - (qv[0] * v[2]));
        const T       t2 = T{2} * ((qv[0] * v[1]) - (qv[1] * v[0]));
        return Vec3<T>{
            v[0] + (w() * t0) + ((qv[1] * t2) - (qv[2] * t1)),
            v[1] + (w() * t1) + ((qv[2] * t0) - (qv[0] * t2)),
            v[2] + (w() * t2) + ((qv[0] * t1) - (qv[1] * t0)),
        };
    }

    /// Convert to a DCM (normalizes first).
    [[nodiscard]] constexpr DCM<T> to_dcm() const {
        const Quaternion qn = normalized();

        const T ww = qn.w();
        const T xx = qn.x();
        const T yy = qn.y();
        const T zz = qn.z();

        return DCM<T>{
            {T{1} - (T{2} * ((yy * yy) + (zz * zz))), T{2} * ((xx * yy) - (zz * ww)), T{2} * ((xx * zz) + (yy * ww))},
            {T{2} * ((xx * yy) + (zz * ww)), T{1} - (T{2} * ((xx * xx) + (zz * zz))), T{2} * ((yy * zz) - (xx * ww))},
            {T{2} * ((xx * zz) - (yy * ww)), T{2} * ((yy * zz) + (xx * ww)), T{1} - (T{2} * ((xx * xx) + (yy * yy)))},
        };
    }

    /// Convert to Euler angles (default aerospace @c ZYX).
    template<EulerOrder Order = EulerOrder::ZYX>
    [[nodiscard]] constexpr Euler<T, Order> to_euler() const {
        return Euler<T, Order>::from_dcm(to_dcm());
    }

    /**
     * @brief Build a unit quaternion from a DCM (Shepperd / largest-component branch)
     * @param R   Rotation matrix
     * @param eps Normalization / branch floor
     * @return Unit quaternion, or nullopt if a branch is singular
     */
    [[nodiscard]] static constexpr damp::optional<Quaternion> from_dcm(const DCM<T>& R, T eps = static_cast<T>(1e-6)) {
        const T half = static_cast<T>(0.5);
        const T trace = R(0, 0) + R(1, 1) + R(2, 2);
        if (trace > T{0}) {
            const T s = damp::sqrt(trace + T{1});
            const T inv_s = half / s;
            return Quaternion{
                half * s,
                (R(2, 1) - R(1, 2)) * inv_s,
                (R(0, 2) - R(2, 0)) * inv_s,
                (R(1, 0) - R(0, 1)) * inv_s,
            }
                .normalized_safe(eps);
        }
        if (R(0, 0) > R(1, 1) && R(0, 0) > R(2, 2)) {
            const T s = damp::sqrt(T{1} + R(0, 0) - R(1, 1) - R(2, 2));
            if (s <= eps) {
                return damp::nullopt;
            }
            const T inv_s = half / s;
            return Quaternion{
                (R(2, 1) - R(1, 2)) * inv_s,
                half * s,
                (R(0, 1) + R(1, 0)) * inv_s,
                (R(0, 2) + R(2, 0)) * inv_s,
            }
                .normalized_safe(eps);
        }
        if (R(1, 1) > R(2, 2)) {
            const T s = damp::sqrt(T{1} + R(1, 1) - R(0, 0) - R(2, 2));
            if (s <= eps) {
                return damp::nullopt;
            }
            const T inv_s = half / s;
            return Quaternion{
                (R(0, 2) - R(2, 0)) * inv_s,
                (R(0, 1) + R(1, 0)) * inv_s,
                half * s,
                (R(1, 2) + R(2, 1)) * inv_s,
            }
                .normalized_safe(eps);
        }
        const T s = damp::sqrt(T{1} + R(2, 2) - R(0, 0) - R(1, 1));
        if (s <= eps) {
            return damp::nullopt;
        }
        const T inv_s = half / s;
        return Quaternion{
            (R(1, 0) - R(0, 1)) * inv_s,
            (R(0, 2) + R(2, 0)) * inv_s,
            (R(1, 2) + R(2, 1)) * inv_s,
            half * s,
        }
            .normalized_safe(eps);
    }

    /// Build from Euler angles (identity if the DCM conversion fails).
    template<EulerOrder Order>
    [[nodiscard]] static constexpr Quaternion from_euler(const Euler<T, Order>& e) {
        auto q_opt = from_dcm(e.to_dcm());
        return q_opt.value_or(identity());
    }

    /**
     * @brief Unit quaternion from axis-angle
     *
     * @p eps is a linear axis-length floor (@f$ \|a\|^2 \le \varepsilon^2 @f$);
     * see @ref DCM::from_axis_angle for the squared-norm rationale.
     *
     * @param axis  Rotation axis (normalized internally)
     * @param angle Rotation angle [rad]
     * @param eps   Linear axis-length floor
     * @return Unit quaternion, or nullopt if the axis is (near) zero
     */
    [[nodiscard]] static constexpr damp::optional<Quaternion> from_axis_angle(const Vec3<T>& axis, T angle, T eps = static_cast<T>(1e-9)) {
        T axis_norm2 = (axis[0] * axis[0]) + (axis[1] * axis[1]) + (axis[2] * axis[2]);
        if (axis_norm2 <= (eps * eps)) {
            return damp::nullopt;
        }
        T inv_axis_norm = T{1} / damp::sqrt(axis_norm2);
        T half = angle * static_cast<T>(0.5);
        const auto [s_half, c_half] = damp::sincos(half);
        T s = s_half * inv_axis_norm;
        return Quaternion{c_half, axis[0] * s, axis[1] * s, axis[2] * s};
    }

    /**
     * @brief Axis-angle extraction (inverse of @ref from_axis_angle)
     *
     * Returns a unit axis and angle in @f$ [0, \pi] @f$ via
     * @f$ 2 \operatorname{atan2}(\|v\|, w) @f$ (no prior normalize required).
     * @c q and @c -q are the same rotation; @c w is forced non-negative first.
     * Near-zero rotation defaults the axis to @f$ (1,0,0) @f$ with angle 0.
     *
     * @param eps Linear floor on the vector-part length
     */
    [[nodiscard]] constexpr damp::pair<Vec3<T>, T> to_axis_angle(T eps = static_cast<T>(1e-9)) const {
        // q and -q are the same rotation; canonicalize to w ≥ 0 so angle ∈ [0, π].
        T qw = w();
        T qx = x();
        T qy = y();
        T qz = z();
        if (qw < T{0}) {
            qw = -qw;
            qx = -qx;
            qy = -qy;
            qz = -qz;
        }
        T v2 = (qx * qx) + (qy * qy) + (qz * qz);
        if (v2 <= (eps * eps)) {
            return {Vec3<T>{T{1}, T{0}, T{0}}, T{0}};
        }
        T vnorm = damp::sqrt(v2);
        T angle = T{2} * damp::atan2(vnorm, qw);
        T inv = T{1} / vnorm;
        return {Vec3<T>{qx * inv, qy * inv, qz * inv}, angle};
    }

    /**
     * @brief Log map to a rotation vector (axis · angle) [rad]
     *
     * Inverse of the axis-angle construction used by @ref from_axis_angle.
     * @f$ \exp(\log q) = \pm q @f$ (same rotation).
     */
    [[nodiscard]] constexpr Vec3<T> log(T eps = static_cast<T>(1e-9)) const {
        auto [axis, angle] = to_axis_angle(eps);
        return Vec3<T>{axis[0] * angle, axis[1] * angle, axis[2] * angle};
    }

    /**
     * @brief First-order body-rate integration: q ← q ⊗ (1, ½ ω Δt)
     * @param omega Body angular rate [rad/s]
     * @param dt    Step [s]
     */
    [[nodiscard]] constexpr Quaternion integrate_body_rates(const Vec3<T>& omega, T dt) const {
        const T half_dt = dt * static_cast<T>(0.5);
        return ((*this) * Quaternion{T{1}, omega[0] * half_dt, omega[1] * half_dt, omega[2] * half_dt})
            .normalized();
    }

    /**
     * @brief Spherical linear interpolation on the unit sphere (shortest arc)
     * @param a Start orientation
     * @param b End orientation
     * @param t Blend in [0,1] (clamped)
     */
    [[nodiscard]] static constexpr Quaternion slerp(const Quaternion& a, const Quaternion& b, T t) {
        T tt = t;
        if (tt < T{0}) {
            tt = T{0};
        }
        if (tt > T{1}) {
            tt = T{1};
        }

        T cos_theta = (a.w() * b.w()) + (a.x() * b.x()) + (a.y() * b.y()) + (a.z() * b.z());

        // Same hemisphere as a (shortest arc)
        const T bw = (cos_theta < T{0}) ? -b.w() : b.w();
        const T bx = (cos_theta < T{0}) ? -b.x() : b.x();
        const T by = (cos_theta < T{0}) ? -b.y() : b.y();
        const T bz = (cos_theta < T{0}) ? -b.z() : b.z();
        if (cos_theta < T{0}) {
            cos_theta = -cos_theta;
        }

        const T k_small = static_cast<T>(1e-6);
        if (cos_theta > (T{1} - k_small)) {
            return Quaternion{
                a.w() + (tt * (bw - a.w())),
                a.x() + (tt * (bx - a.x())),
                a.y() + (tt * (by - a.y())),
                a.z() + (tt * (bz - a.z())),
            }
                .normalized();
        }

        const T theta = damp::acos(cos_theta);
        const T sin_theta = damp::sin(theta);
        const T w1 = damp::sin((T{1} - tt) * theta) / sin_theta;
        const T w2 = damp::sin(tt * theta) / sin_theta;
        return Quaternion{
            (a.w() * w1) + (bw * w2),
            (a.x() * w1) + (bx * w2),
            (a.y() * w1) + (by * w2),
            (a.z() * w1) + (bz * w2),
        }
            .normalized();
    }
};

/**
 * @brief Left scalar multiply for @ref Quaternion
 * @related Quaternion
 */
template<typename T>
[[nodiscard]] constexpr Quaternion<T> operator*(T scalar, const Quaternion<T>& q) {
    return q * scalar;
}

// Out-of-line members defined after all rotation types are complete.

template<typename T>
constexpr damp::optional<Quaternion<T>> DCM<T>::to_quaternion(T eps) const {
    return Quaternion<T>::from_dcm(*this, eps);
}

template<typename T>
template<EulerOrder Order>
constexpr Euler<T, Order> DCM<T>::to_euler() const {
    return Euler<T, Order>::from_dcm(*this);
}

template<typename T>
template<EulerOrder Order>
constexpr DCM<T> DCM<T>::from_euler(const Euler<T, Order>& e) {
    return e.to_dcm();
}

template<typename T>
constexpr DCM<T> DCM<T>::from_quaternion(const Quaternion<T>& q) {
    return q.to_dcm();
}

template<typename T, EulerOrder Order>
constexpr Quaternion<T> Euler<T, Order>::to_quaternion() const {
    return Quaternion<T>::from_euler(*this);
}

template<typename T, EulerOrder Order>
constexpr Euler<T, Order> Euler<T, Order>::from_quaternion(const Quaternion<T>& q) {
    return q.template to_euler<Order>();
}

/**
 * @brief 4×4 homogeneous transform (SE(3) interop / DH export)
 *
 * Not the working pose representation — use @ref Pose (@c Quaternion + translation)
 * for composition. Convert with @c Pose::to_transform4 / @c Pose::from_transform4.
 *
 * @tparam T Floating-point scalar
 */
template<typename T>
struct Transform4 : public Mat4<T> {

    static_assert(std::is_floating_point_v<T>, "Transform4 element type must be floating point");

    constexpr Transform4() : Mat4<T>(Mat4<T>::identity()) {}
    constexpr Transform4(const Transform4&) = default;
    constexpr Transform4& operator=(const Transform4&) = default;
    constexpr Transform4(Transform4&&) = default;
    constexpr Transform4& operator=(Transform4&&) = default;
    constexpr ~Transform4() = default;

    /// Construct from a raw @ref Mat4.
    constexpr explicit Transform4(const Mat4<T>& m) : Mat4<T>(m) {}

    /// Row-nested initializer list (prvalue / RVO-friendly).
    constexpr Transform4(std::initializer_list<std::initializer_list<T>> init) : Mat4<T>(init) {}

    /// Identity transform.
    [[nodiscard]] static constexpr Transform4 identity() { return Transform4{}; }

    /// Build from rotation matrix and translation.
    [[nodiscard]] static constexpr Transform4 from_rotation_translation(const DCM<T>& R, const Vec3<T>& t) {
        return Transform4{
            {R(0, 0), R(0, 1), R(0, 2), t[0]},
            {R(1, 0), R(1, 1), R(1, 2), t[1]},
            {R(2, 0), R(2, 1), R(2, 2), t[2]},
            {T{0}, T{0}, T{0}, T{1}},
        };
    }

    /// Build from quaternion orientation and translation.
    [[nodiscard]] static constexpr Transform4 from_quaternion_translation(const Quaternion<T>& q, const Vec3<T>& t) {
        return from_rotation_translation(q.to_dcm(), t);
    }

    /// Build from Euler orientation and translation.
    template<EulerOrder Order>
    [[nodiscard]] static constexpr Transform4 from_euler_translation(const Euler<T, Order>& e, const Vec3<T>& t) {
        return from_rotation_translation(e.to_dcm(), t);
    }

    /// Upper-left 3×3 rotation block as a @ref DCM.
    [[nodiscard]] constexpr DCM<T> rotation() const {
        return DCM<T>{
            {this->data_[(0 * 4) + 0], this->data_[(0 * 4) + 1], this->data_[(0 * 4) + 2]},
            {this->data_[(1 * 4) + 0], this->data_[(1 * 4) + 1], this->data_[(1 * 4) + 2]},
            {this->data_[(2 * 4) + 0], this->data_[(2 * 4) + 1], this->data_[(2 * 4) + 2]},
        };
    }

    /// Translation column (t_x, t_y, t_z).
    [[nodiscard]] constexpr Vec3<T> translation() const {
        return Vec3<T>{this->data_[(0 * 4) + 3], this->data_[(1 * 4) + 3], this->data_[(2 * 4) + 3]};
    }

    /// Compose transforms: T_this T_rhs.
    [[nodiscard]] constexpr Transform4 operator*(const Transform4& rhs) const {
        return Transform4(static_cast<const Mat4<T>&>(*this) * static_cast<const Mat4<T>&>(rhs));
    }

    /// In-place composition with @p rhs.
    constexpr Transform4& operator*=(const Transform4& rhs) {
        return *this = *this * rhs;
    }

    /// Map a point with rotation and translation (homogeneous; divides by @c w if needed).
    [[nodiscard]] constexpr Vec3<T> transform_point(const Vec3<T>& p) const {
        T x = (this->data_[(0 * 4) + 0] * p[0]) + (this->data_[(0 * 4) + 1] * p[1])
            + (this->data_[(0 * 4) + 2] * p[2]) + this->data_[(0 * 4) + 3];
        T y = (this->data_[(1 * 4) + 0] * p[0]) + (this->data_[(1 * 4) + 1] * p[1])
            + (this->data_[(1 * 4) + 2] * p[2]) + this->data_[(1 * 4) + 3];
        T z = (this->data_[(2 * 4) + 0] * p[0]) + (this->data_[(2 * 4) + 1] * p[1])
            + (this->data_[(2 * 4) + 2] * p[2]) + this->data_[(2 * 4) + 3];
        T w = (this->data_[(3 * 4) + 0] * p[0]) + (this->data_[(3 * 4) + 1] * p[1])
            + (this->data_[(3 * 4) + 2] * p[2]) + this->data_[(3 * 4) + 3];

        if (w != T{1} && w != T{0}) {
            T inv_w = T{1} / w;
            x *= inv_w;
            y *= inv_w;
            z *= inv_w;
        }

        return Vec3<T>{x, y, z};
    }

    /// Map a free vector (rotation only; no translation).
    [[nodiscard]] constexpr Vec3<T> transform_vector(const Vec3<T>& v) const {
        return rotation() * v;
    }

    /**
     * @brief Rigid inverse [R t; 0 1]⁻¹ = [Rᵀ, −Rᵀt; 0 1]
     *
     * Analytic for a rigid transform (no 4×4 LU). Wrapped in @c optional so the
     * signature matches @c Mat4::inverse(); always succeeds for finite input.
     */
    [[nodiscard]] constexpr damp::optional<Transform4> inverse() const {
        const DCM<T> Rt = rotation().transpose();
        return from_rotation_translation(Rt, -(Rt * translation()));
    }

    /// Split into unit quaternion and translation (identity quat if DCM conversion fails).
    [[nodiscard]] constexpr damp::pair<Quaternion<T>, Vec3<T>> to_quaternion_translation() const {
        return {
            rotation().to_quaternion().value_or(Quaternion<T>::identity()),
            translation(),
        };
    }

    /// Split into Euler angles and translation (default aerospace @c ZYX).
    template<EulerOrder Order = EulerOrder::ZYX>
    [[nodiscard]] constexpr damp::pair<Euler<T, Order>, Vec3<T>> to_euler_translation() const {
        return {rotation().template to_euler<Order>(), translation()};
    }

    /// Underlying @ref Mat4 view.
    [[nodiscard]] constexpr const Mat4<T>& matrix() const { return *this; }
    /// Underlying @ref Mat4 view (mutable).
    [[nodiscard]] constexpr Mat4<T>& matrix() { return *this; }
};

using Transform4f = Transform4<float>;  ///< @c Transform4 with @c float
using Transform4d = Transform4<double>; ///< @c Transform4 with @c double

/// Alias: Special Euclidean group SE(3) ≡ @ref Transform4.
template<typename T>
using SE3 = Transform4<T>;
using SE3f = Transform4<float>;  ///< @c SE3 with @c float
using SE3d = Transform4<double>; ///< @c SE3 with @c double

using Quatf = Quaternion<float>;  ///< @c Quaternion with @c float
using Quatd = Quaternion<double>; ///< @c Quaternion with @c double

using Vec3f = Vec3<float>;  ///< @c Vec3 with @c float
using Vec3d = Vec3<double>; ///< @c Vec3 with @c double

using DCMf = DCM<float>;  ///< @c DCM with @c float
using DCMd = DCM<double>; ///< @c DCM with @c double

/// Alias: Special Orthogonal group SO(3) ≡ @ref DCM.
template<typename T>
using SO3 = DCM<T>;
using SO3f = DCM<float>;  ///< @c SO3 with @c float
using SO3d = DCM<double>; ///< @c SO3 with @c double

using EulerZYXf = EulerZYX<float>;  ///< Aerospace YPR, @c float
using EulerZYXd = EulerZYX<double>; ///< Aerospace YPR, @c double

using EulerXYZf = EulerXYZ<float>;  ///< Robotics RPY, @c float
using EulerXYZd = EulerXYZ<double>; ///< Robotics RPY, @c double

} // namespace damp
