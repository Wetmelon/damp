// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <numbers>

#include "damp/backend.hpp"
#include "damp/math/geometry.hpp"
#include "damp/matrix/matrix.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

constexpr float kPi = std::numbers::pi_v<float>;

TEST_SUITE("Quaternion") {
    TEST_CASE("Identity and normalization") {
        Quaternion<float> q_id = Quaternion<float>::identity();
        // Exact float equality is fragile under -ffast-math; Approx is the suite norm.
        CHECK(q_id.w() == doctest::Approx(1.0f));
        CHECK(q_id.x() == doctest::Approx(0.0f));
        CHECK(q_id.y() == doctest::Approx(0.0f));
        CHECK(q_id.z() == doctest::Approx(0.0f));

        auto q_norm = q_id.normalized();
        CHECK(q_norm.w() == doctest::Approx(1.0f));
        CHECK(q_norm.x() == doctest::Approx(0.0f));

        Quaternion<float> q_zero{0.0f, 0.0f, 0.0f, 0.0f};
        auto              maybe_norm = q_zero.normalized_safe();
        CHECK_FALSE(maybe_norm.has_value());
    }

    TEST_CASE("Conjugate and inverse") {
        auto q_opt = Quaternion<float>::from_axis_angle(Vec3<float>{0.0f, 0.0f, 1.0f}, kPi * 0.5f);
        REQUIRE(q_opt.has_value());
        const auto& q = q_opt.value();

        Vec3<float> v{1.0f, 0.0f, 0.0f};
        auto        v_rot = q.rotate(v);
        CHECK(v_rot[0] == doctest::Approx(0.0).epsilon(1e-5));
        CHECK(v_rot[1] == doctest::Approx(1.0).epsilon(1e-5));
        CHECK(v_rot[2] == doctest::Approx(0.0).epsilon(1e-5));

        auto q_inv = q.inverse();
        REQUIRE(q_inv.has_value());
        auto v_restored = q_inv->rotate(v_rot);
        CHECK(v_restored[0] == doctest::Approx(static_cast<double>(v[0])).epsilon(1e-5));
        CHECK(v_restored[1] == doctest::Approx(static_cast<double>(v[1])).epsilon(1e-5));
        CHECK(v_restored[2] == doctest::Approx(static_cast<double>(v[2])).epsilon(1e-5));
    }

    TEST_CASE("Rotation matrix conversion round-trip") {
        // Create quaternion from Euler ZYX (yaw, pitch, roll)
        EulerZYX<float> e_in(kPi * 0.5f, 0.0f, 0.0f); // 90 deg yaw only
        auto            q = e_in.to_quaternion();
        auto            R = q.to_dcm();

        // Expected 90 deg about Z
        CHECK(R(0, 0) == doctest::Approx(0.0).epsilon(1e-6));
        CHECK(R(0, 1) == doctest::Approx(-1.0).epsilon(1e-6));
        CHECK(R(1, 0) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(R(1, 1) == doctest::Approx(0.0).epsilon(1e-6));

        auto q_rt = R.to_quaternion();
        REQUIRE(q_rt.has_value());

        Vec3<float> v{1.0f, 0.0f, 0.0f};
        auto        v_rot_mat = R * v;
        auto        v_rot_q = q_rt->rotate(v);
        CHECK(v_rot_mat[0] == doctest::Approx(static_cast<double>(v_rot_q[0])).epsilon(1e-5));
        CHECK(v_rot_mat[1] == doctest::Approx(static_cast<double>(v_rot_q[1])).epsilon(1e-5));
        CHECK(v_rot_mat[2] == doctest::Approx(static_cast<double>(v_rot_q[2])).epsilon(1e-5));
    }

    TEST_CASE("Integrate body rates (small step)") {
        Quaternion<float> q0 = Quaternion<float>::identity();
        Vec3<float>       omega{0.0f, 0.0f, kPi}; // rad/s about Z
        float             dt = 0.05f;             // small step for first-order approx (~9 deg)
        auto              q1 = q0.integrate_body_rates(omega, dt);

        Vec3<float> v{1.0f, 0.0f, 0.0f};
        auto        v_rot = q1.rotate(v);

        // Expected small-angle rotation using exact trig for reference
        float angle = omega[2] * dt;
        float c = std::cos(angle);
        float s = std::sin(angle);
        CHECK(v_rot[0] == doctest::Approx(static_cast<double>(c)).epsilon(1e-3));
        CHECK(v_rot[1] == doctest::Approx(static_cast<double>(s)).epsilon(1e-3));
        CHECK(v_rot[2] == doctest::Approx(0.0).epsilon(1e-4));
    }

    TEST_CASE("Slerp mid-rotation about Z") {
        auto qa = Quaternion<float>::identity();
        auto qb_opt = Quaternion<float>::from_axis_angle(Vec3<float>{0.0f, 0.0f, 1.0f}, kPi * 0.5f); // 90 deg
        REQUIRE(qb_opt.has_value());
        const auto& qb = qb_opt.value();

        auto        qm = Quaternion<float>::slerp(qa, qb, 0.5f); // ~45 deg about Z
        Vec3<float> v{1.0f, 0.0f, 0.0f};
        auto        v_rot = qm.rotate(v);
        CHECK(v_rot[0] == doctest::Approx(static_cast<double>(std::cos(kPi * 0.25f))).epsilon(1e-5));
        CHECK(v_rot[1] == doctest::Approx(static_cast<double>(std::sin(kPi * 0.25f))).epsilon(1e-5));
        CHECK(v_rot[2] == doctest::Approx(0.0).epsilon(1e-5));
    }

    TEST_CASE("Log / axis-angle round-trips exp (from_axis_angle)") {
        Vec3<double> axis{1.0, 2.0, -0.5};
        double       angle = 0.9; // < π, so log recovers it exactly
        auto         q_opt = Quatd::from_axis_angle(axis, angle);
        REQUIRE(q_opt.has_value());

        // to_axis_angle recovers the (normalized) axis and angle.
        auto [rec_axis, rec_angle] = q_opt.value().to_axis_angle();
        double n = std::sqrt((axis[0] * axis[0]) + (axis[1] * axis[1]) + (axis[2] * axis[2]));
        CHECK(rec_angle == doctest::Approx(angle).epsilon(1e-9));
        CHECK(rec_axis[0] == doctest::Approx(axis[0] / n).epsilon(1e-9));
        CHECK(rec_axis[1] == doctest::Approx(axis[1] / n).epsilon(1e-9));
        CHECK(rec_axis[2] == doctest::Approx(axis[2] / n).epsilon(1e-9));

        // log() = axis·angle; exp(log(q)) == q reproduces the same rotation.
        Vec3<double> rv = q_opt.value().log();
        double       rvn = std::sqrt((rv[0] * rv[0]) + (rv[1] * rv[1]) + (rv[2] * rv[2]));
        auto         q2 = Quatd::from_axis_angle(rv, rvn);
        REQUIRE(q2.has_value());
        Vec3<double> v{0.3, -0.7, 1.1};
        auto         a = q_opt.value().rotate(v);
        auto         b = q2.value().rotate(v);
        CHECK(b[0] == doctest::Approx(a[0]).epsilon(1e-9));
        CHECK(b[1] == doctest::Approx(a[1]).epsilon(1e-9));
        CHECK(b[2] == doctest::Approx(a[2]).epsilon(1e-9));

        // Identity logs to the zero rotation vector.
        Vec3<double> z = Quatd::identity().log();
        CHECK(z[0] == doctest::Approx(0.0));
        CHECK(z[1] == doctest::Approx(0.0));
        CHECK(z[2] == doctest::Approx(0.0));
    }

    TEST_CASE("Euler <-> DCM conversion") {
        EulerZYX<float> e(0.5f, -0.4f, 0.3f); // yaw, pitch, roll
        auto            R = e.to_dcm();
        auto            e_back = EulerZYX<float>::from_dcm(R);

        CHECK(e_back.roll() == doctest::Approx(static_cast<double>(e.roll())).epsilon(1e-5));
        CHECK(e_back.pitch() == doctest::Approx(static_cast<double>(e.pitch())).epsilon(1e-5));
        CHECK(e_back.yaw() == doctest::Approx(static_cast<double>(e.yaw())).epsilon(1e-5));

        // Spot-check a couple matrix entries against direct trig reference
        float cr = std::cos(e.roll());
        float sr = std::sin(e.roll());
        float cp = std::cos(e.pitch());
        float sp = std::sin(e.pitch());
        float cy = std::cos(e.yaw());
        float sy = std::sin(e.yaw());
        CHECK(R(0, 0) == doctest::Approx(static_cast<double>(cy * cp)).epsilon(1e-6));
        CHECK(R(2, 0) == doctest::Approx(static_cast<double>(-sp)).epsilon(1e-6));
        CHECK(R(0, 1) == doctest::Approx(static_cast<double>((cy * sp * sr) - (sy * cr))).epsilon(1e-6));
    }

    TEST_CASE("Quat <-> Euler conversion") {
        EulerZYX<float> e(0.4f, -0.1f, 0.2f); // yaw, pitch, roll
        auto            q = e.to_quaternion();
        auto            e_back = q.to_euler<EulerOrder::ZYX>();

        CHECK(e_back.roll() == doctest::Approx(static_cast<double>(e.roll())).epsilon(1e-5));
        CHECK(e_back.pitch() == doctest::Approx(static_cast<double>(e.pitch())).epsilon(1e-5));
        CHECK(e_back.yaw() == doctest::Approx(static_cast<double>(e.yaw())).epsilon(1e-5));
    }

    TEST_CASE("DCM <-> Quaternion conversion") {
        EulerZYX<float> e(-0.45f, 0.35f, -0.25f); // yaw, pitch, roll
        auto            R = e.to_dcm();
        auto            q_opt = R.to_quaternion();
        REQUIRE(q_opt.has_value());
        auto R_back = q_opt->to_dcm();

        for (size_t r = 0; r < 3; ++r) {
            for (size_t c = 0; c < 3; ++c) {
                CHECK(R_back(r, c) == doctest::Approx(static_cast<double>(R(r, c))).epsilon(1e-5));
            }
        }
    }
}

TEST_SUITE("Transform4") {
    TEST_CASE("Identity and basic construction") {
        Transform4d T_id = Transform4d::identity();
        CHECK(T_id(0, 0) == 1.0);
        CHECK(T_id(1, 1) == 1.0);
        CHECK(T_id(2, 2) == 1.0);
        CHECK(T_id(3, 3) == 1.0);
        CHECK(T_id(0, 3) == 0.0);
        CHECK(T_id(1, 3) == 0.0);
        CHECK(T_id(2, 3) == 0.0);
    }

    TEST_CASE("Construction from rotation and translation") {
        DCMd        R = DCMd::rotate_z(static_cast<double>(kPi) * 0.5); // 90° rotation about Z
        Vec3d       t{1.0, 2.0, 3.0};
        Transform4d T = Transform4d::from_rotation_translation(R, t);

        // Check rotation part
        CHECK(T(0, 0) == doctest::Approx(0.0).epsilon(1e-6));
        CHECK(T(0, 1) == doctest::Approx(-1.0).epsilon(1e-6));
        CHECK(T(1, 0) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(T(1, 1) == doctest::Approx(0.0).epsilon(1e-6));

        // Check translation part
        CHECK(T(0, 3) == 1.0);
        CHECK(T(1, 3) == 2.0);
        CHECK(T(2, 3) == 3.0);

        // Check homogeneous row
        CHECK(T(3, 0) == 0.0);
        CHECK(T(3, 1) == 0.0);
        CHECK(T(3, 2) == 0.0);
        CHECK(T(3, 3) == 1.0);
    }

    TEST_CASE("Point transformation") {
        DCMd        R = DCMd::rotate_z(static_cast<double>(kPi) * 0.5); // 90° rotation about Z
        Vec3d       t{1.0, 0.0, 0.0};
        Transform4d T = Transform4d::from_rotation_translation(R, t);

        Vec3d p{1.0, 0.0, 0.0}; // Point at (1,0,0)
        Vec3d p_transformed = T.transform_point(p);

        // After 90° rotation + translation: (1,0,0) -> (0,1,0) + (1,0,0) = (1,1,0)
        CHECK(p_transformed[0] == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(p_transformed[1] == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(p_transformed[2] == doctest::Approx(0.0).epsilon(1e-6));
    }

    TEST_CASE("Vector transformation (no translation)") {
        DCMd        R = DCMd::rotate_z(static_cast<double>(kPi) * 0.5); // 90° rotation about Z
        Vec3d       t{1.0, 0.0, 0.0};
        Transform4d T = Transform4d::from_rotation_translation(R, t);

        Vec3d v{1.0, 0.0, 0.0}; // Vector along X
        Vec3d v_transformed = T.transform_vector(v);

        // Only rotation: (1,0,0) -> (0,1,0)
        CHECK(v_transformed[0] == doctest::Approx(0.0).epsilon(1e-6));
        CHECK(v_transformed[1] == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(v_transformed[2] == doctest::Approx(0.0).epsilon(1e-6));
    }

    TEST_CASE("Transform composition") {
        // First transform: 90° rotation about Z + translation (1,0,0)
        DCMd        R1 = DCMd::rotate_z(static_cast<double>(kPi) * 0.5);
        Vec3d       t1{1.0, 0.0, 0.0};
        Transform4d T1 = Transform4d::from_rotation_translation(R1, t1);

        // Second transform: translation (0,1,0)
        Transform4d T2 = Transform4d::from_rotation_translation(DCMd::identity(), Vec3d{0.0, 1.0, 0.0});

        // Composed transform: apply T1 then T2 (T2 * T1)
        Transform4d T_composed = T2 * T1;

        Vec3d p{0.0, 0.0, 0.0};
        Vec3d p_transformed = T_composed.transform_point(p);

        // Should be: rotate 90° then translate (1,0,0), then translate (0,1,0)
        // Result: (0,0,0) -> rotate -> (0,0,0) -> translate (1,0,0) -> (1,0,0) -> translate (0,1,0) -> (1,1,0)
        CHECK(p_transformed[0] == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(p_transformed[1] == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(p_transformed[2] == doctest::Approx(0.0).epsilon(1e-6));
    }

    TEST_CASE("Transform inversion") {
        DCMd        R = DCMd::rotate_z(static_cast<double>(kPi) * 0.5);
        Vec3d       t{1.0, 2.0, 3.0};
        Transform4d transform = Transform4d::from_rotation_translation(R, t);

        auto T_inv_opt = transform.inverse();
        REQUIRE(T_inv_opt.has_value());
        const Transform4d& T_inv = *T_inv_opt;

        // T * T_inv should be identity
        Transform4d T_identity = transform * T_inv;
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                double expected = (i == j) ? 1.0 : 0.0;
                CHECK(T_identity(i, j) == doctest::Approx(expected).epsilon(1e-6));
            }
        }
    }

    TEST_CASE("Quaternion and Euler integration") {
        // Create transform from quaternion + translation
        auto q_opt = Quatd::from_axis_angle(Vec3d{0.0, 0.0, 1.0}, static_cast<double>(kPi) * 0.5);
        REQUIRE(q_opt.has_value());

        const Quatd& q = *q_opt;
        Vec3d        t{1.0, 2.0, 3.0};
        Transform4d  T = Transform4d::from_quaternion_translation(q, t);

        // Extract back
        auto [q_back, t_back] = T.to_quaternion_translation();
        CHECK(q_back.w() == doctest::Approx(q.w()).epsilon(1e-6));
        CHECK(q_back.x() == doctest::Approx(q.x()).epsilon(1e-6));
        CHECK(q_back.y() == doctest::Approx(q.y()).epsilon(1e-6));
        CHECK(q_back.z() == doctest::Approx(q.z()).epsilon(1e-6));
        CHECK(t_back[0] == doctest::Approx(t[0]).epsilon(1e-6));
        CHECK(t_back[1] == doctest::Approx(t[1]).epsilon(1e-6));
        CHECK(t_back[2] == doctest::Approx(t[2]).epsilon(1e-6));

        // Test with Euler angles
        EulerZYXd   e{static_cast<double>(kPi) * 0.25, 0.0, 0.0}; // 45° yaw
        Transform4d T_euler = Transform4d::from_euler_translation(e, t);
        auto [e_back, t_e_back] = T_euler.to_euler_translation<EulerOrder::ZYX>();
        CHECK(e_back.yaw() == doctest::Approx(e.yaw()).epsilon(1e-6));
        CHECK(e_back.pitch() == doctest::Approx(e.pitch()).epsilon(1e-6));
        CHECK(e_back.roll() == doctest::Approx(e.roll()).epsilon(1e-6));
        CHECK(t_e_back[0] == doctest::Approx(t[0]).epsilon(1e-6));
        CHECK(t_e_back[1] == doctest::Approx(t[1]).epsilon(1e-6));
        CHECK(t_e_back[2] == doctest::Approx(t[2]).epsilon(1e-6));
    }
}
