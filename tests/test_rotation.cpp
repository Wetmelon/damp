// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>

#include "damp/math/geometry.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {
constexpr double kEps = std::numeric_limits<double>::epsilon();

template<size_t R, size_t C>
double max_abs_diff(const Matrix<R, C, double>& a, const Matrix<R, C, double>& b) {
    double m = 0.0;
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            m = std::fmax(m, std::fabs(a(i, j) - b(i, j)));
        }
    }
    return m;
}

template<EulerOrder Order>
void check_euler_round_trip(double a1, double a2, double a3) {
    const Euler<double, Order> e{a1, a2, a3};
    const DCM<double>          R = e.to_dcm();
    const DCM<double>          R_back = Euler<double, Order>::from_dcm(R).to_dcm();
    CHECK(max_abs_diff<3, 3>(R, R_back) <= 64 * kEps);
}

template<EulerOrder Order>
void check_euler_order() {
    constexpr double half_pi = std::numbers::pi / 2.0;
    // Generic, near gimbal lock (cos θ ≈ 1e-3 and 1e-6), and exactly at it.
    for (const double a2 : {0.4, -1.1, half_pi - 1e-3, -(half_pi - 1e-3), half_pi - 1e-6, half_pi, -half_pi}) {
        check_euler_round_trip<Order>(0.3, a2, 0.7);
        check_euler_round_trip<Order>(-2.5, a2, -0.2);
    }
}

// Round trip of DCMs built through a quaternion, which carry rounding of size ε
// in every entry, including the cos θ₂-sized ones that set θ₃ near lock.
template<typename T, EulerOrder Order>
T max_euler_error_near_lock() {
    const T half_pi = std::numbers::pi_v<T> / T{2};
    T       worst = T{0};
    for (const double c : {1e-2, 1e-3, 3e-4, 1e-4, 1e-5, 1e-6, 1e-8, 1e-9, 1e-12, 0.0}) {
        for (const T sign : {T{1}, T{-1}}) {
            for (const T a1 : {T{0.75}, T{-2.5}}) {
                const Euler<T, Order> e{a1, sign * (half_pi - static_cast<T>(c)), T{-0.375}};
                const DCM<T>          R = e.to_quaternion().to_dcm();
                const DCM<T>          R_back = Euler<T, Order>::from_dcm(R).to_dcm();
                for (size_t r = 0; r < 3; ++r) {
                    for (size_t k = 0; k < 3; ++k) {
                        worst = std::fmax(worst, std::fabs(R(r, k) - R_back(r, k)));
                    }
                }
            }
        }
    }
    return worst;
}

template<EulerOrder Order>
void check_euler_near_lock() {
    CHECK(max_euler_error_near_lock<float, Order>() <= 16.0f * std::numeric_limits<float>::epsilon());
    CHECK(max_euler_error_near_lock<double, Order>() <= 16.0 * kEps);
}

} // namespace

constexpr float kPi = std::numbers::pi_v<float>;

TEST_SUITE("DCM") {
    TEST_CASE("DCM identity") {
        DCM<float> R = DCM<float>::identity();
        CHECK(R(0, 0) == 1.0f);
        CHECK(R(1, 1) == 1.0f);
        CHECK(R(2, 2) == 1.0f);
        CHECK(R(0, 1) == 0.0f);
        CHECK(R(1, 0) == 0.0f);
    }

    TEST_CASE("DCM basic rotations") {
        // 90 deg about Z
        auto        Rz = DCM<float>::rotate_z(kPi / 2);
        Vec3<float> v{1.0f, 0.0f, 0.0f};
        auto        v_rot = Rz * v;
        CHECK(v_rot[0] == doctest::Approx(0.0f).epsilon(1e-6));
        CHECK(v_rot[1] == doctest::Approx(1.0f).epsilon(1e-6));
        CHECK(v_rot[2] == doctest::Approx(0.0f).epsilon(1e-6));

        // 90 deg about X
        auto        Rx = DCM<float>::rotate_x(kPi / 2);
        Vec3<float> vy{0.0f, 1.0f, 0.0f};
        auto        vy_rot = Rx * vy;
        CHECK(vy_rot[0] == doctest::Approx(0.0f).epsilon(1e-6));
        CHECK(vy_rot[1] == doctest::Approx(0.0f).epsilon(1e-6));
        CHECK(vy_rot[2] == doctest::Approx(1.0f).epsilon(1e-6));

        // 90 deg about Y
        auto        Ry = DCM<float>::rotate_y(kPi / 2);
        Vec3<float> vz{0.0f, 0.0f, 1.0f};
        auto        vz_rot = Ry * vz;
        CHECK(vz_rot[0] == doctest::Approx(1.0f).epsilon(1e-6));
        CHECK(vz_rot[1] == doctest::Approx(0.0f).epsilon(1e-6));
        CHECK(vz_rot[2] == doctest::Approx(0.0f).epsilon(1e-6));
    }

    TEST_CASE("DCM composition") {
        auto R1 = DCM<float>::rotate_z(kPi / 4);
        auto R2 = DCM<float>::rotate_z(kPi / 4);
        auto R_combined = R1 * R2; // Should be 90 deg

        Vec3<float> v{1.0f, 0.0f, 0.0f};
        auto        v_rot = R_combined * v;
        CHECK(v_rot[0] == doctest::Approx(0.0f).epsilon(1e-5));
        CHECK(v_rot[1] == doctest::Approx(1.0f).epsilon(1e-5));
    }

    TEST_CASE("DCM transpose is inverse") {
        auto R = DCM<float>::rotate_z(0.5f) * DCM<float>::rotate_x(0.3f);
        auto R_inv = R.transpose();
        auto I = R * R_inv;

        CHECK(I(0, 0) == doctest::Approx(1.0f).epsilon(1e-6));
        CHECK(I(1, 1) == doctest::Approx(1.0f).epsilon(1e-6));
        CHECK(I(2, 2) == doctest::Approx(1.0f).epsilon(1e-6));
        CHECK(I(0, 1) == doctest::Approx(0.0f).epsilon(1e-6));
    }

    TEST_CASE("DCM from axis-angle") {
        auto R_opt = DCM<float>::from_axis_angle(Vec3<float>{0.0f, 0.0f, 1.0f}, kPi / 2);
        REQUIRE(R_opt.has_value());
        auto R = R_opt.value();

        Vec3<float> v{1.0f, 0.0f, 0.0f};
        auto        v_rot = R * v;
        CHECK(v_rot[0] == doctest::Approx(0.0f).epsilon(1e-5));
        CHECK(v_rot[1] == doctest::Approx(1.0f).epsilon(1e-5));
    }

    // Regression: eps is a *linear* axis-length tolerance, so a short (but
    // nonzero) axis — e.g. the small rotation vectors iterative solvers feed in —
    // must still be accepted, not rejected by a squared-norm-vs-eps mismatch.
    TEST_CASE("DCM/Quaternion from a short axis is accepted (eps is linear)") {
        const Vec3<double> tiny{1e-5, 0.0, 0.0}; // |axis|² = 1e-10 < default eps 1e-9
        const double       angle = 0.02;

        auto R_opt = DCM<double>::from_axis_angle(tiny, angle);
        REQUIRE(R_opt.has_value());
        // Direction (axis normalized internally) governs the rotation, not |axis|.
        auto R = R_opt.value();
        auto ref = DCM<double>::rotate_x(angle);
        CHECK(R(1, 1) == doctest::Approx(ref(1, 1)).epsilon(1e-9));
        CHECK(R(2, 1) == doctest::Approx(ref(2, 1)).epsilon(1e-9));

        auto q_opt = Quaternion<double>::from_axis_angle(tiny, angle);
        REQUIRE(q_opt.has_value());
        CHECK(q_opt.value().w() == doctest::Approx(std::cos(angle / 2.0)));

        // A genuinely (near-)zero axis is still rejected.
        CHECK_FALSE(DCM<double>::from_axis_angle(Vec3<double>{0.0, 0.0, 0.0}, angle).has_value());
        CHECK_FALSE(Quaternion<double>::from_axis_angle(Vec3<double>{0.0, 0.0, 0.0}, angle).has_value());
    }
}

TEST_SUITE("Euler") {
    TEST_CASE("EulerZYX construction and accessors") {
        EulerZYX<float> e(0.1f, 0.2f, 0.3f);
        CHECK(e.yaw() == 0.1f);
        CHECK(e.pitch() == 0.2f);
        CHECK(e.roll() == 0.3f);
        CHECK(e.angle1 == 0.1f);
        CHECK(e.angle2 == 0.2f);
        CHECK(e.angle3 == 0.3f);
    }

    TEST_CASE("EulerXYZ construction and angles") {
        EulerXYZ<float> e(0.1f, 0.2f, 0.3f);
        CHECK(e.angle1 == 0.1f);
        CHECK(e.angle2 == 0.2f);
        CHECK(e.angle3 == 0.3f);
    }

    TEST_CASE("EulerZYX to DCM and back") {
        EulerZYX<float> e_orig(0.3f, -0.2f, 0.5f); // yaw, pitch, roll
        auto            R = e_orig.to_dcm();
        auto            e_back = EulerZYX<float>::from_dcm(R);

        CHECK(e_back.yaw() == doctest::Approx(e_orig.yaw()).epsilon(1e-5));
        CHECK(e_back.pitch() == doctest::Approx(e_orig.pitch()).epsilon(1e-5));
        CHECK(e_back.roll() == doctest::Approx(e_orig.roll()).epsilon(1e-5));
    }

    TEST_CASE("EulerXYZ to DCM and back") {
        EulerXYZ<float> e_orig(0.2f, -0.1f, 0.4f); // roll, pitch, yaw
        auto            R = e_orig.to_dcm();
        auto            e_back = EulerXYZ<float>::from_dcm(R);

        CHECK(e_back.angle1 == doctest::Approx(e_orig.angle1).epsilon(1e-5));
        CHECK(e_back.angle2 == doctest::Approx(e_orig.angle2).epsilon(1e-5));
        CHECK(e_back.angle3 == doctest::Approx(e_orig.angle3).epsilon(1e-5));
    }

    TEST_CASE("Euler order produces correct rotation") {
        // For ZYX: rotate Z first, then Y, then X
        float yaw = kPi / 4;
        float pitch = 0.0f;
        float roll = 0.0f;

        EulerZYX<float> e(yaw, pitch, roll);
        auto            R = e.to_dcm();

        // Should be same as rotate_z(yaw)
        auto R_expected = DCM<float>::rotate_z(yaw);
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                CHECK(R(i, j) == doctest::Approx(R_expected(i, j)).epsilon(1e-6));
            }
        }
    }

    TEST_CASE("Euler extraction round-trips near and at gimbal lock, every order") {
        check_euler_order<EulerOrder::ZYX>();
        check_euler_order<EulerOrder::XYZ>();
        check_euler_order<EulerOrder::ZXY>();
        check_euler_order<EulerOrder::YXZ>();
        check_euler_order<EulerOrder::YZX>();
        check_euler_order<EulerOrder::XZY>();
    }

    TEST_CASE("Euler extraction reproduces quaternion-built DCMs to rounding near gimbal lock, float and double") {
        check_euler_near_lock<EulerOrder::ZYX>();
        check_euler_near_lock<EulerOrder::XYZ>();
        check_euler_near_lock<EulerOrder::ZXY>();
        check_euler_near_lock<EulerOrder::YXZ>();
        check_euler_near_lock<EulerOrder::YZX>();
        check_euler_near_lock<EulerOrder::XZY>();
    }

    TEST_CASE("Euler extraction sets theta3 to 0 at exact gimbal lock") {
        // R = Rz(0.5) Ry(π/2) built exactly: R_ij = R_ii = 0 for ZYX.
        const DCM<double> R{{0.0, -std::sin(0.5), std::cos(0.5)}, {0.0, std::cos(0.5), std::sin(0.5)}, {-1.0, 0.0, 0.0}};
        const auto        e = EulerZYX<double>::from_dcm(R);
        CHECK(e.roll() == 0.0);
        CHECK(e.yaw() == doctest::Approx(0.5));
        CHECK(max_abs_diff<3, 3>(R, e.to_dcm()) <= 4 * kEps);
    }

    TEST_CASE("EulerXYZ is intrinsic X-Y'-Z'' and URDF rpy is EulerZYX") {
        const double r = 0.3;
        const double p = -0.2;
        const double y = 0.5;
        const auto   Rx = DCM<double>::rotate_x(r);
        const auto   Ry = DCM<double>::rotate_y(p);
        const auto   Rz = DCM<double>::rotate_z(y);
        CHECK(max_abs_diff<3, 3>(EulerXYZ<double>{r, p, y}.to_dcm(), Rx * Ry * Rz) <= 4 * kEps);
        CHECK(max_abs_diff<3, 3>(EulerZYX<double>{y, p, r}.to_dcm(), Rz * Ry * Rx) <= 4 * kEps);
    }
}

TEST_SUITE("Quaternion (rotation.hpp)") {
    TEST_CASE("Quaternion identity") {
        Quaternion<float> q;
        CHECK(q.w() == 1.0f);
        CHECK(q.x() == 0.0f);
        CHECK(q.y() == 0.0f);
        CHECK(q.z() == 0.0f);

        auto q_id = Quaternion<float>::identity();
        CHECK(q_id.w() == 1.0f);
    }

    TEST_CASE("Quaternion from axis-angle") {
        auto q_opt = Quaternion<float>::from_axis_angle(Vec3<float>{0.0f, 0.0f, 1.0f}, kPi / 2);
        REQUIRE(q_opt.has_value());
        auto q = q_opt.value();

        Vec3<float> v{1.0f, 0.0f, 0.0f};
        auto        v_rot = q.rotate(v);
        CHECK(v_rot[0] == doctest::Approx(0.0f).epsilon(1e-5));
        CHECK(v_rot[1] == doctest::Approx(1.0f).epsilon(1e-5));
        CHECK(v_rot[2] == doctest::Approx(0.0f).epsilon(1e-5));
    }

    TEST_CASE("Quaternion slerp") {
        auto qa = Quaternion<float>::identity();
        auto qb_opt = Quaternion<float>::from_axis_angle(Vec3<float>{0.0f, 0.0f, 1.0f}, kPi / 2);
        REQUIRE(qb_opt.has_value());
        auto qb = qb_opt.value();

        auto        qm = Quaternion<float>::slerp(qa, qb, 0.5f);
        Vec3<float> v{1.0f, 0.0f, 0.0f};
        auto        v_rot = qm.rotate(v);

        // Should be 45 deg rotation
        CHECK(v_rot[0] == doctest::Approx(std::cos(kPi / 4)).epsilon(1e-5));
        CHECK(v_rot[1] == doctest::Approx(std::sin(kPi / 4)).epsilon(1e-5));
    }

    TEST_CASE("Quaternion slerp keeps double precision for small separations") {
        const Quatd a = Quatd::identity();
        for (const double th : {1e-6, 1e-4, 2e-3, 2.8e-3}) {
            const Quatd  b = Quatd::from_axis_angle(Vec3<double>{0.0, 0.0, 1.0}, th).value();
            const Quatd  q = Quatd::slerp(a, b, 0.3);
            const double got = 2.0 * std::atan2(std::hypot(q.x(), q.y(), q.z()), q.w());
            CHECK(std::fabs(got - (0.3 * th)) <= 8 * kEps * 0.3 * th);
        }
    }
}

TEST_SUITE("Cross-conversions") {
    TEST_CASE("Quaternion <-> DCM round-trip") {
        auto q_orig_opt = Quaternion<float>::from_axis_angle(Vec3<float>{1.0f, 2.0f, 3.0f}, 0.7f);
        REQUIRE(q_orig_opt.has_value());
        auto q_orig = q_orig_opt.value();

        auto R = q_orig.to_dcm();
        auto q_back = R.to_quaternion();

        // Quaternions may differ by sign, check rotation is same
        Vec3<float> v{1.0f, 2.0f, 3.0f};
        auto        v1 = q_orig.rotate(v);
        auto        v2 = q_back.rotate(v);

        CHECK(v1[0] == doctest::Approx(v2[0]).epsilon(1e-5));
        CHECK(v1[1] == doctest::Approx(v2[1]).epsilon(1e-5));
        CHECK(v1[2] == doctest::Approx(v2[2]).epsilon(1e-5));
    }

    TEST_CASE("Quaternion <-> EulerZYX round-trip") {
        EulerZYX<float> e_orig(0.3f, -0.2f, 0.5f);
        auto            q = e_orig.to_quaternion();
        auto            e_back = q.to_euler<EulerOrder::ZYX>();

        CHECK(e_back.angle1 == doctest::Approx(e_orig.angle1).epsilon(1e-5));
        CHECK(e_back.angle2 == doctest::Approx(e_orig.angle2).epsilon(1e-5));
        CHECK(e_back.angle3 == doctest::Approx(e_orig.angle3).epsilon(1e-5));
    }

    TEST_CASE("DCM <-> EulerZYX round-trip") {
        EulerZYX<float> e_orig(0.4f, 0.1f, -0.3f);
        auto            R = DCM<float>::from_euler(e_orig);
        auto            e_back = R.to_euler<EulerOrder::ZYX>();

        CHECK(e_back.angle1 == doctest::Approx(e_orig.angle1).epsilon(1e-5));
        CHECK(e_back.angle2 == doctest::Approx(e_orig.angle2).epsilon(1e-5));
        CHECK(e_back.angle3 == doctest::Approx(e_orig.angle3).epsilon(1e-5));
    }

    TEST_CASE("All representations produce same rotation") {
        EulerZYX<float> e(0.25f, -0.15f, 0.35f);
        auto            R = e.to_dcm();
        auto            q = e.to_quaternion();

        Vec3<float> v{1.0f, -2.0f, 0.5f};

        auto v_dcm = R * v;
        auto v_quat = q.rotate(v);

        CHECK(v_dcm[0] == doctest::Approx(v_quat[0]).epsilon(1e-5));
        CHECK(v_dcm[1] == doctest::Approx(v_quat[1]).epsilon(1e-5));
        CHECK(v_dcm[2] == doctest::Approx(v_quat[2]).epsilon(1e-5));
    }

    // Python-generated golden test values
    // import numpy as np
    // from scipy.spatial.transform import Rotation
    // # EulerZYX(yaw, pitch, roll) = intrinsic 'ZYX' [yaw, pitch, roll]
    // #   = extrinsic (lowercase) 'xyz' [roll, pitch, yaw] = URDF / ROS rpy
    // # So EulerZYX(0.3, -0.2, 0.5) means yaw=0.3, pitch=-0.2, roll=0.5
    // r = Rotation.from_euler('xyz', [0.5, -0.2, 0.3])  # roll, pitch, yaw -> extrinsic xyz
    // print("quat (w,x,y,z):", r.as_quat()[[3,0,1,2]])
    // print("dcm:", r.as_matrix())
    TEST_CASE("Golden test: EulerZYX(0.3, -0.2, 0.5) to Quaternion") {
        EulerZYX<double> e(0.3, -0.2, 0.5); // yaw, pitch, roll
        auto             q = e.to_quaternion();

        // scipy from_euler('xyz', [0.5, -0.2, 0.3]): w=0.9496, x=0.2579, y=-0.0589, z=0.1685
        CHECK(q.w() == doctest::Approx(0.949555).epsilon(1e-5));
        CHECK(q.x() == doctest::Approx(0.257859).epsilon(1e-5));
        CHECK(q.y() == doctest::Approx(-0.0588568).epsilon(1e-5));
        CHECK(q.z() == doctest::Approx(0.168491).epsilon(1e-5));
    }

    TEST_CASE("Golden test: EulerZYX(0.3, -0.2, 0.5) to DCM") {
        EulerZYX<double> e(0.3, -0.2, 0.5); // yaw, pitch, roll
        auto             R = e.to_dcm();

        // scipy from_euler('xyz', [0.5, -0.2, 0.3]).as_matrix():
        // [[ 0.93629336, -0.35033617, -0.02488181]
        //  [ 0.28962948,  0.81023922, -0.50953595]
        //  [ 0.19866933,  0.46986869,  0.86008934]]
        CHECK(R(0, 0) == doctest::Approx(0.93629336).epsilon(1e-5));
        CHECK(R(0, 1) == doctest::Approx(-0.35033617).epsilon(1e-5));
        CHECK(R(0, 2) == doctest::Approx(-0.02488181).epsilon(1e-5));
        CHECK(R(1, 0) == doctest::Approx(0.28962948).epsilon(1e-5));
        CHECK(R(1, 1) == doctest::Approx(0.81023922).epsilon(1e-5));
        CHECK(R(1, 2) == doctest::Approx(-0.50953595).epsilon(1e-5));
        CHECK(R(2, 0) == doctest::Approx(0.19866933).epsilon(1e-5));
        CHECK(R(2, 1) == doctest::Approx(0.46986869).epsilon(1e-5));
        CHECK(R(2, 2) == doctest::Approx(0.86008934).epsilon(1e-5));
    }

    // Python golden: rotate vector [1, 2, 3] by ZYX(yaw=0.3, pitch=-0.2, roll=0.5)
    // from scipy.spatial.transform import Rotation
    // r = Rotation.from_euler('xyz', [0.5, -0.2, 0.3])  # roll, pitch, yaw
    // r.apply([1,2,3]) yields the expected rotation
    TEST_CASE("Golden test: rotate vector [1,2,3] by EulerZYX(0.3,-0.2,0.5)") {
        EulerZYX<double> e(0.3, -0.2, 0.5); // yaw, pitch, roll
        auto             q = e.to_quaternion();
        auto             R = e.to_dcm();

        Vec3<double> v{1.0, 2.0, 3.0};
        auto         v_q = q.rotate(v);
        auto         v_R = R * v;

        // Verify quat and DCM produce the same rotation
        CHECK(v_q[0] == doctest::Approx(v_R[0]).epsilon(1e-6));
        CHECK(v_q[1] == doctest::Approx(v_R[1]).epsilon(1e-6));
        CHECK(v_q[2] == doctest::Approx(v_R[2]).epsilon(1e-6));

        // Values from running the implementation
        CHECK(v_R[0] == doctest::Approx(0.161307).epsilon(1e-3));
        CHECK(v_R[1] == doctest::Approx(0.381499).epsilon(1e-3));
        CHECK(v_R[2] == doctest::Approx(3.71868).epsilon(1e-3));
    }
}

TEST_SUITE("Type aliases") {
    TEST_CASE("Convenience aliases compile") {
        Quatf     qf;
        Quatd     qd;
        DCMf      Rf;
        DCMd      Rd;
        EulerZYXf ef;
        EulerZYXd ed;
        EulerXYZf exf;
        EulerXYZd exd;

        CHECK(qf.w() == 1.0f);
        CHECK(qd.w() == 1.0);
        CHECK(Rf(0, 0) == 1.0f);
        CHECK(Rd(0, 0) == 1.0);
    }
}
