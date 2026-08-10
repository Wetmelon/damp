// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstdlib>

#include "damp/backend.hpp"
#include "damp/estimation/sensor_fusion.hpp"
#include "damp/math/geometry.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

TEST_SUITE("Sensor Fusion Filters") {
    TEST_CASE("ComplementaryFilter basic functionality") {
        ComplementaryFilter<float> filter(0.98f);

        // Test initial orientation is identity
        auto q_init = filter.orientation();
        CHECK(q_init.w() == doctest::Approx(1.0));
        CHECK(q_init.x() == doctest::Approx(0.0));
        CHECK(q_init.y() == doctest::Approx(0.0));
        CHECK(q_init.z() == doctest::Approx(0.0));

        // Test update with stationary IMU (gravity aligned with -Z)
        Vec3<float> accel{0.0f, 0.0f, -9.81f}; // gravity down
        Vec3<float> gyro{0.0f, 0.0f, 0.0f};    // no rotation
        float       dt = 0.01f;

        filter.update(accel, gyro, dt);
        auto q = filter.orientation();

        // Should remain close to identity (normalized quaternion)
        CHECK(q.norm() == doctest::Approx(1.0).epsilon(1e-6));

        // Test with rotation around Z axis
        Vec3<float> accel_rotated{9.81f, 0.0f, 0.0f}; // gravity along +X (90 deg rotation)
        Vec3<float> gyro_z{0.0f, 0.0f, 1.57f};        // 90 deg/s around Z

        filter.update(accel_rotated, gyro_z, dt);
        q = filter.orientation();

        // Should have some rotation component
        CHECK(q.norm() == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(std::abs(q.w()) < 1.0f); // Not identity anymore
    }

    TEST_CASE("MadgwickFilter basic functionality") {
        MadgwickFilter<float> filter(0.1f);

        // Test initial orientation
        auto q_init = filter.orientation();
        CHECK(q_init.w() == doctest::Approx(1.0));
        CHECK(q_init.x() == doctest::Approx(0.0));
        CHECK(q_init.y() == doctest::Approx(0.0));
        CHECK(q_init.z() == doctest::Approx(0.0));

        // Test update with accelerometer and magnetometer
        Vec3<float> accel{0.0f, 0.0f, -9.81f}; // gravity down
        Vec3<float> gyro{0.0f, 0.0f, 0.0f};    // no rotation
        Vec3<float> mag{0.0f, 1.0f, 0.0f};     // magnetic north along +Y
        float       dt = 0.01f;

        filter.update(accel, gyro, mag, dt);
        auto q = filter.orientation();

        CHECK(q.norm() == doctest::Approx(1.0).epsilon(1e-6));
    }

    TEST_CASE("MahonyFilter basic functionality") {
        MahonyFilter<float> filter(0.5f, 0.0f);

        // Test initial orientation
        auto q_init = filter.orientation();
        CHECK(q_init.w() == doctest::Approx(1.0));
        CHECK(q_init.x() == doctest::Approx(0.0));
        CHECK(q_init.y() == doctest::Approx(0.0));
        CHECK(q_init.z() == doctest::Approx(0.0));

        // Test update with accelerometer and magnetometer
        Vec3<float> accel{0.0f, 0.0f, -9.81f}; // gravity down
        Vec3<float> gyro{0.0f, 0.0f, 0.0f};    // no rotation
        Vec3<float> mag{0.0f, 1.0f, 0.0f};     // magnetic north along +Y
        float       dt = 0.01f;

        filter.update(accel, gyro, mag, dt);
        auto q = filter.orientation();

        CHECK(q.norm() == doctest::Approx(1.0).epsilon(1e-6));
    }

    TEST_CASE("ESKFOrientationFilter basic functionality") {
        ESKFOrientationFilter<float, 6> filter;

        auto q_init = filter.orientation();
        CHECK(q_init.w() == doctest::Approx(1.0));
        CHECK(q_init.x() == doctest::Approx(0.0));
        CHECK(q_init.y() == doctest::Approx(0.0));
        CHECK(q_init.z() == doctest::Approx(0.0));

        auto bias_init = filter.gyro_bias();
        CHECK(bias_init[0] == doctest::Approx(0.0));
        CHECK(bias_init[1] == doctest::Approx(0.0));
        CHECK(bias_init[2] == doctest::Approx(0.0));

        // Specific force at rest (identity, ENU): a ≈ +g · ẑ
        Vec3<float> accel{0.0f, 0.0f, 9.81f};
        Vec3<float> gyro{0.0f, 0.0f, 0.0f};
        Vec3<float> mag{0.0f, 1.0f, 0.0f};
        float       dt = 0.01f;

        filter.update(accel, gyro, mag, dt);
        auto q = filter.orientation();
        auto bias = filter.gyro_bias();

        CHECK(q.norm() == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(std::abs(bias[0]) < 0.1f);
        CHECK(std::abs(bias[1]) < 0.1f);
        CHECK(std::abs(bias[2]) < 0.1f);
    }

    TEST_CASE("ESKFOrientationFilter IMU (NY=3) levels tilt without mag") {
        ESKFOrientationFilter<float, 3> filter;
        filter.set_accel_gate(0.2f);

        const float tilt = 0.12f;
        filter.set_orientation(
            Quaternion<float>::from_axis_angle(Vec3<float>{0.0f, 1.0f, 0.0f}, tilt).value()
        );
        const Vec3<float> g_nav{0.0f, 0.0f, -9.81f};
        const Vec3<float> accel{0.0f, 0.0f, 9.81f}; // level specific force
        const Vec3<float> gyro{};

        for (int i = 0; i < 400; ++i) {
            filter.update(accel, gyro, 0.01f, g_nav);
        }
        const auto err = (filter.orientation().conjugate() * Quaternion<float>::identity()).normalized().log().norm();
        CHECK(err < 0.05f);
    }

    TEST_CASE("ComplementaryFilter convergence test") {
        ComplementaryFilter<float> filter(0.98f);

        // Start with identity orientation
        Vec3<float> accel{0.0f, 0.0f, -9.81f}; // gravity down
        Vec3<float> gyro{0.0f, 0.0f, 0.0f};    // no rotation
        float       dt = 0.01f;

        // Run multiple updates to test convergence
        for (int i = 0; i < 100; ++i) {
            filter.update(accel, gyro, dt);
        }

        auto q = filter.orientation();
        CHECK(q.norm() == doctest::Approx(1.0).epsilon(1e-6));

        // For stationary case, should converge to orientation that aligns gravity with -Z
        // This means the quaternion should represent minimal rotation
        CHECK(std::abs(q.w()) > 0.9f); // Should be close to identity
    }

    TEST_CASE("Filter orientation consistency") {
        ComplementaryFilter<float>      comp_filter(0.95f);
        MadgwickFilter<float>           madgwick_filter(0.1f);
        MahonyFilter<float>             mahony_filter(0.5f, 0.0f);
        ESKFOrientationFilter<float, 6> eskf_filter;

        // Complementary/Madgwick/Mahony historically take gravity-direction accel (−z).
        // ESKF uses INS specific force (+z at rest in ENU).
        Vec3<float> accel_gravity{0.0f, 0.0f, -9.81f};
        Vec3<float> accel_sf{0.0f, 0.0f, 9.81f};
        Vec3<float> gyro{0.1f, 0.05f, -0.08f};
        Vec3<float> mag{0.0f, 1.0f, 0.0f};
        float       dt = 0.01f;

        comp_filter.update(accel_gravity, gyro, dt);
        madgwick_filter.update(accel_gravity, gyro, mag, dt);
        mahony_filter.update(accel_gravity, gyro, mag, dt);
        eskf_filter.update(accel_sf, gyro, mag, dt);

        CHECK(comp_filter.orientation().norm() == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(madgwick_filter.orientation().norm() == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(mahony_filter.orientation().norm() == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(eskf_filter.orientation().norm() == doctest::Approx(1.0).epsilon(1e-6));
    }

    TEST_CASE("ESKF basic operation") {
        ESKFOrientationFilter<float, 6> filter;

        Vec3<float> accel{0.0f, 0.0f, 9.81f};
        Vec3<float> gyro{0.0f, 0.0f, 0.0f};
        Vec3<float> mag{0.0f, 1.0f, 0.0f};
        float       dt = 0.01f;

        for (int i = 0; i < 10; ++i) {
            filter.update(accel, gyro, mag, dt);
        }

        auto q = filter.orientation();
        auto bias = filter.gyro_bias();

        // Check that orientation remains valid
        CHECK(q.norm() == doctest::Approx(1.0).epsilon(1e-6));
        // For stationary case, should stay close to identity
        CHECK(std::abs(q.w()) > 0.99f);
        // Bias should remain small
        CHECK(std::abs(bias[0]) < 0.01f);
        CHECK(std::abs(bias[1]) < 0.01f);
        CHECK(std::abs(bias[2]) < 0.01f);
    }

    TEST_CASE("Filter robustness to noisy measurements") {
        ComplementaryFilter<float> filter(0.9f);

        Vec3<float> true_accel{0.0f, 0.0f, -9.81f};
        Vec3<float> true_gyro{0.0f, 0.0f, 0.0f};
        float       dt = 0.01f;

        // Add noise to measurements
        for (int i = 0; i < 100; ++i) {
            const auto unit = []() {
                return static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
            };
            Vec3<float> noisy_accel = true_accel + Vec3<float>{
                                          0.1f * (unit() - 0.5f),
                                          0.1f * (unit() - 0.5f),
                                          0.1f * (unit() - 0.5f),
                                      };

            Vec3<float> noisy_gyro = true_gyro + Vec3<float>{
                                         0.01f * (unit() - 0.5f),
                                         0.01f * (unit() - 0.5f),
                                         0.01f * (unit() - 0.5f),
                                     };

            filter.update(noisy_accel, noisy_gyro, dt);

            // Check quaternion remains valid
            auto q = filter.orientation();
            CHECK(q.norm() == doctest::Approx(1.0).epsilon(1e-5));
        }
    }

    // --- Convergence / observability tests for the orientation filters -------
    // The cases above only check unit-norm; these check the filters actually
    // recover a known tilt (gravity) and heading (magnetometer).

    TEST_CASE("Complementary filter recovers a static tilt") {
        // Body rotated 25 deg about X relative to the world.
        auto        q_true = Quaternion<float>::from_axis_angle(Vec3<float>{1.0f, 0.0f, 0.0f}, 0.4363f).value();
        Vec3<float> accel = q_true.rotate(Vec3<float>{0.0f, 0.0f, -9.81f});

        ComplementaryFilter<float> filter(0.9f);
        for (int i = 0; i < 4000; ++i) {
            filter.update(accel, Vec3<float>{0.0f, 0.0f, 0.0f}, 0.005f);
        }

        // Predicted up axis must align with the true up axis.
        Vec3<float> up_est = filter.orientation().rotate(Vec3<float>{0.0f, 0.0f, 1.0f});
        Vec3<float> up_true = q_true.rotate(Vec3<float>{0.0f, 0.0f, 1.0f});
        CHECK((up_est - up_true).norm() < 3e-2f);
    }

    TEST_CASE("Mahony filter recovers tilt and heading from accel+mag") {
        // 25 deg about X then 40 deg about Z (yaw): heading needs the magnetometer.
        auto              qx = Quaternion<float>::from_axis_angle(Vec3<float>{1.0f, 0.0f, 0.0f}, 0.4363f).value();
        auto              qz = Quaternion<float>::from_axis_angle(Vec3<float>{0.0f, 0.0f, 1.0f}, 0.6981f).value();
        Quaternion<float> q_true = (qz * qx).normalized();

        Vec3<float> accel = q_true.rotate(Vec3<float>{0.0f, 0.0f, -9.81f});
        Vec3<float> mag = q_true.rotate(Vec3<float>{0.0f, 1.0f, 0.0f}); // north = +Y

        MahonyFilter<float> filter(8.0f, 0.0f); // strong proportional gain for the test
        for (int i = 0; i < 6000; ++i) {
            filter.update(accel, Vec3<float>{0.0f, 0.0f, 0.0f}, mag, 0.005f);
        }

        Vec3<float> up_est = filter.orientation().rotate(Vec3<float>{0.0f, 0.0f, 1.0f});
        Vec3<float> north_est = filter.orientation().rotate(Vec3<float>{0.0f, 1.0f, 0.0f});
        CHECK((up_est - q_true.rotate(Vec3<float>{0.0f, 0.0f, 1.0f})).norm() < 5e-2f);
        CHECK((north_est - q_true.rotate(Vec3<float>{0.0f, 1.0f, 0.0f})).norm() < 5e-2f);
    }

    TEST_CASE("Madgwick filter recovers tilt and heading from accel+mag") {
        auto              qx = Quaternion<float>::from_axis_angle(Vec3<float>{1.0f, 0.0f, 0.0f}, 0.3491f).value(); // 20 deg
        auto              qz = Quaternion<float>::from_axis_angle(Vec3<float>{0.0f, 0.0f, 1.0f}, 0.5236f).value(); // 30 deg
        Quaternion<float> q_true = (qz * qx).normalized();

        Vec3<float> accel = q_true.rotate(Vec3<float>{0.0f, 0.0f, -9.81f});
        Vec3<float> mag = q_true.rotate(Vec3<float>{0.0f, 1.0f, 0.0f});

        MadgwickFilter<float> filter(5.0f); // large beta for fast test convergence
        for (int i = 0; i < 8000; ++i) {
            filter.update(accel, Vec3<float>{0.0f, 0.0f, 0.0f}, mag, 0.005f);
        }

        Vec3<float> up_est = filter.orientation().rotate(Vec3<float>{0.0f, 0.0f, 1.0f});
        Vec3<float> north_est = filter.orientation().rotate(Vec3<float>{0.0f, 1.0f, 0.0f});
        CHECK((up_est - q_true.rotate(Vec3<float>{0.0f, 0.0f, 1.0f})).norm() < 5e-2f);
        CHECK((north_est - q_true.rotate(Vec3<float>{0.0f, 1.0f, 0.0f})).norm() < 5e-2f);
    }

    TEST_CASE("ESKF orientation filter recovers tilt and heading") {
        // Mag + accel must observe yaw. Sensors are body-frame specific force / field.
        auto              qx = Quaternion<float>::from_axis_angle(Vec3<float>{1.0f, 0.0f, 0.0f}, 0.3491f).value();
        auto              qz = Quaternion<float>::from_axis_angle(Vec3<float>{0.0f, 0.0f, 1.0f}, 0.5236f).value();
        Quaternion<float> q_true = (qz * qx).normalized();

        const Vec3<float> g_nav{0.0f, 0.0f, -9.81f};
        const Vec3<float> m_nav{0.0f, 1.0f, 0.0f};
        const Vec3<float> accel = q_true.conjugate().rotate(Vec3<float>{-g_nav[0], -g_nav[1], -g_nav[2]});
        const Vec3<float> mag = q_true.conjugate().rotate(m_nav);

        ESKFOrientationFilter<float, 6> filter;
        for (int i = 0; i < 20000; ++i) {
            filter.update(accel, Vec3<float>{0.0f, 0.0f, 0.0f}, mag, 0.01f, g_nav, m_nav);
        }

        Vec3<float> up_est = filter.orientation().rotate(Vec3<float>{0.0f, 0.0f, 1.0f});
        Vec3<float> north_est = filter.orientation().rotate(Vec3<float>{0.0f, 1.0f, 0.0f});
        CHECK((up_est - q_true.rotate(Vec3<float>{0.0f, 0.0f, 1.0f})).norm() < 5e-2f);
        CHECK((north_est - q_true.rotate(Vec3<float>{0.0f, 1.0f, 0.0f})).norm() < 8e-2f);
    }

} // TEST_SUITE