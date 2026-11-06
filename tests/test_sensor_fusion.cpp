// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstdlib>
#include <limits>

#include "damp/backend.hpp"
#include "damp/estimation/eskf.hpp"
#include "damp/estimation/ins_mechanization.hpp"
#include "damp/estimation/sensor_fusion.hpp"
#include "damp/math/geometry.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

// A reset that would install an identity covariance in place of the design's P₀.
template<typename F>
concept BareReset = requires(F& f) { f.reset(); };
template<typename F, typename X>
concept StateOnlyReset = requires(F& f, const X& x) { f.reset(x); };

[[nodiscard]] double rotation_angle_between(const Quaternion<double>& a, const Quaternion<double>& b) {
    const double w = damp::abs((a.conjugate() * b).w());
    return 2.0 * std::acos(damp::min(1.0, w));
}

} // namespace

TEST_SUITE("Sensor Fusion Filters") {
    TEST_CASE("ComplementaryFilter basic functionality") {
        ComplementaryFilter<float> filter(0.49f);

        // Test initial orientation is identity
        auto q_init = filter.orientation();
        CHECK(q_init.w() == doctest::Approx(1.0));
        CHECK(q_init.x() == doctest::Approx(0.0));
        CHECK(q_init.y() == doctest::Approx(0.0));
        CHECK(q_init.z() == doctest::Approx(0.0));

        // Test update with stationary IMU (specific force +g along Z)
        Vec3<float> accel{0.0f, 0.0f, 9.81f}; // level at rest
        Vec3<float> gyro{0.0f, 0.0f, 0.0f};   // no rotation
        float       dt = 0.01f;

        filter.step(accel, gyro, dt);
        auto q = filter.orientation();

        // Should remain close to identity (normalized quaternion)
        CHECK(q.norm() == doctest::Approx(1.0).epsilon(1e-6));

        // Test with rotation around Z axis
        Vec3<float> accel_rotated{-9.81f, 0.0f, 0.0f}; // gravity along +X (90 deg rotation)
        Vec3<float> gyro_z{0.0f, 0.0f, 1.57f};         // 90 deg/s around Z

        filter.step(accel_rotated, gyro_z, dt);
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
        Vec3<float> accel{0.0f, 0.0f, 9.81f}; // level at rest (specific force)
        Vec3<float> gyro{0.0f, 0.0f, 0.0f};   // no rotation
        Vec3<float> mag{0.0f, 1.0f, 0.0f};    // magnetic north along +Y
        float       dt = 0.01f;

        filter.step(accel, gyro, mag, dt);
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
        Vec3<float> accel{0.0f, 0.0f, 9.81f}; // level at rest (specific force)
        Vec3<float> gyro{0.0f, 0.0f, 0.0f};   // no rotation
        Vec3<float> mag{0.0f, 1.0f, 0.0f};    // magnetic north along +Y
        float       dt = 0.01f;

        filter.step(accel, gyro, mag, dt);
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

        filter.step(accel, gyro, mag, dt);
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
            filter.step(accel, gyro, 0.01f, g_nav);
        }
        const auto err = (filter.orientation().conjugate() * Quaternion<float>::identity()).normalized().log().norm();
        CHECK(err < 0.05f);
    }

    TEST_CASE("ComplementaryFilter crossover τ does not depend on the sample rate") {
        //! τ = 0.5 s: after 1 s the tilt error is about e⁻² of its start at both 100 Hz and 1 kHz.
        const Vec3<double> accel = Quaternion<double>::from_axis_angle(Vec3<double>{1.0, 0.0, 0.0}, 0.5)
                                       .value()
                                       .conjugate()
                                       .rotate(Vec3<double>{0.0, 0.0, 9.81});
        double             residual[2] = {0.0, 0.0};
        int                i = 0;
        for (const double dt : {0.01, 0.001}) {
            ComplementaryFilter<double> comp(0.5);
            const int                   n = static_cast<int>(1.0 / dt + 0.5);
            for (int k = 0; k < n; ++k) {
                comp.step(accel, Vec3<double>{0.0, 0.0, 0.0}, dt);
            }
            residual[i++] = 0.5 - Euler<double, EulerOrder::ZYX>::from_quaternion(comp.orientation()).roll();
        }
        CHECK(residual[0] == doctest::Approx(0.5 * std::exp(-2.0)).epsilon(0.03));
        CHECK(residual[1] == doctest::Approx(0.5 * std::exp(-2.0)).epsilon(0.03));
    }

    TEST_CASE("ComplementaryFilter convergence test") {
        ComplementaryFilter<float> filter(0.49f);

        // Start with identity orientation
        Vec3<float> accel{0.0f, 0.0f, 9.81f}; // level at rest (specific force)
        Vec3<float> gyro{0.0f, 0.0f, 0.0f};   // no rotation
        float       dt = 0.01f;

        // Run multiple updates to test convergence
        for (int i = 0; i < 100; ++i) {
            filter.step(accel, gyro, dt);
        }

        auto q = filter.orientation();
        CHECK(q.norm() == doctest::Approx(1.0).epsilon(1e-6));

        // For stationary case, should converge to orientation that aligns up with +Z
        // This means the quaternion should represent minimal rotation
        CHECK(std::abs(q.w()) > 0.9f); // Should be close to identity
    }

    TEST_CASE("Filter orientation consistency") {
        ComplementaryFilter<float>      comp_filter(0.19f);
        MadgwickFilter<float>           madgwick_filter(0.1f);
        MahonyFilter<float>             mahony_filter(0.5f, 0.0f);
        ESKFOrientationFilter<float, 6> eskf_filter;

        // Every filter takes specific force (+z at rest in ENU).
        Vec3<float> accel_sf{0.0f, 0.0f, 9.81f};
        Vec3<float> gyro{0.1f, 0.05f, -0.08f};
        Vec3<float> mag{0.0f, 1.0f, 0.0f};
        float       dt = 0.01f;

        comp_filter.step(accel_sf, gyro, dt);
        madgwick_filter.step(accel_sf, gyro, mag, dt);
        mahony_filter.step(accel_sf, gyro, mag, dt);
        eskf_filter.step(accel_sf, gyro, mag, dt);

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
            filter.step(accel, gyro, mag, dt);
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
        ComplementaryFilter<float> filter(0.09f);

        Vec3<float> true_accel{0.0f, 0.0f, 9.81f};
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

            filter.step(noisy_accel, noisy_gyro, dt);

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
        auto q_true = Quaternion<float>::from_axis_angle(Vec3<float>{1.0f, 0.0f, 0.0f}, 0.4363f).value();
        // Body-frame specific force for a body→world attitude.
        Vec3<float> accel = q_true.conjugate().rotate(Vec3<float>{0.0f, 0.0f, 9.81f});

        ComplementaryFilter<float> filter(0.045f);
        for (int i = 0; i < 4000; ++i) {
            filter.step(accel, Vec3<float>{0.0f, 0.0f, 0.0f}, 0.005f);
        }

        // World up seen in the body must match (tilt only; yaw is unobserved).
        Vec3<float> up_est = filter.orientation().conjugate().rotate(Vec3<float>{0.0f, 0.0f, 1.0f});
        Vec3<float> up_true = q_true.conjugate().rotate(Vec3<float>{0.0f, 0.0f, 1.0f});
        CHECK((up_est - up_true).norm() < 3e-2f);
    }

    TEST_CASE("Mahony filter recovers tilt and heading from accel+mag") {
        // 25 deg about X then 40 deg about Z (yaw): heading needs the magnetometer.
        auto              qx = Quaternion<float>::from_axis_angle(Vec3<float>{1.0f, 0.0f, 0.0f}, 0.4363f).value();
        auto              qz = Quaternion<float>::from_axis_angle(Vec3<float>{0.0f, 0.0f, 1.0f}, 0.6981f).value();
        Quaternion<float> q_true = (qz * qx).normalized();

        Vec3<float> accel = q_true.conjugate().rotate(Vec3<float>{0.0f, 0.0f, 9.81f});
        Vec3<float> mag = q_true.conjugate().rotate(Vec3<float>{0.0f, 1.0f, 0.0f}); // north = +Y

        MahonyFilter<float> filter(8.0f, 0.0f); // strong proportional gain for the test
        for (int i = 0; i < 6000; ++i) {
            filter.step(accel, Vec3<float>{0.0f, 0.0f, 0.0f}, mag, 0.005f);
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

        Vec3<float> accel = q_true.conjugate().rotate(Vec3<float>{0.0f, 0.0f, 9.81f});
        Vec3<float> mag = q_true.conjugate().rotate(Vec3<float>{0.0f, 1.0f, 0.0f});

        MadgwickFilter<float> filter(5.0f); // large beta for fast test convergence
        for (int i = 0; i < 8000; ++i) {
            filter.step(accel, Vec3<float>{0.0f, 0.0f, 0.0f}, mag, 0.005f);
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
            filter.step(accel, Vec3<float>{0.0f, 0.0f, 0.0f}, mag, 0.01f, g_nav, m_nav);
        }

        Vec3<float> up_est = filter.orientation().rotate(Vec3<float>{0.0f, 0.0f, 1.0f});
        Vec3<float> north_est = filter.orientation().rotate(Vec3<float>{0.0f, 1.0f, 0.0f});
        CHECK((up_est - q_true.rotate(Vec3<float>{0.0f, 0.0f, 1.0f})).norm() < 5e-2f);
        CHECK((north_est - q_true.rotate(Vec3<float>{0.0f, 1.0f, 0.0f})).norm() < 8e-2f);
    }


    TEST_CASE("ESKFOrientationFilter uses the full accel R block") {
        // Correlated accel noise must reach the update, not just its diagonal.
        auto d = design::eskf_imu(0.003, 0.03, 0.0001, 0.01);
        REQUIRE(d.success);
        Matrix<3, 3> R_corr = d.R;
        R_corr(0, 1) = R_corr(1, 0) = 0.6 * d.R(0, 0);
        ESKFOrientationFilter<double, 3> filt(d);
        filt.set_measurement_noise(R_corr);

        // Level and at rest: zero innovation, so δx = 0 and the reset is G = I.
        const Vec3<double> g_nav = gravity_nav(NavFrame::ENU);
        const Vec3<double> accel{0.0, 0.0, kStandardGravity<double>};
        const double       dt = 0.01;
        filt.step(accel, Vec3<double>{}, dt, g_nav);

        ErrorStateKalmanFilter<6, double> ref{d};
        Matrix<6, 6>                      F = Matrix<6, 6>::identity();
        F(0, 3) = F(1, 4) = F(2, 5) = -dt;
        ref.predict([&](double) { return ErrorStateJacobian<double, 6>{F, Matrix<6, 6>::identity()}; }, dt);
        Matrix<3, 6> H = Matrix<3, 6>::zeros();
        H.block<3, 3>(0, 0) = detail::skew_sym(accel);
        REQUIRE(ref.update_innovation(ColVec<3>{}, H, R_corr));

        for (size_t i = 0; i < 6; ++i) {
            for (size_t j = 0; j < 6; ++j) {
                CHECK(filt.covariance()(i, j) == doctest::Approx(ref.covariance()(i, j)).epsilon(1e-12));
            }
        }
    }

    TEST_CASE("Classic attitude filters track a rotating body (body to world convention)") {
        // Truth is body→world and the sensors read world vectors in the body
        // (q*·v·q), as on hardware. The filters' gyro integration (q ⊗ ω) and
        // their accel/mag models must agree on that convention.
        MahonyFilter<double>        mahony(2.0, 0.0);
        MadgwickFilter<double>      madgwick(0.5);
        ComplementaryFilter<double> comp(0.245);
        const double                dt = 0.005;
        const Vec3<double>          gyro{0.3, -0.2, 0.4};
        Quaternion<double>          q_true = Quaternion<double>::identity();
        for (int k = 0; k < 1200; ++k) {
            q_true = q_true.integrate_body_rates(gyro, dt);
            const Vec3<double> a = q_true.conjugate().rotate(Vec3<double>{0.0, 0.0, 9.81});
            const Vec3<double> m = q_true.conjugate().rotate(Vec3<double>{0.0, 1.0, 0.0});
            mahony.step(a, gyro, m, dt);
            madgwick.step(a, gyro, m, dt);
            comp.step(a, gyro, dt);
        }
        CHECK(rotation_angle_between(mahony.orientation(), q_true) < 1e-2);
        CHECK(rotation_angle_between(madgwick.orientation(), q_true) < 1e-2);
        CHECK(rotation_angle_between(comp.orientation(), q_true) < 1e-2);
    }

    TEST_CASE("Complementary filter keeps gyro-integrated yaw") {
        // Level body yawing: the accelerometer carries no heading, so the blend
        // must not drag yaw back toward zero.
        ComplementaryFilter<double> comp(0.09);
        const double                dt = 0.01;
        const double                rate = 0.5;
        for (int k = 0; k < 200; ++k) {
            comp.step(Vec3<double>{0.0, 0.0, 9.81}, Vec3<double>{0.0, 0.0, rate}, dt);
        }
        const double yaw = Euler<double, EulerOrder::ZYX>::from_quaternion(comp.orientation()).yaw();
        CHECK(yaw == doctest::Approx(rate * 200 * dt).epsilon(1e-4));
    }

    TEST_CASE("Classic attitude filters take specific force, as the ESKF does") {
        // A raw IMU reading tilted 20 deg in roll: every filter must converge to the
        // same attitude from the same accelerometer sample.
        const double             roll = 0.3491;
        const Quaternion<double> q_true = Quaternion<double>::from_axis_angle(Vec3<double>{1.0, 0.0, 0.0}, roll).value();
        const Vec3<double>       f = q_true.conjugate().rotate(Vec3<double>{0.0, 0.0, kStandardGravity<double>});
        const Vec3<double>       m = q_true.conjugate().rotate(Vec3<double>{0.0, 1.0, 0.0});
        const Vec3<double>       w{};

        ComplementaryFilter<double>      comp(0.9);
        MahonyFilter<double>             mahony(4.0, 0.0);
        MadgwickFilter<double>           madgwick(0.2);
        ESKFOrientationFilter<double, 6> eskf{design::eskf_marg(0.003, 0.03, 0.03, 0.0001, 0.01)};
        for (int k = 0; k < 3000; ++k) {
            comp.step(f, w, 0.01);
            mahony.step(f, w, m, 0.01);
            madgwick.step(f, w, m, 0.01);
            eskf.step(f, w, m, 0.01);
        }
        CHECK(rotation_angle_between(comp.orientation(), q_true) < 1e-2);
        CHECK(rotation_angle_between(mahony.orientation(), q_true) < 1e-2);
        CHECK(rotation_angle_between(madgwick.orientation(), q_true) < 1e-2);
        CHECK(rotation_angle_between(eskf.orientation(), q_true) < 1e-2);
    }

    TEST_CASE("Complementary filter skips the blend on a zero or non-finite accel") {
        // A dropped sample reads (0, 0, 0): the tilt must hold, not slerp toward 180 deg roll.
        const Quaternion<double>    q_true = Quaternion<double>::from_axis_angle(Vec3<double>{1.0, 0.0, 0.0}, 0.5).value();
        const Vec3<double>          f = q_true.conjugate().rotate(Vec3<double>{0.0, 0.0, kStandardGravity<double>});
        ComplementaryFilter<double> comp(0.98);
        for (int k = 0; k < 2000; ++k) {
            comp.step(f, Vec3<double>{}, 0.01);
        }
        REQUIRE(rotation_angle_between(comp.orientation(), q_true) < 1e-3);
        for (int k = 0; k < 100; ++k) {
            comp.step(Vec3<double>{}, Vec3<double>{}, 0.01);
        }
        comp.step(Vec3<double>{std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0}, Vec3<double>{}, 0.01);
        CHECK(rotation_angle_between(comp.orientation(), q_true) < 1e-3);
    }

    TEST_CASE("ESKFOrientationFilter::reset states the covariance") {
        // No identity default: a bare reset() would discard the design's P₀ (1 (rad/s)² gyro-bias prior).
        static_assert(!BareReset<ESKFOrientationFilter<double, 3>>);
        static_assert(!StateOnlyReset<ESKFOrientationFilter<double, 3>, Quaternion<double>>);
        const auto                       d = design::eskf_imu(0.003, 0.03, 1e-4, 0.01);
        ESKFOrientationFilter<double, 3> filt{d};
        filt.reset(Quaternion<double>::identity(), Vec3<double>{}, d.P0);
        CHECK(filt.covariance()(3, 3) == doctest::Approx(d.P0(3, 3)));
    }
} // TEST_SUITE
