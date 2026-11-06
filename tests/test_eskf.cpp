// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>
#include <numbers>

#include "damp/estimation/eskf.hpp"
#include "damp/math/geometry.hpp"
#include "damp/matrix/matrix.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {
template<size_t N>
[[nodiscard]] bool exactly_symmetric(const Matrix<N, N, double>& P) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < i; ++j) {
            if (P(i, j) != P(j, i)) {
                return false;
            }
        }
    }
    return true;
}

} // namespace

constexpr double kPi = std::numbers::pi_v<double>;

// Type aliases for ESKF with 6-state error [δθ (3), δb_gyro (3)]
using Mat6d = Matrix<6, 6, double>;
using Mat3x6d = Matrix<3, 6, double>;

// Convenient aliases for the explicit Jacobian types
using PredictJac6 = ErrorStateJacobian<double, 6>;

TEST_SUITE("Error-State Kalman Filter (ESKF)") {
    // ESKF with 6D error state: [attitude(3), gyro_bias(3)]
    // Measurements: accelerometer (3D) + magnetometer (3D)

    TEST_CASE("predict G=I path matches FPF^T+Q") {
        auto                              P0 = Mat6d::identity() * 0.2;
        auto                              Q = Mat6d::identity() * 1e-5;
        ErrorStateKalmanFilter<6, double> eskf{{.Q = Q, .P0 = P0}};
        Mat6d                             F = Mat6d::identity();
        F(0, 0) = 0.99;
        F(0, 3) = -0.01;
        F(3, 3) = 1.0;
        const Mat6d G = Mat6d::identity();
        eskf.predict(
            [F, G](double) -> PredictJac6 {
            return {F, G};
        },
            0.01
        );
        const Mat6d expect = F * P0 * F.t() + Q;
        for (size_t i = 0; i < 6; ++i) {
            for (size_t j = 0; j < 6; ++j) {
                CHECK(eskf.covariance()(i, j) == doctest::Approx(expect(i, j)).epsilon(1e-12));
            }
        }
    }

    TEST_CASE("ESKF initialization") {
        // Initialize with moderate uncertainty
        auto P0 = Mat6d::identity() * 0.1; // Initial covariance
        auto Q = Mat6d::identity() * 1e-6; // Process noise (gyro drift)

        ErrorStateKalmanFilter<6, double> eskf{{.Q = Q, .P0 = P0}};

        CHECK(eskf.covariance()(0, 0) == doctest::Approx(0.1));
        CHECK(eskf.process_noise_covariance()(0, 0) == doctest::Approx(1e-6));
    }

    TEST_CASE("ESKF covariance re-initialization") {
        auto P0 = Mat6d::identity() * 0.1;
        auto Q = Mat6d::identity() * 1e-6;

        ErrorStateKalmanFilter<6, double> eskf{{.Q = Q, .P0 = P0}};

        // Reset covariance (e.g., after significant disturbance)
        auto P_reset = Mat6d::identity() * 0.2;
        eskf.reset(ColVec<6>{}, P_reset);
        CHECK(eskf.covariance()(0, 0) == doctest::Approx(0.2));
        CHECK(eskf.error_state()[0] == 0.0);
        CHECK(eskf.process_noise_covariance()(0, 0) == doctest::Approx(1e-6)); // Q stays
    }

    TEST_CASE("ESKF set_error_state bounds the correction before injection") {
        // The ESKF estimates the correction δx; the nominal state lives outside
        // the filter. Before injecting δx into the nominal quaternion/biases, the
        // caller may need to cap a component (e.g. limit an attitude-error step so
        // a faulty measurement can't slew the quaternion). set_error_state is that
        // hook.
        auto                              P0 = Mat6d::identity() * 0.1;
        auto                              Q = Mat6d::identity() * 1e-6;
        ErrorStateKalmanFilter<6, double> eskf{{.Q = Q, .P0 = P0}};

        // Hand-set a correction with an over-large attitude-error component.
        ColVec<6> dx{};
        dx[0] = 5.0; // unreasonably large attitude error [rad]
        eskf.set_error_state(dx);
        CHECK(static_cast<double>(eskf.error_state()[0]) == doctest::Approx(5.0));

        // Clamp the single component the way a caller would before injection.
        const double cap = 0.1;
        if (eskf.error_state()[0] > cap) {
            dx[0] = cap;
            eskf.set_error_state(dx);
        }
        CHECK(static_cast<double>(eskf.error_state()[0]) == doctest::Approx(cap));

        // reset_error_state still zeroes it afterward.
        eskf.reset_error_state();
        CHECK(static_cast<double>(eskf.error_state()[0]) == doctest::Approx(0.0));
    }

    TEST_CASE("ESKF predict step with proper gyro-bias coupling") {
        auto P0 = Mat6d::identity() * 0.1;
        auto Q = Mat6d::identity() * 1e-6;

        ErrorStateKalmanFilter<6, double> eskf{{.Q = Q, .P0 = P0}};

        double dt = 0.01; // 10 ms

        // Proper F matrix for ESKF error dynamics:
        // δθ_dot ≈ -δb_g  (gyro bias error causes attitude error to grow)
        // δb_g_dot ≈ 0    (bias is approximately constant)
        //
        // Discrete-time: F = I + dt * A, where A has -I in the (attitude, gyro_bias) block
        auto F = Mat6d::identity();
        // Gyro bias error propagates into attitude error
        F(0, 3) = -dt; // δb_gx → δθ_x
        F(1, 4) = -dt; // δb_gy → δθ_y
        F(2, 5) = -dt; // δb_gz → δθ_z

        auto G = Mat6d::identity() * 0.01;

        // Return explicit ErrorStateJacobian type (concept-checked)
        auto propagate = [F, G](double) -> PredictJac6 {
            return {.F = F, .G = G};
        };

        double P_before = eskf.covariance()(0, 0);
        eskf.predict(propagate, dt);
        double P_after = eskf.covariance()(0, 0);

        // Covariance should grow due to process noise
        CHECK(P_after > P_before * 0.99);

        // Cross-covariance between attitude and gyro bias should be non-zero
        // This shows the coupling is working
        double P_03 = eskf.covariance()(0, 3);
        CHECK(P_03 != doctest::Approx(0.0).epsilon(1e-10));
    }

    TEST_CASE("ESKF measurement update (accel + mag correction)") {
        auto P0 = Mat6d::identity() * 0.1;
        auto Q = Mat6d::identity() * 1e-6;
        auto R = Mat6d::identity() * 0.01;

        ErrorStateKalmanFilter<6, double> eskf{{.Q = Q, .P0 = P0}};

        // Measurements: accelerometer and magnetometer
        constexpr double g = 9.81;
        Vec3d            accel_meas{0.0, 0.0, g}; // Gravity pointing down (Z)
        Vec3d            mag_meas{0.0, 1.0, 0.0}; // Magnetic field along Y

        // Measurement Jacobian H (6x6)
        auto H = Mat6d{};
        // Accelerometer part (first 3 rows)
        H(0, 1) = -g; // δθ_y affects accel_x
        H(1, 0) = g;  // δθ_x affects accel_y
        // Magnetometer part (last 3 rows)
        H(3, 1) = -1.0; // δθ_y affects mag_x
        H(4, 0) = 1.0;  // δθ_x affects mag_y

        // Predicted measurements (from nominal state - assuming upright)
        Vec3d     accel_pred{0.0, 0.0, 9.81};
        Vec3d     mag_pred{0.0, 1.0, 0.0};
        ColVec<6> y_pred{accel_pred[0], accel_pred[1], accel_pred[2], mag_pred[0], mag_pred[1], mag_pred[2]};

        // The reading linearized at the nominal state: ν = y − ŷ and H.
        ColVec<6> y{accel_meas[0], accel_meas[1], accel_meas[2], mag_meas[0], mag_meas[1], mag_meas[2]};
        bool      success = eskf.update_innovation(y - y_pred, H, R);
        REQUIRE(success);

        // Error state should remain small (no significant error detected)
        auto err_state = eskf.error_state();
        CHECK(std::abs(static_cast<double>(err_state[0])) < 0.1); // Attitude errors small
        CHECK(std::abs(static_cast<double>(err_state[3])) < 0.1); // Gyro bias errors small

        // Covariance should shrink (measurement reduces uncertainty)
        CHECK(eskf.covariance()(0, 0) < static_cast<double>(P0(0, 0)));
    }

    TEST_CASE("ESKF with attitude error (detects misalignment)") {
        auto P0 = Mat6d::identity() * 0.1;
        auto Q = Mat6d::identity() * 1e-6;
        auto R = Mat6d::identity() * 0.01;

        ErrorStateKalmanFilter<6, double> eskf{{.Q = Q, .P0 = P0}};

        // Simulated scenario: nominal attitude is tilted, sensors say upright
        // This simulates a gyro drift that caused attitude error
        constexpr double g = 9.81;

        // Correct H matrix for accel + mag
        auto H = Mat6d{};
        H(0, 1) = -g;   // δθ_y affects accel_x
        H(1, 0) = g;    // δθ_x affects accel_y
        H(3, 1) = -1.0; // δθ_y affects mag_x
        H(4, 0) = 1.0;  // δθ_x affects mag_y


        // Measurements: upright
        Vec3d accel_meas{0.0, 0.0, 9.81};
        Vec3d mag_meas{0.0, 1.0, 0.0};

        // Predicted: tilted (wrong)
        double tilt_angle = 5.0 * kPi / 180.0; // 5 degrees
        Vec3d  accel_pred{
            9.81 * std::sin(tilt_angle),
            0.0,
            9.81 * std::cos(tilt_angle)
        };
        Vec3d mag_pred{
            -std::sin(tilt_angle),
            1.0,
            0.0
        };
        ColVec<6> y_pred{accel_pred[0], accel_pred[1], accel_pred[2], mag_pred[0], mag_pred[1], mag_pred[2]};

        // The reading linearized at the nominal state: ν = y − ŷ and H.
        ColVec<6> y{accel_meas[0], accel_meas[1], accel_meas[2], mag_meas[0], mag_meas[1], mag_meas[2]};
        bool      success = eskf.update_innovation(y - y_pred, H, R);
        REQUIRE(success);

        // Error state should detect the tilt error
        auto err_state = eskf.error_state();
        CHECK(std::abs(static_cast<double>(err_state[1])) > 0.001); // Should detect pitch/roll error
    }

    TEST_CASE("ESKF cycle: predict then update") {
        auto P0 = Mat6d::identity() * 0.1;
        auto Q = Mat6d::identity() * 1e-7;
        auto R = Mat6d::identity() * 0.001;

        ErrorStateKalmanFilter<6, double> eskf{{.Q = Q, .P0 = P0}};

        double dt = 0.01;

        // Propagation Jacobians (simple model)
        auto F = Mat6d::identity();
        F(0, 0) = 0.99;
        F(1, 1) = 0.99;
        F(2, 2) = 0.99;
        F(3, 3) = 0.999;
        F(4, 4) = 0.999;
        F(5, 5) = 0.999;

        auto G = Mat6d::identity() * 0.1;

        // Return explicit ErrorStateJacobian type (concept-checked)
        auto propagate = [F, G](double) -> PredictJac6 {
            return {.F = F, .G = G};
        };

        // Measurement - correct H for accel + mag
        constexpr double g = 9.81;
        auto             H = Mat6d{};
        H(0, 1) = -g;   // δθ_y affects accel_x
        H(1, 0) = g;    // δθ_x affects accel_y
        H(3, 1) = -1.0; // δθ_y affects mag_x
        H(4, 0) = 1.0;  // δθ_x affects mag_y

        Vec3d     accel_meas{0.0, 0.0, g};
        Vec3d     mag_meas{0.0, 1.0, 0.0};
        Vec3d     accel_pred{0.0, 0.0, g};
        Vec3d     mag_pred{0.0, 1.0, 0.0};
        ColVec<6> y_pred{accel_pred[0], accel_pred[1], accel_pred[2], mag_pred[0], mag_pred[1], mag_pred[2]};

        // The reading linearized at the nominal state: ν = y − ŷ and H.
        ColVec<6> y{accel_meas[0], accel_meas[1], accel_meas[2], mag_meas[0], mag_meas[1], mag_meas[2]};

        // Run several iterations
        for (int i = 0; i < 5; ++i) {
            eskf.predict(propagate, dt);
            [[maybe_unused]] bool success = eskf.update_innovation(y - y_pred, H, R);
            REQUIRE(success);

            // Error state should stay bounded
            auto err = eskf.error_state();
            for (size_t j = 0; j < 6; ++j) {
                CHECK(std::abs(static_cast<double>(err[j])) < 0.5);
            }
        }

        // Covariance should converge to stable value
        double final_P_00 = eskf.covariance()(0, 0);
        CHECK(final_P_00 > 0.0);
        CHECK(final_P_00 < static_cast<double>(P0(0, 0))); // Should have reduced uncertainty
    }

    TEST_CASE("ESKF predict carries an uninjected error mean through F") {
        // δx after an update but before injection is the error mean; predict must
        // propagate it with F, or the next update stacks on a stale correction.
        ErrorStateKalmanFilter<2, double> eskf{{.Q = Matrix<2, 2>::identity() * 1e-3, .P0 = Matrix<2, 2>::identity()}};

        const Matrix<1, 2> H{{0.0, 1.0}};
        REQUIRE(eskf.update_innovation(ColVec<1>{0.5}, H, Matrix<1, 1>{0.1}));
        const ColVec<2> dx = eskf.error_state();
        REQUIRE(dx[1] != 0.0);

        const Matrix<2, 2> F{{1.0, 0.1}, {0.0, 1.0}};
        eskf.predict([&](double) { return ErrorStateJacobian<double, 2>{F, Matrix<2, 2>::identity()}; }, 0.1);
        const ColVec<2> expected = F * dx;
        CHECK(eskf.error_state()[0] == doctest::Approx(expected[0]));
        CHECK(eskf.error_state()[1] == doctest::Approx(expected[1]));
    }

    TEST_CASE("ESKF reset_error_state keeps P exactly symmetric") {
        const Matrix<3, 3>                P{{2.3, 0.37, -0.19}, {0.37, 1.7, 0.29}, {-0.19, 0.29, 0.93}};
        ErrorStateKalmanFilter<3, double> eskf{{.Q = Matrix<3, 3>::identity(), .P0 = P}};
        const Vec3<double>                dth{0.031, -0.017, 0.043};
        Matrix<3, 3>                      G = detail::attitude_reset_jacobian(dth);
        G(0, 2) += 0.013; // generic (non-attitude) reset coupling
        for (int k = 0; k < 20; ++k) {
            eskf.reset_error_state(G);
        }
        CHECK(exactly_symmetric(eskf.covariance()));
    }

    TEST_CASE("ESKF attitude reset Jacobian is I - 1/2 [dtheta]x") {
        const Vec3<double> dth{0.2, -0.1, 0.3};
        const auto         G = detail::attitude_reset_jacobian(dth);
        const Vec3<double> w{0.7, 0.4, -0.5};
        // G w = w − ½ δθ × w
        const Vec3<double> Gw = G * w;
        const Vec3<double> ref = w - (dth.cross(w) * 0.5);
        for (size_t i = 0; i < 3; ++i) {
            CHECK(Gw[i] == doctest::Approx(ref[i]));
        }
    }
}

TEST_CASE("design::eskf_marg matches expected covariances") {
    constexpr double gyro_noise_density = 0.01;
    constexpr double accel_noise_density = 0.1;
    constexpr double mag_noise_density = 0.01;
    constexpr double gyro_bias_rw = 0.001;
    constexpr double dt = 0.01;
    constexpr double initial_attitude_uncertainty = 0.1;
    constexpr double initial_bias_uncertainty = 0.01;

    constexpr auto result = design::eskf_marg(
        gyro_noise_density,
        accel_noise_density,
        mag_noise_density,
        gyro_bias_rw,
        dt,
        initial_attitude_uncertainty,
        initial_bias_uncertainty
    );
    static_assert(result.success);

    constexpr double sg2 = gyro_noise_density * gyro_noise_density;
    constexpr double sbg2 = gyro_bias_rw * gyro_bias_rw;
    constexpr double expected_q_tt = sg2 * dt + sbg2 * (dt * dt * dt) / 3.0;
    constexpr double expected_bias_var = sbg2 * dt;
    constexpr double expected_q_tb = -sbg2 * (dt * dt) / 2.0;

    CHECK(doctest::Approx(result.Q(0, 0)).epsilon(1e-10) == expected_q_tt);
    CHECK(doctest::Approx(result.Q(1, 1)).epsilon(1e-10) == expected_q_tt);
    CHECK(doctest::Approx(result.Q(2, 2)).epsilon(1e-10) == expected_q_tt);
    CHECK(doctest::Approx(result.Q(3, 3)).epsilon(1e-10) == expected_bias_var);
    CHECK(doctest::Approx(result.Q(4, 4)).epsilon(1e-10) == expected_bias_var);
    CHECK(doctest::Approx(result.Q(5, 5)).epsilon(1e-10) == expected_bias_var);
    CHECK(doctest::Approx(result.Q(0, 3)).epsilon(1e-10) == expected_q_tb);

    constexpr double expected_accel_var = (accel_noise_density * accel_noise_density) / dt;
    constexpr double expected_mag_var = (mag_noise_density * mag_noise_density) / dt;
    CHECK(doctest::Approx(result.R(0, 0)).epsilon(1e-10) == expected_accel_var);
    CHECK(doctest::Approx(result.R(1, 1)).epsilon(1e-10) == expected_accel_var);
    CHECK(doctest::Approx(result.R(2, 2)).epsilon(1e-10) == expected_accel_var);
    CHECK(doctest::Approx(result.R(3, 3)).epsilon(1e-10) == expected_mag_var);
    CHECK(doctest::Approx(result.R(4, 4)).epsilon(1e-10) == expected_mag_var);
    CHECK(doctest::Approx(result.R(5, 5)).epsilon(1e-10) == expected_mag_var);

    constexpr double expected_attitude_var = initial_attitude_uncertainty * initial_attitude_uncertainty;
    constexpr double expected_bias_init_var = initial_bias_uncertainty * initial_bias_uncertainty;
    CHECK(doctest::Approx(result.P0(0, 0)).epsilon(1e-10) == expected_attitude_var);
    CHECK(doctest::Approx(result.P0(1, 1)).epsilon(1e-10) == expected_attitude_var);
    CHECK(doctest::Approx(result.P0(2, 2)).epsilon(1e-10) == expected_attitude_var);
    CHECK(doctest::Approx(result.P0(3, 3)).epsilon(1e-10) == expected_bias_init_var);
    CHECK(doctest::Approx(result.P0(4, 4)).epsilon(1e-10) == expected_bias_init_var);
    CHECK(doctest::Approx(result.P0(5, 5)).epsilon(1e-10) == expected_bias_init_var);
}

TEST_CASE("design::eskf_imu matches expected covariances (NY=3)") {
    constexpr double gyro_noise_density = 0.01;
    constexpr double accel_noise_density = 0.1;
    constexpr double gyro_bias_rw = 0.001;
    constexpr double dt = 0.01;

    constexpr auto result = design::eskf_imu(
        gyro_noise_density, accel_noise_density, gyro_bias_rw, dt, 0.1, 0.01
    );
    static_assert(result.success);

    constexpr double sg2 = gyro_noise_density * gyro_noise_density;
    constexpr double sbg2 = gyro_bias_rw * gyro_bias_rw;
    constexpr double expected_q_tt = sg2 * dt + sbg2 * (dt * dt * dt) / 3.0;
    constexpr double expected_bias_var = sbg2 * dt;
    constexpr double expected_accel_var = (accel_noise_density * accel_noise_density) / dt;

    CHECK(doctest::Approx(result.Q(0, 0)).epsilon(1e-10) == expected_q_tt);
    CHECK(doctest::Approx(result.Q(3, 3)).epsilon(1e-10) == expected_bias_var);
    CHECK(doctest::Approx(result.R(0, 0)).epsilon(1e-10) == expected_accel_var);
    CHECK(doctest::Approx(result.R(2, 2)).epsilon(1e-10) == expected_accel_var);
    CHECK(result.R.rows() == 3);
    CHECK(result.R.cols() == 3);

    constexpr auto bad = design::eskf_imu(0.01, 0.1, 0.001, 0.0);
    static_assert(!bad.success);
}

TEST_CASE("design::eskf_imu / eskf_marg refuse a zero measurement density") {
    // R = 0 with the rank-2 attitude Jacobian [a]× leaves S singular: every update would fail.
    CHECK_FALSE(design::eskf_imu(0.003, 0.0, 1e-4, 0.01).success);
    CHECK_FALSE(design::eskf_marg(0.003, 0.0, 0.03, 1e-4, 0.01).success);
    CHECK_FALSE(design::eskf_marg(0.003, 0.03, 0.0, 1e-4, 0.01).success);
    CHECK(design::eskf_imu(0.0, 0.03, 0.0, 0.01).success); // zero process noise is a valid model
}

TEST_CASE("ESKFResult::as<U>() conversion") {
    constexpr auto eskf_d = design::eskf_marg(0.01, 0.1, 0.01, 0.001, 0.01);
    constexpr auto eskf_f = eskf_d.as<float>();
    static_assert(eskf_f.success);
    CHECK(eskf_f.success);
    CHECK(eskf_f.Q(0, 0) > 0.0f);
    CHECK(eskf_f.R(0, 0) > 0.0f);
}
