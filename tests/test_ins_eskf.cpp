// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "damp/backend.hpp"
#include "damp/estimation/ins_eskf.hpp"
#include "damp/estimation/ins_mechanization.hpp"
#include "damp/math/geometry.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

constexpr double g0 = kStandardGravity<double>;

[[nodiscard]] double max_abs_diff(const Vec3<double>& a, const Vec3<double>& b) {
    double m = 0.0;
    for (size_t i = 0; i < 3; ++i) {
        m = damp::max(m, damp::abs(a[i] - b[i]));
    }
    return m;
}

[[nodiscard]] ImuSample<double>
synthetic_imu(const Quaternion<double>& q, const Vec3<double>& a_nav, const Vec3<double>& omega_body, NavFrame frame) {
    const Vec3<double> g_n = gravity_nav(frame);
    const Vec3<double> a_b = q.conjugate().rotate(a_nav - g_n);
    return ImuSample<double>{.gyro = omega_body, .accel = a_b};
}

} // namespace

TEST_SUITE("INS ESKF (#26 PR2)") {

    TEST_CASE("design::ins_eskf_design fills Q R P0 (continuous–discrete)") {
        constexpr double dt = 0.01;
        constexpr auto   r = design::ins_eskf_design(0.01, 0.1, 0.001, 0.01, dt, 2.0);
        static_assert(r.success);
        CHECK(r.success);
        CHECK(r.Q(ins_err::dtheta, ins_err::dtheta) > 0.0);
        CHECK(r.Q(ins_err::dv, ins_err::dv) > 0.0);
        CHECK(r.Q(ins_err::dbg, ins_err::dbg) > 0.0);
        CHECK(r.Q(ins_err::dba, ins_err::dba) > 0.0);
        // Position from double-integrated white accel: σ_a² Δt³/3
        constexpr double sa2 = 0.1 * 0.1;
        CHECK(r.Q(ins_err::dp, ins_err::dp) == doctest::Approx(sa2 * dt * dt * dt / 3.0));
        CHECK(r.Q(ins_err::dtheta, ins_err::dbg) < 0.0); // θ–bias cross
        CHECK(r.R(0, 0) == doctest::Approx(4.0));
        CHECK(r.P0(ins_err::dp, ins_err::dp) > 0.0);

        const auto rf = r.as<float>();
        CHECK(rf.success);
        CHECK(rf.Q(0, 0) > 0.0f);

        constexpr auto bad = design::ins_eskf_design(0.01, 0.1, 0.001, 0.01, 0.0);
        static_assert(!bad.success);
    }

    TEST_CASE("ins_inject reset G is Solà half-angle") {
        InsState<double>             x{};
        ColVec<kInsErrorDim, double> dx{};
        dx[0] = 0.2; // large-ish δθ_x
        const auto G = ins_inject(x, dx);
        // δθ=(0.2,0,0): [δθ]×_{23}=−0.2, [δθ]×_{32}=+0.2 → G=I−½S
        CHECK(G(1, 2) == doctest::Approx(0.5 * 0.2));
        CHECK(G(2, 1) == doctest::Approx(-0.5 * 0.2));
        CHECK(G(1, 1) == doctest::Approx(1.0));
    }

    TEST_CASE("ins_propagate_covariance matches dense F P F^T + Q") {
        const auto ej = ins_error_jacobian(
            Quaternion<double>::identity(),
            Vec3<double>{0.1, -0.05, 0.02},
            Vec3<double>{0.0, 0.0, 9.8},
            0.01
        );
        Matrix<kInsErrorDim, kInsErrorDim, double> P = Matrix<kInsErrorDim, kInsErrorDim, double>::identity() * 0.1;
        Matrix<kInsErrorDim, kInsErrorDim, double> Q = Matrix<kInsErrorDim, kInsErrorDim, double>::identity() * 1e-4;
        Matrix<kInsErrorDim, kInsErrorDim, double> P_dense = ej.F * P * ej.F.t() + Q;
        ins_propagate_covariance(P, ej, Q);
        for (std::size_t i = 0; i < kInsErrorDim; ++i) {
            for (std::size_t j = 0; j < kInsErrorDim; ++j) {
                CHECK(P(i, j) == doctest::Approx(P_dense(i, j)).epsilon(1e-12));
            }
        }
    }

    TEST_CASE("specific_force_at_rest ENU identity is +g Up") {
        const auto a = specific_force_at_rest(Quaternion<double>::identity(), NavFrame::ENU, g0);
        CHECK(a[0] == doctest::Approx(0.0));
        CHECK(a[1] == doctest::Approx(0.0));
        CHECK(a[2] == doctest::Approx(g0)); // +Up, not −g
        const auto a_ned = specific_force_at_rest(Quaternion<double>::identity(), NavFrame::NED, g0);
        CHECK(a_ned[2] == doctest::Approx(-g0));
    }

    TEST_CASE("ins_error_jacobian is identity at dt=0") {
        const auto ej = ins_error_jacobian(
            Quaternion<double>::identity(),
            Vec3<double>{0.1, 0.0, 0.0},
            Vec3<double>{0.0, 0.0, -g0},
            0.0
        );
        for (std::size_t i = 0; i < kInsErrorDim; ++i) {
            for (std::size_t j = 0; j < kInsErrorDim; ++j) {
                const double expect = (i == j) ? 1.0 : 0.0;
                CHECK(ej.F(i, j) == doctest::Approx(expect));
            }
        }
    }

    TEST_CASE("ins_error_jacobian couples gyro bias into attitude and accel into velocity") {
        constexpr double   dt = 0.01;
        const Vec3<double> omega{0.0, 0.0, 0.5};
        const Vec3<double> a_b{1.0, 0.0, -g0};
        const auto         ej = ins_error_jacobian(Quaternion<double>::identity(), omega, a_b, dt);

        // F_θ,bg = -dt I
        CHECK(ej.F(ins_err::dtheta, ins_err::dbg) == doctest::Approx(-dt));
        CHECK(ej.F(ins_err::dtheta + 1, ins_err::dbg + 1) == doctest::Approx(-dt));

        // F_p,v = dt I
        CHECK(ej.F(ins_err::dp, ins_err::dv) == doctest::Approx(dt));

        // F_v,ba = -R dt; R=I ⇒ −dt on diagonal
        CHECK(ej.F(ins_err::dv, ins_err::dba) == doctest::Approx(-dt));

        // Attitude block includes −[ω]× dt: F_θθ(0,1) = −(−ωz) dt = ωz dt for ω=(0,0,ωz)?
        // [ω]× with ωz: row0 col1 = −ωz ⇒ I − [ω]×dt has (0,1) = +ωz dt
        CHECK(ej.F(0, 1) == doctest::Approx(omega[2] * dt));
    }

    TEST_CASE("ins_inject right-multiplies attitude and adds linear states") {
        InsState<double> x{};
        x.q = Quaternion<double>::identity();
        x.v = Vec3<double>{1.0, 2.0, 3.0};
        x.p = Vec3<double>{4.0, 5.0, 6.0};
        x.b_g = Vec3<double>{0.01, 0.0, 0.0};
        x.b_a = Vec3<double>{0.0, 0.02, 0.0};

        ColVec<kInsErrorDim, double> dx{};
        dx[2] = 0.1; // small yaw error [rad]
        dx[3] = 0.5; // δv_x
        dx[6] = 1.0; // δp_x
        dx[9] = 0.001;
        dx[13] = 0.003;

        (void)ins_inject(x, dx);

        CHECK(x.v[0] == doctest::Approx(1.5));
        CHECK(x.p[0] == doctest::Approx(5.0));
        CHECK(x.b_g[0] == doctest::Approx(0.011));
        CHECK(x.b_a[1] == doctest::Approx(0.023));
        // yaw about body z: q ≈ [cos(0.05), 0, 0, sin(0.05)]
        CHECK(x.q.w() == doctest::Approx(damp::cos(0.05)).epsilon(1e-6));
        CHECK(x.q.z() == doctest::Approx(damp::sin(0.05)).epsilon(1e-6));
        CHECK(x.q.norm() == doctest::Approx(1.0).epsilon(1e-9));
    }

    TEST_CASE("position-aided ESKF bounds accel-bias drift vs free-running mechanization") {
        constexpr NavFrame frame = NavFrame::ENU;
        constexpr double   dt = 0.01;
        constexpr double   t_end = 30.0;
        constexpr int      n = static_cast<int>(t_end / dt);
        constexpr int      aid_every = 50; // 2 Hz absolute position
        const Vec3<double> sensor_ba{0.05, 0.0, 0.0};

        auto a_nav = [](double t) -> Vec3<double> {
            if (t < 3.0) {
                return {1.0, 0.0, 0.0};
            }
            return {0.0, 0.0, 0.0};
        };
        auto omega = [](double t) -> Vec3<double> {
            if (t >= 5.0 && t < 8.0) {
                return {0.0, 0.0, 0.2};
            }
            return {0.0, 0.0, 0.0};
        };

        InsState<double> truth{};
        truth.q = Quaternion<double>::identity();

        InsState<double> free_run{};
        free_run.q = Quaternion<double>::identity();

        InsState<double> aided{};
        aided.q = Quaternion<double>::identity();

        // Low velocity process noise so residual accel is attributed to b_a.
        auto design = design::ins_eskf_design(
            1e-3, // gyro nd
            1e-3, // accel nd (small → trust IMU, learn bias)
            1e-4, // bg rw
            1e-3, // ba rw
            dt,
            0.2,  // pos noise std [m]
            0.05, // init att
            0.2,  // init vel
            1.0,  // init pos
            0.01, // init bg
            0.1   // init ba (wide)
        );
        ErrorStateKalmanFilter<kInsErrorDim, kInsMeasDim, double> eskf(design.P0, design.Q, design.R);

        double err_free = 0.0;
        double err_aid = 0.0;

        for (int k = 0; k < n; ++k) {
            const double       t = static_cast<double>(k) * dt;
            const Vec3<double> a_n = a_nav(t);
            const Vec3<double> w_b = omega(t);

            const ImuSample<double> imu_true = synthetic_imu(truth.q, a_n, w_b, frame);
            truth = mechanize_step(truth, imu_true, dt, frame);

            ImuSample<double> imu_biased = synthetic_imu(free_run.q, a_n, w_b, frame);
            imu_biased.accel = imu_biased.accel + sensor_ba;
            free_run = mechanize_step(free_run, imu_biased, dt, frame);

            ImuSample<double> imu_aid = synthetic_imu(aided.q, a_n, w_b, frame);
            imu_aid.accel = imu_aid.accel + sensor_ba;
            ins_predict(eskf, aided, imu_aid, dt, frame);

            if ((k % aid_every) == 0) {
                // Perfect absolute position (mocap-style) — bounds filter
                CHECK(ins_update_position(eskf, aided, truth.p));
            }

            err_free = damp::max(err_free, max_abs_diff(free_run.p, truth.p));
            err_aid = damp::max(err_aid, max_abs_diff(aided.p, truth.p));
        }

        // Free-running with 0.05 m/s² bias over 30 s drifts many metres.
        CHECK(err_free > 10.0);
        // Aided filter stays within a few metres of truth (and well under free-run).
        CHECK(err_aid < 2.0);
        CHECK(err_aid < 0.2 * err_free);

        // Position aids make b_a partially observable; expect the right sign / bulk.
        CHECK(aided.b_a[0] > 0.02);
        CHECK(damp::abs(aided.b_a[0] - sensor_ba[0]) < 0.04);
    }

    TEST_CASE("ZUPT freezes a stationary biased IMU") {
        constexpr NavFrame frame = NavFrame::ENU;
        constexpr double   dt = 0.01;
        const Vec3<double> sensor_ba{0.02, -0.01, 0.0};

        InsState<double> x{};
        x.q = Quaternion<double>::identity();

        const auto                                                design = design::ins_eskf_design(1e-3, 1e-2, 1e-4, 1e-3, dt, 1.0);
        ErrorStateKalmanFilter<kInsErrorDim, kInsMeasDim, double> eskf(design.P0, design.Q, design.R);

        const Matrix<3, 3, double> R_zupt = Matrix<3, 3, double>::identity() * 1e-4;

        for (int k = 0; k < 500; ++k) {
            ImuSample<double> imu{
                .gyro = {0.0, 0.0, 0.0},
                .accel = specific_force_at_rest(x.q, frame) + sensor_ba,
            };
            ins_predict(eskf, x, imu, dt, frame);
            if ((k % 10) == 0) {
                CHECK(ins_update_zupt(eskf, x, R_zupt));
            }
        }

        // Without ZUPT, residual accel would ramp velocity; with ZUPT, |v| stays small.
        CHECK(max_abs_diff(x.v, Vec3<double>{}) < 0.05);
        CHECK(max_abs_diff(x.p, Vec3<double>{}) < 0.5);
    }

    TEST_CASE("float instantiation compiles") {
        const auto                                               design = design::ins_eskf_design(0.01f, 0.1f, 0.001f, 0.01f, 0.01f).as<float>();
        ErrorStateKalmanFilter<kInsErrorDim, kInsMeasDim, float> eskf(design.P0, design.Q, design.R);
        InsState<float>                                          x{};
        x.q = Quaternion<float>::identity();
        const ImuSample<float> imu{
            .gyro = {0.f, 0.f, 0.f},
            .accel = specific_force_at_rest(x.q, NavFrame::ENU),
        };
        ins_predict(eskf, x, imu, 0.01f, NavFrame::ENU);
        CHECK(ins_update_position(eskf, x, Vec3<float>{0.f, 0.f, 0.f}));
        CHECK(x.q.norm() == doctest::Approx(1.0f).epsilon(1e-5f));
    }
}

TEST_SUITE("INS ESKF (#26 PR3 heading / pose)") {

    TEST_CASE("heading_from_baseline_nav matches ENU and NED conventions") {
        // Pointing North: ENU forward = (0,1,0), NED forward = (1,0,0)
        CHECK(heading_from_baseline_nav(Vec3<double>{0.0, 1.0, 0.0}, NavFrame::ENU) == doctest::Approx(0.0).epsilon(1e-12));
        CHECK(heading_from_baseline_nav(Vec3<double>{1.0, 0.0, 0.0}, NavFrame::NED) == doctest::Approx(0.0).epsilon(1e-12));
        // Pointing East: ψ = +π/2
        CHECK(heading_from_baseline_nav(Vec3<double>{1.0, 0.0, 0.0}, NavFrame::ENU) == doctest::Approx(damp::numbers::pi_v<double> / 2.0).epsilon(1e-12));
        CHECK(heading_from_baseline_nav(Vec3<double>{0.0, 1.0, 0.0}, NavFrame::NED) == doctest::Approx(damp::numbers::pi_v<double> / 2.0).epsilon(1e-12));
    }

    TEST_CASE("ins_heading matches heading_from_baseline_nav for any attitude") {
        constexpr NavFrame frame = NavFrame::ENU;
        InsState<double>   x{};
        // Arbitrary level yaw: body x points somewhere in the horizontal plane
        x.q = Quaternion<double>::from_axis_angle(Vec3<double>{0.0, 0.0, 1.0}, 0.7).value();
        const Vec3<double> baseline_nav = x.q.rotate(Vec3<double>{1.0, 0.0, 0.0});
        CHECK(ins_heading(x, frame) == doctest::Approx(heading_from_baseline_nav(baseline_nav, frame)).epsilon(1e-12));
        // NED identity: body +x = North → heading 0
        x.q = Quaternion<double>::identity();
        CHECK(ins_heading(x, NavFrame::NED) == doctest::Approx(0.0).epsilon(1e-12));
    }

    TEST_CASE("heading update corrects a static yaw error") {
        constexpr NavFrame frame = NavFrame::ENU;
        constexpr double   dt = 0.01;
        const double       psi_true = 0.5;

        InsState<double> x{};
        x.q = Quaternion<double>::identity(); // wrong heading

        auto                                                      design = design::ins_eskf_design(1e-3, 1e-2, 1e-4, 1e-3, dt, 1.0, 0.5);
        ErrorStateKalmanFilter<kInsErrorDim, kInsMeasDim, double> eskf(design.P0, design.Q, design.R);
        const Matrix<1, 1, double>                                R_psi{{0.01 * 0.01}};

        for (int k = 0; k < 40; ++k) {
            ImuSample<double> imu{
                .gyro = {0.0, 0.0, 0.0},
                .accel = specific_force_at_rest(x.q, frame),
            };
            ins_predict(eskf, x, imu, dt, frame);
            CHECK(ins_update_heading(eskf, x, psi_true, R_psi, frame));
        }

        CHECK(ins_heading(x, frame) == doctest::Approx(psi_true).epsilon(1e-2));
    }

    TEST_CASE("lever-arm position update corrects for antenna offset") {
        constexpr NavFrame frame = NavFrame::ENU;
        constexpr double   dt = 0.01;
        const Vec3<double> lever{1.0, 0.0, 0.0}; // antenna 1 m forward of IMU

        InsState<double> truth{};
        truth.q = Quaternion<double>::identity();
        truth.p = Vec3<double>{2.0, 3.0, 0.0};

        InsState<double> x = truth;
        x.p = Vec3<double>{0.0, 0.0, 0.0}; // bad position prior

        auto                                                      design = design::ins_eskf_design(1e-3, 1e-2, 1e-4, 1e-3, dt, 0.05, 0.05, 0.5, 5.0);
        ErrorStateKalmanFilter<kInsErrorDim, kInsMeasDim, double> eskf(design.P0, design.Q, design.R);

        const Vec3<double> p_ant = ins_aid_position(truth, lever);
        for (int k = 0; k < 20; ++k) {
            ImuSample<double> imu{
                .gyro = {0.0, 0.0, 0.0},
                .accel = specific_force_at_rest(x.q, frame),
            };
            ins_predict(eskf, x, imu, dt, frame);
            CHECK(ins_update_position(eskf, x, p_ant, design.R, lever));
        }

        CHECK(max_abs_diff(x.p, truth.p) < 0.1);
    }

    TEST_CASE("ins_update_pose (position + orientation) snaps to mocap") {
        constexpr NavFrame frame = NavFrame::ENU;
        constexpr double   dt = 0.01;

        InsState<double> truth{};
        truth.q = Quaternion<double>::from_axis_angle(Vec3<double>{0.0, 0.0, 1.0}, 0.3).value();
        truth.p = Vec3<double>{1.0, -2.0, 0.5};

        InsState<double> x{};
        x.q = Quaternion<double>::identity();

        auto                                                      design = design::ins_eskf_design(1e-3, 1e-2, 1e-4, 1e-3, dt, 0.05, 0.5, 1.0, 5.0);
        ErrorStateKalmanFilter<kInsErrorDim, kInsMeasDim, double> eskf(design.P0, design.Q, design.R);
        const Matrix<3, 3, double>                                R_att = Matrix<3, 3, double>::identity() * (0.02 * 0.02);

        for (int k = 0; k < 15; ++k) {
            ImuSample<double> imu{
                .gyro = {0.0, 0.0, 0.0},
                .accel = specific_force_at_rest(x.q, frame),
            };
            ins_predict(eskf, x, imu, dt, frame);
            CHECK(ins_update_pose(eskf, x, truth.p, truth.q, design.R, R_att));
        }

        CHECK(max_abs_diff(x.p, truth.p) < 0.05);
        CHECK(x.q.w() == doctest::Approx(truth.q.w()).epsilon(1e-2));
        CHECK(x.q.z() == doctest::Approx(truth.q.z()).epsilon(1e-2));
    }
}

TEST_SUITE("INS ESKF (#26 PR4 InsNavigator)") {

    TEST_CASE("InsNavigator DiD path: design → predict → position + heading") {
        constexpr NavFrame frame = NavFrame::ENU;
        constexpr double   dt = 0.01;
        constexpr double   t_end = 15.0;
        constexpr int      n = static_cast<int>(t_end / dt);
        const Vec3<double> sensor_ba{0.04, 0.0, 0.0};

        auto a_nav = [](double t) -> Vec3<double> {
            return (t < 2.0) ? Vec3<double>{0.8, 0.0, 0.0} : Vec3<double>{};
        };
        auto omega = [](double t) -> Vec3<double> {
            return (t >= 3.0 && t < 6.0) ? Vec3<double>{0.0, 0.0, 0.25} : Vec3<double>{};
        };

        InsState<double> truth{};
        truth.q = Quaternion<double>::identity();

        const auto design = design::ins_eskf_design(
            1e-3, 1e-3, 1e-4, 1e-3, dt, 0.2, 0.1, 0.2, 1.0, 0.01, 0.1
        );
        InsNavigator<double> nav(design, frame);

        InsState<double> free_run{};
        free_run.q = Quaternion<double>::identity();

        const Matrix<1, 1, double> R_psi{{0.02 * 0.02}};
        double                     err_free = 0.0;
        double                     err_nav = 0.0;

        for (int k = 0; k < n; ++k) {
            const double       t = static_cast<double>(k) * dt;
            const Vec3<double> a_n = a_nav(t);
            const Vec3<double> w_b = omega(t);

            const ImuSample<double> imu_true = synthetic_imu(truth.q, a_n, w_b, frame);
            truth = mechanize_step(truth, imu_true, dt, frame);

            ImuSample<double> imu_b = synthetic_imu(free_run.q, a_n, w_b, frame);
            imu_b.accel = imu_b.accel + sensor_ba;
            free_run = mechanize_step(free_run, imu_b, dt, frame);

            ImuSample<double> imu_n = synthetic_imu(nav.state().q, a_n, w_b, frame);
            imu_n.accel = imu_n.accel + sensor_ba;
            nav.predict(imu_n, dt);

            if ((k % 50) == 0) {
                CHECK(nav.update_pose_heading(truth.p, ins_heading(truth, frame), design.R, R_psi));
            }

            err_free = damp::max(err_free, max_abs_diff(free_run.p, truth.p));
            err_nav = damp::max(err_nav, max_abs_diff(nav.state().p, truth.p));
        }

        CHECK(err_free > 3.0);
        CHECK(err_nav < 1.0);
        CHECK(err_nav < 0.25 * err_free);
        CHECK(damp::abs(nav.heading() - ins_heading(truth, frame)) < 0.15);
    }

    TEST_CASE("InsNavigator as<float> and constinit-shaped design path") {
        constexpr auto design_d = design::ins_eskf_design(1e-3, 1e-2, 1e-4, 1e-3, 0.01);
        static_assert(design_d.success);
        const auto             design_f = design_d.as<float>();
        InsNavigator<float>    nav(design_f, NavFrame::ENU);
        const ImuSample<float> imu{
            .gyro = {0.f, 0.f, 0.f},
            .accel = specific_force_at_rest(nav.state().q, NavFrame::ENU),
        };
        nav.predict(imu, 0.01f);
        CHECK(nav.update_position(Vec3<float>{0.f, 0.f, 0.f}));
        CHECK(nav.state().q.norm() == doctest::Approx(1.0f).epsilon(1e-5f));

        const auto nav_d = nav.as<double>();
        CHECK(nav_d.frame() == NavFrame::ENU);
    }
}
