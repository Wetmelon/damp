// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Regression tests for the estimation-header review (group g4).

#include <cmath>
#include <type_traits>

#include "damp/estimation/dob.hpp"
#include "damp/estimation/eskf.hpp"
#include "damp/estimation/frequency_response.hpp"
#include "damp/estimation/ins_eskf.hpp"
#include "damp/estimation/ins_mechanization.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/estimation/mhe.hpp"
#include "damp/estimation/relay_autotune.hpp"
#include "damp/estimation/rls.hpp"
#include "damp/estimation/sensor_fusion.hpp"
#include "damp/estimation/ukf.hpp"
#include "damp/math/geometry.hpp"
#include "damp/math/math.hpp"
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

[[nodiscard]] double rotation_angle_between(const Quaternion<double>& a, const Quaternion<double>& b) {
    const double w = damp::abs((a.conjugate() * b).w());
    return 2.0 * std::acos(damp::min(1.0, w));
}

} // namespace

TEST_SUITE("Review g4 - estimation") {
    TEST_CASE("ESKF predict carries an uninjected error mean through F") {
        // δx after an update but before injection is the error mean; predict must
        // propagate it with F, or the next update stacks on a stale correction.
        const Matrix<2, 2>                   P0 = Matrix<2, 2>::identity();
        const Matrix<2, 2>                   Q = Matrix<2, 2>::identity() * 1e-3;
        const Matrix<1, 1>                   R{{0.1}};
        ErrorStateKalmanFilter<2, 1, double> eskf{P0, Q, R};

        const Matrix<1, 2> H{{0.0, 1.0}};
        REQUIRE(eskf.update_innovation(ColVec<1>{{0.5}}, H, R));
        const ColVec<2> dx = eskf.error_state();
        REQUIRE(dx[1] != 0.0);

        const Matrix<2, 2> F{{1.0, 0.1}, {0.0, 1.0}};
        eskf.predict([&](double) { return ErrorStateJacobian<double, 2>{F, Matrix<2, 2>::identity()}; }, 0.1);
        const ColVec<2> expected = F * dx;
        CHECK(eskf.error_state()[0] == doctest::Approx(expected[0]));
        CHECK(eskf.error_state()[1] == doctest::Approx(expected[1]));
    }

    TEST_CASE("ESKF reset_error_state keeps P exactly symmetric") {
        const Matrix<3, 3>                   P{{2.3, 0.37, -0.19}, {0.37, 1.7, 0.29}, {-0.19, 0.29, 0.93}};
        ErrorStateKalmanFilter<3, 1, double> eskf{P, Matrix<3, 3>::identity(), Matrix<1, 1>{{1.0}}};
        const Vec3<double>                   dth{0.031, -0.017, 0.043};
        Matrix<3, 3>                         G = detail::attitude_reset_jacobian(dth);
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

    TEST_CASE("INS predict keeps the error covariance exactly symmetric") {
        const auto d = design::ins_eskf_design(1e-3, 1e-2, 1e-4, 1e-3, 0.01);
        REQUIRE(d.success);
        InsNavigator<double> nav(d, NavFrame::NED);
        InsState<double>     x{};
        x.q = Quaternion<double>::from_axis_angle(Vec3<double>{0.3, -0.5, 0.8}, 0.7).value();
        nav.set_state(x);
        const ImuSample<double> imu{.gyro = {0.11, -0.07, 0.23}, .accel = {0.4, -0.3, -9.7}};
        for (int k = 0; k < 50; ++k) {
            nav.predict(imu, 0.01);
        }
        CHECK(exactly_symmetric(nav.filter().covariance()));
    }

    TEST_CASE("INS heading Jacobian is invariant to the scale of the forward vector") {
        const Vec3<double> f{0.6, 0.8, 0.0};
        const double       s = 1e-13; // tiny but perfectly horizontal
        const Vec3<double> fs = f * s;
        for (const NavFrame frame : {NavFrame::NED, NavFrame::ENU}) {
            const auto J = detail::dheading_dforward(f, frame);
            const auto Js = detail::dheading_dforward(fs, frame);
            CHECK(J(0, 0) != 0.0);
            CHECK(Js(0, 0) * s == doctest::Approx(J(0, 0)));
            CHECK(Js(0, 1) * s == doctest::Approx(J(0, 1)));
        }
        // Exactly vertical: no heading information, zero Jacobian (no division by zero).
        const auto Jv = detail::dheading_dforward(Vec3<double>{0.0, 0.0, 1.0}, NavFrame::NED);
        CHECK(Jv(0, 0) == 0.0);
        CHECK(Jv(0, 1) == 0.0);
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
        filt.update(accel, Vec3<double>{}, dt, g_nav);

        ErrorStateKalmanFilter<6, 3, double> ref{d.P0, d.Q, R_corr};
        Matrix<6, 6>                         F = Matrix<6, 6>::identity();
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
        ComplementaryFilter<double> comp(0.98);
        const double                dt = 0.005;
        const Vec3<double>          gyro{0.3, -0.2, 0.4};
        Quaternion<double>          q_true = Quaternion<double>::identity();
        for (int k = 0; k < 1200; ++k) {
            q_true = q_true.integrate_body_rates(gyro, dt);
            const Vec3<double> a = q_true.conjugate().rotate(Vec3<double>{0.0, 0.0, -9.81});
            const Vec3<double> m = q_true.conjugate().rotate(Vec3<double>{0.0, 1.0, 0.0});
            mahony.update(a, gyro, m, dt);
            madgwick.update(a, gyro, m, dt);
            comp.update(a, gyro, dt);
        }
        CHECK(rotation_angle_between(mahony.orientation(), q_true) < 1e-2);
        CHECK(rotation_angle_between(madgwick.orientation(), q_true) < 1e-2);
        CHECK(rotation_angle_between(comp.orientation(), q_true) < 1e-2);
    }

    TEST_CASE("Complementary filter keeps gyro-integrated yaw") {
        // Level body yawing: the accelerometer carries no heading, so the blend
        // must not drag yaw back toward zero.
        ComplementaryFilter<double> comp(0.9);
        const double                dt = 0.01;
        const double                rate = 0.5;
        for (int k = 0; k < 200; ++k) {
            comp.update(Vec3<double>{0.0, 0.0, -9.81}, Vec3<double>{0.0, 0.0, rate}, dt);
        }
        const double yaw = Euler<double, EulerOrder::ZYX>::from_quaternion(comp.orientation()).yaw();
        CHECK(yaw == doctest::Approx(rate * 200 * dt).epsilon(1e-4));
    }

    TEST_CASE("UKF alias names the same type as UnscentedKalmanFilter") {
        CHECK(std::is_same_v<UKF<2, 1, 1>, UnscentedKalmanFilter<2, 1, 1>>);
    }
    TEST_CASE("RLS forgetting does not wind the covariance up without excitation") {
        // λ < 1 with no information divides P by λ every sample: unbounded growth
        // (inf in float after a few hundred samples, then NaN gains).
        estimation::RlsConfig<float> cfg{};
        cfg.lambda = 0.9f;
        cfg.p0 = 10.0f;
        estimation::Rls<float>          rls(cfg);
        estimation::RlsVector<2, float> rlsv(cfg);
        for (int k = 0; k < 1000; ++k) {
            REQUIRE(rls.update(0.0f, 1.0f));
            REQUIRE(rlsv.update(ColVec<2, float>{0.0f, 0.0f}, 1.0f));
        }
        CHECK(rls.state().covariance <= cfg.p0);
        CHECK(rlsv.state().covariance.trace() <= 2.0f * cfg.p0 * (1.0f + 1e-6f));
        // Still able to learn afterwards.
        for (int k = 0; k < 50; ++k) {
            REQUIRE(rls.update(1.0f, 2.0f));
        }
        CHECK(rls.state().theta == doctest::Approx(2.0f).epsilon(1e-3));
    }

    TEST_CASE("RLS accepts a valid tiny forgetting factor") {
        // λ + φᵀPφ ≥ λ > 0 for any valid config; an absolute threshold rejected it.
        estimation::RlsConfig<double> cfg{};
        cfg.lambda = 1e-13;
        cfg.p0 = 1.0;
        estimation::Rls<double> rls(cfg);
        CHECK(rls.update(0.0, 1.0));
    }
    TEST_CASE("Unconstrained MHE endpoint equals the steady-state Kalman filter") {
        // Arrival cost weight is P∞⁻¹ (a-priori). Its mean must be the matching
        // a-priori prediction, not the previous smoothed estimate (which already
        // used the measurements still inside the window: double counting).
        constexpr double                  Ts = 0.1;
        StateSpace<2, 1, 1, double, 2, 1> plant{};
        plant.A = Matrix<2, 2>{{1.0, Ts}, {0.0, 1.0}};
        plant.B = Matrix<2, 1>{{0.5 * Ts * Ts}, {Ts}};
        plant.C = Matrix<1, 2>{{1.0, 0.0}};
        plant.Ts = Ts;
        const Matrix<2, 2> Q{{1.0e-4, 0.0}, {0.0, 1.0e-3}};
        const Matrix<1, 1> R{{2.5e-3}};

        const auto art = design::mhe<8>(plant, Q, R);
        REQUIRE(art.success);
        MHE        estimator{art};
        const auto kf_design = design::kalman(plant, Q, R);
        REQUIRE(kf_design.success);
        KalmanFilter<2, 1, 1, double, 2, 1> kf{kf_design};

        ColVec<2> x_true{0.2, -0.1};
        double    max_err = 0.0;
        for (size_t k = 0; k < 120; ++k) {
            const double u = 0.3 * damp::sin(0.7 * static_cast<double>(k));
            x_true = plant.A * x_true + plant.B * ColVec<1>{u};
            const ColVec<1> y{x_true(0) + 0.05 * damp::sin(2.3 * static_cast<double>(k))};
            const auto      x_mhe = estimator.update(y, ColVec<1>{u});
            kf.predict(ColVec<1>{u});
            kf.update(y, ColVec<1>{u});
            if (k >= 80) { // priming transient decays geometrically
                max_err = damp::max(max_err, (x_mhe - kf.state()).norm());
            }
        }
        CHECK(max_err < 1e-6);
    }
    TEST_CASE("Asymmetric relay estimates the process static gain y/u") {
        // FOPDT plant, K = 2: the biased relay's mean output over whole limit
        // cycles is K times its mean command, so Ks must come out near 2 (not 1/2).
        constexpr double                    Ts = 0.01;
        constexpr double                    K = 2.0;
        constexpr double                    tau = 1.0;
        constexpr size_t                    delay = 30;
        design::RelayAutotuneConfig<double> cfg{};
        cfg.amplitude = 1.0;
        cfg.amplitude_neg = 0.5;
        cfg.hysteresis = 0.02;
        cfg.warmup_cycles = 4;
        cfg.measure_cycles = 8;
        cfg.max_duration = 200.0;
        RelayAutotuner<double> tuner(design::relay_autotune(cfg), Ts);

        damp::array<double, delay>  buf{};
        size_t                      head = 0;
        double                      y = 0.0;
        RelayAutotuneOutput<double> out{};
        for (int k = 0; k < 20000 && out.status != RelayAutotuneStatus::Done; ++k) {
            out = tuner.step(y);
            const double u_delayed = buf[head];
            buf[head] = out.u;
            head = (head + 1) % delay;
            y += (Ts / tau) * (-y + K * u_delayed);
        }
        REQUIRE(out.status == RelayAutotuneStatus::Done);
        CHECK(out.Ks == doctest::Approx(K).epsilon(0.15));
    }
    TEST_CASE("margins_from_frf finds the gain margin from wrapped measured phase") {
        // L = 2 / (s (s+1) (s+2)): phase −180° at ω = √2 with |L| = 1/3 → GM = 9.54 dB.
        // Measured phase is atan2-wrapped to (−π, π], so it jumps from −π to +π
        // at the crossover instead of passing through −π.
        constexpr std::size_t             Nf = 60;
        damp::array<FrfPoint<double>, Nf> table{};
        const double                      two_pi = 2.0 * damp::numbers::pi_v<double>;
        for (std::size_t i = 0; i < Nf; ++i) {
            const double w = 0.2 + 0.05 * static_cast<double>(i);
            const double mag = 2.0 / (w * damp::sqrt(1.0 + w * w) * damp::sqrt(4.0 + w * w));
            const double ph = -damp::numbers::pi_v<double> / 2.0 - damp::atan(w) - damp::atan(w / 2.0);
            table[i] = {w / two_pi, mag, damp::atan2(damp::sin(ph), damp::cos(ph)), 1.0, true};
        }
        const auto m = design::margins_from_frf(table);
        REQUIRE(m.has_gain_margin);
        CHECK(m.gain_margin_db == doctest::Approx(20.0 * damp::log10(3.0)).epsilon(0.02));
        CHECK(m.phase_crossover_hz * two_pi == doctest::Approx(damp::sqrt(2.0)).epsilon(0.02));
        REQUIRE(m.has_phase_margin);
        CHECK(m.phase_margin_deg > 0.0);
        CHECK(m.phase_margin_deg < 90.0);
    }
    TEST_CASE("classical_dob realizability check is invariant to plant gain units") {
        // Same plant, input in different units: gain 0.1 → 1e-13. Bn[0] ≠ 0 either way.
        const damp::array<double, 2> An{1.0, -0.9};
        const damp::array<double, 1> Qn{0.15};
        const damp::array<double, 2> Qd{1.0, -0.85};
        CHECK(estimation::classical_dob(damp::array<double, 1>{0.1}, An, Qn, Qd).success);
        CHECK(estimation::classical_dob(damp::array<double, 1>{1e-13}, An, Qn, Qd).success);
        // A true input delay (leading zero) is still rejected.
        CHECK_FALSE(estimation::classical_dob(damp::array<double, 2>{0.0, 1e-13}, An, Qn, Qd).success);
    }
}
