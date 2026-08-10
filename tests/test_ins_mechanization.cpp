// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "damp/backend.hpp"
#include "damp/estimation/ins_mechanization.hpp"
#include "damp/math/geometry.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"

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

} // namespace

TEST_SUITE("INS mechanization (#26 PR1)") {

    TEST_CASE("gravity_nav signs for NED and ENU") {
        const auto g_ned = gravity_nav(NavFrame::NED, g0);
        const auto g_enu = gravity_nav(NavFrame::ENU, g0);
        CHECK(g_ned[0] == doctest::Approx(0.0));
        CHECK(g_ned[1] == doctest::Approx(0.0));
        CHECK(g_ned[2] == doctest::Approx(g0));
        CHECK(g_enu[2] == doctest::Approx(-g0));
        // Same physical down direction under the axis map
        CHECK(max_abs_diff(ned_from_enu(g_enu), g_ned) < 1e-15);
        CHECK(max_abs_diff(enu_from_ned(g_ned), g_enu) < 1e-15);
    }

    TEST_CASE("ned_from_enu is an involution with enu_from_ned") {
        const Vec3<double> enu{1.0, 2.0, 3.0};
        const Vec3<double> ned = ned_from_enu(enu);
        CHECK(ned[0] == doctest::Approx(2.0)); // N ← N
        CHECK(ned[1] == doctest::Approx(1.0)); // E ← E
        CHECK(ned[2] == doctest::Approx(-3.0));
        CHECK(max_abs_diff(enu_from_ned(ned), enu) < 1e-15);
    }

    TEST_CASE("static hover: specific force cancels gravity (ENU and NED)") {
        for (const NavFrame frame : {NavFrame::ENU, NavFrame::NED}) {
            InsState<double> x{};
            x.q = Quaternion<double>::identity();
            const ImuSample<double> imu{
                .gyro = {0.0, 0.0, 0.0},
                .accel = specific_force_at_rest(x.q, frame),
            };
            constexpr double dt = 0.01;
            for (int k = 0; k < 1000; ++k) {
                x = mechanize_step(x, imu, dt, frame);
            }
            CHECK(max_abs_diff(x.v, Vec3<double>{}) < 1e-9);
            CHECK(max_abs_diff(x.p, Vec3<double>{}) < 1e-9);
            CHECK(x.q.w() == doctest::Approx(1.0).epsilon(1e-9));
        }
    }

    TEST_CASE("constant body accel: analytic position (ENU, level)") {
        // Level ENU, accel_meas = (0,0,-g) + (1,0,0) ⇒ a_n = (1,0,0)
        InsState<double> x{};
        x.q = Quaternion<double>::identity();
        const ImuSample<double> imu{
            .gyro = {0.0, 0.0, 0.0},
            .accel = specific_force_at_rest(x.q, NavFrame::ENU) + Vec3<double>{1.0, 0.0, 0.0},
        };
        constexpr double dt = 0.001;
        constexpr double T = 1.0;
        const int        n = static_cast<int>(T / dt);
        for (int k = 0; k < n; ++k) {
            x = mechanize_step(x, imu, dt, NavFrame::ENU);
        }
        // v = a T, p = ½ a T²
        CHECK(x.v[0] == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(x.p[0] == doctest::Approx(0.5).epsilon(1e-6));
        CHECK(damp::abs(x.v[1]) < 1e-9);
        CHECK(damp::abs(x.v[2]) < 1e-9);
    }

    TEST_CASE("constant body accel matches in NED after axis map") {
        // Same physical Eastward 1 m/s²: ENU +X = NED +Y
        InsState<double>        x_enu{};
        InsState<double>        x_ned{};
        const auto              a_rest_enu = specific_force_at_rest(Quaternion<double>::identity(), NavFrame::ENU);
        const auto              a_rest_ned = specific_force_at_rest(Quaternion<double>::identity(), NavFrame::NED);
        const ImuSample<double> imu_enu{.gyro = {}, .accel = a_rest_enu + Vec3<double>{1.0, 0.0, 0.0}};
        const ImuSample<double> imu_ned{.gyro = {}, .accel = a_rest_ned + Vec3<double>{0.0, 1.0, 0.0}};

        constexpr double dt = 0.001;
        for (int k = 0; k < 1000; ++k) {
            x_enu = mechanize_step(x_enu, imu_enu, dt, NavFrame::ENU);
            x_ned = mechanize_step(x_ned, imu_ned, dt, NavFrame::NED);
        }
        CHECK(max_abs_diff(ned_from_enu(x_enu.p), x_ned.p) < 1e-6);
        CHECK(max_abs_diff(ned_from_enu(x_enu.v), x_ned.v) < 1e-6);
    }

    TEST_CASE("yaw rate integrates heading (ENU level)") {
        InsState<double> x{};
        x.q = Quaternion<double>::identity();
        // ω_z body = π/2 rad/s for 1 s → 90° about Up
        const ImuSample<double> imu{
            .gyro = {0.0, 0.0, damp::numbers::pi_v<double> / 2.0},
            .accel = specific_force_at_rest(x.q, NavFrame::ENU),
        };
        constexpr double dt = 1e-3;
        for (int k = 0; k < 1000; ++k) {
            x = mechanize_step(x, imu, dt, NavFrame::ENU);
        }
        // Body +X should map toward +Y (East → North after +90° about Up in ENU:
        // right-hand: +Z up, +90° from +X(E) toward +Y(N))
        const Vec3<double> x_body{1.0, 0.0, 0.0};
        const Vec3<double> x_nav = x.q.rotate(x_body);
        CHECK(x_nav[0] == doctest::Approx(0.0).epsilon(1e-3));
        CHECK(x_nav[1] == doctest::Approx(1.0).epsilon(1e-3));
        CHECK(damp::abs(x_nav[2]) < 1e-3);
        CHECK(max_abs_diff(x.v, Vec3<double>{}) < 1e-6);
    }

    TEST_CASE("accel bias causes free-running drift (documented)") {
        // Sensor reports resting specific force + bias; filter bias state is zero → residual a.
        InsState<double> y{};
        y.q = Quaternion<double>::identity();
        const ImuSample<double> imu_biased{
            .gyro = {},
            .accel = specific_force_at_rest(y.q, NavFrame::ENU) + Vec3<double>{0.01, 0.0, 0.0},
        };
        constexpr double dt = 0.01;
        for (int k = 0; k < 1000; ++k) { // 10 s
            y = mechanize_step(y, imu_biased, dt, NavFrame::ENU);
        }
        // p ≈ ½ * 0.01 * 10² = 0.5 m
        CHECK(y.p[0] == doctest::Approx(0.5).epsilon(1e-3));
        CHECK(y.v[0] == doctest::Approx(0.1).epsilon(1e-3));
    }

    TEST_CASE("non-positive dt is a no-op") {
        InsState<double> x{};
        x.p = Vec3<double>{1.0, 2.0, 3.0};
        const auto y = mechanize_step(x, ImuSample<double>{}, 0.0, NavFrame::NED);
        CHECK(max_abs_diff(y.p, x.p) < 1e-15);
    }

    TEST_CASE("constexpr mechanize_step folds") {
        constexpr InsState<double>  x0{};
        constexpr auto              a0 = specific_force_at_rest(Quaternion<double>::identity(), NavFrame::NED);
        constexpr ImuSample<double> imu{.gyro = {}, .accel = a0};
        constexpr auto              x1 = mechanize_step(x0, imu, 0.01, NavFrame::NED);
        static_assert(x1.v[0] == 0.0);
        static_assert(x1.v[1] == 0.0);
        // residual after one step at rest should be tiny in constant evaluation
        CHECK(damp::abs(x1.v[2]) < 1e-12);
        CHECK(damp::abs(x1.p[2]) < 1e-12);
    }
}
