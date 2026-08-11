// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <limits>

#include "damp/controllers/pid.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

/**
 * @file test_pid_mode_control.cpp
 * @brief Runtime tests for PIDController (discrete) and ContinuousPID:
 *        Auto / Tracking mode and bumpless-transfer behavior.
 */

TEST_SUITE("PID Mode Control") {
    TEST_CASE("PI: disable returns clamped u_track from control()") {
        constexpr float Ts = 0.01f;

        PIController<float> pi{design::pid(1.0f, 1.0f, 0.0f, -5.0f, 5.0f).discretize(Ts)};

        REQUIRE(pi.is_enabled());

        pi.disable(2.5f);
        CHECK(!pi.is_enabled());
        CHECK(pi.control(1.0f, 0.0f) == doctest::Approx(2.5f));

        pi.disable(100.0f);
        CHECK(pi.control(1.0f, 0.0f) == doctest::Approx(5.0f));

        pi.disable(-100.0f);
        CHECK(pi.control(1.0f, 0.0f) == doctest::Approx(-5.0f));
    }

    TEST_CASE("PI: enable after disable produces bumpless re-engagement") {
        constexpr float     Ts = 0.01f;
        PIController<float> pi{design::pid(2.0f, 5.0f, 0.0f).discretize(Ts)};

        for (int i = 0; i < 20; ++i) {
            (void)pi.control(1.0f, 0.5f);
        }
        REQUIRE(pi.integral != doctest::Approx(0.0f));

        const float u_track = 3.7f;
        pi.disable(u_track);
        for (int i = 0; i < 10; ++i) {
            const float u = pi.control(1.0f, 0.5f);
            CHECK(u == doctest::Approx(u_track));
        }

        pi.enable();
        const float u_after = pi.control(1.0f, 0.5f);
        CHECK(u_after == doctest::Approx(u_track).epsilon(1e-5f));
    }

    TEST_CASE("PID: bumpless transfer with active derivative term") {
        constexpr float      Ts = 0.01f;
        PIDController<float> pid{design::pid(1.5f, 3.0f, 0.05f).discretize(Ts)};

        for (int i = 0; i < 20; ++i) {
            const float y = 0.5f + (0.01f * static_cast<float>(i));
            (void)pid.control(1.0f, y);
        }

        const float u_track = -1.2f;
        const float r_hold = 1.0f;
        const float y_hold = 0.7f;
        pid.disable(u_track);
        for (int i = 0; i < 20; ++i) {
            const float u = pid.control(r_hold, y_hold);
            CHECK(u == doctest::Approx(u_track));
        }

        pid.enable();
        const float u_after = pid.control(r_hold, y_hold);
        CHECK(u_after == doctest::Approx(u_track).epsilon(1e-5f));
    }

    TEST_CASE("PID: tracking-mode integrator preload respects i_min/i_max") {
        constexpr float Ts = 0.01f;
        // i limits are on the integral term (same units as u); pass through discretize
        PIDController<float> pid{design::pid(
                                     1.0f,
                                     1.0f,
                                     0.0f,
                                     -std::numeric_limits<float>::infinity(),
                                     std::numeric_limits<float>::infinity(),
                                     -2.0f,
                                     2.0f
        )
                                     .discretize(Ts)};

        pid.disable(100.0f);
        (void)pid.control(0.0f, 0.0f);
        CHECK(pid.integral == doctest::Approx(2.0f));

        pid.disable(-100.0f);
        (void)pid.control(0.0f, 0.0f);
        CHECK(pid.integral == doctest::Approx(-2.0f));
    }

    TEST_CASE("P: tracking mode is degenerate but API works for generic code") {
        // P accepts continuous PIDResult (no rate-dependent term)
        PController<float> p{design::pid(2.0f, 0.0f, 0.0f, -10.0f, 10.0f)};

        REQUIRE(p.is_enabled());
        p.disable(4.2f);
        CHECK(!p.is_enabled());

        CHECK(p.control(0.0f, 0.0f) == doctest::Approx(4.2f));
        CHECK(p.control(1.0f, 0.5f) == doctest::Approx(4.2f));

        p.enable();
        CHECK(p.is_enabled());
        CHECK(p.control(1.0f, 0.0f) == doctest::Approx(2.0f));
    }

    TEST_CASE("reset preserves the runtime mode") {
        constexpr float     Ts = 0.01f;
        PIController<float> pi{design::pid(1.0f, 1.0f, 0.0f).discretize(Ts)};

        pi.disable(1.5f);
        REQUIRE(!pi.is_enabled());

        pi.reset();
        CHECK(!pi.is_enabled());
        CHECK(pi.integral == doctest::Approx(0.0f));
    }

    TEST_CASE("back_calculate and mode control compose: disable, then re-engage, then saturate") {
        constexpr float     Ts = 0.01f;
        PIController<float> pi{design::pid(
                                   1.0f,
                                   1.0f,
                                   0.0f,
                                   -1.0f,
                                   1.0f,
                                   -std::numeric_limits<float>::max(),
                                   std::numeric_limits<float>::max(),
                                   0.5f
        )
                                   .discretize(Ts)};

        pi.disable(0.5f);
        (void)pi.control(2.0f, 0.0f);
        pi.enable();
        const float u_first = pi.control(2.0f, 0.0f);
        CHECK(u_first == doctest::Approx(0.5f).epsilon(1e-5f));

        for (int i = 0; i < 50; ++i) {
            (void)pi.control(5.0f, 0.0f);
        }
        const float u_saturated = pi.control(5.0f, 0.0f);
        CHECK(u_saturated == doctest::Approx(1.0f));
    }

    TEST_CASE("ContinuousPID: non-positive Ts holds state and returns the current command") {
        ContinuousPID<double> pid{design::pid(2.0, 5.0, 0.1)};
        (void)pid.control(1.0, 0.0, 0.01);
        const double integ = pid.integral;

        const double u0 = pid.control(1.0, 0.0, 0.0);
        CHECK(pid.integral == doctest::Approx(integ));
        // Integral state is the term (already includes Ki history)
        CHECK(u0 == doctest::Approx((2.0 * 1.0) + integ));

        (void)pid.control(1.0, 0.0, -0.5);
        CHECK(pid.integral == doctest::Approx(integ));

        ContinuousPID<float> pi{design::pid(1.0f, 1.0f, 0.0f)};
        (void)pi.control(1.0f, 0.0f, 0.01f);
        const float pinteg = pi.integral;
        (void)pi.control(1.0f, 0.0f, 0.0f);
        CHECK(pi.integral == doctest::Approx(pinteg));
    }

    TEST_CASE("back_calculate is a no-op when Kbc == 0, bleeds when Kbc > 0") {
        constexpr double      Ts = 0.01;
        PIDController<double> pid{design::pid(1.0, 1.0, 0.0).discretize(Ts)};
        (void)pid.control(1.0, 0.0);
        const double integ = pid.integral;
        pid.back_calculate(5.0, 1.0); // Kbc == 0 → no-op
        CHECK(pid.integral == doctest::Approx(integ));

        constexpr double      inf = std::numeric_limits<double>::infinity();
        PIDController<double> pid2{design::pid(1.0, 1.0, 0.0, -inf, inf, -inf, inf, 0.5).discretize(Ts)};
        (void)pid2.control(1.0, 0.0);
        const double before = pid2.integral;
        pid2.back_calculate(5.0, 1.0); // u_sat < u_unsat → integral down
        CHECK(pid2.integral < before);
    }

    TEST_CASE("cross-precision converting ctor copies gains and live state") {
        constexpr double      Ts = 0.01;
        PIDController<double> d{design::pid(1.5, 2.0, 0.3).discretize(Ts)};
        d.disable(0.4);
        (void)d.control(1.0, 0.2);
        d.enable();
        (void)d.control(1.0, 0.25);

        PIDController<float> f(d);
        CHECK(f.Kp == doctest::Approx(static_cast<float>(d.Kp)));
        CHECK(f.Ki == doctest::Approx(static_cast<float>(d.Ki)));
        CHECK(f.d_a == doctest::Approx(static_cast<float>(d.d_a)));
        CHECK(f.d_b == doctest::Approx(static_cast<float>(d.d_b)));
        CHECK(f.integral == doctest::Approx(static_cast<float>(d.integral)));
        CHECK(f.prev_cr_minus_y == doctest::Approx(static_cast<float>(d.prev_cr_minus_y)));
        CHECK(f.deriv == doctest::Approx(static_cast<float>(d.deriv)));
        CHECK(f.is_enabled() == d.is_enabled());
    }

    TEST_CASE("discrete PID matches ContinuousPID at fixed Ts") {
        constexpr double Ts = 0.01;
        const auto       cont_gains = design::pid(2.0, 5.0, 0.1, -10.0, 10.0, -100.0, 100.0, 2.0, 1.0, 0.0, 0.05);

        PIDController<double> disc{cont_gains.discretize(Ts)};
        ContinuousPID<double> cont{cont_gains};

        for (int k = 0; k < 80; ++k) {
            const double r = 1.0;
            const double y = 0.01 * static_cast<double>(k);
            const double ud = disc.control(r, y);
            const double uc = cont.control(r, y, Ts);
            CHECK(ud == doctest::Approx(uc).epsilon(1e-12));
        }
    }

    TEST_CASE("PID first tick has no derivative kick; Tf attenuates the derivative") {
        constexpr double Ts = 0.01;
        constexpr double inf = std::numeric_limits<double>::infinity();

        PIDController<double> d{design::pid(0.0, 0.0, 1.0).discretize(Ts)};
        CHECK(d.control(0.0, 1.0) == doctest::Approx(0.0));    // seeded: no kick
        CHECK(d.control(0.0, 2.0) == doctest::Approx(-100.0)); // Kd*(-1)/Ts

        auto step = [&](double Tf) {
            PIDController<double> c{
                design::pid(0.0, 0.0, 1.0, -inf, inf, -inf, inf, 0.0, 1.0, 1.0, Tf).discretize(Ts)
            };
            (void)c.control(0.0, 0.0);
            return c.control(0.0, 1.0);
        };
        CHECK(std::abs(step(0.1)) < std::abs(step(0.0)));
    }

    TEST_CASE("discretize scales Ki; integral-term limits pass through") {
        constexpr double Ts = 0.01;
        const auto       cont = design::pid(1.0, 20.0, 0.0, -1.0, 1.0, -0.5, 0.5);
        const auto       disc = cont.discretize(Ts);

        CHECK(disc.Kp == doctest::Approx(1.0));
        CHECK(disc.Ki == doctest::Approx(20.0 * Ts));
        CHECK(disc.i_min == doctest::Approx(-0.5));
        CHECK(disc.i_max == doctest::Approx(0.5));
        CHECK(disc.Ts == doctest::Approx(Ts));
        CHECK(disc.d_a == doctest::Approx(0.0));
        CHECK(disc.d_b == doctest::Approx(0.0));
    }

    TEST_CASE("integral term I is independent of later Ki edits") {
        constexpr double     Ts = 0.01;
        PIController<double> pi{design::pid(0.0, 10.0, 0.0).discretize(Ts)};
        for (int k = 0; k < 50; ++k) {
            (void)pi.control(1.0, 0.0); // e = 1 each tick
        }
        const double I_before = pi.integral;
        REQUIRE(I_before != doctest::Approx(0.0));
        const double u_before = pi.control(0.0, 0.0); // e = 0: u = I (Kp=0)
        pi.Ki *= 2.0;
        const double u_after = pi.control(0.0, 0.0);
        CHECK(pi.integral == doctest::Approx(I_before));
        CHECK(u_after == doctest::Approx(u_before));
        CHECK(u_after == doctest::Approx(I_before));
    }

    TEST_CASE("discretize keeps default unbounded I limits for constinit") {
        constexpr double Ts = 0.01;
        constexpr auto   disc = design::pid(4.0, 8.0, 0.0).discretize(Ts);
        constexpr double lim = std::numeric_limits<double>::max();
        static_assert(disc.i_min == -lim);
        static_assert(disc.i_max == lim);
        // Full design→deploy path is a constant expression.
        constexpr auto ctrl = PIController<float>{design::pid(4.0, 8.0, 0.0).discretize(Ts).as<float>()};
        static_assert(ctrl.Kp == 4.0f);
        CHECK(ctrl.Ki == doctest::Approx(static_cast<float>(8.0 * Ts)));
    }
} // TEST_SUITE
