// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <limits>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/encoder/unwrap.hpp"
#include "damp/math/math.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

/**
 * @file test_pid_mode_control.cpp
 * @brief Runtime tests for PIDController (discrete) and ContinuousPID:
 *        Auto / Tracking mode and bumpless-transfer behavior.
 */

template<typename C>
concept has_integral = requires(C c) { c.integral; };

template<typename C>
concept has_derivative = requires(C c) { c.deriv; };

TEST_SUITE("PID Mode Control") {
    TEST_CASE("PI: disable returns clamped u_track from control()") {
        constexpr float Ts = 0.01f;

        PIController<float> pi{design::PIDGains{.Kp = 1.0f, .Ki = 1.0f, .u_lim = {
                                                                            .lower = -5.0f,
                                                                            .upper = 5.0f,
                                                                        }}
                                   .discretize(Ts)};

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

    TEST_CASE("PID: bumpless with a filtered derivative still decaying at hand-off") {
        constexpr double Ts = 0.01;
        auto             gains = design::pid(1.5, 3.0, 0.05);
        gains.Tf = 0.05; // d_a = Tf/(Ts+Tf) ≈ 0.83: D carries over into the first Auto tick
        PIDController<double> pid{gains.discretize(Ts)};
        ContinuousPID<double> cont{gains};

        pid.disable(-1.2);
        cont.disable(-1.2);
        double y = 0.0;
        for (int i = 0; i < 20; ++i) {
            y = 0.5 + (0.02 * static_cast<double>(i)); // moving y keeps D nonzero
            CHECK(pid.control(1.0, y) == doctest::Approx(-1.2));
            CHECK(cont.control(1.0, y, Ts) == doctest::Approx(-1.2));
        }
        REQUIRE(pid.deriv != doctest::Approx(0.0));

        // Same (r, y) as the last tracking tick: the first Auto tick must emit u_track.
        pid.enable();
        cont.enable();
        CHECK(pid.control(1.0, y) == doctest::Approx(-1.2).epsilon(1e-12));
        CHECK(cont.control(1.0, y, Ts) == doctest::Approx(-1.2).epsilon(1e-12));
    }

    TEST_CASE("PD: P + filtered D, Tracking holds u_track") {
        constexpr double Ts = 0.01;
        auto             gains = design::pid(2.0, 0.0, 0.1);
        gains.Tf = 0.02;
        PDController<double> pd{gains.discretize(Ts)};

        CHECK(pd.control(1.0, 0.0) == doctest::Approx(2.0)); // first tick: no derivative kick
        const double u = pd.control(1.0, 0.1);               // y moved 0.1: D = d_b·(−0.1)
        CHECK(u == doctest::Approx((2.0 * 0.9) - ((0.1 / (Ts + 0.02)) * 0.1)));

        pd.disable(-3.0);
        CHECK(pd.control(1.0, 0.1) == doctest::Approx(-3.0));
    }

    TEST_CASE("each mode stores only its own terms") {
        // RAM is per mode: no integrator state without I, no derivative state without D.
        static_assert(has_integral<PIDController<float>> && has_derivative<PIDController<float>>);
        static_assert(has_integral<PIController<float>> && !has_derivative<PIController<float>>);
        static_assert(!has_integral<PDController<float>> && has_derivative<PDController<float>>);
        static_assert(!has_integral<PController<float>> && !has_derivative<PController<float>>);
        static_assert(sizeof(PController<float>) < sizeof(PIController<float>));
        static_assert(sizeof(PIController<float>) < sizeof(PIDController<float>));
        static_assert(sizeof(PDController<float>) < sizeof(PIDController<float>));
    }

    TEST_CASE("PID: tracking-mode integrator preload respects the integrator limits i") {
        constexpr float Ts = 0.01f;
        // i limits are on the integral term (same units as u); pass through discretize
        PIDController<float> pid{design::PIDGains{.Kp = 1.0f, .Ki = 1.0f, .integral_lim = {
                                                                              .lower = -2.0f,
                                                                              .upper = 2.0f,
                                                                          }}
                                     .discretize(Ts)};

        pid.disable(100.0f);
        (void)pid.control(0.0f, 0.0f);
        CHECK(pid.integral == doctest::Approx(2.0f));

        pid.disable(-100.0f);
        (void)pid.control(0.0f, 0.0f);
        CHECK(pid.integral == doctest::Approx(-2.0f));
    }

    TEST_CASE("P: tracking mode is degenerate but API works for generic code") {
        // P accepts continuous PIDGains (no rate-dependent term)
        PController<float> p{design::PIDGains{.Kp = 2.0f, .u_lim = {
                                                              .lower = -10.0f,
                                                              .upper = 10.0f,
                                                          }}
                                 .discretize(0.01f)};

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
        PIController<float> pi{design::PIDGains{.Kp = 1.0f, .Ki = 1.0f, .u_lim = {
                                                                            .lower = -1.0f,
                                                                            .upper = 1.0f,
                                                                        },
                                                .Kbc = 0.5f}
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

    TEST_CASE("ContinuousPID: back_calculate after a held tick is a no-op; a NaN Ts holds") {
        // Kp = 1, Ki = 10, |u| ≤ 1, Kbc = 1: a negative measured period used to charge
        // back_calculate with Ki·Ts < 0, winding the integral further into the limit.
        ContinuousPID<double> pid{design::PIDGains<double>{.Kp = 1.0, .Ki = 10.0, .u_lim = {.lower = -1.0, .upper = 1.0}, .Kbc = 1.0}};
        for (int k = 0; k < 5; ++k) {
            (void)pid.control(0.5, 0.0, 0.01);
        }
        const double integ = pid.integral;
        (void)pid.control(0.5, 0.0, -0.01);
        pid.back_calculate(pid.last_command() - 0.5);
        CHECK(pid.integral == integ);

        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double u = pid.control(0.5, 0.0, nan);
        CHECK(std::isfinite(u));
        CHECK(std::isfinite(pid.integral));
        CHECK(std::isfinite(pid.deriv));
        CHECK(std::isfinite(pid.control(0.5, 0.0, 0.01)));
    }

    TEST_CASE("discretize refuses a negative Tf and non-finite gains") {
        // Tf = −0.75·Ts gave a derivative pole of −3; Tf = −Ts divided by zero.
        constexpr double Ts = 0.01;
        const auto       neg_tf = design::PIDGains<double>{.Kp = 1.0, .Ki = 2.0, .Kd = 0.1, .Tf = -0.75 * Ts}.discretize(Ts);
        CHECK_FALSE(neg_tf.success);
        CHECK(neg_tf.status == design::DesignStatus::invalid_parameter);
        CHECK_FALSE(design::PIDGains<double>{.Kp = 1.0, .Ki = 2.0, .Kd = 0.1, .Tf = -Ts}.discretize(Ts).success);
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double inf = std::numeric_limits<double>::infinity();
        CHECK_FALSE(design::PIDGains<double>{.Kp = nan, .Ki = 2.0}.discretize(Ts).success);
        CHECK_FALSE(design::PIDGains<double>{.Kp = 1.0, .Ki = inf}.discretize(Ts).success);
        CHECK_FALSE(design::PIDGains<double>{.Kp = 1.0, .Ki = 2.0, .Kd = 0.1, .Tf = nan}.discretize(Ts).success);
        CHECK_FALSE(design::PIDGains<double>{.Kp = 1.0, .Ki = 2.0, .b = inf}.discretize(Ts).success);
        static_assert(!damp::detail::pid_gains_in_domain(1.0, 2.0, 0.1, -0.001, 0.0, 1.0, 1.0));
        static_assert(damp::detail::pid_gains_in_domain(1.0, 2.0, 0.1, 0.0, 0.0, 1.0, 1.0));
    }

    TEST_CASE("back_calculate is a no-op when Kbc == 0, bleeds when Kbc > 0") {
        constexpr double      Ts = 0.01;
        PIDController<double> pid{design::pid(1.0, 1.0, 0.0).discretize(Ts)};
        (void)pid.control(1.0, 0.0);
        const double integ = pid.integral;
        pid.back_calculate(0.5); // Kbc == 0 → no-op
        CHECK(pid.integral == doctest::Approx(integ));

        PIDController<double> pid2{design::PIDGains{.Kp = 1.0, .Ki = 1.0, .Kbc = 0.5}.discretize(Ts)};
        (void)pid2.control(1.0, 0.0);
        const double before = pid2.integral;
        pid2.back_calculate(0.5); // u_sat below the last command → integral down
        CHECK(pid2.integral < before);
    }

    TEST_CASE("a Kbc of the opposite sign to Ki is refused") {
        // T_t = Kbc/Ki < 0: ΔI = Ki·(u_sat − u_unsat)/Kbc would push the integral further into
        // the limit. discretize refuses it; constructing a runtime from such gains directly is
        // refused like a mode mismatch, ill-formed under constant evaluation and a trap at run
        // time, so that half cannot be a CHECK:
        // static constinit PIController<double> bad{design::DiscretePIDGains<double>{.Kp = 1.0, .Ki = 0.01, .Kbc = -1.0}}; // ill-formed
        constexpr auto bad = design::PIDGains{.Kp = 1.0, .Ki = 2.0, .Kbc = -1.0}.discretize(0.01);
        static_assert(!bad.success);
        CHECK(bad.status == design::DesignStatus::invalid_parameter);
        CHECK_FALSE(design::PIDGains{.Kp = 1.0, .Ki = -2.0, .Kbc = 1.0}.discretize(0.01).success);

        // Same sign, Ki = 0 (back-calculation inert) and Kbc = 0 (clamp only) all deploy; a
        // Kp < 0 law takes Kbc > 0 like any other with Ki > 0.
        static_assert(design::PIDGains{.Kp = -1.0, .Ki = 2.0, .Kbc = 1.0}.discretize(0.01).success);
        static_assert(design::PIDGains{.Kp = 1.0, .Ki = -2.0, .Kbc = -1.0}.discretize(0.01).success);
        static_assert(design::PIDGains{.Kp = 1.0, .Ki = 0.0, .Kbc = -1.0}.discretize(0.01).success);
        static_assert(design::PIDGains{.Kp = 1.0, .Ki = 2.0}.discretize(0.01).success);
        static constinit PIController<double> ok{design::PIDGains{.Kp = -1.0, .Ki = 2.0, .Kbc = 1.0}.discretize(0.01)};
        CHECK(ok.Kbc == 1.0);
        static_assert(damp::detail::pid_kbc_winds_up(2.0, -1.0));
        static_assert(damp::detail::pid_kbc_winds_up(-2.0, 1.0));
        static_assert(!damp::detail::pid_kbc_winds_up(-2.0, -1.0));
        static_assert(!damp::detail::pid_kbc_winds_up(0.0, -1.0));
    }

    TEST_CASE("back_calculate(u_sat, u_ff) reproduces the two-point ΔI = gain·(u_sat − u_unsat)/Kbc exactly") {
        // A PI plus feedforward, clamped downstream. The reference applies the formula with
        // u_unsat = the PI's command + u_ff, which is what callers passed before.
        constexpr double      Ts = 0.01;
        const auto            gains = design::PIDGains{.Kp = 2.0, .Ki = 20.0, .Kbc = 2.0}.discretize(Ts);
        PIController<double>  law{gains};
        PIController<double>  ref{gains};
        ContinuousPID<double> claw{design::PIDGains{.Kp = 2.0, .Ki = 20.0, .Kbc = 2.0}};
        ContinuousPID<double> cref{design::PIDGains{.Kp = 2.0, .Ki = 20.0, .Kbc = 2.0}};
        constexpr double      u_ff = 0.3;
        constexpr double      u_lim = 1.0;
        double                y = 0.0;
        int                   clamped = 0;
        for (int k = 0; k < 200; ++k) {
            const double r = (k < 100) ? 2.0 : -2.0;
            const double u = law.control(r, y);
            const double u_ref = ref.control(r, y);
            const double u_sat = damp::clamp(u + u_ff, -u_lim, u_lim);
            if (u_sat != u + u_ff) {
                ++clamped;
                law.back_calculate(u_sat, u_ff);
                damp::detail::pid_wind_back(ref, ref.Ki, u_sat, u_ref + u_ff);

                const double cu = claw.control(r, y, Ts);
                const double cu_ref = cref.control(r, y, Ts);
                const double cu_sat = damp::clamp(cu + u_ff, -u_lim, u_lim);
                claw.back_calculate(cu_sat, u_ff);
                damp::detail::pid_wind_back(cref, cref.Ki * Ts, cu_sat, cu_ref + u_ff);
            }
            CHECK(law.integral == ref.integral);
            CHECK(claw.integral == cref.integral);
            y += Ts * (u_sat - y); // first-order plant
        }
        CHECK(clamped > 0);
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
        CHECK(f.prev_x == doctest::Approx(static_cast<float>(d.prev_x)));
        CHECK(f.deriv == doctest::Approx(static_cast<float>(d.deriv)));
        CHECK(f.is_enabled() == d.is_enabled());
    }

    TEST_CASE("discrete PID matches ContinuousPID at fixed Ts") {
        constexpr double       Ts = 0.01;
        const design::PIDGains cont_gains{
            .Kp = 2.0,
            .Ki = 5.0,
            .Kd = 0.1,
            .u_lim = {
                .lower = -10.0,
                .upper = 10.0,
            },
            .integral_lim = {
                .lower = -100.0,
                .upper = 100.0,
            },
            .Kbc = 2.0,
            .c = 0.0,
            .Tf = 0.05
        };

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

        PIDController<double> d{design::pid(0.0, 0.0, 1.0).discretize(Ts)};
        CHECK(d.control(0.0, 1.0) == doctest::Approx(0.0));    // seeded: no kick
        CHECK(d.control(0.0, 2.0) == doctest::Approx(-100.0)); // Kd*(-1)/Ts

        auto step = [&](double Tf) {
            PIDController<double> c{
                design::PIDGains{.Kd = 1.0, .Tf = Tf}.discretize(Ts)
            };
            (void)c.control(0.0, 0.0);
            return c.control(0.0, 1.0);
        };
        CHECK(std::abs(step(0.1)) < std::abs(step(0.0)));
    }

    TEST_CASE("discretize scales Ki; integral-term limits pass through") {
        constexpr double Ts = 0.01;
        const auto       cont = design::PIDGains{.Kp = 1.0, .Ki = 20.0, .u_lim = {
                                                                            .lower = -1.0,
                                                                            .upper = 1.0,
                                                                        },
                                                 .integral_lim = {
                                                     .lower = -0.5,
                                                     .upper = 0.5,
                                                 }};
        const auto       disc = cont.discretize(Ts);

        CHECK(disc.Kp == doctest::Approx(1.0));
        CHECK(disc.Ki == doctest::Approx(20.0 * Ts));
        CHECK(disc.integral_lim.lower(0) == doctest::Approx(-0.5));
        CHECK(disc.integral_lim.upper(0) == doctest::Approx(0.5));
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
        static_assert(disc.integral_lim.lower(0) == -lim);
        static_assert(disc.integral_lim.upper(0) == lim);
        // Full design→deploy path is a constant expression.
        constexpr auto ctrl = PIController<float>{design::pid(4.0, 8.0, 0.0).discretize(Ts).as<float>()};
        static_assert(ctrl.Kp == 4.0f);
        CHECK(ctrl.Ki == doctest::Approx(static_cast<float>(8.0 * Ts)));
    }

    TEST_CASE("PIDGains designated initializers deduce T and deploy as static constinit") {
        constexpr double                     Ts = 0.01;
        static constinit PIController<float> pi = design::PIDGains{.Kp = 4.0, .Ki = 8.0, .u_lim = {
                                                                                             .lower = -12.0,
                                                                                             .upper = 12.0,
                                                                                         }}
                                                      .discretize(Ts)
                                                      .as<float>();
        static_assert(std::is_same_v<decltype(design::PIDGains{.Kp = 4.0}), design::PIDGains<double>>);
        constexpr design::PIDGains gains{.Kp = 4.0, .Ki = 8.0, .u_lim = {
                                                                   .lower = -12.0,
                                                                   .upper = 12.0,
                                                               }};
        static_assert(gains.Kd == 0.0 && gains.b == 1.0 && gains.c == 1.0 && gains.Tf == 0.0);
        static_assert(gains.integral_lim.upper(0) == std::numeric_limits<double>::max());
        CHECK(pi.u_lim.lower(0) == -12.0f);
        CHECK(pi.u_lim.upper(0) == 12.0f);
        CHECK(pi.control(100.0f, 0.0f) == doctest::Approx(12.0f));
    }

    TEST_CASE("PID aliases deduce the scalar and take their mode from the alias") {
        constexpr double              Ts = 0.01;
        constexpr auto                gains = design::pid(4.0, 8.0, 0.0).discretize(Ts);
        static constinit PIController pi = gains.as<float>(); // scalar deduced, mode from the alias
        static_assert(std::is_same_v<decltype(pi), PIDController<float, PIDMode::PI>>);
        static_assert(std::is_same_v<decltype(PIDController{gains}), PIDController<double, PIDMode::PID>>);
        static_assert(std::is_same_v<decltype(PController{gains}), PIDController<double, PIDMode::P>>);
        static_assert(std::is_same_v<decltype(PDController{gains.as<float>()}), PIDController<float, PIDMode::PD>>);
        CHECK(pi.Kp == 4.0f);
    }

    TEST_CASE("angle input: r′ = y + wrap_pi(r − y) on an unwrapped y takes the short way across ±π") {
        constexpr double       Ts = 0.01;
        constexpr double       turn = 2.0 * damp::numbers::pi_v<double>;
        constexpr double       rate = 0.01; // [rad/tick] sensor slope (1 rad/s)
        constexpr double       lead = 0.05; // [rad] reference ahead of the sensor
        const auto             gains = design::PIDGains{.Kp = 2.0, .Ki = 5.0, .Kd = 0.1, .c = 0.0, .Tf = 0.02}.discretize(Ts);
        PIDController<double>  pid{gains};
        PIDController<double>  wrapped_y{gains}; // fed the sensor angle itself
        PhaseUnwrapper<double> unwrap;

        double d_max = 0.0;
        double wrapped_d_max = 0.0;
        for (int k = 0; k < 40; ++k) {
            const double theta = wrap_pi(3.0 + (rate * static_cast<double>(k))); // sensor wraps +π → −π at k = 15
            const double r = wrap_pi(theta + lead);                              // reference wraps first, at k = 10
            const double y = turn * unwrap.update(theta / turn).value();
            const double u = pid.control(y + wrap_pi(r - y), y);
            CHECK(u - pid.integral - pid.deriv == doctest::Approx(gains.Kp * lead)); // P on the short-way error
            d_max = damp::max(d_max, damp::abs(pid.deriv));

            (void)wrapped_y.control(theta + wrap_pi(r - theta), theta);
            wrapped_d_max = damp::max(wrapped_d_max, damp::abs(wrapped_y.deriv));
        }
        CHECK(pid.integral == doctest::Approx(40.0 * gains.Ki * lead));          // I took +lead every tick
        CHECK(d_max <= doctest::Approx(gains.Kd * rate / Ts));                   // D ≤ Kd·ẏ: no 2π step
        CHECK(wrapped_d_max == doctest::Approx(gains.d_b * turn).epsilon(0.01)); // a wrapped y steps D by d_b·2π
    }
} // TEST_SUITE
