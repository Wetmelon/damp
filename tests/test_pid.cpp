// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <limits>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/design/status.hpp"
#include "damp/encoder/unwrap.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "damp/toolbox/bounds.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

/**
 * @file test_pid.cpp
 * @brief Runtime tests for PIDController: the tick at T_s and at a measured dt, track(u, r, y)
 *        handover, gains and limits at run time, back-calculation, Law<kRes, T>.
 */

namespace {

/// |a − b| within a few ulps of @p scale: two inlined copies of a formula may round apart (FMA contraction).
template<typename T>
bool within_ulps(T a, T b, T scale, T ulps = T{8}) {
    return damp::abs(a - b) <= ulps * eps<T>() * (scale + damp::abs(a) + damp::abs(b));
}

} // namespace

TEST_SUITE("PID handover") {
    TEST_CASE("PI: track then step returns u, clamped into the input limits") {
        constexpr float     Ts = 0.01f;
        PIController<float> pi{design::pid(1.0f, 1.0f, 0.0f).discretize(Ts, {.input = {.lower = -5.0f, .upper = 5.0f}})};

        pi.track(2.5f, 1.0f, 0.0f);
        CHECK(pi.last_command() == 2.5f); // recorded as the command applied
        CHECK(pi.step(1.0f, 0.0f) == doctest::Approx(2.5f));

        pi.track(100.0f, 1.0f, 0.0f);
        CHECK(pi.step(1.0f, 0.0f) == doctest::Approx(5.0f));

        pi.track(-100.0f, 1.0f, 0.0f);
        CHECK(pi.step(1.0f, 0.0f) == doctest::Approx(-5.0f));
        //! The preload is toward u clamped: the integral did not wind up to −100.
        CHECK(pi.integral() == doctest::Approx(-5.0f - 1.0f));
    }

    TEST_CASE("PI: tracking for a while, then stepping, is bumpless") {
        constexpr float     Ts = 0.01f;
        PIController<float> pi{design::pid(2.0f, 5.0f, 0.0f).discretize(Ts)};

        for (int i = 0; i < 20; ++i) {
            (void)pi.step(1.0f, 0.5f);
        }
        REQUIRE(pi.integral() != doctest::Approx(0.0f));

        const float u_track = 3.7f;
        for (int i = 0; i < 10; ++i) {
            pi.track(u_track, 1.0f, 0.5f);
        }
        CHECK(pi.step(1.0f, 0.5f) == doctest::Approx(u_track).epsilon(1e-5f));
    }

    TEST_CASE("PID: bumpless transfer with active derivative term") {
        constexpr float      Ts = 0.01f;
        PIDController<float> pid{design::pid(1.5f, 3.0f, 0.05f).discretize(Ts)};

        for (int i = 0; i < 20; ++i) {
            const float y = 0.5f + (0.01f * static_cast<float>(i));
            (void)pid.step(1.0f, y);
        }

        const float u_track = -1.2f;
        for (int i = 0; i < 20; ++i) {
            pid.track(u_track, 1.0f, 0.7f);
        }
        CHECK(pid.step(1.0f, 0.7f) == doctest::Approx(u_track).epsilon(1e-5f));
    }

    TEST_CASE("PID: bumpless with a filtered derivative still decaying at hand-off, at T_s and at a measured dt") {
        constexpr double Ts = 0.01;
        auto             gains = design::pid(1.5, 3.0, 0.05);
        gains.Tf = 0.05; // d_a = Tf/(Ts+Tf) ≈ 0.83: D carries over into the first step
        PIDController<double> pid{gains.discretize(Ts)};
        PIDController<double> measured{gains.discretize(Ts)};

        double y = 0.0;
        for (int i = 0; i < 20; ++i) {
            y = 0.5 + (0.02 * static_cast<double>(i)); // moving y keeps D nonzero
            pid.track(-1.2, 1.0, y);
            measured.track(-1.2, 1.0, y, Ts);
        }
        REQUIRE(pid.derivative() != doctest::Approx(0.0));

        // Same (r, y) as the last track: the first step emits u_track.
        CHECK(pid.step(1.0, y) == doctest::Approx(-1.2).epsilon(1e-12));
        CHECK(measured.step(1.0, y, Ts) == doctest::Approx(-1.2).epsilon(1e-12));
    }

    TEST_CASE_TEMPLATE("track(u, r, y) then step(r, y) returns u with b ≠ 1 and r ≠ y", T, float, double) {
        //! The old track(u, y) preloaded as if r = y: the first step added Kp·(b·r − y) − Kp·(b − 1)·y
        //! to u, 1.15 here, so the handover was exact only where b·r = y.
        const T                    Ts = static_cast<T>(0.01);
        const design::PIDGains<T>  g{.Kp = T{2}, .Ki = T{8}, .Kd = static_cast<T>(0.05), .b = static_cast<T>(0.3), .c = T{0}, .Tf = static_cast<T>(0.02)};
        const design::PIDResult<T> res = g.discretize(Ts, {.input = {.lower = T{-10}, .upper = T{10}}});
        PIDController<T>           pid{res};
        PIController<T>            pi{design::PIDGains<T>{.Kp = T{2}, .Ki = T{8}, .b = static_cast<T>(0.3)}.discretize(Ts)};
        const T                    r = T{2};
        const T                    u = static_cast<T>(1.75);
        T                          y = T{0};
        for (int k = 0; k < 30; ++k) {
            y = static_cast<T>(0.5) + (static_cast<T>(0.01) * static_cast<T>(k));
            (void)pid.step(r, y);
            (void)pi.step(r, y);
        }
        pid.track(u, r, y);
        pi.track(u, r, y);
        const T scale = damp::abs(g.Kp * ((g.b * r) - y)) + damp::abs(pid.integral()) + damp::abs(pid.derivative());
        CHECK(within_ulps(pid.step(r, y), u, scale));
        CHECK(within_ulps(pi.step(r, y), u, scale));

        //! At a measured dt the same: exact for a next tick of the same dt.
        PIDController<T> measured{res};
        const T          dt = static_cast<T>(0.013);
        (void)measured.step(r, T{0}, dt);
        measured.track(u, r, y, dt);
        CHECK(within_ulps(measured.step(r, y, dt), u, scale));
    }

    TEST_CASE("PD: P + filtered D; track has nothing to preload and keeps the derivative current") {
        constexpr double Ts = 0.01;
        auto             gains = design::pid(2.0, 0.0, 0.1);
        gains.Tf = 0.02;
        PDController<double> pd{gains.discretize(Ts)};
        const double         d_a = 0.02 / (Ts + 0.02);
        const double         d_b = 0.1 / (Ts + 0.02);

        CHECK(pd.step(1.0, 0.0) == doctest::Approx(2.0)); // first tick: no derivative kick
        const double u = pd.step(1.0, 0.1);               // y moved 0.1: D = d_b·(−0.1)
        CHECK(u == doctest::Approx((2.0 * 0.9) - (d_b * 0.1)));

        //! The old Tracking mode held u_track until enable(); now the step after a track is the
        //! law's own command, its derivative decayed through the track at the same y.
        pd.track(-3.0, 1.0, 0.1);
        CHECK(pd.last_command() == -3.0);
        CHECK(pd.step(1.0, 0.1) == doctest::Approx((2.0 * 0.9) - (d_a * d_a * d_b * 0.1)));

        //! A first track seeds the derivative memory: no kick on the step after it.
        PDController<double> fresh{gains.discretize(Ts)};
        fresh.track(0.0, 1.0, 0.4);
        CHECK(fresh.step(1.0, 0.4) == doctest::Approx(2.0 * 0.6));
    }

    TEST_CASE("PID: the track preload respects the integral limits") {
        constexpr float      Ts = 0.01f;
        PIDController<float> pid{design::pid(1.0f, 1.0f, 0.0f).discretize(Ts, {.integral = {.lower = -2.0f, .upper = 2.0f}})};

        pid.track(100.0f, 0.0f, 0.0f);
        CHECK(pid.integral() == doctest::Approx(2.0f));

        pid.track(-100.0f, 0.0f, 0.0f);
        CHECK(pid.integral() == doctest::Approx(-2.0f));
    }

    TEST_CASE("P: track has nothing to preload and records u; the next step is the plain law") {
        PController<float> p{design::PIDGains{.Kp = 2.0f}.discretize(0.01f, {.input = {.lower = -10.0f, .upper = 10.0f}})};

        p.track(4.2f, 0.0f, 0.0f);
        CHECK(p.last_command() == 4.2f);
        CHECK(p.step(1.0f, 0.0f) == doctest::Approx(2.0f));
        CHECK(p.step(100.0f, 0.0f) == doctest::Approx(10.0f));
    }

    TEST_CASE("reset clears the integral, the derivative and the last command, and keeps gains and limits") {
        constexpr float      Ts = 0.01f;
        PIDController<float> pid{design::pid(1.0f, 1.0f, 0.1f).discretize(Ts, {.input = {.lower = -1.0f, .upper = 1.0f}})};
        for (int k = 0; k < 5; ++k) {
            (void)pid.step(1.0f, 0.1f * static_cast<float>(k));
        }
        pid.reset();
        CHECK(pid.integral() == 0.0f);
        CHECK(pid.derivative() == 0.0f);
        CHECK(pid.last_command() == 0.0f);
        CHECK(pid.gains().Kd == 0.1f);
        CHECK(pid.limits().input.upper(0) == 1.0f);
        CHECK(pid.step(0.5f, 0.0f) == doctest::Approx(0.5f + (Ts * 0.5f))); // no derivative kick after reset
    }

    TEST_CASE("back_calculate and track compose: track, then step, then saturate") {
        constexpr float     Ts = 0.01f;
        PIController<float> pi{design::PIDGains{.Kp = 1.0f, .Ki = 1.0f, .tracking_time = 0.5f}.discretize(Ts, {.input = {.lower = -1.0f, .upper = 1.0f}})};

        pi.track(0.5f, 2.0f, 0.0f);
        CHECK(pi.step(2.0f, 0.0f) == doctest::Approx(0.5f).epsilon(1e-5f));

        for (int i = 0; i < 50; ++i) {
            (void)pi.step(5.0f, 0.0f);
        }
        CHECK(pi.step(5.0f, 0.0f) == doctest::Approx(1.0f));
    }
} // TEST_SUITE

TEST_SUITE("PID tick") {
    TEST_CASE("step(r, y, dt): a dt that is not positive holds the state and returns the last command") {
        PIDController<double> pid{design::pid(2.0, 5.0, 0.1).discretize(0.01)};
        const double          u1 = pid.step(1.0, 0.0, 0.01);
        const double          integ = pid.integral();

        const double u0 = pid.step(1.0, 0.0, 0.0);
        CHECK(pid.integral() == doctest::Approx(integ));
        CHECK(u0 == u1);
        // Integral state is the term (already includes Ki history)
        CHECK(u0 == doctest::Approx((2.0 * 1.0) + integ));

        (void)pid.step(1.0, 0.0, -0.5);
        CHECK(pid.integral() == doctest::Approx(integ));

        PIController<float> pi{design::pid(1.0f, 1.0f, 0.0f).discretize(0.01f)};
        (void)pi.step(1.0f, 0.0f, 0.01f);
        const float pinteg = pi.integral();
        (void)pi.step(1.0f, 0.0f, 0.0f);
        CHECK(pi.integral() == doctest::Approx(pinteg));
    }

    TEST_CASE("step(r, y, dt): back_calculate after a held tick is a no-op; a NaN dt holds") {
        // Kp = 1, Ki = 10, |u| ≤ 1, T_t = 0.1: a negative measured period used to charge
        // back_calculate with Ki·Ts < 0, winding the integral further into the limit.
        PIController<double> pid{design::PIDGains<double>{.Kp = 1.0, .Ki = 10.0, .tracking_time = 0.1}.discretize(0.01, {.input = {.lower = -1.0, .upper = 1.0}})};
        for (int k = 0; k < 5; ++k) {
            (void)pid.step(0.5, 0.0, 0.01);
        }
        const double integ = pid.integral();
        (void)pid.step(0.5, 0.0, -0.01);
        pid.back_calculate(pid.last_command() - 0.5);
        CHECK(pid.integral() == integ);

        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double u = pid.step(0.5, 0.0, nan);
        CHECK(std::isfinite(u));
        CHECK(std::isfinite(pid.integral()));
        CHECK(std::isfinite(pid.derivative()));
        CHECK(std::isfinite(pid.step(0.5, 0.0, 0.01)));
    }

    TEST_CASE_TEMPLATE("step(r, y, dt) at dt = T_s matches step(r, y) to a few ulps over a run", T, float, double) {
        //! The fixed tick reads coefficients computed at set_gains; the measured one computes them
        //! from 1/(dt + T_f): the same numbers up to the rounding of that reciprocal.
        const T                   Ts = static_cast<T>(0.01);
        const design::PIDGains<T> g{.Kp = T{2}, .Ki = T{5}, .Kd = static_cast<T>(0.1), .tracking_time = static_cast<T>(0.4), .c = T{0}, .Tf = static_cast<T>(0.05)};
        const auto                res = g.discretize(Ts, {.input = {.lower = T{-10}, .upper = T{10}}, .integral = {.lower = T{-100}, .upper = T{100}}});
        PIDController<T>          fixed{res};
        PIDController<T>          measured{res};
        int                       clamped = 0;
        for (int k = 0; k < 200; ++k) {
            const T r = (k < 100) ? T{8} : T{-8};
            const T y = static_cast<T>(0.01) * static_cast<T>(k % 37);
            const T ud = fixed.step(r, y);
            const T uc = measured.step(r, y, Ts);
            clamped += (damp::abs(ud) >= T{10}) ? 1 : 0;
            CHECK(within_ulps(ud, uc, damp::abs(fixed.integral()) + damp::abs(fixed.derivative()) + T{1}, T{64}));
        }
        CHECK(clamped > 0);
    }

    TEST_CASE("PID first tick has no derivative kick; Tf attenuates the derivative") {
        constexpr double Ts = 0.01;

        PIDController<double> d{design::pid(0.0, 0.0, 1.0).discretize(Ts)};
        CHECK(d.step(0.0, 1.0) == doctest::Approx(0.0));    // seeded: no kick
        CHECK(d.step(0.0, 2.0) == doctest::Approx(-100.0)); // Kd*(-1)/Ts

        auto step = [&](double Tf) {
            PIDController<double> c{
                design::PIDGains{.Kd = 1.0, .Tf = Tf}.discretize(Ts)
            };
            (void)c.step(0.0, 0.0);
            return c.step(0.0, 1.0);
        };
        CHECK(std::abs(step(0.1)) < std::abs(step(0.0)));
    }

    TEST_CASE("discretize carries the gains, the limits and T_s; the runtime scales Ki by T_s") {
        constexpr double Ts = 0.01;
        const auto       disc = design::PIDGains{.Kp = 1.0, .Ki = 20.0}.discretize(Ts, {.input = {.lower = -1.0, .upper = 1.0}, .integral = {.lower = -0.5, .upper = 0.5}});

        CHECK(disc.success);
        CHECK(disc.gains.Kp == 1.0);
        CHECK(disc.gains.Ki == 20.0); // continuous: the runtime takes Ki·Ts
        CHECK(disc.limits.integral.lower(0) == doctest::Approx(-0.5));
        CHECK(disc.limits.integral.upper(0) == doctest::Approx(0.5));
        CHECK(disc.Ts == doctest::Approx(Ts));
        PIController<double> pi{disc};
        (void)pi.step(1.0, 0.9);
        CHECK(pi.integral() == doctest::Approx(20.0 * Ts * 0.1));
        CHECK(pi.sample_time() == Ts);
    }

    TEST_CASE("discretize keeps default unbounded limits for constinit") {
        constexpr double Ts = 0.01;
        constexpr auto   disc = design::pid(4.0, 8.0, 0.0).discretize(Ts);
        constexpr double lim = std::numeric_limits<double>::max();
        static_assert(disc.limits.integral.lower(0) == -lim);
        static_assert(disc.limits.input.upper(0) == lim);
        // Full design→deploy path is a constant expression.
        constexpr auto ctrl = PIController<float>{design::pid(4.0, 8.0, 0.0).discretize(Ts).as<float>()};
        static_assert(ctrl.gains().Kp == 4.0f);
        static_assert(ctrl.gains().Ki == 8.0f); // was the discrete Ki·Ts: gains() returns the continuous gains
        static_assert(ctrl.sample_time() == 0.01f);
    }

    TEST_CASE("PIDGains designated initializers deduce T and deploy as static constinit") {
        constexpr double                     Ts = 0.01;
        static constinit PIController<float> pi = design::PIDGains{.Kp = 4.0, .Ki = 8.0}.discretize(Ts, {.input = {.lower = -12.0, .upper = 12.0}}).as<float>();
        static_assert(std::is_same_v<decltype(design::PIDGains{.Kp = 4.0}), design::PIDGains<double>>);
        constexpr design::PIDGains gains{.Kp = 4.0, .Ki = 8.0};
        static_assert(gains.Kd == 0.0 && gains.b == 1.0 && gains.c == 1.0 && gains.Tf == 0.0 && gains.tracking_time == 0.0);
        static_assert(design::PIDLimits<double>{}.integral.upper(0) == std::numeric_limits<double>::max());
        CHECK(pi.limits().input.lower(0) == -12.0f);
        CHECK(pi.limits().input.upper(0) == 12.0f);
        CHECK(pi.step(100.0f, 0.0f) == doctest::Approx(12.0f));
    }

    TEST_CASE("PID aliases deduce the scalar and take their mode from the alias") {
        constexpr double              Ts = 0.01;
        constexpr auto                gains = design::pid(4.0, 8.0, 0.0).discretize(Ts);
        static constinit PIController pi = gains.as<float>(); // scalar deduced, mode from the alias
        static_assert(std::is_same_v<decltype(pi), PIDController<float, PIDMode::PI>>);
        static_assert(std::is_same_v<decltype(PIDController{gains}), PIDController<double, PIDMode::PID>>);
        static_assert(std::is_same_v<decltype(PController{gains}), PIDController<double, PIDMode::P>>);
        static_assert(std::is_same_v<decltype(PDController{gains.as<float>()}), PIDController<float, PIDMode::PD>>);
        //! One runtime written once: the aliases name its modes.
        static_assert(std::is_same_v<PIController<float>, PIDController<float, PIDMode::PI>>);
        CHECK(pi.gains().Kp == 4.0f);
    }

    TEST_CASE("cross-precision converting ctor copies gains, limits and live state") {
        constexpr double      Ts = 0.01;
        PIDController<double> d{design::pid(1.5, 2.0, 0.3).discretize(Ts, {.input = {.lower = -2.0, .upper = 2.0}})};
        d.track(0.4, 1.0, 0.2);
        (void)d.step(1.0, 0.25);

        PIDController<float> f(d);
        CHECK(f.valid());
        CHECK(f.gains().Kp == doctest::Approx(static_cast<float>(d.gains().Kp)));
        CHECK(f.gains().Ki == doctest::Approx(static_cast<float>(d.gains().Ki)));
        CHECK(f.gains().Kd == doctest::Approx(static_cast<float>(d.gains().Kd)));
        CHECK(f.limits().input.upper(0) == 2.0f);
        CHECK(f.sample_time() == doctest::Approx(0.01f));
        CHECK(f.integral() == doctest::Approx(static_cast<float>(d.integral())));
        CHECK(f.derivative() == doctest::Approx(static_cast<float>(d.derivative())));
        CHECK(f.last_command() == doctest::Approx(static_cast<float>(d.last_command())));
        CHECK(f.step(1.0f, 0.25f) == doctest::Approx(static_cast<float>(d.step(1.0, 0.25))).epsilon(1e-5));
    }

    TEST_CASE("angle input: r′ = y + wrap_pi(r − y) on an unwrapped y takes the short way across ±π") {
        constexpr double       Ts = 0.01;
        constexpr double       turn = 2.0 * damp::numbers::pi_v<double>;
        constexpr double       rate = 0.01; // [rad/tick] sensor slope (1 rad/s)
        constexpr double       lead = 0.05; // [rad] reference ahead of the sensor
        const auto             res = design::PIDGains{.Kp = 2.0, .Ki = 5.0, .Kd = 0.1, .c = 0.0, .Tf = 0.02}.discretize(Ts);
        const double           d_b = res.gains.Kd / (Ts + res.gains.Tf);
        PIDController<double>  pid{res};
        PIDController<double>  wrapped_y{res}; // fed the sensor angle itself
        PhaseUnwrapper<double> unwrap;

        double d_max = 0.0;
        double wrapped_d_max = 0.0;
        for (int k = 0; k < 40; ++k) {
            const double theta = wrap_pi(3.0 + (rate * static_cast<double>(k))); // sensor wraps +π → −π at k = 15
            const double r = wrap_pi(theta + lead);                              // reference wraps first, at k = 10
            const double y = turn * unwrap.update(theta / turn).value();
            const double u = pid.step(y + wrap_pi(r - y), y);
            CHECK(u - pid.integral() - pid.derivative() == doctest::Approx(res.gains.Kp * lead)); // P on the short-way error
            d_max = damp::max(d_max, damp::abs(pid.derivative()));

            (void)wrapped_y.step(theta + wrap_pi(r - theta), theta);
            wrapped_d_max = damp::max(wrapped_d_max, damp::abs(wrapped_y.derivative()));
        }
        CHECK(pid.integral() == doctest::Approx(40.0 * res.gains.Ki * Ts * lead)); // I took +lead every tick
        CHECK(d_max <= doctest::Approx(res.gains.Kd * rate / Ts));                 // D ≤ Kd·ẏ: no 2π step
        CHECK(wrapped_d_max == doctest::Approx(d_b * turn).epsilon(0.01));         // a wrapped y steps D by d_b·2π
    }
} // TEST_SUITE

TEST_SUITE("PID gains and limits at run time") {
    TEST_CASE("get-modify-set: a new Kp, Ki or Kd takes effect at the next tick") {
        constexpr double      Ts = 0.01;
        const auto            res = design::PIDGains{.Kp = 1.0, .Ki = 4.0, .Kd = 0.05, .Tf = 0.02}.discretize(Ts);
        PIDController<double> law{res};
        double                y = 0.0;
        for (int k = 0; k < 50; ++k) {
            y = 0.01 * static_cast<double>(k);
            (void)law.step(1.0, y);
        }

        //! Kp: the next command moves by ΔKp·(b·r − y), nothing else.
        PIDController<double> kp = law;
        auto                  g = kp.gains();
        g.Kp = 2.5;
        REQUIRE(kp.set_gains(g));
        CHECK(kp.gains().Kp == 2.5);
        const double          u_kp = kp.step(1.0, y);
        PIDController<double> same = law;
        const double          u_same = same.step(1.0, y);
        CHECK(u_kp - u_same == doctest::Approx(1.5 * (1.0 - y)).epsilon(1e-12));

        //! Kd: the next derivative takes the new gain on Δ(c·r − y).
        PIDController<double> kd = law;
        g = kd.gains();
        g.Kd = 0.1;
        REQUIRE(kd.set_gains(g));
        PIDController<double> kd_ref = law;
        (void)kd.step(1.0, y + 0.01);
        (void)kd_ref.step(1.0, y + 0.01);
        CHECK(kd.derivative() - kd_ref.derivative() == doctest::Approx((0.05 / (Ts + 0.02)) * -0.01).epsilon(1e-9));
    }

    TEST_CASE("a new Ki does not bump the command: the integral is held as the I term") {
        constexpr double     Ts = 0.01;
        PIController<double> pi{design::pid(0.0, 10.0, 0.0).discretize(Ts)};
        for (int k = 0; k < 50; ++k) {
            (void)pi.step(1.0, 0.0); // e = 1 each tick
        }
        const double I_before = pi.integral();
        REQUIRE(I_before != doctest::Approx(0.0));
        const double u_before = pi.step(0.0, 0.0); // e = 0: u = I (Kp = 0)
        auto         g = pi.gains();
        g.Ki *= 2.0;
        REQUIRE(pi.set_gains(g));
        const double u_after = pi.step(0.0, 0.0);
        CHECK(pi.integral() == I_before);
        CHECK(u_after == u_before);
        CHECK(u_after == I_before);
        //! It changes the rate: the next error integrates at the new Ki.
        (void)pi.step(1.0, 0.0);
        CHECK(pi.integral() - I_before == doctest::Approx(20.0 * Ts).epsilon(1e-9));
    }

    TEST_CASE("set_limits moves the rails at the next tick and pulls the integral into its new limits") {
        PIController<float> pi{design::pid(1.0f, 10.0f, 0.0f).discretize(0.01f)};
        for (int k = 0; k < 100; ++k) {
            (void)pi.step(1.0f, 0.5f);
        }
        REQUIRE(pi.integral() > 2.0f);
        auto lim = pi.limits();
        lim.input = {.lower = -1.0f, .upper = 1.0f};
        lim.integral = {.lower = -0.75f, .upper = 0.75f};
        REQUIRE(pi.set_limits(lim));
        CHECK(pi.integral() == 0.75f); // at once
        CHECK(pi.step(1.0f, 0.5f) == 1.0f);
        CHECK(pi.limits().input.lower(0) == -1.0f);
    }

    TEST_CASE("set_gains refuses a negative Tf, a negative tracking time and a gain the mode does not compute, and changes nothing") {
        constexpr double     Ts = 0.01;
        PIController<double> pi{design::PIDGains{.Kp = 1.0, .Ki = 2.0, .tracking_time = 0.5}.discretize(Ts)};
        (void)pi.step(1.0, 0.0);
        const auto   before = pi.gains();
        const double I = pi.integral();

        auto neg_tf = before;
        neg_tf.Tf = -0.001;
        CHECK_FALSE(pi.set_gains(neg_tf));
        auto neg_tt = before;
        neg_tt.tracking_time = -0.1;
        CHECK_FALSE(pi.set_gains(neg_tt));
        auto with_d = before;
        with_d.Kd = 0.1; // a PI cannot realize the D term
        CHECK_FALSE(pi.set_gains(with_d));
        CHECK_FALSE(pi.set_gains(design::detail::refused<design::PIDGains<double>>(design::DesignStatus::not_placeable)));
        CHECK(pi.gains().Kp == before.Kp);
        CHECK(pi.gains().Ki == before.Ki);
        CHECK(pi.gains().Tf == before.Tf);
        CHECK(pi.gains().tracking_time == before.tracking_time);
        CHECK(pi.integral() == I);

        PDController<double> pd{design::pid(1.0, 0.0, 0.1).discretize(Ts)};
        CHECK_FALSE(pd.set_gains(design::pid(1.0, 0.5, 0.1))); // a PD cannot realize the I term
        CHECK(pd.gains().Ki == 0.0);

        auto inverted = pi.limits();
        inverted.input = {.lower = 1.0, .upper = -1.0};
        CHECK_FALSE(pi.set_limits(inverted));
        CHECK(pi.limits().input.lower(0) == -std::numeric_limits<double>::max());
    }

    TEST_CASE("a default-constructed law ticks at a measured dt with gains and limits set at run time") {
        //! A hand-written C PID's shape: no design, no T_s.
        PIController<float> pi;
        REQUIRE(pi.set_gains({.Kp = 2.0f, .Ki = 10.0f}));
        REQUIRE(pi.set_limits({.input = {.lower = -1.0f, .upper = 1.0f}}));
        CHECK(pi.sample_time() == 0.0f);
        CHECK(pi.step(0.25f, 0.0f, 0.01f) == doctest::Approx(0.5f + 0.025f));
        CHECK(pi.step(5.0f, 0.0f, 0.01f) == 1.0f);
        //! Its step(r, y) has no period: the integral holds.
        const float I = pi.integral();
        (void)pi.step(0.1f, 0.0f);
        CHECK(pi.integral() == I);
    }
} // TEST_SUITE

/// What happens to the integral at a limit, with no anti-windup parameter set: it holds.
TEST_SUITE("PID limits") {
    constexpr double kTs = 0.01;

    TEST_CASE_TEMPLATE("an output limit alone keeps the integral from winding up", T, float, double) {
        //! With Kbc = 0 the integral ran on behind the clamp: Kp 1, Ki 10, ±1, ten seconds on the
        //! rail left I = 500 and the output stayed on the rail for 99.6 s after the error reversed.
        PIController<T> pi{design::PIDGains<double>{.Kp = 1.0, .Ki = 10.0}.discretize(kTs, {.input = {.lower = -1.0, .upper = 1.0}}).template as<T>()};
        for (int k = 0; k < 1000; ++k) {
            CHECK(pi.step(T{5}, T{0}) == T{1});
        }
        CHECK(pi.integral() == T{0});         // P alone is past the limit: nothing of I is needed
        CHECK(pi.step(T{0.5}, T{0}) < T{1});  // off the rail on the first tick the error allows it
        CHECK(pi.step(T{-5}, T{0}) == T{-1}); // and onto the other
        CHECK(pi.integral() >= T{0});         // without the integral having gone the wrong way
    }

    TEST_CASE_TEMPLATE("the integral stops where the output meets the limit", T, float, double) {
        PIController<T> pi{design::PIDGains<double>{.Kp = 1.0, .Ki = 10.0}.discretize(kTs, {.input = {.lower = -1.0, .upper = 1.0}}).template as<T>()};
        T               u{0};
        for (int k = 0; k < 200; ++k) {
            u = pi.step(T{0.5}, T{0}); // P = 0.5, I ramps 0.05 a tick towards 0.5
        }
        CHECK(u == T{1});
        CHECK(damp::abs(pi.integral() - T{0.5}) <= T{8} * eps<T>());
        //! A step of the P term through the limit does not reset what the integral holds.
        (void)pi.step(T{3}, T{0});
        CHECK(damp::abs(pi.integral() - T{0.5}) <= T{8} * eps<T>());
        //! Nor does the integral limit stop applying.
        auto lim = pi.limits();
        lim.integral = {.lower = T{-0.25}, .upper = T{0.25}};
        REQUIRE(pi.set_limits(lim));
        (void)pi.step(static_cast<T>(0.1), T{0});
        CHECK(pi.integral() == T{0.25});
    }

    TEST_CASE_TEMPLATE("an input limit alone is the law clamped outside and told with back_calculate", T, float, double) {
        //! The oracle of the LQ limits: one tick, whichever stage clamps.
        const auto       gains = design::PIDGains<double>{.Kp = 2.0, .Ki = 8.0, .Kd = 0.02, .c = 0.0, .Tf = 0.005};
        const auto       bounded_res = gains.discretize(kTs, {.input = {.lower = -1.0, .upper = 1.0}}).template as<T>();
        PIDController<T> bounded{bounded_res};
        PIDController<T> outside{gains.discretize(kTs).template as<T>()};
        PIDController<T> measured{bounded_res};
        T                ya{0};
        T                yb{0};
        T                yc{0};
        int              saturated = 0;
        for (int k = 0; k < 600; ++k) {
            const T r = (k < 300) ? T{1} : T{-0.5};
            const T ua = bounded.step(r, ya);
            const T ub = damp::clamp(outside.step(r, yb), T{-1}, T{1});
            outside.back_calculate(ub);
            const T uc = measured.step(r, yc, static_cast<T>(kTs));
            saturated += (damp::abs(ua) >= T{1}) ? 1 : 0;
            CHECK(damp::abs(ua - ub) <= T{64} * eps<T>());
            CHECK(damp::abs(ua - uc) <= T{64} * eps<T>());
            ya += static_cast<T>(kTs) * (ua - (T{0.5} * ya)); // a lag with a gain of 2: r = 1 needs u = 0.5
            yb += static_cast<T>(kTs) * (ub - (T{0.5} * yb));
            yc += static_cast<T>(kTs) * (uc - (T{0.5} * yc));
        }
        CHECK(saturated > 10);
        CHECK(static_cast<double>(ya) == doctest::Approx(-0.5).epsilon(0.02));
    }

    TEST_CASE("back-calculation takes at most the whole excess: a tracking time under T_s does not diverge") {
        //! ΔI = Ki·Ts·(u − u_unsat)/Kbc with Kbc = Ki·Ts/10 (T_t = Ts/10) overshot the unsaturated
        //! value tenfold: the integral changed sign and grew ninefold each tick, to ±3.4e38.
        PIController<float> pi{design::PIDGains<double>{.Kp = 1.0, .Ki = 10.0, .tracking_time = 0.001}.discretize(kTs, {.input = {.lower = -1.0, .upper = 1.0}}).as<float>()};
        for (int k = 0; k < 200; ++k) {
            CHECK(pi.step(5.0F, 0.0F) == 1.0F);
            CHECK(damp::abs(pi.integral()) <= 5.0F); // I tracks the limit less P: −4, each tick
        }
    }

    TEST_CASE("back_calculate takes back the last increment with T_t = 0, bleeds with T_t > 0") {
        PIDController<double> pid{design::pid(1.0, 1.0, 0.0).discretize(kTs)};
        (void)pid.step(1.0, 0.0);
        CHECK(pid.integral() == doctest::Approx(0.01));
        pid.back_calculate(0.5); // T_t = 0: the clamp cut 0.51 off, more than the tick added
        CHECK(pid.integral() == 0.0);
        pid.back_calculate(0.5); // nothing left of that tick to give back
        CHECK(pid.integral() == 0.0);

        PIDController<double> pid2{design::PIDGains{.Kp = 1.0, .Ki = 1.0, .tracking_time = 0.5}.discretize(kTs)};
        (void)pid2.step(1.0, 0.0);
        const double before = pid2.integral();
        pid2.back_calculate(0.5); // u_sat below the last command → integral down
        CHECK(pid2.integral() < before);
    }

    TEST_CASE_TEMPLATE("back_calculate(u_sat, u_ff) is ΔI = min(h/T_t, 1)·(u_sat − (u + u_ff)), h the tick's period", T, float, double) {
        // A PI plus feedforward, clamped downstream, at T_s and at a measured dt.
        const T                    Ts = static_cast<T>(0.01);
        const design::PIDResult<T> res = design::PIDGains<T>{.Kp = T{2}, .Ki = T{20}, .tracking_time = static_cast<T>(0.1)}.discretize(Ts);
        PIController<T>            law{res};
        PIController<T>            measured{res};
        const T                    u_ff = static_cast<T>(0.3);
        const T                    u_lim = T{1};
        const T                    dt = static_cast<T>(0.004);
        T                          y = T{0};
        int                        clamped = 0;
        for (int k = 0; k < 200; ++k) {
            const T r = (k < 100) ? T{2} : T{-2};
            const T u = law.step(r, y);
            const T u_sat = damp::clamp(u + u_ff, -u_lim, u_lim);
            const T um = measured.step(r, y, dt);
            const T um_sat = damp::clamp(um + u_ff, -u_lim, u_lim);
            if (u_sat != u + u_ff) {
                ++clamped;
                const T I = law.integral();
                law.back_calculate(u_sat, u_ff);
                CHECK(within_ulps(law.integral() - I, (Ts / static_cast<T>(0.1)) * (u_sat - (u + u_ff)), u_lim));
            }
            if (um_sat != um + u_ff) {
                const T I = measured.integral();
                measured.back_calculate(um_sat, u_ff);
                CHECK(within_ulps(measured.integral() - I, (dt / static_cast<T>(0.1)) * (um_sat - (um + u_ff)), u_lim));
            }
            y += Ts * (u_sat - y); // first-order plant
        }
        CHECK(clamped > 0);
    }

    TEST_CASE_TEMPLATE("tracking_time reproduces the old Kbc on a saturating run: T_t = Kbc/Ki", T, float, double) {
        //! The old PI integrated I += Ki·Ts·e and unwound min(Ki·Ts/Kbc, 1)·(u − u_unsat): the same
        //! numbers, to rounding, as min(Ts/T_t, 1) with T_t = Kbc/Ki.
        const T Ts = static_cast<T>(0.01);
        const T Kp = T{2};
        const T Ki = T{25};
        const T Kbc = static_cast<T>(1.5);
        struct {
            T Ki_d, Kp, Kbc, I{0};
            T tick(T r, T y) {
                I += Ki_d * (r - y);
                const T u_unsat = (Kp * (r - y)) + I;
                const T u = damp::clamp(u_unsat, T{-1}, T{1});
                if (u != u_unsat) {
                    I += damp::min(Ki_d / Kbc, T{1}) * (u - u_unsat);
                }
                return u;
            }
        } old{Ki * Ts, Kp, Kbc};
        PIController<T> pi{design::PIDGains<T>{.Kp = Kp, .Ki = Ki, .tracking_time = Kbc / Ki}.discretize(Ts, {.input = {.lower = T{-1}, .upper = T{1}}})};
        T               y = T{0};
        int             clamped = 0;
        for (int k = 0; k < 400; ++k) {
            const T r = (k % 200 < 100) ? T{3} : T{-3};
            const T u_old = old.tick(r, y);
            const T u = pi.step(r, y);
            clamped += (damp::abs(u) >= T{1}) ? 1 : 0;
            CHECK(within_ulps(u, u_old, damp::abs(old.I) + T{1}, T{64}));
            CHECK(within_ulps(pi.integral(), old.I, T{1}, T{64}));
            y += Ts * (u - y);
        }
        CHECK(clamped > 50);
    }

    TEST_CASE("a runtime built at run time from a result it refuses is invalid and sends 0, in every mode") {
        // They trapped (__builtin_trap, a fault on the target) where every other law is invalid.
        const Bounds<1>           inverted{.lower = 1.0, .upper = -1.0};
        design::PIDResult<double> g{.gains = {.Kp = 1.0}, .limits = {.input = inverted}, .Ts = 0.01};
        design::detail::set_status(g, design::DesignStatus::ok); // hand-built: discretize would refuse it
        PIDController<double> pid{g};
        PIController<double>  pi{g};
        PDController<double>  pd{g};
        PController<double>   p{g};
        auto                  with_i = design::pid(1.0, 0.1, 0.0).discretize(0.01);
        PDController<double>  pd_of_pi{with_i};                                             // a PD cannot realize I
        PController<double>   failed{design::PIDGains<double>{.Kp = 1.0}.discretize(-1.0)}; // a failed design
        auto                  no_ts = design::pid(1.0, 0.0, 0.0).discretize(0.01);
        no_ts.Ts = 0.0;
        PController<double> zero_ts{no_ts}; // no period to tick at
        CHECK_FALSE(pid.valid());
        CHECK_FALSE(pi.valid());
        CHECK_FALSE(pd.valid());
        CHECK_FALSE(p.valid());
        CHECK_FALSE(pd_of_pi.valid());
        CHECK_FALSE(failed.valid());
        CHECK_FALSE(zero_ts.valid());
        CHECK(pid.step(1.0, 0.0) == 0.0);
        CHECK(pi.step(1.0, 0.0) == 0.0);
        CHECK(pd.step(1.0, 0.0) == 0.0);
        CHECK(p.step(1.0, 0.0) == 0.0);
        CHECK(pid.step(1.0, 0.0, 0.01) == 0.0);
        CHECK(failed.step(1.0, 0.0) == 0.0);
        CHECK(PController<double>{design::pid(1.0, 0.0, 0.0).discretize(0.01)}.valid());
        CHECK(PIDController<double>{}.valid()); // a hand-built law: zero gains, valid
    }

    TEST_CASE("a limit that is not a box is refused, as the LQ and MPC designs refuse it") {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const auto   g = design::PIDGains<double>{.Kp = 1.0, .Ki = 1.0};
        CHECK(g.discretize(kTs, {.input = {.lower = 1.0, .upper = -1.0}}).status == design::DesignStatus::invalid_parameter);
        CHECK(g.discretize(kTs, {.input = {.lower = nan, .upper = 1.0}}).status == design::DesignStatus::invalid_parameter);
        CHECK(g.discretize(kTs, {.integral = {.lower = 0.5, .upper = -0.5}}).status == design::DesignStatus::invalid_parameter);
        CHECK(g.discretize(kTs, {.input = {.lower = 0.0, .upper = 0.0}}).success); // pinned, not refused
    }

    TEST_CASE("discretize refuses a negative Tf, a negative tracking time and non-finite gains") {
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

        // T_t < 0 would push the integral further into the limit; it replaces the old Kbc against Ki's sign.
        constexpr auto bad = design::PIDGains{.Kp = 1.0, .Ki = 2.0, .tracking_time = -1.0}.discretize(0.01);
        static_assert(!bad.success);
        CHECK(bad.status == design::DesignStatus::invalid_parameter);
        CHECK_FALSE(design::PIDGains<double>{.Kp = 1.0, .Ki = 2.0, .tracking_time = inf}.discretize(Ts).success);
        // A Kp < 0 law, a Ki < 0 one and Ki = 0 all take a positive tracking time.
        static_assert(design::PIDGains{.Kp = -1.0, .Ki = 2.0, .tracking_time = 0.5}.discretize(0.01).success);
        static_assert(design::PIDGains{.Kp = 1.0, .Ki = -2.0, .tracking_time = 0.5}.discretize(0.01).success);
        static_assert(design::PIDGains{.Kp = 1.0, .Ki = 0.0, .tracking_time = 0.5}.discretize(0.01).success);
        static constinit PIController<double> ok{design::PIDGains{.Kp = -1.0, .Ki = 2.0, .tracking_time = 0.5}.discretize(0.01)};
        CHECK(ok.gains().tracking_time == 0.5);
    }
}

namespace {
inline constexpr auto kPI = design::pid(2.0, 5.0, 0.0).discretize(0.01, {.input = {.lower = -1.0, .upper = 1.0}});
inline constexpr auto kPID = design::PIDGains{.Kp = 2.0, .Ki = 5.0, .Kd = 0.1, .Tf = 0.02}.discretize(0.01);
inline constexpr auto kPD = design::pid(2.0, 0.0, 0.1).discretize(0.01);
inline constexpr auto kP = design::pid(2.0, 0.0, 0.0).discretize(0.01);
inline constexpr auto kFailed = design::pid(2.0, 5.0, 0.0).discretize(-1.0);
const auto            kAtStart = design::pid(2.0, 5.0, 0.0).discretize(0.01); // initialized at start-up

struct Axis {
    Law<kPI, float> law; // a struct member's type, named by its design
    float           offset{0.0f};
};
} // namespace

TEST_SUITE("PID Law") {
    TEST_CASE("Law<kRes, float> deploys a constexpr design in the narrowest mode its gains need") {
        static_assert(std::is_same_v<Law<kPI, float>::Runtime, PIController<float>>);
        static_assert(std::is_same_v<Law<kPID, float>::Runtime, PIDController<float>>);
        static_assert(std::is_same_v<Law<kPD, float>::Runtime, PDController<float>>);
        static_assert(std::is_same_v<Law<kP, float>::Runtime, PController<float>>);
        static_assert(std::is_base_of_v<PIController<float>, Law<kPI, float>>);
        static_assert(std::is_same_v<Law<kAtStart, float>::Runtime, PIDController<float>>); // not constexpr: the full law
        //! A failed design does not compile wherever its Law is declared: the static_assert in the
        //! Law is constant_design_runs, which refuses it.
        static_assert(!detail::constant_design<kAtStart>);
        static_assert(!kFailed.success && kFailed.status == design::DesignStatus::invalid_sample_time);
        // static constinit Law<kFailed, float> refused; // ill-formed: design_refused::invalid_sample_time()

        static constinit Axis axis;
        PIController<float>   ref{kPI.as<float>()};
        for (int k = 0; k < 20; ++k) {
            CHECK(axis.law.step(1.0f, 0.1f) == ref.step(1.0f, 0.1f));
        }
        CHECK(axis.law.limits().input.upper(0) == 1.0f);
        //! Its gains and limits stay settable within its mode.
        auto g = axis.law.gains();
        g.Kp = 3.0f;
        CHECK(axis.law.set_gains(g));
        g.Kd = 0.1f;
        CHECK_FALSE(axis.law.set_gains(g)); // a PI's Law has no D term
        Law<kAtStart, float> held;
        CHECK(held.valid());
        CHECK(held.step(1.0f, 0.0f) == doctest::Approx(2.0f + 0.05f));
    }
}
