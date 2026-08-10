// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>

#include "damp/controllers/pr.hpp"
#include "damp/filters/filters.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

/**
 * @brief Tests for Proportional-Resonant controller (pr.hpp)
 *
 * Validates PR design, discretization, and harmonic tracking.
 */

TEST_SUITE("PR Controller") {
    TEST_CASE("PR design result") {
        constexpr double w0 = 2.0 * std::numbers::pi * 50.0; // 50 Hz
        constexpr double wc = 10.0;
        constexpr double Ts = 1.0 / 10000.0; // 10 kHz

        constexpr auto result = design::pr(1.0, 100.0, w0, wc, Ts);
        static_assert(result.success);
        static_assert(result.Kp == 1.0);
        static_assert(result.Ki == 100.0);
        CHECK(result.w0 == doctest::Approx(w0));
        CHECK(result.wc == doctest::Approx(wc));
    }

    TEST_CASE("PR rejects invalid parameters") {
        static_assert(!design::pr(1.0, 100.0, 0.0, 10.0, 1e-4).success);    // w0 <= 0
        static_assert(!design::pr(1.0, 100.0, 314.0, -1.0, 1e-4).success);  // wc < 0
        static_assert(!design::pr(1.0, 100.0, 314.0, 10.0, -1e-4).success); // Ts < 0
        static_assert(design::pr(1.0, 100.0, 314.0, 0.0, 1e-4).success);    // ideal PR
        static_assert(design::pr(1.0, 100.0, 314.0, 10.0, 0.0).success);    // continuous

        PRController<double> bad(design::pr(1.0, 100.0, 0.0, 10.0, 1e-4));
        CHECK_FALSE(bad.valid());
        CHECK(bad.control(1.0) == 0.0);
    }

    TEST_CASE("Ideal PR (wc=0) has non-zero resonant term") {
        constexpr double w0 = 2.0 * std::numbers::pi * 50.0;
        constexpr double Ts = 1.0 / 10000.0;
        constexpr auto   r = design::pr(1.0, 200.0, w0, 0.0, Ts);
        static_assert(r.success);

        // Continuous ideal law: num = {Kp*w0², Ki, Kp}, den = {w0², 0, 1}
        auto tf = r.to_tf();
        CHECK(tf.num[0] == doctest::Approx(1.0 * w0 * w0));
        CHECK(tf.num[1] == doctest::Approx(200.0));
        CHECK(tf.num[2] == doctest::Approx(1.0));
        CHECK(tf.den[1] == doctest::Approx(0.0));

        auto ss = r.to_ss();
        CHECK(ss.C(0, 1) == doctest::Approx(200.0)); // Ki, not 2*Ki*wc=0
        CHECK(ss.A(1, 1) == doctest::Approx(0.0));   // no damping

        // Runtime: first sample of ideal resonator is non-zero (not pure-P).
        PRController<double> ctrl(r);
        CHECK(ctrl.valid());
        double u = ctrl.control(1.0);
        // Pure-P would give u = Kp * 1 = 1; ideal resonator adds a b0 contribution.
        CHECK(damp::abs(u - 1.0) > 1e-6);
        CHECK(damp::isfinite(u));
    }

    TEST_CASE("PR controller tracks sinusoidal reference") {
        // PR at 50 Hz should achieve zero steady-state error for 50 Hz sine
        double w0 = 2.0 * std::numbers::pi * 50.0;
        double wc = 10.0;
        double Ts = 1.0 / 10000.0;

        auto                 result = design::pr(1.0, 200.0, w0, wc, Ts);
        PRController<double> ctrl(result);

        // Simulate: reference = sin(w0*t), plant = unit gain (identity)
        double y = 0.0;
        double error_sum = 0.0;
        size_t n_steps = 20000; // 2 seconds at 10 kHz

        for (size_t k = 0; k < n_steps; ++k) {
            double t = static_cast<double>(k) * Ts;
            double ref = std::sin(w0 * t);
            double error = ref - y;
            double u = ctrl.control(error);
            y = u; // Unit gain plant

            // Track error in last 0.5 seconds (settled)
            if (k > 15000) {
                error_sum += std::abs(error);
            }
        }

        double avg_error = error_sum / 5000.0;
        // Should have very small tracking error at steady state
        CHECK(avg_error < 0.05);
    }

    TEST_CASE("PR controller type conversion") {
        constexpr auto      result = design::pr(1.0, 100.0, 314.0, 10.0, 0.0001);
        PRController<float> ctrl(result.as<float>());
        CHECK(ctrl.Kp == doctest::Approx(1.0f));
        CHECK(ctrl.w0 == doctest::Approx(314.0f).epsilon(0.1));
    }

    TEST_CASE("PR controller reset") {
        auto                 result = design::pr(1.0, 100.0, 314.0, 10.0, 0.0001);
        PRController<double> ctrl(result);

        // Run some steps
        (void)ctrl.control(1.0);
        (void)ctrl.control(0.5);
        (void)ctrl.control(0.3);

        ctrl.reset();
        // After reset, output should be just Kp * error (no resonant history)
        double u = ctrl.control(1.0);
        // First sample after reset: resonant contribution is b0_r * 1.0
        // The Kp contribution is 1.0
        CHECK(u != 0.0);
    }

    TEST_CASE("PR frequency update") {
        auto                 result = design::pr(1.0, 100.0, 314.0, 10.0, 0.0001);
        PRController<double> ctrl(result);

        // Change to 60 Hz
        double w0_60 = 2.0 * std::numbers::pi * 60.0;
        ctrl.set_frequency(w0_60);
        CHECK(ctrl.w0 == doctest::Approx(w0_60));
    }

    TEST_CASE("PR control(r,y) matches control(error)") {
        auto                 result = design::pr(2.0, 100.0, 314.0, 10.0, 0.0001);
        PRController<double> a(result);
        PRController<double> b(result);
        for (int k = 0; k < 50; ++k) {
            double r = std::sin(0.1 * k);
            double y = 0.3 * r;
            CHECK(a.control(r, y) == doctest::Approx(b.control(r - y)));
        }
    }

    TEST_CASE("PR converting ctor preserves Kbc and state") {
        PRController<double> d(design::pr(1.0, 100.0, 314.0, 10.0, 0.0001));
        d.Kbc = 5.0;
        (void)d.control(0.7);
        (void)d.control(-0.2);

        PRController<float> f(d); // U -> T converting ctor
        CHECK(f.Kbc == doctest::Approx(5.0f));
        // Resonant delay line carried over.
        CHECK(f.resonant.last_output() == doctest::Approx(static_cast<float>(d.resonant.last_output())));
        // Coefficients carried over too: the next step matches within float precision.
        CHECK(static_cast<double>(f.control(0.1f)) == doctest::Approx(d.control(0.1)).epsilon(0.01));
    }

    TEST_CASE("PR back_calculate unwinds resonant state") {
        PRController<double> ctrl(design::pr(1.0, 100.0, 314.0, 10.0, 0.0001));
        ctrl.Kbc = 2.0;
        (void)ctrl.control(1.0);
        double y_before = ctrl.resonant.last_output();

        // Saturation: realizable command below the requested one -> pull output down.
        ctrl.back_calculate(5.0, 3.0, 0.0001);
        CHECK(ctrl.resonant.last_output() < y_before);

        // No-op when Kbc == 0 (back-calculation not configured), even with excess.
        PRController<double> nctrl(design::pr(1.0, 100.0, 314.0, 10.0, 0.0001));
        (void)nctrl.control(1.0);
        double y_noop = nctrl.resonant.last_output();
        nctrl.back_calculate(5.0, 3.0, 0.0001); // Kbc == 0 -> unchanged
        CHECK(nctrl.resonant.last_output() == doctest::Approx(y_noop));
    }

    TEST_CASE("PR Ts<=0 yields proportional-only output") {
        PRController<double> ctrl(design::pr(2.0, 100.0, 314.0, 10.0, 0.0));
        CHECK(ctrl.control(1.0) == doctest::Approx(2.0)); // Kp*error, no resonant
        CHECK(ctrl.control(-3.0) == doctest::Approx(-6.0));
    }

    TEST_CASE("PRResult to_tf / to_ss match the PR law at DC and resonance") {
        constexpr double w0 = 314.0;
        constexpr double wc = 10.0;
        constexpr double Kp = 2.0;
        constexpr double Ki = 100.0;
        auto             tf = design::pr(Kp, Ki, w0, wc, 0.0).to_tf();
        auto             ss = design::pr(Kp, Ki, w0, wc, 0.0).to_ss();

        // DC gain (s=0): num[0]/den[0] = Kp*w0²/w0² = Kp
        CHECK(tf.num[0] / tf.den[0] == doctest::Approx(Kp));
        // State-space D feedthrough is the proportional term.
        CHECK(ss.D(0, 0) == doctest::Approx(Kp));
        CHECK(ss.C(0, 1) == doctest::Approx(2.0 * Ki * wc));
    }

    TEST_CASE("MultiPRController runtime: control and reset") {
        constexpr double                w_fund = 2.0 * std::numbers::pi * 50.0;
        constexpr std::array<size_t, 3> harmonics = {1, 3, 5};
        auto                            results = design::pr_harmonics(1.0, 100.0, w_fund, 10.0, 1.0 / 10000.0, harmonics);

        MultiPRController<3, double> ctrl(results);
        double                       u = ctrl.control(1.0); // first sample: Kp*1 + sum of resonant b0
        CHECK(u != 0.0);

        // control(r,y) forwards the error to control(error).
        MultiPRController<3, double> a(results);
        MultiPRController<3, double> b(results);
        CHECK(a.control(0.5, 0.1) == doctest::Approx(b.control(0.4)));

        // reset clears all resonant history -> first post-reset sample is finite/nonzero.
        ctrl.reset();
        CHECK(ctrl.control(1.0) != 0.0);
    }

    TEST_CASE("MultiPRController converting ctor preserves state") {
        constexpr double                w_fund = 2.0 * std::numbers::pi * 50.0;
        constexpr std::array<size_t, 3> harmonics = {1, 3, 5};
        auto                            results = design::pr_harmonics(1.0, 100.0, w_fund, 10.0, 1.0 / 10000.0, harmonics);

        MultiPRController<3, double> d(results);
        (void)d.control(0.5);

        MultiPRController<3, float> f(d); // U -> T converting ctor
        CHECK(f.Kp == doctest::Approx(1.0f));
        CHECK(f.resonants[0].resonant.last_output() == doctest::Approx(static_cast<float>(d.resonants[0].resonant.last_output())));
    }

    TEST_CASE("Multi-harmonic PR design") {
        constexpr double                w_fund = 2.0 * std::numbers::pi * 50.0;
        constexpr double                wc = 10.0;
        constexpr double                Ts = 1.0 / 10000.0;
        constexpr std::array<size_t, 4> harmonics = {1, 3, 5, 7};

        constexpr auto results = design::pr_harmonics(1.0, 100.0, w_fund, wc, Ts, harmonics);

        // First harmonic should have Kp and Ki_fund
        CHECK(results[0].Kp == doctest::Approx(1.0));
        CHECK(results[0].Ki == doctest::Approx(100.0));
        CHECK(results[0].w0 == doctest::Approx(w_fund));

        // Third harmonic: Ki/3, w0*3
        CHECK(results[1].Kp == doctest::Approx(0.0));
        CHECK(results[1].Ki == doctest::Approx(100.0 / 3.0));
        CHECK(results[1].w0 == doctest::Approx(3.0 * w_fund));

        // 7th harmonic
        CHECK(results[3].w0 == doctest::Approx(7.0 * w_fund));
    }
}
