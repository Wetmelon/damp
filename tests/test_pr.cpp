// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <numbers>

#include "damp/controllers/pr.hpp"
#include "damp/filters/filters.hpp"
#include "damp/math/math.hpp"
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

        constexpr auto result = design::pr({.Kp = 1.0, .Ki = 100.0, .w0 = w0, .wc = wc}, Ts);
        static_assert(result.success);
        static_assert(result.Kp == 1.0);
        static_assert(result.Ki == 100.0);
        CHECK(result.w0 == doctest::Approx(w0));
        CHECK(result.wc == doctest::Approx(wc));
    }

    TEST_CASE("PR rejects invalid parameters") {
        static_assert(!design::pr({.Kp = 1.0, .Ki = 100.0, .w0 = 0.0, .wc = 10.0}, 1e-4).success);    // w0 <= 0
        static_assert(!design::pr({.Kp = 1.0, .Ki = 100.0, .w0 = 314.0, .wc = -1.0}, 1e-4).success);  // wc < 0
        static_assert(!design::pr({.Kp = 1.0, .Ki = 100.0, .w0 = 314.0, .wc = 10.0}, -1e-4).success); // Ts < 0
        static_assert(design::pr({.Kp = 1.0, .Ki = 100.0, .w0 = 314.0, .wc = 0.0}, 1e-4).success);    // ideal PR
        static_assert(design::pr({.Kp = 1.0, .Ki = 100.0, .w0 = 314.0, .wc = 10.0}, 0.0).success);    // continuous

        PRController<double> bad(design::pr({.Kp = 1.0, .Ki = 100.0, .w0 = 0.0, .wc = 10.0}, 1e-4));
        CHECK_FALSE(bad.valid());
        CHECK(bad.control(1.0) == 0.0);
    }

    TEST_CASE("Ideal PR (wc=0) has non-zero resonant term") {
        constexpr double w0 = 2.0 * std::numbers::pi * 50.0;
        constexpr double Ts = 1.0 / 10000.0;
        constexpr auto   r = design::pr({.Kp = 1.0, .Ki = 200.0, .w0 = w0, .wc = 0.0}, Ts);
        static_assert(r.success);

        // Continuous ideal law: num = {Kp, Ki, Kp*w0²}, den = {1, 0, w0²}
        auto tf = r.to_tf();
        CHECK(tf.num[0] == doctest::Approx(1.0));
        CHECK(tf.num[1] == doctest::Approx(200.0));
        CHECK(tf.num[2] == doctest::Approx(1.0 * w0 * w0));
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

        auto                 result = design::pr({.Kp = 1.0, .Ki = 200.0, .w0 = w0, .wc = wc}, Ts);
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
        constexpr auto      result = design::pr({.Kp = 1.0, .Ki = 100.0, .w0 = 314.0, .wc = 10.0}, 0.0001);
        PRController<float> ctrl(result.as<float>());
        CHECK(ctrl.Kp == doctest::Approx(1.0f));
        CHECK(ctrl.w0 == doctest::Approx(314.0f).epsilon(0.1));
    }

    TEST_CASE("PR controller reset") {
        auto                 result = design::pr({.Kp = 1.0, .Ki = 100.0, .w0 = 314.0, .wc = 10.0}, 0.0001);
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
        auto                 result = design::pr({.Kp = 1.0, .Ki = 100.0, .w0 = 314.0, .wc = 10.0}, 0.0001);
        PRController<double> ctrl(result);

        // Change to 60 Hz
        double w0_60 = 2.0 * std::numbers::pi * 60.0;
        ctrl.set_frequency(w0_60);
        CHECK(ctrl.w0 == doctest::Approx(w0_60));
    }

    TEST_CASE("PR control(r,y) matches control(error)") {
        auto                 result = design::pr({.Kp = 2.0, .Ki = 100.0, .w0 = 314.0, .wc = 10.0}, 0.0001);
        PRController<double> a(result);
        PRController<double> b(result);
        for (int k = 0; k < 50; ++k) {
            double r = std::sin(0.1 * k);
            double y = 0.3 * r;
            CHECK(a.control(r, y) == doctest::Approx(b.control(r - y)));
        }
    }

    TEST_CASE("PR converting ctor preserves Tt and state") {
        PRController<double> d(design::pr({.Kp = 1.0, .Ki = 100.0, .w0 = 314.0, .wc = 10.0, .Tt = 5.0}, 0.0001));
        (void)d.control(0.7);
        (void)d.control(-0.2);

        PRController<float> f(d); // U -> T converting ctor
        CHECK(f.Tt == doctest::Approx(5.0f));
        // Resonant delay line carried over.
        CHECK(f.resonant.last_output() == doctest::Approx(static_cast<float>(d.resonant.last_output())));
        // Coefficients carried over too: the next step matches within float precision.
        CHECK(static_cast<double>(f.control(0.1f)) == doctest::Approx(d.control(0.1)).epsilon(0.01));
    }

    TEST_CASE("PR back_calculate unwinds resonant state") {
        PRController<double> ctrl(design::pr({.Kp = 1.0, .Ki = 100.0, .w0 = 314.0, .wc = 10.0, .Tt = 2.0}, 0.0001));
        const double         u = ctrl.control(1.0);
        double               y_before = ctrl.resonant.last_output();

        // No clamp: the applied command is the law's own plus the feedforward -> unchanged.
        ctrl.back_calculate(u + 0.5, 0.5);
        CHECK(ctrl.resonant.last_output() == y_before);

        // Saturation: realizable command 2 below the requested one -> pull output down by Ts·2/Tt.
        ctrl.back_calculate(u - 2.0);
        CHECK(ctrl.resonant.last_output() == doctest::Approx(y_before - (0.0001 * 2.0 / 2.0)));

        // No-op when Tt == 0 (back-calculation not configured), even with excess.
        PRController<double> nctrl(design::pr({.Kp = 1.0, .Ki = 100.0, .w0 = 314.0, .wc = 10.0}, 0.0001));
        const double         un = nctrl.control(1.0);
        double               y_noop = nctrl.resonant.last_output();
        nctrl.back_calculate(un - 2.0); // Tt == 0 -> unchanged
        CHECK(nctrl.resonant.last_output() == doctest::Approx(y_noop));
    }

    TEST_CASE("A continuous (Ts = 0) design leaves PRController invalid") {
        // There is no Tustin map at Ts = 0: the runtime used to drop the resonator and emit Kp·e.
        PRController<double> ctrl(design::pr({.Kp = 2.0, .Ki = 100.0, .w0 = 314.0, .wc = 10.0}, 0.0));
        CHECK_FALSE(ctrl.valid());
        CHECK(ctrl.control(1.0) == 0.0);
        CHECK(ctrl.control(-3.0) == 0.0);
    }

    TEST_CASE("PR refuses a resonance at or above Nyquist; set_frequency ignores one") {
        constexpr double Ts = 1e-4; // fs = 10 kHz, Nyquist 5 kHz
        const auto       above = design::pr({.Kp = 1.0, .Ki = 100.0, .w0 = 2.0 * std::numbers::pi * 6000.0, .wc = 10.0}, Ts);
        CHECK_FALSE(above.success);
        CHECK(above.status == design::DesignStatus::invalid_parameter);
        CHECK_FALSE(design::pr({.Kp = 1.0, .Ki = 100.0, .w0 = std::numbers::pi / Ts, .wc = 10.0}, Ts).success);
        CHECK(design::pr({.Kp = 1.0, .Ki = 100.0, .w0 = 2.0 * std::numbers::pi * 4000.0, .wc = 10.0}, Ts).success);

        PRController<double> ctrl(design::pr({.Kp = 1.0, .Ki = 100.0, .w0 = 314.0, .wc = 10.0}, Ts));
        ctrl.set_frequency(2.0 * std::numbers::pi * 6000.0);
        CHECK(ctrl.w0 == 314.0);

        const double nan = std::numeric_limits<double>::quiet_NaN();
        CHECK(design::pr({.Kp = 1.0, .Ki = 100.0, .w0 = 314.0, .wc = 10.0}, nan).status == design::DesignStatus::invalid_sample_time);
        CHECK(design::pr({.Kp = 1.0, .Ki = 100.0, .w0 = 314.0, .wc = 10.0, .Tt = -1.0}, Ts).status == design::DesignStatus::invalid_parameter);
    }

    TEST_CASE("PRResult to_tf / to_ss match the PR law at DC and resonance") {
        constexpr double w0 = 314.0;
        constexpr double wc = 10.0;
        constexpr double Kp = 2.0;
        constexpr double Ki = 100.0;
        auto             tf = design::pr({.Kp = Kp, .Ki = Ki, .w0 = w0, .wc = wc}, 0.0).to_tf();
        auto             ss = design::pr({.Kp = Kp, .Ki = Ki, .w0 = w0, .wc = wc}, 0.0).to_ss();

        // DC gain (s=0): num[2]/den[2] = Kp*w0²/w0² = Kp
        CHECK(tf.num[2] / tf.den[2] == doctest::Approx(Kp));
        // Resonance (s = jw0), num(s) = num[0] s² + num[1] s + num[2]: |C(jw0)| = Kp + Ki
        const double nre = tf.num[2] - (tf.num[0] * w0 * w0);
        const double dre = tf.den[2] - (tf.den[0] * w0 * w0);
        const double nim = tf.num[1] * w0;
        const double dim = tf.den[1] * w0;
        CHECK(damp::sqrt(((nre * nre) + (nim * nim)) / ((dre * dre) + (dim * dim))) == doctest::Approx(Kp + Ki));
        // State-space D feedthrough is the proportional term.
        CHECK(ss.D(0, 0) == doctest::Approx(Kp));
        CHECK(ss.C(0, 1) == doctest::Approx(2.0 * Ki * wc));
    }

    TEST_CASE("PRResult to_discrete_ss resonates on the 11th and 13th harmonics at 5 kHz") {
        // Frequency response C(zI − A)⁻¹B + D of a 2-state discrete model at f [Hz].
        const auto response = [](const StateSpace<2, 1, 1, double>& ss, double f) {
            const std::complex<double> z = std::polar(1.0, 2.0 * std::numbers::pi * f * ss.Ts);
            const std::complex<double> m00 = z - ss.A(0, 0);
            const std::complex<double> m01 = -ss.A(0, 1);
            const std::complex<double> m10 = -ss.A(1, 0);
            const std::complex<double> m11 = z - ss.A(1, 1);
            const std::complex<double> det = (m00 * m11) - (m01 * m10);
            const std::complex<double> x0 = ((m11 * ss.B(0, 0)) - (m01 * ss.B(1, 0))) / det;
            const std::complex<double> x1 = ((m00 * ss.B(1, 0)) - (m10 * ss.B(0, 0))) / det;
            return (ss.C(0, 0) * x0) + (ss.C(0, 1) * x1) + ss.D(0, 0);
        };

        constexpr double                w_fund = 2.0 * std::numbers::pi * 50.0;
        constexpr double                Ts = 1.0 / 5000.0;
        constexpr std::array<size_t, 3> harmonics = {1, 11, 13};
        constexpr auto                  results = design::pr_harmonics({.Kp = 1.0, .Ki = 100.0, .w0 = w_fund, .wc = 10.0}, Ts, harmonics);

        for (size_t i = 1; i < 3; ++i) {
            const auto   dss = results[i].to_discrete_ss();
            const double f0 = 50.0 * static_cast<double>(harmonics[i]); // 550 Hz, 650 Hz
            CAPTURE(f0);
            CHECK(dss.Ts == Ts);

            // Kp = 0 on the harmonics: the resonator alone, real and equal to Ki on ω₀.
            const std::complex<double> h0 = response(dss, f0);
            CHECK(h0.real() == doctest::Approx(results[i].Ki).epsilon(1e-9));
            CHECK(std::abs(h0.imag()) < 1e-9 * results[i].Ki);

            // Peak of |H| over f0 ± 30 Hz on a 0.01 Hz grid sits on f0 (plain Tustin: 529.6 / 617 Hz).
            double f_peak = 0.0;
            double h_peak = 0.0;
            for (int k = -3000; k <= 3000; ++k) {
                const double f = f0 + (0.01 * static_cast<double>(k));
                const double h = std::abs(response(dss, f));
                if (h > h_peak) {
                    h_peak = h;
                    f_peak = f;
                }
            }
            CHECK(f_peak == doctest::Approx(f0).epsilon(1e-5));
        }
    }

    TEST_CASE("PRResult to_discrete_ss matches PRController tick by tick") {
        constexpr double w0 = 2.0 * std::numbers::pi * 550.0;
        constexpr double Ts = 1.0 / 5000.0;

        for (const double wc : {10.0, 0.0}) { // non-ideal and ideal PR
            CAPTURE(wc);
            const auto r = design::pr({.Kp = 1.5, .Ki = 80.0, .w0 = w0, .wc = wc}, Ts);
            const auto dss = r.to_discrete_ss();

            for (const bool impulse : {false, true}) {
                CAPTURE(impulse);
                PRController<double> ctrl(r);
                Matrix<2, 1, double> x = Matrix<2, 1, double>::zeros();
                for (int k = 0; k < 200; ++k) {
                    const double e = (impulse && k > 0) ? 0.0 : 1.0;
                    const double u_model = (dss.C * x)(0, 0) + (dss.D(0, 0) * e);
                    const double u_runtime = ctrl.control(e);
                    CHECK(u_model == doctest::Approx(u_runtime).epsilon(1e-9).scale(r.Ki));
                    x = (dss.A * x) + (dss.B * Matrix<1, 1, double>{{e}});
                }
            }
        }
    }

    TEST_CASE("Multi-harmonic PR design") {
        constexpr double                w_fund = 2.0 * std::numbers::pi * 50.0;
        constexpr double                wc = 10.0;
        constexpr double                Ts = 1.0 / 10000.0;
        constexpr std::array<size_t, 4> harmonics = {1, 3, 5, 7};

        constexpr auto results = design::pr_harmonics({.Kp = 1.0, .Ki = 100.0, .w0 = w_fund, .wc = wc}, Ts, harmonics);

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
