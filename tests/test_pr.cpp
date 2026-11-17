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

#include "damp/backend.hpp"
#include "damp/controllers/pr.hpp"
#include "damp/design/status.hpp"
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

namespace {

/** @brief sin and cos of ω·k·Ts in double, one rotation per tick. */
struct SineSource {
    double s{0.0};
    double c{1.0};
    double sin_step{};
    double cos_step{};

    SineSource(double w, double Ts) : sin_step(damp::sin(w * Ts)), cos_step(damp::cos(w * Ts)) {}

    /** @brief Advance one tick, pulled back onto the unit circle so the amplitude does not drift over millions of ticks. */
    void step() {
        const double sn = (s * cos_step) + (c * sin_step);
        const double cn = (c * cos_step) - (s * sin_step);
        const double pull = 1.5 - (0.5 * ((sn * sn) + (cn * cn)));
        s = sn * pull;
        c = cn * pull;
    }
};

/** @brief A steady-state response u ≈ gain·amp·sin(ωt + phase), and what the run did on the way. */
struct SineResponse {
    double gain{};
    double phase{};
    double peak{};       ///< Largest |u| of the whole run, the settling included
    bool   finite{true}; ///< No output was NaN or ±Inf
};

/**
 * @brief Steady-state response of @p ctrl to amp·sin(ω·k·Ts)
 *
 * Runs settle + fit ticks and fits u ≈ a·sin + b·cos over the last fit by least squares, in
 * double, so the fit needs no whole number of periods.
 */
template<typename T>
SineResponse sine_response(PRController<T>& ctrl, double w, double Ts, double amp, size_t settle, size_t fit) {
    SineSource   src(w, Ts);
    SineResponse out{};
    double       ss = 0.0;
    double       sc = 0.0;
    double       cc = 0.0;
    double       us = 0.0;
    double       uc = 0.0;
    for (size_t k = 0; k < settle + fit; ++k) {
        const T u = ctrl.control(static_cast<T>(amp * src.s));
        out.finite = out.finite && damp::isfinite(u);
        out.peak = damp::max(out.peak, damp::abs(static_cast<double>(u)));
        if (k >= settle) {
            ss += src.s * src.s;
            sc += src.s * src.c;
            cc += src.c * src.c;
            us += static_cast<double>(u) * src.s;
            uc += static_cast<double>(u) * src.c;
        }
        src.step();
    }
    const double det = (ss * cc) - (sc * sc);
    const double a = ((us * cc) - (uc * sc)) / det;
    const double b = ((uc * ss) - (us * sc)) / det;
    out.gain = damp::hypot(a, b) / amp;
    out.phase = damp::atan2(b, a);
    return out;
}

/**
 * @brief R(jω) of the pre-warped Tustin design, in closed form
 *
 * With t = tan(ωTs/2), β = tan(ω₀Ts/2) and α = ω_c·Ts/2: R = 2·Kᵢ·α·jt / (β² − t² + 2α·jt),
 * which is Kᵢ at ω = ω₀.
 */
std::complex<double> pr_resonant_response(double Ki, double w0, double wc, double Ts, double w) {
    const double t = std::tan(w * Ts / 2.0);
    const double beta = std::tan(w0 * Ts / 2.0);
    const double alpha = wc * Ts / 2.0;
    return std::complex<double>{0.0, 2.0 * Ki * alpha * t} / std::complex<double>{(beta * beta) - (t * t), 2.0 * alpha * t};
}

/**
 * @brief What a resonator in T holds of its design, relative to the response (in rad for a phase)
 *
 * ω₀ is placed to a few ulps and the peak is ω_c wide, so an ulp of ω₀ moves the response at ω₀
 * by ε·Q, Q = ω₀/(2ω_c). The roundoff of one tick, ε of the state, adds up over the 2/c₂ ticks the
 * resonator remembers, c₂ = 2ω_c·Ts: ε/√c₂. @p rate is ω_c, or the decay rate of the loop the
 * resonator is closed in.
 */
template<typename T>
double resonator_tol(double w0, double rate, double Ts) {
    return static_cast<double>(eps<T>()) * ((8.0 * (1.0 + (w0 / (2.0 * rate)))) + (4.0 / damp::sqrt(2.0 * rate * Ts)));
}

/** @brief Ticks after which a transient decaying at @p rate [rad/s] is ε^1.25 of its start: 20 time constants in float, 45 in double. */
template<typename T>
size_t settle_ticks(double rate, double Ts) {
    return static_cast<size_t>(damp::ceil(-1.25 * damp::log(static_cast<double>(eps<T>())) / (rate * Ts)));
}

} // namespace

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

        // Tt == 0: the resonant section holds. It gives back what the tick added (from rest, its
        // whole output), up to the excess, and a second report of the clamp takes nothing more.
        PRController<double> nctrl(design::pr({.Kp = 1.0, .Ki = 100.0, .w0 = 314.0, .wc = 10.0}, 0.0001));
        const double         un = nctrl.control(1.0);
        REQUIRE(nctrl.resonant.last_output() > 0.0);
        nctrl.back_calculate(un - 2.0);
        CHECK(nctrl.resonant.last_output() == 0.0); // was unchanged: no back-calculation configured
        nctrl.back_calculate(un - 2.0);
        CHECK(nctrl.resonant.last_output() == 0.0);
    }

    TEST_CASE_TEMPLATE("u_lim: the resonant section holds against the clamp, as against one reported", T, float, double) {
        const double Ts = 1e-4;
        const double w0 = 2.0 * std::numbers::pi * 50.0;
        auto         res = design::pr({.Kp = 0.0, .Ki = 200.0, .w0 = w0, .wc = 5.0}, Ts);
        res.u_lim = {.lower = -0.1, .upper = 0.1};
        PRController<T> inside{res.template as<T>()};
        auto            open = res;
        open.u_lim = {};
        PRController<T> outside{open.template as<T>()};
        const T         lim = static_cast<T>(0.1);
        T               peak{0};
        for (int k = 0; k < 4000; ++k) { // 20 periods of an error at ω₀, whose gain there is Ki = 200
            const T e = static_cast<T>(std::sin(w0 * Ts * k));
            const T u = inside.control(e);
            const T sent = std::clamp(outside.control(e), -lim, lim);
            outside.back_calculate(sent);
            CHECK(damp::abs(u - sent) <= T{64} * eps<T>() * lim);
            peak = damp::max(peak, damp::abs(inside.resonant.last_output()));
        }
        CHECK(peak <= lim * static_cast<T>(1.001)); // was 200: the resonance wound up behind the clamp
    }

    TEST_CASE("a tracking time under the sample period unwinds at most the whole excess") {
        // ΔR = (Ts/Tt)·excess overshot the realizable command past Tt < Ts and diverged past Tt < Ts/2.
        const double         Ts = 1e-4;
        const double         w0 = 2.0 * std::numbers::pi * 50.0;
        PRController<double> pr{design::pr({.Kp = 0.0, .Ki = 200.0, .w0 = w0, .wc = 5.0, .Tt = Ts / 4.0}, Ts)};
        double               peak = 0.0;
        for (int k = 0; k < 4000; ++k) {
            const double u = pr.control(std::sin(w0 * Ts * k));
            pr.back_calculate(std::clamp(u, -0.1, 0.1));
            peak = std::max(peak, std::abs(pr.resonant.last_output()));
        }
        CHECK(std::isfinite(peak));
        CHECK(peak < 1.0);
    }

    TEST_CASE("a u_lim that is not a box leaves PRController invalid") {
        auto res = design::pr({.Kp = 1.0, .Ki = 100.0, .w0 = 314.0, .wc = 10.0}, 0.0001);
        res.u_lim = {.lower = 1.0, .upper = -1.0};
        PRController<double> pr{res};
        CHECK_FALSE(pr.valid());
        CHECK(pr.control(1.0) == 0.0);
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

    TEST_CASE("PRController<float> holds the design's gain and phase from f0/fs = 5e-5 to 0.1") {
        // As Direct Form I sections in float, the first two ran to Inf and the third had 0.78 of the
        // design's gain at ω₀. Float only: the response is measured in double, which does not
        // resolve what a double resonator holds.
        struct Case {
            double f0;
            double fs;
            double wc;
        };
        for (const Case& c : {
                 Case{.f0 = 1.0, .fs = 20e3, .wc = 5.0},
                 Case{.f0 = 5.0, .fs = 100e3, .wc = 10.0},
                 Case{.f0 = 50.0, .fs = 100e3, .wc = 1.0},
                 Case{.f0 = 50.0, .fs = 10e3, .wc = 10.0},
                 Case{.f0 = 2000.0, .fs = 20e3, .wc = 5.0},
             }) {
            const auto res = design::pr<float>(
                {.Kp = 0.0f, .Ki = 100.0f, .w0 = static_cast<float>(2.0 * std::numbers::pi * c.f0), .wc = static_cast<float>(c.wc)},
                static_cast<float>(1.0 / c.fs)
            );
            REQUIRE(res.success);
            // The design is the float numbers it holds: the oracle takes them as they are.
            const double w0 = static_cast<double>(res.w0);
            const double wc = static_cast<double>(res.wc);
            const double Ts = static_cast<double>(res.Ts);
            const double tol = resonator_tol<float>(w0, wc, Ts);

            for (const double w : {w0, w0 - wc, w0 + wc, w0 / 10.0, w0 * 10.0}) {
                if (w * Ts > 0.9 * std::numbers::pi) {
                    continue; // a decade above 0.1·fs is past Nyquist
                }
                CAPTURE(c.f0);
                CAPTURE(c.fs);
                CAPTURE(w / w0);
                PRController<float>        ctrl(res);
                const size_t               two_periods = static_cast<size_t>(damp::ceil(4.0 * std::numbers::pi / (w * Ts)));
                const SineResponse         got = sine_response(ctrl, w, Ts, 1.0, settle_ticks<float>(wc, Ts), two_periods);
                const std::complex<double> want = pr_resonant_response(100.0, w0, wc, Ts, w);
                CHECK(got.finite);
                CHECK(damp::abs(got.gain - std::abs(want)) <= tol * std::abs(want));
                CHECK(damp::abs(got.phase - std::arg(want)) <= tol);
            }
        }
    }

    TEST_CASE("PRController<float> stays on the double result through the runs that used to diverge") {
        // A 0.01 sine at ω₀. As Direct Form I sections the float outputs reached Inf near tick 710000
        // (1 Hz at 20 kHz, ω_c = 5 rad/s) and near tick 430000 (5 Hz at 100 kHz, ω_c = 10 rad/s).
        struct Case {
            double f0;
            double fs;
            double wc;
            size_t ticks;
        };
        for (const Case& c : {
                 Case{.f0 = 1.0, .fs = 20e3, .wc = 5.0, .ticks = 1000000U},
                 Case{.f0 = 5.0, .fs = 100e3, .wc = 10.0, .ticks = 600000U},
             }) {
            CAPTURE(c.f0);
            const auto res = design::pr<float>(
                {.Kp = 1.0f, .Ki = 100.0f, .w0 = static_cast<float>(2.0 * std::numbers::pi * c.f0), .wc = static_cast<float>(c.wc)},
                static_cast<float>(1.0 / c.fs)
            );
            REQUIRE(res.success);
            PRController<float>  ctrl(res);
            PRController<double> ref(res.as<double>()); // the same design: the float numbers, in double
            const double         w0 = static_cast<double>(res.w0);
            const double         Ts = static_cast<double>(res.Ts);
            const double         tol = resonator_tol<float>(w0, static_cast<double>(res.wc), Ts);

            SineSource src(w0, Ts);
            bool       finite = true;
            double     peak = 0.0;  // largest |u_ref| so far: the scale of the difference
            double     worst = 0.0; // largest |u − u_ref| / peak
            for (size_t k = 0; k < c.ticks; ++k) {
                const float  e = static_cast<float>(0.01 * src.s);
                const float  u = ctrl.control(e);
                const double u_ref = ref.control(static_cast<double>(e));
                finite = finite && damp::isfinite(u);
                peak = damp::max(peak, damp::abs(u_ref));
                if (peak > 0.0) {
                    worst = damp::max(worst, damp::abs(static_cast<double>(u) - u_ref) / peak);
                }
                src.step();
            }
            CHECK(finite);
            CHECK(worst <= tol);
            CHECK(peak < 1.05 * 0.01 * 101.0); // the double run itself: Kp + Ki = 101 at ω₀, and a start-up overshoot
        }
    }

    TEST_CASE_TEMPLATE("set_frequency from 50 Hz to 1 Hz lands on the new resonance", T, float, double) {
        // 1 Hz at 20 kHz is f0/fs = 5e-5: in float the retuned Direct Form I section ran to Inf as
        // the design at 1 Hz did, near tick 710000.
        const T    Ts = static_cast<T>(1.0 / 20e3);
        const T    w_old = static_cast<T>(2.0 * std::numbers::pi * 50.0);
        const T    w_new = static_cast<T>(2.0 * std::numbers::pi * 1.0);
        const auto res = design::pr<T>({.Kp = T{1}, .Ki = T{100}, .w0 = w_old, .wc = T{5}}, Ts);
        REQUIRE(res.success);
        PRController<T> ctrl(res);
        const double    wc = 5.0;
        const double    Tsd = static_cast<double>(res.Ts);

        // Settled on 50 Hz first, so the retune finds a running resonator. Kp + Ki = 101 at ω₀, phase 0.
        const double       w50 = static_cast<double>(w_old);
        const SineResponse before = sine_response(ctrl, w50, Tsd, 1.0, settle_ticks<T>(wc, Tsd), 800U);
        CHECK(damp::abs(before.gain - 101.0) <= resonator_tol<T>(w50, wc, Tsd) * 101.0);

        ctrl.set_frequency(w_new);
        CHECK(ctrl.w0 == w_new);

        // A million ticks at 1 Hz, the last two periods fitted: past the tick that reached Inf, and
        // long after the retune has died out (settle_ticks is 180000 in double).
        const double       w1 = static_cast<double>(w_new);
        const SineResponse after = sine_response(ctrl, w1, Tsd, 1.0, 960000U, 40000U);
        const double       tol = resonator_tol<T>(w1, wc, Tsd);
        CHECK(after.finite);
        // The state it kept is the increment of a 50 Hz output: that of a 1 Hz output 50 times larger.
        CHECK(after.peak <= (50.0 + 2.0) * 101.0);
        CHECK(damp::abs(after.gain - 101.0) <= tol * 101.0);
        CHECK(damp::abs(after.phase) <= tol);
    }

    TEST_CASE_TEMPLATE("PR controller tracks a sinusoidal reference in float and double", T, float, double) {
        // Mean |r − y| over the last ticks of a unit sine at ω₀, around a unit-gain plant y[k] = u[k].
        const auto mean_error = [](T Kp, T Ki, T w0, T wc, T Ts, size_t settle, size_t ticks) {
            const auto res = design::pr<T>({.Kp = Kp, .Ki = Ki, .w0 = w0, .wc = wc}, Ts);
            REQUIRE(res.success);
            PRController<T> ctrl(res);
            SineSource      ref(static_cast<double>(w0), static_cast<double>(Ts));
            T               y = T{0};
            double          sum = 0.0;
            for (size_t k = 0; k < settle + ticks; ++k) {
                const T error = static_cast<T>(ref.s) - y;
                y = ctrl.control(error);
                if (k >= settle) {
                    sum += damp::abs(static_cast<double>(error));
                }
                ref.step();
            }
            return sum / static_cast<double>(ticks);
        };

        // The double test above, in T: 50 Hz at 10 kHz, a non-ideal resonator.
        const T w50 = static_cast<T>(2.0 * std::numbers::pi * 50.0);
        CHECK(mean_error(T{1}, T{200}, w50, T{10}, static_cast<T>(1e-4), 15000U, 5000U) < 0.05);

        // 1 Hz at 20 kHz, f0/fs = 5e-5, with an ideal resonator: no error is left at ω₀ but what T
        // holds of it. The loop (1 + Kp)(s² + ω₀²) + Kᵢ·s decays at Kᵢ/(2(1 + Kp)) = 5 rad/s.
        const T      w1 = static_cast<T>(2.0 * std::numbers::pi * 1.0);
        const T      Ts = static_cast<T>(1.0 / 20e3);
        const double rate = 15.0 / (2.0 * 1.5);
        const double Tsd = static_cast<double>(Ts);
        const double left = mean_error(static_cast<T>(0.5), T{15}, w1, T{0}, Ts, settle_ticks<T>(rate, Tsd), 20000U);
        CHECK(left <= resonator_tol<T>(static_cast<double>(w1), rate, Tsd));
    }
}
