// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "damp/backend.hpp"
#include "damp/filters/filters.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

/**
 * @brief Tests for filter design and runtime implementations
 */

TEST_SUITE("Filter Design") {
    TEST_CASE("First-order low-pass design") {
        const auto coeffs = design::lowpass_1st<double>(10.0, 0.001); // 10 Hz cutoff, 1kHz sample rate
        // Check that coefficients are reasonable (non-zero)
        CHECK(coeffs.b0 != doctest::Approx(0.0));
        CHECK(coeffs.b1 != doctest::Approx(0.0));
        CHECK(coeffs.a1 != doctest::Approx(0.0));

        // Unit DC gain by construction: (b0+b1) == (1+a1) exactly under -ffast-math
        CHECK((coeffs.b0 + coeffs.b1) == doctest::Approx(1.0 + coeffs.a1).epsilon(1e-12));
        const double dc_gain = (coeffs.b0 + coeffs.b1) / (1.0 + coeffs.a1);
        CHECK(dc_gain == doctest::Approx(1.0).epsilon(1e-12));
    }

    TEST_CASE("First-order high-pass has exact DC null and unit Nyquist by construction") {
        const auto c = design::highpass_1st<double>(10.0, 0.001);
        CHECK(c.b1 == doctest::Approx(-c.b0).epsilon(1e-15)); // DC null shape
        const double dc = (c.b0 + c.b1) / (1.0 + c.a1);
        CHECK(damp::abs(dc) < 1e-12);
        const double nyq = (c.b0 - c.b1) / (1.0 - c.a1);
        CHECK(nyq == doctest::Approx(1.0).epsilon(1e-12));
    }

    TEST_CASE("Invalid fc/Ts yield zero first-order coeffs") {
        CHECK(design::lowpass_1st<double>(0.0, 0.001).b0 == doctest::Approx(0.0));
        CHECK(design::lowpass_1st<double>(10.0, 0.0).b0 == doctest::Approx(0.0));
        CHECK(design::lowpass_1st<double>(600.0, 0.001).b0 == doctest::Approx(0.0)); // past Nyquist
        CHECK(design::highpass_1st<double>(-1.0, 0.001).b0 == doctest::Approx(0.0));
    }

    TEST_CASE("Second-order low-pass design") {
        const auto coeffs = design::lowpass_2nd<double>(10.0, 0.001, 0.707); // 10 Hz, 1kHz sample rate, Butterworth
        // Check that coefficients are reasonable (non-zero)
        CHECK(coeffs.b0 != doctest::Approx(0.0));
        CHECK(coeffs.b1 != doctest::Approx(0.0));
        CHECK(coeffs.b2 != doctest::Approx(0.0));
        CHECK(coeffs.a1 != doctest::Approx(0.0));
        CHECK(coeffs.a2 != doctest::Approx(0.0));

        // DC gain should be 1: (b0 + b1 + b2) / (1 + a1 + a2) ≈ 1
        const double dc_gain = (coeffs.b0 + coeffs.b1 + coeffs.b2) / (1.0 + coeffs.a1 + coeffs.a2);
        CHECK(dc_gain == doctest::Approx(1.0).epsilon(0.01));
    }

    TEST_CASE("First-order Pade delay") {
        constexpr double T_delay = 0.01; // 10ms delay
        constexpr auto   tf = design::pade_delay_1st<double>(T_delay);

        // First-order Pade: H(s) = (1 - sT/2) / (1 + sT/2)
        const double half_T = T_delay / 2.0;

        // Numerator: 1, -T/2
        CHECK(tf.num[0] == doctest::Approx(1.0));
        CHECK(tf.num[1] == doctest::Approx(-half_T));

        // Denominator: 1, T/2
        CHECK(tf.den[0] == doctest::Approx(1.0));
        CHECK(tf.den[1] == doctest::Approx(half_T));
    }

    TEST_CASE("Second-order Pade delay") {
        constexpr double T_delay = 0.01; // 10ms delay
        constexpr auto   tf = design::pade_delay_2nd<double>(T_delay);

        // Second-order Pade: H(s) = (1 - sT/2 + (sT)²/12) / (1 + sT/2 + (sT)²/12)
        const double half_T = T_delay / 2.0;
        const double T_sq_12 = T_delay * T_delay / 12.0;

        // Numerator: 1, -T/2, T²/12
        CHECK(tf.num[0] == doctest::Approx(1.0));
        CHECK(tf.num[1] == doctest::Approx(-half_T));
        CHECK(tf.num[2] == doctest::Approx(T_sq_12));

        // Denominator: 1, T/2, T²/12
        CHECK(tf.den[0] == doctest::Approx(1.0));
        CHECK(tf.den[1] == doctest::Approx(half_T));
        CHECK(tf.den[2] == doctest::Approx(T_sq_12));
    }

    TEST_CASE("Continuous first-order low-pass TF is 1/(τs+1)") {
        constexpr double fc = 10.0;
        constexpr auto   tf = design::lowpass_1st<double>(fc);
        const double     omega = 2.0 * damp::numbers::pi_v<double> * fc;
        const double     tau = 1.0 / omega;

        // Ascending powers: num = 1, den = 1 + τs
        CHECK(tf.num[0] == doctest::Approx(1.0));
        CHECK(tf.num[1] == doctest::Approx(0.0));
        CHECK(tf.den[0] == doctest::Approx(1.0));
        CHECK(tf.den[1] == doctest::Approx(tau));

        // Companion SS: pole at −1/τ = −ω, DC gain C(−A)⁻¹B = 1
        const auto sys = tf.to_state_space().value();
        CHECK(sys.A(0, 0) == doctest::Approx(-omega));
        CHECK(sys.C(0, 0) == doctest::Approx(omega));
        CHECK(sys.D(0, 0) == doctest::Approx(0.0));
    }

    TEST_CASE("Continuous second-order low-pass TF is ω²/(s²+2ζωs+ω²)") {
        constexpr double fc = 25.0;
        constexpr double zeta = 0.5;
        constexpr auto   tf = design::lowpass_2nd_continuous<double>(fc, zeta);
        const double     omega = 2.0 * damp::numbers::pi_v<double> * fc;
        const double     w2 = omega * omega;

        CHECK(tf.num[0] == doctest::Approx(w2));
        CHECK(tf.num[1] == doctest::Approx(0.0));
        CHECK(tf.num[2] == doctest::Approx(0.0));
        // den ascending: ω² + 2ζω s + s² (monic in s²)
        CHECK(tf.den[0] == doctest::Approx(w2));
        CHECK(tf.den[1] == doctest::Approx(2.0 * zeta * omega));
        CHECK(tf.den[2] == doctest::Approx(1.0));

        const auto sys = tf.to_state_space().value();
        CHECK(sys.A(1, 0) == doctest::Approx(-w2));
        CHECK(sys.A(1, 1) == doctest::Approx(-2.0 * zeta * omega));
        CHECK(sys.C(0, 0) == doctest::Approx(w2));
        CHECK(sys.C(0, 1) == doctest::Approx(0.0));

        // One-arg convenience matches default-ζ continuous
        const auto def = design::lowpass_2nd<double>(fc);
        const auto exp = design::lowpass_2nd_continuous<double>(fc);
        CHECK(def.den[0] == doctest::Approx(exp.den[0]));
        CHECK(def.den[1] == doctest::Approx(exp.den[1]));
    }

    TEST_CASE("Butterworth continuous factors use correct damping") {
        constexpr double fc = 10.0;
        const double     omega = 2.0 * damp::numbers::pi_v<double> * fc;
        const double     inv_sqrt2 = damp::sin(damp::numbers::pi_v<double> / 4.0);

        // n=2: single quadratic ζ = 1/√2
        {
            const auto tf = design::lowpass_2nd_continuous<double>(fc, inv_sqrt2);
            CHECK(tf.den[1] == doctest::Approx(2.0 * inv_sqrt2 * omega));
            const auto sys = design::butterworth_lowpass<2, double>(fc);
            CHECK(sys.A.rows() == 2);
            CHECK(sys.A(1, 0) == doctest::Approx(-omega * omega));
            CHECK(sys.A(1, 1) == doctest::Approx(-2.0 * inv_sqrt2 * omega));
        }

        // n=3: real pole × quadratic ζ = 1/2 — product den matches cascade
        {
            const auto   tf1 = design::lowpass_1st<double>(fc);
            const auto   tf2 = design::lowpass_2nd_continuous<double>(fc, 0.5);
            const auto   prod = tf1 * tf2;
            const double tau = 1.0 / omega;
            CHECK(prod.den.size() == 4); // degree 3
            // (1 + τ s)(ω² + ω s + s²) → leading τ, constant ω²
            CHECK(prod.den[0] == doctest::Approx(omega * omega));
            CHECK(prod.den[3] == doctest::Approx(tau));
            // Companion monicizes by den[3]: pole sum / leading = den[2]/den[3]
            const auto sys = design::butterworth_lowpass<3, double>(fc);
            CHECK(sys.A.rows() == 3);
            CHECK(sys.A(2, 2) == doctest::Approx(-prod.den[2] / prod.den[3]));
        }

        // n=4: ζ = sin(π/8), sin(3π/8) — not two copies of 0.707
        {
            const double z1 = damp::sin(damp::numbers::pi_v<double> / 8.0);
            const double z2 = damp::sin(3.0 * damp::numbers::pi_v<double> / 8.0);
            CHECK(z1 == doctest::Approx(0.38268343236).epsilon(1e-9));
            CHECK(z2 == doctest::Approx(0.92387953251).epsilon(1e-9));
            CHECK(z1 != doctest::Approx(inv_sqrt2).epsilon(1e-3));

            const auto a = design::lowpass_2nd_continuous<double>(fc, z1);
            const auto b = design::lowpass_2nd_continuous<double>(fc, z2);
            const auto prod = a * b;
            const auto sys = design::butterworth_lowpass<4, double>(fc);
            CHECK(sys.A.rows() == 4);
            // Cascade product is monic degree 4 with constant ω⁴
            CHECK(prod.den[0] == doctest::Approx(omega * omega * omega * omega));
            CHECK(prod.den[4] == doctest::Approx(1.0));
            // Companion last-row constant term matches −den[0]/den[4]
            CHECK(sys.A(3, 0) == doctest::Approx(-prod.den[0]));
        }
    }

    TEST_CASE("Butterworth discrete Tustin has unit DC gain") {
        constexpr double fc = 20.0;
        constexpr double Ts = 0.001;
        const auto       sys = design::butterworth_lowpass<3, double>(fc, Ts);
        CHECK(sys.is_discrete());
        // Steady-state gain of discrete SS for constant u: C(I−A)⁻¹B + D
        const auto I = Matrix<3, 3, double>::identity();
        const auto solved = mat::solve(I - sys.A, sys.B);
        REQUIRE(solved.has_value());
        const auto y_ss = sys.C * (*solved) + sys.D;
        CHECK(y_ss(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
    }

} // TEST_SUITE

TEST_SUITE("Runtime Filters") {
    TEST_CASE("LowPass runtime") {
        // Hardcoded coeffs for 10 Hz cutoff, 1kHz sample rate
        design::FirstOrderCoeffs<float> coeffs{.b0 = 0.030418f, .b1 = 0.030418f, .a1 = -0.938874f};
        LowPass<1, float>               lpf(std::array<float, 2>{coeffs.b0, coeffs.b1}, std::array<float, 1>{coeffs.a1});

        // Test step response
        lpf.reset();
        float output = lpf(1.0f); // Step input

        // Should start responding
        CHECK(output > 0.0f);
        CHECK(output < 1.0f);

        // Let it settle
        for (int i = 0; i < 1000; ++i) {
            output = lpf(1.0f);
        }

        // Should be close to 1.0 (DC gain = 1)
        CHECK(output == doctest::Approx(1.0).epsilon(0.01));
    }

    TEST_CASE("LowPass2nd runtime") {
        auto              coeffs = design::lowpass_2nd(10.0f, 0.001f, 0.707f).as<float>(); // 10 Hz cutoff, Butterworth, 1kHz sample rate
        LowPass<2, float> lpf({coeffs.b0, coeffs.b1, coeffs.b2}, {coeffs.a1, coeffs.a2});

        lpf.reset();
        float output = lpf(1.0f); // Step input

        // The first sample of a step is the direct feedthrough b0 — a proper
        // 2nd-order low-pass barely responds on sample one. (The old buggy
        // coefficients made b0 ≈ 1, i.e. a near-passthrough; this guards against
        // regressing to that.)
        CHECK(output > 0.0f);
        CHECK(output == doctest::Approx(coeffs.b0));
        CHECK(output < 0.01f); // nowhere near settled yet

        // Let it settle
        for (int i = 0; i < 2000; ++i) {
            output = lpf(1.0f);
        }

        // Should settle at 1.0. lowpass_2nd derives its numerator taps from the
        // unit-DC-gain identity (b0 + b1 + b2) = (1 + a1 + a2) by construction, so
        // unity DC gain holds exactly even under the test runner's -ffast-math (the
        // taps are dc_sum·{1/4, 1/2, 1/4}, exact powers of two).
        CHECK(output == doctest::Approx(1.0).epsilon(0.01));
    }

    TEST_CASE("Discrete delay runtime") {
        Delay<10, float> delay; // Max 10 sample delay
        delay.init(3);          // 3 sample delay

        delay.reset();

        // Test impulse response
        float output = delay(1.0f);            // Input impulse at t=0
        CHECK(output == doctest::Approx(0.0)); // Should output 0 (no history)

        output = delay(0.0f);                  // t=1
        CHECK(output == doctest::Approx(0.0)); // Still no output

        output = delay(0.0f);                  // t=2
        CHECK(output == doctest::Approx(0.0)); // Still no output

        output = delay(0.0f);                  // t=3
        CHECK(output == doctest::Approx(1.0)); // Should output the impulse from t=0

        output = delay(0.0f);                  // t=4
        CHECK(output == doctest::Approx(0.0)); // Back to zero

        // Test reset
        delay.reset();
        output = delay(2.0f); // New impulse
        CHECK(output == doctest::Approx(0.0));

        output = delay(0.0f);
        CHECK(output == doctest::Approx(0.0));

        output = delay(0.0f);
        CHECK(output == doctest::Approx(0.0));

        output = delay(0.0f);
        CHECK(output == doctest::Approx(2.0)); // Should output the reset impulse
    }

} // TEST_SUITE

TEST_SUITE("Discrete designer guards") {
    TEST_CASE("lowpass_1st invalid Ts/fc returns zero coeffs") {
        const auto z_ts = design::lowpass_1st<double>(10.0, 0.0);
        CHECK(z_ts.b0 == doctest::Approx(0.0));
        CHECK(z_ts.b1 == doctest::Approx(0.0));
        CHECK(z_ts.a1 == doctest::Approx(0.0));

        const auto z_fc = design::lowpass_1st<double>(0.0, 0.001);
        CHECK(z_fc.b0 == doctest::Approx(0.0));

        // fc at/above Nyquist (fs=1000 → Nyquist 500)
        const auto z_nyq = design::lowpass_1st<double>(500.0, 0.001);
        CHECK(z_nyq.b0 == doctest::Approx(0.0));
    }

    TEST_CASE("lowpass_2nd invalid zeta returns zero coeffs") {
        const auto z = design::lowpass_2nd<double>(10.0, 0.001, 0.0);
        CHECK(z.b0 == doctest::Approx(0.0));
        CHECK(z.a1 == doctest::Approx(0.0));
    }

    TEST_CASE("highpass_1st matches Tustin and rejects invalid") {
        const auto ok = design::highpass_1st<double>(10.0, 0.001);
        // DC gain of high-pass: (b0+b1)/(1+a1) = 0
        const double dc = (ok.b0 + ok.b1) / (1.0 + ok.a1);
        CHECK(dc == doctest::Approx(0.0).epsilon(1e-12));
        CHECK(ok.b0 > 0.0);

        const auto bad = design::highpass_1st<double>(-1.0, 0.001);
        CHECK(bad.b0 == doctest::Approx(0.0));
    }

    TEST_CASE("RBJ designers zero on bad Q or Nyquist") {
        const auto n = design::notch<double>(50.0, 0.0, 0.001);
        CHECK(n.b0 == doctest::Approx(0.0));
        const auto bp = design::bandpass<double>(600.0, 5.0, 0.001); // above Nyquist
        CHECK(bp.b0 == doctest::Approx(0.0));
    }

    TEST_CASE("LowPass(tf,Ts) N=1 uses to_coeffs path") {
        const auto         tf = design::lowpass_1st<double>(10.0);
        LowPass<1, double> via_tf(tf, 0.001);
        const auto         c = design::to_coeffs(tf, 0.001);
        LowPass<1, double> via_c(c);

        // Same first sample on a unit step (both start from zero state)
        CHECK(via_tf(1.0) == doctest::Approx(via_c(1.0)));
    }
}

TEST_SUITE("Coeff <-> StateSpace conversion") {
    TEST_CASE("First-order coeffs round-trip through a discrete state space") {
        const design::FirstOrderCoeffs<double> c{0.3, 0.2, -0.5};
        const auto                             sys = c.to_state_space(0.001);
        CHECK(sys.is_discrete());

        const auto c2 = design::to_coeffs(sys); // sys already discrete: no re-discretization
        CHECK(c2.b0 == doctest::Approx(c.b0));
        CHECK(c2.b1 == doctest::Approx(c.b1));
        CHECK(c2.a1 == doctest::Approx(c.a1));
    }

    TEST_CASE("Second-order coeffs round-trip with D != 0") {
        // notch has b0 != 0, so the D-dependent numerator terms (b1's −D·trace and
        // the full b2) must be correct — the case the old to_coeffs got wrong.
        const auto c = design::notch<double>(50.0, 5.0, 1.0 / 1000.0);
        REQUIRE(c.b0 != doctest::Approx(0.0));

        const auto sys = c.to_state_space(0.001);
        const auto c2 = design::to_coeffs(sys);
        CHECK(c2.b0 == doctest::Approx(c.b0));
        CHECK(c2.b1 == doctest::Approx(c.b1));
        CHECK(c2.b2 == doctest::Approx(c.b2));
        CHECK(c2.a1 == doctest::Approx(c.a1));
        CHECK(c2.a2 == doctest::Approx(c.a2));
    }

    TEST_CASE("Biquad coeffs and their state-space realization share an impulse response") {
        const double   Ts = 1.0 / 2000.0;
        const auto     c = design::peaking<double>(100.0, 2.0, 6.0, Ts); // D != 0
        Biquad<double> bq{c};
        const auto     sys = c.to_state_space(Ts);

        ColVec<2, double> x{};
        for (int k = 0; k < 32; ++k) {
            const double u = (k == 0) ? 1.0 : 0.0; // impulse
            const double y_bq = bq(u);
            const double y_ss = (sys.C * x)(0, 0) + (sys.D(0, 0) * u);
            x = (sys.A * x) + (sys.B * u);
            CHECK(y_ss == doctest::Approx(y_bq).epsilon(1e-9));
        }
    }
}