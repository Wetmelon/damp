// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <complex>
#include <numbers>
#include <type_traits>
#include <vector>

#include "damp/analysis/analysis.hpp"
#include "damp/backend.hpp"
#include "damp/design/status.hpp"
#include "damp/filters/filters.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {
constexpr double pi = std::numbers::pi_v<double>;

double gain1(const design::FirstOrderCoeffs<double>& c, double w) {
    const std::complex<double> z = std::exp(std::complex<double>(0.0, -w));
    return std::abs((c.b0 + c.b1 * z) / (1.0 + c.a1 * z));
}

double gain2(const design::SecondOrderCoeffs<double>& c, double w) {
    const std::complex<double> z = std::exp(std::complex<double>(0.0, -w));
    return std::abs((c.b0 + c.b1 * z + c.b2 * z * z) / (1.0 + c.a1 * z + c.a2 * z * z));
}

} // namespace

namespace {

// |H(jω)| for a first-order analog TF (highest power first: {b1, b0}).
[[nodiscard]] double mag_jw(const TransferFunction<2, 2, double>& tf, double w) {
    const double nr = tf.num[1];
    const double ni = tf.num[0] * w;
    const double dr = tf.den[1];
    const double di = tf.den[0] * w;
    return damp::hypot(nr, ni) / damp::hypot(dr, di);
}

// |H(jω)| for a second-order analog TF (highest power first: {b2, b1, b0}).
[[nodiscard]] double mag_jw(const TransferFunction<3, 3, double>& tf, double w) {
    const double w2 = w * w;
    const double nr = tf.num[2] - (tf.num[0] * w2);
    const double ni = tf.num[1] * w;
    const double dr = tf.den[2] - (tf.den[0] * w2);
    const double di = tf.den[1] * w;
    return damp::hypot(nr, ni) / damp::hypot(dr, di);
}

template<size_t Nnum, size_t Nden>
void check_tf_finite(const TransferFunction<Nnum, Nden, double>& tf) {
    for (size_t i = 0; i < Nnum; ++i) {
        CHECK(damp::isfinite(tf.num[i]));
    }
    for (size_t i = 0; i < Nden; ++i) {
        CHECK(damp::isfinite(tf.den[i]));
    }
}

} // namespace

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

    TEST_CASE("Invalid fc/Ts are refused with zero first-order coeffs") {
        using design::DesignStatus;
        const auto fc0 = design::lowpass_1st<double>(0.0, 0.001);
        CHECK_FALSE(fc0.success);
        CHECK(fc0.status == DesignStatus::invalid_parameter);
        CHECK(fc0.b0 == doctest::Approx(0.0));
        CHECK(design::lowpass_1st<double>(10.0, 0.0).status == DesignStatus::invalid_sample_time);
        CHECK(design::lowpass_1st<double>(600.0, 0.001).status == DesignStatus::invalid_parameter); // past Nyquist
        CHECK(design::highpass_1st<double>(-1.0, 0.001).status == DesignStatus::invalid_parameter);
        CHECK(design::lowpass_2nd<double>(10.0, 0.0, 0.001).status == DesignStatus::invalid_parameter);
        CHECK(design::lowpass_1st<double>(10.0, 0.001).success);
    }

    TEST_CASE("Second-order low-pass design") {
        const auto coeffs = design::lowpass_2nd<double>(10.0, 0.707, 0.001); // 10 Hz, Butterworth, 1kHz sample rate
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
        constexpr auto   tf = design::pade_delay_1st<double>(T_delay).tf;

        // First-order Pade: H(s) = (1 - sT/2) / (1 + sT/2)
        const double half_T = T_delay / 2.0;

        // Numerator, highest power first: -T/2, 1
        CHECK(tf.num[0] == doctest::Approx(-half_T));
        CHECK(tf.num[1] == doctest::Approx(1.0));

        // Denominator: T/2, 1
        CHECK(tf.den[0] == doctest::Approx(half_T));
        CHECK(tf.den[1] == doctest::Approx(1.0));
    }

    TEST_CASE("Second-order Pade delay") {
        constexpr double T_delay = 0.01; // 10ms delay
        constexpr auto   tf = design::pade_delay_2nd<double>(T_delay).tf;

        // Second-order Pade: H(s) = (1 - sT/2 + (sT)²/12) / (1 + sT/2 + (sT)²/12)
        const double half_T = T_delay / 2.0;
        const double T_sq_12 = T_delay * T_delay / 12.0;

        // Numerator, highest power first: T²/12, -T/2, 1
        CHECK(tf.num[0] == doctest::Approx(T_sq_12));
        CHECK(tf.num[1] == doctest::Approx(-half_T));
        CHECK(tf.num[2] == doctest::Approx(1.0));

        // Denominator: T²/12, T/2, 1
        CHECK(tf.den[0] == doctest::Approx(T_sq_12));
        CHECK(tf.den[1] == doctest::Approx(half_T));
        CHECK(tf.den[2] == doctest::Approx(1.0));
    }

    TEST_CASE("Tustin Padé delays keep their all-pass coefficients in float") {
        //! T = Ts: H(z) = (z²/7 + 4z/7 + 1)/(z² + 4z/7 + 1/7); T = 5·Ts: H(z) = (1 − 2z/3)/(z − 2/3).
        //! Both have unity gain at z = 1, which a numerator from moved zeros loses.
        const auto second = design::pade_delay_2nd<float>(1e-4F, 1e-4F);
        REQUIRE(second.success);
        CHECK(static_cast<double>(second.tf.num[0]) == doctest::Approx(1.0 / 7.0).epsilon(1e-5));
        CHECK(static_cast<double>(second.tf.num[1]) == doctest::Approx(4.0 / 7.0).epsilon(1e-5));
        CHECK(static_cast<double>(second.tf.num[2]) == doctest::Approx(1.0).epsilon(1e-5));
        CHECK(static_cast<double>(second.tf.den[1]) == doctest::Approx(4.0 / 7.0).epsilon(1e-5));
        CHECK(static_cast<double>(second.tf.den[2]) == doctest::Approx(1.0 / 7.0).epsilon(1e-5));

        const auto first = design::pade_delay_1st<float>(5e-4F, 1e-4F);
        REQUIRE(first.success);
        CHECK(static_cast<double>(first.tf.num[0]) == doctest::Approx(-2.0 / 3.0).epsilon(1e-5));
        CHECK(static_cast<double>(first.tf.num[1]) == doctest::Approx(1.0).epsilon(1e-5));
        CHECK(static_cast<double>(first.tf.den[1]) == doctest::Approx(-2.0 / 3.0).epsilon(1e-5));
    }

    TEST_CASE("Continuous first-order low-pass TF is 1/(τs+1)") {
        constexpr double fc = 10.0;
        constexpr auto   tf = design::lowpass_1st<double>(fc).tf;
        const double     omega = 2.0 * damp::numbers::pi_v<double> * fc;
        const double     tau = 1.0 / omega;

        // Highest power first: num = 1, den = τs + 1
        CHECK(tf.num[0] == doctest::Approx(0.0));
        CHECK(tf.num[1] == doctest::Approx(1.0));
        CHECK(tf.den[0] == doctest::Approx(tau));
        CHECK(tf.den[1] == doctest::Approx(1.0));

        // Companion SS: pole at −1/τ = −ω, DC gain C(−A)⁻¹B = 1
        const auto sys = tf.to_state_space().value();
        CHECK(sys.A(0, 0) == doctest::Approx(-omega));
        CHECK(sys.C(0, 0) == doctest::Approx(omega));
        CHECK(sys.D(0, 0) == doctest::Approx(0.0));
    }

    TEST_CASE("Continuous second-order low-pass TF is ω²/(s²+(ω/Q)s+ω²)") {
        constexpr double fc = 25.0;
        constexpr double Q = 0.5;
        constexpr auto   tf = design::lowpass_2nd<double>(fc, Q).tf;
        const double     omega = 2.0 * damp::numbers::pi_v<double> * fc;
        const double     w2 = omega * omega;

        CHECK(tf.num[0] == doctest::Approx(0.0));
        CHECK(tf.num[1] == doctest::Approx(0.0));
        CHECK(tf.num[2] == doctest::Approx(w2));
        // den highest power first: s² + (ω/Q) s + ω² (monic)
        CHECK(tf.den[0] == doctest::Approx(1.0));
        CHECK(tf.den[1] == doctest::Approx(omega / Q));
        CHECK(tf.den[2] == doctest::Approx(w2));

        const auto sys = tf.to_state_space().value();
        CHECK(sys.A(1, 0) == doctest::Approx(-w2));
        CHECK(sys.A(1, 1) == doctest::Approx(-omega / Q));
        CHECK(sys.C(0, 0) == doctest::Approx(w2));
        CHECK(sys.C(0, 1) == doctest::Approx(0.0));

        // The default Q is the Butterworth 1/√2
        const auto def = design::lowpass_2nd<double>(fc).tf;
        const auto exp = design::lowpass_2nd<double>(fc, damp::numbers::inv_sqrt2_v<double>).tf;
        CHECK(def.den[0] == doctest::Approx(exp.den[0]));
        CHECK(def.den[1] == doctest::Approx(exp.den[1]));
        CHECK(def.den[2] == doctest::Approx(exp.den[2]));
    }

    TEST_CASE("Butterworth continuous factors use correct damping") {
        constexpr double fc = 10.0;
        const double     omega = 2.0 * damp::numbers::pi_v<double> * fc;
        const double     inv_sqrt2 = damp::sin(damp::numbers::pi_v<double> / 4.0);

        // n=2: single quadratic Q = 1/√2
        {
            const auto tf = design::lowpass_2nd<double>(fc, inv_sqrt2).tf;
            CHECK(tf.den[1] == doctest::Approx(2.0 * inv_sqrt2 * omega));
            const auto sys = design::butterworth_lowpass<2, double>(fc).sys;
            CHECK(sys.A.rows() == 2);
            CHECK(sys.A(1, 0) == doctest::Approx(-omega * omega));
            CHECK(sys.A(1, 1) == doctest::Approx(-2.0 * inv_sqrt2 * omega));
        }

        // n=3: real pole × quadratic Q = 1 — product den matches cascade
        {
            const auto   tf1 = design::lowpass_1st<double>(fc).tf;
            const auto   tf2 = design::lowpass_2nd<double>(fc, 1.0).tf;
            const auto   prod = tf1 * tf2;
            const double tau = 1.0 / omega;
            CHECK(prod.den.size() == 4); // degree 3
            // (τ s + 1)(s² + ω s + ω²) → leading τ, constant ω²
            CHECK(prod.den[0] == doctest::Approx(tau));
            CHECK(prod.den[3] == doctest::Approx(omega * omega));
            // Companion monicizes by den[0]: pole sum / leading = den[1]/den[0]
            const auto sys = design::butterworth_lowpass<3, double>(fc).sys;
            CHECK(sys.A.rows() == 3);
            CHECK(sys.A(2, 2) == doctest::Approx(-prod.den[1] / prod.den[0]));
        }

        // n=4: ζ = sin(π/8), sin(3π/8) — not two copies of 0.707
        {
            const double z1 = damp::sin(damp::numbers::pi_v<double> / 8.0);
            const double z2 = damp::sin(3.0 * damp::numbers::pi_v<double> / 8.0);
            CHECK(z1 == doctest::Approx(0.38268343236).epsilon(1e-9));
            CHECK(z2 == doctest::Approx(0.92387953251).epsilon(1e-9));
            CHECK(z1 != doctest::Approx(inv_sqrt2).epsilon(1e-3));

            const auto a = design::lowpass_2nd<double>(fc, 1.0 / (2.0 * z1)).tf;
            const auto b = design::lowpass_2nd<double>(fc, 1.0 / (2.0 * z2)).tf;
            const auto prod = a * b;
            const auto sys = design::butterworth_lowpass<4, double>(fc).sys;
            CHECK(sys.A.rows() == 4);
            // Cascade product is monic degree 4 with constant ω⁴
            CHECK(prod.den[0] == doctest::Approx(1.0));
            CHECK(prod.den[4] == doctest::Approx(omega * omega * omega * omega));
            // Companion last-row constant term matches −den[4]/den[0]
            CHECK(sys.A(3, 0) == doctest::Approx(-prod.den[4]));
        }
    }

    TEST_CASE("Butterworth discrete Tustin has unit DC gain") {
        constexpr double fc = 20.0;
        constexpr double Ts = 0.001;
        const auto       sys = design::butterworth_lowpass<3, double>(fc, Ts).sys;
        CHECK(sys.is_discrete());
        // Steady-state gain of discrete SS for constant u: C(I−A)⁻¹B + D
        const auto I = Matrix<3, 3, double>::identity();
        const auto solved = mat::solve(I - sys.A, sys.B);
        REQUIRE(solved.has_value());
        const auto y_ss = sys.C * (*solved) + sys.D;
        CHECK(y_ss(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
    }

    TEST_CASE("Continuous first-order high-pass TF is s/(s+ωc)") {
        constexpr double fc = 10.0;
        constexpr auto   tf = design::highpass_1st<double>(fc).tf;
        const double     omega = 2.0 * damp::numbers::pi_v<double> * fc;

        // Highest power first: num = s, den = s + ωc
        CHECK(tf.num[0] == doctest::Approx(1.0));
        CHECK(tf.num[1] == doctest::Approx(0.0));
        CHECK(tf.den[0] == doctest::Approx(1.0));
        CHECK(tf.den[1] == doctest::Approx(omega));

        CHECK(mag_jw(tf, 0.0) == doctest::Approx(0.0).epsilon(1e-12));          // DC = 0
        CHECK(mag_jw(tf, 100.0 * omega) == doctest::Approx(1.0).epsilon(1e-4)); // HF → 1
        CHECK(mag_jw(tf, omega) == doctest::Approx(damp::numbers::inv_sqrt2_v<double>).epsilon(1e-12));
    }

    TEST_CASE("Continuous second-order high-pass TF is s²/(s²+(ω/Q)s+ω²)") {
        constexpr double fc = 25.0;
        constexpr double Q = 2.0;
        constexpr auto   tf = design::highpass_2nd<double>(fc, Q).tf;
        const double     omega = 2.0 * damp::numbers::pi_v<double> * fc;
        const double     w2 = omega * omega;

        CHECK(tf.num[0] == doctest::Approx(1.0));
        CHECK(tf.num[1] == doctest::Approx(0.0));
        CHECK(tf.num[2] == doctest::Approx(0.0));
        CHECK(tf.den[0] == doctest::Approx(1.0));
        CHECK(tf.den[1] == doctest::Approx(omega / Q));
        CHECK(tf.den[2] == doctest::Approx(w2));

        CHECK(mag_jw(tf, 0.0) == doctest::Approx(0.0).epsilon(1e-12));          // DC = 0
        CHECK(mag_jw(tf, 100.0 * omega) == doctest::Approx(1.0).epsilon(1e-4)); // HF → 1

        // One-arg convenience matches Butterworth Q = 1/√2
        const auto def = design::highpass_2nd<double>(fc).tf;
        const auto exp = design::highpass_2nd<double>(fc, damp::numbers::inv_sqrt2_v<double>).tf;
        CHECK(def.den[2] == doctest::Approx(exp.den[2]));
        CHECK(def.den[1] == doctest::Approx(exp.den[1]));
        CHECK(def.num[0] == doctest::Approx(exp.num[0]));
    }

    TEST_CASE("Continuous band-pass TF peaks at unity at ω0") {
        constexpr double f0 = 40.0;
        constexpr double Q = 5.0;
        constexpr auto   tf = design::bandpass<double>(f0, Q).tf;
        const double     omega = 2.0 * damp::numbers::pi_v<double> * f0;

        CHECK(tf.num[0] == doctest::Approx(0.0));
        CHECK(tf.num[1] == doctest::Approx(omega / Q));
        CHECK(tf.num[2] == doctest::Approx(0.0));
        CHECK(tf.den[0] == doctest::Approx(1.0));

        CHECK(mag_jw(tf, 0.0) == doctest::Approx(0.0).epsilon(1e-12));   // DC = 0
        CHECK(mag_jw(tf, omega) == doctest::Approx(1.0).epsilon(1e-12)); // peak
        CHECK(mag_jw(tf, 100.0 * omega) < 0.01);                         // HF → 0
    }

    TEST_CASE("Continuous notch TF nulls ω0 and is unity at DC/HF") {
        constexpr double f0 = 50.0;
        constexpr double Q = 8.0;
        constexpr auto   tf = design::notch<double>(f0, Q).tf;
        const double     omega = 2.0 * damp::numbers::pi_v<double> * f0;
        const double     w2 = omega * omega;

        CHECK(tf.num[0] == doctest::Approx(1.0));
        CHECK(tf.num[1] == doctest::Approx(0.0));
        CHECK(tf.num[2] == doctest::Approx(w2));
        CHECK(tf.den[0] == doctest::Approx(1.0));
        CHECK(tf.den[1] == doctest::Approx(omega / Q));
        CHECK(tf.den[2] == doctest::Approx(w2));

        CHECK(mag_jw(tf, 0.0) == doctest::Approx(1.0).epsilon(1e-12));          // DC = 1
        CHECK(mag_jw(tf, omega) == doctest::Approx(0.0).epsilon(1e-12));        // |H(jω0)| ≈ 0
        CHECK(mag_jw(tf, 100.0 * omega) == doctest::Approx(1.0).epsilon(1e-4)); // HF → 1
    }

    TEST_CASE("Continuous all-pass TF has |H(jω)| = 1") {
        constexpr double f0 = 30.0;
        constexpr double Q = 0.7;
        constexpr auto   tf = design::allpass<double>(f0, Q).tf;
        const double     omega = 2.0 * damp::numbers::pi_v<double> * f0;

        CHECK(tf.num[0] == doctest::Approx(tf.den[0]));
        CHECK(tf.num[1] == doctest::Approx(-tf.den[1]));
        CHECK(tf.num[2] == doctest::Approx(tf.den[2]));
        CHECK(tf.den[0] == doctest::Approx(1.0));

        CHECK(mag_jw(tf, 0.0) == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(mag_jw(tf, omega) == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(mag_jw(tf, 0.3 * omega) == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(mag_jw(tf, 10.0 * omega) == doctest::Approx(1.0).epsilon(1e-12));
    }

    TEST_CASE("Invalid analog designers are refused") {
        using design::DesignStatus;
        const auto hp1 = design::highpass_1st<double>(0.0);
        CHECK_FALSE(hp1.success);
        CHECK(hp1.status == DesignStatus::invalid_parameter);
        check_tf_finite(hp1.tf);

        CHECK(design::highpass_1st<double>(-5.0).status == DesignStatus::invalid_parameter);
        CHECK(design::lowpass_1st<double>(0.0).status == DesignStatus::invalid_parameter);
        CHECK(design::highpass_2nd<double>(10.0, 0.0).status == DesignStatus::invalid_parameter);
        CHECK(design::bandpass<double>(-1.0, 2.0).status == DesignStatus::invalid_parameter);
        CHECK(design::notch<double>(50.0, -1.0).status == DesignStatus::invalid_parameter);
        CHECK(design::allpass<double>(0.0, 1.0).status == DesignStatus::invalid_parameter);
        CHECK(design::pade_delay_1st<double>(-0.01).status == DesignStatus::invalid_parameter);
        CHECK(design::pade_delay_2nd<double>(0.01, 0.0).status == DesignStatus::invalid_sample_time);
        CHECK(design::butterworth_lowpass<3, double>(0.0).status == DesignStatus::invalid_parameter);
        CHECK(design::butterworth_lowpass<3, double>(600.0, 0.001).status == DesignStatus::invalid_parameter); // past Nyquist
        CHECK(design::butterworth_lowpass<3, double>(10.0, 0.0).status == DesignStatus::invalid_sample_time);
        CHECK(design::lowpass_2nd<double>(10.0).success);
    }


    TEST_CASE("First/second-order bilinear designs put the -3 dB corner at fc (prewarped)") {
        const double fc = 250.0;
        const double Ts = 1e-3; // fc = fs/4: unwarped Tustin gave |H(fc)| = 0.618
        const double w = 2.0 * pi * fc * Ts;
        CHECK(gain1(design::lowpass_1st<double>(fc, Ts), w) == doctest::Approx(1.0 / std::sqrt(2.0)).epsilon(1e-9));
        CHECK(gain1(design::highpass_1st<double>(fc, Ts), w) == doctest::Approx(1.0 / std::sqrt(2.0)).epsilon(1e-9));
        CHECK(gain2(design::lowpass_2nd<double>(fc, damp::numbers::inv_sqrt2_v<double>, Ts), w) == doctest::Approx(1.0 / std::sqrt(2.0)).epsilon(1e-9));
        // Same prewarped bilinear as the RBJ cookbook high-pass: |H(fc)| equal.
        CHECK(gain2(design::highpass_2nd<double>(fc, damp::numbers::inv_sqrt2_v<double>, Ts), w) == doctest::Approx(1.0 / std::sqrt(2.0)).epsilon(1e-9));
    }

    TEST_CASE("Discrete Butterworth corner lands on fc") {
        const double fc = 200.0;
        const double Ts = 1e-3;
        const auto   sys = design::butterworth_lowpass<3, double>(fc, Ts).sys;
        REQUIRE(sys.is_discrete());
        const auto b = analysis::bode(sys, std::vector<double>{2.0 * pi * fc});
        REQUIRE(b.points.size() == 1);
        CHECK(b.points[0].magnitude == doctest::Approx(1.0 / std::sqrt(2.0)).epsilon(1e-9));
    }
} // TEST_SUITE

TEST_SUITE("Runtime Filters") {
    TEST_CASE("LowPass runtime") {
        // 10 Hz cutoff, 1 kHz sample rate, designed in double and deployed in float
        LowPass<1, float> lpf{design::lowpass<1>(10.0, 0.001)};

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

    TEST_CASE("LowPass(fc, Ts) is LowPass(design::lowpass<N>(fc, Ts))") {
        static constinit LowPass<2, float> shorthand{10.0f, 0.001f};
        static constinit LowPass<2, float> designed{design::lowpass<2>(10.0f, 0.001f)};
        for (int k = 0; k < 50; ++k) {
            const float x = (k % 7 == 0) ? 1.0f : -0.25f;
            CHECK(shorthand(x) == designed(x));
        }
        LowPass<1, double> refused{50.0, 0.01}; // fc above Nyquist (50 Hz at 100 Hz): inert at run time
        CHECK(refused(1.0) == 0.0);
    }

    TEST_CASE("LowPass2nd runtime") {
        const auto        coeffs = design::lowpass<2>(10.0f, 0.001f); // 10 Hz cutoff, Butterworth, 1kHz sample rate
        LowPass<2, float> lpf{coeffs};

        lpf.reset();
        float output = lpf(1.0f); // Step input

        // The first sample of a step is the direct feedthrough b0 — a proper
        // 2nd-order low-pass barely responds on sample one. (The old buggy
        // coefficients made b0 ≈ 1, i.e. a near-passthrough; this guards against
        // regressing to that.)
        CHECK(output > 0.0f);
        CHECK(output == doctest::Approx(coeffs.b[0]));
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


    TEST_CASE("Delay<1> reports the delay it actually realizes") {
        Delay<1, double> d;
        CHECK(d.get_delay_samples() == 0);
        CHECK(d(5.0) == doctest::Approx(5.0));
        Delay<3, double> d3;
        CHECK(d3.get_delay_samples() == 1);
        CHECK(d3(5.0) == doctest::Approx(0.0));
        CHECK(d3(6.0) == doctest::Approx(5.0));
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

    TEST_CASE("lowpass_2nd refuses Q = 0 with invalid_parameter and zero coeffs") {
        const auto z = design::lowpass_2nd<double>(10.0, 0.0, 0.001);
        CHECK_FALSE(z.success);
        CHECK(z.status == design::DesignStatus::invalid_parameter);
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

    TEST_CASE("discrete allpass puts −180° at f0, keeps |H| = 1 in float, and rejects invalid") {
        const auto c = design::allpass<double>(50.0, 0.7, 0.001);
        CHECK(damp::isfinite(c.b0));
        CHECK(damp::isfinite(c.a1));
        const double dc = (c.b0 + c.b1 + c.b2) / (1.0 + c.a1 + c.a2);
        CHECK(dc == doctest::Approx(1.0).epsilon(1e-9));

        // Plain Tustin put the −180° point at 286 Hz for f0 = 400 Hz (Ts = 1 ms, Q = 0.7).
        for (const double f0 : {50.0, 150.0, 250.0, 400.0}) {
            const auto                 d = design::allpass<double>(f0, 0.7, 0.001);
            const std::complex<double> z = std::exp(std::complex<double>(0.0, -2.0 * pi * f0 * 0.001));
            const std::complex<double> h = (d.b0 + d.b1 * z + d.b2 * z * z) / (1.0 + d.a1 * z + d.a2 * z * z);
            CHECK(std::abs(h) == doctest::Approx(1.0).epsilon(1e-12));
            CHECK(std::abs(std::arg(h)) == doctest::Approx(pi).epsilon(1e-9));

            const auto                              f = design::allpass<float>(static_cast<float>(f0), 0.7f, 0.001f);
            const design::SecondOrderCoeffs<double> fd{
                .b0 = static_cast<double>(f.b0),
                .b1 = static_cast<double>(f.b1),
                .b2 = static_cast<double>(f.b2),
                .a1 = static_cast<double>(f.a1),
                .a2 = static_cast<double>(f.a2),
            };
            for (const double w : {0.01, 0.5, 1.5, 3.0}) {
                CHECK(gain2(fd, w) == doctest::Approx(1.0).epsilon(1e-12)); // mirrored taps: exact under rounding
            }
        }

        const auto bad = design::allpass<double>(50.0, 0.0, 0.001);
        CHECK(bad.b0 == doctest::Approx(0.0));
        CHECK(bad.a1 == doctest::Approx(0.0));
        const auto nyq = design::allpass<double>(600.0, 1.0, 0.001);
        CHECK(nyq.b0 == doctest::Approx(0.0));
    }

    TEST_CASE("design::lowpass<3/4> in float is stable with unit DC gain at fc/fs = 0.01") {
        //! The impulse-response (Prony) fit returned |p| = 1.02 for lowpass<4, float>(10, 1e-3) with success.
        const auto f3 = design::lowpass<3>(10.0f, 1e-3f);
        const auto f4 = design::lowpass<4>(10.0f, 1e-3f);
        REQUIRE(f3.success);
        REQUIRE(f4.success);
        //! The float numerator moved 20 % per tap with compiler flags; (1 + z⁻¹)^N pins Σb = 1 + Σa.
        const auto dc_gain = [](const auto& r) {
            float num = 0.0f;
            float den = 1.0f;
            for (const float b : r.b) {
                num += b;
            }
            for (const float a : r.a) {
                den += a;
            }
            return static_cast<double>(num / den);
        };
        CHECK(dc_gain(f3) == doctest::Approx(1.0).epsilon(1e-5));
        CHECK(dc_gain(f4) == doctest::Approx(1.0).epsilon(1e-5));
        CHECK(static_cast<double>(f4.b[1]) == doctest::Approx(4.0 * static_cast<double>(f4.b[0])).epsilon(1e-6));
        CHECK(static_cast<double>(f4.b[2]) == doctest::Approx(6.0 * static_cast<double>(f4.b[0])).epsilon(1e-6));

        LowPass<3, float> lp3{f3};
        LowPass<4, float> lp4{f4};
        float             y3 = 0.0f;
        float             y4 = 0.0f;
        for (int k = 0; k < 4000; ++k) {
            y3 = lp3(1.0f);
            y4 = lp4(1.0f);
        }
        // Direct form I in float wanders about 1 % around the DC gain at this pole radius.
        CHECK(damp::abs(y3 - 1.0f) < 0.03f);
        CHECK(damp::abs(y4 - 1.0f) < 0.03f);

        const auto d4 = design::lowpass<4>(10.0, 1e-3);
        REQUIRE(d4.success);
        double num = 0.0;
        double den = 1.0;
        for (size_t i = 0; i < 5; ++i) {
            num += d4.b[i];
        }
        for (size_t i = 0; i < 4; ++i) {
            den += d4.a[i];
        }
        CHECK(num / den == doctest::Approx(1.0).epsilon(1e-9));
    }

    TEST_CASE("design::lowpass refuses a direct form whose rounded poles leave the unit circle") {
        //! lowpass<4, float>(1, 1e-3), (499, 1e-3) and lowpass<3, float>(0.1, 1e-3) reported success and diverged.
        CHECK(design::lowpass<4>(1.0f, 1e-3f).status == design::DesignStatus::invalid_parameter);
        CHECK(design::lowpass<4>(499.0f, 1e-3f).status == design::DesignStatus::invalid_parameter);
        CHECK(design::lowpass<3>(0.1f, 1e-3f).status == design::DesignStatus::invalid_parameter);
        CHECK(design::lowpass<4>(1.0, 1e-3).success); // double holds these poles
        CHECK(design::lowpass<4>(499.0, 1e-3).success);
        // An unstable TF is refused on both conversion paths.
        const TransferFunction<2, 2, double> unstable_1st{{0.0, 1.0}, {1.0, -1.0}};
        const TransferFunction<4, 4, double> unstable_3rd{{0.0, 0.0, 0.0, 1.0}, {1.0, 2.0, -1.0, -2.0}}; // (s−1)(s+1)(s+2)
        CHECK(design::lowpass(unstable_1st, 0.01).status == design::DesignStatus::invalid_parameter);
        CHECK(design::lowpass(unstable_3rd, 0.01).status == design::DesignStatus::invalid_parameter);
    }

    TEST_CASE("design::lowpass(tf, Ts) accepts a non-minimal TF") {
        // (s+1)/((s+1)(s+2)(s+3)) = 1/((s+2)(s+3)): DC gain 1/6; its Hankel matrix is rank-deficient.
        const TransferFunction<4, 4, double> tf{{0.0, 0.0, 1.0, 1.0}, {1.0, 6.0, 11.0, 6.0}};
        const auto                           r = design::lowpass(tf, 0.01);
        REQUIRE(r.success);
        double num = 0.0;
        double den = 1.0;
        for (size_t i = 0; i < 4; ++i) {
            num += r.b[i];
        }
        for (size_t i = 0; i < 3; ++i) {
            den += r.a[i];
        }
        CHECK(num / den == doctest::Approx(1.0 / 6.0).epsilon(1e-9));
    }

    TEST_CASE("Padé delays refuse T_delay = 0 or non-finite as invalid_parameter") {
        CHECK(design::pade_delay_1st<double>(0.0).status == design::DesignStatus::invalid_parameter);
        CHECK(design::pade_delay_2nd<double>(0.0).status == design::DesignStatus::invalid_parameter);
        CHECK(design::pade_delay_1st<double>(0.0, 1e-3).status == design::DesignStatus::invalid_parameter);
        CHECK(design::pade_delay_2nd<double>(0.0, 1e-3).status == design::DesignStatus::invalid_parameter);
    }

    TEST_CASE("design::lowpass(tf, Ts) N=1 uses the to_coeffs path") {
        const auto tf = design::lowpass_1st<double>(10.0).tf;
        const auto res = design::lowpass(tf, 0.001);
        const auto c = design::to_coeffs(tf, 0.001);
        REQUIRE(res.success);
        CHECK(res.b[0] == doctest::Approx(c.b0));
        CHECK(res.b[1] == doctest::Approx(c.b1));
        CHECK(res.a[0] == doctest::Approx(c.a1));
        LowPass via_tf{res};
        static_assert(std::is_same_v<decltype(via_tf), LowPass<1, double>>);
        CHECK(via_tf(1.0) == doctest::Approx(c.b0));
    }

    TEST_CASE("design::lowpass<N>(fc, Ts) refuses an invalid cutoff or sample time") {
        CHECK(design::lowpass<1>(0.0, 0.001).status == design::DesignStatus::invalid_parameter);
        CHECK(design::lowpass<2>(600.0, 0.001).status == design::DesignStatus::invalid_parameter);
        CHECK(design::lowpass<3>(10.0, 0.0).status == design::DesignStatus::invalid_sample_time);
        CHECK(design::lowpass(design::lowpass_1st<double>(10.0).tf, -1.0).status == design::DesignStatus::invalid_sample_time);
        // Third order: unit DC gain, Σb / (1 + Σa) = 1.
        const auto r3 = design::lowpass<3>(20.0, 0.001);
        REQUIRE(r3.success);
        CHECK((r3.b[0] + r3.b[1] + r3.b[2] + r3.b[3]) / (1.0 + r3.a[0] + r3.a[1] + r3.a[2]) == doctest::Approx(1.0).epsilon(1e-9));
    }

    TEST_CASE("LowPass N=3 matches Tustin state space on a step") {
        // 1/(s+1)^3. Orders above 2 recover b from the impulse response.
        TransferFunction<4, 4> tf{.num = {0.0, 0.0, 0.0, 1.0}, .den = {1.0, 3.0, 3.0, 1.0}};
        const double           Ts = 0.01;
        const auto             res = design::lowpass(tf, Ts);
        REQUIRE(res.success);
        // Tustin maps the constant numerator to b ∝ (1 + z⁻¹)³ = {1, 3, 3, 1} (s³ would give {1, −3, 3, −1})
        CHECK(res.b[1] == doctest::Approx(3.0 * res.b[0]).epsilon(1e-9));
        CHECK(res.b[3] == doctest::Approx(res.b[0]).epsilon(1e-9));
        LowPass<3, double> lpf{res};
        const auto         sd = *discretize(tf.to_state_space().value(), Ts, DiscretizationMethod::Tustin);
        ColVec<3>          x{};
        const ColVec<1>    u{{1.0}};
        for (int k = 0; k < 30; ++k) {
            const double y_lpf = lpf(1.0);
            const double y_ss = (sd.C * x + sd.D * u)(0, 0);
            CHECK(y_lpf == doctest::Approx(y_ss).epsilon(1e-6));
            x = sd.A * x + sd.B * u;
        }
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