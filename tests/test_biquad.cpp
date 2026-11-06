// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/design/status.hpp"
#include "damp/filters/filters.hpp"
#include "damp/math/complex.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/systems/zpk.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

// |H(e^{jω})| for a normalized biquad at frequency f (Hz), sample time Ts.
double biquad_mag(const design::SecondOrderCoeffs<double>& c, double f, double Ts) {
    const double               w = 2.0 * std::numbers::pi * f * Ts;
    const std::complex<double> z1 = std::exp(std::complex<double>(0.0, -w));
    const std::complex<double> num = c.b0 + c.b1 * z1 + c.b2 * z1 * z1;
    const std::complex<double> den = 1.0 + c.a1 * z1 + c.a2 * z1 * z1;
    return std::abs(num / den);
}

constexpr double Ts = 1.0 / 1000.0; // 1 kHz

// Closed-form gains at the band edges, evaluated from the *stored* coefficients
// the way the difference equation does (z = +1 at DC, z = −1 at Nyquist). These
// are the exact identities a fast-math reassociation could disturb.
double dc_gain(const design::SecondOrderCoeffs<double>& c) {
    return (c.b0 + c.b1 + c.b2) / (1.0 + c.a1 + c.a2);
}
double nyquist_gain(const design::SecondOrderCoeffs<double>& c) {
    return (c.b0 - c.b1 + c.b2) / (1.0 - c.a1 + c.a2);
}

} // namespace

TEST_SUITE("Biquad Designs") {
    TEST_CASE("notch rejects f0 and passes the rest") {
        const auto c = design::notch<double>(50.0, 5.0, Ts);
        CHECK(biquad_mag(c, 50.0, Ts) == doctest::Approx(0.0).epsilon(1e-6)); // deep null at f0
        CHECK(biquad_mag(c, 10.0, Ts) == doctest::Approx(1.0).epsilon(0.05)); // passband
        CHECK(biquad_mag(c, 200.0, Ts) == doctest::Approx(1.0).epsilon(0.05));
    }

    TEST_CASE("bandpass peaks at f0 with unity gain") {
        const auto c = design::bandpass<double>(50.0, 5.0, Ts);
        CHECK(biquad_mag(c, 50.0, Ts) == doctest::Approx(1.0).epsilon(1e-6)); // 0 dB peak
        CHECK(biquad_mag(c, 5.0, Ts) < 0.3);                                  // rejects DC side
        CHECK(biquad_mag(c, 400.0, Ts) < 0.3);                                // rejects HF side
    }

    TEST_CASE("highpass blocks DC and passes high frequencies") {
        const auto c = design::highpass_2nd<double>(100.0, damp::numbers::inv_sqrt2_v<double>, Ts);
        CHECK(biquad_mag(c, 0.0, Ts) == doctest::Approx(0.0).epsilon(1e-9)); // exact zero at DC
        CHECK(biquad_mag(c, 480.0, Ts) == doctest::Approx(1.0).epsilon(0.05));
    }

    TEST_CASE("peaking boosts f0 by the requested dB") {
        const auto   c = design::peaking<double>(50.0, 5.0, 6.0, Ts);
        const double expected = std::pow(10.0, 6.0 / 20.0); // +6 dB ≈ 1.995
        CHECK(biquad_mag(c, 50.0, Ts) == doctest::Approx(expected).epsilon(1e-3));
        CHECK(biquad_mag(c, 5.0, Ts) == doctest::Approx(1.0).epsilon(0.02)); // unity away from f0
    }

    TEST_CASE("low-shelf boosts below corner, unity above") {
        const auto   c = design::lowshelf<double>(100.0, 6.0, Ts);
        const double expected = std::pow(10.0, 6.0 / 20.0);
        CHECK(biquad_mag(c, 1.0, Ts) == doctest::Approx(expected).epsilon(0.02)); // DC boosted
        CHECK(biquad_mag(c, 490.0, Ts) == doctest::Approx(1.0).epsilon(0.05));    // HF unity
    }

    TEST_CASE("high-shelf boosts above corner, unity below") {
        const auto   c = design::highshelf<double>(100.0, 6.0, Ts);
        const double expected = std::pow(10.0, 6.0 / 20.0);
        CHECK(biquad_mag(c, 490.0, Ts) == doctest::Approx(expected).epsilon(0.02)); // HF boosted
        CHECK(biquad_mag(c, 1.0, Ts) == doctest::Approx(1.0).epsilon(0.02));        // DC unity
    }

    // The test runner builds with -ffast-math (see tests/Tupfile.lua). These check
    // each designer's exact band-edge gain identity straight from the stored
    // coefficients, tightly — guarding against an associative-math reordering
    // quietly breaking a cancellation the property relies on (roadmap #17 sweep).
    TEST_CASE("band-edge gain identities hold tightly under -ffast-math") {
        // notch: unity passband at both DC and Nyquist.
        const auto n = design::notch<double>(50.0, 5.0, Ts);
        CHECK(dc_gain(n) == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(nyquist_gain(n) == doctest::Approx(1.0).epsilon(1e-9));

        // bandpass: blocks DC and Nyquist exactly (numerator {α, 0, −α}).
        const auto bp = design::bandpass<double>(50.0, 5.0, Ts);
        CHECK(std::abs(dc_gain(bp)) < 1e-12);
        CHECK(std::abs(nyquist_gain(bp)) < 1e-12);

        // highpass: zero gain at DC.
        const auto hp = design::highpass_2nd<double>(100.0, damp::numbers::inv_sqrt2_v<double>, Ts);
        CHECK(std::abs(dc_gain(hp)) < 1e-9);

        // lowpass_2nd: unity DC gain by construction (the lead #17 fix).
        const auto lp = design::lowpass_2nd<double>(100.0, 0.707, Ts); // (fc, zeta, Ts) coeffs overload
        CHECK(dc_gain(lp) == doctest::Approx(1.0).epsilon(1e-12));

        // peaking & shelves: unity gain away from the affected band.
        const auto pk = design::peaking<double>(50.0, 5.0, 6.0, Ts);
        CHECK(dc_gain(pk) == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(nyquist_gain(pk) == doctest::Approx(1.0).epsilon(1e-9));

        const auto ls = design::lowshelf<double>(100.0, 6.0, Ts);
        CHECK(nyquist_gain(ls) == doctest::Approx(1.0).epsilon(1e-9)); // unity above corner

        const auto hs = design::highshelf<double>(100.0, 6.0, Ts);
        CHECK(dc_gain(hs) == doctest::Approx(1.0).epsilon(1e-9)); // unity below corner
    }

    TEST_CASE("notch coefficients are symmetric and constexpr") {
        constexpr auto c = design::notch<double>(50.0, 5.0, Ts);
        static_assert(c.b0 == c.b2, "notch numerator is symmetric");
        static_assert(c.b1 == c.a1, "notch shares the −2cosω₀ term");
        CHECK(c.b0 == doctest::Approx(c.b2));
    }

    TEST_CASE("float and double designs agree") {
        const auto cf = design::notch<float>(50.0f, 5.0f, 1.0f / 1000.0f);
        const auto cd = design::notch<double>(50.0, 5.0, Ts);
        CHECK(cf.b0 == doctest::Approx(static_cast<float>(cd.b0)).epsilon(1e-4));
        CHECK(cf.a2 == doctest::Approx(static_cast<float>(cd.a2)).epsilon(1e-4));
    }
}

TEST_SUITE("Biquad Runtime") {
    TEST_CASE("Biquad impulse response starts at b0") {
        const auto     c = design::notch<double>(50.0, 5.0, Ts);
        Biquad<double> bq(c);
        CHECK(bq(1.0) == doctest::Approx(c.b0)); // first impulse sample = b0
    }

    TEST_CASE("Biquad attenuates a tone at the notch frequency") {
        const auto     c = design::notch<double>(50.0, 10.0, Ts);
        Biquad<double> bq(c);

        // Drive a 50 Hz sine; after settling the output amplitude should be tiny.
        double peak = 0.0;
        for (int n = 0; n < 2000; ++n) {
            const double t = n * Ts;
            const double y = bq(std::sin(2.0 * std::numbers::pi * 50.0 * t));
            if (n > 1500 && std::abs(y) > peak) {
                peak = std::abs(y);
            }
        }
        CHECK(peak < 0.1); // 50 Hz tone strongly rejected
    }

    TEST_CASE("BiquadCascade chains sections (two notches)") {
        const std::array<design::SecondOrderCoeffs<double>, 2> sections{
            design::notch<double>(50.0, 8.0, Ts),
            design::notch<double>(150.0, 8.0, Ts)
        };
        BiquadCascade<2, double> cascade(sections);

        // Cascade gain = product of section gains.
        const double g50 = biquad_mag(sections[0], 50.0, Ts) * biquad_mag(sections[1], 50.0, Ts);
        const double g150 = biquad_mag(sections[0], 150.0, Ts) * biquad_mag(sections[1], 150.0, Ts);
        CHECK(g50 < 0.05);
        CHECK(g150 < 0.05);

        // Runtime impulse first sample = product of the two b0 terms.
        CHECK(cascade(1.0) == doctest::Approx(sections[0].b0 * sections[1].b0));
    }
}

namespace {

// |H(e^{jω})| of a cascade: the product of its sections' magnitudes.
template<size_t NS>
double cascade_mag(const design::SecondOrderSectionsResult<NS, double>& sos, double f, double period) {
    double m = 1.0;
    for (const auto& c : sos.sections) {
        m *= biquad_mag(c, f, period);
    }
    return m;
}

// Impulse response of a transfer function in z run as one direct form (the oracle).
template<size_t K, size_t Nnum, size_t Nden>
std::array<double, K> direct_impulse(const TransferFunction<Nnum, Nden, double>& tf) {
    constexpr size_t         M = Nden - 1;
    std::array<double, Nden> a{};
    std::array<double, Nden> b{}; // z⁻¹ taps: a shorter numerator is a delay
    for (size_t i = 0; i < Nden; ++i) {
        a[i] = tf.den[i] / tf.den[0];
    }
    for (size_t i = 0; i < Nnum; ++i) {
        b[i + (Nden - Nnum)] = tf.num[i] / tf.den[0];
    }
    std::array<double, K> x{};
    std::array<double, K> y{};
    x[0] = 1.0;
    for (size_t n = 0; n < K; ++n) {
        double acc = 0.0;
        for (size_t i = 0; i <= M && i <= n; ++i) {
            acc += b[i] * x[n - i];
        }
        for (size_t i = 1; i <= M && i <= n; ++i) {
            acc -= a[i] * y[n - i];
        }
        y[n] = acc;
    }
    return y;
}

template<size_t K, size_t NS>
std::array<double, K> cascade_impulse(const design::SecondOrderSectionsResult<NS, double>& sos) {
    BiquadCascade         run{sos};
    std::array<double, K> y{};
    for (size_t n = 0; n < K; ++n) {
        y[n] = run(n == 0 ? 1.0 : 0.0);
    }
    return y;
}

// A 5th-order filter in z with complex and real poles and zeros, from its factors.
constexpr ZPK<5, 5, double> kZpk{
    .zeros = {damp::complex<double>{-1.0, 0.0}, damp::complex<double>{0.3, 0.8}, damp::complex<double>{0.3, -0.8}, damp::complex<double>{-0.5, 0.0}, damp::complex<double>{0.9, 0.0}},
    .poles = {damp::complex<double>{0.95, 0.2}, damp::complex<double>{0.95, -0.2}, damp::complex<double>{0.6, 0.0}, damp::complex<double>{-0.4, 0.5}, damp::complex<double>{-0.4, -0.5}},
    .gain = 0.25,
};

} // namespace

TEST_SUITE("Second-order sections") {
    TEST_CASE("lowpass_sections is the bilinear Butterworth at every order") {
        const double fc = 50.0;
        auto         check = [&](const auto& sos, size_t N) {
            REQUIRE(sos.success);
            for (const double f : {5.0, 25.0, 50.0, 80.0, 200.0, 450.0}) {
                const double r = std::tan(std::numbers::pi * f * Ts) / std::tan(std::numbers::pi * fc * Ts);
                const double expected = 1.0 / std::sqrt(1.0 + std::pow(r, 2.0 * static_cast<double>(N)));
                CHECK(cascade_mag(sos, f, Ts) == doctest::Approx(expected).epsilon(1e-9));
            }
        };
        check(design::lowpass_sections<1>(fc, Ts), 1);
        check(design::lowpass_sections<2>(fc, Ts), 2);
        check(design::lowpass_sections<3>(fc, Ts), 3);
        check(design::lowpass_sections<6>(fc, Ts), 6);
        check(design::lowpass_sections<7>(fc, Ts), 7);
    }

    TEST_CASE("highpass_sections is the bilinear Butterworth high-pass") {
        const double fc = 50.0;
        const auto   sos = design::highpass_sections<5>(fc, Ts);
        REQUIRE(sos.success);
        for (const double f : {5.0, 50.0, 120.0, 400.0}) {
            const double r = std::tan(std::numbers::pi * fc * Ts) / std::tan(std::numbers::pi * f * Ts);
            CHECK(cascade_mag(sos, f, Ts) == doctest::Approx(1.0 / std::sqrt(1.0 + std::pow(r, 10.0))).epsilon(1e-9));
        }
    }

    TEST_CASE("lowpass_sections<4> matches the direct-form design::lowpass<4>") {
        const auto sos = design::lowpass_sections<4>(30.0, Ts);
        const auto df = design::lowpass<4>(30.0, Ts);
        REQUIRE(sos.success);
        REQUIRE(df.success);
        BiquadCascade      cascade{sos};
        LowPass<4, double> direct{df};
        for (size_t n = 0; n < 400; ++n) {
            const double x = (n < 200) ? 1.0 : -0.5;
            CHECK(cascade(x) == doctest::Approx(direct(x)).epsilon(1e-6)); // the direct form is a Prony fit
        }
    }

    TEST_CASE("Butterworth sections run from the lowest Q to the highest") {
        const auto sos = design::lowpass_sections<6>(50.0, Ts);
        REQUIRE(sos.success);
        // A complex pair's pole radius is √a₂: farthest from the unit circle first (MATLAB's order up).
        CHECK(std::sqrt(sos.sections[0].a2) < std::sqrt(sos.sections[1].a2));
        CHECK(std::sqrt(sos.sections[1].a2) < std::sqrt(sos.sections[2].a2));
        const auto odd = design::lowpass_sections<5>(50.0, Ts);
        CHECK(odd.sections[0].a2 == 0.0); // the first-order section first
        CHECK(odd.sections[0].b2 == 0.0);
    }

    TEST_CASE("second_order_sections of a TF runs as the TF does") {
        const auto tf = kZpk.to_transfer_function();
        const auto sos = design::second_order_sections(tf);
        REQUIRE(sos.success);
        static_assert(std::is_same_v<decltype(sos), const design::SecondOrderSectionsResult<3, double>>);
        const auto want = direct_impulse<200>(tf);
        const auto got = cascade_impulse<200>(sos);
        for (size_t n = 0; n < 200; ++n) {
            CHECK(got[n] == doctest::Approx(want[n]).epsilon(1e-9).scale(1.0));
        }
        // The poles closest to the unit circle (|p| ≈ 0.97) are in the last section, with the zeros
        // nearest them: the real pair −1, 0.9, so (z + 1)(z − 0.9) = z² + 0.1z − 0.9.
        CHECK(sos.sections[2].a2 == doctest::Approx((0.95 * 0.95) + (0.2 * 0.2)));
        CHECK(sos.sections[2].b1 / sos.sections[2].b0 == doctest::Approx(0.1));
        CHECK(sos.sections[2].b2 / sos.sections[2].b0 == doctest::Approx(-0.9));
    }

    TEST_CASE("second_order_sections of the ZPK gives the TF's sections") {
        const auto from_zpk = design::second_order_sections(kZpk);
        const auto from_tf = design::second_order_sections(kZpk.to_transfer_function());
        REQUIRE(from_zpk.success);
        for (size_t k = 0; k < 3; ++k) {
            CHECK(from_zpk.sections[k].b0 == doctest::Approx(from_tf.sections[k].b0).epsilon(1e-9));
            CHECK(from_zpk.sections[k].b1 == doctest::Approx(from_tf.sections[k].b1).epsilon(1e-9));
            CHECK(from_zpk.sections[k].a1 == doctest::Approx(from_tf.sections[k].a1).epsilon(1e-9));
            CHECK(from_zpk.sections[k].a2 == doctest::Approx(from_tf.sections[k].a2).epsilon(1e-9));
        }
    }

    TEST_CASE("a numerator of lower degree keeps its delay") {
        // (z − 0.5) / ((z − 0.9)(z + 0.3)): one sample of delay, h[0] = 0
        const TransferFunction<2, 3, double> tf{.num = {1.0, -0.5}, .den = {1.0, -0.6, -0.27}};
        const auto                           sos = design::second_order_sections(tf);
        REQUIRE(sos.success);
        const auto want = direct_impulse<50>(tf);
        const auto got = cascade_impulse<50>(sos);
        CHECK(got[0] == 0.0);
        for (size_t n = 0; n < 50; ++n) {
            CHECK(got[n] == doctest::Approx(want[n]).epsilon(1e-12).scale(1.0));
        }
    }

    TEST_CASE("Thiran fractional delay runs as sections") {
        const auto delay = fractional_delay<5>(4.6e-3, 1e-3);
        REQUIRE(delay.has_value());
        const auto sos = design::second_order_sections(*delay);
        REQUIRE(sos.success);
        const auto want = direct_impulse<120>(*delay);
        const auto got = cascade_impulse<120>(sos);
        for (size_t n = 0; n < 120; ++n) {
            CHECK(got[n] == doctest::Approx(want[n]).epsilon(1e-9).scale(1.0));
        }
        // All-pass: unit gain at any frequency
        CHECK(cascade_mag(sos, 37.0, 1e-3) == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(cascade_mag(sos, 410.0, 1e-3) == doctest::Approx(1.0).epsilon(1e-9));
    }

    TEST_CASE("a whole-sample Thiran delay is an exact pure delay") {
        const auto delay = fractional_delay<4>(1.0, 0.25); // D = 4 samples exactly: z⁻⁴
        REQUIRE(delay.has_value());
        const auto sos = design::second_order_sections(*delay);
        REQUIRE(sos.success);
        BiquadCascade run{sos};
        for (int n = 0; n < 12; ++n) {
            const double x = static_cast<double>(n + 1);
            CHECK(run(x) == ((n >= 4) ? static_cast<double>(n - 3) : 0.0));
        }
    }

    TEST_CASE("a first-order TF is one first-order section") {
        const TransferFunction<2, 2, double> tf{.num = {1.0, -0.5}, .den = {1.0, -0.9}};
        const auto                           sos = design::second_order_sections(tf);
        REQUIRE(sos.success);
        CHECK(sos.sections[0].a2 == 0.0);
        CHECK(sos.sections[0].b2 == 0.0);
        CHECK(sos.sections[0].a1 == doctest::Approx(-0.9));
        CHECK(sos.sections[0].b1 / sos.sections[0].b0 == doctest::Approx(-0.5));
        const auto want = direct_impulse<40>(tf);
        const auto got = cascade_impulse<40>(sos);
        for (size_t n = 0; n < 40; ++n) {
            CHECK(got[n] == doctest::Approx(want[n]).epsilon(1e-12).scale(1.0));
        }
        const auto from_zpk = design::second_order_sections(ZPK<1, 1, double>{
            .zeros = {damp::complex<double>{0.5, 0.0}},
            .poles = {damp::complex<double>{0.9, 0.0}},
            .gain = 1.0,
        });
        REQUIRE(from_zpk.success);
        CHECK(from_zpk.sections[0].a1 == doctest::Approx(-0.9));
        CHECK(from_zpk.sections[0].b1 / from_zpk.sections[0].b0 == doctest::Approx(-0.5));
    }

    TEST_CASE("an odd order puts one first-order section first") {
        const ZPK<3, 3, double> zpk{
            .zeros = {damp::complex<double>{-1.0, 0.0}, damp::complex<double>{-1.0, 0.0}, damp::complex<double>{-1.0, 0.0}},
            .poles = {damp::complex<double>{0.2, 0.0}, damp::complex<double>{0.7, 0.4}, damp::complex<double>{0.7, -0.4}},
            .gain = 0.1,
        };
        const auto sos = design::second_order_sections(zpk.to_transfer_function());
        REQUIRE(sos.success);
        CHECK(sos.sections[0].a2 == 0.0);
        CHECK(sos.sections[0].b2 == 0.0);
        CHECK(sos.sections[0].a1 == doctest::Approx(-0.2));
        CHECK(sos.sections[1].a2 == doctest::Approx((0.7 * 0.7) + (0.4 * 0.4)));
    }

    TEST_CASE("refusals") {
        using design::DesignStatus;
        // The denominator's leading zero leaves more zeros than poles: not causal
        CHECK(design::second_order_sections(TransferFunction<3, 3, double>{.num = {1.0, 2.0, 3.0}, .den = {0.0, 1.0, 0.5}}).status == DesignStatus::invalid_parameter);
        CHECK(design::second_order_sections(TransferFunction<2, 2, double>{.num = {1.0, 0.0}, .den = {0.0, 0.0}}).status == DesignStatus::invalid_parameter); // zero denominator
        CHECK(design::second_order_sections(TransferFunction<2, 2, double>{.num = {std::nan(""), 0.0}, .den = {1.0, 0.5}}).status == DesignStatus::invalid_parameter);
        CHECK(design::lowpass_sections<4>(600.0, Ts).status == DesignStatus::invalid_parameter); // above Nyquist
        CHECK(design::lowpass_sections<4>(50.0, 0.0).status == DesignStatus::invalid_sample_time);
        CHECK(design::highpass_sections<3>(-1.0, Ts).status == DesignStatus::invalid_parameter);

        BiquadCascade run{design::lowpass_sections<4>(600.0, Ts)}; // a failed design outputs 0 at run time
        CHECK(run(1.0) == 0.0);
        CHECK(run(1.0) == 0.0);
    }

    TEST_CASE("designs and cascades are constexpr; a cascade deduces from its design") {
        static constexpr auto sos = design::lowpass_sections<6>(100.0, 1e-3);
        static_assert(sos.success);
        static constinit BiquadCascade<3, float> lpf{sos};
        static constexpr auto                    factored = design::second_order_sections(kZpk.to_transfer_function());
        static_assert(factored.success);
        BiquadCascade deduced{factored};
        static_assert(std::is_same_v<decltype(deduced), BiquadCascade<3, double>>);
        float y = 0.0f;
        for (int n = 0; n < 2000; ++n) {
            y = lpf(1.0f);
        }
        CHECK(y == doctest::Approx(1.0f).epsilon(1e-4)); // unity DC gain in float
    }
}
