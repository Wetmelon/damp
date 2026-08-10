// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <array>
#include <cmath>
#include <complex>
#include <numbers>

#include "damp/filters/filters.hpp"

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
        const auto c = design::highpass_2nd<double>(100.0, Ts);
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
        const auto hp = design::highpass_2nd<double>(100.0, Ts);
        CHECK(std::abs(dc_gain(hp)) < 1e-9);

        // lowpass_2nd: unity DC gain by construction (the lead #17 fix).
        const auto lp = design::lowpass_2nd<double>(100.0, Ts, 0.707); // (fc, Ts, zeta) coeffs overload
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
