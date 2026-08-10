// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>

#include "damp/filters/fir.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

// Sum of coefficients = DC gain for a linear FIR
template<size_t MaxTaps>
double dc_gain(const design::FirDesignResult<MaxTaps, double>& des) {
    double s = 0.0;
    for (size_t i = 0; i < des.n_taps; ++i) {
        s += des.b[i];
    }
    return s;
}

// Alternating sum ≈ Nyquist gain
template<size_t MaxTaps>
double nyquist_gain(const design::FirDesignResult<MaxTaps, double>& des) {
    double s = 0.0;
    for (size_t i = 0; i < des.n_taps; ++i) {
        s += ((i % 2) == 0 ? 1.0 : -1.0) * des.b[i];
    }
    return s;
}

} // namespace

TEST_SUITE("FIR") {
    TEST_CASE("lowpass DC gain ≈ 1") {
        const auto des = design::fir1<64, double>(31, 0.1, design::FirType::Lowpass, design::FirWindow::Hamming);
        REQUIRE(des.success);
        CHECK(des.n_taps == 31);
        CHECK(dc_gain(des) == doctest::Approx(1.0).epsilon(1e-9));
    }

    TEST_CASE("highpass rejects DC and passes Nyquist") {
        const auto des = design::fir1<64, double>(31, 0.15, design::FirType::Highpass, design::FirWindow::Hann);
        REQUIRE(des.success);
        CHECK(des.n_taps % 2 == 1); // Type I
        CHECK(std::abs(dc_gain(des)) < 1e-12);
        // Zero-phase gain at Nyquist is +1; alternating sum equals (−1)^α times that
        // gain, so only the magnitude is checked here.
        CHECK(std::abs(nyquist_gain(des)) == doctest::Approx(1.0).epsilon(1e-6));
    }

    TEST_CASE("impulse response equals coefficients") {
        const auto des = design::fir1<32, double>(9, 0.2, design::FirType::Lowpass, design::FirWindow::Rectangular);
        REQUIRE(des.success);

        FirFilter<32, double> fir(des);
        // Impulse at t=0, then zeros: y[k] must equal b[k]
        for (size_t k = 0; k < des.n_taps; ++k) {
            const double x = (k == 0) ? 1.0 : 0.0;
            const double y = fir(x);
            CHECK(y == doctest::Approx(des.b[k]).epsilon(1e-12));
        }
        // After the impulse fully exits the delay line, output is ~0
        CHECK(std::abs(fir(0.0)) < 1e-12);
    }

    TEST_CASE("runtime lowpass settles to DC input") {
        const auto des = design::fir1<64, double>(21, 0.05, design::FirType::Lowpass);
        REQUIRE(des.success);
        FirFilter<64, double> fir(des);

        double y = 0.0;
        for (int k = 0; k < 200; ++k) {
            y = fir(3.0);
        }
        CHECK(y == doctest::Approx(3.0).epsilon(1e-6));
    }

    TEST_CASE("highpass runtime rejects constant") {
        const auto des = design::fir1<64, double>(33, 0.1, design::FirType::Highpass);
        REQUIRE(des.success);
        FirFilter<64, double> fir(des);

        double y = 0.0;
        for (int k = 0; k < 300; ++k) {
            y = fir(2.5);
        }
        CHECK(std::abs(y) < 1e-6);
    }

    TEST_CASE("bandpass and bandstop design succeed") {
        const auto bp = design::fir1<64, double>(41, 0.1, 0.25, design::FirType::Bandpass, design::FirWindow::Blackman);
        REQUIRE(bp.success);
        CHECK(std::abs(dc_gain(bp)) < 0.15); // stopband at DC (window leakage ok)

        const auto bs = design::fir1<64, double>(41, 0.1, 0.25, design::FirType::Bandstop, design::FirWindow::Hamming);
        REQUIRE(bs.success);
        CHECK(dc_gain(bs) == doctest::Approx(1.0).epsilon(1e-6));
    }

    TEST_CASE("invalid fc → success=false") {
        CHECK_FALSE(design::fir1<32, double>(15, 0.0).success);
        CHECK_FALSE(design::fir1<32, double>(15, 0.5).success);
        CHECK_FALSE(design::fir1<32, double>(15, -0.1).success);
        CHECK_FALSE(design::fir1<32, double>(15, 0.6).success);
        CHECK_FALSE(design::fir1<32, double>(0, 0.1).success);
        CHECK_FALSE(design::fir1<32, double>(64, 0.1).success); // > MaxTaps
        CHECK_FALSE(design::fir1<32, double>(15, 0.3, 0.1, design::FirType::Bandpass).success);
        CHECK_FALSE(design::fir1<32, double>(15, 0.1, design::FirType::Bandpass).success);
    }

    TEST_CASE("windows: rectangular, hann, hamming, blackman") {
        for (auto w : {design::FirWindow::Rectangular, design::FirWindow::Hann, design::FirWindow::Hamming, design::FirWindow::Blackman}) {
            const auto des = design::fir1<32, double>(17, 0.12, design::FirType::Lowpass, w);
            REQUIRE(des.success);
            CHECK(dc_gain(des) == doctest::Approx(1.0).epsilon(1e-9));
        }
    }

    TEST_CASE("highpass even length is adjusted to odd") {
        const auto des = design::fir1<64, double>(20, 0.2, design::FirType::Highpass);
        REQUIRE(des.success);
        CHECK(des.n_taps == 21); // bumped to odd
    }

    TEST_CASE("fir1_hz matches normalized fir1") {
        const auto a = design::fir1<32, double>(15, 100.0 / 1000.0);
        const auto b = design::fir1_hz<32, double>(15, 100.0, 1000.0);
        REQUIRE(a.success);
        REQUIRE(b.success);
        REQUIRE(a.n_taps == b.n_taps);
        for (size_t i = 0; i < a.n_taps; ++i) {
            CHECK(a.b[i] == doctest::Approx(b.b[i]).epsilon(1e-12));
        }
    }

    TEST_CASE("fir_window alias and as<float>") {
        constexpr auto des = design::fir_window<16, double>(9, 0.2);
        static_assert(des.success);
        const auto f = des.as<float>();
        CHECK(f.success);
        CHECK(f.n_taps == 9);
        CHECK(static_cast<double>(f.b[0]) == doctest::Approx(des.b[0]).epsilon(1e-6));
    }

    TEST_CASE("constexpr design smoke") {
        constexpr auto des = design::fir1<16, double>(7, 0.25, design::FirType::Lowpass, design::FirWindow::Rectangular);
        static_assert(des.success);
        static_assert(des.n_taps == 7);
        // DC gain ≈ 1 at compile time
        constexpr double g = des.b[0] + des.b[1] + des.b[2] + des.b[3] + des.b[4] + des.b[5]
                           + des.b[6];
        static_assert(g > 0.999 && g < 1.001);
        CHECK(g == doctest::Approx(1.0).epsilon(1e-9));
    }

    TEST_CASE("FirFilter reset and identity default") {
        FirFilter<8, double> fir;
        CHECK(fir.n_taps() == 1);
        CHECK(fir(4.0) == doctest::Approx(4.0));

        damp::array<double, 8> b{};
        b[0] = 0.5;
        b[1] = 0.5;
        REQUIRE(fir.init(b, 2));
        fir(1.0);
        fir.reset();
        // After reset delay is zero: first sample → 0.5 * x
        CHECK(fir(2.0) == doctest::Approx(1.0));
    }

    TEST_CASE("NLMS identifies a short FIR plant") {
        // Plant: y = 0.5 x[n] + 0.3 x[n-1] + 0.1 x[n-2]
        const double b_true[3] = {0.5, 0.3, 0.1};
        double       x_hist[3] = {0.0, 0.0, 0.0};

        NlmsFilter<8, double> nlms;
        nlms.init(3, 1.0, 1e-12);

        // LCG white ±1 excitation (period 2^31 — persistently exciting for length 3)
        unsigned state = 1u;
        auto     next = [&]() -> double {
            state = state * 1103515245u + 12345u;
            return (state & 0x10000u) ? 1.0 : -1.0;
        };

        double err2 = 0.0;
        for (int k = 0; k < 5000; ++k) {
            const double x = next();
            x_hist[2] = x_hist[1];
            x_hist[1] = x_hist[0];
            x_hist[0] = x;
            const double d = b_true[0] * x_hist[0] + b_true[1] * x_hist[1] + b_true[2] * x_hist[2];
            nlms(x, d);
            if (k >= 4000) {
                err2 += nlms.last_error() * nlms.last_error();
            }
        }

        CHECK(nlms.coefficient(0) == doctest::Approx(b_true[0]).epsilon(0.02));
        CHECK(nlms.coefficient(1) == doctest::Approx(b_true[1]).epsilon(0.02));
        CHECK(nlms.coefficient(2) == doctest::Approx(b_true[2]).epsilon(0.02));
        CHECK(err2 / 1000.0 < 1e-6); // mean-square residual over last 1000 samples
    }
}
