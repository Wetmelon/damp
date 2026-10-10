// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Review regressions (group g5): filters, spectral, PLL/SOGI, excitation, analysis.

#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <vector>

#include "damp/analysis/analysis.hpp"
#include "damp/estimation/excitation/chirp.hpp"
#include "damp/estimation/excitation/prbs.hpp"
#include "damp/filters/filters.hpp"
#include "damp/filters/pll.hpp"
#include "damp/filters/sogi.hpp"
#include "damp/filters/spectral.hpp"

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

TEST_SUITE("Review g5 - IIR design") {
    TEST_CASE("First/second-order bilinear designs put the -3 dB corner at fc (prewarped)") {
        const double fc = 250.0;
        const double Ts = 1e-3; // fc = fs/4: unwarped Tustin gave |H(fc)| = 0.618
        const double w = 2.0 * pi * fc * Ts;
        CHECK(gain1(design::lowpass_1st<double>(fc, Ts), w) == doctest::Approx(1.0 / std::sqrt(2.0)).epsilon(1e-9));
        CHECK(gain1(design::highpass_1st<double>(fc, Ts), w) == doctest::Approx(1.0 / std::sqrt(2.0)).epsilon(1e-9));
        CHECK(gain2(design::lowpass_2nd<double>(fc, Ts), w) == doctest::Approx(1.0 / std::sqrt(2.0)).epsilon(1e-9));
        // Same prewarped bilinear as the RBJ cookbook high-pass: |H(fc)| equal.
        CHECK(gain2(design::highpass_2nd<double>(fc, Ts), w) == doctest::Approx(1.0 / std::sqrt(2.0)).epsilon(1e-9));
    }

    TEST_CASE("Discrete Butterworth corner lands on fc") {
        const double fc = 200.0;
        const double Ts = 1e-3;
        const auto   sys = design::butterworth_lowpass<3, double>(fc, Ts);
        REQUIRE(sys.is_discrete());
        const auto b = analysis::bode(sys, std::vector<double>{2.0 * pi * fc});
        REQUIRE(b.points.size() == 1);
        CHECK(b.points[0].magnitude == doctest::Approx(1.0 / std::sqrt(2.0)).epsilon(1e-9));
    }
}

TEST_SUITE("Review g5 - spectral") {
    TEST_CASE("Goertzel phase is the cosine phase of the block (no w(N-1) rotation)") {
        Goertzel<double> g(50.0, 1000.0, 200);
        for (int k = 0; k < 200; ++k) {
            (void)g.push(std::cos(2.0 * pi * 50.0 * k / 1000.0 + 0.5));
        }
        CHECK(g.phase() == doctest::Approx(0.5).epsilon(1e-9));
        CHECK(g.amplitude() == doctest::Approx(1.0).epsilon(1e-9));
    }

    TEST_CASE("Goertzel DC bin amplitude is the DC level") {
        Goertzel<double> g(0.0, 1000.0, 200);
        for (int k = 0; k < 200; ++k) {
            (void)g.push(3.0);
        }
        CHECK(g.amplitude() == doctest::Approx(3.0));
    }

    TEST_CASE("HarmonicAnalyzer completes when its top harmonic is above Nyquist") {
        HarmonicAnalyzer<5, double> ha(50.0, 400.0, 80); // 5*50 Hz > fs/2
        bool                        any = false;
        for (int k = 0; k < 400; ++k) {
            any = ha.push(std::sin(2.0 * pi * 50.0 * k / 400.0)) || any;
        }
        CHECK(any);
        CHECK(ha.amplitude(1) == doctest::Approx(1.0).epsilon(1e-9));
    }

    TEST_CASE("THD is scale invariant") {
        for (const double A : {1.0, 1e-9}) {
            HarmonicAnalyzer<3, double> ha(50.0, 1000.0, 200);
            for (int k = 0; k < 200; ++k) {
                const double th = 2.0 * pi * 50.0 * k / 1000.0;
                (void)ha.push(A * (std::sin(th) + 0.1 * std::sin(3.0 * th)));
            }
            CHECK(ha.thd() == doctest::Approx(0.1).epsilon(1e-6));
        }
    }
}

TEST_SUITE("Review g5 - PLL / SOGI") {
    TEST_CASE("SinglePhasePLL locks phase, not just frequency") {
        const double Ts = 1e-4;
        for (const double ph0 : {0.0, 1.5, -2.5}) {
            SinglePhasePLL<double> pll(50.0);
            double                 th = 0.0;
            for (int k = 0; k < 20000; ++k) {
                th = 2.0 * pi * 50.0 * k * Ts + ph0;
                pll.step(std::cos(th), Ts);
            }
            CHECK(std::cos(pll.phase() - th) > 0.99);
            CHECK(pll.frequency() == doctest::Approx(50.0).epsilon(1e-3));
        }
    }

    TEST_CASE("SogiFll lock is independent of input amplitude") {
        double f[2];
        int    i = 0;
        for (const double A : {1.0, 1e-6}) {
            SogiFll<double> fll(45.0, 1e-4);
            for (int k = 0; k < 2000; ++k) {
                fll.update(A * std::sin(2.0 * pi * 50.0 * k * 1e-4));
            }
            f[i++] = fll.frequency_hz();
        }
        CHECK(f[1] == doctest::Approx(f[0]).epsilon(1e-9));
    }
}

TEST_SUITE("Review g5 - runtime filters") {
    TEST_CASE("FirFilter from a failed design stays identity") {
        FirFilter<4, double> f(design::FirDesignResult<4, double>{});
        CHECK(f(1.0) == doctest::Approx(1.0));
        CHECK(f(2.0) == doctest::Approx(2.0));
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

    TEST_CASE("FIR design normalization is scale-free (float lowpass still succeeds)") {
        const auto des = design::fir1<31, float>(31, 0.1f, design::FirType::Highpass, design::FirWindow::Blackman);
        REQUIRE(des.success);
        float g = 0.0f;
        for (size_t n = 0; n < des.n_taps; ++n) {
            g += des.b[n];
        }
        CHECK(std::abs(g) < 1e-5f); // DC null
    }
}

TEST_SUITE("Review g5 - excitation") {
    TEST_CASE("Log chirp phase matches the closed-form integral") {
        design::ChirpConfig<double> cfg{};
        cfg.f_start_hz = 1.0;
        cfg.f_end_hz = 100.0;
        cfg.duration_s = 2.0;
        cfg.mode = design::ChirpMode::Log;
        const Chirp<double> c(design::chirp(cfg));
        const double        beta = std::log(100.0) / 2.0;
        const double        t = 1.3;
        const double        ph = 2.0 * pi * (std::exp(beta * t) - 1.0) / beta;
        CHECK(c.step(t) == doctest::Approx(std::sin(ph)).epsilon(1e-9));
    }

    TEST_CASE("PRBS step(t) at a huge time stays periodic and in range") {
        design::PRBSConfig<double> cfg{};
        cfg.lfsr_order = 5;
        cfg.clock_period_s = 1.0;
        const PRBS<double> p(design::prbs(cfg), 1.0);
        const double       period = 31.0;
        CHECK(p.step(3.0) == p.step(3.0 + period * 1e6));
        CHECK(std::abs(p.step(1e30)) == doctest::Approx(1.0));
    }

    TEST_CASE("PRBS feedback masks are maximal length") {
        for (std::size_t order = 2; order <= 16; ++order) {
            std::uint32_t s = 1u;
            std::size_t   n = 0;
            do {
                s = damp::detail::prbs_advance(s, order);
                ++n;
            } while (s != 1u && n < (std::size_t{1} << order));
            CHECK(n == (std::size_t{1} << order) - 1);
        }
    }
}

TEST_SUITE("Review g5 - analysis") {
    TEST_CASE("Gain margin finds a +180 deg phase crossover") {
        // L = (s+1)^2 / s^3: phase -270 -> -90, crosses -180 (wrapped +180) at w = 1, |L| = 2.
        TransferFunction<3, 4, double> tf{{1.0, 2.0, 1.0}, {0.0, 0.0, 0.0, 1.0}};
        const auto                     b = analysis::bode(tf, analysis::logspace(0.01, 100.0, 2000));
        const auto                     gm = b.gain_margin();
        REQUIRE(gm.has_value());
        CHECK(gm->first == doctest::Approx(-20.0 * std::log10(2.0)).epsilon(1e-3));
        CHECK(gm->second == doctest::Approx(1.0).epsilon(1e-2));
    }

    TEST_CASE("TF bode is invariant to scaling the denominator") {
        TransferFunction<1, 2, double> tf{{1e-15}, {1e-15, 1e-15}};
        const auto                     b = analysis::bode(tf, std::vector<double>{1.0});
        REQUIRE(b.points.size() == 1);
        CHECK(b.points[0].magnitude == doctest::Approx(1.0 / std::sqrt(2.0)));
    }

    TEST_CASE("Middlebrook from impedance data pairs samples by frequency") {
        analysis::ImpedanceResult<double> zs, zl;
        for (const double w : {1.0, 2.0, 3.0}) {
            zs.points.push_back({w, {1.0, 0.0}, 1.0, 0.0, 0.0});
        }
        // Load sweep lost its w = 2 sample (contour pole skipped).
        zl.points.push_back({1.0, {10.0, 0.0}, 10.0, 20.0, 0.0});
        zl.points.push_back({3.0, {4.0, 0.0}, 4.0, 12.0, 0.0});
        const auto mb = analysis::middlebrook(zs, zl);
        REQUIRE(mb.minor_loop_gain.points.size() == 2);
        CHECK(mb.minor_loop_gain.points[1].omega == doctest::Approx(3.0));
        CHECK(mb.minor_loop_gain.points[1].magnitude == doctest::Approx(0.25));
    }

    TEST_CASE("Discrete H-infinity norm includes DC") {
        StateSpace<1, 1, 1, double> s{
            .A = Matrix<1, 1>{{0.999}},
            .B = Matrix<1, 1>{{0.001}},
            .C = Matrix<1, 1>{{1.0}},
            .D = Matrix<1, 1>{{0.0}},
            .Ts = 0.01,
        };
        CHECK(analysis::norm_hinf(s) == doctest::Approx(1.0).epsilon(1e-9));
    }

    TEST_CASE("pzmap ignores trailing-zero padding of a TF numerator") {
        const auto pz = analysis::pzmap(design::lowpass_2nd_continuous<double>(1.0));
        CHECK(pz.zeros.empty());
        REQUIRE(pz.poles.size() == 2);
        CHECK(pz.poles[0].abs() == doctest::Approx(2.0 * pi));

        TransferFunction<3, 3, double> tf{{2.0, 1.0, 0.0}, {2.0, 3.0, 1.0}};
        const auto                     pz2 = analysis::pzmap(tf);
        REQUIRE(pz2.zeros.size() == 1);
        CHECK(pz2.zeros[0].real() == doctest::Approx(-2.0));
    }

    TEST_CASE("Discrete impulse is the unit-pulse response") {
        StateSpace<1, 1, 1, double> s{
            .A = Matrix<1, 1>{{0.5}},
            .B = Matrix<1, 1>{{1.0}},
            .C = Matrix<1, 1>{{1.0}},
            .D = Matrix<1, 1>{{0.25}},
            .Ts = 1.0,
        };
        const auto r = analysis::impulse(s, analysis::linspace(0.0, 3.0, 4));
        REQUIRE(r.y.size() == 4);
        CHECK(r.y[0](0, 0) == doctest::Approx(0.25)); // D
        CHECK(r.y[1](0, 0) == doctest::Approx(1.0));  // CB
        CHECK(r.y[2](0, 0) == doctest::Approx(0.5));  // CAB
        CHECK(r.y[3](0, 0) == doctest::Approx(0.25));
    }

    TEST_CASE("Non-uniform time grids are rejected") {
        TransferFunction<1, 2, double> tf{{1.0}, {1.0, 1.0}};
        const std::vector<double>      t{0.0, 0.1, 0.2, 0.5, 0.6};
        CHECK(analysis::step(tf, t).t.empty());
        // float linspace with a long horizon still counts as uniform
        TransferFunction<1, 2, float> tff{{1.0f}, {1.0f, 1.0f}};
        CHECK(analysis::step(tff, analysis::linspace(0.0f, 100.0f, 100001)).t.size() == 100001);
    }

    TEST_CASE("stepinfo flags a response that never reaches 90%") {
        const std::vector<double> t{0.0, 1.0, 2.0, 3.0};
        const std::vector<double> y{0.0, 0.2, 0.5, 0.6};
        const auto                info = analysis::stepinfo(y, t, 1.0);
        CHECK_FALSE(info.risen);
        CHECK(info.rise_time == doctest::Approx(0.0));
        CHECK(analysis::stepinfo(y, std::vector<double>{0.0, 1.0}, 1.0).peak == doctest::Approx(0.0));
    }

    TEST_CASE("lsiminfo settles a decay to zero") {
        std::vector<double> t, y;
        for (int k = 0; k <= 1000; ++k) {
            t.push_back(k * 0.01);
            y.push_back(std::exp(-k * 0.01));
        }
        const auto info = analysis::lsiminfo(y, t, 0.0);
        CHECK(info.settling_time == doctest::Approx(std::log(50.0)).epsilon(1e-2));
    }
}
