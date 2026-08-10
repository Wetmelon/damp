// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>

#include "damp/backend.hpp"
#include "damp/estimation/excitation.hpp"
#include "damp/estimation/frequency_response.hpp"
#include "damp/filters/spectral.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

// Discrete FRF of a SISO plant at frequency f [Hz]: G(z) with z = e^{j 2π f Ts}.
template<size_t NX>
void discrete_siso_frf(
    const StateSpace<NX, 1, 1, double>& sys_d, double f_hz, double& mag, double& phase
) {
    const double                w = 2.0 * damp::numbers::pi_v<double> * f_hz;
    const double                th = w * sys_d.Ts;
    const damp::complex<double> z{std::cos(th), std::sin(th)};
    const auto                  G = *eval_frf(sys_d, z);
    mag = damp::abs(G(0, 0));
    phase = damp::arg(G(0, 0));
}

// Continuous second-order: G(s) = wn² / (s² + 2 ζ wn s + wn²) — for mode-extractor table only.
void second_order_analytic(
    double wn, double zeta, double f_hz, double& mag, double& phase
) {
    const double w = 2.0 * damp::numbers::pi_v<double> * f_hz;
    const double re = (wn * wn) - (w * w);
    const double im = 2.0 * zeta * wn * w;
    const double den = (re * re) + (im * im);
    const double gr = (wn * wn) * re / den;
    const double gi = -(wn * wn) * im / den;
    mag = std::hypot(gr, gi);
    phase = std::atan2(gi, gr);
}

} // namespace

TEST_SUITE("Frequency response estimator") {
    TEST_CASE("coherent_block_size rounds cycles * fs / f") {
        CHECK(FrequencyResponseEstimator<4, double>::coherent_block_size(50.0, 1000.0, 10) == 200);
        CHECK(FrequencyResponseEstimator<4, double>::coherent_block_size(25.0, 1000.0, 4) == 160);
    }

    TEST_CASE("H1 lock-in matches analytic first-order plant") {
        // Plant: G = 2 / (0.05 s + 1). fs = 1 kHz, pure sine in → discrete ZOH sim.
        constexpr double Ts = 0.001;
        constexpr double fs = 1.0 / Ts;
        constexpr double K = 2.0;
        constexpr double tau = 0.05;

        TransferFunction<1, 2, double> tf{
            .num = {K},
            .den = {1.0, tau},
        };
        const auto sys_c = tf.to_state_space().value();
        const auto sys_d = *discretize(sys_c, Ts, DiscretizationMethod::ZOH);

        constexpr std::size_t   Nf = 5;
        damp::array<double, Nf> freqs{5.0, 10.0, 20.0, 40.0, 80.0};

        FrequencyResponseEstimator<Nf, double> frf(fs);
        constexpr std::size_t                  measure_cycles = 20;
        constexpr std::size_t                  settle_cycles = 10;
        constexpr double                       A = 1.0;

        for (std::size_t i = 0; i < Nf; ++i) {
            const double      f = freqs[i];
            const auto        N = frf.coherent_block_size(f, measure_cycles);
            const auto        Nset = frf.coherent_block_size(f, settle_cycles);
            ColVec<1, double> x{};

            frf.begin_bin(i, f, N, /*num_blocks=*/2, /*settle=*/Nset);
            // Drive until the bin completes (settle + 2 blocks).
            for (int guard = 0; guard < 100000; ++guard) {
                const double t = static_cast<double>(guard) * Ts;
                const double u = A * std::sin(2.0 * damp::numbers::pi_v<double> * f * t);
                // y = C x + D u; x+ = A x + B u
                const auto   yv = sys_d.C * x + sys_d.D * ColVec<1, double>{u};
                const double y = yv(0, 0);
                x = sys_d.A * x + sys_d.B * ColVec<1, double>{u};
                if (frf.push(u, y)) {
                    break;
                }
            }

            REQUIRE(frf.point(i).valid);
            double mag_a = 0.0;
            double ph_a = 0.0;
            discrete_siso_frf(sys_d, f, mag_a, ph_a);

            CHECK(frf.point(i).magnitude == doctest::Approx(mag_a).epsilon(0.03));
            CHECK(frf.point(i).phase_rad == doctest::Approx(ph_a).epsilon(0.05));
            CHECK(frf.point(i).coherence == doctest::Approx(1.0).epsilon(1e-6));
        }
    }

    TEST_CASE("push_bank recovers multi-sine FRF bins") {
        constexpr double               Ts = 0.001;
        constexpr double               fs = 1000.0;
        TransferFunction<1, 2, double> tf{.num = {1.0}, .den = {1.0, 0.01}};
        const auto                     sys_d = *discretize(*tf.to_state_space(), Ts, DiscretizationMethod::ZOH);

        constexpr std::size_t   Nf = 3;
        damp::array<double, Nf> freqs{10.0, 20.0, 40.0};
        // Coherent for all three: N = 200 → 2, 4, 8 cycles at 10/20/40 Hz.
        constexpr std::size_t N = 200;

        FrequencyResponseEstimator<Nf, double> frf(fs);
        for (std::size_t i = 0; i < Nf; ++i) {
            frf.configure_bin(i, freqs[i], N);
        }

        // Discard plant transient over a full block, then measure one block.
        ColVec<1, double> x{};
        auto              drive = [&](std::size_t n_samples, bool measure) {
            bool done = false;
            for (std::size_t k = 0; k < n_samples; ++k) {
                const double t = static_cast<double>(k) * Ts;
                double       u = 0.0;
                for (std::size_t i = 0; i < Nf; ++i) {
                    u += std::sin(2.0 * damp::numbers::pi_v<double> * freqs[i] * t);
                }
                const auto   yv = sys_d.C * x + sys_d.D * ColVec<1, double>{u};
                const double y = yv(0, 0);
                x = sys_d.A * x + sys_d.B * ColVec<1, double>{u};
                if (measure) {
                    done = frf.push_bank(u, y);
                }
            }
            return done;
        };
        drive(N, false);
        // Restart bank Goertzels after settle (configure again clears accumulators).
        for (std::size_t i = 0; i < Nf; ++i) {
            frf.configure_bin(i, freqs[i], N);
        }
        const bool done = drive(N, true);
        REQUIRE(done);
        CHECK(frf.valid_count() == Nf);
        for (std::size_t i = 0; i < Nf; ++i) {
            double mag_a = 0.0;
            double ph_a = 0.0;
            discrete_siso_frf(sys_d, freqs[i], mag_a, ph_a);
            CHECK(frf.point(i).valid);
            CHECK(frf.point(i).magnitude == doctest::Approx(mag_a).epsilon(0.05));
            CHECK(frf.point(i).phase_rad == doctest::Approx(ph_a).epsilon(0.08));
            (void)ph_a;
        }
    }

    TEST_CASE("SteppedSine walks the frequency table") {
        constexpr design::SteppedSineConfig<3, double> cfg{
            .frequencies_hz = {10.0, 20.0, 40.0},
            .amplitude = 1.5,
            .cycles_per_freq = 2,
        };
        static_assert(cfg.valid());
        constexpr auto design = design::stepped_sine(cfg);
        static_assert(design.success);

        SteppedSine<3, double> exc(design, 0.001);
        CHECK(exc.frequency_index() == 0);
        CHECK(exc.frequency_hz() == doctest::Approx(10.0));

        std::size_t last_idx = 0;
        int         steps = 0;
        while (!exc.done() && steps < 100000) {
            (void)exc.step();
            if (exc.frequency_index() != last_idx) {
                last_idx = exc.frequency_index();
            }
            ++steps;
        }
        CHECK(exc.done());
        CHECK(last_idx == 2);
        CHECK(steps > 100);
    }
}

TEST_SUITE("Mode extractor") {
    TEST_CASE("recovers underdamped second-order resonance from analytic FRF") {
        // G = wn² / (s² + 2ζwn s + wn²), fn = 25 Hz, ζ = 0.05
        constexpr double fn = 25.0;
        constexpr double wn = 2.0 * damp::numbers::pi_v<double> * fn;
        constexpr double zeta = 0.05;

        constexpr std::size_t             Nf = 41;
        damp::array<FrfPoint<double>, Nf> table{};
        // Dense log-ish grid around the resonance: 5..80 Hz
        for (std::size_t i = 0; i < Nf; ++i) {
            const double f = 5.0 + static_cast<double>(i) * (75.0 / static_cast<double>(Nf - 1));
            double       mag = 0.0;
            double       ph = 0.0;
            second_order_analytic(wn, zeta, f, mag, ph);
            table[i] = FrfPoint<double>{f, mag, ph, 1.0, true};
        }

        const auto modes = design::extract_modes<Nf, 2>(table);
        REQUIRE(modes.success);
        REQUIRE(modes.count >= 1);
        CHECK(modes.modes[0].frequency_hz == doctest::Approx(fn).epsilon(0.08));
        // Half-power on a discrete grid is approximate for light damping.
        CHECK(modes.modes[0].zeta == doctest::Approx(zeta).epsilon(0.35));
        CHECK(modes.modes[0].valid);
        CHECK(modes.modes[0].peak_gain > 1.0);
    }

    TEST_CASE("extract_modes is constexpr on a tiny table") {
        constexpr bool ok = []() consteval {
            damp::array<FrfPoint<double>, 5> table{};
            // Synthetic triangle peak at 30 Hz.
            table[0] = {10.0, 0.5, 0.0, 1.0, true};
            table[1] = {20.0, 1.0, 0.0, 1.0, true};
            table[2] = {30.0, 2.0, 0.0, 1.0, true};
            table[3] = {40.0, 1.0, 0.0, 1.0, true};
            table[4] = {50.0, 0.5, 0.0, 1.0, true};
            ModeExtractorConfig<double> cfg{};
            cfg.require_valid_zeta = false; // coarse 5-point table
            const auto modes = design::extract_modes<5, 1>(table, cfg);
            return modes.success && modes.modes[0].frequency_hz == 30.0
                && modes.modes[0].kind == FrfFeatureKind::Peak;
        }();
        static_assert(ok);
        CHECK(ok);
    }

    TEST_CASE("extract_valleys finds a synthetic anti-resonance") {
        damp::array<FrfPoint<double>, 7> table{};
        // Magnitude dips at 40 Hz.
        const double freqs[] = {10, 20, 30, 40, 50, 60, 70};
        const double mags[] = {1.0, 0.9, 0.5, 0.1, 0.5, 0.9, 1.0};
        for (std::size_t i = 0; i < 7; ++i) {
            table[i] = {freqs[i], mags[i], 0.0, 1.0, true};
        }
        ModeExtractorConfig<double> cfg{};
        cfg.require_valid_zeta = false;
        cfg.min_prominence = 0.05;
        const auto valleys = design::extract_valleys<7, 2>(table, cfg);
        REQUIRE(valleys.success);
        CHECK(valleys.modes[0].frequency_hz == doctest::Approx(40.0));
        CHECK(valleys.modes[0].kind == FrfFeatureKind::Valley);
    }

    TEST_CASE("margins_from_frf recovers PM of a first-order loop shape") {
        // L = K / (j w tau) style table: |L|=K/(w tau), phase = -90°
        // Gain crossover when K/(w tau)=1 → w = K/tau
        constexpr double                  K = 10.0;
        constexpr double                  tau = 0.1;
        constexpr std::size_t             Nf = 20;
        damp::array<FrfPoint<double>, Nf> table{};
        for (std::size_t i = 0; i < Nf; ++i) {
            const double f = 1.0 + static_cast<double>(i) * 2.0; // 1..39 Hz
            const double w = 2.0 * damp::numbers::pi_v<double> * f;
            const double mag = K / (w * tau);
            table[i] = {f, mag, -damp::numbers::pi_v<double> / 2.0, 1.0, true};
        }
        const auto m = design::margins_from_frf(table);
        REQUIRE(m.has_phase_margin);
        // PM = 180 + (-90) = 90°
        CHECK(m.phase_margin_deg == doctest::Approx(90.0).epsilon(0.05));
        const double f_c_expect = (K / tau) / (2.0 * damp::numbers::pi_v<double>);
        CHECK(m.gain_crossover_hz == doctest::Approx(f_c_expect).epsilon(0.15));
    }
}

TEST_SUITE("Goertzel FRF accessors") {
    TEST_CASE("real/imag match magnitude/phase") {
        const double      fs = 1000.0;
        const double      f = 50.0;
        const std::size_t N = 200;
        Goertzel<double>  g(f, fs, N);
        for (std::size_t k = 0; k < N; ++k) {
            g.push(std::cos(2.0 * damp::numbers::pi_v<double> * f * static_cast<double>(k) / fs));
        }
        REQUIRE(g.complete());
        CHECK(std::hypot(g.real(), g.imag()) == doctest::Approx(g.magnitude()).epsilon(1e-12));
        CHECK(std::atan2(g.imag(), g.real()) == doctest::Approx(g.phase()).epsilon(1e-12));
    }
}
