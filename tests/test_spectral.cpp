// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>
#include <utility>

#include "damp/filters/spectral.hpp"
#include "damp/math/math.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {
constexpr double pi = 3.14159265358979323846;

// Feed `n` samples of `gen(k)` into a streaming spectral block; returns the
// block-complete flag from the final sample.
template<typename Block, typename Gen>
bool feed(Block& b, size_t n, Gen gen) {
    bool done = false;
    for (size_t k = 0; k < n; ++k) {
        done = b.push(gen(k));
    }
    return done;
}
} // namespace

TEST_SUITE("spectral") {
    TEST_CASE("Goertzel recovers a sinusoid's amplitude on a coherent block") {
        // fs = 1000 Hz, f = 50 Hz, N = 200 -> exactly 10 cycles (coherent).
        const double fs = 1000.0;
        const double f = 50.0;
        const size_t N = 200;
        const double A = 2.0;

        Goertzel<double> g(f, fs, N);
        const bool       done = feed(g, N, [&](size_t k) {
            return A * std::cos(2.0 * pi * f * static_cast<double>(k) / fs);
        });

        CHECK(done);
        CHECK(g.complete());
        CHECK(g.amplitude() == doctest::Approx(A).epsilon(1e-9));
        CHECK(g.magnitude() == doctest::Approx(A * N / 2.0).epsilon(1e-9));
    }

    TEST_CASE("Goertzel rejects a frequency that is not its bin") {
        const double fs = 1000.0;
        const size_t N = 200;
        // Detector tuned to 50 Hz, signal is a (coherent) 70 Hz tone -> ~0.
        Goertzel<double> g(50.0, fs, N);
        feed(g, N, [&](size_t k) { return std::cos(2.0 * pi * 70.0 * static_cast<double>(k) / fs); });
        CHECK(g.amplitude() < 1e-9);
    }

    TEST_CASE("Goertzel streams back-to-back blocks (auto-restart)") {
        const double     fs = 1000.0;
        const double     f = 50.0;
        const size_t     N = 200;
        Goertzel<double> g(f, fs, N);
        // Two consecutive blocks of the same tone both report the amplitude.
        for (int block = 0; block < 2; ++block) {
            const bool done = feed(g, N, [&](size_t k) {
                return 1.5 * std::cos(2.0 * pi * f * static_cast<double>(k) / fs);
            });
            CHECK(done);
            CHECK(g.amplitude() == doctest::Approx(1.5).epsilon(1e-9));
        }
    }

    TEST_CASE("HarmonicAnalyzer separates harmonics and computes THD") {
        // f0 = 50 Hz, fs = 1000 Hz, N = 200 (coherent for f0 and all harmonics).
        const double fs = 1000.0;
        const double f0 = 50.0;
        const size_t N = 200;

        // x = 1.0·H1 + 0.3·H3 + 0.1·H5  (odd harmonics, typical of a power stage).
        auto gen = [&](size_t k) {
            const double t = static_cast<double>(k) / fs;
            return 1.0 * std::cos(2.0 * pi * (1 * f0) * t) + 0.3 * std::cos(2.0 * pi * (3 * f0) * t)
                 + 0.1 * std::cos(2.0 * pi * (5 * f0) * t);
        };

        HarmonicAnalyzer<5, double> ha(f0, fs, N); // tracks 50,100,150,200,250 Hz
        const bool                  done = feed(ha, N, gen);
        REQUIRE(done);

        CHECK(ha.amplitude(1) == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(ha.amplitude(2) == doctest::Approx(0.0).epsilon(1e-9));
        CHECK(ha.amplitude(3) == doctest::Approx(0.3).epsilon(1e-9));
        CHECK(ha.amplitude(4) == doctest::Approx(0.0).epsilon(1e-9));
        CHECK(ha.amplitude(5) == doctest::Approx(0.1).epsilon(1e-9));

        // THD = sqrt(0.3² + 0.1²) / 1.0 = sqrt(0.10) ≈ 0.31623
        CHECK(ha.thd() == doctest::Approx(std::sqrt(0.10)).epsilon(1e-9));
        // Fundamental RMS = 1/√2; total harmonic RMS = sqrt(0.10)/√2.
        CHECK(ha.rms() == doctest::Approx(1.0 / std::sqrt(2.0)).epsilon(1e-9));
        CHECK(ha.total_harmonic_rms() == doctest::Approx(std::sqrt(0.10) / std::sqrt(2.0)).epsilon(1e-9));
    }

    TEST_CASE("HarmonicAnalyzer reports zero THD for a clean fundamental") {
        const double                fs = 1000.0;
        const double                f0 = 50.0;
        const size_t                N = 200;
        HarmonicAnalyzer<5, double> ha(f0, fs, N);
        feed(ha, N, [&](size_t k) { return std::cos(2.0 * pi * f0 * static_cast<double>(k) / fs); });
        CHECK(ha.thd() == doctest::Approx(0.0).epsilon(1e-9));
    }

    TEST_CASE("float specialization") {
        const float     fs = 1000.0f;
        const float     f = 50.0f;
        const size_t    N = 200;
        Goertzel<float> g(f, fs, N);
        feed(g, N, [&](size_t k) { return 1.0f * std::cos(2.0f * 3.14159265f * f * static_cast<float>(k) / fs); });
        CHECK(g.amplitude() == doctest::Approx(1.0f).epsilon(1e-4));
    }

    TEST_CASE("Goertzel is constexpr-evaluable") {
        constexpr bool ok = []() consteval {
            Goertzel<double> g(50.0, 1000.0, 200);
            bool             done = false;
            for (size_t k = 0; k < 200; ++k) {
                done = g.push(damp::cos(2.0 * 3.14159265358979323846 * 50.0 * static_cast<double>(k) / 1000.0));
            }
            return done && damp::abs(g.amplitude() - 1.0) < 1e-9;
        }();
        static_assert(ok, "Goertzel must work at compile time");
        CHECK(ok);
    }

    TEST_CASE("Goertzel rejects invalid block_size / fs / freq") {
        Goertzel<double> zero_n(50.0, 1000.0, 0);
        CHECK_FALSE(zero_n.push(1.0));
        CHECK(zero_n.amplitude() == doctest::Approx(0.0));
        CHECK(zero_n.block_size() == 0);

        Goertzel<double> bad_fs(50.0, 0.0, 200);
        CHECK_FALSE(bad_fs.push(1.0));
        CHECK(bad_fs.amplitude() == doctest::Approx(0.0));

        Goertzel<double> above_nyquist(600.0, 1000.0, 200);
        CHECK_FALSE(above_nyquist.push(1.0));
        CHECK(above_nyquist.amplitude() == doctest::Approx(0.0));
    }

    TEST_CASE("HarmonicAnalyzer amplitude/phase out-of-range are zero") {
        HarmonicAnalyzer<3, double> ha(50.0, 1000.0, 200);
        feed(ha, 200, [&](size_t k) { return std::cos(2.0 * pi * 50.0 * static_cast<double>(k) / 1000.0); });
        CHECK(ha.amplitude(0) == doctest::Approx(0.0));
        CHECK(ha.amplitude(4) == doctest::Approx(0.0));
        CHECK(ha.phase(0) == doctest::Approx(0.0));
        CHECK(ha.phase(99) == doctest::Approx(0.0));
        CHECK(ha.amplitude(1) == doctest::Approx(1.0).epsilon(1e-9));
    }
}
