// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>

#include "damp/backend.hpp"
#include "damp/encoder/axis_position.hpp"
#include "damp/encoder/quadrature.hpp"
#include "damp/encoder/tracker.hpp"
#include "damp/encoder/unwrap.hpp"
#include "damp/math/split_real.hpp"
#include "damp/matrix/colvec.hpp"
#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {
// Quadrature cycles in (A, B), each starting and ending at state 00 so repeated
// cycles chain cleanly from the decoder's initial 00 state.
// Forward: 00 -> 01 -> 11 -> 10 -> 00.
constexpr std::array<damp::pair<bool, bool>, 4> kForward{{{false, true}, {true, true}, {true, false}, {false, false}}};
// Reverse: 00 -> 10 -> 11 -> 01 -> 00.
constexpr std::array<damp::pair<bool, bool>, 4> kReverse{{{true, false}, {true, true}, {false, true}, {false, false}}};

template<typename Dec>
constexpr void drive(Dec& d, const std::array<damp::pair<bool, bool>, 4>& seq, int cycles) {
    for (int c = 0; c < cycles; ++c) {
        for (const auto& s : seq) {
            d.update(s.first, s.second);
        }
    }
}
} // namespace

TEST_SUITE("Quadrature encoder & speed") {
    TEST_CASE("X4 counts every edge, both directions") {
        QuadratureDecoder fwd{QuadMode::X4};
        drive(fwd, kForward, 3);
        CHECK(fwd.position() == 12); // 4 edges * 3 cycles

        QuadratureDecoder rev{QuadMode::X4};
        drive(rev, kReverse, 3);
        CHECK(rev.position() == -12);
    }

    TEST_CASE("X1 counts one detent per cycle") {
        QuadratureDecoder x1{QuadMode::X1};
        drive(x1, kForward, 5);
        CHECK(x1.position() == 5);
    }

    TEST_CASE("X2 counts two edges per cycle") {
        QuadratureDecoder x2{QuadMode::X2};
        drive(x2, kForward, 5);
        CHECK(x2.position() == 10);
    }

    TEST_CASE("illegal (double-bit) transition is ignored") {
        QuadratureDecoder d{QuadMode::X4};
        d.update(false, false); // state 00
        d.update(true, true);   // 00 -> 11 : illegal, no count
        CHECK(d.position() == 0);
    }

    TEST_CASE("index resets position to zero") {
        QuadratureDecoder d{QuadMode::X4};
        drive(d, kForward, 2);
        CHECK(d.position() == 8);
        d.index(true);
        CHECK(d.position() == 0);
    }

    TEST_CASE("wrapped_delta handles counter rollover") {
        CHECK(wrapped_delta<uint16_t>(0xFFFE, 0x0002) == 4);  // wrap up
        CHECK(wrapped_delta<uint16_t>(0x0002, 0xFFFE) == -4); // wrap down
        CHECK(wrapped_delta<uint16_t>(100, 150) == 50);
    }

    TEST_CASE("PhaseUnwrapper accumulates continuously across the wrap seam") {
        PhaseUnwrapper<double> unwrap;                              // phase in turns, [−0.5, 0.5)
        CHECK(unwrap.update(0.4).value() == doctest::Approx(0.4));  // first sample: from 0
        CHECK(unwrap.update(-0.4).value() == doctest::Approx(0.6)); // +0.2 forward across the seam, not −0.8
        CHECK(unwrap.update(0.4).value() == doctest::Approx(0.4));  // −0.2 back across the seam
        CHECK(unwrap.update(0.3).value() == doctest::Approx(0.3));  // plain −0.1, no wrap
    }

    TEST_CASE("PhaseUnwrapper takes radians as a fraction of a turn, several turns") {
        PhaseUnwrapper<double> unwrap;
        const double           two_pi = 2.0 * std::numbers::pi_v<double>;
        const double           step = 2.0; // rad/sample, under half a turn (π)
        for (int i = 1; i <= 10; ++i) {
            (void)unwrap.update(std::remainder(static_cast<double>(i) * step, two_pi) / two_pi);
        }
        CHECK(unwrap.position().value() * two_pi == doctest::Approx(10.0 * step)); // 20 rad, several wraps
    }

    TEST_CASE("QuadratureDecoder X1 does not drift while dithering on an A edge") {
        QuadratureDecoder dec{QuadMode::X1};
        dec.update(false, true); // state 01
        for (int i = 0; i < 10; ++i) {
            dec.update(true, true);  // 01 -> 11 (A rises, forward)
            dec.update(false, true); // 11 -> 01 (A falls, reverse)
        }
        CHECK(dec.position() == 0);
    }
}

TEST_SUITE("EncoderTracker") {
    TEST_CASE("ctor seeds omega from x0 like reset") {
        const ColVec<2, double> x0{{0.25}, {12.5}}; // [counts, counts/s]
        EncoderTracker<double>  obs{2000.0, x0};
        CHECK(obs.position().value() == doctest::Approx(0.25));
        CHECK(obs.omega() == doctest::Approx(12.5));

        // Default x0 leaves omega at 0
        EncoderTracker<double> zeroed{1000.0};
        CHECK(zeroed.omega() == doctest::Approx(0.0));

        // Config path and reset() agree with ctor seed
        EncoderTrackerConfig<double> cfg{.bandwidth = 1500.0, .x0 = ColVec<2, double>{{-0.1}, {-3.0}}};
        EncoderTracker<double>       from_cfg{cfg};
        CHECK(from_cfg.position().value() == doctest::Approx(-0.1));
        CHECK(from_cfg.omega() == doctest::Approx(-3.0));

        zeroed.reset(x0);
        CHECK(zeroed.position().value() == doctest::Approx(0.25));
        CHECK(zeroed.omega() == doctest::Approx(12.5));
    }

    TEST_CASE("predict advances the position with seeded omega") {
        EncoderTracker<double> obs{2000.0, ColVec<2, double>{{0.0}, {10.0}}}; // 10 counts/s
        obs.predict(0.01);                                                    // 0.1 counts
        CHECK(obs.position().value() == doctest::Approx(0.1));
        CHECK(obs.omega() == doctest::Approx(10.0));
    }

    TEST_CASE("update tracks constant speed") {
        // Critically-damped tracker: predict (free run) then update (encoder sample).
        EncoderTracker<double> obs{100.0};
        const double           dt = 0.001;
        const double           omega = 0.5; // counts/s
        double                 true_pos = 0.0;

        for (int k = 0; k < 3000; ++k) {
            true_pos += omega * dt;
            obs.predict(dt);
            obs.update(SplitReal<double>(true_pos), dt);
        }
        CHECK(obs.omega() == doctest::Approx(omega).epsilon(0.05));
        CHECK(obs.position().value() == doctest::Approx(true_pos).epsilon(0.02));
    }
    TEST_CASE("position is the reading exactly; only speed is smoothed") {
        EncoderTracker<double> obs{50.0}; // slow velocity loop: ω is far from settled early on
        const double           dt = 1e-3;
        double                 true_pos = 0.0;
        int                    mismatches = 0;
        for (int k = 0; k < 200; ++k) {
            true_pos += 2.0 * dt; // 2 turns/s
            const auto reading = SplitReal<double>(true_pos);
            obs.predict(dt);
            obs.update(reading, dt);
            mismatches += (obs.position() == reading) ? 0 : 1;
        }
        CHECK(mismatches == 0);
        CHECK(obs.omega() < 2.0); // still converging — and the position never waited for it
    }

    TEST_CASE("extrapolates the last reading by ω between readings") {
        EncoderTracker<double> obs{2000.0, ColVec<2, double>{{0.0}, {10.0}}};
        obs.update(SplitReal<double>(0.25), 1e-4);
        const double w = obs.omega();
        for (int k = 0; k < 4; ++k) { // a 40 kHz current loop between 10 kHz readings
            obs.predict(2.5e-5);
        }
        CHECK((obs.position() - SplitReal<double>(0.25)) == doctest::Approx(w * 1e-4));
    }

    TEST_CASE("speed from quantized counts is much smoother than differencing") {
        constexpr double       cpr = 4096.0;
        const double           dt = 1e-4;
        const double           speed = 10.0; // turns/s
        EncoderTracker<double> obs{300.0, ColVec<2, double>{{0.0}, {speed * cpr}}};
        double                 prev = 0.0;
        double                 worst_diff = 0.0;
        double                 worst_obs = 0.0;
        for (int k = 1; k <= 20000; ++k) {
            const double counts = std::floor(speed * k * dt * cpr + 0.3); // quantized reading
            const double turns = counts / cpr;
            obs.predict(dt);
            obs.update(SplitReal<double>(counts), dt);
            if (k > 5000) { // settled
                worst_diff = std::max(worst_diff, std::abs(((turns - prev) / dt) - speed));
                worst_obs = std::max(worst_obs, std::abs((obs.omega() / cpr) - speed));
            }
            prev = turns;
        }
        CHECK(worst_diff > 1.0);             // one count per sample = 2.4 turns/s of jitter
        CHECK(worst_obs < 0.1 * worst_diff); // the velocity loop removes >90% of it
    }

    TEST_CASE("interpolate_position: the loop's estimate moves within a count, the reading steps") {
        // 100 counts/s read every 0.1 ms: a count every 100 samples. At ω_bw = 50 rad/s the loop
        // filters the reading's staircase (a sawtooth at 628 rad/s) to a few hundredths of a
        // count, so the interpolated position follows the ramp less the −½ count cell offset.
        const double            dt = 1e-4;
        const double            speed = 100.0; // counts/s
        const ColVec<2, double> x0{0.0, speed};
        EncoderTracker<double>  raw{EncoderTrackerConfig<double>{.bandwidth = 50.0, .x0 = x0}};
        EncoderTracker<double>  interp{EncoderTrackerConfig<double>{.bandwidth = 50.0, .x0 = x0, .interpolate_position = true}};
        double                  worst_interp = 0.0;
        double                  worst_step = 0.0;
        double                  prev_interp = 0.0;
        bool                    raw_in_cell = true;
        for (int k = 1; k <= 20000; ++k) {
            const double truth = speed * k * dt;
            const auto   reading = SplitReal<double>(std::floor(truth)); // quadrature: count n is [n, n + 1)
            for (auto* obs : {&raw, &interp}) {
                obs->predict(dt);
                obs->update(reading, dt);
            }
            const double p = interp.position().value();
            if (k > 10000) { // settled (≈ 10/ω_bw)
                worst_interp = std::max(worst_interp, std::abs(p - (truth - 0.5)));
                worst_step = std::max(worst_step, std::abs(p - prev_interp));
                const double e = truth - raw.position().value();
                raw_in_cell = raw_in_cell && e >= 0.0 && e < 1.0;
            }
            prev_interp = p;
        }
        CHECK(raw_in_cell);         // the default is the reading itself
        CHECK(worst_interp < 0.15); // within a sixth of a count of the cell-centred ramp
        CHECK(worst_step < 0.05);   // no one-count steps (the ramp moves 0.01 count per sample)
    }

    TEST_CASE("float tracker in counts: 40M count/rev speed noise is the same everywhere in the turn") {
        // A float angle in turns cannot resolve every count of a 40M count/rev encoder near
        // half a turn, and the speed loop's sub-count corrections round away there: speed noise
        // was ~10x larger near ±½ turn than near 0. Counted in counts it is uniform.
        constexpr std::int64_t cpr = 40'000'000;
        constexpr double       dt = 1.0 / 20000.0;
        constexpr double       rps = 0.3712345678; // not a whole number of counts per sample
        EncoderTracker<float>  obs{EncoderTrackerConfig<float>{.x0 = {0.0F, static_cast<float>(rps * static_cast<double>(cpr))}}};
        double                 sq_near_zero = 0.0;
        double                 sq_near_half = 0.0;
        int                    n_zero = 0;
        int                    n_half = 0;
        for (int k = 0; k < 60000; ++k) {
            const auto c = static_cast<std::int64_t>(std::floor(rps * (k * dt) * static_cast<double>(cpr)));
            obs.predict(static_cast<float>(dt));
            obs.update(SplitReal<float>{c, 0.0F}, static_cast<float>(dt));
            if (k < 10000) {
                continue; // settle
            }
            const double err = static_cast<double>(obs.omega()) - (rps * static_cast<double>(cpr)); // [counts/s]
            const double phase = std::abs(static_cast<double>(AxisPosition<float>{obs.position(), cpr}.phase()));
            if (phase < 0.125) {
                sq_near_zero += err * err;
                ++n_zero;
            } else if (phase > 0.375) {
                sq_near_half += err * err;
                ++n_half;
            }
        }
        REQUIRE(n_zero > 1000);
        REQUIRE(n_half > 1000);
        const double rms_zero = std::sqrt(sq_near_zero / n_zero);
        const double rms_half = std::sqrt(sq_near_half / n_half);
        CHECK(rms_half < 1.5 * rms_zero);
        CHECK(rms_half < 200.0); // the encoder's own quantization floor here is ~73 counts/s
        // The exact count comes back out, near half a turn included.
        CHECK(obs.position().frac == 0.0F);
    }

    TEST_CASE("AxisPosition: phase and distance_to work from the exact count") {
        constexpr std::int64_t    cpr = 40'000'000;
        const AxisPosition<float> p{SplitReal<float>{(7 * cpr) + (cpr / 2) - 1000, 0.0F}, cpr}; // 1000 counts short of 7½ turns
        CHECK(p.phase() == doctest::Approx(0.5F - (1000.0F / static_cast<float>(cpr))).epsilon(1e-7));
        CHECK(p.distance_to(p.counts + 1.0F) == doctest::Approx(1.0F / static_cast<float>(cpr)));
        CHECK(AxisPosition<float>{SplitReal<float>{-3 * cpr, 0.25F}, cpr}.phase() == doctest::Approx(0.25F / static_cast<float>(cpr)));
    }

    TEST_CASE("AxisPosition: the same type for belt travel and a glass scale") {
        // A 23-bit motor encoder on a 40 mm/rev pulley: one period is one turn, 0.040 m of belt.
        constexpr std::int64_t    cpr = std::int64_t{1} << 23;
        const AxisPosition<float> belt{SplitReal<float>{(2 * cpr) + (cpr / 4), 0.0F}, cpr, 0.040F};
        CHECK(belt.value() == doctest::Approx(0.090F)); // 2¼ turns of belt
        CHECK(belt.distance_to(belt.counts + 1.0F) == doctest::Approx(0.040F / static_cast<float>(cpr)));
        CHECK(belt.phase() == doctest::Approx(0.25F)); // a quarter turn in
        // A glass scale with a 20 µm grating interpolated ×400 (50 nm counts), 1.6 m out.
        const AxisPosition<float> scale{SplitReal<float>{32'000'000, 0.0F}, 400, 20e-6F};
        CHECK(scale.distance_to(scale.counts + 3.0F) == doctest::Approx(150e-9F)); // exact 3 counts, far from zero
        CHECK(scale.distance_to(SplitReal<float>{31'999'999, 0.0F}) == doctest::Approx(-50e-9F));
    }

    TEST_CASE("EncoderTracker: default construction matches the default config") {
        EncoderTracker<double> a{};
        EncoderTracker<double> b{EncoderTrackerConfig<double>{}};
        a.update(SplitReal<double>(0.1), 1e-4);
        b.update(SplitReal<double>(0.1), 1e-4);
        CHECK(a.position().value() == doctest::Approx(b.position().value()));
        CHECK(a.omega() == doctest::Approx(b.omega()));
    }

    TEST_CASE("EncoderTrackerConfig keeps the initial state across as<U>()") {
        EncoderTrackerConfig<double> tc{};
        tc.x0 = ColVec<2, double>{0.25, 3.0};
        CHECK(tc.as<float>().x0[0] == doctest::Approx(0.25F));
        CHECK(tc.as<float>().x0[1] == doctest::Approx(3.0F));
    }
}
