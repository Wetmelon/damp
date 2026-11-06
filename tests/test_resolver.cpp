// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstdint>
#include <limits>

#include "damp/backend.hpp"
#include "damp/encoder/sincos.hpp"
#include "damp/encoder/tracker.hpp"
#include "damp/encoder/unwrap.hpp"
#include "damp/math/split_real.hpp"
#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {
constexpr double two_pi = 2.0 * damp::numbers::pi_v<double>;
} // namespace

TEST_SUITE("SplitReal") {

    TEST_CASE("splits an exact integer part from a remainder in [-0.5, 0.5)") {
        constexpr auto p = SplitReal<double>(2.75);
        static_assert(p.whole == 3);
        static_assert(p.frac == -0.25);
        static_assert(SplitReal<double>(0.5).whole == 1);
        static_assert(SplitReal<double>(0.5).frac == -0.5);
        static_assert(SplitReal<double>(-0.5).whole == 0);
        static_assert(SplitReal<double>(-3.25).whole == -3);
    }

    TEST_CASE("is affine: position ± displacement, position − position") {
        constexpr auto a = SplitReal<double>(10.4);
        constexpr auto b = a + 0.3; // crosses a cycle boundary
        static_assert(b.whole == 11);
        CHECK(b.frac == doctest::Approx(-0.3));
        CHECK((b - a) == doctest::Approx(0.3));
        CHECK(((b - 0.3) - a) == doctest::Approx(0.0));
    }

    TEST_CASE("keeps float resolution far from zero") {
        // Ten million turns: a plain float position is quantized to a whole turn there.
        const SplitReal<float> base{10'000'000, 0.0F};
        const SplitReal<float> moved = base + 1e-4F;
        CHECK((moved - base) == doctest::Approx(1e-4F).epsilon(1e-3));
        CHECK(static_cast<float>(10'000'000.0001) == 10'000'000.0F); // what a plain float keeps
    }

    TEST_CASE("meters decades of energy at a fixed step without stalling") {
        // A 10 kW inverter integrating P·dt at 20 kHz adds 0.5 J per tick. A float total stops
        // moving once it passes 2^24 · 0.5 J ≈ 8 MJ — about 14 minutes in. Run past that point.
        const float      step = 0.5F;                                        // [J]
        SplitReal<float> energy = SplitReal<float>{9'500'000'000'000, 0.0F}; // ≈ 30 years at 10 kW
        float            naive = 9.5e12F;
        const auto       start = energy;
        for (int i = 0; i < 100'000; ++i) {
            energy += step;
            naive += step;
        }
        CHECK((energy - start) == doctest::Approx(50'000.0F)); // every joule counted
        CHECK(naive == 9.5e12F);                               // a float total never moved
    }

    TEST_CASE("normalizing never moves the value and keeps the remainder in [-0.5, 0.5)") {
        // Odd integer in [2^23, 2^24): frac + 0.5 would tie up to the next even number.
        constexpr auto big = SplitReal<float>(8388609.0F);
        static_assert(big.whole == 8388609);
        static_assert(big.frac == 0.0F);
        // Just below 0.5: frac + 0.5 would round to 1 and move the value by an ulp.
        constexpr auto below_half = SplitReal<float>(0.49999997F);
        static_assert(below_half.whole == 0);
        static_assert(below_half.frac == 0.49999997F);
        constexpr auto below_half_d = SplitReal<double>(0.49999999999999994);
        static_assert(below_half_d.whole == 0);
        static_assert(below_half_d.frac == 0.49999999999999994);
        CHECK(SplitReal<float>(8388609.0F) == SplitReal<float>(8388609, 0.0F));
    }

    TEST_CASE("a negative remainder is kept exactly: small negative steps accumulate") {
        // frac − floor(frac) would round frac + 1 onto the 2⁻²⁴ grid of [0.5, 1).
        constexpr auto tiny_neg = SplitReal<float>(-1e-9F);
        static_assert(tiny_neg.whole == 0);
        static_assert(tiny_neg.frac == -1e-9F);
        constexpr auto neg_third = SplitReal<float>(-0.3F);
        static_assert(neg_third.whole == 0);
        static_assert(neg_third.frac == -0.3F);
        constexpr auto below_neg_half = SplitReal<float>(-0.7F); // folds to −1 + 0.3, exactly
        static_assert(below_neg_half.whole == -1);
        static_assert((below_neg_half - SplitReal<float>{}) == -0.7F);
        SplitReal<float> down{};
        SplitReal<float> up{};
        for (int i = 0; i < 1'000'000; ++i) {
            down -= 3e-8F;
            up += 3e-8F;
        }
        CHECK((down - SplitReal<float>{}) == doctest::Approx(-0.03).epsilon(0.01));
        CHECK((down - SplitReal<float>{}) == -(up - SplitReal<float>{})); // no direction bias
    }

    TEST_CASE("the integer part wraps instead of overflowing past 2^63") {
        constexpr std::int64_t top = std::numeric_limits<std::int64_t>::max();
        constexpr auto         p = SplitReal<double>(top, 1.0); // constant evaluation rejects signed overflow
        static_assert(p.whole == std::numeric_limits<std::int64_t>::min());
        static_assert((p - SplitReal<double>(top, 0.0)) == 1.0);
    }

    TEST_CASE("a non-finite step leaves the whole count alone (no UB cast)") {
        SplitReal<float> p{7, 0.0F};
        p += std::numeric_limits<float>::infinity();
        CHECK(p.whole == 7);
    }
}

TEST_SUITE("Sin/cos decoding") {

    TEST_CASE("decode_sincos: phase is atan2 in periods, amplitude is the trimmed radius") {
        const double theta = 0.7;
        const auto   s = decode_sincos(std::sin(theta), std::cos(theta));
        CHECK(s.frac == doctest::Approx(theta / two_pi));
        CHECK(s.amplitude == doctest::Approx(1.0));
        // atan2's +π end maps to the bottom of the range.
        CHECK(decode_sincos(0.0, -1.0).frac == doctest::Approx(-0.5));
    }

    TEST_CASE("decode_sincos: calibration trims recover the phase from a skewed front end") {
        const SinCosCalibration<double> cal{.sin_gain = 1.0 / 2.0, .cos_gain = 1.0 / 1.5, .sin_offset = 0.3, .cos_offset = -0.1};
        const double                    theta = 0.9;
        const auto                      s = decode_sincos(0.3 + (2.0 * std::sin(theta)), -0.1 + (1.5 * std::cos(theta)), cal);
        CHECK(s.frac == doctest::Approx(theta / two_pi));
        CHECK(s.amplitude == doctest::Approx(1.0));
    }

    TEST_CASE("combine_quadrature: exact period over ±5 periods, with edges a count early or late") {
        int failures = 0;
        for (int k = -2000; k <= 2000; ++k) {
            const double       p = (k * 0.0025) + 0.0003; // [periods]
            const double       frac = p - std::floor(p + 0.5);
            const std::int64_t q = static_cast<std::int64_t>(std::floor(4.0 * p));
            for (std::int64_t jitter = -1; jitter <= 1; ++jitter) {
                const auto pos = combine_quadrature(q + jitter, frac);
                failures += (std::abs(pos.value() - p) > 1e-9) ? 1 : 0;
            }
        }
        CHECK(failures == 0);
    }

    TEST_CASE("combine_quadrature: exact period with edges 0.3 period early or late") {
        // Regression: the rounding was centred on q/4 instead of on q = ⌊4p⌋'s mean lag of
        // 1/8, so a late edge past 1/4 period gave a position one period off.
        int failures = 0;
        for (const double offset : {-0.3, 0.3}) {
            for (int k = -2000; k <= 2000; ++k) {
                const double       p = (k * 0.0025) + 0.0003; // [periods]
                const double       frac = p - std::floor(p + 0.5);
                const std::int64_t q = static_cast<std::int64_t>(std::floor(4.0 * (p - offset)));
                failures += (std::abs(combine_quadrature(q, frac).value() - p) > 1e-9) ? 1 : 0;
            }
        }
        CHECK(failures == 0);
    }

    TEST_CASE("PhaseUnwrapper counts resolver turns across the seam, both directions") {
        PhaseUnwrapper<double> unwrap;
        double                 theta = 0.0;
        for (int i = 0; i < 4000; ++i) { // 0.1 turn per sample, 400 turns forward
            theta += 0.1 * two_pi;
            unwrap.update(decode_sincos(std::sin(theta), std::cos(theta)).frac);
        }
        CHECK(unwrap.position().value() == doctest::Approx(400.0));
        for (int i = 0; i < 1500; ++i) { // and 150 back
            theta -= 0.1 * two_pi;
            unwrap.update(decode_sincos(std::sin(theta), std::cos(theta)).frac);
        }
        CHECK(unwrap.position().value() == doctest::Approx(250.0));
    }
}

TEST_SUITE("EncoderTracker on decoded sin/cos") {

    TEST_CASE("converges to a static resolver angle") {
        EncoderTracker<double> tracker{1000.0}; // bw·dt = 0.1, the stable tracking regime
        PhaseUnwrapper<double> unwrap;
        const double           theta = -1.2;
        for (int i = 0; i < 5000; ++i) {
            tracker.predict(1e-4);
            tracker.update(unwrap.update(decode_sincos(std::sin(theta), std::cos(theta)).frac), 1e-4);
        }
        CHECK(tracker.position().value() * two_pi == doctest::Approx(theta).epsilon(1e-3));
        CHECK(std::abs(tracker.omega()) < 0.1);
    }

    TEST_CASE("tracks a constant-speed ramp: speed and multi-turn position") {
        EncoderTracker<double> tracker{1000.0};
        PhaseUnwrapper<double> unwrap;
        const double           w = 300.0; // rad/s
        const double           dt = 1e-4;
        double                 theta = 0.0;
        for (int i = 0; i < 8000; ++i) { // 0.8 s to settle the loop
            theta += w * dt;
            tracker.predict(dt);
            tracker.update(unwrap.update(decode_sincos(std::sin(theta), std::cos(theta)).frac), dt);
        }
        CHECK(tracker.omega() * two_pi == doctest::Approx(w).epsilon(0.02));
        CHECK(tracker.position().value() == doctest::Approx(theta / two_pi).epsilon(0.02));
    }

    TEST_CASE("absolute updates track a multi-period encoder moving several periods per sample") {
        // 2048-line encoder at 3000 rpm sampled at 10 kHz: 10.24 periods per sample. Shortest-arc
        // unwrapping aliases here; counter + fine phase does not.
        const double           periods_per_s = 2048.0 * 50.0;
        const double           dt = 1e-4;
        EncoderTracker<double> tracker{2000.0};
        double                 p = 0.0; // [periods]
        for (int i = 0; i < 4000; ++i) {
            p += periods_per_s * dt;
            const double phase = p * two_pi;
            const auto   s = decode_sincos(std::sin(phase), std::cos(phase));
            const auto   q = static_cast<std::int64_t>(std::floor(4.0 * p));
            tracker.predict(dt);
            tracker.update(combine_quadrature(q, s.frac), dt);
        }
        CHECK(tracker.omega() == doctest::Approx(periods_per_s).epsilon(1e-3));
        CHECK(std::abs(tracker.position() - SplitReal<double>(p)) < 1e-3);
    }

    TEST_CASE("float tracker keeps sub-micro-turn resolution ten million turns out") {
        EncoderTracker<float> tracker{1000.0F};
        const auto            start = SplitReal<float>{10'000'000, 0.0F};
        tracker.reset(start, 0.0F);
        // Step 1e-4 turns and hold: the tracker must see the step, not round it away.
        for (int i = 0; i < 20000; ++i) {
            tracker.predict(1e-4F);
            tracker.update(start + 1e-4F, 1e-4F);
        }
        CHECK((tracker.position() - start) == doctest::Approx(1e-4F).epsilon(1e-2));
    }
}
