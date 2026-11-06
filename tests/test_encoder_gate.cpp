// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <bit>
#include <cstdint>
#include <type_traits>

#include "damp/encoder/gate.hpp"
#include "damp/encoder/tracker.hpp"
#include "damp/math/split_real.hpp"
#include "doctest.h"

using namespace damp;

namespace {

constexpr double dt = 1e-4;       // [s] encoder frame period
constexpr double speed = 20000.0; // [counts/s] 2 counts per frame

// Built from the bit pattern: under -ffast-math numeric_limits may not give a real NaN / ∞.
const double not_a_number = std::bit_cast<double>(std::uint64_t{0x7ff8000000000000});
const double infinite = std::bit_cast<double>(std::uint64_t{0x7ff0000000000000});

constexpr EncoderGateConfig<double> cfg{.tolerance = 4.0, .max_accel = 1e8, .stale_time = 1e-3, .resync_count = 3};

/// [counts] the true position at frame k, moving at constant speed from 1e9 counts out.
SplitReal<double> truth(int k) {
    return SplitReal<double>{1'000'000'000, 0.0} + (speed * dt * k);
}

/// A gate and tracker locked onto the axis, k frames in.
struct Locked {
    EncoderGate<double>    gate{cfg};
    EncoderTracker<double> tracker{500.0};
    int                    k{0};

    Locked() {
        for (; k < 200; ++k) {
            tracker.predict(dt);
            (void)gate.update(tracker, truth(k), dt);
        }
    }
    EncoderVerdict feed(const SplitReal<double>& reading, bool valid = true) {
        tracker.predict(dt);
        const EncoderVerdict v = gate.update(tracker, reading, dt, valid);
        ++k;
        return v;
    }
};

} // namespace

TEST_SUITE("Encoder gate") {

    TEST_CASE("start-up: locks after resync_count consistent readings, then accepts") {
        EncoderGate<double>    gate{cfg};
        EncoderTracker<double> tracker{500.0};
        CHECK_FALSE(gate.locked());
        CHECK(gate.update(tracker, truth(0), dt) == EncoderVerdict::Jump);
        CHECK(gate.update(tracker, truth(1), dt) == EncoderVerdict::Jump); // 2 counts: inside the window at ω = 0
        CHECK(gate.update(tracker, truth(2), dt) == EncoderVerdict::Resync);
        CHECK(gate.locked());
        CHECK(tracker.position() == truth(2)); // restarted at the confirmed reading, not walked there
        CHECK(gate.update(tracker, truth(3), dt) == EncoderVerdict::Accepted);
        CHECK(tracker.position() == truth(3));

        EncoderGate<double> eager{{.tolerance = 4.0, .resync_count = 1}};
        CHECK(eager.check(truth(0), 0.0, dt) == EncoderVerdict::Resync);
    }

    TEST_CASE("consistent readings pass straight through to the tracker") {
        Locked l;
        REQUIRE(l.gate.locked());
        for (int i = 0; i < 100; ++i) {
            CHECK(l.feed(truth(l.k)) == EncoderVerdict::Accepted);
            CHECK(l.tracker.position() == truth(l.k - 1));
        }
        CHECK(l.tracker.omega() == doctest::Approx(speed).epsilon(1e-3));
        CHECK(l.gate.interval() == doctest::Approx(dt));
    }

    TEST_CASE("CRC failures and non-finite readings are dropped") {
        Locked                  l;
        const SplitReal<double> before = l.tracker.position();
        CHECK(l.feed(truth(l.k), /*valid=*/false) == EncoderVerdict::Invalid);
        CHECK(l.feed(SplitReal<double>{0, not_a_number}) == EncoderVerdict::Invalid);
        CHECK(l.feed(SplitReal<double>{0, infinite}) == EncoderVerdict::Invalid);
        CHECK((l.tracker.position() - before) == doctest::Approx(l.tracker.omega() * 3.0 * dt)); // only extrapolated
        // The next good reading is accepted against a prediction three frames out.
        CHECK(l.feed(truth(l.k)) == EncoderVerdict::Accepted);
        CHECK(l.gate.interval() == doctest::Approx(4 * dt));
        CHECK(l.tracker.position() == truth(l.k - 1));
    }

    TEST_CASE("a glitch is dropped and never reaches the tracker") {
        Locked l;
        for (const double glitch : {1000.0, -37.0, 4.0e6}) {
            const SplitReal<double> last = l.gate.last();
            CHECK(l.feed(truth(l.k) + glitch) == EncoderVerdict::Jump);
            CHECK(l.gate.last() == last);
            CHECK(l.feed(truth(l.k)) == EncoderVerdict::Accepted);
            CHECK(l.tracker.position() == truth(l.k - 1));
        }
        CHECK(l.tracker.omega() == doctest::Approx(speed).epsilon(1e-2)); // speed undisturbed
    }

    TEST_CASE("glitches that do not agree with each other never resync") {
        Locked l;
        for (int i = 0; i < 8; ++i) { // (fewer than stale_time's worth)
            const double glitch = (i % 2 == 0) ? 500.0 : -800.0;
            CHECK(l.feed(truth(l.k) + glitch) == EncoderVerdict::Jump);
        }
        CHECK(l.gate.locked());
        CHECK(l.feed(truth(l.k)) == EncoderVerdict::Accepted);
    }

    TEST_CASE("the window widens with the time since the last accepted reading") {
        // A 30-count deviation is too much one frame out, within reach after dropped frames (½·a·t²).
        Locked l;
        CHECK(l.feed(truth(l.k) + 30.0) == EncoderVerdict::Jump);
        for (int i = 0; i < 6; ++i) {
            CHECK(l.feed(truth(l.k), /*valid=*/false) == EncoderVerdict::Invalid);
        }
        CHECK(l.feed(truth(l.k) + 30.0) == EncoderVerdict::Accepted); // t = 8 frames: ½·1e8·(8e-4)² = 32
    }

    TEST_CASE("a persistent jump resyncs after resync_count readings, keeping the speed") {
        Locked       l;
        const double w = l.tracker.omega();
        const double offset = 5000.0; // the encoder was re-zeroed / slipped
        CHECK(l.feed(truth(l.k) + offset) == EncoderVerdict::Jump);
        CHECK(l.feed(truth(l.k) + offset) == EncoderVerdict::Jump);
        const SplitReal<double> confirmed = truth(l.k) + offset;
        CHECK(l.feed(confirmed) == EncoderVerdict::Resync);
        CHECK(l.tracker.position() == confirmed);
        CHECK(l.tracker.omega() == w);
        CHECK(l.feed(truth(l.k) + offset) == EncoderVerdict::Accepted);
    }

    TEST_CASE("no accepted reading for stale_time raises the stale fault until a resync") {
        Locked l;
        CHECK_FALSE(l.gate.stale());
        for (int i = 0; i < 9; ++i) { // 0.9 ms of failed frames
            (void)l.feed(truth(l.k), /*valid=*/false);
        }
        CHECK_FALSE(l.gate.stale());
        (void)l.feed(truth(l.k), /*valid=*/false);
        (void)l.feed(truth(l.k), /*valid=*/false); // 1.1 ms
        CHECK(l.gate.stale());
        CHECK_FALSE(l.gate.locked());
        // Good frames again: the position must be re-confirmed, not taken from a stale prediction.
        CHECK(l.feed(truth(l.k)) == EncoderVerdict::Jump);
        CHECK(l.feed(truth(l.k)) == EncoderVerdict::Jump);
        CHECK(l.feed(truth(l.k)) == EncoderVerdict::Resync);
        CHECK_FALSE(l.gate.stale());
        CHECK(l.feed(truth(l.k)) == EncoderVerdict::Accepted);

        EncoderGate<double> silent{cfg}; // never locked: stale too
        for (int i = 0; i < 11; ++i) {
            CHECK(silent.check(truth(0), 0.0, dt, /*valid=*/false) == EncoderVerdict::Invalid);
        }
        CHECK(silent.stale());
    }

    TEST_CASE("bad dt and speed inputs never poison the gate") {
        Locked l;
        CHECK(l.gate.check(truth(l.k), not_a_number, dt) == EncoderVerdict::Accepted);
        CHECK(l.gate.check(truth(l.k + 1), speed, -1.0) != EncoderVerdict::Invalid);
        CHECK(l.gate.interval() >= 0.0);
    }

    TEST_CASE("float, as<>, and reset") {
        EncoderGate<float> gate{cfg.as<float>()};
        CHECK(gate.check(SplitReal<float>{1'000'000'000, 0.0F}, 0.0F, 1e-4F) == EncoderVerdict::Jump);
        gate.reset();
        CHECK_FALSE(gate.locked());
        CHECK(cfg.as<float>().resync_count == 3);
    }

    TEST_CASE("EncoderGate<> and EncoderTracker<> compose with their default scalar") {
        // Regression: EncoderTracker defaulted to double while the gate defaulted to float, so
        // gate.update(tracker, ...) did not compile with both written with their defaults.
        EncoderGate<>    gate{};
        EncoderTracker<> tracker{};
        static_assert(std::is_same_v<EncoderTracker<>, EncoderTracker<float>>);
        const EncoderVerdict v = gate.update(tracker, SplitReal<float>(1.0F), 1e-3F);
        CHECK(v != EncoderVerdict::Accepted); // the first reading cannot pass an unlocked gate
    }
}
