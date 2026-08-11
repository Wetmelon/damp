// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file toolbox/timing.hpp
 * @brief Non-blocking software timers.
 *
 * The plain scheduling idioms firmware actually uses — "how long since…",
 * "has it been N ticks yet", "do this every N ticks" — without blocking.
 * Each is fed the elapsed amount `dt` per call (so they are clock-source
 * agnostic and unit-test friendly), `constexpr`, and allocation-free.
 *
 * The time type `T` defaults to `uint32_t` — the common embedded case of a
 * free-running tick / millisecond counter — but any arithmetic type works; use
 * `float`/`double` for `dt` in seconds.
 *
 * These complement the IEC 61131-3 `TON`/`TOF`/`TP` PLC-scan timers in
 * `iec61131.hpp`; use those for ladder-style logic, these for general scheduling.
 *
 * @see iec61131.hpp for the standard PLC timers.
 */

#include <cstdint>

namespace damp {

/**
 * @brief Free-running elapsed-time accumulator.
 *
 * Accumulates `dt` each @ref tick and reports the total since the last
 * @ref reset. A stopwatch for elapsed-time measurements; reset it when you
 * start timing something. With the default `uint32_t` the total wraps modulo
 * 2³² like a hardware tick counter — reset it, or take wrapped differences, for
 * spans longer than the counter range.
 *
 * @tparam T Time/count type (default: uint32_t)
 */
template<typename T = uint32_t>
class Stopwatch {
public:
    constexpr Stopwatch() = default;

    /// Advance the clock by @p dt; returns the new elapsed total.
    constexpr T tick(T dt) {
        elapsed_ += dt;
        return elapsed_;
    }

    [[nodiscard]] constexpr T elapsed() const { return elapsed_; }

    constexpr void reset(T t0 = T{0}) { elapsed_ = t0; }

private:
    T elapsed_{0};
};

/**
 * @brief One-shot timeout.
 *
 * Counts elapsed time toward a @ref duration and latches @ref expired once it
 * is reached. The accumulator saturates at `duration` (it never grows
 * unbounded). Re-arm with @ref reset.
 *
 * @code
 * Timeout<> comms{50};          // 50 ticks (e.g. ms); or Timeout<float>{0.050f} for seconds
 * // each tick:
 * if (comms.tick(dt)) { enter_comms_fault(); }
 * if (rx_ok) { comms.reset(); }
 * @endcode
 *
 * @tparam T Time/count type (default: uint32_t)
 */
template<typename T = uint32_t>
class Timeout {
public:
    T duration{0}; ///< Time until expiry, in `dt` units (ticks or seconds).

    constexpr Timeout() = default;
    constexpr explicit Timeout(T duration_) : duration(duration_) {}

    /// Advance by @p dt; returns whether the timeout has now expired.
    constexpr bool tick(T dt) {
        if (elapsed_ < duration) {
            elapsed_ += dt;
            if (elapsed_ > duration) {
                elapsed_ = duration; // clamp: bounded accumulator
            }
        }
        return expired();
    }

    [[nodiscard]] constexpr bool expired() const { return elapsed_ >= duration; }
    [[nodiscard]] constexpr T    remaining() const {
        return elapsed_ >= duration ? T{0} : duration - elapsed_;
    }

    constexpr void reset() { elapsed_ = T{0}; }

private:
    T elapsed_{0};
};

/**
 * @brief Periodic trigger — fires once per elapsed @ref period.
 *
 * The "do X every N seconds" idiom: feed it `dt` and it returns `true` on the
 * ticks where a full period has elapsed, carrying the remainder forward so the
 * average rate has no drift. If a single `dt` spans multiple periods the
 * remainder accumulates and it fires again on following calls (gradual
 * catch-up); the stored phase stays bounded by one period. A non-positive
 * @ref period never fires (avoids the `period == 0` infinite-fire trap).
 *
 * @code
 * Periodic<> telemetry{100};       // every 100 ticks; or Periodic<float>{0.1f} for 10 Hz
 * if (telemetry(dt)) { send_telemetry(); }
 * @endcode
 *
 * @tparam T Time/count type (default: uint32_t)
 */
template<typename T = uint32_t>
class Periodic {
public:
    T period{0}; ///< Trigger interval, in `dt` units (ticks or seconds). Must be > 0 to fire.

    constexpr Periodic() = default;
    constexpr explicit Periodic(T period_) : period(period_) {}

    /// Advance by @p dt; returns true once per elapsed period (never when period ≤ 0).
    constexpr bool operator()(T dt) {
        if (!(period > T{0})) {
            return false;
        }
        phase_ += dt;
        if (phase_ >= period) {
            phase_ -= period;
            return true;
        }
        return false;
    }

    constexpr void reset() { phase_ = T{0}; }

private:
    T phase_{0};
};

} // namespace damp
