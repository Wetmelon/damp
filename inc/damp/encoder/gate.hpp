// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file damp/encoder/gate.hpp
 * @brief Encoder reading validation between the raw reading and the tracker.
 *
 * EncoderTracker is the core algorithm and assumes good readings. EncoderGate decides
 * which readings are good: it drops readings the interface flags as bad (CRC, frame
 * error) or that are not finite, drops jumps further from the prediction than the axis
 * can move, raises a stale fault when no reading has been accepted for too long, and
 * re-locks onto a new position once it has held for a few consecutive readings. A
 * dropped reading never reaches the tracker.
 *
 * `constexpr`, allocation-free, embeddable.
 */

#include "damp/backend.hpp"
#include "damp/encoder/tracker.hpp" // EncoderTracker
#include "damp/math/math.hpp"       // damp::abs, damp::isfinite
#include "damp/math/split_real.hpp" // SplitReal

namespace damp {

/** @brief What EncoderGate made of a reading. */
enum class EncoderVerdict : uint8_t {
    Accepted, ///< consistent with the prediction: pass it to the estimator
    Invalid,  ///< flagged bad by the interface (CRC, frame) or not finite: dropped
    Jump,     ///< further from the prediction than the axis can move: dropped
    Resync,   ///< a new position held for EncoderGateConfig::resync_count readings: restart the estimator there
};

/**
 * @brief Configuration for EncoderGate. Positions in counts, as the tracker's.
 * @tparam T Scalar type.
 */
template<typename T = float>
struct EncoderGateConfig {
    /**
     * @brief [counts] deviation from the prediction always allowed: quantization, noise, and the predictor's
     * speed error over a reading interval. ≤ 0 turns the jump check off.
     */
    T tolerance{T{0}};
    /**
     * @brief [counts/s²] fastest the axis can accelerate; the window widens by ½·a·t² with the time since the
     * last accepted reading, so it opens up across dropped readings.
     */
    T max_accel{T{0}};
    /**
     * @brief [s] no accepted reading for longer than this raises the stale fault, and the gate re-locks
     * (≤ 0: never stale).
     */
    T stale_time{static_cast<T>(0.01)};
    /**
     * @brief consecutive mutually consistent readings that confirm a new position (after a persistent jump,
     * at start-up, after staleness). 1 accepts the first reading.
     */
    uint16_t resync_count{4};

    template<typename U>
    [[nodiscard]] constexpr EncoderGateConfig<U> as() const {
        return {static_cast<U>(tolerance), static_cast<U>(max_accel), static_cast<U>(stale_time), resync_count};
    }
};

/**
 * @brief Validates encoder readings before they reach the tracker.
 *
 * Call once per expected reading with the time since the previous call, passing
 * @p valid = false when the frame failed its CRC or did not arrive, so that time still
 * counts towards staleness. Each reading is judged against the prediction
 * @f$ \hat p = p_{last} + \omega\,t @f$ from the last accepted reading, with @f$ t @f$ the
 * time since it and @f$ \omega @f$ the caller's speed estimate:
 * @f[
 *   |p - \hat p| \le \text{tolerance} + \tfrac12\,a_{max}\,t^2 .
 * @f]
 * - Accepted: inside the window. Feed it to the estimator (for EncoderTracker, with
 *   @ref interval as its dt).
 * - Invalid / Jump: dropped. A jump starts (or extends) a candidate track: readings that
 *   stay consistent with each other under the same window and speed.
 * - Resync: a candidate track reached EncoderGateConfig::resync_count readings: the
 *   position really moved (a slipped count, a re-zeroed encoder). Restart the estimator at
 *   the reading, keeping its speed. A position offset leaves the speed alone, which is
 *   why the candidates are predicted with the same speed.
 *
 * The gate starts unlocked and locks with a Resync, so the first readings are
 * confirmed too (they report Jump until then). Once no reading has been accepted (or
 * none ever was) for EncoderGateConfig::stale_time, @ref stale is raised (the fault) and
 * the gate unlocks; it clears at the next Resync.
 *
 * Gate before anything stateful: a PhaseUnwrapper fed a glitch would count it as a turn.
 * For a single-turn sensor, place the reading on the turn nearest the prediction first.
 *
 * @code
 * EncoderGate<float> gate{{.tolerance = 4.0F, .max_accel = 2e9F, .stale_time = 1e-3F}};
 * // each encoder frame:
 * gate.update(tracker, SplitReal<float>{count, 0.0F}, dt, crc_ok);
 * if (gate.stale()) { fault(); }
 * @endcode
 *
 * @tparam T Scalar type.
 */
template<typename T = float>
class EncoderGate {
public:
    constexpr EncoderGate() noexcept = default;
    constexpr explicit EncoderGate(const EncoderGateConfig<T>& config) noexcept : cfg_(config) {}

    /**
     * @brief Judge a reading.
     * @param reading [counts] absolute position
     * @param speed   [counts/s] speed estimate the prediction uses (the tracker's ω)
     * @param dt      [s] time since the previous call
     * @param valid   false when the interface flags the frame as bad (CRC, timeout)
     */
    [[nodiscard]] constexpr EncoderVerdict check(const SplitReal<T>& reading, T speed, T dt, bool valid = true) {
        const T step = (damp::isfinite(dt) && dt > T{0}) ? dt : T{0};
        const T w = damp::isfinite(speed) ? speed : T{0};
        since_ = saturated(since_ + step);
        cand_since_ = saturated(cand_since_ + step);
        if (cfg_.stale_time > T{0} && since_ > cfg_.stale_time) {
            stale_ = true;
            locked_ = false;
        }
        if (!valid || !damp::isfinite(reading.frac)) {
            return EncoderVerdict::Invalid;
        }
        if (locked_ && within(reading - (last_ + (w * since_)), since_)) {
            interval_ = since_;
            last_ = reading;
            since_ = T{0};
            candidates_ = 0;
            return EncoderVerdict::Accepted;
        }
        return candidate(reading, w);
    }

    /**
     * @brief Gate a reading into @p tracker: accepted readings update it, a Resync restarts it at
     *        the reading with its speed, dropped readings never reach it.
     * @return the verdict, as @ref check
     */
    constexpr EncoderVerdict update(EncoderTracker<T>& tracker, const SplitReal<T>& reading, T dt, bool valid = true) {
        const EncoderVerdict v = check(reading, tracker.omega(), dt, valid);
        if (v == EncoderVerdict::Accepted) {
            tracker.update(reading, interval_);
        } else if (v == EncoderVerdict::Resync) {
            tracker.reset(reading, tracker.omega());
        }
        return v;
    }

    [[nodiscard]] constexpr bool                locked() const { return locked_; }     ///< a position is confirmed and readings are being accepted
    [[nodiscard]] constexpr bool                stale() const { return stale_; }       ///< fault: no reading accepted for stale_time (clears at the next Resync)
    [[nodiscard]] constexpr T                   interval() const { return interval_; } ///< [s] between the last two accepted readings (the tracker's dt)
    [[nodiscard]] constexpr const SplitReal<T>& last() const { return last_; }         ///< [counts] last accepted reading

    /** @brief Forget everything: unlocked, not stale; the next readings must confirm a position. */
    constexpr void reset() {
        last_ = {};
        cand_ = {};
        since_ = T{0};
        cand_since_ = T{0};
        interval_ = T{0};
        candidates_ = 0;
        locked_ = false;
        stale_ = false;
    }

private:
    static constexpr T max_time = static_cast<T>(1e6); ///< [s] time accumulators saturate here (finite window)

    [[nodiscard]] static constexpr T saturated(T t) { return (t < max_time) ? t : max_time; }

    /** @brief Inside the window after @p t seconds? Off when tolerance ≤ 0. */
    [[nodiscard]] constexpr bool within(T deviation, T t) const {
        if (!(cfg_.tolerance > T{0})) {
            return true;
        }
        const T accel = (cfg_.max_accel > T{0}) ? cfg_.max_accel : T{0};
        return damp::abs(deviation) <= cfg_.tolerance + (static_cast<T>(0.5) * accel * t * t);
    }

    /** @brief Extend (or start) the candidate track with a reading the locked window did not take. */
    constexpr EncoderVerdict candidate(const SplitReal<T>& reading, T w) {
        const bool extends = candidates_ > 0 && within(reading - (cand_ + (w * cand_since_)), cand_since_);
        candidates_ = extends ? static_cast<uint16_t>(candidates_ + 1U) : uint16_t{1};
        cand_ = reading;
        cand_since_ = T{0};
        if (candidates_ < cfg_.resync_count) {
            return EncoderVerdict::Jump;
        }
        last_ = reading;
        interval_ = T{0};
        since_ = T{0};
        candidates_ = 0;
        locked_ = true;
        stale_ = false;
        return EncoderVerdict::Resync;
    }

    EncoderGateConfig<T> cfg_{};
    SplitReal<T>         last_{};           ///< [counts] last accepted reading
    SplitReal<T>         cand_{};           ///< [counts] latest reading of the candidate track
    T                    since_{T{0}};      ///< [s] since the last accepted reading
    T                    cand_since_{T{0}}; ///< [s] since the latest candidate reading
    T                    interval_{T{0}};   ///< [s] between the last two accepted readings
    uint16_t             candidates_{0};    ///< length of the candidate track
    bool                 locked_{false};    ///< a position is confirmed
    bool                 stale_{false};     ///< no reading accepted for stale_time
};

} // namespace damp
