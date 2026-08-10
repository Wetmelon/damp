// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file toolbox/logic.hpp
 * @brief Modern C++ discrete logic & timing blocks — edge detect, latch, timers,
 *        debounce, toggle, counter.
 *
 * Value-returning equivalents of the classic PLC function blocks: `operator()`
 * returns the primary output directly (no `.Q`/`.CV` member poking), state is
 * private with named accessors, everything is `constexpr` and allocation-free.
 *
 * Timers take the elapsed time per call (`operator()(in, dt)`), matching the
 * rest of Damp (`PIDController` etc.) — jitter-tolerant and multi-rate. Template
 * the timer on the time type: `OnDelayTimer<float>` counts seconds; an integral
 * `OnDelayTimer<int>` with `dt = 1` counts ticks exactly, same code.
 *
 * @see toolbox/iec61131.hpp (namespace `damp::plc`) for the strict, spec-faithful
 *      IEC 61131-3 surface (`.Q`-member API) these are the modern counterpart to.
 * @see toolbox/conditioning.hpp for the analog signal-shaping siblings.
 */

namespace damp {

/// Rising-edge detector: true on the tick @p x goes false → true.
class RisingEdge {
public:
    constexpr bool operator()(bool x) {
        const bool edge = x && !prev_;
        prev_ = x;
        return edge;
    }
    constexpr void reset() { prev_ = false; }

private:
    bool prev_{false};
};

/// Falling-edge detector: true on the tick @p x goes true → false.
class FallingEdge {
public:
    constexpr bool operator()(bool x) {
        const bool edge = !x && prev_;
        prev_ = x;
        return edge;
    }
    constexpr void reset() { prev_ = false; }

private:
    bool prev_{false};
};

/**
 * @brief Set/reset latch. @tparam SetDominant which input wins when both are
 *        asserted (default: set-dominant, e.g. a trip overriding a clear).
 */
template<bool SetDominant = true>
class Latch {
public:
    constexpr bool operator()(bool set, bool reset) {
        if constexpr (SetDominant) {
            if (set) {
                q_ = true;
            } else if (reset) {
                q_ = false;
            }
        } else {
            if (reset) {
                q_ = false;
            } else if (set) {
                q_ = true;
            }
        }
        return q_;
    }
    [[nodiscard]] constexpr bool value() const { return q_; }
    constexpr void               reset() { q_ = false; }

private:
    bool q_{false};
};

/// Reset-dominant latch (clear wins over set).
using ResetDominantLatch = Latch<false>;

/**
 * @brief On-delay timer: output goes true once @p in has been held true
 *        continuously for @p delay; drops immediately when @p in goes false.
 * @tparam T time type — float for seconds, integral with `dt = 1` for ticks.
 */
template<typename T = float>
class OnDelayTimer {
public:
    constexpr OnDelayTimer() = default;
    constexpr explicit OnDelayTimer(T delay) : delay_(delay) {}

    constexpr bool operator()(bool in, T dt) {
        if (!in) {
            elapsed_ = T{0};
            out_ = false;
        } else if (!out_) {
            elapsed_ += dt;
            if (elapsed_ >= delay_) {
                out_ = true; // latch; stop accumulating so integral T can't overflow
            }
        }
        return out_;
    }

    [[nodiscard]] constexpr bool value() const { return out_; }
    [[nodiscard]] constexpr T    elapsed() const { return elapsed_; }
    constexpr void               reset() {
        elapsed_ = T{0};
        out_ = false;
    }

private:
    T    delay_{};
    T    elapsed_{};
    bool out_{false};
};

/**
 * @brief Off-delay timer: output goes true immediately when @p in is true and
 *        stays true until @p in has been false continuously for @p delay.
 * @tparam T time type (see @ref OnDelayTimer).
 */
template<typename T = float>
class OffDelayTimer {
public:
    constexpr OffDelayTimer() = default;
    constexpr explicit OffDelayTimer(T delay) : delay_(delay) {}

    constexpr bool operator()(bool in, T dt) {
        if (in) {
            out_ = true;
            elapsed_ = T{0};
        } else if (out_) {
            elapsed_ += dt;
            if (elapsed_ >= delay_) {
                out_ = false;
            }
        }
        return out_;
    }

    [[nodiscard]] constexpr bool value() const { return out_; }
    [[nodiscard]] constexpr T    elapsed() const { return elapsed_; }
    constexpr void               reset() {
        elapsed_ = T{0};
        out_ = false;
    }

private:
    T    delay_{};
    T    elapsed_{};
    bool out_{false};
};

/**
 * @brief Pulse timer (non-retriggerable): a rising edge of @p in emits a fixed
 *        @p width output pulse; edges during the pulse are ignored.
 * @tparam T time type (see @ref OnDelayTimer).
 */
template<typename T = float>
class PulseTimer {
public:
    constexpr PulseTimer() = default;
    constexpr explicit PulseTimer(T width) : width_(width) {}

    constexpr bool operator()(bool in, T dt) {
        const bool rising = in && !prev_;
        prev_ = in;
        if (rising && !out_) {
            out_ = true;
            elapsed_ = T{0};
        }
        if (out_) {
            elapsed_ += dt;
            if (elapsed_ >= width_) {
                out_ = false;
            }
        }
        return out_;
    }

    [[nodiscard]] constexpr bool value() const { return out_; }
    [[nodiscard]] constexpr T    elapsed() const { return elapsed_; }
    constexpr void               reset() {
        elapsed_ = T{0};
        out_ = false;
        prev_ = false;
    }

private:
    T    width_{};
    T    elapsed_{};
    bool out_{false};
    bool prev_{false};
};

/**
 * @brief Debounce: the output adopts @p in only after @p in differs from the
 *        current output continuously for @p stable_time. Rejects contact bounce
 *        and brief glitches. (Not an IEC block — the one everyone hand-rolls.)
 * @tparam T time type (see @ref OnDelayTimer).
 */
template<typename T = float>
class Debounce {
public:
    constexpr Debounce() = default;
    constexpr explicit Debounce(T stable_time, bool initial = false)
        : stable_(stable_time), state_(initial) {}

    constexpr bool operator()(bool in, T dt) {
        if (in == state_) {
            elapsed_ = T{0}; // no pending change
        } else {
            elapsed_ += dt;
            if (elapsed_ >= stable_) {
                state_ = in;
                elapsed_ = T{0};
            }
        }
        return state_;
    }

    [[nodiscard]] constexpr bool value() const { return state_; }
    constexpr void               reset(bool state = false) {
        state_ = state;
        elapsed_ = T{0};
    }

private:
    T    stable_{};
    T    elapsed_{};
    bool state_{false};
};

/// Toggle (T flip-flop): output flips on each rising edge of @p in.
class Toggle {
public:
    constexpr explicit Toggle(bool initial = false) : q_(initial) {}

    constexpr bool operator()(bool in) {
        if (in && !prev_) {
            q_ = !q_;
        }
        prev_ = in;
        return q_;
    }
    [[nodiscard]] constexpr bool value() const { return q_; }
    constexpr void               reset(bool state = false) {
        q_ = state;
        prev_ = false;
    }

private:
    bool q_{false};
    bool prev_{false};
};

/**
 * @brief Edge-counting up/down counter: increments on each rising edge of @p up,
 *        decrements on each rising edge of @p down. Returns the running count.
 * @tparam T integer count type.
 */
template<typename T = int>
class Counter {
public:
    constexpr Counter() = default;
    constexpr explicit Counter(T initial) : count_(initial) {}

    constexpr T operator()(bool up, bool down) {
        if (up && !prev_up_) {
            ++count_;
        }
        if (down && !prev_down_) {
            --count_;
        }
        prev_up_ = up;
        prev_down_ = down;
        return count_;
    }

    [[nodiscard]] constexpr T value() const { return count_; }
    constexpr void            reset(T value = T{0}) {
        count_ = value;
        prev_up_ = false;
        prev_down_ = false;
    }

private:
    T    count_{};
    bool prev_up_{false};
    bool prev_down_{false};
};

} // namespace damp
