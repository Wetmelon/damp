// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file toolbox/io.hpp
 * @brief Tier-2 I/O appliances — batteries-included field-I/O channels.
 *
 * Where toolbox/conditioning.hpp, toolbox/logic.hpp, and toolbox/scaling.hpp
 * provide the primitives (range monitor, debounce, edge detect, affine cal),
 * this header composes them into the constructs you actually wire up: a
 * calibrated+validated analog input, a debounced push-button with edges and
 * hold, a debounced maintained switch. Live in `damp::io` to mark them as the
 * higher-level composites.
 *
 * @see toolbox/conditioning.hpp, toolbox/logic.hpp, toolbox/scaling.hpp (tier 1).
 */

#include "damp/toolbox/conditioning.hpp" // RangeMonitor, SignalStatus, is_valid/is_fault
#include "damp/toolbox/logic.hpp"        // Debounce, OnDelayTimer
#include "damp/toolbox/scaling.hpp"      // AffineCal

namespace damp::io {

/**
 * @brief A single analog input: range/fault check on the raw reading, then
 *        affine calibration to engineering units.
 *
 * The everyday "read one sensor" appliance. The raw value is classified first
 * (a wire fault is a property of the raw V/mA, not the scaled value), then
 * scaled. `update()` returns the engineering value; `status()`/`faulted()`
 * report the channel health from the same tick.
 *
 * @code
 * io::AnalogInput<float> level{
 *     two_point_cal(0.5f, 0.0f, 4.5f, 100.0f),   // 0.5–4.5 V → 0–100 %
 *     RangeMonitor<float>{0.25f, 0.5f, 4.5f, 4.75f}}; // NE43 fault bands
 * float pct = level.update(volts);
 * if (level.faulted()) flag_wire_break();
 * @endcode
 */
template<typename T = float>
class AnalogInput {
public:
    constexpr AnalogInput() = default;
    constexpr AnalogInput(const AffineCal<T>& cal, const RangeMonitor<T>& monitor)
        : cal_(cal), monitor_(monitor) {}

    /// Validate the raw reading, scale it, and store both results.
    constexpr T update(T raw) {
        status_ = monitor_(raw);
        value_ = cal_.apply(raw);
        return value_;
    }

    [[nodiscard]] constexpr T            value() const { return value_; }
    [[nodiscard]] constexpr SignalStatus status() const { return status_; }
    [[nodiscard]] constexpr bool         valid() const { return is_valid(status_); }
    [[nodiscard]] constexpr bool         faulted() const { return is_fault(status_); }

private:
    AffineCal<T>    cal_{};
    RangeMonitor<T> monitor_{};
    T               value_{};
    SignalStatus    status_{SignalStatus::Valid};
};

/**
 * @brief Dual-channel analog input with same-slope or opposite-slope cross-check.
 *
 * Each channel is an AnalogInput (range/fault on raw, then affine to engineering
 * units). Each tick both channels update; @ref cross_check_analog votes them with
 * @ref AnalogCrossMode. Hydraulic and safety I/O use this for redundant pressure /
 * position / joystick sensors.
 *
 * @code
 * io::RedundantAnalogPair<float> rail{
 *     {two_point_cal(0.5f, 0.f, 4.5f, 250.f), RangeMonitor<float>{0.25f, 0.5f, 4.5f, 4.75f}},
 *     {two_point_cal(0.5f, 0.f, 4.5f, 250.f), RangeMonitor<float>{0.25f, 0.5f, 4.5f, 4.75f}},
 *     AnalogCrossMode::SameSlope,
 *     5.f}; // 5 bar agreement
 * auto r = rail.update(v_a, v_b);
 * if (r.agree()) set_pressure(r.value);
 * @endcode
 */
template<typename T = float>
class RedundantAnalogPair {
public:
    constexpr RedundantAnalogPair() = default;

    constexpr RedundantAnalogPair(
        AnalogInput<T>  channel_a,
        AnalogInput<T>  channel_b,
        AnalogCrossMode mode,
        T               tol,
        T               span = T{0}
    )
        : a_(channel_a), b_(channel_b), mode_(mode), tol_(tol), span_(span) {}

    /// Update both raw readings; return the cross-check result.
    constexpr AnalogCrossResult<T> update(T raw_a, T raw_b) {
        a_.update(raw_a);
        b_.update(raw_b);
        last_ = cross_check_analog(a_.value(), a_.status(), b_.value(), b_.status(), mode_, tol_, span_);
        return last_;
    }

    [[nodiscard]] constexpr const AnalogCrossResult<T>& last() const { return last_; }
    [[nodiscard]] constexpr const AnalogInput<T>&       channel_a() const { return a_; }
    [[nodiscard]] constexpr const AnalogInput<T>&       channel_b() const { return b_; }

private:
    AnalogInput<T>       a_{};
    AnalogInput<T>       b_{};
    AnalogCrossMode      mode_{AnalogCrossMode::SameSlope};
    T                    tol_{};
    T                    span_{};
    AnalogCrossResult<T> last_{};
};

/**
 * @brief Affine map from a normalized command to a proportional current (mA).
 *
 * Hydraulic valves and solenoid drivers are typically commanded as a percent (or
 * @f$[0,1]@f$) that maps linearly onto a coil current window
 * @f$ [I_{\min}, I_{\max}] @f$ (often with a deadband or pilot current at zero).
 * This is that map only — PWM current regulation lives in the driver hardware or
 * a PI current loop on the same stack.
 *
 * @code
 * io::ProportionalCurrentMap<float> valve{200.f, 800.f}; // pilot 200 mA, full 800 mA
 * float mA = valve(0.5f); // mid command → 500 mA
 * @endcode
 */
template<typename T = float>
class ProportionalCurrentMap {
public:
    constexpr ProportionalCurrentMap() = default;

    /// @param i_min Current at command 0 [mA]. @param i_max Current at command 1 [mA].
    constexpr ProportionalCurrentMap(T i_min, T i_max) : i_min_(i_min), i_max_(i_max) {}

    /// Map command in @f$[0,1]@f$ (clamped) to mA.
    [[nodiscard]] constexpr T operator()(T command) const {
        const T u = damp::clamp(command, T{0}, T{1});
        return i_min_ + ((i_max_ - i_min_) * u);
    }

    /// Inverse: mA → command in @f$[0,1]@f$ (clamped). Requires @p i_max ≠ @p i_min.
    [[nodiscard]] constexpr T invert(T milliamps) const {
        if (i_max_ == i_min_) {
            return T{0};
        }
        return damp::clamp((milliamps - i_min_) / (i_max_ - i_min_), T{0}, T{1});
    }

    [[nodiscard]] constexpr T i_min() const { return i_min_; }
    [[nodiscard]] constexpr T i_max() const { return i_max_; }

private:
    T i_min_{T{0}};
    T i_max_{T{1}};
};

/**
 * @brief Operator-axis conditioning chain (joystick / RC stick → command).
 *
 * The teleop counterpart to AnalogInput: instead of validating a sensor, it
 * *shapes an operator input*. Each `update(raw)` runs the standard chain —
 * affine cal to normalized `[-1, 1]` → center scaled dead zone → exponential
 * response curve → output scale (with optional inversion):
 *
 *   1. `cal.apply(raw)`        raw counts/volts → ~`[-1, 1]` (clamped before shaping)
 *   2. `scaled_deadband(·, dz)` ignore slop around center, full range preserved
 *   3. `expo(·, k)`            soften near center for fine control
 *   4. `· × scale` (± invert)  to the actuator command range
 *
 * @code
 * io::AxisInput<float> steer{
 *     two_point_cal(0.0f, -1.0f, 4095.0f, 1.0f), // 12-bit ADC → [-1, 1]
 *     0.05f,   // 5% center dead zone
 *     0.6f,    // expo
 *     1.0f};   // output scale
 * float cmd = steer.update(adc);
 * @endcode
 */
template<typename T = float>
class AxisInput {
public:
    constexpr AxisInput() = default;
    constexpr AxisInput(const AffineCal<T>& cal, T deadzone = T{0}, T expo_k = T{0}, T scale = T{1}, bool invert = false)
        : cal_(cal), deadzone_(deadzone), expo_(expo_k), scale_(scale), invert_(invert) {}

    constexpr T update(T raw) {
        T x = damp::clamp(cal_.apply(raw), T{-1}, T{1}); // normalize, guard cal overshoot
        x = scaled_deadband(x, deadzone_);
        x = expo(x, expo_);
        if (invert_) {
            x = -x;
        }
        value_ = scale_ * x;
        return value_;
    }

    [[nodiscard]] constexpr T value() const { return value_; }

    constexpr void               set_deadzone(T dz) { deadzone_ = dz; }
    constexpr void               set_expo(T k) { expo_ = k; }
    constexpr void               set_scale(T s) { scale_ = s; }
    constexpr void               set_invert(bool inv) { invert_ = inv; }
    [[nodiscard]] constexpr T    deadzone() const { return deadzone_; }
    [[nodiscard]] constexpr T    expo_k() const { return expo_; } ///< expo factor (not damp::expo)
    [[nodiscard]] constexpr T    scale() const { return scale_; }
    [[nodiscard]] constexpr bool invert() const { return invert_; }

private:
    AffineCal<T> cal_{};
    T            deadzone_{T{0}};
    T            expo_{T{0}};
    T            scale_{T{1}};
    bool         invert_{false};
    T            value_{T{0}};
};

/**
 * @brief Debounced momentary push-button with edge and long-press detection.
 *
 * `down()` is the debounced level; `pressed()`/`released()` are one-tick edges;
 * `held()` is true once the button has been down continuously past the hold
 * time. A zero hold time (default) makes `held()` track `down()`.
 *
 * @tparam T time type (float seconds, or integral with `dt = 1` for ticks).
 */
template<typename T = float>
class Button {
public:
    constexpr Button() = default;
    constexpr explicit Button(T debounce_time, T hold_time = T{0})
        : db_(debounce_time), hold_(hold_time) {}

    /// Advance one tick with the raw electrical state and elapsed time.
    constexpr void update(bool raw, T dt) {
        const bool s = db_(raw, dt);
        pressed_ = s && !prev_;
        released_ = !s && prev_;
        prev_ = s;
        held_ = hold_(s, dt);
    }

    [[nodiscard]] constexpr bool down() const { return prev_; }         ///< debounced level
    [[nodiscard]] constexpr bool pressed() const { return pressed_; }   ///< rising edge (one tick)
    [[nodiscard]] constexpr bool released() const { return released_; } ///< falling edge (one tick)
    [[nodiscard]] constexpr bool held() const { return held_; }         ///< down past the hold time

    constexpr void reset() {
        db_.reset();
        hold_.reset();
        prev_ = pressed_ = released_ = held_ = false;
    }

private:
    Debounce<T>     db_{};
    OnDelayTimer<T> hold_{};
    bool            prev_{false};
    bool            pressed_{false};
    bool            released_{false};
    bool            held_{false};
};

/**
 * @brief Debounced maintained switch (toggle/selector contact) with change flag.
 *
 * `on()` is the debounced level; `changed()` is true the tick the debounced
 * state flips. The smoothed contact for a panel switch or limit/proximity sensor.
 *
 * @tparam T time type (see @ref Button).
 */
template<typename T = float>
class Switch {
public:
    constexpr Switch() = default;
    constexpr explicit Switch(T debounce_time, bool initial = false)
        : db_(debounce_time, initial), prev_(initial) {}

    /// Advance one tick; returns the debounced state.
    constexpr bool update(bool raw, T dt) {
        const bool s = db_(raw, dt);
        changed_ = (s != prev_);
        prev_ = s;
        return s;
    }

    [[nodiscard]] constexpr bool on() const { return prev_; }
    [[nodiscard]] constexpr bool changed() const { return changed_; }

    constexpr void reset(bool initial = false) {
        db_.reset(initial);
        prev_ = initial;
        changed_ = false;
    }

private:
    Debounce<T> db_{};
    bool        prev_{false};
    bool        changed_{false};
};

/**
 * @brief Left/right actuator pair (differential/arcade drive output).
 */
template<typename T = float>
struct DrivePair {
    T left{};
    T right{};
};

/**
 * @brief Differential (tank) drive mixer: throttle + turn → left/right.
 *
 * The tank-drive kinematics
 * @f[ \text{left} = v + \omega,\qquad \text{right} = v - \omega @f]
 * with anti-clip normalization: if either wheel command would exceed
 * @f$[-1, 1]@f$, both are scaled down by the overshoot so the turn/throttle
 * *ratio* is preserved instead of hard-clamping one wheel (which distorts the
 * commanded curvature). Inputs are assumed already in @f$[-1, 1]@f$ (e.g. from
 * @ref AxisInput). Dual-path hydrostatic drives are the headline user of the
 * ratio-preserving clip.
 *
 * @param throttle forward/back command in [-1, 1]
 * @param turn     turn-rate command in [-1, 1] (positive → turn one way)
 *
 * @code
 * auto [l, r] = io::differential_drive(0.8f, 0.5f); // preserves 0.8:0.5 split
 * @endcode
 */
template<typename T>
[[nodiscard]] constexpr DrivePair<T> differential_drive(T throttle, T turn) {
    T       left = throttle + turn;
    T       right = throttle - turn;
    const T peak = damp::max(damp::abs(left), damp::abs(right));
    if (peak > T{1}) {
        left /= peak;
        right /= peak;
    }
    return {left, right};
}

/**
 * @brief Arcade (single-stick) drive mixer: x (steer) + y (throttle) → left/right.
 *
 * Identical kinematics to @ref differential_drive with single-stick argument
 * ordering (horizontal steer axis, vertical throttle axis).
 *
 * @param x steering axis in [-1, 1]
 * @param y throttle axis in [-1, 1]
 */
template<typename T>
[[nodiscard]] constexpr DrivePair<T> arcade_drive(T x, T y) {
    return differential_drive(y, x);
}

/**
 * @brief H-pattern (dual-lever / tank) drive: one lever per track, no mixing.
 *
 * The classic dozer/skid-steer lever control — each stick drives its own track
 * directly (both forward → straight, split → spin-in-place). No cross-mixing,
 * so nothing can exceed range from combination; each lever is only clamped to
 * @f$[-1, 1]@f$ to guard input overshoot.
 *
 * @param left_lever  left-track command in [-1, 1]
 * @param right_lever right-track command in [-1, 1]
 */
template<typename T>
[[nodiscard]] constexpr DrivePair<T> h_pattern(T left_lever, T right_lever) {
    return {damp::clamp(left_lever, T{-1}, T{1}), damp::clamp(right_lever, T{-1}, T{1})};
}

/**
 * @brief ISO-S (single-stick, speed-summed) drive — arcade travel stick.
 *
 * The ISO travel stick with speed-summed steering: @f$\text{left}=v+\omega@f$,
 * @f$\text{right}=v-\omega@f$ with ratio-preserving anti-clip. Identical to
 * @ref differential_drive; the turn is a fixed track *difference*, so a turn
 * command at zero throttle spins in place. Contrast @ref iso_c_drive.
 *
 * @param throttle forward/back command in [-1, 1]
 * @param turn     steer command in [-1, 1]
 */
template<typename T>
[[nodiscard]] constexpr DrivePair<T> iso_s_drive(T throttle, T turn) {
    return differential_drive(throttle, turn);
}

/**
 * @brief ISO-C (single-stick, coordinated / curvature) drive.
 *
 * The coordinated ISO travel stick: the steer axis commands a *curvature*, not
 * a fixed track difference — the turn effect scales with throttle
 * @f[ \text{left} = v\,(1+\omega),\qquad \text{right} = v\,(1-\omega), @f]
 * so @f$\omega@f$ sets the turn radius (full steer pivots about the inner
 * track) and there is no spin-in-place at zero throttle. Ratio-preserving
 * anti-clip as in @ref differential_drive. Smoother for travel/road use where
 * ISO-S's zero-speed spin is undesirable.
 *
 * @param throttle forward/back command in [-1, 1]
 * @param turn     curvature command in [-1, 1]
 */
template<typename T>
[[nodiscard]] constexpr DrivePair<T> iso_c_drive(T throttle, T turn) {
    T       left = throttle * (T{1} + turn);
    T       right = throttle * (T{1} - turn);
    const T peak = damp::max(damp::abs(left), damp::abs(right));
    if (peak > T{1}) {
        left /= peak;
        right /= peak;
    }
    return {left, right};
}

/**
 * @brief Curve-driven steering map — arbitrary `(throttle, turn)` geometry.
 *
 * Where @ref iso_s_drive / @ref iso_c_drive fix the surface in closed form, this
 * reads the inner-track command from a caller-supplied 2-D curve over
 * @f$(|throttle|, |turn|)@f$: the outer track runs at the commanded throttle,
 * the inner track takes the curve value, the sign of `turn` picks which side is
 * inner, and the sign of `throttle` sets travel direction (car-like in reverse,
 * as in @ref iso_c_drive). This expresses geometries neither closed form can —
 * e.g. pivot about the inner track at low speed yet hold a fixed inner-track
 * ratio at high speed (a 3×3 `Lut2D` calibration).
 *
 * @p inner_curve is *any* callable @f$(|throttle|, |turn|) \to@f$ inner-track
 * command, so the interpolation is the caller's choice — a bilinear @ref Lut2D,
 * a spline/expo-shaped surface, or a closed-form lambda. Outputs are clamped to
 * @f$[-1, 1]@f$ to guard a mis-authored table.
 *
 * @param inner_curve inner-track response surface, indexed by (|throttle|, |turn|)
 * @param throttle    forward/back command in [-1, 1]
 * @param turn        steer command in [-1, 1]
 *
 * @code
 * // Pivot at half throttle, inner track holds 0.5 at full throttle.
 * Lut2D<3, 3, float> inner{{0,0.5f,1}, {0,0.5f,1},
 *     {{0,0,0}, {0.5f,0.25f,0}, {1,0.75f,0.5f}}};
 * auto [l, r] = io::steer_map(inner, throttle, turn);
 * @endcode
 */
template<typename Surface, typename T>
[[nodiscard]] constexpr DrivePair<T> steer_map(const Surface& inner_curve, T throttle, T turn) {
    const T sign = (throttle < T{0}) ? T{-1} : T{1};
    const T outer = damp::clamp(throttle, T{-1}, T{1});
    const T inner = damp::clamp(sign * inner_curve(damp::abs(throttle), damp::abs(turn)), T{-1}, T{1});
    return (turn < T{0}) ? DrivePair<T>{inner, outer} : DrivePair<T>{outer, inner};
}

/**
 * @brief Four-wheel actuator set (mecanum/holonomic drive output), X-config.
 */
template<typename T = float>
struct HolonomicOutput {
    T front_left{};
    T front_right{};
    T rear_left{};
    T rear_right{};
};

/**
 * @brief Four-wheel holonomic drive mixer (mecanum / omni).
 *
 * Maps a body-frame velocity command to the four wheel speeds:
 * @f[
 *   \begin{aligned}
 *   \text{FL} &= y + x + \omega, & \text{FR} &= y - x - \omega, \\
 *   \text{RL} &= y - x + \omega, & \text{RR} &= y + x - \omega,
 *   \end{aligned}
 * @f]
 * where @f$y@f$ is forward, @f$x@f$ is strafe-right, and @f$\omega@f$ is yaw. The
 * @f$\pm x@f$ handedness follows the diagonal (X) roller pattern (`FL, RR` one lay;
 * `FR, RL` the other) — the only pattern that strafes; a left/right split would
 * spin instead. Mecanum-H and omni-X reduce to this same law, so they are aliases
 * (@ref mecanum_drive, @ref omniwheel_drive) rather than separate kinematics.
 *
 * Anti-clip normalization scales all four wheels by the common overshoot so the
 * translation/rotation blend is preserved when a combined command saturates. Pure
 * forward/strafe/yaw inputs recover the canonical holonomic wheel patterns.
 *
 * @param x        strafe-right command in [-1, 1]
 * @param y        forward command in [-1, 1]
 * @param rotation yaw command in [-1, 1]
 */
template<typename T>
[[nodiscard]] constexpr HolonomicOutput<T> holonomic_drive(T x, T y, T rotation) {
    T fl = y + x + rotation;
    T fr = y - x - rotation;
    T rl = y - x + rotation;
    T rr = y + x - rotation;

    const T peak = damp::max({damp::abs(fl), damp::abs(fr), damp::abs(rl), damp::abs(rr)});
    if (peak > T{1}) {
        fl /= peak;
        fr /= peak;
        rl /= peak;
        rr /= peak;
    }

    return {fl, fr, rl, rr};
}

/// Mecanum (H-layout, 45° rollers) mixer — alias for @ref holonomic_drive.
template<typename T>
[[nodiscard]] constexpr HolonomicOutput<T> mecanum_drive(T x, T y, T rotation) {
    return holonomic_drive(x, y, rotation);
}

/// Omni-wheel (X-layout) mixer — alias for @ref holonomic_drive.
template<typename T>
[[nodiscard]] constexpr HolonomicOutput<T> omniwheel_drive(T x, T y, T rotation) {
    return holonomic_drive(x, y, rotation);
}

} // namespace damp::io
