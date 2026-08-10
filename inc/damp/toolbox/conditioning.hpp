// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file toolbox/conditioning.hpp
 * @brief Nonlinear signal-shaping primitives for real control loops.
 *
 * The everyday command/feedback conditioning blocks that aren't filters and
 * aren't calibration math: dead zones, their inverse (hard step, soft
 * sine / tanh-like, Coulomb–viscous / Stribeck friction maps), rate limiting,
 * and a Schmitt-trigger hysteresis comparator. All `constexpr` and
 * allocation-free — drop them straight into an ISR.
 *
 * Forward deadband/expo live on the stick path (reject noise). inverse
 * deadband / friction live on the actuator command after the online motion
 * planner (valve overlap, stiction, Coulomb) — not on the stick.
 *
 * @see toolbox/scaling.hpp for stateless range/affine calibration helpers.
 * @see filters/filters.hpp for the linear filter runtimes.
 */

#include <cstdint>
#include <limits>

#include "damp/math/math.hpp" // damp::abs, damp::copysign, damp::sin, damp::sqrt, damp::exp

namespace damp {

/**
 * @brief Dead zone over `[lower, upper]`, matching Simulink®'s Dead Zone block.
 *
 * Output is 0 inside the band and the signed distance past the edge outside it
 * (`x − lower` below, `x − upper` above) — continuous through both edges. The
 * standard joystick / valve-command dead zone that suppresses noise around
 * neutral without introducing a step.
 *
 * @param x     Input.
 * @param lower Lower edge of the dead zone.
 * @param upper Upper edge of the dead zone (`upper ≥ lower`).
 * @return 0 within `[lower, upper]`, else `x − lower` / `x − upper`.
 */
template<typename T>
[[nodiscard]] constexpr T deadband(T x, T lower, T upper) {
    if (x > upper) {
        return x - upper;
    }
    if (x < lower) {
        return x - lower;
    }
    return T{0};
}

/// Symmetric dead zone of half-width @p band (i.e. `[−band, band]`).
template<typename T>
[[nodiscard]] constexpr T deadband(T x, T band) {
    return deadband(x, -band, band);
}

/**
 * @brief Inverse dead zone: add an offset to overcome a physical dead zone
 *        (valve overlap, static friction, motor stiction), with independent
 *        negative/positive offsets.
 *
 * A positive command is boosted by @p upper, a negative one by @p lower (use a
 * negative @p lower for a symmetric-feeling boost), so the actuator actually
 * starts moving. Commands within ±@p threshold are forced to 0 so sensor noise
 * doesn't cause jitter/dither. This is the compensating inverse of @ref deadband.
 *
 * @param x         Input command.
 * @param lower     Offset added to negative commands (typically ≤ 0).
 * @param upper     Offset added to positive commands (typically ≥ 0).
 * @param threshold Magnitude below which the output is forced to 0 (default 0).
 */
template<typename T>
[[nodiscard]] constexpr T inverse_deadband(T x, T lower, T upper, T threshold = T{0}) {
    if (damp::abs(x) <= threshold) {
        return T{0};
    }
    return x + (x > T{0} ? upper : lower);
}

/// Symmetric inverse dead zone: boost by @p band in the command's direction.
template<typename T>
[[nodiscard]] constexpr T inverse_deadband(T x, T band) {
    return inverse_deadband(x, -band, band, T{0});
}

/**
 * @brief Smooth unit direction @f$ x / \sqrt{x^{2}+\varepsilon^{2}} \in (-1,1) @f$
 *
 * Algebraic soft-sign (C∞, no transcendental overflow). As @f$\varepsilon \to 0@f$
 * this approaches @f$\mathrm{sgn}(x)@f$; larger @p eps softens the transition.
 * Used for smooth Coulomb / stiction maps (tanh-like shape).
 *
 * @param x   Input (velocity or command).
 * @param eps Softness scale (same units as @p x); use a small positive value.
 */
template<typename T>
[[nodiscard]] constexpr T soft_sign(T x, T eps) {
    // Floor eps so we never divide by zero under finite-math builds
    const T e = damp::max(damp::abs(eps), std::numeric_limits<T>::epsilon());
    return x / damp::sqrt((x * x) + (e * e));
}

/**
 * @brief Smooth inverse dead zone via soft-sign (tanh-like Coulomb boost)
 *
 * @f[
 *   y = x + b\,\mathrm{soft\_sign}(x,\varepsilon)
 * @f]
 * with optional center hold: @f$|x|\le\mathrm{threshold}\Rightarrow y=0@f$.
 * As @p eps → 0 this matches @ref inverse_deadband (hard ±band). Larger @p eps
 * eases the boost in near zero (less dither, softer breakaway).
 *
 * Apply on the actuator command after OTG, not on the stick.
 *
 * @param x         Desired command (after planner).
 * @param band      Asymptotic boost magnitude (same units as @p x).
 * @param eps       Softness of the sign (same units); @f$\sim 0.01\cdot\mathrm{scale}@f$.
 * @param threshold Magnitude below which output is forced to 0 (default 0).
 * @see soft_sign, inverse_deadband, inverse_deadband_sine
 */
template<typename T>
[[nodiscard]] constexpr T inverse_deadband_soft(T x, T band, T eps, T threshold = T{0}) {
    if (damp::abs(x) <= threshold) {
        return T{0};
    }
    return x + (damp::abs(band) * soft_sign(x, eps));
}

/**
 * @brief Sine-shaped inverse dead zone (full boost once @f$|x|\ge w@f$)
 *
 * @f[
 *   y = x + b\,\sin\bigl(\tfrac{\pi}{2}\,\mathrm{sat}(x/w)\bigr)
 * @f]
 * Offset ramps smoothly for @f$|x| < w@f$ and equals @f$\pm b@f$ outside.
 * Continuous and odd; derivative is continuous at @f$\pm w@f$ in the sin sense
 * (flat slope of the offset at the ends of the ramp).
 *
 * @param x         Desired command.
 * @param band      Full boost once past @p width.
 * @param width     Half-width of the sine ramp (must be > 0).
 * @param threshold Center hold; output 0 when @f$|x|\le@f$ threshold.
 * @see inverse_deadband_soft, inverse_deadband
 */
template<typename T>
[[nodiscard]] constexpr T inverse_deadband_sine(T x, T band, T width, T threshold = T{0}) {
    if (damp::abs(x) <= threshold) {
        return T{0};
    }
    const T w = damp::max(damp::abs(width), std::numeric_limits<T>::epsilon());
    const T s = damp::clamp(x / w, T{-1}, T{1});
    const T offset = damp::abs(band) * damp::sin(static_cast<T>(0.5) * damp::numbers::pi_v<T> * s);
    return x + offset;
}

/**
 * @brief Coulomb + viscous friction compensation (static inverse, massless)
 *
 * @f[
 *   u = u^{\star} + F_{c}\,\mathrm{soft\_sign}(u^{\star},\varepsilon) + b\,u^{\star}
 * @f]
 * Boosts the command so a plant with Coulomb stiction and linear viscous drag
 * can track @f$u^{\star}@f$ in steady state. Purely algebraic — no internal state.
 *
 * @param cmd      Desired command / velocity proxy.
 * @param coulomb  Coulomb magnitude @f$F_c \ge 0@f$.
 * @param viscous  Viscous coefficient @f$b@f$ (default 0).
 * @param eps      Soft-sign scale (default: machine-ish from @p coulomb).
 * @see compensate_stribeck, soft_sign, inverse_deadband_soft
 */
template<typename T>
[[nodiscard]] constexpr T compensate_coulomb_viscous(
    T cmd, T coulomb, T viscous = T{0}, T eps = T{0}
) {
    const T fc = damp::max(coulomb, T{0});
    const T e = (eps > T{0}) ? eps : (static_cast<T>(1e-3) * (T{1} + fc));
    return cmd + (fc * soft_sign(cmd, e)) + (viscous * cmd);
}

/**
 * @brief Stribeck-style friction compensation (static + Coulomb + viscous)
 *
 * @f[
 *   F(v) = \Bigl(F_{c} + (F_{s}-F_{c})\,e^{-(v/v_{s})^{2}}\Bigr)
 *          \mathrm{soft\_sign}(v,\varepsilon) + b\,v
 * @f]
 * @f$ u = u^{\star} + F(u^{\star}) @f$. Captures higher breakaway (@p static_f)
 * that falls to Coulomb as speed rises. Still massless / memoryless.
 *
 * @param cmd       Desired command / velocity proxy @f$v@f$.
 * @param coulomb   Sliding Coulomb @f$F_c \ge 0@f$.
 * @param static_f  Breakaway / static level @f$F_s \ge F_c@f$.
 * @param v_stribeck Stribeck velocity scale @f$v_s > 0@f$.
 * @param viscous   Viscous @f$b@f$ (default 0).
 * @param eps       Soft-sign scale (default from Coulomb).
 * @see compensate_coulomb_viscous
 */
template<typename T>
[[nodiscard]] constexpr T compensate_stribeck(
    T cmd,
    T coulomb,
    T static_f,
    T v_stribeck,
    T viscous = T{0},
    T eps = T{0}
) {
    const T fc = damp::max(coulomb, T{0});
    const T fs = damp::max(static_f, fc);
    const T vs = damp::max(damp::abs(v_stribeck), std::numeric_limits<T>::epsilon());
    const T e = (eps > T{0}) ? eps : (static_cast<T>(1e-3) * (T{1} + fc));
    const T r = cmd / vs;
    const T stribeck = fc + ((fs - fc) * damp::exp(-(r * r)));
    return cmd + (stribeck * soft_sign(cmd, e)) + (viscous * cmd);
}

/**
 * @brief Center dead zone that rescales the surviving range back to full span.
 *
 * Unlike @ref deadband (which returns the raw distance past the edge), this maps
 * the dead-zone edge to 0 and the input extreme ±1 back to ±1 — the joystick/RC
 * convention, so removing the dead zone doesn't shrink the usable range.
 *
 * @param x        Normalized input, expected in `[-1, 1]`.
 * @param deadzone Dead-zone half-width as a fraction in `[0, 1)`.
 */
template<typename T>
[[nodiscard]] constexpr T scaled_deadband(T x, T deadzone) {
    const T m = damp::abs(x);
    if (m <= deadzone) {
        return T{0};
    }
    return damp::copysign((m - deadzone) / (T{1} - deadzone), x);
}

/**
 * @brief Exponential response curve `y = (1−k)·x + k·x³` (RC "expo").
 *
 * Softens response near center for fine control while keeping full authority at
 * the extremes: `k = 0` is linear, `k = 1` is fully cubic. Endpoints (±1) and
 * sign are preserved and the curve is monotonic for `k ∈ [0, 1]`.
 *
 * @param x Normalized input, expected in `[-1, 1]`.
 * @param k Expo factor in `[0, 1]`.
 */
template<typename T>
[[nodiscard]] constexpr T expo(T x, T k) {
    return ((T{1} - k) * x) + (k * x * x * x);
}

/**
 * @brief Slew-rate limiter: bound how fast the output may follow the target.
 *
 * Caps the per-step change at `rate·dt`, with independent up/down rates so you
 * can, e.g., ramp a command on slowly but allow a fast retreat. Rates are
 * positive magnitudes in units/second; the default (∞) is pass-through. The
 * first sample seeds the state to the target (no start-up ramp). A non-positive
 * `dt` holds the output (no reverse step from a negative interval).
 *
 * @tparam T Scalar type.
 */
template<typename T = float>
class SlewLimiter {
public:
    constexpr SlewLimiter() = default;

    /// Symmetric rate limit [units/s].
    constexpr explicit SlewLimiter(T rate) : rate_up_(rate), rate_down_(rate) {}

    /// Independent rising/falling rate limits [units/s].
    constexpr SlewLimiter(T rate_up, T rate_down) : rate_up_(rate_up), rate_down_(rate_down) {}

    /// Advance one step toward @p target over @p dt seconds; holds when dt ≤ 0.
    constexpr T operator()(T target, T dt) {
        if (!init_) {
            y_ = target;
            init_ = true;
            return y_;
        }
        if (!(dt > T{0})) {
            return y_; // non-positive step: hold output (no reverse / NaN rate·dt)
        }
        const T max_up = rate_up_ * dt;
        const T max_down = rate_down_ * dt;
        const T delta = target - y_;
        if (delta > max_up) {
            y_ += max_up;
        } else if (delta < -max_down) {
            y_ -= max_down;
        } else {
            y_ = target;
        }
        return y_;
    }

    [[nodiscard]] constexpr T value() const { return y_; }

    constexpr void reset() {
        y_ = T{0};
        init_ = false;
    }

    /// Pre-seed the output (and mark initialized) for bumpless start.
    constexpr void reset(T y) {
        y_ = y;
        init_ = true;
    }

private:
    T    rate_up_{std::numeric_limits<T>::max()};
    T    rate_down_{std::numeric_limits<T>::max()};
    T    y_{T{0}};
    bool init_{false};
};

/**
 * @brief Hysteresis comparator (Schmitt trigger): bool output with separate
 *        on/off thresholds to reject chatter.
 *
 * Output latches true once the input rises above @p high and stays true until
 * it falls below @p low (`low ≤ high`); between the thresholds it holds. The
 * standard cure for relay chatter on a noisy threshold crossing.
 *
 * @tparam T Scalar type.
 */
template<typename T = float>
class Hysteresis {
public:
    constexpr Hysteresis() = default;

    /// @param low Falling (release) threshold. @param high Rising (trip) threshold.
    constexpr Hysteresis(T low, T high) : low_(low), high_(high) {}

    /// Update with one sample, return the latched state.
    constexpr bool operator()(T x) {
        if (state_) {
            if (x < low_) {
                state_ = false;
            }
        } else if (x > high_) {
            state_ = true;
        }
        return state_;
    }

    [[nodiscard]] constexpr bool state() const { return state_; }

    constexpr void reset(bool state = false) { state_ = state; }

private:
    T    low_{0};
    T    high_{0};
    bool state_{false};
};

/**
 * @brief Classification of an analog input against its valid/fault bands.
 *
 * Models NAMUR NE43-style signal-range monitoring: the live-zero margins of a
 * ratiometric (0.5–4.5 V) or 4–20 mA sensor let a genuine reading be told apart
 * from a wire fault.
 */
enum class SignalStatus : std::uint8_t {
    Valid,      ///< within the nominal measuring span
    UnderRange, ///< below span but not a wire fault (saturated low)
    OverRange,  ///< above span but not a wire fault (saturated high)
    FaultLow,   ///< below the fault floor — open circuit / short to ground
    FaultHigh,  ///< above the fault ceiling — short to supply
};

/// True only for the in-span status.
[[nodiscard]] constexpr bool is_valid(SignalStatus s) {
    return s == SignalStatus::Valid;
}

/// True for a wire fault (FaultLow/FaultHigh) — i.e. not a real reading at all.
[[nodiscard]] constexpr bool is_fault(SignalStatus s) {
    return s == SignalStatus::FaultLow || s == SignalStatus::FaultHigh;
}

/**
 * @brief Classify @p x against the four band edges `[fault_lo (valid_lo,
 *        valid_hi) fault_hi]` (assumed ordered, non-decreasing).
 *
 * Bands: `x < fault_lo` → FaultLow; `[fault_lo, valid_lo)` → UnderRange;
 * `[valid_lo, valid_hi]` → Valid; `(valid_hi, fault_hi]` → OverRange;
 * `x > fault_hi` → FaultHigh. Pure/stateless.
 */
template<typename T>
[[nodiscard]] constexpr SignalStatus
classify_range(T x, T fault_lo, T valid_lo, T valid_hi, T fault_hi) {
    if (x < fault_lo) {
        return SignalStatus::FaultLow;
    }
    if (x > fault_hi) {
        return SignalStatus::FaultHigh;
    }
    if (x < valid_lo) {
        return SignalStatus::UnderRange;
    }
    if (x > valid_hi) {
        return SignalStatus::OverRange;
    }
    return SignalStatus::Valid;
}

/**
 * @brief Analog-input range/fault monitor (NAMUR NE43 pattern).
 *
 * Classifies a reading into valid / out-of-range / wire-fault bands so a broken
 * wire (short to ground/supply) is distinguished from a real low/high reading —
 * the standard front-end check before scaling a 0.5–4.5 V or 4–20 mA sensor.
 *
 * With a nonzero @p hysteresis the current band is widened by that margin, so a
 * signal hovering on an edge must move `hysteresis` past it to re-classify (no
 * boundary chatter). `hysteresis = 0` (default) is a plain per-call classifier.
 * For *time*-based fault qualification, compose `faulted()` with a `Debounce`.
 *
 * @code
 * RangeMonitor<float> ai{0.25f, 0.5f, 4.5f, 4.75f}; // [0.25 (0.5 4.5) 4.75]
 * if (ai.faulted(v)) trip();   // broken wire, not a low reading
 * @endcode
 */
template<typename T = float>
class RangeMonitor {
public:
    constexpr RangeMonitor() = default;

    /// @param fault_lo,valid_lo,valid_hi,fault_hi The four band edges.
    /// @param hysteresis Re-classification margin (≥ 0; 0 = none).
    constexpr RangeMonitor(T fault_lo, T valid_lo, T valid_hi, T fault_hi, T hysteresis = T{0})
        : fault_lo_(fault_lo), valid_lo_(valid_lo), valid_hi_(valid_hi), fault_hi_(fault_hi), hyst_(hysteresis) {}

    /// Classify one sample (applies hysteresis around the previous status).
    constexpr SignalStatus operator()(T x) {
        // Bias the current band's edges outward by hyst_ so leaving it costs an
        // extra margin; entering bands keep their nominal edges.
        T lo_f = fault_lo_;
        T lo_v = valid_lo_;
        T hi_v = valid_hi_;
        T hi_f = fault_hi_;
        switch (status_) {
            case SignalStatus::Valid:
                lo_v -= hyst_;
                hi_v += hyst_;
                break;
            case SignalStatus::UnderRange:
                lo_f -= hyst_;
                lo_v += hyst_;
                break;
            case SignalStatus::OverRange:
                hi_v -= hyst_;
                hi_f += hyst_;
                break;
            case SignalStatus::FaultLow:
                lo_f += hyst_;
                break;
            case SignalStatus::FaultHigh:
                hi_f -= hyst_;
                break;
        }
        status_ = classify_range(x, lo_f, lo_v, hi_v, hi_f);
        return status_;
    }

    /// Convenience predicates that also advance the monitor with @p x.
    constexpr bool valid(T x) { return is_valid((*this)(x)); }
    constexpr bool faulted(T x) { return is_fault((*this)(x)); }

    [[nodiscard]] constexpr SignalStatus value() const { return status_; }
    constexpr void                       reset() { status_ = SignalStatus::Valid; }

private:
    T            fault_lo_{};
    T            valid_lo_{};
    T            valid_hi_{};
    T            fault_hi_{};
    T            hyst_{T{0}};
    SignalStatus status_{SignalStatus::Valid};
};

/**
 * @brief How two redundant analog channels should relate in engineering units.
 *
 * After each channel is calibrated, the pair is either:
 * - SameSlope — both read the same quantity (@f$ a \approx b @f$), or
 * - OppositeSlope — complementary reads that sum to a known span
 *   (@f$ a + b \approx S @f$, e.g. dual 0–100% sensors wired with opposite slope).
 *
 * Opposite-slope sensors that are already inverted in their AffineCal (so both
 * report the same engineering direction) use SameSlope after calibration.
 */
enum class AnalogCrossMode : std::uint8_t {
    SameSlope,     ///< @f$ a \approx b @f$
    OppositeSlope, ///< @f$ a + b \approx \mathrm{span} @f$
};

/**
 * @brief Result of comparing two redundant analog channels.
 */
enum class AnalogCrossStatus : std::uint8_t {
    Agree,        ///< both usable and within tolerance
    Disagree,     ///< both usable but disagree (no silent average as “truth”)
    ChannelAOnly, ///< A usable, B not — value from A (degraded)
    ChannelBOnly, ///< B usable, A not — value from B (degraded)
    BothUnusable, ///< neither channel is Valid/over/under usable
};

/**
 * @brief Output of @ref cross_check_analog.
 *
 * @tparam T Scalar type.
 */
template<typename T = float>
struct AnalogCrossResult {
    T                 value{};          ///< fused or single-channel engineering value when @ref has_value
    bool              has_value{false}; ///< true when a usable value is provided
    AnalogCrossStatus status{AnalogCrossStatus::BothUnusable};
    SignalStatus      status_a{SignalStatus::Valid};
    SignalStatus      status_b{SignalStatus::Valid};

    [[nodiscard]] constexpr bool agree() const { return status == AnalogCrossStatus::Agree; }
    [[nodiscard]] constexpr bool degraded() const {
        return status == AnalogCrossStatus::ChannelAOnly || status == AnalogCrossStatus::ChannelBOnly;
    }
};

/**
 * @brief True when a channel reading may still be used (in-span or saturated).
 *
 * Wire faults are excluded. Over/under range are treated as usable-but-saturated
 * for dual-channel voting (the peer can still disagree on the engineering value).
 */
[[nodiscard]] constexpr bool is_usable(SignalStatus s) {
    return s == SignalStatus::Valid || s == SignalStatus::UnderRange || s == SignalStatus::OverRange;
}

/**
 * @brief Cross-check two calibrated analog readings for redundant sensors.
 *
 * Validity is decided first (per-channel @ref SignalStatus); only then is the
 * same-slope or opposite-slope relationship tested. When both channels are usable
 * and agree, the returned value is the average (same-slope) or channel A (opposite
 * slope, equivalent to @f$ (a + (S-b))/2 @f$ when they agree). When they disagree,
 * @ref AnalogCrossResult::has_value is false — the caller must not treat a mid-point
 * of disagreeing sensors as truth. When only one channel is usable, that channel's
 * value is returned with a degraded status.
 *
 * @param eng_a   Channel A in engineering units (after affine cal).
 * @param status_a Channel A range/fault class (from @ref RangeMonitor / AnalogInput).
 * @param eng_b   Channel B in engineering units.
 * @param status_b Channel B range/fault class.
 * @param mode    @ref AnalogCrossMode::SameSlope or @ref AnalogCrossMode::OppositeSlope.
 * @param tol     Allowed absolute mismatch (@f$ |a-b| @f$ or @f$ |a+b-S| @f$).
 * @param span    Expected @f$ a+b @f$ for opposite-slope mode (ignored for same-slope).
 * @return Cross-check result with optional fused value.
 *
 * @code
 * // Two pressure sensors, same slope, 2 bar agreement window
 * auto r = cross_check_analog(p_a, sa, p_b, sb, AnalogCrossMode::SameSlope, 2.0f);
 * if (r.agree()) use(r.value);
 *
 * // Complementary 0–100 % pair (opposite slope after cal): a + b ≈ 100
 * auto r2 = cross_check_analog(pct_a, sa, pct_b, sb, AnalogCrossMode::OppositeSlope, 1.0f, 100.0f);
 * @endcode
 */
template<typename T>
[[nodiscard]] constexpr AnalogCrossResult<T> cross_check_analog(
    T               eng_a,
    SignalStatus    status_a,
    T               eng_b,
    SignalStatus    status_b,
    AnalogCrossMode mode,
    T               tol,
    T               span = T{0}
) {
    const bool a_ok = is_usable(status_a);
    const bool b_ok = is_usable(status_b);

    if (!a_ok && !b_ok) {
        return AnalogCrossResult<T>{
            .status = AnalogCrossStatus::BothUnusable,
            .status_a = status_a,
            .status_b = status_b,
        };
    }
    if (a_ok && !b_ok) {
        return AnalogCrossResult<T>{
            .value = eng_a,
            .has_value = true,
            .status = AnalogCrossStatus::ChannelAOnly,
            .status_a = status_a,
            .status_b = status_b,
        };
    }
    if (!a_ok && b_ok) {
        return AnalogCrossResult<T>{
            .value = eng_b,
            .has_value = true,
            .status = AnalogCrossStatus::ChannelBOnly,
            .status_a = status_a,
            .status_b = status_b,
        };
    }

    // Both usable — test the geometric relationship.
    T mismatch{};
    if (mode == AnalogCrossMode::SameSlope) {
        mismatch = damp::abs(eng_a - eng_b);
    } else {
        mismatch = damp::abs((eng_a + eng_b) - span);
    }

    if (mismatch <= tol) {
        const T value = (mode == AnalogCrossMode::SameSlope)
                          ? ((eng_a + eng_b) * static_cast<T>(0.5))
                          // A and (span − B) both estimate the quantity; average them.
                          : ((eng_a + (span - eng_b)) * static_cast<T>(0.5));
        return AnalogCrossResult<T>{
            .value = value,
            .has_value = true,
            .status = AnalogCrossStatus::Agree,
            .status_a = status_a,
            .status_b = status_b,
        };
    }

    return AnalogCrossResult<T>{
        .has_value = false,
        .status = AnalogCrossStatus::Disagree,
        .status_a = status_a,
        .status_b = status_b,
    };
}

} // namespace damp
