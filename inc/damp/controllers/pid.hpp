// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file pid.hpp
 * @brief PID design and runtime controllers
 */


#include "damp/backend.hpp"
#include "damp/controllers/detail/anti_windup.hpp"
#include "damp/controllers/law.hpp" // IWYU pragma: export
#include "damp/design/status.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/toolbox/bounds.hpp"

// Arduino.h defines PI as a macro, which would rewrite PIDMode::PI: set aside for this header.
#pragma push_macro("PI")
#undef PI

namespace damp {

/**
 * @brief Compile-time selection of the PID control-law structure.
 *
 * Picks which terms PIDController computes: the P term always, the integral with PI and PID,
 * the filtered derivative with PD and PID. A design with a gain the mode does not compute is
 * refused.
 */
enum class PIDMode : uint8_t {
    P,
    PI,
    PD,
    PID,
};

namespace detail {

/**
 * @brief True when a PID gain set has no non-finite entry, a non-negative derivative filter and a
 *        non-negative tracking time
 *
 * A Tf < 0 puts the filtered-derivative pole Tf/(Ts + Tf) outside the unit circle for
 * Tf < −Ts/2, and divides by zero at Tf = −Ts. A tracking time T_t < 0 would push the integral
 * further into the limit it is to unwind from.
 */
template<typename T>
[[nodiscard]] constexpr bool pid_gains_in_domain(T Kp, T Ki, T Kd, T Tf, T tracking_time, T b, T c) {
    return damp::isfinite(Kp) && damp::isfinite(Ki) && damp::isfinite(Kd) && damp::finite_non_negative(Tf)
        && damp::finite_non_negative(tracking_time) && damp::isfinite(b) && damp::isfinite(c);
}

/** @brief Whether @p Mode computes every term the gains ask for: no Kᵢ without an integral, no K_d without a derivative. */
template<PIDMode Mode, typename T>
[[nodiscard]] constexpr bool pid_mode_realizes(T Ki, T Kd) {
    constexpr bool has_i = (Mode == PIDMode::PI) || (Mode == PIDMode::PID);
    constexpr bool has_d = (Mode == PIDMode::PD) || (Mode == PIDMode::PID);
    return (has_i || Ki == T{0}) && (has_d || Kd == T{0});
}

/**
 * @brief The per-tick coefficients of a PID at a fixed sample period, computed once when the gains
 *        or the period are set, so the fixed-rate tick divides by nothing
 */
template<typename T>
struct PIDCoefficients {
    T ki_ts{T{0}};  ///< Kᵢ·T_s: the integral's increment per unit error
    T d_a{T{0}};    ///< T_f/(T_s + T_f): the filtered derivative's pole
    T d_b{T{0}};    ///< K_d/(T_s + T_f): its gain on Δ(c·r − y)
    T inv_tt{T{0}}; ///< 1/T_t, 0 for T_t = 0 (conditional integration)

    /** @brief Every coefficient finite (a gain times T_s, or 1/T_t, can overflow). */
    [[nodiscard]] constexpr bool finite() const {
        return damp::isfinite(ki_ts) && damp::isfinite(d_a) && damp::isfinite(d_b) && damp::isfinite(inv_tt);
    }
};

/**
 * @brief The coefficients of the gains @p Ki, @p Kd, @p Tf, @p tracking_time at the period @p Ts
 *
 * Backward-Euler integral, I ← I + Kᵢ·T_s·e, and filtered derivative,
 * D ← T_f/(T_s + T_f)·D + K_d/(T_s + T_f)·Δx. With K_d = 0 the derivative is off (both zero).
 * A @p Ts that is not positive (a law with no fixed period, ticked by step(r, y, dt) only) gives
 * no increment and no derivative.
 */
template<typename T>
[[nodiscard]] constexpr PIDCoefficients<T> pid_coefficients(T Ki, T Kd, T Tf, T tracking_time, T Ts) {
    PIDCoefficients<T> k{};
    if (tracking_time > T{0}) {
        k.inv_tt = T{1} / tracking_time;
    }
    if (!damp::finite_positive(Ts)) {
        return k;
    }
    k.ki_ts = Ki * Ts;
    if (Kd != T{0}) {
        const T den = Ts + Tf;
        k.d_a = Tf / den;
        k.d_b = Kd / den;
    }
    return k;
}

/** @brief A PIDController's integral memory: its limits and I (absent on a P or a PD). */
template<typename T>
struct PIDIntegral {
    Bounds<1, T> lim{};       ///< The limits on I (u's units)
    T            value{T{0}}; ///< The integral term I (u's units)
    T            step{T{0}};  ///< What the last tick added to I: what a clamp takes back when T_t = 0
};

/** @brief A PIDController's derivative memory (absent on a P or a PI). */
template<typename T>
struct PIDDerivative {
    T    prev_x{T{0}}; ///< Previous c·r − y
    T    value{T{0}};  ///< The filtered derivative term D
    bool first{true};  ///< The next tick seeds prev_x (no derivative kick)
};

/**
 * @brief A PIDController's gains, its fixed sample period and the coefficients at it: the members
 *        a PIDController holds, or the constants of a `constexpr` design's Law
 *
 * Only the terms @p Mode computes have storage.
 */
template<typename T, PIDMode Mode>
struct PIDGainSet {
    static constexpr bool constant = false; ///< Members of the runtime, settable

    /** @brief The integral term's gain, tracking time and their coefficients. */
    struct Integral {
        T Ki{T{0}};     ///< Continuous integral gain
        T tt{T{0}};     ///< [s] Tracking time T_t (0: conditional integration)
        T ki_ts{T{0}};  ///< Kᵢ·T_s: the increment per unit error at the fixed rate
        T inv_tt{T{0}}; ///< 1/T_t, 0 for T_t = 0
    };

    /** @brief The derivative term's gain, filter, setpoint weight and coefficients. */
    struct Derivative {
        T Kd{T{0}}; ///< Derivative gain
        T Tf{T{0}}; ///< [s] Filter time constant
        T c{T{1}};  ///< Derivative setpoint weight
        T a{T{0}};  ///< T_f/(T_s + T_f): the filter's pole at the fixed rate
        T b{T{0}};  ///< K_d/(T_s + T_f): its gain on Δ(c·r − y) at the fixed rate
    };

    static constexpr bool kI = (Mode == PIDMode::PI) || (Mode == PIDMode::PID); ///< The mode integrates
    static constexpr bool kD = (Mode == PIDMode::PD) || (Mode == PIDMode::PID); ///< The mode differentiates

    T                                                                          Kp{T{0}}; ///< Proportional gain
    T                                                                          b{T{1}};  ///< Proportional setpoint weight
    T                                                                          Ts{T{0}}; ///< [s] The fixed-rate tick's period (0: none)
    [[no_unique_address]] maybe<kI, Integral, OffSlot::pid_integral_gains>     i{};      ///< The integral's gains: none on a P or a PD
    [[no_unique_address]] maybe<kD, Derivative, OffSlot::pid_derivative_gains> d{};      ///< The derivative's gains: none on a P or a PI
};

} // namespace detail

namespace design {

/**
 * @struct PIDLimits
 * @brief The limits of a PID law: on its command and on its integral term
 *
 * Unbounded by default. design::PIDGains::discretize takes the ones a law starts with, and the law
 * reads and changes them at run time with PIDController::limits and PIDController::set_limits.
 *
 * @p integral clamps the integral term I, which enters u with unit gain, so it is in u's units.
 * That differs from StateFeedbackLimits::integral, which bounds the LQ integrators ξ ≈ ∫e dt in
 * the units of e times seconds.
 */
template<typename T = double>
struct PIDLimits {
    Bounds<1, T> input{};    ///< The command u [lower, upper]
    Bounds<1, T> integral{}; ///< The integral term I, in u's units

    /** @brief The same limits in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr PIDLimits<U> as() const {
        return PIDLimits<U>{input.template as<U>(), integral.template as<U>()};
    }
};

template<typename T>
struct PIDResult;

/**
 * @struct PIDGains
 * @brief 2-DOF continuous-time PID gains, setpoint weights, derivative filter and tracking time
 *
 * An aggregate: every field has a safe default (no back-calculation, b = c = 1, unfiltered
 * derivative), so name only the knobs you set — `design::PIDGains{.Kp = 4.0, .Ki = 8.0}`. The
 * limits are not gains: discretize takes them, `.discretize(0.01, {.input = {.lower = -12.0, .upper = 12.0}})`,
 * and the law changes them at run time with set_limits. An output limit is enough for
 * anti-windup: with tracking_time = 0 the integral holds against the clamp (conditional
 * integration), its own or one reported with back_calculate. A tracking_time T_t > 0 unwinds it
 * by back-calculation instead, a share min(T_s/T_t, 1) of the clamp's excess each tick.
 * @ref design::pid covers the plain (Kp, Ki, Kd) case. Deploy via discretize into a
 * @ref PIDResult for PIDController, which ticks at that period with step(r, y) or at a measured
 * one with step(r, y, dt), and takes new gains at run time with set_gains.
 *
 * Special cases:
 * - b=1, c=1: Standard PID (P and D on error)
 * - b=1, c=0: PI-D (D on measurement - no derivative kick)
 * - b=0, c=0: I-PD (P and D on measurement - no setpoint kick)
 *
 * @see Astrom & Hagglund, "Advanced PID Control" (2006), Sec. 4.4
 */
template<typename T = double>
struct [[nodiscard]] PIDGains {
    T Kp{};
    T Ki{};
    T Kd{};
    T tracking_time{T{0}}; ///< [s] Back-calculation tracking time T_t: a clamp's excess unwinds the integral by min(T_s/T_t, 1) of it each tick.
                           ///< 0 (default): no unwinding, I holds against a clamp (conditional integration). T_t = T_i = Kp/Ki is the usual choice.
    T b = T{1};            ///< Proportional setpoint weight (0=I-PD, 1=standard PID)
    T c = T{1};            ///< Derivative setpoint weight   (0=PI-D, 1=standard PID)
    T Tf = T{0};           ///< Derivative filter time constant (0 = unfiltered)

    bool         success{true};            ///< status == DesignStatus::ok (a hand-built set is a valid design)
    DesignStatus status{DesignStatus::ok}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The same gains in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return PIDGains<U>{
            static_cast<U>(Kp),
            static_cast<U>(Ki),
            static_cast<U>(Kd),
            static_cast<U>(tracking_time),
            static_cast<U>(b),
            static_cast<U>(c),
            static_cast<U>(Tf),
            success,
            status,
        };
    }

    /**
     * @brief The deploy form: these gains at the sample period @p ts, with the limits the law starts with
     *
     * @code
     * constinit PIController controller =
     *     design::pid(4.0, 8.0, 0.0).discretize(0.01).as<float>();
     * @endcode
     *
     * A failed design stays failed (its status is carried over), a @p ts that is not
     * positive and finite is refused with DesignStatus::invalid_sample_time, and a
     * non-finite gain, a Tf < 0, a tracking_time that is negative or not finite, or a limit
     * that is not finite or whose lower end is above its upper one is refused with
     * DesignStatus::invalid_parameter.
     *
     * @param ts     Sample period [s]
     * @param limits The command's and the integral term's limits (unbounded by default)
     * @return PIDResult for PIDController
     */
    [[nodiscard]] constexpr PIDResult<T> discretize(T ts, const PIDLimits<T>& limits = {}) const;

    /**
     * @brief Error-channel C(s) = Kp + Ki/s + Kd s/(1 + Tf s)
     *
     * Feedback law (y → −u at r = 0). Bode / series / feedback of L use this
     * map. Setpoint weights b, c do not appear. Tf = 0 is the ideal
     * (Kd s² + Kp s + Ki)/s.
     */
    [[nodiscard]] constexpr TransferFunction<3, 3, T> to_tf() const {
        return TransferFunction<3, 3, T>{
            .num = {(Kp * Tf) + Kd, Kp + (Ki * Tf), Ki}, // (Kp Tf + Kd) s² + (Kp + Ki Tf) s + Ki
            .den = {Tf, T{1}, T{0}},                     // Tf s² + s
        };
    }

    /**
     * @brief Setpoint-channel C_ff(s) = b Kp + Ki/s + c Kd s/(1 + Tf s)
     *
     * Tracking is T_yr = C_ff P / (1 + C_fb P). Equals @ref to_tf when b = c = 1.
     */
    [[nodiscard]] constexpr TransferFunction<3, 3, T> to_tf_ff() const {
        return TransferFunction<3, 3, T>{
            .num = {(b * Kp * Tf) + (c * Kd), (b * Kp) + (Ki * Tf), Ki},
            .den = {Tf, T{1}, T{0}},
        };
    }

    /**
     * @brief Error-channel C_fb(jω) for L = C_fb P
     */
    [[nodiscard]] constexpr damp::complex<T> eval_fb(T w) const {
        return eval_channel(w, T{1}, T{1});
    }

    /**
     * @brief Setpoint-channel C_ff(jω) (P weight b, D weight c)
     */
    [[nodiscard]] constexpr damp::complex<T> eval_ff(T w) const {
        return eval_channel(w, b, c);
    }

    /**
     * @brief Error-channel realization C(s) = Kₚ + Kᵢ/s + K_d s/(1 + T_f s), the map of @ref to_tf
     *
     * States: the integral of e, and the derivative filter state, as in @ref to_ss_2dof
     * (its error column). With K_d = 0 the filter state has no input and no output.
     *
     * @return The realization, or nullopt when K_d ≠ 0 and T_f ≤ 0 (an unfiltered derivative is not proper)
     */
    [[nodiscard]] constexpr damp::optional<StateSpace<2, 1, 1, T>> to_ss() const {
        const auto two = to_ss_2dof();
        if (!two) {
            return damp::nullopt;
        }
        // e = r − y enters as −y at r = 0: negate the y column.
        StateSpace<2, 1, 1, T> out{};
        out.A = two->A;
        out.B(0, 0) = -two->B(0, 1);
        out.B(1, 0) = -two->B(1, 1);
        out.C = two->C;
        out.D(0, 0) = -two->D(0, 1);
        return out;
    }

    /**
     * @brief Two-degree-of-freedom realization u = K[r; y] with the setpoint weights, for analysis::gang_of_seven
     *
     * States: the integral of r − y, and the derivative filter state of c·r − y through
     * 1/(1 + T_f s). So u = b·Kₚ r − Kₚ y + Kᵢ x_I + (K_d/T_f)(c·r − y − x_D), the setpoint
     * channel @ref to_tf_ff from r and the error channel @ref to_tf from −y. With K_d = 0 the
     * filter state has no input and no output; it sits at −1/T_f (−1 when T_f = 0) and cancels
     * from every transfer function.
     *
     * @return The realization, or nullopt when K_d ≠ 0 and T_f ≤ 0 (an unfiltered derivative is not proper)
     */
    [[nodiscard]] constexpr damp::optional<StateSpace<2, 2, 1, T>> to_ss_2dof() const {
        if ((Kd != T{0}) && !damp::finite_positive(Tf)) {
            return damp::nullopt;
        }
        StateSpace<2, 2, 1, T> out{};
        out.B(0, 0) = T{1};
        out.B(0, 1) = T{-1};
        out.C(0, 0) = Ki;
        out.D(0, 0) = b * Kp;
        out.D(0, 1) = -Kp;
        if (Kd == T{0}) {
            out.A(1, 1) = Tf > T{0} ? -T{1} / Tf : T{-1};
            return out;
        }
        const T a = T{1} / Tf;
        out.A(1, 1) = -a;
        out.B(1, 0) = c * a;
        out.B(1, 1) = -a;
        out.C(0, 1) = -Kd * a;
        out.D(0, 0) += c * Kd * a;
        out.D(0, 1) -= Kd * a;
        return out;
    }

private:
    [[nodiscard]] constexpr damp::complex<T> eval_channel(T w, T bp, T cd) const {
        using Cplx = damp::complex<T>;
        Cplx out{bp * Kp, T{0}};
        if ((w > T{0}) && (Ki != T{0})) {
            out = out + (Cplx{Ki, T{0}} / Cplx{T{0}, w});
        }
        if ((w > T{0}) && (Kd != T{0}) && (cd != T{0})) {
            const Cplx s{T{0}, w};
            out = out + ((cd * Kd) * s) / (Cplx{T{1}, T{0}} + (Tf * s));
        }
        return out;
    }
};

/**
 * @struct PIDResult
 * @brief A PID design deployed at a sample period: the continuous gains, the limits the law starts with, and T_s
 *
 * Produced by PIDGains::discretize. The runtime computes its per-tick coefficients from the
 * continuous gains and T_s (backward-Euler integral, filtered derivative):
 * - @f$I \leftarrow I + K_i T_s\,e@f$, @f$u_I = I@f$ (same units as @f$u@f$), clamped into
 *   limits.integral
 * - @f$D \leftarrow a_d D + b_d\,\Delta(c r - y)@f$ with
 *   @f$a_d = T_f/(T_s+T_f)@f$, @f$b_d = K_d/(T_s+T_f)@f$
 *   (@f$T_f=0@f$ ⇒ raw @f$K_d/T_s@f$ difference)
 * - Back-calculation: @f$\Delta I = \min(T_s/T_t, 1)\,(u_{\mathrm{sat}}-u_{\mathrm{unsat}})@f$
 *
 * @see PIDGains::discretize
 * @see Åström & Hägglund, "Advanced PID Control" (2006), Sec. 3.3
 */
template<typename T = double>
struct [[nodiscard]] PIDResult {
    PIDGains<T>  gains{};                             ///< The continuous gains
    PIDLimits<T> limits{};                            ///< The limits the law starts with
    T            Ts{T{0}};                            ///< [s] Sample period of the fixed-rate tick
    bool         success{false};                      ///< status == DesignStatus::ok
    DesignStatus status{DesignStatus::not_converged}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The same design in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr PIDResult<U> as() const {
        return PIDResult<U>{gains.template as<U>(), limits.template as<U>(), static_cast<U>(Ts), success, status};
    }
};

namespace detail {
using ::damp::detail::pid_gains_in_domain;
} // namespace detail

template<typename T>
constexpr PIDResult<T> PIDGains<T>::discretize(T ts, const PIDLimits<T>& limits) const {
    if (!success) {
        return detail::refused<PIDResult<T>>(status);
    }
    if (!damp::finite_positive(ts)) {
        return detail::refused<PIDResult<T>>(DesignStatus::invalid_sample_time);
    }
    if (!detail::pid_gains_in_domain(Kp, Ki, Kd, Tf, tracking_time, b, c)) {
        return detail::refused<PIDResult<T>>(DesignStatus::invalid_parameter);
    }
    if (!limits.input.valid() || !limits.integral.valid()) {
        return detail::refused<PIDResult<T>>(DesignStatus::invalid_parameter); // as the LQ and MPC designs
    }
    PIDResult<T> out{.gains = *this, .limits = limits, .Ts = ts};
    detail::set_status(out, DesignStatus::ok);
    return out;
}

/**
 * @brief Continuous PID gains Kₚ, Kᵢ, K_d with every other knob at its default
 *
 * Short form for the common case: no back-calculation, b = c = 1, unfiltered derivative. For
 * setpoint weights, a derivative filter or a tracking time, build @ref PIDGains directly with
 * designated initializers. Deploy with PIDGains::discretize(Ts, limits) for PIDController.
 *
 * The control law is:
 *
 *     u = Kp(b*r - y) + I + Kd*d/dt(c*r - y)
 *     dI/dt = Ki·(r - y)
 *
 * Preferred deploy construction:
 *
 * @code
 * constinit PIController controller =
 *     design::pid(4.0, 8.0, 0.0).discretize(0.01).as<float>();
 *
 * constinit PIController clamped =
 *     design::pid(4.0, 8.0, 0.0).discretize(0.01, {.input = {.lower = -12.0, .upper = 12.0}}).as<float>();
 * @endcode
 *
 * @param Kp Proportional gain
 * @param Ki Integral gain (continuous)
 * @param Kd Derivative gain
 *
 * @return PIDGains{Kp, Ki, Kd} with default weights, filter and tracking time
 */
template<typename T = double>
[[nodiscard]] constexpr PIDGains<T> pid(T Kp, T Ki, T Kd) {
    return PIDGains<T>{.Kp = Kp, .Ki = Ki, .Kd = Kd};
}

} // namespace design

namespace detail {

/**
 * @brief The gain set of @p g at the fixed period @p Ts for @p Mode, or nullopt when a gain is out
 *        of its domain, the mode does not compute it (K_d ≠ 0 on a PI, Kᵢ ≠ 0 on a P or a PD), or
 *        a coefficient overflows
 */
template<PIDMode Mode, typename T>
[[nodiscard]] constexpr damp::optional<PIDGainSet<T, Mode>> pid_gain_set(const design::PIDGains<T>& g, T Ts) {
    if (!pid_gains_in_domain(g.Kp, g.Ki, g.Kd, g.Tf, g.tracking_time, g.b, g.c) || !pid_mode_realizes<Mode>(g.Ki, g.Kd)) {
        return damp::nullopt;
    }
    const auto k = pid_coefficients(g.Ki, g.Kd, g.Tf, g.tracking_time, Ts);
    if (!k.finite()) {
        return damp::nullopt;
    }
    PIDGainSet<T, Mode> set{};
    set.Kp = g.Kp;
    set.b = g.b;
    set.Ts = Ts;
    if constexpr (PIDGainSet<T, Mode>::kI) {
        set.i = {.Ki = g.Ki, .tt = g.tracking_time, .ki_ts = k.ki_ts, .inv_tt = k.inv_tt};
    }
    if constexpr (PIDGainSet<T, Mode>::kD) {
        set.d = {.Kd = g.Kd, .Tf = g.Tf, .c = g.c, .a = k.d_a, .b = k.d_b};
    }
    return set;
}

} // namespace detail

/**
 * @ingroup discrete_controllers
 * @brief 2-DOF PID runtime: fixed-rate step(r, y), measured-period step(r, y, dt), gains and limits settable at run time
 *
 * One runtime for the four laws: @ref PIDMode picks the terms the tick computes (P, PI, PD, PID),
 * and the aliases PController, PIController and PDController name the narrower ones. A design
 * with a gain its mode does not compute (K_d ≠ 0 on a PI, Kᵢ ≠ 0 on a P or a PD), a failed
 * design, a non-finite gain, a T_f or tracking time below 0, a T_s that is not positive or a
 * limit that is not a box is refused at construction: ill-formed under constant evaluation
 * (`constinit` deploy); at run time it leaves the controller invalid (valid() false, output 0),
 * as every law does.
 *
 * | Member | What it is |
 * | --- | --- |
 * | step(r, y) | one tick at the design's T_s, on coefficients computed when the gains were set |
 * | step(r, y, dt) | one tick over a measured period dt [s] (jitter, a variable-rate loop) |
 * | track(u, r, y), track(u, r, y, dt) | in place of each step while another source drives the plant |
 * | back_calculate(u_sat, u_ff) | a clamp after the law: the input the plant received |
 * | gains(), set_gains(g) | the continuous gains in force (design::PIDGains); change them at run time |
 * | limits(), set_limits(lim) | the command's and the integral's limits (design::PIDLimits); change them at run time |
 * | integral(), derivative(), set_integral(v) | the integral term I and the filtered derivative D, both in u's units |
 * | last_command(), sample_time(), valid(), reset() | the last command; T_s; whether the design ran; clear I and D |
 *
 * Control law, at T_s or at dt in its place:
 * @f[
 *   I \leftarrow I + K_i T_s\,e,\quad
 *   D \leftarrow \tfrac{T_f}{T_s + T_f} D + \tfrac{K_d}{T_s + T_f}\,\Delta(c r - y),\quad
 *   u = \mathrm{sat}\big(K_p(b r - y) + I + D\big)
 * @f]
 * The first tick seeds Δ with zero (no derivative kick). The integral is clamped into
 * limits().integral, the command into limits().input, and the integral holds against that clamp
 * (it gives back what the tick added, up to what the clamp cut off: conditional integration), or
 * with a tracking time T_t > 0 unwinds by back-calculation, ΔI = min(T_s/T_t, 1)·(u − u_unsat).
 * The fixed-rate tick divides by nothing: its coefficients are computed when the gains are set,
 * and 1/T_t with them. The measured tick takes one reciprocal, 1/(dt + T_f), for the derivative,
 * and none without one.
 *
 * Gains and limits are tunable at run time, get-modify-set, and take effect at the next tick:
 * @code
 * auto g = pid.gains();
 * g.Kp = 2.5f;
 * if (!pid.set_gains(g)) { } // refused: nothing changed
 * auto lim = pid.limits();
 * lim.input = {.lower = -v_bus, .upper = v_bus};
 * if (!pid.set_limits(lim)) { }
 * @endcode
 * The integral is held as the I term, in u's units, not as ∫e, so a new Kᵢ changes its rate
 * and not its value: retuning does not bump the command.
 *
 * While another source drives the plant (an operator, another law), track(u, r, y) takes the
 * place of step(r, y): the tick's derivative memory, then the integral preloaded so that
 * step(r, y) at this r and y returns u clamped into the input limits, and u recorded as the
 * command applied. There is no mode: the law tracks for as long as the caller calls track, and the
 * first step after it hands over. A P or a PD has nothing to preload; track keeps its derivative
 * memory current, so the first step has no kick.
 *
 * @code
 * constinit auto c = PIDController{
 *     design::pid(2.0, 5.0, 0.1).discretize(0.01).as<float>()};
 * float u = c.step(r, y);
 * @endcode
 *
 * A law only ever ticked by step(r, y, dt) needs no T_s: default-construct it and set_gains,
 * set_limits. Its step(r, y) then has no period to integrate or differentiate over: the
 * integral and the derivative hold.
 *
 * Angle input: pass r′ = y + wrap_pi(r − y), the reference nearest the measurement,
 * `pid.step(y + wrap_pi(r - y), y)`, so P and I see the short-way error. With c ≠ 1
 * (derivative on measurement) D also reads y, so unwrap y first with @ref PhaseUnwrapper.
 * Keep b = 1 on an angle axis: with b ≠ 1, P carries −Kₚ(1 − b)·y, which depends on the
 * absolute angle and the turn count rather than on the error.
 *
 * MATLAB® equivalent: design with `pid(...)` then `c2d`; runtime is the discrete law.
 * A profiled PID takes its setpoint from a motion profile each tick (the motion pack's
 * trapezoidal profile, not yet public).
 * @see Åström & Hägglund, "Advanced PID Control" (2006), Sec. 3.3–4.4, §3.5 (anti-windup, bumpless transfer)
 *
 * @note Gains from a continuous-time designer (ziegler_nichols, cohen_coon,
 *       simc, lambda_tuning, pidtune) are Ts-agnostic until discretize(Ts).
 *       Gains from the sampled design::pid_pole_placement(plant, poles, Ts) already
 *       embed that Ts in continuous form — discretize with the same Ts.
 *
 * @tparam T     Scalar
 * @tparam Mode  The terms the tick computes
 * @tparam Gains Where the gains live: void, the default, for members (detail::PIDGainSet, settable), or the
 *               constants of a `constexpr` design's Law (detail::ConstantPIDGains)
 */
// The default is void, not detail::PIDGainSet<T, Mode>: a default that names T and Mode crashes clang 23 deducing through the PIController alias.
template<typename T = float, PIDMode Mode = PIDMode::PID, typename Gains = void>
class PIDController {
    static constexpr bool kI = (Mode == PIDMode::PI) || (Mode == PIDMode::PID);                            ///< The tick integrates
    static constexpr bool kD = (Mode == PIDMode::PD) || (Mode == PIDMode::PID);                            ///< The tick differentiates
    using GainSet = damp::conditional_t<damp::is_same_v<Gains, void>, detail::PIDGainSet<T, Mode>, Gains>; ///< The gains' storage
    static constexpr bool kHeld = !GainSet::constant;                                                      ///< The gains are members: set_gains changes them

    template<typename, PIDMode, typename>
    friend class PIDController;

public:
    constexpr PIDController() noexcept = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr PIDController(const design::PIDResult<U>& r) // NOLINT(google-explicit-constructor): assign-from-design
        : PIDController(r.template as<T>()) {}

    /**
     * @brief From a design result. Non-explicit for assign-from-design:
     * `PIController c = design::pid(...).discretize(Ts).as<float>();`
     */
    constexpr PIDController(const design::PIDResult<T>& r) // NOLINT(google-explicit-constructor): assign-from-design
        : u_lim_(r.limits.input), i_(integral_within(r.limits.integral)) {
        if (!r.success && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(r.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
        if (!detail::pid_mode_realizes<Mode>(r.gains.Ki, r.gains.Kd) && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(); // the mode cannot realize a gain: a PI's D term, a P's or a PD's I term
        }
        if (!damp::finite_positive(r.Ts) && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(design::DesignStatus::invalid_sample_time);
        }
        valid_ = r.success && damp::finite_positive(r.Ts) && r.limits.input.valid() && r.limits.integral.valid() && take(r.gains, r.Ts);
        if (!valid_ && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(design::DesignStatus::invalid_parameter); // a non-finite gain, Tf or tracking_time below 0, or a limit that is not a box
        }
    }

    /** @brief The same law, with its live state, in another scalar (or off a Law's constant gains, onto members). */
    template<typename U, typename G>
        requires kHeld
    constexpr explicit PIDController(const PIDController<U, Mode, G>& o)
        : u_lim_(o.u_lim_.template as<T>()),
          u_cmd_(static_cast<T>(o.u_cmd_)),
          h_(static_cast<T>(o.h_)),
          i_(integral_within(o.limits().integral.template as<T>())) {
        valid_ = o.valid_ && take(o.gains().template as<T>(), static_cast<T>(o.sample_time()));
        if constexpr (kI) {
            i_.value = static_cast<T>(o.i_.value);
            i_.step = static_cast<T>(o.i_.step);
        }
        if constexpr (kD) {
            d_.prev_x = static_cast<T>(o.d_.prev_x);
            d_.value = static_cast<T>(o.d_.value);
            d_.first = o.d_.first;
        }
    }

    /** @brief False when built from a design it refuses: step() returns 0. */
    [[nodiscard]] constexpr bool valid() const { return valid_; }

    /** @brief One fixed-rate tick at sample_time(): u from the reference @p r and the measurement @p y. */
    [[nodiscard]] constexpr T step(T r, T y) {
        if (!valid_) {
            return T{0};
        }
        return tick(r, y, ki_per(k_.Ts, true), fixed_derivative(), k_.Ts);
    }

    /**
     * @brief One tick over the measured period @p dt [s]
     *
     * Kᵢ·dt and the derivative's coefficients for this dt, one reciprocal 1/(dt + T_f) shared by
     * both (none without a derivative). A @p dt that is not finite and positive holds the state
     * and returns the last command; a back_calculate after it is a no-op. A back_calculate after
     * a measured tick charges that tick's dt.
     */
    [[nodiscard]] constexpr T step(T r, T y, T dt) {
        if (!valid_) {
            return T{0};
        }
        if (!damp::finite_positive(dt)) {
            h_ = T{0};
            return u_cmd_;
        }
        return tick(r, y, ki_per(dt, false), derivative_at(dt), dt);
    }

    /**
     * @brief In place of step(@p r, @p y) while another source drives the plant: the memory that starts the next command from @p u, the input applied
     *
     * The derivative memory moves as the tick's, then the integral is preloaded so that
     * step(@p r, @p y) at this same r and y returns @p u clamped into the input limits (the
     * setpoint weight b included), then clamped into the integral limits, which leaves the command
     * short of u by what that clamp took. @p u is recorded as the command applied, so a
     * back_calculate of it finds no excess. A P or a PD has nothing to preload: its derivative
     * memory is kept current, so the first step has no kick.
     */
    constexpr void track(T u, T r, T y) {
        if (!valid_) {
            return;
        }
        follow(u, r, y, ki_per(k_.Ts, true), fixed_derivative(), k_.Ts);
    }

    /**
     * @brief In place of step(@p r, @p y, @p dt): track(@p u, @p r, @p y) for the next measured tick
     *
     * The preload is exact for a next tick of the same @p dt. A @p dt that is not finite and
     * positive holds the state and records @p u.
     */
    constexpr void track(T u, T r, T y, T dt) {
        if (!valid_) {
            return;
        }
        if (!damp::finite_positive(dt)) {
            h_ = T{0};
            u_cmd_ = u;
            return;
        }
        follow(u, r, y, ki_per(dt, false), derivative_at(dt), dt);
    }

    /**
     * @brief Anti-windup after a downstream clamp: the integral holds against the excess
     *        u_sat − (u + u_ff), or with a tracking time unwinds by min(h/T_t, 1) of it
     *
     * Call after the stage that clamps the summed command, in the same tick as step(): @p u_sat
     * is what the plant received, @p u_ff whatever was added to this law's command before the
     * clamp, u the command step() returned (@ref last_command) and h that tick's period, T_s or
     * the measured dt. With T_t = 0 the integral gives back what the last tick added, up to the
     * excess. A P or a PD has nothing to wind up: a no-op.
     */
    constexpr void back_calculate(T u_sat, T u_ff = T{0}) {
        if constexpr (kI) {
            if (valid_) {
                wind_back(u_sat, u_cmd_ + u_ff);
            }
        } else {
            static_cast<void>(u_sat);
            static_cast<void>(u_ff);
        }
    }

    /** @brief The command the last step() returned, or the input the last track() recorded. */
    [[nodiscard]] constexpr T last_command() const { return u_cmd_; }

    /** @brief Clear the integral, the derivative and the last command; keeps the gains and the limits. */
    constexpr void reset() {
        if constexpr (kI) {
            i_.value = T{0};
            i_.step = T{0};
        }
        if constexpr (kD) {
            d_.prev_x = T{0};
            d_.value = T{0};
            d_.first = true;
        }
        u_cmd_ = T{0};
        h_ = T{0};
    }

    /**
     * @brief The continuous gains in force: the design's, or the last set_gains
     *
     * A term the mode does not compute reads at its default: Kᵢ = 0 and T_t = 0 on a P or a PD,
     * K_d = 0, T_f = 0 and c = 1 on a P or a PI.
     */
    [[nodiscard]] constexpr design::PIDGains<T> gains() const {
        design::PIDGains<T> g{};
        g.Kp = k_.Kp;
        g.b = k_.b;
        if constexpr (kI) {
            g.Ki = k_.i.Ki;
            g.tracking_time = k_.i.tt;
        }
        if constexpr (kD) {
            g.Kd = k_.d.Kd;
            g.Tf = k_.d.Tf;
            g.c = k_.d.c;
        }
        return g;
    }

    /**
     * @brief Change the gains at run time; the next tick uses them
     *
     * Recomputes the fixed-rate coefficients at sample_time() and keeps the integral I: it is the
     * I term, in u's units, so a new Kᵢ changes its rate and not its value, and the command does
     * not jump. Refused (false, nothing changed) for a failed @p g, a non-finite gain, a Tf or a
     * tracking_time below 0 or not finite, a gain the mode does not compute (K_d ≠ 0 on a PI,
     * Kᵢ ≠ 0 on a P or a PD), or a coefficient that overflows at sample_time().
     *
     * @code
     * auto g = pid.gains();
     * g.Kp = 2.5f;
     * if (!pid.set_gains(g)) { } // refused: nothing changed
     * @endcode
     */
    [[nodiscard]] constexpr bool set_gains(const design::PIDGains<T>& g)
        requires kHeld
    {
        return g.success && take(g, k_.Ts);
    }

    /** @brief The limits in force: the design's, or the last set_limits (integral unbounded on a P or a PD). */
    [[nodiscard]] constexpr design::PIDLimits<T> limits() const {
        design::PIDLimits<T> lim{};
        lim.input = u_lim_;
        if constexpr (kI) {
            lim.integral = i_.lim;
        }
        return lim;
    }

    /**
     * @brief Change the limits at run time; the next tick clamps into them
     *
     * The integral is clamped into the new integral limits at once. Refused (false, nothing
     * changed) when an end is NaN or ±Inf, or a lower end is above its upper one. A P or a PD
     * has no integral: its integral limits are checked and not kept.
     */
    [[nodiscard]] constexpr bool set_limits(const design::PIDLimits<T>& lim) {
        if (!lim.input.valid() || !lim.integral.valid()) {
            return false;
        }
        u_lim_ = Bounds<1, T>(lim.input);
        if constexpr (kI) {
            i_.lim = Bounds<1, T>(lim.integral);
            const T held = i_.value;
            i_.value = i_.lim.saturate(i_.value);
            if (i_.value != held) {
                i_.step = T{0}; // the clamp took it: nothing of the last step left to give back
            }
        }
        return true;
    }

    /** @brief The integral term I (u's units; 0 on a P or a PD). */
    [[nodiscard]] constexpr T integral() const {
        if constexpr (kI) {
            return i_.value;
        } else {
            return T{0};
        }
    }

    /** @brief I ← @p value clamped into the integral limits: a preload, or a leak an outer loop applies. */
    constexpr void set_integral(T value)
        requires kI
    {
        i_.value = i_.lim.saturate(value);
        i_.step = T{0};
    }

    /** @brief The filtered derivative term D (u's units; 0 on a P or a PI). */
    [[nodiscard]] constexpr T derivative() const {
        if constexpr (kD) {
            return d_.value;
        } else {
            return T{0};
        }
    }

    /** @brief [s] The fixed-rate tick's sample period T_s (0 for a law with no fixed period). */
    [[nodiscard]] constexpr T sample_time() const { return k_.Ts; }

private:
    /**
     * @brief Store @p g at the period @p Ts when it is in the domain, the mode computes it and no
     *        coefficient overflows; on a constexpr design's Law, whose gains are the type's, only check it
     */
    [[nodiscard]] constexpr bool take(const design::PIDGains<T>& g, T Ts) {
        const auto set = detail::pid_gain_set<Mode>(g, Ts);
        if (!set) {
            return false;
        }
        if constexpr (kHeld) {
            k_ = *set;
        }
        return true;
    }

    /** @brief The integral term at rest within the limits @p lim (nothing on a P or a PD): built, not assigned, so a constant initializer can copy the design's limits. */
    [[nodiscard]] static constexpr auto integral_within(const Bounds<1, T>& lim) {
        if constexpr (kI) {
            return detail::PIDIntegral<T>{.lim = lim};
        } else {
            static_cast<void>(lim);
            return detail::maybe<false, detail::PIDIntegral<T>, detail::OffSlot::pid_integral>{};
        }
    }

    /** @brief The filtered derivative's coefficients over one period. */
    struct Derivative {
        T a{T{0}}; ///< T_f/(dt + T_f): the pole
        T b{T{0}}; ///< K_d/(dt + T_f): the gain on Δ(c·r − y)
    };

    /** @brief The integral's increment per unit error over @p h: the stored Kᵢ·T_s at the fixed rate, Kᵢ·h at a measured one, 0 without an integral. */
    [[nodiscard]] constexpr T ki_per(T h, bool fixed) const {
        if constexpr (kI) {
            return fixed ? k_.i.ki_ts : k_.i.Ki * h;
        } else {
            static_cast<void>(h);
            static_cast<void>(fixed);
            return T{0};
        }
    }

    /** @brief The derivative's pole and gain at the fixed rate, computed when the gains were set. */
    [[nodiscard]] constexpr Derivative fixed_derivative() const {
        if constexpr (kD) {
            return Derivative{k_.d.a, k_.d.b};
        } else {
            return Derivative{};
        }
    }

    /** @brief The derivative's pole and gain over @p dt: one reciprocal shared by both, none without a derivative. */
    [[nodiscard]] constexpr Derivative derivative_at(T dt) const {
        if constexpr (kD) {
            if (k_.d.Kd != T{0}) {
                const T inv = T{1} / (dt + k_.d.Tf);
                return Derivative{k_.d.Tf * inv, k_.d.Kd * inv};
            }
        }
        static_cast<void>(dt);
        return Derivative{};
    }

    /** @brief The derivative term over one period: D ← a·D + b·Δ(c·r − y), Δ zero on the first tick (no derivative kick). */
    constexpr void differentiate(T r, T y, const Derivative& k)
        requires kD
    {
        const T x = (k_.d.c * r) - y;
        if (d_.first) {
            d_.prev_x = x;
            d_.first = false;
        }
        d_.value = (k.a * d_.value) + (k.b * (x - d_.prev_x));
        d_.prev_x = x;
    }

    /** @brief One tick over a period @p h: integral increment @p ki_h per unit error, derivative coefficients @p k. */
    [[nodiscard]] constexpr T tick(T r, T y, T ki_h, const Derivative& k, T h) {
        T u_unsat = k_.Kp * ((k_.b * r) - y);
        if constexpr (kD) {
            differentiate(r, y, k);
            u_unsat = u_unsat + d_.value;
        }
        if constexpr (kI) {
            const T held = i_.value;
            i_.value = i_.lim.saturate(i_.value + (ki_h * (r - y)));
            i_.step = i_.value - held;
            u_unsat = u_unsat + i_.value;
        }
        static_cast<void>(ki_h);
        static_cast<void>(k);
        const T u = u_lim_.saturate(u_unsat);
        h_ = h;
        u_cmd_ = u;
        if constexpr (kI) {
            wind_back(u, u_unsat);
        }
        return u;
    }

    /** @brief track's memory: the tick's derivative, then I for a next tick at (r, y) of the same period that returns u clamped. */
    constexpr void follow(T u, T r, T y, T ki_h, const Derivative& k, T h) {
        T next = k_.Kp * ((k_.b * r) - y); // the next tick's P + D at this (r, y)
        if constexpr (kD) {
            differentiate(r, y, k);
            next = next + (k.a * d_.value); // Δ(c·r − y) = 0 at the next tick: D only decays
        }
        if constexpr (kI) {
            i_.value = i_.lim.saturate(u_lim_.saturate(u) - next - (ki_h * (r - y)));
            i_.step = T{0};
        }
        static_cast<void>(ki_h);
        static_cast<void>(k);
        h_ = h;
        u_cmd_ = u;
    }

    /** @brief Hold I against @p u_sat − @p u_unsat (conditional integration), or with T_t > 0 unwind min(h/T_t, 1) of it. */
    constexpr void wind_back(T u_sat, T u_unsat)
        requires kI
    {
        if (k_.i.Ki == T{0} || h_ == T{0} || u_sat == u_unsat) {
            return;
        }
        if (k_.i.inv_tt == T{0}) {
            detail::hold_against_clamp(i_.value, i_.step, u_sat - u_unsat);
            return;
        }
        // The whole excess at most: past it an overshoot of the unsaturated value, past twice divergence.
        i_.value = i_.lim.saturate(i_.value + (damp::min(h_ * k_.i.inv_tt, T{1}) * (u_sat - u_unsat)));
    }

    using Integral = detail::maybe<kI, detail::PIDIntegral<T>, detail::OffSlot::pid_integral>;
    using Deriv = detail::maybe<kD, detail::PIDDerivative<T>, detail::OffSlot::pid_derivative>;

    [[no_unique_address]] GainSet  k_{};         ///< The gains, T_s and the coefficients at it: members, or no storage on a constexpr design's Law
    Bounds<1, T>                   u_lim_{};     ///< The command's limits
    T                              u_cmd_{T{0}}; ///< Last command step() returned, or the input track() recorded
    T                              h_{T{0}};     ///< [s] The last tick's period, which back_calculate charges (0: none)
    [[no_unique_address]] Integral i_{};         ///< The integral term: no storage on a P or a PD
    [[no_unique_address]] Deriv    d_{};         ///< The derivative term: no storage on a P or a PI
    bool                           valid_{true}; ///< False when built from a design it refuses: step() returns 0
};

/**
 * @brief Deduce the runtime from its design, keeping its scalar: PIDController pid{res} is a PID;
 * through an alias the mode comes from the alias: PIController pi{res};
 */
template<typename T, PIDMode Mode = PIDMode::PID>
PIDController(const design::PIDResult<T>&) -> PIDController<T, Mode>;

template<typename T = float>
using PController = PIDController<T, PIDMode::P>;

template<typename T = float>
using PIController = PIDController<T, PIDMode::PI>;

template<typename T = float>
using PDController = PIDController<T, PIDMode::PD>;

namespace detail {

/** @brief Whether @p R is a PID design result. */
template<typename R>
inline constexpr bool is_pid_result = false;

/** @brief design::PIDGains::discretize's result. */
template<typename U>
inline constexpr bool is_pid_result<design::PIDResult<U>> = true;

/** @brief A PID design (design::PIDGains::discretize). */
template<const auto& Res>
concept pid_design = is_pid_result<design_t<Res>>;

/** @brief The narrowest mode that computes every term of @p g: no integral at Kᵢ = 0, no derivative at K_d = 0. */
template<typename U>
[[nodiscard]] constexpr PIDMode narrowest_pid_mode(const design::PIDGains<U>& g) {
    const bool has_i = g.Ki != U{0};
    const bool has_d = g.Kd != U{0};
    if (has_i) {
        return has_d ? PIDMode::PID : PIDMode::PI;
    }
    return has_d ? PIDMode::PD : PIDMode::P;
}

} // namespace detail

/**
 * @brief Law<kRes, T> of a PID design initialized at start-up: the PIDController it runs as. A runtime of the same shape assigned replaces it.
 */
template<const auto& Res, typename T>
    requires detail::pid_design<Res>
struct Law<Res, T> : PIDController<T> {
    using Runtime = PIDController<T>; ///< The runtime the design runs as

    /** @brief The design's law. */
    constexpr Law() : Runtime{Res} {}

    using Runtime::operator=;
};

namespace detail {

/**
 * @brief A `constexpr` PID design's gains as constants of the type: what a PIDController reads in
 *        place of its gain members, in the runtime's scalar, at the mode the gains need
 *
 * Each member is static, so the Law stores none of them: the gains sit in flash and fold into
 * the tick (a zero or unit gain is no code). A result whose gains its mode refuses has no set,
 * which the Law's static_assert reports.
 */
template<const auto& Res, typename T, PIDMode Mode>
struct ConstantPIDGains {
    static constexpr bool constant = true; ///< The gains are the type's: set_gains is absent

    static constexpr auto set = [] {
        const auto s = pid_gain_set<Mode>(Res.gains.template as<T>(), static_cast<T>(Res.Ts));
        return s ? *s : PIDGainSet<T, Mode>{};
    }(); ///< The gains at the design's T_s

    static constexpr bool ok = pid_gain_set<Mode>(Res.gains.template as<T>(), static_cast<T>(Res.Ts)).has_value(); ///< The gains are in the mode's domain at T_s

    static constexpr T    Kp = set.Kp; ///< Proportional gain
    static constexpr T    b = set.b;   ///< Proportional setpoint weight
    static constexpr T    Ts = set.Ts; ///< [s] The fixed-rate tick's period
    static constexpr auto i = set.i;   ///< The integral's gains (empty on a P or a PD)
    static constexpr auto d = set.d;   ///< The derivative's gains (empty on a P or a PI)
};

} // namespace detail

/**
 * @brief Law<kRes, T> of a `constexpr` PID design: the PIDController of the narrowest mode its
 *        gains need, with those gains constants of the type
 *
 * K_d = 0 is a PI, Kᵢ = 0 a PD, both a P, so the tick computes no term the design has not and
 * stores none of its memory. The gains, T_s and the coefficients at it are the type's
 * (detail::ConstantPIDGains): in flash, folded into the tick, no RAM. The limits stay members,
 * settable at run time; there is no set_gains. A law retuned at run time is a PIDController,
 * `Law<kRes, T>::Runtime law{kRes.as<T>()};`. A failed design, one with no positive T_s or with a
 * gain its mode refuses does not compile wherever the Law is declared.
 */
template<const auto& Res, typename T>
    requires(detail::pid_design<Res> && detail::constant_design<Res>)
struct Law<Res, T> : PIDController<T, detail::narrowest_pid_mode(Res.gains), detail::ConstantPIDGains<Res, T, detail::narrowest_pid_mode(Res.gains)>> {
    static constexpr PIDMode kMode = detail::narrowest_pid_mode(Res.gains); ///< The terms the tick computes
    using Runtime = PIDController<T, kMode>;                                ///< The runtime with members: the type of a member a redesign assigns
    using Base = PIDController<T, kMode, detail::ConstantPIDGains<Res, T, kMode>>;

    static_assert(detail::constant_design_runs(Res.success, Res.status, Res.Ts));
    static_assert(detail::ConstantPIDGains<Res, T, kMode>::ok, "Law<kRes, T>: a gain of the PID design is not finite, a Tf or tracking_time is below 0, or a coefficient overflows at its Ts");

    /** @brief The design's law: its gains the type's, its limits members. */
    constexpr Law() : Base{Res} {}

    /** @brief There is nothing to replace: the gains are in the type. */
    template<typename Other>
        requires(!damp::is_same_v<damp::remove_cvref_t<Other>, Law>)
    constexpr Law& operator=(const Other& other) {
        static_assert(sizeof(Other) == 0, "This Law is a constexpr PID design's: its gains are constants of the type, so there is nothing to replace. A member a redesign assigns, or a law retuned with set_gains, is Law<kRes, T>::Runtime law{kRes.as<T>()};");
        static_cast<void>(other);
        return *this;
    }
};

} // namespace damp

#pragma pop_macro("PI")
