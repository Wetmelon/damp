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
#include "damp/design/status.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/toolbox/bounds.hpp"

namespace damp {
namespace detail {

/**
 * @brief True when Kbc and Ki have opposite signs
 *
 * Back-calculation adds Ki·(u_sat − u_unsat)/Kbc to the integral term, which enters u with
 * unit gain, so it opposes windup only when the tracking time T_t = Kbc/Ki is positive.
 */
template<typename T>
[[nodiscard]] constexpr bool pid_kbc_winds_up(T Ki, T Kbc) {
    return ((Ki > T{0}) && (Kbc < T{0})) || ((Ki < T{0}) && (Kbc > T{0}));
}

/**
 * @brief True when a PID gain set has no non-finite entry and a non-negative derivative filter
 *
 * A Tf < 0 puts the filtered-derivative pole Tf/(Ts + Tf) outside the unit circle for
 * Tf < −Ts/2, and divides by zero at Tf = −Ts.
 */
template<typename T>
[[nodiscard]] constexpr bool pid_gains_in_domain(T Kp, T Ki, T Kd, T Tf, T Kbc, T b, T c) {
    return damp::isfinite(Kp) && damp::isfinite(Ki) && damp::isfinite(Kd) && damp::finite_non_negative(Tf)
        && damp::isfinite(Kbc) && damp::isfinite(b) && damp::isfinite(c);
}

} // namespace detail

namespace design {

/**
 * @struct DiscretePIDGains
 * @brief Fixed-rate discrete PID coefficients (canonical deploy form)
 *
 * Produced by PIDGains::discretize. Sample time is baked into the
 * integral gain and derivative filter coeffs so the runtime controller is a
 * pure `control(r, y)` tick with no per-call `Ts`.
 *
 * Discretization (backward Euler integrator, filtered derivative):
 * - Integral term: @f$K_{i,d} = K_i T_s@f$, @f$I \leftarrow I + K_{i,d}\,e@f$,
 *   @f$u_I = I@f$ (same units as @f$u@f$)
 * - @f$D = a_d D + b_d\,\Delta x@f$ with
 *   @f$a_d = T_f/(T_s+T_f)@f$, @f$b_d = K_d/(T_s+T_f)@f$
 *   (@f$T_f=0@f$ ⇒ raw @f$K_d/T_s@f$ difference)
 * - Integrator limits clamp @f$I@f$ (same units as @f$u@f$), copied from the
 *   continuous design result
 * - Back-calculation: @f$\Delta I = K_{i,d}(u-u_{\mathrm{unsat}})/K_{bc}@f$
 *
 * @see PIDGains::discretize
 * @see Åström & Hägglund, "Advanced PID Control" (2006), Sec. 3.3
 */
template<typename T = double>
struct [[nodiscard]] DiscretePIDGains {
    T            Kp{};
    T            Ki{};           ///< Discrete integral gain Kᵢ T_s (I += Ki·e)
    T            Kd{};           ///< Continuous K_d (documentation / re-discretize); runtime uses @ref d_a / @ref d_b
    T            d_a{};          ///< Filtered-derivative pole: D ← d_a D + d_b Δx
    T            d_b{};          ///< Filtered-derivative input gain
    Bounds<1, T> u_lim{};        ///< Output limits [lower, upper] (default unbounded)
    Bounds<1, T> integral_lim{}; ///< Clamp on integral term I (same units as u; default unbounded)
    T            Kbc = T{0};     ///< Back-calculation coefficient, same sign as Ki; 0 = output clamp only
    T            b = T{1};       ///< Proportional setpoint weight
    T            c = T{1};       ///< Derivative setpoint weight
    T            Tf = T{0};      ///< Continuous derivative filter time constant (documentation)
    T            Ts = T{0};      ///< Sample period used for discretization [s]

    bool         success{true};            ///< status == DesignStatus::ok (a hand-built set is a valid design)
    DesignStatus status{DesignStatus::ok}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return DiscretePIDGains<U>{
            static_cast<U>(Kp),
            static_cast<U>(Ki),
            static_cast<U>(Kd),
            static_cast<U>(d_a),
            static_cast<U>(d_b),
            u_lim.template as<U>(),
            integral_lim.template as<U>(),
            static_cast<U>(Kbc),
            static_cast<U>(b),
            static_cast<U>(c),
            static_cast<U>(Tf),
            static_cast<U>(Ts),
            success,
            status,
        };
    }
};

/**
 * @struct PIDGains
 * @brief 2-DOF continuous-time PID gains, limits and setpoint weights
 *
 * An aggregate: every field has a safe default (unbounded limits, no
 * back-calculation, b = c = 1, unfiltered derivative), so name only the knobs
 * you set — `design::PIDGains{.Kp = 4.0, .Ki = 8.0, .u_lim = {.lower = -12.0, .upper = 12.0}}`.
 * @ref design::pid covers the plain (Kp, Ki, Kd) case. Deploy via
 * discretize into a @ref DiscretePIDGains for the fixed-rate
 * PIDController, or hand the continuous gains to @ref ContinuousPID
 * when the sample period is measured each tick.
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
    T            Kp{};
    T            Ki{};
    T            Kd{};
    Bounds<1, T> u_lim{};        ///< Output limits [lower, upper] (default unbounded)
    Bounds<1, T> integral_lim{}; ///< Clamp on integral term I (same units as u; default unbounded)
    T            Kbc = T{0};     ///< Back-calculation coefficient: ΔI += Ki·(Ts/Kbc)(u−u_unsat).
                                 ///< Tracking time T_t = Kbc/Ki (parallel form); Kbc = Kp ⇒ T_t = T_i.
                                 ///< Same sign as Ki (T_t > 0). Larger |Kbc| → slower unwind; 0 = output clamp only.
    T            b = T{1};       ///< Proportional setpoint weight (0=I-PD, 1=standard PID)
    T            c = T{1};       ///< Derivative setpoint weight   (0=PI-D, 1=standard PID)
    T            Tf = T{0};      ///< Derivative filter time constant (0 = unfiltered)

    bool         success{true};            ///< status == DesignStatus::ok (a hand-built set is a valid design)
    DesignStatus status{DesignStatus::ok}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return PIDGains<U>{
            static_cast<U>(Kp),
            static_cast<U>(Ki),
            static_cast<U>(Kd),
            u_lim.template as<U>(),
            integral_lim.template as<U>(),
            static_cast<U>(Kbc),
            static_cast<U>(b),
            static_cast<U>(c),
            static_cast<U>(Tf),
            success,
            status,
        };
    }

    /**
     * @brief Bake sample period into discrete coefficients (canonical deploy path)
     *
     * @code
     * constinit PIController controller =
     *     design::pid(4.0, 8.0, 0.0).discretize(0.01).as<float>();
     * @endcode
     *
     * A failed design stays failed (its status is carried over), a @p ts that is not
     * positive and finite is refused with DesignStatus::invalid_sample_time, and a
     * non-finite gain, a Tf < 0, or a Kbc whose sign differs from Ki's (T_t = Kbc/Ki < 0,
     * back-calculation would wind the integrator up) is refused with
     * DesignStatus::invalid_parameter.
     *
     * @param ts Sample period [s]
     * @return DiscretePIDGains for PIDController
     */
    [[nodiscard]] constexpr DiscretePIDGains<T> discretize(T ts) const {
        if (!success) {
            return ::damp::detail::refused<DiscretePIDGains<T>>(status);
        }
        if (!(ts > T{0}) || !damp::isfinite(ts)) {
            return ::damp::detail::refused<DiscretePIDGains<T>>(DesignStatus::invalid_sample_time);
        }
        if (!::damp::detail::pid_gains_in_domain(Kp, Ki, Kd, Tf, Kbc, b, c) || ::damp::detail::pid_kbc_winds_up(Ki, Kbc)) {
            return ::damp::detail::refused<DiscretePIDGains<T>>(DesignStatus::invalid_parameter);
        }
        DiscretePIDGains<T> d{};
        d.Kp = Kp;
        d.Kd = Kd;
        d.u_lim = u_lim;
        d.Kbc = Kbc;
        d.b = b;
        d.c = c;
        d.Tf = Tf;
        d.Ts = ts;
        d.integral_lim = integral_lim;
        d.Ki = Ki * ts;
        const T den = ts + Tf;
        d.d_a = Tf / den;
        d.d_b = Kd / den;
        return d;
    }

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
        if ((Kd != T{0}) && !(Tf > T{0})) {
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
 * @brief Continuous PID gains Kₚ, Kᵢ, K_d with every other knob at its default
 *
 * Short form for the common case: unbounded output and integrator, no
 * back-calculation, b = c = 1, unfiltered derivative. For limits, setpoint
 * weights, or a derivative filter, build @ref PIDGains directly with
 * designated initializers. Deploy with PIDGains::discretize for fixed-rate
 * PIDController, or construct @ref ContinuousPID for variable-rate ticks.
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
 *     design::PIDGains{.Kp = 4.0, .Ki = 8.0, .u_lim = {.lower = -12.0, .upper = 12.0}}
 *         .discretize(0.01)
 *         .as<float>();
 * @endcode
 *
 * @param Kp Proportional gain
 * @param Ki Integral gain (continuous)
 * @param Kd Derivative gain
 *
 * @return PIDGains{Kp, Ki, Kd} with default limits, weights and filter
 */
template<typename T = double>
[[nodiscard]] constexpr PIDGains<T> pid(T Kp, T Ki, T Kd) {
    return PIDGains<T>{.Kp = Kp, .Ki = Ki, .Kd = Kd};
}

} // namespace design

/**
 * @brief Compile-time selection of the PID control-law structure.
 *
 * Picks which specialization of PIDController is instantiated, so each mode
 * stores and computes only its own terms. Distinct from @ref PIDRuntimeMode,
 * which selects runtime behavior (Auto vs Tracking) at each tick.
 */
enum class PIDMode : uint8_t {
    P,
    PI,
    PD,
    PID,
};

/**
 * @brief Runtime operating mode for PIDController / @ref ContinuousPID.
 *
 * - `Auto`:     standard control law with back-calculation anti-windup. The
 *               controller follows the reference.
 * - `Tracking`: the controller's output is `u_track` clamped into the output limits `u`
 *               and on every tick the integrator (and derivative state) is
 *               pre-loaded so that the next `enable()` call resumes in `Auto`
 *               *without a bump* in command. This is the standard ISA-PID
 *               pattern for operator manual mode, leader/follower handoff,
 *               and multi-controller bumpless transfer.
 *
 * Semantically `u_track` is the same quantity as the `u_sat` argument of
 * back_calculate — both represent the value actually being applied to the
 * plant. Use `back_calculate` when the controller stays in `Auto` but its
 * output was clamped downstream; use `disable(u_track)` when the controller is
 * being temporarily replaced by another command source.
 */
enum class PIDRuntimeMode : uint8_t {
    Auto,
    Tracking,
};

namespace detail {

// Integral + output stage shared by PIDController and ContinuousPID. Free
// functions over the controller's public fields, so neither type needs a base.

/** @brief I ← v clamped into the integrator limits i */
template<typename C, typename T>
constexpr void pid_set_integral(C& c, T v) {
    c.integral = c.integral_lim.saturate(v);
}

/** @brief Back-calculation ΔI = gain·(u_sat − u_unsat)/Kbc. No-op when Kbc or gain is zero. */
template<typename C, typename T>
constexpr void pid_wind_back(C& c, T gain, T u_sat, T u_unsat) {
    if (c.Kbc == T{0} || gain == T{0} || u_sat == u_unsat) {
        return;
    }
    pid_set_integral(c, c.integral + (gain * ((u_sat - u_unsat) / c.Kbc)));
}

/**
 * @brief One tick: Tracking holds u_track, Auto integrates and saturates.
 *
 * @param dI       Integral increment this tick
 * @param pd       Proportional + derivative terms this tick
 * @param pd_next  P + D the next tick would produce at the same (r, y)
 * @param gain     Back-calculation gain (Ki, or Ki·Ts for a measured period)
 *
 * Tracking preloads I = u_track − pd_next − dI, so the next Auto tick at the
 * same (r, y) emits exactly u_track.
 */
template<typename C, typename T>
[[nodiscard]] constexpr T pid_step(C& c, T dI, T pd, T pd_next, T gain) {
    if (c.runtime_mode == PIDRuntimeMode::Tracking) {
        pid_set_integral(c, c.u_track - pd_next - dI);
        c.u_cmd = c.u_lim.saturate(c.u_track);
        return c.u_cmd;
    }
    pid_set_integral(c, c.integral + dI);
    const T u_unsat = pd + c.integral;
    const T u = c.u_lim.saturate(u_unsat);
    pid_wind_back(c, gain, u, u_unsat);
    c.u_cmd = u;
    return u;
}

/** @brief Δ(c·r − y) since the last tick; zero on the first (no derivative kick). */
template<typename C, typename T>
[[nodiscard]] constexpr T pid_delta(C& c, T x) {
    if (c.first) {
        c.prev_x = x;
        c.first = false;
    }
    const T dx = x - c.prev_x;
    c.prev_x = x;
    return dx;
}

} // namespace detail

template<typename T = float, PIDMode Mode = PIDMode::PID>
struct PIDController;

/**
 * @ingroup discrete_controllers
 * @brief Fixed-rate discrete 2-DOF PID (canonical runtime)
 *
 * Sample time is baked into the gains at design time via
 * @ref design::PIDGains::discretize. One tick is `control(r, y)` — no `Ts`
 * argument. Matches the library-wide pre-discretized controller model
 * (LQR, Kalman, PR, lead-lag).
 *
 * @ref PIDMode picks the specialization (P, PI, PD, PID); each stores and
 * computes only its own terms. A design result whose gains the mode cannot
 * realize, or whose Kbc has the opposite sign to Ki, is refused at
 * construction: ill-formed under constant evaluation (`constinit` deploy), a
 * trap at runtime.
 *
 * Control law:
 * @f[
 *   I \leftarrow I + K_{i,d}\,e,\quad
 *   D \leftarrow a_d D + b_d\,\Delta(c r - y),\quad
 *   u = K_p(b r - y) + I + D
 * @f]
 *
 * @code
 * constinit auto c = PIDController{
 *     design::pid(2.0, 5.0, 0.1).discretize(0.01).as<float>()};
 * float u = c.control(r, y);
 * @endcode
 *
 * For continuous gains with a measured period each tick, use @ref ContinuousPID.
 *
 * Angle input: pass r′ = y + wrap_pi(r − y), the reference nearest the measurement,
 * `pid.control(y + wrap_pi(r - y), y)`, so P and I see the short-way error. With c ≠ 1
 * (derivative on measurement) D also reads y, so unwrap y first with @ref PhaseUnwrapper.
 * Keep b = 1 on an angle axis: with b ≠ 1, P carries −Kₚ(1 − b)·y, which depends on the
 * absolute angle and the turn count rather than on the error.
 *
 * MATLAB® equivalent: design with `pid(...)` then `c2d`; runtime is the discrete law.
 * A profiled PID takes its setpoint from a motion profile each tick (the motion pack's
 * trapezoidal profile, not yet public).
 * @see Åström & Hägglund, "Advanced PID Control" (2006), Sec. 3.3–4.4
 *
 * @note Gains from a continuous-time designer (ziegler_nichols, cohen_coon,
 *       simc, lambda_tuning, pidtune) are Ts-agnostic until discretize(Ts).
 *       Gains from the sampled design::pid_pole_placement(plant, poles, Ts) already
 *       embed that Ts in continuous form — discretize with the same Ts.
 */
template<typename T>
struct PIDController<T, PIDMode::PID> {
    T            Kp{};
    T            Ki{};  ///< Discrete integral gain Kᵢ T_s (I += Ki·e)
    T            Kd{};  ///< Continuous Kd retained for inspection; tick uses d_a / d_b
    T            d_a{}; ///< D ← d_a D + d_b Δx
    T            d_b{};
    Bounds<1, T> u_lim{};        ///< Output limits [lower, upper] (default unbounded)
    Bounds<1, T> integral_lim{}; ///< Clamp on integral term I (same units as u; default unbounded)
    T            Kbc = T{0};     ///< Back-calculation: ΔI += Ki·(u−u_unsat)/Kbc; 0 = output clamp only
    T            b = T{1};       ///< Proportional setpoint weight
    T            c = T{1};       ///< Derivative setpoint weight
    T            Ts = T{0};      ///< Sample period from discretize (documentation)

    T    integral = T{0}; ///< Integral term I (same units as u)
    T    prev_x = T{0};   ///< Previous c·r − y
    T    deriv = T{0};    ///< Filtered derivative term D
    bool first = true;    ///< Next tick seeds prev_x (no derivative kick)

    PIDRuntimeMode runtime_mode{PIDRuntimeMode::Auto};
    T              u_track{T{0}};
    T              u_cmd{T{0}}; ///< Last command control() returned: the u in back_calculate's u_sat − (u + u_ff)

    constexpr PIDController() noexcept = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr PIDController(const design::DiscretePIDGains<U>& r) // NOLINT(google-explicit-constructor): assign-from-design
        : PIDController(r.template as<T>()) {}

    /**
     * @brief From discrete design result. Non-explicit for assign-from-design:
     * `PIDController c = design::pid(...).discretize(Ts).as<float>();`
     */
    constexpr PIDController(const design::DiscretePIDGains<T>& r) // NOLINT(google-explicit-constructor): assign-from-design
        : Kp(r.Kp),
          Ki(r.Ki),
          Kd(r.Kd),
          d_a(r.d_a),
          d_b(r.d_b),
          u_lim(r.u_lim),
          integral_lim(r.integral_lim),
          Kbc(r.Kbc),
          b(r.b),
          c(r.c),
          Ts(r.Ts) {
        if (!r.success && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(r.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
        if (detail::pid_kbc_winds_up(r.Ki, r.Kbc)) {
            detail::refuse_failed_design(design::DesignStatus::invalid_parameter); // Kbc against Ki's sign winds the integrator up
        }
    }

    template<typename U>
    constexpr explicit PIDController(const PIDController<U, PIDMode::PID>& o)
        : Kp(static_cast<T>(o.Kp)),
          Ki(static_cast<T>(o.Ki)),
          Kd(static_cast<T>(o.Kd)),
          d_a(static_cast<T>(o.d_a)),
          d_b(static_cast<T>(o.d_b)),
          u_lim(o.u_lim.template as<T>()),
          integral_lim(o.integral_lim.template as<T>()),
          Kbc(static_cast<T>(o.Kbc)),
          b(static_cast<T>(o.b)),
          c(static_cast<T>(o.c)),
          Ts(static_cast<T>(o.Ts)),
          integral(static_cast<T>(o.integral)),
          prev_x(static_cast<T>(o.prev_x)),
          deriv(static_cast<T>(o.deriv)),
          first(o.first),
          runtime_mode(o.runtime_mode),
          u_track(static_cast<T>(o.u_track)),
          u_cmd(static_cast<T>(o.u_cmd)) {}

    /**
     * @brief One fixed-rate discrete PID tick
     *
     * In Auto mode runs the 2-DOF law with back-calculation anti-windup.
     * In Tracking mode the output is `u_track` clamped into the output limits `u` and the
     * integrator is pre-loaded for bumpless re-enable.
     */
    [[nodiscard]] constexpr T control(T r, T y) {
        deriv = (d_a * deriv) + (d_b * detail::pid_delta(*this, (c * r) - y));
        const T p = Kp * ((b * r) - y);
        return detail::pid_step(*this, Ki * (r - y), p + deriv, p + (d_a * deriv), Ki);
    }

    /**
     * @brief Hold in Tracking mode, preloading the integrator to emit @p u_track bumplessly
     * @param u_track_ Output value to follow
     * @param y        Current measurement (preload assumes r = y, zero error)
     */
    constexpr void track(T u_track_, T y) {
        disable(u_track_);
        (void)control(y, y);
    }

    /**
     * @brief Anti-windup after a downstream clamp: ΔI = Ki·(u_sat − (u + u_ff))/Kbc
     *
     * Call after the stage that clamps the summed command. @p u_sat is what the plant
     * received, @p u_ff whatever was added to this law's command before the clamp; u is
     * the command control() returned (@ref last_command).
     */
    constexpr void back_calculate(T u_sat, T u_ff = T{0}) { detail::pid_wind_back(*this, Ki, u_sat, u_cmd + u_ff); }

    /** @brief The command the last control() returned. */
    [[nodiscard]] constexpr T last_command() const { return u_cmd; }

    /** @brief Clear integral and derivative state; keeps gains and runtime mode. */
    constexpr void reset() {
        integral = T{0};
        prev_x = T{0};
        deriv = T{0};
        first = true;
    }

    constexpr void enable() { runtime_mode = PIDRuntimeMode::Auto; }

    constexpr void disable(T track) {
        runtime_mode = PIDRuntimeMode::Tracking;
        u_track = track;
    }

    [[nodiscard]] constexpr bool is_enabled() const { return runtime_mode == PIDRuntimeMode::Auto; }
};

/**
 * @ingroup discrete_controllers
 * @brief Fixed-rate discrete 2-DOF PI controller specialization
 *
 * Constructing from a design with a derivative term (Kd ≠ 0), or with a Kbc of
 * the opposite sign to Ki, is refused: a compile error under constant
 * evaluation (`constinit` deploy), a trap at runtime.
 */
template<typename T>
struct PIDController<T, PIDMode::PI> {
    T            Kp{};
    T            Ki{};           ///< Discrete integral gain Kᵢ T_s (I += Ki·e)
    Bounds<1, T> u_lim{};        ///< Output limits [lower, upper] (default unbounded)
    Bounds<1, T> integral_lim{}; ///< Clamp on integral term I (same units as u; default unbounded)
    T            Kbc = T{0};     ///< Back-calculation: ΔI += Ki·(u−u_unsat)/Kbc; 0 = output clamp only
    T            b = T{1};       ///< Proportional setpoint weight
    T            Ts = T{0};      ///< Sample period from discretize (documentation)

    T integral = T{0}; ///< Integral term I (same units as u)

    PIDRuntimeMode runtime_mode{PIDRuntimeMode::Auto};
    T              u_track{T{0}};
    T              u_cmd{T{0}}; ///< Last command control() returned: the u in back_calculate's u_sat − (u + u_ff)

    constexpr PIDController() noexcept = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr PIDController(const design::DiscretePIDGains<U>& r) // NOLINT(google-explicit-constructor): assign-from-design
        : PIDController(r.template as<T>()) {}

    constexpr PIDController(const design::DiscretePIDGains<T>& r) // NOLINT(google-explicit-constructor): assign-from-design
        : Kp(r.Kp),
          Ki(r.Ki),
          u_lim(r.u_lim),
          integral_lim(r.integral_lim),
          Kbc(r.Kbc),
          b(r.b),
          Ts(r.Ts) {
        if (!r.success && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(r.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
        if (r.Kd != T{0}) {
            detail::refuse_failed_design(); // a PI cannot realize the D term
        }
        if (detail::pid_kbc_winds_up(r.Ki, r.Kbc)) {
            detail::refuse_failed_design(design::DesignStatus::invalid_parameter); // Kbc against Ki's sign winds the integrator up
        }
    }

    template<typename U>
    constexpr explicit PIDController(const PIDController<U, PIDMode::PI>& o)
        : Kp(static_cast<T>(o.Kp)),
          Ki(static_cast<T>(o.Ki)),
          u_lim(o.u_lim.template as<T>()),
          integral_lim(o.integral_lim.template as<T>()),
          Kbc(static_cast<T>(o.Kbc)),
          b(static_cast<T>(o.b)),
          Ts(static_cast<T>(o.Ts)),
          integral(static_cast<T>(o.integral)),
          runtime_mode(o.runtime_mode),
          u_track(static_cast<T>(o.u_track)),
          u_cmd(static_cast<T>(o.u_cmd)) {}

    [[nodiscard]] constexpr T control(T r, T y) {
        const T p = Kp * ((b * r) - y);
        return detail::pid_step(*this, Ki * (r - y), p, p, Ki);
    }

    /**
     * @brief Hold in Tracking mode, preloading the integrator to emit @p u_track bumplessly
     * @param u_track_ Output value to follow
     * @param y        Current measurement (preload assumes r = y, zero error)
     */
    constexpr void track(T u_track_, T y) {
        disable(u_track_);
        (void)control(y, y);
    }

    /**
     * @brief Anti-windup after a downstream clamp: ΔI = Ki·(u_sat − (u + u_ff))/Kbc
     *
     * Call after the stage that clamps the summed command. @p u_sat is what the plant
     * received, @p u_ff whatever was added to this law's command before the clamp; u is
     * the command control() returned (@ref last_command).
     */
    constexpr void back_calculate(T u_sat, T u_ff = T{0}) { detail::pid_wind_back(*this, Ki, u_sat, u_cmd + u_ff); }

    /** @brief The command the last control() returned. */
    [[nodiscard]] constexpr T last_command() const { return u_cmd; }

    /** @brief Clear the integral; keeps gains and runtime mode. */
    constexpr void reset() { integral = T{0}; }

    constexpr void enable() { runtime_mode = PIDRuntimeMode::Auto; }

    constexpr void disable(T track) {
        runtime_mode = PIDRuntimeMode::Tracking;
        u_track = track;
    }

    [[nodiscard]] constexpr bool is_enabled() const { return runtime_mode == PIDRuntimeMode::Auto; }
};

/**
 * @ingroup discrete_controllers
 * @brief Fixed-rate discrete 2-DOF PD controller specialization
 *
 * P plus filtered derivative on (c·r − y). With no integrator there is nothing
 * to preload, so Tracking holds u_track and re-enable resumes the plain law.
 * Constructing from a design with an integral term (Ki ≠ 0) is refused.
 */
template<typename T>
struct PIDController<T, PIDMode::PD> {
    T            Kp{};
    T            Kd{};  ///< Continuous Kd retained for inspection; tick uses d_a / d_b
    T            d_a{}; ///< D ← d_a D + d_b Δx
    T            d_b{};
    Bounds<1, T> u_lim{};   ///< Output limits [lower, upper] (default unbounded)
    T            b = T{1};  ///< Proportional setpoint weight
    T            c = T{1};  ///< Derivative setpoint weight
    T            Ts = T{0}; ///< Sample period from discretize (documentation)

    T    prev_x = T{0}; ///< Previous c·r − y
    T    deriv = T{0};  ///< Filtered derivative term D
    bool first = true;  ///< Next tick seeds prev_x (no derivative kick)

    PIDRuntimeMode runtime_mode{PIDRuntimeMode::Auto};
    T              u_track{T{0}};

    constexpr PIDController() noexcept = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr PIDController(const design::DiscretePIDGains<U>& r) // NOLINT(google-explicit-constructor): assign-from-design
        : PIDController(r.template as<T>()) {}

    constexpr PIDController(const design::DiscretePIDGains<T>& r) // NOLINT(google-explicit-constructor): assign-from-design
        : Kp(r.Kp),
          Kd(r.Kd),
          d_a(r.d_a),
          d_b(r.d_b),
          u_lim(r.u_lim),
          b(r.b),
          c(r.c),
          Ts(r.Ts) {
        if (!r.success && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(r.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
        if (r.Ki != T{0}) {
            detail::refuse_failed_design(); // a PD cannot realize the I term
        }
    }

    template<typename U>
    constexpr explicit PIDController(const PIDController<U, PIDMode::PD>& o)
        : Kp(static_cast<T>(o.Kp)),
          Kd(static_cast<T>(o.Kd)),
          d_a(static_cast<T>(o.d_a)),
          d_b(static_cast<T>(o.d_b)),
          u_lim(o.u_lim.template as<T>()),
          b(static_cast<T>(o.b)),
          c(static_cast<T>(o.c)),
          Ts(static_cast<T>(o.Ts)),
          prev_x(static_cast<T>(o.prev_x)),
          deriv(static_cast<T>(o.deriv)),
          first(o.first),
          runtime_mode(o.runtime_mode),
          u_track(static_cast<T>(o.u_track)) {}

    [[nodiscard]] constexpr T control(T r, T y) {
        deriv = (d_a * deriv) + (d_b * detail::pid_delta(*this, (c * r) - y));
        const T cmd = (runtime_mode == PIDRuntimeMode::Tracking) ? u_track : (Kp * ((b * r) - y)) + deriv;
        return u_lim.saturate(cmd);
    }

    /** @brief Clear derivative state; keeps gains and runtime mode. */
    constexpr void reset() {
        prev_x = T{0};
        deriv = T{0};
        first = true;
    }

    constexpr void enable() { runtime_mode = PIDRuntimeMode::Auto; }

    constexpr void disable(T track) {
        runtime_mode = PIDRuntimeMode::Tracking;
        u_track = track;
    }

    [[nodiscard]] constexpr bool is_enabled() const { return runtime_mode == PIDRuntimeMode::Auto; }
};

/**
 * @ingroup discrete_controllers
 * @brief Discrete proportional controller specialization
 *
 * Constructing from a design with an integral or derivative term is refused.
 */
template<typename T>
struct PIDController<T, PIDMode::P> {
    T            Kp{};
    Bounds<1, T> u_lim{};  ///< Output limits [lower, upper] (default unbounded)
    T            b = T{1}; ///< Proportional setpoint weight

    PIDRuntimeMode runtime_mode{PIDRuntimeMode::Auto};
    T              u_track{T{0}};

    constexpr PIDController() noexcept = default;

    constexpr explicit PIDController(T Kp_) noexcept
        : Kp(Kp_) {}

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr PIDController(const design::DiscretePIDGains<U>& r) // NOLINT(google-explicit-constructor): assign-from-design
        : PIDController(r.template as<T>()) {}

    constexpr PIDController(const design::DiscretePIDGains<T>& r) // NOLINT(google-explicit-constructor): assign-from-design
        : Kp(r.Kp), u_lim(r.u_lim), b(r.b) {
        if (!r.success && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(r.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
        if (r.Ki != T{0} || r.Kd != T{0}) {
            detail::refuse_failed_design(); // a P cannot realize the I or D term
        }
    }

    template<typename U>
    constexpr explicit PIDController(const PIDController<U, PIDMode::P>& o)
        : Kp(static_cast<T>(o.Kp)),
          u_lim(o.u_lim.template as<T>()),
          b(static_cast<T>(o.b)),
          runtime_mode(o.runtime_mode),
          u_track(static_cast<T>(o.u_track)) {}

    [[nodiscard]] constexpr T control(T r, T y) {
        const T cmd = (runtime_mode == PIDRuntimeMode::Tracking) ? u_track : Kp * ((b * r) - y);
        return u_lim.saturate(cmd);
    }

    constexpr void reset() {}

    constexpr void enable() { runtime_mode = PIDRuntimeMode::Auto; }

    constexpr void disable(T track) {
        runtime_mode = PIDRuntimeMode::Tracking;
        u_track = track;
    }

    [[nodiscard]] constexpr bool is_enabled() const { return runtime_mode == PIDRuntimeMode::Auto; }
};

/**
 * @brief Deduce the runtime from its gains, keeping their scalar: PIDController pid{gains} is a PID;
 * through an alias the mode comes from the alias: PIController pi{gains};
 */
template<typename T, PIDMode Mode = PIDMode::PID>
PIDController(const design::DiscretePIDGains<T>&) -> PIDController<T, Mode>;

template<typename T = float>
using PController = PIDController<T, PIDMode::P>;

template<typename T = float>
using PIController = PIDController<T, PIDMode::PI>;

template<typename T = float>
using PDController = PIDController<T, PIDMode::PD>;

/**
 * @brief Continuous-gain PID with per-tick sample time (variable-rate secondary form)
 *
 * Stores continuous-time @f$K_p, K_i, K_d@f$ and applies backward-Euler
 * integration / filtered differentiation using the measured @p Ts each call:
 * `control(r, y, Ts)`. Use this when the loop period is not fixed at design
 * time (jitter, multi-rate motor loops, PLL step periods).
 *
 * Integral term @f$I@f$ with @f$\dot I = K_i e@f$ and @f$u_I = I@f$, matching
 * the discrete form after @ref design::PIDGains::discretize.
 *
 * For fixed-rate embedded deploy, prefer
 * `design::pid(...).discretize(Ts)` → PIDController.
 *
 * Gains with a non-finite entry, a Tf < 0, or a Kbc of the opposite sign to Ki are
 * refused at construction: ill-formed under constant evaluation, a trap at runtime.
 *
 * @code
 * ContinuousPID<float> c{design::pid(2.0f, 5.0f, 0.1f)};
 * float u = c.control(r, y, measured_dt);
 * @endcode
 *
 * @see PIDController for the canonical fixed-rate path
 * @see Åström & Hägglund, "Advanced PID Control" (2006), Sec. 3.3–4.4
 */
template<typename T = float>
struct ContinuousPID {
    T            Kp{};
    T            Ki{}; ///< Continuous integral gain (I += Ki·e·Ts)
    T            Kd{};
    T            Tf = T{0};      ///< Derivative filter time constant (0 = unfiltered)
    Bounds<1, T> u_lim{};        ///< Output limits [lower, upper] (default unbounded)
    Bounds<1, T> integral_lim{}; ///< Clamp on integral term I (same units as u; default unbounded)
    T            Kbc = T{0};     ///< Back-calculation: ΔI += Ki·(Ts/Kbc)(u−u_unsat); 0 = output clamp only
    T            b = T{1};       ///< Proportional setpoint weight
    T            c = T{1};       ///< Derivative setpoint weight

    T    integral = T{0}; ///< Integral term I (same units as u)
    T    prev_x = T{0};   ///< Previous c·r − y
    T    deriv = T{0};    ///< Filtered derivative term D
    bool first = true;    ///< Next tick seeds prev_x (no derivative kick)

    PIDRuntimeMode runtime_mode{PIDRuntimeMode::Auto};
    T              u_track{T{0}};
    T              u_cmd{T{0}};   ///< Last command control() returned: the u in back_calculate's u_sat − (u + u_ff)
    T              Ts_last{T{0}}; ///< [s] Period of the last control() tick, which back_calculate charges

    constexpr ContinuousPID() noexcept = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr ContinuousPID(const design::PIDGains<U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : ContinuousPID(result.template as<T>()) {}

    constexpr ContinuousPID(const design::PIDGains<T>& r) // NOLINT(google-explicit-constructor): assign-from-design
        : Kp(r.Kp),
          Ki(r.Ki),
          Kd(r.Kd),
          Tf(r.Tf),
          u_lim(r.u_lim),
          integral_lim(r.integral_lim),
          Kbc(r.Kbc),
          b(r.b),
          c(r.c) {
        if (!r.success && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(r.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
        if (!detail::pid_gains_in_domain(r.Kp, r.Ki, r.Kd, r.Tf, r.Kbc, r.b, r.c) || detail::pid_kbc_winds_up(r.Ki, r.Kbc)) {
            detail::refuse_failed_design(design::DesignStatus::invalid_parameter); // a non-finite gain, Tf < 0, or Kbc against Ki's sign
        }
    }

    template<typename U>
    constexpr explicit ContinuousPID(const ContinuousPID<U>& o)
        : Kp(static_cast<T>(o.Kp)),
          Ki(static_cast<T>(o.Ki)),
          Kd(static_cast<T>(o.Kd)),
          Tf(static_cast<T>(o.Tf)),
          u_lim(o.u_lim.template as<T>()),
          integral_lim(o.integral_lim.template as<T>()),
          Kbc(static_cast<T>(o.Kbc)),
          b(static_cast<T>(o.b)),
          c(static_cast<T>(o.c)),
          integral(static_cast<T>(o.integral)),
          prev_x(static_cast<T>(o.prev_x)),
          deriv(static_cast<T>(o.deriv)),
          first(o.first),
          runtime_mode(o.runtime_mode),
          u_track(static_cast<T>(o.u_track)),
          u_cmd(static_cast<T>(o.u_cmd)),
          Ts_last(static_cast<T>(o.Ts_last)) {}

    /**
     * @brief Compute 2-DOF PID control output with measured sample period
     *
     * A @p Ts that is not positive and finite holds state and emits the current command (no
     * division by zero on the derivative path); a back_calculate after such a tick is a no-op.
     */
    [[nodiscard]] constexpr T control(T r, T y, T Ts) {
        const T p = Kp * ((b * r) - y);
        if (!damp::finite_positive(Ts)) {
            Ts_last = T{0};
            u_cmd = u_lim.saturate(runtime_mode == PIDRuntimeMode::Tracking ? u_track : p + integral + deriv);
            return u_cmd;
        }
        Ts_last = Ts;

        // Filtered derivative D = Kd*s/(1 + Tf*s) on (c*r − y), backward Euler
        const T decay = Tf / (Ts + Tf);
        deriv = (decay * deriv) + ((Kd / (Ts + Tf)) * detail::pid_delta(*this, (c * r) - y));
        return detail::pid_step(*this, Ki * (r - y) * Ts, p + deriv, p + (decay * deriv), Ki * Ts);
    }

    /**
     * @brief Hold in Tracking mode, preloading the integrator to emit @p u_track bumplessly
     *
     * The preload assumes the next tick has the same @p Ts.
     */
    constexpr void track(T u_track_, T y, T Ts) {
        disable(u_track_);
        (void)control(y, y, Ts);
    }

    /**
     * @brief Anti-windup after a downstream clamp: ΔI = Ki·Ts·(u_sat − (u + u_ff))/Kbc
     *
     * Call after the stage that clamps the summed command, in the same tick as control():
     * Ts is that tick's period and u the command it returned (@ref last_command).
     * @p u_ff is whatever was added to this law's command before the clamp.
     */
    constexpr void back_calculate(T u_sat, T u_ff = T{0}) { detail::pid_wind_back(*this, Ki * Ts_last, u_sat, u_cmd + u_ff); }

    /** @brief The command the last control() returned. */
    [[nodiscard]] constexpr T last_command() const { return u_cmd; }

    /** @brief Clear integral and derivative state; keeps gains and runtime mode. */
    constexpr void reset() {
        integral = T{0};
        prev_x = T{0};
        deriv = T{0};
        first = true;
    }

    constexpr void enable() { runtime_mode = PIDRuntimeMode::Auto; }

    constexpr void disable(T track) {
        runtime_mode = PIDRuntimeMode::Tracking;
        u_track = track;
    }

    [[nodiscard]] constexpr bool is_enabled() const { return runtime_mode == PIDRuntimeMode::Auto; }
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<typename T>
ContinuousPID(const design::PIDGains<T>&) -> ContinuousPID<T>;

} // namespace damp
