// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file pid.hpp
 * @brief PID design and runtime controllers
 */

#include <cstdint>
#include <limits>

#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

namespace damp {
namespace design {

/**
 * @struct DiscretePIDResult
 * @brief Fixed-rate discrete PID coefficients (canonical deploy form)
 *
 * Produced by PIDResult::discretize. Sample time is baked into the
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
 * @see PIDResult::discretize
 * @see Åström & Hägglund, "Advanced PID Control" (2006), Sec. 3.3
 */
template<typename T = double>
struct DiscretePIDResult {
    T Kp{};
    T Ki{};  ///< Discrete integral gain Kᵢ T_s (I += Ki·e)
    T Kd{};  ///< Continuous K_d (documentation / re-discretize); runtime uses @ref d_a / @ref d_b
    T d_a{}; ///< Filtered-derivative pole: D ← d_a D + d_b Δx
    T d_b{}; ///< Filtered-derivative input gain
    T u_min = -std::numeric_limits<T>::max();
    T u_max = std::numeric_limits<T>::max();
    T i_min = -std::numeric_limits<T>::max(); ///< Clamp on integral term I (same units as u)
    T i_max = std::numeric_limits<T>::max();
    T Kbc = T{0}; ///< Back-calculation coefficient; 0 = output clamp only
    T b = T{1};   ///< Proportional setpoint weight
    T c = T{1};   ///< Derivative setpoint weight
    T Tf = T{0};  ///< Continuous derivative filter time constant (documentation)
    T Ts = T{0};  ///< Sample period used for discretization [s]

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return DiscretePIDResult<U>{
            static_cast<U>(Kp),
            static_cast<U>(Ki),
            static_cast<U>(Kd),
            static_cast<U>(d_a),
            static_cast<U>(d_b),
            static_cast<U>(u_min),
            static_cast<U>(u_max),
            static_cast<U>(i_min),
            static_cast<U>(i_max),
            static_cast<U>(Kbc),
            static_cast<U>(b),
            static_cast<U>(c),
            static_cast<U>(Tf),
            static_cast<U>(Ts),
        };
    }
};

/**
 * @struct PIDResult
 * @brief 2-DOF continuous-time PID controller design result
 *
 * Contains continuous gains and setpoint weights. Deploy via
 * discretize into a @ref DiscretePIDResult for the fixed-rate
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
struct PIDResult {
    T Kp{};
    T Ki{};
    T Kd{};
    T u_min = -std::numeric_limits<T>::max();
    T u_max = std::numeric_limits<T>::max();
    T i_min = -std::numeric_limits<T>::max(); ///< Clamp on integral term I (same units as u)
    T i_max = std::numeric_limits<T>::max();
    T Kbc = T{0}; ///< Back-calculation coefficient: ΔI += Ki·(Ts/Kbc)(u−u_unsat).
                  ///< Tracking time T_t = Kbc/Ki (parallel form); Kbc = Kp ⇒ T_t = T_i.
                  ///< Larger Kbc → slower unwind; 0 = output clamp only.
    T b = T{1};   ///< Proportional setpoint weight (0=I-PD, 1=standard PID)
    T c = T{1};   ///< Derivative setpoint weight   (0=PI-D, 1=standard PID)
    T Tf = T{0};  ///< Derivative filter time constant (0 = unfiltered)

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return PIDResult<U>{
            static_cast<U>(Kp),
            static_cast<U>(Ki),
            static_cast<U>(Kd),
            static_cast<U>(u_min),
            static_cast<U>(u_max),
            static_cast<U>(i_min),
            static_cast<U>(i_max),
            static_cast<U>(Kbc),
            static_cast<U>(b),
            static_cast<U>(c),
            static_cast<U>(Tf),
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
     * @param ts Sample period [s]; must be positive for a useful integral / derivative map
     * @return DiscretePIDResult for PIDController
     */
    [[nodiscard]] constexpr DiscretePIDResult<T> discretize(T ts) const {
        DiscretePIDResult<T> d{};
        d.Kp = Kp;
        d.Kd = Kd;
        d.u_min = u_min;
        d.u_max = u_max;
        d.Kbc = Kbc;
        d.b = b;
        d.c = c;
        d.Tf = Tf;
        d.Ts = ts;
        d.i_min = i_min;
        d.i_max = i_max;
        if (ts > T{0}) {
            d.Ki = Ki * ts;
            const T den = ts + Tf;
            d.d_a = Tf / den;
            d.d_b = Kd / den;
        } else {
            d.Ki = T{0};
            d.d_a = T{0};
            d.d_b = T{0};
        }
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
            .num = {Ki, Kp + (Ki * Tf), (Kp * Tf) + Kd},
            .den = {T{0}, T{1}, Tf},
        };
    }

    /**
     * @brief Setpoint-channel C_ff(s) = b Kp + Ki/s + c Kd s/(1 + Tf s)
     *
     * Tracking is T_yr = C_ff P / (1 + C_fb P). Equals @ref to_tf when b = c = 1.
     */
    [[nodiscard]] constexpr TransferFunction<3, 3, T> to_tf_ff() const {
        return TransferFunction<3, 3, T>{
            .num = {Ki, (b * Kp) + (Ki * Tf), (b * Kp * Tf) + (c * Kd)},
            .den = {T{0}, T{1}, Tf},
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
     * @brief Error-channel PI realization C(s) = Kₚ + Kᵢ/s
     *
     * One integrator. @f$ K_d @f$ / @f$ T_f @f$ are not in this map — use
     * @ref to_tf for the full (possibly improper) C(s). This is the LTI
     * composed into cascade and feedforward plants.
     */
    [[nodiscard]] constexpr StateSpace<1, 1, 1, T> to_ss() const {
        return StateSpace<1, 1, 1, T>{
            .A = {{T{0}}},
            .B = {{T{1}}},
            .C = {{Ki}},
            .D = {{Kp}},
        };
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
 * @brief 2-DOF continuous PID controller design
 *
 * Constructs a PIDResult with the given continuous-time gains and optional
 * setpoint weights. Deploy with PIDResult::discretize for fixed-rate
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
 * @endcode
 *
 * @param Kp    Proportional gain
 * @param Ki    Integral gain (continuous)
 * @param Kd    Derivative gain
 * @param u_min Minimum control output
 * @param u_max Maximum control output
 * @param i_min Minimum integral term (same units as u)
 * @param i_max Maximum integral term
 * @param Kbc   Back-calculation coefficient: ΔI += Ki·(Ts/Kbc)(u−u_unsat); T_t = Kbc/Ki
 *              (set Kbc = Kp for T_t = T_i). 0 = clamping only
 * @param b     Proportional setpoint weight (default: 1 - standard PID)
 * @param c     Derivative setpoint weight   (default: 1 - standard PID)
 * @param Tf    Derivative filter time constant (default: 0 - unfiltered)
 *
 * @return PIDResult with the specified parameters
 */
template<typename T = double>
[[nodiscard]] constexpr PIDResult<T> pid(
    T Kp,
    T Ki,
    T Kd,
    T u_min = -std::numeric_limits<T>::max(),
    T u_max = std::numeric_limits<T>::max(),
    T i_min = -std::numeric_limits<T>::max(),
    T i_max = std::numeric_limits<T>::max(),
    T Kbc = T{0},
    T b = T{1},
    T c = T{1},
    T Tf = T{0}
) {
    return PIDResult<T>{Kp, Ki, Kd, u_min, u_max, i_min, i_max, Kbc, b, c, Tf};
}

} // namespace design

/**
 * @brief Compile-time selection of the PID control-law structure.
 *
 * Picks which specialization of PIDController is instantiated, so each mode
 * stores and computes only its own terms. Distinct from @ref PIDRuntimeMode,
 * which selects runtime behavior (Auto vs Tracking) at each tick.
 */
enum class PIDMode : std::uint8_t {
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
 * - `Tracking`: the controller's output is `clamp(u_track, u_min, u_max)`
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
enum class PIDRuntimeMode : std::uint8_t {
    Auto,
    Tracking,
};

namespace detail {

// Integral + output stage shared by PIDController and ContinuousPID. Free
// functions over the controller's public fields, so neither type needs a base.

/// I ← clamp(v, i_min, i_max)
template<typename C, typename T>
constexpr void pid_set_integral(C& c, T v) {
    c.integral = damp::clamp(v, c.i_min, c.i_max);
}

/// Back-calculation ΔI = gain·(u_sat − u_unsat)/Kbc. No-op when Kbc or gain is zero.
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
        return damp::clamp(c.u_track, c.u_min, c.u_max);
    }
    pid_set_integral(c, c.integral + dI);
    const T u_unsat = pd + c.integral;
    const T u = damp::clamp(u_unsat, c.u_min, c.u_max);
    pid_wind_back(c, gain, u, u_unsat);
    return u;
}

/// Δ(c·r − y) since the last tick; zero on the first (no derivative kick).
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
 * @ref design::PIDResult::discretize. One tick is `control(r, y)` — no `Ts`
 * argument. Matches the library-wide pre-discretized controller model
 * (LQR, Kalman, PR, lead-lag).
 *
 * @ref PIDMode picks the specialization (P, PI, PD, PID); each stores and
 * computes only its own terms. A design result whose gains the mode cannot
 * realize is refused at construction: ill-formed under constant evaluation
 * (`constinit` deploy), a trap at runtime.
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
 * MATLAB® equivalent: design with `pid(...)` then `c2d`; runtime is the discrete law.
 * @see Åström & Hägglund, "Advanced PID Control" (2006), Sec. 3.3–4.4
 *
 * @note Gains from a continuous-time designer (ziegler_nichols, cohen_coon,
 *       simc, lambda_tuning, pid_from_bandwidth) are Ts-agnostic until
 *       discretize(Ts). Gains from design::pid_pole_placement already embed a
 *       design Ts in continuous form — discretize with that same Ts.
 */
template<typename T>
struct PIDController<T, PIDMode::PID> {
    T Kp{};
    T Ki{};  ///< Discrete integral gain Kᵢ T_s (I += Ki·e)
    T Kd{};  ///< Continuous Kd retained for inspection; tick uses d_a / d_b
    T d_a{}; ///< D ← d_a D + d_b Δx
    T d_b{};
    T u_min = -std::numeric_limits<T>::max();
    T u_max = std::numeric_limits<T>::max();
    T i_min = -std::numeric_limits<T>::max(); ///< Clamp on integral term I (same units as u)
    T i_max = std::numeric_limits<T>::max();
    T Kbc = T{0}; ///< Back-calculation: ΔI += Ki·(u−u_unsat)/Kbc; 0 = output clamp only
    T b = T{1};   ///< Proportional setpoint weight
    T c = T{1};   ///< Derivative setpoint weight
    T Ts = T{0};  ///< Sample period from discretize (documentation)

    T    integral = T{0}; ///< Integral term I (same units as u)
    T    prev_x = T{0};   ///< Previous c·r − y
    T    deriv = T{0};    ///< Filtered derivative term D
    bool first = true;    ///< Next tick seeds prev_x (no derivative kick)

    PIDRuntimeMode runtime_mode{PIDRuntimeMode::Auto};
    T              u_track{T{0}};

    constexpr PIDController() noexcept = default;

    /// From discrete design result. Non-explicit for assign-from-design:
    /// `PIDController c = design::pid(...).discretize(Ts).as<float>();`
    constexpr PIDController(const design::DiscretePIDResult<T>& r)
        : Kp(r.Kp),
          Ki(r.Ki),
          Kd(r.Kd),
          d_a(r.d_a),
          d_b(r.d_b),
          u_min(r.u_min),
          u_max(r.u_max),
          i_min(r.i_min),
          i_max(r.i_max),
          Kbc(r.Kbc),
          b(r.b),
          c(r.c),
          Ts(r.Ts) {}

    template<typename U>
    constexpr explicit PIDController(const PIDController<U, PIDMode::PID>& o)
        : Kp(static_cast<T>(o.Kp)),
          Ki(static_cast<T>(o.Ki)),
          Kd(static_cast<T>(o.Kd)),
          d_a(static_cast<T>(o.d_a)),
          d_b(static_cast<T>(o.d_b)),
          u_min(static_cast<T>(o.u_min)),
          u_max(static_cast<T>(o.u_max)),
          i_min(static_cast<T>(o.i_min)),
          i_max(static_cast<T>(o.i_max)),
          Kbc(static_cast<T>(o.Kbc)),
          b(static_cast<T>(o.b)),
          c(static_cast<T>(o.c)),
          Ts(static_cast<T>(o.Ts)),
          integral(static_cast<T>(o.integral)),
          prev_x(static_cast<T>(o.prev_x)),
          deriv(static_cast<T>(o.deriv)),
          first(o.first),
          runtime_mode(o.runtime_mode),
          u_track(static_cast<T>(o.u_track)) {}

    /**
     * @brief One fixed-rate discrete PID tick
     *
     * In Auto mode runs the 2-DOF law with back-calculation anti-windup.
     * In Tracking mode the output is `clamp(u_track, u_min, u_max)` and the
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

    /// Anti-windup hook driven by a downstream stage: ΔI = Ki·(u_sat − u_unsat)/Kbc.
    constexpr void back_calculate(T u_unsat, T u_sat) { detail::pid_wind_back(*this, Ki, u_sat, u_unsat); }

    /// Clear integral and derivative state; keeps gains and runtime mode.
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
 * Constructing from a design with a derivative term (Kd ≠ 0) is refused: a
 * compile error under constant evaluation (`constinit` deploy), a trap at
 * runtime.
 */
template<typename T>
struct PIDController<T, PIDMode::PI> {
    T Kp{};
    T Ki{}; ///< Discrete integral gain Kᵢ T_s (I += Ki·e)
    T u_min = -std::numeric_limits<T>::max();
    T u_max = std::numeric_limits<T>::max();
    T i_min = -std::numeric_limits<T>::max(); ///< Clamp on integral term I (same units as u)
    T i_max = std::numeric_limits<T>::max();
    T Kbc = T{0}; ///< Back-calculation: ΔI += Ki·(u−u_unsat)/Kbc; 0 = output clamp only
    T b = T{1};   ///< Proportional setpoint weight
    T Ts = T{0};  ///< Sample period from discretize (documentation)

    T integral = T{0}; ///< Integral term I (same units as u)

    PIDRuntimeMode runtime_mode{PIDRuntimeMode::Auto};
    T              u_track{T{0}};

    constexpr PIDController() noexcept = default;

    constexpr PIDController(const design::DiscretePIDResult<T>& r)
        : Kp(r.Kp),
          Ki(r.Ki),
          u_min(r.u_min),
          u_max(r.u_max),
          i_min(r.i_min),
          i_max(r.i_max),
          Kbc(r.Kbc),
          b(r.b),
          Ts(r.Ts) {
        if (r.Kd != T{0}) {
            detail::refuse_failed_design(); // a PI cannot realize the D term
        }
    }

    template<typename U>
    constexpr explicit PIDController(const PIDController<U, PIDMode::PI>& o)
        : Kp(static_cast<T>(o.Kp)),
          Ki(static_cast<T>(o.Ki)),
          u_min(static_cast<T>(o.u_min)),
          u_max(static_cast<T>(o.u_max)),
          i_min(static_cast<T>(o.i_min)),
          i_max(static_cast<T>(o.i_max)),
          Kbc(static_cast<T>(o.Kbc)),
          b(static_cast<T>(o.b)),
          Ts(static_cast<T>(o.Ts)),
          integral(static_cast<T>(o.integral)),
          runtime_mode(o.runtime_mode),
          u_track(static_cast<T>(o.u_track)) {}

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

    /// Anti-windup hook driven by a downstream stage: ΔI = Ki·(u_sat − u_unsat)/Kbc.
    constexpr void back_calculate(T u_unsat, T u_sat) { detail::pid_wind_back(*this, Ki, u_sat, u_unsat); }

    /// Clear the integral; keeps gains and runtime mode.
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
    T Kp{};
    T Kd{};  ///< Continuous Kd retained for inspection; tick uses d_a / d_b
    T d_a{}; ///< D ← d_a D + d_b Δx
    T d_b{};
    T u_min = -std::numeric_limits<T>::max();
    T u_max = std::numeric_limits<T>::max();
    T b = T{1};  ///< Proportional setpoint weight
    T c = T{1};  ///< Derivative setpoint weight
    T Ts = T{0}; ///< Sample period from discretize (documentation)

    T    prev_x = T{0}; ///< Previous c·r − y
    T    deriv = T{0};  ///< Filtered derivative term D
    bool first = true;  ///< Next tick seeds prev_x (no derivative kick)

    PIDRuntimeMode runtime_mode{PIDRuntimeMode::Auto};
    T              u_track{T{0}};

    constexpr PIDController() noexcept = default;

    constexpr PIDController(const design::DiscretePIDResult<T>& r)
        : Kp(r.Kp),
          Kd(r.Kd),
          d_a(r.d_a),
          d_b(r.d_b),
          u_min(r.u_min),
          u_max(r.u_max),
          b(r.b),
          c(r.c),
          Ts(r.Ts) {
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
          u_min(static_cast<T>(o.u_min)),
          u_max(static_cast<T>(o.u_max)),
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
        const T u = (runtime_mode == PIDRuntimeMode::Tracking) ? u_track : (Kp * ((b * r) - y)) + deriv;
        return damp::clamp(u, u_min, u_max);
    }

    /// Clear derivative state; keeps gains and runtime mode.
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
    T Kp{};
    T u_min = -std::numeric_limits<T>::max();
    T u_max = std::numeric_limits<T>::max();
    T b = T{1}; ///< Proportional setpoint weight

    PIDRuntimeMode runtime_mode{PIDRuntimeMode::Auto};
    T              u_track{T{0}};

    constexpr PIDController() noexcept = default;

    constexpr explicit PIDController(T Kp_) noexcept
        : Kp(Kp_) {}

    constexpr PIDController(const design::DiscretePIDResult<T>& r)
        : Kp(r.Kp), u_min(r.u_min), u_max(r.u_max), b(r.b) {
        if (r.Ki != T{0} || r.Kd != T{0}) {
            detail::refuse_failed_design(); // a P cannot realize the I or D term
        }
    }

    template<typename U>
    constexpr explicit PIDController(const PIDController<U, PIDMode::P>& o)
        : Kp(static_cast<T>(o.Kp)),
          u_min(static_cast<T>(o.u_min)),
          u_max(static_cast<T>(o.u_max)),
          b(static_cast<T>(o.b)),
          runtime_mode(o.runtime_mode),
          u_track(static_cast<T>(o.u_track)) {}

    [[nodiscard]] constexpr T control(T r, T y) {
        const T u = (runtime_mode == PIDRuntimeMode::Tracking) ? u_track : Kp * ((b * r) - y);
        return damp::clamp(u, u_min, u_max);
    }

    constexpr void reset() {}

    constexpr void enable() { runtime_mode = PIDRuntimeMode::Auto; }

    constexpr void disable(T track) {
        runtime_mode = PIDRuntimeMode::Tracking;
        u_track = track;
    }

    [[nodiscard]] constexpr bool is_enabled() const { return runtime_mode == PIDRuntimeMode::Auto; }
};

template<typename T>
PIDController(const design::DiscretePIDResult<T>&) -> PIDController<T, PIDMode::PID>;

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
 * the discrete form after @ref design::PIDResult::discretize.
 *
 * For fixed-rate embedded deploy, prefer
 * `design::pid(...).discretize(Ts)` → PIDController.
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
    T Kp{};
    T Ki{}; ///< Continuous integral gain (I += Ki·e·Ts)
    T Kd{};
    T Tf = T{0}; ///< Derivative filter time constant (0 = unfiltered)
    T u_min = -std::numeric_limits<T>::max();
    T u_max = std::numeric_limits<T>::max();
    T i_min = -std::numeric_limits<T>::max(); ///< Clamp on integral term I (same units as u)
    T i_max = std::numeric_limits<T>::max();
    T Kbc = T{0}; ///< Back-calculation: ΔI += Ki·(Ts/Kbc)(u−u_unsat); 0 = output clamp only
    T b = T{1};   ///< Proportional setpoint weight
    T c = T{1};   ///< Derivative setpoint weight

    T    integral = T{0}; ///< Integral term I (same units as u)
    T    prev_x = T{0};   ///< Previous c·r − y
    T    deriv = T{0};    ///< Filtered derivative term D
    bool first = true;    ///< Next tick seeds prev_x (no derivative kick)

    PIDRuntimeMode runtime_mode{PIDRuntimeMode::Auto};
    T              u_track{T{0}};

    constexpr ContinuousPID() noexcept = default;

    constexpr ContinuousPID(const design::PIDResult<T>& r) noexcept
        : Kp(r.Kp),
          Ki(r.Ki),
          Kd(r.Kd),
          Tf(r.Tf),
          u_min(r.u_min),
          u_max(r.u_max),
          i_min(r.i_min),
          i_max(r.i_max),
          Kbc(r.Kbc),
          b(r.b),
          c(r.c) {}

    template<typename U>
    constexpr explicit ContinuousPID(const ContinuousPID<U>& o)
        : Kp(static_cast<T>(o.Kp)),
          Ki(static_cast<T>(o.Ki)),
          Kd(static_cast<T>(o.Kd)),
          Tf(static_cast<T>(o.Tf)),
          u_min(static_cast<T>(o.u_min)),
          u_max(static_cast<T>(o.u_max)),
          i_min(static_cast<T>(o.i_min)),
          i_max(static_cast<T>(o.i_max)),
          Kbc(static_cast<T>(o.Kbc)),
          b(static_cast<T>(o.b)),
          c(static_cast<T>(o.c)),
          integral(static_cast<T>(o.integral)),
          prev_x(static_cast<T>(o.prev_x)),
          deriv(static_cast<T>(o.deriv)),
          first(o.first),
          runtime_mode(o.runtime_mode),
          u_track(static_cast<T>(o.u_track)) {}

    /**
     * @brief Compute 2-DOF PID control output with measured sample period
     *
     * Non-positive @p Ts holds state and emits the current command (no
     * division by zero on the derivative path).
     */
    [[nodiscard]] constexpr T control(T r, T y, T Ts) {
        const T p = Kp * ((b * r) - y);
        if (Ts <= T{0}) {
            return damp::clamp(runtime_mode == PIDRuntimeMode::Tracking ? u_track : p + integral + deriv, u_min, u_max);
        }

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

    /// Anti-windup hook: ΔI = Ki·Ts·(u_sat − u_unsat)/Kbc.
    constexpr void back_calculate(T u_unsat, T u_sat, T Ts) { detail::pid_wind_back(*this, Ki * Ts, u_sat, u_unsat); }

    /// Clear integral and derivative state; keeps gains and runtime mode.
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

} // namespace damp
