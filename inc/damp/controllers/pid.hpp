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
     * @brief Error-channel PI realization @f$ C(s) = K_p + K_i/s @f$
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
 * Picks which specialization of PIDController is instantiated. Distinct
 * from @ref PIDRuntimeMode, which selects runtime behavior (Auto vs Tracking)
 * at each tick.
 */
enum class PIDMode : std::uint8_t {
    P,
    PI,
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
    T i_min = -std::numeric_limits<T>::max();
    T i_max = std::numeric_limits<T>::max();
    T Kbc = T{0}; ///< Back-calculation: ΔI += Ki·(u−u_unsat)/Kbc
    T b = T{1};
    T c = T{1};
    T Ts = T{0}; ///< Sample period from discretize (documentation)

    T    integral = T{0};        ///< Integral term I (same units as u)
    T    prev_cr_minus_y = T{0}; ///< Previous (c*r − y) for derivative
    T    deriv = T{0};           ///< Filtered derivative term
    bool first_ = true;

    PIDRuntimeMode runtime_mode{PIDRuntimeMode::Auto};
    T              u_track{T{0}};

    constexpr PIDController() = default;

    /// From discrete design result. Non-explicit for assign-from-design:
    /// `PIController c = design::pid(...).discretize(Ts).as<float>();`
    constexpr PIDController(const design::DiscretePIDResult<T>& result)
        : Kp(result.Kp),
          Ki(result.Ki),
          Kd(result.Kd),
          d_a(result.d_a),
          d_b(result.d_b),
          u_min(result.u_min),
          u_max(result.u_max),
          i_min(result.i_min),
          i_max(result.i_max),
          Kbc(result.Kbc),
          b(result.b),
          c(result.c),
          Ts(result.Ts) {}

    template<typename U>
    constexpr explicit PIDController(const PIDController<U, PIDMode::PID>& other)
        : Kp(static_cast<T>(other.Kp)),
          Ki(static_cast<T>(other.Ki)),
          Kd(static_cast<T>(other.Kd)),
          d_a(static_cast<T>(other.d_a)),
          d_b(static_cast<T>(other.d_b)),
          u_min(static_cast<T>(other.u_min)),
          u_max(static_cast<T>(other.u_max)),
          i_min(static_cast<T>(other.i_min)),
          i_max(static_cast<T>(other.i_max)),
          Kbc(static_cast<T>(other.Kbc)),
          b(static_cast<T>(other.b)),
          c(static_cast<T>(other.c)),
          Ts(static_cast<T>(other.Ts)),
          integral(static_cast<T>(other.integral)),
          prev_cr_minus_y(static_cast<T>(other.prev_cr_minus_y)),
          deriv(static_cast<T>(other.deriv)),
          first_(other.first_),
          runtime_mode(other.runtime_mode),
          u_track(static_cast<T>(other.u_track)) {}

    /**
     * @brief One fixed-rate discrete PID tick
     *
     * In Auto mode runs the 2-DOF law with back-calculation anti-windup.
     * In Tracking mode the output is `clamp(u_track, u_min, u_max)` and the
     * integrator (plus derivative state) is pre-loaded for bumpless re-enable.
     */
    [[nodiscard]] constexpr T control(T r, T y) {
        const T cr_minus_y = (c * r) - y;
        if (first_) {
            prev_cr_minus_y = cr_minus_y;
            first_ = false;
        }
        const T dX = cr_minus_y - prev_cr_minus_y;
        prev_cr_minus_y = cr_minus_y;

        // Filtered (or raw when d_a = 0) derivative on (c*r − y)
        deriv = (d_a * deriv) + (d_b * dX);

        const T e = r - y;

        if (runtime_mode == PIDRuntimeMode::Tracking) {
            // Preload I so next Auto tick with the same (r,y) yields u_track.
            const T target = u_track - (Kp * ((b * r) - y)) - deriv - (Ki * e);
            integral = damp::clamp(target, i_min, i_max);
            return damp::clamp(u_track, u_min, u_max);
        }

        integral += Ki * e;
        integral = damp::clamp(integral, i_min, i_max);
        const T u_unsat = (Kp * ((b * r) - y)) + integral + deriv;
        const T u = damp::clamp(u_unsat, u_min, u_max);

        // Back-calculation on the term state: ΔI = Ki·(u − u_unsat)/Kbc
        if (Kbc != T{0} && Ki != T{0}) {
            integral += Ki * ((u - u_unsat) / Kbc);
            integral = damp::clamp(integral, i_min, i_max);
        }

        return u;
    }

    constexpr void reset() {
        integral = T{0};
        prev_cr_minus_y = T{0};
        deriv = T{0};
        first_ = true;
    }

    constexpr void enable() { runtime_mode = PIDRuntimeMode::Auto; }

    constexpr void disable(T track) {
        runtime_mode = PIDRuntimeMode::Tracking;
        u_track = track;
    }

    [[nodiscard]] constexpr bool is_enabled() const { return runtime_mode == PIDRuntimeMode::Auto; }

    /**
     * @brief Anti-windup hook driven by a downstream stage
     *
     * Winds the integral term by `Ki * (u_sat - u_unsat) / Kbc`. No-op when
     * `Kbc == 0` or `Ki == 0`.
     */
    constexpr void back_calculate(T u_unsat, T u_sat) {
        if (Kbc == T{0} || Ki == T{0} || u_unsat == u_sat) {
            return;
        }
        integral += Ki * ((u_sat - u_unsat) / Kbc);
        integral = damp::clamp(integral, i_min, i_max);
    }
};

/**
 * @ingroup discrete_controllers
 * @brief Fixed-rate discrete 2-DOF PI controller specialization
 */
template<typename T>
struct PIDController<T, PIDMode::PI> {
    T Kp{};
    T Ki{}; ///< Discrete integral gain Kᵢ T_s (I += Ki·e)
    T u_min = -std::numeric_limits<T>::max();
    T u_max = std::numeric_limits<T>::max();
    T i_min = -std::numeric_limits<T>::max();
    T i_max = std::numeric_limits<T>::max();
    T Kbc = T{0};
    T b = T{1};
    T Ts = T{0};

    T integral = T{0}; ///< Integral term I (same units as u)

    PIDRuntimeMode runtime_mode{PIDRuntimeMode::Auto};
    T              u_track{T{0}};

    constexpr PIDController() = default;

    constexpr PIDController(const design::DiscretePIDResult<T>& result)
        : Kp(result.Kp),
          Ki(result.Ki),
          u_min(result.u_min),
          u_max(result.u_max),
          i_min(result.i_min),
          i_max(result.i_max),
          Kbc(result.Kbc),
          b(result.b),
          Ts(result.Ts) {}

    template<typename U>
    constexpr explicit PIDController(const PIDController<U, PIDMode::PI>& other)
        : Kp(static_cast<T>(other.Kp)),
          Ki(static_cast<T>(other.Ki)),
          u_min(static_cast<T>(other.u_min)),
          u_max(static_cast<T>(other.u_max)),
          i_min(static_cast<T>(other.i_min)),
          i_max(static_cast<T>(other.i_max)),
          Kbc(static_cast<T>(other.Kbc)),
          b(static_cast<T>(other.b)),
          Ts(static_cast<T>(other.Ts)),
          integral(static_cast<T>(other.integral)),
          runtime_mode(other.runtime_mode),
          u_track(static_cast<T>(other.u_track)) {}

    [[nodiscard]] constexpr T control(T r, T y) {
        const T e = r - y;

        if (runtime_mode == PIDRuntimeMode::Tracking) {
            const T target = u_track - (Kp * ((b * r) - y)) - (Ki * e);
            integral = damp::clamp(target, i_min, i_max);
            return damp::clamp(u_track, u_min, u_max);
        }

        integral += Ki * e;
        integral = damp::clamp(integral, i_min, i_max);
        const T u_unsat = (Kp * ((b * r) - y)) + integral;
        const T u = damp::clamp(u_unsat, u_min, u_max);

        if (Kbc != T{0} && Ki != T{0}) {
            integral += Ki * ((u - u_unsat) / Kbc);
            integral = damp::clamp(integral, i_min, i_max);
        }
        return u;
    }

    constexpr void reset() { integral = T{0}; }

    constexpr void enable() { runtime_mode = PIDRuntimeMode::Auto; }
    constexpr void disable(T track) {
        runtime_mode = PIDRuntimeMode::Tracking;
        u_track = track;
    }

    /**
     * @brief Hold in Tracking mode, preloading the integrator to emit @p u_track bumplessly
     * @param u_track_ Output value to follow
     * @param y       Current measurement (preload assumes r = y, zero error)
     */
    constexpr void track(T u_track_, T y) {
        disable(u_track_);
        (void)control(y, y);
    }

    constexpr void back_calculate(T u_unsat, T u_sat) {
        if (Kbc == T{0} || Ki == T{0} || u_unsat == u_sat) {
            return;
        }
        integral += Ki * ((u_sat - u_unsat) / Kbc);
        integral = damp::clamp(integral, i_min, i_max);
    }

    [[nodiscard]] constexpr bool is_enabled() const { return runtime_mode == PIDRuntimeMode::Auto; }
};

/**
 * @ingroup discrete_controllers
 * @brief Discrete proportional controller specialization
 */
template<typename T>
struct PIDController<T, PIDMode::P> {
    T Kp{};
    T u_min = -std::numeric_limits<T>::max();
    T u_max = std::numeric_limits<T>::max();
    T b = T{1};

    PIDRuntimeMode runtime_mode{PIDRuntimeMode::Auto};
    T              u_track{T{0}};

    constexpr PIDController() = default;

    constexpr explicit PIDController(T Kp_)
        : Kp(Kp_) {}

    constexpr PIDController(const design::DiscretePIDResult<T>& result)
        : Kp(result.Kp), u_min(result.u_min), u_max(result.u_max), b(result.b) {}

    /// From continuous design (P has no rate-dependent term)
    constexpr PIDController(const design::PIDResult<T>& result)
        : Kp(result.Kp), u_min(result.u_min), u_max(result.u_max), b(result.b) {}

    template<typename U>
    constexpr explicit PIDController(const PIDController<U, PIDMode::P>& other)
        : Kp(static_cast<T>(other.Kp)),
          u_min(static_cast<T>(other.u_min)),
          u_max(static_cast<T>(other.u_max)),
          b(static_cast<T>(other.b)),
          runtime_mode(other.runtime_mode),
          u_track(static_cast<T>(other.u_track)) {}

    [[nodiscard]] constexpr T control(T r, T y) {
        if (runtime_mode == PIDRuntimeMode::Tracking) {
            return damp::clamp(u_track, u_min, u_max);
        }
        return damp::clamp(Kp * ((b * r) - y), u_min, u_max);
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
    T Ki{};
    T Kd{};
    T u_min = -std::numeric_limits<T>::max();
    T u_max = std::numeric_limits<T>::max();
    T i_min = -std::numeric_limits<T>::max();
    T i_max = std::numeric_limits<T>::max();
    T Kbc = T{0}; ///< Back-calculation: ΔI += Ki·(Ts/Kbc)(u−u_unsat)
    T b = T{1};
    T c = T{1};
    T Tf = T{0};

    T    integral = T{0}; ///< Integral term I (same units as u)
    T    prev_cr_minus_y = T{0};
    T    deriv = T{0};
    bool first_ = true;

    PIDRuntimeMode runtime_mode{PIDRuntimeMode::Auto};
    T              u_track{T{0}};

    constexpr ContinuousPID() = default;

    constexpr ContinuousPID(const design::PIDResult<T>& result)
        : Kp(result.Kp),
          Ki(result.Ki),
          Kd(result.Kd),
          u_min(result.u_min),
          u_max(result.u_max),
          i_min(result.i_min),
          i_max(result.i_max),
          Kbc(result.Kbc),
          b(result.b),
          c(result.c),
          Tf(result.Tf) {}

    template<typename U>
    constexpr explicit ContinuousPID(const ContinuousPID<U>& other)
        : Kp(static_cast<T>(other.Kp)),
          Ki(static_cast<T>(other.Ki)),
          Kd(static_cast<T>(other.Kd)),
          u_min(static_cast<T>(other.u_min)),
          u_max(static_cast<T>(other.u_max)),
          i_min(static_cast<T>(other.i_min)),
          i_max(static_cast<T>(other.i_max)),
          Kbc(static_cast<T>(other.Kbc)),
          b(static_cast<T>(other.b)),
          c(static_cast<T>(other.c)),
          Tf(static_cast<T>(other.Tf)),
          integral(static_cast<T>(other.integral)),
          prev_cr_minus_y(static_cast<T>(other.prev_cr_minus_y)),
          deriv(static_cast<T>(other.deriv)),
          first_(other.first_),
          runtime_mode(other.runtime_mode),
          u_track(static_cast<T>(other.u_track)) {}

    /**
     * @brief Compute 2-DOF PID control output with measured sample period
     *
     * Non-positive @p Ts holds state and emits the current command (no
     * division by zero on the derivative path).
     */
    [[nodiscard]] constexpr T control(T r, T y, T Ts) {
        if (Ts <= T{0}) {
            if (runtime_mode == PIDRuntimeMode::Tracking) {
                return damp::clamp(u_track, u_min, u_max);
            }
            return damp::clamp((Kp * ((b * r) - y)) + integral + deriv, u_min, u_max);
        }

        const T cr_minus_y = (c * r) - y;
        if (first_) {
            prev_cr_minus_y = cr_minus_y;
            first_ = false;
        }
        const T dX = cr_minus_y - prev_cr_minus_y;
        prev_cr_minus_y = cr_minus_y;

        // Filtered derivative D = Kd*s/(1 + Tf*s) on (c*r − y), backward Euler
        deriv = ((Tf * deriv) + (Kd * dX)) / (Ts + Tf);

        const T e = r - y;

        if (runtime_mode == PIDRuntimeMode::Tracking) {
            const T target = u_track - (Kp * ((b * r) - y)) - deriv - (Ki * e * Ts);
            integral = damp::clamp(target, i_min, i_max);
            return damp::clamp(u_track, u_min, u_max);
        }

        integral += Ki * e * Ts;
        integral = damp::clamp(integral, i_min, i_max);
        const T u_unsat = (Kp * ((b * r) - y)) + integral + deriv;
        const T u = damp::clamp(u_unsat, u_min, u_max);

        if (Kbc != T{0} && Ki != T{0}) {
            integral += Ki * (Ts / Kbc) * (u - u_unsat);
            integral = damp::clamp(integral, i_min, i_max);
        }

        return u;
    }

    constexpr void reset() {
        integral = T{0};
        prev_cr_minus_y = T{0};
        deriv = T{0};
        first_ = true;
    }

    constexpr void enable() { runtime_mode = PIDRuntimeMode::Auto; }

    constexpr void disable(T track) {
        runtime_mode = PIDRuntimeMode::Tracking;
        u_track = track;
    }

    /**
     * @brief Hold in Tracking mode, preloading the integrator to emit @p u_track bumplessly
     */
    constexpr void track(T u_track_, T y, T Ts) {
        disable(u_track_);
        (void)control(y, y, Ts);
    }

    [[nodiscard]] constexpr bool is_enabled() const { return runtime_mode == PIDRuntimeMode::Auto; }

    /**
     * @brief Anti-windup hook: ΔI += Ki · (u_sat − u_unsat) · Ts / Kbc
     */
    constexpr void back_calculate(T u_unsat, T u_sat, T Ts) {
        if (Kbc == T{0} || Ki == T{0} || u_unsat == u_sat) {
            return;
        }
        integral += Ki * Ts * ((u_sat - u_unsat) / Kbc);
        integral = damp::clamp(integral, i_min, i_max);
    }
};

} // namespace damp
