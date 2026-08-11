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
 * - @f$K_{i,d} = K_i T_s@f$ with error-sum state @f$I \leftarrow I + e@f$,
 *   contribution @f$K_{i,d}\,I@f$
 * - @f$D = a_d D + b_d\,\Delta x@f$ with
 *   @f$a_d = T_f/(T_s+T_f)@f$, @f$b_d = K_d/(T_s+T_f)@f$
 *   (@f$T_f=0@f$ ⇒ raw @f$K_d/T_s@f$ difference)
 * - Integrator limits scale as @f$i_{\min}/T_s@f$ (same units as the error sum)
 * - Back-calculation: @f$\Delta I = (u-u_{\mathrm{unsat}})/K_{bc}@f$ (Ts already in I)
 *
 * @see PIDResult::discretize
 * @see Åström & Hägglund, "Advanced PID Control" (2006), Sec. 3.3
 */
template<typename T = double>
struct DiscretePIDResult {
    T Kp{};
    T Ki{};  ///< Discrete integral gain Kᵢ T_s; I += e, u_I = Kᵢ I
    T Kd{};  ///< Continuous K_d (documentation / re-discretize); runtime uses @ref d_a / @ref d_b
    T d_a{}; ///< Filtered-derivative pole: D ← d_a D + d_b Δx
    T d_b{}; ///< Filtered-derivative input gain
    T u_min = -std::numeric_limits<T>::max();
    T u_max = std::numeric_limits<T>::max();
    T i_min = -std::numeric_limits<T>::max(); ///< Limit on discrete error-sum state
    T i_max = std::numeric_limits<T>::max();
    T Kbc = T{0}; ///< Back-calculation coefficient (continuous units); 0 = clamp-only
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
    T i_min = -std::numeric_limits<T>::max(); ///< Continuous integrator limit (∫e)
    T i_max = std::numeric_limits<T>::max();
    T Kbc = T{0}; ///< Back-calculation coefficient (not a pure time): ΔI += (Ts/Kbc)(u−u_unsat).
                  ///< Tracking time constant T_t = Kbc/Ki (parallel form). Kbc = Kp ⇒ T_t = T_i.
                  ///< Larger Kbc → slower unwind; 0 = clamp-only (no back-calculation).
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
        if (ts > T{0}) {
            d.Ki = Ki * ts;
            // Scale finite integrator limits only. ±max stays ±max so constinit
            // deploy paths (default unbounded I) remain constant-expression-safe
            // (max/Ts is not representable).
            constexpr T lim = std::numeric_limits<T>::max();
            d.i_min = (i_min > -lim) ? (i_min / ts) : -lim;
            d.i_max = (i_max < lim) ? (i_max / ts) : lim;
            const T den = ts + Tf;
            d.d_a = Tf / den;
            d.d_b = Kd / den;
        } else {
            d.Ki = T{0};
            d.i_min = i_min;
            d.i_max = i_max;
            d.d_a = T{0};
            d.d_b = T{0};
        }
        return d;
    }

    /**
     * @brief Continuous-time controller transfer function C(s).
     *
     * Returns @f$C(s) = K_p + K_i/s + K_d s/(1 + T_f s)@f$ as a 2nd-order TF in
     * ascending powers of s, so PID drops into the analysis tooling (Bode,
     * `series`/`feedback`, `discretize`) like the lead-lag / PR design results.
     * `Tf = 0` gives the ideal form @f$(K_d s^2 + K_p s + K_i)/s@f$.
     */
    [[nodiscard]] constexpr TransferFunction<3, 3, T> to_tf() const {
        return TransferFunction<3, 3, T>{
            .num = {Ki, Kp + (Ki * Tf), (Kp * Tf) + Kd},
            .den = {T{0}, T{1}, Tf},
        };
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
 *     u = Kp(b*r - y) + Ki*integral(r-y)dt + Kd*d/dt(c*r - y)
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
 * @param i_min Minimum continuous integrator value (∫e)
 * @param i_max Maximum continuous integrator value
 * @param Kbc   Back-calculation coefficient: ΔI += (Ts/Kbc)(u−u_unsat); T_t = Kbc/Ki
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
 * Control law (error-sum integrator, filtered derivative):
 * @f[
 *   I \leftarrow I + e,\quad
 *   D \leftarrow a_d D + b_d\,\Delta(c r - y),\quad
 *   u = K_p(b r - y) + K_{i,d} I + D
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
    T Ki{};  ///< Discrete integral gain Kᵢ T_s
    T Kd{};  ///< Continuous Kd (unused in the tick; retained for inspection)
    T d_a{}; ///< D ← d_a D + d_b Δx
    T d_b{};
    T u_min = -std::numeric_limits<T>::max();
    T u_max = std::numeric_limits<T>::max();
    T i_min = -std::numeric_limits<T>::max();
    T i_max = std::numeric_limits<T>::max();
    T Kbc = T{0}; ///< Back-calculation coefficient; ΔI += (u−u_unsat)/Kbc
    T b = T{1};
    T c = T{1};
    T Ts = T{0}; ///< Sample period baked at discretize (documentation)

    T    integral = T{0};        ///< Discrete error-sum state (∑e)
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
            // Auto does I += e before the output; store the pre-integration value.
            if (Ki != T{0}) {
                const T target = ((u_track - (Kp * ((b * r) - y)) - deriv) / Ki) - e;
                integral = damp::clamp(target, i_min, i_max);
            }
            return damp::clamp(u_track, u_min, u_max);
        }

        integral += e;
        integral = damp::clamp(integral, i_min, i_max);
        const T u_unsat = (Kp * ((b * r) - y)) + (Ki * integral) + deriv;
        const T u = damp::clamp(u_unsat, u_min, u_max);

        // Back-calculation: ΔI += (u − u_unsat)/Kbc (Ts already in discrete I)
        if (Kbc != T{0}) {
            integral += (u - u_unsat) / Kbc;
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
     * Winds the discrete integrator by `(u_sat - u_unsat) / Kbc`. No-op when
     * `Kbc == 0`. Sample time is already baked into the integrator units.
     */
    constexpr void back_calculate(T u_unsat, T u_sat) {
        if (Kbc == T{0} || u_unsat == u_sat) {
            return;
        }
        integral += (u_sat - u_unsat) / Kbc;
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
    T Ki{}; ///< Discrete integral gain Kᵢ T_s
    T u_min = -std::numeric_limits<T>::max();
    T u_max = std::numeric_limits<T>::max();
    T i_min = -std::numeric_limits<T>::max();
    T i_max = std::numeric_limits<T>::max();
    T Kbc = T{0};
    T b = T{1};
    T Ts = T{0};

    T integral = T{0}; ///< Discrete error-sum state (∑e)

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
            if (Ki != T{0}) {
                const T target = ((u_track - (Kp * ((b * r) - y))) / Ki) - e;
                integral = damp::clamp(target, i_min, i_max);
            }
            return damp::clamp(u_track, u_min, u_max);
        }

        integral += e;
        integral = damp::clamp(integral, i_min, i_max);
        const T u_unsat = (Kp * ((b * r) - y)) + (Ki * integral);
        const T u = damp::clamp(u_unsat, u_min, u_max);

        if (Kbc != T{0}) {
            integral += (u - u_unsat) / Kbc;
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
        if (Kbc == T{0} || u_unsat == u_sat) {
            return;
        }
        integral += (u_sat - u_unsat) / Kbc;
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
    T Kbc = T{0}; ///< Back-calculation: ΔI += (Ts/Kbc)(u−u_unsat); T_t = Kbc/Ki
    T b = T{1};
    T c = T{1};
    T Tf = T{0};

    T    integral = T{0}; ///< Continuous integrator state (∫e dt)
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
            return damp::clamp((Kp * ((b * r) - y)) + (Ki * integral), u_min, u_max);
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
            if (Ki != T{0}) {
                const T target = ((u_track - (Kp * ((b * r) - y)) - deriv) / Ki) - (e * Ts);
                integral = damp::clamp(target, i_min, i_max);
            }
            return damp::clamp(u_track, u_min, u_max);
        }

        integral += e * Ts;
        integral = damp::clamp(integral, i_min, i_max);
        const T u_unsat = (Kp * ((b * r) - y)) + (Ki * integral) + deriv;
        const T u = damp::clamp(u_unsat, u_min, u_max);

        if (Kbc != T{0}) {
            integral += (Ts / Kbc) * (u - u_unsat);
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
     * @brief Anti-windup hook: ΔI += (u_sat − u_unsat) * Ts / Kbc
     */
    constexpr void back_calculate(T u_unsat, T u_sat, T Ts) {
        if (Kbc == T{0} || u_unsat == u_sat) {
            return;
        }
        integral += Ts * ((u_sat - u_unsat) / Kbc);
        integral = damp::clamp(integral, i_min, i_max);
    }
};

} // namespace damp
