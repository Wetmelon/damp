// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file smc.hpp
 * @brief Sliding-mode control design and runtime
 */

#include "damp/math/math.hpp"

namespace damp {

namespace design {

/**
 * @struct SMCResult
 * @brief Tuning parameters for a first-order sliding-mode controller.
 *
 * Just the three numbers @ref SMCController needs. SMC has no gain *formula* the
 * way LQR or pole-placement do — you pick the sliding-surface slope @p lambda and
 * the switching gain @p k directly, so this struct is a plain parameter bundle
 * rather than the output of a synthesis step. Use `.as\<U\>()` to down-cast for an
 * embedded build (e.g. `double` design → `float` on an MCU).
 *
 * @tparam T Scalar type (default: double)
 */
template<typename T = double>
struct SMCResult {
    T    lambda{};       ///< Sliding-surface slope λ in s = λ·e + ė [1/s]. Larger = faster but noisier.
    T    k{};            ///< Switching gain (must exceed the disturbance/uncertainty bound) [output units].
    T    b0{};           ///< Control effectiveness ṡ = b0·u + … ; the nominal plant input gain.
    bool success{false}; ///< true if λ, k, b0 are all > 0 (b0 = 0 would divide-by-zero at runtime).

    /// Re-cast the parameters to a different scalar type (e.g. double → float).
    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return SMCResult<U>{static_cast<U>(lambda), static_cast<U>(k), static_cast<U>(b0), success};
    }
};

/**
 * @brief Bundle hand-picked SMC parameters into an @ref SMCResult.
 *
 * Convenience factory so call sites read `design::smc(λ, k, b0)` symmetrically
 * with the other controllers' design functions. There is no computation: SMC is
 * tuned by choosing @p lambda and @p k directly (see @ref SMCController for how
 * to pick them).
 *
 * @param lambda Sliding-surface slope λ [1/s] (> 0).
 * @param k      Switching gain (> the worst-case lumped disturbance).
 * @param b0     Control effectiveness / input gain (> 0).
 * @return SMCResult carrying the three parameters.
 */
template<typename T = double>
[[nodiscard]] constexpr SMCResult<T> smc(T lambda, T k, T b0) {
    // reject λ/k ≤ 0 or b0 ≤ 0 rather than producing a divide-by-zero controller
    if (!(lambda > T{0} && k > T{0} && b0 > T{0})) {
        return SMCResult<T>{};
    }
    return SMCResult<T>{
        .lambda = lambda,
        .k = k,
        .b0 = b0,
        .success = true,
    };
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief First-order sliding-mode controller (SMC) for a SISO plant.
 *
 * Sliding-mode control is a discontinuous nonlinear law: it forces the tracking
 * error onto a chosen line in the error/error-rate plane (the sliding surface)
 * and keeps it there under bounded disturbances and model error (finite-time
 * reachability of the surface when the switching gain bounds the lumped
 * disturbance).
 *
 * How it works, step by step. With error `e = r − y` and its rate `ė`, define
 * the sliding surface
 *
 * @f[ s = \lambda\,e + \dot e @f]
 *
 * `s = 0` is a line through the origin with slope `−λ`; on it the error decays as
 * `e(t) = e(0)\,\mathrm{e}^{-\lambda t}`, so λ sets the closed-loop speed once
 * sliding. The control drives `s` to zero with a relay (bang-bang) term:
 *
 * @f[ u = -\frac{k}{b_0}\,\operatorname{sign}(s) @f]
 *
 * As long as the switching gain `k` exceeds the worst-case lumped disturbance,
 * `s` reaches 0 in finite time and stays there. The price of the hard `sign()`
 * is chattering — buzzing as `u` slams between ±k/b₀ every sample. The
 * @p phi boundary layer fixes that: inside a band `|s| < phi` the relay is
 * replaced by the saturation `sat(s/phi)` (a linear ramp through the layer, ±1
 * outside, so full switching gain is restored beyond it), trading a little
 * steady-state accuracy for a continuous command. For *continuous* control with no such
 * trade-off, see the super-twisting controller (@ref STSMCController),
 * which hides the discontinuity under an integrator.
 *
 * Picking the knobs (Arduino®-friendly):
 *  - `λ` — start near your desired bandwidth in rad/s; bigger is faster but
 *    amplifies measurement noise through `ė`.
 *  - `k` — increase until disturbances are rejected; if it buzzes, that's
 *    chattering — add a `phi` boundary layer rather than dropping `k`.
 *  - `b0` — your plant's input gain (how much `ṡ` moves per unit `u`); if unsure,
 *    set 1 and fold the scaling into `k`.
 *
 * @note `ė` is estimated by backward difference `(e − e_prev)/Ts`, so a noisy `y`
 *       feeds noise straight into the surface — filter the measurement (or its
 *       derivative) on a real sensor. `Ts` is a per-call argument, so
 *       the same controller works at any/varying loop rate.
 *
 * Example — regulate a position with a boundary layer to avoid actuator buzz:
 * @code
 * using namespace damp;
 * constexpr auto gains = design::smc(20.0f, 5.0f, 1.0f); // λ, k, b0
 * SMCController<float> ctrl(gains);
 * // In the loop: Ts = 1 ms, phi = 0.05 boundary layer to avoid chatter.
 * float u = ctrl.control(r, y, 1e-3f, 0.05f);
 * @endcode
 *
 * @tparam T Scalar type (default: float)
 *
 * @see STSMCController for the chatter-free (second-order) sibling.
 * @see Utkin, "Variable Structure Systems with Sliding Modes," IEEE TAC 22(2),
 *      1977, https://doi.org/10.1109/TAC.1977.1101446
 * @see Slotine & Li, *Applied Nonlinear Control*, Prentice Hall, 1991 (ch. 7) —
 *      the standard textbook treatment of the surface/boundary-layer design.
 */
template<typename T = float>
class SMCController {
    template<typename>
    friend class SMCController; // cross-precision converting ctor

    T lambda{}; ///< Sliding-surface slope λ [1/s].
    T k{};      ///< Switching gain (> disturbance bound).
    T b0{};     ///< Control effectiveness / input gain.

    T    error_prev{};  ///< Previous error, for the backward-difference rate ė.
    T    Ts_{};         ///< Stored sample time for the 2-arg control(r, y) form [s].
    T    phi_{};        ///< Stored boundary-layer thickness for the 2-arg form.
    bool first_{true};  ///< Seed error_prev on the first control() tick (no derivative kick).
    bool valid_{false}; ///< From the design's success flag; gates control().

public:
    constexpr SMCController() = default;

    constexpr SMCController(const design::SMCResult<T>& result, T Ts = T{0}, T phi = T{0})
        : lambda(result.lambda), k(result.k), b0(result.b0), Ts_(Ts), phi_(phi), valid_(result.success) {}

    template<typename U>
    constexpr SMCController(const SMCController<U>& other)
        : lambda(static_cast<T>(other.lambda)), k(static_cast<T>(other.k)), b0(static_cast<T>(other.b0)), error_prev(static_cast<T>(other.error_prev)), Ts_(static_cast<T>(other.Ts_)), phi_(static_cast<T>(other.phi_)), first_(other.first_), valid_(other.valid_) {}

    /**
     * @brief Run one control step: u = −(k/b0)·sign(λ·e + ė).
     *
     * @param r   Reference (setpoint).
     * @param y   Measurement (plant output).
     * @param Ts  Sample time [s] since the last call; used for the backward-
     *            difference rate `ė = (e − e_prev)/Ts`. Pass your actual loop
     *            period — it need not be constant.
     * @param phi Boundary-layer thickness (default 0 = hard `sign()`). Set
     *            `phi > 0` to soften the relay and suppress chattering; bigger
     *            `phi` = smoother command but larger steady-state error.
     *
     * @return Control command `u`.
     */
    [[nodiscard]] constexpr T control(T r, T y, T Ts, T phi = T{0}) {
        // A non-positive Ts (uninitialized / stalled dt) would divide-by-zero in the
        // backward-difference rate; an invalid design has a zero/garbage b0. Either
        // way emit no correction and hold the state.
        if (!valid_ || Ts <= T{0}) {
            return T{0};
        }

        const T error = r - y;
        if (first_) {
            error_prev = error; // seed: this tick's rate is 0, no derivative kick
            first_ = false;
        }
        const T dot_e = (error - error_prev) / Ts;
        error_prev = error;

        // Sliding surface: λ weights the error against its rate.
        const T s = (lambda * error) + dot_e;

        // phi <= 0: hard relay. phi > 0: saturation boundary layer sat(s/phi) — a
        // linear ramp through |s| < phi, full ±1 switching gain restored beyond it.
        if (phi <= T{0}) {
            return -(k / b0) * static_cast<T>(damp::sgn(s));
        }
        return -(k / b0) * damp::clamp(s / phi, T{-1}, T{1});
    }

    /**
     * @brief 2-arg fixed-rate form using the stored sample time and boundary layer.
     *
     * The SisoController concept shape: construct with Ts (and optionally phi)
     * and call control(r, y) each tick. With no stored Ts it emits zero and
     * holds state, exactly like the 4-arg form with Ts = 0.
     */
    [[nodiscard]] constexpr T control(T r, T y) { return control(r, y, Ts_, phi_); }

    [[nodiscard]] constexpr bool valid() const { return valid_; }

    constexpr void reset() {
        error_prev = T{};
        first_ = true;
    }
};

} // namespace damp