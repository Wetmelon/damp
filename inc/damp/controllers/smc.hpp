// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once


/**
 * @file smc.hpp
 * @brief Sliding-mode control design and runtime
 */

#include "damp/backend.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/systems/transfer_function.hpp"

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
    T            lambda{};                            ///< Sliding-surface slope λ in s = λ·e + ė [1/s]. Larger = faster but noisier.
    T            k{};                                 ///< Switching gain (must exceed the lumped disturbance bound) [units of ÿ].
    T            b0{};                                ///< Nominal plant input gain, ÿ = f + b0·u (> 0).
    T            phi{};                               ///< Boundary-layer thickness (0 = hard sign()); > 0 softens the relay to sat(s/φ).
    bool         success{false};                      ///< status == DesignStatus::ok
    DesignStatus status{DesignStatus::not_converged}; ///< ok, or why the design failed
    T            Ts{0};                               ///< [s] Sample time of the law: 0 continuous (analysis), > 0 the runtime's tick

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief Re-cast the parameters to a different scalar type (e.g. double → float). */
    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return SMCResult<U>{static_cast<U>(lambda), static_cast<U>(k), static_cast<U>(b0), static_cast<U>(phi), success, status, static_cast<U>(Ts)};
    }

    /**
     * @brief Boundary-layer linearization C(s) = (k / (b₀ φ)) (s + λ)
     *
     * Inside |s| < φ the relay is sat(s/φ), so u = (k/b₀φ)(λ e + ė) with e = r − y.
     * Improper (relative degree −1). Requires the design's φ > 0; otherwise a zero TF.
     */
    [[nodiscard]] constexpr TransferFunction<2, 1, T> to_tf() const {
        TransferFunction<2, 1, T> tf{};
        if (!success || !(phi > T{0}) || !(b0 > T{0})) {
            return tf;
        }
        const T g = k / (b0 * phi);
        tf.num = {g, g * lambda}; // g s + g λ
        tf.den = {T{1}};
        return tf;
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
 * @param k      Switching gain (> the worst-case lumped disturbance on ÿ) [units of ÿ].
 * @param b0     Input gain of the relative-degree-2 plant ÿ = f + b0·u (> 0).
 * @param Ts     [s] Sample time of the runtime's backward-difference rate (0 = continuous, analysis only).
 * @param phi    Boundary-layer thickness (≥ 0; 0 = hard sign()).
 * @return SMCResult carrying the parameters; refused with invalid_parameter for λ, k, b0 ≤ 0 or φ < 0,
 *         and with invalid_sample_time for a negative or non-finite Ts.
 */
template<typename T = double>
[[nodiscard]] constexpr SMCResult<T> smc(T lambda, T k, T b0, T Ts = T{0}, T phi = T{0}) {
    // reject λ/k ≤ 0 or b0 ≤ 0 rather than producing a divide-by-zero controller
    if (!(damp::finite_positive(lambda) && damp::finite_positive(k) && damp::finite_positive(b0) && damp::finite_non_negative(phi))) {
        return ::damp::detail::refused<SMCResult<T>>(DesignStatus::invalid_parameter);
    }
    if (!damp::finite_non_negative(Ts)) {
        return ::damp::detail::refused<SMCResult<T>>(DesignStatus::invalid_sample_time);
    }
    return SMCResult<T>{
        .lambda = lambda,
        .k = k,
        .b0 = b0,
        .phi = phi,
        .success = true,
        .status = DesignStatus::ok,
        .Ts = Ts,
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
 * sliding. The surface @f$ (d/dt + \lambda)\,e @f$ is the relative-degree-2 case of
 * Slotine's @f$ (d/dt + \lambda)^{n-1} e @f$: for the plant @f$ \ddot y = f + b_0 u @f$
 * with input gain @f$ b_0 > 0 @f$ it obeys @f$ \dot s = -b_0 u + \ldots @f$, so the
 * relay (bang-bang) term that drives `s` to zero is
 *
 * @f[ u = +\frac{k}{b_0}\,\operatorname{sign}(s) @f]
 *
 * — positive error gives a positive command, the same sign as a PID on @f$ r - y @f$.
 * On a first-order plant @f$ \dot y = f + b_0 u @f$ the surface itself contains
 * @f$ -b_0 u @f$ and the relay settles off the reference; use @ref ADRCController
 * (NX = 1) or a PID there.
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
 *  - `b0` — your plant's input gain (ÿ = f + b0·u, > 0); if unsure,
 *    set 1 and fold the scaling into `k`.
 *
 * @note `ė` is estimated by backward difference `(e − e_prev)/Ts`, so a noisy `y`
 *       feeds noise straight into the surface — filter the measurement (or its
 *       derivative) on a real sensor. Ts comes from the design result: design
 *       at the loop's period.
 *
 * Example — regulate a position with a boundary layer to avoid actuator buzz:
 * @code
 * using namespace damp;
 * // λ, k, b0, Ts = 1 ms, phi = 0.05 boundary layer to avoid chatter
 * constexpr auto gains = design::smc(20.0f, 5.0f, 1.0f, 1e-3f, 0.05f);
 * static constinit SMCController ctrl{gains};
 * // In the loop:
 * float u = ctrl.control(r, y);
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
    T    Ts_{};         ///< [s] Sample time, from the design result.
    T    phi_{};        ///< Boundary-layer thickness, from the design result.
    bool first_{true};  ///< Seed error_prev on the first control() tick (no derivative kick).
    bool valid_{false}; ///< From the design's success flag; gates control().

public:
    constexpr SMCController() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr SMCController(const design::SMCResult<U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : SMCController(result.template as<T>()) {}

    /**
     * @brief From the design result, which carries the sample time and boundary layer.
     *
     * A failed or continuous (Ts = 0) design is ill-formed under constant evaluation; at run
     * time it leaves the controller invalid (valid() false, output 0).
     */
    constexpr SMCController(const design::SMCResult<T>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : lambda(result.lambda), k(result.k), b0(result.b0), Ts_(result.Ts), phi_(result.phi), valid_(result.success && result.Ts > T{0}) {
        if ((!result.success || !(result.Ts > T{0})) && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(result.success ? design::DesignStatus::invalid_sample_time : result.status); // constant-evaluated (constinit / constexpr): a failed or continuous design does not compile
        }
    }

    template<typename U>
    constexpr SMCController(const SMCController<U>& other)
        : lambda(static_cast<T>(other.lambda)), k(static_cast<T>(other.k)), b0(static_cast<T>(other.b0)), error_prev(static_cast<T>(other.error_prev)), Ts_(static_cast<T>(other.Ts_)), phi_(static_cast<T>(other.phi_)), first_(other.first_), valid_(other.valid_) {}

    /**
     * @brief Run one control step: u = (k/b0)·sign(λ·e + ė), e = r − y.
     *
     * The rate is the backward difference `ė = (e − e_prev)/Ts` over the design's Ts; the
     * design's boundary layer φ > 0 softens the relay to sat(s/φ) (smoother command, larger
     * steady-state error).
     *
     * @param r   Reference (setpoint).
     * @param y   Measurement (plant output).
     * @return Control command `u`.
     */
    [[nodiscard]] constexpr T control(T r, T y) {
        // A failed or continuous design (no b0, no sampled rate): emit no correction, hold the state.
        if (!valid_) {
            return T{0};
        }

        const T error = r - y;
        if (first_) {
            error_prev = error; // seed: this tick's rate is 0, no derivative kick
            first_ = false;
        }
        const T dot_e = (error - error_prev) / Ts_;
        error_prev = error;

        // Sliding surface: λ weights the error against its rate.
        const T s = (lambda * error) + dot_e;

        // phi <= 0: hard relay. phi > 0: saturation boundary layer sat(s/phi) — a
        // linear ramp through |s| < phi, full ±1 switching gain restored beyond it.
        if (phi_ <= T{0}) {
            return (k / b0) * static_cast<T>(damp::sgn(s));
        }
        return (k / b0) * damp::clamp(s / phi_, T{-1}, T{1});
    }

    [[nodiscard]] constexpr bool valid() const { return valid_; }

    constexpr void reset() {
        error_prev = T{};
        first_ = true;
    }
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<typename T>
SMCController(const design::SMCResult<T>&) -> SMCController<T>;

} // namespace damp