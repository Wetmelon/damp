// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once


/**
 * @file stsmc.hpp
 * @brief Super-twisting (second-order sliding-mode) controller — continuous
 *        finite-time rejection of a Lipschitz matched disturbance.
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/anti_windup.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/toolbox/bounds.hpp"
namespace damp {

namespace design {

/**
 * @struct STSMCResult
 * @brief Super-twisting (second-order sliding-mode) controller design result.
 *
 * Holds the two super-twisting gains plus the optional shaping parameters. Use
 * .as\<U\>() to convert for embedded deployment.
 *
 * @see "Strict Lyapunov Functions for the Super-Twisting Algorithm" (Moreno &
 *      Osorio, 2012), IEEE TAC
 */
template<typename T = double>
struct [[nodiscard]] STSMCResult {
    T            k1{};                                ///< Continuous gain on |s|^½·sign(s)
    T            k2{};                                ///< Integral gain on sign(s)
    T            lambda{};                            ///< Sliding-surface slope for the (r, y) convenience form (σ = −(λe + ė), e = r − y)
    T            k_lin{};                             ///< Generalized-STA linear gain (0 = classic super-twisting)
    T            phi{};                               ///< Boundary-layer half-width φ of the smoothed sign (0 = true sign)
    bool         success{false};                      ///< status == DesignStatus::ok
    DesignStatus status{DesignStatus::not_converged}; ///< ok, or why the design failed
    T            Ts{0};                               ///< [s] Sample time of the law: 0 continuous (analysis), > 0 the runtime's tick
    Bounds<1, T> u_lim{};                             ///< Limits on the command (default unbounded): set here, or on the runtime

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return STSMCResult<U>{
            static_cast<U>(k1), static_cast<U>(k2), static_cast<U>(lambda), static_cast<U>(k_lin),
            static_cast<U>(phi), success, status, static_cast<U>(Ts), u_lim.template as<U>()
        };
    }
};

/**
 * @brief Synthesize super-twisting gains from a disturbance-derivative bound.
 *
 * The super-twisting algorithm rejects a matched disturbance d(t) acting on the
 * sliding dynamics ṡ = u + d, provided the disturbance is Lipschitz with a known
 * bound on its rate, |ḋ| ≤ L. A standard sufficient gain choice (Levant /
 * Moreno-Osorio) is
 *
 *     k₁ = 1.5·√L,   k₂ = 1.1·L
 *
 * scaled by an optional @p gain_margin ≥ 1 for robustness to model error. The
 * unique feature versus first-order SMC: the control is *continuous* (the
 * discontinuity is hidden under the integrator), so it rejects the disturbance
 * with finite-time convergence without the boundary-layer chattering of
 * `SMCController`.
 *
 * For a noisy measurement or a higher-order/under-modelled actuator (e.g. a
 * hydraulic swash-plate stage), set @p k_lin > 0 to use the generalized
 * super-twisting algorithm, which adds linear damping terms (`+ k_lin·s`) that
 * stiffen convergence and de-sensitize the loop to the `|s|^½` term amplifying
 * sensor noise — the usual cure when textbook super-twisting still chatters on
 * real hardware.
 *
 * @note Discretized with explicit Euler on the integral state. Keep gains matched
 *       to the sample rate; if `s` is differentiated from a noisy signal, filter
 *       it (or supply `s` directly) — noise on `s` is the dominant real-world
 *       chattering source.
 *
 * @tparam T Scalar type
 * @param disturbance_bound L, the bound on |ḋ| (must be > 0)
 * @param Ts                [s] Sample time of the runtime's explicit-Euler integral (0 = continuous, analysis only)
 * @param lambda            Sliding-surface slope for the (r, y) form (≥ 0; 0 if you supply s directly)
 * @param k_lin             Generalized-STA linear gain (≥ 0; 0 = classic super-twisting)
 * @param phi               Boundary-layer half-width φ of the smoothed sign (≥ 0; 0 = true sign):
 *                          the sign is exact for |s| ≥ φ, so |s| settles within about φ
 * @param gain_margin       Robustness multiplier on both gains (≥ 1)
 * @return STSMCResult with computed k₁, k₂ and the shaping parameters; refused with
 *         invalid_parameter for an out-of-range input, invalid_sample_time for a negative Ts
 *
 * @see "Sliding order and sliding accuracy in sliding mode control" (Levant, 2003)
 * @see smc() for the first-order (boundary-layer) sliding-mode controller
 */
template<typename T = double>
[[nodiscard]] constexpr STSMCResult<T> stsmc(
    T disturbance_bound,
    T Ts = T{0},
    T lambda = T{0},
    T k_lin = T{0},
    T phi = T{0},
    T gain_margin = T{1}
) {
    if (!(damp::finite_positive(disturbance_bound) && damp::isfinite(gain_margin) && gain_margin >= T{1}
          && damp::finite_non_negative(k_lin) && damp::finite_non_negative(phi) && damp::finite_non_negative(lambda))) {
        return ::damp::detail::refused<STSMCResult<T>>(DesignStatus::invalid_parameter);
    }
    if (!damp::finite_non_negative(Ts)) {
        return ::damp::detail::refused<STSMCResult<T>>(DesignStatus::invalid_sample_time);
    }
    const T L = disturbance_bound;
    return STSMCResult<T>{
        .k1 = gain_margin * static_cast<T>(1.5) * damp::sqrt(L),
        .k2 = gain_margin * static_cast<T>(1.1) * L,
        .lambda = lambda,
        .k_lin = k_lin,
        .phi = phi,
        .success = true,
        .status = DesignStatus::ok,
        .Ts = Ts,
    };
}

/**
 * @brief Super-twisting controller from gains you specify directly.
 *
 * Bypasses the disturbance-bound formula for users who tune k₁/k₂ themselves.
 * Separate name from @ref stsmc so a two-arg call resolves cleanly against that
 * API's defaulted `(L, Ts)` form.
 *
 * @param k1 Continuous gain on |s|^½·sign(s) (> 0)
 * @param k2 Integral gain on sign(s) (> 0)
 * @param Ts      [s] Sample time of the runtime (0 = continuous, analysis only)
 * @param lambda  Sliding-surface slope for the (r, y) form (≥ 0)
 * @param k_lin   Generalized-STA linear gain (≥ 0)
 * @param phi     Boundary-layer half-width φ of the smoothed sign (≥ 0)
 */
template<typename T = double>
[[nodiscard]] constexpr STSMCResult<T>
stsmc_gains(T k1, T k2, T Ts = T{0}, T lambda = T{0}, T k_lin = T{0}, T phi = T{0}) {
    if (!(damp::finite_positive(k1) && damp::finite_positive(k2) && damp::finite_non_negative(k_lin)
          && damp::finite_non_negative(phi) && damp::finite_non_negative(lambda))) {
        return ::damp::detail::refused<STSMCResult<T>>(DesignStatus::invalid_parameter);
    }
    if (!damp::finite_non_negative(Ts)) {
        return ::damp::detail::refused<STSMCResult<T>>(DesignStatus::invalid_sample_time);
    }
    return STSMCResult<T>{
        .k1 = k1,
        .k2 = k2,
        .lambda = lambda,
        .k_lin = k_lin,
        .phi = phi,
        .success = true,
        .status = DesignStatus::ok,
        .Ts = Ts,
    };
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Super-twisting controller (second-order sliding mode).
 *
 * Continuous-control sliding mode for systems of relative degree 1 in the
 * sliding variable (ṡ = u + d). The control law is
 *
 *     u   = −k₁·(|s|^½·sign(s) + k_lin·s) + v
 *     v̇   = −k₂·(sign(s) + 3·k_lin·|s|^½·sign(s) + 2·k_lin²·s)
 *
 * With k_lin = 0 this is the classic super-twisting algorithm (Levant):
 * `u = −k₁·|s|^½·sign(s) + v`, `v̇ = −k₂·sign(s)`. The integral state v is what
 * makes u continuous — finite-time rejection of a Lipschitz matched disturbance
 * with none of the first-order boundary-layer chattering. With k_lin > 0 the
 * generalized form adds linear damping for noisy / higher-order plants.
 *
 * Two entry points per step:
 *  - `control_sliding(s)` — canonical: you form the sliding variable s yourself.
 *  - `control(r, y)`      — convenience: builds σ = −(λ·e + ė) on e = r − y (ė by
 *                          backward difference) — `SMCController`'s surface,
 *                          signed so that σ̇ = u + d for ÿ = u + d.
 *
 * The integral state is updated by explicit Euler. The output is the
 * sliding-mode control contribution; add any equivalent/feedforward control
 * separately.
 *
 * Example: reject a sinusoidal load disturbance on an integrator plant
 * @code
 * using namespace damp;
 * // |ḋ| ≤ 2, 1 kHz loop; supply s directly.
 * constexpr auto art = design::stsmc(2.0, 1e-3);
 * static constinit STSMCController u_st{art.as<float>()}; // a failed or continuous design does not compile
 * // In the loop: float u = u_st.control_sliding(s);  // s = your sliding variable
 * @endcode
 *
 * @tparam T Scalar type (default: float)
 * @see "Strict Lyapunov Functions for the Super-Twisting Algorithm" (Moreno &
 *      Osorio, 2012), IEEE TAC, https://doi.org/10.1109/TAC.2012.2186179
 */
template<typename T = float>
class STSMCController {
public:
    Bounds<1, T> u_lim{}; ///< Limits on the command (the result's; default unbounded): the integral state holds against them

    constexpr STSMCController() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr STSMCController(const design::STSMCResult<U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : STSMCController(result.template as<T>()) {}

    /**
     * @brief From the design result, which carries the sample time.
     *
     * A failed or continuous (Ts = 0) design, or one whose u_lim is not a box, is ill-formed
     * under constant evaluation; at run time it leaves the controller invalid (valid() false,
     * output 0 from both entry points).
     */
    constexpr STSMCController(const design::STSMCResult<T>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : u_lim(result.u_lim),
          k1_(result.k1),
          k2_(result.k2),
          lambda_(result.lambda),
          k_lin_(result.k_lin),
          phi_(result.phi),
          Ts_(result.Ts),
          valid_(result.success && result.Ts > T{0} && result.u_lim.valid()) {
        if (result.success && !result.u_lim.valid() && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(design::DesignStatus::invalid_parameter); // a limit that is not finite, or lower above upper
        }
        if ((!result.success || !(result.Ts > T{0})) && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(result.success ? design::DesignStatus::invalid_sample_time : result.status); // constant-evaluated (constinit / constexpr): a failed or continuous design does not compile
        }
    }

    template<typename U>
    constexpr STSMCController(const STSMCController<U>& other)
        : u_lim(other.u_lim.template as<T>()),
          k1_(static_cast<T>(other.k1_)),
          k2_(static_cast<T>(other.k2_)),
          lambda_(static_cast<T>(other.lambda_)),
          k_lin_(static_cast<T>(other.k_lin_)),
          phi_(static_cast<T>(other.phi_)),
          Ts_(static_cast<T>(other.Ts_)),
          v_(static_cast<T>(other.v_)),
          v_step_(static_cast<T>(other.v_step_)),
          u_cmd_(static_cast<T>(other.u_cmd_)),
          e_prev_(static_cast<T>(other.e_prev_)),
          first_(other.first_),
          valid_(other.valid_) {}

    /**
     * @brief Super-twisting control from the sliding variable @p s (canonical form).
     * Named distinctly from control(r, y) because both would otherwise be
     * two-scalar overloads.
     * The integral state advances by one explicit-Euler step of the design's Ts. The command
     * is clamped into u_lim, and a step of the integral state that pushes the command further
     * past a limit it is on is taken back (conditional integration): v holds against the clamp
     * where it ran on while s could not reach 0.
     * @param s  Sliding variable
     */
    [[nodiscard]] constexpr T control_sliding(T s) {
        if (!valid_) {
            return T{0};
        }
        const T abs_s = damp::abs(s);
        // sign(s), optionally smoothed inside |s| < φ to trade a little accuracy for less
        // numerical chatter near s = 0; exact ±1 outside, so the full k₂ acts beyond the layer.
        const T sign_s = damp::smooth_sign(s, phi_);
        const T sqrt_s = damp::sqrt(abs_s);

        const T phi1 = (sqrt_s * sign_s) + (k_lin_ * s);
        const T u_unsat = (-k1_ * phi1) + v_;
        const T u = u_lim.saturate(u_unsat);

        // φ₂ = sign(s) + 3·k_lin·|s|^½·sign(s) + 2·k_lin²·s  (= classic sign(s) when k_lin = 0,
        // keeping the k₁ = 1.5√L / k₂ = 1.1L formulas valid). Explicit-Euler integral.
        const T phi2 = sign_s + (T{3} * k_lin_ * sqrt_s * sign_s) + (T{2} * k_lin_ * k_lin_ * s);
        v_step_ = -k2_ * phi2 * Ts_; // explicit-Euler integral
        v_ += v_step_;
        detail::hold_against_clamp(v_, v_step_, u - u_unsat);

        u_cmd_ = u;
        return u;
    }

    /**
     * @brief Record the input the plant received, @p u_sat, when a stage after this law clamped
     *        its command, and the feedforward @p u_ff added to the command before that clamp
     *
     * The integral state gives back what the last tick added, as far as the clamp cut the
     * command off: u_sat − (u + u_ff), u the command the last tick returned.
     */
    constexpr void back_calculate(T u_sat, T u_ff = T{0}) { detail::hold_against_clamp(v_, v_step_, u_sat - (u_cmd_ + u_ff)); }

    /** @brief The command the last tick returned. */
    [[nodiscard]] constexpr T last_command() const { return u_cmd_; }

    /**
     * @brief Convenience: build σ = −(λ·e + ė), e = r − y (ė via backward difference),
     * then apply STA to σ.
     * @note For a unit-gain plant ÿ = u + d the surface on the output error
     *       y − r obeys σ̇ = u + (d − λė − r̈), the canonical sign control_sliding
     *       assumes — so positive tracking error gives a positive command, as
     *       for PID / SMC. Scale u for a non-unit input gain, or use
     *       `control_sliding(s)` with a surface you sign yourself.
     */
    [[nodiscard]] constexpr T control(T r, T y) {
        if (!valid_) {
            return T{0};
        }
        const T e = r - y;
        if (first_) {
            e_prev_ = e; // seed: no first-tick derivative kick into the surface
            first_ = false;
        }
        const T de = (e - e_prev_) / Ts_;
        e_prev_ = e;
        return control_sliding(-((lambda_ * e) + de));
    }

    /** @brief Current value of the integral state v (the disturbance estimate, −d̂). */
    [[nodiscard]] constexpr T integral_state() const { return v_; }

    [[nodiscard]] constexpr bool valid() const { return valid_; }

    constexpr void reset() {
        v_ = T{0};
        v_step_ = T{0};
        u_cmd_ = T{0};
        e_prev_ = T{0};
        first_ = true;
    }

private:
    template<typename>
    friend class STSMCController; // cross-precision converting ctor

    T    k1_{T{0}};
    T    k2_{T{0}};
    T    lambda_{T{0}};
    T    k_lin_{T{0}};
    T    phi_{T{0}}; ///< Boundary-layer half-width of the smoothed sign
    T    Ts_{T{0}};  ///< [s] Sample time, from the design result
    T    v_{T{0}};
    T    v_step_{T{0}}; ///< What the last tick added to v: what a clamp takes back
    T    u_cmd_{T{0}};  ///< The command the last tick returned
    T    e_prev_{T{0}};
    bool first_{true};
    bool valid_{false};
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<typename T>
STSMCController(const design::STSMCResult<T>&) -> STSMCController<T>;

} // namespace damp
