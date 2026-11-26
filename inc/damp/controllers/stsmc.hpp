// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file stsmc.hpp
 * @brief Super-twisting (second-order sliding-mode) control: design::stsmc(spec, Ts, opts...) →
 *        STSMCResult → STSMCController, continuous finite-time rejection of a Lipschitz matched
 *        disturbance
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/anti_windup.hpp"
#include "damp/controllers/detail/sliding_mode.hpp"
#include "damp/controllers/law.hpp" // IWYU pragma: export
#include "damp/design/detail/siso_input_limits.hpp"
#include "damp/design/input_limits.hpp" // IWYU pragma: export
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/toolbox/bounds.hpp"

namespace damp {

namespace design::detail {
using ::damp::detail::stsmc_gain_set;
} // namespace design::detail

namespace design {

/**
 * @struct STSMCGains
 * @brief The gains of a super-twisting law: k₁, k₂, the (r, y) form's surface slope λ, the
 *        generalized algorithm's linear gain and the smoothed sign's half-width
 *
 * design::stsmc's input when you tune k₁ and k₂ yourself, the result's gains, and what
 * STSMCController::gains and set_gains read and take:
 * `design::stsmc({.k1 = 2.0, .k2 = 2.2, .lambda = 5.0}, Ts)`.
 *
 * @tparam T Scalar type (deduced from k1 when the type is named)
 */
template<typename T = double>
struct STSMCGains {
    T                        k1{};     ///< Gain on |s|^½·sign(s) (> 0)
    damp::type_identity_t<T> k2{};     ///< Integral gain on sign(s) (> 0)
    damp::type_identity_t<T> lambda{}; ///< [1/s] Slope of the (r, y) form's surface σ = −(λe + ė) (≥ 0; 0 when you supply s yourself)
    damp::type_identity_t<T> k_lin{};  ///< Generalized-STA linear gain (≥ 0; 0: the classic algorithm)
    damp::type_identity_t<T> phi{};    ///< Half-width φ of the smoothed sign (≥ 0; 0: the true sign): exact for |s| ≥ φ, so |s| settles within about φ

    /** @brief The same gains in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr STSMCGains<U> as() const {
        return STSMCGains<U>{static_cast<U>(k1), static_cast<U>(k2), static_cast<U>(lambda), static_cast<U>(k_lin), static_cast<U>(phi)};
    }
};

/**
 * @struct STSMCSpec
 * @brief A super-twisting design from the bound L on the matched disturbance's rate, |ḋ| ≤ L
 *
 * design::stsmc's input when the gains come from Levant's and Moreno and Osorio's sufficient
 * choice, k₁ = 1.5√L, k₂ = 1.1L, both scaled by the margin:
 * `design::stsmc({.disturbance_bound = 2.0, .lambda = 5.0}, Ts)`.
 *
 * @tparam T Scalar type (deduced from disturbance_bound when the type is named)
 */
template<typename T = double>
struct STSMCSpec {
    T                        disturbance_bound{}; ///< L, the bound on |ḋ| (> 0)
    damp::type_identity_t<T> lambda{};            ///< [1/s] Slope of the (r, y) form's surface (≥ 0; 0 when you supply s yourself)
    damp::type_identity_t<T> k_lin{};             ///< Generalized-STA linear gain (≥ 0; 0: the classic algorithm)
    damp::type_identity_t<T> phi{};               ///< Half-width φ of the smoothed sign (≥ 0; 0: the true sign)
    damp::type_identity_t<T> gain_margin{T{1}};   ///< Robustness multiplier on both gains (≥ 1)
};

/**
 * @struct STSMCLimits
 * @brief The limits of a super-twisting law: on its command
 *
 * Unbounded by default. design::stsmc sets them from an InputLimits among its options, and the
 * law reads and changes them at run time with STSMCController::limits and set_limits.
 */
template<typename T = double>
struct STSMCLimits {
    Bounds<1, T> input{}; ///< The command u [lower, upper]: the integral state holds against it

    /** @brief The same limits in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr STSMCLimits<U> as() const {
        return STSMCLimits<U>{input.template as<U>()};
    }
};

/**
 * @struct STSMCResult
 * @brief A super-twisting design deployed at a sample period: the gains, the limits the law
 *        starts with, and T_s
 *
 * @see "Strict Lyapunov Functions for the Super-Twisting Algorithm" (Moreno & Osorio, 2012),
 *      IEEE TAC
 */
template<typename T = double>
struct [[nodiscard]] STSMCResult {
    STSMCGains<T>  gains{};                             ///< k₁, k₂, λ, k_lin, φ
    STSMCLimits<T> limits{};                            ///< The limits the law starts with
    T              Ts{T{0}};                            ///< [s] Sample period of the tick
    bool           success{false};                      ///< status == DesignStatus::ok
    DesignStatus   status{DesignStatus::not_converged}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The same design in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr STSMCResult<U> as() const {
        return STSMCResult<U>{gains.template as<U>(), limits.template as<U>(), static_cast<U>(Ts), success, status};
    }
};

namespace detail {

/** @brief The gains deployed at @p Ts with the InputLimits among @p opts: refused as design::stsmc documents. */
template<typename T, typename... Opts>
[[nodiscard]] constexpr STSMCResult<T> stsmc_deployed(const STSMCGains<T>& gains, T Ts, const Opts&... opts) {
    check_siso_input_limits<Opts...>();
    const STSMCLimits<T> limits{siso_input_limits_of<T>(opts...)};
    if (!damp::finite_positive(Ts)) {
        return refused<STSMCResult<T>>(DesignStatus::invalid_sample_time);
    }
    if (!limits.input.valid() || !stsmc_gain_set(gains.k1, gains.k2, gains.lambda, gains.k_lin, gains.phi, Ts)) {
        return refused<STSMCResult<T>>(DesignStatus::invalid_parameter);
    }
    STSMCResult<T> out{.gains = gains, .limits = limits, .Ts = Ts};
    set_status(out, DesignStatus::ok);
    return out;
}

} // namespace detail

/**
 * @brief Super-twisting gains from a bound on the disturbance's rate, deployed at @p Ts
 *
 * The super-twisting algorithm rejects a matched disturbance d(t) acting on the
 * sliding dynamics ṡ = u + d, provided the disturbance is Lipschitz with a known
 * bound on its rate, |ḋ| ≤ L. A standard sufficient gain choice (Levant /
 * Moreno-Osorio) is
 *
 *     k₁ = 1.5·√L,   k₂ = 1.1·L
 *
 * scaled by the spec's gain_margin ≥ 1 for robustness to model error. The
 * unique feature versus first-order SMC: the control is *continuous* (the
 * discontinuity is hidden under the integrator), so it rejects the disturbance
 * with finite-time convergence without the boundary-layer chattering of
 * `SMCController`.
 *
 * For a noisy measurement or a higher-order/under-modelled actuator (e.g. a
 * hydraulic swash-plate stage), set k_lin > 0 to use the generalized
 * super-twisting algorithm, which adds linear damping terms (`+ k_lin·s`) that
 * stiffen convergence and de-sensitize the loop to the `|s|^½` term amplifying
 * sensor noise — the usual cure when textbook super-twisting still chatters on
 * real hardware.
 *
 * @code
 * constexpr auto kSt = design::stsmc({.disturbance_bound = 2.0, .lambda = 5.0}, 1e-3,
 *                                    design::InputLimits{.lower = {-1.0}, .upper = {1.0}});
 * static constinit Law<kSt, float> st; // each tick: u = st.step(r, y);
 * @endcode
 *
 * @note Discretized with explicit Euler on the integral state. Keep gains matched
 *       to the sample rate; if `s` is differentiated from a noisy signal, filter
 *       it (or supply `s` directly) — noise on `s` is the dominant real-world
 *       chattering source.
 *
 * Refused with invalid_sample_time for a Ts that is not finite and positive, and with
 * invalid_parameter for L not finite and positive, a gain_margin below 1 or not finite, a λ,
 * k_lin or φ below 0 or not finite, or an InputLimits that is not a box. An option other than one
 * InputLimits on one input does not compile.
 *
 * @param spec L, λ, k_lin, φ and the margin
 * @param Ts   [s] The law's sample period (> 0)
 * @param opts InputLimits{.lower = {lo}, .upper = {hi}}
 * @return The design: the gains, the limits and Ts, for STSMCController or Law
 *
 * @see "Sliding order and sliding accuracy in sliding mode control" (Levant, 2003)
 * @see smc() for the first-order (boundary-layer) sliding-mode controller
 */
template<typename T, typename... Opts>
[[nodiscard]] constexpr STSMCResult<T> stsmc(const STSMCSpec<T>& spec, T Ts, const Opts&... opts) {
    const T    L = spec.disturbance_bound;
    const bool in_domain = damp::finite_positive(L) && damp::isfinite(spec.gain_margin) && spec.gain_margin >= T{1};
    const auto gains = in_domain ? STSMCGains<T>{
                                       .k1 = spec.gain_margin * static_cast<T>(1.5) * damp::sqrt(L),
                                       .k2 = spec.gain_margin * static_cast<T>(1.1) * L,
                                       .lambda = spec.lambda,
                                       .k_lin = spec.k_lin,
                                       .phi = spec.phi,
                                   }
                                 : STSMCGains<T>{}; // refused below: k₁ = k₂ = 0
    return detail::stsmc_deployed(gains, Ts, opts...);
}

/**
 * @brief Super-twisting gains you specify directly, deployed at @p Ts
 *
 * Bypasses the disturbance-bound formula for users who tune k₁ and k₂ themselves. Refused as the
 * spec's form is, and with invalid_parameter for k₁ or k₂ not finite and positive.
 *
 * @param gains k₁, k₂, λ, k_lin, φ
 * @param Ts    [s] The law's sample period (> 0)
 * @param opts  InputLimits{.lower = {lo}, .upper = {hi}}
 * @return The design: the gains, the limits and Ts, for STSMCController or Law
 */
template<typename T, typename... Opts>
[[nodiscard]] constexpr STSMCResult<T> stsmc(const STSMCGains<T>& gains, T Ts, const Opts&... opts) {
    return detail::stsmc_deployed(gains, Ts, opts...);
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
 * Two entry points per step, each with its handover twin:
 *  - `step_sliding(s)`, `track_sliding(u, s)` — canonical: you form the sliding variable s yourself.
 *  - `step(r, y)`, `track(u, r, y)` — convenience: builds σ = −(λ·e + ė) on e = r − y (ė by
 *    backward difference) — `SMCController`'s surface, signed so that σ̇ = u + d for ÿ = u + d.
 *
 * | Member | What it is |
 * | --- | --- |
 * | step(r, y), step_sliding(s) | one tick at the design's T_s |
 * | track(u, r, y), track_sliding(u, s) | in place of each while another source drives the plant |
 * | back_calculate(u_sat, u_ff) | a clamp after the law: v holds against it |
 * | gains(), set_gains(g) | k₁, k₂, λ, k_lin, φ in force (design::STSMCGains); change them at run time |
 * | limits(), set_limits(lim) | the command's limits (design::STSMCLimits); change them at run time |
 * | integral_state(), last_command(), sample_time(), valid(), reset() | v; the last command; T_s; whether the design ran; clear v and the rate's memory |
 *
 * The integral state is updated by explicit Euler. The output is the
 * sliding-mode control contribution; add any equivalent/feedforward control
 * separately. The command is clamped into limits().input, and v holds against that clamp
 * (conditional integration), as against one reported with back_calculate.
 *
 * While another source drives the plant, track(u, r, y) takes the place of step(r, y): the rate's
 * memory moves as the tick's, and v is preloaded so that step(r, y) at this r and y returns u
 * clamped into limits().input, u recorded as the command applied. track_sliding(u, s) does the
 * same for step_sliding(s).
 *
 * Example: reject a sinusoidal load disturbance on an integrator plant
 * @code
 * using namespace damp;
 * // |ḋ| ≤ 2, 1 kHz loop; supply s directly.
 * constexpr auto art = design::stsmc({.disturbance_bound = 2.0}, 1e-3);
 * static constinit STSMCController u_st{art.as<float>()}; // a failed design does not compile
 * // In the loop: float u = u_st.step_sliding(s);  // s = your sliding variable
 * @endcode
 *
 * @tparam T     Scalar type (default: float)
 * @tparam Gains Where the gains live: void, the default, for members (settable), or the constants
 *               of a `constexpr` design's Law (detail::ConstantSTSMCGains)
 * @see "Strict Lyapunov Functions for the Super-Twisting Algorithm" (Moreno &
 *      Osorio, 2012), IEEE TAC, https://doi.org/10.1109/TAC.2012.2186179
 */
template<typename T = float, typename Gains = void>
class STSMCController {
    using GainSet = damp::conditional_t<damp::is_same_v<Gains, void>, detail::STSMCGainSet<T>, Gains>; ///< The gains' storage
    static constexpr bool kHeld = !GainSet::constant;                                                  ///< The gains are members: set_gains changes them

    template<typename, typename>
    friend class STSMCController; // cross-precision converting ctor

public:
    constexpr STSMCController() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr STSMCController(const design::STSMCResult<U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : STSMCController(result.template as<T>()) {}

    /**
     * @brief From the design result, which carries the gains, the limits and the sample time.
     *
     * A failed design, a T_s that is not finite and positive, a gain outside its domain or a limit
     * that is not a box is ill-formed under constant evaluation; at run time it leaves the
     * controller invalid (valid() false, output 0 from both entry points).
     */
    constexpr STSMCController(const design::STSMCResult<T>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : u_lim_(result.limits.input) {
        if (!result.success && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(result.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
        if (!damp::finite_positive(result.Ts) && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(design::DesignStatus::invalid_sample_time);
        }
        valid_ = result.success && result.limits.input.valid() && take(result.gains, result.Ts);
        if (!valid_ && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(design::DesignStatus::invalid_parameter); // a gain outside its domain, or a limit that is not a box
        }
    }

    /** @brief The same law, with its live state, in another scalar (or off a Law's constant gains, onto members). */
    template<typename U, typename G>
        requires kHeld
    constexpr explicit STSMCController(const STSMCController<U, G>& other)
        : u_lim_(other.u_lim_.template as<T>()),
          v_(static_cast<T>(other.v_)),
          v_step_(static_cast<T>(other.v_step_)),
          u_cmd_(static_cast<T>(other.u_cmd_)),
          e_prev_(static_cast<T>(other.e_prev_)),
          first_(other.first_) {
        valid_ = other.valid_ && take(other.gains().template as<T>(), static_cast<T>(other.sample_time()));
    }

    /**
     * @brief Super-twisting control from the sliding variable @p s (canonical form).
     *
     * Named distinctly from step(r, y) because both would otherwise be two-scalar overloads.
     * The integral state advances by one explicit-Euler step of the design's Ts. The command
     * is clamped into limits().input, and a step of the integral state that pushes the command
     * further past a limit it is on is taken back (conditional integration): v holds against the
     * clamp where it ran on while s could not reach 0. The step enters the next command, so v
     * holds at most one step, k₂·Ts, past the limit, as it does for a clamp reported with
     * back_calculate.
     *
     * @param s Sliding variable
     */
    [[nodiscard]] constexpr T step_sliding(T s) {
        if (!valid_) {
            return T{0};
        }
        const T abs_s = damp::abs(s);
        // sign(s), optionally smoothed inside |s| < φ to trade a little accuracy for less
        // numerical chatter near s = 0; exact ±1 outside, so the full k₂ acts beyond the layer.
        const T sign_s = damp::smooth_sign(s, k_.phi);
        const T sqrt_s = damp::sqrt(abs_s);

        const T u_unsat = -(k_.k1 * proportional(s, sign_s, sqrt_s)) + v_;
        const T u = u_lim_.saturate(u_unsat);

        // φ₂ = sign(s) + 3·k_lin·|s|^½·sign(s) + 2·k_lin²·s  (= classic sign(s) when k_lin = 0,
        // keeping the k₁ = 1.5√L / k₂ = 1.1L formulas valid). Explicit-Euler integral.
        const T phi2 = sign_s + (T{3} * k_.k_lin * sqrt_s * sign_s) + (T{2} * k_.k_lin * k_.k_lin * s);
        v_step_ = -k_.k2 * phi2 * k_.Ts; // explicit-Euler integral
        v_ += v_step_;
        detail::hold_against_clamp(v_, v_step_, u - u_unsat);

        u_cmd_ = u;
        return u;
    }

    /**
     * @brief Convenience: build σ = −(λ·e + ė), e = r − y (ė via backward difference),
     * then apply STA to σ.
     * @note For a unit-gain plant ÿ = u + d the surface on the output error
     *       y − r obeys σ̇ = u + (d − λė − r̈), the canonical sign step_sliding
     *       assumes — so positive tracking error gives a positive command, as
     *       for PID / SMC. Scale u for a non-unit input gain, or use
     *       `step_sliding(s)` with a surface you sign yourself.
     */
    [[nodiscard]] constexpr T step(T r, T y) {
        if (!valid_) {
            return T{0};
        }
        return step_sliding(surface(r - y));
    }

    /**
     * @brief In place of step_sliding(@p s) while another source drives the plant, @p u the input applied
     *
     * Preloads v so that step_sliding(@p s) at this same s returns @p u clamped into
     * limits().input: v = sat(u) + k₁·(|s|^½·sign(s) + k_lin·s). @p u is recorded as the command
     * applied, so a back_calculate of it finds no excess.
     */
    constexpr void track_sliding(T u, T s) {
        if (!valid_) {
            return;
        }
        const T sign_s = damp::smooth_sign(s, k_.phi);
        v_ = u_lim_.saturate(u) + (k_.k1 * proportional(s, sign_s, damp::sqrt(damp::abs(s))));
        v_step_ = T{0};
        u_cmd_ = u;
    }

    /**
     * @brief In place of step(@p r, @p y) while another source drives the plant, @p u the input applied
     *
     * The rate's memory moves as the tick's, then v is preloaded so that step(@p r, @p y) at this
     * same r and y, where ė = 0, returns @p u clamped into limits().input, as track_sliding does.
     */
    constexpr void track(T u, T r, T y) {
        if (!valid_) {
            return;
        }
        static_cast<void>(surface(r - y));
        track_sliding(u, -(k_.lambda * (r - y))); // the next step's surface at this r and y: ė = 0
    }

    /**
     * @brief Record the input the plant received, @p u_sat, when a stage after this law clamped
     *        its command, and the feedforward @p u_ff added to the command before that clamp
     *
     * The integral state gives back what the last tick added, as far as the clamp cut the
     * command off: u_sat − (u + u_ff), u the command the last tick returned.
     */
    constexpr void back_calculate(T u_sat, T u_ff = T{0}) { detail::hold_against_clamp(v_, v_step_, u_sat - (u_cmd_ + u_ff)); }

    /** @brief The command the last tick returned, or the input the last track recorded. */
    [[nodiscard]] constexpr T last_command() const { return u_cmd_; }

    /** @brief Current value of the integral state v (the disturbance estimate, −d̂). */
    [[nodiscard]] constexpr T integral_state() const { return v_; }

    /** @brief False when built from a design it refuses: both entry points return 0. */
    [[nodiscard]] constexpr bool valid() const { return valid_; }

    /** @brief Clear v, the rate's memory and the last command; keeps the gains and the limits. */
    constexpr void reset() {
        v_ = T{0};
        v_step_ = T{0};
        u_cmd_ = T{0};
        e_prev_ = T{0};
        first_ = true;
    }

    /** @brief The gains in force: the design's, or the last set_gains. */
    [[nodiscard]] constexpr design::STSMCGains<T> gains() const { return {k_.k1, k_.k2, k_.lambda, k_.k_lin, k_.phi}; }

    /**
     * @brief Change the gains at run time; the next tick uses them
     *
     * Keeps v, so the command does not jump by the integral. Refused (false, nothing changed) for
     * k₁ or k₂ not finite and positive, or a λ, k_lin or φ below 0 or not finite.
     */
    [[nodiscard]] constexpr bool set_gains(const design::STSMCGains<T>& g)
        requires kHeld
    {
        return take(g, k_.Ts);
    }

    /** @brief The limits in force: the design's, or the last set_limits. */
    [[nodiscard]] constexpr design::STSMCLimits<T> limits() const { return {u_lim_}; }

    /** @brief Change the limits at run time; the next tick clamps into them. Refused (false, nothing changed) when an end is NaN or ±Inf, or the lower end is above the upper one. */
    [[nodiscard]] constexpr bool set_limits(const design::STSMCLimits<T>& lim) {
        if (!lim.input.valid()) {
            return false;
        }
        u_lim_ = lim.input;
        return true;
    }

    /** @brief [s] The tick's sample period T_s. */
    [[nodiscard]] constexpr T sample_time() const { return k_.Ts; }

private:
    /** @brief Store @p g at the period @p Ts when it is in its domain; on a constexpr design's Law, whose gains are the type's, only check it. */
    [[nodiscard]] constexpr bool take(const design::STSMCGains<T>& g, T Ts) {
        const auto set = detail::stsmc_gain_set(g.k1, g.k2, g.lambda, g.k_lin, g.phi, Ts);
        if (!set) {
            return false;
        }
        if constexpr (kHeld) {
            k_ = *set;
        }
        return true;
    }

    /** @brief φ₁ = |s|^½·sign(s) + k_lin·s, the term k₁ multiplies. */
    [[nodiscard]] constexpr T proportional(T s, T sign_s, T sqrt_s) const { return (sqrt_s * sign_s) + (k_.k_lin * s); }

    /** @brief σ = −(λe + ė) on the error @p e, ė its backward difference, 0 on the first tick (no derivative kick). */
    [[nodiscard]] constexpr T surface(T e) {
        if (first_) {
            e_prev_ = e; // seed: no first-tick derivative kick into the surface
            first_ = false;
        }
        const T de = (e - e_prev_) * k_.inv_ts;
        e_prev_ = e;
        return -((k_.lambda * e) + de);
    }

    [[no_unique_address]] GainSet k_{};          ///< The gains, T_s and 1/T_s: members, or no storage on a constexpr design's Law
    Bounds<1, T>                  u_lim_{};      ///< The command's limits: the integral state holds against them
    T                             v_{T{0}};      ///< The integral state v
    T                             v_step_{T{0}}; ///< What the last tick added to v: what a clamp takes back
    T                             u_cmd_{T{0}};  ///< The command the last tick returned, or the input track recorded
    T                             e_prev_{T{0}}; ///< Previous error, for the (r, y) form's backward difference
    bool                          first_{true};  ///< Seed e_prev_ on the first (r, y) tick
    bool                          valid_{false}; ///< From the design's success flag; gates both entry points
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<typename T>
STSMCController(const design::STSMCResult<T>&) -> STSMCController<T>;

namespace detail {

/** @brief Whether @p R is a super-twisting design result. */
template<typename R>
inline constexpr bool is_stsmc_result = false;

/** @brief design::stsmc's result. */
template<typename U>
inline constexpr bool is_stsmc_result<design::STSMCResult<U>> = true;

/** @brief A super-twisting design (design::stsmc). */
template<const auto& Res>
concept stsmc_design = is_stsmc_result<design_t<Res>>;

} // namespace detail

/** @brief Law<kRes, T> of a super-twisting design initialized at start-up: the STSMCController it runs as. A runtime assigned replaces it. */
template<const auto& Res, typename T>
    requires detail::stsmc_design<Res>
struct Law<Res, T> : STSMCController<T> {
    using Runtime = STSMCController<T>; ///< The runtime the design runs as

    /** @brief The design's law. */
    constexpr Law() : Runtime{Res} {}

    using Runtime::operator=;
};

/**
 * @brief Law<kRes, T> of a `constexpr` super-twisting design: the STSMCController with its gains
 *        constants of the type (detail::ConstantSTSMCGains)
 *
 * The gains and T_s are in flash and fold into the tick; the limits stay members, settable at
 * run time; there is no set_gains. A law retuned at run time is
 * `Law<kRes, T>::Runtime law{kRes.as<T>()};`. A failed design, or one whose gains or T_s its
 * runtime refuses, does not compile wherever the Law is declared.
 */
template<const auto& Res, typename T>
    requires(detail::stsmc_design<Res> && detail::constant_design<Res>)
struct Law<Res, T> : STSMCController<T, detail::ConstantSTSMCGains<Res, T>> {
    using Runtime = STSMCController<T>; ///< The runtime with members: the type of a member a redesign assigns
    using Base = STSMCController<T, detail::ConstantSTSMCGains<Res, T>>;

    static_assert(detail::constant_design_runs(Res.success, Res.status, Res.Ts));
    static_assert(detail::ConstantSTSMCGains<Res, T>::ok, "Law<kRes, T>: a gain of the super-twisting design is outside its domain at its Ts");

    /** @brief The design's law: its gains the type's, its limits members. */
    constexpr Law() : Base{Res} {}

    /** @brief There is nothing to replace: the gains are in the type. */
    template<typename Other>
        requires(!damp::is_same_v<damp::remove_cvref_t<Other>, Law>)
    constexpr Law& operator=(const Other& other) {
        static_assert(sizeof(Other) == 0, "This Law is a constexpr super-twisting design's: its gains are constants of the type, so there is nothing to replace. A member a redesign assigns, or a law retuned with set_gains, is Law<kRes, T>::Runtime law{kRes.as<T>()};");
        static_cast<void>(other);
        return *this;
    }
};

} // namespace damp
