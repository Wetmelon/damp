// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file smc.hpp
 * @brief Sliding-mode control: design::smc(gains, Ts, opts...) → SMCResult → SMCController
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/sliding_mode.hpp"
#include "damp/controllers/law.hpp" // IWYU pragma: export
#include "damp/design/detail/siso_input_limits.hpp"
#include "damp/design/input_limits.hpp" // IWYU pragma: export
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/toolbox/bounds.hpp"

namespace damp {

namespace design::detail {
using ::damp::detail::smc_gain_set;
} // namespace design::detail

namespace design {

/**
 * @struct SMCGains
 * @brief The gains of a first-order sliding-mode law: the surface slope λ, the switching gain k,
 *        the plant's input gain b₀ and the boundary layer φ
 *
 * design::smc's input, the result's gains, and what SMCController::gains and set_gains read and
 * take. SMC has no gain formula: you pick λ and k (SMCController says how). The caller names
 * each number, `design::smc({.lambda = 20.0, .k = 5.0, .b0 = 1.0, .phi = 0.05}, Ts)`.
 *
 * @tparam T Scalar type (deduced from lambda when the type is named)
 */
template<typename T = double>
struct SMCGains {
    T                        lambda{}; ///< [1/s] Sliding-surface slope λ in s = λe + ė (> 0): larger is faster and noisier
    damp::type_identity_t<T> k{};      ///< Switching gain (> 0), above the lumped disturbance's bound on ÿ
    damp::type_identity_t<T> b0{};     ///< Input gain of the plant ÿ = f + b₀u (> 0)
    damp::type_identity_t<T> phi{};    ///< Boundary-layer thickness φ (≥ 0; 0: the hard relay sign(s), > 0: sat(s/φ))

    /** @brief The same gains in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr SMCGains<U> as() const {
        return SMCGains<U>{static_cast<U>(lambda), static_cast<U>(k), static_cast<U>(b0), static_cast<U>(phi)};
    }
};

/**
 * @struct SMCLimits
 * @brief The limits of an SMC law: on its command
 *
 * Unbounded by default. design::smc sets them from an InputLimits among its options, and the law
 * reads and changes them at run time with SMCController::limits and set_limits.
 */
template<typename T = double>
struct SMCLimits {
    Bounds<1, T> input{}; ///< The command u [lower, upper]

    /** @brief The same limits in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr SMCLimits<U> as() const {
        return SMCLimits<U>{input.template as<U>()};
    }
};

/**
 * @struct SMCResult
 * @brief A sliding-mode design deployed at a sample period: the gains, the limits the law starts
 *        with, and T_s
 *
 * The runtime differentiates the error by backward difference over T_s, and computes k/b₀, 1/T_s
 * and 1/φ when the gains are set. Use `.as<U>()` to deploy in another scalar.
 *
 * @tparam T Scalar type (default: double)
 */
template<typename T = double>
struct [[nodiscard]] SMCResult {
    SMCGains<T>  gains{};                             ///< λ, k, b₀, φ
    SMCLimits<T> limits{};                            ///< The limits the law starts with
    T            Ts{T{0}};                            ///< [s] Sample period of the tick
    bool         success{false};                      ///< status == DesignStatus::ok
    DesignStatus status{DesignStatus::not_converged}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The same design in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr SMCResult<U> as() const {
        return SMCResult<U>{gains.template as<U>(), limits.template as<U>(), static_cast<U>(Ts), success, status};
    }

    /**
     * @brief Boundary-layer linearization C(s) = (k/(b₀φ))(s + λ), for analysis
     *
     * Inside |s| < φ the relay is sat(s/φ), so u = (k/b₀φ)(λe + ė) with e = r − y. Improper
     * (relative degree −1). The zero TF for a failed design or one without a boundary layer.
     */
    [[nodiscard]] constexpr TransferFunction<2, 1, T> to_tf() const {
        TransferFunction<2, 1, T> tf{};
        if (!success || !damp::finite_positive(gains.phi) || !damp::finite_positive(gains.b0)) {
            return tf;
        }
        const T g = gains.k / (gains.b0 * gains.phi);
        tf.num = {g, g * gains.lambda}; // g s + g λ
        tf.den = {T{1}};
        return tf;
    }
};

/**
 * @brief A first-order sliding-mode law from the gains you pick, deployed at @p Ts
 *
 * There is no synthesis: SMC is tuned by choosing λ and k (SMCController says how). The design
 * checks them and deploys them at the loop's period, with the limits among the options.
 *
 * @code
 * constexpr auto kSmc = design::smc({.lambda = 20.0, .k = 5.0, .b0 = 1.0, .phi = 0.05}, 1e-3,
 *                                   design::InputLimits{.lower = {-2.0}, .upper = {2.0}});
 * static constinit Law<kSmc, float> smc; // each tick: u = smc.step(r, y);
 * @endcode
 *
 * Refused with invalid_sample_time for a Ts that is not finite and positive, and with
 * invalid_parameter for λ, k or b₀ not finite and positive, φ below 0 or not finite, a
 * coefficient (k/b₀, 1/T_s, 1/φ) that overflows, or an InputLimits that is not a box. An option
 * other than one InputLimits on one input does not compile.
 *
 * @param gains λ, k, b₀, φ
 * @param Ts    [s] The law's sample period (> 0)
 * @param opts  InputLimits{.lower = {lo}, .upper = {hi}}
 * @return The design: the gains, the limits and Ts, for SMCController or Law
 */
template<typename T, typename... Opts>
[[nodiscard]] constexpr SMCResult<T> smc(const SMCGains<T>& gains, T Ts, const Opts&... opts) {
    detail::check_siso_input_limits<Opts...>();
    const SMCLimits<T> limits{detail::siso_input_limits_of<T>(opts...)};
    if (!damp::finite_positive(Ts)) {
        return detail::refused<SMCResult<T>>(DesignStatus::invalid_sample_time);
    }
    if (!limits.input.valid() || !detail::smc_gain_set(gains.lambda, gains.k, gains.b0, gains.phi, Ts)) {
        return detail::refused<SMCResult<T>>(DesignStatus::invalid_parameter);
    }
    SMCResult<T> out{.gains = gains, .limits = limits, .Ts = Ts};
    detail::set_status(out, DesignStatus::ok);
    return out;
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
 * | Member | What it is |
 * | --- | --- |
 * | step(r, y) | one tick at the design's T_s |
 * | track(u, r, y) | in place of step while another source drives the plant |
 * | back_calculate(u_sat, u_ff) | a clamp after the law: a no-op, nothing winds up |
 * | gains(), set_gains(g) | λ, k, b₀, φ in force (design::SMCGains); change them at run time |
 * | limits(), set_limits(lim) | the command's limits (design::SMCLimits); change them at run time |
 * | last_command(), sample_time(), valid(), reset() | the last command; T_s; whether the design ran; clear the rate's memory |
 *
 * The command is clamped into limits().input. The law keeps only the last error, for the
 * backward difference ė, so a clamp after it winds nothing up (back_calculate is a no-op), and
 * there is nothing to preload: track(u, r, y) keeps the last error current, so the step after it
 * at this r and y has ė = 0, no derivative kick, and records u as the command applied.
 *
 * @note `ė` is estimated by backward difference `(e − e_prev)/Ts`, so a noisy `y`
 *       feeds noise straight into the surface — filter the measurement (or its
 *       derivative) on a real sensor. Ts comes from the design result: design
 *       at the loop's period.
 *
 * Example — regulate a position with a boundary layer to avoid actuator buzz:
 * @code
 * using namespace damp;
 * constexpr auto kSmc = design::smc({.lambda = 20.0, .k = 5.0, .b0 = 1.0, .phi = 0.05}, 1e-3);
 * static constinit SMCController ctrl{kSmc.as<float>()};
 * // In the loop:
 * float u = ctrl.step(r, y);
 * @endcode
 *
 * @tparam T     Scalar type (default: float)
 * @tparam Gains Where the gains live: void, the default, for members (settable), or the constants
 *               of a `constexpr` design's Law (detail::ConstantSMCGains)
 *
 * @see STSMCController for the chatter-free (second-order) sibling.
 * @see Utkin, "Variable Structure Systems with Sliding Modes," IEEE TAC 22(2),
 *      1977, https://doi.org/10.1109/TAC.1977.1101446
 * @see Slotine & Li, *Applied Nonlinear Control*, Prentice Hall, 1991 (ch. 7) —
 *      the standard textbook treatment of the surface/boundary-layer design.
 */
template<typename T = float, typename Gains = void>
class SMCController {
    using GainSet = damp::conditional_t<damp::is_same_v<Gains, void>, detail::SMCGainSet<T>, Gains>; ///< The gains' storage
    static constexpr bool kHeld = !GainSet::constant;                                                ///< The gains are members: set_gains changes them

    template<typename, typename>
    friend class SMCController; // cross-precision converting ctor

public:
    constexpr SMCController() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr SMCController(const design::SMCResult<U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : SMCController(result.template as<T>()) {}

    /**
     * @brief From the design result, which carries the gains, the limits and the sample time.
     *
     * A failed design, a T_s that is not finite and positive, a gain outside its domain or a limit
     * that is not a box is ill-formed under constant evaluation; at run time it leaves the
     * controller invalid (valid() false, output 0).
     */
    constexpr SMCController(const design::SMCResult<T>& result) // NOLINT(google-explicit-constructor): assign-from-design
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
    constexpr explicit SMCController(const SMCController<U, G>& other)
        : u_lim_(other.u_lim_.template as<T>()),
          e_prev_(static_cast<T>(other.e_prev_)),
          u_cmd_(static_cast<T>(other.u_cmd_)),
          first_(other.first_) {
        valid_ = other.valid_ && take(other.gains().template as<T>(), static_cast<T>(other.sample_time()));
    }

    /**
     * @brief One tick: u = (k/b₀)·sign(λe + ė), e = r − y, clamped into limits().input
     *
     * The rate is the backward difference `ė = (e − e_prev)/Ts` over the design's Ts, zero on the
     * first tick (no derivative kick); the boundary layer φ > 0 softens the relay to sat(s/φ)
     * (smoother command, larger steady-state error).
     *
     * @param r Reference (setpoint)
     * @param y Measurement (plant output)
     * @return The command u, 0 when the law is invalid
     */
    [[nodiscard]] constexpr T step(T r, T y) {
        if (!valid_) {
            return T{0};
        }
        const T error = r - y;
        const T s = (k_.lambda * error) + rate(error); // λ weights the error against its rate
        // φ = 0: the hard relay. φ > 0: the saturation sat(s/φ), a linear ramp through |s| < φ.
        const T relay = (k_.inv_phi == T{0}) ? static_cast<T>(damp::sgn(s)) : damp::clamp(s * k_.inv_phi, T{-1}, T{1});
        u_cmd_ = u_lim_.saturate(k_.relay * relay);
        return u_cmd_;
    }

    /**
     * @brief In place of step(@p r, @p y) while another source drives the plant, @p u the input applied
     *
     * The relay has no memory to preload: the last error is kept current, so the step after it at
     * this r and y has ė = 0 and no derivative kick, and @p u is recorded as the command applied.
     */
    constexpr void track(T u, T r, T y) {
        if (!valid_) {
            return;
        }
        static_cast<void>(rate(r - y));
        u_cmd_ = u;
    }

    /** @brief A clamp downstream: the law keeps only the last error, which no clamp winds up, so nothing to take back. */
    constexpr void back_calculate(T u_sat, T u_ff = T{0}) {
        static_cast<void>(u_sat);
        static_cast<void>(u_ff);
    }

    /** @brief The command the last step() returned, or the input the last track() recorded. */
    [[nodiscard]] constexpr T last_command() const { return u_cmd_; }

    /** @brief False when built from a design it refuses: step() returns 0. */
    [[nodiscard]] constexpr bool valid() const { return valid_; }

    /** @brief Clear the rate's memory and the last command; keeps the gains and the limits. */
    constexpr void reset() {
        e_prev_ = T{0};
        u_cmd_ = T{0};
        first_ = true;
    }

    /** @brief The gains in force: the design's, or the last set_gains. */
    [[nodiscard]] constexpr design::SMCGains<T> gains() const { return {k_.lambda, k_.k, k_.b0, k_.phi}; }

    /**
     * @brief Change the gains at run time; the next tick uses them
     *
     * Recomputes k/b₀ and 1/φ and keeps the last error. Refused (false, nothing changed) for λ, k
     * or b₀ not finite and positive, φ below 0 or not finite, or a coefficient that overflows.
     */
    [[nodiscard]] constexpr bool set_gains(const design::SMCGains<T>& g)
        requires kHeld
    {
        return take(g, k_.Ts);
    }

    /** @brief The limits in force: the design's, or the last set_limits. */
    [[nodiscard]] constexpr design::SMCLimits<T> limits() const { return {u_lim_}; }

    /** @brief Change the limits at run time; the next tick clamps into them. Refused (false, nothing changed) when an end is NaN or ±Inf, or the lower end is above the upper one. */
    [[nodiscard]] constexpr bool set_limits(const design::SMCLimits<T>& lim) {
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
    [[nodiscard]] constexpr bool take(const design::SMCGains<T>& g, T Ts) {
        const auto set = detail::smc_gain_set(g.lambda, g.k, g.b0, g.phi, Ts);
        if (!set) {
            return false;
        }
        if constexpr (kHeld) {
            k_ = *set;
        }
        return true;
    }

    /** @brief The backward difference of the error @p e, 0 on the first tick; keeps @p e for the next. */
    [[nodiscard]] constexpr T rate(T e) {
        if (first_) {
            e_prev_ = e; // seed: this tick's rate is 0, no derivative kick
            first_ = false;
        }
        const T de = (e - e_prev_) * k_.inv_ts;
        e_prev_ = e;
        return de;
    }

    [[no_unique_address]] GainSet k_{};          ///< The gains, T_s and the coefficients at it: members, or no storage on a constexpr design's Law
    Bounds<1, T>                  u_lim_{};      ///< The command's limits
    T                             e_prev_{};     ///< Previous error, for the backward-difference rate ė
    T                             u_cmd_{};      ///< The command the last step() returned, or the input track() recorded
    bool                          first_{true};  ///< Seed e_prev_ on the first tick (no derivative kick)
    bool                          valid_{false}; ///< From the design's success flag; gates step()
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<typename T>
SMCController(const design::SMCResult<T>&) -> SMCController<T>;

namespace detail {

/** @brief Whether @p R is an SMC design result. */
template<typename R>
inline constexpr bool is_smc_result = false;

/** @brief design::smc's result. */
template<typename U>
inline constexpr bool is_smc_result<design::SMCResult<U>> = true;

/** @brief An SMC design (design::smc). */
template<const auto& Res>
concept smc_design = is_smc_result<design_t<Res>>;

} // namespace detail

/** @brief Law<kRes, T> of an SMC design initialized at start-up: the SMCController it runs as. A runtime assigned replaces it. */
template<const auto& Res, typename T>
    requires detail::smc_design<Res>
struct Law<Res, T> : SMCController<T> {
    using Runtime = SMCController<T>; ///< The runtime the design runs as

    /** @brief The design's law. */
    constexpr Law() : Runtime{Res} {}

    using Runtime::operator=;
};

/**
 * @brief Law<kRes, T> of a `constexpr` SMC design: the SMCController with its gains constants of
 *        the type (detail::ConstantSMCGains)
 *
 * The gains, T_s and k/b₀, 1/T_s, 1/φ are in flash and fold into the tick; the limits stay
 * members, settable at run time; there is no set_gains. A law retuned at run time is
 * `Law<kRes, T>::Runtime law{kRes.as<T>()};`. A failed design, or one whose gains or T_s its
 * runtime refuses, does not compile wherever the Law is declared.
 */
template<const auto& Res, typename T>
    requires(detail::smc_design<Res> && detail::constant_design<Res>)
struct Law<Res, T> : SMCController<T, detail::ConstantSMCGains<Res, T>> {
    using Runtime = SMCController<T>; ///< The runtime with members: the type of a member a redesign assigns
    using Base = SMCController<T, detail::ConstantSMCGains<Res, T>>;

    static_assert(detail::constant_design_runs(Res.success, Res.status, Res.Ts));
    static_assert(detail::ConstantSMCGains<Res, T>::ok, "Law<kRes, T>: a gain of the SMC design is outside its domain, or a coefficient overflows at its Ts");

    /** @brief The design's law: its gains the type's, its limits members. */
    constexpr Law() : Base{Res} {}

    /** @brief There is nothing to replace: the gains are in the type. */
    template<typename Other>
        requires(!damp::is_same_v<damp::remove_cvref_t<Other>, Law>)
    constexpr Law& operator=(const Other& other) {
        static_assert(sizeof(Other) == 0, "This Law is a constexpr SMC design's: its gains are constants of the type, so there is nothing to replace. A member a redesign assigns, or a law retuned with set_gains, is Law<kRes, T>::Runtime law{kRes.as<T>()};");
        static_cast<void>(other);
        return *this;
    }
};

} // namespace damp
