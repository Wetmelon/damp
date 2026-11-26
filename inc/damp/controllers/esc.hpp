// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file esc.hpp
 * @brief Extremum-seeking control (ESC) — model-free real-time optimization of a
 *        measured objective by dither-and-demodulate gradient estimation.
 *
 * ESC drives an input û toward the extremum of an *unknown* input→objective map
 * J(u) using only the measured J — no model, no derivative. It adds a small
 * sinusoidal dither, correlates the resulting objective ripple with that dither
 * to estimate the local gradient, and integrates the gradient to climb (or
 * descend) the map. It is exactly the principled form of "perturb and observe":
 *
 *     u    = û + a·sin(ωt)                 // perturbed input
 *     J    = objective(u)                  // measured (plant + map)
 *     ξ    = HPF(J) · sin(ωt − φ)          // demodulated gradient estimate
 *     û̇    = ±k·LPF(ξ)                      // integrate (+ maximize, − minimize)
 *
 * The high-pass removes the operating-point DC so only the dither-induced ripple
 * remains; multiplying by the dither and averaging (the low-pass) yields a signal
 * whose mean is ∝ ∂J/∂u·a/2. Frequency separation makes it work: the dither ω is
 * fast, the high-pass corner sits below ω, and the integrator (gain k) is slow.
 *
 * Targets: PV/wind MPPT under changing irradiance, thermal/drivetrain efficiency
 * optimization, source seeking, online auto-tuning. Allocation-free, constexpr.
 *
 * @see Y. Tan, D. Nešić, I. Mareels, "On non-local stability properties of
 *      extremum seeking control," Automatica, 2006; K. B. Ariyur & M. Krstić,
 *      "Real-Time Optimization by Extremum-Seeking Control," Wiley, 2003.
 * @see Y. Tan et al., "Extremum seeking control for discrete-time systems," IEEE
 *      TAC, 2002, https://doi.org/10.1109/9.983370
 */


#include "damp/backend.hpp"
#include "damp/controllers/law.hpp" // IWYU pragma: export
#include "damp/design/detail/siso_input_limits.hpp"
#include "damp/design/input_limits.hpp" // IWYU pragma: export
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"      // damp::sin/cos, damp::clamp, numbers
#include "damp/toolbox/bounds.hpp" // Bounds

namespace damp {

/** @brief Whether ESC climbs to a maximum or descends to a minimum of the objective. */
enum class ExtremumType {
    Maximize,
    Minimize
};

namespace design {

/**
 * @brief The extremum-seeking design input: dither, gain and filter corners; Ts is passed beside it
 *
 * The design input of @ref esc and @ref esc_mppt, its frequencies in rad/s as every design's
 * are. The caller names every number:
 * `design::esc({.dither_amplitude = 0.05, .dither_omega = 2π·40, .gain = 60.0, .lp_omega = 2π·8}, Ts)`.
 * The design turns it into an @ref ESCResult (filter coefficients), the runtime realization; the
 * two are distinct types. The limits are an InputLimits among the design's options.
 *
 * dither_amplitude, dither_omega and gain have no usable default: left at 0 the design is
 * refused.
 *
 * @tparam T Scalar type
 */
template<typename T = double>
struct ESCSpec {
    T            dither_amplitude{};           ///< a, perturbation amplitude in input units (> 0)
    T            dither_omega{};               ///< [rad/s] Perturbation frequency ω (> 0, at least 8 samples per period: ω·Ts ≤ π/4)
    T            gain{};                       ///< k, adaptation (integrator) gain (> 0)
    ExtremumType type{ExtremumType::Maximize}; ///< Maximize (climb) or Minimize (descend)
    T            hp_omega{T{0}};               ///< [rad/s] High-pass corner, below the dither (≥ 0; 0 = dither_omega/5)
    T            lp_omega{T{0}};               ///< [rad/s] Gradient low-pass corner (≥ 0; 0 = no low-pass)
    T            u_init{T{0}};                 ///< Initial operating point û₀
    T            demod_phase{T{0}};            ///< Demodulation phase offset φ [rad] for plant lag
};

/**
 * @struct ESCLimits
 * @brief The limits of an extremum-seeking law: on the input it applies, û + dither, and on its
 *        operating point û
 *
 * Unbounded by default. design::esc sets `input` from an InputLimits among its options, and the
 * law reads and changes both at run time with ExtremumSeekingController::limits and set_limits.
 * û stays inside both, so with `input` alone it stays inside the input's limits, as the applied
 * input does; `operating_point` narrows it further, e.g. by the dither's amplitude on each side
 * so the dither is never clipped. Where the two do not meet, the input's limits win.
 */
template<typename T = double>
struct ESCLimits {
    Bounds<1, T> input{};           ///< The applied input û + a·sin φ [lower, upper]
    Bounds<1, T> operating_point{}; ///< The operating point û, inside the input's limits as well

    /** @brief The same limits in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr ESCLimits<U> as() const {
        return ESCLimits<U>{input.template as<U>(), operating_point.template as<U>()};
    }

    /** @brief Both are boxes: every end finite and no lower end above its upper one. */
    [[nodiscard]] constexpr bool valid() const { return input.valid() && operating_point.valid(); }
};

/**
 * @brief Extremum-seeking design result: the discrete realization the runtime runs
 *
 * Built by @ref esc from an @ref ESCSpec (Hz → rad/s, filter corners → coefficients).
 *
 * @tparam T Scalar type
 */
template<typename T = double>
struct [[nodiscard]] ESCResult {
    T            dither_amplitude{};                  ///< a, perturbation amplitude (input units)
    T            dither_omega{};                      ///< ω, perturbation frequency [rad/s]
    T            gain{};                              ///< k, integrator (adaptation) gain
    T            hp_alpha{};                          ///< High-pass coefficient (objective DC removal)
    T            lp_alpha{T{1}};                      ///< Gradient low-pass coefficient (1 = no LPF)
    T            demod_phase{};                       ///< φ, demodulation phase offset [rad] (plant lag)
    T            direction{T{1}};                     ///< +1 maximize, −1 minimize
    T            u_init{};                            ///< Initial operating point û₀
    ESCLimits<T> limits{};                            ///< The limits the law starts with
    T            Ts{};                                ///< Sample time [s]
    bool         success{false};                      ///< status == DesignStatus::ok
    DesignStatus status{DesignStatus::not_converged}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr ESCResult<damp::remove_const_t<U>> as() const {
        using out_t = damp::remove_const_t<U>;
        return ESCResult<out_t>{
            .dither_amplitude = static_cast<out_t>(dither_amplitude),
            .dither_omega = static_cast<out_t>(dither_omega),
            .gain = static_cast<out_t>(gain),
            .hp_alpha = static_cast<out_t>(hp_alpha),
            .lp_alpha = static_cast<out_t>(lp_alpha),
            .demod_phase = static_cast<out_t>(demod_phase),
            .direction = static_cast<out_t>(direction),
            .u_init = static_cast<out_t>(u_init),
            .limits = limits.template as<out_t>(),
            .Ts = static_cast<out_t>(Ts),
            .success = success,
            .status = status,
        };
    }
};

namespace detail {

/** @brief First-order discrete low-pass coefficient for corner @p wc [rad/s] at @p Ts; 1 (pass-through) for wc = 0. */
template<typename T>
[[nodiscard]] constexpr T esc_lpf_alpha(T wc, T Ts) {
    if (!(wc > T{0})) {
        return T{1};
    }
    return (wc * Ts) / (T{1} + (wc * Ts));
}
} // namespace detail

/**
 * @brief Synthesize an extremum-seeking controller.
 *
 * The scalar type follows @p Ts, so a braced spec needs no template argument.
 *
 * The dither must be well resolved by the sample rate, not merely below Nyquist: ESC needs the
 * perturbation cleanly sampled to demodulate it, so the design requires ≥ 8 samples per dither
 * period (ω·Ts ≤ π/4). Below that the gradient estimate degrades and the scale separation
 * (HPF < ω < Nyquist) collapses.
 *
 * @code
 * constexpr auto kMppt = design::esc({.dither_amplitude = 0.02, .dither_omega = 628.0, .gain = 200.0, .u_init = 0.5}, 1e-4,
 *                                    design::InputLimits{.lower = {0.2}, .upper = {0.8}});
 * static constinit Law<kMppt, float> mppt; // each tick: u = mppt.step(J);
 * @endcode
 *
 * @param spec Dither amplitude, frequency and gain (each > 0), direction, filter corners (≥ 0),
 *             initial operating point and demodulation phase
 * @param Ts   Sample time [s] (> 0)
 * @param opts InputLimits{.lower = {lo}, .upper = {hi}}: the limits on the applied input, and on û
 * @return ESCResult; refused with invalid_sample_time for a Ts that is not finite and positive,
 *         and with invalid_parameter for a spec entry outside its domain, an InputLimits that is
 *         not a box (an end not finite, or lower above upper) or a dither faster than Ts allows.
 *         An option other than one InputLimits on one input does not compile.
 */
template<typename T, typename... Opts>
[[nodiscard]] constexpr ESCResult<T> esc(const ESCSpec<T>& spec, T Ts, const Opts&... opts) {
    detail::check_siso_input_limits<Opts...>();
    const ESCLimits<T> limits{.input = detail::siso_input_limits_of<T>(opts...)};
    if (!damp::finite_positive(Ts)) {
        return detail::refused<ESCResult<T>>(DesignStatus::invalid_sample_time);
    }
    const bool in_domain = damp::finite_positive(spec.dither_amplitude) && damp::finite_positive(spec.dither_omega)
                        && damp::finite_positive(spec.gain) && damp::finite_non_negative(spec.hp_omega)
                        && damp::finite_non_negative(spec.lp_omega) && damp::isfinite(spec.u_init)
                        && damp::isfinite(spec.demod_phase) && limits.valid();
    if (!in_domain) {
        return detail::refused<ESCResult<T>>(DesignStatus::invalid_parameter);
    }
    const T      omega = spec.dither_omega;
    const T      hp = (spec.hp_omega > T{0}) ? spec.hp_omega : omega / T{5};
    ESCResult<T> out{
        .dither_amplitude = spec.dither_amplitude,
        .dither_omega = omega,
        .gain = spec.gain,
        .hp_alpha = detail::esc_lpf_alpha(hp, Ts),
        .lp_alpha = detail::esc_lpf_alpha(spec.lp_omega, Ts),
        .demod_phase = spec.demod_phase,
        .direction = (spec.type == ExtremumType::Maximize) ? T{1} : T{-1},
        .u_init = spec.u_init,
        .limits = limits,
        .Ts = Ts,
    };
    const bool resolved = (omega * Ts <= damp::numbers::pi_v<T> / T{4}) && damp::isfinite(out.hp_alpha)
                       && damp::isfinite(out.lp_alpha);
    detail::set_status(out, resolved ? DesignStatus::ok : DesignStatus::invalid_parameter);
    return out;
}

/**
 * @brief MPPT-flavored ESC: maximize a power measurement by perturbing the
 *        operating point (e.g. converter duty or reference voltage).
 *
 * @ref esc with two MPPT defaults: it always maximizes (spec.type is ignored), and the gradient
 * low-pass sits at dither_omega/5 unless spec.lp_omega sets one. Set u_init in the spec and the
 * band in an InputLimits — the principled form of perturb-and-observe.
 *
 * @param spec Dither, gain and initial operating point, as for @ref esc
 * @param Ts   Sample time [s] (> 0)
 * @param opts InputLimits{.lower = {lo}, .upper = {hi}}, the band
 */
template<typename T, typename... Opts>
[[nodiscard]] constexpr ESCResult<T> esc_mppt(const ESCSpec<T>& spec, T Ts, const Opts&... opts) {
    ESCSpec<T> mppt = spec;
    mppt.type = ExtremumType::Maximize;
    if (damp::isfinite(mppt.lp_omega) && mppt.lp_omega == T{0}) {
        mppt.lp_omega = spec.dither_omega / T{5};
    }
    return esc(mppt, Ts, opts...);
}

} // namespace design

/**
 * @ingroup controllers
 * @brief Extremum-seeking controller runtime (model-free online optimizer).
 *
 * Usage each control tick: apply the controller's input to the plant, measure the
 * objective, and feed it back to get the next input:
 *
 * @code
 * auto esc = ExtremumSeekingController(design::esc({.dither_amplitude = 0.05, .dither_omega = 251.0, .gain = 60.0}, Ts));
 * double u = esc.input();                 // first input (= û₀)
 * for (;;) {
 *     double J = measure_objective(u);    // apply u, read the objective
 *     u = esc.step(J);                    // -> next perturbed input
 * }
 * double u_opt = esc.operating_point();   // converged operating point û
 * @endcode
 *
 * `step(J, measurement_valid)` *freezes* the estimator when the measurement is
 * flagged bad (sensor fault / transient): û and the high-/low-pass states hold
 * (the bad J never enters them) while the dither keeps running — the "freeze on
 * degraded measurement" safety behavior.
 *
 * | Member | What it is |
 * | --- | --- |
 * | step(J, valid) | one tick: the objective J measured at the last input, the next input returned |
 * | track(u) | in place of step while another source drives the plant: û starts from u |
 * | back_calculate(u_sat, u_ff) | a clamp after the law: û follows the input applied |
 * | limits(), set_limits(lim) | the applied input's and û's limits (design::ESCLimits); change them at run time |
 * | input(), operating_point(), gradient() | the input to apply; û; the filtered gradient estimate |
 * | last_command(), valid(), reset() | the last input; whether the design ran; back to û₀ |
 *
 * It is an optimizer, not a SISO law: it has no reference, and its tick takes the objective, so
 * its handover is track(u). While an operator or another law drives the plant at u, track(u)
 * takes the place of step: û ← u (inside the limits), and the gradient's filters restart (the
 * high-pass seeds on the next objective, the low-pass from 0), so no stale gradient moves û when
 * it resumes. The dither holds its phase. The step after it returns u plus the dither, and seeks
 * from there.
 *
 * @tparam T Scalar type (float or double)
 */
template<typename T = float>
class ExtremumSeekingController {
public:
    constexpr ExtremumSeekingController() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr ExtremumSeekingController(const design::ESCResult<U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : ExtremumSeekingController(result.template as<T>()) {}

    /**
     * @brief From the design result, which carries the dither, the filters, the limits and T_s.
     *
     * A failed design, or one whose limits are not boxes, is ill-formed under constant
     * evaluation; at run time it leaves the controller invalid (valid() false, û held).
     */
    constexpr ExtremumSeekingController(const design::ESCResult<T>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : config_(result),
          lim_(result.limits),
          uhat_(within(result.limits, result.u_init)),
          valid_(result.success && result.limits.valid()) {
        if (!result.success && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(result.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
        if (!valid_ && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(design::DesignStatus::invalid_parameter); // a limit that is not a box
        }
        init_phasors();
        u_cmd_ = input();
    }

    /**
     * @brief Current perturbed input to apply, û + a·sin φ.
     *
     * The *applied* command is clamped into limits().input: the estimate û is bounded, but
     * û + dither could otherwise ride a full dither amplitude past a hard actuator limit. Near a
     * bound this clips the perturbation asymmetrically (a small, benign gradient bias) in
     * exchange for never commanding past the limit; limits().operating_point can keep û a dither
     * amplitude inside.
     */
    [[nodiscard]] constexpr T input() const { return lim_.input.saturate(uhat_ + (config_.dither_amplitude * dither_s_)); }

    /**
     * @brief Feed back the measured objective; returns the next input to apply.
     * @param objective         Measured J at the previously applied input.
     * @param measurement_valid If false, ignore @p objective: hold û and the
     *                          filter states, advance only the dither.
     */
    [[nodiscard]] constexpr T step(T objective, bool measurement_valid = true) {
        if (!valid_) {
            return uhat_;
        }
        // A flagged sample is not evidence: keep it out of the filter states too,
        // or it reappears as a spurious gradient once measurements recover.
        if (!measurement_valid) {
            advance_dither();
            u_cmd_ = input();
            return u_cmd_;
        }
        // Seed the DC estimate to the first objective so the high-pass starts at 0
        // instead of dumping the full operating-point level through as a spurious
        // gradient on the first few ticks.
        if (!primed_) {
            hp_state_ = objective;
            primed_ = true;
        }
        // High-pass the objective: hp = J − LPF(J), removing the operating-point DC.
        const T hp = objective - hp_state_;
        hp_state_ += config_.hp_alpha * (objective - hp_state_);

        // Demodulate with the dither that produced this objective:
        // sin(φ − d) = sin φ cos d − cos φ sin d.
        const T demod = (dither_s_ * demod_c_) - (dither_c_ * demod_s_);
        const T raw_grad = hp * demod;
        grad_state_ += config_.lp_alpha * (raw_grad - grad_state_);
        gradient_ = grad_state_;

        // Integrate the gradient toward the extremum.
        uhat_ = within(lim_, uhat_ + (config_.direction * config_.gain * gradient_ * config_.Ts));

        // Advance the dither phase and form the next perturbed input.
        advance_dither();
        u_cmd_ = input();
        return u_cmd_;
    }

    /** @brief Not a SISO law: step(r, y) would read y as the validity flag, so it does not compile. */
    constexpr T step(T, T) = delete;

    /**
     * @brief In place of step while another source drives the plant at @p u, the input applied
     *
     * û ← @p u, inside the limits, and the gradient's filters restart: the high-pass seeds on
     * the next objective and the low-pass starts from 0, so the step after it seeks from u with no
     * stale gradient. The dither holds its phase. @p u is recorded as the last command.
     */
    constexpr void track(T u) {
        if (!valid_) {
            return;
        }
        uhat_ = within(lim_, u);
        grad_state_ = T{0};
        gradient_ = T{0};
        primed_ = false;
        u_cmd_ = u;
    }

    /**
     * @brief Record the input the plant received, @p u_sat, when a stage after this law clamped
     *        or replaced its command, and the feedforward @p u_ff added to the command before it
     *
     * The operating point moves by what the clamp cut off, û ← û + u_sat − (u + u_ff) with u the
     * input step() returned, then into the limits: û stays where the applied input was, less the
     * dither, and does not integrate on behind the clamp. An input applied as commanded changes
     * nothing.
     */
    constexpr void back_calculate(T u_sat, T u_ff = T{0}) {
        const T asked = input() + u_ff;
        const T diff = u_sat - asked;
        const T tol = T{4} * eps<T>() * (damp::abs(u_sat) + damp::abs(asked));
        if (valid_ && damp::abs(diff) > tol) {
            uhat_ = within(lim_, uhat_ + diff);
        }
    }

    /** @brief The input the last step() returned, or the one track() recorded; before the first step, input(). */
    [[nodiscard]] constexpr T last_command() const { return u_cmd_; }

    /** @brief The limits in force: the design's, or the last set_limits. */
    [[nodiscard]] constexpr design::ESCLimits<T> limits() const { return lim_; }

    /**
     * @brief Change the limits at run time; the next tick clamps into them, and û moves inside them at once
     *
     * Refused (false, nothing changed) when an end of either is NaN or ±Inf, or a lower end is
     * above its upper one.
     */
    [[nodiscard]] constexpr bool set_limits(const design::ESCLimits<T>& lim) {
        if (!lim.valid()) {
            return false;
        }
        lim_ = lim;
        uhat_ = within(lim_, uhat_);
        return true;
    }

    /** @brief The operating point û: the current estimate of the extremum's input, without the dither. */
    [[nodiscard]] constexpr T operating_point() const { return uhat_; }

    /** @brief The filtered gradient estimate the last step integrated. */
    [[nodiscard]] constexpr T gradient() const { return gradient_; }

    /** @brief False when built from a design it refuses: step() holds û. */
    [[nodiscard]] constexpr bool valid() const { return valid_; }

    /** @brief Back to û₀ (inside the limits as they stand), the dither at phase 0 and the filters cleared. */
    constexpr void reset() {
        uhat_ = within(lim_, config_.u_init);
        dither_s_ = T{0};
        dither_c_ = T{1};
        hp_state_ = T{0};
        grad_state_ = T{0};
        gradient_ = T{0};
        primed_ = false;
        u_cmd_ = input();
    }

private:
    /** @brief @p x inside the operating point's limits and the input's: the input's win where the two do not meet. */
    [[nodiscard]] static constexpr T within(const design::ESCLimits<T>& lim, T x) { return lim.input.saturate(lim.operating_point.saturate(x)); }

    // The dither phase φ advances by a fixed ωTs each tick, so it is carried as a
    // unit phasor (sin φ, cos φ) rotated by precomputed cos/sin(ωTs): the tick
    // path needs no trig and no fmod. Trig runs once, at construction.
    constexpr void init_phasors() {
        const T step = config_.dither_omega * config_.Ts;
        rot_c_ = damp::cos(step);
        rot_s_ = damp::sin(step);
        demod_c_ = damp::cos(config_.demod_phase);
        demod_s_ = damp::sin(config_.demod_phase);
    }

    // φ ← φ + ωTs, then pull |phasor| back to 1 with one Newton step of
    // 1/√(s² + c²) (no sqrt; holds the norm at 1 to roundoff indefinitely).
    constexpr void advance_dither() {
        const T s = (dither_s_ * rot_c_) + (dither_c_ * rot_s_);
        const T c = (dither_c_ * rot_c_) - (dither_s_ * rot_s_);
        const T k = (T{3} - ((s * s) + (c * c))) * static_cast<T>(0.5);
        dither_s_ = s * k;
        dither_c_ = c * k;
    }

    design::ESCResult<T> config_{}; ///< The design's realization
    design::ESCLimits<T> lim_{};    ///< The limits in force

    T    uhat_{T{0}};     // Converged operating point estimate û (the optimizer's answer).
    T    u_cmd_{T{0}};    // The input the last step() returned, or the one track() recorded.
    T    dither_s_{T{0}}; // sin φ of the current dither phase
    T    dither_c_{T{1}}; // cos φ
    T    rot_c_{T{1}};    // cos(ωTs), sin(ωTs): one-tick phase rotation
    T    rot_s_{T{0}};
    T    demod_c_{T{1}}; // cos d, sin d of the demodulation phase offset
    T    demod_s_{T{0}};
    T    hp_state_{T{0}};
    T    grad_state_{T{0}};
    T    gradient_{T{0}}; // Most recent (filtered) gradient estimate.
    bool valid_{false};
    bool primed_{false}; // Whether hp_state_ has been seeded to the first objective.
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<typename T>
ExtremumSeekingController(const design::ESCResult<T>&) -> ExtremumSeekingController<T>;

namespace detail {

/** @brief Whether @p R is an extremum-seeking design result. */
template<typename R>
inline constexpr bool is_esc_result = false;

/** @brief design::esc's result. */
template<typename U>
inline constexpr bool is_esc_result<design::ESCResult<U>> = true;

/** @brief An extremum-seeking design (design::esc, design::esc_mppt). */
template<const auto& Res>
concept esc_design = is_esc_result<design_t<Res>>;

} // namespace detail

/** @brief Law<kRes, T> of an extremum-seeking design initialized at start-up: the ExtremumSeekingController it runs as. A runtime assigned replaces it. */
template<const auto& Res, typename T>
    requires detail::esc_design<Res>
struct Law<Res, T> : ExtremumSeekingController<T> {
    using Runtime = ExtremumSeekingController<T>; ///< The runtime the design runs as

    /** @brief The design's law. */
    constexpr Law() : Runtime{Res} {}

    using Runtime::operator=;
};

/**
 * @brief Law<kRes, T> of a `constexpr` extremum-seeking design: the same runtime, a failed design
 *        refused wherever the Law is declared
 */
// ponytail: a constexpr design's coefficients (the rotation, the filters, k·T_s) could be constants of the type, as a PID Law's gains are; they stay members here.
template<const auto& Res, typename T>
    requires(detail::esc_design<Res> && detail::constant_design<Res>)
struct Law<Res, T> : ExtremumSeekingController<T> {
    using Runtime = ExtremumSeekingController<T>; ///< The runtime the design runs as

    static_assert(detail::constant_design_runs(Res.success, Res.status, Res.Ts));

    /** @brief The design's law. */
    constexpr Law() : Runtime{Res} {}

    using Runtime::operator=;
};

} // namespace damp
