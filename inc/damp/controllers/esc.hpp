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
 * @brief The extremum-seeking design input: dither, gain, filter corners and clamp; Ts is passed beside it
 *
 * The design input of @ref esc and @ref esc_mppt, in the caller's units (Hz, rad). The caller
 * names every number:
 * `design::esc({.dither_amplitude = 0.05, .dither_freq_hz = 40.0, .gain = 60.0, .lp_cutoff_hz = 8.0}, Ts)`.
 * The design turns it into an @ref ESCResult (rad/s, filter coefficients), the runtime
 * realization; the two are distinct types.
 *
 * dither_amplitude, dither_freq_hz and gain have no usable default: left at 0 the design is
 * refused.
 *
 * @tparam T Scalar type
 */
template<typename T = double>
struct ESCSpec {
    T            dither_amplitude{};           ///< a, perturbation amplitude in input units (> 0)
    T            dither_freq_hz{};             ///< Perturbation frequency [Hz] (> 0, at least 8 samples per period)
    T            gain{};                       ///< k, adaptation (integrator) gain (> 0)
    ExtremumType type{ExtremumType::Maximize}; ///< Maximize (climb) or Minimize (descend)
    T            hp_cutoff_hz{T{0}};           ///< High-pass corner [Hz], below the dither; 0 = dither_freq_hz/5
    T            lp_cutoff_hz{T{0}};           ///< Gradient low-pass corner [Hz]; 0 = no low-pass
    T            u_init{T{0}};                 ///< Initial operating point û₀
    Bounds<1, T> u_lim{};                      ///< Clamp on û (default unbounded; off when lower ≥ upper)
    T            demod_phase{T{0}};            ///< Demodulation phase offset φ [rad] for plant lag
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
    Bounds<1, T> u_lim{};                             ///< Clamp on û (default unbounded; off when lower ≥ upper)
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
            .u_lim = u_lim.template as<out_t>(),
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
 * @param spec Dither amplitude, frequency and gain (each > 0), direction, filter corners (≥ 0),
 *             initial operating point, clamp and demodulation phase
 * @param Ts   Sample time [s] (> 0)
 * @return ESCResult; refused with invalid_sample_time for a non-positive or non-finite Ts, and
 *         with invalid_parameter for a spec entry outside its domain or a dither faster than Ts allows
 */
template<typename T = double>
[[nodiscard]] constexpr ESCResult<T> esc(const ESCSpec<T>& spec, T Ts) {
    if (!damp::finite_positive(Ts)) {
        return ::damp::detail::refused<ESCResult<T>>(DesignStatus::invalid_sample_time);
    }
    const bool in_domain = damp::finite_positive(spec.dither_amplitude) && damp::finite_positive(spec.dither_freq_hz)
                        && damp::finite_positive(spec.gain) && damp::finite_non_negative(spec.hp_cutoff_hz)
                        && damp::finite_non_negative(spec.lp_cutoff_hz) && damp::isfinite(spec.u_init)
                        && damp::isfinite(spec.demod_phase);
    if (!in_domain) {
        return ::damp::detail::refused<ESCResult<T>>(DesignStatus::invalid_parameter);
    }
    const T      two_pi = T{2} * damp::numbers::pi_v<T>;
    const T      omega = two_pi * spec.dither_freq_hz;
    const T      hp_hz = (spec.hp_cutoff_hz > T{0}) ? spec.hp_cutoff_hz : spec.dither_freq_hz / T{5};
    ESCResult<T> out{
        .dither_amplitude = spec.dither_amplitude,
        .dither_omega = omega,
        .gain = spec.gain,
        .hp_alpha = detail::esc_lpf_alpha(two_pi * hp_hz, Ts),
        .lp_alpha = detail::esc_lpf_alpha(two_pi * spec.lp_cutoff_hz, Ts),
        .demod_phase = spec.demod_phase,
        .direction = (spec.type == ExtremumType::Maximize) ? T{1} : T{-1},
        .u_init = spec.u_init,
        .u_lim = spec.u_lim,
        .Ts = Ts,
    };
    const bool resolved = (omega * Ts <= damp::numbers::pi_v<T> / T{4}) && damp::isfinite(out.hp_alpha)
                       && damp::isfinite(out.lp_alpha);
    ::damp::detail::set_status(out, resolved ? DesignStatus::ok : DesignStatus::invalid_parameter);
    return out;
}

/**
 * @brief MPPT-flavored ESC: maximize a power measurement by perturbing the
 *        operating point (e.g. converter duty or reference voltage).
 *
 * @ref esc with two MPPT defaults: it always maximizes (spec.type is ignored), and the gradient
 * low-pass sits at dither_freq_hz/5 unless spec.lp_cutoff_hz sets one. Set u_init and the clamp
 * band u in the spec — the principled form of perturb-and-observe.
 *
 * @param spec Dither, gain, initial operating point and clamp band, as for @ref esc
 * @param Ts   Sample time [s] (> 0)
 */
template<typename T = double>
[[nodiscard]] constexpr ESCResult<T> esc_mppt(const ESCSpec<T>& spec, T Ts) {
    ESCSpec<T> mppt = spec;
    mppt.type = ExtremumType::Maximize;
    if (mppt.lp_cutoff_hz == T{0}) {
        mppt.lp_cutoff_hz = spec.dither_freq_hz / T{5};
    }
    return esc(mppt, Ts);
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
 * auto esc = ExtremumSeekingController(design::esc({.dither_amplitude = 0.05, .dither_freq_hz = 40.0, .gain = 60.0}, Ts));
 * double u = esc.input();                 // first input (= û₀)
 * for (;;) {
 *     double J = measure_objective(u);    // apply u, read the objective
 *     u = esc.control(J);                 // -> next perturbed input
 * }
 * double u_opt = esc.operating_point();   // converged operating point û
 * @endcode
 *
 * `control(J, measurement_valid)` *freezes* the estimator when the measurement is
 * flagged bad (sensor fault / transient): û and the high-/low-pass states hold
 * (the bad J never enters them) while the dither keeps running — the "freeze on
 * degraded measurement" safety behavior.
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

    constexpr ExtremumSeekingController(const design::ESCResult<T>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : config_(result), uhat_(result.u_init), valid_(result.success) {
        if (!result.success && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(result.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
        init_phasors();
    }

    /**
     * @brief Current perturbed input to apply, û + a·sin φ.
     *
     * When the clamp band is active (`u_lim.lower < u_lim.upper`) the *applied* command is
     * clamped to it: the estimate û is bounded, but û + dither could otherwise
     * ride a full dither amplitude past a hard actuator limit. Near a bound this
     * clips the perturbation asymmetrically (a small, benign gradient bias) in
     * exchange for never commanding past the limit.
     */
    [[nodiscard]] constexpr T input() const {
        const T u = uhat_ + (config_.dither_amplitude * dither_s_);
        if (config_.u_lim.lower(0) < config_.u_lim.upper(0)) {
            return config_.u_lim.saturate(u);
        }
        return u;
    }

    /**
     * @brief Feed back the measured objective; returns the next input to apply.
     * @param objective         Measured J at the previously applied input.
     * @param measurement_valid If false, ignore @p objective: hold û and the
     *                          filter states, advance only the dither.
     */
    [[nodiscard]] constexpr T control(T objective, bool measurement_valid = true) {
        if (!valid_) {
            return uhat_;
        }
        // A flagged sample is not evidence: keep it out of the filter states too,
        // or it reappears as a spurious gradient once measurements recover.
        if (!measurement_valid) {
            advance_dither();
            return input();
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
        uhat_ += config_.direction * config_.gain * gradient_ * config_.Ts;
        if (config_.u_lim.lower(0) < config_.u_lim.upper(0)) {
            uhat_ = config_.u_lim.saturate(uhat_);
        }

        // Advance the dither phase and form the next perturbed input.
        advance_dither();
        return input();
    }

    /** @brief Not a SISO law: control(r, y) would read y as the validity flag, so it does not compile. */
    constexpr T control(T, T) = delete;

    /** @brief The operating point û: the current estimate of the extremum's input, without the dither. */
    [[nodiscard]] constexpr T    operating_point() const { return uhat_; }
    [[nodiscard]] constexpr T    gradient() const { return gradient_; }
    [[nodiscard]] constexpr bool valid() const { return valid_; }

    constexpr void reset() {
        uhat_ = config_.u_init;
        dither_s_ = T{0};
        dither_c_ = T{1};
        hp_state_ = T{0};
        grad_state_ = T{0};
        gradient_ = T{0};
        primed_ = false;
    }

private:
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

    T    uhat_{T{0}};     // Converged operating point estimate û (the optimizer's answer).
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

} // namespace damp
