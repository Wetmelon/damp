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

#include <type_traits>

#include "damp/math/math.hpp"            // damp::sin, damp::fmod, damp::clamp, numbers
#include "damp/matrix/matrix_traits.hpp" // scalar_type_t, default_tol

namespace damp {

/// Whether ESC climbs to a maximum or descends to a minimum of the objective.
enum class ExtremumType {
    Maximize,
    Minimize
};

namespace design {

/**
 * @brief Extremum-seeking controller configuration (discrete realization).
 * @tparam T Scalar type
 */
template<typename T = double>
struct ESCConfig {

    using scalar = scalar_type_t<T>;

    scalar dither_amplitude{static_cast<scalar>(0.01)}; ///< a, perturbation amplitude (input units)
    scalar dither_omega{scalar{0}};                     ///< ω, perturbation frequency [rad/s]
    scalar gain{scalar{0}};                             ///< k, integrator (adaptation) gain
    scalar hp_alpha{scalar{0}};                         ///< high-pass coefficient (objective DC removal)
    scalar lp_alpha{scalar{1}};                         ///< gradient low-pass coefficient (1 = no LPF)
    scalar demod_phase{scalar{0}};                      ///< φ, demodulation phase offset [rad] (plant lag)
    scalar direction{scalar{1}};                        ///< +1 maximize, −1 minimize
    scalar Ts{scalar{0}};                               ///< sample time [s]
    scalar u_init{scalar{0}};                           ///< initial operating point û₀
    scalar u_min{scalar{0}};                            ///< lower clamp on û (active iff u_min < u_max)
    scalar u_max{scalar{0}};                            ///< upper clamp on û

    [[nodiscard]] constexpr bool valid() const {
        if (dither_amplitude <= scalar{0} || dither_omega <= scalar{0} || Ts <= scalar{0}) {
            return false;
        }
        if (gain <= scalar{0}) {
            return false;
        }
        if (hp_alpha <= scalar{0} || hp_alpha > scalar{1}) {
            return false;
        }
        if (lp_alpha <= scalar{0} || lp_alpha > scalar{1}) {
            return false;
        }
        // The dither must be well-resolved by the sample rate — not merely below
        // Nyquist. ESC needs the perturbation cleanly sampled to demodulate it, so
        // require ≥ 8 samples per dither period (ω·Ts ≤ π/4). Below that the gradient
        // estimate degrades and the scale separation (HPF < ω < Nyquist) collapses.
        if (dither_omega * Ts > damp::numbers::pi_v<scalar> / scalar{4}) {
            return false;
        }
        return true;
    }
};

/**
 * @brief Design result for the extremum-seeking controller.
 */
template<typename T = double>
struct ESCResult {

    ESCConfig<T> config{};
    bool         success{false};

    template<typename U>
    [[nodiscard]] constexpr ESCResult<std::remove_const_t<U>> as() const {
        using out_t = std::remove_const_t<U>;
        using os = scalar_type_t<out_t>;
        ESCResult<out_t> out{};
        out.config.dither_amplitude = static_cast<os>(config.dither_amplitude);
        out.config.dither_omega = static_cast<os>(config.dither_omega);
        out.config.gain = static_cast<os>(config.gain);
        out.config.hp_alpha = static_cast<os>(config.hp_alpha);
        out.config.lp_alpha = static_cast<os>(config.lp_alpha);
        out.config.demod_phase = static_cast<os>(config.demod_phase);
        out.config.direction = static_cast<os>(config.direction);
        out.config.Ts = static_cast<os>(config.Ts);
        out.config.u_init = static_cast<os>(config.u_init);
        out.config.u_min = static_cast<os>(config.u_min);
        out.config.u_max = static_cast<os>(config.u_max);
        out.success = success;
        return out;
    }
};

/// First-order discrete low-pass coefficient for corner @p wc [rad/s] at @p Ts.
template<typename T>
[[nodiscard]] constexpr T esc_lpf_alpha(T wc, T Ts) {
    if (wc <= T{0}) {
        return T{1}; // no filtering (pass-through)
    }
    return (wc * Ts) / (T{1} + (wc * Ts));
}

/**
 * @brief Synthesize an extremum-seeking controller.
 *
 * @param dither_amplitude  a, perturbation amplitude in input units (> 0)
 * @param dither_freq_hz    perturbation frequency [Hz] (> 0, below Nyquist)
 * @param gain              k, adaptation/integrator gain (> 0)
 * @param Ts                sample time [s] (> 0)
 * @param type              Maximize (climb) or Minimize (descend)
 * @param hp_cutoff_hz      high-pass corner [Hz] (0 → dither_freq/5); below the dither
 * @param lp_cutoff_hz      gradient low-pass corner [Hz] (0 → no LPF)
 * @param u_init            initial operating point û₀
 * @param u_min,u_max       clamp on û (active iff u_min < u_max)
 * @param demod_phase       demodulation phase offset [rad] for plant lag (default 0)
 */
template<typename T = double>
[[nodiscard]] constexpr ESCResult<T> esc(
    T            dither_amplitude,
    T            dither_freq_hz,
    T            gain,
    T            Ts,
    ExtremumType type = ExtremumType::Maximize,
    T            hp_cutoff_hz = T{0},
    T            lp_cutoff_hz = T{0},
    T            u_init = T{0},
    T            u_min = T{0},
    T            u_max = T{0},
    T            demod_phase = T{0}
) {
    if (dither_freq_hz <= T{0} || Ts <= T{0}) {
        return ESCResult<T>{};
    }
    const T two_pi = T{2} * damp::numbers::pi_v<T>;
    const T hp_hz = (hp_cutoff_hz > T{0}) ? hp_cutoff_hz : dither_freq_hz / T{5};

    using scalar = scalar_type_t<T>;
    const ESCConfig<T> config{
        .dither_amplitude = static_cast<scalar>(dither_amplitude),
        .dither_omega = static_cast<scalar>(two_pi * dither_freq_hz),
        .gain = static_cast<scalar>(gain),
        .hp_alpha = static_cast<scalar>(esc_lpf_alpha(two_pi * hp_hz, Ts)),
        .lp_alpha = static_cast<scalar>(esc_lpf_alpha(two_pi * lp_cutoff_hz, Ts)),
        .demod_phase = static_cast<scalar>(demod_phase),
        .direction = static_cast<scalar>((type == ExtremumType::Maximize) ? T{1} : T{-1}),
        .Ts = static_cast<scalar>(Ts),
        .u_init = static_cast<scalar>(u_init),
        .u_min = static_cast<scalar>(u_min),
        .u_max = static_cast<scalar>(u_max),
    };
    return ESCResult<T>{
        .config = config,
        .success = config.valid(),
    };
}

/**
 * @brief MPPT-flavored ESC: maximize a power measurement by perturbing the
 *        operating point (e.g. converter duty or reference voltage).
 *
 * Thin wrapper over @ref esc with Maximize and a clamp band — the
 * principled form of perturb-and-observe.
 */
template<typename T = double>
[[nodiscard]] constexpr ESCResult<T> esc_mppt(
    T dither_amplitude,
    T dither_freq_hz,
    T gain,
    T Ts,
    T u_init,
    T u_min,
    T u_max
) {
    return esc<T>(dither_amplitude, dither_freq_hz, gain, Ts, ExtremumType::Maximize, T{0}, dither_freq_hz / T{5}, u_init, u_min, u_max);
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
 * auto esc = ExtremumSeekingController(design::esc(...));
 * double u = esc.input();                 // first input (= û₀)
 * for (;;) {
 *     double J = measure_objective(u);    // apply u, read the objective
 *     u = esc.step(J);                    // -> next perturbed input
 * }
 * double u_opt = esc.estimate();          // converged operating point û
 * @endcode
 *
 * `step(J, measurement_valid)` *freezes* the integrator when the measurement is
 * flagged bad (sensor fault / transient), holding û while still dithering — the
 * "freeze on degraded measurement" safety behavior.
 *
 * @tparam T Scalar type (float or double)
 */
template<typename T = float>
class ExtremumSeekingController {
public:
    constexpr ExtremumSeekingController() = default;

    constexpr explicit ExtremumSeekingController(const design::ESCConfig<T>& config)
        : config_(config), uhat_(config.u_init), valid_(config.valid()) {}

    constexpr explicit ExtremumSeekingController(const design::ESCResult<T>& design)
        : config_(design.config), uhat_(design.config.u_init), valid_(design.success && design.config.valid()) {}

    /// Current perturbed input to apply, û + a·sin(phase).
    ///
    /// When a clamp band is configured (`u_min < u_max`) the *applied* command is
    /// clamped to it: the estimate û is bounded, but û + dither could otherwise
    /// ride a full dither amplitude past a hard actuator limit. Near a bound this
    /// clips the perturbation asymmetrically (a small, benign gradient bias) in
    /// exchange for never commanding past the limit.
    [[nodiscard]] constexpr T input() const {
        const T u = uhat_ + (config_.dither_amplitude * damp::sin(phase_));
        if (config_.u_min < config_.u_max) {
            return damp::clamp(u, config_.u_min, config_.u_max);
        }
        return u;
    }

    /**
     * @brief Feed back the measured objective; returns the next input to apply.
     * @param objective         Measured J at the previously applied input.
     * @param measurement_valid If false, freeze the integrator (hold û) this tick.
     */
    [[nodiscard]] constexpr T step(T objective, bool measurement_valid = true) {
        if (!valid_) {
            return uhat_;
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

        // Demodulate with the dither that produced this objective (phase_).
        const T demod = damp::sin(phase_ - config_.demod_phase);
        const T raw_grad = hp * demod;
        grad_state_ += config_.lp_alpha * (raw_grad - grad_state_);
        gradient_ = grad_state_;

        // Integrate the gradient toward the extremum (frozen on bad measurement).
        if (measurement_valid) {
            uhat_ += config_.direction * config_.gain * gradient_ * config_.Ts;
            if (config_.u_min < config_.u_max) {
                uhat_ = damp::clamp(uhat_, config_.u_min, config_.u_max);
            }
        }

        // Advance the dither phase and form the next perturbed input.
        phase_ = damp::fmod(phase_ + (config_.dither_omega * config_.Ts), two_pi());
        return input();
    }

    [[nodiscard]] constexpr T    estimate() const { return uhat_; }
    [[nodiscard]] constexpr T    gradient() const { return gradient_; }
    [[nodiscard]] constexpr bool valid() const { return valid_; }

    constexpr void reset() {
        uhat_ = config_.u_init;
        phase_ = T{0};
        hp_state_ = T{0};
        grad_state_ = T{0};
        gradient_ = T{0};
        primed_ = false;
    }

private:
    static constexpr T two_pi() { return T{2} * damp::numbers::pi_v<T>; }

    design::ESCConfig<T> config_{};

    T    uhat_{T{0}}; // Converged operating point estimate û (the optimizer's answer).
    T    phase_{T{0}};
    T    hp_state_{T{0}};
    T    grad_state_{T{0}};
    T    gradient_{T{0}}; // Most recent (filtered) gradient estimate.
    bool valid_{false};
    bool primed_{false}; // Whether hp_state_ has been seeded to the first objective.
};

} // namespace damp
