// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file pr.hpp
 * @brief Proportional-resonant controller
 *
 * @defgroup pr_controller Proportional-Resonant Controller
 * @brief Proportional-resonant controller for AC tracking (multi-harmonic: HarmonicSuppressor)
 *
 * PR controllers provide infinite gain at a specific resonant frequency,
 * achieving zero steady-state error for sinusoidal references. They are
 * the AC equivalent of PI controllers for DC references.
 *
 * Common applications:
 * - Grid-tied inverter current control
 * - Active power filter harmonic compensation
 * - Rotating reference frame alternatives (no Park transform needed)
 *
 * Transfer function (non-ideal PR, @f$\omega_c > 0@f$):
 * @f[
 *   C(s) = K_p + \frac{2 K_i \omega_c\, s}{s^2 + 2 \omega_c s + \omega_0^2}
 * @f]
 * where @f$K_p@f$ is the proportional gain, @f$K_i@f$ the resonant (integral)
 * gain, @f$\omega_0@f$ the resonant frequency (e.g. @f$2\pi\cdot 50@f$ for a
 * 50 Hz grid) and @f$\omega_c@f$ the resonant-term bandwidth.
 *
 * Ideal PR (@f$\omega_c = 0@f$):
 * @f[
 *   C(s) = K_p + \frac{K_i s}{s^2 + \omega_0^2}
 * @f]
 * Infinite gain exactly at @f$\omega_0@f$; fragile to grid-frequency drift, so
 * non-ideal PR with a few rad/s of @f$\omega_c@f$ is preferred in practice.
 *
 * MATLAB® equivalent (non-ideal):
 * @code{.m}
 *   s = tf('s');
 *   C = Kp + 2*Ki*wc*s / (s^2 + 2*wc*s + w0^2);
 * @endcode
 *
 * @see R. Teodorescu et al., "Proportional-resonant controllers and filters for
 *      grid-connected voltage-source converters," IEE Proc. Electr. Power Appl.,
 *      2006. DOI: 10.1049/ip-epa:20060008
 * @see D. N. Zmood and D. G. Holmes, "Stationary frame current regulation of PWM
 *      inverters with zero steady-state error," IEEE Trans. Power Electron.,
 *      2003. DOI: 10.1109/TPEL.2003.810852
 */

#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/design/status.hpp"
#include "damp/filters/filters.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/core.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

namespace damp {

namespace design {

namespace detail {

/**
 * @brief Pre-warped resonant frequency (2/Ts)·tan(ω₀Ts/2) that PRController discretizes at
 *
 * The bilinear map s → (2/Ts)(z−1)/(z+1) of a resonator rebuilt at this frequency places the
 * digital peak exactly on ω₀. Returns ω₀ unchanged when Ts ≤ 0 (continuous) or when ω₀ is at or
 * above Nyquist (tan has no finite value there).
 */
template<typename T>
[[nodiscard]] constexpr T pr_prewarped_w0(T w0, T Ts) {
    const T half = w0 * Ts / T{2};
    if (!(Ts > T{0}) || !(half < damp::numbers::pi_v<T> / T{2})) {
        return w0;
    }
    return (T{2} / Ts) * damp::tan(half);
}

} // namespace detail

/**
 * @struct PRSpec
 * @brief The gains and frequencies of a proportional-resonant design; Ts is passed beside it
 *
 * The caller names every number: `design::pr({.Kp = 1.0, .Ki = 200.0, .w0 = w, .wc = 10.0}, Ts)`.
 * For @ref pr_harmonics and @ref harmonic_suppressor it describes the fundamental: Ki and w0 are
 * the fundamental's, and each harmonic h gets Ki/h at h·w0.
 */
template<typename T = double>
struct PRSpec {
    T Kp{}; ///< Proportional gain (finite)
    T Ki{}; ///< Resonant (integral) gain (finite)
    T w0{}; ///< Resonant frequency [rad/s], e.g. 2π·50 for a 50 Hz grid (> 0)
    T wc{}; ///< Resonant-term bandwidth [rad/s] (≥ 0); 0 = ideal PR, 5–15 rad/s typical non-ideal
    T Tt{}; ///< Back-calculation tracking time [s] (≥ 0): ΔR = (Ts/Tt)(u_sat − (u + u_ff)); 0 = none
};

/**
 * @struct PRResult
 * @brief Proportional-Resonant controller design result
 */
template<typename T = double>
struct [[nodiscard]] PRResult {
    T            Kp{};                                ///< Proportional gain
    T            Ki{};                                ///< Resonant (integral) gain
    T            w0{};                                ///< Resonant frequency (rad/s)
    T            wc{};                                ///< Cutoff bandwidth of resonant term (rad/s); 0 = ideal PR
    T            Ts{};                                ///< Sampling time (seconds); 0 = continuous-time design
    T            Tt{};                                ///< Back-calculation tracking time [s]; 0 = none
    bool         success{false};                      ///< status == DesignStatus::ok
    DesignStatus status{DesignStatus::not_converged}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return PRResult<U>{
            static_cast<U>(Kp), static_cast<U>(Ki), static_cast<U>(w0), static_cast<U>(wc),
            static_cast<U>(Ts), static_cast<U>(Tt), success, status
        };
    }

    /**
     * @brief Convert to continuous-time transfer function
     *
     * Non-ideal (@f$\omega_c > 0@f$):
     *   C(s) = Kp + 2*Ki*wc*s / (s² + 2*wc*s + w0²)
     *   num = {Kp, 2*wc*(Kp+Ki), Kp*w0²}, den = {1, 2*wc, w0²}
     *
     * Ideal (@f$\omega_c = 0@f$):
     *   C(s) = Kp + Ki*s / (s² + w0²)
     *   num = {Kp, Ki, Kp*w0²}, den = {1, 0, w0²}
     *
     * Numerator and denominator highest power of s first.
     */
    [[nodiscard]] constexpr TransferFunction<3, 3, T> to_tf() const {
        if (wc > T{0}) {
            return TransferFunction<3, 3, T>{
                .num = {Kp, T{2} * wc * (Kp + Ki), Kp * w0 * w0},
                .den = {T{1}, T{2} * wc, w0 * w0},
            };
        }
        // Ideal PR: C(s) = Kp + Ki*s/(s² + w0²)
        return TransferFunction<3, 3, T>{
            .num = {Kp, Ki, Kp * w0 * w0},
            .den = {T{1}, T{0}, w0 * w0},
        };
    }

    /**
     * @brief Convert to continuous-time state-space (2nd order SISO)
     *
     * Controllable canonical form of C(s); D = Kp carries the proportional
     * feedthrough. Resonant C row is [0, 2*Ki*wc] (non-ideal) or [0, Ki] (ideal).
     */
    [[nodiscard]] constexpr StateSpace<2, 1, 1, T> to_ss() const {
        const T c_res = (wc > T{0}) ? (T{2} * Ki * wc) : Ki;
        return StateSpace<2, 1, 1, T>{
            .A = Matrix<2, 2, T>{{T{0}, T{1}}, {-w0 * w0, -T{2} * wc}},
            .B = Matrix<2, 1, T>{{T{0}}, {T{1}}},
            .C = Matrix<1, 2, T>{{T{0}, c_res}},
            .D = Matrix<1, 1, T>{{Kp}},
        };
    }

    /**
     * @brief Discrete-time state-space of the PRController runtime (Tustin, ω₀ pre-warped)
     *
     * Discretizes the way @ref PRController does: the resonator is rebuilt at the pre-warped
     * frequency @f$\omega_0' = (2/T_s)\tan(\omega_0 T_s/2)@f$ and mapped with
     * @f$s = (2/T_s)(z-1)/(z+1)@f$, so the digital resonance sits exactly on @f$\omega_0@f$
     * and the model's response equals the runtime's tick for tick. A resonance at or above
     * Nyquist is not pre-warped (as in the runtime). Ts = 0 (continuous design) returns
     * @ref to_ss. For another discretization use `discretize(to_ss(), Ts, method)`.
     */
    [[nodiscard]] constexpr StateSpace<2, 1, 1, T> to_discrete_ss() const {
        PRResult warped = *this;
        warped.w0 = detail::pr_prewarped_w0(w0, Ts);
        return (Ts > T{0}) ? *discretize(warped.to_ss(), Ts, DiscretizationMethod::Tustin) : to_ss();
    }
};

/**
 * @brief Design a Proportional-Resonant controller
 *
 * The scalar type follows @p Ts, so a braced spec needs no template argument:
 * `design::pr({.Kp = 1.0, .Ki = 200.0, .w0 = w, .wc = 10.0}, 1e-4)`.
 *
 * @param spec Kp and Ki finite, w0 > 0 (below Nyquist π/Ts when Ts > 0), wc ≥ 0 (0 = ideal PR),
 *             Tt ≥ 0 (0 = no back-calculation)
 * @param Ts   Sampling time [s]; 0 = continuous-time design (analysis only: the runtime refuses it). Must be ≥ 0.
 * @return PRResult with design parameters; refused with invalid_sample_time for a negative or
 *         non-finite Ts, and with invalid_parameter for any other input outside its domain
 */
template<typename T = double>
[[nodiscard]] constexpr PRResult<T> pr(const PRSpec<T>& spec, T Ts) {
    if (!damp::finite_non_negative(Ts)) {
        return ::damp::detail::refused<PRResult<T>>(DesignStatus::invalid_sample_time);
    }
    const bool ok = damp::isfinite(spec.Kp) && damp::isfinite(spec.Ki) && damp::finite_positive(spec.w0)
                 && damp::finite_non_negative(spec.wc) && damp::finite_non_negative(spec.Tt);
    // At or above Nyquist the bilinear map folds the resonance onto another frequency.
    if (!ok || ((Ts > T{0}) && !(spec.w0 * Ts < damp::numbers::pi_v<T>))) {
        return ::damp::detail::refused<PRResult<T>>(DesignStatus::invalid_parameter);
    }
    return PRResult<T>{
        .Kp = spec.Kp,
        .Ki = spec.Ki,
        .w0 = spec.w0,
        .wc = spec.wc,
        .Ts = Ts,
        .Tt = spec.Tt,
        .success = true,
        .status = DesignStatus::ok,
    };
}

/**
 * @brief Design multiple-harmonic PR controller gains
 *
 * For harmonic compensation, returns an array of PRResult for harmonics
 * 1, 3, 5, 7, ... (or user-specified harmonic numbers), each refused on its own (e.g. a harmonic
 * at or above Nyquist). @ref design::harmonic_suppressor builds the deployable bank from the
 * same spec.
 *
 * @param fundamental The fundamental's resonator: Kp (carried on the first entry only), Ki and w0
 *                    of the fundamental, wc and Tt shared by every harmonic
 * @param Ts          Sampling time [s]
 * @param harmonics   Array of harmonic numbers (e.g. {1, 3, 5, 7}); harmonic h gets Ki/h at h·w0
 * @return Array of PRResult, one per harmonic
 */
template<size_t N, typename T = double>
[[nodiscard]] constexpr damp::array<PRResult<T>, N>
pr_harmonics(const PRSpec<T>& fundamental, T Ts, const damp::array<size_t, N>& harmonics) {
    damp::array<PRResult<T>, N> results{};
    for (size_t i = 0; i < N; ++i) {
        const T order = static_cast<T>(harmonics[i]);
        // Typically reduce Ki for higher harmonics
        results[i] = pr(
            PRSpec<T>{
                .Kp = (i == 0) ? fundamental.Kp : T{0},
                .Ki = fundamental.Ki / order,
                .w0 = fundamental.w0 * order,
                .wc = fundamental.wc,
                .Tt = fundamental.Tt,
            },
            Ts
        );
    }
    return results;
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Discrete Proportional-Resonant Controller
 *
 * PR controller realized as a proportional path plus a Tustin-discretized
 * resonant @ref Biquad. The resonant frequency is pre-warped so the bilinear
 * map lands the digital peak exactly on w0, preserving zero steady-state error
 * at the target frequency (the unwarped map places the peak at the
 * bilinear-warped frequency, worst for high harmonics / low fs).
 *
 * Non-ideal (@f$\omega_c > 0@f$): @f$R(s) = 2 K_i \omega_c s/(s^2 + 2\omega_c s + \omega_0^2)@f$.
 * Ideal (@f$\omega_c = 0@f$): @f$R(s) = K_i s/(s^2 + \omega_0^2)@f$.
 * Both use @f$s = (2/T_s)(z-1)/(z+1)@f$.
 *
 * A failed or continuous (Ts = 0) design does not compile under constant evaluation; at run time
 * it leaves the controller invalid (valid() false, output 0).
 *
 * @code
 * // 50 Hz grid current regulator at 10 kHz
 * auto                d = design::pr({.Kp = 1.0, .Ki = 200.0, .w0 = 2 * std::numbers::pi * 50, .wc = 10.0}, 1e-4);
 * PRController<double> pr(d);
 * double              u = pr.control(i_ref, i_meas); // (r, y) form
 * @endcode
 *
 * @tparam T Scalar type (default: float)
 */
template<typename T = float>
struct PRController {
    T    Kp{};
    T    Ki{};
    T    w0{};
    T    wc{};
    T    Ts{};
    T    Tt{T{0}};      ///< Back-calculation tracking time [s]: ΔR = (Ts/Tt)(u_sat − (u + u_ff)); 0 = none
    T    u_cmd{T{0}};   ///< Last command control() returned: the u in back_calculate's u_sat − (u + u_ff)
    bool valid_{false}; ///< From the design's success flag; gates control()

    //! Discretized resonant term (Direct Form I biquad, inert until designed).
    Biquad<T> resonant{design::SecondOrderCoeffs<T>{}};

    constexpr PRController() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!std::is_same_v<U, T>)
    constexpr PRController(const design::PRResult<U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : PRController(result.template as<T>()) {}

    constexpr PRController(const design::PRResult<T>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : Kp(result.Kp),
          Ki(result.Ki),
          w0(result.w0),
          wc(result.wc),
          Ts(result.Ts),
          Tt(result.Tt),
          valid_(result.success && result.Ts > T{0}) {
        if ((!result.success || !(result.Ts > T{0})) && std::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(result.success ? design::DesignStatus::invalid_sample_time : result.status); // constant-evaluated (constinit / constexpr): a failed or continuous design does not compile
        }
        if (valid_) {
            compute_coefficients();
        }
    }

    template<typename U>
    constexpr explicit PRController(const PRController<U>& other)
        : Kp(static_cast<T>(other.Kp)),
          Ki(static_cast<T>(other.Ki)),
          w0(static_cast<T>(other.w0)),
          wc(static_cast<T>(other.wc)),
          Ts(static_cast<T>(other.Ts)),
          Tt(static_cast<T>(other.Tt)),
          u_cmd(static_cast<T>(other.u_cmd)),
          valid_(other.valid_),
          resonant(other.resonant) {}

    /**
     * @brief Compute control output
     *
     * @param error Current error (reference - measurement)
     * @return Control output u
     */
    [[nodiscard]] constexpr T control(T error) {
        if (!valid_) {
            return T{0};
        }
        u_cmd = Kp * error + resonant(error);
        return u_cmd;
    }

    /**
     * @brief Reference-tracking overload satisfying SISOController.
     *
     * Computes the error internally and forwards to the error-form `control()`.
     * Lets this controller be used directly in `Cascade<Outer, PRController>`
     * and in tuning harnesses that speak the (r, y) protocol.
     */
    [[nodiscard]] constexpr T control(T r, T y) {
        return control(r - y);
    }

    /** @brief Clear the resonant filter's delay line and the last command; gains and config are preserved. */
    constexpr void reset() {
        resonant.reset();
        u_cmd = T{0};
    }

    [[nodiscard]] constexpr bool valid() const { return valid_; }

    /**
     * @brief Anti-windup after a downstream clamp: ΔR = (Ts/Tt)·(u_sat − (u + u_ff))
     *
     * Damped unwind: adds the excess of the applied command over this law's last command
     * u (plus any feedforward @p u_ff added downstream), times Ts / Tt, to the most recent
     * resonant output, so the next tick's resonant contribution is pulled back toward the
     * realizable command. No Ki factor: the resonant output is already in command units, so
     * Tt is a tracking time [s] (the PID laws' Kbc is a gain ratio instead). No-op when `Tt == 0`.
     *
     * @param u_sat Command actually applied after downstream clamping.
     * @param u_ff  Feedforward added to this law's command before the clamp (default 0).
     */
    constexpr void back_calculate(T u_sat, T u_ff = T{0}) {
        const T u_unsat = u_cmd + u_ff;
        if (!valid_ || Tt == T{0} || u_unsat == u_sat) {
            return;
        }
        resonant.set_last_output(resonant.last_output() + Ts * ((u_sat - u_unsat) / Tt));
    }

    /**
     * @brief Update resonant frequency (for grid frequency adaptation)
     *
     * Ignores a non-finite @p new_w0 or one outside (0, π/Ts), which the bilinear map cannot place.
     *
     * @param new_w0 New resonant frequency (rad/s)
     */
    constexpr void set_frequency(T new_w0) {
        if (!damp::finite_positive(new_w0) || !(new_w0 * Ts < damp::numbers::pi_v<T>)) {
            return;
        }
        w0 = new_w0;
        if (valid_) {
            compute_coefficients();
        }
    }

private:
    constexpr void compute_coefficients() {
        const T k = T{2} / Ts; // Tustin substitution factor 2/Ts

        // Pre-warp the resonant frequency so the bilinear map places the digital
        // peak exactly on w0 — w_pw = (2/Ts)*tan(w0*Ts/2). Without this the peak
        // lands at the warped frequency and zero-steady-state error at w0 erodes
        // (negligible at 50 Hz/10 kHz, real for high harmonics / low fs). Skipped
        // when the resonance is at/above Nyquist (tan undefined past pi/2).
        // PRResult::to_discrete_ss() uses the same frequency.
        const T w0_d = design::detail::pr_prewarped_w0(w0, Ts);

        // Tustin of R(s), s = k*(z-1)/(z+1), multiplied through by (z+1)².
        // Non-ideal: R(s) = 2*Ki*wc*s / (s² + 2*wc*s + w0_d²)
        // Ideal:     R(s) = Ki*s / (s² + w0_d²)
        const T a0 = k * k + T{2} * wc * k + w0_d * w0_d;
        const T a1 = T{2} * (w0_d * w0_d - k * k);
        const T a2 = k * k - T{2} * wc * k + w0_d * w0_d;
        const T b_scale = (wc > T{0}) ? (T{2} * Ki * wc * k) : (Ki * k);
        const T b0 = b_scale;
        const T b2 = -b_scale;

        resonant.set_coefficients(design::detail::normalize_biquad<T>(b0, T{0}, b2, a0, a1, a2));
    }
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<typename T>
PRController(const design::PRResult<T>&) -> PRController<T>;

} // namespace damp
