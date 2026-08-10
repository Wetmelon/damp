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
 * @brief Proportional-resonant (and multi-harmonic) controllers for AC tracking
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

#include "damp/backend.hpp"
#include "damp/filters/filters.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/core.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

namespace damp {

namespace design {

/**
 * @struct PRResult
 * @brief Proportional-Resonant controller design result
 */
template<typename T = double>
struct PRResult {
    T    Kp{};           ///< Proportional gain
    T    Ki{};           ///< Resonant (integral) gain
    T    w0{};           ///< Resonant frequency (rad/s)
    T    wc{};           ///< Cutoff bandwidth of resonant term (rad/s); 0 = ideal PR
    T    Ts{};           ///< Sampling time (seconds); 0 = continuous-time design
    bool success{false}; ///< true if w0 > 0, wc ≥ 0, Ts ≥ 0, and all parameters finite

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return PRResult<U>{
            static_cast<U>(Kp), static_cast<U>(Ki), static_cast<U>(w0), static_cast<U>(wc), static_cast<U>(Ts), success
        };
    }

    /**
     * @brief Convert to continuous-time transfer function
     *
     * Non-ideal (@f$\omega_c > 0@f$):
     *   C(s) = Kp + 2*Ki*wc*s / (s² + 2*wc*s + w0²)
     *   num = {Kp*w0², 2*wc*(Kp+Ki), Kp}, den = {w0², 2*wc, 1}
     *
     * Ideal (@f$\omega_c = 0@f$):
     *   C(s) = Kp + Ki*s / (s² + w0²)
     *   num = {Kp*w0², Ki, Kp}, den = {w0², 0, 1}
     *
     * Numerator and denominator in ascending powers of s.
     */
    [[nodiscard]] constexpr TransferFunction<3, 3, T> to_tf() const {
        if (wc > T{0}) {
            return TransferFunction<3, 3, T>{
                .num = {Kp * w0 * w0, T{2} * wc * (Kp + Ki), Kp},
                .den = {w0 * w0, T{2} * wc, T{1}},
            };
        }
        // Ideal PR: C(s) = Kp + Ki*s/(s² + w0²)
        return TransferFunction<3, 3, T>{
            .num = {Kp * w0 * w0, Ki, Kp},
            .den = {w0 * w0, T{0}, T{1}},
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
     * @brief Convert to discrete-time state-space
     *
     * @param method Discretization method (default: Tustin)
     */
    [[nodiscard]] constexpr StateSpace<2, 1, 1, T>
    to_discrete_ss(DiscretizationMethod method = DiscretizationMethod::Tustin) const {
        return *discretize(to_ss(), Ts, method);
    }
};

/**
 * @brief Design a Proportional-Resonant controller
 *
 * @param Kp  Proportional gain (must be finite)
 * @param Ki  Resonant gain (must be finite)
 * @param w0  Resonant frequency (rad/s, e.g. 2*pi*50 for 50Hz); must be > 0
 * @param wc  Cutoff bandwidth of resonant term (rad/s). Use 0 for ideal PR.
 *            Typical values: 5–15 rad/s for non-ideal PR. Must be ≥ 0.
 * @param Ts  Sampling time (seconds); 0 = continuous-time design. Must be ≥ 0.
 * @return PRResult with design parameters; success=false if any constraint fails
 */
template<typename T = double>
[[nodiscard]] constexpr PRResult<T> pr(T Kp, T Ki, T w0, T wc, T Ts) {
    const bool ok = damp::isfinite(Kp) && damp::isfinite(Ki) && damp::isfinite(w0) && damp::isfinite(wc)
                 && damp::isfinite(Ts) && (w0 > T{0}) && (wc >= T{0}) && (Ts >= T{0});
    if (!ok) {
        return PRResult<T>{};
    }
    return PRResult<T>{
        .Kp = Kp,
        .Ki = Ki,
        .w0 = w0,
        .wc = wc,
        .Ts = Ts,
        .success = true,
    };
}

/**
 * @brief Design multiple-harmonic PR controller gains
 *
 * For harmonic compensation, returns an array of PRResult for harmonics
 * 1, 3, 5, 7, ... (or user-specified harmonic numbers).
 *
 * @param Kp          Proportional gain (shared across all harmonics)
 * @param Ki_fund     Resonant gain for fundamental
 * @param w_fund      Fundamental frequency (rad/s)
 * @param wc          Cutoff bandwidth (rad/s)
 * @param Ts          Sampling time
 * @param harmonics   Array of harmonic numbers (e.g. {1, 3, 5, 7})
 * @return Array of PRResult, one per harmonic
 */
template<size_t N, typename T = double>
[[nodiscard]] constexpr damp::array<PRResult<T>, N>
pr_harmonics(T Kp, T Ki_fund, T w_fund, T wc, T Ts, const damp::array<size_t, N>& harmonics) {
    damp::array<PRResult<T>, N> results{};
    for (size_t i = 0; i < N; ++i) {
        T w_h = w_fund * static_cast<T>(harmonics[i]);
        // Typically reduce Ki for higher harmonics
        T Ki_h = Ki_fund / static_cast<T>(harmonics[i]);
        results[i] = pr((i == 0) ? Kp : T{0}, Ki_h, w_h, wc, Ts);
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
 * @code
 * // 50 Hz grid current regulator at 10 kHz
 * auto                d = design::pr(1.0, 200.0, 2 * std::numbers::pi * 50, 10.0, 1e-4);
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
    T    Kbc{T{0}};     ///< Back-calculation coefficient (same law as PID); larger = slower unwind (0 = unit-rate fallback)
    bool valid_{false}; ///< From the design's success flag; gates control()

    //! Discretized resonant term (Direct Form I biquad, inert until designed).
    Biquad<T> resonant{design::SecondOrderCoeffs<T>{}};

    constexpr PRController() = default;

    constexpr PRController(const design::PRResult<T>& result)
        : Kp(result.Kp), Ki(result.Ki), w0(result.w0), wc(result.wc), Ts(result.Ts), valid_(result.success) {
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
          Kbc(static_cast<T>(other.Kbc)),
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
        return Kp * error + resonant(error);
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

    /// Clear the resonant filter's delay line; gains and config are preserved.
    constexpr void reset() {
        resonant.reset();
    }

    [[nodiscard]] constexpr bool valid() const { return valid_; }

    /**
     * @brief Anti-windup hook for cascade-level saturation propagation.
     *
     * Damped unwind: adds `(u_sat - u_unsat) * Ts / Kbc` to the most recent
     * resonant output so the next tick's resonant contribution is pulled back
     * toward the realizable command. Same shape and sign convention as
     * `ContinuousPID::back_calculate` (Ts per call). Discrete
     * `PIDController::back_calculate` bakes Ts into the integrator units and
     * takes only (u_unsat, u_sat). No-op when `Kbc == 0`.
     *
     * @param u_unsat Command this controller requested.
     * @param u_sat   Command actually applied after downstream clamping.
     * @param sample_time Sample time (s).
     */
    constexpr void back_calculate(T u_unsat, T u_sat, T sample_time) {
        if (!valid_ || Kbc == T{0} || u_unsat == u_sat) {
            return;
        }
        resonant.set_last_output(resonant.last_output() + sample_time * ((u_sat - u_unsat) / Kbc));
    }

    /**
     * @brief Update resonant frequency (for grid frequency adaptation)
     *
     * @param new_w0 New resonant frequency (rad/s)
     */
    constexpr void set_frequency(T new_w0) {
        if (!(new_w0 > T{0}) || !damp::isfinite(new_w0)) {
            return;
        }
        w0 = new_w0;
        if (valid_) {
            compute_coefficients();
        }
    }

private:
    constexpr void compute_coefficients() {
        // Ts <= 0 means a continuous-time design (e.g. a design factory's default
        // Ts = 0) was handed to a discrete controller — no valid Tustin mapping.
        // Zero the resonant term so control() falls back to pure proportional.
        if (!(Ts > T{0})) {
            resonant.set_coefficients(design::SecondOrderCoeffs<T>{});
            return;
        }

        const T k = T{2} / Ts; // Tustin substitution factor 2/Ts

        // Pre-warp the resonant frequency so the bilinear map places the digital
        // peak exactly on w0 — w_pw = (2/Ts)*tan(w0*Ts/2). Without this the peak
        // lands at the warped frequency and zero-steady-state error at w0 erodes
        // (negligible at 50 Hz/10 kHz, real for high harmonics / low fs). Skip
        // when the resonance is at/above Nyquist (tan undefined past pi/2).
        T       w0_d = w0;
        const T half = w0 * Ts / T{2};
        if (half < damp::numbers::pi_v<T> / T{2}) {
            w0_d = k * damp::tan(half);
        }

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

/**
 * @ingroup discrete_controllers
 * @brief Multi-harmonic PR Controller
 *
 * Combines proportional gain with multiple resonant terms for harmonic
 * compensation. Each resonant term tracks one harmonic frequency.
 *
 * @tparam N Number of harmonic terms
 * @tparam T Scalar type (default: float)
 */
template<size_t N, typename T = float>
struct MultiPRController {
    T                               Kp{};
    damp::array<PRController<T>, N> resonants{};
    bool                            valid_{false};

    constexpr MultiPRController() = default;

    template<size_t M>
    constexpr MultiPRController(const damp::array<design::PRResult<T>, M>& results)
        requires(M == N)
        : Kp(results[0].Kp), valid_(results[0].success) {
        for (size_t i = 0; i < N; ++i) {
            auto r = results[i];
            r.Kp = T{0}; // Kp is applied once, not per-harmonic
            resonants[i] = PRController<T>(r);
            valid_ = valid_ && r.success;
        }
    }

    template<typename U>
    constexpr explicit MultiPRController(const MultiPRController<N, U>& other)
        : Kp(static_cast<T>(other.Kp)), valid_(other.valid_) {
        for (size_t i = 0; i < N; ++i) {
            resonants[i] = PRController<T>(other.resonants[i]);
        }
    }

    /**
     * @brief Compute control output (sum of all resonant terms + Kp)
     */
    [[nodiscard]] constexpr T control(T error) {
        if (!valid_) {
            return T{0};
        }
        T u = Kp * error;
        for (size_t i = 0; i < N; ++i) {
            u += resonants[i].control(error);
        }
        return u;
    }

    /**
     * @brief Reference-tracking overload satisfying SISOController.
     */
    [[nodiscard]] constexpr T control(T r, T y) {
        return control(r - y);
    }

    constexpr void reset() {
        for (size_t i = 0; i < N; ++i) {
            resonants[i].reset();
        }
    }

    [[nodiscard]] constexpr bool valid() const { return valid_; }
};

} // namespace damp
