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


#include "damp/backend.hpp"
#include "damp/controllers/detail/anti_windup.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/core.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/toolbox/bounds.hpp"

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

/**
 * @brief β = tan(ω₀Ts/2): the pre-warped resonant frequency of pr_prewarped_w0 over the Tustin factor 2/Ts
 *
 * What PRController builds its coefficients from: the frequency of @ref pr_prewarped_w0 without
 * the factor 2/Ts, whose square would swamp it in single precision. At or above Nyquist it is
 * ω₀Ts/2 itself, the frequency left unwarped there.
 */
template<typename T>
[[nodiscard]] constexpr T pr_tan_half_angle(T w0, T Ts) {
    const T half = (w0 * Ts) / T{2};
    return (half < (damp::numbers::pi_v<T> / T{2})) ? damp::tan(half) : half;
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
    T Kp{};            ///< Proportional gain (finite)
    T Ki{};            ///< Resonant (integral) gain (finite)
    T w0{};            ///< Resonant frequency [rad/s], e.g. 2π·50 for a 50 Hz grid (> 0)
    T wc{};            ///< Resonant-term bandwidth [rad/s] (≥ 0); 0 = ideal PR, 5–15 rad/s typical non-ideal
    T tracking_time{}; ///< Back-calculation tracking time T_t [s] (≥ 0): ΔR = min(Ts/T_t, 1)(u_sat − (u + u_ff)); 0: the resonant section holds against a clamp
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
    T            tracking_time{};                     ///< Back-calculation tracking time T_t [s]; 0: the resonant section holds against a clamp
    bool         success{false};                      ///< status == DesignStatus::ok
    DesignStatus status{DesignStatus::not_converged}; ///< ok, or why the design failed
    Bounds<1, T> u_lim{};                             ///< Limits on the command (default unbounded): set here, or on the runtime

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return PRResult<U>{
            static_cast<U>(Kp), static_cast<U>(Ki), static_cast<U>(w0), static_cast<U>(wc),
            static_cast<U>(Ts), static_cast<U>(tracking_time), success, status, u_lim.template as<U>()
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
 *             tracking_time ≥ 0 (0 = no back-calculation)
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
                 && damp::finite_non_negative(spec.wc) && damp::finite_non_negative(spec.tracking_time);
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
        .tracking_time = spec.tracking_time,
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
 *                    of the fundamental, wc and tracking_time shared by every harmonic
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
                .tracking_time = fundamental.tracking_time,
            },
            Ts
        );
    }
    return results;
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Resonant section of PRController, run on its output y and the output's increment d
 *
 * The Tustin map of the resonant term is a second-order section with zeros at z = ±1,
 * @f[
 *   \frac{Y(z)}{X(z)} = \frac{\bar b_0\,(1 - z^{-2})}{1 + \bar a_1 z^{-1} + \bar a_2 z^{-2}},
 *   \qquad \bar a_1 = \kappa + c_2 - 2, \quad \bar a_2 = 1 - c_2 .
 * @f]
 * It is held by κ = 1 + ā₁ + ā₂ and c₂ = 1 − ā₂ and run on the output and its increment
 * @f$d[n] = y[n] - y[n-1]@f$:
 * @f[
 *   d[n] = d[n-1] + \bigl(\bar b_0\,(x[n] - x[n-2]) - c_2\, d[n-1] - \kappa\, y[n-1]\bigr),
 *   \qquad y[n] = y[n-1] + d[n] .
 * @f]
 * Substituting d gives the direct form
 * @f$y[n] = -\bar a_1 y[n-1] - \bar a_2 y[n-2] + \bar b_0 (x[n] - x[n-2])@f$ term for term, so
 * the transfer function and the meaning of the state (the last output, and through d the one
 * before it) are those of a Direct Form I section; only the rounding differs. It is the
 * section written around z = 1, an accumulator in place of each delay.
 *
 * What it guards against. With θ = ω₀Ts and a narrow resonance, the poles sit next to z = 1:
 * ā₁ ≈ −2 + θ² + 2ω_cTs and ā₂ ≈ 1 − 2ω_cTs. Single precision stores ā₁ to 1.2e-7, and θ² is
 * 9.9e-8 at 1 Hz and 20 kHz, so the direct form cannot hold the resonance: a float Direct Form I
 * section at 1 Hz and 20 kHz with ω_c = 5 rad/s ran to Inf in 36 s of a 0.01 sine at ω₀, and at
 * 50 Hz and 100 kHz with ω_c = 1 rad/s its gain at ω₀ was 0.78 of the design's. Here κ ≈ θ² and
 * c₂ ≈ 2ω_cTs are stored with their own relative precision, and each product in the recursion
 * is small next to the state it corrects, so its rounding is too.
 *
 * Limits. The resonance is placed to a few ulps of ω₀ and its peak is ω_c/(1 + tan²(θ/2)) wide
 * (the bilinear map narrows it toward Nyquist), so the gain and phase at ω₀ are the design's to
 * a few ε·Q, Q = ω₀(1 + tan²(θ/2))/(2ω_c): a float resonance with Q above about 1e6 is not
 * resolved. Roundoff in the recursion is about ε/√c₂ of the output, and the ideal form
 * (c₂ = 0), which is marginally stable by design, accumulates it. Next to Nyquist the poles
 * approach z = −1, where this form has the direct form's cancellation: the stability margin
 * 4 − κ − 2c₂ is 4/(1 + ω_cTs + tan²(θ/2)) and falls below rounding within about fs·√ε/(2π)
 * of Nyquist (5.5e-5·fs in float).
 *
 * @see "New Recursive Digital Filter Structures Having Very Low Sensitivity and Roundoff Noise"
 *      (Agarwal & Burrus, 1975), IEEE Trans. Circuits and Systems 22(12): second-order sections
 *      for poles near z = 1, whose coefficient sensitivity and roundoff noise do not grow with
 *      the sampling rate
 * @see "Improved Finite Word Length Characteristics in Digital Control Using Delta Operators"
 *      (Middleton & Goodwin, 1986), IEEE TAC 31(11): the delta operator, the increment in place
 *      of the shift, for better coefficient representation and less roundoff noise
 *
 * @tparam T Scalar type (default: float)
 */
template<typename T = float>
class PRResonator {
public:
    constexpr PRResonator() = default;

    /** @brief Convert across scalar precision, preserving coefficients and state. */
    template<typename U>
    constexpr explicit PRResonator(const PRResonator<U>& o)
        : b0_(static_cast<T>(o.b0_)),
          c2_(static_cast<T>(o.c2_)),
          kappa_(static_cast<T>(o.kappa_)),
          x1_(static_cast<T>(o.x1_)),
          x2_(static_cast<T>(o.x2_)),
          y_(static_cast<T>(o.y_)),
          d_(static_cast<T>(o.d_)) {}

    /** @brief Process one sample: d += b̄₀·(x[n] − x[n−2]) − c₂·d − κ·y, then y += d. */
    constexpr T operator()(T x) {
        // The three small terms are summed first, so d takes one rounding at its own scale. Behind
        // a barrier: -ffast-math regroups d − c₂·d into (1 − c₂)·d, the direct form's ā₂ again
        // (GCC 14 with FMA: 1.3e-3 off the gain at ω₀ at 5 Hz and 100 kHz, ω_c = 1 rad/s).
        d_ += detail::assoc_barrier(((b0_ * (x - x2_)) - (c2_ * d_)) - (kappa_ * y_));
        y_ += d_;
        x2_ = x1_;
        x1_ = x;
        return y_;
    }

    /** @brief Clear the inputs, the output and its increment; the coefficients are kept. */
    constexpr void reset() {
        x1_ = x2_ = y_ = d_ = T{0};
    }

    /**
     * @brief Replace the coefficients in place, keeping the state, so a retune keeps the output and its increment
     *
     * @param b0    b̄₀, the gain on x[n] − x[n−2]
     * @param c2    c₂ = 1 − ā₂ (≥ 0; 0 = undamped)
     * @param kappa κ = 1 + ā₁ + ā₂ (> 0)
     */
    constexpr void set_coefficients(T b0, T c2, T kappa) {
        b0_ = b0;
        c2_ = c2;
        kappa_ = kappa;
    }

    /** @brief Most recent output y[n−1]. */
    [[nodiscard]] constexpr T last_output() const { return y_; }

    /** @brief Add @p dy to the most recent output y[n−1] and leave y[n−2] in place: the increment moves with it. */
    constexpr void shift_last_output(T dy) {
        y_ += dy;
        d_ += dy;
    }

    /**
     * @brief Take back what the last sample added to the output, d[n−1], as far as a clamp cut the
     *        command off (@p excess = u_sat − u_asked), and no further (detail::hold_against_clamp)
     * @return How far the output moved
     */
    [[nodiscard]] constexpr T hold_against_clamp(T excess) {
        const T held = y_;
        detail::hold_against_clamp(y_, d_, excess);
        return y_ - held;
    }

private:
    template<typename>
    friend class PRResonator;

    T b0_{0};    ///< b̄₀; with c₂ = κ = 1 the section is inert (y = 0) until designed
    T c2_{1};    ///< c₂ = 1 − ā₂
    T kappa_{1}; ///< κ = 1 + ā₁ + ā₂
    T x1_{0};    ///< x[n−1]
    T x2_{0};    ///< x[n−2]
    T y_{0};     ///< y[n−1]
    T d_{0};     ///< d[n−1] = y[n−1] − y[n−2]
};

/**
 * @ingroup discrete_controllers
 * @brief Discrete Proportional-Resonant Controller
 *
 * PR controller realized as a proportional path plus a Tustin-discretized
 * resonant section, @ref PRResonator. The resonant frequency is pre-warped so the bilinear
 * map lands the digital peak exactly on w0, preserving zero steady-state error
 * at the target frequency (the unwarped map places the peak at the
 * bilinear-warped frequency, worst for high harmonics / low fs).
 *
 * Non-ideal (@f$\omega_c > 0@f$): @f$R(s) = 2 K_i \omega_c s/(s^2 + 2\omega_c s + \omega_0^2)@f$.
 * Ideal (@f$\omega_c = 0@f$): @f$R(s) = K_i s/(s^2 + \omega_0^2)@f$.
 * Both use @f$s = (2/T_s)(z-1)/(z+1)@f$.
 *
 * A failed or continuous (Ts = 0) design does not compile under constant evaluation; at run time
 * it leaves the controller invalid (valid() false, output 0). So does a u_lim that is not a box.
 *
 * The command is clamped into u_lim. Against that clamp, and one reported with back_calculate, the
 * resonant section gives back what the tick added to its output, as far as the clamp cut the
 * command off (conditional integration, tracking_time = 0), or unwinds by back-calculation with the tracking
 * time T_t, tracking_time.
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
    T            Kp{};
    T            Ki{};
    T            w0{};
    T            wc{};
    T            Ts{};
    T            tracking_time{T{0}}; ///< Back-calculation tracking time T_t [s]: ΔR = min(Ts/T_t, 1)(u_sat − (u + u_ff)); 0: the resonant section holds
    Bounds<1, T> u_lim{};             ///< Limits on the command (the result's; default unbounded)
    T            u_cmd{T{0}};         ///< Last command control() returned: the u in back_calculate's u_sat − (u + u_ff)
    bool         valid_{false};       ///< From the design's success flag; gates control()

    //! Discretized resonant term (inert until designed).
    PRResonator<T> resonant{};

    constexpr PRController() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr PRController(const design::PRResult<U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : PRController(result.template as<T>()) {}

    constexpr PRController(const design::PRResult<T>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : Kp(result.Kp),
          Ki(result.Ki),
          w0(result.w0),
          wc(result.wc),
          Ts(result.Ts),
          tracking_time(result.tracking_time),
          u_lim(result.u_lim),
          valid_(result.success && result.Ts > T{0} && result.u_lim.valid()) {
        if (result.success && !result.u_lim.valid() && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(design::DesignStatus::invalid_parameter); // a limit that is not finite, or lower above upper
        }
        if ((!result.success || !(result.Ts > T{0})) && damp::is_constant_evaluated()) {
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
          tracking_time(static_cast<T>(other.tracking_time)),
          u_lim(other.u_lim.template as<T>()),
          u_cmd(static_cast<T>(other.u_cmd)),
          valid_(other.valid_),
          resonant(other.resonant) {}

    /**
     * @brief Compute control output
     *
     * @param error Current error (reference - measurement)
     * @return Control output u, clamped into u_lim (zero if the design was invalid)
     */
    [[nodiscard]] constexpr T control(T error) {
        if (!valid_) {
            return T{0};
        }
        const T u = (Kp * error) + resonant(error);
        u_cmd = u_lim.saturate(u);
        static_cast<void>(wind_back(u_cmd - u));
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
     * @brief Anti-windup after a downstream clamp: the resonant section holds (tracking_time = 0), or
     *        ΔR = min(Ts/T_t, 1)·(u_sat − (u + u_ff))
     *
     * With tracking_time = 0 the resonant output gives back what the last tick added to it, as far as the
     * clamp cut the command off, and no further, as against u_lim. With tracking_time > 0, a damped unwind:
     * the excess of the applied command over this law's last command u (plus any feedforward
     * @p u_ff added downstream), times Ts / T_t, is added to the most recent resonant output, so
     * the next tick's resonant contribution is pulled back toward the realizable command. No Ki
     * factor: the resonant output is already in command units, so T_t is a time [s] (the
     * PID laws' Kbc is a gain ratio instead). It takes at most the whole excess per tick: past
     * that it would overshoot, and past twice it diverges.
     *
     * @param u_sat Command actually applied after downstream clamping.
     * @param u_ff  Feedforward added to this law's command before the clamp (default 0).
     */
    constexpr void back_calculate(T u_sat, T u_ff = T{0}) {
        if (valid_) {
            static_cast<void>(wind_back(u_sat - (u_cmd + u_ff)));
        }
    }

    /**
     * @brief Take a clamp's @p excess (u_sat − u_asked) back from the resonant section: hold
     *        (tracking_time = 0) or a share min(Ts/T_t, 1) of it
     * @return How far the resonant output moved
     */
    [[nodiscard]] constexpr T wind_back(T excess) {
        if (tracking_time == T{0}) {
            return resonant.hold_against_clamp(excess);
        }
        const T shift = damp::min(Ts / tracking_time, T{1}) * excess;
        resonant.shift_last_output(shift);
        return shift;
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
    /**
     * @brief Coefficients of the resonant section, each a ratio of positive terms
     *
     * Tustin of R(s), s = k(z−1)/(z+1) with k = 2/Ts, at the pre-warped frequency
     * @f$\omega_0' = k\beta@f$, @f$\beta = \tan(\omega_0 T_s/2)@f$, so that the digital peak sits
     * exactly on ω₀ (@ref design::PRResult::to_discrete_ss uses the same frequency). Divided
     * through by k², with @f$\alpha = \omega_c T_s/2@f$, the denominator is
     * @f$(z-1)^2 + 2\alpha(z^2-1) + \beta^2(z+1)^2@f$, and normalizing by
     * @f$\mathrm{den} = 1 + 2\alpha + \beta^2@f$ gives the coefficients of @ref PRResonator:
     * @f[
     *   \kappa = \frac{4\beta^2}{\mathrm{den}}, \qquad c_2 = \frac{4\alpha}{\mathrm{den}}, \qquad
     *   \bar b_0 = \frac{2 K_i \alpha}{\mathrm{den}} = \frac{K_i c_2}{2}
     *   \quad\Bigl(\text{ideal, } \alpha = 0:\ \bar b_0 = \frac{K_i T_s}{2\,\mathrm{den}}\Bigr).
     * @f]
     * No k² is formed and nothing cancels: in float the sum k² + 2ω_c·k + ω₀'² drops ω₀'² once
     * it is below the ulp of k² (39 against 128 at 1 Hz and 20 kHz), and −2 + κ + c₂ cannot be
     * stored either (see @ref PRResonator). κ and c₂ are within a few ulps of their values. The
     * poles have |z|² = 1 − c₂ and are inside the unit circle while κ + 2c₂ < 4, which holds by
     * 4/den. A den that is not positive and finite (a hand-built result) leaves the section inert.
     */
    constexpr void compute_coefficients() {
        const T alpha = (wc * Ts) / T{2};
        const T beta = design::detail::pr_tan_half_angle(w0, Ts);
        const T beta2 = beta * beta;
        const T den = T{1} + (T{2} * alpha) + beta2;
        if (!damp::finite_positive(den)) {
            resonant.set_coefficients(T{0}, T{1}, T{1});
            return;
        }
        const T c2 = T{4} * (alpha / den);
        const T kappa = T{4} * (beta2 / den);
        const T b0 = (wc > T{0}) ? (Ki * (c2 / T{2})) : (Ki * (Ts / (T{2} * den)));
        resonant.set_coefficients(b0, c2, kappa);
    }
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<typename T>
PRController(const design::PRResult<T>&) -> PRController<T>;

} // namespace damp
