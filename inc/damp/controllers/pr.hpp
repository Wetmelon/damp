// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file pr.hpp
 * @brief The proportional-resonant law: design::pr(gains, Ts, opts...) → design::PRResult →
 *        PRController, one resonance or a bank of them at harmonics of a fundamental
 *
 * @defgroup pr_controller Proportional-Resonant Controller
 * @brief Proportional-resonant control of AC references and disturbances, at one frequency or at
 *        a sparse set of harmonics (design::Harmonics)
 *
 * PR controllers provide infinite gain at a specific resonant frequency, achieving zero
 * steady-state error for sinusoidal references. They are the AC equivalent of PI controllers for
 * DC references.
 *
 * Common applications:
 * - Grid-tied inverter current control
 * - Active power filter harmonic compensation (a bank: Harmonics{{1, 5, 7, 11, 13}})
 * - Rotating reference frame alternatives (no Park transform needed)
 * - Motor torque-ripple cancellation
 *
 * Transfer function (non-ideal PR, @f$\omega_c > 0@f$):
 * @f[
 *   C(s) = K_p + \frac{2 K_i \omega_c\, s}{s^2 + 2 \omega_c s + \omega_0^2}
 * @f]
 * where @f$K_p@f$ is the proportional gain, @f$K_i@f$ the resonant (integral) gain,
 * @f$\omega_0@f$ the resonant frequency (e.g. @f$2\pi\cdot 50@f$ for a 50 Hz grid) and
 * @f$\omega_c@f$ the resonant-term bandwidth.
 *
 * Ideal PR (@f$\omega_c = 0@f$):
 * @f[
 *   C(s) = K_p + \frac{K_i s}{s^2 + \omega_0^2}
 * @f]
 * Infinite gain exactly at @f$\omega_0@f$; fragile to grid-frequency drift, so non-ideal PR with
 * a few rad/s of @f$\omega_c@f$ is preferred in practice.
 *
 * A bank, `design::pr(gains, Ts, design::Harmonics{{1, 5, 7}})`, adds one resonant term per
 * harmonic h, of gain @f$K_i/h@f$ at @f$h\omega_0@f$, with the proportional term once:
 * @f$C(s) = K_p + \sum_h R_h(s)@f$. It targets a specific, sparse set of harmonics with
 * independent resonators. Repetitive control (controllers/repetitive.hpp) rejects every harmonic
 * of a period with one delay loop instead.
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
 * @see filters/spectral.hpp (Goertzel, HarmonicAnalyzer) for detecting which harmonics to target
 */


#include "damp/backend.hpp"
#include "damp/controllers/detail/pr_resonator.hpp"
#include "damp/controllers/law.hpp"     // IWYU pragma: export
#include "damp/design/input_limits.hpp" // IWYU pragma: export
#include "damp/design/options.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/systems/zpk.hpp"
#include "damp/toolbox/bounds.hpp"

namespace damp {

namespace detail {

/**
 * @brief The coefficients of every resonator of a PR bank: the h-th harmonic's resonant term,
 *        Kᵢ/h at h·ω₀, for each order h of @p orders, at the period @p Ts
 * @return The coefficients, or nullopt when Kp or Kᵢ is not finite, ω_c or the tracking time is
 *         below 0 or not finite, ω₀ or Ts is not finite and positive, an order is 0, a resonance
 *         is at or above Nyquist, or a coefficient overflows
 */
template<typename T, size_t NH>
[[nodiscard]] constexpr damp::optional<damp::array<PRCoefficients<T>, NH>>
pr_bank_coefficients(T Kp, T Ki, T w0, T wc, T tracking_time, T Ts, const damp::array<size_t, NH>& orders) {
    if (!damp::isfinite(Kp) || !damp::isfinite(Ki) || !damp::finite_positive(w0) || !damp::finite_non_negative(wc)
        || !damp::finite_non_negative(tracking_time) || !damp::finite_positive(Ts)) {
        return damp::nullopt;
    }
    damp::array<PRCoefficients<T>, NH> out{};
    for (size_t i = 0; i < NH; ++i) {
        if (orders[i] == 0) {
            return damp::nullopt; // the 0th harmonic (DC) is not a resonant target
        }
        const T order = static_cast<T>(orders[i]);
        if (!pr_resonance_ok(order * w0, Ts)) {
            return damp::nullopt;
        }
        out[i] = pr_coefficients(Ki / order, order * w0, wc, Ts);
        if (!out[i].finite()) {
            return damp::nullopt;
        }
    }
    return out;
}

} // namespace detail

namespace design {

namespace detail {
using ::damp::detail::pr_bank_coefficients;
using ::damp::detail::pr_prewarped_w0;
} // namespace detail

/**
 * @struct PRGains
 * @brief The gains and frequency of a proportional-resonant law: what design::pr takes, the result
 *        holds and the runtime retunes (PRController::gains, PRController::set_gains)
 *
 * The caller names every number: `design::pr({.Kp = 1.0, .Ki = 200.0, .w0 = w, .wc = 10.0}, Ts)`.
 * With design::Harmonics among the options they describe the fundamental: each harmonic h gets
 * Kᵢ/h at h·ω₀ (@ref harmonic), and Kp is applied once.
 */
template<typename T = double>
struct PRGains {
    T Kp{};            ///< Proportional gain (finite)
    T Ki{};            ///< Resonant (integral) gain (finite): the resonant term's gain at ω₀
    T w0{};            ///< [rad/s] Resonant frequency, the fundamental's in a bank, e.g. 2π·50 for a 50 Hz grid (> 0)
    T wc{};            ///< [rad/s] Resonant-term bandwidth (≥ 0); 0 = ideal PR, 5–15 rad/s typical non-ideal
    T tracking_time{}; ///< [s] Back-calculation tracking time T_t (≥ 0): ΔR = min(Ts/T_t, 1)(u_sat − (u + u_ff)); 0: the resonant terms hold against a clamp

    /** @brief The same gains in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr PRGains<U> as() const {
        return PRGains<U>{static_cast<U>(Kp), static_cast<U>(Ki), static_cast<U>(w0), static_cast<U>(wc), static_cast<U>(tracking_time)};
    }

    /** @brief The resonant term at the harmonic @p h of ω₀: Kᵢ/h at h·ω₀, no proportional term (the zero law for h = 0). */
    [[nodiscard]] constexpr PRGains harmonic(size_t h) const {
        if (h == 0) {
            return PRGains{};
        }
        const T order = static_cast<T>(h);
        return PRGains{.Kp = T{0}, .Ki = Ki / order, .w0 = w0 * order, .wc = wc, .tracking_time = tracking_time};
    }

    /**
     * @brief Continuous-time transfer function C(s) = Kp + R(s), highest power of s first
     *
     * Non-ideal (ω_c > 0): num = {Kp, 2ω_c(Kp + Kᵢ), Kp·ω₀²}, den = {1, 2ω_c, ω₀²}.
     * Ideal (ω_c = 0): num = {Kp, Kᵢ, Kp·ω₀²}, den = {1, 0, ω₀²}.
     */
    [[nodiscard]] constexpr TransferFunction<3, 3, T> to_tf() const {
        if (wc > T{0}) {
            return TransferFunction<3, 3, T>{
                .num = {Kp, T{2} * wc * (Kp + Ki), Kp * w0 * w0},
                .den = {T{1}, T{2} * wc, w0 * w0},
            };
        }
        return TransferFunction<3, 3, T>{
            .num = {Kp, Ki, Kp * w0 * w0},
            .den = {T{1}, T{0}, w0 * w0},
        };
    }

    /**
     * @brief Continuous-time state-space of C(s) (2nd order SISO)
     *
     * Controllable canonical form; D = Kp carries the proportional feedthrough. The resonant C
     * row is [0, 2Kᵢω_c] (non-ideal) or [0, Kᵢ] (ideal).
     */
    [[nodiscard]] constexpr StateSpace<2, 1, 1, T> to_ss() const {
        const T c_res = (wc > T{0}) ? (T{2} * Ki * wc) : Ki;
        return StateSpace{
            .A = Matrix<2, 2, T>{{T{0}, T{1}}, {-w0 * w0, -T{2} * wc}},
            .B = Matrix<2, 1, T>{{T{0}}, {T{1}}},
            .C = Matrix<1, 2, T>{{T{0}, c_res}},
            .D = Matrix<1, 1, T>{{Kp}},
        };
    }
};

/**
 * @struct Harmonics
 * @brief The harmonics of ω₀ a PR law resonates at, among design::pr's options: `Harmonics{{1, 5, 7}}`
 *
 * One resonant term per order h listed, of gain Kᵢ/h at h·ω₀, run in the order listed; the law's
 * proportional term is applied once. Left out, the law has one resonance, at ω₀. An order of 0 (DC
 * is not a resonant target) or a resonance at or above Nyquist is refused with
 * DesignStatus::invalid_parameter. The list is an array, so its length is deduced from the braces.
 *
 * @tparam N The resonators: the orders listed (≥ 1)
 */
template<size_t N>
struct Harmonics {
    static_assert(N >= 1, "Harmonics lists at least one order: Harmonics{{1, 5, 7}}.");
    static constexpr size_t listed = N;  ///< The orders listed: the law's resonators
    using option_kind = kind::Harmonics; ///< Its kind of design option
    size_t orders[N]{};                  ///< The harmonic orders h ≥ 1, the fundamental 1
};

/**
 * @struct PRLimits
 * @brief The limits of a PR law: on its command
 *
 * Unbounded by default. design::pr sets the ones a law starts with from its options (InputLimits),
 * and the law reads and changes them at run time with PRController::limits and
 * PRController::set_limits.
 */
template<typename T = double>
struct PRLimits {
    Bounds<1, T> input{}; ///< The command u [lower, upper]

    /** @brief The same limits in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr PRLimits<U> as() const {
        return PRLimits<U>{input.template as<U>()};
    }
};

/**
 * @struct PRResult
 * @brief A proportional-resonant design deployed at a sample period: the gains, the harmonics it
 *        resonates at, the limits the law starts with, and T_s
 *
 * Produced by design::pr. The runtime, PRController, discretizes each resonant term by Tustin at
 * its pre-warped frequency (@ref to_discrete_ss is the same model).
 *
 * @tparam T  Scalar
 * @tparam NH Resonators (1 without design::Harmonics)
 */
template<typename T = double, size_t NH = 1>
struct [[nodiscard]] PRResult {
    static_assert(NH >= 1, "A PR law has at least one resonator.");
    static constexpr size_t resonators = NH; ///< The resonant terms

    PRGains<T>              gains{};                             ///< The gains, the fundamental's in a bank
    damp::array<size_t, NH> harmonics{};                         ///< The order h of each resonator, Kᵢ/h at h·ω₀ ({1} without Harmonics)
    PRLimits<T>             limits{};                            ///< The limits the law starts with
    T                       Ts{};                                ///< [s] Sample period of the fixed-rate tick
    bool                    success{false};                      ///< status == DesignStatus::ok
    DesignStatus            status{DesignStatus::not_converged}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The same design in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr PRResult<U, NH> as() const {
        return PRResult<U, NH>{gains.template as<U>(), harmonics, limits.template as<U>(), static_cast<U>(Ts), success, status};
    }

    /** @brief The continuous law of resonator @p i: its harmonic's resonant term, and Kp on the first. */
    [[nodiscard]] constexpr PRGains<T> resonator(size_t i) const {
        PRGains<T> g = gains.harmonic(harmonics[i]);
        g.Kp = (i == 0) ? gains.Kp : T{0};
        return g;
    }

    /**
     * @brief Continuous-time state-space of the law, Kp + Σ R_h(s): the resonators in parallel
     *
     * State dimension 2·NH, resonator i in states 2i and 2i + 1. A failed design is the zero system.
     */
    [[nodiscard]] constexpr StateSpace<2 * NH, 1, 1, T> to_ss() const {
        if (!success) {
            return {};
        }
        return bank<NH - 1>([this](size_t i) { return resonator(i).to_ss(); });
    }

    /** @brief Continuous-time transfer function of the law: closed form for one resonator, of @ref to_ss for a bank. */
    [[nodiscard]] constexpr TransferFunction<(2 * NH) + 1, (2 * NH) + 1, T> to_tf() const {
        if constexpr (NH == 1) {
            return success ? resonator(0).to_tf() : TransferFunction<3, 3, T>{};
        } else {
            return ::damp::to_transfer_function(to_ss());
        }
    }

    /**
     * @brief Discrete-time state-space of the PRController runtime (Tustin, each ω pre-warped)
     *
     * Discretizes the way @ref PRController does: each resonator is rebuilt at its pre-warped
     * frequency @f$\omega' = (2/T_s)\tan(\omega T_s/2)@f$ and mapped with
     * @f$s = (2/T_s)(z-1)/(z+1)@f$, so each digital resonance sits exactly on its ω and the
     * model's response equals the runtime's tick for tick. For another discretization use
     * `discretize(to_ss(), Ts, method)`. A failed design is the zero system.
     */
    [[nodiscard]] constexpr StateSpace<2 * NH, 1, 1, T> to_discrete_ss() const {
        if (!success) {
            return {};
        }
        return bank<NH - 1>([this](size_t i) {
            PRGains<T> warped = resonator(i);
            warped.w0 = detail::pr_prewarped_w0(warped.w0, Ts);
            return *discretize(warped.to_ss(), Ts, DiscretizationMethod::Tustin);
        });
    }

private:
    /** @brief The resonators 0 … I of @p section in parallel. */
    template<size_t I, typename Section>
    [[nodiscard]] static constexpr StateSpace<2 * (I + 1), 1, 1, T> bank(const Section& section) {
        if constexpr (I == 0) {
            return section(0);
        } else {
            return *::damp::parallel(bank<I - 1>(section), section(I));
        }
    }
};

namespace detail {

/** @brief The options design::pr takes. */
using pr_options = option_set<kind::InputLimits, kind::Harmonics>;

/** @brief The resonators a PR design's options give: as many as a Harmonics lists, else 1. */
template<typename... Opts>
inline constexpr size_t pr_resonators = (option_count<kind::Harmonics, Opts...> > 0) ? listed_by<kind::Harmonics, Opts...> : size_t{1};

/** @brief The compile-time check of design::pr's options: their kinds, one of each, and one input. */
template<typename... Opts>
consteval void check_pr_options() {
    pr_options::check<Opts...>();
    if constexpr (option_count<kind::InputLimits, Opts...> > 0) {
        static_assert(option_t<kind::InputLimits, void, Opts...>::inputs == 1, "A PR law has one input: InputLimits{.lower = {lo}, .upper = {hi}}, one entry at each end.");
    }
}

/** @brief The harmonic orders among @p opts: a Harmonics' list, else the fundamental alone. */
template<size_t NH, typename... Opts>
[[nodiscard]] constexpr damp::array<size_t, NH> pr_harmonics_of(const Opts&... opts) {
    damp::array<size_t, NH> out{};
    if constexpr (option_count<kind::Harmonics, Opts...> > 0) {
        const auto given = the<kind::Harmonics>(opts...);
        for (size_t i = 0; i < NH; ++i) {
            out[i] = given.orders[i];
        }
    } else {
        (static_cast<void>(opts), ...);
        out[0] = 1;
    }
    return out;
}

/** @brief The InputLimits among @p opts in the scalar @p T (none, or one of more inputs, which check_pr_options refuses: unbounded, with no second message). */
template<typename T, typename... Opts>
[[nodiscard]] constexpr Bounds<1, T> pr_input_limits(const Opts&... opts) {
    if constexpr (option_t<kind::InputLimits, InputLimits<1, T>, Opts...>::inputs == 1) {
        return input_limits_of<1, T>(opts...); // unbounded without one
    } else {
        (static_cast<void>(opts), ...);
        return Bounds<1, T>{};
    }
}

} // namespace detail

/**
 * @brief A proportional-resonant law from its gains, deployed at @p Ts
 *
 * The scalar follows @p Ts, so a braced spec needs no template argument:
 * `design::pr({.Kp = 1.0, .Ki = 200.0, .w0 = w, .wc = 10.0}, 1e-4)`. The options trail in any
 * order, at most one of each kind: InputLimits{.lower, .upper} (the command, one input) and
 * Harmonics{{h, …}} (a bank at harmonics of ω₀).
 *
 * @code
 * // 50 Hz grid current regulator at 10 kHz, clamped to the bus
 * constexpr auto kCurrent = design::pr({.Kp = 1.0, .Ki = 200.0, .w0 = 2.0 * pi * 50.0, .wc = 10.0}, 1e-4,
 *                                      design::InputLimits{.lower = {-400.0}, .upper = {400.0}});
 * // The fundamental and the 5th and 7th line harmonics, Kᵢ/h at h·ω₀
 * constexpr auto kBank = design::pr({.Kp = 0.5, .Ki = 200.0, .w0 = 2.0 * pi * 50.0, .wc = 5.0}, 1e-4,
 *                                   design::Harmonics{{1, 5, 7}});
 * @endcode
 *
 * @param gains Kp and Kᵢ finite, ω₀ > 0, ω_c ≥ 0 (0 = ideal PR), tracking_time ≥ 0 (0 = the
 *              resonant terms hold against a clamp)
 * @param Ts    [s] The law's sample period (finite and positive)
 * @param opts  InputLimits, Harmonics
 * @return The design, for PRController or Law; refused, with zero gains, with invalid_sample_time
 *         for a Ts that is not finite and positive, and with invalid_parameter for a gain out of
 *         its domain, an order of 0, a resonance at or above Nyquist π/Ts, or a limit that is not a box
 */
template<typename T = double, typename... Opts>
[[nodiscard]] constexpr PRResult<T, detail::pr_resonators<Opts...>> pr(const PRGains<T>& gains, T Ts, const Opts&... opts) {
    detail::check_pr_options<Opts...>();
    using Result = PRResult<T, detail::pr_resonators<Opts...>>;
    if (!damp::finite_positive(Ts)) {
        return detail::refused<Result>(DesignStatus::invalid_sample_time);
    }
    Result out{
        .gains = gains,
        .harmonics = detail::pr_harmonics_of<Result::resonators>(opts...),
        .limits = {.input = detail::pr_input_limits<T>(opts...)},
        .Ts = Ts,
    };
    if (!out.limits.input.valid() || !detail::pr_bank_coefficients(gains.Kp, gains.Ki, gains.w0, gains.wc, gains.tracking_time, Ts, out.harmonics)) {
        return detail::refused<Result>(DesignStatus::invalid_parameter);
    }
    detail::set_status(out, DesignStatus::ok);
    return out;
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Proportional-resonant runtime: Kp·e plus NH resonant terms, its gains, frequency and
 *        limits settable at run time
 *
 * A proportional path plus a Tustin-discretized resonant section per resonance (detail::PRResonator),
 * each at its pre-warped frequency so the bilinear map lands the digital peak exactly on it,
 * preserving zero steady-state error there (the unwarped map places the peak at the
 * bilinear-warped frequency, worst for high harmonics / low fs). NH = 1 is the PR law at ω₀; a bank
 * (design::Harmonics) runs one resonator per harmonic h, of gain Kᵢ/h at h·ω₀.
 *
 * Non-ideal (@f$\omega_c > 0@f$): @f$R(s) = 2 K_i \omega_c s/(s^2 + 2\omega_c s + \omega_0^2)@f$.
 * Ideal (@f$\omega_c = 0@f$): @f$R(s) = K_i s/(s^2 + \omega_0^2)@f$.
 * Both use @f$s = (2/T_s)(z-1)/(z+1)@f$.
 *
 * | Member | What it is |
 * | --- | --- |
 * | step(e), step(r, y) | one tick at the design's T_s on the error e = r − y |
 * | track(u, e), track(u, r, y) | in place of each step while another source drives the plant |
 * | back_calculate(u_sat, u_ff) | a clamp after the law: the input the plant received |
 * | gains(), set_gains(g), set_frequency(w0) | the gains in force (design::PRGains); change them, or ω₀ alone, at run time |
 * | limits(), set_limits(lim) | the command's limits (design::PRLimits); change them at run time |
 * | resonant(), harmonics(), sample_time() | the resonant terms' sum in the last command; the orders; T_s |
 * | last_command(), valid(), reset() | the last command; whether the design ran; clear the resonators |
 *
 * A failed design, a T_s that is not positive, a gain out of its domain or a limit that is not a
 * box is refused at construction: ill-formed under constant evaluation (`constinit` deploy); at
 * run time it leaves the controller invalid (valid() false, output 0).
 *
 * The command is clamped into limits().input. Against that clamp, and one reported with
 * back_calculate, the resonators give back what the tick added to their outputs, in the order
 * listed, until the command meets the clamp (conditional integration, tracking_time = 0), or unwind
 * by back-calculation with the tracking time T_t, a share min(T_s/T_t, 1) of the excess: never more
 * in all than the clamp cut off.
 *
 * Gains, ω₀ and limits are tunable at run time and take effect at the next tick; each setter
 * refuses a value out of its domain and changes nothing. A retune keeps each resonator's output
 * and its increment, so the command does not jump:
 * @code
 * if (!pr.set_frequency(w_grid)) { } // a PLL's estimate: refused at or above Nyquist, nothing changed
 * auto lim = pr.limits();
 * lim.input = {.lower = -v_bus, .upper = v_bus};
 * if (!pr.set_limits(lim)) { }
 * @endcode
 *
 * @code
 * static constinit PRController<float> pr{design::pr({.Kp = 1.0, .Ki = 200.0, .w0 = 2.0 * pi * 50.0, .wc = 10.0}, 1e-4)};
 * float u = pr.step(i_ref, i_meas); // (r, y) form
 * @endcode
 *
 * @tparam T  Scalar type (default: float)
 * @tparam NH Resonators (1: the PR law at ω₀)
 */
template<typename T = float, size_t NH = 1>
class PRController {
    static_assert(NH >= 1, "A PR law has at least one resonator.");

    template<typename, size_t>
    friend class PRController;

public:
    constexpr PRController() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr PRController(const design::PRResult<U, NH>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : PRController(result.template as<T>()) {}

    /** @brief From a design result. Non-explicit for assign-from-design: `PRController<float> pr = design::pr(gains, Ts).as<float>();` */
    constexpr PRController(const design::PRResult<T, NH>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : gains_(result.gains), orders_(result.harmonics), u_lim_(result.limits.input), Ts_(result.Ts) {
        if (!result.success && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(result.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
        if (!damp::finite_positive(result.Ts) && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(design::DesignStatus::invalid_sample_time);
        }
        valid_ = result.success && damp::finite_positive(result.Ts) && result.limits.input.valid() && retune(result.gains);
        if (!valid_ && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(design::DesignStatus::invalid_parameter); // a gain out of its domain, or a limit that is not a box
        }
    }

    /** @brief The same law, with its live state, in another scalar. */
    template<typename U>
    constexpr explicit PRController(const PRController<U, NH>& o)
        : gains_(o.gains_.template as<T>()),
          orders_(o.orders_),
          u_lim_(o.u_lim_.template as<T>()),
          Ts_(static_cast<T>(o.Ts_)),
          unwind_(static_cast<T>(o.unwind_)),
          u_cmd_(static_cast<T>(o.u_cmd_)),
          valid_(o.valid_) {
        for (size_t i = 0; i < NH; ++i) {
            res_[i] = detail::PRResonator<T>(o.res_[i]);
        }
    }

    /** @brief False when built from a design it refuses: step() returns 0. */
    [[nodiscard]] constexpr bool valid() const { return valid_; }

    /** @brief One tick on the error @p e = r − y: Kp·e plus the resonant terms, clamped into limits().input (0 if invalid). */
    [[nodiscard]] constexpr T step(T e) {
        if (!valid_) {
            return T{0};
        }
        T u = gains_.Kp * e;
        for (size_t i = 0; i < NH; ++i) {
            u += res_[i](e);
        }
        u_cmd_ = u_lim_.saturate(u);
        wind_back(u_cmd_ - u);
        return u_cmd_;
    }

    /** @brief One tick on the reference @p r and the measurement @p y: step(r − y). */
    [[nodiscard]] constexpr T step(T r, T y) { return step(r - y); }

    /**
     * @brief In place of step(@p e) while another source drives the plant: the memory that holds @p u, the input applied
     *
     * The resonators run on e, so their input memory is current, then the first one's output is
     * shifted so that this tick's command, Kp·e plus the resonant terms, is u clamped into the input
     * limits: the resonant memory is the applied input's history less the proportional term. @p u
     * is recorded as the command applied, so a back_calculate of it finds no excess.
     *
     * A resonator holds no constant (it has zeros at z = ±1), so this is not a preload of the next
     * command: the step after a track continues from that memory. For a u at the resonance, the
     * resonator's own oscillation, it is the input's next sample; for a constant u it is u short
     * by κ·(u − Kp·e) on the first step, κ = 4β²/(1 + 2α + β²) ≈ (ω₀T_s)², then the resonance rings
     * down at ω_c. With Kᵢ = 0 the law is a P and there is nothing to preload.
     */
    constexpr void track(T u, T e) {
        if (!valid_) {
            return;
        }
        T v = gains_.Kp * e;
        for (size_t i = 0; i < NH; ++i) {
            v += res_[i](e);
        }
        if (gains_.Ki != T{0}) {
            res_[0].shift_last_output(u_lim_.saturate(u) - v);
        }
        u_cmd_ = u;
    }

    /** @brief In place of step(@p r, @p y): track(@p u, @p r − @p y). */
    constexpr void track(T u, T r, T y) { track(u, r - y); }

    /**
     * @brief Anti-windup after a downstream clamp: the resonant terms hold against the excess
     *        u_sat − (u + u_ff), or with a tracking time unwind by min(T_s/T_t, 1) of it
     *
     * Call after the stage that clamps the summed command, in the same tick as step(): @p u_sat is
     * what the plant received, @p u_ff whatever was added to this law's command before the clamp, u
     * the command step() returned (@ref last_command). With tracking_time = 0 each resonator gives
     * back what the last tick added to its output, as far as the clamp cut the command off, and no
     * further, as against limits().input. With tracking_time > 0, a damped unwind: the excess times
     * min(T_s/T_t, 1) is added to the most recent resonant output, so the next tick's resonant
     * contribution is pulled back toward the realizable command. No Kᵢ factor: the resonant output
     * is already in command units, so T_t is a time [s], as the PID laws' tracking_time is. It takes
     * at most the whole excess per tick: past that it would overshoot, and past twice it diverges.
     *
     * @param u_sat Command actually applied after downstream clamping.
     * @param u_ff  Feedforward added to this law's command before the clamp (default 0).
     */
    constexpr void back_calculate(T u_sat, T u_ff = T{0}) {
        if (valid_) {
            wind_back(u_sat - (u_cmd_ + u_ff));
        }
    }

    /** @brief The command the last step() returned, or the input the last track() recorded. */
    [[nodiscard]] constexpr T last_command() const { return u_cmd_; }

    /** @brief Clear the resonators' memory and the last command; keeps the gains and the limits. */
    constexpr void reset() {
        for (size_t i = 0; i < NH; ++i) {
            res_[i].reset();
        }
        u_cmd_ = T{0};
    }

    /** @brief The gains in force: the design's, or the last set_gains or set_frequency. */
    [[nodiscard]] constexpr design::PRGains<T> gains() const { return gains_; }

    /**
     * @brief Change the gains at run time; the next tick uses them
     *
     * Recomputes every resonator's coefficients at sample_time() and keeps its output and
     * increment, so the command does not jump. Refused (false, nothing changed) for a Kp or Kᵢ that
     * is not finite, an ω_c or a tracking time below 0 or not finite, an ω₀ that is not finite and
     * positive or that puts a resonance at or above Nyquist, or a coefficient that overflows.
     */
    [[nodiscard]] constexpr bool set_gains(const design::PRGains<T>& g) { return retune(g); }

    /**
     * @brief Move the resonance (in a bank, the fundamental and every harmonic with it) to @p w0
     *        [rad/s]: grid-frequency adaptation
     *
     * set_gains with ω₀ alone changed. Refused (false, nothing changed) for an @p w0 that is not
     * finite and positive, or one that would put a resonance at or above Nyquist: the bank is
     * retuned whole or not at all.
     */
    [[nodiscard]] constexpr bool set_frequency(T w0) {
        design::PRGains<T> g = gains_;
        g.w0 = w0;
        return retune(g);
    }

    /** @brief The limits in force: the design's, or the last set_limits. */
    [[nodiscard]] constexpr design::PRLimits<T> limits() const { return design::PRLimits<T>{u_lim_}; }

    /**
     * @brief Change the limits at run time; the next tick clamps into them
     *
     * Refused (false, nothing changed) when an end is NaN or ±Inf, or the lower end is above the upper one.
     */
    [[nodiscard]] constexpr bool set_limits(const design::PRLimits<T>& lim) {
        if (!lim.input.valid()) {
            return false;
        }
        u_lim_ = Bounds<1, T>(lim.input);
        return true;
    }

    /** @brief The resonant terms' sum in the last command (u's units): the resonators' most recent outputs. */
    [[nodiscard]] constexpr T resonant() const {
        T sum{0};
        for (size_t i = 0; i < NH; ++i) {
            sum += res_[i].last_output();
        }
        return sum;
    }

    /** @brief The harmonic order of each resonator, Kᵢ/h at h·ω₀. */
    [[nodiscard]] constexpr const damp::array<size_t, NH>& harmonics() const { return orders_; }

    /** @brief [s] The tick's sample period T_s. */
    [[nodiscard]] constexpr T sample_time() const { return Ts_; }

private:
    /** @brief Take @p g at sample_time() when every resonator's coefficients are in their domain, keeping each resonator's state; nothing changed otherwise. */
    [[nodiscard]] constexpr bool retune(const design::PRGains<T>& g) {
        const auto k = detail::pr_bank_coefficients(g.Kp, g.Ki, g.w0, g.wc, g.tracking_time, Ts_, orders_);
        if (!k) {
            return false;
        }
        for (size_t i = 0; i < NH; ++i) {
            res_[i].set_coefficients((*k)[i]);
        }
        gains_ = g;
        unwind_ = (g.tracking_time > T{0}) ? damp::min(Ts_ / g.tracking_time, T{1}) : T{0};
        return true;
    }

    /** @brief Take a clamp's @p excess (u_sat − u_asked) back from the resonators in turn, each from what is left: hold (tracking_time = 0) or a share min(T_s/T_t, 1). */
    constexpr void wind_back(T excess) {
        for (size_t i = 0; i < NH; ++i) {
            if (gains_.tracking_time == T{0}) {
                excess -= res_[i].hold_against_clamp(excess);
            } else {
                const T shift = unwind_ * excess;
                res_[i].shift_last_output(shift);
                excess -= shift;
            }
        }
    }

    damp::array<detail::PRResonator<T>, NH> res_{};        ///< One resonant section per harmonic (inert until designed)
    design::PRGains<T>                      gains_{};      ///< The gains in force
    damp::array<size_t, NH>                 orders_{};     ///< The harmonic order of each resonator
    Bounds<1, T>                            u_lim_{};      ///< The command's limits
    T                                       Ts_{T{0}};     ///< [s] The tick's period
    T                                       unwind_{T{0}}; ///< min(T_s/T_t, 1), 0 for T_t = 0 (hold)
    T                                       u_cmd_{T{0}};  ///< Last command step() returned, or the input track() recorded
    bool                                    valid_{false}; ///< False when built from a design it refuses: step() returns 0
};

/** @brief Deduce the runtime from its design, keeping the result's scalar and resonators. */
template<typename T, size_t NH>
PRController(const design::PRResult<T, NH>&) -> PRController<T, NH>;

namespace detail {

/** @brief Whether @p R is a PR design result. */
template<typename R>
inline constexpr bool is_pr_result = false;

/** @brief design::pr's result. */
template<typename U, size_t NH>
inline constexpr bool is_pr_result<design::PRResult<U, NH>> = true;

/** @brief A PR design (design::pr). */
template<const auto& Res>
concept pr_design = is_pr_result<design_t<Res>>;

} // namespace detail

/**
 * @brief Law<kRes, T> of a PR design: the PRController it runs as, with the design's resonators
 *
 * A `constexpr` design that failed, or has no positive T_s, does not compile wherever the Law is
 * declared. The gains stay members, settable at run time.
 */
// ponytail: a constexpr design's coefficients could be constants of the type, as PID's
// ConstantPIDGains are; PRResonator would read them from a gain source in place of members.
template<const auto& Res, typename T>
    requires detail::pr_design<Res>
struct Law<Res, T> : PRController<T, detail::design_t<Res>::resonators> {
    using Runtime = PRController<T, detail::design_t<Res>::resonators>; ///< The runtime the design runs as

    static_assert(detail::constant_design_checked<Res>());

    /** @brief The design's law. */
    constexpr Law() : Runtime{Res} {}

    using Runtime::operator=;
};

} // namespace damp
