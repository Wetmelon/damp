// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file adrc.hpp
 * @brief Active Disturbance Rejection Control: design::adrc<NX>(gains, Ts, opts...) → ADRCResult
 *        → ADRCController
 *
 * Linear ADRC: an extended-state observer (ESO) estimates the plant output, its derivatives
 * and the total disturbance f on the model y⁽ⁿ⁾ = f + b₀u, and a PD law on the estimate cancels
 * f̂. The observer is @ref design::eso / @ref ESO; this header is that observer plus the PD law.
 * The runtime's ESO is the discrete one of Miklosovic, Radke and Gao: the ESO's linear
 * dynamics discretized by ZOH, in current-estimator form, with its poles at the image
 * @f$e^{-\omega_o T_s}@f$ of the continuous design's poles at @f$-\omega_o@f$.
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/adrc_gains.hpp"
#include "damp/controllers/law.hpp" // IWYU pragma: export
#include "damp/design/detail/siso_input_limits.hpp"
#include "damp/design/input_limits.hpp" // IWYU pragma: export
#include "damp/design/status.hpp"
#include "damp/estimation/eso.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/systems/zpk.hpp"
#include "damp/toolbox/bounds.hpp"

namespace damp {

namespace design {

/**
 * @struct ADRCGains
 * @brief The tuning of a linear ADRC law: the controller bandwidth ω_c, the observer bandwidth
 *        ω_o and the plant's gain b₀
 *
 * design::adrc's input, `design::adrc<2>({.wc = 10.0, .wo = 50.0, .b0 = 1.0}, Ts)`, the result's
 * gains, and what ADRCController::gains and set_gains read and take.
 *
 * @tparam T Scalar type (deduced from wc when the type is named)
 */
template<typename T = double>
struct ADRCGains {
    T                        wc{}; ///< [rad/s] Controller bandwidth ω_c (> 0)
    damp::type_identity_t<T> wo{}; ///< [rad/s] Observer bandwidth ω_o (> 0), as design::eso's
    damp::type_identity_t<T> b0{}; ///< Plant gain b₀ in y⁽ⁿ⁾ = f + b₀u (> 0)

    /** @brief The same gains in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr ADRCGains<U> as() const {
        return ADRCGains<U>{static_cast<U>(wc), static_cast<U>(wo), static_cast<U>(b0)};
    }
};

/**
 * @struct ADRCLimits
 * @brief The limits of an ADRC law: on its command
 *
 * Unbounded by default. design::adrc sets them from an InputLimits among its options, and the law
 * reads and changes them at run time with ADRCController::limits and set_limits.
 */
template<typename T = double>
struct ADRCLimits {
    Bounds<1, T> input{}; ///< The command u [lower, upper]: the ESO runs on the clamped command

    /** @brief The same limits in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr ADRCLimits<U> as() const {
        return ADRCLimits<U>{input.template as<U>()};
    }
};

/**
 * @struct ADRCResult
 * @brief Active Disturbance Rejection Control design result
 *
 * The gains, the observer and controller gains they give, and the discrete ESO the runtime
 * ticks, z̄_k = Φ ẑ_{k−1} + Γ u_{k−1}, ẑ_k = z̄_k + L (y_k − z̄_k,₀). The continuous ESO gains β
 * are for analysis (to_ss, to_tf). Use .as\<U\>() to convert for type conversion (e.g., double
 * to float).
 *
 * @see "From PID to Active Disturbance Rejection Control" (Han, 1998)
 * @see "Discrete implementation and generalization of the extended state observer"
 *      (Miklosovic, Radke & Gao, 2006), ACC
 */
template<size_t NX, typename T = double>
struct [[nodiscard]] ADRCResult {
    static_assert(NX == 1 || NX == 2, "ADRC supports 1st- and 2nd-order plants only");

    ADRCGains<T> gains{}; ///< ω_c, ω_o, b₀

    ColVec<NX + 1, T> beta{}; ///< Continuous ESO gains (poles at −ω_o)

    T Kp{}; ///< Unscaled proportional gain (wc or wc²); tick divides by b0
    T Kd{}; ///< Unscaled derivative gain (0 or 2·wc); tick divides by b0

    Matrix<NX + 1, NX + 1, T> Phi{};   ///< ZOH state map Φ of the ESO
    ColVec<NX + 1, T>         Gamma{}; ///< ZOH input map Γ on u_{k−1}
    ColVec<NX + 1, T>         L{};     ///< Current-estimator gain on the innovation y_k − z̄_k,₀

    ADRCLimits<T> limits{};                            ///< The limits the law starts with
    T             Ts{0};                               ///< [s] Sample period of the tick
    bool          success{false};                      ///< status == DesignStatus::ok
    DesignStatus  status{DesignStatus::not_converged}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The same design in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr ADRCResult<NX, U> as() const {
        return ADRCResult<NX, U>{
            .gains = gains.template as<U>(),
            .beta = beta.template as<U>(),
            .Kp = static_cast<U>(Kp),
            .Kd = static_cast<U>(Kd),
            .Phi = Phi.template as<U>(),
            .Gamma = ColVec<NX + 1, U>(Gamma.template as<U>()),
            .L = ColVec<NX + 1, U>(L.template as<U>()),
            .limits = limits.template as<U>(),
            .Ts = static_cast<U>(Ts),
            .success = success,
            .status = status,
        };
    }

    /** @brief The observer this law ticks: @ref design::eso of the same ω_o, b₀ and T_s. */
    [[nodiscard]] constexpr ESOResult<NX, T> as_eso() const {
        return {
            .wo = gains.wo,
            .b0 = gains.b0,
            .beta = beta,
            .Phi = Phi,
            .Gamma = Gamma,
            .L = L,
            .success = success,
            .status = status,
            .Ts = Ts,
        };
    }

    /**
     * @brief Linear ADRC as a 2-DOF compensator: u from [r; y]
     *
     * Continuous ESO + PD (unsaturated). State is the ESO [ẑ₁, (ẑ₂,) f̂].
     * Failed designs return a zero system.
     */
    [[nodiscard]] constexpr StateSpace<NX + 1, 2, 1, T> to_ss_2dof() const {
        StateSpace<NX + 1, 2, 1, T> sys{};
        const T                     b0 = gains.b0;
        if (!success || !damp::finite_positive(b0)) {
            return sys;
        }
        if constexpr (NX == 1) {
            sys.A(0, 0) = -(beta[0] + Kp);
            sys.A(1, 0) = -beta[1];
            sys.B(0, 0) = Kp;
            sys.B(0, 1) = beta[0];
            sys.B(1, 1) = beta[1];
            sys.C(0, 0) = -Kp / b0;
            sys.C(0, 1) = -T{1} / b0;
            sys.D(0, 0) = Kp / b0;
        } else {
            sys.A(0, 0) = -beta[0];
            sys.A(0, 1) = T{1};
            sys.A(1, 0) = -(beta[1] + Kp);
            sys.A(1, 1) = -Kd;
            sys.A(2, 0) = -beta[2];
            sys.B(1, 0) = Kp;
            sys.B(0, 1) = beta[0];
            sys.B(1, 1) = beta[1];
            sys.B(2, 1) = beta[2];
            sys.C(0, 0) = -Kp / b0;
            sys.C(0, 1) = -Kd / b0;
            sys.C(0, 2) = -T{1} / b0;
            sys.D(0, 0) = Kp / b0;
        }
        return sys;
    }

    /**
     * @brief SISO loop map C(s) = −U/Y at r = 0
     *
     * Unity-negative-feedback convention (same as PID / lead / PR `to_ss`):
     * u = C(s)(r − y). Bode / `series` / `feedback` consume this map.
     */
    [[nodiscard]] constexpr StateSpace<NX + 1, 1, 1, T> to_ss() const {
        const auto                  full = to_ss_2dof();
        StateSpace<NX + 1, 1, 1, T> sys{};
        sys.A = full.A;
        for (size_t i = 0; i < NX + 1; ++i) {
            sys.B(i, 0) = full.B(i, 1);
        }
        for (size_t j = 0; j < NX + 1; ++j) {
            sys.C(0, j) = -full.C(0, j);
        }
        sys.D(0, 0) = -full.D(0, 1);
        return sys;
    }

    /**
     * @brief Transfer function of @ref to_ss (Leverrier)
     *
     * NX = 1: C(s) = (K_p β₂ + (K_p β₁ + β₂)s) / (b₀ s (s + β₁ + K_p)).
     */
    [[nodiscard]] constexpr TransferFunction<NX + 2, NX + 2, T> to_tf() const {
        return ::damp::to_transfer_function(to_ss());
    }
};

/**
 * @brief Active Disturbance Rejection Control design, deployed at @p Ts
 *
 * The observer is @ref design::eso (poles at −ω_o, plant gain b₀, the discrete ESO at Ts).
 * Controller gains are unscaled (b0 is applied only in the runtime law u = (u0 − f̂)/b0):
 * - 1st order: Kp = wc, Kd = 0
 * - 2nd order: Kp = wc², Kd = 2·wc
 *
 * Only plant orders NX = 1 and NX = 2 are supported (same surface as typical
 * linear ADRC toolboxes / MATLAB® discrete ADRC). Higher-order ADRC is out of
 * scope.
 *
 * @code
 * constexpr auto kAdrc = design::adrc<2>({.wc = 10.0, .wo = 50.0, .b0 = 1.0}, 1e-3,
 *                                       design::InputLimits{.lower = {-2.5}, .upper = {2.5}});
 * static constinit Law<kAdrc, float> adrc; // each tick: u = adrc.step(r, y);
 * @endcode
 *
 * Refused with invalid_sample_time for a Ts that is not finite and positive, and with
 * invalid_parameter for ω_c, ω_o or b₀ not finite and positive or an InputLimits that is not a box
 * (with zero gains, so the runtime cannot divide by an invalid b0), and as design::eso refuses its
 * observer. An option other than one InputLimits on one input does not compile.
 *
 * @param gains ω_c, ω_o and b₀
 * @param Ts    [s] Sample period of the runtime (> 0)
 * @param opts  InputLimits{.lower = {lo}, .upper = {hi}}
 *
 * @return ADRCResult with computed gains
 *
 * @see design::eso — the observer this law ticks
 * @see "From PID to Active Disturbance Rejection Control" (Han, 1998)
 * @see "Discrete implementation and generalization of the extended state observer"
 *      (Miklosovic, Radke & Gao, 2006), ACC: the current discrete ESO and its gains
 */
template<size_t NX, typename T, typename... Opts>
[[nodiscard]] constexpr ADRCResult<NX, T> adrc(const ADRCGains<T>& gains, T Ts, const Opts&... opts) {
    static_assert(NX == 1 || NX == 2, "ADRC supports 1st- and 2nd-order plants only");
    detail::check_siso_input_limits<Opts...>();
    const ADRCLimits<T> limits{detail::siso_input_limits_of<T>(opts...)};
    if (!damp::finite_positive(Ts)) {
        return detail::refused<ADRCResult<NX, T>>(DesignStatus::invalid_sample_time);
    }
    const T    wc = gains.wc;
    const bool ok = damp::finite_positive(wc) && damp::finite_positive(gains.wo) && damp::finite_positive(gains.b0) && limits.input.valid();
    if (!ok) {
        return detail::refused<ADRCResult<NX, T>>(DesignStatus::invalid_parameter);
    }

    const auto e = eso<NX>(ESOGains<T>{.wo = gains.wo, .b0 = gains.b0}, Ts);
    if (!e.success) {
        return detail::refused<ADRCResult<NX, T>>(e.status);
    }

    return ADRCResult<NX, T>{
        .gains = gains,
        .beta = e.beta,
        .Kp = (NX == 1) ? wc : (wc * wc),
        .Kd = (NX == 1) ? T{0} : (T{2} * wc),
        .Phi = e.Phi,
        .Gamma = e.Gamma,
        .L = e.L,
        .limits = limits,
        .Ts = e.Ts,
        .success = true,
        .status = DesignStatus::ok,
    };
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Active Disturbance Rejection Control (ADRC)
 *
 * Discrete ADRC with an extended-state observer (ESO) for SISO systems.
 * Supports both 1st-order and 2nd-order plant models.
 * Based on the work of Jingqing Han and Zhiqiang Gao. The observer is @ref ESO
 * (design::eso): one mat-vec per tick; this runtime adds the PD law that cancels f̂.
 *
 * | Member | What it is |
 * | --- | --- |
 * | step(r, y) | one tick at the design's T_s: the ESO, then the PD law on its estimate |
 * | track(u, r, y) | in place of step while another source drives the plant: the ESO runs on the input applied |
 * | back_calculate(u_sat, u_ff) | a clamp after the law: the ESO's next tick takes the input applied |
 * | gains(), set_gains(g) | ω_c, ω_o, b₀ in force (design::ADRCGains); change them at run time |
 * | limits(), set_limits(lim) | the command's limits (design::ADRCLimits); change them at run time |
 * | observer(), applied(), last_command(), sample_time(), valid(), reset() | the ESO state; the input it takes as applied; the last command; T_s; whether the design ran; clear the ESO |
 *
 * While another source drives the plant, track(u, r, y) takes the place of step(r, y): the ESO
 * ticks on y with the input applied over the last period, then holds u as applied over the next,
 * so its estimate follows the plant under the other source. Nothing is preloaded: f̂ stays the
 * estimate of the total disturbance. The first step after it returns the law's own command on
 * that current estimate, sat((K_p(r − ẑ₁) − K_d ẑ₂ − f̂)/b₀), with no kick from a stale one; it
 * is u when the plant rests at r under u (f̂ = −b₀u there).
 *
 * References:
 * - Han, J. (1998). From PID to Active Disturbance Rejection Control. IEEE Transactions on Industrial Electronics, 45(5), 900-906.
 * - Gao, Z. (2006). Active Disturbance Rejection Control: A Paradigm Shift in Feedback Control System Design. American Control Conference.
 *
 * @see "Discrete implementation and generalization of the extended state observer"
 *      (Miklosovic, Radke & Gao, 2006), ACC
 * @tparam NX    Plant order: 1 (1st-order) or 2 (2nd-order) only
 * @tparam T     Scalar type (default: float)
 * @tparam Gains Where the gains and the ESO's maps live: void, the default, for members (detail::ADRCGainSet,
 *               settable), or the constants of a `constexpr` design's Law (detail::ConstantADRCGains)
 */
// The default is void, not detail::ADRCGainSet<NX, T>: a default that names the parameters crashed clang 23 deducing through an alias (pid.hpp).
template<size_t NX, typename T = float, typename Gains = void>
class ADRCController {
    static_assert(NX == 1 || NX == 2, "ADRC supports 1st- and 2nd-order plants only");
    using GainSet = damp::conditional_t<damp::is_same_v<Gains, void>, detail::ADRCGainSet<NX, T>, Gains>; ///< The gains' storage
    using EsoStore = typename GainSet::Eso;                                                               ///< Where the ESO's maps live
    static constexpr bool kHeld = !GainSet::constant;                                                     ///< The gains are members: set_gains changes them

    template<size_t, typename, typename>
    friend class ADRCController; // cross-precision converting ctor

    /** @brief The runtime's observer of @p result: on its matrices, or the type's constants. */
    [[nodiscard]] static constexpr ESO<NX, T, EsoStore> observer_of(const design::ADRCResult<NX, T>& result) {
        if constexpr (EsoStore::constant) {
            static_cast<void>(result);
            return ESO<NX, T, EsoStore>{};
        } else {
            return ESO<NX, T, EsoStore>{result.as_eso()};
        }
    }

    /** @brief The gains' storage for @p result: its gains, PD gains and T_s, or nothing on a constexpr design's Law. */
    [[nodiscard]] static constexpr GainSet gains_of(const design::ADRCResult<NX, T>& result) {
        if constexpr (kHeld) {
            return GainSet{.wc = result.gains.wc, .wo = result.gains.wo, .b0 = result.gains.b0, .Kp = result.Kp, .Kd = result.Kd, .Ts = result.Ts};
        } else {
            static_cast<void>(result);
            return GainSet{};
        }
    }

public:
    constexpr ADRCController() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr ADRCController(const design::ADRCResult<NX, U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : ADRCController(result.template as<T>()) {}

    /**
     * @brief From the design result, which carries the gains, the limits and the sample time.
     *
     * A failed design, one whose Ts is not finite and positive or whose limit is not a box is
     * ill-formed under constant evaluation; at run time it leaves the controller invalid
     * (valid() false, output 0).
     */
    constexpr ADRCController(const design::ADRCResult<NX, T>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : k_(gains_of(result)),
          u_lim_(result.limits.input),
          eso_(observer_of(result)),
          valid_(result.success && damp::finite_positive(result.Ts) && damp::finite_positive(result.gains.b0) && result.limits.input.valid()) {
        if ((!result.success || !damp::finite_positive(result.Ts)) && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(result.success ? design::DesignStatus::invalid_sample_time : result.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
        if (!valid_ && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(design::DesignStatus::invalid_parameter); // b₀ not finite and positive, or a limit that is not a box
        }
    }

    /** @brief The same law, with its live state, in another scalar (or off a Law's constant gains and ESO maps, onto members). */
    template<typename U, typename G>
        requires kHeld
    constexpr explicit ADRCController(const ADRCController<NX, U, G>& other)
        : k_{.wc = static_cast<T>(other.k_.wc),
             .wo = static_cast<T>(other.k_.wo),
             .b0 = static_cast<T>(other.k_.b0),
             .Kp = static_cast<T>(other.k_.Kp),
             .Kd = static_cast<T>(other.k_.Kd),
             .Ts = static_cast<T>(other.k_.Ts)},
          u_lim_(other.u_lim_.template as<T>()),
          eso_(other.eso_),
          u_cmd_(static_cast<T>(other.u_cmd_)),
          valid_(other.valid_) {}

    /**
     * @brief Reference-tracking control step (Gao's linear ADRC).
     *
     * Runs the discrete extended-state observer one step (predict with the previously-applied
     * command, correct with y_k), then forms the control law
     *
     *     u = (Kp·(r − ẑ₁) − Kd·ẑ₂ − f̂) / b0
     *
     * where `f̂ = z_{NX+1}` is the ESO's total-disturbance estimate. Cancelling
     * `f̂` is what gives ADRC its disturbance rejection and (integral-free) zero
     * steady-state error. The two design knobs are the controller and observer
     * bandwidths `wc`/`wo` (see @ref design::adrc); there is no separate
     * integrator gain.
     *
     * The command is clamped into limits().input, and the ESO takes that clamped command as the
     * input applied, so its disturbance estimate does not absorb what the limit cut off. If a
     * stage after this law clamps further, follow up with `back_calculate(u_sat)` so the next ESO
     * tick uses what was actually applied.
     *
     * A refused design emits zero and holds the ESO state, so the control law never divides by an
     * invalid b0 or runs without the observer.
     *
     * @param r  Reference (setpoint).
     * @param y  Measurement (plant output).
     * @return Control command `u`.
     */
    [[nodiscard]] constexpr T step(T r, T y) {
        if (!valid()) {
            return T{0};
        }

        const auto& z = eso_.step(y); // current estimator on the last applied command

        // --- Control law: PD on the estimated state, minus the disturbance.
        const T u0 = (k_.Kp * (r - z[0])) - (k_.Kd * z[1]); // Kd == 0 for NX == 1
        const T u = u_lim_.saturate((u0 - z[NX]) / k_.b0);  // cancel f̂ = z_{NX+1}

        eso_.back_calculate(u);
        u_cmd_ = u;
        return u;
    }

    /**
     * @brief In place of step(@p r, @p y) while another source drives the plant, @p u the input applied
     *
     * The ESO ticks on @p y with the input applied over the last period, then takes @p u as the
     * input applied over the next, so its estimate stays current. Nothing is preloaded: the first
     * step after it returns the law's own command on that estimate (see the class).
     * @p u is recorded as the command applied.
     */
    constexpr void track(T u, T r, T y) {
        if (!valid()) {
            return;
        }
        static_cast<void>(r);
        static_cast<void>(eso_.step(y));
        eso_.back_calculate(u);
        u_cmd_ = u;
    }

    /**
     * @brief Anti-windup hook for cascade-level saturation propagation.
     *
     * Overrides the stored command with the actually-realized value so
     * the next ESO tick reflects what the plant received, not what the
     * controller wanted. This is the natural ADRC anti-windup behavior:
     * the ESO's disturbance estimate self-corrects once it sees the true
     * input, without an explicit integrator unwind step.
     *
     * Same signature as the other laws' back_calculate. The ESO needs only the applied
     * input, and a feedforward added downstream is already part of @p u_sat, so @p u_ff
     * does not enter.
     *
     * @param u_sat Command actually applied after downstream clamping.
     * @param u_ff  Feedforward added to this law's command before the clamp (unused here).
     */
    constexpr void back_calculate(T u_sat, [[maybe_unused]] T u_ff = T{0}) { eso_.back_calculate(u_sat, u_ff); }

    /** @brief A successful design with b0 > 0 and a limit that is a box; gates step(). */
    [[nodiscard]] constexpr bool valid() const { return valid_ && eso_.valid(); }

    /** @brief ESO state [ŷ, (ŷ̇,) f̂] — plot these when commissioning. */
    [[nodiscard]] constexpr const ColVec<NX + 1, T>& observer() const { return eso_.state(); }

    /** @brief The command the last step() returned, or the input the last track() recorded. */
    [[nodiscard]] constexpr T last_command() const { return u_cmd_; }

    /** @brief The input the next ESO tick takes as applied: the last command, or what back_calculate recorded since. */
    [[nodiscard]] constexpr T applied() const { return eso_.last_command(); }

    /** @brief The gains in force: the design's, or the last set_gains. */
    [[nodiscard]] constexpr design::ADRCGains<T> gains() const { return design::ADRCGains<T>{.wc = k_.wc, .wo = k_.wo, .b0 = k_.b0}; }

    /**
     * @brief Change the gains at run time; the next tick uses them
     *
     * Recomputes the PD gains and the ESO's maps at sample_time() (design::adrc), and keeps the
     * estimate: a new b₀ moves f̂ by −Δb₀·u, u the input the ESO holds as applied, so the estimate
     * of y⁽ⁿ⁾ = f + b₀u and the next command at the same estimate do not jump. Refused (false,
     * nothing changed) for ω_c, ω_o or b₀ not finite and positive, or an observer design::eso
     * refuses at sample_time(). A `constexpr` design's Law, whose gains and ESO maps are the
     * type's, has no set_gains.
     */
    [[nodiscard]] constexpr bool set_gains(const design::ADRCGains<T>& g)
        requires kHeld
    {
        const auto res = design::adrc<NX>(g, k_.Ts);
        if (!res.success) {
            return false;
        }
        const T u = eso_.last_command();
        auto    z = eso_.state();
        z[NX] = z[NX] - ((res.gains.b0 - k_.b0) * u); // f = y⁽ⁿ⁾ − b₀u: the same y⁽ⁿ⁾ under the new b₀
        ESO<NX, T, EsoStore> retuned{res.as_eso()};
        retuned.reset(z);
        retuned.back_calculate(u);
        eso_ = retuned;
        k_ = gains_of(res);
        return true;
    }

    /** @brief The limits in force: the design's, or the last set_limits. */
    [[nodiscard]] constexpr design::ADRCLimits<T> limits() const { return {u_lim_}; }

    /** @brief Change the limits at run time; the next tick clamps into them. Refused (false, nothing changed) when an end is NaN or ±Inf, or the lower end is above the upper one. */
    [[nodiscard]] constexpr bool set_limits(const design::ADRCLimits<T>& lim) {
        if (!lim.input.valid()) {
            return false;
        }
        u_lim_ = lim.input;
        return true;
    }

    /** @brief [s] The tick's sample period T_s. */
    [[nodiscard]] constexpr T sample_time() const { return k_.Ts; }

    /** @brief Zero the ESO state and the held input. */
    constexpr void reset() {
        eso_.reset();
        u_cmd_ = T{0};
    }

    /** @brief Start the ESO from @p z0 with no held input: e.g. [y₀, 0, 0] avoids the first-command kick from ẑ = 0. */
    constexpr void reset(const ColVec<NX + 1, T>& z0) {
        eso_.reset(z0);
        u_cmd_ = T{0};
    }

private:
    [[no_unique_address]] GainSet k_{};          ///< ω_c, ω_o, b₀ (the tick divides by b₀), K_p, K_d and T_s: members, or no storage on a constexpr design's Law
    Bounds<1, T>                  u_lim_{};      ///< The command's limits
    ESO<NX, T, EsoStore>          eso_{};        ///< Discrete ESO (Miklosovic–Radke–Gao current estimator)
    T                             u_cmd_{};      ///< The command the last step() returned, or the input track() recorded
    bool                          valid_{false}; ///< Built from a design it runs
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<size_t NX, typename T>
ADRCController(const design::ADRCResult<NX, T>&) -> ADRCController<NX, T>;

namespace detail {

/** @brief The runtime an ADRC design runs as: none for anything else. */
template<const auto& Res, typename T, typename R = detail::design_t<Res>>
struct adrc_runtime {};

/** @brief design::adrc: with members, or its gains and ESO on the design's constants. */
template<const auto& Res, typename T, size_t NX, typename U>
struct adrc_runtime<Res, T, design::ADRCResult<NX, U>> {
    using type = ADRCController<NX, T>;                                    ///< With members
    using on_constants = ADRCController<NX, T, ConstantADRCGains<Res, T>>; ///< Its gains and ESO on the design's constants
};

/** @brief An ADRC design. */
template<const auto& Res, typename T>
concept adrc_design = requires { typename adrc_runtime<Res, T>::type; };

} // namespace detail

/** @brief Law<kRes, T> of a design::adrc initialized at start-up: the ADRCController it runs as. A runtime of the same shape assigned replaces it. */
template<const auto& Res, typename T>
    requires detail::adrc_design<Res, T>
struct Law<Res, T> : detail::adrc_runtime<Res, T>::type {
    using Runtime = detail::adrc_runtime<Res, T>::type; ///< The runtime the design runs as

    /** @brief The design's law. */
    constexpr Law() : Runtime{Res} {}

    using Runtime::operator=;
};

/**
 * @brief Law<kRes, T> of a `constexpr` design::adrc: the same law with its gains and its ESO's
 *        maps constants of the type (detail::ConstantADRCGains, detail::ConstantEso)
 *
 * K_p, K_d, b₀, T_s and the ESO's Φ, Γ and L are in flash and fold into the tick; the ESO's
 * estimate and the limits stay members, the limits settable at run time; there is no set_gains.
 * A law retuned at run time is `Law<kRes, T>::Runtime law{kRes.as<T>()};`. A failed design does
 * not compile wherever the Law is declared.
 */
template<const auto& Res, typename T>
    requires(detail::adrc_design<Res, T> && detail::constant_design<Res>)
struct Law<Res, T> : detail::adrc_runtime<Res, T>::on_constants {
    using Runtime = detail::adrc_runtime<Res, T>::type; ///< The runtime with members: the type of a member a redesign assigns
    using Base = detail::adrc_runtime<Res, T>::on_constants;

    static_assert(detail::constant_design_runs(Res.success, Res.status, Res.Ts));
    static_assert(Res.limits.input.valid(), "The design's limits are not a box: an end is not finite, or the lower one is above the upper one.");
    static_assert(detail::ConstantADRCGains<Res, T>::ok, "Law<kRes, T>: the ADRC design's b0 is not finite and positive, or a PD gain is not finite, in the runtime's scalar");

    /** @brief The design's law: its gains and ESO maps the type's, its limits members. */
    constexpr Law() : Base(Res.template as<T>()) {}

    /** @brief There is nothing to replace: the gains are in the type. */
    template<typename Other>
        requires(!damp::is_same_v<damp::remove_cvref_t<Other>, Law>)
    constexpr Law& operator=(const Other& other) {
        static_assert(sizeof(Other) == 0, "This Law is a constexpr ADRC design's: its gains and its ESO's maps are constants of the type, so there is nothing to replace. A member a redesign assigns, or a law retuned with set_gains, is Law<kRes, T>::Runtime law{kRes.as<T>()};");
        static_cast<void>(other);
        return *this;
    }
};
} // namespace damp
