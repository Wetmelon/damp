// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file pid_design.hpp
 * @brief design::pid: a PID from a tuning rule on a model, or from gains you have, deployed at Ts
 *
 * @defgroup pid_design PID Design
 * @brief One design call: design::pid(model, Ts, rule, opts...) → PIDResult → PIDController
 *
 * The model is what you know of the plant: a first-order-plus-dead-time fit (FOPDT), a relay
 * test's ultimate point (Ultimate), or a continuous StateSpace or TransferFunction. The rule turns
 * it into gains, and Ts deploys them. In place of the model and the rule, gains you already have:
 * Parallel{.Kp, .Ki, .Kd} or Standard{.Kp, .Ti, .Td}.
 *
 * | Rule | Models | What it does |
 * | --- | --- | --- |
 * | SIMC{.tau_c, .mode} | FOPDT | Skogestad's IMC rule, one number: the closed-loop time constant |
 * | AMIGO{} | FOPDT; Ultimate with .K | Åström and Hägglund's robust PI (M_s ≈ 1.4) |
 * | LoopShape{.wc, .phase_margin, .mode} | StateSpace, TransferFunction; Ultimate | a crossover and a phase margin, the closed loop certified on a plant; the relay autotuner's phase-margin design at the ultimate point |
 * | PolePlacement{.omega, .zeta, .omega_i} | first- or second-order StateSpace / TransferFunction | a PI on the first, a PD or PID with its filter on the second, in closed form |
 * | Poles{p…} | first-order StateSpace / TransferFunction | the s-plane poles of the law sampled at Ts, exact: two for a PI, three for a PID |
 *
 * The options trail in any order, at most one of each kind (design/options.hpp), and set what the
 * rule does not: InputLimits{.lower, .upper} (the command), DerivativeFilter{.Tf},
 * SetpointWeights{.b, .c} and AntiWindup{.tracking_time, .limits}. A rule that cannot use a model,
 * an option twice or of another design, and a rule among the options do not compile.
 *
 * @code
 * constexpr auto kTemp = design::pid(design::FOPDT{.K = 2.5, .L = 4.0, .tau = 30.0}, 0.5, design::SIMC{.tau_c = 4.0},
 *                                    design::InputLimits{.lower = {0.0}, .upper = {100.0}});
 * static constinit Law<kTemp, float> heater; // each tick: u = heater.step(r, y);
 * @endcode
 *
 * Every design reports through PIDResult::success and status, with zero gains when it fails:
 * invalid_sample_time for a Ts that is not finite and positive or a discrete StateSpace, and
 * invalid_parameter for a model, rule or option number outside its domain; LoopShape and the
 * placements add singular, not_placeable and not_stabilizable.
 */

#include "damp/backend.hpp"
#include "damp/controllers/pid.hpp" // IWYU pragma: export
#include "damp/design/detail/pid_rules.hpp"
#include "damp/design/input_limits.hpp" // IWYU pragma: export
#include "damp/design/options.hpp"
#include "damp/design/poles.hpp" // IWYU pragma: export
#include "damp/design/status.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/toolbox/bounds.hpp"

// Arduino.h defines PI as a macro, which would rewrite PIDMode::PI: set aside for this header.
#pragma push_macro("PI")
#undef PI

namespace damp {
namespace design {

// ============================================================================
// Models
// ============================================================================

/**
 * @struct FOPDT
 * @brief A first-order-plus-dead-time model K·e^{−Ls}/(τs + 1): what a step test fits, and what SIMC and AMIGO tune
 *
 * `FOPDT{.K = 2.5, .L = 4.0, .tau = 30.0}`: the scalar is K's (the others convert to it). A plant
 * of higher order reduces to one by Skogestad's half rule, or by a fit to its step response.
 *
 * @tparam T The scalar (deduced from K)
 */
template<typename T = double>
struct FOPDT {
    static_assert(damp::is_floating_point_v<T>, "FOPDT's numbers are floating-point: FOPDT{.K = 2.0, .L = 0.5, .tau = 3.0}.");
    T                        K{};   ///< Static gain (finite and nonzero; < 0 for a reverse-acting plant)
    damp::type_identity_t<T> L{};   ///< [s] Dead time (≥ 0; AMIGO needs > 0)
    damp::type_identity_t<T> tau{}; ///< [s] Time constant (> 0)
};

/**
 * @struct Ultimate
 * @brief A relay test's ultimate point: the gain K_u and period T_u at which the loop oscillates, and the static gain K
 *
 * What design::RelayAutotuner measures, `Ultimate{.Ku = out.Ku, .Tu = out.Tu}`: the plant's
 * response at ω_u = 2π/T_u is P(jω_u) = −1/K_u. LoopShape fits the law there; AMIGO also needs
 * the static gain (a biased relay's Ks, or a step test's), `.K = out.Ks`, and refuses K = 0, not
 * known.
 *
 * @tparam T The scalar (deduced from Ku)
 */
template<typename T = double>
struct Ultimate {
    static_assert(damp::is_floating_point_v<T>, "Ultimate's numbers are floating-point: Ultimate{.Ku = 4.0, .Tu = 2.0}.");
    T                        Ku{}; ///< Ultimate gain (> 0)
    damp::type_identity_t<T> Tu{}; ///< [s] Ultimate period (> 0)
    damp::type_identity_t<T> K{};  ///< Static gain, which only AMIGO reads (0: not known)
};

/**
 * @struct Parallel
 * @brief Gains you already have, in the parallel form Kp + Kᵢ/s + K_d s: design::pid's model and rule in one
 *
 * `design::pid(Parallel{.Kp = 4.0, .Ki = 8.0}, Ts, opts...)`: the scalar is Kp's. The options add
 * the derivative filter, the setpoint weights, the anti-windup and the limits.
 *
 * @tparam T The scalar (deduced from Kp)
 */
template<typename T = double>
struct Parallel {
    static_assert(damp::is_floating_point_v<T>, "Parallel's gains are floating-point: Parallel{.Kp = 4.0, .Ki = 8.0}.");
    T                        Kp{}; ///< Proportional gain
    damp::type_identity_t<T> Ki{}; ///< Integral gain [1/s]
    damp::type_identity_t<T> Kd{}; ///< Derivative gain [s]
};

/**
 * @struct Standard
 * @brief Gains you already have, in the standard (ISA) form Kp(1 + 1/(Tᵢ s) + T_d s): Kᵢ = Kp/Tᵢ, K_d = Kp·T_d
 *
 * `design::pid(Standard{.Kp = 4.0, .Ti = 0.5, .Td = 0.0025}, Ts, opts...)`: Tᵢ = 0 means no
 * integral, and a Tᵢ or T_d below 0 is refused with DesignStatus::invalid_parameter.
 *
 * @note Compare with MATLAB®'s pidstd(Kp, Ti, Td) (matlab::pidstd), whose filter N is a
 *       DerivativeFilter{.Tf = Td/N} here.
 * @tparam T The scalar (deduced from Kp)
 */
template<typename T = double>
struct Standard {
    static_assert(damp::is_floating_point_v<T>, "Standard's gains are floating-point: Standard{.Kp = 4.0, .Ti = 0.5}.");
    T                        Kp{}; ///< Proportional gain
    damp::type_identity_t<T> Ti{}; ///< [s] Integral time (0: no integral)
    damp::type_identity_t<T> Td{}; ///< [s] Derivative time (0: no derivative)
};

// ============================================================================
// Rules
// ============================================================================

/**
 * @struct SIMC
 * @brief Skogestad's SIMC rule on an FOPDT: Kp = τ/(K(τ_c + L)), Tᵢ = min(τ, 4(τ_c + L))
 *
 * One number, τ_c, the closed-loop time constant: τ_c = L is the rule's own choice, tight
 * control with good robustness, and a larger τ_c is slower and more robust. The mode picks the
 * law, a PI by default; P drops the integral, and the derivative of a PD or PID, K_d = Kp·L/2, is
 * a dead-time heuristic, not SIMC's (design::detail::simc_gains).
 *
 * @tparam T The scalar of τ_c (converted to the model's)
 */
template<typename T = double>
struct SIMC {
    using option_kind = kind::Rule; ///< A rule: the argument after the model and Ts, never an option
    T       tau_c{};                ///< [s] Closed-loop time constant τ_c (τ_c + L > 0)
    PIDMode mode{PIDMode::PI};      ///< P, PI, PD or PID
};

/**
 * @brief Åström and Hägglund's AMIGO PI: the robust PI (M_s ≈ 1.4) on an FOPDT, or on a relay test's ultimate point and static gain
 *
 * On an Ultimate it first recovers the FOPDT with that ultimate point and static gain, so it needs
 * `.K`. A PI only: it takes no DerivativeFilter (design::detail::amigo_pi,
 * design::detail::amigo_from_ultimate).
 */
struct AMIGO {
    using option_kind = kind::Rule; ///< A rule: the argument after the model and Ts, never an option
};

/**
 * @struct LoopShape
 * @brief A crossover ωc and a phase margin: the loop C·P crosses the unit circle at ωc with that margin
 *
 * On a continuous StateSpace or TransferFunction the law is fitted to the plant's response at
 * .wc, the filter of a DerivativeFilter inside the fit, and the closed loop certified (every pole
 * in the open left half plane, else not_stabilizable). With SetpointWeights other than 1 the
 * crossover is raised so the reference response keeps the 1-DOF design's bandwidth. A PID has
 * Tᵢ = 4T_d. On a relay test's Ultimate the crossover is the ultimate frequency, so .wc stays 0
 * (another is refused with invalid_parameter): the design moves the measured point onto the unit
 * circle, the weights do not move it, and with one point of the plant the loop is not certified.
 * There a PD or a PID; a PI cannot add the lead the point needs, and a P is K_u itself
 * (design::detail::loop_shape, design::detail::loop_shape_ultimate).
 *
 * @note Compare with MATLAB®'s pidtune(sys, type, wc, pidtuneOptions('PhaseMargin', pm)).
 * @tparam T The scalar of ωc (converted to the model's)
 */
template<typename T = double>
struct LoopShape {
    using option_kind = kind::Rule;               ///< A rule: the argument after the model and Ts, never an option
    T                        wc{};                ///< [rad/s] Loop crossover ωc (> 0; 0 on an Ultimate: its ultimate frequency)
    damp::type_identity_t<T> phase_margin{T{60}}; ///< [deg] Phase margin at the crossover, in (0, 180)
    PIDMode                  mode{PIDMode::PID};  ///< P, PI, PD or PID
};

/**
 * @struct PolePlacement
 * @brief The closed-loop poles of a first- or second-order plant in closed form: the pair s² + 2ζω s + ω², and the integrator's at −ωᵢ
 *
 * On a first-order plant b₀/(a₁ s + a₀) a PI places the pair (ωᵢ = 0). On a second-order plant
 * b₀/(a₂ s² + a₁ s + a₀) a PD places the pair and its derivative filter's pole at −10ω, and a PID
 * with ωᵢ > 0 the integrator's pole too; the derivative is on the measurement (c = 0), unless a
 * SetpointWeights sets both weights. The rule sets the filter, so it takes no DerivativeFilter.
 * Keep ω about a tenth of the sample rate or below for the continuous placement to hold at Ts
 * (design::detail::pole_pair_pi, design::detail::pole_placement_pid).
 *
 * @tparam T The scalar of ω (converted to the model's)
 */
template<typename T = double>
struct PolePlacement {
    using option_kind = kind::Rule;         ///< A rule: the argument after the model and Ts, never an option
    T                        omega{};       ///< [rad/s] Frequency ω of the placed pair (> 0)
    damp::type_identity_t<T> zeta{T{1}};    ///< Damping ratio ζ of the pair (> 0); 1 = critically damped
    damp::type_identity_t<T> omega_i{T{0}}; ///< [rad/s] Integrator pole at −ωᵢ, second-order plant only (0: a PD)
};

// ============================================================================
// Options
// ============================================================================

/**
 * @struct DerivativeFilter
 * @brief The derivative's filter K_d s/(1 + T_f s), among design::pid's options
 *
 * With gains given (Parallel, Standard) and with SIMC it is the law's filter; LoopShape fits the
 * law with it, so the delivered loop meets the crossover and the phase margin. AMIGO (a PI),
 * PolePlacement (which places its own) and Poles (exact for an unfiltered derivative) refuse it at
 * compile time. A T_f below 0 or not finite is refused with DesignStatus::invalid_parameter.
 *
 * @tparam T The scalar (converted to the design's)
 */
template<typename T = double>
struct DerivativeFilter {
    using option_kind = kind::DerivativeFilter; ///< Its kind of design option
    T Tf{};                                     ///< [s] Filter time constant T_f (≥ 0; 0: unfiltered)
};

/**
 * @struct SetpointWeights
 * @brief The setpoint weights of a 2-DOF PID, u = Kp(b·r − y) + I + K_d d/dt(c·r − y), among design::pid's options
 *
 * The weights shape the response to r and leave the loop as it is: b = 0 is an I-P (no
 * proportional kick), c = 0 a PI-D (no derivative kick). Given, both are the law's, so PolePlacement's
 * c = 0 becomes this c. LoopShape raises the crossover on a plant so that the reference response
 * keeps the 1-DOF design's bandwidth. A weight that is not finite is refused with
 * DesignStatus::invalid_parameter.
 *
 * @see "Advanced PID Control" (Åström & Hägglund, 2006), §4.4 (setpoint weighting)
 * @tparam T The scalar (converted to the design's)
 */
template<typename T = double>
struct SetpointWeights {
    using option_kind = kind::SetpointWeights; ///< Its kind of design option
    T                        b{T{1}};          ///< Proportional setpoint weight (1: on the error, 0: I-P)
    damp::type_identity_t<T> c{T{1}};          ///< Derivative setpoint weight (1: on the error, 0: on the measurement)
};

/**
 * @struct AntiWindup
 * @brief How the integral meets a clamp: the tracking time, and the integral term's own limits, among design::pid's options
 *
 * Without it, as with a tracking time of 0, the integral holds where integrating would push a
 * clamped command further (conditional integration): it keeps what it carries, which is why that
 * is the default, as on the LQ laws (IntegralAction). A tracking time T_t > 0 unwinds it by
 * back-calculation instead, min(T_s/T_t, 1) of the clamp's excess each tick. `.limits` clamps the
 * integral term I itself.
 *
 * Not IntegralAction, which an LQ design takes: its `.limits` bound ξ ≈ ∫e dt, which a gain row
 * turns into u, while a PID holds its integral as the I term, already in u's units, so a clamp on
 * it is in u's units too (PIDLimits::integral). A tracking time below 0 or not finite, or limits
 * that are not a box, are refused with DesignStatus::invalid_parameter.
 *
 * The same T_t governs the law's own clamp and a clamp after it, reported with back_calculate(u_sat,
 * u_ff), and holding suits both: a clamp is a clamp wherever it sits. An input from another source
 * (an operator, another law) is not a clamp; that is what track(u, r, y) is for.
 *
 * @note Compare with MATLAB®'s PID Controller block, whose back-calculation coefficient Kb (its own
 *       saturation) and tracking coefficient Kt (the TR input, which back_calculate takes here) are
 *       both ≈ 1/T_t. Damp has one T_t for the two, and its 0 is the block's clamping method.
 *
 * @see "Advanced PID Control" (Åström & Hägglund, 2006), §3.5: back-calculation and conditional
 *      integration
 * @tparam T The scalar (converted to the design's)
 */
template<typename T = double>
struct AntiWindup {
    using option_kind = kind::AntiWindup;                ///< Its kind of design option
    T                                   tracking_time{}; ///< [s] Back-calculation tracking time T_t (0: conditional integration)
    damp::type_identity_t<Bounds<1, T>> limits{};        ///< The integral term I's limits, in u's units (default unbounded)
};

namespace detail {

/** @brief The options design::pid takes. A rule or Poles among them is the rule given twice. */
using pid_options = option_set<kind::InputLimits, kind::DerivativeFilter, kind::SetpointWeights, kind::AntiWindup, kind::Rule, kind::Weights>;

/** @brief A PID design's options as it reads them, in its scalar @p T: left out, their defaults. */
template<typename T>
struct PIDOptions {
    PIDLimits<T> limits{};        ///< InputLimits and AntiWindup::limits
    T            Tf{T{0}};        ///< DerivativeFilter::Tf
    bool         filtered{false}; ///< A DerivativeFilter is given
    T            b{T{1}};         ///< SetpointWeights::b
    T            c{T{1}};         ///< SetpointWeights::c
    bool         weighted{false}; ///< A SetpointWeights is given
    T            tracking_time{}; ///< AntiWindup::tracking_time

    /** @brief Every number in its domain: Tf and the tracking time finite and ≥ 0, the weights finite, the limits boxes. */
    [[nodiscard]] constexpr bool valid() const {
        return damp::finite_non_negative(Tf) && damp::isfinite(b) && damp::isfinite(c) && damp::finite_non_negative(tracking_time)
            && limits.input.valid() && limits.integral.valid();
    }
};

/** @brief The options among @p opts, read by kind into the design's scalar @p T. */
template<typename T, typename... Opts>
[[nodiscard]] constexpr PIDOptions<T> pid_options_of(const Opts&... opts) {
    PIDOptions<T> o{};
    if constexpr (option_count<kind::InputLimits, Opts...> > 0) {
        if constexpr (option_t<kind::InputLimits, void, Opts...>::inputs == 1) { // more is refused at compile time
            o.limits.input = input_limits_of<1, T>(opts...);
        }
    }
    if constexpr (option_count<kind::DerivativeFilter, Opts...> > 0) {
        o.Tf = static_cast<T>(the<kind::DerivativeFilter>(opts...).Tf);
        o.filtered = true;
    }
    if constexpr (option_count<kind::SetpointWeights, Opts...> > 0) {
        const auto w = the<kind::SetpointWeights>(opts...);
        o.b = static_cast<T>(w.b);
        o.c = static_cast<T>(w.c);
        o.weighted = true;
    }
    if constexpr (option_count<kind::AntiWindup, Opts...> > 0) {
        const auto aw = the<kind::AntiWindup>(opts...);
        o.tracking_time = static_cast<T>(aw.tracking_time);
        o.limits.integral = aw.limits.template as<T>();
    }
    (static_cast<void>(opts), ...);
    return o;
}

/**
 * @brief The compile-time check of design::pid's options: their kinds, one of each, one input, and
 *        no rule among them (@p kGains: the call took gains, Parallel or Standard, which take none)
 */
template<bool kGains, typename... Opts>
consteval void check_pid_options() {
    static_assert(
        option_count<kind::IntegralAction, Opts...> == 0,
        "A PID's integral is its I term, in u's units: AntiWindup{.tracking_time = ..., .limits = {...}} sets how it meets a clamp. IntegralAction is an LQ design's."
    );
    pid_options::check<Opts...>();
    constexpr bool kRule = ((OptionOf<Opts, kind::Rule> || OptionOf<Opts, kind::Weights>) || ...);
    if constexpr (kGains) {
        static_assert(
            !kRule,
            "Gains take no rule: Parallel{...} and Standard{...} are the gains, pid(gains, Ts, opts...); a rule tunes a model, pid(model, Ts, rule, opts...)."
        );
    } else {
        static_assert(
            !kRule,
            "The rule given twice: SIMC{...}, AMIGO{}, LoopShape{...}, PolePlacement{...} or Poles{...} is the one argument after the model and Ts, pid(model, Ts, rule, opts...)."
        );
    }
    if constexpr (option_count<kind::InputLimits, Opts...> > 0) {
        static_assert(option_t<kind::InputLimits, void, Opts...>::inputs == 1, "A PID has one input: InputLimits{.lower = {lo}, .upper = {hi}}, one entry at each end.");
    }
}

/** @brief Whether @p R is SIMC. */
template<typename R>
inline constexpr bool is_simc = false;

/** @brief SIMC is. */
template<typename U>
inline constexpr bool is_simc<SIMC<U>> = true;

/** @brief Whether @p R is LoopShape. */
template<typename R>
inline constexpr bool is_loop_shape = false;

/** @brief LoopShape is. */
template<typename U>
inline constexpr bool is_loop_shape<LoopShape<U>> = true;

/** @brief Whether @p R is PolePlacement. */
template<typename R>
inline constexpr bool is_pole_placement = false;

/** @brief PolePlacement is. */
template<typename U>
inline constexpr bool is_pole_placement<PolePlacement<U>> = true;

/** @brief Whether @p R is Poles. */
template<typename R>
inline constexpr bool is_poles = false;

/** @brief Poles is. */
template<size_t N, typename U>
inline constexpr bool is_poles<Poles<N, U>> = true;

/** @brief What a model is, as a rule reads it: none of design::pid's. */
template<typename M>
struct pid_model {
    static constexpr bool   fopdt = false;    ///< An FOPDT
    static constexpr bool   ultimate = false; ///< An Ultimate
    static constexpr bool   plant = false;    ///< A continuous StateSpace or TransferFunction
    static constexpr size_t order = 0;        ///< A plant's order (its states, its denominator's degree)
    static constexpr bool   all_pole = false; ///< A plant whose form has no zero (a StateSpace's is checked at run time)
};

/** @brief An FOPDT. */
template<typename U>
struct pid_model<FOPDT<U>> : pid_model<void> {
    static constexpr bool fopdt = true; ///< An FOPDT
};

/** @brief An Ultimate. */
template<typename U>
struct pid_model<Ultimate<U>> : pid_model<void> {
    static constexpr bool ultimate = true; ///< An Ultimate
};

/** @brief A SISO StateSpace of NX states. */
template<size_t NX, typename U, size_t NW, size_t NV>
struct pid_model<StateSpace<NX, 1, 1, U, NW, NV>> : pid_model<void> {
    static constexpr bool   plant = true;    ///< A plant
    static constexpr size_t order = NX;      ///< Its states
    static constexpr bool   all_pole = true; ///< Checked at run time (all_pole_tf)
};

/** @brief A TransferFunction of denominator degree ND − 1, all-pole when its numerator is a constant. */
template<size_t NN, size_t ND, typename U>
struct pid_model<TransferFunction<NN, ND, U>> : pid_model<void> {
    static constexpr bool   plant = true;       ///< A plant
    static constexpr size_t order = ND - 1;     ///< Its denominator's degree
    static constexpr bool   all_pole = NN == 1; ///< num = {b₀}
};

/** @brief Whether the rule @p Rule takes the model @p Model; refused at compile time, with the models it takes, when it does not. */
template<typename Rule, typename Model>
consteval bool pid_rule_fits() {
    using M = pid_model<Model>;
    if constexpr (is_simc<Rule>) {
        static_assert(M::fopdt, "SIMC tunes an FOPDT model: pid(FOPDT{.K = ..., .L = ..., .tau = ...}, Ts, SIMC{.tau_c = ...}).");
        return M::fopdt;
    } else if constexpr (damp::is_same_v<Rule, AMIGO>) {
        static_assert(M::fopdt || M::ultimate, "AMIGO tunes an FOPDT model, FOPDT{.K, .L, .tau}, or a relay test's Ultimate{.Ku, .Tu, .K} with its static gain.");
        return M::fopdt || M::ultimate;
    } else if constexpr (is_loop_shape<Rule>) {
        constexpr bool ok = (M::plant && M::order >= 1) || M::ultimate;
        static_assert(ok, "LoopShape fits the frequency response of a continuous StateSpace or TransferFunction with a pole, or a relay test's Ultimate{.Ku, .Tu}.");
        return ok;
    } else if constexpr (is_pole_placement<Rule>) {
        constexpr bool ok = M::plant && M::all_pole && (M::order == 1 || M::order == 2);
        static_assert(
            ok,
            "PolePlacement places the poles of a first- or second-order all-pole plant: a StateSpace of one or two states, or TransferFunction{.num = {b0}, .den = {a1, a0}} or {.num = {b0}, .den = {a2, a1, a0}}."
        );
        return ok;
    } else if constexpr (is_poles<Rule>) {
        constexpr bool ok = M::plant && M::all_pole && (M::order == 1);
        constexpr bool count = (listed_by<kind::Weights, Rule> == 2) || (listed_by<kind::Weights, Rule> == 3);
        static_assert(ok, "Poles places a first-order plant sampled at Ts: a StateSpace of one state, or TransferFunction{.num = {b0}, .den = {a1, a0}}.");
        static_assert(!ok || count, "Poles on a first-order plant: two for a PI, three for a PID.");
        return ok && count;
    } else {
        static_assert(
            !damp::is_same_v<Rule, Rule>,
            "The rule after the model and Ts: SIMC{...}, AMIGO{}, LoopShape{...}, PolePlacement{...} or Poles{...}, pid(model, Ts, rule, opts...)."
        );
        return false;
    }
}

/** @brief Whether the options @p Opts fit the rule @p Rule: refused at compile time with why when they do not. */
template<typename Rule, typename... Opts>
consteval bool pid_options_fit() {
    constexpr bool kFilter = option_count<kind::DerivativeFilter, Opts...> > 0;
    if constexpr (damp::is_same_v<Rule, AMIGO>) {
        static_assert(!kFilter, "AMIGO is a PI rule: there is no derivative for a DerivativeFilter to filter.");
    } else if constexpr (is_pole_placement<Rule>) {
        static_assert(!kFilter, "PolePlacement places the derivative filter's pole itself, at 10 omega, and a first-order plant's PI has no derivative: no DerivativeFilter.");
    } else if constexpr (is_poles<Rule>) {
        static_assert(!kFilter, "Poles places the sampled law with an unfiltered derivative, exact at Ts: no DerivativeFilter.");
    }
    return !kFilter || !(damp::is_same_v<Rule, AMIGO> || is_pole_placement<Rule> || is_poles<Rule>);
}

/** @brief LoopShape's spec with the filter and the weights among the options. */
template<typename T, typename U>
[[nodiscard]] constexpr LoopSpec<T> loop_spec_of(const LoopShape<U>& rule, const PIDOptions<T>& o) {
    return LoopSpec<T>{
        .wc = static_cast<T>(rule.wc),
        .phase_margin_deg = static_cast<T>(rule.phase_margin),
        .mode = rule.mode,
        .Tf = o.Tf,
        .b = o.b,
        .c = o.c,
    };
}

/** @brief SIMC on an FOPDT. */
template<typename T, typename U>
[[nodiscard]] constexpr PIDTuning<T> tune(const FOPDT<T>& m, const SIMC<U>& rule, T /*Ts*/, const PIDOptions<T>& /*o*/) {
    return simc_gains(m.K, m.L, m.tau, static_cast<T>(rule.tau_c), rule.mode);
}

/** @brief AMIGO on an FOPDT. */
template<typename T>
[[nodiscard]] constexpr PIDTuning<T> tune(const FOPDT<T>& m, const AMIGO& /*rule*/, T /*Ts*/, const PIDOptions<T>& /*o*/) {
    return amigo_pi(m.K, m.L, m.tau);
}

/** @brief AMIGO on an ultimate point and its static gain. */
template<typename T>
[[nodiscard]] constexpr PIDTuning<T> tune(const Ultimate<T>& m, const AMIGO& /*rule*/, T /*Ts*/, const PIDOptions<T>& /*o*/) {
    return amigo_from_ultimate(m.K, m.Ku, m.Tu);
}

/** @brief LoopShape at an ultimate point. */
template<typename T, typename U>
[[nodiscard]] constexpr PIDTuning<T> tune(const Ultimate<T>& m, const LoopShape<U>& rule, T /*Ts*/, const PIDOptions<T>& o) {
    if (!damp::isfinite(m.K)) {
        return pid_invalid<T>();
    }
    return loop_shape_ultimate(m.Ku, m.Tu, loop_spec_of(rule, o));
}

/** @brief LoopShape on a StateSpace. */
template<size_t NX, typename T, size_t NW, size_t NV, typename U>
[[nodiscard]] constexpr PIDTuning<T> tune(const StateSpace<NX, 1, 1, T, NW, NV>& sys, const LoopShape<U>& rule, T /*Ts*/, const PIDOptions<T>& o) {
    return loop_shape(sys, loop_spec_of(rule, o));
}

/** @brief LoopShape on a TransferFunction: its controllable realization. */
template<size_t NN, size_t ND, typename T, typename U>
[[nodiscard]] constexpr PIDTuning<T> tune(const TransferFunction<NN, ND, T>& tf, const LoopShape<U>& rule, T /*Ts*/, const PIDOptions<T>& o) {
    const auto sys = tf.to_state_space();
    if (!sys) {
        return pid_invalid<T>(); // den[0] = 0
    }
    return loop_shape(*sys, loop_spec_of(rule, o));
}

/** @brief PolePlacement on a first-order transfer function: the PI. */
template<typename T, typename U>
[[nodiscard]] constexpr PIDTuning<T> tune(const TransferFunction<1, 2, T>& tf, const PolePlacement<U>& rule, T /*Ts*/, const PIDOptions<T>& /*o*/) {
    return pole_pair_pi(tf, static_cast<T>(rule.omega), static_cast<T>(rule.zeta), static_cast<T>(rule.omega_i));
}

/** @brief PolePlacement on a second-order transfer function: the PD or PID. */
template<typename T, typename U>
[[nodiscard]] constexpr PIDTuning<T> tune(const TransferFunction<1, 3, T>& tf, const PolePlacement<U>& rule, T /*Ts*/, const PIDOptions<T>& /*o*/) {
    return pole_placement_pid(tf, static_cast<T>(rule.omega), static_cast<T>(rule.zeta), static_cast<T>(rule.omega_i));
}

/** @brief PolePlacement on a StateSpace: its all-pole transfer function. */
template<size_t NX, typename T, size_t NW, size_t NV, typename U>
[[nodiscard]] constexpr PIDTuning<T> tune(const StateSpace<NX, 1, 1, T, NW, NV>& sys, const PolePlacement<U>& rule, T Ts, const PIDOptions<T>& o) {
    const auto all_pole = all_pole_tf(sys);
    if (all_pole.status != DesignStatus::ok) {
        return refused<PIDTuning<T>>(all_pole.status);
    }
    return tune(all_pole.tf, rule, Ts, o);
}

/** @brief Poles on a first-order transfer function, sampled at Ts. */
template<typename T, size_t N, typename U>
[[nodiscard]] constexpr PIDTuning<T> tune(const TransferFunction<1, 2, T>& tf, const Poles<N, U>& rule, T Ts, const PIDOptions<T>& /*o*/) {
    return sampled_placement<N>(tf, pole_list<T>(rule.poles), Ts);
}

/** @brief Poles on a first-order StateSpace: its all-pole transfer function. */
template<typename T, size_t NW, size_t NV, size_t N, typename U>
[[nodiscard]] constexpr PIDTuning<T> tune(const StateSpace<1, 1, 1, T, NW, NV>& sys, const Poles<N, U>& rule, T Ts, const PIDOptions<T>& o) {
    const auto all_pole = all_pole_tf(sys);
    if (all_pole.status != DesignStatus::ok) {
        return refused<PIDTuning<T>>(all_pole.status);
    }
    return tune(all_pole.tf, rule, Ts, o);
}

/**
 * @brief A rule's gains deployed at @p Ts with the options @p o: the filter and weights they give,
 *        their tracking time and limits; the rule's status when it failed
 *
 * Refused with invalid_parameter for a gain that is not finite, or one whose per-tick coefficient
 * (Kᵢ·Ts, K_d/(Ts + T_f), 1/T_t) overflows, which the runtime would refuse.
 */
template<typename T>
[[nodiscard]] constexpr PIDResult<T> pid_deployed(const PIDTuning<T>& tuning, T Ts, const PIDOptions<T>& o) {
    if (!tuning.success) {
        return refused<PIDResult<T>>(tuning.status);
    }
    PIDGains<T> g = tuning.gains;
    if (o.filtered) {
        g.Tf = o.Tf;
    }
    if (o.weighted) {
        g.b = o.b;
        g.c = o.c;
    }
    g.tracking_time = o.tracking_time;
    if (!pid_gain_set<PIDMode::PID>(g, Ts)) { // the runtime's own check: every gain in its domain, every coefficient finite
        return refused<PIDResult<T>>(DesignStatus::invalid_parameter);
    }
    PIDResult<T> out{.gains = g, .limits = o.limits, .Ts = Ts};
    set_status(out, DesignStatus::ok);
    return out;
}

/** @brief design::pid once its call compiles: Ts, the options, then the gains @p tune_at gives (from the options) deployed at Ts. */
template<typename T, typename Tune, typename... Opts>
[[nodiscard]] constexpr PIDResult<T> pid_with(T Ts, const Tune& tune_at, const Opts&... opts) {
    const PIDOptions<T> o = pid_options_of<T>(opts...);
    if (!damp::finite_positive(Ts)) {
        return refused<PIDResult<T>>(DesignStatus::invalid_sample_time);
    }
    if (!o.valid()) {
        return refused<PIDResult<T>>(DesignStatus::invalid_parameter);
    }
    return pid_deployed(tune_at(o), Ts, o);
}

/** @brief design::pid of a rule on a model: the rule and the options checked at compile time, then designed. */
template<typename T, typename Model, typename Rule, typename... Opts>
[[nodiscard]] constexpr PIDResult<T> pid_of(const Model& model, T Ts, const Rule& rule, const Opts&... opts) {
    check_pid_options<false, Opts...>();
    constexpr bool kFits = pid_rule_fits<Rule, Model>() && pid_options_fit<Rule, Opts...>();
    if constexpr (!kFits) {
        static_cast<void>(model);
        static_cast<void>(rule);
        static_cast<void>(Ts);
        (static_cast<void>(opts), ...);
        return PIDResult<T>{}; // the message says which: nothing further
    } else {
        return pid_with(Ts, [&](const PIDOptions<T>& o) { return tune(model, rule, Ts, o); }, opts...);
    }
}

} // namespace detail

/**
 * @brief A PID from a rule on an FOPDT model, deployed at @p Ts: SIMC or AMIGO
 *
 * @code
 * constexpr auto kTemp = design::pid(design::FOPDT{.K = 2.5, .L = 4.0, .tau = 30.0}, 0.5, design::SIMC{.tau_c = 4.0},
 *                                    design::InputLimits{.lower = {0.0}, .upper = {100.0}});
 * constexpr auto kAmigo = design::pid(design::FOPDT{.K = 2.5, .L = 4.0, .tau = 30.0}, 0.5, design::AMIGO{});
 * @endcode
 *
 * Refused with invalid_sample_time for a Ts that is not finite and positive, and with
 * invalid_parameter for K = 0, L < 0 (AMIGO: L ≤ 0), τ ≤ 0, a rule's number outside its domain, or
 * an option's.
 *
 * @param model The plant K·e^{−Ls}/(τs + 1)
 * @param Ts    [s] The law's sample period
 * @param rule  SIMC{.tau_c, .mode} or AMIGO{}
 * @param opts  InputLimits, DerivativeFilter (SIMC), SetpointWeights, AntiWindup
 * @return The design: the gains, the limits and Ts, for PIDController or Law
 */
template<typename T, typename Rule, typename... Opts>
[[nodiscard]] constexpr PIDResult<T> pid(const FOPDT<T>& model, damp::type_identity_t<T> Ts, const Rule& rule, const Opts&... opts) {
    return detail::pid_of(model, Ts, rule, opts...);
}

/**
 * @brief A PID from a relay test's ultimate point, deployed at @p Ts: LoopShape (Åström and Hägglund's phase-margin design) or AMIGO
 *
 * @code
 * const auto out = tuner.step(y); // design::RelayAutotuner, Done
 * const auto res = design::pid(design::Ultimate{.Ku = out.Ku, .Tu = out.Tu}, Ts, design::LoopShape{.phase_margin = 45.0});
 * @endcode
 *
 * LoopShape moves the measured point P(jω_u) = −1/K_u onto the unit circle at −180° + PM, the
 * crossover at the ultimate frequency (.wc stays 0), with a PD or a PID. The loop is not certified:
 * one point is all of the plant the design has. AMIGO recovers the FOPDT with that point and the
 * static gain `.K` (> 0) and applies its PI.
 *
 * Refused with invalid_sample_time for a Ts that is not finite and positive; invalid_parameter for
 * K_u or T_u not finite and positive, a LoopShape .wc other than 0, AMIGO without K or with
 * K·K_u ≤ 1, or an option out of its domain; not_placeable for a PI under LoopShape, and
 * not_stabilizable for a P (K_u itself).
 *
 * @param model The ultimate point (and the static gain for AMIGO)
 * @param Ts    [s] The law's sample period
 * @param rule  LoopShape{.phase_margin, .mode} or AMIGO{}
 * @param opts  InputLimits, DerivativeFilter (LoopShape), SetpointWeights, AntiWindup
 * @return The design: the gains, the limits and Ts, for PIDController or Law
 */
template<typename T, typename Rule, typename... Opts>
[[nodiscard]] constexpr PIDResult<T> pid(const Ultimate<T>& model, damp::type_identity_t<T> Ts, const Rule& rule, const Opts&... opts) {
    return detail::pid_of(model, Ts, rule, opts...);
}

/**
 * @brief A PID from a rule on a continuous SISO state-space plant, deployed at @p Ts: LoopShape, PolePlacement or Poles
 *
 * @code
 * constexpr auto kLoop = design::pid(plant, 1e-3, design::LoopShape{.wc = 200.0});                 // certified closed loop
 * constexpr auto kPair = design::pid(plant, 1e-3, design::PolePlacement{.omega = 200.0});          // one or two states
 * constexpr auto kZ = design::pid(plant, 1e-3, design::Poles{-200.0, -200.0});                     // one state, exact at Ts
 * @endcode
 *
 * LoopShape fits the plant's frequency response and certifies the closed loop. PolePlacement and
 * Poles take a plant of the order they place and with no finite zero and D = 0, refused with
 * invalid_parameter otherwise; noise channels (G, H) do not enter. A discrete plant is refused
 * with invalid_sample_time, and so is a Ts that is not finite and positive.
 *
 * @note Compare with MATLAB®'s pidtune(sys, type, wc) for LoopShape (matlab::pidtune).
 * @param plant The continuous plant (one input, one output)
 * @param Ts    [s] The law's sample period
 * @param rule  LoopShape{.wc, .phase_margin, .mode}, PolePlacement{.omega, .zeta, .omega_i} or Poles{p…}
 * @param opts  InputLimits, DerivativeFilter (LoopShape), SetpointWeights, AntiWindup
 * @return The design: the gains, the limits and Ts, for PIDController or Law
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Rule, typename... Opts>
[[nodiscard]] constexpr PIDResult<T>
pid(const StateSpace<NX, NU, NY, T, NW, NV>& plant, damp::type_identity_t<T> Ts, const Rule& rule, const Opts&... opts) {
    static_assert(NU == 1 && NY == 1, "A PID has one input and one output: a StateSpace<NX, 1, 1>.");
    if constexpr (NU != 1 || NY != 1) {
        static_cast<void>(plant);
        static_cast<void>(Ts);
        static_cast<void>(rule);
        (static_cast<void>(opts), ...);
        return PIDResult<T>{};
    } else {
        return detail::pid_of(plant, Ts, rule, opts...);
    }
}

/**
 * @brief A PID from a rule on a transfer function, deployed at @p Ts: LoopShape, PolePlacement or Poles
 *
 * As the state-space form: LoopShape on any proper transfer function with a pole (its
 * controllable realization), PolePlacement on `{.num = {b₀}, .den = {a₁, a₀}}` or
 * `{.num = {b₀}, .den = {a₂, a₁, a₀}}`, Poles on the first. A transfer function has no sample time:
 * it is read as continuous.
 *
 * @code
 * constexpr auto kInertia = design::pid(TransferFunction{.num = {1.0}, .den = {0.05, 0.0, 0.0}}, 1e-3,
 *                                       design::PolePlacement{.omega = 8.0, .omega_i = 1.0});
 * @endcode
 *
 * @param plant The continuous plant num(s)/den(s)
 * @param Ts    [s] The law's sample period
 * @param rule  LoopShape{.wc, .phase_margin, .mode}, PolePlacement{.omega, .zeta, .omega_i} or Poles{p…}
 * @param opts  InputLimits, DerivativeFilter (LoopShape), SetpointWeights, AntiWindup
 * @return The design: the gains, the limits and Ts, for PIDController or Law
 */
template<size_t NN, size_t ND, typename T, typename Rule, typename... Opts>
[[nodiscard]] constexpr PIDResult<T> pid(const TransferFunction<NN, ND, T>& plant, damp::type_identity_t<T> Ts, const Rule& rule, const Opts&... opts) {
    return detail::pid_of(plant, Ts, rule, opts...);
}

/**
 * @brief A PID from gains you already have, in the parallel form Kp + Kᵢ/s + K_d s, deployed at @p Ts
 *
 * @code
 * constexpr auto kSpeed = design::pid(design::Parallel{.Kp = 4.0, .Ki = 8.0, .Kd = 0.02}, 1e-3,
 *                                     design::DerivativeFilter{.Tf = 1e-3}, design::InputLimits{.lower = {-12.0}, .upper = {12.0}});
 * static constinit Law<kSpeed, float> speed;
 * @endcode
 *
 * Refused with invalid_sample_time for a Ts that is not finite and positive, and with
 * invalid_parameter for a gain that is not finite or an option out of its domain.
 *
 * @note Compare with MATLAB®'s c2d(pid(Kp, Ki, Kd, Tf), Ts) (matlab::pid, matlab::c2d).
 * @param gains The gains
 * @param Ts    [s] The law's sample period
 * @param opts  InputLimits, DerivativeFilter, SetpointWeights, AntiWindup
 * @return The design: the gains, the limits and Ts, for PIDController or Law
 */
template<typename T, typename... Opts>
[[nodiscard]] constexpr PIDResult<T> pid(const Parallel<T>& gains, damp::type_identity_t<T> Ts, const Opts&... opts) {
    detail::check_pid_options<true, Opts...>();
    return detail::pid_with(
        Ts,
        [&](const detail::PIDOptions<T>& /*o*/) { return detail::tuned(PIDGains<T>{.Kp = gains.Kp, .Ki = gains.Ki, .Kd = gains.Kd}); },
        opts...
    );
}

/**
 * @brief A PID from gains you already have, in the standard (ISA) form Kp(1 + 1/(Tᵢ s) + T_d s), deployed at @p Ts
 *
 * Kᵢ = Kp/Tᵢ (none at Tᵢ = 0) and K_d = Kp·T_d, then as the parallel form.
 *
 * @code
 * constexpr auto kSpeed = design::pid(design::Standard{.Kp = 4.0, .Ti = 0.5, .Td = 0.0025}, 1e-3);
 * @endcode
 *
 * Refused with invalid_sample_time for a Ts that is not finite and positive, and with
 * invalid_parameter for Tᵢ or T_d below 0, a number that is not finite, or an option out of its
 * domain.
 *
 * @note Compare with MATLAB®'s c2d(pidstd(Kp, Ti, Td), Ts) (matlab::pidstd, matlab::c2d).
 * @param gains The gains
 * @param Ts    [s] The law's sample period
 * @param opts  InputLimits, DerivativeFilter, SetpointWeights, AntiWindup
 * @return The design: the gains, the limits and Ts, for PIDController or Law
 */
template<typename T, typename... Opts>
[[nodiscard]] constexpr PIDResult<T> pid(const Standard<T>& gains, damp::type_identity_t<T> Ts, const Opts&... opts) {
    detail::check_pid_options<true, Opts...>();
    return detail::pid_with(
        Ts,
        [&](const detail::PIDOptions<T>& /*o*/) {
        if (!damp::finite_non_negative(gains.Ti) || !damp::finite_non_negative(gains.Td) || !damp::isfinite(gains.Kp)) {
            return detail::pid_invalid<T>();
        }
        const T Ki = (gains.Ti > T{0}) ? (gains.Kp / gains.Ti) : T{0};
        return detail::tuned(PIDGains<T>{.Kp = gains.Kp, .Ki = Ki, .Kd = gains.Kp * gains.Td});
    },
        opts...
    );
}

/**
 * @brief Map percent overshoot target to equivalent damping ratio
 *
 * Uses the standard second-order relation:
 *
 *   OS = exp(-zeta*pi/sqrt(1-zeta^2))
 *
 * with OS as a unit fraction (e.g. 10% -> 0.10).
 *
 * @param overshoot_percent Desired percent overshoot in [0, 100)
 * @return Estimated damping ratio zeta
 */
template<typename T = double>
[[nodiscard]] constexpr T damping_ratio_from_overshoot_percent(T overshoot_percent) {
    constexpr T pi = damp::numbers::pi_v<T>;

    if (overshoot_percent <= T{0}) {
        return T{1};
    }

    T os = overshoot_percent / T{100};
    if (os >= T{1}) {
        os = static_cast<T>(0.999);
    }

    const T log_os = damp::log(os);
    return -log_os / damp::sqrt((pi * pi) + (log_os * log_os));
}

/**
 * @brief Approximate phase margin from damping ratio
 *
 * Exact for the standard second-order loop L(s) = ωₙ²/(s(s + 2ζωₙ)):
 *
 *   PM = atan(2*zeta / sqrt(sqrt(1 + 4*zeta^4) - 2*zeta^2))
 *
 * valid for every ζ > 0 (PM → 90° as ζ → ∞). The inner difference is evaluated
 * as 1/(√(1+4ζ⁴) + 2ζ²) to avoid cancellation at large ζ.
 *
 * @param zeta Damping ratio
 * @return Phase margin in degrees; 30° fallback for ζ ≤ 0 or not finite
 */
template<typename T = double>
[[nodiscard]] constexpr T phase_margin_from_damping_ratio(T zeta) {
    constexpr T pi = damp::numbers::pi_v<T>;

    if (!damp::finite_positive(zeta)) {
        return T{30};
    }

    const T z2 = zeta * zeta;
    // √(1+4ζ⁴) − 2ζ² rewritten without the subtraction.
    const T inner = T{1} / (damp::sqrt(T{1} + (T{4} * z2 * z2)) + (T{2} * z2));
    const T pm_rad = damp::atan2(T{2} * zeta, damp::sqrt(inner));
    return pm_rad * T{180} / pi;
}

} // namespace design
} // namespace damp

#pragma pop_macro("PI")
