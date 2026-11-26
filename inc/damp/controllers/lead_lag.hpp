// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once


/**
 * @file lead_lag.hpp
 * @brief Lead-lag compensators: design::lead, lag, lead_lag and lead_lag_direct → design::LeadLagResult
 *        → LeadLagController, one first-order section or a lead and a lag in series
 *
 * @code
 * using namespace damp;
 * constexpr auto kLead = design::lead({.phase_max = 0.5, .wc = 100.0}, 1e-3); // φ_max [rad], ωc [rad/s], Ts [s]
 * constexpr auto tf = kLead.to_tf();                                             // C(s) = K (s + z)/(s + p)
 * static constinit LeadLagController<float> c{kLead};
 * constexpr auto kBoth = design::lead_lag({.phase_max = 0.5, .wc = 100.0, .dc_gain_boost = 10.0}, 1e-3,
 *                                         design::InputLimits{.lower = {-1.0}, .upper = {1.0}});
 * static constinit LeadLagController<float, 2> both{kBoth};                     // the lead, then the lag
 * @endcode
 *
 * @defgroup lead_lag Lead-Lag Compensator
 * @brief Lead/lag compensator design and discrete runtime for classical loop shaping
 *
 * A lead-lag section has the transfer function:
 *
 *   C(s) = K * (s + z) / (s + p)
 *
 * - Lead (z < p): Adds phase at a target frequency. Used to increase
 *   phase margin and improve transient response. Analogous to adding
 *   derivative action in a bounded frequency range.
 *
 * - Lag (z > p): Adds low-frequency gain without disturbing crossover.
 *   Used to reduce steady-state error. Analogous to adding integral
 *   action in a bounded frequency range.
 *
 * Design functions compute zero/pole locations and gain from intuitive
 * specifications (desired phase boost, crossover frequency, etc.); lead_lag puts a lead and a lag
 * in series, C(s) = Π K_i(s + z_i)/(s + p_i), which the same runtime runs section by section.
 *
 * Common applications:
 * - Voltage/current loop compensation in DC-DC converters
 * - Phase margin improvement in feedback loops
 * - Gain shaping for disturbance rejection
 * - Cascade compensator design (lead + lag sections)
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/lead_lag_stage.hpp"
#include "damp/controllers/law.hpp" // IWYU pragma: export
#include "damp/design/detail/lead_lag_rules.hpp"
#include "damp/design/input_limits.hpp" // IWYU pragma: export
#include "damp/design/options.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/toolbox/bounds.hpp"

namespace damp {

namespace design {

/**
 * @struct LeadLagGains
 * @brief One section of a lead-lag law, C(s) = K(s + z)/(s + p): what design::lead_lag_direct takes,
 *        a result holds per section and the runtime retunes (LeadLagController::gains, set_gains)
 */
template<typename T = double>
struct LeadLagGains {
    T K{}; ///< High-frequency gain (finite); the DC gain is K·z/p
    T z{}; ///< [rad/s] Zero location, positive: the zero is at s = −z
    T p{}; ///< [rad/s] Pole location, positive: the pole is at s = −p

    /** @brief The same section in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr LeadLagGains<U> as() const {
        return LeadLagGains<U>{static_cast<U>(K), static_cast<U>(z), static_cast<U>(p)};
    }

    /** @brief Continuous-time transfer function C(s) = K(s + z)/(s + p): num = {K, K·z}, den = {1, p}. */
    [[nodiscard]] constexpr TransferFunction<2, 2, T> to_tf() const {
        return TransferFunction<2, 2, T>{
            .num = {K, K * z},
            .den = {T{1}, p},
        };
    }

    /** @brief Continuous-time state-space, controllable canonical form: A = −p, B = 1, C = K(z − p), D = K. */
    [[nodiscard]] constexpr StateSpace<1, 1, 1, T> to_ss() const {
        return StateSpace{
            .A = Matrix<1, 1, T>{{-p}},
            .B = Matrix<1, 1, T>{{T{1}}},
            .C = Matrix<1, 1, T>{{K * (z - p)}},
            .D = Matrix<1, 1, T>{{K}},
        };
    }
};

/** @brief design::lead's spec: the phase to add and where. */
template<typename T = double>
struct LeadSpec {
    T phase_max{}; ///< [rad] Maximum phase boost, in (0, π/2)
    T wc{};        ///< [rad/s] Crossover frequency where the boost is (> 0)
};

/** @brief design::lag's spec: the DC gain to add, and how far below the crossover. */
template<typename T = double>
struct LagSpec {
    T dc_gain_boost{};      ///< DC gain increase factor (> 1, e.g. 10 for +20 dB)
    T wc{};                 ///< [rad/s] Crossover frequency to stay below (> 0)
    T margin_factor{T{10}}; ///< How far below ωc the section's zero sits (> 0); larger reduces the phase lag at ωc
};

/** @brief design::lead_lag's spec: a lead's and a lag's, at the same crossover. */
template<typename T = double>
struct LeadLagSpec {
    T phase_max{};          ///< [rad] The lead section's maximum phase boost, in (0, π/2)
    T wc{};                 ///< [rad/s] Crossover frequency (> 0)
    T dc_gain_boost{};      ///< The lag section's DC gain increase factor (> 1)
    T margin_factor{T{10}}; ///< How far below ωc the lag section's zero sits (> 0)
};

/**
 * @struct LeadLagLimits
 * @brief The limits of a lead-lag law: on its output
 *
 * Unbounded by default. The design sets the ones a law starts with from its options (InputLimits),
 * and the law reads and changes them at run time with LeadLagController::limits and
 * LeadLagController::set_limits.
 */
template<typename T = double>
struct LeadLagLimits {
    Bounds<1, T> input{}; ///< The output u [lower, upper]

    /** @brief The same limits in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr LeadLagLimits<U> as() const {
        return LeadLagLimits<U>{input.template as<U>()};
    }
};

/**
 * @struct LeadLagResult
 * @brief A lead-lag design deployed at a sample period: its sections in series, the limits the law
 *        starts with, and T_s
 *
 * C(s) = Π K_i(s + z_i)/(s + p_i), section 0 first. LeadLagController runs each section by Tustin
 * at T_s, which is @ref to_discrete_ss: the map of a product is the product of the maps.
 *
 * @tparam T  Scalar
 * @tparam NS Sections in series (1: a lead or a lag; 2: design::lead_lag's lead, then its lag)
 */
template<typename T = double, size_t NS = 1>
struct [[nodiscard]] LeadLagResult {
    static_assert(NS >= 1, "A lead-lag law has at least one section.");
    static constexpr size_t sections = NS; ///< The first-order sections in series

    damp::array<LeadLagGains<T>, NS> gains{};                             ///< Each section's K, z, p, in series from the input
    LeadLagLimits<T>                 limits{};                            ///< The limits the law starts with
    T                                Ts{};                                ///< [s] Sample period of the fixed-rate tick
    bool                             success{false};                      ///< status == DesignStatus::ok
    DesignStatus                     status{DesignStatus::not_converged}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The same design in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr LeadLagResult<U, NS> as() const {
        LeadLagResult<U, NS> out{.limits = limits.template as<U>(), .Ts = static_cast<U>(Ts), .success = success, .status = status};
        for (size_t i = 0; i < NS; ++i) {
            out.gains[i] = gains[i].template as<U>();
        }
        return out;
    }

    /** @brief Continuous-time transfer function Π K_i(s + z_i)/(s + p_i), highest power of s first. */
    [[nodiscard]] constexpr TransferFunction<NS + 1, NS + 1, T> to_tf() const {
        TransferFunction<NS + 1, NS + 1, T> tf{};
        tf.num[0] = T{1};
        tf.den[0] = T{1};
        for (size_t i = 0; i < NS; ++i) { // multiply by K_i s + K_i z_i and by s + p_i, i + 1 coefficients so far
            for (size_t j = i + 1; j > 0; --j) {
                tf.num[j] = (gains[i].K * tf.num[j]) + (gains[i].K * gains[i].z * tf.num[j - 1]);
                tf.den[j] = tf.den[j] + (gains[i].p * tf.den[j - 1]);
            }
            tf.num[0] = gains[i].K * tf.num[0];
        }
        return tf;
    }

    /** @brief Continuous-time state-space of the series, section 0 first (1 state per section). */
    [[nodiscard]] constexpr StateSpace<NS, 1, 1, T> to_ss() const { return chain<NS - 1>(); }

    /**
     * @brief Discrete-time state-space: the continuous series mapped at T_s
     *
     * Tustin (the default) is the runtime's model tick for tick. A failed design is the zero system.
     *
     * @param method Discretization method (default: Tustin)
     */
    [[nodiscard]] constexpr StateSpace<NS, 1, 1, T> to_discrete_ss(DiscretizationMethod method = DiscretizationMethod::Tustin) const {
        return *discretize(to_ss(), Ts, method);
    }

private:
    /** @brief Sections 0 … I in series. */
    template<size_t I>
    [[nodiscard]] constexpr StateSpace<I + 1, 1, 1, T> chain() const {
        if constexpr (I == 0) {
            return gains[0].to_ss();
        } else {
            return *::damp::series(chain<I - 1>(), gains[I].to_ss());
        }
    }
};

namespace detail {

using ::damp::detail::lead_lag_coefficients;
using ::damp::detail::lead_lag_section_ok;

/** @brief The options a lead-lag design takes. */
using lead_lag_options = option_set<kind::InputLimits>;

/**
 * @brief A lead-lag design once its sections are placed: the options checked, then @p sections
 *        deployed at @p Ts, refused when a section or a limit is out of its domain
 */
template<typename T, size_t NS, typename... Opts>
[[nodiscard]] constexpr LeadLagResult<T, NS> lead_lag_deployed(const damp::array<ZeroPoleGain<T>, NS>& sections, T Ts, const Opts&... opts) {
    lead_lag_options::check<Opts...>();
    using Result = LeadLagResult<T, NS>;
    if constexpr (option_count<kind::InputLimits, Opts...> > 0) {
        static_assert(option_t<kind::InputLimits, void, Opts...>::inputs == 1, "A lead-lag law has one output: InputLimits{.lower = {lo}, .upper = {hi}}, one entry at each end.");
    }
    if (!damp::finite_positive(Ts)) {
        return refused<Result>(DesignStatus::invalid_sample_time);
    }
    Result out{.Ts = Ts};
    if constexpr (option_t<kind::InputLimits, InputLimits<1, T>, Opts...>::inputs == 1) {
        out.limits.input = input_limits_of<1, T>(opts...); // unbounded without one
    }
    for (size_t i = 0; i < NS; ++i) {
        out.gains[i] = LeadLagGains<T>{.K = sections[i].K, .z = sections[i].z, .p = sections[i].p};
        const bool in_domain = sections[i].ok && lead_lag_section_ok(sections[i].K, sections[i].z, sections[i].p);
        if (!in_domain || !lead_lag_coefficients(sections[i].K, sections[i].z, sections[i].p, Ts).finite()) {
            return refused<Result>(DesignStatus::invalid_parameter);
        }
    }
    if (!out.limits.input.valid()) {
        return refused<Result>(DesignStatus::invalid_parameter);
    }
    set_status(out, DesignStatus::ok);
    return out;
}

} // namespace detail

/**
 * @brief A lead section from its phase boost at a crossover, deployed at @p Ts
 *
 * The maximum phase advance φ_max is at ωc, where the gain is 1 (detail::lead_section):
 * `design::lead({.phase_max = 0.6, .wc = 20.0}, 1e-3)`.
 *
 * @see "Feedback Control of Dynamic Systems" (Franklin, Powell & Emami-Naeini), lead compensation
 * @see lag() for the lag dual; lead_lag() for a lead and a lag in series
 *
 * @param spec φ_max in (0, π/2) and ωc > 0
 * @param Ts   [s] The law's sample period (finite and positive)
 * @param opts InputLimits
 * @return The design; refused with invalid_sample_time for a Ts that is not finite and positive,
 *         and with invalid_parameter for a spec or a limit out of its domain
 */
template<typename T = double, typename... Opts>
[[nodiscard]] constexpr LeadLagResult<T, 1> lead(const LeadSpec<T>& spec, T Ts, const Opts&... opts) {
    return detail::lead_lag_deployed<T, 1>({detail::lead_section(spec.phase_max, spec.wc)}, Ts, opts...);
}

/**
 * @brief A lag section from its DC gain boost below a crossover, deployed at @p Ts
 *
 * z = ωc/margin_factor, p = z/dc_gain_boost, K = 1 (detail::lag_section):
 * `design::lag({.dc_gain_boost = 10.0, .wc = 20.0}, 1e-3)`.
 *
 * @see "Feedback Control of Dynamic Systems" (Franklin, Powell & Emami-Naeini), lag compensation
 * @see lead() for phase boost at crossover
 *
 * @param spec A boost above 1, ωc > 0 and a margin factor > 0 (default 10)
 * @param Ts   [s] The law's sample period (finite and positive)
 * @param opts InputLimits
 * @return The design; refused with invalid_sample_time for a Ts that is not finite and positive,
 *         and with invalid_parameter for a spec or a limit out of its domain
 */
template<typename T = double, typename... Opts>
[[nodiscard]] constexpr LeadLagResult<T, 1> lag(const LagSpec<T>& spec, T Ts, const Opts&... opts) {
    return detail::lead_lag_deployed<T, 1>({detail::lag_section(spec.dc_gain_boost, spec.wc, spec.margin_factor)}, Ts, opts...);
}

/**
 * @brief A lead and a lag in series at one crossover, deployed at @p Ts
 *
 *   C(s) = K_lead (s + z_lead)/(s + p_lead) · (s + z_lag)/(s + p_lag)
 *
 * Section 0 is lead() of the phase boost, section 1 lag() of the DC boost; the runtime clamps the
 * lag's output, whose slow pole is the one that would wind up behind a clamp.
 *
 * @code
 * constexpr auto kLoop = design::lead_lag({.phase_max = 0.6, .wc = 20.0, .dc_gain_boost = 10.0}, 1e-3);
 * static constinit LeadLagController<float, 2> loop{kLoop};
 * @endcode
 *
 * @param spec φ_max in (0, π/2), ωc > 0, a DC boost above 1 and a margin factor > 0 (default 10)
 * @param Ts   [s] The law's sample period (finite and positive)
 * @param opts InputLimits
 * @return The design; refused with invalid_sample_time for a Ts that is not finite and positive,
 *         and with invalid_parameter for a spec or a limit out of its domain
 */
template<typename T = double, typename... Opts>
[[nodiscard]] constexpr LeadLagResult<T, 2> lead_lag(const LeadLagSpec<T>& spec, T Ts, const Opts&... opts) {
    return detail::lead_lag_deployed<T, 2>(
        {detail::lead_section(spec.phase_max, spec.wc), detail::lag_section(spec.dc_gain_boost, spec.wc, spec.margin_factor)}, Ts, opts...
    );
}

/**
 * @brief A section from the zero, pole and gain you already have, deployed at @p Ts
 *
 * `design::lead_lag_direct({.K = 2.0, .z = 100.0, .p = 500.0}, 1e-3)`.
 *
 * @param gains K finite, z and p finite and positive
 * @param Ts    [s] The law's sample period (finite and positive)
 * @param opts  InputLimits
 * @return The design; refused with invalid_sample_time for a Ts that is not finite and positive,
 *         and with invalid_parameter for a section or a limit out of its domain
 */
template<typename T = double, typename... Opts>
[[nodiscard]] constexpr LeadLagResult<T, 1> lead_lag_direct(const LeadLagGains<T>& gains, T Ts, const Opts&... opts) {
    return detail::lead_lag_deployed<T, 1>({detail::ZeroPoleGain<T>{.K = gains.K, .z = gains.z, .p = gains.p, .ok = true}}, Ts, opts...);
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Discrete lead-lag compensator: NS first-order sections in series, each K(s + z)/(s + p) by
 *        Tustin at the design's T_s, the gains and limits settable at run time
 *
 * | Member | What it is |
 * | --- | --- |
 * | step(e), step(r, y) | one tick at the design's T_s on the error e = r − y |
 * | track(u, e), track(u, r, y) | in place of each step while another source drives the plant |
 * | back_calculate(u_sat, u_ff) | a clamp after the law: the output that went out |
 * | gains(), set_gains(g) | each section's K, z, p (design::LeadLagGains); change them at run time |
 * | limits(), set_limits(lim) | the output's limits (design::LeadLagLimits); change them at run time |
 * | last_command(), sample_time(), valid(), reset() | the last output; T_s; whether the design ran; clear the memory |
 *
 * A failed design, a T_s that is not positive, a section out of its domain or a limit that is not a
 * box is refused at construction: ill-formed under constant evaluation; at run time it constructs
 * an invalid controller that emits zero and holds state (fail-closed).
 *
 * The output is clamped into limits().input, and the last section's recursion runs on that clamped
 * output: a lag section's slow pole then holds against the clamp where it kept charging towards
 * its DC gain behind it. For a stage that clamps after this one, back_calculate(u_sat) puts the
 * recursion on the output that went out. The sections before the last run on their own outputs.
 *
 * Gains and limits are tunable at run time, get-modify-set, and take effect at the next tick; a
 * retune keeps each section's memory:
 * @code
 * auto g = ll.gains();
 * g[0].z = 12.0f;
 * if (!ll.set_gains(g)) { } // refused: nothing changed
 * @endcode
 *
 * @see "Computer-Controlled Systems" (Åström & Wittenmark, 1997): the controller's recursion
 *      driven by the saturated output
 *
 * @tparam T     Scalar type (default: float)
 * @tparam NS    Sections in series (1: a lead or a lag; 2: design::lead_lag's)
 * @tparam Gains Where the sections live: void, the default, for members (detail::LeadLagGainSet, settable), or
 *               the constants of a `constexpr` design's Law (detail::ConstantLeadLagGains)
 */
template<typename T = float, size_t NS = 1, typename Gains = void>
class LeadLagController {
    static_assert(NS >= 1, "A lead-lag law has at least one section.");
    using GainSet = damp::conditional_t<damp::is_same_v<Gains, void>, detail::LeadLagGainSet<T, NS>, Gains>; ///< The sections' storage
    static constexpr bool kHeld = !GainSet::constant;                                                        ///< The sections are members: set_gains changes them

    template<typename, size_t, typename>
    friend class LeadLagController;

public:
    constexpr LeadLagController() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr LeadLagController(const design::LeadLagResult<U, NS>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : LeadLagController(result.template as<T>()) {}

    /** @brief From a design result. Non-explicit for assign-from-design: `LeadLagController<float> c = design::lead(spec, Ts).as<float>();` */
    constexpr LeadLagController(const design::LeadLagResult<T, NS>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : u_lim_(result.limits.input) {
        if (!result.success && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(result.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
        if (!damp::finite_positive(result.Ts) && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(design::DesignStatus::invalid_sample_time);
        }
        valid_ = result.success && damp::finite_positive(result.Ts) && result.limits.input.valid() && take(result.gains, result.Ts);
        if (!valid_ && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(design::DesignStatus::invalid_parameter); // a section out of its domain, or a limit that is not a box
        }
    }

    /** @brief The same law, with its live memory, in another scalar (or off a Law's constant sections, onto members). */
    template<typename U, typename G>
        requires kHeld
    constexpr explicit LeadLagController(const LeadLagController<U, NS, G>& o) : u_lim_(o.u_lim_.template as<T>()), u_cmd_(static_cast<T>(o.u_cmd_)) {
        damp::array<design::LeadLagGains<T>, NS> g{};
        for (size_t i = 0; i < NS; ++i) {
            stages_[i] = detail::LeadLagStage<T>(o.stages_[i]);
            g[i] = o.gains()[i].template as<T>();
        }
        valid_ = o.valid_ && take(g, static_cast<T>(o.sample_time()));
    }

    /** @brief False when built from a design it refuses: step() returns 0. */
    [[nodiscard]] constexpr bool valid() const { return valid_; }

    /**
     * @brief One tick on the compensator's input @p e (typically the tracking error)
     * @return The output, clamped into limits().input (zero if the design was invalid)
     */
    [[nodiscard]] constexpr T step(T e) {
        if (!valid_) {
            return T{0};
        }
        const T v = into_last(e);
        auto&   last = stages_[NS - 1];
        const T u = u_lim_.saturate(last.raw(v, k_.sections[NS - 1].k));
        last.e_prev = v;
        last.u_prev = u;
        u_cmd_ = u;
        return u;
    }

    /** @brief One tick on the reference @p r and the measurement @p y: step(r − y), so the law drops into a Cascade or any (r, y) consumer. */
    [[nodiscard]] constexpr T step(T r, T y) { return step(r - y); }

    /**
     * @brief In place of step(@p e) while another source drives the plant: the memory that starts the next output from @p u, the output applied
     *
     * The sections before the last run on e, so their memory is current; the last one's memory is
     * preloaded so that step(@p e) at this same e returns @p u clamped into limits().input: its
     * previous input is this tick's, and its previous output the value that, fed back, makes the
     * next output u. A last section with its pole at z = 0 (p = 2/T_s) has no memory of its output
     * to preload: the next step is the law's own output. @p u is recorded as the output applied, so
     * a back_calculate of it finds no clamp.
     */
    constexpr void track(T u, T e) {
        if (!valid_) {
            return;
        }
        T v = e;     // the last section's input on this tick
        T ahead = e; // and on the next one, at the same e
        for (size_t i = 0; i + 1 < NS; ++i) {
            v = stages_[i](v, k_.sections[i].k);
            ahead = stages_[i].raw(ahead, k_.sections[i].k);
        }
        stages_[NS - 1].preload(v, ahead, u_lim_.saturate(u), k_.sections[NS - 1].k);
        u_cmd_ = u;
    }

    /** @brief In place of step(@p r, @p y): track(@p u, @p r − @p y). */
    constexpr void track(T u, T r, T y) { track(u, r - y); }

    /**
     * @brief Record the output that went out, @p u_sat, when a stage after this one clamped it,
     *        and the feedforward @p u_ff added to this output before that clamp
     *
     * When u_sat − u_ff is not this law's last output (beyond rounding), the last section's
     * recursion continues from u_sat − u_ff, this compensator's share of what was applied, so it
     * holds against that clamp as against its own. Without a clamp nothing changes.
     */
    constexpr void back_calculate(T u_sat, T u_ff = T{0}) {
        if (!valid_) {
            return;
        }
        const T applied = u_sat - u_ff;
        if (damp::abs(applied - u_cmd_) > (T{4} * eps<T>() * (damp::abs(applied) + damp::abs(u_cmd_)))) {
            stages_[NS - 1].u_prev = applied;
        }
    }

    /** @brief The output the last step() returned, or the output the last track() recorded. */
    [[nodiscard]] constexpr T last_command() const { return u_cmd_; }

    /** @brief Clear every section's memory and the last output; keeps the gains and the limits. */
    constexpr void reset() {
        for (size_t i = 0; i < NS; ++i) {
            stages_[i].reset();
        }
        u_cmd_ = T{0};
    }

    /** @brief Each section's K, z, p in force: the design's, or the last set_gains. */
    [[nodiscard]] constexpr damp::array<design::LeadLagGains<T>, NS> gains() const {
        damp::array<design::LeadLagGains<T>, NS> g{};
        for (size_t i = 0; i < NS; ++i) {
            g[i] = design::LeadLagGains<T>{.K = k_.sections[i].K, .z = k_.sections[i].z, .p = k_.sections[i].p};
        }
        return g;
    }

    /**
     * @brief Change the sections at run time; the next tick uses them
     *
     * Recomputes the Tustin coefficients at sample_time() and keeps each section's memory.
     * Refused (false, nothing changed) for a K that is not finite, a z or p that is not finite and
     * positive, or a coefficient that overflows. A `constexpr` design's Law, whose sections are the
     * type's, has none.
     */
    [[nodiscard]] constexpr bool set_gains(const damp::array<design::LeadLagGains<T>, NS>& g)
        requires kHeld
    {
        return take(g, k_.Ts);
    }

    /** @brief The limits in force: the design's, or the last set_limits. */
    [[nodiscard]] constexpr design::LeadLagLimits<T> limits() const { return design::LeadLagLimits<T>{u_lim_}; }

    /**
     * @brief Change the limits at run time; the next tick clamps into them
     *
     * Refused (false, nothing changed) when an end is NaN or ±Inf, or the lower end is above the upper one.
     */
    [[nodiscard]] constexpr bool set_limits(const design::LeadLagLimits<T>& lim) {
        if (!lim.input.valid()) {
            return false;
        }
        u_lim_ = lim.input;
        return true;
    }

    /** @brief [s] The tick's sample period T_s. */
    [[nodiscard]] constexpr T sample_time() const { return k_.Ts; }

private:
    /** @brief The sections before the last, one tick on @p e: the last section's input. */
    [[nodiscard]] constexpr T into_last(T e) {
        T v = e;
        for (size_t i = 0; i + 1 < NS; ++i) {
            v = stages_[i](v, k_.sections[i].k);
        }
        return v;
    }

    /**
     * @brief Store @p g at the period @p Ts when every section is in its domain, keeping each
     *        section's memory; nothing changed otherwise. On a `constexpr` design's Law, whose
     *        sections are the type's, whether they are.
     */
    [[nodiscard]] constexpr bool take(const damp::array<design::LeadLagGains<T>, NS>& g, T Ts) {
        if constexpr (kHeld) {
            const auto set = detail::lead_lag_gain_set(g, Ts);
            if (!set) {
                return false;
            }
            k_ = *set;
            return true;
        } else {
            static_cast<void>(g);
            static_cast<void>(Ts);
            return GainSet::ok; // checked once, when the type was formed
        }
    }

    damp::array<detail::LeadLagStage<T>, NS> stages_{};     ///< Each section's memory, in series
    [[no_unique_address]] GainSet            k_{};          ///< Each section's K, z, p, T_s and the coefficients at it: members, or no storage on a constexpr design's Law
    Bounds<1, T>                             u_lim_{};      ///< The output's limits
    T                                        u_cmd_{T{0}};  ///< The output step() returned, or the one track() recorded
    bool                                     valid_{false}; ///< False when built from a design it refuses: step() returns 0
};

/** @brief Deduce the runtime from its design, keeping the result's scalar and sections. */
template<typename T, size_t NS>
LeadLagController(const design::LeadLagResult<T, NS>&) -> LeadLagController<T, NS>;

namespace detail {

/** @brief Whether @p R is a lead-lag design result. */
template<typename R>
inline constexpr bool is_lead_lag_result = false;

/** @brief design::lead's, lag's, lead_lag's and lead_lag_direct's result. */
template<typename U, size_t NS>
inline constexpr bool is_lead_lag_result<design::LeadLagResult<U, NS>> = true;

/** @brief A lead-lag design. */
template<const auto& Res>
concept lead_lag_design = is_lead_lag_result<design_t<Res>>;

} // namespace detail

/** @brief Law<kRes, T> of a lead-lag design initialized at start-up: the LeadLagController it runs as, with the design's sections. A runtime of the same shape assigned replaces it. */
template<const auto& Res, typename T>
    requires detail::lead_lag_design<Res>
struct Law<Res, T> : LeadLagController<T, detail::design_t<Res>::sections> {
    using Runtime = LeadLagController<T, detail::design_t<Res>::sections>; ///< The runtime the design runs as

    /** @brief The design's law. */
    constexpr Law() : Runtime{Res} {}

    using Runtime::operator=;
};

/**
 * @brief Law<kRes, T> of a `constexpr` lead-lag design: the LeadLagController of its sections with
 *        their gains constants of the type (detail::ConstantLeadLagGains)
 *
 * Each section's K, z, p, T_s and its Tustin coefficients b₀, b₁, a₁ are in flash and fold into
 * the tick; the sections' memory and the limits stay members, the limits settable at run time;
 * there is no set_gains. A law retuned at run time is `Law<kRes, T>::Runtime law{kRes.as<T>()};`.
 * A failed design, or one whose sections its runtime refuses at T_s, does not compile wherever the
 * Law is declared.
 */
template<const auto& Res, typename T>
    requires(detail::lead_lag_design<Res> && detail::constant_design<Res>)
struct Law<Res, T> : LeadLagController<T, detail::design_t<Res>::sections, detail::ConstantLeadLagGains<Res, T>> {
    using Runtime = LeadLagController<T, detail::design_t<Res>::sections>; ///< The runtime with members: the type of a member a redesign assigns
    using Base = LeadLagController<T, detail::design_t<Res>::sections, detail::ConstantLeadLagGains<Res, T>>;

    static_assert(detail::constant_design_runs(Res.success, Res.status, Res.Ts));
    static_assert(detail::ConstantLeadLagGains<Res, T>::ok, "Law<kRes, T>: a section of the lead-lag design is outside its domain, or a coefficient overflows at its Ts");

    /** @brief The design's law: its sections the type's, its limits members. */
    constexpr Law() : Base{Res} {}

    /** @brief There is nothing to replace: the sections are in the type. */
    template<typename Other>
        requires(!damp::is_same_v<damp::remove_cvref_t<Other>, Law>)
    constexpr Law& operator=(const Other& other) {
        static_assert(sizeof(Other) == 0, "This Law is a constexpr lead-lag design's: its sections are constants of the type, so there is nothing to replace. A member a redesign assigns, or a law retuned with set_gains, is Law<kRes, T>::Runtime law{kRes.as<T>()};");
        static_cast<void>(other);
        return *this;
    }
};

} // namespace damp
