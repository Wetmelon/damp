// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once


/**
 * @file lead_lag.hpp
 * @brief Lead-lag compensator design and runtime
 *
 * @code
 * using namespace damp;
 * constexpr auto lead = design::lead(0.5, 100.0, 0.0); // phi_max [rad], wc [rad/s], Ts [s] (0 = continuous)
 * static_assert(lead.success);
 * constexpr auto tf = lead.to_tf(); // C(s) = K (s+z)/(s+p)
 * static constinit LeadLagController<float> c{design::lead(0.5, 100.0, 1e-3).as<float>()}; // the runtime needs Ts > 0
 * @endcode
 *
 * @defgroup lead_lag Lead-Lag Compensator
 * @brief Lead/lag compensator design and discrete runtime for classical loop shaping
 *
 * A lead-lag compensator has the transfer function:
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
 * specifications (desired phase boost, crossover frequency, etc.).
 *
 * Common applications:
 * - Voltage/current loop compensation in DC-DC converters
 * - Phase margin improvement in feedback loops
 * - Gain shaping for disturbance rejection
 * - Cascade compensator design (lead + lag sections)
 */

#include "damp/backend.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

namespace damp {

namespace design {

/**
 * @struct LeadLagResult
 * @brief Lead-lag compensator design result
 */
template<typename T = double>
struct [[nodiscard]] LeadLagResult {
    T            K{};                                 ///< High-frequency gain (DC gain is K·z/p)
    T            z{};                                 ///< Zero location (rad/s, positive value; actual zero at s = -z)
    T            p{};                                 ///< Pole location (rad/s, positive value; actual pole at s = -p)
    T            Ts{};                                ///< Sampling time (seconds), 0 = continuous
    bool         success{false};                      ///< status == DesignStatus::ok
    DesignStatus status{DesignStatus::not_converged}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr LeadLagResult<U> as() const {
        return {static_cast<U>(K), static_cast<U>(z), static_cast<U>(p), static_cast<U>(Ts), success, status};
    }

    /**
     * @brief Convert to continuous-time transfer function
     *
     * Returns C(s) = K * (s + z) / (s + p)
     * Numerator and denominator highest power of s first: {s, const}
     */
    [[nodiscard]] constexpr TransferFunction<2, 2, T> to_tf() const {
        // C(s) = K * (s + z) / (s + p)
        // Numerator:   K*s + K*z   = {K, K*z}
        // Denominator: s + p       = {1, p}
        return TransferFunction<2, 2, T>{
            .num = {K, K * z},
            .den = {T{1}, p},
        };
    }

    /**
     * @brief Convert to continuous-time state-space (1st order SISO)
     */
    [[nodiscard]] constexpr StateSpace<1, 1, 1, T> to_ss() const {
        // From C(s) = K*(s+z)/(s+p), controllable canonical form:
        // A = [-p], B = [1], C = [K*(z-p)], D = [K]
        return StateSpace<1, 1, 1, T>{
            .A = Matrix<1, 1, T>{{-p}},
            .B = Matrix<1, 1, T>{{T{1}}},
            .C = Matrix<1, 1, T>{{K * (z - p)}},
            .D = Matrix<1, 1, T>{{K}},
        };
    }

    /**
     * @brief Convert to discrete-time state-space
     *
     * Ts = 0 (continuous design) returns @ref to_ss.
     *
     * @param method Discretization method (default: Tustin)
     */
    [[nodiscard]] constexpr StateSpace<1, 1, 1, T>
    to_discrete_ss(DiscretizationMethod method = DiscretizationMethod::Tustin) const {
        return (Ts > T{0}) ? *discretize(to_ss(), Ts, method) : to_ss();
    }
};

/**
 * @brief Design a lead compensator from desired phase boost at a target frequency
 *
 * Places zero and pole symmetrically (in log-frequency) around wc so that
 * maximum phase advance occurs exactly at wc:
 *
 *     α = (1 − sin φ_max) / (1 + sin φ_max),   z = wc √α,   p = wc / √α
 *
 * with α = z/p < 1. Gain K = 1/√α so |C(j wc)| = 1 (phase boost without
 * shifting the gain-crossover frequency).
 *
 * @see Franklin, Powell & Emami-Naeini, "Feedback Control of Dynamic Systems"
 * @see lag() for the lag dual; lead_lag() for cascade lead+lag
 *
 * @param phi_max  Desired maximum phase boost (radians, 0 < phi_max < π/2)
 * @param wc       Target crossover frequency (rad/s), must be > 0
 * @param Ts       Sampling time (seconds); 0 = continuous-time design (analysis only: the runtime refuses it)
 * @return LeadLagResult with K, z, p; refused with invalid_parameter for φ_max outside (0, π/2) or a
 *         non-positive or non-finite wc, and with invalid_sample_time for a negative or non-finite Ts
 */
template<typename T = double>
[[nodiscard]] constexpr LeadLagResult<T> lead(T phi_max, T wc, T Ts) {
    // Valid lead requires 0 < phi_max < pi/2 (else alpha <= 0) and wc > 0.
    if (!damp::finite_positive(phi_max) || !(phi_max < damp::numbers::pi_v<T> / T{2}) || !damp::finite_positive(wc)) {
        return ::damp::detail::refused<LeadLagResult<T>>(DesignStatus::invalid_parameter);
    }
    if (!damp::finite_non_negative(Ts)) {
        return ::damp::detail::refused<LeadLagResult<T>>(DesignStatus::invalid_sample_time);
    }

    // alpha = (1 - sin(phi)) / (1 + sin(phi))
    // For a lead: alpha < 1
    T s = damp::sin(phi_max);
    T alpha = (T{1} - s) / (T{1} + s);

    // Zero and pole placed symmetrically around wc in log-frequency:
    // z = wc * sqrt(alpha),  p = wc / sqrt(alpha)
    // So geometric mean = wc, and z < p (lead).
    T sqrt_alpha = damp::sqrt(alpha);
    T z = wc * sqrt_alpha;
    T p = wc / sqrt_alpha;

    // Gain set for unity magnitude at wc: |C(j*wc)| = 1
    // |C(j*wc)| = K * sqrt(wc²+z²) / sqrt(wc²+p²)
    //           = K * sqrt(1 + alpha) / sqrt(1 + 1/alpha)
    //           = K * sqrt(alpha)
    // So K = 1/sqrt(alpha) for unity gain at crossover.
    T K = T{1} / sqrt_alpha;

    return LeadLagResult<T>{K, z, p, Ts, true, DesignStatus::ok};
}

/**
 * @brief Design a lag compensator from desired low-frequency gain boost
 *
 * Places the lag section a factor @p margin_factor below crossover so DC gain
 * rises by @p dc_gain_boost while phase lag at wc stays small:
 *
 *     z = wc / margin_factor,   p = z / dc_gain_boost,   K = 1
 *
 * so DC gain of K(s+z)/(s+p) equals dc_gain_boost (z > p).
 *
 * @see Franklin, Powell & Emami-Naeini, "Feedback Control of Dynamic Systems"
 * @see lead() for phase boost at crossover
 *
 * @param dc_gain_boost  Desired DC gain increase factor (> 1, e.g. 10 for +20 dB)
 * @param wc             Crossover frequency to stay below (rad/s), must be > 0
 * @param Ts             Sampling time (seconds); 0 = continuous (analysis only: the runtime refuses it)
 * @param margin_factor  How far below wc to place the compensator (default: 10).
 *                       Larger values reduce phase lag at wc.
 * @return LeadLagResult with K, z, p; refused with invalid_parameter for an out-of-range or
 *         non-finite spec, and with invalid_sample_time for a negative or non-finite Ts
 */
template<typename T = double>
[[nodiscard]] constexpr LeadLagResult<T> lag(T dc_gain_boost, T wc, T Ts, T margin_factor = T{10}) {
    // Valid lag requires a real DC boost, a positive crossover, and a positive margin.
    if (!damp::isfinite(dc_gain_boost) || !(dc_gain_boost > T{1}) || !damp::finite_positive(wc)
        || !damp::finite_positive(margin_factor)) {
        return ::damp::detail::refused<LeadLagResult<T>>(DesignStatus::invalid_parameter);
    }
    if (!damp::finite_non_negative(Ts)) {
        return ::damp::detail::refused<LeadLagResult<T>>(DesignStatus::invalid_sample_time);
    }

    // Place zero at wc / margin_factor
    // Place pole at zero / dc_gain_boost (further below, so z > p → lag)
    T z = wc / margin_factor;
    T p = z / dc_gain_boost;

    // DC gain of (s+z)/(s+p) → z/p = dc_gain_boost
    // Set K = 1 so the total DC gain is exactly dc_gain_boost
    T K = T{1};

    return LeadLagResult<T>{K, z, p, Ts, true, DesignStatus::ok};
}

/**
 * @struct LeadLagSeriesResult
 * @brief Cascaded lead+lag design result (2nd-order StateSpace + success)
 */
template<typename T = double>
struct [[nodiscard]] LeadLagSeriesResult {
    StateSpace<2, 1, 1, T> sys{};                               ///< Cascaded compensator (zero if design failed)
    bool                   success{false};                      ///< status == DesignStatus::ok
    DesignStatus           status{DesignStatus::not_converged}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr LeadLagSeriesResult<U> as() const {
        return {sys.template as<U>(), success, status};
    }
};

/**
 * @brief Design a lead-lag compensator (cascade of lead + lag sections)
 *
 * Returns a 2nd-order system combining a lead section (for phase margin)
 * and a lag section (for DC gain). The two sections are cascaded:
 *
 *   C(s) = K_lead * (s+z_lead)/(s+p_lead) * (s+z_lag)/(s+p_lag)
 *
 * Propagates section status: if either lead() or lag() fails, the result is
 * refused with that section's status and sys is the zero StateSpace.
 *
 * @param phi_max        Desired phase boost from lead section (radians)
 * @param wc             Target crossover frequency (rad/s)
 * @param dc_gain_boost  Desired DC gain increase from lag section (> 1)
 * @param Ts             Sampling time (seconds), 0 = continuous
 * @param margin_factor  How far below wc to place lag section (default: 10)
 * @return LeadLagSeriesResult with the 2nd-order SISO compensator
 */
template<typename T = double>
[[nodiscard]] constexpr LeadLagSeriesResult<T>
lead_lag(T phi_max, T wc, T dc_gain_boost, T Ts, T margin_factor = T{10}) {
    const auto lead_r = lead(phi_max, wc, Ts);
    const auto lag_r = lag(dc_gain_boost, wc, Ts, margin_factor);

    if (!lead_r.success) {
        return ::damp::detail::refused<LeadLagSeriesResult<T>>(lead_r.status);
    }
    if (!lag_r.success) {
        return ::damp::detail::refused<LeadLagSeriesResult<T>>(lag_r.status);
    }

    // Cascade: series connection of lead and lag state-space representations
    const auto combined = *series(lead_r.to_ss(), lag_r.to_ss());
    const auto sys = (Ts > T{0}) ? *discretize(combined, Ts, DiscretizationMethod::Tustin) : combined;

    return LeadLagSeriesResult<T>{
        .sys = sys,
        .success = true,
        .status = DesignStatus::ok,
    };
}

/**
 * @brief Direct lead-lag specification from zero/pole locations
 *
 * For users who already know their zero and pole locations.
 *
 * @param K   Gain
 * @param z   Zero location (positive; actual zero at s = -z)
 * @param p   Pole location (positive; actual pole at s = -p)
 * @param Ts  Sampling time (seconds), 0 = continuous
 * @return LeadLagResult; refused with invalid_parameter for a non-finite K or a non-positive or
 *         non-finite z or p, and with invalid_sample_time for a negative or non-finite Ts
 */
template<typename T = double>
[[nodiscard]] constexpr LeadLagResult<T> lead_lag_direct(T K, T z, T p, T Ts) {
    if (!damp::isfinite(K) || !damp::finite_positive(z) || !damp::finite_positive(p)) {
        return ::damp::detail::refused<LeadLagResult<T>>(DesignStatus::invalid_parameter);
    }
    if (!damp::finite_non_negative(Ts)) {
        return ::damp::detail::refused<LeadLagResult<T>>(DesignStatus::invalid_sample_time);
    }
    return LeadLagResult<T>{K, z, p, Ts, true, DesignStatus::ok};
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Discrete lead-lag compensator
 *
 * Implements @f$ C(s) = K(s+z)/(s+p) @f$ discretized via Tustin. A failed or continuous
 * (Ts = 0) design does not compile under constant evaluation; at run time it constructs an
 * invalid controller that emits zero and holds state (fail-closed).
 *
 * @tparam T Scalar type (default: float)
 */
template<typename T = float>
struct LeadLagController {
    T K{};
    T z{};
    T p{};

    // Pre-computed IIR coefficients (1st-order: u[n] = b0*e[n] + b1*e[n-1] - a1*u[n-1])
    T b0{0};
    T b1{0};
    T a1{0};

    // State
    T e_prev{0}; ///< previous input e[n-1]
    T u_prev{0}; ///< previous output u[n-1]

    constexpr LeadLagController() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr LeadLagController(const design::LeadLagResult<U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : LeadLagController(result.template as<T>()) {}

    constexpr LeadLagController(const design::LeadLagResult<T>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : K(result.K),
          z(result.z),
          p(result.p),
          valid_(result.success && result.Ts > T{0} && result.z > T{0} && result.p > T{0}) {
        if ((!result.success || !(result.Ts > T{0})) && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(result.success ? design::DesignStatus::invalid_sample_time : result.status); // constant-evaluated (constinit / constexpr): a failed or continuous design does not compile
        }
        if (valid_) {
            compute_coefficients(result.Ts);
        }
    }

    template<typename U>
    constexpr LeadLagController(const LeadLagController<U>& other)
        : K(static_cast<T>(other.K)),
          z(static_cast<T>(other.z)),
          p(static_cast<T>(other.p)),
          b0(static_cast<T>(other.b0)),
          b1(static_cast<T>(other.b1)),
          a1(static_cast<T>(other.a1)),
          e_prev(static_cast<T>(other.e_prev)),
          u_prev(static_cast<T>(other.u_prev)),
          valid_(other.valid()) {}

    /**
     * @brief Compute compensator output
     *
     * @param e Compensator input (typically the tracking error)
     * @return Compensated output (zero if the design was invalid)
     */
    [[nodiscard]] constexpr T control(T e) {
        if (!valid_) {
            return T{0};
        }
        const T u = (b0 * e) + (b1 * e_prev) - (a1 * u_prev);
        e_prev = e;
        u_prev = u;
        return u;
    }

    /**
     * @brief Reference-tracking overload satisfying SISOController.
     *
     * Lead-lag is structurally a SISO filter (frequency-domain compensator),
     * but it's most commonly applied to the tracking error in a feedback
     * loop. This overload computes the error internally so the controller
     * can drop into `Cascade<Outer, LeadLagController>` and any other
     * SISOController-shaped consumer.
     *
     * The single-argument `control(T e)` form is still available for use as
     * a standalone filter inside a larger composite controller.
     */
    [[nodiscard]] constexpr T control(T r, T y) {
        return control(r - y);
    }

    constexpr void reset() {
        e_prev = T{0};
        u_prev = T{0};
    }

    [[nodiscard]] constexpr bool valid() const { return valid_; }

private:
    bool valid_{true}; // default ctor leaves a zero-gain filter usable as pass-through placeholder
    /**
     * @brief Compute Tustin-discretized IIR coefficients
     *
     * C(s) = K*(s+z)/(s+p)
     * Tustin: s = (2/Ts)*(z-1)/(z+1)
     *
     * Numerator:   K * ((2/Ts)*(z-1)/(z+1) + z_loc)
     *            = K * ((2/Ts + z_loc)*z + (z_loc - 2/Ts)) / (z+1)
     *
     * Denominator: (2/Ts)*(z-1)/(z+1) + p_loc
     *            = ((2/Ts + p_loc)*z + (p_loc - 2/Ts)) / (z+1)
     */
    constexpr void compute_coefficients(T Ts) {
        T k = T{2} / Ts; // Tustin substitution factor

        // Denominator coefficients (before normalization)
        T a0_raw = k + p;
        T a1_raw = p - k;

        // Numerator coefficients (before normalization)
        T b0_raw = K * (k + z);
        T b1_raw = K * (z - k);

        // Normalize by a0
        b0 = b0_raw / a0_raw;
        b1 = b1_raw / a0_raw;
        a1 = a1_raw / a0_raw;
    }
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<typename T>
LeadLagController(const design::LeadLagResult<T>&) -> LeadLagController<T>;

} // namespace damp
