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
 * constexpr auto lead = design::lead(0.5, 100.0); // phi_max [rad], wc [rad/s]
 * static_assert(lead.success);
 * constexpr auto tf = lead.to_tf(); // C(s) = K (s+z)/(s+p)
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
struct LeadLagResult {
    T    K{};            ///< DC gain of the compensator
    T    z{};            ///< Zero location (rad/s, positive value; actual zero at s = -z)
    T    p{};            ///< Pole location (rad/s, positive value; actual pole at s = -p)
    T    Ts{};           ///< Sampling time (seconds), 0 = continuous
    bool success{false}; ///< true if the spec produced valid (finite, positive) z, p

    template<typename U>
    [[nodiscard]] constexpr LeadLagResult<U> as() const {
        return {static_cast<U>(K), static_cast<U>(z), static_cast<U>(p), static_cast<U>(Ts), success};
    }

    /**
     * @brief Convert to continuous-time transfer function
     *
     * Returns C(s) = K * (s + z) / (s + p)
     * Numerator and denominator in ascending powers of s: {const, s}
     */
    [[nodiscard]] constexpr TransferFunction<2, 2, T> to_tf() const {
        // C(s) = K * (s + z) / (s + p)
        // Numerator:   K*z + K*s   = {K*z, K}
        // Denominator: p + s       = {p, 1}
        return TransferFunction<2, 2, T>{
            .num = {K * z, K},
            .den = {p, T{1}},
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
     * @param method Discretization method (default: Tustin)
     */
    [[nodiscard]] constexpr StateSpace<1, 1, 1, T>
    to_discrete_ss(DiscretizationMethod method = DiscretizationMethod::Tustin) const {
        return *discretize(to_ss(), Ts, method);
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
 * @param Ts       Sampling time (seconds), 0 = continuous-time design
 * @return LeadLagResult with K, z, p; success false if the spec is out of range
 */
template<typename T = double>
[[nodiscard]] constexpr LeadLagResult<T> lead(T phi_max, T wc, T Ts = T{0}) {
    // Valid lead requires 0 < phi_max < pi/2 (else alpha <= 0) and wc > 0.
    if (!(phi_max > T{0}) || !(phi_max < damp::numbers::pi_v<T> / T{2}) || !(wc > T{0})) {
        return LeadLagResult<T>{};
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

    return LeadLagResult<T>{K, z, p, Ts, true};
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
 * @param margin_factor  How far below wc to place the compensator (default: 10).
 *                       Larger values reduce phase lag at wc.
 * @param Ts             Sampling time (seconds), 0 = continuous
 * @return LeadLagResult with K, z, p; success false if the spec is out of range
 */
template<typename T = double>
[[nodiscard]] constexpr LeadLagResult<T> lag(T dc_gain_boost, T wc, T margin_factor = T{10}, T Ts = T{0}) {
    // Valid lag requires a real DC boost, a positive crossover, and a positive margin.
    if (!(dc_gain_boost > T{1}) || !(wc > T{0}) || !(margin_factor > T{0})) {
        return LeadLagResult<T>{};
    }

    // Place zero at wc / margin_factor
    // Place pole at zero / dc_gain_boost (further below, so z > p → lag)
    T z = wc / margin_factor;
    T p = z / dc_gain_boost;

    // DC gain of (s+z)/(s+p) → z/p = dc_gain_boost
    // Set K = 1 so the total DC gain is exactly dc_gain_boost
    T K = T{1};

    return LeadLagResult<T>{K, z, p, Ts, true};
}

/**
 * @struct LeadLagSeriesResult
 * @brief Cascaded lead+lag design result (2nd-order StateSpace + success)
 */
template<typename T = double>
struct LeadLagSeriesResult {
    StateSpace<2, 1, 1, T> sys{};          ///< Cascaded compensator (zero if design failed)
    bool                   success{false}; ///< true if both lead and lag sections succeeded

    template<typename U>
    [[nodiscard]] constexpr LeadLagSeriesResult<U> as() const {
        return {sys.template as<U>(), success};
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
 * Propagates section success: if either lead() or lag() fails, success is
 * false and sys is the zero StateSpace.
 *
 * @param phi_max        Desired phase boost from lead section (radians)
 * @param wc             Target crossover frequency (rad/s)
 * @param dc_gain_boost  Desired DC gain increase from lag section (> 1)
 * @param margin_factor  How far below wc to place lag section (default: 10)
 * @param Ts             Sampling time (seconds), 0 = continuous
 * @return LeadLagSeriesResult with the 2nd-order SISO compensator
 */
template<typename T = double>
[[nodiscard]] constexpr LeadLagSeriesResult<T>
lead_lag(T phi_max, T wc, T dc_gain_boost, T margin_factor = T{10}, T Ts = T{0}) {
    const auto lead_r = lead(phi_max, wc, Ts);
    const auto lag_r = lag(dc_gain_boost, wc, margin_factor, Ts);

    if (!lead_r.success || !lag_r.success) {
        return LeadLagSeriesResult<T>{};
    }

    // Cascade: series connection of lead and lag state-space representations
    const auto combined = *series(lead_r.to_ss(), lag_r.to_ss());
    const auto sys = (Ts > T{0}) ? *discretize(combined, Ts, DiscretizationMethod::Tustin) : combined;

    return LeadLagSeriesResult<T>{
        .sys = sys,
        .success = true,
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
 * @return LeadLagResult
 */
template<typename T = double>
[[nodiscard]] constexpr LeadLagResult<T> lead_lag_direct(T K, T z, T p, T Ts = T{0}) {
    return LeadLagResult<T>{K, z, p, Ts, z > T{0} && p > T{0}};
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Discrete lead-lag compensator
 *
 * Implements @f$ C(s) = K(s+z)/(s+p) @f$ discretized via Tustin. A failed design
 * (@c success = false) constructs an invalid controller that emits zero and holds
 * state (fail-closed).
 *
 * @tparam T Scalar type (default: float)
 */
template<typename T = float>
struct LeadLagController {
    T K{};
    T z{};
    T p{};

    // Pre-computed IIR coefficients (1st-order: y[n] = b0*u[n] + b1*u[n-1] - a1*y[n-1])
    T b0{0};
    T b1{0};
    T a1{0};

    // State
    T u_prev{0};
    T y_prev{0};

    constexpr LeadLagController() = default;

    constexpr LeadLagController(const design::LeadLagResult<T>& result)
        : K(result.K), z(result.z), p(result.p), valid_(result.success && result.z > T{0} && result.p > T{0}) {
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
          u_prev(static_cast<T>(other.u_prev)),
          y_prev(static_cast<T>(other.y_prev)),
          valid_(other.valid()) {}

    /**
     * @brief Compute compensator output
     *
     * @param u Input signal (typically error signal)
     * @return Compensated output (zero if the design was invalid)
     */
    [[nodiscard]] constexpr T control(T u) {
        if (!valid_) {
            return T{0};
        }
        T y = b0 * u + b1 * u_prev - a1 * y_prev;
        u_prev = u;
        y_prev = y;
        return y;
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
     * The single-argument `control(T u)` form is still available for use as
     * a standalone filter inside a larger composite controller.
     */
    [[nodiscard]] constexpr T control(T r, T y) {
        return control(r - y);
    }

    constexpr void reset() {
        u_prev = T{0};
        y_prev = T{0};
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
        // Ts <= 0 means a continuous-time design (e.g. the design factories' default
        // Ts = 0) was handed to a discrete controller. There is no valid Tustin
        // mapping.
        if (!(Ts > T{0})) {
            b0 = K;
            b1 = T{0};
            a1 = T{0};
            return;
        }

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

} // namespace damp
