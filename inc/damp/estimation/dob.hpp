// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file dob.hpp
 * @brief Disturbance-observer estimation primitives and lightweight runtime.
 */

#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp" // damp::array
#include "damp/math/math.hpp"
#include "damp/matrix/matrix_traits.hpp"

namespace damp::estimation {

/**
 * @brief Configuration for a first-order disturbance observer.
 *
 * The runtime update is:
 *
 *     innovation = y_measured - y_predicted
 *     d_hat[k+1] = (1 - leak) * d_hat[k] + gain * innovation
 */
template<typename T = double>
struct DOBConfig {

    using scalar_type = scalar_type_t<T>;

    scalar_type gain{scalar_type{0.1}};
    scalar_type leak{scalar_type{0}};
    scalar_type innovation_deadband{scalar_type{0}};
    scalar_type max_disturbance_magnitude{scalar_type{0}};
    bool        clamp_enabled{false};

    [[nodiscard]] constexpr bool valid() const {
        if (gain < scalar_type{0}) {
            return false;
        }
        if (gain > scalar_type{1}) {
            return false;
        }
        if (leak < scalar_type{0}) {
            return false;
        }
        if (leak >= scalar_type{1}) {
            return false;
        }
        if (innovation_deadband < scalar_type{0}) {
            return false;
        }
        if (max_disturbance_magnitude < scalar_type{0}) {
            return false;
        }
        return true;
    }
};

template<typename T = double>
struct DOBResult {

    using scalar_type = scalar_type_t<T>;

    DOBConfig<T> config{};

    /// Constant-innovation steady-state gain d_hat_ss/innovation = gain/leak.
    /// Sentinel 0 means "no finite steady state": leak == 0 is a pure integrator
    /// whose estimate grows without bound under constant innovation.
    T    steady_state_gain{};
    bool success{false};

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        using out_t = std::remove_const_t<U>;
        return DOBResult<out_t>{
            DOBConfig<out_t>{
                static_cast<scalar_type_t<out_t>>(config.gain),
                static_cast<scalar_type_t<out_t>>(config.leak),
                static_cast<scalar_type_t<out_t>>(config.innovation_deadband),
                static_cast<scalar_type_t<out_t>>(config.max_disturbance_magnitude),
                config.clamp_enabled,
            },
            static_cast<out_t>(steady_state_gain),
            success,
        };
    }
};

template<typename T = float>
struct DOBState {

    T    disturbance_hat{};
    T    innovation{};
    bool initialized{false};
};

namespace detail {

template<typename T>
[[nodiscard]] constexpr T limit_magnitude(const T& value, scalar_type_t<T> max_magnitude) {
    if (max_magnitude <= scalar_type_t<T>{0}) {
        return value;
    }

    const auto magnitude = damp::abs(value);
    if (magnitude <= max_magnitude) {
        return value;
    }
    if (magnitude <= default_tol<T>()) {
        return T{};
    }

    const auto scale = max_magnitude / magnitude;
    return value * scale;
}

} // namespace detail

/**
 * @brief Validate and package DOB configuration into a runtime-ready design result.
 */
template<typename T = double>
[[nodiscard]] constexpr DOBResult<T> dob(
    const DOBConfig<T>& config
) {
    if (!config.valid()) {
        return DOBResult<T>{.config = config, .success = false};
    }

    const auto one_minus_leak = scalar_type_t<T>{1} - config.leak;
    if (one_minus_leak <= default_tol<T>()) {
        return DOBResult<T>{.config = config, .steady_state_gain = T{}, .success = true};
    }

    // Constant-innovation steady-state gain: d_hat_ss = gain / leak * innovation.
    const T steady_state_gain = (config.leak <= default_tol<T>()) ? T{} : static_cast<T>(config.gain / config.leak);
    return DOBResult<T>{
        .config = config,
        .steady_state_gain = steady_state_gain,
        .success = true,
    };
}

/**
 * @brief Lightweight SISO disturbance observer runtime.
 */
template<typename T = float>
class DOB {
public:
    using scalar_type = scalar_type_t<T>;

    constexpr DOB() = default;

    constexpr explicit DOB(const DOBConfig<T>& config)
        : config_(config) {}

    constexpr explicit DOB(const DOBResult<T>& design)
        : config_(design.config), valid_(design.success) {}

    /**
     * @brief Update disturbance estimate from predicted and measured outputs.
     * @return true when update succeeded.
     */
    [[nodiscard]] constexpr bool update(T y_predicted, T y_measured) {
        if (!config_.valid()) {
            valid_ = false;
            return false;
        }

        const scalar_type innovation_mag_deadband = config_.innovation_deadband;
        T                 innovation = y_measured - y_predicted;

        if (innovation_mag_deadband > scalar_type{0} && damp::abs(innovation) < innovation_mag_deadband) {
            innovation = T{};
        }

        const scalar_type alpha = scalar_type{1} - config_.leak;
        state_.disturbance_hat = (static_cast<T>(alpha) * state_.disturbance_hat) + (static_cast<T>(config_.gain) * innovation);

        if (config_.clamp_enabled && config_.max_disturbance_magnitude > scalar_type{0}) {
            state_.disturbance_hat = detail::limit_magnitude(state_.disturbance_hat, config_.max_disturbance_magnitude);
        }

        state_.innovation = innovation;
        state_.initialized = true;
        valid_ = true;
        return true;
    }

    /**
     * @brief Disturbance-compensated command (u = u_nominal - d_hat).
     */
    [[nodiscard]] constexpr T compensate(const T& u_nominal) const {
        return u_nominal - state_.disturbance_hat;
    }

    constexpr void reset() {
        state_ = DOBState<T>{};
        valid_ = true;
    }

    [[nodiscard]] constexpr const auto& config() const { return config_; }
    [[nodiscard]] constexpr const auto& state() const { return state_; }
    [[nodiscard]] constexpr bool        valid() const { return valid_; }

private:
    DOBConfig<T> config_{};
    DOBState<T>  state_{};
    bool         valid_{true};
};

// ===========================================================================
// Classical Pn^-1 * Q disturbance observer (Ohnishi DOB) — bolt-on for an
// existing controller. Estimates the input-referred lumped disturbance from a
// nominal plant model and a low-pass Q-filter, and subtracts it from the
// command. Complements the scalar innovation-based observer above (richer plant
// model) and ADRC's ESO (this is a feed-around bolt-on, not a full controller).
//
//   plant:        y = P·(u + d)          (d = input-referred disturbance)
//   estimate:     d_hat = Q·(Pn^-1·y − u) ≈ Q·d   (exact at DC when P = Pn)
//   compensate:   u = u_command − d_hat
//
// All polynomials are discrete z^-1 digital filters (c0 + c1·z^-1 + …), the
// natural form for an embedded loop. Realizability needs the leading numerator
// of Pn and the leading denominators nonzero (Pn must be causally invertible —
// no pure input delay in the nominal model).
//
// @see "Disturbance Observer-Based Control" (S. Li et al., CRC Press, 2016)
// @see "Robust Motion Control by Disturbance Observer" (Ohnishi et al., 1996)
// ===========================================================================

namespace detail {

/// Polynomial (z^-1) convolution = digital-filter series multiplication.
template<size_t Na, size_t Nb, typename T>
[[nodiscard]] constexpr damp::array<T, Na + Nb - 1> poly_conv(const damp::array<T, Na>& a, const damp::array<T, Nb>& b) {
    damp::array<T, Na + Nb - 1> r{};
    for (size_t i = 0; i < Na; ++i) {
        for (size_t j = 0; j < Nb; ++j) {
            r[i + j] += a[i] * b[j];
        }
    }
    return r;
}

/// Direct-Form-II IIR for a digital filter H(z^-1) = num / den (den[0] ≠ 0).
template<size_t Nnum, size_t Nden, typename T>
struct IirDF2 {
    static constexpr size_t Nstate = ((Nnum > Nden) ? Nnum : Nden) - 1;

    damp::array<T, Nnum>   num{};
    damp::array<T, Nden>   den{};
    damp::array<T, Nstate> w{}; // w[0] = most recent state

    constexpr T step(T x) {
        T w0 = x;
        for (size_t i = 1; i < Nden; ++i) {
            w0 -= den[i] * w[i - 1];
        }
        w0 /= den[0];
        T y = num[0] * w0;
        for (size_t i = 1; i < Nnum; ++i) {
            y += num[i] * w[i - 1];
        }
        for (size_t i = Nstate; i-- > 1;) {
            w[i] = w[i - 1];
        }
        if constexpr (Nstate > 0) {
            w[0] = w0;
        }
        return y;
    }

    constexpr void reset() { w = damp::array<T, Nstate>{}; }
};

} // namespace detail

/**
 * @brief Design result for the classical Pn^-1·Q disturbance observer.
 *
 * Holds the two realized digital filters: Fy = Q·Pn^-1 (applied to the
 * measurement y) and Fu = Q (applied to the applied input u). The disturbance
 * estimate is d_hat = Fy(y) − Fu(u).
 *
 * @tparam NBn,NAn Nominal plant Pn = Bn/An sizes (z^-1)
 * @tparam NQn,NQd Q-filter = Qn/Qd sizes (z^-1)
 */
template<size_t NBn, size_t NAn, size_t NQn, size_t NQd, typename T = double>
struct ClassicalDobResult {

    static constexpr size_t NFyNum = NQn + NAn - 1; // Qn · An
    static constexpr size_t NFyDen = NQd + NBn - 1; // Qd · Bn

    damp::array<T, NFyNum> fy_num{};
    damp::array<T, NFyDen> fy_den{};
    damp::array<T, NQn>    fu_num{}; // = Qn
    damp::array<T, NQd>    fu_den{}; // = Qd
    bool                   success{false};

    template<typename U>
    [[nodiscard]] constexpr ClassicalDobResult<NBn, NAn, NQn, NQd, std::remove_const_t<U>> as() const {
        ClassicalDobResult<NBn, NAn, NQn, NQd, std::remove_const_t<U>> out{};
        using O = std::remove_const_t<U>;
        for (size_t i = 0; i < NFyNum; ++i) {
            out.fy_num[i] = static_cast<O>(fy_num[i]);
        }
        for (size_t i = 0; i < NFyDen; ++i) {
            out.fy_den[i] = static_cast<O>(fy_den[i]);
        }
        for (size_t i = 0; i < NQn; ++i) {
            out.fu_num[i] = static_cast<O>(fu_num[i]);
        }
        for (size_t i = 0; i < NQd; ++i) {
            out.fu_den[i] = static_cast<O>(fu_den[i]);
        }
        out.success = success;
        return out;
    }
};

/**
 * @brief Synthesize a classical disturbance observer from a nominal plant and Q-filter.
 *
 * Forms Fy = Q·Pn^-1 = (Qn·An)/(Qd·Bn) and Fu = Q = Qn/Qd. Validates causal
 * realizability: the leading coefficients Bn[0], An[0], Qd[0] must be nonzero
 * (the nominal plant must be causally invertible — no pure input delay).
 *
 * @param Bn,An Nominal plant Pn = Bn/An as z^-1 polynomials
 * @param Qn,Qd Low-pass Q-filter = Qn/Qd as z^-1 polynomials (DC gain ~1)
 */
template<size_t NBn, size_t NAn, size_t NQn, size_t NQd, typename T = double>
[[nodiscard]] constexpr ClassicalDobResult<NBn, NAn, NQn, NQd, T> classical_dob(
    const damp::array<T, NBn>& Bn,
    const damp::array<T, NAn>& An,
    const damp::array<T, NQn>& Qn,
    const damp::array<T, NQd>& Qd
) {
    const T tol = default_tol<T>();
    if (damp::abs(Bn[0]) <= tol || damp::abs(An[0]) <= tol || damp::abs(Qd[0]) <= tol) {
        return ClassicalDobResult<NBn, NAn, NQn, NQd, T>{}; // not causally realizable
    }
    return ClassicalDobResult<NBn, NAn, NQn, NQd, T>{
        .fy_num = detail::poly_conv(Qn, An), // Q·Pn^-1 numerator
        .fy_den = detail::poly_conv(Qd, Bn), // Q·Pn^-1 denominator
        .fu_num = Qn,
        .fu_den = Qd,
        .success = true,
    };
}

/**
 * @ingroup estimators
 * @brief Classical Pn^-1·Q disturbance observer runtime (bolt-on compensator).
 *
 * Drop it around an existing controller: feed it the measurement and your
 * controller's command each tick, and it returns the disturbance-compensated
 * command. The one-sample delay on the applied input breaks the algebraic loop
 * (standard discrete-DOB practice) and is exact at DC.
 *
 * @tparam NBn,NAn,NQn,NQd Plant/Q sizes (must match the design result)
 * @tparam T Scalar type (float or double)
 */
template<size_t NBn, size_t NAn, size_t NQn, size_t NQd, typename T = float>
class ClassicalDOB {
public:
    using result_type = ClassicalDobResult<NBn, NAn, NQn, NQd, T>;

    constexpr ClassicalDOB() = default;

    constexpr explicit ClassicalDOB(const result_type& design)
        : valid_(design.success) {
        fy_.num = design.fy_num;
        fy_.den = design.fy_den;
        fu_.num = design.fu_num;
        fu_.den = design.fu_den;
    }

    /// Estimate the input-referred disturbance from measurement @p y and the
    /// applied input @p u (advances the internal filters one step).
    constexpr T estimate(T y, T u) {
        if (!valid_) {
            return T{0};
        }
        d_hat_ = fy_.step(y) - fu_.step(u);
        return d_hat_;
    }

    /// Bolt-on: return the disturbance-compensated command u = u_command − d_hat,
    /// using the previously applied command in the Q·u path (breaks the loop).
    constexpr T compensate(T u_command, T y) {
        if (!valid_) {
            return u_command;
        }
        d_hat_ = fy_.step(y) - fu_.step(u_prev_);
        const T u = u_command - d_hat_;
        u_prev_ = u;
        return u;
    }

    [[nodiscard]] constexpr T    disturbance() const { return d_hat_; }
    [[nodiscard]] constexpr bool valid() const { return valid_; }

    constexpr void reset() {
        fy_.reset();
        fu_.reset();
        d_hat_ = T{0};
        u_prev_ = T{0};
    }

private:
    detail::IirDF2<ClassicalDobResult<NBn, NAn, NQn, NQd, T>::NFyNum, ClassicalDobResult<NBn, NAn, NQn, NQd, T>::NFyDen, T>
                                fy_{};
    detail::IirDF2<NQn, NQd, T> fu_{};
    T                           d_hat_{T{0}};
    T                           u_prev_{T{0}};
    bool                        valid_{false};
};

} // namespace damp::estimation
