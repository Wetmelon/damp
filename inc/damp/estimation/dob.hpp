// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file dob.hpp
 * @brief Disturbance-observer estimation primitives and lightweight runtime.
 */


#include "damp/backend.hpp" // damp::array
#include "damp/math/math.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "damp/systems/transfer_function.hpp"

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

    /**
     * @brief Constant-innovation steady-state gain d_hat_ss/innovation = gain/leak.
     * Sentinel 0 means "no finite steady state": leak == 0 is a pure integrator
     * whose estimate grows without bound under constant innovation.
     */
    T    steady_state_gain{};
    bool success{false};

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        using out_t = damp::remove_const_t<U>;
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
    const auto scale = max_magnitude / magnitude;
    return value * scale;
}

/**
 * @brief True if p[0] is zero to rounding relative to the polynomial's largest
 * coefficient (scale-invariant: a plant gain in other units is not "zero").
 */
template<size_t N, typename T>
[[nodiscard]] constexpr bool leading_coefficient_negligible(const damp::array<T, N>& p) {
    scalar_type_t<T> scale{0};
    for (size_t i = 0; i < N; ++i) {
        scale = damp::max(scale, damp::abs(p[i]));
    }
    return !(damp::abs(p[0]) > static_cast<scalar_type_t<T>>(N) * eps<T>() * scale);
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

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr DOB(const DOBResult<U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : DOB(result.template as<T>()) {}

    constexpr DOB(const DOBResult<T>& design) // NOLINT(google-explicit-constructor): assign-from-design
        : config_(design.config), valid_(design.success) {
        if (!design.success && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    /**
     * @brief Update disturbance estimate from predicted and measured outputs.
     * @return true when update succeeded.
     */
    [[nodiscard]] constexpr bool update(T y_predicted, T y_measured) {
        if (!valid_) {
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
        return true;
    }

    /**
     * @brief Disturbance-compensated command (u = u_nominal - d_hat).
     */
    [[nodiscard]] constexpr T compensate(const T& u_nominal) const {
        return u_nominal - state_.disturbance_hat;
    }

    constexpr void reset() { state_ = DOBState<T>{}; }

    [[nodiscard]] constexpr const auto& config() const { return config_; }
    [[nodiscard]] constexpr const auto& state() const { return state_; }
    [[nodiscard]] constexpr bool        valid() const { return valid_; }

private:
    DOBConfig<T> config_{};
    DOBState<T>  state_{};
    bool         valid_{true};
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<typename T>
DOB(const DOBResult<T>&) -> DOB<T>;

// ===========================================================================
// Classical Pn^-1 * Q disturbance observer (Ohnishi DOB) — bolt-on for an
// existing controller. Estimates the input-referred lumped disturbance from a
// nominal plant model and a low-pass Q-filter, and subtracts it from the
// command. Complements the scalar innovation-based observer above (richer plant
// model) and the extended-state observer (`design::eso` / `ESO`; ADRC is that
// observer plus PD). This is a feed-around bolt-on, not a full controller.
//
//   plant:        y = P·(u + d)          (d = input-referred disturbance)
//   estimate:     d_hat = Q·(Pn^-1·y − u) ≈ Q·d   (exact at DC when P = Pn)
//   compensate:   u = u_command − d_hat
//
// All polynomials are discrete z^-1 digital filters (c0 + c1·z^-1 + …), the
// natural form for an embedded loop. Pn is the sampled plant with its
// one-sample ZOH delay removed (y[k+1] = a·y[k] + b·u[k] gives Bn = {b},
// An = {1, −a}); the runtime supplies that delay by feeding the previous
// command to Q. Pn⁻¹ must be causal and stable: Bn[0] ≠ 0 and every zero of Bn
// strictly inside the unit circle.
//
// @see "Disturbance Observer-Based Control" (S. Li et al., CRC Press, 2016)
// @see "Robust Motion Control by Disturbance Observer" (Ohnishi et al., 1996)
// ===========================================================================

namespace detail {

/** @brief Polynomial (z^-1) convolution = digital-filter series multiplication. */
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

/**
 * @brief The z⁻¹ filter (c₀ + c₁z⁻¹ + …)/(d₀ + d₁z⁻¹ + …) as a TransferFunction in z
 *
 * Multiplying both by z^(N−1), N the longer length, keeps each array's order: highest power of z
 * first, as TransferFunction stores it. The shorter one is zero-padded at the end (its missing
 * z⁻ᵏ terms are the lowest powers of z).
 */
template<size_t NN, size_t ND, typename T>
[[nodiscard]] constexpr TransferFunction<(NN > ND ? NN : ND), (NN > ND ? NN : ND), T>
zinv_to_tf(const damp::array<T, NN>& num, const damp::array<T, ND>& den) {
    TransferFunction<(NN > ND ? NN : ND), (NN > ND ? NN : ND), T> tf{};
    for (size_t i = 0; i < NN; ++i) {
        tf.num[i] = num[i];
    }
    for (size_t i = 0; i < ND; ++i) {
        tf.den[i] = den[i];
    }
    return tf;
}

/** @brief Direct-Form-II IIR for a digital filter H(z^-1) = num / den (den[0] ≠ 0). */
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

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr ClassicalDobResult<NBn, NAn, NQn, NQd, damp::remove_const_t<U>> as() const {
        ClassicalDobResult<NBn, NAn, NQn, NQd, damp::remove_const_t<U>> out{};
        using O = damp::remove_const_t<U>;
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

    /**
     * @brief Q-filter Fu = Q as a discrete TF in z (the IIR's coefficients, the shorter array zero-padded)
     */
    [[nodiscard]] constexpr auto q_tf() const { return detail::zinv_to_tf(fu_num, fu_den); }

    /**
     * @brief Measurement filter Fy = Q·Pₙ⁻¹ as a discrete TF in z (as @ref q_tf)
     */
    [[nodiscard]] constexpr auto fy_tf() const { return detail::zinv_to_tf(fy_num, fy_den); }
};

/**
 * @brief Synthesize a classical disturbance observer from a nominal plant and Q-filter.
 *
 * Forms Fy = Q·Pn^-1 = (Qn·An)/(Qd·Bn) and Fu = Q = Qn/Qd. Pn is the sampled
 * plant with its one-sample ZOH delay removed: y[k+1] = a·y[k] + b·u[k] is
 * Bn = {b}, An = {1, −a}.
 *
 * Refused (success = false) unless both filters are causal and stable: Bn[0],
 * An[0] and Qd[0] nonzero, and Bn and Qd Schur-stable (every zero of Bn becomes
 * a pole of Fy). A ZOH plant of relative degree two has a sampling zero near
 * z = −1: on it for a pure double integrator (refused), just inside for
 * 1/(s(s+1)) (−0.997 at Ts = 0.01, passes). Relative degree three or more puts
 * one outside (refused). A zero z₀ just inside the circle passes, and Fy then
 * amplifies noise near Nyquist by about 1/(1−|z₀|). Such a plant needs a stable
 * approximate inverse in Bn, e.g. the single coefficient {Σbᵢ} (same DC gain).
 *
 * @param Bn,An Nominal plant Pn = Bn/An as z^-1 polynomials (delay removed)
 * @param Qn,Qd Low-pass Q-filter = Qn/Qd as z^-1 polynomials (DC gain ~1)
 */
template<size_t NBn, size_t NAn, size_t NQn, size_t NQd, typename T = double>
[[nodiscard]] constexpr ClassicalDobResult<NBn, NAn, NQn, NQd, T> classical_dob(
    const damp::array<T, NBn>& Bn,
    const damp::array<T, NAn>& An,
    const damp::array<T, NQn>& Qn,
    const damp::array<T, NQd>& Qd
) {
    if (detail::leading_coefficient_negligible(Bn) || detail::leading_coefficient_negligible(An)
        || detail::leading_coefficient_negligible(Qd)) {
        return ClassicalDobResult<NBn, NAn, NQn, NQd, T>{}; // not causally realizable
    }
    if (!tf_detail::schur_stable(Bn) || !tf_detail::schur_stable(Qd)) {
        return ClassicalDobResult<NBn, NAn, NQn, NQd, T>{}; // Fy or Fu unstable
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
 * command to apply. The Q·u path uses the command returned on the previous
 * tick, which is the one-sample plant delay classical_dob removes from Pn.
 *
 * @tparam NBn,NAn,NQn,NQd Plant/Q sizes (must match the design result)
 * @tparam T Scalar type (float or double)
 */
template<size_t NBn, size_t NAn, size_t NQn, size_t NQd, typename T = float>
class ClassicalDOB {
public:
    using result_type = ClassicalDobResult<NBn, NAn, NQn, NQd, T>;

    constexpr ClassicalDOB() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr ClassicalDOB(const ClassicalDobResult<NBn, NAn, NQn, NQd, U>& design) // NOLINT(google-explicit-constructor): assign-from-design
        : ClassicalDOB(design.template as<T>()) {}

    constexpr ClassicalDOB(const result_type& design) // NOLINT(google-explicit-constructor): assign-from-design
        : fy_{.num = design.fy_num, .den = design.fy_den}, fu_{.num = design.fu_num, .den = design.fu_den}, valid_(design.success) {
        if (!design.success && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    /**
     * @brief One tick: estimate d̂ from measurement @p y and return the command to apply,
     * u = u_command − d̂.
     *
     * The Q·u path takes the command returned on the previous tick, so apply the returned
     * value. d̂ is read back with @ref disturbance.
     *
     * @param y         Measurement this tick
     * @param u_command Controller command this tick
     * @return u_command − d̂ (u_command unchanged when the design failed)
     */
    [[nodiscard]] constexpr T step(T y, T u_command) {
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

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<size_t NBn, size_t NAn, size_t NQn, size_t NQd, typename T>
ClassicalDOB(const ClassicalDobResult<NBn, NAn, NQn, NQd, T>&) -> ClassicalDOB<NBn, NAn, NQn, NQd, T>;

} // namespace damp::estimation
