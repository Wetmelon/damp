// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file damp/estimation/state_variable_filter.hpp
 * @brief Regression front end: filtered signals and their derivatives without differentiating.
 *
 * Physical parameters usually sit in an equation with derivatives in it, e.g. a belt
 * @f$ F = k\,\delta + c\,\dot\delta @f$ or a mass @f$ F = m\,\ddot x + b\,\dot x @f$. It is
 * linear in the parameters, so recursive least squares (RlsVector) can fit it online,
 * but the derivatives of quantized encoder signals are mostly noise.
 *
 * The state-variable filter passes every signal of the regression through the same
 * stable filter
 * @f[
 *   F(s) = \frac{\omega^3}{(s + \omega)^3},
 * @f]
 * whose three internal states are the filtered signal and its first and second
 * derivatives. No numerical differentiation happens anywhere, and because the second
 * derivative is a state (relative degree 3), none of the outputs feeds the raw sample
 * through. Filtering both sides of @f$ y = \phi^\top \theta @f$ with the same linear
 * filter leaves θ unchanged, so the regression runs on the filtered quantities as is:
 *
 * @code
 * // F = k δ + c δ̇: filter δ and F, regress F_f on [δ_f, δ̇_f].
 * constexpr auto svf = estimation::state_variable_filter(2.0 * pi * 20.0, Ts); // 20 Hz corner
 * StateVariableFilter<2, float> f{svf};
 * f.step({delta_um, force});                       // δ in µm: see the note on scaling
 * (void)rls.update(ColVec<2, float>{f[0].x, f[0].dx}, f[1].x);
 * @endcode
 *
 * Choose ω above the dynamics being identified and below where the sensor noise lives.
 * The filter is discretized by zero-order hold.
 *
 * @note Scale the regression, and size p₀ to it. RLS starts from a prior θ₀ with covariance
 *       p₀·I. Too small and the data never outweigh it (they need Σφ² ≫ 1/p₀: a 10⁻⁵ m belt
 *       stretch in metres is pulled to a quarter of its true stiffness at p₀ = 10⁶). Too large
 *       and φᵀPφ dwarfs λ, so the covariance update cancels away in float (a 10⁵ µm/s²
 *       acceleration with p₀ = 10⁶ loses the damping term). Pick units that put the regressors
 *       near 1–10⁵, and p₀ a few orders above the expected θ².
 *
 * `constexpr`, allocation-free, embeddable.
 *
 * @see "Parameter estimation for continuous-time models — a survey" (Young, 1981),
 *      Automatica 17(1), §3 (the state-variable filter method)
 * @see "System Identification: Theory for the User" (Ljung, 1999), §11.2 (RLS)
 * @see rls.hpp (RlsVector), parameter_estimation.hpp (ParameterDriftMonitor, AdaptationGate)
 */


#include "damp/backend.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

namespace damp::estimation {

/** @brief One signal after the state-variable filter: its filtered value and first two derivatives. */
template<typename T = float>
struct Filtered {
    T x{};   ///< filtered signal
    T dx{};  ///< its first derivative
    T ddx{}; ///< its second derivative
};

/**
 * @brief Discretized state-variable filter (shared by every channel).
 * @tparam T Scalar type.
 */
template<typename T = double>
struct [[nodiscard]] StateVariableFilterResult {
    Matrix<3, 3, T> Ad{};    ///< state update over one sample, states [x_f, ẋ_f, ẍ_f]
    ColVec<3, T>    Bd{};    ///< input coupling over one sample
    T               omega{}; ///< [rad/s] filter corner
    bool            success{false};

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr StateVariableFilterResult<U> as() const {
        return {Ad.template as<U>(), Bd.template as<U>(), static_cast<U>(omega), success};
    }
};

/**
 * @brief Design the filter F(s) = ω³/(s + ω)³ sampled at @p Ts.
 * @param omega [rad/s] corner (triple pole)
 * @param Ts    [s] sample period
 */
template<typename T = double>
[[nodiscard]] constexpr StateVariableFilterResult<T> state_variable_filter(T omega, T Ts) {
    StateVariableFilterResult<T> r{.omega = omega};
    if (!(omega > T{0}) || !(Ts > T{0})) {
        return r;
    }
    // x⃛_f = ω³(u − x_f) − 3ω² ẋ_f − 3ω ẍ_f in companion form.
    StateSpace<3, 1, 1, T> cont{};
    cont.A(0, 1) = T{1};
    cont.A(1, 2) = T{1};
    cont.A(2, 0) = -(omega * omega * omega);
    cont.A(2, 1) = -T{3} * omega * omega;
    cont.A(2, 2) = -T{3} * omega;
    cont.B(2, 0) = omega * omega * omega;
    cont.C(0, 0) = T{1};
    const auto disc = discretize(cont, Ts, DiscretizationMethod::ZOH);
    if (!disc) {
        return r;
    }
    r.Ad = disc->A;
    r.Bd = ColVec<3, T>{disc->B(0, 0), disc->B(1, 0), disc->B(2, 0)};
    r.success = true;
    return r;
}

/**
 * @brief State-variable filter over @p NS signals: filtered values and derivatives, no differentiation.
 * @tparam NS Number of signals (every signal of the regression, both sides).
 * @tparam T  Scalar type.
 */
template<size_t NS, typename T = float>
class StateVariableFilter {
public:
    constexpr StateVariableFilter() = default;

    template<typename U>
    constexpr StateVariableFilter(const StateVariableFilterResult<U>& design) // NOLINT(google-explicit-constructor): assign-from-design
        : Ad_(design.Ad.template as<T>()), Bd_(design.Bd.template as<T>()), valid_(design.success) {
        if (!design.success && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    /** @brief Advance one sample with the raw signals @p u. */
    constexpr void step(const ColVec<NS, T>& u) {
        for (size_t i = 0; i < NS; ++i) {
            const T s0 = s_[i].x;
            const T s1 = s_[i].dx;
            const T s2 = s_[i].ddx;
            s_[i].x = (Ad_(0, 0) * s0) + (Ad_(0, 1) * s1) + (Ad_(0, 2) * s2) + (Bd_[0] * u[i]);
            s_[i].dx = (Ad_(1, 0) * s0) + (Ad_(1, 1) * s1) + (Ad_(1, 2) * s2) + (Bd_[1] * u[i]);
            s_[i].ddx = (Ad_(2, 0) * s0) + (Ad_(2, 1) * s1) + (Ad_(2, 2) * s2) + (Bd_[2] * u[i]);
        }
    }

    /** @brief Signal @p i after the filter. */
    [[nodiscard]] constexpr const Filtered<T>& operator[](size_t i) const { return s_[i]; }

    [[nodiscard]] constexpr bool valid() const { return valid_; } ///< false if the design failed

    /** @brief Restart every channel at rest at @p u0 (avoids a start-up transient in the regression). */
    constexpr void reset(const ColVec<NS, T>& u0 = ColVec<NS, T>{}) {
        for (size_t i = 0; i < NS; ++i) {
            s_[i] = Filtered<T>{.x = u0[i]};
        }
    }

private:
    Matrix<3, 3, T> Ad_{};
    ColVec<3, T>    Bd_{};
    bool            valid_{false};
    Filtered<T>     s_[NS]{};
};

} // namespace damp::estimation
