// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file rls.hpp
 * @brief Recursive Least Squares (RLS) primitives for online identification.
 *
 * RLS for linear-in-parameters regression @f$ y = \phi^\top \theta + e @f$
 * with NP parameters; NP = 1 is the scalar case (a ColVec<1, T> regressor
 * converts from T). Updates fail closed (return false) when the design failed
 * or @f$ \lambda + \phi^\top P \phi \le 0 @f$ (P has lost positive
 * semi-definiteness). A zero regressor is a valid sample that carries no
 * information: θ is unchanged.
 *
 * Forgetting (λ < 1) divides P by λ every sample, so without persistent
 * excitation P grows without bound (estimator windup). The covariance is
 * therefore capped at its prior, tr P ≤ NP·p₀, by uniform rescaling (the
 * bounded-trace algorithm).
 *
 * Persistent excitation (PE) is the caller's responsibility: without rich
 * @f$\phi@f$ the covariance stays large and @c valid() alone is not a
 * confidence claim. Pair with FirstOrderPlantEstimator::confidence or an
 * external PE gate before acting on @c theta.
 *
 * @note Compare with MATLAB®'s recursiveLS / dsp.RLSFilter (forgetting form).
 * @see "System Identification: Theory for the User" (Ljung, 1999), §11.2
 * @see parameter_estimation.hpp for grey-box extraction + AdaptationGate
 */


#include "damp/backend.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "damp/toolbox/bounds.hpp"

namespace damp::estimation {

namespace detail {

/** @brief Real part of a real or complex scalar. */
template<typename T>
[[nodiscard]] constexpr scalar_type_t<T> real_part(const T& x) {
    if constexpr (is_complex_v<T>) {
        return x.real();
    } else {
        return x;
    }
}

} // namespace detail

/**
 * @brief Common RLS configuration.
 * @tparam T  Scalar type (real or complex)
 * @tparam NP Number of parameters
 */
template<typename T = double, size_t NP = 1>
struct RlsConfig {

    using scalar_type = scalar_type_t<T>;

    scalar_type             lambda{scalar_type{1}};
    scalar_type             p0{scalar_type{1}};
    bool                    projection_enabled{false};
    Bounds<NP, scalar_type> theta_lim{}; ///< Per-parameter projection box (real T only; default unbounded)

    [[nodiscard]] constexpr bool valid() const {
        if (lambda <= scalar_type{0}) {
            return false;
        }
        if (lambda > scalar_type{1}) {
            return false;
        }
        if (p0 <= scalar_type{0}) {
            return false;
        }
        return true;
    }
};

/**
 * @brief Vector RLS runtime state for N parameters.
 */
template<size_t NP, typename T = float>
struct RlsVectorState {

    ColVec<NP, T>     theta{};
    Matrix<NP, NP, T> covariance{Matrix<NP, NP, T>::identity()};
    T                 predicted_output{};
    T                 residual{};
    ColVec<NP, T>     gain{};
    bool              initialized{false};
};

/**
 * @brief Vector RLS design payload for N parameters.
 */
template<size_t NP, typename T = double>
struct [[nodiscard]] RlsVectorResult {

    RlsConfig<T, NP> config{}; ///< The runtime starts from P₀ = p₀·I
    ColVec<NP, T>    theta0{};
    bool             success{false};

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        using out_t = damp::remove_const_t<U>;
        return RlsVectorResult<NP, out_t>{
            RlsConfig<out_t, NP>{
                static_cast<scalar_type_t<out_t>>(config.lambda),
                static_cast<scalar_type_t<out_t>>(config.p0),
                config.projection_enabled,
                config.theta_lim.template as<scalar_type_t<out_t>>(),
            },
            theta0.template as<out_t>(),
            success,
        };
    }
};

/**
 * @brief Build vector RLS design payload.
 */
template<size_t NP, typename T = double>
[[nodiscard]] constexpr RlsVectorResult<NP, T> rls_vector(
    const RlsConfig<T, NP>& config,
    const ColVec<NP, T>&    theta0 = ColVec<NP, T>{}
) {
    return RlsVectorResult<NP, T>{
        .config = config,
        .theta0 = theta0,
        .success = config.valid(),
    };
}

/**
 * @brief Vector runtime RLS estimator (NP parameters).
 */
template<size_t NP, typename T = float>
class RlsVector {
public:
    using scalar_type = scalar_type_t<T>;

    constexpr RlsVector() = default;

    /** @brief From a design result in another scalar: converts it with .as\<T\>() and builds from that. */
    template<typename U>
        requires(!damp::is_same_v<U, T>)
    constexpr RlsVector(const RlsVectorResult<NP, U>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : RlsVector(result.template as<T>()) {}

    constexpr RlsVector(const RlsVectorResult<NP, T>& design) // NOLINT(google-explicit-constructor): assign-from-design
        : config_(design.config), state_{.theta = design.theta0, .covariance = prior_covariance(design.config)}, valid_(design.success) {
        if (!design.success && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    /**
     * @brief Update estimate using one sample y = phi^H * theta + noise.
     *
     * @return false if the design failed or information
     *         λ + φᵀ P φ is degenerate. On false the prior
     *         (theta, P) is left unchanged — not a PE claim of success.
     */
    [[nodiscard]] constexpr bool update(const ColVec<NP, T>& phi, T y) {
        if (!valid_) {
            return false;
        }

        const ColVec<NP, T> p_phi = state_.covariance * phi;
        const T             denom = static_cast<T>(config_.lambda) + dot(phi, p_phi);
        // denom ≥ λ > 0 while P ≥ 0; otherwise P is corrupt: skip without
        // poisoning P or theta.
        if (!(detail::real_part(denom) > scalar_type{0})) {
            return false;
        }

        state_.gain = p_phi / denom;
        state_.predicted_output = dot(phi, state_.theta);
        state_.residual = y - state_.predicted_output;
        state_.theta += state_.gain * state_.residual;

        const auto phi_h = phi.conjugate_transpose();
        const auto one_over_lambda = T{1} / static_cast<T>(config_.lambda);
        state_.covariance = (state_.covariance - (state_.gain * phi_h * state_.covariance)) * one_over_lambda;

        if constexpr (!is_complex_v<T>) {
            if (config_.projection_enabled) {
                state_.theta = config_.theta_lim.saturate(state_.theta);
            }
        }

        // Keep covariance numerically symmetric/Hermitian after finite-precision updates.
        state_.covariance = (state_.covariance + state_.covariance.conjugate_transpose()) * static_cast<T>(0.5);

        // Anti-windup (bounded trace): forgetting without excitation must not
        // let P grow past its prior, tr P ≤ N·p₀.
        const scalar_type tr = detail::real_part(state_.covariance.trace());
        const scalar_type tr_max = static_cast<scalar_type>(NP) * config_.p0;
        if (tr > tr_max) {
            state_.covariance = state_.covariance * static_cast<T>(tr_max / tr);
        }
        state_.initialized = true;
        return true;
    }

    [[nodiscard]] constexpr T predict(const ColVec<NP, T>& phi) const {
        return dot(phi, state_.theta);
    }

    constexpr void reset(const ColVec<NP, T>& theta0 = ColVec<NP, T>{}) {
        state_ = RlsVectorState<NP, T>{.theta = theta0, .covariance = prior_covariance(config_)};
    }

    [[nodiscard]] constexpr const auto& config() const { return config_; }
    [[nodiscard]] constexpr const auto& state() const { return state_; }
    [[nodiscard]] constexpr bool        valid() const { return valid_; }

private:
    [[nodiscard]] static constexpr Matrix<NP, NP, T> prior_covariance(const RlsConfig<T, NP>& config) {
        return Matrix<NP, NP, T>::identity() * static_cast<T>(config.p0);
    }

    RlsConfig<T, NP>      config_{};
    RlsVectorState<NP, T> state_{};
    bool                  valid_{true};
};

/** @brief Deduce the runtime from its design, keeping the result's scalar. */
template<size_t NP, typename T>
RlsVector(const RlsVectorResult<NP, T>&) -> RlsVector<NP, T>;

} // namespace damp::estimation
