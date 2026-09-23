// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file rls.hpp
 * @brief Recursive Least Squares (RLS) primitives for online identification.
 *
 * Scalar and vector RLS for linear-in-parameters regression
 * @f$ y = \phi^\top \theta + e @f$. Updates fail closed (return false) when
 * the configuration is invalid or the instantaneous information is degenerate
 * (@f$ \lambda + \phi^\top P \phi \approx 0 @f$ — zero / collinear regressor).
 *
 * Persistent excitation (PE) is the caller's responsibility: without rich
 * @f$\phi@f$ the covariance stays large and @c valid() alone is not a
 * confidence claim. Pair with FirstOrderPlantEstimator::confidence or an
 * external PE gate before acting on @c theta.
 *
 * @note Compare with MATLAB®'s recursiveLS / dsp.RLSFilter (forgetting form).
 * @see Ljung, "System Identification," 2nd ed., 1999, ch. 11
 * @see parameter_estimation.hpp for grey-box extraction + AdaptationGate
 */

#include <cstddef>
#include <limits>
#include <type_traits>

#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"

namespace damp::estimation {

/**
 * @brief Common RLS configuration.
 */
template<typename T = double>
struct RlsConfig {

    using scalar_type = scalar_type_t<T>;

    scalar_type lambda{scalar_type{1}};
    scalar_type p0{scalar_type{1}};
    bool        projection_enabled{false};
    T           theta_min{};
    T           theta_max{};

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
 * @brief Scalar RLS runtime state.
 */
template<typename T = float>
struct RlsState {

    T    theta{};
    T    covariance{T{1}};
    T    predicted_output{};
    T    residual{};
    T    gain{};
    bool initialized{false};
};

/**
 * @brief Scalar RLS design payload.
 */
template<typename T = double>
struct RlsResult {

    RlsConfig<T> config{};
    T            theta0{};
    T            covariance0{T{1}};
    bool         success{false};

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        using out_t = std::remove_const_t<U>;
        return RlsResult<out_t>{
            RlsConfig<out_t>{
                static_cast<scalar_type_t<out_t>>(config.lambda),
                static_cast<scalar_type_t<out_t>>(config.p0),
                config.projection_enabled,
                static_cast<out_t>(config.theta_min),
                static_cast<out_t>(config.theta_max),
            },
            static_cast<out_t>(theta0),
            static_cast<out_t>(covariance0),
            success,
        };
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
struct RlsVectorResult {

    RlsConfig<T>      config{};
    ColVec<NP, T>     theta0{};
    Matrix<NP, NP, T> covariance0{Matrix<NP, NP, T>::identity()};
    bool              success{false};

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        using out_t = std::remove_const_t<U>;
        return RlsVectorResult<NP, out_t>{
            RlsConfig<out_t>{
                static_cast<scalar_type_t<out_t>>(config.lambda),
                static_cast<scalar_type_t<out_t>>(config.p0),
                config.projection_enabled,
                static_cast<out_t>(config.theta_min),
                static_cast<out_t>(config.theta_max),
            },
            theta0.template as<out_t>(),
            covariance0.template as<out_t>(),
            success,
        };
    }
};

/**
 * @brief Build scalar RLS design payload.
 */
template<typename T = double>
[[nodiscard]] constexpr RlsResult<T> rls(
    const RlsConfig<T>& config,
    T                   theta0 = T{}
) {
    return RlsResult<T>{
        .config = config,
        .theta0 = theta0,
        .covariance0 = static_cast<T>(config.p0),
        .success = config.valid(),
    };
}

/**
 * @brief Build vector RLS design payload.
 */
template<size_t NP, typename T = double>
[[nodiscard]] constexpr RlsVectorResult<NP, T> rls_vector(
    const RlsConfig<T>&  config,
    const ColVec<NP, T>& theta0 = ColVec<NP, T>{}
) {
    if (!config.valid()) {
        return RlsVectorResult<NP, T>{
            .config = config,
            .theta0 = theta0,
            .success = false,
        };
    }
    return RlsVectorResult<NP, T>{
        .config = config,
        .theta0 = theta0,
        .covariance0 = Matrix<NP, NP, T>::identity() * static_cast<T>(config.p0),
        .success = true,
    };
}

/**
 * @brief Scalar runtime RLS estimator.
 */
template<typename T = float>
class Rls {
public:
    using scalar_type = scalar_type_t<T>;

    constexpr Rls() = default;

    constexpr explicit Rls(const RlsConfig<T>& config)
        : config_(config) {
        state_.covariance = static_cast<T>(config_.p0);
    }

    constexpr explicit Rls(const RlsResult<T>& design)
        : config_(design.config), valid_(design.success) {
        state_.theta = design.theta0;
        state_.covariance = design.covariance0;
    }

    /**
     * @brief Update estimate using one regression sample y = phi * theta + noise.
     *
     * @return false if config invalid or information λ + φᴴ P φ
     *         is degenerate (no PE this sample). State is left unchanged on false.
     */
    [[nodiscard]] constexpr bool update(T phi, T y) {
        if (!config_.valid()) {
            valid_ = false;
            return false;
        }

        // Hermitian form (conj(phi)) so the complex scalar case matches the
        // vector path's phi^H and keeps the denominator real; conj is the identity
        // for real T.
        const T phi_h = damp::conj(phi);
        const T denom = static_cast<T>(config_.lambda) + phi_h * state_.covariance * phi;
        // Degenerate information: zero regressor with tiny P, or λ≈0 misconfig.
        // Do not touch theta/P — honest "no update" rather than a NaN gain.
        if (damp::abs(denom) <= default_tol<T>()) {
            return false;
        }

        state_.gain = (state_.covariance * phi) / denom;
        state_.predicted_output = phi_h * state_.theta;
        state_.residual = y - state_.predicted_output;
        state_.theta = state_.theta + state_.gain * state_.residual;

        const T one_over_lambda = T{1} / static_cast<T>(config_.lambda);
        state_.covariance = (state_.covariance - state_.gain * phi_h * state_.covariance) * one_over_lambda;

        if constexpr (!is_complex_v<T>) {
            if (config_.projection_enabled) {
                if (state_.theta < config_.theta_min) {
                    state_.theta = config_.theta_min;
                }
                if (state_.theta > config_.theta_max) {
                    state_.theta = config_.theta_max;
                }
            }

            // Floor above zero. A covariance of exactly 0 freezes the gain at 0.
            const T cov_floor = std::numeric_limits<T>::epsilon() * (damp::abs(static_cast<T>(config_.p0)) + T{1});
            if (state_.covariance < cov_floor) {
                state_.covariance = cov_floor;
            }
        }

        state_.initialized = true;
        valid_ = true;
        return true;
    }

    [[nodiscard]] constexpr T predict(T phi) const {
        return damp::conj(phi) * state_.theta;
    }

    constexpr void reset(T theta0 = T{}) {
        state_ = RlsState<T>{};
        state_.theta = theta0;
        state_.covariance = static_cast<T>(config_.p0);
        valid_ = true;
    }

    [[nodiscard]] constexpr const auto& config() const { return config_; }
    [[nodiscard]] constexpr const auto& state() const { return state_; }
    [[nodiscard]] constexpr bool        valid() const { return valid_; }

private:
    RlsConfig<T> config_{};
    RlsState<T>  state_{};
    bool         valid_{true};
};

/**
 * @brief Vector runtime RLS estimator (NP parameters).
 */
template<size_t NP, typename T = float>
class RlsVector {
public:
    constexpr RlsVector() = default;

    constexpr explicit RlsVector(const RlsConfig<T>& config)
        : config_(config) {
        state_.covariance = Matrix<NP, NP, T>::identity() * static_cast<T>(config_.p0);
    }

    constexpr explicit RlsVector(const RlsVectorResult<NP, T>& design)
        : config_(design.config), valid_(design.success) {
        state_.theta = design.theta0;
        state_.covariance = design.covariance0;
    }

    /**
     * @brief Update estimate using one sample y = phi^H * theta + noise.
     *
     * @return false if config invalid or information
     *         λ + φᵀ P φ is degenerate. On false the prior
     *         (theta, P) is left unchanged — not a PE claim of success.
     */
    [[nodiscard]] constexpr bool update(const ColVec<NP, T>& phi, T y) {
        if (!config_.valid()) {
            valid_ = false;
            return false;
        }

        const ColVec<NP, T> p_phi = state_.covariance * phi;
        const T             denom = static_cast<T>(config_.lambda) + dot(phi, p_phi);
        // No instantaneous PE (e.g. phi = 0): skip without poisoning P or theta.
        if (damp::abs(denom) <= default_tol<T>()) {
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
                for (size_t i = 0; i < NP; ++i) {
                    if (state_.theta[i] < config_.theta_min) {
                        state_.theta[i] = config_.theta_min;
                    }
                    if (state_.theta[i] > config_.theta_max) {
                        state_.theta[i] = config_.theta_max;
                    }
                }
            }
        }

        // Keep covariance numerically symmetric/Hermitian after finite-precision updates.
        state_.covariance = (state_.covariance + state_.covariance.conjugate_transpose()) * static_cast<T>(0.5);
        state_.initialized = true;
        valid_ = true;
        return true;
    }

    [[nodiscard]] constexpr T predict(const ColVec<NP, T>& phi) const {
        return dot(phi, state_.theta);
    }

    constexpr void reset(const ColVec<NP, T>& theta0 = ColVec<NP, T>{}) {
        state_ = RlsVectorState<NP, T>{};
        state_.theta = theta0;
        state_.covariance = Matrix<NP, NP, T>::identity() * static_cast<T>(config_.p0);
        valid_ = true;
    }

    [[nodiscard]] constexpr const auto& config() const { return config_; }
    [[nodiscard]] constexpr const auto& state() const { return state_; }
    [[nodiscard]] constexpr bool        valid() const { return valid_; }

private:
    RlsConfig<T>          config_{};
    RlsVectorState<NP, T> state_{};
    bool                  valid_{true};
};

} // namespace damp::estimation
