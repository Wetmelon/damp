// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file ekf.hpp
 * @brief Extended Kalman filter on a continuous nonlinear model
 *
 * The model is the continuous Dynamics ẋ = f(x, u) of concepts.hpp, a struct or a captureless
 * lambda: the same object design::linearize and the sim::simulate* plants take. Each reading
 * comes with its Measurement y = h(x, u) and its R. The filter integrates f over its sample
 * time with one RK4 step and linearizes f and h with the jacobian(x, u) the model provides, or
 * by central differences when it has none.
 *
 * @see "Optimal State Estimation" (Simon, 2006), §13.2.3
 */

#include "damp/backend.hpp"
#include "damp/concepts.hpp"
#include "damp/estimation/sqrt_covariance.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/dynamics.hpp"

namespace damp {

/**
 * @brief Configuration of an ExtendedKalmanFilter: its sample time, the process noise and the
 *        initial estimate
 *
 * An aggregate; name the fields you set. The model is the filter's first constructor argument:
 * `ExtendedKalmanFilter ekf{f, {.Ts = 0.01F, .Q = Q, .P0 = P0}}`.
 *
 * @tparam Model Dynamics ẋ = f(x, u), a struct or a lambda whose call operator gives NX and NU
 * @tparam T     Scalar type (default: the model's)
 */
template<typename Model, typename T = detail::model_shape<Model>::scalar>
    requires detail::shaped_model<Model>
struct ExtendedKalmanFilterConfig {
    static constexpr size_t NX = detail::model_shape<Model>::states; ///< States

    T                 Ts{};                              ///< Sample time [s] predict integrates over, finite and positive
    Matrix<NX, NX, T> Q{Matrix<NX, NX, T>::identity()};  ///< Process-noise covariance per step, on the full state
    ColVec<NX, T>     x0{};                              ///< Initial state estimate
    Matrix<NX, NX, T> P0{Matrix<NX, NX, T>::identity()}; ///< Initial estimate covariance
};

/**
 * @brief Extended Kalman filter: continuous nonlinear dynamics, sampled measurements
 *
 * Estimates the state of
 *
 *     ẋ      = f(x, u)                       (Dynamics)
 *     x[k+1] = φ(x[k], u[k]) + w[k],   w ~ N(0, Q)
 *     y[k]   = h(x[k], u[k]) + v[k],   v ~ N(0, R)   (Measurement, per reading)
 *
 * where φ is one classical RK4 step of f over the sample time Ts with u held (zero-order
 * hold). The filter holds its model, Ts and Q; each tick is predict(u) with the input applied
 * now and update(h, y, R) per reading, in any order and at any rate per sensor, each reading of
 * any size. The input h sees is the one the last predict held, the input the plant holds when
 * the reading is taken.
 *
 * Transition Jacobian. The covariance propagates with the RK4 stability function of
 * @f$ M = A\,T_s @f$, where @f$ A = \partial f/\partial x @f$ at @f$ (\hat x[k|k], u[k]) @f$:
 * @f[
 *   \Phi = I + M + \tfrac{1}{2}M^2 + \tfrac{1}{6}M^3 + \tfrac{1}{24}M^4 ,
 * @f]
 * the exact Jacobian ∂φ/∂x of the RK4 step with A held across it. On a linear model the filter is
 * then exactly the Kalman filter of the RK4-discretized plant: mean and covariance propagate
 * with the same Φ. On a nonlinear model it leaves out the change of A across the step, a term
 * of order @f$ T_s^2 @f$; the truncated series @f$ I + M + M^2/2 @f$ leaves out the same term,
 * and the stage-by-stage RK4 Jacobian keeps it at four Jacobian evaluations per tick instead
 * of one.
 *
 * Noise. Q is the covariance of the state increment over one step, on the full state (G = I).
 * A white-noise density @f$ Q_c @f$ on ẋ gives @f$ Q \approx Q_c T_s @f$, or exactly
 * design::discretize_process_noise(A, I, Q_c, Ts). Each reading's R is the covariance of an
 * additive measurement noise, positive semidefinite: a row with no noise is a perfect
 * measurement, taken at the rounding of its own prior variance.
 *
 * Square-root form. The filter carries the lower-triangular factor S of the covariance,
 * P = S Sᵀ, as KalmanFilter does (estimation/sqrt_covariance.hpp): predict triangularizes
 * [Φ S, Q½] by Householder reflections, and an update whitens the reading with the factor of
 * its R and applies Carlson's triangular update one whitened row at a time. A prior
 * much wider than a reading's noise, which a covariance stored as P cannot hold in float, costs
 * the factor nothing. covariance() forms P = S Sᵀ on request; covariance_factor() is S. P₀ and
 * Q are factored once, at construction and in reset(), and may be singular; one that is not
 * positive semidefinite or not finite is refused under constant evaluation, and at run time P₀
 * becomes the identity and such a Q adds no process noise.
 *
 * Central differences, for a model with no jacobian(x, u), step state i by
 * hᵢ = ∛eps·max(|x̂ᵢ|, 1), as linearize does. A state much smaller than one unit (a gap of a few
 * millimetres in metres) is stepped by more than its own size, the more so in float
 * (∛eps ≈ 4.9e-3). State it in units where it is of order one, or give the model an analytic
 * jacobian(x, u). The step does not follow P: a large P₀ would sample f far from x̂, and a
 * small P would put the step below the roundoff of f.
 *
 * Cost per tick, NX states: predict evaluates f 4 times (RK4) plus jacobian(x, u) once, or f
 * 2·NX more times by central differences; an update evaluates h once plus jacobian(x, u) once,
 * or h 2·NX more times, and factors R once (Cholesky). The counts are fixed, so the worst-case
 * time is bounded.
 *
 * @code
 * // Damped pendulum, θ from the upright: θ̈ = (g/l)·sin θ − c·θ̇ + u
 * constexpr auto pendulum = [](const ColVec<2, float>& x, const ColVec<1, float>& u) {
 *     return ColVec<2, float>{x[1], (9.81F * damp::sin(x[0])) - (0.5F * x[1]) + u[0]};
 * };
 * constexpr auto angle = [](const ColVec<2, float>& x, const ColVec<1, float>&) { return ColVec<1, float>{x[0]}; };
 * static constinit ExtendedKalmanFilter ekf{pendulum, {.Ts = 0.01F, .Q = Q, .P0 = P0}};
 * // each tick:
 * (void)ekf.update(angle, y, 1e-6F); // y[k], with the input the plant holds
 * u = law.step(ekf.state());
 * ekf.predict(u);                    // to k + 1 with the input applied now
 * @endcode
 *
 * A struct member names a model struct: `ExtendedKalmanFilter<Pendulum> ekf_{Pendulum{}, {.Ts = 0.01F}};`.
 *
 * @see "Optimal State Estimation" (Simon, 2006), §13.2.3 (discrete-time EKF)
 * @see "Solving Ordinary Differential Equations II" (Hairer & Wanner, 1996), §IV.2 (RK4 stability function)
 * @see "Fast triangular formulation of the square root filter" (Carlson, 1973), AIAA Journal 11(9)
 *
 * @tparam Model Dynamics ẋ = f(x, u), a struct or a lambda whose call operator gives NX and NU
 * @tparam T     Scalar type (default: the model's)
 */
template<typename Model, typename T = detail::model_shape<Model>::scalar>
    requires detail::shaped_model<Model> && Dynamics<Model, detail::model_shape<Model>::states, detail::model_shape<Model>::inputs, T>
struct ExtendedKalmanFilter {
    static constexpr size_t NX = detail::model_shape<Model>::states; ///< States
    static constexpr size_t NU = detail::model_shape<Model>::inputs; ///< Inputs

    constexpr ExtendedKalmanFilter() = default;

    /**
     * @brief Filter from its process model and its configuration
     *
     * A sample time that is not finite and positive is refused under constant evaluation; at run
     * time it becomes 0 and predict holds the state. A P0 or a Q that is not positive semidefinite
     * is refused under constant evaluation too; at run time P0 becomes the identity and Q zero.
     */
    constexpr ExtendedKalmanFilter(const Model& model, const ExtendedKalmanFilterConfig<Model, T>& config)
        : model_(model),
          x(config.x0),
          S(detail::covariance_factor_or_identity(config.P0)),
          Q(config.Q),
          Gamma(detail::process_noise_factor(Matrix<NX, NX, T>::identity(), config.Q)),
          Ts(detail::filter_sample_time(config.Ts)) {}

    /**
     * @brief Time update over one sample: x ← RK4 step of f, P ← Φ·P·Φᵀ + Q on its factor
     *
     *     A        = ∂f/∂x at (x[k|k], u)
     *     Φ        = I + M + M²/2 + M³/6 + M⁴/24,  M = A·Ts
     *     x[k+1|k] = RK4 step of f from x[k|k] over Ts, u held
     *     S[k+1|k] = the triangular factor of [Φ·S[k|k], Q½], so P[k+1|k] = Φ·P[k|k]·Φᵀ + Q
     *
     * @param u Input held over the step; the next updates' h sees it too
     */
    constexpr void predict(const ColVec<NU, T>& u) {
        const Matrix<NX, NX, T> A = detail::dynamics_state_jacobian<NX, NU, T>(model_, x, u, detail::default_jacobian_step<T>());
        const Matrix<NX, NX, T> Phi = detail::rk4_transition(A, Ts);
        x = detail::rk4_step<NX, NU, T>(model_, x, u, Ts);
        u_ = u;
        detail::sqrt_time_update<true>(S, Phi, Gamma);
    }

    /**
     * @brief Measurement update with one reading y = h(x, u) + v of any size, u the input of the
     *        last predict
     *
     *     C = ∂h/∂x at (x[k|k−1], u)
     *     R = L·Lᵀ                                  (Cholesky: R is symmetric positive semidefinite)
     *     ν̃ = L⁻¹·(y − h(x[k|k−1], u)),  C̃ = L⁻¹·C   (unit, uncorrelated noise)
     *     each row c̃ of C̃ in turn:  f = Sᵀ·c̃ᵀ,  α = 1 + fᵀf,  x ← x + S·f·ν̃/α,  S ← S·T   [Carlson]
     *
     * which is K = P·Cᵀ·(C·P·Cᵀ + R)⁻¹, x ← x + K·ν and P ← (I − K·C)·P, on the factor.
     *
     * @param h Measurement y = h(x, u) (Measurement), a struct or a lambda, with an optional
     *          analytic jacobian(x, u). A jacobian that returns C as a PatternMatrix
     *          (PatternMeasurement) has no term for its zeros and no multiply for its ones in C·S.
     * @param y Reading (NY rows)
     * @param R Covariance of the reading's noise, positive semidefinite: a zero-variance row is a perfect measurement
     * @return false if R is not positive semidefinite, C·P·Cᵀ + R is singular or the reading is not finite: the estimate and the covariance are unchanged
     */
    template<typename H, size_t NY>
        requires Measurement<H, NX, NU, NY, T>
    [[nodiscard]] constexpr bool update(const H& h, const ColVec<NY, T>& y, const damp::type_identity_t<Matrix<NY, NY, T>>& R) {
        const ColVec<NY, T> nu = y - h(x, u_);
        const auto          noise = detail::reading_noise_factor(R); // an R that is not a covariance is refused below
        if constexpr (PatternMeasurement<H, NX, NU, NY, T>) {
            return detail::sqrt_measurement_update(S, x, nu, h.jacobian(x, u_), noise);
        } else {
            const Matrix<NY, NX, T> C = detail::measurement_state_jacobian<NX, NU, NY, T>(h, x, u_, detail::default_jacobian_step<T>());
            return detail::sqrt_measurement_update(S, x, nu, C, noise);
        }
    }

    /** @brief The state estimate x̂. */
    [[nodiscard]] constexpr const ColVec<NX, T>& state() const { return x; }
    /** @brief The estimate covariance P = S Sᵀ, formed from the factor on each call: exactly symmetric, rounded to eps·‖P‖. */
    [[nodiscard]] constexpr Matrix<NX, NX, T> covariance() const { return detail::covariance_of_factor(S); }
    /** @brief The factor S the filter carries, P = S Sᵀ: lower triangular with a non-negative diagonal. */
    [[nodiscard]] constexpr const Matrix<NX, NX, T>& covariance_factor() const { return S; }
    /** @brief The process-noise covariance Q. */
    [[nodiscard]] constexpr const Matrix<NX, NX, T>& process_noise_covariance() const { return Q; }
    /** @brief The process model. */
    [[nodiscard]] constexpr const Model& model() const { return model_; }
    /** @brief Sample time [s] that predict integrates over. */
    [[nodiscard]] constexpr T sample_time() const { return Ts; }

    /**
     * @brief Overwrite the estimate between updates (a clamp to a physical bound, an angle wrap);
     *        P and the held input are kept
     *
     * With sequential scalar readings, clamping between the updates enforces the bound each
     * reading sees, which a clamp after one fused vector update cannot.
     */
    constexpr void set_state(const ColVec<NX, T>& x_new) { x = x_new; }

    /**
     * @brief Re-initialize the estimate and its covariance and forget the held input. The model, Ts and Q stay.
     *
     * @p P0 is factored here and may be singular. One that is not positive semidefinite or not
     * finite does not compile when constant-evaluated; at run time the covariance becomes the
     * identity.
     */
    constexpr void reset(const ColVec<NX, T>& x0, const Matrix<NX, NX, T>& P0) {
        x = x0;
        S = detail::covariance_factor_or_identity(P0);
        u_ = ColVec<NU, T>{};
    }

private:
    Model             model_{}; ///< Continuous dynamics ẋ = f(x, u)
    ColVec<NX, T>     x{};      ///< Estimate
    Matrix<NX, NX, T> S{};      ///< Lower-triangular factor of the estimate covariance, P = S Sᵀ
    Matrix<NX, NX, T> Q{};      ///< Process-noise covariance per step
    Matrix<NX, NX, T> Gamma{};  ///< Lower-triangular factor of Q, factored once
    T                 Ts{};     ///< Sample time [s]; 0 (default-constructed) holds the state in predict
    ColVec<NU, T>     u_{};     ///< The input the last predict held, for the readings' h
};

/** @brief Deduce the filter from its model, a struct or a captureless lambda: `ExtendedKalmanFilter ekf{f, {.Ts = Ts, .Q = Q}}`. */
template<typename Model>
ExtendedKalmanFilter(const Model&, const ExtendedKalmanFilterConfig<Model>&) -> ExtendedKalmanFilter<Model>;

} // namespace damp
