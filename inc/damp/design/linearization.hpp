// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file linearization.hpp
 * @brief Nonlinear model linearization (design-time)
 *
 * @defgroup linearization Nonlinear Model Linearization
 * @brief Jacobian linearization of a Dynamics / Measurement model to continuous A/B/C/D at an operating point
 *
 * Linearizes the continuous model of concepts.hpp,
 *
 *   ẋ = f(x, u)   (Dynamics)
 *   y = h(x, u)   (Measurement)
 *
 * around an operating point (x_op, u_op), returning the Jacobians
 *
 *   A = ∂f/∂x,  B = ∂f/∂u,  C = ∂h/∂x,  D = ∂h/∂u
 *
 * from the model's own jacobian(x, u) when it has one (AnalyticDynamics, AnalyticMeasurement),
 * else by central differences. The same f and h are the model of ExtendedKalmanFilter and
 * UnscentedKalmanFilter and the plant of the continuous sim::simulate* functions.
 *
 * Primary API: @c damp::design::linearize (same pattern as @c design::pid, @c design::place).
 */

#include <cstddef>

#include "damp/concepts.hpp"
#include "damp/design/status.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/dynamics.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace design {

/**
 * @brief Result of nonlinear operating-point linearization
 */
template<size_t NX, size_t NU, size_t NY, typename T = double>
struct [[nodiscard]] LinearizationResult {
    Matrix<NX, NX, T> A{}; ///< State Jacobian A = ∂f/∂x
    Matrix<NX, NU, T> B{}; ///< Input Jacobian B = ∂f/∂u
    Matrix<NY, NX, T> C{}; ///< Output Jacobian C = ∂h/∂x
    Matrix<NY, NU, T> D{}; ///< Feedthrough Jacobian D = ∂h/∂u

    ColVec<NX, T> x_op{}; ///< State operating point
    ColVec<NU, T> u_op{}; ///< Input operating point
    ColVec<NY, T> y_op{}; ///< Output at operating point

    T epsilon{}; ///< Relative central-difference step (hᵢ = ε·max(|pᵢ|, 1)); unused for Jacobians the model supplies

    bool         success{false};                      ///< status == DesignStatus::ok
    DesignStatus status{DesignStatus::not_converged}; ///< ok, or why the linearization is refused

    /** @brief True when the linearization succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The same linearization in the scalar @p U. */
    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return LinearizationResult<NX, NU, NY, U>{
            A.template as<U>(),
            B.template as<U>(),
            C.template as<U>(),
            D.template as<U>(),
            x_op.template as<U>(),
            u_op.template as<U>(),
            y_op.template as<U>(),
            static_cast<U>(epsilon),
            success,
            status
        };
    }

    /**
     * @brief Convert Jacobians into a continuous-time state-space model
     */
    [[nodiscard]] constexpr StateSpace<NX, NU, NY, T> to_state_space() const {
        return StateSpace<NX, NU, NY, T>{
            .A = A,
            .B = B,
            .C = C,
            .D = D,
            .Ts = T{0}
        };
    }
};

/**
 * @brief Linearize a Dynamics / Measurement model about an operating point
 *
 * A and B come from the dynamics' jacobian(x, u) when it has one, C and D from the
 * measurement's; the rest by central differences.
 *
 * The central-difference step hᵢ = ε·max(|pᵢ|, 1) has a floor of one unit of each coordinate,
 * so a coordinate much smaller than one unit at the operating point (a gap of a few millimetres
 * in metres) is stepped by more than its own size and its column comes out wrong, the more so
 * in float, where ∛eps ≈ 4.9e-3. State such coordinates in units where they are of order one,
 * or give the model an analytic jacobian(x, u), which no step affects.
 *
 * @param dynamics Continuous dynamics ẋ = f(x, u) (Dynamics)
 * @param output   Measurement y = h(x, u) (Measurement)
 * @param x_op     State operating point
 * @param u_op     Input operating point
 * @param epsilon  Relative central-difference step: hᵢ = ε·max(|pᵢ|, 1) per coordinate (default ∛eps; ≤ 0 → default)
 * @return LinearizationResult with A, B, C, D Jacobians; refused with DesignStatus::invalid_parameter
 *         for a non-finite @p epsilon and DesignStatus::bad_weights when the operating point, the
 *         output there or a Jacobian entry is not finite
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, typename F, typename H>
    requires Dynamics<F, NX, NU, T> && Measurement<H, NX, NU, NY, T>
[[nodiscard]] constexpr LinearizationResult<NX, NU, NY, T> linearize(
    const F&             dynamics,
    const H&             output,
    const ColVec<NX, T>& x_op,
    const ColVec<NU, T>& u_op,
    T                    epsilon = ::damp::detail::default_jacobian_step<T>()
) {
    using Result = LinearizationResult<NX, NU, NY, T>;
    if (!damp::isfinite(epsilon)) {
        return ::damp::detail::refused<Result>(DesignStatus::invalid_parameter);
    }
    if (epsilon <= T{0}) {
        epsilon = ::damp::detail::default_jacobian_step<T>();
    }
    const auto finite = [](T v) { return damp::isfinite(v); };
    if (!mat::all_of(x_op, finite) || !mat::all_of(u_op, finite)) {
        return ::damp::detail::refused<Result>(DesignStatus::bad_weights);
    }
    const DynamicsJacobian<NX, NU, T>        fj = ::damp::detail::dynamics_jacobian<NX, NU, T>(dynamics, x_op, u_op, epsilon);
    const MeasurementJacobian<NX, NU, NY, T> hj = ::damp::detail::measurement_jacobian<NX, NU, NY, T>(output, x_op, u_op, epsilon);
    Result                                   out{
        .A = fj.A,
        .B = fj.B,
        .C = hj.C,
        .D = hj.D,
        .x_op = x_op,
        .u_op = u_op,
        .y_op = output(x_op, u_op),
        .epsilon = epsilon,
    };
    const bool finite_result = mat::all_of(out.A, finite) && mat::all_of(out.B, finite) && mat::all_of(out.C, finite)
                            && mat::all_of(out.D, finite) && mat::all_of(out.y_op, finite);
    return ::damp::detail::set_status(out, finite_result ? DesignStatus::ok : DesignStatus::bad_weights);
}

/**
 * @brief Linearize nonlinear dynamics with identity output y = x
 *
 * Convenience overload for state-feedback workflows.
 */
template<size_t NX, size_t NU, typename T = double, typename F>
    requires Dynamics<F, NX, NU, T>
[[nodiscard]] constexpr LinearizationResult<NX, NU, NX, T> linearize(
    const F&             dynamics,
    const ColVec<NX, T>& x_op,
    const ColVec<NU, T>& u_op,
    T                    epsilon = ::damp::detail::default_jacobian_step<T>()
) {
    const auto output = [](const ColVec<NX, T>& x, const ColVec<NU, T>& /*u*/) -> ColVec<NX, T> {
        return x;
    };

    return linearize<NX, NU, NX, T>(dynamics, output, x_op, u_op, epsilon);
}

} // namespace design
} // namespace damp
