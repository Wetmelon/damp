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
 * @ingroup design
 * @brief Finite-difference Jacobian linearization to continuous A/B/C/D at an operating point
 *
 * Provides numerical linearization for nonlinear plants:
 *
 *   x_dot = f(x, u) or f(t, x, u)
 *   y     = h(x) or h(x, u)
 *
 * around an operating point (x_op, u_op), returning Jacobians:
 *
 *   A = ∂f/∂x,  B = ∂f/∂u,  C = ∂h/∂x,  D = ∂h/∂u
 *
 * using central finite differences.
 *
 * Primary API: @c damp::design::linearize (same pattern as @c design::pid, @c design::place).
 */

#include <concepts>
#include <cstddef>
#include <limits>

#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace design {

/**
 * @brief Result of nonlinear operating-point linearization
 */
template<size_t NX, size_t NU, size_t NY, typename T = double>
struct LinearizationResult {
    Matrix<NX, NX, T> A{}; ///< State Jacobian A = ∂f/∂x
    Matrix<NX, NU, T> B{}; ///< Input Jacobian B = ∂f/∂u
    Matrix<NY, NX, T> C{}; ///< Output Jacobian C = ∂h/∂x
    Matrix<NY, NU, T> D{}; ///< Feedthrough Jacobian D = ∂h/∂u

    ColVec<NX, T> x_op{}; ///< State operating point
    ColVec<NU, T> u_op{}; ///< Input operating point
    ColVec<NY, T> y_op{}; ///< Output at operating point

    T epsilon{}; ///< Finite-difference perturbation used

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
            static_cast<U>(epsilon)
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

namespace detail {

template<size_t NX, size_t NU, typename T, typename Dynamics>
[[nodiscard]] constexpr ColVec<NX, T>
eval_dynamics(const Dynamics& dynamics, const ColVec<NX, T>& x, const ColVec<NU, T>& u) {
    if constexpr (requires { { dynamics(x, u) } -> std::convertible_to<ColVec<NX, T>>; }) {
        return dynamics(x, u);
    } else {
        static_assert(
            requires { { dynamics(T{0}, x, u) } -> std::convertible_to<ColVec<NX, T>>; },
            "Dynamics must be callable as f(x, u) or f(t, x, u)."
        );
        return dynamics(T{0}, x, u);
    }
}

template<size_t NX, size_t NU, size_t NY, typename T, typename Output>
[[nodiscard]] constexpr ColVec<NY, T>
eval_output(const Output& output, const ColVec<NX, T>& x, const ColVec<NU, T>& u) {
    if constexpr (requires { { output(x, u) } -> std::convertible_to<ColVec<NY, T>>; }) {
        return output(x, u);
    } else {
        static_assert(
            requires { { output(x) } -> std::convertible_to<ColVec<NY, T>>; },
            "Output must be callable as h(x, u) or h(x)."
        );
        return output(x);
    }
}

template<typename T>
[[nodiscard]] constexpr T default_linearization_epsilon() {
    // Central differences: optimal step ~ ε^(1/3) balances the O(h²) truncation
    // error against round-off, unlike the ε^(1/2) of one-sided differences.
    return damp::cbrt(std::numeric_limits<T>::epsilon());
}

} // namespace detail

/**
 * @brief Linearize nonlinear dynamics and output maps about an operating point
 *
 * Computes central finite-difference Jacobians for:
 * - dynamics f: x_dot = f(x, u) or f(t, x, u)
 * - output   h: y = h(x) or h(x, u)
 *
 * @param dynamics Nonlinear dynamics callable
 * @param output   Nonlinear output callable
 * @param x_op     State operating point
 * @param u_op     Input operating point
 * @param epsilon  Finite-difference perturbation
 * @return LinearizationResult with A, B, C, D Jacobians
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, typename Dynamics, typename Output>
[[nodiscard]] constexpr LinearizationResult<NX, NU, NY, T> linearize(
    const Dynamics&      dynamics,
    const Output&        output,
    const ColVec<NX, T>& x_op,
    const ColVec<NU, T>& u_op,
    T                    epsilon = detail::default_linearization_epsilon<T>()
) {
    if (epsilon <= T{0}) {
        epsilon = detail::default_linearization_epsilon<T>();
    }

    LinearizationResult<NX, NU, NY, T> result{};
    result.x_op = x_op;
    result.u_op = u_op;
    result.y_op = detail::eval_output<NX, NU, NY, T>(output, x_op, u_op);
    result.epsilon = epsilon;

    const T inv_2eps = T{1} / (T{2} * epsilon);

    for (size_t i = 0; i < NX; ++i) {
        ColVec<NX, T> x_plus = x_op;
        ColVec<NX, T> x_minus = x_op;
        x_plus(i, 0) += epsilon;
        x_minus(i, 0) -= epsilon;

        const auto f_plus = detail::eval_dynamics<NX, NU, T>(dynamics, x_plus, u_op);
        const auto f_minus = detail::eval_dynamics<NX, NU, T>(dynamics, x_minus, u_op);

        const auto y_plus = detail::eval_output<NX, NU, NY, T>(output, x_plus, u_op);
        const auto y_minus = detail::eval_output<NX, NU, NY, T>(output, x_minus, u_op);

        for (size_t r = 0; r < NX; ++r) {
            result.A(r, i) = (f_plus(r, 0) - f_minus(r, 0)) * inv_2eps;
        }
        for (size_t r = 0; r < NY; ++r) {
            result.C(r, i) = (y_plus(r, 0) - y_minus(r, 0)) * inv_2eps;
        }
    }

    for (size_t j = 0; j < NU; ++j) {
        ColVec<NU, T> u_plus = u_op;
        ColVec<NU, T> u_minus = u_op;
        u_plus(j, 0) += epsilon;
        u_minus(j, 0) -= epsilon;

        const auto f_plus = detail::eval_dynamics<NX, NU, T>(dynamics, x_op, u_plus);
        const auto f_minus = detail::eval_dynamics<NX, NU, T>(dynamics, x_op, u_minus);

        const auto y_plus = detail::eval_output<NX, NU, NY, T>(output, x_op, u_plus);
        const auto y_minus = detail::eval_output<NX, NU, NY, T>(output, x_op, u_minus);

        for (size_t r = 0; r < NX; ++r) {
            result.B(r, j) = (f_plus(r, 0) - f_minus(r, 0)) * inv_2eps;
        }
        for (size_t r = 0; r < NY; ++r) {
            result.D(r, j) = (y_plus(r, 0) - y_minus(r, 0)) * inv_2eps;
        }
    }

    return result;
}

/**
 * @brief Linearize nonlinear dynamics with identity output y = x
 *
 * Convenience overload for state-feedback workflows.
 */
template<size_t NX, size_t NU, typename T = double, typename Dynamics>
[[nodiscard]] constexpr LinearizationResult<NX, NU, NX, T> linearize(
    const Dynamics&      dynamics,
    const ColVec<NX, T>& x_op,
    const ColVec<NU, T>& u_op,
    T                    epsilon = detail::default_linearization_epsilon<T>()
) {
    auto output = [](const ColVec<NX, T>& x) -> ColVec<NX, T> {
        return x;
    };

    return linearize<NX, NU, NX, T>(dynamics, output, x_op, u_op, epsilon);
}

} // namespace design
} // namespace damp
