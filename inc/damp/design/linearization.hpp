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

    T epsilon{}; ///< Relative finite-difference step used (hᵢ = ε·max(|pᵢ|, 1))

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
    return damp::cbrt(eps<T>());
}

/// Central-difference step for coordinate value v: rel · max(|v|, 1), so the
/// perturbation stays resolvable at large operating points.
template<typename T>
[[nodiscard]] constexpr T fd_step(T v, T rel) {
    const T mag = damp::abs(v);
    return rel * ((mag > T{1}) ? mag : T{1});
}

/// Central difference of every output of g along coordinate i of p.
/// The divisor is the realized step (p+h) − (p−h), not 2h, so rounding of
/// p ± h does not bias the quotient.
template<size_t NP, size_t NOUT, typename T, typename Eval>
constexpr void central_column(
    Matrix<NOUT, NP, T>& J,
    const ColVec<NP, T>& p,
    size_t               i,
    T                    rel,
    const Eval&          g
) {
    ColVec<NP, T> plus = p;
    ColVec<NP, T> minus = p;
    const T       h = fd_step(p(i, 0), rel);
    plus(i, 0) += h;
    minus(i, 0) -= h;
    const T    inv_span = T{1} / (plus(i, 0) - minus(i, 0));
    const auto g_plus = g(plus);
    const auto g_minus = g(minus);
    for (size_t r = 0; r < NOUT; ++r) {
        J(r, i) = (g_plus(r, 0) - g_minus(r, 0)) * inv_span;
    }
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
 * @param epsilon  Relative central-difference step: hᵢ = ε·max(|pᵢ|, 1) per coordinate (default ∛eps; ≤ 0 → default)
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

    const auto f_of_x = [&](const ColVec<NX, T>& x) { return detail::eval_dynamics<NX, NU, T>(dynamics, x, u_op); };
    const auto h_of_x = [&](const ColVec<NX, T>& x) { return detail::eval_output<NX, NU, NY, T>(output, x, u_op); };
    const auto f_of_u = [&](const ColVec<NU, T>& u) { return detail::eval_dynamics<NX, NU, T>(dynamics, x_op, u); };
    const auto h_of_u = [&](const ColVec<NU, T>& u) { return detail::eval_output<NX, NU, NY, T>(output, x_op, u); };

    for (size_t i = 0; i < NX; ++i) {
        detail::central_column(result.A, x_op, i, epsilon, f_of_x);
        detail::central_column(result.C, x_op, i, epsilon, h_of_x);
    }
    for (size_t j = 0; j < NU; ++j) {
        detail::central_column(result.B, u_op, j, epsilon, f_of_u);
        detail::central_column(result.D, u_op, j, epsilon, h_of_u);
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
