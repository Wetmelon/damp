// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file integrator.hpp
 * @brief Fixed-step ODE integrators for simulation
 */

#include <concepts>
#include <cstddef>
#include <type_traits>

#include "damp/matrix/colvec.hpp"
#include "damp/matrix/functions.hpp"
#include "damp/matrix/matrix.hpp"

namespace damp::sim {

/**
 * @struct IntegrationResult
 * @brief Result of an integration step
 *
 * @tparam NX Number of states
 * @tparam T  Scalar type
 */
template<size_t NX, typename T = double>
struct IntegrationResult {
    ColVec<NX, T> x{};     ///< State at the end of the step
    ColVec<NX, T> error{}; ///< Unscaled embedded local-error vector (zeros for fixed-step methods)
};

/**
 * @brief Integrator that reports a genuine embedded local-error estimate
 *
 * Every integrator's `evolve` returns IntegrationResult with an `error`
 * vector; fixed-step methods leave that vector zero. Wrapping those in
 * AdaptiveStepSolver would silently freeze the step size at the initial
 * guess (scale stays 1 forever). Opt-in traits pin the real pairs:
 *
 * - `provides_embedded_error = true`
 * - `error_order` — local error scales as O(h^{error_order}); the adaptive
 *   controller uses exponent `1/error_order`
 *
 * The error vector is the unscaled embedded difference (e.g. x5 − x̂4 for
 * DP45). AdaptiveStepSolver applies SciPy-style atol/rtol weighting.
 *
 * | Method   | Role                         | MATLAB® |
 * |----------|------------------------------|--------|
 * | DP45     | Non-stiff default (order 5)  | ode45  |
 * | RK23     | Cheap non-stiff (order 3)    | ode23  |
 * | TRBDF2   | Stiff / L-stable (order 2)   | ode23tb|
 *
 * Prefer adaptive_solve for one-liners; AdaptiveStepSolver when you
 * need events, zero-crossings, or callbacks.
 *
 * @tparam I  Integrator type
 * @tparam NX Number of states
 * @tparam T  Scalar type
 */
template<typename I, size_t NX, typename T>
concept AdaptiveStepIntegrator = requires {
    requires I::provides_embedded_error;
    { I::error_order } -> std::convertible_to<int>;
    requires I::error_order > 0;
} && requires(const I& i, ColVec<NX, T> (*f)(T, const ColVec<NX, T>&), const ColVec<NX, T>& x, T t, T h) {
    { i.evolve(f, x, t, h) } -> std::convertible_to<IntegrationResult<NX, T>>;
};

/**
 * @brief Discrete-time integrator (no integration, just one step)
 *
 * Applies an already-discretized state-space update directly; use this when
 * A and B are discrete-time matrices rather than a continuous plant.
 *
 * @tparam NX Number of states
 * @tparam T  Scalar type
 */
template<size_t NX, typename T = double>
struct Discrete {
    /**
     * @brief Advance an already-discrete system one step: x_{n+1} = A x_n + B u_n
     *
     * No integration is performed; A and B are treated as discrete-time
     * matrices (e.g. already produced by @ref damp::discretize).
     *
     * @param A     Discrete system matrix
     * @param B     Discrete input matrix
     * @param x     Current state
     * @param u     Current input
     *
     * @return      Next state and error estimate (always 0)
     */
    template<size_t NU>
    constexpr IntegrationResult<NX, T> evolve(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B, const ColVec<NX, T>& x, const ColVec<NU, T>& u) const {
        return {A * x + B * u, {}};
    }
};

/**
 * @brief Exact integrator for LTI systems
 *
 * Uses the matrix exponential of the augmented [A B; 0 0] matrix to advance
 * the state with zero truncation error, assuming u is held constant over h
 * (zero-order hold). Exact but more expensive per step than the explicit
 * methods below.
 *
 * @note Compare with MATLAB®'s expm() applied to the augmented matrix, or c2d(sys, h, 'zoh') followed by one discrete step.
 *
 * @tparam NX Number of states
 * @tparam T  Scalar type
 */
template<size_t NX, typename T = double>
struct Exact {

    /**
     * @brief Evolve LTI system state exactly over step h
     *
     * @param A     System matrix
     * @param B     Input matrix
     * @param x     Current state
     * @param u     Input (held constant over the duration of the step)
     * @param h     [s] Step size
     *
     * @return      Next state and error estimate (always 0)
     */
    template<size_t NU>
    constexpr IntegrationResult<NX, T> evolve(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B, const ColVec<NX, T>& x, const ColVec<NU, T>& u, T h) const {
        // Augmented matrix
        Matrix<NX + NU, NX + NU, T> M = Matrix<NX + NU, NX + NU, T>::zeros();

        M.template block<NX, NX>(0, 0) = A;
        M.template block<NX, NU>(0, NX) = B;

        // Matrix exponential
        auto expM = mat::expm(M * h);

        // Extract blocks
        auto expAh = expM.template block<NX, NX>(0, 0);
        auto Gamma = expM.template block<NX, NU>(0, NX);

        // Exact update
        ColVec<NX, T> x_next = expAh * x + Gamma * u;

        return {x_next, {}};
    }
};

/**
 * @brief Forward Euler integrator
 *
 * Fixed-step, explicit, 1st-order accurate:
 *
 *     x_{n+1} = x_n + h f(t_n, x_n)
 *
 * Cheapest method available (one function evaluation per step) but least
 * accurate and only conditionally stable; use small step sizes, especially
 * for stiff systems.
 *
 * @note Compare with MATLAB®/Simulink®'s ode1() fixed-step (explicit Euler) solver.
 *
 * @tparam NX Number of states
 * @tparam T  Scalar type
 */
template<size_t NX, typename T = double>
struct ForwardEuler {

    /**
     * @brief Evolve nonlinear system state using Forward Euler
     *
     * @param f     Right-hand side function: dx/dt = f(t, x)
     * @param x     Current state
     * @param t     [s] Current time (for passing to f)
     * @param h     [s] Step size
     *
     * @return      Next state and estimated error (always 0)
     */
    constexpr IntegrationResult<NX, T> evolve(const auto& f, const ColVec<NX, T>& x, T t, T h) const {
        k1 = f(t, x);
        return {x + h * k1, {}};
    }

    /**
     * @brief Evolve LTI system state using Forward Euler: x_{n+1} = x_n + h(A x_n + B u_n)
     *
     * @param A     System matrix
     * @param B     Input matrix
     * @param x     Current state
     * @param u     Current input
     * @param h     [s] Step size
     *
     * @return      Next state and estimated error (always 0)
     */
    template<size_t NU>
    constexpr IntegrationResult<NX, T> evolve(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B, const ColVec<NX, T>& x, const ColVec<NU, T>& u, T h) const {
        x_next = x + h * (A * x + B * u);
        return {x_next, {}};
    }

private:
    mutable ColVec<NX, T> k1, x_next;
};

/**
 * @brief Semi-implicit (symplectic) Euler for mechanical systems
 *
 * Fixed-step, explicit, 1st-order accurate. State is stacked as
 * `x = [q; v]` with `NQ = NX/2` configuration coordinates and `NQ` velocities.
 * One step is
 *
 *     v_{n+1} = v_n + h a(t_n, q_n, v_n)
 *     q_{n+1} = q_n + h v_{n+1}
 *
 * where `a` is the lower half of `f(t, x)`. Position is advanced with the
 * *updated* velocity, so the method is symplectic on separable Hamiltonians
 * and keeps energy error bounded over long horizons (unlike Forward Euler,
 * which drifts monotonically). Requires even `NX` and the mechanical
 * partitioning above — not a general ODE solver.
 *
 * @see "Geometric Numerical Integration" (Hairer, Lubich & Wanner, 2006), §1.1
 *
 * @tparam NX Number of states (must be even: NQ + NQ)
 * @tparam T  Scalar type
 */
template<size_t NX, typename T = double>
struct SymplecticEuler {
    static_assert(NX % 2 == 0, "SymplecticEuler requires even NX (state = [q; v])");
    static constexpr size_t NQ = NX / 2;

    /**
     * @brief Evolve LTI mechanical system one symplectic-Euler step
     *
     * @param A     System matrix
     * @param B     Input matrix
     * @param x     Current state [q; v]
     * @param u     Input (held constant over the step)
     * @param h     [s] Step size
     *
     * @return      Next state and error estimate (always 0)
     */
    template<size_t NU>
    constexpr IntegrationResult<NX, T> evolve(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B, const ColVec<NX, T>& x, const ColVec<NU, T>& u, T h) const {
        auto f = [&](T, const ColVec<NX, T>& xx) { return A * xx + B * u; };
        return evolve(f, x, T{0}, h);
    }

    /**
     * @brief Evolve nonlinear mechanical system one symplectic-Euler step
     *
     * @param f     Right-hand side: dx/dt = f(t, x) = [v_expected; a]; only `a` is used
     * @param x     Current state [q; v]
     * @param t     [s] Current time
     * @param h     [s] Step size
     *
     * @return      Next state and error estimate (always 0)
     */
    template<typename F>
    constexpr IntegrationResult<NX, T> evolve(const F& f, const ColVec<NX, T>& x, T t, T h) const {
        fx = f(t, x);
        x_next = x;
        for (size_t i = 0; i < NQ; ++i) {
            x_next[NQ + i] = x[NQ + i] + h * fx[NQ + i]; // v ← v + h a
        }
        for (size_t i = 0; i < NQ; ++i) {
            x_next[i] = x[i] + h * x_next[NQ + i]; // q ← q + h v_new
        }
        return {x_next, {}};
    }

private:
    mutable ColVec<NX, T> fx, x_next;
};

/**
 * @brief Backward Euler integrator
 *
 * Implicit, 1st-order accurate:
 *
 *     x_{n+1} = x_n + h f(t_{n+1}, x_{n+1})
 *
 * A-stable, so it remains stable at large step sizes for stiff systems where
 * Forward Euler would diverge, at the cost of a linear solve (LTI case) or
 * fixed-point iteration (nonlinear case) per step.
 *
 * @note Compare with MATLAB®/Simulink®'s ode1be() / implicit Euler fixed-step solver.
 * @see "Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems" (Hairer & Wanner, 1996)
 *
 * @tparam NX Number of states
 * @tparam T  Scalar type
 */
template<size_t NX, typename T = double>
struct BackwardEuler {

    /**
     * @brief Evolve LTI system state using Backward Euler: (I − hA)x_{n+1} = x_n + hBu_n
     *
     * @param A     System matrix
     * @param B     Input matrix
     * @param x     Current state
     * @param u     Current input
     * @param h     [s] Step size
     *
     * @return      Next state and estimated error (always 0)
     */
    template<size_t NU>
    constexpr IntegrationResult<NX, T> evolve(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B, const ColVec<NX, T>& x, const ColVec<NU, T>& u, T h) const {
        Matrix I = Matrix<NX, NX, T>::identity();
        // (I − hA)·x_next = x + hBu  (Backward Euler implicit step).
        x_next = *mat::solve(I - h * A, x + h * B * u);
        return {x_next, {}};
    }

    /**
     * @brief Evolve nonlinear system state using Backward Euler, solved by
     *        fixed-point iteration on y_{n+1} = x_n + h f(t_{n+1}, y_{n+1})
     *
     * @param f     Right-hand side function: dx/dt = f(t, x)
     * @param x     Current state
     * @param t     [s] Current time (for passing to f)
     * @param h     [s] Step size
     *
     * @return      Next state and estimated error (always 0)
     */
    constexpr IntegrationResult<NX, T> evolve(const auto& f, const ColVec<NX, T>& x, T t, T h) const {
        const size_t max_iter = 10;
        const double tol = 1e-10;

        // Initial guess: explicit Euler
        ColVec<NX, T> y = x + h * f(t, x);
        for (size_t i = 0; i < max_iter; ++i) {
            y_next = x + h * f(t + h, y);
            if ((y_next - y).norm() <= tol) {
                y = y_next;
                break;
            }
            y = y_next;
        }
        return {y, {}};
    }

private:
    mutable ColVec<NX, T> y_next, x_next;
};

/**
 * @brief Backward Differentiation Formula 2 (BDF2) integrator
 *
 * Implicit, 2-step, 2nd-order accurate:
 *
 *     x_{n+1} = (4/3)x_n − (1/3)x_{n-1} + (2/3)h f(t_{n+1}, x_{n+1})
 *
 * A-stable multistep method; more accurate than Backward Euler for the same
 * step size on stiff problems, at the cost of keeping one step of history.
 * The first step has no x_{n-1} yet and falls back to Backward Euler.
 *
 * @note Compare with MATLAB®'s ode15s() in fixed-order-2, fixed-step mode.
 * @see "Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems" (Hairer & Wanner, 1996)
 *
 * @tparam NX Number of states
 * @tparam T  Scalar type
 */
template<size_t NX, typename T = double>
struct BDF2 {
    BDF2() : first_step(true) {}

    // BDF2 coefficients: x_{n+1} = (4/3)*x_n - (1/3)*x_{n-1} + (2/3)*h*f(t_{n+1}, x_{n+1})
    static constexpr double c0 = 4.0 / 3.0;  // Coefficient for x_n
    static constexpr double c1 = -1.0 / 3.0; // Coefficient for x_{n-1}
    static constexpr double c2 = 2.0 / 3.0;  // Coefficient for h*f

    /**
     * @brief Evolve LTI system state using BDF2
     *
     * First step falls back to Backward Euler (no x_{n-1} history yet);
     * subsequent steps use the BDF2 formula above, solved implicitly since
     * x_{n+1} appears on both sides.
     *
     * @param A     System matrix
     * @param B     Input matrix
     * @param x     Current state
     * @param u     Current input
     * @param h     [s] Step size
     *
     * @return      Next state and estimated error (always 0)
     */
    template<size_t NU>
    constexpr IntegrationResult<NX, T> evolve(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B, const ColVec<NX, T>& x, const ColVec<NU, T>& u, T h) const {
        Matrix I = Matrix<NX, NX, T>::identity();

        if (first_step) {
            // Use Backward Euler for first step
            ColVec<NX, T> x_next = *mat::solve(I - h * A, x + h * B * u);
            x_prev = x;
            first_step = false;
            return {x_next, {}};
        }

        // BDF2: (I - (2/3)*h*A)*x_{n+1} = (4/3)*x_n - (1/3)*x_{n-1} + (2/3)*h*B*u
        Matrix        lhs = I - c2 * h * A;
        ColVec        rhs = c0 * x + c1 * x_prev + c2 * h * B * u;
        ColVec<NX, T> x_next = *mat::solve(lhs, rhs);

        x_prev = x;
        return {x_next, {}};
    }

    /**
     * @brief Evolve nonlinear system state using BDF2, solved by fixed-point
     *        iteration on y_{n+1} = (4/3)x_n − (1/3)x_{n-1} + (2/3)h f(t_{n+1}, y_{n+1})
     *
     * First step falls back to Backward Euler (no x_{n-1} history yet).
     *
     * @param f     Right-hand side function: dx/dt = f(t, x)
     * @param x     Current state
     * @param t     [s] Current time (for passing to f)
     * @param h     [s] Step size
     *
     * @return      Next state and estimated error (always 0)
     */
    template<typename F>
    constexpr IntegrationResult<NX, T> evolve(const F& f, const ColVec<NX, T>& x, T t, T h) const {
        constexpr size_t max_iter = 20;
        constexpr double tol = 1e-10;

        if (first_step) {
            // Use Backward Euler for first step
            ColVec y = x + h * f(t, x); // Initial guess
            for (size_t i = 0; i < max_iter; ++i) {
                ColVec y_next = x + h * f(t + h, y);
                if ((y_next - y).norm() <= tol) {
                    y = y_next;
                    break;
                }
                y = y_next;
            }
            x_prev = x;
            first_step = false;
            return {y, {}};
        }

        // BDF2: x_{n+1} = (4/3)*x_n - (1/3)*x_{n-1} + (2/3)*h*f(t_{n+1}, x_{n+1})
        // Solve using fixed-point iteration
        ColVec y = c0 * x + c1 * x_prev; // Initial guess without implicit term

        for (size_t i = 0; i < max_iter; ++i) {
            ColVec y_next = c0 * x + c1 * x_prev + c2 * h * f(t + h, y);
            if ((y_next - y).norm() <= tol) {
                y = y_next;
                break;
            }
            y = y_next;
        }

        x_prev = x;
        return {y, {}};
    }

    /**
     * @brief Discard multistep history and restart from Backward Euler on
     *        the next call to evolve() (e.g. after a discontinuity in x)
     */
    void reset() const {
        first_step = true;
    }

private:
    mutable ColVec<NX, T> x_prev;     // Previous state (for multistep)
    mutable bool          first_step; // Track if this is the first step
};

/**
 * @brief TR-BDF2 composite integrator — the stiff adaptive pair (ode23tb)
 *
 * One-step, 2nd-order, L-stable (γ = 2 − √2). Stage 1 advances to t_n + γh by
 * the trapezoidal rule; stage 2 finishes with BDF2 on {x_n, x_γ, x_{n+1}}.
 * Both stages share the iteration matrix I − (γh/2)A (two LTI solves, same
 * LHS). No multistep history — no `reset()` after discontinuities (unlike
 * BDF2). Prefer this over Trapezoidal on very stiff modes (Trapezoidal only
 * A-stable and can ring) and over BDF2 when the plant switches.
 *
 * Embedded 3rd-order companion (Hosea–Shampine): reuses f at {t_n, t_γ, t_{n+1}}
 * so the error estimate is free of extra stages; local error ~ O(h³). Usable
 * fixed-step by ignoring `error`, or under AdaptiveStepSolver /
 * adaptive_solve for stiff adaptive integration.
 *
 * @note Compare with MATLAB®'s ode23tb().
 * @see "Analysis and implementation of TR-BDF2" (Hosea & Shampine, 1996)
 * @see Bank et al., "Transient simulation of silicon devices and circuits" (1985)
 *
 * @tparam NX Number of states
 * @tparam T  Scalar type
 */
template<size_t NX, typename T = double>
struct TRBDF2 {
    /// @brief Opt-in for @ref AdaptiveStepIntegrator
    static constexpr bool provides_embedded_error = true;
    /// @brief Local error scales as O(h³) (embedded 2/3 pair)
    static constexpr int error_order = 3;

    /**
     * @brief Evolve LTI system state using TR-BDF2
     *
     * @param A     System matrix
     * @param B     Input matrix
     * @param x     Current state
     * @param u     Input (held constant over the step)
     * @param h     [s] Step size
     *
     * @return      2nd-order next state and embedded 3rd-order error estimate
     */
    template<size_t NU>
    constexpr IntegrationResult<NX, T> evolve(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B, const ColVec<NX, T>& x, const ColVec<NU, T>& u, T h) const {
        // γ = 2 − √2 ⇒ L-stable; both stages share M = I − (γh/2) A.
        const T g = T{2} - damp::sqrt(T{2});
        const T gh2 = g * h * static_cast<T>(0.5);
        const T c_g = T{1} / (g * (T{2} - g));
        const T c_n = -((T{1} - g) * (T{1} - g)) / (g * (T{2} - g));

        const Matrix I = Matrix<NX, NX, T>::identity();
        const Matrix M = I - gh2 * A;

        // Stage 1 — trapezoidal over [t, t+γh]
        const ColVec<NX, T> x_gamma = *mat::solve(M, (I + gh2 * A) * x + (g * h) * (B * u));

        // Stage 2 — BDF2 to t+h
        const ColVec<NX, T> x_next = *mat::solve(M, c_g * x_gamma + c_n * x + gh2 * (B * u));

        const ColVec<NX, T> f0 = A * x + B * u;
        const ColVec<NX, T> fg = A * x_gamma + B * u;
        const ColVec<NX, T> f1 = A * x_next + B * u;
        return IntegrationResult<NX, T>{x_next, embedded_error(h, g, c_g, f0, fg, f1)};
    }

    /**
     * @brief Evolve nonlinear system state using TR-BDF2 with embedded error
     *
     * @param f     Right-hand side function: dx/dt = f(t, x)
     * @param x     Current state
     * @param t     [s] Current time
     * @param h     [s] Step size
     *
     * @return      2nd-order next state and unscaled embedded error x₂ − x̂₃
     */
    template<typename F>
    constexpr IntegrationResult<NX, T> evolve(const F& f, const ColVec<NX, T>& x, T t, T h) const {
        constexpr size_t max_iter = 20;
        constexpr T      tol = static_cast<T>(1e-10);

        const T g = T{2} - damp::sqrt(T{2});
        const T gh2 = g * h * static_cast<T>(0.5);
        const T c_g = T{1} / (g * (T{2} - g));
        const T c_n = -((T{1} - g) * (T{1} - g)) / (g * (T{2} - g));
        const T t_g = t + g * h;

        // Working vectors are locals (see DP45) so evolve() stays re-entrant.
        const ColVec<NX, T> f0 = f(t, x);
        // Stage 1 — trapezoidal fixed-point: y = x + (γh/2)(f(t,x) + f(t_γ,y))
        ColVec<NX, T> y = x + (g * h) * f0; // explicit-Euler predictor over γh
        for (size_t i = 0; i < max_iter; ++i) {
            const ColVec<NX, T> y_next = x + gh2 * (f0 + f(t_g, y));
            if ((y_next - y).norm() <= tol) {
                y = y_next;
                break;
            }
            y = y_next;
        }
        const ColVec<NX, T> x_gamma = y;

        // Stage 2 — BDF2 fixed-point: y = c_g x_γ + c_n x + (γh/2) f(t+h, y)
        y = c_g * x_gamma + c_n * x;
        for (size_t i = 0; i < max_iter; ++i) {
            const ColVec<NX, T> y_next = c_g * x_gamma + c_n * x + gh2 * f(t + h, y);
            if ((y_next - y).norm() <= tol) {
                y = y_next;
                break;
            }
            y = y_next;
        }
        const ColVec<NX, T> x_next = y;

        const ColVec<NX, T> fg = f(t_g, x_gamma);
        const ColVec<NX, T> f1 = f(t + h, x_next);
        return IntegrationResult<NX, T>{x_next, embedded_error(h, g, c_g, f0, fg, f1)};
    }

private:
    /// Hosea–Shampine embedded 3rd-order companion: e = β − b on stages {0, γ, 1}.
    [[nodiscard]] static constexpr ColVec<NX, T> embedded_error(T h, T g, T c_g, const ColVec<NX, T>& f0, const ColVec<NX, T>& fg, const ColVec<NX, T>& f1) {
        // Order conditions Σb=1, Σb c=1/2, Σb c²=1/3 at c ∈ {0, γ, 1}.
        const T b2 = (T{-1} / T{6}) / (g * (g - T{1}));
        const T b3 = static_cast<T>(0.5) - b2 * g;
        const T b1 = T{1} - b2 - b3;

        // TR-BDF2 weights: x₂ = x + h Σ βᵢ fᵢ with β₁ = β₂ = c_g·γ/2, β₃ = γ/2.
        const T beta1 = c_g * (g * static_cast<T>(0.5));
        const T beta2 = beta1;
        const T beta3 = g * static_cast<T>(0.5);
        return h * ((beta1 - b1) * f0 + (beta2 - b2) * fg + (beta3 - b3) * f1);
    }
};

/**
 * @brief Trapezoidal (Tustin) integrator
 *
 * Implicit, 2nd-order accurate, A-stable:
 *
 *     x_{n+1} = x_n + (h/2)(f(t_n, x_n) + f(t_{n+1}, x_{n+1}))
 *
 * Averages the slope at the start and end of the step; good accuracy and
 * stability across a wide range of problems, including stiff ones, without
 * the numerical damping Backward Euler introduces.
 *
 * @note Compare with MATLAB®'s c2d(sys, h, 'tustin') applied per step.
 * @see "Feedback Control of Dynamic Systems" (Franklin et al., 2015), §8.3
 *
 * @tparam NX Number of states
 * @tparam T  Scalar type
 */
template<size_t NX, typename T = double>
struct Trapezoidal {
    /**
     * @brief Evolve LTI system state using the trapezoidal rule:
     *        (I − ½hA)x_{n+1} = (I + ½hA)x_n + hBu_n
     *
     * @param A     System matrix
     * @param B     Input matrix
     * @param x     Current state
     * @param u     Current input
     * @param h     [s] Step size
     *
     * @return      Next state and estimated error (always 0)
     */
    template<size_t NU>
    constexpr IntegrationResult<NX, T> evolve(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B, const ColVec<NX, T>& x, const ColVec<NU, T>& u, T h) const {
        Matrix I = Matrix<NX, NX, T>::identity();
        // (I − ½hA)·x_next = (I + ½hA)·x + hBu  (trapezoidal / Tustin implicit step).
        x_next = *mat::solve(I - 0.5 * h * A, (I + 0.5 * h * A) * x + h * B * u);
        return {x_next, {}};
    }

    /**
     * @brief Evolve nonlinear system state using the trapezoidal rule
     *
     * Predictor-corrector form (equivalent to Heun's method applied as a
     * fixed-point step rather than solved implicitly):
     *
     *     x_{n+1} = x_n + (h/2)(f(t_n, x_n) + f(t_n+h, x_n + h f(t_n, x_n)))
     *
     * @param f     Right-hand side function: dx/dt = f(t, x)
     * @param x     Current state
     * @param t     [s] Current time (for passing to f)
     * @param h     [s] Step size
     *
     * @return      Next state and estimated error (always 0)
     */
    template<typename F>
    constexpr IntegrationResult<NX, T> evolve(const F& f, const ColVec<NX, T>& x, T t, T h) const {
        k1 = f(t, x);
        xp = x + h * k1;
        k2 = f(t + h, xp);
        return {x + 0.5 * h * (k1 + k2), {}};
    }

private:
    mutable ColVec<NX, T> k1, k2, xp, x_next;
};

/**
 * @brief Classical 4th-order Runge-Kutta (RK4) integrator
 *
 * Fixed-step, explicit, 4th-order accurate. The standard general-purpose
 * choice when no error estimate is needed.
 *
 * @note Compare with MATLAB®/Simulink®'s ode4() fixed-step solver, not the adaptive ode45.
 *
 * @tparam NX Number of states
 * @tparam T  Scalar type
 */
template<size_t NX, typename T = double>
struct RK4 {
    /**
     * @brief Evolve LTI system state using RK4 (wraps the ODE form below)
     *
     * @param A     System matrix
     * @param B     Input matrix
     * @param x     Current state
     * @param u     Input (held constant over the duration of the step)
     * @param h     [s] Step size
     *
     * @return      Next state and estimated error (always 0)
     */
    template<size_t NU>
    constexpr IntegrationResult<NX, T> evolve(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B, const ColVec<NX, T>& x, const ColVec<NU, T>& u, T h) const {
        // RK4 for LTI system
        auto f = [&](double, const ColVec<NX, T>& xk) { return A * xk + B * u; };
        return evolve(f, x, 0.0, h);
    }

    /**
     * @brief Evolve nonlinear system state using RK4
     *
     *     k1 = f(t, x)
     *     k2 = f(t+h/2, x+h/2 k1)
     *     k3 = f(t+h/2, x+h/2 k2)
     *     k4 = f(t+h, x+h k3)
     *     x_{n+1} = x + (h/6)(k1 + 2k2 + 2k3 + k4)
     *
     * @param f     Right-hand side function: dx/dt = f(t, x)
     * @param x     Current state
     * @param t     [s] Current time (for passing to f)
     * @param h     [s] Step size
     *
     * @return      Next state and estimated error (always 0)
     */
    template<typename F>
    constexpr IntegrationResult<NX, T> evolve(const F& f, const ColVec<NX, T>& x, T t, T h) const {
        k1 = f(t, x);
        k2 = f(t + 0.5 * h, x + 0.5 * h * k1);
        k3 = f(t + 0.5 * h, x + 0.5 * h * k2);
        k4 = f(t + h, x + h * k3);

        x_next = x + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
        return {x_next, {}};
    }

private:
    mutable ColVec<NX, T> k1, k2, k3, k4, x_next;
};

/**
 * @brief Dormand-Prince 5(4) adaptive integrator — the ode45 pair
 *
 * Embedded 5th/4th-order explicit pair. Advances with the 5th-order solution
 * (local extrapolation) and returns the unscaled embedded difference
 *
 *     err = x5 − x̂4 = h Σᵢ (bᵢ − b̂ᵢ) kᵢ
 *
 * for step-size control (weighted by atol/rtol in AdaptiveStepSolver).
 * Seven stages; the pair is FSAL (the 7th stage equals the next step's first
 * stage), which this implementation does not cache — each evolve() is
 * stateless so rejected or discontinuous steps need no bookkeeping, at the
 * cost of one extra function evaluation per step. Also usable at a fixed
 * step (the error output is simply ignored), where it is 5th-order accurate.
 *
 * @note Compare with MATLAB®'s ode45(), which uses this exact pair.
 * @see "A family of embedded Runge-Kutta formulae" (Dormand & Prince, 1980),
 *      https://doi.org/10.1016/0771-050X(80)90013-3
 * @see "The MATLAB® ODE Suite" (Shampine & Reichelt, 1997), §2
 *
 * @tparam NX Number of states
 * @tparam T  Scalar type
 */
template<size_t NX, typename T = double>
struct DP45 {
    /// @brief Opt-in for @ref AdaptiveStepIntegrator
    static constexpr bool provides_embedded_error = true;
    /// @brief Local error scales as O(h⁵) (embedded 5(4) pair)
    static constexpr int error_order = 5;

    /**
     * @brief Evolve LTI system state using DP45 (wraps the ODE form below)
     *
     * @param A     System matrix
     * @param B     Input matrix
     * @param x     Current state
     * @param u     Input (held constant over the duration of the step)
     * @param h     [s] Step size
     *
     * @return      5th-order next state and embedded error estimate
     */
    template<size_t NU>
    constexpr IntegrationResult<NX, T> evolve(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B, const ColVec<NX, T>& x, const ColVec<NU, T>& u, T h) const {
        // Evolution of LTI system is a special case of generic ODE solver
        auto f = [&](T, const ColVec<NX, T>& xk) { return A * xk + B * u; };
        return evolve(f, x, 0.0, h);
    }

    /**
     * @brief Evolve nonlinear system state using the Dormand-Prince 5(4) pair,
     *        returning the 5th-order solution and the embedded error estimate
     *
     * @param f     Right-hand side function: dx/dt = f(t, x)
     * @param x     Current state
     * @param t     [s] Current time (for passing to f)
     * @param h     [s] Step size
     *
     * @return      5th-order next state and unscaled embedded error x5 − x̂4
     */
    template<typename F>
    constexpr IntegrationResult<NX, T> evolve(const F& f, const ColVec<NX, T>& x, T t, T h) const {
        // Dormand-Prince 5(4) tableau (precomputed to avoid divisions in hot loop)
        constexpr T c2 = 1.0 / 5.0, c3 = 3.0 / 10.0, c4 = 4.0 / 5.0, c5 = 8.0 / 9.0;
        constexpr T a21 = 1.0 / 5.0;
        constexpr T a31 = 3.0 / 40.0, a32 = 9.0 / 40.0;
        constexpr T a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0;
        constexpr T a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0, a53 = 64448.0 / 6561.0, a54 = -212.0 / 729.0;
        constexpr T a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0, a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0, a65 = -5103.0 / 18656.0;

        // 5th-order weights (b2 = b7 = 0); the 7th stage sits at the solution (FSAL).
        constexpr T b1 = 35.0 / 384.0, b3 = 500.0 / 1113.0, b4 = 125.0 / 192.0, b5 = -2187.0 / 6784.0, b6 = 11.0 / 84.0;

        // Error weights e = b − b̂ against the embedded 4th-order solution (e2 = 0).
        constexpr T e1 = 71.0 / 57600.0, e3 = -71.0 / 16695.0, e4 = 71.0 / 1920.0,
                    e5 = -17253.0 / 339200.0, e6 = 22.0 / 525.0, e7 = -1.0 / 40.0;

        // Stages are locals (not mutable members). Returning views of member
        // stages was implicated in intermittent AVs under g++ 14 MinGW -O3 when
        // multiple adaptive_solve instantiations shared a TU; locals keep each
        // evolve() fully re-entrant and copy-clean.
        const ColVec<NX, T> k1 = f(t, x);
        const ColVec<NX, T> k2 = f(t + c2 * h, x + h * a21 * k1);
        const ColVec<NX, T> k3 = f(t + c3 * h, x + h * (a31 * k1 + a32 * k2));
        const ColVec<NX, T> k4 = f(t + c4 * h, x + h * (a41 * k1 + a42 * k2 + a43 * k3));
        const ColVec<NX, T> k5 = f(t + c5 * h, x + h * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4));
        const ColVec<NX, T> k6 = f(t + h, x + h * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5));

        const ColVec<NX, T> x5 = x + h * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6);
        const ColVec<NX, T> k7 = f(t + h, x5); // FSAL stage — completes the embedded 4th-order comparison

        const ColVec<NX, T> err = h * (e1 * k1 + e3 * k3 + e4 * k4 + e5 * k5 + e6 * k6 + e7 * k7);
        return IntegrationResult<NX, T>{x5, err};
    }
};

/**
 * @brief Heun's method (Improved Euler, RK2) integrator
 *
 * Fixed-step, explicit, 2nd-order accurate predictor-corrector: one Euler
 * predictor step followed by averaging its slope with the corrector slope.
 *
 * @note Compare with MATLAB®'s ode2() fixed-step solver.
 *
 * @tparam NX Number of states
 * @tparam T  Scalar type
 */
template<size_t NX, typename T = double>
struct Heun {
    /**
     * @brief Evolve LTI system state using Heun's method (wraps the ODE form below)
     *
     * @param A     System matrix
     * @param B     Input matrix
     * @param x     Current state
     * @param u     Input (held constant over the duration of the step)
     * @param h     [s] Step size
     *
     * @return      Next state and estimated error (always 0)
     */
    template<size_t NU>
    constexpr IntegrationResult<NX, T> evolve(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B, const ColVec<NX, T>& x, const ColVec<NU, T>& u, T h) const {
        auto f = [&](T, const ColVec<NX, T>& xk) { return A * xk + B * u; };
        return evolve(f, x, 0.0, h);
    }

    /**
     * @brief Evolve nonlinear system state using Heun's method
     *
     *     k1 = f(t, x)
     *     k2 = f(t+h, x+h k1)
     *     x_{n+1} = x + (h/2)(k1 + k2)
     *
     * @param f     Right-hand side function: dx/dt = f(t, x)
     * @param x     Current state
     * @param t     [s] Current time (for passing to f)
     * @param h     [s] Step size
     *
     * @return      Next state and estimated error (always 0)
     */
    template<typename F>
    constexpr IntegrationResult<NX, T> evolve(const F& f, const ColVec<NX, T>& x, T t, T h) const {
        // Heun's method (Improved Euler, RK2)
        // Predictor-corrector: first estimate then correct
        constexpr double c2 = 1.0;
        constexpr double b1 = 0.5, b2 = 0.5;

        k1 = f(t, x);
        k2 = f(t + c2 * h, x + h * k1);
        return {x + h * (b1 * k1 + b2 * k2), {}};
    }

private:
    mutable ColVec<NX, T> k1, k2;
};

/**
 * @brief Bogacki-Shampine 2(3) adaptive integrator
 *
 * Embedded 2nd/3rd-order pair; returns the 3rd-order solution with the
 * difference against the 2nd-order solution as an error estimate. Cheaper
 * per step than DP45 (3 stages plus a reusable FSAL-style evaluation).
 *
 * @note Compare with MATLAB®'s ode23().
 * @see "A 3(2) pair of Runge-Kutta formulas" (Bogacki & Shampine, 1989)
 *
 * @tparam NX Number of states
 * @tparam T  Scalar type
 */
template<size_t NX, typename T = double>
struct RK23 {
    /// @brief Opt-in for @ref AdaptiveStepIntegrator
    static constexpr bool provides_embedded_error = true;
    /// @brief Local error scales as O(h³) (embedded 2(3) pair)
    static constexpr int error_order = 3;

    /**
     * @brief Evolve LTI system state using RK23 (wraps the ODE form below)
     *
     * @param A     System matrix
     * @param B     Input matrix
     * @param x     Current state
     * @param u     Input (held constant over the duration of the step)
     * @param h     [s] Step size
     *
     * @return      Next state and unscaled embedded error x3 − x2
     */
    template<size_t NU>
    constexpr IntegrationResult<NX, T> evolve(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B, const ColVec<NX, T>& x, const ColVec<NU, T>& u, T h) const {
        auto f = [&](double, const ColVec<NX, T>& xk) { return A * xk + B * u; };
        return evolve(f, x, 0.0, h);
    }

    /**
     * @brief Evolve nonlinear system state using RK23, returning both the
     *        3rd-order solution and a 2nd/3rd-order error estimate
     *
     * @param f     Right-hand side function: dx/dt = f(t, x)
     * @param x     Current state
     * @param t     [s] Current time (for passing to f)
     * @param h     [s] Step size
     *
     * @return      3rd-order next state and unscaled embedded error x3 − x2
     */
    template<typename F>
    constexpr IntegrationResult<NX, T> evolve(const F& f, const ColVec<NX, T>& x, T t, T h) const {
        // Bogacki-Shampine (2,3) pair — stage vectors are locals (see DP45).
        constexpr T c2 = static_cast<T>(0.5), c3 = static_cast<T>(0.75);
        constexpr T a21 = static_cast<T>(0.5);
        constexpr T a32 = static_cast<T>(0.75);
        constexpr T b1 = T{2} / T{9}, b2 = T{1} / T{3}, b3 = T{4} / T{9};
        constexpr T e1 = T{7} / T{24}, e2 = static_cast<T>(0.25), e3 = T{1} / T{3}, e4 = static_cast<T>(0.125);

        const ColVec<NX, T> k1 = f(t, x);
        const ColVec<NX, T> k2 = f(t + c2 * h, x + h * a21 * k1);
        const ColVec<NX, T> k3 = f(t + c3 * h, x + h * a32 * k2);

        const ColVec<NX, T> x3 = x + h * (b1 * k1 + b2 * k2 + b3 * k3);
        const ColVec<NX, T> k4 = f(t + h, x3);

        // 2nd order solution (for error estimate)
        const ColVec<NX, T> x2 = x + h * (e1 * k1 + e2 * k2 + e3 * k3 + e4 * k4);

        return IntegrationResult<NX, T>{x3, x3 - x2};
    }
};

/**
 * @brief Classical 3rd-order Runge-Kutta (RK3) integrator
 *
 * Fixed-step, explicit, 3rd-order accurate. Cheaper than RK4 (3 stages vs
 * 4) when 3rd-order accuracy is sufficient.
 *
 * @tparam NX Number of states
 * @tparam T  Scalar type
 */
template<size_t NX, typename T = double>
struct RK3 {
    /**
     * @brief Evolve LTI system state using RK3 (wraps the ODE form below)
     *
     * @param A     System matrix
     * @param B     Input matrix
     * @param x     Current state
     * @param u     Input (held constant over the duration of the step)
     * @param h     [s] Step size
     *
     * @return      Next state and estimated error (always 0)
     */
    template<size_t NU>
    constexpr IntegrationResult<NX, T> evolve(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B, const ColVec<NX, T>& x, const ColVec<NU, T>& u, T h) const {
        auto f = [&](double, const ColVec<NX, T>& xk) { return A * xk + B * u; };
        return evolve(f, x, 0.0, h);
    }

    /**
     * @brief Evolve nonlinear system state using classical RK3
     *
     *     k1 = f(t, x)
     *     k2 = f(t+h/2, x+h/2 k1)
     *     k3 = f(t+h, x−h k1+2h k2)
     *     x_{n+1} = x + (h/6)(k1 + 4k2 + k3)
     *
     * @param f     Right-hand side function: dx/dt = f(t, x)
     * @param x     Current state
     * @param t     [s] Current time (for passing to f)
     * @param h     [s] Step size
     *
     * @return      Next state and estimated error (always 0)
     */
    template<typename F>
    constexpr IntegrationResult<NX, T> evolve(const F& f, const ColVec<NX, T>& x, T t, T h) const {
        // Classical 3rd order Runge-Kutta
        constexpr double c2 = 0.5;
        constexpr double a21 = 0.5;
        constexpr double a31 = -1.0, a32 = 2.0;
        constexpr double b1 = 1.0 / 6.0, b2 = 4.0 / 6.0, b3 = 1.0 / 6.0;

        k1 = f(t, x);
        k2 = f(t + c2 * h, x + h * a21 * k1);
        k3 = f(t + h, x + h * (a31 * k1 + a32 * k2));
        return {x + h * (b1 * k1 + b2 * k2 + b3 * k3), {}};
    }

private:
    mutable ColVec<NX, T> k1, k2, k3;
};

} // namespace damp::sim
