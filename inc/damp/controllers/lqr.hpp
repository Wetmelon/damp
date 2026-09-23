// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lqr.hpp
 * @brief Linear-Quadratic Regulator design and runtime controller
 *
 * @code
 * using namespace damp;
 * constexpr Matrix<2,2> A{{1.0, 0.01}, {0.0, 1.0}};
 * constexpr Matrix<2,1> B{{0.0}, {0.01}};
 * constexpr auto res = design::discrete_lqr(
 *     A, B, Matrix<2,2>::identity(), Matrix<1,1>{{0.1}});
 * static_assert(res.success);
 * LQR ctrl(res.as<float>());  // runtime defaults to float
 * auto u = ctrl.control(x);   // u = -K x
 * @endcode
 *
 * Short MATLAB®-style aliases (embeddable): design::dlqr, design::lqrd.
 * Design results default to @c double; runtime LQR / StateFeedback default to @c float.
 *
 * @see examples/control/cart_pole/ for Design Is Deploy (deploy.hpp + sketch + sil)
 * @see examples/control/pendulum/ for a 2-state pendulum design-is-deploy example
 */

#include <cstddef>

#include "damp/design/riccati.hpp"
#include "damp/design/stability.hpp"
#include "damp/math/complex.hpp"
#include "damp/matrix/block.hpp"
#include "damp/matrix/functions.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace design {

/**
 * @struct LQRResult
 * @brief Linear-Quadratic Regulator design result
 *
 * Contains the optimal gain K, Riccati solution S, and closed-loop poles.
 * Use .as<float>() to convert for embedded deployment.
 *
 * @note Compare with MATLAB®'s [K,S,P] = dlqr(...) output structure.
 * @see "Optimal Control" (Anderson & Moore, 1990), §4.3
 */
template<size_t NX, size_t NU, typename T = double>
struct LQRResult {
    Matrix<NU, NX, T>            K{};            ///< Optimal gain: u = −Kx
    Matrix<NX, NX, T>            S{};            ///< DARE solution (positive semidefinite)
    ColVec<NX, damp::complex<T>> e{};            ///< Closed-loop poles (eigenvalues of A − BK)
    bool                         success{false}; ///< true if DARE converged
    bool                         discrete{true}; ///< false: @c e are s-plane poles (continuous_lqr)

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return LQRResult<NX, NU, U>{
            K.template as<U>(),
            S.template as<U>(),
            e.template as<damp::complex<U>>(),
            success,
            discrete
        };
    }

    /**
     * @brief Check if the closed-loop system is stable
     *
     * Stability is determined by checking if all closed-loop poles lie within the unit circle.
     *
     * @return true if stable, false otherwise
     */
    [[nodiscard]] constexpr bool is_stable() const {
        for (size_t i = 0; i < NX; ++i) {
            if (discrete) {
                if (e[i].abs() >= static_cast<T>(1.0)) {
                    return false;
                }
            } else if (e[i].real() >= static_cast<T>(0)) {
                return false;
            }
        }
        return true;
    }
};

/**
 * @brief Discrete-time Linear-Quadratic Regulator design
 *
 * Computes the optimal state-feedback gain @f$ K @f$ that minimizes
 * @f[
 *   J = \sum_k \bigl( x_k^\top Q x_k + u_k^\top R u_k + 2 x_k^\top N u_k \bigr)
 * @f]
 * subject to the discrete-time dynamics @f$ x_{k+1} = A x_k + B u_k @f$.
 *
 * The gain is applied as @f$ u = -Kx @f$. The solution is found via the Discrete
 * Algebraic Riccati Equation (DARE).
 *
 * @note Compare with MATLAB®'s dlqr(A, B, Q, R, N).
 *
 * @see dare() for the underlying Riccati solver
 * @see discrete_lqr_from_continuous() to design from a continuous-time system
 * @see "Optimal Control" (Anderson & Moore, 1990), Chapter 4
 *
 * @param A  State transition matrix (NX × NX)
 * @param B  Control input matrix (NX × NU)
 * @param Q  State cost matrix (NX × NX, positive semidefinite)
 * @param R  Input cost matrix (NU × NU, positive definite)
 * @param N  Cross-term cost matrix (NX × NU, default: zero)
 * @return LQRResult with gain K, Riccati solution S, and closed-loop poles
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr LQRResult<NX, NU, T> discrete_lqr(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N = Matrix<NX, NU, T>{}
) {
    const auto dare_opt = dare(A, B, Q, R, N);
    if (!dare_opt) {
        return LQRResult<NX, NU, T>{};
    }
    const Matrix<NX, NX, T> S = dare_opt.value();

    const auto K_opt = lqr_gain(A, B, S, R, N);
    if (!K_opt) {
        return LQRResult<NX, NU, T>{};
    }

    const Matrix<NU, NX, T> K = K_opt.value();
    return LQRResult<NX, NU, T>{K, S, stability::closed_loop_poles(A, B, K), true};
}

/**
 * @struct LQRCost
 * @brief Discretized LQR cost weights (Q, R, N) for a sampled-data problem
 */
template<size_t NX, size_t NU, typename T = double>
struct LQRCost {
    Matrix<NX, NX, T> Q{}; ///< Discrete state cost
    Matrix<NU, NU, T> R{}; ///< Discrete input cost
    Matrix<NX, NU, T> N{}; ///< Discrete cross-term cost
};

/**
 * @brief Discretize a continuous LQR cost integral over one sample (Van Loan)
 *
 * Maps the continuous running cost
 * @f[
 *   J = \int_0^\infty \big( x^\top Q x + 2 x^\top N u + u^\top R u \big)\, dt
 * @f]
 * to its exact discrete equivalent @f$ \sum (x^\top Q_d x + 2 x^\top N_d u + u^\top R_d u) @f$
 * for a zero-order-hold input. Naively reusing the continuous @f$ (Q,R,N) @f$ on
 * the discretized dynamics is only first-order accurate in @f$ T_s @f$; this is exact.
 *
 * Augmenting the held input as constant states @f$ \bar A = [A\;B;\,0\;0] @f$ with
 * weight @f$ \bar Q = [Q\;N;\,N^\top\;R] @f$, the discrete weights follow from a
 * single matrix exponential (Van Loan, 1978):
 * @f[
 *   \exp\!\left( \begin{bmatrix} -\bar A^\top & \bar Q \\ 0 & \bar A \end{bmatrix} T_s \right)
 *     = \begin{bmatrix} M_{11} & M_{12} \\ 0 & M_{22} \end{bmatrix}, \quad
 *   \begin{bmatrix} Q_d & N_d \\ N_d^\top & R_d \end{bmatrix} = M_{22}^\top M_{12}.
 * @f]
 *
 * @note This is the cost-discretization step of MATLAB®'s lqrd(A, B, Q, R, Ts).
 * @see "Computing Integrals Involving the Matrix Exponential" (Van Loan, 1978),
 *      IEEE TAC 23(3), https://doi.org/10.1109/TAC.1978.1101743
 *
 * @param A   State transition matrix (continuous-time, NX × NX)
 * @param B   Control input matrix (continuous-time, NX × NU)
 * @param Q   Continuous state cost (NX × NX, positive semidefinite)
 * @param R   Continuous input cost (NU × NU, positive definite)
 * @param Ts  Sampling time [s]
 * @param N   Continuous cross-term cost (NX × NU, default: zero)
 * @return Discrete cost weights {Q_d, R_d, N_d}
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr LQRCost<NX, NU, T> discretize_lqr_cost(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    T                        Ts,
    const Matrix<NX, NU, T>& N = Matrix<NX, NU, T>{}
) {
    constexpr size_t NA = NX + NU; // augmented (state + held input) dimension

    //! Augmented dynamics Ā = [A B; 0 0] and weight Q̄ = [Q N; Nᵀ R]
    const Matrix<NU, NX, T> Nt = N.t();
    Matrix<NA, NA, T>       Abar{};
    Abar.template block<NX, NX>(0, 0) = A;
    Abar.template block<NX, NU>(0, NX) = B;

    Matrix<NA, NA, T> Qbar{};
    Qbar.template block<NX, NX>(0, 0) = Q;
    Qbar.template block<NX, NU>(0, NX) = N;
    Qbar.template block<NU, NX>(NX, 0) = Nt;
    Qbar.template block<NU, NU>(NX, NX) = R;

    //! Van Loan block Z = [-Āᵀ Q̄; 0 Ā]·Ts, then exp(Z) = [M11 M12; 0 M22]
    const Matrix<NA, NA, T>   AbarT = Abar.t();
    const Matrix<NA, NA, T>   neg_AbarT = AbarT * (-Ts);
    const Matrix<NA, NA, T>   Qbar_Ts = Qbar * Ts;
    const Matrix<NA, NA, T>   Abar_Ts = Abar * Ts;
    Matrix<2 * NA, 2 * NA, T> Z{};
    Z.template block<NA, NA>(0, 0) = neg_AbarT;
    Z.template block<NA, NA>(0, NA) = Qbar_Ts;
    Z.template block<NA, NA>(NA, NA) = Abar_Ts;

    const Matrix<2 * NA, 2 * NA, T> G = mat::expm(Z);
    const Matrix<NA, NA, T>         M12 = G.template block<NA, NA>(0, NA);
    const Matrix<NA, NA, T>         M22 = G.template block<NA, NA>(NA, NA);
    const Matrix<NA, NA, T>         M22t = M22.t();

    //! Q̄_d = M22ᵀ·M12; symmetrize to scrub round-off asymmetry before slicing
    Matrix<NA, NA, T>       Qbar_d = M22t * M12;
    const Matrix<NA, NA, T> Qbar_dt = Qbar_d.t();
    Qbar_d = (Qbar_d + Qbar_dt) * static_cast<T>(0.5);

    return LQRCost<NX, NU, T>{
        Qbar_d.template block<NX, NX>(0, 0),
        Qbar_d.template block<NU, NU>(NX, NX),
        Qbar_d.template block<NX, NU>(0, NX)
    };
}

/**
 * @brief Design discrete LQR from continuous-time system via discretization
 *
 * Discretizes both the dynamics (ZOH) and the cost integral (Van Loan), then
 * solves the discrete LQR problem. Equivalent to MATLAB®'s lqrd — discretizing
 * the cost is what distinguishes this from feeding continuous Q, R, N into the
 * sampled dynamics, which is only first-order accurate in Ts.
 *
 * @note Compare with MATLAB®'s lqrd(A, B, Q, R, Ts).
 *
 * @see discrete_lqr() for the discrete-time design
 * @see discretize() for the ZOH dynamics discretization step
 * @see discretize_lqr_cost() for the Van Loan cost discretization step
 *
 * @param A   State transition matrix (continuous-time, NX × NX)
 * @param B   Control input matrix (continuous-time, NX × NU)
 * @param Q   State cost matrix (NX × NX, positive semidefinite)
 * @param R   Input cost matrix (NU × NU, positive definite)
 * @param Ts  Sampling time [s]
 * @param N   Cross-term cost matrix (NX × NU, default: zero)
 * @return LQRResult with gain K, Riccati solution S, and closed-loop poles
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr LQRResult<NX, NU, T> discrete_lqr_from_continuous(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    T                        Ts,
    const Matrix<NX, NU, T>& N = Matrix<NX, NU, T>{}
) {
    StateSpace<NX, NU, NX, T, NX, NX> sys_c{A, B, Matrix<NX, NX, T>::identity()};
    const auto                        sys_d = *discretize(sys_c, Ts, DiscretizationMethod::ZOH);
    const auto                        cost = discretize_lqr_cost(A, B, Q, R, Ts, N);
    return discrete_lqr(sys_d.A, sys_d.B, cost.Q, cost.R, cost.N);
}

/**
 * @brief LQR from continuous state-space system via discretization
 *
 * @param sys State-space system (continuous-time)
 * @param Q   State cost matrix (NX × NX, positive semidefinite)
 * @param R   Input cost matrix (NU × NU, positive definite)
 * @param Ts  Sampling time [s]
 * @param N   Cross-term cost matrix (NX × NU, default: zero)
 * @return LQRResult with gain K, Riccati solution S, and closed-loop poles
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr LQRResult<NX, NU, T> discrete_lqr_from_continuous(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NX, NX, T>&                 Q,
    const Matrix<NU, NU, T>&                 R,
    T                                        Ts,
    const Matrix<NX, NU, T>&                 N = Matrix<NX, NU, T>{}
) {
    return discrete_lqr_from_continuous(sys.A, sys.B, Q, R, Ts, N);
}

/**
 * @brief Discrete-time LQR design (MATLAB®-style short name)
 *
 * Thin alias of @ref discrete_lqr. Prefer the descriptive name in new prose;
 * use this when matching MATLAB® / textbooks that write `dlqr`.
 *
 * @note Compare with MATLAB®'s dlqr(A, B, Q, R, N).
 * @see discrete_lqr()
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr LQRResult<NX, NU, T> dlqr(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N = Matrix<NX, NU, T>{}
) {
    return discrete_lqr(A, B, Q, R, N);
}

/**
 * @brief Sampled-data LQR from continuous plant (MATLAB®-style short name)
 *
 * Thin alias of @ref discrete_lqr_from_continuous (ZOH dynamics + Van Loan cost).
 *
 * @note Compare with MATLAB®'s lqrd(A, B, Q, R, Ts).
 * @see discrete_lqr_from_continuous()
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr LQRResult<NX, NU, T> lqrd(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    T                        Ts,
    const Matrix<NX, NU, T>& N = Matrix<NX, NU, T>{}
) {
    return discrete_lqr_from_continuous(A, B, Q, R, Ts, N);
}

/**
 * @brief Sampled-data LQR from continuous StateSpace (MATLAB®-style short name)
 *
 * Thin alias of @ref discrete_lqr_from_continuous.
 *
 * @note Compare with MATLAB®'s lqrd(sys, Q, R, Ts).
 * @see discrete_lqr_from_continuous()
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr LQRResult<NX, NU, T> lqrd(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NX, NX, T>&                 Q,
    const Matrix<NU, NU, T>&                 R,
    T                                        Ts,
    const Matrix<NX, NU, T>&                 N = Matrix<NX, NU, T>{}
) {
    return discrete_lqr_from_continuous(sys, Q, R, Ts, N);
}

/**
 * @brief Continuous-time Linear-Quadratic Regulator design
 *
 * Computes the optimal state-feedback gain @f$ K @f$ that minimizes
 * @f[
 *   J = \int_0^\infty \bigl( x^\top Q x + u^\top R u + 2 x^\top N u \bigr)\, dt
 * @f]
 * subject to the continuous-time dynamics @f$ \dot x = A x + B u @f$, applied as
 * @f$ u = -Kx @f$. Solves the Continuous Algebraic Riccati Equation (CARE) for
 * the stabilizing @f$ S @f$, then forms the gain via
 * @f$ R K^\top = B^\top S + N^\top @f$ (solve, not explicit @f$ R^{-1} @f$).
 *
 * @note Compare with MATLAB®'s [K,S,e] = lqr(A, B, Q, R, N).
 * @note Unlike @ref discrete_lqr, the returned `e` are continuous-time
 *       (s-plane) closed-loop poles of @f$ (A - BK) @f$ and @c discrete is
 *       false, so @c is_stable tests @f$ \mathrm{Re}(e) < 0 @f$.
 *
 * @see care() for the underlying Riccati solver
 * @see discrete_lqr_from_continuous() for the sampled-data (digital) design
 * @see "Optimal Control" (Anderson & Moore, 1990), §3.3
 *
 * @param A  State matrix (NX × NX)
 * @param B  Control input matrix (NX × NU)
 * @param Q  State cost matrix (NX × NX, positive semidefinite)
 * @param R  Input cost matrix (NU × NU, positive definite)
 * @param N  Cross-term cost matrix (NX × NU, default: zero)
 * @return LQRResult with gain K, CARE solution S, and continuous-time poles
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr LQRResult<NX, NU, T> continuous_lqr(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N = Matrix<NX, NU, T>{}
) {
    const auto care_opt = care(A, B, Q, R, N);
    if (!care_opt) {
        return LQRResult<NX, NU, T>{};
    }
    const Matrix<NX, NX, T> S = care_opt.value();

    // Continuous-time optimal gain K = R⁻¹(BᵀS + Nᵀ), via a decomposition solve
    // of R·K = BᵀS + Nᵀ rather than forming R⁻¹ explicitly (better conditioned
    // for NU > 1; matches care_schur's lu_solve usage).
    const auto K_opt = mat::lu_solve(R, (B.transpose() * S) + N.transpose());
    if (!K_opt) {
        return LQRResult<NX, NU, T>{};
    }
    const Matrix<NU, NX, T> K = K_opt.value();
    return LQRResult<NX, NU, T>{K, S, stability::closed_loop_poles(A, B, K), true, false};
}

/// @brief Continuous LQR from a continuous-time StateSpace (uses A, B).
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr LQRResult<NX, NU, T> continuous_lqr(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NX, NX, T>&                 Q,
    const Matrix<NU, NU, T>&                 R,
    const Matrix<NX, NU, T>&                 N = Matrix<NX, NU, T>{}
) {
    return continuous_lqr(sys.A, sys.B, Q, R, N);
}

} // namespace design

namespace detail {

/// Failed design:: result passed to a runtime law. Not constexpr, so a
/// constant-evaluated constructor becomes ill-formed; at runtime it traps.
[[noreturn]] inline void refuse_failed_design() noexcept {
#if defined(__GNUC__) || defined(__clang__)
    __builtin_trap();
#elif defined(_MSC_VER)
    __fastfail(7);
#else
    for (;;) {
    }
#endif
}

} // namespace detail

/**
 * @ingroup discrete_controllers
 * @brief Runtime full-state feedback law u = −Kx
 *
 * Stores only the gain matrix K. One matrix-vector multiply per call — suitable for
 * ISR. Gains may come from pole placement (design::place), LQR
 * (@ref design::discrete_lqr / design::dlqr), or any other design that produces K.
 *
 * LQR is an alias of this type (MATLAB®-shaped name for the same runtime).
 *
 * @see design::place(), design::discrete_lqr(), design::dlqr()
 * @see "Optimal Control" (Anderson & Moore, 1990), §4.1
 *
 * @tparam NX Number of states
 * @tparam NU Number of control inputs
 * @tparam T  Scalar type (default: float — embedded runtime)
 */
template<size_t NX, size_t NU, typename T = float>
struct StateFeedback {
    Matrix<NU, NX, T> K{};

    constexpr StateFeedback() = default;
    constexpr explicit StateFeedback(const Matrix<NU, NX, T>& K_) : K(K_) {}

    /// From a successful design result (any scalar); converts K via @c .as\<T\>().
    ///
    /// A failed design must not become a law. Constant evaluation of a failed
    /// result is ill-formed. A runtime failed result does not return.
    template<typename U>
    constexpr StateFeedback(const design::LQRResult<NX, NU, U>& result) // NOLINT
        : K{} {
        if (result.success) {
            K = result.K.template as<T>();
        } else {
            detail::refuse_failed_design();
        }
    }

    template<typename U>
    constexpr StateFeedback(const StateFeedback<NX, NU, U>& other) : K(other.getK().template as<T>()) {} // NOLINT

    /**
     * @brief Compute regulator control law
     *
     * Drives state to zero: u = -K*x
     *
     * @param x Current state vector
     * @return Control input vector u
     */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NX, T>& x) const {
        return ColVec<NU, T>(-K * x);
    }

    /**
     * @brief Servo control law: state-error feedback u = -K*(x - x_ref)
     *
     * The feedback half of a 2-DOF tracking law. The complete law adds the feedforward input
     * u_ff that makes x_ref a feasible trajectory (x_ref[k+1] = A*x_ref[k] + B*u_ff[k]):
     * u = u_ff - K*(x - x_ref). With no integrator to supply it, omitting u_ff leaves a
     * steady-state offset (droop) — the classic reason the regulator-tracking form is paired
     * with the Nu*r input feedforward (the Nx/Nu precompensator). @p x_ref must be a feasible
     * state, not an arbitrary target.
     *
     * Reference-first argument order, matching the library-wide control(r, feedback...)
     * convention.
     *
     * @param x_ref Reference state vector
     * @param x     Current state vector
     * @return Control input vector u
     */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NX, T>& x_ref, const ColVec<NX, T>& x) const {
        return ColVec<NU, T>(-K * (x - x_ref));
    }

    /// No internal state; provided so StateFeedback models the StateFeedbackController concept.
    constexpr void reset() {}

    [[nodiscard]] constexpr const Matrix<NU, NX, T>& getK() const { return K; }
};

/**
 * @brief Alias for StateFeedback — MATLAB®-shaped name for the same runtime.
 *
 * Prefer StateFeedback when the gain came from pole placement or another
 * non-LQR design; LQR remains the familiar name when the gain is from
 * @ref design::discrete_lqr.
 */
template<size_t NX, size_t NU, typename T = float>
using LQR = StateFeedback<NX, NU, T>;

} // namespace damp
