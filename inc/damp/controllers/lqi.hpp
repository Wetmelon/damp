// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lqi.hpp
 * @brief Linear-Quadratic-Integral design and runtime controller
 *
 * @code
 * using namespace damp;
 * // Discrete plant with Ts already set on sys
 * constexpr auto res = design::discrete_lqi(
 *     sys,
 *     Matrix<NX + NY, NX + NY>::identity(), // Q on [x; xi]
 *     Matrix<NU, NU>::identity());          // R
 * static_assert(res.success);
 * LQI ctrl(res);
 * auto u = ctrl.control(r, y, x); // advances internal xi
 * @endcode
 */

#include <cstddef>

#include "damp/design/riccati.hpp"
#include "damp/design/stability.hpp"
#include "damp/math/complex.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace design {

/**
 * @struct LQIResult
 * @brief LQI design result
 *
 * Optimal gain on the augmented state [x; xi] for discrete LQI tracking,
 * with Riccati solution S and closed-loop poles of the augmented system.
 * Use .as<float>() for deploy.
 *
 * @see discrete_lqi()
 * @see "Optimal Control" (Anderson & Moore, 1990), §4 (LQR on the augmented plant)
 */
template<size_t NX, size_t NU, size_t NY, typename T = double>
struct LQIResult {
    Matrix<NU, NX + NY, T>            K{};            ///< Optimal gain: u = -K*[x; xi]
    Matrix<NX + NY, NX + NY, T>       S{};            ///< Riccati equation solution
    ColVec<NX + NY, damp::complex<T>> e{};            ///< Closed-loop poles
    bool                              success{false}; ///< true if DARE converged

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return LQIResult<NX, NU, NY, U>{
            K.template as<U>(),
            S.template as<U>(),
            e.template as<damp::complex<U>>(),
            success
        };
    }
};

/**
 * @brief Discrete Linear-Quadratic-Integral (LQI) design for output tracking
 *
 * Augments the discrete plant with an output-error integrator and solves DARE
 * LQR on the augmented pair (A_aug, B_aug):
 *
 *     A_aug = |  A   0 |      B_aug = | B |
 *             | −C   I |              | 0 |
 *
 *     xi[k+1] = xi[k] + (r[k] − y[k])
 *     u[k]    = −K [x[k]; xi[k]]
 *
 * Q is the cost on [x; xi] (size (NX+NY)×(NX+NY)); R is the input cost.
 * Constant references and constant output disturbances yield zero steady-state
 * tracking error when the design succeeds and the loop is stable.
 *
 * @note Compare with MATLAB®'s lqi(sys, Q, R) — exposed as the lqi() alias in matlab.hpp.
 *
 * @see dare() for the Riccati solve
 * @see discrete_lqr() for pure regulation without integral action
 * @see "Optimal Control" (Anderson & Moore, 1990), §4
 *
 * @param sys  Discrete-time state-space plant (uses A, B, C)
 * @param Q    Augmented state cost ((NX+NY)×(NX+NY), positive semidefinite)
 * @param R    Input cost (NU×NU, positive definite)
 * @return LQIResult with K, S, closed-loop poles, and success
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr LQIResult<NX, NU, NY, T> discrete_lqi(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NX + NY, NX + NY, T>&       Q,
    const Matrix<NU, NU, T>&                 R
) {
    // Build augmented system
    Matrix<NX + NY, NX + NY, T> A_aug{};
    Matrix<NX + NY, NU, T>      B_aug{};

    /* Canonical form for LQI augmentation:
     *
     * A_aug = |  A      0 |
     *         | -C      I |
     *
     * B_aug = | B |
     *         | 0 |
     */

    // Top-left block: original A
    A_aug.template block<NX, NX>(0, 0) = sys.A;

    // Top-right block: C matrix for integral action
    A_aug.template block<NY, NX>(NX, 0) = -sys.C;

    // Top block of B_aug: original B
    B_aug.template block<NX, NU>(0, 0) = sys.B;

    // Bottom-right block: integrator memory
    A_aug.template block<NY, NY>(NX, NX) = Matrix<NY, NY, T>::identity();

    // Solve DARE for augmented system
    const auto dare_opt = dare(A_aug, B_aug, Q, R);
    if (!dare_opt) {
        return LQIResult<NX, NU, NY, T>{};
    }
    Matrix<NX + NY, NX + NY, T> P_aug = dare_opt.value();

    const auto K_opt = lqr_gain(A_aug, B_aug, P_aug, R);
    if (!K_opt) {
        return LQIResult<NX, NU, NY, T>{};
    }
    Matrix<NU, NX + NY, T> K_aug = K_opt.value();

    ColVec<NX + NY, damp::complex<T>> poles = stability::closed_loop_poles(A_aug, B_aug, K_aug);
    return LQIResult<NX, NU, NY, T>{K_aug, P_aug, poles, true};
}
} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Linear-Quadratic-Integral (LQI) controller
 *
 * Output tracking controller with integral action: u = -K * [x; xi]
 * where xi integrates the output error (r - y).
 * Provides zero steady-state error for constant references and disturbances.
 *
 * @tparam NX Number of states
 * @tparam NU Number of control inputs
 * @tparam NY Number of outputs
 * @tparam T  Scalar type (default: float — embedded runtime)
 */
template<size_t NX, size_t NU, size_t NY, typename T = float>
struct LQI {
    Matrix<NU, NX + NY, T> K{};  ///< Optimal gain: u = -K*[x; xi]
    ColVec<NY, T>          xi{}; ///< Integral of tracking error, owned by the controller

    constexpr LQI() = default;
    constexpr LQI(const Matrix<NU, NX + NY, T>& K_) : K(K_) {} // NOLINT

    /// From design result (any scalar); converts K via @c .as\<T\>().
    template<typename U>
    constexpr LQI(const design::LQIResult<NX, NU, NY, U>& result) // NOLINT
        : K(result.K.template as<T>()) {}

    template<typename U>
    constexpr LQI(const LQI<NX, NU, NY, U>& other) : K(other.K.template as<T>()), xi(other.xi.template as<T>()) {} // NOLINT

    /**
     * @brief Compute control from a caller-supplied augmented state
     *
     * Computes u = -K * x_aug where x_aug = [x; xi]. Use this when you maintain
     * the integral state yourself; it does not touch the internal integrator.
     *
     * @param x_aug Augmented state vector [x; xi]
     * @return Control input vector u
     */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NX + NY, T>& x_aug) const {
        return ColVec<NU, T>(-K * x_aug);
    }

    /**
     * @brief Compute control with the controller's own integral action
     *
     * Uses the internal integrator state: u = -[Kx Ki]·[x; xi] with the current
     * xi, then advances xi[k+1] = xi[k] + (r - y) (matching the −C augmentation
     * sign convention). Gives zero steady-state error to constant references
     * without the caller tracking the integral.
     *
     * Reference-first argument order (reference, measurement, state), matching
     * the library-wide control(r, feedback...) convention.
     *
     * @param r Output reference
     * @param y Measured output
     * @param x Current (estimated or measured) plant state
     * @return Control input vector u
     */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NY, T>& r, const ColVec<NY, T>& y, const ColVec<NX, T>& x) {
        return control(r, y, x, ColVec<NX, T>{});
    }

    /**
     * @brief Compute control against a full state reference (trajectory tracking)
     *
     * u = -Kx·(x − x_ref) − Ki·xi, then xi[k+1] = xi[k] + (r − y). The state-feedback term is
     * the feedback half of a 2-DOF law (see LQR::control(const ColVec<NX,T>&, const
     * ColVec<NX,T>&) for the u_ff input feedforward it omits). Without @p x_ref the feedback
     * regulates the whole state toward the origin, so on a moving reference it opposes the
     * commanded motion and the integrator alone drags the plant along — a velocity-proportional
     * lag plus the overshoot that windup dumps at arrival. Supplying a feasible reference state
     * (e.g. the profile's position and velocity, consistent with C·x_ref = r) makes the feedback
     * act on trajectory error; the integrator then supplies the quasi-static feedforward, so
     * unlike the pure regulator there is no droop — but a fast-changing u_ff (the acceleration
     * phase of a move) still lags unless fed forward on the input.
     *
     * Reference-first argument order, matching the library-wide control(r, feedback...)
     * convention.
     *
     * @param r     Output reference (drives the integrator)
     * @param y     Measured output
     * @param x     Current (estimated or measured) plant state
     * @param x_ref State reference the feedback regulates toward (must be feasible)
     * @return Control input vector u
     */
    [[nodiscard]] constexpr ColVec<NU, T>
    control(const ColVec<NY, T>& r, const ColVec<NY, T>& y, const ColVec<NX, T>& x, const ColVec<NX, T>& x_ref) {
        const auto    Kx = K.template block<NU, NX>(0, 0);
        const auto    Ki = K.template block<NU, NY>(0, NX);
        ColVec<NU, T> u = ColVec<NU, T>(-((Kx * (x - x_ref)) + (Ki * xi)));
        xi = ColVec<NY, T>(xi + (r - y));
        return u;
    }

    /// Clear the integral state.
    constexpr void reset() { xi = ColVec<NY, T>{}; }

    [[nodiscard]] constexpr const Matrix<NU, NX + NY, T>& getK() const { return K; }
};

} // namespace damp