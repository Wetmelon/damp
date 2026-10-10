// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file discretization.hpp
 * @brief Continuous-to-discrete state-space discretization
 */

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/matrix/matrix.hpp"
#include "state_space.hpp"

namespace damp {

/**
 * @brief Discretization methods for continuous-time state-space systems
 *
 * @see "Feedback Control of Dynamic Systems" (Franklin et al., 2015), Chapter 8
 */
enum class DiscretizationMethod {
    ForwardEuler, ///< Explicit Euler: Ad = I + ATs (first-order, simple, low overhead)
    ZOH,          ///< Zero-Order Hold: Ad = e^(ATs), exact for piecewise-constant inputs
    Tustin,       ///< Bilinear transform: s → (2/Ts)(z−1)/(z+1), preserves stability
};

namespace detail {

/**
 * @brief Discretize using Forward Euler (explicit Euler)
 *
 *     A_d = I + ATs,  B_d = BTs,  C_d = C,  D_d = D
 *
 * First-order approximation of the matrix exponential. Simple and fast but
 * only accurate when ‖A‖·Ts ≪ 1. Does not preserve stability for stiff systems.
 *
 * @see "Feedback Control of Dynamic Systems" (Franklin et al., 2015), §8.3
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr StateSpace<NX, NU, NY, T, NW, NV> discretize_forward_euler_impl(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        sampling_time
) {
    // A_d = I + A·Ts, B_d = B·Ts; C/D unchanged. Scale G like B when NW > 0.
    return StateSpace<NX, NU, NY, T, NW, NV>{
        .A = Matrix<NX, NX, T>::identity() + (sys.A * sampling_time),
        .B = sys.B * sampling_time,
        .C = sys.C,
        .D = sys.D,
        .G = sys.G * sampling_time,
        .H = sys.H,
        .Ts = sampling_time,
    };
}

/// Copy the NX × NC block of M starting at column c0 (rows 0 … NX−1).
template<size_t NX, size_t NC, size_t N, typename T>
[[nodiscard]] constexpr Matrix<NX, NC, T> leading_rows_block(const Matrix<N, N, T>& M, size_t c0) {
    Matrix<NX, NC, T> out{};
    for (size_t i = 0; i < NX; ++i) {
        for (size_t j = 0; j < NC; ++j) {
            out(i, j) = M(i, c0 + j);
        }
    }
    return out;
}

/**
 * @brief Discretize using Zero-Order Hold (ZOH)
 *
 *     A_d = e^(ATs)
 *     B_d = ∫₀^{Ts} e^{Aτ} B dτ
 *     G_d = ∫₀^{Ts} e^{Aτ} G dτ   (same map as B; process-noise input)
 *
 * Exact discretization assuming piecewise-constant input between samples. All
 * three come from one exponential of the augmented matrix (Van Loan, 1978):
 * @f[
 *   \exp\left(\begin{bmatrix} A & B & G \\ 0 & 0 & 0 \end{bmatrix} T_s\right)
 *   = \begin{bmatrix} A_d & B_d & G_d \\ 0 & I & 0 \\ 0 & 0 & I \end{bmatrix}.
 * @f]
 * Unlike @f$ A^{-1}(e^{AT_s} - I)B @f$ this has no cancellation when @f$ A @f$ is
 * singular or nearly so (integrators, very slow modes).
 *
 * @note Compare with MATLAB®'s c2d(sys, Ts, 'zoh').
 * @see "Computing integrals involving the matrix exponential" (Van Loan, 1978)
 * @see "Feedback Control of Dynamic Systems" (Franklin et al., 2015), §8.3
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr StateSpace<NX, NU, NY, T, NW, NV> discretize_zoh_impl(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        sampling_time
) {
    constexpr size_t N = NX + NU + NW;
    Matrix<N, N, T>  M{};
    for (size_t i = 0; i < NX; ++i) {
        for (size_t j = 0; j < NX; ++j) {
            M(i, j) = sys.A(i, j) * sampling_time;
        }
        for (size_t j = 0; j < NU; ++j) {
            M(i, NX + j) = sys.B(i, j) * sampling_time;
        }
        for (size_t j = 0; j < NW; ++j) {
            M(i, NX + NU + j) = sys.G(i, j) * sampling_time;
        }
    }
    const Matrix<N, N, T> E = mat::expm(M);

    return StateSpace<NX, NU, NY, T, NW, NV>{
        .A = leading_rows_block<NX, NX>(E, 0),
        .B = leading_rows_block<NX, NU>(E, NX),
        .C = sys.C,
        .D = sys.D,
        .G = leading_rows_block<NX, NW>(E, NX + NU),
        .H = sys.H,
        .Ts = sampling_time,
    };
}

/**
 * @brief Discretize a continuous-time state-space system using Tustin method
 *
 * Tustin (Bilinear Transform) implementation
 *    Maps: s → (2/Ts) · (z − 1) / (z + 1)
 *
 *    Let L = (I − A·Ts/2). All results are computed via LU solve
 *    against L rather than forming L⁻¹ explicitly:
 *      A_d:  solve L · A_d = (I + A·Ts/2)
 *      B_d:  solve L · X = B,  then B_d = Ts · X
 *      C_d:  solve Lᵀ · X = Cᵀ, then C_d = Xᵀ
 *      D_d = D + (Ts/2) · C_d · B
 *
 * @note Compare with MATLAB®'s c2d(sys, Ts, 'tustin').
 * @see "Feedback Control of Dynamic Systems" (Franklin et al., 2015), §8.6
 *
 * @param sys             Continuous-time state-space model
 * @param sampling_time   Desired sampling period for discrete system
 * @return Discretized state-space system
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr damp::optional<StateSpace<NX, NU, NY, T, NW, NV>> discretize_tustin_impl(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        sampling_time
) {
    const T      ts_half = sampling_time / T{2};
    const Matrix I = Matrix<NX, NX, T>::identity();
    const Matrix L = I - sys.A * ts_half;

    // A_d: solve L · A_d = (I + A·Ts/2)
    const auto A_d_opt = mat::lu_solve(L, I + sys.A * ts_half);
    if (!A_d_opt) {
        return damp::nullopt;
    }
    const Matrix A_d = A_d_opt.value();

    // B_d: solve L · X = B, then B_d = Ts · X
    const auto B_d_opt = mat::lu_solve(L, sys.B);
    if (!B_d_opt) {
        return damp::nullopt;
    }
    const Matrix B_d = B_d_opt.value() * sampling_time;

    // C_d = C · L⁻¹ → solve Lᵀ · X = Cᵀ, then C_d = Xᵀ
    const auto C_d_t_opt = mat::lu_solve(L.transpose(), sys.C.transpose());
    if (!C_d_t_opt) {
        return damp::nullopt;
    }
    const Matrix C_d = C_d_t_opt.value().transpose();

    // D_d = D + (Ts/2) · C_d · B
    const Matrix D_d = sys.D + C_d * sys.B * ts_half;

    // G_d: same map as B
    const auto G_d_opt = mat::lu_solve(L, sys.G);
    if (!G_d_opt) {
        return damp::nullopt;
    }
    const Matrix G_d = G_d_opt.value() * sampling_time;
    const Matrix H_d = sys.H;

    return StateSpace{A_d, B_d, C_d, D_d, G_d, H_d, sampling_time};
}

} // namespace detail

/**
 * @brief Discretize a continuous-time state-space system
 *
 * Converts @f$ \dot x = Ax + Bu @f$ to @f$ x[k+1] = A_d x[k] + B_d u[k] @f$.
 * Already-discrete models and non-positive @p sampling_time return @p sys unchanged.
 * Tustin returns @c nullopt if the bilinear factor @f$ I - A T_s/2 @f$ is singular
 * (no silent fallback to ZOH).
 *
 * @note Compare with MATLAB®'s c2d(sys, Ts, method).
 * @see "Feedback Control of Dynamic Systems" (Franklin et al., 2015), Chapter 8
 *
 * @param sys           Continuous-time state-space model (@c Ts = 0)
 * @param sampling_time Desired sampling period [s]
 * @param method        Discretization method (ForwardEuler, ZOH, or Tustin)
 * @return Discrete system, or @c nullopt if the method cannot be applied
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr damp::optional<StateSpace<NX, NU, NY, T, NW, NV>> discretize(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        sampling_time,
    DiscretizationMethod                     method = DiscretizationMethod::ZOH
) {
    if (sampling_time <= T{0}) {
        return sys;
    }

    if (sys.Ts > T{0}) {
        return sys;
    }

    switch (method) {
        case DiscretizationMethod::ForwardEuler:
            return detail::discretize_forward_euler_impl(sys, sampling_time);
        case DiscretizationMethod::ZOH:
            return detail::discretize_zoh_impl(sys, sampling_time);
        case DiscretizationMethod::Tustin:
            return detail::discretize_tustin_impl(sys, sampling_time);
        default:
            return detail::discretize_forward_euler_impl(sys, sampling_time);
    }
}
} // namespace damp