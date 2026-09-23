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
#include <limits>

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
 * only accurate when ‖A‖��Ts ≪ 1. Does not preserve stability for stiff systems.
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

/**
 * @brief φ₁(Z) = ∫₀¹ e^{Z s} ds by scaling and squaring
 *
 * @return nullopt if the scaled Taylor series does not converge
 */
template<size_t NX, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>> phi1_expm(const Matrix<NX, NX, T>& Z) {
    using real_t = scalar_type_t<T>;
    const Matrix<NX, NX, T> I = Matrix<NX, NX, T>::identity();
    const real_t            norm = mat::infinity_norm(Z);
    size_t                  s = 0;
    real_t                  scaled = norm;
    constexpr real_t        theta = real_t(0.5);
    while (scaled > theta) {
        scaled *= real_t(0.5);
        ++s;
        if (s > 16) {
            return damp::nullopt;
        }
    }

    T scale = T{1};
    for (size_t i = 0; i < s; ++i) {
        scale *= T{0.5};
    }
    const Matrix<NX, NX, T> Zs = Z * scale;

    Matrix<NX, NX, T> E = I;
    Matrix<NX, NX, T> P = I;
    Matrix<NX, NX, T> term_e = I;
    const T           eps = static_cast<T>(std::numeric_limits<real_t>::epsilon()) * T{32};
    bool              converged = (norm == real_t{0});
    for (size_t k = 1; k <= 24 && !converged; ++k) {
        term_e = (term_e * Zs) * (T{1} / static_cast<T>(k));
        const Matrix<NX, NX, T> term_p = term_e * (T{1} / static_cast<T>(k + 1));
        E = E + term_e;
        P = P + term_p;
        if (mat::infinity_norm(term_p) <= eps) {
            converged = true;
        }
    }
    if (!converged) {
        return damp::nullopt;
    }

    for (size_t i = 0; i < s; ++i) {
        P = (P * (I + E)) * T{0.5};
        E = E * E;
    }
    return P;
}

/**
 * @brief ZOH integral map M ↦ ∫₀^{Ts} e^{Aτ} M dτ
 *
 * Equals A⁻¹(e^{A Ts} − I)M when A is invertible. A singular A (an integrator)
 * uses φ₁(A Ts) with scaling and squaring — the same family as @ref mat::expm,
 * not a fixed-length Taylor polynomial. Shared by B_d and G_d.
 *
 * @return nullopt if the φ₁ series does not converge
 */
template<size_t NX, size_t NC, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NC, T>> zoh_integrate_input(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NC, T>& M,
    const Matrix<NX, NX, T>& exp_A_Ts,
    T                        sampling_time
) {
    if constexpr (NC == 0) {
        return Matrix<NX, NC, T>{};
    }

    const Matrix I = Matrix<NX, NX, T>::identity();
    const Matrix rhs = (exp_A_Ts - I) * M;
    const auto   X_opt = mat::lu_solve(A, rhs);

    if (X_opt) {
        return X_opt.value();
    }

    const auto phi = phi1_expm(A * sampling_time);
    if (!phi) {
        return damp::nullopt;
    }
    return ((*phi) * sampling_time) * M;
}

/**
 * @brief Discretize using Zero-Order Hold (ZOH)
 *
 *     A_d = e^(ATs)
 *     B_d = ∫₀^{Ts} e^{Aτ} B dτ = A⁻¹(e^(ATs) − I)B   [φ₁ if A singular]
 *     G_d = ∫₀^{Ts} e^{Aτ} G dτ   (same map as B; process-noise input)
 *
 * Exact discretization assuming piecewise-constant input between samples.
 * Uses matrix exponential for A_d and LU solve for B_d/G_d (avoids forming A⁻¹).
 *
 * @note Compare with MATLAB®'s c2d(sys, Ts, 'zoh').
 * @see "Feedback Control of Dynamic Systems" (Franklin et al., 2015), §8.3
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr damp::optional<StateSpace<NX, NU, NY, T, NW, NV>> discretize_zoh_impl(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        sampling_time
) {
    // Compute A_d = exp(A * Ts)
    const Matrix A_scaled = sys.A * sampling_time;
    const Matrix exp_A_Ts = mat::expm(A_scaled);

    const auto B_d = zoh_integrate_input(sys.A, sys.B, exp_A_Ts, sampling_time);
    const auto G_d = zoh_integrate_input(sys.A, sys.G, exp_A_Ts, sampling_time);
    if (!B_d || !G_d) {
        return damp::nullopt;
    }

    //! C_d = C, D_d = D (output equation is the same); H_d = H (direct)
    const Matrix C_d = sys.C;
    const Matrix D_d = sys.D;
    const Matrix H_d = sys.H;

    return StateSpace{exp_A_Ts, B_d.value(), C_d, D_d, G_d.value(), H_d, sampling_time};
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