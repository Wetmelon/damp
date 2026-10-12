// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file stability.hpp
 * @brief Stability analysis helpers
 *
 * @defgroup stability_analysis Stability Analysis
 * @brief Controllability, observability, and related structural stability checks
 *
 * For continuous systems: stable if all eigenvalues have Re(λ) < 0 (left half plane)
 * For discrete systems: stable if all eigenvalues have |λ| < 1 (inside unit circle)
 *
 * @note The eigenvalue-based routines (is_stable_discrete, stability_margin_*,
 *       closed_loop_poles) use mat::compute_eigenvalues, which resolves
 *       fully-complex eigenvalues for any N (closed form for N ≤ 2, Francis QR
 *       beyond).
 */
#include <cstddef>

#include "damp/backend.hpp"
#include "damp/design/lyapunov.hpp"
#include "damp/math/complex.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/svd.hpp"

namespace damp {
namespace stability {

// ============================================================================
// Structural Analysis (Controllability / Observability)
// ============================================================================

/**
 * @brief Compute the controllability matrix [B, AB, A²B, ..., A^(N-1)B]
 *
 * The system (A, B) is controllable iff controllability_matrix(A, B) has full row rank.
 *
 * @tparam NX Number of states
 * @tparam NU Number of inputs
 * @tparam T  Scalar type
 * @param A   State matrix
 * @param B   Input matrix
 * @return Matrix<NX, NX*NU, T> controllability matrix
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr Matrix<NX, NX * NU, T>
controllability_matrix(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B) noexcept {
    Matrix<NX, NX * NU, T> Co{};
    Matrix<NX, NU, T>      AB = B;
    for (size_t k = 0; k < NX; ++k) {
        for (size_t r = 0; r < NX; ++r) {
            for (size_t c = 0; c < NU; ++c) {
                Co(r, (k * NU) + c) = AB(r, c);
            }
        }
        if (k + 1 < NX) {
            AB = A * AB;
        }
    }
    return Co;
}

/**
 * @brief Compute the observability matrix [C; CA; CA²; ...; CA^(N-1)]
 *
 * The system (A, C) is observable iff observability_matrix(A, C) has full column rank.
 *
 * @tparam NX Number of states
 * @tparam NY Number of outputs
 * @tparam T  Scalar type
 * @param A   State matrix
 * @param C   Output matrix
 * @return Matrix<NX*NY, NX, T> observability matrix
 */
template<size_t NX, size_t NY, typename T = double>
[[nodiscard]] constexpr Matrix<NX * NY, NX, T>
observability_matrix(const Matrix<NX, NX, T>& A, const Matrix<NY, NX, T>& C) noexcept {
    Matrix<NX * NY, NX, T> Ob{};
    Matrix<NY, NX, T>      CA = C;
    for (size_t k = 0; k < NX; ++k) {
        for (size_t r = 0; r < NY; ++r) {
            for (size_t c = 0; c < NX; ++c) {
                Ob((k * NY) + r, c) = CA(r, c);
            }
        }
        if (k + 1 < NX) {
            CA = CA * A;
        }
    }
    return Ob;
}

namespace detail {

/// Gramian Lyapunov solve, refused unless A is Hurwitz (continuous) / Schur (discrete).
/// For an unstable A the Lyapunov equation can still have a unique solution, but it
/// is not the (divergent) Gramian integral/sum.
template<size_t NX, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>>
gramian_lyapunov(const Matrix<NX, NX, T>& A, const Matrix<NX, NX, T>& Q, bool discrete) {
    const auto eig = mat::compute_eigenvalues(A);
    if (!eig.converged) {
        return damp::nullopt;
    }
    for (size_t i = 0; i < NX; ++i) {
        const bool stable = discrete ? (damp::abs(eig.values[i]) < T{1}) : (eig.values[i].real() < T{0});
        if (!stable) {
            return damp::nullopt;
        }
    }
    return discrete ? dlyap(A, Q) : lyap(A, Q);
}

/// Orthonormal basis of the controllable subspace of (A, B) and its dimension.
template<size_t NX, typename T>
struct ControllableSubspace {
    Matrix<NX, NX, T> basis{};         ///< leading `dim` columns: an orthonormal basis of the subspace
    size_t            dim{0};          ///< its dimension (NX when (A, B) is controllable)
    bool              converged{true}; ///< every SVD converged
};

/**
 * @brief Controllable subspace by the block-Arnoldi staircase (the ctrbf construction).
 *
 * Stair 0 is range(B). Stair k adds the directions of A·(stair k−1) that are new: the
 * block is projected against the basis found so far (twice, so it is orthogonal to working
 * precision) and its left singular vectors above the tolerance are appended. A is scaled by
 * 1/‖A‖∞ and the basis columns have unit norm, so the test σ > @p tol is unit-free: a
 * direction counts when it is reached with a coupling of at least @p tol relative to the
 * plant's fastest rate; stair 0 is tested against tol·σₘₐₓ(B).
 *
 * The rank of the Krylov matrix [B, AB, …, Aⁿ⁻¹B] answers the same question in exact
 * arithmetic, but its columns align with the dominant mode and the matrix loses rank
 * numerically from n ≈ 10, or sooner when the modes spread over decades. Each stair here is
 * tested on its own, so the verdict does not degrade with n.
 *
 * @see Van Dooren, "The generalized eigenstructure problem in linear system theory,"
 *      IEEE TAC 26(1), 1981
 * @note Compare with MATLAB®'s ctrbf.
 */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr ControllableSubspace<NX, T>
controllable_subspace(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B, T tol) {
    constexpr size_t            K = (NX < NU) ? NX : NU;
    ControllableSubspace<NX, T> out{};
    const T                     a_norm = mat::infinity_norm(A);
    const Matrix<NX, NX, T>     As = (a_norm > T{0}) ? Matrix<NX, NX, T>(A * (T{1} / a_norm)) : A;
    Matrix<NX, NU, T>           W = B;
    for (size_t stair = 0; stair < NX && out.dim < NX; ++stair) {
        for (int pass = 0; pass < 2; ++pass) { // classical Gram–Schmidt, twice
            for (size_t j = 0; j < out.dim; ++j) {
                for (size_t c = 0; c < NU; ++c) {
                    T dot{0};
                    for (size_t i = 0; i < NX; ++i) {
                        dot += out.basis(i, j) * W(i, c);
                    }
                    for (size_t i = 0; i < NX; ++i) {
                        W(i, c) -= dot * out.basis(i, j);
                    }
                }
            }
        }
        const auto s = mat::svd(W);
        if (!s.converged) {
            out.converged = false;
            return out;
        }
        const T floor = tol * ((stair == 0) ? s.singular_values[0] : T{1});
        size_t  r = 0;
        for (size_t k = 0; k < K; ++k) {
            if (s.singular_values[k] > floor) {
                ++r;
            }
        }
        if (r == 0) {
            break;
        }
        r = damp::min(r, NX - out.dim);
        for (size_t k = 0; k < r; ++k) {
            for (size_t i = 0; i < NX; ++i) {
                out.basis(i, out.dim + k) = s.singular_U(i, k);
            }
        }
        Matrix<NX, NU, T> next{};
        for (size_t k = 0; k < r; ++k) {
            for (size_t i = 0; i < NX; ++i) {
                T acc{0};
                for (size_t m = 0; m < NX; ++m) {
                    acc += As(i, m) * out.basis(m, out.dim + k);
                }
                next(i, k) = acc;
            }
        }
        W = next;
        out.dim += r;
    }
    return out;
}

} // namespace detail

/**
 * @brief Continuous/discrete controllability Gramian W_c.
 *
 * Solves the Lyapunov equation whose solution measures how strongly each state
 * direction is excited by the input:
 * - continuous: @f$ A W_c + W_c A^\top + B B^\top = 0 @f$
 * - discrete:   @f$ A W_c A^\top - W_c + B B^\top = 0 @f$
 *
 * @f$ W_c \succ 0 @f$ iff @f$ (A,B) @f$ is controllable. The Gramian is finite only
 * for a stable @f$ A @f$ (Hurwitz / Schur); an unstable @f$ A @f$ returns
 * damp::nullopt even when the Lyapunov equation itself is solvable.
 *
 * @note Compare with MATLAB®'s @c gram(sys,'c').
 * @see observability_gramian(), damp::lyap(), damp::dlyap().
 *
 * @param A        State matrix (NX × NX).
 * @param B        Input matrix (NX × NU).
 * @param discrete false → continuous-time, true → discrete-time.
 * @return W_c (NX × NX), or damp::nullopt if A is not stable or the solve is singular.
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>> controllability_gramian(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    bool                     discrete = false
) {
    return detail::gramian_lyapunov(A, Matrix<NX, NX, T>{B * B.transpose()}, discrete);
}

/**
 * @brief Continuous/discrete observability Gramian W_o.
 *
 * Solves the dual Lyapunov equation measuring how strongly each state direction
 * shows up in the output:
 * - continuous: @f$ A^\top W_o + W_o A + C^\top C = 0 @f$
 * - discrete:   @f$ A^\top W_o A - W_o + C^\top C = 0 @f$
 *
 * @f$ W_o \succ 0 @f$ iff @f$ (A,C) @f$ is observable. The Gramian is finite only
 * for a stable @f$ A @f$; an unstable @f$ A @f$ returns damp::nullopt.
 *
 * @note Compare with MATLAB®'s @c gram(sys,'o').
 * @see controllability_gramian(), damp::lyap(), damp::dlyap().
 *
 * @param A        State matrix (NX × NX).
 * @param C        Output matrix (NY × NX).
 * @param discrete false → continuous-time, true → discrete-time.
 * @return W_o (NX × NX), or damp::nullopt if A is not stable or the solve is singular.
 */
template<size_t NX, size_t NY, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>> observability_gramian(
    const Matrix<NX, NX, T>& A,
    const Matrix<NY, NX, T>& C,
    bool                     discrete = false
) {
    return detail::gramian_lyapunov(A.transpose(), Matrix<NX, NX, T>{C.transpose() * C}, discrete);
}

/**
 * @brief Numerical rank of a matrix (forwards to mat::rank)
 *
 * @param M   Input matrix
 * @param tol Relative tolerance: singular values ≤ tol·σ_max count as zero
 *            (default: mat::rank's eps-based default)
 * @return size_t rank of the matrix
 */
template<size_t R, size_t C, typename T>
[[nodiscard]] constexpr size_t rank(const Matrix<R, C, T>& M) {
    return mat::rank(M);
}

/// @overload explicit relative tolerance
template<size_t R, size_t C, typename T>
[[nodiscard]] constexpr size_t rank(const Matrix<R, C, T>& M, T tol) {
    return mat::rank(M, tol);
}

/**
 * @brief Check if a system is controllable
 *
 * The controllable subspace (detail::controllable_subspace, the staircase construction)
 * has dimension NX. The tolerance NX·64·eps is relative to the plant's fastest rate and
 * to B, so the answer is invariant to the signal and time units.
 *
 * @param A State matrix
 * @param B Input matrix
 * @return true if (A, B) is controllable
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr bool is_controllable(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B) {
    return detail::controllable_subspace(A, B, static_cast<T>(NX) * T{64} * eps<T>()).dim == NX;
}

/**
 * @brief Check if a system is observable
 *
 * Dual of is_controllable(): the controllable subspace of (Aᵀ, Cᵀ) has dimension NX.
 *
 * @param A State matrix
 * @param C Output matrix
 * @return true if (A, C) is observable
 */
template<size_t NX, size_t NY, typename T = double>
[[nodiscard]] constexpr bool is_observable(const Matrix<NX, NX, T>& A, const Matrix<NY, NX, T>& C) {
    return detail::controllable_subspace(Matrix<NX, NX, T>(A.transpose()), Matrix<NX, NY, T>(C.transpose()), static_cast<T>(NX) * T{64} * eps<T>()).dim == NX;
}

/**
 * @brief Check if a discrete-time system matrix A is stable
 *
 * A discrete system is stable if all eigenvalues have magnitude less than 1 (inside unit circle).
 *
 * @tparam N   Number of states
 * @tparam T   Numeric type (default: double)
 * @param A    State matrix to check
 *
 * @return true if all eigenvalues satisfy |λ| < 1, false otherwise
 */
template<size_t N, typename T = double>
[[nodiscard]] constexpr bool is_stable_discrete(const Matrix<N, N, T>& A) {
    auto eigen = mat::compute_eigenvalues(A);
    if (!eigen.converged) {
        return false;
    }

    for (size_t i = 0; i < N; ++i) {
        T magnitude = damp::abs(eigen.values[i]);
        if (magnitude >= T{1}) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Check closed-loop stability for discrete system with state feedback
 *
 * Checks stability of the closed-loop system A_cl = A - B*K with feedback u = -K*x.
 *
 * @tparam NX  Number of states
 * @tparam NU  Number of inputs
 * @tparam T   Numeric type (default: double)
 * @param A    State matrix
 * @param B    Control input matrix
 * @param K    State feedback gain matrix (u = -K*x)
 *
 * @return true if closed-loop system is stable, false otherwise
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr bool is_closed_loop_stable_discrete(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NU, NX, T>& K
) {
    Matrix<NX, NX, T> A_cl = A - B * K;
    return is_stable_discrete(A_cl);
}

/**
 * @brief Compute stability margin for continuous system
 *
 * Returns the distance to the stability boundary (imaginary axis).
 * Computed as the negative of the most positive real eigenvalue part.
 * Positive values indicate stability; larger values indicate more stability margin.
 *
 * @tparam N   Number of states
 * @tparam T   Numeric type (default: double)
 * @param A    State matrix
 *
 * @return Stability margin (positive = stable, negative = unstable)
 */
template<size_t N, typename T = double>
[[nodiscard]] constexpr T stability_margin_continuous(const Matrix<N, N, T>& A) {
    auto eigen = mat::compute_eigenvalues(A);
    if (!eigen.converged) {
        return T{-1}; // Negative margin signals "not provably stable"
    }

    T max_real = eigen.values[0].real();
    for (size_t i = 1; i < N; ++i) {
        if (eigen.values[i].real() > max_real) {
            max_real = eigen.values[i].real();
        }
    }
    return -max_real; // Positive means stable, larger is more stable
}

/**
 * @brief Compute stability margin for discrete system
 *
 * Returns the distance to the stability boundary (unit circle).
 * Computed as 1 - (maximum magnitude eigenvalue).
 * Positive values indicate stability; larger values indicate more stability margin.
 *
 * @tparam N   Number of states
 * @tparam T   Numeric type (default: double)
 * @param A    State matrix
 *
 * @return Stability margin (positive = stable, negative = unstable)
 */
template<size_t N, typename T = double>
[[nodiscard]] constexpr T stability_margin_discrete(const Matrix<N, N, T>& A) {
    auto eigen = mat::compute_eigenvalues(A);
    if (!eigen.converged) {
        return T{-1}; // Return unstable indicator
    }

    T max_mag = T{0};
    for (size_t i = 0; i < N; ++i) {
        const T magnitude = eigen.values[i].abs();
        if (magnitude > max_mag) {
            max_mag = magnitude;
        }
    }
    return T{1} - max_mag; // Positive means stable, larger is more stable
}

/**
 * @brief Compute closed-loop poles (eigenvalues) with state feedback
 *
 * Computes the eigenvalues of the closed-loop state matrix A_cl = A - B*K.
 * These poles determine the closed-loop system dynamics.
 *
 * If the eigen solver does not converge, returns a zero vector rather than
 * uninitialized/garbage pole locations. Callers that set design @c success from
 * the Riccati/gain path alone (e.g. LQR) keep @c success true when K is valid;
 * a zeroed pole vector signals that pole reporting failed independently of K.
 *
 * @tparam NX  Number of states
 * @tparam NU  Number of inputs
 * @tparam T   Numeric type (default: double)
 * @param A    State matrix
 * @param B    Control input matrix
 * @param K    State feedback gain matrix (u = -K*x)
 *
 * @return Vector of closed-loop pole locations (eigenvalues as complex numbers);
 *         all zeros if eigenvalue computation does not converge
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr ColVec<NX, damp::complex<T>> closed_loop_poles(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NU, NX, T>& K
) {
    Matrix<NX, NX, T> A_cl = A - B * K;
    auto              eigen = mat::compute_eigenvalues(A_cl);
    if (!eigen.converged) {
        return ColVec<NX, damp::complex<T>>{};
    }
    return eigen.values;
}

} // namespace stability
} // namespace damp
