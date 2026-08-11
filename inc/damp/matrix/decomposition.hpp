// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file decomposition.hpp
 * @brief Matrix factorizations: Cholesky, LU (partial pivoting), and QR
 *        (modified Gram-Schmidt), plus the symmetry/Hermitian predicate used to
 *        choose between them. Linear-system solvers built on these live in
 *        solve.hpp; the eigensolver in eigen.hpp.
 */

#include <cstddef>
#include <type_traits>

#include "core.hpp"
#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "matrix_traits.hpp"

namespace damp {

namespace mat {

template<size_t N, typename T>
constexpr bool is_symmetric_or_hermitian(const Matrix<N, N, T>& A) {
    if constexpr (std::is_floating_point_v<T>) {
        constexpr auto tol = default_tol<T>();
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = i + 1; j < N; ++j) {
                if (damp::abs(A(i, j) - A(j, i)) > tol) {
                    return false;
                }
            }
        }
        return true;
    } else {
        constexpr auto tol = default_tol<T>();
        constexpr auto tol_sq = tol * tol;
        // Hermitian requires real diagonal elements
        for (size_t i = 0; i < N; ++i) {
            auto imag_diag = damp::imag(A(i, i));
            if (imag_diag * imag_diag > tol_sq) {
                return false;
            }
        }
        // Off-diagonal: A(i,j) == conj(A(j,i))
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = i + 1; j < N; ++j) {
                auto diff = A(i, j) - damp::conj(A(j, i));
                auto mag_sq = damp::real(diff * damp::conj(diff));
                if (mag_sq > tol_sq) {
                    return false;
                }
            }
        }
        return true;
    }
}

/**
 * @brief Cholesky decomposition for positive-definite matrices
 *
 * For real matrices, computes L such that A = LLᵀ (A must be symmetric).
 * For complex matrices, computes L such that A = LLᴴ (A must be Hermitian).
 *
 * Reads only the lower triangle and assumes the upper mirrors it; symmetry is
 * not checked. Callers that may pass a non-symmetric matrix should gate on
 * is_symmetric_or_hermitian() first, as solve() does. The pivot-positivity test
 * is the PD guard, so non-PD input is still rejected via nullopt.
 *
 * @note Compare with MATLAB®'s chol(A, 'lower').
 * @see Golub & Van Loan, "Matrix Computations" (4th ed., 2013), §4.2
 *
 * @param A Symmetric positive-definite (real) or Hermitian positive-definite (complex) matrix
 * @return Lower-triangular L, or damp::nullopt if A is not positive definite
 */
template<size_t N, typename T>
constexpr damp::optional<Matrix<N, N, T>> cholesky(const Matrix<N, N, T>& A) {
    Matrix<N, N, T> L = Matrix<N, N, T>::zeros();

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j <= i; ++j) {

            T sum = A(i, j);
            for (size_t k = 0; k < j; ++k) {
                sum -= L(i, k) * damp::conj(L(j, k));
            }

            if (i == j) {
                if constexpr (std::is_floating_point_v<T>) {
                    auto diag_val = sum;
                    if (diag_val <= 0) {
                        return damp::nullopt;
                    }
                    L(i, j) = damp::sqrt(diag_val);
                } else {
                    auto diag_val = damp::real(sum);
                    if (diag_val <= 0) {
                        return damp::nullopt;
                    }
                    L(i, j) = damp::sqrt(diag_val);
                }
            } else {
                L(i, j) = sum / L(j, j);
            }
        }
    }

    return L;
}

/**
 * @brief LU decomposition with partial pivoting
 *
 * Factors @f$ P A = L U @f$ where @f$ L @f$ is unit lower-triangular,
 * @f$ U @f$ is upper-triangular, and @f$ P @f$ is encoded by the pivot vector
 * (@c piv[i] is the original row moved into position @c i).
 *
 * Returns damp::nullopt if a pivot magnitude falls below default_tol
 * (numerically singular).
 *
 * @note Compare with MATLAB®'s [L,U,P] = lu(A) (row-permutation form).
 * @see Golub & Van Loan, "Matrix Computations" (4th ed., 2013), §3.2
 */
template<size_t N, typename T>
constexpr damp::optional<damp::tuple<Matrix<N, N, T>, Matrix<N, N, T>, damp::array<size_t, N>>>
lu_decomposition(const Matrix<N, N, T>& A) {
    Matrix<N, N, T>        L = Matrix<N, N, T>::identity();
    Matrix<N, N, T>        U = A;
    damp::array<size_t, N> piv;
    for (size_t i = 0; i < N; ++i) {
        piv[i] = i;
    }

    for (size_t i = 0; i < N; ++i) {
        // Partial pivot
        size_t max_row = i;
        auto   max_val = damp::abs(U(i, i));
        for (size_t r = i + 1; r < N; ++r) {
            auto val = damp::abs(U(r, i));
            if (val > max_val) {
                max_val = val;
                max_row = r;
            }
        }

        if (max_val < default_tol<T>()) {
            return damp::nullopt;
        }

        if (max_row != i) {
            // Swap rows in U
            for (size_t col = 0; col < N; ++col) {
                damp::swap(U(i, col), U(max_row, col));
            }
            // Swap previous columns in L
            for (size_t col = 0; col < i; ++col) {
                damp::swap(L(i, col), L(max_row, col));
            }
            // Track pivot
            damp::swap(piv[i], piv[max_row]);
        }

        // Elimination
        for (size_t j = i + 1; j < N; ++j) {
            auto factor = U(j, i) / U(i, i);
            L(j, i) = factor;
            for (size_t k = i; k < N; ++k) {
                U(j, k) -= factor * U(i, k);
            }
        }
    }

    return damp::make_tuple(L, U, piv);
}

/**
 * @brief Thin QR factorization A = QR (modified Gram–Schmidt)
 *
 * @f$ Q @f$ is @f$ N \times M @f$ with orthonormal columns (when full column
 * rank); @f$ R @f$ is @f$ M \times M @f$ upper-triangular.
 *
 * @see Golub & Van Loan, "Matrix Computations" (4th ed., 2013), §5.2.8
 * @see full_qr() for the Householder full unitary factor
 */
template<typename T, size_t N, size_t M>
struct QRDecomposition {
    Matrix<N, M, T> Q{}; ///< Orthonormal columns (thin Q)
    Matrix<M, M, T> R{}; ///< Upper-triangular factor

    /**
     * @brief True if every diagonal entry of R exceeds default_tol
     */
    [[nodiscard]] constexpr bool is_valid() const {
        constexpr auto tol = default_tol<T>();
        for (size_t i = 0; i < M; ++i) {
            if (damp::abs(R(i, i)) < tol) {
                return false;
            }
        }
        return true;
    }
};

/**
 * @brief Thin QR via modified Gram–Schmidt
 *
 * @note Compare with MATLAB®'s [Q,R] = qr(A, 0) (economy size).
 * @see Golub & Van Loan, "Matrix Computations" (4th ed., 2013), §5.2.8
 *
 * @tparam T   Element type
 * @tparam N   Rows
 * @tparam M   Columns
 * @param A    Input matrix
 * @param eps  Column-norm floor treated as rank loss (default_tol by default)
 * @return Thin QR factors
 */
template<typename T, size_t N, size_t M>
[[nodiscard]] constexpr QRDecomposition<T, N, M> qr_decompose(
    const Matrix<N, M, T>& A,
    scalar_type_t<T>       eps = default_tol<T>()
) {
    using real_t = scalar_type_t<T>;
    QRDecomposition<T, N, M> result;
    result.Q = A;

    for (size_t j = 0; j < M; ++j) {
        real_t col_norm_sq = real_t{0};
        for (size_t i = 0; i < N; ++i) {
            const real_t a = damp::abs(result.Q(i, j));
            col_norm_sq += a * a;
        }
        const real_t col_norm = damp::sqrt(col_norm_sq);
        result.R(j, j) = static_cast<T>(col_norm);

        if (col_norm > eps) {
            for (size_t i = 0; i < N; ++i) {
                result.Q(i, j) /= result.R(j, j);
            }

            for (size_t k = j + 1; k < M; ++k) {
                // Hermitian inner product for complex T
                T dot_prod = T{0};
                for (size_t i = 0; i < N; ++i) {
                    dot_prod += damp::conj(result.Q(i, j)) * result.Q(i, k);
                }
                result.R(j, k) = dot_prod;

                for (size_t i = 0; i < N; ++i) {
                    result.Q(i, k) -= result.R(j, k) * result.Q(i, j);
                }
            }
        }
    }

    return result;
}

/**
 * @brief Result of a full (complete) QR factorization.
 *
 * A = Q·R with Q an N×N unitary (orthogonal, for real T) matrix and R N×M
 * upper-triangular. Unlike qr_decompose (thin modified-Gram-Schmidt, N×M Q),
 * this retains the *full* unitary factor: its leading rank columns span the
 * range of A and its trailing columns span the orthogonal complement (the kernel
 * of Aᴴ). Those complement columns are what robust pole placement needs.
 */
template<typename T, size_t N, size_t M>
struct FullQR {
    Matrix<N, N, T> Q{}; ///< Unitary factor (full N×N).
    Matrix<N, M, T> R{}; ///< Upper-triangular factor (N×M).
};

/**
 * @brief Full QR factorization via Householder reflections (real or complex T).
 *
 * Numerically robust (unitary reflections, not Gram-Schmidt). For A with full
 * column rank M ≤ N, columns 0..M−1 of Q form an orthonormal basis of range(A)
 * and columns M..N−1 an orthonormal basis of its complement (Qᴴ_⊥·A = 0). The
 * Householder vector uses a complex phase so v[k] avoids cancellation; for real
 * T this reduces to the usual sign choice and the reflector to I − 2·v·vᵀ/(vᵀv).
 *
 * @see Golub & Van Loan, "Matrix Computations" (4th ed., 2013), §5.2 and §5.1.13
 *      (complex Householder)
 */
template<typename T, size_t N, size_t M>
[[nodiscard]] constexpr FullQR<T, N, M> full_qr(const Matrix<N, M, T>& A) {
    using real_t = scalar_type_t<T>;
    FullQR<T, N, M> out;
    out.Q = Matrix<N, N, T>::identity();
    out.R = A;

    constexpr size_t steps = (M < N) ? M : (N - 1); // columns to reflect
    for (size_t k = 0; k < steps; ++k) {
        // Householder vector v that zeroes R(k+1.., k). The reflector is the
        // Hermitian unitary H = I − β·v·vᴴ with real β = 2/‖v‖².
        real_t norm_sq = real_t{0};
        for (size_t i = k; i < N; ++i) {
            const real_t a = damp::abs(out.R(i, k));
            norm_sq += a * a;
        }
        const real_t nrm = damp::sqrt(norm_sq);
        if (nrm == real_t{0}) {
            continue; // column already zero below the diagonal
        }
        // Reflect R(k,k) onto α·e_k with |α| = ‖x‖; the phase −R(k,k)/|R(k,k)|
        // maximizes |v[k]| (avoids cancellation). For real T this is the sign.
        const real_t r0 = damp::abs(out.R(k, k));
        const T      phase = (r0 > real_t{0}) ? (out.R(k, k) / r0) : T{1};
        const T      alpha = -(phase * nrm);

        damp::array<T, N> v{};
        v[k] = out.R(k, k) - alpha;
        for (size_t i = k + 1; i < N; ++i) {
            v[i] = out.R(i, k);
        }
        real_t vhv = real_t{0};
        for (size_t i = k; i < N; ++i) {
            const real_t a = damp::abs(v[i]);
            vhv += a * a;
        }
        if (vhv == real_t{0}) {
            continue;
        }
        const real_t beta = real_t{2} / vhv;

        // Apply H from the left: R −= β·v·(vᴴR).
        for (size_t j = 0; j < M; ++j) {
            T dot = T{0};
            for (size_t i = k; i < N; ++i) {
                dot += damp::conj(v[i]) * out.R(i, j);
            }
            const T s = beta * dot;
            for (size_t i = k; i < N; ++i) {
                out.R(i, j) -= s * v[i];
            }
        }
        // Accumulate Q = Q·H from the right: Q −= β·(Q·v)·vᴴ.
        for (size_t i = 0; i < N; ++i) {
            T dot = T{0};
            for (size_t j = k; j < N; ++j) {
                dot += out.Q(i, j) * v[j];
            }
            const T s = beta * dot;
            for (size_t j = k; j < N; ++j) {
                out.Q(i, j) -= s * damp::conj(v[j]);
            }
        }
    }
    return out;
}

} // namespace mat
} // namespace damp
