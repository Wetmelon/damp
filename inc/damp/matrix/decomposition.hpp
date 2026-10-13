// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file decomposition.hpp
 * @brief Matrix factorizations: Cholesky, LDLᵀ, LU (partial pivoting), and QR
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

/**
 * @brief True if A equals its transpose (real) or conjugate transpose (complex)
 *
 * Each pair may differ by N·64·eps times its own scale
 * max(|aᵢⱼ|, |aⱼᵢ|, √|aᵢᵢ aⱼⱼ|), so the test is invariant to scaling A and to
 * graded (diagonally scaled) matrices.
 */
template<size_t N, typename T>
constexpr bool is_symmetric_or_hermitian(const Matrix<N, N, T>& A) {
    using real_t = scalar_type_t<T>;
    constexpr real_t rel = static_cast<real_t>(N) * real_t{64} * eps<T>();
    for (size_t i = 0; i < N; ++i) {
        if constexpr (!std::is_floating_point_v<T>) {
            // Hermitian requires a real diagonal
            if (damp::abs(damp::imag(A(i, i))) > rel * damp::abs(A(i, i))) {
                return false;
            }
        }
        for (size_t j = i + 1; j < N; ++j) {
            const real_t scale = damp::max({damp::abs(A(i, j)), damp::abs(A(j, i)), damp::sqrt(damp::abs(A(i, i)) * damp::abs(A(j, j)))});
            if (damp::abs(A(i, j) - damp::conj(A(j, i))) > rel * scale) {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief Cholesky decomposition for positive-definite matrices
 *
 * For real matrices, computes L such that A = LLᵀ (A must be symmetric).
 * For complex matrices, computes L such that A = LLᴴ (A must be Hermitian).
 *
 * Reads only the lower triangle and assumes the upper mirrors it; symmetry is
 * not checked. Callers that may pass a non-symmetric matrix should gate on
 * is_symmetric_or_hermitian() first, as solve() does.
 *
 * Returns damp::nullopt when a pivot aᵢᵢ − Σₖ|lᵢₖ|² does not exceed 16·(i+1)·eps times
 * |aᵢᵢ| + Σₖ|lᵢₖ|², the magnitude of the i+1 terms that formed it. A singular
 * positive-semidefinite matrix leaves a roundoff-sized last pivot of either sign, and a
 * positive one would factor into entries of 1/√eps. The floor grows with each pivot's own
 * term count, not with N: the early pivots of a large graded matrix (an MPC Hessian, whose
 * first-move diagonal is 10⁴ times the last) carry only a few terms' roundoff, and an
 * N-wide floor refused them in float. No floor separates singular from positive definite
 * exactly, since the last pivot of a singular matrix also carries the forward error of the
 * leading block: about 3 % of random singular Gram matrices clear this one, against half
 * for a positivity test. Decide rank with mat::rank. The test is invariant to scaling A.
 * LAPACK's xPOTRF tests only positivity.
 *
 * @note Compare with MATLAB®'s chol(A, 'lower').
 * @see Golub & Van Loan, "Matrix Computations" (4th ed., 2013), §4.2
 *
 * @param A Symmetric positive-definite (real) or Hermitian positive-definite (complex) matrix
 * @return Lower-triangular L, or damp::nullopt if A is not positive definite
 */
template<size_t N, typename T>
constexpr damp::optional<Matrix<N, N, T>> cholesky(const Matrix<N, N, T>& A) {
    using real_t = scalar_type_t<T>;
    Matrix<N, N, T> L = Matrix<N, N, T>::zeros();

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j <= i; ++j) {

            T sum = A(i, j);
            for (size_t k = 0; k < j; ++k) {
                sum -= L(i, k) * damp::conj(L(j, k));
            }

            if (i == j) {
                real_t diag_val{};
                real_t mag = damp::abs(A(i, i)); // size of the terms that formed the pivot
                for (size_t k = 0; k < j; ++k) {
                    const real_t l = damp::abs(L(i, k));
                    mag += l * l;
                }
                if constexpr (std::is_floating_point_v<T>) {
                    diag_val = sum;
                } else {
                    diag_val = damp::real(sum);
                }
                if (!(diag_val > real_t{16} * static_cast<real_t>(i + 1) * eps<T>() * mag)) {
                    return damp::nullopt;
                }
                L(i, j) = damp::sqrt(diag_val);
            } else {
                L(i, j) = sum / L(j, j);
            }
        }
    }

    return L;
}

/// LDLᴴ factors: A = L D Lᴴ with L unit lower-triangular and D real diagonal.
template<size_t N, typename T>
struct LDLT {
    Matrix<N, N, T>                  L{}; ///< Unit lower-triangular factor
    damp::array<scalar_type_t<T>, N> d{}; ///< Diagonal of D (real; may carry negative entries)
};

/**
 * @brief LDLᴴ decomposition of a symmetric (Hermitian) matrix, no pivoting
 *
 * A = L D Lᴴ with L unit lower-triangular and D real diagonal. No square roots, and
 * D may carry negative entries, so this covers the symmetric indefinite matrices whose
 * leading principal minors are all nonzero: in particular the quasi-definite KKT
 * matrices @f$ \begin{bmatrix} H & A^\top \\ A & -E \end{bmatrix} @f$ with
 * @f$ H \succ 0,\; E \succeq 0 @f$, which factor in any symmetric ordering
 * (Vanderbei). A general indefinite matrix needs Bunch–Kaufman pivoting, which this does
 * not do: a pivot |dⱼ| at or below 16·(j+1)·eps of the terms that formed it (the floor
 * cholesky() uses) returns damp::nullopt rather than dividing by roundoff.
 *
 * Reads only the lower triangle and assumes the upper mirrors it.
 *
 * @note Compare with MATLAB®'s ldl(A) (which pivots).
 * @see Golub & Van Loan, "Matrix Computations" (4th ed., 2013), §4.1.2
 * @see Vanderbei, "Symmetric Quasidefinite Matrices," SIAM J. Optim. 5(1), 1995
 * @see ldlt_solve(), cholesky()
 *
 * @param A Symmetric (real) or Hermitian (complex) matrix
 * @return L and d, or damp::nullopt on a vanishing pivot
 */
template<size_t N, typename T>
constexpr damp::optional<LDLT<N, T>> ldlt(const Matrix<N, N, T>& A) {
    using real_t = scalar_type_t<T>;
    const auto real_part = [](const T& x) -> real_t {
        if constexpr (std::is_floating_point_v<T>) {
            return x;
        } else {
            return damp::real(x);
        }
    };

    LDLT<N, T> out{};
    out.L = Matrix<N, N, T>::identity();
    for (size_t j = 0; j < N; ++j) {
        // dⱼ = aⱼⱼ − Σₖ |lⱼₖ|² dₖ, tested against the size of those terms
        real_t dj = real_part(A(j, j));
        real_t mag = damp::abs(dj);
        for (size_t k = 0; k < j; ++k) {
            const real_t l = damp::abs(out.L(j, k));
            dj -= l * l * out.d[k];
            mag += l * l * damp::abs(out.d[k]);
        }
        if (!(damp::abs(dj) > real_t{16} * static_cast<real_t>(j + 1) * eps<T>() * mag)) {
            return damp::nullopt;
        }
        out.d[j] = dj;
        for (size_t i = j + 1; i < N; ++i) {
            T sum = A(i, j);
            for (size_t k = 0; k < j; ++k) {
                sum -= out.L(i, k) * damp::conj(out.L(j, k)) * out.d[k];
            }
            out.L(i, j) = sum / dj;
        }
    }
    return out;
}

/**
 * @brief LU decomposition with partial pivoting
 *
 * Factors @f$ P A = L U @f$ where @f$ L @f$ is unit lower-triangular,
 * @f$ U @f$ is upper-triangular, and @f$ P @f$ is encoded by the pivot vector
 * (@c piv[i] is the original row moved into position @c i).
 *
 * Returns damp::nullopt when the chosen pivot is zero or lost to cancellation:
 * below N·64·eps times the magnitude |aᵢⱼ| + Σₖ|lᵢₖ||uₖⱼ| of the terms that
 * formed it. The test is invariant to row and column scaling, so graded but
 * well-posed systems such as [[1, 10¹⁶], [0, 1]] still factor.
 *
 * @note Compare with MATLAB®'s [L,U,P] = lu(A) (row-permutation form).
 * @see Golub & Van Loan, "Matrix Computations" (4th ed., 2013), §3.2, §3.3.1
 */
template<size_t N, typename T>
constexpr damp::optional<damp::tuple<Matrix<N, N, T>, Matrix<N, N, T>, damp::array<size_t, N>>>
lu_decomposition(const Matrix<N, N, T>& A) {
    using real_t = scalar_type_t<T>;
    constexpr real_t rel = static_cast<real_t>(N) * real_t{64} * eps<T>();

    Matrix<N, N, T>        L = Matrix<N, N, T>::identity();
    Matrix<N, N, T>        U = A;
    Matrix<N, N, real_t>   mag; // |A| + Σ|l||u|: the size of the terms behind each U entry
    damp::array<size_t, N> piv;
    for (size_t i = 0; i < N; ++i) {
        piv[i] = i;
        for (size_t j = 0; j < N; ++j) {
            mag(i, j) = damp::abs(A(i, j));
        }
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

        if (!(max_val > rel * mag(max_row, i))) {
            return damp::nullopt;
        }

        if (max_row != i) {
            // Swap rows in U and its magnitude bound
            for (size_t col = 0; col < N; ++col) {
                damp::swap(U(i, col), U(max_row, col));
                damp::swap(mag(i, col), mag(max_row, col));
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
                mag(j, k) += damp::abs(factor) * damp::abs(U(i, k));
            }
        }
    }

    return damp::make_tuple(L, U, piv);
}

namespace detail {

/**
 * @brief Householder vector that zeroes X(k+1…N−1, k) when applied as H = I − β·v·vᴴ.
 *
 * Reflects X(k…, k) onto α·e_k with |α| = ‖X(k…, k)‖ and the phase of α opposite to that
 * of X(k, k), so v[k] = X(k, k) − α never cancels; for real T that is the usual sign choice.
 * Entries of @p v below k are left as they were (the reflector acts on rows k…N−1 only).
 *
 * @see "Matrix Computations" (Golub & Van Loan, 2013), §5.1.2–5.1.3 and §5.1.13 (complex)
 *
 * @return β = 2/‖v‖², or 0 when the column is already zero from row k down (no reflection).
 */
template<typename T, size_t N, size_t M>
[[nodiscard]] constexpr scalar_type_t<T> householder_vector(const Matrix<N, M, T>& X, size_t k, damp::array<T, N>& v) {
    using real_t = scalar_type_t<T>;
    real_t norm_sq = real_t{0};
    for (size_t i = k; i < N; ++i) {
        const real_t a = damp::abs(X(i, k));
        norm_sq += a * a;
    }
    const real_t nrm = damp::sqrt(norm_sq);
    if (nrm == real_t{0}) {
        return real_t{0};
    }
    const real_t r0 = damp::abs(X(k, k));
    const T      phase = (r0 > real_t{0}) ? (X(k, k) / r0) : T{1};
    const T      alpha = -(phase * nrm);
    v[k] = X(k, k) - alpha;
    for (size_t i = k + 1; i < N; ++i) {
        v[i] = X(i, k);
    }
    real_t vhv = real_t{0};
    for (size_t i = k; i < N; ++i) {
        const real_t a = damp::abs(v[i]);
        vhv += a * a;
    }
    return (vhv == real_t{0}) ? real_t{0} : (real_t{2} / vhv);
}

/// X ← (I − β·v·vᴴ)·X on rows k…N−1, columns c0…C−1.
template<typename T, size_t N, size_t C>
constexpr void apply_householder_left(Matrix<N, C, T>& X, const damp::array<T, N>& v, scalar_type_t<T> beta, size_t k, size_t c0) {
    for (size_t j = c0; j < C; ++j) {
        T dot = T{0};
        for (size_t i = k; i < N; ++i) {
            dot += damp::conj(v[i]) * X(i, j);
        }
        const T s = beta * dot;
        for (size_t i = k; i < N; ++i) {
            X(i, j) -= s * v[i];
        }
    }
}

} // namespace detail

/**
 * @brief Thin QR factorization A = QR (Householder)
 *
 * @f$ Q @f$ is @f$ N \times M @f$ with orthonormal columns; @f$ R @f$ is
 * @f$ M \times M @f$ upper-triangular with a real, non-negative diagonal.
 *
 * @see Golub & Van Loan, "Matrix Computations" (4th ed., 2013), §5.2.1–5.2.2
 * @see full_qr() for the full unitary factor
 */
template<typename T, size_t N, size_t M>
struct QRDecomposition {
    Matrix<N, M, T> Q{}; ///< Orthonormal columns (thin Q)
    Matrix<M, M, T> R{}; ///< Upper-triangular factor

    /**
     * @brief True if A had full column rank (no diagonal of R was zeroed as rank loss)
     */
    [[nodiscard]] constexpr bool is_valid() const {
        for (size_t i = 0; i < M; ++i) {
            if (R(i, i) == T{0}) {
                return false;
            }
        }
        return true;
    }
};

/**
 * @brief Thin QR via Householder reflections
 *
 * Reduces A to upper-triangular form by reflections H₀…H_{s−1}, s = min(M, N−1), and forms
 * the thin Q = H₀⋯H_{s−1}·[I; 0] by backward accumulation. Q is orthonormal to working
 * precision whatever the conditioning of A: ‖QᴴQ − I‖ is a modest multiple of ε. Modified
 * Gram–Schmidt, used here before, loses orthogonality in proportion to κ(A)·ε, so a
 * conditioning of 10⁸ left columns of Q 10⁻⁸ away from orthogonal.
 *
 * Each column of Q and row of R is then rotated by the phase of its diagonal entry, so R has
 * a real, non-negative diagonal (the convention Gram–Schmidt produces) and QR is unchanged.
 *
 * Rank loss is scale-invariant per column: when |R(j, j)|, the part of column j orthogonal to
 * the columns before it, is at most rel_tol times ‖A(:, j)‖, it is roundoff and R(j, j) is
 * stored as zero, which is what is_valid() reports. Q keeps an orthonormal column there, so
 * QR reconstructs A to within that dropped roundoff. Columns past N (M > N) have R(j, j) = 0
 * and a zero Q column.
 *
 * @note Compare with MATLAB®'s [Q,R] = qr(A, 0) (economy size).
 * @see Golub & Van Loan, "Matrix Computations" (4th ed., 2013), §5.1.6 (backward accumulation) and §5.2.2
 * @see "Accuracy and Stability of Numerical Algorithms" (Higham, 2002), §19.3 (Householder QR) and §19.8 (MGS orthogonality loss)
 *
 * @tparam T   Element type
 * @tparam N   Rows
 * @tparam M   Columns
 * @param A        Input matrix
 * @param rel_tol  Relative column-norm floor treated as rank loss (default N·64·eps)
 * @return Thin QR factors
 */
template<typename T, size_t N, size_t M>
[[nodiscard]] constexpr QRDecomposition<T, N, M> qr_decompose(
    const Matrix<N, M, T>& A,
    scalar_type_t<T>       rel_tol = static_cast<scalar_type_t<T>>(N) * scalar_type_t<T>{64} * eps<T>()
) {
    using real_t = scalar_type_t<T>;
    constexpr size_t steps = (M < N) ? M : (N - 1); // reflections; the last row needs none
    constexpr size_t K = (M < N) ? M : N;           // diagonal entries R can carry

    Matrix<N, M, T>                                         W = A;
    damp::array<damp::array<T, N>, (steps > 0 ? steps : 1)> vs{};
    damp::array<real_t, (steps > 0 ? steps : 1)>            betas{};
    for (size_t k = 0; k < steps; ++k) {
        betas[k] = detail::householder_vector(W, k, vs[k]);
        if (betas[k] != real_t{0}) {
            detail::apply_householder_left(W, vs[k], betas[k], k, k);
        }
    }

    // Thin Q = H₀⋯H_{s−1}·[I; 0], applying the reflections from the last to the first.
    QRDecomposition<T, N, M> result;
    for (size_t j = 0; j < K; ++j) {
        result.Q(j, j) = T{1};
    }
    for (size_t kk = steps; kk-- > 0;) {
        if (betas[kk] != real_t{0}) {
            detail::apply_householder_left(result.Q, vs[kk], betas[kk], kk, kk);
        }
    }

    for (size_t i = 0; i < K; ++i) {
        for (size_t j = i; j < M; ++j) {
            result.R(i, j) = W(i, j);
        }
    }

    // Real, non-negative diagonal: Q(:, i) ·= φ, R(i, :) ·= φ̄ with φ the phase of R(i, i).
    for (size_t i = 0; i < K; ++i) {
        const real_t r = damp::abs(result.R(i, i));
        if (r == real_t{0}) {
            continue;
        }
        const T phase = result.R(i, i) / r;
        for (size_t j = i; j < M; ++j) {
            result.R(i, j) = damp::conj(phase) * result.R(i, j);
        }
        result.R(i, i) = T{r};
        for (size_t row = 0; row < N; ++row) {
            result.Q(row, i) = result.Q(row, i) * phase;
        }
    }

    // Rank loss, relative to each column of A.
    for (size_t j = 0; j < K; ++j) {
        real_t sq = real_t{0};
        for (size_t i = 0; i < N; ++i) {
            const real_t a = damp::abs(A(i, j));
            sq += a * a;
        }
        if (!(damp::abs(result.R(j, j)) > rel_tol * damp::sqrt(sq))) {
            result.R(j, j) = T{0};
        }
    }
    return result;
}

/**
 * @brief Result of a full (complete) QR factorization.
 *
 * A = Q·R with Q an N×N unitary (orthogonal, for real T) matrix and R N×M
 * upper-triangular. Unlike qr_decompose (thin, N×M Q), this retains the *full*
 * unitary factor: its leading rank columns span the
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
        // Householder vector v that zeroes R(k+1.., k); H = I − β·v·vᴴ is Hermitian unitary.
        damp::array<T, N> v{};
        const real_t      beta = detail::householder_vector(out.R, k, v);
        if (beta == real_t{0}) {
            continue; // column already zero below the diagonal
        }
        // Apply H from the left: R −= β·v·(vᴴR).
        detail::apply_householder_left(out.R, v, beta, k, 0);
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
