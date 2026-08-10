// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file svd.hpp
 * @brief Singular value decomposition via the one-sided Jacobi method, plus the
 *        two routines it underpins: the Moore–Penrose pseudoinverse and an
 *        orthonormal null-space basis. Works for real and complex matrices of
 *        any shape (tall, wide, square, rank-deficient).
 *
 * The one-sided Jacobi SVD is chosen over the bidiagonal Golub–Kahan SVD because
 * it is simple, constexpr-friendly, and computes small singular values with high
 * relative accuracy — exactly the property that makes rank, null-space, and
 * pseudoinverse robust. For the matrix sizes typical of control design (n ≲ 20)
 * its O(n³) per sweep cost is negligible.
 *
 * @see Golub & Van Loan, "Matrix Computations" (4th ed., 2013), §8.6 (Jacobi SVD)
 * @see J. Demmel, K. Veselić, "Jacobi's method is more accurate than QR,"
 *      SIAM J. Matrix Anal. Appl. 13(4), 1992,
 *      https://doi.org/10.1137/0613074
 */

#include <cstddef>

#include "core.hpp"
#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "matrix_traits.hpp"

namespace damp {

namespace mat {

namespace detail {

/// Result of a one-sided Jacobi SVD of a tall matrix (P ≥ Q).
template<size_t P, size_t Q, typename T>
struct TallSVD {
    Matrix<P, P, T>                  U{}; ///< Left singular vectors (P×P unitary).
    damp::array<scalar_type_t<T>, Q> s{}; ///< Singular values, descending, ≥ 0.
    Matrix<Q, Q, T>                  V{}; ///< Right singular vectors (Q×Q unitary).
    bool                             converged = false;
};

/**
 * @brief One-sided Jacobi SVD of a tall/square matrix A (P×Q, P ≥ Q).
 *
 * Orthogonalizes the columns of A by a sequence of two-sided rotations; the
 * orthogonalized column norms are the singular values, the normalized columns
 * the (economy) left vectors, and the accumulated rotations the right vectors.
 * The left factor is then completed to a full P×P unitary by Gram–Schmidt so
 * the result is a genuine full SVD even when A is rank-deficient.
 */
template<size_t P, size_t Q, typename T>
[[nodiscard]] constexpr TallSVD<P, Q, T> jacobi_svd_tall(const Matrix<P, Q, T>& A) {
    static_assert(P >= Q, "jacobi_svd_tall requires P >= Q (transpose wide inputs)");
    using real_t = scalar_type_t<T>;
    constexpr real_t eps = default_tol<T>();

    Matrix<P, Q, T> W = A;
    Matrix<Q, Q, T> V = Matrix<Q, Q, T>::identity();

    // Sweep over column pairs, rotating each pair to mutual orthogonality.
    constexpr size_t max_sweeps = 60;
    bool             converged = false;
    for (size_t sweep = 0; sweep < max_sweeps; ++sweep) {
        real_t off = real_t{0};
        for (size_t p = 0; p < Q; ++p) {
            for (size_t q = p + 1; q < Q; ++q) {
                // 2×2 Hermitian Gram block of columns p, q.
                real_t alpha = real_t{0};
                real_t beta = real_t{0};
                T      gamma = T{0};
                for (size_t i = 0; i < P; ++i) {
                    const T      wip = W(i, p);
                    const T      wiq = W(i, q);
                    const real_t aip = damp::abs(wip);
                    const real_t aiq = damp::abs(wiq);
                    alpha += aip * aip;
                    beta += aiq * aiq;
                    gamma += damp::conj(wip) * wiq;
                }
                const real_t gmag = damp::abs(gamma);
                const real_t scale = damp::sqrt(alpha * beta);
                if (scale == real_t{0} || gmag <= eps * scale) {
                    continue; // already orthogonal
                }
                const real_t rel = gmag / scale;
                if (rel > off) {
                    off = rel;
                }

                // Real Jacobi angle that zeros the off-diagonal of the
                // phase-aligned block [[alpha, gmag], [gmag, beta]]. The smaller
                // root t = -sign(zeta)/(|zeta|+sqrt(zeta²+1)) keeps the rotation
                // near identity for stability.
                const real_t zeta = (beta - alpha) / (real_t{2} * gmag);
                const real_t zsign = (zeta >= real_t{0}) ? real_t{1} : real_t{-1};
                const real_t t = -zsign / (damp::abs(zeta) + damp::sqrt((zeta * zeta) + real_t{1}));
                const real_t c = real_t{1} / damp::sqrt(real_t{1} + (t * t));
                const real_t s = c * t;
                // ph = e^{-i·arg(gamma)} folds the complex phase into the
                // rotation; for real T it is simply ±1.
                const T ph = damp::conj(gamma) / gmag;
                const T phs = ph * s;
                const T phc = ph * c;

                for (size_t i = 0; i < P; ++i) {
                    const T wip = W(i, p);
                    const T wiq = W(i, q);
                    W(i, p) = (c * wip) + (phs * wiq);
                    W(i, q) = (-s * wip) + (phc * wiq);
                }
                for (size_t i = 0; i < Q; ++i) {
                    const T vip = V(i, p);
                    const T viq = V(i, q);
                    V(i, p) = (c * vip) + (phs * viq);
                    V(i, q) = (-s * vip) + (phc * viq);
                }
            }
        }
        if (off <= eps) {
            converged = true;
            break;
        }
    }

    // Singular values = orthogonalized column norms; economy U = W normalized.
    TallSVD<P, Q, T> out;
    out.V = V;
    out.converged = converged;
    real_t smax = real_t{0};
    for (size_t j = 0; j < Q; ++j) {
        real_t nrm_sq = real_t{0};
        for (size_t i = 0; i < P; ++i) {
            const real_t a = damp::abs(W(i, j));
            nrm_sq += a * a;
        }
        out.s[j] = damp::sqrt(nrm_sq);
        if (out.s[j] > smax) {
            smax = out.s[j];
        }
    }

    // Sort singular values descending, permuting Ue (= W) and V columns with them.
    for (size_t i = 0; i < Q; ++i) {
        size_t best = i;
        for (size_t j = i + 1; j < Q; ++j) {
            if (out.s[j] > out.s[best]) {
                best = j;
            }
        }
        if (best != i) {
            damp::swap(out.s[i], out.s[best]);
            for (size_t r = 0; r < P; ++r) {
                damp::swap(W(r, i), W(r, best));
            }
            for (size_t r = 0; r < Q; ++r) {
                damp::swap(out.V(r, i), out.V(r, best));
            }
        }
    }

    // Normalize the reliable (nonzero) singular columns into U.
    const real_t sv_tol = smax * eps;
    size_t       rank = 0;
    for (size_t j = 0; j < Q; ++j) {
        if (out.s[j] > sv_tol) {
            for (size_t i = 0; i < P; ++i) {
                out.U(i, j) = W(i, j) / out.s[j];
            }
            ++rank;
        }
    }

    // Complete U to a full P×P unitary: repeatedly add the standard basis vector
    // with the largest residual against the columns chosen so far. The leading
    // `rank` columns (the reliable singular vectors) are kept fixed.
    for (size_t filled = rank; filled < P; ++filled) {
        damp::array<T, P> best_v{};
        real_t            best_norm = real_t{0};
        for (size_t k = 0; k < P; ++k) {
            damp::array<T, P> v{};
            v[k] = T{1};
            for (size_t j = 0; j < filled; ++j) {
                T proj = T{0};
                for (size_t i = 0; i < P; ++i) {
                    proj += damp::conj(out.U(i, j)) * v[i];
                }
                for (size_t i = 0; i < P; ++i) {
                    v[i] -= proj * out.U(i, j);
                }
            }
            real_t nrm_sq = real_t{0};
            for (size_t i = 0; i < P; ++i) {
                const real_t a = damp::abs(v[i]);
                nrm_sq += a * a;
            }
            const real_t nrm = damp::sqrt(nrm_sq);
            if (nrm > best_norm) {
                best_norm = nrm;
                best_v = v;
            }
        }
        for (size_t i = 0; i < P; ++i) {
            out.U(i, filled) = best_v[i] / best_norm;
        }
    }

    return out;
}

} // namespace detail

/**
 * @brief Result of a full singular value decomposition A = U·Σ·Vᴴ.
 *
 * U (M×M) and V (N×N) are unitary; singular_values holds the min(M,N) diagonal
 * entries of Σ in non-increasing order. For real T, "unitary" means orthogonal
 * and Vᴴ = Vᵀ.
 */
template<size_t M, size_t N, typename T>
struct SVDResult {
    Matrix<M, M, T>                                singular_U{};      ///< Left singular vectors.
    damp::array<scalar_type_t<T>, (M < N ? M : N)> singular_values{}; ///< Descending, ≥ 0.
    Matrix<N, N, T>                                singular_V{};      ///< Right singular vectors.
    bool                                           converged = false;
};

/**
 * @brief Full singular value decomposition A = U·Σ·Vᴴ (one-sided Jacobi).
 *
 * @tparam M Rows, @tparam N Columns, @tparam T float/double or damp::complex.
 * @param A Input matrix (any shape).
 * @return SVDResult with unitary U (M×M), descending singular values, unitary
 *         V (N×N). Wide inputs (M < N) are handled by decomposing Aᴴ internally.
 *
 * @note Compare with MATLAB®'s [U, S, V] = svd(A).
 */
template<size_t M, size_t N, typename T>
[[nodiscard]] constexpr SVDResult<M, N, T> svd(const Matrix<M, N, T>& A) {
    SVDResult<M, N, T> out;
    if constexpr (M >= N) {
        const auto r = detail::jacobi_svd_tall<M, N, T>(A);
        out.singular_U = r.U;
        out.singular_V = r.V;
        out.singular_values = r.s;
        out.converged = r.converged;
    } else {
        // A = (Aᴴ)ᴴ = (Uᵣ·Σ·Vᵣᴴ)ᴴ = Vᵣ·Σ·Uᵣᴴ ⇒ swap the factors.
        const auto r = detail::jacobi_svd_tall<N, M, T>(A.conjugate_transpose());
        out.singular_U = r.V;
        out.singular_V = r.U;
        out.singular_values = r.s;
        out.converged = r.converged;
    }
    return out;
}

/**
 * @brief Numerical rank from a precomputed SVD result.
 * @param result SVD of the matrix.
 * @param rel_tol Singular values ≤ rel_tol·σ_max are treated as zero.
 */
template<size_t M, size_t N, typename T>
[[nodiscard]] constexpr size_t rank_from_svd(
    const SVDResult<M, N, T>& result,
    scalar_type_t<T>          rel_tol = default_tol<T>()
) {
    constexpr size_t       K = (M < N) ? M : N;
    const scalar_type_t<T> smax = (K > 0) ? result.singular_values[0] : scalar_type_t<T>{0};
    const scalar_type_t<T> tol = smax * rel_tol;
    size_t                 r = 0;
    for (size_t k = 0; k < K; ++k) {
        if (result.singular_values[k] > tol) {
            ++r;
        }
    }
    return r;
}

/**
 * @brief Moore–Penrose pseudoinverse A⁺ via SVD.
 *
 * A⁺ = V·Σ⁺·Uᴴ, where Σ⁺ inverts each singular value above the tolerance and
 * zeros the rest. Gives the least-squares minimum-norm solution: for an
 * overdetermined system x = A⁺·b minimizes ‖A·x − b‖; for an underdetermined
 * one it returns the smallest-norm solution. Reduces to A⁻¹ when A is square and
 * nonsingular.
 *
 * @param A Input matrix (any shape).
 * @param rel_tol Singular values ≤ rel_tol·σ_max are treated as zero.
 * @return Pseudoinverse (N×M).
 *
 * @note Compare with MATLAB®'s pinv(A).
 */
template<size_t M, size_t N, typename T>
[[nodiscard]] constexpr Matrix<N, M, T> pseudo_inverse(
    const Matrix<M, N, T>& A,
    scalar_type_t<T>       rel_tol = default_tol<T>()
) {
    using real_t = scalar_type_t<T>;
    constexpr size_t K = (M < N) ? M : N;

    const auto   s = svd(A);
    const real_t smax = (K > 0) ? s.singular_values[0] : real_t{0};
    const real_t tol = smax * rel_tol;

    // Sp_Uh = Σ⁺·Uᴴ : an N×M matrix whose row k is (1/σ_k)·(Uᴴ row k) for the
    // first K rows with σ_k > tol, and zero elsewhere.
    const Matrix<M, M, T> Uh = s.singular_U.conjugate_transpose();
    Matrix<N, M, T>       Sp_Uh = Matrix<N, M, T>::zeros();
    for (size_t k = 0; k < K; ++k) {
        if (s.singular_values[k] > tol) {
            const T inv = T{1} / T{s.singular_values[k]};
            for (size_t c = 0; c < M; ++c) {
                Sp_Uh(k, c) = inv * Uh(k, c);
            }
        }
    }
    return Matrix<N, M, T>(s.singular_V * Sp_Uh);
}

/// Orthonormal basis for the null space (kernel) of a matrix.
template<size_t N, typename T>
struct NullSpace {
    Matrix<N, N, T> vectors{}; ///< Kernel basis = the last `dim` columns.
    size_t          dim = 0;   ///< Nullity (number of basis vectors).
};

/**
 * @brief Orthonormal basis for the null space {x : A·x = 0} via SVD.
 *
 * The right singular vectors with zero singular value span the kernel. Because
 * the singular values are sorted descending, those vectors are the trailing
 * columns of V — so the basis is the last `dim` columns of the returned matrix.
 *
 * @param A Input matrix (any shape).
 * @param rel_tol Singular values ≤ rel_tol·σ_max are treated as zero.
 * @return NullSpace whose last `dim` columns form an orthonormal kernel basis.
 *
 * @note Compare with MATLAB®'s null(A) (which returns those columns directly).
 */
template<size_t M, size_t N, typename T>
[[nodiscard]] constexpr NullSpace<N, T> null_space(
    const Matrix<M, N, T>& A,
    scalar_type_t<T>       rel_tol = default_tol<T>()
) {
    const auto      s = svd(A);
    const size_t    r = rank_from_svd(s, rel_tol);
    NullSpace<N, T> out;
    out.vectors = s.singular_V;
    out.dim = N - r;
    return out;
}

} // namespace mat

} // namespace damp
