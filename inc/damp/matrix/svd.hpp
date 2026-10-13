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
#include <limits>

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
    Matrix<P, P, T>                  U{};               ///< Left singular vectors (P×P unitary).
    damp::array<scalar_type_t<T>, Q> s{};               ///< Singular values, descending, ≥ 0.
    Matrix<Q, Q, T>                  V{};               ///< Right singular vectors (Q×Q unitary).
    bool                             converged = false; ///< A sweep finished with no rotation.
};

/// Default relative rank threshold max(M, N)·eps (MATLAB® rank / pinv / null).
template<size_t M, size_t N, typename T>
[[nodiscard]] constexpr scalar_type_t<T> svd_default_rel_tol() {
    return static_cast<scalar_type_t<T>>(M > N ? M : N) * damp::eps<T>();
}

/**
 * @brief One-sided Jacobi SVD of a tall/square matrix A (P×Q, P ≥ Q).
 *
 * Orthogonalizes the columns of A by a sequence of two-sided rotations; the
 * orthogonalized column norms are the singular values, the normalized columns
 * the (economy) left vectors, and the accumulated rotations the right vectors.
 * The left factor is then completed to a full P×P unitary by Gram–Schmidt so
 * the result is a genuine full SVD even when A is rank-deficient.
 *
 * A pair of columns is rotated while its cosine |⟨a_p, a_q⟩| / (‖a_p‖‖a_q‖)
 * exceeds P·eps (Demmel–Veselić), so every nonzero column — however small —
 * ends orthogonal to working precision and becomes its own left vector.
 *
 * Nothing is squared. Each column is held as a unit vector times its norm (a scaled norm,
 * damp::detail::scaled_norm2), the cosine is the inner product of the unit vectors, and the
 * Jacobi angle comes from the norm ratio ρ = min/max and the cosine c: @f$ |\zeta| = (1-\rho^2)/(2\rho|c|) @f$,
 * with @f$ t = 1/(|\zeta| + \sqrt{\zeta^2+1}) @f$ for |ζ| ≤ 1 and
 * @f$ t = u/(1 + \sqrt{1+u^2}) @f$, @f$ u = 1/|\zeta| @f$, beyond. A Gram sum ‖a‖² flushes
 * a column 10⁻²⁰ below the largest to zero in float, which cost a graded matrix its small
 * singular values. A is rescaled by its largest entry only when every entry is below one
 * or the column norms could overflow, so graded entries stay where they are; a column
 * whose norm falls below min/eps (absolute) is roundoff residue and is reported as zero.
 *
 * @see "Jacobi's method is more accurate than QR" (Demmel & Veselić, 1992), SIAM J. Matrix Anal. Appl. 13(4)
 * @see "New fast and accurate Jacobi SVD algorithm I" (Drmač & Veselić, 2008), SIAM J. Matrix Anal. Appl. 29(4)
 */
template<size_t P, size_t Q, typename T>
[[nodiscard]] constexpr TallSVD<P, Q, T> jacobi_svd_tall(const Matrix<P, Q, T>& A) {
    static_assert(P >= Q, "jacobi_svd_tall requires P >= Q (transpose wide inputs)");
    using real_t = scalar_type_t<T>;
    const real_t tol = static_cast<real_t>(P) * damp::eps<T>();
    // Below min/eps a column norm has lost relative precision: its direction is noise.
    const real_t tiny = std::numeric_limits<real_t>::min() / damp::eps<T>();
    // A column norm can reach √P·max|a|; past this, rescale before taking norms.
    const real_t big = std::numeric_limits<real_t>::max() / (real_t{2} * static_cast<real_t>(P));

    real_t amax = real_t{0};
    for (size_t i = 0; i < P; ++i) {
        for (size_t j = 0; j < Q; ++j) {
            amax = damp::max(amax, damp::abs(A(i, j)));
        }
    }
    // Rescale only when needed: up when every entry is below one (none is pushed toward
    // underflow), down when a norm could overflow. Otherwise a graded A keeps its entries.
    const bool      rescale = (amax > real_t{0}) && (amax < real_t{1} || amax > big);
    const real_t    inv_amax = rescale ? real_t{1} / amax : real_t{1};
    Matrix<P, Q, T> W = A * T{inv_amax};
    Matrix<Q, Q, T> V = Matrix<Q, Q, T>::identity();

    // Column j of the working matrix is W(:, j)·d[j] with W(:, j) a unit vector (Drmač–Veselić):
    // cosines are plain inner products of unit vectors, with no reciprocal norms that fast-math
    // could hoist out of the sum as one overflowing or flushing product.
    damp::array<real_t, Q> d{};
    const auto             renormalize = [&](size_t j, real_t scale) {
        const real_t n = damp::detail::scaled_norm2<real_t>(P, [&](size_t i) { return damp::abs(W(i, j)); });
        const real_t dj = scale * n;
        if (!(n > real_t{0}) || dj < tiny) {
            for (size_t i = 0; i < P; ++i) {
                W(i, j) = T{0};
            }
            d[j] = real_t{0};
            return;
        }
        const real_t inv = real_t{1} / n;
        for (size_t i = 0; i < P; ++i) {
            W(i, j) = W(i, j) * T{inv};
        }
        d[j] = dj;
    };
    for (size_t j = 0; j < Q; ++j) {
        renormalize(j, real_t{1});
    }

    // Sweep over column pairs, rotating each pair to mutual orthogonality. A
    // single column has no pairs and is already orthogonal; skipping the sweep
    // at compile time also keeps its out-of-range W(i, q) out of the instance.
    constexpr size_t max_sweeps = 60;
    bool             converged = (Q < 2);
    if constexpr (Q > 1) {
        for (size_t sweep = 0; sweep < max_sweeps; ++sweep) {
            bool rotated = false;
            for (size_t p = 0; p < Q; ++p) {
                for (size_t q = p + 1; q < Q; ++q) {
                    // A zero column (or roundoff residue below min/eps) has no direction.
                    if (d[p] == real_t{0} || d[q] == real_t{0}) {
                        continue;
                    }
                    T cosine = T{0};
                    for (size_t i = 0; i < P; ++i) {
                        cosine += damp::conj(W(i, p)) * W(i, q);
                    }
                    const real_t gmag = damp::abs(cosine);
                    if (gmag <= tol) {
                        continue; // already orthogonal
                    }
                    const bool   p_big = d[p] > d[q];
                    const real_t rho = p_big ? (d[q] / d[p]) : (d[p] / d[q]);
                    if (rho == real_t{0}) {
                        continue; // norm ratio below the range: no representable rotation moves either column
                    }
                    rotated = true;

                    // Real Jacobi angle that zeros the off-diagonal of the phase-aligned
                    // Gram block [[d_p², d_p d_q |c|], [d_p d_q |c|, d_q²]]:
                    // ζ = (d_q² − d_p²)/(2 d_p d_q |c|), smaller root t = −sign(ζ)/(|ζ| + √(ζ²+1)),
                    // written in ρ = min/max so nothing is squared or overflows.
                    const real_t zsign = p_big ? real_t{-1} : real_t{1};
                    const real_t gap = (real_t{1} - rho) * (real_t{1} + rho); // 1 − ρ²
                    const real_t two_rc = real_t{2} * rho * gmag;
                    real_t       t_abs = real_t{0};
                    if (gap <= two_rc) {
                        const real_t az = gap / two_rc; // |ζ| ≤ 1
                        t_abs = real_t{1} / (az + damp::sqrt((az * az) + real_t{1}));
                    } else {
                        const real_t u = two_rc / gap; // 1/|ζ| < 1
                        t_abs = u / (real_t{1} + damp::sqrt(real_t{1} + (u * u)));
                    }
                    const real_t t = -zsign * t_abs;
                    const real_t c = real_t{1} / damp::sqrt(real_t{1} + (t * t));
                    const real_t sn = c * t;
                    // ph = e^{-i·arg(cosine)} folds the complex phase into the rotation;
                    // for real T it is simply ±1.
                    const T ph = damp::conj(cosine) / gmag;

                    // w_p ← c·w_p + ph·s·w_q, w_q ← −s·w_p + ph·c·w_q on w = W·d: in unit-vector
                    // form the cross terms carry s·d_q/d_p and s·d_p/d_q, both bounded by |c|·…
                    // since |s| ≈ ρ|c| whenever ρ is small.
                    const real_t sq_over_p = p_big ? (sn * rho) : (sn / rho); // s·d_q/d_p
                    const real_t sp_over_q = p_big ? (sn / rho) : (sn * rho); // s·d_p/d_q
                    const T      phs = ph * T{sq_over_p};
                    const T      phc = ph * T{c};
                    for (size_t i = 0; i < P; ++i) {
                        const T wip = W(i, p);
                        const T wiq = W(i, q);
                        W(i, p) = (c * wip) + (phs * wiq);
                        W(i, q) = (-sp_over_q * wip) + (phc * wiq);
                    }
                    renormalize(p, d[p]);
                    renormalize(q, d[q]);
                    for (size_t i = 0; i < Q; ++i) {
                        const T vip = V(i, p);
                        const T viq = V(i, q);
                        V(i, p) = (c * vip) + ((ph * T{sn}) * viq);
                        V(i, q) = (-sn * vip) + (phc * viq);
                    }
                }
            }
            if (!rotated) {
                converged = true;
                break;
            }
        }
    }

    // Singular values = column norms; economy U = the unit columns.
    TallSVD<P, Q, T> out;
    out.V = V;
    out.converged = converged;
    for (size_t j = 0; j < Q; ++j) {
        out.s[j] = d[j];
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

    // Every nonzero column is orthogonal to the others in the relative sense,
    // so it normalizes into a left singular vector (sorted: zeros trail).
    // Residue columns (norm below min/eps) were zeroed by renormalize and are not
    // orthogonal to the others; the completion below supplies their left vectors.
    size_t rank = 0;
    for (size_t j = 0; j < Q; ++j) {
        if (out.s[j] > real_t{0}) {
            for (size_t i = 0; i < P; ++i) {
                out.U(i, j) = W(i, j);
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

    if (rescale) {
        for (size_t j = 0; j < Q; ++j) {
            out.s[j] *= amax;
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
    bool                                           converged = false; ///< Jacobi sweeps converged.
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
 * @param rel_tol Singular values ≤ rel_tol·σ_max are treated as zero (default max(M, N)·eps, as MATLAB®).
 */
template<size_t M, size_t N, typename T>
[[nodiscard]] constexpr size_t rank_from_svd(
    const SVDResult<M, N, T>& result,
    scalar_type_t<T>          rel_tol = detail::svd_default_rel_tol<M, N, T>()
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
 * @param rel_tol Singular values ≤ rel_tol·σ_max are treated as zero (default max(M, N)·eps, as MATLAB®).
 * @return Pseudoinverse (N×M).
 *
 * @note Compare with MATLAB®'s pinv(A).
 */
template<size_t M, size_t N, typename T>
[[nodiscard]] constexpr Matrix<N, M, T> pseudo_inverse(
    const Matrix<M, N, T>& A,
    scalar_type_t<T>       rel_tol = detail::svd_default_rel_tol<M, N, T>()
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
 * @param rel_tol Singular values ≤ rel_tol·σ_max are treated as zero (default max(M, N)·eps, as MATLAB®).
 * @return NullSpace whose last `dim` columns form an orthonormal kernel basis.
 *
 * @note Compare with MATLAB®'s null(A) (which returns those columns directly).
 */
template<size_t M, size_t N, typename T>
[[nodiscard]] constexpr NullSpace<N, T> null_space(
    const Matrix<M, N, T>& A,
    scalar_type_t<T>       rel_tol = detail::svd_default_rel_tol<M, N, T>()
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
