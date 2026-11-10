// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file decomposition.hpp
 * @brief Matrix factorizations: Cholesky, LDLᵀ (rook pivoting), LU (partial pivoting), and thin
 *        and full QR (Householder), plus the symmetry/Hermitian predicate used to
 *        choose between them. Linear-system solvers built on these live in
 *        solve.hpp; the eigensolver in eigen.hpp.
 */


#include "core.hpp"
#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "matrix_traits.hpp"

namespace damp {

namespace mat {

namespace detail {

/**
 * @brief True if @p v is finite (both parts of a damp::complex), by the bit-pattern test.
 *
 * The parts are tested directly: damp::abs of a damp::complex with a NaN part and a zero
 * part returns 0.
 */
template<typename T>
[[nodiscard]] constexpr bool is_finite_entry(const T& v) {
    if constexpr (damp::is_floating_point_v<T>) {
        return damp::isfinite(v);
    } else {
        return damp::isfinite(v.real()) && damp::isfinite(v.imag());
    }
}

/**
 * @brief True if each pair aᵢⱼ, conj(aⱼᵢ) agrees to @p rel times max(|aᵢⱼ|, |aⱼᵢ|, √|aᵢᵢ aⱼⱼ|); false for a non-finite entry.
 *
 * A NaN or ±Inf fails every comparison against the bound (and makes the bound itself
 * infinite), so finiteness is tested first: a matrix with a non-finite strict-upper entry
 * is not Hermitian, and cholesky(), which reads only the lower triangle, never sees it.
 */
template<size_t N, typename T>
[[nodiscard]] constexpr bool is_hermitian_within(const Matrix<N, N, T>& A, scalar_type_t<T> rel) {
    using real_t = scalar_type_t<T>;
    if (!mat::all_of(A, [](const T& v) { return is_finite_entry(v); })) {
        return false;
    }
    for (size_t i = 0; i < N; ++i) {
        if constexpr (!damp::is_floating_point_v<T>) {
            // Hermitian requires a real diagonal
            if (damp::abs(damp::imag(A(i, i))) > rel * damp::abs(A(i, i))) {
                return false;
            }
        }
        for (size_t j = i + 1; j < N; ++j) {
            // √(|aᵢᵢ||aⱼⱼ|) as hi·√(lo/hi): the product overflows past √max (1.8e19 in
            // float), and -ffast-math folds √x·√y back into √(xy).
            const real_t dhi = damp::max(damp::abs(A(i, i)), damp::abs(A(j, j)));
            const real_t dlo = damp::min(damp::abs(A(i, i)), damp::abs(A(j, j)));
            const real_t gmean = (dhi > real_t{0}) ? dhi * damp::sqrt(dlo / dhi) : real_t{0};
            const real_t scale = damp::max({damp::abs(A(i, j)), damp::abs(A(j, i)), gmean});
            if (damp::abs(A(i, j) - damp::conj(A(j, i))) > rel * scale) {
                return false;
            }
        }
    }
    return true;
}

} // namespace detail

/**
 * @brief True if A equals its transpose (real) or conjugate transpose (complex)
 *
 * Each pair may differ by N·64·eps times its own scale
 * max(|aᵢⱼ|, |aⱼᵢ|, √|aᵢᵢ aⱼⱼ|), so the test is invariant to scaling A and to
 * graded (diagonally scaled) matrices. A matrix with a NaN or ±Inf entry is not
 * symmetric.
 */
template<size_t N, typename T>
[[nodiscard]] constexpr bool is_symmetric_or_hermitian(const Matrix<N, N, T>& A) {
    return detail::is_hermitian_within(A, static_cast<scalar_type_t<T>>(N) * scalar_type_t<T>{64} * eps<T>());
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
 * A NaN or ±Inf in the lower triangle reaches a pivot (an off-diagonal one through lᵢⱼ and
 * |lᵢⱼ|²), and a pivot that is not finite is refused by its bit pattern (damp::isfinite)
 * before the floor is read. The floor test alone does not refuse it: under
 * -ffinite-math-only GCC folds the negated comparison into ≤, which a NaN fails, so a NaN
 * pivot factored and cholesky_solve() returned NaN.
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
                if constexpr (damp::is_floating_point_v<T>) {
                    diag_val = sum;
                } else {
                    diag_val = damp::real(sum);
                }
                // Finiteness first: -ffinite-math-only folds !(d > floor) into d ≤ floor, false for a NaN d.
                if (!detail::is_finite_entry(sum) || !(diag_val > real_t{16} * static_cast<real_t>(i + 1) * eps<T>() * mag)) {
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

/**
 * @brief LDLᴴ factors with symmetric pivoting: P·A·Pᵀ = L·D·Lᴴ, D block diagonal with 1×1 and 2×2 blocks
 *
 * D is Hermitian with real diagonal d. A 2×2 block on rows k, k+1 keeps its off-diagonal
 * D(k+1, k) in e[k], the block's largest entry in magnitude and never zero; every other
 * entry of e is zero. (P·A·Pᵀ)(i, j) = A(perm[i], perm[j]), the convention of
 * lu_decomposition()'s pivot vector.
 */
template<size_t N, typename T>
struct [[nodiscard]] LDLT {
    Matrix<N, N, T>                  L{};    ///< Unit lower-triangular factor, |lᵢⱼ| ≤ 1/(1 − α) ≈ 2.78
    damp::array<scalar_type_t<T>, N> d{};    ///< Diagonal of D (real)
    damp::array<T, N>                e{};    ///< Subdiagonal of D: e[k] = D(k+1, k), nonzero only on the first row of a 2×2 block
    damp::array<size_t, N>           perm{}; ///< Row i of P·A·Pᵀ is row perm[i] of A
};

namespace detail {

/** @brief Bunch and Kaufman's α = (1 + √17)/8, which equalizes the growth bounds of a 1×1 and a 2×2 step. */
template<typename R>
inline constexpr R ldlt_alpha = static_cast<R>(0.64038820320220756872767623199676);

/** @brief Real part of a diagonal entry (the entry itself for real T). */
template<typename T>
[[nodiscard]] constexpr scalar_type_t<T> real_entry(const T& x) {
    if constexpr (damp::is_floating_point_v<T>) {
        return x;
    } else {
        return damp::real(x);
    }
}

/**
 * @brief A 2×2 block E = [[a, w̄], [w, c]] of ldlt()'s D, with |a|, |c| < α|w|, solved scaled by |w|
 *
 * With â = a/|w| and ĉ = c/|w|, det E = |w|²·(âĉ − 1) and |âĉ| < α² ≈ 0.41, so δ = âĉ − 1 lies
 * in (−1.41, −0.59): the determinant never cancels, and |w|², which overflows past 1.8·10¹⁹
 * and flushes below 10⁻¹⁹ in float, is never formed. κ(E) ≤ (1 + α)/(1 − α) ≈ 4.6.
 */
template<typename T>
struct LDLTBlock {
    using real_t = scalar_type_t<T>;
    real_t m{};     ///< |w|
    real_t a{};     ///< â = a/|w|
    real_t c{};     ///< ĉ = c/|w|
    real_t delta{}; ///< δ = âĉ − 1 = det E/|w|²
    T      u{};     ///< w/|w| (±1 for real T)

    /** @brief The block with diagonal a0, c0 and subdiagonal w. */
    [[nodiscard]] static constexpr LDLTBlock of(real_t a0, real_t c0, const T& w) {
        LDLTBlock E{};
        E.m = damp::abs(w);
        E.a = a0 / E.m;
        E.c = c0 / E.m;
        E.delta = (E.a * E.c) - real_t{1};
        E.u = w / E.m;
        return E;
    }

    /** @brief (y₁, y₂) ← E⁻¹(y₁, y₂) = ((ĉy₁ − ūy₂), (ây₂ − uy₁)) / (|w|δ) */
    constexpr void solve(T& y1, T& y2) const {
        const T s1 = y1 / m;
        const T s2 = y2 / m;
        y1 = ((c * s1) - (damp::conj(u) * s2)) / delta;
        y2 = ((a * s2) - (u * s1)) / delta;
    }
};

/**
 * @brief ldlt()'s state: the active matrix S (lower triangle), the size M of the terms behind
 *        each entry of S, the factors so far, and k, the number of rows eliminated.
 */
template<size_t N, typename T>
struct LDLTWork {
    using real_t = scalar_type_t<T>;

    /** @brief Largest non-roundoff off-diagonal magnitude of a column and its row (N if none). */
    struct ColMax {
        real_t w{};           ///< ωᵢ
        size_t row = N;       ///< Its row, or N when every off-diagonal entry is roundoff
        bool   finite = true; ///< False when an entry of the column (or its M) is not finite
    };

    LDLT<N, T>           f{};   ///< Factors so far
    Matrix<N, N, T>      S{};   ///< Active matrix, rows and columns k…N−1 of the lower triangle
    Matrix<N, N, real_t> M{};   ///< |A| + Σ|l||s|: the size of the terms behind each entry of S
    size_t               k = 0; ///< Rows eliminated

    /** @brief Loads the lower triangle of A; false for a non-finite entry. */
    [[nodiscard]] constexpr bool load(const Matrix<N, N, T>& A) {
        f.L = Matrix<N, N, T>::identity();
        for (size_t i = 0; i < N; ++i) {
            f.perm[i] = i;
            for (size_t j = 0; j <= i; ++j) {
                if (!is_finite_entry(A(i, j))) {
                    return false;
                }
                S(i, j) = (i == j) ? T{real_entry(A(i, i))} : A(i, j);
                M(i, j) = damp::abs(S(i, j));
            }
        }
        return true;
    }

    /** @brief Entry (i, j), i ≥ j, is cancellation roundoff: it does not clear 16·(k+1)·eps of its k+1 terms. */
    [[nodiscard]] constexpr bool roundoff(size_t i, size_t j) const {
        return !(damp::abs(S(i, j)) > real_t{16} * static_cast<real_t>(k + 1) * eps<T>() * M(i, j));
    }

    /** @brief |S(i, i)|, or 0 when it is roundoff. */
    [[nodiscard]] constexpr real_t diag(size_t i) const { return roundoff(i, i) ? real_t{0} : damp::abs(S(i, i)); }

    /** @brief ωᵢ over the active rows of column i, read from the lower triangle; roundoff entries count as zero. */
    [[nodiscard]] constexpr ColMax col_max(size_t i) const {
        ColMax out{};
        for (size_t x = k; x < N; ++x) {
            const size_t r = (x > i) ? x : i;
            const size_t c = (x > i) ? i : x;
            if (!is_finite_entry(S(r, c)) || !damp::isfinite(M(r, c))) {
                out.finite = false;
                return out;
            }
            if (x != i && !roundoff(r, c) && (out.row == N || damp::abs(S(r, c)) > out.w)) {
                out.w = damp::abs(S(r, c));
                out.row = x;
            }
        }
        return out;
    }

    /**
     * @brief Rook search from column k: {i, N} for a 1×1 pivot on i, {i, r} for a 2×2 pivot on
     *        (i, r), {N, N} to refuse (a column that is roundoff throughout, or a non-finite entry).
     */
    [[nodiscard]] constexpr damp::array<size_t, 2> rook_pivot() const {
        constexpr real_t alpha = ldlt_alpha<real_t>;
        size_t           i = k;
        ColMax           wi = col_max(i);
        if (!wi.finite) {
            return {N, N};
        }
        if (wi.row == N) {
            // Off-diagonal roundoff throughout: the diagonal pivots, or the column is lost.
            return {roundoff(i, i) ? N : i, N};
        }
        if (diag(i) >= alpha * wi.w) {
            return {i, N};
        }
        // ω rises strictly from column to column, so none repeats: at most N − k − 1 passes.
        for (size_t pass = 0; pass < N; ++pass) {
            const size_t r = wi.row;
            const ColMax wr = col_max(r);
            if (!wr.finite) {
                return {N, N};
            }
            if (diag(r) >= alpha * wr.w) {
                return {r, N};
            }
            if (!(wr.w > wi.w)) {
                return {i, r}; // aᵣᵢ is the largest off-diagonal entry of both columns
            }
            i = r;
            wi = wr;
        }
        return {N, N};
    }

    /** @brief Symmetric interchange of active indices q ≤ p: S and M in place, the rows of L so far, perm. */
    constexpr void swap(size_t q, size_t p) {
        if (q == p) {
            return;
        }
        for (size_t x = k; x < q; ++x) {
            damp::swap(S(q, x), S(p, x));
            damp::swap(M(q, x), M(p, x));
        }
        damp::swap(S(q, q), S(p, p));
        damp::swap(M(q, q), M(p, p));
        for (size_t x = q + 1; x < p; ++x) {
            const T s = S(x, q);
            S(x, q) = damp::conj(S(p, x));
            S(p, x) = damp::conj(s);
            damp::swap(M(x, q), M(p, x));
        }
        S(p, q) = damp::conj(S(p, q));
        for (size_t x = p + 1; x < N; ++x) {
            damp::swap(S(x, q), S(x, p));
            damp::swap(M(x, q), M(x, p));
        }
        for (size_t x = 0; x < k; ++x) {
            damp::swap(f.L(q, x), f.L(p, x));
        }
        damp::swap(f.perm[q], f.perm[p]);
    }

    /** @brief S(i, j) −= Σ L(i, c)·conj(S(j, c)) over the pivot columns c ∈ [k, k + s), for i ≥ j ≥ k + s. */
    constexpr void update(size_t s) {
        for (size_t j = k + s; j < N; ++j) {
            for (size_t i = j; i < N; ++i) {
                for (size_t c = k; c < k + s; ++c) {
                    S(i, j) -= f.L(i, c) * damp::conj(S(j, c));
                    M(i, j) += damp::abs(f.L(i, c)) * damp::abs(S(j, c));
                }
            }
            S(j, j) = T{real_entry(S(j, j))};
        }
        k += s;
    }

    /** @brief Eliminates the 1×1 pivot at k. */
    constexpr void eliminate_1x1() {
        const real_t p = real_entry(S(k, k));
        f.d[k] = p;
        for (size_t i = k + 1; i < N; ++i) {
            // A roundoff entry the pivot test did not bound is zero within its rounding bound.
            if (roundoff(i, k) && ldlt_alpha<real_t> * damp::abs(S(i, k)) > damp::abs(p)) {
                S(i, k) = T{0};
            }
            f.L(i, k) = S(i, k) / p;
        }
        update(1);
    }

    /** @brief Eliminates the 2×2 pivot on rows k, k+1: Lᵢᴴ = E⁻¹ (S(i, k), S(i, k+1))ᴴ. */
    constexpr void eliminate_2x2() {
        // The search saw a roundoff diagonal as zero; storing it so keeps |a|, |c| < α|w|.
        for (size_t c = k; c < k + 2; ++c) {
            if (roundoff(c, c)) {
                S(c, c) = T{0};
            }
        }
        f.d[k] = real_entry(S(k, k));
        f.d[k + 1] = real_entry(S(k + 1, k + 1));
        f.e[k] = S(k + 1, k);
        const LDLTBlock<T> E = LDLTBlock<T>::of(f.d[k], f.d[k + 1], f.e[k]);
        for (size_t i = k + 2; i < N; ++i) {
            for (size_t c = k; c < k + 2; ++c) {
                if (roundoff(i, c) && damp::abs(S(i, c)) > E.m) {
                    S(i, c) = T{0};
                }
            }
            T l1 = damp::conj(S(i, k));
            T l2 = damp::conj(S(i, k + 1));
            E.solve(l1, l2);
            f.L(i, k) = damp::conj(l1);
            f.L(i, k + 1) = damp::conj(l2);
        }
        update(2);
    }
};

} // namespace detail

/**
 * @brief Symmetric indefinite LDLᴴ with rook pivoting: P·A·Pᵀ = L·D·Lᴴ, D block diagonal (1×1 and 2×2)
 *
 * Factors any symmetric (Hermitian) matrix that is nonsingular to working precision, definite
 * or not, without square roots. Each step takes a 1×1 or a 2×2 pivot by the rook search of
 * Ashcraft, Grimes and Lewis (their bounded Bunch–Kaufman) with Bunch and Kaufman's
 * @f$ \alpha = (1 + \sqrt{17})/8 \approx 0.64 @f$. With @f$ \omega_i @f$ the largest
 * off-diagonal magnitude in column i of the active matrix, the search starts at the next
 * column i. If @f$ |a_{ii}| \ge \alpha\,\omega_i @f$, i is a 1×1 pivot; otherwise the search
 * moves to the row r of @f$ \omega_i @f$ and stops on the 2×2 pivot (i, r) once
 * @f$ a_{ri} @f$ is the largest off-diagonal entry of column r too. The ω it meets rise
 * strictly, so no column is visited twice: a step scans at most N − k columns (under two on
 * average), and the work stays O(N³).
 *
 * The pivot test alone bounds both factors, for every matrix:
 * - A 1×1 pivot has multipliers @f$ |a_{ri}|/|a_{ii}| \le 1/\alpha @f$ and grows the active
 *   entries by at most @f$ 1 + 1/\alpha \approx 2.56 @f$.
 * - A 2×2 pivot @f$ E = \begin{bmatrix} a & \bar w \\ w & c \end{bmatrix} @f$ has
 *   @f$ |a|, |c| < \alpha|w| @f$ with |w| the largest entry of both its columns, so
 *   @f$ |\det E| > (1 - \alpha^2)|w|^2 @f$, its multipliers are at most
 *   @f$ (1 + \alpha)/(1 - \alpha^2) = 1/(1 - \alpha) \approx 2.78 @f$, and it grows the active
 *   entries by at most @f$ 1 + 2/(1 - \alpha) = (1 + 1/\alpha)^2 @f$, the balance that fixes α.
 *
 * The computed factors satisfy @f$ PAP^\top + \Delta A = LDL^H @f$ with
 * @f$ |\Delta A| \le p(n)\,\varepsilon\,(|A| + |L||D||L^H|) @f$ (Higham, §11.1). With |L| at
 * most 2.78, each entry of the second term is at most @f$ 2 \cdot 2.78^2\, n \max|D| @f$, so the
 * factors and ldlt_solve() are backward stable up to the growth of D alone, the guarantee of
 * LU with partial pivoting. Bunch–Kaufman partial pivoting (LAPACK's xSYTRF) stops the search
 * at the second column: D grows no faster, but |L| is unbounded. On
 * @f$ \begin{bmatrix} \delta^2 & \delta & 0 \\ \delta & 0 & 1 \\ 0 & 1 & 1 \end{bmatrix} @f$
 * it takes the 1×1 pivot δ² and l = 1/δ, and Higham (1997) bounds its |L||D||Lᵀ| in norm only,
 * through small entries of D offsetting large ones of L. On random indefinite 7×7 matrices
 * in float, the rook search kept max |L| at 1.8 against Bunch–Kaufman's 8.0 and
 * ‖|L||D||Lᵀ|‖/‖A‖ at 6.0 against 14, scanning 1.7 columns per step against 1.5.
 *
 * An entry of the active matrix that does not exceed 16·(k+1)·eps of the magnitude of the
 * k+1 terms that formed it is cancellation roundoff (the floor of cholesky() and
 * lu_decomposition()): the search counts it as zero, and one too large for the bound on L is
 * eliminated as zero, a perturbation within its rounding bound. Returns damp::nullopt when a
 * column of the active matrix is roundoff throughout (A is singular to working precision), or
 * for a non-finite entry of the lower triangle or of the active matrix (overflow); finiteness
 * is tested on the bit pattern, so a NaN is refused under -ffast-math. The tests compare
 * entries to the terms that formed them and so are invariant to scaling A; the pivot choice
 * compares absolute sizes, so a graded matrix pivots on its large-unit entries first.
 *
 * A 2×2 block is solved scaled by its off-diagonal (detail::LDLTBlock): neither its
 * determinant nor |w|² is formed. Every 2×2 block has det E < 0, so A has as many negative
 * eigenvalues as D has negative 1×1 pivots plus 2×2 blocks (Sylvester's law of inertia).
 *
 * Reads only the lower triangle and assumes the upper mirrors it.
 *
 * @note Compare with MATLAB®'s [L, D, P] = ldl(A), whose P'·A·P = L·D·L' makes its P the Pᵀ here.
 * @see "Accurate Symmetric Indefinite Linear Equation Solvers" (Ashcraft, Grimes & Lewis, 1998), §2–3
 * @see "Some Stable Methods for Calculating Inertia and Solving Symmetric Linear Systems" (Bunch & Kaufman, 1977)
 * @see "Accuracy and Stability of Numerical Algorithms" (Higham, 2002), §11.1
 * @see "Stability of the Diagonal Pivoting Method with Partial Pivoting" (Higham, 1997)
 * @see "Matrix Computations" (Golub & Van Loan, 2013), §4.4
 * @see ldlt_solve(), cholesky()
 *
 * @param A Symmetric (real) or Hermitian (complex) matrix
 * @return L, D (d, e) and P (perm), or damp::nullopt for a singular or non-finite A
 */
template<size_t N, typename T>
constexpr damp::optional<LDLT<N, T>> ldlt(const Matrix<N, N, T>& A) {
    detail::LDLTWork<N, T> w{};
    if (!w.load(A)) {
        return damp::nullopt;
    }
    // k rises by 1 or 2 per pass
    for (size_t pass = 0; pass < N && w.k < N; ++pass) {
        const auto piv = w.rook_pivot();
        if (piv[0] == N) {
            return damp::nullopt;
        }
        if (piv[1] == N) {
            w.swap(w.k, piv[0]);
            w.eliminate_1x1();
        } else {
            w.swap(w.k, damp::min(piv[0], piv[1]));
            w.swap(w.k + 1, damp::max(piv[0], piv[1]));
            w.eliminate_2x2();
        }
    }
    return w.f;
}

/**
 * @brief LU decomposition with partial pivoting
 *
 * Factors @f$ P A = L U @f$ where @f$ L @f$ is unit lower-triangular,
 * @f$ U @f$ is upper-triangular, and @f$ P @f$ is encoded by the pivot vector
 * (@c piv[i] is the original row moved into position @c i).
 *
 * Each entry of the active column at step i is measured against the magnitude
 * mᵢⱼ = |aᵢⱼ| + Σₖ|lᵢₖ||uₖⱼ| of the i+1 terms that formed it. An entry that does not exceed
 * 16·(i+1)·eps·mᵢⱼ (the rounding bound of those terms, the floor cholesky() uses) is
 * cancellation roundoff: it is not a pivot candidate, and when it is larger than the pivot
 * it is taken as zero. The pivot is the largest entry that is not roundoff, eliminated with
 * |lᵢⱼ| ≤ 1, and damp::nullopt is returned when it does not exceed N·64·eps·mᵢⱼ (the column
 * is lost to cancellation, or an entry that still carries digits is too small to pivot on),
 * for a non-finite entry of A, or for a pivot that the elimination overflowed to a non-finite
 * value.
 * The roundoff test is invariant to row and column scaling, so graded but well-posed systems
 * such as [[1, 10¹⁶], [0, 1]] still factor when their cancellation residues stay within the
 * 16·(i+1)·eps bound, and a row whose largest entry is cancellation roundoff does not hide a
 * genuine pivot below it. The pivot choice by absolute size is not: in a row-graded matrix,
 * an entry that carries a few digits and outranks a clean pivot in a smaller row is refused.
 * U is stored exactly upper-triangular.
 *
 * @note Compare with MATLAB®'s [L,U,P] = lu(A) (row-permutation form).
 * @see Golub & Van Loan, "Matrix Computations" (4th ed., 2013), §3.2, §3.3.1
 */
template<size_t N, typename T>
constexpr damp::optional<damp::tuple<Matrix<N, N, T>, Matrix<N, N, T>, damp::array<size_t, N>>>
lu_decomposition(const Matrix<N, N, T>& A) {
    using real_t = scalar_type_t<T>;
    constexpr real_t rel = static_cast<real_t>(N) * real_t{64} * eps<T>();
    if (!mat::all_of(A, [](const T& v) { return detail::is_finite_entry(v); })) {
        return damp::nullopt; // an infinite entry would read as roundoff of its infinite magnitude and be zeroed
    }

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

    // Entry (r, i) is cancellation roundoff when it does not clear the rounding bound of its
    // i+1 terms, 16·(i+1)·eps·mag(r, i).
    const auto roundoff = [&](size_t r, size_t i) {
        return !(damp::abs(U(r, i)) > real_t{16} * static_cast<real_t>(i + 1) * eps<T>() * mag(r, i));
    };

    for (size_t i = 0; i < N; ++i) {
        // Partial pivot among the entries that are not roundoff; a row-graded matrix can
        // put the largest entry of the column in a row where it is pure cancellation.
        size_t max_row = N;
        real_t max_val = real_t{0};
        for (size_t r = i; r < N; ++r) {
            const real_t val = damp::abs(U(r, i));
            if (!roundoff(r, i) && (max_row == N || val > max_val)) {
                max_val = val;
                max_row = r;
            }
        }

        // The pivot must clear the N·64·eps floor. One below it either is lost to cancellation
        // or still carries digits while being too small; skipping it for a smaller entry would
        // perturb A by more than the pivot. A finite A can still overflow to Inf − Inf = NaN in
        // the elimination, and under -ffinite-math-only the roundoff and floor tests can pass a
        // NaN, so finiteness is tested on its bits.
        if (max_row == N || !detail::is_finite_entry(U(max_row, i)) || !(max_val > rel * mag(max_row, i))) {
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

        // Elimination. A roundoff entry larger than the pivot (one classical partial pivoting
        // would have chosen) is taken as zero, a perturbation within its rounding bound that
        // keeps |lⱼᵢ| ≤ 1; every other entry is eliminated as usual. Column i below the pivot
        // is stored as the exact zero it is in U.
        for (size_t j = i + 1; j < N; ++j) {
            const bool zero = roundoff(j, i) && damp::abs(U(j, i)) > damp::abs(U(i, i));
            const T    factor = zero ? T{0} : (U(j, i) / U(i, i));
            U(j, i) = T{0};
            L(j, i) = factor;
            if (zero) {
                continue;
            }
            for (size_t k = i + 1; k < N; ++k) {
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
 * v is normalized to v[k] = 1 (LAPACK xLARFG), so |v[i]| ≤ 1 and β = 1 + |X(k, k)|/‖X(k…, k)‖
 * lies in [1, 2]. The norm is scaled (damp::detail::scaled_norm2) and no entry is squared, so
 * columns anywhere in the finite range reflect: a plain sum of squares overflows past √max
 * (1.8e19 in float) and flushes below √min.
 *
 * @see "Matrix Computations" (Golub & Van Loan, 2013), §5.1.2–5.1.3 and §5.1.13 (complex)
 * @see LAPACK xLARFG
 *
 * @return β, or 0 when the column is already zero from row k down (no reflection).
 */
template<typename T, size_t N, size_t M>
[[nodiscard]] constexpr scalar_type_t<T> householder_vector(const Matrix<N, M, T>& X, size_t k, damp::array<T, N>& v) {
    using real_t = scalar_type_t<T>;
    const real_t nrm = damp::detail::scaled_norm2<real_t>(N - k, [&](size_t i) { return real_t{damp::abs(X(k + i, k))}; });
    if (!(nrm > real_t{0})) {
        return real_t{0};
    }
    const real_t r0 = damp::abs(X(k, k));
    const T      phase = (r0 > real_t{0}) ? (X(k, k) / r0) : T{1};
    const T      vk = phase * (r0 + nrm); // X(k, k) − α with α = −phase·‖x‖: no cancellation
    v[k] = T{1};
    for (size_t i = k + 1; i < N; ++i) {
        v[i] = X(i, k) / vk;
    }
    return real_t{1} + (r0 / nrm);
}

/** @brief X ← (I − β·v·vᴴ)·X on rows k…N−1, columns c0…C−1. */
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
struct [[nodiscard]] QRDecomposition {
    Matrix<N, M, T> Q{}; ///< Orthonormal columns (thin Q)
    Matrix<M, M, T> R{}; ///< Upper-triangular factor
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
 * The roundoff test is scale-invariant per column: when |R(j, j)|, the part of column j
 * orthogonal to the columns before it, is at most rel_tol times ‖A(:, j)‖, it is roundoff and
 * R(j, j) is stored as zero, so a triangular solve on R refuses instead of dividing by it. Q keeps an
 * orthonormal column there, so QR reconstructs A to within that dropped roundoff. Columns past
 * N (M > N) have R(j, j) = 0 and a zero Q column.
 *
 * That zero is not a rank test. Without column pivoting each column is measured against the
 * ones before it only, so a near-dependency can hide in the column order: for
 * A = [a₀, a₀ + 10⁻⁴·a₂, a₂] with σ₃/σ₁ = 10⁻¹⁷ every |R(j, j)| clears its floor, while the same
 * columns in the order [a₀, a₂, a₀ + 10⁻⁴·a₂] zero the last one. Decide rank with mat::rank.
 *
 * @note Compare with MATLAB®'s [Q,R] = qr(A, 0) (economy size).
 * @see Golub & Van Loan, "Matrix Computations" (4th ed., 2013), §5.1.6 (backward accumulation) and §5.2.2
 * @see "Accuracy and Stability of Numerical Algorithms" (Higham, 2002), §19.3 (Householder QR) and §19.8 (MGS orthogonality loss)
 *
 * @tparam T   Element type
 * @tparam N   Rows
 * @tparam M   Columns
 * @param A        Input matrix
 * @param rel_tol  Relative column-norm floor below which R(j, j) is roundoff (default N·64·eps)
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
        const real_t col = damp::detail::scaled_norm2<real_t>(N, [&](size_t i) { return real_t{damp::abs(A(i, j))}; });
        if (!(damp::abs(result.R(j, j)) > rel_tol * col)) {
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
struct [[nodiscard]] FullQR {
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
        // Apply H from the left: R −= β·v·(vᴴR). Column k below the diagonal is zero by
        // construction; store it as such.
        detail::apply_householder_left(out.R, v, beta, k, 0);
        for (size_t i = k + 1; i < N; ++i) {
            out.R(i, k) = T{0};
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
