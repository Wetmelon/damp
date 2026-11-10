// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file eigen.hpp
 * @brief Eigenvalue algorithms
 */


#include "core.hpp"
#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"

namespace damp {
namespace mat {

/**
 * @brief Eigenvalue computation result.
 *
 * Returned by compute_eigenvalues(). Eigenvalues are fully resolved complex
 * values, a list in the order of the diagonal of the real Schur form; a complex
 * pair occupies two consecutive entries, positive imaginary part first.
 * @c schur_vectors is not a matrix of eigenvectors: it is the real Schur
 * basis D·Z of the balanced matrix mapped back to A, with unit columns. Its
 * leading k columns span an invariant subspace of A whenever k does not split a
 * 2×2 block (a complex pair); it is orthogonal only when A needed no balancing.
 * Column 0 is an eigenvector only when values[0] is real, and every column is an
 * eigenvector when A is symmetric (the Schur form is then diagonal).
 *
 * @tparam N Number of states
 * @tparam T Numeric type (default: double)
 */
template<size_t N, typename T = double>
struct [[nodiscard]] EigenResult {
    damp::array<damp::complex<T>, N> values{};        ///< Eigenvalues, in Schur-diagonal order
    Matrix<N, N, T>                  schur_vectors{}; ///< Real Schur basis D·Z with unit columns (nested invariant subspaces)
    bool                             success{false};  ///< Every eigenvalue deflated within the budget, from finite input

    /** @brief True when every eigenvalue deflated, from finite input. */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }
};

namespace detail {

/**
 * @brief Multiply @p x by 2^@p exp.
 *
 * Balancing scales are exact powers of two (LAPACK xGEBAL). Repeated doubling
 * keeps the factor exact inside the exponent range of @c T.
 */
template<typename T>
[[nodiscard]] constexpr T scale_by_pow2(T x, int exp) {
    if (exp > 0) {
        for (int k = 0; k < exp; ++k) {
            x *= T{2};
        }
    } else if (exp < 0) {
        for (int k = 0; k < -exp; ++k) {
            x *= static_cast<T>(0.5);
        }
    }
    return x;
}

/**
 * @brief One Parlett–Reinsch radix step on a row/column pair.
 *
 * @p c and @p r are the 2-norms of the column and the row, @p ca and @p ra the
 * matching max-norms. The returned exponent is log2(F) for the diagonal
 * similarity that multiplies the column by F and the row by 1/F. @c apply is
 * false when the step would not shrink c + r by the LAPACK factor 0.95.
 *
 * @see LAPACK xGEBAL (scale-only), Parlett & Reinsch (1969)
 */
template<typename T>
struct GebalStep {
    int  exp{0};
    bool apply{false};
};

template<typename T>
[[nodiscard]] constexpr GebalStep<T> gebal_step(T c, T r, T ca, T ra) {
    GebalStep<T> step;
    if (!(c > T{0}) || !(r > T{0}) || !damp::isfinite(c + ca + r + ra)) {
        return step;
    }
    const T s = c + r;
    const T sfmin1 = damp::numeric_limits<T>::min() / eps<T>();
    const T sfmin2 = sfmin1 * T{2};
    const T sfmax2 = T{1} / sfmin2;

    T   f = T{1};
    T   g = r * static_cast<T>(0.5);
    int f_exp = 0;
    int guard = 0;
    while (c < g && damp::max(damp::max(f, c), ca) < sfmax2 && damp::min(damp::min(r, g), ra) > sfmin2
           && guard < 4096) {
        f *= T{2};
        c *= T{2};
        ca *= T{2};
        r *= static_cast<T>(0.5);
        g *= static_cast<T>(0.5);
        ra *= static_cast<T>(0.5);
        ++f_exp;
        ++guard;
    }
    g = c * static_cast<T>(0.5);
    guard = 0;
    while (g >= r && damp::max(r, ra) < sfmax2 && damp::min(damp::min(f, c), damp::min(g, ca)) > sfmin2
           && guard < 4096) {
        f *= static_cast<T>(0.5);
        c *= static_cast<T>(0.5);
        g *= static_cast<T>(0.5);
        ca *= static_cast<T>(0.5);
        r *= T{2};
        ra *= T{2};
        --f_exp;
        ++guard;
    }
    //! (c + r) still within 5% of its original value: xGEBAL skips the update.
    if ((c + r) >= static_cast<T>(0.95) * s || f_exp == 0) {
        return step;
    }
    step.exp = f_exp;
    step.apply = true;
    return step;
}

/**
 * @brief log2 of the xGEBAL scale factors of @p A.
 *
 * @p A is overwritten by D⁻¹ A D with D = diag(2^exp). The row and column norms
 * include whatever diagonal the caller leaves: balance_exponents keeps it, as
 * LAPACK ≥ 3.5.0 does, and the Hamiltonian balance clears it, as Benner's
 * symplectic recipe and SciPy's CARE do. Row and column norms are scaled
 * (damp::detail::scaled_norm2), so entries anywhere in the finite range balance:
 * a plain sum of squares overflows past √max and flushes small rows to zero.
 *
 * Only the indices in [@p lo, @p hi] are scaled, and their norms count only that block, as after
 * permute_and_scale's isolation: the rest of A is triangular and its eigenvalues are already exact.
 * The max-norms that guard against overflow take the whole row and column, which the similarity
 * scales.
 *
 * @return false, with @p A and @p exp untouched, if @p A holds a non-finite entry.
 * @see LAPACK xGEBAL, JOB='S'
 * @see "Balancing a Matrix for Calculation of Eigenvalues and Eigenvectors" (Parlett & Reinsch, 1969)
 */
template<size_t M, typename T>
[[nodiscard]] constexpr bool gebal_log2_scale(Matrix<M, M, T>& A, damp::array<int, M>& exp, size_t lo = 0, size_t hi = M - 1) {
    if (!mat::all_of(A, [](T v) { return damp::isfinite(v); })) {
        return false;
    }
    if constexpr (M == 0) {
        return true; // hi = M − 1 wraps: nothing to scale
    }
    const T      sfmin1 = damp::numeric_limits<T>::min() / eps<T>();
    const T      sfmax1 = T{1} / sfmin1;
    const size_t width = hi - lo + 1;

    for (int sweep = 0; sweep < 64; ++sweep) {
        bool noconv = false;
        for (size_t i = lo; i <= hi; ++i) {
            T ca = T{0};
            T ra = T{0};
            for (size_t j = 0; j < M; ++j) {
                ra = damp::max(ra, damp::abs(A(i, j)));
                ca = damp::max(ca, damp::abs(A(j, i)));
            }
            const T    c = damp::detail::scaled_norm2<T>(width, [&](size_t k) { return damp::abs(A(lo + k, i)); });
            const T    r = damp::detail::scaled_norm2<T>(width, [&](size_t k) { return damp::abs(A(i, lo + k)); });
            const auto step = gebal_step(c, r, ca, ra);
            if (!step.apply) {
                continue;
            }
            const int f_exp = step.exp;
            if (f_exp < 0 && exp[i] < 0) {
                const T scale_now = scale_by_pow2(T{1}, exp[i]);
                const T f = scale_by_pow2(T{1}, f_exp);
                if (f * scale_now <= sfmin1) {
                    continue;
                }
            }
            if (f_exp > 0 && exp[i] > 0) {
                const T scale_now = scale_by_pow2(T{1}, exp[i]);
                const T f = scale_by_pow2(T{1}, f_exp);
                if (scale_now >= sfmax1 / f) {
                    continue;
                }
            }
            exp[i] += f_exp;
            noconv = true;
            for (size_t j = 0; j < M; ++j) {
                if (j == i) {
                    continue;
                }
                A(i, j) = scale_by_pow2(A(i, j), -f_exp);
                A(j, i) = scale_by_pow2(A(j, i), f_exp);
            }
        }
        if (!noconv) {
            break;
        }
    }
    return true;
}

/**
 * @brief xGEBAL scaling exponents of @p A: D = diag(2^e) with D⁻¹ A D balanced.
 *
 * Row and column 2-norms of the balanced matrix agree to within a factor of two. A
 * similarity transform, so the eigenvalues are unchanged, but the QR iteration then
 * computes them to a relative accuracy that does not depend on the units of each state.
 *
 * The norms include the diagonal, as xGEBAL's have since LAPACK 3.5.0. Without it a matrix
 * whose diagonal dominates is over-balanced. A result mapped back by D (eigenvectors, a
 * Schur basis, √A, log A) has a backward error bounded only by κ(D)·‖D⁻¹AD‖_F/‖A‖_F
 * units of roundoff. A similarity of a 3×3 Jordan block took κ(D) = 16 for a 1.16× smaller
 * norm, and its square root lost up to 9× in accuracy; with the diagonal, D = I.
 *
 * Scaling only, over the whole matrix: the conditioning step for rank and pencil tests that do
 * not reduce A itself. A Schur reduction balances with permute_and_scale (real_schur), which
 * first isolates the eigenvalues a permutation exposes; with nothing to isolate the two agree.
 *
 * @see "Balancing a Matrix for Calculation of Eigenvalues and Eigenvectors" (Parlett & Reinsch, 1969)
 * @see "On matrix balancing and eigenvector computation" (James, Langou & Lowery, 2014), §2–4
 * @see LAPACK xGEBAL (JOB = 'S'); MATLAB®'s balance
 */
template<size_t N, typename T>
[[nodiscard]] constexpr damp::array<int, N> balance_exponents(const Matrix<N, N, T>& A) {
    Matrix<N, N, T> mag{};
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            mag(i, j) = damp::abs(A(i, j));
        }
    }
    damp::array<int, N> e{};
    if (!gebal_log2_scale(mag, e)) {
        return damp::array<int, N>{}; // a non-finite entry: leave A unscaled
    }
    return e;
}

/** @brief D⁻¹ A D for D = diag(2^e); exact (powers of two). */
template<size_t N, typename T>
[[nodiscard]] constexpr Matrix<N, N, T> balanced(const Matrix<N, N, T>& A, const damp::array<int, N>& e) {
    Matrix<N, N, T> out{};
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            out(i, j) = scale_by_pow2(A(i, j), e[j] - e[i]);
        }
    }
    return out;
}

/**
 * @brief A permutation P and a diagonal D = diag(2^e) with B = D⁻¹·Pᵀ·A·P·D balanced.
 *
 * B is upper triangular outside the block [lo, hi]: its rows below hi and its columns left of
 * lo hold an eigenvalue of A on the diagonal, exactly. D scales only the block (e is 0 outside).
 */
template<size_t N>
struct Balancing {
    damp::array<size_t, N> perm{}; ///< Row and column k of B are row and column perm[k] of A
    damp::array<int, N>    e{};    ///< log2 of D, by the index of B
    size_t                 lo{0};  ///< First index of the block the QR iteration must reduce
    size_t                 hi{0};  ///< Last index of that block
};

/** @brief Exchange indices i and j of the similarity: rows and columns of @p M, and @p perm. */
template<size_t N, typename T>
constexpr void swap_indices(Matrix<N, N, T>& M, damp::array<size_t, N>& perm, size_t i, size_t j) {
    if (i == j) {
        return;
    }
    for (size_t k = 0; k < N; ++k) {
        const T t = M(i, k);
        M(i, k) = M(j, k);
        M(j, k) = t;
    }
    for (size_t k = 0; k < N; ++k) {
        const T t = M(k, i);
        M(k, i) = M(k, j);
        M(k, j) = t;
    }
    const size_t t = perm[i];
    perm[i] = perm[j];
    perm[j] = t;
}

/** @brief True when row (@p column false) or column @p k of @p M is zero on [lo, hi] off the diagonal. */
template<size_t N, typename T>
[[nodiscard]] constexpr bool isolated_index(const Matrix<N, N, T>& M, size_t k, size_t lo, size_t hi, bool column) {
    for (size_t j = lo; j <= hi; ++j) {
        if (j != k && (column ? M(j, k) : M(k, j)) != T{0}) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Isolate the eigenvalues a permutation exposes, then scale the rest (xGEBAL, JOB = 'B').
 *
 * A row whose entries left of column hi are zero off the diagonal holds an eigenvalue on its
 * diagonal: it is exchanged to position hi and the search repeats on the leading hi × hi part.
 * Then a column whose entries in rows [lo, hi] are zero off the diagonal moves to position lo, and
 * the search repeats on [lo + 1, hi]. No column found then can turn a row into an isolated one: an
 * isolated column has no nonzero in a row of the block. Each search moves one index or ends, so
 * each takes at most N passes. The block [lo, hi] is then scaled (gebal_log2_scale).
 *
 * The exchanges are exact, so a triangular A, or one a permutation makes triangular (a cascade,
 * an integrator appended to a plant), keeps every eigenvalue exact: the Hessenberg reduction and
 * the QR iteration never touch those rows. Scaling alone mixes them, and a graded triangular
 * matrix then loses its small eigenvalues (and the sign of one).
 *
 * @see "Balancing a Matrix for Calculation of Eigenvalues and Eigenvectors" (Parlett & Reinsch, 1969), §2 (permutation)
 * @see EISPACK `balanc`; LAPACK xGEBAL (JOB = 'B')
 */
template<size_t N, typename T>
[[nodiscard]] constexpr Balancing<N> permute_and_scale(const Matrix<N, N, T>& A) {
    Balancing<N> out{};
    if constexpr (N == 0) {
        return out;
    }
    Matrix<N, N, T> mag{};
    for (size_t i = 0; i < N; ++i) {
        out.perm[i] = i;
        for (size_t j = 0; j < N; ++j) {
            mag(i, j) = damp::abs(A(i, j));
        }
    }
    size_t hi = N - 1;
    for (size_t pass = 0; pass < N && hi > 0; ++pass) {
        size_t found = N;
        for (size_t i = hi + 1; i-- > 0;) {
            if (isolated_index(mag, i, 0, hi, false)) {
                found = i;
                break;
            }
        }
        if (found == N) {
            break;
        }
        swap_indices(mag, out.perm, found, hi);
        --hi;
    }
    size_t lo = 0;
    for (size_t pass = 0; pass < N && lo < hi; ++pass) {
        size_t found = N;
        for (size_t j = lo; j <= hi; ++j) {
            if (isolated_index(mag, j, lo, hi, true)) {
                found = j;
                break;
            }
        }
        if (found == N) {
            break;
        }
        swap_indices(mag, out.perm, found, lo);
        ++lo;
    }
    out.lo = lo;
    out.hi = hi;
    if (lo < hi && !gebal_log2_scale(mag, out.e, lo, hi)) {
        out.e = damp::array<int, N>{}; // a non-finite entry: leave the block unscaled
    }
    return out;
}

/**
 * @brief Reduce a square matrix to upper Hessenberg form by Householder reflections.
 *
 * On return @p H is upper Hessenberg and @p Z holds the orthogonal similarity
 * transform with Zᵀ·A·Z = H (eigenvalues preserved). @p Z is seeded with the
 * identity and accumulates the reflections, so the subsequent Francis sweeps can
 * continue accumulating Schur vectors into the same matrix.
 *
 * @see EISPACK `orthes` / `ortran`
 */
template<typename T, size_t N>
constexpr void hessenberg_reduce(Matrix<N, N, T>& H, Matrix<N, N, T>& Z) {
    Z = Matrix<N, N, T>::identity();
    if constexpr (N > 2) {
        damp::array<T, N> ort{};
        for (size_t m = 1; m + 1 < N; ++m) {
            T scale = T{0};
            for (size_t i = m; i < N; ++i) {
                scale += damp::abs(H(i, m - 1));
            }
            if (scale == T{0}) {
                continue;
            }
            // Build the Householder vector for column m-1 below the subdiagonal.
            T h = T{0};
            for (size_t i = N; i-- > m;) {
                ort[i] = H(i, m - 1) / scale;
                h += ort[i] * ort[i];
            }
            T g = damp::sqrt(h);
            if (ort[m] > T{0}) {
                g = -g;
            }
            h -= ort[m] * g;
            ort[m] -= g;
            // H ← (I − v·vᵀ/h)·H   (left)
            for (size_t j = m; j < N; ++j) {
                T f = T{0};
                for (size_t i = N; i-- > m;) {
                    f += ort[i] * H(i, j);
                }
                f /= h;
                for (size_t i = m; i < N; ++i) {
                    H(i, j) -= f * ort[i];
                }
            }
            // H ← H·(I − v·vᵀ/h)   (right)
            for (size_t i = 0; i < N; ++i) {
                T f = T{0};
                for (size_t j = N; j-- > m;) {
                    f += ort[j] * H(i, j);
                }
                f /= h;
                for (size_t j = m; j < N; ++j) {
                    H(i, j) -= f * ort[j];
                }
            }
            // Z ← Z·(I − v·vᵀ/h)   (accumulate Schur vectors)
            for (size_t i = 0; i < N; ++i) {
                T f = T{0};
                for (size_t j = N; j-- > m;) {
                    f += ort[j] * Z(i, j);
                }
                f /= h;
                for (size_t j = m; j < N; ++j) {
                    Z(i, j) -= f * ort[j];
                }
            }
            H(m, m - 1) = scale * g;
            for (size_t i = m + 1; i < N; ++i) {
                H(i, m - 1) = T{0};
            }
        }
    }
}

/** @brief Eigenvalues of a standardized 2×2 Schur block (wi = 0 for a real pair). */
template<typename T>
struct SchurBlockEigenvalues {
    T wr0{T{0}}; ///< Real part, leading position
    T wi0{T{0}}; ///< Imaginary part, leading position (positive for a complex pair)
    T wr1{T{0}}; ///< Real part, trailing position
    T wi1{T{0}}; ///< Imaginary part, trailing position
};

/**
 * @brief Apply the rotation G = [[c, −s], [s, c]] as Gᵀ·H·G on rows/columns i, i+1.
 *
 * H is upper quasi-triangular outside the 2×2 block at i, so the rows are rotated from
 * column i on and the columns down to row i+1. The Schur vectors take Z·G.
 */
template<typename T, size_t N>
constexpr void rotate_schur_pair(Matrix<N, N, T>& H, Matrix<N, N, T>& Z, size_t i, T c, T s) {
    for (size_t j = i; j < N; ++j) {
        const T h0 = H(i, j);
        const T h1 = H(i + 1, j);
        H(i, j) = (c * h0) + (s * h1);
        H(i + 1, j) = (c * h1) - (s * h0);
    }
    for (size_t r = 0; r <= i + 1; ++r) {
        const T h0 = H(r, i);
        const T h1 = H(r, i + 1);
        H(r, i) = (c * h0) + (s * h1);
        H(r, i + 1) = (c * h1) - (s * h0);
    }
    for (size_t r = 0; r < N; ++r) {
        const T z0 = Z(r, i);
        const T z1 = Z(r, i + 1);
        Z(r, i) = (c * z0) + (s * z1);
        Z(r, i + 1) = (c * z1) - (s * z0);
    }
}

/**
 * @brief √(|b|·|c|) without forming the product: max·√(min/max).
 *
 * The imaginary part of a standardized 2×2 block, @f$ \sqrt{|b|}\sqrt{|c|} @f$, is no safer
 * written that way: -ffast-math rewrites it into @f$ \sqrt{|bc|} @f$, whose product overflows past
 * √max and flushes to zero below √min (about 1.8·10¹⁹ and 1.1·10⁻¹⁹ in float). The quotient here
 * is at most 1 and underflows only when the result is below max·√min, negligible beside max.
 */
template<typename T>
[[nodiscard]] constexpr T sqrt_abs_product(T b, T c) {
    const T hi = damp::max(damp::abs(b), damp::abs(c));
    const T lo = damp::min(damp::abs(b), damp::abs(c));
    if (!(hi > T{0})) {
        return T{0};
    }
    return hi * damp::sqrt(lo / hi);
}

/**
 * @brief Bring the 2×2 diagonal block at @p i to the standard real Schur form.
 *
 * A block with real eigenvalues is triangularized; one with a complex pair is rotated to
 * equal diagonal entries with off-diagonal entries of opposite sign, so its eigenvalues are
 * @f$ m \pm i\sqrt{|b|}\sqrt{|c|} @f$. Both are orthogonal similarities on rows and columns
 * i, i+1, accumulated into @p Z.
 *
 * For @f$ \bigl[\begin{smallmatrix} a & b \\ c & d \end{smallmatrix}\bigr] @f$ the
 * eigenvalues are real when @f$ p^2 + bc \ge 0 @f$, @f$ p = (a-d)/2 @f$, evaluated on the
 * entries divided by @f$ \max(|p|, |b|, |c|) @f$ so it cannot overflow. A real pair comes
 * from the root @f$ z = p + \operatorname{sign}(p)\sqrt{p^2+bc} @f$, which never cancels:
 * @f$ \lambda_0 = d + z @f$, @f$ \lambda_1 = d - (b/z)c @f$, and the rotation whose first
 * column is the eigenvector @f$ [z, c] @f$ zeroes c. Taking @f$ (a+d)/2 \pm \sqrt{\cdot} @f$
 * instead would lose a small eigenvalue next to a large one. A complex pair takes the
 * rotation by θ with @f$ (\cos 2\theta, \sin 2\theta) \propto (b+c,\, d-a) @f$, chosen with
 * @f$ \cos 2\theta \ge 0 @f$ so @f$ \cos\theta \ge 1/\sqrt 2 @f$; it equalizes the diagonal
 * and leaves @f$ b - c @f$, the skew part, unchanged.
 *
 * A lower triangular block (b = 0) is exchanged by the exact permutation G = [[0, −1], [1, 0]]
 * to @f$ \bigl[\begin{smallmatrix} d & -c \\ 0 & a \end{smallmatrix}\bigr] @f$: the real-root
 * path would form a as @f$ d + (a - d) @f$, which cancels when |a| ≪ |d|. No entry is squared:
 * the rotations take damp::hypot, the imaginary part sqrt_abs_product, so the block's
 * eigenvalues do not overflow or flush to zero with its scale.
 *
 * @see "Matrix Computations" (Golub & Van Loan, 2013), §7.4.1 (real Schur form)
 * @see "On swapping diagonal blocks in real Schur form" (Bai & Demmel, 1993), §2 (standardized blocks)
 *
 * @return The block's eigenvalues; a complex pair has wi0 > 0 and wi1 = −wi0.
 */
template<typename T, size_t N>
constexpr SchurBlockEigenvalues<T> standardize_schur_block(Matrix<N, N, T>& H, Matrix<N, N, T>& Z, size_t i) {
    SchurBlockEigenvalues<T> ev{};
    const T                  a = H(i, i);
    const T                  b = H(i, i + 1);
    const T                  c = H(i + 1, i);
    const T                  d = H(i + 1, i + 1);
    if (c == T{0}) {
        ev.wr0 = a;
        ev.wr1 = d;
        return ev;
    }
    if (b == T{0}) {
        rotate_schur_pair(H, Z, i, T{0}, T{1});
        H(i, i) = d;
        H(i + 1, i + 1) = a;
        H(i + 1, i) = T{0};
        ev.wr0 = d;
        ev.wr1 = a;
        return ev;
    }
    const T p = static_cast<T>(0.5) * (a - d);
    const T scale = damp::max(damp::abs(p), damp::max(damp::abs(b), damp::abs(c)));
    const T ps = p / scale;
    const T disc = (ps * ps) + ((b / scale) * (c / scale)); // (p² + bc) / scale²
    if (disc >= T{0}) {
        const T root = scale * damp::sqrt(disc);
        const T z = p + damp::copysign(root, p);
        const T l0 = d + z;
        const T l1 = (z != T{0}) ? (d - ((b / z) * c)) : a; // z = 0 only when a = d and bc/scale² flushes
        const T r = damp::hypot(z, c);
        rotate_schur_pair(H, Z, i, z / r, c / r);
        H(i, i) = l0;
        H(i + 1, i + 1) = l1;
        H(i + 1, i) = T{0};
        ev.wr0 = l0;
        ev.wr1 = l1;
        return ev;
    }
    // Complex pair: equalize the diagonal.
    T       c2 = b + c;
    T       s2 = d - a;
    const T r2 = damp::hypot(c2, s2);
    if (r2 > T{0}) {
        if (c2 < T{0}) {
            c2 = -c2;
            s2 = -s2;
        }
        const T cs = damp::sqrt(static_cast<T>(0.5) * (T{1} + (c2 / r2)));
        const T sn = (s2 / r2) / (T{2} * cs);
        rotate_schur_pair(H, Z, i, cs, sn);
    }
    const T m = static_cast<T>(0.5) * (H(i, i) + H(i + 1, i + 1));
    H(i, i) = m;
    H(i + 1, i + 1) = m;
    const T b1 = H(i, i + 1);
    const T c1 = H(i + 1, i);
    const T g = sqrt_abs_product(b1, c1); // √|b₁|·√|c₁|
    if ((b1 > T{0} && c1 < T{0}) || (b1 < T{0} && c1 > T{0})) {
        ev.wr0 = m;
        ev.wr1 = m;
        ev.wi0 = g;
        ev.wi1 = -g;
        return ev;
    }
    // Rounding left a (nearly) double real eigenvalue: triangularize about m with the
    // eigenvector [√|b₁|, sign(c₁)√|c₁|] of m + √|b₁c₁|, of norm √(|b₁| + |c₁|).
    const T sb = damp::sqrt(damp::abs(b1));
    const T sc = damp::sqrt(damp::abs(c1));
    const T r = damp::sqrt(damp::abs(b1) + damp::abs(c1));
    if (r > T{0}) {
        rotate_schur_pair(H, Z, i, sb / r, damp::copysign(sc, c1) / r);
    }
    H(i, i) = m + g;
    H(i + 1, i + 1) = m - g;
    H(i + 1, i) = T{0};
    ev.wr0 = H(i, i);
    ev.wr1 = H(i + 1, i + 1);
    return ev;
}

/**
 * @brief Francis double-shift QR on an upper Hessenberg matrix.
 *
 * Reduces @p H to real Schur (quasi-triangular) form in place, accumulating the
 * orthogonal transforms into @p Z, and writes the eigenvalues to @p wr (real
 * parts) and @p wi (imaginary parts). Complex eigenvalues emerge as conjugate
 * pairs from the surviving 2×2 diagonal blocks. The implicit double shift makes
 * this robust for complex and equal-magnitude spectra that unshifted QR cannot
 * deflate.
 *
 * A subdiagonal entry is taken as zero by the Ahues–Tisseur criterion: the usual
 * @f$ |h_{k,k-1}| \le \epsilon(|h_{k-1,k-1}| + |h_{kk}|) @f$, and then also
 * @f$ |h_{k,k-1}|\,|h_{k-1,k}| \le \epsilon\,|h_{kk}|\,|h_{k-1,k-1} - h_{kk}| @f$, the
 * product that actually moves the eigenvalue. The first test alone deflates too early on a
 * graded matrix and costs a small eigenvalue all its digits. Exceptional shifts come every
 * ten sweeps on the current eigenvalue (the count restarts at each deflation); a separate
 * count bounds the total work. A converged 2×2 block is standardized
 * (standardize_schur_block).
 *
 * Outside the block [@p lo, @p hi] H must be upper triangular (permute_and_scale isolates such
 * rows): their eigenvalues are read off the diagonal, and the iteration and its shifts touch only
 * the block, so they stay exact. The shift products are of entries in the block, which
 * schur_reduce brings into a range where a product of two cannot overflow; under -ffast-math a
 * local scale factor would not help, since the compiler folds (x·k)(y·k) back into x·y·k².
 *
 * @return true if every eigenvalue deflated within the iteration budget.
 * @see EISPACK `hqr2`; Golub & Van Loan, "Matrix Computations" §7.5
 * @see "A New Deflation Criterion for the QR Algorithm" (Ahues & Tisseur, 1997), LAPACK Working Note 122
 *
 * Indexing: EISPACK/Golub–Van Loan use a signed active-block cursor that
 * counts down through −1 as the exit condition. Accesses cast to size_t; the
 * conversion is intentional and bounded (nn starts at hi and only decreases).
 */
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
template<typename T, size_t N>
constexpr bool francis_qr(
    Matrix<N, N, T>&   H,
    Matrix<N, N, T>&   Z,
    damp::array<T, N>& wr,
    damp::array<T, N>& wi,
    size_t             lo = 0,
    size_t             hi = N - 1
) {
    constexpr T  eps = damp::eps<T>();
    constexpr T  smlnum = damp::numeric_limits<T>::min() * (static_cast<T>(N) / eps);
    constexpr T  quotient_guard = T{8} / damp::numeric_limits<T>::max(); // |x|·guard < |y| ⇒ |x/y| < max/8
    const size_t max_its = (30 * N) + 30;
    const int    first = static_cast<int>(lo);

    T anorm = T{0};
    for (size_t i = 0; i < N; ++i) {
        if (i < lo || i > hi) {
            wr[i] = H(i, i); // isolated: the eigenvalue is on the diagonal
            wi[i] = T{0};
            continue;
        }
        for (size_t j = ((i > lo) ? i - 1 : lo); j <= hi; ++j) {
            anorm += damp::abs(H(i, j));
        }
    }

    int    nn = static_cast<int>(hi);
    T      t = T{0}; // accumulated shift
    bool   converged = true;
    size_t its = 0;      // total sweeps: the budget
    size_t its_here = 0; // sweeps on the current eigenvalue: the exceptional-shift schedule

    // Ahues–Tisseur: may H(l, l−1) be taken as zero? Diagonal entries carry −t.
    const auto negligible = [&](int l, int last) {
        const T h = damp::abs(H(l, l - 1));
        if (h <= smlnum) {
            return true;
        }
        const T hll = H(l, l) + t;
        const T hpp = H(l - 1, l - 1) + t;
        T       tst = damp::abs(hpp) + damp::abs(hll);
        if constexpr (N >= 3) {
            if (tst == T{0}) {
                if (l >= first + 2) {
                    tst += damp::abs(H(l - 1, l - 2));
                }
                if (l + 1 <= last) {
                    tst += damp::abs(H(l + 1, l));
                }
            }
        }
        if (tst == T{0}) {
            tst = anorm;
        }
        if (h > eps * tst) {
            return false;
        }
        const T up = damp::abs(H(l - 1, l));
        const T ab = damp::max(h, up);
        const T ba = damp::min(h, up);
        const T gap = damp::abs(hpp - hll);
        const T aa = damp::max(damp::abs(hll), gap);
        const T bb = damp::min(damp::abs(hll), gap);
        const T sc = aa + ab;
        return ba * (ab / sc) <= damp::max(smlnum, eps * (bb * (aa / sc)));
    };

    while (nn >= first) {
        bool deflated = false;
        while (!deflated) {
            // Locate the first small subdiagonal element (deflation point).
            int l = nn;
            while (l > first) {
                if (negligible(l, nn)) {
                    H(l, l - 1) = T{0};
                    break;
                }
                --l;
            }

            T x = H(nn, nn);
            if (l == nn) {
                // One real eigenvalue has converged. Restore the accumulated
                // shift to the diagonal so H stays a valid Schur form (the
                // recorded eigenvalue already carries t; the matrix must too).
                H(nn, nn) = x + t;
                wr[nn] = x + t;
                wi[nn] = T{0};
                --nn;
                its_here = 0;
                deflated = true;
                break;
            }

            T y = H(nn - 1, nn - 1);
            T w = H(nn, nn - 1) * H(nn - 1, nn);
            if (l == nn - 1) {
                // A 2×2 block has converged. Restore the accumulated shift to its diagonal
                // (off-diagonals are shift-invariant), then bring it to standard form.
                H(nn - 1, nn - 1) = y + t;
                H(nn, nn) = x + t;
                const auto ev = standardize_schur_block(H, Z, static_cast<size_t>(nn - 1));
                wr[nn - 1] = ev.wr0;
                wi[nn - 1] = ev.wi0;
                wr[nn] = ev.wr1;
                wi[nn] = ev.wi1;
                nn -= 2;
                its_here = 0;
                deflated = true;
                break;
            }

            // 1×1 and 2×2 active blocks always deflate via the two branches
            // above; only a 3×3-or-larger active window ever reaches the bulge
            // chase below. Guarding it with `if constexpr (N >= 3)` keeps the
            // H(nn-1, nn-2)-style indices out of the N<3 instantiations entirely
            // — they are unreachable there, but the optimizer cannot prove it
            // and would emit a -Warray-bounds false positive for N == 2.
            if constexpr (N >= 3) {
                if (its >= max_its) {
                    // Give up on this eigenvalue but keep going (honest convergence flag).
                    converged = false;
                    H(nn, nn) = x + t;
                    wr[nn] = x + t;
                    wi[nn] = T{0};
                    --nn;
                    deflated = true;
                    break;
                }

                // Exceptional shift every 10 sweeps on this eigenvalue to break out of cycles.
                if (its_here > 0 && its_here % 10 == 0) {
                    t += x;
                    for (int i = first; i <= nn; ++i) {
                        H(i, i) -= x;
                    }
                    T s = damp::abs(H(nn, nn - 1)) + damp::abs(H(nn - 1, nn - 2));
                    x = static_cast<T>(0.75) * s;
                    y = x;
                    w = static_cast<T>(-0.4375) * s * s;
                }
                ++its;
                ++its_here;

                // Find two consecutive small subdiagonals to start the bulge.
                int m = nn - 2;
                T   p = T{0};
                T   q = T{0};
                T   r = T{0};
                while (m >= l) {
                    T z2 = H(m, m);
                    r = x - z2;
                    T s = y - z2;
                    //! First column of (H − σ₁)(H − σ₂) from row m, divided by h₍ₘ₊₁,ₘ₎. Beside a
                    //! tiny h₍ₘ₊₁,ₘ₎ the quotient could overflow: then the column is kept undivided,
                    //! the same direction.
                    const T h21 = H(m + 1, m);
                    const T num = (r * s) - w;
                    q = H(m + 1, m + 1) - z2 - r - s;
                    r = H(m + 2, m + 1);
                    if (damp::abs(num) * quotient_guard < damp::abs(h21)) {
                        p = (num / h21) + H(m, m + 1);
                    } else {
                        p = num + (H(m, m + 1) * h21);
                        q *= h21;
                        r *= h21;
                    }
                    T s2 = damp::abs(p) + damp::abs(q) + damp::abs(r);
                    p /= s2;
                    q /= s2;
                    r /= s2;
                    if (m == l) {
                        break;
                    }
                    T u = damp::abs(H(m, m - 1)) * (damp::abs(q) + damp::abs(r));
                    T v = damp::abs(p) * (damp::abs(H(m - 1, m - 1)) + damp::abs(z2) + damp::abs(H(m + 1, m + 1)));
                    if (u <= eps * v) {
                        break;
                    }
                    --m;
                }
                for (int i = m + 2; i <= nn; ++i) {
                    H(i, i - 2) = T{0};
                    if (i != m + 2) {
                        H(i, i - 3) = T{0};
                    }
                }

                // Double-shift QR sweep: chase the bulge down the band, applying
                // the 3×3 (and trailing 2×2) Householder reflections to H rows, H
                // columns, and the Schur-vector accumulator Z.
                for (int k = m; k <= nn - 1; ++k) {
                    if (k != m) {
                        p = H(k, k - 1);
                        q = H(k + 1, k - 1);
                        r = (k != nn - 1) ? H(k + 2, k - 1) : T{0};
                        x = damp::abs(p) + damp::abs(q) + damp::abs(r);
                        if (x != T{0}) {
                            p /= x;
                            q /= x;
                            r /= x;
                        }
                    }
                    T s = damp::copysign(damp::sqrt((p * p) + (q * q) + (r * r)), p);
                    if (s == T{0}) {
                        continue;
                    }
                    if (k == m) {
                        if (l != m) {
                            H(k, k - 1) = -H(k, k - 1);
                        }
                    } else {
                        H(k, k - 1) = -s * x;
                    }
                    //! The reflector I − τ·v·vᵀ, v = [1, q/p, r/p], τ = p/s, with p ← p + s: the row
                    //! factor (τ, τ·v₁, τ·v₂) is formed from v, not as [p, q, r]/s, so the
                    //! transform stays symmetric when -ffast-math turns x/s into x·(1/s). Divided
                    //! separately, square roots of 3×3 Jordan blocks lost about a quarter of their
                    //! accuracy in float and double.
                    p += s;
                    const T px = p / s;
                    const T qq = q / p;
                    const T rr = r / p;
                    const T qx = qq * px;
                    const T rx = rr * px;
                    // Row modification.
                    for (int j = k; j < static_cast<int>(N); ++j) {
                        p = H(k, j) + (qq * H(k + 1, j));
                        if (k != nn - 1) {
                            p += rr * H(k + 2, j);
                            H(k + 2, j) -= p * rx;
                        }
                        H(k + 1, j) -= p * qx;
                        H(k, j) -= p * px;
                    }
                    const int row_max = (nn < k + 3) ? nn : (k + 3);
                    // Column modification.
                    for (int i = 0; i <= row_max; ++i) {
                        p = (px * H(i, k)) + (qx * H(i, k + 1));
                        if (k != nn - 1) {
                            p += rx * H(i, k + 2);
                            H(i, k + 2) -= p * rr;
                        }
                        H(i, k + 1) -= p * qq;
                        H(i, k) -= p;
                    }
                    // Schur-vector accumulation.
                    for (int i = 0; i < static_cast<int>(N); ++i) {
                        p = (px * Z(i, k)) + (qx * Z(i, k + 1));
                        if (k != nn - 1) {
                            p += rx * Z(i, k + 2);
                            Z(i, k + 2) -= p * rr;
                        }
                        Z(i, k + 1) -= p * qq;
                        Z(i, k) -= p;
                    }
                }
            }
        }
    }
    // The bulge chase leaves entries below the subdiagonal that the reflections annihilated
    // in exact arithmetic but never stored as zero (EISPACK clears them only at the start of
    // the next sweep). The Schur form has structural zeros there.
    for (size_t i = 2; i < N; ++i) {
        for (size_t j = 0; j + 1 < i; ++j) {
            H(i, j) = T{0};
        }
    }
    return converged;
}
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

/**
 * @brief Even exponent σ with 2^−σ·max |hᵢⱼ| over the block [lo, hi] in [√min/ε, ε/√min].
 *
 * The range is LAPACK xGEEV's: a product of two entries there neither overflows nor flushes to
 * zero (about [9·10⁻¹³, 1.1·10¹²] in float, [7·10⁻¹³⁹, 1.5·10¹³⁸] in double). An even σ keeps
 * every square root exact under the scaling, so a matrix already inside (σ = 0) and one brought
 * there round alike.
 */
template<typename T, size_t N>
[[nodiscard]] constexpr int schur_prescale_exponent(const Matrix<N, N, T>& H, size_t lo, size_t hi) {
    constexpr T small = damp::sqrt(damp::numeric_limits<T>::min()) / damp::eps<T>();
    constexpr T big = T{1} / small;
    T           amax = T{0};
    for (size_t i = lo; i <= hi; ++i) {
        for (size_t j = lo; j <= hi; ++j) {
            amax = damp::max(amax, damp::abs(H(i, j)));
        }
    }
    int sigma = 0;
    for (int k = 0; k < 1024 && amax > big; ++k) {
        amax *= static_cast<T>(0.25);
        sigma += 2;
    }
    for (int k = 0; k < 1024 && amax > T{0} && amax < small; ++k) {
        amax *= T{4};
        sigma -= 2;
    }
    return sigma;
}

/**
 * @brief Real Schur form of @p H in place: Hessenberg reduction and Francis QR, at a safe scale.
 *
 * The block [lo, hi] is multiplied by 2^−σ (schur_prescale_exponent) before the reduction and its
 * Schur form, eigenvalues included, by 2^σ after, both exact. Only the block: outside it H must be
 * upper triangular, and the coupling rows and columns take orthogonal transforms computed from the
 * block alone, which do not depend on its scale. On return @p H is the standardized real Schur
 * form, Zᵀ·H₀·Z = H with @p Z orthogonal, and @p wr, @p wi hold the eigenvalues in its diagonal
 * order. Without the prescale a block whose entries pass √max (1.8·10¹⁹ in float) returns an
 * infinite imaginary part, and one below √min a complex pair as a real double eigenvalue.
 *
 * @return true if every eigenvalue deflated within the iteration budget.
 * @see LAPACK xGEEV (scaling A into [SMLNUM, BIGNUM] before xGEBAL and xHSEQR)
 */
template<typename T, size_t N>
[[nodiscard]] constexpr bool schur_reduce(Matrix<N, N, T>& H, Matrix<N, N, T>& Z, damp::array<T, N>& wr, damp::array<T, N>& wi, size_t lo = 0, size_t hi = N - 1) {
    Z = Matrix<N, N, T>::identity();
    if constexpr (N == 0) {
        return true;
    } else if constexpr (N == 1) {
        wr[0] = H(0, 0);
        wi[0] = T{0};
        return true;
    } else {
        const int sigma = schur_prescale_exponent(H, lo, hi);
        for (size_t i = lo; sigma != 0 && i <= hi; ++i) {
            for (size_t j = lo; j <= hi; ++j) {
                H(i, j) = scale_by_pow2(H(i, j), -sigma);
            }
        }
        hessenberg_reduce(H, Z);
        const bool ok = francis_qr(H, Z, wr, wi, lo, hi);
        for (size_t i = lo; sigma != 0 && i <= hi; ++i) {
            for (size_t j = lo; j <= hi; ++j) {
                H(i, j) = scale_by_pow2(H(i, j), sigma);
            }
            wr[i] = scale_by_pow2(wr[i], sigma);
            wi[i] = scale_by_pow2(wi[i], sigma);
        }
        return ok;
    }
}

/**
 * @brief Real Schur form of a balanced A: A = D·Z·T·Zᵀ·D⁻¹, D = diag(2^e), Z orthogonal.
 *
 * Every row of Z, and every entry of e, is indexed as A's; the permutation of the balance is
 * folded into Z, so Z stays orthogonal and D diagonal. T is in A's units.
 */
template<size_t N, typename T>
struct RealSchur {
    Matrix<N, N, T>     T_form{};         ///< Standardized real Schur form, quasi upper triangular
    Matrix<N, N, T>     Z{};              ///< Orthogonal Schur vectors of D⁻¹AD
    damp::array<int, N> e{};              ///< log2 of the balancing D
    damp::array<T, N>   wr{};             ///< Real parts of the eigenvalues, in T's diagonal order
    damp::array<T, N>   wi{};             ///< Imaginary parts, positive first within a pair
    bool                converged{false}; ///< Finite A, every eigenvalue deflated, finite T
};

/**
 * @brief Balance (permute_and_scale), reduce (schur_reduce) and map back, for every real Schur caller.
 *
 * The eigenvalues the permutation isolates are exact, the rest are accurate relative to the
 * balanced block, and neither depends on the overall scale of A from the smallest normal number
 * to the largest. A non-finite entry, a QR iteration that does not deflate, or a Schur form that
 * leaves the finite range on the way back (a true eigenvalue past the largest number) returns
 * converged = false.
 *
 * @see "Balancing a Matrix for Calculation of Eigenvalues and Eigenvectors" (Parlett & Reinsch, 1969)
 * @see LAPACK xGEEV (balance, scale, reduce, then undo the scaling and the balance)
 */
template<size_t N, typename T>
[[nodiscard]] constexpr RealSchur<N, T> real_schur(const Matrix<N, N, T>& A) {
    RealSchur<N, T> out{};
    if (!mat::all_of(A, [](T v) { return damp::isfinite(v); })) {
        return out;
    }
    const Balancing<N> bal = permute_and_scale(A);
    Matrix<N, N, T>    H{};
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            H(i, j) = scale_by_pow2(A(bal.perm[i], bal.perm[j]), bal.e[j] - bal.e[i]);
        }
    }
    Matrix<N, N, T> Zb{};
    const bool      ok = schur_reduce(H, Zb, out.wr, out.wi, bal.lo, bal.hi);
    out.T_form = H;
    for (size_t i = 0; i < N; ++i) {
        out.e[bal.perm[i]] = bal.e[i];
        for (size_t j = 0; j < N; ++j) {
            out.Z(bal.perm[i], j) = Zb(i, j);
        }
    }
    out.converged = ok && mat::all_of(H, [](T v) { return damp::isfinite(v); });
    for (size_t i = 0; i < N; ++i) {
        out.converged = out.converged && damp::isfinite(out.wr[i]) && damp::isfinite(out.wi[i]);
    }
    return out;
}

} // namespace detail

/**
 * @brief Compute the eigenvalues (and Schur vectors) of a real square matrix.
 *
 * Resolves both real and complex-conjugate eigenvalues for any N, returned as
 * fully-resolved complex values. 1×1 is trivial; larger matrices go through
 * detail::real_schur: balanced (detail::permute_and_scale, a permutation and a
 * diagonal similarity by powers of two), brought to a safe scale, then reduced by
 * Hessenberg + Francis double-shift QR (a 2×2 deflates at once via the
 * cancellation-free quadratic). An eigenvalue the permutation isolates (a
 * triangular A, a cascade, an appended integrator) is the diagonal entry itself,
 * exactly. Balancing makes the others accurate relative to the balanced norm, not
 * to the largest entry of A, so a state in metres next to one in counts does not
 * cost the small eigenvalues their digits, and the prescale makes them independent
 * of A's overall scale across the whole floating-point range. The @c schur_vectors
 * field holds the real Schur basis of the balanced matrix mapped back (x = D x_b)
 * with unit columns: nested invariant subspaces of A, orthogonal when A needed no
 * scaling, and the eigenvectors for symmetric A. It does not hold eigenvectors of a
 * general A.
 *
 * @note Compare with MATLAB®'s eig(A) — the algorithm is selected internally.
 * @warning QR is iterative (data-dependent iteration count). Suitable for
 *          design-time analysis; avoid in hard-real-time loops where WCET matters.
 *
 * The iteration budget scales with N internally and a subdiagonal entry deflates
 * by the local Ahues–Tisseur criterion (see detail::francis_qr), so no tolerance
 * or iteration cap is taken. A non-finite entry in A, or an eigenvalue past the
 * largest finite number, returns success = false.
 *
 * @param A input matrix
 * @return EigenResult with complex eigenvalues, Schur vectors, and an honest
 *         convergence flag.
 */
template<size_t N, typename T = double>
    requires damp::is_floating_point_v<T> // the real Schur form: real A only
[[nodiscard]] constexpr EigenResult<N, T> compute_eigenvalues(const Matrix<N, N, T>& A) {
    using Cplx = damp::complex<T>;
    EigenResult<N, T> result;

    if (!mat::all_of(A, [](T v) { return damp::isfinite(v); })) {
        result.success = false;
        return result;
    }

    if constexpr (N == 1) {
        // Scalar: the lone entry is the eigenvalue.
        result.values[0] = Cplx{A(0, 0), T{0}};
        result.schur_vectors(0, 0) = T{1};
        result.success = true;
    } else {
        // N == 2 deflates immediately through the stable 2×2 branch of
        // francis_qr (its bulge chase is compiled only for N ≥ 3).
        const auto rs = detail::real_schur(A);
        result.success = rs.converged;

        for (size_t i = 0; i < N; ++i) {
            result.values[i] = Cplx{rs.wr[i], rs.wi[i]};
            // Schur vectors of A are D·Z; renormalize each column (D may span many decades).
            for (size_t j = 0; j < N; ++j) {
                result.schur_vectors(j, i) = detail::scale_by_pow2(rs.Z(j, i), rs.e[j]);
            }
            const T nrm = damp::detail::scaled_norm2<T>(N, [&](size_t j) { return damp::abs(result.schur_vectors(j, i)); });
            const T inv = (nrm > T{0}) ? (T{1} / nrm) : T{1};
            for (size_t j = 0; j < N; ++j) {
                result.schur_vectors(j, i) *= inv;
            }
        }
    }
    return result;
}
} // namespace mat
} // namespace damp