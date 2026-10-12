// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file eigen.hpp
 * @brief Eigenvalue algorithms
 */

#include <cstddef>
#include <limits>

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
 * values; @c vectors holds the orthogonal real Schur basis as columns. For
 * symmetric A the Schur form is diagonal, so these are the eigenvectors; for a
 * general A only the first column is an eigenvector (complex pairs span a
 * 2-D invariant subspace).
 *
 * @tparam N Number of states
 * @tparam T Numeric type (default: double)
 */
template<size_t N, typename T = double>
struct EigenResult {
    ColVec<N, damp::complex<T>>    values{};        ///< Eigenvalues
    Matrix<N, N, damp::complex<T>> vectors{};       ///< Eigenvectors (columns)
    bool                           converged{true}; ///< Every eigenvalue deflated within the budget
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
    const T sfmin1 = std::numeric_limits<T>::min() / eps<T>();
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
 * @p A is overwritten by D⁻¹ A D with D = diag(2^exp). Diagonal entries are
 * ignored by the caller (they are zero), matching the Benner / SciPy recipe
 * of clearing diag(|H|) before balancing. Returns false if a norm overflows.
 *
 * @see LAPACK xGEBAL, JOB='S'
 */
template<size_t M, typename T>
[[nodiscard]] constexpr bool gebal_log2_scale(Matrix<M, M, T>& A, damp::array<int, M>& exp) {
    const T sfmin1 = std::numeric_limits<T>::min() / eps<T>();
    const T sfmax1 = T{1} / sfmin1;

    for (int sweep = 0; sweep < 64; ++sweep) {
        bool noconv = false;
        for (size_t i = 0; i < M; ++i) {
            T c2 = T{0};
            T r2 = T{0};
            T ca = T{0};
            T ra = T{0};
            for (size_t j = 0; j < M; ++j) {
                const T aij = damp::abs(A(i, j));
                const T aji = damp::abs(A(j, i));
                r2 += aij * aij;
                c2 += aji * aji;
                if (aij > ra) {
                    ra = aij;
                }
                if (aji > ca) {
                    ca = aji;
                }
            }
            const T    c = damp::sqrt(c2);
            const T    r = damp::sqrt(r2);
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
 * Row and column norms of the balanced matrix agree to within a factor of two
 * (the diagonal is ignored, as in LAPACK). A similarity transform, so the
 * eigenvalues are unchanged, but the QR iteration then computes them to a
 * relative accuracy that does not depend on the units of each state.
 *
 * @see "Balancing a Matrix for Calculation of Eigenvalues and Eigenvectors" (Parlett & Reinsch, 1969)
 * @see LAPACK xGEBAL (JOB = 'S'); MATLAB®'s balance
 */
template<size_t N, typename T>
[[nodiscard]] constexpr damp::array<int, N> balance_exponents(const Matrix<N, N, T>& A) {
    Matrix<N, N, T> mag = Matrix<N, N, T>::zeros();
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            if (i != j) {
                mag(i, j) = damp::abs(A(i, j));
            }
        }
    }
    damp::array<int, N> e{};
    if (!gebal_log2_scale(mag, e)) {
        return damp::array<int, N>{}; // a norm overflowed: leave A unscaled
    }
    return e;
}

/// D⁻¹ A D for D = diag(2^e); exact (powers of two).
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

/**
 * @brief Rotate a deflated 2×2 block with real eigenvalues to upper-triangular form.
 *
 * Applies the Givens rotation of EISPACK `hqr2` to rows/columns @p n−1, @p n of
 * @p H and to the Schur vectors @p Z, so @c H(n,n−1) → 0 and the leading
 * eigenvalue moves to @c H(n−1,n−1). @p z is the stable root offset
 * @f$ p + \operatorname{sign}(p)\sqrt{p^2 + w} @f$ of the block.
 */
template<typename T, size_t N>
constexpr void triangularize_real_block(Matrix<N, N, T>& H, Matrix<N, N, T>& Z, size_t n, T z) {
    const T sub = H(n, n - 1);
    const T s = damp::abs(sub) + damp::abs(z);
    if (s == T{0}) {
        return;
    }
    T       p = sub / s;
    T       q = z / s;
    const T r = damp::sqrt((p * p) + (q * q));
    p /= r;
    q /= r;
    for (size_t j = n - 1; j < N; ++j) {
        const T h = H(n - 1, j);
        H(n - 1, j) = (q * h) + (p * H(n, j));
        H(n, j) = (q * H(n, j)) - (p * h);
    }
    for (size_t i = 0; i <= n; ++i) {
        const T h = H(i, n - 1);
        H(i, n - 1) = (q * h) + (p * H(i, n));
        H(i, n) = (q * H(i, n)) - (p * h);
    }
    for (size_t i = 0; i < N; ++i) {
        const T v = Z(i, n - 1);
        Z(i, n - 1) = (q * v) + (p * Z(i, n));
        Z(i, n) = (q * Z(i, n)) - (p * v);
    }
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
 * @return true if every eigenvalue deflated within the iteration budget.
 * @see EISPACK `hqr2`; Golub & Van Loan, "Matrix Computations" §7.5
 *
 * Indexing: EISPACK/Golub–Van Loan use a signed active-block cursor that
 * counts down through −1 as the exit condition. Accesses cast to size_t; the
 * conversion is intentional and bounded (nn starts at N−1 and only decreases).
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
    damp::array<T, N>& wi
) {
    constexpr T  eps = damp::eps<T>();
    const size_t max_its = (30 * N) + 30;

    T anorm = T{0};
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = (i ? i - 1 : 0); j < N; ++j) {
            anorm += damp::abs(H(i, j));
        }
    }

    int    nn = static_cast<int>(N) - 1;
    T      t = T{0}; // accumulated shift
    bool   converged = true;
    size_t its = 0;

    while (nn >= 0) {
        bool deflated = false;
        while (!deflated) {
            // Locate the first small subdiagonal element (deflation point).
            int l = nn;
            while (l > 0) {
                T s = damp::abs(H(l - 1, l - 1)) + damp::abs(H(l, l));
                if (s == T{0}) {
                    s = anorm;
                }
                if (damp::abs(H(l, l - 1)) <= eps * s) {
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
                deflated = true;
                break;
            }

            T y = H(nn - 1, nn - 1);
            T w = H(nn, nn - 1) * H(nn - 1, nn);
            if (l == nn - 1) {
                // A 2×2 block has converged: solve its characteristic quadratic.
                T p = static_cast<T>(0.5) * (y - x);
                T q = (p * p) + w;
                T z = damp::sqrt(damp::abs(q));
                x += t;
                // Restore the accumulated shift to the block diagonal so H stays
                // a valid real Schur form (off-diagonals are shift-invariant).
                H(nn - 1, nn - 1) = y + t;
                H(nn, nn) = x;
                if (q >= T{0}) {
                    // Real eigenvalue pair.
                    z = p + damp::copysign(z, p);
                    wr[nn - 1] = x + z;
                    wr[nn] = (z != T{0}) ? (x - (w / z)) : (x + z);
                    wi[nn - 1] = T{0};
                    wi[nn] = T{0};
                    triangularize_real_block(H, Z, static_cast<size_t>(nn), z);
                } else {
                    // Complex conjugate pair.
                    wr[nn - 1] = x + p;
                    wr[nn] = x + p;
                    wi[nn - 1] = z;
                    wi[nn] = -z;
                }
                nn -= 2;
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

                // Exceptional shift every 10 iterations to break out of cycles.
                if ((its % 10 == 9) && (its != 0)) {
                    t += x;
                    for (int i = 0; i <= nn; ++i) {
                        H(i, i) -= x;
                    }
                    T s = damp::abs(H(nn, nn - 1)) + damp::abs(H(nn - 1, nn - 2));
                    x = static_cast<T>(0.75) * s;
                    y = x;
                    w = static_cast<T>(-0.4375) * s * s;
                }
                ++its;

                // Find two consecutive small subdiagonals to start the bulge.
                int m = nn - 2;
                T   p = T{0};
                T   q = T{0};
                T   r = T{0};
                while (m >= l) {
                    T z2 = H(m, m);
                    r = x - z2;
                    T s = y - z2;
                    p = (((r * s) - w) / H(m + 1, m)) + H(m, m + 1);
                    q = H(m + 1, m + 1) - z2 - r - s;
                    r = H(m + 2, m + 1);
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
                    p += s;
                    const T px = p / s;
                    const T qx = q / s;
                    const T rx = r / s;
                    const T qq = q / p;
                    const T rr = r / p;
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
    return converged;
}
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

} // namespace detail

/**
 * @brief Compute the eigenvalues (and Schur vectors) of a real square matrix.
 *
 * Resolves both real and complex-conjugate eigenvalues for any N, returned as
 * fully-resolved complex values. 1×1 is trivial; larger matrices are balanced
 * (detail::balance_exponents, a diagonal similarity by powers of two) and then
 * reduced by Hessenberg + Francis double-shift QR (a 2×2 deflates at once via
 * the cancellation-free quadratic). Balancing makes the computed eigenvalues
 * accurate relative to the balanced norm, not to the largest entry of A, so a
 * state in metres next to one in counts does not cost the small eigenvalues
 * their digits. The @c vectors field holds the real Schur basis of the balanced
 * matrix mapped back (x = D x_b) with unit columns: nested invariant subspaces of
 * A, orthogonal when A needed no balancing, and the eigenvectors for symmetric A.
 *
 * @note Compare with MATLAB®'s eig(A) — the algorithm is selected internally.
 * @warning QR is iterative (data-dependent iteration count). Suitable for
 *          design-time analysis; avoid in hard-real-time loops where WCET matters.
 *
 * The iteration budget scales with N internally and deflation uses the machine
 * epsilon relative to the matrix norm, so no tolerance or iteration cap is taken.
 *
 * @param A input matrix
 * @return EigenResult with complex eigenvalues, Schur vectors, and an honest
 *         convergence flag.
 */
template<size_t N, typename T = double>
[[nodiscard]] constexpr EigenResult<N, T> compute_eigenvalues(const Matrix<N, N, T>& A) {
    using Cplx = damp::complex<T>;
    EigenResult<N, T> result;

    if constexpr (N == 1) {
        // Scalar: the lone entry is the eigenvalue.
        result.values[0] = Cplx{A(0, 0), T{0}};
        result.vectors(0, 0) = Cplx{T{1}, T{0}};
        result.converged = true;
    } else {
        // N == 2 deflates immediately through the stable 2×2 branch of
        // francis_qr (its bulge chase is compiled only for N ≥ 3).
        const auto      e = detail::balance_exponents(A);
        Matrix<N, N, T> H = detail::balanced(A, e);
        Matrix<N, N, T> Z;
        detail::hessenberg_reduce(H, Z);

        damp::array<T, N> wr{};
        damp::array<T, N> wi{};
        result.converged = detail::francis_qr(H, Z, wr, wi);

        for (size_t i = 0; i < N; ++i) {
            result.values[i] = Cplx{wr[i], wi[i]};
            // Schur vectors of A are D·Z; renormalize each column.
            T nrm_sq = T{0};
            for (size_t j = 0; j < N; ++j) {
                const T v = detail::scale_by_pow2(Z(j, i), e[j]);
                result.vectors(j, i) = Cplx{v, T{0}};
                nrm_sq += v * v;
            }
            const T inv = (nrm_sq > T{0}) ? (T{1} / damp::sqrt(nrm_sq)) : T{1};
            for (size_t j = 0; j < N; ++j) {
                result.vectors(j, i) = result.vectors(j, i) * inv;
            }
        }
    }
    return result;
}
} // namespace mat
} // namespace damp