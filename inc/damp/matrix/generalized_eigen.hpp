// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file generalized_eigen.hpp
 * @brief Generalized eigenvalues of a real square pencil A − λB by the QZ algorithm
 *
 * The pencil is reduced to Hessenberg–triangular form (A upper Hessenberg, B upper triangular)
 * by Householder reflections and plane rotations, then the implicit double-shift QZ iteration
 * drives A to quasi-triangular form while B stays triangular. Every transformation is
 * orthogonal, so the computed eigenvalues are exact for a pencil within a few ε of (A, B) in
 * norm. The arithmetic is real, as in mat::compute_eigenvalues: a complex conjugate pair comes
 * out of a 2×2 diagonal block as an exact conjugate pair.
 *
 * A zero (or negligible) diagonal entry of B is an infinite eigenvalue. It is chased to the
 * bottom of the active block by rotations and deflated there, so the finite eigenvalues are
 * computed from a pencil whose B part is safely nonsingular. A diagonal entry of B counts as
 * zero when it is at most n·ε·‖B‖_F: a perturbation of that size, which the backward error of
 * the reduction already allows, makes it exactly zero. The converse does not hold: rounding in
 * the reduction and the sweeps can leave the diagonal entry of an exactly singular B some ten
 * to a hundred times ε·‖B‖_F. That infinite eigenvalue then comes out finite, with |λ| of order
 * ‖A‖/(ε‖B‖) or more, as with LAPACK's xGGEV.
 *
 * @code
 * constexpr Matrix<3, 3> A{{1.0, 2.0, 0.0}, {0.0, 3.0, 1.0}, {1.0, 0.0, 2.0}};
 * constexpr Matrix<3, 3> B{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}};
 * constexpr auto         ev = mat::generalized_eigenvalues(A, B);
 * static_assert(ev.success && ev.n_finite == 2); // B has rank 2: one eigenvalue at infinity
 * @endcode
 *
 * @see "An Algorithm for Generalized Matrix Eigenvalue Problems" (Moler & Stewart, 1973), SIAM J.
 *      Numer. Anal. 10(2)
 * @see "Matrix Computations" (Golub & Van Loan, 2013), §7.7 (the QZ method)
 */


#include "core.hpp"
#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "eigen.hpp" // scale_by_pow2

namespace damp {
namespace mat {

/**
 * @brief Generalized eigenvalues of A − λB: the finite ones, and how many there are
 *
 * The finite eigenvalues fill @c values[0..n_finite) in the order of their diagonal position in
 * the generalized Schur form, top first; the other N − n_finite eigenvalues are infinite (B
 * singular). A complex pair occupies two consecutive entries, positive imaginary part first.
 *
 * @tparam N Pencil dimension
 * @tparam T Floating-point scalar type
 */
template<size_t N, typename T = double>
struct [[nodiscard]] GeneralizedEigenResult {
    damp::array<damp::complex<T>, N> values{};       ///< Finite eigenvalues λ with det(A − λB) = 0, the first n_finite entries
    size_t                           n_finite{0};    ///< Number of finite eigenvalues; the other N − n_finite are infinite
    bool                             success{false}; ///< The pencil is finite and every eigenvalue deflated within the iteration budget

    /** @brief True when the eigenvalues were found (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The same eigenvalues in the scalar @p U. */
    template<typename U>
    [[nodiscard]] constexpr GeneralizedEigenResult<N, U> as() const {
        GeneralizedEigenResult<N, U> out{};
        for (size_t i = 0; i < N; ++i) {
            out.values[i] = damp::complex<U>{static_cast<U>(values[i].real()), static_cast<U>(values[i].imag())};
        }
        out.n_finite = n_finite;
        out.success = success;
        return out;
    }
};

namespace detail {

/** @brief A plane rotation [c s; −s c]. */
template<typename T>
struct PlaneRotation {
    T c{1}; ///< Cosine
    T s{0}; ///< Sine
};

/** @brief The rotation G with G·[f; g] = [r; 0], r = ‖(f, g)‖ (the identity when both are zero). */
template<typename T>
[[nodiscard]] constexpr PlaneRotation<T> rotation_onto_first(T f, T g) {
    const T scale = damp::max(damp::abs(f), damp::abs(g));
    if (!(scale > T{0})) {
        return {};
    }
    const T fs = f / scale;
    const T gs = g / scale;
    const T r = damp::sqrt((fs * fs) + (gs * gs));
    return {fs / r, gs / r};
}

/** @brief Rows i, i+1 of M on the columns [c0, c1) ← G·(rows i, i+1). */
template<size_t N, typename T>
constexpr void rotate_rows(Matrix<N, N, T>& M, size_t i, size_t c0, size_t c1, const PlaneRotation<T>& g) {
    for (size_t j = c0; j < c1; ++j) {
        const T x = M(i, j);
        const T y = M(i + 1, j);
        M(i, j) = (g.c * x) + (g.s * y);
        M(i + 1, j) = (g.c * y) - (g.s * x);
    }
}

/**
 * @brief Columns j, j+1 of M on the rows [r0, r1) ← (columns j, j+1)·Gᵀ
 *
 * With g = rotation_onto_first(q, p) this maps a row's pair [p, q] to [0, r]: it zeroes the
 * entry in column j.
 */
template<size_t N, typename T>
constexpr void rotate_cols(Matrix<N, N, T>& M, size_t j, size_t r0, size_t r1, const PlaneRotation<T>& g) {
    for (size_t i = r0; i < r1; ++i) {
        const T x = M(i, j);
        const T y = M(i, j + 1);
        M(i, j) = (g.c * x) - (g.s * y);
        M(i, j + 1) = (g.s * x) + (g.c * y);
    }
}

/** @brief A Householder reflector H = I − β·u·uᵀ on three consecutive indices. */
template<typename T>
struct Reflector3 {
    damp::array<T, 3> u{};     ///< Householder vector, scaled to its largest entry
    T                 beta{0}; ///< 2/(uᵀu), or 0 for the identity
};

/** @brief The reflector with H·(x₀, x₁, x₂) = ∓‖x‖·e_t, the sign chosen so u_t does not cancel. */
template<typename T>
[[nodiscard]] constexpr Reflector3<T> reflector3(T x0, T x1, T x2, size_t t) {
    Reflector3<T> h{};
    const T       scale = damp::max(damp::abs(x0), damp::max(damp::abs(x1), damp::abs(x2)));
    if (!(scale > T{0})) {
        return h;
    }
    h.u = {x0 / scale, x1 / scale, x2 / scale};
    const T nrm = damp::sqrt((h.u[0] * h.u[0]) + (h.u[1] * h.u[1]) + (h.u[2] * h.u[2]));
    const T ut = h.u[t];
    h.u[t] = ut + damp::copysign(nrm, ut);
    h.beta = T{1} / (nrm * (nrm + damp::abs(ut)));
    return h;
}

/** @brief Rows k … k+2 of M on the columns [c0, c1) ← H·(rows k … k+2). */
template<size_t N, typename T>
constexpr void reflect_rows3(Matrix<N, N, T>& M, const Reflector3<T>& h, size_t k, size_t c0, size_t c1) {
    for (size_t j = c0; j < c1; ++j) {
        const T w = h.beta * ((h.u[0] * M(k, j)) + (h.u[1] * M(k + 1, j)) + (h.u[2] * M(k + 2, j)));
        M(k, j) -= w * h.u[0];
        M(k + 1, j) -= w * h.u[1];
        M(k + 2, j) -= w * h.u[2];
    }
}

/** @brief Columns k … k+2 of M on the rows [r0, r1) ← (columns k … k+2)·H. */
template<size_t N, typename T>
constexpr void reflect_cols3(Matrix<N, N, T>& M, const Reflector3<T>& h, size_t k, size_t r0, size_t r1) {
    for (size_t i = r0; i < r1; ++i) {
        const T w = h.beta * ((M(i, k) * h.u[0]) + (M(i, k + 1) * h.u[1]) + (M(i, k + 2) * h.u[2]));
        M(i, k) -= w * h.u[0];
        M(i, k + 1) -= w * h.u[1];
        M(i, k + 2) -= w * h.u[2];
    }
}

/**
 * @brief Reduce the leading n×n block of the pencil (A, B) to Hessenberg–triangular form
 *
 * B = QR by Householder reflections (Qᵀ applied to A), then for each column of A the entries
 * below the subdiagonal are rotated out from the bottom up. Each row rotation fills in one
 * entry under B's diagonal, which a column rotation removes again.
 *
 * @see "Matrix Computations" (Golub & Van Loan, 2013), Algorithm 7.7.1
 */
template<size_t N, typename T>
constexpr void hessenberg_triangular(Matrix<N, N, T>& A, Matrix<N, N, T>& B, size_t n) {
    damp::array<T, N> u{};
    for (size_t j = 0; j + 1 < n; ++j) {
        T scale = T{0};
        for (size_t i = j; i < n; ++i) {
            scale = damp::max(scale, damp::abs(B(i, j)));
        }
        if (!(scale > T{0})) {
            continue;
        }
        T nrm_sq = T{0};
        for (size_t i = j; i < n; ++i) {
            u[i] = B(i, j) / scale;
            nrm_sq += u[i] * u[i];
        }
        const T nrm = damp::sqrt(nrm_sq);
        const T ut = u[j];
        u[j] = ut + damp::copysign(nrm, ut);
        const T    beta = T{1} / (nrm * (nrm + damp::abs(ut)));
        const auto reflect = [&](Matrix<N, N, T>& M, size_t c0) {
            for (size_t c = c0; c < n; ++c) {
                T w = T{0};
                for (size_t i = j; i < n; ++i) {
                    w += u[i] * M(i, c);
                }
                w *= beta;
                for (size_t i = j; i < n; ++i) {
                    M(i, c) -= w * u[i];
                }
            }
        };
        reflect(B, j);
        reflect(A, 0);
        for (size_t i = j + 1; i < n; ++i) {
            B(i, j) = T{0};
        }
    }
    //! A 2×2 pencil is already Hessenberg; compiling the rotations only for N ≥ 3 keeps their
    //! indices out of instantiations where they cannot run (a -Warray-bounds false positive).
    if constexpr (N >= 3) {
        for (size_t j = 0; j + 2 < n; ++j) {
            for (size_t i = n - 1; i >= j + 2; --i) {
                const auto row = rotation_onto_first(A(i - 1, j), A(i, j));
                rotate_rows(A, i - 1, j, n, row);
                rotate_rows(B, i - 1, i - 1, n, row);
                A(i, j) = T{0};
                const auto col = rotation_onto_first(B(i, i), B(i, i - 1));
                rotate_cols(B, i - 1, 0, i + 1, col);
                rotate_cols(A, i - 1, 0, n, col);
                B(i, i - 1) = T{0};
            }
        }
    }
}

/** @brief The two eigenvalues of a 2×2 block of a Hessenberg–triangular pencil. */
template<typename T>
struct PencilPair {
    damp::complex<T> first{};  ///< Positive imaginary part for a complex pair
    damp::complex<T> second{}; ///< Conjugate of first for a complex pair
};

/**
 * @brief Roots of det(A₂ − λB₂) for A₂ = [a₀₀ a₀₁; a₁₀ a₁₁], B₂ = [b₀₀ b₀₁; 0 b₁₁], b₀₀b₁₁ ≠ 0
 *
 * A₂ and B₂ are divided by their largest entries, and the roots are found relative to a shift
 * σ: the diagonal ratio s₀ = a₀₀/b₀₀ or s₁ = a₁₁/b₁₁ nearer zero, so that A₂ − σB₂ rounds least.
 * For σ = s₀ the (0, 0) entry of A₂ − σB₂ is zero up to a rounding of a₀₀, which is dropped as
 * a backward error, and μ = λ − σ solves
 * @f[ \mu^2 - 2p\mu - q = 0,\quad p = \tfrac12\bigl(s_1 - s_0 - \hat s\, b_{01}\bigr),\quad
 *     q = \hat s\,(a_{01} - \sigma b_{01}),\quad \hat s = \frac{a_{10}}{b_{00} b_{11}} @f]
 * (for σ = s₁, exchange s₀ and s₁ in p). p is half the gap between the two ratios, less
 * ŝb₀₁/2, so the discriminant @f$ p^2 + q @f$ is a square plus a product: no difference of
 * large squares. Two close eigenvalues keep their digits, and whether they are real is
 * decided by comparing |p| with √(−q), not by the sign of a quantity that cancelled. The
 * unshifted quadratic in λ loses ε·λ²/gap instead, enough to report a complex pair with a
 * small imaginary part as real.
 *
 * Neither p² nor q is formed, so a small block does not underflow in float:
 * @f$ \sqrt{|q|} = \sqrt{|\hat s|}\,\sqrt{|a_{01} - \sigma b_{01}|} @f$, and
 * @f$ \sqrt{p^2 + q} @f$ is the scaled norm of (p, √q) for q ≥ 0 and
 * @f$ \sqrt{|p| - \sqrt{|q|}}\,\sqrt{|p| + \sqrt{|q|}} @f$ for q < 0. A real pair takes
 * @f$ \mu = p + \operatorname{sign}(p)\sqrt{p^2 + q} @f$, which does not cancel, and the other
 * root as @f$ -q/\mu @f$. Adding σ back can still cancel in the eigenvalue of smaller magnitude;
 * when it is less than half the other, it is taken as det A₂/(b₀₀b₁₁) over the other instead,
 * so a small b₀₀ or b₁₁ gives one large eigenvalue without costing the other its digits.
 *
 * @see LAPACK xLAG2, whose shifted formulation it credits to C. Van Loan
 */
template<typename T>
[[nodiscard]] constexpr PencilPair<T> pencil_2x2_eigenvalues(T a00, T a01, T a10, T a11, T b00, T b01, T b11) {
    const T a_max = damp::max(damp::max(damp::abs(a00), damp::abs(a01)), damp::max(damp::abs(a10), damp::abs(a11)));
    const T b_max = damp::max(damp::max(damp::abs(b00), damp::abs(b01)), damp::abs(b11));
    if (!(a_max > T{0})) {
        return {};
    }
    //! Divisors clamped to [min, 1/min] keep their reciprocals normal under -freciprocal-math.
    constexpr T lo = damp::numeric_limits<T>::min();
    const T     sa = damp::clamp(a_max, lo, T{1} / lo);
    const T     sb = damp::clamp(b_max, lo, T{1} / lo);
    a00 /= sa;
    a01 /= sa;
    a10 /= sa;
    a11 /= sa;
    b00 /= sb;
    b01 /= sb;
    b11 /= sb;
    const T    s0 = a00 / b00;
    const T    s1 = a11 / b11;
    const T    s_hat = (a10 / b00) / b11;
    const bool shift_by_s0 = damp::abs(s0) <= damp::abs(s1);
    const T    sigma = shift_by_s0 ? s0 : s1;
    const T    other = shift_by_s0 ? s1 : s0;
    const T    p = static_cast<T>(0.5) * ((other - sigma) - (s_hat * b01));
    const T    a01_shifted = a01 - (sigma * b01);
    const T    root_q = damp::sqrt(damp::abs(s_hat)) * damp::sqrt(damp::abs(a01_shifted)); // √|q|
    const bool q_negative = (root_q > T{0}) && ((s_hat < T{0}) != (a01_shifted < T{0}));
    const T    abs_p = damp::abs(p);
    const T    unit = sa / sb;
    if (q_negative && abs_p < root_q) {
        const T re = sigma + p;
        const T im = damp::sqrt(root_q - abs_p) * damp::sqrt(root_q + abs_p);
        return {damp::complex<T>{re * unit, im * unit}, damp::complex<T>{re * unit, -im * unit}};
    }
    const T r = q_negative ? (damp::sqrt(abs_p - root_q) * damp::sqrt(abs_p + root_q))
                           : damp::detail::scaled_norm2<T>(2, [&](size_t k) { return (k == 0) ? abs_p : root_q; });
    const T mu_far = p + damp::copysign(r, p);
    T       mu_near = T{0}; // −q/μ_far; μ_far is zero only when p and q both are
    if (mu_far != T{0}) {
        const T q_over_mu = root_q * (root_q / mu_far);
        mu_near = q_negative ? q_over_mu : -q_over_mu;
    }
    T big = sigma + mu_far;
    T small = sigma + mu_near;
    if (damp::abs(small) > damp::abs(big)) {
        damp::swap(big, small);
    }
    if (damp::abs(small) < static_cast<T>(0.5) * damp::abs(big)) {
        small = ((((a00 * a11) - (a01 * a10)) / b00) / b11) / big;
    }
    return {damp::complex<T>{big * unit, T{0}}, damp::complex<T>{small * unit, T{0}}};
}

/**
 * @brief Move the zero at B(j, j) to B(h, h) and deflate it as an infinite eigenvalue at h
 *
 * A row rotation on rows k, k+1 zeroes B(k+1, k+1); the entry it fills in at A(k+1, k−1) is
 * removed by a column rotation on columns k−1, k, which also restores B(k−1, k−1), the zero the
 * step before started from. At the bottom a column rotation zeroes A(h, h−1) and restores
 * B(h−1, h−1), leaving the eigenvalue (A(h, h), 0) split off. The rotations act on the window
 * [l, h] only.
 *
 * @see "Matrix Computations" (Golub & Van Loan, 2013), §7.7.5 (deflation of a zero bᵢᵢ)
 */
template<size_t N, typename T>
constexpr void chase_infinite_eigenvalue(Matrix<N, N, T>& A, Matrix<N, N, T>& B, size_t l, size_t h, size_t j) {
    B(j, j) = T{0};
    for (size_t k = j; k < h; ++k) {
        const auto row = rotation_onto_first(B(k, k + 1), B(k + 1, k + 1));
        rotate_rows(B, k, k + 1, h + 1, row);
        rotate_rows(A, k, (k > l) ? k - 1 : l, h + 1, row);
        B(k + 1, k + 1) = T{0};
        if (k > l) {
            const auto col = rotation_onto_first(A(k + 1, k), A(k + 1, k - 1));
            rotate_cols(A, k - 1, l, k + 2, col);
            rotate_cols(B, k - 1, l, k + 1, col);
            A(k + 1, k - 1) = T{0};
        }
    }
    if (h > l) {
        const auto col = rotation_onto_first(A(h, h), A(h, h - 1));
        rotate_cols(A, h - 1, l, h + 1, col);
        rotate_cols(B, h - 1, l, h + 1, col);
        A(h, h - 1) = T{0};
    }
}

/**
 * @brief One implicit double-shift QZ sweep on the unreduced window [l, h] (h ≥ l + 2)
 *
 * The shifts are the eigenvalues of the trailing 2×2 block of M = A·B⁻¹, and the sweep starts
 * from the first column of (M − σ₁I)(M − σ₂I), which has three nonzero entries. M is never
 * formed: the entries needed follow from M·B = A by back substitution on B's leading 2×2 and
 * trailing 3×3 blocks. An exceptional sweep takes the pair c ± i·0.66e instead, c = m_hh + 0.75e,
 * e = |m_{h,h−1}| + |m_{h−1,h−2}|. The bulge is chased down by a reflector on three rows of the
 * pencil, then a reflector and a rotation on columns that restore B's triangle.
 *
 * @see "Matrix Computations" (Golub & Van Loan, 2013), Algorithm 7.7.2
 */
template<size_t N, typename T>
constexpr void qz_sweep(Matrix<N, N, T>& A, Matrix<N, N, T>& B, size_t l, size_t h, bool exceptional) {
    // Leading entries of M = A·B⁻¹.
    const T m00 = A(l, l) / B(l, l);
    const T m10 = A(l + 1, l) / B(l, l);
    const T m01 = (A(l, l + 1) - (m00 * B(l, l + 1))) / B(l + 1, l + 1);
    const T m11 = (A(l + 1, l + 1) - (m10 * B(l, l + 1))) / B(l + 1, l + 1);
    const T m21 = A(l + 2, l + 1) / B(l + 1, l + 1);
    // Trailing 2×2 of M (rows and columns h−1, h) and the entry left of it.
    const T mpq = A(h - 1, h - 2) / B(h - 2, h - 2);
    const T mpp = (A(h - 1, h - 1) - (mpq * B(h - 2, h - 1))) / B(h - 1, h - 1);
    const T mpr = (A(h - 1, h) - (mpq * B(h - 2, h)) - (mpp * B(h - 1, h))) / B(h, h);
    const T mrp = A(h, h - 1) / B(h - 1, h - 1);
    const T mrr = (A(h, h) - (mrp * B(h - 1, h))) / B(h, h);
    T       x{0};
    T       y{0};
    T       z = m21 * m10;
    if (exceptional) {
        const T e = damp::abs(mrp) + damp::abs(mpq);
        const T c = mrr + (static_cast<T>(0.75) * e);
        const T s = T{2} * c;
        const T t = (c * c) + (static_cast<T>(0.4375) * e * e);
        x = (m00 * (m00 - s)) + t + (m01 * m10);
        y = m10 * (m00 + m11 - s);
    } else {
        x = ((m00 - mpp) * (m00 - mrr)) - (mpr * mrp) + (m01 * m10);
        y = m10 * ((m00 - mpp) + (m11 - mrr));
    }
    for (size_t k = l; k + 2 <= h; ++k) {
        if (k > l) {
            x = A(k, k - 1);
            y = A(k + 1, k - 1);
            z = A(k + 2, k - 1);
        }
        const auto q = reflector3(x, y, z, 0);
        reflect_rows3(A, q, k, (k > l) ? k - 1 : l, h + 1);
        reflect_rows3(B, q, k, k, h + 1);
        if (k > l) {
            A(k + 1, k - 1) = T{0};
            A(k + 2, k - 1) = T{0};
        }
        const size_t a_rows = (k + 4 < h + 1) ? k + 4 : h + 1;
        const auto   z1 = reflector3(B(k + 2, k), B(k + 2, k + 1), B(k + 2, k + 2), 2);
        reflect_cols3(A, z1, k, l, a_rows);
        reflect_cols3(B, z1, k, l, k + 3);
        B(k + 2, k) = T{0};
        B(k + 2, k + 1) = T{0};
        const auto z2 = rotation_onto_first(B(k + 1, k + 1), B(k + 1, k));
        rotate_cols(A, k, l, a_rows, z2);
        rotate_cols(B, k, l, k + 2, z2);
        B(k + 1, k) = T{0};
    }
    const auto q = rotation_onto_first(A(h - 1, h - 2), A(h, h - 2));
    rotate_rows(A, h - 1, h - 2, h + 1, q);
    rotate_rows(B, h - 1, h - 1, h + 1, q);
    A(h, h - 2) = T{0};
    const auto z2 = rotation_onto_first(B(h, h), B(h, h - 1));
    rotate_cols(A, h - 1, l, h + 1, z2);
    rotate_cols(B, h - 1, l, h + 1, z2);
    B(h, h - 1) = T{0};
}

/**
 * @brief Divide the leading n×n block of @p M by 2^e, e the binary exponent of its largest entry
 *
 * The largest entry then lies in [1, 2). Powers of two scale exactly, short of the subnormal range.
 *
 * @return e (0 for a zero block)
 */
template<size_t N, typename T>
constexpr int scale_to_unit_max(Matrix<N, N, T>& M, size_t n) {
    T m_max = T{0};
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            m_max = damp::max(m_max, damp::abs(M(i, j)));
        }
    }
    const int e = damp::detail::binary_exponent(m_max);
    if (e != 0) {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                M(i, j) = scale_by_pow2(M(i, j), -e);
            }
        }
    }
    return e;
}

/**
 * @brief QZ eigenvalues of the leading n×n block (2 ≤ n ≤ N) of the pencil (A, B), overwritten
 *
 * Hessenberg–triangular reduction, then QZ sweeps on the active window. A subdiagonal entry of
 * A deflates when it is at most ε(|a_{k−1,k−1}| + |a_kk|) (ε‖A‖_F when both are zero); a
 * diagonal entry of B at most n·ε·‖B‖_F is an infinite eigenvalue (chase_infinite_eigenvalue).
 * A converged 1×1 block gives a_hh/b_hh, a 2×2 block pencil_2x2_eigenvalues. Every tenth sweep
 * on the same eigenvalue is exceptional, and 30n + 30 sweeps in total bound the work. Only the
 * window is transformed, since the eigenvalues are all that is kept.
 *
 * @param lambda Eigenvalue at each position, valid where @p finite is set
 * @param finite Whether the eigenvalue at each position is finite
 * @return true if every eigenvalue deflated within the budget
 */
template<size_t N, typename T>
constexpr bool qz_eigenvalues(
    Matrix<N, N, T>&                  A,
    Matrix<N, N, T>&                  B,
    size_t                            n,
    damp::array<damp::complex<T>, N>& lambda,
    damp::array<bool, N>&             finite
) {
    hessenberg_triangular(A, B, n);
    const T      a_norm = A.norm();
    const T      b_tol = static_cast<T>(n) * eps<T>() * B.norm();
    const size_t max_its = (30 * n) + 30;
    size_t       its = 0;
    size_t       its_here = 0;
    size_t       active = n; // eigenvalues at positions [active, n) are found
    while (active > 0) {
        const size_t h = active - 1;
        size_t       l = h;
        while (l > 0) {
            T local = damp::abs(A(l - 1, l - 1)) + damp::abs(A(l, l));
            if (local == T{0}) {
                local = a_norm;
            }
            if (damp::abs(A(l, l - 1)) <= eps<T>() * local) {
                A(l, l - 1) = T{0};
                break;
            }
            --l;
        }
        size_t zero_at = n;
        for (size_t j = l; j <= h && zero_at == n; ++j) {
            if (damp::abs(B(j, j)) <= b_tol) {
                zero_at = j;
            }
        }
        if (zero_at < n) {
            chase_infinite_eigenvalue(A, B, l, h, zero_at);
            finite[h] = false;
            active = h;
            its_here = 0;
        } else if (l == h) {
            lambda[h] = damp::complex<T>{A(h, h) / B(h, h), T{0}};
            finite[h] = true;
            active = h;
            its_here = 0;
        } else if (l + 1 == h) {
            const auto pair = pencil_2x2_eigenvalues(A(l, l), A(l, h), A(h, l), A(h, h), B(l, l), B(l, h), B(h, h));
            lambda[l] = pair.first;
            lambda[h] = pair.second;
            finite[l] = true;
            finite[h] = true;
            active = l;
            its_here = 0;
        } else if constexpr (N >= 3) { // a window of three or more: only N ≥ 3 has one
            if (its >= max_its) {
                return false;
            }
            ++its;
            ++its_here;
            qz_sweep(A, B, l, h, its_here % 10 == 0);
        }
    }
    return true;
}

/**
 * @brief Eigenvalues of the leading n×n block of the pencil (A, B), overwritten by the iteration
 *
 * A and B are first divided by the powers of two that bring their largest entries into [1, 2)
 * (scale_to_unit_max), and every eigenvalue is multiplied back by 2^(e_A − e_B). The scaling is
 * exact, so it changes no digit of a pencil that was safe already. It is needed because the
 * first column of the double-shift polynomial is quadratic in A·B⁻¹ (qz_sweep): on the pencil
 * as given it overflowed or underflowed once ‖A‖/‖B‖ left about [2⁻⁶⁰, 2⁶⁰] in float
 * ([2⁻⁴⁹⁷, 2⁵⁰⁸] in double), and products of small entries flushed near the bottom of the
 * range. The eigenvalues are then found by qz_eigenvalues (n ≥ 2) or as a₀₀/b₀₀ (n = 1).
 *
 * @param lambda Eigenvalue at each position, valid where @p finite is set
 * @param finite Whether the eigenvalue at each position is finite
 * @return true if every eigenvalue deflated within the budget and every finite one is within the
 *         range of T (an eigenvalue of magnitude above its largest value is lost, not clamped)
 */
template<size_t N, typename T>
constexpr bool pencil_eigenvalues(
    Matrix<N, N, T>&                  A,
    Matrix<N, N, T>&                  B,
    size_t                            n,
    damp::array<damp::complex<T>, N>& lambda,
    damp::array<bool, N>&             finite
) {
    n = (n < N) ? n : N;
    const int e_a = scale_to_unit_max(A, n);
    const int e_b = scale_to_unit_max(B, n);
    bool      converged = true;
    if constexpr (N == 1) { // a 1×1 pencil: no iteration, and no 2×2 block indices to compile
        if (n == 1 && B(0, 0) != T{0}) {
            lambda[0] = damp::complex<T>{A(0, 0) / B(0, 0), T{0}};
            finite[0] = true;
        }
    } else {
        converged = qz_eigenvalues(A, B, n, lambda, finite);
    }
    bool representable = true;
    for (size_t i = 0; i < n; ++i) {
        if (finite[i]) {
            const T re = scale_by_pow2(lambda[i].real(), e_a - e_b);
            const T im = scale_by_pow2(lambda[i].imag(), e_a - e_b);
            lambda[i] = damp::complex<T>{re, im};
            representable = representable && damp::isfinite(re) && damp::isfinite(im);
        }
    }
    return converged && representable;
}

} // namespace detail

/**
 * @brief Generalized eigenvalues of the real pencil A − λB: the λ with det(A − λB) = 0
 *
 * Hessenberg–triangular reduction and the double-shift QZ iteration (see the file
 * description). The finite eigenvalues are returned with their count; the rest are infinite:
 * generically one for each rank B loses, though an infinite eigenvalue can come out as a
 * finite one of magnitude near ‖A‖/(ε‖B‖) (see the file description). With B = I they are the
 * eigenvalues of A. The pencil must be regular (det(A − λB) not identically zero): for a
 * singular pencil every λ makes A − λB singular, and the values returned are not determined
 * by it. A and B are each scaled by a power of two (exact), so any finite pencil whose
 * eigenvalues fit in T can be solved, whatever ‖A‖/‖B‖. No balancing is applied, so the
 * accuracy is relative to ‖A‖ and ‖B‖.
 *
 * @note Compare with MATLAB®'s eig(A, B) and its 'qz' algorithm (MATLAB® reports an infinite
 *       eigenvalue as Inf; here it is counted, not stored).
 * @warning QZ is iterative (data-dependent iteration count). Suitable for design-time
 *          analysis; avoid in hard-real-time loops where WCET matters.
 * @see "An Algorithm for Generalized Matrix Eigenvalue Problems" (Moler & Stewart, 1973), SIAM J.
 *      Numer. Anal. 10(2)
 * @see "Matrix Computations" (Golub & Van Loan, 2013), §7.7
 *
 * @param A Left matrix of the pencil
 * @param B Right matrix of the pencil
 * @return GeneralizedEigenResult; success is false for a non-finite entry, when the iteration
 *         budget runs out (the eigenvalues found by then are kept), or when a finite eigenvalue
 *         is too large for T
 */
template<size_t N, typename T = double>
    requires damp::is_floating_point_v<T> // the real QZ: real A and B only
[[nodiscard]] constexpr GeneralizedEigenResult<N, T> generalized_eigenvalues(const Matrix<N, N, T>& A, const Matrix<N, N, T>& B) {
    GeneralizedEigenResult<N, T> out{};
    const auto                   finite_entry = [](T v) { return damp::isfinite(v); };
    if (!mat::all_of(A, finite_entry) || !mat::all_of(B, finite_entry)) {
        return out;
    }
    Matrix<N, N, T>                  H = A;
    Matrix<N, N, T>                  R = B;
    damp::array<damp::complex<T>, N> lambda{};
    damp::array<bool, N>             finite{};
    out.success = detail::pencil_eigenvalues(H, R, N, lambda, finite);
    for (size_t i = 0; i < N; ++i) {
        if (finite[i]) {
            out.values[out.n_finite] = lambda[i];
            ++out.n_finite;
        }
    }
    return out;
}

} // namespace mat
} // namespace damp
