// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lyapunov.hpp
 * @brief Lyapunov equation solvers: Bartels–Stewart for X, Hammarling for its Cholesky factor
 *
 * Both reduce A once to the complex Schur form A = D·U·S·Uᴴ·D⁻¹ (D the balancing diagonal,
 * U unitary, S upper triangular) and then solve column by column on S, so every pivot is a
 * scalar. The cost is O(NX³) time and O(NX²) storage.
 */

#include <cstddef>
#include <limits>

#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/decomposition.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/matrix.hpp"

namespace damp {

namespace detail {

/// Complex Schur form of a balanced A: A = D·U·S·Uᴴ·D⁻¹ with D = diag(2^e).
template<size_t N, typename T>
struct ComplexSchur {
    Matrix<N, N, damp::complex<T>> S{};       ///< Upper triangular, the eigenvalues on its diagonal
    Matrix<N, N, damp::complex<T>> U{};       ///< Unitary Schur vectors of the balanced matrix
    damp::array<int, N>            e{};       ///< Balancing exponents
    T                              smax{};    ///< max |sᵢⱼ|, the scale for singularity tests
    bool                           ok{false}; ///< The Francis iteration converged
};

/**
 * @brief Complex Schur form from the standardized real Schur form.
 *
 * A standardized complex block @f$ \bigl[\begin{smallmatrix} m & b \\ c & m \end{smallmatrix}\bigr] @f$
 * (bc < 0) has the eigenvector @f$ [b,\; i\omega] @f$ for @f$ \lambda = m + i\omega @f$,
 * @f$ \omega = \sqrt{|b|}\sqrt{|c|} @f$; the unitary G with that first column triangularizes it.
 *
 * @see "Matrix Computations" (Golub & Van Loan, 2013), §7.4.1 (real and complex Schur forms)
 */
template<size_t N, typename T>
[[nodiscard]] constexpr ComplexSchur<N, T> complex_schur(const Matrix<N, N, T>& A) {
    using C = damp::complex<T>;
    ComplexSchur<N, T> out{};
    out.e = mat::detail::balance_exponents(A);
    Matrix<N, N, T> H = mat::detail::balanced(A, out.e);
    Matrix<N, N, T> Z = Matrix<N, N, T>::identity();
    if constexpr (N > 1) {
        mat::detail::hessenberg_reduce(H, Z);
        damp::array<T, N> wr{};
        damp::array<T, N> wi{};
        if (!mat::detail::francis_qr(H, Z, wr, wi)) {
            return out;
        }
    }
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            out.S(i, j) = C{H(i, j), T{0}};
            out.U(i, j) = C{Z(i, j), T{0}};
        }
    }
    for (size_t k = 0; k + 1 < N; ++k) {
        if (H(k + 1, k) == T{0}) {
            continue;
        }
        const T m = H(k, k);
        const T b = H(k, k + 1);
        const T w = damp::sqrt(damp::abs(b)) * damp::sqrt(damp::abs(H(k + 1, k)));
        const T nv = damp::sqrt((b * b) + (w * w));
        const C g00{b / nv, T{0}}; // G = [[g00, g01], [g10, g11]], first column the eigenvector
        const C g10{T{0}, w / nv};
        const C g01 = -damp::conj(g10);
        const C g11 = damp::conj(g00);
        for (size_t j = k; j < N; ++j) { // S ← Gᴴ S
            const C r0 = out.S(k, j);
            const C r1 = out.S(k + 1, j);
            out.S(k, j) = (damp::conj(g00) * r0) + (damp::conj(g10) * r1);
            out.S(k + 1, j) = (damp::conj(g01) * r0) + (damp::conj(g11) * r1);
        }
        for (size_t i = 0; i <= k + 1; ++i) { // S ← S G
            const C c0 = out.S(i, k);
            const C c1 = out.S(i, k + 1);
            out.S(i, k) = (c0 * g00) + (c1 * g10);
            out.S(i, k + 1) = (c0 * g01) + (c1 * g11);
        }
        for (size_t i = 0; i < N; ++i) { // U ← U G
            const C c0 = out.U(i, k);
            const C c1 = out.U(i, k + 1);
            out.U(i, k) = (c0 * g00) + (c1 * g10);
            out.U(i, k + 1) = (c0 * g01) + (c1 * g11);
        }
        out.S(k, k) = C{m, w};
        out.S(k + 1, k + 1) = C{m, -w};
        out.S(k + 1, k) = C{T{0}, T{0}};
        ++k;
    }
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i; j < N; ++j) {
            out.smax = damp::max(out.smax, damp::abs(out.S(i, j)));
        }
    }
    out.ok = true;
    return out;
}

/// Uᴴ·D⁻¹·M·(D⁻¹ on the right when @p both): a real matrix into the Schur coordinates.
template<size_t N, size_t M, typename T>
[[nodiscard]] constexpr Matrix<N, M, damp::complex<T>> to_schur_coordinates(const ComplexSchur<N, T>& cs, const Matrix<N, M, T>& X, bool both) {
    using C = damp::complex<T>;
    Matrix<N, M, C> Xb{};
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            int ex = -cs.e[i];
            if (both) {
                ex -= cs.e[j];
            }
            Xb(i, j) = C{mat::detail::scale_by_pow2(X(i, j), ex), T{0}};
        }
    }
    Matrix<N, M, C> out{};
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            C acc{};
            for (size_t k = 0; k < N; ++k) {
                acc += damp::conj(cs.U(k, i)) * Xb(k, j);
            }
            out(i, j) = acc;
        }
    }
    if (both) {
        if constexpr (N == M) {
            Matrix<N, M, C> right{};
            for (size_t i = 0; i < N; ++i) {
                for (size_t j = 0; j < M; ++j) {
                    C acc{};
                    for (size_t k = 0; k < N; ++k) {
                        acc += out(i, k) * cs.U(k, j);
                    }
                    right(i, j) = acc;
                }
            }
            return right;
        }
    }
    return out;
}

/**
 * @brief Bartels–Stewart on the complex Schur form: S Y + Y Sᴴ = −C or S Y Sᴴ − Y = −C.
 *
 * Column j of Y, last to first, solves an upper triangular system: continuous
 * @f$ (S + \bar\lambda_j I)\,y_j = -c_j - \sum_{k>j} \bar s_{jk}\, y_k @f$; discrete
 * @f$ (\bar\lambda_j S - I)\,y_j = -c_j - S\sum_{k>j} \bar s_{jk}\, y_k @f$. A pivot
 * @f$ \lambda_i + \bar\lambda_j @f$ (or @f$ \lambda_i\bar\lambda_j - 1 @f$) at or below
 * 64·N·ε of its scale means no unique solution.
 */
template<size_t N, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, damp::complex<T>>>
schur_lyapunov(const ComplexSchur<N, T>& cs, const Matrix<N, N, damp::complex<T>>& Cm, bool discrete) {
    using C = damp::complex<T>;
    const T         floor = static_cast<T>(64 * N) * eps<T>();
    const T         tol = discrete ? floor * ((cs.smax * cs.smax) + T{1}) : floor * cs.smax;
    Matrix<N, N, C> Y{};
    for (size_t jj = N; jj-- > 0;) {
        const C           lj = damp::conj(cs.S(jj, jj));
        damp::array<C, N> rhs{};
        damp::array<C, N> w{};
        for (size_t i = 0; i < N; ++i) {
            rhs[i] = -Cm(i, jj);
            for (size_t k = jj + 1; k < N; ++k) {
                w[i] += damp::conj(cs.S(jj, k)) * Y(i, k);
            }
        }
        if (discrete) {
            for (size_t i = 0; i < N; ++i) {
                C sw{};
                for (size_t l = i; l < N; ++l) {
                    sw += cs.S(i, l) * w[l];
                }
                rhs[i] -= sw;
            }
        } else {
            for (size_t i = 0; i < N; ++i) {
                rhs[i] -= w[i];
            }
        }
        for (size_t ii = N; ii-- > 0;) {
            C acc{};
            for (size_t l = ii + 1; l < N; ++l) {
                acc += cs.S(ii, l) * Y(l, jj);
            }
            const C piv = discrete ? ((lj * cs.S(ii, ii)) - C{T{1}, T{0}}) : (cs.S(ii, ii) + lj);
            if (!(damp::abs(piv) > tol)) {
                return damp::nullopt;
            }
            Y(ii, jj) = (rhs[ii] - (discrete ? (lj * acc) : acc)) / piv;
        }
    }
    return Y;
}

/// X = D·Re(U Y Uᴴ)·D, back from the Schur coordinates.
template<size_t N, typename T>
[[nodiscard]] constexpr Matrix<N, N, T> from_schur_coordinates(const ComplexSchur<N, T>& cs, const Matrix<N, N, damp::complex<T>>& Y) {
    using C = damp::complex<T>;
    Matrix<N, N, C> UY{};
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            C acc{};
            for (size_t k = 0; k < N; ++k) {
                acc += cs.U(i, k) * Y(k, j);
            }
            UY(i, j) = acc;
        }
    }
    Matrix<N, N, T> X{};
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            C acc{};
            for (size_t k = 0; k < N; ++k) {
                acc += UY(i, k) * damp::conj(cs.U(j, k));
            }
            X(i, j) = mat::detail::scale_by_pow2(acc.real(), cs.e[i] + cs.e[j]);
        }
    }
    return X;
}

template<size_t N, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> bartels_stewart(const Matrix<N, N, T>& A, const Matrix<N, N, T>& Q, bool discrete) {
    const auto cs = complex_schur(A);
    if (!cs.ok) {
        return damp::nullopt;
    }
    const auto Y = schur_lyapunov(cs, to_schur_coordinates(cs, Q, true), discrete);
    if (!Y) {
        return damp::nullopt;
    }
    return from_schur_coordinates(cs, *Y);
}

/**
 * @brief Hammarling's method on the complex Schur form: the upper triangular R̃ with X̃ = R̃ R̃ᴴ.
 *
 * With S = [[S₁, a], [0, λ]], the last row of the current B written bᴴ, and R̃ = [[R₁, u], [0, μ]]:
 * - continuous: @f$ \mu = \|b\|/\sqrt{-2\,\mathrm{Re}\,\lambda} @f$,
 *   @f$ (S_1 + \bar\lambda I)u = -a\mu - B_1\beta @f$ with @f$ \beta = b/\mu @f$, and the
 *   leading problem continues with @f$ B_1 - u\beta^H @f$;
 * - discrete: @f$ \mu = \|b\|/\sqrt{1-|\lambda|^2} @f$,
 *   @f$ (\bar\lambda S_1 - I)u = -a\bar\lambda\mu - B_1\beta @f$, and with
 *   @f$ y = S_1 u + a\mu @f$, @f$ s = \|\beta\| @f$, @f$ v = B_1\beta/s @f$ the leading
 *   problem continues with @f$ [\,s y - \lambda v,\; B_1 H_{:,2:m}\,] @f$, H the Householder
 *   reflector that takes β/s to a multiple of e₁ (its other columns span β⊥).
 * β is formed as (b/‖b‖)·√(−2 Re λ) (or ·√(1−|λ|²)), bounded however small ‖b‖ is.
 *
 * @see "Numerical Solution of the Stable, Non-negative Definite Lyapunov Equation" (Hammarling, 1982), IMA J. Numer. Anal. 2(3), §5–6
 */
template<size_t N, size_t M, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, damp::complex<T>>>
hammarling(const ComplexSchur<N, T>& cs, Matrix<N, M, damp::complex<T>> Bw, bool discrete) {
    using C = damp::complex<T>;
    Matrix<N, N, C> R{};
    for (size_t kk = N; kk-- > 0;) {
        const C lam = cs.S(kk, kk);
        T       denom{};
        if (discrete) {
            const T mag = damp::abs(lam);
            if (!(mag < T{1})) {
                return damp::nullopt;
            }
            denom = damp::sqrt((T{1} - mag) * (T{1} + mag));
        } else {
            if (!(lam.real() < T{0})) {
                return damp::nullopt;
            }
            denom = damp::sqrt(T{-2} * lam.real());
        }
        T bmax = T{0};
        for (size_t j = 0; j < M; ++j) {
            bmax = damp::max(bmax, damp::abs(Bw(kk, j)));
        }
        if (bmax == T{0}) {
            continue; // μ = 0, u = 0: the leading problem keeps B₁
        }
        T ssq = T{0};
        for (size_t j = 0; j < M; ++j) {
            const T r = damp::abs(Bw(kk, j)) / bmax;
            ssq += r * r;
        }
        const T bnorm = bmax * damp::sqrt(ssq);
        const T mu = bnorm / denom;
        R(kk, kk) = C{mu, T{0}};
        damp::array<C, M> brow{}; // βᴴ
        for (size_t j = 0; j < M; ++j) {
            brow[j] = Bw(kk, j) * (denom / bnorm);
        }
        damp::array<C, N> bb{}; // B₁β
        for (size_t i = 0; i < kk; ++i) {
            for (size_t j = 0; j < M; ++j) {
                bb[i] += Bw(i, j) * damp::conj(brow[j]);
            }
        }
        const C           lbar = damp::conj(lam);
        damp::array<C, N> u{};
        for (size_t ii = kk; ii-- > 0;) {
            C acc{};
            for (size_t l = ii + 1; l < kk; ++l) {
                acc += cs.S(ii, l) * u[l];
            }
            if (discrete) {
                const C rhs = -(cs.S(ii, kk) * lbar * mu) - bb[ii];
                u[ii] = (rhs - (lbar * acc)) / ((lbar * cs.S(ii, ii)) - C{T{1}, T{0}});
            } else {
                const C rhs = -(cs.S(ii, kk) * mu) - bb[ii];
                u[ii] = (rhs - acc) / (cs.S(ii, ii) + lbar);
            }
            R(ii, kk) = u[ii];
        }
        if (!discrete) {
            for (size_t i = 0; i < kk; ++i) {
                for (size_t j = 0; j < M; ++j) {
                    Bw(i, j) -= u[i] * brow[j];
                }
            }
            continue;
        }
        // Discrete: B₁ ← [s y − λ v, B₁ H(:, 2:m)].
        const T           s = denom;
        damp::array<C, N> wcol{};
        for (size_t i = 0; i < kk; ++i) {
            C y = cs.S(i, kk) * mu;
            for (size_t l = i; l < kk; ++l) {
                y += cs.S(i, l) * u[l];
            }
            wcol[i] = (y * s) - (lam * (bb[i] / s));
        }
        if constexpr (M > 1) {
            damp::array<C, M> h{}; // H = I − 2hhᴴ/(hᴴh), H·β̂ ∝ e₁, β̂ = β/s
            for (size_t j = 0; j < M; ++j) {
                h[j] = damp::conj(brow[j]) / s;
            }
            const T x0 = damp::abs(h[0]);
            const C alpha = (x0 > T{0}) ? (-(h[0] / x0)) : C{T{-1}, T{0}};
            h[0] -= alpha;
            T hh = T{0};
            for (size_t j = 0; j < M; ++j) {
                hh += damp::abs(h[j]) * damp::abs(h[j]);
            }
            for (size_t i = 0; i < kk; ++i) {
                C t{};
                for (size_t j = 0; j < M; ++j) {
                    t += Bw(i, j) * h[j];
                }
                const C f = t * (T{2} / hh);
                for (size_t j = 1; j < M; ++j) {
                    Bw(i, j) -= f * damp::conj(h[j]);
                }
            }
        }
        for (size_t i = 0; i < kk; ++i) {
            Bw(i, 0) = wcol[i];
        }
    }
    return R;
}

/// The real upper triangular R with RᵀR = G Gᴴ, G = D·U·R̃, by QR of [Re Gᵀ; Im Gᵀ].
template<size_t N, typename T>
[[nodiscard]] constexpr Matrix<N, N, T> real_cholesky_factor(const ComplexSchur<N, T>& cs, const Matrix<N, N, damp::complex<T>>& Rt) {
    using C = damp::complex<T>;
    Matrix<2 * N, N, T> F{};
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            C acc{};
            for (size_t k = 0; k <= j; ++k) {
                acc += cs.U(i, k) * Rt(k, j);
            }
            const T sr = mat::detail::scale_by_pow2(acc.real(), cs.e[i]);
            const T si = mat::detail::scale_by_pow2(acc.imag(), cs.e[i]);
            F(j, i) = sr; // Gᵀ: row j of Gᵀ is column j of G
            F(N + j, i) = si;
        }
    }
    const auto      qr = mat::full_qr(F);
    Matrix<N, N, T> R{};
    for (size_t i = 0; i < N; ++i) {
        const T sign = (qr.R(i, i) < T{0}) ? T{-1} : T{1};
        for (size_t j = i; j < N; ++j) {
            R(i, j) = sign * qr.R(i, j);
        }
    }
    return R;
}

template<size_t N, size_t M, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> lyapunov_cholesky_impl(const Matrix<N, N, T>& A, const Matrix<N, M, T>& B, bool discrete) {
    const auto cs = complex_schur(A);
    if (!cs.ok) {
        return damp::nullopt;
    }
    const auto Rt = hammarling(cs, to_schur_coordinates(cs, B, false), discrete);
    if (!Rt) {
        return damp::nullopt;
    }
    return real_cholesky_factor(cs, *Rt);
}

} // namespace detail

/**
 * @brief Solve the continuous-time Lyapunov equation A X + X Aᵀ + Q = 0.
 *
 * The unique solution exists iff @f$ A @f$ and @f$ -A @f$ share no eigenvalue
 * (i.e. @f$ \lambda_i + \lambda_j \neq 0 @f$ for all @f$ i,j @f$). For a Hurwitz
 * @f$ A @f$ and symmetric positive-semidefinite @f$ Q @f$, @f$ X @f$ is symmetric
 * positive semidefinite and equals @f$ \int_0^\infty e^{At} Q\, e^{A^\top t}\,dt @f$
 * — the continuous controllability/observability Gramian when @f$ Q = BB^\top @f$
 * / @f$ C^\top C @f$ (see stability::controllability_gramian). Q need not be symmetric.
 *
 * Bartels–Stewart: A is balanced and reduced to complex Schur form, the equation is
 * solved column by column on the triangular factor, and X is mapped back. O(NX³).
 *
 * @note Compare with MATLAB®'s @c X=lyap(A,Q) (which solves @f$ AX+XA^\top+Q=0 @f$).
 *
 * @see dlyap() — the discrete-time counterpart; lyapunov_cholesky() for the factor of X.
 * @see "Solution of the Matrix Equation AX + XB = C" (Bartels & Stewart, 1972), Comm. ACM 15(9)
 * @see "Matrix Computations" (Golub & Van Loan, 2013), §7.6.3.
 *
 * @param A State matrix (NX × NX).
 * @param Q Right-hand side (NX × NX).
 * @return Solution X (NX × NX), or damp::nullopt if A and −A share an eigenvalue (to
 *         64·NX·ε of the Schur form's largest entry) or the Schur reduction fails.
 */
template<size_t NX, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>>
lyap(const Matrix<NX, NX, T>& A, const Matrix<NX, NX, T>& Q) {
    return detail::bartels_stewart(A, Q, false);
}

/**
 * @brief Solve the discrete-time Lyapunov (Stein) equation A X Aᵀ − X + Q = 0.
 *
 * The unique solution exists iff no product of eigenvalues satisfies
 * @f$ \lambda_i \lambda_j = 1 @f$. For a Schur-stable @f$ A @f$ (all
 * @f$ |\lambda| < 1 @f$) and symmetric positive-semidefinite @f$ Q @f$,
 * @f$ X = \sum_{k=0}^{\infty} A^k Q (A^\top)^k @f$ — the discrete
 * controllability/observability Gramian when @f$ Q = BB^\top @f$ / @f$ C^\top C @f$.
 *
 * Bartels–Stewart on the complex Schur form, as lyap(). O(NX³).
 *
 * @note Compare with MATLAB®'s @c X=dlyap(A,Q) (which solves @f$ AXA^\top-X+Q=0 @f$).
 *
 * @see lyap() — the continuous-time counterpart; discrete_lyapunov_cholesky() for the factor of X.
 * @see "Solution of the Matrix Equation AX + XB = C" (Bartels & Stewart, 1972), Comm. ACM 15(9)
 *
 * @param A State transition matrix (NX × NX).
 * @param Q Right-hand side (NX × NX).
 * @return Solution X (NX × NX), or damp::nullopt if some λᵢλⱼ = 1 or the Schur reduction fails.
 */
template<size_t NX, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>>
dlyap(const Matrix<NX, NX, T>& A, const Matrix<NX, NX, T>& Q) {
    return detail::bartels_stewart(A, Q, true);
}

/**
 * @brief Cholesky factor of the continuous Lyapunov solution: A X + X Aᵀ + B Bᵀ = 0, X = RᵀR.
 *
 * Hammarling's method computes R directly, never forming X, for a Hurwitz A. X is a
 * Gramian when B is the input matrix (or Cᵀ with Aᵀ for the observability Gramian), and a
 * Gramian of a nearly uncontrollable system is too ill-conditioned for cholesky(lyap(…))
 * in floating point: its small eigenvalues are lost in X, while R carries them. R is upper
 * triangular with a nonnegative diagonal; B may have any number of columns.
 *
 * @note Compare with MATLAB®'s @c R=lyapchol(A,B).
 * @see "Numerical Solution of the Stable, Non-negative Definite Lyapunov Equation" (Hammarling, 1982), IMA J. Numer. Anal. 2(3)
 *
 * @param A State matrix (NX × NX), Hurwitz.
 * @param B Input matrix (NX × NU).
 * @return R (NX × NX), or damp::nullopt if A is not Hurwitz or the Schur reduction fails.
 */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>>
lyapunov_cholesky(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B) {
    return detail::lyapunov_cholesky_impl(A, B, false);
}

/**
 * @brief Cholesky factor of the discrete Lyapunov solution: A X Aᵀ − X + B Bᵀ = 0, X = RᵀR.
 *
 * Hammarling's method for a Schur-stable A (all |λ| < 1), as lyapunov_cholesky().
 *
 * @note Compare with MATLAB®'s @c R=dlyapchol(A,B).
 * @see "Numerical Solution of the Stable, Non-negative Definite Lyapunov Equation" (Hammarling, 1982), IMA J. Numer. Anal. 2(3), §6
 *
 * @param A State transition matrix (NX × NX), Schur stable.
 * @param B Input matrix (NX × NU).
 * @return R (NX × NX), or damp::nullopt if A is not Schur stable or the Schur reduction fails.
 */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>>
discrete_lyapunov_cholesky(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B) {
    return detail::lyapunov_cholesky_impl(A, B, true);
}

} // namespace damp
