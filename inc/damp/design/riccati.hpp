// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file riccati.hpp
 * @brief Algebraic Riccati equation solvers (DARE/CARE)
 */

#include <cstddef>
#include <cstdint>
#include <limits>

#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/matrix.hpp"

namespace damp {

/**
 * @brief Check if (A, B) is a stabilizable pair
 *
 * (A, B) is stabilizable iff every uncontrollable eigenvalue of A lies strictly
 * inside the unit circle. An eigenvalue λ is uncontrollable when
 * rank([λI − A, B]) < n.
 *
 * Stabilizability is weaker than controllability — it permits uncontrollable
 * modes as long as they are already stable (|λ| < 1).
 *
 * @see "Optimal Control" (Anderson & Moore, 1990), §2.4
 *
 * @tparam NX Number of states
 * @tparam NU Number of inputs
 * @tparam T  Scalar type
 * @param A   State matrix (NX × NX)
 * @param B   Input matrix (NX × NU)
 * @return true if (A, B) is stabilizable
 */
template<size_t NX, size_t NU, typename T = double>
constexpr bool is_stabilizable(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B
) {
    using Cplx = damp::complex<T>;

    auto eigen = mat::compute_eigenvalues(A);
    if (!eigen.converged) {
        return false;
    }

    const T tol = std::is_same_v<T, float> ? static_cast<T>(1e-5) : static_cast<T>(1e-10);

    for (size_t i = 0; i < NX; ++i) {
        // Only check unstable eigenvalues (|λ| >= 1)
        if (eigen.values[i].abs() < T{1}) {
            continue;
        }

        // Form [λI - A, B] in complex arithmetic and check rank.
        Matrix<NX, NX + NU, Cplx> test_mat{};

        // λI - A
        for (size_t r = 0; r < NX; ++r) {
            for (size_t c = 0; c < NX; ++c) {
                test_mat(r, c) = Cplx(-A(r, c), T{0});
            }
            test_mat(r, r) = test_mat(r, r) + eigen.values[i];
        }

        // B
        for (size_t r = 0; r < NX; ++r) {
            for (size_t c = 0; c < NU; ++c) {
                test_mat(r, NX + c) = Cplx(B(r, c), T{0});
            }
        }

        // Uncontrollable unstable mode ⇒ not stabilizable.
        if (mat::rank(test_mat, tol) < NX) {
            return false;
        }
    }

    return true;
}

namespace detail {

/**
 * @brief Solve DARE via Structure-Preserving Doubling Algorithm (SDA)
 *
 * Solves AᵀXA − X − (AᵀXB + N)(R + BᵀXB)⁻¹(BᵀXA + Nᵀ) + Q = 0
 *
 * Quadratic convergence (doubles correct digits each iteration).
 * Requires R positive definite (needs R⁻¹ for cross-term reduction).
 * No precondition checks — use dare() for validated entry point.
 *
 * @see Chu et al., "Structure-Preserving Algorithms for Periodic DRE" (2004)
 * @see "Optimal Control" (Anderson & Moore, 1990), §4.3
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>> dare_sda(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N = Matrix<NX, NU, T>{}
) {
    //! Handle cross-term N by reducing to standard form:
    //!   A_bar = A − B R⁻¹ Nᵀ,  Q_bar = Q − N R⁻¹ Nᵀ
    //! Solve via R * X = Nᵀ  and  R * X = Bᵀ  to avoid explicit R⁻¹
    const auto R_inv_Nt_opt = mat::lu_solve(R, N.transpose());
    if (!R_inv_Nt_opt) {
        return damp::nullopt;
    }
    const Matrix R_inv_Nt = R_inv_Nt_opt.value();

    const Matrix A_eff = A - B * R_inv_Nt;
    const Matrix Q_eff = Q - N * R_inv_Nt;

    //! Compute Gk = B R⁻¹ Bᵀ via solve: R * X = Bᵀ → X = R⁻¹Bᵀ → Gk = B * X
    const auto R_inv_Bt_opt = mat::lu_solve(R, B.transpose());
    if (!R_inv_Bt_opt) {
        return damp::nullopt;
    }

    //! Structure-preserving Doubling Algorithm (SDA) for DARE:
    //!   Initialize:  Ak = A,  Gk = B R⁻¹ Bᵀ,  Hk = Q
    //!   Iterate:
    //!     Vk = (I + Gk Hk)⁻¹  (solved via LU)
    //!     Ak₊₁ = Ak Vk Ak
    //!     Gk₊₁ = Gk + Ak Vk Gk Akᵀ
    //!     Hk₊₁ = Hk + Akᵀ Hk Vk Ak
    //!   Converges quadratically: Hk → X for any stabilizable/detectable pair.

    Matrix<NX, NX, T> Ak = A_eff;
    Matrix<NX, NX, T> Gk = B * R_inv_Bt_opt.value();
    Matrix<NX, NX, T> Hk = Q_eff;

    const T   tol = default_tol<T>(); // 1e-6f / 1e-12 — float DARE must be able to succeed
    const int max_iter = 100;         //! Quadratic convergence needs far fewer iterations
    bool      converged = false;

    for (int iter = 0; iter < max_iter; ++iter) {
        //! Solve (I + Gk Hk) Vk = I via LU decomposition (recommended by SDA paper)
        const Matrix Vk_arg = Matrix<NX, NX, T>::identity() + Gk * Hk;
        const auto   Vk_opt = mat::lu_solve(Vk_arg, Matrix<NX, NX, T>::identity());
        if (!Vk_opt) {
            return damp::nullopt;
        }
        const Matrix Vk = Vk_opt.value();

        const Matrix Ak_next = Ak * Vk * Ak;
        const Matrix Gk_next = Gk + Ak * Vk * Gk * Ak.t();
        const Matrix Hk_next = Hk + Ak.t() * Hk * Vk * Ak;

        const T diff = (Hk_next - Hk).norm();
        Ak = Ak_next;
        Gk = Gk_next;
        Hk = Hk_next;

        if (diff < tol * damp::max(T{1}, Hk.norm())) {
            converged = true;
            break;
        }
    }

    if (!converged) {
        return damp::nullopt;
    }

    //! Symmetrize for numerical cleanup
    return (Hk + Hk.t()) * static_cast<T>(0.5);
}

/**
 * @brief Solve DARE via Riccati Difference Equation (RDE) iteration
 *
 * Iterates the discrete Riccati recursion to steady state:
 *
 *     X[k+1] = AᵀX[k]A + Q − AᵀX[k]B(R + BᵀX[k]B)⁻¹BᵀX[k]A
 *
 * Linear convergence. Handles R ≥ 0 (only requires R + BᵀXB invertible at
 * each step, not R itself). Useful when R is singular (e.g., cheap-control
 * problems or minimum-energy estimation).
 *
 * No precondition checks — use dare() for validated entry point.
 *
 * @see "Optimal Control" (Anderson & Moore, 1990), §4.2
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>> dare_rde(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N = Matrix<NX, NU, T>{}
) {
    //! Handle cross-term N by reducing to standard form (requires R invertible)
    Matrix<NX, NX, T> A_eff = A;
    Matrix<NX, NX, T> Q_eff = Q;

    const T n_tol = default_tol<T>();
    bool    n_is_zero = true;
    for (size_t i = 0; i < NX && n_is_zero; ++i) {
        for (size_t j = 0; j < NU; ++j) {
            if (damp::abs(N(i, j)) > n_tol) {
                n_is_zero = false;
                break;
            }
        }
    }

    if (!n_is_zero) {
        const auto R_inv_Nt_opt = mat::lu_solve(R, N.transpose());
        if (!R_inv_Nt_opt) {
            return damp::nullopt; // singular R with non-zero N is unsupported
        }
        const Matrix R_inv_Nt = R_inv_Nt_opt.value();
        A_eff = A - B * R_inv_Nt;
        Q_eff = Q - N * R_inv_Nt;
    }

    //! Riccati Difference Equation iteration
    Matrix<NX, NX, T> X = Q_eff;
    const T           tol = default_tol<T>(); // 1e-6f / 1e-12
    const int         max_iter = 500;
    // Type-safe ceiling: float cannot represent 1e150 (narrowing under -Wnarrowing).
    const T guard = std::numeric_limits<T>::max() / T{100};

    bool converged = false;

    for (int iter = 0; iter < max_iter; ++iter) {
        const Matrix BtX = B.t() * X;
        const Matrix BtXA = BtX * A_eff;
        const Matrix S = R + BtX * B;

        const auto M_opt = mat::lu_solve(S, BtXA);
        if (!M_opt) {
            //! First-iteration kick: if S = R + B'QB is singular, try X₀ = Q + αI
            if (iter == 0) {
                T trace_q = T{0};
                for (size_t i = 0; i < NX; ++i) {
                    trace_q += damp::abs(Q_eff(i, i));
                }
                X = Q_eff + Matrix<NX, NX, T>::identity() * (trace_q / T(NX) + T{1});
                continue;
            }
            return damp::nullopt;
        }
        const Matrix M = M_opt.value();

        const Matrix X_next = A_eff.t() * X * A_eff + Q_eff
                            - A_eff.t() * X * B * M;

        //! Divergence guard (matches dlyap pattern)
        bool diverged = false;
        for (size_t i = 0; i < NX && !diverged; ++i) {
            for (size_t j = 0; j < NX; ++j) {
                if (!damp::isfinite(X_next(i, j)) || damp::abs(X_next(i, j)) > guard) {
                    diverged = true;
                    break;
                }
            }
        }
        if (diverged) {
            return damp::nullopt;
        }

        //! Convergence check (Frobenius norm with guard clamping)
        T diff_norm_sq = T{0};
        T x_norm_sq = T{0};
        for (size_t i = 0; i < NX; ++i) {
            for (size_t j = 0; j < NX; ++j) {
                T diff = X_next(i, j) - X(i, j);
                if (diff > guard) {
                    diff = guard;
                } else if (diff < -guard) {
                    diff = -guard;
                }
                diff_norm_sq += diff * diff;

                T xval = X_next(i, j);
                if (xval > guard) {
                    xval = guard;
                } else if (xval < -guard) {
                    xval = -guard;
                }
                x_norm_sq += xval * xval;
            }
        }

        X = X_next;

        if (damp::sqrt(diff_norm_sq) < tol * damp::max(T{1}, damp::sqrt(x_norm_sq))) {
            converged = true;
            break;
        }
    }

    if (!converged) {
        return damp::nullopt;
    }

    return (X + X.t()) * static_cast<T>(0.5);
}

/**
 * @brief Split a real-eigenvalue 2×2 Schur block into two 1×1 blocks.
 *
 * Francis QR leaves any 2×2 diagonal block with @e real eigenvalues
 * untriangularized (it only records the eigenvalues). This applies the Givens
 * rotation whose first column is the dominant eigenvector, driving the
 * subdiagonal to zero and accumulating the transform into the Schur vectors
 * @p Z. Genuine complex-conjugate pairs (negative discriminant) are left intact.
 *
 * The block is @f$ \bigl[\begin{smallmatrix} a & b \\ c & d \end{smallmatrix}\bigr] @f$
 * at diagonal offset @p i; its eigenvalues are real iff
 * @f$ ((a-d)/2)^2 + bc \ge 0 @f$.
 *
 * @see LAPACK dlanv2 (the standardization step of the real Schur form)
 */
template<size_t M, typename T>
constexpr void split_real_2x2(Matrix<M, M, T>& Tm, Matrix<M, M, T>& Z, size_t i) {
    const T a = Tm(i, i);
    const T b = Tm(i, i + 1);
    const T c = Tm(i + 1, i);
    const T d = Tm(i + 1, i + 1);
    const T p = static_cast<T>(0.5) * (a - d);
    const T disc = (p * p) + (b * c);
    if (disc <= T{0}) {
        return; //! Complex conjugate pair — keep the 2×2 block.
    }
    //! Dominant eigenvalue and its eigenvector v ∝ [b, λ−a] (fall back to [λ−d, c]).
    const T lambda = (static_cast<T>(0.5) * (a + d)) + damp::copysign(damp::sqrt(disc), p);
    T       v0 = b;
    T       v1 = lambda - a;
    if ((damp::abs(v0) + damp::abs(v1)) == T{0}) {
        v0 = lambda - d;
        v1 = c;
    }
    const T nrm = damp::sqrt((v0 * v0) + (v1 * v1));
    if (nrm == T{0}) {
        return;
    }
    const T cs = v0 / nrm;
    const T sn = v1 / nrm;
    //! Givens G = [[cs,−sn],[sn,cs]] with first column = eigenvector ⇒ GᵀTmG upper-triangular.
    for (size_t k = 0; k < M; ++k) {
        const T r0 = Tm(i, k);
        const T r1 = Tm(i + 1, k);
        Tm(i, k) = (cs * r0) + (sn * r1);
        Tm(i + 1, k) = (-sn * r0) + (cs * r1);
    }
    for (size_t k = 0; k < M; ++k) {
        const T c0 = Tm(k, i);
        const T c1 = Tm(k, i + 1);
        Tm(k, i) = (c0 * cs) + (c1 * sn);
        Tm(k, i + 1) = (-c0 * sn) + (c1 * cs);
        const T z0 = Z(k, i);
        const T z1 = Z(k, i + 1);
        Z(k, i) = (z0 * cs) + (z1 * sn);
        Z(k, i + 1) = (-z0 * sn) + (z1 * cs);
    }
    Tm(i + 1, i) = T{0};
}

/**
 * @brief Swap two adjacent diagonal blocks of a real Schur form.
 *
 * Exchanges the P×P block A immediately above-left of the Q×Q block C at diagonal
 * offset @p j (so C's eigenvalues end up first), via an orthogonal similarity that
 * is accumulated into the Schur vectors @p Z. The C-invariant subspace within the
 * window is @f$ \mathrm{span}\bigl(\bigl[\begin{smallmatrix} X \\ I \end{smallmatrix}\bigr]\bigr) @f$,
 * where X solves the Sylvester equation A·X − X·C = −B (B the P×Q coupling block);
 * orthonormalizing it by QR and using it as the leading columns performs the swap.
 * The Sylvester system is tiny (≤ 4×4 via a Kronecker expansion). P,Q ∈ {1,2}.
 *
 * @return false if the Sylvester solve is singular (the blocks share an eigenvalue).
 * @see Golub & Van Loan §7.6.2; Bai & Demmel, "On swapping diagonal blocks" (1993)
 */
template<size_t P, size_t Q, size_t M, typename T>
constexpr bool swap_schur_blocks(Matrix<M, M, T>& Tm, Matrix<M, M, T>& Z, size_t j) {
    constexpr size_t S = P + Q;
    constexpr size_t PQ = P * Q;

    const Matrix<P, P, T> A = Tm.template block<P, P>(j, j);
    const Matrix<Q, Q, T> C = Tm.template block<Q, Q>(j + P, j + P);
    const Matrix<P, Q, T> B = Tm.template block<P, Q>(j, j + P);

    //! Sylvester A·X − X·C = −B via Kronecker: [(I_Q⊗A) − (Cᵀ⊗I_P)]·vec(X) = −vec(B)
    //! (column-major vec). PQ ≤ 4 — a tiny dense LU solve.
    Matrix<PQ, PQ, T> K = Matrix<PQ, PQ, T>::zeros();
    Matrix<PQ, 1, T>  rhs{};
    for (size_t cc = 0; cc < Q; ++cc) {
        for (size_t rr = 0; rr < P; ++rr) {
            const size_t out = rr + (cc * P);
            rhs(out, 0) = -B(rr, cc);
            for (size_t k = 0; k < P; ++k) {
                K(out, k + (cc * P)) += A(rr, k);
            }
            for (size_t k = 0; k < Q; ++k) {
                K(out, rr + (k * P)) -= C(k, cc);
            }
        }
    }
    const auto x_opt = mat::lu_solve(K, rhs);
    if (!x_opt) {
        return false;
    }

    //! Orthonormal basis of the C-invariant subspace [[X];[I_Q]] (S×Q) via QR.
    Matrix<S, Q, T> Mq{};
    for (size_t cc = 0; cc < Q; ++cc) {
        for (size_t rr = 0; rr < P; ++rr) {
            Mq(rr, cc) = x_opt.value()(rr + (cc * P), 0);
        }
        Mq(P + cc, cc) = T{1};
    }
    const Matrix<S, S, T> G = mat::full_qr(Mq).Q;

    //! Apply the similarity Gᵀ(·)G to the affected rows/columns, and G to Z.
    const Matrix<S, M, T> row_blk = Tm.template block<S, M>(j, 0).to_matrix();
    Tm.template block<S, M>(j, 0) = G.transpose() * row_blk;
    const Matrix<M, S, T> col_blk = Tm.template block<M, S>(0, j).to_matrix();
    Tm.template block<M, S>(0, j) = col_blk * G;
    const Matrix<M, S, T> z_blk = Z.template block<M, S>(0, j).to_matrix();
    Z.template block<M, S>(0, j) = z_blk * G;
    return true;
}

/**
 * @brief Reorder a real Schur form so eigenvalues satisfying @p in_front lead.
 *
 * First standardizes the form (split_real_2x2) so every remaining 2×2 block is a
 * genuine complex pair, then bubbles each block whose eigenvalue real part
 * satisfies the predicate to the top-left via adjacent-block swaps, keeping the
 * Schur vectors @p Z orthogonal. For the CARE Hamiltonian the predicate selects
 * the stable spectrum (Re λ < 0), collecting the stabilizing invariant subspace
 * into the leading columns. Best-effort: a singular swap is skipped (the caller's
 * subsequent solve detects the resulting rank deficiency).
 *
 * @see LAPACK dtrsen / dtrexc
 */
template<size_t M, typename T, typename Pred>
constexpr void reorder_schur(Matrix<M, M, T>& Tm, Matrix<M, M, T>& Z, Pred in_front) {
    constexpr T eps = std::numeric_limits<T>::epsilon();
    T           anorm = T{0};
    for (size_t i = 0; i < M; ++i) {
        anorm += damp::abs(Tm(i, i));
    }
    const T tol = eps * damp::max(anorm, T{1});

    const auto is_2x2 = [&](size_t i) {
        return (i + 1 < M) && (damp::abs(Tm(i + 1, i)) > tol);
    };

    //! Pass 1: standardize — split real 2×2 blocks into 1×1 blocks.
    for (size_t i = 0; i + 1 < M;) {
        if (is_2x2(i)) {
            split_real_2x2(Tm, Z, i);
            i += is_2x2(i) ? size_t{2} : size_t{1};
        } else {
            ++i;
        }
    }

    //! Pass 2: bubble selected eigenvalues to the front. Invariant: [top, i) holds
    //! the already-passed unselected blocks.
    size_t top = 0;
    size_t i = 0;
    while (i < M) {
        const size_t b = is_2x2(i) ? 2 : 1;
        const T      re = (b == 2) ? (static_cast<T>(0.5) * (Tm(i, i) + Tm(i + 1, i + 1))) : Tm(i, i);
        if (in_front(re)) {
            size_t cur = i;
            while (cur > top) {
                const size_t pa = (cur >= 2 && damp::abs(Tm(cur - 1, cur - 2)) > tol) ? 2 : 1;
                const size_t jj = cur - pa;
                bool         ok = true;
                if (pa == 1 && b == 1) {
                    ok = swap_schur_blocks<1, 1>(Tm, Z, jj);
                } else if (pa == 1 && b == 2) {
                    ok = swap_schur_blocks<1, 2>(Tm, Z, jj);
                } else if (pa == 2 && b == 1) {
                    ok = swap_schur_blocks<2, 1>(Tm, Z, jj);
                } else {
                    ok = swap_schur_blocks<2, 2>(Tm, Z, jj);
                }
                if (!ok) {
                    break;
                }
                cur = jj;
            }
            top += b;
        }
        i += b;
    }
}

/**
 * @brief Solve CARE via the ordered real-Schur method (Laub's method).
 *
 * Solves AᵀX + XA − (XB + N)R⁻¹(BᵀX + Nᵀ) + Q = 0 from the stable invariant
 * subspace of the Hamiltonian
 *
 *     H = ⎡  A   −G  ⎤,   G = B R⁻¹ Bᵀ   (cross-term N folded into A, Q first)
 *         ⎣ −Q   −Aᵀ ⎦
 *
 * H is reduced to real Schur form Zᵀ H Z = T (Hessenberg + Francis double-shift
 * QR), the spectrum is reordered so the NX stable eigenvalues (Re λ < 0) lead, and
 * the leading NX Schur vectors [U₁₁; U₂₁] span the stabilizing subspace. The
 * solution is X = U₂₁ U₁₁⁻¹ (obtained as the transpose of U₁₁ᵀ Xᵀ = U₂₁ᵀ).
 *
 * Built entirely from orthogonal transforms, so unlike the matrix-sign-function
 * iteration it stays accurate near the imaginary axis. No precondition checks —
 * use care() for the validated entry point. Returns nullopt if the QR iteration
 * fails to converge or the stabilizing subspace is rank-deficient (e.g. H has
 * eigenvalues on the imaginary axis, i.e. no stabilizing solution exists).
 *
 * @see Laub, "A Schur method for solving algebraic Riccati equations," IEEE TAC
 *      1979, https://doi.org/10.1109/TAC.1979.1102178
 * @see Golub & Van Loan, "Matrix Computations" §7.6 (ordered Schur form)
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>> care_schur(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N = Matrix<NX, NU, T>{}
) {
    //! G = B R⁻¹ Bᵀ and cross-term reduction A_eff = A − B R⁻¹ Nᵀ,
    //! Q_eff = Q − N R⁻¹ Nᵀ — solved against R rather than forming R⁻¹.
    const auto Rinv_Bt_opt = mat::lu_solve(R, B.transpose());
    if (!Rinv_Bt_opt) {
        return damp::nullopt;
    }
    const auto Rinv_Nt_opt = mat::lu_solve(R, N.transpose());
    if (!Rinv_Nt_opt) {
        return damp::nullopt;
    }
    const Matrix<NX, NX, T> G = B * Rinv_Bt_opt.value();
    const Matrix<NX, NX, T> A_eff = A - B * Rinv_Nt_opt.value();
    const Matrix<NX, NX, T> Q_eff = Q - N * Rinv_Nt_opt.value();

    //! Assemble the 2NX×2NX Hamiltonian H = [[A_eff, −G], [−Q_eff, −A_effᵀ]].
    constexpr size_t M = 2 * NX;
    Matrix<M, M, T>  H = Matrix<M, M, T>::zeros();
    H.template block<NX, NX>(0, 0) = A_eff;
    H.template block<NX, NX>(0, NX) = G * T{-1};
    H.template block<NX, NX>(NX, 0) = Q_eff * T{-1};
    H.template block<NX, NX>(NX, NX) = A_eff.transpose() * T{-1};

    //! Real Schur form Zᵀ H Z = T, accumulating Schur vectors in Z.
    Matrix<M, M, T> T_schur = H;
    Matrix<M, M, T> Z;
    mat::detail::hessenberg_reduce(T_schur, Z);
    damp::array<T, M> wr{};
    damp::array<T, M> wi{};
    if (!mat::detail::francis_qr(T_schur, Z, wr, wi)) {
        return damp::nullopt;
    }

    //! Reorder the NX stable eigenvalues (Re λ < 0) into the leading block, so the
    //! first NX Schur vectors span the stabilizing invariant subspace.
    reorder_schur(T_schur, Z, [](T re) { return re < T{0}; });

    const Matrix<NX, NX, T> U11 = Z.template block<NX, NX>(0, 0);
    const Matrix<NX, NX, T> U21 = Z.template block<NX, NX>(NX, 0);

    //! X = U₂₁ U₁₁⁻¹, solved as the transpose of U₁₁ᵀ Xᵀ = U₂₁ᵀ (no explicit inverse).
    const auto Xt_opt = mat::lu_solve(U11.transpose(), U21.transpose());
    if (!Xt_opt) {
        return damp::nullopt;
    }
    const Matrix<NX, NX, T> X = Xt_opt.value().transpose();

    //! Symmetrize for numerical cleanup.
    return (X + X.t()) * static_cast<T>(0.5);
}

} // namespace detail

enum class DareMethod : uint8_t { Auto,
                                  SDA,
                                  RDE };

/**
 * @brief Solve the Discrete Algebraic Riccati Equation (DARE)
 *
 * Finds the unique stabilizing solution X to the control-form equation:
 *
 *     AᵀXA − X − (AᵀXB + N)(R + BᵀXB)⁻¹(BᵀXA + Nᵀ) + Q = 0
 *
 * Used directly by discrete LQR (X = S, B = control input matrix).
 *
 * Filter / estimator DARE (dual). The steady-state Kalman covariance P
 * satisfies the filter Riccati. It is the same solver on the dual plant:
 *
 *     dare(Aᵀ, Cᵀ, Q_eff, R_eff)  →  P
 *
 * with Q_eff = G Q Gᵀ (process) and R_eff = H R Hᵀ (measurement). Prefer
 * design::kalman when the deliverable is the estimator gain L as well as P;
 * call this dual form only when the Riccati matrix itself is the goal.
 *
 * Preconditions (checked internally):
 * - Q symmetric positive semidefinite
 * - R symmetric (positive definite for SDA, positive semidefinite for RDE)
 * - (A, B) stabilizable
 *
 * @note Compare with MATLAB®'s idare(A, B, Q, R, N).
 *
 * @see dare_sda() — quadratic convergence, requires R > 0
 * @see dare_rde() — linear convergence, handles R ≥ 0
 * @see design::kalman — filter DARE + steady-state gain L (uses this dual)
 * @see design::discrete_lqr — control DARE + gain K
 * @see "Optimal Control" (Anderson & Moore, 1990), Chapter 4
 *
 * @param A       State transition matrix (NX × NX)
 * @param B       Input matrix (NX × NU); for the filter dual, pass Cᵀ
 * @param Q       State cost / process weight (NX × NX, positive semidefinite)
 * @param R       Input cost / measurement weight (NU × NU; for dual, NY × NY)
 * @param N       Cross-term matrix (NX × NU, default: zero)
 * @param method  DareMethod::Auto (default) selects SDA when R > 0, RDE when R ≥ 0.
 *                DareMethod::SDA forces SDA (requires R positive definite).
 *                DareMethod::RDE forces RDE (handles R positive semidefinite).
 * @return Solution matrix X (NX × NX, positive semidefinite) or damp::nullopt on failure
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>> dare(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N = Matrix<NX, NU, T>{},
    DareMethod               method = DareMethod::Auto
) {
    //! Check R is symmetric
    if (!mat::is_symmetric_or_hermitian(R)) {
        return damp::nullopt;
    }

    //! Check Q is symmetric
    if (!mat::is_symmetric_or_hermitian(Q)) {
        return damp::nullopt;
    }

    //! Check Q is positive semidefinite via Cholesky of (Q + εI)
    //! If Q is PSD, Q + εI is PD for any ε > 0, so Cholesky succeeds.
    //! If Q has a negative eigenvalue, Q + εI will still fail for small ε.
    const T eps = default_tol<T>();
    {
        const Matrix<NX, NX, T> Q_shifted = Q + Matrix<NX, NX, T>::identity() * eps;
        if (!mat::cholesky(Q_shifted)) {
            return damp::nullopt;
        }
    }

    //! Check (A, B) is stabilizable
    if (!is_stabilizable(A, B)) {
        return damp::nullopt;
    }

    const bool r_is_pd = mat::cholesky(R).has_value();
    const bool r_is_psd = r_is_pd || mat::cholesky(R + Matrix<NU, NU, T>::identity() * eps).has_value();

    if (!r_is_psd) {
        return damp::nullopt;
    }

    switch (method) {
        case DareMethod::SDA:
            if (!r_is_pd) {
                return damp::nullopt;
            }
            return detail::dare_sda(A, B, Q, R, N);
        case DareMethod::RDE:
            return detail::dare_rde(A, B, Q, R, N);
        case DareMethod::Auto:
        default:
            if (r_is_pd) {
                return detail::dare_sda(A, B, Q, R, N);
            }
            return detail::dare_rde(A, B, Q, R, N);
    }
}

/**
 * @brief Optimal LQR state-feedback gain from a Riccati solution
 *
 * Given the DARE solution @f$ S @f$, computes the gain @f$ K @f$ for @f$ u = -Kx @f$:
 * @f[
 *   K = (R + B^\top S B)^{-1} (B^\top S A + N^\top).
 * @f]
 * @f$ R + B^\top S B @f$ is symmetric positive definite, so the system is solved
 * by Cholesky factorization rather than forming an explicit inverse — the
 * numerically stabler and cheaper route. Shared by dlqr/lqi/lqg gain synthesis.
 *
 * @see dare() — produces the Riccati solution S
 * @see "Optimal Control" (Anderson & Moore, 1990), §4.3
 *
 * @param A  State transition matrix (NX × NX)
 * @param B  Input matrix (NX × NU)
 * @param S  DARE solution (NX × NX, symmetric positive semidefinite)
 * @param R  Input cost matrix (NU × NU, positive definite)
 * @param N  Cross-term cost matrix (NX × NU, default: zero)
 * @return Optimal gain K (NU × NX) or damp::nullopt if the Cholesky solve fails
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NU, NX, T>> lqr_gain(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& S,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N = Matrix<NX, NU, T>{}
) {
    const Matrix<NU, NU, T> denom = R + B.t() * S * B;
    const Matrix<NU, NX, T> rhs = B.t() * S * A + N.t();
    return mat::cholesky_solve(denom, rhs);
}

/**
 * @brief Solve the Continuous-time Algebraic Riccati Equation (CARE)
 *
 * Finds the unique stabilizing solution X to:
 *
 *     AᵀX + XA − (XB + N)R⁻¹(BᵀX + Nᵀ) + Q = 0
 *
 * This is the continuous-time counterpart of dare(). It underpins continuous
 * LQR/LQG design that has not been discretized first; the discrete pipeline
 * (LQR/LQI/LQG/LQGI) discretizes the plant and uses dare() instead.
 *
 * Preconditions (checked internally):
 * - Q symmetric positive semidefinite
 * - R symmetric positive definite (CARE forms G = B R⁻¹ Bᵀ)
 *
 * Existence of the stabilizing solution additionally requires (A, B)
 * stabilizable and (A, Q) detectable; those are not pre-screened here (the
 * continuous stabilizability/detectability tests differ from the discrete
 * is_stabilizable() used by dare()). Instead, infeasibility surfaces as the
 * Schur reduction failing to converge or the stabilizing subspace being
 * rank-deficient, in which case care() returns damp::nullopt.
 *
 * @note Compare with MATLAB®'s icare(A, B, Q, R) / care(A, B, Q, R).
 *
 * @see care_schur() — the underlying ordered-Schur (Laub's method) solver
 * @see dare() — the discrete-time counterpart
 * @see "Optimal Control" (Anderson & Moore, 1990), §3.3
 *
 * @param A  State matrix (NX × NX)
 * @param B  Input matrix (NX × NU)
 * @param Q  State cost matrix (NX × NX, positive semidefinite)
 * @param R  Input cost matrix (NU × NU, positive definite)
 * @param N  Cross-term matrix (NX × NU, default: zero)
 * @return Stabilizing solution X (NX × NX, symmetric positive semidefinite) or
 *         damp::nullopt on failure
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>> care(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N = Matrix<NX, NU, T>{}
) {
    //! R and Q must be symmetric (CARE assumes symmetric weights).
    if (!mat::is_symmetric_or_hermitian(R)) {
        return damp::nullopt;
    }
    if (!mat::is_symmetric_or_hermitian(Q)) {
        return damp::nullopt;
    }

    //! Q positive semidefinite via Cholesky of (Q + εI) — same trick as dare().
    const T eps = default_tol<T>();
    {
        const Matrix<NX, NX, T> Q_shifted = Q + Matrix<NX, NX, T>::identity() * eps;
        if (!mat::cholesky(Q_shifted)) {
            return damp::nullopt;
        }
    }

    //! R must be positive definite (G = B R⁻¹ Bᵀ).
    if (!mat::cholesky(R)) {
        return damp::nullopt;
    }

    return detail::care_schur(A, B, Q, R, N);
}

} // namespace damp