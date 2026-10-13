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
#include "damp/design/lyapunov.hpp"
#include "damp/math/complex.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/matrix.hpp"

namespace damp {

namespace detail {

using mat::detail::gebal_log2_scale;
using mat::detail::scale_by_pow2;

/// (A, B) in balanced coordinates: (D⁻¹ A D, D⁻¹ B), D = diag(2^e) from mat::detail::balance_exponents(A).
template<size_t NX, size_t NU, typename T>
struct BalancedPair {
    damp::array<int, NX> e{}; ///< the scaling exponents
    Matrix<NX, NX, T>    A{}; ///< D⁻¹ A D
    Matrix<NX, NU, T>    B{}; ///< D⁻¹ B
};

/// @copybrief BalancedPair
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr BalancedPair<NX, NU, T> balanced_pair(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B) {
    BalancedPair<NX, NU, T> out{.e = mat::detail::balance_exponents(A)};
    out.A = mat::detail::balanced(A, out.e);
    for (size_t r = 0; r < NX; ++r) {
        for (size_t c = 0; c < NU; ++c) {
            out.B(r, c) = scale_by_pow2(B(r, c), -out.e[r]);
        }
    }
    return out;
}

/**
 * @brief PBH test at one eigenvalue of a balanced pair: true when rank([λI − A, B]) < n
 *        (λ is uncontrollable).
 *
 * @p lambda must come from mat::compute_eigenvalues (which balances the same way): a
 * diagonal similarity does not move an eigenvalue, but it does rescale the roundoff in a
 * computed one, and that roundoff is what the rank tolerance has to absorb.
 *
 * The rank tolerance is relative to the largest singular value, so the test matrix is
 * conditioned first, by operations that do not change its rank. Balancing (a delay
 * register that copies a 10⁵-gain accelerometer row into A otherwise puts that row next
 * to rows of 10⁻³, and the tolerance reads the small ones as zero). Then every B column is
 * brought to A's scale: a large-magnitude input (or, for the dual, a large-gain sensor row)
 * cannot hide the others. Then every row is brought to unit size. Row scaling alone, on an
 * unbalanced pair, amplifies the eigenvalue's roundoff on the small rows by the scale
 * spread and calls an uncontrollable mode controllable. The tolerance is a relative √eps:
 * an uncontrollable eigenvalue is only known to O(√eps) when defective, so a tighter
 * floor would call it controllable.
 */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr bool uncontrollable_at(const BalancedPair<NX, NU, T>& pair, damp::complex<T> lambda) {
    using Cplx = damp::complex<T>;
    const auto& A = pair.A;
    const auto& B = pair.B;
    T           a_scale{0};
    for (size_t r = 0; r < NX; ++r) {
        for (size_t c = 0; c < NX; ++c) {
            a_scale = damp::max(a_scale, damp::abs(A(r, c)));
        }
    }
    a_scale = damp::max(a_scale, T{1});

    Matrix<NX, NX + NU, Cplx> test_mat{};
    for (size_t r = 0; r < NX; ++r) {
        for (size_t c = 0; c < NX; ++c) {
            test_mat(r, c) = Cplx(-A(r, c), T{0});
        }
        test_mat(r, r) = test_mat(r, r) + lambda;
    }
    for (size_t c = 0; c < NU; ++c) {
        T col{0};
        for (size_t r = 0; r < NX; ++r) {
            col = damp::max(col, damp::abs(B(r, c)));
        }
        const T gain = (col > T{0}) ? (a_scale / col) : T{1};
        for (size_t r = 0; r < NX; ++r) {
            test_mat(r, NX + c) = Cplx(B(r, c) * gain, T{0});
        }
    }
    for (size_t r = 0; r < NX; ++r) {
        T row{0};
        for (size_t c = 0; c < NX + NU; ++c) {
            row = damp::max(row, test_mat(r, c).abs());
        }
        if (row > T{0}) {
            for (size_t c = 0; c < NX + NU; ++c) {
                test_mat(r, c) = test_mat(r, c) / Cplx(row, T{0});
            }
        }
    }
    return mat::rank(test_mat, damp::sqrt(eps<T>())) < NX;
}

} // namespace detail

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
 * The pair is balanced first (mat::detail::balance_exponents), so the verdict does not
 * depend on the units of the states. An eigenvalue within √eps of the unit circle counts as
 * unstable: one on the circle (an integrator, a random-walk disturbance) is computed to
 * roundoff on either side, and the DARE certificate uses the same band.
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
    const auto eigen = mat::compute_eigenvalues(A); // balanced internally
    const auto bal = detail::balanced_pair(A, B);
    if (!eigen.converged) {
        return false;
    }
    const T slack = damp::sqrt(eps<T>());
    for (size_t i = 0; i < NX; ++i) {
        // Only unstable eigenvalues (|λ| ≥ 1 − √eps) matter. The conjugate of a complex
        // eigenvalue gives the same rank, so test one of each pair.
        const auto lambda = eigen.values[i];
        if (lambda.abs() < T{1} - slack || lambda.imag() < T{0}) {
            continue;
        }
        if (detail::uncontrollable_at(bal, lambda)) {
            return false; // uncontrollable unstable mode
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

    const T   tol = default_tol<T>(); // relative: 1e-6f / 1e-12 of ‖X‖, so a covariance in m² converges like one in counts²
    const int max_iter = 100;         //! Quadratic convergence needs far fewer iterations
    bool      converged = false;

    for (int iter = 0; iter < max_iter; ++iter) {
        //! Wk = I + Gk Hk, solved for Wk⁻¹Ak and Wk⁻¹Gk — never Wk⁻¹ itself
        //! (six products per step, not ten).
        const Matrix<NX, NX, T> Wk = Matrix<NX, NX, T>::identity() + (Gk * Hk);
        const auto              WinvA_opt = mat::lu_solve(Wk, Ak);
        const auto              WinvG_opt = mat::lu_solve(Wk, Gk);
        if (!WinvA_opt || !WinvG_opt) {
            return damp::nullopt;
        }
        const Matrix<NX, NX, T>& WinvA = WinvA_opt.value();
        const Matrix<NX, NX, T>& WinvG = WinvG_opt.value();
        const Matrix<NX, NX, T>  AkT = Ak.t();

        const Matrix Ak_next = Ak * WinvA;
        const Matrix Gk_next = Gk + ((Ak * WinvG) * AkT);
        const Matrix Hk_next = Hk + ((AkT * Hk) * WinvA);

        const T diff = (Hk_next - Hk).norm();
        Ak = Ak_next;
        Gk = Gk_next;
        Hk = Hk_next;

        if (diff <= tol * Hk.norm()) {
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

    bool n_is_zero = true;
    for (size_t i = 0; i < NX && n_is_zero; ++i) {
        for (size_t j = 0; j < NU; ++j) {
            if (N(i, j) != T{0}) {
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
    const T           tol = default_tol<T>(); // relative to ‖X‖ (see dare_sda)
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

        if (damp::sqrt(diff_norm_sq) <= tol * damp::sqrt(x_norm_sq)) {
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
 * Two 1×1 blocks swap by the Givens rotation whose first column is [b, c − a], which is
 * unconditionally stable. A swap involving a 2×2 block is computed on a copy of the S×S
 * window and kept only if it passes Bai and Demmel's two stability tests: the block that
 * must vanish, rows Q… by columns …Q of GᵀWG, is at most τ (weak), and so is the window
 * rebuilt from the swapped form with that block zeroed, ‖W − G·W′·Gᵀ‖ (strong), entrywise.
 * τ = max(10·S·ε·max|W|, smallest normal): forming G·W′·Gᵀ to check it costs about
 * 2S·ε·max|W| by itself, so a 10·ε threshold would refuse sound swaps on its own rounding.
 * Nearly equal eigenvalues with a large coupling make X huge and the swap inaccurate; it is
 * then refused and the form is left as it was. An accepted swap stores the vanished block
 * as zero and re-standardizes its 2×2 blocks (mat::detail::standardize_schur_block).
 *
 * @return false if the Sylvester solve is singular (the blocks share an eigenvalue) or the
 *         swap fails a stability test; @p Tm and @p Z are then unchanged.
 * @see Golub & Van Loan §7.6.2
 * @see "On swapping diagonal blocks in real Schur form" (Bai & Demmel, 1993), §3 (stability tests)
 */
template<size_t P, size_t Q, size_t M, typename T>
constexpr bool swap_schur_blocks(Matrix<M, M, T>& Tm, Matrix<M, M, T>& Z, size_t j) {
    constexpr size_t S = P + Q;
    constexpr size_t PQ = P * Q;
    if constexpr (P == 1 && Q == 1) {
        //! [[a, b], [0, c]] → [[c, b], [0, a]] by the rotation with first column ∝ [b, c − a].
        const T a = Tm(j, j);
        const T b = Tm(j, j + 1);
        const T c = Tm(j + 1, j + 1);
        if (a == c) {
            return false; //! equal eigenvalues: nothing to order
        }
        const T x0 = b;
        const T x1 = c - a;
        const T r = damp::sqrt((x0 * x0) + (x1 * x1));
        if (!(r > T{0}) || !damp::isfinite(r)) {
            return false;
        }
        const T cs = x0 / r;
        const T sn = x1 / r;
        for (size_t k = j; k < M; ++k) {
            const T r0 = Tm(j, k);
            const T r1 = Tm(j + 1, k);
            Tm(j, k) = (cs * r0) + (sn * r1);
            Tm(j + 1, k) = (cs * r1) - (sn * r0);
        }
        for (size_t k = 0; k <= j + 1; ++k) {
            const T c0 = Tm(k, j);
            const T c1 = Tm(k, j + 1);
            Tm(k, j) = (cs * c0) + (sn * c1);
            Tm(k, j + 1) = (cs * c1) - (sn * c0);
        }
        for (size_t k = 0; k < M; ++k) {
            const T z0 = Z(k, j);
            const T z1 = Z(k, j + 1);
            Z(k, j) = (cs * z0) + (sn * z1);
            Z(k, j + 1) = (cs * z1) - (sn * z0);
        }
        Tm(j, j) = c;
        Tm(j + 1, j + 1) = a;
        Tm(j + 1, j) = T{0};
        return true;
    }
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
    for (size_t k = 0; k < PQ; ++k) {
        if (!damp::isfinite(x_opt.value()(k, 0))) {
            return false;
        }
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

    //! Bai–Demmel: try the swap on the window before touching Tm and Z.
    const Matrix<S, S, T> W = Tm.template block<S, S>(j, j).to_matrix();
    T                     wmax = T{0};
    for (size_t r = 0; r < S; ++r) {
        for (size_t c = 0; c < S; ++c) {
            wmax = damp::max(wmax, damp::abs(W(r, c)));
        }
    }
    const T         tau = damp::max(T{10} * static_cast<T>(S) * eps<T>() * wmax, std::numeric_limits<T>::min());
    Matrix<S, S, T> Ws = G.transpose() * W * G;
    for (size_t r = Q; r < S; ++r) {
        for (size_t c = 0; c < Q; ++c) {
            if (!(damp::abs(Ws(r, c)) <= tau)) {
                return false; //! weak stability: the swapped form is not block upper triangular
            }
            Ws(r, c) = T{0};
        }
    }
    const Matrix<S, S, T> back = G * Ws * G.transpose();
    for (size_t r = 0; r < S; ++r) {
        for (size_t c = 0; c < S; ++c) {
            if (!(damp::abs(back(r, c) - W(r, c)) <= tau)) {
                return false; //! strong stability: the zeroed form is not similar to the window
            }
        }
    }

    //! Apply the similarity Gᵀ(·)G to the affected rows/columns, and G to Z.
    const Matrix<S, M, T> row_blk = Tm.template block<S, M>(j, 0).to_matrix();
    Tm.template block<S, M>(j, 0) = G.transpose() * row_blk;
    const Matrix<M, S, T> col_blk = Tm.template block<M, S>(0, j).to_matrix();
    Tm.template block<M, S>(0, j) = col_blk * G;
    const Matrix<M, S, T> z_blk = Z.template block<M, S>(0, j).to_matrix();
    Z.template block<M, S>(0, j) = z_blk * G;
    for (size_t r = Q; r < S; ++r) {
        for (size_t c = 0; c < Q; ++c) {
            Tm(j + r, j + c) = T{0};
        }
    }
    if constexpr (Q == 2) {
        static_cast<void>(mat::detail::standardize_schur_block(Tm, Z, j));
    }
    if constexpr (P == 2) {
        static_cast<void>(mat::detail::standardize_schur_block(Tm, Z, j + Q));
    }
    return true;
}

/**
 * @brief True when T(i+1, i) is a genuine Schur subdiagonal, not roundoff.
 *
 * Compares against the block's own entries, so the answer does not depend on the units:
 * the Francis iteration, the 2×2 standardization and the block swap all store a vanished
 * subdiagonal as an exact zero, and what remains is either zero or part of a 2×2 block.
 */
template<size_t M, typename T>
constexpr bool schur_subdiagonal(const Matrix<M, M, T>& Tm, size_t i) {
    if (i + 1 >= M) {
        return false;
    }
    const T sub = damp::abs(Tm(i + 1, i));
    const T scale = damp::abs(Tm(i, i)) + damp::abs(Tm(i + 1, i + 1)) + damp::abs(Tm(i, i + 1));
    return sub > eps<T>() * scale || (scale == T{0} && sub > T{0});
}

/**
 * @brief Reorder a real Schur form so eigenvalues satisfying @p in_front lead.
 *
 * First standardizes the form (mat::detail::standardize_schur_block) so every remaining
 * 2×2 block is a genuine complex pair, then bubbles each block whose eigenvalue real part
 * satisfies the predicate to the top-left via adjacent-block swaps, keeping the
 * Schur vectors @p Z orthogonal. For the CARE Hamiltonian the predicate selects
 * the stable spectrum (Re λ < 0), collecting the stabilizing invariant subspace
 * into the leading columns.
 *
 * @return Number of eigenvalues (counting a complex pair as two) that actually
 *         reached the leading block. A swap that fails, or that would reach back
 *         into the already settled prefix, leaves the block where it is and does
 *         not count it — the caller must not treat a short prefix as the invariant
 *         subspace.
 *
 * @see LAPACK dtrsen / dtrexc
 */
template<size_t M, typename T, typename Pred>
constexpr size_t reorder_schur(Matrix<M, M, T>& Tm, Matrix<M, M, T>& Z, Pred in_front) {
    //! Pass 1: standardize — real 2×2 blocks become two 1×1 blocks.
    for (size_t i = 0; i + 1 < M;) {
        if (schur_subdiagonal(Tm, i)) {
            static_cast<void>(mat::detail::standardize_schur_block(Tm, Z, i));
            i += schur_subdiagonal(Tm, i) ? size_t{2} : size_t{1};
        } else {
            ++i;
        }
    }

    //! Pass 2: bubble selected eigenvalues to the front. Invariant: [0, top) holds
    //! the selected blocks that reached the front; [top, i) holds blocks already
    //! passed that were not selected.
    size_t top = 0;
    size_t i = 0;
    while (i < M) {
        const size_t b = schur_subdiagonal(Tm, i) ? 2 : 1;
        const T      re = (b == 2) ? (static_cast<T>(0.5) * (Tm(i, i) + Tm(i + 1, i + 1))) : Tm(i, i);
        if (in_front(re)) {
            size_t cur = i;
            bool   reached = true;
            while (cur > top) {
                const size_t pa = (cur >= 2 && schur_subdiagonal(Tm, cur - 2)) ? 2 : 1;
                const size_t jj = cur - pa;
                if (jj < top) {
                    reached = false;
                    break;
                }
                bool ok = true;
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
                    reached = false;
                    break;
                }
                cur = jj;
            }
            if (reached) {
                top += b;
            }
        }
        i += b;
    }
    return top;
}

/**
 * @brief Symplectic diagonal balance of a 2n×2n Hamiltonian.
 *
 * Scale-only xGEBAL on |H| with a zero diagonal, then the Benner projection
 * sᵢ = round((e_{n+i} − eᵢ) / 2) so the similarity is D = diag(2^s, 2^{−s})
 * and stays symplectic. H is replaced by D H D⁻¹. @p s receives the state
 * exponents; the original CARE solution is recovered by Xᵢⱼ *= 2^{sᵢ+sⱼ}.
 * A non-finite scaled entry leaves H untouched and @p s zero.
 *
 * @see Benner, "Symplectic Balancing of Hamiltonian Matrices," SIAM J. Sci.
 *      Comput. 22(5), 2001, https://doi.org/10.1137/S1064827500367993
 */
template<size_t NX, typename T>
constexpr void balance_hamiltonian(Matrix<2 * NX, 2 * NX, T>& H, damp::array<int, NX>& s) {
    constexpr size_t M = 2 * NX;
    Matrix<M, M, T>  Mag = Matrix<M, M, T>::zeros();
    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < M; ++j) {
            if (i != j) {
                Mag(i, j) = damp::abs(H(i, j));
            }
        }
    }
    damp::array<int, M> e{};
    if (!gebal_log2_scale(Mag, e)) {
        return;
    }
    constexpr int kCap = 60;
    bool          any = false;
    for (size_t i = 0; i < NX; ++i) {
        T half = static_cast<T>(e[NX + i] - e[i]) * static_cast<T>(0.5);
        if (half > T{60}) {
            half = T{60};
        } else if (half < T{-60}) {
            half = T{-60};
        }
        int si = static_cast<int>(damp::nearbyint(half));
        if (si > kCap) {
            si = kCap;
        } else if (si < -kCap) {
            si = -kCap;
        }
        s[i] = si;
        if (si != 0) {
            any = true;
        }
    }
    if (!any) {
        return;
    }
    Matrix<M, M, T> Hb = H;
    for (size_t i = 0; i < M; ++i) {
        const int er = (i < NX) ? s[i] : -s[i - NX];
        for (size_t j = 0; j < M; ++j) {
            const int ec = (j < NX) ? s[j] : -s[j - NX];
            const T   v = scale_by_pow2(H(i, j), er - ec);
            if (!damp::isfinite(v)) {
                for (size_t k = 0; k < NX; ++k) {
                    s[k] = 0;
                }
                return;
            }
            Hb(i, j) = v;
        }
    }
    H = Hb;
}

/**
 * @brief True when the leading @p n eigenvalues of a real Schur form are stable.
 *
 * A block that crosses the cut, or whose real part is positive by more than a few
 * ulps of the largest entry of the form, means the stable subspace was not isolated.
 */
template<size_t M, typename T>
[[nodiscard]] constexpr bool schur_prefix_stable(const Matrix<M, M, T>& Tm, size_t n) {
    T tmax = T{0};
    for (size_t r = 0; r < M; ++r) {
        for (size_t c = 0; c < M; ++c) {
            tmax = damp::max(tmax, damp::abs(Tm(r, c)));
        }
    }
    const T margin = (eps<T>() * T{1024}) * tmax;
    size_t  i = 0;
    while (i < n) {
        const size_t b = schur_subdiagonal(Tm, i) ? size_t{2} : size_t{1};
        if (i + b > n) {
            return false;
        }
        const T re = (b == 2) ? (static_cast<T>(0.5) * (Tm(i, i) + Tm(i + 1, i + 1))) : Tm(i, i);
        if (re > margin) {
            return false;
        }
        i += b;
    }
    return true;
}

/**
 * @brief Laub's symmetry test on the stabilizing basis: U₁₁ᵀ U₂₁ ≈ symmetric.
 *
 * A large skew part means the selected subspace is not Lagrangian (the stable
 * eigenvalues were not separated from their unstable mirrors). The threshold
 * matches SciPy's: max(1000 ε, 0.1 ‖U₁₁ᵀ U₂₁‖), with the Frobenius norm in
 * place of the 1-norm.
 *
 * @see Laub (1979), Theorem 5; SciPy linalg.solve_continuous_are
 */
template<size_t NX, typename T>
[[nodiscard]] constexpr bool care_basis_symmetric(const Matrix<NX, NX, T>& U11, const Matrix<NX, NX, T>& U21) {
    const Matrix<NX, NX, T> gram = U11.transpose() * U21;
    const Matrix<NX, NX, T> skew = gram - gram.transpose();
    const T                 tol = damp::max(eps<T>() * T{1000}, static_cast<T>(0.1) * gram.norm());
    return skew.norm() <= tol;
}

/**
 * @brief X = U₂₁ U₁₁⁻¹ with rows of U₁₁ᵀ equilibrated before the LU solve.
 *
 * lu_decomposition rejects pivots below an absolute tolerance. A stabilizing
 * basis whose entries are all ~1e−14 can be well conditioned and still fail
 * that test. Row scaling leaves X unchanged and brings a uniform basis up to
 * unit row size; a genuinely singular U₁₁ still has a zero row.
 */
template<size_t NX, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>>
care_basis_solve(const Matrix<NX, NX, T>& U11, const Matrix<NX, NX, T>& U21) {
    Matrix<NX, NX, T> Ut = U11.transpose();
    Matrix<NX, NX, T> Rt = U21.transpose();
    for (size_t i = 0; i < NX; ++i) {
        T row_max = T{0};
        for (size_t j = 0; j < NX; ++j) {
            row_max = damp::max(row_max, damp::abs(Ut(i, j)));
        }
        if (!(row_max > T{0}) || !damp::isfinite(row_max)) {
            return damp::nullopt;
        }
        const T inv = T{1} / row_max;
        for (size_t j = 0; j < NX; ++j) {
            Ut(i, j) *= inv;
            Rt(i, j) *= inv;
            if (!damp::isfinite(Ut(i, j)) || !damp::isfinite(Rt(i, j))) {
                return damp::nullopt;
            }
        }
    }
    const auto Xt = mat::lu_solve(Ut, Rt);
    if (!Xt) {
        return damp::nullopt;
    }
    return Xt.value().transpose();
}

/**
 * @brief ‖AᵀX + XA − XGX + Q‖_F / (‖AᵀX‖_F + ‖XA‖_F + ‖XGX‖_F + ‖Q‖_F).
 */
template<size_t NX, typename T>
[[nodiscard]] constexpr T care_relative_residual(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NX, T>& G,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NX, NX, T>& X
) {
    const Matrix<NX, NX, T> AtX = A.transpose() * X;
    const Matrix<NX, NX, T> XA = X * A;
    const Matrix<NX, NX, T> XGX = X * G * X;
    const Matrix<NX, NX, T> Res = AtX + XA - XGX + Q;
    const T                 den = AtX.norm() + XA.norm() + XGX.norm() + Q.norm();
    if (!(den > T{0}) || !damp::isfinite(den)) {
        return Res.norm();
    }
    return Res.norm() / den;
}

/**
 * @brief Map a balanced CARE solution back to the original coordinates.
 *
 * X = D₁ X_b D₁ with D₁ = diag(2^s), then symmetrized.
 */
template<size_t NX, typename T>
[[nodiscard]] constexpr Matrix<NX, NX, T> care_unscale(Matrix<NX, NX, T> X, const damp::array<int, NX>& s) {
    for (size_t i = 0; i < NX; ++i) {
        for (size_t j = 0; j < NX; ++j) {
            X(i, j) = scale_by_pow2(X(i, j), s[i] + s[j]);
        }
    }
    return (X + X.transpose()) * static_cast<T>(0.5);
}

/**
 * @brief Relative-residual certificate for a CARE solution.
 *
 * About 100 √ε (≈ 1.5×10⁻⁶ in double). This is a residual of the Riccati
 * equation, so it scales with √ε rather than with default_tol (an absolute
 * pivot floor). A larger residual is a failed solve.
 */
template<typename T>
[[nodiscard]] constexpr T care_accept_tol() {
    return T{100} * damp::sqrt(eps<T>());
}

/**
 * @brief Up to two Newton steps while the residual is moderate.
 *
 * The step solves (A − GX)ᵀ dX + dX (A − GX) = −R(X) via lyap. It runs only
 * when the subspace checks have already passed and the residual is below 10⁻²,
 * so it polishes a stabilizing solution and does not iterate a corrupted one.
 * A step that fails to decrease the residual is discarded.
 */
template<size_t NX, typename T>
[[nodiscard]] constexpr Matrix<NX, NX, T> care_newton_polish(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NX, T>& G,
    const Matrix<NX, NX, T>& Q,
    Matrix<NX, NX, T>        X
) {
    const T accept = care_accept_tol<T>();
    const T polish_below = static_cast<T>(1e-2);
    for (int step = 0; step < 2; ++step) {
        const T rel = care_relative_residual(A, G, Q, X);
        if (!damp::isfinite(rel) || !(rel > accept) || !(rel < polish_below)) {
            break;
        }
        const Matrix<NX, NX, T> Acl = A - (G * X);
        const Matrix<NX, NX, T> Res = (A.transpose() * X) + (X * A) - (X * G * X) + Q;
        const auto              dX = lyap(Acl.transpose(), Res);
        if (!dX) {
            break;
        }
        const Matrix<NX, NX, T> stepped = X + dX.value();
        const Matrix<NX, NX, T> trial = (stepped + stepped.transpose()) * static_cast<T>(0.5);
        const T                 trial_rel = care_relative_residual(A, G, Q, trial);
        if (!damp::isfinite(trial_rel) || !(trial_rel < rel)) {
            break;
        }
        X = trial;
    }
    return X;
}

/**
 * @brief True when every eigenvalue of A − GX has real part within a few ulps of ≤ 0.
 *
 * The anti-stabilizing Riccati solution also has a tiny residual; its closed-loop
 * poles are the negatives of the stabilizing ones and fail this test.
 */
template<size_t NX, typename T>
[[nodiscard]] constexpr bool care_closed_loop_hurwitz(const Matrix<NX, NX, T>& Acl) {
    const auto ev = mat::compute_eigenvalues(Acl);
    if (!ev.converged) {
        return false;
    }
    const T scale = damp::max(T{1}, Acl.norm());
    //! √ε · ‖A‖: roundoff on a near-imaginary pole stays inside; an
    //! anti-stabilizing pole (order-1 positive real part) does not. (10⁶ ε
    //! was 0.12‖A‖ in float.)
    const T margin = damp::sqrt(eps<T>()) * scale;
    for (size_t i = 0; i < NX; ++i) {
        if (ev.values[i].real() > margin) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Solve CARE via an ordered real Schur method (Laub), symplectically balanced.
 *
 * Solves AᵀX + XA − (XB + N)R⁻¹(BᵀX + Nᵀ) + Q = 0 from the stable invariant
 * subspace of the Hamiltonian
 *
 *     H = ⎡  A   −G  ⎤,   G = B R⁻¹ Bᵀ   (cross-term N folded into A, Q first)
 *         ⎣ −Q   −Aᵀ ⎦
 *
 * H is symplectically balanced, reduced to real Schur form Zᵀ H Z = T
 * (Hessenberg + Francis double-shift QR), and reordered so the NX stable
 * eigenvalues (Re λ < 0) lead. The leading Schur vectors [U₁₁; U₂₁] span the
 * stabilizing subspace and X = U₂₁ U₁₁⁻¹, mapped back from the balanced scaling.
 *
 * Returns damp::nullopt unless the reorder isolates exactly NX stable
 * eigenvalues, U₁₁ᵀ U₂₁ passes Laub's symmetry test, U₁₁ is nonsingular, the
 * closed loop A − GX is Hurwitz, and the relative residual is at most about
 * 100 √ε. Q and R definiteness are checked by care().
 *
 * @see Laub, "A Schur method for solving algebraic Riccati equations," IEEE TAC
 *      1979, https://doi.org/10.1109/TAC.1979.1102178
 * @see Benner, "Symplectic Balancing of Hamiltonian Matrices," SIAM J. Sci.
 *      Comput. 22(5), 2001, https://doi.org/10.1137/S1064827500367993
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

    //! Assemble the 2NX×2NX Hamiltonian H = [[A_eff, −G], [−Q_eff, −A_effᵀ]]
    //! and balance it before the Schur reduction.
    constexpr size_t     M = 2 * NX;
    Matrix<M, M, T>      H = Matrix<M, M, T>::zeros();
    damp::array<int, NX> scale{};
    H.template block<NX, NX>(0, 0) = A_eff;
    H.template block<NX, NX>(0, NX) = G * T{-1};
    H.template block<NX, NX>(NX, 0) = Q_eff * T{-1};
    H.template block<NX, NX>(NX, NX) = A_eff.transpose() * T{-1};
    balance_hamiltonian(H, scale);

    //! Real Schur form Zᵀ H Z = T, accumulating Schur vectors in Z.
    Matrix<M, M, T>   T_schur = H;
    Matrix<M, M, T>   Z;
    damp::array<T, M> wr{};
    damp::array<T, M> wi{};
    mat::detail::hessenberg_reduce(T_schur, Z);
    if (!mat::detail::francis_qr(T_schur, Z, wr, wi)) {
        return damp::nullopt;
    }

    //! Reorder the NX stable eigenvalues (Re λ < 0) into the leading block. A
    //! short count, or a leading block whose real part is positive, means the
    //! stable spectrum was not separated.
    if (reorder_schur(T_schur, Z, [](T re) { return re < T{0}; }) != NX) {
        return damp::nullopt;
    }
    if (!schur_prefix_stable(T_schur, NX)) {
        return damp::nullopt;
    }

    const Matrix<NX, NX, T> U11 = Z.template block<NX, NX>(0, 0);
    const Matrix<NX, NX, T> U21 = Z.template block<NX, NX>(NX, 0);
    if (!care_basis_symmetric(U11, U21)) {
        return damp::nullopt;
    }
    const auto Xb = care_basis_solve(U11, U21);
    if (!Xb) {
        return damp::nullopt;
    }

    //! Undo the symplectic scaling, polish a moderate residual, then certify.
    Matrix<NX, NX, T> X = care_unscale(Xb.value(), scale);
    X = care_newton_polish(A_eff, G, Q_eff, X);
    const T rel = care_relative_residual(A_eff, G, Q_eff, X);
    if (!damp::isfinite(rel) || rel > care_accept_tol<T>()) {
        return damp::nullopt;
    }
    if (!care_closed_loop_hurwitz(A_eff - (G * X))) {
        return damp::nullopt;
    }
    return X;
}

} // namespace detail

namespace detail {

/**
 * @brief Positive semidefinite test: Cholesky of M + δI with δ = N·64·eps·max|mᵢᵢ|.
 *
 * The shift is relative to M's own scale, so a uniformly tiny M is judged the
 * same as a unit one; the floor at the smallest normal keeps M = 0 PSD.
 */
template<size_t N, typename T>
[[nodiscard]] constexpr bool is_psd(const Matrix<N, N, T>& M) {
    T diag_max = T{0};
    for (size_t i = 0; i < N; ++i) {
        diag_max = damp::max(diag_max, damp::abs(M(i, i)));
    }
    //! Above cholesky's pivot floor (16·(i+1)·eps of the pivot's terms, i < N), so a shifted zero pivot clears it.
    const T shift = (static_cast<T>(N) * T{256} * eps<T>() * diag_max) + std::numeric_limits<T>::min();
    return mat::cholesky(M + Matrix<N, N, T>::identity() * shift).has_value();
}

/**
 * @brief True when A − BK is Schur stable for K = (R + BᵀXB)⁻¹(BᵀXA + Nᵀ).
 *
 * The certificate that X is the stabilizing DARE solution. SDA and RDE both
 * converge to a non-stabilizing fixed point when (A, Q) has an unstable
 * unobservable mode (e.g. Q = 0 with unstable A gives X = 0); the residual
 * alone cannot tell them apart.
 */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr bool dare_closed_loop_schur(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& X,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N
) {
    const auto K = mat::solve(Matrix<NU, NU, T>(R + B.t() * X * B), Matrix<NU, NX, T>(B.t() * X * A + N.t()));
    if (!K) {
        return false;
    }
    //! compute_eigenvalues balances the closed loop, so the poles are computed to a relative
    //! accuracy independent of the states' units; the PBH test below runs on the pair
    //! balanced the same way.
    const auto ev = mat::compute_eigenvalues(Matrix<NX, NX, T>(A - B * K.value()));
    if (!ev.converged) {
        return false;
    }
    //! √eps slack keeps a lightly damped pole near |z| = 1 from failing on roundoff.
    //! An uncontrollable mode is invariant under feedback, so one on the unit circle
    //! shows up here; inside the slack band the PBH test tells it from a slow
    //! controllable pole.
    const T    slack = damp::sqrt(eps<T>());
    const auto pair = detail::balanced_pair(A, B);
    for (size_t i = 0; i < NX; ++i) {
        const T mag = ev.values[i].abs();
        if (!(mag < T{1} + slack)) {
            return false;
        }
        if (mag >= T{1} - slack && detail::uncontrollable_at(pair, ev.values[i])) {
            return false;
        }
    }
    return true;
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
 * - Q symmetric positive semidefinite (Q − N R⁻¹ Nᵀ when R > 0)
 * - R symmetric (positive definite for SDA, positive semidefinite for RDE)
 * - (A, B) stabilizable
 *
 * The result is certified: A − BK must be Schur stable, so a non-stabilizing
 * fixed point (e.g. Q = 0 with unstable A) returns damp::nullopt. A Schur-stable
 * A − BK is itself the proof that (A, B) is stabilizable; only a closed-loop pole
 * on the unit circle (within √eps) gets the PBH rank test, so an uncontrollable
 * marginal mode is still rejected without an eigen/rank pass per open-loop pole.
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
    //! R and Q must be symmetric
    if (!mat::is_symmetric_or_hermitian(R) || !mat::is_symmetric_or_hermitian(Q)) {
        return damp::nullopt;
    }
    if (!detail::is_psd(Q)) {
        return damp::nullopt;
    }

    const bool r_is_pd = mat::cholesky(R).has_value();
    if (!r_is_pd && !detail::is_psd(R)) {
        return damp::nullopt;
    }

    //! With a cross term the cost is PSD only if Q − N R⁻¹ Nᵀ is (Schur
    //! complement of the joint weight [[Q, N], [Nᵀ, R]]).
    if (r_is_pd) {
        const auto Rinv_Nt = mat::cholesky_solve(R, N.transpose());
        if (!Rinv_Nt || !detail::is_psd(Matrix<NX, NX, T>(Q - N * Rinv_Nt.value()))) {
            return damp::nullopt;
        }
    }

    damp::optional<Matrix<NX, NX, T>> X;
    switch (method) {
        case DareMethod::SDA:
            if (r_is_pd) {
                X = detail::dare_sda(A, B, Q, R, N);
            }
            break;
        case DareMethod::RDE:
            X = detail::dare_rde(A, B, Q, R, N);
            break;
        case DareMethod::Auto:
        default:
            X = r_is_pd ? detail::dare_sda(A, B, Q, R, N) : detail::dare_rde(A, B, Q, R, N);
            break;
    }

    if (!X || !detail::dare_closed_loop_schur(A, B, X.value(), R, N)) {
        return damp::nullopt;
    }
    return X;
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
 * is_stabilizable() used by dare()). care() returns damp::nullopt when the
 * Schur reduction does not converge, the stable subspace cannot be isolated,
 * or the computed X fails the residual or closed-loop certificate. A matrix
 * that does not solve the equation is not returned.
 *
 * @note Compare with MATLAB®'s icare(A, B, Q, R) / care(A, B, Q, R).
 *
 * @see care_schur() — ordered Schur (Laub) with symplectic balancing (Benner)
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

    //! Q positive semidefinite — same scale-free test as dare().
    if (!detail::is_psd(Q)) {
        return damp::nullopt;
    }

    //! R must be positive definite (G = B R⁻¹ Bᵀ).
    if (!mat::cholesky(R)) {
        return damp::nullopt;
    }

    return detail::care_schur(A, B, Q, R, N);
}

} // namespace damp