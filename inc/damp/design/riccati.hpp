// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file riccati.hpp
 * @brief Algebraic Riccati equation solvers (DARE/CARE) and the LQR built on them:
 *        discrete_lqr, continuous_lqr and the sampled-data cost (discretize_lqr_cost)
 *
 * The LQ designs on a plant (design::lqr, lqi, lqg, lqgi, oflqg, each in its runtime's header in
 * controllers/) solve through these.
 */


#include "damp/backend.hpp"
#include "damp/design/lyapunov.hpp"
#include "damp/design/stability.hpp"
#include "damp/design/status.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/decomposition.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/functions.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace detail {

using mat::detail::gebal_log2_scale;
using mat::detail::scale_by_pow2;

/** @brief A Riccati solution and why it failed: X is meaningful only when status is ok. */
template<size_t NX, typename T>
struct [[nodiscard]] RiccatiSolution {
    Matrix<NX, NX, T>    X{};                                         ///< Stabilizing solution (status ok)
    design::DesignStatus status{design::DesignStatus::not_converged}; ///< ok, or the cause class of the failure
};

/** @brief The residual of a Riccati equation at X, the gain and closed loop it implies, and its relative size. */
template<size_t NX, size_t NU, typename T>
struct RiccatiResidual {
    Matrix<NX, NX, T> res{};        ///< R(X), symmetrized
    Matrix<NU, NX, T> K{};          ///< the gain at X: R⁻¹(BᵀX + Nᵀ) or (R + BᵀXB)⁻¹(BᵀXA + Nᵀ)
    Matrix<NX, NX, T> Acl{};        ///< A − BK
    T                 rel{T{1}};    ///< ‖R(X)‖_F / Σ ‖term‖_F
    T                 floor{T{1}};  ///< max(1, ρ): how far rounding X alone can lift rel above n·ε
    bool              valid{false}; ///< the gain solve succeeded and rel is finite
};

/**
 * @brief Riccati residual in the factored form, from the original data.
 *
 * CARE: R(X) = AᵀX + XA − (XB + N) K + Q,         K = R⁻¹(BᵀX + Nᵀ).
 * DARE: R(X) = AᵀXA − X − (AᵀXB + N) K + Q,      K = (R + BᵀXB)⁻¹(BᵀXA + Nᵀ).
 *
 * The quadratic term is formed through the gain, never as X·(B R⁻¹ Bᵀ)·X: with a small R
 * (cheap control) G = B R⁻¹ Bᵀ is huge while XB is tiny, and X G X then cancels to an
 * O(1) result from O(‖X‖²‖G‖) terms, losing every digit the residual is meant to measure.
 * rel is the residual over the sum of the terms' norms, so it is unit-free.
 *
 * floor is the reach of rounding: perturbing X by ε|X| moves R(X) by the closed-loop
 * Lyapunov (Stein) operator, at most ε·‖|A_cl|ᵀ|X| + |X||A_cl|‖ (‖|A_cl|ᵀ|X||A_cl| + |X|‖),
 * so ρ, that over Σ‖terms‖, is how many n·ε of residual even the correctly rounded solution
 * may carry. It is about 1 unless ‖A_cl‖ dwarfs the terms, as with cheap control, whose fast
 * closed-loop poles sit near −1/√R.
 */
template<bool Discrete, size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr RiccatiResidual<NX, NU, T> riccati_residual(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N,
    const Matrix<NX, NX, T>& X
) {
    RiccatiResidual<NX, NU, T>        out{};
    damp::array<Matrix<NX, NX, T>, 4> terms{};
    if constexpr (Discrete) {
        const Matrix<NX, NX, T> XA = X * A;
        const Matrix<NX, NU, T> XB = X * B;
        const auto              K = mat::solve(Matrix<NU, NU, T>(R + (B.t() * XB)), Matrix<NU, NX, T>((B.t() * XA) + N.t()));
        if (!K) {
            return out;
        }
        out.K = K.value();
        terms = {Matrix<NX, NX, T>(A.t() * XA), X * T{-1}, Matrix<NX, NX, T>(((A.t() * XB) + N) * out.K) * T{-1}, Q};
    } else {
        const Matrix<NX, NU, T> L = (X * B) + N;
        const auto              K = mat::solve(R, Matrix<NU, NX, T>(L.t()));
        if (!K) {
            return out;
        }
        out.K = K.value();
        terms = {Matrix<NX, NX, T>(A.t() * X), Matrix<NX, NX, T>(X * A), Matrix<NX, NX, T>(L * out.K) * T{-1}, Q};
    }
    out.Acl = A - (B * out.K);
    Matrix<NX, NX, T> res = Matrix<NX, NX, T>::zeros();
    T                 den{0};
    for (const auto& t : terms) {
        res = res + t;
        den += t.norm();
    }
    const Matrix<NX, NX, T> absX = mat::transform(X, [](T v) { return damp::abs(v); });
    const Matrix<NX, NX, T> absAcl = mat::transform(out.Acl, [](T v) { return damp::abs(v); });
    const Matrix<NX, NX, T> reach = Discrete ? Matrix<NX, NX, T>((absAcl.t() * absX * absAcl) + absX)
                                             : Matrix<NX, NX, T>((absAcl.t() * absX) + (absX * absAcl));
    out.res = (res + res.t()) * static_cast<T>(0.5);
    out.rel = (den > T{0}) ? (res.norm() / den) : res.norm();
    out.floor = (den > T{0}) ? damp::max(T{1}, reach.norm() / den) : T{1};
    out.valid = damp::isfinite(out.rel) && damp::isfinite(out.floor);
    return out;
}

/**
 * @brief Relative residual a Riccati solution is polished to: 4·n·ε.
 *
 * Forming each term of the residual is a length-n inner product, so the computed relative
 * residual of the exact solution is itself of order n·ε; below that a Newton step cannot be
 * confirmed to help. 4·n·ε stops on the step that reaches the rounding floor, 8× inside the
 * 32·n·ε the oracle suite (tests/test_oracle_*.cpp) holds the designs to.
 */
template<size_t NX, typename T>
[[nodiscard]] constexpr T riccati_polish_tol() {
    return T{4} * static_cast<T>(NX) * eps<T>();
}

/**
 * @brief The residual certificate: a relative residual above 32·n·ε·floor is a failed solve.
 *
 * After the Newton polish the residual of a solvable problem sits at its rounding floor
 * (RiccatiResidual::floor n·ε) whatever the conditioning of X: the residual is a backward
 * error. One that stays 8× above the polish target, at that floor, did not converge.
 */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr bool riccati_residual_certified(const RiccatiResidual<NX, NU, T>& r) {
    return r.valid && r.rel <= T{8} * riccati_polish_tol<NX, T>() * r.floor;
}

/** @brief A polished Riccati solution and its residual. */
template<size_t NX, size_t NU, typename T>
struct PolishedRiccati {
    Matrix<NX, NX, T>          X{}; ///< the solution after the polish
    RiccatiResidual<NX, NU, T> r{}; ///< its residual, gain and closed loop
};

/**
 * @brief Newton (Kleinman / Hewer) steps on a Riccati solution while its residual is above n·ε.
 *
 * One step linearizes the equation at X about the closed loop A_cl = A − BK(X) and solves
 * for the correction, A_clᵀ dX + dX A_cl = −R(X) (lyap) or A_clᵀ dX A_cl − dX = −R(X)
 * (dlyap), so the residual of X + dX is quadratic in dX. From the Schur or SDA answer, which
 * already lies in the region where Newton converges quadratically, it reaches
 * riccati_polish_tol() in one or two steps; there (residual below 10⁻²) a step that does not
 * decrease the residual is discarded and ends the polish.
 *
 * From a start further out whose gain is stabilizing (an RDE iterate stopped at its cap), the
 * steps are Hewer's (Kleinman's) iteration: every iterate stays stabilizing and X decreases
 * monotonically to the stabilizing solution, while the residual need not, so those steps are
 * taken whatever the residual does. Any other start (residual above 10⁻², gain not
 * stabilizing) is returned as is. At most 30 steps are taken.
 *
 * @see "On an iterative technique for Riccati equation computations" (Kleinman, 1968)
 * @see "An iterative technique for the computation of the steady state gains for the discrete
 *      optimal regulator" (Hewer, 1971)
 */
template<bool Discrete, size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr PolishedRiccati<NX, NU, T> riccati_newton_polish(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N,
    const Matrix<NX, NX, T>& X0
) {
    PolishedRiccati<NX, NU, T> out{.X = X0, .r = riccati_residual<Discrete>(A, B, Q, R, N, X0)};
    const T                    target = riccati_polish_tol<NX, T>();
    const T                    polish_below = static_cast<T>(1e-2);
    constexpr int              kMaxSteps = 30;
    const auto                 stabilizing = [](const Matrix<NX, NX, T>& Acl) {
        const auto ev = mat::compute_eigenvalues(Acl);
        for (size_t i = 0; i < NX && ev.success; ++i) {
            if (Discrete ? !(ev.values[i].abs() < T{1}) : !(ev.values[i].real() < T{0})) {
                return false;
            }
        }
        return ev.success;
    };
    const bool hewer = out.r.valid && !(out.r.rel < polish_below) && stabilizing(out.r.Acl);
    for (int step = 0; step < kMaxSteps; ++step) {
        const bool far = !(out.r.rel < polish_below);
        if (!out.r.valid || !(out.r.rel > target) || (far && !hewer)) {
            break;
        }
        const Matrix<NX, NX, T>           Acl_t = out.r.Acl.t();
        damp::optional<Matrix<NX, NX, T>> dX{};
        if constexpr (Discrete) {
            dX = dlyap(Acl_t, out.r.res);
        } else {
            dX = lyap(Acl_t, out.r.res);
        }
        if (!dX) {
            break;
        }
        const Matrix<NX, NX, T> stepped = out.X + dX.value();
        const Matrix<NX, NX, T> trial = (stepped + stepped.t()) * static_cast<T>(0.5);
        const auto              next = riccati_residual<Discrete>(A, B, Q, R, N, trial);
        if (!next.valid || (!far && !(next.rel < out.r.rel))) {
            break;
        }
        out.X = trial;
        out.r = next;
    }
    return out;
}

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

    const Matrix            A_eff = A - B * R_inv_Nt;
    const Matrix<NX, NX, T> Q_fold = Q - N * R_inv_Nt;
    const Matrix<NX, NX, T> Q_eff = (Q_fold + Q_fold.t()) * T{0.5}; // the product is symmetric only to roundoff

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
    Gk = (Gk + Gk.t()) * T{0.5};
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
 * Linear convergence, by ρ(A − BK)² per step. Handles R ≥ 0 (only requires R + BᵀXB
 * invertible at each step, not R itself). Useful when R is singular (e.g., cheap-control
 * problems or minimum-energy estimation). With a slow closed loop the 500-step cap comes
 * before the step tolerance; the last iterate is then returned, and dare() finishes it with
 * Hewer's iteration (riccati_newton_polish()) once its gain is stabilizing. Only a singular
 * R + BᵀXB or a non-finite iterate returns damp::nullopt.
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
        Q_eff = (Q_eff + Q_eff.t()) * T{0.5}; // the product is symmetric only to roundoff
    }

    //! Riccati Difference Equation iteration
    Matrix<NX, NX, T> X = Q_eff;
    const T           tol = default_tol<T>(); // relative to ‖X‖ (see dare_sda)
    const int         max_iter = 500;
    // Type-safe ceiling: float cannot represent 1e150 (narrowing under -Wnarrowing).
    const T guard = damp::numeric_limits<T>::max() / T{100};

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
        if (mat::any_of(X_next, [guard](T v) { return !damp::isfinite(v) || damp::abs(v) > guard; })) {
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
            break;
        }
    }
    return (X + X.t()) * static_cast<T>(0.5); // converged, or the iterate at the cap (dare() certifies it)
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
    if constexpr (S > M) {
        return false; //! the pair does not fit in M; reorder_schur's run-time dispatch never asks
    } else if constexpr (P == 1 && Q == 1) {
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
    } else {
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
        if (!x_opt || !mat::all_of(*x_opt, [](T v) { return damp::isfinite(v); })) {
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

        //! Bai–Demmel: try the swap on the window before touching Tm and Z.
        const Matrix<S, S, T> W = Tm.template block<S, S>(j, j).to_matrix();
        T                     wmax = T{0};
        for (size_t r = 0; r < S; ++r) {
            for (size_t c = 0; c < S; ++c) {
                wmax = damp::max(wmax, damp::abs(W(r, c)));
            }
        }
        const T         tau = damp::max(T{10} * static_cast<T>(S) * eps<T>() * wmax, damp::numeric_limits<T>::min());
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
 * 2×2 block is a genuine complex pair, then bubbles each block whose eigenvalue satisfies
 * the predicate to the top-left via adjacent-block swaps, keeping the Schur vectors @p Z
 * orthogonal. The predicate takes the block's eigenvalue as a damp::complex (the one with
 * Im λ ≥ 0 for a pair). For the CARE Hamiltonian it selects the stable spectrum (Re λ < 0),
 * collecting the stabilizing invariant subspace into the leading columns; for a discrete
 * split it selects |λ| < 1.
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
        const T      im = (b == 2) ? (damp::sqrt(damp::abs(Tm(i, i + 1))) * damp::sqrt(damp::abs(Tm(i + 1, i)))) : T{0};
        if (in_front(damp::complex<T>{re, im})) {
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
        if (half > static_cast<T>(kCap)) {
            half = static_cast<T>(kCap);
        } else if (half < -static_cast<T>(kCap)) {
            half = -static_cast<T>(kCap);
        }
        const int si = static_cast<int>(damp::nearbyint(half)); // |half| ≤ kCap: the cast is in range
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
 * Row scaling leaves X unchanged and makes partial pivoting choose pivots by
 * their size relative to their row (scaled partial pivoting), so a basis with
 * rows of very different magnitude is not pivoted on its largest row alone; a
 * genuinely singular U₁₁ still has a zero row.
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

/** @brief The cost with its cross term folded in: A − B R⁻¹ Nᵀ and Q − N R⁻¹ Nᵀ (A and Q when R is singular, where N = 0). */
template<size_t NX, typename T>
struct CostPair {
    Matrix<NX, NX, T> A{}; ///< A − B R⁻¹ Nᵀ
    Matrix<NX, NX, T> Q{}; ///< Q − N R⁻¹ Nᵀ
};

/**
 * @copybrief CostPair
 *
 * Both differences cancel to roundoff for an output cost (Q = CᵀC, N = CᵀD, R = DᵀD), and the
 * PBH test balances the pair and scales its rows, which would lift that roundoff to a weight
 * or a coupling: a boundary mode of the zero dynamics then reads as seen. So an entry within
 * n·256·ε of the scale of its terms is set to zero: √(sᵢsⱼ) with sᵢ = max(|qᵢᵢ|, |(N R⁻¹ Nᵀ)ᵢᵢ|)
 * for Q (the tolerance riccati_weights_valid() allows), max(|aᵢⱼ|, |(B R⁻¹ Nᵀ)ᵢⱼ|) for A. The
 * scales follow the states' units, and with N = 0 nothing changes.
 */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr CostPair<NX, T> cost_pair(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N
) {
    CostPair<NX, T> out{.A = A, .Q = Q};
    if (const auto Rinv_Nt = mat::cholesky_solve(R, N.transpose())) {
        const Matrix<NX, NX, T> BRN = B * Rinv_Nt.value();
        const Matrix<NX, NX, T> NRN = N * Rinv_Nt.value();
        out.A = A - BRN;
        out.Q = Q - NRN;
        out.Q = (out.Q + out.Q.t()) * T{0.5}; // the product is symmetric only to roundoff
        const T tol = static_cast<T>(NX) * T{256} * eps<T>();
        for (size_t i = 0; i < NX; ++i) {
            for (size_t j = 0; j < NX; ++j) {
                const T si = damp::max(damp::abs(Q(i, i)), damp::abs(NRN(i, i)));
                const T sj = damp::max(damp::abs(Q(j, j)), damp::abs(NRN(j, j)));
                if (damp::abs(out.Q(i, j)) <= tol * damp::sqrt(si * sj)) {
                    out.Q(i, j) = T{0};
                }
                if (damp::abs(out.A(i, j)) <= tol * damp::max(damp::abs(A(i, j)), damp::abs(BRN(i, j)))) {
                    out.A(i, j) = T{0};
                }
            }
        }
    }
    return out;
}

/**
 * @brief True when A − BK is stable with every pole off the stability boundary by more than
 *        roundoff: Schur (@p Discrete) or Hurwitz.
 *
 * The certificate that X is the stabilizing Riccati solution. SDA and RDE can converge to a
 * non-stabilizing fixed point (Q = 0 with unstable A gives X = 0), and the anti-stabilizing
 * CARE solution has a tiny residual too; the residual alone cannot tell them apart. Two bands
 * on the distance inside the boundary (1 − |λ| or −Re λ), with the norm of the balanced closed
 * loop ‖D⁻¹(A − BK)D‖_F (D the balancing compute_eigenvalues applies) as the scale where
 * one is needed, so the verdict does not depend on the states' units or the time unit:
 * - within n·ε·max(1, ‖D⁻¹(A − BK)D‖_F) (discrete) or n·ε·‖D⁻¹(A − BK)D‖_F (continuous),
 *   the pole is on the boundary to within the roundoff of computing it: refused. No
 *   stabilizing solution has one.
 * - within √ε (discrete: |λ| has no units) or 4·√ε·‖D⁻¹(A − BK)D‖_F (continuous), the pole
 *   must be controllable, PBH on (A, B), and seen by the cost, the same test on
 *   (A_effᵀ, Q_eff) (cost_pair()). An uncontrollable mode is invariant under feedback, so one
 *   on the boundary stays there whatever X is. The maximal solution of a problem with a
 *   controllable boundary mode that Q does not see (a position loop that weights only
 *   velocity) leaves that mode's pole on the boundary too. Either is computed only to
 *   roundoff, amplified by the pole's condition number, or to O(√ε) when defective. The CARE's
 *   boundary pole is a defective double root of the Hamiltonian and lands up to about 2× outside
 *   √ε·‖D⁻¹(A − BK)D‖_F, hence the wider continuous band; the rank tests tell both cases from a
 *   lightly damped pole, so a wider band refuses only poles that fail them.
 * - with a pole within 16× that band, a boundary mode of A_eff that Q_eff does not see (the
 *   same test on the open loop's eigenvalues) is refused: no stabilizing solution exists, and
 *   the solve can leave that pole O(√ε) inside, up to 2× outside the band.
 *
 * compute_eigenvalues balances the closed loop, so the poles are computed to a relative
 * accuracy independent of the states' units; the PBH tests run on pairs balanced the same way.
 *
 * @see "Controllability and observability conditions of linear autonomous systems" (Hautus, 1969)
 */
template<bool Discrete, size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr bool closed_loop_certified(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N,
    const Matrix<NX, NX, T>& Acl
) {
    const auto ev = mat::compute_eigenvalues(Acl);
    if (!ev.success) {
        return false;
    }
    const T            bal_norm = mat::detail::balanced(Acl, mat::detail::balance_exponents(Acl)).norm(); // the scale of the poles' roundoff
    const T            on_boundary = static_cast<T>(NX) * eps<T>() * (Discrete ? damp::max(T{1}, bal_norm) : bal_norm);
    const T            band = damp::sqrt(eps<T>()) * (Discrete ? T{1} : T{4} * bal_norm);
    damp::array<T, NX> margin{};
    bool               in_band = false;
    bool               near = false;
    for (size_t i = 0; i < NX; ++i) {
        margin[i] = Discrete ? (T{1} - ev.values[i].abs()) : -ev.values[i].real(); // distance inside the boundary
        if (!damp::isfinite(margin[i]) || !(margin[i] > on_boundary)) {
            return false;
        }
        in_band = in_band || (margin[i] <= band);
        near = near || (margin[i] <= T{16} * band);
    }
    if (!near) {
        return true;
    }
    //! A boundary mode the cost does not see leaves no stabilizing solution, but the solve puts
    //! its pole only O(√ε) inside, at times outside the band: the open loop's is tested too.
    const auto              cost = cost_pair(A, B, Q, R, N);
    const Matrix<NX, NX, T> cost_At = cost.A.transpose();
    if (stability::detail::uncontrollable_unstable_mode<Discrete>(cost_At, cost.Q, true)) {
        return false;
    }
    if (!in_band) {
        return true;
    }
    const auto pair = stability::detail::balanced_pair(A, B);
    const auto seen = stability::detail::balanced_pair(cost_At, cost.Q);
    for (size_t i = 0; i < NX; ++i) {
        if (margin[i] <= band && (stability::detail::uncontrollable_at(pair, ev.values[i]) || stability::detail::uncontrollable_at(seen, ev.values[i]))) {
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
 * X is then polished by Newton steps (riccati_newton_polish()) until its relative
 * residual reaches the n·ε rounding floor: the Schur solution of an ill-conditioned
 * problem (cond U₁₁ large) is accurate only to κ·ε and needs it.
 *
 * Fails with not_converged unless the reorder isolates exactly NX stable eigenvalues,
 * U₁₁ᵀ U₂₁ passes Laub's symmetry test, U₁₁ is nonsingular and the relative residual is at
 * most 32·n·ε times its rounding floor (riccati_residual_certified()), and with
 * not_stabilizable unless the closed loop A − BK is Hurwitz with no pole on the imaginary
 * axis (closed_loop_certified()). solve_care() checks the weights and classifies a failure
 * by the problem's structure (riccati_failure()).
 *
 * @see Laub, "A Schur method for solving algebraic Riccati equations," IEEE TAC
 *      1979, https://doi.org/10.1109/TAC.1979.1102178
 * @see Benner, "Symplectic Balancing of Hamiltonian Matrices," SIAM J. Sci.
 *      Comput. 22(5), 2001, https://doi.org/10.1137/S1064827500367993
 * @see Golub & Van Loan, "Matrix Computations" §7.6 (ordered Schur form)
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr RiccatiSolution<NX, T> care_schur(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N = Matrix<NX, NU, T>{}
) {
    //! G = B R⁻¹ Bᵀ and cross-term reduction A_eff = A − B R⁻¹ Nᵀ,
    //! Q_eff = Q − N R⁻¹ Nᵀ — solved against R rather than forming R⁻¹.
    using design::DesignStatus;
    RiccatiSolution<NX, T> out{};
    const auto             Rinv_Bt_opt = mat::lu_solve(R, B.transpose());
    const auto             Rinv_Nt_opt = mat::lu_solve(R, N.transpose());
    if (!Rinv_Bt_opt || !Rinv_Nt_opt) {
        out.status = DesignStatus::bad_weights; // R singular
        return out;
    }
    //! G and Q_eff are symmetric only to roundoff as products; symmetrized, H is Hamiltonian
    //! to the last bit and Laub's symmetry test reads only the solve's error.
    const Matrix<NX, NX, T> G_prod = B * Rinv_Bt_opt.value();
    const Matrix<NX, NX, T> Q_fold = Q - N * Rinv_Nt_opt.value();
    const Matrix<NX, NX, T> G = (G_prod + G_prod.t()) * T{0.5};
    const Matrix<NX, NX, T> A_eff = A - B * Rinv_Nt_opt.value();
    const Matrix<NX, NX, T> Q_eff = (Q_fold + Q_fold.t()) * T{0.5};

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
        return out; // not_converged
    }

    //! Reorder the NX stable eigenvalues (Re λ < 0) into the leading block. A
    //! short count, or a leading block whose real part is positive, means the
    //! stable spectrum was not separated.
    if (reorder_schur(T_schur, Z, [](const damp::complex<T>& lambda) { return lambda.real() < T{0}; }) != NX || !schur_prefix_stable(T_schur, NX)) {
        return out; // not_converged
    }

    const Matrix<NX, NX, T> U11 = Z.template block<NX, NX>(0, 0);
    const Matrix<NX, NX, T> U21 = Z.template block<NX, NX>(NX, 0);
    if (!care_basis_symmetric(U11, U21)) {
        return out; // not_converged
    }
    //! A singular U₁₁ means the stable subspace is not a graph over x (an unstabilizable
    //! (A, B), which solve_care classifies).
    const auto Xb = care_basis_solve(U11, U21);
    if (!Xb) {
        return out; // not_converged
    }

    //! Undo the symplectic scaling, polish to the rounding floor, then certify: the residual
    //! (on the original data, in the factored form) and the closed loop A − BK.
    const auto sol = riccati_newton_polish<false>(A, B, Q, R, N, care_unscale(Xb.value(), scale));
    if (!riccati_residual_certified(sol.r)) {
        return out; // not_converged
    }
    if (!closed_loop_certified<false>(A, B, Q, R, N, sol.r.Acl)) {
        out.status = DesignStatus::not_stabilizable;
        return out;
    }
    out.X = sol.X;
    out.status = DesignStatus::ok;
    return out;
}

} // namespace detail

namespace detail {

/** @brief The largest |mᵢᵢ|: the scale of a symmetric PSD matrix. */
template<size_t N, typename T>
[[nodiscard]] constexpr T diag_max(const Matrix<N, N, T>& M) {
    T out = T{0};
    for (size_t i = 0; i < N; ++i) {
        out = damp::max(out, damp::abs(M(i, i)));
    }
    return out;
}

/**
 * @brief Positive semidefinite test: Cholesky of M + δI with δ = N·256·eps·scale, the scale of the
 *        data M was formed from.
 *
 * The shift is relative to the data's scale, so a uniformly tiny M is judged the same as a
 * unit one; the floor at the smallest normal keeps M = 0 PSD. A matrix formed by cancellation
 * (a Schur complement) passes the scale of its terms: its own diagonal is then roundoff.
 */
template<size_t N, typename T>
[[nodiscard]] constexpr bool is_psd(const Matrix<N, N, T>& M, T scale) {
    //! Above cholesky's pivot floor (16·(i+1)·eps of the pivot's terms, i < N), so a shifted zero pivot clears it.
    const T shift = (static_cast<T>(N) * T{256} * eps<T>() * scale) + damp::numeric_limits<T>::min();
    return mat::cholesky(M + Matrix<N, N, T>::identity() * shift).has_value();
}

/** @brief Positive semidefinite test with the shift scaled to M's own diagonal, max|mᵢᵢ|. */
template<size_t N, typename T>
[[nodiscard]] constexpr bool is_psd(const Matrix<N, N, T>& M) {
    return is_psd(M, diag_max(M));
}

/**
 * @brief The weights of a Riccati equation are consistent: Q and R symmetric, Q PSD, and the
 *        joint weight [[Q, N], [Nᵀ, R]] PSD.
 *
 * With R positive definite the joint weight is PSD iff its Schur complement Q − N R⁻¹ Nᵀ is.
 * A singular R (R PSD, only when @p singular_r) is taken with N = 0 alone: the reduction of N
 * needs R⁻¹.
 */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr bool riccati_weights_valid(
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N,
    bool                     singular_r
) {
    if (!mat::is_symmetric_or_hermitian(R) || !mat::is_symmetric_or_hermitian(Q) || !is_psd(Q)) {
        return false;
    }
    if (!mat::cholesky(R).has_value()) {
        return singular_r && is_psd(R) && mat::all_of(N, [](T v) { return v == T{0}; });
    }
    const auto Rinv_Nt = mat::cholesky_solve(R, N.transpose());
    if (!Rinv_Nt.has_value()) {
        return false;
    }
    //! The complement cancels to roundoff for a cost on z = Cx + Du (Q = CᵀC, N = CᵀD, R = DᵀD),
    //! so its shift scales with its terms, not with its own diagonal.
    const Matrix<NX, NX, T> NRN = N * Rinv_Nt.value();
    const Matrix<NX, NX, T> S = Q - NRN;
    return is_psd(Matrix<NX, NX, T>((S + S.t()) * T{0.5}), damp::max(diag_max(Q), diag_max(NRN)));
}

/**
 * @brief Why a Riccati solve with valid weights failed, from the problem's structure.
 *
 * A stabilizing solution exists iff (A, B) is stabilizable and no mode on the stability
 * boundary is unseen by the cost, i.e. (A − B R⁻¹ Nᵀ, Q − N R⁻¹ Nᵀ) has no unobservable
 * eigenvalue on it (with a singular R, N = 0 and the pair is (A, Q)). So: not_stabilizable
 * when (A, B) is not stabilizable, bad_weights when such a mode is unseen (Q is what to
 * change), else @p otherwise: a solver that did not converge, or a certificate that failed.
 *
 * @see "Optimal Control" (Anderson & Moore, 1990), §3.3 and §4.3
 */
template<bool Discrete, size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr design::DesignStatus riccati_failure(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N,
    design::DesignStatus     otherwise
) {
    if (stability::detail::uncontrollable_unstable_mode<Discrete>(A, B, false)) {
        return design::DesignStatus::not_stabilizable;
    }
    const auto cost = cost_pair(A, B, Q, R, N);
    if (stability::detail::uncontrollable_unstable_mode<Discrete>(Matrix<NX, NX, T>(cost.A.transpose()), cost.Q, true)) {
        return design::DesignStatus::bad_weights; // a boundary mode the cost does not see
    }
    return otherwise;
}

} // namespace detail

/** @brief DARE solver selection for dare(). */
enum class DareMethod : uint8_t {
    Auto, ///< SDA when R > 0, else RDE
    SDA,  ///< Structure-preserving doubling: needs R > 0
    RDE,  ///< Riccati difference iteration: R ≥ 0
};

namespace detail {

/**
 * @brief dare() with the reason it failed: bad_weights, not_stabilizable or not_converged.
 *
 * Inconsistent weights (riccati_weights_valid(); R must be positive definite for SDA, and with
 * a singular R, N must be zero) are bad_weights. A solve that fails is classified once, on the
 * failure path, by riccati_failure(): an (A, B) that is not stabilizable is not_stabilizable,
 * a boundary mode the cost does not see is bad_weights, and otherwise a solver that stopped
 * without a certified residual is not_converged and a closed loop that fails its certificate
 * not_stabilizable.
 */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr RiccatiSolution<NX, T> solve_dare(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N,
    DareMethod               method
) {
    using design::DesignStatus;
    RiccatiSolution<NX, T> out{};
    out.status = DesignStatus::bad_weights;
    if (!riccati_weights_valid(Q, R, N, method != DareMethod::SDA)) {
        return out;
    }
    const bool r_is_pd = mat::cholesky(R).has_value();

    //! Bit-pattern finiteness: under -ffast-math a NaN iterate can pass every comparison
    //! downstream (an unstabilizable pair makes SDA blow up), so it is caught here.
    const auto finite = [](const damp::optional<Matrix<NX, NX, T>>& M) {
        return M.has_value() && mat::all_of(*M, [](T v) { return damp::isfinite(v); });
    };
    const bool use_sda = method == DareMethod::SDA || (method == DareMethod::Auto && r_is_pd);
    auto       X = use_sda ? dare_sda(A, B, Q, R, N) : dare_rde(A, B, Q, R, N);
    //! SDA works with G = B R⁻¹ Bᵀ, and I + GₖHₖ has condition about ‖G‖‖X‖: a tiny R (cheap
    //! control) costs SDA that many digits, all of them in float (R = 10⁻⁸ against
    //! 1/ε_float ≈ 10⁷). RDE works with R + BᵀXB instead, which stays the size of BᵀXB,
    //! so Auto falls back to it when SDA fails on a stabilizable pair (on an unstabilizable
    //! one both diverge). The Newton polish below restores whatever digits either loses.
    if (!finite(X) && use_sda && method == DareMethod::Auto && !stability::detail::uncontrollable_unstable_mode<true>(A, B, false)) {
        X = dare_rde(A, B, Q, R, N);
    }
    const auto failed = [&](DesignStatus otherwise) {
        out.status = riccati_failure<true>(A, B, Q, R, N, otherwise);
        return out;
    };
    if (!finite(X)) {
        return failed(DesignStatus::not_converged);
    }
    //! Polish to the rounding floor, then certify: the residual (SDA and RDE stop on the size
    //! of their step, which says nothing about it), and the closed loop A − BK — a fixed point
    //! whose A − BK is not Schur stable is not the stabilizing solution.
    const auto sol = riccati_newton_polish<true>(A, B, Q, R, N, X.value());
    if (!riccati_residual_certified(sol.r)) {
        return failed(DesignStatus::not_converged);
    }
    if (!closed_loop_certified<true>(A, B, Q, R, N, sol.r.Acl)) {
        return failed(DesignStatus::not_stabilizable);
    }
    out.X = sol.X;
    out.status = DesignStatus::ok;
    return out;
}

} // namespace detail

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
 * - R symmetric (positive definite for SDA, positive semidefinite for RDE); a nonzero N
 *   needs R positive definite
 * - (A, B) stabilizable, and no mode on the unit circle unseen by Q
 *
 * SDA and RDE stop on the size of their step, so their answer is Newton-polished
 * (riccati_newton_polish(), Hewer's iteration) to the n·ε rounding floor of the
 * residual; an RDE stopped at its cap with a stabilizing gain is finished the same way.
 * With Auto, an SDA that fails on a stabilizable pair (a tiny R makes I + GₖHₖ singular
 * to working precision) is retried with RDE.
 *
 * The result is certified: the relative residual must be at most 32·n·ε times its
 * rounding floor (riccati_residual_certified()), and A − BK must be Schur stable with no
 * pole on the unit circle to within roundoff (closed_loop_certified()), so a
 * non-stabilizing fixed point (e.g. Q = 0 with unstable A, or the maximal solution when Q
 * does not see a marginal mode) returns damp::nullopt. A Schur-stable A − BK is itself the
 * proof that (A, B) is stabilizable; only a closed-loop pole within √eps of the unit circle
 * gets the PBH rank tests, on (A, B) and on the cost's dual, so an uncontrollable marginal
 * mode, or one Q does not see, is rejected wherever roundoff puts its pole in that band,
 * without an eigen/rank pass per open-loop pole.
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
    const auto sol = detail::solve_dare(A, B, Q, R, N, method);
    if (sol.status != design::DesignStatus::ok) {
        return damp::nullopt;
    }
    return sol.X;
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
 * numerically stabler and cheaper route. Shared by discrete_lqr / lqi gain synthesis.
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

namespace detail {

/**
 * @brief care() with the reason it failed.
 *
 * Q or R not symmetric, Q not PSD, R not positive definite or the joint weight
 * [[Q, N], [Nᵀ, R]] not PSD is bad_weights (riccati_weights_valid(): the weights are what to
 * change). A failed Schur solve is classified as the DARE's is (riccati_failure()), so the
 * same defect gets the same status in either domain.
 */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr RiccatiSolution<NX, T> solve_care(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N
) {
    if (!riccati_weights_valid(Q, R, N, false)) {
        RiccatiSolution<NX, T> out{};
        out.status = design::DesignStatus::bad_weights;
        return out;
    }
    auto sol = care_schur(A, B, Q, R, N);
    if (sol.status != design::DesignStatus::ok) {
        sol.status = riccati_failure<false>(A, B, Q, R, N, sol.status);
    }
    return sol;
}

/**
 * @brief ∫₀^Ts e^{Fᵀτ} S e^{Fτ} dτ (Van Loan), by scaling and squaring on the pair (e^{Fh}, W(h)).
 *
 * The integral of a cost (F = Ā, S = Q̄) or of process noise (F = Aᵀ, S = G Q_c Gᵀ) over one
 * period. A single exponential of @f$ [-F^\top\; S;\, 0\; F]\,T_s @f$ gives
 * @f$ W = M_{22}^\top M_{12} @f$, but a fast stable mode makes @f$ e^{-F^\top T_s} @f$ grow
 * like @f$ e^{|\lambda| T_s} @f$, so @f$ M_{12} @f$ is that large and the product carries an
 * absolute error of @f$ \varepsilon\, e^{|\lambda| T_s} \|S\| @f$ (from |λ|·Ts ≈ 35 in
 * double, 15 in float). The block is taken instead at @f$ h = T_s/2^s @f$ with
 * @f$ \|F\|_F h \le 1/2 @f$, where nothing cancels, and doubled s times:
 * @f$ W(2h) = W(h) + e^{F^\top h} W(h) e^{F h} @f$, @f$ e^{2Fh} = (e^{Fh})^2 @f$.
 *
 * @return W (symmetric), or damp::nullopt on a non-finite entry or result
 * @see "Computing Integrals Involving the Matrix Exponential" (Van Loan, 1978), §3
 */
template<size_t N, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> van_loan_integral(const Matrix<N, N, T>& F, const Matrix<N, N, T>& S, T Ts) {
    T a = F.norm() * damp::abs(Ts);
    if (!damp::isfinite(a)) {
        return damp::nullopt;
    }
    //! a is finite, so it reaches ½ within max_exponent + 1 halvings.
    constexpr int max_halvings = damp::numeric_limits<T>::max_exponent + 1;
    int           s = 0;
    T             h = Ts;
    while (a > T{0.5} && s < max_halvings) {
        a *= T{0.5};
        h *= T{0.5};
        ++s;
    }
    Matrix<2 * N, 2 * N, T> Z{};
    Z.template block<N, N>(0, 0) = F.t() * (-h);
    Z.template block<N, N>(0, N) = S * h;
    Z.template block<N, N>(N, N) = F * h;
    const auto E = mat::expm(Z);
    if (!E) {
        return damp::nullopt;
    }
    Matrix<N, N, T> Phi = E->template block<N, N>(N, N);
    Matrix<N, N, T> W = Phi.t() * E->template block<N, N>(0, N);
    for (int k = 0; k < s; ++k) {
        W = W + ((Phi.t() * W) * Phi);
        Phi = Phi * Phi;
    }
    W = (W + W.t()) * T{0.5};
    if (!mat::all_of(W, [](T v) { return damp::isfinite(v); })) {
        return damp::nullopt; // an unstable mode overflowed over the period
    }
    return W;
}

} // namespace detail

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
 * - Q symmetric positive semidefinite, and Q − N R⁻¹ Nᵀ positive semidefinite
 * - R symmetric positive definite (CARE forms G = B R⁻¹ Bᵀ)
 *
 * Existence of the stabilizing solution additionally requires (A, B)
 * stabilizable and no mode on the imaginary axis unseen by Q; those are not
 * pre-screened, only used to classify a failure. care() returns damp::nullopt when
 * the Schur reduction does not converge, the stable subspace cannot be isolated,
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
    const auto sol = detail::solve_care(A, B, Q, R, N);
    if (sol.status != design::DesignStatus::ok) {
        return damp::nullopt;
    }
    return sol.X;
}

namespace design {

/**
 * @struct LQRResult
 * @brief The LQR gain: design::lqr on a plant, discrete_lqr and continuous_lqr on given matrices
 *
 * Contains the optimal gain K, Riccati solution S, and closed-loop poles. A discrete one runs
 * as an LQR (`LQR law{res}`).
 *
 * @note Compare with MATLAB®'s [K,S,P] = dlqr(...) output structure.
 * @see "Optimal Control" (Anderson & Moore, 1990), §4.3
 */
template<size_t NX, size_t NU, typename T = double>
struct [[nodiscard]] LQRResult {
    Matrix<NU, NX, T>                 K{};                                 ///< Optimal gain: u = −Kx
    Matrix<NX, NX, T>                 S{};                                 ///< Riccati solution (positive semidefinite)
    damp::array<damp::complex<T>, NX> e{};                                 ///< Closed-loop poles (eigenvalues of A − BK)
    bool                              success{false};                      ///< status == DesignStatus::ok
    DesignStatus                      status{DesignStatus::not_converged}; ///< ok, or why the design failed
    T                                 Ts{0};                               ///< Sample time of the law [s]: the plant's (design::lqr); 0 when it has none (continuous, or discrete_lqr on bare matrices)
    bool                              discrete{false};                     ///< A discrete law (e in the z-plane): discrete_lqr, or design::lqr on a discrete plant or with a Ts

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The result in scalar U. */
    template<typename U>
    [[nodiscard]] constexpr LQRResult<NX, NU, U> as() const {
        return {K.template as<U>(), S.template as<U>(), array_as<damp::complex<U>>(e), success, status, static_cast<U>(Ts), discrete};
    }

    /**
     * @brief Convert the LQR law to a static state-space block: no states, input x, output u = −Kx
     *
     * @return StateSpace with 0 states, NX inputs (x), NU outputs (u)
     */
    [[nodiscard]] constexpr StateSpace<0, NX, NU, T> to_ss() const {
        return StateSpace<0, NX, NU, T>{
            .A = Matrix<0, 0, T>{},
            .B = Matrix<0, NX, T>{},
            .C = Matrix<NU, 0, T>{},
            .D = -K,
            .Ts = Ts,
        };
    }

    /**
     * @brief Check if the closed-loop system is stable
     *
     * Discrete: every pole strictly inside the unit circle. Continuous: every pole strictly in the
     * left half plane. False for a failed design.
     */
    [[nodiscard]] constexpr bool is_stable() const {
        if (!success) {
            return false; // a failed design has no closed loop (e is zero)
        }
        for (size_t i = 0; i < NX; ++i) {
            if (discrete ? (e[i].abs() >= T{1}) : (e[i].real() >= T{0})) {
                return false;
            }
        }
        return true;
    }
};

/**
 * @struct LQRCost
 * @brief A quadratic cost on (x, u): state cost Q, input cost R, cross term N
 *
 * What discretize_lqr_cost returns: the discrete weights of a continuous cost integral.
 */
template<size_t NX, size_t NU, typename T = double>
struct LQRCost {
    Matrix<NX, NX, T> Q{}; ///< State cost
    Matrix<NU, NU, T> R{}; ///< Input cost
    Matrix<NX, NU, T> N{}; ///< Cross-term cost
};

/**
 * @brief Discrete-time Linear-Quadratic Regulator design
 *
 * Computes the optimal state-feedback gain @f$ K @f$ that minimizes
 * @f[
 *   J = \sum_k \bigl( x_k^\top Q x_k + u_k^\top R u_k + 2 x_k^\top N u_k \bigr)
 * @f]
 * subject to the discrete-time dynamics @f$ x_{k+1} = A x_k + B u_k @f$.
 *
 * The gain is applied as @f$ u = -Kx @f$. The solution is found via the Discrete
 * Algebraic Riccati Equation (DARE). Bare matrices carry no sample time: the result is
 * @c discrete with @c Ts = 0; design::lqr on a plant sets the plant's Ts.
 *
 * @note Compare with MATLAB®'s dlqr(A, B, Q, R, N).
 *
 * @see dare() for the underlying Riccati solver
 * @see lqr() to design from a plant, continuous or discrete
 * @see "Optimal Control" (Anderson & Moore, 1990), Chapter 4
 *
 * @param A  State transition matrix (NX × NX)
 * @param B  Control input matrix (NX × NU)
 * @param Q  State cost matrix (NX × NX, positive semidefinite)
 * @param R  Input cost matrix (NU × NU, positive definite)
 * @param N  Cross-term cost matrix (NX × NU, default: zero)
 * @return LQRResult with gain K, Riccati solution S, closed-loop poles and status
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr LQRResult<NX, NU, T> discrete_lqr(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N = Matrix<NX, NU, T>{}
) {
    LQRResult<NX, NU, T> out{.discrete = true};
    if (!mat::cholesky(R).has_value()) {
        return damp::detail::set_status(out, DesignStatus::bad_weights); // R not positive definite, as continuous_lqr
    }
    const auto sol = damp::detail::solve_dare(A, B, Q, R, N, DareMethod::Auto);
    if (sol.status != DesignStatus::ok) {
        return damp::detail::set_status(out, sol.status);
    }
    const auto K = lqr_gain(A, B, sol.X, R, N);
    if (!K) {
        return damp::detail::set_status(out, DesignStatus::singular); // R + BᵀSB not positive definite
    }
    out.K = K.value();
    out.S = sol.X;
    out.e = stability::closed_loop_poles(A, B, out.K);
    return damp::detail::set_status(out, DesignStatus::ok);
}

/** @brief discrete_lqr on a cost: the DARE with @p cost's Q, R and N. */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr LQRResult<NX, NU, T> discrete_lqr(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B, const LQRCost<NX, NU, T>& cost) {
    return discrete_lqr(A, B, cost.Q, cost.R, cost.N);
}

/**
 * @brief Discretize a continuous LQR cost integral over one sample (Van Loan)
 *
 * Maps the continuous running cost
 * @f[
 *   J = \int_0^\infty \big( x^\top Q x + 2 x^\top N u + u^\top R u \big)\, dt
 * @f]
 * to its exact discrete equivalent @f$ \sum (x^\top Q_d x + 2 x^\top N_d u + u^\top R_d u) @f$
 * for a zero-order-hold input. Naively reusing the continuous @f$ (Q,R,N) @f$ on
 * the discretized dynamics is only first-order accurate in @f$ T_s @f$; this is exact.
 *
 * Augmenting the held input as constant states @f$ \bar A = [A\;B;\,0\;0] @f$ with
 * weight @f$ \bar Q = [Q\;N;\,N^\top\;R] @f$, the discrete weights follow from a
 * single matrix exponential (Van Loan, 1978):
 * @f[
 *   \exp\!\left( \begin{bmatrix} -\bar A^\top & \bar Q \\ 0 & \bar A \end{bmatrix} T_s \right)
 *     = \begin{bmatrix} M_{11} & M_{12} \\ 0 & M_{22} \end{bmatrix}, \quad
 *   \begin{bmatrix} Q_d & N_d \\ N_d^\top & R_d \end{bmatrix} = M_{22}^\top M_{12}.
 * @f]
 * The exponential is taken over a fraction of the period short enough that the product does
 * not cancel, then doubled up to @f$ T_s @f$ (detail::van_loan_integral()), so a model that
 * keeps a fast stable mode (an electrical pole under a slowly sampled outer loop) gets its
 * weights to working precision.
 *
 * @note This is the cost-discretization step of MATLAB®'s lqrd(A, B, Q, R, Ts).
 * @see "Computing Integrals Involving the Matrix Exponential" (Van Loan, 1978),
 *      IEEE TAC 23(3), https://doi.org/10.1109/TAC.1978.1101743
 *
 * @param A   State transition matrix (continuous-time, NX × NX)
 * @param B   Control input matrix (continuous-time, NX × NU)
 * @param Q   Continuous state cost (NX × NX, positive semidefinite)
 * @param R   Continuous input cost (NU × NU, positive definite)
 * @param Ts  Sampling time [s]
 * @param N   Continuous cross-term cost (NX × NU, default: zero)
 * @return Discrete cost weights {Q_d, R_d, N_d}, or damp::nullopt on a non-finite entry or
 *         result
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr damp::optional<LQRCost<NX, NU, T>> discretize_lqr_cost(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    T                        Ts,
    const Matrix<NX, NU, T>& N = Matrix<NX, NU, T>{}
) {
    constexpr size_t NA = NX + NU; // augmented (state + held input) dimension

    //! Augmented dynamics Ā = [A B; 0 0] and weight Q̄ = [Q N; Nᵀ R]
    const Matrix<NU, NX, T> Nt = N.t();
    Matrix<NA, NA, T>       Abar{};
    Abar.template block<NX, NX>(0, 0) = A;
    Abar.template block<NX, NU>(0, NX) = B;

    Matrix<NA, NA, T> Qbar{};
    Qbar.template block<NX, NX>(0, 0) = Q;
    Qbar.template block<NX, NU>(0, NX) = N;
    Qbar.template block<NU, NX>(NX, 0) = Nt;
    Qbar.template block<NU, NU>(NX, NX) = R;

    //! Q̄_d = ∫₀^Ts e^{Āᵀτ} Q̄ e^{Āτ} dτ, symmetric
    const auto W = ::damp::detail::van_loan_integral(Abar, Qbar, Ts);
    if (!W) {
        return damp::nullopt;
    }
    const Matrix<NA, NA, T>& Qbar_d = *W;

    return LQRCost<NX, NU, T>{
        Qbar_d.template block<NX, NX>(0, 0),
        Qbar_d.template block<NU, NU>(NX, NX),
        Qbar_d.template block<NX, NU>(0, NX),
    };
}

/**
 * @brief Continuous-time Linear-Quadratic Regulator design
 *
 * Computes the optimal state-feedback gain @f$ K @f$ that minimizes
 * @f[
 *   J = \int_0^\infty \bigl( x^\top Q x + u^\top R u + 2 x^\top N u \bigr)\, dt
 * @f]
 * subject to the continuous-time dynamics @f$ \dot x = A x + B u @f$, applied as
 * @f$ u = -Kx @f$. Solves the Continuous Algebraic Riccati Equation (CARE) for
 * the stabilizing @f$ S @f$, then forms the gain via
 * @f$ R K^\top = B^\top S + N^\top @f$ (solve, not explicit @f$ R^{-1} @f$).
 *
 * @note Compare with MATLAB®'s [K,S,e] = lqr(A, B, Q, R, N).
 * @note Unlike @ref discrete_lqr, the returned `e` are continuous-time
 *       (s-plane) closed-loop poles of @f$ (A - BK) @f$ and @c Ts is 0,
 *       so @c is_stable tests @f$ \mathrm{Re}(e) < 0 @f$.
 *
 * @see care() for the underlying Riccati solver
 * @see lqr() for the plant form and the sampled-data (digital) design
 * @see "Optimal Control" (Anderson & Moore, 1990), §3.3
 *
 * @param A  State matrix (NX × NX)
 * @param B  Control input matrix (NX × NU)
 * @param Q  State cost matrix (NX × NX, positive semidefinite)
 * @param R  Input cost matrix (NU × NU, positive definite)
 * @param N  Cross-term cost matrix (NX × NU, default: zero)
 * @return LQRResult with gain K, CARE solution S, continuous-time poles and status
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr LQRResult<NX, NU, T> continuous_lqr(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N = Matrix<NX, NU, T>{}
) {
    LQRResult<NX, NU, T> out{};
    const auto           sol = damp::detail::solve_care(A, B, Q, R, N);
    if (sol.status != DesignStatus::ok) {
        return damp::detail::set_status(out, sol.status);
    }
    // K = R⁻¹(BᵀS + Nᵀ) as the solve R·K = BᵀS + Nᵀ, not an explicit R⁻¹
    const auto K = mat::lu_solve(R, (B.transpose() * sol.X) + N.transpose());
    if (!K) {
        return damp::detail::set_status(out, DesignStatus::bad_weights); // R singular
    }
    out.K = K.value();
    out.S = sol.X;
    out.e = stability::closed_loop_poles(A, B, out.K);
    return damp::detail::set_status(out, DesignStatus::ok);
}

} // namespace design

} // namespace damp
