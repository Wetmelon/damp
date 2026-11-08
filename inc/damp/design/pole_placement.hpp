// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file pole_placement.hpp
 * @brief State-feedback pole placement: design::place and design::place_jordan
 *
 * @defgroup pole_placement Robust pole placement (place)
 * @brief Multi-input place (Kautsky–Nichols–Van Dooren) and arbitrary Jordan-structure placement
 *
 * For a controllable pair (A, B) and a desired closed-loop spectrum, find a
 * state-feedback gain K such that the eigenvalues of (A − B·K) equal the desired
 * poles. A multi-input system has extra freedom beyond merely assigning the
 * spectrum; the Kautsky–Nichols–Van Dooren (KNV) method spends that freedom to
 * make the eigenvector matrix X as well-conditioned as possible, which minimizes
 * the sensitivity of the placed poles to perturbations in A, B, and K.
 *
 * Method 0 (orthogonal) is implemented: each eigenvector is repeatedly re-chosen
 * within its admissible subspace to be maximally orthogonal to the span of the
 * others, driving κ(X) down. The spectrum is assigned *exactly* regardless of
 * how far the conditioning sweep runs (placement correctness depends only on X
 * being invertible); the sweep only improves robustness.
 *
 * ## Time domain
 *
 * Both forms take the plant as a StateSpace, and the domain follows it:
 *
 * | Call | Plant | Poles | Gain |
 * |------|-------|-------|------|
 * | `place(sys, poles)` | continuous (`sys.Ts == 0`) | s-plane | continuous K for ẋ = (A − BK)x |
 * | `place(sys, poles)` | discrete (`sys.Ts > 0`) | z-plane | discrete K for x⁺ = (A − BK)x |
 * | `place(sys, poles, Ts)` | continuous, sampled at Ts | s-plane | discrete K on the ZOH model, poles mapped by z = e^{s·Ts} |
 *
 * A trailing `Ts` only means "this continuous plant, sampled at Ts". A discrete plant
 * with a trailing `Ts` is refused (DesignStatus::invalid_sample_time), so a discrete
 * model is never discretized a second time.
 *
 * The result is a PlaceResult, a gain design that LQR runs: `LQR law{res}`
 * deploys it, and a continuous result is refused there under constant evaluation.
 *
 * Observer gains (A − LC) are the dual: design::place_observer in estimation/fixed_gain_observer.hpp.
 *
 * @see J. Kautsky, N. K. Nichols, P. Van Dooren, "Robust pole assignment in
 *      linear state feedback," Int. J. Control 41(5), 1985,
 *      https://doi.org/10.1080/00207178508933420
 */

#include "damp/backend.hpp"
#include "damp/design/stability.hpp"
#include "damp/design/status.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/decomposition.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/functions.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace design {

namespace detail {

/** @brief Relative tolerance of the pole tests below: 64·ε, so a verdict does not depend on the units. */
template<typename T>
[[nodiscard]] constexpr T pole_rel_tol() {
    return T{64} * eps<T>();
}

/** @brief a and b are the same requested eigenvalue: |a − b| ≤ 64·ε·max(|a|, |b|). */
template<typename T>
[[nodiscard]] constexpr bool same_eigenvalue(damp::complex<T> a, damp::complex<T> b) {
    return damp::abs(a - b) <= (pole_rel_tol<T>() * damp::max(damp::abs(a), damp::abs(b)));
}

/**
 * @brief a and b are one cluster: |a − b| ≤ √ε·max(|a|, |b|)
 *
 * Closer than this, a single-input KNV eigenvector basis loses half the digits, and the computed
 * closed-loop eigenvalues spread as a Jordan chain's do.
 */
template<typename T>
[[nodiscard]] constexpr bool same_cluster(damp::complex<T> a, damp::complex<T> b) {
    return damp::abs(a - b) <= (damp::sqrt(eps<T>()) * damp::max(damp::abs(a), damp::abs(b)));
}

/** @brief λ is a real pole: |Im λ| ≤ 64·ε·|λ|. */
template<typename T>
[[nodiscard]] constexpr bool is_real_pole(damp::complex<T> lam) {
    return damp::abs(lam.imag()) <= (pole_rel_tol<T>() * damp::abs(lam));
}

/** @brief Every pole is real (is_real_pole). */
template<size_t N, typename T>
[[nodiscard]] constexpr bool all_real_poles(const damp::array<damp::complex<T>, N>& poles) {
    for (size_t i = 0; i < N; ++i) {
        if (!is_real_pole(poles[i])) {
            return false;
        }
    }
    return true;
}

/**
 * @brief KNV Method 0 placement of a real spectrum.
 *
 * @return Gain K (NU×NX), or nullopt if B is rank-deficient or the assigned
 *         eigenvectors are linearly dependent (a pole repeated more than NU times).
 */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NU, NX, T>> knv_place(
    const Matrix<NX, NX, T>&  A,
    const Matrix<NX, NU, T>&  B,
    const damp::array<T, NX>& poles
) {
    static_assert(NU <= NX, "place requires NU <= NX (no more inputs than states)");

    // Build the closed-loop matrix A − B·K = X·Λ·X⁻¹ from a chosen eigenvector
    // basis X, then recover K. The eigenvector matrix is assembled below.
    Matrix<NX, NX, T> X;

    if constexpr (NU == NX) {
        // Square, full-rank B: any spectrum is assignable with X = I. Each
        // standard basis vector is a valid eigenvector, so A − B·K = diag(poles).
        X = Matrix<NX, NX, T>::identity();
    } else {
        // QR of B: U0 = range(B), U1 = its orthogonal complement (U1ᵀ·B = 0).
        const auto qrB = mat::full_qr(B);

        // U1 = last (NX − NU) columns of Q.
        constexpr size_t  NC = NX - NU; // complement dimension
        Matrix<NX, NC, T> U1;
        for (size_t i = 0; i < NX; ++i) {
            for (size_t c = 0; c < NC; ++c) {
                U1(i, c) = qrB.Q(i, NU + c);
            }
        }

        // For each eigenvalue λ_j, the admissible eigenvectors are the null space
        // of U1ᵀ(A − λ_j I): vectors x with (A − λ_j I)x ∈ range(B). A basis S_j
        // (NX×NU) is the trailing NU columns of the full QR of (A − λ_j I)ᵀ·U1.
        damp::array<Matrix<NX, NU, T>, NX> S{};
        for (size_t j = 0; j < NX; ++j) {
            // Mtʲ = (A − λ_j I)ᵀ · U1   (NX×NC), whose left null space (the
            // trailing NU columns of its full-Q) is the admissible space S_j.
            Matrix<NX, NC, T> Mt;
            for (size_t r = 0; r < NX; ++r) {
                for (size_t c = 0; c < NC; ++c) {
                    T acc = T{0};
                    for (size_t k = 0; k < NX; ++k) {
                        const T aki = (A(k, r) - (k == r ? poles[j] : T{0})); // (A − λI)ᵀ at (r,k)
                        acc += aki * U1(k, c);
                    }
                    Mt(r, c) = acc;
                }
            }
            const auto qrM = mat::full_qr(Mt);
            for (size_t i = 0; i < NX; ++i) {
                for (size_t s = 0; s < NU; ++s) {
                    S[j](i, s) = qrM.Q(i, NC + s);
                }
            }
        }

        // Initialize each eigenvector as the first admissible basis vector.
        for (size_t j = 0; j < NX; ++j) {
            for (size_t i = 0; i < NX; ++i) {
                X(i, j) = S[j](i, 0);
            }
        }

        // KNV Method 0: sweep, re-choosing each x_j ∈ S_j to be maximally
        // orthogonal to the span of the other columns. q⊥ = the 1-D orthogonal
        // complement of the others (trailing column of their full-Q); the new
        // x_j is q⊥ projected onto S_j and normalized.
        constexpr size_t max_sweeps = 40;
        const T          tol = T{64} * eps<T>();
        for (size_t sweep = 0; sweep < max_sweeps; ++sweep) {
            T max_change = T{0};
            for (size_t j = 0; j < NX; ++j) {
                // Others = X with column j removed (NX × (NX−1)).
                Matrix<NX, NX - 1, T> others;
                for (size_t i = 0; i < NX; ++i) {
                    size_t cc = 0;
                    for (size_t c = 0; c < NX; ++c) {
                        if (c == j) {
                            continue;
                        }
                        others(i, cc) = X(i, c);
                        ++cc;
                    }
                }
                const auto qrO = mat::full_qr(others);
                // q⊥ = trailing column of Q (orthogonal to span(others)).
                damp::array<T, NX> qperp{};
                for (size_t i = 0; i < NX; ++i) {
                    qperp[i] = qrO.Q(i, NX - 1);
                }
                // Project q⊥ onto S_j: x_new = S_j (S_jᵀ q⊥).
                damp::array<T, NU> coeffs{};
                for (size_t s = 0; s < NU; ++s) {
                    T acc = T{0};
                    for (size_t i = 0; i < NX; ++i) {
                        acc += S[j](i, s) * qperp[i];
                    }
                    coeffs[s] = acc;
                }
                damp::array<T, NX> xnew{};
                T                  nrm_sq = T{0};
                for (size_t i = 0; i < NX; ++i) {
                    T acc = T{0};
                    for (size_t s = 0; s < NU; ++s) {
                        acc += S[j](i, s) * coeffs[s];
                    }
                    xnew[i] = acc;
                    nrm_sq += acc * acc;
                }
                const T nrm = damp::sqrt(nrm_sq);
                if (nrm < tol) {
                    continue; // degenerate; keep the current vector
                }
                // Track change (sign-agnostic) and commit the normalized vector.
                T diff_pos = T{0};
                T diff_neg = T{0};
                for (size_t i = 0; i < NX; ++i) {
                    const T u = xnew[i] / nrm;
                    diff_pos += (u - X(i, j)) * (u - X(i, j));
                    diff_neg += (u + X(i, j)) * (u + X(i, j));
                    xnew[i] = u;
                }
                const T change = damp::sqrt(diff_pos < diff_neg ? diff_pos : diff_neg);
                if (change > max_change) {
                    max_change = change;
                }
                for (size_t i = 0; i < NX; ++i) {
                    X(i, j) = xnew[i];
                }
            }
            if (max_change < tol) {
                break;
            }
        }
    }

    // X has unit-norm columns, so an absolute eps floor is a relative one: a
    // dependent basis (uncontrollable mode, multiplicity > NU) shows up as a
    // roundoff-sized pivot that LU alone would accept.
    if (mat::rank(X, static_cast<T>(NX) * T{64} * eps<T>()) < NX) {
        return damp::nullopt;
    }

    // Recover the gain. A − B·K = X·Λ·X⁻¹, so B·K = A − X·Λ·X⁻¹.
    // M = X·Λ·X⁻¹ ⇔ M·X = X·Λ: solve Xᵀ Mᵀ = (X·Λ)ᵀ rather than forming X⁻¹.
    Matrix<NX, NX, T> XL;
    for (size_t i = 0; i < NX; ++i) {
        for (size_t j = 0; j < NX; ++j) {
            XL(i, j) = X(i, j) * poles[j];
        }
    }
    const auto Mt_opt = mat::solve(X.transpose(), XL.transpose());
    if (!Mt_opt) {
        return damp::nullopt; // eigenvectors dependent (e.g. multiplicity > NU)
    }
    const Matrix<NX, NX, T> M = Mt_opt.value().transpose();
    const Matrix<NX, NX, T> A_minus_M = A - M;

    if constexpr (NU == NX) {
        // K = B⁻¹·(A − M) by solving B·K = (A − M).
        const auto K = mat::solve(B, A_minus_M);
        if (!K) {
            return damp::nullopt; // B not invertible
        }
        return *K;
    } else {
        // B = U0·Z with Z the top NU×NU block of R, so U0ᵀ·B = Z and
        // K = Z⁻¹·U0ᵀ·(A − M).
        const auto        qrB = mat::full_qr(B);
        Matrix<NX, NU, T> U0;
        for (size_t i = 0; i < NX; ++i) {
            for (size_t c = 0; c < NU; ++c) {
                U0(i, c) = qrB.Q(i, c);
            }
        }
        Matrix<NU, NU, T> Z;
        for (size_t r = 0; r < NU; ++r) {
            for (size_t c = 0; c < NU; ++c) {
                Z(r, c) = qrB.R(r, c);
            }
        }
        // K = Z⁻¹·U0ᵀ·(A − M) by solving Z·K = U0ᵀ·(A − M).
        const auto K = mat::solve(Z, Matrix<NU, NX, T>(U0.transpose() * A_minus_M));
        if (!K) {
            return damp::nullopt; // B rank-deficient
        }
        return *K;
    }
}

/**
 * @brief KNV placement of a spectrum with complex-conjugate pairs.
 *
 * Each conjugate pair (σ ± jω) is assigned in real arithmetic via a real
 * eigenvector pair (Re v, Im v) and a 2×2 real block [[σ, ω], [−ω, σ]] in the
 * closed-loop matrix, so K stays real. Complex poles must be supplied as adjacent
 * conjugate pairs. An all-real spectrum forwards to the real-pole routine (which
 * also runs the Method-0 conditioning sweep); the complex path assigns each
 * eigenvector from the orthonormal admissible basis without the sweep.
 *
 * @return Gain K, or nullopt if B is rank-deficient, the poles are not in valid
 *         conjugate pairs, or the eigenvectors are dependent.
 */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NU, NX, T>> knv_place(
    const Matrix<NX, NX, T>&                 A,
    const Matrix<NX, NU, T>&                 B,
    const damp::array<damp::complex<T>, NX>& poles
) {
    static_assert(NU <= NX, "place requires NU <= NX (no more inputs than states)");

    // All-real spectrum → reuse the conditioned real-pole routine.
    if (all_real_poles(poles)) {
        damp::array<T, NX> rp{};
        for (size_t i = 0; i < NX; ++i) {
            rp[i] = poles[i].real();
        }
        return knv_place(A, B, rp);
    }

    // Build the real block-diagonal Λ from the (paired) spectrum: a real pole is
    // a 1×1 entry, a conjugate pair (σ ± jω) a 2×2 block [[σ, ω], [−ω, σ]].
    // Returns false if the poles are not validly paired.
    Matrix<NX, NX, T> Lambda = Matrix<NX, NX, T>::zeros();
    {
        size_t j = 0;
        while (j < NX) {
            if (is_real_pole(poles[j])) {
                Lambda(j, j) = poles[j].real();
                ++j;
            } else {
                if (j + 1 >= NX) {
                    return damp::nullopt; // dangling complex pole
                }
                const T sg = poles[j].real();
                const T w = poles[j].imag();
                if (!same_eigenvalue(poles[j + 1], damp::conj(poles[j]))) {
                    return damp::nullopt; // not a conjugate pair
                }
                Lambda(j, j) = sg;
                Lambda(j, j + 1) = w;
                Lambda(j + 1, j) = -w;
                Lambda(j + 1, j + 1) = sg;
                j += 2;
            }
        }
    }

    if constexpr (NU == NX) {
        // X = I ⇒ A − B·K = Λ, so K = B⁻¹·(A − Λ) by solving B·K = (A − Λ).
        const auto K = mat::solve(B, Matrix<NX, NX, T>(A - Lambda));
        if (!K) {
            return damp::nullopt;
        }
        return *K;
    } else {
        const auto        qrB = mat::full_qr(B);
        constexpr size_t  NC = NX - NU;
        Matrix<NX, NC, T> U1;
        Matrix<NX, NU, T> U0;
        for (size_t i = 0; i < NX; ++i) {
            for (size_t c = 0; c < NU; ++c) {
                U0(i, c) = qrB.Q(i, c);
            }
            for (size_t c = 0; c < NC; ++c) {
                U1(i, c) = qrB.Q(i, NU + c);
            }
        }
        Matrix<NU, NU, T> Z;
        for (size_t r = 0; r < NU; ++r) {
            for (size_t c = 0; c < NU; ++c) {
                Z(r, c) = qrB.R(r, c);
            }
        }

        // U1ᵀ(A − σI) as an NC×NX block (function of σ).
        const auto P_block = [&](size_t r, size_t c, T sigma) {
            T acc = T{0};
            for (size_t k = 0; k < NX; ++k) {
                acc += U1(k, r) * (A(k, c) - (k == c ? sigma : T{0}));
            }
            return acc;
        };

        Matrix<NX, NX, T> X;
        size_t            j = 0;
        while (j < NX) {
            if (is_real_pole(poles[j])) {
                // Real pole: x_j = first admissible vector (null space of
                // U1ᵀ(A − σI) = trailing NU columns of full_qr((A − σI)ᵀU1)).
                const T           sigma = poles[j].real();
                Matrix<NX, NC, T> Mt;
                for (size_t r = 0; r < NX; ++r) {
                    for (size_t c = 0; c < NC; ++c) {
                        Mt(r, c) = P_block(c, r, sigma);
                    }
                }
                const auto qrM = mat::full_qr(Mt);
                for (size_t i = 0; i < NX; ++i) {
                    X(i, j) = qrM.Q(i, NC);
                }
                ++j;
            } else {
                // Complex pair: admissible (Re v, Im v) from the real-stacked
                // null space of [[U1ᵀ(A−σI), ω U1ᵀ], [−ω U1ᵀ, U1ᵀ(A−σI)]].
                const T sg = poles[j].real();
                const T w = poles[j].imag();

                Matrix<2 * NX, 2 * NC, T> Ht; // (2NX × 2NC) = transpose of the stacked H
                for (size_t a = 0; a < 2 * NX; ++a) {
                    const size_t ar = a % NX;
                    const bool   a_imag = a >= NX;
                    for (size_t b = 0; b < 2 * NC; ++b) {
                        const size_t br = b % NC;
                        const bool   b_bot = b >= NC;
                        T            h;
                        if (!b_bot) {
                            h = a_imag ? (w * U1(ar, br)) : P_block(br, ar, sg);
                        } else {
                            h = a_imag ? P_block(br, ar, sg) : (-w * U1(ar, br));
                        }
                        Ht(a, b) = h;
                    }
                }
                const auto   qrH = mat::full_qr(Ht);
                const size_t col0 = (2 * NX) - (2 * NU); // first trailing (null-space) column
                for (size_t i = 0; i < NX; ++i) {
                    X(i, j) = qrH.Q(i, col0);          // Re v
                    X(i, j + 1) = qrH.Q(NX + i, col0); // Im v
                }
                j += 2;
            }
        }

        // Columns of X come from orthonormal QR factors, so an absolute eps floor
        // is a relative one: a dependent basis (unobservable / uncontrollable
        // mode) shows up as a roundoff-sized pivot that LU alone would accept.
        if (mat::rank(X, static_cast<T>(NX) * T{64} * eps<T>()) < NX) {
            return damp::nullopt; // eigenvectors dependent
        }

        // M = X·Λ·X⁻¹ ⇔ M·X = X·Λ: solve Xᵀ Mᵀ = (X·Λ)ᵀ rather than forming X⁻¹.
        const Matrix<NX, NX, T> XL = X * Lambda;
        const auto              Mt_opt = mat::solve(X.transpose(), XL.transpose());
        if (!Mt_opt) {
            return damp::nullopt; // eigenvectors dependent
        }
        const Matrix<NX, NX, T> M = Mt_opt.value().transpose();
        // K = Z⁻¹·U0ᵀ·(A − M) by solving Z·K = U0ᵀ·(A − M).
        const auto K = mat::solve(Z, Matrix<NU, NX, T>(U0.transpose() * (A - M)));
        if (!K) {
            return damp::nullopt; // B rank-deficient
        }
        return *K;
    }
}

} // namespace detail

/**
 * @brief One Jordan mini-block of a desired closed-loop spectrum.
 *
 * A block contributes `size` consecutive eigenvalues equal to `eigenvalue`,
 * coupled into a single Jordan chain of that order. `size == 1` is an ordinary
 * (semisimple) eigenvalue; `size > 1` is a defective block requiring generalized
 * eigenvectors. Complex eigenvalues must be supplied as conjugate-pair blocks of
 * equal size (e.g. one block at σ+jω and one at σ−jω), with the same block orders at
 * λ and λ̄.
 */
template<typename T = double>
struct JordanBlock {
    damp::complex<T> eigenvalue{}; ///< The eigenvalue λ this block places.
    size_t           size{0};      ///< Mini-block order pᵢₖ (chain length); 0 = unused slot.
};

namespace detail {

/**
 * @brief Precomputed, K-independent data for the Klein–Moore construction.
 *
 * For each distinct processed eigenvalue (reals and the +imaginary member of each
 * conjugate pair) it stores the kernel basis N and pseudoinverse M of the pencil
 * [A − λI | B] — neither depends on the free parameter K — plus the mini-block
 * layout. assemble_vw() then turns any K into the eigenvector/input matrices,
 * which lets an optimizer sweep K without recomputing any SVDs.
 */
template<size_t NX, size_t NU, size_t NB, typename T>
struct JordanPlan {
    using Cplx = damp::complex<T>;
    static constexpr size_t NS = NX + NU;

    size_t                                   ndistinct = 0; ///< Processed eigenvalues.
    size_t                                   total_pos = 0; ///< Chain positions (used K columns).
    damp::array<bool, NB>                    is_real{};
    damp::array<size_t, NB>                  g{};          ///< Mini-blocks per eigenvalue.
    damp::array<damp::array<size_t, NU>, NB> orders{};     ///< Block orders.
    damp::array<size_t, NB>                  pos_offset{}; ///< First K column for eigenvalue.
    damp::array<size_t, NB>                  npos{};       ///< Positions (= multiplicity).
    damp::array<Matrix<NS, NU, Cplx>, NB>    N{};          ///< Kernel bases.
    damp::array<Matrix<NS, NX, Cplx>, NB>    M{};          ///< Pseudoinverses.
};

/** @brief The first n entries of a and b are equal as multisets (block orders in any order). */
template<size_t N>
[[nodiscard]] constexpr bool same_block_orders(damp::array<size_t, N> a, damp::array<size_t, N> b, size_t n) {
    const size_t m = (n < N) ? n : N;
    if constexpr (N > 1) {
        for (size_t i = 1; i < m; ++i) { // insertion sort, m ≤ N
            for (size_t j = i; j > 0 && a[j - 1] > a[j]; --j) {
                const size_t t = a[j - 1];
                a[j - 1] = a[j];
                a[j] = t;
            }
            for (size_t j = i; j > 0 && b[j - 1] > b[j]; --j) {
                const size_t t = b[j - 1];
                b[j - 1] = b[j];
                b[j] = t;
            }
        }
    }
    for (size_t i = 0; i < m; ++i) {
        if (a[i] != b[i]) {
            return false;
        }
    }
    return true;
}

/** @brief Build the K-independent plan, or nullopt if the requested structure is inadmissible. */
template<size_t NX, size_t NU, size_t NB, typename T>
[[nodiscard]] constexpr damp::optional<JordanPlan<NX, NU, NB, T>> prepare_jordan_plan(
    const Matrix<NX, NX, T>&               A,
    const Matrix<NX, NU, T>&               B,
    const damp::array<JordanBlock<T>, NB>& blocks
) {
    using Cplx = damp::complex<T>;
    constexpr size_t NS = NX + NU;

    size_t total = 0;
    for (size_t b = 0; b < NB; ++b) {
        total += blocks[b].size;
    }
    if (total != NX) {
        return damp::nullopt; // block sizes must place exactly NX eigenvalues
    }

    const Matrix<NX, NX, Cplx> Ac = A.template as<Cplx>();
    const Matrix<NX, NU, Cplx> Bc = B.template as<Cplx>();

    JordanPlan<NX, NU, NB, T> plan;
    size_t                    col = 0;
    damp::array<bool, NB>     consumed{};
    for (size_t bi = 0; bi < NB; ++bi) {
        if (consumed[bi] || blocks[bi].size == 0) {
            continue; // already grouped, or an unused (size-0) slot
        }
        const Cplx lam = blocks[bi].eigenvalue;
        const bool is_real = is_real_pole(lam);
        if (!is_real && lam.imag() < T{0}) {
            continue; // negative-imaginary member is realified with its partner
        }

        // Count then gather this eigenvalue's mini-blocks.
        size_t cnt = 0;
        for (size_t bj = 0; bj < NB; ++bj) {
            if (!consumed[bj] && blocks[bj].size > 0 && same_eigenvalue(blocks[bj].eigenvalue, lam)) {
                ++cnt;
            }
        }
        if (cnt > NU) {
            return damp::nullopt; // more chains than the kernel dimension allows
        }
        const size_t            e = plan.ndistinct;
        damp::array<size_t, NU> ords{};
        size_t                  g = 0;
        size_t                  mult = 0;
        for (size_t bj = 0; bj < NB; ++bj) {
            if (!consumed[bj] && blocks[bj].size > 0 && same_eigenvalue(blocks[bj].eigenvalue, lam)) {
                ords[g] = blocks[bj].size;
                mult += blocks[bj].size;
                ++g;
                consumed[bj] = true;
            }
        }
        // A real K mirrors the Jordan structure at λ onto λ̄: the conjugate's blocks must have
        // the same orders.
        if (!is_real) {
            damp::array<size_t, NU> conj_ords{};
            size_t                  gc = 0;
            for (size_t bj = 0; bj < NB; ++bj) {
                if (!consumed[bj] && blocks[bj].size > 0 && same_eigenvalue(blocks[bj].eigenvalue, damp::conj(lam))) {
                    if (gc >= NU) {
                        return damp::nullopt; // more chains at λ̄ than at λ
                    }
                    conj_ords[gc] = blocks[bj].size;
                    ++gc;
                    consumed[bj] = true;
                }
            }
            if (gc != g || !same_block_orders(ords, conj_ords, g)) {
                return damp::nullopt; // unmatched conjugate pair
            }
        }

        // Pencil S = [A − λI | B] (NX × NS). Its kernel basis N and a right
        // inverse M come from one QR of Sᴴ. Reachability makes S full row rank
        // for *every* λ (PBH test), so this single path also covers λ ∈ spec(A)
        // with no special case — unlike a factorization of (A − λI) alone.
        Matrix<NX, NS, Cplx> S;
        for (size_t i = 0; i < NX; ++i) {
            for (size_t c = 0; c < NX; ++c) {
                S(i, c) = Ac(i, c) - (i == c ? lam : Cplx{0});
            }
            for (size_t c = 0; c < NU; ++c) {
                S(i, NX + c) = Bc(i, c);
            }
        }
        // Sᴴ = Q·R: Q (NS×NS) unitary, R (NS×NX) upper-triangular with leading
        // NX×NX block R1 invertible iff S is full row rank.
        const auto           qr = mat::full_qr(Matrix<NS, NX, Cplx>(S.conjugate_transpose()));
        Matrix<NX, NX, Cplx> R1;
        for (size_t i = 0; i < NX; ++i) {
            for (size_t j = 0; j < NX; ++j) {
                R1(i, j) = qr.R(i, j);
            }
        }
        // Kernel N = trailing NU columns of Q (orthonormal basis of ker S).
        Matrix<NS, NU, Cplx> N;
        // Q1 = leading NX columns of Q. Right inverse of S is S⁺ = Q1·(R1⁻¹)ᴴ;
        // M = Q1·(R1ᴴ)⁻¹ ⇔ R1·Mᴴ = Q1ᴴ, solved without forming R1⁻¹.
        Matrix<NS, NX, Cplx> Q1;
        for (size_t i = 0; i < NS; ++i) {
            for (size_t j = 0; j < NX; ++j) {
                Q1(i, j) = qr.Q(i, j);
            }
            for (size_t k = 0; k < NU; ++k) {
                N(i, k) = qr.Q(i, NX + k);
            }
        }
        const auto MH_opt = mat::solve(R1, Q1.conjugate_transpose());
        if (!MH_opt) {
            return damp::nullopt; // (A,B) not reachable at λ
        }

        plan.is_real[e] = is_real;
        plan.g[e] = g;
        plan.orders[e] = ords;
        plan.N[e] = N;
        plan.M[e] = Matrix<NS, NX, Cplx>(MH_opt.value().conjugate_transpose());
        plan.pos_offset[e] = col;
        plan.npos[e] = mult;
        col += mult;
        ++plan.ndistinct;
    }
    plan.total_pos = col;
    if (col == 0) {
        return damp::nullopt;
    }
    return plan;
}

/** @brief Eigenvector matrix V (state parts, real) and input matrix W (input parts) for a given K. */
template<size_t NX, size_t NU, size_t NB, typename T>
constexpr void assemble_vw(
    const JordanPlan<NX, NU, NB, T>&        plan,
    const Matrix<NU, NX, damp::complex<T>>& K,
    Matrix<NX, NX, T>&                      V,
    Matrix<NU, NX, T>&                      W
) {
    using Cplx = damp::complex<T>;
    constexpr size_t NS = NX + NU;
    V = Matrix<NX, NX, T>::zeros();
    W = Matrix<NU, NX, T>::zeros();
    size_t col = 0;

    for (size_t e = 0; e < plan.ndistinct; ++e) {
        const bool           is_real = plan.is_real[e];
        const auto&          N = plan.N[e];
        const auto&          M = plan.M[e];
        size_t               pcol = plan.pos_offset[e]; // K-column cursor
        Matrix<NS, NX, Cplx> chains{};
        size_t               nchain = 0;
        for (size_t k = 0; k < plan.g[e]; ++k) {
            damp::array<Cplx, NS> prev{};
            for (size_t l = 0; l < plan.orders[e][k]; ++l) {
                damp::array<Cplx, NS> h{};
                for (size_t i = 0; i < NS; ++i) {
                    Cplx acc{0};
                    for (size_t r = 0; r < NU; ++r) {
                        acc += N(i, r) * K(r, pcol); // N · K(l)
                    }
                    h[i] = acc;
                }
                if (l > 0) {
                    for (size_t i = 0; i < NS; ++i) {
                        Cplx acc{0};
                        for (size_t r = 0; r < NX; ++r) {
                            acc += M(i, r) * prev[r]; // + M · overp(prev)
                        }
                        h[i] += acc;
                    }
                }
                ++pcol;
                if (is_real) {
                    for (size_t i = 0; i < NX; ++i) {
                        V(i, col) = h[i].real();
                    }
                    for (size_t i = 0; i < NU; ++i) {
                        W(i, col) = h[NX + i].real();
                    }
                    ++col;
                } else {
                    for (size_t i = 0; i < NS; ++i) {
                        chains(i, nchain) = h[i];
                    }
                    ++nchain;
                }
                prev = h;
            }
        }
        if (!is_real) {
            // Realify: real parts then imaginary parts of the chain columns.
            for (size_t t = 0; t < nchain; ++t) {
                for (size_t i = 0; i < NX; ++i) {
                    V(i, col) = chains(i, t).real();
                }
                for (size_t i = 0; i < NU; ++i) {
                    W(i, col) = chains(NX + i, t).real();
                }
                ++col;
            }
            for (size_t t = 0; t < nchain; ++t) {
                for (size_t i = 0; i < NX; ++i) {
                    V(i, col) = chains(i, t).imag();
                }
                for (size_t i = 0; i < NU; ++i) {
                    W(i, col) = chains(NX + i, t).imag();
                }
                ++col;
            }
        }
    }
}

/** @brief The canonical (minimum-norm) parameter: chain k starts from kernel column k. */
template<size_t NX, size_t NU, size_t NB, typename T>
[[nodiscard]] constexpr Matrix<NU, NX, damp::complex<T>> canonical_kparams(
    const JordanPlan<NX, NU, NB, T>& plan
) {
    using Cplx = damp::complex<T>;
    Matrix<NU, NX, Cplx> K = Matrix<NU, NX, Cplx>::zeros();
    for (size_t e = 0; e < plan.ndistinct; ++e) {
        size_t pcol = plan.pos_offset[e];
        for (size_t k = 0; k < plan.g[e]; ++k) {
            K(k, pcol) = Cplx{1}; // h(1) = N·e_k = kernel column k
            pcol += plan.orders[e][k];
        }
    }
    return K;
}

/** @brief Gain from the canonical (minimum-norm) parameter: K = −W·V⁻¹, or nullopt if V is singular. */
template<size_t NX, size_t NU, size_t NB, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NU, NX, T>> jordan_canonical_gain(const JordanPlan<NX, NU, NB, T>& plan) {
    Matrix<NX, NX, T> V;
    Matrix<NU, NX, T> W;
    assemble_vw(plan, canonical_kparams(plan), V, W);
    // K = −W·V⁻¹ ⇔ K·V = −W: solve Vᵀ Kᵀ = −Wᵀ rather than forming V⁻¹.
    const auto Kt_opt = mat::solve(V.transpose(), Matrix<NX, NU, T>(T{-1} * W.transpose()));
    if (!Kt_opt) {
        return damp::nullopt; // canonical parameter gave dependent eigenvectors
    }
    return Matrix<NU, NX, T>(Kt_opt.value().transpose());
}

/**
 * @brief Placement of a spectrum with repeated poles by the Klein–Moore canonical construction.
 *
 * With one input every closed-loop eigenvalue has one eigenvector, so a pole of
 * multiplicity m forms one Jordan chain of order m and the gain is unique (the
 * Ackermann gain). With NU ≥ 2 each pole is its own 1×1 block, so a pole repeated
 * m ≤ NU times gets m independent eigenvectors (the chain heads are distinct kernel
 * columns), and one repeated more than NU times is refused.
 */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NU, NX, T>> repeated_pole_place(
    const Matrix<NX, NX, T>&                 A,
    const Matrix<NX, NU, T>&                 B,
    const damp::array<damp::complex<T>, NX>& poles
) {
    damp::array<JordanBlock<T>, NX> blocks{};
    if constexpr (NU == 1) {
        damp::array<bool, NX> grouped{};
        size_t                nb = 0;
        for (size_t i = 0; i < NX; ++i) {
            if (grouped[i]) {
                continue;
            }
            size_t m = 0;
            for (size_t j = i; j < NX; ++j) {
                if (!grouped[j] && same_eigenvalue(poles[j], poles[i])) {
                    grouped[j] = true;
                    ++m;
                }
            }
            blocks[nb] = JordanBlock<T>{poles[i], m};
            ++nb;
        }
    } else {
        for (size_t i = 0; i < NX; ++i) {
            blocks[i] = JordanBlock<T>{poles[i], 1};
        }
    }
    const auto plan = prepare_jordan_plan(A, B, blocks);
    if (!plan) {
        return damp::nullopt;
    }
    return jordan_canonical_gain(*plan);
}

/** @brief true if some pole appears more than once (same_eigenvalue). */
template<size_t NX, typename T>
[[nodiscard]] constexpr bool has_repeated_pole(const damp::array<damp::complex<T>, NX>& poles) {
    for (size_t i = 0; i < NX; ++i) {
        for (size_t j = i + 1; j < NX; ++j) {
            if (same_eigenvalue(poles[j], poles[i])) {
                return true;
            }
        }
    }
    return false;
}

/** @brief true if two requested poles are one cluster (same_cluster). */
template<size_t NX, typename T>
[[nodiscard]] constexpr bool has_clustered_poles(const damp::array<damp::complex<T>, NX>& poles) {
    for (size_t i = 0; i < NX; ++i) {
        for (size_t j = i + 1; j < NX; ++j) {
            if (same_cluster(poles[j], poles[i])) {
                return true;
            }
        }
    }
    return false;
}

/** @brief Every requested pole is finite. */
template<size_t N, typename T>
[[nodiscard]] constexpr bool finite_poles(const damp::array<damp::complex<T>, N>& poles) {
    for (size_t i = 0; i < N; ++i) {
        if (!damp::isfinite(poles[i].real()) || !damp::isfinite(poles[i].imag())) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Certificate: eig(A − B·K) matches the requested poles
 *
 * Each requested pole claims the nearest unclaimed eigenvalue, which must lie within
 * 64·NX·ε^{1/(m+1)}·(‖A‖ + ‖B‖‖K‖), m the number of requested poles in its cluster (same_cluster,
 * itself included): a repeat, or poles too close for the eigensolver to tell apart from one. The
 * exponent leaves room for the ε^{1/m} spread of a Jordan chain of order m and for an
 * ill-conditioned eigenvector basis; a gain that misses the spectrum, an eigensolver that does
 * not converge, or a non-finite gain fails.
 */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr bool placement_matches(
    const Matrix<NX, NX, T>&                 A,
    const Matrix<NX, NU, T>&                 B,
    const Matrix<NU, NX, T>&                 K,
    const damp::array<damp::complex<T>, NX>& poles
) {
    const auto ev = mat::compute_eigenvalues(Matrix<NX, NX, T>(A - (B * K)));
    if (!ev.success) {
        return false;
    }
    const T size = A.norm() + (B.norm() * K.norm());
    if (!damp::isfinite(size)) {
        return false;
    }
    damp::array<bool, NX> claimed{};
    for (size_t i = 0; i < NX; ++i) {
        size_t m = 0;
        for (size_t j = 0; j < NX; ++j) {
            if (same_cluster(poles[j], poles[i])) {
                ++m;
            }
        }
        const T tol = T{64} * static_cast<T>(NX) * size * damp::exp(damp::log(eps<T>()) / static_cast<T>(m + 1));
        size_t  best = NX;
        T       best_d = T{0};
        for (size_t j = 0; j < NX; ++j) {
            const T d = damp::abs(damp::complex<T>(ev.values[j]) - poles[i]);
            if (!claimed[j] && (best == NX || d < best_d)) {
                best = j;
                best_d = d;
            }
        }
        if (best == NX || !(best_d <= tol)) {
            return false;
        }
        claimed[best] = true;
    }
    return true;
}

/**
 * @brief Single-input gain from the requested characteristic polynomial (Ackermann)
 *
 * K = e_nᵀ·𝒞⁻¹·φ(A), φ(s) = Π(s − pᵢ) with its real coefficients, 𝒞 the controllability
 * matrix; e_nᵀ𝒞⁻¹ is a solve with 𝒞ᵀ. The formula is evaluated in the time unit where ‖A‖ and
 * every |pᵢ| are at most 1 (A′ = A/α, pᵢ′ = pᵢ/α, K = α·K′), so the rank test on 𝒞′ does not
 * depend on the time unit. nullopt if 𝒞′ is numerically rank-deficient (mat::rank): a mode
 * that B reaches only below roundoff would need a gain of order 1/ε. The imaginary parts of φ's
 * coefficients are dropped, so a set that is not conjugate-closed gives a gain that misses it
 * (and fails placement_matches).
 *
 * @see "Der Entwurf linearer Regelungssysteme im Zustandsraum" (Ackermann, 1972)
 */
template<size_t NX, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<1, NX, T>> ackermann_gain(
    const Matrix<NX, NX, T>&                 A,
    const Matrix<NX, 1, T>&                  B,
    const damp::array<damp::complex<T>, NX>& poles
) {
    T alpha = A.norm();
    for (size_t i = 0; i < NX; ++i) {
        alpha = damp::max(alpha, damp::abs(poles[i]));
    }
    if (!(alpha > T{0}) || !damp::isfinite(alpha)) {
        alpha = T{1};
    }
    const T                               inv_alpha = T{1} / alpha;
    const Matrix<NX, NX, T>               As = A * inv_alpha;
    damp::array<damp::complex<T>, NX + 1> c{}; // monic φ′(s), highest power first
    c[0] = damp::complex<T>{T{1}, T{0}};
    for (size_t i = 0; i < NX; ++i) {
        const damp::complex<T> p = poles[i] * inv_alpha;
        for (size_t k = i + 1; k > 0; --k) {
            c[k] = c[k] - (p * c[k - 1]);
        }
    }
    Matrix<NX, NX, T> phi = Matrix<NX, NX, T>::identity(); // φ′(A′) by Horner
    for (size_t k = 1; k <= NX; ++k) {
        phi = Matrix<NX, NX, T>((As * phi) + (Matrix<NX, NX, T>::identity() * c[k].real()));
    }
    const Matrix<NX, NX, T> Co = stability::controllability_matrix(As, B);
    if (mat::rank(Co) < NX) {
        return damp::nullopt;
    }
    Matrix<NX, 1, T> en{};
    en(NX - 1, 0) = T{1};
    const auto y = mat::solve(Matrix<NX, NX, T>(Co.transpose()), en); // yᵀ = e_nᵀ𝒞′⁻¹
    if (!y) {
        return damp::nullopt;
    }
    return Matrix<1, NX, T>((y->transpose() * phi) * alpha);
}

/**
 * @brief The placement kernel behind place() and place_observer(): K with eig(A − B·K) = poles.
 *
 * KNV for a non-defective spectrum; a single-input pair with repeated poles, or a repeated
 * pole in a spectrum with complex pairs, goes through the Klein–Moore construction instead.
 * With one input the gain is unique, so a cluster (same_cluster), whose KNV eigenvectors are
 * nearly parallel, takes the Ackermann gain, and so does a request the other paths miss. The
 * result is certified against the request (placement_matches). Writes K and returns the
 * status: invalid_parameter for a non-finite pole, singular for a rank-deficient B,
 * not_placeable otherwise.
 */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr DesignStatus place_gain(
    const Matrix<NX, NX, T>&                 A,
    const Matrix<NX, NU, T>&                 B,
    const damp::array<damp::complex<T>, NX>& poles,
    Matrix<NU, NX, T>&                       K
) {
    static_assert(NU <= NX, "place requires NU <= NX (no more inputs than states)");
    if (!finite_poles(poles)) {
        return DesignStatus::invalid_parameter;
    }
    if (mat::rank(B) < NU) {
        return DesignStatus::singular; // rank-deficient B (C for an observer)
    }
    // The all-real KNV sweep separates the eigenvectors of a repeat; the complex path and a
    // single input cannot, so those repeats take the Jordan construction.
    const bool                        chains = has_repeated_pole(poles) && ((NU == 1) || !all_real_poles(poles));
    damp::optional<Matrix<NU, NX, T>> gain{};
    if (chains) {
        gain = repeated_pole_place(A, B, poles);
    } else if ((NU > 1) || !has_clustered_poles(poles)) {
        gain = knv_place(A, B, poles);
    }
    if constexpr (NU == 1) {
        // One input: the gain is unique, from the requested polynomial.
        if (!gain || !placement_matches(A, B, *gain, poles)) {
            gain = ackermann_gain(A, B, poles);
        }
    }
    if (!gain || !placement_matches(A, B, *gain, poles)) {
        return DesignStatus::not_placeable;
    }
    K = *gain;
    return DesignStatus::ok;
}

/** @brief z = e^{s·Ts} for one complex pole. */
template<typename T>
[[nodiscard]] constexpr damp::complex<T> pole_s_to_z(damp::complex<T> s, T Ts) {
    const T er = damp::exp(s.real() * Ts);
    const auto [sn, cs] = damp::sincos(s.imag() * Ts);
    return damp::complex<T>{er * cs, er * sn};
}

/** @brief Map each s-plane pole to z = e^{s·Ts}. */
template<size_t N, typename T>
[[nodiscard]] constexpr damp::array<damp::complex<T>, N> poles_s_to_z(const damp::array<damp::complex<T>, N>& poles_s, T Ts) {
    damp::array<damp::complex<T>, N> z{};
    for (size_t i = 0; i < N; ++i) {
        z[i] = pole_s_to_z(poles_s[i], Ts);
    }
    return z;
}

/**
 * @brief The ZOH model of a continuous plant at Ts, for the trailing-Ts design forms.
 *
 * Refuses a discrete plant (it is already sampled), and a Ts that is not positive
 * and finite, with DesignStatus::invalid_sample_time. Writes the sampled plant to
 * @p out on success.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr DesignStatus sample_plant(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        Ts,
    StateSpace<NX, NU, NY, T, NW, NV>&       out
) {
    if (sys.Ts != T{0} || !(Ts > T{0}) || !damp::isfinite(Ts)) {
        return DesignStatus::invalid_sample_time;
    }
    const auto sampled = discretize(sys, Ts, DiscretizationMethod::ZOH);
    if (!sampled) {
        return DesignStatus::singular;
    }
    out = *sampled;
    return DesignStatus::ok;
}

} // namespace detail

/**
 * @struct PlaceResult
 * @brief Pole-placement design result: the gain K, the achieved closed-loop poles e, and Ts
 *
 * A gain design, so it runs as LQR:
 * `LQR law{res.as<float>()}`, or `static constinit LQR law{res}`
 * for a discrete design. `e` holds eig(A − BK) computed from the gain, not the requested
 * poles; `Ts` is the law's sample time (0 for a continuous gain).
 *
 * @tparam NX Number of states
 * @tparam NU Number of inputs
 * @tparam T  Scalar type
 */
template<size_t NX, size_t NU, typename T = double>
struct [[nodiscard]] PlaceResult {
    Matrix<NU, NX, T>                 K{};                                 ///< The gain: u = −Kx
    damp::array<damp::complex<T>, NX> e{};                                 ///< Achieved closed-loop poles, eig(A − BK)
    bool                              success{false};                      ///< status == DesignStatus::ok
    DesignStatus                      status{DesignStatus::not_converged}; ///< ok, or why the design failed
    T                                 Ts{0};                               ///< Sample time of the law: 0 continuous, > 0 discrete

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief Whether the design succeeded and the closed loop is stable: every pole inside the unit circle (Ts > 0) or in the left half plane (Ts = 0). */
    [[nodiscard]] constexpr bool is_stable() const {
        if (!success) {
            return false;
        }
        for (size_t i = 0; i < NX; ++i) {
            if ((Ts > T{0}) ? (e[i].abs() >= T{1}) : (e[i].real() >= T{0})) {
                return false;
            }
        }
        return true;
    }

    template<typename U>
    [[nodiscard]] constexpr PlaceResult<NX, NU, U> as() const {
        return PlaceResult<NX, NU, U>{K.template as<U>(), array_as<damp::complex<U>>(e), success, status, static_cast<U>(Ts)};
    }
};

namespace detail {

/** @brief place() on a plant whose domain is already the law's: K, e and Ts from sys. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr PlaceResult<NX, NU, T> place_on(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const damp::array<damp::complex<T>, NX>& poles
) {
    PlaceResult<NX, NU, T> out{};
    out.Ts = sys.Ts;
    const DesignStatus status = place_gain(sys.A, sys.B, poles, out.K);
    if (status == DesignStatus::ok) {
        out.e = stability::closed_loop_poles(sys.A, sys.B, out.K);
    }
    ::damp::detail::set_status(out, status);
    return out;
}

} // namespace detail

/**
 * @ingroup pole_placement
 * @brief Robust pole placement (KNV): the gain K with eig(A − B·K) at @p poles, in the domain of @p sys
 *
 * A continuous plant takes s-plane poles and gives a continuous gain; a discrete plant
 * takes z-plane poles and gives a discrete gain. Multi-input (NU ≤ NX): the freedom
 * beyond the spectrum goes to eigenvector conditioning (KNV Method 0 for real spectra).
 * Complex poles come as adjacent conjugate pairs, named as complex values
 * (`design::place(sys, {damp::complex{-1.0, 2.0}, damp::complex{-1.0, -2.0}})`); a real pole
 * converts implicitly (`design::place(sys, {-1.0, -2.0})`). A pole repeated more than NU times has no
 * independent eigenvectors and is refused, except for a single-input plant, where the
 * repeat forms one Jordan chain and the gain is the unique (Ackermann) gain. Other
 * repeated structures go through place_jordan.
 *
 * Two poles are the same pole when they agree to 64·ε relative to their magnitude, so
 * distinct poles are never merged, whatever their units. For a single input, nearly
 * repeated poles (agreeing to √ε relative) take the Ackermann gain. The achieved
 * eig(A − BK) is checked against the request before the design is reported as a success.
 *
 * Status: invalid_parameter if a pole is non-finite; not_placeable if the spectrum cannot
 * be assigned (uncontrollable mode, bad conjugate pairing, a repeat beyond NU) or the
 * achieved poles miss the request; singular if B is rank-deficient.
 *
 * @param sys   Plant; only A and B are used
 * @param poles Desired closed-loop eigenvalues, adjacent conjugate pairs
 * @return PlaceResult with K, e = eig(A − BK) and Ts = sys.Ts
 *
 * @note Compare with MATLAB®'s K = place(A, B, p), and K = acker(A, B, p) for a single input.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr PlaceResult<NX, NU, T> place(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const damp::array<damp::complex<T>, NX>& poles
) {
    return detail::place_on(sys, poles);
}

/**
 * @ingroup pole_placement
 * @brief Pole placement for a continuous plant sampled at @p Ts: s-plane poles, discrete gain
 *
 * ZOH-discretizes @p sys at @p Ts, maps each pole to z = e^{s·Ts}, and places on the
 * sampled pair, so the gain deploys at Ts with the closed-loop poles chosen in rad/s.
 * A discrete @p sys is refused with DesignStatus::invalid_sample_time (it is already
 * sampled; pass z-plane poles to place(sys, poles) instead), as is Ts ≤ 0 or non-finite.
 *
 * @param sys   Continuous plant (sys.Ts == 0)
 * @param poles Desired continuous closed-loop poles (s-plane), adjacent conjugate pairs
 * @param Ts    [s] Sample period of the law
 * @return PlaceResult with the discrete K, e = eig(A_d − B_d K) (z-plane) and Ts
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr PlaceResult<NX, NU, T> place(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const damp::array<damp::complex<T>, NX>& poles,
    T                                        Ts
) {
    StateSpace<NX, NU, NY, T, NW, NV> sampled{};
    const DesignStatus                status = detail::sample_plant(sys, Ts, sampled);
    if (status != DesignStatus::ok) {
        return ::damp::detail::refused<PlaceResult<NX, NU, T>>(status);
    }
    return detail::place_on(sampled, detail::poles_s_to_z(poles, Ts));
}

/** @brief How place_jordan spends the free parameter of the Klein–Moore family (Schmid et al.). */
enum class JordanObjective : uint8_t {
    Canonical,              ///< No search: the canonical minimum-norm parameter.
    ConditionNumber,        ///< Method 1: minimize the Frobenius condition number of V.
    DepartureFromNormality, ///< Method 2: minimize the departure from normality of A − B·K.
};

/**
 * @brief Options for place_jordan
 *
 * The default places the structure with the canonical parameter and runs no search.
 * ConditionNumber / DepartureFromNormality minimize α·R(K) + (1 − α)·‖K‖²_F by BFGS.
 */
template<typename T = double>
struct JordanOptions {
    JordanObjective objective{JordanObjective::Canonical}; ///< Robustness measure R, or Canonical for no search
    T               alpha{T{1}};                           ///< Robustness/gain weight in [0,1]: 1 = robust, 0 = minimum gain
    size_t          max_iter{200};                         ///< BFGS iteration cap
};

/**
 * @struct JordanPlaceResult
 * @brief place_jordan result: the PlaceResult fields plus the achieved robustness metrics
 *
 * Every parameter in the Klein–Moore family places the structure exactly, so success
 * reports the placement; @c converged reports whether the search reached a stationary
 * point (true when no search was asked).
 */
template<size_t NX, size_t NU, typename T = double>
struct [[nodiscard]] JordanPlaceResult : PlaceResult<NX, NU, T> {
    T      cond_fro{};       ///< Achieved κ_F(V) = ‖V‖F · ‖V⁻¹‖F (eigenvalue robustness)
    T      gain_fro{};       ///< Achieved ‖K‖F (control effort)
    T      departure_fro{};  ///< Achieved δ_F(A − B·K) (departure from normality)
    size_t iterations{};     ///< BFGS iterations taken
    bool   converged{false}; ///< Search reached a stationary point (true for Canonical)

    template<typename U>
    [[nodiscard]] constexpr JordanPlaceResult<NX, NU, U> as() const {
        return JordanPlaceResult<NX, NU, U>{
            PlaceResult<NX, NU, T>::template as<U>(),
            static_cast<U>(cond_fro),
            static_cast<U>(gain_fro),
            static_cast<U>(departure_fro),
            iterations,
            converged,
        };
    }
};

namespace detail {

/** @brief The Jordan placement on (A, B): gain, metrics and search state, or nullopt if inadmissible. */
template<size_t NX, size_t NU, size_t NB, typename T>
[[nodiscard]] constexpr damp::optional<JordanPlaceResult<NX, NU, T>> jordan_place(
    const Matrix<NX, NX, T>&               A,
    const Matrix<NX, NU, T>&               B,
    const damp::array<JordanBlock<T>, NB>& blocks,
    const JordanOptions<T>&                opts
) {
    static_assert(NU <= NX, "place_jordan requires NU <= NX");
    const T               alpha = opts.alpha;
    const JordanObjective objective = opts.objective;
    const bool            search = objective != JordanObjective::Canonical;
    const size_t          max_iter = search ? opts.max_iter : size_t{0};
    using Cplx = damp::complex<T>;

    const auto plan_opt = prepare_jordan_plan(A, B, blocks);
    if (!plan_opt) {
        return damp::nullopt;
    }
    const auto& plan = plan_opt.value();

    // Σ|λ|² over all eigenvalues, for the departure-from-normality measure.
    T sum_lambda_sq = T{0};
    for (size_t b = 0; b < NB; ++b) {
        const Cplx lam = blocks[b].eigenvalue;
        sum_lambda_sq += static_cast<T>(blocks[b].size) * ((lam.real() * lam.real()) + (lam.imag() * lam.imag()));
    }

    // Real degrees of freedom: NU per real-eigenvalue chain position, 2·NU per
    // complex one (the conjugate is determined). Laid out into a flat θ vector.
    constexpr size_t MaxDof = 2 * NU * NX;
    size_t           dof = 0;
    for (size_t e = 0; e < plan.ndistinct; ++e) {
        dof += plan.npos[e] * (plan.is_real[e] ? NU : (2 * NU));
    }

    const auto theta_to_k = [&](const damp::array<T, MaxDof>& th) {
        Matrix<NU, NX, Cplx> K = Matrix<NU, NX, Cplx>::zeros();
        size_t               idx = 0;
        for (size_t e = 0; e < plan.ndistinct; ++e) {
            size_t pcol = plan.pos_offset[e];
            for (size_t p = 0; p < plan.npos[e]; ++p) {
                for (size_t r = 0; r < NU; ++r) {
                    if (plan.is_real[e]) {
                        K(r, pcol) = Cplx{th[idx], T{0}};
                        idx += 1;
                    } else {
                        K(r, pcol) = Cplx{th[idx], th[idx + 1]};
                        idx += 2;
                    }
                }
                ++pcol;
            }
        }
        return K;
    };

    // Squared Frobenius norm (sum of squares — no sqrt round-trip).
    const auto fro2 = [](const auto& Mx) {
        T s = T{0};
        for (size_t i = 0; i < Mx.rows(); ++i) {
            for (size_t j = 0; j < Mx.cols(); ++j) {
                s += Mx(i, j) * Mx(i, j);
            }
        }
        return s;
    };

    // assemble_vw is *linear* in K, so ∂V/∂θⱼ and ∂W/∂θⱼ are the constant matrices
    // produced by a unit parameter at component j. Precompute them once; the
    // analytic gradient is then a Frobenius product with these directions.
    damp::array<Matrix<NX, NX, T>, MaxDof> dV{};
    damp::array<Matrix<NU, NX, T>, MaxDof> dW{};
    if (search) {
        size_t j = 0;
        for (size_t e = 0; e < plan.ndistinct; ++e) {
            size_t pcol = plan.pos_offset[e];
            for (size_t p = 0; p < plan.npos[e]; ++p) {
                for (size_t r = 0; r < NU; ++r) {
                    Matrix<NU, NX, Cplx> Kr = Matrix<NU, NX, Cplx>::zeros();
                    Kr(r, pcol) = Cplx{T{1}, T{0}};
                    assemble_vw(plan, Kr, dV[j], dW[j]);
                    ++j;
                    if (!plan.is_real[e]) {
                        Matrix<NU, NX, Cplx> Ki = Matrix<NU, NX, Cplx>::zeros();
                        Ki(r, pcol) = Cplx{T{0}, T{1}};
                        assemble_vw(plan, Ki, dV[j], dW[j]);
                        ++j;
                    }
                }
                ++pcol;
            }
        }
    }

    // Assemble V, W and form P = V⁻¹ at θ (the inverse is reused for the
    // objective, analytic gradient, and κ_F metric). Returns false when V is
    // singular. Split from the objective/gradient so the gradient can reuse the
    // line search's accepted inverse instead of recomputing it.
    const auto assemble_and_invert = [&](const damp::array<T, MaxDof>& th, Matrix<NX, NX, T>& V,
                                         Matrix<NU, NX, T>& W, Matrix<NX, NX, T>& P) -> bool {
        assemble_vw(plan, theta_to_k(th), V, W);
        // Inverse itself is the deliverable (reused many times); form via solve(V, I).
        const auto P_opt = mat::solve(V, Matrix<NX, NX, T>::identity());
        if (!P_opt) {
            return false;
        }
        P = P_opt.value();
        return true;
    };

    // Objective value from a solved (V, W, P = V⁻¹). Squared Frobenius (sqrt-free).
    const auto objective_value = [&](const Matrix<NX, NX, T>& V, const Matrix<NU, NX, T>& W,
                                     const Matrix<NX, NX, T>& P) -> T {
        const Matrix<NU, NX, T> F = Matrix<NU, NX, T>(W * P);
        const T                 gain2 = fro2(F);
        const T                 ga = T{1} - alpha;
        if (objective == JordanObjective::ConditionNumber) {
            return (alpha * (fro2(V) + fro2(P))) + (ga * gain2);
        }
        const Matrix<NX, NX, T> M = Matrix<NX, NX, T>(A + (B * F));
        T                       dep2 = fro2(M) - sum_lambda_sq;
        if (dep2 < T{0}) {
            dep2 = T{0};
        }
        return (alpha * dep2) + (ga * gain2);
    };

    // Analytic gradient from a solved (V, W, P) via the adjoint matrices
    // GV = ∂f/∂V, GW = ∂f/∂W; gradⱼ = ⟨GV, ∂V/∂θⱼ⟩ + ⟨GW, ∂W/∂θⱼ⟩.
    const auto gradient_at = [&](const Matrix<NX, NX, T>& V, const Matrix<NU, NX, T>& W,
                                 const Matrix<NX, NX, T>& P, damp::array<T, MaxDof>& grad_out) {
        const Matrix<NX, NX, T> Pt = Matrix<NX, NX, T>(P.transpose());
        const Matrix<NU, NX, T> F = Matrix<NU, NX, T>(W * P); // paper's F = −K
        const Matrix<NX, NU, T> Ft = Matrix<NX, NU, T>(F.transpose());
        const T                 ga = T{1} - alpha;
        Matrix<NX, NX, T>       GV;
        Matrix<NU, NX, T>       GW;
        if (objective == JordanObjective::ConditionNumber) {
            const Matrix<NX, NX, T> Gb = Matrix<NX, NX, T>(Matrix<NX, NX, T>(Pt * P) * Pt);    // Pᵀ P Pᵀ
            const Matrix<NX, NX, T> FtFPt = Matrix<NX, NX, T>(Matrix<NX, NX, T>(Ft * F) * Pt); // Fᵀ F Pᵀ
            GV = Matrix<NX, NX, T>(((T{2} * alpha) * V) - ((T{2} * alpha) * Gb) - ((T{2} * ga) * FtFPt));
            GW = Matrix<NU, NX, T>((T{2} * ga) * Matrix<NU, NX, T>(F * Pt));
        } else {
            const Matrix<NX, NX, T> M = Matrix<NX, NX, T>(A + (B * F));                                                      // closed loop A + B·F
            const Matrix<NU, NX, T> BtMPt = Matrix<NU, NX, T>(Matrix<NU, NX, T>(Matrix<NU, NX, T>(B.transpose()) * M) * Pt); // Bᵀ M Pᵀ
            const Matrix<NX, NX, T> FtBtMPt = Matrix<NX, NX, T>(Ft * BtMPt);                                                 // Fᵀ Bᵀ M Pᵀ
            const Matrix<NX, NX, T> FtFPt = Matrix<NX, NX, T>(Matrix<NX, NX, T>(Ft * F) * Pt);                               // Fᵀ F Pᵀ
            GV = Matrix<NX, NX, T>(((T{-2} * alpha) * FtBtMPt) - ((T{2} * ga) * FtFPt));
            GW = Matrix<NU, NX, T>(((T{2} * alpha) * BtMPt) + ((T{2} * ga) * Matrix<NU, NX, T>(F * Pt)));
        }
        for (size_t j = 0; j < dof; ++j) {
            grad_out[j] = GV.dot(dV[j]) + GW.dot(dW[j]);
        }
    };

    // Start from the canonical parameter.
    damp::array<T, MaxDof> theta{};
    {
        const auto Kc = canonical_kparams(plan);
        size_t     idx = 0;
        for (size_t e = 0; e < plan.ndistinct; ++e) {
            size_t pcol = plan.pos_offset[e];
            for (size_t p = 0; p < plan.npos[e]; ++p) {
                for (size_t r = 0; r < NU; ++r) {
                    theta[idx] = Kc(r, pcol).real();
                    idx += plan.is_real[e] ? size_t{1} : size_t{2};
                    if (!plan.is_real[e]) {
                        theta[idx - 1] = Kc(r, pcol).imag();
                    }
                }
                ++pcol;
            }
        }
    }

    damp::array<T, MaxDof> grad{};
    Matrix<NX, NX, T>      Vcur;
    Matrix<NU, NX, T>      Wcur;
    Matrix<NX, NX, T>      Pcur;
    if (!assemble_and_invert(theta, Vcur, Wcur, Pcur)) {
        return damp::nullopt; // canonical parameter already gives a singular V
    }
    T f0 = objective_value(Vcur, Wcur, Pcur);
    if (search) {
        gradient_at(Vcur, Wcur, Pcur, grad);
    }

    // BFGS with an Armijo backtracking line search warm-started at the unit
    // (quasi-Newton) step. H approximates the inverse Hessian. The accepted
    // line-search inverse is reused for the gradient (one V⁻¹ per step).
    damp::array<damp::array<T, MaxDof>, MaxDof> H{};
    for (size_t i = 0; i < dof; ++i) {
        H[i][i] = T{1};
    }
    const T c1 = static_cast<T>(1e-4);
    const T ftol = static_cast<T>(1e-12);
    const T gtol = damp::sqrt(eps<T>()); // relative to 1 + |f|
    bool    converged = !search;
    size_t  iter = 0;
    for (; iter < max_iter; ++iter) {
        T gnorm2 = T{0};
        for (size_t i = 0; i < dof; ++i) {
            gnorm2 += grad[i] * grad[i];
        }
        if (damp::sqrt(gnorm2) < (gtol * (T{1} + damp::abs(f0)))) {
            converged = true;
            break;
        }

        // Search direction p = −H·g; fall back to steepest descent if not a descent.
        damp::array<T, MaxDof> p{};
        T                      gp = T{0};
        for (size_t i = 0; i < dof; ++i) {
            T acc = T{0};
            for (size_t k = 0; k < dof; ++k) {
                acc += H[i][k] * grad[k];
            }
            p[i] = -acc;
            gp += grad[i] * p[i];
        }
        if (gp >= T{0}) {
            for (size_t i = 0; i < dof; ++i) {
                for (size_t k = 0; k < dof; ++k) {
                    H[i][k] = (i == k) ? T{1} : T{0};
                }
                p[i] = -grad[i];
            }
            gp = -gnorm2;
        }

        // Line search; cache the accepted (V, W, P) for the gradient.
        T                      step = T{1};
        bool                   ok = false;
        damp::array<T, MaxDof> theta_new{};
        T                      f_new = f0;
        Matrix<NX, NX, T>      Va;
        Matrix<NU, NX, T>      Wa;
        Matrix<NX, NX, T>      Pa;
        for (size_t ls = 0; ls < 40; ++ls) {
            for (size_t i = 0; i < dof; ++i) {
                theta_new[i] = theta[i] + (step * p[i]);
            }
            Matrix<NX, NX, T> Vt;
            Matrix<NU, NX, T> Wt;
            Matrix<NX, NX, T> Pt2;
            if (assemble_and_invert(theta_new, Vt, Wt, Pt2)) {
                const T ft = objective_value(Vt, Wt, Pt2);
                if (ft <= (f0 + (c1 * step * gp))) {
                    Va = Vt;
                    Wa = Wt;
                    Pa = Pt2;
                    f_new = ft;
                    ok = true;
                    break;
                }
            }
            step *= static_cast<T>(0.5);
        }
        if (!ok) {
            break; // line search stalled before the gradient test passed: not converged
        }

        // Gradient at the accepted point reuses the cached inverse (no new solve).
        damp::array<T, MaxDof> grad_new{};
        gradient_at(Va, Wa, Pa, grad_new);

        // BFGS inverse-Hessian update.
        damp::array<T, MaxDof> sv{};
        damp::array<T, MaxDof> yv{};
        T                      sy = T{0};
        T                      ss = T{0};
        T                      yy = T{0};
        for (size_t i = 0; i < dof; ++i) {
            sv[i] = theta_new[i] - theta[i];
            yv[i] = grad_new[i] - grad[i];
            sy += sv[i] * yv[i];
            ss += sv[i] * sv[i];
            yy += yv[i] * yv[i];
        }
        if (sy > (eps<T>() * damp::sqrt(ss * yy))) { // curvature condition, relative to |s||y|
            damp::array<T, MaxDof> Hy{};
            T                      yHy = T{0};
            for (size_t i = 0; i < dof; ++i) {
                T acc = T{0};
                for (size_t k = 0; k < dof; ++k) {
                    acc += H[i][k] * yv[k];
                }
                Hy[i] = acc;
                yHy += yv[i] * Hy[i];
            }
            const T rho = T{1} / sy;
            const T coef = rho * (T{1} + (rho * yHy));
            for (size_t i = 0; i < dof; ++i) {
                for (size_t k = 0; k < dof; ++k) {
                    H[i][k] += (coef * sv[i] * sv[k]) - (rho * ((Hy[i] * sv[k]) + (sv[i] * Hy[k])));
                }
            }
        }

        const T f_prev = f0;
        theta = theta_new;
        grad = grad_new;
        f0 = f_new;
        if ((f_prev - f0) <= (ftol * (damp::abs(f_prev) + T{1}))) {
            converged = true;
            break;
        }
    }

    // Assemble the final gain and metrics (reported as true norms).
    // P = V⁻¹ is reused for F = W·P and κ_F(V) = ‖V‖_F·‖P‖_F.
    Matrix<NX, NX, T> V;
    Matrix<NU, NX, T> W;
    Matrix<NX, NX, T> P;
    if (!assemble_and_invert(theta, V, W, P)) {
        return damp::nullopt;
    }
    const Matrix<NU, NX, T> F = Matrix<NU, NX, T>(W * P);
    const Matrix<NX, NX, T> M = Matrix<NX, NX, T>(A + (B * F));
    T                       dep2 = (M.norm() * M.norm()) - sum_lambda_sq;
    if (dep2 < T{0}) {
        dep2 = T{0};
    }
    JordanPlaceResult<NX, NU, T> result{};
    result.K = Matrix<NU, NX, T>(T{-1} * F);
    result.cond_fro = V.norm() * P.norm();
    result.gain_fro = F.norm();
    result.departure_fro = damp::sqrt(dep2);
    result.iterations = iter;
    result.converged = converged;
    return result;
}

/** @brief The NX eigenvalues of a Jordan structure, each block's eigenvalue repeated by its size (the plan checked the sizes sum to NX). */
template<size_t NX, typename T, size_t NB>
[[nodiscard]] constexpr damp::array<damp::complex<T>, NX> block_spectrum(const damp::array<JordanBlock<T>, NB>& blocks) {
    damp::array<damp::complex<T>, NX> poles{};
    size_t                            k = 0;
    for (size_t b = 0; b < NB; ++b) {
        for (size_t l = 0; l < blocks[b].size && k < NX; ++l) {
            poles[k] = blocks[b].eigenvalue;
            ++k;
        }
    }
    return poles;
}

/** @brief place_jordan on a plant whose domain is already the law's: refuses inadmissible options and structures. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NB>
[[nodiscard]] constexpr JordanPlaceResult<NX, NU, T> place_jordan_on(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const damp::array<JordanBlock<T>, NB>&   blocks,
    const JordanOptions<T>&                  opts
) {
    using Result = JordanPlaceResult<NX, NU, T>;
    if (!(opts.alpha >= T{0}) || !(opts.alpha <= T{1})) {
        return ::damp::detail::refused<Result>(DesignStatus::bad_weights);
    }
    for (size_t b = 0; b < NB; ++b) {
        if (!damp::isfinite(blocks[b].eigenvalue.real()) || !damp::isfinite(blocks[b].eigenvalue.imag())) {
            return ::damp::detail::refused<Result>(DesignStatus::invalid_parameter);
        }
    }
    if (mat::rank(sys.B) < NU) {
        return ::damp::detail::refused<Result>(DesignStatus::singular);
    }
    const auto placed = jordan_place(sys.A, sys.B, blocks, opts);
    if (!placed || !placement_matches(sys.A, sys.B, placed->K, block_spectrum<NX>(blocks))) {
        return ::damp::detail::refused<Result>(DesignStatus::not_placeable);
    }
    Result out = *placed;
    out.Ts = sys.Ts;
    out.e = stability::closed_loop_poles(sys.A, sys.B, out.K);
    ::damp::detail::set_status(out, DesignStatus::ok);
    return out;
}

} // namespace detail

/**
 * @ingroup pole_placement
 * @brief Exact pole placement with an arbitrary Jordan structure
 *        (Schmid–Ntogramatzidis–Nguyen–Pandey / Klein–Moore parametric form)
 *
 * Unlike place(), which assigns a non-defective spectrum (each eigenvalue with
 * independent eigenvectors), this assigns *any* admissible eigenstructure: any
 * eigenvalues with any algebraic multiplicities and any Jordan mini-block
 * orders, including fully defective blocks and closed-loop poles that coincide
 * with open-loop ones. It computes a real gain K such that A − B·K has exactly
 * the requested Jordan structure (same A − B·K convention and domain rules as place()).
 *
 * The method builds each Jordan chain from the kernel and Moore–Penrose
 * pseudoinverse of the matrix pencil S(λ) = [A − λI | B]: the chain head is a
 * kernel vector (a closed-loop eigenvector) and each successor solves
 * S(λ)·h = v_prev for the generalized eigenvector, which the pseudoinverse does
 * directly since S(λ) has full row rank for a reachable (A, B). Stacking the
 * state parts gives the eigenvector matrix V and the input parts W; the gain is
 * K = −W·V⁻¹. Conjugate eigenpairs are realified (Re/Im columns) so K is real.
 *
 * Every K in the family places the structure exactly (Theorem 2.1), so @p opts only
 * trades robustness against gain. With JordanObjective::Canonical (the default) the
 * canonical minimum-norm parameter is used. Otherwise the objective
 * f(K) = α·R(K) + (1 − α)·‖K‖²_F is minimized, with R selected by the objective:
 *  - ConditionNumber (Method 1): R = ‖V‖²_F + ‖V⁻¹‖²_F, the Byers–Nash proxy for
 *    the Frobenius condition number κ_F(V) (the eigenvalue sensitivity).
 *  - DepartureFromNormality (Method 2): R = δ²_F(A − B·K) = ‖A − B·K‖²_F − Σ|λᵢ|²,
 *    a closed form since the closed-loop eigenvalues are exactly the targets.
 *
 * α = 1 is the pure robust problem (REPP), α = 0 the pure minimum-gain problem
 * (MGEPP). The nonconvex objective is minimized by BFGS with an Armijo line search
 * from the canonical parameter; the result is a local minimum dependent on that
 * start. The eigenvector/input matrices are linear in K, so the gradient is analytic
 * (one V⁻¹ per iteration) and the kernel/pseudoinverse factors are computed once.
 * Constant evaluation of larger systems may need `-fconstexpr-ops-limit` raised.
 *
 * Status: invalid_parameter if an eigenvalue is non-finite; not_placeable if the structure
 * is inadmissible (block sizes do not sum to NX, an eigenvalue is asked for more mini-blocks
 * than NU, a conjugate pair's block orders differ, (A, B) is not reachable at some λ, or the
 * parameter yields a singular eigenvector matrix) or the achieved eig(A − BK) misses the
 * requested eigenvalues; singular if B is rank-deficient, bad_weights if α is outside [0, 1].
 *
 * @tparam NB Number of Jordan blocks supplied
 * @param sys    Plant; only A and B are used; (A, B) reachable
 * @param blocks Desired Jordan structure in the domain of @p sys; block sizes sum to
 *               NX, complex eigenvalues given as equal-size conjugate pairs
 * @param opts   Search objective, weight α and iteration cap
 * @return JordanPlaceResult with K, e, Ts and the achieved metrics
 *
 * @see "A unified method for optimal arbitrary pole placement" (Schmid, Ntogramatzidis,
 *      Nguyen, Pandey, 2014), Automatica 50(8), https://doi.org/10.1016/j.automatica.2014.05.020
 * @see "Eigenvalue-generalized eigenvector assignment with state feedback" (Klein, Moore,
 *      1977), IEEE TAC 22(1)
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NB>
[[nodiscard]] constexpr JordanPlaceResult<NX, NU, T> place_jordan(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const damp::array<JordanBlock<T>, NB>&   blocks,
    const JordanOptions<T>&                  opts = {}
) {
    return detail::place_jordan_on(sys, blocks, opts);
}

/**
 * @ingroup pole_placement
 * @brief place_jordan for a continuous plant sampled at @p Ts: s-plane eigenvalues, discrete gain
 *
 * ZOH-discretizes @p sys and maps each block's eigenvalue to z = e^{λ·Ts}; block sizes are
 * kept. Refuses a discrete @p sys or a Ts that is not positive and finite with
 * DesignStatus::invalid_sample_time, as place(sys, poles, Ts) does.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NB>
[[nodiscard]] constexpr JordanPlaceResult<NX, NU, T> place_jordan(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const damp::array<JordanBlock<T>, NB>&   blocks,
    T                                        Ts,
    const JordanOptions<T>&                  opts = {}
) {
    StateSpace<NX, NU, NY, T, NW, NV> sampled{};
    const DesignStatus                status = detail::sample_plant(sys, Ts, sampled);
    if (status != DesignStatus::ok) {
        return ::damp::detail::refused<JordanPlaceResult<NX, NU, T>>(status);
    }
    damp::array<JordanBlock<T>, NB> mapped = blocks;
    for (size_t b = 0; b < NB; ++b) {
        mapped[b].eigenvalue = detail::pole_s_to_z(blocks[b].eigenvalue, Ts);
    }
    return detail::place_jordan_on(sampled, mapped, opts);
}

} // namespace design

} // namespace damp
