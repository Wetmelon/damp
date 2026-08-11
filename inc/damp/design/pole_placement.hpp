// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file pole_placement.hpp
 * @brief Pole placement (place/acker)
 *
 * @defgroup pole_placement Robust pole placement (place)
 * @brief Multi-input place (Kautsky–Nichols–Van Dooren) plus Ackermann and observer duals
 *
 * For a controllable pair (A, B) and a desired closed-loop spectrum, find a
 * state-feedback gain K such that the eigenvalues of (A − B·K) equal the desired
 * poles. A multi-input system has extra freedom beyond merely assigning the
 * spectrum; the Kautsky–Nichols–Van Dooren (KNV) method spends that freedom to
 * make the eigenvector matrix X as well-conditioned as possible, which minimizes
 * the sensitivity of the placed poles to perturbations in A, B, and K — the
 * robustness that distinguishes `place` from the single-input Ackermann formula
 * (`matlab::acker`).
 *
 * Method 0 (orthogonal) is implemented: each eigenvector is repeatedly re-chosen
 * within its admissible subspace to be maximally orthogonal to the span of the
 * others, driving κ(X) down. The spectrum is assigned *exactly* regardless of
 * how far the conditioning sweep runs (placement correctness depends only on X
 * being invertible); the sweep only improves robustness.
 *
 * ## Overload contract (state feedback)
 *
 * | Call | Plant | Poles | Action |
 * |------|-------|-------|--------|
 * | `place(A, B, poles)` | continuous \(A,B\) | s-plane | continuous \(K\) for \(A-BK\) |
 * | `place(A, B, poles, Ts)` | continuous \(A,B\) | s-plane | ZOH plant + \(z=e^{sT_s}\) + discrete place → discrete \(K\) |
 * | `place(sys, poles)` | from `sys` | domain of `sys` | `place(sys.A, sys.B, poles)` — no pole mapping |
 * | `place_discrete(Ad, Bd, Zi)` | already discrete | z-plane | discrete place |
 * | `place_discrete(Ad, Bd, poles, Ts)` | already discrete | s-plane | map \(z=e^{sT_s}\) only (no second c2d), then place |
 *
 * Rule: an extra `Ts` argument always means poles are continuous (s-plane).
 * `place(sys, poles)` never maps poles — domain follows `sys.Ts` only.
 *
 * Luenberger duals: `place_observer` / `place_observer_discrete` with the same shape
 * (L for A-LC via the transpose dual of place).
 *
 * @see J. Kautsky, N. K. Nichols, P. Van Dooren, "Robust pole assignment in
 *      linear state feedback," Int. J. Control 41(5), 1985,
 *      https://doi.org/10.1080/00207178508933420
 * @see matlab::acker for the single-input Ackermann path
 */

#include <cstddef>
#include <cstdint>
#include <limits>

#include "damp/backend.hpp"
#include "damp/design/stability.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/decomposition.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace design {

/**
 * @brief Robust multi-input pole placement (Kautsky–Nichols–Van Dooren, real poles).
 *
 * Computes the state-feedback gain K placing the eigenvalues of (A − B·K) at the
 * requested real poles, using KNV Method 0 to minimize the conditioning of the
 * eigenvector basis.
 *
 * @tparam NX Number of states
 * @tparam NU Number of inputs (NU ≤ NX)
 * @param A      State matrix (NX×NX)
 * @param B      Input matrix (NX×NU), assumed full column rank (controllable)
 * @param poles  Desired closed-loop eigenvalues (real)
 * @return Gain K (NU×NX), or damp::nullopt if B is rank-deficient or the assigned
 *         eigenvectors are linearly dependent (e.g. a pole repeated more than NU
 *         times — not assignable with independent eigenvectors).
 *
 * @note Compare with MATLAB®'s K = place(A, B, poles).
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NU, NX, T>> place(
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
        const T          tol = T{64} * std::numeric_limits<T>::epsilon();
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
 * @brief Robust pole placement with complex-conjugate poles (KNV).
 *
 * Generalizes place() to complex spectra. Each conjugate pair (σ ± jω) is
 * assigned in real arithmetic via a real eigenvector pair (Re v, Im v) and a
 * 2×2 real block [[σ, ω], [−ω, σ]] in the closed-loop matrix, so K stays real.
 * Complex poles must be supplied as adjacent conjugate pairs.
 *
 * An all-real spectrum forwards to the real-pole overload (which additionally
 * runs the Method-0 conditioning sweep). The complex path assigns each
 * eigenvector from the orthonormal admissible basis without the extra sweep —
 * the spectrum is still placed exactly; the conditioning is good but not sweep-
 * optimized.
 *
 * @param A     Open-loop state matrix (NX × NX)
 * @param B     Input matrix (NX × NU)
 * @param poles Desired closed-loop eigenvalues as adjacent conjugate pairs
 *              (real poles may appear anywhere).
 * @return Gain K (NU×NX), or damp::nullopt if B is rank-deficient, the poles are
 *         not in valid conjugate pairs, or the eigenvectors are dependent.
 *
 * @note Compare with MATLAB®'s K = place(A, B, poles) for complex spectra.
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NU, NX, T>> place(
    const Matrix<NX, NX, T>&                 A,
    const Matrix<NX, NU, T>&                 B,
    const damp::array<damp::complex<T>, NX>& poles
) {
    static_assert(NU <= NX, "place requires NU <= NX (no more inputs than states)");
    constexpr T tol_im = static_cast<T>(1e-9);

    // All-real spectrum → reuse the conditioned real-pole routine.
    bool all_real = true;
    for (size_t i = 0; i < NX; ++i) {
        if (damp::abs(poles[i].imag()) >= tol_im) {
            all_real = false;
            break;
        }
    }
    if (all_real) {
        damp::array<T, NX> rp{};
        for (size_t i = 0; i < NX; ++i) {
            rp[i] = poles[i].real();
        }
        return place(A, B, rp);
    }

    // Build the real block-diagonal Λ from the (paired) spectrum: a real pole is
    // a 1×1 entry, a conjugate pair (σ ± jω) a 2×2 block [[σ, ω], [−ω, σ]].
    // Returns false if the poles are not validly paired.
    Matrix<NX, NX, T> Lambda = Matrix<NX, NX, T>::zeros();
    {
        size_t j = 0;
        while (j < NX) {
            if (damp::abs(poles[j].imag()) < tol_im) {
                Lambda(j, j) = poles[j].real();
                ++j;
            } else {
                if (j + 1 >= NX) {
                    return damp::nullopt; // dangling complex pole
                }
                const T sg = poles[j].real();
                const T w = poles[j].imag();
                if (damp::abs(poles[j + 1].real() - sg) > static_cast<T>(1e-7) || damp::abs(poles[j + 1].imag() + w) > static_cast<T>(1e-7)) {
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
            if (damp::abs(poles[j].imag()) < tol_im) {
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

// -----------------------------------------------------------------------------
// Domain helpers: continuous poles + Ts → discrete z-poles (and StateSpace place)
// -----------------------------------------------------------------------------

namespace detail {

/// \(z = e^{s T_s}\) for a real pole.
template<typename T>
[[nodiscard]] constexpr T pole_s_to_z(T s, T Ts) {
    return damp::exp(s * Ts);
}

/// Continuous-to-discrete complex pole via exp and sincos.
template<typename T>
[[nodiscard]] constexpr damp::complex<T> pole_s_to_z(damp::complex<T> s, T Ts) {
    const T er = damp::exp(s.real() * Ts);
    const T wt = s.imag() * Ts;
    const auto [sn, cs] = damp::sincos(wt);
    return damp::complex<T>{er * cs, er * sn};
}

template<size_t NX, typename T>
[[nodiscard]] constexpr damp::array<T, NX> poles_s_to_z(const damp::array<T, NX>& poles_s, T Ts) {
    damp::array<T, NX> z{};
    for (size_t i = 0; i < NX; ++i) {
        z[i] = pole_s_to_z(poles_s[i], Ts);
    }
    return z;
}

template<size_t NX, typename T>
[[nodiscard]] constexpr damp::array<damp::complex<T>, NX>
poles_s_to_z(const damp::array<damp::complex<T>, NX>& poles_s, T Ts) {
    damp::array<damp::complex<T>, NX> z{};
    for (size_t i = 0; i < NX; ++i) {
        z[i] = pole_s_to_z(poles_s[i], Ts);
    }
    return z;
}

} // namespace detail

/**
 * @brief Continuous plant + continuous poles + sample period → discrete gain \(K\).
 *
 * ZOH-discretizes \((A,B)\), maps poles \(z_i = e^{p_i T_s}\), then places on the
 * discrete pair. Use this when you think in rad/s but ship \(u[k] = -K x[k]\).
 *
 * @param A      Continuous state matrix
 * @param B      Continuous input matrix
 * @param poles  Desired continuous closed-loop poles (s-plane)
 * @param Ts     Sample period [s] (must be > 0)
 * @return Discrete gain K, or nullopt if Ts <= 0 or placement fails
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NU, NX, T>> place(
    const Matrix<NX, NX, T>&  A,
    const Matrix<NX, NU, T>&  B,
    const damp::array<T, NX>& poles,
    T                         Ts
) {
    if (Ts <= T{0}) {
        return damp::nullopt;
    }
    StateSpace<NX, NU, NX, T> sys_c{
        .A = A,
        .B = B,
        .C = Matrix<NX, NX, T>::identity(),
        .Ts = T{0},
    };
    const auto sys_d = *discretize(sys_c, Ts, DiscretizationMethod::ZOH);
    return place(sys_d.A, sys_d.B, detail::poles_s_to_z(poles, Ts));
}

/**
 * @brief Continuous plant + complex continuous poles + \(T_s\) → discrete \(K\).
 * @see place(A, B, poles, Ts) for the real-pole form
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NU, NX, T>> place(
    const Matrix<NX, NX, T>&                 A,
    const Matrix<NX, NU, T>&                 B,
    const damp::array<damp::complex<T>, NX>& poles,
    T                                        Ts
) {
    if (Ts <= T{0}) {
        return damp::nullopt;
    }
    StateSpace<NX, NU, NX, T> sys_c{
        .A = A,
        .B = B,
        .C = Matrix<NX, NX, T>::identity(),
        .Ts = T{0},
    };
    const auto sys_d = *discretize(sys_c, Ts, DiscretizationMethod::ZOH);
    return place(sys_d.A, sys_d.B, detail::poles_s_to_z(poles, Ts));
}

/**
 * @brief Place poles of a StateSpace plant (domain follows `sys.Ts`).
 *
 * - Continuous (`Ts == 0`): poles are s-plane; returns continuous \(K\).
 * - Discrete (`Ts > 0`): poles are z-plane; returns discrete \(K\).
 *
 * Does not map poles or re-discretize. For continuous physics + digital gain
 * with s-plane poles, call `place(sys.A, sys.B, poles_s, Ts)` instead.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr damp::optional<Matrix<NU, NX, T>> place(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const damp::array<T, NX>&                poles
) {
    return place(sys.A, sys.B, poles);
}

/**
 * @brief Place complex poles of a StateSpace plant (domain follows `sys.Ts`).
 * @see place(sys, poles) for the real-pole form
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr damp::optional<Matrix<NU, NX, T>> place(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const damp::array<damp::complex<T>, NX>& poles
) {
    return place(sys.A, sys.B, poles);
}

/**
 * @brief Discrete plant + z-plane poles → discrete \(K\).
 *
 * Plant is already discrete; poles are already in z. Same math as
 * `place(Ad, Bd, Zi)` — the name documents the domain.
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NU, NX, T>> place_discrete(
    const Matrix<NX, NX, T>&  Ad,
    const Matrix<NX, NU, T>&  Bd,
    const damp::array<T, NX>& Zi
) {
    return place(Ad, Bd, Zi);
}

/**
 * @brief Discrete plant + complex z-plane poles → discrete \(K\).
 * @see place_discrete(Ad, Bd, Zi)
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NU, NX, T>> place_discrete(
    const Matrix<NX, NX, T>&                 Ad,
    const Matrix<NX, NU, T>&                 Bd,
    const damp::array<damp::complex<T>, NX>& Zi
) {
    return place(Ad, Bd, Zi);
}

/**
 * @brief Discrete plant + continuous poles + \(T_s\) → discrete \(K\).
 *
 * Maps \(z_i = e^{p_i T_s}\) only; does not discretize the plant (caller
 * already has \(A_d, B_d\)). Use when the model is discrete but you still pick
 * poles in rad/s.
 *
 * @param Ad     Discrete state matrix
 * @param Bd     Discrete input matrix
 * @param poles  Desired continuous poles (s-plane)
 * @param Ts     Sample period used for the pole map (must be > 0)
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NU, NX, T>> place_discrete(
    const Matrix<NX, NX, T>&  Ad,
    const Matrix<NX, NU, T>&  Bd,
    const damp::array<T, NX>& poles,
    T                         Ts
) {
    if (Ts <= T{0}) {
        return damp::nullopt;
    }
    return place(Ad, Bd, detail::poles_s_to_z(poles, Ts));
}

/**
 * @brief Discrete plant + complex continuous poles + \(T_s\) → discrete \(K\).
 * @see place_discrete(Ad, Bd, poles, Ts)
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NU, NX, T>> place_discrete(
    const Matrix<NX, NX, T>&                 Ad,
    const Matrix<NX, NU, T>&                 Bd,
    const damp::array<damp::complex<T>, NX>& poles,
    T                                        Ts
) {
    if (Ts <= T{0}) {
        return damp::nullopt;
    }
    return place(Ad, Bd, detail::poles_s_to_z(poles, Ts));
}

// -----------------------------------------------------------------------------
// Luenberger duals: L for (A − L C), same domain rules as place
// -----------------------------------------------------------------------------

/**
 * @brief Continuous observer gain: place eigenvalues of \(A - LC\) (s-plane poles).
 *
 * Dual of place: L is recovered from place(A^T, C^T, p)^T.
 *
 * @param A      Continuous state matrix
 * @param C      Output matrix (NY × NX), NY ≤ NX
 * @param poles  Desired continuous error poles
 * @return Luenberger gain \(L\) (NX × NY), or nullopt if unplaceable
 */
template<size_t NX, size_t NY, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NY, T>> place_observer(
    const Matrix<NX, NX, T>&  A,
    const Matrix<NY, NX, T>&  C,
    const damp::array<T, NX>& poles
) {
    static_assert(NY <= NX, "place_observer requires NY <= NX");
    const auto Lt = place<NX, NY, T>(A.transpose(), C.transpose(), poles);
    if (!Lt) {
        return damp::nullopt;
    }
    return Lt->transpose();
}

/**
 * @brief Continuous observer gain with complex s-plane poles.
 * @see place_observer(A, C, poles)
 */
template<size_t NX, size_t NY, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NY, T>> place_observer(
    const Matrix<NX, NX, T>&                 A,
    const Matrix<NY, NX, T>&                 C,
    const damp::array<damp::complex<T>, NX>& poles
) {
    static_assert(NY <= NX, "place_observer requires NY <= NX");
    const auto Lt = place<NX, NY, T>(A.transpose(), C.transpose(), poles);
    if (!Lt) {
        return damp::nullopt;
    }
    return Lt->transpose();
}

/**
 * @brief Continuous plant + s-plane poles + \(T_s\) → discrete observer gain \(L\).
 *
 * ZOH-discretizes the plant (dummy \(B\); \(C_d = C\)), maps poles \(z = e^{s T_s}\),
 * then places on the discrete dual. Does not ZOH the dual pair (A^T, C^T)
 * — that would not equal (Ad^T, Cd^T).
 */
template<size_t NX, size_t NY, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NY, T>> place_observer(
    const Matrix<NX, NX, T>&  A,
    const Matrix<NY, NX, T>&  C,
    const damp::array<T, NX>& poles,
    T                         Ts
) {
    static_assert(NY <= NX, "place_observer requires NY <= NX");
    if (Ts <= T{0}) {
        return damp::nullopt;
    }
    StateSpace<NX, 1, NY, T> sys_c{
        .A = A,
        .B = Matrix<NX, 1, T>{},
        .C = C,
        .Ts = T{0},
    };
    const auto sys_d = *discretize(sys_c, Ts, DiscretizationMethod::ZOH);
    return place_observer(sys_d.A, sys_d.C, detail::poles_s_to_z(poles, Ts));
}

/**
 * @brief Continuous plant + complex s-plane poles + \(T_s\) → discrete \(L\).
 * @see place_observer(A, C, poles, Ts)
 */
template<size_t NX, size_t NY, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NY, T>> place_observer(
    const Matrix<NX, NX, T>&                 A,
    const Matrix<NY, NX, T>&                 C,
    const damp::array<damp::complex<T>, NX>& poles,
    T                                        Ts
) {
    static_assert(NY <= NX, "place_observer requires NY <= NX");
    if (Ts <= T{0}) {
        return damp::nullopt;
    }
    StateSpace<NX, 1, NY, T> sys_c{
        .A = A,
        .B = Matrix<NX, 1, T>{},
        .C = C,
        .Ts = T{0},
    };
    const auto sys_d = *discretize(sys_c, Ts, DiscretizationMethod::ZOH);
    return place_observer(sys_d.A, sys_d.C, detail::poles_s_to_z(poles, Ts));
}

/**
 * @brief Luenberger place on a StateSpace (domain follows `sys.Ts`; no pole mapping).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NY, T>> place_observer(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const damp::array<T, NX>&                poles
) {
    return place_observer(sys.A, sys.C, poles);
}

/**
 * @brief Luenberger place on a StateSpace with complex poles (domain follows `sys.Ts`).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NY, T>> place_observer(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const damp::array<damp::complex<T>, NX>& poles
) {
    return place_observer(sys.A, sys.C, poles);
}

/**
 * @brief Discrete plant + z-plane poles → discrete observer gain \(L\).
 */
template<size_t NX, size_t NY, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NY, T>> place_observer_discrete(
    const Matrix<NX, NX, T>&  Ad,
    const Matrix<NY, NX, T>&  Cd,
    const damp::array<T, NX>& Zi
) {
    return place_observer(Ad, Cd, Zi);
}

/**
 * @brief Discrete plant + complex z-plane poles → discrete \(L\).
 */
template<size_t NX, size_t NY, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NY, T>> place_observer_discrete(
    const Matrix<NX, NX, T>&                 Ad,
    const Matrix<NY, NX, T>&                 Cd,
    const damp::array<damp::complex<T>, NX>& Zi
) {
    return place_observer(Ad, Cd, Zi);
}

/**
 * @brief Discrete plant + continuous poles + \(T_s\) → discrete \(L\) (map poles only).
 */
template<size_t NX, size_t NY, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NY, T>> place_observer_discrete(
    const Matrix<NX, NX, T>&  Ad,
    const Matrix<NY, NX, T>&  Cd,
    const damp::array<T, NX>& poles,
    T                         Ts
) {
    if (Ts <= T{0}) {
        return damp::nullopt;
    }
    return place_observer(Ad, Cd, detail::poles_s_to_z(poles, Ts));
}

/**
 * @brief Discrete plant + complex continuous poles + \(T_s\) → discrete \(L\).
 */
template<size_t NX, size_t NY, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NY, T>> place_observer_discrete(
    const Matrix<NX, NX, T>&                 Ad,
    const Matrix<NY, NX, T>&                 Cd,
    const damp::array<damp::complex<T>, NX>& poles,
    T                                        Ts
) {
    if (Ts <= T{0}) {
        return damp::nullopt;
    }
    return place_observer(Ad, Cd, detail::poles_s_to_z(poles, Ts));
}

/**
 * @brief One Jordan mini-block of a desired closed-loop spectrum.
 *
 * A block contributes `size` consecutive eigenvalues equal to `eigenvalue`,
 * coupled into a single Jordan chain of that order. `size == 1` is an ordinary
 * (semisimple) eigenvalue; `size > 1` is a defective block requiring generalized
 * eigenvectors. Complex eigenvalues must be supplied as conjugate-pair blocks of
 * equal size (e.g. one block at σ+jω and one at σ−jω).
 */
template<typename T = double>
struct JordanBlock {
    damp::complex<T> eigenvalue; ///< The eigenvalue λ this block places.
    size_t           size;       ///< Mini-block order pᵢₖ (chain length).
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

/// Build the K-independent plan, or nullopt if the requested structure is inadmissible.
template<size_t NX, size_t NU, size_t NB, typename T>
[[nodiscard]] constexpr damp::optional<JordanPlan<NX, NU, NB, T>> prepare_jordan_plan(
    const Matrix<NX, NX, T>&               A,
    const Matrix<NX, NU, T>&               B,
    const damp::array<JordanBlock<T>, NB>& blocks
) {
    using Cplx = damp::complex<T>;
    constexpr size_t NS = NX + NU;
    const T          tol_eq = static_cast<T>(1e-7);

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
        if (consumed[bi]) {
            continue;
        }
        const Cplx lam = blocks[bi].eigenvalue;
        const bool is_real = (damp::abs(lam.imag()) <= tol_eq);
        if (!is_real && lam.imag() < T{0}) {
            continue; // negative-imaginary member is realified with its partner
        }

        // Count then gather this eigenvalue's mini-blocks.
        size_t cnt = 0;
        for (size_t bj = 0; bj < NB; ++bj) {
            if (!consumed[bj] && damp::abs(blocks[bj].eigenvalue - lam) <= tol_eq) {
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
            if (!consumed[bj] && damp::abs(blocks[bj].eigenvalue - lam) <= tol_eq) {
                ords[g] = blocks[bj].size;
                mult += blocks[bj].size;
                ++g;
                consumed[bj] = true;
            }
        }
        // A complex eigenvalue needs a matching conjugate of equal multiplicity.
        if (!is_real) {
            size_t conj_mult = 0;
            for (size_t bj = 0; bj < NB; ++bj) {
                if (!consumed[bj] && damp::abs(blocks[bj].eigenvalue - damp::conj(lam)) <= tol_eq) {
                    conj_mult += blocks[bj].size;
                    consumed[bj] = true;
                }
            }
            if (conj_mult != mult) {
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

/// Eigenvector matrix V (state parts, real) and input matrix W (input parts) for a given K.
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

/// The canonical (minimum-norm) parameter: chain k starts from kernel column k.
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

} // namespace detail

/**
 * @ingroup pole_placement
 * @brief Exact pole placement with an arbitrary Jordan structure
 *        (Schmid–Ntogramatzidis–Nguyen–Pandey / Klein–Moore parametric form).
 *
 * Unlike place(), which assigns a non-defective spectrum (each eigenvalue with
 * independent eigenvectors), this assigns *any* admissible eigenstructure: any
 * eigenvalues with any algebraic multiplicities and any Jordan mini-block
 * orders — including fully defective blocks and closed-loop poles that coincide
 * with open-loop ones. It computes a real gain K such that A − B·K has exactly
 * the requested Jordan structure (same A − B·K convention as place()).
 *
 * The method builds each Jordan chain from the kernel and Moore–Penrose
 * pseudoinverse of the matrix pencil S(λ) = [A − λI | B]: the chain head is a
 * kernel vector (a closed-loop eigenvector) and each successor solves
 * S(λ)·h = v_prev for the generalized eigenvector, which the pseudoinverse does
 * directly since S(λ) has full row rank for a reachable (A, B). Stacking the
 * state parts gives the eigenvector matrix V and the input parts W; the gain is
 * K = −W·V⁻¹. Conjugate eigenpairs are realified (Re/Im columns) so K is real.
 *
 * This uses the canonical minimum-norm chain parameter, which places the
 * structure exactly but does not optimize robustness or gain. For the paper's
 * robust/minimum-gain selection of the free parameter, use place_jordan_optimal.
 *
 * @tparam NX Number of states
 * @tparam NU Number of inputs (NU ≤ NX)
 * @tparam NB Number of Jordan blocks supplied
 * @param A      State matrix (NX×NX)
 * @param B      Input matrix (NX×NU), assumed full column rank and (A,B) reachable
 * @param blocks Desired Jordan structure; block sizes must sum to NX, complex
 *               eigenvalues given as equal-size conjugate pairs
 * @return Gain K (NU×NX) with A − B·K in the requested Jordan form, or
 *         damp::nullopt if the structure is inadmissible (block sizes do not sum
 *         to NX, an eigenvalue is asked for more mini-blocks than NU, (A,B) is
 *         not reachable at some λ, or the canonical parameter yields a singular
 *         eigenvector matrix).
 *
 * @see R. Schmid, L. Ntogramatzidis, T. Nguyen, A. Pandey, "A unified method for
 *      optimal arbitrary pole placement," Automatica 50(8), 2014,
 *      https://doi.org/10.1016/j.automatica.2014.05.020
 * @see G. Klein, B. C. Moore, "Eigenvalue-generalized eigenvector assignment
 *      with state feedback," IEEE TAC 22(1), 1977.
 * @see place for the non-defective (distinct/semisimple) robust path.
 */
template<size_t NX, size_t NU, size_t NB, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NU, NX, T>> place_jordan(
    const Matrix<NX, NX, T>&               A,
    const Matrix<NX, NU, T>&               B,
    const damp::array<JordanBlock<T>, NB>& blocks
) {
    static_assert(NU <= NX, "place_jordan requires NU <= NX (no more inputs than states)");

    const auto plan_opt = detail::prepare_jordan_plan(A, B, blocks);
    if (!plan_opt) {
        return damp::nullopt;
    }
    const auto&       plan = plan_opt.value();
    const auto        K = detail::canonical_kparams(plan);
    Matrix<NX, NX, T> V;
    Matrix<NU, NX, T> W;
    detail::assemble_vw(plan, K, V, W);

    // K = −W·V⁻¹ ⇔ K·V = −W: solve Vᵀ Kᵀ = −Wᵀ rather than forming V⁻¹.
    const auto Kt_opt = mat::solve(V.transpose(), Matrix<NX, NU, T>(T{-1} * W.transpose()));
    if (!Kt_opt) {
        return damp::nullopt; // canonical parameter gave dependent eigenvectors
    }
    // A − B·K = V·Λ·V⁻¹ ⇒ closed loop has the Jordan form.
    return Matrix<NU, NX, T>(Kt_opt.value().transpose());
}

/// Robustness objective for place_jordan_optimal (the paper's two methods).
enum class JordanObjective : std::uint8_t {
    ConditionNumber,        ///< Method 1: minimize the Frobenius condition number of V.
    DepartureFromNormality, ///< Method 2: minimize the departure from normality of A − B·K.
};

/// Result of optimized arbitrary pole placement (place_jordan_optimal).
template<size_t NU, size_t NX, typename T = double>
struct OptimalJordanPlacement {
    Matrix<NU, NX, T> gain{};           ///< K, with the A − B·K convention.
    T                 cond_fro{};       ///< Achieved κ_F(V) = ‖V‖F · ‖V⁻¹‖F (eigenvalue robustness).
    T                 gain_fro{};       ///< Achieved ‖K‖F (control effort).
    T                 departure_fro{};  ///< Achieved δ_F(A − B·K) (departure from normality).
    size_t            iterations{};     ///< Gradient-descent iterations taken.
    bool              converged{false}; ///< Search reached a stationary point.

    template<typename U>
    [[nodiscard]] constexpr OptimalJordanPlacement<NU, NX, U> as() const {
        return OptimalJordanPlacement<NU, NX, U>{
            gain.template as<U>(),
            static_cast<U>(cond_fro),
            static_cast<U>(gain_fro),
            static_cast<U>(departure_fro),
            iterations,
            converged,
        };
    }
};

/**
 * @ingroup pole_placement
 * @brief Robust / minimum-gain arbitrary pole placement (Schmid et al., Methods 1–2).
 *
 * Places the same arbitrary Jordan structure as place_jordan, but spends the free
 * parameter K of the Klein–Moore parameterization to optimize a weighted blend of
 * eigenvalue robustness and control effort. Every K in the family places the
 * structure *exactly* (Theorem 2.1), so the search only trades robustness against
 * gain — the assigned poles never move.
 *
 * The objective is f(K) = α·R(K) + (1 − α)·‖K‖²_F, with the robustness term R
 * selected by @p objective:
 *  - ConditionNumber (Method 1): R = ‖V‖²_F + ‖V⁻¹‖²_F, the Byers–Nash proxy for
 *    the Frobenius condition number κ_F(V) (the eigenvalue sensitivity).
 *  - DepartureFromNormality (Method 2): R = δ²_F(A − B·K) = ‖A − B·K‖²_F − Σ|λᵢ|²,
 *    a closed form since the closed-loop eigenvalues are exactly the targets.
 *
 * α = 1 is the pure robust problem (REPP), α = 0 the pure minimum-gain problem
 * (MGEPP). The unconstrained nonconvex objective is minimized by BFGS with an
 * Armijo line search from the canonical parameter; as the paper notes, the result
 * is a local minimum dependent on that start. Because the eigenvector/input
 * matrices are *linear* in K, the gradient is analytic (one V⁻¹ per iteration via
 * the adjoint, reused from the line search, not a finite-difference sweep), and
 * the kernel/pseudoinverse factors are precomputed once — so each iteration is
 * cheap. This keeps it `constexpr`-evaluable for small systems; constant
 * evaluation of larger ones may need `-fconstexpr-ops-limit` raised.
 *
 * @tparam NX Number of states
 * @tparam NU Number of inputs (NU ≤ NX)
 * @tparam NB Number of Jordan blocks supplied
 * @param A         State matrix (NX×NX)
 * @param B         Input matrix (NX×NU), full column rank and (A,B) reachable
 * @param blocks    Desired Jordan structure (see place_jordan)
 * @param alpha     Robustness/gain weight in [0,1]; 1 = robust, 0 = minimum gain
 * @param objective Which robustness measure to minimize (Method 1 or 2)
 * @param max_iter  Maximum gradient-descent iterations
 * @return OptimalJordanPlacement with the gain and achieved metrics, or
 *         damp::nullopt if the structure is inadmissible (as in place_jordan).
 *
 * @see place_jordan for the unoptimized (canonical-parameter) placement.
 * @see R. Schmid et al., "A unified method for optimal arbitrary pole placement,"
 *      Automatica 50(8), 2014.
 */
template<size_t NX, size_t NU, size_t NB, typename T = double>
[[nodiscard]] constexpr damp::optional<OptimalJordanPlacement<NU, NX, T>> place_jordan_optimal(
    const Matrix<NX, NX, T>&               A,
    const Matrix<NX, NU, T>&               B,
    const damp::array<JordanBlock<T>, NB>& blocks,
    T                                      alpha = T{1},
    JordanObjective                        objective = JordanObjective::ConditionNumber,
    size_t                                 max_iter = 200
) {
    static_assert(NU <= NX, "place_jordan_optimal requires NU <= NX");
    using Cplx = damp::complex<T>;

    const auto plan_opt = detail::prepare_jordan_plan(A, B, blocks);
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
    {
        size_t j = 0;
        for (size_t e = 0; e < plan.ndistinct; ++e) {
            size_t pcol = plan.pos_offset[e];
            for (size_t p = 0; p < plan.npos[e]; ++p) {
                for (size_t r = 0; r < NU; ++r) {
                    Matrix<NU, NX, Cplx> Kr = Matrix<NU, NX, Cplx>::zeros();
                    Kr(r, pcol) = Cplx{T{1}, T{0}};
                    detail::assemble_vw(plan, Kr, dV[j], dW[j]);
                    ++j;
                    if (!plan.is_real[e]) {
                        Matrix<NU, NX, Cplx> Ki = Matrix<NU, NX, Cplx>::zeros();
                        Ki(r, pcol) = Cplx{T{0}, T{1}};
                        detail::assemble_vw(plan, Ki, dV[j], dW[j]);
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
        detail::assemble_vw(plan, theta_to_k(th), V, W);
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
        const auto Kc = detail::canonical_kparams(plan);
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
    gradient_at(Vcur, Wcur, Pcur, grad);

    // BFGS with an Armijo backtracking line search warm-started at the unit
    // (quasi-Newton) step. H approximates the inverse Hessian. The accepted
    // line-search inverse is reused for the gradient (one V⁻¹ per step).
    damp::array<damp::array<T, MaxDof>, MaxDof> H{};
    for (size_t i = 0; i < dof; ++i) {
        H[i][i] = T{1};
    }
    const T c1 = static_cast<T>(1e-4);
    const T ftol = static_cast<T>(1e-12);
    const T gtol = static_cast<T>(1e-8);
    bool    converged = false;
    size_t  iter = 0;
    for (; iter < max_iter; ++iter) {
        T gnorm2 = T{0};
        for (size_t i = 0; i < dof; ++i) {
            gnorm2 += grad[i] * grad[i];
        }
        if (damp::sqrt(gnorm2) < gtol) {
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
            converged = true; // stalled at a stationary point
            break;
        }

        // Gradient at the accepted point reuses the cached inverse (no new solve).
        damp::array<T, MaxDof> grad_new{};
        gradient_at(Va, Wa, Pa, grad_new);

        // BFGS inverse-Hessian update.
        damp::array<T, MaxDof> sv{};
        damp::array<T, MaxDof> yv{};
        T                      sy = T{0};
        for (size_t i = 0; i < dof; ++i) {
            sv[i] = theta_new[i] - theta[i];
            yv[i] = grad_new[i] - grad[i];
            sy += sv[i] * yv[i];
        }
        if (sy > static_cast<T>(1e-12)) {
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
    OptimalJordanPlacement<NU, NX, T> result;
    result.gain = Matrix<NU, NX, T>(T{-1} * F);
    result.cond_fro = V.norm() * P.norm();
    result.gain_fro = F.norm();
    result.departure_fro = damp::sqrt(dep2);
    result.iterations = iter;
    result.converged = converged;
    return result;
}

/**
 * @brief Single-input pole placement via Ackermann's formula.
 *
 * Computes the state-feedback gain K placing the eigenvalues of (A − B·K) at the
 * requested poles. Unlike place(), Ackermann works through the characteristic
 * polynomial rather than eigenvector assignment, so it assigns repeated/defective
 * spectra fine — but it is single-input only and numerically weaker for large NX.
 * Use place() for multi-input systems or when robustness matters.
 *
 * @tparam NX Number of states
 * @param A      State matrix (NX×NX)
 * @param B      Input vector (NX×1)
 * @param poles  Desired closed-loop eigenvalues (complex; conjugate pairs cancel
 *               to a real characteristic polynomial)
 * @return Gain K (1×NX), or damp::nullopt if (A, B) is uncontrollable.
 *
 * @note Compare with MATLAB®'s K = acker(A, B, p).
 * @see place for the robust multi-input path.
 */
template<size_t NX, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<1, NX, T>> ackermann(
    const Matrix<NX, NX, T>&                 A,
    const Matrix<NX, 1, T>&                  B,
    const damp::array<damp::complex<T>, NX>& poles
) {
    // Controllability matrix Co = [B, AB, …, A^{NX-1}B] (NU = 1 ⇒ NX×NX).
    const Matrix<NX, NX, T> Co = stability::controllability_matrix(A, B);

    // Desired characteristic polynomial φ(s) = Π(s − pᵢ), built in complex so
    // conjugate pairs cancel to real coefficients.
    damp::array<damp::complex<T>, NX + 1> cc{};
    cc[0] = damp::complex<T>{T{1}, T{0}};
    for (size_t i = 0; i < NX; ++i) {
        const damp::complex<T> root = poles[i];
        damp::complex<T>       carry = cc[0];
        cc[0] = damp::complex<T>{T{0}, T{0}} - (root * cc[0]);
        for (size_t j = 1; j <= NX; ++j) {
            const damp::complex<T> next = cc[j];
            cc[j] = carry - (root * cc[j]);
            carry = next;
        }
    }

    // φ(A) = Σ Re(coeffs[k]) · Aᵏ
    Matrix<NX, NX, T> phi_A = Matrix<NX, NX, T>::zeros();
    Matrix<NX, NX, T> A_power = Matrix<NX, NX, T>::identity();
    for (size_t k = 0; k <= NX; ++k) {
        phi_A = phi_A + (cc[k].real() * A_power);
        A_power = A_power * A;
    }

    // K = e_Nᵀ · Co⁻¹ · φ(A). Solve Co·X = φ(A) instead of forming Co⁻¹;
    // a singular Co (uncontrollable) makes solve() fail.
    Matrix<1, NX, T> e_N{};
    e_N(0, NX - 1) = T{1};
    const auto CoInv_phi = mat::solve(Co, phi_A);
    if (!CoInv_phi) {
        return damp::nullopt;
    }
    return Matrix<1, NX, T>(e_N * (*CoInv_phi));
}

} // namespace design

} // namespace damp
