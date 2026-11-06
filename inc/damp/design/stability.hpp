// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file stability.hpp
 * @brief Stability analysis helpers
 *
 * @defgroup stability_analysis Stability Analysis
 * @brief Controllability, observability, and related structural stability checks
 *
 * For continuous systems: stable if all eigenvalues have Re(λ) < 0 (left half plane)
 * For discrete systems: stable if all eigenvalues have |λ| < 1 (inside unit circle)
 *
 * @note The eigenvalue-based routines (is_stable, stability_margin,
 *       closed_loop_poles) use mat::compute_eigenvalues, which resolves
 *       fully-complex eigenvalues for any N (closed form for N ≤ 2, Francis QR
 *       beyond).
 */
#include <cstddef>
#include <limits>

#include "damp/backend.hpp"
#include "damp/design/lyapunov.hpp"
#include "damp/math/complex.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/svd.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace stability {

// ============================================================================
// Structural Analysis (Controllability / Observability)
// ============================================================================

/**
 * @brief Compute the controllability matrix [B, AB, A²B, ..., A^(N-1)B]
 *
 * The system (A, B) is controllable iff controllability_matrix(A, B) has full row rank.
 *
 * @tparam NX Number of states
 * @tparam NU Number of inputs
 * @tparam T  Scalar type
 * @param A   State matrix
 * @param B   Input matrix
 * @return Matrix<NX, NX*NU, T> controllability matrix
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr Matrix<NX, NX * NU, T>
controllability_matrix(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B) noexcept {
    Matrix<NX, NX * NU, T> Co{};
    Matrix<NX, NU, T>      AB = B;
    for (size_t k = 0; k < NX; ++k) {
        for (size_t r = 0; r < NX; ++r) {
            for (size_t c = 0; c < NU; ++c) {
                Co(r, (k * NU) + c) = AB(r, c);
            }
        }
        if (k + 1 < NX) {
            AB = A * AB;
        }
    }
    return Co;
}

/**
 * @brief Compute the observability matrix [C; CA; CA²; ...; CA^(N-1)]
 *
 * The system (A, C) is observable iff observability_matrix(A, C) has full column rank.
 *
 * @tparam NX Number of states
 * @tparam NY Number of outputs
 * @tparam T  Scalar type
 * @param A   State matrix
 * @param C   Output matrix
 * @return Matrix<NX*NY, NX, T> observability matrix
 */
template<size_t NX, size_t NY, typename T = double>
[[nodiscard]] constexpr Matrix<NX * NY, NX, T>
observability_matrix(const Matrix<NX, NX, T>& A, const Matrix<NY, NX, T>& C) noexcept {
    Matrix<NX * NY, NX, T> Ob{};
    Matrix<NY, NX, T>      CA = C;
    for (size_t k = 0; k < NX; ++k) {
        for (size_t r = 0; r < NY; ++r) {
            for (size_t c = 0; c < NX; ++c) {
                Ob((k * NY) + r, c) = CA(r, c);
            }
        }
        if (k + 1 < NX) {
            CA = CA * A;
        }
    }
    return Ob;
}

namespace detail {

/**
 * @brief Gramian Lyapunov solve, refused unless A is Hurwitz (continuous) / Schur (discrete).
 * For an unstable A the Lyapunov equation can still have a unique solution, but it
 * is not the (divergent) Gramian integral/sum.
 */
template<size_t NX, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>>
gramian_lyapunov(const Matrix<NX, NX, T>& A, const Matrix<NX, NX, T>& Q, bool discrete) {
    const auto eig = mat::compute_eigenvalues(A);
    if (!eig.success) {
        return damp::nullopt;
    }
    for (size_t i = 0; i < NX; ++i) {
        const bool stable = discrete ? (damp::abs(eig.values[i]) < T{1}) : (eig.values[i].real() < T{0});
        if (!stable) {
            return damp::nullopt;
        }
    }
    return discrete ? dlyap(A, Q) : lyap(A, Q);
}

/** @brief Orthonormal basis of the controllable subspace of (A, B) and its dimension. */
template<size_t NX, typename T>
struct ControllableSubspace {
    Matrix<NX, NX, T> basis{};         ///< leading `dim` columns: an orthonormal basis of the subspace
    size_t            dim{0};          ///< its dimension (NX when (A, B) is controllable)
    bool              converged{true}; ///< every SVD converged
};

/**
 * @brief Controllable subspace by the block-Arnoldi staircase (the ctrbf construction).
 *
 * Stair 0 is range(B). Stair k adds the directions of A·(stair k−1) that are new: the
 * block is projected against the basis found so far (twice, so it is orthogonal to working
 * precision) and its left singular vectors above the tolerance are appended. A is scaled by
 * 1/‖A‖∞ and the basis columns have unit norm, so the test σ > @p tol is unit-free: a
 * direction counts when it is reached with a coupling of at least @p tol relative to the
 * plant's fastest rate; stair 0 is tested against tol·σₘₐₓ(B).
 *
 * The rank of the Krylov matrix [B, AB, …, Aⁿ⁻¹B] answers the same question in exact
 * arithmetic, but its columns align with the dominant mode and the matrix loses rank
 * numerically from n ≈ 10, or sooner when the modes spread over decades. Each stair here is
 * tested on its own, so the verdict does not degrade with n.
 *
 * @see Van Dooren, "The generalized eigenstructure problem in linear system theory,"
 *      IEEE TAC 26(1), 1981
 * @note Compare with MATLAB®'s ctrbf.
 */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr ControllableSubspace<NX, T>
controllable_subspace(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B, T tol) {
    constexpr size_t            K = (NX < NU) ? NX : NU;
    ControllableSubspace<NX, T> out{};
    const T                     a_norm = mat::infinity_norm(A);
    const Matrix<NX, NX, T>     As = (a_norm > T{0}) ? Matrix<NX, NX, T>(A * (T{1} / a_norm)) : A;
    Matrix<NX, NU, T>           W = B;
    for (size_t stair = 0; stair < NX && out.dim < NX; ++stair) {
        for (int pass = 0; pass < 2; ++pass) { // classical Gram–Schmidt, twice
            for (size_t j = 0; j < out.dim; ++j) {
                for (size_t c = 0; c < NU; ++c) {
                    T dot{0};
                    for (size_t i = 0; i < NX; ++i) {
                        dot += out.basis(i, j) * W(i, c);
                    }
                    for (size_t i = 0; i < NX; ++i) {
                        W(i, c) -= dot * out.basis(i, j);
                    }
                }
            }
        }
        const auto s = mat::svd(W);
        if (!s.success) {
            out.converged = false;
            return out;
        }
        const T floor = tol * ((stair == 0) ? s.singular_values[0] : T{1});
        size_t  r = 0;
        for (size_t k = 0; k < K; ++k) {
            if (s.singular_values[k] > floor) {
                ++r;
            }
        }
        if (r == 0) {
            break;
        }
        r = damp::min(r, NX - out.dim);
        for (size_t k = 0; k < r; ++k) {
            for (size_t i = 0; i < NX; ++i) {
                out.basis(i, out.dim + k) = s.singular_U(i, k);
            }
        }
        Matrix<NX, NU, T> next{};
        for (size_t k = 0; k < r; ++k) {
            for (size_t i = 0; i < NX; ++i) {
                T acc{0};
                for (size_t m = 0; m < NX; ++m) {
                    acc += As(i, m) * out.basis(m, out.dim + k);
                }
                next(i, k) = acc;
            }
        }
        W = next;
        out.dim += r;
    }
    return out;
}

/** @brief (A, B) in balanced coordinates: (D⁻¹ A D, D⁻¹ B), D = diag(2^e) from mat::detail::balance_exponents(A). */
template<size_t NX, size_t NU, typename T>
struct BalancedPair {
    damp::array<int, NX> e{}; ///< the scaling exponents
    Matrix<NX, NX, T>    A{}; ///< D⁻¹ A D
    Matrix<NX, NU, T>    B{}; ///< D⁻¹ B
};

/** @copybrief BalancedPair */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr BalancedPair<NX, NU, T> balanced_pair(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B) {
    BalancedPair<NX, NU, T> out{.e = mat::detail::balance_exponents(A)};
    out.A = mat::detail::balanced(A, out.e);
    for (size_t r = 0; r < NX; ++r) {
        for (size_t c = 0; c < NU; ++c) {
            out.B(r, c) = mat::detail::scale_by_pow2(B(r, c), -out.e[r]);
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
 * cannot hide the others. Then every row is brought to unit size, so a decoupled state's
 * input gain counts whatever the state's unit. The exception is a row whose off-diagonal and
 * input entries are all exactly zero (a mode with no input and no coupling): there λ − a_rr
 * cancels to λ's roundoff, and dividing by that would make it a unit entry and call the mode
 * controllable, so that row is divided by max(‖A‖_max, |λ|, 1) instead. Row scaling alone, on
 * an unbalanced pair, amplifies the eigenvalue's roundoff on the small rows by the scale
 * spread and calls an uncontrollable mode controllable. The tolerance is a relative √eps:
 * an uncontrollable eigenvalue is only known to O(√eps) when defective, so a tighter
 * tolerance would call it controllable.
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
    //! Each row to unit size. A row whose only entry is λ − a_rr (no coupling, no input) holds
    //! λ's roundoff (~ε) there, which must not become a unit entry: it is divided by the
    //! data's scale instead. a_scale ≥ 1, so no row is divided by 0.
    for (size_t r = 0; r < NX; ++r) {
        T data{0};
        for (size_t c = 0; c < NX + NU; ++c) {
            if (c != r) {
                data = damp::max(data, test_mat(r, c).abs());
            }
        }
        const T row = (data > T{0}) ? damp::max(data, test_mat(r, r).abs()) : damp::max(a_scale, lambda.abs());
        for (size_t c = 0; c < NX + NU; ++c) {
            test_mat(r, c) = test_mat(r, c) / Cplx(row, T{0});
        }
    }
    return mat::rank(test_mat, damp::sqrt(eps<T>())) < NX;
}


/**
 * @brief True when (A, B) has an uncontrollable eigenvalue at or beyond the stability boundary,
 *        or, with @p boundary_only, one on the boundary.
 *
 * The boundary bands are relative, so the verdict does not depend on the states' units or the
 * time unit: |λ| ≥ 1 − √ε (Discrete), Re λ ≥ −√ε·‖A‖ (continuous, A balanced); with
 * @p boundary_only also |λ| ≤ 1 + √ε, or Re λ ≤ √ε·‖A‖. An eigenvalue on the boundary (an
 * integrator, a random-walk disturbance) is computed only to roundoff on either side, and a
 * defective one to O(√ε) (a 2×2 Jordan block), so it is tested whichever side it lands on. A
 * row with no coupling and no input, eᵣᵀA = aᵣᵣ·eᵣᵀ and eᵣᵀB = 0, makes aᵣᵣ an uncontrollable
 * eigenvalue exactly; it is tested at aᵣᵣ, so a larger Jordan block, computed only to
 * O(ε^{1/k}), is still found there (a triple integrator whose Q does not see position). The
 * dual, (Aᵀ, Q), asks whether Q sees each such mode. An eigenvalue computation that does not
 * converge counts as such a mode.
 *
 * @see "Controllability and observability conditions of linear autonomous systems" (Hautus, 1969)
 */
template<bool Discrete, size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr bool uncontrollable_unstable_mode(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B, bool boundary_only) {
    const auto eigen = mat::compute_eigenvalues(A); // balanced internally
    if (!eigen.success) {
        return true;
    }
    const auto bal = balanced_pair(A, B);
    const T    slack = damp::sqrt(eps<T>());
    const T    band = Discrete ? slack : slack * bal.A.norm();
    //! The band is closed: a continuous A = 0 has band = 0 and its eigenvalue exactly on the axis.
    const auto skipped = [band, boundary_only](const damp::complex<T>& lambda) {
        const T dist = Discrete ? (lambda.abs() - T{1}) : lambda.real(); // signed distance past the boundary
        return dist < -band || (boundary_only && dist > band);
    };
    for (size_t r = 0; r < NX; ++r) {
        bool isolated = true;
        for (size_t c = 0; c < NX; ++c) {
            isolated = isolated && (c == r || A(r, c) == T{0});
        }
        for (size_t c = 0; c < NU; ++c) {
            isolated = isolated && (B(r, c) == T{0});
        }
        if (isolated && !skipped(damp::complex<T>(A(r, r), T{0}))) {
            return true; // aᵣᵣ is an uncontrollable eigenvalue exactly
        }
    }
    for (size_t i = 0; i < NX; ++i) {
        // The conjugate of a complex eigenvalue gives the same rank, so test one of each pair.
        const auto lambda = eigen.values[i];
        if (skipped(lambda) || lambda.imag() < T{0}) {
            continue;
        }
        if (uncontrollable_at(bal, lambda)) {
            return true;
        }
    }
    return false;
}

} // namespace detail

/**
 * @brief Continuous/discrete controllability Gramian W_c.
 *
 * Solves the Lyapunov equation whose solution measures how strongly each state
 * direction is excited by the input:
 * - continuous: @f$ A W_c + W_c A^\top + B B^\top = 0 @f$
 * - discrete:   @f$ A W_c A^\top - W_c + B B^\top = 0 @f$
 *
 * @f$ W_c \succ 0 @f$ iff @f$ (A,B) @f$ is controllable. The Gramian is finite only
 * for a stable @f$ A @f$ (Hurwitz / Schur); an unstable @f$ A @f$ returns
 * damp::nullopt even when the Lyapunov equation itself is solvable.
 *
 * @note Compare with MATLAB®'s @c gram(sys,'c').
 * @see observability_gramian(), damp::lyap(), damp::dlyap().
 *
 * @param A        State matrix (NX × NX).
 * @param B        Input matrix (NX × NU).
 * @param discrete false → continuous-time, true → discrete-time.
 * @return W_c (NX × NX), or damp::nullopt if A is not stable or the solve is singular.
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>> controllability_gramian(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    bool                     discrete = false
) {
    return detail::gramian_lyapunov(A, Matrix<NX, NX, T>{B * B.transpose()}, discrete);
}

/**
 * @brief Continuous/discrete observability Gramian W_o.
 *
 * Solves the dual Lyapunov equation measuring how strongly each state direction
 * shows up in the output:
 * - continuous: @f$ A^\top W_o + W_o A + C^\top C = 0 @f$
 * - discrete:   @f$ A^\top W_o A - W_o + C^\top C = 0 @f$
 *
 * @f$ W_o \succ 0 @f$ iff @f$ (A,C) @f$ is observable. The Gramian is finite only
 * for a stable @f$ A @f$; an unstable @f$ A @f$ returns damp::nullopt.
 *
 * @note Compare with MATLAB®'s @c gram(sys,'o').
 * @see controllability_gramian(), damp::lyap(), damp::dlyap().
 *
 * @param A        State matrix (NX × NX).
 * @param C        Output matrix (NY × NX).
 * @param discrete false → continuous-time, true → discrete-time.
 * @return W_o (NX × NX), or damp::nullopt if A is not stable or the solve is singular.
 */
template<size_t NX, size_t NY, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>> observability_gramian(
    const Matrix<NX, NX, T>& A,
    const Matrix<NY, NX, T>& C,
    bool                     discrete = false
) {
    return detail::gramian_lyapunov(A.transpose(), Matrix<NX, NX, T>{C.transpose() * C}, discrete);
}

/**
 * @brief Numerical rank of a matrix (forwards to mat::rank)
 *
 * @param M   Input matrix
 * @return size_t rank of the matrix
 */
template<size_t R, size_t C, typename T>
[[nodiscard]] constexpr size_t rank(const Matrix<R, C, T>& M) {
    return mat::rank(M);
}

/** @brief Numerical rank with an explicit relative tolerance: singular values ≤ tol·σ_max count as zero. */
template<size_t R, size_t C, typename T>
[[nodiscard]] constexpr size_t rank(const Matrix<R, C, T>& M, T tol) {
    return mat::rank(M, tol);
}

/**
 * @brief Check if a system is controllable
 *
 * The controllable subspace (detail::controllable_subspace, the staircase construction)
 * has dimension NX. The tolerance NX·64·eps is relative to the plant's fastest rate and
 * to B, so the answer is invariant to the signal and time units.
 *
 * @param A State matrix
 * @param B Input matrix
 * @return true if (A, B) is controllable
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr bool is_controllable(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B) {
    return detail::controllable_subspace(A, B, static_cast<T>(NX) * T{64} * eps<T>()).dim == NX;
}

/**
 * @brief Check if a system is observable
 *
 * Dual of is_controllable(): the controllable subspace of (Aᵀ, Cᵀ) has dimension NX.
 *
 * @param A State matrix
 * @param C Output matrix
 * @return true if (A, C) is observable
 */
template<size_t NX, size_t NY, typename T = double>
[[nodiscard]] constexpr bool is_observable(const Matrix<NX, NX, T>& A, const Matrix<NY, NX, T>& C) {
    return detail::controllable_subspace(Matrix<NX, NX, T>(A.transpose()), Matrix<NX, NY, T>(C.transpose()), static_cast<T>(NX) * T{64} * eps<T>()).dim == NX;
}

/**
 * @brief Check if (A, B) of @p sys is a stabilizable pair, in the domain of @p sys
 *
 * (A, B) is stabilizable iff every uncontrollable eigenvalue of A is stable: |λ| < 1
 * (discrete, sys.Ts > 0) or Re λ < 0 (continuous, sys.Ts = 0). An eigenvalue λ is
 * uncontrollable when rank([λI − A, B]) < n (the PBH test).
 *
 * The pair is balanced first (mat::detail::balance_exponents), so the verdict does not
 * depend on the units of the states. An eigenvalue within √ε of the boundary (relative to
 * ‖A‖ when continuous) counts as unstable: one on it (an integrator, a random-walk
 * disturbance) is computed to roundoff on either side, and the Riccati certificates use the
 * same band.
 *
 * @see "Optimal Control" (Anderson & Moore, 1990), §2.4
 *
 * @param sys State-space model; only A, B and Ts are used
 * @return true if (A, B) is stabilizable
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr bool is_stabilizable(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    return sys.Ts > T{0} ? !detail::uncontrollable_unstable_mode<true>(sys.A, sys.B, false)
                         : !detail::uncontrollable_unstable_mode<false>(sys.A, sys.B, false);
}

/**
 * @brief Stability of a state-space model, in the domain of @p sys
 *
 * Continuous (sys.Ts = 0): every eigenvalue of A has Re(λ) < 0. Discrete (sys.Ts > 0):
 * every eigenvalue has |λ| < 1. An eigenvalue computation that does not converge
 * reports unstable.
 *
 * @param sys State-space model; only A and Ts are used
 * @return true if every pole is inside the domain's stability region
 *
 * @note Compare with MATLAB®'s isstable(sys).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr bool is_stable(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    const auto eigen = mat::compute_eigenvalues(sys.A);
    if (!eigen.success) {
        return false;
    }
    for (size_t i = 0; i < NX; ++i) {
        if (sys.Ts > T{0} ? !(damp::abs(eigen.values[i]) < T{1}) : !(eigen.values[i].real() < T{0})) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Distance from the poles of @p sys to its domain's stability boundary
 *
 * Continuous (sys.Ts = 0): −max Re(λ), the distance to the imaginary axis. Discrete
 * (sys.Ts > 0): 1 − max |λ|, the distance to the unit circle. Positive is stable; a
 * non-converged eigenvalue computation returns −1. A system with no states (a static gain)
 * returns 1 when discrete and the largest finite T when continuous.
 *
 * @param sys State-space model; only A and Ts are used
 * @return Stability margin (positive = stable, negative = unstable)
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr T stability_margin(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    const auto eigen = mat::compute_eigenvalues(sys.A);
    if (!eigen.success) {
        return T{-1}; // negative margin signals "not provably stable"
    }
    T worst = sys.Ts > T{0} ? T{0} : -std::numeric_limits<T>::max(); // finite sentinel: NX = 0 has no eigenvalue to start from
    for (size_t i = 0; i < NX; ++i) {
        const T v = sys.Ts > T{0} ? eigen.values[i].abs() : eigen.values[i].real();
        if (v > worst) {
            worst = v;
        }
    }
    return sys.Ts > T{0} ? (T{1} - worst) : -worst;
}

/**
 * @brief Compute closed-loop poles (eigenvalues) with state feedback
 *
 * Computes the eigenvalues of the closed-loop state matrix A_cl = A - B*K.
 * These poles determine the closed-loop system dynamics.
 *
 * If the eigen solver does not converge, returns a zero vector rather than
 * uninitialized/garbage pole locations. Callers that set design @c success from
 * the Riccati/gain path alone (e.g. LQR) keep @c success true when K is valid;
 * a zeroed pole vector signals that pole reporting failed independently of K.
 *
 * @tparam NX  Number of states
 * @tparam NU  Number of inputs
 * @tparam T   Numeric type (default: double)
 * @param A    State matrix
 * @param B    Control input matrix
 * @param K    State feedback gain matrix (u = -K*x)
 *
 * @return Vector of closed-loop pole locations (eigenvalues as complex numbers);
 *         all zeros if eigenvalue computation does not converge
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr damp::array<damp::complex<T>, NX> closed_loop_poles(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NU, NX, T>& K
) {
    Matrix<NX, NX, T>                 A_cl = A - B * K;
    auto                              eigen = mat::compute_eigenvalues(A_cl);
    damp::array<damp::complex<T>, NX> poles{};
    if (eigen.success) {
        for (size_t i = 0; i < NX; ++i) {
            poles[i] = eigen.values[i];
        }
    }
    return poles;
}

} // namespace stability
} // namespace damp
