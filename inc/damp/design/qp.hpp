// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file qp.hpp
 * @brief Dense convex quadratic programming: Goldfarb–Idnani dual active-set
 *        solver over the library's fixed-size matrices.
 *
 * Solves the strictly convex QP
 * @f[
 *   \min_x \tfrac{1}{2} x^\top H x + f^\top x
 *   \quad \text{s.t.} \quad A x \le b
 * @f]
 * with @f$ H \succ 0 @f$, allocation-free and `constexpr`-capable. This is the
 * shared small dense QP used by constrained MPC (roadmap #14) and the
 * constraint-enforcement safety filters (roadmap #34).
 *
 * Rows of @f$ A x \le b @f$ whose bound is the unbounded sentinel
 * (@f$ b_i \ge \texttt{max}/4 @f$, see unbounded_bound) are ignored, so a
 * fixed constraint block can carry optional limits without rebuilding the
 * matrices (the library builds with -ffinite-math-only; sentinels are
 * `numeric_limits` extremes, never infinity).
 *
 * @note Compare with MATLAB®'s quadprog(H, f, A, b).
 * @see Goldfarb & Idnani, "A numerically stable dual method for solving
 *      strictly convex quadratic programs," Mathematical Programming 27, 1983,
 *      https://doi.org/10.1007/BF02591962
 * @see Nocedal & Wright, "Numerical Optimization," 2nd ed., 2006, §16.5
 */

#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/solve.hpp"

namespace damp {

namespace design {

/**
 * @brief Termination status of a QP solve
 */
enum class QPStatus : uint8_t {
    Success,             ///< KKT point found: primal feasible, dual feasible, complementary
    MaxIterations,       ///< Iteration budget exhausted — x is NOT guaranteed primal feasible (dual method)
    Infeasible,          ///< Constraints admit no solution
    Degenerate,          ///< Active-set normals became numerically dependent
    NotPositiveDefinite, ///< H failed its Cholesky factorization
    Unbounded,           ///< (ADMM) The objective decreases without bound on the feasible set: H d = 0, A d ≤ 0, fᵀd < 0
    NonFinite,           ///< The data or an iterate holds a NaN or an infinity: no solution is reported
};

/**
 * @brief Sentinel bound treated as "no constraint" on that row
 *
 * Any row with @f$ b_i \ge @f$ unbounded_bound()/4 is skipped entirely, so
 * the arithmetic @f$ n_i^\top x - b_i @f$ is never evaluated on sentinel rows
 * (finite-math safe).
 */
template<typename T>
[[nodiscard]] constexpr T unbounded_bound() {
    return std::numeric_limits<T>::max();
}

/**
 * @struct QPResult
 * @brief Result of a dense QP solve
 *
 * Contains the minimizer x, the Lagrange multipliers of the inequality
 * constraints, the objective value, and the iteration count actually used.
 * Use .as<float>() to convert for embedded deployment.
 *
 * @note Compare with MATLAB®'s [x, fval, exitflag, ~, lambda] = quadprog(...) outputs.
 */
template<size_t NV, size_t NI, typename T = double>
struct QPResult {
    ColVec<NV, T> x{};                             ///< Minimizer
    ColVec<NI, T> lambda{};                        ///< Multipliers: λᵢ ≥ 0, λᵢ·(Ax − b)ᵢ = 0
    T             objective{};                     ///< ½xᵀHx + fᵀx at x
    size_t        iterations{0};                   ///< Active-set changes performed
    QPStatus      status{QPStatus::MaxIterations}; ///< Termination status
    bool          success{false};                  ///< true iff status == Success

    template<typename U>
    [[nodiscard]] constexpr QPResult<NV, NI, U> as() const {
        return QPResult<NV, NI, U>{
            x.template as<U>(),
            lambda.template as<U>(),
            static_cast<U>(objective),
            iterations,
            status,
            success,
        };
    }
};

namespace detail {

/** @brief A ColVec with every entry set to value (bound-vector construction helper) */
template<size_t N, typename T>
[[nodiscard]] constexpr ColVec<N, T> filled_vector(T value) {
    ColVec<N, T> out{};
    for (size_t i = 0; i < N; ++i) {
        out(i) = value;
    }
    return out;
}

/**
 * @brief Bound-vector scalar conversion that maps unbounded sentinels
 *
 * A plain static_cast of numeric_limits<double>::max() to float overflows
 * (undefined under -ffinite-math-only), so sentinel entries are re-issued in
 * the target type instead of cast.
 */
template<typename U, size_t N, typename T>
[[nodiscard]] constexpr ColVec<N, U> convert_bounds(const ColVec<N, T>& v) {
    constexpr T  hi = std::numeric_limits<T>::max() / T{4};
    ColVec<N, U> out{};
    for (size_t i = 0; i < N; ++i) {
        if (v(i) >= hi) {
            out(i) = std::numeric_limits<U>::max();
        } else if (v(i) <= -hi) {
            out(i) = -std::numeric_limits<U>::max();
        } else {
            out(i) = static_cast<U>(v(i));
        }
    }
    return out;
}

/**
 * @brief constexpr floor for branch-and-bound bounds. A value at or past 1/ε in magnitude is
 *        already an integer and is returned as is, as is a NaN: the long long cast never sees them.
 */
template<typename T>
[[nodiscard]] constexpr T floor_scalar(T v) {
    if (!(damp::abs(v) < T{1} / eps<T>())) {
        return v;
    }
    T fl = static_cast<T>(static_cast<long long>(v));
    if (fl > v) {
        fl -= T{1};
    }
    return fl;
}

/** @brief constexpr round-to-nearest for integrality checks */
template<typename T>
[[nodiscard]] constexpr T round_scalar(T v) {
    return floor_scalar(v + static_cast<T>(0.5));
}

/**
 * @brief Relative pivot floor for the active-set Gram factor: a pivot below this
 * fraction of its diagonal entry means the new normal is (numerically) in the
 * span of the active ones. Relative, so it is invariant to scaling H, f, A, b.
 */
template<size_t KMAX, typename T>
[[nodiscard]] constexpr T gram_pivot_tol() {
    return static_cast<T>(KMAX + 1) * T{64} * eps<T>();
}

/**
 * @brief Cholesky-factor the leading k×k block of M into L
 *
 * Returns false if a pivot falls below gram_pivot_tol() of its diagonal entry
 * (numerically dependent normals in the active-set method).
 */
template<size_t KMAX, typename T>
constexpr bool cholesky_leading(const Matrix<KMAX, KMAX, T>& M, Matrix<KMAX, KMAX, T>& L, size_t k) {
    constexpr T tol = gram_pivot_tol<KMAX, T>();

    k = (k < KMAX) ? k : KMAX; // provable index bound for the optimizer
    for (size_t i = 0; i < k; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            T sum = M(i, j);
            for (size_t c = 0; c < j; ++c) {
                sum -= L(i, c) * L(j, c);
            }
            if (i == j) {
                if (!(sum > tol * M(i, i))) {
                    return false;
                }
                L(i, i) = damp::sqrt(sum);
            } else {
                L(i, j) = sum / L(j, j);
            }
        }
    }
    return true;
}

/** @brief Solve L Lᵀ r = rhs in place on the leading k×k block of a cholesky_leading factor */
template<size_t KMAX, typename T>
constexpr void factored_solve_leading(const Matrix<KMAX, KMAX, T>& L, ColVec<KMAX, T>& rhs, size_t k) {
    k = (k < KMAX) ? k : KMAX; // provable index bound for the optimizer
    for (size_t i = 0; i < k; ++i) {
        T sum = rhs(i);
        for (size_t j = 0; j < i; ++j) {
            sum -= L(i, j) * rhs(j);
        }
        rhs(i) = sum / L(i, i);
    }
    for (size_t ii = k; ii > 0; --ii) {
        const size_t i = ii - 1;
        T            sum = rhs(i);
        for (size_t j = i + 1; j < k; ++j) {
            sum -= L(j, i) * rhs(j);
        }
        rhs(i) = sum / L(i, i);
    }
}

/**
 * @brief Grow a Cholesky factor by one row: O(k²) instead of an O(k³) refactor
 *
 * L holds the factor of the leading k×k block; col carries the new matrix
 * column (off-diagonal entries in col(0..k−1), diagonal in col(k)). Returns
 * false if the new pivot falls below gram_pivot_tol() of col(k) (dependent normal).
 */
template<size_t KMAX, typename T>
constexpr bool cholesky_append(Matrix<KMAX, KMAX, T>& L, const ColVec<KMAX, T>& col, size_t k) {
    constexpr T tol = gram_pivot_tol<KMAX, T>();

    if (k >= KMAX) {
        return false; // row k would fall outside the fixed-size factor
    }
    // Forward-substitute the new off-diagonal row: L(0:k,0:k) l = col(0:k).
    T diag_sq = col(k);
    for (size_t j = 0; j < k; ++j) {
        T sum = col(j);
        for (size_t c = 0; c < j; ++c) {
            sum -= L(j, c) * L(k, c);
        }
        L(k, j) = sum / L(j, j);
        diag_sq -= L(k, j) * L(k, j);
    }
    if (!(diag_sq > tol * col(k))) {
        return false;
    }
    L(k, k) = damp::sqrt(diag_sq);
    return true;
}

/**
 * @brief Goldfarb–Idnani dual active-set core with a pre-factored Hessian
 *
 * Shared by solve_qp() (cold: factors H per call, empty hint) and
 * WarmStartActiveSetSolver (caches the factor and seeds the hint from the
 * previous solve). The Gram matrix M = NᵀH⁻¹N of the active normals and its
 * Cholesky factor are maintained incrementally: an O(k²) append per activation
 * (the new column falls out of quantities the step computation already needs),
 * with a full O(k³) refactor only on the rare constraint drop.
 *
 * The hint restricts each constraint pick to the most-violated *hinted* row
 * while any hinted row is violated, falling back to the full most-violated
 * scan — with a good hint (the previous solve's active set) the pick costs
 * O(hint·NV) instead of O(NI·NV). Any violated constraint is a valid pick, so
 * a stale or wrong hint affects speed, never correctness.
 */
template<size_t NV, size_t NI, typename T>
[[nodiscard]] constexpr QPResult<NV, NI, T> gi_active_set(
    const Matrix<NV, NV, T>&                        H,
    const Matrix<NV, NV, T>&                        L,
    const ColVec<NV, T>&                            f,
    const Matrix<NI, NV, T>&                        A,
    const ColVec<NI, T>&                            b,
    size_t                                          max_iterations,
    const damp::array<size_t, (NV < NI) ? NV : NI>& hint,
    size_t                                          hint_count
) {
    constexpr size_t KMAX = (NV < NI) ? NV : NI;
    // All tests below are relative to the magnitudes that form the tested
    // quantity, so scaling H, f, or a constraint row (aᵢ, bᵢ) never changes the
    // active-set path.
    constexpr T pivot_tol = gram_pivot_tol<KMAX, T>();
    constexpr T slack_tol = static_cast<T>(NV + 1) * T{64} * eps<T>();
    constexpr T big = unbounded_bound<T>();
    constexpr T skip_threshold = big / T{4};

    QPResult<NV, NI, T> result{};

    // x = −H⁻¹ f via the supplied factor
    const auto h_solve = [&L](const ColVec<NV, T>& v) {
        return mat::backward_substitute_transpose(L, mat::forward_substitute(L, v));
    };
    ColVec<NV, T> x = -h_solve(f);

    // Internal ≥ convention: nᵢ = −A(i,:)ᵀ, βᵢ = −bᵢ, slack sᵢ = nᵢᵀx − βᵢ ≥ 0.
    const auto normal = [&A](size_t i) {
        ColVec<NV, T> n{};
        for (size_t c = 0; c < NV; ++c) {
            n(c) = -A(i, c);
        }
        return n;
    };
    const auto slack = [&](size_t i) {
        T s = b(i);
        for (size_t c = 0; c < NV; ++c) {
            s -= A(i, c) * x(c);
        }
        return s; // = nᵢᵀx − βᵢ
    };
    bool       non_finite = false; // a NaN or infinite slack: a non-finite bᵢ, row of A, or x
    const auto violated = [&](size_t i, T& s_out) {
        if (b(i) >= skip_threshold) {
            return false;
        }
        // Roundoff scale of sᵢ = bᵢ − aᵢᵀx: slack_tol·(|bᵢ| + Σ|A(i,c) x(c)|) for the row's own
        // terms, plus (NV+1)·ε·‖aᵢ‖₁‖x‖∞ for the error the whole solve leaves in x, so a state
        // pinned at 0 by an equal pair of limits does not read its paired row as violated. The
        // second term stays ε-sized: slack_tol·‖x‖∞ accepted real violations in float on a small
        // variable beside a large one.
        T mag = damp::abs(b(i));
        T a_norm{0};
        T x_max{0};
        for (size_t c = 0; c < NV; ++c) {
            mag += damp::abs(A(i, c) * x(c));
            a_norm += damp::abs(A(i, c));
            x_max = damp::max(x_max, damp::abs(x(c)));
        }
        const T tol = (slack_tol * mag) + (static_cast<T>(NV + 1) * eps<T>() * a_norm * x_max);
        s_out = slack(i);
        if (!damp::isfinite(s_out) || !damp::isfinite(tol)) {
            non_finite = true; // a NaN compares false: the row would read as satisfied
            return false;
        }
        return s_out < -tol;
    };

    // Active set: indices, multipliers, cached wⱼ = H⁻¹nⱼ, and the Gram matrix
    // M = NᵀH⁻¹N with its Cholesky factor LM (leading k×k blocks valid).
    damp::array<size_t, KMAX>        active{};
    damp::array<T, KMAX>             lam{};
    damp::array<ColVec<NV, T>, KMAX> w{};
    Matrix<KMAX, KMAX, T>            M{};
    Matrix<KMAX, KMAX, T>            LM{};
    damp::array<bool, NI>            in_active{};
    size_t                           k = 0;

    const auto finish = [&](QPStatus status) {
        for (size_t a = 0; a < k; ++a) {
            result.lambda(active[a]) = lam[a];
        }
        result.x = x;
        T obj = T{0};
        for (size_t r = 0; r < NV; ++r) {
            T Hx = T{0};
            for (size_t c = 0; c < NV; ++c) {
                Hx += H(r, c) * x(c);
            }
            obj += x(r) * (static_cast<T>(0.5) * Hx + f(r));
        }
        result.objective = obj;
        result.status = status;
        result.success = (status == QPStatus::Success);
        return result;
    };

    for (size_t v = 0; v < NV; ++v) {
        if (!damp::isfinite(x(v))) {
            return finish(QPStatus::NonFinite); // a non-finite H or f
        }
    }

    while (true) {
        // Pick a violated constraint: most-violated among the hinted rows first
        // (an O(hint·NV) scan — at steady state the hint is the whole active
        // set), falling back to the full O(NI·NV) most-violated scan.
        size_t p = NI;
        T      s_i = T{0};
        T      worst = T{0};
        for (size_t h = 0; h < hint_count; ++h) {
            const size_t i = hint[h];
            if (i >= NI || in_active[i]) {
                continue;
            }
            if (violated(i, s_i) && s_i < worst) {
                worst = s_i;
                p = i;
            }
        }
        if (p == NI) {
            for (size_t i = 0; i < NI; ++i) {
                if (in_active[i]) {
                    continue;
                }
                if (violated(i, s_i) && s_i < worst) {
                    worst = s_i;
                    p = i;
                }
            }
        }
        if (non_finite) {
            return finish(QPStatus::NonFinite);
        }
        if (p == NI) {
            return finish(QPStatus::Success);
        }

        const ColVec<NV, T> n_p = normal(p);
        const ColVec<NV, T> u_p = h_solve(n_p); // H⁻¹nₚ — fixed while p is the candidate
        T                   npup = T{0};
        for (size_t v = 0; v < NV; ++v) {
            npup += n_p(v) * u_p(v);
        }
        // nu(a) = n_aᵀuₚ — doubles as the new Gram column when p is activated.
        ColVec<KMAX, T> nu{};
        for (size_t a = 0; a < k; ++a) {
            T dot = T{0};
            for (size_t v = 0; v < NV; ++v) {
                dot += A(active[a], v) * u_p(v); // n_a = −A(a,:) …
            }
            nu(a) = -dot; // … so negate once
        }
        T lam_p = T{0};

        // Inner loop: step toward constraint p, dropping blockers as needed.
        while (true) {
            // r = M⁻¹ nu via the maintained factor; d = uₚ − Σ r_a w_a.
            ColVec<KMAX, T> r = nu;
            detail::factored_solve_leading(LM, r, k);

            ColVec<NV, T> d = u_p;
            for (size_t a = 0; a < k; ++a) {
                for (size_t v = 0; v < NV; ++v) {
                    d(v) -= r(a) * w[a](v);
                }
            }

            T kappa = T{0};
            for (size_t v = 0; v < NV; ++v) {
                kappa += n_p(v) * d(v);
            }

            // κ = nₚᵀ(I − projection)H⁻¹nₚ is the Schur complement of the new
            // Gram column: ≈ 0 relative to nₚᵀH⁻¹nₚ means nₚ is dependent on the
            // active normals, so only a dual (drop) step is possible. With KMAX
            // normals active they span ℝ^NV and κ is zero up to roundoff.
            const T    s_p = slack(p);
            const bool dependent = !(k < KMAX && kappa > T{2} * pivot_tol * npup);
            const T    t2 = dependent ? big : (-s_p / kappa);

            T      t1 = big;
            size_t drop = KMAX;
            for (size_t a = 0; a < k; ++a) {
                // r_a carries units ‖nₚ‖/‖n_a‖ (H⁻¹ metric); compare in those units.
                const T r_floor = pivot_tol * damp::sqrt(npup / M(a, a));
                if (r(a) > r_floor && lam[a] / r(a) < t1) {
                    t1 = lam[a] / r(a);
                    drop = a;
                }
            }

            if (t1 >= big && t2 >= big) {
                // Infeasible needs a Farkas certificate: y = (−r, 1) ≥ 0 on (active, p) with
                // Σ yᵢaᵢ ≈ 0 and Σ yᵢbᵢ < 0. Without one, nₚ is only nearly dependent in the
                // H⁻¹ metric, where κ fell to roundoff (an ill-conditioned H in float), and the
                // active set cannot grow: Degenerate.
                T res = T{0};
                T res_mag = T{0};
                for (size_t v = 0; v < NV; ++v) {
                    T rv = A(p, v);
                    T mv = damp::abs(A(p, v));
                    for (size_t a = 0; a < k; ++a) {
                        rv -= r(a) * A(active[a], v);
                        mv += damp::abs(r(a) * A(active[a], v));
                    }
                    res = damp::max(res, damp::abs(rv));
                    res_mag = damp::max(res_mag, mv);
                }
                T gap = b(p);
                T gap_mag = damp::abs(b(p));
                for (size_t a = 0; a < k; ++a) {
                    gap -= r(a) * b(active[a]);
                    gap_mag += damp::abs(r(a) * b(active[a]));
                }
                const bool certified = (res <= damp::sqrt(pivot_tol) * res_mag) && (gap < -(slack_tol * gap_mag));
                return finish(certified ? QPStatus::Infeasible : QPStatus::Degenerate);
            }

            // A dependent nₚ has primal direction d = 0, so the step is in the dual alone
            // (Goldfarb & Idnani's partial step with z = 0): the computed d is roundoff, and
            // the step length, a multiplier as large as the cost's weights, would carry it
            // into x.
            const T t = (t1 < t2) ? t1 : t2;
            if (!dependent) {
                for (size_t v = 0; v < NV; ++v) {
                    x(v) += t * d(v);
                }
            }
            for (size_t a = 0; a < k; ++a) {
                lam[a] -= t * r(a);
            }
            lam_p += t;

            ++result.iterations;
            if (result.iterations >= max_iterations) {
                if (k < KMAX && t2 <= t1) {
                    active[k] = p; // multiplier bookkeeping only — returning immediately
                    lam[k] = lam_p;
                    w[k] = u_p;
                    ++k;
                }
                return finish(QPStatus::MaxIterations);
            }

            if (t2 <= t1) {
                // Full step: p joins the active set. The Gram column is exactly
                // (nu, npup), so the factor grows in O(k²).
                if (k >= KMAX) {
                    return finish(QPStatus::Degenerate); // unreachable (t2 is finite only with k < KMAX): the index bound for the optimizer
                }
                for (size_t a = 0; a < k; ++a) {
                    M(a, k) = nu(a);
                    M(k, a) = nu(a);
                }
                M(k, k) = npup;
                ColVec<KMAX, T> col = nu;
                col(k) = npup;
                if (!detail::cholesky_append(LM, col, k)) {
                    return finish(QPStatus::Degenerate);
                }
                active[k] = p;
                lam[k] = lam_p;
                w[k] = u_p;
                in_active[p] = true;
                ++k;
                break;
            }

            // Partial (dual) step: drop the blocking constraint (swap-remove +
            // symmetric permutation of M), refactor — drops are rare, adds are
            // the O(k²) fast path.
            in_active[active[drop]] = false;
            const size_t last = k - 1;
            if (drop != last) {
                active[drop] = active[last];
                lam[drop] = lam[last];
                w[drop] = w[last];
                nu(drop) = nu(last);
                for (size_t j = 0; j < k; ++j) {
                    damp::swap(M(drop, j), M(last, j));
                }
                for (size_t j = 0; j < k; ++j) {
                    damp::swap(M(j, drop), M(j, last));
                }
            }
            k = last;
            if (!detail::cholesky_leading(M, LM, k)) {
                return finish(QPStatus::Degenerate);
            }
        }
    }
}

} // namespace detail

/**
 * @brief Solve a strictly convex inequality-constrained QP (dual active-set)
 *
 * Minimizes
 * @f[
 *   \tfrac{1}{2} x^\top H x + f^\top x \quad \text{s.t.} \quad A x \le b
 * @f]
 * by the Goldfarb–Idnani dual active-set method: start at the unconstrained
 * minimum @f$ x = -H^{-1} f @f$ and add violated constraints one at a time,
 * taking dual steps (dropping blocking constraints) as needed. Iterates are
 * dual feasible throughout; primal feasibility holds only at Success — a
 * MaxIterations result may violate constraints and callers must clamp or
 * reject it.
 *
 * Each active-set change costs one @f$ O(N_V^2) @f$ solve against the cached
 * Cholesky factor of H plus one @f$ O(k^3) @f$ dense solve in the number of
 * active constraints k.
 *
 * @note Compare with MATLAB®'s quadprog(H, f, A, b) ('active-set' algorithm).
 * @see QPResult, unbounded_bound()
 * @see Goldfarb & Idnani, Mathematical Programming 27, 1983,
 *      https://doi.org/10.1007/BF02591962
 *
 * @code
 * // MATLAB® quadprog documentation example: minimum at x = [2/3, 4/3]
 * constexpr Matrix<2, 2> H{{1.0, -1.0}, {-1.0, 2.0}};
 * constexpr ColVec<2>    f{-2.0, -6.0};
 * constexpr Matrix<3, 2> A{{1.0, 1.0}, {-1.0, 2.0}, {2.0, 1.0}};
 * constexpr ColVec<3>    b{2.0, 2.0, 3.0};
 * constexpr auto res = design::solve_qp(H, f, A, b);
 * static_assert(res.success);
 * @endcode
 *
 * @param H               Hessian (NV × NV, symmetric positive definite)
 * @param f               Linear cost term (NV)
 * @param A               Inequality constraint matrix (NI × NV), rows aᵢᵀx ≤ bᵢ
 * @param b               Inequality bounds (NI); entries ≥ unbounded_bound()/4 disable the row
 * @param max_iterations  Active-set change budget (default 10·(NV + NI))
 * @return QPResult with minimizer, multipliers, objective, and status
 */
template<size_t NV, size_t NI, typename T = double>
    requires std::is_floating_point_v<T>
[[nodiscard]] constexpr QPResult<NV, NI, T> solve_qp(
    const Matrix<NV, NV, T>& H,
    const ColVec<NV, T>&     f,
    const Matrix<NI, NV, T>& A,
    const ColVec<NI, T>&     b,
    size_t                   max_iterations = 10 * (NV + NI)
) {
    QPResult<NV, NI, T> result{};

    const auto L_opt = mat::cholesky(H);
    if (!L_opt) {
        result.status = QPStatus::NotPositiveDefinite;
        return result;
    }

    constexpr damp::array<size_t, (NV < NI) ? NV : NI> no_hint{};
    return detail::gi_active_set(H, L_opt.value(), f, A, b, max_iterations, no_hint, 0);
}

/**
 * @brief Solve a strictly convex QP with equality and inequality constraints
 *
 * Minimizes
 * @f[
 *   \tfrac{1}{2} x^\top H x + f^\top x \quad \text{s.t.} \quad
 *   A x \le b,\;\; A_{eq} x = b_{eq}
 * @f]
 * by encoding each equality as the paired inequalities
 * @f$ a^\top x \le b @f$ and @f$ -a^\top x \le -b @f$ and running the same
 * Goldfarb–Idnani solver. The twins are mutually exclusive in the active set
 * (violating one means satisfying the other), so the pairing is safe; it costs
 * two constraint rows per equality and lets the equality's multiplier change
 * sign by swapping which twin is active.
 * // ponytail: paired-row encoding, no native equality phase; a dedicated
 * // never-drop equality pass would save a few iterations if it ever matters.
 *
 * The multiplier of equality e is lambda(NI+e) − lambda(NI+NE+e) in the
 * returned result (upper twin minus lower twin). Inconsistent equalities are
 * reported as Infeasible.
 *
 * @note Compare with MATLAB®'s quadprog(H, f, A, b, Aeq, beq).
 * @see solve_qp
 *
 * @param H               Hessian (NV × NV, symmetric positive definite)
 * @param f               Linear cost term (NV)
 * @param A               Inequality constraint matrix (NI × NV), rows aᵢᵀx ≤ bᵢ
 * @param b               Inequality bounds; sentinel rows (≥ max/4) are skipped
 * @param Aeq             Equality constraint matrix (NE × NV)
 * @param beq             Equality values (NE)
 * @param max_iterations  Active-set change budget
 * @return QPResult over NI + 2·NE rows (equalities appended as twin pairs)
 */
template<size_t NV, size_t NI, size_t NE, typename T = double>
    requires std::is_floating_point_v<T>
[[nodiscard]] constexpr QPResult<NV, NI + (2 * NE), T> solve_qp(
    const Matrix<NV, NV, T>& H,
    const ColVec<NV, T>&     f,
    const Matrix<NI, NV, T>& A,
    const ColVec<NI, T>&     b,
    const Matrix<NE, NV, T>& Aeq,
    const ColVec<NE, T>&     beq,
    size_t                   max_iterations = 10 * (NV + NI + (2 * NE))
) {
    Matrix<NI + (2 * NE), NV, T> A_ext{};
    ColVec<NI + (2 * NE), T>     b_ext{};
    for (size_t i = 0; i < NI; ++i) {
        for (size_t v = 0; v < NV; ++v) {
            A_ext(i, v) = A(i, v);
        }
        b_ext(i) = b(i);
    }
    for (size_t e = 0; e < NE; ++e) {
        for (size_t v = 0; v < NV; ++v) {
            A_ext(NI + e, v) = Aeq(e, v);
            A_ext(NI + NE + e, v) = -Aeq(e, v);
        }
        b_ext(NI + e) = beq(e);
        b_ext(NI + NE + e) = -beq(e);
    }
    return solve_qp(H, f, A_ext, b_ext, max_iterations);
}

/**
 * @brief Default QP solver policy: the Goldfarb–Idnani active-set solve_qp()
 *
 * damp::MPC is templated on a solver policy with this call shape so alternative
 * solvers (ADMM, interior-point, MIQP) can be swapped in without touching the
 * controller. Policies are objects, not free functions — a solver may carry
 * state (warm-start caches, factorizations) across calls. This default is
 * stateless.
 *
 * @see solve_qp, damp::MPC
 */
struct ActiveSetSolver {
    template<size_t NV, size_t NI, typename T>
    [[nodiscard]] constexpr QPResult<NV, NI, T> operator()(
        const Matrix<NV, NV, T>& H,
        const ColVec<NV, T>&     f,
        const Matrix<NI, NV, T>& A,
        const ColVec<NI, T>&     b,
        size_t                   max_iterations
    ) const {
        return solve_qp(H, f, A, b, max_iterations);
    }
};

/**
 * @brief Warm-started active-set solver policy — the damp::MPC default
 *
 * Same Goldfarb–Idnani algorithm as solve_qp(), with the two cross-call
 * savings that matter on a small microcontroller:
 *
 * - the Cholesky factor of H is computed once and cached, removing the
 *   O(NV³) factorization from every tick, and
 * - the active set of the previous solve seeds the next one, so each
 *   constraint pick costs O(NV) against the hint instead of an O(NI·NV)
 *   most-violated scan, and add/drop churn from picking the wrong
 *   constraints first is avoided.
 *
 * A hint is only ever a hint — a stale or wrong one costs iterations, never
 * correctness. The cache assumes H and A are constant across calls (true for
 * damp::MPC, whose QP matrices are fixed at synthesis); call reset() if they
 * change. RAM cost over the stateless policy: one NV×NV factor plus the hint
 * indices.
 *
 * @see solve_qp, ActiveSetSolver, damp::MPC
 */
template<size_t NV, size_t NI, typename T = double>
    requires std::is_floating_point_v<T>
class WarmStartActiveSetSolver {
public:
    static constexpr size_t KMAX = (NV < NI) ? NV : NI;

    constexpr WarmStartActiveSetSolver() = default;

    [[nodiscard]] constexpr QPResult<NV, NI, T> operator()(
        const Matrix<NV, NV, T>& H,
        const ColVec<NV, T>&     f,
        const Matrix<NI, NV, T>& A,
        const ColVec<NI, T>&     b,
        size_t                   max_iterations
    ) {
        if (!factored_) {
            const auto L_opt = mat::cholesky(H);
            if (!L_opt) {
                QPResult<NV, NI, T> failed{};
                failed.status = QPStatus::NotPositiveDefinite;
                return failed;
            }
            L_ = L_opt.value();
            factored_ = true;
        }

        const auto result = detail::gi_active_set(H, L_, f, A, b, max_iterations, hint_, hint_count_);

        // Harvest the active set as the next call's hint.
        hint_count_ = 0;
        for (size_t i = 0; i < NI && hint_count_ < KMAX; ++i) {
            if (result.lambda(i) > T{0}) {
                hint_[hint_count_] = i;
                ++hint_count_;
            }
        }
        return result;
    }

    /** @brief Drop the warm start and the cached factorization (call when H or A change) */
    constexpr void reset() {
        factored_ = false;
        hint_count_ = 0;
    }

private:
    Matrix<NV, NV, T>         L_{};    // cached chol(H)
    damp::array<size_t, KMAX> hint_{}; // previous solve's active set
    size_t                    hint_count_{0};
    bool                      factored_{false};
};

/**
 * @struct AdmmSettings
 * @brief Tuning parameters for the ADMM QP solver
 *
 * The defaults follow OSQP's, with tighter tolerances so results agree with
 * the exact active-set solver at test precision: eps_abs and eps_rel default to
 * √ε (ε the machine epsilon of T), 1.5e-8 in double and 3.4e-4 in float, which
 * float iterates can reach. Loosen them (e.g. to 1e-4) when a cheaper
 * moderate-accuracy solve is acceptable.
 *
 * eps_infeasible is relative: each infeasibility certificate test compares a
 * residual with the magnitude of the terms that form it (see AdmmSolver). Its
 * default √ε keeps half the digits of T.
 */
template<typename T = double>
struct AdmmSettings {
    T rho{static_cast<T>(0.1)};             ///< Constraint penalty ρ
    T sigma{static_cast<T>(1.0e-6)};        ///< Primal regularization σ
    T alpha{static_cast<T>(1.6)};           ///< Over-relaxation α ∈ (0, 2)
    T eps_abs{damp::sqrt(eps<T>())};        ///< Absolute residual tolerance
    T eps_rel{damp::sqrt(eps<T>())};        ///< Relative residual tolerance
    T eps_infeasible{damp::sqrt(eps<T>())}; ///< Relative tolerance of the primal and dual infeasibility certificates
};

/**
 * @brief ADMM (OSQP-style) QP solver policy — warm-started, fixed-cost iterations
 *
 * Operator-splitting solver for the same problem as solve_qp():
 * @f$ \min \tfrac{1}{2}x^\top H x + f^\top x @f$ s.t. @f$ Ax \le b @f$. One
 * Cholesky factorization of @f$ K = H + \sigma I + \rho A^\top A @f$ is cached
 * on the first call; each iteration then costs a fixed pair of triangular
 * solves plus a projection — the deterministic-per-tick embedded workhorse.
 * The primal/dual iterates persist across calls (warm start), so a receding-
 * horizon problem that changes little between ticks converges in a handful of
 * iterations.
 *
 * The cached factorization assumes H and A are constant across calls (true for
 * damp::MPC, whose QP matrices are fixed at synthesis); call reset() if they
 * change. Compared with the active-set solver: solutions are accurate to the
 * configured tolerance rather than exact. Because σ > 0 keeps K positive
 * definite, H may be positive semidefinite (H = 0 is an LP).
 *
 * Infeasibility is detected from the successive iterate differences
 * @f$ \delta y^k = y^k - y^{k-1} @f$ and @f$ \delta x^k = x^k - x^{k-1} @f$,
 * which converge on every convex QP; a nonzero limit is a certificate. With
 * @f$ \hat y = \Pi(\delta y) @f$, δy with negative entries and sentinel rows
 * zeroed (the projection onto the polar of the constraint set's recession
 * cone), each iteration that fails the convergence test checks
 * @f[
 *   \text{primal: } \|A^\top \hat y\|_\infty \le \epsilon \textstyle\sum_i \|a_i\|_\infty \hat y_i,
 *   \quad b^\top \hat y < -\epsilon \textstyle\sum_i |b_i| \hat y_i ,
 * @f]
 * @f[
 *   \text{dual: } \|H\delta x\|_\infty \le \epsilon \textstyle\sum_j \|h_j\|_\infty |\delta x_j|,
 *   \quad a_i^\top \delta x \le \epsilon \|a_i\|_\infty \|\delta x\|_\infty,
 *   \quad f^\top \delta x < -\epsilon \textstyle\sum_j |f_j \delta x_j| ,
 * @f]
 * the middle dual test on bounded rows only, with ε = eps_infeasible. The
 * primal test runs first. A primal certificate (no x satisfies Ax ≤ b) is
 * reported as Infeasible. A dual certificate (a direction along which H is
 * singular and the objective decreases without bound on the feasible set) is
 * reported as Unbounded; solve_qp() and solve_qp_interior_point() report
 * NotPositiveDefinite for the same problem, because they cannot factor H. Each
 * test is relative to the constraint rows and the cost, so scaling a row
 * (aᵢ, bᵢ) or the whole problem does not change its outcome for given
 * iterates; the variables are assumed comparably scaled (damp::MPC scales its
 * moves). The iteration itself is not row-scale invariant (one ρ for every
 * row): a badly scaled row slows detection as it slows convergence. A
 * certificate also clears the warm start (the iterates diverge on an
 * infeasible problem) and keeps the factorization.
 *
 * Worst-case time: an iteration costs four passes over A, one over H and the
 * two triangular solves, about 4·NI·NV + 2·NV² multiply-adds. In an iteration
 * that fails the convergence test, the certificates add at most two passes
 * over A and one over H (2·NI·NV + NV²), each behind an O(NI) or O(NV) sign
 * test that usually skips it. Budget the worst case at 1.5× the iteration
 * cost times max_iterations.
 *
 * Sized on the problem (state carries fixed-size iterates):
 * @code
 * using Ctl = MPC<2, 1, 1, 15, 5, 0, double>;
 * MPC<2, 1, 1, 15, 5, 0, double, design::AdmmSolver<Ctl::NZ1, Ctl::NI1, double>> controller{art};
 * @endcode
 *
 * @see solve_qp_admm() for one-shot use, ActiveSetSolver for the exact default
 * @see Stellato et al., "OSQP: an operator splitting solver for quadratic
 *      programs," Mathematical Programming Computation 12, 2020,
 *      https://doi.org/10.1007/s12532-020-00179-2
 * @see "Infeasibility detection in the alternating direction method of
 *      multipliers for convex optimization" (Banjac, Goulart, Stellato & Boyd,
 *      2019), §3.2 Prop. 3.1 (certificates), §5.1 Thm. 5.1 (limits of δx, δy),
 *      §5.2 (termination criteria, Remark 5.3 projection), JOTA 183(2),
 *      https://doi.org/10.1007/s10957-019-01575-y
 */
template<size_t NV, size_t NI, typename T = double>
    requires std::is_floating_point_v<T>
class AdmmSolver {
public:
    constexpr AdmmSolver() = default;
    constexpr explicit AdmmSolver(const AdmmSettings<T>& settings) : settings_(settings) {}

    [[nodiscard]] constexpr QPResult<NV, NI, T> operator()(
        const Matrix<NV, NV, T>& H,
        const ColVec<NV, T>&     f,
        const Matrix<NI, NV, T>& A,
        const ColVec<NI, T>&     b,
        size_t                   max_iterations
    ) {
        constexpr T big = unbounded_bound<T>();
        constexpr T skip_threshold = big / T{4};

        QPResult<NV, NI, T> result{};
        for (size_t i = 0; i < NI; ++i) {
            if (!damp::isfinite(b(i))) {
                result.status = QPStatus::NonFinite; // a NaN bound is no sentinel; -ffast-math would drop its row
                return result;
            }
        }

        if (!factored_) {
            // K = H + σI + ρAᵀA (sentinel rows included — harmless regularization,
            // and it keeps the factorization independent of which bounds are active).
            Matrix<NV, NV, T> K = H;
            for (size_t r = 0; r < NV; ++r) {
                K(r, r) += settings_.sigma;
                for (size_t c = 0; c < NV; ++c) {
                    T acc = T{0};
                    for (size_t i = 0; i < NI; ++i) {
                        acc += A(i, r) * A(i, c);
                    }
                    K(r, c) += settings_.rho * acc;
                }
            }
            const auto L_opt = mat::cholesky(K);
            if (!L_opt) {
                result.status = QPStatus::NotPositiveDefinite;
                return result;
            }
            L_ = L_opt.value();
            factored_ = true;
        }

        const T rho = settings_.rho;
        const T alpha = settings_.alpha;

        size_t   it = 0;
        QPStatus status = QPStatus::MaxIterations;
        for (; it < max_iterations; ++it) {
            // x̃ = K⁻¹(σx − f + Aᵀ(ρz − y))
            ColVec<NV, T> rhs{};
            for (size_t v = 0; v < NV; ++v) {
                T acc = settings_.sigma * x_(v) - f(v);
                for (size_t i = 0; i < NI; ++i) {
                    acc += A(i, v) * (rho * z_(i) - y_(i));
                }
                rhs(v) = acc;
            }
            const ColVec<NV, T> x_tilde = mat::backward_substitute_transpose(L_, mat::forward_substitute(L_, rhs));

            // Relaxed updates with projection onto (−∞, b].
            ColVec<NV, T> x_next{};
            for (size_t v = 0; v < NV; ++v) {
                x_next(v) = alpha * x_tilde(v) + (T{1} - alpha) * x_(v);
            }
            ColVec<NI, T> z_next{};
            ColVec<NI, T> dy{}; // δy = y^k − y^{k−1}
            for (size_t i = 0; i < NI; ++i) {
                T ax = T{0};
                for (size_t v = 0; v < NV; ++v) {
                    ax += A(i, v) * x_tilde(v);
                }
                const T v_relaxed = (alpha * ax) + ((T{1} - alpha) * z_(i));
                const T candidate = v_relaxed + (y_(i) / rho);
                if (b(i) >= skip_threshold || candidate < b(i)) {
                    z_next(i) = candidate;
                } else {
                    z_next(i) = b(i);
                }
                dy(i) = rho * (v_relaxed - z_next(i));
                y_(i) += dy(i);
            }
            const ColVec<NV, T> dx = x_next - x_; // δx = x^k − x^{k−1}
            x_ = x_next;
            z_ = z_next;

            // Termination: ‖Ax − z‖∞ and ‖Hx + f + Aᵀy‖∞ under (abs + rel) tolerances.
            // The maxima drop a NaN, so the sum of the residuals carries it.
            T r_sum = T{0};
            T r_prim = T{0};
            T ax_norm = T{0};
            T z_norm = T{0};
            for (size_t i = 0; i < NI; ++i) {
                T ax = T{0};
                for (size_t v = 0; v < NV; ++v) {
                    ax += A(i, v) * x_(v);
                }
                const T rp = damp::abs(ax - z_(i));
                r_sum += rp;
                r_prim = (rp > r_prim) ? rp : r_prim;
                ax_norm = (damp::abs(ax) > ax_norm) ? damp::abs(ax) : ax_norm;
                z_norm = (damp::abs(z_(i)) > z_norm) ? damp::abs(z_(i)) : z_norm;
            }
            T r_dual = T{0};
            T scale_dual = T{0};
            for (size_t v = 0; v < NV; ++v) {
                T hx = f(v);
                for (size_t c = 0; c < NV; ++c) {
                    hx += H(v, c) * x_(c);
                }
                T aty = T{0};
                for (size_t i = 0; i < NI; ++i) {
                    aty += A(i, v) * y_(i);
                }
                const T rd = damp::abs(hx + aty);
                r_sum += rd;
                r_dual = (rd > r_dual) ? rd : r_dual;
                const T sd = damp::abs(hx - f(v)) + damp::abs(aty) + damp::abs(f(v));
                scale_dual = (sd > scale_dual) ? sd : scale_dual;
            }
            const T eps_prim = settings_.eps_abs + settings_.eps_rel * ((ax_norm > z_norm) ? ax_norm : z_norm);
            const T eps_dual = settings_.eps_abs + settings_.eps_rel * scale_dual;
            if (!damp::isfinite(r_sum)) {
                status = QPStatus::NonFinite; // a non-finite H, f, A or b
            } else if (r_prim <= eps_prim && r_dual <= eps_dual) {
                status = QPStatus::Success;
            } else if (primal_infeasible(A, b, dy)) {
                status = QPStatus::Infeasible;
            } else if (dual_infeasible(H, f, A, b, dx)) {
                status = QPStatus::Unbounded;
            }
            if (status != QPStatus::MaxIterations) {
                ++it;
                break;
            }
        }

        result.x = x_;
        for (size_t i = 0; i < NI; ++i) {
            result.lambda(i) = (y_(i) > T{0}) ? y_(i) : T{0};
        }
        result.iterations = it;
        result.status = status;
        result.success = (status == QPStatus::Success);
        T obj = T{0};
        for (size_t r = 0; r < NV; ++r) {
            T hx = T{0};
            for (size_t c = 0; c < NV; ++c) {
                hx += H(r, c) * x_(c);
            }
            obj += x_(r) * (static_cast<T>(0.5) * hx + f(r));
        }
        result.objective = obj;
        if (status == QPStatus::Infeasible || status == QPStatus::Unbounded || status == QPStatus::NonFinite) {
            clear_iterates(); // diverging or non-finite iterates are no warm start for the next call
        }
        return result;
    }

    /** @brief Drop the warm start and the cached factorization (call when H or A change) */
    constexpr void reset() {
        clear_iterates();
        factored_ = false;
    }

private:
    constexpr void clear_iterates() {
        x_ = ColVec<NV, T>{};
        z_ = ColVec<NI, T>{};
        y_ = ColVec<NI, T>{};
    }

    /**
     * @brief Primal infeasibility certificate from δy (Banjac et al. 2019, (30) with Remark 5.3)
     *
     * ŷ is δy projected onto the polar of the recession cone of {z : z ≤ b}: negative
     * entries and sentinel rows zeroed. True when Aᵀŷ cancels and the support bᵀŷ is
     * negative, each relative to the terms that form it. The O(NI) support test runs first.
     */
    [[nodiscard]] constexpr bool primal_infeasible(const Matrix<NI, NV, T>& A, const ColVec<NI, T>& b, const ColVec<NI, T>& dy) const {
        constexpr T   skip_threshold = unbounded_bound<T>() / T{4};
        const T       tol = settings_.eps_infeasible;
        ColVec<NI, T> y_hat{};
        T             support = T{0};
        T             support_scale = T{0};
        for (size_t i = 0; i < NI; ++i) {
            if (b(i) < skip_threshold && dy(i) > T{0}) {
                y_hat(i) = dy(i);
                support += b(i) * dy(i);
                support_scale += damp::abs(b(i)) * dy(i);
            }
        }
        if (!(support < -(tol * support_scale))) {
            return false;
        }
        ColVec<NV, T> aty{};
        T             aty_scale = T{0}; // Σᵢ ‖aᵢ‖∞ ŷᵢ
        for (size_t i = 0; i < NI; ++i) {
            T row_norm = T{0};
            for (size_t v = 0; v < NV; ++v) {
                aty(v) += A(i, v) * y_hat(i);
                row_norm = damp::max(row_norm, damp::abs(A(i, v)));
            }
            aty_scale += row_norm * y_hat(i);
        }
        for (size_t v = 0; v < NV; ++v) {
            if (damp::abs(aty(v)) > tol * aty_scale) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Dual infeasibility certificate from δx (Banjac et al. 2019, (31))
     *
     * True when δx is a descent direction (fᵀδx < 0) in the null space of H that stays in
     * the recession cone of the constraints (aᵢᵀδx ≤ 0 on bounded rows), each relative to
     * the terms that form it. The O(NV) descent test runs first. H is taken symmetric, so
     * its row norms are its column norms.
     */
    [[nodiscard]] constexpr bool dual_infeasible(
        const Matrix<NV, NV, T>& H, const ColVec<NV, T>& f, const Matrix<NI, NV, T>& A, const ColVec<NI, T>& b, const ColVec<NV, T>& dx
    ) const {
        constexpr T skip_threshold = unbounded_bound<T>() / T{4};
        const T     tol = settings_.eps_infeasible;
        T           descent = T{0};
        T           descent_scale = T{0};
        T           dx_norm = T{0};
        for (size_t v = 0; v < NV; ++v) {
            descent += f(v) * dx(v);
            descent_scale += damp::abs(f(v) * dx(v));
            dx_norm = damp::max(dx_norm, damp::abs(dx(v)));
        }
        if (!(descent < -(tol * descent_scale))) {
            return false;
        }
        T hdx_norm = T{0};
        T hdx_scale = T{0}; // Σⱼ ‖hⱼ‖∞ |δxⱼ|
        for (size_t r = 0; r < NV; ++r) {
            T hdx = T{0};
            T row_norm = T{0};
            for (size_t c = 0; c < NV; ++c) {
                hdx += H(r, c) * dx(c);
                row_norm = damp::max(row_norm, damp::abs(H(r, c)));
            }
            hdx_norm = damp::max(hdx_norm, damp::abs(hdx));
            hdx_scale += row_norm * damp::abs(dx(r));
        }
        if (hdx_norm > tol * hdx_scale) {
            return false;
        }
        for (size_t i = 0; i < NI; ++i) {
            if (b(i) >= skip_threshold) {
                continue; // a free row bounds no direction
            }
            T adx = T{0};
            T row_norm = T{0};
            for (size_t v = 0; v < NV; ++v) {
                adx += A(i, v) * dx(v);
                row_norm = damp::max(row_norm, damp::abs(A(i, v)));
            }
            if (adx > tol * row_norm * dx_norm) {
                return false;
            }
        }
        return true;
    }

    Matrix<NV, NV, T> L_{}; // Cholesky factor of H + σI + ρAᵀA
    ColVec<NV, T>     x_{}; // warm-started primal iterate
    ColVec<NI, T>     z_{}; // warm-started auxiliary iterate
    ColVec<NI, T>     y_{}; // warm-started dual iterate
    AdmmSettings<T>   settings_{};
    bool              factored_{false};
};

/**
 * @brief One-shot ADMM QP solve (cold start)
 *
 * Convenience wrapper constructing a local AdmmSolver. For receding-horizon
 * use, hold an AdmmSolver instance instead and keep its warm start.
 *
 * @see AdmmSolver, solve_qp
 */
template<size_t NV, size_t NI, typename T = double>
[[nodiscard]] constexpr QPResult<NV, NI, T> solve_qp_admm(
    const Matrix<NV, NV, T>& H,
    const ColVec<NV, T>&     f,
    const Matrix<NI, NV, T>& A,
    const ColVec<NI, T>&     b,
    size_t                   max_iterations = 4000,
    const AdmmSettings<T>&   settings = {}
) {
    AdmmSolver<NV, NI, T> solver{settings};
    return solver(H, f, A, b, max_iterations);
}

/**
 * @brief Interior-point QP solver (primal-dual path following)
 *
 * Solves the same problem as solve_qp() by damped Newton steps on the
 * perturbed KKT system with slacks s > 0 and multipliers λ > 0:
 * @f[
 *   Hx + f + A^\top\lambda = 0,\quad Ax + s = b,\quad s_i\lambda_i = \sigma\mu,
 * @f]
 * with fixed centering σ = 0.1 and fraction-to-boundary step τ = 0.995. Each
 * iteration solves one NV×NV SPD system @f$ (H + A^\top \mathrm{diag}(\lambda/s) A)\,\Delta x @f$;
 * convergence is superlinear near the solution and typically takes 10–30
 * iterations regardless of how many constraints are active — the preferred
 * profile when many constraints bind at once. Iterates approach feasibility
 * from the first few steps, so a budget-exhausted stop is near-feasible
 * (unlike the dual active-set method). Infeasible problems are detected by
 * step-length collapse (heuristic, not a certificate).
 *
 * @note Compare with MATLAB®'s mpcInteriorPointSolver / quadprog
 *       'interior-point-convex'.
 * @see solve_qp, InteriorPointSolver
 * @see Nocedal & Wright, "Numerical Optimization," 2nd ed., 2006, ch. 16.6
 *      (predictor-corrector variants) and ch. 19 (interior methods)
 *
 * @param H               Hessian (NV × NV, symmetric positive definite)
 * @param f               Linear cost term (NV)
 * @param A               Inequality constraint matrix (NI × NV), rows aᵢᵀx ≤ bᵢ
 * @param b               Inequality bounds; sentinel rows (≥ max/4) are skipped
 * @param max_iterations  Newton-step budget (default 50)
 * @return QPResult with minimizer, multipliers, objective, and status
 */
template<size_t NV, size_t NI, typename T = double>
    requires std::is_floating_point_v<T>
[[nodiscard]] constexpr QPResult<NV, NI, T> solve_qp_interior_point(
    const Matrix<NV, NV, T>& H,
    const ColVec<NV, T>&     f,
    const Matrix<NI, NV, T>& A,
    const ColVec<NI, T>&     b,
    size_t                   max_iterations = 50
) {
    constexpr T big = unbounded_bound<T>();
    constexpr T skip_threshold = big / T{4};
    constexpr T sigma = static_cast<T>(0.1); // centering
    constexpr T tau = static_cast<T>(0.995); // fraction to boundary
    // Residuals are tested relative to the terms that form them (scale-invariant).
    const T tol = damp::sqrt(eps<T>());

    QPResult<NV, NI, T> result{};

    // Active-row mask (sentinel rows are absent from the problem; a NaN bound is not one).
    damp::array<bool, NI> active{};
    size_t                m = 0;
    for (size_t i = 0; i < NI; ++i) {
        if (!damp::isfinite(b(i))) {
            result.status = QPStatus::NonFinite;
            return result;
        }
        active[i] = (b(i) < skip_threshold);
        m += active[i] ? 1U : 0U;
    }

    const auto L_opt = mat::cholesky(H);
    if (!L_opt) {
        result.status = QPStatus::NotPositiveDefinite;
        return result;
    }
    ColVec<NV, T> x = -mat::backward_substitute_transpose(L_opt.value(), mat::forward_substitute(L_opt.value(), f));
    for (size_t v = 0; v < NV; ++v) {
        if (!damp::isfinite(x(v))) {
            result.status = QPStatus::NonFinite; // a non-finite H or f
            return result;
        }
    }

    if (m == 0) {
        result.x = x;
        result.status = QPStatus::Success;
        result.success = true;
        T obj = T{0};
        for (size_t r = 0; r < NV; ++r) {
            T hx = T{0};
            for (size_t c = 0; c < NV; ++c) {
                hx += H(r, c) * x(c);
            }
            obj += x(r) * (static_cast<T>(0.5) * hx + f(r));
        }
        result.objective = obj;
        return result;
    }

    // Strictly interior start.
    ColVec<NI, T> s{};
    ColVec<NI, T> lam{};
    for (size_t i = 0; i < NI; ++i) {
        if (!active[i]) {
            continue;
        }
        T slack = b(i);
        for (size_t v = 0; v < NV; ++v) {
            slack -= A(i, v) * x(v);
        }
        s(i) = (slack > T{1}) ? slack : T{1};
        lam(i) = T{1};
    }

    size_t stalls = 0;
    for (size_t it = 0; it < max_iterations; ++it) {
        // Residuals and duality measure.
        ColVec<NV, T> r_d = f;
        bool          dual_ok = true;
        T             gap_scale = T{0}; // magnitude of the objective and of λᵀ(b − Ax)
        for (size_t v = 0; v < NV; ++v) {
            T mag = damp::abs(f(v));
            T hx = T{0};
            for (size_t c = 0; c < NV; ++c) {
                hx += H(v, c) * x(c);
                mag += damp::abs(H(v, c) * x(c));
            }
            r_d(v) += hx;
            for (size_t i = 0; i < NI; ++i) {
                if (active[i]) {
                    r_d(v) += A(i, v) * lam(i);
                    mag += damp::abs(A(i, v) * lam(i));
                }
            }
            dual_ok = dual_ok && (damp::abs(r_d(v)) <= tol * mag);
            gap_scale += damp::abs(x(v) * hx) + damp::abs(f(v) * x(v));
        }
        ColVec<NI, T> r_p{};
        T             mu = T{0};
        bool          primal_ok = true;
        for (size_t i = 0; i < NI; ++i) {
            if (!active[i]) {
                continue;
            }
            T ax = T{0};
            T mag = damp::abs(s(i)) + damp::abs(b(i));
            for (size_t v = 0; v < NV; ++v) {
                ax += A(i, v) * x(v);
                mag += damp::abs(A(i, v) * x(v));
            }
            r_p(i) = ax + s(i) - b(i);
            primal_ok = primal_ok && (damp::abs(r_p(i)) <= tol * mag);
            gap_scale += lam(i) * mag; // λᵢ(|bᵢ| + |aᵢᵀx| + sᵢ) bounds the gap term λᵢsᵢ
            mu += s(i) * lam(i);
        }
        mu /= static_cast<T>(m);

        if (dual_ok && primal_ok && (static_cast<T>(m) * mu <= tol * gap_scale)) {
            result.status = QPStatus::Success;
            result.success = true;
            break;
        }

        // Newton system on the condensed normal equations.
        Matrix<NV, NV, T> M = H;
        ColVec<NV, T>     rhs{};
        for (size_t v = 0; v < NV; ++v) {
            rhs(v) = -r_d(v);
        }
        for (size_t i = 0; i < NI; ++i) {
            if (!active[i]) {
                continue;
            }
            const T ratio = lam(i) / s(i);
            for (size_t r = 0; r < NV; ++r) {
                for (size_t c = 0; c < NV; ++c) {
                    M(r, c) += ratio * A(i, r) * A(i, c);
                }
            }
            // r_c_i = s_i λ_i − σμ;  rhs −= aᵢ (λ_i r_p_i − r_c_i)/s_i
            const T coeff = (lam(i) * r_p(i) - (s(i) * lam(i) - sigma * mu)) / s(i);
            for (size_t v = 0; v < NV; ++v) {
                rhs(v) -= A(i, v) * coeff;
            }
        }
        const auto M_chol = mat::cholesky(M);
        if (!M_chol) {
            result.status = QPStatus::Degenerate;
            break;
        }
        const ColVec<NV, T> dx = mat::backward_substitute_transpose(M_chol.value(), mat::forward_substitute(M_chol.value(), rhs));

        // Recover Δs, Δλ and the fraction-to-boundary step length.
        ColVec<NI, T> ds{};
        ColVec<NI, T> dlam{};
        T             step = T{1};
        for (size_t i = 0; i < NI; ++i) {
            if (!active[i]) {
                continue;
            }
            T adx = T{0};
            for (size_t v = 0; v < NV; ++v) {
                adx += A(i, v) * dx(v);
            }
            ds(i) = -r_p(i) - adx;
            dlam(i) = (-(s(i) * lam(i) - sigma * mu) - lam(i) * ds(i)) / s(i);
            if (ds(i) < T{0}) {
                const T limit = -tau * s(i) / ds(i);
                step = (limit < step) ? limit : step;
            }
            if (dlam(i) < T{0}) {
                const T limit = -tau * lam(i) / dlam(i);
                step = (limit < step) ? limit : step;
            }
        }

        if (step < static_cast<T>(1.0e-10)) {
            ++stalls;
            if (stalls >= 3) {
                result.status = QPStatus::Infeasible;
                break;
            }
        } else {
            stalls = 0;
        }

        for (size_t v = 0; v < NV; ++v) {
            x(v) += step * dx(v);
        }
        for (size_t i = 0; i < NI; ++i) {
            if (active[i]) {
                s(i) += step * ds(i);
                lam(i) += step * dlam(i);
            }
        }
        result.iterations = it + 1;
    }

    for (size_t v = 0; v < NV; ++v) {
        if (!damp::isfinite(x(v))) {
            result.status = QPStatus::NonFinite; // a non-finite row of A
            result.success = false;
        }
    }
    result.x = x;
    for (size_t i = 0; i < NI; ++i) {
        result.lambda(i) = active[i] ? lam(i) : T{0};
    }
    T obj = T{0};
    for (size_t r = 0; r < NV; ++r) {
        T hx = T{0};
        for (size_t c = 0; c < NV; ++c) {
            hx += H(r, c) * x(c);
        }
        obj += x(r) * (static_cast<T>(0.5) * hx + f(r));
    }
    result.objective = obj;
    return result;
}

/**
 * @brief Interior-point QP solver policy for damp::MPC
 *
 * Stateless (interior-point methods do not warm start well); each call runs
 * solve_qp_interior_point(). Prefer this over the active-set default when many
 * constraints bind at once, or when a budget-exhausted stop must still be
 * near-feasible.
 *
 * @see solve_qp_interior_point, ActiveSetSolver, AdmmSolver
 */
struct InteriorPointSolver {
    template<size_t NV, size_t NI, typename T>
    [[nodiscard]] constexpr QPResult<NV, NI, T> operator()(
        const Matrix<NV, NV, T>& H,
        const ColVec<NV, T>&     f,
        const Matrix<NI, NV, T>& A,
        const ColVec<NI, T>&     b,
        size_t                   max_iterations
    ) const {
        return solve_qp_interior_point(H, f, A, b, (max_iterations < 50) ? max_iterations : size_t{50});
    }
};

/**
 * @brief Mixed-integer QP by branch and bound over the active-set relaxation
 *
 * Solves @f$ \min \tfrac{1}{2}x^\top H x + f^\top x @f$ s.t. @f$ Ax \le b @f$
 * with selected variables restricted to integers. Depth-first branch and
 * bound: each node's continuous relaxation (solve_qp with the node's integer
 * bounds appended as box rows) provides the lower bound; branching splits the
 * most fractional integer variable into floor/ceil children; nodes whose
 * relaxation cannot beat the incumbent are pruned.
 *
 * The node stack is fixed-size (MaxNodes); if it overflows or the node budget
 * is exhausted, the best integer-feasible solution found so far is returned
 * with status MaxIterations (anytime behaviour). `iterations` reports the
 * number of relaxations solved.
 *
 * For integer manipulated variables in damp::MPC, restrict the Δu decision
 * variables through a custom solver policy wrapping this function, and keep
 * the input scale factors at 1 (integrality is not preserved by rescaling) —
 * i.e. unbounded u limits or a constraint span of exactly 1.
 *
 * @note Compare with MATLAB®'s mpc integer MV support (BranchBound options).
 * @see solve_qp
 * @see Floudas, "Nonlinear and Mixed-Integer Optimization," Oxford, 1995, ch. 5
 *      (branch and bound); Bemporad & Morari, "Control of systems integrating
 *      logic, dynamics, and constraints," Automatica 35(3), 1999 (hybrid MPC)
 *
 * @tparam MaxNodes   Node-stack capacity (compile-time; default 64)
 * @param H           Hessian (NV × NV, symmetric positive definite)
 * @param f           Linear cost term (NV)
 * @param A           Inequality constraint matrix (NI × NV)
 * @param b           Inequality bounds; sentinel rows are skipped
 * @param is_integer  Per-variable integrality mask
 * @param max_nodes   Relaxation budget (default: explore up to 16·NV nodes)
 * @param integer_tol Distance from an integer treated as integral, plus 8·(NV+1)·ε of the
 *                    value (the relaxation's roundoff at that magnitude), capped at 0.25
 * @return QPResult over the original variables (lambda reports the incumbent's
 *         relaxation multipliers for the original NI rows); NonFinite for non-finite data
 */
template<size_t MaxNodes = 64, size_t NV, size_t NI, typename T>
    requires std::is_floating_point_v<T>
[[nodiscard]] constexpr QPResult<NV, NI, T> solve_miqp(
    const Matrix<NV, NV, T>&     H,
    const ColVec<NV, T>&         f,
    const Matrix<NI, NV, T>&     A,
    const ColVec<NI, T>&         b,
    const damp::array<bool, NV>& is_integer,
    size_t                       max_nodes = 16 * NV,
    T                            integer_tol = static_cast<T>(1.0e-6)
) {
    constexpr T      big = unbounded_bound<T>();
    constexpr size_t NIE = NI + (2 * NV); // relaxation rows: original + integer box bounds
    // Integral within integer_tol plus the roundoff of the value (a float relaxation on a bound
    // is off by its ulp). The cap keeps a half-integer fractional at any magnitude.
    constexpr T roundoff = static_cast<T>(8 * (NV + 1)) * eps<T>();
    const auto  fraction = [&](T value) {
        const T frac = damp::abs(value - detail::round_scalar(value));
        return (frac > integer_tol + damp::min(roundoff * damp::abs(value), T{0.25})) ? frac : T{0};
    };

    QPResult<NV, NI, T> result{};

    //! The search runs about x₀, the root relaxation rounded on the integer variables: each
    //! relaxation then sees an integer bound at the scale of its offset from x₀ rather than of x.
    //! In float a bound at 1e4 lies inside the active-set solver's slack tolerance ((NV+1)·64·ε of
    //! the row's magnitude) and would never bind.
    const auto root = solve_qp(H, f, A, b);
    if (root.status == QPStatus::NonFinite || root.status == QPStatus::Infeasible) {
        result.status = root.status;
        result.iterations = 1;
        return result;
    }
    ColVec<NV, T> x0{};
    bool          integral = root.success;
    for (size_t v = 0; v < NV; ++v) {
        if (is_integer[v] && root.success) {
            x0(v) = detail::round_scalar(root.x(v));
            integral = integral && !(fraction(root.x(v)) > T{0});
        }
    }
    if (integral) {
        QPResult<NV, NI, T> out = root;
        out.iterations = 1; // the root relaxation is integral
        return out;
    }
    const ColVec<NV, T> fs{f + (H * x0)};
    ColVec<NI, T>       bs{};
    for (size_t i = 0; i < NI; ++i) {
        T ax = T{0};
        for (size_t v = 0; v < NV; ++v) {
            ax += A(i, v) * x0(v);
        }
        bs(i) = (b(i) >= big / T{4}) ? b(i) : (b(i) - ax);
    }
    result.iterations = 1;

    // Extended constraint matrix: [A; I; −I]. Bound rows are sentinel-disabled
    // until branching tightens them.
    Matrix<NIE, NV, T> A_ext{};
    for (size_t i = 0; i < NI; ++i) {
        for (size_t v = 0; v < NV; ++v) {
            A_ext(i, v) = A(i, v);
        }
    }
    for (size_t v = 0; v < NV; ++v) {
        A_ext(NI + v, v) = T{1};
        A_ext(NI + NV + v, v) = T{-1};
    }

    struct Node {
        ColVec<NV, T> lb{};
        ColVec<NV, T> ub{};
    };
    damp::array<Node, MaxNodes> stack{};
    for (size_t v = 0; v < NV; ++v) {
        stack[0].lb(v) = -big;
        stack[0].ub(v) = big;
    }
    size_t depth = 1;

    T    best_obj = big;
    bool truncated = false;

    while (depth > 0 && result.iterations < max_nodes) {
        const Node node = stack[--depth];

        ColVec<NIE, T> b_ext{};
        for (size_t i = 0; i < NI; ++i) {
            b_ext(i) = bs(i);
        }
        for (size_t v = 0; v < NV; ++v) {
            b_ext(NI + v) = node.ub(v);
            b_ext(NI + NV + v) = (node.lb(v) <= -big / T{4}) ? big : -node.lb(v);
        }

        const auto relaxed = solve_qp(H, fs, A_ext, b_ext);
        ++result.iterations;
        if (relaxed.status == QPStatus::NonFinite) {
            const size_t solved = result.iterations;
            result = QPResult<NV, NI, T>{};
            result.iterations = solved;
            result.status = QPStatus::NonFinite;
            return result;
        }
        if (!relaxed.success) {
            // Only a proven-infeasible relaxation prunes soundly; a budget or
            // degeneracy stop leaves the subtree unexplored, so the result can no
            // longer claim optimality (or infeasibility).
            truncated = truncated || (relaxed.status != QPStatus::Infeasible);
            continue;
        }
        if (relaxed.objective >= best_obj) {
            continue; // bound — cannot beat the incumbent
        }

        // Most fractional integer variable (offsets from the integral x₀ keep integrality).
        size_t branch_var = NV;
        T      worst_frac = T{0};
        for (size_t v = 0; v < NV; ++v) {
            if (!is_integer[v]) {
                continue;
            }
            const T frac = fraction(relaxed.x(v));
            if (frac > worst_frac) {
                worst_frac = frac;
                branch_var = v;
            }
        }

        if (branch_var == NV) {
            // Integer feasible: new incumbent.
            best_obj = relaxed.objective;
            result.x = relaxed.x;
            for (size_t i = 0; i < NI; ++i) {
                result.lambda(i) = relaxed.lambda(i);
            }
            result.objective = relaxed.objective;
            result.success = true;
            continue;
        }

        if (depth + 2 > MaxNodes) {
            truncated = true;
            continue;
        }
        const T value = relaxed.x(branch_var);
        const T fl = detail::floor_scalar(value);
        // Down child (x ≤ floor) and up child (x ≥ ceil); explore the nearer first.
        Node down = node;
        down.ub(branch_var) = fl;
        Node up = node;
        up.lb(branch_var) = fl + T{1};
        if (value - fl < static_cast<T>(0.5)) {
            stack[depth++] = up;
            stack[depth++] = down; // popped first
        } else {
            stack[depth++] = down;
            stack[depth++] = up;
        }
    }

    truncated = truncated || (depth > 0);
    if (result.success) {
        // Anytime: the incumbent is a valid integer-feasible point even when
        // the search was truncated (it just may not be proven optimal).
        result.status = truncated ? QPStatus::MaxIterations : QPStatus::Success;
        result.x = result.x + x0; // back from the offsets
        T obj = T{0};
        for (size_t r = 0; r < NV; ++r) {
            T hx = T{0};
            for (size_t c = 0; c < NV; ++c) {
                hx += H(r, c) * result.x(c);
            }
            obj += result.x(r) * ((static_cast<T>(0.5) * hx) + f(r));
        }
        result.objective = obj;
    } else if (truncated) {
        result.status = QPStatus::MaxIterations; // budget ran out before any integer point was found
    } else {
        result.status = QPStatus::Infeasible;
    }
    return result;
}

} // namespace design

} // namespace damp
