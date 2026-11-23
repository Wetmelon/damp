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


#include "damp/backend.hpp"
#include "damp/math/math.hpp"
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
    LimitsRefused,       ///< (MPC runtime) A limit assigned since the last tick was refused, not a box or moved under an InvariantSet: the tick solved on the limits in force
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
    return damp::numeric_limits<T>::max();
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
struct [[nodiscard]] QPResult {
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
    constexpr T  hi = damp::numeric_limits<T>::max() / T{4};
    ColVec<N, U> out{};
    for (size_t i = 0; i < N; ++i) {
        if (v(i) >= hi) {
            out(i) = damp::numeric_limits<U>::max();
        } else if (v(i) <= -hi) {
            out(i) = -damp::numeric_limits<U>::max();
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
    if (!damp::isfinite(v) || !(damp::abs(v) < T{1} / eps<T>())) {
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
 * @brief Floor on the sine of the angle between a constraint normal and the span of the active
 *        ones, in the metric of H⁻¹: at or below it the normal counts as dependent on them.
 *
 * 8·(KMAX+1)·ε. Rotating a normal ñ into the basis of the active set leaves a roundoff of at
 * most about (KMAX+1)·ε·‖ñ‖ in each coordinate of its free part, and far less in practice: the
 * dependent rows of a soft-constrained MPC (11 variables, single precision) read 6ε, a sixteenth
 * of the floor. The floor decides between a primal step and a dual one; the point of a full
 * step is recomputed from the factors and does not use the length it tested (gi_active_set).
 * Relative to ‖ñ‖, so it does not change with the scale of H, f, A or b.
 */
template<size_t KMAX, typename T>
[[nodiscard]] constexpr T dependence_tol() {
    return static_cast<T>(KMAX + 1) * T{8} * eps<T>();
}

/**
 * @brief Orthogonal factorization of the active constraint normals in the metric of H⁻¹
 *
 * With H = LLᵀ and y = Lᵀx, a constraint nᵀx ≥ β reads ñᵀy ≥ β with ñ = L⁻¹n, and
 * nᵢᵀH⁻¹nⱼ = ñᵢᵀñⱼ: in y the QP is a least-distance problem. The k active ñ are the columns of
 * @f[
 *   \tilde N = Q \begin{bmatrix} R \\ 0 \end{bmatrix},
 * @f]
 * Q orthogonal and R upper triangular: the first k columns of Q span the active normals, the
 * others the directions still free. For a candidate ñ, c = Qᵀñ splits into its coefficients on
 * the active normals (r = R⁻¹c₁) and its free part c₂, whose length is ‖ñ‖ times the sine of the
 * angle to their span. A normal joins by one Householder reflection of the free columns and
 * leaves by a sweep of Givens rotations, O(NV²) each. Q is kept apart from L (the paper holds
 * J = L⁻ᵀQ): a normal costs one triangular solve, and no inverse of L is formed.
 *
 * Nothing is formed from the Gram matrix ÑᵀÑ = NᵀH⁻¹N. Its Cholesky factor sees the square of
 * that sine, so it tells a dependent normal from an independent one only down to an angle of
 * √ε, 3·10⁻⁴ in single precision; the orthogonal factors resolve ε.
 *
 * @see "A numerically stable dual method for solving strictly convex quadratic programs"
 *      (Goldfarb & Idnani, 1983), Mathematical Programming 27: the factors J = L⁻ᵀQ and R of
 *      the active set, and their updates when a constraint is added or dropped
 * @see "Methods for Modifying Matrix Factorizations" (Gill, Golub, Murray & Saunders, 1974),
 *      Mathematics of Computation 28(126): a column added to and removed from a QR factorization
 */
template<size_t NV, size_t KMAX, typename T>
struct ActiveNormals {
    Matrix<NV, NV, T>     Q = Matrix<NV, NV, T>::identity(); ///< Orthogonal: columns 0 … k−1 span the active normals
    Matrix<KMAX, KMAX, T> R{};                               ///< Upper triangular on its leading k × k block
    size_t                k{0};                              ///< Active normals, ≤ KMAX

    /** @brief c = Qᵀv: the first k entries along the active normals, the rest in the free directions */
    [[nodiscard]] constexpr ColVec<NV, T> rotate(const ColVec<NV, T>& v) const {
        ColVec<NV, T> c{};
        for (size_t j = 0; j < NV; ++j) {
            T dot{0};
            for (size_t i = 0; i < NV; ++i) {
                dot += Q(i, j) * v(i);
            }
            c(j) = dot;
        }
        return c;
    }

    /** @brief ‖c₂‖²: the squared length of the free part of a rotated normal, never negative */
    [[nodiscard]] constexpr T free_norm_sq(const ColVec<NV, T>& c) const {
        T sum{0};
        for (size_t j = k; j < NV; ++j) {
            sum += c(j) * c(j);
        }
        return sum;
    }

    /** @brief Q₂c₂: the free part of a rotated normal, back in y (the primal step direction LᵀΔx) */
    [[nodiscard]] constexpr ColVec<NV, T> free_part(const ColVec<NV, T>& c) const {
        ColVec<NV, T> z{};
        for (size_t i = 0; i < NV; ++i) {
            T sum{0};
            for (size_t j = k; j < NV; ++j) {
                sum += Q(i, j) * c(j);
            }
            z(i) = sum;
        }
        return z;
    }

    /** @brief r = R⁻¹c₁: the least-squares coefficients of a rotated normal on the active ones */
    [[nodiscard]] constexpr ColVec<KMAX, T> coefficients(const ColVec<NV, T>& c) const {
        const size_t    n = (k < KMAX) ? k : KMAX; // provable index bound for the optimizer
        ColVec<KMAX, T> r{};
        for (size_t ii = n; ii > 0; --ii) {
            const size_t i = ii - 1;
            T            sum = c(i);
            for (size_t j = i + 1; j < n; ++j) {
                sum -= R(i, j) * r(j);
            }
            r(i) = sum / R(i, i);
        }
        return r;
    }

    /**
     * @brief The normal with rotated coordinates @p c joins the active set
     *
     * One Householder reflection P = I − u·uᵀ/(‖c₂‖·(‖c₂‖ + |c_k|)), u = c₂ − σ·e₁, of the free
     * columns sends c₂ to (σ, 0, …), σ = ∓‖c₂‖ with the sign that adds in u's first entry. The new
     * column of R is (c₁, σ). The caller has checked that c₂ is not zero.
     */
    constexpr void append(const ColVec<NV, T>& c) {
        if (k >= KMAX) {
            return; // the column would fall outside the fixed-size factor
        }
        const T       len = damp::sqrt(free_norm_sq(c));
        const T       alpha = c(k);
        const T       sigma = (alpha > T{0}) ? -len : len;
        const T       denom = len * (len + damp::abs(alpha));
        ColVec<NV, T> u{};
        for (size_t j = k + 1; j < NV; ++j) {
            u(j) = c(j);
        }
        u(k) = alpha - sigma;
        if (denom > T{0}) {
            for (size_t i = 0; i < NV; ++i) {
                T w{0};
                for (size_t j = k; j < NV; ++j) {
                    w += Q(i, j) * u(j);
                }
                w /= denom;
                for (size_t j = k; j < NV; ++j) {
                    Q(i, j) -= w * u(j);
                }
            }
        }
        for (size_t i = 0; i < k; ++i) {
            R(i, k) = c(i);
        }
        R(k, k) = sigma;
        ++k;
    }

    /**
     * @brief Active normal @p drop leaves: the columns after it shift left, and one Givens rotation
     *        per shifted column restores the triangle. @p c, a rotated candidate, turns with Q.
     */
    constexpr void remove(size_t drop, ColVec<NV, T>& c) {
        if (k == 0 || k > KMAX || drop >= k) {
            return;
        }
        const size_t last = k - 1;
        for (size_t j = drop; j < last; ++j) {
            for (size_t i = 0; i <= j + 1; ++i) {
                R(i, j) = R(i, j + 1);
            }
        }
        for (size_t j = drop; j < last; ++j) {
            const T h = damp::hypot(R(j, j), R(j + 1, j));
            if (!(h > T{0})) {
                continue; // nothing below the diagonal to remove
            }
            const T cs = R(j, j) / h;
            const T sn = R(j + 1, j) / h;
            R(j, j) = h;
            R(j + 1, j) = T{0};
            for (size_t col = j + 1; col < last; ++col) {
                const T top = R(j, col);
                const T bottom = R(j + 1, col);
                R(j, col) = (cs * top) + (sn * bottom);
                R(j + 1, col) = (cs * bottom) - (sn * top);
            }
            for (size_t i = 0; i < NV; ++i) {
                const T left = Q(i, j);
                const T right = Q(i, j + 1);
                Q(i, j) = (cs * left) + (sn * right);
                Q(i, j + 1) = (cs * right) - (sn * left);
            }
            const T upper = c(j);
            const T lower = c(j + 1);
            c(j) = (cs * upper) + (sn * lower);
            c(j + 1) = (cs * lower) - (sn * upper);
        }
        k = last;
    }

    /**
     * @brief Q₁·R⁻ᵀ·rhs: the shortest y with ñₐᵀy = rhsₐ on every active normal
     *
     * The forward substitution is backward stable: the active constraints are met to the
     * roundoff of the terms of Rᵀ(R⁻ᵀ·rhs), however close to dependent the normals are.
     */
    [[nodiscard]] constexpr ColVec<NV, T> range_point(const ColVec<KMAX, T>& rhs) const {
        const size_t    n = (k < KMAX) ? k : KMAX; // provable index bound for the optimizer
        ColVec<KMAX, T> w{};
        for (size_t i = 0; i < n; ++i) {
            T sum = rhs(i);
            for (size_t j = 0; j < i; ++j) {
                sum -= R(j, i) * w(j);
            }
            w(i) = sum / R(i, i);
        }
        ColVec<NV, T> y{};
        for (size_t i = 0; i < NV; ++i) {
            T sum{0};
            for (size_t j = 0; j < n; ++j) {
                sum += Q(i, j) * w(j);
            }
            y(i) = sum;
        }
        return y;
    }

    /**
     * @brief The point of the active face nearest y₀: y = Q₁·R⁻ᵀβ + Q₂·Q₂ᵀy₀
     *
     * The minimizer of the QP with the active constraints ñₐᵀy = βₐ held as equalities, from the
     * factors alone: the part along the active normals that meets them, and y₀ in the free
     * directions.
     *
     * @see "Numerical Optimization" (Nocedal & Wright, 2006), §16.2: the null-space method
     */
    [[nodiscard]] constexpr ColVec<NV, T> nearest(const ColVec<NV, T>& y0, const ColVec<KMAX, T>& beta) const {
        ColVec<NV, T> y = range_point(beta);
        for (size_t j = k; j < NV; ++j) {
            T dot{0};
            for (size_t i = 0; i < NV; ++i) {
                dot += Q(i, j) * y0(i);
            }
            for (size_t i = 0; i < NV; ++i) {
                y(i) += Q(i, j) * dot;
            }
        }
        return y;
    }
};

/**
 * @brief Sine above which a normal is independent of the active ones in any metric: √(64·(KMAX+1)·ε)
 *
 * The square of a sine at this level is what a Gram factor resolves. Between dependence_tol()
 * and this, the metric of H⁻¹ alone does not tell a row that a heavy weight hides from the
 * roundoff of a dependent row, which grows with the conditioning of L; the row itself does
 * (row_stands_clear).
 *
 * Its limit: row_stands_clear takes this as its tolerance too, so a row that reaches a heavily
 * weighted variable only through a coefficient below it, against the row's largest, reads as
 * dependent. In float that is about 5·10⁻³ (more with more active rows), and a feasible problem
 * that can only meet such a row through that variable can end Infeasible. In double the level
 * is about 10⁻⁷.
 */
template<size_t KMAX, typename T>
[[nodiscard]] constexpr T clear_tol() {
    return damp::sqrt(static_cast<T>(KMAX + 1) * T{64} * eps<T>());
}

/**
 * @brief Whether row @p p stands clear of the first @p k active rows as a row: the largest
 *        entry of aₚ − Σ rₐaₐ against @p tol of the largest column's terms
 *
 * Linear dependence does not depend on the metric. With r the coefficients of nₚ on the active
 * normals in the metric of H⁻¹, a row dependent on them leaves only roundoff here, formed from
 * the rows as given; a row that reaches a variable the active rows do not leaves that entry
 * whole, however heavily H weights the variable.
 */
template<size_t NV, size_t NI, size_t KMAX, typename T>
[[nodiscard]] constexpr bool row_stands_clear(
    const Matrix<NI, NV, T>&         A,
    const damp::array<size_t, KMAX>& active,
    const ColVec<KMAX, T>&           r,
    size_t                           k,
    size_t                           p,
    T                                tol
) {
    k = (k < KMAX) ? k : KMAX; // provable index bound for the optimizer
    T res{0};
    T res_mag{0};
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
    return res > tol * res_mag;
}

/**
 * @brief Farkas certificate that candidate row @p p cannot be met with the active rows
 *
 * With r the coefficients of nₚ on the first @p k active normals, all at or below the ratio
 * test's floor, y = (−r, 1) ≥ 0 on (active, p). The rows are infeasible when Σ yᵢaᵢ ≈ 0 (the
 * row does not stand clear of the active ones by clear_tol(), row_stands_clear) and
 * Σ yᵢbᵢ < 0 past the roundoff of its terms. Without it the normal is only close to dependent in the metric of
 * H⁻¹, and the problem may have a solution.
 */
template<size_t NV, size_t NI, size_t KMAX, typename T>
[[nodiscard]] constexpr bool farkas_infeasible(
    const Matrix<NI, NV, T>&         A,
    const ColVec<NI, T>&             b,
    const damp::array<size_t, KMAX>& active,
    const ColVec<KMAX, T>&           r,
    size_t                           k,
    size_t                           p
) {
    constexpr T gap_tol = static_cast<T>(NV + 1) * T{64} * eps<T>();
    k = (k < KMAX) ? k : KMAX; // provable index bound for the optimizer
    T gap = b(p);
    T gap_mag = damp::abs(b(p));
    for (size_t a = 0; a < k; ++a) {
        gap -= r(a) * b(active[a]);
        gap_mag += damp::abs(r(a) * b(active[a]));
    }
    return !row_stands_clear(A, active, r, k, p, clear_tol<KMAX, T>()) && (gap < -(gap_tol * gap_mag));
}

/**
 * @brief Goldfarb–Idnani dual active-set core with a pre-factored Hessian
 *
 * Shared by solve_qp() (cold: factors H per call, empty hint) and
 * WarmStartActiveSetSolver (caches the factor and seeds the hint from the
 * previous solve).
 *
 * The method. From the unconstrained minimizer, each iteration picks a violated constraint p
 * and steps toward it along the part of H⁻¹nₚ that keeps the active constraints tight, in the
 * primal and in the multipliers together. The step stops where p becomes tight (a full step: p
 * joins the active set) or, sooner, where a multiplier reaches zero (a partial step: that
 * constraint leaves, and the step goes on). The multipliers stay nonnegative throughout; the
 * point is primal feasible only at Success.
 *
 * The linear algebra is in y = Lᵀx (H = LLᵀ), where the metric of H⁻¹ is the Euclidean one, on
 * an orthogonal factorization of the active normals (ActiveNormals):
 *
 * - The dependence test is on the sine of the angle between ñₚ = L⁻¹nₚ and the span of the
 *   active normals, ‖c₂‖/‖ñₚ‖. Above clear_tol() the normal is independent, at or below
 *   dependence_tol() dependent. In between the row decides (row_stands_clear): aₚ − Σ rₐaₐ,
 *   formed from the rows as given, is roundoff for a dependent row whatever the conditioning
 *   of L, and whole for a row that reaches a variable the active rows do not. A dependent
 *   normal has no primal direction, and its step is in the multipliers alone.
 * - The step length to p is −sₚ/‖c₂‖²: a sum of squares, never a difference of two products.
 * - After a full step the point is not the sum of the steps taken. It is recomputed as the
 *   minimizer on the active face (ActiveNormals::nearest), then corrected once on the active
 *   rows evaluated in x (one step of iterative refinement, which removes the cond(L)·ε the
 *   back substitution to x leaves). Each active constraint is then met to the roundoff of its
 *   own terms, and the point carries no error from one iteration to the next. The multipliers
 *   are accumulated along the steps, as in the paper: each step adds a multiple of
 *   ñₚ − Ñ·r to their combination, so stationarity holds to ε of its largest term.
 *
 * What it guards against. A constraint that can only be met by moving a heavily weighted
 * variable has a normal close to the span of the active ones in the metric of H⁻¹, though far
 * from it as a row. A soft-constrained MPC is the case: its slack carries the weight ρ = 10⁵
 * against move weights of 10⁻², and once the moves are pinned by their limits a further row is
 * met through the slack alone. On a double integrator (NP = 30, NC = 10, 11 variables) the sine
 * of such a normal is 2·10⁻⁴ to 10⁻², its square 5·10⁻⁸ to 10⁻⁴. A Cholesky factor of the Gram
 * matrix NᵀH⁻¹N holds the square, computed by cancellation, against a pivot floor of 2·10⁻⁴ in
 * single precision: the float solve ended Degenerate on 17 to 35 of that loop's 300 ticks
 * (the count depends on the optimization flags), and from a start 13 times farther out on 295
 * of them, or reported Infeasible on all 300. The orthogonal factors hold the sine itself,
 * against a floor of 10⁻⁵.
 *
 * Limits. A normal within dependence_tol() of the active span is still taken as dependent:
 * in single precision at 11 variables a sine of 1.1·10⁻⁵. The solve then ends Degenerate when no
 * active constraint can leave and there is no Farkas certificate. In the case above the sine
 * falls as 1/√ρ: a slack weight of 10⁷ is on the floor (it still solves on every tick), and at
 * 10⁸ some ticks end Degenerate. The minimizer on a face whose normals meet
 * at the sine s is known to about ε/s of the constraint data, and its multipliers to ε/s of the
 * largest, so a nearly dependent active set costs digits in proportion. ñ = L⁻¹n carries the
 * forward error of the triangular solve, which grows with the conditioning of L: a dependent
 * normal then reads above dependence_tol(), and the row test is what refuses it. Against an
 * enumeration of the active sets on QPs of 4 variables with exactly dependent rows, every
 * answer is right in both precisions with L graded by up to 10; graded by 100, a double solve
 * in 3000 returns Success on an infeasible problem, and single precision is past cond(H)·ε.
 * Success means every row met, active or not, by one test (the active ones refined, then
 * tested like the others); the multipliers are nonnegative by the ratio test and are not
 * tested again.
 *
 * At the budget: a solve whose last allowed iteration completes it is Success. A MaxIterations
 * point is the minimizer on the last complete active set, moved along the partial steps taken
 * since.
 *
 * The hint restricts each constraint pick to the most-violated *hinted* row
 * while any hinted row is violated, falling back to the full most-violated
 * scan — with a good hint (the previous solve's active set) the pick costs
 * O(hint·NV) instead of O(NI·NV). Any violated constraint is a valid pick, so
 * a stale or wrong hint affects speed, never correctness.
 *
 * @see "A numerically stable dual method for solving strictly convex quadratic programs"
 *      (Goldfarb & Idnani, 1983), Mathematical Programming 27: the dual algorithm, and its
 *      implementation on the orthogonal factors of the active set
 * @see "Numerical Optimization" (Nocedal & Wright, 2006), §16.2: the null-space method, the
 *      minimizer on a set of active constraints
 * @see "Accuracy and Stability of Numerical Algorithms" (Higham, 2002), ch. 12: iterative
 *      refinement in fixed precision
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
    constexpr T dep_tol = dependence_tol<KMAX, T>();
    constexpr T clear = clear_tol<KMAX, T>();
    constexpr T row_tol = static_cast<T>(KMAX + 1) * T{64} * eps<T>();
    constexpr T slack_tol = static_cast<T>(NV + 1) * T{64} * eps<T>();
    constexpr T big = unbounded_bound<T>();
    constexpr T skip_threshold = big / T{4};

    QPResult<NV, NI, T> result{};

    // y = Lᵀx. The unconstrained minimizer is y₀ = −L⁻¹f, x₀ = L⁻ᵀy₀.
    const ColVec<NV, T> y0 = -mat::forward_substitute(L, f);
    ColVec<NV, T>       y = y0;
    ColVec<NV, T>       x = mat::backward_substitute_transpose(L, y0);

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

    // Active set: the factors of its normals, and per constraint its row, multiplier, ‖ñₐ‖ and βₐ.
    ActiveNormals<NV, KMAX, T> normals{};
    damp::array<size_t, KMAX>  active{};
    damp::array<T, KMAX>       lam{};
    damp::array<T, KMAX>       length{};
    ColVec<KMAX, T>            beta{};
    damp::array<bool, NI>      in_active{};

    const auto finish = [&](QPStatus status) {
        for (size_t a = 0; a < normals.k && a < KMAX; ++a) {
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
        if (p == NI) {
            // No inactive row is violated. Success is every row met: the active ones, each
            // refined after the last step, take the same test. One that fails it was stepped to
            // on a free part at roundoff, and the active set cannot be held at this precision.
            for (size_t a = 0; a < normals.k && a < KMAX; ++a) {
                if (violated(active[a], s_i)) {
                    return finish(non_finite ? QPStatus::NonFinite : QPStatus::Degenerate);
                }
            }
        }
        if (non_finite) {
            return finish(QPStatus::NonFinite);
        }
        if (p == NI) {
            return finish(QPStatus::Success);
        }
        if (result.iterations >= max_iterations) {
            return finish(QPStatus::MaxIterations); // x: the minimizer on the last active set
        }

        const ColVec<NV, T> v_p = mat::forward_substitute(L, normal(p)); // ñₚ = L⁻¹nₚ
        T                   npup = T{0};                                 // ‖ñₚ‖² = nₚᵀH⁻¹nₚ
        for (size_t v = 0; v < NV; ++v) {
            npup += v_p(v) * v_p(v);
        }
        const T       v_len = damp::sqrt(npup);
        ColVec<NV, T> c = normals.rotate(v_p);
        T             s_p = worst; // the slack of p, < 0: followed along the partial steps
        T             lam_p = T{0};

        // Inner loop: step toward constraint p, dropping blockers as needed.
        while (true) {
            const size_t          k = (normals.k < KMAX) ? normals.k : KMAX;
            const ColVec<KMAX, T> r = normals.coefficients(c);
            // κ = ‖c₂‖² = nₚᵀ(I − projection)H⁻¹nₚ: at or below dep_tol²·‖ñₚ‖², nₚ is dependent
            // on the active normals and only a dual (drop) step is possible. With KMAX normals
            // active they span ℝ^NV and nothing is free.
            const T kappa = normals.free_norm_sq(c);
            bool    dependent = !(k < KMAX && kappa > dep_tol * dep_tol * npup);
            if (!dependent && !(kappa > clear * clear * npup)) {
                // A small free part is a row a heavy weight of H hides, or the roundoff of a
                // dependent row, which grows with the conditioning of L. The row tells which.
                dependent = !row_stands_clear(A, active, r, k, p, row_tol);
            }
            const T t2 = dependent ? big : (-s_p / kappa);

            T      t1 = big;
            size_t drop = KMAX;
            for (size_t a = 0; a < k; ++a) {
                // r_a carries units ‖ñₚ‖/‖ñ_a‖; compare in those units.
                const T r_floor = dep_tol * v_len / length[a];
                if (r(a) > r_floor && lam[a] / r(a) < t1) {
                    t1 = lam[a] / r(a);
                    drop = a;
                }
            }

            if (t1 >= big && t2 >= big) {
                // Dependent and nothing can leave: Infeasible with a Farkas certificate, else
                // nₚ is only nearly dependent in the metric of H⁻¹ and the active set cannot
                // grow at this precision: Degenerate.
                return finish(farkas_infeasible(A, b, active, r, k, p) ? QPStatus::Infeasible : QPStatus::Degenerate);
            }

            // The step in the multipliers, of the shorter length. The ratio test keeps every one
            // at or above zero; the clamp takes the roundoff of the one that reaches it.
            const T t = (t1 < t2) ? t1 : t2;
            for (size_t a = 0; a < k; ++a) {
                lam[a] = damp::max(lam[a] - (t * r(a)), T{0});
            }
            lam_p += t;
            ++result.iterations;

            if (t2 <= t1) {
                // Full step: p joins the active set, and x is the minimizer on the new face.
                if (k >= KMAX) {
                    return finish(QPStatus::Degenerate); // unreachable (t2 is finite only with k < KMAX): the index bound for the optimizer
                }
                normals.append(c);
                active[k] = p;
                lam[k] = lam_p;
                length[k] = v_len;
                beta(k) = -b(p);
                in_active[p] = true;
                y = normals.nearest(y0, beta);
                x = mat::backward_substitute_transpose(L, y);
                // One step of refinement on the active rows, in x: the back substitution left
                // each of them off by cond(L)·ε of x, and the correction, added to x, leaves
                // the roundoff of the row's own terms. In single precision the first is past the
                // tolerance the inactive rows are held to.
                ColVec<KMAX, T> miss{};
                for (size_t a = 0; a <= k; ++a) {
                    miss(a) = -slack(active[a]);
                }
                const ColVec<NV, T> dy = normals.range_point(miss);
                const ColVec<NV, T> dx = mat::backward_substitute_transpose(L, dy);
                for (size_t v = 0; v < NV; ++v) {
                    y(v) += dy(v);
                    x(v) += dx(v);
                }
                break;
            }

            // Partial step: the blocking constraint leaves. A dependent nₚ has no primal
            // direction (Goldfarb & Idnani's step in the dual alone): its computed free part is
            // roundoff, and the step length, a multiplier as large as the cost's weights, would
            // carry it into the point.
            if (!dependent) {
                const ColVec<NV, T> z = normals.free_part(c);
                for (size_t v = 0; v < NV; ++v) {
                    y(v) += t * z(v);
                }
                s_p += t * kappa;
            }
            in_active[active[drop]] = false;
            for (size_t a = drop; a + 1 < k; ++a) {
                active[a] = active[a + 1];
                lam[a] = lam[a + 1];
                length[a] = length[a + 1];
                beta(a) = beta(a + 1);
            }
            normals.remove(drop, c);
            if (result.iterations >= max_iterations) {
                x = mat::backward_substitute_transpose(L, y);
                return finish(QPStatus::MaxIterations);
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
 * Each active-set change costs @f$ O(N_V^2) @f$: triangular solves against the
 * Cholesky factor of H, and the update of the orthogonal factors of the active
 * normals, which resolve a nearly dependent constraint in single precision
 * (detail::gi_active_set has the method, what it guards against and its limits).
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
    requires damp::is_floating_point_v<T>
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
    requires damp::is_floating_point_v<T>
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
    requires damp::is_floating_point_v<T>
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
    requires damp::is_floating_point_v<T>
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
    requires damp::is_floating_point_v<T>
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
    requires damp::is_floating_point_v<T>
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
