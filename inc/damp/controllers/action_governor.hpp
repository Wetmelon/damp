// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file action_governor.hpp
 * @brief Command projection and action governors
 *
 * A minimally-invasive layer between any controller and the actuator: given a
 * desired command @f$u_{\mathrm{des}}@f$, return the nearest admissible
 * @f$u@f$ that satisfies the stated constraints. The filter is the identity
 * when @f$u_{\mathrm{des}}@f$ already satisfies the constraints.
 *
 * Two projection modes:
 * - Box @f$\text{lower} \le u \le \text{upper}@f$ (a Bounds) — closed-form per-channel clamp
 *   (no solver; bounded per-tick cost).
 * - Affine @f$A u \le b@f$ — Euclidean projection
 *   @f$\min_u \tfrac12\|u - u_{\mathrm{des}}\|^2@f$ subject to the inequalities,
 *   solved by design::solve_qp (Goldfarb–Idnani). Intended for small
 *   @f$N_U@f$ (1–4) on-target.
 *
 * A relative-degree-1 control barrier function (CBF) condition
 * @f$\dot h \ge -\alpha h@f$ is encoded as one affine row via
 * design::cbf_relative_degree_1.
 *
 * @code
 * using namespace damp;
 *
 * // Box: clamp a 2-input command into [-1, 1] × [-2, 2]
 * constexpr Bounds<2, float> box{.lower = {-1.f, -2.f}, .upper = {1.f, 2.f},};
 * auto u = project_box(ColVec<2, float>{1.5f, -3.f}, box); // → {1, -2}
 *
 * // Affine half-space: u₀ + u₁ ≤ 1
 * constexpr Matrix<1, 2> A{{1.0, 1.0}};
 * constexpr ColVec<1>    b{1.0};
 * constexpr auto res = design::project_affine(ColVec<2>{2.0, 2.0}, A, b);
 * static_assert(res.success);
 * // res.u ≈ {0.5, 0.5}
 * @endcode
 *
 * @see design::solve_qp
 * @see Ames et al., "Control Barrier Functions: Theory and Applications,"
 *      ECC 2019, https://doi.org/10.23919/ECC.2019.8796030
 * @see Gurriet et al., realizable command filters / action governors
 */


#include "damp/backend.hpp"
#include "damp/design/qp.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/toolbox/bounds.hpp"

namespace damp {

// =============================================================================
// Box projection (closed-form)
// =============================================================================

/**
 * @brief Euclidean projection of @p u onto a Bounds box
 *
 * Per-channel clamp: @f$u_i \leftarrow \mathrm{clamp}(u_i, \text{lower}_i, \text{upper}_i)@f$.
 * Identity when @p u already lies in the box.
 *
 * @param u      Desired command (NU)
 * @param bounds Box [lower, upper] (NU)
 * @return Projected command
 */
template<size_t NU, typename T = double>
    requires damp::is_floating_point_v<T>
[[nodiscard]] constexpr ColVec<NU, T> project_box(const ColVec<NU, T>& u, const Bounds<NU, T>& bounds) {
    return bounds.saturate(u);
}

// =============================================================================
// Affine projection (QP)
// =============================================================================

/**
 * @brief Result of an affine command projection
 *
 * @tparam NU Number of inputs
 * @tparam T  Scalar type
 */
template<size_t NU, typename T = double>
struct CommandProjectionResult {
    ColVec<NU, T>    u{};                                     ///< Projected command; on failure, @c u_des from project_affine and the last safe command from ActionGovernor::filter
    bool             success{false};                          ///< true if the QP found a feasible minimizer
    bool             modified{false};                         ///< true if @c u differs from the desired command
    design::QPStatus status{design::QPStatus::MaxIterations}; ///< Underlying QP status

    template<typename U>
    [[nodiscard]] constexpr CommandProjectionResult<NU, U> as() const {
        return CommandProjectionResult<NU, U>{
            u.template as<U>(),
            success,
            modified,
            status,
        };
    }
};

namespace design {

/**
 * @brief Project @p u_des onto the polyhedron A u ≤ b (Euclidean)
 *
 * Solves
 * @f[
 *   \min_u \tfrac12 \|u - u_{\mathrm{des}}\|^2
 *   \quad\text{s.t.}\quad A u \le b
 * @f]
 * via @ref solve_qp with @f$H = I@f$, @f$f = -u_{\mathrm{des}}@f$. Identity
 * when @p u_des already satisfies the inequalities. Rows of @p b at the
 * unbounded sentinel are ignored (see unbounded_bound).
 *
 * On failure (infeasible / non-convergence), @c u is left as @p u_des and
 * @c success is false. ActionGovernor's filter() and operator() are fail-closed:
 * they return the last successful projection (zero until the first success).
 *
 * @note Compare with a one-step action governor / command filter.
 * @see solve_qp, Ames et al. (ECC 2019), Gurriet et al.
 *
 * @param u_des  Desired command (NU)
 * @param A      Inequality matrix (NI × NU), rows aᵢᵀ u ≤ bᵢ
 * @param b      Inequality bounds (NI)
 * @param max_iterations  Active-set budget (default 10·(NU+NI))
 * @return CommandProjectionResult with projected @c u and status
 */
template<size_t NU, size_t NI, typename T = double>
    requires damp::is_floating_point_v<T>
[[nodiscard]] constexpr CommandProjectionResult<NU, T> project_affine(
    const ColVec<NU, T>&     u_des,
    const Matrix<NI, NU, T>& A,
    const ColVec<NI, T>&     b,
    size_t                   max_iterations = 10 * (NU + NI)
) {
    CommandProjectionResult<NU, T> out{};
    out.u = u_des;

    // min ½‖u − u_des‖² = ½ uᵀu − u_desᵀu + const  →  H = I, f = −u_des
    const Matrix<NU, NU, T> H = Matrix<NU, NU, T>::identity();
    ColVec<NU, T>           f{};
    for (size_t i = 0; i < NU; ++i) {
        f(i) = -u_des(i);
    }

    const auto qp = solve_qp(H, f, A, b, max_iterations);
    out.status = qp.status;
    out.success = qp.success;
    if (qp.success) {
        out.u = qp.x;
        constexpr T tol = damp::default_tol<T>();
        for (size_t i = 0; i < NU; ++i) {
            if (damp::abs(qp.x(i) - u_des(i)) > tol) {
                out.modified = true;
                break;
            }
        }
    }
    return out;
}

/**
 * @brief One affine row from a relative-degree-1 CBF condition
 *
 * For @f$\dot x = f(x) + g(x) u@f$ and barrier @f$h(x)@f$ (safe set
 * @f$\{h \ge 0\}@f$) of relative degree one, the CBF inequality
 * @f[
 *   \dot h(x,u) = L_f h(x) + L_g h(x)\, u \;\ge\; -\alpha\, h(x)
 * @f]
 * rearranges to the affine constraint
 * @f[
 *   \bigl(-L_g h\bigr)\, u \;\le\; L_f h + \alpha\, h.
 * @f]
 *
 * @see Ames et al., "Control Barrier Functions: Theory and Applications,"
 *      ECC 2019, https://doi.org/10.23919/ECC.2019.8796030
 *
 * @tparam NU Number of inputs
 * @tparam T  Scalar type
 */
template<size_t NU, typename T = double>
struct CBFConstraint {
    Matrix<1, NU, T> A{}; ///< Single row: (−L_g h) u ≤ b
    ColVec<1, T>     b{}; ///< L_f h + α h
};

/**
 * @brief Build a relative-degree-1 CBF inequality row
 *
 * @param Lg_h   Control Lie derivative L_g h (NU); (L_g h)ᵢ = ∇hᵀ g_·i
 * @param Lf_h   Drift Lie derivative L_f h = ∇hᵀ f
 * @param h      Barrier value h(x) (safe when ≥ 0)
 * @param alpha  Class-𝒦 gain α > 0 (default 1)
 * @return CBFConstraint with @c A and @c b ready for @ref project_affine
 */
template<size_t NU, typename T = double>
    requires damp::is_floating_point_v<T>
[[nodiscard]] constexpr CBFConstraint<NU, T> cbf_relative_degree_1(
    const ColVec<NU, T>& Lg_h,
    T                    Lf_h,
    T                    h,
    T                    alpha = T{1}
) {
    CBFConstraint<NU, T> row{};
    for (size_t i = 0; i < NU; ++i) {
        row.A(0, i) = -Lg_h(i);
    }
    row.b = Lf_h + (alpha * h);
    return row;
}

} // namespace design

/**
 * @brief Runtime action governor — QP projection of u_des onto A u ≤ b
 *
 * Stores a fixed (or tick-updated) inequality set and projects each call via
 * @ref design::project_affine. Default scalar type is @c float for embedded
 * deployment; for design-time checks prefer @c double via @ref design::project_affine.
 *
 * @tparam NU Number of inputs (keep small: 1–4)
 * @tparam NI Number of inequality rows
 * @tparam T  Scalar type (default float)
 *
 * @see design::project_affine, design::cbf_relative_degree_1
 * @see Ames et al., ECC 2019; Gurriet et al., action governors
 */
template<size_t NU, size_t NI, typename T = float>
    requires damp::is_floating_point_v<T>
class ActionGovernor {
public:
    constexpr ActionGovernor() = default;

    constexpr ActionGovernor(const Matrix<NI, NU, T>& A, const ColVec<NI, T>& b) : A_(A), b_(b) {}

    constexpr void set_constraints(const Matrix<NI, NU, T>& A, const ColVec<NI, T>& b) {
        A_ = A;
        b_ = b;
    }

    /** @brief Update only the bound vector (state-dependent b(x), fixed normals). */
    constexpr void set_bounds(const ColVec<NI, T>& b) { b_ = b; }

    [[nodiscard]] constexpr const Matrix<NI, NU, T>& A() const { return A_; }
    [[nodiscard]] constexpr const ColVec<NI, T>&     b() const { return b_; }

    /**
     * @brief Project @p u_des onto the stored polyhedron
     * @return CommandProjectionResult; on failure @c success is false and @c u is the last safe
     *         command (@ref last_safe, zero until the first success), never the unprojected @p u_des
     *
     * Updates the fail-closed hold (@ref last_safe) when @c success is true.
     */
    [[nodiscard]] constexpr CommandProjectionResult<NU, T> filter(const ColVec<NU, T>& u_des) {
        return hold(design::project_affine(u_des, A_, b_));
    }

    /**
     * @brief Project @p u_des with a one-shot bound vector @p b
     *
     * Does not overwrite the stored @c b_. Useful when @f$b = b(x)@f$ changes
     * every tick (e.g. CBF) while @f$A@f$ is fixed or also passed via
     * @ref set_constraints. Fail-closed like filter(u_des).
     */
    [[nodiscard]] constexpr CommandProjectionResult<NU, T> filter(
        const ColVec<NU, T>& u_des,
        const ColVec<NI, T>& b
    ) {
        return hold(design::project_affine(u_des, A_, b));
    }

    /**
     * @brief Project with full (A,b) for this tick only; fail-closed like filter(u_des)
     */
    [[nodiscard]] constexpr CommandProjectionResult<NU, T> filter(
        const ColVec<NU, T>&     u_des,
        const Matrix<NI, NU, T>& A,
        const ColVec<NI, T>&     b
    ) {
        return hold(design::project_affine(u_des, A, b));
    }

    /**
     * @brief Fail-closed projection: the command of @ref filter (last successful @c u, else zero if none yet)
     *
     * Prefer @ref filter when the caller must observe QP failure.
     */
    [[nodiscard]] constexpr ColVec<NU, T> operator()(const ColVec<NU, T>& u_des) { return filter(u_des).u; }

    /** @brief Last successful projection (zero until the first success). */
    [[nodiscard]] constexpr const ColVec<NU, T>& last_safe() const { return last_safe_; }

    /**
     * @brief Convert constraint storage; fail-closed hold restarts empty (design-time convert)
     *
     * Disabled rows (b at or past the unbounded sentinel) map to the target type's sentinel
     * rather than overflowing it.
     */
    template<typename U>
    [[nodiscard]] constexpr ActionGovernor<NU, NI, U> as() const {
        return ActionGovernor<NU, NI, U>{A_.template as<U>(), design::detail::convert_bounds<U>(b_)};
    }

private:
    /** @brief Fail closed: record a successful projection, or replace a failed one's command with the hold. */
    [[nodiscard]] constexpr CommandProjectionResult<NU, T> hold(CommandProjectionResult<NU, T> res) {
        if (res.success) {
            last_safe_ = res.u;
        } else {
            res.u = last_safe_;
            res.modified = true;
        }
        return res;
    }

    Matrix<NI, NU, T> A_{};
    ColVec<NI, T>     b_{};
    ColVec<NU, T>     last_safe_{}; ///< fail-closed hold
};

} // namespace damp
