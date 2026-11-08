// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file invariant_set.hpp
 * @brief Maximal positively invariant polytopes of x⁺ = A x under G₀ x ≤ h₀
 *
 * Gilbert & Tan iteration: start from the finite rows of G₀ x ≤ h₀ and add
 * G₀ Aᵏ rows until a full step is redundant. Storage is a fixed-capacity
 * matrix; overflow or a non-converged iteration is `not_converged`, never an
 * unfinished success.
 *
 * @see "Linear Systems with State and Control Constraints: The Theory and
 *      Application of Maximal Output Admissible Sets" (Gilbert & Tan, 1991)
 */


#include "damp/backend.hpp"
#include "damp/design/qp.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"

namespace damp {
namespace design {

/**
 * @brief Maximal positively invariant polytope G e ≤ h of x⁺ = A x
 *
 * `n_rows` of `G` / `h` are the live facets; unused capacity is the unbounded
 * sentinel so a fixed-size QP can skip it.
 */
template<size_t NT_max, size_t NX, size_t NG, typename T = double, size_t K_max = 32>
struct [[nodiscard]] InvariantSetResult {
    Matrix<NT_max, NX, T> G{};                                 ///< Facet normals; first n_rows used
    ColVec<NT_max, T>     h{};                                 ///< Facet offsets
    size_t                n_rows{0};                           ///< Live row count
    bool                  success{false};                      ///< status == DesignStatus::ok
    DesignStatus          status{DesignStatus::not_converged}; ///< ok, or why the set failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return InvariantSetResult<NT_max, NX, NG, U, K_max>{
            G.template as<U>(),
            detail::convert_bounds<U>(h),
            n_rows,
            success,
            status,
        };
    }
};

namespace detail {

template<typename T>
[[nodiscard]] constexpr bool bound_is_unbounded(T h) {
    return h >= unbounded_bound<T>() / T{4};
}

template<size_t NX, typename T>
[[nodiscard]] constexpr T ridge_weight() {
    return static_cast<T>(NX + 1) * damp::sqrt(eps<T>());
}

template<size_t N, typename T>
[[nodiscard]] constexpr Matrix<N, N, T> ridge_hessian() {
    Matrix<N, N, T> H{};
    const T         w = ridge_weight<N, T>();
    for (size_t i = 0; i < N; ++i) {
        H(i, i) = w;
    }
    return H;
}

template<size_t NX, typename T>
[[nodiscard]] constexpr T row_abs_sum(const ColVec<NX, T>& a) {
    T s{0};
    for (size_t j = 0; j < NX; ++j) {
        s += damp::abs(a(j));
    }
    return s;
}

/** @brief A row of G₀ that is zero to roundoff relative to the largest row, @p g_scale (its ‖·‖₁). */
template<size_t NX, typename T>
[[nodiscard]] constexpr bool row_is_zero(const ColVec<NX, T>& a, T g_scale) {
    const T mag = row_abs_sum(a);
    return mag <= static_cast<T>(NX) * eps<T>() * g_scale;
}

template<size_t NR, size_t NX, typename T>
[[nodiscard]] constexpr ColVec<NX, T> matrix_row(const Matrix<NR, NX, T>& G, size_t i) {
    ColVec<NX, T> a{};
    for (size_t j = 0; j < NX; ++j) {
        a(j) = G(i, j);
    }
    return a;
}

template<size_t N, typename T>
constexpr void set_unbounded_offsets(ColVec<N, T>& h) {
    for (size_t i = 0; i < N; ++i) {
        h(i) = unbounded_bound<T>();
    }
}

/** @brief Candidate aᵀe ≤ β is redundant iff max aᵀe on G e ≤ h is ≤ β (data-scaled tol). */
enum class RowRedundancy : uint8_t { redundant,
                                     needed,
                                     empty };

/** @brief The size of G e ≤ h: its largest facet distance |hᵢ|/‖Gᵢ‖ (1 when every hᵢ is 0). */
template<size_t NT_max, size_t NX, typename T>
[[nodiscard]] constexpr T polytope_scale(const Matrix<NT_max, NX, T>& G, const ColVec<NT_max, T>& h) {
    T rho{0};
    for (size_t i = 0; i < NT_max; ++i) {
        if (bound_is_unbounded(h(i))) {
            continue;
        }
        const T nrm = matrix_row(G, i).norm();
        if (nrm > T{0}) {
            rho = damp::max(rho, damp::abs(h(i)) / nrm);
        }
    }
    return (rho > T{0}) ? rho : T{1};
}

/**
 * @brief Whether aᵀe ≤ β is redundant on G e ≤ h, from the ridge QP min ½w‖ẽ‖² − âᵀẽ in units
 *        of the polytope (ẽ = e/ρ, â = a/‖a‖).
 *
 * The ridge stands in for the LP max âᵀẽ. In a direction where the polytope is unbounded it
 * returns âᵀẽ ≈ 1/w instead of +∞, so the test is posed in units where the polytope has size 1
 * and β/(ρ‖a‖) is compared with it: scaling h by any factor scales ρ and β alike, and the
 * outcome does not change.
 */
template<size_t NT_max, size_t NX, typename T>
[[nodiscard]] constexpr RowRedundancy row_redundancy(
    const Matrix<NT_max, NX, T>& G,
    const ColVec<NT_max, T>&     h,
    const ColVec<NX, T>&         a,
    T                            beta
) {
    if constexpr (NX == 0) {
        return (beta >= T{0}) ? RowRedundancy::redundant : RowRedundancy::empty;
    }
    const T a_norm = a.norm();
    if (!(a_norm > T{0})) {
        return (beta >= T{0}) ? RowRedundancy::redundant : RowRedundancy::empty; // 0 ≤ β (A singular)
    }
    const T           rho = polytope_scale(G, h);
    ColVec<NT_max, T> h_unit{};
    for (size_t i = 0; i < NT_max; ++i) {
        h_unit(i) = bound_is_unbounded(h(i)) ? h(i) : (h(i) / rho);
    }
    const ColVec<NX, T> a_unit{a * (T{1} / a_norm)};
    const T             beta_unit = beta / (rho * a_norm);
    const auto          qp = solve_qp(ridge_hessian<NX, T>(), ColVec<NX, T>(-a_unit), G, h_unit);
    if (qp.status == QPStatus::Infeasible) {
        return RowRedundancy::empty;
    }
    T max_dot{0};
    for (size_t j = 0; j < NX; ++j) {
        max_dot += a_unit(j) * qp.x(j);
    }
    const T mag = row_abs_sum(a_unit) + damp::abs(beta_unit) + damp::abs(max_dot);
    const T tol = static_cast<T>(8 * (NX + 1)) * eps<T>() * mag;
    if (qp.status == QPStatus::Success && max_dot <= beta_unit + tol) {
        return RowRedundancy::redundant;
    }
    return RowRedundancy::needed;
}

/** @brief Chebyshev radius of G e ≤ h: max r s.t. Gᵢ e + r ‖Gᵢ‖₂ ≤ hᵢ. r ≤ 0 is empty. */
template<size_t NT_max, size_t NX, typename T>
[[nodiscard]] constexpr T chebyshev_radius(const Matrix<NT_max, NX, T>& G, const ColVec<NT_max, T>& h, size_t n_rows) {
    if (n_rows == 0) {
        return T{1};
    }
    constexpr size_t      NV = NX + 1;
    Matrix<NT_max, NV, T> A{};
    ColVec<NV, T>         f{};
    f(NX) = T{-1};
    for (size_t i = 0; i < n_rows; ++i) {
        const auto a = matrix_row(G, i);
        const T    nrm = a.norm();
        for (size_t j = 0; j < NX; ++j) {
            A(i, j) = a(j);
        }
        A(i, NX) = nrm;
    }
    const auto qp = solve_qp(ridge_hessian<NV, T>(), f, A, h);
    if (qp.status != QPStatus::Success) {
        return T{-1};
    }
    return qp.x(NX);
}

template<size_t NT_max, size_t NX, typename T>
[[nodiscard]] constexpr bool append_facet(
    Matrix<NT_max, NX, T>& G,
    ColVec<NT_max, T>&     h,
    size_t&                n_rows,
    const ColVec<NX, T>&   a,
    T                      beta
) {
    if (n_rows >= NT_max) {
        return false;
    }
    for (size_t j = 0; j < NX; ++j) {
        G(n_rows, j) = a(j);
    }
    h(n_rows) = beta;
    ++n_rows;
    return true;
}

} // namespace detail

/**
 * @brief Maximal positively invariant set of x⁺ = A x under G₀ x ≤ h₀
 *
 * Skips unbounded G₀ rows (hᵢ ≥ unbounded_bound()/4). A zero row with a
 * negative offset is empty. Iteration cap is @p K_max (default 32); storage
 * cap is @p NT_max. The returned polytope is certified non-empty by a
 * Chebyshev-center LP (r > 0).
 *
 * @note Compare with MATLAB®'s mpc TerminalSet / setTerminalSet.
 * @see "Linear Systems with State and Control Constraints: The Theory and
 *      Application of Maximal Output Admissible Sets" (Gilbert & Tan, 1991)
 *
 * @tparam NT_max Facet capacity
 * @tparam K_max  Maximum propagation steps (default 32)
 * @param A   Discrete dynamics (should be Schur; not checked here)
 * @param G0  Initial facet normals G₀
 * @param h0  Initial facet offsets h₀
 * @return G, h, n_rows; success only when the iteration closed and the set is non-empty
 */
template<size_t NT_max, size_t K_max = 32, size_t NX, size_t NG, typename T>
[[nodiscard]] constexpr InvariantSetResult<NT_max, NX, NG, T, K_max>
maximal_invariant_set(const Matrix<NX, NX, T>& A, const Matrix<NG, NX, T>& G0, const ColVec<NG, T>& h0) {
    using Result = InvariantSetResult<NT_max, NX, NG, T, K_max>;
    Result out{};
    detail::set_unbounded_offsets(out.h);

    T g_scale{0}; // the largest live row of G₀ (‖·‖₁): a zero row is zero relative to it
    for (size_t i = 0; i < NG; ++i) {
        if (!detail::bound_is_unbounded(h0(i))) {
            g_scale = damp::max(g_scale, detail::row_abs_sum(detail::matrix_row(G0, i)));
        }
    }
    size_t n_usable = 0;
    for (size_t i = 0; i < NG; ++i) {
        if (detail::bound_is_unbounded(h0(i))) {
            continue;
        }
        const auto a = detail::matrix_row(G0, i);
        if (detail::row_is_zero(a, g_scale)) {
            if (h0(i) < T{0}) {
                return damp::detail::set_status(out, DesignStatus::not_converged);
            }
            continue;
        }
        ++n_usable;
    }
    if (n_usable > NT_max) {
        return damp::detail::set_status(out, DesignStatus::not_converged);
    }

    for (size_t i = 0; i < NG; ++i) {
        if (detail::bound_is_unbounded(h0(i))) {
            continue;
        }
        const auto a = detail::matrix_row(G0, i);
        if (detail::row_is_zero(a, g_scale)) {
            continue;
        }
        (void)detail::append_facet(out.G, out.h, out.n_rows, a, h0(i));
    }

    if (!(detail::chebyshev_radius(out.G, out.h, out.n_rows) > T{0})) {
        return damp::detail::set_status(out, DesignStatus::not_converged);
    }
    if (out.n_rows == 0) {
        return damp::detail::set_status(out, DesignStatus::ok);
    }

    Matrix<NX, NX, T> A_pow = Matrix<NX, NX, T>::identity();
    for (size_t k = 0; k < K_max; ++k) {
        A_pow = A_pow * A;
        Matrix<NG, NX, T> pending_G{};
        ColVec<NG, T>     pending_h{};
        size_t            n_pending = 0;
        for (size_t i = 0; i < NG; ++i) {
            if (detail::bound_is_unbounded(h0(i))) {
                continue;
            }
            const auto a0 = detail::matrix_row(G0, i);
            if (detail::row_is_zero(a0, g_scale)) {
                continue;
            }
            const ColVec<NX, T> a{A_pow.transpose() * a0};
            if (a == a0) {
                continue; // an A-invariant row (e.g. a constant target's): already in the set since k = 0
            }
            const auto red = detail::row_redundancy(out.G, out.h, a, h0(i));
            if (red == detail::RowRedundancy::empty) {
                return damp::detail::set_status(out, DesignStatus::not_converged);
            }
            if (red == detail::RowRedundancy::redundant) {
                continue;
            }
            for (size_t j = 0; j < NX; ++j) {
                pending_G(n_pending, j) = a(j);
            }
            pending_h(n_pending) = h0(i);
            ++n_pending;
        }
        if (n_pending == 0) {
            if (!(detail::chebyshev_radius(out.G, out.h, out.n_rows) > T{0})) {
                return damp::detail::set_status(out, DesignStatus::not_converged);
            }
            return damp::detail::set_status(out, DesignStatus::ok);
        }
        for (size_t p = 0; p < n_pending; ++p) {
            if (!detail::append_facet(out.G, out.h, out.n_rows, detail::matrix_row(pending_G, p), pending_h(p))) {
                return damp::detail::set_status(out, DesignStatus::not_converged);
            }
        }
    }
    return damp::detail::set_status(out, DesignStatus::not_converged);
}

} // namespace design
} // namespace damp
