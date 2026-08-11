// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file solve.hpp
 * @brief Linear-system solvers built on the factorizations in decomposition.hpp:
 *        triangular substitution, Cholesky/LU solves, the SPD-aware generic
 *        solve(), and the Matrix::inverse() definition.
 */

#include <cstddef>
#include <type_traits>

#include "core.hpp"
#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "decomposition.hpp"
#include "matrix_traits.hpp"

namespace damp {

namespace mat {

/**
 * @brief Forward substitution for Lx = b
 *
 * Assumes a nonsingular lower-triangular @f$ L @f$ (unit diagonal from LU, or
 * positive diagonal from Cholesky). Does not check pivots — use the
 * optional-returning solve(const LowerTriangle&) overload when the factor
 * may be singular. Call only after a successful factorization.
 */
template<size_t N, typename T>
constexpr ColVec<N, T> forward_substitute(const Matrix<N, N, T>& L, const ColVec<N, T>& b) {
    ColVec<N, T> x;

    for (size_t i = 0; i < N; ++i) {
        T sum = b(i);
        for (size_t j = 0; j < i; ++j) {
            sum -= L(i, j) * x(j);
        }

        x(i) = sum / L(i, i); // unit-diagonal LU (1/1) or Cholesky diagonal
    }

    return x;
}

/**
 * @brief Backward substitution for Lᵀx = b (real) / Lᴴx = b layout
 *
 * Solves using the lower factor stored in @p L by walking the transpose pattern
 * (@c L(j,i) for @c j > i). Caller must ensure nonzero diagonal (post-Cholesky).
 */
template<size_t N, typename T>
constexpr ColVec<N, T> backward_substitute_transpose(const Matrix<N, N, T>& L, const ColVec<N, T>& b) {
    ColVec<N, T> x;

    for (size_t i = N; i-- > 0;) {
        T sum = b(i);
        for (size_t j = i + 1; j < N; ++j) {
            sum -= L(j, i) * x(j);
        }

        x(i) = sum / L(i, i);
    }

    return x;
}

/**
 * @brief Solve lower-triangular system LX = B via forward substitution
 *
 * Requires nonsingular lower-triangular @p L. Returns damp::nullopt if any
 * diagonal magnitude is below default_tol.
 */
template<size_t N, size_t M, typename T>
constexpr damp::optional<Matrix<N, M, std::remove_const_t<T>>>
solve(const LowerTriangle<N, T>& L, const Matrix<N, M, std::remove_const_t<T>>& B) {
    using VT = std::remove_const_t<T>;
    constexpr auto tol = default_tol<VT>();

    Matrix<N, M, VT> X;
    for (size_t col = 0; col < M; ++col) {
        ColVec<N, VT> b;
        for (size_t i = 0; i < N; ++i) {
            b(i) = B(i, col);
        }

        ColVec<N, VT> x;
        for (size_t i = 0; i < N; ++i) {
            if (damp::abs(L(i, i)) < tol) {
                return damp::nullopt;
            }
            VT sum = b(i);
            for (size_t j = 0; j < i; ++j) {
                sum -= L(i, j) * x(j);
            }
            x(i) = sum / L(i, i);
        }

        for (size_t i = 0; i < N; ++i) {
            X(i, col) = x(i);
        }
    }
    return X;
}

/**
 * @brief Solve upper-triangular system UX = B via backward substitution
 *
 * Returns damp::nullopt if any diagonal magnitude is below default_tol.
 */
template<size_t N, size_t M, typename T>
constexpr damp::optional<Matrix<N, M, std::remove_const_t<T>>>
solve(const UpperTriangle<N, T>& U, const Matrix<N, M, std::remove_const_t<T>>& B) {
    using VT = std::remove_const_t<T>;
    constexpr auto tol = default_tol<VT>();

    Matrix<N, M, VT> X;
    for (size_t col = 0; col < M; ++col) {
        ColVec<N, VT> b;
        for (size_t i = 0; i < N; ++i) {
            b(i) = B(i, col);
        }

        ColVec<N, VT> x;
        for (size_t i = N; i-- > 0;) {
            if (damp::abs(U(i, i)) < tol) {
                return damp::nullopt;
            }
            VT sum = b(i);
            for (size_t j = i + 1; j < N; ++j) {
                sum -= U(i, j) * x(j);
            }
            x(i) = sum / U(i, i);
        }

        for (size_t i = 0; i < N; ++i) {
            X(i, col) = x(i);
        }
    }
    return X;
}

/**
 * @brief Solve AX = B via Cholesky (A = LLᴴ)
 *
 * @p A must be Hermitian positive definite. Returns damp::nullopt if Cholesky
 * fails or a triangular solve hits a zero diagonal.
 *
 * @note Prefer this over forming @f$ A^{-1} @f$ then multiplying.
 * @see cholesky(), solve()
 */
template<size_t N, size_t M, typename T>
constexpr damp::optional<Matrix<N, M, T>> cholesky_solve(const Matrix<N, N, T>& A, const Matrix<N, M, T>& B) {
    auto Lopt = cholesky(A);
    if (!Lopt) {
        return damp::nullopt;
    }

    const auto& L = Lopt.value();

    // Forward substitution: solve L * Y = B
    auto Y_opt = solve(L.lower_triangle(), B);
    if (!Y_opt) {
        return damp::nullopt;
    }

    // Backward substitution: solve Lᴴ * X = Y
    const auto Lh = L.conjugate_transpose();
    return solve(Lh.upper_triangle(), Y_opt.value());
}

/**
 * @brief Solve AX = B via LU with partial pivoting
 *
 * Returns damp::nullopt if the factorization or a triangular solve fails.
 *
 * @note Prefer this over forming @f$ A^{-1} @f$ then multiplying.
 * @see lu_decomposition(), solve()
 */
template<size_t N, size_t M, typename T>
constexpr damp::optional<Matrix<N, M, T>> lu_solve(const Matrix<N, N, T>& A, const Matrix<N, M, T>& B) {
    auto lu_opt = lu_decomposition(A);
    if (!lu_opt) {
        return damp::nullopt;
    }

    const auto& [L, U, piv] = lu_opt.value();

    // Apply row permutation to B
    Matrix<N, M, T> B_perm;
    for (size_t r = 0; r < N; ++r) {
        for (size_t c = 0; c < M; ++c) {
            B_perm(r, c) = B(piv[r], c);
        }
    }

    // Forward substitution: solve L * Y = P*B
    auto Y_opt = solve(L.lower_triangle(), B_perm);
    if (!Y_opt) {
        return damp::nullopt;
    }

    // Backward substitution: solve U * X = Y
    return solve(U.upper_triangle(), Y_opt.value());
}

/**
 * @brief Solve AX = B — Cholesky when Hermitian, else LU
 *
 * Primary entry for linear solves. Prefer this over Matrix::inverse when
 * the goal is @f$ A^{-1} B @f$ rather than the inverse matrix itself.
 *
 * @return Solution X, or damp::nullopt if @p A is singular / not PD under Cholesky
 * @note Compare with MATLAB®'s A \\ B.
 * @see cholesky_solve(), lu_solve()
 */
template<size_t N, size_t M, typename T>
constexpr damp::optional<Matrix<N, M, T>> solve(const Matrix<N, N, T>& A, const Matrix<N, M, T>& B) {
    if (is_symmetric_or_hermitian<N, T>(A)) {
        auto Xopt = cholesky_solve(A, B);
        if (Xopt) {
            return Xopt;
        }
    }

    return lu_solve(A, B);
}

} // namespace mat

/**
 * @brief Matrix inverse A⁻¹ via mat::solve(A, I)
 *
 * Use only when the inverse matrix itself is the deliverable. For
 * @f$ A^{-1} b @f$ / @f$ A^{-1} B @f$, call mat::solve directly.
 *
 * @return Inverse if nonsingular, damp::nullopt otherwise
 * @note Compare with MATLAB®'s inv(A).
 * @see mat::solve(), mat::cholesky_solve(), mat::lu_solve()
 */
template<size_t Rows, size_t Cols, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<Rows, Cols, T>> Matrix<Rows, Cols, T>::inverse() const
    requires(Rows == Cols)
{
    return mat::solve(*this, Matrix<Rows, Cols, T>::identity());
}

} // namespace damp
