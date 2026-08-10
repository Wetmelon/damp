// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file lyapunov.hpp
 * @brief Lyapunov equation solvers
 */

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/solve.hpp"

namespace damp {

namespace detail {

/**
 * @brief Solve a linear matrix equation L(X) + Q = 0 by Kronecker vectorization.
 *
 * The Lyapunov operators are linear in X, so @f$ L(X) = -Q @f$ is an ordinary
 * @f$ n^2 \times n^2 @f$ linear system in @c vec(X). The system matrix is built
 * one column at a time by applying @p apply to each basis matrix @f$ E_{ij} @f$
 * (row-major vectorization), then solved by dense LU.
 *
 * @param Q     Right-hand side (NX × NX).
 * @param apply Callable returning @f$ L(E) @f$ for a basis matrix @f$ E @f$.
 * @return X (NX × NX), or damp::nullopt if the system is singular (the operator
 *         has a zero eigenvalue — e.g. @f$ \lambda_i + \lambda_j = 0 @f$ for lyap,
 *         @f$ \lambda_i \lambda_j = 1 @f$ for dlyap — so no unique solution exists).
 *
 * @see Golub & Van Loan, "Matrix Computations" §12.3 (Kronecker / vec form).
 *
 * @note ponytail: O(NX^6) dense solve. Fine for the small fixed NX of an embedded
 *       plant. Swap to Bartels–Stewart (real Schur form, reusing the QR machinery
 *       already in design/riccati.hpp) only if NX grows large enough to matter.
 */
template<size_t NX, typename T, typename Op>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>>
solve_lyapunov_kron(const Matrix<NX, NX, T>& Q, Op apply) {
    constexpr size_t N2 = NX * NX;

    Matrix<N2, N2, T> K = Matrix<N2, N2, T>::zeros();
    Matrix<N2, 1, T>  rhs{};

    for (size_t i = 0; i < NX; ++i) {
        for (size_t j = 0; j < NX; ++j) {
            const size_t col = (i * NX) + j;

            Matrix<NX, NX, T> E = Matrix<NX, NX, T>::zeros();
            E(i, j) = T{1};
            const Matrix<NX, NX, T> L = apply(E);

            for (size_t r = 0; r < NX; ++r) {
                for (size_t c = 0; c < NX; ++c) {
                    K((r * NX) + c, col) = L(r, c);
                }
            }
            rhs(col, 0) = -Q(i, j);
        }
    }

    const auto x = mat::lu_solve(K, rhs);
    if (!x) {
        return damp::nullopt;
    }

    Matrix<NX, NX, T> X{};
    for (size_t r = 0; r < NX; ++r) {
        for (size_t c = 0; c < NX; ++c) {
            X(r, c) = x.value()((r * NX) + c, 0);
        }
    }
    return X;
}

} // namespace detail

/**
 * @brief Solve the continuous-time Lyapunov equation @f$ A X + X A^\top + Q = 0 @f$.
 *
 * The unique solution exists iff @f$ A @f$ and @f$ -A @f$ share no eigenvalue
 * (i.e. @f$ \lambda_i + \lambda_j \neq 0 @f$ for all @f$ i,j @f$). For a Hurwitz
 * @f$ A @f$ and symmetric positive-semidefinite @f$ Q @f$, @f$ X @f$ is symmetric
 * positive semidefinite and equals @f$ \int_0^\infty e^{At} Q\, e^{A^\top t}\,dt @f$
 * — the continuous controllability/observability Gramian when @f$ Q = BB^\top @f$
 * / @f$ C^\top C @f$ (see stability::controllability_gramian).
 *
 * @note Compare with MATLAB®'s @c X=lyap(A,Q) (which solves @f$ AX+XA^\top+Q=0 @f$).
 *
 * @see dlyap() — the discrete-time counterpart.
 * @see "Matrix Computations" (Golub & Van Loan), §12.3.
 *
 * @param A State matrix (NX × NX).
 * @param Q Symmetric right-hand side (NX × NX).
 * @return Solution X (NX × NX), or damp::nullopt if A and -A share an eigenvalue.
 */
template<size_t NX, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>>
lyap(const Matrix<NX, NX, T>& A, const Matrix<NX, NX, T>& Q) {
    return detail::solve_lyapunov_kron(
        Q, [&](const Matrix<NX, NX, T>& E) { return A * E + E * A.transpose(); }
    );
}

/**
 * @brief Solve the discrete-time Lyapunov (Stein) equation @f$ A X A^\top - X + Q = 0 @f$.
 *
 * The unique solution exists iff no product of eigenvalues satisfies
 * @f$ \lambda_i \lambda_j = 1 @f$. For a Schur-stable @f$ A @f$ (all
 * @f$ |\lambda| < 1 @f$) and symmetric positive-semidefinite @f$ Q @f$,
 * @f$ X = \sum_{k=0}^{\infty} A^k Q (A^\top)^k @f$ — the discrete
 * controllability/observability Gramian when @f$ Q = BB^\top @f$ / @f$ C^\top C @f$.
 *
 * @note Compare with MATLAB®'s @c X=dlyap(A,Q) (which solves @f$ AXA^\top-X+Q=0 @f$).
 *
 * @see lyap() — the continuous-time counterpart.
 *
 * @param A State transition matrix (NX × NX).
 * @param Q Symmetric right-hand side (NX × NX).
 * @return Solution X (NX × NX), or damp::nullopt if some @f$ \lambda_i\lambda_j = 1 @f$.
 */
template<size_t NX, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>>
dlyap(const Matrix<NX, NX, T>& A, const Matrix<NX, NX, T>& Q) {
    return detail::solve_lyapunov_kron(
        Q, [&](const Matrix<NX, NX, T>& E) { return A * E * A.transpose() - E; }
    );
}

} // namespace damp
