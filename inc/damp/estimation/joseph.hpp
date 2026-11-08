// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file joseph.hpp
 * @brief The Joseph-form covariance update for a measurement matrix whose zeros and ones are known
 *        at compile time, shared by the filters that run on one
 */

#include "damp/backend.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/pattern_matrix.hpp"
#include "damp/matrix/solve.hpp"

namespace damp::detail {

/**
 * @brief P ← (I − K H) P (I − K H)ᵀ + K R Kᵀ with K = P Hᵀ S⁻¹ and S = H P Hᵀ + R, for an H whose
 *        zeros and ones are known at compile time
 *
 * The sum of two positive semidefinite terms, not the shorter P − K H P, which can lose symmetry and
 * positive definiteness to rounding. S is solved by Cholesky. Only the entries of (I − K H) P that
 * the result reads are computed: its lower triangle, and past the diagonal the columns H reads. The
 * result is exactly symmetric.
 *
 * @param P Covariance, updated in place
 * @param H Measurement matrix on its pattern; its values may be constants or runtime values
 * @param R Covariance of the reading's noise
 * @param K The gain, every entry written on success (it may start unfilled)
 * @return false if S is not positive definite: P and K are unchanged
 *
 * @see "Optimal State Estimation" (Simon, 2006), §5.1, §7.6.2 (Joseph form)
 */
template<const auto& Pattern, size_t NX, size_t NY, typename T>
[[nodiscard]] constexpr bool joseph_update(Matrix<NX, NX, T>& P, const PatternMatrix<Pattern, T>& H, const Matrix<NY, NY, T>& R, Matrix<NX, NY, T>& K) {
    static_assert(PatternMatrix<Pattern, T>::NR == NY && PatternMatrix<Pattern, T>::NC == NX, "H is NY × NX.");
    const Matrix<NY, NX, T> HP = H * P;
    const Matrix<NY, NY, T> S = times_transpose(HP, H) + R;
    const auto              K_opt = mat::cholesky_solve(S, HP);
    if (!K_opt) {
        return false;
    }
    const auto& Kt = K_opt.value();
    for_each_entry<NX, NY>([&](size_t i, size_t m) { K(i, m) = Kt(m, i); }); // the transpose, in place: no temporary to copy

    const Matrix<NX, NY, T> KR = K * R;
    constexpr auto&         read_by_H = pattern_read_columns<Pattern>;
    // A P with A = I − K H, row by row: its lower triangle, and past the diagonal only the columns H reads.
    T AP[NX][NX]; // written before it is read, only where the result reads it
    for (size_t i = 0; i < NX; ++i) {
        const auto ap_entry = [&](size_t l) __attribute__((always_inline)) {
            T s = P(i, l);
            [&]<size_t... M>(damp::index_sequence<M...>) __attribute__((always_inline)) {
                ((s -= K(i, M) * HP(M, l)), ...);
            }(damp::make_index_sequence<NY>{});
            AP[i][l] = s;
        };
        for (size_t l = 0; l <= i; ++l) {
            ap_entry(l);
        }
        [&]<size_t... C>(damp::index_sequence<C...>) __attribute__((always_inline)) {
            ((read_by_H[C] && C > i ? ap_entry(C) : void()), ...);
        }(damp::make_index_sequence<NX>{});
    }
    // P ← A P Aᵀ + K R Kᵀ, row by row: the part of each row that does not depend on the column once.
    for (size_t i = 0; i < NX; ++i) {
        T d[NY];
        [&]<size_t... M>(damp::index_sequence<M...>) __attribute__((always_inline)) {
            ((d[M] = KR(i, M) - H.template row_dot<M>([&](auto c) __attribute__((always_inline)) { return AP[i][decltype(c)::value]; })), ...);
        }(damp::make_index_sequence<NY>{});
        for (size_t j = 0; j <= i; ++j) {
            T s = AP[i][j];
            [&]<size_t... M>(damp::index_sequence<M...>) __attribute__((always_inline)) {
                ((s += K(j, M) * d[M]), ...);
            }(damp::make_index_sequence<NY>{});
            P(i, j) = s;
            P(j, i) = s;
        }
    }
    return true;
}

} // namespace damp::detail
