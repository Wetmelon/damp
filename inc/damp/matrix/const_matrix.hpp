// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file const_matrix.hpp
 * @brief A matrix known at compile time: its product with a vector is its nonzero terms and nothing else
 */

#include <type_traits>

#include "damp/backend.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/core.hpp"

namespace damp {

namespace detail {

/**
 * @brief The columns of row @p R of @p M whose entry is not its type's zero: their count and their
 *        indices. @p M is anything with rows(), cols() and a constexpr operator(): a matrix, or a
 *        sparsity pattern.
 */
template<const auto& M, size_t R>
inline constexpr auto nonzero_columns = [] {
    using Stored = std::remove_cvref_t<decltype(M)>;
    using Value = std::remove_cvref_t<decltype(M(R, size_t{0}))>;
    struct Columns {
        size_t                              n{0};
        damp::array<size_t, Stored::cols()> col{};
    };
    Columns out{};
    for (size_t c = 0; c < Stored::cols(); ++c) {
        if (M(R, c) != Value{}) {
            out.col[out.n] = c;
            ++out.n;
        }
    }
    return out;
}();

/** @brief One entry of @p M, a constant in any build. */
template<const auto& M, size_t R, size_t C>
inline constexpr auto const_entry = M(R, C);

/**
 * @brief Row @p R of @p M times @p x: a sum over the nonzero entries, with no term for the others.
 *
 * The fold has no initial zero, so no zero is added under IEEE rules, and its arity is one row's
 * nonzero count, never Rows·Cols (Matrix's copy folds per row for the same reason).
 */
template<const auto& M, size_t R, typename T, size_t N>
[[nodiscard, gnu::always_inline]] constexpr T const_row_times(const ColVec<N, T>& x) {
    return [&]<size_t... I>(damp::index_sequence<I...>) -> T {
        if constexpr (sizeof...(I) == 0) {
            return T{0};
        } else {
            return (... + (const_entry<M, R, nonzero_columns<M, R>.col[I]> * x[nonzero_columns<M, R>.col[I]]));
        }
    }(damp::make_index_sequence<nonzero_columns<M, R>.n>{});
}

} // namespace detail

/**
 * @brief A matrix that is a constant of its type: no storage, and a product with a vector that has
 *        no term for an entry that is exactly 0.
 *
 * @p M is an object with static storage, usable in constant expressions. A runtime that holds this
 * in place of a Matrix reads the same, `A * x`, and the compiler sees every coefficient: a 0 is no
 * code, a 1 needs no multiply. Only exact zeros go; an entry of 1e-20 from a matrix exponential is
 * a term like any other. There is deliberately no `value_type`, so the dense `operator*` for
 * MatrixLike operands cannot pick this up: a product is this one or it does not compile.
 *
 * @code
 * static constexpr Matrix<2, 3, float> kA{{1.0f, 0.5f, 0.0f}, {0.0f, 0.0f, 1.0f}};
 * ConstMatrix<kA> A;
 * ColVec<2, float> y = A * x; // y₀ = 1·x₀ + 0.5·x₁, y₁ = x₂
 * @endcode
 *
 * @tparam M The matrix
 */
template<const auto& M>
struct ConstMatrix {
    using Stored = std::remove_cvref_t<decltype(M)>; ///< The matrix's type
    using Scalar = typename Stored::value_type;      ///< Its scalar

    /** @brief Rows. */
    [[nodiscard]] static constexpr size_t rows() { return Stored::rows(); }

    /** @brief Columns. */
    [[nodiscard]] static constexpr size_t cols() { return Stored::cols(); }

    /** @brief One entry. */
    [[nodiscard]] constexpr Scalar operator()(size_t r, size_t c) const { return M(r, c); }

    /** @brief The matrix itself, for anything that reads a Matrix. */
    constexpr operator const Stored&() const { return M; } // NOLINT(google-explicit-constructor): reads as the matrix it is

    /** @brief The nonzero terms of each row, and nothing else. Each store has a constant index, so a small result stays in registers. */
    [[nodiscard, gnu::always_inline]] constexpr ColVec<Stored::rows(), Scalar> operator*(const ColVec<Stored::cols(), Scalar>& x) const {
        ColVec<Stored::rows(), Scalar> out{};
        [&]<size_t... R>(damp::index_sequence<R...>) {
            ((out[R] = detail::const_row_times<M, R>(x)), ...);
        }(damp::make_index_sequence<Stored::rows()>{});
        return out;
    }
};

} // namespace damp
