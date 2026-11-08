// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file pattern_matrix.hpp
 * @brief A matrix whose zeros and ones are known at compile time and whose other entries are runtime values
 */

#include <cstddef>
#include <cstdint>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/const_matrix.hpp"
#include "damp/matrix/core.hpp"

namespace damp {

/** @brief What an entry of a sparsity pattern is known to be at compile time. */
enum class PatternEntry : std::uint8_t {
    zero,  ///< Exactly 0: no term
    one,   ///< Exactly 1: no multiply
    value, ///< Known only at run time: one number in the matrix's storage
};

namespace detail {

/** @brief Refuse a block of a pattern that does not fit inside it: no constant expression at compile time, a trap at run time. */
constexpr void pattern_block_inside(size_t r0, size_t c0, size_t nr, size_t nc, size_t rows, size_t cols) {
    if (r0 > rows || nr > rows - r0 || c0 > cols || nc > cols - c0) {
        ::damp::detail::refuse_initializer_overflow();
    }
}

} // namespace detail

/**
 * @brief Which entries of an R × C matrix are exactly 0, exactly 1, or known only at run time
 *
 * Built at compile time by chaining, starting from all zeros or the identity:
 * @code
 * inline constexpr auto kF = SparsityPattern<4, 4>::identity().value_diagonal(0, 2, 2).values(2, 0, 2, 2);
 * @endcode
 *
 * @tparam R Rows
 * @tparam C Columns
 */
template<size_t R, size_t C>
struct SparsityPattern {
    damp::array<PatternEntry, R * C> entries{}; ///< By row

    /** @brief Rows. */
    [[nodiscard]] static constexpr size_t rows() { return R; }
    /** @brief Columns. */
    [[nodiscard]] static constexpr size_t cols() { return C; }

    /** @brief One entry. */
    [[nodiscard]] constexpr PatternEntry operator()(size_t r, size_t c) const { return entries[(r * C) + c]; }

    /** @brief Ones on the main diagonal, zeros elsewhere. */
    [[nodiscard]] static constexpr SparsityPattern identity() { return SparsityPattern{}.unit_diagonal(0, 0, (R < C) ? R : C); }

    /** @brief Exact ones on @p n entries of a diagonal starting at (@p r0, @p c0). */
    [[nodiscard]] constexpr SparsityPattern unit_diagonal(size_t r0, size_t c0, size_t n) const {
        detail::pattern_block_inside(r0, c0, n, n, R, C);
        SparsityPattern p = *this;
        for (size_t i = 0; i < n; ++i) {
            p.entries[((r0 + i) * C) + c0 + i] = PatternEntry::one;
        }
        return p;
    }

    /** @brief Runtime values on @p n entries of a diagonal starting at (@p r0, @p c0). */
    [[nodiscard]] constexpr SparsityPattern value_diagonal(size_t r0, size_t c0, size_t n) const {
        detail::pattern_block_inside(r0, c0, n, n, R, C);
        SparsityPattern p = *this;
        for (size_t i = 0; i < n; ++i) {
            p.entries[((r0 + i) * C) + c0 + i] = PatternEntry::value;
        }
        return p;
    }

    /** @brief Runtime values on a whole @p nr × @p nc block at (@p r0, @p c0). */
    [[nodiscard]] constexpr SparsityPattern values(size_t r0, size_t c0, size_t nr, size_t nc) const {
        detail::pattern_block_inside(r0, c0, nr, nc, R, C);
        SparsityPattern p = *this;
        for (size_t i = 0; i < nr; ++i) {
            for (size_t j = 0; j < nc; ++j) {
                p.entries[((r0 + i) * C) + c0 + j] = PatternEntry::value;
            }
        }
        return p;
    }

    /** @brief Runtime values off the diagonal of an @p n × @p n block at (@p r0, @p c0), its diagonal left as it is. */
    [[nodiscard]] constexpr SparsityPattern values_off_diagonal(size_t r0, size_t c0, size_t n) const {
        detail::pattern_block_inside(r0, c0, n, n, R, C);
        SparsityPattern p = *this;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if (i != j) {
                    p.entries[((r0 + i) * C) + c0 + j] = PatternEntry::value;
                }
            }
        }
        return p;
    }

    /** @brief How many entries are runtime values. */
    [[nodiscard]] constexpr size_t value_count() const {
        size_t n = 0;
        for (size_t i = 0; i < R * C; ++i) {
            n += (entries[i] == PatternEntry::value) ? size_t{1} : size_t{0};
        }
        return n;
    }

    /** @brief Whether row @p r is exactly the unit vector eᵣ: that row of F·P is that row of P. */
    [[nodiscard]] constexpr bool unit_row(size_t r) const {
        for (size_t c = 0; c < C; ++c) {
            if ((*this)(r, c) != ((c == r) ? PatternEntry::one : PatternEntry::zero)) {
                return false;
            }
        }
        return true;
    }
};

namespace detail {

/** @brief For each entry of @p P, its place among the runtime values. */
template<const auto& P>
inline constexpr auto pattern_slots = [] {
    using Pattern = std::remove_cvref_t<decltype(P)>;
    damp::array<size_t, Pattern::rows() * Pattern::cols()> slots{};
    size_t                                                 n = 0;
    for (size_t i = 0; i < Pattern::rows() * Pattern::cols(); ++i) {
        slots[i] = n;
        n += (P.entries[i] == PatternEntry::value) ? size_t{1} : size_t{0};
    }
    return slots;
}();

/** @brief Whether each row of a square @p P is a unit row. */
template<const auto& P>
inline constexpr auto pattern_unit_rows = [] {
    using Pattern = std::remove_cvref_t<decltype(P)>;
    damp::array<bool, Pattern::rows()> unit{};
    for (size_t r = 0; r < Pattern::rows(); ++r) {
        unit[r] = P.unit_row(r);
    }
    return unit;
}();

/** @brief Whether any row of @p P reads each column. */
template<const auto& P>
inline constexpr auto pattern_read_columns = [] {
    using Pattern = std::remove_cvref_t<decltype(P)>;
    damp::array<bool, Pattern::cols()> read{};
    for (size_t r = 0; r < Pattern::rows(); ++r) {
        for (size_t c = 0; c < Pattern::cols(); ++c) {
            read[c] = read[c] || (P(r, c) != PatternEntry::zero);
        }
    }
    return read;
}();

/**
 * @brief The entries of W = F·S that F S Fᵀ reads, by row: W(i, l) for a unit row l (it is the
 *        result's entry (l, i)), and for the columns of every other row j ≤ i. A unit row of F has no W row.
 */
template<const auto& P>
inline constexpr auto pattern_sandwich_reads = [] {
    constexpr size_t         N = std::remove_cvref_t<decltype(P)>::rows();
    constexpr auto&          unit = pattern_unit_rows<P>;
    damp::array<bool, N * N> reads{};
    for (size_t i = 0; i < N; ++i) {
        if (unit[i]) {
            continue;
        }
        for (size_t l = 0; l < N; ++l) {
            bool read = unit[l];
            for (size_t j = 0; j <= i; ++j) {
                read = read || (!unit[j] && (P(j, l) != PatternEntry::zero));
            }
            reads[(i * N) + l] = read;
        }
    }
    return reads;
}();

} // namespace detail

/**
 * @brief The runtime values of a sparsity pattern, and nothing else: an entry the pattern knows has
 *        no storage, and in a product it has no term (a 0) or no multiply (a 1)
 *
 * For a matrix whose numbers change every tick but whose shape does not: an error-state Jacobian
 * that depends on the attitude, a measurement Jacobian. A matrix whose numbers are known as well
 * is ConstMatrix. There is deliberately no `value_type`, so the dense `operator*` for MatrixLike
 * operands cannot pick this up: a product is this one or it does not compile.
 *
 * @code
 * inline constexpr auto kF = SparsityPattern<2, 2>::identity().value_diagonal(0, 1, 1);
 * PatternMatrix<kF, float> F;
 * F.at<0, 1>() = dt;
 * ColVec<2, float> y = F * x; // y₀ = x₀ + dt·x₁, y₁ = x₁
 * @endcode
 *
 * @tparam P The pattern, with static storage and usable in constant expressions
 * @tparam T Scalar
 */
template<const auto& P, typename T = float>
struct PatternMatrix {
    using Pattern = std::remove_cvref_t<decltype(P)>; ///< The pattern's type
    using Scalar = T;                                 ///< The scalar

    static constexpr size_t NR = Pattern::rows(); ///< Rows
    static constexpr size_t NC = Pattern::cols(); ///< Columns
    static constexpr size_t NV = P.value_count(); ///< Runtime values

    damp::array<T, (NV > 0) ? NV : 1> values{}; ///< The runtime values, by row (one unused number when there are none)

    /** @brief The storage of entry (@p R, @p C), which has to be a runtime value. */
    template<size_t R, size_t C>
    [[nodiscard]] constexpr T& at() {
        static_assert(R < NR && C < NC, "That entry is outside the pattern.");
        static_assert(P(R, C) == PatternEntry::value, "The pattern knows that entry: it has no storage.");
        return values[detail::pattern_slots<P>[(R * NC) + C]];
    }

    /** @brief One entry, read at run time: for checks and display, not for a tick. */
    [[nodiscard]] constexpr T operator()(size_t r, size_t c) const {
        switch (P(r, c)) {
            case PatternEntry::zero:
                return T{0};
            case PatternEntry::one:
                return T{1};
            default:
                return values[detail::pattern_slots<P>[(r * NC) + c]];
        }
    }

    /** @brief The whole matrix, dense: for checks and display, not for a tick. */
    [[nodiscard]] constexpr Matrix<NR, NC, T> to_matrix() const {
        Matrix<NR, NC, T> m{};
        for (size_t r = 0; r < NR; ++r) {
            for (size_t c = 0; c < NC; ++c) {
                m(r, c) = (*this)(r, c);
            }
        }
        return m;
    }

    /** @brief Write the runtime values of a block at (@p R0, @p C0) from @p block: entries the pattern knows are left out. */
    template<size_t R0, size_t C0, size_t BR, size_t BC>
    [[gnu::always_inline]] constexpr void set_block(const Matrix<BR, BC, T>& block) {
        static_assert(R0 + BR <= NR && C0 + BC <= NC, "The block does not fit inside the pattern.");
        detail::for_index_pairs<BR, BC>([&]<size_t I, size_t J>() __attribute__((always_inline)) {
            if constexpr (P(R0 + I, C0 + J) == PatternEntry::value) {
                values[detail::pattern_slots<P>[((R0 + I) * NC) + C0 + J]] = block(I, J);
            }
        });
    }

    /** @brief Entry (@p R, @p C) times @p v: no multiply for a one, and 0 for a zero. */
    template<size_t R, size_t C>
    [[nodiscard, gnu::always_inline]] constexpr T times(T v) const {
        static_assert(R < NR && C < NC, "That entry is outside the pattern.");
        if constexpr (P(R, C) == PatternEntry::zero) {
            return T{0};
        } else if constexpr (P(R, C) == PatternEntry::one) {
            return v;
        } else {
            return values[detail::pattern_slots<P>[(R * NC) + C]] * v;
        }
    }

    /**
     * @brief Row @p R times a vector read through @p read, called with each column as an
     *        integral_constant: a sum of the row's nonzero terms and nothing else
     */
    template<size_t R, typename Read>
    [[nodiscard, gnu::always_inline]] constexpr T row_dot(const Read& read) const {
        return [&]<size_t... I>(damp::index_sequence<I...>) __attribute__((always_inline)) -> T {
            if constexpr (sizeof...(I) == 0) {
                return T{0};
            } else {
                return (... + times<R, detail::nonzero_columns<P, R>.col[I]>(read(std::integral_constant<size_t, detail::nonzero_columns<P, R>.col[I]>{})));
            }
        }(damp::make_index_sequence<detail::nonzero_columns<P, R>.n>{});
    }

    /** @brief This matrix times a vector. */
    [[nodiscard, gnu::always_inline]] constexpr ColVec<NR, T> operator*(const ColVec<NC, T>& x) const {
        ColVec<NR, T> out{};
        [&]<size_t... R>(damp::index_sequence<R...>) __attribute__((always_inline)) {
            ((out[R] = row_dot<R>([&](auto c) __attribute__((always_inline)) { return x[decltype(c)::value]; })), ...);
        }(damp::make_index_sequence<NR>{});
        return out;
    }

    /** @brief This matrix times a dense one, column by column. */
    template<size_t K>
    [[nodiscard, gnu::always_inline]] constexpr Matrix<NR, K, T> operator*(const Matrix<NC, K, T>& m) const {
        Matrix<NR, K, T> out{};
        detail::for_index_pairs<NR, K>([&]<size_t R, size_t L>() __attribute__((always_inline)) {
            out(R, L) = row_dot<R>([&](auto c) __attribute__((always_inline)) { return m(decltype(c)::value, L); });
        });
        return out;
    }
};

/** @brief A dense matrix times Hᵀ for a sparse H: A Hᵀ. */
template<const auto& P, typename T, size_t K, size_t NC>
[[nodiscard, gnu::always_inline]] constexpr Matrix<K, PatternMatrix<P, T>::NR, T> times_transpose(const Matrix<K, NC, T>& a, const PatternMatrix<P, T>& h) {
    static_assert(NC == PatternMatrix<P, T>::NC, "A Hᵀ needs A with as many columns as H.");
    constexpr size_t NR = PatternMatrix<P, T>::NR;
    Matrix<K, NR, T> out{};
    detail::for_index_pairs<K, NR>([&]<size_t I, size_t M>() __attribute__((always_inline)) {
        out(I, M) = h.template row_dot<M>([&](auto c) __attribute__((always_inline)) { return a(I, decltype(c)::value); });
    });
    return out;
}

/**
 * @brief @p S ← F S Fᵀ in place, for a symmetric S: every multiply the pattern can skip, skipped
 *
 * S has to be exactly symmetric, as quadratic_form's X: entries are read from either triangle.
 *
 * W = F·S is computed first, only where it is read, and only then is S written, so no copy of S
 * is needed. A unit row of F has no row of W: it is that row of S. Only the lower triangle is
 * computed, and the result is exactly symmetric. An entry of the result in a unit row l and another
 * row i is W(i, l), which other rows usually read anyway, and an entry in two unit rows is the one
 * already there. That sharing is most of what a symbolic common-subexpression pass finds.
 *
 * @tparam P F's pattern, square
 */
template<const auto& P, size_t N, typename T>
[[gnu::always_inline]] constexpr void quadratic_form_in_place(const PatternMatrix<P, T>& F, Matrix<N, N, T>& S) {
    static_assert(P.rows() == N && P.cols() == N, "F S Fᵀ takes a square F of S's size.");
    constexpr auto& unit = detail::pattern_unit_rows<P>;
    T               w[N][N]; // written before it is read, only where the result reads it
    detail::for_index_pairs<N, N>([&]<size_t I, size_t L>() __attribute__((always_inline)) {
        if constexpr (detail::pattern_sandwich_reads<P>[(I * N) + L]) {
            w[I][L] = F.template row_dot<I>([&](auto k) __attribute__((always_inline)) { return S(decltype(k)::value, L); });
        }
    });
    detail::for_index_pairs<N, N>([&]<size_t I, size_t J>() __attribute__((always_inline)) {
        if constexpr (J <= I && !(unit[I] && unit[J])) {
            T v{};
            if constexpr (unit[I]) {
                v = w[J][I];
            } else if constexpr (unit[J]) {
                v = w[I][J];
            } else {
                v = F.template row_dot<J>([&](auto l) __attribute__((always_inline)) { return w[I][decltype(l)::value]; });
            }
            S(I, J) = v;
            if constexpr (J != I) {
                S(J, I) = v;
            }
        }
    });
}

namespace detail {

/** @brief The sparsity pattern of a constant matrix @p M: its exact zeros and ones, and every other entry a value. */
template<const auto& M>
inline constexpr auto pattern_of = [] {
    using Stored = std::remove_cvref_t<decltype(M)>;
    using Value = typename Stored::value_type;
    SparsityPattern<Stored::rows(), Stored::cols()> p{};
    for (size_t r = 0; r < Stored::rows(); ++r) {
        for (size_t c = 0; c < Stored::cols(); ++c) {
            const Value v = M(r, c);
            p.entries[(r * Stored::cols()) + c] = (v == Value{0}) ? PatternEntry::zero : ((v == Value{1}) ? PatternEntry::one : PatternEntry::value);
        }
    }
    return p;
}();

/**
 * @brief A constant matrix @p M as a PatternMatrix whose values are constants too: the pattern
 *        products and quadratic forms then skip its zeros and ones and read its other entries as literals
 */
template<const auto& M>
inline constexpr auto as_pattern_matrix = [] {
    using Stored = std::remove_cvref_t<decltype(M)>;
    PatternMatrix<pattern_of<M>, typename Stored::value_type> out{};
    for (size_t r = 0; r < Stored::rows(); ++r) {
        for (size_t c = 0; c < Stored::cols(); ++c) {
            if (pattern_of<M>(r, c) == PatternEntry::value) {
                out.values[pattern_slots<pattern_of<M>>[(r * Stored::cols()) + c]] = M(r, c);
            }
        }
    }
    return out;
}();

} // namespace detail

namespace detail {

/** @brief Whether M is a PatternMatrix. */
template<typename M>
inline constexpr bool is_pattern_matrix = false;

template<const auto& P, typename T>
inline constexpr bool is_pattern_matrix<PatternMatrix<P, T>> = true;

} // namespace detail

/** @brief F S Fᵀ for a symmetric S (quadratic_form_in_place, on a copy). */
template<const auto& P, size_t N, typename T>
[[nodiscard, gnu::always_inline]] constexpr Matrix<N, N, T> quadratic_form(const PatternMatrix<P, T>& F, const Matrix<N, N, T>& S) {
    Matrix<N, N, T> out = S;
    quadratic_form_in_place(F, out);
    return out;
}

} // namespace damp
