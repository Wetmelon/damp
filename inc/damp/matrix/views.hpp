// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file views.hpp
 * @brief Matrix view types (transpose, diagonal, triangles)
 */

#include "core.hpp"
#include "damp/backend.hpp"
#include "damp/matrix/matrix_traits.hpp"

namespace damp {

/**
 * @brief Diagonal view of a square matrix
 * @ingroup linear_algebra
 *
 * Provides non-owning access to the diagonal elements of a matrix.
 *
 * @tparam N Matrix dimension
 * @tparam T Element type
 */
template<size_t N, typename T>
struct Diagonal {
private:
    using Parent = damp::conditional_t<damp::is_const_v<T>, const Matrix<N, N, damp::remove_const_t<T>>, Matrix<N, N, damp::remove_const_t<T>>>; ///< const T: a read-only view
    Parent*               mat_ptr{};
    static constexpr bool is_const = damp::is_const_v<T>;

public:
    using value_type = damp::remove_const_t<T>;

    constexpr explicit Diagonal(Parent& mat) : mat_ptr(&mat) {}

    constexpr Diagonal(const Diagonal&) = default;
    constexpr Diagonal(Diagonal&&) noexcept = default;
    constexpr ~Diagonal() = default;

    /** @brief Copies the other diagonal's elements into this one (writes through; never rebinds). */
    constexpr Diagonal& operator=(const Diagonal& other) {
        static_assert(!is_const, "a const view is read-only");
        if (this == &other) {
            return *this;
        }
        for (size_t i = 0; i < N; ++i) {
            (*this)(i) = other(i);
        }
        return *this;
    }

    /** @brief Writes through, as the copy assignment does. */
    constexpr Diagonal& operator=(Diagonal&& other) noexcept {
        *this = static_cast<const Diagonal&>(other);
        return *this;
    }

    [[nodiscard]] static constexpr size_t rows() { return N; }
    [[nodiscard]] static constexpr size_t cols() { return N; }

    [[nodiscard]] constexpr auto& operator()(size_t i) {
        return mat_ptr->data()[(i * N) + i];
    }

    [[nodiscard]] constexpr const value_type& operator()(size_t i) const {
        return mat_ptr->data()[(i * N) + i];
    }

    /**
     * @brief Two-arg access for MatrixLike compliance; returns the element at (r,c)
     * Off-diagonal elements return zero
     */
    [[nodiscard]] constexpr value_type operator()(size_t r, size_t c) const {
        if (r == c) {
            return static_cast<value_type>(mat_ptr->data()[(r * N) + r]);
        }
        return value_type{};
    }

    [[nodiscard]] constexpr const value_type* storage_first() const { return mat_ptr->data(); }
    [[nodiscard]] constexpr const value_type* storage_last() const {
        if constexpr (N == 0) {
            return mat_ptr->data();
        }
        return mat_ptr->data() + ((N * N) - 1);
    }
    [[nodiscard]] static constexpr size_t size() { return N; }
    // Convert to ColVec
    [[nodiscard]] constexpr ColVec<N, damp::remove_const_t<T>> to_vector() const {
        ColVec<N, damp::remove_const_t<T>> result;
        for (size_t i = 0; i < N; ++i) {
            result(i) = (*this)(i);
        }
        return result;
    }

    // Assignment from any MatrixLike with cols()==1 (column vector shape)
    template<MatrixLike M>
        requires(M::rows() == N && M::cols() == 1 && ConvertsExactly<typename M::value_type, damp::remove_const_t<T>>)
    constexpr Diagonal& operator=(const M& vec) {
        for (size_t i = 0; i < N; ++i) {
            (*this)(i) = static_cast<damp::remove_const_t<T>>(vec(i, 0));
        }
        return *this;
    }
};

/**
 * @brief Upper triangular view of a square matrix
 * @ingroup linear_algebra
 *
 * Provides non-owning access to the upper triangular elements of a matrix.
 *
 * @tparam N Matrix dimension
 * @tparam T Element type
 */
template<size_t N, typename T>
struct UpperTriangle {
private:
    const Matrix<N, N, damp::remove_const_t<T>>* mat_ptr{};

public:
    using value_type = damp::remove_const_t<T>;

    constexpr explicit UpperTriangle(const Matrix<N, N, damp::remove_const_t<T>>& mat) : mat_ptr(&mat) {}

    constexpr UpperTriangle(const UpperTriangle&) = default;
    constexpr UpperTriangle(UpperTriangle&&) noexcept = default;
    UpperTriangle& operator=(const UpperTriangle&) = delete; ///< Read-only view: neither rebinds nor writes through
    UpperTriangle& operator=(UpperTriangle&&) = delete;      ///< Read-only view: neither rebinds nor writes through
    constexpr ~UpperTriangle() = default;

    [[nodiscard]] constexpr value_type operator()(size_t r, size_t c) const {
        return r <= c ? mat_ptr->data()[(r * N) + c] : value_type{0};
    }

    [[nodiscard]] static constexpr size_t rows() { return N; }
    [[nodiscard]] static constexpr size_t cols() { return N; }

    [[nodiscard]] constexpr const value_type* storage_first() const { return mat_ptr->data(); }
    [[nodiscard]] constexpr const value_type* storage_last() const {
        if constexpr (N == 0) {
            return mat_ptr->data();
        }
        return mat_ptr->data() + ((N * N) - 1);
    }
};

/**
 * @brief Lower triangular view of a square matrix
 * @ingroup linear_algebra
 *
 * Provides non-owning access to the lower triangular elements of a matrix.
 *
 * @tparam N Matrix dimension
 * @tparam T Element type
 */
template<size_t N, typename T>
struct LowerTriangle {
private:
    const Matrix<N, N, damp::remove_const_t<T>>* mat_ptr{};

public:
    using value_type = damp::remove_const_t<T>;

    constexpr explicit LowerTriangle(const Matrix<N, N, damp::remove_const_t<T>>& mat) : mat_ptr(&mat) {}

    constexpr LowerTriangle(const LowerTriangle&) = default;
    constexpr LowerTriangle(LowerTriangle&&) noexcept = default;
    LowerTriangle& operator=(const LowerTriangle&) = delete; ///< Read-only view: neither rebinds nor writes through
    LowerTriangle& operator=(LowerTriangle&&) = delete;      ///< Read-only view: neither rebinds nor writes through
    constexpr ~LowerTriangle() = default;

    [[nodiscard]] constexpr value_type operator()(size_t r, size_t c) const {
        return r >= c ? mat_ptr->data()[(r * N) + c] : value_type{0};
    }

    [[nodiscard]] static constexpr size_t rows() { return N; }
    [[nodiscard]] static constexpr size_t cols() { return N; }

    [[nodiscard]] constexpr const value_type* storage_first() const { return mat_ptr->data(); }
    [[nodiscard]] constexpr const value_type* storage_last() const {
        if constexpr (N == 0) {
            return mat_ptr->data();
        }
        return mat_ptr->data() + ((N * N) - 1);
    }
};

/**
 * @brief Non-owning row view of a matrix
 *
 * Lightweight slice: pointer to parent + row index. Contiguous in row-major
 * storage. @c MatrixLike (1 × Cols).
 */
template<size_t Rows, size_t Cols, typename T>
struct RowView {
private:
    using Parent = damp::conditional_t<damp::is_const_v<T>, const Matrix<Rows, Cols, damp::remove_const_t<T>>, Matrix<Rows, Cols, damp::remove_const_t<T>>>; ///< const T: a read-only view
    Parent*               mat_ptr{};
    size_t                row_index{};
    static constexpr bool is_const = damp::is_const_v<T>;

public:
    using value_type = damp::remove_const_t<T>;

    constexpr explicit RowView(Parent& mat, size_t row_idx)
        : mat_ptr(&mat), row_index(row_idx) {}

    constexpr RowView(const RowView&) = default;
    constexpr RowView(RowView&&) noexcept = default;
    constexpr ~RowView() = default;

    /** @brief Copies the other row's elements into this one (writes through; never rebinds). */
    constexpr RowView& operator=(const RowView& other) {
        static_assert(!is_const, "a const view is read-only");
        if (this == &other) {
            return *this;
        }
        this->template operator= <RowView>(other);
        return *this;
    }

    /** @brief Writes through, as the copy assignment does. */
    constexpr RowView& operator=(RowView&& other) noexcept {
        *this = static_cast<const RowView&>(other);
        return *this;
    }

    [[nodiscard]] constexpr auto& operator()(size_t c) {
        return mat_ptr->data()[(row_index * Cols) + c];
    }

    [[nodiscard]] constexpr const value_type& operator()(size_t c) const {
        return mat_ptr->data()[(row_index * Cols) + c];
    }

    [[nodiscard]] constexpr auto& operator[](size_t c) { return (*this)(c); }
    [[nodiscard]] constexpr auto& operator[](size_t c) const { return (*this)(c); }

    [[nodiscard]] constexpr auto& operator()(size_t /*r*/, size_t c) {
        return operator()(c);
    }

    [[nodiscard]] constexpr auto& operator()(size_t /*r*/, size_t c) const {
        return operator()(c);
    }

    [[nodiscard]] constexpr auto*    data() { return mat_ptr->data() + (row_index * Cols); }
    [[nodiscard]] constexpr const T* data() const { return mat_ptr->data() + (row_index * Cols); }

    [[nodiscard]] constexpr const value_type* storage_first() const { return mat_ptr->data() + (row_index * Cols); }
    [[nodiscard]] constexpr const value_type* storage_last() const {
        if constexpr (Cols == 0) {
            return storage_first();
        }
        return storage_first() + (Cols - 1);
    }
    [[nodiscard]] static constexpr size_t size() { return Cols; }
    [[nodiscard]] static constexpr size_t rows() { return 1; }
    [[nodiscard]] static constexpr size_t cols() { return Cols; }

    /** @brief Owning @c RowVec copy of this row */
    [[nodiscard]] constexpr RowVec<Cols, value_type> to_vector() const {
        RowVec<Cols, value_type> result;
        for (size_t i = 0; i < Cols; ++i) {
            result(i) = (*this)(i);
        }
        return result;
    }

    /** @brief Owning column of the same length (for matvec that want N×1) */
    [[nodiscard]] constexpr ColVec<Cols, value_type> to_col_vector() const {
        ColVec<Cols, value_type> result;
        for (size_t i = 0; i < Cols; ++i) {
            result[i] = (*this)(i);
        }
        return result;
    }

    /** @brief Owning 1×Cols matrix copy */
    [[nodiscard]] constexpr Matrix<1, Cols, value_type> to_matrix() const {
        Matrix<1, Cols, value_type> result;
        for (size_t i = 0; i < Cols; ++i) {
            result(0, i) = (*this)(i);
        }
        return result;
    }

    /** @brief Assign from a 1×Cols row-shaped @c MatrixLike */
    template<MatrixLike M>
        requires(M::rows() == 1 && M::cols() == Cols && ConvertsExactly<typename M::value_type, value_type>)
    constexpr RowView& operator=(const M& vec) {
        if (detail::may_alias(storage_first(), storage_last(), vec)) {
            *this = Matrix<1, Cols, value_type>(vec);
            return *this;
        }
        for (size_t i = 0; i < Cols; ++i) {
            (*this)(i) = static_cast<value_type>(vec(0, i));
        }
        return *this;
    }

    /** @brief Assign from a Cols×1 column-shaped @c MatrixLike (layout transfer into this row) */
    template<MatrixLike M>
        requires(M::rows() == Cols && M::cols() == 1 && !(M::rows() == 1 && M::cols() == Cols) && ConvertsExactly<typename M::value_type, value_type>)
    constexpr RowView& operator=(const M& col) {
        if (detail::may_alias(storage_first(), storage_last(), col)) {
            *this = Matrix<Cols, 1, value_type>(col); // e.g. A.row(2) = A.col(0)
            return *this;
        }
        for (size_t i = 0; i < Cols; ++i) {
            (*this)(i) = static_cast<value_type>(col(i, 0));
        }
        return *this;
    }
};

/**
 * @brief Non-owning column view of a matrix
 *
 * Lightweight slice: pointer to parent + column index. Strided in row-major
 * storage (strided / non-contiguous). @c MatrixLike (Rows × 1).
 */
template<size_t Rows, size_t Cols, typename T>
struct ColView {
private:
    using Parent = damp::conditional_t<damp::is_const_v<T>, const Matrix<Rows, Cols, damp::remove_const_t<T>>, Matrix<Rows, Cols, damp::remove_const_t<T>>>; ///< const T: a read-only view
    Parent*               mat_ptr{};
    size_t                col_index{};
    static constexpr bool is_const = damp::is_const_v<T>;

public:
    using value_type = damp::remove_const_t<T>;

    constexpr explicit ColView(Parent& mat, size_t col_idx)
        : mat_ptr(&mat), col_index(col_idx) {}

    constexpr ColView(const ColView&) = default;
    constexpr ColView(ColView&&) noexcept = default;
    constexpr ~ColView() = default;

    /** @brief Copies the other column's elements into this one (writes through; never rebinds). */
    constexpr ColView& operator=(const ColView& other) {
        static_assert(!is_const, "a const view is read-only");
        if (this == &other) {
            return *this;
        }
        this->template operator= <ColView>(other);
        return *this;
    }

    /** @brief Writes through, as the copy assignment does. */
    constexpr ColView& operator=(ColView&& other) noexcept {
        *this = static_cast<const ColView&>(other);
        return *this;
    }

    [[nodiscard]] constexpr auto& operator()(size_t r) {
        return mat_ptr->data()[(r * Cols) + col_index];
    }

    [[nodiscard]] constexpr const value_type& operator()(size_t r) const {
        return mat_ptr->data()[(r * Cols) + col_index];
    }

    [[nodiscard]] constexpr auto& operator[](size_t r) { return (*this)(r); }
    [[nodiscard]] constexpr auto& operator[](size_t r) const { return (*this)(r); }

    [[nodiscard]] constexpr auto& operator()(size_t r, size_t /*c*/) {
        return operator()(r);
    }

    [[nodiscard]] constexpr auto& operator()(size_t r, size_t /*c*/) const {
        return operator()(r);
    }

    [[nodiscard]] static constexpr size_t size() { return Rows; }
    [[nodiscard]] static constexpr size_t rows() { return Rows; }
    [[nodiscard]] static constexpr size_t cols() { return 1; }

    [[nodiscard]] constexpr const value_type* storage_first() const { return mat_ptr->data() + col_index; }
    [[nodiscard]] constexpr const value_type* storage_last() const {
        if constexpr (Rows == 0) {
            return storage_first();
        }
        return mat_ptr->data() + ((Rows - 1) * Cols) + col_index;
    }

    /** @brief Owning @c ColVec copy of this column */
    [[nodiscard]] constexpr ColVec<Rows, value_type> to_vector() const {
        ColVec<Rows, value_type> result;
        for (size_t i = 0; i < Rows; ++i) {
            result(i) = (*this)(i);
        }
        return result;
    }

    /** @brief Owning Rows×1 matrix copy */
    [[nodiscard]] constexpr Matrix<Rows, 1, value_type> to_matrix() const {
        Matrix<Rows, 1, value_type> result;
        for (size_t i = 0; i < Rows; ++i) {
            result(i, 0) = (*this)(i);
        }
        return result;
    }

    /** @brief Assign from a Rows×1 column-shaped @c MatrixLike */
    template<MatrixLike M>
        requires(M::rows() == Rows && M::cols() == 1 && ConvertsExactly<typename M::value_type, value_type>)
    constexpr ColView& operator=(const M& vec) {
        if (detail::may_alias(storage_first(), storage_last(), vec)) {
            *this = Matrix<Rows, 1, value_type>(vec);
            return *this;
        }
        for (size_t i = 0; i < Rows; ++i) {
            (*this)(i) = static_cast<value_type>(vec(i, 0));
        }
        return *this;
    }

    /** @brief Assign from a 1×Rows row-shaped @c MatrixLike (layout transfer into this column) */
    template<MatrixLike M>
        requires(M::rows() == 1 && M::cols() == Rows && !(M::rows() == Rows && M::cols() == 1) && ConvertsExactly<typename M::value_type, value_type>)
    constexpr ColView& operator=(const M& row) {
        if (detail::may_alias(storage_first(), storage_last(), row)) {
            *this = Matrix<1, Rows, value_type>(row); // e.g. A.col(2) = A.row(0)
            return *this;
        }
        for (size_t i = 0; i < Rows; ++i) {
            (*this)(i) = static_cast<value_type>(row(0, i));
        }
        return *this;
    }
};

/**
 * @brief Non-owning transpose view of a matrix (zero-copy)
 * @ingroup linear_algebra
 *
 * Provides non-owning access to matrix elements with swapped row/column
 * indices, effectively presenting an M×N matrix as N×M without copying.
 * Satisfies the MatrixLike concept for use with free operators.
 *
 * @tparam Rows Rows of the *parent* matrix
 * @tparam Cols Cols of the *parent* matrix
 * @tparam T    Element type (const-qualified for const access)
 */
template<size_t Rows, size_t Cols, typename T>
struct TransposeView {
private:
    using Parent = damp::conditional_t<damp::is_const_v<T>, const Matrix<Rows, Cols, damp::remove_const_t<T>>, Matrix<Rows, Cols, damp::remove_const_t<T>>>; ///< const T: a read-only view
    Parent*               mat_ptr{};
    static constexpr bool is_const = damp::is_const_v<T>;

public:
    using value_type = damp::remove_const_t<T>;

    constexpr explicit TransposeView(Parent& mat) : mat_ptr(&mat) {}

    constexpr TransposeView(const TransposeView&) = default;
    constexpr TransposeView(TransposeView&&) noexcept = default;
    constexpr ~TransposeView() = default;

    /** @brief Copies the other transpose's elements into this one (writes through; never rebinds). */
    constexpr TransposeView& operator=(const TransposeView& other) {
        static_assert(!is_const, "a const view is read-only");
        if (this == &other) {
            return *this;
        }
        this->template operator= <TransposeView>(other);
        return *this;
    }

    /** @brief Writes through, as the copy assignment does. */
    constexpr TransposeView& operator=(TransposeView&& other) noexcept {
        *this = static_cast<const TransposeView&>(other);
        return *this;
    }

    /** @brief Transposed dimensions: parent Rows×Cols becomes Cols×Rows */
    [[nodiscard]] static constexpr size_t rows() { return Cols; }
    [[nodiscard]] static constexpr size_t cols() { return Rows; }

    /** @brief Access element (r, c) of the transposed view → parent element (c, r) */
    [[nodiscard]] constexpr auto& operator()(size_t r, size_t c) {
        return mat_ptr->data()[(c * Cols) + r];
    }

    [[nodiscard]] constexpr const value_type& operator()(size_t r, size_t c) const {
        return mat_ptr->data()[(c * Cols) + r];
    }

    /** @brief Convert to owning transposed Matrix */
    [[nodiscard]] constexpr Matrix<Cols, Rows, value_type> to_matrix() const {
        Matrix<Cols, Rows, value_type> result;
        for (size_t r = 0; r < Cols; ++r) {
            for (size_t c = 0; c < Rows; ++c) {
                result(r, c) = (*this)(r, c);
            }
        }
        return result;
    }

    [[nodiscard]] constexpr const value_type* storage_first() const { return mat_ptr->data(); }
    [[nodiscard]] constexpr const value_type* storage_last() const {
        if constexpr (Rows == 0 || Cols == 0) {
            return mat_ptr->data();
        }
        return mat_ptr->data() + ((Rows * Cols) - 1);
    }

    /** @brief Assignment from any MatrixLike with transposed dimensions */
    template<MatrixLike M>
        requires(M::rows() == Cols && M::cols() == Rows && ConvertsExactly<typename M::value_type, value_type>)
    constexpr TransposeView& operator=(const M& mat) {
        if constexpr (Rows == 0 || Cols == 0) {
            return *this;
        }
        // A.t() = B writes through. A.t() = A shares this storage, so snapshot.
        const value_type* const dest_first = mat_ptr->data();
        const value_type* const dest_last = storage_last();
        if (detail::storage_overlaps(dest_first, dest_last, mat)) {
            Matrix<Cols, Rows, value_type> tmp{};
            for (size_t r = 0; r < Cols; ++r) {
                for (size_t c = 0; c < Rows; ++c) {
                    tmp(r, c) = static_cast<value_type>(mat(r, c));
                }
            }
            for (size_t r = 0; r < Cols; ++r) {
                for (size_t c = 0; c < Rows; ++c) {
                    (*this)(r, c) = tmp(r, c);
                }
            }
        } else {
            for (size_t r = 0; r < Cols; ++r) {
                for (size_t c = 0; c < Rows; ++c) {
                    (*this)(r, c) = static_cast<value_type>(mat(r, c));
                }
            }
        }
        return *this;
    }
};

} // namespace damp