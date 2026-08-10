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

namespace damp {

namespace detail {
/// Helper for views: returns mutable or const reference depending on IsConst
template<bool IsConst, typename T>
[[nodiscard]] constexpr auto& view_element(const T* data, size_t idx) {
    if constexpr (IsConst) {
        return static_cast<const T&>(data[idx]);
    } else {
        return const_cast<T&>(data[idx]);
    }
}

/// Helper for views: returns mutable or const pointer depending on IsConst
template<bool IsConst, typename T>
[[nodiscard]] constexpr auto* view_data(const T* data) {
    if constexpr (IsConst) {
        return static_cast<const T*>(data);
    } else {
        return const_cast<T*>(data);
    }
}
} // namespace detail

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
    const Matrix<N, N, std::remove_const_t<T>>* mat_ptr;
    static constexpr bool                       is_const = std::is_const_v<T>;

public:
    using value_type = std::remove_const_t<T>;

    constexpr explicit Diagonal(const Matrix<N, N, std::remove_const_t<T>>& mat) : mat_ptr(&mat) {}

    [[nodiscard]] static constexpr size_t rows() { return N; }
    [[nodiscard]] static constexpr size_t cols() { return N; }

    [[nodiscard]] constexpr auto& operator()(size_t i) {
        return detail::view_element<is_const>(mat_ptr->data(), i * N + i);
    }

    [[nodiscard]] constexpr auto& operator()(size_t i) const {
        return static_cast<const T&>(mat_ptr->data()[i * N + i]);
    }

    /// Two-arg access for MatrixLike compliance; returns the element at (r,c)
    /// Off-diagonal elements return zero
    [[nodiscard]] constexpr value_type operator()(size_t r, size_t c) const {
        if (r == c) {
            return static_cast<value_type>(mat_ptr->data()[r * N + r]);
        }
        return value_type{};
    }

    [[nodiscard]] constexpr auto* data() {
        return detail::view_data<is_const>(mat_ptr->data());
    }
    [[nodiscard]] constexpr const T*      data() const { return mat_ptr->data(); }
    [[nodiscard]] static constexpr size_t size() { return N; }
    // Convert to ColVec
    [[nodiscard]] constexpr ColVec<N, std::remove_const_t<T>> to_vector() const {
        ColVec<N, std::remove_const_t<T>> result;
        for (size_t i = 0; i < N; ++i) {
            result(i) = (*this)(i);
        }
        return result;
    }

    // Assignment from any MatrixLike with cols()==1 (column vector shape)
    template<MatrixLike M>
        requires(M::rows() == N && M::cols() == 1)
    constexpr Diagonal& operator=(const M& vec) {
        for (size_t i = 0; i < N; ++i) {
            (*this)(i) = static_cast<std::remove_const_t<T>>(vec(i, 0));
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
    const Matrix<N, N, std::remove_const_t<T>>* mat_ptr;

public:
    using value_type = std::remove_const_t<T>;

    constexpr explicit UpperTriangle(const Matrix<N, N, std::remove_const_t<T>>& mat) : mat_ptr(&mat) {}

    [[nodiscard]] constexpr value_type operator()(size_t r, size_t c) const {
        return r <= c ? mat_ptr->data()[r * N + c] : value_type{0};
    }

    [[nodiscard]] static constexpr size_t rows() { return N; }
    [[nodiscard]] static constexpr size_t cols() { return N; }
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
    const Matrix<N, N, std::remove_const_t<T>>* mat_ptr;

public:
    using value_type = std::remove_const_t<T>;

    constexpr explicit LowerTriangle(const Matrix<N, N, std::remove_const_t<T>>& mat) : mat_ptr(&mat) {}

    [[nodiscard]] constexpr value_type operator()(size_t r, size_t c) const {
        return r >= c ? mat_ptr->data()[r * N + c] : value_type{0};
    }

    [[nodiscard]] static constexpr size_t rows() { return N; }
    [[nodiscard]] static constexpr size_t cols() { return N; }
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
    const Matrix<Rows, Cols, std::remove_const_t<T>>* mat_ptr;
    size_t const                                      row_index;
    static constexpr bool                             is_const = std::is_const_v<T>;

public:
    using value_type = std::remove_const_t<T>;

    constexpr explicit RowView(const Matrix<Rows, Cols, std::remove_const_t<T>>& mat, size_t row_idx)
        : mat_ptr(&mat), row_index(row_idx) {}

    [[nodiscard]] constexpr auto& operator()(size_t c) {
        return detail::view_element<is_const>(mat_ptr->data(), row_index * Cols + c);
    }

    [[nodiscard]] constexpr auto& operator()(size_t c) const {
        return static_cast<const T&>(mat_ptr->data()[row_index * Cols + c]);
    }

    [[nodiscard]] constexpr auto& operator[](size_t c) { return (*this)(c); }
    [[nodiscard]] constexpr auto& operator[](size_t c) const { return (*this)(c); }

    [[nodiscard]] constexpr auto& operator()(size_t /*r*/, size_t c) {
        return operator()(c);
    }

    [[nodiscard]] constexpr auto& operator()(size_t /*r*/, size_t c) const {
        return operator()(c);
    }

    [[nodiscard]] constexpr auto*         data() { return detail::view_data<is_const>(mat_ptr->data() + row_index * Cols); }
    [[nodiscard]] constexpr const T*      data() const { return mat_ptr->data() + row_index * Cols; }
    [[nodiscard]] static constexpr size_t size() { return Cols; }
    [[nodiscard]] static constexpr size_t rows() { return 1; }
    [[nodiscard]] static constexpr size_t cols() { return Cols; }

    /// Owning @c RowVec copy of this row
    [[nodiscard]] constexpr RowVec<Cols, value_type> to_vector() const {
        RowVec<Cols, value_type> result;
        for (size_t i = 0; i < Cols; ++i) {
            result(i) = (*this)(i);
        }
        return result;
    }

    /// Owning column of the same length (for matvec that want N×1)
    [[nodiscard]] constexpr ColVec<Cols, value_type> to_col_vector() const {
        ColVec<Cols, value_type> result;
        for (size_t i = 0; i < Cols; ++i) {
            result[i] = (*this)(i);
        }
        return result;
    }

    /// Owning 1×Cols matrix copy
    [[nodiscard]] constexpr Matrix<1, Cols, value_type> to_matrix() const {
        Matrix<1, Cols, value_type> result;
        for (size_t i = 0; i < Cols; ++i) {
            result(0, i) = (*this)(i);
        }
        return result;
    }

    /// Assign from a 1×Cols row-shaped @c MatrixLike
    template<MatrixLike M>
        requires(M::rows() == 1 && M::cols() == Cols)
    constexpr RowView& operator=(const M& vec) {
        for (size_t i = 0; i < Cols; ++i) {
            (*this)(i) = static_cast<value_type>(vec(0, i));
        }
        return *this;
    }

    /// Assign from a Cols×1 column-shaped @c MatrixLike (layout transfer into this row)
    template<MatrixLike M>
        requires(M::rows() == Cols && M::cols() == 1 && !(M::rows() == 1 && M::cols() == Cols))
    constexpr RowView& operator=(const M& col) {
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
 * storage (not a contiguous span). @c MatrixLike (Rows × 1).
 */
template<size_t Rows, size_t Cols, typename T>
struct ColView {
private:
    const Matrix<Rows, Cols, std::remove_const_t<T>>* mat_ptr;
    size_t const                                      col_index;
    static constexpr bool                             is_const = std::is_const_v<T>;

public:
    using value_type = std::remove_const_t<T>;

    constexpr explicit ColView(const Matrix<Rows, Cols, std::remove_const_t<T>>& mat, size_t col_idx)
        : mat_ptr(&mat), col_index(col_idx) {}

    [[nodiscard]] constexpr auto& operator()(size_t r) {
        return detail::view_element<is_const>(mat_ptr->data(), r * Cols + col_index);
    }

    [[nodiscard]] constexpr auto& operator()(size_t r) const {
        return static_cast<const T&>(mat_ptr->data()[r * Cols + col_index]);
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

    /// Owning @c ColVec copy of this column
    [[nodiscard]] constexpr ColVec<Rows, value_type> to_vector() const {
        ColVec<Rows, value_type> result;
        for (size_t i = 0; i < Rows; ++i) {
            result(i) = (*this)(i);
        }
        return result;
    }

    /// Owning Rows×1 matrix copy
    [[nodiscard]] constexpr Matrix<Rows, 1, value_type> to_matrix() const {
        Matrix<Rows, 1, value_type> result;
        for (size_t i = 0; i < Rows; ++i) {
            result(i, 0) = (*this)(i);
        }
        return result;
    }

    /// Assign from a Rows×1 column-shaped @c MatrixLike
    template<MatrixLike M>
        requires(M::rows() == Rows && M::cols() == 1)
    constexpr ColView& operator=(const M& vec) {
        for (size_t i = 0; i < Rows; ++i) {
            (*this)(i) = static_cast<value_type>(vec(i, 0));
        }
        return *this;
    }

    /// Assign from a 1×Rows row-shaped @c MatrixLike (layout transfer into this column)
    template<MatrixLike M>
        requires(M::rows() == 1 && M::cols() == Rows && !(M::rows() == Rows && M::cols() == 1))
    constexpr ColView& operator=(const M& row) {
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
    const Matrix<Rows, Cols, std::remove_const_t<T>>* mat_ptr;
    static constexpr bool                             is_const = std::is_const_v<T>;

public:
    using value_type = std::remove_const_t<T>;

    constexpr explicit TransposeView(const Matrix<Rows, Cols, std::remove_const_t<T>>& mat) : mat_ptr(&mat) {}

    /// Transposed dimensions: parent Rows×Cols becomes Cols×Rows
    [[nodiscard]] static constexpr size_t rows() { return Cols; }
    [[nodiscard]] static constexpr size_t cols() { return Rows; }

    /// Access element (r, c) of the transposed view → parent element (c, r)
    [[nodiscard]] constexpr auto& operator()(size_t r, size_t c) {
        return detail::view_element<is_const>(mat_ptr->data(), c * Cols + r);
    }

    [[nodiscard]] constexpr auto& operator()(size_t r, size_t c) const {
        return static_cast<const T&>(mat_ptr->data()[c * Cols + r]);
    }

    /// Convert to owning transposed Matrix
    [[nodiscard]] constexpr Matrix<Cols, Rows, value_type> to_matrix() const {
        Matrix<Cols, Rows, value_type> result;
        for (size_t r = 0; r < Cols; ++r) {
            for (size_t c = 0; c < Rows; ++c) {
                result(r, c) = (*this)(r, c);
            }
        }
        return result;
    }

    /// Assignment from any MatrixLike with transposed dimensions
    template<MatrixLike M>
        requires(M::rows() == Cols && M::cols() == Rows)
    constexpr TransposeView& operator=(const M& mat) {
        for (size_t r = 0; r < Cols; ++r) {
            for (size_t c = 0; c < Rows; ++c) {
                (*this)(r, c) = static_cast<value_type>(mat(r, c));
            }
        }
        return *this;
    }
};

} // namespace damp