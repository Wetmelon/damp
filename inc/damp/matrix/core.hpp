// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file core.hpp
 * @brief Fixed-size stack-allocated matrix type
 */

#include <cstddef>
#include <initializer_list>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "matrix_traits.hpp"

namespace damp {

template<typename T>
struct is_matrix_type : std::false_type {};

/* Concrete (owning) types */
template<size_t Rows, size_t Cols, typename T>
struct Matrix;

template<size_t N, typename T>
struct ColVec;

template<size_t N, typename T>
struct RowVec;

/* Non-owning Views */
template<size_t N, typename T>
struct Diagonal;

template<size_t N, typename T>
struct UpperTriangle;

template<size_t N, typename T>
struct LowerTriangle;

template<size_t Rows, size_t Cols, size_t ParentCols, typename T>
struct Block;

template<size_t Rows, size_t Cols, typename T>
struct RowView;

template<size_t Rows, size_t Cols, typename T>
struct ColView;

template<size_t Rows, size_t Cols, typename T>
struct TransposeView;

/**
 * @ingroup linear_algebra
 * @brief Fixed-size, stack-allocated matrix for linear algebra operations
 *
 * Provides efficient matrix arithmetic, multiplication, transposition, trace, and inversion.
 * Supports constexpr evaluation for compile-time matrix computations.
 *
 * @tparam Rows Number of rows
 * @tparam Cols Number of columns
 * @tparam T Element type — only @c float, @c double, or @c damp::complex of those
 *         (see @ref is_matrix_element_v). Nested @c Matrix / block types are
 *         not elements: use a flat matrix + @c .block<>() or multi-axis
 *         @c Lut1D/@c Lut2D/@c Lut3D for multi-parameter tables.
 */
template<size_t Rows, size_t Cols, typename T = double>
struct Matrix {
protected:
    damp::array<T, Rows * Cols> data_{};

    template<size_t, size_t, typename>
    friend struct Matrix;

    /// One row of constant-index stores (fold arity = Cols, not Rows*Cols).
    template<size_t R>
    constexpr void copy_row_from(const Matrix& other) {
        [&]<size_t... C>(damp::index_sequence<C...>) {
            ((data_[(R * Cols) + C] = other.data_[(R * Cols) + C]), ...);
        }(damp::make_index_sequence<Cols>{});
    }

    /**
     * @brief Element-wise copy via per-row folds
     *
     * Each store uses a compile-time constant index so GCC can SROA small local
     * matrices. Fold arity is Cols (one row), with an outer pack over Rows —
     * not a single pack of length Rows*Cols.
     */
    constexpr void copy_from(const Matrix& other) {
        [&]<size_t... R>(damp::index_sequence<R...>) {
            (copy_row_from<R>(other), ...);
        }(damp::make_index_sequence<Rows>{});
    }

public:
    typedef T value_type;

    static_assert(
        is_matrix_element_v<T>,
        "Matrix element type must be float, double, or damp::complex<float|double> "
        "(not nested Matrix or other types — use flat Matrix + block views, "
        "or Lut1D/Lut2D/Lut3D for multi-parameter lookup tables)"
    );

    /**
     * @brief Default constructor, initializes all elements to zero
     */
    constexpr Matrix() = default;

    /**
     * @brief Copy constructor and assignment operator
     *
     * Copies per row with a Cols-wide index_sequence fold (constant indices for
     * GCC SROA of small locals). Arity is Cols rather than Rows*Cols, which
     * keeps g++ frontend cost reasonable for typical matrix shapes.
     */
    constexpr Matrix(const Matrix& other) { copy_from(other); }

    constexpr Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            copy_from(other);
        }
        return *this;
    }

    /**
     * @brief Move constructor and assignment operator
     */
    constexpr Matrix(Matrix&&) = default;
    constexpr Matrix& operator=(Matrix&&) = default;

    /**
     * @brief Destructor
     */
    constexpr ~Matrix() = default;

    /**
     * @brief Type conversion assignment operator (copy from different element type)
     */
    template<typename U>
    constexpr Matrix& operator=(const Matrix<Rows, Cols, U>& other) {
        for (size_t r = 0; r < Rows; ++r) {
            for (size_t c = 0; c < Cols; ++c) {
                data_[(r * Cols) + c] = static_cast<T>(other.data_[(r * Cols) + c]);
            }
        }
        return *this;
    }

    /**
     * @brief Type conversion assignment operator (move from different element type)
     */
    template<typename U>
    constexpr Matrix& operator=(Matrix<Rows, Cols, U>&& other) {
        for (size_t r = 0; r < Rows; ++r) {
            for (size_t c = 0; c < Cols; ++c) {
                data_[(r * Cols) + c] = static_cast<T>(other.data_[(r * Cols) + c]);
            }
        }
        return *this;
    }

    /**
     * @brief Assignment from nested initializer list
     */
    constexpr Matrix& operator=(std::initializer_list<std::initializer_list<T>> init) {
        size_t r = 0;
        for (const auto& row : init) {
            size_t c = 0;
            for (const auto& val : row) {
                if (r < Rows && c < Cols) {
                    data_[(r * Cols) + c] = val;
                }
                ++c;
            }
            ++r;
        }
        return *this;
    }

    /**
     * @brief Type conversion constructor from another matrix with different element type
     */
    template<typename U>
    explicit constexpr Matrix(const Matrix<Rows, Cols, U>& other) : Matrix() {
        for (size_t r = 0; r < Rows; ++r) {
            for (size_t c = 0; c < Cols; ++c) {
                data_[(r * Cols) + c] = static_cast<T>(other.data_[(r * Cols) + c]);
            }
        }
    }

    /**
     * @brief Constructor from nested initializer list
     */
    constexpr Matrix(std::initializer_list<std::initializer_list<T>> init) : Matrix() {
        size_t r = 0;
        for (const auto& row : init) {
            size_t c = 0;
            for (const auto& val : row) {
                if (r < Rows && c < Cols) {
                    data_[(r * Cols) + c] = val;
                }
                ++c;
            }
            ++r;
        }
    }

    /**
     * @brief Constructor from damp::array, enables class template argument deduction
     */
    constexpr Matrix(const damp::array<T, Rows * Cols>& arr) : Matrix() {
        for (size_t i = 0; i < Rows * Cols; ++i) {
            data_[i] = arr[i];
        }
    }

    /**
     * @brief Returns a view over a sub-block of the matrix
     *
     * @param start_row Starting row index of the block
     * @param start_col Starting column index of the block
     * @tparam Brows Number of rows in the block
     * @tparam Bcols Number of columns in the block
     *
     * @return a Block view of the specified submatrix
     */
    template<size_t Brows, size_t Bcols>
    constexpr Block<Brows, Bcols, Cols, T> block(size_t start_row, size_t start_col) {
        if constexpr (Brows == 0 || Bcols == 0) {
            return Block<Brows, Bcols, Cols, T>{nullptr, 0};
        } else {
            return Block<Brows, Bcols, Cols, T>{data_.data(), (start_row * Cols) + start_col};
        }
    }

    template<size_t Brows, size_t Bcols>
    constexpr Block<Brows, Bcols, Cols, const T> block(size_t start_row, size_t start_col) const {
        if constexpr (Brows == 0 || Bcols == 0) {
            return Block<Brows, Bcols, Cols, const T>{nullptr, 0};
        } else {
            return Block<Brows, Bcols, Cols, const T>{data_.data(), (start_row * Cols) + start_col};
        }
    }

    /**
     * @brief Returns a view of the diagonal elements
     *
     * @return Diagonal view of the matrix
     */
    constexpr Diagonal<Rows, T> diagonal()
        requires(Rows == Cols)
    {
        return Diagonal<Rows, T>{*this};
    }

    constexpr Diagonal<Rows, const T> diagonal() const
        requires(Rows == Cols)
    {
        return Diagonal<Rows, const T>{*this};
    }

    /**
     * @brief Returns a view of the upper triangular elements
     *
     * @return UpperTriangle view of the matrix
     */
    constexpr UpperTriangle<Rows, T> upper_triangle() const
        requires(Rows == Cols)
    {
        return UpperTriangle<Rows, T>{*this};
    }

    /**
     * @brief Returns a view of the lower triangular elements
     *
     * @return LowerTriangle view of the matrix
     */
    constexpr LowerTriangle<Rows, T> lower_triangle() const
        requires(Rows == Cols)
    {
        return LowerTriangle<Rows, T>{*this};
    }

    /**
     * @brief Returns a view of the specified row
     *
     * @param row_index The index of the row to view
     * @return RowView of the specified row
     */
    constexpr RowView<Rows, Cols, T> row(size_t row_index) {
        return RowView<Rows, Cols, T>{*this, row_index};
    }

    constexpr RowView<Rows, Cols, const T> row(size_t row_index) const {
        return RowView<Rows, Cols, const T>{*this, row_index};
    }

    /**
     * @brief Returns a view of the specified column
     *
     * @param col_index The index of the column to view
     * @return ColView of the specified column
     */
    constexpr ColView<Rows, Cols, T> col(size_t col_index) {
        return ColView<Rows, Cols, T>{*this, col_index};
    }

    constexpr ColView<Rows, Cols, const T> col(size_t col_index) const {
        return ColView<Rows, Cols, const T>{*this, col_index};
    }

    /**
     * @brief Constructor from a compatible Block (copies data)
     */
    template<size_t ParentCols, typename U>
    constexpr Matrix(const Block<Rows, Cols, ParentCols, U>& block) : Matrix() {
        for (size_t r = 0; r < Rows; ++r) {
            for (size_t c = 0; c < Cols; ++c) {
                data_[(r * Cols) + c] = static_cast<T>(block(r, c));
            }
        }
    }

    /**
     * @brief Constructor from any MatrixLike type with matching dimensions (copies data)
     */
    template<MatrixLike M>
        requires(M::rows() == Rows && M::cols() == Cols && !std::is_same_v<std::remove_cvref_t<M>, Matrix<Rows, Cols, T>>)
    constexpr Matrix(const M& other) : Matrix() {
        for (size_t r = 0; r < Rows; ++r) {
            for (size_t c = 0; c < Cols; ++c) {
                data_[(r * Cols) + c] = static_cast<T>(other(r, c));
            }
        }
    }

    /**
     * @brief Assignment from any MatrixLike type with matching dimensions
     */
    template<MatrixLike M>
        requires(M::rows() == Rows && M::cols() == Cols)
    constexpr Matrix& operator=(const M& other) {
        for (size_t r = 0; r < Rows; ++r) {
            for (size_t c = 0; c < Cols; ++c) {
                data_[(r * Cols) + c] = static_cast<T>(other(r, c));
            }
        }
        return *this;
    }

    /**
     * @brief Create a zero matrix
     */
    [[nodiscard]] static constexpr Matrix zeros() { return Matrix{}; }

    /**
     * @brief Create a matrix with all elements set to a constant value
     */
    [[nodiscard]] static constexpr Matrix constant(const T& value) {
        Matrix result;
        result.data_.fill(value);
        return result;
    }

    /**
     * @brief  Create an identity matrix (square matrix with ones on the diagonal)
     */
    [[nodiscard]] static constexpr Matrix identity()
        requires(Rows == Cols)
    {
        Matrix result{};
        for (size_t r = 0; r < Rows; ++r) {
            result.data_[(r * Cols) + r] = T{1};
        }
        return result;
    }

    /**
     * @brief Create a diagonal matrix from an array of diagonal elements
     *
     * @param diag Array containing the diagonal elements
     * @return Diagonal matrix with the specified diagonal elements
     */
    [[nodiscard]] static constexpr Matrix diagonal(const damp::array<T, Rows>& diag)
        requires(Rows == Cols)
    {
        Matrix result{};
        for (size_t r = 0; r < Rows; ++r) {
            result.data_[(r * Cols) + r] = diag[r];
        }
        return result;
    }

    /**
     * @brief Fill all elements with a constant value
     */
    constexpr void fill(const T& value) {
        data_.fill(value);
    }

    /**
     * @brief Type conversion method for converting to a different element type
     */
    template<typename U>
    [[nodiscard]] constexpr Matrix<Rows, Cols, U> as() const {
        return Matrix<Rows, Cols, U>(*this);
    }

    /**
     * @brief Get number of rows
     */
    [[nodiscard]] static constexpr size_t rows() { return Rows; }

    /**
     * @brief Get number of columns
     */
    [[nodiscard]] static constexpr size_t cols() { return Cols; }

    /**
     * @brief Get const pointer to data in row-major order
     */
    [[nodiscard]] constexpr const T* data() const { return data_.data(); }

    /**
     * @brief Get pointer to data in row-major order
     */
    [[nodiscard]] constexpr T* data() { return data_.data(); }

    /**
     * @brief Iterator begin for range-based access
     */
    [[nodiscard]] constexpr const T* begin() const { return data(); }

    /**
     * @brief Iterator end for range-based access
     */
    [[nodiscard]] constexpr const T* end() const { return data() + (Rows * Cols); }

    /**
     * @brief Iterator begin for range-based access
     */
    [[nodiscard]] constexpr T* begin() { return data(); }

    /**
     * @brief Iterator end for range-based access
     */
    [[nodiscard]] constexpr T* end() { return data() + (Rows * Cols); }

    /**
     * @brief Element access operator
     * @param row Row index
     * @param col Column index
     * @return Value of element at (row, col)
     */
    constexpr T operator()(size_t row, size_t col) const { return data_[(row * Cols) + col]; }

    /**
     * @brief Element access operator
     * @param row Row index
     * @param col Column index
     * @return Reference to element at (row, col)
     */
    constexpr T& operator()(size_t row, size_t col) { return data_[(row * Cols) + col]; }

    /**
     * @brief Equality comparison operator
     */
    [[nodiscard]] constexpr bool operator==(const Matrix& other) const {
        for (size_t r = 0; r < Rows; ++r) {
            for (size_t c = 0; c < Cols; ++c) {
                if (data_[(r * Cols) + c] != other.data_[(r * Cols) + c]) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * @brief Inequality comparison operator
     */
    [[nodiscard]] constexpr bool operator!=(const Matrix& other) const { return !(*this == other); }

    /**
     * @brief In-place addition with any MatrixLike type
     */
    template<MatrixLike M>
        requires(M::rows() == Rows && M::cols() == Cols)
    constexpr Matrix& operator+=(const M& other) {
        for (size_t r = 0; r < Rows; ++r) {
            for (size_t c = 0; c < Cols; ++c) {
                data_[(r * Cols) + c] += static_cast<T>(other(r, c));
            }
        }
        return *this;
    }

    /**
     * @brief In-place subtraction with any MatrixLike type
     */
    template<MatrixLike M>
        requires(M::rows() == Rows && M::cols() == Cols)
    constexpr Matrix& operator-=(const M& other) {
        for (size_t r = 0; r < Rows; ++r) {
            for (size_t c = 0; c < Cols; ++c) {
                data_[(r * Cols) + c] -= static_cast<T>(other(r, c));
            }
        }
        return *this;
    }

    /**
     * @brief In-place addition with scalar
     */
    constexpr Matrix& operator+=(T scalar) {
        for (size_t i = 0; i < Rows * Cols; ++i) {
            data_[i] += scalar;
        }
        return *this;
    }

    /**
     * @brief In-place subtraction with scalar
     */
    constexpr Matrix& operator-=(T scalar) {
        for (size_t i = 0; i < Rows * Cols; ++i) {
            data_[i] -= scalar;
        }
        return *this;
    }

    /**
     * @brief In-place multiplication with scalar
     */
    constexpr Matrix& operator*=(T scalar) {
        for (size_t i = 0; i < Rows * Cols; ++i) {
            data_[i] *= scalar;
        }
        return *this;
    }

    /**
     * @brief In-place division with scalar
     */
    constexpr Matrix& operator/=(T scalar) {
        for (size_t i = 0; i < Rows * Cols; ++i) {
            data_[i] /= scalar;
        }
        return *this;
    }

    /**
     * @brief Matrix transpose (concrete copy)
     * @return Transposed matrix (Cols × Rows)
     * @see t() for a zero-copy TransposeView
     */
    [[nodiscard]] constexpr Matrix<Cols, Rows, T> transpose() const {
        Matrix<Cols, Rows, T> result;
        for (size_t r = 0; r < Rows; ++r) {
            for (size_t c = 0; c < Cols; ++c) {
                result.data_[(c * Rows) + r] = data_[(r * Cols) + c];
            }
        }
        return result;
    }

    /**
     * @brief Non-owning transpose view (zero-copy)
     * @return TransposeView that swaps row/column indexing
     */
    [[nodiscard]] constexpr TransposeView<Rows, Cols, T> t() {
        return TransposeView<Rows, Cols, T>{*this};
    }

    [[nodiscard]] constexpr TransposeView<Rows, Cols, const T> t() const {
        return TransposeView<Rows, Cols, const T>{*this};
    }

    /**
     * @brief Conjugate transpose (Hermitian adjoint, concrete copy)
     *
     * Returns Aᴴ where Aᴴ(i,j) = conj(A(j,i)). For real types, identical
     * to transpose().
     *
     * @return Conjugate-transposed matrix (Cols × Rows)
     */
    [[nodiscard]] constexpr Matrix<Cols, Rows, T> conjugate_transpose() const {
        Matrix<Cols, Rows, T> result;
        for (size_t r = 0; r < Rows; ++r) {
            for (size_t c = 0; c < Cols; ++c) {
                result.data_[(c * Rows) + r] = damp::conj(data_[(r * Cols) + c]);
            }
        }
        return result;
    }

    /**
     * @brief Reshape the matrix to new dimensions (total elements must match)
     * @tparam NewRows New number of rows
     * @tparam NewCols New number of columns
     * @return Reshaped matrix
     */
    template<size_t NewRows, size_t NewCols>
    [[nodiscard]] constexpr Matrix<NewRows, NewCols, T> reshape() const
        requires(NewRows* NewCols == Rows * Cols)
    {
        Matrix<NewRows, NewCols, T> result;
        for (size_t i = 0; i < Rows * Cols; ++i) {
            size_t r_new = i / NewCols;
            size_t c_new = i % NewCols;
            size_t r_old = i / Cols;
            size_t c_old = i % Cols;
            result(r_new, c_new) = (*this)(r_old, c_old);
        }
        return result;
    }

    /**
     * @brief Sum of diagonal elements (square matrices only)
     *
     * @return Trace of the matrix
     */
    [[nodiscard]] constexpr T trace() const
        requires(Rows == Cols)
    {
        T accum = T{0};
        for (size_t i = 0; i < Rows; ++i) {
            accum += data_[(i * Cols) + i];
        }
        return accum;
    }

    /**
     * @brief Frobenius norm (generalized "length") of the matrix
     *
     * The Frobenius norm is the square root of the sum of the absolute squares
     * of every element: ||A||_F = sqrt( sum |a_ij|^2 ).
     *
     * For a column vector this is the Euclidean length. For a matrix it is
     * a measure of overall magnitude, analogous to the Euclidean norm but
     * applied to all entries.
     *
     * Always returns a real value, even for complex-valued matrices.
     *
     * @return Frobenius norm (real scalar)
     */
    [[nodiscard]] constexpr scalar_type_t<T> norm() const {
        using real_t = scalar_type_t<T>;
        real_t sum_sq = real_t{0};
        for (size_t r = 0; r < Rows; ++r) {
            for (size_t c = 0; c < Cols; ++c) {
                auto abs_val = damp::abs(data_[(r * Cols) + c]);
                sum_sq += abs_val * abs_val;
            }
        }
        return damp::sqrt(sum_sq);
    }

    /**
     * @brief Frobenius inner product with another matrix
     *
     * For real matrices:   dot(A,B) = sum( a_ij * b_ij )
     * For complex matrices: dot(A,B) = sum( conj(a_ij) * b_ij )  (=trace(A^H B))
     *
     * @param other Matrix to compute inner product with
     * @return Scalar inner product
     */
    [[nodiscard]] constexpr T dot(const Matrix& other) const {
        T result = T{0};
        for (size_t r = 0; r < Rows; ++r) {
            for (size_t c = 0; c < Cols; ++c) {
                result += damp::conj(data_[(r * Cols) + c]) * other.data_[(r * Cols) + c];
            }
        }
        return result;
    }

    /**
     * @brief Sum of all elements in the matrix
     * @return Sum of all elements
     */
    [[nodiscard]] constexpr T sum() const {
        T result = T{0};
        for (size_t r = 0; r < Rows; ++r) {
            for (size_t c = 0; c < Cols; ++c) {
                result += data_[(r * Cols) + c];
            }
        }
        return result;
    }

    /**
     * @brief Mean (average) of all elements in the matrix
     * @return Mean of all elements
     */
    [[nodiscard]] constexpr T mean() const {
        return sum() / static_cast<T>(Rows * Cols);
    }

    /**
     * @brief Matrix inverse @f$ A^{-1} @f$ via mat::solve(@f$ A, I @f$)
     *
     * Solves @f$ A X = I @f$ (Cholesky when Hermitian PD, else LU with partial
     * pivoting). Prefer mat::solve / @ref mat::cholesky_solve when the
     * inverse is only an intermediate (@f$ A^{-1} b @f$) — this method is for
     * when the inverse matrix itself is the deliverable.
     *
     * @return Inverse if nonsingular, damp::nullopt if singular
     * @note Compare with MATLAB®'s inv(A).
     * @see mat::solve()
     */
    [[nodiscard]] constexpr damp::optional<Matrix> inverse() const
        requires(Rows == Cols);

    /**
     * @brief Extract a row as a 1×Cols matrix
     * @param r_idx Row index
     * @return Row vector
     */
    [[nodiscard]] constexpr Matrix<1, Cols, T> row_vector(size_t r_idx) const {
        Matrix<1, Cols, T> result;
        for (size_t c = 0; c < Cols; ++c) {
            result.data_[c] = data_[(r_idx * Cols) + c];
        }
        return result;
    }

    /**
     * @brief Extract a column as a Rows×1 matrix
     * @param c_idx Column index
     * @return Column vector
     */
    [[nodiscard]] constexpr Matrix<Rows, 1, T> col_vector(size_t c_idx) const {
        Matrix<Rows, 1, T> result;
        for (size_t r = 0; r < Rows; ++r) {
            result.data_[r] = data_[(r * Cols) + c_idx];
        }
        return result;
    }

    /**
     * @brief Extract the first NewRows rows as a thin view (Block)
     * @tparam NewRows Number of rows to extract (must be <= Rows)
     * @return Block view with the first NewRows rows
     */
    template<size_t NewRows>
    [[nodiscard]] constexpr Block<NewRows, Cols, Cols, T> head() {
        static_assert(NewRows <= Rows, "head<NewRows> called with NewRows > Rows");
        return Block<NewRows, Cols, Cols, T>(data_.data(), 0);
    }

    /**
     * @brief Extract the first NewRows rows as a thin const view (Block)
     * @tparam NewRows Number of rows to extract (must be <= Rows)
     * @return Const Block view with the first NewRows rows
     */
    template<size_t NewRows>
    [[nodiscard]] constexpr Block<NewRows, Cols, Cols, const T> head() const {
        static_assert(NewRows <= Rows, "head<NewRows> called with NewRows > Rows");
        return Block<NewRows, Cols, Cols, const T>(data_.data(), 0);
    }

    /**
     * @brief Extract the last NewRows rows as a thin view (Block)
     * @tparam NewRows Number of rows to extract (must be <= Rows)
     * @return Block view with the last NewRows rows
     */
    template<size_t NewRows>
    [[nodiscard]] constexpr Block<NewRows, Cols, Cols, T> tail() {
        static_assert(NewRows <= Rows, "tail<NewRows> called with NewRows > Rows");
        return Block<NewRows, Cols, Cols, T>(&data_[(Rows - NewRows) * Cols]);
    }

    /**
     * @brief Extract the last NewRows rows as a thin const view (Block)
     * @tparam NewRows Number of rows to extract (must be <= Rows)
     * @return Const Block view with the last NewRows rows
     */
    template<size_t NewRows>
    [[nodiscard]] constexpr Block<NewRows, Cols, Cols, const T> tail() const {
        static_assert(NewRows <= Rows, "tail<NewRows> called with NewRows > Rows");
        return Block<NewRows, Cols, Cols, const T>(&data_[(Rows - NewRows) * Cols]);
    }
};

template<typename T, size_t N, size_t M>
struct is_matrix_type<Matrix<N, M, T>> : std::true_type {};

// Common small-matrix aliases
template<typename T>
using Mat2 = Matrix<2, 2, T>;

template<typename T>
using Mat3 = Matrix<3, 3, T>;

template<typename T>
using Mat4 = Matrix<4, 4, T>;

template<typename T>
using Mat3x4 = Matrix<3, 4, T>;

template<typename T>
using Mat4x3 = Matrix<4, 3, T>;

// =============================================================================
// Free operators for MatrixLike types (Matrix, Block, RowView, ColView, etc.)
// =============================================================================

/**
 * @brief Addition of two MatrixLike types (with broadcasting support)
 */
template<MatrixLike A, MatrixLike B>
[[nodiscard]] constexpr auto operator+(const A& a, const B& b) {
    constexpr size_t ResRows = damp::max(A::rows(), B::rows());
    constexpr size_t ResCols = damp::max(A::cols(), B::cols());
    static_assert(
        (A::rows() == 1 || A::rows() == ResRows) && (B::rows() == 1 || B::rows() == ResRows),
        "Incompatible row dimensions for broadcasting"
    );
    static_assert(
        (A::cols() == 1 || A::cols() == ResCols) && (B::cols() == 1 || B::cols() == ResCols),
        "Incompatible column dimensions for broadcasting"
    );
    using T = typename A::value_type;
    Matrix<ResRows, ResCols, T> result;
    for (size_t i = 0; i < ResRows; ++i) {
        for (size_t j = 0; j < ResCols; ++j) {
            size_t ai = (A::rows() == 1) ? 0 : i;
            size_t aj = (A::cols() == 1) ? 0 : j;
            size_t bi = (B::rows() == 1) ? 0 : i;
            size_t bj = (B::cols() == 1) ? 0 : j;
            result(i, j) = static_cast<T>(a(ai, aj)) + static_cast<T>(b(bi, bj));
        }
    }
    if constexpr (std::is_same_v<A, B> && std::is_constructible_v<A, decltype(result)>) {
        return A(result);
    } else {
        return result;
    }
}

/**
 * @brief Subtraction of two MatrixLike types (with broadcasting support)
 */
template<MatrixLike A, MatrixLike B>
[[nodiscard]] constexpr auto operator-(const A& a, const B& b) {
    constexpr size_t ResRows = damp::max(A::rows(), B::rows());
    constexpr size_t ResCols = damp::max(A::cols(), B::cols());
    static_assert(
        (A::rows() == 1 || A::rows() == ResRows) && (B::rows() == 1 || B::rows() == ResRows),
        "Incompatible row dimensions for broadcasting"
    );
    static_assert(
        (A::cols() == 1 || A::cols() == ResCols) && (B::cols() == 1 || B::cols() == ResCols),
        "Incompatible column dimensions for broadcasting"
    );
    using T = typename A::value_type;
    Matrix<ResRows, ResCols, T> result;
    for (size_t i = 0; i < ResRows; ++i) {
        for (size_t j = 0; j < ResCols; ++j) {
            size_t ai = (A::rows() == 1) ? 0 : i;
            size_t aj = (A::cols() == 1) ? 0 : j;
            size_t bi = (B::rows() == 1) ? 0 : i;
            size_t bj = (B::cols() == 1) ? 0 : j;
            result(i, j) = static_cast<T>(a(ai, aj)) - static_cast<T>(b(bi, bj));
        }
    }
    if constexpr (std::is_same_v<A, B> && std::is_constructible_v<A, decltype(result)>) {
        return A(result);
    } else {
        return result;
    }
}

/**
 * @brief Matrix multiplication of two MatrixLike types
 */
template<MatrixLike A, MatrixLike B>
    requires(A::cols() == B::rows())
[[nodiscard]] constexpr auto operator*(const A& a, const B& b) {
    using T = typename A::value_type;
    Matrix<A::rows(), B::cols(), T> result;
    for (size_t i = 0; i < A::rows(); ++i) {
        for (size_t k = 0; k < A::cols(); ++k) {
            T aik = static_cast<T>(a(i, k));
            for (size_t j = 0; j < B::cols(); ++j) {
                result(i, j) += aik * static_cast<T>(b(k, j));
            }
        }
    }
    return result;
}

/**
 * @brief Symmetric congruence (quadratic) form  S = M X Mᵀ
 *
 * Computes the congruence transform of a symmetric matrix X by M. Because X is
 * symmetric the result is symmetric, so only the lower triangle is evaluated and
 * mirrored: the output is symmetric to the last bit (no rounding-induced
 * asymmetry) and the outer product costs ~half a general triple product. This is
 * the standard building block for covariance propagation — A P Aᵀ, C P Cᵀ, the
 * Joseph-form update — where keeping the covariance exactly symmetric is what
 * preserves its positive-definiteness across iterations.
 *
 * Precondition: X is symmetric. If it is not, the result is built from the lower
 * triangle of M X Mᵀ and is therefore not the true product. (For complex T this
 * computes M X Mᵀ with the plain transpose, not the Hermitian adjoint.)
 *
 * @tparam Mat MatrixLike, p×n
 * @tparam Sym MatrixLike, n×n (symmetric)
 * @return p×p symmetric Matrix  M X Mᵀ
 */
template<MatrixLike Mat, MatrixLike Sym>
    requires(Mat::cols() == Sym::rows() && Sym::rows() == Sym::cols())
[[nodiscard]] constexpr auto quadratic_form(const Mat& m, const Sym& x) {
    using T = typename Mat::value_type;
    constexpr size_t P = Mat::rows();
    constexpr size_t N = Mat::cols();
    const auto       mx = m * x; // p×n
    Matrix<P, P, T>  result{};
    for (size_t i = 0; i < P; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            T acc{};
            for (size_t k = 0; k < N; ++k) {
                acc += static_cast<T>(mx(i, k)) * static_cast<T>(m(j, k));
            }
            result(i, j) = acc;
            result(j, i) = acc; // mirror → exact symmetry
        }
    }
    return result;
}

/**
 * @brief Unary negation for any MatrixLike type
 */
template<MatrixLike A>
[[nodiscard]] constexpr auto operator-(const A& a) {
    using T = typename A::value_type;
    Matrix<A::rows(), A::cols(), T> result;
    for (size_t i = 0; i < A::rows(); ++i) {
        for (size_t j = 0; j < A::cols(); ++j) {
            result(i, j) = -static_cast<T>(a(i, j));
        }
    }
    if constexpr (std::is_constructible_v<A, decltype(result)>) {
        return A(result);
    } else {
        return result;
    }
}

// MatrixLike + Scalar
template<MatrixLike A, typename Scalar>
    requires(std::is_arithmetic_v<Scalar> || is_matrix_element_v<Scalar>)
[[nodiscard]] constexpr auto operator+(const A& a, Scalar scalar) {
    using T = typename A::value_type;
    Matrix<A::rows(), A::cols(), T> result;
    for (size_t i = 0; i < A::rows(); ++i) {
        for (size_t j = 0; j < A::cols(); ++j) {
            result(i, j) = static_cast<T>(a(i, j)) + static_cast<T>(scalar);
        }
    }
    if constexpr (std::is_constructible_v<A, decltype(result)>) {
        return A(result);
    } else {
        return result;
    }
}

// MatrixLike - Scalar
template<MatrixLike A, typename Scalar>
    requires(std::is_arithmetic_v<Scalar> || is_matrix_element_v<Scalar>)
[[nodiscard]] constexpr auto operator-(const A& a, Scalar scalar) {
    using T = typename A::value_type;
    Matrix<A::rows(), A::cols(), T> result;
    for (size_t i = 0; i < A::rows(); ++i) {
        for (size_t j = 0; j < A::cols(); ++j) {
            result(i, j) = static_cast<T>(a(i, j)) - static_cast<T>(scalar);
        }
    }
    if constexpr (std::is_constructible_v<A, decltype(result)>) {
        return A(result);
    } else {
        return result;
    }
}

// Scalar - MatrixLike
template<typename Scalar, MatrixLike A>
    requires(std::is_arithmetic_v<Scalar> || is_matrix_element_v<Scalar>)
[[nodiscard]] constexpr auto operator-(Scalar scalar, const A& a) {
    using T = typename A::value_type;
    Matrix<A::rows(), A::cols(), T> result;
    for (size_t i = 0; i < A::rows(); ++i) {
        for (size_t j = 0; j < A::cols(); ++j) {
            result(i, j) = static_cast<T>(scalar) - static_cast<T>(a(i, j));
        }
    }
    if constexpr (std::is_constructible_v<A, decltype(result)>) {
        return A(result);
    } else {
        return result;
    }
}

// MatrixLike * Scalar
template<MatrixLike A, typename Scalar>
    requires(std::is_arithmetic_v<Scalar> || is_matrix_element_v<Scalar>)
[[nodiscard]] constexpr auto operator*(const A& a, Scalar scalar) {
    using T = typename A::value_type;
    Matrix<A::rows(), A::cols(), T> result;
    for (size_t i = 0; i < A::rows(); ++i) {
        for (size_t j = 0; j < A::cols(); ++j) {
            result(i, j) = static_cast<T>(a(i, j)) * static_cast<T>(scalar);
        }
    }
    if constexpr (std::is_constructible_v<A, decltype(result)>) {
        return A(result);
    } else {
        return result;
    }
}

// Scalar * MatrixLike
template<typename Scalar, MatrixLike A>
    requires(std::is_arithmetic_v<Scalar> || is_matrix_element_v<Scalar>)
[[nodiscard]] constexpr auto operator*(Scalar scalar, const A& a) {
    return a * scalar;
}

// MatrixLike / Scalar
template<MatrixLike A, typename Scalar>
    requires(std::is_arithmetic_v<Scalar> || is_matrix_element_v<Scalar>)
[[nodiscard]] constexpr auto operator/(const A& a, Scalar scalar) {
    using T = typename A::value_type;
    Matrix<A::rows(), A::cols(), T> result;
    for (size_t i = 0; i < A::rows(); ++i) {
        for (size_t j = 0; j < A::cols(); ++j) {
            result(i, j) = static_cast<T>(a(i, j)) / static_cast<T>(scalar);
        }
    }
    if constexpr (std::is_constructible_v<A, decltype(result)>) {
        return A(result);
    } else {
        return result;
    }
}

} // namespace damp
