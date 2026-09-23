// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file block.hpp
 * @brief Non-owning matrix block views
 */

#include <cstddef>
#include <type_traits>

#include "core.hpp"
#include "matrix_traits.hpp"

namespace damp {

// Forward declarations for to_vector()
template<size_t N, typename T>
struct ColVec;
template<size_t N, typename T>
struct RowVec;

/**
 * @brief  Block view (non-owning) into a parent matrix
 * @ingroup linear_algebra
 *
 * @tparam Rows        Number of rows in the block
 * @tparam Cols        Number of columns in the block
 * @tparam ParentCols  Number of columns in the parent matrix
 * @tparam T           Element type
 *
 */
template<size_t Rows, size_t Cols, size_t ParentCols, typename T>
struct Block {
private:
    T* const     data_ptr;
    size_t const offset;

public:
    typedef std::remove_const_t<T> value_type;
    // Create a block from a pointer to the parent data and an offset
    constexpr explicit Block(T* data_ptr_, size_t offset_ = 0)
        : data_ptr(data_ptr_), offset(offset_) {}

    constexpr Block(const Block& other) : data_ptr(other.data_ptr), offset(other.offset) {}
    constexpr Block& operator=(const Block& other) {
        if constexpr (Rows == 0 || Cols == 0) {
            return *this;
        }
        // Direct write when the rectangles do not share an element.
        // Overlap (A.block(0,0) = A.block(1,0)) reads through a snapshot.
        if (detail::rects_share_element(
                storage_first(), Rows, Cols, ParentCols, other.storage_first(), Rows, Cols, ParentCols
            )) {
            Matrix<Rows, Cols, value_type> tmp{};
            for (size_t i = 0; i < Rows; ++i) {
                for (size_t j = 0; j < Cols; ++j) {
                    tmp(i, j) = other(i, j);
                }
            }
            for (size_t i = 0; i < Rows; ++i) {
                for (size_t j = 0; j < Cols; ++j) {
                    (*this)(i, j) = tmp(i, j);
                }
            }
        } else {
            for (size_t i = 0; i < Rows; ++i) {
                for (size_t j = 0; j < Cols; ++j) {
                    (*this)(i, j) = other(i, j);
                }
            }
        }
        return *this;
    }

    constexpr Block(Block&& other) noexcept : data_ptr(other.data_ptr), offset(other.offset) {}
    constexpr Block& operator=(Block&& other) noexcept {
        for (size_t i = 0; i < Rows; ++i) {
            for (size_t j = 0; j < Cols; ++j) {
                (*this)(i, j) = other(i, j);
            }
        }
        return *this;
    }

    [[nodiscard]] constexpr T& operator()(size_t r, size_t c) {
        return data_ptr[offset + (r * ParentCols) + c];
    }

    [[nodiscard]] constexpr const T& operator()(size_t r, size_t c) const {
        return data_ptr[offset + (r * ParentCols) + c];
    }

    constexpr Block& operator=(const Matrix<Rows, Cols, T>& mat) {
        for (size_t i = 0; i < Rows; ++i) {
            for (size_t j = 0; j < Cols; ++j) {
                (*this)(i, j) = mat(i, j);
            }
        }
        return *this;
    }

    /**
     * @brief Assignment from any MatrixLike type with matching dimensions
     */
    template<MatrixLike M>
        requires(M::rows() == Rows && M::cols() == Cols)
    constexpr Block& operator=(const M& other) {
        if constexpr (Rows == 0 || Cols == 0) {
            return *this;
        } else if (detail::storage_overlaps(storage_first(), storage_last(), other)) {
            Matrix<Rows, Cols, value_type> tmp{};
            for (size_t i = 0; i < Rows; ++i) {
                for (size_t j = 0; j < Cols; ++j) {
                    tmp(i, j) = static_cast<value_type>(other(i, j));
                }
            }
            return *this = tmp;
        } else {
            for (size_t i = 0; i < Rows; ++i) {
                for (size_t j = 0; j < Cols; ++j) {
                    (*this)(i, j) = static_cast<value_type>(other(i, j));
                }
            }
            return *this;
        }
    }

    constexpr Block& operator=(T scalar) {
        for (size_t i = 0; i < Rows; ++i) {
            for (size_t j = 0; j < Cols; ++j) {
                (*this)(i, j) = scalar;
            }
        }
        return *this;
    }

    [[nodiscard]] constexpr T*       data() { return data_ptr + offset; }
    [[nodiscard]] constexpr const T* data() const { return data_ptr + offset; }

    /// First element of this rectangle.
    [[nodiscard]] constexpr const value_type* storage_first() const { return data_ptr + offset; }

    /// Last element of this rectangle (inclusive). The span includes gaps when
    /// Cols < ParentCols.
    [[nodiscard]] constexpr const value_type* storage_last() const {
        if constexpr (Rows == 0 || Cols == 0) {
            return data_ptr + offset;
        }
        return data_ptr + offset + ((Rows - 1) * ParentCols) + (Cols - 1);
    }
    [[nodiscard]] static constexpr size_t size() { return Rows * Cols; }
    [[nodiscard]] static constexpr size_t rows() { return Rows; }
    [[nodiscard]] static constexpr size_t cols() { return Cols; }

    // Sub-block access
    template<size_t Brows, size_t Bcols>
    constexpr Block<Brows, Bcols, ParentCols, T> block(size_t start_row, size_t start_col) {
        return Block<Brows, Bcols, ParentCols, T>(data_ptr, offset + (start_row * ParentCols) + start_col);
    }

    template<size_t Brows, size_t Bcols>
    constexpr Block<Brows, Bcols, ParentCols, const T> block(size_t start_row, size_t start_col) const {
        return Block<Brows, Bcols, ParentCols, const T>(
            const_cast<T*>(data_ptr), offset + (start_row * ParentCols) + start_col
        );
    }

    /**
     * @brief Convert to owning Matrix
     */
    [[nodiscard]] constexpr Matrix<Rows, Cols, value_type> to_matrix() const {
        Matrix<Rows, Cols, value_type> result;
        for (size_t i = 0; i < Rows; ++i) {
            for (size_t j = 0; j < Cols; ++j) {
                result(i, j) = (*this)(i, j);
            }
        }
        return result;
    }

    /**
     * @brief Compound addition assignment from any MatrixLike
     */
    template<MatrixLike M>
        requires(M::rows() == Rows && M::cols() == Cols)
    constexpr Block& operator+=(const M& other) {
        for (size_t i = 0; i < Rows; ++i) {
            for (size_t j = 0; j < Cols; ++j) {
                (*this)(i, j) += static_cast<T>(other(i, j));
            }
        }
        return *this;
    }

    /**
     * @brief Compound subtraction assignment from any MatrixLike
     */
    template<MatrixLike M>
        requires(M::rows() == Rows && M::cols() == Cols)
    constexpr Block& operator-=(const M& other) {
        for (size_t i = 0; i < Rows; ++i) {
            for (size_t j = 0; j < Cols; ++j) {
                (*this)(i, j) -= static_cast<T>(other(i, j));
            }
        }
        return *this;
    }

    constexpr Block& operator+=(T scalar) {
        for (size_t i = 0; i < Rows; ++i) {
            for (size_t j = 0; j < Cols; ++j) {
                (*this)(i, j) += scalar;
            }
        }
        return *this;
    }

    constexpr Block& operator-=(T scalar) {
        for (size_t i = 0; i < Rows; ++i) {
            for (size_t j = 0; j < Cols; ++j) {
                (*this)(i, j) -= scalar;
            }
        }
        return *this;
    }

    constexpr Block& operator*=(T scalar) {
        for (size_t i = 0; i < Rows; ++i) {
            for (size_t j = 0; j < Cols; ++j) {
                (*this)(i, j) *= scalar;
            }
        }
        return *this;
    }

    constexpr Block& operator/=(T scalar) {
        for (size_t i = 0; i < Rows; ++i) {
            for (size_t j = 0; j < Cols; ++j) {
                (*this)(i, j) /= scalar;
            }
        }
        return *this;
    }

    /**
     * @brief Convert column block to ColVec
     */
    template<size_t R = Rows, size_t C = Cols>
        requires(C == 1)
    [[nodiscard]] constexpr ColVec<R, value_type> to_vector() const {
        ColVec<R, value_type> result;
        for (size_t i = 0; i < R; ++i) {
            result[i] = (*this)(i, 0);
        }
        return result;
    }

    /**
     * @brief Convert row block to RowVec
     */
    template<size_t R = Rows, size_t C = Cols>
        requires(R == 1 && C > 1)
    [[nodiscard]] constexpr RowVec<C, value_type> to_vector() const {
        RowVec<C, value_type> result;
        for (size_t i = 0; i < C; ++i) {
            result[i] = (*this)(0, i);
        }
        return result;
    }
};

} // namespace damp