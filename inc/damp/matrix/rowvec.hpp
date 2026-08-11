// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file rowvec.hpp
 * @brief Row vector type
 */

#include <cstddef>
#include <initializer_list>

#include "core.hpp"
#include "damp/backend.hpp"
#include "matrix_traits.hpp"

namespace damp {

/**
 * @brief Row vector specialization of Matrix<1, N, T>
 * @ingroup linear_algebra
 * @tparam N Vector dimension
 * @tparam T Element type
 */
template<size_t N, typename T = double>
struct RowVec : public Matrix<1, N, T> {
    /**
     * @brief Default constructor
     */
    constexpr RowVec() = default;
    /**
     * @brief Copy constructor and assignment
     */
    constexpr RowVec(const RowVec&) = default;
    constexpr RowVec& operator=(const RowVec&) = default;

    /**
     * @brief Move constructor and assignment
     */
    constexpr RowVec(RowVec&&) = default;
    constexpr RowVec& operator=(RowVec&&) = default;
    /**
     * @brief Destructor
     */
    constexpr ~RowVec() = default;

    /**
     * @brief Constructor from initializer list
     */
    constexpr RowVec(std::initializer_list<T> values) : Matrix<1, N, T>() {
        size_t i = 0;
        for (const auto& val : values) {
            if (i < N) {
                this->data_[i] = val;
            }
            ++i;
        }
    }

    // Same-precision wrap of a 1×N matrix is a no-op (identical shape and
    // storage), so it is implicit. Cross-precision conversions stay explicit
    // below, matching the library's "no silent float↔double" rule (use .as\<U\>()).
    constexpr RowVec(const Matrix<1, N, T>& other) : Matrix<1, N, T>(other) {}

    template<typename U>
    explicit constexpr RowVec(const RowVec<N, U>& other) : Matrix<1, N, T>(other) {}

    template<typename U>
    explicit constexpr RowVec(const Matrix<1, N, U>& other) : Matrix<1, N, T>(other) {}

    template<typename U>
    constexpr RowVec& operator=(const Matrix<1, N, U>& other) {
        Matrix<1, N, T>::operator=(other);
        return *this;
    }

    // Bring base class operator() into scope (single-arg overload hides it otherwise)
    using Matrix<1, N, T>::operator();

    /**
     * @brief Element access operator
     * @param idx Vector index
     * @return Const reference to element
     */
    constexpr const T& operator[](size_t idx) const { return this->data_[idx]; }
    constexpr const T& operator()(size_t idx) const { return this->data_[idx]; }

    /**
     * @brief Element access operator
     * @param idx Vector index
     * @return Reference to element
     */
    constexpr T& operator[](size_t idx) { return this->data_[idx]; }
    constexpr T& operator()(size_t idx) { return this->data_[idx]; }

    /**
     * @brief Inner product of two row vectors
     *
     * For real vectors:    dot(u,v) = sum( u_i * v_i )
     * For complex vectors: dot(u,v) = sum( conj(u_i) * v_i )
     *
     * @param vec1 First vector (conjugated)
     * @param vec2 Second vector
     * @return Scalar inner product
     */
    [[nodiscard]] friend constexpr T dot(const RowVec<N, T>& vec1, const RowVec<N, T>& vec2) {
        T result = T{0};
        for (std::size_t i = 0; i < N; ++i) {
            result += damp::conj(vec1.data_[i]) * vec2.data_[i];
        }
        return result;
    }

    // Cross product for 3D vectors
    template<size_t D = N>
    /**
     * @brief 3D cross product with another 3D row vector
     * @param other The other 3D row vector
     * @return 3D cross product result
     */
    [[nodiscard]] constexpr RowVec<3, T> cross(const RowVec<3, T>& other) const
        requires(D == 3)
    {
        return RowVec<3, T>{
            this->data_[1] * other.data_[2] - this->data_[2] * other.data_[1],
            this->data_[2] * other.data_[0] - this->data_[0] * other.data_[2],
            this->data_[0] * other.data_[1] - this->data_[1] * other.data_[0]
        };
    }

    /**
     * @brief Euclidean norm (magnitude) of the vector
     *
     * ||v|| = sqrt( sum |v_i|^2 )
     *
     * Always returns a real value, even for complex-valued vectors.
     *
     * @return Vector magnitude (real scalar)
     */
    [[nodiscard]] constexpr scalar_type_t<T> norm() const {
        using real_t = scalar_type_t<T>;
        real_t sum_sq = real_t{0};
        for (size_t i = 0; i < N; ++i) {
            auto abs_val = damp::abs(this->data_[i]);
            sum_sq += abs_val * abs_val;
        }
        return damp::sqrt(sum_sq);
    }

    /**
     * @brief Normalized vector (unit length)
     * @return Normalized vector
     */
    [[nodiscard]] constexpr RowVec normalized() const {
        auto n = norm();
        if (n == scalar_type_t<T>{0}) {
            return *this;
        }
        return *this * (T{1} / static_cast<T>(n));
    }

    /**
     * @brief Contiguous subvector view (Eigen-style segment)
     *
     * Non-owning @c Block over @p Len elements starting at @p start. Equivalent
     * to @c block\<1, Len\>(0, start). Precondition: @c start + Len <= N.
     *
     * @tparam Len Segment length (compile-time)
     * @param start First index of the segment
     * @return Mutable block view of shape 1 × Len
     */
    template<size_t Len>
    [[nodiscard]] constexpr Block<1, Len, N, T> segment(size_t start) {
        static_assert(Len <= N, "segment length exceeds vector dimension");
        return Block<1, Len, N, T>(this->data_.data(), start);
    }

    /**
     * @brief Contiguous subvector view (const; Eigen-style segment)
     *
     * @tparam Len Segment length (compile-time)
     * @param start First index of the segment
     * @return Const block view of shape 1 × Len
     */
    template<size_t Len>
    [[nodiscard]] constexpr Block<1, Len, N, const T> segment(size_t start) const {
        static_assert(Len <= N, "segment length exceeds vector dimension");
        return Block<1, Len, N, const T>(this->data_.data(), start);
    }

    // Type-preserving compound assignment (base returns Matrix&)
    template<MatrixLike M>
    constexpr RowVec& operator+=(const M& other) {
        Matrix<1, N, T>::operator+=(other);
        return *this;
    }

    template<MatrixLike M>
    constexpr RowVec& operator-=(const M& other) {
        Matrix<1, N, T>::operator-=(other);
        return *this;
    }

    constexpr RowVec& operator+=(T scalar) {
        Matrix<1, N, T>::operator+=(scalar);
        return *this;
    }

    constexpr RowVec& operator-=(T scalar) {
        Matrix<1, N, T>::operator-=(scalar);
        return *this;
    }

    constexpr RowVec& operator*=(T scalar) {
        Matrix<1, N, T>::operator*=(scalar);
        return *this;
    }

    constexpr RowVec& operator/=(T scalar) {
        Matrix<1, N, T>::operator/=(scalar);
        return *this;
    }
};

// Type-preserving operators: RowVec op RowVec → RowVec
template<size_t N, typename T>
[[nodiscard]] constexpr RowVec<N, T> operator+(const RowVec<N, T>& a, const RowVec<N, T>& b) {
    RowVec<N, T> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

template<size_t N, typename T>
[[nodiscard]] constexpr RowVec<N, T> operator-(const RowVec<N, T>& a, const RowVec<N, T>& b) {
    RowVec<N, T> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

template<size_t N, typename T>
[[nodiscard]] constexpr RowVec<N, T> operator-(const RowVec<N, T>& a) {
    RowVec<N, T> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = -a[i];
    }
    return result;
}

template<size_t N, typename T, typename Scalar>
    requires(std::is_arithmetic_v<Scalar> || is_matrix_element_v<Scalar>)
[[nodiscard]] constexpr RowVec<N, T> operator*(const RowVec<N, T>& a, Scalar s) {
    RowVec<N, T> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = a[i] * static_cast<T>(s);
    }
    return result;
}

template<typename Scalar, size_t N, typename T>
    requires(std::is_arithmetic_v<Scalar> || is_matrix_element_v<Scalar>)
[[nodiscard]] constexpr RowVec<N, T> operator*(Scalar s, const RowVec<N, T>& a) {
    return a * s;
}

template<size_t N, typename T, typename Scalar>
    requires(std::is_arithmetic_v<Scalar> || is_matrix_element_v<Scalar>)
[[nodiscard]] constexpr RowVec<N, T> operator/(const RowVec<N, T>& a, Scalar s) {
    RowVec<N, T> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = a[i] / static_cast<T>(s);
    }
    return result;
}

// RowVec * Matrix returns RowVec
template<size_t N, size_t Cols, typename T, typename U>
[[nodiscard]] constexpr RowVec<Cols, T> operator*(const RowVec<N, T>& vec, const Matrix<N, Cols, U>& mat) {
    RowVec<Cols, T> result;
    for (size_t j = 0; j < Cols; ++j) {
        T sum = 0;
        for (size_t k = 0; k < N; ++k) {
            sum += vec[k] * static_cast<T>(mat(k, j));
        }
        result[j] = sum;
    }
    return result;
}

// Deduce RowVec<N, T> from variadic constructor arguments, e.g. RowVec vec{1.0f, 2.0f, 3.0f};
// Integer literals deduce to double
template<typename T, typename... Args>
    requires(std::is_integral_v<T> && (std::is_integral_v<Args> && ...))
RowVec(T, Args...) -> RowVec<1 + sizeof...(Args), double>;

template<typename T, typename... Args>
    requires(!std::is_integral_v<T> || !(std::is_integral_v<Args> && ...))
RowVec(T, Args...) -> RowVec<1 + sizeof...(Args), T>;

template<typename U, size_t M>
RowVec(Matrix<1, M, U>) -> RowVec<M, U>;

template<typename T, size_t N>
RowVec(const damp::array<T, N>&) -> RowVec<N, T>;

template<typename T>
using RowVec2 = RowVec<2, T>;

template<typename T>
using RowVec3 = RowVec<3, T>;

template<typename T>
using RowVec4 = RowVec<4, T>;

} // namespace damp