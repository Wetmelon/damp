// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file concat.hpp
 * @brief Size-deducing matrix vocabulary: diag, blkdiag, hcat, vcat.
 *
 * These are the assembly words for weights and augmented plants. They live
 * here (reachable from control.hpp) so firmware does not pull workbench.hpp.
 */

#include <cstddef>
#include <type_traits>

#include "colvec.hpp"
#include "core.hpp"
#include "damp/backend.hpp"
#include "matrix_traits.hpp"
#include "rowvec.hpp"

namespace damp {
namespace detail {

template<typename T>
concept diag_scalar = is_matrix_element_v<T> || std::is_integral_v<T>;

template<typename Common>
using diag_scalar_t = std::conditional_t<is_matrix_element_v<Common>, Common, double>;

template<size_t R, size_t C, typename T, MatrixLike M>
constexpr void copy_block(Matrix<R, C, T>& dest, size_t r0, size_t c0, const M& src) {
    for (size_t i = 0; i < M::rows(); ++i) {
        for (size_t j = 0; j < M::cols(); ++j) {
            dest(r0 + i, c0 + j) = static_cast<T>(src(i, j));
        }
    }
}

} // namespace detail

/**
 * @brief Square diagonal matrix from scalars: diag(10.0, 0.1, 100.0, 10.0)
 * @ingroup linear_algebra
 *
 * Integral arguments become double so `diag(10, 1, 100)` is a valid Q.
 *
 * @note Compare with MATLAB®'s diag([a b …]).
 */
template<detail::diag_scalar... Ts>
    requires(sizeof...(Ts) >= 1)
[[nodiscard]] constexpr auto diag(Ts... xs) noexcept {
    using T = detail::diag_scalar_t<std::common_type_t<Ts...>>;
    constexpr size_t N = sizeof...(Ts);
    Matrix<N, N, T>  result{};
    const T          vals[N] = {static_cast<T>(xs)...};
    for (size_t i = 0; i < N; ++i) {
        result(i, i) = vals[i];
    }
    return result;
}

/**
 * @brief Square diagonal matrix from an array of diagonal entries
 * @ingroup linear_algebra
 *
 * @note Compare with MATLAB®'s diag(v).
 */
template<size_t N, typename T>
    requires is_matrix_element_v<T>
[[nodiscard]] constexpr Matrix<N, N, T> diag(const damp::array<T, N>& d) noexcept {
    return Matrix<N, N, T>::diagonal(d);
}

/**
 * @brief Square diagonal matrix from a column vector (N ≠ 1; 1×1 uses extract)
 * @ingroup linear_algebra
 */
template<MatrixLike M>
    requires(M::cols() == 1 && M::rows() != 1)
[[nodiscard]] constexpr auto diag(const M& v) noexcept {
    using T = M::value_type;
    constexpr size_t N = M::rows();
    Matrix<N, N, T>  result{};
    for (size_t i = 0; i < N; ++i) {
        result(i, i) = v(i, 0);
    }
    return result;
}

/**
 * @brief Square diagonal matrix from a row vector (N ≠ 1)
 * @ingroup linear_algebra
 */
template<MatrixLike M>
    requires(M::rows() == 1 && M::cols() != 1)
[[nodiscard]] constexpr auto diag(const M& v) noexcept {
    using T = M::value_type;
    constexpr size_t N = M::cols();
    Matrix<N, N, T>  result{};
    for (size_t i = 0; i < N; ++i) {
        result(i, i) = v(0, i);
    }
    return result;
}

/**
 * @brief Main diagonal of a square matrix as a column vector
 * @ingroup linear_algebra
 *
 * @note Compare with MATLAB®'s diag(A).
 */
template<MatrixLike M>
    requires(M::rows() == M::cols())
[[nodiscard]] constexpr auto diag(const M& A) noexcept {
    constexpr size_t                  N = M::rows();
    ColVec<N, typename M::value_type> result{};
    for (size_t i = 0; i < N; ++i) {
        result[i] = A(i, i);
    }
    return result;
}

/**
 * @brief Block-diagonal assembly of MatrixLike arguments
 * @ingroup linear_algebra
 *
 * @note Compare with MATLAB®'s blkdiag.
 */
template<MatrixLike... Ms>
    requires(sizeof...(Ms) >= 1)
[[nodiscard]] constexpr auto blkdiag(const Ms&... ms) noexcept {
    constexpr size_t R = (Ms::rows() + ...);
    constexpr size_t C = (Ms::cols() + ...);
    using T = std::common_type_t<typename Ms::value_type...>;
    Matrix<R, C, T> result{};
    size_t          r0 = 0;
    size_t          c0 = 0;
    auto            put = [&](const auto& m) {
        detail::copy_block(result, r0, c0, m);
        r0 += m.rows();
        c0 += m.cols();
    };
    (put(ms), ...);
    return result;
}

/**
 * @brief Horizontal concatenation (same row count)
 * @ingroup linear_algebra
 *
 * @note Compare with MATLAB®'s [A B] / horzcat.
 */
template<MatrixLike M0, MatrixLike... Ms>
    requires((Ms::rows() == M0::rows()) && ...)
[[nodiscard]] constexpr auto hcat(const M0& m0, const Ms&... ms) noexcept {
    constexpr size_t R = M0::rows();
    constexpr size_t C = (M0::cols() + ... + Ms::cols());
    using T = std::common_type_t<typename M0::value_type, typename Ms::value_type...>;
    Matrix<R, C, T> result{};
    size_t          c0 = 0;
    auto            put = [&](const auto& m) {
        detail::copy_block(result, 0, c0, m);
        c0 += m.cols();
    };
    put(m0);
    (put(ms), ...);
    return result;
}

/**
 * @brief Vertical concatenation (same column count)
 * @ingroup linear_algebra
 *
 * @note Compare with MATLAB®'s [A; B] / vertcat.
 */
template<MatrixLike M0, MatrixLike... Ms>
    requires((Ms::cols() == M0::cols()) && ...)
[[nodiscard]] constexpr auto vcat(const M0& m0, const Ms&... ms) noexcept {
    constexpr size_t R = (M0::rows() + ... + Ms::rows());
    constexpr size_t C = M0::cols();
    using T = std::common_type_t<typename M0::value_type, typename Ms::value_type...>;
    Matrix<R, C, T> result{};
    size_t          r0 = 0;
    auto            put = [&](const auto& m) {
        detail::copy_block(result, r0, 0, m);
        r0 += m.rows();
    };
    put(m0);
    (put(ms), ...);
    return result;
}

} // namespace damp
