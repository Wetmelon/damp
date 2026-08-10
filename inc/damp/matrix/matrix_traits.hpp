// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file matrix_traits.hpp
 * @brief Matrix element traits and tolerances
 */

#include <concepts>
#include <cstddef>
#include <type_traits>

#include "damp/math/complex.hpp"

namespace damp {

// True if T is damp::complex\<U\> for some U
template<typename T>
struct is_complex : std::false_type {};

template<typename T>
struct is_complex<damp::complex<T>> : std::true_type {};

template<typename T>
inline constexpr bool is_complex_v = is_complex<T>::value;

/**
 * @brief True if @p T may be a Matrix element type
 *
 * Allowed: @c float, @c double, @c damp::complex&lt;float&gt;,
 * @c damp::complex&lt;double&gt;, and const-qualified forms.
 *
 * Rejected on purpose: nested @c Matrix, integers, @c long double, pointers, …
 * Multi-parameter maps use @c Lut* ; block structure uses views on a flat matrix.
 */
template<typename T>
struct is_matrix_element : std::bool_constant<
                               std::is_same_v<std::remove_const_t<T>, float> || std::is_same_v<std::remove_const_t<T>, double> || is_complex<std::remove_const_t<T>>::value> {};

template<typename T>
inline constexpr bool is_matrix_element_v = is_matrix_element<T>::value;

// Extracts the underlying scalar type (float/double) from T or damp::complex\<T\>
template<typename T>
struct scalar_type {
    using type = T;
};

template<typename T>
struct scalar_type<damp::complex<T>> {
    using type = T;
};

template<typename T>
using scalar_type_t = typename scalar_type<T>::type;

/**
 * @brief Concept for any type that provides 2D matrix-like element access
 *
 * A MatrixLike type must provide:
 * - Static `rows()` and `cols()` returning dimensions
 * - `operator()(size_t, size_t)` for element access
 * - A `value_type` typedef
 */
template<typename M>
concept MatrixLike = requires(const M& m, size_t i, size_t j) {
    { M::rows() } -> std::convertible_to<size_t>;
    { M::cols() } -> std::convertible_to<size_t>;
    { m(i, j) } -> std::convertible_to<typename M::value_type>;
};

/**
 * @brief Concept for a MatrixLike type with specific dimensions
 */
template<typename M, size_t R, size_t C>
concept MatrixLikeOf = MatrixLike<M> && (M::rows() == R) && (M::cols() == C);

/// Type-appropriate default tolerance for floating-point comparisons.
/// float  ~7 decimal digits  → 1e-6
/// double ~15 decimal digits → 1e-12
template<typename T>
constexpr scalar_type_t<std::remove_const_t<T>> default_tol() {
    using real_t = scalar_type_t<std::remove_const_t<T>>;
    if constexpr (std::is_same_v<real_t, float>) {
        return real_t{1e-6f};
    } else {
        return real_t{1e-12};
    }
}

} // namespace damp
