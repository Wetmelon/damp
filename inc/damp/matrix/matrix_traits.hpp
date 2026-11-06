// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file matrix_traits.hpp
 * @brief Matrix element traits and tolerances
 */

#include "damp/backend.hpp"
#include "damp/math/complex.hpp"

namespace damp {

// True if T is damp::complex\<U\> for some U
template<typename T>
struct is_complex : damp::false_type {};

template<typename T>
struct is_complex<damp::complex<T>> : damp::true_type {};

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
struct is_matrix_element : damp::bool_constant<
                               damp::is_same_v<damp::remove_const_t<T>, float> || damp::is_same_v<damp::remove_const_t<T>, double> || is_complex<damp::remove_const_t<T>>::value> {};

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

namespace detail {
template<typename From, typename To>
inline constexpr bool converts_exactly_v = damp::is_same_v<damp::remove_cv_t<From>, damp::remove_cv_t<To>>
                                        || (is_complex_v<damp::remove_cv_t<To>> && damp::is_same_v<damp::remove_cv_t<From>, scalar_type_t<damp::remove_cv_t<To>>>);

/** @brief Element type of a mixed operation: the shared type, or the complex one when the other side is its real. */
template<typename A, typename B>
using common_element_t = damp::conditional_t<is_complex_v<damp::remove_cv_t<B>>, damp::remove_cv_t<B>, damp::remove_cv_t<A>>;
} // namespace detail

/**
 * @brief Matrix elements of type @p From convert to @p To without rounding: the same type, or a real
 * into the complex of its own precision
 *
 * Implicit matrix conversions, assignments and matrix–matrix operators require it. A change of
 * precision (double ↔ float) is written `.as<U>()` or with the explicit converting constructor.
 */
template<typename From, typename To>
concept ConvertsExactly = detail::converts_exactly_v<From, To>; // change precision with .as<U>()

/**
 * @brief Two matrix element types mix in one operator without rounding: equal, or a real and the
 * complex of its own precision
 */
template<typename A, typename B>
concept SamePrecision = ConvertsExactly<A, B> || ConvertsExactly<B, A>; // change precision with .as<U>()

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
    { M::rows() } -> damp::convertible_to<size_t>;
    { M::cols() } -> damp::convertible_to<size_t>;
    { m(i, j) } -> damp::convertible_to<typename M::value_type>;
};

/**
 * @brief Concept for a MatrixLike type with specific dimensions
 */
template<typename M, size_t R, size_t C>
concept MatrixLikeOf = MatrixLike<M> && (M::rows() == R) && (M::cols() == C);

/**
 * @brief Machine epsilon of T's real scalar (MATLAB® eps). Scale it by problem size and
 * magnitude to build a relative tolerance.
 */
template<typename T>
constexpr scalar_type_t<damp::remove_const_t<T>> eps() {
    return damp::numeric_limits<scalar_type_t<damp::remove_const_t<T>>>::epsilon();
}

/**
 * @brief Type-appropriate default tolerance for floating-point comparisons.
 * float  ~7 decimal digits  → 1e-6
 * double ~15 decimal digits → 1e-12
 */
template<typename T>
constexpr scalar_type_t<damp::remove_const_t<T>> default_tol() {
    using real_t = scalar_type_t<damp::remove_const_t<T>>;
    if constexpr (damp::is_same_v<real_t, float>) {
        return real_t{1e-6f};
    } else {
        return real_t{1e-12};
    }
}

} // namespace damp
