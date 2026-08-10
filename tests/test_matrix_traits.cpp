// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <type_traits>

#include "damp/math/complex.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

/**
 * @file test_matrix_traits.cpp
 * @brief Pins the matrix type traits and concepts: is_complex_v,
 *        is_matrix_element_v, scalar_type_t, default_tol, and the MatrixLike /
 *        MatrixLikeOf concepts (both acceptance and rejection). These are
 *        compile-time contracts exercised transitively everywhere; this file
 *        asserts them directly so a regression fails loudly and locally.
 */

// -------- is_complex_v --------
static_assert(is_complex_v<damp::complex<float>>);
static_assert(is_complex_v<damp::complex<double>>);
static_assert(!is_complex_v<float>);
static_assert(!is_complex_v<double>);
static_assert(!is_complex_v<int>);

// -------- is_matrix_element_v (accept) --------
static_assert(is_matrix_element_v<float>);
static_assert(is_matrix_element_v<double>);
static_assert(is_matrix_element_v<damp::complex<float>>);
static_assert(is_matrix_element_v<damp::complex<double>>);
static_assert(is_matrix_element_v<const double>); // const-qualified allowed

// -------- is_matrix_element_v (reject) --------
static_assert(!is_matrix_element_v<int>);
static_assert(!is_matrix_element_v<long double>);
static_assert(!is_matrix_element_v<bool>);
static_assert(!is_matrix_element_v<float*>);
// Nested Matrix is not a scalar element (would be a block-matrix type, not Lut n-D).
static_assert(!is_matrix_element_v<Matrix<2, 2, double>>);
static_assert(!is_matrix_element_v<ColVec<3, float>>);

// -------- scalar_type_t --------
static_assert(std::is_same_v<scalar_type_t<float>, float>);
static_assert(std::is_same_v<scalar_type_t<double>, double>);
static_assert(std::is_same_v<scalar_type_t<damp::complex<float>>, float>);
static_assert(std::is_same_v<scalar_type_t<damp::complex<double>>, double>);

// -------- MatrixLike concept (acceptance) --------
static_assert(MatrixLike<Matrix<2, 3, double>>);
static_assert(MatrixLike<Matrix<4, 4, damp::complex<double>>>);
static_assert(MatrixLike<ColVec<3, double>>);
static_assert(MatrixLike<RowVec<3, double>>);

// -------- MatrixLike concept (rejection) --------
static_assert(!MatrixLike<double>);
static_assert(!MatrixLike<int>);
struct NotAMatrix {
    int x;
};
static_assert(!MatrixLike<NotAMatrix>);

// -------- MatrixLikeOf (dimension-constrained) --------
static_assert(MatrixLikeOf<Matrix<2, 3, double>, 2, 3>);
static_assert(!MatrixLikeOf<Matrix<2, 3, double>, 3, 2>); // wrong dims
static_assert(!MatrixLikeOf<double, 1, 1>);

TEST_SUITE("matrix_traits") {
    TEST_CASE("default_tol is type-appropriate") {
        CHECK(default_tol<float>() == doctest::Approx(1e-6f));
        CHECK(default_tol<double>() == doctest::Approx(1e-12));
        CHECK(default_tol<damp::complex<float>>() == doctest::Approx(1e-6f));
        CHECK(default_tol<damp::complex<double>>() == doctest::Approx(1e-12));
        CHECK(default_tol<const double>() == doctest::Approx(1e-12));

        // Return type is the real scalar type, not the complex element type.
        static_assert(std::is_same_v<decltype(default_tol<damp::complex<double>>()), double>);
        static_assert(std::is_same_v<decltype(default_tol<float>()), float>);
    }

    TEST_CASE("is_matrix_type trait identifies Matrix specializations") {
        static_assert(is_matrix_type<Matrix<3, 3, double>>::value);
        static_assert(!is_matrix_type<double>::value);
        static_assert(!is_matrix_type<int>::value);
        CHECK(true); // anchor a runtime assertion for the case
    }
}
