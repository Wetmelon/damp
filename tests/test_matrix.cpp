// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>
#include <random>
#include <vector>

#include "damp/matrix/matrix.hpp"
#include "fmt/core.h"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

TEST_CASE("Matrix basic construction and access") {
    Matrix<2, 3, float> mat;
    // Default initialized to 0
    CHECK(mat(0, 0) == 0.0f);
    CHECK(mat(1, 1) == 0.0f);
    CHECK(mat(1, 2) == 0.0f);

    // Set values
    mat(0, 0) = 1.0f;
    mat(0, 1) = 2.0f;
    mat(0, 2) = 3.0f;
    mat(1, 0) = 4.0f;
    mat(1, 1) = 5.0f;
    mat(1, 2) = 6.0f;

    CHECK(mat(0, 0) == 1.0f);
    CHECK(mat(0, 1) == 2.0f);
    CHECK(mat(0, 2) == 3.0f);
    CHECK(mat(1, 0) == 4.0f);
    CHECK(mat(1, 1) == 5.0f);
    CHECK(mat(1, 2) == 6.0f);
}

TEST_CASE("Matrix initializer list constructor") {
    Matrix<2, 2> mat = {{1, 2}, {3, 4}};

    CHECK(mat(0, 0) == 1);
    CHECK(mat(0, 1) == 2);
    CHECK(mat(1, 0) == 3);
    CHECK(mat(1, 1) == 4);
}

TEST_CASE("Matrix copy and assignment") {
    Matrix<2, 2, float> mat1 = {{1.0f, 2.0f}, {3.0f, 4.0f}};

    Matrix mat2 = mat1;
    CHECK(mat2(0, 0) == 1.0f);
    CHECK(mat2(1, 1) == 4.0f);

    Matrix<2, 2, float> mat3;
    mat3 = mat1;
    CHECK(mat3(0, 1) == 2.0f);
}

TEST_CASE("Matrix type conversion") {
    Matrix<2, 2> mat_int = {{1, 2}, {3, 4}};

    Matrix<2, 2, float> mat_float(mat_int);
    CHECK(mat_float(0, 0) == 1.0f);
    CHECK(mat_float(1, 1) == 4.0f);

    Matrix<2, 2, float> mat_assign;
    mat_assign = mat_int;
    CHECK(mat_assign(0, 1) == 2.0f);
}

TEST_CASE("Matrix addition and subtraction") {
    Matrix<2, 2> mat1 = {{1, 2}, {3, 4}};

    Matrix<2, 2> mat2 = {{5, 6}, {7, 8}};

    auto sum = mat1 + mat2;
    CHECK(sum(0, 0) == 6);
    CHECK(sum(1, 1) == 12);

    auto diff = mat2 - mat1;
    CHECK(diff(0, 0) == 4);
    CHECK(diff(1, 1) == 4);

    mat1 += mat2;
    CHECK(mat1(0, 0) == 6);

    mat1 -= mat2;
    CHECK(mat1(0, 0) == 1);
}

TEST_CASE("Matrix scalar operations") {
    Matrix<2, 2> mat = {{1, 2}, {3, 4}};

    Matrix<2, 2> scaled = mat * 2;
    CHECK(scaled(0, 0) == 2);
    CHECK(scaled(1, 1) == 8);

    scaled = 3 * mat;
    CHECK(scaled(0, 1) == 6);

    scaled = mat / 2;
    CHECK(scaled(0, 0) == doctest::Approx(0.5));
    CHECK(scaled(1, 0) == doctest::Approx(1.5));

    mat *= 2;
    CHECK(mat(0, 0) == 2);

    mat /= 2;
    CHECK(mat(0, 0) == 1);
}

TEST_CASE("Matrix equality") {
    Matrix<2, 2> mat1 = {{1, 2}, {3, 4}};

    Matrix<2, 2> mat2 = {{1, 2}, {3, 4}};

    Matrix<2, 2> mat3 = {{1, 2}, {3, 5}};

    CHECK(mat1 == mat2);
    CHECK(mat1 != mat3);
}

TEST_CASE("Matrix unary negation") {
    Matrix<2, 2> mat = {{1, -2}, {3, 4}};

    auto neg = -mat;
    CHECK(neg(0, 0) == -1);
    CHECK(neg(0, 1) == 2);
    CHECK(neg(1, 0) == -3);
    CHECK(neg(1, 1) == -4);
}

TEST_CASE("Matrix transpose") {
    Matrix<2, 3> mat = {{1, 2, 3}, {4, 5, 6}};

    auto trans = mat.transpose();
    CHECK(trans(0, 0) == 1);
    CHECK(trans(0, 1) == 4);
    CHECK(trans(1, 0) == 2);
    CHECK(trans(1, 1) == 5);
    CHECK(trans(2, 0) == 3);
    CHECK(trans(2, 1) == 6);
}

TEST_CASE("Matrix multiplication") {
    Matrix<2, 2> mat1 = {{1, 2}, {3, 4}};

    Matrix<2, 2> mat2 = {{5, 6}, {7, 8}};

    auto prod = mat1 * mat2;
    CHECK(prod(0, 0) == 19); // 1*5 + 2*7
    CHECK(prod(0, 1) == 22); // 1*6 + 2*8
    CHECK(prod(1, 0) == 43); // 3*5 + 4*7
    CHECK(prod(1, 1) == 50); // 3*6 + 4*8
}

TEST_CASE("Matrix constexpr") {
    constexpr Matrix<2, 2> mat = {{1, 2}, {3, 4}};

    CHECK(mat(0, 0) == 1);
    CHECK(mat(1, 1) == 4);
}

TEST_CASE("Matrix vector operations") {
    // Using Matrix<N, 1, T> as column vectors
    Matrix<3, 1> vec1 = {{1}, {2}, {3}};

    Matrix<3, 1> vec2 = {{4}, {5}, {6}};

    // Vector addition
    auto vec_sum = vec1 + vec2;
    CHECK(vec_sum(0, 0) == 5);
    CHECK(vec_sum(1, 0) == 7);
    CHECK(vec_sum(2, 0) == 9);

    // Scalar multiplication
    auto vec_scaled = vec1 * 2;
    CHECK(vec_scaled(0, 0) == 2);
    CHECK(vec_scaled(1, 0) == 4);
    CHECK(vec_scaled(2, 0) == 6);

    // Vector subtraction
    auto vec_diff = vec2 - vec1;
    CHECK(vec_diff(0, 0) == 3);
    CHECK(vec_diff(1, 0) == 3);
    CHECK(vec_diff(2, 0) == 3);

    // Unary negation
    auto vec_neg = -vec1;
    CHECK(vec_neg(0, 0) == -1);
    CHECK(vec_neg(1, 0) == -2);
    CHECK(vec_neg(2, 0) == -3);

    // Dot product via matrix multiplication (vec1^T * vec2)
    Matrix<1, 3> vec1_T = vec1.transpose();
    auto         dot_product = vec1_T * vec2;
    CHECK(dot_product(0, 0) == 32); // 1*4 + 2*5 + 3*6

    // Test Matrix::dot member function
    auto direct_dot = vec1.dot(vec2);
    CHECK(direct_dot == 32); // 1*4 + 2*5 + 3*6
}

TEST_CASE("Matrix zeros() static method") {
    auto zero_mat = Matrix<3, 2, float>::zeros();
    CHECK(zero_mat(0, 0) == 0.0f);
    CHECK(zero_mat(1, 0) == 0.0f);
    CHECK(zero_mat(2, 0) == 0.0f);
    CHECK(zero_mat(0, 1) == 0.0f);
    CHECK(zero_mat(1, 1) == 0.0f);
    CHECK(zero_mat(2, 1) == 0.0f);

    // Test with different types
    auto zero_int = Matrix<2, 2, double>::zeros();
    CHECK(zero_int(0, 0) == 0);
    CHECK(zero_int(1, 1) == 0);
}

TEST_CASE("Matrix identity() static method") {
    auto identity_mat = Matrix<3, 3, float>::identity();
    CHECK(identity_mat(0, 0) == 1.0f);
    CHECK(identity_mat(1, 1) == 1.0f);
    CHECK(identity_mat(2, 2) == 1.0f);
    CHECK(identity_mat(0, 1) == 0.0f);
    CHECK(identity_mat(0, 2) == 0.0f);
    CHECK(identity_mat(1, 0) == 0.0f);
    CHECK(identity_mat(1, 2) == 0.0f);
    CHECK(identity_mat(2, 0) == 0.0f);
    CHECK(identity_mat(2, 1) == 0.0f);

    // Test with different types
    auto identity_int = Matrix<2, 2, double>::identity();
    CHECK(identity_int(0, 0) == 1);
    CHECK(identity_int(1, 1) == 1);
    CHECK(identity_int(0, 1) == 0);
    CHECK(identity_int(1, 0) == 0);
}

TEST_CASE("Matrix diagonal() static method") {
    auto diag_mat = Matrix<3, 3, float>::diagonal({1.0f, 2.0f, 3.0f});
    CHECK(diag_mat(0, 0) == 1.0f);
    CHECK(diag_mat(1, 1) == 2.0f);
    CHECK(diag_mat(2, 2) == 3.0f);
    CHECK(diag_mat(0, 1) == 0.0f);
    CHECK(diag_mat(0, 2) == 0.0f);
    CHECK(diag_mat(1, 0) == 0.0f);
    CHECK(diag_mat(1, 2) == 0.0f);
    CHECK(diag_mat(2, 0) == 0.0f);
    CHECK(diag_mat(2, 1) == 0.0f);

    // Test with different types
    auto diag_int = Matrix<2, 2, double>::diagonal({5.0, 7.0});
    CHECK(diag_int(0, 0) == 5.0);
    CHECK(diag_int(1, 1) == 7.0);
    CHECK(diag_int(0, 1) == 0.0);
    CHECK(diag_int(1, 0) == 0.0);

    // Test compile-time usage
    static_assert(Matrix<2, 2, double>::diagonal({1.0, 2.0})(0, 0) == 1.0);
    static_assert(Matrix<2, 2, double>::diagonal({1.0, 2.0})(1, 1) == 2.0);
    static_assert(Matrix<2, 2, double>::diagonal({1.0, 2.0})(0, 1) == 0.0);
}

TEST_CASE("Matrix fill, data, dims") {
    Matrix<2, 3, double> mat;
    mat.fill(7);
    CHECK(mat(0, 0) == 7);
    CHECK(mat(1, 2) == 7);

    auto* ptr = mat.data();
    CHECK(ptr[0] == 7);
    CHECK(ptr[5] == 7);

    static_assert(Matrix<2, 3, double>::rows() == 2);
    static_assert(Matrix<2, 3, double>::cols() == 3);
}

TEST_CASE("Matrix inverse() - 1x1 matrix") {
    Matrix<1, 1, float> mat = {{2.0f}};
    auto                inv = mat.inverse();
    REQUIRE(inv.has_value());
    CHECK(inv->operator()(0, 0) == 0.5f);

    // Test identity property: mat * inv = identity
    auto product = mat * inv.value();
    CHECK(product(0, 0) == doctest::Approx(1.0));
}

TEST_CASE("Matrix inverse() - 2x2 matrix") {
    Matrix<2, 2, float> mat = {{4.0f, 2.0f}, {7.0f, 6.0f}};
    auto                inv = mat.inverse();
    REQUIRE(inv.has_value());

    // Expected inverse of [[4, 2], [7, 6]] is [[0.6, -0.2], [-0.7, 0.4]]
    CHECK(inv->operator()(0, 0) == doctest::Approx(0.6));
    CHECK(inv->operator()(0, 1) == doctest::Approx(-0.2));
    CHECK(inv->operator()(1, 0) == doctest::Approx(-0.7));
    CHECK(inv->operator()(1, 1) == doctest::Approx(0.4));

    // Test identity property: mat * inv = identity
    auto product = mat * inv.value();
    CHECK(product(0, 0) == doctest::Approx(1.0));
    CHECK(product(0, 1) == doctest::Approx(0.0));
    CHECK(product(1, 0) == doctest::Approx(0.0));
    CHECK(product(1, 1) == doctest::Approx(1.0));
}

TEST_CASE("Matrix inverse() - singular matrix returns nullopt") {
    Matrix<2, 2, float> singular = {
        {1.0f, 2.0f},
        {2.0f, 4.0f} // determinant = 0
    };
    auto inv = singular.inverse();
    CHECK_FALSE(inv.has_value());
}

TEST_CASE("Matrix-matrix multiplication - various dimensions") {
    // 2x3 * 3x2 = 2x2
    Matrix<2, 3> mat1 = {{1, 2, 3}, {4, 5, 6}};
    Matrix<3, 2> mat2 = {{7, 8}, {9, 10}, {11, 12}};
    auto         result = mat1 * mat2;
    CHECK(result(0, 0) == 58);  // 1*7 + 2*9 + 3*11
    CHECK(result(0, 1) == 64);  // 1*8 + 2*10 + 3*12
    CHECK(result(1, 0) == 139); // 4*7 + 5*9 + 6*11
    CHECK(result(1, 1) == 154); // 4*8 + 5*10 + 6*12

    // 3x1 * 1x3 = 3x3 (outer product)
    Matrix<3, 1> vec1 = {{1}, {2}, {3}};
    Matrix<1, 3> vec2 = {{4, 5, 6}};
    auto         outer = vec1 * vec2;
    CHECK(outer(0, 0) == 4);
    CHECK(outer(0, 1) == 5);
    CHECK(outer(0, 2) == 6);
    CHECK(outer(1, 0) == 8);
    CHECK(outer(1, 1) == 10);
    CHECK(outer(1, 2) == 12);
    CHECK(outer(2, 0) == 12);
    CHECK(outer(2, 1) == 15);
    CHECK(outer(2, 2) == 18);
}

TEST_CASE("is_matrix_type type trait") {
    // Test that Matrix types are detected
    CHECK(is_matrix_type<Matrix<2, 2, double>>::value == true);
    CHECK(is_matrix_type<Matrix<3, 1, float>>::value == true);
    CHECK(is_matrix_type<Matrix<1, 5, double>>::value == true);

    // Test that non-Matrix types are not detected
    CHECK(is_matrix_type<int>::value == false);
    CHECK(is_matrix_type<float>::value == false);
    CHECK(is_matrix_type<std::vector<int>>::value == false);

    // Test with typedefs (ColVec, RowVec)
    using ColVec2f = Matrix<2, 1, float>;
    using RowVec2f = Matrix<1, 2, float>;
    CHECK(is_matrix_type<ColVec2f>::value == true);
    CHECK(is_matrix_type<RowVec2f>::value == true);
}

TEST_CASE("Matrix operations with floating point precision") {
    // Test inverse with floating point precision
    Matrix<2, 2> mat = {{1.0, 0.5}, {0.5, 1.0}};
    auto         inv = mat.inverse();
    CHECK(inv.has_value());

    // Expected inverse of [[1, 0.5], [0.5, 1]] is [[1.333..., -0.666...], [-0.666..., 1.333...]]
    CHECK(inv->operator()(0, 0) == doctest::Approx(1.333333333).epsilon(1e-6));
    CHECK(inv->operator()(0, 1) == doctest::Approx(-0.666666667).epsilon(1e-6));
    CHECK(inv->operator()(1, 0) == doctest::Approx(-0.666666667).epsilon(1e-6));
    CHECK(inv->operator()(1, 1) == doctest::Approx(1.333333333).epsilon(1e-6));

    // Test identity property
    auto product = mat * inv.value();
    CHECK(product(0, 0) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(product(0, 1) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(product(1, 0) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(product(1, 1) == doctest::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("Matrix trace") {
    Matrix<2, 2> mat = {{3, 1}, {2, 5}};
    CHECK(mat.trace() == 8);
}

TEST_CASE("Block Assignment") {
    Matrix<4, 4> mat = Matrix<4, 4>::zeros();
    Matrix<2, 2> block = {{1, 2}, {3, 4}};

    mat.block<2, 2>(0, 0) = block;

    CHECK(mat(0, 0) == 1);
    CHECK(mat(0, 1) == 2);
    CHECK(mat(1, 0) == 3);
    CHECK(mat(1, 1) == 4);
}

TEST_CASE("Block view operations") {
    Matrix<4, 4, double> mat = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};

    // Create a 2x2 block starting at (1,1)
    auto block = mat.block<2, 2>(1, 1);

    // Check access
    CHECK(block(0, 0) == 6);
    CHECK(block(0, 1) == 7);
    CHECK(block(1, 0) == 10);
    CHECK(block(1, 1) == 11);

    // Modify through block
    block(0, 0) = 99;
    CHECK(mat(1, 1) == 99);

    // Assign a matrix to block
    Matrix<2, 2, double> sub = {{100, 101}, {102, 103}};
    block = sub;
    CHECK(mat(1, 1) == 100);
    CHECK(mat(1, 2) == 101);
    CHECK(mat(2, 1) == 102);
    CHECK(mat(2, 2) == 103);

    // Test copy construction
    auto block2 = block;
    block2(0, 0) = 200;
    CHECK(mat(1, 1) == 200);

    // Test copy assignment
    Matrix<4, 4, double> mat2 = {{21, 22, 23, 24}, {25, 26, 27, 28}, {29, 30, 31, 32}, {33, 34, 35, 36}};
    auto                 block3 = mat2.block<2, 2>(1, 1);
    block3 = block; // Assign block to block3
    CHECK(mat2(1, 1) == 200);
    CHECK(mat2(1, 2) == 101);
    CHECK(mat2(2, 1) == 102);
    CHECK(mat2(2, 2) == 103);

    // Test const block
    const auto& const_mat = mat;
    auto        const_block = const_mat.block<2, 2>(1, 1);
    CHECK(const_block(0, 0) == 200);
    // const_block(0, 0) = 300; // This should not compile if const-correct

    // Test move operations (basic check)
    auto block4 = damp::move(block3);
    CHECK(block4(0, 0) == 200);
}

TEST_SUITE("Matrix Inversion Tests") {
    /**
     * @brief Helper function that generates a random matrix of given size
     *
     * @param rng Random number generator
     *
     * @return Randomly generated matrix (rows x cols) [-10.0, 10.0]
     */
    template<size_t rows, size_t cols>
    auto make_random_matrix(std::mt19937 & rng) {
        Matrix<rows, cols> mat;
        for (size_t r = 0; r < rows; ++r) {
            for (size_t c = 0; c < cols; ++c) {
                mat(r, c) = std::uniform_real_distribution<double>(-10.0, 10.0)(rng);
            }
        }
        return mat;
    };

    template<size_t N>
    void test_inverse_for_size() {
        std::mt19937 rng(0x1234); // Fixed seed for reproducibility

        constexpr size_t num_tests = 10;

        for (size_t i = 0; i < num_tests; ++i) {
            // Generate a random matrix
            auto mat = make_random_matrix<N, N>(rng);

            // Take its inverse
            auto inv = mat.inverse();
            REQUIRE(inv.has_value());

            // Verify that mat * inv = identity
            auto ident = mat * inv.value();
            for (size_t r = 0; r < N; ++r) {
                for (size_t c = 0; c < N; ++c) {
                    if (r == c) {
                        CHECK(ident(r, c) == doctest::Approx(1.0).epsilon(1e-10));
                    } else {
                        CHECK(ident(r, c) == doctest::Approx(0.0).epsilon(1e-10));
                    }
                }
            }
        }
    }

    template<size_t... I>
    void test_all_sizes(damp::index_sequence<I...>) {
        // I = 0..8  →  N = 1..9
        (test_inverse_for_size<I + 1>(), ...);
    }

    TEST_CASE("Matrix Inversion for sizes 1 to 9") {
        test_all_sizes(damp::make_index_sequence<9>{});
    }
}

TEST_CASE("Matrix sum and mean") {
    Matrix<2, 3> mat = {{1, 2, 3}, {4, 5, 6}};

    CHECK(mat.sum() == 21);
    CHECK(mat.mean() == doctest::Approx(3.5));

    Matrix<1, 1> scalar_mat = {{5.0}};
    CHECK(scalar_mat.sum() == 5.0);
    CHECK(scalar_mat.mean() == 5.0);
}

TEST_CASE("Block scalar assignment") {
    Matrix<3, 3> mat = Matrix<3, 3>::zeros();

    auto block = mat.block<2, 2>(0, 0);
    block = 7.0;

    CHECK(mat(0, 0) == 7);
    CHECK(mat(0, 1) == 7);
    CHECK(mat(1, 0) == 7);
    CHECK(mat(1, 1) == 7);
    CHECK(mat(2, 2) == 0); // Unchanged
}
TEST_CASE("Broadcasting addition") {
    Matrix<1, 3> row = {{1, 2, 3}};
    Matrix<2, 3> mat = {{4, 5, 6}, {7, 8, 9}};

    auto result = mat + row; // Should broadcast row to both rows of mat

    CHECK(result.rows() == 2);
    CHECK(result.cols() == 3);
    CHECK(result(0, 0) == 5);  // 4 + 1
    CHECK(result(0, 1) == 7);  // 5 + 2
    CHECK(result(0, 2) == 9);  // 6 + 3
    CHECK(result(1, 0) == 8);  // 7 + 1
    CHECK(result(1, 1) == 10); // 8 + 2
    CHECK(result(1, 2) == 12); // 9 + 3
}

TEST_CASE("Matrix reshape") {
    Matrix<2, 3> mat = {{1, 2, 3}, {4, 5, 6}};

    auto reshaped = mat.reshape<3, 2>();

    CHECK(reshaped.rows() == 3);
    CHECK(reshaped.cols() == 2);
    CHECK(reshaped(0, 0) == 1);
    CHECK(reshaped(0, 1) == 2);
    CHECK(reshaped(1, 0) == 3);
    CHECK(reshaped(1, 1) == 4);
    CHECK(reshaped(2, 0) == 5);
    CHECK(reshaped(2, 1) == 6);
}

TEST_CASE("Transpose assignment does not drop the lower triangle") {
    Matrix<2, 2> A{{1.0, 2.0}, {3.0, 4.0}};
    A.t() = A;
    CHECK(A(0, 0) == doctest::Approx(1.0));
    CHECK(A(0, 1) == doctest::Approx(3.0));
    CHECK(A(1, 0) == doctest::Approx(2.0));
    CHECK(A(1, 1) == doctest::Approx(4.0));
}

TEST_CASE("Overlapping block assignment copies through a temporary") {
    Matrix<3, 3> A{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
    A.block<2, 2>(0, 0) = A.block<2, 2>(1, 0);
    CHECK(A(0, 0) == doctest::Approx(4.0));
    CHECK(A(0, 1) == doctest::Approx(5.0));
    CHECK(A(1, 0) == doctest::Approx(7.0));
    CHECK(A(1, 1) == doctest::Approx(8.0));
}

TEST_CASE("Disjoint blocks of one matrix assign in place") {
    Matrix<2, 4> A{{1.0, 2.0, 3.0, 4.0}, {5.0, 6.0, 7.0, 8.0}};
    A.block<2, 2>(0, 0) = A.block<2, 2>(0, 2);
    CHECK(A(0, 0) == doctest::Approx(3.0));
    CHECK(A(0, 1) == doctest::Approx(4.0));
    CHECK(A(1, 0) == doctest::Approx(7.0));
    CHECK(A(1, 1) == doctest::Approx(8.0));
    CHECK(A(0, 2) == doctest::Approx(3.0));
    CHECK(A(1, 3) == doctest::Approx(8.0));
}

TEST_CASE("View of a temporary owns its elements") {
    Matrix<2, 2> A{{1.0, 2.0}, {3.0, 4.0}};
    Matrix<2, 2> B{{5.0, 6.0}, {7.0, 8.0}};
    const auto   block = (A + B).block<2, 2>(0, 0);
    CHECK(block(0, 0) == doctest::Approx(6.0));
    CHECK(block(1, 1) == doctest::Approx(12.0));
    const auto row = (A + B).row(1);
    CHECK(row(0, 0) == doctest::Approx(10.0));
    CHECK(row(0, 1) == doctest::Approx(12.0));
}

TEST_CASE("Block span-like methods") {
    Matrix<3, 3> mat = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

    auto block = mat.block<2, 2>(1, 1);

    CHECK(block.size() == 4);
    CHECK(block.rows() == 2);
    CHECK(block.cols() == 2);
    CHECK(block.data() == &mat(1, 1));
}

TEST_CASE("Block arithmetic operators") {
    Matrix<3, 3> mat1 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    Matrix<3, 3> mat2 = {{9, 8, 7}, {6, 5, 4}, {3, 2, 1}};

    auto block1 = mat1.block<2, 2>(0, 0);
    auto block2 = mat2.block<2, 2>(0, 0);

    auto sum = block1 + block2;
    CHECK(sum.rows() == 2);
    CHECK(sum.cols() == 2);
    CHECK(sum(0, 0) == 1 + 9);
    CHECK(sum(0, 1) == 2 + 8);
    CHECK(sum(1, 0) == 4 + 6);
    CHECK(sum(1, 1) == 5 + 5);

    auto diff = block1 - block2;
    CHECK(diff(0, 0) == 1 - 9);
    CHECK(diff(1, 1) == 5 - 5);

    auto prod = block1 * block2; // Matrix multiplication
    CHECK(prod(0, 0) == 21);     // 1*9 + 2*6
    CHECK(prod(1, 1) == 57);     // 4*8 + 5*5

    auto scalar_add = block1 + 10.0;
    CHECK(scalar_add(0, 0) == 1 + 10);
    CHECK(scalar_add(1, 1) == 5 + 10);

    auto scalar_mul = block1 * 2.0;
    CHECK(scalar_mul(0, 0) == 1 * 2);
    CHECK(scalar_mul(1, 1) == 5 * 2);
}

TEST_CASE("Block from block") {
    Matrix<4, 4> mat = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}};

    auto outer_block = mat.block<3, 3>(0, 0);
    auto inner_block = outer_block.block<2, 2>(1, 1);

    CHECK(inner_block.rows() == 2);
    CHECK(inner_block.cols() == 2);
    CHECK(inner_block(0, 0) == 6);  // mat(1,1)
    CHECK(inner_block(0, 1) == 7);  // mat(1,2)
    CHECK(inner_block(1, 0) == 10); // mat(2,1)
    CHECK(inner_block(1, 1) == 11); // mat(2,2)

    // Modify through inner block
    inner_block = 99.0;
    CHECK(mat(1, 1) == 99);
    CHECK(mat(2, 2) == 99);
}

TEST_CASE("Block compound assignments") {
    Matrix<3, 3> mat1 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    Matrix<3, 3> mat2 = {{1, 1, 1}, {2, 2, 2}, {3, 3, 3}};

    auto block1 = mat1.block<2, 2>(0, 0);
    auto block2 = mat2.block<2, 2>(0, 0);

    // Test += with another block (direct view arithmetic)
    block1 = block1 + block2;
    CHECK(mat1(0, 0) == 1 + 1);
    CHECK(mat1(0, 1) == 2 + 1);
    CHECK(mat1(1, 0) == 4 + 2);
    CHECK(mat1(1, 1) == 5 + 2);

    // Test -= with scalar
    block1 = block1 - 1.0;
    CHECK(mat1(0, 0) == 2 - 1);
    CHECK(mat1(1, 1) == 7 - 1);

    // Test *= with scalar
    block1 = block1 * 2.0;
    CHECK(mat1(0, 0) == 1 * 2);
    CHECK(mat1(1, 1) == 6 * 2);

    // Test /= with scalar
    block1 = block1 / 2.0;
    CHECK(mat1(0, 0) == 1);
    CHECK(mat1(1, 1) == 6);
}

TEST_CASE("Matrix inverse") {
    // Test with a simple 2x2 matrix
    Matrix<2, 2, double> A = {
        {2.0, 1.0},
        {1.0, 1.0},
    };

    auto A_inv_opt = A.inverse();
    REQUIRE(A_inv_opt.has_value());
    auto A_inv = *A_inv_opt;

    // Verify A * A_inv ≈ I
    auto product = A * A_inv;
    auto identity = Matrix<2, 2, double>::identity();

    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            CHECK(doctest::Approx(product(i, j)).epsilon(1e-12) == identity(i, j));
        }
    }

    // Test with singular matrix (should return nullopt)
    Matrix<2, 2, double> singular = {
        {1.0, 2.0},
        {2.0, 4.0} // rank deficient
    };

    auto singular_inv = singular.inverse();
    CHECK(!singular_inv.has_value());
}

// 3x3 matrix test
TEST_CASE("Matrix inverse 3x3") {
    using T = double;
    Matrix<3, 3, T> A{
        {1.496714153011233e+00, -1.382643011711847e-01, 6.476885381006925e-01},
        {1.523029856408025e+00, 7.658466252766640e-01, -2.341369569491805e-01},
        {1.579212815507391e+00, 7.674347291529088e-01, 5.305256140650478e-01},
    };

    Matrix<3, 3, T> expected_inv{
        {5.781271633034201e-01, 5.627613896905597e-01, -4.574389841541137e-01},
        {-1.161962990634306e+00, -2.257244021189988e-01, 1.318955517996783e+00},
        {-4.006417732725574e-02, -1.348646011238505e+00, 1.338637035691151e+00},
    };

    auto actual_inv_opt = A.inverse();
    REQUIRE(actual_inv_opt.has_value());

    Matrix<3, 3, T> actual_inv = actual_inv_opt.value();
    Matrix<3, 3, T> product = A * actual_inv;
    Matrix<3, 3, T> identity = Matrix<3, 3, T>::identity();
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            CHECK(product(i, j) == doctest::Approx(identity(i, j)).epsilon(1e-10));
        }
    }

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            CHECK(actual_inv(i, j) == doctest::Approx(expected_inv(i, j)).epsilon(1e-10));
        }
    }
}

// 4x4 matrix test
TEST_CASE("Matrix inverse 4x4") {
    using T = double;
    Matrix<4, 4, T> A{
        {1.542560043585965e+00, -4.634176928124623e-01, -4.657297535702569e-01, 2.419622715660341e-01},
        {-1.913280244657798e+00, -7.249178325130328e-01, -5.622875292409727e-01, -1.012831120334424e+00},
        {3.142473325952739e-01, -9.080240755212109e-01, -4.123037013352915e-01, 1.465648768921554e+00},
        {-2.257763004865357e-01, 6.752820468792384e-02, -1.424748186213457e+00, 4.556172754748173e-01},
    };

    Matrix<4, 4, T> expected_inv{
        {5.037582970690051e-01, -1.340823238222598e-01, -1.550277395898095e-01, -6.689174905130418e-02},
        {-4.731028803935471e-01, -5.065999237407754e-01, -4.199895986719887e-01, 4.761235563953082e-01},
        {-2.533124665129921e-01, -1.032239124411711e-01, 1.651250224224445e-01, -6.261211597034416e-01},
        {-4.723745476684392e-01, -3.141472402757916e-01, 5.017832785754724e-01, 1.331834327233001e-01},
    };

    auto actual_inv_opt = A.inverse();
    REQUIRE(actual_inv_opt.has_value());

    Matrix<4, 4, T> actual_inv = actual_inv_opt.value();
    Matrix<4, 4, T> product = A * actual_inv;
    Matrix<4, 4, T> identity = Matrix<4, 4, T>::identity();
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            CHECK(product(i, j) == doctest::Approx(identity(i, j)).epsilon(1e-10));
        }
    }
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            CHECK(actual_inv(i, j) == doctest::Approx(expected_inv(i, j)).epsilon(1e-10));
        }
    }
}

// 5x5 matrix test
TEST_CASE("Matrix inverse 5x5") {
    using T = double;
    Matrix<5, 5, T> A{
        {1.110922589709866e+00, -1.150993577422303e+00, 3.756980183456720e-01, -6.006386899188050e-01, -2.916937497932768e-01},
        {-6.017066122293969e-01, 2.852278184508938e+00, -1.349722473793392e-02, -1.057710928955900e+00, 8.225449121031890e-01},
        {-1.220843649971022e+00, 2.088635950047554e-01, -9.596701238797756e-01, -1.328186048898431e+00, 1.968612358691235e-01},
        {7.384665799954104e-01, 1.713682811899705e-01, -1.156482823882405e-01, 6.988963044107113e-01, -1.478521990367427e+00},
        {-7.198442083947086e-01, -4.606387709597875e-01, 1.057122226218916e+00, 3.436182895684614e-01, -7.630401553627340e-01},
    };

    Matrix<5, 5, T> expected_inv{
        {3.552909350106445e-01, 9.853931163143488e-02, -2.609975195988257e-01, 1.294609088927118e-01, -3.477854065332936e-01},
        {-9.920046853805239e-02, 3.012035005445296e-01, -1.033366003429114e-01, 1.908400986567891e-01, -3.383136477268624e-02},
        {2.927238385928023e-01, 2.418046371159175e-01, -3.133630255020795e-01, -2.406511667581695e-01, 5.342162749473734e-01},
        {-5.725951358265908e-01, -2.239482798413569e-01, -3.444302598780314e-01, 3.203817369922038e-04, -1.120042678596126e-01},
        {-1.276051129676750e-01, -4.064613301355580e-02, -2.806368143768866e-01, -5.705959758948312e-01, -2.723571841992745e-01},
    };

    auto actual_inv_opt = A.inverse();
    REQUIRE(actual_inv_opt.has_value());

    Matrix<5, 5, T> actual_inv = actual_inv_opt.value();
    Matrix<5, 5, T> product = A * actual_inv;
    Matrix<5, 5, T> identity = Matrix<5, 5, T>::identity();
    for (size_t i = 0; i < 5; ++i) {
        for (size_t j = 0; j < 5; ++j) {
            CHECK(product(i, j) == doctest::Approx(identity(i, j)).epsilon(1e-10));
        }
    }

    for (size_t i = 0; i < 5; ++i) {
        for (size_t j = 0; j < 5; ++j) {
            CHECK(actual_inv(i, j) == doctest::Approx(expected_inv(i, j)).epsilon(1e-10));
        }
    }
}

// 6x6 matrix test
TEST_CASE("Matrix inverse 6x6") {
    using T = double;
    Matrix<6, 6, T> A{
        {1.324083969394795e+00, -3.850822804163165e-01, -6.769220003059587e-01, 6.116762888408679e-01, 1.030999522495951e+00, 9.312801191161986e-01},
        {-8.392175232226385e-01, 6.907876241487854e-01, 3.312634314035640e-01, 9.755451271223592e-01, -4.791742378452900e-01, -1.856589766638171e-01},
        {-1.106334974006028e+00, -1.196206624080671e+00, 1.812525822394198e+00, 1.356240028570823e+00, -7.201012158033385e-02, 1.003532897892024e+00},
        {3.616360250476341e-01, -6.451197546051243e-01, 3.613956055084139e-01, 2.538036566465969e+00, -3.582603910995154e-02, 1.564643655814006e+00},
        {-2.619745104089744e+00, 8.219025043752238e-01, 8.704706823817122e-02, -2.990073504658675e-01, 1.091760776535502e+00, -1.987568914600893e+00},
        {-2.196718878375119e-01, 3.571125715117464e-01, 1.477894044741516e+00, -5.182702182736474e-01, -8.084936028931876e-01, 4.982429564154635e-01},
    };

    Matrix<6, 6, T> expected_inv{
        {9.219593018189169e+00, 1.583755605643359e+01, 8.344834809301510e+00, -1.461169029797190e+01, -5.863293654960798e+00, -5.642964213805874e+00},
        {-8.348702197662898e-01, -1.800590441808834e+00, -1.647832822540016e+00, 2.190624502374752e+00, 1.031908610057714e+00, 1.445679147918631e+00},
        {8.234790794532385e+00, 1.400916099252166e+01, 7.654966074641130e+00, -1.298647356251552e+01, -4.964382008258302e+00, -4.611863761711620e+00},
        {4.527505650232546e+00, 8.276275058376543e+00, 4.258789160878783e+00, -7.125420735593321e+00, -2.885783067825769e+00, -3.092055035706307e+00},
        {2.082805881009056e+00, 2.255382098288948e+00, 1.332351492490935e+00, -2.262233966949203e+00, -4.941536023906656e-01, -6.032732461210786e-01},
        {-1.167364951804336e+01, -2.101217341412075e+01, -1.125402782479168e+01, 1.942558268741152e+01, 7.597070438874137e+00, 7.967423400536967e+00},
    };

    auto actual_inv_opt = A.inverse();
    REQUIRE(actual_inv_opt.has_value());

    Matrix<6, 6, T> actual_inv = actual_inv_opt.value();
    Matrix<6, 6, T> product = A * actual_inv;
    Matrix<6, 6, T> identity = Matrix<6, 6, T>::identity();
    for (size_t i = 0; i < 6; ++i) {
        for (size_t j = 0; j < 6; ++j) {
            CHECK(product(i, j) == doctest::Approx(identity(i, j)).epsilon(1e-10));
        }
    }

    for (size_t i = 0; i < 6; ++i) {
        for (size_t j = 0; j < 6; ++j) {
            CHECK(actual_inv(i, j) == doctest::Approx(expected_inv(i, j)).epsilon(1e-10));
        }
    }
}

// 7x7 matrix test
TEST_CASE("Matrix inverse 7x7") {
    using T = double;
    Matrix<7, 7, T> A{
        {1.915402117702074e+00, 3.287511096596845e-01, -5.297602037670388e-01, 5.132674331133561e-01, 9.707754934804039e-02, 9.686449905328892e-01, -7.020530938773524e-01},
        {-3.276621465977682e-01, 6.078918468678424e-01, -1.463514948132119e+00, 2.961202770645761e-01, 2.610552721798893e-01, 5.113456642460890e-03, -2.345871333751469e-01},
        {-1.415370742050414e+00, -4.206453227653590e-01, 6.572854834732305e-01, -8.022772692216189e-01, -1.612857116660091e-01, 4.040508568145384e-01, 1.886185901210530e+00},
        {1.745778128318390e-01, 2.575503907227644e-01, -7.444591576616721e-02, -9.187712152990415e-01, -2.651387544921688e-02, 6.023020994102644e-02, 2.463242112485286e+00},
        {-1.923609647811225e-01, 3.015473423336125e-01, -3.471176970524331e-02, -1.168678037619532e+00, 2.142822814515021e+00, 7.519330326867741e-01, 7.910319470430469e-01},
        {-9.093874547947389e-01, 1.402794310936099e+00, -1.401851062792281e+00, 5.868570938002703e-01, 2.190455625809979e+00, 9.463674869311656e-03, -5.662977296027719e-01},
        {9.965136508764122e-02, -5.034756541161992e-01, -1.550663431066133e+00, 6.856297480602733e-02, -1.062303713726105e+00, 4.735924306351816e-01, 8.057576576619685e-02},
    };

    Matrix<7, 7, T> expected_inv{
        {5.524303285598250e-02, -7.247828353655729e-01, -5.130477785074442e-01, 3.688362568598772e-01, 6.863328820953631e-03, 1.754466692520696e-01, 2.712383595586376e-01},
        {6.178292293116449e-01, 3.718797764864832e+00, 7.427923947154822e-01, -3.522761277566802e-01, 1.049127107557167e-01, -1.506878632217408e+00, -2.029187965552086e+00},
        {2.286591796168871e-01, 2.764205690394361e-01, 3.443624640835058e-01, -1.818301330450062e-01, -7.893126721933388e-02, -2.391385945259077e-01, -6.112368764452979e-01},
        {1.174207832263735e-01, -3.092740696609128e+00, -4.420411967405580e-02, 4.519144559671758e-01, -9.176196053829432e-01, 1.842116964953159e+00, 1.193584635242779e+00},
        {-2.714553174267275e-01, -2.077338632640303e+00, -4.803392557998336e-01, 2.978042885828822e-01, 3.764915904409664e-02, 1.068704492755904e+00, 8.684416022294921e-01},
        {7.747341127711851e-01, 1.059613730053016e+00, 9.670732195025066e-01, -5.463505279535646e-01, 1.300017537626713e-01, -4.355117251329946e-01, -4.377849397937827e-01},
        {-3.967150438433407e-02, -1.530945170711477e+00, -7.620017678490418e-02, 5.962915524434225e-01, -3.588802287843774e-01, 8.471419376657615e-01, 6.397205918165388e-01},
    };

    auto actual_inv_opt = A.inverse();
    REQUIRE(actual_inv_opt.has_value());

    Matrix<7, 7, T> actual_inv = actual_inv_opt.value();
    Matrix<7, 7, T> product = A * actual_inv;
    Matrix<7, 7, T> identity = Matrix<7, 7, T>::identity();
    for (size_t i = 0; i < 7; ++i) {
        for (size_t j = 0; j < 7; ++j) {
            CHECK(product(i, j) == doctest::Approx(identity(i, j)).epsilon(1e-10));
        }
    }

    for (size_t i = 0; i < 7; ++i) {
        for (size_t j = 0; j < 7; ++j) {
            CHECK(actual_inv(i, j) == doctest::Approx(expected_inv(i, j)).epsilon(1e-10));
        }
    }
}

// 8x8 matrix test
TEST_CASE("Matrix inverse 8x8") {
    using T = double;
    Matrix<8, 8, T> A{
        {2.549934405017539e+00, -7.832532923362371e-01, -3.220615162056756e-01, 8.135172173696698e-01, -1.230864316433955e+00, 2.274599346041294e-01, 1.307142754282428e+00, -1.607483234561228e+00},
        {1.846338585323042e-01, 1.259882794248424e+00, 7.818228717773104e-01, -1.236950710878082e+00, -1.320456613084276e+00, 5.219415656168976e-01, 2.969846732331861e-01, 2.504928503458765e-01},
        {3.464482094969757e-01, -6.800247215784908e-01, 1.232253697161004e+00, 2.93072473986812e-01, -7.143514180263678e-01, 1.865774511144757e+00, 4.738329209117875e-01, -1.191303497202649e+00},
        {6.565536086338297e-01, -9.746816702273214e-01, 7.870846037424520e-01, 2.158595579007404e+00, -8.206823183517105e-01, 9.633761292443218e-01, 4.127809269364983e-01, 8.220601599944900e-01},
        {1.896792982653947e+00, -2.453881160028705e-01, -7.537361643574896e-01, -8.895144296255233e-01, 1.841897150345617e-01, -7.710170941410420e-02, 3.411519748166439e-01, 2.766907993300191e-01},
        {8.271832490360238e-01, 1.300189187790702e-02, 1.453534077157317e+00, -2.64656833379561e-01, 2.720169166589619e+00, 1.625667347765006e+00, -8.571575564162826e-01, -1.070892498061112e+00},
        {4.824724152431853e-01, -2.234627853258509e-01, 7.140004940920920e-01, 4.732376245735448e-01, -7.282891265687277e-02, -8.467937180684050e-01, -5.148472246858646e-01, -4.465149520670211e-01},
        {8.563987943234723e-01, 2.140937441302040e-01, -1.245738778711988e+00, 1.731809258511820e-01, 3.853173797288368e-01, -8.838574362011330e-01, 1.537251059455279e-01, 1.058208718446000e+00},
    };

    Matrix<8, 8, T> expected_inv{
        {-3.881975840554259e-01, 6.987217242906382e-02, 1.723443446787736e+00, -4.990672124881113e-01, -2.058551469082437e-01, -9.270709121974034e-02, 8.704617316613201e-01, 2.048969099829658e+00},
        {-2.721599759028239e-01, 4.242559724911424e-01, 1.909540449104829e+00, -7.948106822998497e-01, -1.237125770821691e+00, 3.635207536917446e-02, 6.567547327528404e-01, 2.890677390398517e+00},
        {1.657263400466631e+00, 4.304663259122828e-01, -6.500776357643364e+00, 2.288318972703170e+00, 1.561320386880343e+00, 7.965117954841815e-01, -2.127808992028737e+00, -7.180494691408549e+00},
        {-2.987597906239056e-01, -5.845506414338723e-02, 1.807878521999630e+00, -4.748822444734350e-01, -9.968200606400034e-01, -9.010980649554717e-02, 7.064512416818532e-01, 2.431712214936735e+00},
        {1.270551177971183e+00, 5.601692534403702e-02, -4.668518559017953e+00, 1.468991165220061e+00, 9.447263687182463e-01, 7.481887866892959e-01, -1.905850305741251e+00, -4.774125959881475e+00},
        {-1.375649735214174e+00, -1.480908024776835e-01, 5.099083977440352e+00, -1.535582615125854e+00, -1.127252078334652e+00, -4.586394096888384e-01, 1.494726694664588e+00, 5.339990270101679e+00},
        {3.665864952442371e+00, 5.273273854714353e-01, -1.261227198816661e+01, 4.136008252164986e+00, 2.667318811673531e+00, 1.431851310151671e+00, -5.249078387970734e+00, -1.343101989812065e+01},
        {2.249058165060260e-01, 1.532438643522488e-01, -1.938736014870986e+00, 9.179479619102422e-01, 7.450317307629486e-01, 1.565749057052497e-01, -7.528853401498666e-01, -1.999349149170280e+00},
    };

    auto actual_inv_opt = A.inverse();
    REQUIRE(actual_inv_opt.has_value());

    Matrix<8, 8, T> actual_inv = actual_inv_opt.value();
    Matrix<8, 8, T> product = A * actual_inv;
    Matrix<8, 8, T> identity = Matrix<8, 8, T>::identity();
    for (size_t i = 0; i < 8; ++i) {
        for (size_t j = 0; j < 8; ++j) {
            CHECK(product(i, j) == doctest::Approx(identity(i, j)).epsilon(1e-10));
        }
    }

    for (size_t i = 0; i < 8; ++i) {
        for (size_t j = 0; j < 8; ++j) {
            CHECK(actual_inv(i, j) == doctest::Approx(expected_inv(i, j)).epsilon(1e-8));
        }
    }
}

// 9x9 matrix test
TEST_CASE("Matrix inverse 9x9") {
    using T = double;
    Matrix<9, 9, T> A{{
        {-1.429702978306231e-01, 3.577873603482833e-01, 5.607845263682344e-01, 1.083051243175277e+00, 1.053802052034903e+00, -1.377669367957091e+00, -9.378250399151228e-01, 5.150352672086598e-01, 5.137859509122088e-01},
        {5.150476863060479e-01, 4.852731490654721e+00, 5.708905106931670e-01, 1.135565640180599e+00, 9.540017634932023e-01, 6.513912513057980e-01, -3.152692446403456e-01, 7.589692204932674e-01, -7.728252145375718e-01},
        {-2.368186067400089e-01, -4.853635478291035e-01, 1.081874139386322e+00, 2.314658566673509e+00, -1.867265192591748e+00, 6.862601903745135e-01, -1.612715871189652e+00, -4.719318657894335e-01, 1.088950596967366e+00},
        {6.428001909546277e-02, -1.077744777929306e+00, -7.153037092599682e-01, 1.679597748934676e+00, -7.303666317171367e-01, 2.164585895819749e-01, 4.557183990381378e-02, -6.516003476058171e-01, 2.143944089325326e+00},
        {6.339190223180112e-01, -2.025142586657607e+00, 1.864543147694276e-01, -6.617864647683880e-01, 1.852433334796224e+00, -7.925207384327007e-01, -1.147364414668990e-01, 5.049872789804571e-01, 8.657551941701215e-01},
        {-1.200296407055776e+00, -3.345012358409484e-01, -4.749453111609562e-01, -6.533292325737119e-01, 1.765454240281097e+00, 1.404981710960955e+00, -1.260883954335045e+00, 9.178619470547761e-01, 2.122156197012633e+00},
        {1.032465260551147e+00, -1.519369965954013e+00, -4.842340728662514e-01, 1.266911149186623e+00, -7.076694656187807e-01, 4.438194281462284e-01, 1.774634053429337e+00, -9.269304715780829e-01, -5.952535606180008e-02},
        {-3.241267340069073e+00, -1.024387641334290e+00, -2.525681513931603e-01, -1.247783181964849e+00, 1.632411303931635e+00, -1.430141377960633e+00, -4.400444866969838e-01, 1.130740577286091e+00, 1.441273289066116e+00},
        {-1.435862151179439e+00, 1.163163752154960e+00, 1.023306101958705e-02, -9.815086510479509e-01, 4.621034742632708e-01, 1.990596955734700e-01, -6.002168771587947e-01, 6.980208499001891e-02, 6.146864031382397e-01},
    }};
    Matrix<9, 9, T> expected_inv{{
        {-3.279775389485940e-01, 1.421681342010094e-01, -4.870437916932439e-02, 3.263446906053169e-01, 3.838478043616686e-01, -1.700305930954809e-01, -4.208776327749045e-01, -1.404173327423720e-01, -2.642105468546092e-01},
        {-1.635200269324444e-01, 2.013946169135109e-01, -4.321188902424791e-02, 2.163199244507130e-01, 9.193598783775307e-02, -1.777386409953094e-01, -1.593776986493337e-01, 8.581131294870656e-03, 1.605294160060367e-01},
        {-8.671620357996345e-01, 5.184209675126222e-01, 7.014208288052837e-01, -1.764594540709505e-01, 1.063020802536212e+00, -5.745135214647223e-01, 1.726318266551771e-01, 3.436161232004729e-01, 4.467545042941501e-01},
        {5.974414194682195e-01, -9.016285735363481e-02, -1.562671271079270e-02, -1.940649871867750e-01, -4.366117158674283e-01, 2.872552623939815e-01, 3.993551454746352e-01, -2.989640668465251e-02, -1.761824586071995e-01},
        {8.934890791031886e-01, -2.522541656973538e-01, -1.961940716409480e-01, -4.599443936446771e-01, -2.275586400228906e-01, 4.183090226502261e-01, 6.585514305739710e-01, -3.180813452105098e-01, 5.737329520127740e-01},
        {-1.127146202836684e-01, 2.098645783744640e-02, 1.405073606890025e-01, -3.030250472876553e-01, -8.527985780003769e-02, 3.491577623343067e-01, 3.160076388137164e-01, -6.691252360819307e-02, 3.075894923381855e-02},
        {-9.304330238029551e-01, 5.437917734732418e-01, 2.096158910802638e-01, 2.142563163349185e-01, 6.314663587661783e-01, -4.696384145712628e-01, 1.555847115594764e-01, 5.018686185376598e-01, -8.692775688583117e-02},
        {-1.288271269098638e+00, 6.497627032451500e-01, 2.335893054324069e-01, 4.102100439628416e-01, 4.168132801779376e-01, -2.275395365868569e-01, -8.589862222447939e-01, 7.636957351330317e-01, -1.626187723528615e+00},
        {-8.857319326276011e-01, 4.384454293852911e-01, 2.115174119645866e-01, 6.525741674084081e-01, 7.757671890040314e-01, -4.528915057168739e-01, -3.946981929263851e-01, 2.664245051373069e-01, 7.564952431955441e-02},
    }};

    auto actual_inv_opt = A.inverse();
    REQUIRE(actual_inv_opt.has_value());

    Matrix<9, 9, T> actual_inv = actual_inv_opt.value();
    Matrix<9, 9, T> product = A * actual_inv;
    Matrix<9, 9, T> identity = Matrix<9, 9, T>::identity();
    for (size_t i = 0; i < 9; ++i) {
        for (size_t j = 0; j < 9; ++j) {
            CHECK(product(i, j) == doctest::Approx(identity(i, j)).epsilon(1e-10));
        }
    }

    for (size_t i = 0; i < 9; ++i) {
        for (size_t j = 0; j < 9; ++j) {
            CHECK(actual_inv(i, j) == doctest::Approx(expected_inv(i, j)).epsilon(1e-10));
        }
    }
}
