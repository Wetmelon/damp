// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cstddef>
#include <type_traits>

#include "damp/math/complex.hpp"
#include "damp/matrix/matrix.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

TEST_CASE("Block view basic functionality") {
    Matrix<3, 4, double> mat = {
        {1, 2, 3, 4},
        {5, 6, 7, 8},
        {9, 10, 11, 12}
    };

    // Test block view
    auto block = mat.block<2, 2>(1, 1);
    CHECK(block(0, 0) == 6);
    CHECK(block(0, 1) == 7);
    CHECK(block(1, 0) == 10);
    CHECK(block(1, 1) == 11);

    // Test block modification
    block(0, 0) = 99;
    CHECK(mat(1, 1) == 99);
    CHECK(block(0, 0) == 99);

    // Test block arithmetic
    Matrix<3, 4, double> other_mat = {
        {100, 101, 102, 103},
        {104, 105, 106, 107},
        {108, 109, 110, 111}
    };
    auto other_block = other_mat.block<2, 2>(1, 1);
    auto sum = block + other_block;
    CHECK(sum(0, 0) == 99 + 105);
    CHECK(sum(0, 1) == 7 + 106);
    CHECK(sum(1, 0) == 10 + 109);
    CHECK(sum(1, 1) == 11 + 110);
}

TEST_CASE("Block view const correctness") {
    const Matrix<3, 4, double> mat = {
        {1, 2, 3, 4},
        {5, 6, 7, 8},
        {9, 10, 11, 12}
    };

    auto block = mat.block<2, 2>(1, 1);
    CHECK(block(0, 0) == 6);
    CHECK(block(0, 1) == 7);
    CHECK(block(1, 0) == 10);
    CHECK(block(1, 1) == 11);
}

TEST_CASE("Diagonal view basic functionality") {
    Matrix<3, 3, double> mat = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    };

    // Test diagonal view
    auto diag = mat.diagonal();
    CHECK(diag(0) == 1);
    CHECK(diag(1) == 5);
    CHECK(diag(2) == 9);

    // Test diagonal modification
    diag(1) = 99;
    CHECK(mat(1, 1) == 99);
    CHECK(diag(1) == 99);

    // Test conversion to vector
    auto vec = diag.to_vector();
    CHECK(vec(0) == 1);
    CHECK(vec(1) == 99);
    CHECK(vec(2) == 9);

    // Test assignment from vector
    ColVec<3, double> new_vec = {10, 20, 30};
    diag = new_vec;
    CHECK(mat(0, 0) == 10);
    CHECK(mat(1, 1) == 20);
    CHECK(mat(2, 2) == 30);
    CHECK(diag(0) == 10);
    CHECK(diag(1) == 20);
    CHECK(diag(2) == 30);
}

TEST_CASE("Diagonal view const correctness") {
    const Matrix<3, 3, double> mat = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    };

    auto diag = mat.diagonal();
    CHECK(diag(0) == 1);
    CHECK(diag(1) == 5);
    CHECK(diag(2) == 9);

    // Test conversion to vector
    auto vec = diag.to_vector();
    CHECK(vec(0) == 1);
    CHECK(vec(1) == 5);
    CHECK(vec(2) == 9);
}

TEST_CASE("Upper triangle view basic functionality") {
    Matrix<3, 3, double> mat = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    };

    auto upper = mat.upper_triangle();

    // On-triangle elements return stored values
    CHECK(upper(0, 0) == 1);
    CHECK(upper(0, 1) == 2);
    CHECK(upper(0, 2) == 3);
    CHECK(upper(1, 1) == 5);
    CHECK(upper(1, 2) == 6);
    CHECK(upper(2, 2) == 9);

    // Off-triangle elements return zero
    CHECK(upper(1, 0) == 0);
    CHECK(upper(2, 0) == 0);
    CHECK(upper(2, 1) == 0);
}

TEST_CASE("Lower triangle view basic functionality") {
    Matrix<3, 3, double> mat = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    };

    auto lower = mat.lower_triangle();

    // On-triangle elements return stored values
    CHECK(lower(0, 0) == 1);
    CHECK(lower(1, 0) == 4);
    CHECK(lower(1, 1) == 5);
    CHECK(lower(2, 0) == 7);
    CHECK(lower(2, 1) == 8);
    CHECK(lower(2, 2) == 9);

    // Off-triangle elements return zero
    CHECK(lower(0, 1) == 0);
    CHECK(lower(0, 2) == 0);
    CHECK(lower(1, 2) == 0);
}

TEST_CASE("Lower triangle view const correctness") {
    const Matrix<3, 3, double> mat = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    };

    auto lower = mat.lower_triangle();
    CHECK(lower(0, 0) == 1);
    CHECK(lower(1, 0) == 4);
    CHECK(lower(1, 1) == 5);
    CHECK(lower(2, 0) == 7);
    CHECK(lower(2, 1) == 8);
    CHECK(lower(2, 2) == 9);
}

TEST_CASE("RowView basic functionality") {
    Matrix<3, 4, double> mat = {
        {1, 2, 3, 4},
        {5, 6, 7, 8},
        {9, 10, 11, 12}
    };

    // Test row view
    auto row1 = mat.row(1);
    CHECK(row1(0) == 5);
    CHECK(row1(1) == 6);
    CHECK(row1(2) == 7);
    CHECK(row1(3) == 8);

    // Test row modification
    row1(0) = 99;
    CHECK(mat(1, 0) == 99);
    CHECK(row1(0) == 99);

    // Test to_vector
    auto vec = row1.to_vector();
    CHECK(vec(0) == 99);
    CHECK(vec(1) == 6);
    CHECK(vec(2) == 7);
    CHECK(vec(3) == 8);

    // Test assignment from RowVec
    RowVec<4, double> new_row{100, 101, 102, 103};
    row1 = new_row;
    CHECK(mat(1, 0) == 100);
    CHECK(mat(1, 1) == 101);
    CHECK(mat(1, 2) == 102);
    CHECK(mat(1, 3) == 103);
}

TEST_CASE("RowView const correctness") {
    const Matrix<3, 4, double> mat = {
        {1, 2, 3, 4},
        {5, 6, 7, 8},
        {9, 10, 11, 12}
    };

    auto row1 = mat.row(1);
    CHECK(row1(0) == 5);
    CHECK(row1(1) == 6);
    CHECK(row1(2) == 7);
    CHECK(row1(3) == 8);

    auto vec = row1.to_vector();
    CHECK(vec(0) == 5);
}

TEST_CASE("ColView basic functionality") {
    Matrix<3, 4, double> mat = {
        {1, 2, 3, 4},
        {5, 6, 7, 8},
        {9, 10, 11, 12}
    };

    // Test col view
    auto col2 = mat.col(2);
    CHECK(col2(0) == 3);
    CHECK(col2(1) == 7);
    CHECK(col2(2) == 11);
    CHECK(col2[1] == 7);

    // Test col modification
    col2(0) = 99;
    CHECK(mat(0, 2) == 99);
    CHECK(col2(0) == 99);

    // Test to_vector
    auto vec = col2.to_vector();
    CHECK(vec(0) == 99);
    CHECK(vec(1) == 7);
    CHECK(vec(2) == 11);

    // Test assignment from ColVec
    ColVec<3, double> new_col{100, 101, 102};
    col2 = new_col;
    CHECK(mat(0, 2) == 100);
    CHECK(mat(1, 2) == 101);
    CHECK(mat(2, 2) == 102);
}

TEST_CASE("ColView / RowView write-back without pack loops") {
    Matrix<3, 3, double> A = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    };
    Matrix<3, 3, double> B = Matrix<3, 3, double>::zeros();

    // Column view assign from matmul result (owning ColVec)
    const Matrix<3, 3, double> I = Matrix<3, 3, double>::identity();
    B.col(1) = I * A.col(1).to_vector();
    CHECK(B(0, 1) == doctest::Approx(2.0));
    CHECK(B(1, 1) == doctest::Approx(5.0));
    CHECK(B(2, 1) == doctest::Approx(8.0));
    CHECK(B(0, 0) == doctest::Approx(0.0));

    // Row view assign from column-shaped result (layout transfer)
    B.row(0) = ColVec<3, double>{10.0, 20.0, 30.0};
    CHECK(B(0, 0) == doctest::Approx(10.0));
    CHECK(B(0, 1) == doctest::Approx(20.0));
    CHECK(B(0, 2) == doctest::Approx(30.0));

    // Column view assign from row-shaped result
    B.col(2) = RowVec<3, double>{1.0, 2.0, 3.0};
    CHECK(B(0, 2) == doctest::Approx(1.0));
    CHECK(B(1, 2) == doctest::Approx(2.0));
    CHECK(B(2, 2) == doctest::Approx(3.0));

    // Row → column owning transfer
    const auto c0 = A.row(2).to_col_vector();
    CHECK(c0[0] == doctest::Approx(7.0));
    CHECK(c0[1] == doctest::Approx(8.0));
    CHECK(c0[2] == doctest::Approx(9.0));
}

TEST_CASE("ColView const correctness") {
    const Matrix<3, 4, double> mat = {
        {1, 2, 3, 4},
        {5, 6, 7, 8},
        {9, 10, 11, 12}
    };

    auto col2 = mat.col(2);
    CHECK(col2(0) == 3);
    CHECK(col2(1) == 7);
    CHECK(col2(2) == 11);

    auto vec = col2.to_vector();
    CHECK(vec(0) == 3);
}

TEST_CASE("Views with complex numbers") {
    using Cplx = damp::complex<double>;
    Matrix<2, 2, Cplx> mat = {
        {Cplx{1, 1}, Cplx{2, 2}},
        {Cplx{3, 3}, Cplx{4, 4}}
    };

    auto diag = mat.diagonal();
    CHECK(diag(0).real() == 1);
    CHECK(diag(0).imag() == 1);
    CHECK(diag(1).real() == 4);
    CHECK(diag(1).imag() == 4);

    auto upper = mat.upper_triangle();
    CHECK(upper(0, 0).real() == 1);
    CHECK(upper(0, 0).imag() == 1);
    CHECK(upper(0, 1).real() == 2);
    CHECK(upper(0, 1).imag() == 2);
    CHECK(upper(1, 1).real() == 4);
    CHECK(upper(1, 1).imag() == 4);

    auto lower = mat.lower_triangle();
    CHECK(lower(0, 0).real() == 1);
    CHECK(lower(0, 0).imag() == 1);
    CHECK(lower(1, 0).real() == 3);
    CHECK(lower(1, 0).imag() == 3);
    CHECK(lower(1, 1).real() == 4);
    CHECK(lower(1, 1).imag() == 4);

    auto row0 = mat.row(0);
    CHECK(row0(0).real() == 1);
    CHECK(row0(1).real() == 2);

    auto col0 = mat.col(0);
    CHECK(col0(0).real() == 1);
    CHECK(col0(1).real() == 3);
}

TEST_CASE("Views constexpr compatibility") {
    static constexpr Matrix<2, 2, double> mat = {
        {1, 2},
        {3, 4}
    };

    constexpr auto diag = mat.diagonal();
    static_assert(diag(0) == 1);
    static_assert(diag(1) == 4);

    constexpr auto upper = mat.upper_triangle();
    static_assert(upper(0, 0) == 1);
    static_assert(upper(0, 1) == 2);
    static_assert(upper(1, 1) == 4);

    constexpr auto lower = mat.lower_triangle();
    static_assert(lower(0, 0) == 1);
    static_assert(lower(1, 0) == 3);
    static_assert(lower(1, 1) == 4);

    constexpr auto row0 = mat.row(0);
    static_assert(row0(0) == 1);
    static_assert(row0(1) == 2);

    constexpr auto col0 = mat.col(0);
    static_assert(col0(0) == 1);
    static_assert(col0(1) == 3);
}

TEST_CASE("RowView arithmetic operations") {
    Matrix<3, 4, double> mat = {
        {1, 2, 3, 4},
        {5, 6, 7, 8},
        {9, 10, 11, 12}
    };

    auto row1 = mat.row(1); // [5, 6, 7, 8]
    auto row2 = mat.row(2); // [9, 10, 11, 12]

    // RowView + RowView
    auto sum = row1 + row2;
    CHECK(sum(0, 0) == 5 + 9);
    CHECK(sum(0, 1) == 6 + 10);
    CHECK(sum(0, 2) == 7 + 11);
    CHECK(sum(0, 3) == 8 + 12);

    // RowView - RowView
    auto diff = row1 - row2;
    CHECK(diff(0, 0) == 5 - 9);
    CHECK(diff(0, 1) == 6 - 10);

    // RowView + scalar
    auto row_plus_scalar = row1 + 10.0;
    CHECK(row_plus_scalar(0, 0) == 5 + 10);
    CHECK(row_plus_scalar(0, 1) == 6 + 10);

    // RowView * scalar
    auto row_times_scalar = row1 * 2.0;
    CHECK(row_times_scalar(0, 0) == 5 * 2);
    CHECK(row_times_scalar(0, 1) == 6 * 2);

    // RowView + Matrix<1, 4>
    Matrix<1, 4, double> row_mat{{1, 2, 3, 4}};
    auto                 row_plus_mat = row1 + row_mat;
    CHECK(row_plus_mat(0, 0) == 5 + 1);
    CHECK(row_plus_mat(0, 1) == 6 + 2);

    // RowView + RowVec
    RowVec<4, double> row_vec{1, 2, 3, 4};
    auto              row_plus_vec = row1 + row_vec;
    CHECK(row_plus_vec(0, 0) == 5 + 1);
    CHECK(row_plus_vec(0, 1) == 6 + 2);
}

TEST_CASE("ColView arithmetic operations") {
    Matrix<3, 4, double> mat = {
        {1, 2, 3, 4},
        {5, 6, 7, 8},
        {9, 10, 11, 12}
    };

    auto col1 = mat.col(1); // [2, 6, 10]
    auto col2 = mat.col(2); // [3, 7, 11]

    // ColView + ColView
    auto sum = col1 + col2;
    CHECK(sum(0, 0) == 2 + 3);
    CHECK(sum(1, 0) == 6 + 7);
    CHECK(sum(2, 0) == 10 + 11);

    // ColView - ColView
    auto diff = col1 - col2;
    CHECK(diff(0, 0) == 2 - 3);
    CHECK(diff(1, 0) == 6 - 7);

    // ColView + scalar
    auto col_plus_scalar = col1 + 10.0;
    CHECK(col_plus_scalar(0, 0) == 2 + 10);
    CHECK(col_plus_scalar(1, 0) == 6 + 10);

    // ColView * scalar
    auto col_times_scalar = col1 * 2.0;
    CHECK(col_times_scalar(0, 0) == 2 * 2);
    CHECK(col_times_scalar(1, 0) == 6 * 2);

    // ColView + Matrix<3, 1>
    Matrix<3, 1, double> col_mat{{1}, {2}, {3}};
    auto                 col_plus_mat = col1 + col_mat;
    CHECK(col_plus_mat(0, 0) == 2 + 1);
    CHECK(col_plus_mat(1, 0) == 6 + 2);

    // ColView + ColVec
    ColVec<3, double> col_vec{1, 2, 3};
    auto              col_plus_vec = col1 + col_vec;
    CHECK(col_plus_vec(0, 0) == 2 + 1);
    CHECK(col_plus_vec(1, 0) == 6 + 2);
}

TEST_CASE("Cross-type view arithmetic operations") {
    Matrix<3, 4, double> mat = {
        {1, 2, 3, 4},
        {5, 6, 7, 8},
        {9, 10, 11, 12}
    };

    auto row1 = mat.row(1); // [5, 6, 7, 8]
    auto col1 = mat.col(1); // [2, 6, 10]

    // RowView + ColView (broadcasting to 3x4 matrix)
    auto row_plus_col = row1 + col1;
    CHECK(row_plus_col.rows() == 3);
    CHECK(row_plus_col.cols() == 4);
    CHECK(row_plus_col(0, 0) == 5 + 2);
    CHECK(row_plus_col(0, 1) == 6 + 2);
    CHECK(row_plus_col(0, 2) == 7 + 2);
    CHECK(row_plus_col(0, 3) == 8 + 2);
    CHECK(row_plus_col(1, 0) == 5 + 6);
    CHECK(row_plus_col(1, 1) == 6 + 6);
    CHECK(row_plus_col(2, 0) == 5 + 10);

    // Block + RowView
    auto block = mat.block<2, 4>(1, 0); // rows 1-2, cols 0-3
    auto block_plus_row = block + row1;
    CHECK(block_plus_row.rows() == 2);
    CHECK(block_plus_row.cols() == 4);
    CHECK(block_plus_row(0, 0) == 5 + 5);
    CHECK(block_plus_row(0, 1) == 6 + 6);
    CHECK(block_plus_row(1, 0) == 9 + 5);
    CHECK(block_plus_row(1, 1) == 10 + 6);
}

TEST_CASE("TransposeView basic functionality") {
    Matrix<2, 3, double> mat = {
        {1, 2, 3},
        {4, 5, 6}
    };

    auto tv = mat.t();

    // Dimensions are swapped
    CHECK(tv.rows() == 3);
    CHECK(tv.cols() == 2);

    // Element access: tv(r,c) == mat(c,r)
    CHECK(tv(0, 0) == 1);
    CHECK(tv(0, 1) == 4);
    CHECK(tv(1, 0) == 2);
    CHECK(tv(1, 1) == 5);
    CHECK(tv(2, 0) == 3);
    CHECK(tv(2, 1) == 6);

    // to_matrix produces an owning copy
    auto copied = tv.to_matrix();
    CHECK(copied(0, 0) == 1);
    CHECK(copied(0, 1) == 4);
    CHECK(copied(2, 1) == 6);

    // Equivalent to .transpose()
    auto transposed = mat.transpose();
    for (size_t r = 0; r < 3; ++r) {
        for (size_t c = 0; c < 2; ++c) {
            CHECK(tv(r, c) == transposed(r, c));
        }
    }
}

TEST_CASE("TransposeView modification through view") {
    Matrix<2, 3, double> mat = {
        {1, 2, 3},
        {4, 5, 6}
    };

    auto tv = mat.t();

    // Write through the view
    tv(0, 1) = 99.0; // should modify mat(1, 0)
    CHECK(mat(1, 0) == 99.0);

    tv(2, 0) = -7.0; // should modify mat(0, 2)
    CHECK(mat(0, 2) == -7.0);
}

TEST_CASE("TransposeView const correctness") {
    const Matrix<2, 3, double> mat = {
        {1, 2, 3},
        {4, 5, 6}
    };

    auto tv = mat.t();

    // Should be readable
    CHECK(tv(0, 0) == 1);
    CHECK(tv(1, 1) == 5);

    // const view: tv(r,c) = X should not compile
    static_assert(std::is_const_v<std::remove_reference_t<decltype(tv(0, 0))>>);
}

TEST_CASE("TransposeView arithmetic via MatrixLike") {
    Matrix<2, 3, double> A = {
        {1, 2, 3},
        {4, 5, 6}
    };

    Matrix<3, 2, double> B = {
        {10, 20},
        {30, 40},
        {50, 60}
    };

    // TransposeView + Matrix (element-wise)
    auto sum = A.t() + B;
    CHECK(sum(0, 0) == 1 + 10);
    CHECK(sum(0, 1) == 4 + 20);
    CHECK(sum(1, 0) == 2 + 30);
    CHECK(sum(2, 1) == 6 + 60);

    // TransposeView - Matrix (element-wise)
    auto diff = A.t() - B;
    CHECK(diff(0, 0) == 1 - 10);
    CHECK(diff(2, 1) == 6 - 60);

    // Matrix * TransposeView (matmul)
    // A is 2x3, A.t() is 3x2 → A * A.t() is 2x2
    auto product = A * A.t();
    // A * A^T = [1*1+2*2+3*3,  1*4+2*5+3*6]
    //           [4*1+5*2+6*3,  4*4+5*5+6*6]
    CHECK(product(0, 0) == doctest::Approx(14));
    CHECK(product(0, 1) == doctest::Approx(32));
    CHECK(product(1, 0) == doctest::Approx(32));
    CHECK(product(1, 1) == doctest::Approx(77));

    // Unary minus
    auto neg = -A.t();
    CHECK(neg(0, 0) == -1);
    CHECK(neg(2, 1) == -6);
}

TEST_CASE("TransposeView square matrix (A.t() * P * A pattern)") {
    // Classic Kalman-style pattern: A^T * P * A
    Matrix<3, 3, double> A = {
        {1, 2, 0},
        {0, 1, 3},
        {4, 0, 1}
    };

    Matrix<3, 3, double> P = Matrix<3, 3, double>::identity();

    // A.t() * P * A  — uses TransposeView for A.t(), no copy until multiplication
    auto result = A.t() * P * A;
    // With P = I, this is just A^T * A
    // A^T * A = [1  0  4] [1  2  0]   [17  2  4 ]
    //           [2  1  0] [0  1  3] = [2   5  3 ]
    //           [0  3  1] [4  0  1]   [4   3  10]
    CHECK(result(0, 0) == doctest::Approx(17));
    CHECK(result(0, 1) == doctest::Approx(2));
    CHECK(result(0, 2) == doctest::Approx(4));
    CHECK(result(1, 0) == doctest::Approx(2));
    CHECK(result(1, 1) == doctest::Approx(5));
    CHECK(result(1, 2) == doctest::Approx(3));
    CHECK(result(2, 0) == doctest::Approx(4));
    CHECK(result(2, 1) == doctest::Approx(3));
    CHECK(result(2, 2) == doctest::Approx(10));
}

TEST_CASE("TransposeView constexpr compatibility") {
    static constexpr Matrix<2, 3, double> mat = {
        {1, 2, 3},
        {4, 5, 6}
    };

    static constexpr auto tv = mat.t();
    static_assert(tv.rows() == 3);
    static_assert(tv.cols() == 2);
    static_assert(tv(0, 0) == 1);
    static_assert(tv(0, 1) == 4);
    static_assert(tv(2, 0) == 3);
    static_assert(tv(2, 1) == 6);

    // to_matrix at compile time
    constexpr auto copied = tv.to_matrix();
    static_assert(copied(0, 0) == 1);
    static_assert(copied(2, 1) == 6);
}

TEST_CASE("TransposeView assignment") {
    Matrix<2, 3, double> mat = {
        {1, 2, 3},
        {4, 5, 6}
    };

    Matrix<3, 2, double> replacement = {
        {10, 40},
        {20, 50},
        {30, 60}
    };

    auto tv = mat.t();
    tv = replacement;

    // mat should now be the transpose of replacement
    CHECK(mat(0, 0) == 10);
    CHECK(mat(0, 1) == 20);
    CHECK(mat(0, 2) == 30);
    CHECK(mat(1, 0) == 40);
    CHECK(mat(1, 1) == 50);
    CHECK(mat(1, 2) == 60);
}