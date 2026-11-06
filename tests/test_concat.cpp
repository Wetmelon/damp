// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>
#include <type_traits>

#include "damp/control.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

TEST_SUITE("diag / blkdiag / hcat / vcat") {

    TEST_CASE("constexpr diag scalars are the north-star Q") {
        constexpr auto Q = diag(10.0, 0.1, 100.0, 10.0);
        static_assert(Q.rows() == 4);
        static_assert(Q.cols() == 4);
        static_assert(std::is_same_v<decltype(Q)::value_type, double>);
        static_assert(Q(0, 0) == 10.0);
        static_assert(Q(1, 1) == 0.1);
        static_assert(Q(2, 2) == 100.0);
        static_assert(Q(3, 3) == 10.0);
        static_assert(Q(0, 1) == 0.0);
        static_assert(Q(1, 0) == 0.0);
        static_assert(Q(2, 3) == 0.0);

        CHECK(Q(0, 0) == 10.0);
        CHECK(Q(1, 1) == 0.1);
        CHECK(Q(2, 2) == 100.0);
        CHECK(Q(3, 3) == 10.0);
        for (size_t r = 0; r < Q.rows(); ++r) {
            for (size_t c = 0; c < Q.cols(); ++c) {
                if (r != c) {
                    CHECK(Q(r, c) == 0.0);
                }
            }
        }
    }

    TEST_CASE("diag of ints promotes to double") {
        constexpr auto Q = diag(10, 1, 100);
        static_assert(Q.rows() == 3);
        static_assert(Q.cols() == 3);
        static_assert(std::is_same_v<decltype(Q)::value_type, double>);
        static_assert(Q(0, 0) == 10.0);
        static_assert(Q(1, 1) == 1.0);
        static_assert(Q(2, 2) == 100.0);
        CHECK(Q(0, 1) == 0.0);
        CHECK(Q(2, 0) == 0.0);
    }

    TEST_CASE("1x1 diag from a scalar") {
        constexpr auto s = diag(3.0);
        static_assert(s.rows() == 1);
        static_assert(s.cols() == 1);
        static_assert(std::is_same_v<decltype(diag(3.0)), Matrix<1, 1, double>>);
        static_assert(s(0, 0) == 3.0);
        CHECK(s(0, 0) == 3.0);
    }

    TEST_CASE("extract diagonal of a square matrix") {
        constexpr auto d = diag(Matrix<2, 2>{{1.0, 2.0}, {3.0, 4.0}});
        static_assert(std::is_same_v<decltype(diag(Matrix<2, 2>{{1.0, 2.0}, {3.0, 4.0}})), ColVec<2, double>>);
        static_assert(d.rows() == 2);
        static_assert(d.cols() == 1);
        CHECK(d[0] == 1.0);
        CHECK(d[1] == 4.0);
    }

    TEST_CASE("diag from ColVec constructs a square matrix") {
        constexpr ColVec<3> v{1.0, 2.0, 3.0};
        constexpr auto      D = diag(v);
        static_assert(D.rows() == 3);
        static_assert(D.cols() == 3);
        static_assert(std::is_same_v<decltype(D), const Matrix<3, 3, double>>);
        CHECK(D(0, 0) == 1.0);
        CHECK(D(1, 1) == 2.0);
        CHECK(D(2, 2) == 3.0);
        CHECK(D(0, 1) == 0.0);
        CHECK(D(1, 0) == 0.0);
        CHECK(D(2, 0) == 0.0);
    }

    TEST_CASE("1x1 ColVec uses extract, not construct") {
        constexpr ColVec<1> c{3.0};
        constexpr auto      d = diag(c);
        static_assert(std::is_same_v<decltype(diag(ColVec<1>{3.0})), ColVec<1, double>>);
        static_assert(d.rows() == 1);
        static_assert(d.cols() == 1);
        CHECK(d[0] == 3.0);
    }

    TEST_CASE("blkdiag of 2x2, 3x3, 1x1") {
        constexpr Matrix<2, 2> mat1 = {{1, 2}, {3, 4}};
        constexpr Matrix<3, 3> mat2 = {{5, 6, 7}, {8, 9, 10}, {11, 12, 13}};
        constexpr Matrix<1, 1> mat3 = {{14}};

        constexpr auto block_diag = blkdiag(mat1, mat2, mat3);
        static_assert(block_diag.rows() == 6);
        static_assert(block_diag.cols() == 6);

        CHECK(block_diag(0, 0) == 1);
        CHECK(block_diag(0, 1) == 2);
        CHECK(block_diag(1, 0) == 3);
        CHECK(block_diag(1, 1) == 4);

        CHECK(block_diag(2, 2) == 5);
        CHECK(block_diag(2, 3) == 6);
        CHECK(block_diag(2, 4) == 7);
        CHECK(block_diag(3, 2) == 8);
        CHECK(block_diag(3, 3) == 9);
        CHECK(block_diag(3, 4) == 10);
        CHECK(block_diag(4, 2) == 11);
        CHECK(block_diag(4, 3) == 12);
        CHECK(block_diag(4, 4) == 13);

        CHECK(block_diag(5, 5) == 14);

        for (size_t r = 0; r < block_diag.rows(); ++r) {
            for (size_t c = 0; c < block_diag.cols(); ++c) {
                if ((r < 2 && c < 2) || (r >= 2 && r < 5 && c >= 2 && c < 5) || (r == 5 && c == 5)) {
                    continue;
                }
                CHECK(block_diag(r, c) == 0);
            }
        }
    }

    TEST_CASE("blkdiag of Matrix and ColVec") {
        constexpr Matrix<2, 2> A{{1.0, 2.0}, {3.0, 4.0}};
        constexpr ColVec<3>    v{5.0, 6.0, 7.0};
        constexpr auto         B = blkdiag(A, v);
        static_assert(B.rows() == 5);
        static_assert(B.cols() == 3);
        CHECK(B(0, 0) == 1.0);
        CHECK(B(0, 1) == 2.0);
        CHECK(B(1, 0) == 3.0);
        CHECK(B(1, 1) == 4.0);
        CHECK(B(2, 2) == 5.0);
        CHECK(B(3, 2) == 6.0);
        CHECK(B(4, 2) == 7.0);
        CHECK(B(0, 2) == 0.0);
        CHECK(B(2, 0) == 0.0);
        CHECK(B(4, 0) == 0.0);
    }

    TEST_CASE("hcat two 2x1 is 2x2") {
        constexpr Matrix<2, 1> a{{1.0}, {2.0}};
        constexpr Matrix<2, 1> b{{3.0}, {4.0}};
        constexpr auto         H = hcat(a, b);
        static_assert(H.rows() == 2);
        static_assert(H.cols() == 2);
        CHECK(H(0, 0) == 1.0);
        CHECK(H(1, 0) == 2.0);
        CHECK(H(0, 1) == 3.0);
        CHECK(H(1, 1) == 4.0);
    }

    TEST_CASE("vcat two 1x2 is 2x2") {
        constexpr Matrix<1, 2> a{{1.0, 2.0}};
        constexpr Matrix<1, 2> b{{3.0, 4.0}};
        constexpr auto         V = vcat(a, b);
        static_assert(V.rows() == 2);
        static_assert(V.cols() == 2);
        CHECK(V(0, 0) == 1.0);
        CHECK(V(0, 1) == 2.0);
        CHECK(V(1, 0) == 3.0);
        CHECK(V(1, 1) == 4.0);
    }

    TEST_CASE("mixed ColVec / Matrix hcat") {
        constexpr ColVec<2>    c{1.0, 2.0};
        constexpr Matrix<2, 2> M{{3.0, 4.0}, {5.0, 6.0}};
        constexpr auto         H = hcat(c, M);
        static_assert(H.rows() == 2);
        static_assert(H.cols() == 3);
        CHECK(H(0, 0) == 1.0);
        CHECK(H(1, 0) == 2.0);
        CHECK(H(0, 1) == 3.0);
        CHECK(H(0, 2) == 4.0);
        CHECK(H(1, 1) == 5.0);
        CHECK(H(1, 2) == 6.0);
    }
}
