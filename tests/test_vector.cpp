// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <type_traits>

#include "damp/matrix/matrix.hpp"


#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

TEST_SUITE("ColVec") {
    TEST_CASE("ColVec basic construction and access") {
        ColVec<3, float> vec;
        // Default initialized to 0
        CHECK(vec[0] == 0.0f);
        CHECK(vec[1] == 0.0f);
        CHECK(vec[2] == 0.0f);

        // Set values
        vec[0] = 1.0f;
        vec[1] = 2.0f;
        vec[2] = 3.0f;

        CHECK(vec[0] == 1.0f);
        CHECK(vec[1] == 2.0f);
        CHECK(vec[2] == 3.0f);
    }

    TEST_CASE("ColVec C array reference constructor") {
        const float      abc[3] = {1.0f, 2.0f, 3.0f};
        ColVec<3, float> v(abc);
        CHECK(v[0] == doctest::Approx(1.0f));
        CHECK(v[1] == doctest::Approx(2.0f));
        CHECK(v[2] == doctest::Approx(3.0f));

        // CTAD from T[N]
        ColVec w = abc;
        static_assert(std::is_same_v<decltype(w), ColVec<3, float>>);
        CHECK(w[0] == doctest::Approx(1.0f));

        // Implicit conversion at call boundary
        auto first = [](const ColVec<3, float>& x) { return x[0]; };
        CHECK(first(abc) == doctest::Approx(1.0f));
    }

    TEST_CASE("ColVec initializer list constructor") {
        ColVec<3> vec = {1, 2, 3};

        CHECK(vec[0] == 1);
        CHECK(vec[1] == 2);
        CHECK(vec[2] == 3);
    }

    TEST_CASE("ColVec copy and assignment") {
        ColVec vec1 = {1.0f, 2.0f, 3.0f};

        ColVec vec2 = vec1;
        CHECK(vec2[0] == 1.0f);
        CHECK(vec2[2] == 3.0f);

        ColVec<3, float> vec3;
        vec3 = vec1;
        CHECK(vec3[1] == 2.0f);
    }

    TEST_CASE("ColVec addition and subtraction") {
        ColVec vec1 = {1, 2, 3};
        ColVec vec2 = {4, 5, 6};

        auto vec_sum = vec1 + vec2;
        CHECK(vec_sum[0] == 5);
        CHECK(vec_sum[1] == 7);
        CHECK(vec_sum[2] == 9);

        auto vec_diff = vec2 - vec1;
        CHECK(vec_diff[0] == 3);
        CHECK(vec_diff[1] == 3);
        CHECK(vec_diff[2] == 3);

        vec1 += vec2;
        CHECK(vec1[0] == 5);

        vec1 -= vec2;
        CHECK(vec1[0] == 1);
    }

    TEST_CASE("ColVec scalar operations") {
        ColVec<3> vec = {1, 2, 3};

        auto vec_scaled = vec * 2;
        CHECK(vec_scaled[0] == 2);
        CHECK(vec_scaled[1] == 4);
        CHECK(vec_scaled[2] == 6);

        vec_scaled = 3 * vec;
        CHECK(vec_scaled[0] == 3);
        CHECK(vec_scaled[1] == 6);
        CHECK(vec_scaled[2] == 9);

        vec_scaled = vec / 2;
        CHECK(vec_scaled[0] == doctest::Approx(0.5));
        CHECK(vec_scaled[1] == doctest::Approx(1.0));
        CHECK(vec_scaled[2] == doctest::Approx(1.5));

        vec *= 2;
        CHECK(vec[0] == 2);

        vec /= 2;
        CHECK(vec[0] == 1);
    }

    TEST_CASE("ColVec equality") {
        ColVec<3> vec1 = {1, 2, 3};
        ColVec<3> vec2 = {1, 2, 3};
        ColVec<3> vec3 = {1, 2, 4};

        CHECK(vec1 == vec2);
        CHECK(vec1 != vec3);
    }

    TEST_CASE("ColVec unary negation") {
        ColVec<3> vec = {1, -2, 3};

        auto vec_neg = -vec;
        CHECK(vec_neg[0] == -1);
        CHECK(vec_neg[1] == 2);
        CHECK(vec_neg[2] == -3);
    }

    TEST_CASE("ColVec dot product") {
        ColVec<3> vec1 = {1, 2, 3};
        ColVec<3> vec2 = {4, 5, 6};

        auto dot_prod = dot(vec1, vec2);
        CHECK(dot_prod == 32); // 1*4 + 2*5 + 3*6
    }

    TEST_CASE("ColVec cross product") {
        ColVec vec1 = {1, 2, 3};
        ColVec vec2 = {4, 5, 6};

        auto vec_cross = vec1.cross(vec2);
        CHECK(vec_cross[0] == -3); // 2*6 - 3*5
        CHECK(vec_cross[1] == 6);  // 3*4 - 1*6
        CHECK(vec_cross[2] == -3); // 1*5 - 2*4
    }

    TEST_CASE("ColVec norm and normalized") {
        ColVec vec = {3.0f, 4.0f, 0.0f};

        auto n = vec.norm();
        CHECK(n == 5.0f); // sqrt(9 + 16 + 0)

        auto vec_norm = vec.normalized();
        CHECK(vec_norm[0] == 0.6f); // 3/5
        CHECK(vec_norm[1] == 0.8f); // 4/5
        CHECK(vec_norm[2] == 0.0f);
    }

    TEST_CASE("ColVec constexpr") {
        constexpr ColVec<3> vec = {1, 2, 3};

        static_assert(vec[0] == 1);
        static_assert(vec[2] == 3);
    }

    TEST_CASE("ColVec segment view") {
        ColVec<6> x = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};

        // Read middle triad
        auto mid = x.template segment<3>(2);
        CHECK(mid(0, 0) == doctest::Approx(2.0));
        CHECK(mid(1, 0) == doctest::Approx(3.0));
        CHECK(mid(2, 0) == doctest::Approx(4.0));

        // Write-back into parent
        mid = ColVec<3>{10.0, 20.0, 30.0};
        CHECK(x[2] == doctest::Approx(10.0));
        CHECK(x[3] == doctest::Approx(20.0));
        CHECK(x[4] == doctest::Approx(30.0));
        CHECK(x[0] == doctest::Approx(0.0));
        CHECK(x[5] == doctest::Approx(5.0));

        // Owning copy via to_vector
        const ColVec<6> y = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        const ColVec<3> head3 = y.template segment<3>(0).to_vector();
        CHECK(head3[0] == doctest::Approx(1.0));
        CHECK(head3[1] == doctest::Approx(2.0));
        CHECK(head3[2] == doctest::Approx(3.0));

        // Mat3 * segment → usable as MatrixLike product
        const Matrix<3, 3> I = Matrix<3, 3>::identity();
        const auto         prod = I * y.template segment<3>(3).to_vector();
        CHECK(prod[0] == doctest::Approx(4.0));
        CHECK(prod[1] == doctest::Approx(5.0));
        CHECK(prod[2] == doctest::Approx(6.0));

        // segment itself is constexpr-callable (Block is a pointer view, so the
        // returned object is not a usable constant-expression value).
        constexpr auto sum_mid = []() {
            ColVec<5>  z = {1.0, 2.0, 3.0, 4.0, 5.0};
            const auto s = z.template segment<2>(1);
            return s(0, 0) + s(1, 0);
        }();
        CHECK(sum_mid == doctest::Approx(5.0));
    }
}

TEST_SUITE("RowVec") {
    TEST_CASE("RowVec basic construction and access") {
        RowVec<3, float> vec;
        CHECK(vec[0] == 0.0f);
        CHECK(vec[1] == 0.0f);
        CHECK(vec[2] == 0.0f);

        vec[0] = 1.0f;
        vec[1] = 2.0f;
        vec[2] = 3.0f;

        CHECK(vec[0] == 1.0f);
        CHECK(vec[1] == 2.0f);
        CHECK(vec[2] == 3.0f);
    }

    TEST_CASE("RowVec initializer list constructor") {
        RowVec<3> vec = {1, 2, 3};

        CHECK(vec[0] == 1);
        CHECK(vec[1] == 2);
        CHECK(vec[2] == 3);
    }

    TEST_CASE("RowVec dot product") {
        RowVec<3> vec1 = {1, 2, 3};
        RowVec<3> vec2 = {4, 5, 6};

        auto dot_prod = dot(vec1, vec2);
        CHECK(dot_prod == 32);
    }

    TEST_CASE("RowVec cross product") {
        RowVec<3> vec1 = {1, 2, 3};
        RowVec<3> vec2 = {4, 5, 6};

        auto vec_cross = vec1.cross(vec2);
        CHECK(vec_cross[0] == -3);
        CHECK(vec_cross[1] == 6);
        CHECK(vec_cross[2] == -3);
    }

    TEST_CASE("RowVec norm and normalized") {
        RowVec vec = {3.0f, 4.0f, 0.0f};

        auto n = vec.norm();
        CHECK(n == 5.0f);

        auto vec_norm = vec.normalized();
        CHECK(vec_norm[0] == 0.6f);
        CHECK(vec_norm[1] == 0.8f);
        CHECK(vec_norm[2] == 0.0f);
    }

    TEST_CASE("RowVec segment view") {
        RowVec<5> r = {0.0, 1.0, 2.0, 3.0, 4.0};

        auto mid = r.template segment<3>(1);
        CHECK(mid(0, 0) == doctest::Approx(1.0));
        CHECK(mid(0, 1) == doctest::Approx(2.0));
        CHECK(mid(0, 2) == doctest::Approx(3.0));

        mid = RowVec<3>{7.0, 8.0, 9.0};
        CHECK(r[0] == doctest::Approx(0.0));
        CHECK(r[1] == doctest::Approx(7.0));
        CHECK(r[2] == doctest::Approx(8.0));
        CHECK(r[3] == doctest::Approx(9.0));
        CHECK(r[4] == doctest::Approx(4.0));
    }
}
