// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Oracle conformance: fixed-seed property tests (checks in oracle_checks.hpp, cases in oracle_cases.hpp).

#include "oracle_checks.hpp"

using namespace oracle_checks;

TEST_SUITE("Oracle conformance") {
    TEST_CASE("property: QR and SVD reconstruct random matrices with orthogonal factors") {
        auto shapes = []<typename T>() {
            qr_property<2, 2, T>();
            qr_property<3, 3, T>();
            qr_property<4, 3, T>();
            qr_property<5, 5, T>();
            qr_property<6, 4, T>();
            qr_property<6, 6, T>();
            svd_property<2, 2, T>();
            svd_property<3, 5, T>();
            svd_property<4, 4, T>();
            svd_property<5, 3, T>();
            svd_property<6, 2, T>();
            svd_property<6, 6, T>();
        };
        shapes.template operator()<double>();
        shapes.template operator()<float>();
    }
    TEST_CASE("property: LU and Cholesky solves are backward stable on random systems") {
        auto sizes = []<typename T>() {
            solve_property<2, T>();
            solve_property<3, T>();
            solve_property<4, T>();
            solve_property<5, T>();
            solve_property<6, T>();
        };
        sizes.template operator()<double>();
        sizes.template operator()<float>();
    }
    TEST_CASE("property: lyap and dlyap residuals on random stable systems") {
        auto sizes = []<typename T>() {
            lyapunov_property<2, T>();
            lyapunov_property<3, T>();
            lyapunov_property<4, T>();
            lyapunov_property<5, T>();
            lyapunov_property<6, T>();
        };
        sizes.template operator()<double>();
        sizes.template operator()<float>();
    }
    TEST_CASE("property: ZOH resampled through the logarithm equals direct discretization") {
        auto sizes = []<typename T>() {
            resample_property<2, 1, T>();
            resample_property<3, 2, T>();
            resample_property<4, 1, T>();
            resample_property<5, 2, T>();
        };
        sizes.template operator()<double>();
        sizes.template operator()<float>();
    }
}
