// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Oracle conformance: LU, QR, Cholesky, LDLᵀ and SVD (checks in oracle_checks.hpp, cases in oracle_cases.hpp).

#include "oracle_checks.hpp"

using namespace oracle_checks;

TEST_SUITE("Oracle conformance") {
    TEST_CASE("lu factors, determinant and solve against exact arithmetic") {
        run_cases<oracle::lu_cases>([]<typename C, typename T>() { check_lu<C, T>(); });
    }
    TEST_CASE("qr is orthogonal, reconstructs A and matches scipy's R") {
        run_cases<oracle::qr_cases>([]<typename C, typename T>() { check_qr<C, T>(); });
    }
    TEST_CASE("cholesky factor and solve against scipy and exact arithmetic") {
        run_cases<oracle::chol_cases>([]<typename C, typename T>() { check_cholesky<C, T>(); });
        run_cases<oracle::chol_fail_cases>([]<typename C, typename T>() { check_cholesky_fails<C, T>(); });
    }
    TEST_CASE("ldlt factors, inertia and solve against exact arithmetic") {
        run_cases<oracle::ldlt_cases>([]<typename C, typename T>() { check_ldlt<C, T>(); });
        run_cases<oracle::ldlt_fail_cases>([]<typename C, typename T>() { check_ldlt_fails<C, T>(); });
    }
    TEST_CASE("svd reconstructs A and matches scipy's singular values") {
        run_cases<oracle::svd_cases>([]<typename C, typename T>() { check_svd<C, T>(); });
    }
}
