// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Oracle conformance: eigenvalues, the Schur basis, generalized eigenvalues (QZ), logm and sqrtm
// (checks in oracle_checks.hpp, cases in oracle_cases.hpp).

#include "oracle_checks.hpp"

using namespace oracle_checks;

TEST_SUITE("Oracle conformance") {
    TEST_CASE("eigenvalues match scipy within their condition numbers") {
        run_cases<oracle::eig_cases>([]<typename C, typename T>() { check_eig<C, T>(); });
    }
    TEST_CASE("the Schur basis is orthogonal and triangularizes the balanced matrix") {
        run_cases<oracle::eig_cases>([]<typename C, typename T>() { check_schur<C, T>(); });
    }
    TEST_CASE("generalized eigenvalues match scipy within their condition numbers and count the infinite ones") {
        run_cases<oracle::qz_cases>([]<typename C, typename T>() { check_qz<C, T>(); });
    }
    TEST_CASE("logm matches scipy within its condition number") {
        run_cases<oracle::logm_cases>([]<typename C, typename T>() { check_logm<C, T>(); });
    }
    TEST_CASE("sqrtm squares back to A and matches scipy") {
        run_cases<oracle::sqrtm_cases>([]<typename C, typename T>() { check_sqrtm<C, T>(); });
    }
    TEST_CASE("sqrtm and logm balance a graded matrix, not one of a single scale") {
        // Balancing sqrtm_jordan_3 (κ(D) = 16 for a 1.16× smaller norm) cost up to 9× in accuracy, which
        // only some instruction schedules showed: check the decision itself, on every host.
        using J = oracle::sqrtm_jordan_3;
        using G = oracle::sqrtm_graded_4;
        const auto A = load<J::n, J::n, double>(J::A);
        const auto B = load<G::n, G::n, double>(G::A);
        CHECK(mat::detail::balance_exponents(A) != damp::array<int, J::n>{});
        CHECK(mat::detail::function_balance_exponents(A) == damp::array<int, J::n>{});
        CHECK(mat::detail::function_balance_exponents(B) == mat::detail::balance_exponents(B));
        CHECK(mat::detail::balance_exponents(B) != damp::array<int, G::n>{});
    }
}
