// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Oracle conformance: expm, CARE and DARE (checks in oracle_checks.hpp, cases in oracle_cases.hpp).

#include "oracle_checks.hpp"

using namespace oracle_checks;

TEST_SUITE("Oracle conformance") {
    TEST_CASE("expm matches scipy within its condition number") {
        for_each_case<oracle::expm_cases>([]<typename C>() {
            check_expm<C, double>();
            if constexpr (C::float_ok) {
                check_expm<C, float>();
            }
        });
    }
    TEST_CASE("care solutions certify against the Riccati equation") {
        check_riccati_cases<false, oracle::care_cases, false>();
    }
    TEST_CASE("dare solutions certify against the Riccati equation") {
        check_riccati_cases<true, oracle::dare_cases, false>();
    }
    TEST_CASE("care refuses a problem with no stabilizing solution") {
        check_riccati_fail_cases<false, oracle::care_fail_cases, false>();
    }
    TEST_CASE("dare refuses a problem with no stabilizing solution") {
        check_riccati_fail_cases<true, oracle::dare_fail_cases, false>();
    }
    TEST_CASE("known library bugs found by the oracle suite" * doctest::may_fail()) {
        check_riccati_cases<false, oracle::care_cases, true>();
        check_riccati_cases<true, oracle::dare_cases, true>();
        check_riccati_fail_cases<false, oracle::care_fail_cases, true>();
        check_riccati_fail_cases<true, oracle::dare_fail_cases, true>();
    }
    TEST_CASE("expm at compile time agrees with run time") {
        check_expm_constexpr<oracle::expm_random_3, double>();
        check_expm_constexpr<oracle::expm_random_3, float>();
        check_expm_constexpr<oracle::expm_nonnormal_rotated_2, double>();
        check_expm_constexpr<oracle::expm_stiff_3, double>();
        check_expm_constexpr<oracle::expm_stiff_3, float>();
    }
    TEST_CASE("care at compile time agrees with run time") {
        check_riccati_constexpr<false, oracle::care_random_3x1, double>();
        check_riccati_constexpr<false, oracle::care_random_3x1, float>();
        check_riccati_constexpr<false, oracle::care_lightly_damped_2, double>();
        static_assert(!run_riccati<false, oracle::care_fail_unstabilizable_2, double>(load<2, 2, double>(oracle::care_fail_unstabilizable_2::A)).has_value());
    }
    TEST_CASE("dare at compile time agrees with run time") {
        check_riccati_constexpr<true, oracle::dare_random_3x1, double>();
        check_riccati_constexpr<true, oracle::dare_random_3x1, float>();
        check_riccati_constexpr<true, oracle::dare_lightly_damped_2, double>();
        check_riccati_constexpr<true, oracle::dare_cross_term_3, double>();
    }
}
