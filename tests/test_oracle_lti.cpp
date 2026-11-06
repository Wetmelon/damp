// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Oracle conformance: lyap / dlyap, discretization, the SISO conversions, invariant zeros and place
// (checks in oracle_checks.hpp, cases in oracle_cases.hpp).

#include "oracle_checks.hpp"

using namespace oracle_checks;

TEST_SUITE("Oracle conformance") {
    TEST_CASE("lyap and dlyap certify against the Lyapunov equation and match scipy") {
        run_cases<oracle::lyap_cases>([]<typename C, typename T>() { check_lyapunov<false, C, T>(); });
        run_cases<oracle::dlyap_cases>([]<typename C, typename T>() { check_lyapunov<true, C, T>(); });
        run_cases<oracle::lyap_fail_cases>([]<typename C, typename T>() { check_lyapunov_fails<false, C, T>(); });
        run_cases<oracle::dlyap_fail_cases>([]<typename C, typename T>() { check_lyapunov_fails<true, C, T>(); });
    }
    TEST_CASE("discretize ZOH and Tustin match scipy") {
        run_cases<oracle::zoh_cases>([]<typename C, typename T>() { check_zoh<C, T>(); });
        run_cases<oracle::tustin_cases>([]<typename C, typename T>() { check_tustin<C, T>(); });
        run_cases<oracle::tustin_fail_cases>([]<typename C, typename T>() { check_tustin_fails<C, T>(); });
    }
    TEST_CASE("ss2tf and ss2zpk match exact Leverrier coefficients and their roots") {
        run_cases<oracle::ss2tf_cases>([]<typename C, typename T>() {
            check_ss2tf<C, T>();
            check_ss2zpk<C, T>();
        });
    }
    TEST_CASE("tf2zpk and zpk2tf match the factored and expanded forms") {
        run_cases<oracle::tf2zpk_cases>([]<typename C, typename T>() { check_tf2zpk<C, T>(); });
        run_cases<oracle::zpk2tf_cases>([]<typename C, typename T>() { check_zpk2tf<C, T>(); });
    }
    TEST_CASE("transmission_zeros matches the exact SISO zeros and the MIMO zeros of the system pencil") {
        run_cases<oracle::ss2tf_cases>([]<typename C, typename T>() { check_ss2tf_tzero<C, T>(); });
        run_cases<oracle::tzero_cases>([]<typename C, typename T>() { check_tzero<C, T>(); });
    }
    TEST_CASE("tf2ss realizes the transfer function's Markov parameters") {
        run_cases<oracle::tf2ss_cases>([]<typename C, typename T>() { check_tf2ss<C, T>(); });
    }
    TEST_CASE("place puts the closed-loop eigenvalues on the requested poles") {
        run_cases<oracle::place_cases>([]<typename C, typename T>() { check_place<C, T>(); });
    }
}
