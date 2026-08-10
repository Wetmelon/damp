// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>

#include "damp/backend.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

/**
 * @brief Tests for discretization methods (ZOH, Tustin)
 *
 * These tests verify the accuracy and correctness of continuous-to-discrete
 * system transformations.
 */

TEST_SUITE("Control Design: Discretization") {
    // Test 1: Simple 1st-order system discretization via ZOH
    // System: dx/dt = -x + u (time constant = 1 second)
    // Analytical solution: x[k+1] = exp(-Ts)*x[k] + (1-exp(-Ts))*u[k]
    TEST_CASE("ZOH Discretization: 1st-order system") {
        double     Ts = 0.1; // 100ms sampling time
        StateSpace sys{
            Matrix<1, 1>{{-1.0}}, // A: dx/dt = -x
            Matrix<1, 1>{{1.0}},  // B: input gain
            Matrix<1, 1>{{1.0}}   // C: output is state
        };

        // Use new unified discretize() function
        StateSpace sys_d = *discretize(sys, Ts, DiscretizationMethod::ZOH);

        // Verify A_d ≈ exp(-Ts)
        double A_d_expected = std::exp(-Ts);
        double A_d_actual = sys_d.A(0, 0);
        CHECK(doctest::Approx(A_d_actual).epsilon(1e-4) == A_d_expected);

        // Verify B_d ≈ 1 - exp(-Ts)
        double B_d_expected = 1.0 - std::exp(-Ts);
        double B_d_actual = sys_d.B(0, 0);
        CHECK(doctest::Approx(B_d_actual).epsilon(1e-4) == B_d_expected);
    }

    // Test 2: Tustin Discretization
    TEST_CASE("Tustin Discretization: 1st-order system") {
        double     Ts = 0.1;
        StateSpace sys{
            Matrix<1, 1>{{-1.0}},
            Matrix<1, 1>{{1.0}},
            Matrix<1, 1>{{1.0}}
        };

        // Use new unified discretize() function
        StateSpace sys_d = *discretize(sys, Ts, DiscretizationMethod::Tustin);

        // For Tustin with A=-1, Ts=0.1:
        // M = (I - A*Ts/2)^{-1} = (1 + 0.05)^{-1} = 1/1.05
        // A_d = M * (I + A*Ts/2) = (1/1.05) * 0.95 = (2-Ts)/(2+Ts)
        // B_d = Ts * M * B = 0.1/1.05
        // C_d = C * M = 1/1.05
        // D_d = D + (Ts/2) * C_d * B = 0 + 0.05/1.05
        double expected_A_d = (2.0 - Ts) / (2.0 + Ts);       // 0.904762...
        double expected_B_d = Ts / (1.0 + Ts / 2.0);         // 0.095238...
        double expected_C_d = 1.0 / (1.0 + Ts / 2.0);        // 0.952381...
        double expected_D_d = (Ts / 2.0) / (1.0 + Ts / 2.0); // 0.047619...

        CHECK(sys_d.A(0, 0) == doctest::Approx(expected_A_d).epsilon(1e-12));
        CHECK(sys_d.B(0, 0) == doctest::Approx(expected_B_d).epsilon(1e-12));
        CHECK(sys_d.C(0, 0) == doctest::Approx(expected_C_d).epsilon(1e-12));
        CHECK(sys_d.D(0, 0) == doctest::Approx(expected_D_d).epsilon(1e-12));

        // Discrete pole should be inside unit circle (stable system stays stable)
        CHECK(std::abs(sys_d.A(0, 0)) < 1.0);
    }

    // Test 3: ZOH vs Tustin comparison using unified API
    TEST_CASE("ZOH and Tustin produce different discretizations") {
        double     Ts = 0.5;
        StateSpace sys{
            Matrix<1, 1>{{-2.0}},
            Matrix<1, 1>{{1.0}},
            Matrix<1, 1>{{1.0}}
        };

        StateSpace sys_d_zoh = *discretize(sys, Ts, DiscretizationMethod::ZOH);
        StateSpace sys_d_tustin = *discretize(sys, Ts, DiscretizationMethod::Tustin);

        // They should produce different results
        CHECK_NE(sys_d_zoh.A(0, 0), doctest::Approx(sys_d_tustin.A(0, 0)));
    }

    // Test 4: Default method is ZOH
    TEST_CASE("eval_frf returns nullopt at a continuous pole") {
        StateSpace<1, 1, 1> sys{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}};
        // s = -1 is a pole of 1/(s+1)
        CHECK_FALSE(eval_frf(sys, damp::complex<double>{-1.0, 0.0}).has_value());
        // s = j is fine
        CHECK(eval_frf(sys, damp::complex<double>{0.0, 1.0}).has_value());
    }

    TEST_CASE("Default discretization method is ZOH") {
        double     Ts = 0.1;
        StateSpace sys{
            Matrix<1, 1>{{-1.0}},
            Matrix<1, 1>{{1.0}},
            Matrix<1, 1>{{1.0}}
        };

        StateSpace sys_d_default = *discretize(sys, Ts);
        StateSpace sys_d_zoh = *discretize(sys, Ts, DiscretizationMethod::ZOH);

        CHECK(doctest::Approx(sys_d_default.A(0, 0)) == sys_d_zoh.A(0, 0));
        CHECK(doctest::Approx(sys_d_default.B(0, 0)) == sys_d_zoh.B(0, 0));
    }

    // ZOH maps process-noise input G with the same ∫e^{Aτ}G dτ as B
    TEST_CASE("ZOH G_d matches B_d integral for identical columns") {
        double                            Ts = 0.1;
        StateSpace<1, 1, 1, double, 1, 0> sys{
            .A = {{-1.0}},
            .B = {{1.0}},
            .C = {{1.0}},
            .D = {{0.0}},
            .G = {{1.0}}, // same column as B
        };

        auto sys_d = *discretize(sys, Ts, DiscretizationMethod::ZOH);

        // G_d must equal B_d, not the Forward-Euler G*Ts
        CHECK(sys_d.G(0, 0) == doctest::Approx(sys_d.B(0, 0)).epsilon(1e-12));
        CHECK(sys_d.G(0, 0) == doctest::Approx(1.0 - std::exp(-Ts)).epsilon(1e-10));
        CHECK(std::abs(sys_d.G(0, 0) - Ts) > 1e-6); // not FE
    }

    TEST_CASE("ForwardEuler: Ad = I + A Ts, Bd = B Ts") {
        const double Ts = 0.1;
        StateSpace   sys{
              .A = Matrix<1, 1>{{-1.0}},
              .B = Matrix<1, 1>{{1.0}},
              .C = Matrix<1, 1>{{1.0}},
        };
        const auto sys_d = *discretize(sys, Ts, DiscretizationMethod::ForwardEuler);
        CHECK(sys_d.A(0, 0) == doctest::Approx(1.0 - Ts)); // 0.9
        CHECK(sys_d.B(0, 0) == doctest::Approx(Ts));
        CHECK(sys_d.Ts == doctest::Approx(Ts));
    }

    TEST_CASE("discretize Ts<=0 or already-discrete returns input unchanged") {
        StateSpace cont{
            .A = Matrix<1, 1>{{-1.0}},
            .B = Matrix<1, 1>{{1.0}},
            .C = Matrix<1, 1>{{1.0}},
            .Ts = 0.0,
        };
        const auto same = *discretize(cont, 0.0, DiscretizationMethod::ZOH);
        CHECK(same.Ts == doctest::Approx(0.0));
        CHECK(same.A(0, 0) == doctest::Approx(-1.0));

        StateSpace disc{
            .A = Matrix<1, 1>{{-0.5}},
            .B = Matrix<1, 1>{{0.1}},
            .C = Matrix<1, 1>{{1.0}},
            .Ts = 0.01,
        };
        const auto passthrough = *discretize(disc, 0.1, DiscretizationMethod::ZOH);
        CHECK(passthrough.Ts == doctest::Approx(0.01));
        CHECK(passthrough.A(0, 0) == doctest::Approx(-0.5));
        CHECK(passthrough.B(0, 0) == doctest::Approx(0.1));
    }
}
