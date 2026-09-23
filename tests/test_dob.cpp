// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <array>
#include <cmath>

#include "damp/estimation/dob.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using estimation::DOB;
using estimation::DOBConfig;

// The disturbance observer is a lightweight SISO estimator:
//   innovation = y_measured − y_predicted
//   d_hat[k+1] = (1 − leak)·d_hat[k] + gain·innovation
// with optional innovation deadband and output-magnitude clamp. It is the
// estimator core for the planned DOB control law. These tests cover config
// validation, the recursion's steady state, the deadband, the clamp, and the
// compensate()/reset() runtime surface.
//
// @see Li et al., "Disturbance Observer-Based Control" (CRC Press, 2016).

TEST_SUITE("Disturbance Observer") {
    TEST_CASE("config validation rejects out-of-range parameters") {
        // Valid baseline.
        CHECK(DOBConfig<double>{0.1, 0.0, 0.0, 0.0, false}.valid());

        // gain must be in [0, 1].
        CHECK_FALSE(DOBConfig<double>{-0.1}.valid());
        CHECK_FALSE(DOBConfig<double>{1.5}.valid());

        // leak must be in [0, 1).
        CHECK_FALSE(DOBConfig<double>{0.1, -0.1}.valid());
        CHECK_FALSE(DOBConfig<double>{0.1, 1.0}.valid());

        // negative deadband / magnitude are invalid.
        CHECK_FALSE(DOBConfig<double>{0.1, 0.0, -1.0}.valid());
        CHECK_FALSE(DOBConfig<double>{0.1, 0.0, 0.0, -1.0}.valid());
    }

    TEST_CASE("synthesize reports steady-state gain and rejects bad config") {
        constexpr DOBConfig<double> good{0.2, 0.1, 0.0, 0.0, false};
        constexpr auto              res = estimation::dob(good);
        static_assert(res.success);
        // Constant-innovation steady state: d_hat_ss = gain/leak · innovation.
        CHECK(res.steady_state_gain == doctest::Approx(0.2 / 0.1));

        constexpr DOBConfig<double> bad{2.0}; // gain > 1
        constexpr auto              res_bad = estimation::dob(bad);
        static_assert(!res_bad.success);
    }

    TEST_CASE("leaky integrator settles at (gain/leak)·innovation") {
        // The recursion d_hat[k+1] = (1−leak)·d_hat[k] + gain·innovation is a
        // leaky integrator: with a constant innovation it settles at
        // (gain/leak)·innovation. Choosing gain = leak gives unity DC gain, so
        // the estimate converges to the innovation itself. (With leak = 0 it is a
        // pure integrator and would diverge on a constant innovation — that is
        // the intended physics, not a settling estimator.)
        DOB<double> dob{DOBConfig<double>{0.25, 0.25, 0.0, 0.0, false}};

        const double d_true = 3.0;
        for (int k = 0; k < 500; ++k) {
            // y_predicted = 0, y_measured = d_true → innovation = d_true each tick.
            REQUIRE(dob.update(0.0, d_true));
        }
        CHECK(dob.state().disturbance_hat == doctest::Approx(d_true).epsilon(1e-3));
    }

    TEST_CASE("compensate subtracts the disturbance estimate from the command") {
        // gain = leak ⇒ unity-DC-gain leaky integrator settles d_hat at the
        // constant innovation (here 2.0).
        DOB<double> dob{DOBConfig<double>{0.5, 0.5, 0.0, 0.0, false}};
        for (int k = 0; k < 200; ++k) {
            REQUIRE(dob.update(0.0, 2.0)); // drive estimate toward +2
        }
        // u_compensated = u_nominal − d_hat ≈ 10 − 2 = 8.
        CHECK(dob.compensate(10.0) == doctest::Approx(8.0).epsilon(1e-2));
    }

    TEST_CASE("innovation deadband freezes the estimate for small errors") {
        DOB<double> dob{DOBConfig<double>{0.5, 0.0, 0.1, 0.0, false}};
        // Innovation of 0.05 is below the 0.1 deadband → estimate stays put.
        for (int k = 0; k < 50; ++k) {
            REQUIRE(dob.update(0.0, 0.05));
        }
        CHECK(dob.state().disturbance_hat == doctest::Approx(0.0));

        // Innovation above the deadband moves the estimate.
        REQUIRE(dob.update(0.0, 1.0));
        CHECK(dob.state().disturbance_hat > 0.0);
    }

    TEST_CASE("magnitude clamp bounds the estimate") {
        DOB<double> dob{DOBConfig<double>{0.5, 0.0, 0.0, 1.0, true}};
        for (int k = 0; k < 200; ++k) {
            REQUIRE(dob.update(0.0, 100.0)); // would blow past the clamp
        }
        CHECK(dob.state().disturbance_hat <= doctest::Approx(1.0));
        CHECK(dob.state().disturbance_hat > 0.0);
    }

    TEST_CASE("reset clears state") {
        DOB<double> dob{DOBConfig<double>{0.5, 0.0, 0.0, 0.0, false}};
        REQUIRE(dob.update(0.0, 5.0));
        REQUIRE(dob.state().disturbance_hat != 0.0);

        dob.reset();
        CHECK(dob.state().disturbance_hat == doctest::Approx(0.0));
        CHECK_FALSE(dob.state().initialized);
    }

    TEST_CASE("constructing from a design result carries the config") {
        constexpr DOBConfig<double> cfg{0.3, 0.05, 0.0, 0.0, false};
        const auto                  design = estimation::dob(cfg);
        REQUIRE(design.success);

        DOB<double> dob{design};
        CHECK(dob.valid());
        CHECK(dob.config().gain == doctest::Approx(0.3));
        CHECK(dob.config().leak == doctest::Approx(0.05));
    }

    TEST_CASE("as<float>() converts the design result") {
        constexpr DOBConfig<double> cfg{0.2, 0.1, 0.01, 5.0, true};
        const auto                  res = estimation::dob(cfg);
        const auto                  rf = res.as<float>();
        CHECK(rf.success == res.success);
        CHECK(rf.config.gain == doctest::Approx(0.2));
        CHECK(rf.config.clamp_enabled);
    }
}

// ---------------------------------------------------------------------------
// Classical Pn^-1 * Q disturbance observer
// ---------------------------------------------------------------------------

using estimation::classical_dob;
using estimation::ClassicalDOB;

namespace {
// First-order discrete plant y[k] = alpha*y[k-1] + beta*(u[k] + d[k]).
// As a z^-1 transfer function from (u+d) to y:  Pn = beta / (1 - alpha z^-1).
constexpr double alpha = 0.9;
constexpr double beta = 0.1;
// 1st-order low-pass Q = (1 - rho) / (1 - rho z^-1)  (unity DC gain).
constexpr double rho = 0.85;

constexpr auto make_dob() {
    const damp::array<double, 1> Bn{beta};
    const damp::array<double, 2> An{1.0, -alpha};
    const damp::array<double, 1> Qn{1.0 - rho};
    const damp::array<double, 2> Qd{1.0, -rho};
    return classical_dob(Bn, An, Qn, Qd);
}
} // namespace

TEST_SUITE("Classical DOB") {
    TEST_CASE("synthesize forms Q*Pn^-1 and validates realizability") {
        const auto d = make_dob();
        REQUIRE(d.success);
        // Fy = Q*Pn^-1 = (Qn * An)/(Qd * Bn). Numerator = (1-rho)*[1,-alpha].
        CHECK(d.fy_num[0] == doctest::Approx((1.0 - rho)));
        CHECK(d.fy_num[1] == doctest::Approx(-(1.0 - rho) * alpha));
        CHECK(d.fy_den[0] == doctest::Approx(beta));        // Qd[0]*Bn[0] = 1*beta
        CHECK(d.fy_den[1] == doctest::Approx(-rho * beta)); // Qd[1]*Bn[0]
        CHECK(d.fu_num[0] == doctest::Approx(1.0 - rho));   // Fu = Q
        CHECK(d.q_tf().num[0] == doctest::Approx(d.fu_num[0]));
        CHECK(d.q_tf().den[1] == doctest::Approx(d.fu_den[1]));
        CHECK(d.fy_tf().num[0] == doctest::Approx(d.fy_num[0]));
        CHECK(d.fy_tf().den[0] == doctest::Approx(d.fy_den[0]));

        // Non-invertible nominal plant (Bn[0] == 0) is rejected.
        const damp::array<double, 2> Bn_delay{0.0, beta};
        const damp::array<double, 2> An{1.0, -alpha};
        const damp::array<double, 1> Qn{1.0 - rho};
        const damp::array<double, 2> Qd{1.0, -rho};
        CHECK_FALSE(classical_dob(Bn_delay, An, Qn, Qd).success);
    }

    TEST_CASE("rejects a step load disturbance (nominal model matches plant)") {
        ClassicalDOB<1, 2, 1, 2, double> dob(make_dob());
        REQUIRE(dob.valid());

        // Closed loop: base command 0, constant input-referred disturbance d = 1.
        // The DOB should drive y back to 0 and estimate d_hat -> 1.
        double       y = 0.0;
        const double d = 1.0;
        double       u = 0.0;
        for (int k = 0; k < 400; ++k) {
            u = dob.compensate(0.0, y);     // u = -d_hat
            y = alpha * y + beta * (u + d); // plant
        }
        CHECK(std::abs(y) < 1e-3);                                      // disturbance rejected
        CHECK(dob.disturbance() == doctest::Approx(1.0).epsilon(1e-2)); // d_hat -> d
    }

    TEST_CASE("still rejects a DC disturbance under model mismatch") {
        // Actual plant gain 20% higher than the nominal model used to design.
        ClassicalDOB<1, 2, 1, 2, double> dob(make_dob());
        const double                     alpha_act = 0.88;
        const double                     beta_act = 0.12;
        double                           y = 0.0;
        const double                     d = 1.0;
        for (int k = 0; k < 1500; ++k) {
            const double u = dob.compensate(0.0, y);
            y = alpha_act * y + beta_act * (u + d);
        }
        CHECK(std::abs(y) < 5e-2); // DC rejection survives the mismatch
    }

    TEST_CASE("estimate() recovers a known disturbance open-loop") {
        ClassicalDOB<1, 2, 1, 2, double> dob(make_dob());
        const double                     d = 0.5;
        double                           y = 0.0;
        double                           last = 0.0;
        for (int k = 0; k < 400; ++k) {
            const double u = 0.2;           // fixed input
            y = alpha * y + beta * (u + d); // plant with disturbance
            last = dob.estimate(y, u);
        }
        CHECK(last == doctest::Approx(0.5).epsilon(1e-2));
    }

    TEST_CASE("invalid design is inert; reset clears filter state") {
        const damp::array<double, 1>     Bad{0.0}; // Bn[0] == 0 -> invalid
        const damp::array<double, 2>     An{1.0, -alpha};
        const damp::array<double, 1>     Qn{1.0 - rho};
        const damp::array<double, 2>     Qd{1.0, -rho};
        ClassicalDOB<1, 2, 1, 2, double> bad(classical_dob(Bad, An, Qn, Qd));
        CHECK_FALSE(bad.valid());
        CHECK(bad.compensate(3.0, 1.0) == doctest::Approx(3.0)); // pass-through

        ClassicalDOB<1, 2, 1, 2, double> dob(make_dob());
        (void)dob.compensate(0.0, 1.0);
        dob.reset();
        CHECK(dob.disturbance() == doctest::Approx(0.0));
    }

    TEST_CASE("float deployment via as<float>()") {
        ClassicalDOB<1, 2, 1, 2, float> dob(make_dob().as<float>());
        REQUIRE(dob.valid());
        float y = 0.0f;
        for (int k = 0; k < 400; ++k) {
            const float u = dob.compensate(0.0f, y);
            y = 0.9f * y + 0.1f * (u + 1.0f);
        }
        CHECK(std::abs(y) < 1e-2f);
    }

    TEST_CASE("classical DOB is constexpr-evaluable") {
        constexpr double y_final = []() consteval {
            ClassicalDOB<1, 2, 1, 2, double> dob(make_dob());
            double                           y = 0.0;
            for (int k = 0; k < 400; ++k) {
                const double u = dob.compensate(0.0, y);
                y = 0.9 * y + 0.1 * (u + 1.0);
            }
            return y;
        }();
        static_assert(y_final < 1e-3 && y_final > -1e-3, "DOB must reject the DC load at compile time");
        CHECK(y_final == doctest::Approx(0.0).epsilon(1e-3));
    }
}
