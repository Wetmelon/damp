// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <numbers>

#include "damp/backend.hpp"
#include "damp/controllers/lead_lag.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

/**
 * @brief Tests for Lead-Lag compensator (lead_lag.hpp)
 *
 * Validates lead/lag design, discretization, and step response.
 * Reference values generated with scipy.signal.
 */

TEST_SUITE("Lead-Lag Compensator") {

    TEST_CASE("Lead design - 45 degree phase boost at wc=1000") {
        constexpr double phi = std::numbers::pi / 4.0; // 45 degrees
        constexpr double wc = 1000.0;
        constexpr auto   result = design::lead(phi, wc);

        // Check zero, pole, gain against analytical values
        // alpha = (1-sin(45))/(1+sin(45)) = 0.17157...
        // sqrt(alpha) = 0.41421...
        static_assert(result.z < result.p); // Lead: zero < pole
        CHECK(result.z == doctest::Approx(414.2135623730950).epsilon(1e-12));
        CHECK(result.p == doctest::Approx(2414.213562373095).epsilon(1e-12));
        CHECK(result.K == doctest::Approx(2.414213562373095).epsilon(1e-12));
    }

    TEST_CASE("Lead design - phase at wc is exactly phi_max") {
        // C(jw) = K*(jw + z)/(jw + p)
        // Phase at wc should be phi_max, magnitude should be 1.0
        constexpr double phi = std::numbers::pi / 4.0;
        constexpr double wc = 1000.0;
        constexpr auto   r = design::lead(phi, wc);

        // Phase = atan(wc/z) - atan(wc/p)
        double phase_z = std::atan(wc / r.z);
        double phase_p = std::atan(wc / r.p);
        double phase_deg = (phase_z - phase_p) * 180.0 / std::numbers::pi;
        CHECK(phase_deg == doctest::Approx(45.0).epsilon(1e-10));

        // Magnitude = K * sqrt(wc²+z²) / sqrt(wc²+p²)
        double mag = r.K * std::sqrt((wc * wc) + (r.z * r.z)) / std::sqrt((wc * wc) + (r.p * r.p));
        CHECK(mag == doctest::Approx(1.0).epsilon(1e-12));
    }

    TEST_CASE("Lead design - 30 degree boost at wc=500") {
        constexpr double phi = 30.0 * std::numbers::pi / 180.0;
        constexpr double wc = 500.0;
        constexpr auto   r = design::lead(phi, wc);

        CHECK(r.z == doctest::Approx(288.6751345948128).epsilon(1e-12));
        CHECK(r.p == doctest::Approx(866.0254037844387).epsilon(1e-12));
        CHECK(r.K == doctest::Approx(1.732050807568877).epsilon(1e-12));

        // Verify phase and magnitude at wc
        double phase_z = std::atan(wc / r.z);
        double phase_p = std::atan(wc / r.p);
        double phase_deg = (phase_z - phase_p) * 180.0 / std::numbers::pi;
        CHECK(phase_deg == doctest::Approx(30.0).epsilon(1e-10));

        double mag = r.K * std::sqrt((wc * wc) + (r.z * r.z)) / std::sqrt((wc * wc) + (r.p * r.p));
        CHECK(mag == doctest::Approx(1.0).epsilon(1e-12));
    }

    TEST_CASE("Lag design - 10x DC gain boost") {
        constexpr double dc_boost = 10.0;
        constexpr double wc = 1000.0;
        constexpr auto   r = design::lag(dc_boost, wc);

        // z = wc/margin_factor = 1000/10 = 100
        // p = z/dc_gain_boost = 100/10 = 10
        static_assert(r.z > r.p); // Lag: zero > pole
        CHECK(r.z == doctest::Approx(100.0).epsilon(1e-12));
        CHECK(r.p == doctest::Approx(10.0).epsilon(1e-12));
        CHECK(r.K == doctest::Approx(1.0).epsilon(1e-12));

        // DC gain = K * z / p = 10
        double dc_gain = r.K * r.z / r.p;
        CHECK(dc_gain == doctest::Approx(10.0).epsilon(1e-12));

        // Phase at wc should be small (< 6 degrees lag)
        double phase_z = std::atan(wc / r.z);
        double phase_p = std::atan(wc / r.p);
        double phase_deg = (phase_z - phase_p) * 180.0 / std::numbers::pi;
        CHECK(phase_deg == doctest::Approx(-5.137654439816157).epsilon(1e-6));
    }

    TEST_CASE("Lead-lag direct specification") {
        constexpr auto r = design::lead_lag_direct(2.0, 100.0, 500.0);
        static_assert(r.K == 2.0);
        static_assert(r.z == 100.0);
        static_assert(r.p == 500.0);
        static_assert(r.Ts == 0.0);
    }

    TEST_CASE("LeadLagResult to_tf") {
        constexpr auto r = design::lead_lag_direct(2.0, 100.0, 500.0);
        constexpr auto tf = r.to_tf();

        // C(s) = 2*(s+100)/(s+500) = (200 + 2s)/(500 + s)
        CHECK(tf.num[0] == doctest::Approx(200.0).epsilon(1e-12)); // K*z
        CHECK(tf.num[1] == doctest::Approx(2.0).epsilon(1e-12));   // K
        CHECK(tf.den[0] == doctest::Approx(500.0).epsilon(1e-12)); // p
        CHECK(tf.den[1] == doctest::Approx(1.0).epsilon(1e-12));   // 1
    }

    TEST_CASE("LeadLagResult to_ss") {
        constexpr auto r = design::lead_lag_direct(2.0, 100.0, 500.0);
        constexpr auto ss = r.to_ss();

        // A = [-p] = [-500], B = [1], C = [K*(z-p)] = [2*(100-500)] = [-800], D = [K] = [2]
        CHECK(ss.A(0, 0) == doctest::Approx(-500.0).epsilon(1e-12));
        CHECK(ss.B(0, 0) == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(ss.C(0, 0) == doctest::Approx(-800.0).epsilon(1e-12));
        CHECK(ss.D(0, 0) == doctest::Approx(2.0).epsilon(1e-12));
    }

    TEST_CASE("Tustin discretization matches scipy") {
        constexpr double phi = std::numbers::pi / 4.0;
        constexpr double wc = 1000.0;
        constexpr double Ts = 0.001;

        constexpr auto            r = design::lead(phi, wc, Ts);
        LeadLagController<double> ctrl(r);

        // Coefficients verified against scipy.signal.cont2discrete bilinear
        CHECK(ctrl.b0 == doctest::Approx(1.320377241017041).epsilon(1e-12));
        CHECK(ctrl.b1 == doctest::Approx(-0.8672954016950677).epsilon(1e-12));
        CHECK(ctrl.a1 == doctest::Approx(0.09383632135605431).epsilon(1e-12));
    }

    TEST_CASE("Discrete step response matches scipy") {
        constexpr double phi = std::numbers::pi / 4.0;
        constexpr double wc = 1000.0;
        constexpr double Ts = 0.001;

        constexpr auto            r = design::lead(phi, wc, Ts);
        LeadLagController<double> ctrl(r);

        // Step response reference from scipy simulation
        constexpr double expected[] = {
            1.320377241017041e+00,
            3.291824962226776e-01,
            4.221925648216337e-01,
            4.134648421352334e-01,
            4.142838195259410e-01,
        };

        for (int i = 0; i < 5; ++i) {
            double y = ctrl.control(1.0);
            CHECK(y == doctest::Approx(expected[i]).epsilon(1e-12));
        }
    }

    TEST_CASE("DC gain consistency") {
        // Lead: DC gain = K*z/p < 1 (attenuates DC)
        constexpr double phi = std::numbers::pi / 4.0;
        constexpr double wc = 1000.0;
        constexpr auto   lead_r = design::lead(phi, wc);
        double           dc_lead = lead_r.K * lead_r.z / lead_r.p;
        CHECK(dc_lead == doctest::Approx(0.4142135623730950).epsilon(1e-12));
        CHECK(dc_lead < 1.0); // Lead attenuates DC

        // Lag: DC gain = K*z/p > 1 (boosts DC)
        constexpr auto lag_r = design::lag(10.0, wc);
        double         dc_lag = lag_r.K * lag_r.z / lag_r.p;
        CHECK(dc_lag == doctest::Approx(10.0).epsilon(1e-12));
        CHECK(dc_lag > 1.0); // Lag boosts DC
    }

    TEST_CASE("Controller reset clears state") {
        constexpr double          phi = std::numbers::pi / 4.0;
        constexpr auto            r = design::lead(phi, 1000.0, 0.001);
        LeadLagController<double> ctrl(r);

        // Run some steps
        (void)ctrl.control(1.0);
        (void)ctrl.control(1.0);
        (void)ctrl.control(1.0);

        ctrl.reset();

        // After reset, first step should be same as fresh controller
        LeadLagController<double> fresh(r);
        CHECK(ctrl.control(1.0) == doctest::Approx(fresh.control(1.0)).epsilon(1e-15));
    }

    TEST_CASE("Online lead design matches design") {
        constexpr double phi = std::numbers::pi / 4.0;
        constexpr double wc = 1000.0;
        constexpr auto   design_r = design::lead(phi, wc);
        auto             online_r = design::lead(phi, wc);

        CHECK(online_r.K == doctest::Approx(design_r.K).epsilon(1e-12));
        CHECK(online_r.z == doctest::Approx(design_r.z).epsilon(1e-12));
        CHECK(online_r.p == doctest::Approx(design_r.p).epsilon(1e-12));
    }

    TEST_CASE("Online lag design matches design") {
        constexpr double wc = 1000.0;
        constexpr auto   design_r = design::lag(10.0, wc);
        auto             online_r = design::lag(10.0, wc);

        CHECK(online_r.K == doctest::Approx(design_r.K).epsilon(1e-12));
        CHECK(online_r.z == doctest::Approx(design_r.z).epsilon(1e-12));
        CHECK(online_r.p == doctest::Approx(design_r.p).epsilon(1e-12));
    }

    TEST_CASE("Lead-lag combined (2nd order)") {
        // 45-degree lead + 10x lag at wc=1000, continuous
        constexpr double phi = std::numbers::pi / 4.0;
        constexpr double wc = 1000.0;
        constexpr auto   result = design::lead_lag(phi, wc, 10.0);
        static_assert(result.success);
        constexpr auto ss = result.sys;

        // Should be a 2x2 state-space system (series of two 1st-order)
        // Check dimensions compile
        static_assert(ss.A.rows() == 2);
        static_assert(ss.A.cols() == 2);
        static_assert(ss.B.rows() == 2);
        static_assert(ss.B.cols() == 1);
        static_assert(ss.C.rows() == 1);
        static_assert(ss.C.cols() == 2);

        // DC gain of combined = lead_dc * lag_dc = 0.41421... * 10 = 4.1421...
        // Compute via state-space: -C * A^-1 * B + D
        auto A_inv = ss.A.inverse();
        REQUIRE(A_inv.has_value());
        auto dc = (ss.C * (-A_inv.value()) * ss.B + ss.D);
        CHECK(dc(0, 0) == doctest::Approx(0.4142135623730950 * 10.0).epsilon(1e-10));
    }

    TEST_CASE("Lead-lag series propagates section failure") {
        // Invalid lead phi -> success false, zero system
        static_assert(!design::lead_lag(0.0, 1000.0, 10.0).success);
        static_assert(!design::lead_lag(std::numbers::pi / 4.0, 1000.0, 1.0).success); // lag needs boost > 1
        constexpr auto bad = design::lead_lag(-1.0, 1000.0, 10.0);
        static_assert(!bad.success);
        CHECK(bad.sys.A(0, 0) == doctest::Approx(0.0));
        CHECK(bad.sys.D(0, 0) == doctest::Approx(0.0));
    }

    TEST_CASE("Type conversion via .as<float>()") {
        constexpr auto r = design::lead(std::numbers::pi / 4.0, 1000.0, 0.001);
        constexpr auto rf = r.as<float>();

        LeadLagController<float> ctrl(rf);
        float                    y = ctrl.control(1.0f);
        CHECK(y == doctest::Approx(1.320377f).epsilon(1e-4));
    }

    TEST_CASE("control(r, y) equals control(r - y)") {
        constexpr double          phi = std::numbers::pi / 4.0;
        constexpr auto            r = design::lead(phi, 1000.0, 0.001);
        LeadLagController<double> a(r);
        LeadLagController<double> b(r);

        // Feed the same error through both forms; outputs must match step for step.
        const double refs[] = {1.0, 1.0, 0.5, -0.3, 0.0};
        const double meas[] = {0.0, 0.4, 0.6, 0.1, -0.2};
        for (int i = 0; i < 5; ++i) {
            double ya = a.control(refs[i], meas[i]);
            double yb = b.control(refs[i] - meas[i]);
            CHECK(ya == doctest::Approx(yb).epsilon(1e-15));
        }
    }

    TEST_CASE("Invalid specs report success = false") {
        // Lead requires 0 < phi < pi/2 and wc > 0.
        static_assert(!design::lead(0.0, 1000.0).success);
        static_assert(!design::lead(std::numbers::pi / 2.0, 1000.0).success);
        static_assert(!design::lead(std::numbers::pi / 4.0, -1.0).success);
        static_assert(design::lead(std::numbers::pi / 4.0, 1000.0).success);

        // Lag requires dc_gain_boost > 1, wc > 0, margin > 0.
        static_assert(!design::lag(1.0, 1000.0).success);
        static_assert(!design::lag(10.0, 0.0).success);
        static_assert(design::lag(10.0, 1000.0).success);

        static_assert(!design::lead_lag_direct(1.0, -1.0, 5.0).success);
        static_assert(design::lead_lag_direct(1.0, 100.0, 500.0).success);
    }

    TEST_CASE("Failed design constructs fail-closed controller (emits zero)") {
        constexpr auto bad = design::lead(0.0, 1000.0, 0.001);
        static_assert(!bad.success);
        LeadLagController<double> ctrl(bad);
        CHECK_FALSE(ctrl.valid());
        CHECK(ctrl.control(1.0) == doctest::Approx(0.0));
        CHECK(ctrl.control(2.0, 1.0) == doctest::Approx(0.0));
        ctrl.reset();
        CHECK(ctrl.control(-3.0) == doctest::Approx(0.0));
    }

    TEST_CASE("to_discrete_ss roundtrip (dimensions and sample time)") {
        constexpr double phi = std::numbers::pi / 4.0;
        constexpr double wc = 1000.0;
        constexpr double Ts = 0.001;
        constexpr auto   r = design::lead(phi, wc, Ts);

        // Get discrete SS
        constexpr auto dss = r.to_discrete_ss();

        // Verify dimensions
        static_assert(dss.A.rows() == 1);
        static_assert(dss.A.cols() == 1);
        static_assert(dss.B.rows() == 1);
        static_assert(dss.B.cols() == 1);
        static_assert(dss.C.rows() == 1);
        static_assert(dss.C.cols() == 1);

        // Verify sampling time set
        CHECK(dss.Ts == doctest::Approx(Ts).epsilon(1e-15));
    }

    TEST_CASE("to_discrete_ss matches controller step response") {
        constexpr double phi = std::numbers::pi / 4.0;
        constexpr double wc = 1000.0;
        constexpr double Ts = 0.001;
        constexpr auto   r = design::lead(phi, wc, Ts);

        // Discrete SS via Tustin
        constexpr auto dss = r.to_discrete_ss();

        // Verify SS matrices against scipy (bilinear)
        CHECK(dss.A(0, 0) == doctest::Approx(-9.383632135605437e-02).epsilon(1e-12));
        CHECK(dss.B(0, 0) == doctest::Approx(4.530818393219728e-04).epsilon(1e-12));
        CHECK(dss.C(0, 0) == doctest::Approx(-2.187672642712108e+03).epsilon(1e-6));
        CHECK(dss.D(0, 0) == doctest::Approx(1.320377241017041e+00).epsilon(1e-12));

        // Step response from state-space should match IIR controller
        LeadLagController<double> ctrl(r);

        Matrix<1, 1, double> x = Matrix<1, 1, double>::zeros();
        Matrix<1, 1, double> u{{1.0}};

        for (int i = 0; i < 5; ++i) {
            auto   y_ss = dss.C * x + dss.D * u;
            double y_ctrl = ctrl.control(1.0);
            CHECK(y_ss(0, 0) == doctest::Approx(y_ctrl).epsilon(1e-10));
            x = dss.A * x + dss.B * u;
        }
    }
}
