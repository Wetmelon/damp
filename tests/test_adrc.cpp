// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/controllers/adrc.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

TEST_SUITE("Active Disturbance Rejection Control (ADRC)") {
    TEST_CASE("1st-order ADRC gain computation") {
        constexpr size_t NX = 1;

        double wc = 5.0;  // Controller bandwidth
        double wo = 20.0; // Luenberger bandwidth
        double b0 = 1.0;  // Plant gain

        auto adrc_result = design::adrc<NX>(wc, wo, b0);
        CHECK(adrc_result.success);

        // Expected ESO gains for 1st-order system with poles at -wo
        damp::array<double, NX + 1> expected_beta = {
            2 * wo,
            wo * wo,
        };

        CHECK(adrc_result.beta.size() == expected_beta.size());
        for (size_t i = 0; i < expected_beta.size(); ++i) {
            CHECK(adrc_result.beta[i] == doctest::Approx(expected_beta[i]).epsilon(1e-6));
        }

        // Check Kp and Kd gains for 1st order
        double expected_Kp = wc / b0;
        double expected_Kd = 0.0;

        CHECK(adrc_result.Kp == doctest::Approx(expected_Kp).epsilon(1e-6));
        CHECK(adrc_result.Kd == doctest::Approx(expected_Kd).epsilon(1e-6));
    }

    TEST_CASE("2nd-order ADRC gain computation") {
        constexpr size_t NX = 2;

        double wc = 5.0;  // Controller bandwidth
        double wo = 20.0; // Luenberger bandwidth
        double b0 = 1.0;  // Plant gain

        auto adrc_result = design::adrc<NX>(wc, wo, b0);
        CHECK(adrc_result.success);

        // Expected ESO gains for 2nd-order system with poles at -wo
        damp::array<double, NX + 1> expected_beta = {
            3 * wo,
            3 * wo * wo,
            wo * wo * wo,
        };

        CHECK(adrc_result.beta.size() == expected_beta.size());
        for (size_t i = 0; i < expected_beta.size(); ++i) {
            CHECK(adrc_result.beta[i] == doctest::Approx(expected_beta[i]).epsilon(1e-6));
        }

        // Check Kp and Kd gains
        double expected_Kp = (wc * wc) / b0;
        double expected_Kd = (2 * wc) / b0;

        CHECK(adrc_result.Kp == doctest::Approx(expected_Kp).epsilon(1e-6));
        CHECK(adrc_result.Kd == doctest::Approx(expected_Kd).epsilon(1e-6));
    }

    TEST_CASE("invalid b0 / bandwidths report success=false and zero gains") {
        auto bad_b0 = design::adrc<1>(5.0, 20.0, 0.0);
        CHECK_FALSE(bad_b0.success);
        CHECK(bad_b0.Kp == 0.0);
        CHECK(bad_b0.b0 == 0.0);

        auto bad_wc = design::adrc<1>(-1.0, 20.0, 1.0);
        CHECK_FALSE(bad_wc.success);

        auto bad_wo = design::adrc<2>(5.0, 0.0, 1.0);
        CHECK_FALSE(bad_wo.success);

        // Runtime from a failed design must not produce Inf/NaN.
        ADRCController<1, double> ctrl(bad_b0);
        CHECK_FALSE(ctrl.valid());
        CHECK(ctrl.control(1.0, 0.0, 1e-3) == 0.0);
    }

    TEST_CASE("converting ctor preserves gains and ESO state") {
        ADRCController<1, double> d(design::adrc<1>(20.0, 100.0, 1.0), 1e-3);
        (void)d.control(1.0, 0.2, 1e-3);
        (void)d.control(1.0, 0.4, 1e-3);

        ADRCController<1, float> f(d);
        CHECK(f.valid());
        // Next step matches within float precision (state + gains carried over).
        CHECK(static_cast<double>(f.control(1.0f, 0.5f, 1e-3f)) == doctest::Approx(d.control(1.0, 0.5, 1e-3)).epsilon(1e-5));
    }

    TEST_CASE("ADRCResult::as propagates success") {
        constexpr auto ok = design::adrc<1>(5.0, 20.0, 1.0);
        static_assert(ok.success);
        constexpr auto ok_f = ok.as<float>();
        static_assert(ok_f.success);
        CHECK(ok_f.Kp == doctest::Approx(static_cast<float>(ok.Kp)));

        constexpr auto bad = design::adrc<1>(5.0, 20.0, -1.0);
        static_assert(!bad.success);
        static_assert(!bad.as<float>().success);
    }

    TEST_CASE("1st-order runtime: tracks setpoint and rejects a constant disturbance") {
        const double Ts = 1e-3;
        auto         ctrl = ADRCController<1, double>(design::adrc<1>(20.0, 100.0, 1.0));

        // Plant: ẏ = −a·y + b·u + d. The −a·y term and the constant d are both
        // unmodeled, lumped into the ESO's total-disturbance estimate.
        const double a = 2.0, b = 1.0, d = 0.5, r = 1.0;
        double       y = 0.0;
        for (int k = 0; k < 20000; ++k) {
            const double u = ctrl.control(r, y, Ts);
            y += (-a * y + b * u + d) * Ts;
        }
        CHECK(y == doctest::Approx(r).epsilon(0.01)); // zero steady-state error
    }

    TEST_CASE("2nd-order runtime: double integrator tracks with disturbance rejection") {
        const double Ts = 1e-3;
        auto         ctrl = ADRCController<2, double>(design::adrc<2>(10.0, 50.0, 1.0));

        // Plant: ÿ = b·u + d (double integrator + constant load disturbance).
        const double b = 1.0, d = 0.3, r = 1.0;
        double       y = 0.0, v = 0.0;
        for (int k = 0; k < 30000; ++k) {
            const double u = ctrl.control(r, y, Ts);
            v += (b * u + d) * Ts;
            y += v * Ts;
        }
        CHECK(y == doctest::Approx(r).epsilon(0.01));   // converged to setpoint
        CHECK(v == doctest::Approx(0.0).epsilon(1e-3)); // and settled (no residual velocity)
    }

    TEST_CASE("reset clears the observer state") {
        const double Ts = 1e-3;
        auto         ctrl = ADRCController<1, double>(design::adrc<1>(20.0, 100.0, 1.0));
        for (int k = 0; k < 100; ++k) {
            (void)ctrl.control(1.0, 0.0, Ts);
        }
        ctrl.reset();
        // After reset the ESO state is zero, so the first command equals the
        // fresh-controller command for the same inputs.
        ADRCController<1, double> fresh(design::adrc<1>(20.0, 100.0, 1.0));
        CHECK(ctrl.control(1.0, 0.0, Ts) == doctest::Approx(fresh.control(1.0, 0.0, Ts)));
    }

    TEST_CASE("2-arg control uses stored Ts; back_calculate overrides u_prev") {
        const double Ts = 1e-3;
        auto         res = design::adrc<1>(20.0, 100.0, 1.0);
        REQUIRE(res.success);
        ADRCController<1, double> ctrl{res, Ts}; // Siso form stores Ts

        const double u3 = ctrl.control(1.0, 0.0, Ts);
        ctrl.reset();
        const double u2 = ctrl.control(1.0, 0.0);
        CHECK(u2 == doctest::Approx(u3));

        // Warm the ESO with a nonzero measurement so u_prev enters the update.
        ctrl.reset();
        for (int k = 0; k < 20; ++k) {
            (void)ctrl.control(1.0, 0.5);
        }
        const double u_free = ctrl.control(1.0, 0.5);
        ctrl.back_calculate(u_free, 0.0); // plant saturated to 0
        // Same inputs; ESO next tick uses u_prev=0 instead of u_free
        const double u_after = ctrl.control(1.0, 0.5);
        CHECK(u_after != doctest::Approx(u_free));
    }
}
