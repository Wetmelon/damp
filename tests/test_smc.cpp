// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)


#include <cmath>
#include <limits>

#include "damp/controllers/smc.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

TEST_SUITE("Sliding Mode Control (SMC)") {
    TEST_CASE("SMC gain computation") {
        float lambda = 10.0f; // Sliding surface parameter
        float k = 5.0f;       // Switching gain
        float b0 = 1.0f;      // Plant gain

        auto smc_result = design::smc(lambda, k, b0);

        CHECK(smc_result.lambda == doctest::Approx(lambda).epsilon(1e-6));
        CHECK(smc_result.k == doctest::Approx(k).epsilon(1e-6));
        CHECK(smc_result.b0 == doctest::Approx(b0).epsilon(1e-6));
    }

    TEST_CASE("SMC controller basic functionality") {
        float lambda = 10.0f;
        float k = 5.0f;
        float b0 = 1.0f;
        float Ts = 0.01f;

        auto          smc_result = design::smc(lambda, k, b0, Ts);
        SMCController controller(smc_result);

        // Test control with zero error
        float r = 0.0f;
        float y = 0.0f;
        float u = controller.control(r, y);
        CHECK(u == 0.0f); // s = 0, sign(0) = 0, u = 0

        // Test control with positive error
        r = 1.0f;
        y = 0.0f;
        u = controller.control(r, y);
        CHECK(u == doctest::Approx(k / b0).epsilon(1e-6)); // s > 0, u = +k/b0

        // Test control with negative error
        r = 0.0f;
        y = 1.0f;
        u = controller.control(r, y);
        CHECK(u == doctest::Approx(-k / b0).epsilon(1e-6)); // s < 0, u = -k/b0
    }

    TEST_CASE("SMC boundary layer: continuous, saturates to the relay magnitude") {
        const float          Ts = 0.01f;
        const float          phi = 0.5f;
        const auto           g = design::smc(10.0f, 5.0f, 1.0f, Ts, phi);
        SMCController<float> ctrl(g);

        (void)ctrl.control(0.0f, 0.0f); // first tick seeds (error 0)

        // Large surface -> sat() clamps to the full relay magnitude k/b0.
        const float u_big = ctrl.control(100.0f, 0.0f);
        CHECK(u_big == doctest::Approx(5.0f).epsilon(1e-4));

        // Small surface -> proportional region, |u| strictly below the relay, but nonzero.
        SMCController<float> c2(g);
        (void)c2.control(0.0f, 0.0f);
        const float u_small = c2.control(1e-4f, 0.0f);
        CHECK(damp::abs(u_small) < 5.0f);
        CHECK(u_small != 0.0f);
    }

    TEST_CASE("SMC reset clears the rate history and re-seeds the first tick") {
        SMCController<double> ctrl(design::smc(10.0, 5.0, 1.0, 0.01));
        (void)ctrl.control(1.0, 0.0);
        (void)ctrl.control(0.5, 0.0);
        ctrl.reset();
        // First post-reset tick re-seeds: e_dot = 0, s = lambda*error = 10 > 0, so
        // the relay gives exactly +k/b0 with no derivative kick carried from history.
        CHECK(ctrl.control(1.0, 0.0) == doctest::Approx(5.0));
    }

    TEST_CASE("SMC: a continuous or invalid design is inert (no divide-by-zero)") {
        SMCController<double> continuous(design::smc(10.0, 5.0, 1.0)); // Ts = 0
        CHECK_FALSE(continuous.valid());                               // used to report valid() and emit nothing
        CHECK(continuous.control(1.0, 0.0) == doctest::Approx(0.0));   // no sampled rate: hold
        CHECK_FALSE(design::smc(std::numeric_limits<double>::quiet_NaN(), 5.0, 1.0, 0.01).success);

        SMCController<double> bad(design::smc(10.0, 5.0, 0.0, 0.01)); // b0 = 0 -> invalid
        CHECK_FALSE(bad.valid());
        CHECK(bad.control(1.0, 0.0) == doctest::Approx(0.0)); // no k/b0 divide-by-zero

        CHECK(design::smc(10.0, 5.0, 1.0, -0.01).status == design::DesignStatus::invalid_sample_time);
        CHECK(design::smc(10.0, 5.0, 1.0, 0.01, -0.1).status == design::DesignStatus::invalid_parameter);
    }

    TEST_CASE("SMC converting ctor preserves gains and rate state") {
        SMCController<double> d(design::smc(10.0, 5.0, 1.0, 0.01, 100.0));
        (void)d.control(1.0, 0.0); // error_prev now 1, first_ cleared
        SMCController<float> f(d); // converting ctor (exercises friend access)
        CHECK(f.valid());

        // In the boundary layer the output depends on s (hence error_prev); the
        // copied state reproduces the double's next command within float precision.
        const double u_d = d.control(1.2, 0.0);
        const float  u_f = f.control(1.2f, 0.0f);
        CHECK(static_cast<double>(u_f) == doctest::Approx(u_d).epsilon(1e-3));
        CHECK(u_d != doctest::Approx(0.0));
    }

    TEST_CASE("SMC as<float>() down-casts the design via static_cast") {
        constexpr auto g = design::smc(10.0, 5.0, 1.0, 1e-3, 0.05);
        constexpr auto gf = g.as<float>();
        static_assert(gf.success);
        CHECK(gf.lambda == doctest::Approx(10.0f));
        CHECK(gf.k == doctest::Approx(5.0f));
        CHECK(gf.b0 == doctest::Approx(1.0f));
        CHECK(gf.phi == doctest::Approx(0.05f));
        CHECK(gf.Ts == doctest::Approx(1e-3f));
    }

    TEST_CASE("boundary-layer to_tf is (k/(b0 φ))(s + λ)") {
        const double   phi = 0.05;
        constexpr auto g = design::smc(20.0, 5.0, 2.0, 0.0, 0.05);
        const auto     tf = g.to_tf();
        const double   gn = 5.0 / (2.0 * phi);
        CHECK(tf.num[0] == doctest::Approx(gn));        // s
        CHECK(tf.num[1] == doctest::Approx(gn * 20.0)); // λ
        CHECK(tf.den[0] == doctest::Approx(1.0));
        CHECK(design::smc(20.0, 5.0, 2.0).to_tf().num[1] == doctest::Approx(0.0));
    }

    TEST_CASE("SMC (r, y) regulates a positive-gain double integrator") {
        // ÿ = b·u with b = 2 > 0: SMC must push y toward r, not away.
        const double          b = 2.0;
        const double          Ts = 1e-3;
        SMCController<double> ctrl(design::smc(5.0, 5.0, b, Ts, 0.05));
        double                y = 0.0;
        double                v = 0.0;
        for (int k = 0; k < 6000; ++k) {
            const double u = ctrl.control(1.0, y);
            v += Ts * b * u;
            y += Ts * v;
        }
        CHECK(std::abs(y - 1.0) < 0.05);
    }
}
