// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>

#include "damp/controllers/smc.hpp"
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

        auto          smc_result = design::smc(lambda, k, b0);
        SMCController controller(smc_result);

        // Test control with zero error
        float r = 0.0f;
        float y = 0.0f;
        float u = controller.control(r, y, Ts);
        CHECK(u == 0.0f); // s = 0, sign(0) = 0, u = 0

        // Test control with positive error
        r = 1.0f;
        y = 0.0f;
        u = controller.control(r, y, Ts);
        CHECK(u == doctest::Approx(-k / b0).epsilon(1e-6)); // s > 0, u = -k/b0

        // Test control with negative error
        r = 0.0f;
        y = 1.0f;
        u = controller.control(r, y, Ts);
        CHECK(u == doctest::Approx(k / b0).epsilon(1e-6)); // s < 0, u = k/b0
    }

    TEST_CASE("SMC boundary layer: continuous, saturates to the relay magnitude") {
        const auto           g = design::smc(10.0f, 5.0f, 1.0f);
        SMCController<float> ctrl(g);
        const float          Ts = 0.01f;
        const float          phi = 0.5f;

        (void)ctrl.control(0.0f, 0.0f, Ts, phi); // first tick seeds (error 0)

        // Large surface -> sat() clamps to the full relay magnitude k/b0.
        const float u_big = ctrl.control(100.0f, 0.0f, Ts, phi);
        CHECK(u_big == doctest::Approx(-5.0f).epsilon(1e-4));

        // Small surface -> proportional region, |u| strictly below the relay, but nonzero.
        SMCController<float> c2(g);
        (void)c2.control(0.0f, 0.0f, Ts, phi);
        const float u_small = c2.control(1e-4f, 0.0f, Ts, phi);
        CHECK(damp::abs(u_small) < 5.0f);
        CHECK(u_small != 0.0f);
    }

    TEST_CASE("SMC reset clears the rate history and re-seeds the first tick") {
        SMCController<double> ctrl(design::smc(10.0, 5.0, 1.0));
        (void)ctrl.control(1.0, 0.0, 0.01);
        (void)ctrl.control(0.5, 0.0, 0.01);
        ctrl.reset();
        // First post-reset tick re-seeds: e_dot = 0, s = lambda*error = 10 > 0, so
        // the relay gives exactly -k/b0 with no derivative kick carried from history.
        CHECK(ctrl.control(1.0, 0.0, 0.01) == doctest::Approx(-5.0));
    }

    TEST_CASE("SMC: non-positive Ts and invalid design are inert (no divide-by-zero)") {
        SMCController<double> good(design::smc(10.0, 5.0, 1.0));
        CHECK(good.valid());
        CHECK(good.control(1.0, 0.0, 0.0) == doctest::Approx(0.0));   // Ts = 0 -> hold
        CHECK(good.control(1.0, 0.0, -0.01) == doctest::Approx(0.0)); // Ts < 0 -> hold

        SMCController<double> bad(design::smc(10.0, 5.0, 0.0)); // b0 = 0 -> invalid
        CHECK_FALSE(bad.valid());
        CHECK(bad.control(1.0, 0.0, 0.01) == doctest::Approx(0.0)); // no k/b0 divide-by-zero
    }

    TEST_CASE("SMC converting ctor preserves gains and rate state") {
        SMCController<double> d(design::smc(10.0, 5.0, 1.0));
        (void)d.control(1.0, 0.0, 0.01); // error_prev now 1, first_ cleared
        SMCController<float> f(d);       // converting ctor (exercises friend access)
        CHECK(f.valid());

        // In the boundary layer the output depends on s (hence error_prev); the
        // copied state reproduces the double's next command within float precision.
        const double u_d = d.control(1.2, 0.0, 0.01, 100.0);
        const float  u_f = f.control(1.2f, 0.0f, 0.01f, 100.0f);
        CHECK(static_cast<double>(u_f) == doctest::Approx(u_d).epsilon(1e-3));
        CHECK(u_d != doctest::Approx(0.0));
    }

    TEST_CASE("SMC as<float>() down-casts the design via static_cast") {
        constexpr auto g = design::smc(10.0, 5.0, 1.0);
        constexpr auto gf = g.as<float>();
        static_assert(gf.success);
        CHECK(gf.lambda == doctest::Approx(10.0f));
        CHECK(gf.k == doctest::Approx(5.0f));
        CHECK(gf.b0 == doctest::Approx(1.0f));
    }

    TEST_CASE("boundary-layer to_tf is (k/(b0 φ))(s + λ)") {
        constexpr auto g = design::smc(20.0, 5.0, 2.0);
        const double   phi = 0.05;
        const auto     tf = g.to_tf(phi);
        const double   gn = 5.0 / (2.0 * phi);
        CHECK(tf.num[0] == doctest::Approx(gn * 20.0));
        CHECK(tf.num[1] == doctest::Approx(gn));
        CHECK(tf.den[0] == doctest::Approx(1.0));
        CHECK(g.to_tf(0.0).num[1] == doctest::Approx(0.0));
    }
}
