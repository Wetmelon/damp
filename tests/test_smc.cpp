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
    TEST_CASE("SMC InputLimits clamps the relay; limits that are not a box are refused") {
        const auto res = design::smc({.lambda = 20.0, .k = 5.0, .b0 = 1.0}, 1e-3, design::InputLimits{.lower = {-1.0}, .upper = {1.0}});
        REQUIRE(res.success);
        SMCController<double> law{res};
        CHECK(law.step(1.0, 0.0) == 1.0); // k/b₀ = 5 past the limit
        CHECK(law.step(-1.0, 0.0) == -1.0);
        CHECK(design::smc({.lambda = 20.0, .k = 5.0, .b0 = 1.0}, 1e-3, design::InputLimits{.lower = {1.0}, .upper = {-1.0}}).status == design::DesignStatus::invalid_parameter);
        auto by_hand = res;
        by_hand.limits.input = {.lower = 1.0, .upper = -1.0};
        SMCController<double> bad{by_hand};
        CHECK_FALSE(bad.valid());
        CHECK(bad.step(1.0, 0.0) == 0.0);
    }

    TEST_CASE("SMC gain computation") {
        float lambda = 10.0f; // Sliding surface parameter
        float k = 5.0f;       // Switching gain
        float b0 = 1.0f;      // Plant gain

        auto smc_result = design::smc({.lambda = lambda, .k = k, .b0 = b0}, 0.01f);

        CHECK(smc_result.gains.lambda == doctest::Approx(lambda).epsilon(1e-6));
        CHECK(smc_result.gains.k == doctest::Approx(k).epsilon(1e-6));
        CHECK(smc_result.gains.b0 == doctest::Approx(b0).epsilon(1e-6));
    }

    TEST_CASE("SMC controller basic functionality") {
        float lambda = 10.0f;
        float k = 5.0f;
        float b0 = 1.0f;
        float Ts = 0.01f;

        auto          smc_result = design::smc({.lambda = lambda, .k = k, .b0 = b0}, Ts);
        SMCController controller(smc_result);

        // Test control with zero error
        float r = 0.0f;
        float y = 0.0f;
        float u = controller.step(r, y);
        CHECK(u == 0.0f); // s = 0, sign(0) = 0, u = 0

        // Test control with positive error
        r = 1.0f;
        y = 0.0f;
        u = controller.step(r, y);
        CHECK(u == doctest::Approx(k / b0).epsilon(1e-6)); // s > 0, u = +k/b0

        // Test control with negative error
        r = 0.0f;
        y = 1.0f;
        u = controller.step(r, y);
        CHECK(u == doctest::Approx(-k / b0).epsilon(1e-6)); // s < 0, u = -k/b0
        CHECK(controller.last_command() == u);
    }

    TEST_CASE("SMC boundary layer: continuous, saturates to the relay magnitude") {
        const float          Ts = 0.01f;
        const float          phi = 0.5f;
        const auto           g = design::smc({.lambda = 10.0f, .k = 5.0f, .b0 = 1.0f, .phi = phi}, Ts);
        SMCController<float> ctrl(g);

        (void)ctrl.step(0.0f, 0.0f); // first tick seeds (error 0)

        // Large surface -> sat() clamps to the full relay magnitude k/b0.
        const float u_big = ctrl.step(100.0f, 0.0f);
        CHECK(u_big == doctest::Approx(5.0f).epsilon(1e-4));

        // Small surface -> proportional region, |u| strictly below the relay, but nonzero.
        SMCController<float> c2(g);
        (void)c2.step(0.0f, 0.0f);
        const float u_small = c2.step(1e-4f, 0.0f);
        CHECK(damp::abs(u_small) < 5.0f);
        CHECK(u_small != 0.0f);
    }

    TEST_CASE("SMC reset clears the rate history and re-seeds the first tick") {
        SMCController<double> ctrl(design::smc({.lambda = 10.0, .k = 5.0, .b0 = 1.0}, 0.01));
        (void)ctrl.step(1.0, 0.0);
        (void)ctrl.step(0.5, 0.0);
        ctrl.reset();
        // First post-reset tick re-seeds: e_dot = 0, s = lambda*error = 10 > 0, so
        // the relay gives exactly +k/b0 with no derivative kick carried from history.
        CHECK(ctrl.step(1.0, 0.0) == doctest::Approx(5.0));
        ctrl.reset();
        CHECK(ctrl.last_command() == 0.0);
    }

    TEST_CASE("SMC: an invalid design is inert (no divide-by-zero); Ts must be finite and positive") {
        //! Ts = 0 used to be the continuous design, for analysis; the runtime refused it. Every Ts
        //! that is not finite and positive is refused now, and to_tf needs none.
        CHECK(design::smc({.lambda = 10.0, .k = 5.0, .b0 = 1.0}, 0.0).status == design::DesignStatus::invalid_sample_time);
        CHECK(design::smc({.lambda = 10.0, .k = 5.0, .b0 = 1.0}, -0.01).status == design::DesignStatus::invalid_sample_time);
        CHECK(design::smc({.lambda = 10.0, .k = 5.0, .b0 = 1.0}, std::numeric_limits<double>::quiet_NaN()).status == design::DesignStatus::invalid_sample_time);
        SMCController<double> continuous(design::smc({.lambda = 10.0, .k = 5.0, .b0 = 1.0}, 0.0));
        CHECK_FALSE(continuous.valid());
        CHECK(continuous.step(1.0, 0.0) == doctest::Approx(0.0)); // no sampled rate: hold
        CHECK_FALSE(design::smc({.lambda = std::numeric_limits<double>::quiet_NaN(), .k = 5.0, .b0 = 1.0}, 0.01).success);

        SMCController<double> bad(design::smc({.lambda = 10.0, .k = 5.0, .b0 = 0.0}, 0.01)); // b0 = 0 -> invalid
        CHECK_FALSE(bad.valid());
        CHECK(bad.step(1.0, 0.0) == doctest::Approx(0.0)); // no k/b0 divide-by-zero

        CHECK(design::smc({.lambda = 10.0, .k = 5.0, .b0 = 1.0, .phi = -0.1}, 0.01).status == design::DesignStatus::invalid_parameter);
    }

    TEST_CASE("SMC converting ctor preserves gains and rate state") {
        SMCController<double> d(design::smc({.lambda = 10.0, .k = 5.0, .b0 = 1.0, .phi = 100.0}, 0.01));
        (void)d.step(1.0, 0.0);    // error_prev now 1, first_ cleared
        SMCController<float> f(d); // converting ctor (exercises friend access)
        CHECK(f.valid());

        // In the boundary layer the output depends on s (hence error_prev); the
        // copied state reproduces the double's next command within float precision.
        const double u_d = d.step(1.2, 0.0);
        const float  u_f = f.step(1.2f, 0.0f);
        CHECK(static_cast<double>(u_f) == doctest::Approx(u_d).epsilon(1e-3));
        CHECK(u_d != doctest::Approx(0.0));
    }

    TEST_CASE("SMC as<float>() down-casts the design via static_cast") {
        constexpr auto g = design::smc({.lambda = 10.0, .k = 5.0, .b0 = 1.0, .phi = 0.05}, 1e-3);
        constexpr auto gf = g.as<float>();
        static_assert(gf.success);
        CHECK(gf.gains.lambda == doctest::Approx(10.0f));
        CHECK(gf.gains.k == doctest::Approx(5.0f));
        CHECK(gf.gains.b0 == doctest::Approx(1.0f));
        CHECK(gf.gains.phi == doctest::Approx(0.05f));
        CHECK(gf.Ts == doctest::Approx(1e-3f));
    }

    TEST_CASE("boundary-layer to_tf is (k/(b0 φ))(s + λ)") {
        const double   phi = 0.05;
        constexpr auto g = design::smc({.lambda = 20.0, .k = 5.0, .b0 = 2.0, .phi = 0.05}, 1e-3);
        const auto     tf = g.to_tf();
        const double   gn = 5.0 / (2.0 * phi);
        CHECK(tf.num[0] == doctest::Approx(gn));        // s
        CHECK(tf.num[1] == doctest::Approx(gn * 20.0)); // λ
        CHECK(tf.den[0] == doctest::Approx(1.0));
        CHECK(design::smc({.lambda = 20.0, .k = 5.0, .b0 = 2.0}, 1e-3).to_tf().num[1] == doctest::Approx(0.0));
    }

    TEST_CASE("SMC (r, y) regulates a positive-gain double integrator") {
        // ÿ = b·u with b = 2 > 0: SMC must push y toward r, not away.
        const double          b = 2.0;
        const double          Ts = 1e-3;
        SMCController<double> ctrl(design::smc({.lambda = 5.0, .k = 5.0, .b0 = b, .phi = 0.05}, Ts));
        double                y = 0.0;
        double                v = 0.0;
        for (int k = 0; k < 6000; ++k) {
            const double u = ctrl.step(1.0, y);
            v += Ts * b * u;
            y += Ts * v;
        }
        CHECK(std::abs(y - 1.0) < 0.05);
    }

    TEST_CASE("SMC track keeps the rate current: the step after it has no derivative kick") {
        //! The relay has no memory to preload; what track must not leave behind is a stale last
        //! error, which the next step would difference into a kick of (Δe)/Ts.
        const auto            res = design::smc({.lambda = 10.0, .k = 5.0, .b0 = 1.0, .phi = 100.0}, 0.01);
        SMCController<double> tracked{res};
        SMCController<double> fresh{res};
        (void)tracked.step(0.0, 0.0);
        for (int k = 0; k < 20; ++k) {
            tracked.track(2.5, 1.0, 0.05 * static_cast<double>(k)); // an operator drives the plant
        }
        CHECK(tracked.last_command() == 2.5);
        const double y = 0.05 * 19.0;
        CHECK(tracked.step(1.0, y) == doctest::Approx(fresh.step(1.0, y)).epsilon(1e-12)); // ė = 0, as on a first tick
    }

    TEST_CASE("SMC set_gains retunes the relay and set_limits moves the clamp; a bad value is refused") {
        SMCController<double> law{design::smc({.lambda = 10.0, .k = 5.0, .b0 = 1.0}, 0.01)};
        CHECK(law.step(1.0, 0.0) == 5.0);
        auto g = law.gains();
        g.k = 2.0;
        g.b0 = 4.0;
        REQUIRE(law.set_gains(g));
        CHECK(law.step(1.0, 0.0) == doctest::Approx(0.5)); // k/b₀
        auto lim = law.limits();
        lim.input = {.lower = -0.25, .upper = 0.25};
        REQUIRE(law.set_limits(lim));
        CHECK(law.step(1.0, 0.0) == 0.25);
        auto neg = g;
        neg.lambda = -1.0;
        CHECK_FALSE(law.set_gains(neg));
        auto thin = g;
        thin.phi = -0.1;
        CHECK_FALSE(law.set_gains(thin));
        CHECK(law.gains().k == 2.0);
        CHECK(law.sample_time() == 0.01);
    }

    TEST_CASE("a constexpr SMC design's Law: its gains constants of the type, its limits settable") {
        static constexpr auto kSmc = design::smc({.lambda = 10.0, .k = 5.0, .b0 = 1.0}, 0.01, design::InputLimits{.lower = {-2.0}, .upper = {2.0}});
        Law<kSmc, float>      law;
        SMCController<float>  held{kSmc};
        static_assert(sizeof(Law<kSmc, float>) < sizeof(SMCController<float>));
        CHECK(law.step(1.0f, 0.0f) == held.step(1.0f, 0.0f));
        CHECK(law.step(1.0f, 0.0f) == 2.0f);
        CHECK(law.gains().k == 5.0f);
        auto lim = law.limits();
        lim.input.upper(0) = 1.0f;
        REQUIRE(law.set_limits(lim));
        CHECK(law.step(1.0f, 0.0f) == 1.0f);
        Law<kSmc, float>::Runtime retunable{kSmc.as<float>()};
        CHECK(retunable.set_gains(retunable.gains()));
    }
}
