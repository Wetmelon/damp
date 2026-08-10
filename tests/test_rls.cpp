// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>

#include "damp/estimation/rls.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using namespace damp::estimation;

TEST_SUITE("Recursive Least Squares") {
    TEST_CASE("Scalar RLS converges to the true parameter") {
        RlsConfig<double> cfg{};
        cfg.lambda = 1.0;
        cfg.p0 = 1000.0; // large prior covariance: trust the data

        Rls<double> rls(cfg);

        // y = theta_true * phi, theta_true = 3
        constexpr double theta_true = 3.0;
        for (int k = 1; k <= 50; ++k) {
            const double phi = 0.1 * k;
            REQUIRE(rls.update(phi, theta_true * phi));
        }

        CHECK(rls.state().theta == doctest::Approx(theta_true).epsilon(1e-4));
        CHECK(rls.predict(2.0) == doctest::Approx(2.0 * theta_true).epsilon(1e-4));
        CHECK(rls.valid());
    }

    TEST_CASE("Scalar RLS rejects an invalid configuration") {
        RlsConfig<double> cfg{};
        cfg.lambda = 1.5; // > 1 is invalid
        CHECK_FALSE(cfg.valid());

        Rls<double> rls(cfg);
        CHECK_FALSE(rls.update(1.0, 1.0));
        CHECK_FALSE(rls.valid());
    }

    TEST_CASE("Scalar RLS projection clamps the estimate") {
        RlsConfig<double> cfg{};
        cfg.lambda = 1.0;
        cfg.p0 = 1000.0;
        cfg.projection_enabled = true;
        cfg.theta_min = 0.0;
        cfg.theta_max = 1.0;

        Rls<double> rls(cfg);
        // True parameter (5) is well above theta_max → estimate must be clamped.
        for (int k = 1; k <= 50; ++k) {
            const double phi = 0.1 * k;
            REQUIRE(rls.update(phi, 5.0 * phi));
        }
        CHECK(rls.state().theta <= 1.0);
        CHECK(rls.state().theta >= 0.0);
    }

    TEST_CASE("Scalar RLS reset restores the prior") {
        RlsConfig<double> cfg{};
        cfg.p0 = 10.0;
        Rls<double> rls(cfg);

        REQUIRE(rls.update(1.0, 2.0));
        rls.reset(0.5);
        CHECK(rls.state().theta == doctest::Approx(0.5));
        CHECK(rls.state().covariance == doctest::Approx(10.0));
    }

    TEST_CASE("Vector RLS converges to the true parameter vector") {
        RlsConfig<double> cfg{};
        cfg.lambda = 1.0;
        cfg.p0 = 1000.0;

        RlsVector<2, double> rls(cfg);

        const ColVec<2, double> theta_true{2.0, -1.0};
        // Persistently exciting regressors.
        for (int k = 0; k < 100; ++k) {
            const double      a = damp::sin(0.3 * k);
            const double      b = damp::cos(0.17 * k);
            ColVec<2, double> phi{a, b};
            const double      y = dot(phi, theta_true);
            REQUIRE(rls.update(phi, y));
        }

        CHECK(rls.state().theta[0] == doctest::Approx(2.0).epsilon(1e-4));
        CHECK(rls.state().theta[1] == doctest::Approx(-1.0).epsilon(1e-4));

        // Covariance stays symmetric after the update.
        const auto& P = rls.state().covariance;
        CHECK(P(0, 1) == doctest::Approx(P(1, 0)).epsilon(1e-9));
    }

    TEST_CASE("Vector RLS with forgetting tracks a parameter jump") {
        RlsConfig<double> cfg{};
        cfg.lambda = 0.9; // forgetting factor: discount old data
        cfg.p0 = 100.0;

        RlsVector<2, double> rls(cfg);

        ColVec<2, double> theta{1.0, 1.0};
        for (int k = 0; k < 400; ++k) {
            if (k == 200) {
                theta = ColVec<2, double>{-2.0, 3.0}; // parameters jump
            }
            ColVec<2, double> phi{damp::sin(0.3 * k), damp::cos(0.11 * k)};
            REQUIRE(rls.update(phi, dot(phi, theta)));
        }

        // After the jump, forgetting lets the estimate catch up.
        CHECK(rls.state().theta[0] == doctest::Approx(-2.0).epsilon(1e-3));
        CHECK(rls.state().theta[1] == doctest::Approx(3.0).epsilon(1e-3));
    }

    TEST_CASE("RLS design payload round-trips through as<float>()") {
        RlsConfig<double> cfg{};
        cfg.lambda = 0.98;
        cfg.p0 = 50.0;
        const auto design = rls(cfg, 1.25);
        REQUIRE(design.success);

        const auto design_f = design.as<float>();
        Rls<float> rls(design_f);
        CHECK(rls.state().theta == doctest::Approx(1.25f));
        CHECK(rls.valid());
    }

    TEST_CASE("Scalar RLS zero regressor does not move theta") {
        RlsConfig<double> cfg{};
        cfg.lambda = 1.0;
        cfg.p0 = 1.0;
        Rls<double> rls(cfg);
        rls.reset(3.0);
        // phi = 0 ⇒ gain = 0; denom = λ > 0 so the sample is accepted but is not PE.
        REQUIRE(rls.update(0.0, 99.0));
        CHECK(rls.state().theta == doctest::Approx(3.0));
        CHECK(rls.state().gain == doctest::Approx(0.0));
    }

    TEST_CASE("Vector RLS zero regressor leaves theta unchanged") {
        RlsConfig<double> cfg{};
        cfg.lambda = 1.0;
        cfg.p0 = 10.0;
        RlsVector<2, double> rls(cfg);
        REQUIRE(rls.update(ColVec<2, double>{1.0, 0.0}, 1.0));
        const auto theta0 = rls.state().theta;
        REQUIRE(rls.update(ColVec<2, double>{0.0, 0.0}, 99.0)); // no PE this sample
        CHECK(rls.state().theta[0] == doctest::Approx(theta0[0]));
        CHECK(rls.state().theta[1] == doctest::Approx(theta0[1]));
    }
}
