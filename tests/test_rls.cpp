// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)


#include "damp/estimation/rls.hpp"
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

        RlsVector<1, double> rls(rls_vector(cfg));

        // y = theta_true * phi, theta_true = 3
        constexpr double theta_true = 3.0;
        for (int k = 1; k <= 50; ++k) {
            const double phi = 0.1 * k;
            REQUIRE(rls.update(phi, theta_true * phi));
        }

        CHECK(rls.state().theta[0] == doctest::Approx(theta_true).epsilon(1e-4));
        CHECK(rls.predict(2.0) == doctest::Approx(2.0 * theta_true).epsilon(1e-4));
        CHECK(rls.valid());
    }

    TEST_CASE("Scalar RLS rejects an invalid configuration") {
        RlsConfig<double> cfg{};
        cfg.lambda = 1.5; // > 1 is invalid
        CHECK_FALSE(cfg.valid());

        RlsVector<1, double> rls(rls_vector(cfg));
        CHECK_FALSE(rls.valid());
        CHECK_FALSE(rls.update(1.0, 1.0));
        rls.reset();
        CHECK_FALSE(rls.valid()); // reset keeps the failed design inert
        CHECK_FALSE(rls.update(1.0, 1.0));
    }

    TEST_CASE("Scalar RLS projection clamps the estimate") {
        RlsConfig<double> cfg{};
        cfg.lambda = 1.0;
        cfg.p0 = 1000.0;
        cfg.projection_enabled = true;
        cfg.theta_lim = {
            .lower = 0.0,
            .upper = 1.0,
        };

        RlsVector<1, double> rls(rls_vector(cfg));
        // True parameter (5) is well above theta.upper → estimate must be clamped.
        for (int k = 1; k <= 50; ++k) {
            const double phi = 0.1 * k;
            REQUIRE(rls.update(phi, 5.0 * phi));
        }
        CHECK(rls.state().theta[0] <= 1.0);
        CHECK(rls.state().theta[0] >= 0.0);
    }

    TEST_CASE("Scalar RLS reset restores the prior") {
        RlsConfig<double> cfg{};
        cfg.p0 = 10.0;
        RlsVector<1, double> rls(rls_vector(cfg));

        REQUIRE(rls.update(1.0, 2.0));
        rls.reset(0.5);
        CHECK(rls.state().theta[0] == doctest::Approx(0.5));
        CHECK(rls.state().covariance(0, 0) == doctest::Approx(10.0));
    }

    TEST_CASE("Vector RLS converges to the true parameter vector") {
        RlsConfig<double, 2> cfg{};
        cfg.lambda = 1.0;
        cfg.p0 = 1000.0;

        RlsVector<2, double> rls(rls_vector(cfg));

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
        RlsConfig<double, 2> cfg{};
        cfg.lambda = 0.9; // forgetting factor: discount old data
        cfg.p0 = 100.0;

        RlsVector<2, double> rls(rls_vector(cfg));

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
        const auto design = rls_vector(cfg, ColVec<1, double>{1.25});
        REQUIRE(design.success);

        const auto          design_f = design.as<float>();
        RlsVector<1, float> rls(design_f);
        CHECK(rls.state().theta[0] == doctest::Approx(1.25f));
        CHECK(rls.state().covariance(0, 0) == doctest::Approx(50.0f));
        CHECK(rls.valid());
    }

    TEST_CASE("Scalar RLS zero regressor does not move theta") {
        RlsConfig<double> cfg{};
        cfg.lambda = 1.0;
        cfg.p0 = 1.0;
        RlsVector<1, double> rls(rls_vector(cfg));
        rls.reset(3.0);
        // phi = 0 ⇒ gain = 0; denom = λ > 0 so the sample is accepted but is not PE.
        REQUIRE(rls.update(0.0, 99.0));
        CHECK(rls.state().theta[0] == doctest::Approx(3.0));
        CHECK(rls.state().gain[0] == doctest::Approx(0.0));
    }

    TEST_CASE("Vector RLS zero regressor leaves theta unchanged") {
        RlsConfig<double, 2> cfg{};
        cfg.lambda = 1.0;
        cfg.p0 = 10.0;
        RlsVector<2, double> rls(rls_vector(cfg));
        REQUIRE(rls.update(ColVec<2, double>{1.0, 0.0}, 1.0));
        const auto theta0 = rls.state().theta;
        REQUIRE(rls.update(ColVec<2, double>{0.0, 0.0}, 99.0)); // no PE this sample
        CHECK(rls.state().theta[0] == doctest::Approx(theta0[0]));
        CHECK(rls.state().theta[1] == doctest::Approx(theta0[1]));
    }

    TEST_CASE("RLS forgetting does not wind the covariance up without excitation") {
        // λ < 1 with no information divides P by λ every sample: unbounded growth
        // (inf in float after a few hundred samples, then NaN gains).
        estimation::RlsConfig<float> cfg{};
        cfg.lambda = 0.9f;
        cfg.p0 = 10.0f;
        estimation::RlsVector<1, float> rls(rls_vector(cfg));
        estimation::RlsVector<2, float> rlsv(rls_vector(RlsConfig<float, 2>{.lambda = cfg.lambda, .p0 = cfg.p0}));
        for (int k = 0; k < 1000; ++k) {
            REQUIRE(rls.update(0.0f, 1.0f));
            REQUIRE(rlsv.update(ColVec<2, float>{0.0f, 0.0f}, 1.0f));
        }
        CHECK(rls.state().covariance(0, 0) <= cfg.p0 * (1.0f + 1e-6f));
        CHECK(rlsv.state().covariance.trace() <= 2.0f * cfg.p0 * (1.0f + 1e-6f));
        // Still able to learn afterwards.
        for (int k = 0; k < 50; ++k) {
            REQUIRE(rls.update(1.0f, 2.0f));
        }
        CHECK(rls.state().theta[0] == doctest::Approx(2.0f).epsilon(1e-3));
    }

    TEST_CASE("RLS accepts a valid tiny forgetting factor") {
        // λ + φᵀPφ ≥ λ > 0 for any valid config; an absolute threshold rejected it.
        estimation::RlsConfig<double> cfg{};
        cfg.lambda = 1e-13;
        cfg.p0 = 1.0;
        estimation::RlsVector<1, double> rls(rls_vector(cfg));
        CHECK(rls.update(0.0, 1.0));
    }

    TEST_CASE("RlsVector deploys with static constinit") {
        // The state is built in the member-initializer list: an assignment there compared the
        // constinit object's address with the design's, which g++ rejects as non-constant.
        static constexpr auto                             design = rls_vector(RlsConfig<double, 2>{.lambda = 0.99, .p0 = 100.0});
        static constinit estimation::RlsVector<2, double> r{design};
        static constinit estimation::RlsVector<2, float>  rf{design};
        CHECK(r.valid());
        CHECK(r.state().covariance(1, 1) == doctest::Approx(100.0));
        CHECK(rf.state().covariance(0, 0) == doctest::Approx(100.0f));
    }
}
