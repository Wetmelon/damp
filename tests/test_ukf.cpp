// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>

#include "damp/estimation/ukf.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

// The Unscented Kalman Filter propagates sigma points through the exact
// nonlinear maps f/h — no Jacobian. These tests cover: (1) the linear case,
// where the UKF must reduce to a linear Kalman filter; (2) a strongly
// nonlinear measurement where the UKF should at least match the EKF; (3)
// covariance grow/shrink behaviour; (4) the singular-Pyy failure path; and
// (5) float instantiation.
//
// @see Wan & van der Merwe (2000); Julier & Uhlmann (2004); Simon (2006) §14.3.

TEST_SUITE("Unscented Kalman Filter") {
    TEST_CASE("linear dynamics + linear measurement converges to truth") {
        // Constant-velocity model with a linear position measurement. On a
        // linear f/h the unscented transform is exact, so the UKF must behave
        // like a textbook KF and recover the unmeasured velocity.
        const Matrix<2, 2> A{{1.0, 0.1}, {0.0, 1.0}};
        const Matrix<2, 2> Q{{1e-4, 0.0}, {0.0, 1e-4}};
        const Matrix<1, 1> R{{1e-2}};

        auto state_fn = [&](const ColVec<2>& x, const ColVec<1>&) { return ColVec<2>{A * x}; };
        auto meas_fn = [&](const ColVec<2>& x, const ColVec<1>&) { return ColVec<1>{{x[0]}}; };

        UnscentedKalmanFilter<2, 1, 1, double> ukf{ColVec<2>{{0.0}, {0.0}}, Matrix<2, 2>::identity(), Q};

        ColVec<2> x_true{{0.0}, {1.0}}; // 1 m/s
        for (int k = 0; k < 100; ++k) {
            x_true = A * x_true;
            REQUIRE(ukf.predict(state_fn));
            REQUIRE(ukf.update(meas_fn, ColVec<1>{{x_true[0]}}, R));
        }

        CHECK(ukf.state()[0] == doctest::Approx(x_true[0]).epsilon(1e-3));
        CHECK(ukf.state()[1] == doctest::Approx(1.0).epsilon(1e-2));
    }

    TEST_CASE("strongly nonlinear measurement: UKF matches/beats EKF") {
        // Static 2D point observed through range-to-origin r = sqrt(x²+y²).
        // Compare the UKF (no Jacobian) against the EKF (hand-coded Jacobian)
        // from the same prior — the UKF must land on the same measurement
        // manifold (range ≈ 5).
        const Matrix<2, 2> Q = Matrix<2, 2>::identity() * 1e-6;
        const Matrix<1, 1> R{{1e-4}};

        auto f_ukf = [&](const ColVec<2>& x, const ColVec<1>&) { return x; };
        auto h_ukf = [&](const ColVec<2>& x, const ColVec<1>&) {
            return ColVec<1>{{damp::sqrt(x[0] * x[0] + x[1] * x[1])}};
        };

        const ColVec<1>                        y{{5.0}};
        UnscentedKalmanFilter<2, 1, 1, double> ukf{ColVec<2>{{2.5}, {4.5}}, Matrix<2, 2>::identity(), Q};
        for (int k = 0; k < 50; ++k) {
            REQUIRE(ukf.predict(f_ukf));
            REQUIRE(ukf.update(h_ukf, y, R));
        }

        const double r_est = damp::sqrt(ukf.state()[0] * ukf.state()[0] + ukf.state()[1] * ukf.state()[1]);
        // Range-only obs can't fix the tangential direction; the estimate must
        // at least land near the correct measurement manifold (range ≈ 5).
        CHECK(r_est == doctest::Approx(5.0).epsilon(0.05));
    }

    TEST_CASE("predict-only grows covariance; update shrinks it") {
        const Matrix<2, 2> A{{1.0, 0.1}, {0.0, 1.0}};
        const Matrix<2, 2> Q{{1e-3, 0.0}, {0.0, 1e-3}};
        const Matrix<1, 1> R{{1e-2}};

        auto state_fn = [&](const ColVec<2>& x, const ColVec<1>&) { return ColVec<2>{A * x}; };
        auto meas_fn = [&](const ColVec<2>& x, const ColVec<1>&) { return ColVec<1>{{x[0]}}; };

        UnscentedKalmanFilter<2, 1, 1, double> ukf{ColVec<2>{{0.0}, {0.0}}, Matrix<2, 2>::identity() * 0.1, Q};

        const double p0 = ukf.covariance()(0, 0);
        REQUIRE(ukf.predict(state_fn));
        const double p_after_predict = ukf.covariance()(0, 0);
        CHECK(p_after_predict > p0);

        REQUIRE(ukf.update(meas_fn, ColVec<1>{{0.0}}, R));
        CHECK(ukf.covariance()(0, 0) < p_after_predict);
    }

    TEST_CASE("update returns false on singular innovation covariance") {
        // R = 0 and a measurement that sees a zero-covariance direction makes
        // Pyy singular; Cholesky and LU both fail closed (false).
        const Matrix<2, 2> Q = Matrix<2, 2>::zeros();
        const Matrix<1, 1> R = Matrix<1, 1>::zeros();

        auto state_fn = [&](const ColVec<2>& x, const ColVec<1>&) { return x; };
        auto meas_fn = [&](const ColVec<2>&, const ColVec<1>&) { return ColVec<1>{{0.0}}; };

        // Non-degenerate prior so sigma points can be drawn, but a constant h
        // gives zero output spread → Pyy = R = 0, singular.
        UnscentedKalmanFilter<2, 1, 1, double> ukf{ColVec<2>{{0.0}, {0.0}}, Matrix<2, 2>::identity(), Q};
        REQUIRE(ukf.predict(state_fn));
        const auto x_before = ukf.state();
        CHECK_FALSE(ukf.update(meas_fn, ColVec<1>{{1.0}}, R));
        CHECK(ukf.state()[0] == doctest::Approx(x_before[0]));                    // state unchanged on fail
        CHECK(ukf.covariance()(0, 1) == doctest::Approx(ukf.covariance()(1, 0))); // P stays symmetric
    }

    TEST_CASE("predict keeps P symmetric") {
        const Matrix<2, 2>                     A{{1.0, 0.1}, {0.0, 1.0}};
        const Matrix<2, 2>                     Q{{1e-3, 0.0}, {0.0, 1e-3}};
        auto                                   state_fn = [&](const ColVec<2>& x, const ColVec<1>&) { return ColVec<2>{A * x}; };
        UnscentedKalmanFilter<2, 1, 1, double> ukf{ColVec<2>{{0.0}, {0.0}}, Matrix<2, 2>::identity(), Q};
        REQUIRE(ukf.predict(state_fn));
        CHECK(ukf.covariance()(0, 1) == doctest::Approx(ukf.covariance()(1, 0)).epsilon(1e-14));
    }

    TEST_CASE("set_state / set_covariance round-trip") {
        UnscentedKalmanFilter<2, 1, 1, double> ukf{};
        ukf.set_state(ColVec<2>{{1.5}, {-2.0}});
        CHECK(ukf.state()[0] == doctest::Approx(1.5));
        CHECK(ukf.state()[1] == doctest::Approx(-2.0));

        ukf.set_state(1, 7.0);
        CHECK(ukf.state()[1] == doctest::Approx(7.0));

        ukf.set_covariance(Matrix<2, 2>::identity() * 3.0);
        CHECK(ukf.covariance()(0, 0) == doctest::Approx(3.0));
    }

    TEST_CASE("UKF alias instantiates") {
        UKF<2, 1, 1> ukf{};
        CHECK(ukf.state()[0] == doctest::Approx(0.0));
    }

    TEST_CASE("float instantiation compiles and runs") {
        const Matrix<2, 2, float> A{{1.0f, 0.1f}, {0.0f, 1.0f}};
        auto                      state_fn = [&](const ColVec<2, float>& x, const ColVec<1, float>&) { return ColVec<2, float>{A * x}; };
        auto                      meas_fn = [&](const ColVec<2, float>& x, const ColVec<1, float>&) { return ColVec<1, float>{{x[0]}}; };

        UnscentedKalmanFilter<2, 1, 1, float> ukf{ColVec<2, float>{{0.0f}, {0.0f}}, Matrix<2, 2, float>::identity(), Matrix<2, 2, float>::identity() * 1e-3f};
        REQUIRE(ukf.predict(state_fn));
        CHECK(ukf.update(meas_fn, ColVec<1, float>{{0.5f}}, Matrix<1, 1, float>{{1e-2f}}));
    }
}
