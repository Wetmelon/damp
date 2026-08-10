// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>

#include "damp/estimation/ekf.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

// The Extended Kalman Filter handles nonlinear dynamics/measurement by
// linearizing about the current estimate. The user supplies callables that
// return both the function value and its Jacobian in one shot. These tests
// cover: (1) the linear case, where the EKF must reduce exactly to a linear
// Kalman filter; (2) a genuinely nonlinear measurement (range/bearing-style);
// and (3) the singular-innovation failure path.
//
// @see "Optimal State Estimation" (Simon, 2006), §13.3.

TEST_SUITE("Extended Kalman Filter") {
    TEST_CASE("linear dynamics + linear measurement converges to truth") {
        // x[k+1] = A x[k] (constant-velocity), y = position. A linear f/h with
        // constant Jacobians is the degenerate EKF — it should behave like a KF.
        const Matrix<2, 2> A{{1.0, 0.1}, {0.0, 1.0}};
        const Matrix<2, 2> Q{{1e-4, 0.0}, {0.0, 1e-4}};
        const Matrix<1, 1> R{{1e-2}};

        auto state_fn = [&](const ColVec<2>& x, const ColVec<1>&) {
            return StateJacobian<double, 2>{A * x, A, Matrix<2, 2>::identity()};
        };
        auto meas_fn = [&](const ColVec<2>& x, const ColVec<1>&) {
            return MeasJacobian<double, 1, 2>{ColVec<1>{{x[0]}}, Matrix<1, 2>{{1.0, 0.0}}, Matrix<1, 1>::identity()};
        };

        ExtendedKalmanFilter<2, 1, 1, double> ekf{ColVec<2>{{0.0}, {0.0}}, Matrix<2, 2>::identity(), Q};

        // Truth propagates with the same A; feed noiseless position measurements.
        ColVec<2> x_true{{0.0}, {1.0}}; // moving at 1 m/s
        for (int k = 0; k < 100; ++k) {
            x_true = A * x_true;
            ekf.predict(state_fn);
            const ColVec<1> y{{x_true[0]}};
            REQUIRE(ekf.update(meas_fn, y, R));
        }

        // Estimator must recover both the measured position and the *unmeasured*
        // velocity (observable through the dynamics coupling).
        CHECK(ekf.state()[0] == doctest::Approx(x_true[0]).epsilon(1e-3));
        CHECK(ekf.state()[1] == doctest::Approx(1.0).epsilon(1e-2));
    }

    TEST_CASE("nonlinear measurement (range to origin) is tracked") {
        // State = 2D position, static. Measurement is range r = sqrt(x² + y²),
        // a nonlinear h whose Jacobian is [x/r, y/r]. The EKF must localize the
        // point from range alone given a decent prior.
        const Matrix<2, 2> Q = Matrix<2, 2>::identity() * 1e-6;
        const Matrix<1, 1> R{{1e-4}};

        auto state_fn = [&](const ColVec<2>& x, const ColVec<1>&) {
            // Static point: f(x) = x, F = I.
            return StateJacobian<double, 2>{x, Matrix<2, 2>::identity(), Matrix<2, 2>::identity()};
        };
        auto meas_fn = [&](const ColVec<2>& x, const ColVec<1>&) {
            const double r = damp::sqrt(x[0] * x[0] + x[1] * x[1]);
            const double rs = r > 1e-9 ? r : 1e-9;
            return MeasJacobian<double, 1, 2>{
                ColVec<1>{{r}},
                Matrix<1, 2>{{x[0] / rs, x[1] / rs}},
                Matrix<1, 1>::identity()
            };
        };

        const ColVec<2>                       x_true{{3.0}, {4.0}}; // true range = 5
        ExtendedKalmanFilter<2, 1, 1, double> ekf{ColVec<2>{{2.5}, {4.5}}, Matrix<2, 2>::identity(), Q};

        const ColVec<1> y{{5.0}};
        for (int k = 0; k < 50; ++k) {
            ekf.predict(state_fn);
            REQUIRE(ekf.update(meas_fn, y, R));
        }

        // Range-only obs can't fix the tangential direction, but the estimate
        // must at least land on the correct measurement manifold (range ≈ 5).
        const double r_est = damp::sqrt(ekf.state()[0] * ekf.state()[0] + ekf.state()[1] * ekf.state()[1]);
        CHECK(r_est == doctest::Approx(5.0).epsilon(1e-2));
    }

    TEST_CASE("predict-only grows covariance; update shrinks it") {
        const Matrix<2, 2> A{{1.0, 0.1}, {0.0, 1.0}};
        const Matrix<2, 2> Q{{1e-3, 0.0}, {0.0, 1e-3}};
        const Matrix<1, 1> R{{1e-2}};

        auto state_fn = [&](const ColVec<2>& x, const ColVec<1>&) {
            return StateJacobian<double, 2>{A * x, A, Matrix<2, 2>::identity()};
        };
        auto meas_fn = [&](const ColVec<2>& x, const ColVec<1>&) {
            return MeasJacobian<double, 1, 2>{ColVec<1>{{x[0]}}, Matrix<1, 2>{{1.0, 0.0}}, Matrix<1, 1>::identity()};
        };

        ExtendedKalmanFilter<2, 1, 1, double> ekf{ColVec<2>{{0.0}, {0.0}}, Matrix<2, 2>::identity() * 0.1, Q};

        const double p0 = ekf.covariance()(0, 0);
        ekf.predict(state_fn);
        const double p_after_predict = ekf.covariance()(0, 0);
        CHECK(p_after_predict > p0); // prediction adds uncertainty

        REQUIRE(ekf.update(meas_fn, ColVec<1>{{0.0}}, R));
        const double p_after_update = ekf.covariance()(0, 0);
        CHECK(p_after_update < p_after_predict); // measurement removes uncertainty
    }

    TEST_CASE("update returns false on singular innovation covariance") {
        // R = 0 and a measurement Jacobian that sees a zero-covariance direction
        // can make S singular; the Cholesky solve must fail gracefully (false),
        // not crash.
        const Matrix<2, 2> Q = Matrix<2, 2>::zeros();
        const Matrix<1, 1> R = Matrix<1, 1>::zeros();

        auto state_fn = [&](const ColVec<2>& x, const ColVec<1>&) {
            return StateJacobian<double, 2>{x, Matrix<2, 2>::identity(), Matrix<2, 2>::identity()};
        };
        // H = 0 → S = H P Hᵀ + R = 0, singular.
        auto meas_fn = [&](const ColVec<2>&, const ColVec<1>&) {
            return MeasJacobian<double, 1, 2>{ColVec<1>{{0.0}}, Matrix<1, 2>{{0.0, 0.0}}, Matrix<1, 1>::identity()};
        };

        ExtendedKalmanFilter<2, 1, 1, double> ekf{ColVec<2>{{0.0}, {0.0}}, Matrix<2, 2>::zeros(), Q};
        ekf.predict(state_fn);
        CHECK_FALSE(ekf.update(meas_fn, ColVec<1>{{1.0}}, R));
    }

    TEST_CASE("sequential scalar updates with inter-measurement state clamping") {
        // The intended workflow for constrained estimation: instantiate with
        // NY == 1, fold in measurements one scalar at a time, and clamp the
        // affected state between updates via set_state. Clamping a fused vector
        // update can't enforce the constraint each measurement sees; clamping
        // between scalar updates can. Here state[0] is physically bounded to
        // [0, 1] (e.g. a duty cycle / SoC) and must never leave that box even
        // though a noisy measurement pushes the raw estimate past 1.
        const Matrix<2, 2> Q = Matrix<2, 2>::identity() * 1e-3;

        auto state_fn = [&](const ColVec<2>& x, const ColVec<1>&) {
            return StateJacobian<double, 2>{x, Matrix<2, 2>::identity(), Matrix<2, 2>::identity()};
        };
        // Two separate scalar measurements of the two states.
        auto meas0 = [&](const ColVec<2>& x, const ColVec<1>&) {
            return MeasJacobian<double, 1, 2>{ColVec<1>{{x[0]}}, Matrix<1, 2>{{1.0, 0.0}}, Matrix<1, 1>::identity()};
        };
        auto meas1 = [&](const ColVec<2>& x, const ColVec<1>&) {
            return MeasJacobian<double, 1, 2>{ColVec<1>{{x[1]}}, Matrix<1, 2>{{0.0, 1.0}}, Matrix<1, 1>::identity()};
        };

        const Matrix<1, 1>                    R{{1e-2}};
        ExtendedKalmanFilter<2, 1, 1, double> ekf{ColVec<2>{{0.5}, {0.5}}, Matrix<2, 2>::identity(), Q};

        for (int k = 0; k < 200; ++k) {
            ekf.predict(state_fn);

            // Scalar update #1: measurement of state[0] sits above its physical
            // bound (1.2 > 1.0). Fuse it, then clamp before the next update.
            REQUIRE(ekf.update(meas0, ColVec<1>{{1.2}}, R));
            if (ekf.state()[0] > 1.0) {
                ekf.set_state(0, 1.0);
            } else if (ekf.state()[0] < 0.0) {
                ekf.set_state(0, 0.0);
            }

            // The clamped estimate is what the second scalar update sees.
            CHECK(ekf.state()[0] <= 1.0 + 1e-12);
            REQUIRE(ekf.update(meas1, ColVec<1>{{0.3}}, R));
        }

        CHECK(ekf.state()[0] == doctest::Approx(1.0)); // pinned at the bound
        CHECK(ekf.state()[1] == doctest::Approx(0.3).epsilon(0.05));
    }

    TEST_CASE("set_state / set_covariance round-trip") {
        ExtendedKalmanFilter<2, 1, 1, double> ekf{};
        ekf.set_state(ColVec<2>{{1.5}, {-2.0}});
        CHECK(ekf.state()[0] == doctest::Approx(1.5));
        CHECK(ekf.state()[1] == doctest::Approx(-2.0));

        ekf.set_state(1, 7.0);
        CHECK(ekf.state()[1] == doctest::Approx(7.0));

        ekf.set_covariance(Matrix<2, 2>::identity() * 3.0);
        CHECK(ekf.covariance()(0, 0) == doctest::Approx(3.0));
        CHECK(ekf.covariance()(1, 1) == doctest::Approx(3.0));
    }

    TEST_CASE("float instantiation compiles and runs") {
        const Matrix<2, 2, float> A{{1.0f, 0.1f}, {0.0f, 1.0f}};
        auto                      state_fn = [&](const ColVec<2, float>& x, const ColVec<1, float>&) {
            return StateJacobian<float, 2>{A * x, A, Matrix<2, 2, float>::identity()};
        };
        auto meas_fn = [&](const ColVec<2, float>& x, const ColVec<1, float>&) {
            return MeasJacobian<float, 1, 2>{ColVec<1, float>{{x[0]}}, Matrix<1, 2, float>{{1.0f, 0.0f}}, Matrix<1, 1, float>::identity()};
        };

        ExtendedKalmanFilter<2, 1, 1, float> ekf{ColVec<2, float>{{0.0f}, {0.0f}}, Matrix<2, 2, float>::identity(), Matrix<2, 2, float>::identity() * 1e-3f};
        ekf.predict(state_fn);
        CHECK(ekf.update(meas_fn, ColVec<1, float>{{0.5f}}, Matrix<1, 1, float>{{1e-2f}}));
    }
}
