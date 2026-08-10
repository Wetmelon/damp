// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/design/qp.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/estimation/mhe.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

// Deterministic "measurement noise": a fixed repeating pattern so the tests
// are reproducible while still exercising the estimator like real noise does.
constexpr damp::array<double, 12> noise_pattern{0.030, -0.050, 0.010, -0.040, 0.050, -0.020, 0.045, -0.035, 0.005, -0.045, 0.025, -0.015};
constexpr double                  noise(size_t k) {
    return noise_pattern[k % noise_pattern.size()];
}

} // namespace

TEST_CASE("example: tank-level estimator that cannot report a negative level") {
    // ------------------------------------------------------------------
    // The MHE workflow, end to end. Scenario: a draining tank whose level
    // is measured with a noisy sensor. Near empty, sensor noise makes a
    // Kalman filter report physically impossible negative levels; the MHE
    // gets the constraint x >= 0 and cannot.
    // ------------------------------------------------------------------

    // 1. Model the plant, exactly as for a Kalman filter design:
    //    level drains 10%/tick, inflow u adds directly; sensor reads level.
    constexpr StateSpace<1, 1, 1> tank{
        .A = {{0.90}},
        .B = {{0.10}},
        .C = {{1.0}},
        .Ts = 0.1,
    };

    // 2. Pick the noise covariances - the SAME numbers a Kalman design uses.
    //    Small Q: the drain model is trusted. R matches the sensor noise
    //    (std ~0.05 -> R = 0.0025).
    const Matrix<1, 1> Q{{1.0e-4}};
    const Matrix<1, 1> R{{2.5e-3}};

    // 3. State the physics the estimate must respect: a level is never
    //    negative. (Defaults are unbounded - this is the only extra input
    //    MHE needs over a Kalman filter.)
    design::MHEConstraints<1> limits{};
    limits.x_min = ColVec<1>{0.0};

    // 4. Synthesize with an 8-sample window and construct the runtime.
    const auto art = design::mhe<8>(tank, Q, R, limits);
    REQUIRE(art.success);
    MHE estimator{art};
    estimator.set_initial_state(ColVec<1>{0.5}); // prior for the first window

    // Reference: the equivalent (unconstrained) steady-state Kalman filter.
    StateSpace<1, 1, 1, double, 1, 1> tank_noise{};
    tank_noise.A = tank.A;
    tank_noise.B = tank.B;
    tank_noise.C = tank.C;
    tank_noise.Ts = tank.Ts;
    const auto kf_design = design::kalman(tank_noise, Q, R);
    REQUIRE(kf_design.success);
    KalmanFilter<1, 1, 1, double, 1, 1> kf{kf_design};
    kf.set_state(ColVec<1>{0.5});

    // 5. Run: truth drains from 0.5 toward empty; each tick both estimators
    //    receive the same noisy level reading and the same (zero) inflow.
    double x_true = 0.5;
    double kf_min = 1.0;
    double mhe_min = 1.0;
    double mhe_final = 0.0;
    for (size_t k = 0; k < 80; ++k) {
        const double u = 0.0;
        x_true = (0.90 * x_true) + (0.10 * u);
        const ColVec<1> y{x_true + noise(k)};

        const auto x_mhe = estimator.update(y, ColVec<1>{u});
        CHECK(estimator.last_status() == design::QPStatus::Success);

        kf.predict(ColVec<1>{u});
        kf.update(y, ColVec<1>{u});

        mhe_min = (x_mhe(0) < mhe_min) ? x_mhe(0) : mhe_min;
        kf_min = (kf.state()(0) < kf_min) ? kf.state()(0) : kf_min;
        mhe_final = x_mhe(0);
    }

    // 6. The payoff: the Kalman filter went negative near empty; the MHE
    //    stayed physical, and still converged to the (near-empty) truth.
    CHECK(kf_min < 0.0);
    CHECK(mhe_min >= -1e-9);
    CHECK(mhe_final == doctest::Approx(x_true).epsilon(0.05));
}

TEST_CASE("unconstrained MHE matches the steady-state Kalman filter") {
    // Roadmap #15 acceptance: with no active constraints, the window endpoint
    // reproduces the steady-state Kalman estimate (the arrival cost is the
    // steady-state Kalman covariance, so the two solve the same problem).
    constexpr double Ts = 0.1;

    StateSpace<2, 1, 1, double, 2, 1> plant{};
    plant.A = Matrix<2, 2>{{1.0, Ts}, {0.0, 1.0}};
    plant.B = Matrix<2, 1>{{0.5 * Ts * Ts}, {Ts}};
    plant.C = Matrix<1, 2>{{1.0, 0.0}};
    plant.Ts = Ts;

    const Matrix<2, 2> Q{{1.0e-4, 0.0}, {0.0, 1.0e-3}};
    const Matrix<1, 1> R{{2.5e-3}};

    const auto art = design::mhe<8>(plant, Q, R);
    REQUIRE(art.success);
    MHE estimator{art};

    const auto kf_design = design::kalman(plant, Q, R);
    REQUIRE(kf_design.success);
    KalmanFilter<2, 1, 1, double, 2, 1> kf{kf_design};

    ColVec<2> x_true{0.2, -0.1};
    for (size_t k = 0; k < 60; ++k) {
        const double u = 0.3 * noise(k + 3); // wiggle the input deterministically
        x_true = ColVec<2>{
            x_true(0) + (Ts * x_true(1)) + (0.5 * Ts * Ts * u),
            x_true(1) + (Ts * u),
        };
        const ColVec<1> y{x_true(0) + (0.3 * noise(k))};

        const auto x_mhe = estimator.update(y, ColVec<1>{u});
        kf.predict(ColVec<1>{u});
        kf.update(y, ColVec<1>{u});

        if (k >= 20) { // after the window has flushed its priming samples
            CHECK(x_mhe(0) == doctest::Approx(kf.state()(0)).epsilon(0.02));
            CHECK(x_mhe(1) == doctest::Approx(kf.state()(1)).epsilon(0.05));
        }
    }
}

TEST_CASE("MHE warm start settles to a few iterations per tick") {
    constexpr StateSpace<1, 1, 1> tank{
        .A = {{0.90}},
        .B = {{0.10}},
        .C = {{1.0}},
        .Ts = 0.1,
    };
    const auto art = design::mhe<8>(tank, Matrix<1, 1>{{1e-4}}, Matrix<1, 1>{{2.5e-3}});
    REQUIRE(art.success);
    MHE estimator{art};

    size_t settled_iterations = 0;
    for (size_t k = 0; k < 40; ++k) {
        (void)estimator.update(ColVec<1>{0.2 + (0.05 * noise(k))});
        if (k >= 30) {
            settled_iterations += estimator.last_iterations();
        }
    }
    // Unconstrained window at steady state: the QP is solved at its
    // unconstrained minimum, no active-set changes at all.
    CHECK(settled_iterations == 0);
}

TEST_CASE("mhe rejects invalid inputs") {
    StateSpace<1, 1, 1> continuous{};
    continuous.A = Matrix<1, 1>{{0.9}};
    continuous.B = Matrix<1, 1>{{0.1}};
    continuous.C = Matrix<1, 1>{{1.0}};
    continuous.Ts = 0.0; // not discrete
    CHECK(!design::mhe<4>(continuous, Matrix<1, 1>{{1e-4}}, Matrix<1, 1>{{1e-3}}).success);

    StateSpace<1, 1, 1> feedthrough = continuous;
    feedthrough.Ts = 0.1;
    feedthrough.D = Matrix<1, 1>{{1.0}};
    CHECK(!design::mhe<4>(feedthrough, Matrix<1, 1>{{1e-4}}, Matrix<1, 1>{{1e-3}}).success);

    StateSpace<1, 1, 1> ok = feedthrough;
    ok.D = Matrix<1, 1>{};
    CHECK(!design::mhe<4>(ok, Matrix<1, 1>{{0.0}}, Matrix<1, 1>{{1e-3}}).success); // singular Q
}

TEST_CASE("mhe is constexpr") {
    constexpr StateSpace<1, 1, 1> tank{
        .A = {{0.90}},
        .B = {{0.10}},
        .C = {{1.0}},
        .Ts = 0.1,
    };
    constexpr auto art = design::mhe<3>(tank, Matrix<1, 1>{{1e-4}}, Matrix<1, 1>{{2.5e-3}});
    static_assert(art.success);
}

TEST_CASE("MHE deploys in float") {
    constexpr StateSpace<1, 1, 1> tank{
        .A = {{0.90}},
        .B = {{0.10}},
        .C = {{1.0}},
        .Ts = 0.1,
    };
    design::MHEConstraints<1> limits{};
    limits.x_min = ColVec<1>{0.0};

    const auto art = design::mhe<6>(tank, Matrix<1, 1>{{1e-4}}, Matrix<1, 1>{{2.5e-3}}, limits);
    REQUIRE(art.success);

    MHE   estimator{art.as<float>()};
    float x_true = 0.4F;
    float x_hat = 0.0F;
    for (size_t k = 0; k < 60; ++k) {
        x_true *= 0.90F;
        const auto e = estimator.update(ColVec<1, float>{x_true + (0.05F * static_cast<float>(noise(k)))});
        x_hat = e(0);
        CHECK(x_hat >= -1e-6F);
    }
    CHECK(x_hat == doctest::Approx(x_true).epsilon(0.1));
}
