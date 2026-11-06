// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)


#include <type_traits>

#include "damp/estimation/ukf.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

// The Unscented Kalman Filter propagates sigma points through the exact
// nonlinear maps f/h — no Jacobian. The model is the continuous Dynamics
// ẋ = f(x, u), a struct or a lambda the filter holds, carried over the sample time by
// one RK4 step per sigma point, and each reading comes with its
// Measurement y = h(x, u) and R. These tests cover: (1) the linear case, where the UKF
// must reduce to a linear Kalman filter; (2) a strongly nonlinear measurement
// where the UKF should at least match the EKF; (3) covariance grow/shrink
// behaviour; (4) the singular-Pyy failure path; and (5) float instantiation.
//
// @see Wan & van der Merwe (2000); Julier & Uhlmann (2004); Simon (2006) §14.3.

namespace {

// Constant velocity: ṗ = v, v̇ = 0. One RK4 step is exact (A is nilpotent).
struct ConstantVelocity {
    [[nodiscard]] ColVec<2> operator()(const ColVec<2>& x, const ColVec<1>& /*u*/) const { return ColVec<2>{x[1], 0.0}; }
};

// A static point: ẋ = 0.
struct Stationary {
    [[nodiscard]] ColVec<2> operator()(const ColVec<2>& /*x*/, const ColVec<1>& /*u*/) const { return ColVec<2>{}; }
};

struct Position {
    [[nodiscard]] ColVec<1> operator()(const ColVec<2>& x, const ColVec<1>& /*u*/) const { return ColVec<1>{x[0]}; }
};

struct Range {
    [[nodiscard]] ColVec<1> operator()(const ColVec<2>& x, const ColVec<1>& /*u*/) const { return ColVec<1>{damp::sqrt((x[0] * x[0]) + (x[1] * x[1]))}; }
};

struct Blind {
    [[nodiscard]] ColVec<1> operator()(const ColVec<2>& /*x*/, const ColVec<1>& /*u*/) const { return ColVec<1>{0.0}; }
};

struct ConstantVelocityF {
    [[nodiscard]] ColVec<2, float> operator()(const ColVec<2, float>& x, const ColVec<1, float>& /*u*/) const { return ColVec<2, float>{x[1], 0.0F}; }
};

struct PositionF {
    [[nodiscard]] ColVec<1, float> operator()(const ColVec<2, float>& x, const ColVec<1, float>& /*u*/) const { return ColVec<1, float>{x[0]}; }
};

const ColVec<1> kNoInput{};

// The model as a captureless lambda: the filter deduces from it, the scalar from its ColVec.
constexpr auto                         kConstantVelocityLambda = [](const ColVec<2, float>& x, const ColVec<1, float>& /*u*/) { return ColVec<2, float>{x[1], 0.0F}; };
static constinit UnscentedKalmanFilter kLambdaUKF{kConstantVelocityLambda, {.Ts = 0.1F, .Q = Matrix<2, 2, float>::identity() * 1e-3F}};
static_assert(std::is_same_v<decltype(kLambdaUKF), UnscentedKalmanFilter<std::remove_cvref_t<decltype(kConstantVelocityLambda)>, float>>);

} // namespace

TEST_SUITE("Unscented Kalman Filter") {
    TEST_CASE("linear dynamics + linear measurement converges to truth") {
        // Constant-velocity model with a linear position measurement. On a
        // linear f/h the unscented transform is exact, so the UKF must behave
        // like a textbook KF and recover the unmeasured velocity.
        const double       Ts = 0.1;
        const Matrix<2, 2> A{{1.0, Ts}, {0.0, 1.0}};

        UnscentedKalmanFilter<ConstantVelocity, double> ukf{ConstantVelocity{}, {.Ts = Ts, .Q = Matrix<2, 2>{{1e-4, 0.0}, {0.0, 1e-4}}}};
        const Matrix<1, 1>                              R{1e-2};

        ColVec<2> x_true{{0.0}, {1.0}}; // 1 m/s
        for (int k = 0; k < 100; ++k) {
            x_true = A * x_true;
            REQUIRE(ukf.predict(kNoInput));
            REQUIRE(ukf.update(Position{}, ColVec<1>{{x_true[0]}}, R));
        }

        CHECK(ukf.state()[0] == doctest::Approx(x_true[0]).epsilon(1e-3));
        CHECK(ukf.state()[1] == doctest::Approx(1.0).epsilon(1e-2));
    }

    TEST_CASE("strongly nonlinear measurement: UKF matches/beats EKF") {
        // Static 2D point observed through range-to-origin r = sqrt(x²+y²).
        // The UKF (no Jacobian) must land on the same measurement manifold as
        // the EKF does from the same prior (range ≈ 5).
        const ColVec<1>                           y{{5.0}};
        UnscentedKalmanFilter<Stationary, double> ukf{Stationary{}, {.Ts = 0.01, .Q = Matrix<2, 2>::identity() * 1e-6, .x0 = ColVec<2>{{2.5}, {4.5}}}};
        for (int k = 0; k < 50; ++k) {
            REQUIRE(ukf.predict(kNoInput));
            REQUIRE(ukf.update(Range{}, y, 1e-4));
        }

        const double r_est = damp::sqrt((ukf.state()[0] * ukf.state()[0]) + (ukf.state()[1] * ukf.state()[1]));
        // Range-only obs can't fix the tangential direction; the estimate must
        // at least land near the correct measurement manifold (range ≈ 5).
        CHECK(r_est == doctest::Approx(5.0).epsilon(0.05));
    }

    TEST_CASE("predict-only grows covariance; update shrinks it") {
        UnscentedKalmanFilter<ConstantVelocity, double> ukf{ConstantVelocity{}, {.Ts = 0.1, .Q = Matrix<2, 2>::identity() * 1e-3, .P0 = Matrix<2, 2>::identity() * 0.1}};

        const double p0 = ukf.covariance()(0, 0);
        REQUIRE(ukf.predict(kNoInput));
        const double p_after_predict = ukf.covariance()(0, 0);
        CHECK(p_after_predict > p0);

        REQUIRE(ukf.update(Position{}, ColVec<1>{{0.0}}, 1e-2));
        CHECK(ukf.covariance()(0, 0) < p_after_predict);
    }

    TEST_CASE("update returns false on singular innovation covariance") {
        // R = 0 and a constant h: the sigma points can be drawn, but the output spread is zero,
        // so Pyy = R = 0 has no Cholesky factor and the update fails closed (false).
        UnscentedKalmanFilter<Stationary, double> ukf{Stationary{}, {.Ts = 0.1, .Q = Matrix<2, 2>::zeros()}};
        REQUIRE(ukf.predict(kNoInput));
        const auto x_before = ukf.state();
        CHECK_FALSE(ukf.update(Blind{}, ColVec<1>{{1.0}}, 0.0));
        CHECK(ukf.state()[0] == doctest::Approx(x_before[0]));                    // state unchanged on fail
        CHECK(ukf.covariance()(0, 1) == doctest::Approx(ukf.covariance()(1, 0))); // P stays symmetric
    }

    TEST_CASE("a semidefinite P (a state known exactly) draws sigma points; an indefinite one does not") {
        //! P0 = diag(0, 1): Cholesky refused it, predict returned before adding Q, and the filter
        //! stayed frozen at x0.
        UnscentedKalmanFilter<ConstantVelocity, double> ukf{ConstantVelocity{}, {.Ts = 0.1, .Q = Matrix<2, 2>::identity() * 1e-3, .P0 = Matrix<2, 2>{{0.0, 0.0}, {0.0, 1.0}}}};
        for (int k = 0; k < 20; ++k) {
            REQUIRE(ukf.predict(kNoInput));
            REQUIRE(ukf.update(Position{}, ColVec<1>{{0.5}}, 1e-2));
        }
        CHECK(ukf.covariance()(0, 0) > 0.0);
        CHECK(ukf.state()[0] == doctest::Approx(0.5).epsilon(0.05));

        UnscentedKalmanFilter<ConstantVelocity, double> indefinite{ConstantVelocity{}, {.Ts = 0.1, .Q = Matrix<2, 2>::identity() * 1e-3, .P0 = Matrix<2, 2>{{-1.0, 0.0}, {0.0, 1.0}}}};
        CHECK_FALSE(indefinite.predict(kNoInput));
        UnscentedKalmanFilter<ConstantVelocity, double> coupled{ConstantVelocity{}, {.Ts = 0.1, .Q = Matrix<2, 2>::identity() * 1e-3, .P0 = Matrix<2, 2>{{0.0, 1.0}, {1.0, 1.0}}}};
        CHECK_FALSE(coupled.predict(kNoInput)); // a zero variance with a nonzero covariance
    }

    TEST_CASE("predict keeps P symmetric") {
        UnscentedKalmanFilter<ConstantVelocity, double> ukf{ConstantVelocity{}, {.Ts = 0.1, .Q = Matrix<2, 2>::identity() * 1e-3}};
        REQUIRE(ukf.predict(kNoInput));
        CHECK(ukf.covariance()(0, 1) == doctest::Approx(ukf.covariance()(1, 0)).epsilon(1e-14));
    }

    TEST_CASE("set_state / reset") {
        UnscentedKalmanFilter<Stationary, double> ukf{};
        ukf.set_state(ColVec<2>{{1.5}, {-2.0}});
        CHECK(ukf.state()[0] == doctest::Approx(1.5));
        CHECK(ukf.state()[1] == doctest::Approx(-2.0));

        ukf.reset(ColVec<2>{0.5, 7.0}, Matrix<2, 2>::identity() * 3.0);
        CHECK(ukf.state()[1] == doctest::Approx(7.0));
        CHECK(ukf.covariance()(0, 0) == doctest::Approx(3.0));
    }

    TEST_CASE("a lambda model in static constinit, with a lambda and a struct measurement") {
        constexpr auto by_lambda = [](const ColVec<2, float>& x, const ColVec<1, float>& /*u*/) { return ColVec<1, float>{x[0]}; };
        // P = I, R = 1: a position reading 0.5 moves p by 0.5/2 = 0.25 whichever form h takes.
        auto a = kLambdaUKF;
        auto b = kLambdaUKF;
        REQUIRE(a.update(by_lambda, ColVec<1, float>{0.5F}, 1.0F));
        REQUIRE(b.update(PositionF{}, ColVec<1, float>{0.5F}, 1.0F));
        CHECK(a.state()[0] == doctest::Approx(0.25).epsilon(1e-5));
        CHECK(b.state()[0] == doctest::Approx(0.25).epsilon(1e-5));
    }

    TEST_CASE("the default scalar is the model's") {
        static_assert(std::is_same_v<UnscentedKalmanFilter<ConstantVelocityF>, UnscentedKalmanFilter<ConstantVelocityF, float>>);
        static_assert(std::is_same_v<UnscentedKalmanFilter<ConstantVelocity>, UnscentedKalmanFilter<ConstantVelocity, double>>);
        UnscentedKalmanFilter<ConstantVelocityF> ukf{ConstantVelocityF{}, {.Ts = 0.1F, .Q = Matrix<2, 2, float>::identity() * 1e-3F}};
        REQUIRE(ukf.predict({}));
        CHECK(ukf.update(PositionF{}, ColVec<1, float>{{0.5F}}, 1e-2F));
        CHECK(ukf.sample_time() == doctest::Approx(0.1));
    }

    TEST_CASE("UKF alias names the same type as UnscentedKalmanFilter") {
        CHECK(std::is_same_v<UKF<ConstantVelocityF>, UnscentedKalmanFilter<ConstantVelocityF, float>>);
        UKF<ConstantVelocityF> ukf{};
        CHECK(ukf.state()[0] == doctest::Approx(0.0));
    }
}
