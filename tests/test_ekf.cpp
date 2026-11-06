// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)


#include <cstddef>
#include <type_traits>

#include "damp/concepts.hpp"
#include "damp/estimation/ekf.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

// The Extended Kalman Filter handles nonlinear dynamics/measurement by
// linearizing about the current estimate. The model is the continuous
// Dynamics ẋ = f(x, u), a struct or a lambda the filter holds, and each reading
// comes with its Measurement y = h(x, u) and R; the filter integrates f
// over its sample time with one RK4 step, and takes the model's own Jacobians
// when it provides them, central differences otherwise. These tests cover:
// (1) the linear case, where the EKF must reduce to the linear Kalman filter of
// the RK4-discretized plant; (2) a genuinely nonlinear measurement
// (range-to-origin) with analytic and numerical Jacobians; (3) two readings of
// different sizes on one filter; and (4) the singular-innovation failure path.
//
// @see "Optimal State Estimation" (Simon, 2006), §13.3.

namespace {

// Constant velocity: ṗ = v, v̇ = 0. One RK4 step is exact (A is nilpotent).
struct ConstantVelocity {
    [[nodiscard]] ColVec<2> operator()(const ColVec<2>& x, const ColVec<1>& /*u*/) const { return ColVec<2>{x[1], 0.0}; }
};

// A static point: ẋ = 0.
struct Stationary {
    [[nodiscard]] ColVec<2> operator()(const ColVec<2>& /*x*/, const ColVec<1>& /*u*/) const { return ColVec<2>{}; }
};

// Damped oscillator ẋ = A x + B u.
struct Oscillator {
    Matrix<2, 2> A{{0.0, 1.0}, {-4.0, -0.4}};
    Matrix<2, 1> B{{0.0}, {1.0}};

    [[nodiscard]] ColVec<2> operator()(const ColVec<2>& x, const ColVec<1>& u) const { return (A * x) + (B * u); }
};

// One state measured, by central differences.
template<size_t I>
struct Component {
    [[nodiscard]] ColVec<1> operator()(const ColVec<2>& x, const ColVec<1>& /*u*/) const { return ColVec<1>{x[I]}; }
};

// Both states measured.
struct BothStates {
    [[nodiscard]] ColVec<2> operator()(const ColVec<2>& x, const ColVec<1>& /*u*/) const { return x; }
};

// A constant h: its central-difference Jacobian is C = 0.
struct Blind {
    [[nodiscard]] ColVec<1> operator()(const ColVec<2>& /*x*/, const ColVec<1>& /*u*/) const { return ColVec<1>{0.0}; }
};

// Range to the origin r = √(x² + y²), with its Jacobian [x/r, y/r].
struct RangeToOrigin {
    [[nodiscard]] ColVec<1> operator()(const ColVec<2>& x, const ColVec<1>& /*u*/) const {
        return ColVec<1>{damp::sqrt((x[0] * x[0]) + (x[1] * x[1]))};
    }
    [[nodiscard]] MeasurementJacobian<2, 1, 1> jacobian(const ColVec<2>& x, const ColVec<1>& u) const {
        const double r = (*this)(x, u)[0];
        const double rs = r > 1e-9 ? r : 1e-9;
        return {.C = Matrix<1, 2>{{x[0] / rs, x[1] / rs}}};
    }
};
static_assert(AnalyticMeasurement<RangeToOrigin, 2, 1, 1, double>);

// The same range without the Jacobian: central differences.
struct RangeNumeric {
    [[nodiscard]] ColVec<1> operator()(const ColVec<2>& x, const ColVec<1>& u) const { return RangeToOrigin{}(x, u); }
};
static_assert(!AnalyticMeasurement<RangeNumeric, 2, 1, 1, double>);

struct ConstantVelocityF {
    [[nodiscard]] ColVec<2, float> operator()(const ColVec<2, float>& x, const ColVec<1, float>& /*u*/) const { return ColVec<2, float>{x[1], 0.0F}; }
};

struct PositionF {
    [[nodiscard]] ColVec<1, float> operator()(const ColVec<2, float>& x, const ColVec<1, float>& /*u*/) const { return ColVec<1, float>{x[0]}; }
};

const ColVec<1> kNoInput{};

// Magnetic levitation, ẍ = g − (k/m)·(i/x)², the gap in metres (per_metre = 1) or millimetres
// (per_metre = 1000). At the equilibrium current ∂ẍ/∂x = 2g/x, near 3924 /s² at a 5 mm gap in
// either unit.
struct MaglevF {
    float i{};
    float per_metre{1.0F};

    [[nodiscard]] ColVec<2, float> operator()(const ColVec<2, float>& x, const ColVec<1, float>& /*u*/) const {
        const float ratio = i * per_metre / x[0]; // i / gap in metres
        return ColVec<2, float>{x[1], per_metre * (9.81F - (5e-4F * ratio * ratio))};
    }
};

// The same plant with its Jacobian, ∂ẍ/∂x = 2(k/m)·(i/x)²/x.
struct MaglevAnalyticF : MaglevF {
    [[nodiscard]] DynamicsJacobian<2, 1, float> jacobian(const ColVec<2, float>& x, const ColVec<1, float>& /*u*/) const {
        const float ratio = i * per_metre / x[0];
        return {.A = Matrix<2, 2, float>{{0.0F, 1.0F}, {per_metre * 2.0F * 5e-4F * ratio * ratio / x[0], 0.0F}}};
    }
};
static_assert(AnalyticDynamics<MaglevAnalyticF, 2, 1, float>);

// The pendulum of ExtendedKalmanFilter's documentation, θ from the upright: θ̈ = 9.81·sin θ − 0.5·θ̇ + u.
struct PendulumF {
    [[nodiscard]] ColVec<2, float> operator()(const ColVec<2, float>& x, const ColVec<1, float>& u) const {
        return ColVec<2, float>{x[1], (9.81F * damp::sin(x[0])) - (0.5F * x[1]) + u[0]};
    }
};

// A heading and a gyro bias: θ̇ = ω − b with ω the gyro reading (the input), ḃ = 0; ∂θ̇/∂b = −1.
// The volatile store rounds ω − b to float as a model compiled apart would; inlined under
// -ffast-math, (ω − b⁺) − (ω − b⁻) would fold to b⁻ − b⁺ and hide the step's roundoff.
struct HeadingBiasF {
    [[nodiscard]] ColVec<2, float> operator()(const ColVec<2, float>& x, const ColVec<1, float>& u) const {
        volatile float rate = u[0] - x[1];
        return ColVec<2, float>{rate, 0.0F};
    }
};

// Φ P₀ Φᵀ, Φ the RK4 stability function of A·Ts: one predict's covariance with Q = 0.
Matrix<2, 2> rk4_covariance(const Matrix<2, 2>& A, double Ts, const Matrix<2, 2>& P0) {
    const Matrix<2, 2> M = A * Ts;
    const Matrix<2, 2> M2 = M * M;
    const Matrix<2, 2> M3 = M2 * M;
    const Matrix<2, 2> Phi = Matrix<2, 2>::identity() + M + (M2 * 0.5) + (M3 * (1.0 / 6.0)) + ((M3 * M) * (1.0 / 24.0));
    return Phi * P0 * Phi.transpose();
}

// The model as a captureless lambda: the filter deduces from it, the scalar from its ColVec.
constexpr auto                        kConstantVelocityLambda = [](const ColVec<2, float>& x, const ColVec<1, float>& /*u*/) { return ColVec<2, float>{x[1], 0.0F}; };
static constinit ExtendedKalmanFilter kLambdaEKF{kConstantVelocityLambda, {.Ts = 0.1F, .Q = Matrix<2, 2, float>::identity() * 1e-3F}};
static_assert(std::is_same_v<decltype(kLambdaEKF), ExtendedKalmanFilter<std::remove_cvref_t<decltype(kConstantVelocityLambda)>, float>>);

// A struct model deduces too, and its scalar is the default.
static constinit ExtendedKalmanFilter kStructEKF{ConstantVelocity{}, {.Ts = 0.1, .Q = Matrix<2, 2>::identity() * 1e-3}};
static_assert(std::is_same_v<decltype(kStructEKF), ExtendedKalmanFilter<ConstantVelocity, double>>);

// A struct member names the model type.
struct HoldsFilter {
    ExtendedKalmanFilter<ConstantVelocityF> ekf_{ConstantVelocityF{}, {.Ts = 0.1F}};
};

} // namespace

TEST_SUITE("Extended Kalman Filter") {
    TEST_CASE("linear dynamics + linear measurement converges to truth") {
        // Constant velocity with a position measurement: a linear f/h is the
        // degenerate EKF — it should behave like a KF.
        const double       Ts = 0.1;
        const Matrix<2, 2> A{{1.0, Ts}, {0.0, 1.0}};
        const Matrix<2, 2> Q{{1e-4, 0.0}, {0.0, 1e-4}};
        const Matrix<1, 1> R{1e-2};

        ExtendedKalmanFilter<ConstantVelocity, double> ekf{ConstantVelocity{}, {.Ts = Ts, .Q = Q}};

        // Truth propagates with the same A; feed noiseless position measurements.
        ColVec<2> x_true{{0.0}, {1.0}}; // moving at 1 m/s
        for (int k = 0; k < 100; ++k) {
            x_true = A * x_true;
            ekf.predict(kNoInput);
            const ColVec<1> y{{x_true[0]}};
            REQUIRE(ekf.update(Component<0>{}, y, R));
        }

        // Estimator must recover both the measured position and the *unmeasured*
        // velocity (observable through the dynamics coupling).
        CHECK(ekf.state()[0] == doctest::Approx(x_true[0]).epsilon(1e-3));
        CHECK(ekf.state()[1] == doctest::Approx(1.0).epsilon(1e-2));
    }

    TEST_CASE("on a linear model predict is the Kalman time update of the RK4 discretization") {
        // Damped oscillator ẋ = A x + B u. One RK4 step with u held is
        // x⁺ = Φ x + Ts·S·B u, Φ = I + M + M²/2 + M³/6 + M⁴/24 and
        // S = I + M/2 + M²/6 + M³/24 (M = A·Ts): the mean and the covariance must
        // both propagate with that Φ.
        const Oscillator   osc{};
        const double       Ts = 0.05;
        const Matrix<2, 2> Q = Matrix<2, 2>::identity() * 1e-3;
        const Matrix<2, 2> P0{{0.5, 0.1}, {0.1, 0.2}};
        const ColVec<2>    x0{0.3, -0.2};
        const ColVec<1>    u{0.7};

        const Matrix<2, 2> I = Matrix<2, 2>::identity();
        const Matrix<2, 2> M = osc.A * Ts;
        const Matrix<2, 2> M2 = M * M;
        const Matrix<2, 2> M3 = M2 * M;
        const Matrix<2, 2> Phi = I + M + (M2 * 0.5) + (M3 * (1.0 / 6.0)) + ((M3 * M) * (1.0 / 24.0));
        const Matrix<2, 2> S = I + (M * 0.5) + (M2 * (1.0 / 6.0)) + (M3 * (1.0 / 24.0));
        const ColVec<2>    x_expected = (Phi * x0) + ((S * (osc.B * u)) * Ts);
        const Matrix<2, 2> P_expected = (Phi * P0 * Phi.transpose()) + Q;

        ExtendedKalmanFilter<Oscillator, double> ekf{osc, {.Ts = Ts, .Q = Q, .x0 = x0, .P0 = P0}};
        ekf.predict(u);

        for (size_t i = 0; i < 2; ++i) {
            CHECK(ekf.state()[i] == doctest::Approx(x_expected[i]).epsilon(1e-12));
            for (size_t j = 0; j < 2; ++j) {
                CHECK(ekf.covariance()(i, j) == doctest::Approx(P_expected(i, j)).epsilon(1e-12));
            }
        }
    }

    TEST_CASE("a millimetre-scale state in metres takes an analytic jacobian(x, u) or millimetre units, in float") {
        //! Central differences step the gap by ∛eps·max(|x|, 1) ≈ 4.9e-3 m in float, as large as
        //! the 5 mm gap itself (the documented limitation); the model's own Jacobian, or the gap
        //! stated in millimetres, gives the A₂₁ column. P(1, 0) ≈ Ts·A₂₁·P₀₀ + P₀₁ reads it.
        const float        gap = 0.005F;
        const float        current = gap * damp::sqrt(9.81F / 5e-4F);
        const double       Ts = 1e-4;
        const Matrix<2, 2> P0{{1e-6, 0.0}, {0.0, 1e-4}};
        const Matrix<2, 2> P_expected = rk4_covariance(Matrix<2, 2>{{0.0, 1.0}, {2.0 * 9.81 / static_cast<double>(gap), 0.0}}, Ts, P0);
        const auto         p10 = [&]<typename Model>(const Model& model) {
            const float                        unit = model.per_metre;
            ExtendedKalmanFilter<Model, float> ekf{model, {.Ts = static_cast<float>(Ts), .Q = Matrix<2, 2, float>{}, .x0 = ColVec<2, float>{gap * unit, 0.0F}, .P0 = P0.as<float>() * (unit * unit)}};
            ekf.predict(ColVec<1, float>{});
            return static_cast<double>(ekf.covariance()(1, 0) / (unit * unit)) / P_expected(1, 0);
        };
        CHECK(p10(MaglevAnalyticF{MaglevF{.i = current}}) == doctest::Approx(1.0).epsilon(1e-2));
        CHECK(p10(MaglevF{.i = current, .per_metre = 1000.0F}) == doctest::Approx(1.0).epsilon(1e-2));
        CHECK(damp::abs(p10(MaglevF{.i = current}) - 1.0) > 0.5); // central differences in metres
    }

    TEST_CASE("central differences keep the Jacobian for a large P₀ and a small P, in float") {
        //! The step is ∛eps·max(|x̂ᵢ|, 1), not tied to √Pᵢᵢ: a step of √P₀ = 1e3 rad would sample
        //! sin θ far from 0, and a step of a few 1e-8 would fall below the roundoff of f.
        SUBCASE("documented pendulum, P₀ = 1e6·I at θ = 0") {
            const Matrix<2, 2>                     P0 = Matrix<2, 2>::identity() * 1e6;
            const Matrix<2, 2>                     P_expected = rk4_covariance(Matrix<2, 2>{{0.0, 1.0}, {9.81, -0.5}}, 0.01, P0);
            ExtendedKalmanFilter<PendulumF, float> ekf{PendulumF{}, {.Ts = 0.01F, .Q = Matrix<2, 2, float>{}, .P0 = P0.as<float>()}};
            ekf.predict(ColVec<1, float>{});
            CHECK(static_cast<double>(ekf.covariance()(1, 0)) / P_expected(1, 0) == doctest::Approx(1.0).epsilon(1e-2));
        }
        SUBCASE("gyro bias with std(b) = 1e-5 rad/s at b = 0, ω = 3 rad/s") {
            const Matrix<2, 2>                        P0{{1e-4, 0.0}, {0.0, 1e-10}};
            const Matrix<2, 2>                        P_expected = rk4_covariance(Matrix<2, 2>{{0.0, -1.0}, {0.0, 0.0}}, 0.01, P0);
            ExtendedKalmanFilter<HeadingBiasF, float> ekf{HeadingBiasF{}, {.Ts = 0.01F, .Q = Matrix<2, 2, float>{}, .P0 = P0.as<float>()}};
            ekf.predict(ColVec<1, float>{3.0F});
            CHECK(static_cast<double>(ekf.covariance()(0, 1)) / P_expected(0, 1) == doctest::Approx(1.0).epsilon(1e-2)); // −Ts·P₁₁: the ∂θ̇/∂b column
        }
    }

    TEST_CASE("nonlinear range measurement, analytic and central-difference Jacobians") {
        // A static point observed through range-to-origin r = √(x² + y²):
        // a nonlinear h whose Jacobian is [x/r, y/r]. The EKF must localize the
        // point from range alone given a decent prior, with the analytic
        // Jacobian and with central differences alike.
        const ExtendedKalmanFilterConfig<Stationary, double> config{.Ts = 0.01, .Q = Matrix<2, 2>::identity() * 1e-6, .x0 = ColVec<2>{{2.5}, {4.5}}};
        const Matrix<1, 1>                                   R{1e-4};

        ExtendedKalmanFilter<Stationary, double> analytic{Stationary{}, config};
        ExtendedKalmanFilter<Stationary, double> numeric{Stationary{}, config};

        const ColVec<1> y{{5.0}}; // true point (3, 4)
        for (int k = 0; k < 50; ++k) {
            analytic.predict(kNoInput);
            numeric.predict(kNoInput);
            REQUIRE(analytic.update(RangeToOrigin{}, y, R));
            REQUIRE(numeric.update(RangeNumeric{}, y, R));
        }

        // Range-only obs can't fix the tangential direction, but the estimate
        // must at least land on the correct measurement manifold (range ≈ 5).
        const double r_est = damp::sqrt(analytic.state()[0] * analytic.state()[0] + analytic.state()[1] * analytic.state()[1]);
        CHECK(r_est == doctest::Approx(5.0).epsilon(1e-2));

        // A central difference carries a relative error of order ∛eps² ≈ 4e-11.
        const double fd_error = damp::cbrt(eps<double>()) * damp::cbrt(eps<double>());
        CHECK(numeric.state()[0] == doctest::Approx(analytic.state()[0]).epsilon(1e3 * fd_error));
        CHECK(numeric.state()[1] == doctest::Approx(analytic.state()[1]).epsilon(1e3 * fd_error));
    }

    TEST_CASE("two readings of different sizes on one filter") {
        // A scalar position reading and a 2-axis full-state reading update the same estimate,
        // each with its own (h, R). P = I, R₁ = 1, C₁ = [1 0]: the scalar reading y₁ = 2 moves
        // x₀ by P₀₀/(P₀₀ + R₁)·2 = 1 and leaves P = diag(0.5, 1). The 2-axis reading
        // y = (4, 2) with R = I then moves x₀ by 0.5/1.5·(4 − 1) = 1 and x₁ by 1/2·2 = 1, and
        // leaves P = diag(1/3, 1/2).
        ExtendedKalmanFilter<Stationary, double> ekf{Stationary{}, {.Ts = 0.1, .Q = Matrix<2, 2>::zeros()}};

        REQUIRE(ekf.update(Component<0>{}, ColVec<1>{2.0}, 1.0));
        CHECK(ekf.state()[0] == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(ekf.state()[1] == doctest::Approx(0.0));
        CHECK(ekf.covariance()(0, 0) == doctest::Approx(0.5).epsilon(1e-12));

        REQUIRE(ekf.update(BothStates{}, ColVec<2>{4.0, 2.0}, Matrix<2, 2>::identity()));
        CHECK(ekf.state()[0] == doctest::Approx(2.0).epsilon(1e-12));
        CHECK(ekf.state()[1] == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(ekf.covariance()(0, 0) == doctest::Approx(1.0 / 3.0).epsilon(1e-12));
        CHECK(ekf.covariance()(1, 1) == doctest::Approx(0.5).epsilon(1e-12));
    }

    TEST_CASE("predict-only grows covariance; update shrinks it") {
        ExtendedKalmanFilter<ConstantVelocity, double> ekf{ConstantVelocity{}, {.Ts = 0.1, .Q = Matrix<2, 2>::identity() * 1e-3, .P0 = Matrix<2, 2>::identity() * 0.1}};

        const double p0 = ekf.covariance()(0, 0);
        ekf.predict(kNoInput);
        const double p_after_predict = ekf.covariance()(0, 0);
        CHECK(p_after_predict > p0); // prediction adds uncertainty

        REQUIRE(ekf.update(Component<0>{}, ColVec<1>{{0.0}}, 1e-2));
        const double p_after_update = ekf.covariance()(0, 0);
        CHECK(p_after_update < p_after_predict); // measurement removes uncertainty
    }

    TEST_CASE("update returns false on singular innovation covariance") {
        // R = 0 and a measurement that sees nothing makes S = C P Cᵀ + R = 0,
        // not positive definite; the update must fail gracefully (false), not crash.
        ExtendedKalmanFilter<Stationary, double> ekf{Stationary{}, {.Ts = 0.1, .Q = Matrix<2, 2>::zeros(), .P0 = Matrix<2, 2>::zeros()}};
        ekf.predict(kNoInput);
        CHECK_FALSE(ekf.update(Blind{}, ColVec<1>{{1.0}}, 0.0));
        CHECK(ekf.state()[0] == 0.0); // the estimate is unchanged
    }

    TEST_CASE("sequential scalar updates with inter-measurement state clamping") {
        // The intended workflow for constrained estimation: fold in measurements one scalar
        // reading at a time, and clamp the affected state between updates via set_state.
        // Clamping a fused vector update can't enforce the constraint each measurement sees;
        // clamping between scalar updates can. Here state[0] is physically bounded to
        // [0, 1] (e.g. a duty cycle / SoC) and must never leave that box even
        // though a noisy measurement pushes the raw estimate past 1.
        ExtendedKalmanFilter<Stationary, double> ekf{Stationary{}, {.Ts = 0.01, .Q = Matrix<2, 2>::identity() * 1e-3, .x0 = ColVec<2>{{0.5}, {0.5}}}};
        const Matrix<1, 1>                       R{1e-2};

        for (int k = 0; k < 200; ++k) {
            ekf.predict(kNoInput);

            // Scalar update #1: measurement of state[0] sits above its physical
            // bound (1.2 > 1.0). Fuse it, then clamp before the next update.
            REQUIRE(ekf.update(Component<0>{}, ColVec<1>{{1.2}}, R));
            ColVec<2> x = ekf.state();
            x[0] = damp::clamp(x[0], 0.0, 1.0);
            ekf.set_state(x);

            // The clamped estimate is what the second scalar update sees.
            CHECK(ekf.state()[0] <= 1.0 + 1e-12);
            REQUIRE(ekf.update(Component<1>{}, ColVec<1>{{0.3}}, R));
        }

        CHECK(ekf.state()[0] == doctest::Approx(1.0)); // pinned at the bound
        CHECK(ekf.state()[1] == doctest::Approx(0.3).epsilon(0.05));
    }

    TEST_CASE("set_state / reset") {
        ExtendedKalmanFilter<Stationary, double> ekf{};
        ekf.set_state(ColVec<2>{{1.5}, {-2.0}});
        CHECK(ekf.state()[0] == doctest::Approx(1.5));
        CHECK(ekf.state()[1] == doctest::Approx(-2.0));

        ekf.reset(ColVec<2>{0.5, 7.0}, Matrix<2, 2>::identity() * 3.0);
        CHECK(ekf.state()[1] == doctest::Approx(7.0));
        CHECK(ekf.covariance()(0, 0) == doctest::Approx(3.0));
        CHECK(ekf.covariance()(1, 1) == doctest::Approx(3.0));
    }

    TEST_CASE("a non-positive Ts holds the state in predict") {
        const Matrix<2, 2> Q = Matrix<2, 2>::identity() * 1e-3;
        const ColVec<2>    x0{1.0, 2.0};
        for (const double Ts : {0.0, -0.1}) {
            ExtendedKalmanFilter<ConstantVelocity, double> ekf{ConstantVelocity{}, {.Ts = Ts, .Q = Q, .x0 = x0}};
            CHECK(ekf.sample_time() == 0.0);
            ekf.predict(kNoInput);
            CHECK(ekf.state()[0] == 1.0);
            CHECK(ekf.state()[1] == 2.0);
            CHECK(ekf.covariance()(0, 0) == doctest::Approx(1.0 + 1e-3));
        }
    }

    TEST_CASE("a lambda model in static constinit, with a lambda and a struct measurement") {
        // Constant velocity from x = 0, v = 0, P = I: a position reading y = 0.5 with R = 1 moves
        // p by P₀₀/(P₀₀ + R)·0.5 = 0.25.
        HoldsFilter    holder{};
        constexpr auto position = [](const ColVec<2, float>& x, const ColVec<1, float>& /*u*/) { return ColVec<1, float>{x[0]}; };
        auto           ekf = kLambdaEKF;
        REQUIRE(ekf.update(position, ColVec<1, float>{0.5F}, 1.0F));
        CHECK(ekf.state()[0] == doctest::Approx(0.25).epsilon(1e-5));
        REQUIRE(holder.ekf_.update(PositionF{}, ColVec<1, float>{0.5F}, 1.0F));
        CHECK(holder.ekf_.state()[0] == doctest::Approx(0.25).epsilon(1e-6));
        CHECK(kStructEKF.sample_time() == doctest::Approx(0.1));
    }

    TEST_CASE("the default scalar is the model's") {
        static_assert(std::is_same_v<ExtendedKalmanFilter<ConstantVelocityF>, ExtendedKalmanFilter<ConstantVelocityF, float>>);
        static_assert(std::is_same_v<ExtendedKalmanFilter<ConstantVelocity>, ExtendedKalmanFilter<ConstantVelocity, double>>);
        ExtendedKalmanFilter<ConstantVelocityF> ekf{ConstantVelocityF{}, {.Ts = 0.1F, .Q = Matrix<2, 2, float>::identity() * 1e-3F}};
        ekf.predict({});
        CHECK(ekf.update(PositionF{}, ColVec<1, float>{{0.5F}}, 1e-2F));
        CHECK(ekf.sample_time() == doctest::Approx(0.1));
    }
}
