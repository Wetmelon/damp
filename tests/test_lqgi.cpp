// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/fixed_gain_observer.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

// LQGI = LQI (state feedback with integral action) + Kalman estimator. The
// integral state drives steady-state output error to zero for constant
// references/disturbances. design::lqgi() was only covered indirectly;
// this exercises the design and the runtime tracking loop, including the
// integrator augmentation [x; xi].
//
// @see "Optimal Control" (Anderson & Moore, 1990); integral augmentation is the
//      standard servo/Type-1 construction.

namespace {
// Discrete double integrator at Ts = 0.1 s, measuring position.
constexpr StateSpace<2, 1, 1, double, 2, 1> make_plant() {
    return StateSpace<2, 1, 1, double, 2, 1>{
        Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
        Matrix<2, 1>{{0.005}, {0.1}},
        Matrix<1, 2>{{1.0, 0.0}},
        Matrix<1, 1>::zeros(),
        Matrix<2, 2>::identity(),
        Matrix<1, 1>::identity(),
        0.1 // Ts
    };
}

// The cost: position, velocity, the input, and the integral of the tracking error ∫e, weighted to pull the steady-state error to 0.
constexpr design::LQRConfig<2, 1, 1> kW{.Q = Matrix<2, 2>::identity(), .R = 0.1, .Q_integral = 1000.0};
} // namespace

TEST_SUITE("LQGI") {
    TEST_CASE("design succeeds with integral augmentation") {
        constexpr auto         sys = make_plant();
        constexpr Matrix<1, 1> R{{0.1}};
        constexpr Matrix<2, 2> Q_kf{{0.01, 0.0}, {0.0, 0.01}};
        constexpr Matrix<1, 1> R_kf{{0.1}};

        constexpr auto result = design::lqgi(sys, {.Q = kW.Q, .R = R, .Q_integral = kW.Q_integral}, design::KalmanNoise{.process = Q_kf, .measurement = R_kf});
        static_assert(result.success);
        static_assert(result.kalman.success);

        // Gain has a column per augmented state: [x(2) ; xi(1)] = 3.
        CHECK(result.K(0, 0) != 0.0);
        CHECK(result.K(0, 2) != 0.0); // integral gain is non-trivial
    }

    TEST_CASE("runtime LQGI tracks a constant reference with zero steady-state error") {
        const auto         sys = make_plant();
        const Matrix<1, 1> R{{0.1}};
        const Matrix<2, 2> Q_kf{{0.01, 0.0}, {0.0, 0.01}};
        const Matrix<1, 1> R_kf{{0.1}};

        const auto result = design::lqgi(sys, {.Q = kW.Q, .R = R, .Q_integral = kW.Q_integral}, design::KalmanNoise{.process = Q_kf, .measurement = R_kf});
        REQUIRE(result.success);

        auto controller = LQGI{result}; // exercises kf(result.kalman)

        const double r = 1.0;     // position setpoint
        ColVec<2>    x{0.0, 0.0}; // true plant state
        double       xi = 0.0;    // integral of (r - y)

        for (int k = 0; k < 600; ++k) {
            const ColVec<1> y{x[0]};
            controller.update(y);

            // Augmented state fed to LQI is [x_hat; xi].
            const auto      xhat = controller.estimator.state();
            const ColVec<3> x_aug{xhat[0], xhat[1], xi};
            const ColVec<1> u = -(controller.K * x_aug);

            x = sys.A * x + sys.B * u;
            controller.back_calculate(u);
            controller.predict();
            xi += sys.Ts * (r - x[0]); // integrate the tracking error, xi += Ts·e (matches the −Ts·C augmentation)
        }

        CHECK(x[0] == doctest::Approx(r).epsilon(0.02)); // output reaches the reference
    }

    TEST_CASE("runtime LQGI tracks via the controller's own integrator (update, control(r, x_ref, u_ref), predict)") {
        const auto sys = make_plant();
        const auto result = design::lqgi(sys, kW, {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = 0.1});
        REQUIRE(result.success);

        auto            controller = LQGI{result};
        const ColVec<1> r{1.0};
        ColVec<2>       x{0.0, 0.0};
        for (int k = 0; k < 600; ++k) {
            const ColVec<1> y{x[0]};
            controller.update(y);
            const ColVec<1> u = controller.control(r); // pulls x̂, the y from update, and owns xi
            x = sys.A * x + sys.B * u;
            controller.back_calculate(u);
            controller.predict();
        }
        CHECK(x[0] == doctest::Approx(r[0]).epsilon(0.02));

        controller.reset();
        CHECK(controller.integral.xi[0] == doctest::Approx(0.0));
    }

    TEST_CASE("LQGIResult::to_ss compensator tracks a reference") {
        // Validates the [r; y] -> u realization (steady-state L) by closing the loop.
        const auto sys = make_plant();
        const auto result = design::lqgi(sys, kW, {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = 0.1});
        REQUIRE(result.success);

        const auto ss = result.to_ss();            // [x̂⁻; xi; u_prev], in [r; y], out u
        CHECK(ss.D(0, 0) == doctest::Approx(0.0)); // no direct term from r
        CHECK(ss.Ts == doctest::Approx(sys.Ts));

        auto            runtime = LQGI{result};
        ColVec<4>       xc{}; // NX + NY + NU
        ColVec<2>       xp{0.0, 0.0};
        ColVec<2>       xp_rt = xp;
        const ColVec<1> r{1.0};
        for (int k = 0; k < 8; ++k) {
            const ColVec<2> in{r[0], xp[0]};
            const ColVec<1> u = ss.C * xc + ss.D * in;
            const ColVec<1> u_rt = runtime.step(r, ColVec<1>{{xp_rt[0]}});
            CHECK(u[0] == doctest::Approx(u_rt[0]));
            xp = sys.A * xp + sys.B * u;
            xp_rt = sys.A * xp_rt + sys.B * u_rt;
            xc = ss.A * xc + ss.B * in;
        }

        for (int k = 0; k < 600; ++k) {
            const ColVec<2> in{r[0], xp[0]};
            const ColVec<1> u = ss.C * xc + ss.D * in;
            xp = sys.A * xp + sys.B * u;
            xc = ss.A * xc + ss.B * in;
        }
        CHECK(xp[0] == doctest::Approx(r[0]).epsilon(0.02));
    }

    TEST_CASE("LQGIResult::as<float>() preserves the design") {
        const auto sys = make_plant();
        const auto result = design::lqgi(sys, kW, {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = 0.1});
        REQUIRE(result.success);

        const auto rf = result.as<float>();
        CHECK(rf.success);
        CHECK(rf.K(0, 0) == doctest::Approx(static_cast<float>(result.K(0, 0))));
        CHECK(rf.kalman.M(0, 0) == doctest::Approx(static_cast<float>(result.kalman.M(0, 0))));
    }

    TEST_CASE("integrating a subset of the sensors: same argument shape as lqi") {
        // Position and velocity sensors; integrate the position only.
        constexpr StateSpace<2, 1, 2, double, 2, 2> sys2{
            Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}}, Matrix<2, 1>{{0.005}, {0.1}},
            Matrix<2, 2>::identity(), Matrix<2, 1>{}, Matrix<2, 2>::identity(), Matrix<2, 2>::identity(), 0.1
        };
        const Matrix<2, 2> Q_kf{{0.01, 0.0}, {0.0, 0.01}};
        const Matrix<2, 2> R_kf{{0.1, 0.0}, {0.0, 0.1}};
        const auto         res = design::lqgi(sys2, kW, design::KalmanNoise{.process = Q_kf, .measurement = R_kf}, design::IntegralAction<1>{{0}});
        REQUIRE(res.success);
        const auto lqi_only = design::lqi(sys2, kW, design::IntegralAction<1>{{0}});
        CHECK(res.K(0, 2) == doctest::Approx(lqi_only.K(0, 2)));

        // The compensator realization takes [r (NI); y (NY)] and matches the runtime.
        const auto      ss = res.to_ss();
        LQGI            runtime{res.as<float>()};
        ColVec<4>       xc{}; // NX + NI + NU
        ColVec<2>       xp{};
        const ColVec<1> r{1.0};
        for (int k = 0; k < 600; ++k) {
            const ColVec<2> y = sys2.C * xp;
            const ColVec<3> in{r[0], y[0], y[1]};
            const ColVec<1> u = ss.C * xc + ss.D * in;
            const ColVec<1> u_rt = runtime.step(r.as<float>(), y.as<float>()).as<double>();
            if (k < 8) {
                CHECK(u[0] == doctest::Approx(u_rt[0]).epsilon(1e-4));
            }
            xp = sys2.A * xp + sys2.B * u;
            xc = ss.A * xc + ss.B * in;
        }
        CHECK(xp[0] == doctest::Approx(1.0).epsilon(0.02));

        const auto bad = design::lqgi(sys2, kW, design::KalmanNoise{.process = Q_kf, .measurement = R_kf}, design::IntegralAction<1>{{5}});
        CHECK_FALSE(bad.success);
        CHECK(bad.status == design::DesignStatus::bad_weights);
    }

    TEST_CASE("predict holds the integrator against the clamp") {
        const auto result = design::lqgi(make_plant(), kW, {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = 0.1});
        REQUIRE(result.success);
        const auto   sys = make_plant();
        const double u_max = 0.2;
        const auto   run = [&](bool aware) {
            auto      ctrl = LQGI{result};
            ColVec<2> x{};
            double    peak = 0.0;
            for (int k = 0; k < 2000; ++k) {
                const ColVec<1> y{{x[0]}};
                ctrl.update(y);
                const ColVec<1> u_req = ctrl.control(ColVec<1>{{1.0}});
                const ColVec<1> u{{damp::clamp(u_req[0], -u_max, u_max)}};
                // The unaware loop hides the clamp from the integrator (it still feeds the estimator u).
                if (aware) {
                    ctrl.back_calculate(u);
                    ctrl.predict();
                } else {
                    ctrl.back_calculate(u_req); // the integrator told the request, the filter the clamped input
                    ctrl.predict();
                    ctrl.estimator.set_state(ctrl.estimator.state() + (ctrl.estimator.model().B * (u - u_req)));
                }
                x = sys.A * x + sys.B * u;
                peak = damp::max(peak, x[0]);
            }
            return damp::array<double, 2>{peak, x[0]};
        };
        const auto naive = run(false);
        const auto aware = run(true);
        CHECK(aware[0] < naive[0]);
        CHECK(aware[1] == doctest::Approx(1.0).epsilon(1e-3));
    }

    TEST_CASE("LQGI deduces sizes and scalar from its design") {
        static constexpr auto res = design::lqgi(make_plant(), kW, {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = 0.1});
        static constinit LQGI law{res.as<float>()}; // sizes and scalar deduced from the result
        static_assert(std::is_same_v<decltype(law), LQGI<2, 1, 1, 1, float, FixedGainObserver<2, 1, 1, float>>>);
        // A double design does not narrow to the float default.
        static_assert(std::is_same_v<decltype(LQGI{res}), LQGI<2, 1, 1, 1, double, FixedGainObserver<2, 1, 1, double>>>);
        CHECK(law.K(0, 0) == doctest::Approx(static_cast<float>(res.K(0, 0))));
    }
}

TEST_SUITE("LQGI status") {
    TEST_CASE("a continuous plant is refused by both halves") {
        auto continuous = make_plant();
        continuous.Ts = 0.0;
        const auto refused = design::lqgi(continuous, kW, {});
        CHECK_FALSE(refused.success);
        CHECK(refused.status == design::DesignStatus::invalid_sample_time);
        CHECK(refused.kalman.status == design::DesignStatus::invalid_sample_time);

        constexpr auto ok = design::lqgi(make_plant(), kW, {});
        static_assert(ok.success && ok.status == design::DesignStatus::ok && ok.Ts == 0.1);
        CHECK(ok.as<float>().status == design::DesignStatus::ok);
    }
}

namespace {
// Continuous double integrator, position measured, process noise on the velocity.
constexpr StateSpace<2, 1, 1, double, 1, 1> make_continuous_plant() {
    return StateSpace<2, 1, 1, double, 1, 1>{
        .A = Matrix<2, 2>{{0.0, 1.0}, {0.0, 0.0}},
        .B = Matrix<2, 1>{{0.0}, {1.0}},
        .C = Matrix<1, 2>{{1.0, 0.0}},
        .G = Matrix<2, 1>{{0.0}, {1.0}},
        .H = Matrix<1, 1>{{1.0}},
    };
}
constexpr design::LQRConfig<2, 1, 1> kContinuousWeights{.Q = Matrix<2, 2>::diagonal({10.0, 1.0}), .R = 1.0, .Q_integral = 100.0};
constexpr design::KalmanNoise<1, 1>  kNoise{.process = 0.5, .measurement = 1e-4}; // density; per-reading variance
constexpr double                     kTs = 1e-3;
} // namespace

TEST_SUITE("LQGI sampled-data (continuous plant, Ts)") {
    TEST_CASE("each half is its own sampled-data design; the measurement variance goes in unscaled") {
        static constexpr auto res = design::lqgi(make_continuous_plant(), kTs, kContinuousWeights, kNoise);
        static_assert(res.success && res.Ts == kTs && res.kalman.sys.Ts == kTs);
        static_assert(res.integral.step == kTs);
        static_assert(res.K == design::lqi(make_continuous_plant(), kTs, kContinuousWeights).K);
        static_assert(res.kalman.M == design::kalman(make_continuous_plant(), kTs, kNoise).M);
        static_assert(res.kalman.R == kNoise.measurement);

        //! Braced aggregates deduce at the call, and the index form matches.
        const auto braced = design::lqgi(make_continuous_plant(), kTs, {.Q = kContinuousWeights.Q, .R = kContinuousWeights.R, .Q_integral = kContinuousWeights.Q_integral}, design::KalmanNoise{.process = kNoise.process, .measurement = kNoise.measurement}, design::IntegralAction<1>{{0}});
        REQUIRE(braced.success);
        for (size_t j = 0; j < 3; ++j) {
            CHECK(braced.K(0, j) == doctest::Approx(res.K(0, j)).epsilon(1e-9)); // run time vs constant-evaluated
        }
        for (size_t i = 0; i < 2; ++i) {
            CHECK(braced.kalman.M(i, 0) == doctest::Approx(res.kalman.M(i, 0)).epsilon(1e-9));
        }
    }

    TEST_CASE("it deploys constinit in float and tracks r with zero steady-state error in SIL") {
        static constexpr auto res = design::lqgi(make_continuous_plant(), kTs, kContinuousWeights, kNoise);
        static constinit auto law = LQGI{res.as<float>()};
        CHECK(law.integral.step == static_cast<float>(kTs));

        const auto      plant_d = res.kalman.sys; // the ZOH plant: the continuous plant sampled exactly
        const ColVec<1> r{1.0};
        ColVec<2>       x{};
        for (int k = 0; k < 20000; ++k) {
            const ColVec<1> y = plant_d.C * x;
            const ColVec<1> u = law.step(r.as<float>(), y.as<float>()).as<double>();
            x = plant_d.A * x + plant_d.B * (u + ColVec<1>{0.5}); // constant input disturbance
        }
        CHECK(x[0] == doctest::Approx(1.0).epsilon(1e-4));
        CHECK(damp::abs(x[1]) < 1e-4);
    }

    TEST_CASE("to_ss realizes the Ts-form runtime, integrator step included") {
        const auto result = design::lqgi(make_continuous_plant(), kTs, kContinuousWeights, kNoise);
        REQUIRE(result.success);
        const auto      ss = result.to_ss();
        auto            runtime = LQGI{result};
        const auto&     sys = result.kalman.sys;
        ColVec<4>       xc{};
        ColVec<2>       xp{};
        ColVec<2>       xp_rt{};
        const ColVec<1> r{1.0};
        for (int k = 0; k < 200; ++k) {
            const ColVec<2> in{r[0], xp[0]};
            const ColVec<1> u = ss.C * xc + ss.D * in;
            const ColVec<1> u_rt = runtime.step(r, ColVec<1>{{xp_rt[0]}});
            CHECK(u[0] == doctest::Approx(u_rt[0]).epsilon(1e-9));
            // xi = ∫e dt in both, through the previous tick: the step holds this tick's error for its next predict.
            CHECK(xc[2] == doctest::Approx(runtime.integral.xi[0]).epsilon(1e-9));
            xp = sys.A * xp + sys.B * u;
            xp_rt = sys.A * xp_rt + sys.B * u_rt;
            xc = ss.A * xc + ss.B * in;
        }
    }

    TEST_CASE("Ts form refuses a discrete plant and a bad Ts") {
        auto discrete = make_continuous_plant();
        discrete.Ts = kTs;
        const auto twice = design::lqgi(discrete, kTs, kContinuousWeights, kNoise);
        CHECK(twice.status == design::DesignStatus::invalid_sample_time);
        CHECK(twice.kalman.status == design::DesignStatus::invalid_sample_time);
        CHECK(design::lqgi(make_continuous_plant(), 0.0, kContinuousWeights, kNoise).status == design::DesignStatus::invalid_sample_time);
        CHECK(design::lqgi(make_continuous_plant(), -kTs, kContinuousWeights, kNoise).status == design::DesignStatus::invalid_sample_time);
        CHECK(design::lqgi(make_continuous_plant(), kTs, kContinuousWeights, kNoise, design::IntegralAction<1>{{2}}).status == design::DesignStatus::bad_weights);
    }
}

TEST_SUITE("LQGI on a plant with feedthrough") {
    TEST_CASE("the design and the runtime integrate the same output: y holds r, and to_ss realizes the runtime") {
        // ẋ = −x + u, y = x + D·u; the reading at k carries the input held over the period before it.
        for (const double D : {0.5, -2.0}) {
            const double                            a = damp::exp(-0.1);
            const double                            b = 1.0 - a;
            const StateSpace<1, 1, 1, double, 1, 1> pd{.A = Matrix<1, 1>{{a}}, .B = Matrix<1, 1>{{b}}, .C = Matrix<1, 1>{{1.0}}, .D = Matrix<1, 1>{{D}}, .G = Matrix<1, 1>{{1.0}}, .H = Matrix<1, 1>{{1.0}}, .Ts = 0.1};
            StateSpace<1, 1, 1, double, 1, 1>       pc = pd;
            pc.A = Matrix<1, 1>{{-1.0}};
            pc.B = Matrix<1, 1>{{1.0}};
            pc.Ts = 0.0;
            const auto w = design::LQRConfig<1, 1, 1>{.Q = 1.0, .R = 1.0, .Q_integral = 10.0};
            const auto n = design::KalmanNoise<1, 1, double>{.process = 1.0, .measurement = 0.01};
            for (const bool sampled : {false, true}) {
                const auto res = sampled ? design::lqgi(pc, 0.1, w, n) : design::lqgi(pd, w, n);
                REQUIRE(res.success);
                auto      law = LQGI{res};
                auto      runtime = LQGI{res};
                auto      ss = res.to_ss();
                double    x = 0.0;
                double    u_prev = 0.0;
                double    y = 0.0;
                double    xs = 0.0;
                double    us_prev = 0.0;
                ColVec<3> s{};
                for (int k = 0; k < 2000; ++k) {
                    const double u = law.step(ColVec<1>{1.0}, ColVec<1>{x + (D * u_prev)})[0];
                    y = x + (D * u);
                    x = (a * x) + (b * u);
                    u_prev = u;
                    if (k < 200) {
                        const ColVec<2> in{1.0, xs + (D * us_prev)};
                        const double    u_ss = (ss.C * s + ss.D * in)[0];
                        const double    u_rt = runtime.step(ColVec<1>{1.0}, ColVec<1>{in[1]})[0];
                        CHECK(u_ss == doctest::Approx(u_rt).epsilon(1e-9));
                        s = ss.A * s + ss.B * in;
                        xs = (a * xs) + (b * u_rt);
                        us_prev = u_rt;
                    }
                }
                CHECK(y == doctest::Approx(1.0).epsilon(1e-9));
            }
        }
    }
}
