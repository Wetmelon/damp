// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "damp/controllers/oflqg.hpp"
#include "damp/controllers/state_feedback_types.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/status.hpp"
#include "damp/design/steady_state_target.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {
constexpr StateSpace<2, 1, 1, double, 2, 1> make_plant() {
    return StateSpace<2, 1, 1, double, 2, 1>{
        Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
        Matrix<2, 1>{{0.005}, {0.1}},
        Matrix<1, 2>{{1.0, 0.0}},
        Matrix<1, 1>::zeros(),
        Matrix<2, 2>::identity(),
        Matrix<1, 1>::identity(),
        0.1
    };
}
} // namespace

TEST_SUITE("Offset-free LQG: oflqg with a disturbance model") {
    TEST_CASE("design succeeds") {
        constexpr auto sys = make_plant();
        constexpr auto result = design::oflqg(sys, {.Q = Matrix<2, 2>::identity(), .R = 0.1}, design::KalmanNoise{Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, Matrix<1, 1>{{0.1}}, Matrix<1, 1>{{1.0}}}, design::InputDisturbance{});
        static_assert(result.success);
        static_assert(result.status == design::DesignStatus::ok);
        static_assert(result.Ts == 0.1 && result.kalman.sys.Ts == 0.1); // the plant's, both halves
        CHECK(result.success);
        CHECK(result.K(0, 0) != 0.0);
        CHECK(result.kalman.M(0, 0) != 0.0);
        CHECK(result.as<float>().status == design::DesignStatus::ok);

        //! A velocity held at r: [A−I, B; C, D] is singular (the first column is zero), so there is no unique target.
        auto velocity = sys;
        velocity.C = Matrix<1, 2>{{0.0, 1.0}};
        CHECK(design::steady_state_target(velocity).status == design::DesignStatus::singular);
    }

    TEST_CASE("a continuous plant is refused, as by lqg; continuous plant + Ts designs at Ts") {
        constexpr design::LQRConfig<2, 1> w{.Q = Matrix<2, 2>::identity(), .R = 0.1};
        auto                              continuous = make_plant();
        continuous.Ts = 0.0;
        const auto refused = design::oflqg(continuous, w, design::KalmanNoise<2, 1, double, 1>{}, design::InputDisturbance{});
        CHECK_FALSE(refused.success);
        CHECK(refused.status == design::DesignStatus::invalid_sample_time);
        const auto refused_model = design::oflqg(continuous, w, design::KalmanNoise<2, 1, double, 1>{}, design::DisturbanceModel<2, 1, 1>{.Bd = continuous.B});
        CHECK(refused_model.status == design::DesignStatus::invalid_sample_time);

        //! The sample-time form takes the continuous plant only: a discrete one would be discretized twice.
        constexpr StateSpace<2, 1, 1> di{.A = Matrix<2, 2>{{0.0, 1.0}, {0.0, 0.0}}, .B = Matrix<2, 1>{{0.0}, {1.0}}, .C = Matrix<1, 2>{{1.0, 0.0}}};
        constexpr auto                sampled = design::oflqg(di, 0.1, w, design::KalmanNoise<0, 0, double, 1>{}, design::InputDisturbance{});
        static_assert(sampled.success && sampled.Ts == 0.1);
        CHECK(design::oflqg(*discretize(di, 0.1), 0.1, w, design::KalmanNoise<0, 0, double, 1>{}, design::InputDisturbance{}).status == design::DesignStatus::invalid_sample_time);
        CHECK(design::oflqg(di, 0.0, w, design::KalmanNoise<0, 0, double, 1>{}, design::InputDisturbance{}).status == design::DesignStatus::invalid_sample_time);
    }

    TEST_CASE("tracks a setpoint under a constant input load") {
        const auto sys = make_plant();
        const auto result = design::oflqg(sys, {.Q = Matrix<2, 2>::identity(), .R = 0.1}, design::KalmanNoise{Matrix<2, 2>{{0.05, 0.0}, {0.0, 0.05}}, Matrix<1, 1>{{0.05}}, Matrix<1, 1>{{1.0}}}, design::InputDisturbance{});
        REQUIRE(result.success);

        auto            ctrl = OFLQG{result};
        ColVec<2>       x{0.0, 0.0};
        const ColVec<1> r{{1.0}};
        const double    w_true = 0.4;
        for (int k = 0; k < 400; ++k) {
            const ColVec<1> y{{x[0]}};
            const ColVec<1> u = ctrl.step(r, y);
            x = sys.A * x + sys.B * (u + ColVec<1>{{w_true}});
        }
        CHECK(x[0] == doctest::Approx(1.0).epsilon(0.05));
        CHECK(ctrl.disturbance_estimate()(0) == doctest::Approx(w_true).epsilon(0.15));
    }

    TEST_CASE("as<float> preserves success") {
        constexpr auto result = design::oflqg(make_plant(), {.Q = Matrix<2, 2>::identity(), .R = 0.1}, design::KalmanNoise{Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, Matrix<1, 1>{{0.1}}, Matrix<1, 1>{{1.0}}}, design::InputDisturbance{});
        constexpr auto f = result.as<float>();
        static_assert(f.success);
        CHECK(f.success);
    }

    TEST_CASE("more sensors than inputs: a load disturbance, one controlled output") {
        // Mass on a spring behind a compliant drive: x = [p_motor, v_motor, p_load, v_load].
        // Measured: motor position, load position, load acceleration (which feels the load
        // force directly). The disturbance is a force on the load; only the load is held at r.
        constexpr double Ts = 0.001;
        constexpr double k = 400.0;
        constexpr double m = 1.0;
        Matrix<4, 4>     Ac{};
        Ac(0, 1) = 1.0;
        Ac(1, 0) = -k / m;
        Ac(1, 2) = k / m;
        Ac(1, 1) = -2.0;
        Ac(2, 3) = 1.0;
        Ac(3, 0) = k / m;
        Ac(3, 2) = -k / m;
        Ac(3, 3) = -2.0;
        StateSpace<4, 2, 3> cont{.A = Ac};
        cont.B(1, 0) = 1.0 / m; // motor force
        cont.B(3, 1) = 1.0 / m; // load force (the disturbance channel)
        cont.C(0, 0) = 1.0;
        cont.C(1, 2) = 1.0;
        for (std::size_t c = 0; c < 4; ++c) {
            cont.C(2, c) = Ac(3, c);
        }
        const auto disc = discretize(cont, Ts, DiscretizationMethod::ZOH);
        REQUIRE(disc.has_value());
        StateSpace<4, 1, 3, double, 4, 3> sys{.A = disc->A, .C = disc->C}; // the process noise on every state (G = I)
        sys.B = disc->B.template block<4, 1>(0, 0);
        sys.Ts = Ts;

        design::DisturbanceModel<4, 3, 1, 1> model{.outputs = {{1}}}; // hold the load position, output 1
        model.Bd = disc->B.template block<4, 1>(0, 1);
        model.Cd(2, 0) = 1.0 / m;

        Matrix<4, 4> Q = Matrix<4, 4>::identity();
        Q(2, 2) = 1.0e4;
        const auto res = design::oflqg(sys, {.Q = Q, .R = 1.0e-4}, design::KalmanNoise{Matrix<4, 4>::identity() * 1.0e-8, Matrix<3, 3>{{1.0e-8, 0.0, 0.0}, {0.0, 1.0e-8, 0.0}, {0.0, 0.0, 1.0e-2}}, Matrix<1, 1>{{1.0}}}, model);
        REQUIRE(res.success);

        OFLQG           ctrl{res};
        ColVec<4>       x{};
        const ColVec<1> r{{0.01}};
        const double    load_force = -0.8;
        for (int i = 0; i < 4000; ++i) {
            const ColVec<3> y = (sys.C * x) + (model.Cd * ColVec<1>{{load_force}});
            const ColVec<1> u = ctrl.step(r, y);
            x = (sys.A * x) + (sys.B * u) + (model.Bd * ColVec<1>{{load_force}});
        }
        CHECK(x[2] == doctest::Approx(0.01).epsilon(1e-3));                                 // load held at r
        CHECK(ctrl.disturbance_estimate()(0) == doctest::Approx(load_force).epsilon(0.02)); // and the force found
    }

    TEST_CASE("tracking a constant reference is the setpoint form") {
        // Same design, same plant and load: step({}, x*, u*, y) with x* = [1, 0], u* = 0 must
        // produce exactly the inputs step(r = 1, {}, {}, y) does (re-anchoring makes the two filters equal).
        const auto sys = make_plant();
        const auto result = design::oflqg(sys, {.Q = Matrix<2, 2>::identity(), .R = 0.1}, design::KalmanNoise{Matrix<2, 2>{{0.05, 0.0}, {0.0, 0.05}}, Matrix<1, 1>{{0.05}}, Matrix<1, 1>{{1.0}}}, design::InputDisturbance{});
        REQUIRE(result.success);
        OFLQG     setpoint{result};
        OFLQG     tracking{result};
        ColVec<2> xs{};
        ColVec<2> xt{};
        double    worst = 0.0;
        for (int k = 0; k < 300; ++k) {
            const ColVec<1> us = setpoint.step(ColVec<1>{{1.0}}, ColVec<1>{{xs[0]}});
            const ColVec<1> ut = tracking.step({}, ColVec<2>{1.0, 0.0}, ColVec<1>{}, ColVec<1>{{xt[0]}});
            worst = std::max(worst, std::abs(us[0] - ut[0]));
            xs = sys.A * xs + sys.B * (us + ColVec<1>{{0.4}});
            xt = sys.A * xt + sys.B * (ut + ColVec<1>{{0.4}});
        }
        CHECK(worst < 1e-9);
    }

    TEST_CASE("tracks a moving reference under a constant load") {
        const auto sys = make_plant();
        const auto result = design::oflqg(sys, {.Q = Matrix<2, 2>{{100.0, 0.0}, {0.0, 1.0}}, .R = 0.01}, design::KalmanNoise{Matrix<2, 2>{{0.05, 0.0}, {0.0, 0.05}}, Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{1.0}}}, design::InputDisturbance{});
        REQUIRE(result.success);
        OFLQG        ctrl{result};
        ColVec<2>    x{};
        const double w = 0.3; // [rad/s]
        const double Ts = 0.1;
        double       worst = 0.0;
        for (int k = 0; k < 2000; ++k) {
            const double    t = k * Ts;
            const ColVec<2> x_ref{std::sin(w * t), w * std::cos(w * t)};
            const ColVec<1> u_ref{{-w * w * std::sin(w * t)}}; // the reference's own acceleration
            const ColVec<1> u = ctrl.step({}, x_ref, u_ref, ColVec<1>{{x[0]}});
            x = sys.A * x + sys.B * (u + ColVec<1>{{0.4}});
            if (t > 100.0) {
                worst = std::max(worst, std::abs(x[0] - std::sin(w * (t + Ts))));
            }
        }
        CHECK(worst < 0.01);                                                         // on the sine, not lagging it
        CHECK(ctrl.disturbance_estimate()(0) == doctest::Approx(0.4).epsilon(0.05)); // and the load found
    }

    TEST_CASE("the offset-free LQG is offset-free with feedthrough D != 0") {
        const StateSpace<2, 1, 1, double, 2, 1> sys{
            Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
            Matrix<2, 1>{{0.005}, {0.1}},
            Matrix<1, 2>{{1.0, 0.0}},
            Matrix<1, 1>{{0.5}},
            Matrix<2, 2>::identity(),
            Matrix<1, 1>::identity(),
            0.1
        };
        const auto result = design::oflqg(
            sys, {.Q = Matrix<2, 2>::identity(), .R = 0.1},
            design::KalmanNoise{Matrix<2, 2>{{0.05, 0.0}, {0.0, 0.05}}, Matrix<1, 1>{{0.05}}, Matrix<1, 1>{{1.0}}},
            design::DisturbanceModel<2, 1, 1>{.Bd = sys.B} // the load at the input; the sensor reads the command alone
        );
        REQUIRE(result.success);

        OFLQG           ctrl{result};
        ColVec<2>       x{0.0, 0.0};
        ColVec<1>       u_held{};
        const ColVec<1> r{{1.0}};
        double          y = 0.0;
        for (int k = 0; k < 600; ++k) {
            // Sampled before this tick's update: the held input is still u[k-1].
            y = x[0] + (0.5 * u_held(0));
            u_held = ctrl.step(r, ColVec<1>{{y}});
            x = sys.A * x + sys.B * (u_held + ColVec<1>{{0.4}});
        }
        CHECK(y == doctest::Approx(1.0).epsilon(1e-3));
    }

    TEST_CASE("a disturbance that reaches the controlled output (C_d ≠ 0) leaves no offset") {
        // An output bias (C_d = I), and an input load on a plant with feedthrough (C_d = D). The
        // target holds C_z·x + D_z·u + C_{d,z}·d̂ at r; without the C_d term y settled at r + C_d·d.
        const auto run = [](double D, const auto& spec, bool output_bias) {
            const StateSpace<2, 1, 1, double, 2, 1> sys{
                Matrix<2, 2>{{0.9, 0.1}, {0.0, 0.8}}, Matrix<2, 1>{{0.0}, {0.2}}, Matrix<1, 2>{{1.0, 0.0}}, Matrix<1, 1>{{D}}, Matrix<2, 2>::identity(), Matrix<1, 1>::identity(), 0.01
            };
            const design::KalmanNoise<2, 1, double, 1> noise{.process = Matrix<2, 2>::identity() * 0.01, .measurement = 0.01, .disturbance = 0.01};
            const design::LQRConfig<2, 1>              w{.Q = Matrix<2, 2>::identity(), .R = 1.0};
            const auto                                 result = design::oflqg(sys, w, noise, spec);
            REQUIRE(result.success);
            OFLQG        law{result};
            ColVec<2>    x{};
            double       u = 0.0;
            double       y = 0.0;
            const double d = 0.3;
            for (int k = 0; k < 6000; ++k) {
                const double d_in = output_bias ? 0.0 : d;
                y = (sys.C * x)(0) + (D * (u + d_in)) + (output_bias ? d : 0.0); // the held input is u[k−1]
                u = law.step(ColVec<1>{{1.0}}, ColVec<1>{{y}})(0);
                x = (sys.A * x) + (sys.B * ColVec<1>{{u + d_in}});
            }
            CHECK(law.disturbance_estimate()(0) == doctest::Approx(d).epsilon(1e-4));
            return y;
        };
        CHECK(run(0.0, design::OutputDisturbance{}, true) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(run(0.5, design::OutputDisturbance{}, true) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(run(0.5, design::InputDisturbance{}, false) == doctest::Approx(1.0).epsilon(1e-6));
    }
}
