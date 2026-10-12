// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "damp/backend.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {
// Discrete double integrator at Ts = 0.1 s, position measured. NX=2, NU=1, NY=1.
constexpr StateSpace<2, 1, 1, double, 1, 1> make_plant() {
    return StateSpace<2, 1, 1, double, 1, 1>{
        Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}}, // A
        Matrix<2, 1>{{0.005}, {0.1}},         // B
        Matrix<1, 2>{{1.0, 0.0}},             // C: measure position
        Matrix<1, 1>::zeros()                 // D
    };
}

constexpr auto design_lqi() {
    return design::discrete_lqi(
        make_plant(),
        Matrix<3, 3>{{10.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 5.0}}, // Q on [x; xi]
        Matrix<1, 1>{{1.0}}                                               // R
    );
}
} // namespace

TEST_SUITE("LQI") {
    TEST_CASE("LQIResult::as<U>() conversion") {
        constexpr auto lqi_d = design::discrete_lqi(
            StateSpace<1, 1, 1, double, 1, 1>{Matrix<1, 1>{{1.0}}, Matrix<1, 1>{{1.0}}, Matrix<1, 1>{{1.0}}, Matrix<1, 1>::zeros()},
            Matrix<2, 2>{{1.0, 0.0}, {0.0, 1.0}}, Matrix<1, 1>{{1.0}}
        );

        constexpr auto lqi_f = lqi_d.as<float>();

        static_assert(lqi_f.success);
        CHECK(lqi_f.success);
        CHECK(lqi_f.K(0, 0) != 0.0f);
        CHECK(lqi_f.K(0, 1) != 0.0f);
    }

    TEST_CASE("runtime control(x_aug) applies u = -K*[x; xi]") {
        constexpr auto r = design_lqi();
        REQUIRE(r.success);

        LQI<2, 1, 1, double> ctrl{r};
        const ColVec<3>      x_aug{{0.5}, {-0.2}, {1.3}};
        const ColVec<1>      u = ctrl.control(x_aug);
        const ColVec<1>      expected = -r.K * x_aug;
        CHECK(u[0] == doctest::Approx(expected[0]).epsilon(1e-12));
    }

    TEST_CASE("own integrator matches manual augmented state") {
        constexpr auto       r = design_lqi();
        LQI<2, 1, 1, double> a{r}; // uses internal integrator
        LQI<2, 1, 1, double> b{r}; // fed a manual augmented vector

        ColVec<1>       xi{{0.0}};
        const ColVec<2> states[] = {{{0.5}, {0.0}}, {{0.4}, {-0.1}}, {{0.2}, {0.05}}};
        const double    refs[] = {1.0, 1.0, 1.0};
        for (int k = 0; k < 3; ++k) {
            const ColVec<1> y{{states[k][0]}};
            const ColVec<1> rf{{refs[k]}};

            const ColVec<3> x_aug{{states[k][0]}, {states[k][1]}, {xi[0]}};
            const ColVec<1> ub = b.control(x_aug);            // manual
            const ColVec<1> ua = a.control(rf, y, states[k]); // owns xi
            CHECK(ua[0] == doctest::Approx(ub[0]).epsilon(1e-12));

            xi[0] += refs[k] - y[0]; // mirror the controller's update
            CHECK(a.xi[0] == doctest::Approx(xi[0]).epsilon(1e-12));
        }
    }

    TEST_CASE("closed loop drives output to the reference with zero steady-state error") {
        const auto sys = make_plant();
        const auto r = design_lqi();
        REQUIRE(r.success);

        LQI<2, 1, 1, double> ctrl{r};
        ColVec<2>            x{{0.0}, {0.0}};
        const ColVec<1>      ref{{1.0}};
        for (int k = 0; k < 500; ++k) {
            const ColVec<1> y{{x[0]}};
            const ColVec<1> u = ctrl.control(ref, y, x);
            x = sys.A * x + sys.B * u;
        }
        CHECK(x[0] == doctest::Approx(1.0).epsilon(1e-3)); // tracks position
        CHECK(x[1] == doctest::Approx(0.0).epsilon(1e-3)); // at rest

        ctrl.reset();
        CHECK(ctrl.xi[0] == doctest::Approx(0.0));
    }

    TEST_CASE("integrating a subset of outputs matches a hand-built single-output model") {
        // Same double integrator, two sensors (position, velocity); integrate position only.
        constexpr StateSpace<2, 1, 2, double> sys2{
            Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}}, Matrix<2, 1>{{0.005}, {0.1}},
            Matrix<2, 2>::identity(), Matrix<2, 1>{}
        };
        constexpr Matrix<3, 3> Q{{10.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 5.0}};
        constexpr auto         picked = design::discrete_lqi(sys2, Q, Matrix<1, 1>{{1.0}}, {0});
        constexpr auto         by_hand = design_lqi(); // C = [1 0]
        static_assert(picked.success);
        for (size_t j = 0; j < 3; ++j) {
            CHECK(picked.K(0, j) == doctest::Approx(by_hand.K(0, j)).epsilon(1e-12));
        }
        CHECK(picked.Ci(0, 0) == 1.0);
        CHECK(picked.Ci(0, 1) == 0.0);

        // The selector form takes the reference for the integrated output and the full measurement.
        LQI             ctrl{picked.as<float>()};
        ColVec<2>       x{};
        const ColVec<1> r{{1.0}};
        for (int k = 0; k < 500; ++k) {
            const ColVec<2, float> y{static_cast<float>(x[0]), static_cast<float>(x[1])};
            const ColVec<1, float> u = ctrl.control(r.as<float>(), y, y);
            x = sys2.A * x + sys2.B * u.as<double>();
        }
        CHECK(x[0] == doctest::Approx(1.0).epsilon(1e-3));

        // An index past the last output is an invalid configuration.
        constexpr auto bad = design::discrete_lqi(sys2, Q, Matrix<1, 1>{{1.0}}, {2});
        static_assert(!bad.success);
        // Integrating both outputs of a single-input plant cannot be held apart: not stabilizable.
        const auto both = design::discrete_lqi(sys2, Matrix<4, 4>::identity(), Matrix<1, 1>{{1.0}}, {0, 1});
        CHECK_FALSE(both.success);
        LQI<2, 1, 2, double, 2> idle{both};
        CHECK(idle.control(ColVec<2>{{1.0}, {1.0}}, ColVec<2>{}, ColVec<2>{})[0] == 0.0); // failed design: zero output
    }

    TEST_CASE("feedback / commit hold the integrator while the actuator is clamped") {
        const auto sys = make_plant();
        const auto r = design_lqi();
        REQUIRE(r.success);
        const double u_max = 0.2;
        const auto   clamp_u = [&](const ColVec<1>& u) { return ColVec<1>{{damp::clamp(u[0], -u_max, u_max)}}; };

        const auto run = [&](bool aware) {
            LQI<2, 1, 1, double> ctrl{r};
            ColVec<2>            x{};
            double               peak = 0.0;
            double               xi_peak = 0.0;
            for (int k = 0; k < 2000; ++k) {
                const ColVec<1> y{{x[0]}};
                const ColVec<1> u_req = ctrl.feedback(ColVec<1>{{1.0}}, y, x);
                const ColVec<1> u = clamp_u(u_req);
                ctrl.commit(aware ? u : u_req); // the unaware loop reports the unclamped command
                x = sys.A * x + sys.B * u;
                peak = damp::max(peak, x[0]);
                xi_peak = damp::max(xi_peak, damp::abs(ctrl.xi[0]));
            }
            return damp::array<double, 3>{peak, xi_peak, x[0]};
        };
        const auto naive = run(false);
        const auto aware = run(true);
        CHECK(aware[0] < naive[0]);                            // less overshoot
        CHECK(aware[1] < naive[1]);                            // the integral stays smaller
        CHECK(aware[2] == doctest::Approx(1.0).epsilon(1e-3)); // still offset-free

        // Feedforward before the clamp is not mistaken for saturation.
        LQI<2, 1, 1, double> a{r};
        LQI<2, 1, 1, double> b{r};
        const ColVec<1>      y{{0.1}};
        const ColVec<2>      x{{0.1}, {0.0}};
        const ColVec<1>      u_ff{{0.05}};
        (void)a.control(ColVec<1>{{0.2}}, y, x);
        const ColVec<1> ub = b.feedback(ColVec<1>{{0.2}}, y, x);
        b.commit(ub + u_ff, u_ff);
        CHECK(b.xi[0] == doctest::Approx(a.xi[0]));
    }
}
