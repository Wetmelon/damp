// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "damp/backend.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/math/math.hpp"
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
}
