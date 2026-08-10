// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

// LQGI = LQI (state feedback with integral action) + Kalman estimator. The
// integral state drives steady-state output error to zero for constant
// references/disturbances. design::discrete_lqgi() was only covered indirectly;
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
        Matrix<1, 1>::identity()
    };
}

// Augmented-state cost: penalize position, velocity, and the integral state.
constexpr Matrix<3, 3> make_Q_aug() {
    Matrix<3, 3> Q{};
    Q(0, 0) = 1.0;  // position
    Q(1, 1) = 1.0;  // velocity
    Q(2, 2) = 10.0; // integral of tracking error — weight it to pull SS error to 0
    return Q;
}
} // namespace

TEST_SUITE("LQGI") {
    TEST_CASE("design succeeds with integral augmentation") {
        constexpr auto         sys = make_plant();
        constexpr auto         Q_aug = make_Q_aug();
        constexpr Matrix<1, 1> R{{0.1}};
        constexpr Matrix<2, 2> Q_kf{{0.01, 0.0}, {0.0, 0.01}};
        constexpr Matrix<1, 1> R_kf{{0.1}};

        constexpr auto result = design::discrete_lqgi(sys, Q_aug, R, Q_kf, R_kf);
        static_assert(result.success);
        static_assert(result.lqi.success);
        static_assert(result.kalman.success);

        // Gain has a column per augmented state: [x(2) ; xi(1)] = 3.
        CHECK(result.lqi.K(0, 0) != 0.0);
        CHECK(result.lqi.K(0, 2) != 0.0); // integral gain is non-trivial
    }

    TEST_CASE("runtime LQGI tracks a constant reference with zero steady-state error") {
        const auto         sys = make_plant();
        const auto         Q_aug = make_Q_aug();
        const Matrix<1, 1> R{{0.1}};
        const Matrix<2, 2> Q_kf{{0.01, 0.0}, {0.0, 0.01}};
        const Matrix<1, 1> R_kf{{0.1}};

        const auto result = design::discrete_lqgi(sys, Q_aug, R, Q_kf, R_kf);
        REQUIRE(result.success);

        LQGI<2, 1, 1, double, 2, 1> controller{result}; // exercises kf(result.kalman)

        const double r = 1.0;     // position setpoint
        ColVec<2>    x{0.0, 0.0}; // true plant state
        double       xi = 0.0;    // integral of (r - y)

        for (int k = 0; k < 600; ++k) {
            const ColVec<1> y{x[0]};
            controller.update(y);

            // Augmented state fed to LQI is [x_hat; xi].
            const auto      xhat = controller.kf.state();
            const ColVec<3> x_aug{xhat[0], xhat[1], xi};
            const ColVec<1> u = controller.control(x_aug);

            x = sys.A * x + sys.B * u;
            controller.predict(u);
            xi += (r - x[0]); // accumulate tracking error (matches -C augmentation sign convention)
        }

        CHECK(x[0] == doctest::Approx(r).epsilon(0.02)); // output reaches the reference
    }

    TEST_CASE("runtime LQGI tracks via the controller's own integrator (control(r, y))") {
        const auto sys = make_plant();
        const auto result = design::discrete_lqgi(sys, make_Q_aug(), Matrix<1, 1>{{0.1}}, Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, Matrix<1, 1>{{0.1}});
        REQUIRE(result.success);

        LQGI<2, 1, 1, double, 2, 1> controller{result};
        const ColVec<1>             r{1.0};
        ColVec<2>                   x{0.0, 0.0};
        for (int k = 0; k < 600; ++k) {
            const ColVec<1> y{x[0]};
            controller.update(y);
            const ColVec<1> u = controller.control(r, y); // pulls x̂ and owns xi
            x = sys.A * x + sys.B * u;
            controller.predict(u);
        }
        CHECK(x[0] == doctest::Approx(r[0]).epsilon(0.02));

        controller.reset();
        CHECK(controller.lqi.xi[0] == doctest::Approx(0.0));
    }

    TEST_CASE("LQGIResult::to_ss compensator tracks a reference") {
        // Validates the [r; y] -> u realization (steady-state L) by closing the loop.
        const auto sys = make_plant();
        const auto result = design::discrete_lqgi(sys, make_Q_aug(), Matrix<1, 1>{{0.1}}, Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, Matrix<1, 1>{{0.1}});
        REQUIRE(result.success);

        const auto ss = result.to_ss(); // StateSpace<3, 2, 1>: in [r;y], out u, state [x̂;xi]
        CHECK(ss.D(0, 0) == doctest::Approx(0.0));
        CHECK(ss.Ts == doctest::Approx(sys.Ts));

        ColVec<3>    xc{0.0, 0.0, 0.0}; // [x̂; xi]
        ColVec<2>    xp{0.0, 0.0};      // plant
        const double ref = 1.0;
        for (int k = 0; k < 600; ++k) {
            const ColVec<2> in{ref, xp[0]}; // [r; y]
            const ColVec<1> u = ss.C * xc + ss.D * in;
            xp = sys.A * xp + sys.B * u;
            xc = ss.A * xc + ss.B * in;
        }
        CHECK(xp[0] == doctest::Approx(ref).epsilon(0.02));
    }

    TEST_CASE("LQGIResult::as<float>() preserves the design") {
        const auto sys = make_plant();
        const auto result = design::discrete_lqgi(sys, make_Q_aug(), Matrix<1, 1>{{0.1}}, Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, Matrix<1, 1>{{0.1}});
        REQUIRE(result.success);

        const auto rf = result.as<float>();
        CHECK(rf.success);
        CHECK(rf.lqi.K(0, 0) == doctest::Approx(static_cast<float>(result.lqi.K(0, 0))));
        CHECK(rf.kalman.L(0, 0) == doctest::Approx(static_cast<float>(result.kalman.L(0, 0))));
    }
}
