// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "damp/controllers/offset_free_lqg.hpp"
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

TEST_SUITE("OFLQG") {
    TEST_CASE("design succeeds") {
        constexpr auto sys = make_plant();
        constexpr auto result = design::discrete_oflqg(
            sys, Matrix<2, 2>::identity(), Matrix<1, 1>{{0.1}},
            {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = Matrix<1, 1>{{0.1}}}
        );
        static_assert(result.success);
        CHECK(result.success);
        CHECK(result.lqr.K(0, 0) != 0.0);
        CHECK(result.kalman.L(0, 0) != 0.0);
    }

    TEST_CASE("tracks a setpoint under a constant input load") {
        const auto sys = make_plant();
        const auto result = design::discrete_oflqg(
            sys, Matrix<2, 2>::identity(), Matrix<1, 1>{{0.1}},
            design::OffsetFreeNoise{Matrix<2, 2>{{0.05, 0.0}, {0.0, 0.05}}, Matrix<1, 1>{{0.05}}, Matrix<1, 1>{{1.0}}}
        );
        REQUIRE(result.success);

        OFLQG<2, 1, 1, double> ctrl{result};
        ColVec<2>              x{0.0, 0.0};
        const ColVec<1>        r{{1.0}};
        const double           w_true = 0.4;
        for (int k = 0; k < 400; ++k) {
            const ColVec<1> y{{x[0]}};
            const ColVec<1> u = ctrl.control(r, y);
            x = sys.A * x + sys.B * (u + ColVec<1>{{w_true}});
        }
        CHECK(x[0] == doctest::Approx(1.0).epsilon(0.05));
        CHECK(ctrl.disturbance_estimate()(0) == doctest::Approx(w_true).epsilon(0.15));
    }

    TEST_CASE("as<float> preserves success") {
        constexpr auto result = design::discrete_oflqg(
            make_plant(), Matrix<2, 2>::identity(), Matrix<1, 1>{{0.1}},
            {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = Matrix<1, 1>{{0.1}}}
        );
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
        StateSpace<4, 1, 3> sys{.A = disc->A, .C = disc->C};
        sys.B = disc->B.template block<4, 1>(0, 0);
        sys.Ts = Ts;

        design::OffsetFreeModel<4, 1, 3, 1, 1> model{};
        model.Bd = disc->B.template block<4, 1>(0, 1);
        model.Cd(2, 0) = 1.0 / m;
        model.Cz(0, 2) = 1.0; // hold the load position

        Matrix<4, 4> Q = Matrix<4, 4>::identity();
        Q(2, 2) = 1.0e4;
        const auto res = design::discrete_oflqg(
            sys, model, Q, Matrix<1, 1>{{1.0e-4}}, design::OffsetFreeNoise{Matrix<4, 4>::identity() * 1.0e-8, Matrix<3, 3>{{1.0e-8, 0.0, 0.0}, {0.0, 1.0e-8, 0.0}, {0.0, 0.0, 1.0e-2}}, Matrix<1, 1>{{1.0}}}
        );
        REQUIRE(res.success);

        OFLQG<4, 1, 3, double, 1, 1> ctrl{res};
        ColVec<4>                    x{};
        const ColVec<1>              r{{0.01}};
        const double                 load_force = -0.8;
        for (int i = 0; i < 4000; ++i) {
            const ColVec<3> y = (sys.C * x) + (model.Cd * ColVec<1>{{load_force}});
            const ColVec<1> u = ctrl.control(r, y);
            x = (sys.A * x) + (sys.B * u) + (model.Bd * ColVec<1>{{load_force}});
        }
        CHECK(x[2] == doctest::Approx(0.01).epsilon(1e-3));                                 // load held at r
        CHECK(ctrl.disturbance_estimate()(0) == doctest::Approx(load_force).epsilon(0.02)); // and the force found
    }

    TEST_CASE("tracking a constant reference is the setpoint form") {
        // Same design, same plant and load: control(y, x*, u*) with x* = [1, 0], u* = 0 must
        // produce exactly the inputs control(r = 1, y) does (re-anchoring makes the two filters equal).
        const auto sys = make_plant();
        const auto result = design::discrete_oflqg(
            sys, Matrix<2, 2>::identity(), Matrix<1, 1>{{0.1}}, design::OffsetFreeNoise{Matrix<2, 2>{{0.05, 0.0}, {0.0, 0.05}}, Matrix<1, 1>{{0.05}}, Matrix<1, 1>{{1.0}}}
        );
        REQUIRE(result.success);
        OFLQG<2, 1, 1, double> setpoint{result};
        OFLQG<2, 1, 1, double> tracking{result};
        ColVec<2>              xs{};
        ColVec<2>              xt{};
        double                 worst = 0.0;
        for (int k = 0; k < 300; ++k) {
            const ColVec<1> us = setpoint.control(ColVec<1>{{1.0}}, ColVec<1>{{xs[0]}});
            const ColVec<1> ut = tracking.control(ColVec<1>{{xt[0]}}, ColVec<2>{1.0, 0.0}, ColVec<1>{});
            worst = std::max(worst, std::abs(us[0] - ut[0]));
            xs = sys.A * xs + sys.B * (us + ColVec<1>{{0.4}});
            xt = sys.A * xt + sys.B * (ut + ColVec<1>{{0.4}});
        }
        CHECK(worst < 1e-9);
    }

    TEST_CASE("tracks a moving reference under a constant load") {
        const auto sys = make_plant();
        const auto result = design::discrete_oflqg(
            sys, Matrix<2, 2>{{100.0, 0.0}, {0.0, 1.0}}, Matrix<1, 1>{{0.01}}, design::OffsetFreeNoise{Matrix<2, 2>{{0.05, 0.0}, {0.0, 0.05}}, Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{1.0}}}
        );
        REQUIRE(result.success);
        OFLQG<2, 1, 1, double> ctrl{result};
        ColVec<2>              x{};
        const double           w = 0.3; // [rad/s]
        const double           Ts = 0.1;
        double                 worst = 0.0;
        for (int k = 0; k < 2000; ++k) {
            const double    t = k * Ts;
            const ColVec<2> x_ref{std::sin(w * t), w * std::cos(w * t)};
            const ColVec<1> u_ref{{-w * w * std::sin(w * t)}}; // the reference's own acceleration
            const ColVec<1> u = ctrl.control(ColVec<1>{{x[0]}}, x_ref, u_ref);
            x = sys.A * x + sys.B * (u + ColVec<1>{{0.4}});
            if (t > 100.0) {
                worst = std::max(worst, std::abs(x[0] - std::sin(w * (t + Ts))));
            }
        }
        CHECK(worst < 0.01);                                                         // on the sine, not lagging it
        CHECK(ctrl.disturbance_estimate()(0) == doctest::Approx(0.4).epsilon(0.05)); // and the load found
    }
}
