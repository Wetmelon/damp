// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <numbers>

#include "damp/backend.hpp"
#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/pr.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/synthesis.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

TEST_SUITE("Workflow Glue") {
    TEST_CASE("lqg_bundle builds design analysis and runtime bundles") {
        StateSpace<2, 1, 1, double, 2, 1> sys{
            .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
            .B = Matrix<2, 1>{{0.005}, {0.1}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .D = Matrix<1, 1>::zeros(),
            .G = Matrix<2, 2>::identity(),
            .H = Matrix<1, 1>::identity(),
            .Ts = 0.1
        };

        const Matrix<2, 2> Q_lqr = Matrix<2, 2>::identity();
        const Matrix<1, 1> R_lqr{{0.1}};
        const Matrix<2, 2> Q_kf{{0.01, 0.0}, {0.0, 0.01}};
        const Matrix<1, 1> R_kf{{0.1}};

        const auto artifacts = design::lqg_bundle(sys, Q_lqr, R_lqr, Q_kf, R_kf);

        CHECK(artifacts.success);
        CHECK(artifacts.design.success);
        CHECK(artifacts.models.state_feedback_closed_loop.A(0, 0) != doctest::Approx(sys.A(0, 0)));
        CHECK(artifacts.models.observer_error_dynamics(0, 0) != doctest::Approx(sys.A(0, 0)));

        auto                   runtime = artifacts.runtime;
        const ColVec<1, float> y{1.0f};
        const auto             u = runtime.step(y);
        CHECK(u(0, 0) != doctest::Approx(0.0f));
    }

    TEST_CASE("lqg_pr_bundle adds SISO PR internal model") {
        StateSpace<2, 1, 1, double, 2, 1> sys{
            .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
            .B = Matrix<2, 1>{{0.005}, {0.1}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .D = Matrix<1, 1>::zeros(),
            .G = Matrix<2, 2>::identity(),
            .H = Matrix<1, 1>::identity(),
            .Ts = 0.001
        };

        const Matrix<2, 2> Q_lqr = Matrix<2, 2>::identity();
        const Matrix<1, 1> R_lqr{{0.1}};
        const Matrix<2, 2> Q_kf{{0.01, 0.0}, {0.0, 0.01}};
        const Matrix<1, 1> R_kf{{0.1}};

        const auto pr = design::pr(
            0.0,
            10.0,
            2.0 * std::numbers::pi,
            5.0,
            0.001
        );

        const auto artifacts = design::lqg_pr_bundle(sys, Q_lqr, R_lqr, Q_kf, R_kf, pr);

        CHECK(artifacts.success);
        CHECK(artifacts.runtime_pr.pr_design.Ki == doctest::Approx(10.0f));

        auto       runtime = artifacts.runtime_pr;
        const auto u = runtime.step(1.0f, 0.0f);
        CHECK(u(0, 0) != doctest::Approx(0.0f));
    }

    TEST_CASE("lqi_bundle builds servo artifacts and a ready-to-run LQI") {
        StateSpace<2, 1, 1> sys{
            .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
            .B = Matrix<2, 1>{{0.005}, {0.1}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .D = Matrix<1, 1>::zeros(),
            .Ts = 0.1
        };

        Matrix<3, 3> Q_aug{};
        Q_aug(0, 0) = 1.0;
        Q_aug(1, 1) = 1.0;
        Q_aug(2, 2) = 10.0;
        const Matrix<1, 1> R{{0.1}};

        const auto artifacts = design::lqi_bundle(sys, Q_aug, R);

        CHECK(artifacts.success);
        CHECK(artifacts.design.success);

        auto                   runtime = artifacts.runtime;
        const ColVec<2, float> x{1.0f, 0.0f};
        const ColVec<1, float> y{1.0f};
        const ColVec<1, float> r{0.0f};
        const auto             u = runtime.control(r, y, x);

        CHECK(u(0, 0) != doctest::Approx(0.0f));
    }

    TEST_CASE("lqgi_bundle builds servo observer artifacts and a ready-to-run LQGI") {
        StateSpace<2, 1, 1, double, 2, 1> sys{
            .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
            .B = Matrix<2, 1>{{0.005}, {0.1}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .D = Matrix<1, 1>::zeros(),
            .G = Matrix<2, 2>::identity(),
            .H = Matrix<1, 1>::identity(),
            .Ts = 0.1
        };

        Matrix<3, 3> Q_aug{};
        Q_aug(0, 0) = 1.0;
        Q_aug(1, 1) = 1.0;
        Q_aug(2, 2) = 10.0;

        const Matrix<1, 1> R{{0.1}};
        const Matrix<2, 2> Q_kf{{0.01, 0.0}, {0.0, 0.01}};
        const Matrix<1, 1> R_kf{{0.1}};

        const auto artifacts = design::lqgi_bundle(sys, Q_aug, R, Q_kf, R_kf);

        CHECK(artifacts.success);
        CHECK(artifacts.design.success);
        CHECK(artifacts.models.observer_error_dynamics(0, 0) != doctest::Approx(sys.A(0, 0)));

        auto                   runtime = artifacts.runtime;
        const ColVec<1, float> y{0.0f};
        const ColVec<1, float> r{1.0f};

        // u[k] = -K·[x̂[k]; xi[k]] uses the *current* integral state (the DARE augmentation
        // convention), so the first tick from a zero estimate and zero xi commands nothing;
        // xi picks up the error and the second tick drives the plant.
        CHECK(runtime.step(r, y)(0, 0) == doctest::Approx(0.0f));
        CHECK(runtime.step(r, y)(0, 0) != doctest::Approx(0.0f));
    }

    TEST_CASE("lqi_bundle drives steady-state tracking error to zero") {
        StateSpace<2, 1, 1> sys{
            .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
            .B = Matrix<2, 1>{{0.005}, {0.1}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .D = Matrix<1, 1>::zeros(),
            .Ts = 0.1
        };

        Matrix<3, 3> Q_aug{};
        Q_aug(0, 0) = 1.0;
        Q_aug(1, 1) = 1.0;
        Q_aug(2, 2) = 10.0;
        const Matrix<1, 1> R{{0.1}};

        const auto artifacts = design::lqi_bundle(sys, Q_aug, R);
        REQUIRE(artifacts.success);

        // Closed-loop simulation against the (linear, discrete) design plant.
        // The integral state uses a unit discrete integrator, so the runtime
        // must accumulate xi += (r - y) with no Ts scaling; if it did not,
        // integral action would be ~Ts-times too weak and y would never reach r.
        const auto Af = sys.A.as<float>();
        const auto Bf = sys.B.as<float>();
        const auto Cf = sys.C.as<float>();

        auto             runtime = artifacts.runtime;
        ColVec<2, float> x{0.0f, 0.0f};
        const float      r = 1.0f;
        float            y = 0.0f;

        for (int k = 0; k < 5000; ++k) {
            y = (Cf * x)(0, 0);
            const auto u = runtime.control(ColVec<1, float>{r}, ColVec<1, float>{y}, x);
            x = Af * x + Bf * u;
        }

        CHECK(y == doctest::Approx(1.0f).epsilon(0.01));
    }

    TEST_CASE("workflow artifacts synthesize at compile time") {
        constexpr auto artifacts = [] {
            StateSpace<2, 1, 1> sys{
                .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
                .B = Matrix<2, 1>{{0.005}, {0.1}},
                .C = Matrix<1, 2>{{1.0, 0.0}},
                .D = Matrix<1, 1>::zeros(),
                .Ts = 0.1
            };
            Matrix<3, 3> Q_aug{};
            Q_aug(0, 0) = 1.0;
            Q_aug(1, 1) = 1.0;
            Q_aug(2, 2) = 10.0;
            const Matrix<1, 1> R{{0.1}};
            return design::lqi_bundle(sys, Q_aug, R);
        }();

        static_assert(artifacts.success, "LQI workflow must synthesize at compile time");
        CHECK(artifacts.success);
    }
}
