// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>
#include <numbers>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/pr.hpp"
#include "damp/design/synthesis.hpp"
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

    TEST_CASE("lqi_bundle / lqgi_bundle integrate a subset of the outputs, same argument as discrete_lqi") {
        // Position and velocity sensors on a double integrator; integrate the position only.
        constexpr StateSpace<2, 1, 2, double, 2, 2> sys{
            .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
            .B = Matrix<2, 1>{{0.005}, {0.1}},
            .C = Matrix<2, 2>::identity(),
            .G = Matrix<2, 2>::identity(),
            .H = Matrix<2, 2>::identity(),
            .Ts = 0.1,
        };
        constexpr Matrix<3, 3> Q_aug{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 10.0}};
        constexpr Matrix<1, 1> R{{0.1}};
        constexpr Matrix<1, 2> pick_position{{1.0, 0.0}};

        constexpr auto by_index = design::lqi_bundle(sys, Q_aug, R, {0});
        static_assert(by_index.success);
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(by_index.runtime)>, LQI<2, 1, 2, float, 1>>);
        const auto by_selector = design::lqi_bundle(sys, Q_aug, R, pick_position);
        const auto direct = design::discrete_lqi(sys, Q_aug, R, {0});
        REQUIRE(by_selector.success);
        for (size_t j = 0; j < 3; ++j) {
            CHECK(by_index.design.K(0, j) == doctest::Approx(direct.K(0, j)));
            CHECK(by_selector.design.K(0, j) == doctest::Approx(direct.K(0, j)));
        }
        // The servo model integrates r − y₀ only: the integrator row reads −C₀ = [−1, 0].
        const auto& servo = by_index.models.augmented_servo_closed_loop;
        CHECK(servo.A(2, 0) == doctest::Approx(-1.0));
        CHECK(servo.A(2, 1) == doctest::Approx(0.0));
        CHECK(servo.A(2, 2) == doctest::Approx(1.0));

        // The bundled runtime takes one reference and both measurements, and holds y₀ at r.
        auto             lqi = by_index.runtime;
        ColVec<2, float> x{};
        for (int k = 0; k < 600; ++k) {
            const auto u = lqi.control(ColVec<1, float>{1.0F}, x, x);
            x = sys.A.as<float>() * x + sys.B.as<float>() * u;
        }
        CHECK(x(0) == doctest::Approx(1.0F).epsilon(1e-3));

        // LQGI: the same argument, an estimator underneath.
        const Matrix<2, 2> Q_kf{{0.01, 0.0}, {0.0, 0.01}};
        const Matrix<2, 2> R_kf{{0.1, 0.0}, {0.0, 0.1}};
        const auto         lqgi_index = design::lqgi_bundle(sys, Q_aug, R, Q_kf, R_kf, {0});
        const auto         lqgi_selector = design::lqgi_bundle(sys, Q_aug, R, Q_kf, R_kf, pick_position);
        REQUIRE(lqgi_index.success);
        REQUIRE(lqgi_selector.success);
        CHECK(lqgi_index.design.lqi.K(0, 2) == doctest::Approx(direct.K(0, 2)));
        CHECK(lqgi_selector.design.lqi.K(0, 2) == doctest::Approx(direct.K(0, 2)));
        CHECK(lqgi_index.models.augmented_servo_closed_loop.A(2, 1) == doctest::Approx(0.0));
        auto      lqgi = lqgi_index.runtime;
        ColVec<2> xp{};
        for (int k = 0; k < 600; ++k) {
            const auto u = lqgi.step(ColVec<1, float>{1.0F}, (sys.C * xp).as<float>());
            xp = sys.A * xp + sys.B * u.as<double>();
        }
        CHECK(xp(0) == doctest::Approx(1.0).epsilon(0.02));

        // An index past the last output is an invalid configuration, in either bundle.
        CHECK_FALSE(design::lqi_bundle(sys, Q_aug, R, {2}).success);
        CHECK_FALSE(design::lqgi_bundle(sys, Q_aug, R, Q_kf, R_kf, {2}).success);
    }
}
