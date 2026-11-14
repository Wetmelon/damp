// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/fixed_gain_observer.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

// LQG = LQR (state feedback) + Kalman filter (state estimation), tied together
// by the separation principle. design::lqg() was previously only covered
// indirectly via test_design/test_api; this exercises the full Tier 2 → Tier 3
// path including the runtime LQG controller built straight from the result
// (which feeds the Kalman filter from result.kalman.sys — the constructor
// dogfooded by lqg.hpp).
//
// @see "Optimal Control" (Anderson & Moore, 1990), §8 (separation principle).

namespace {
// Discrete double integrator at Ts = 0.1 s: x = [position, velocity].
constexpr StateSpace<2, 1, 1, double, 2, 1> make_plant() {
    return StateSpace<2, 1, 1, double, 2, 1>{
        Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}}, // A
        Matrix<2, 1>{{0.005}, {0.1}},         // B
        Matrix<1, 2>{{1.0, 0.0}},             // C: measure position only
        Matrix<1, 1>::zeros(),                // D
        Matrix<2, 2>::identity(),             // G: process noise on both states
        Matrix<1, 1>::identity(),             // H
        0.1                                   // Ts
    };
}
} // namespace

TEST_SUITE("LQG") {
    TEST_CASE("design succeeds and reports stable closed loop") {
        constexpr auto         sys = make_plant();
        constexpr auto         Q_lqr = Matrix<2, 2>::identity();
        constexpr Matrix<1, 1> R_lqr{{0.1}};
        constexpr Matrix<2, 2> Q_kf{{0.01, 0.0}, {0.0, 0.01}};
        constexpr Matrix<1, 1> R_kf{{0.1}};

        constexpr auto result = design::lqg(sys, {.Q = Q_lqr, .R = R_lqr}, design::KalmanNoise{.process = Q_kf, .measurement = R_kf});
        static_assert(result.success, "LQG design must converge at compile time");
        static_assert(result.kalman.success);

        // LQR closed loop A − BK must be inside the unit circle.
        CHECK(result.is_stable());
        // Kalman error covariance is positive on the diagonal.
        CHECK(result.kalman.P(0, 0) > 0.0);
        CHECK(result.kalman.P(1, 1) > 0.0);
    }

    TEST_CASE("runtime LQG regulates output-feedback plant to zero") {
        // The controller only sees the measured position y = x[0]; it must
        // reconstruct velocity through the Kalman filter and still drive the
        // full state to zero. This is the separation principle working end to end.
        const auto         sys = make_plant();
        const auto         Q_lqr = Matrix<2, 2>::identity();
        const Matrix<1, 1> R_lqr{{0.1}};
        const Matrix<2, 2> Q_kf{{0.01, 0.0}, {0.0, 0.01}};
        const Matrix<1, 1> R_kf{{0.1}};

        const auto result = design::lqg(sys, {.Q = Q_lqr, .R = R_lqr}, design::KalmanNoise{.process = Q_kf, .measurement = R_kf});
        REQUIRE(result.success);

        // Built straight from the result — exercises kf(result.kalman).
        auto controller = LQG{result};

        ColVec<2> x{1.0, 0.0}; // start displaced by 1 m
        for (int k = 0; k < 200; ++k) {
            const ColVec<1> y{{x[0]}}; // measure position
            controller.update(y);      // Kalman correct
            const ColVec<1> u = controller.control();
            x = sys.A * x + sys.B * u; // advance true plant
            controller.back_calculate(u);
            controller.predict(); // Kalman predict
        }

        CHECK(x[0] == doctest::Approx(0.0).epsilon(0.02)); // position regulated
        CHECK(x[1] == doctest::Approx(0.0).epsilon(0.02)); // velocity regulated
    }

    TEST_CASE("preferred step(x_ref, u_ref, y) path regulates plant to zero") {
        // Preferred self-contained tick: predict with the last command, update, control.
        const auto sys = make_plant();
        const auto result = design::lqg(
            sys, {.Q = Matrix<2, 2>::identity(), .R = 0.1},
            {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = 0.1}
        );
        REQUIRE(result.success);

        auto      controller = LQG{result};
        ColVec<2> x{1.0, 0.0};
        ColVec<1> last_u{};
        for (int k = 0; k < 200; ++k) {
            const ColVec<1> y{{x[0]}};
            last_u = controller.step(y);
            x = sys.A * x + sys.B * last_u;
        }

        CHECK(x[0] == doctest::Approx(0.0).epsilon(0.02));
        CHECK(x[1] == doctest::Approx(0.0).epsilon(0.02));
        // The tick recorded its command; the next step predicts with it.
        CHECK(controller.last_command()[0] == last_u[0]);
    }

    TEST_CASE("step with the applied input: saturation bookkeeping uses what the plant received") {
        // Clamp the raw command and hand the realized input to the next step so
        // the estimator tracks what the plant received.
        const auto sys = make_plant();
        const auto result = design::lqg(
            sys, {.Q = Matrix<2, 2>::identity(), .R = 0.1},
            {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = 0.1}
        );
        REQUIRE(result.success);

        auto         controller = LQG{result};
        ColVec<2>    x{1.0, 0.0};
        const double u_max = 0.05; // tight clamp forces saturation early on
        ColVec<1>    u_applied{};
        for (int k = 0; k < 50; ++k) {
            const ColVec<1> y{{x[0]}};
            controller.back_calculate(u_applied);
            ColVec<1> u = controller.step(y);
            CHECK(controller.applied()[0] == u_applied[0]); // predicted with last tick's clamped input
            if (u[0] > u_max) {
                u[0] = u_max;
            }
            if (u[0] < -u_max) {
                u[0] = -u_max;
            }
            u_applied = u;
            x = sys.A * x + sys.B * u;
        }
        // Still makes progress toward the origin under the clamp.
        CHECK(damp::abs(x[0]) < 1.0);
    }

    TEST_CASE("step is predict, update, control in that order; the no-input step predicts with the last command") {
        const auto sys = make_plant();
        const auto result = design::lqg(
            sys, {.Q = Matrix<2, 2>::identity(), .R = 0.1},
            {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = 0.1}
        );
        REQUIRE(result.success);

        // A clamped loop: the step and the three building blocks agree to rounding (fast-math may contract either differently).
        auto            fused = LQG{result};
        auto            split = LQG{result};
        ColVec<2>       x{1.0, 0.0};
        ColVec<1>       u_applied{};
        const ColVec<2> x_ref{0.2, 0.0};
        for (int k = 0; k < 40; ++k) {
            const ColVec<1> y{{x[0]}};
            fused.back_calculate(u_applied);
            const ColVec<1> u_fused = fused.step(x_ref, {}, y);
            split.back_calculate(u_applied);
            split.predict();
            split.update(y);
            const ColVec<1> u_split = split.control(x_ref, {});
            CHECK(u_fused[0] == doctest::Approx(u_split[0]).epsilon(1e-12));
            CHECK(fused.estimator.state()[0] == doctest::Approx(split.estimator.state()[0]).epsilon(1e-12));
            CHECK(fused.estimator.state()[1] == doctest::Approx(split.estimator.state()[1]).epsilon(1e-12));
            u_applied = ColVec<1>{{damp::clamp(u_fused[0], -0.05, 0.05)}};
            x = sys.A * x + sys.B * u_applied;
        }

        // Without an applied input the step predicts with the command it returned last.
        auto      implicit = LQG{result};
        auto      explicit_u = LQG{result};
        ColVec<2> xi{1.0, 0.0};
        ColVec<2> xe{1.0, 0.0};
        for (int k = 0; k < 40; ++k) {
            const ColVec<1> ui = implicit.step(ColVec<1>{{xi[0]}});
            explicit_u.back_calculate(explicit_u.last_command());
            const ColVec<1> ue = explicit_u.step(ColVec<1>{{xe[0]}});
            CHECK(ui[0] == doctest::Approx(ue[0]).epsilon(1e-12));
            CHECK(implicit.last_command()[0] == ui[0]);
            xi = sys.A * xi + sys.B * ui;
            xe = sys.A * xe + sys.B * ue;
        }
    }

    TEST_CASE("LQGResult::as<float>() preserves the design") {
        const auto sys = make_plant();
        const auto result = design::lqg(
            sys, {.Q = Matrix<2, 2>::identity(), .R = 0.1},
            {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = 0.1}
        );
        REQUIRE(result.success);

        const auto rf = result.as<float>();
        CHECK(rf.success);
        CHECK(rf.K(0, 0) == doctest::Approx(static_cast<float>(result.K(0, 0))));
        CHECK(rf.kalman.M(0, 0) == doctest::Approx(static_cast<float>(result.kalman.M(0, 0))));
    }

    TEST_CASE("lqgreg combines independent Kalman and LQR designs") {
        const auto sys = make_plant();
        const auto kf = design::kalman(sys, {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = 0.1});
        const auto lqr = design::lqr(sys, {.Q = Matrix<2, 2>::identity(), .R = 0.1});
        REQUIRE(kf.success);
        REQUIRE(lqr.success);

        const auto combined = design::lqg_from_parts(kf, lqr);
        CHECK(combined.success);
        CHECK(combined.status == design::DesignStatus::ok);
        // The combined result reuses the inputs verbatim.
        CHECK(combined.K(0, 0) == lqr.K(0, 0));
        CHECK(combined.kalman.M(0, 0) == kf.M(0, 0));
        CHECK(combined.Ts == kf.sys.Ts);
    }

    TEST_CASE("LQG uses designed steady-state L (not time-varying P)") {
        const auto sys = make_plant();
        const auto result = design::lqg(
            sys, {.Q = Matrix<2, 2>::identity(), .R = 0.1},
            {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = 0.1}
        );
        REQUIRE(result.success);

        auto controller = LQG{result};
        // Fixed gain matches design (the law's FixedGainObserver holds M).
        CHECK(controller.estimator.gain()(0, 0) == doctest::Approx(result.kalman.M(0, 0)));
        CHECK(controller.estimator.gain()(1, 0) == doctest::Approx(result.kalman.M(1, 0)));

        // One update applies that L: x̂⁺ = x̂ + L (y − C x̂) from zero state.
        const ColVec<1> y{{1.0}};
        controller.update(y);
        const ColVec<2> xhat = controller.estimator.state();
        CHECK(xhat[0] == doctest::Approx(result.kalman.M(0, 0) * 1.0));
        CHECK(xhat[1] == doctest::Approx(result.kalman.M(1, 0) * 1.0));
    }

    TEST_CASE("LQGResult::to_ss matches step and regulates") {
        // Current-estimator realization: one ss step equals LQG::step, and the
        // closed loop regulates a displaced plant.
        const auto sys = make_plant();
        const auto result = design::lqg(
            sys, {.Q = Matrix<2, 2>::identity(), .R = 0.1},
            {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = 0.1}
        );
        REQUIRE(result.success);

        const auto ss = result.to_ss();
        const auto direct = -(result.K * result.kalman.M);
        CHECK(ss.D(0, 0) == doctest::Approx(direct(0, 0)));
        CHECK(ss.Ts == doctest::Approx(sys.Ts));

        auto      runtime = LQG{result};
        ColVec<2> xc{{0.0}, {0.0}};
        ColVec<2> xp{{1.0}, {0.0}};
        ColVec<2> xp_rt = xp;
        for (int k = 0; k < 8; ++k) {
            const ColVec<1> y{{xp[0]}};
            const ColVec<1> u = ss.C * xc + ss.D * y;
            const ColVec<1> u_rt = runtime.step(ColVec<1>{{xp_rt[0]}});
            CHECK(u[0] == doctest::Approx(u_rt[0]));
            xp = sys.A * xp + sys.B * u;
            xp_rt = sys.A * xp_rt + sys.B * u_rt;
            xc = ss.A * xc + ss.B * y;
        }

        for (int k = 0; k < 200; ++k) {
            const ColVec<1> y{{xp[0]}};
            const ColVec<1> u = ss.C * xc + ss.D * y;
            xp = sys.A * xp + sys.B * u;
            xc = ss.A * xc + ss.B * y;
        }
        CHECK(xp[0] == doctest::Approx(0.0).epsilon(0.02));
        CHECK(xp[1] == doctest::Approx(0.0).epsilon(0.02));
    }

    TEST_CASE("design::lqg converges with nonzero gains (smoke)") {
        constexpr Matrix<2, 2> A{{0.9, 0.1}, {0.0, 0.8}};
        constexpr Matrix<2, 1> B{{0.1}, {0.2}};
        constexpr Matrix<1, 2> C{{1.0, 0.0}};
        constexpr Matrix<2, 2> Q_lqr{{1.0, 0.0}, {0.0, 1.0}};
        constexpr Matrix<1, 1> R_lqr{{1.0}};
        constexpr Matrix<2, 2> Q_kf{{0.1, 0.0}, {0.0, 0.1}};
        constexpr Matrix<1, 1> R_kf{{0.5}};

        constexpr auto result = design::lqg(
            StateSpace<2, 1, 1, double, 2, 1>{
                A, B, C, Matrix<1, 1>::zeros(), Matrix<2, 2>::identity(), Matrix<1, 1>::identity(), 1.0
            },
            {.Q = Q_lqr, .R = R_lqr}, design::KalmanNoise{.process = Q_kf, .measurement = R_kf}
        );
        static_assert(result.success);
        static_assert(result.kalman.success);
        CHECK(result.K(0, 0) != 0.0);
        CHECK(result.K(0, 1) != 0.0);
        CHECK(result.kalman.M(0, 0) != 0.0);
        CHECK(result.kalman.M(1, 0) != 0.0);
    }

    TEST_CASE("LQGResult::as<U>() conversion") {
        constexpr auto lqg_d = design::lqg(
            StateSpace<1, 1, 1, double, 1, 1>{
                Matrix<1, 1>{{1.0}},
                Matrix<1, 1>{{1.0}},
                Matrix<1, 1>{{1.0}},
                Matrix<1, 1>::zeros(),
                Matrix<1, 1>::identity(),
                Matrix<1, 1>::identity(),
                1.0
            },
            {.Q = 1.0, .R = 1.0},
            {.process = 0.1, .measurement = 0.5}
        );
        constexpr auto lqg_f = lqg_d.as<float>();
        static_assert(lqg_f.success);
        CHECK(lqg_f.success);
        CHECK(lqg_f.K(0, 0) != 0.0f);
        CHECK(lqg_f.kalman.M(0, 0) != 0.0f);
    }

    TEST_CASE("LQG deduces sizes and scalar from its design") {
        static constexpr auto res = design::lqg(
            make_plant(), {.Q = Matrix<2, 2>::identity(), .R = 0.1},
            {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = 0.1}
        );
        static constinit LQG law{res.as<float>()}; // sizes and scalar deduced from the result
        static_assert(std::is_same_v<decltype(law), LQG<2, 1, 1, float, FixedGainObserver<2, 1, 1, float>>>);
        // A double design does not narrow to the float default.
        static_assert(std::is_same_v<decltype(LQG{res}), LQG<2, 1, 1, double, FixedGainObserver<2, 1, 1, double>>>);
        CHECK(law.K(0, 0) == doctest::Approx(static_cast<float>(res.K(0, 0))));
    }

    TEST_CASE("design::lqg with TrackedOutputs holds the output at r") {
        // A damped plant (no integrator): position measured, held at r through the target.
        StateSpace<2, 1, 1, double, 2, 1> sys = make_plant();
        sys.A = Matrix<2, 2>{{1.0, 0.1}, {-0.2, 0.95}};
        const auto res = design::lqg(
            sys, {.Q = Matrix<2, 2>::identity(), .R = 0.1}, {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = 0.1},
            design::TrackedOutputs{}
        );
        REQUIRE(res.success);
        LQG law{res};
        static_assert(std::is_same_v<decltype(law), LQG<2, 1, 1, double, FixedGainObserver<2, 1, 1, double>, 1>>);
        ColVec<2>       x{};
        const ColVec<1> r{0.5};
        for (int k = 0; k < 600; ++k) {
            const ColVec<1> u = law.step(r, ColVec<1>{{x[0]}});
            x = (sys.A * x) + (sys.B * u);
        }
        CHECK(x[0] == doctest::Approx(0.5).epsilon(1e-9));
    }

    TEST_CASE("Law names a member's type by its design") {
        static constexpr auto res = design::lqg(
            make_plant(), {.Q = Matrix<2, 2>::identity(), .R = 0.1},
            {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = 0.1}
        );
        struct Axis {
            Law<res, float> law; // the design's law
            float           offset{0.0f};
        };
        static constinit Axis axis{};
        static_assert(std::is_same_v<Law<res, float>::Runtime, decltype(LQG{res.as<float>()})>);
        const LQG deduced{res.as<float>()};
        CHECK(axis.law.K(0, 0) == deduced.K(0, 0));
        CHECK(axis.law.K(0, 1) == deduced.K(0, 1));

        // A redesign at run time replaces a member of the runtime's type: a constexpr design's Law has its matrices in the type.
        Law<res, float>::Runtime replaceable{res.as<float>()};
        replaceable = LQG{design::lqg(
                              make_plant(), {.Q = Matrix<2, 2>::identity() * 10.0, .R = 0.1},
                              {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = 0.1}
        )
                              .as<float>()};
        CHECK(replaceable.K(0, 0) > deduced.K(0, 0));
    }
}

TEST_SUITE("LQG status") {
    TEST_CASE("a continuous plant is refused; the first failure is the result's status") {
        auto continuous = make_plant();
        continuous.Ts = 0.0;
        const auto refused = design::lqg(continuous, {.Q = Matrix<2, 2>::identity(), .R = 1.0}, {});
        CHECK_FALSE(refused.success);
        CHECK(refused.status == design::DesignStatus::invalid_sample_time);
        CHECK(refused.kalman.status == design::DesignStatus::invalid_sample_time);

        //! Measuring velocity only leaves the position (a mode at z = 1) undetectable.
        auto blind = make_plant();
        blind.C = Matrix<1, 2>{{0.0, 1.0}};
        const auto undetectable = design::lqg(blind, {.Q = Matrix<2, 2>::identity(), .R = 1.0}, {});
        CHECK_FALSE(undetectable.is_stable()); // a failed design reports no stable closed loop, even when its gain half succeeded
        CHECK_FALSE(undetectable.success);
        CHECK(undetectable.status == design::DesignStatus::not_detectable);
        CHECK(design::lqg_from_parts(undetectable.kalman, design::lqr(blind, {.Q = Matrix<2, 2>::identity(), .R = 1.0})).status == design::DesignStatus::not_detectable);

        constexpr auto ok = design::lqg(make_plant(), {.Q = Matrix<2, 2>::identity(), .R = 1.0}, {});
        static_assert(ok.success && ok.status == design::DesignStatus::ok && ok.Ts == 0.1);
        CHECK(ok.as<float>().status == design::DesignStatus::ok);
    }

    TEST_CASE("lqg(sys, Ts, w, noise): a continuous plant with continuous cost and a noise density") {
        constexpr StateSpace<2, 1, 1, double, 1, 1> continuous{
            .A = Matrix<2, 2>{{0.0, 1.0}, {0.0, 0.0}},
            .B = Matrix<2, 1>{{0.0}, {1.0}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .G = Matrix<2, 1>{{0.0}, {1.0}},
            .H = Matrix<1, 1>{{1.0}},
        };
        constexpr design::LQRConfig<2, 1>   w{.Q = Matrix<2, 2>::identity(), .R = 0.1};
        constexpr design::KalmanNoise<1, 1> noise{.process = 0.5, .measurement = 0.01};
        constexpr double                    Ts = 0.05;

        static constexpr auto res = design::lqg(continuous, Ts, w, noise);
        static_assert(res.success && res.Ts == Ts && res.kalman.sys.Ts == Ts);

        //! Each half is its own sampled-data design.
        constexpr auto lqr_d = design::lqr(continuous, Ts, w);
        constexpr auto kf_d = design::kalman(continuous, Ts, noise);
        static_assert(res.kalman.R == noise.measurement); // per-reading variance, no R/Ts
        static_assert(res.K == lqr_d.K);
        static_assert(res.kalman.M == kf_d.M);
        CHECK(res.is_stable());

        //! It deploys (float, constant-initialized), and regulates the discretized plant.
        static constinit auto law_f = LQG{res.as<float>()};
        CHECK(law_f.getK()(0, 0) == doctest::Approx(res.K(0, 0)).epsilon(1e-5));
        auto      law = LQG{res};
        ColVec<2> x{1.0, 0.0};
        for (int k = 0; k < 400; ++k) {
            const ColVec<1> u = law.step(ColVec<1>{{x[0]}});
            x = res.kalman.sys.A * x + res.kalman.sys.B * u;
        }
        CHECK(damp::abs(x[0]) < 1e-3);
        CHECK(damp::abs(x[1]) < 1e-3);

        //! Rows that refuse: a discrete plant (even at the same Ts), a bad Ts.
        auto discrete = continuous;
        discrete.Ts = Ts;
        const auto twice = design::lqg(discrete, Ts, w, noise);
        CHECK(twice.status == design::DesignStatus::invalid_sample_time);
        CHECK(design::lqg(continuous, 0.0, w, noise).status == design::DesignStatus::invalid_sample_time);
        CHECK(design::lqg(continuous, -Ts, w, noise).status == design::DesignStatus::invalid_sample_time);
    }

    TEST_CASE("lqg_from_parts refuses a continuous LQR or a Ts mismatch") {
        const auto sys = make_plant();
        const auto kf = design::kalman(sys, {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = 0.1});
        const auto lqr_d = design::lqr(sys, {.Q = Matrix<2, 2>::identity(), .R = 1.0});
        REQUIRE(kf.success);
        REQUIRE(lqr_d.success);

        const auto ok = design::lqg_from_parts(kf, lqr_d);
        CHECK(ok.success);
        CHECK(ok.status == design::DesignStatus::ok);

        constexpr StateSpace<2, 1, 1> continuous{
            .A = Matrix<2, 2>{{0.0, 1.0}, {0.0, 0.0}},
            .B = Matrix<2, 1>{{0.0}, {1.0}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
        };
        const auto lqr_c = design::lqr(continuous, {.Q = Matrix<2, 2>::identity(), .R = 1.0});
        REQUIRE(lqr_c.success); // CARE can succeed; LQG still cannot take it
        CHECK(lqr_c.Ts == 0.0);
        const auto refused_c = design::lqg_from_parts(kf, lqr_c);
        CHECK_FALSE(refused_c.success);
        CHECK(refused_c.status == design::DesignStatus::invalid_sample_time);

        auto lqr_mis = lqr_d;
        lqr_mis.Ts = sys.Ts * 2.0;
        const auto refused_m = design::lqg_from_parts(kf, lqr_mis);
        CHECK_FALSE(refused_m.success);
        CHECK(refused_m.status == design::DesignStatus::invalid_sample_time);
    }
}
