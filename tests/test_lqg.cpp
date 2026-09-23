// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

// LQG = LQR (state feedback) + Kalman filter (state estimation), tied together
// by the separation principle. design::discrete_lqg() was previously only covered
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

        constexpr auto result = design::discrete_lqg(sys, Q_lqr, R_lqr, Q_kf, R_kf);
        static_assert(result.success, "LQG design must converge at compile time");
        static_assert(result.lqr.success);
        static_assert(result.kalman.success);

        // LQR closed loop A − BK must be inside the unit circle.
        CHECK(result.lqr.is_stable());
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

        const auto result = design::discrete_lqg(sys, Q_lqr, R_lqr, Q_kf, R_kf);
        REQUIRE(result.success);

        // Built straight from the result — exercises kf(result.kalman).
        LQG<2, 1, 1, double, 2, 1> controller{result};

        ColVec<2> x{1.0, 0.0}; // start displaced by 1 m
        for (int k = 0; k < 200; ++k) {
            const ColVec<1> y{{x[0]}}; // measure position
            controller.update(y);      // Kalman correct
            const ColVec<1> u = controller.control();
            x = sys.A * x + sys.B * u; // advance true plant
            controller.predict(u);     // Kalman predict
        }

        CHECK(x[0] == doctest::Approx(0.0).epsilon(0.02)); // position regulated
        CHECK(x[1] == doctest::Approx(0.0).epsilon(0.02)); // velocity regulated
    }

    TEST_CASE("preferred step(y) path regulates plant to zero") {
        // Preferred self-contained tick: feedback(y) then commit(u).
        const auto sys = make_plant();
        const auto result = design::discrete_lqg(
            sys, Matrix<2, 2>::identity(), Matrix<1, 1>{{0.1}},
            Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, Matrix<1, 1>{{0.1}}
        );
        REQUIRE(result.success);

        LQG<2, 1, 1, double, 2, 1> controller{result};
        ColVec<2>                  x{1.0, 0.0};
        ColVec<1>                  last_u{};
        for (int k = 0; k < 200; ++k) {
            const ColVec<1> y{{x[0]}};
            last_u = controller.step(y);
            x = sys.A * x + sys.B * last_u;
        }

        CHECK(x[0] == doctest::Approx(0.0).epsilon(0.02));
        CHECK(x[1] == doctest::Approx(0.0).epsilon(0.02));
        // commit stored the last applied (unsaturated) input in u_prev.
        CHECK(controller.u_prev[0] == doctest::Approx(last_u[0]));
    }

    TEST_CASE("feedback/commit saturation bookkeeping uses applied input") {
        // Split path: clamp the raw command and commit the realized input so
        // the estimator tracks what the plant received.
        const auto sys = make_plant();
        const auto result = design::discrete_lqg(
            sys, Matrix<2, 2>::identity(), Matrix<1, 1>{{0.1}},
            Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, Matrix<1, 1>{{0.1}}
        );
        REQUIRE(result.success);

        LQG<2, 1, 1, double, 2, 1> controller{result};
        ColVec<2>                  x{1.0, 0.0};
        const double               u_max = 0.05; // tight clamp forces saturation early on
        for (int k = 0; k < 50; ++k) {
            const ColVec<1> y{{x[0]}};
            ColVec<1>       u = controller.feedback(y);
            if (u[0] > u_max) {
                u[0] = u_max;
            }
            if (u[0] < -u_max) {
                u[0] = -u_max;
            }
            controller.commit(u);
            CHECK(controller.u_prev[0] == doctest::Approx(u[0]));
            x = sys.A * x + sys.B * u;
        }
        // Still makes progress toward the origin under the clamp.
        CHECK(damp::abs(x[0]) < 1.0);
    }

    TEST_CASE("LQGResult::as<float>() preserves the design") {
        const auto sys = make_plant();
        const auto result = design::discrete_lqg(
            sys, Matrix<2, 2>::identity(), Matrix<1, 1>{{0.1}},
            Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, Matrix<1, 1>{{0.1}}
        );
        REQUIRE(result.success);

        const auto rf = result.as<float>();
        CHECK(rf.success);
        CHECK(rf.lqr.K(0, 0) == doctest::Approx(static_cast<float>(result.lqr.K(0, 0))));
        CHECK(rf.kalman.L(0, 0) == doctest::Approx(static_cast<float>(result.kalman.L(0, 0))));
    }

    TEST_CASE("lqgreg combines independent Kalman and LQR designs") {
        const auto sys = make_plant();
        const auto kf = design::kalman(sys, Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, Matrix<1, 1>{{0.1}});
        const auto lqr = design::discrete_lqr(sys.A, sys.B, Matrix<2, 2>::identity(), Matrix<1, 1>{{0.1}});
        REQUIRE(kf.success);
        REQUIRE(lqr.success);

        const auto combined = design::lqg_from_parts(kf, lqr);
        CHECK(combined.success);
        // The combined result reuses the inputs verbatim.
        CHECK(combined.lqr.K(0, 0) == lqr.K(0, 0));
        CHECK(combined.kalman.L(0, 0) == kf.L(0, 0));
    }

    TEST_CASE("LQG uses designed steady-state L (not time-varying P)") {
        const auto sys = make_plant();
        const auto result = design::discrete_lqg(
            sys, Matrix<2, 2>::identity(), Matrix<1, 1>{{0.1}},
            Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, Matrix<1, 1>{{0.1}}
        );
        REQUIRE(result.success);

        LQG<2, 1, 1, double, 2, 1> controller{result};
        // Fixed gain matches design (LQG stores SteadyStateKalmanFilter).
        CHECK(controller.kf.gain()(0, 0) == doctest::Approx(result.kalman.L(0, 0)));
        CHECK(controller.kf.gain()(1, 0) == doctest::Approx(result.kalman.L(1, 0)));

        // One update applies that L: x̂⁺ = x̂ + L (y − C x̂) from zero state.
        const ColVec<1> y{{1.0}};
        controller.update(y);
        const ColVec<2> xhat = controller.kf.state();
        CHECK(xhat[0] == doctest::Approx(result.kalman.L(0, 0) * 1.0));
        CHECK(xhat[1] == doctest::Approx(result.kalman.L(1, 0) * 1.0));
    }

    TEST_CASE("LQGResult::to_ss matches step and regulates") {
        // Current-estimator realization: one ss step equals LQG::step, and the
        // closed loop regulates a displaced plant.
        const auto sys = make_plant();
        const auto result = design::discrete_lqg(
            sys, Matrix<2, 2>::identity(), Matrix<1, 1>{{0.1}},
            Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, Matrix<1, 1>{{0.1}}
        );
        REQUIRE(result.success);

        const auto ss = result.to_ss();
        const auto direct = -(result.lqr.K * result.kalman.L);
        CHECK(ss.D(0, 0) == doctest::Approx(direct(0, 0)));
        CHECK(ss.Ts == doctest::Approx(sys.Ts));

        LQG<2, 1, 1, double, 2, 1> runtime{result};
        ColVec<2>                  xc{{0.0}, {0.0}};
        ColVec<2>                  xp{{1.0}, {0.0}};
        ColVec<2>                  xp_rt = xp;
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

        constexpr auto result = design::discrete_lqg(
            StateSpace<2, 1, 1, double, 2, 1>{
                A, B, C, Matrix<1, 1>::zeros(), Matrix<2, 2>::identity(), Matrix<1, 1>::identity(), 1.0
            },
            Q_lqr, R_lqr, Q_kf, R_kf
        );
        static_assert(result.success);
        static_assert(result.lqr.success);
        static_assert(result.kalman.success);
        CHECK(result.lqr.K(0, 0) != 0.0);
        CHECK(result.lqr.K(0, 1) != 0.0);
        CHECK(result.kalman.L(0, 0) != 0.0);
        CHECK(result.kalman.L(1, 0) != 0.0);
    }

    TEST_CASE("LQGResult::as<U>() conversion") {
        constexpr auto lqg_d = design::discrete_lqg(
            StateSpace<1, 1, 1, double, 1, 1>{
                Matrix<1, 1>{{1.0}},
                Matrix<1, 1>{{1.0}},
                Matrix<1, 1>{{1.0}},
                Matrix<1, 1>::zeros(),
                Matrix<1, 1>::identity(),
                Matrix<1, 1>::identity(),
                1.0
            },
            Matrix<1, 1>{{1.0}},
            Matrix<1, 1>{{1.0}},
            Matrix<1, 1>{{0.1}},
            Matrix<1, 1>{{0.5}}
        );
        constexpr auto lqg_f = lqg_d.as<float>();
        static_assert(lqg_f.success);
        CHECK(lqg_f.success);
        CHECK(lqg_f.lqr.K(0, 0) != 0.0f);
        CHECK(lqg_f.kalman.L(0, 0) != 0.0f);
    }
}
