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

        constexpr auto result = design::discrete_lqgi(sys, Q_aug, R, {.process = Q_kf, .measurement = R_kf});
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

        const auto result = design::discrete_lqgi(sys, Q_aug, R, {.process = Q_kf, .measurement = R_kf});
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
        const auto result = design::discrete_lqgi(sys, make_Q_aug(), Matrix<1, 1>{{0.1}}, {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = Matrix<1, 1>{{0.1}}});
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
        const auto result = design::discrete_lqgi(sys, make_Q_aug(), Matrix<1, 1>{{0.1}}, {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = Matrix<1, 1>{{0.1}}});
        REQUIRE(result.success);

        const auto ss = result.to_ss();            // [x̂⁻; xi; u_prev], in [r; y], out u
        CHECK(ss.D(0, 0) == doctest::Approx(0.0)); // no direct term from r
        CHECK(ss.Ts == doctest::Approx(sys.Ts));

        LQGI<2, 1, 1, double, 2, 1> runtime{result};
        ColVec<4>                   xc{}; // NX + NY + NU
        ColVec<2>                   xp{0.0, 0.0};
        ColVec<2>                   xp_rt = xp;
        const ColVec<1>             r{1.0};
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
        const auto result = design::discrete_lqgi(sys, make_Q_aug(), Matrix<1, 1>{{0.1}}, {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = Matrix<1, 1>{{0.1}}});
        REQUIRE(result.success);

        const auto rf = result.as<float>();
        CHECK(rf.success);
        CHECK(rf.lqi.K(0, 0) == doctest::Approx(static_cast<float>(result.lqi.K(0, 0))));
        CHECK(rf.kalman.L(0, 0) == doctest::Approx(static_cast<float>(result.kalman.L(0, 0))));
    }

    TEST_CASE("integrating a subset of the sensors: same argument shape as discrete_lqi") {
        // Position and velocity sensors; integrate the position only.
        constexpr StateSpace<2, 1, 2, double, 2, 2> sys2{
            Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}}, Matrix<2, 1>{{0.005}, {0.1}},
            Matrix<2, 2>::identity(), Matrix<2, 1>{}, Matrix<2, 2>::identity(), Matrix<2, 2>::identity()
        };
        const Matrix<2, 2> Q_kf{{0.01, 0.0}, {0.0, 0.01}};
        const Matrix<2, 2> R_kf{{0.1, 0.0}, {0.0, 0.1}};
        const auto         res = design::discrete_lqgi(sys2, make_Q_aug(), Matrix<1, 1>{{0.1}}, {.process = Q_kf, .measurement = R_kf}, {0});
        REQUIRE(res.success);
        const auto lqi_only = design::discrete_lqi(sys2, make_Q_aug(), Matrix<1, 1>{{0.1}}, {0});
        CHECK(res.lqi.K(0, 2) == doctest::Approx(lqi_only.K(0, 2)));

        // The compensator realization takes [r (NI); y (NY)] and matches the runtime.
        const auto      ss = res.to_ss();
        LQGI            runtime{res};
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

        const auto bad = design::discrete_lqgi(sys2, make_Q_aug(), Matrix<1, 1>{{0.1}}, {.process = Q_kf, .measurement = R_kf}, {5});
        CHECK_FALSE(bad.success);
    }

    TEST_CASE("commit holds the integrator against the clamp") {
        const auto result = design::discrete_lqgi(make_plant(), make_Q_aug(), Matrix<1, 1>{{0.1}}, {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = Matrix<1, 1>{{0.1}}});
        REQUIRE(result.success);
        const auto   sys = make_plant();
        const double u_max = 0.2;
        const auto   run = [&](bool aware) {
            LQGI<2, 1, 1, double, 2, 1> ctrl{result};
            ColVec<2>                   x{};
            double                      peak = 0.0;
            for (int k = 0; k < 2000; ++k) {
                const ColVec<1> y{{x[0]}};
                const ColVec<1> u_req = ctrl.feedback(ColVec<1>{{1.0}}, y);
                const ColVec<1> u{{damp::clamp(u_req[0], -u_max, u_max)}};
                // The unaware loop hides the clamp from the integrator (it still feeds the estimator u).
                if (aware) {
                    ctrl.commit(u);
                } else {
                    ctrl.lqi.commit(u_req);
                    ctrl.kf.predict(u);
                    ctrl.u_prev = u;
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
}
