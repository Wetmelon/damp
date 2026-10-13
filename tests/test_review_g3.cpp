// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Regression tests for the g3 review (controllers + motor modulation).

#include <cmath>
#include <numbers>

#include "damp/controllers/esc.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/controllers/offset_free_lqg.hpp"
#include "damp/controllers/repetitive.hpp"
#include "damp/controllers/smc.hpp"
#include "damp/controllers/stsmc.hpp"
#include "damp/math/transforms.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/motor/modulation.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using namespace damp::motor;

TEST_SUITE("Review g3") {
    TEST_CASE("SMC (r, y) regulates a positive-gain double integrator") {
        // ÿ = b·u with b = 2 > 0: SMC must push y toward r, not away.
        const double          b = 2.0;
        const double          Ts = 1e-3;
        SMCController<double> ctrl(design::smc(5.0, 5.0, b), Ts, 0.05);
        double                y = 0.0;
        double                v = 0.0;
        for (int k = 0; k < 6000; ++k) {
            const double u = ctrl.control(1.0, y);
            v += Ts * b * u;
            y += Ts * v;
        }
        CHECK(std::abs(y - 1.0) < 0.05);
    }

    TEST_CASE("STSMC (r, y) regulates a unit-gain double integrator under a Lipschitz load") {
        const double            Ts = 1e-3;
        STSMCController<double> ctrl(design::stsmc(2.0, 5.0), Ts);
        double                  y = 0.0;
        double                  v = 0.0;
        for (int k = 0; k < 10000; ++k) {
            const double u = ctrl.control(1.0, y);
            v += Ts * (u + (0.5 * std::sin(k * Ts)));
            y += Ts * v;
        }
        CHECK(std::abs(y - 1.0) < 0.02);
    }

    TEST_CASE("THIPWM stays linear up to the SVPWM limit Vdc/sqrt(3)") {
        const double Vdc = 100.0;
        const double Vpk = 0.57 * Vdc; // > Vdc/2 (SPWM limit), < Vdc/√3
        for (int deg = 0; deg < 360; ++deg) {
            const double            th = deg * std::numbers::pi / 180.0;
            const AlphaBeta<double> v = {.alpha = Vpk * std::cos(th), .beta = Vpk * std::sin(th)};
            CHECK_FALSE(modulation_duty_cycles<double>(v, Vdc, PwmScheme::Thipwm).is_clipped);
        }
    }

    TEST_CASE("DPWM3 clamps each phase in two 30-degree windows off its peak; DPWM0 leads, DPWM2 lags") {
        const double Vdc = 100.0;
        const double Vpk = 40.0;
        int          high_deg = 0;
        for (int deg = -180; deg < 180; ++deg) {
            const double            th = (deg + 0.5) * std::numbers::pi / 180.0;
            const AlphaBeta<double> v = {.alpha = Vpk * std::cos(th), .beta = Vpk * std::sin(th)};
            const auto              d = modulation_duty_cycles<double>(v, Vdc, PwmScheme::Dpwm3);
            const bool              a_high = d.duties[0] > 1.0 - 1e-12;
            const bool              expected = (deg >= 30 && deg < 60) || (deg >= -60 && deg < -30);
            CHECK(a_high == expected);
            high_deg += a_high ? 1 : 0;
        }
        CHECK(high_deg == 60);

        auto a_high_at = [&](PwmScheme s, double deg) {
            const double            th = deg * std::numbers::pi / 180.0;
            const AlphaBeta<double> v = {.alpha = Vpk * std::cos(th), .beta = Vpk * std::sin(th)};
            return modulation_duty_cycles<double>(v, Vdc, s).duties[0] > 1.0 - 1e-12;
        };
        CHECK(a_high_at(PwmScheme::Dpwm0, -30.0)); // window [-60, 0]: before the peak
        CHECK_FALSE(a_high_at(PwmScheme::Dpwm0, 30.0));
        CHECK(a_high_at(PwmScheme::Dpwm2, 30.0)); // window [0, 60]: after the peak
        CHECK_FALSE(a_high_at(PwmScheme::Dpwm2, -30.0));
    }

    TEST_CASE("OFLQG is offset-free with feedthrough D != 0") {
        const StateSpace<2, 1, 1, double, 2, 1> sys{
            Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
            Matrix<2, 1>{{0.005}, {0.1}},
            Matrix<1, 2>{{1.0, 0.0}},
            Matrix<1, 1>{{0.5}},
            Matrix<2, 2>::identity(),
            Matrix<1, 1>::identity(),
            0.1
        };
        const auto result = design::discrete_oflqg(
            sys, Matrix<2, 2>::identity(), Matrix<1, 1>{{0.1}}, design::OffsetFreeNoise{Matrix<2, 2>{{0.05, 0.0}, {0.0, 0.05}}, Matrix<1, 1>{{0.05}}, Matrix<1, 1>{{1.0}}}
        );
        REQUIRE(result.success);

        OFLQG<2, 1, 1, double> ctrl{result};
        ColVec<2>              x{0.0, 0.0};
        ColVec<1>              u_held{};
        const ColVec<1>        r{{1.0}};
        double                 y = 0.0;
        for (int k = 0; k < 600; ++k) {
            // Sampled before this tick's update: the held input is still u[k-1].
            y = x[0] + (0.5 * u_held(0));
            u_held = ctrl.control(r, ColVec<1>{{y}});
            x = sys.A * x + sys.B * (u_held + ColVec<1>{{0.4}});
        }
        CHECK(y == doctest::Approx(1.0).epsilon(1e-3));
    }

    TEST_CASE("state_mpc rejects any nonzero feedthrough, whatever its units") {
        StateSpace<2, 1, 1> sys{
            .A = {{1.0, 0.1}, {0.0, 1.0}},
            .B = {{0.005}, {0.1}},
            .C = {{1.0, 0.0}},
            .D = {{1e-13}}, // e.g. y in V, u in MW: real feedthrough, tiny number
            .Ts = 0.1,
        };
        CHECK_FALSE(design::state_mpc<5, 2>(sys).success);
        sys.D = Matrix<1, 1>{};
        CHECK(design::state_mpc<5, 2>(sys).success);
    }

    TEST_CASE("suggest_mpc_horizon does not read a defective integrator as a slow pole") {
        // Double integrator in a non-orthogonal basis with a large off-diagonal:
        // the computed eigenvalues sit O(√eps) off 1.
        const Matrix<2, 2> J{{1.0, 100.0}, {0.0, 1.0}};
        const Matrix<2, 2> Tm{{1.0, 0.3}, {0.7, 1.9}};
        const auto         Ti = Tm.inverse();
        REQUIRE(Ti.has_value());
        const StateSpace<2, 1, 1> sys{.A = Tm * J * *Ti, .B = {{0.0}, {1.0}}, .C = {{1.0, 0.0}}, .Ts = 0.1};
        CHECK_FALSE(design::suggest_mpc_horizon(sys).success);

        // A genuinely stable pole still yields a suggestion.
        const StateSpace<1, 1, 1> stable{.A = {{0.9}}, .B = {{1.0}}, .C = {{1.0}}, .Ts = 0.1};
        CHECK(design::suggest_mpc_horizon(stable).success);
        CHECK_FALSE(design::suggest_mpc_horizon(1e300, 1e-300).success); // NP not representable
    }

    TEST_CASE("discrete_lqr_from_continuous fails on a non-positive sample time") {
        const Matrix<2, 2> A{{0.0, 1.0}, {0.0, 0.0}};
        const Matrix<2, 1> B{{0.0}, {1.0}};
        CHECK_FALSE(design::discrete_lqr_from_continuous(A, B, Matrix<2, 2>::identity(), Matrix<1, 1>{{1.0}}, 0.0).success);
        CHECK(design::discrete_lqr_from_continuous(A, B, Matrix<2, 2>::identity(), Matrix<1, 1>{{1.0}}, 0.01).success);
    }

    TEST_CASE("ESC: a flagged-bad objective does not leak into the estimate") {
        const auto art = design::esc(0.1, 5.0, 2.0, 1e-3);
        REQUIRE(art.success);
        ExtremumSeekingController<double> a(art);
        ExtremumSeekingController<double> b(art);
        auto                              J = [](double u) { return -((u - 1.0) * (u - 1.0)); };
        double                            ua = a.input();
        double                            ub = b.input();
        for (int k = 0; k < 2000; ++k) {
            const bool fault = (k >= 500 && k < 510);
            ua = a.step(fault ? 1e6 : J(ua), !fault); // sensor spike, flagged
            ub = b.step(J(ub), !fault);               // same flag, sane value
        }
        CHECK(a.estimate() == doctest::Approx(b.estimate()).epsilon(1e-12));
    }

    TEST_CASE("repetitive design rejects rates whose period overflows size_t") {
        CHECK_FALSE(design::repetitive(1e300, 1e-300).success);
        CHECK(design::repetitive(10000.0, 50.0).config.period == 200);
        CHECK(design::repetitive_binomial<2>(10000.0f, 50.0f).success); // float taps still sum to ≤ 1
    }
}
