// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file test_hybrid_lti.cpp
 * @brief Piecewise-LTI Exact hybrid simulation (Phase 1)
 *
 * No plotlypp — keep incremental test builds light. CHECKs only.
 */

#include <cmath>

#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/simulation/hybrid.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using namespace damp::sim;

TEST_SUITE("Hybrid piecewise LTI") {

    TEST_CASE("Exact LTI step matches ZOH for scalar integrator ẋ = u") {
        StateSpace sys{
            .A = Matrix<1, 1>{{0.0}},
            .B = Matrix<1, 1>{{1.0}},
            .C = Matrix<1, 1>{{1.0}},
            .Ts = 0.0,
        };
        ColVec<1>    x0{2.0};
        ColVec<1>    u{3.0};
        const double h = 0.05;
        const auto   x1 = exact_lti_step(sys, x0, u, h);
        CHECK(x1(0) == doctest::Approx(2.0 + 3.0 * h).epsilon(1e-12));
    }

    TEST_CASE("Exact LTI step: stable scalar ẋ = -a x") {
        const double a = 2.0;
        StateSpace   sys{
              .A = Matrix<1, 1>{{-a}},
              .B = Matrix<1, 1>{{0.0}},
              .C = Matrix<1, 1>{{1.0}},
              .Ts = 0.0,
        };
        ColVec<1>    x0{1.0};
        ColVec<1>    u{0.0};
        const double h = 0.1;
        const auto   x1 = exact_lti_step(sys, x0, u, h);
        CHECK(x1(0) == doctest::Approx(std::exp(-a * h)).epsilon(1e-10));
    }

    TEST_CASE("Two-mode RL chopper: open-loop PWM events only") {
        // Series RL: L=1e-3, R=1, Vin=10 ON / 0 OFF
        constexpr double L = 1e-3;
        constexpr double R = 1.0;
        constexpr double Vin = 10.0;
        constexpr double Tp = 100e-6; // 10 kHz
        constexpr double duty = 0.4;

        StateSpace on{
            .A = Matrix<1, 1>{{-R / L}},
            .B = Matrix<1, 1>{{1.0 / L}},
            .C = Matrix<1, 1>{{1.0}},
            .Ts = 0.0,
        };
        StateSpace off = on;

        PiecewiseLTI<2, 1, 1, 1> plant{.modes = {on, off}};

        // Explicit edge calendar: ON [kTp, kTp+d), OFF [kTp+d, (k+1)Tp)
        double next_edge = duty * Tp; // first falling (start in ON)
        bool   high = true;

        auto schedule = [&](double /*t*/, std::size_t, const ColVec<1>&, const ColVec<1>&) {
            return next_edge;
        };

        auto on_event = [&](double /*t*/, std::size_t, const ColVec<1>& x, const ColVec<1>&, const ColVec<1>&) {
            HybridEventAction<1, 1, 1> act{};
            act.x = x;
            if (high) {
                // falling edge → OFF
                high = false;
                act.mode = 1;
                act.u = ColVec<1>{0.0};
                // next rising at ceil to next period
                const double k = std::floor((next_edge + 1e-18) / Tp);
                next_edge = (k + 1.0) * Tp;
            } else {
                // rising edge → ON
                high = true;
                act.mode = 0;
                act.u = ColVec<1>{Vin};
                next_edge = next_edge + duty * Tp;
            }
            return act;
        };

        ColVec<1>        x0{0.0};
        ColVec<1>        u0{Vin};
        constexpr double t_end = 50 * Tp;

        auto tr = simulate_piecewise_lti(plant, 0, x0, u0, {0.0, t_end}, schedule, on_event);

        CHECK(tr.event_count >= 90); // ~2 edges per period × 50
        CHECK(tr.t.size() == tr.x.size());
        CHECK(tr.mode.size() == tr.t.size());
        CHECK_FALSE(tr.truncated);
        CHECK(tr.t.back() == doctest::Approx(t_end).epsilon(1e-9));

        double       i_sum = 0.0;
        std::size_t  n = 0;
        const double t_mid = 0.5 * t_end;
        for (std::size_t k = 0; k < tr.t.size(); ++k) {
            if (tr.t[k] >= t_mid) {
                i_sum += tr.x[k](0);
                ++n;
            }
        }
        REQUIRE(n > 0);
        const double i_avg = i_sum / static_cast<double>(n);
        const double i_ss = duty * Vin / R;
        CHECK(i_avg == doctest::Approx(i_ss).epsilon(0.15));
    }

    TEST_CASE("FixedStride denser than EventsOnly on same PWM") {
        constexpr double L = 1e-3;
        constexpr double R = 1.0;
        constexpr double Vin = 10.0;
        constexpr double Tp = 100e-6;
        constexpr double duty = 0.5;

        StateSpace on{
            .A = Matrix<1, 1>{{-R / L}},
            .B = Matrix<1, 1>{{1.0 / L}},
            .C = Matrix<1, 1>{{1.0}},
            .Ts = 0.0,
        };
        PiecewiseLTI<2, 1, 1, 1> plant{.modes = {on, on}};

        double next_edge = duty * Tp;
        bool   high = true;
        auto   schedule = [&](double, std::size_t, const ColVec<1>&, const ColVec<1>&) {
            return next_edge;
        };
        auto on_event = [&](double, std::size_t, const ColVec<1>& x, const ColVec<1>&, const ColVec<1>&) {
            HybridEventAction<1, 1, 1> act{};
            act.x = x;
            if (high) {
                high = false;
                act.mode = 1;
                act.u = ColVec<1>{0.0};
                const double k = std::floor((next_edge + 1e-18) / Tp);
                next_edge = (k + 1.0) * Tp;
            } else {
                high = true;
                act.mode = 0;
                act.u = ColVec<1>{Vin};
                next_edge = next_edge + duty * Tp;
            }
            return act;
        };

        EventSimConfig<> cfg_ev{.log = LogPolicy::EventsOnly};
        EventSimConfig<> cfg_dense{.log = LogPolicy::FixedStride, .h_fine = Tp / 20.0};

        ColVec<1>        x0{0.0};
        ColVec<1>        u0{Vin};
        constexpr double t_end = 5 * Tp;

        auto a = simulate_piecewise_lti(plant, 0, x0, u0, {0.0, t_end}, schedule, on_event, cfg_ev);

        next_edge = duty * Tp;
        high = true;
        auto b = simulate_piecewise_lti(plant, 0, x0, u0, {0.0, t_end}, schedule, on_event, cfg_dense);

        CHECK(b.t.size() > a.t.size());
        CHECK(a.x.back()(0) == doctest::Approx(b.x.back()(0)).epsilon(1e-9));
    }

    TEST_CASE("locate_zero_crossing_exact finds i=0 on freewheel segment") {
        // OFF: di/dt = -v/L with v fixed via algebraic trick: use 1-state  di/dt = -a
        // Simple: ẋ = -1, x(0)=0.1 → zero at t=0.1
        StateSpace sys{
            .A = Matrix<1, 1>{{0.0}},
            .B = Matrix<1, 1>{{-1.0}},
            .C = Matrix<1, 1>{{1.0}},
            .Ts = 0.0,
        };
        ColVec<1>    x0{0.1};
        ColVec<1>    u{1.0}; // ẋ = -u
        auto         g = [](const ColVec<1>& x) { return x(0); };
        const double t_zc = locate_zero_crossing_exact(sys, x0, u, 0.0, 1.0, g, 1e-14);
        CHECK(t_zc == doctest::Approx(0.1).epsilon(1e-9));
    }
}

// NPY export (host) — no plotly
#include <cstdio>

#include "damp/simulation/npy_export.hpp"

TEST_SUITE("NPY export") {
    TEST_CASE("write_npy_f64 round-trips header magic") {
        std::vector<double> data = {1.0, 2.0, 3.0, 4.0};
        const char*         path = "tests/build/test_trace.npy";
        REQUIRE(write_npy_f64(path, {2, 2}, data));
        FILE* f = std::fopen(path, "rb");
        REQUIRE(f != nullptr);
        char magic[6]{};
        REQUIRE(std::fread(magic, 1, 6, f) == 6);
        CHECK(magic[0] == '\x93');
        CHECK(magic[1] == 'N');
        CHECK(magic[2] == 'U');
        CHECK(magic[3] == 'M');
        CHECK(magic[4] == 'P');
        CHECK(magic[5] == 'Y');
        std::fclose(f);
    }
}
