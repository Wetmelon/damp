// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file test_multirate.cpp
 * @brief Multi-rate simulation harness (Phase 1 / roadmap #18)
 *
 * No plotlypp — CHECKs only.
 */

#include <cmath>
#include <cstddef>

#include "damp/backend.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/design/pid_design.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/simulation/integrator.hpp"
#include "damp/simulation/multirate.hpp"
#include "damp/simulation/simulate.hpp"
#include "damp/simulation/solver.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using namespace damp::sim;

TEST_SUITE("Multi-rate harness") {

    TEST_CASE("fires_at: tick 0 and integer divisors") {
        CHECK(fires_at(0, 1));
        CHECK(fires_at(0, 10));
        CHECK(fires_at(10, 10));
        CHECK_FALSE(fires_at(1, 10));
        CHECK_FALSE(fires_at(9, 10));
        CHECK_FALSE(fires_at(0, 0));
        CHECK_FALSE(fires_at(5, 0));
    }

    TEST_CASE("single-rate simulate_multirate matches simulate_sampled") {
        // Plant: ẋ = -x + u, y = x
        auto plant = [](double, const ColVec<1>& x, const ColVec<1>& u) -> ColVec<1> {
            return ColVec<1>{-x[0] + u[0]};
        };
        auto output = [](const ColVec<1>& x) -> ColVec<1> { return x; };

        constexpr double Ts = 0.01;
        constexpr double r = 1.0;
        auto             ctrl = [&](double, const ColVec<1>& y) -> ColVec<1> {
            // Proportional: u = r - y
            return ColVec<1>{r - y[0]};
        };
        auto step = [&](std::size_t, double t, const ColVec<1>& y) -> ColVec<1> {
            return ctrl(t, y);
        };

        RK4<1>          rk4;
        FixedStepSolver fine(rk4, 0.001);
        ColVec<1>       x0{0.0};

        auto a = simulate_sampled<1, 1, 1, double>(plant, output, ctrl, fine, Ts, x0, {0.0, 0.5});
        auto b = simulate_multirate<1, 1, 1, double>(plant, output, step, fine, Ts, x0, {0.0, 0.5});

        REQUIRE(a.t.size() == b.t.size());
        REQUIRE(a.t.size() > 10);
        for (std::size_t i = 0; i < a.t.size(); ++i) {
            CHECK(a.t[i] == doctest::Approx(b.t[i]).epsilon(1e-14));
            CHECK(a.x[i][0] == doctest::Approx(b.x[i][0]).epsilon(1e-12));
            CHECK(a.u[i][0] == doctest::Approx(b.u[i][0]).epsilon(1e-12));
        }
        // ẋ = r − 2x → x(t) = (r/2)(1 − e^{−2t}); t = 0.5 → ≈ 0.316
        CHECK(b.x.back()[0] == doctest::Approx(0.5 * (1.0 - std::exp(-1.0))).epsilon(1e-3));
    }

    TEST_CASE("two-rate: outer reference held between outer firings") {
        // Integrator ẋ = u. Outer posts a constant r = 1 every outer period;
        // inner proportional u = r − y every base tick. r must stay 1 between
        // outer fires; plant tracks to 1 under u = 1 − y (τ = 1 s).
        auto plant = [](double, const ColVec<1>&, const ColVec<1>& u) -> ColVec<1> { return u; };
        auto output = [](const ColVec<1>& x) -> ColVec<1> { return x; };

        constexpr std::size_t outer_div = 5;
        constexpr double      Ts_base = 0.01;
        std::size_t           outer_fires = 0;

        auto outer = [&](double, const ColVec<1>&) -> ColVec<1> {
            ++outer_fires;
            return ColVec<1>{1.0};
        };
        auto inner = [](double, const ColVec<1>& r, const ColVec<1>& y) -> ColVec<1> {
            return ColVec<1>{r[0] - y[0]};
        };

        RK4<1>          rk4;
        FixedStepSolver fine(rk4, 0.001);
        ColVec<1>       x0{0.0};

        auto tr = simulate_two_rate<1, 1, 1, 1, double>(
            plant, output, outer, outer_div, inner, 1, fine, Ts_base, x0, {0.0, 3.0}
        );

        REQUIRE_FALSE(tr.t.empty());
        REQUIRE(tr.r.size() == tr.t.size());
        REQUIRE(tr.base_tick_count == 300); // 3.0 / 0.01
        CHECK(outer_fires == 60);           // ticks 0,5,...,295

        for (const auto& rv : tr.r) {
            CHECK(rv[0] == doctest::Approx(1.0));
        }
        CHECK(tr.x.back()[0] == doctest::Approx(1.0).epsilon(0.05));
    }

    TEST_CASE("two-rate: outer holds last sample while plant evolves") {
        // Free response ẋ = −x (inner always applies u = 0). Outer latches y into
        // r only on its firings — between firings r is ZOH while y decays.
        auto plant = [](double, const ColVec<1>& x, const ColVec<1>&) -> ColVec<1> {
            return ColVec<1>{-x[0]};
        };
        auto output = [](const ColVec<1>& x) -> ColVec<1> { return x; };

        constexpr std::size_t outer_div = 4;
        constexpr double      Ts_base = 0.05;

        auto outer = [](double, const ColVec<1>& y) -> ColVec<1> { return y; };
        auto inner = [](double, const ColVec<1>&, const ColVec<1>&) -> ColVec<1> {
            return ColVec<1>{0.0};
        };

        RK4<1>          rk4;
        FixedStepSolver fine(rk4, 0.005);
        ColVec<1>       x0{1.0};

        MultiRateConfig<double> cfg{.Ts_base = Ts_base, .log_stride = 1, .max_base_ticks = 12};
        auto                    tr = simulate_two_rate<1, 1, 1, 1, double>(
            plant, output, outer, outer_div, inner, 1, fine, x0, {0.0, 10.0}, cfg
        );

        REQUIRE(tr.base_tick_count == 12);
        REQUIRE(tr.truncated);

        const double r0 = tr.r.front()[0];
        CHECK(r0 == doctest::Approx(1.0).epsilon(1e-12));

        // Mid first outer period [0, 4·Ts): y has decayed, r still r0
        bool found_mid = false;
        for (std::size_t i = 1; i < tr.t.size(); ++i) {
            if (tr.t[i] > 0.0 && tr.t[i] < 3.5 * Ts_base) {
                if (std::abs(tr.y[i][0] - r0) > 1e-3) {
                    CHECK(tr.r[i][0] == doctest::Approx(r0).epsilon(1e-12));
                    found_mid = true;
                    break;
                }
            }
        }
        CHECK(found_mid);

        // After second outer fire (tick 4), r steps from 1 → y(4·Ts) ≈ e^{−0.2}.
        // Fine samples that end the previous base period still carry the old r;
        // find the first sample where r has left r0.
        bool found_second = false;
        for (std::size_t i = 0; i < tr.t.size(); ++i) {
            if (std::abs(tr.r[i][0] - r0) > 1e-9) {
                const double r1 = tr.r[i][0];
                CHECK(r1 == doctest::Approx(std::exp(-4.0 * Ts_base)).epsilon(1e-2));
                // Held through the rest of that outer period while y keeps decaying
                for (std::size_t j = i + 1; j < tr.t.size() && tr.t[j] < 7.5 * Ts_base; ++j) {
                    CHECK(tr.r[j][0] == doctest::Approx(r1).epsilon(1e-12));
                    if (std::abs(tr.y[j][0] - r1) > 1e-3) {
                        found_second = true;
                    }
                }
                break;
            }
        }
        CHECK(found_second);
    }

    TEST_CASE("two-rate: slow outer step ref, fast inner PI tracks") {
        // Outer posts a step reference r = 1 only on its grid (open-loop outer).
        // Inner PI at base rate tracks that held r on G = 1/(s+1).
        constexpr double      Ts_base = 0.001; // 1 kHz inner
        constexpr std::size_t outer_div = 10;  // 100 Hz outer
        constexpr double      r_cmd = 1.0;

        StateSpace sys{
            .A = Matrix<1, 1>{{-1.0}},
            .B = Matrix<1, 1>{{1.0}},
            .C = Matrix<1, 1>{{1.0}},
            .Ts = 0.0,
        };

        const auto           inner_d = design::pi_pole_placement_first_order(1.0, 1.0, 20.0).discretize(Ts_base);
        PIController<double> inner_pi{inner_d};

        auto outer = [&](double, const ColVec<1>&) -> ColVec<1> { return ColVec<1>{r_cmd}; };
        auto inner = [&](double, const ColVec<1>& r_i, const ColVec<1>& y) -> ColVec<1> {
            return ColVec<1>{inner_pi.control(r_i[0], y[0])};
        };

        RK4<1>          rk4;
        FixedStepSolver fine(rk4, 0.0002);
        ColVec<1>       x0{0.0};

        auto tr = simulate_two_rate_lti<1, 1, 1, 1>(
            sys, outer, outer_div, inner, 1, fine, Ts_base, x0, {0.0, 1.5}
        );

        REQUIRE(tr.x.size() > 100);
        CHECK(tr.r.size() == tr.t.size());
        for (const auto& rv : tr.r) {
            CHECK(rv[0] == doctest::Approx(r_cmd));
        }
        CHECK(tr.x.back()[0] == doctest::Approx(r_cmd).epsilon(0.05));
    }

    TEST_CASE("two-rate cascade: outer position PI + inner velocity PI on double integrator") {
        // Plant: x = [pos, vel],  ṗ = v, v̇ = u. Outer (slow) sets v* from position
        // error; inner (base rate) sets u from velocity error — classic cascade rates.
        constexpr double      Ts_base = 0.001;
        constexpr std::size_t outer_div = 10;
        constexpr double      p_ref = 1.0;

        auto plant = [](double, const ColVec<2>& x, const ColVec<1>& u) -> ColVec<2> {
            return ColVec<2>{x[1], u[0]};
        };
        auto output = [](const ColVec<2>& x) -> ColVec<2> { return x; }; // y = [pos, vel]

        // Treat velocity plant as 1/s (a1=1, a0=0) and position under tight velocity as 1/s
        const auto           outer_d = design::pi_pole_placement_first_order(1.0, 0.0, 4.0).discretize(Ts_base * static_cast<double>(outer_div));
        const auto           inner_d = design::pi_pole_placement_first_order(1.0, 0.0, 30.0).discretize(Ts_base);
        PIController<double> outer_pi{outer_d};
        PIController<double> inner_pi{inner_d};

        auto outer = [&](double, const ColVec<2>& y) -> ColVec<1> {
            return ColVec<1>{outer_pi.control(p_ref, y[0])}; // v* from position
        };
        auto inner = [&](double, const ColVec<1>& r_v, const ColVec<2>& y) -> ColVec<1> {
            return ColVec<1>{inner_pi.control(r_v[0], y[1])}; // u from velocity
        };

        RK4<2>          rk4;
        FixedStepSolver fine(rk4, 0.0002);
        ColVec<2>       x0{0.0, 0.0};

        auto tr = simulate_two_rate<2, 1, 2, 1, double>(
            plant, output, outer, outer_div, inner, 1, fine, Ts_base, x0, {0.0, 3.0}
        );

        REQUIRE(tr.x.size() > 100);
        CHECK(tr.x.back()[0] == doctest::Approx(p_ref).epsilon(0.08));
        CHECK(std::abs(tr.x.back()[1]) < 0.15); // velocity near rest
    }

    TEST_CASE("max_base_ticks truncates") {
        auto plant = [](double, const ColVec<1>& x, const ColVec<1>&) -> ColVec<1> { return -x; };
        auto output = [](const ColVec<1>& x) -> ColVec<1> { return x; };
        auto step = [](std::size_t, double, const ColVec<1>&) -> ColVec<1> { return ColVec<1>{0.0}; };

        RK4<1>                  rk4;
        FixedStepSolver         fine(rk4, 0.01);
        MultiRateConfig<double> cfg{.Ts_base = 0.1, .max_base_ticks = 3};
        auto                    tr = simulate_multirate<1, 1, 1, double>(
            plant, output, step, fine, ColVec<1>{1.0}, {0.0, 10.0}, cfg
        );
        // 3 base periods of 0.1 → last time ≈ 0.3
        CHECK(tr.t.back() == doctest::Approx(0.3).epsilon(1e-9));
    }
}
