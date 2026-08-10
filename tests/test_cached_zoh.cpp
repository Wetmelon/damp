// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file test_cached_zoh.cpp
 * @brief Cached ZOH discrete maps (expm once, Discrete step many times)
 */

#include <cmath>
#include <vector>

#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/simulation/cached_zoh.hpp"
#include "damp/simulation/hybrid.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using namespace damp::sim;

TEST_SUITE("Cached ZOH") {

    TEST_CASE("CachedZoh step matches Exact for scalar decay") {
        const double a = 2.0;
        const double h = 0.1;
        StateSpace   sys{
              .A = Matrix<1, 1>{{-a}},
              .B = Matrix<1, 1>{{0.0}},
              .C = Matrix<1, 1>{{1.0}},
              .Ts = 0.0,
        };
        ColVec<1> x0{1.0};
        ColVec<1> u{0.0};

        const auto zoh = CachedZoh<1, 1>::from(sys, h);
        const auto x_cache = zoh.step(x0, u);
        const auto x_exact = exact_lti_step(sys, x0, u, h);

        CHECK(x_cache(0) == doctest::Approx(x_exact(0)).epsilon(1e-14));
        CHECK(x_cache(0) == doctest::Approx(std::exp(-a * h)).epsilon(1e-12));
    }

    TEST_CASE("CachedZoh integrator plant ẋ = u matches ramp") {
        StateSpace sys{
            .A = Matrix<1, 1>{{0.0}},
            .B = Matrix<1, 1>{{1.0}},
            .C = Matrix<1, 1>{{1.0}},
            .Ts = 0.0,
        };
        const double h = 0.05;
        ColVec<1>    x{2.0};
        ColVec<1>    u{3.0};
        const auto   zoh = CachedZoh<1, 1>::from(sys.A, sys.B, h);

        x = zoh.step(x, u);
        CHECK(x(0) == doctest::Approx(2.0 + 3.0 * h).epsilon(1e-12));
        x = zoh.step(x, u);
        CHECK(x(0) == doctest::Approx(2.0 + 2.0 * 3.0 * h).epsilon(1e-12));
    }

    TEST_CASE("simulate_cached_zoh held step input on first-order lag") {
        // ẋ = -x + u, u = 1, x0 = 0 → x(t) = 1 - e^{-t}
        StateSpace sys{
            .A = Matrix<1, 1>{{-1.0}},
            .B = Matrix<1, 1>{{1.0}},
            .C = Matrix<1, 1>{{1.0}},
            .Ts = 0.0,
        };
        const double h = 0.01;
        const size_t n = 100;
        ColVec<1>    x0{0.0};
        ColVec<1>    u{1.0};

        auto sim = simulate_cached_zoh(sys, x0, u, h, n);
        CHECK(sim.t.size() == n + 1);
        CHECK(sim.t.back() == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(sim.x.back()(0) == doctest::Approx(1.0 - std::exp(-1.0)).epsilon(1e-6));
        CHECK(sim.y.back()(0) == doctest::Approx(sim.x.back()(0)).epsilon(1e-14));
        CHECK(sim.u.back()(0) == doctest::Approx(1.0));
    }

    TEST_CASE("simulate_cached_zoh_multirate holds u across fine steps") {
        StateSpace sys{
            .A = Matrix<1, 1>{{0.0}},
            .B = Matrix<1, 1>{{1.0}},
            .C = Matrix<1, 1>{{1.0}},
            .Ts = 0.0,
        };
        // Two control holds: u=1 for 3 fine steps, then u=2 for 3 fine steps; h=0.1
        std::vector<ColVec<1>> u_ctrl{ColVec<1>{1.0}, ColVec<1>{2.0}};
        auto                   sim = simulate_cached_zoh_multirate(sys, ColVec<1>{0.0}, u_ctrl, 0.1, 3);

        // 1 + 2*3 samples
        CHECK(sim.t.size() == 7);
        CHECK(sim.x[3](0) == doctest::Approx(0.3).epsilon(1e-12)); // after 3×0.1 at u=1
        CHECK(sim.x.back()(0) == doctest::Approx(0.3 + 0.6).epsilon(1e-12));
        CHECK(sim.u[1](0) == doctest::Approx(1.0));
        CHECK(sim.u[4](0) == doctest::Approx(2.0));
    }

    TEST_CASE("as_discrete_ss Ts and one-step map match cache") {
        StateSpace cont{
            .A = Matrix<1, 1>{{-1.0}},
            .B = Matrix<1, 1>{{1.0}},
            .C = Matrix<1, 1>{{1.0}},
            .D = Matrix<1, 1>{{0.0}},
            .Ts = 0.0,
        };
        const double h = 0.02;
        const auto   zoh = CachedZoh<1, 1>::from(cont, h);
        const auto   dsys = zoh.as_discrete_ss(cont.C, cont.D);
        CHECK(dsys.Ts == doctest::Approx(h));
        ColVec<1> x{1.0};
        ColVec<1> u{0.5};
        CHECK(zoh.step(x, u)(0) == doctest::Approx((dsys.A * x + dsys.B * u)(0)).epsilon(1e-14));
    }
}
