// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Canary: the embeddable umbrella must be self-contained and allocation-free.
//
// This TU includes ONLY "damp/control.hpp" (plus doctest). If it compiles, the
// umbrella is usable on its own for real controller code. The companion
// `make embedded-check` target additionally proves no <vector> is reachable
// from this header, which is the embedded (no-heap) guarantee.
#include "damp/backend.hpp"
#include "damp/control.hpp"
#include "damp/math/math.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

TEST_SUITE("Embedded Umbrella") {
    TEST_CASE("damp/control.hpp alone supports compile-time synthesis + runtime control") {
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
            return damp::design::lqi_bundle(sys, Q_aug, R);
        }();

        static_assert(artifacts.success, "synthesis must succeed at compile time via control.hpp alone");

        auto                   runtime = artifacts.runtime;
        const ColVec<2, float> x{1.0f, 0.0f};
        const ColVec<1, float> y{1.0f};
        const ColVec<1, float> r{0.0f};
        const auto             u = runtime.control(r, y, x);
        CHECK(u(0, 0) != doctest::Approx(0.0f));
    }
}
