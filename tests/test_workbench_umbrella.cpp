// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Canary: the host workbench umbrella must be self-contained for design +
// analysis + simulation on a workstation (heap OK). Complements
// test_embedded_umbrella.cpp for damp/control.hpp.
#include "damp/backend.hpp"
#include "damp/workbench.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

TEST_SUITE("Workbench Umbrella") {
    TEST_CASE("damp/workbench.hpp alone supports FRF helpers + design + sim types") {
        // Analysis sweeps allocate on the host — that is intentional here.
        const auto omega = analysis::logspace(1.0, 100.0, 8);
        CHECK(omega.size() == 8);
        CHECK(omega.front() == doctest::Approx(1.0));
        CHECK(omega.back() == doctest::Approx(100.0));

        StateSpace sys{
            .A = Matrix<1, 1>{{-1.0}},
            .B = Matrix<1, 1>{{1.0}},
            .C = Matrix<1, 1>{{1.0}},
        };
        const auto disc = *discretize(sys, 0.01, DiscretizationMethod::ZOH);
        CHECK(disc.Ts == doctest::Approx(0.01));

        // MATLAB®-style alias reachable through workbench
        const auto K = matlab::dlqr(
            disc.A, disc.B, Matrix<1, 1>{{1.0}}, Matrix<1, 1>{{1.0}}
        );
        CHECK(K.success);
        CHECK(K.K(0, 0) != 0.0);
    }
}
