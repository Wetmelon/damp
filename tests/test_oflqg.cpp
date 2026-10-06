// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "damp/controllers/offset_free_lqg.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {
constexpr StateSpace<2, 1, 1, double, 2, 1> make_plant() {
    return StateSpace<2, 1, 1, double, 2, 1>{
        Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
        Matrix<2, 1>{{0.005}, {0.1}},
        Matrix<1, 2>{{1.0, 0.0}},
        Matrix<1, 1>::zeros(),
        Matrix<2, 2>::identity(),
        Matrix<1, 1>::identity(),
        0.1
    };
}
} // namespace

TEST_SUITE("OFLQG") {
    TEST_CASE("design succeeds") {
        constexpr auto sys = make_plant();
        constexpr auto result = design::discrete_oflqg(
            sys, Matrix<2, 2>::identity(), Matrix<1, 1>{{0.1}},
            Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, Matrix<1, 1>{{0.1}}
        );
        static_assert(result.success);
        CHECK(result.success);
        CHECK(result.lqr.K(0, 0) != 0.0);
        CHECK(result.kalman.L(0, 0) != 0.0);
    }

    TEST_CASE("tracks a setpoint under a constant input load") {
        const auto sys = make_plant();
        const auto result = design::discrete_oflqg(
            sys, Matrix<2, 2>::identity(), Matrix<1, 1>{{0.1}},
            Matrix<2, 2>{{0.05, 0.0}, {0.0, 0.05}}, Matrix<1, 1>{{0.05}},
            Matrix<1, 1>{{1.0}}
        );
        REQUIRE(result.success);

        OFLQG<2, 1, 1, double> ctrl{result};
        ColVec<2>              x{0.0, 0.0};
        const ColVec<1>        r{{1.0}};
        const double           w_true = 0.4;
        for (int k = 0; k < 400; ++k) {
            const ColVec<1> y{{x[0]}};
            const ColVec<1> u = ctrl.control(r, y);
            x = sys.A * x + sys.B * (u + ColVec<1>{{w_true}});
        }
        CHECK(x[0] == doctest::Approx(1.0).epsilon(0.05));
        CHECK(ctrl.disturbance_estimate()(0) == doctest::Approx(w_true).epsilon(0.15));
    }

    TEST_CASE("as<float> preserves success") {
        constexpr auto result = design::discrete_oflqg(
            make_plant(), Matrix<2, 2>::identity(), Matrix<1, 1>{{0.1}},
            Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, Matrix<1, 1>{{0.1}}
        );
        constexpr auto f = result.as<float>();
        static_assert(f.success);
        CHECK(f.success);
    }
}
