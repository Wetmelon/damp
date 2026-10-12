// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <type_traits>

#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "doctest.h"

using namespace damp;

namespace {

/// Double integrator (position, speed), force input, position output.
constexpr StateSpace<2, 1, 1> double_integrator() {
    StateSpace<2, 1, 1> s{};
    s.A(0, 1) = 1.0;
    s.B(1, 0) = 1.0;
    s.C(0, 0) = 1.0;
    return s;
}

/// A sampled input sequence with some structure.
double u_at(int k) {
    return (k < 0) ? 0.0 : std::sin(0.37 * k) + ((k % 7 == 0) ? 1.0 : 0.0);
}

} // namespace

TEST_SUITE("Input delay") {
    TEST_CASE("input_delay_steps rounds a fractional delay up and keeps a whole one") {
        static_assert(input_delay_steps(0.0, 1e-3) == 0);
        static_assert(input_delay_steps(2e-3, 1e-3) == 2);
        static_assert(input_delay_steps(1.4e-3, 1e-3) == 2);
        static_assert(input_delay_steps(0.3e-3, 1e-3) == 1);
        static_assert(input_delay_steps(3e-3, 1e-3) == 3); // 3e-3/1e-3 is 2.9999999999999996 in double
    }

    TEST_CASE("a delayed discrete model outputs the undelayed output D samples later") {
        const auto sys = discretize(double_integrator(), 0.01, DiscretizationMethod::ZOH).value();
        const auto delayed = discretize(sys, 0.01, InputDelay<3>{}).value(); // three whole samples
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(delayed)>, StateSpace<5, 1, 1>>);
        ColVec<2> x{};
        ColVec<5> xa{};
        double    y_hist[200]{};
        double    worst = 0.0;
        for (int k = 0; k < 200; ++k) {
            y_hist[k] = (sys.C * x)(0, 0);
            const double ya = (delayed.C * xa)(0, 0);
            if (k >= 3) {
                worst = std::max(worst, std::abs(ya - y_hist[k - 3]));
            }
            x = sys.A * x + sys.B * ColVec<1>{{u_at(k)}};
            xa = delayed.A * xa + delayed.B * ColVec<1>{{u_at(k)}};
        }
        CHECK(worst < 1e-12);
        CHECK(discretize(sys, 0.01, InputDelay<0>{}).value().A == sys.A); // no delay: the same model
    }

    TEST_CASE("a whole-sample delay is the same from the continuous or the discrete model") {
        constexpr double Ts = 0.01;
        const auto       a = discretize(double_integrator(), Ts, InputDelay<2>{2.0 * Ts});
        REQUIRE(a.has_value());
        const auto b = discretize(discretize(double_integrator(), Ts, DiscretizationMethod::ZOH).value(), Ts, InputDelay<2>{2.0 * Ts}).value();
        CHECK((a->A - b.A).norm() < 1e-12);
        CHECK((a->B - b.B).norm() < 1e-12);
        CHECK((a->C - b.C).norm() < 1e-12);
    }

    TEST_CASE("a fractional delay is exact against the continuous plant") {
        // 1.4 samples of delay: during each sample the plant holds u(k−2) for 0.4·Ts, then u(k−1).
        constexpr double Ts = 0.01;
        constexpr double delay = 1.4 * Ts;
        const auto       model = discretize(double_integrator(), Ts, InputDelay<input_delay_steps(delay, Ts)>{delay});
        REQUIRE(model.has_value());
        ColVec<4> xa{};
        double    p = 0.0;
        double    v = 0.0;
        double    worst = 0.0;
        // Exact double-integrator step over h at constant acceleration a.
        const auto hold = [&](double a, double h) {
            p += (v * h) + (0.5 * a * h * h);
            v += a * h;
        };
        for (int k = 0; k < 300; ++k) {
            worst = std::max({worst, std::abs(xa[0] - p), std::abs(xa[1] - v)});
            hold(u_at(k - 2), 0.4 * Ts);
            hold(u_at(k - 1), 0.6 * Ts);
            xa = model->A * xa + model->B * ColVec<1>{{u_at(k)}};
        }
        CHECK(worst < 1e-12 * 1e3); // states grow to ~10³ over the run
    }

    TEST_CASE("the delay states are known: no process noise reaches them") {
        const auto sys = discretize(double_integrator(), 0.01, DiscretizationMethod::ZOH).value();
        const auto delayed = discretize(StateSpace<2, 1, 1, double, 2, 1>{.A = sys.A, .B = sys.B, .C = sys.C, .G = Matrix<2, 2>::identity(), .Ts = 0.01}, 0.01, InputDelay<2>{0.02}).value();
        for (size_t r = 2; r < 4; ++r) {
            CHECK(delayed.G(r, 0) == 0.0);
            CHECK(delayed.G(r, 1) == 0.0);
        }
        CHECK(delayed.G(0, 0) == 1.0);
    }

    TEST_CASE("a delay that cannot be honoured is refused, not silently wrong") {
        CHECK_FALSE(discretize(double_integrator(), 0.01, InputDelay<1>{0.014}).has_value()); // needs 2
        CHECK_FALSE(discretize(double_integrator(), 0.01, InputDelay<2>{-0.001}).has_value());
        CHECK_FALSE(discretize(double_integrator(), 0.01, InputDelay<2>{0.014}, DiscretizationMethod::Tustin).has_value()); // fractional needs ZOH
        const auto sys_d = discretize(double_integrator(), 0.01, DiscretizationMethod::ZOH).value();
        CHECK_FALSE(discretize(sys_d, 0.01, InputDelay<2>{0.014}).has_value());                                      // fractional on a discrete model
        CHECK(discretize(double_integrator(), 0.01, InputDelay<2>{0.02}, DiscretizationMethod::Tustin).has_value()); // whole: any method
    }
}
