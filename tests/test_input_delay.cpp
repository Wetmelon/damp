// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

#include "damp/math/math.hpp"
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
    TEST_CASE("the register model is exact against the continuous plant for a delay of up to one sample") {
        // θ = 0.4·Ts: during each sample the plant holds u(k−1) for 0.4·Ts, then u(k); the register
        // holds u(k−1). θ = Ts: it holds u(k−1) for the whole sample.
        constexpr double Ts = 0.01;
        for (const double frac : {0.4, 1.0}) {
            CAPTURE(frac);
            const auto model = discretize(double_integrator(), Ts, InputDelay{frac * Ts});
            REQUIRE(model.has_value());
            static_assert(std::is_same_v<std::remove_cvref_t<decltype(*model)>, StateSpace<3, 1, 1>>);
            ColVec<3> xa{};
            double    p = 0.0;
            double    v = 0.0;
            double    worst = 0.0;
            // Exact double-integrator step over h at constant acceleration a.
            const auto hold = [&](double a, double h) {
                p += (v * h) + (0.5 * a * h * h);
                v += a * h;
            };
            for (int k = 0; k < 300; ++k) {
                worst = std::max({worst, std::abs(xa[0] - p), std::abs(xa[1] - v), std::abs(xa[2] - u_at(k - 1))});
                hold(u_at(k - 1), frac * Ts);
                hold(u_at(k), (1.0 - frac) * Ts);
                xa = model->A * xa + model->B * ColVec<1>{{u_at(k)}};
            }
            CHECK(worst < 1e-12 * 1e3); // states grow to ~10³ over the run
        }
    }

    TEST_CASE("boundaries 0, 16εTs, Ts−16εTs, Ts: snap, ordinary ZOH on the plant, Γ continuous") {
        constexpr double Ts = 0.01;
        const auto       plain = discretize(double_integrator(), Ts, DiscretizationMethod::ZOH).value();
        const double     snap = 16.0 * eps<double>() * Ts;
        const auto       delayed = [&](double theta) {
            const auto m = discretize(double_integrator(), Ts, InputDelay{theta});
            REQUIRE(m.has_value());
            return *m;
        };
        const auto gamma0 = [](const StateSpace<3, 1, 1>& m) { return Matrix<2, 1>(m.B.block<2, 1>(0, 0)); };
        const auto gamma1 = [](const StateSpace<3, 1, 1>& m) { return Matrix<2, 1>(m.A.block<2, 1>(0, 2)); };

        const auto zero = delayed(0.0);
        CHECK((zero.A.block<2, 2>(0, 0) - plain.A).norm() < 1e-15);
        CHECK((gamma0(zero) - plain.B).norm() < 1e-15);
        CHECK(gamma1(zero).norm() < 1e-15);
        CHECK((delayed(0.5 * snap).A - zero.A).norm() == 0.0); // snaps to no delay

        const auto whole = delayed(Ts);
        CHECK(gamma0(whole).norm() < 1e-15);
        CHECK((gamma1(whole) - plain.B).norm() < 1e-15);
        CHECK((delayed(Ts - 0.5 * snap).A - whole.A).norm() == 0.0); // snaps to one sample

        const double gscale = 1.0 + plain.B.norm();
        const auto   tiny = delayed(2.0 * snap); // just outside the zero snap
        CHECK((gamma0(tiny) - plain.B).norm() < 1e-9 * gscale);
        CHECK(gamma1(tiny).norm() < 1e-9 * gscale);
        const auto almost = delayed(Ts - 2.0 * snap);
        CHECK(gamma0(almost).norm() < 1e-9 * gscale);
        CHECK((gamma1(almost) - plain.B).norm() < 1e-9 * gscale);

        // Interior delays stay exact against the held continuous plant, including the ends.
        for (const double frac : {0.0, 2.0 * snap / Ts, 0.5, 1.0 - 2.0 * snap / Ts, 1.0}) {
            CAPTURE(frac);
            const auto model = delayed(frac * Ts);
            ColVec<3>  xa{};
            double     p = 0.0;
            double     v = 0.0;
            double     worst = 0.0;
            const auto hold = [&](double a, double h) {
                p += (v * h) + (0.5 * a * h * h);
                v += a * h;
            };
            for (int k = 0; k < 300; ++k) {
                worst = std::max({worst, std::abs(xa[0] - p), std::abs(xa[1] - v), std::abs(xa[2] - u_at(k - 1))});
                hold(u_at(k - 1), frac * Ts);
                hold(u_at(k), (1.0 - frac) * Ts);
                xa = model.A * xa + model.B * ColVec<1>{{u_at(k)}};
            }
            CHECK(worst < 1e-12 * 1e3);
        }
    }

    TEST_CASE("no delay leaves the register empty; a whole sample is Γ₀ = 0, Γ₁ = Γ; continuous and discrete models agree") {
        constexpr double Ts = 0.01;
        const auto       plain = discretize(double_integrator(), Ts, DiscretizationMethod::ZOH).value();
        for (const double theta : {0.0, Ts}) {
            CAPTURE(theta);
            const auto a = discretize(double_integrator(), Ts, InputDelay{theta});
            const auto b = discretize(plain, Ts, InputDelay{theta});
            REQUIRE(a.has_value());
            REQUIRE(b.has_value());
            CHECK((a->A - b->A).norm() < 1e-12);
            CHECK((a->B - b->B).norm() < 1e-12);
            CHECK((a->C - b->C).norm() < 1e-12);
            CHECK((a->A.block<2, 2>(0, 0) - plain.A).norm() < 1e-15);
            const Matrix<2, 1> gamma0 = a->B.block<2, 1>(0, 0);
            const Matrix<2, 1> gamma1 = a->A.block<2, 1>(0, 2);
            const Matrix<2, 1> expect0 = (theta == 0.0) ? plain.B : Matrix<2, 1>{};
            const Matrix<2, 1> expect1 = (theta == 0.0) ? Matrix<2, 1>{} : plain.B;
            CHECK((gamma0 - expect0).norm() < 1e-15);
            CHECK((gamma1 - expect1).norm() < 1e-15);
            CHECK(a->B(2, 0) == 1.0); // the register loads u(k)
            CHECK(a->A(2, 2) == 0.0);
        }
        // A whole sample with another method: the register feeds that method's Γ.
        const auto tustin = discretize(double_integrator(), Ts, DiscretizationMethod::Tustin).value();
        const auto late = discretize(double_integrator(), Ts, InputDelay{Ts}, DiscretizationMethod::Tustin);
        REQUIRE(late.has_value());
        CHECK((late->A.block<2, 1>(0, 2) - tustin.B).norm() < 1e-15);
        CHECK(Matrix<2, 1>(late->B.block<2, 1>(0, 0)).norm() == 0.0);
    }

    TEST_CASE("the feedthrough reads the input the plant holds at the sample instant") {
        // y = x + 0.5·u on a first-order lag: with no delay the plant holds u(k) at the sample instant,
        // with any delay in (0, Ts] it still holds u(k−1), the register.
        StateSpace<1, 1, 1> lag{};
        lag.A(0, 0) = -1.0;
        lag.B(0, 0) = 1.0;
        lag.C(0, 0) = 1.0;
        lag.D(0, 0) = 0.5;
        constexpr double Ts = 0.1;
        const auto       now = discretize(lag, Ts, InputDelay{0.0}).value();
        CHECK(now.D(0, 0) == 0.5);
        CHECK(now.C(0, 1) == 0.0);
        for (const double theta : {0.4 * Ts, Ts}) {
            CAPTURE(theta);
            const auto held = discretize(lag, Ts, InputDelay{theta}).value();
            CHECK(held.D(0, 0) == 0.0);
            CHECK(held.C(0, 1) == 0.5);
        }
    }

    TEST_CASE("the register is a command already sent: no process noise reaches it") {
        const StateSpace<2, 1, 1, double, 2, 1> plant{.A = double_integrator().A, .B = double_integrator().B, .C = double_integrator().C, .G = Matrix<2, 2>::identity()};
        for (const double theta : {0.004, 0.01}) {
            CAPTURE(theta);
            const auto delayed = discretize(plant, 0.01, InputDelay{theta}).value();
            CHECK(delayed.G(2, 0) == 0.0);
            CHECK(delayed.G(2, 1) == 0.0);
            CHECK(delayed.G(0, 0) != 0.0);
        }
    }

    TEST_CASE("the register gets no process noise when NW == NX + NU (G is not the square default identity)") {
        const StateSpace<2, 1, 1, double, 3> plant{.A = double_integrator().A, .B = double_integrator().B, .C = double_integrator().C, .G = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}}};
        for (const double theta : {0.0, 0.004, 0.01}) {
            CAPTURE(theta);
            const auto delayed = discretize(plant, 0.01, InputDelay{theta}).value();
            CHECK(delayed.G(2, 0) == 0.0);
            CHECK(delayed.G(2, 1) == 0.0);
            CHECK(delayed.G(2, 2) == 0.0);
            CHECK(delayed.G(0, 0) != 0.0);
        }
    }

    TEST_CASE("a delay that cannot be honoured is refused, not silently wrong") {
        constexpr double nan = std::numeric_limits<double>::quiet_NaN();
        CHECK_FALSE(discretize(double_integrator(), 0.01, InputDelay{0.014}).has_value()); // more than one sample
        CHECK_FALSE(discretize(double_integrator(), 0.01, InputDelay{0.02}).has_value());  // two whole samples
        CHECK_FALSE(discretize(double_integrator(), 0.01, InputDelay{-0.001}).has_value());
        CHECK_FALSE(discretize(double_integrator(), 0.01, InputDelay{nan}).has_value());
        CHECK_FALSE(discretize(double_integrator(), 0.01, InputDelay{0.004}, DiscretizationMethod::Tustin).has_value()); // a fraction needs ZOH
        const auto sys_d = discretize(double_integrator(), 0.01, DiscretizationMethod::ZOH).value();
        CHECK_FALSE(discretize(sys_d, 0.01, InputDelay{0.004}).has_value()); // a fraction on a discrete model
        CHECK_FALSE(discretize(sys_d, 0.02, InputDelay{0.02}).has_value());  // a discrete model at another Ts
        CHECK(discretize(sys_d, 0.01, InputDelay{0.01}).has_value());        // a whole sample on a discrete model
        for (const auto method : {DiscretizationMethod::ForwardEuler, DiscretizationMethod::ZOH, DiscretizationMethod::Tustin}) {
            CHECK(discretize(double_integrator(), 0.01, InputDelay{0.01}, method).has_value()); // a whole sample: any method
        }
    }
}
