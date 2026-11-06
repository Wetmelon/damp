// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "damp/analysis/norms.hpp"
#include "damp/analysis/sensitivity.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/math/complex.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "doctest.h"

using namespace damp;

namespace {
using Cplx = damp::complex<double>;

// P(s) = 1/(s + 1).
constexpr StateSpace<1, 1, 1> kLag{
    .A = Matrix<1, 1>{{-1.0}},
    .B = Matrix<1, 1>{{1.0}},
    .C = Matrix<1, 1>{{1.0}},
    .D = Matrix<1, 1>{{0.0}},
};

constexpr design::PIDGains<double> kPI{.Kp = 2.0, .Ki = 3.0};

// A coupled 2×2 plant and a static 2×2 controller.
constexpr StateSpace<2, 2, 2> kPlant2{
    .A = Matrix<2, 2>{{-1.0, 0.4}, {0.0, -2.0}},
    .B = Matrix<2, 2>{{1.0, 0.0}, {0.3, 1.0}},
    .C = Matrix<2, 2>{{1.0, 0.5}, {0.2, 1.0}},
    .D = Matrix<2, 2>{{0.0, 0.1}, {0.0, 0.0}},
};
constexpr StateSpace<0, 2, 2> kGain2{.D = Matrix<2, 2>{{2.0, 0.3}, {0.1, 1.5}}};

template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
Matrix<NY, NU, Cplx> frf(const StateSpace<NX, NU, NY, T, NW, NV>& sys, Cplx s) {
    const auto G = eval_frf(sys, s);
    REQUIRE(G.has_value());
    return *G;
}

void check_close(Cplx got, Cplx want, double rel = 1e-12) {
    CHECK(damp::abs(got - want) <= rel * (1.0 + damp::abs(want)));
}

template<size_t R, size_t C>
void check_close(const Matrix<R, C, Cplx>& got, const Matrix<R, C, Cplx>& want, double rel = 1e-12) {
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            check_close(got(i, j), want(i, j), rel);
        }
    }
}
} // namespace

TEST_SUITE("Gang of four") {
    TEST_CASE("SISO PI on a first-order lag matches S, T, PS and CS") {
        const auto g = analysis::gang_of_four(kLag, *kPI.to_ss());
        REQUIRE(g.success);
        CHECK(g.stable);
        for (const double w : {0.1, 1.0, 10.0}) {
            const Cplx s{0.0, w};
            const Cplx P = 1.0 / (s + 1.0);
            const Cplx C = 2.0 + (3.0 / s);
            const Cplx den = 1.0 + (P * C);
            check_close(frf(g.sensitivity, s)(0, 0), 1.0 / den);
            check_close(frf(g.complementary_sensitivity, s)(0, 0), (P * C) / den);
            check_close(frf(g.load_sensitivity, s)(0, 0), P / den);
            check_close(frf(g.noise_sensitivity, s)(0, 0), C / den);
        }
    }

    TEST_CASE("MIMO members satisfy S + T = I, PS = S·P and CS = C·S") {
        const auto g = analysis::gang_of_four(kPlant2, kGain2);
        REQUIRE(g.success);
        CHECK(g.stable);
        const Cplx s{0.0, 0.7};
        const auto P = frf(kPlant2, s);
        const auto C = frf(kGain2, s);
        const auto S = frf(g.sensitivity, s);
        const auto I = Matrix<2, 2, Cplx>::identity();
        check_close(S + frf(g.complementary_sensitivity, s), I);
        check_close(S * (I + (P * C)), I);
        check_close(frf(g.load_sensitivity, s), S * P);
        check_close(frf(g.noise_sensitivity, s), C * S);
    }

    TEST_CASE("Discrete loop: S = 1/(1 + P(z)C) on the unit circle") {
        const auto Pd = discretize(kLag, 0.1);
        REQUIRE(Pd.success);
        auto Cd = ss_gain(0.5);
        Cd.Ts = 0.1;
        const auto g = analysis::gang_of_four(*Pd, Cd);
        REQUIRE(g.success);
        CHECK(g.stable);
        const Cplx z{std::cos(0.8), std::sin(0.8)};
        const Cplx P = frf(*Pd, z)(0, 0);
        check_close(frf(g.sensitivity, z)(0, 0), 1.0 / (1.0 + (0.5 * P)));
        check_close(frf(g.complementary_sensitivity, z)(0, 0), (0.5 * P) / (1.0 + (0.5 * P)));
    }

    TEST_CASE("Peak sensitivity from norm_hinf bounds a dense grid from above") {
        // P = 1/(s + 1)², C = 20: closed-loop ζ = 1/√21, so |S| peaks above 1 at a finite frequency.
        constexpr StateSpace<2, 1, 1> lag2{
            .A = Matrix<2, 2>{{-1.0, 1.0}, {0.0, -1.0}},
            .B = Matrix<2, 1>{{0.0}, {1.0}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .D = Matrix<1, 1>{{0.0}},
        };
        const auto g = analysis::gang_of_four(lag2, ss_gain(20.0));
        const auto Ms = analysis::norm_hinf(g.sensitivity);
        REQUIRE(Ms.success);
        CHECK(Ms.gain > 1.1);
        double grid_max = 0.0;
        for (int k = 0; k <= 200000; ++k) {
            const double w = std::pow(10.0, -3.0 + (6.0 * k / 200000.0));
            grid_max = std::max(grid_max, damp::abs(frf(g.sensitivity, Cplx{0.0, w})(0, 0)));
        }
        CHECK(Ms.gain >= grid_max * (1.0 - 1e-12));
        CHECK(Ms.gain <= grid_max * (1.0 + 1e-6));
    }

    TEST_CASE("Unstable closed loop: success, not stable, and norm_hinf refuses") {
        constexpr StateSpace<1, 1, 1> unstable{
            .A = Matrix<1, 1>{{1.0}},
            .B = Matrix<1, 1>{{1.0}},
            .C = Matrix<1, 1>{{1.0}},
            .D = Matrix<1, 1>{{0.0}},
        };
        const auto g = analysis::gang_of_four(unstable, ss_gain(0.5)); // closed-loop pole at +0.5
        CHECK(g.success);
        CHECK_FALSE(g.stable);
        CHECK_FALSE(analysis::norm_hinf(g.sensitivity).success);
    }

    TEST_CASE("Refused: sample times differ, or the loop is not well posed") {
        auto Cd = ss_gain(0.5);
        Cd.Ts = 0.1;
        CHECK_FALSE(analysis::gang_of_four(kLag, Cd).success);

        const StateSpace<0, 1, 1> unit{.D = Matrix<1, 1>{{1.0}}};
        const auto                ill_posed = analysis::gang_of_four(unit, ss_gain(-1.0)); // 1 + PC = 0
        CHECK_FALSE(ill_posed.success);
        CHECK_FALSE(ill_posed.stable);
    }

    TEST_CASE("Evaluates at compile time") {
        static constexpr auto g = analysis::gang_of_four(kLag, ss_gain(2.0));
        static_assert(g.success && g.stable);
        static_assert(g.sensitivity.A(0, 0) == -3.0);
        CHECK(g.sensitivity.D(0, 0) == 1.0);
    }
}

TEST_SUITE("Gang of seven") {
    TEST_CASE("PIDGains::to_ss_2dof realizes the setpoint and error channels") {
        constexpr design::PIDGains<double> pid{.Kp = 2.0, .Ki = 3.0, .Kd = 0.5, .b = 0.4, .c = 0.0, .Tf = 0.05};
        const auto                         K = pid.to_ss_2dof();
        REQUIRE(K.has_value());
        for (const double w : {0.1, 1.0, 30.0}) {
            const auto Kw = frf(*K, Cplx{0.0, w});
            check_close(Kw(0, 0), pid.eval_ff(w));
            check_close(Kw(0, 1), -pid.eval_fb(w));
        }

        const auto Kpi = kPI.to_ss_2dof(); // K_d = 0, T_f = 0: the filter state is disconnected
        REQUIRE(Kpi.has_value());
        CHECK(Kpi->A(1, 1) < 0.0);
        check_close(frf(*Kpi, Cplx{0.0, 2.0})(0, 1), -kPI.eval_fb(2.0));

        constexpr design::PIDGains<double> unfiltered{.Kp = 1.0, .Kd = 0.5};
        CHECK_FALSE(unfiltered.to_ss_2dof().has_value());
    }

    TEST_CASE("Setpoint weights change the reference responses, not the four") {
        constexpr design::PIDGains<double> pid{.Kp = 2.0, .Ki = 3.0, .Kd = 0.5, .b = 0.4, .c = 0.0, .Tf = 0.05};
        const auto                         g = analysis::gang_of_seven(kLag, *pid.to_ss_2dof());
        REQUIRE(g.success);
        CHECK(g.stable);
        for (const double w : {0.1, 1.0, 10.0}) {
            const Cplx s{0.0, w};
            const Cplx P = 1.0 / (s + 1.0);
            const Cplx Cff = pid.eval_ff(w);
            const Cplx Cfb = pid.eval_fb(w);
            const Cplx den = 1.0 + (P * Cfb);
            check_close(frf(g.sensitivity, s)(0, 0), 1.0 / den);
            check_close(frf(g.complementary_sensitivity, s)(0, 0), (P * Cfb) / den);
            check_close(frf(g.load_sensitivity, s)(0, 0), P / den);
            check_close(frf(g.noise_sensitivity, s)(0, 0), Cfb / den);
            check_close(frf(g.reference_to_output, s)(0, 0), (P * Cff) / den);
            check_close(frf(g.reference_to_control, s)(0, 0), Cff / den);
            check_close(frf(g.reference_to_error, s)(0, 0), 1.0 - ((P * Cff) / den));
        }
    }

    TEST_CASE("One-degree-of-freedom controller: reference responses equal T, CS and S") {
        constexpr design::PIDGains<double> pid{.Kp = 2.0, .Ki = 3.0, .Kd = 0.5, .Tf = 0.05};
        const auto                         g = analysis::gang_of_seven(kLag, *pid.to_ss_2dof());
        REQUIRE(g.success);
        const Cplx s{0.0, 1.3};
        check_close(frf(g.reference_to_output, s)(0, 0), frf(g.complementary_sensitivity, s)(0, 0));
        check_close(frf(g.reference_to_control, s)(0, 0), frf(g.noise_sensitivity, s)(0, 0));
        check_close(frf(g.reference_to_error, s)(0, 0), frf(g.sensitivity, s)(0, 0));
    }
}
