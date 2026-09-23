// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>

#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

// Scalar type T sits before noise dims NW/NV so float (and noise) need no dummy 0,0:
//   StateSpace<2,1,1>              — double, no noise
//   StateSpace<2,1,1,float>        — float, no noise
//   StateSpace<2,1,1,double,2,1>   — double with G/H channels
static_assert(std::is_same_v<decltype(StateSpace<2, 1, 1>{}.Ts), double>);
static_assert(std::is_same_v<decltype(StateSpace<2, 1, 1, float>{}.Ts), float>);
static_assert(StateSpace<2, 1, 1, double, 2, 1>{}.G.cols() == 2);
static_assert(StateSpace<2, 1, 1, double, 2, 1>{}.H.cols() == 1);

TEST_CASE("StateSpace Series Connection") {
    // Simple first-order systems: x' = -a*x + u, y = x
    StateSpace<1, 1, 1, double, 1, 1> sys1{
        .A = {{-1.0}},
        .B = {{1.0}},
        .C = {{1.0}},
        .D = {{0.0}},
    };

    StateSpace<1, 1, 1, double, 1, 1> sys2{
        .A = {{-2.0}},
        .B = {{1.0}},
        .C = {{1.0}},
        .D = {{0.0}},
    };

    auto result = *series(sys1, sys2);

    // Expected dimensions
    CHECK(result.A.rows() == 2);
    CHECK(result.A.cols() == 2);
    CHECK(result.B.rows() == 2);
    CHECK(result.B.cols() == 1);
    CHECK(result.C.rows() == 1);
    CHECK(result.C.cols() == 2);
    CHECK(result.D.rows() == 1);
    CHECK(result.D.cols() == 1);

    // Check A matrix structure
    CHECK(result.A(0, 0) == doctest::Approx(-1.0));
    CHECK(result.A(1, 0) == doctest::Approx(1.0)); // B2*C1 = 1*1
    CHECK(result.A(1, 1) == doctest::Approx(-2.0));

    // Check B matrix
    CHECK(result.B(0, 0) == doctest::Approx(1.0));
    CHECK(result.B(1, 0) == doctest::Approx(0.0)); // B2*D1 = 1*0

    // Check C matrix
    CHECK(result.C(0, 0) == doctest::Approx(0.0)); // D2*C1 = 0*1
    CHECK(result.C(0, 1) == doctest::Approx(1.0)); // C2

    // Check D matrix
    CHECK(result.D(0, 0) == doctest::Approx(0.0)); // D2*D1 = 0*0
}

TEST_CASE("StateSpace Series via Operator*") {
    StateSpace<1, 1, 1, double, 1, 1> sys1{
        .A = {{-1.0}},
        .B = {{1.0}},
        .C = {{1.0}},
        .D = {{0.0}},
    };

    StateSpace<1, 1, 1, double, 1, 1> sys2{
        .A = {{-2.0}},
        .B = {{1.0}},
        .C = {{1.0}},
        .D = {{0.0}},
    };

    // sys1 * sys2 should be equivalent to series(sys2, sys1)
    auto result = *(sys1 * sys2);

    CHECK(result.A.rows() == 2);
    CHECK(result.A(0, 0) == doctest::Approx(-2.0));
    CHECK(result.A(1, 0) == doctest::Approx(1.0));
}

TEST_CASE("StateSpace Parallel Connection") {
    StateSpace<1, 1, 1, double, 1, 1> sys1{
        .A = {{-1.0}},
        .B = {{1.0}},
        .C = {{2.0}},
        .D = {{0.5}},
    };

    StateSpace<1, 1, 1, double, 1, 1> sys2{
        .A = {{-3.0}},
        .B = {{1.0}},
        .C = {{1.0}},
        .D = {{0.3}},
    };

    auto result = *parallel(sys1, sys2);

    // Expected dimensions
    CHECK(result.A.rows() == 2);
    CHECK(result.A.cols() == 2);
    CHECK(result.B.rows() == 2);
    CHECK(result.B.cols() == 1);
    CHECK(result.C.rows() == 1);
    CHECK(result.C.cols() == 2);
    CHECK(result.D.rows() == 1);
    CHECK(result.D.cols() == 1);

    // Check A matrix (block diagonal)
    CHECK(result.A(0, 0) == doctest::Approx(-1.0));
    CHECK(result.A(0, 1) == doctest::Approx(0.0));
    CHECK(result.A(1, 0) == doctest::Approx(0.0));
    CHECK(result.A(1, 1) == doctest::Approx(-3.0));

    // Check B matrix (stacked)
    CHECK(result.B(0, 0) == doctest::Approx(1.0));
    CHECK(result.B(1, 0) == doctest::Approx(1.0));

    // Check C matrix (side-by-side)
    CHECK(result.C(0, 0) == doctest::Approx(2.0));
    CHECK(result.C(0, 1) == doctest::Approx(1.0));

    // Check D matrix (summed)
    CHECK(result.D(0, 0) == doctest::Approx(0.8)); // 0.5 + 0.3
}

TEST_CASE("StateSpace Parallel via Operator+") {
    StateSpace<1, 1, 1, double, 1, 1> sys1{
        .A = {{-1.0}},
        .B = {{1.0}},
        .C = {{2.0}},
        .D = {{0.5}},
    };

    StateSpace<1, 1, 1, double, 1, 1> sys2{
        .A = {{-3.0}},
        .B = {{1.0}},
        .C = {{1.0}},
        .D = {{0.3}},
    };

    auto result = *(sys1 + sys2);

    CHECK(result.A.rows() == 2);
    CHECK(result.D(0, 0) == doctest::Approx(0.8));
}

TEST_CASE("StateSpace Subtraction Connection") {
    StateSpace<1, 1, 1, double, 1, 1> sys1{
        .A = {{-1.0}},
        .B = {{1.0}},
        .C = {{2.0}},
        .D = {{0.5}},
    };

    StateSpace<1, 1, 1, double, 1, 1> sys2{
        .A = {{-3.0}},
        .B = {{1.0}},
        .C = {{1.0}},
        .D = {{0.3}},
    };

    auto result = *subtract(sys1, sys2);

    // Expected dimensions
    CHECK(result.A.rows() == 2);
    CHECK(result.A.cols() == 2);
    CHECK(result.B.rows() == 2);
    CHECK(result.B.cols() == 1);
    CHECK(result.C.rows() == 1);
    CHECK(result.C.cols() == 2);
    CHECK(result.D.rows() == 1);
    CHECK(result.D.cols() == 1);

    // Check A matrix (block diagonal)
    CHECK(result.A(0, 0) == doctest::Approx(-1.0));
    CHECK(result.A(0, 1) == doctest::Approx(0.0));
    CHECK(result.A(1, 0) == doctest::Approx(0.0));
    CHECK(result.A(1, 1) == doctest::Approx(-3.0));

    // Check B matrix (stacked)
    CHECK(result.B(0, 0) == doctest::Approx(1.0));
    CHECK(result.B(1, 0) == doctest::Approx(1.0));

    // Check C matrix (sys1 minus sys2)
    CHECK(result.C(0, 0) == doctest::Approx(2.0));
    CHECK(result.C(0, 1) == doctest::Approx(-1.0)); // -C2

    // Check D matrix (subtracted)
    CHECK(result.D(0, 0) == doctest::Approx(0.2)); // 0.5 - 0.3
}

TEST_CASE("StateSpace Subtraction via Operator-") {
    StateSpace<1, 1, 1, double, 1, 1> sys1{
        .A = {{-1.0}},
        .B = {{1.0}},
        .C = {{2.0}},
        .D = {{0.5}},
    };

    StateSpace<1, 1, 1, double, 1, 1> sys2{
        .A = {{-3.0}},
        .B = {{1.0}},
        .C = {{1.0}},
        .D = {{0.3}},
    };

    auto result = *(sys1 - sys2);

    CHECK(result.A.rows() == 2);
    CHECK(result.C(0, 0) == doctest::Approx(2.0));
    CHECK(result.C(0, 1) == doctest::Approx(-1.0));
    CHECK(result.D(0, 0) == doctest::Approx(0.2));
}

TEST_CASE("StateSpace Negative Feedback") {
    // Simple plant: G(s) = 1/(s+1)
    StateSpace<1, 1, 1, double, 1, 1> plant{
        .A = {{-1.0}},
        .B = {{1.0}},
        .C = {{1.0}},
        .D = {{0.0}},
    };

    // Dynamic feedback with B2 ≠ 1 (catches the old B = B1*B2 bug)
    StateSpace<1, 1, 1, double, 1, 1> controller{
        .A = {{0.0}},
        .B = {{3.0}},
        .C = {{2.0}},
        .D = {{0.0}},
    };

    auto result = *feedback(plant, controller);

    // Expected dimensions
    CHECK(result.A.rows() == 2);
    CHECK(result.A.cols() == 2);
    CHECK(result.B.rows() == 2);
    CHECK(result.B.cols() == 1);
    CHECK(result.C.rows() == 1);
    CHECK(result.C.cols() == 2);
    CHECK(result.D.rows() == 1);
    CHECK(result.D.cols() == 1);

    // D1=D2=0: A = [A1, -B1*C2; B2*C1, A2] = [-1, -2; 3, 0]
    CHECK(result.A(0, 0) == doctest::Approx(-1.0));
    CHECK(result.A(0, 1) == doctest::Approx(-2.0));
    CHECK(result.A(1, 0) == doctest::Approx(3.0));
    CHECK(result.A(1, 1) == doctest::Approx(0.0));

    // B = [B1; 0] — independent of B2 (must be 1, not B1*B2=3)
    CHECK(result.B(0, 0) == doctest::Approx(1.0));
    CHECK(result.B(1, 0) == doctest::Approx(0.0));

    // C = [C1, 0], D = 0
    CHECK(result.C(0, 0) == doctest::Approx(1.0));
    CHECK(result.C(0, 1) == doctest::Approx(0.0));
    CHECK(result.D(0, 0) == doctest::Approx(0.0));
}

TEST_CASE("StateSpace feedback static gain K=2 closes pole at -3") {
    // Plant G(s) = 1/(s+1); negative feedback of static gain K=2
    // Closed loop: G/(1+KG) = 1/(s+3), pole at -3
    StateSpace<1, 1, 1> plant{
        .A = {{-1.0}},
        .B = {{1.0}},
        .C = {{1.0}},
        .D = {{0.0}},
    };

    // Static gain as D-only system (no states)
    StateSpace<0, 1, 1> K{
        .D = {{2.0}},
    };

    auto cl = *feedback(plant, K);

    CHECK(cl.A.rows() == 1);
    CHECK(cl.A(0, 0) == doctest::Approx(-3.0));
    CHECK(cl.B(0, 0) == doctest::Approx(1.0));
    CHECK(cl.C(0, 0) == doctest::Approx(1.0));
    CHECK(cl.D(0, 0) == doctest::Approx(0.0));
}

TEST_CASE("StateSpace Feedback via Operator/") {
    StateSpace<1, 1, 1, double, 1, 1> plant{
        .A = {{-1.0}},
        .B = {{1.0}},
        .C = {{1.0}},
        .D = {{0.0}},
    };

    StateSpace<1, 1, 1, double, 1, 1> controller{
        .A = {{0.0}},
        .B = {{3.0}},
        .C = {{2.0}},
        .D = {{0.0}},
    };

    auto result = *(plant / controller);

    CHECK(result.A.rows() == 2);
    CHECK(result.A(0, 0) == doctest::Approx(-1.0));
    CHECK(result.A(0, 1) == doctest::Approx(-2.0));
    CHECK(result.B(0, 0) == doctest::Approx(1.0)); // not B1*B2
}

TEST_CASE("StateSpace Constexpr") {
    constexpr StateSpace<1, 1, 1, double, 1, 1> sys1{
        .A = {{-1.0}},
        .B = {{1.0}},
        .C = {{1.0}},
        .D = {{0.0}},
        .G = {},
        .H = {},
        .Ts = 0.0
    };

    constexpr StateSpace<1, 1, 1, double, 1, 1> sys2{
        .A = {{-2.0}},
        .B = {{1.0}},
        .C = {{1.0}},
        .D = {{0.0}},
        .G = {},
        .H = {},
        .Ts = 0.0
    };

    // Should compile as constexpr
    constexpr auto result_series = *series(sys1, sys2);
    constexpr auto result_parallel = *parallel(sys1, sys2);
    constexpr auto result_feedback = *feedback(sys1, sys2);

    static_assert(result_series.A.rows() == 2);
    static_assert(result_parallel.A.rows() == 2);
    static_assert(result_feedback.A.rows() == 2);
}

TEST_CASE("StateSpace interconnect rejects mismatched sample times") {
    StateSpace<1, 1, 1> cont{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}, .Ts = 0.0};
    StateSpace<1, 1, 1> disc{.A = {{-0.9}}, .B = {{1.0}}, .C = {{1.0}}, .Ts = 0.01};
    CHECK_FALSE(series(cont, disc).has_value());
    CHECK_FALSE(parallel(cont, disc).has_value());
    CHECK_FALSE(feedback(cont, disc).has_value());
}

TEST_CASE("StateSpace feedback rejects singular algebraic loop") {
    // D1 = 1, D2 = -1 ⇒ I + D2 D1 = 0 (singular)
    StateSpace<1, 1, 1> sys1{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}, .D = {{1.0}}};
    StateSpace<1, 1, 1> sys2{.A = {{-2.0}}, .B = {{1.0}}, .C = {{1.0}}, .D = {{-1.0}}};
    CHECK_FALSE(feedback(sys1, sys2).has_value());
}

TEST_CASE("StateSpace * matches MATLAB series order") {
    // MATLAB: sys1*sys2 ≡ series(sys2, sys1) ≡ u → sys2 → sys1
    StateSpace<1, 1, 1> a{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}};
    StateSpace<1, 1, 1> b{.A = {{-2.0}}, .B = {{1.0}}, .C = {{1.0}}};
    auto                via_star = *(a * b);
    auto                via_series = *series(b, a);
    CHECK(via_star.A(0, 0) == doctest::Approx(via_series.A(0, 0)));
    CHECK(via_star.A(1, 1) == doctest::Approx(via_series.A(1, 1)));
}

TEST_CASE("Statespace 2x2 with noise") {
    StateSpace<2, 1, 1, double, 2, 1> sys = {
        .A = {{0.0, 1.0}, {-2.0, -3.0}},
        .B = {{0.0}, {1.0}},
        .C = {{1.0, 0.0}},
        .D = {{0.0}},
        .G = {{0.0, 0.0}, {1.0, 0.0}}, // Process noise input
        .H = {{0.1}},                  // Measurement noise input
        .Ts = 0.0
    };

    CHECK(sys.A.rows() == 2);
    CHECK(sys.B.cols() == 1);
    CHECK(sys.G.cols() == 2);
    CHECK(sys.H.cols() == 1);

    CHECK(sys.A(1, 0) == doctest::Approx(-2.0));
    CHECK(sys.G(1, 0) == doctest::Approx(1.0));
    CHECK(sys.H(0, 0) == doctest::Approx(0.1));
}

TEST_CASE("series and parallel preserve process/measurement noise G/H") {
    // Two SISO first-order plants with independent noise dims.
    StateSpace<1, 1, 1, double, 1, 1> a{
        .A = {{-1.0}},
        .B = {{1.0}},
        .C = {{1.0}},
        .D = {{0.0}},
        .G = {{2.0}},
        .H = {{0.3}},
        .Ts = 0.0,
    };
    StateSpace<1, 1, 1, double, 1, 1> b{
        .A = {{-2.0}},
        .B = {{1.0}},
        .C = {{1.0}},
        .D = {{0.0}},
        .G = {{4.0}},
        .H = {{0.5}},
        .Ts = 0.0,
    };

    const auto ser = *series(a, b); // u → a → b → y; NW=2, NV=2
    CHECK(ser.G.rows() == 2);
    CHECK(ser.G.cols() == 2);
    CHECK(ser.H.rows() == 1);
    CHECK(ser.H.cols() == 2);
    // G block-diag path: sys1.G on top-left; sys2.G bottom-right
    CHECK(ser.G(0, 0) == doctest::Approx(2.0));
    CHECK(ser.G(1, 1) == doctest::Approx(4.0));
    // H: [sys2.D * sys1.H | sys2.H] = [0 | 0.5] when D2=0
    CHECK(ser.H(0, 0) == doctest::Approx(0.0));
    CHECK(ser.H(0, 1) == doctest::Approx(0.5));

    const auto par = *parallel(a, b);
    CHECK(par.G(0, 0) == doctest::Approx(2.0));
    CHECK(par.G(1, 1) == doctest::Approx(4.0));
    CHECK(par.H(0, 0) == doctest::Approx(0.3));
    CHECK(par.H(0, 1) == doctest::Approx(0.5));
}

TEST_CASE("select slices a MIMO plant") {
    constexpr StateSpace<1, 1, 2> P{
        .A = {{-1.0}},
        .B = {{1.0}},
        .C = {{1.0}, {2.0}},
        .D = {{0.0}, {0.5}},
    };
    constexpr auto G0 = *select(P, 0, 0);
    constexpr auto G1 = *select(P, 1, 0);
    static_assert(G0.C(0, 0) == 1.0);
    static_assert(G1.D(0, 0) == 0.5);
    CHECK(G1.C(0, 0) == doctest::Approx(2.0));
    CHECK_FALSE(select(P, 2, 0).has_value());
}

TEST_CASE("mix of y into u is G/(1-G)") {
    // G = 1/(s+1), two identical outputs. v = u + y0 ⇒ y0/u = 1/s.
    constexpr StateSpace<1, 1, 2> P{
        .A = {{-1.0}},
        .B = {{1.0}},
        .C = {{1.0}, {1.0}},
        .D = {{0.0}, {0.0}},
    };
    constexpr damp::array<SsFeed<double>, 1> ff{SsFeed<double>{0, 1.0}};
    constexpr auto                           mixed = *mix(P, 0, 1.0, ff);
    constexpr auto                           Gid = *select(mixed, 0, 0);
    static_assert(Gid.A(0, 0) == 0.0);
    CHECK(Gid.B(0, 0) == doctest::Approx(1.0));
    CHECK(Gid.C(0, 0) == doctest::Approx(1.0));
    CHECK(Gid.D(0, 0) == doctest::Approx(0.0));
}

TEST_CASE("ss_gain / ss_integrator / mix_outputs are constexpr") {
    constexpr auto K = ss_gain(3.0);
    constexpr auto I = ss_integrator(2.0);
    static_assert(K.D(0, 0) == 3.0);
    static_assert(I.C(0, 0) == 2.0);
    constexpr StateSpace<1, 1, 2> P{
        .A = {{-1.0}},
        .B = {{1.0}},
        .C = {{1.0}, {4.0}},
        .D = {{0.0}, {0.0}},
    };
    constexpr damp::array<SsFeed<double>, 2> terms{
        SsFeed<double>{0, 1.0},
        SsFeed<double>{1, -0.25},
    };
    constexpr auto y = *mix_outputs(P, terms);
    CHECK(y.C(0, 0) == doctest::Approx(0.0));
}