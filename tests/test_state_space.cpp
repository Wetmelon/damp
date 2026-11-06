// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>
#include <limits>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
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

TEST_CASE("subtract negates the measurement noise of the subtracted system") {
    const StateSpace<1, 1, 1, double, 1, 1> sys1{.A = {{-1.0}}, .B = {{1.0}}, .C = {{2.0}}, .H = {{0.5}}};
    const StateSpace<1, 1, 1, double, 1, 1> sys2{.A = {{-3.0}}, .B = {{1.0}}, .C = {{1.0}}, .H = {{0.25}}};
    const auto                              diff = subtract(sys1, sys2);
    REQUIRE(diff.has_value());
    CHECK(diff->H(0, 0) == 0.5);
    CHECK(diff->H(0, 1) == -0.25); // y₁ − y₂ = … − H₂v₂
    CHECK(diff->G(0, 0) == 1.0);
    CHECK(diff->G(1, 1) == 1.0);
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

TEST_CASE("StateSpace Feedback with a dynamic controller") {
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

    auto result = *feedback(plant, controller);

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

TEST_CASE("feedback refuses an algebraic loop lost to cancellation, relative to the data") {
    const double eps = std::numeric_limits<double>::epsilon();
    for (const double scale : {1.0, 1048576.0}) {
        CAPTURE(scale);
        const StateSpace<1, 1, 1> sys1{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}, .D = {{scale}}};
        CHECK_FALSE(feedback(sys1, ss_gain(-(1.0 + (4.0 * eps)) / scale)).has_value()); // 1 + D₂D₁ = −4·eps (LU alone accepts it)
        CHECK_FALSE(feedback(sys1, ss_gain((-1.0 + 1e-15) / scale)).has_value());
        const auto ok = feedback(sys1, ss_gain(-0.25 / scale)); // well posed: 1 + D₂D₁ = 3/4
        REQUIRE(ok.has_value());
        CHECK(ok->D(0, 0) == doctest::Approx(scale * 4.0 / 3.0));
    }
}

TEST_CASE("feedback accepts MIMO loops whose channels do not interact, whatever their gains") {
    //! The scale is the size of the terms D₂·D₁ is formed from, not ‖D₂‖·‖D₁‖.
    const auto loop = []<typename T>(T d1a, T d1b, T d2a, T d2b) {
        const StateSpace<1, 2, 2, T> plant{.A = {{T{-1}}}, .B = {{T{1}, T{1}}}, .C = {{T{1}}, {T{1}}}, .D = {{d1a, T{0}}, {T{0}, d1b}}};
        const StateSpace<1, 2, 2, T> ctrl{.A = {{T{-2}}}, .B = {{T{1}, T{1}}}, .C = {{T{1}}, {T{1}}}, .D = {{d2a, T{0}}, {T{0}, d2b}}};
        return feedback(plant, ctrl);
    };
    CHECK(loop(1e7, 0.0, 0.0, 1e7).has_value());                                                                 // D₂D₁ = 0
    CHECK(loop(10.0F, 0.0F, 1.0F, 1e4F).has_value());                                                            // I + D₂D₁ = diag(11, 1)
    CHECK(loop(1e3F, 1e-3F, 1e-3F, 1e3F).has_value());                                                           // I + D₂D₁ = 2I, units apart
    CHECK_FALSE(loop(1e7, 0.0, (-1.0 - (4.0 * std::numeric_limits<double>::epsilon())) / 1e7, 0.0).has_value()); // 1 + D₂D₁ = −4·eps on channel 1
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

TEST_CASE("mix gives the exact output map of the measurement noise through the loop") {
    //! y = x + 0.5·v + n, v = u + 0.4·y: y = … + (1 + 0.5·0.4/(1 − 0.4·0.5))·n = … + 1.25·n.
    const StateSpace<1, 1, 1, double, 0, 1> P{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}, .D = {{0.5}}, .H = {{1.0}}};
    const auto                              mixed = mix(P, 0, 1.0, damp::array<SsFeed<double>, 1>{SsFeed<double>{0, 0.4}});
    REQUIRE(mixed.has_value());
    CHECK(mixed->H(0, 0) == doctest::Approx(1.25));
    CHECK(mixed->D(0, 0) == doctest::Approx(0.625));
}

TEST_CASE("mix_outputs with one measurement-noise input starts H from zero") {
    const StateSpace<2, 1, 2, double, 0, 1> P{
        .A = {{-1.0, 0.0}, {0.0, -2.0}},
        .B = {{1.0}, {1.0}},
        .C = {{1.0, 0.0}, {0.0, 1.0}},
        .H = {{0.5}, {0.25}},
    };
    const auto y = mix_outputs(P, damp::array<SsFeed<double>, 1>{SsFeed<double>{1, 2.0}});
    REQUIRE(y.has_value());
    CHECK(y->H(0, 0) == 0.5); // 2·0.25, not 1 + 2·0.25
    CHECK(y->C(0, 1) == 2.0);
}

TEST_CASE("feedback_read closes the loop on G_loop and reads G_out, noise map included") {
    //! u = r − k·y_loop with k = 2, 1 + k·D_loop = 2:
    //! ẋ = −2x + 0.5r, y_out = 1.75x + 0.125r + (0.3 − 0.25)v.
    const StateSpace<1, 1, 1, double, 0, 1> G_loop{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}, .D = {{0.5}}, .H = {{1.0}}};
    const StateSpace<1, 1, 1, double, 0, 1> G_out{.A = {{-1.0}}, .B = {{1.0}}, .C = {{2.0}}, .D = {{0.25}}, .H = {{0.3}}};
    const auto                              cl = feedback_read(G_loop, G_out, 2.0);
    REQUIRE(cl.has_value());
    CHECK(cl->A(0, 0) == doctest::Approx(-2.0));
    CHECK(cl->B(0, 0) == doctest::Approx(0.5));
    CHECK(cl->C(0, 0) == doctest::Approx(1.75));
    CHECK(cl->D(0, 0) == doctest::Approx(0.125));
    CHECK(cl->H(0, 0) == doctest::Approx(0.05));

    //! Not one realization: another A, another B, another sample time.
    StateSpace<1, 1, 1, double, 0, 1> other = G_out;
    other.A(0, 0) = -5.0;
    CHECK_FALSE(feedback_read(G_loop, other, 2.0).has_value());
    other = G_out;
    other.B(0, 0) = 3.0;
    CHECK_FALSE(feedback_read(G_loop, other, 2.0).has_value());
    other = G_out;
    other.Ts = 0.1;
    CHECK_FALSE(feedback_read(G_loop, other, 2.0).has_value());
    CHECK_FALSE(feedback_read(G_loop, G_out, -2.0).has_value()); // 1 + k·D_loop = 0
}

TEST_CASE("ss_gain / mix_outputs are constexpr") {
    constexpr auto K = ss_gain(3.0);
    static_assert(K.D(0, 0) == 3.0);
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
TEST_CASE("StateSpace deduces its sizes and scalar from braces") {
    constexpr StateSpace first_order{.A = {{-1.0 / 30.0}}, .B = {{6.0 / 30.0}}, .C = {{1.0}}};
    static_assert(std::is_same_v<decltype(first_order), const StateSpace<1, 1, 1, double>>);
    static_assert(first_order.B(0, 0) == 6.0 / 30.0);

    constexpr StateSpace with_d{.A = {{0.9, 0.1}, {0.0, 0.8}}, .B = {{1.0, 0.0}, {0.0, 1.0}}, .C = {{1.0, 0.0}}, .D = {{0.0, 0.5}}};
    static_assert(std::is_same_v<decltype(with_d), const StateSpace<2, 2, 1, double>>);
    static_assert(with_d.D(0, 1) == 0.5 && with_d.Ts == 0.0);

    constexpr StateSpace discrete{.A = {{0.9f}}, .B = {{1.0f}}, .C = {{1.0f}}, .Ts = 0.01f};
    static_assert(std::is_same_v<decltype(discrete), const StateSpace<1, 1, 1, float>>);
    static_assert(discrete.Ts == 0.01f);

    constexpr StateSpace all{.A = {{0.5}}, .B = {{1.0, 2.0}}, .C = {{1.0}, {3.0}}, .D = {{0.0, 0.0}, {0.0, 1.0}}, .Ts = 0.1};
    static_assert(std::is_same_v<decltype(all), const StateSpace<1, 2, 2, double>>);
    static_assert(all.C(1, 0) == 3.0 && all.D(1, 1) == 1.0 && all.Ts == 0.1);
    CHECK(all.is_discrete());
}

namespace {
// A lightly damped plant and a lead controller, both with feedthrough so the algebraic loop matters.
constexpr StateSpace<2, 1, 1, double, 1, 1> lft_plant() {
    return StateSpace<2, 1, 1, double, 1, 1>{
        .A = {{0.0, 1.0}, {-4.0, -0.4}},
        .B = {{0.0}, {1.0}},
        .C = {{1.0, 0.25}},
        .D = {{0.2}},
        .G = {{0.0}, {0.5}},
        .H = {{0.1}},
    };
}

constexpr StateSpace<1, 1, 1, double, 1, 1> lft_controller() {
    return StateSpace<1, 1, 1, double, 1, 1>{.A = {{-5.0}}, .B = {{2.0}}, .C = {{-3.0}}, .D = {{1.5}}, .G = {{0.7}}, .H = {{0.05}}};
}

// −K: the same controller with its output negated (lft closes u = K y, feedback u = r − K y).
template<size_t NX, size_t NU, size_t NY, size_t NW, size_t NV>
constexpr StateSpace<NX, NU, NY, double, NW, NV> negated(StateSpace<NX, NU, NY, double, NW, NV> K) {
    K.C = K.C * -1.0;
    K.D = K.D * -1.0;
    K.H = K.H * -1.0;
    return K;
}

template<size_t R, size_t C>
void check_matrix_close(const Matrix<R, C>& a, const Matrix<R, C>& b) {
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            CHECK(a(i, j) == doctest::Approx(b(i, j)).epsilon(1e-12).scale(1.0));
        }
    }
}
} // namespace

TEST_SUITE("lft (lower linear fractional transformation)") {
    TEST_CASE("with no w / z channels, lft(G, K) closes the loop feedback(G, −K) closes") {
        constexpr auto G = lft_plant();
        constexpr auto K = lft_controller();
        const auto     closed = lft(G, K);
        const auto     fb = feedback(G, negated(K));
        REQUIRE(closed.has_value());
        REQUIRE(fb.has_value());
        static_assert(std::is_same_v<decltype(closed), const damp::optional<StateSpace<3, 0, 0, double, 2, 2>>>);
        check_matrix_close(closed->A, fb->A);
        check_matrix_close(closed->G, fb->G);
    }

    TEST_CASE("P = [G G; G G] with −K is feedback(G, K), noise channels included") {
        constexpr auto                    G = lft_plant();
        constexpr auto                    K = lft_controller();
        StateSpace<2, 2, 2, double, 1, 1> P{};
        P.A = G.A;
        P.B = Matrix<2, 2>{{G.B(0, 0), G.B(0, 0)}, {G.B(1, 0), G.B(1, 0)}}; // [w; u] both enter as G's input
        P.C = Matrix<2, 2>{{G.C(0, 0), G.C(0, 1)}, {G.C(0, 0), G.C(0, 1)}}; // z = y = G's output
        P.D = Matrix<2, 2>{{G.D(0, 0), G.D(0, 0)}, {G.D(0, 0), G.D(0, 0)}};
        P.G = G.G;
        P.H = Matrix<2, 1>{{G.H(0, 0)}, {G.H(0, 0)}};

        const auto closed = lft(P, negated(K)); // u = −K y: G's input is w − K y
        const auto fb = feedback(G, K);
        REQUIRE(closed.has_value());
        REQUIRE(fb.has_value());
        check_matrix_close(closed->A, fb->A);
        check_matrix_close(closed->B, fb->B);
        check_matrix_close(closed->C, fb->C);
        check_matrix_close(closed->D, fb->D);
        check_matrix_close(closed->G, fb->G);
        check_matrix_close(closed->H, fb->H);
    }

    TEST_CASE("a 2-input / 2-output plant matches python-control's lft (MATLAB's lft convention)") {
        //! P: [w; u] → [z; y]; K: y → u with feedthrough, so W = 1/(1 − D22·DK) = 1/1.28.
        constexpr StateSpace<2, 2, 2> P{
            .A = {{0.0, 1.0}, {-2.0, -0.5}},
            .B = {{0.0, 0.0}, {1.0, 2.0}},
            .C = {{1.0, 0.0}, {0.5, 1.0}},
            .D = {{0.1, 0.2}, {0.3, 0.4}},
        };
        constexpr StateSpace<1, 1, 1> K{.A = {{-3.0}}, .B = {{1.5}}, .C = {{-2.0}}, .D = {{-0.7}}};
        constexpr auto                closed = lft(P, K);
        static_assert(closed.has_value() && closed->Ts == 0.0);

        //! control.ss(A, B, C, D).lft(control.ss(AK, BK, CK, DK)), python-control 0.10.2.
        check_matrix_close(closed->A, Matrix<3, 3>{{0.0, 1.0, 0.0}, {-2.546875, -1.59375, -3.125}, {0.5859375, 1.171875, -3.9375}});
        check_matrix_close(closed->B, Matrix<3, 1>{{0.0}, {0.671875}, {0.3515625}});
        check_matrix_close(closed->C, Matrix<1, 3>{{0.9453125, -0.109375, -0.3125}});
        check_matrix_close(closed->D, Matrix<1, 1>{{0.0671875}});

        //! And the frequency response is P11 + P12 K (1 − P22 K)⁻¹ P21: positive feedback.
        for (const double w : {0.3, 1.0, 7.0}) {
            const damp::complex<double> s{0.0, w};
            const auto                  Ps = *eval_frf(P, s);
            const auto                  Ks = (*eval_frf(K, s))(0, 0);
            const auto                  expected = Ps(0, 0) + (Ps(0, 1) * Ks * Ps(1, 0) / (damp::complex<double>{1.0, 0.0} - (Ps(1, 1) * Ks)));
            const auto                  got = (*eval_frf(*closed, s))(0, 0);
            CHECK((got - expected).abs() <= 1e-12 * expected.abs());
        }
    }

    TEST_CASE("lft refuses an ill-posed loop and mismatched sample times") {
        //! w and u both drive the one state; z = y = x, and y also reads u through D22.
        const auto check_loops = [](double d22_scale) {
            StateSpace<1, 2, 2> P{.A = {{-1.0}}, .B = {{1.0, 1.0}}, .C = {{1.0}, {1.0}}, .D = {{0.0, 0.0}, {0.0, d22_scale}}};
            StateSpace<0, 1, 1> K{};
            const double        eps = std::numeric_limits<double>::epsilon();

            K.D(0, 0) = 1.0 / d22_scale; // I − D22·DK = 0
            CHECK_FALSE(lft(P, K).has_value());
            K.D(0, 0) = (1.0 + (4.0 * eps)) / d22_scale; // −4·eps, lost to cancellation (LU alone accepts it)
            CHECK_FALSE(lft(P, K).has_value());
            K.D(0, 0) = 0.25 / d22_scale; // well posed: W = 4/3
            const auto ok = lft(P, K);
            REQUIRE(ok.has_value());
            CHECK(ok->A(0, 0) == doctest::Approx(-1.0 + (K.D(0, 0) * 4.0 / 3.0))); // A + B2·DK·W·C2
            CHECK(ok->B(0, 0) == doctest::Approx(1.0));
        };
        check_loops(1.0);
        //! D22 scaled by 2²⁰ and DK by 2⁻²⁰: the same verdicts, the test is relative to the data.
        check_loops(1048576.0);

        const StateSpace<1, 2, 2> P{.A = {{-1.0}}, .B = {{1.0, 1.0}}, .C = {{1.0}, {1.0}}};
        const StateSpace<1, 1, 1> K_discrete{.A = {{0.5}}, .B = {{1.0}}, .C = {{1.0}}, .Ts = 0.01};
        CHECK_FALSE(lft(P, K_discrete).has_value());
    }
}

TEST_CASE("series does not route measurement noise into process-noise columns") {
    // u → a → b. a's measurement noise v₁ is not a process-noise input w₁.
    const StateSpace<1, 1, 1, double, 1, 2> a{
        .A = {{-1.0}},
        .B = {{1.0}},
        .C = {{1.0}},
        .G = {{2.0}},
        .H = {{0.3, 0.7}},
    };
    const StateSpace<1, 1, 1, double, 1, 1> b{
        .A = {{-2.0}},
        .B = {{1.0}},
        .C = {{1.0}},
        .G = {{4.0}},
        .H = {{0.5}},
    };
    const auto ser = series(a, b);
    REQUIRE(ser.has_value());
    CHECK(ser->G(0, 0) == 2.0);
    CHECK(ser->G(1, 0) == 0.0); // w₁ only drives x₁
    CHECK(ser->G(1, 1) == 4.0);

    const auto fb = feedback(a, b);
    REQUIRE(fb.has_value());
    CHECK(fb->G(1, 0) == 0.0);
    CHECK(fb->H(0, 0) == doctest::Approx(0.3)); // D₁ = 0: y = C₁x₁ + H₁v₁
    CHECK(fb->H(0, 2) == 0.0);                  // v₂ reaches y only through x₁
}

TEST_CASE("sample-time match is relative") {
    const StateSpace<1, 1, 1> a{.A = {{0.5}}, .B = {{1.0}}, .C = {{1.0}}, .Ts = 1e-13};
    const StateSpace<1, 1, 1> b{.A = {{0.5}}, .B = {{1.0}}, .C = {{1.0}}, .Ts = 2e-13};
    CHECK_FALSE(series(a, b).has_value());
}

TEST_CASE("mix accepts a high but finite algebraic-loop gain") {
    // v = u + y with D = 1 − 1e-10: 1 − D is small but not lost to roundoff.
    const StateSpace<1, 1, 1> P{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}, .D = {{1.0 - 1e-13}}};
    const auto                m = mix(P, 0, 1.0, damp::array<SsFeed<double>, 1>{SsFeed<double>{0, 1.0}});
    CHECK(m.has_value());
}
