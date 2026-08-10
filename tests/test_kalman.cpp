// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>

#include "damp/estimation/kalman.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

TEST_CASE("Kalman predict/update 1D") {
    StateSpace<1, 1, 1, double, 1, 1> sys{};
    sys.A(0, 0) = 1.0;
    sys.B(0, 0) = 1.0;
    sys.C(0, 0) = 1.0;
    sys.G(0, 0) = 1.0;
    sys.H(0, 0) = 1.0;
    sys.Ts = 1.0; // discrete system

    KalmanFilter<1, 1, 1, double, 1, 1> kf(sys, Matrix<1, 1, double>{{0.1}}, Matrix<1, 1, double>{{0.25}}, ColVec<1>{0.0}, Matrix<1, 1, double>::identity());

    ColVec<1> u = {0.0};
    kf.predict(u);

    ColVec<1> y = {1.2};
    bool      ok = kf.update(y, u);
    CHECK(ok);
    CHECK(kf.state()[0] == doctest::Approx(0.978).epsilon(1e-3));
    CHECK(kf.covariance()(0, 0) == doctest::Approx(0.204).epsilon(1e-3));
}

TEST_CASE("KalmanFilter from KalmanResult + fused estimate()") {
    StateSpace<1, 1, 1, double, 1, 1> sys{};
    sys.A(0, 0) = 0.95;
    sys.B(0, 0) = 1.0;
    sys.C(0, 0) = 1.0;
    sys.G(0, 0) = 1.0;
    sys.H(0, 0) = 1.0;
    sys.Ts = 1.0;

    const auto result = design::kalman(sys, Matrix<1, 1>{{0.1}}, Matrix<1, 1>{{0.25}});
    REQUIRE(result.success);

    KalmanFilter<1, 1, 1, double, 1, 1> kf{result};
    CHECK(kf.covariance()(0, 0) == doctest::Approx(result.P(0, 0)));
    CHECK(kf.model().A(0, 0) == doctest::Approx(0.95));

    const ColVec<1> y{1.0};
    const ColVec<1> u{0.0};
    const auto&     x = kf.estimate(y, u);
    CHECK(std::isfinite(x[0]));
    // After a measurement of 1 from a near-zero seed, estimate moves toward y
    CHECK(x[0] > 0.0);
    CHECK(x[0] < 1.5);
}

TEST_CASE("SteadyStateKalmanFilter uses designed L and set_gain/reset") {
    StateSpace<1, 1, 1, double, 1, 1> sys{};
    sys.A(0, 0) = 0.9;
    sys.B(0, 0) = 1.0;
    sys.C(0, 0) = 1.0;
    sys.G(0, 0) = 1.0;
    sys.H(0, 0) = 1.0;
    sys.Ts = 1.0;

    const auto result = design::kalman(sys, Matrix<1, 1>{{0.1}}, Matrix<1, 1>{{0.1}});
    REQUIRE(result.success);

    SteadyStateKalmanFilter<1, 1, 1, double, 1, 1> sskf{result};
    CHECK(sskf.gain()(0, 0) == doctest::Approx(result.L(0, 0)));
    CHECK(sskf.covariance()(0, 0) == doctest::Approx(result.P(0, 0)));

    const auto& x = sskf.estimate(ColVec<1>{2.0}, ColVec<1>{0.0});
    CHECK(x[0] != 0.0);

    sskf.set_gain(Matrix<1, 1>{{0.0}}); // open-loop (no correction)
    sskf.reset(ColVec<1>{0.0});
    sskf.predict(ColVec<1>{0.0});
    REQUIRE(sskf.update(ColVec<1>{5.0}));
    CHECK(sskf.state()[0] == doctest::Approx(0.0)); // L=0 → state unchanged by y
}

TEST_CASE("Kalman design R!=0, square C (NY=NX=2)") {
    // Double integrator, dt=0.1: position+velocity state, both measured directly.
    // Reference values generated with scipy.linalg.solve_discrete_are.
    constexpr double                  dt = 0.1;
    StateSpace<2, 1, 2, double, 2, 2> sys{};
    sys.A = {{1.0, dt}, {0.0, 1.0}};
    sys.B = {{0.5 * dt * dt}, {dt}};
    sys.C = Matrix<2, 2, double>::identity();
    sys.G = Matrix<2, 2, double>::identity();
    sys.H = Matrix<2, 2, double>::identity();
    sys.Ts = dt;

    Matrix<2, 2> Q_proc = {{1e-4, 0.0}, {0.0, 1e-2}};
    Matrix<2, 2> R_meas = {{1e-2, 0.0}, {0.0, 5e-2}};

    auto result = design::kalman(sys, Q_proc, R_meas);
    REQUIRE(result.success);

    // scipy: solve_discrete_are(A.T, C.T, Q, R)
    CHECK(result.P(0, 0) == doctest::Approx(0.002606406390998).epsilon(1e-8));
    CHECK(result.P(0, 1) == doctest::Approx(0.003583502847688).epsilon(1e-8));
    CHECK(result.P(1, 0) == doctest::Approx(0.003583502847688).epsilon(1e-8));
    CHECK(result.P(1, 1) == doctest::Approx(0.027171128304431).epsilon(1e-8));

    CHECK(result.L(0, 0) == doctest::Approx(0.196141710450478).epsilon(1e-8));
    CHECK(result.L(0, 1) == doctest::Approx(0.037327800344898).epsilon(1e-8));
    CHECK(result.L(1, 0) == doctest::Approx(0.186639001724491).epsilon(1e-8));
    CHECK(result.L(1, 1) == doctest::Approx(0.343422566088615).epsilon(1e-8));

    // compile-time version must agree
    constexpr auto dresult = [&]() {
        StateSpace<2, 1, 2, double, 2, 2> s{};
        s.A = {{1.0, dt}, {0.0, 1.0}};
        s.B = {{0.5 * dt * dt}, {dt}};
        s.C = Matrix<2, 2, double>::identity();
        s.G = Matrix<2, 2, double>::identity();
        s.H = Matrix<2, 2, double>::identity();
        s.Ts = dt;
        Matrix<2, 2> Qp = {{1e-4, 0.0}, {0.0, 1e-2}};
        Matrix<2, 2> Rm = {{1e-2, 0.0}, {0.0, 5e-2}};
        return design::kalman(s, Qp, Rm);
    }();
    static_assert(dresult.success);
    CHECK(dresult.P(0, 0) == doctest::Approx(0.002606406390998).epsilon(1e-8));
    CHECK(dresult.L(1, 1) == doctest::Approx(0.343422566088615).epsilon(1e-8));
}

TEST_CASE("Kalman design R!=0, non-square C (NY=1, NX=2)") {
    // Same double integrator, only position is measured (velocity is hidden state).
    // Reference values generated with scipy.linalg.solve_discrete_are.
    constexpr double                  dt = 0.1;
    StateSpace<2, 1, 1, double, 2, 1> sys{};
    sys.A = {{1.0, dt}, {0.0, 1.0}};
    sys.B = {{0.5 * dt * dt}, {dt}};
    sys.C = {{1.0, 0.0}}; // 1x2: observe position only
    sys.G = Matrix<2, 2, double>::identity();
    sys.H = Matrix<1, 1, double>::identity();
    sys.Ts = dt;

    Matrix<2, 2> Q_proc = {{1e-4, 0.0}, {0.0, 1e-2}};
    Matrix<1, 1> R_meas = {{5e-2}};

    auto result = design::kalman(sys, Q_proc, R_meas);
    REQUIRE(result.success);

    // scipy: solve_discrete_are(A.T, C.T, Q, R) with C = [[1,0]], R = [[0.05]]
    CHECK(result.P(0, 0) == doctest::Approx(0.017691003582802).epsilon(1e-8));
    CHECK(result.P(0, 1) == doctest::Approx(0.026017494803075).epsilon(1e-8));
    CHECK(result.P(1, 0) == doctest::Approx(0.026017494803075).epsilon(1e-8));
    CHECK(result.P(1, 1) == doctest::Approx(0.077996568142721).epsilon(1e-8));

    // L is 2x1: Kalman gain for position measurement
    CHECK(result.L(0, 0) == doctest::Approx(0.261349406072277).epsilon(1e-8));
    CHECK(result.L(1, 0) == doctest::Approx(0.384356759776049).epsilon(1e-8));

    // compile-time version must agree
    constexpr auto dresult = [&]() {
        StateSpace<2, 1, 1, double, 2, 1> s{};
        s.A = {{1.0, dt}, {0.0, 1.0}};
        s.B = {{0.5 * dt * dt}, {dt}};
        s.C = {{1.0, 0.0}};
        s.G = Matrix<2, 2, double>::identity();
        s.H = Matrix<1, 1, double>::identity();
        s.Ts = dt;
        Matrix<2, 2> Qp = {{1e-4, 0.0}, {0.0, 1e-2}};
        Matrix<1, 1> Rm = {{5e-2}};
        return design::kalman(s, Qp, Rm);
    }();
    static_assert(dresult.success);
    CHECK(dresult.P(0, 0) == doctest::Approx(0.017691003582802).epsilon(1e-8));
    CHECK(dresult.L(1, 0) == doctest::Approx(0.384356759776049).epsilon(1e-8));
}

TEST_CASE("Kalman design with R=0 (perfect measurements)") {
    // With R=0 (noiseless measurements) the standard DARE cannot be used (R
    // must be positive definite).  kalman() detects this, and for a square
    // invertible C solves analytically: P_ss = 0, L = C^{-1}.

    StateSpace<2, 1, 2, double, 2, 2> sys{};
    sys.A = {{0.9, 0.1}, {0.0, 0.8}};
    sys.B = {{0.0}, {1.0}};
    sys.C = Matrix<2, 2, double>::identity();
    sys.G = Matrix<2, 2, double>::identity();
    sys.H = Matrix<2, 2, double>::identity();
    sys.Ts = 0.01;

    Matrix<2, 2> Q_proc = {{0.01, 0.0}, {0.0, 0.01}};
    Matrix<2, 2> R_zero = Matrix<2, 2, double>::zeros();

    // kalman() should now succeed for the R=0 + square-C case
    auto result = design::kalman(sys, Q_proc, R_zero);
    CHECK(result.success);

    // P_ss = Q_eff = G*Q*G' = Q_proc (G=I here).
    // Scipy confirms: as R->0, solve_discrete_are(A', C', Q, R) -> Q.
    // Intuition: each predict step injects Q; the L=C^{-1}=I update fully
    // removes it, so the steady-state *prior* covariance equals Q.
    CHECK(result.P(0, 0) == doctest::Approx(Q_proc(0, 0)).epsilon(1e-12));
    CHECK(result.P(1, 1) == doctest::Approx(Q_proc(1, 1)).epsilon(1e-12));
    CHECK(result.P(0, 1) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(result.P(1, 0) == doctest::Approx(0.0).epsilon(1e-12));

    // L = C^{-1} = I (since C = I here): trust the measurement completely
    CHECK(result.L(0, 0) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(result.L(1, 1) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(result.L(0, 1) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(result.L(1, 0) == doctest::Approx(0.0).epsilon(1e-12));

    // Non-square C with R=0 now works via RDE fallback
    StateSpace<2, 1, 1, double, 2, 1> sys1d{};
    sys1d.A = {{0.9, 0.1}, {0.0, 0.8}};
    sys1d.B = {{0.0}, {1.0}};
    sys1d.C = {{1.0, 0.0}};
    sys1d.G = Matrix<2, 2, double>::identity();
    sys1d.H = Matrix<1, 1, double>::identity();
    sys1d.Ts = 0.01;

    Matrix<2, 2> Q2 = {{0.01, 0.0}, {0.0, 0.01}};
    Matrix<1, 1> R1_zero = Matrix<1, 1, double>::zeros();
    auto         result1d = design::kalman(sys1d, Q2, R1_zero);
    CHECK(result1d.success);
}

TEST_CASE("Kalman design R!=0, non-square C (NY=2, NX=4)") {
    // 4-state integrator chain, 2 measurements observing states 0 and 2.
    // Typical of a 2-axis inertial system where position is measured but not velocity.
    // A = chain-of-integrators with dt=0.1, C picks out states 0 and 2.
    // Reference values generated with scipy.linalg.solve_discrete_are.
    StateSpace<4, 1, 2, double, 4, 2> sys{};

    sys.A = {
        {1.0, 0.1, 0.0, 0.0},
        {0.0, 1.0, 0.1, 0.0},
        {0.0, 0.0, 1.0, 0.1},
        {0.0, 0.0, 0.0, 1.0},
    };

    // Single input on the last state
    sys.B = {
        {0.0},
        {0.0},
        {0.0},
        {1.0},
    };

    // 2x4: observe states 0 and 2
    sys.C = {
        {1.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0},
    };

    sys.G = Matrix<4, 4, double>::identity();
    sys.H = Matrix<2, 2, double>::identity();
    sys.Ts = 0.1;

    Matrix<4, 4> Q_proc = Matrix<4, 4, double>::identity();
    for (size_t i = 0; i < 4; ++i) {
        Q_proc(i, i) = 0.01;
    }

    Matrix<2, 2> R_meas = Matrix<2, 2, double>::identity();
    for (size_t i = 0; i < 2; ++i) {
        R_meas(i, i) = 0.1;
    }

    auto result = design::kalman(sys, Q_proc, R_meas);
    REQUIRE(result.success);

    // P is 4x4 symmetric
    CHECK(result.P(0, 0) == doctest::Approx(0.050147732129045).epsilon(1e-8));
    CHECK(result.P(0, 1) == doctest::Approx(0.040639610088462).epsilon(1e-8));
    CHECK(result.P(0, 2) == doctest::Approx(0.001193670164291).epsilon(1e-8));
    CHECK(result.P(0, 3) == doctest::Approx(-0.000207024981456).epsilon(1e-8));
    CHECK(result.P(1, 1) == doctest::Approx(0.146757944186451).epsilon(1e-8));
    CHECK(result.P(1, 2) == doctest::Approx(0.010539815102521).epsilon(1e-8));
    CHECK(result.P(1, 3) == doctest::Approx(0.003949387027125).epsilon(1e-8));
    CHECK(result.P(2, 2) == doctest::Approx(0.049606917603299).epsilon(1e-8));
    CHECK(result.P(2, 3) == doctest::Approx(0.038675628845935).epsilon(1e-8));
    CHECK(result.P(3, 3) == doctest::Approx(0.138213910986039).epsilon(1e-8));

    // L is 4x2: Kalman gain
    CHECK(result.L(0, 0) == doctest::Approx(0.333947026467608).epsilon(1e-8));
    CHECK(result.L(0, 1) == doctest::Approx(0.005314243318957).epsilon(1e-8));
    CHECK(result.L(1, 0) == doctest::Approx(0.270121220006257).epsilon(1e-8));
    CHECK(result.L(1, 1) == doctest::Approx(0.068294833054111).epsilon(1e-8));
    CHECK(result.L(2, 0) == doctest::Approx(0.005314243318957).epsilon(1e-8));
    CHECK(result.L(2, 1) == doctest::Approx(0.331539309439724).epsilon(1e-8));
    CHECK(result.L(3, 0) == doctest::Approx(-0.003434213066953).epsilon(1e-8));
    CHECK(result.L(3, 1) == doctest::Approx(0.258542377473313).epsilon(1e-8));

    // compile-time version must agree
    constexpr auto dresult = []() {
        StateSpace<4, 1, 2, double, 4, 2> s{};
        s.A = {{1.0, 0.1, 0.0, 0.0}, {0.0, 1.0, 0.1, 0.0}, {0.0, 0.0, 1.0, 0.1}, {0.0, 0.0, 0.0, 1.0}};
        s.B = {{0.0}, {0.0}, {0.0}, {1.0}};
        s.C = {{1.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}};
        s.G = Matrix<4, 4, double>::identity();
        s.H = Matrix<2, 2, double>::identity();
        s.Ts = 0.1;
        Matrix<4, 4> Qp = Matrix<4, 4, double>::identity();
        for (size_t i = 0; i < 4; ++i) {
            Qp(i, i) = 0.01;
        }
        Matrix<2, 2> Rm = Matrix<2, 2, double>::identity();
        for (size_t i = 0; i < 2; ++i) {
            Rm(i, i) = 0.1;
        }
        return design::kalman(s, Qp, Rm);
    }();

    static_assert(dresult.success);
    CHECK(dresult.P(0, 0) == doctest::Approx(0.050147732129045).epsilon(1e-8));
    CHECK(dresult.L(3, 1) == doctest::Approx(0.258542377473313).epsilon(1e-8));
}

TEST_CASE("Kalman design R=0, non-square C (NY=2, NX=4)") {
    // scipy reference (solve_discrete_are(A.T, C.T, Q, zeros(2,2))):
    StateSpace<4, 1, 2, double, 4, 2> sys{};

    sys.A = {
        {1.0, 0.1, 0.0, 0.0},
        {0.0, 1.0, 0.1, 0.0},
        {0.0, 0.0, 1.0, 0.1},
        {0.0, 0.0, 0.0, 1.0},
    };

    sys.B = {
        {0.0},
        {0.0},
        {0.0},
        {1.0},
    };

    sys.C = {
        {1.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0},
    };

    sys.G = Matrix<4, 4, double>::identity();
    sys.H = Matrix<2, 2, double>::identity();
    sys.Ts = 0.1;

    Matrix<4, 4> Q_proc = Matrix<4, 4, double>::identity();
    for (size_t i = 0; i < 4; ++i) {
        Q_proc(i, i) = 0.01;
    }

    Matrix<2, 2> R_zero = Matrix<2, 2, double>::zeros();

    auto result = design::kalman(sys, Q_proc, R_zero);
    REQUIRE(result.success);
    CHECK(result.P(0, 0) == doctest::Approx(1.105124921972504e-02).epsilon(1e-8));
    CHECK(result.P(0, 1) == doctest::Approx(1.051249219725040e-02).epsilon(1e-8));
    CHECK(result.P(0, 2) == doctest::Approx(0.0).epsilon(1e-8));
    CHECK(result.P(0, 3) == doctest::Approx(0.0).epsilon(1e-8));
    CHECK(result.P(1, 1) == doctest::Approx(1.151249219725040e-01).epsilon(1e-8));
    CHECK(result.P(1, 2) == doctest::Approx(0.0).epsilon(1e-8));
    CHECK(result.P(1, 3) == doctest::Approx(0.0).epsilon(1e-8));
    CHECK(result.P(2, 2) == doctest::Approx(1.105124921972504e-02).epsilon(1e-8));
    CHECK(result.P(2, 3) == doctest::Approx(1.051249219725037e-02).epsilon(1e-8));
    CHECK(result.P(3, 3) == doctest::Approx(1.151249219725037e-01).epsilon(1e-8));
    CHECK(result.L(0, 0) == doctest::Approx(1.0).epsilon(1e-8));
    CHECK(result.L(0, 1) == doctest::Approx(0.0).epsilon(1e-8));
    CHECK(result.L(1, 0) == doctest::Approx(9.512492197250400e-01).epsilon(1e-8));
    CHECK(result.L(1, 1) == doctest::Approx(0.0).epsilon(1e-8));
    CHECK(result.L(2, 0) == doctest::Approx(0.0).epsilon(1e-8));
    CHECK(result.L(2, 1) == doctest::Approx(1.0).epsilon(1e-8));
    CHECK(result.L(3, 0) == doctest::Approx(0.0).epsilon(1e-8));
    CHECK(result.L(3, 1) == doctest::Approx(9.512492197250375e-01).epsilon(1e-8));
}

TEST_CASE("Kalman design R=0 (NY=1, NX=4)") {
    StateSpace<4, 1, 1, double, 4, 1> sys{};

    sys.A = {
        {1.0000000000000000e+00, 1.2500000000000003e-04, 0.0, 0.0},
        {0.0, 1.0000000000000000e+00, 0.0, 0.0},
        {1.2044286656514691e-01, 7.6886479655468900e-06, 8.7955713343485309e-01, 0.0},
        {7.5633576541031576e-03, 3.2193761045041167e-07, 1.1287950891104376e-01, 8.7955713343485309e-01},
    };

    sys.B = {
        {7.8125000000000030e-09},
        {1.2500000000000003e-04},
        {3.2375688155733029e-10},
        {1.0189648978629236e-11},
    };

    sys.C = {{0.0, 0.0, 2.0, -1.0}};

    sys.G = Matrix<4, 4, double>::identity();
    sys.H = Matrix<1, 1, double>::identity();
    sys.Ts = 1.25e-4;

    // Qd from Van Loan discretization of continuous-time noise model
    Matrix<4, 4> Qd = {
        {6.5104166666666704e-13, 7.8125000000000046e-09, 3.0287979498172528e-14, 1.0167935310197761e-15},
        {7.8125000000000030e-09, 1.2500000000000003e-04, 3.2375688155733029e-10, 1.0189648978629238e-11},
        {3.0287979498172516e-14, 3.2375688155733013e-10, 2.7948205317517660e-10, 1.7167518071887303e-11},
        {1.0167935310197763e-15, 1.0189648978629233e-11, 1.7167518071887290e-11, 1.4368706442600565e-12},
    };

    Matrix<1, 1> R_zero = Matrix<1, 1, double>::zeros();

    auto result = design::kalman(sys, Qd, R_zero);
    REQUIRE(result.success);

    // Kalman gain (4x1): scipy reference
    CHECK(result.L(0, 0) == doctest::Approx(5.860848699642386e-01).epsilon(1e-4));
    CHECK(result.L(1, 0) == doctest::Approx(3.213133860108700e+02).epsilon(1e-4));
    CHECK(result.L(2, 0) == doctest::Approx(5.158533783244187e-01).epsilon(1e-4));
    CHECK(result.L(3, 0) == doctest::Approx(3.170675664883725e-02).epsilon(1e-4));

    // Covariance P (spot checks on diagonal and key off-diagonals)
    CHECK(result.P(0, 0) == doctest::Approx(3.1493750780952886e-09).epsilon(1e-4));
    CHECK(result.P(1, 1) == doctest::Approx(1.8865288001709108e-03).epsilon(1e-4));
    CHECK(result.P(2, 2) == doctest::Approx(3.2303061178439298e-10).epsilon(1e-4));
    CHECK(result.P(3, 3) == doctest::Approx(4.6007682900265547e-12).epsilon(1e-4));
    CHECK(result.P(0, 1) == doctest::Approx(1.7788444152705028e-06).epsilon(1e-4));
}

TEST_CASE("Kalman set_state clamps the estimate to a physical bound") {
    // 1D filter; the state represents a physically non-negative quantity (e.g.
    // a concentration / SoC). A measurement pulls the estimate negative; the
    // caller clamps it back to 0 via set_state, and the next predict/update
    // proceeds from the constrained value rather than the non-physical one.
    StateSpace<1, 1, 1, double, 1, 1> sys{};
    sys.A(0, 0) = 1.0;
    sys.C(0, 0) = 1.0;
    sys.G(0, 0) = 1.0;
    sys.H(0, 0) = 1.0;
    sys.Ts = 1.0;

    KalmanFilter<1, 1, 1, double, 1, 1> kf(sys, Matrix<1, 1>{{0.1}}, Matrix<1, 1>{{0.25}}, ColVec<1>{0.0}, Matrix<1, 1>::identity());

    kf.predict();
    REQUIRE(kf.update(ColVec<1>{-5.0})); // measurement drives estimate negative
    REQUIRE(kf.state()[0] < 0.0);

    if (kf.state()[0] < 0.0) {
        kf.set_state(0, 0.0); // clamp to physical floor
    }
    CHECK(kf.state()[0] == doctest::Approx(0.0));

    // Covariance setter likewise lets the caller re-inflate uncertainty after a
    // hard reset of the estimate.
    kf.set_covariance(Matrix<1, 1>{{2.0}});
    CHECK(kf.covariance()(0, 0) == doctest::Approx(2.0));

    // Vector setter overload.
    kf.set_state(ColVec<1>{0.7});
    CHECK(kf.state()[0] == doctest::Approx(0.7));
}

TEST_CASE("design::kalman rejects continuous plant (Ts=0)") {
    StateSpace<1, 1, 1, double, 1, 1> sys{};
    sys.A(0, 0) = 0.0;
    sys.B(0, 0) = 1.0;
    sys.C(0, 0) = 1.0;
    sys.G(0, 0) = 1.0;
    sys.H(0, 0) = 1.0;
    sys.Ts = 0.0; // continuous
    auto result = design::kalman(sys, Matrix<1, 1>{{0.1}}, Matrix<1, 1>{{0.1}});
    CHECK_FALSE(result.success);
}

TEST_CASE("Kalman Joseph update keeps P symmetric; singular S fails closed") {
    StateSpace<2, 1, 1, double, 2, 1> sys{};
    sys.A = {{1.0, 0.1}, {0.0, 1.0}};
    sys.B = {{0.005}, {0.1}};
    sys.C = {{1.0, 0.0}};
    sys.G = Matrix<2, 2, double>::identity();
    sys.H = Matrix<1, 1, double>::identity();
    sys.Ts = 0.1;

    KalmanFilter<2, 1, 1, double, 2, 1> kf(
        sys, Matrix<2, 2, double>::identity() * 0.01, Matrix<1, 1, double>{{0.1}}, ColVec<2>{0.0, 0.0},
        Matrix<2, 2, double>::identity()
    );

    kf.predict(ColVec<1>{0.0});
    REQUIRE(kf.update(ColVec<1>{1.0}, ColVec<1>{0.0}));
    const auto& P = kf.covariance();
    CHECK(P(0, 1) == doctest::Approx(P(1, 0)).epsilon(1e-14));
    CHECK(P(0, 0) > 0.0);
    CHECK(P(1, 1) > 0.0);

    // Zero R and zero prior on the measured channel → S singular → update rejects.
    KalmanFilter<2, 1, 1, double, 2, 1> kf_sing(
        sys, Matrix<2, 2, double>::zeros(), Matrix<1, 1, double>{{0.0}}, ColVec<2>{}, Matrix<2, 2, double>::zeros()
    );
    kf_sing.predict(ColVec<1>{0.0});
    CHECK_FALSE(kf_sing.update(ColVec<1>{1.0}, ColVec<1>{0.0}));
}

TEST_CASE("design::kalman R≈0 square C uses LU solve path (L = C^{-1})") {
    StateSpace<2, 1, 2, double, 2, 2> sys{};
    sys.A = {{0.9, 0.1}, {0.0, 0.8}};
    sys.B = {{0.0}, {1.0}};
    sys.C = Matrix<2, 2, double>::identity();
    sys.G = Matrix<2, 2, double>::identity();
    sys.H = Matrix<2, 2, double>::identity();
    sys.Ts = 0.01;
    // Near-zero R → analytical L = C^{-1} = I
    const auto result = design::kalman(sys, Matrix<2, 2, double>::identity() * 0.01, Matrix<2, 2, double>::zeros());
    REQUIRE(result.success);
    CHECK(result.L(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(result.L(1, 1) == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(result.P(0, 1) == doctest::Approx(result.P(1, 0)).epsilon(1e-14));
}

TEST_CASE("design::discrete_lqe matches design::kalman (sys form)") {
    // Dual-of-LQR spelling is a pure alias — bit-identical L and P.
    StateSpace<2, 1, 1, double, 2, 1> sys{};
    sys.A = Matrix<2, 2>{{0.9, 0.1}, {0.0, 0.8}};
    sys.B = Matrix<2, 1>{{0.0}, {1.0}};
    sys.C = Matrix<1, 2>{{1.0, 0.0}};
    sys.G = Matrix<2, 2>::identity();
    sys.H = Matrix<1, 1>{{1.0}};
    sys.Ts = 0.01;

    const Matrix<2, 2> Q = Matrix<2, 2>::diagonal({0.1, 0.1});
    const Matrix<1, 1> R{{0.5}};

    const auto via_kalman = design::kalman(sys, Q, R);
    const auto via_lqe = design::discrete_lqe(sys, Q, R);
    REQUIRE(via_kalman.success);
    REQUIRE(via_lqe.success);
    CHECK(via_lqe.L(0, 0) == doctest::Approx(via_kalman.L(0, 0)));
    CHECK(via_lqe.L(1, 0) == doctest::Approx(via_kalman.L(1, 0)));
    CHECK(via_lqe.P(0, 0) == doctest::Approx(via_kalman.P(0, 0)));
    CHECK(via_lqe.P(1, 1) == doctest::Approx(via_kalman.P(1, 1)));
}

TEST_CASE("design::discrete_lqe matrix form matches sys form with G and H=I") {
    const Matrix<2, 2> A{{0.9, 0.1}, {0.0, 0.8}};
    const Matrix<2, 2> G = Matrix<2, 2>::identity();
    const Matrix<1, 2> C{{1.0, 0.0}};
    const Matrix<2, 2> Q = Matrix<2, 2>::diagonal({0.1, 0.1});
    const Matrix<1, 1> R{{0.5}};

    StateSpace<2, 1, 1, double, 2, 1> sys{};
    sys.A = A;
    sys.G = G;
    sys.C = C;
    sys.H = Matrix<1, 1>{{1.0}};
    sys.Ts = 1.0;

    const auto via_sys = design::discrete_lqe(sys, Q, R);
    const auto via_mat = design::discrete_lqe(A, G, C, Q, R);
    REQUIRE(via_sys.success);
    REQUIRE(via_mat.success);
    CHECK(via_mat.L(0, 0) == doctest::Approx(via_sys.L(0, 0)));
    CHECK(via_mat.L(1, 0) == doctest::Approx(via_sys.L(1, 0)));
    CHECK(via_mat.P(0, 0) == doctest::Approx(via_sys.P(0, 0)));
    CHECK(via_mat.P(1, 1) == doctest::Approx(via_sys.P(1, 1)));
}

TEST_CASE("design::kalman matches scipy/control golden data") {
    constexpr Matrix<2, 2> A{{0.9, 0.1}, {0.0, 0.8}};
    constexpr Matrix<2, 1> B{{0.1}, {0.2}};
    constexpr Matrix<1, 2> C{{1.0, 0.0}};
    constexpr Matrix<2, 2> Q{{0.1, 0.0}, {0.0, 0.1}};
    constexpr Matrix<1, 1> R{{0.5}};

    constexpr auto result = design::kalman(
        StateSpace<2, 1, 1, double, 2, 1>{
            A, B, C, Matrix<1, 1>::zeros(), Matrix<2, 2>::identity(), Matrix<1, 1>::identity(), 1.0
        },
        Q, R
    );
    static_assert(result.success);
    CHECK(doctest::Approx(result.L(0, 0)).epsilon(1e-3) == 0.323);
    CHECK(doctest::Approx(result.L(1, 0)).epsilon(1e-3) == 0.057);
    CHECK(result.P(0, 0) > 0.0);
    CHECK(result.P(1, 1) > 0.0);
}

TEST_CASE("KalmanResult::as<U>() conversion") {
    constexpr auto kalman_d = design::kalman(
        StateSpace<1, 1, 1, double, 1, 1>{
            Matrix<1, 1>{{1.0}},
            Matrix<1, 1>{{1.0}},
            Matrix<1, 1>{{1.0}},
            Matrix<1, 1>::zeros(),
            Matrix<1, 1>::identity(),
            Matrix<1, 1>::identity(),
            1.0
        },
        Matrix<1, 1>{{0.1}},
        Matrix<1, 1>{{0.5}}
    );
    constexpr auto kalman_f = kalman_d.as<float>();
    static_assert(kalman_f.success);
    CHECK(kalman_f.success);
    CHECK(kalman_f.L(0, 0) != 0.0f);
}
