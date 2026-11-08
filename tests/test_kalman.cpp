// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>

#include "damp/design/status.hpp"
#include "damp/estimation/fixed_gain_observer.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
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

    KalmanFilter kf{sys, {.Q = 0.1, .R = 0.25}};

    ColVec<1> u = {0.0};
    kf.predict(u);

    ColVec<1> y = {1.2};
    bool      ok = kf.update(y);
    CHECK(ok);
    CHECK(kf.state()[0] == doctest::Approx(0.978).epsilon(1e-3));
    CHECK(kf.covariance()(0, 0) == doctest::Approx(0.204).epsilon(1e-3));
}

TEST_CASE("KalmanFilter from KalmanResult + fused step()") {
    StateSpace<1, 1, 1, double, 1, 1> sys{};
    sys.A(0, 0) = 0.95;
    sys.B(0, 0) = 1.0;
    sys.C(0, 0) = 1.0;
    sys.G(0, 0) = 1.0;
    sys.H(0, 0) = 1.0;
    sys.Ts = 1.0;

    const auto result = design::kalman(sys, {.process = 0.1, .measurement = 0.25});
    REQUIRE(result.success);

    KalmanFilter<1, 1, 1, double, 1, 1> kf{result};
    CHECK(kf.covariance()(0, 0) == doctest::Approx(result.P(0, 0)));
    CHECK(kf.model().A(0, 0) == doctest::Approx(0.95));

    const ColVec<1> y{1.0};
    const ColVec<1> u{0.0};
    const auto&     x = kf.step(y, u);
    CHECK(std::isfinite(x[0]));
    // After a measurement of 1 from a near-zero seed, estimate moves toward y
    CHECK(x[0] > 0.0);
    CHECK(x[0] < 1.5);
}

TEST_CASE("FixedGainObserver runs design::kalman's fixed gain M; set_gain and reset") {
    StateSpace<1, 1, 1, double, 1, 1> sys{};
    sys.A(0, 0) = 0.9;
    sys.B(0, 0) = 1.0;
    sys.C(0, 0) = 1.0;
    sys.G(0, 0) = 1.0;
    sys.H(0, 0) = 1.0;
    sys.Ts = 1.0;

    const auto result = design::kalman(sys, {.process = 0.1, .measurement = 0.1});
    REQUIRE(result.success);

    FixedGainObserver sskf{result};
    static_assert(std::is_same_v<decltype(sskf), FixedGainObserver<1, 1, 1, double>>);
    CHECK(sskf.valid());
    CHECK(sskf.gain()(0, 0) == doctest::Approx(result.M(0, 0)));

    const auto& x = sskf.step(ColVec<1>{2.0}, ColVec<1>{0.0});
    CHECK(x[0] != 0.0);

    sskf.set_gain(Matrix<1, 1>{{0.0}}); // open-loop (no correction)
    sskf.reset(ColVec<1>{0.0});
    sskf.predict(ColVec<1>{0.0});
    sskf.update(ColVec<1>{5.0});
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

    auto result = design::kalman(sys, {.process = Q_proc, .measurement = R_meas});
    REQUIRE(result.success);

    // scipy: solve_discrete_are(A.T, C.T, Q, R)
    CHECK(result.P(0, 0) == doctest::Approx(0.002606406390998).epsilon(1e-8));
    CHECK(result.P(0, 1) == doctest::Approx(0.003583502847688).epsilon(1e-8));
    CHECK(result.P(1, 0) == doctest::Approx(0.003583502847688).epsilon(1e-8));
    CHECK(result.P(1, 1) == doctest::Approx(0.027171128304431).epsilon(1e-8));

    CHECK(result.M(0, 0) == doctest::Approx(0.196141710450478).epsilon(1e-8));
    CHECK(result.M(0, 1) == doctest::Approx(0.037327800344898).epsilon(1e-8));
    CHECK(result.M(1, 0) == doctest::Approx(0.186639001724491).epsilon(1e-8));
    CHECK(result.M(1, 1) == doctest::Approx(0.343422566088615).epsilon(1e-8));

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
        return design::kalman(s, {.process = Qp, .measurement = Rm});
    }();
    static_assert(dresult.success);
    CHECK(dresult.P(0, 0) == doctest::Approx(0.002606406390998).epsilon(1e-8));
    CHECK(dresult.M(1, 1) == doctest::Approx(0.343422566088615).epsilon(1e-8));
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

    auto result = design::kalman(sys, {.process = Q_proc, .measurement = R_meas});
    REQUIRE(result.success);

    // scipy: solve_discrete_are(A.T, C.T, Q, R) with C = [[1,0]], R = [[0.05]]
    CHECK(result.P(0, 0) == doctest::Approx(0.017691003582802).epsilon(1e-8));
    CHECK(result.P(0, 1) == doctest::Approx(0.026017494803075).epsilon(1e-8));
    CHECK(result.P(1, 0) == doctest::Approx(0.026017494803075).epsilon(1e-8));
    CHECK(result.P(1, 1) == doctest::Approx(0.077996568142721).epsilon(1e-8));

    // L is 2x1: Kalman gain for position measurement
    CHECK(result.M(0, 0) == doctest::Approx(0.261349406072277).epsilon(1e-8));
    CHECK(result.M(1, 0) == doctest::Approx(0.384356759776049).epsilon(1e-8));

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
        return design::kalman(s, {.process = Qp, .measurement = Rm});
    }();
    static_assert(dresult.success);
    CHECK(dresult.P(0, 0) == doctest::Approx(0.017691003582802).epsilon(1e-8));
    CHECK(dresult.M(1, 0) == doctest::Approx(0.384356759776049).epsilon(1e-8));
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
    auto result = design::kalman(sys, {.process = Q_proc, .measurement = R_zero});
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
    CHECK(result.M(0, 0) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(result.M(1, 1) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(result.M(0, 1) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(result.M(1, 0) == doctest::Approx(0.0).epsilon(1e-12));

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
    auto         result1d = design::kalman(sys1d, {.process = Q2, .measurement = R1_zero});
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

    auto result = design::kalman(sys, {.process = Q_proc, .measurement = R_meas});
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
    CHECK(result.M(0, 0) == doctest::Approx(0.333947026467608).epsilon(1e-8));
    CHECK(result.M(0, 1) == doctest::Approx(0.005314243318957).epsilon(1e-8));
    CHECK(result.M(1, 0) == doctest::Approx(0.270121220006257).epsilon(1e-8));
    CHECK(result.M(1, 1) == doctest::Approx(0.068294833054111).epsilon(1e-8));
    CHECK(result.M(2, 0) == doctest::Approx(0.005314243318957).epsilon(1e-8));
    CHECK(result.M(2, 1) == doctest::Approx(0.331539309439724).epsilon(1e-8));
    CHECK(result.M(3, 0) == doctest::Approx(-0.003434213066953).epsilon(1e-8));
    CHECK(result.M(3, 1) == doctest::Approx(0.258542377473313).epsilon(1e-8));

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
        return design::kalman(s, {.process = Qp, .measurement = Rm});
    }();

    static_assert(dresult.success);
    CHECK(dresult.P(0, 0) == doctest::Approx(0.050147732129045).epsilon(1e-8));
    CHECK(dresult.M(3, 1) == doctest::Approx(0.258542377473313).epsilon(1e-8));
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

    auto result = design::kalman(sys, {.process = Q_proc, .measurement = R_zero});
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
    CHECK(result.M(0, 0) == doctest::Approx(1.0).epsilon(1e-8));
    CHECK(result.M(0, 1) == doctest::Approx(0.0).epsilon(1e-8));
    CHECK(result.M(1, 0) == doctest::Approx(9.512492197250400e-01).epsilon(1e-8));
    CHECK(result.M(1, 1) == doctest::Approx(0.0).epsilon(1e-8));
    CHECK(result.M(2, 0) == doctest::Approx(0.0).epsilon(1e-8));
    CHECK(result.M(2, 1) == doctest::Approx(1.0).epsilon(1e-8));
    CHECK(result.M(3, 0) == doctest::Approx(0.0).epsilon(1e-8));
    CHECK(result.M(3, 1) == doctest::Approx(9.512492197250375e-01).epsilon(1e-8));
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

    auto result = design::kalman(sys, {.process = Qd, .measurement = R_zero});
    REQUIRE(result.success);

    // Kalman gain (4x1): scipy reference
    CHECK(result.M(0, 0) == doctest::Approx(5.860848699642386e-01).epsilon(1e-4));
    CHECK(result.M(1, 0) == doctest::Approx(3.213133860108700e+02).epsilon(1e-4));
    CHECK(result.M(2, 0) == doctest::Approx(5.158533783244187e-01).epsilon(1e-4));
    CHECK(result.M(3, 0) == doctest::Approx(3.170675664883725e-02).epsilon(1e-4));

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
    // proceeds from the constrained value rather than the non-physical one;
    // reset re-initializes the estimate and re-inflates its covariance.
    StateSpace<1, 1, 1, double, 1, 1> sys{};
    sys.A(0, 0) = 1.0;
    sys.C(0, 0) = 1.0;
    sys.G(0, 0) = 1.0;
    sys.H(0, 0) = 1.0;
    sys.Ts = 1.0;

    KalmanFilter kf{sys, {.Q = 0.1, .R = 0.25}};

    kf.predict({});
    REQUIRE(kf.update(ColVec<1>{-5.0})); // measurement drives estimate negative
    REQUIRE(kf.state()[0] < 0.0);

    if (kf.state()[0] < 0.0) {
        kf.set_state(ColVec<1>{0.0}); // clamp to physical floor
    }
    CHECK(kf.state()[0] == doctest::Approx(0.0));

    kf.reset(ColVec<1>{0.7}, Matrix<1, 1>{{2.0}});
    CHECK(kf.state()[0] == doctest::Approx(0.7));
    CHECK(kf.covariance()(0, 0) == doctest::Approx(2.0));
}

TEST_CASE("design::kalman rejects continuous plant (Ts=0)") {
    StateSpace<1, 1, 1, double, 1, 1> sys{};
    sys.A(0, 0) = 0.0;
    sys.B(0, 0) = 1.0;
    sys.C(0, 0) = 1.0;
    sys.G(0, 0) = 1.0;
    sys.H(0, 0) = 1.0;
    sys.Ts = 0.0; // continuous
    auto result = design::kalman(sys, {.process = 0.1, .measurement = 0.1});
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

    KalmanFilter kf{sys, {.Q = Matrix<2, 2, double>::identity() * 0.01, .R = 0.1}};

    kf.predict(ColVec<1>{0.0});
    REQUIRE(kf.update(ColVec<1>{1.0}));
    const auto& P = kf.covariance();
    CHECK(P(0, 1) == doctest::Approx(P(1, 0)).epsilon(1e-14));
    CHECK(P(0, 0) > 0.0);
    CHECK(P(1, 1) > 0.0);

    // Zero R and zero prior on the measured channel → S singular → update rejects.
    KalmanFilter kf_sing{sys, {.Q = Matrix<2, 2, double>::zeros(), .R = 0.0, .P0 = Matrix<2, 2, double>::zeros()}};
    kf_sing.predict(ColVec<1>{0.0});
    CHECK_FALSE(kf_sing.update(ColVec<1>{1.0}));
}

TEST_CASE("KalmanFilter fuses a second sensor of another size: update(y, C, R) and update(y, C, D, R)") {
    // A = I, B = [0; 1], the model's own reading y = x₀ (C = [1 0], H = 1, R = 1), Q = 0,
    // x̂₀ = 0, P₀ = I. Every P stays diagonal, so each gain is P_ii / (P_ii + r) per channel.
    StateSpace<2, 1, 1, double, 2, 1> sys{};
    sys.A = Matrix<2, 2, double>::identity();
    sys.B = {{0.0}, {1.0}};
    sys.C = {{1.0, 0.0}};
    sys.G = Matrix<2, 2, double>::identity();
    sys.H = Matrix<1, 1, double>::identity();
    sys.Ts = 1.0;

    KalmanFilter kf{sys, {.Q = Matrix<2, 2, double>::zeros(), .R = 1.0}};
    static_assert(std::is_same_v<decltype(kf), KalmanFilter<2, 1, 1, double, 2, 1>>);
    kf.predict(ColVec<1>{0.0});

    // A 2-axis reading of the full state, y = I x + v, R = diag(1, 3): S = diag(2, 4),
    // K = diag(1/2, 1/4), x̂ = K (2, 4) = (1, 1), P = (I − K)² + K² R = diag(1/2, 3/4).
    REQUIRE(kf.update(ColVec<2>{2.0, 4.0}, Matrix<2, 2, double>::identity(), Matrix<2, 2, double>::diagonal({1.0, 3.0})));
    CHECK(kf.state()[0] == doctest::Approx(1.0).epsilon(1e-14));
    CHECK(kf.state()[1] == doctest::Approx(1.0).epsilon(1e-14));
    CHECK(kf.covariance()(0, 0) == doctest::Approx(0.5).epsilon(1e-14));
    CHECK(kf.covariance()(1, 1) == doctest::Approx(0.75).epsilon(1e-14));
    CHECK(kf.covariance()(0, 1) == doctest::Approx(0.0));

    // The model's own reading y = 1.5: S = 1/2 + 1, K = (1/3, 0), x̂₀ = 1 + (1.5 − 1)/3 = 7/6, P₀₀ = 1/3.
    REQUIRE(kf.update(ColVec<1>{1.5}));
    CHECK(kf.state()[0] == doctest::Approx(7.0 / 6.0).epsilon(1e-14));
    CHECK(kf.state()[1] == doctest::Approx(1.0).epsilon(1e-14));
    CHECK(kf.covariance()(0, 0) == doctest::Approx(1.0 / 3.0).epsilon(1e-14));

    // predict(u = 2) moves x̂₁ to 3 and leaves P. A scalar reading of x₁ with feedthrough on the
    // held input, y = x₁ + 0.5 u, R = 3/4: ŷ = 3 + 1 = 4, y = 4.5, S = 3/4 + 3/4, K = (0, 1/2),
    // x̂₁ = 3 + 0.5·0.5 = 3.25, P₁₁ = 3/4 − (3/4)·(1/2) = 3/8.
    kf.predict(ColVec<1>{2.0});
    REQUIRE(kf.update(ColVec<1>{4.5}, Matrix<1, 2, double>{{0.0, 1.0}}, Matrix<1, 1, double>{0.5}, 0.75));
    CHECK(kf.state()[0] == doctest::Approx(7.0 / 6.0).epsilon(1e-14));
    CHECK(kf.state()[1] == doctest::Approx(3.25).epsilon(1e-14));
    CHECK(kf.covariance()(0, 0) == doctest::Approx(1.0 / 3.0).epsilon(1e-14));
    CHECK(kf.covariance()(1, 1) == doctest::Approx(0.375).epsilon(1e-14));

    // A reading whose S = P + R is not positive definite is refused; the estimate and P stay.
    CHECK_FALSE(kf.update(ColVec<2>{0.0, 0.0}, Matrix<2, 2, double>::identity(), Matrix<2, 2, double>::identity() * -1.0));
    CHECK(kf.state()[1] == doctest::Approx(3.25).epsilon(1e-14));
    CHECK(kf.covariance()(1, 1) == doctest::Approx(0.375).epsilon(1e-14));
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
    const auto result = design::kalman(sys, {.process = Matrix<2, 2, double>::identity() * 0.01, .measurement = Matrix<2, 2, double>::zeros()});
    REQUIRE(result.success);
    CHECK(result.M(0, 0) == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(result.M(1, 1) == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(result.P(0, 1) == doctest::Approx(result.P(1, 0)).epsilon(1e-14));
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
        {.process = Q, .measurement = R}
    );
    static_assert(result.success);
    CHECK(doctest::Approx(result.M(0, 0)).epsilon(1e-3) == 0.323);
    CHECK(doctest::Approx(result.M(1, 0)).epsilon(1e-3) == 0.057);
    CHECK(result.P(0, 0) > 0.0);
    CHECK(result.P(1, 1) > 0.0);
}

TEST_CASE("KalmanResult carries MATLAB's two gains: the innovation gain M and the predictor gain L = A M") {
    //! The plain filter: M = P Cᵀ (C P Cᵀ + H R Hᵀ)⁻¹, L = A M (MATLAB®'s Mx and L)
    StateSpace<2, 1, 1, double, 2, 1> sys{};
    sys.A = {{0.9, 0.1}, {0.0, 0.8}};
    sys.B = {{0.1}, {0.2}};
    sys.C = {{1.0, 0.0}};
    sys.Ts = 1.0;
    constexpr Matrix<2, 2> Q{{0.1, 0.0}, {0.0, 0.1}};
    constexpr Matrix<1, 1> R{{0.5}};
    const auto             res = design::kalman(sys, {.process = Q, .measurement = R});
    REQUIRE(res.success);
    const Matrix<1, 1> S = (sys.C * res.P * sys.C.t()) + R;
    CHECK((res.M * S - (res.P * sys.C.t())).norm() < 1e-14);
    CHECK((res.L - (sys.A * res.M)).norm() == 0.0);

    //! The predictor's error A − L C and the current form's (I − M C) A share their poles.
    const auto   pred = mat::compute_eigenvalues<2, double>(sys.A - (res.L * sys.C)).values;
    const auto   curr = mat::compute_eigenvalues<2, double>((Matrix<2, 2>::identity() - (res.M * sys.C)) * sys.A).values;
    const double p0 = damp::min(pred[0].real(), pred[1].real());
    const double c0 = damp::min(curr[0].real(), curr[1].real());
    CHECK(p0 == doctest::Approx(c0).epsilon(1e-12));

    //! A design on [x; d]: both gains on the augmented model
    const auto aug = design::kalman(sys, {.process = Q, .measurement = R, .disturbance = Matrix<1, 1>{{1e-3}}}, design::InputDisturbance{});
    REQUIRE(aug.success);
    CHECK((aug.L - (aug.sys.A * aug.M)).norm() == 0.0);

    //! as<U>() carries both
    const auto res_f = res.as<float>();
    CHECK(res_f.L(0, 0) == doctest::Approx(static_cast<float>(res.L(0, 0))));
    CHECK(res_f.M(1, 0) == doctest::Approx(static_cast<float>(res.M(1, 0))));
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
        {.process = 0.1, .measurement = 0.5}
    );
    constexpr auto kalman_f = kalman_d.as<float>();
    static_assert(kalman_f.success);
    CHECK(kalman_f.success);
    CHECK(kalman_f.M(0, 0) != 0.0f);
}

TEST_CASE("design::kalman status: undetectable, continuous, ok") {
    //! The mode at 1.5 is unstable and the sensor does not see it: (A, C) is not detectable.
    StateSpace<2, 1, 1, double, 2, 1> sys{};
    sys.A = Matrix<2, 2>{{1.5, 0.0}, {0.0, 0.5}};
    sys.B = Matrix<2, 1>{{0.0}, {1.0}};
    sys.C = Matrix<1, 2>{{0.0, 1.0}};
    sys.G = Matrix<2, 2>::identity();
    sys.H = Matrix<1, 1>{{1.0}};
    sys.Ts = 0.1;
    const Matrix<2, 2> Q = Matrix<2, 2>::identity();
    const Matrix<1, 1> R{{1.0}};

    const auto undetectable = design::kalman(sys, {.process = Q, .measurement = R});
    CHECK_FALSE(undetectable.success);
    CHECK(undetectable.status == design::DesignStatus::not_detectable);

    //! No continuous-time Kalman design: a continuous plant is refused, not solved as if discrete.
    auto continuous = sys;
    continuous.Ts = 0.0;
    continuous.C = Matrix<1, 2>{{1.0, 1.0}};
    const auto refused = design::kalman(continuous, {.process = Q, .measurement = R});
    CHECK_FALSE(refused.success);
    CHECK(refused.status == design::DesignStatus::invalid_sample_time);

    auto detectable = sys;
    detectable.C = Matrix<1, 2>{{1.0, 1.0}};
    const auto ok = design::kalman(detectable, {.process = Q, .measurement = R});
    CHECK(ok.success);
    CHECK(ok.status == design::DesignStatus::ok);
    CHECK(ok.as<float>().status == design::DesignStatus::ok);

    //! A non-finite model entry is refused, also one the filter DARE never reads (B, D): the
    //! predictor x̂ = A x̂ + B u would run on it.
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (int entry = 0; entry < 3; ++entry) {
        auto bad = detectable;
        if (entry == 0) {
            bad.B(1, 0) = nan;
        } else if (entry == 1) {
            bad.C(0, 0) = nan;
        } else {
            bad.D(0, 0) = nan;
        }
        const auto plain = design::kalman(bad, {.process = Q, .measurement = R});
        CHECK_FALSE(plain.success);
        CHECK(plain.status == design::DesignStatus::bad_weights);
        const auto with_disturbance = design::kalman(bad, {.process = Q, .measurement = R, .disturbance = Matrix<1, 1>{{1e-3}}}, design::InputDisturbance{});
        CHECK(with_disturbance.status == design::DesignStatus::bad_weights);
    }
}

TEST_CASE("design::kalman refuses process noise that leaves a unit-circle mode unexcited") {
    //! A constant bias (z = 1) with no noise on it: no filter makes its error decay.
    StateSpace<2, 1, 1, double, 1, 1> sys{};
    sys.A = Matrix<2, 2>{{0.9, 0.0}, {0.0, 1.0}};
    sys.B = Matrix<2, 1>{{1.0}, {0.0}};
    sys.C = Matrix<1, 2>{{1.0, 1.0}};
    sys.G = Matrix<2, 1>{{1.0}, {0.0}};
    sys.H = Matrix<1, 1>{{1.0}};
    sys.Ts = 0.1;
    const auto bias = design::kalman(sys, {.process = 1.0, .measurement = 1.0});
    CHECK_FALSE(bias.success);
    CHECK(bias.status == design::DesignStatus::bad_weights);

    //! The same for a disturbance state with no walk; any positive walk designs.
    StateSpace<1, 1, 1, double, 1, 1> plant{};
    plant.A = Matrix<1, 1>{{0.9}};
    plant.B = Matrix<1, 1>{{1.0}};
    plant.C = Matrix<1, 1>{{1.0}};
    plant.G = Matrix<1, 1>{{1.0}};
    plant.H = Matrix<1, 1>{{1.0}};
    plant.Ts = 0.1;
    const auto no_walk = design::kalman(plant, {.process = 1.0, .measurement = 1.0, .disturbance = 0.0}, design::InputDisturbance{});
    CHECK(no_walk.status == design::DesignStatus::bad_weights);
    CHECK(design::kalman(plant, {.process = 1.0, .measurement = 1.0, .disturbance = 1e-12}, design::InputDisturbance{}).success);
}

namespace {
// Continuous damped oscillator, position sensor, process noise on the velocity.
constexpr StateSpace<2, 1, 1, double, 1, 1> make_continuous_oscillator() {
    return StateSpace<2, 1, 1, double, 1, 1>{
        .A = Matrix<2, 2>{{0.0, 1.0}, {-2.0, -0.5}},
        .B = Matrix<2, 1>{{0.0}, {1.0}},
        .C = Matrix<1, 2>{{1.0, 0.0}},
        .G = Matrix<2, 1>{{0.0}, {1.0}},
        .H = Matrix<1, 1>{{1.0}},
    };
}
} // namespace

TEST_CASE("discretize_process_noise: Van Loan against SciPy, the double integrator, and Ts → 0") {
    constexpr auto sys = make_continuous_oscillator();
    constexpr auto Qd = design::discretize_process_noise(sys.A, sys.G, Matrix<1, 1>{{0.3}}, 0.1).value();
    static_assert(Qd(0, 1) == Qd(1, 0)); // symmetrized

    //! scipy.linalg.expm on [−A GQGᵀ; 0 Aᵀ]·Ts
    CHECK(Qd(0, 0) == doctest::Approx(9.5952957935106529e-05).epsilon(1e-12));
    CHECK(Qd(0, 1) == doctest::Approx(1.4176528742916635e-03).epsilon(1e-12));
    CHECK(Qd(1, 1) == doctest::Approx(2.8361612946555551e-02).epsilon(1e-12));

    //! Double integrator: Qd = q·[Ts³/3 Ts²/2; Ts²/2 Ts] exactly.
    const Matrix<2, 2> A_di{{0.0, 1.0}, {0.0, 0.0}};
    const Matrix<2, 1> G_di{{0.0}, {1.0}};
    const double       Ts = 0.05;
    const double       q = 4.0;
    const auto         Qd_di = design::discretize_process_noise(A_di, G_di, Matrix<1, 1>{{q}}, Ts).value();
    CHECK(Qd_di(0, 0) == doctest::Approx(q * Ts * Ts * Ts / 3.0).epsilon(1e-12));
    CHECK(Qd_di(0, 1) == doctest::Approx(q * Ts * Ts / 2.0).epsilon(1e-12));
    CHECK(Qd_di(1, 1) == doctest::Approx(q * Ts).epsilon(1e-12));

    //! Short period: Qd → G Qc Gᵀ Ts, the error first order in Ts.
    const double Ts_small = 1e-5;
    const auto   Qd_small = design::discretize_process_noise(sys.A, sys.G, Matrix<1, 1>{{0.3}}, Ts_small).value();
    CHECK(Qd_small(1, 1) == doctest::Approx(0.3 * Ts_small).epsilon(1e-4));
    CHECK(damp::abs(Qd_small(0, 0)) < 1e-4 * 0.3 * Ts_small);

    //! Float: the same construction at single precision.
    constexpr auto Qd_f = design::discretize_process_noise(sys.A.as<float>(), sys.G.as<float>(), Matrix<1, 1, float>{{0.3f}}, 0.1f).value();
    CHECK(Qd_f(1, 1) == doctest::Approx(2.8361612946555551e-02).epsilon(1e-5));
}

TEST_CASE_TEMPLATE("discretize_process_noise keeps a fast stable mode", T, float, double) {
    //! A = [[−1, 1], [0, −40]], G = Qc = I, Ts = 1. One exponential of the Van Loan block holds
    //! e^{40} and its product cancels: Qd₀₁ came out −0.0057 in double and −1048576 in float.
    //! Reference: SciPy quad_vec of ∫₀^Ts e^{At} e^{Aᵀt} dt.
    const Matrix<2, 2, T> A{{T{-1}, T{1}}, {T{0}, T{-40}}};
    const auto            Qd = design::discretize_process_noise(A, Matrix<2, 2, T>::identity(), Matrix<2, 2, T>::identity(), T{1});
    REQUIRE(Qd.has_value());
    const double tol = std::is_same_v<T, float> ? 1e-4 : 1e-10;
    CHECK(static_cast<double>((*Qd)(0, 0)) == doctest::Approx(0.432592747514223).epsilon(tol));
    CHECK(static_cast<double>((*Qd)(0, 1)) == doctest::Approx(0.000304878048780494).epsilon(tol));
    CHECK(static_cast<double>((*Qd)(1, 1)) == doctest::Approx(0.0125).epsilon(tol));
}

TEST_CASE("kalman(sys, Ts, noise): a continuous plant and a noise density, against SciPy") {
    constexpr auto         sys = make_continuous_oscillator();
    constexpr Matrix<1, 1> Qc{{0.3}}; // process-noise density
    constexpr Matrix<1, 1> Rd{{0.2}}; // variance of each reading
    constexpr double       Ts = 0.1;

    constexpr auto res = design::kalman(sys, Ts, {.process = Qc, .measurement = Rd});
    static_assert(res.success && res.status == design::DesignStatus::ok);
    static_assert(res.sys.Ts == Ts);

    //! The model it deploys: ZOH plant, G = I with Q = Qd, and the measurement variance as given.
    CHECK(res.sys.A(0, 0) == doctest::Approx(0.990180930582829).epsilon(1e-12));
    CHECK(res.sys.A(1, 0) == doctest::Approx(-0.19443270467639365).epsilon(1e-12));
    CHECK(res.sys.B(1, 0) == doctest::Approx(0.09721635233819682).epsilon(1e-12));
    CHECK(res.sys.G == Matrix<2, 2>::identity());
    CHECK(res.Q(0, 1) == doctest::Approx(1.4176528742916635e-03).epsilon(1e-12));
    CHECK(res.R(0, 0) == doctest::Approx(0.2).epsilon(1e-14));

    //! scipy.linalg.solve_discrete_are on the discretized problem
    CHECK(res.M(0, 0) == doctest::Approx(0.15876599438771183).epsilon(1e-9));
    CHECK(res.M(1, 0) == doctest::Approx(0.13651425878842288).epsilon(1e-9));
    CHECK(res.P(0, 0) == doctest::Approx(0.0377459763462973).epsilon(1e-9));
    CHECK(res.P(0, 1) == doctest::Approx(0.03245571574084469).epsilon(1e-9));
    CHECK(res.P(1, 1) == doctest::Approx(0.13784428598747336).epsilon(1e-9));

    //! It deploys like any Kalman design, in float: the fixed gain on FixedGainObserver.
    constexpr FixedGainObserver<2, 1, 1, float> filt{res.as<float>()};
    CHECK(filt.gain()(0, 0) == doctest::Approx(0.15876599438771183).epsilon(1e-5));

    //! A float plant designs in float.
    constexpr auto res_fd = design::kalman(sys.as<float>(), 0.1f, {.process = Qc.as<float>(), .measurement = Rd.as<float>()});
    static_assert(res_fd.success);
    CHECK(res_fd.M(0, 0) == doctest::Approx(0.15876599438771183).epsilon(1e-4));
}

TEST_CASE("kalman(sys, Ts, noise) is the discrete design on the ZOH plant, the Van Loan Qd and the unscaled variance") {
    constexpr auto         sys = make_continuous_oscillator();
    constexpr Matrix<1, 1> Qc{{0.3}};
    constexpr Matrix<1, 1> sigma2{{0.05}};

    for (const double Ts : {1e-3, 1e-2, 0.1}) {
        const auto res = design::kalman(sys, Ts, {.process = Qc, .measurement = sigma2});
        REQUIRE(res.success);
        CHECK(res.R == sigma2); // no R/Ts: the per-reading variance goes in as given

        const auto                              zoh = discretize(StateSpace<2, 1, 1>{.A = sys.A, .B = sys.B, .C = sys.C}, Ts, DiscretizationMethod::ZOH);
        const StateSpace<2, 1, 1, double, 2, 1> sys_d{.A = zoh->A, .B = zoh->B, .C = sys.C, .D = sys.D, .G = Matrix<2, 2>::identity(), .H = sys.H, .Ts = Ts};
        const auto                              by_hand = design::kalman(sys_d, {.process = design::discretize_process_noise(sys.A, sys.G, Qc, Ts).value(), .measurement = sigma2});
        REQUIRE(by_hand.success);
        for (size_t i = 0; i < 2; ++i) {
            CHECK(res.M(i, 0) == doctest::Approx(by_hand.M(i, 0)).epsilon(1e-12));
        }
    }

    //! Constant-evaluated, in double and float.
    constexpr auto res_c = design::kalman(sys, 0.01, {.process = Qc, .measurement = sigma2});
    static_assert(res_c.success && res_c.R(0, 0) == 0.05);
    constexpr auto res_f = design::kalman(sys.as<float>(), 0.01f, {.process = Qc.as<float>(), .measurement = sigma2.as<float>()});
    static_assert(res_f.success);
    CHECK(res_f.M(0, 0) == doctest::Approx(res_c.M(0, 0)).epsilon(1e-4));
}

TEST_CASE("kalman(sys, Ts, noise) dispatch: one discretization, at a valid Ts") {
    constexpr auto         sys = make_continuous_oscillator();
    constexpr Matrix<1, 1> Qc{{0.3}};
    constexpr Matrix<1, 1> Rc{{0.02}};

    auto discrete = sys;
    discrete.Ts = 0.1;
    const auto twice = design::kalman(discrete, 0.1, {.process = Qc, .measurement = Rc}); // even at the same Ts
    CHECK_FALSE(twice.success);
    CHECK(twice.status == design::DesignStatus::invalid_sample_time);

    for (const double bad : {0.0, -0.1, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
        const auto refused = design::kalman(sys, bad, {.process = Qc, .measurement = Rc});
        CHECK_FALSE(refused.success);
        CHECK(refused.status == design::DesignStatus::invalid_sample_time);
    }

    //! Continuous plant without Ts: still refused (there is no Kalman–Bucy design).
    CHECK(design::kalman(sys, {.process = Qc, .measurement = Rc}).status == design::DesignStatus::invalid_sample_time);

    //! The plant takes discretize's status: a non-finite entry, or a model past the range of T.
    auto bad = sys;
    bad.B(1, 0) = std::numeric_limits<double>::quiet_NaN();
    CHECK(design::kalman(bad, 0.1, {.process = Qc, .measurement = Rc}).status == design::DesignStatus::bad_weights);
    bad = sys;
    bad.G(1, 0) = std::numeric_limits<double>::quiet_NaN();
    CHECK(design::kalman(bad, 0.1, {.process = Qc, .measurement = Rc}).status == design::DesignStatus::bad_weights);
    bad = sys;
    bad.C(0, 0) = std::numeric_limits<double>::quiet_NaN(); // not sampled: refused before the DARE
    CHECK(design::kalman(bad, 0.1, {.process = Qc, .measurement = Rc}).status == design::DesignStatus::bad_weights);
    bad = sys;
    bad.A = Matrix<2, 2>{{1000.0, 1.0}, {0.0, 1000.0}};
    CHECK(design::kalman(bad, 10.0, {.process = Qc, .measurement = Rc}).status == design::DesignStatus::not_converged);
}
