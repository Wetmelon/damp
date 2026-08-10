// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <type_traits>
#include <vector>

#include "damp/analysis/analysis.hpp"
#include "damp/design/stability.hpp"
#include "damp/math/complex.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

TEST_CASE("poles() handles systems larger than 4 states") {
    // Diagonal A with known eigenvalues -1..-6 (6 states > old NX<=4 limit).
    Matrix<6, 6, double> A{};
    const double         expected[6] = {-1.0, -2.0, -3.0, -4.0, -5.0, -6.0};
    for (size_t i = 0; i < 6; ++i) {
        A(i, i) = expected[i];
    }
    auto p = analysis::poles(A);
    for (size_t i = 0; i < 6; ++i) {
        CHECK(p[i].real() == doctest::Approx(expected[i]).epsilon(1e-9));
        CHECK(p[i].imag() == doctest::Approx(0.0).epsilon(1e-9));
    }
    CHECK(analysis::is_stable_continuous(A));
}

TEST_CASE("step/impulse/initial response of 1/(s+1)") {
    // 1/(s+1): A=-1, B=1, C=1, D=0
    TransferFunction<1, 2, double> tf{{1.0}, {1.0, 1.0}};
    auto                           t = analysis::linspace(0.0, 5.0, 501);

    auto s = analysis::step(tf, t);
    auto i = analysis::impulse(tf, t);
    for (size_t k = 0; k < t.size(); ++k) {
        CHECK(s.y[k](0, 0) == doctest::Approx(1.0 - std::exp(-t[k])).epsilon(1e-3));
        CHECK(i.y[k](0, 0) == doctest::Approx(std::exp(-t[k])).epsilon(1e-3));
    }

    // initial response from x0=2: y(t) = 2 e^{-t}
    auto      ss = tf.to_state_space().value();
    ColVec<1> x0{2.0};
    auto      in = analysis::initial(ss, x0, t);
    for (size_t k = 0; k < t.size(); ++k) {
        CHECK(in.y[k][0] == doctest::Approx(2.0 * std::exp(-t[k])).epsilon(1e-3));
    }
}

TEST_CASE("step response of a 2x2 MIMO system separates input channels") {
    // Two decoupled first-order plants: A=diag(-1,-2), B=C=I, D=0.
    // Step on input j must drive only output j: y[k](i,j)=0 for i!=j.
    Matrix<2, 2, double> A{{-1.0, 0.0}, {0.0, -2.0}};
    Matrix<2, 2, double> B{{1.0, 0.0}, {0.0, 1.0}};
    Matrix<2, 2, double> C{{1.0, 0.0}, {0.0, 1.0}};
    const auto           sys = StateSpace{A, B, C};

    const auto s = analysis::step(sys, analysis::linspace(0.0, 5.0, 501));
    for (size_t k = 0; k < s.t.size(); ++k) {
        CHECK(s.y[k](0, 0) == doctest::Approx(1.0 - std::exp(-s.t[k])).epsilon(1e-3));
        CHECK(s.y[k](1, 1) == doctest::Approx(0.5 * (1.0 - std::exp(-2.0 * s.t[k]))).epsilon(1e-3));
        CHECK(s.y[k](0, 1) == doctest::Approx(0.0).epsilon(1e-9));
        CHECK(s.y[k](1, 0) == doctest::Approx(0.0).epsilon(1e-9));
    }
}

TEST_CASE("lsim drives a system with an arbitrary input (MIMO output)") {
    // 1/(s+1): step input via lsim should match the closed-form 1 - e^{-t}.
    // Two-output C = [1; 2] checks the MIMO result shape (y[k] is a ColVec).
    Matrix<1, 1, double> A{{-1.0}};
    Matrix<1, 1, double> B{{1.0}};
    Matrix<2, 1, double> C{{1.0}, {2.0}};
    const auto           sys = StateSpace{A, B, C};

    const auto                t = analysis::linspace(0.0, 5.0, 501);
    const std::vector<double> u(t.size(), 1.0); // unit step
    const auto                res = analysis::lsim(sys, u, t);

    REQUIRE(res.y.size() == t.size());
    for (size_t k = 0; k < t.size(); ++k) {
        const double expected = 1.0 - std::exp(-t[k]);
        CHECK(res.y[k][0] == doctest::Approx(expected).epsilon(1e-3));
        CHECK(res.y[k][1] == doctest::Approx(2.0 * expected).epsilon(1e-3));
    }
}

TEST_CASE("lsim returns empty result when u and time lengths disagree") {
    Matrix<1, 1, double> A{{-1.0}};
    Matrix<1, 1, double> B{{1.0}};
    Matrix<1, 1, double> C{{1.0}};
    const auto           sys = StateSpace{A, B, C};

    const auto                t = analysis::linspace(0.0, 1.0, 11);
    const std::vector<double> u_short(5, 1.0);
    const auto                bad = analysis::lsim(sys, u_short, t);
    CHECK(bad.t.empty());
    CHECK(bad.y.empty());
    CHECK(bad.x.empty());
}

TEST_CASE("discrete step/lsim empty when host grid Ts disagrees with sys.Ts") {
    // Discrete plant at 100 Hz; grid at 50 Hz must not silently re-label steps.
    StateSpace<1, 1, 1> sys{
        .A = Matrix<1, 1>{{0.9}},
        .B = Matrix<1, 1>{{0.1}},
        .C = Matrix<1, 1>{{1.0}},
        .D = Matrix<1, 1>{{0.0}},
        .Ts = 0.01,
    };
    const auto t_bad = analysis::linspace(0.0, 0.1, 6); // dt = 0.02 ≠ 0.01
    CHECK(analysis::step(sys, t_bad).t.empty());
    CHECK(analysis::impulse(sys, t_bad).t.empty());

    const std::vector<double> u(t_bad.size(), 1.0);
    CHECK(analysis::lsim(sys, u, t_bad).t.empty());

    // Matching grid works
    const auto t_ok = analysis::linspace(0.0, 0.1, 11); // dt = 0.01
    const auto s = analysis::step(sys, t_ok);
    REQUIRE(s.t.size() == t_ok.size());
    REQUIRE_FALSE(s.y.empty());
}

TEST_CASE("TF step empty when companion realization fails") {
    TransferFunction<1, 2, double> bad{{1.0}, {1.0, 0.0}}; // leading den = 0
    const auto                     t = analysis::linspace(0.0, 1.0, 11);
    CHECK(analysis::step(bad, t).t.empty());
    CHECK(analysis::impulse(bad, t).t.empty());
}

TEST_CASE("stepinfo on a first-order system: no overshoot, known rise/settling") {
    // 1/(s+1): y=1-e^{-t}. Rise(10-90%)=ln(9)≈2.197, settle(2%)=ln(50)≈3.912.
    TransferFunction<1, 2, double> tf{{1.0}, {1.0, 1.0}};
    const auto                     info = analysis::stepinfo(*tf.to_state_space(), analysis::linspace(0.0, 10.0, 2001));

    CHECK(info.overshoot == doctest::Approx(0.0).epsilon(1e-6));
    CHECK(info.rise_time == doctest::Approx(2.197).epsilon(2e-2));
    CHECK(info.settling_time == doctest::Approx(3.912).epsilon(2e-2));
}

TEST_CASE("stepinfo on an underdamped 2nd-order system matches the overshoot formula") {
    // 1/(s^2+s+1): wn=1, zeta=0.5 → overshoot = exp(-pi*zeta/sqrt(1-zeta^2))*100 ≈ 16.3%.
    TransferFunction<1, 3, double> tf{{1.0}, {1.0, 1.0, 1.0}};
    const auto                     y = analysis::step(*tf.to_state_space(), analysis::linspace(0.0, 20.0, 4001));

    std::vector<double> sig, t = y.t;
    for (const auto& yk : y.y) {
        sig.push_back(yk(0, 0));
    }
    const auto info = analysis::stepinfo(sig, t, sig.back());

    const double expected_os = std::exp(-std::numbers::pi_v<double> * 0.5 / std::sqrt(1.0 - 0.25)) * 100.0;
    CHECK(info.overshoot == doctest::Approx(expected_os).epsilon(2e-2));
    CHECK(info.peak > 1.0); // overshoots past the steady-state value of 1
}

TEST_CASE("lsiminfo reports extremes and settling of a signal") {
    // 1/(s+1) step: min at t=0 (=0), max →1, settles within 2% at ln(50)≈3.912.
    TransferFunction<1, 2, double> tf{{1.0}, {1.0, 1.0}};
    const auto                     resp = analysis::step(*tf.to_state_space(), analysis::linspace(0.0, 10.0, 2001));

    std::vector<double> sig, t = resp.t;
    for (const auto& yk : resp.y) {
        sig.push_back(yk(0, 0));
    }
    REQUIRE_FALSE(sig.empty());
    const auto info = analysis::lsiminfo(sig, t, sig.back());

    CHECK(info.min == doctest::Approx(0.0).epsilon(1e-9));
    CHECK(info.min_time == doctest::Approx(0.0).epsilon(1e-9));
    CHECK(info.max == doctest::Approx(1.0).epsilon(1e-2));
    CHECK(info.settling_time == doctest::Approx(3.912).epsilon(2e-2));
}

TEST_CASE("pzmap of a transfer function finds poles and zeros") {
    // (s+2) / (s^2 + 3s + 2) = (s+2)/((s+1)(s+2)): zero at -2, poles at -1,-2.
    TransferFunction<2, 3, double> tf{{2.0, 1.0}, {2.0, 3.0, 1.0}};
    const auto                     pz = analysis::pzmap(tf);

    REQUIRE(pz.zeros.size() == 1);
    REQUIRE(pz.poles.size() == 2);
    CHECK(pz.zeros[0].real() == doctest::Approx(-2.0).epsilon(1e-9));

    // Poles are {-1, -2} in some order.
    double pmin = std::min(pz.poles[0].real(), pz.poles[1].real());
    double pmax = std::max(pz.poles[0].real(), pz.poles[1].real());
    CHECK(pmin == doctest::Approx(-2.0).epsilon(1e-9));
    CHECK(pmax == doctest::Approx(-1.0).epsilon(1e-9));
}

TEST_CASE("pzmap of a state-space system returns eigenvalue poles") {
    Matrix<2, 2, double> A{{0.0, 1.0}, {-2.0, -3.0}}; // char poly s^2+3s+2 → -1,-2
    const auto           pz = analysis::pzmap(A);

    REQUIRE(pz.poles.size() == 2);
    CHECK(pz.zeros.empty()); // transmission zeros not computed
    double pmin = std::min(pz.poles[0].real(), pz.poles[1].real());
    double pmax = std::max(pz.poles[0].real(), pz.poles[1].real());
    CHECK(pmin == doctest::Approx(-2.0).epsilon(1e-9));
    CHECK(pmax == doctest::Approx(-1.0).epsilon(1e-9));
}

TEST_CASE("linspace and logspace support float") {
    const auto lin = analysis::linspace(0.0f, 1.0f, 5);
    const auto log = analysis::logspace(1.0f, 100.0f, 3);

    CHECK(std::is_same_v<typename decltype(lin)::value_type, float>);
    CHECK(std::is_same_v<typename decltype(log)::value_type, float>);

    CHECK(lin.front() == doctest::Approx(0.0f));
    CHECK(lin.back() == doctest::Approx(1.0f));
    CHECK(log.front() == doctest::Approx(1.0f));
    CHECK(log.back() == doctest::Approx(100.0f));
}

TEST_CASE("geomspace matches geometric progression and endpoints") {
    const auto g = analysis::geomspace(1.0, 1000.0, 4);
    REQUIRE(g.size() == 4);
    CHECK(g[0] == doctest::Approx(1.0));
    CHECK(g[1] == doctest::Approx(10.0).epsilon(1e-12));
    CHECK(g[2] == doctest::Approx(100.0).epsilon(1e-12));
    CHECK(g[3] == doctest::Approx(1000.0));

    // Same span as damp logspace(value endpoints) for positive bases.
    const auto lg = analysis::logspace(1.0, 1000.0, 4);
    REQUIRE(lg.size() == g.size());
    for (size_t i = 0; i < g.size(); ++i) {
        CHECK(g[i] == doctest::Approx(lg[i]).epsilon(1e-12));
    }

    CHECK(analysis::geomspace(-1.0, 1.0, 5).empty()); // mixed signs
    CHECK(analysis::geomspace(0.0, 10.0, 5).empty()); // zero endpoint
}

TEST_CASE("arange is half-open like NumPy") {
    const auto a = analysis::arange(0.0, 5.0, 1.0);
    REQUIRE(a.size() == 5);
    CHECK(a.front() == doctest::Approx(0.0));
    CHECK(a.back() == doctest::Approx(4.0));

    const auto b = analysis::arange(1.0, 4.0); // step 1
    REQUIRE(b.size() == 3);
    CHECK(b[0] == doctest::Approx(1.0));
    CHECK(b[2] == doctest::Approx(3.0));

    const auto c = analysis::arange(5.0); // 0..5
    REQUIRE(c.size() == 5);
    CHECK(c.back() == doctest::Approx(4.0));

    const auto d = analysis::arange(5.0, 0.0, -2.0);
    REQUIRE(d.size() == 3);
    CHECK(d[0] == doctest::Approx(5.0));
    CHECK(d[2] == doctest::Approx(1.0));

    CHECK(analysis::arange(0.0, 5.0, 0.0).empty());
    CHECK(analysis::arange(5.0, 0.0, 1.0).empty());
}

TEST_SUITE("Stability Analysis") {
    TEST_CASE("Discrete stability check - stable system") {
        // Stable discrete system: eigenvalues inside unit circle
        Matrix<2, 2> A{{0.5, 0.0}, {0.0, 0.8}};
        CHECK(stability::is_stable_discrete(A) == true);
        CHECK(stability::stability_margin_discrete(A) > 0.0);
    }

    TEST_CASE("Discrete stability check - unstable system") {
        // Unstable: eigenvalue outside unit circle
        Matrix<2, 2> A{{1.5, 0.0}, {0.0, 0.5}};
        CHECK(stability::is_stable_discrete(A) == false);
        CHECK(stability::stability_margin_discrete(A) < 0.0);
    }

    TEST_CASE("closed_loop_poles returns finite poles when eigen converges") {
        // Diagonal A, zero K ⇒ poles are the open-loop eigenvalues (well-conditioned).
        // On non-convergence the API zeros the vector — here we assert the happy path.
        Matrix<2, 2> A{{0.5, 0.0}, {0.0, 0.8}};
        Matrix<2, 1> B{{1.0}, {0.0}};
        Matrix<1, 2> K{{0.0, 0.0}};

        const auto poles = stability::closed_loop_poles(A, B, K);
        CHECK(std::isfinite(poles[0].real()));
        CHECK(std::isfinite(poles[0].imag()));
        CHECK(std::isfinite(poles[1].real()));
        CHECK(std::isfinite(poles[1].imag()));
        // Order is not guaranteed; both diagonal entries must appear.
        const double r0 = poles[0].real();
        const double r1 = poles[1].real();
        const double lo = (r0 < r1) ? r0 : r1;
        const double hi = (r0 < r1) ? r1 : r0;
        CHECK(lo == doctest::Approx(0.5));
        CHECK(hi == doctest::Approx(0.8));
    }
}
