// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Every design family on a delay-augmented model: the design on discretize(..., InputDelay /
// OutputDelay) holds the loop, the same design on the undelayed model does not.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <type_traits>

#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/controllers/offset_free_lqg.hpp"
#include "damp/controllers/offset_free_mpc.hpp"
#include "damp/design/synthesis.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/estimation/mhe.hpp"
#include "damp/matlab.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "doctest.h"

using namespace damp;

namespace {

constexpr double kTs = 0.01;

/// Lightly damped mass-spring (ω = 2 rad/s, ζ = 0.1): position, speed; force in; position out.
/// Process noise on both states, measurement noise on the reading.
constexpr StateSpace<2, 1, 1, double, 2, 1> mass_spring() {
    StateSpace<2, 1, 1, double, 2, 1> s{};
    s.A(0, 1) = 1.0;
    s.A(1, 0) = -4.0;
    s.A(1, 1) = -0.4;
    s.B(1, 0) = 1.0;
    s.C(0, 0) = 1.0;
    return s; // G = I, H = I by default
}

/// The same mass-spring with a speed sensor too: position and speed out.
constexpr StateSpace<2, 1, 2, double, 2, 2> mass_spring_two_sensors() {
    StateSpace<2, 1, 2, double, 2, 2> s{};
    s.A = mass_spring().A;
    s.B = mass_spring().B;
    s.C = Matrix<2, 2>::identity();
    return s;
}

/// A structured input sequence.
double u_at(int k) {
    return (k < 0) ? 0.0 : std::sin(0.037 * k) + ((k / 150) % 2 == 0 ? 0.5 : -0.5);
}

const auto kUndelayed = discretize(mass_spring(), kTs).value();
const auto kDelay2 = discretize(mass_spring(), kTs, InputDelay<2>{}).value();           // two whole samples
const auto kDelay14 = discretize(mass_spring(), kTs, InputDelay<2>{1.4 * kTs}).value(); // 1.4 samples
const auto kBoth = discretize(mass_spring(), kTs, InputDelay<2>{1.4 * kTs}, OutputDelay<1>{}).value();
const auto kBothTwoSensors = discretize(mass_spring_two_sensors(), kTs, InputDelay<2>{1.4 * kTs}, OutputDelay<1>{}).value();

/// The true plant: a delay-augmented model (exact; see test_input_delay.cpp).
template<size_t N>
struct Truth {
    StateSpace<N, 1, 1, double, 2, 1> m{};
    ColVec<N>                         x{};
    [[nodiscard]] double              y() const { return (m.C * x)(0, 0); }
    [[nodiscard]] double              position() const { return x[0]; }
    void                              step(double u, double load = 0.0) {
        x = m.A * x + m.B * ColVec<1>{{u + load}};
    }
};

/// Largest |position| over the last quarter of a run: ≈ 0 when the loop holds.
struct Tail {
    double worst = 0.0;
    void   add(int k, int n, double v) {
        if (4 * k >= 3 * n) {
            worst = std::max(worst, std::abs(v));
        }
    }
};

} // namespace

TEST_SUITE("Delay-aware design") {
    TEST_CASE("an output delay reads the plant output D samples late, feedthrough included") {
        auto plant = mass_spring();
        plant.D(0, 0) = 0.5;
        const auto base = discretize(plant, kTs).value();
        const auto late = discretize(plant, kTs, OutputDelay<3>{}).value();
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(late)>, StateSpace<5, 1, 1, double, 2, 1>>);
        ColVec<2> x{};
        ColVec<5> xa{};
        double    y_hist[300]{};
        double    worst = 0.0;
        for (int k = 0; k < 300; ++k) {
            const ColVec<1> u{{u_at(k)}};
            y_hist[k] = (base.C * x + base.D * u)(0, 0);
            const double ya = (late.C * xa + late.D * u)(0, 0);
            worst = std::max(worst, std::abs(ya - ((k >= 3) ? y_hist[k - 3] : 0.0)));
            x = base.A * x + base.B * u;
            xa = late.A * xa + late.B * u;
        }
        CHECK(worst < 1e-12);
        for (size_t r = 2; r < 5; ++r) { // the register is deterministic
            CHECK(late.G(r, 0) == 0.0);
            CHECK(late.G(r, 1) == 0.0);
        }
        CHECK(late.H(0, 0) == 1.0); // measurement noise stays on the reading
        CHECK(discretize(plant, kTs, OutputDelay<0>{}).value().A == base.A);
        CHECK_FALSE(discretize(base, 2.0 * kTs, OutputDelay<1>{}).has_value()); // a discrete model at another Ts
        CHECK(discretize(base, kTs, OutputDelay<1>{}).has_value());
    }

    TEST_CASE("input and output delay together: the input-delayed output, one sample late") {
        const auto in_only = kDelay14;
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(kBoth)>, StateSpace<5, 1, 1, double, 2, 1>>);
        ColVec<4> x{};
        ColVec<5> xb{};
        double    prev = 0.0;
        double    worst = 0.0;
        for (int k = 0; k < 300; ++k) {
            worst = std::max(worst, std::abs((kBoth.C * xb)(0, 0) - prev));
            prev = (in_only.C * x)(0, 0);
            x = in_only.A * x + in_only.B * ColVec<1>{{u_at(k)}};
            xb = kBoth.A * xb + kBoth.B * ColVec<1>{{u_at(k)}};
        }
        CHECK(worst < 1e-12);
    }

    TEST_CASE("LQR on the delayed model holds a loop that the undelayed design loses") {
        const Matrix<2, 2> Q{{1e4, 0.0}, {0.0, 0.0}};
        const Matrix<1, 1> R{{1e-4}};
        const auto         naive = design::discrete_lqr(kUndelayed.A, kUndelayed.B, Q, R);
        const auto         aware = design::discrete_lqr(kDelay14.A, kDelay14.B, Q, R); // Q sized to the plant
        REQUIRE(naive.success);
        REQUIRE(aware.success);
        CHECK(aware.is_stable());

        constexpr int               n = 600;
        Truth<4>                    a{kDelay14, ColVec<4>{1.0, 0.0, 0.0, 0.0}};
        Truth<4>                    b = a;
        StateFeedback<4, 1, double> ka{aware};
        StateFeedback<2, 1, double> kb{naive};
        Tail                        ta;
        Tail                        tb;
        for (int k = 0; k < n; ++k) {
            // State feedback on the delayed model: the plant state plus the commands already sent.
            a.step(ka.control(a.x)[0]);
            b.step(std::clamp(kb.control(ColVec<2>{b.x[0], b.x[1]})[0], -1e6, 1e6));
            ta.add(k, n, a.position());
            tb.add(k, n, b.position());
        }
        CHECK(ta.worst < 1e-6);
        CHECK(tb.worst > 1.0); // unstable through 1.4 samples of delay
    }

    TEST_CASE("Kalman on the delayed model tracks the true state; the undelayed one lags") {
        const Matrix<2, 2> Qn{{1e-6, 0.0}, {0.0, 1e-4}};
        const Matrix<1, 1> Rn{{1e-6}};
        SUBCASE("input delay") {
            SteadyStateKalmanFilter<4, 1, 1, double, 2, 1> aware{design::kalman(kDelay2, Qn, Rn)};
            SteadyStateKalmanFilter<2, 1, 1, double, 2, 1> naive{design::kalman(kUndelayed, Qn, Rn)};
            Truth<4>                                       t{kDelay2};
            Tail                                           ea;
            Tail                                           eb;
            constexpr int                                  n = 1200;
            for (int k = 0; k < n; ++k) {
                const ColVec<1> y{{t.y()}};
                const ColVec<1> u{{u_at(k)}};
                aware.update(y);
                naive.update(y);
                ea.add(k, n, aware.state()[1] - t.x[1]);
                eb.add(k, n, naive.state()[1] - t.x[1]);
                aware.predict(u);
                naive.predict(u);
                t.step(u_at(k));
            }
            CHECK(ea.worst < 1e-9);
            CHECK(eb.worst > 1e-3);
        }
        SUBCASE("output delay") {
            const auto                                     late = discretize(mass_spring(), kTs, OutputDelay<2>{}).value();
            SteadyStateKalmanFilter<4, 1, 1, double, 2, 1> aware{design::kalman(late, Qn, Rn)};
            SteadyStateKalmanFilter<2, 1, 1, double, 2, 1> naive{design::kalman(kUndelayed, Qn, Rn)};
            Truth<4>                                       t{late};
            Tail                                           ea;
            Tail                                           eb;
            constexpr int                                  n = 1200;
            for (int k = 0; k < n; ++k) {
                const ColVec<1> y{{t.y()}};
                const ColVec<1> u{{u_at(k)}};
                aware.update(y);
                naive.update(y);
                ea.add(k, n, aware.state()[0] - t.x[0]); // the current position, from a stale reading
                eb.add(k, n, naive.state()[0] - t.x[0]);
                aware.predict(u);
                naive.predict(u);
                t.step(u_at(k));
            }
            CHECK(ea.worst < 1e-9);
            CHECK(eb.worst > 1e-3);
        }
    }

    TEST_CASE("LQG regulates through input and output delay; the undelayed LQG goes unstable") {
        const Matrix<2, 2> Q{{1e4, 0.0}, {0.0, 0.0}};
        const Matrix<1, 1> R{{1e-4}};
        const Matrix<2, 2> Qn = Matrix<2, 2>::identity();
        const Matrix<1, 1> Rn{{1e-6}};
        const auto         aware = design::discrete_lqg(kBoth, Q, R, {.process = Qn, .measurement = Rn});
        const auto         naive = design::discrete_lqg(kUndelayed, Q, R, {.process = Qn, .measurement = Rn});
        REQUIRE(aware.success);
        REQUIRE(naive.success);
        LQG<5, 1, 1, double, 2, 1> ca{aware};
        LQG<2, 1, 1, double, 2, 1> cb{naive};
        Truth<5>                   a{kBoth, ColVec<5>{0.5, 0.0, 0.0, 0.0, 0.0}};
        Truth<5>                   b = a;
        Tail                       ta;
        Tail                       tb;
        constexpr int              n = 2000;
        for (int k = 0; k < n; ++k) {
            a.step(ca.step(ColVec<1>{{a.y()}})[0]);
            b.step(std::clamp(cb.step(ColVec<1>{{b.y()}})[0], -1e6, 1e6));
            ta.add(k, n, a.position());
            tb.add(k, n, b.position());
        }
        CHECK(ta.worst < 1e-6);
        CHECK(tb.worst > 1.0);
    }

    TEST_CASE("LQGI tracks a step through the delay, with a plant-sized state reference") {
        Matrix<3, 3> Qp{}; // [x; ∫e], sized to the plant: the delay states get no cost
        Qp(0, 0) = 1e3;
        Qp(2, 2) = 1e2;
        const Matrix<1, 1> R{{1e-4}};
        const Matrix<2, 2> Qn = Matrix<2, 2>::identity();
        const Matrix<1, 1> Rn{{1e-6}};
        const auto         aware = design::discrete_lqgi(kDelay2, Qp, R, {.process = Qn, .measurement = Rn});
        const auto         naive = design::discrete_lqgi(kUndelayed, Qp, R, {.process = Qn, .measurement = Rn});
        REQUIRE(aware.success);
        REQUIRE(naive.success);
        LQGI<4, 1, 1, double, 2, 1> ca{aware};
        LQGI<2, 1, 1, double, 2, 1> cb{naive};
        Truth<4>                    a{kDelay2};
        Truth<4>                    b{kDelay2};
        const ColVec<1>             r{{0.3}};
        const ColVec<2>             x_ref{0.3, 0.0}; // the plant states only
        Tail                        ta;
        Tail                        tb;
        constexpr int               n = 2000;
        for (int k = 0; k < n; ++k) {
            a.step(ca.step(r, ColVec<1>{{a.y()}}, x_ref)[0], 0.2); // with an input load
            b.step(std::clamp(cb.step(r, ColVec<1>{{b.y()}}, x_ref)[0], -1e6, 1e6), 0.2);
            ta.add(k, n, a.position() - 0.3);
            tb.add(k, n, b.position() - 0.3);
        }
        CHECK(ta.worst < 1e-4);
        CHECK(tb.worst > 1e-2);
    }

    TEST_CASE("LQI / LQGI integrate a subset of the outputs on a delayed model, weight sized to the plant") {
        // Position and speed sensors, both delayed; integrate the position only. The weight on
        // [x_plant; ∫e_position] is 3×3; the model has 2 plant + 2 input-delay + 2 output-delay states.
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(kBothTwoSensors)>, StateSpace<6, 1, 2, double, 2, 2>>);
        Matrix<3, 3> Qp{};
        Qp(0, 0) = 1e3;
        Qp(2, 2) = 1e2;
        Matrix<7, 7> Qp_full{}; // the same weight, padded by hand
        Qp_full(0, 0) = 1e3;
        Qp_full(6, 6) = 1e2;
        const Matrix<1, 1> R{{1e-4}};
        const Matrix<1, 2> pick_position{{1.0, 0.0}};
        const Matrix<2, 2> Qn = Matrix<2, 2>::identity();
        const Matrix<2, 2> Rn{{1e-6, 0.0}, {0.0, 1e-4}};

        const auto by_hand = design::discrete_lqi(kBothTwoSensors, Qp_full, R, {0});
        REQUIRE(by_hand.success);
        const auto exact = [&](const auto& K) { return (K - by_hand.K).norm() == 0.0; };
        CHECK(exact(design::discrete_lqi(kBothTwoSensors, Qp, R, {0}).K));
        CHECK(exact(design::discrete_lqi(kBothTwoSensors, Qp, R, pick_position).K));
        CHECK(exact(design::discrete_lqgi(kBothTwoSensors, Qp, R, {.process = Qn, .measurement = Rn}, {0}).lqi.K));
        CHECK(exact(design::discrete_lqgi(kBothTwoSensors, Qp, R, {.process = Qn, .measurement = Rn}, pick_position).lqi.K));
        CHECK(exact(design::lqi_bundle(kBothTwoSensors, Qp, R, {0}).design.K));
        CHECK(exact(design::lqi_bundle(kBothTwoSensors, Qp, R, pick_position).design.K));
        CHECK(exact(design::lqgi_bundle(kBothTwoSensors, Qp, R, {.process = Qn, .measurement = Rn}, {0}).design.lqi.K));
        CHECK(exact(design::lqgi_bundle(kBothTwoSensors, Qp, R, {.process = Qn, .measurement = Rn}, pick_position).design.lqi.K));
        CHECK(exact(matlab::lqi(kBothTwoSensors, Qp, R, {0}).K));
        CHECK(exact(matlab::lqi(kBothTwoSensors, Qp, R, pick_position).K));
        CHECK(exact(matlab::lqgtrack(kBothTwoSensors, Qp, R, Qn, Rn, {0}).lqi.K));
        CHECK(exact(matlab::lqgtrack(kBothTwoSensors, Qp, R, Qn, Rn, pick_position).lqi.K));

        // From the two delayed readings alone, under an input load, the position settles on r.
        const auto aware = design::discrete_lqgi(kBothTwoSensors, Qp, R, {.process = Qn, .measurement = Rn}, {0});
        REQUIRE(aware.success);
        LQGI<6, 1, 2, double, 2, 2, 1> ctrl{aware};
        ColVec<6>                      x{};
        const ColVec<1>                r{{0.3}};
        const ColVec<2>                x_ref{0.3, 0.0}; // the plant states only
        Tail                           tail;
        constexpr int                  n = 2000;
        for (int k = 0; k < n; ++k) {
            const auto u = ctrl.step(r, kBothTwoSensors.C * x, x_ref);
            x = kBothTwoSensors.A * x + kBothTwoSensors.B * ColVec<1>{{u[0] + 0.2}};
            tail.add(k, n, x[0] - 0.3);
        }
        CHECK(tail.worst < 1e-4);
    }

    TEST_CASE("MPC on the delayed model tracks; the undelayed MPC rings") {
        design::MPCWeights<1, 1> w{};
        w.Qy = Matrix<1, 1>{{1e4}};
        w.Rdu = Matrix<1, 1>{{1e-4}};
        const auto aware = design::state_mpc<40, 10>(kDelay14, w);
        const auto naive = design::state_mpc<40, 10>(kUndelayed, w);
        REQUIRE(aware.success);
        REQUIRE(naive.success);
        MPC             ca{aware};
        MPC             cb{naive};
        Truth<4>        a{kDelay14};
        Truth<4>        b{kDelay14};
        const ColVec<1> r{{0.2}};
        Tail            ta;
        Tail            tb;
        constexpr int   n = 800;
        for (int k = 0; k < n; ++k) {
            a.step(ca.control(r, a.x)[0]); // the commands already sent are part of the state
            b.step(std::clamp(cb.control(r, ColVec<2>{b.x[0], b.x[1]})[0], -1e6, 1e6));
            ta.add(k, n, a.position() - 0.2);
            tb.add(k, n, b.position() - 0.2);
        }
        CHECK(ta.worst < 1e-6);
        CHECK(tb.worst > 1e-2);
    }

    TEST_CASE("OffsetFreeMPC on the delayed model rejects a load through the delay") {
        design::MPCWeights<1, 1> w{};
        w.Qy = Matrix<1, 1>{{1e3}};
        w.Rdu = Matrix<1, 1>{{1e-7}};
        const Matrix<2, 2> Qp{{1e-4, 0.0}, {0.0, 1e-2}};
        const Matrix<1, 1> Qd{{1e-2}};
        const Matrix<1, 1> Rn{{1e-6}};
        const auto         aware = design::mpc<40, 10>(kBoth, w, {}, design::OffsetFreeNoise{Qp, Rn, Qd}); // Q_process sized to the plant
        const auto         naive = design::mpc<40, 10>(kUndelayed, w, {}, design::OffsetFreeNoise{Qp, Rn, Qd});
        REQUIRE(aware.success);
        REQUIRE(naive.success);
        OffsetFreeMPC   ca{aware};
        OffsetFreeMPC   cb{naive};
        Truth<5>        a{kBoth};
        Truth<5>        b{kBoth};
        const ColVec<1> r{{0.2}};
        Tail            ta;
        Tail            tb;
        constexpr int   n = 2000;
        for (int k = 0; k < n; ++k) {
            const double load = (k > n / 4) ? 0.5 : 0.0;
            a.step(ca.control(r, ColVec<1>{{a.y()}})[0], load);
            b.step(std::clamp(cb.control(r, ColVec<1>{{b.y()}})[0], -1e6, 1e6), load);
            ta.add(k, n, a.position() - 0.2);
            tb.add(k, n, b.position() - 0.2);
        }
        CHECK(ta.worst < 1e-5);
        CHECK(tb.worst > 1e-2);
    }

    TEST_CASE("OffsetFreeMPC with a plant-sized disturbance model on the delayed model") {
        // The same OffsetFreeModel as the OFLQG case: a load at the plant input, after the
        // delay, the position held. Sized to the plant; the delay states get none of it.
        design::OffsetFreeModel<2, 1, 1, 1, 1> model{};
        model.Bd = kUndelayed.B;
        model.Cz(0, 0) = 1.0;
        design::MPCWeights<1, 1> w{};
        w.Qy = Matrix<1, 1>{{1e3}};
        w.Rdu = Matrix<1, 1>{{1e-7}};
        const Matrix<2, 2> Qp{{1e-4, 0.0}, {0.0, 1e-2}};
        const Matrix<1, 1> Qd{{1e-2}};
        const Matrix<1, 1> Rn{{1e-6}};
        const auto         aware = design::mpc<40, 10>(kBoth, model, w, {}, design::OffsetFreeNoise{Qp, Rn, Qd});
        const auto         naive = design::mpc<40, 10>(kUndelayed, model, w, {}, design::OffsetFreeNoise{Qp, Rn, Qd});
        REQUIRE(aware.success);
        REQUIRE(naive.success);
        OffsetFreeMPC   ca{aware};
        OffsetFreeMPC   cb{naive};
        Truth<5>        a{kBoth};
        Truth<5>        b{kBoth};
        const ColVec<1> r{{0.2}};
        Tail            ta;
        Tail            tb;
        constexpr int   n = 2000;
        for (int k = 0; k < n; ++k) {
            const double load = (k > n / 4) ? 0.5 : 0.0;
            a.step(ca.control(r, ColVec<1>{{a.y()}})[0], load);
            b.step(std::clamp(cb.control(r, ColVec<1>{{b.y()}})[0], -1e6, 1e6), load);
            ta.add(k, n, a.position() - 0.2);
            tb.add(k, n, b.position() - 0.2);
        }
        CHECK(ta.worst < 1e-5);
        CHECK(tb.worst > 1e-2);
        CHECK(ca.disturbance_estimate()(0) == doctest::Approx(0.5).epsilon(1e-3));
    }

    TEST_CASE("OFLQG tracks a moving reference through both delays from the plant-state reference alone") {
        // The model, weights and reference are sized to the plant; only the model is delayed.
        const auto                             plant_d = discretize(mass_spring(), kTs).value();
        design::OffsetFreeModel<2, 1, 1, 1, 1> model{};
        model.Bd = plant_d.B; // a load at the plant input, after the delay
        model.Cz(0, 0) = 1.0;
        const Matrix<2, 2> Q{{1e4, 0.0}, {0.0, 0.0}};
        const Matrix<1, 1> R{{1e-4}};
        const Matrix<2, 2> Qn{{1e-4, 0.0}, {0.0, 1e-2}};
        const Matrix<1, 1> Rn{{1e-6}};
        const Matrix<1, 1> Qd{{1e-2}};
        const auto         aware = design::discrete_oflqg(kBoth, model, Q, R, design::OffsetFreeNoise{Qn, Rn, Qd});
        const auto         naive = design::discrete_oflqg(kUndelayed, model, Q, R, design::OffsetFreeNoise{Qn, Rn, Qd});
        REQUIRE(aware.success);
        REQUIRE(naive.success);

        OFLQG<5, 1, 1, double, 1, 1> ca{aware};
        OFLQG<5, 1, 1, double, 1, 1> by_hand{aware}; // the same law, fed the full delayed reference
        OFLQG<2, 1, 1, double, 1, 1> cb{naive};
        Truth<5>                     a{kBoth};
        Truth<5>                     h{kBoth};
        Truth<5>                     b{kBoth};
        ColVec<5>                    x_full{}; // [x*; delay states], propagated through the model
        constexpr double             w = 1.5;
        Tail                         ta;
        Tail                         tb;
        double                       same = 0.0;
        constexpr int                n = 2000;
        for (int k = 0; k < n; ++k) {
            const double    t = k * kTs;
            const ColVec<2> x_ref{0.2 * std::sin(w * t), 0.2 * w * std::cos(w * t)};
            const ColVec<1> u_ref{{0.2 * (4.0 - w * w) * std::sin(w * t) + 0.4 * 0.2 * w * std::cos(w * t)}};
            x_full[0] = x_ref[0];
            x_full[1] = x_ref[1];
            const double ua = ca.control(ColVec<1>{{a.y()}}, x_ref, u_ref)[0];
            const double uh = by_hand.control(ColVec<1>{{h.y()}}, x_full, u_ref)[0];
            const double ub = std::clamp(cb.control(ColVec<1>{{b.y()}}, x_ref, u_ref)[0], -1e6, 1e6);
            x_full = kBoth.A * x_full + kBoth.B * u_ref;
            same = std::max(same, std::abs(ua - uh));
            a.step(ua, 0.3);
            h.step(uh, 0.3);
            b.step(ub, 0.3);
            ta.add(k, n, a.position() - x_ref[0]);
            tb.add(k, n, b.position() - x_ref[0]);
        }
        CHECK(same < 1e-9);     // the history kept inside is the one a caller would build
        CHECK(ta.worst < 5e-3); // within a few % of the 0.2 amplitude, under the load
        CHECK(tb.worst > 0.1);
    }

    TEST_CASE("MHE on the delayed model, Q sized to the plant, tracks the true state") {
        const Matrix<2, 2> Qn{{1e-6, 0.0}, {0.0, 1e-4}};
        const auto         art = design::mhe<10>(kDelay2, Qn, Matrix<1, 1>{{1e-6}});
        REQUIRE(art.success);
        MHE<4, 1, 1, 10, double> mhe{art};
        Truth<4>                 t{kDelay2};
        Tail                     e;
        constexpr int            n = 1200;
        for (int k = 0; k < n; ++k) {
            const auto x_hat = mhe.update(ColVec<1>{{t.y()}}, ColVec<1>{{u_at(k - 1)}});
            e.add(k, n, x_hat[1] - t.x[1]);
            t.step(u_at(k));
        }
        CHECK(e.worst < 1e-9);
    }

    TEST_CASE("one rule for a plant-sized state reference across LQR, LQI, LQG and LQGI") {
        // Laws without a feedforward input: the past-input states' reference is zero, and the
        // past-output states carry no LQR gain, so LQG's model-filled tail acts like zeros.
        const Matrix<2, 2> Q{{1e4, 0.0}, {0.0, 0.0}};
        const Matrix<1, 1> R{{1e-4}};
        const auto         lqr = design::discrete_lqr(kBoth.A, kBoth.B, Q, R);
        REQUIRE(lqr.success);
        CHECK(std::abs(lqr.K(0, 4)) < 1e-9 * std::abs(lqr.K(0, 0))); // the past output steers nothing
        const StateFeedback<5, 1, double> sf{lqr};
        const ColVec<5>                   x{0.1, -0.2, 0.3, 0.05, 0.02};
        const ColVec<2>                   x_ref{0.2, 0.1};
        const ColVec<5>                   x_ref_full{0.2, 0.1, 0.0, 0.0, 0.0};
        CHECK(sf.control(x_ref, x)[0] == sf.control(x_ref_full, x)[0]);

        Matrix<3, 3> Qi{};
        Qi(0, 0) = 1e3;
        Qi(2, 2) = 1e2;
        Matrix<5, 5> Qi_full{}; // the same weight, padded by hand
        Qi_full(0, 0) = 1e3;
        Qi_full(4, 4) = 1e2;
        const auto lqi_design = design::discrete_lqi(kDelay2, Qi, R);
        CHECK((lqi_design.K - design::discrete_lqi(kDelay2, Qi_full, R).K).norm() == 0.0);
        LQI<4, 1, 1, double> lqi_a{lqi_design};
        LQI<4, 1, 1, double> lqi_b = lqi_a;
        const ColVec<4>      xi{0.1, -0.2, 0.3, 0.05};
        CHECK(lqi_a.control(ColVec<1>{{0.2}}, ColVec<1>{{0.1}}, xi, x_ref)[0] == lqi_b.control(ColVec<1>{{0.2}}, ColVec<1>{{0.1}}, xi, ColVec<4>{0.2, 0.1, 0.0, 0.0})[0]);

        const Matrix<2, 2>         Qn = Matrix<2, 2>::identity();
        const auto                 lqg = design::discrete_lqg(kBoth, Q, R, {.process = Qn, .measurement = Matrix<1, 1>{{1e-6}}});
        LQG<5, 1, 1, double, 2, 1> g_plant{lqg};
        LQG<5, 1, 1, double, 2, 1> g_full{lqg};
        double                     worst = 0.0;
        for (int k = 0; k < 50; ++k) {
            const ColVec<1> y{{0.01 * k}};
            const ColVec<2> r{0.2 * std::sin(0.1 * k), 0.02 * std::cos(0.1 * k)};
            worst = std::max(worst, std::abs(g_plant.step(y, r)[0] - g_full.step(y, ColVec<5>{r[0], r[1], 0.0, 0.0, 0.0})[0]));
        }
        CHECK(worst < 1e-6 * std::abs(lqr.K(0, 0)));
    }

    TEST_CASE("a disturbance in G is discretized with the plant and kept off the delay states") {
        auto plant = mass_spring();
        plant.G = Matrix<2, 2>{{0.0, 0.0}, {0.0, 1.0}}; // a force on the mass, the second column
        const auto plain = discretize(plant, kTs).value();
        const auto delayed = discretize(plant, kTs, InputDelay<2>{1.4 * kTs}, OutputDelay<1>{}).value();
        for (size_t i = 0; i < 5; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                CHECK(delayed.G(i, j) == ((i < 2) ? plain.G(i, j) : 0.0)); // not delayed with the command
            }
        }
        CHECK(std::abs(plain.G(0, 1) - plain.B(0, 0)) < 1e-15); // the same ZOH as the input
        CHECK(std::abs(plain.G(1, 1) - plain.B(1, 0)) < 1e-15);
    }
}
