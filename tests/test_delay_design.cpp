// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Design families on a plant with an input delay. The delay is part of the model,
// discretize(plant, Ts, InputDelay<N>{θ}): NX + N·NU states, the registers v(k−1) … v(k−N) last.
// Every design takes that model like any other plant, with zero weight on the register. A
// measured-state law reads [x; u(k−1)], the last applied input appended by the caller; a Kalman
// filter estimates the register like any other state. The design on the delayed model holds
// the loop; the same design on the delay-blind model does not. The true plants are the
// register form (tests/register_form.hpp), built independently of discretize.

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/controllers/oflqg.hpp"
#include "damp/controllers/ofmpc.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/riccati.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/estimation/luenberger.hpp"
#include "damp/estimation/mhe.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "doctest.h"
#include "register_form.hpp"

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

/// A structured input sequence.
double u_at(int k) {
    return (k < 0) ? 0.0 : std::sin(0.037 * k) + ((k / 150) % 2 == 0 ? 0.5 : -0.5);
}

const auto kUndelayed = discretize(mass_spring(), kTs).value();
const auto kDelayed = discretize(mass_spring(), kTs, InputDelay{kTs}).value(); // the design model: one whole period
const auto kReg1 = test::register_form<1>(mass_spring(), kTs, kTs);            // the true plant: one whole period
const auto kReg2 = test::register_form<2>(mass_spring(), kTs, 2.0 * kTs);      // two periods, a model built by hand

/// The weight on the position alone, on the delayed model: none on the speed or the register.
constexpr Matrix<3, 3> position_weight(double q) {
    return Matrix<3, 3>{{q, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
}

/// The process noise on the plant states as given (G = I), as the offset-free designs model it;
/// the register of the delayed model gets none.
const auto kUndelayedG = [] {
    auto s = kUndelayed;
    s.G = Matrix<2, 2>::identity();
    return s;
}();
const auto kDelayedG = [] {
    auto s = kDelayed;
    s.G = Matrix<3, 2>{{1.0, 0.0}, {0.0, 1.0}, {0.0, 0.0}};
    return s;
}();

/// A load at the plant input, not delayed with the command: B_d = Γ on the plant states, 0 on the register.
const Matrix<3, 1> kLoadOnPlant{{kUndelayed.B(0, 0)}, {kUndelayed.B(1, 0)}, {0.0}};

/// The true plant: a delay-augmented model (exact; see test_input_delay.cpp). A load adds to the
/// command, so it is delayed with it.
template<size_t N>
struct Truth {
    StateSpace<N, 1, 1, double, 2, 1> m{};
    ColVec<N>                         x{};
    [[nodiscard]] ColVec<1>           y() const { return ColVec<1>(m.C * x); }
    [[nodiscard]] double              position() const { return x[0]; }
    [[nodiscard]] ColVec<2>           plant() const { return ColVec<2>{x[0], x[1]}; }
    void                              step(double u, double load = 0.0) { x = m.A * x + m.B * ColVec<1>{{u + load}}; }
};

/// The delayed model's state as a measured-state law reads it: the plant state and the input applied last period.
ColVec<3> with_last_input(const ColVec<2>& x, double u_last) {
    return ColVec<3>{x[0], x[1], u_last};
}

/// Largest |value| over the last quarter of a run: ≈ 0 when the loop holds.
struct Tail {
    double worst = 0.0;
    void   add(int k, int n, double v) {
        if (4 * k >= 3 * n) {
            worst = std::max(worst, std::abs(v));
        }
    }
};

/// The inputs a law applied, newest first: what the caller of a measured-state law keeps for the registers.
template<size_t N>
struct AppliedInputs {
    ColVec<N> u{};

    void push(double v) {
        for (size_t j = N - 1; j > 0; --j) {
            u[j] = u[j - 1];
        }
        u[0] = v;
    }

    [[nodiscard]] ColVec<2 + N> after(const ColVec<2>& x) const {
        ColVec<2 + N> out{};
        out[0] = x[0];
        out[1] = x[1];
        for (size_t j = 0; j < N; ++j) {
            out[2 + j] = u[j];
        }
        return out;
    }
};

/// The position weight on N registers' model: none on the speed or the registers.
template<size_t N>
Matrix<2 + N, 2 + N> position_weight_on(double q) {
    return blkdiag(Matrix<2, 2>{{q, 0.0}, {0.0, 0.0}}, Matrix<N, N>{});
}

/// The measured-state laws, each on the true plant with the true delay.
template<size_t N>
void check_measured_laws(double theta) {
    constexpr size_t NS = 2 + N;
    const auto       model = discretize(mass_spring(), kTs, InputDelay<N>{theta}).value();
    const auto       truth = test::register_form<N>(mass_spring(), kTs, theta);
    ColVec<NS>       x0{};
    x0[0] = 1.0;

    const auto lqr = design::lqr(model, {.Q = position_weight_on<N>(1e4), .R = Matrix<1, 1>{{1e-4}}});
    REQUIRE(lqr.success);
    CHECK(lqr.is_stable());
    LQR              regulator{lqr};
    Truth<NS>        a{truth, x0};
    AppliedInputs<N> applied_a;
    Tail             ta;
    for (int k = 0; k < 1000; ++k) {
        const double u = regulator.step(applied_a.after(a.plant()))[0];
        a.step(u);
        applied_a.push(u);
        ta.add(k, 1000, a.position());
    }
    CHECK(ta.worst < 1e-6);

    const auto lqi = design::lqi(model, {.Q = position_weight_on<N>(1e3), .R = Matrix<1, 1>{{1e-7}}, .Q_integral = Matrix<1, 1>{{1e6}}});
    REQUIRE(lqi.success);
    LQI              integrating{lqi};
    Truth<NS>        b{truth};
    AppliedInputs<N> applied_b;
    Tail             tb;
    for (int k = 0; k < 2000; ++k) {
        const double u = integrating.step(ColVec<1>{{0.3}}, applied_b.after(b.plant()))[0];
        b.step(u, 0.2);
        applied_b.push(u);
        tb.add(k, 2000, b.position() - 0.3);
    }
    CHECK(tb.worst < 1e-4);

    design::MPCWeights<1, 1> w{};
    w.Qy = 1e4;
    w.Rdu = 1e-4;
    const auto mpc = design::mpc<20 + N, 5>(model, w, {});
    REQUIRE(mpc.success);
    MPC              planner{mpc};
    Truth<NS>        c{truth};
    AppliedInputs<N> applied_c;
    Tail             tc;
    for (int k = 0; k < 1500; ++k) {
        const double u = planner.step(ColVec<1>{{0.2}}, applied_c.after(c.plant()))[0];
        c.step(u);
        applied_c.push(u);
        tc.add(k, 1500, c.position() - 0.2);
    }
    CHECK(tc.worst < 1e-4);
}

/// The laws that read the output, offset-free under a load the command carries.
template<size_t N>
void check_output_laws(double theta) {
    constexpr size_t NS = 2 + N;
    const auto       model = discretize(mass_spring(), kTs, InputDelay<N>{theta}).value();
    const auto       truth = test::register_form<N>(mass_spring(), kTs, theta);

    const design::KalmanNoise lqgi_noise{.process = Matrix<2, 2>::identity(), .measurement = Matrix<1, 1>{{1e-6}}};
    const auto                lqgi = design::lqgi(model, {.Q = position_weight_on<N>(1e3), .R = Matrix<1, 1>{{1e-7}}, .Q_integral = Matrix<1, 1>{{1e6}}}, lqgi_noise);
    REQUIRE(lqgi.success);
    LQGI      integrating{lqgi};
    Truth<NS> a{truth};
    Tail      ta;
    for (int k = 0; k < 2000; ++k) {
        a.step(integrating.step(ColVec<1>{{0.3}}, a.y())[0], 0.2);
        ta.add(k, 2000, a.position() - 0.3);
    }
    CHECK(ta.worst < 1e-4);

    const design::KalmanNoise noise{Matrix<2, 2>{{1e-4, 0.0}, {0.0, 1e-2}}, Matrix<1, 1>{{1e-6}}, Matrix<1, 1>{{1e-2}}};
    const auto                oflqg = design::oflqg(model, {.Q = position_weight_on<N>(1e3), .R = Matrix<1, 1>{{1e-4}}}, noise, design::InputDisturbance{});
    REQUIRE(oflqg.success);
    OFLQG     holding{oflqg};
    Truth<NS> b{truth};
    Tail      tb;
    for (int k = 0; k < 2000; ++k) {
        b.step(holding.step(ColVec<1>{{0.3}}, b.y())[0], (k > 500) ? 0.5 : 0.0);
        tb.add(k, 2000, b.position() - 0.3);
    }
    CHECK(tb.worst < 1e-4);

    design::MPCWeights<1, 1> w{};
    w.Qy = 1e3;
    w.Rdu = 1e-4;
    const auto ofmpc = design::ofmpc<20 + N, 5>(model, w, {}, noise, design::InputDisturbance{});
    REQUIRE(ofmpc.success);
    OFMPC     planner{ofmpc};
    Truth<NS> c{truth};
    Tail      tc;
    for (int k = 0; k < 2000; ++k) {
        c.step(planner.step(ColVec<1>{{0.2}}, c.y())[0], (k > 500) ? 0.5 : 0.0);
        tc.add(k, 2000, c.position() - 0.2);
    }
    CHECK(tc.worst < 1e-4);

    // Deployed in float: the same offset-free law, to float's precision.
    OFLQG     single{oflqg.template as<float>()};
    Truth<NS> d{truth};
    Tail      td;
    for (int k = 0; k < 2000; ++k) {
        const ColVec<1, float> y{{static_cast<float>(d.y()[0])}};
        d.step(static_cast<double>(single.step(ColVec<1, float>{{0.3F}}, y)[0]), (k > 500) ? 0.5 : 0.0);
        td.add(k, 2000, d.position() - 0.3);
    }
    CHECK(td.worst < 1e-3);
}

/// design::luenberger on N registers: the plant's poles as asked, the registers' at 0, and the error gone on the true plant.
template<size_t N>
void check_placed_observer(double Ts) {
    constexpr size_t NS = 2 + N;
    CAPTURE(N);
    const auto                       model = discretize(mass_spring(), Ts, InputDelay<N>{static_cast<double>(N) * Ts}).value();
    damp::array<complex<double>, NS> poles{};
    poles[0] = complex<double>{0.9, 0.0};
    poles[1] = complex<double>{0.92, 0.0};
    const auto placed = design::luenberger(model, poles);
    REQUIRE(placed.success);
    CHECK(placed.is_stable());
    Luenberger observer{placed};
    ColVec<NS> x0{};
    x0[0] = 1.0;
    Truth<NS> t{test::register_form<N>(mass_spring(), Ts, static_cast<double>(N) * Ts), x0};
    Tail      e;
    for (int k = 0; k < 600; ++k) {
        observer.update(t.y());
        e.add(k, 600, std::abs(observer.state_estimate()[0] - t.x[0]) + std::abs(observer.state_estimate()[1] - t.x[1]));
        observer.predict(ColVec<1>{{u_at(k)}});
        t.step(u_at(k));
    }
    CHECK(e.worst < 1e-9);
}

/// Every estimator, from a wrong initial state: the registers known, the plant's estimate the undelayed filter's.
template<size_t N>
void check_estimators(double theta) {
    constexpr size_t   NS = 2 + N;
    constexpr int      n = 3000;
    const auto         model = discretize(mass_spring(), kTs, InputDelay<N>{theta}).value();
    const auto         truth = test::register_form<N>(mass_spring(), kTs, theta);
    const Matrix<2, 2> Qn{{1e-6, 0.0}, {0.0, 1e-4}};
    ColVec<NS>         x0{};
    x0[0] = 1.0;

    const design::KalmanNoise<2, 1> noise{.process = Qn, .measurement = 1e-6};
    const auto                      filter = design::kalman(model, noise);
    const auto                      undelayed = design::kalman(kUndelayed, noise);
    REQUIRE(filter.success);
    REQUIRE(undelayed.success);
    CHECK((Matrix<2, 1>(filter.M.template block<2, 1>(0, 0)) - undelayed.M).norm() <= 1e-9 * undelayed.M.norm());
    CHECK(Matrix<N, 1>(filter.M.template block<N, 1>(2, 0)).norm() == 0.0); // a command already sent is known

    const auto run = [&](auto& est, const char* name) {
        CAPTURE(name);
        Truth<NS> t{truth, x0};
        Tail      e;
        for (int k = 0; k < n; ++k) {
            est.update(t.y());
            const auto x_hat = [&] {
                if constexpr (requires { est.state_estimate(); }) {
                    return est.state_estimate();
                } else {
                    return est.state();
                }
            }();
            e.add(k, n, std::abs(x_hat[0] - t.x[0]) + std::abs(x_hat[1] - t.x[1]));
            est.predict(ColVec<1>{{u_at(k)}});
            t.step(u_at(k));
        }
        CHECK(e.worst < 1e-6);
    };
    Luenberger fixed{filter};
    run(fixed, "Kalman gain");
    KalmanFilter recursive{model, {.Q = Qn, .R = 1e-6}};
    run(recursive, "recursive Kalman");

    const auto art = design::mhe<6>(model, {.process = Qn, .measurement = 1e-6});
    REQUIRE(art.success);
    MHE<NS, 1, 1, 6, double> mhe{art};
    Truth<NS>                t{truth, x0};
    Tail                     e;
    for (int k = 0; k < n; ++k) {
        const auto x_hat = mhe.step(t.y(), ColVec<1>{{u_at(k - 1)}});
        e.add(k, n, x_hat[1] - t.x[1]);
        t.step(u_at(k));
    }
    CHECK(e.worst < 1e-6);
}

} // namespace

TEST_SUITE("Design on the delayed model") {
    TEST_CASE("LQR on the delayed model holds a loop that the delay-blind design loses") {
        // The per-step weights Q/R = 1e8 put the delay-blind loop's poles near 100 rad/s
        // (|z| ≈ 0.5, ω·Ts ≈ 1); one whole sample of delay moves two of them outside the unit
        // circle (|z| ≈ 1.22). The design on the delayed model plans with the command landing a
        // period late; it measures the last applied input as its third state.
        const Matrix<1, 1> R{{1e-4}};
        const auto         naive = design::discrete_lqr(kUndelayed.A, kUndelayed.B, Matrix<2, 2>{{1e4, 0.0}, {0.0, 0.0}}, R);
        const auto         aware = design::lqr(kDelayed, {.Q = position_weight(1e4), .R = R});
        REQUIRE(naive.success);
        REQUIRE(aware.success);
        CHECK(aware.is_stable());

        constexpr int     n = 600;
        Truth<3>          a{kReg1, ColVec<3>{1.0, 0.0, 0.0}};
        Truth<3>          b = a;
        LQR               ka{aware};
        LQR<2, 1, double> kb{naive};
        double            u_last = 0.0;
        Tail              ta;
        Tail              tb;
        for (int k = 0; k < n; ++k) {
            const double u = ka.step(with_last_input(a.plant(), u_last))[0];
            a.step(u);
            u_last = u;
            b.step(std::clamp(kb.step(b.plant())[0], -1e6, 1e6));
            ta.add(k, n, a.position());
            tb.add(k, n, b.position());
        }
        CHECK(ta.worst < 1e-6);
        CHECK(tb.worst > 1.0); // unstable through one sample of delay
    }

    TEST_CASE("a Kalman filter on the delayed model, fed the command, converges on the true state; the delay-blind filter lags") {
        // The register is a state of the filter's model: the filter takes the command when it is
        // sent and carries it until the plant receives it. The delay-blind filter takes the
        // command to act at once. A fraction of a period and a whole one; the plant starts away
        // from the estimate.
        const Matrix<2, 2> Qn{{1e-6, 0.0}, {0.0, 1e-4}};
        const Matrix<1, 1> Rn{{1e-6}};
        for (const double theta : {0.4 * kTs, kTs}) {
            const auto                  model = discretize(mass_spring(), kTs, InputDelay{theta}).value();
            Luenberger<3, 1, 1, double> aware{design::kalman(model, {.process = Qn, .measurement = Rn})};
            Luenberger<2, 1, 1, double> naive{design::kalman(kUndelayed, {.process = Qn, .measurement = Rn})};
            Truth<3>                    t{test::register_form<1>(mass_spring(), kTs, theta), ColVec<3>{0.5, -0.2, 0.0}};
            Tail                        ea;
            Tail                        eb;
            constexpr int               n = 2000;
            for (int k = 0; k < n; ++k) {
                const ColVec<1> y = t.y();
                const ColVec<1> u{{u_at(k)}};
                aware.update(y);
                naive.update(y);
                ea.add(k, n, aware.state()[1] - t.x[1]);
                eb.add(k, n, naive.state()[1] - t.x[1]);
                aware.predict(u); // the command as sent
                naive.predict(u);
                t.step(u[0]);
            }
            CHECK(ea.worst < 1e-9);
            CHECK(eb.worst > 1e-3);
        }
    }

    TEST_CASE("LQGI on the delayed model tracks a step under a load about its reference point; the delay-blind LQGI goes unstable") {
        // One sample of input delay. The reference point is on the model's states: the plant at
        // rest at 0.3 and, in the register, the input that holds it there, u* = 4·0.3 (the
        // spring), which is also the input reference. The integrator takes up the load.
        const Matrix<1, 1>        Qi{{1e6}}; // on ∫e
        const Matrix<1, 1>        R{{1e-7}};
        const design::KalmanNoise noise{.process = Matrix<2, 2>::identity(), .measurement = Matrix<1, 1>{{1e-6}}};
        const auto                aware = design::lqgi(kDelayed, {.Q = position_weight(1e3), .R = R, .Q_integral = Qi}, noise);
        const auto                naive = design::lqgi(kUndelayed, {.Q = Matrix<2, 2>{{1e3, 0.0}, {0.0, 0.0}}, .R = R, .Q_integral = Qi}, noise);
        REQUIRE(aware.success);
        REQUIRE(naive.success);
        LQGI            ca{aware};
        LQGI            cb{naive};
        Truth<3>        a{kReg1};
        Truth<3>        b{kReg1};
        const double    u_star = 4.0 * 0.3;
        const ColVec<1> r{{0.3}};
        const ColVec<1> u_ref{{u_star}};
        Tail            ta;
        Tail            tb;
        constexpr int   n = 2000;
        for (int k = 0; k < n; ++k) {
            a.step(ca.step(r, ColVec<3>{0.3, 0.0, u_star}, u_ref, a.y())[0], 0.2);
            b.step(std::clamp(cb.step(r, ColVec<2>{0.3, 0.0}, u_ref, b.y())[0], -1e6, 1e6), 0.2);
            ta.add(k, n, a.position() - 0.3);
            tb.add(k, n, b.position() - 0.3);
        }
        CHECK(ta.worst < 1e-4);
        CHECK(tb.worst > 1e-2);
    }

    TEST_CASE("LQI on [x; u(k−1)] and LQGI on the reading, a fraction of a period and a whole one, offset-free under a load") {
        // The delayed model of each delay; the true plant is the register form of each. The LQI
        // reads the plant state and the input it applied last period; the LQGI the position alone.
        const Matrix<1, 1> Qi{{1e3}}; // on ∫e dt
        const Matrix<1, 1> R{{1e-4}};
        for (const double theta : {0.4 * kTs, kTs}) {
            const auto model = discretize(mass_spring(), kTs, InputDelay{theta}).value();
            const auto lqi = design::lqi(model, {.Q = position_weight(1e3), .R = R, .Q_integral = Qi});
            const auto lqgi = design::lqgi(model, {.Q = position_weight(1e3), .R = R, .Q_integral = Qi}, {.process = Matrix<2, 2>::identity(), .measurement = 1e-6});
            REQUIRE(lqi.success);
            REQUIRE(lqgi.success);
            LQI             li{lqi};
            LQGI            lg{lqgi};
            const auto      truth = test::register_form<1>(mass_spring(), kTs, theta);
            Truth<3>        a{truth};
            Truth<3>        b{truth};
            const ColVec<1> r{{0.3}};
            double          u_last = 0.0;
            Tail            ta;
            Tail            tb;
            constexpr int   n = 3000;
            for (int k = 0; k < n; ++k) {
                const double u = li.step(r, with_last_input(a.plant(), u_last))[0];
                a.step(u, 0.2);
                u_last = u;
                b.step(lg.step(r, b.y())[0], 0.2);
                ta.add(k, n, a.position() - 0.3);
                tb.add(k, n, b.position() - 0.3);
            }
            CHECK(ta.worst < 1e-4);
            CHECK(tb.worst < 1e-4);
        }
    }

    TEST_CASE("MPC on the delayed model tracks; the delay-blind MPC rings") {
        // The state measured, the last applied input appended; the position the controlled
        // output. The same per-step weights on both models. A switch inside the period (0.4·Ts)
        // and a whole period; the true plant is the register form of each.
        design::MPCWeights<1, 1> w{};
        w.Qy = 1e4;
        w.Rdu = 1e-4;
        const auto naive = design::mpc<40, 10>(kUndelayed, w, {});
        REQUIRE(naive.success);
        for (const double theta : {0.4 * kTs, kTs}) {
            const auto aware = design::mpc<40, 10>(discretize(mass_spring(), kTs, InputDelay{theta}).value(), w, {});
            REQUIRE(aware.success);
            MPC             ca{aware};
            MPC             cb{naive};
            const auto      truth = test::register_form<1>(mass_spring(), kTs, theta);
            Truth<3>        a{truth};
            Truth<3>        b{truth};
            const ColVec<1> r{0.2};
            double          u_last = 0.0;
            Tail            ta;
            Tail            tb;
            constexpr int   n = 800;
            for (int k = 0; k < n; ++k) {
                const double u = ca.step(r, with_last_input(a.plant(), u_last))[0];
                a.step(u);
                u_last = u;
                b.step(std::clamp(cb.step(r, b.plant())[0], -1e6, 1e6));
                ta.add(k, n, a.position() - 0.2);
                tb.add(k, n, b.position() - 0.2);
            }
            CHECK(ta.worst < 1e-6);
            if (theta == kTs) { // at 0.4·Ts the delay-blind MPC still holds; a whole sample it does not
                CHECK(tb.worst > 1e-2);
            }
        }
    }

    TEST_CASE("offset-free MPC on the delayed model rejects a load; the delay-blind one does not") {
        // One sample of input delay, the position read. The load adds to the command, so it is
        // delayed with it; in steady state it acts as Γ·d on the plant either way, so both
        // disturbance models estimate d̂ = the load. At RΔu = 1e-8 the delay-blind MPC moves
        // fast enough that the delay destabilizes it.
        design::MPCWeights<1, 1> w{};
        w.Qy = 1e3;
        w.Rdu = 1e-8;
        const design::KalmanNoise noise{Matrix<2, 2>{{1e-4, 0.0}, {0.0, 1e-2}}, Matrix<1, 1>{{1e-6}}, Matrix<1, 1>{{1e-2}}};
        const auto                naive = design::ofmpc<40, 10>(kUndelayedG, w, {}, noise, design::DisturbanceModel<2, 1, 1>{.Bd = kUndelayedG.B});
        REQUIRE(naive.success);
        const auto run = [&](const auto& aware) {
            REQUIRE(aware.success);
            OFMPC           ca{aware};
            OFMPC           cb{naive};
            Truth<3>        a{kReg1};
            Truth<3>        b{kReg1};
            const ColVec<1> r{{0.2}};
            Tail            ta;
            Tail            tb;
            constexpr int   n = 2000;
            for (int k = 0; k < n; ++k) {
                const double load = (k > n / 4) ? 0.5 : 0.0;
                a.step(ca.step(r, a.y())[0], load);
                b.step(std::clamp(cb.step(r, b.y())[0], -1e6, 1e6), load);
                ta.add(k, n, a.position() - 0.2);
                tb.add(k, n, b.position() - 0.2);
            }
            CHECK(ta.worst < 1e-5);
            CHECK(tb.worst > 1e-2);
            CHECK(ca.disturbance_estimate()(0) == doctest::Approx(0.5).epsilon(1e-3));
        };
        SUBCASE("a load on the plant states, not delayed (B_d zero on the register)") {
            run(design::ofmpc<40, 10>(kDelayedG, w, {}, noise, design::DisturbanceModel<3, 1, 1>{.Bd = kLoadOnPlant}));
        }
        SUBCASE("an input disturbance, delayed with the command through the register") {
            run(design::ofmpc<40, 10>(kDelayedG, w, {}, noise, design::InputDisturbance{}));
        }
    }

    TEST_CASE("the offset-free LQG on the delayed model tracks a moving reference, its feedforward one period ahead") {
        // The register holds the input the plant receives this period, so the reference point
        // for the delayed model is [x*(t); u*(t)] with the input reference u*(t + Ts): the
        // command leads the plant by the delay, and x_ref(k+1) = A·x_ref(k) + B·u_ref(k) holds
        // in the register row. u* is the exact force that moves the mass-spring along
        // x* = 0.2·sin(ωt); the plant holds each sample of it for a period, the residual the
        // stiff loop leaves of that is below a micrometre. At R = 1e-7 the delay-blind law is
        // fast enough that the delay destabilizes it.
        const Matrix<1, 1>        R{{1e-7}};
        const design::KalmanNoise noise{Matrix<2, 2>{{1e-4, 0.0}, {0.0, 1e-2}}, Matrix<1, 1>{{1e-6}}, Matrix<1, 1>{{1e-2}}};
        const auto                aware = design::oflqg(kDelayedG, {.Q = position_weight(1e4), .R = R}, noise, design::DisturbanceModel<3, 1, 1>{.Bd = kLoadOnPlant});
        const auto                naive = design::oflqg(kUndelayedG, {.Q = Matrix<2, 2>{{1e4, 0.0}, {0.0, 0.0}}, .R = R}, noise, design::DisturbanceModel<2, 1, 1>{.Bd = kUndelayedG.B});
        REQUIRE(aware.success);
        REQUIRE(naive.success);

        OFLQG            ca{aware};
        OFLQG            cc{aware}; // the same law, its feedforward not led
        OFLQG            cb{naive};
        Truth<3>         a{kReg1};
        Truth<3>         c{kReg1};
        Truth<3>         b{kReg1};
        constexpr double w = 1.5;
        const auto       x_star = [](double t) { return ColVec<2>{0.2 * std::sin(w * t), 0.2 * w * std::cos(w * t)}; };
        const auto       u_star = [](double t) { return (0.2 * (4.0 - w * w) * std::sin(w * t)) + (0.4 * 0.2 * w * std::cos(w * t)); };
        Tail             ta;
        Tail             tc;
        Tail             tb;
        constexpr int    n = 2000;
        for (int k = 0; k < n; ++k) {
            const double    t = k * kTs;
            const ColVec<2> xs = x_star(t);
            ta.add(k, n, a.position() - xs[0]); // x(k) against x*(t_k)
            tc.add(k, n, c.position() - xs[0]);
            tb.add(k, n, b.position() - xs[0]);
            const ColVec<3> x_ref{xs[0], xs[1], u_star(t)};
            a.step(ca.step({}, x_ref, ColVec<1>{{u_star(t + kTs)}}, a.y())[0], 0.3);
            c.step(cc.step({}, x_ref, ColVec<1>{{u_star(t)}}, c.y())[0], 0.3);
            b.step(std::clamp(cb.step({}, xs, ColVec<1>{{u_star(t)}}, b.y())[0], -1e6, 1e6), 0.3);
        }
        CHECK(ta.worst < 2e-6); // of the 0.2 amplitude, under the load
        CHECK(ta.worst < tc.worst);
        CHECK(ca.disturbance_estimate()(0) == doctest::Approx(0.3).epsilon(1e-3));
        CHECK(tb.worst > 0.1);
    }

    TEST_CASE("MHE on a hand-built delay-augmented model, Q sized to the plant, tracks the true state") {
        const Matrix<2, 2> Qn{{1e-6, 0.0}, {0.0, 1e-4}};
        const auto         art = design::mhe<10>(kReg2, {.process = Qn, .measurement = 1e-6}); // G leaves the registers noise-free
        REQUIRE(art.success);
        MHE<4, 1, 1, 10, double> mhe{art};
        Truth<4>                 t{kReg2};
        Tail                     e;
        constexpr int            n = 1200;
        for (int k = 0; k < n; ++k) {
            const auto x_hat = mhe.step(t.y(), ColVec<1>{{u_at(k - 1)}});
            e.add(k, n, x_hat[1] - t.x[1]);
            t.step(u_at(k));
        }
        CHECK(e.worst < 1e-9);
    }

    TEST_CASE("the measured-state laws take N registers: whole samples, a fraction, spare capacity") {
        check_measured_laws<2>(2.0 * kTs);
        check_measured_laws<5>(5.0 * kTs);
        check_measured_laws<3>(2.4 * kTs);
        check_measured_laws<5>(2.4 * kTs);
    }

    TEST_CASE("the laws that read the output take N registers: whole samples, a fraction, spare capacity") {
        check_output_laws<2>(2.0 * kTs);
        check_output_laws<5>(5.0 * kTs);
        check_output_laws<3>(2.4 * kTs);
        check_output_laws<5>(2.4 * kTs);
    }

    TEST_CASE("the estimators take N registers: whole samples, a fraction, spare capacity") {
        check_estimators<2>(2.0 * kTs);
        check_estimators<5>(5.0 * kTs);
        check_estimators<3>(2.4 * kTs);
        check_estimators<5>(2.4 * kTs);
    }

    TEST_CASE("LQR on N registers is the delay-free LQR on the state the command will meet") {
        // K = K₀·[Aᴺ, B, A·B, …, Aᴺ⁻¹·B]: the registers are newest first, and the newest has the
        // longest to travel, so x(k+N) = Aᴺx + Σⱼ Aʲ⁻¹·B·v(k−j).
        constexpr size_t N = 5;
        const auto       model = discretize(mass_spring(), kTs, InputDelay<N>{N * kTs}).value();
        const auto       free = design::lqr(kUndelayed, {.Q = Matrix<2, 2>{{1e4, 0.0}, {0.0, 0.0}}, .R = Matrix<1, 1>{{1e-4}}});
        const auto       held = design::lqr(model, {.Q = position_weight_on<N>(1e4), .R = Matrix<1, 1>{{1e-4}}});
        REQUIRE(free.success);
        REQUIRE(held.success);
        Matrix<1, 2 + N> expected{};
        Matrix<2, 2>     power = Matrix<2, 2>::identity();
        for (size_t j = 0; j < N; ++j) {
            expected(0, 2 + j) = Matrix<1, 1>(free.K * Matrix<2, 1>(power * kUndelayed.B))(0, 0);
            power = power * kUndelayed.A;
        }
        const Matrix<1, 2> on_plant = free.K * power;
        expected(0, 0) = on_plant(0, 0);
        expected(0, 1) = on_plant(0, 1);
        CHECK((held.K - expected).norm() <= 1e-9 * expected.norm());
    }

    TEST_CASE("a clamped command: the filter's registers hold what back_calculate recorded") {
        constexpr size_t N = 3;
        constexpr size_t NS = 2 + N;
        const auto       model = discretize(mass_spring(), kTs, InputDelay<N>{N * kTs}).value();
        const auto       truth = test::register_form<N>(mass_spring(), kTs, N * kTs);
        const auto       res = design::lqg(model, {.Q = position_weight_on<N>(1e4), .R = Matrix<1, 1>{{1e-4}}}, design::KalmanNoise<2, 1>{.process = Matrix<2, 2>::identity(), .measurement = 1e-6});
        REQUIRE(res.success);
        LQG        law{res};
        ColVec<NS> x0{};
        x0[0] = 1.0;
        Truth<NS>        t{truth, x0};
        AppliedInputs<N> applied;
        double           worst = 0.0;
        int              clamped = 0;
        for (int k = 0; k < 300; ++k) {
            law.update(t.y());
            const double requested = law.control()[0];
            const double u = std::clamp(requested, -0.5, 0.5);
            clamped += (u != requested) ? 1 : 0;
            law.back_calculate(ColVec<1>{{u}});
            law.predict();
            applied.push(u);
            t.step(u);
            for (size_t j = 0; j < N; ++j) {
                worst = std::max(worst, std::abs(law.state_estimate()[2 + j] - applied.u[j]));
            }
        }
        CHECK(clamped > 10);
        CHECK(worst < 1e-12);
    }

    TEST_CASE("a placed observer on N registers: their poles at 0, and every register read") {
        //! Two or more registers repeat the pole at 0, and the placement on (A, C A) then turned on
        //! roundoff: refused at some Ts, and at Ts = 0.01 only under -O3.
        for (const double Ts : {0.001, 0.01, 0.1, 0.5}) {
            CAPTURE(Ts);
            check_placed_observer<1>(Ts);
            check_placed_observer<2>(Ts);
            check_placed_observer<3>(Ts);
            check_placed_observer<5>(Ts);
        }

        // Registers the plant never reads cannot be observed: refused, where the Kalman gain takes them.
        const auto                      spare = discretize(mass_spring(), kTs, InputDelay<3>{1.0 * kTs}).value();
        damp::array<complex<double>, 5> poles{};
        poles[0] = complex<double>{0.9, 0.0};
        poles[1] = complex<double>{0.92, 0.0};
        CHECK(design::luenberger(spare, poles).status == design::DesignStatus::not_placeable);

        // A pole the pinned mode cannot have is still refused, with the fallback in the way.
        const auto                      whole = discretize(mass_spring(), kTs, InputDelay<2>{2.0 * kTs}).value();
        damp::array<complex<double>, 4> slow{};
        slow[0] = complex<double>{0.9, 0.0};
        slow[1] = complex<double>{0.92, 0.0};
        slow[2] = complex<double>{0.5, 0.0};
        slow[3] = complex<double>{0.6, 0.0};
        CHECK(design::luenberger(whole, slow).status == design::DesignStatus::not_placeable);
    }

    TEST_CASE("two inputs on two registers each: offset-free LQG and MPC on the true plant") {
        const StateSpace<2, 2, 2, double, 2, 2> plant{.A = {{-1.0, 0.2}, {0.1, -2.0}}, .B = {{1.0, 0.0}, {0.0, 1.0}}, .C = {{1.0, 0.0}, {0.0, 1.0}}};
        constexpr double                        theta = 2.0 * kTs;
        const auto                              model = discretize(plant, kTs, InputDelay<2>{theta}).value();
        const auto                              truth = test::register_form<2>(plant, kTs, theta);
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(model)>, StateSpace<6, 2, 2, double, 2, 2>>);
        const design::KalmanNoise noise{Matrix<2, 2>::identity() * 1e-4, Matrix<2, 2>::identity() * 1e-6, Matrix<2, 2>::identity() * 1e-2};
        const ColVec<2>           r{0.3, -0.2};
        const ColVec<2>           load{0.4, 0.1};

        const auto run = [&](auto& law) {
            ColVec<6> x{};
            Tail      e;
            for (int k = 0; k < 2000; ++k) {
                const ColVec<2> y = truth.C * x;
                const ColVec<2> u = law.step(r, y);
                e.add(k, 2000, std::abs(y[0] - r[0]) + std::abs(y[1] - r[1]));
                x = truth.A * x + truth.B * ColVec<2>((k > 500) ? ColVec<2>(u + load) : u);
            }
            CHECK(e.worst < 1e-4);
        };

        const auto lq = design::oflqg(model, {.Q = blkdiag(Matrix<2, 2>::identity() * 1e3, Matrix<4, 4>{}), .R = Matrix<2, 2>::identity() * 1e-2}, noise, design::InputDisturbance{});
        REQUIRE(lq.success);
        OFLQG holding{lq};
        run(holding);

        design::MPCWeights<2, 2> w{};
        w.Qy = Matrix<2, 2>::identity() * 1e3;
        w.Qy_terminal = w.Qy;
        w.Rdu = Matrix<2, 2>::identity() * 1e-2;
        const auto mpc = design::ofmpc<30, 5>(model, w, {}, noise, design::InputDisturbance{});
        REQUIRE(mpc.success);
        OFMPC planning{mpc};
        run(planning);
    }

    TEST_CASE("a law on N registers is designed at compile time and deployed in float") {
        static constexpr StateSpace<1, 1, 1, double, 1, 1> lag{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}};
        static constexpr auto                              model = *discretize(lag, 0.05, InputDelay<2>{0.075});
        static constexpr auto                              res = design::lqg(model, {.Q = blkdiag(Matrix<1, 1>{{1.0}}, Matrix<2, 2>{}), .R = 0.1}, {.process = 0.05, .measurement = 0.001});
        static_assert(res.success);
        static constinit LQG law{res.as<float>()};
        CHECK(std::isfinite(law.step(ColVec<1, float>{{1.0F}})[0]));
    }

    TEST_CASE("a disturbance in G is discretized with the plant and kept off the register") {
        auto plant = mass_spring();
        plant.G = Matrix<2, 2>{{0.0, 0.0}, {0.0, 1.0}}; // a force on the mass, the second column
        const auto plain = discretize(plant, kTs).value();
        const auto delayed = discretize(plant, kTs, InputDelay{0.4 * kTs}).value();
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                CHECK(delayed.G(i, j) == ((i < 2) ? plain.G(i, j) : 0.0)); // not delayed with the command
            }
        }
        CHECK(std::abs(plain.G(0, 1) - plain.B(0, 0)) < 1e-15); // the same ZOH as the input
        CHECK(std::abs(plain.G(1, 1) - plain.B(1, 0)) < 1e-15);
    }
}
