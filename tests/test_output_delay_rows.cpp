// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// OutputDelay on some outputs only (a sensor read over a field bus next to a local one): the
// register is exact, bad selections are refused, and every design family works on the result
// the way it does on a model whose outputs are all late.

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
#include "damp/estimation/kalman.hpp"
#include "damp/estimation/mhe.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "doctest.h"

using namespace damp;

namespace {

constexpr double kTs = 0.01;

/// Two masses on springs, coupled: [x₁, v₁, x₂, v₂], a force on each, both positions measured.
/// y₀ = x₁ is local; y₁ = x₂ comes over a bus (the late one). Process noise on every state.
constexpr StateSpace<4, 2, 2, double, 4, 2> two_masses() {
    StateSpace<4, 2, 2, double, 4, 2> s{};
    s.A(0, 1) = 1.0;
    s.A(1, 0) = -5.0;
    s.A(1, 1) = -0.4;
    s.A(1, 2) = 4.0;
    s.A(2, 3) = 1.0;
    s.A(3, 0) = 4.0;
    s.A(3, 2) = -6.0;
    s.A(3, 3) = -0.3;
    s.B(1, 0) = 1.0;
    s.B(3, 1) = 1.0;
    s.C(0, 0) = 1.0;
    s.C(1, 2) = 1.0;
    return s; // G = I, H = I by default
}

ColVec<2> u_at(int k) {
    return (k < 0) ? ColVec<2>{} : ColVec<2>{std::sin(0.037 * k) + 0.3, std::cos(0.021 * k) - ((k / 150) % 2 == 0 ? 0.4 : -0.4)};
}

const auto kPlain = discretize(two_masses(), kTs).value();
const auto kLate = discretize(two_masses(), kTs, OutputDelay<2, 1>{{1}}).value(); // y₁ two samples late
const auto kBoth = discretize(two_masses(), kTs, InputDelay<2>{1.4 * kTs}, OutputDelay<2, 1>{{1}}).value();

/// The true plant: the delayed model itself (exact).
template<size_t N>
struct Truth {
    StateSpace<N, 2, 2, double, 4, 2> m{};
    ColVec<N>                         x{};
    [[nodiscard]] ColVec<2>           y() const { return m.C * x; }
    void                              step(const ColVec<2>& u, const ColVec<2>& load = ColVec<2>{}) { x = m.A * x + m.B * (u + load); }
};

/// Largest |v| over the last quarter of a run.
struct Tail {
    double worst = 0.0;
    void   add(int k, int n, double v) {
        if (4 * k >= 3 * n) {
            worst = std::max(worst, std::abs(v));
        }
    }
};

} // namespace

TEST_SUITE("Output delay on selected outputs") {
    TEST_CASE("only the listed outputs are late; the others stay current, feedthrough included") {
        auto plant = two_masses();
        plant.D(0, 0) = 0.5;
        plant.D(1, 1) = 0.3;
        const auto base = discretize(plant, kTs).value();
        const auto late = discretize(plant, kTs, OutputDelay<2, 1>{{1}}).value();
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(late)>, StateSpace<6, 2, 2, double, 4, 2>>); // one row × two samples
        ColVec<4> x{};
        ColVec<6> xa{};
        double    y1_hist[300]{};
        double    worst = 0.0;
        for (int k = 0; k < 300; ++k) {
            const ColVec<2> u = u_at(k);
            const ColVec<2> y = base.C * x + base.D * u;
            const ColVec<2> ya = late.C * xa + late.D * u;
            y1_hist[k] = y[1];
            worst = std::max({worst, std::abs(ya[0] - y[0]), std::abs(ya[1] - ((k >= 2) ? y1_hist[k - 2] : 0.0))});
            x = base.A * x + base.B * u;
            xa = late.A * xa + late.B * u;
        }
        CHECK(worst < 1e-12);
        for (size_t r = 4; r < 6; ++r) { // the register is deterministic
            for (size_t j = 0; j < 4; ++j) {
                CHECK(late.G(r, j) == 0.0);
            }
        }
        CHECK(late.D(1, 1) == 0.0); // the late reading's feedthrough went into the register
        CHECK(late.D(0, 0) == 0.5); // the current one keeps it
    }

    TEST_CASE("listing every output is the all-outputs delay; bad lists are refused") {
        const auto all = discretize(two_masses(), kTs, OutputDelay<1>{}).value();
        const auto listed = discretize(two_masses(), kTs, OutputDelay<1, 2>{{0, 1}}).value();
        CHECK((all.A - listed.A).norm() < 1e-12);
        CHECK((all.B - listed.B).norm() < 1e-12);
        CHECK((all.C - listed.C).norm() < 1e-12);
        CHECK((all.D - listed.D).norm() < 1e-12);
        CHECK_FALSE(discretize(two_masses(), kTs, OutputDelay<1, 2>{{1, 1}}).has_value()); // listed twice
        CHECK_FALSE(discretize(two_masses(), kTs, OutputDelay<1, 1>{{2}}).has_value());    // no such output
        CHECK_FALSE(discretize(two_masses(), kTs, InputDelay<1>{}, OutputDelay<1, 1>{{5}}).has_value());
        CHECK((discretize(two_masses(), kTs, OutputDelay<0, 1>{{1}}).value().A - kPlain.A).norm() < 1e-12); // no delay: the plain model
        CHECK(discretize(kPlain, kTs, OutputDelay<1, 1>{{1}}).has_value());                                 // an already-discrete model, same Ts
    }

    TEST_CASE("with an input delay: the input-delayed model, its listed output late") {
        const auto in_only = discretize(two_masses(), kTs, InputDelay<2>{1.4 * kTs}).value();
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(kBoth)>, StateSpace<10, 2, 2, double, 4, 2>>);
        ColVec<8>  x{};
        ColVec<10> xb{};
        double     hist[300]{};
        double     worst = 0.0;
        for (int k = 0; k < 300; ++k) {
            const ColVec<2> y = in_only.C * x;
            const ColVec<2> yb = kBoth.C * xb;
            hist[k] = y[1];
            worst = std::max({worst, std::abs(yb[0] - y[0]), std::abs(yb[1] - ((k >= 2) ? hist[k - 2] : 0.0))});
            x = in_only.A * x + in_only.B * u_at(k);
            xb = kBoth.A * xb + kBoth.B * u_at(k);
        }
        CHECK(worst < 1e-12);
    }

    TEST_CASE("Kalman and MHE reconstruct the current state; the undelayed Kalman lags") {
        const Matrix<4, 4>                             Qn = Matrix<4, 4>::identity() * 1e-5;
        const Matrix<2, 2>                             Rn = Matrix<2, 2>::identity() * 1e-6;
        SteadyStateKalmanFilter<6, 2, 2, double, 4, 2> aware{design::kalman(kLate, Qn, Rn)};
        SteadyStateKalmanFilter<4, 2, 2, double, 4, 2> naive{design::kalman(kPlain, Qn, Rn)};
        const auto                                     mhe_art = design::mhe<10>(kLate, Qn, Rn); // Q sized to the plant
        REQUIRE(mhe_art.success);
        MHE<6, 2, 2, 10, double> mhe{mhe_art};
        Truth<6>                 t{kLate};
        Tail                     ea;
        Tail                     eb;
        Tail                     em;
        constexpr int            n = 1200;
        for (int k = 0; k < n; ++k) {
            const ColVec<2> y = t.y();
            aware.update(y);
            naive.update(y);
            const auto x_mhe = mhe.update(y, u_at(k - 1));
            ea.add(k, n, aware.state()[2] - t.x[2]); // the remote mass, from a stale reading
            eb.add(k, n, naive.state()[2] - t.x[2]);
            em.add(k, n, x_mhe[2] - t.x[2]);
            aware.predict(u_at(k));
            naive.predict(u_at(k));
            t.step(u_at(k));
        }
        CHECK(ea.worst < 1e-9);
        CHECK(em.worst < 1e-8);
        CHECK(eb.worst > 1e-4);
    }

    TEST_CASE("state-feedback families on the delayed model, from plant-sized weights") {
        Matrix<4, 4> Q{};
        Q(0, 0) = 1e3;
        Q(2, 2) = 1e3;
        const Matrix<2, 2> R = Matrix<2, 2>::identity() * 1e-3;
        const auto         lqr = design::discrete_lqr(kBoth.A, kBoth.B, Q, R);
        REQUIRE(lqr.success);
        CHECK(lqr.is_stable());
        Matrix<6, 6> Qi{}; // [x; ∫e], plant-sized
        Qi(0, 0) = 1e3;
        Qi(2, 2) = 1e3;
        Qi(4, 4) = 1e2;
        Qi(5, 5) = 1e2;
        const auto lqi = design::discrete_lqi(kBoth, Qi, R);
        REQUIRE(lqi.success);
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 8; j < 10; ++j) { // the past readings steer nothing in LQR (LQI integrates them)
                CHECK(std::abs(lqr.K(i, j)) < 1e-9 * std::abs(lqr.K(i, 0)));
            }
        }
        LQI<10, 2, 2, double> lqi_c{lqi};
        Truth<10>             ti{kBoth};
        Tail                  ei;
        const ColVec<2>       r{0.2, -0.1};
        constexpr int         n = 2500;
        for (int k = 0; k < n; ++k) {
            ti.step(lqi_c.control(r, ti.y(), ti.x), ColVec<2>{0.3, -0.2}); // under a load: the integrators take it
            ei.add(k, n, std::max(std::abs(ti.x[0] - r[0]), std::abs(ti.x[2] - r[1])));
        }
        CHECK(ei.worst < 1e-4);
        design::MPCWeights<2, 2> w{};
        w.Qy = Matrix<2, 2>::identity() * 1e4;
        w.Rdu = Matrix<2, 2>::identity() * 1e-4;
        const auto mpc_art = design::state_mpc<40, 10>(kBoth, w);
        REQUIRE(mpc_art.success);
        MPC       mpc{mpc_art};
        Truth<10> t{kBoth};
        Tail      e;
        for (int k = 0; k < n; ++k) {
            t.step(mpc.control(r, t.x));
            e.add(k, n, std::max(std::abs(t.x[0] - r[0]), std::abs(t.x[2] - r[1])));
        }
        CHECK(e.worst < 1e-6);
    }

    TEST_CASE("output-feedback families hold the loop through both delays") {
        Matrix<4, 4> Q{};
        Q(0, 0) = 1e3;
        Q(2, 2) = 1e3;
        const Matrix<2, 2> R = Matrix<2, 2>::identity() * 1e-3;
        const Matrix<4, 4> Qn = Matrix<4, 4>::identity() * 1e-4;
        const Matrix<2, 2> Rn = Matrix<2, 2>::identity() * 1e-6;
        const ColVec<2>    load{0.3, -0.2};
        constexpr int      n = 2500;

        SUBCASE("LQG regulates") {
            const auto art = design::discrete_lqg(kBoth, Q, R, {.process = Qn, .measurement = Rn});
            REQUIRE(art.success);
            LQG<10, 2, 2, double, 4, 2> c{art};
            Truth<10>                   t{kBoth, ColVec<10>{0.3, 0.0, -0.2}};
            Tail                        e;
            for (int k = 0; k < n; ++k) {
                t.step(c.step(t.y()));
                e.add(k, n, std::max(std::abs(t.x[0]), std::abs(t.x[2])));
            }
            CHECK(e.worst < 1e-6);
        }
        SUBCASE("LQGI tracks a step under a load, from a plant-state reference") {
            Matrix<6, 6> Qi{};
            Qi(0, 0) = 1e3;
            Qi(2, 2) = 1e3;
            Qi(4, 4) = 1e2;
            Qi(5, 5) = 1e2;
            const auto art = design::discrete_lqgi(kBoth, Qi, R, {.process = Qn, .measurement = Rn});
            REQUIRE(art.success);
            LQGI<10, 2, 2, double, 4, 2> c{art};
            Truth<10>                    t{kBoth};
            const ColVec<2>              r{0.2, -0.1};
            const ColVec<4>              x_ref{0.2, 0.0, -0.1, 0.0};
            Tail                         e;
            for (int k = 0; k < n; ++k) {
                t.step(c.step(r, t.y(), x_ref), load);
                e.add(k, n, std::max(std::abs(t.x[0] - r[0]), std::abs(t.x[2] - r[1])));
            }
            CHECK(e.worst < 1e-4);
        }
        SUBCASE("OFLQG tracks a moving reference under a load") {
            design::OffsetFreeModel<4, 2, 2, 2, 2> model{};
            model.Bd = kPlain.B; // a load at the plant input, after the delay
            model.Cz = kPlain.C;
            const auto art = design::discrete_oflqg(kBoth, model, Q, R, design::OffsetFreeNoise{Qn, Rn, Matrix<2, 2>::identity() * 1e-2});
            REQUIRE(art.success);
            OFLQG<10, 2, 2, double, 2, 2> c{art};
            Truth<10>                     t{kBoth};
            Tail                          e;
            for (int k = 0; k < n; ++k) {
                const double    s = 0.1 * std::sin(0.5 * k * kTs);
                const ColVec<4> x_ref{s, 0.05 * std::cos(0.5 * k * kTs), -s, -0.05 * std::cos(0.5 * k * kTs)};
                const ColVec<2> u_ref = ColVec<2>{(5.0 - 4.0) * s, (6.0 - 4.0) * -s}; // quasi-static: the springs' force
                t.step(c.control(t.y(), x_ref, u_ref), load);
                e.add(k, n, std::max(std::abs(t.x[0] - x_ref[0]), std::abs(t.x[2] - x_ref[2])));
            }
            CHECK(e.worst < 5e-3);
        }
        SUBCASE("OffsetFreeMPC rejects the load") {
            design::MPCWeights<2, 2> w{};
            w.Qy = Matrix<2, 2>::identity() * 1e3;
            w.Rdu = Matrix<2, 2>::identity() * 1e-6;
            const auto art = design::mpc<40, 10>(kBoth, w, {}, design::OffsetFreeNoise{Qn, Rn, Matrix<2, 2>::identity() * 1e-2}); // Q_process, R_measurement, Q_disturbance
            REQUIRE(art.success);
            OffsetFreeMPC   c{art};
            Truth<10>       t{kBoth};
            const ColVec<2> r{0.2, -0.1};
            Tail            e;
            for (int k = 0; k < n; ++k) {
                t.step(c.control(r, t.y()), (k > n / 4) ? load : ColVec<2>{});
                e.add(k, n, std::max(std::abs(t.x[0] - r[0]), std::abs(t.x[2] - r[1])));
            }
            CHECK(e.worst < 1e-5);
        }
    }
}
