// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file test_mpc_feedforward.cpp
 * @brief The MPC runtimes' feedforward calls, the LQ laws' spelling: control(r, x_ref, u_ref) and
 *        step(r, x_ref, u_ref, y) pose the QP about the target for r plus (x_ref, u_ref), the
 *        delay states at u_ref; with a measured disturbance, d after the feedforward.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "damp/backend.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/controllers/oflqg.hpp"
#include "damp/controllers/state_feedback.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/input_limits.hpp"
#include "damp/design/riccati.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "doctest.h"

using namespace damp;

namespace {

constexpr double kTs = 0.1;

/// Two lags in a row, both read: x₁' = −x₁ + u, x₂' = x₁ − x₂; at rest x₂ = x₁ = u.
constexpr StateSpace<2, 1, 2> kLagsC{.A = {{-1.0, 0.0}, {1.0, -1.0}}, .B = {{1.0}, {0.0}}, .C = Matrix<2, 2>::identity()};

/// Sampled at kTs.
const auto kLags = discretize(kLagsC, kTs).value();

/// A load on the second lag, per step: x₂ takes kTs·w more each period.
constexpr Matrix<2, 1> kLoadOnX2{{0.0}, {kTs}};

/// The plant a step on: the lags under @p u, the load @p w on x₂.
ColVec<2> lags_step(const ColVec<2>& x, double u, double w = 0.0) {
    return (kLags.A * x) + (kLags.B * ColVec<1>{u}) + (kLoadOnX2 * ColVec<1>{w});
}

/// Both outputs weighed, the second the more; Ru damps the input.
design::MPCWeights<1, 2> both_weighed() {
    design::MPCWeights<1, 2> w{.Qy = diag(1.0, 10.0), .Rdu = 0.1, .Ru = 0.5};
    w.Qy_terminal = w.Qy;
    return w;
}

/// A move limit: the same in the plant's frame and in one shifted by the planned point.
constexpr design::InputLimits kMoveLim{.du_lim = {.lower = -0.2, .upper = 0.2}};

/// The planned state at tick @p k: a slow swing about (0.4, 0.4), the lags apart.
ColVec<2> planned(size_t k) {
    const double s = std::sin(0.05 * static_cast<double>(k));
    return ColVec<2>{0.4 + (0.3 * s), 0.4 + (0.2 * s)};
}

/// The worst gap between two runs' commands, relative to the command.
struct Gap {
    double worst{0.0};

    void add(double a, double b) { worst = std::max(worst, std::abs(a - b) / (1.0 + std::abs(a))); }
};

/**
 * The feedforward law on the lags against the same design handed the shifted problem by hand: its
 * readings less x_ref, the input it holds less u_ref (back_calculate), its command plus u_ref.
 * @p ref_at gives r at tick k; the worst gap of the commands.
 */
template<typename Law, typename RefAt>
double against_shifted(Law fed, Law hand, RefAt ref_at, double u_ref) {
    ColVec<2> x{};
    double    u_prev = 0.0;
    Gap       gap;
    for (size_t k = 0; k < 200; ++k) {
        const auto      r = ref_at(k);
        const ColVec<2> x_ref = planned(k);
        hand.back_calculate(ColVec<1>{u_prev - u_ref});
        const double u = fed.step(r, x_ref, ColVec<1>{u_ref}, x)(0);
        gap.add(u, hand.step(r, ColVec<2>(x - x_ref))(0) + u_ref);
        x = lags_step(x, u);
        u_prev = u;
    }
    return gap.worst;
}

/// A preview of a ramp from step 40 on each of the NREF held outputs, seen from tick @p k.
template<typename Preview, size_t NREF>
Preview ramp_preview(size_t k) {
    Preview p{};
    for (size_t j = 0; j < p.size(); ++j) {
        const double v = 0.01 * std::max(0.0, static_cast<double>(k + j) - 40.0);
        for (size_t i = 0; i < NREF; ++i) {
            p[j](i) = v;
        }
    }
    return p;
}

} // namespace

TEST_SUITE("MPC feedforward") {
    TEST_CASE("control(r, x_ref, u_ref) is the MPC handed the shifted problem: the readings less x_ref, the input less u_ref") {
        //! The QP is posed about the target for r plus (x_ref, u_ref): the law in the frame of the
        //! planned point, to roundoff. u_ref is the input that holds the point, so the shifted law
        //! holds u − u_ref and keeps its default input target, the target's input.
        const double u_ref = 0.4;
        SUBCASE("TrackedOutputs: r on the first output, the second where the target puts it") {
            const auto art = design::mpc<20, 5>(kLags, both_weighed(), kMoveLim, design::TrackedOutputs{{0}});
            REQUIRE(art.success);
            CHECK(against_shifted(MPC{art}, MPC{art}, [](size_t k) { return (k < 60) ? 0.0 : 0.2; }, u_ref) < 1e-9);
        }
        SUBCASE("TrackedOutputs: r on the second output, the first where the target puts it") {
            const auto art = design::mpc<20, 5>(kLags, both_weighed(), kMoveLim, design::TrackedOutputs{{1}});
            REQUIRE(art.success);
            CHECK(against_shifted(MPC{art}, MPC{art}, [](size_t k) { return (k < 60) ? 0.0 : 0.3; }, u_ref) < 1e-9);
        }
        SUBCASE("a Preview of r: each step's target from the preview, moved by the feedforward") {
            const auto first = design::mpc<20, 5>(kLags, both_weighed(), kMoveLim, design::TrackedOutputs{{0}});
            const auto tracked = design::mpc<20, 5>(kLags, both_weighed(), kMoveLim, design::TrackedOutputs{{1}});
            REQUIRE(first.success);
            REQUIRE(tracked.success);
            using First = decltype(MPC{first});
            using Tracked = decltype(MPC{tracked});
            CHECK(against_shifted(MPC{first}, MPC{first}, ramp_preview<First::Preview, 1>, u_ref) < 1e-9);
            CHECK(against_shifted(MPC{tracked}, MPC{tracked}, ramp_preview<Tracked::Preview, 1>, u_ref) < 1e-9);
        }
    }

    TEST_CASE("an unconstrained MPC with the DARE's terminal weight: control(r, x_ref, u_ref) is the LQR's, about the same point") {
        //! A free move every step and the LQR's cost make the MPC the LQR law. With feedforward both
        //! are u = u_t + u_ref − K·(x − x_t − x_ref), the target for r the same.
        constexpr StateSpace<2, 1, 2> di{.A = {{1.0, kTs}, {0.0, 1.0}}, .B = {{0.5 * kTs * kTs}, {kTs}}, .C = Matrix<2, 2>::identity(), .Ts = kTs};
        const Matrix<2, 2>            Q = diag(1.0, 0.1);
        const Matrix<1, 1>            R{{0.5}};
        const auto                    lq = design::state_feedback(di, {.Q = Q, .R = R}, design::TrackedOutputs{{0}});
        REQUIRE(lq.success);
        design::MPCWeights<1, 2> w{.Qy = Q, .Rdu = Matrix<1, 1>{}, .Ru = R};
        w.Qy_terminal = design::discrete_lqr(di.A, di.B, Q, R).S;
        const auto art = design::mpc<8, 8>(di, w, design::TrackedOutputs{{0}});
        REQUIRE(art.success);

        LQR       lqr{lq};
        MPC       mpc{art};
        ColVec<2> x{1.0, -0.5};
        for (size_t k = 0; k < 40; ++k) {
            const double    r = (k < 20) ? 0.0 : 0.3;
            const ColVec<2> x_ref{0.2 * std::sin(0.2 * static_cast<double>(k)), 0.1};
            const ColVec<1> u_ref{0.05};
            lqr.update(x);
            const double u_lqr = lqr.control(r, x_ref, u_ref)(0);
            const double u_mpc = mpc.step(r, x_ref, u_ref, x)(0);
            CHECK(u_mpc == doctest::Approx(u_lqr).epsilon(1e-6));
            CHECK(mpc.input_target()(0) == doctest::Approx(lq.target.target_r(2, 0) * r + 0.05).epsilon(1e-12)); // Ru's anchor: u_t + u_ref
            x = (di.A * x) + (di.B * ColVec<1>{u_mpc});
        }
    }

    TEST_CASE("an offset-free MPC tracks a ramp given through x_ref, the speed its feedforward, as the OFLQG does") {
        //! A mass with friction, the position read: v' = −c·v + b·u. The ramp at v_ref is a
        //! trajectory of the model under u_ref = c·v_ref/b: handed both, the law follows it with no
        //! lag and d̂ stays 0, the filter carrying the point's motion between ticks.
        constexpr double              c = 1.0;
        constexpr double              b = 2.0;
        constexpr double              v_ref = 0.5;
        constexpr StateSpace<2, 1, 1> cart{.A = {{0.0, 1.0}, {0.0, -c}}, .B = {{0.0}, {b}}, .C = {{1.0, 0.0}}};
        const auto                    truth = discretize(cart, kTs).value();
        const design::KalmanNoise     noise{.process = diag(1e-4, 1e-2), .measurement = Matrix<1, 1>{{1e-6}}, .disturbance = Matrix<1, 1>{{1e-1}}};
        design::MPCWeights<1, 1>      w{.Qy = 100.0, .Rdu = 0.01, .Ru = 0.1};
        w.Qy_terminal = w.Qy;
        const auto ofmpc = design::mpc<20, 5>(cart, kTs, w, noise, design::InputDisturbance{});
        const auto oflqg = design::state_feedback(cart, kTs, {.Q = diag(100.0, 1.0), .R = 0.1}, noise, design::InputDisturbance{});
        REQUIRE(ofmpc.success);
        REQUIRE(oflqg.success);

        struct Run {
            double lag{0.0};   ///< The worst |p − p_ref| over the last third
            double d_hat{0.0}; ///< d̂ at the end
        };
        const auto run = [&](auto law, double u_ff, bool fed) {
            ColVec<2> x{};
            Run       out{};
            for (size_t k = 0; k < 600; ++k) {
                const double p_ref = v_ref * kTs * static_cast<double>(k);
                if (k >= 400) {
                    out.lag = std::max(out.lag, std::abs(x(0) - p_ref));
                }
                const ColVec<1> y{x(0)};
                const double    u = fed ? law.step(0.0, ColVec<2>{p_ref, v_ref}, ColVec<1>{u_ff}, y)(0) : law.step(p_ref, y)(0);
                x = (truth.A * x) + (truth.B * ColVec<1>{u});
            }
            out.d_hat = law.disturbance_estimate()(0);
            return out;
        };
        const double u_hold = c * v_ref / b;
        const Run    mpc_fed = run(OFMPC{ofmpc}, u_hold, true);
        const Run    lqg_fed = run(OFLQG{oflqg}, u_hold, true);
        CHECK(mpc_fed.lag < 1e-9);
        CHECK(std::abs(mpc_fed.d_hat) < 1e-9);
        CHECK(lqg_fed.lag < 1e-9);
        CHECK(std::abs(lqg_fed.d_hat) < 1e-9);

        //! The plant carries no load, so the filter finds none either way. Without u_ref the model
        //! carries the point under no input, slowing it, and both laws lag the ramp; with r alone
        //! the target is a standstill at p_ref, and they lag more.
        const Run mpc_speed_only = run(OFMPC{ofmpc}, 0.0, true);
        const Run lqg_speed_only = run(OFLQG{oflqg}, 0.0, true);
        CHECK(mpc_speed_only.lag > 1e-2);
        CHECK(lqg_speed_only.lag > 1e-2);
        CHECK(std::abs(mpc_speed_only.d_hat) < 1e-9);
        CHECK(run(OFMPC{ofmpc}, 0.0, false).lag > 5.0 * mpc_speed_only.lag);
        CHECK(run(OFLQG{oflqg}, 0.0, false).lag > 5.0 * lqg_speed_only.lag);
    }

    TEST_CASE("a measured-state MPC with an InputDelay pads x_ref with u_ref on the delay states it holds") {
        //! At the planned point each input on its way is u_ref, as the LQ laws pad it. The shifted
        //! law holds its own inputs, u − u_ref: the same QP once its delay states start there.
        const auto art = design::mpc<20, 5>(kLagsC, kTs, both_weighed(), kMoveLim, design::TrackedOutputs{{1}}, InputDelay<2>{1.5 * kTs});
        REQUIRE(art.success);
        const auto   truth = discretize(kLagsC, kTs, InputDelay<2>{1.5 * kTs}).value(); // [x; v(k−1); v(k−2)]
        const double u_ref = 0.4;
        MPC          fed{art};
        MPC          hand{art};
        for (size_t i = 0; i < 2; ++i) {
            hand.back_calculate(ColVec<1>{-u_ref}); // its delay states at its own frame's zero input
            hand.predict();
        }
        ColVec<4> x{};
        double    u_prev = 0.0;
        Gap       gap;
        for (size_t k = 0; k < 200; ++k) {
            const ColVec<2> plant{x(0), x(1)};
            const ColVec<2> x_ref = planned(k);
            const double    r = (k < 60) ? 0.0 : 0.3;
            hand.back_calculate(ColVec<1>{u_prev - u_ref});
            const double u = fed.step(r, x_ref, ColVec<1>{u_ref}, plant)(0);
            gap.add(u, hand.step(r, plant - x_ref)(0) + u_ref);
            x = (truth.A * x) + (truth.B * ColVec<1>{u});
            u_prev = u;
        }
        CHECK(gap.worst < 1e-9);
    }

    TEST_CASE("an OFMPC pads them through its filter: an InputDelay's delay states at u_ref, on a late reading") {
        //! At a point the model holds, x_ref = (a, a) under u_ref = a, the law is the shifted one:
        //! readings less C·x_ref, the input less u_ref, both filters started at the point. The QP
        //! reads the input's delay states; the late readings' are the filter's alone, at 0 in the
        //! point as the LQ laws leave them, the anchor's motion carrying the reading of x_ref.
        constexpr StateSpace<2, 1, 2, double, 2, 2> noisy{.A = kLagsC.A, .B = kLagsC.B, .C = kLagsC.C};
        const design::KalmanNoise                   noise{.process = diag(1e-3, 1e-3), .measurement = diag(1e-6, 1e-6), .disturbance = Matrix<1, 1>{{1e-1}}};
        const design::DisturbanceModel              load{.Bd = Matrix<2, 1>{{0.0}, {1.0}}, .outputs = {1}};
        const InputDelay                            in{0.4 * kTs};
        const OutputDelay<2>                        out{1.4 * kTs};
        const auto                                  art = design::mpc<20, 5>(noisy, kTs, both_weighed(), kMoveLim, noise, load, in, out);
        REQUIRE(art.success);
        auto truth = discretize(noisy, kTs, in, out).value(); // [x; v(k−1); the readings on their way]
        static_assert(decltype(truth.A)::rows() == 7);

        const double    a = 0.6;
        const ColVec<2> x_ref{a, a};
        const ColVec<1> u_ref{a};
        OFMPC           fed{art};
        OFMPC           hand{art};
        ColVec<7>       x{};
        for (size_t i = 0; i < 7; ++i) {
            x(i) = a; // at rest at the point: each lag, the input on its way and every reading at a
        }
        fed.estimator.reset(x);
        fed.back_calculate(u_ref);
        double u_prev = a;
        Gap    gap;
        for (size_t k = 0; k < 300; ++k) {
            const double    r = (k < 50) ? 0.0 : 0.2;
            const double    w = (k < 150) ? 0.0 : 0.3;
            const ColVec<2> y = truth.C * x;
            hand.back_calculate(ColVec<1>{u_prev - a});
            const double u = fed.step(r, x_ref, u_ref, y)(0);
            gap.add(u, hand.step(r, y - x_ref)(0) + a);
            x = (truth.A * x) + (truth.B * ColVec<1>{u});
            x(1) += kTs * w; // a load on the second lag
            u_prev = u;
        }
        CHECK(gap.worst < 1e-9);
        CHECK(fed.disturbance_estimate()(0) == doctest::Approx(hand.disturbance_estimate()(0)).epsilon(1e-9));
    }

    TEST_CASE("with a measured disturbance: control(r, x_ref, u_ref, d) and step(r, x_ref, u_ref, d, x), the reference and its feedforward, then what was read") {
        const auto art = design::mpc<20, 5>(kLags, both_weighed(), kMoveLim, design::MeasuredDisturbance{kLoadOnX2}, design::TrackedOutputs{{1}});
        REQUIRE(art.success);
        const double u_ref = 0.4;
        MPC          fed{art};
        MPC          split{art};
        MPC          hand{art};
        ColVec<2>    x{};
        double       u_prev = 0.0;
        Gap          gap;
        Gap          ticks;
        for (size_t k = 0; k < 200; ++k) {
            const ColVec<1> d{(k < 80) ? 0.0 : 0.4};
            const ColVec<2> x_ref = planned(k);
            hand.back_calculate(ColVec<1>{u_prev - u_ref});
            const double u = fed.step(0.3, x_ref, ColVec<1>{u_ref}, d, x)(0);
            split.predict();
            split.update(x);
            ticks.add(u, split.control(0.3, x_ref, ColVec<1>{u_ref}, d)(0));
            gap.add(u, hand.step(0.3, d, x - x_ref)(0) + u_ref);
            x = lags_step(x, u, d(0));
            u_prev = u;
        }
        CHECK(gap.worst < 1e-9);
        CHECK(ticks.worst < 1e-12); // step is predict, update, control
    }

    TEST_CASE("the output limits stay absolute: x_ref moves the cost's point, not the bounds") {
        //! The planned point puts x₂ at 1; a hard x₂ ≤ 0.6 holds it at 0.6, on the output itself, as
        //! MATLAB's limits are on the plant's outputs. A bound that rides the point is moved at run time.
        const auto bounded = design::mpc<20, 5>(kLags, both_weighed(), design::TrackedOutputs{{1}}, design::OutputLimits{.outputs = {1}, .lim = {.upper = 0.6}, .upper_ecr = {0.0}});
        const auto unbounded = design::mpc<20, 5>(kLags, both_weighed(), design::TrackedOutputs{{1}});
        REQUIRE(bounded.success);
        REQUIRE(unbounded.success);
        const auto settle = [](auto law) {
            ColVec<2> x{};
            double    peak = 0.0;
            for (size_t k = 0; k < 300; ++k) {
                x = lags_step(x, law.step(0.0, ColVec<2>{1.0, 1.0}, ColVec<1>{1.0}, x)(0));
                peak = std::max(peak, x(1));
            }
            return ColVec<2>{x(1), peak};
        };
        const ColVec<2> b = settle(MPC{bounded});
        CHECK(b(1) <= 0.6 + 1e-9);
        CHECK(b(0) == doctest::Approx(0.6).epsilon(1e-6));
        CHECK(settle(MPC{unbounded})(0) == doctest::Approx(1.0).epsilon(1e-9)); // unbounded: at the point
    }
}
