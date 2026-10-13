// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file test_delay_predictor.cpp
 * @brief The predictor form of an input delay equals the register form: OFLQG and OffsetFreeMPC
 *        designed on the short model and run through InputDelayPredictor give the same commands,
 *        to roundoff, as the same designs on discretize(sys, Ts, InputDelay<D>).
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <type_traits>

#include "damp/controllers/delay_predictor.hpp"
#include "damp/controllers/offset_free_lqg.hpp"
#include "damp/controllers/offset_free_mpc.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "doctest.h"

using namespace damp;

namespace {

/// The heater of examples/control/heater_ofmpc: K = 6 °C/W, τ = 30 s, 5 s dead time at Ts = 2 s.
constexpr double kTs = 2.0;
constexpr double kDead = 5.0;

constexpr StateSpace<1, 1, 1> heater_block() {
    return StateSpace<1, 1, 1>{.A = {{-1.0 / 30.0}}, .B = {{6.0 / 30.0}}, .C = {{1.0}}};
}

/// Register form: the whole delay in the model (NX = 1 + 3).
constexpr auto kRegister = discretize(heater_block(), kTs, InputDelay<3>{kDead});
/// Predictor form: the fractional second as one register state (NX = 2), two whole samples predicted.
constexpr auto kSplit = split_input_delay(kDead, kTs);
constexpr auto kShort = discretize(heater_block(), kTs, InputDelay<1>{kSplit.fraction});

/// The plant with the true delay: the register model itself, driven by u + w (an input disturbance).
struct Plant {
    ColVec<4> x{};
    double    y() const { return (kRegister->C * x)(0, 0); }
    void      step(double u, double w) { x = kRegister->A * x + kRegister->B * ColVec<1>{u + w}; }
};

} // namespace

TEST_SUITE("Input delay predictor") {
    TEST_CASE("split_input_delay: whole samples and the remainder") {
        static_assert(kSplit.whole == 2);
        static_assert(kSplit.fraction == 1.0);
        static_assert(kShort.has_value());
        static_assert(kRegister.has_value());
        constexpr auto whole = split_input_delay(6.0, kTs);
        static_assert(whole.whole == 3);
        static_assert(whole.fraction == 0.0);
        constexpr auto none = split_input_delay(0.0, kTs);
        static_assert(none.whole == 0);
        constexpr auto can = split_input_delay(1.4e-3, 1e-3);
        static_assert(can.whole == 1);
        CHECK(can.fraction == doctest::Approx(0.4e-3));
    }

    TEST_CASE("PredictedInputDelay: one register state for the remainder, none for a whole delay") {
        // 5 s at 2 s: the same model as InputDelay<1>{1.0}
        constexpr auto pred = discretize(heater_block(), kTs, PredictedInputDelay<2>{kDead});
        static_assert(pred.has_value());
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                CHECK(pred->A(i, j) == kShort->A(i, j));
            }
            CHECK(pred->B(i, 0) == kShort->B(i, 0));
        }
        // 6 s at 2 s: the tap carries nothing
        constexpr auto whole = discretize(heater_block(), kTs, PredictedInputDelay<3>{6.0});
        constexpr auto plain = discretize(heater_block(), kTs);
        static_assert(whole.has_value());
        CHECK(whole->A(0, 0) == plain->A(0, 0));
        CHECK(whole->A(0, 1) == 0.0); // Γ₁ = 0
        CHECK(whole->B(0, 0) == plain->B(0, 0));
        CHECK(whole->B(1, 0) == 1.0);
        CHECK_FALSE(discretize(heater_block(), kTs, PredictedInputDelay<0>{-1.0}).has_value());
        // D is checked against the seconds; no seconds is exactly D samples
        static_assert(!discretize(heater_block(), kTs, PredictedInputDelay<3>{kDead}).has_value());
        constexpr auto three = discretize(heater_block(), kTs, PredictedInputDelay<3>{});
        static_assert(three.has_value());
        CHECK(three->B(0, 0) == whole->B(0, 0));
    }

    TEST_CASE("design from the continuous plant and the delay in periods: the discretize route, scalar ticks") {
        constexpr design::OffsetFreeNoise<1, 1, 1> noise{.process = {{0.01}}, .measurement = {{0.01}}, .disturbance = {{0.05}}};
        constexpr design::MPCWeights<1, 1>         w{.Rdu = {{10.0}}};
        constexpr design::MPCConstraints<1, 1>     c{.u_min = {0.0}, .u_max = {40.0}};
        // 5 s at 2 s: 2.5 periods, two predicted and half a period as one register state
        constexpr auto mpc_c = design::mpc<30, 4, kDead / kTs>(heater_block(), kTs, w, c, noise);
        constexpr auto mpc_d = design::mpc<30, 4, 2.0>(*kShort, w, c, noise);
        constexpr auto lqg_c = design::discrete_oflqg<kDead / kTs>(heater_block(), kTs, Matrix<1, 1>{{1.0}}, Matrix<1, 1>{{0.02}}, noise);
        constexpr auto lqg_d = design::discrete_oflqg<2.0>(*kShort, Matrix<1, 1>{{1.0}}, Matrix<1, 1>{{0.02}}, noise);
        static_assert(mpc_c.success && mpc_d.success && lqg_c.success && lqg_d.success);
        static_assert(std::is_same_v<decltype(OffsetFreeMPC{mpc_c}), decltype(OffsetFreeMPC{mpc_d})>);
        static_assert(std::is_same_v<decltype(OFLQG{lqg_c}), OFLQG<2, 1, 1, double, 1, 1, 2>>);
        // a whole number of periods needs no register state; no delay is the plain discretization
        static_assert(std::is_same_v<decltype(OffsetFreeMPC{design::mpc<30, 4, 3.0>(heater_block(), kTs)}), OffsetFreeMPC<1, 1, 1, 30, 4, double, 1, 1, 3>>);
        static_assert(std::is_same_v<decltype(OFLQG{design::discrete_oflqg(heater_block(), kTs, Matrix<1, 1>{{1.0}}, Matrix<1, 1>{{0.02}})}), OFLQG<1, 1, 1, double, 1, 1, 0>>);

        OffsetFreeMPC mc{mpc_c};
        OffsetFreeMPC md{mpc_d};
        OFLQG         lc{lqg_c};
        OFLQG         ld{lqg_d};
        for (int k = 0; k < 40; ++k) {
            const double y = 0.5 * k;
            CHECK(mc.control(185.0, y) == md.control(ColVec<1>{185.0}, ColVec<1>{y})[0]);
            const double u = lc.feedback(185.0, y);
            CHECK(u == ld.feedback(ColVec<1>{185.0}, ColVec<1>{y})[0]);
            lc.commit(damp::clamp(u, 0.0, 40.0));
            ld.commit(ColVec<1>{damp::clamp(u, 0.0, 40.0)});
        }
    }

    TEST_CASE("predict is the closed form through the buffered commands") {
        const Matrix<2, 2>                   A{{0.9, 0.1}, {0.0, 0.8}};
        const Matrix<2, 1>                   B{{0.0}, {1.0}};
        InputDelayPredictor<2, 1, 3, double> p{StateSpace<2, 1, 1>{.A = A, .B = B}};
        // commands u_0 … u_3 pushed; at tick 4 the buffer holds u_0 … u_3 and predict uses u_1 … u_3
        for (int k = 0; k < 4; ++k) {
            p.push(ColVec<1>{1.0 + k});
        }
        CHECK(p.applied_last()(0) == 1.0); // u_0 = u_{k-D-1}
        CHECK(p.applied_now()(0) == 2.0);  // u_1 = u_{k-D}
        const ColVec<2> x0{1.0, -0.5};
        ColVec<2>       expect = x0;
        for (int j = 1; j <= 3; ++j) {
            expect = A * expect + B * ColVec<1>{1.0 + j};
        }
        const auto got = p.predict(x0);
        CHECK(got(0) == doctest::Approx(expect(0)));
        CHECK(got(1) == doctest::Approx(expect(1)));
    }

    TEST_CASE("OFLQG: predictor form equals the register form to roundoff") {
        // the same weights, sized to the plant
        const auto reg = design::discrete_oflqg(
            *kRegister, Matrix<4, 4>{{1.0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}, Matrix<1, 1>{{0.02}},
            design::OffsetFreeNoise{Matrix<4, 4>{{0.01, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}, Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{0.05}}}
        );
        const auto pre = design::discrete_oflqg<2.0>(
            *kShort, Matrix<2, 2>{{1.0, 0.0}, {0.0, 0.0}}, Matrix<1, 1>{{0.02}}, design::OffsetFreeNoise{Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.0}}, Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{0.05}}}
        );
        REQUIRE(reg.success);
        REQUIRE(pre.success);
        OFLQG law_reg{reg};
        OFLQG law_pre{pre};
        static_assert(std::is_same_v<decltype(law_pre), OFLQG<2, 1, 1, double, 1, 1, 2>>); // D from the result

        Plant  a;
        Plant  b;
        double worst = 0.0;
        double scale = 0.0;
        for (int k = 0; k < 400; ++k) {
            const double w = (k >= 200) ? 2.0 : 0.0; // an unmeasured 2 W loss from 400 s
            const double ua = law_reg.control(ColVec<1>{100.0}, ColVec<1>{a.y()})[0];
            const double ub = law_pre.control(ColVec<1>{100.0}, ColVec<1>{b.y()})[0];
            worst = std::max(worst, std::abs(ua - ub));
            scale = std::max(scale, std::abs(ua));
            a.step(ua, w);
            b.step(ub, w);
        }
        CHECK(worst <= 1e-8 * scale);
        CHECK(std::abs(a.y() - 100.0) < 1e-3); // and both are offset-free under the loss
    }

    TEST_CASE("OffsetFreeMPC: predictor form equals the register form, constraints and preview included") {
        // The register form spends its first D prediction steps on outputs no decision can
        // change, so its usable horizon is NP − D; the predictor plans from the state the
        // command lands on and uses all NP. The two are the same QP when the register horizon
        // is D steps longer: 32 against 30 here.
        design::MPCWeights<1, 1> w{};
        w.Qy = Matrix<1, 1>{{1.0}};
        w.Rdu = Matrix<1, 1>{{10.0}};
        w.Qy_terminal = w.Qy;
        design::MPCConstraints<1, 1> c{};
        c.u_min = ColVec<1>{0.0};
        c.u_max = ColVec<1>{40.0};
        const auto reg = design::mpc<32, 4>(*kRegister, w, c, design::OffsetFreeNoise{Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{0.05}}});
        const auto pre = design::mpc<30, 4, 2.0>(*kShort, w, c, design::OffsetFreeNoise{Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{0.05}}});
        REQUIRE(reg.success);
        REQUIRE(pre.success);
        OffsetFreeMPC law_reg{reg};
        OffsetFreeMPC law_pre{pre};
        static_assert(std::is_same_v<decltype(law_pre), OffsetFreeMPC<2, 1, 1, 30, 4, double, 1, 1, 2>>); // D from the result

        Plant  a;
        Plant  b;
        double worst = 0.0;
        int    saturated = 0;
        for (int k = 0; k < 400; ++k) {
            const double w_in = (k >= 200) ? 2.0 : 0.0;
            const double ua = law_reg.control(ColVec<1>{185.0}, ColVec<1>{a.y()})[0]; // the heater's 210 °C setpoint above ambient
            const double ub = law_pre.control(ColVec<1>{185.0}, ColVec<1>{b.y()})[0];
            worst = std::max(worst, std::abs(ua - ub));
            saturated += (ua >= 40.0 - 1e-9) ? 1 : 0;
            a.step(ua, w_in);
            b.step(ub, w_in);
        }
        CHECK(saturated > 10); // the heat-up runs on the clamp, so the QP's constraints were active
        CHECK(worst <= 1e-6 * 40.0);
        CHECK(std::abs(a.y() - 185.0) < 1e-3);

        // With a preview the two must still agree: the predictor form sees the reference D
        // periods further on, since its plan starts where the command lands.
        Plant      pc;
        Plant      pe;
        double     worst_p = 0.0;
        const auto sp = [](int k) { return (k < 300) ? 185.0 : 165.0; };
        for (int k = 0; k < 400; ++k) {
            decltype(law_reg)::Preview ahead_reg{};
            decltype(law_pre)::Preview ahead_pre{};
            for (std::size_t j = 0; j < 32; ++j) {
                // the same information as the 30-entry preview below: the predictor form holds its
                // last value over the D steps its plan reaches beyond the preview
                ahead_reg[j] = ColVec<1>{sp(k + 1 + static_cast<int>(std::min<std::size_t>(j, 29)))};
            }
            for (std::size_t j = 0; j < 30; ++j) {
                ahead_pre[j] = ColVec<1>{sp(k + 1 + static_cast<int>(j))}; // shifted by D inside: the plan starts at k + 2
            }
            const double ua = law_reg.control(ahead_reg, ColVec<1>{pc.y()})[0];
            const double ub = law_pre.control(ahead_pre, ColVec<1>{pe.y()})[0];
            worst_p = std::max(worst_p, std::abs(ua - ub));
            pc.step(ua, 0.0);
            pe.step(ub, 0.0);
        }
        CHECK(worst_p <= 1e-6 * 40.0);
        CHECK(std::abs(pc.y() - 165.0) < 1e-3);
    }

    TEST_CASE("a 50-sample delay stays at the plant's size and is offset-free") {
        // 10 kHz loop, 5 ms transport delay: the register form would be 51 states.
        constexpr double Ts = 1e-4;
        constexpr auto   split = split_input_delay(5e-3, Ts);
        static_assert(split.whole == 50);
        static_assert(split.fraction == 0.0);
        const StateSpace<1, 1, 1> block{.A = {{-1.0 / 0.02}}, .B = {{6.0 / 0.02}}, .C = {{1.0}}}; // τ = 20 ms
        const auto                model = *discretize(block, Ts);
        const auto                res = design::discrete_oflqg<50.0>(model, Matrix<1, 1>{{1.0}}, Matrix<1, 1>{{0.1}}, design::OffsetFreeNoise{Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{0.05}}});
        REQUIRE(res.success);
        OFLQG law{res};
        static_assert(sizeof(law.delay.buffer) == 51 * sizeof(ColVec<1, double>));

        // the true plant with the 50-sample delay: a ring of past commands
        double                  x = 0.0;
        damp::array<double, 50> pipe{}; // u_{k-50} reaches the plant at k
        size_t                  head = 0;
        double                  worst_over = 0.0;
        for (int k = 0; k < 4000; ++k) {
            const double y = x;
            const double u = law.control(ColVec<1>{10.0}, ColVec<1>{y})[0];
            const double delayed = pipe[head];
            pipe[head] = u;
            head = (head + 1) % pipe.size();
            const double w = (k >= 2000) ? 0.5 : 0.0;
            x = model.A(0, 0) * x + model.B(0, 0) * (delayed + w);
            if (k < 2000) {
                worst_over = std::max(worst_over, y - 10.0); // heat-up only; the load step after it is a separate transient
            }
        }
        CHECK(std::abs(x - 10.0) < 1e-6);
        CHECK(worst_over < 0.05); // the predictor sees the heat in flight: no overshoot from the delay
    }
}
