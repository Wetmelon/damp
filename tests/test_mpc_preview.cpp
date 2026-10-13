// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file test_mpc_preview.cpp
 * @brief Reference preview: a held preview is the single-reference tick exactly; a step known
 *        ahead is acted on before it happens; per-move input targets reach the cost.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "damp/controllers/mpc.hpp"
#include "damp/controllers/offset_free_mpc.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "doctest.h"

using namespace damp;

namespace {

/// First-order heater block, K = 6, τ = 30 s, Ts = 2 s (no delay: the preview is the point here).
constexpr auto        kModel = discretize(StateSpace<1, 1, 1>{.A = {{-1.0 / 30.0}}, .B = {{6.0 / 30.0}}, .C = {{1.0}}}, 2.0);
constexpr std::size_t NP = 20;
constexpr std::size_t NC = 4;

design::MPCArtifacts<1, 1, 1, NP, NC, 0, double> plain_design(double ru = 0.0) {
    design::MPCWeights<1, 1> w{};
    w.Qy = Matrix<1, 1>{{1.0}};
    w.Rdu = Matrix<1, 1>{{1.0}};
    w.Ru = Matrix<1, 1>{{ru}};
    w.Qy_terminal = w.Qy;
    design::MPCConstraints<1, 1> c{};
    c.u_min = ColVec<1>{0.0};
    c.u_max = ColVec<1>{40.0};
    const auto art = design::state_mpc<NP, NC>(*kModel, w, c);
    REQUIRE(art.success);
    return art;
}

/// The setpoint profile: 100 until step kStep, then 80.
constexpr int kStep = 60;
double        setpoint_at(int k) {
    return (k < kStep) ? 100.0 : 80.0;
}

} // namespace

TEST_SUITE("MPC reference preview") {
    TEST_CASE("a held preview is the single-reference tick, exactly") {
        const auto                      art = plain_design();
        MPC<1, 1, 1, NP, NC, 0, double> a{art};
        MPC<1, 1, 1, NP, NC, 0, double> b{art};
        double                          xa = 0.0;
        double                          xb = 0.0;
        for (int k = 0; k < 120; ++k) {
            const double r = setpoint_at(k);
            const double ua = a.control(ColVec<1>{r}, ColVec<1>{xa})[0];
            const double ub = b.control(MPC<1, 1, 1, NP, NC, 0, double>::held(ColVec<1>{r}), ColVec<1>{xb})[0];
            CHECK(ua == ub);
            xa = kModel->A(0, 0) * xa + kModel->B(0, 0) * ua;
            xb = kModel->A(0, 0) * xb + kModel->B(0, 0) * ub;
        }
    }

    TEST_CASE("a step known over the horizon is acted on before it happens") {
        const auto                      art = plain_design();
        MPC<1, 1, 1, NP, NC, 0, double> with{art};
        MPC<1, 1, 1, NP, NC, 0, double> without{art};
        double                          x_with = 0.0;
        double                          x_without = 0.0;
        double                          err_with = 0.0;
        double                          err_without = 0.0;
        double                          u_with_before = 0.0;    // the move two steps before the step, told
        double                          u_without_before = 0.0; // and not told
        for (int k = 0; k < 120; ++k) {
            MPC<1, 1, 1, NP, NC, 0, double>::Preview ahead{};
            for (std::size_t j = 0; j < NP; ++j) {
                ahead[j] = ColVec<1>{setpoint_at(k + 1 + static_cast<int>(j))}; // the next NP setpoints
            }
            const double u_a = with.control(ahead, ColVec<1>{x_with})[0];
            const double u_b = without.control(ColVec<1>{setpoint_at(k)}, ColVec<1>{x_without})[0];
            if (k == kStep - 2) {
                u_with_before = u_a;
                u_without_before = u_b;
            }
            x_with = kModel->A(0, 0) * x_with + kModel->B(0, 0) * u_a;
            x_without = kModel->A(0, 0) * x_without + kModel->B(0, 0) * u_b;
            if (k >= kStep - static_cast<int>(NP) && k < 120) {
                err_with += std::abs(x_with - setpoint_at(k + 1));
                err_without += std::abs(x_without - setpoint_at(k + 1));
            }
        }
        CHECK(u_with_before < u_without_before - 1.0); // power already coming down before the step; the held law has not moved
        CHECK(err_with < 0.5 * err_without);           // and the error over the transition is less than half
        CHECK(std::abs(x_with - 80.0) < 1e-6);
        CHECK(std::abs(x_without - 80.0) < 1e-6);
    }

    TEST_CASE("per-move input targets: a constant preview equals the single target, and a ramp differs") {
        const auto                      art = plain_design(0.5);
        MPC<1, 1, 1, NP, NC, 0, double> one{art};
        MPC<1, 1, 1, NP, NC, 0, double> same{art};
        MPC<1, 1, 1, NP, NC, 0, double> ramp{art};
        one.set_input_target(ColVec<1>{10.0});
        MPC<1, 1, 1, NP, NC, 0, double>::InputPreview flat{};
        MPC<1, 1, 1, NP, NC, 0, double>::InputPreview rising{};
        for (std::size_t j = 0; j < NC; ++j) {
            flat[j] = ColVec<1>{10.0};
            rising[j] = ColVec<1>{10.0 + (5.0 * static_cast<double>(j))};
        }
        same.set_input_target(flat);
        ramp.set_input_target(rising);
        CHECK(same.input_target()(0) == 10.0);
        const ColVec<1> x{50.0};
        const double    u_one = one.feedback(ColVec<1>{60.0}, x)[0];
        const double    u_same = same.feedback(ColVec<1>{60.0}, x)[0];
        const double    u_ramp = ramp.feedback(ColVec<1>{60.0}, x)[0];
        CHECK(u_one == u_same);
        CHECK(u_ramp != u_one); // the later targets reach the first move through the shared Δu
        // The map itself: move i's column collects 2·Ru from every move j ≥ i.
        for (std::size_t i = 0; i < NC; ++i) {
            for (std::size_t j = 0; j < NC; ++j) {
                CHECK(art.F_ut(i, j) == ((j >= i) ? 2.0 * 0.5 : 0.0));
            }
        }
    }

    TEST_CASE("OffsetFreeMPC: a held preview is the single-reference tick through the estimator") {
        design::MPCWeights<1, 1> w{};
        w.Qy = Matrix<1, 1>{{1.0}};
        w.Rdu = Matrix<1, 1>{{1.0}};
        w.Qy_terminal = w.Qy;
        design::MPCConstraints<1, 1> c{};
        c.u_min = ColVec<1>{0.0};
        c.u_max = ColVec<1>{40.0};
        const auto art = design::mpc<NP, NC>(*kModel, w, c, design::OffsetFreeNoise{Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{0.05}}});
        REQUIRE(art.success);
        OffsetFreeMPC<1, 1, 1, NP, NC, double> a{art};
        OffsetFreeMPC<1, 1, 1, NP, NC, double> b{art};
        double                                 xa = 0.0;
        double                                 xb = 0.0;
        for (int k = 0; k < 120; ++k) {
            const double r = setpoint_at(k);
            const double w_in = (k >= 90) ? 1.0 : 0.0;
            const double ua = a.control(ColVec<1>{r}, ColVec<1>{xa})[0];
            const double ub = b.control(MPC<1, 1, 1, NP, NC, 1, double>::held(ColVec<1>{r}), ColVec<1>{xb})[0];
            CHECK(ua == ub);
            xa = kModel->A(0, 0) * xa + kModel->B(0, 0) * (ua + w_in);
            xb = kModel->A(0, 0) * xb + kModel->B(0, 0) * (ub + w_in);
        }
        CHECK(std::abs(xa - 80.0) < 1e-6);
    }
}
