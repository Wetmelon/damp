// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file test_mpc_held_outputs.cpp
 * @brief The outputs an MPC's target holds at r (TrackedOutputs, or a disturbance spec's list) apart
 *        from the outputs its QP weighs, every one; and Ru measured from the target's input.
 */

#include <cmath>
#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/controllers/state_feedback.hpp"
#include "damp/design/config.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/qp.hpp"
#include "damp/design/status.hpp"
#include "damp/design/steady_state_target.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "doctest.h"

using namespace damp;

namespace {

constexpr double kTs = 0.1;

/// Two lags in a row, both read: x₁' = −x₁ + u, x₂' = x₁ − x₂; at rest x₂ = x₁ = u. Sampled at kTs.
constexpr auto kLags = discretize(StateSpace<2, 1, 2>{.A = {{-1.0, 0.0}, {1.0, -1.0}}, .B = {{1.0}, {0.0}}, .C = Matrix<2, 2>::identity()}, kTs);

/// The same with noise on both states and both readings, for a filter.
constexpr StateSpace<2, 1, 2, double, 2, 2> kNoisyLags{.A = kLags->A, .B = kLags->B, .C = kLags->C, .Ts = kTs};

/// A load on the second lag, per step: x₂ takes kTs·w more each period.
constexpr Matrix<2, 1> kLoadOnX2{{0.0}, {kTs}};

/// Both outputs weighed, the second the more; Ru damps the input.
design::MPCWeights<1, 2> both_weighed(double ru = 0.5) {
    design::MPCWeights<1, 2> w{.Qy = diag(1.0, 10.0), .Rdu = 0.1, .Ru = ru};
    w.Qy_terminal = w.Qy;
    return w;
}

/// The plant a step on: the lags under @p u, the load @p w on x₂.
ColVec<2> lags_step(const ColVec<2>& x, double u, double w = 0.0) {
    return (kLags->A * x) + (kLags->B * ColVec<1>{u}) + (kLoadOnX2 * ColVec<1>{w});
}

/// Whether two maps agree to roundoff: different instantiations of one formula may round apart.
template<size_t R, size_t C>
bool same_map(const Matrix<R, C>& a, const Matrix<R, C>& b) {
    return Matrix<R, C>(a - b).norm() <= 1e-12 * (1.0 + a.norm());
}

/// x_t, u_t of a design's target for @p r under @p d: [x_t; u_t] = N_d·d + N_r·r.
template<typename Target>
ColVec<3> target_point(const Target& t, double r, double d = 0.0) {
    ColVec<3> s{};
    for (size_t i = 0; i < 3; ++i) {
        s(i) = t.target_r(i, 0) * r;
        if constexpr (std::remove_cvref_t<decltype(t.target_d)>::cols() > 0) {
            s(i) += t.target_d(i, 0) * d;
        }
    }
    return s;
}

} // namespace

TEST_SUITE("MPC held outputs") {
    TEST_CASE("an offset-free MPC holds the output its spec lists and weighs both: step(r, y) alone, the other output where the target puts it") {
        //! The load on x₂ moves the steady state: x₂ = r needs x₁ = r − w, so the first output
        //! cannot sit at r. DisturbanceModel{.outputs = {1}} holds the second at r; the QP weighs
        //! both, the first toward the target's x₁, and Ru toward the target's input, by default.
        const design::KalmanNoise noise{.process = Matrix<2, 2>::identity() * 1e-4, .measurement = Matrix<2, 2>::identity() * 1e-6, .disturbance = 1e-2};
        const auto                art = design::mpc<20, 5>(kNoisyLags, both_weighed(), noise, design::DisturbanceModel{.Bd = kLoadOnX2, .outputs = {1}});
        REQUIRE(art.success);
        static_assert(decltype(OFMPC{art})::NDE == 1);
        static_assert(std::is_same_v<decltype(OFMPC{art})::Preview, damp::array<ColVec<1>, 20>>); // r on the one held output

        OFMPC        law{art};
        ColVec<2>    x{};
        const double r = 1.0;
        const double w = 0.4;
        for (size_t k = 0; k < 600; ++k) {
            x = lags_step(x, law.step(r, kLags->C * x)(0), w);
        }
        const ColVec<3> t = target_point(art.target, r, w);
        CHECK(x(1) == doctest::Approx(r).epsilon(1e-6));                          // held at r, offset-free
        CHECK(x(0) == doctest::Approx(t(0)).epsilon(1e-6));                       // the other where the target puts it
        CHECK(std::abs(x(0) - r) > 0.3);                                          // which is not r
        CHECK(law.disturbance_estimate()(0) == doctest::Approx(w).epsilon(1e-6)); // the load found
        CHECK(law.applied()(0) == doctest::Approx(t(2)).epsilon(1e-6));           // and the input the target's
        CHECK(law.input_target()(0) == doctest::Approx(t(2)).epsilon(1e-6));      // Ru's anchor, followed
        CHECK(law.last_status() == design::QPStatus::Success);
    }

    TEST_CASE("Ru measures the input from the target's: no offset by default; set_input_target's leaves one, clear_input_target returns") {
        //! One lag, y = x, gain 2: r = 1 takes u = 0.5. Ru on u − 0 would pull the plant short of r.
        constexpr auto           lag = discretize(StateSpace<1, 1, 1>{.A = {{-1.0}}, .B = {{2.0}}, .C = {{1.0}}}, kTs);
        design::MPCWeights<1, 1> w{.Qy = 1.0, .Rdu = 0.1, .Ru = 1.0};
        w.Qy_terminal = w.Qy;
        const auto art = design::mpc<20, 5>(*lag, w);
        REQUIRE(art.success);
        const double u_t = art.target.target_r(1, 0);
        CHECK(u_t == doctest::Approx(0.5).epsilon(1e-12));

        MPC        law{art};
        MPC        by_hand{art};
        double     x = 0.0;
        double     x_hand = 0.0;
        const auto settle = [&](auto& m, double& state) {
            for (size_t k = 0; k < 400; ++k) {
                state = (lag->A(0, 0) * state) + (lag->B(0, 0) * m.step(1.0, ColVec<1>{state})(0));
            }
        };
        by_hand.set_input_target(u_t); // the target's input written out: the same law
        for (size_t k = 0; k < 400; ++k) {
            const double u = law.step(1.0, ColVec<1>{x})(0);
            CHECK(u == doctest::Approx(by_hand.step(1.0, ColVec<1>{x_hand})(0)).epsilon(1e-12));
            x = (lag->A(0, 0) * x) + (lag->B(0, 0) * u);
            x_hand = x;
        }
        CHECK(x == doctest::Approx(1.0).epsilon(1e-9)); // no offset
        CHECK(law.input_target()(0) == doctest::Approx(u_t).epsilon(1e-12));

        law.set_input_target(0.0); // an operating point the application prefers: Ru trades it against r
        settle(law, x);
        CHECK(x < 0.95);
        CHECK(law.input_target()(0) == 0.0);

        law.clear_input_target(); // back to the target's input
        settle(law, x);
        CHECK(x == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(law.input_target()(0) == doctest::Approx(u_t).epsilon(1e-12));
    }

    TEST_CASE("the measured-state MPC with TrackedOutputs: r on the listed output, every output weighed") {
        const auto art = design::mpc<20, 5>(*kLags, both_weighed(), design::TrackedOutputs{{1}});
        REQUIRE(art.success);
        static_assert(std::is_same_v<decltype(MPC{art})::Preview, damp::array<ColVec<1>, 20>>);
        //! The target is the LQR's with the same TrackedOutputs: one path, steady_state_target's.
        const auto lqr = design::state_feedback(*kLags, {.Q = Matrix<2, 2>::identity(), .R = 1.0}, design::TrackedOutputs{{1}});
        REQUIRE(lqr.success);
        CHECK(same_map(art.target.target_r, lqr.target.target_r));
        CHECK(art.target.target_r(1, 0) == doctest::Approx(1.0).epsilon(1e-12));

        MPC       law{art};
        ColVec<2> x{};
        for (size_t k = 0; k < 400; ++k) {
            x = lags_step(x, law.step(0.7, x)(0));
        }
        CHECK(x(1) == doctest::Approx(0.7).epsilon(1e-9));
        CHECK(x(0) == doctest::Approx(0.7).epsilon(1e-9)); // the target's x₁: at rest both lags agree

        SUBCASE("under a measured disturbance, the target under it") {
            const auto md = design::mpc<20, 5>(*kLags, both_weighed(), design::MeasuredDisturbance{kLoadOnX2}, design::TrackedOutputs{{1}});
            REQUIRE(md.success);
            const auto same = design::steady_state_target(*kLags, design::DisturbanceModel{.Bd = kLoadOnX2, .outputs = {1}});
            CHECK(same_map(md.target.target_d, same.target_d));
            CHECK(same_map(md.target.target_r, same.target_r));
            MPC          ml{md};
            ColVec<2>    z{};
            const double w = 0.4;
            for (size_t k = 0; k < 400; ++k) {
                z = lags_step(z, ml.step(0.7, ColVec<1>{w}, z)(0), w);
            }
            const ColVec<3> t = target_point(md.target, 0.7, w);
            CHECK(z(1) == doctest::Approx(0.7).epsilon(1e-9));
            CHECK(z(0) == doctest::Approx(t(0)).epsilon(1e-9));
            CHECK(ml.applied()(0) == doctest::Approx(t(2)).epsilon(1e-9));
        }
    }

    TEST_CASE("the reference preview drives each step's target: the held output's preview, the other at the target's") {
        //! A law that holds x₂ with r previewed equals one that holds both outputs, handed the
        //! target's outputs for that preview: each step's reference on x₁ and each move's input
        //! target are the steady state of that step's r.
        const auto tracked = design::mpc<20, 5>(*kLags, both_weighed(0.2), design::TrackedOutputs{{1}});
        const auto every = design::mpc<20, 5>(*kLags, both_weighed(0.2));
        REQUIRE(tracked.success);
        REQUIRE(every.success);
        MPC       held{tracked};
        MPC       given{every};
        MPC       blind{tracked};
        ColVec<2> xa{};
        ColVec<2> xb{};
        ColVec<2> xc{};
        double    gap = 0.0;
        double    lead = 0.0;
        for (int k = 0; k < 120; ++k) {
            decltype(held)::Preview  r{};
            decltype(given)::Preview y_r{};
            for (size_t j = 0; j < 20; ++j) {
                const double rk = (k + static_cast<int>(j) < 40) ? 0.0 : 0.02 * static_cast<double>(k + static_cast<int>(j) - 40); // a ramp from step 40
                r[j] = ColVec<1>{rk};
                for (size_t i = 0; i < 2; ++i) {
                    y_r[j](i) = tracked.target.target_r(i, 0) * rk; // C = I: the target's outputs
                }
            }
            const double ua = held.step(r, xa)(0);
            const double ub = given.step(y_r, xb)(0);
            const double uc = blind.step(r[0], xc)(0); // the same law, told r now only
            gap = damp::max(gap, std::abs(ua - ub));
            lead = (k == 35) ? ua - uc : lead;
            xa = lags_step(xa, ua);
            xb = lags_step(xb, ub);
            xc = lags_step(xc, uc);
        }
        CHECK(gap < 1e-9);
        CHECK(lead > 0.01); // the preview moves the input before the ramp starts
    }

    TEST_CASE("TrackedOutputs on a filtered MPC without a disturbance: the target holds the listed output") {
        const design::KalmanNoise noise{.process = Matrix<2, 2>::identity() * 1e-4, .measurement = Matrix<2, 2>::identity() * 1e-6};
        const auto                art = design::mpc<20, 5>(kNoisyLags, both_weighed(), noise, design::TrackedOutputs{{1}});
        REQUIRE(art.success);
        static_assert(std::is_same_v<decltype(art.target.target_r), Matrix<3, 1>>); // [x_t; u_t] per unit r on the one held output
        const auto lqg = design::state_feedback(kNoisyLags, {.Q = Matrix<2, 2>::identity(), .R = 1.0}, noise, design::TrackedOutputs{{1}});
        REQUIRE(lqg.success);
        CHECK(same_map(art.target.target_r, lqg.target.target_r)); // the LQG's with the same TrackedOutputs
    }

    TEST_CASE("MPCConfig::tracked passes TrackedOutputs on; TrackedOutputs{} is every output, as left out") {
        const auto pack = design::mpc<10, 4>(*kLags, both_weighed(), design::TrackedOutputs{{1}});
        const auto named = design::mpc<10, 4>(*kLags, design::MPCConfig{.weights = both_weighed(), .tracked = design::TrackedOutputs{{1}}});
        static_assert(std::is_same_v<decltype(pack), decltype(named)>);
        CHECK(pack.target.target_r == named.target.target_r);
        CHECK(pack.mpc.H == named.mpc.H);
        static_assert(std::is_same_v<decltype(design::mpc<10, 4>(*kLags, both_weighed(), design::TrackedOutputs{})), decltype(design::mpc<10, 4>(*kLags, both_weighed()))>);
        CHECK(same_map(design::mpc<10, 4>(*kLags, both_weighed(), design::TrackedOutputs{}).target.target_r, design::mpc<10, 4>(*kLags, both_weighed()).target.target_r));
    }

    TEST_CASE("a held output past the plant's or listed twice is refused, as a disturbance spec's is") {
        CHECK(design::mpc<10, 4>(*kLags, both_weighed(), design::TrackedOutputs{{2}}).status == design::DesignStatus::bad_weights);
        CHECK(design::mpc<10, 4>(*kLags, both_weighed(), design::TrackedOutputs{{1, 1}}).status == design::DesignStatus::bad_weights);
        const design::KalmanNoise noise{.process = Matrix<2, 2>::identity() * 1e-4, .measurement = Matrix<2, 2>::identity() * 1e-6};
        CHECK(design::mpc<10, 4>(kNoisyLags, both_weighed(), noise, design::TrackedOutputs{{3}}).status == design::DesignStatus::bad_weights);
        const design::KalmanNoise load{.process = Matrix<2, 2>::identity() * 1e-4, .measurement = Matrix<2, 2>::identity() * 1e-6, .disturbance = 1e-2};
        CHECK(design::mpc<10, 4>(kNoisyLags, both_weighed(), load, design::DisturbanceModel{.Bd = kLoadOnX2, .outputs = {2}}).status == design::DesignStatus::bad_weights);
        //! A refused design still deploys and holds its input.
        MPC refused{design::mpc<10, 4>(*kLags, both_weighed(), design::TrackedOutputs{{2}})};
        CHECK_FALSE(refused.valid());
        CHECK(refused.step(1.0, ColVec<2>{})(0) == 0.0);
    }

    TEST_CASE("a reading with feedthrough is measured, not weighed: refused only where the QP weighs or bounds it") {
        //! The QP predicts y = C·x. Before every output was the QP's, a spec's list left such a reading out.
        auto through = *kLags;
        through.D(1, 0) = 0.5;
        design::MPCWeights<1, 2> first_only{.Qy = diag(1.0, 0.0), .Rdu = 0.1};
        first_only.Qy_terminal = first_only.Qy;
        CHECK(design::mpc<10, 4>(through, first_only, design::TrackedOutputs{{0}}).success);
        CHECK(design::mpc<10, 4>(through, both_weighed(), design::TrackedOutputs{{0}}).status == design::DesignStatus::bad_weights);
        CHECK(design::mpc<10, 4>(through, first_only, design::OutputLimits{.outputs = {1}, .lim = {.upper = 1.0}}).status == design::DesignStatus::bad_weights);
    }
}
