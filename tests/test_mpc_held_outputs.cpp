// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file test_mpc_held_outputs.cpp
 * @brief The outputs an MPC's target holds at r (TrackedOutputs, or a disturbance spec's list) apart
 *        from the outputs its QP weighs, every one; Ru measured from the target's input; at most
 *        one held output per input, in every design with a target; and with neither, the
 *        regulator, as an LQR is.
 */

#include <cmath>
#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/controllers/state_feedback.hpp"
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

/// The lags with an input on each, both read: as many outputs as inputs. Sampled at kTs.
constexpr auto kTwoInputs = discretize(StateSpace<2, 2, 2>{.A = {{-1.0, 0.0}, {1.0, -1.0}}, .B = Matrix<2, 2>::identity(), .C = Matrix<2, 2>::identity()}, kTs);

/// The same with noise on both states and both readings, for a filter.
constexpr StateSpace<2, 2, 2, double, 2, 2> kNoisyTwo{.A = kTwoInputs->A, .B = kTwoInputs->B, .C = kTwoInputs->C, .Ts = kTs};

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

/// Both readings of the two-input lags weighed, the second the more; a move and an input weight on each input.
design::MPCWeights<2, 2> two_weighed() {
    design::MPCWeights<2, 2> w{.Qy = diag(1.0, 10.0), .Rdu = Matrix<2, 2>::identity() * 0.1, .Ru = Matrix<2, 2>::identity() * 0.2};
    w.Qy_terminal = w.Qy;
    return w;
}

/// The two-input lags a step on, the load @p w on x₂.
ColVec<2> two_step(const ColVec<2>& x, const ColVec<2>& u, double w = 0.0) {
    return (kTwoInputs->A * x) + (kTwoInputs->B * u) + (kLoadOnX2 * ColVec<1>{w});
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
        const auto art = design::mpc<20, 5>(*lag, w, design::TrackedOutputs{});
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
            const auto same = design::steady_state_target(*kLags, design::TrackedOutputs{{1}}, design::MeasuredDisturbance{kLoadOnX2});
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
        //! A law that holds x₂ with r previewed equals one that holds x₁, handed the target's x₁
        //! for that preview: each step's reference on the output r does not hold and each move's
        //! input target are the steady state of that step's r.
        const auto tracked = design::mpc<20, 5>(*kLags, both_weighed(0.2), design::TrackedOutputs{{1}});
        const auto other = design::mpc<20, 5>(*kLags, both_weighed(0.2), design::TrackedOutputs{{0}});
        REQUIRE(tracked.success);
        REQUIRE(other.success);
        MPC       held{tracked};
        MPC       given{other};
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
                y_r[j] = ColVec<1>{tracked.target.target_r(0, 0) * rk}; // C = I: the target's x₁
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

    TEST_CASE("TrackedOutputs{} is every output, the list of all of them in order; out of order it is refused") {
        //! With every output held the runtime takes r itself as the QP's reference on each, so a
        //! full list is read in the plant's order: {{1, 0}} would hand each r to the other output.
        const auto every = design::mpc<10, 4>(*kTwoInputs, {}, design::TrackedOutputs{});
        const auto listed = design::mpc<10, 4>(*kTwoInputs, {}, design::TrackedOutputs{{0, 1}});
        static_assert(std::is_same_v<decltype(every), decltype(listed)>);
        REQUIRE(every.success);
        REQUIRE(listed.success);
        CHECK(same_map(every.target.target_r, listed.target.target_r));
        CHECK(design::mpc<10, 4>(*kTwoInputs, {}, design::TrackedOutputs{{1, 0}}).status == design::DesignStatus::bad_weights);
        CHECK(design::state_feedback(*kTwoInputs, {.Q = Matrix<2, 2>::identity(), .R = Matrix<2, 2>::identity()}, design::TrackedOutputs{{1, 0}}).success); // an LQ law has no QP to read r
    }

    TEST_CASE("a target holds at most one output per input: the count each design's options hold") {
        //! What check_held_outputs counts, and refuses past NU, in state_feedback, mpc and
        //! steady_state_target alike: a spec's list, else a TrackedOutputs', every output when
        //! the list is empty; with neither the design is a regulator, with no target, or under a
        //! measured d every output held at 0, or an LQI's integrated ones. Past NU it does not
        //! compile, a static_assert naming TrackedOutputs{{i, ...}} and the spec's .outputs; there
        //! is no compile-fail harness, so each refusal (a regulator under d, a TrackedOutputs{} or
        //! {{0, 1}}, an InputDisturbance{}, an LQI of both outputs under d) was checked by hand, on
        //! g++ and clang.
        using design::detail::held_count;
        using design::detail::reference_outputs;
        using M = design::MeasuredDisturbance<2, 1>;
        static_assert(held_count<1, 2, 0>() == 0);    // a regulator with no target: design::mpc(*kLags, w) compiles
        static_assert(held_count<1, 2, 0, M>() == 2); // a regulator under d: every output at 0, refused
        static_assert(held_count<1, 2, 1, M>() == 1); // an LQI under d: its one integrated output at 0
        static_assert(held_count<1, 2, 0, design::TrackedOutputs<>>() == 2);
        static_assert(held_count<1, 2, 0, design::TrackedOutputs<1>, M>() == 1);
        static_assert(held_count<1, 2, 0, design::InputDisturbance<>>() == 2);
        static_assert(held_count<1, 2, 0, design::InputDisturbance<1>, M>() == 1);
        static_assert(held_count<1, 2, 0, design::DisturbanceModel<2, 0, 1, 1>>() == 1);
        static_assert(held_count<2, 2, 0, design::TrackedOutputs<>>() == 2); // as many as the inputs: kTwoInputs
        //! The outputs r is on: none for a regulator.
        static_assert(reference_outputs<design::NoDisturbance, 2> == 0);
        static_assert(reference_outputs<design::TrackedOutputs<>, 2> == 2 && reference_outputs<design::TrackedOutputs<1>, 2> == 1);
        static_assert(reference_outputs<design::InputDisturbance<>, 2> == 2);

        //! The designs that compile on the lags, one input and two readings: one output listed.
        const design::KalmanNoise load{.process = Matrix<2, 2>::identity() * 1e-4, .measurement = Matrix<2, 2>::identity() * 1e-6, .disturbance = 1e-2};
        const M                   d{kLoadOnX2};
        static_assert(std::is_same_v<decltype(design::mpc<10, 4>(*kLags, both_weighed(), design::TrackedOutputs{{1}}).target), design::SteadyStateTargetResult<2, 1, 0, 1, double>>);
        static_assert(std::is_same_v<decltype(design::mpc<10, 4>(kNoisyLags, both_weighed(), load, design::InputDisturbance{{0}}, d).target), design::SteadyStateTargetResult<2, 1, 2, 1, double>>);
        static_assert(std::is_same_v<decltype(design::state_feedback(*kLags, {.Q = Matrix<2, 2>::identity(), .R = 1.0}, design::TrackedOutputs{{0}}, d).target), design::SteadyStateTargetResult<2, 1, 1, 1, double>>);
        static_assert(std::is_same_v<decltype(design::state_feedback(*kLags, {.Q = Matrix<2, 2>::identity(), .R = 1.0}, design::IntegralAction{.Q = 1.0, .outputs = {1}}, d).target), design::SteadyStateTargetResult<2, 1, 1, 0, double>>);
        static_assert(std::is_same_v<decltype(design::steady_state_target(*kLags, design::TrackedOutputs{{1}}, d)), design::SteadyStateTargetResult<2, 1, 1, 1, double>>);
        //! A regulator with no measured d has no target: every output read, nothing held.
        static_assert(std::is_same_v<decltype(design::state_feedback(*kLags, {.Q = Matrix<2, 2>::identity(), .R = 1.0})), design::StateFeedbackResult<2, 1, double>>);
        static_assert(std::is_same_v<decltype(design::mpc<10, 4>(*kLags, both_weighed()).target), design::SteadyStateTargetResult<2, 1, 0, 0, double>>);
    }

    TEST_CASE("as many held outputs as inputs: every output held at r by TrackedOutputs{}, in all three designs") {
        //! kTwoInputs: two inputs, two readings, so the target holds both, square and unique.
        const auto tg = design::steady_state_target(*kTwoInputs, design::TrackedOutputs{});
        const auto mpc = design::mpc<20, 5>(*kTwoInputs, {}, design::TrackedOutputs{});
        const auto lqr = design::state_feedback(*kTwoInputs, {.Q = Matrix<2, 2>::identity(), .R = Matrix<2, 2>::identity()}, design::TrackedOutputs{});
        REQUIRE(tg.success);
        REQUIRE(mpc.success);
        REQUIRE(lqr.success);
        static_assert(std::is_same_v<decltype(tg), const design::SteadyStateTargetResult<2, 2, 0, 2, double>>);
        CHECK(same_map(tg.target_r, mpc.target.target_r));
        CHECK(same_map(tg.target_r, lqr.target.target_r));
        const ColVec<2> r{0.7, -0.2};
        const ColVec<4> s = tg.target_r * r;
        CHECK(s(0) == doctest::Approx(0.7).epsilon(1e-12)); // both readings at their r
        CHECK(s(1) == doctest::Approx(-0.2).epsilon(1e-12));

        MPC       law{mpc};
        ColVec<2> x{};
        for (size_t k = 0; k < 400; ++k) {
            x = (kTwoInputs->A * x) + (kTwoInputs->B * law.step(r, x));
        }
        CHECK(x(0) == doctest::Approx(0.7).epsilon(1e-9));
        CHECK(x(1) == doctest::Approx(-0.2).epsilon(1e-9));
    }

    TEST_CASE("at r = 0 the list does not move an MPC: the commands it gave holding both outputs") {
        //! Regulated to 0 the target is 0 whichever outputs it holds, so listing one gives the
        //! commands of the design that held both, recorded on the base before it was refused.
        const auto first = design::mpc<20, 5>(*kLags, both_weighed(), design::TrackedOutputs{{0}});
        const auto second = design::mpc<20, 5>(*kLags, both_weighed(), design::TrackedOutputs{{1}});
        REQUIRE(first.success);
        REQUIRE(second.success);
        MPC       a{first};
        MPC       b{second};
        ColVec<2> x{1.0, -0.5};
        for (size_t k = 0; k < 30; ++k) {
            const double u = a.step(0.0, x)(0);
            CHECK(b.step(0.0, x)(0) == doctest::Approx(u).epsilon(1e-12)); // the same QP: r = 0 holds every output at 0
            if (k == 0) {
                CHECK(u == doctest::Approx(-0.52604853604519486).epsilon(1e-12));
            } else if (k == 5) {
                CHECK(u == doctest::Approx(-0.4879430919665545).epsilon(1e-12));
            } else if (k == 20) {
                CHECK(u == doctest::Approx(-0.0084323571676990922).epsilon(1e-10));
            }
            x = lags_step(x, u);
        }
        CHECK(x(0) == doctest::Approx(-0.0040284944597127779).epsilon(1e-9));
        CHECK(x(1) == doctest::Approx(0.00074001389230396908).epsilon(1e-9));
    }

    TEST_CASE("a held output past the plant's or listed twice is refused, as a disturbance spec's is") {
        CHECK(design::mpc<10, 4>(*kLags, both_weighed(), design::TrackedOutputs{{2}}).status == design::DesignStatus::bad_weights);
        CHECK(design::mpc<10, 4>(*kTwoInputs, {}, design::TrackedOutputs{{1, 1}}).status == design::DesignStatus::bad_weights); // two listed: a plant with two inputs
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
        CHECK(design::mpc<10, 4>(through, first_only, design::TrackedOutputs{{0}}, design::OutputLimits{.outputs = {1}, .upper = {1.0}}).status == design::DesignStatus::bad_weights);
    }
}

TEST_SUITE("MPC regulator") {
    TEST_CASE("an MPC with no TrackedOutputs is a regulator: no r, and the commands of TrackedOutputs{} at r = 0") {
        //! As an LQR: nothing is held at r, so the runtime takes none, and the QP is posed about
        //! x = 0, u = 0, the target TrackedOutputs{} gives at r = 0.
        const auto reg = design::mpc<20, 5>(*kTwoInputs, two_weighed());
        const auto every = design::mpc<20, 5>(*kTwoInputs, two_weighed(), design::TrackedOutputs{});
        REQUIRE(reg.success);
        REQUIRE(every.success);
        static_assert(std::is_same_v<decltype(reg.target), design::SteadyStateTargetResult<2, 2, 0, 0, double>>); // no map at all
        MPC       a{reg};
        MPC       b{every};
        MPC       fed_a{reg};
        MPC       fed_b{every};
        ColVec<2> x{1.0, -0.5};
        ColVec<2> xf{1.0, -0.5};
        for (size_t k = 0; k < 40; ++k) {
            const ColVec<2> u = a.step(x);
            const ColVec<2> v = b.step(ColVec<2>{}, x);
            CHECK(u(0) == doctest::Approx(v(0)).epsilon(1e-12));
            CHECK(u(1) == doctest::Approx(v(1)).epsilon(1e-12));
            x = two_step(x, u);
            //! About a planned point, the feedforward form without r.
            const ColVec<2> x_ref{0.2, 0.1};
            const ColVec<2> u_ref{0.2, -0.1};
            const ColVec<2> uf = fed_a.step(x_ref, u_ref, xf);
            const ColVec<2> vf = fed_b.step(ColVec<2>{}, x_ref, u_ref, xf);
            CHECK(uf(0) == doctest::Approx(vf(0)).epsilon(1e-12));
            CHECK(uf(1) == doctest::Approx(vf(1)).epsilon(1e-12));
            xf = two_step(xf, uf);
        }
        CHECK(Matrix<2, 1>(x).norm() < 1e-3); // regulated to 0
        a.update(x);
        b.update(x);
        CHECK(a.control()(1) == doctest::Approx(b.control(ColVec<2>{})(1)).epsilon(1e-12));

        //! It stores no map of r, and no target: the (x + u) × y entries of TrackedOutputs{}'s go.
        if constexpr (sizeof(size_t) == 8) {
            CHECK(sizeof(a) == 9592);
            CHECK(sizeof(b) == 9664);
        }
    }

    TEST_CASE("a regulator under a measured disturbance: the steady state that holds every output at 0, TrackedOutputs{}'s at r = 0") {
        const design::MeasuredDisturbance load{kLoadOnX2};
        const auto                        reg = design::mpc<20, 5>(*kTwoInputs, two_weighed(), load);
        const auto                        every = design::mpc<20, 5>(*kTwoInputs, two_weighed(), design::TrackedOutputs{}, load);
        REQUIRE(reg.success);
        REQUIRE(every.success);
        static_assert(std::is_same_v<decltype(reg.target), design::SteadyStateTargetResult<2, 2, 1, 0, double>>); // d's map, no r's
        //! The same target as steady_state_target's, the LQR regulator's, and TrackedOutputs{}'s at r = 0.
        const auto tg = design::steady_state_target(*kTwoInputs, load);
        const auto lqr = design::state_feedback(*kTwoInputs, {.Q = Matrix<2, 2>::identity(), .R = Matrix<2, 2>::identity()}, load);
        REQUIRE(tg.success);
        REQUIRE(lqr.success);
        CHECK(same_map(reg.target.target_d, tg.target_d));
        CHECK(same_map(reg.target.target_d, lqr.target.target_d));
        CHECK(same_map(reg.target.target_d, every.target.target_d));

        MPC          a{reg};
        MPC          b{every};
        ColVec<2>    x{0.4, -0.3};
        const double w = 0.5;
        for (size_t k = 0; k < 400; ++k) {
            const ColVec<2> u = a.step(ColVec<1>{w}, x);
            const ColVec<2> v = b.step(ColVec<2>{}, ColVec<1>{w}, x);
            CHECK(u(0) == doctest::Approx(v(0)).epsilon(1e-12));
            CHECK(u(1) == doctest::Approx(v(1)).epsilon(1e-12));
            x = two_step(x, u, w);
        }
        CHECK(Matrix<2, 1>(x).norm() < 1e-9); // both outputs at 0 under the load
        a.update(x);
        b.update(x);
        CHECK(a.control(ColVec<1>{w})(0) == doctest::Approx(b.control(ColVec<2>{}, ColVec<1>{w})(0)).epsilon(1e-12));
        CHECK(a.control(ColVec<2>{}, ColVec<2>{}, ColVec<1>{w})(1) == doctest::Approx(b.control(ColVec<2>{}, ColVec<2>{}, ColVec<2>{}, ColVec<1>{w})(1)).epsilon(1e-12));
        if constexpr (sizeof(size_t) == 8) {
            CHECK(sizeof(a) == 9880);
            CHECK(sizeof(b) == 9944);
        }
    }

    TEST_CASE("on a plant with more outputs than inputs a regulator needs no list") {
        //! The lags read both states on one input: nothing held, so nothing refused. Under a
        //! measured d it would hold both at 0 and does not compile: TrackedOutputs{{i}} with r = 0.
        const auto reg = design::mpc<20, 5>(*kLags, both_weighed());
        REQUIRE(reg.success);
        MPC       law{reg};
        ColVec<2> x{1.0, -0.5};
        for (size_t k = 0; k < 300; ++k) {
            x = lags_step(x, law.step(x)(0));
        }
        CHECK(Matrix<2, 1>(x).norm() < 1e-6);
    }

    TEST_CASE("through a filter: the regulator's design is TrackedOutputs{}'s QP and filter, with no map of r") {
        //! design::mpc with a KalmanNoise and no spec: the QP and the filter do not depend on the
        //! target, so they are those of the design that holds every output. Neither has a runtime:
        //! OFMPC runs with a disturbance spec, whose outputs it holds at r.
        const design::KalmanNoise noise{.process = Matrix<2, 2>::identity() * 1e-4, .measurement = Matrix<2, 2>::identity() * 1e-6};
        const auto                reg = design::mpc<10, 4>(kNoisyTwo, two_weighed(), noise);
        const auto                every = design::mpc<10, 4>(kNoisyTwo, two_weighed(), noise, design::TrackedOutputs{});
        REQUIRE(reg.success);
        REQUIRE(every.success);
        static_assert(std::is_same_v<decltype(reg.target), design::SteadyStateTargetResult<2, 2, 0, 0, double>>);
        CHECK(same_map(reg.mpc.H, every.mpc.H));
        CHECK(same_map(reg.kalman.M, every.kalman.M));
        const design::MeasuredDisturbance load{kLoadOnX2};
        const auto                        reg_d = design::mpc<10, 4>(kNoisyTwo, two_weighed(), noise, load);
        const auto                        every_d = design::mpc<10, 4>(kNoisyTwo, two_weighed(), noise, design::TrackedOutputs{}, load);
        REQUIRE(reg_d.success);
        REQUIRE(every_d.success);
        static_assert(std::is_same_v<decltype(reg_d.target), design::SteadyStateTargetResult<2, 2, 1, 0, double>>);
        CHECK(same_map(reg_d.target.target_d, every_d.target.target_d));
        CHECK(same_map(reg_d.target.target_d, design::steady_state_target(kNoisyTwo, load).target_d));
    }
}
