// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file test_measured_disturbance.cpp
 * @brief MeasuredDisturbance in both families: the LQ laws' target under the measured d, the
 *        filter that takes it as a known input after u, and a measured and an estimated
 *        disturbance together, in the OFLQG and the OFMPC; and design::steady_state_target
 *        under it, the target those laws hold.
 */

#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/law.hpp"
#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/controllers/oflqg.hpp"
#include "damp/controllers/state_feedback.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/status.hpp"
#include "damp/design/steady_state_target.hpp"
#include "damp/estimation/fixed_gain_observer.hpp"
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
constexpr auto kLags = discretize(kLagsC, kTs);

/// The same with noise on both states and both readings, for a filter.
constexpr StateSpace<2, 1, 2, double, 2, 2> kNoisyLags{.A = kLags->A, .B = kLags->B, .C = kLags->C, .Ts = kTs};

/// The lags read at x₂ alone, sampled at kTs: one output per input, so a regulator holds every output.
constexpr StateSpace<2, 1, 1> kSecond{.A = kLags->A, .B = kLags->B, .C = {{0.0, 1.0}}, .Ts = kTs};

/// The same, continuous.
constexpr StateSpace<2, 1, 1> kSecondC{.A = kLagsC.A, .B = kLagsC.B, .C = {{0.0, 1.0}}};

/// A load on the second lag, per step: x₂ takes kTs·w more each period.
constexpr Matrix<2, 1> kLoadOnX2{{0.0}, {kTs}};

/// The plant a step on: the lags under @p u, the load @p w on x₂.
ColVec<2> lags_step(const ColVec<2>& x, double u, double w) {
    return (kLags->A * x) + (kLags->B * ColVec<1>{u}) + (kLoadOnX2 * ColVec<1>{w});
}

constexpr design::LQRWeights<2, 1> kWeights{.Q = Matrix<2, 2>::identity(), .R = {1.0}};

constexpr design::KalmanNoise kNoise{.process = Matrix<2, 2>::identity() * 1e-4, .measurement = Matrix<2, 2>::identity() * 1e-6};
constexpr design::KalmanNoise kLoadNoise{.process = Matrix<2, 2>::identity() * 1e-4, .measurement = Matrix<2, 2>::identity() * 1e-6, .disturbance = 1e-2};

/// Both outputs weighed, the second the more; Ru damps the input.
constexpr design::MPCWeights<1, 2> both_weighed() {
    design::MPCWeights<1, 2> w{.Qy = diag(1.0, 10.0), .Rdu = 0.1, .Ru = 0.5};
    w.Qy_terminal = w.Qy;
    return w;
}

/// Whether two maps agree to roundoff: different instantiations of one formula may round apart.
template<size_t R, size_t C>
bool same_map(const Matrix<R, C>& a, const Matrix<R, C>& b, double tol = 1e-12) {
    return Matrix<R, C>(a - b).norm() <= tol * (1.0 + a.norm());
}

/// A load on x₂, continuous: x₂' gains w.
constexpr Matrix kLoadOnX2C{{0.0}, {1.0}};

/// A quiet NaN from its bits, read through a volatile so the compiler cannot fold it.
double runtime_nan() {
    volatile std::uint64_t bits = 0x7ff8000000000000ull;
    return std::bit_cast<double>(static_cast<std::uint64_t>(bits));
}

/// Whether two targets hold the same maps of d and r.
template<typename A, typename B>
bool same_target(const A& a, const B& b) {
    return same_map(a.target_d, b.target_d) && same_map(a.target_r, b.target_r);
}

} // namespace

TEST_SUITE("MeasuredDisturbance") {
    TEST_CASE("an LQR holds its output at r under a measured load with no offset; the same LQR without it does not") {
        const auto md = design::state_feedback(*kLags, kWeights, design::TrackedOutputs{{1}}, design::MeasuredDisturbance{kLoadOnX2});
        const auto plain = design::state_feedback(*kLags, kWeights, design::TrackedOutputs{{1}});
        REQUIRE(md.success);
        REQUIRE(plain.success);
        CHECK(same_map(md.K, plain.K)); // the gain does not depend on d
        static_assert(std::is_same_v<decltype(md.target.target_d), Matrix<3, 1>>);

        //! The target is the MPC's: steady_state_target's under d, on the same held output.
        const auto same = design::steady_state_target(*kLags, design::TrackedOutputs{{1}}, design::MeasuredDisturbance{kLoadOnX2});
        CHECK(same_map(md.target.target_d, same.target_d));
        CHECK(same_map(md.target.target_r, same.target_r));
        const auto mpc = design::mpc<20, 5>(*kLags, both_weighed(), design::MeasuredDisturbance{kLoadOnX2}, design::TrackedOutputs{{1}});
        REQUIRE(mpc.success);
        CHECK(same_map(md.target.target_d, mpc.target.target_d));

        LQR          law{md};
        LQR          blind{plain};
        ColVec<2>    x{};
        ColVec<2>    xb{};
        const double w = 0.4;
        for (size_t k = 0; k < 400; ++k) {
            x = lags_step(x, law.step(0.7, ColVec<1>{w}, x)(0), w);
            xb = lags_step(xb, blind.step(0.7, xb)(0), w);
        }
        CHECK(x(1) == doctest::Approx(0.7).epsilon(1e-9));
        const double u_t = (md.target.target_r(2, 0) * 0.7) + (md.target.target_d(2, 0) * w);
        CHECK(law.control(0.7, ColVec<1>{w})(0) == doctest::Approx(u_t).epsilon(1e-9)); // the steady input is the target's
        CHECK(std::abs(xb(1) - 0.7) > 0.1);                                             // without it: the load's offset
    }

    TEST_CASE("the regulator holds every output at 0 under d: control(d) and step(d, x), and its feedforward form") {
        //! One reading, x₂: the regulator's target holds it at 0 under the load.
        const auto md = design::state_feedback(kSecond, kWeights, design::MeasuredDisturbance{kLoadOnX2});
        REQUIRE(md.success);
        static_assert(std::is_same_v<decltype(md.target.target_r), Matrix<3, 0>>); // no r
        LQR       law{md};
        LQR       blind{design::state_feedback(kSecond, kWeights)};
        ColVec<2> x{};
        ColVec<2> xb{};
        for (size_t k = 0; k < 400; ++k) {
            x = lags_step(x, law.step(ColVec{0.4}, x)(0), 0.4);
            xb = lags_step(xb, blind.step(xb)(0), 0.4);
        }
        CHECK(std::abs(x(1)) < 1e-9);
        CHECK(std::abs(xb(1)) > 0.05);

        //! About a planned point: u = u_ref + u_t − K·(x − x_ref − x_t).
        const ColVec x_ref{0.1, -0.2};
        const ColVec u_ref{0.3};
        const ColVec d{0.25};
        law.update(x);
        const ColVec<3> t = md.target.target_d * d;
        const ColVec<2> dev{x(0) - x_ref(0) - t(0), x(1) - x_ref(1) - t(1)};
        CHECK(law.control(x_ref, u_ref, d)(0) == doctest::Approx(u_ref(0) + t(2) - (md.K * dev)(0)).epsilon(1e-12));
    }

    TEST_CASE("a regulator with more outputs than inputs lists the one it holds under d: TrackedOutputs with r = 0") {
        //! At rest x₂ = x₁ + d: no steady state holds both readings at 0. Holding both took the
        //! least squares, x = (−0.5·d, 0.5·d), which held neither, and no longer compiles.
        //! Listed, x₂ is held at r = 0, and x₁ sits at −d.
        constexpr Matrix load_c{{0.0}, {1.0}};
        const auto       one = design::state_feedback(kLagsC, kTs, kWeights, design::TrackedOutputs{{1}}, design::MeasuredDisturbance{load_c});
        REQUIRE(one.success);
        CHECK(one.target.target_d(0, 0) == doctest::Approx(-1.0).epsilon(1e-9));
        CHECK(std::abs(one.target.target_d(1, 0)) < 1e-9);
        const auto load = discretize(StateSpace<2, 1, 1>{.A = kLagsC.A, .B = load_c}, kTs).value();
        LQR        law{one};
        ColVec<2>  x{};
        for (size_t k = 0; k < 400; ++k) {
            x = (kLags->A * x) + (kLags->B * law.step(0.0, ColVec{0.4}, x)) + (load.B * ColVec{0.4});
        }
        CHECK(std::abs(x(1)) < 1e-9);                       // the listed reading at 0
        CHECK(x(0) == doctest::Approx(-0.4).epsilon(1e-9)); // the other where the steady state puts it
    }

    TEST_CASE("an LQI feeds the measured load forward: the same steady state with its integrators where they were") {
        const design::IntegralAction action{.Q = 10.0, .outputs = {1}};
        const auto                   md = design::state_feedback(*kLags, kWeights, action, design::MeasuredDisturbance{kLoadOnX2});
        const auto                   plain = design::state_feedback(*kLags, kWeights, action);
        REQUIRE(md.success);
        REQUIRE(plain.success);
        CHECK(same_map(md.K, plain.K));
        LQI          law{md};
        LQI          blind{plain};
        ColVec<2>    x{};
        ColVec<2>    xb{};
        double       xi_before = 0.0;
        double       peak = 0.0;
        double       peak_blind = 0.0;
        const double r = 0.7;
        for (size_t k = 0; k < 900; ++k) {
            const double w = (k < 300) ? 0.0 : 0.4;
            if (k == 300) {
                xi_before = law.integral.xi(0);
            }
            x = lags_step(x, law.step(r, ColVec<1>{w}, x)(0), w);
            xb = lags_step(xb, blind.step(r, xb)(0), w);
            if (k >= 300) {
                peak = std::max(peak, std::abs(x(1) - r));
                peak_blind = std::max(peak_blind, std::abs(xb(1) - r));
            }
        }
        CHECK(x(1) == doctest::Approx(r).epsilon(1e-9));
        CHECK(xb(1) == doctest::Approx(r).epsilon(1e-6));                      // integral action holds it too, in time
        CHECK(law.integral.xi(0) == doctest::Approx(xi_before).epsilon(1e-9)); // the feedforward carried the load, not the integrators
        CHECK(peak < 0.5 * peak_blind);
    }

    TEST_CASE("an LQG takes the measured d as its filter's input: no bias in the estimate, no offset") {
        const auto md = design::state_feedback(kNoisyLags, kWeights, kNoise, design::TrackedOutputs{{1}}, design::MeasuredDisturbance{kLoadOnX2});
        const auto plain = design::state_feedback(kNoisyLags, kWeights, kNoise, design::TrackedOutputs{{1}});
        REQUIRE(md.success);
        REQUIRE(plain.success);
        //! The filter's model: B = [B, B_d], D = [D, 0]; its gain the one without d, which is not noise.
        static_assert(std::is_same_v<decltype(md.kalman.sys.B), Matrix<2, 2>>);
        CHECK(same_map(md.kalman.M, plain.kalman.M));
        CHECK(md.kalman.sys.B(1, 1) == kLoadOnX2(1, 0));
        CHECK(md.kalman.sys.B(0, 0) == plain.kalman.sys.B(0, 0));
        LQG          law{md};
        LQG          blind{plain};
        ColVec<2>    x{};
        ColVec<2>    xb{};
        const double w = 0.4;
        for (size_t k = 0; k < 600; ++k) {
            x = lags_step(x, law.step(0.7, ColVec<1>{w}, x)(0), w);
            xb = lags_step(xb, blind.step(0.7, xb)(0), w);
        }
        CHECK(x(1) == doctest::Approx(0.7).epsilon(1e-9));
        CHECK(Matrix<2, 1>(law.state_estimate() - x).norm() < 1e-9);
        CHECK(std::abs(xb(1) - 0.7) > 0.05);
        CHECK(Matrix<2, 1>(blind.state_estimate() - xb).norm() > 1e-4); // the load read as a model error
    }

    TEST_CASE("an OFLQG with a measured and an estimated disturbance: d̂ settles at the unmeasured part only") {
        const design::DisturbanceModel spec{.Bd = kLoadOnX2, .outputs = {1}};
        const auto                     md = design::state_feedback(kNoisyLags, kWeights, kLoadNoise, spec, design::MeasuredDisturbance{kLoadOnX2});
        const auto                     plain = design::state_feedback(kNoisyLags, kWeights, kLoadNoise, spec);
        REQUIRE(md.success);
        REQUIRE(plain.success);
        //! The target is under [d̂; d]: the estimated column the OFLQG's own, the measured one the same channel's.
        static_assert(std::is_same_v<decltype(md.target.target_d), Matrix<3, 2>>);
        const Matrix<3, 1> est = md.target.target_d.block<3, 1>(0, 0).to_matrix();
        const Matrix<3, 1> mea = md.target.target_d.block<3, 1>(0, 1).to_matrix();
        CHECK(same_map(est, plain.target.target_d));
        CHECK(same_map(mea, plain.target.target_d));

        OFLQG        law{md};
        OFLQG        blind{plain};
        ColVec<2>    x{};
        ColVec<2>    xb{};
        const double w_measured = 0.4;
        const double w_unmeasured = 0.15;
        for (size_t k = 0; k < 800; ++k) {
            const double w = w_measured + w_unmeasured;
            x = lags_step(x, law.step(0.7, ColVec<1>{w_measured}, x)(0), w);
            xb = lags_step(xb, blind.step(0.7, xb)(0), w);
        }
        CHECK(x(1) == doctest::Approx(0.7).epsilon(1e-9));
        CHECK(law.disturbance_estimate()(0) == doctest::Approx(w_unmeasured).epsilon(1e-7));
        CHECK(blind.disturbance_estimate()(0) == doctest::Approx(w_measured + w_unmeasured).epsilon(1e-7));
        CHECK(law.last_command()(0) == doctest::Approx(blind.last_command()(0)).epsilon(1e-7)); // the same steady input
    }

    TEST_CASE("an MPC with a KalmanNoise and a measured d, refused before, settles where the measured-state MPC does") {
        const design::DisturbanceModel spec{.Bd = kLoadOnX2, .outputs = {1}};
        const auto                     of = design::mpc<20, 5>(kNoisyLags, both_weighed(), kLoadNoise, spec, design::MeasuredDisturbance{kLoadOnX2});
        const auto                     ms = design::mpc<20, 5>(*kLags, both_weighed(), design::MeasuredDisturbance{kLoadOnX2}, design::TrackedOutputs{{1}});
        REQUIRE(of.success);
        REQUIRE(ms.success);
        static_assert(std::is_same_v<decltype(of.target.target_d), Matrix<3, 2>>); // [d̂; d]
        CHECK(same_map(of.target.target_d.block<3, 1>(0, 1).to_matrix(), ms.target.target_d));
        CHECK(same_map(of.target.target_r, ms.target.target_r));

        OFMPC filtered{of};
        MPC   measured{ms};
        static_assert(decltype(filtered)::NM == 1);
        ColVec<2> xf{};
        ColVec<2> xm{};
        for (size_t k = 0; k < 400; ++k) {
            const double w = (k < 100) ? 0.0 : 0.4;
            xf = lags_step(xf, filtered.step(0.7, ColVec<1>{w}, xf)(0), w);
            xm = lags_step(xm, measured.step(0.7, ColVec<1>{w}, xm)(0), w);
        }
        CHECK(xf(1) == doctest::Approx(0.7).epsilon(1e-9));
        CHECK(Matrix<2, 1>(xf - xm).norm() < 1e-9);
        CHECK(filtered.applied()(0) == doctest::Approx(measured.applied()(0)).epsilon(1e-9));
        CHECK(std::abs(filtered.disturbance_estimate()(0)) < 1e-9); // every load measured: nothing left to estimate

        //! With a part it does not measure, d̂ takes that part.
        OFMPC     partial{of};
        ColVec<2> xp{};
        for (size_t k = 0; k < 600; ++k) {
            xp = lags_step(xp, partial.step(0.7, ColVec{0.4}, xp)(0), 0.55);
        }
        CHECK(xp(1) == doctest::Approx(0.7).epsilon(1e-9));
        CHECK(partial.disturbance_estimate()(0) == doctest::Approx(0.15).epsilon(1e-7));
    }

    TEST_CASE("an unconstrained MPC with the DARE's terminal weight is the LQR under the same measured d") {
        //! Both are u = u_t(r, d) − K·(x − x_t(r, d)), the MPC holding d over its horizon: for a d
        //! held constant the target alone is the LQ-optimal feedforward.
        constexpr StateSpace<2, 1, 2>     di{.A = {{1.0, kTs}, {0.0, 1.0}}, .B = {{0.5 * kTs * kTs}, {kTs}}, .C = Matrix<2, 2>::identity(), .Ts = kTs};
        const Matrix<2, 2>                Q = diag(1.0, 0.1);
        const Matrix                      R{{0.5}};
        const design::MeasuredDisturbance load{.Bd = di.B};
        const auto                        lq = design::state_feedback(di, {.Q = Q, .R = R}, design::TrackedOutputs{{0}}, load);
        REQUIRE(lq.success);
        design::MPCWeights<1, 2> w{.Qy = Q, .Rdu = Matrix<1, 1>{}, .Ru = R};
        w.Qy_terminal = design::discrete_lqr(di.A, di.B, Q, R).S;
        const auto art = design::mpc<8, 8>(di, w, design::TrackedOutputs{{0}}, load);
        REQUIRE(art.success);
        LQR    lqr{lq};
        MPC    mpc{art};
        ColVec x{1.0, -0.5};
        for (size_t k = 0; k < 40; ++k) {
            const double    r = (k < 20) ? 0.0 : 0.3;
            const ColVec<1> d{(k < 10) ? 0.0 : 0.2};
            const ColVec<2> x_ref{0.05 * std::sin(0.2 * static_cast<double>(k)), 0.0};
            const ColVec    u_ref{0.02};
            lqr.update(x);
            const double u_lqr = lqr.control(r, x_ref, u_ref, d)(0);
            const double u_mpc = mpc.step(r, x_ref, u_ref, d, x)(0);
            CHECK(u_mpc == doctest::Approx(u_lqr).epsilon(1e-6));
            x = (di.A * x) + (di.B * ColVec<1>{u_mpc + d(0)});
        }
    }

    TEST_CASE("the measured forms are predict, update and control, and the filter takes [u; d] by hand") {
        const auto md = design::state_feedback(kNoisyLags, kWeights, kNoise, design::TrackedOutputs{{1}}, design::MeasuredDisturbance{kLoadOnX2});
        REQUIRE(md.success);
        LQG                                law{md};
        LQG                                split{md};
        FixedGainObserver<2, 2, 2, double> obs{md.kalman};
        ColVec                             x{0.2, -0.1};
        double                             u_prev = 0.0;
        double                             d_prev = 0.0;
        for (size_t k = 0; k < 30; ++k) {
            const double    r = (k < 10) ? 0.0 : 0.5;
            const ColVec<1> d{0.3 * std::sin(0.3 * static_cast<double>(k))};
            const double    u = law.step(r, d, x)(0);
            split.predict();
            split.update(x);
            CHECK(split.control(r, d)(0) == doctest::Approx(u).epsilon(1e-12));
            //! By hand: predict with the input and the d the plant held, update, then the law about the target under d.
            obs.predict(ColVec<2>{u_prev, d_prev});
            obs.update(x);
            const ColVec<3> t = (md.target.target_r * ColVec<1>{r}) + (md.target.target_d * d);
            const ColVec<2> dev{obs.state()(0) - t(0), obs.state()(1) - t(1)};
            CHECK(u == doctest::Approx(t(2) - (md.K * dev)(0)).epsilon(1e-12));
            u_prev = u;
            d_prev = d(0);
            x = lags_step(x, u, d(0));
        }

        //! The OFLQG and the OFMPC the same way.
        const design::DisturbanceModel spec{.Bd = kLoadOnX2, .outputs = {1}};
        const auto                     of = design::state_feedback(kNoisyLags, kWeights, kLoadNoise, spec, design::MeasuredDisturbance{kLoadOnX2});
        const auto                     ofm = design::mpc<10, 4>(kNoisyLags, both_weighed(), kLoadNoise, spec, design::MeasuredDisturbance{kLoadOnX2});
        REQUIRE(of.success);
        REQUIRE(ofm.success);
        OFLQG     a{of};
        OFLQG     b{of};
        OFMPC     c{ofm};
        OFMPC     e{ofm};
        ColVec<2> z{};
        for (size_t k = 0; k < 30; ++k) {
            const ColVec<1> d{0.2 + (0.01 * static_cast<double>(k))};
            const ColVec<2> x_ref{0.0, 0.01 * static_cast<double>(k)};
            const ColVec    u_ref{0.01};
            const double    u = a.step(0.6, x_ref, u_ref, d, z)(0);
            b.predict();
            b.update(z);
            CHECK(b.control(0.6, x_ref, u_ref, d)(0) == doctest::Approx(u).epsilon(1e-12));
            const double v = c.step(0.6, d, z)(0);
            e.predict();
            e.update(z);
            CHECK(e.control(0.6, d)(0) == doctest::Approx(v).epsilon(1e-12));
            z = lags_step(z, u, d(0) + 0.05);
        }
    }

    TEST_CASE("to_ss is the law at d = 0, its filter on the plant's inputs: the block of the design without it") {
        const auto                     lqg = design::state_feedback(kNoisyLags, kWeights, kNoise, design::TrackedOutputs{{1}}, design::MeasuredDisturbance{kLoadOnX2});
        const auto                     lqg_plain = design::state_feedback(kNoisyLags, kWeights, kNoise, design::TrackedOutputs{{1}});
        const design::DisturbanceModel spec{.Bd = kLoadOnX2, .outputs = {1}};
        const auto                     of = design::state_feedback(kNoisyLags, kWeights, kLoadNoise, spec, design::MeasuredDisturbance{kLoadOnX2});
        const auto                     of_plain = design::state_feedback(kNoisyLags, kWeights, kLoadNoise, spec);
        const auto                     lqgi = design::state_feedback(kNoisyLags, kWeights, kNoise, design::IntegralAction{.Q = 10.0, .outputs = {1}}, design::MeasuredDisturbance{kLoadOnX2});
        const auto                     lqgi_plain = design::state_feedback(kNoisyLags, kWeights, kNoise, design::IntegralAction{.Q = 10.0, .outputs = {1}});
        REQUIRE((lqg.success && of.success && lqgi.success));
        const auto a = lqg.to_ss();
        const auto b = lqg_plain.to_ss();
        CHECK((same_map(a.A, b.A) && same_map(a.B, b.B) && same_map(a.C, b.C) && same_map(a.D, b.D)));
        const auto c = of.to_ss();
        const auto e = of_plain.to_ss();
        CHECK(same_map(c.A, e.A));
        CHECK(same_map(c.C, e.C));
        const auto g = lqgi.to_ss();
        const auto h = lqgi_plain.to_ss();
        CHECK((same_map(g.A, h.A) && same_map(g.B, h.B) && same_map(g.C, h.C) && same_map(g.D, h.D)));
    }

    TEST_CASE("with an InputDelay d enters the plant, not the delay states, and a late reading as the plant moved under it") {
        constexpr Matrix load_c{{0.0}, {1.0}}; // continuous: x₂' gains w
        const auto       d = discretize(StateSpace<2, 1, 1>{.A = kLagsC.A, .B = load_c}, kTs).value();
        const auto       lqg = design::state_feedback(kLagsC, kTs, kWeights, kNoise, design::TrackedOutputs{{1}}, design::MeasuredDisturbance{load_c}, InputDelay<2>{0.15});
        REQUIRE(lqg.success);
        static_assert(std::is_same_v<decltype(lqg.kalman.sys.B), Matrix<4, 2>>); // [x; v(k−1); v(k−2)], then [u, d]
        CHECK(lqg.kalman.sys.B(0, 1) == doctest::Approx(d.B(0, 0)).epsilon(1e-14));
        CHECK(lqg.kalman.sys.B(1, 1) == doctest::Approx(d.B(1, 0)).epsilon(1e-14));
        CHECK(lqg.kalman.sys.B(2, 1) == 0.0);
        CHECK(lqg.kalman.sys.B(3, 1) == 0.0);
        //! At the target under d each delay state holds u_t.
        CHECK(lqg.target.target_d(2, 0) == doctest::Approx(lqg.target.target_d(4, 0)).epsilon(1e-12));
        CHECK(lqg.target.target_d(3, 0) == doctest::Approx(lqg.target.target_d(4, 0)).epsilon(1e-12));

        //! A measured-state LQR on the delayed plant holds the output under the load.
        const auto lqr = design::state_feedback(kLagsC, kTs, kWeights, design::TrackedOutputs{{1}}, design::MeasuredDisturbance{load_c}, InputDelay<2>{0.15});
        REQUIRE(lqr.success);
        const auto model = discretize(kLagsC, kTs, InputDelay<2>{0.15}).value();
        LQR        law{lqr};
        ColVec<4>  s{};
        for (size_t k = 0; k < 600; ++k) {
            const double u = law.step(0.7, ColVec{0.3}, ColVec<2>{s(0), s(1)})(0);
            ColVec<4>    load{};
            load(0) = d.B(0, 0) * 0.3;
            load(1) = d.B(1, 0) * 0.3;
            s = (model.A * s) + (model.B * ColVec<1>{u}) + load;
        }
        CHECK(s(1) == doctest::Approx(0.7).epsilon(1e-9));

        //! A late reading: its delay state at the target is the steady output, C·x_t.
        const auto late = design::state_feedback(kLagsC, kTs, kWeights, kNoise, design::TrackedOutputs{{1}}, design::MeasuredDisturbance{load_c}, OutputDelay{0.04});
        REQUIRE(late.success);
        static_assert(std::is_same_v<decltype(late.target.target_d), Matrix<5, 1>>); // [x; z] then u
        CHECK(late.target.target_d(2, 0) == doctest::Approx(late.target.target_d(0, 0)).epsilon(1e-12));
        CHECK(late.target.target_d(3, 0) == doctest::Approx(late.target.target_d(1, 0)).epsilon(1e-12));
    }

    TEST_CASE("design::kalman takes a MeasuredDisturbance as known inputs after u, its gain unchanged") {
        const auto kf = design::kalman(kNoisyLags, kNoise, design::MeasuredDisturbance{kLoadOnX2});
        const auto plain = design::kalman(kNoisyLags, kNoise);
        REQUIRE(kf.success);
        static_assert(std::is_same_v<decltype(kf), const design::KalmanResult<2, 2, 2, double, 2, 2, 0>>);
        CHECK(same_map(kf.M, plain.M));
        CHECK(same_map(kf.P, plain.P));
        CHECK(kf.sys.B(1, 1) == kLoadOnX2(1, 0));
        CHECK(kf.sys.D(0, 1) == 0.0);
        //! With a disturbance spec: d's rows of the measured column are zero, d̂ does not move with d.
        const auto both = design::kalman(kNoisyLags, kLoadNoise, design::InputDisturbance{}, design::MeasuredDisturbance{kLoadOnX2});
        REQUIRE(both.success);
        static_assert(std::is_same_v<decltype(both.sys.B), Matrix<3, 2>>);
        CHECK(both.sys.B(2, 1) == 0.0);
        CHECK(both.sys.B(1, 1) == kLoadOnX2(1, 0));
    }

    TEST_CASE("a non-finite B_d is refused, as a non-finite channel of a disturbance spec is") {
        Matrix<2, 1> bad = kLoadOnX2;
        bad(1, 0) = std::numeric_limits<double>::quiet_NaN();
        CHECK(design::state_feedback(*kLags, kWeights, design::TrackedOutputs{{1}}, design::MeasuredDisturbance{bad}).status == design::DesignStatus::bad_weights);
        CHECK(design::state_feedback(kNoisyLags, kWeights, kNoise, design::TrackedOutputs{{1}}, design::MeasuredDisturbance{bad}).status == design::DesignStatus::bad_weights);
        CHECK(design::kalman(kNoisyLags, kNoise, design::MeasuredDisturbance{bad}).status == design::DesignStatus::bad_weights);
        CHECK(design::mpc<10, 4>(kNoisyLags, both_weighed(), kLoadNoise, design::InputDisturbance{{1}}, design::MeasuredDisturbance{bad}).status == design::DesignStatus::bad_weights);
        CHECK(design::mpc<10, 4>(*kLags, both_weighed(), design::TrackedOutputs{{1}}, design::MeasuredDisturbance{bad}).status == design::DesignStatus::bad_weights);
    }

    TEST_CASE("Law<kRes, T> of a constexpr design with a measured disturbance runs on its constants as the runtime does") {
        static constexpr auto kOflqg = design::state_feedback(kNoisyLags, kWeights, kLoadNoise, design::DisturbanceModel{.Bd = kLoadOnX2, .outputs = {1}}, design::MeasuredDisturbance{kLoadOnX2});
        static constexpr auto kLqr = design::state_feedback(kSecond, kWeights, design::MeasuredDisturbance{kLoadOnX2});
        static constexpr auto kOfmpc = design::mpc<6, 2>(kNoisyLags, both_weighed(), kLoadNoise, design::InputDisturbance{{1}}, design::MeasuredDisturbance{kLoadOnX2});
        static_assert(kOflqg.success && kLqr.success && kOfmpc.success);
        Law<kOfmpc, float> f;
        OFMPC              g{kOfmpc.as<float>()};
        Law<kOflqg, float> a;
        OFLQG              b{kOflqg.as<float>()};
        Law<kLqr, float>   c;
        LQR                e{kLqr.as<float>()};
        ColVec<2, float>   x{};
        for (size_t k = 0; k < 50; ++k) {
            const ColVec d{0.3f};
            const float  u = a.step(0.5f, d, x)(0);
            CHECK(static_cast<double>(u) == doctest::Approx(static_cast<double>(b.step(0.5f, d, x)(0))).epsilon(1e-5));
            CHECK(static_cast<double>(c.step(d, x)(0)) == doctest::Approx(static_cast<double>(e.step(d, x)(0))).epsilon(1e-5));
            CHECK(static_cast<double>(f.step(0.5f, d, x)(0)) == doctest::Approx(static_cast<double>(g.step(0.5f, d, x)(0))).epsilon(1e-4));
            const ColVec<2> next = lags_step(ColVec<2>{static_cast<double>(x(0)), static_cast<double>(x(1))}, static_cast<double>(u), 0.4);
            x = ColVec<2, float>{static_cast<float>(next(0)), static_cast<float>(next(1))};
        }
    }

    TEST_CASE("a law without a MeasuredDisturbance is the one it was: its type, its filter's inputs and its size") {
        const auto lqr = design::state_feedback(*kLags, kWeights);
        const auto lqrt = design::state_feedback(*kLags, kWeights, design::TrackedOutputs{{1}});
        const auto lqi = design::state_feedback(*kLags, kWeights, design::IntegralAction{.Q = Matrix<2, 2>::identity()});
        const auto lqg = design::state_feedback(kNoisyLags, kWeights, kNoise);
        const auto lqgt = design::state_feedback(kNoisyLags, kWeights, kNoise, design::TrackedOutputs{{1}});
        const auto lqgi = design::state_feedback(kNoisyLags, kWeights, kNoise, design::IntegralAction{.Q = Matrix<2, 2>::identity()});
        const auto oflqg = design::state_feedback(kNoisyLags, kWeights, kLoadNoise, design::DisturbanceModel{.Bd = kLoadOnX2, .outputs = {1}});
        const auto mpc = design::mpc<10, 4>(*kLags, {}, design::TrackedOutputs{{1}});
        const auto mpcm = design::mpc<10, 4>(*kLags, {}, design::TrackedOutputs{{1}}, design::MeasuredDisturbance{kLoadOnX2});
        const auto ofmpc = design::mpc<10, 4>(kNoisyLags, {}, kLoadNoise, design::DisturbanceModel{.Bd = kLoadOnX2, .outputs = {1}});
        static_assert(std::is_same_v<decltype(lqr), const design::StateFeedbackResult<2, 1, double>>);
        static_assert(std::is_same_v<decltype(LQR{lqr}), LQR<2, 1, double>>);
        static_assert(std::is_same_v<decltype(LQG{lqg}), LQG<2, 1, 2, double, FixedGainObserver<2, 1, 2, double>>>); // no input column for d
        static_assert(std::is_same_v<decltype(oflqg.kalman), design::KalmanResult<3, 1, 2, double, 3, 2, 1>>);
        static_assert(decltype(OFMPC{ofmpc})::NM == 0);
        //! The sizes on the base this was built on (x86-64): no member, no padding for d. The MPCs'
        //! are those without Δu rows (2992, 3184 and 3504 with them, before InputRateLimits), and the
        //! results' hold their limits as one StateFeedbackLimits, with a two-byte empty box for each
        //! limit the law lacks (144 and 752 before).
        if constexpr (sizeof(size_t) == 8) {
            CHECK(sizeof(LQR{lqr}) == 88);
            CHECK(sizeof(LQR{lqr.as<float>()}) == 48);
            CHECK(sizeof(LQR{lqrt}) == 120);
            CHECK(sizeof(LQI{lqi}) == 288);
            CHECK(sizeof(LQG{lqg}) == 288);
            CHECK(sizeof(LQG{lqgt}) == 320);
            CHECK(sizeof(LQGI{lqgi}) == 472);
            CHECK(sizeof(OFLQG{oflqg}) == 488);
            CHECK(sizeof(OFLQG{oflqg.as<float>()}) == 256);
            CHECK(sizeof(MPC{mpc}) == 2368); // one output held, TrackedOutputs{{1}}: 2376 holding both
            CHECK(sizeof(MPC{mpcm}) == 2496);
            CHECK(sizeof(OFMPC{ofmpc}) == 2824);
            CHECK(sizeof(oflqg) == 776);
            CHECK(sizeof(lqr) == 160);
        }
        //! A measured d on a measured state adds the target's map of it and nothing to keep between ticks.
        const auto lqrm = design::state_feedback(*kLags, kWeights, design::TrackedOutputs{{1}}, design::MeasuredDisturbance{kLoadOnX2});
        CHECK(sizeof(LQR{lqrm}) <= sizeof(LQR{lqrt}) + (3 * sizeof(double)));
    }

    TEST_CASE("steady_state_target holds the outputs a TrackedOutputs lists: the LQR's and the MPC's target, with and without a MeasuredDisturbance") {
        const design::TrackedOutputs      held{{1}}; // one output of the two
        const design::MeasuredDisturbance load{kLoadOnX2};
        const auto                        tg = design::steady_state_target(*kLags, held);
        const auto                        lqr = design::state_feedback(*kLags, kWeights, held);
        const auto                        mpc = design::mpc<20, 5>(*kLags, both_weighed(), held);
        REQUIRE(tg.success);
        REQUIRE(lqr.success);
        REQUIRE(mpc.success);
        static_assert(std::is_same_v<decltype(tg), const design::SteadyStateTargetResult<2, 1, 0, 1, double>>);
        static_assert(std::is_same_v<decltype(tg), const decltype(lqr.target)>); // the laws' own target
        static_assert(std::is_same_v<decltype(tg), const decltype(mpc.target)>);
        CHECK(same_target(tg, lqr.target));
        CHECK(same_target(tg, mpc.target));
        CHECK(tg.target_r(1, 0) == doctest::Approx(1.0).epsilon(1e-12)); // x₂ at r: the listed output

        //! Under a measured d: d's column after the (empty) spec's, r's map the same.
        const auto tgm = design::steady_state_target(*kLags, held, load);
        const auto lqrm = design::state_feedback(*kLags, kWeights, held, load);
        const auto mpcm = design::mpc<20, 5>(*kLags, both_weighed(), held, load);
        REQUIRE(tgm.success);
        REQUIRE(lqrm.success);
        REQUIRE(mpcm.success);
        static_assert(std::is_same_v<decltype(tgm), const decltype(lqrm.target)>);
        static_assert(std::is_same_v<decltype(tgm), const decltype(mpcm.target)>);
        CHECK(same_target(tgm, lqrm.target));
        CHECK(same_target(tgm, mpcm.target));
        CHECK(same_map(tgm.target_r, tg.target_r));

        //! TrackedOutputs{}, every output, on a plant with no more outputs than inputs: the list of
        //! all of them, and the LQR's and the MPC's with TrackedOutputs{}.
        const auto every = design::steady_state_target(kSecond, design::TrackedOutputs{});
        static_assert(std::is_same_v<decltype(every), const decltype(design::steady_state_target(kSecond, design::TrackedOutputs{{0}}))>);
        CHECK(same_target(every, design::steady_state_target(kSecond, design::TrackedOutputs{{0}})));
        CHECK(same_target(every, design::state_feedback(kSecond, kWeights, design::TrackedOutputs{}).target));
        CHECK(same_target(every, design::mpc<10, 4>(kSecond, {}, design::TrackedOutputs{}).target));

        //! The steady state it gives under d is the one both laws hold, and their command there.
        const SteadyStateTarget target{tgm};
        const double            r = 0.7;
        const double            w = 0.4;
        const auto              p = target.equilibrium(ColVec<1>{r}, ColVec<1>{w});
        LQR                     law{lqrm};
        MPC                     ctl{mpcm};
        ColVec<2>               x{};
        ColVec<2>               xm{};
        for (size_t k = 0; k < 400; ++k) {
            x = lags_step(x, law.step(r, ColVec<1>{w}, x)(0), w);
            xm = lags_step(xm, ctl.step(r, ColVec<1>{w}, xm)(0), w);
        }
        CHECK(Matrix<2, 1>(x - p.x).norm() < 1e-9);
        CHECK(Matrix<2, 1>(xm - p.x).norm() < 1e-9);
        CHECK(law.control(r, ColVec<1>{w})(0) == doctest::Approx(p.u(0)).epsilon(1e-9));
        CHECK(ctl.applied()(0) == doctest::Approx(p.u(0)).epsilon(1e-9));
    }

    TEST_CASE("with a disturbance spec, steady_state_target is under [d̂; d], the OFLQG's and the OFMPC's target") {
        const design::DisturbanceModel    spec{.Bd = kLoadOnX2, .outputs = {1}};
        const design::MeasuredDisturbance load{kLoadOnX2};
        const auto                        tg = design::steady_state_target(kNoisyLags, spec, load);
        const auto                        of = design::state_feedback(kNoisyLags, kWeights, kLoadNoise, spec, load);
        const auto                        ofm = design::mpc<20, 5>(kNoisyLags, both_weighed(), kLoadNoise, spec, load);
        REQUIRE(tg.success);
        REQUIRE(of.success);
        REQUIRE(ofm.success);
        static_assert(std::is_same_v<decltype(tg.target_d), Matrix<3, 2>>); // d̂'s column, then d's
        CHECK(same_target(tg, of.target));
        CHECK(same_target(tg, ofm.target));
        //! Its measured column is the target's without the spec, on the output the spec holds, and a spec's own on that channel.
        const Matrix<3, 1> measured = tg.target_d.block<3, 1>(0, 1).to_matrix();
        CHECK(same_map(measured, design::steady_state_target(*kLags, design::TrackedOutputs{{1}}, load).target_d));
        CHECK(same_map(measured, design::steady_state_target(*kLags, design::DisturbanceModel{.Bd = kLoadOnX2, .outputs = {1}}).target_d));

        //! The OFLQG settles at its steady state under [d̂; d], d̂ the part not measured.
        const SteadyStateTarget target{tg};
        OFLQG                   law{of};
        ColVec<2>               x{};
        for (size_t k = 0; k < 800; ++k) {
            x = lags_step(x, law.step(0.7, ColVec{0.4}, x)(0), 0.55);
        }
        const auto p = target.equilibrium(ColVec{0.7}, ColVec<2>{law.disturbance_estimate()(0), 0.4});
        CHECK(Matrix<2, 1>(x - p.x).norm() < 1e-9);
        CHECK(law.last_command()(0) == doctest::Approx(p.u(0)).epsilon(1e-9));
    }

    TEST_CASE("steady_state_target with an InputDelay: the laws' target on the delayed model, d on the plant, not the delay states") {
        const design::TrackedOutputs      held{{1}};
        const design::MeasuredDisturbance load{kLoadOnX2C};
        const InputDelay<2>               late{0.15};
        const auto                        plain = design::steady_state_target(kLagsC, kTs, held, late);
        REQUIRE(plain.success);
        CHECK(same_target(plain, design::state_feedback(kLagsC, kTs, kWeights, held, late).target));
        CHECK(same_target(plain, design::mpc<20, 5>(kLagsC, kTs, both_weighed(), held, late).target));

        const auto tg = design::steady_state_target(kLagsC, kTs, held, load, late);
        const auto lqr = design::state_feedback(kLagsC, kTs, kWeights, held, load, late);
        const auto mpc = design::mpc<20, 5>(kLagsC, kTs, both_weighed(), held, load, late);
        REQUIRE(tg.success);
        REQUIRE(lqr.success);
        REQUIRE(mpc.success);
        static_assert(std::is_same_v<decltype(tg.target_d), Matrix<5, 1>>); // [x; v(k−1); v(k−2)] then u
        static_assert(std::is_same_v<decltype(tg), const decltype(lqr.target)>);
        CHECK(same_target(tg, lqr.target));
        CHECK(same_target(tg, mpc.target));
        //! Each delay state holds u_t under d: d does not travel through them.
        CHECK(tg.target_d(2, 0) == doctest::Approx(tg.target_d(4, 0)).epsilon(1e-12));
        CHECK(tg.target_d(3, 0) == doctest::Approx(tg.target_d(4, 0)).epsilon(1e-12));

        //! With an input disturbance estimated beside it, holding the same output: the OFLQG's and the OFMPC's.
        const design::InputDisturbance spec{{1}};
        const auto                     spec_tg = design::steady_state_target(kLagsC, kTs, spec, load, late);
        const auto                     of = design::state_feedback(kLagsC, kTs, kWeights, kLoadNoise, spec, load, late);
        const auto                     ofm = design::mpc<20, 5>(kLagsC, kTs, both_weighed(), kLoadNoise, spec, load, late);
        REQUIRE(spec_tg.success);
        REQUIRE(of.success);
        REQUIRE(ofm.success);
        static_assert(std::is_same_v<decltype(spec_tg.target_d), Matrix<5, 2>>);
        CHECK(same_target(spec_tg, of.target));
        CHECK(same_target(spec_tg, ofm.target));
    }

    TEST_CASE("steady_state_target with an OutputDelay: the filter's model, d on the late readings' delay states") {
        const design::TrackedOutputs      held{{1}};
        const design::MeasuredDisturbance load{kLoadOnX2C};
        const OutputDelay<1>              reading{0.04};
        const auto                        tg = design::steady_state_target(kLagsC, kTs, held, load, reading);
        const auto                        lqg = design::state_feedback(kLagsC, kTs, kWeights, kNoise, held, load, reading);
        const auto                        mpc = design::mpc<20, 5>(kLagsC, kTs, both_weighed(), kNoise, held, load, reading);
        REQUIRE(tg.success);
        REQUIRE(lqg.success);
        REQUIRE(mpc.success);
        static_assert(std::is_same_v<decltype(tg.target_d), Matrix<5, 1>>); // [x; both readings on their way] then u
        CHECK(same_target(tg, lqg.target));
        CHECK(same_target(tg, mpc.target));
        //! Each late reading holds the steady output, C·x_t = x_t.
        CHECK(tg.target_d(2, 0) == doctest::Approx(tg.target_d(0, 0)).epsilon(1e-12));
        CHECK(tg.target_d(3, 0) == doctest::Approx(tg.target_d(1, 0)).epsilon(1e-12));

        //! Both delays, and an input disturbance estimated beside d.
        const design::InputDisturbance spec{{1}};
        const InputDelay<1>            in{0.05};
        const auto                     spec_tg = design::steady_state_target(kLagsC, kTs, spec, load, in, reading);
        const auto                     of = design::state_feedback(kLagsC, kTs, kWeights, kLoadNoise, spec, load, in, reading);
        const auto                     ofm = design::mpc<20, 5>(kLagsC, kTs, both_weighed(), kLoadNoise, spec, load, in, reading);
        REQUIRE(spec_tg.success);
        REQUIRE(of.success);
        REQUIRE(ofm.success);
        static_assert(std::is_same_v<decltype(spec_tg.target_d), Matrix<6, 2>>); // [x; v(k−1); the readings] then u
        CHECK(same_target(spec_tg, of.target));
        CHECK(same_target(spec_tg, ofm.target));
    }

    TEST_CASE("steady_state_target refuses what the designs refuse: an option twice, a TrackedOutputs beside a spec, a bad list, a non-finite B_d") {
        using M = design::MeasuredDisturbance<2, 1>;
        using Held = design::TrackedOutputs<1>;
        using Targets = design::detail::target_options;
        static_assert(Targets::known<M, Held, InputDelay<>, OutputDelay<>>);
        static_assert(!Targets::unique<M, M> && !Targets::unique<Held, Held>);
        static_assert(!design::detail::state_feedback_options::unique<M, M> && !design::detail::mpc_options::unique<M, M> && !design::detail::kalman_options::unique<M, M>);
        //! The held outputs in one place, in every design that holds r: a spec's .outputs, or a TrackedOutputs without one.
        static_assert(!design::detail::held_in_one_place<design::InputDisturbance<1>, Held>);
        static_assert(!design::detail::held_in_one_place<Held, design::DisturbanceModel<2, 0, 1, 1>, M>);
        static_assert(design::detail::held_in_one_place<Held, M, InputDelay<>> && design::detail::held_in_one_place<design::NoDisturbance, Held>);

        //! An output out of range or listed twice, as the LQR refuses it.
        CHECK(design::steady_state_target(*kLags, design::TrackedOutputs{{2}}).status == design::DesignStatus::bad_weights);
        constexpr StateSpace<2, 2, 2> two{.A = kLagsC.A, .B = Matrix<2, 2>::identity(), .C = Matrix<2, 2>::identity()}; // two listed: two inputs
        CHECK(design::steady_state_target(two, kTs, design::TrackedOutputs{{1, 1}}).status == design::DesignStatus::bad_weights);
        CHECK(design::state_feedback(*kLags, kWeights, design::TrackedOutputs{{2}}).status == design::DesignStatus::bad_weights);

        Matrix<2, 1> bad = kLoadOnX2;
        bad(1, 0) = runtime_nan();
        const design::MeasuredDisturbance nan{bad};
        const design::TrackedOutputs      held{{1}};
        CHECK(design::steady_state_target(*kLags, held, nan).status == design::DesignStatus::bad_weights);
        CHECK(design::steady_state_target(kNoisyLags, design::InputDisturbance{{1}}, nan).status == design::DesignStatus::bad_weights);
        CHECK(design::steady_state_target(kLagsC, kTs, held, nan).status == design::DesignStatus::bad_weights);
        CHECK(design::steady_state_target(kLagsC, kTs, held, nan, InputDelay<2>{0.15}, OutputDelay{0.04}).status == design::DesignStatus::bad_weights);
        CHECK(design::state_feedback(*kLags, kWeights, held, nan).status == design::DesignStatus::bad_weights);
        //! A sample time the plant refuses is refused first, as the designs refuse it.
        const design::MeasuredDisturbance load{kLoadOnX2C};
        const auto                        sf = design::state_feedback(kLagsC, 0.0, kWeights, held, load);
        CHECK(sf.status == design::DesignStatus::invalid_sample_time);
        CHECK(design::steady_state_target(kLagsC, 0.0, held, load).status == sf.status);
        CHECK(design::steady_state_target(kLagsC, 0.0, held, nan).status == sf.status);
    }

    TEST_CASE("a target without a MeasuredDisturbance is the one it was: its type and its size") {
        const auto a = design::steady_state_target(kSecond, design::TrackedOutputs{});
        const auto b = design::steady_state_target(*kLags, design::DisturbanceModel{.Bd = kLoadOnX2, .outputs = {1}});
        const auto c = design::steady_state_target(kLagsC, kTs, design::InputDisturbance{{1}}, InputDelay<2>{0.15});
        const auto e = design::steady_state_target(kLagsC, kTs, design::InputDisturbance{{1}}, InputDelay{0.05}, OutputDelay{0.04});
        const auto f = design::steady_state_target(kSecondC, kTs, design::TrackedOutputs{});
        static_assert(std::is_same_v<decltype(a), const design::SteadyStateTargetResult<2, 1, 0, 1, double>>);
        static_assert(std::is_same_v<decltype(b), const design::SteadyStateTargetResult<2, 1, 1, 1, double>>);
        static_assert(std::is_same_v<decltype(c), const design::SteadyStateTargetResult<4, 1, 1, 1, double>>);
        static_assert(std::is_same_v<decltype(e), const design::SteadyStateTargetResult<5, 1, 1, 1, double>>);
        static_assert(std::is_same_v<decltype(f), const design::SteadyStateTargetResult<2, 1, 0, 1, double>>);
        static_assert(std::is_same_v<decltype(SteadyStateTarget{a}), SteadyStateTarget<2, 1, 0, 1, double>>);
        //! The sizes on the base this was built on (x86-64): no column, no member for either.
        if constexpr (sizeof(size_t) == 8) {
            CHECK(sizeof(a) == 40);
            CHECK(sizeof(b) == 56);
            CHECK(sizeof(c) == 88);
            CHECK(sizeof(e) == 104);
            CHECK(sizeof(f) == 40);
            CHECK(sizeof(SteadyStateTarget{a}) == 32);
            CHECK(sizeof(SteadyStateTarget{b}) == 48);
        }
        //! A TrackedOutputs sizes r's map by its list, and a MeasuredDisturbance adds d's NM columns after the spec's.
        static_assert(std::is_same_v<decltype(design::steady_state_target(kLagsC, kTs, design::TrackedOutputs{{1}})), design::SteadyStateTargetResult<2, 1, 0, 1, double>>);
        static_assert(std::is_same_v<decltype(design::steady_state_target(kSecond, design::TrackedOutputs{}, design::MeasuredDisturbance{kLoadOnX2})), design::SteadyStateTargetResult<2, 1, 1, 1, double>>);
        static_assert(std::is_same_v<decltype(design::steady_state_target(kSecond, design::MeasuredDisturbance{kLoadOnX2})), design::SteadyStateTargetResult<2, 1, 1, 0, double>>); // the regulator's: no r
        static_assert(std::is_same_v<decltype(design::steady_state_target(*kLags, design::InputDisturbance{{1}}, design::MeasuredDisturbance{kLoadOnX2})), design::SteadyStateTargetResult<2, 1, 2, 1, double>>);
    }
}
