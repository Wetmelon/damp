// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <bit>
#include <cstddef>
#include <cstdint>

#include "damp/backend.hpp"
#include "damp/controllers/law.hpp"
#include "damp/design/qp.hpp"
#include "damp/design/state_limits.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/estimation/mhe.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

// Deterministic "measurement noise": a fixed repeating pattern so the tests
// are reproducible while still exercising the estimator like real noise does.
constexpr damp::array<double, 12> noise_pattern{0.030, -0.050, 0.010, -0.040, 0.050, -0.020, 0.045, -0.035, 0.005, -0.045, 0.025, -0.015};
constexpr double                  noise(size_t k) {
    return noise_pattern[k % noise_pattern.size()];
}

/// A quiet NaN from its bits, read through a volatile so the compiler cannot fold it (-ffast-math).
template<typename T>
T runtime_nan() {
    if constexpr (sizeof(T) == 4) {
        volatile std::uint32_t bits = 0x7fc00000u;
        return std::bit_cast<float>(static_cast<std::uint32_t>(bits));
    } else {
        volatile std::uint64_t bits = 0x7ff8000000000000ull;
        return std::bit_cast<double>(static_cast<std::uint64_t>(bits));
    }
}

/// +Inf from its bits, read through a volatile.
template<typename T>
T runtime_inf() {
    if constexpr (sizeof(T) == 4) {
        volatile std::uint32_t bits = 0x7f800000u;
        return std::bit_cast<float>(static_cast<std::uint32_t>(bits));
    } else {
        volatile std::uint64_t bits = 0x7ff0000000000000ull;
        return std::bit_cast<double>(static_cast<std::uint64_t>(bits));
    }
}

/// The draining tank of the tests: the level drops 10 % a tick, the inflow adds to it, the sensor reads it.
constexpr StateSpace<1, 1, 1, double, 1, 1> kTank{.A = {{0.90}}, .B = {{0.10}}, .C = {{1.0}}, .Ts = 0.1};

/// A level that is never negative, hard: the estimate never leaves it.
constexpr design::StateLimits kHardLevel{.states = {0}, .lower = {0.0}, .lower_ecr = {0.0}};

/// The same lower end, soft (the default ECR of 1).
constexpr design::StateLimits kSoftLevel{.states = {0}, .lower = {0.0}};

/// A design at namespace scope, for Law<kRes, T>.
inline constexpr auto kLawRes = design::mhe<4>(kTank, {.process = 1e-4, .measurement = 2.5e-3}, kHardLevel);

} // namespace

TEST_CASE("example: tank-level estimator that cannot report a negative level") {
    // ------------------------------------------------------------------
    // The MHE workflow, end to end. Scenario: a draining tank whose level
    // is measured with a noisy sensor. Near empty, sensor noise makes a
    // Kalman filter report physically impossible negative levels; the MHE
    // gets the constraint x >= 0 and cannot.
    // ------------------------------------------------------------------

    // 1. Model the plant, exactly as for a Kalman filter design:
    //    level drains 10%/tick, inflow u adds directly; sensor reads level.
    constexpr StateSpace<1, 1, 1, double, 1, 1> tank{
        .A = {{0.90}},
        .B = {{0.10}},
        .C = {{1.0}},
        .Ts = 0.1,
    };

    // 2. Pick the noise covariances - the SAME numbers a Kalman design uses.
    //    Small Q: the drain model is trusted. R matches the sensor noise
    //    (std ~0.05 -> R = 0.0025).
    const Matrix Q{{1.0e-4}};
    const Matrix R{{2.5e-3}};

    // 3. State the physics the estimate must respect: a level is never
    //    negative, hard (lower_ecr 0). (Left out, every state is unbounded -
    //    this is the only extra input MHE needs over a Kalman filter.)
    const design::StateLimits limits{.states = {0}, .lower = {0.0}, .lower_ecr = {0.0}};

    // 4. Synthesize with an 8-sample window and construct the runtime.
    const auto art = design::mhe<8>(tank, {.process = Q, .measurement = R}, limits);
    REQUIRE(art.success);
    MHE estimator{art};
    estimator.set_initial_state(ColVec{0.5}); // prior for the first window

    // Reference: the equivalent (unconstrained) steady-state Kalman filter.
    StateSpace<1, 1, 1, double, 1, 1> tank_noise{};
    tank_noise.A = tank.A;
    tank_noise.B = tank.B;
    tank_noise.C = tank.C;
    tank_noise.Ts = tank.Ts;
    const auto kf_design = design::kalman(tank_noise, design::KalmanNoise{.process = Q, .measurement = R});
    REQUIRE(kf_design.success);
    KalmanFilter<1, 1, 1, double, 1, 1> kf{kf_design};
    kf.set_state(ColVec{0.5});

    // 5. Run: truth drains from 0.5 toward empty; each tick both estimators
    //    receive the same noisy level reading and the same (zero) inflow.
    double x_true = 0.5;
    double kf_min = 1.0;
    double mhe_min = 1.0;
    double mhe_final = 0.0;
    for (size_t k = 0; k < 80; ++k) {
        const double u = 0.0;
        x_true = (0.90 * x_true) + (0.10 * u);
        const ColVec<1> y{x_true + noise(k)};

        const auto x_mhe = estimator.step(y, ColVec<1>{u});
        CHECK(estimator.last_status() == design::QPStatus::Success);

        kf.predict(ColVec<1>{u});
        REQUIRE(kf.update(y));

        mhe_min = (x_mhe(0) < mhe_min) ? x_mhe(0) : mhe_min;
        kf_min = (kf.state()(0) < kf_min) ? kf.state()(0) : kf_min;
        mhe_final = x_mhe(0);
    }

    // 6. The payoff: the Kalman filter went negative near empty; the MHE
    //    stayed physical, and still converged to the (near-empty) truth.
    CHECK(kf_min < 0.0);
    CHECK(mhe_min >= -1e-9);
    CHECK(mhe_final == doctest::Approx(x_true).epsilon(0.05));
}

TEST_CASE("unconstrained MHE matches the steady-state Kalman filter") {
    // Roadmap #15 acceptance: with no active constraints, the window endpoint
    // reproduces the steady-state Kalman estimate (the arrival cost is the
    // steady-state Kalman covariance, so the two solve the same problem).
    constexpr double Ts = 0.1;

    StateSpace<2, 1, 1, double, 2, 1> plant{};
    plant.A = Matrix<2, 2>{{1.0, Ts}, {0.0, 1.0}};
    plant.B = Matrix<2, 1>{{0.5 * Ts * Ts}, {Ts}};
    plant.C = Matrix{{1.0, 0.0}};
    plant.Ts = Ts;

    const Matrix Q{{1.0e-4, 0.0}, {0.0, 1.0e-3}};
    const Matrix R{{2.5e-3}};

    const auto art = design::mhe<8>(plant, {.process = Q, .measurement = R});
    REQUIRE(art.success);
    MHE estimator{art};

    const auto kf_design = design::kalman(plant, design::KalmanNoise{.process = Q, .measurement = R});
    REQUIRE(kf_design.success);
    KalmanFilter<2, 1, 1, double, 2, 1> kf{kf_design};

    ColVec x_true{0.2, -0.1};
    for (size_t k = 0; k < 60; ++k) {
        const double u = 0.3 * noise(k + 3); // wiggle the input deterministically
        x_true = ColVec<2>{
            x_true(0) + (Ts * x_true(1)) + (0.5 * Ts * Ts * u),
            x_true(1) + (Ts * u),
        };
        const ColVec<1> y{x_true(0) + (0.3 * noise(k))};

        const auto x_mhe = estimator.step(y, ColVec<1>{u});
        kf.predict(ColVec<1>{u});
        REQUIRE(kf.update(y));

        if (k >= 20) { // after the window has flushed its priming samples
            CHECK(x_mhe(0) == doctest::Approx(kf.state()(0)).epsilon(0.02));
            CHECK(x_mhe(1) == doctest::Approx(kf.state()(1)).epsilon(0.05));
        }
    }
}

TEST_CASE("MHE warm start settles to a few iterations per tick") {
    constexpr StateSpace<1, 1, 1, double, 1, 1> tank{
        .A = {{0.90}},
        .B = {{0.10}},
        .C = {{1.0}},
        .Ts = 0.1,
    };
    const auto art = design::mhe<8>(tank, {.process = 1e-4, .measurement = 2.5e-3});
    REQUIRE(art.success);
    MHE estimator{art};

    size_t settled_iterations = 0;
    for (size_t k = 0; k < 40; ++k) {
        (void)estimator.step(ColVec<1>{0.2 + (0.05 * noise(k))}, {});
        if (k >= 30) {
            settled_iterations += estimator.last_iterations();
        }
    }
    // Unconstrained window at steady state: the QP is solved at its
    // unconstrained minimum, no active-set changes at all.
    CHECK(settled_iterations == 0);
}

TEST_CASE("mhe rejects invalid inputs") {
    StateSpace<1, 1, 1, double, 1, 1> continuous{};
    continuous.A = Matrix{{0.9}};
    continuous.B = Matrix{{0.1}};
    continuous.C = Matrix{{1.0}};
    continuous.Ts = 0.0; // not discrete
    const design::KalmanNoise<1, 1> noise{.process = 1e-4, .measurement = 1e-3};
    CHECK(design::mhe<4>(continuous, noise).status == design::DesignStatus::invalid_sample_time);

    StateSpace<1, 1, 1, double, 1, 1> feedthrough = continuous;
    feedthrough.Ts = 0.1;
    feedthrough.D = Matrix{{1.0}};
    CHECK(design::mhe<4>(feedthrough, noise).status == design::DesignStatus::bad_weights);

    StateSpace<1, 1, 1, double, 1, 1> ok = feedthrough;
    ok.D = Matrix<1, 1>{};
    CHECK(design::mhe<4>(ok, {.process = 0.0, .measurement = 1e-3}).status == design::DesignStatus::bad_weights); // zero Q

    StateSpace<1, 1, 1, double, 1, 1> nan_b = ok;
    nan_b.B(0, 0) = runtime_nan<double>();
    CHECK(design::mhe<4>(nan_b, noise).status == design::DesignStatus::bad_weights);

    //! Limits no window can meet: every tick would be infeasible and return a stale estimate.
    CHECK(design::mhe<4>(ok, noise, design::StateLimits{.states = {0}, .lower = {1.0}, .upper = {0.0}}).status == design::DesignStatus::invalid_parameter);
    design::StateLimits<1> nan_limit{.states = {0}};
    nan_limit.upper(0) = runtime_nan<double>();
    CHECK(design::mhe<4>(ok, noise, nan_limit).status == design::DesignStatus::invalid_parameter);
    design::StateLimits<1> inf_limit{.states = {0}};
    inf_limit.lower(0) = -runtime_inf<double>();
    CHECK(design::mhe<4>(ok, noise, inf_limit).status == design::DesignStatus::invalid_parameter);

    //! A state past the plant's, or one listed twice, as design::mpc refuses them.
    CHECK(design::mhe<4>(ok, noise, design::StateLimits{.states = {1}, .lower = {0.0}}).status == design::DesignStatus::invalid_parameter);
    CHECK(design::mhe<4>(ok, noise, design::StateLimits{.states = {0, 0}}).status == design::DesignStatus::invalid_parameter);

    //! A softness below 0, or not a number, as design::mpc's ECR.
    CHECK(design::mhe<4>(ok, noise, design::StateLimits{.states = {0}, .lower = {0.0}, .lower_ecr = {-1.0}}).status == design::DesignStatus::bad_weights);
    design::StateLimits<1> nan_ecr{.states = {0}, .lower = {0.0}};
    nan_ecr.upper_ecr(0) = runtime_nan<double>();
    CHECK(design::mhe<4>(ok, noise, nan_ecr).status == design::DesignStatus::bad_weights);
}

TEST_CASE("MHE takes the noise through the plant's G and H, as design::kalman") {
    //! Noise on the velocity only, G = [0; 1]: Q is 1 × 1. The unconstrained window endpoint is the
    //! Kalman filter of the same (G, Q, R); it was 0.12 off with Q read as the position's.
    constexpr double                  Ts = 0.1;
    StateSpace<2, 1, 1, double, 1, 1> plant{};
    plant.A = Matrix<2, 2>{{1.0, Ts}, {0.0, 1.0}};
    plant.B = Matrix<2, 1>{{0.5 * Ts * Ts}, {Ts}};
    plant.C = Matrix{{1.0, 0.0}};
    plant.G = Matrix{{0.0}, {1.0}};
    plant.Ts = Ts;
    const design::KalmanNoise<1, 1> noise{.process = 1e-3, .measurement = 2.5e-3};

    const auto art = design::mhe<8>(plant, noise);
    REQUIRE(art.success);
    MHE        estimator{art};
    const auto kf_design = design::kalman(plant, noise);
    REQUIRE(kf_design.success);
    KalmanFilter<2, 1, 1, double, 1, 1> kf{kf_design};

    ColVec x_true{0.2, -0.1};
    double max_err = 0.0;
    for (size_t k = 0; k < 120; ++k) {
        const double u = 0.3 * damp::sin(0.7 * static_cast<double>(k));
        x_true = plant.A * x_true + plant.B * ColVec<1>{u};
        const ColVec<1> y{x_true(0) + 0.05 * damp::sin(2.3 * static_cast<double>(k))};
        const auto      x_mhe = estimator.step(y, ColVec<1>{u});
        kf.predict(ColVec<1>{u});
        REQUIRE(kf.update(y));
        if (k >= 80) {
            max_err = damp::max(max_err, (x_mhe - kf.state()).norm());
        }
    }
    CHECK(max_err < 1e-6);
}

TEST_CASE("MHE on a tight QP budget still returns an estimate inside the box") {
    //! The dual active-set iterate at MaxIterations is not primal feasible; it went to −5.7e-3.
    constexpr StateSpace<1, 1, 1, double, 1, 1> tank{
        .A = {{0.90}},
        .B = {{0.10}},
        .C = {{1.0}},
        .Ts = 0.1,
    };
    const auto art = design::mhe<8>(tank, {.process = 1e-4, .measurement = 2.5e-3}, kHardLevel, design::QPIterations{2});
    REQUIRE(art.success);
    CHECK(art.max_qp_iterations == 2);
    MHE estimator{art};

    double lowest = 1.0;
    for (size_t k = 0; k < 60; ++k) {
        const double y = (k < 20) ? 0.3 : -0.05 + (0.1 * damp::sin(2.1 * static_cast<double>(k)));
        (void)estimator.step(ColVec<1>{y}, {});
        for (size_t w = 0; w <= 8; ++w) {
            lowest = damp::min(lowest, estimator.window_state(w)(0));
        }
    }
    CHECK(lowest >= 0.0);
}

TEST_CASE("MHE takes its QP budget as QPIterations, in any order with the limits: left out the default, 0 refused") {
    constexpr StateSpace<1, 1, 1, double, 1, 1> tank{.A = {{0.90}}, .B = {{0.10}}, .C = {{1.0}}, .Ts = 0.1};
    const auto                                  plain = design::mhe<8>(tank, {.process = 1e-4, .measurement = 2.5e-3});
    REQUIRE(plain.success);
    CHECK(plain.max_qp_iterations == 10 * 9); // 10·(NZ + NI): 9 window states, no rows
    const auto bounded = design::mhe<8>(tank, {.process = 1e-4, .measurement = 2.5e-3}, kHardLevel);
    REQUIRE(bounded.success);
    CHECK(bounded.max_qp_iterations == 10 * (10 + 19)); // the 9 window states and ε; 18 box rows and ε ≥ 0

    //! The same window QP, its budget the option's.
    const auto budget = design::mhe<8>(tank, {.process = 1e-4, .measurement = 2.5e-3}, design::QPIterations{20});
    REQUIRE(budget.success);
    CHECK(budget.max_qp_iterations == 20);
    CHECK(budget.H == plain.H);
    static_assert(design::mhe<3>(tank, {.process = 1e-4, .measurement = 2.5e-3}, design::QPIterations{7}).max_qp_iterations == 7);
    const auto either = design::mhe<8>(tank, {.process = 1e-4, .measurement = 2.5e-3}, design::QPIterations{20}, kHardLevel);
    REQUIRE(either.success);
    CHECK(either.max_qp_iterations == 20);
    CHECK(either.H == bounded.H);
    CHECK(either.b_bound == bounded.b_bound);

    //! Not one active-set step, as design::mpc refuses it.
    const auto none = design::mhe<8>(tank, {.process = 1e-4, .measurement = 2.5e-3}, design::QPIterations{0});
    CHECK_FALSE(none.success);
    CHECK(none.status == design::DesignStatus::invalid_parameter);
}

TEST_CASE("MHE with a state pinned at zero by equal limits solves every tick") {
    //! Rows x ≤ 0 and −x ≤ 0: once one is active the other reads x ≈ ±1e-20 against a roundoff
    //! scale of |x| itself, counted as violated, and the solve went Infeasible on 38 of 40 ticks.
    constexpr double                  Ts = 0.1;
    StateSpace<2, 1, 1, double, 2, 1> plant{};
    plant.A = Matrix<2, 2>{{1.0, Ts}, {0.0, 1.0}};
    plant.B = Matrix<2, 1>{{0.5 * Ts * Ts}, {Ts}};
    plant.C = Matrix{{1.0, 0.0}};
    plant.Ts = Ts;
    const design::StateLimits limits{.states = {1}, .lower = {0.0}, .upper = {0.0}, .lower_ecr = {0.0}, .upper_ecr = {0.0}};
    const auto                art = design::mhe<4>(plant, {.process = Matrix{{1e-4, 0.0}, {0.0, 1e-3}}, .measurement = 2.5e-3}, limits);
    REQUIRE(art.success);
    MHE estimator{art};

    size_t failed = 0;
    for (size_t k = 0; k < 40; ++k) {
        const auto x_hat = estimator.step(ColVec<1>{0.5 + (0.3 * noise(k))}, {});
        if (estimator.last_status() != design::QPStatus::Success) {
            ++failed;
        }
        CHECK(damp::abs(x_hat(1)) < 1e-12);
    }
    CHECK(failed == 0);
}

TEST_CASE("mhe is constexpr") {
    constexpr StateSpace<1, 1, 1, double, 1, 1> tank{
        .A = {{0.90}},
        .B = {{0.10}},
        .C = {{1.0}},
        .Ts = 0.1,
    };
    constexpr auto art = design::mhe<3>(tank, {.process = 1e-4, .measurement = 2.5e-3});
    static_assert(art.success);
    const MHE<1, 1, 1, 3, double> estimator = art; // assign-from-design
    CHECK(estimator.valid());
}

TEST_CASE("MHE deploys in float") {
    constexpr StateSpace<1, 1, 1, double, 1, 1> tank{
        .A = {{0.90}},
        .B = {{0.10}},
        .C = {{1.0}},
        .Ts = 0.1,
    };
    const auto art = design::mhe<6>(tank, {.process = 1e-4, .measurement = 2.5e-3}, kHardLevel);
    REQUIRE(art.success);

    MHE   estimator{art.as<float>()};
    float x_true = 0.4F;
    float x_hat = 0.0F;
    for (size_t k = 0; k < 60; ++k) {
        x_true *= 0.90F;
        const auto e = estimator.step(ColVec<1, float>{x_true + (0.05F * static_cast<float>(noise(k)))}, {});
        x_hat = e(0);
        CHECK(x_hat >= -1e-6F);
    }
    CHECK(x_hat == doctest::Approx(x_true).epsilon(0.1));
}

TEST_CASE("MHE refuses a float deployment whose window Hessian does not factor in float") {
    //! Noise on the velocity only (G = [0; 1]) weights the noise-free direction by 1/√ε of the
    //! design scalar: diag(H) spans 2e3 to 2.7e11. as<float>() reported success and every float
    //! tick failed with NotPositiveDefinite, holding x̂ = 0.
    constexpr double                  Ts = 0.1;
    StateSpace<2, 1, 1, double, 1, 1> plant{};
    plant.A = Matrix<2, 2>{{1.0, Ts}, {0.0, 1.0}};
    plant.B = Matrix<2, 1>{{0.5 * Ts * Ts}, {Ts}};
    plant.C = Matrix{{1.0, 0.0}};
    plant.G = Matrix{{0.0}, {1.0}};
    plant.Ts = Ts;
    const auto art = design::mhe<8>(plant, {.process = 1e-3, .measurement = 2.5e-3});
    REQUIRE(art.success);

    const auto deployed = art.as<float>();
    CHECK_FALSE(deployed.success);
    CHECK(deployed.status == design::DesignStatus::singular);
    const MHE estimator{deployed};
    CHECK_FALSE(estimator.valid());
}

TEST_CASE("Unconstrained MHE endpoint equals the steady-state Kalman filter") {
    // Arrival cost weight is P∞⁻¹ (a-priori). Its mean must be the matching
    // a-priori prediction, not the previous smoothed estimate (which already
    // used the measurements still inside the window: double counting).
    constexpr double                  Ts = 0.1;
    StateSpace<2, 1, 1, double, 2, 1> plant{};
    plant.A = Matrix<2, 2>{{1.0, Ts}, {0.0, 1.0}};
    plant.B = Matrix<2, 1>{{0.5 * Ts * Ts}, {Ts}};
    plant.C = Matrix{{1.0, 0.0}};
    plant.Ts = Ts;
    const Matrix Q{{1.0e-4, 0.0}, {0.0, 1.0e-3}};
    const Matrix R{{2.5e-3}};

    const auto art = design::mhe<8>(plant, {.process = Q, .measurement = R});
    REQUIRE(art.success);
    MHE        estimator{art};
    const auto kf_design = design::kalman(plant, design::KalmanNoise{.process = Q, .measurement = R});
    REQUIRE(kf_design.success);
    KalmanFilter<2, 1, 1, double, 2, 1> kf{kf_design};

    ColVec x_true{0.2, -0.1};
    double max_err = 0.0;
    for (size_t k = 0; k < 120; ++k) {
        const double u = 0.3 * damp::sin(0.7 * static_cast<double>(k));
        x_true = plant.A * x_true + plant.B * ColVec<1>{u};
        const ColVec<1> y{x_true(0) + 0.05 * damp::sin(2.3 * static_cast<double>(k))};
        const auto      x_mhe = estimator.step(y, ColVec<1>{u});
        kf.predict(ColVec<1>{u});
        REQUIRE(kf.update(y));
        if (k >= 80) { // priming transient decays geometrically
            max_err = damp::max(max_err, (x_mhe - kf.state()).norm());
        }
    }
    CHECK(max_err < 1e-6);
}

TEST_CASE("mhe reads a plant with no noise inputs as noise on every state and reading") {
    constexpr double                            Ts = 0.1;
    constexpr StateSpace                        bare{.A = {{1.0, Ts}, {0.0, 1.0}}, .B = {{0.5 * Ts * Ts}, {Ts}}, .C = {{1.0, 0.0}}, .Ts = Ts};
    constexpr StateSpace<2, 1, 1, double, 2, 1> full{.A = bare.A, .B = bare.B, .C = bare.C, .Ts = Ts};
    const Matrix                                Q{{1.0e-4, 0.0}, {0.0, 1.0e-3}};
    const Matrix                                R{{2.5e-3}};

    const auto by_bare = design::mhe<4>(bare, {.process = Q, .measurement = R});
    const auto by_full = design::mhe<4>(full, {.process = Q, .measurement = R});
    REQUIRE(by_bare.success);
    CHECK((by_bare.H - by_full.H).norm() <= 1e-12 * by_full.H.norm());
}

TEST_CASE("A StateLimits bounds only the states it lists, and the QP carries rows for those alone") {
    constexpr double                            Ts = 0.1;
    constexpr double                            big = damp::numeric_limits<double>::max();
    constexpr StateSpace<2, 1, 1, double, 2, 1> plant{.A = {{1.0, Ts}, {0.0, 1.0}}, .B = {{0.5 * Ts * Ts}, {Ts}}, .C = {{1.0, 0.0}}, .Ts = Ts};
    const design::KalmanNoise<2, 1>             cov{.process = Matrix{{1e-4, 0.0}, {0.0, 1e-3}}, .measurement = 2.5e-3};

    //! The speed hard at ±0.05, the position free; then both listed, the position unbounded: the
    //! rows a box on every state carried.
    const design::StateLimits speed{.states = {1}, .lower = {-0.05}, .upper = {0.05}, .lower_ecr = {0.0}, .upper_ecr = {0.0}};
    const design::StateLimits both{.states = {0, 1}, .lower = {-big, -0.05}, .upper = {big, 0.05}, .lower_ecr = {0.0, 0.0}, .upper_ecr = {0.0, 0.0}};
    const auto                one = design::mhe<6>(plant, cov, speed);
    const auto                every = design::mhe<6>(plant, cov, both);
    const auto                free = design::mhe<6>(plant, cov);
    REQUIRE(one.success);
    REQUIRE(every.success);
    REQUIRE(free.success);

    using One = decltype(one);
    using Every = decltype(every);
    using Free = decltype(free);
    static_assert(One::NV == 14 && One::NZ == 15 && One::NIB == 14 && One::NI == 15);         // 7 window states of 2; ε; the speed's 2 rows at each, then ε ≥ 0
    static_assert(Every::NV == 14 && Every::NZ == 15 && Every::NIB == 28 && Every::NI == 29); // twice the box rows
    static_assert(Free::NZ == 14 && Free::NI == 0);                                           // no rows, no slack

    MHE          est_one{one};
    MHE          est_every{every};
    MHE          est_free{free};
    ColVec<2>    x_true{0.0, 0.1}; // twice the speed limit
    double       fastest = 0.0;
    double       worst_gap = 0.0;
    const double u = 0.0;
    for (size_t k = 0; k < 60; ++k) {
        x_true = plant.A * x_true + plant.B * ColVec<1>{u};
        const ColVec<1> y{x_true(0) + (0.3 * noise(k))};
        const auto      a = est_one.step(y, ColVec<1>{u});
        const auto      b = est_every.step(y, ColVec<1>{u});
        (void)est_free.step(y, ColVec<1>{u});
        CHECK(est_one.last_status() == design::QPStatus::Success);
        for (size_t w = 0; w <= 6; ++w) {
            fastest = damp::max(fastest, damp::abs(est_one.window_state(w)(1)));
        }
        worst_gap = damp::max(worst_gap, (a - b).norm());
    }
    CHECK(fastest <= 0.05);                             // the listed state held, at every window state
    CHECK(est_free.state()(1) > 0.08);                  // where the readings put it without the limit
    CHECK(est_one.state()(0) > 0.9 * x_true(0));        // the position follows the readings (0.6): no limit on it
    CHECK(worst_gap <= 1e-12 * est_one.state().norm()); // an unbounded state's rows change nothing
}

TEST_CASE("A soft StateLimits gives way to a strong pull of the readings, a hard one does not") {
    //! The level reads 0.2 below empty, tick after tick: a sensor offset the limit disagrees with.
    const design::KalmanNoise<1, 1> noise{.process = 1e-4, .measurement = 1e-3};
    const auto                      hard = design::mhe<8>(kTank, noise, kHardLevel);
    const auto                      soft = design::mhe<8>(kTank, noise, kSoftLevel);
    const auto                      softer = design::mhe<8>(kTank, noise, design::StateLimits{.states = {0}, .lower = {0.0}, .lower_ecr = {10.0}});
    const auto                      free = design::mhe<8>(kTank, noise);
    REQUIRE(hard.success);
    REQUIRE(soft.success);
    REQUIRE(softer.success);
    REQUIRE(free.success);
    static_assert(decltype(hard)::NI == decltype(soft)::NI); // the same rows: the softness is a coefficient on ε

    MHE    h{hard};
    MHE    s{soft};
    MHE    s10{softer};
    MHE    f{free};
    double lowest_hard = 1.0;
    for (size_t k = 0; k < 40; ++k) {
        const ColVec<1> y{-0.2};
        (void)h.step(y, {});
        (void)s.step(y, {});
        (void)s10.step(y, {});
        (void)f.step(y, {});
        for (size_t w = 0; w <= 8; ++w) {
            lowest_hard = damp::min(lowest_hard, h.window_state(w)(0));
        }
    }
    CHECK(lowest_hard >= 0.0);                // hard: never below, at any window state
    CHECK(h.last_slack() == 0.0);             // and no slack taken
    CHECK(s.state()(0) < -1e-3);              // soft: exceeded
    CHECK(s.last_slack() > 0.0);              // through the slack
    CHECK(s.state()(0) > f.state()(0) + 0.1); // still held well above where the readings pull
    CHECK(s10.state()(0) < s.state()(0));     // a larger ECR gives way more readily
    CHECK(s10.state()(0) > f.state()(0));
}

TEST_CASE("mhe(plant, Ts, ...) is mhe(discretize(plant, Ts), ...) with the process-noise density discretized over the period") {
    constexpr double          Ts = 0.1;
    constexpr StateSpace      plant{.A = {{0.0, 1.0}, {0.0, -0.5}}, .B = {{0.0}, {1.0}}, .C = {{1.0, 0.0}}}; // continuous
    const Matrix              Qc{{1e-4, 0.0}, {0.0, 1e-2}};                                                  // a density
    constexpr double          R = 2.5e-3;                                                                    // the variance of one reading, either form
    const design::StateLimits speed{.states = {1}, .lower = {-0.3}, .upper = {0.3}};

    const auto sampled = design::mhe<6>(plant, Ts, {.process = Qc, .measurement = R}, speed);
    REQUIRE(sampled.success);
    const auto d = discretize(plant, Ts);
    REQUIRE(d);
    const auto Qd = design::discretize_process_noise(plant.A, Matrix<2, 2>::identity(), Qc, Ts);
    REQUIRE(Qd);
    const auto discrete = design::mhe<6>(*d, {.process = *Qd, .measurement = R}, speed);
    REQUIRE(discrete.success);

    //! A few ulps of each matrix: the same window on the same model and noise.
    const auto close = [](const auto& a, const auto& b) { return (a - b).norm() <= 16.0 * eps<double>() * b.norm(); };
    CHECK(close(sampled.H, discrete.H));
    CHECK(close(sampled.M_a, discrete.M_a));
    CHECK(close(sampled.M_cy, discrete.M_cy));
    CHECK(close(sampled.M_qu, discrete.M_qu));
    CHECK(close(sampled.M_aqu, discrete.M_aqu));
    CHECK(close(sampled.Ap_x, discrete.Ap_x));
    CHECK(close(sampled.Ap_y, discrete.Ap_y));
    CHECK(close(sampled.Ap_u, discrete.Ap_u));
    CHECK(sampled.A_con == discrete.A_con);
    CHECK(sampled.b_bound == discrete.b_bound);

    //! Unbounded, its window endpoint is design::kalman(plant, Ts, …)'s filter on the same density.
    const auto unbounded = design::mhe<8>(plant, Ts, {.process = Qc, .measurement = R});
    REQUIRE(unbounded.success);
    const auto kf_design = design::kalman(plant, Ts, design::KalmanNoise{.process = Qc, .measurement = R});
    REQUIRE(kf_design.success);
    MHE                                 estimator{unbounded};
    KalmanFilter<2, 1, 1, double, 2, 1> kf{kf_design};
    ColVec                              x_true{0.2, -0.1};
    double                              max_err = 0.0;
    for (size_t k = 0; k < 120; ++k) {
        const double u = 0.3 * damp::sin(0.7 * static_cast<double>(k));
        x_true = d->A * x_true + d->B * ColVec<1>{u};
        const ColVec<1> y{x_true(0) + 0.05 * damp::sin(2.3 * static_cast<double>(k))};
        const auto      x_mhe = estimator.step(y, ColVec<1>{u});
        kf.predict(ColVec<1>{u});
        REQUIRE(kf.update(y));
        if (k >= 80) { // priming transient decays geometrically
            max_err = damp::max(max_err, (x_mhe - kf.state()).norm());
        }
    }
    CHECK(max_err < 1e-6);

    //! One time domain per form: a discrete plant takes no Ts, a continuous one a finite positive one.
    CHECK(design::mhe<4>(kTank, 0.1, {.process = 1e-4, .measurement = 1e-3}).status == design::DesignStatus::invalid_sample_time);
    CHECK(design::mhe<4>(plant, {.process = Qc, .measurement = R}).status == design::DesignStatus::invalid_sample_time);
    CHECK(design::mhe<4>(plant, 0.0, {.process = Qc, .measurement = R}).status == design::DesignStatus::invalid_sample_time);
    CHECK(design::mhe<4>(plant, runtime_nan<double>(), {.process = Qc, .measurement = R}).status == design::DesignStatus::invalid_sample_time);
    Matrix<2, 2> Qc_nan = Qc;
    Qc_nan(1, 1) = runtime_nan<double>();
    CHECK(design::mhe<4>(plant, Ts, {.process = Qc_nan, .measurement = R}).status == design::DesignStatus::bad_weights);
}

TEST_CASE_TEMPLATE("MHE set_limits moves a bound at run time and refuses one that is not a box", T, double, float) {
    const auto art = design::mhe<6>(kTank, {.process = 1e-4, .measurement = 2.5e-3}, kHardLevel).template as<T>();
    REQUIRE(art.success);
    MHE est{art};
    CHECK(est.limits().state.lower(0) == T{0});
    CHECK(est.limits().state.upper(0) == damp::numeric_limits<T>::max());
    const ColVec<1, T> y{static_cast<T>(-0.1)};
    for (size_t k = 0; k < 30; ++k) {
        (void)est.step(y, {});
    }
    CHECK(est.state()(0) >= T{0});
    CHECK(est.state()(0) < static_cast<T>(1e-6));

    //! Not a box: refused, and nothing changes.
    const auto before = est.limits();
    auto       nan_end = before;
    nan_end.state.lower(0) = runtime_nan<T>();
    CHECK_FALSE(est.set_limits(nan_end));
    auto inf_end = before;
    inf_end.state.upper(0) = runtime_inf<T>();
    CHECK_FALSE(est.set_limits(inf_end));
    auto inverted = before;
    inverted.state.lower(0) = T{1};
    inverted.state.upper(0) = static_cast<T>(0.5);
    CHECK_FALSE(est.set_limits(inverted));
    CHECK(est.limits().state.lower(0) == before.state.lower(0));
    CHECK(est.limits().state.upper(0) == before.state.upper(0));
    (void)est.step(y, {});
    CHECK(est.state()(0) >= T{0});

    //! The floor moved to −0.03, above where the readings pull the level (−0.06): it holds there.
    auto lower = before;
    lower.state.lower(0) = static_cast<T>(-0.03);
    REQUIRE(est.set_limits(lower));
    CHECK(est.limits().state.lower(0) == static_cast<T>(-0.03));
    for (size_t k = 0; k < 30; ++k) {
        (void)est.step(y, {});
    }
    CHECK(est.state()(0) >= static_cast<T>(-0.03));
    CHECK(est.state()(0) < static_cast<T>(-0.029));

    //! Unbounded: the estimate goes where the readings put it.
    auto open = before;
    open.state.lower(0) = -damp::numeric_limits<T>::max();
    REQUIRE(est.set_limits(open));
    for (size_t k = 0; k < 40; ++k) {
        (void)est.step(y, {});
    }
    CHECK(est.state()(0) < static_cast<T>(-0.04));

    //! A design without a StateLimits has no state limits to move.
    MHE free{design::mhe<6>(kTank, {.process = 1e-4, .measurement = 2.5e-3}).template as<T>()};
    static_assert(damp::is_same_v<typename decltype(free)::Limits, design::MHELimits<0, T>>);
    CHECK(free.set_limits({}));
}

TEST_CASE("Law<kRes, T> of a constexpr design::mhe is its MHE, sized by the design") {
    using Res = damp::remove_cvref_t<decltype(kLawRes)>;
    static_assert(Res::nx == 1 && Res::nu == 1 && Res::ny == 1 && Res::nh == 4 && Res::nxc == 1);
    static_assert(Res::NV == 5 && Res::NZ == 6 && Res::NIB == 10 && Res::NI == 11);
    static_assert(damp::is_same_v<Law<kLawRes, float>::Runtime, MHE<1, 1, 1, 4, float, 1>>);
    static_assert(damp::is_same_v<Law<kLawRes, double>::Runtime, MHE<1, 1, 1, 4, double, 1>>);

    static constinit MHE constant{kLawRes}; // a constexpr design runs constinit, sized by deduction
    Law<kLawRes, float>  law;
    MHE                  plain{kLawRes.as<float>()};
    CHECK(law.valid());
    CHECK(constant.valid());
    for (size_t k = 0; k < 30; ++k) {
        const float y = -0.05F + (0.1F * static_cast<float>(damp::sin(0.9 * static_cast<double>(k))));
        const auto  x_law = law.step(ColVec<1, float>{y}, {});
        const auto  x_plain = plain.step(ColVec<1, float>{y}, {});
        (void)constant.step(ColVec<1>{static_cast<double>(y)}, {});
        CHECK(x_law(0) == doctest::Approx(x_plain(0)).epsilon(1e-5)); // the same runtime on the same design (to rounding: two inlined copies)
        CHECK(x_law(0) >= 0.0F);
    }
    CHECK(law.limits().state.lower(0) == 0.0F);
}
