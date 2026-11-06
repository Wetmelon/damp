// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file test_offset_free_mpc.cpp
 * @brief Offset-free MPC: design::ofmpc with KalmanNoise and a disturbance spec, run as OFMPC from the
 *        output measurement alone.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <type_traits>

#include "damp/controllers/mpc.hpp"
#include "damp/controllers/ofmpc.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/qp.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/estimation/luenberger.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"
#include "register_form.hpp"

using namespace damp;

namespace {

constexpr double Ts = 0.1;

// Double integrator, ZOH-discretized at Ts = 0.1, position output only; process noise on both
// states (G = I), measurement noise on the reading (H = I).
constexpr StateSpace<2, 1, 1, double, 2, 1> double_integrator{
    .A = {{1.0, Ts}, {0.0, 1.0}},
    .B = {{0.5 * Ts * Ts}, {Ts}},
    .C = {{1.0, 0.0}},
    .Ts = Ts,
};

/// The estimator's noise: unit process, disturbance and measurement noise.
constexpr design::KalmanNoise<2, 1, double, 1> kUnitNoise{};

constexpr ColVec<2> plant_step(const ColVec<2>& x, double u_effective) {
    return ColVec<2>{
        x(0) + (Ts * x(1)) + (0.5 * Ts * Ts * u_effective),
        x(1) + (Ts * u_effective),
    };
}

/// An offset-free design on the double integrator: an input disturbance, the position held.
template<size_t NP, size_t NC, size_t NT_max = 0>
constexpr auto offset_free(const design::MPCWeights<1, 1>& w, const design::MPCConstraints<1, 1>& c = {}, const design::KalmanNoise<2, 1, double, 1>& noise = kUnitNoise) {
    return design::ofmpc<NP, NC, NT_max>(double_integrator, w, c, noise, design::InputDisturbance{});
}

design::MPCWeights<1, 1> move_weight(double rdu) {
    design::MPCWeights<1, 1> w{};
    w.Rdu = rdu;
    return w;
}

} // namespace

TEST_CASE("offset-free MPC rejects an unknown input disturbance from output feedback") {
    // The controller sees only y = position; the plant carries a constant
    // unmodeled load d. The disturbance estimate must converge to d and the
    // output must settle on the reference with no offset — the case the bare
    // velocity-form MPC cannot handle.
    const auto art = offset_free<15, 5>(move_weight(0.05), {}, design::KalmanNoise{Matrix<2, 2>::identity(), Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{1.0}}});
    REQUIRE(art.success);

    OFMPC           controller{art};
    ColVec<2>       x{};
    const ColVec<1> r{0.5};
    const double    d_load = 0.3;

    for (size_t k = 0; k < 600; ++k) {
        const auto u = controller.step(r, ColVec<1>{x(0)});
        x = plant_step(x, u(0) + d_load);
    }

    CHECK(x(0) == doctest::Approx(0.5).epsilon(1e-3));
    CHECK(x(1) == doctest::Approx(0.0).epsilon(1e-3));
    CHECK(controller.disturbance_estimate()(0) == doctest::Approx(d_load).epsilon(1e-3));
    // Steady state: the applied input cancels the load.
    CHECK(controller.applied()(0) == doctest::Approx(-d_load).epsilon(1e-2));
}

TEST_CASE("offset-free MPC removes steady-state error under plant gain mismatch") {
    // The real plant has 30% more input gain than the model. The disturbance
    // integrator absorbs the mismatch at steady state.
    const auto art = offset_free<15, 5>(move_weight(0.05));
    REQUIRE(art.success);

    OFMPC           controller{art};
    ColVec<2>       x{};
    const ColVec<1> r{1.0};

    for (size_t k = 0; k < 800; ++k) {
        const auto u = controller.step(r, ColVec<1>{x(0)});
        x = plant_step(x, 1.3 * u(0)); // plant gain != model gain
    }

    CHECK(x(0) == doctest::Approx(1.0).epsilon(5e-3));
    CHECK(x(1) == doctest::Approx(0.0).epsilon(5e-3));
}

TEST_CASE("offset-free MPC removes a sensor bias (OutputDisturbance, C_d = I)") {
    // The bias reaches the controlled output, not the state: the target and the QP's reference
    // both hold C·x + C_d·d̂ at r. Without the C_d term the reading settled at r + d.
    constexpr StateSpace<2, 1, 1, double, 2, 1> plant{.A = {{0.9, 0.1}, {0.0, 0.8}}, .B = {{0.0}, {0.2}}, .C = {{1.0, 0.0}}, .Ts = 0.01};
    const design::KalmanNoise<2, 1, double, 1>  noise{.process = Matrix<2, 2>::identity() * 0.01, .measurement = 0.01, .disturbance = 0.01};
    const auto                                  art = design::ofmpc<15, 5>(plant, move_weight(0.05), {}, noise, design::OutputDisturbance{});
    REQUIRE(art.success);

    OFMPC        controller{art};
    ColVec<2>    x{};
    const double bias = 0.3;
    double       y = 0.0;
    for (int k = 0; k < 6000; ++k) {
        y = (plant.C * x)(0) + bias;
        const auto u = controller.step(ColVec<1>{1.0}, ColVec<1>{y});
        x = (plant.A * x) + (plant.B * u);
    }
    CHECK(y == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(controller.disturbance_estimate()(0) == doctest::Approx(bias).epsilon(1e-4));

    // A finite y_lim scales the outputs by its span, C_d·d̂ included: unscaled, the limit row
    // counted the bias 4 times, so r = 1.5 settled at 1.09 and r = 3 under a −1 bias at 3.0.
    design::MPCConstraints<1, 1> limits{};
    limits.y_lim = {.lower = -2.0, .upper = 2.0};
    const auto limited = design::ofmpc<15, 5>(plant, move_weight(0.05), limits, noise, design::OutputDisturbance{});
    REQUIRE(limited.success);
    const auto settle = [&](double b, double r) {
        OFMPC     ctl{limited};
        ColVec<2> xs{};
        double    ys = 0.0;
        for (int k = 0; k < 6000; ++k) {
            ys = (plant.C * xs)(0) + b;
            const auto u = ctl.step(ColVec<1>{r}, ColVec<1>{ys});
            xs = (plant.A * xs) + (plant.B * u);
        }
        return ys;
    };
    CHECK(settle(0.4, 1.5) == doctest::Approx(1.5).epsilon(1e-4));
    CHECK(settle(-1.0, 3.0) == doctest::Approx(2.0).epsilon(1e-3));
}

TEST_CASE("offset-free MPC respects input limits while rejecting a load") {
    design::MPCConstraints<1, 1> limits{};
    limits.u_lim = {
        .lower = -0.6,
        .upper = 0.6,
    };

    const auto art = offset_free<15, 5>(move_weight(0.05), limits);
    REQUIRE(art.success);

    OFMPC           controller{art};
    ColVec<2>       x{};
    const ColVec<1> r{0.5};

    for (size_t k = 0; k < 800; ++k) {
        const auto u = controller.step(r, ColVec<1>{x(0)});
        CHECK(u(0) <= 0.6 + 1e-9);
        CHECK(u(0) >= -0.6 - 1e-9);
        x = plant_step(x, u(0) + 0.3); // load leaves ±0.3 of authority
    }
    CHECK(x(0) == doctest::Approx(0.5).epsilon(5e-3));
}

TEST_CASE("caller-sequenced predict/update/control matches the fused tick") {
    const auto art = offset_free<10, 4>(move_weight(0.05));
    REQUIRE(art.success);

    OFMPC fused{art};
    OFMPC sequenced{art};

    ColVec<2>       xf{0.2, -0.1};
    ColVec<2>       xs{0.2, -0.1};
    const ColVec<1> r{0.4};

    for (size_t k = 0; k < 50; ++k) {
        const auto uf = fused.step(r, ColVec<1>{xf(0)});

        sequenced.back_calculate(sequenced.last_command());
        sequenced.predict();
        sequenced.update(ColVec<1>{xs(0)});
        const auto us = sequenced.control(r);

        CHECK(uf(0) == doctest::Approx(us(0)).epsilon(1e-12));
        xf = plant_step(xf, uf(0));
        xs = plant_step(xs, us(0));
    }
}

TEST_CASE("ofmpc rejects undetectable disturbance models") {
    // B = 0: the input disturbance never reaches the output — the
    // Pannocchia–Rawlings rank condition fails.
    StateSpace<1, 1, 1, double, 1, 1> plant{};
    plant.A = Matrix<1, 1>{{0.9}};
    plant.B = Matrix<1, 1>{{0.0}};
    plant.C = Matrix<1, 1>{{1.0}};
    plant.Ts = 0.1;

    const auto art = design::ofmpc<5, 5>(plant, {}, {}, design::KalmanNoise<1, 1, double, 1>{}, design::InputDisturbance{});
    CHECK_FALSE(art.success);
    CHECK(art.status == design::DesignStatus::not_detectable);
}

TEST_CASE("offset-free mpc is constexpr") {
    constexpr auto art = offset_free<3, 3>({});
    static_assert(art.success);
    static_assert(art.status == design::DesignStatus::ok);
    static_assert(art.kalman.sys.Ts == Ts); // the plant's, not an invented sample time
}

TEST_CASE("a continuous plant is refused by the discrete form; continuous plant + Ts designs at Ts") {
    auto continuous = double_integrator;
    continuous.Ts = 0.0;
    const auto refused = design::ofmpc<15, 5>(continuous, {}, {}, kUnitNoise, design::InputDisturbance{});
    CHECK_FALSE(refused.success);
    CHECK(refused.status == design::DesignStatus::invalid_sample_time);
    const auto refused_model = design::ofmpc<15, 5>(continuous, {}, {}, kUnitNoise, design::DisturbanceModel<2, 1, 1>{.Bd = continuous.B});
    CHECK_FALSE(refused_model.success);
    CHECK(refused_model.status == design::DesignStatus::invalid_sample_time);

    //! The sample-time form takes the continuous plant only: a discrete one would be discretized twice.
    constexpr StateSpace<2, 1, 1, double, 2, 1> di{.A = Matrix<2, 2>{{0.0, 1.0}, {0.0, 0.0}}, .B = Matrix<2, 1>{{0.0}, {1.0}}, .C = Matrix<1, 2>{{1.0, 0.0}}};
    const auto                                  sampled = design::ofmpc<10, 4>(di, 0.1, {}, {}, kUnitNoise, design::InputDisturbance{});
    CHECK(sampled.success);
    CHECK(sampled.status == design::DesignStatus::ok);
    CHECK(sampled.kalman.sys.Ts == 0.1);
    StateSpace<2, 1, 1, double, 2, 1> discrete{};
    discrete.A = discretize(di, 0.1)->A;
    discrete.B = discretize(di, 0.1)->B;
    discrete.C = di.C;
    discrete.Ts = 0.1;
    CHECK(design::ofmpc<10, 4>(discrete, 0.1, {}, {}, kUnitNoise, design::InputDisturbance{}).status == design::DesignStatus::invalid_sample_time);
    CHECK(design::ofmpc<10, 4>(di, 0.0, {}, {}, kUnitNoise, design::InputDisturbance{}).status == design::DesignStatus::invalid_sample_time);
}

TEST_CASE("offset-free MPC deploys in float") {
    const auto art = offset_free<10, 4>(move_weight(0.05));
    REQUIRE(art.success);

    OFMPC                  controller{art.as<float>()};
    ColVec<2, float>       x{};
    const ColVec<1, float> r{0.5F};

    for (size_t k = 0; k < 600; ++k) {
        const auto u = controller.step(r, ColVec<1, float>{x(0)});
        x = ColVec<2, float>{
            x(0) + (0.1F * x(1)) + (0.005F * (u(0) + 0.3F)),
            x(1) + (0.1F * (u(0) + 0.3F)),
        };
    }
    CHECK(x(0) == doctest::Approx(0.5).epsilon(0.01));
}

TEST_CASE("step with the applied input: a downstream clamp reaches the held input and the estimator") {
    const auto art = offset_free<15, 5>(move_weight(0.05), {}, design::KalmanNoise{Matrix<2, 2>::identity(), Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{1.0}}});
    REQUIRE(art.success);

    OFMPC           fused{art};
    OFMPC           split{art};
    const ColVec<1> r{0.5};
    ColVec<2>       x{};
    for (size_t k = 0; k < 5; ++k) {
        const ColVec<1> y{x(0)};
        const auto      u_fused = fused.step(r, y);
        split.back_calculate(split.last_command());
        const auto u_req = split.step(r, y);
        CHECK(u_req(0) == doctest::Approx(u_fused(0)));
        CHECK(split.applied()(0) != doctest::Approx(u_req(0))); // not held until the next predict
        x = plant_step(x, u_fused(0));
    }

    // A tighter downstream clamp: the next step predicts with the applied input, and its move starts from it.
    const ColVec<1> y{x(0)};
    split.back_calculate(ColVec<1>{0.01});
    (void)split.step(ColVec<1>{5.0}, y);
    CHECK(split.applied()(0) == 0.01);
}

namespace {

/// Motor and load masses on a compliant transmission: x = [p_motor, v_motor, p_load, v_load].
/// The motor force is the input; a force on the load is the disturbance. Sensors: the motor
/// encoder and an accelerometer on the load, which feels the load force directly. Nothing
/// measures the load position, the one output that matters: it is a third output row whose
/// reading the filter all but ignores (a variance of 10⁴ m²), listed as the controlled output.
struct CompliantDrive {
    static constexpr double Ts = 0.01;
    static constexpr double k = 400.0; ///< [N/m] transmission stiffness
    static constexpr double m = 1.0;   ///< [kg] each mass

    StateSpace<4, 1, 3, double, 4, 3> sys{};     ///< motor force in; encoder, accelerometer, load position out
    Matrix<4, 1>                      G_load{};  ///< discrete entry of the load force
    Matrix<3, 1>                      Cd_load{}; ///< the accelerometer feels it directly

    CompliantDrive() {
        Matrix<4, 4> Ac{};
        Ac(0, 1) = 1.0;
        Ac(1, 0) = -k / m;
        Ac(1, 1) = -2.0;
        Ac(1, 2) = k / m;
        Ac(2, 3) = 1.0;
        Ac(3, 0) = k / m;
        Ac(3, 2) = -k / m;
        Ac(3, 3) = -2.0;
        StateSpace<4, 2, 3> cont{.A = Ac};
        cont.B(1, 0) = 1.0 / m; // motor force
        cont.B(3, 1) = 1.0 / m; // load force
        cont.C(0, 0) = 1.0;
        for (std::size_t c = 0; c < 4; ++c) {
            cont.C(1, c) = Ac(3, c);
        }
        cont.C(2, 2) = 1.0;
        const auto disc = discretize(cont, Ts, DiscretizationMethod::ZOH).value();
        sys.A = disc.A;
        sys.B = disc.B.block<4, 1>(0, 0);
        sys.C = disc.C;
        sys.Ts = Ts;
        G_load = disc.B.block<4, 1>(0, 1);
        Cd_load(1, 0) = 1.0 / m;
    }

    [[nodiscard]] ColVec<3> measure(const ColVec<4>& x, double load_force) const {
        return (sys.C * x) + (Cd_load * ColVec<1>{{load_force}});
    }
    [[nodiscard]] ColVec<4> step(const ColVec<4>& x, double u, double load_force) const {
        return (sys.A * x) + (sys.B * ColVec<1>{{u}}) + (G_load * ColVec<1>{{load_force}});
    }
};

design::MPCWeights<1, 1> load_weights() {
    design::MPCWeights<1, 1> w{};
    w.Qy = 1.0e6;
    w.Rdu = 1.0;
    return w;
}

/// The load force where it acts, felt by the accelerometer; the load position (row 2) held.
design::DisturbanceModel<4, 3, 1, 1> load_force_model(const CompliantDrive& plant) {
    return design::DisturbanceModel<4, 3, 1, 1>{.Bd = plant.G_load, .Cd = plant.Cd_load, .outputs = {{2}}};
}

/// The estimator's noise on the drive: the encoder exact, the accelerometer noisy, the load position unread.
design::KalmanNoise<4, 3, double, 1> drive_noise() {
    return design::KalmanNoise<4, 3, double, 1>{
        .process = Matrix<4, 4>::identity() * 1.0e-8,
        .measurement = Matrix<3, 3>{{1.0e-8, 0.0, 0.0}, {0.0, 1.0e-2, 0.0}, {0.0, 0.0, 1.0e4}},
        .disturbance = 1.0,
    };
}

} // namespace

TEST_CASE("offset-free MPC: a load force sensed by an accelerometer, the load held at r") {
    // Two sensors, one input, a disturbance where it acts, and a controlled output nobody
    // measures. The load-force model holds the load exactly at r. The input-disturbance
    // model, given the same controlled output, explains the load force as a motor force,
    // so it believes the spring relaxed and leaves the load off by F/k.
    const CompliantDrive plant;
    const auto           art_load = design::ofmpc<40, 8>(plant.sys, load_weights(), {}, drive_noise(), load_force_model(plant));
    const auto           art_input = design::ofmpc<40, 8>(plant.sys, load_weights(), {}, drive_noise(), design::InputDisturbance<1>{{2}});
    REQUIRE(art_load.success);
    REQUIRE(art_input.success);
    static_assert(decltype(OFMPC{art_load})::NDE == 1);
    static_assert(decltype(OFMPC{art_load})::NYM == 3);

    OFMPC           c_load{art_load};
    OFMPC           c_input{art_input};
    ColVec<4>       x_load{};
    ColVec<4>       x_input{};
    const ColVec<1> r{{0.01}};
    const double    force = -0.8;
    for (int i = 0; i < 1500; ++i) {
        x_load = plant.step(x_load, c_load.step(r, plant.measure(x_load, force))(0), force);
        x_input = plant.step(x_input, c_input.step(r, plant.measure(x_input, force))(0), force);
    }
    const double deflection = std::abs(force) / CompliantDrive::k;
    CHECK(x_load(2) == doctest::Approx(0.01).epsilon(1e-4));                         // offset-free
    CHECK(c_load.disturbance_estimate()(0) == doctest::Approx(force).epsilon(1e-3)); // and the force found
    CHECK(std::abs(x_input(2) - 0.01) > 0.5 * deflection);                           // ≈ F/k off
    CHECK(std::abs(c_input.state_estimate()(0) - c_input.state_estimate()(2)) < 0.1 * deflection);
}

TEST_CASE("offset-free MPC with a general disturbance model honours the input and rate limits") {
    const CompliantDrive         plant;
    design::MPCConstraints<1, 1> limits{};
    limits.u_lim = {
        .lower = -1.0,
        .upper = 1.0,
    };
    limits.du_lim = {
        .lower = -0.05,
        .upper = 0.05,
    };
    const auto art = design::ofmpc<40, 8>(plant.sys, load_weights(), limits, drive_noise(), load_force_model(plant));
    REQUIRE(art.success);

    OFMPC           ctrl{art};
    ColVec<4>       x{};
    const ColVec<1> r{{0.01}};
    const double    force = -0.8; // holding it takes u = 0.8 of the ±1 authority
    double          u_prev = 0.0;
    double          worst_u = 0.0;
    double          worst_du = 0.0;
    for (int i = 0; i < 2000; ++i) {
        const double u = ctrl.step(r, plant.measure(x, force))(0);
        worst_u = std::max(worst_u, std::abs(u));
        worst_du = std::max(worst_du, std::abs(u - u_prev));
        u_prev = u;
        x = plant.step(x, u, force);
    }
    CHECK(worst_u <= 1.0 + 1e-9);
    CHECK(worst_du <= 0.05 + 1e-9);
    CHECK(worst_du > 0.04); // the rate limit was reached, not idle
    CHECK(x(2) == doctest::Approx(0.01).epsilon(1e-4));
}

TEST_CASE("offset-free MPC on the ADMM policy reports an infeasible tick like the active-set default") {
    // A held input of 2.0 outside |u| <= 0.5, which |Δu| <= 0.1 cannot reach in one move.
    design::MPCConstraints<1, 1> limits{};
    limits.u_lim = {
        .lower = -0.5,
        .upper = 0.5,
    };
    limits.du_lim = {
        .lower = -0.1,
        .upper = 0.1,
    };
    limits.max_qp_iterations = 20000;
    const auto art = offset_free<15, 5>(move_weight(0.05), limits);
    REQUIRE(art.success);

    using Ctl = decltype(OFMPC{art});
    OFMPC           exact{art};
    OFMPC           admm{art, design::AdmmSolver<Ctl::NZ1, Ctl::NI1, double>{}};
    const ColVec<1> r{1.0};
    const ColVec<1> y{0.0};

    exact.back_calculate(ColVec<1>{2.0});
    admm.back_calculate(ColVec<1>{2.0});
    const auto ue = exact.step(r, y);
    const auto ua = admm.step(r, y);
    CHECK(exact.last_status() == design::QPStatus::Infeasible);
    CHECK(admm.last_status() == design::QPStatus::Infeasible);
    CHECK(admm.last_iterations() < 1000);
    CHECK(ue(0) == doctest::Approx(0.5));
    CHECK(ua(0) == doctest::Approx(0.5));

    exact.back_calculate(ColVec<1>{0.3});
    admm.back_calculate(ColVec<1>{0.3});
    const auto ue2 = exact.step(r, y);
    const auto ua2 = admm.step(r, y);
    CHECK(exact.last_status() == design::QPStatus::Success);
    CHECK(admm.last_status() == design::QPStatus::Success);
    CHECK(ua2(0) == doctest::Approx(ue2(0)).epsilon(1e-4));
}

TEST_CASE("offset-free MPC: InputDisturbance is the explicit model B_d = B, z = y") {
    const auto by_spec = offset_free<10, 4>(move_weight(0.05));
    const auto by_model = design::ofmpc<10, 4>(double_integrator, move_weight(0.05), {}, kUnitNoise, design::DisturbanceModel<2, 1, 1>{.Bd = double_integrator.B});
    static_assert(std::is_same_v<decltype(by_spec), decltype(by_model)>);
    REQUIRE(by_spec.success);
    REQUIRE(by_model.success);
    CHECK((by_spec.mpc.H - by_model.mpc.H).norm() == 0.0);
    CHECK((by_spec.mpc.F_x - by_model.mpc.F_x).norm() == 0.0);
    CHECK((by_spec.kalman.M - by_model.kalman.M).norm() == 0.0);
}

TEST_CASE("offset-free MPC rejects an undetectable or unsupported general disturbance model") {
    // Two disturbances that both enter exactly as the input does: the sensors cannot tell
    // them apart (rank [A−I, Bd; C, Cd] < NX + ND).
    const CompliantDrive                 plant;
    design::DisturbanceModel<4, 3, 2, 1> twins{.outputs = {{2}}};
    twins.Bd.block<4, 1>(0, 0) = plant.sys.B;
    twins.Bd.block<4, 1>(0, 1) = plant.sys.B;
    const auto twins_art = design::ofmpc<10, 4>(plant.sys, load_weights(), {}, design::KalmanNoise<4, 3, double, 2>{}, twins);
    CHECK_FALSE(twins_art.success);
    CHECK(twins_art.status == design::DesignStatus::not_detectable);

    // A controlled output with feedthrough is outside the condensed MPC (D_z = 0 only).
    auto feedthrough = plant.sys;
    feedthrough.D(2, 0) = 1.0;
    const auto refused = design::ofmpc<10, 4>(feedthrough, load_weights(), {}, drive_noise(), load_force_model(plant));
    CHECK_FALSE(refused.success);

    // A failed design still deploys: zero estimate, finite (zero) move.
    OFMPC      idle{refused};
    const auto u = idle.step(ColVec<1>{{0.01}}, ColVec<3>{});
    CHECK(damp::isfinite(u(0)));

    // A controlled output out of range is refused.
    CHECK(design::ofmpc<10, 4>(plant.sys, load_weights(), {}, drive_noise(), design::InputDisturbance<1>{{3}}).status == design::DesignStatus::bad_weights);
}

namespace {

/// The heater of examples/control/heater_ofmpc with its 2.5-sample dead time in the model, the
/// process noise on the block's temperature (G = e₀).
constexpr StateSpace<1, 1, 1, double, 1, 1> kBlock{.A = {{-1.0 / 30.0}}, .B = {{6.0 / 30.0}}, .C = {{1.0}}};

constexpr StateSpace<4, 1, 1, double, 1, 1> heater_model() {
    auto m = test::register_form<3>(kBlock, 2.0, 5.0);
    m.G = Matrix<4, 1>{{1.0}, {0.0}, {0.0}, {0.0}};
    return m;
}

constexpr auto kHeater = heater_model();

/// The heater's offset-free design: the power in [0, @p u_max], an input disturbance.
template<size_t NC>
auto heater_design(double u_max, bool terminal_weight = true) {
    design::MPCWeights<1, 1> w{};
    w.Qy = 1.0;
    w.Rdu = 10.0;
    if (terminal_weight) {
        w.Qy_terminal = w.Qy;
    }
    design::MPCConstraints<1, 1> c{};
    c.u_lim = {
        .lower = 0.0,
        .upper = u_max,
    };
    const auto art = design::ofmpc<30, NC>(kHeater, w, c, design::KalmanNoise{Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{0.05}}}, design::InputDisturbance{});
    REQUIRE(art.success);
    return art;
}

} // namespace

TEST_SUITE("offset-free MPC float runtime") {
    // At a 185-unit operating point the QP gradient F_x·x − F_r·r is a small difference of
    // large terms, and in float the lost digits left the heater example's power chattering
    // by 0.3 W every period at steady state. Posed about the steady-state target (x_t, u_t)
    // the gradient is a sum of small terms and the held input is constant.
    TEST_CASE("the held input does not chatter at a large operating point") {
        OFMPC mpc{heater_design<4>(40.0).as<float>()};

        // the model itself as the plant, in double
        ColVec<4>   x{};
        const float r = 185.0f;
        float       u_lo = 1e9f;
        float       u_hi = -1e9f;
        for (int k = 0; k < 400; ++k) {
            const float y = static_cast<float>((kHeater.C * x)(0, 0));
            const float u = mpc.step(ColVec<1, float>{r}, ColVec<1, float>{y})[0];
            x = kHeater.A * x + kHeater.B * ColVec<1>{static_cast<double>(u)};
            if (k >= 300) {
                u_lo = std::min(u_lo, u);
                u_hi = std::max(u_hi, u);
            }
        }
        CHECK(u_hi - u_lo < 5e-3f); // was ≈ 0.3 about the origin
        CHECK(u_hi == doctest::Approx(185.0 / 6.0).epsilon(1e-3));
    }

    // A Δu-form Hessian is graded: the first move drives every later output, so its diagonal
    // is 10⁴ times the last move's. Cholesky's pivot floor once scaled with the matrix size
    // (N·64·eps); with 30 moves in float it refused the early pivots, which carry only a few
    // terms' roundoff, every tick reported NotPositiveDefinite and the heater stayed at 0 W.
    TEST_CASE("a float runtime factors the Hessian of a 30-move control horizon") {
        const auto art = heater_design<30>(40.0, false);
        CHECK(mat::cholesky(art.mpc.H.as<float>()).has_value());

        OFMPC     mpc{art.as<float>()};
        ColVec<4> x{};
        float     u_peak = 0.0f;
        float     u = 0.0f;
        int       failed_solves = 0;
        for (int k = 0; k < 300; ++k) {
            const float y = static_cast<float>((kHeater.C * x)(0, 0));
            u = mpc.step(ColVec<1, float>{185.0f}, ColVec<1, float>{y})[0];
            failed_solves += (mpc.last_status() != design::QPStatus::Success) ? 1 : 0;
            u_peak = std::max(u_peak, u);
            x = kHeater.A * x + kHeater.B * ColVec<1>{static_cast<double>(u)};
        }
        CHECK(failed_solves == 0);
        CHECK(u_peak == doctest::Approx(40.0)); // full power during the heat-up
        CHECK(u == doctest::Approx(185.0 / 6.0).epsilon(1e-3));
    }
}

TEST_SUITE("offset-free MPC error-state estimator") {
    /// The tick in absolute coordinates, from public parts: the steady-state filter on [x; d] and
    /// an MPC on the measured state that takes the estimate as its state and d̂ as a measured
    /// disturbance (the same QP).
    template<typename T>
    struct AbsoluteReplica {
        static constexpr auto law_design(double u_max) {
            const StateSpace<4, 1, 1> plant{.A = kHeater.A, .B = kHeater.B, .C = kHeater.C, .Ts = kHeater.Ts};
            design::MPCWeights<1, 1>  w{};
            w.Rdu = 10.0;
            design::MPCConstraints<1, 1> c{};
            c.u_lim = {
                .lower = 0.0,
                .upper = u_max,
            };
            return design::mpc<30, 4>(plant, w, c, design::MeasuredDisturbance{kHeater.B});
        }
        using Law = decltype(MPC{law_design(1.0).template as<T>()});

        Law                    mpc;
        Luenberger<5, 1, 1, T> kf;
        ColVec<1, T>           u_prev{};
        template<typename Art>
        AbsoluteReplica(const Art& a, double u_max)
            : mpc(law_design(u_max).template as<T>()), kf(a.kalman) {}
        T tick(T r, T y) {
            kf.predict(u_prev);
            kf.update(ColVec<1, T>{y}, u_prev);
            ColVec<4, T> x{};
            for (size_t i = 0; i < 4; ++i) {
                x(i) = kf.state()(i);
            }
            const ColVec<1, T> d{kf.state()(4)};
            const ColVec<1, T> u = mpc.step(ColVec<1, T>{r}, d, x); // the estimate as the measured state, d̂ as a measured disturbance
            u_prev = u;
            return u(0);
        }
    };

    TEST_CASE("equal to the absolute-coordinate tick in double, through saturation and a load step") {
        const auto              art = heater_design<4>(40.0);
        OFMPC                   esk{art};
        AbsoluteReplica<double> abs{art, 40.0};
        ColVec<4>               xa{};
        ColVec<4>               xb{};
        double                  worst = 0.0;
        for (int k = 0; k < 400; ++k) {
            const double w = (k >= 200) ? -0.8 : 0.0; // an unmeasured heat loss
            const double ua = esk.step(ColVec<1>{185.0}, ColVec<1>{(kHeater.C * xa)(0, 0)})[0];
            const double ub = abs.tick(185.0, (kHeater.C * xb)(0, 0));
            worst = std::max(worst, std::abs(ua - ub));
            xa = kHeater.A * xa + kHeater.B * ColVec<1>{ua + w};
            xb = kHeater.A * xb + kHeater.B * ColVec<1>{ub + w};
        }
        CHECK(worst <= 1e-9 * 40.0);
        CHECK(std::abs((kHeater.C * xa)(0, 0) - 185.0) < 1e-6);
        // the public estimates are absolute
        CHECK(esk.state_estimate()(0) == doctest::Approx(185.0).epsilon(1e-6));
        CHECK(esk.disturbance_estimate()(0) == doctest::Approx(-0.8).epsilon(1e-4));
    }

    TEST_CASE("float at a large operating point: the deviation form holds the measurement's floor") {
        // Setpoint 1e5 (in whatever unit), input bound 5e4. A float holds 1e5 to 8e-3, and the
        // measurement y arrives as a float, so that resolution is the floor of any form. The
        // deviation form removes the filter recursion's and the gradient's own rounding at that
        // magnitude: the held input stays within a few ulp(y)·gain of 16667.
        OFMPC     esk{heater_design<4>(5.0e4).as<float>()};
        ColVec<4> xa{};
        float     lo_e = 1e30f;
        float     hi_e = -1e30f;
        for (int k = 0; k < 600; ++k) {
            const float ua = esk.step(ColVec<1, float>{1.0e5f}, ColVec<1, float>{static_cast<float>((kHeater.C * xa)(0, 0))})[0];
            xa = kHeater.A * xa + kHeater.B * ColVec<1>{static_cast<double>(ua)};
            if (k >= 500) {
                lo_e = std::min(lo_e, ua);
                hi_e = std::max(hi_e, ua);
            }
        }
        const float held = 1.0e5f / 6.0f;
        CHECK(hi_e - lo_e <= 3e-5f * held); // within a few ulp(y)·gain: the measurement's floor
        CHECK(std::abs((kHeater.C * xa)(0, 0) - 1.0e5) < 1.0);
    }
}

TEST_SUITE("offset-free MPC terminal") {
    TEST_CASE("default design: no terminal rows") {
        const auto art = offset_free<8, 4>({});
        REQUIRE(art.success);
        CHECK(art.mpc.n_terminal == 0);
        CHECK(art.mpc.terminal == design::MPCTerminal::none);
        static_assert(decltype(OFMPC{art})::NYM == 1); // the measurement is the reading y, not the state: a filter reconstructs x
        static_assert(decltype(OFMPC{art})::NDE == 1);
    }

    TEST_CASE("dare_cost fills a PD terminal Hessian") {
        design::MPCConstraints<1, 1> c{};
        c.terminal = design::MPCTerminal::dare_cost;
        const auto art = offset_free<8, 4>({}, c);
        REQUIRE(art.success);
        CHECK(art.mpc.n_terminal == 0);
        CHECK(art.mpc.terminal == design::MPCTerminal::dare_cost);
        CHECK(art.mpc.P(0, 0) != doctest::Approx(0.0));
        CHECK(mat::cholesky(art.mpc.P).has_value());
    }

    TEST_CASE("invariant_set with NT_max: the set moves with the estimated disturbance's target") {
        design::MPCConstraints<1, 1> limits{};
        limits.u_lim = {
            .lower = -1.0,
            .upper = 1.0,
        };
        limits.y_lim = {
            .lower = -1.0,
            .upper = 1.0,
        };
        limits.terminal = design::MPCTerminal::invariant_set;

        const auto art = offset_free<4, 4, 64>(move_weight(0.01), limits);
        REQUIRE(art.success);
        REQUIRE(art.mpc.n_terminal > 0);
        CHECK(art.mpc.terminal == design::MPCTerminal::invariant_set);

        OFMPC law{art};
        static_assert(decltype(law)::NI1 == (4 * 4) + (2 * 4) + 64 + 1);

        // Regulation from the reading under an unmeasured load: the terminal set sits about the
        // target the load estimate moves (u_t = −d̂), so every solve stays feasible.
        const ColVec<1> r{0.0};
        ColVec<2>       x{0.02, 0.0};
        for (size_t k = 0; k < 300; ++k) {
            const auto u = law.step(r, ColVec<1>{x(0)});
            CHECK(law.last_status() == design::QPStatus::Success);
            CHECK(u(0) <= 1.0 + 1e-9);
            CHECK(u(0) >= -1.0 - 1e-9);
            x = plant_step(x, u(0) + 0.2);
        }
        CHECK(x(0) == doctest::Approx(0.0).scale(1.0).epsilon(1e-3));
        CHECK(law.disturbance_estimate()(0) == doctest::Approx(0.2).epsilon(1e-3));
    }

    TEST_CASE("invariant_set without NT_max is refused as bad_weights") {
        design::MPCConstraints<1, 1> limits{};
        limits.terminal = design::MPCTerminal::invariant_set;
        const auto art = offset_free<4, 2>({}, limits);
        CHECK_FALSE(art.success);
        CHECK(art.status == design::DesignStatus::bad_weights);
        CHECK(art.mpc.status == design::DesignStatus::bad_weights);
    }
}

TEST_CASE("more controlled outputs than inputs: the target is an equilibrium, so d̂ stays 0 on an exact model") {
    //! The stacked least squares put the target off the equilibria: d̂ = −0.019 with no
    //! disturbance and y = (0.334, 0.334) against the measured-state MPC's (0.493, 0.493).
    const StateSpace<2, 1, 2, double, 2, 2> plant{
        .A = {{0.9, 0.0}, {0.0, 0.8}},
        .B = {{0.1}, {0.2}},
        .C = {{1.0, 0.0}, {0.0, 1.0}},
        .G = Matrix<2, 2>::identity(),
        .H = Matrix<2, 2>::identity(),
        .Ts = 0.1,
    };
    design::MPCWeights<1, 2> w{};
    w.Rdu = 0.1;
    const design::KalmanNoise<2, 2, double, 1> noise{.process = Matrix<2, 2>::identity() * 1e-4, .measurement = Matrix<2, 2>::identity() * 1e-4, .disturbance = 1e-4};
    const auto                                 art = design::ofmpc<20, 5>(plant, w, design::MPCConstraints<1, 2>{}, noise, design::InputDisturbance{});
    REQUIRE(art.success);
    const auto nominal = design::mpc<20, 5>(StateSpace<2, 1, 2>{.A = plant.A, .B = plant.B, .C = plant.C, .Ts = 0.1}, w, {});
    REQUIRE(nominal.success);
    OFMPC     law{art};
    MPC       reference{nominal};
    ColVec<2> x{};
    ColVec<2> x_ref{};
    ColVec<1> u{};
    for (size_t k = 0; k < 600; ++k) {
        u = law.step(ColVec<2>{1.0, 0.0}, x);
        x = ColVec<2>{(0.9 * x(0)) + (0.1 * u(0)), (0.8 * x(1)) + (0.2 * u(0))};
        const auto u_ref = reference.step(ColVec<2>{1.0, 0.0}, x_ref);
        x_ref = ColVec<2>{(0.9 * x_ref(0)) + (0.1 * u_ref(0)), (0.8 * x_ref(1)) + (0.2 * u_ref(0))};
    }
    CHECK(law.disturbance_estimate()(0) == doctest::Approx(0.0).scale(1.0).epsilon(1e-9));
    CHECK(x(0) == doctest::Approx(x_ref(0)).epsilon(1e-6));
    CHECK(x(1) == doctest::Approx(x_ref(1)).epsilon(1e-6));
}

TEST_CASE("back_calculate: a clamp is the applied input, the next move starts from it") {
    design::MPCConstraints<1, 1> limits{};
    limits.u_lim = {.lower = -0.3, .upper = 0.3};
    const auto art = offset_free<10, 4>(move_weight(0.05), limits);
    REQUIRE(art.success);
    OFMPC     ctrl{art};
    ColVec<2> x{};
    double    last = 0.0;
    for (int k = 0; k < 80; ++k) {
        const auto   u = ctrl.step(ColVec<1>{5.0}, ColVec<1>{x[0]});
        const double applied = damp::clamp(u[0], -0.3, 0.3);
        ctrl.back_calculate(ColVec<1>{applied});
        CHECK(applied >= -0.3);
        CHECK(applied <= 0.3);
        x = plant_step(x, applied);
        last = applied;
    }
    CHECK(damp::abs(last) == doctest::Approx(0.3));
}
