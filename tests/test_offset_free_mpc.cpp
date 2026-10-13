// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <type_traits>

#include "damp/controllers/mpc.hpp"
#include "damp/controllers/offset_free_lqg.hpp"
#include "damp/controllers/offset_free_mpc.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

constexpr double Ts = 0.1;

// Double integrator, ZOH-discretized at Ts = 0.1, position output only.
constexpr StateSpace<2, 1, 1> double_integrator{
    .A = {{1.0, Ts}, {0.0, 1.0}},
    .B = {{0.5 * Ts * Ts}, {Ts}},
    .C = {{1.0, 0.0}},
    .Ts = Ts,
};

constexpr ColVec<2> plant_step(const ColVec<2>& x, double u_effective) {
    return ColVec<2>{
        x(0) + (Ts * x(1)) + (0.5 * Ts * Ts * u_effective),
        x(1) + (Ts * u_effective),
    };
}

} // namespace

TEST_CASE("offset-free MPC rejects an unknown input disturbance from output feedback") {
    // The controller sees only y = position; the plant carries a constant
    // unmodeled load d. The disturbance estimate must converge to d and the
    // output must settle on the reference with no offset — the case the bare
    // velocity-form MPC cannot handle.
    design::MPCWeights<1, 1> weights{};
    weights.Rdu = Matrix<1, 1>{{0.05}};

    const auto art = design::mpc<15, 5>(
        double_integrator, weights, design::MPCConstraints<1, 1>{},
        design::OffsetFreeNoise{Matrix<2, 2>::identity(), Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{1.0}}}
    );
    REQUIRE(art.success);

    OffsetFreeMPC   controller{art};
    ColVec<2>       x{};
    const ColVec<1> r{0.5};
    const double    d_load = 0.3;

    for (size_t k = 0; k < 600; ++k) {
        const auto u = controller.control(r, ColVec<1>{x(0)});
        x = plant_step(x, u(0) + d_load);
    }

    CHECK(x(0) == doctest::Approx(0.5).epsilon(1e-3));
    CHECK(x(1) == doctest::Approx(0.0).epsilon(1e-3));
    CHECK(controller.disturbance_estimate()(0) == doctest::Approx(d_load).epsilon(1e-3));
    // Steady state: the applied input cancels the load.
    CHECK(controller.mpc.previous_control()(0) == doctest::Approx(-d_load).epsilon(1e-2));
}

TEST_CASE("offset-free MPC removes steady-state error under plant gain mismatch") {
    // The real plant has 30% more input gain than the model. The disturbance
    // integrator absorbs the mismatch at steady state.
    design::MPCWeights<1, 1> weights{};
    weights.Rdu = Matrix<1, 1>{{0.05}};

    const auto art = design::mpc<15, 5>(double_integrator, weights);
    REQUIRE(art.success);

    OffsetFreeMPC   controller{art};
    ColVec<2>       x{};
    const ColVec<1> r{1.0};

    for (size_t k = 0; k < 800; ++k) {
        const auto u = controller.control(r, ColVec<1>{x(0)});
        x = plant_step(x, 1.3 * u(0)); // plant gain != model gain
    }

    CHECK(x(0) == doctest::Approx(1.0).epsilon(5e-3));
    CHECK(x(1) == doctest::Approx(0.0).epsilon(5e-3));
}

TEST_CASE("offset-free MPC respects input limits while rejecting a load") {
    design::MPCWeights<1, 1> weights{};
    weights.Rdu = Matrix<1, 1>{{0.05}};

    design::MPCConstraints<1, 1> limits{};
    limits.u_min = ColVec<1>{-0.6};
    limits.u_max = ColVec<1>{0.6};

    const auto art = design::mpc<15, 5>(double_integrator, weights, limits);
    REQUIRE(art.success);

    OffsetFreeMPC   controller{art};
    ColVec<2>       x{};
    const ColVec<1> r{0.5};

    for (size_t k = 0; k < 800; ++k) {
        const auto u = controller.control(r, ColVec<1>{x(0)});
        CHECK(u(0) <= 0.6 + 1e-9);
        CHECK(u(0) >= -0.6 - 1e-9);
        x = plant_step(x, u(0) + 0.3); // load leaves ±0.3 of authority
    }
    CHECK(x(0) == doctest::Approx(0.5).epsilon(5e-3));
}

TEST_CASE("caller-sequenced predict/update/control matches the fused tick") {
    design::MPCWeights<1, 1> weights{};
    weights.Rdu = Matrix<1, 1>{{0.05}};

    const auto art = design::mpc<10, 4>(double_integrator, weights);
    REQUIRE(art.success);

    OffsetFreeMPC fused{art};
    OffsetFreeMPC sequenced{art};

    ColVec<2>       xf{0.2, -0.1};
    ColVec<2>       xs{0.2, -0.1};
    const ColVec<1> r{0.4};

    for (size_t k = 0; k < 50; ++k) {
        const auto uf = fused.control(r, ColVec<1>{xf(0)});

        sequenced.predict(sequenced.mpc.previous_control());
        sequenced.update(ColVec<1>{xs(0)}, sequenced.mpc.previous_control());
        const auto us = sequenced.control(r);

        CHECK(uf(0) == doctest::Approx(us(0)).epsilon(1e-12));
        xf = plant_step(xf, uf(0));
        xs = plant_step(xs, us(0));
    }
}

TEST_CASE("mpc rejects undetectable disturbance models") {
    // B = 0: the input disturbance never reaches the output — the
    // Pannocchia–Rawlings rank condition fails.
    StateSpace<1, 1, 1> plant{};
    plant.A = Matrix<1, 1>{{0.9}};
    plant.B = Matrix<1, 1>{{0.0}};
    plant.C = Matrix<1, 1>{{1.0}};
    plant.Ts = 0.1;

    CHECK(!design::mpc<5, 5>(plant).success);
}

TEST_CASE("mpc is constexpr") {
    constexpr auto art = design::mpc<3, 3>(double_integrator);
    static_assert(art.success);
}

TEST_CASE("offset-free MPC deploys in float") {
    design::MPCWeights<1, 1> weights{};
    weights.Rdu = Matrix<1, 1>{{0.05}};

    const auto art = design::mpc<10, 4>(double_integrator, weights);
    REQUIRE(art.success);

    OffsetFreeMPC          controller{art.as<float>()};
    ColVec<2, float>       x{};
    const ColVec<1, float> r{0.5F};

    for (size_t k = 0; k < 600; ++k) {
        const auto u = controller.control(r, ColVec<1, float>{x(0)});
        x = ColVec<2, float>{
            x(0) + (0.1F * x(1)) + (0.005F * (u(0) + 0.3F)),
            x(1) + (0.1F * (u(0) + 0.3F)),
        };
    }
    CHECK(x(0) == doctest::Approx(0.5).epsilon(0.01));
}

TEST_CASE("feedback / commit: a downstream clamp reaches the held input and the estimator") {
    design::MPCWeights<1, 1> weights{};
    weights.Rdu = Matrix<1, 1>{{0.05}};
    const auto art = design::mpc<15, 5>(
        double_integrator, weights, design::MPCConstraints<1, 1>{},
        design::OffsetFreeNoise{Matrix<2, 2>::identity(), Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{1.0}}}
    );
    REQUIRE(art.success);

    OffsetFreeMPC   fused{art};
    OffsetFreeMPC   split{art};
    const ColVec<1> r{0.5};
    ColVec<2>       x{};
    for (size_t k = 0; k < 5; ++k) {
        const ColVec<1> y{x(0)};
        const auto      u_fused = fused.control(r, y);
        const auto      u_req = split.feedback(r, y);
        CHECK(u_req(0) == doctest::Approx(u_fused(0)));
        CHECK(split.mpc.previous_control()(0) != doctest::Approx(u_req(0))); // not held until commit
        split.commit(u_req);
        CHECK(split.mpc.previous_control()(0) == doctest::Approx(u_req(0)));
        x = plant_step(x, u_fused(0));
    }

    // A tighter downstream clamp: commit the applied input, and the next move starts from it.
    const ColVec<1> y{x(0)};
    (void)split.feedback(ColVec<1>{5.0}, y);
    split.commit(ColVec<1>{0.01});
    CHECK(split.mpc.previous_control()(0) == doctest::Approx(0.01));
}

namespace {

/// Motor and load masses on a compliant transmission: x = [p_motor, v_motor, p_load, v_load].
/// The motor force is the input; a force on the load is the disturbance. Sensors: the motor
/// encoder and an accelerometer on the load, which feels the load force directly. Nothing
/// measures the load position, the one output that matters.
struct CompliantDrive {
    static constexpr double Ts = 0.01;
    static constexpr double k = 400.0; ///< [N/m] transmission stiffness
    static constexpr double m = 1.0;   ///< [kg] each mass

    StateSpace<4, 1, 2> sys{};     ///< motor force in; encoder, accelerometer out
    Matrix<4, 1>        G_load{};  ///< discrete entry of the load force
    Matrix<2, 1>        Cd_load{}; ///< the accelerometer feels it directly

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
        StateSpace<4, 2, 2> cont{.A = Ac};
        cont.B(1, 0) = 1.0 / m; // motor force
        cont.B(3, 1) = 1.0 / m; // load force
        cont.C(0, 0) = 1.0;
        for (std::size_t c = 0; c < 4; ++c) {
            cont.C(1, c) = Ac(3, c);
        }
        const auto disc = discretize(cont, Ts, DiscretizationMethod::ZOH).value();
        sys.A = disc.A;
        sys.B = disc.B.block<4, 1>(0, 0);
        sys.C = disc.C;
        sys.Ts = Ts;
        G_load = disc.B.block<4, 1>(0, 1);
        Cd_load(1, 0) = 1.0 / m;
    }

    /// Hold the load position (state 2).
    [[nodiscard]] static Matrix<1, 4> load_position() {
        Matrix<1, 4> Cz{};
        Cz(0, 2) = 1.0;
        return Cz;
    }

    [[nodiscard]] ColVec<2> measure(const ColVec<4>& x, double load_force) const {
        return (sys.C * x) + (Cd_load * ColVec<1>{{load_force}});
    }
    [[nodiscard]] ColVec<4> step(const ColVec<4>& x, double u, double load_force) const {
        return (sys.A * x) + (sys.B * ColVec<1>{{u}}) + (G_load * ColVec<1>{{load_force}});
    }
};

design::MPCWeights<1, 1> load_weights() {
    design::MPCWeights<1, 1> w{};
    w.Qy = Matrix<1, 1>{{1.0e6}};
    w.Rdu = Matrix<1, 1>{{1.0}};
    return w;
}

design::OffsetFreeModel<4, 1, 2, 1, 1> load_force_model(const CompliantDrive& plant) {
    design::OffsetFreeModel<4, 1, 2, 1, 1> model{};
    model.Bd = plant.G_load;
    model.Cd = plant.Cd_load;
    model.Cz = CompliantDrive::load_position();
    return model;
}

} // namespace

TEST_CASE("offset-free MPC: a load force sensed by an accelerometer, the load held at r") {
    // Two sensors, one input, a disturbance where it acts, and a controlled output nobody
    // measures. The load-force model holds the load exactly at r. The input-disturbance
    // model, given the same controlled output, explains the load force as a motor force,
    // so it believes the spring relaxed and leaves the load off by F/k.
    const CompliantDrive plant;
    const Matrix<4, 4>   Qn = Matrix<4, 4>::identity() * 1.0e-8;
    const Matrix<1, 1>   Qd{{1.0}};
    const Matrix<2, 2>   Rn{{1.0e-8, 0.0}, {0.0, 1.0e-2}};

    design::OffsetFreeModel<4, 1, 2, 1, 1> at_input{};
    at_input.Bd = plant.sys.B;
    at_input.Cz = CompliantDrive::load_position();

    const auto art_load = design::mpc<40, 8>(plant.sys, load_force_model(plant), load_weights(), {}, design::OffsetFreeNoise{Qn, Rn, Qd});
    const auto art_input = design::mpc<40, 8>(plant.sys, at_input, load_weights(), {}, design::OffsetFreeNoise{Qn, Rn, Qd});
    REQUIRE(art_load.success);
    REQUIRE(art_input.success);
    static_assert(std::is_same_v<decltype(OffsetFreeMPC{art_load}), OffsetFreeMPC<4, 1, 2, 40, 8, double, 1, 1>>);

    OffsetFreeMPC   c_load{art_load};
    OffsetFreeMPC   c_input{art_input};
    ColVec<4>       x_load{};
    ColVec<4>       x_input{};
    const ColVec<1> r{{0.01}};
    const double    force = -0.8;
    for (int i = 0; i < 1500; ++i) {
        x_load = plant.step(x_load, c_load.control(r, plant.measure(x_load, force))(0), force);
        x_input = plant.step(x_input, c_input.control(r, plant.measure(x_input, force))(0), force);
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
    limits.u_min = ColVec<1>{-1.0};
    limits.u_max = ColVec<1>{1.0};
    limits.du_min = ColVec<1>{-0.05};
    limits.du_max = ColVec<1>{0.05};
    const auto art = design::mpc<40, 8>(
        plant.sys, load_force_model(plant), load_weights(), limits, design::OffsetFreeNoise{Matrix<4, 4>::identity() * 1.0e-8, Matrix<2, 2>{{1.0e-8, 0.0}, {0.0, 1.0e-2}}, Matrix<1, 1>{{1.0}}}
    );
    REQUIRE(art.success);

    OffsetFreeMPC   ctrl{art};
    ColVec<4>       x{};
    const ColVec<1> r{{0.01}};
    const double    force = -0.8; // holding it takes u = 0.8 of the ±1 authority
    double          u_prev = 0.0;
    double          worst_u = 0.0;
    double          worst_du = 0.0;
    for (int i = 0; i < 2000; ++i) {
        const double u = ctrl.control(r, plant.measure(x, force))(0);
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

TEST_CASE("offset-free MPC: the default form is the input-disturbance model with z = y") {
    design::MPCWeights<1, 1> weights{};
    weights.Rdu = Matrix<1, 1>{{0.05}};
    const auto                             by_default = design::mpc<10, 4>(double_integrator, weights);
    design::OffsetFreeModel<2, 1, 1, 1, 1> input{};
    input.Bd = double_integrator.B;
    input.Cz = double_integrator.C;
    const auto by_model = design::mpc<10, 4>(double_integrator, input, weights);
    static_assert(std::is_same_v<decltype(by_default), decltype(by_model)>);
    REQUIRE(by_default.success);
    REQUIRE(by_model.success);
    CHECK((by_default.mpc.H - by_model.mpc.H).norm() == 0.0);
    CHECK((by_default.mpc.F_x - by_model.mpc.F_x).norm() == 0.0);
    CHECK((by_default.kalman.L - by_model.kalman.L).norm() == 0.0);
}

TEST_CASE("offset-free MPC rejects an undetectable or unsupported general disturbance model") {
    // Two disturbances that both enter exactly as the input does: the sensors cannot tell
    // them apart (rank [A−I, Bd; C, Cd] < NX + ND).
    const CompliantDrive                   plant;
    design::OffsetFreeModel<4, 1, 2, 2, 1> twins{};
    twins.Bd.block<4, 1>(0, 0) = plant.sys.B;
    twins.Bd.block<4, 1>(0, 1) = plant.sys.B;
    twins.Cz = CompliantDrive::load_position();
    CHECK_FALSE(design::mpc<10, 4>(plant.sys, twins).success);

    // A controlled output with feedthrough is outside the condensed MPC (D_z = 0 only).
    auto feedthrough = load_force_model(plant);
    feedthrough.Dz = Matrix<1, 1>{{1.0}};
    CHECK_FALSE(design::mpc<10, 4>(plant.sys, feedthrough).success);

    // A failed design still deploys: zero estimate, finite (zero) move.
    OffsetFreeMPC<4, 1, 2, 10, 4, double, 1, 1> idle{design::mpc<10, 4>(plant.sys, feedthrough)};
    const auto                                  u = idle.control(ColVec<1>{{0.01}}, ColVec<2>{{0.0}, {0.0}});
    CHECK(damp::isfinite(u(0)));
}

TEST_SUITE("OffsetFreeMPC float runtime") {
    // At a 185-unit operating point the QP gradient F_x·x − F_r·r is a small difference of
    // large terms, and in float the lost digits left the heater example's power chattering
    // by 0.3 W every period at steady state. Posed about the steady-state target (x_t, u_t)
    // the gradient is a sum of small terms and the held input is constant.
    TEST_CASE("the held input does not chatter at a large operating point") {
        // the heater example's plant: K = 6, τ = 30 s, 5 s dead time at Ts = 2 s
        const StateSpace<1, 1, 1> block{.A = {{-1.0 / 30.0}}, .B = {{6.0 / 30.0}}, .C = {{1.0}}};
        const auto                model = discretize(block, 2.0, InputDelay<3>{5.0});
        REQUIRE(model.has_value());
        design::MPCWeights<1, 1> w{};
        w.Qy = Matrix<1, 1>{{1.0}};
        w.Rdu = Matrix<1, 1>{{10.0}};
        w.Qy_terminal = w.Qy;
        design::MPCConstraints<1, 1> c{};
        c.u_min = ColVec<1>{0.0};
        c.u_max = ColVec<1>{40.0};
        const auto art = design::mpc<30, 4>(*model, w, c, design::OffsetFreeNoise{Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{0.05}}});
        REQUIRE(art.success);
        OffsetFreeMPC<4, 1, 1, 30, 4, float> mpc{art.as<float>()};

        // the model itself as the plant, in double
        ColVec<4>   x{};
        const float r = 185.0f;
        float       u_lo = 1e9f;
        float       u_hi = -1e9f;
        for (int k = 0; k < 400; ++k) {
            const float y = static_cast<float>((model->C * x)(0, 0));
            const float u = mpc.control(ColVec<1, float>{r}, ColVec<1, float>{y})[0];
            x = model->A * x + model->B * ColVec<1>{static_cast<double>(u)};
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
        const StateSpace<1, 1, 1> block{.A = {{-1.0 / 30.0}}, .B = {{6.0 / 30.0}}, .C = {{1.0}}};
        const auto                model = discretize(block, 2.0, InputDelay<3>{5.0});
        REQUIRE(model.has_value());
        design::MPCWeights<1, 1> w{};
        w.Qy = Matrix<1, 1>{{1.0}};
        w.Rdu = Matrix<1, 1>{{10.0}};
        design::MPCConstraints<1, 1> c{};
        c.u_min = ColVec<1>{0.0};
        c.u_max = ColVec<1>{40.0};
        const auto art = design::mpc<30, 30>(*model, w, c, design::OffsetFreeNoise{Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{0.05}}});
        REQUIRE(art.success);
        CHECK(mat::cholesky(art.mpc.H.as<float>()).has_value());

        OffsetFreeMPC mpc{art.as<float>()};
        ColVec<4>     x{};
        float         u_peak = 0.0f;
        float         u = 0.0f;
        int           failed_solves = 0;
        for (int k = 0; k < 300; ++k) {
            const float y = static_cast<float>((model->C * x)(0, 0));
            u = mpc.control(ColVec<1, float>{185.0f}, ColVec<1, float>{y})[0];
            failed_solves += (mpc.mpc.last_status() != design::QPStatus::Success) ? 1 : 0;
            u_peak = std::max(u_peak, u);
            x = model->A * x + model->B * ColVec<1>{static_cast<double>(u)};
        }
        CHECK(failed_solves == 0);
        CHECK(u_peak == doctest::Approx(40.0)); // full power during the heat-up
        CHECK(u == doctest::Approx(185.0 / 6.0).epsilon(1e-3));
    }
}

TEST_SUITE("OffsetFreeMPC error-state estimator") {
    // The heater of examples/control/heater_ofmpc with its 2.5-sample dead time in the model.
    static constexpr StateSpace<1, 1, 1> block{.A = {{-1.0 / 30.0}}, .B = {{6.0 / 30.0}}, .C = {{1.0}}};
    static constexpr auto                model = discretize(block, 2.0, InputDelay<3>{5.0});

    template<typename T>
    design::OffsetFreeMPCArtifacts<4, 1, 1, 30, 4, T> heater_design(double u_max) {
        design::MPCWeights<1, 1> w{};
        w.Qy = Matrix<1, 1>{{1.0}};
        w.Rdu = Matrix<1, 1>{{10.0}};
        w.Qy_terminal = w.Qy;
        design::MPCConstraints<1, 1> c{};
        c.u_min = ColVec<1>{0.0};
        c.u_max = ColVec<1>{u_max};
        const auto art = design::mpc<30, 4>(*model, w, c, design::OffsetFreeNoise{Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{0.01}}, Matrix<1, 1>{{0.05}}});
        REQUIRE(art.success);
        return art.template as<T>();
    }

    /// The absolute-coordinate tick this replaces, from public parts: steady-state filter on
    /// [x; d], target from the maps, MPC about that target.
    template<typename T>
    struct AbsoluteReplica {
        MPC<4, 1, 1, 30, 4, 1, T>                 mpc;
        SteadyStateKalmanFilter<5, 1, 1, T, 5, 1> kf;
        Matrix<5, 1, T>                           target_d;
        Matrix<5, 1, T>                           target_r;
        ColVec<1, T>                              u_prev{};
        explicit AbsoluteReplica(const design::OffsetFreeMPCArtifacts<4, 1, 1, 30, 4, T>& a)
            : mpc(a.mpc), kf(a.kalman), target_d(a.target_d), target_r(a.target_r) {}
        T tick(T r, T y) {
            kf.predict(u_prev);
            kf.update(ColVec<1, T>{y}, u_prev);
            ColVec<4, T> x{};
            for (size_t i = 0; i < 4; ++i) {
                x(i) = kf.state()(i);
            }
            const ColVec<1, T> d{kf.state()(4)};
            const ColVec<5, T> s = (target_d * d) + (target_r * ColVec<1, T>{r});
            ColVec<4, T>       x_t{};
            for (size_t i = 0; i < 4; ++i) {
                x_t(i) = s(i);
            }
            const ColVec<1, T> u_t{s(4)};
            const ColVec<1, T> u = mpc.feedback(ColVec<1, T>{r}, x, d, x_t, u_t);
            mpc.commit(u);
            u_prev = u;
            return u(0);
        }
    };

    TEST_CASE("equal to the absolute-coordinate tick in double, through saturation and a load step") {
        const auto                            art = heater_design<double>(40.0);
        OffsetFreeMPC<4, 1, 1, 30, 4, double> esk{art};
        AbsoluteReplica<double>               abs{art};
        ColVec<4>                             xa{};
        ColVec<4>                             xb{};
        double                                worst = 0.0;
        for (int k = 0; k < 400; ++k) {
            const double w = (k >= 200) ? -0.8 : 0.0; // an unmeasured heat loss
            const double ua = esk.control(ColVec<1>{185.0}, ColVec<1>{(model->C * xa)(0, 0)})[0];
            const double ub = abs.tick(185.0, (model->C * xb)(0, 0));
            worst = std::max(worst, std::abs(ua - ub));
            xa = model->A * xa + model->B * ColVec<1>{ua + w};
            xb = model->A * xb + model->B * ColVec<1>{ub + w};
        }
        CHECK(worst <= 1e-9 * 40.0);
        CHECK(std::abs((model->C * xa)(0, 0) - 185.0) < 1e-6);
        // the public estimates are absolute
        CHECK(esk.state_estimate()(0) == doctest::Approx(185.0).epsilon(1e-6));
        CHECK(esk.disturbance_estimate()(0) == doctest::Approx(-0.8).epsilon(1e-4));
    }

    TEST_CASE("float at a large operating point: the deviation form adds no jitter of its own") {
        // Setpoint 1e5 (in whatever unit), input bound 5e4. A float holds 1e5 to 8e-3, and the
        // measurement y arrives as a float, so that resolution is a floor both forms share. The
        // deviation form removes the filter recursion's and the gradient's own rounding at that
        // magnitude; the absolute form adds them on top. Measured: 0.21 vs 0.27 on a held input
        // of 16667 (1.3e-5 relative). The gap widens with a measurement path that carries more
        // than float, which the stored deviation then uses in full.
        const auto                           art = heater_design<float>(5.0e4);
        OffsetFreeMPC<4, 1, 1, 30, 4, float> esk{art};
        AbsoluteReplica<float>               abs{art};
        ColVec<4>                            xa{};
        ColVec<4>                            xb{};
        float                                lo_e = 1e30f, hi_e = -1e30f, lo_a = 1e30f, hi_a = -1e30f;
        for (int k = 0; k < 600; ++k) {
            const float ua = esk.control(ColVec<1, float>{1.0e5f}, ColVec<1, float>{static_cast<float>((model->C * xa)(0, 0))})[0];
            const float ub = abs.tick(1.0e5f, static_cast<float>((model->C * xb)(0, 0)));
            xa = model->A * xa + model->B * ColVec<1>{static_cast<double>(ua)};
            xb = model->A * xb + model->B * ColVec<1>{static_cast<double>(ub)};
            if (k >= 500) {
                lo_e = std::min(lo_e, ua);
                hi_e = std::max(hi_e, ua);
                lo_a = std::min(lo_a, ub);
                hi_a = std::max(hi_a, ub);
            }
        }
        const float held = 1.0e5f / 6.0f;
        CHECK(hi_e - lo_e <= 3e-5f * held); // within a few ulp(y)·gain: the measurement's floor
        CHECK(hi_e - lo_e <= hi_a - lo_a);  // and no worse than the absolute form
        CHECK(std::abs((model->C * xa)(0, 0) - 1.0e5) < 1.0);
    }
}
