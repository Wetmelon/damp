// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>

#include "damp/controllers/mpc.hpp"
#include "damp/controllers/offset_free_mpc.hpp"
#include "damp/matrix/matrix.hpp"
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
        Matrix<2, 2>::identity(), Matrix<1, 1>{{1.0}}, Matrix<1, 1>{{0.01}}
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

    CHECK(!design::mpc<5>(plant).success);
}

TEST_CASE("mpc is constexpr") {
    constexpr auto art = design::mpc<3>(double_integrator);
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
