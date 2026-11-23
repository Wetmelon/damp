// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// An MPC's run-time limits are members, u_lim, du_lim, y_lim and x_lim, as the LQ laws' are: an
// assignment is taken into the QP's bound rows by the next tick, which rebuilds nothing when none
// was made and refuses one that is not a box or that moves an InvariantSet's limits. The recorded
// commands are those set_limits() gave for the same changes, before the members replaced it.

#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/controllers/state_feedback.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/input_limits.hpp"
#include "damp/design/qp.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/toolbox/bounds.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

constexpr double Ts = 0.1;

// Double integrator, ZOH-discretized at Ts = 0.1, the state measured, the position the output.
constexpr StateSpace double_integrator{
    .A = {{1.0, Ts}, {0.0, 1.0}},
    .B = {{0.5 * Ts * Ts}, {Ts}},
    .C = {{1.0, 0.0}},
    .Ts = Ts,
};

// The same plant read through its position, with process noise on both states and noise on the reading.
constexpr StateSpace<2, 1, 1, double, 2, 1> measured{
    .A = {{1.0, Ts}, {0.0, 1.0}},
    .B = {{0.5 * Ts * Ts}, {Ts}},
    .C = {{1.0, 0.0}},
    .Ts = Ts,
};

design::MPCWeights<1, 1> move_weight(double rdu) {
    design::MPCWeights<1, 1> w{};
    w.Rdu = rdu;
    return w;
}

constexpr ColVec<2> plant_step(const ColVec<2>& x, double u) {
    return ColVec<2>{
        x(0) + (Ts * x(1)) + (0.5 * Ts * Ts * u),
        x(1) + (Ts * u),
    };
}

constexpr design::InputLimits<1> kHalf{.u_lim = {.lower = -0.5, .upper = 0.5}};

/// 1 while the limit member @p lim is assigned and not yet taken in, else 0.
template<typename Lim>
unsigned assigned(const Lim& lim) {
    return detail::MPCBoundsAccess::assigned(lim);
}

template<typename L>
concept writes_an_entry = requires(L& law) { law.u_lim.lower(0) = 1.0; };

template<typename L>
concept writes_an_end = requires(L& law) { law.u_lim.lower = law.u_lim.upper; };

template<typename L>
concept assigns_a_box = requires(L& law) { law.u_lim = Bounds<1>{}; };

template<typename L>
concept reads_output_limits = requires(const L& law) { law.y_lim.lower(0); };

template<typename L>
concept reads_state_limits = requires(const L& law) { law.x_lim.lower(0); };

// An offset-free law on the design's constants, its output bound soft.
constexpr design::OutputLimits kSoftCap{.outputs = {0}, .lim = {.upper = 0.8}};
constexpr auto                 kOnConstants = design::mpc<3, 3>(measured, {}, design::KalmanNoise<2, 1, double, 1>{}, design::InputDisturbance{}, kSoftCap);

} // namespace

TEST_CASE("an MPC's limits are members: read as a Bounds, assigned whole, never an entry alone") {
    const auto art = design::mpc<20, 8>(double_integrator, move_weight(0.05), kHalf);
    REQUIRE(art.success);
    MPC law{art};
    using Runtime = decltype(law);
    static_assert(assigns_a_box<Runtime> && !writes_an_entry<Runtime> && !writes_an_end<Runtime>);
    static_assert(sizeof(detail::MPCBounds<1, float>) <= sizeof(Bounds<1, float>) + alignof(float)); // the box and one byte
    static_assert(sizeof(detail::MPCBounds<3, double>) <= sizeof(Bounds<3, double>) + alignof(double));

    //! The design's limits, read as a Bounds' ends, entry by entry, as the vector, or whole.
    CHECK(law.u_lim.lower(0) == -0.5);
    CHECK(law.u_lim.upper(0) == 0.5);
    const ColVec<1>& upper = law.u_lim.upper;
    CHECK(upper(0) == 0.5);
    const Bounds<1> box = law.u_lim;
    CHECK(box.lower(0) == -0.5);
    CHECK(law.du_lim.upper(0) == std::numeric_limits<double>::max());
    CHECK(assigned(law.u_lim) == 0U);

    //! Braces as a Bounds takes them, an end left out unbounded, and `= {}` both.
    law.u_lim = {.upper = 0.3};
    CHECK(law.u_lim.lower(0) == -std::numeric_limits<double>::max());
    CHECK(law.u_lim.upper(0) == 0.3);
    CHECK(assigned(law.u_lim) == 1U);
    law.du_lim = {.lower = {-0.1}, .upper = {0.1}};
    CHECK(law.du_lim.upper(0) == 0.1);
    law.du_lim = {};
    CHECK(law.du_lim.upper(0) == std::numeric_limits<double>::max());

    //! The next tick takes them in, and the mark is cleared.
    static_cast<void>(law.step(ColVec<1>{1.0}, ColVec<2>{}));
    CHECK(law.last_status() == design::QPStatus::Success);
    CHECK(assigned(law.u_lim) == 0U);
    CHECK(assigned(law.du_lim) == 0U);
}

TEST_CASE("assigning u_lim gives the commands set_limits gave") {
    const auto art = design::mpc<20, 8>(double_integrator, move_weight(0.05), kHalf);
    REQUIRE(art.success);
    MPC law{art};
    law.u_lim = {.lower = -0.5, .upper = 0.2};
    ColVec<2> x{};
    for (size_t k = 0; k < 400; ++k) {
        const auto u = law.step(ColVec<1>{1.0}, x);
        REQUIRE(law.last_status() == design::QPStatus::Success);
        CHECK(u(0) <= 0.2 + 1e-9);
        if (k == 0 || k == 1 || k == 5 || k == 20) {
            CHECK(u(0) == doctest::Approx(0.2).epsilon(1e-12));
        }
        if (k == 100) {
            CHECK(u(0) == doctest::Approx(-1.3812749484849882e-06).epsilon(1e-8));
        }
        x = plant_step(x, u(0));
    }
    CHECK(x(0) == doctest::Approx(1.0).epsilon(1e-9));
}

TEST_CASE("assigning u_lim and du_lim while running gives the commands set_limits gave, from the next tick") {
    const auto art = design::mpc<20, 8>(double_integrator, move_weight(0.05), kHalf);
    REQUIRE(art.success);
    MPC       law{art};
    ColVec<2> x{};
    for (size_t k = 0; k < 90; ++k) {
        if (k == 30) {
            law.u_lim = {.lower = -0.3, .upper = 0.3};
            law.du_lim = {.lower = -0.05, .upper = 0.05};
        }
        const auto u = law.step(ColVec<1>{1.0}, x);
        //! At k = 30 no move within Δu's 0.05 brings −0.445 into ±0.3: infeasible, held, and clamped there, as with set_limits.
        CHECK(law.last_status() == ((k == 30) ? design::QPStatus::Infeasible : design::QPStatus::Success));
        if (k == 29) {
            CHECK(u(0) == doctest::Approx(-0.44544704979984895).epsilon(1e-10)); // the design's limits
        }
        if (k == 30) {
            CHECK(u(0) == doctest::Approx(-0.3).epsilon(1e-12)); // the tick after the assignment: the new ones
        }
        if (k == 31) {
            CHECK(u(0) == doctest::Approx(-0.25).epsilon(1e-12)); // a move of 0.05
        }
        if (k == 40) {
            CHECK(u(0) == doctest::Approx(0.2).epsilon(1e-12));
        }
        if (k == 89) {
            CHECK(u(0) == doctest::Approx(0.00016597569284359376).epsilon(1e-6));
        }
        x = plant_step(x, u(0));
    }
    CHECK(x(0) == doctest::Approx(0.9999907080932916).epsilon(1e-10));
    CHECK(x(1) == doctest::Approx(-4.7474449907281548e-06).epsilon(1e-5));
}

TEST_CASE_TEMPLATE("assigning x_lim gives the commands set_limits gave", T, double, float) {
    const design::InputLimits<1> limits{.u_lim = {.lower = -1.0, .upper = 1.0}};
    design::StateLimits          speed{.states = {1}, .lim = {.lower = -0.3, .upper = 0.3}};
    speed.upper_ecr(0) = 0.0;
    const auto art = design::mpc<20, 8>(double_integrator, move_weight(0.01), limits, speed);
    REQUIRE(art.success);
    MPC law{art.template as<T>()};
    CHECK(law.x_lim.upper(0) == T(0.3));
    law.x_lim = {.lower = T(-0.3), .upper = T(0.15)};

    struct Recorded {
        size_t k;
        double u_double;
        double u_float;
    };
    constexpr Recorded recorded[] = {
        {0, 1.0, 1.0},
        {1, 0.49999999999999989, 0.49999994039535522},
        {65, -0.17390959825319985, -0.17387717962265015},
        {66, -0.31881961938157355, -0.31884825229644775},
        {67, -0.36821262888244399, -0.36826425790786743},
        {70, -0.13672324230895161, -0.13671565055847168},
        {75, 0.031794247184645026, 0.03181752935051918},
        {80, -0.0066096021613713841, -0.0066298586316406727},
    };
    constexpr bool is_float = std::is_same_v<T, float>;
    ColVec<2>      x{};
    double         peak = 0.0;
    size_t         next = 0;
    for (size_t k = 0; k < 300; ++k) {
        const auto u = law.step(ColVec<1, T>{T{1}}, x.template as<T>());
        CHECK(law.last_status() == design::QPStatus::Success);
        if (next < std::size(recorded) && recorded[next].k == k) {
            CHECK(static_cast<double>(u(0)) == doctest::Approx(is_float ? recorded[next].u_float : recorded[next].u_double).epsilon(is_float ? 1e-4 : 1e-9)); // float: its rounding, fused or not
            ++next;
        }
        x = plant_step(x, static_cast<double>(u(0)));
        peak = std::fmax(peak, x(1));
    }
    CHECK(next == std::size(recorded));
    CHECK(peak <= 0.15 + (is_float ? 1e-5 : 1e-9)); // hard: plant == model
    CHECK(x(0) == doctest::Approx(1.0).epsilon(is_float ? 1e-6 : 1e-12));
}

TEST_CASE("assigning y_lim gives the commands set_limits gave, on the outputs OutputLimits lists in its order") {
    constexpr double              big = damp::numeric_limits<double>::max();
    constexpr StateSpace<4, 2, 4> two_axes{
        .A = {{1.0, Ts, 0.0, 0.0}, {0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, Ts}, {0.0, 0.0, 0.0, 1.0}},
        .B = {{0.5 * Ts * Ts, 0.0}, {Ts, 0.0}, {0.0, 0.5 * Ts * Ts}, {0.0, Ts}},
        .C = Matrix<4, 4>::identity(),
        .Ts = Ts,
    };
    design::MPCWeights<2, 4> w{};
    w.Qy = diag(1.0, 0.0, 1.0, 0.0);
    w.Qy_terminal = w.Qy;
    w.Rdu = Matrix<2, 2>::identity() * 0.1;
    const design::InputLimits<2> limits{.u_lim = {.lower = {-1.0, -1.0}, .upper = {1.0, 1.0}}};
    design::OutputLimits         subset{.outputs = {1, 2}, .lim = {.lower = {-0.3, -big}, .upper = {0.3, 0.4}}};
    subset.upper_ecr(1) = 0.0;
    const auto art = design::mpc<12, 4>(two_axes, w, limits, subset);
    REQUIRE(art.success);
    MPC law{art};
    static_assert(reads_output_limits<decltype(law)> && !reads_state_limits<decltype(law)>);
    law.y_lim = {.lower = law.y_lim.lower, .upper = {0.25, 0.35}}; // the first axis's speed, then the second's position

    ColVec<4>       x{};
    const ColVec<4> r{1.0, 0.0, 0.5, 0.0};
    double          v_peak = 0.0;
    double          p_peak = 0.0;
    for (size_t k = 0; k < 60; ++k) {
        const auto u = law.step(r, x);
        REQUIRE(law.last_status() == design::QPStatus::Success);
        if (k == 3) {
            CHECK(u(0) == doctest::Approx(-1.975050111013843e-06).epsilon(1e-6));
            CHECK(u(1) == doctest::Approx(1.0).epsilon(1e-12));
        }
        if (k == 6) {
            CHECK(u(1) == doctest::Approx(-0.57449642841255633).epsilon(1e-9));
        }
        if (k == 10) {
            CHECK(u(1) == doctest::Approx(-0.56529296992136513).epsilon(1e-9));
        }
        if (k == 59) {
            CHECK(u(0) == doctest::Approx(-0.0096292158593440352).epsilon(1e-8));
        }
        x = (two_axes.A * x) + (two_axes.B * u);
        v_peak = damp::max(v_peak, x(1));
        p_peak = damp::max(p_peak, x(2));
    }
    CHECK(v_peak == doctest::Approx(0.25001996662729153).epsilon(1e-9)); // soft: about the new bound
    CHECK(p_peak == doctest::Approx(0.35).epsilon(1e-12));               // hard: on it
}

TEST_CASE("an OFMPC's u_lim assigned while running gives the commands set_limits gave") {
    const auto art = design::mpc<15, 5>(measured, move_weight(0.05), design::InputLimits{.u_lim = {.lower = -0.6, .upper = 0.6}}, design::KalmanNoise<2, 1, double, 1>{}, design::InputDisturbance{});
    REQUIRE(art.success);
    OFMPC     law{art};
    ColVec<2> x{};
    for (size_t k = 0; k < 300; ++k) {
        if (k == 100) {
            law.u_lim = {.lower = -0.45, .upper = 0.45}; // with the load, ±0.15 left to move on
        }
        const auto u = law.step(ColVec<1>{(k < 100) ? 0.5 : 2.0}, ColVec<1>{x(0)});
        REQUIRE(law.last_status() == design::QPStatus::Success);
        if (k == 99) {
            CHECK(u(0) == doctest::Approx(-0.3000841481877598).epsilon(1e-10));
        }
        if (k == 100 || k == 101 || k == 105) {
            CHECK(u(0) == doctest::Approx(0.45).epsilon(1e-12));
        }
        if (k == 120 || k == 150) {
            CHECK(u(0) == doctest::Approx(-0.45).epsilon(1e-12));
        }
        if (k == 299) {
            CHECK(u(0) == doctest::Approx(-0.30000094623013929).epsilon(1e-10));
        }
        x = plant_step(x, u(0) + 0.3);
    }
    CHECK(x(0) == doctest::Approx(2.0000000069048185).epsilon(1e-10));
}

TEST_CASE("Law<kRes, T> of an OFMPC on the design's constants takes its limit members as OFMPC does") {
    Law<kOnConstants, double> law{};
    OFMPC                     plain{kOnConstants};
    law.u_lim = {.lower = -0.4, .upper = 0.4};
    law.y_lim = {.upper = 0.5};
    plain.u_lim = law.u_lim; // one member from another: assigned too
    plain.y_lim = law.y_lim;
    CHECK(assigned(plain.u_lim) == 1U);

    ColVec<2> x{};
    for (size_t k = 0; k < 50; ++k) {
        const auto u = law.step(ColVec<1>{0.6}, ColVec<1>{x(0)});
        REQUIRE(law.last_status() == design::QPStatus::Success);
        CHECK(plain.step(ColVec<1>{0.6}, ColVec<1>{x(0)})(0) == doctest::Approx(u(0)).epsilon(1e-12));
        if (k == 0) {
            CHECK(u(0) == doctest::Approx(0.041881765682602101).epsilon(1e-10));
        }
        if (k == 5) {
            CHECK(u(0) == doctest::Approx(0.24562769697553119).epsilon(1e-10));
        }
        if (k == 10) {
            CHECK(u(0) == doctest::Approx(0.4).epsilon(1e-12));
        }
        if (k == 20) {
            CHECK(u(0) == doctest::Approx(-0.4).epsilon(1e-12));
        }
        if (k == 49) {
            CHECK(u(0) == doctest::Approx(-0.17746303709160807).epsilon(1e-8));
        }
        x = plant_step(x, u(0));
    }
    CHECK(x(0) == doctest::Approx(0.086350121397227309).epsilon(1e-8));
    CHECK(law.y_lim.upper(0) == 0.5);
}

TEST_CASE("a tick with nothing assigned leaves the QP's bound rows as built") {
    //! The result's u rows say 0.2 where its limits say 0.5. Ticks with nothing assigned plan to
    //! the rows: rebuilt from the members, they would plan to 0.5. Assigning the very box the law
    //! holds has the next tick rewrite them.
    auto art = design::mpc<20, 8>(double_integrator, move_weight(0.05), kHalf);
    REQUIRE(art.success);
    constexpr size_t NZ = 8; // NC moves of one input
    REQUIRE(art.mpc.scale_u(0) == 1.0);
    for (size_t i = 0; i < NZ; ++i) {
        REQUIRE(art.mpc.b_bound((2 * NZ) + i) == 0.5); // the u ≤ upper rows, after Δu's
        art.mpc.b_bound((2 * NZ) + i) = 0.2;
    }
    MPC       law{art};
    ColVec<2> x{};
    for (size_t k = 0; k < 10; ++k) {
        CHECK(assigned(law.u_lim) + assigned(law.du_lim) == 0U);
        if (k == 5) {
            law.u_lim = Bounds<1>(law.u_lim);
        }
        const auto u = law.step(ColVec<1>{1.0}, x);
        CHECK(law.last_status() == design::QPStatus::Success);
        CHECK(u(0) == doctest::Approx((k < 5) ? 0.2 : 0.5).epsilon(1e-12));
        x = plant_step(x, u(0));
    }
}

TEST_CASE("a limit that is not a box is refused at the next tick: reported once, the limits in force kept and read back") {
    const auto art = design::mpc<20, 8>(double_integrator, move_weight(0.05), kHalf, design::OutputLimits{.outputs = {0}, .lim = {.upper = 1.02}});
    REQUIRE(art.success);
    MPC        law{art};
    MPC        twin{art}; // never assigned
    ColVec<2>  x{};
    const auto tick = [&] {
        const auto u = law.step(ColVec<1>{1.0}, x);
        CHECK(twin.step(ColVec<1>{1.0}, x)(0) == doctest::Approx(u(0)).epsilon(1e-12)); // the same QP: the rows kept the limits in force
        x = plant_step(x, u(0));
    };
    const double nan = std::numeric_limits<double>::quiet_NaN();

    //! An inverted u_lim with a good du_lim: the tick takes neither, as one set_limits call.
    law.u_lim = {.lower = 0.3, .upper = -0.3};
    law.du_lim = {.lower = -0.1, .upper = 0.1};
    tick();
    CHECK(law.last_status() == design::QPStatus::LimitsRefused);
    CHECK(law.u_lim.lower(0) == -0.5);
    CHECK(law.u_lim.upper(0) == 0.5);
    CHECK(law.du_lim.upper(0) == std::numeric_limits<double>::max());
    CHECK(assigned(law.u_lim) + assigned(law.du_lim) == 0U);
    tick();
    CHECK(law.last_status() == design::QPStatus::Success); // reported on that tick alone

    //! A NaN end, on an output limit.
    law.y_lim = {.upper = nan};
    tick();
    CHECK(law.last_status() == design::QPStatus::LimitsRefused);
    CHECK(law.y_lim.upper(0) == 1.02);
    for (size_t k = 3; k < 40; ++k) {
        tick();
        CHECK(law.last_status() == design::QPStatus::Success);
    }
}

TEST_CASE("with an InvariantSet a limit that moves is refused; the same box, or a copy of the whole law, is not") {
    design::OutputLimits<1> position{.outputs = {0}, .lim = {.lower = -1.0, .upper = 1.0}};
    position.lower_ecr(0) = 0.0;
    position.upper_ecr(0) = 0.0;
    const auto art = design::mpc<4, 4, 64>(double_integrator, move_weight(0.01), design::InputLimits{.u_lim = {.lower = -1.0, .upper = 1.0}}, position, design::InvariantSet{});
    REQUIRE(art.success);
    REQUIRE(art.mpc.n_terminal > 0);
    MPC             law{art};
    MPC             twin{art};
    const ColVec<2> x0{0.5, 0.8}; // a fast approach the set keeps feasible
    ColVec<2>       x = x0;
    const ColVec<1> r{0.9};
    const auto      tick = [&] {
        const auto u = law.step(r, x);
        CHECK(twin.step(r, x)(0) == doctest::Approx(u(0)).epsilon(1e-12));
        x = plant_step(x, u(0));
    };

    law.u_lim = {.lower = -0.8, .upper = 0.8};
    tick();
    CHECK(law.last_status() == design::QPStatus::LimitsRefused);
    CHECK(law.u_lim.upper(0) == 1.0);
    law.y_lim = {.lower = -1.0, .upper = 0.9};
    tick();
    CHECK(law.last_status() == design::QPStatus::LimitsRefused);
    CHECK(law.y_lim.upper(0) == 1.0);

    law.u_lim = {.lower = -1.0, .upper = 1.0}; // the box in force
    tick();
    CHECK(law.last_status() == design::QPStatus::Success);
    law = MPC{art}; // a whole law, its limits with it
    twin = MPC{art};
    x = x0;
    for (size_t k = 0; k < 30; ++k) {
        tick();
        CHECK(law.last_status() == design::QPStatus::Success);
    }
}

TEST_CASE("y_lim and x_lim take no storage without their options") {
    const auto plain = design::mpc<10, 4>(double_integrator, move_weight(0.05), kHalf);
    const auto bounded = design::mpc<10, 4>(double_integrator, move_weight(0.05), kHalf, design::StateLimits{.states = {1}, .lim = {.lower = -0.3, .upper = 0.3}});
    REQUIRE(plain.success);
    REQUIRE(bounded.success);
    MPC a{plain};
    MPC b{bounded};
    static_assert(std::is_empty_v<decltype(a.y_lim)> && std::is_empty_v<decltype(a.x_lim)>);
    static_assert(!reads_output_limits<decltype(a)> && !reads_state_limits<decltype(a)>);
    static_assert(std::is_empty_v<decltype(b.y_lim)> && reads_state_limits<decltype(b)>);
    CHECK(b.x_lim.upper(0) == 0.3);

    Law<kOnConstants, double> on_constants{};
    static_assert(reads_output_limits<decltype(on_constants)> && !reads_state_limits<decltype(on_constants)>);
    static_assert(std::is_empty_v<decltype(on_constants.x_lim)>);
}

TEST_CASE("an LQ law's u_lim and an MPC's are assigned the same way, and both bind from the next tick") {
    const auto lq = design::state_feedback(double_integrator, design::LQRWeights{.Q = Matrix<2, 2>::identity(), .R = {0.01}});
    const auto mpc = design::mpc<20, 8>(double_integrator, move_weight(0.05));
    REQUIRE(lq.success);
    REQUIRE(mpc.success);
    LQR             lqr{lq};
    MPC             law{mpc};
    const ColVec<2> far{-1.0, 0.0};
    lqr.update(far);
    law.update(far);
    REQUIRE(lqr.control()(0) > 0.1);
    REQUIRE(law.control(ColVec<1>{0.0})(0) > 0.1);

    lqr.u_lim = {.lower = -0.1, .upper = 0.1};
    law.u_lim = {.lower = -0.1, .upper = 0.1};
    CHECK(lqr.control()(0) == 0.1);
    CHECK(law.control(ColVec<1>{0.0})(0) == doctest::Approx(0.1).epsilon(1e-12));
    CHECK(law.last_status() == design::QPStatus::Success);
}
