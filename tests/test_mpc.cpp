// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/controllers/state_feedback_types.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/input_limits.hpp"
#include "damp/design/qp.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/fixed_gain_observer.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

constexpr double Ts = 0.1;

// Double integrator, ZOH-discretized at Ts = 0.1, the state measured, the position the controlled output.
constexpr StateSpace double_integrator{
    .A = {{1.0, Ts}, {0.0, 1.0}},
    .B = {{0.5 * Ts * Ts}, {Ts}},
    .C = {{1.0, 0.0}},
    .Ts = Ts,
};

// Same plant with every state a controlled output (C = I), for the LQR-equivalence tests.
constexpr StateSpace double_integrator_full{
    .A = {{1.0, Ts}, {0.0, 1.0}},
    .B = {{0.5 * Ts * Ts}, {Ts}},
    .C = {{1.0, 0.0}, {0.0, 1.0}},
    .Ts = Ts,
};

/// The move weight @p rdu, the position weighted by the default Qy = 1.
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

/// A solver policy that keeps the QP of the last tick and its result: the data for a check of the
/// solve that owes nothing to the solver.
template<size_t NV, size_t NI, typename T>
struct RecordingSolver {
    design::WarmStartActiveSetSolver<NV, NI, T> solver{};
    Matrix<NV, NV, T>                           H{};
    ColVec<NV, T>                               f{};
    Matrix<NI, NV, T>                           A{};
    ColVec<NI, T>                               b{};
    design::QPResult<NV, NI, T>                 result{};

    design::QPResult<NV, NI, T> operator()(
        const Matrix<NV, NV, T>& H_in, const ColVec<NV, T>& f_in, const Matrix<NI, NV, T>& A_in, const ColVec<NI, T>& b_in, size_t max_iterations
    ) {
        H = H_in;
        f = f_in;
        A = A_in;
        b = b_in;
        result = solver(H_in, f_in, A_in, b_in, max_iterations);
        return result;
    }
};

/// The KKT conditions of a recorded solve, formed in double from its data. Each residual is over
/// the terms that form it, so it reads against the ε of the solve's own scalar.
struct KktResiduals {
    double primal{0.0};       ///< Largest aᵢᵀx − bᵢ, over |bᵢ| + Σ|aᵢⱼxⱼ| + ‖aᵢ‖₁‖x‖∞
    double tight{0.0};        ///< Largest |aᵢᵀx − bᵢ| on a row with λᵢ > 0, over the same
    double stationarity{0.0}; ///< ‖Hx + f + Aᵀλ‖∞ over the largest sum of one entry's terms
    double lambda_min{0.0};   ///< Smallest multiplier, 0 when none is negative

    void take_worst(const KktResiduals& other) {
        primal = damp::max(primal, other.primal);
        tight = damp::max(tight, other.tight);
        stationarity = damp::max(stationarity, other.stationarity);
        lambda_min = damp::min(lambda_min, other.lambda_min);
    }
};

template<size_t NV, size_t NI, typename T>
KktResiduals kkt_residuals(const RecordingSolver<NV, NI, T>& qp) {
    const auto                         at = [](T v) { return static_cast<double>(v); };
    const double                       disabled = at(design::unbounded_bound<T>()) / 4.0;
    const design::QPResult<NV, NI, T>& res = qp.result;
    KktResiduals                       out{};
    double                             x_max = 0.0;
    for (size_t v = 0; v < NV; ++v) {
        x_max = damp::max(x_max, damp::abs(at(res.x(v))));
    }
    for (size_t i = 0; i < NI; ++i) {
        if (at(qp.b(i)) >= disabled) {
            continue;
        }
        double ax = 0.0;
        double terms = damp::abs(at(qp.b(i)));
        double a_norm = 0.0;
        for (size_t v = 0; v < NV; ++v) {
            const double term = at(qp.A(i, v)) * at(res.x(v));
            ax += term;
            terms += damp::abs(term);
            a_norm += damp::abs(at(qp.A(i, v)));
        }
        const double scale = terms + (a_norm * x_max);
        const double miss = ax - at(qp.b(i));
        if (scale > 0.0) {
            out.primal = damp::max(out.primal, miss / scale);
            if (res.lambda(i) > T{0}) {
                out.tight = damp::max(out.tight, damp::abs(miss) / scale);
            }
        }
        out.lambda_min = damp::min(out.lambda_min, at(res.lambda(i)));
    }
    double gradient = 0.0;
    double largest = 0.0;
    for (size_t v = 0; v < NV; ++v) {
        double g = at(qp.f(v));
        double terms = damp::abs(g);
        for (size_t c = 0; c < NV; ++c) {
            const double term = at(qp.H(v, c)) * at(res.x(c));
            g += term;
            terms += damp::abs(term);
        }
        for (size_t i = 0; i < NI; ++i) {
            const double term = at(qp.A(i, v)) * at(res.lambda(i));
            g += term;
            terms += damp::abs(term);
        }
        gradient = damp::max(gradient, damp::abs(g));
        largest = damp::max(largest, terms);
    }
    //! Below √min the products of T underflow: a loop at rest (x ≈ 1e-33 in float) has no digits to test.
    if (largest > damp::sqrt(at(damp::numeric_limits<T>::min()))) {
        out.stationarity = gradient / largest;
    }
    return out;
}

/// The soft-bounded loop of the float regression: |u| ≤ 1 (soft_box), |Δu| ≤ 0.2 (soft_moves), soft |y| ≤ 0.5 (soft_y), r = 0.
constexpr design::InputLimits<1>     soft_box{.lower = {-1.0}, .upper = {1.0}};
constexpr design::InputRateLimits<1> soft_moves{.lower = {-0.2}, .upper = {0.2}};

/// Its soft |y| ≤ 0.5 on the position.
constexpr design::OutputLimits soft_y{.outputs = {0}, .lower = {-0.5}, .upper = {0.5}};

/// Whether @p got is the type @p Want: what a limit option deduced.
template<typename Want, typename Got>
constexpr bool deduces_to(const Got& got) {
    static_cast<void>(got);
    return std::is_same_v<Want, Got>;
}

/// The limit options deduce inside a function template as at namespace scope: a braced end deduces
/// nothing (one input, double), a typed end and an index list do. On clang a ColVec literal set as the
/// end itself deduces nothing in a template (llvm/llvm-project#114870), so the typed ends are locals.
template<size_t N>
constexpr bool deduced_in_template() {
    const ColVec lower{-1.0, -2.0};
    const ColVec upper{1.0f, 2.0f};
    return deduces_to<design::InputLimits<1>>(design::InputLimits{.upper = {1.0}})
        && deduces_to<design::InputLimits<2>>(design::InputLimits{.lower = lower, .upper = {1.0, 2.0}})
        && deduces_to<design::InputRateLimits<1>>(design::InputRateLimits{.lower = {-0.2}, .upper = {0.2}})
        && deduces_to<design::InputRateLimits<2, float>>(design::InputRateLimits{.upper = upper})
        && deduces_to<design::OutputLimits<2>>(design::OutputLimits{.outputs = {1, 2}, .lower = {-0.4, -0.26}, .upper = {0.4, 0.26}})
        && deduces_to<design::StateLimits<1>>(design::StateLimits{.states = {2}, .lower = {-0.26}, .upper = {0.26}})
        && deduces_to<DemandLimit<3>>(DemandLimit{.state = 2, .outer = {0, 1, 4}, .lower = -0.26, .upper = 0.26}) && N == 0;
}

} // namespace

TEST_CASE("unconstrained MPC with a DARE terminal weight reproduces the LQR law") {
    // Finite-horizon LQ with terminal cost S (the DARE solution) gives
    // u0 = -K x for any horizon; the MPC must match move for move.
    const Matrix Q{{1.0, 0.0}, {0.0, 0.1}};
    const Matrix R{{0.5}};

    const auto lqr = design::discrete_lqr(double_integrator_full.A, double_integrator_full.B, Q, R);
    REQUIRE(lqr.success);

    design::MPCWeights<1, 2> weights{};
    weights.Qy = Q;
    weights.Qy_terminal = lqr.S;
    weights.Rdu = Matrix<1, 1>{}; // no move penalty: pure LQ cost
    weights.Ru = R;

    const auto art = design::mpc<8, 8>(double_integrator_full, weights); // a free move every step: equals LQR; no TrackedOutputs, a regulator
    REQUIRE(art.success);

    MPC    controller{art};
    ColVec x{1.0, -0.5};

    for (size_t k = 0; k < 30; ++k) {
        const auto   u_mpc = controller.step(x);
        const double u_lqr = -((lqr.K(0, 0) * x(0)) + (lqr.K(0, 1) * x(1)));
        CHECK(u_mpc(0) == doctest::Approx(u_lqr).epsilon(1e-6));
        CHECK(controller.last_status() == design::QPStatus::Success);
        CHECK(controller.last_slack() == doctest::Approx(0.0));
        x = plant_step(x, u_mpc(0));
    }
}

TEST_CASE("constrained MPC tracks a step within u and du limits") {
    auto weights = move_weight(0.05);

    constexpr design::InputLimits<1>     limits{.lower = {-0.5}, .upper = {0.5}};
    constexpr design::InputRateLimits<1> moves{.lower = {-0.2}, .upper = {0.2}};

    const auto art = design::mpc<20, 8>(double_integrator, weights, limits, moves, design::TrackedOutputs{});
    REQUIRE(art.success);

    MPC          controller{art};
    ColVec<2>    x{};
    const ColVec r{1.0};

    double u_last = 0.0;
    double u_peak = 0.0;
    for (size_t k = 0; k < 200; ++k) {
        const auto u = controller.step(r, x);
        CHECK(controller.last_status() == design::QPStatus::Success);
        CHECK(controller.last_iterations() <= art.mpc.max_qp_iterations);

        // Hard limits hold at every tick.
        CHECK(u(0) <= 0.5 + 1e-9);
        CHECK(u(0) >= -0.5 - 1e-9);
        CHECK(u(0) - u_last <= 0.2 + 1e-9);
        CHECK(u(0) - u_last >= -0.2 - 1e-9);

        u_peak = (damp::abs(u(0)) > u_peak) ? damp::abs(u(0)) : u_peak;
        u_last = u(0);
        x = plant_step(x, u(0));
    }

    // The step was aggressive enough that the input limit actually bound...
    CHECK(u_peak == doctest::Approx(0.5));
    // ...and the loop still converged to the reference.
    CHECK(x(0) == doctest::Approx(1.0).epsilon(0.02));
    CHECK(x(1) == doctest::Approx(0.0).epsilon(0.02));
}

TEST_CASE("InputLimits and InputRateLimits are the QP and the commands MPCConstraints gave; left out, the inputs are unbounded") {
    //! Recorded with the positional MPCConstraints<1>{.u_lim, .du_lim} before InputLimits replaced
    //! it, and again with InputLimits{.u_lim, .du_lim} before InputRateLimits took the move limit.
    const design::InputLimits     limits{.lower = {-0.5}, .upper = {0.5}};
    const design::InputRateLimits moves{.lower = {-0.2}, .upper = {0.2}};
    const auto                    w = move_weight(0.05);
    const auto                    art = design::mpc<10, 4>(double_integrator, w, limits, moves, design::TrackedOutputs{});
    REQUIRE(art.success);
    const auto& a = art.mpc;
    CHECK(a.H(0, 0) == doctest::Approx(1.3666500000000006).epsilon(1e-12));
    CHECK(a.H(0, 1) == doctest::Approx(0.9834000000000005).epsilon(1e-12));
    CHECK(a.H(3, 3) == doctest::Approx(0.3338000000000001).epsilon(1e-12));
    CHECK(a.F_x(0, 0) == doctest::Approx(3.8500000000000005).epsilon(1e-12));
    CHECK(a.F_x(0, 1) == doctest::Approx(3.0250000000000004).epsilon(1e-12));
    CHECK(a.F_x(0, 2) == doctest::Approx(1.2666500000000005).epsilon(1e-12));
    //! The first rows of Δu (0, NZ) and of u (2·NZ, 3·NZ), NZ = 4, on the input's scale, its range 1.
    CHECK(a.scale_u(0) == 1.0);
    CHECK(a.b_bound(0) == doctest::Approx(0.2).epsilon(1e-15));
    CHECK(a.b_bound(4) == doctest::Approx(0.2).epsilon(1e-15));
    CHECK(a.b_bound(8) == doctest::Approx(0.5).epsilon(1e-15));
    CHECK(a.b_bound(12) == doctest::Approx(0.5).epsilon(1e-15));
    CHECK(a.B_map(8, 2) == -1.0);
    CHECK(a.A_con(8, 0) == 1.0);
    CHECK(a.A_con(9, 0) == 1.0);
    CHECK(a.max_qp_iterations == 220);
    //! The commands at ticks 0, 1, 2, 5, 10, 20 and 39 of a step, then the state at the end.
    const auto run = [](const auto& res) {
        MPC                    law{res};
        ColVec<2>              x{};
        damp::array<double, 9> out{};
        size_t                 next = 0;
        for (int k = 0; k < 40; ++k) {
            const auto u = law.step(ColVec{1.0}, x);
            if (k == 0 || k == 1 || k == 2 || k == 5 || k == 10 || k == 20 || k == 39) {
                out[next++] = u(0);
            }
            x = plant_step(x, u(0));
        }
        out[next++] = x(0);
        out[next] = x(1);
        return out;
    };
    const auto got = run(art);
    CHECK(got[0] == doctest::Approx(0.2).epsilon(1e-12)); // the move limit
    CHECK(got[1] == doctest::Approx(0.4).epsilon(1e-12));
    CHECK(got[2] == doctest::Approx(0.5).epsilon(1e-12)); // the input limit
    CHECK(got[3] == doctest::Approx(0.5).epsilon(1e-12));
    CHECK(got[4] == doctest::Approx(0.5).epsilon(1e-12));
    CHECK(got[5] == doctest::Approx(-0.5).epsilon(1e-12));
    CHECK(got[6] == doctest::Approx(-0.30805135702196029).epsilon(1e-8));
    CHECK(got[7] == doctest::Approx(1.1791913358321511).epsilon(1e-8));
    CHECK(got[8] == doctest::Approx(-0.32435392995393764).epsilon(1e-8));

    //! The same design with the limits after another option.
    const auto same = [&](const auto& other) {
        REQUIRE(other.success);
        const auto& b = other.mpc;
        CHECK((b.H - a.H).norm() <= 4.0 * eps<double>() * a.H.norm());
        CHECK((b.F_x - a.F_x).norm() <= 4.0 * eps<double>() * a.F_x.norm());
        CHECK((b.A_con - a.A_con).norm() <= 4.0 * eps<double>() * a.A_con.norm());
        CHECK((b.B_map - a.B_map).norm() <= 4.0 * eps<double>() * a.B_map.norm());
        for (size_t i = 0; i < 16; ++i) {
            CHECK(b.b_bound(i) == doctest::Approx(a.b_bound(i)).epsilon(1e-15));
        }
        const auto commands = run(other);
        for (size_t i = 0; i < commands.size(); ++i) {
            CHECK(commands[i] == doctest::Approx(got[i]).epsilon(1e-9));
        }
    };
    same(design::mpc<10, 4>(double_integrator, w, design::QPIterations{220}, moves, limits, design::TrackedOutputs{}));

    //! Limits in another scalar are read in the plant's.
    const auto in_float = design::mpc<10, 4>(double_integrator, w, design::InputLimits<1, float>{.lower = {-0.5f}, .upper = {0.5f}}, design::InputRateLimits<1, float>{.lower = {-0.25f}, .upper = {0.25f}}, design::TrackedOutputs{});
    REQUIRE(in_float.success);
    CHECK(in_float.mpc.b_bound(0) == 0.25);
    CHECK(in_float.mpc.b_bound(8) == 0.5);

    //! Left out, u is unbounded, sentinel rows, and Δu has no rows: the cost of the bounded design (its scale is 1 too).
    const auto free = design::mpc<10, 4>(double_integrator, w, design::TrackedOutputs{});
    REQUIRE(free.success);
    static_assert(std::is_same_v<decltype(free.mpc.constraints.limits.input_rate), Bounds<0>> && decltype(free.mpc)::NI == 2 * 4);
    CHECK(free.mpc.constraints.limits.input.upper(0) == damp::numeric_limits<double>::max());
    for (size_t i = 0; i < 8; ++i) {
        CHECK(free.mpc.b_bound(i) >= design::unbounded_bound<double>() / 4.0);
    }
    CHECK((free.mpc.H - a.H).norm() <= 4.0 * eps<double>() * a.H.norm());
    const auto spelled = design::mpc<10, 4>(double_integrator, w, design::InputLimits<1>{}, design::InputRateLimits<1>{}, design::TrackedOutputs{});
    REQUIRE(spelled.success);
    for (size_t i = 0; i < 16; ++i) {
        CHECK(spelled.mpc.b_bound(i) >= design::unbounded_bound<double>() / 4.0);
    }
}

TEST_CASE("an OFMPC on the sampled form with InputLimits is the QP and the commands MPCConstraints gave") {
    //! Recorded with MPCConstraints<1>{.u_lim = ±1.05, .du_lim = ±0.1} before InputLimits replaced it: a
    //! mass on a spring sampled at 0.05 s, its reading filtered, offset-free under a load at the input.
    constexpr StateSpace<2, 1, 1, double, 2, 1> spring{
        .A = Matrix{{0.0, 1.0}, {-1.0, -0.5}},
        .B = Matrix{{0.0}, {1.0}},
        .C = Matrix{{1.0, 0.0}},
        .G = Matrix<2, 2>::identity(),
        .H = Matrix<1, 1>::identity(),
    };
    const design::KalmanNoise<2, 1, double, 1> noise{.process = Matrix{{0.01, 0.0}, {0.0, 0.1}}, .measurement = Matrix{{1e-4}}, .disturbance = Matrix{{1.0}}};
    const design::InputLimits                  limits{.lower = {-1.05}, .upper = {1.05}};
    const design::InputRateLimits              moves{.lower = {-0.1}, .upper = {0.1}};
    const auto                                 of = design::mpc<20, 5>(spring, 0.05, move_weight(0.1), limits, moves, noise, design::InputDisturbance{});
    REQUIRE(of.success);
    CHECK(of.mpc.H(0, 0) == doctest::Approx(4.3012663890329419).epsilon(1e-10));
    CHECK(of.mpc.b_bound(0) == doctest::Approx(0.047619047619047616).epsilon(1e-14)); // 0.1 over the range 2.1
    CHECK(of.mpc.b_bound(10) == doctest::Approx(0.5).epsilon(1e-14));
    const auto d = discretize(spring, 0.05).value();
    OFMPC      law{of};
    ColVec<2>  z{};
    for (int k = 0; k < 120; ++k) {
        const auto u = law.step(ColVec{1.0}, ColVec<1>{z(0)});
        if (k == 0) {
            CHECK(u(0) == doctest::Approx(0.1).epsilon(1e-12));
        }
        if (k == 1) {
            CHECK(u(0) == doctest::Approx(0.2).epsilon(1e-12));
        }
        if (k == 5) {
            CHECK(u(0) == doctest::Approx(0.6).epsilon(1e-12));
        }
        if (k == 20) {
            CHECK(u(0) == doctest::Approx(1.0430532763888376).epsilon(1e-8));
        }
        if (k == 119) {
            CHECK(u(0) == doctest::Approx(0.39438686032214743).epsilon(1e-8));
        }
        z = (d.A * z) + (d.B * (u + ColVec{0.3}));
    }
    CHECK(z(0) == doctest::Approx(1.0915662238798631).epsilon(1e-8));
}

TEST_CASE("set_limits changes the input limit at run time; a bad limit is refused") {
    design::InputLimits<1> limits{.lower = {-0.5}, .upper = {0.5}};
    const auto             art = design::mpc<20, 8>(double_integrator, move_weight(0.05), limits, design::TrackedOutputs{});
    REQUIRE(art.success);

    MPC  controller{art};
    auto lim = controller.limits();
    CHECK(lim.input.upper(0) == 0.5);
    static_assert(std::is_same_v<decltype(lim.input_rate), Bounds<0>>); // no InputRateLimits: no move limit to set

    lim.input.upper(0) = 0.2;
    REQUIRE(controller.set_limits(lim));
    CHECK(controller.limits().input.upper(0) == 0.2);

    ColVec<2>    x{};
    const ColVec r{1.0};
    double       u_peak = 0.0;
    for (size_t k = 0; k < 400; ++k) {
        const auto u = controller.step(r, x);
        REQUIRE(controller.last_status() == design::QPStatus::Success);
        CHECK(u(0) <= 0.2 + 1e-9); // the new limit, not the design's 0.5
        u_peak = (u(0) > u_peak) ? u(0) : u_peak;
        x = plant_step(x, u(0));
    }
    CHECK(u_peak == doctest::Approx(0.2));
    CHECK(x(0) == doctest::Approx(1.0).epsilon(0.02));

    auto bad = controller.limits();
    bad.input.lower(0) = 0.3; // above u.upper
    CHECK_FALSE(controller.set_limits(bad));
    CHECK(controller.limits().input.lower(0) == -0.5); // unchanged
}

TEST_CASE("hard output constraint caps the predicted response") {
    auto weights = move_weight(0.01);

    design::InputLimits<1> limits{.lower = {-1.0}, .upper = {1.0}};
    // An overshoot limiter on the position, hard: no slack allowed.
    const design::OutputLimits cap{.outputs = {0}, .upper = {1.02}, .upper_ecr = {0.0}};

    const auto art = design::mpc<20, 8>(double_integrator, weights, limits, cap, design::TrackedOutputs{});
    REQUIRE(art.success);

    MPC          controller{art};
    ColVec<2>    x{};
    const ColVec r{1.0};

    for (size_t k = 0; k < 300; ++k) {
        const auto u = controller.step(r, x);
        x = plant_step(x, u(0));
        CHECK(x(0) <= 1.02 + 1e-6); // plant == model, so the predicted cap holds exactly
    }
    CHECK(x(0) == doctest::Approx(1.0).epsilon(0.02));
}

TEST_CASE_TEMPLATE("StateLimits bound a state the MPC does not track, hard or soft, and move at run time", T, double, float) {
    // The position is tracked, the speed only bounded: a step of 1 at full input would run past 0.3.
    const design::InputLimits<1> limits{.lower = {-1.0}, .upper = {1.0}};
    design::StateLimits          speed{.states = {1}, .lower = {-0.3}, .upper = {0.3}};
    speed.upper_ecr(0) = 0.0; // hard above
    const auto art = design::mpc<20, 8>(double_integrator, move_weight(0.01), limits, speed, design::TrackedOutputs{});
    REQUIRE(art.success);
    const auto free_art = design::mpc<20, 8>(double_integrator, move_weight(0.01), limits, design::TrackedOutputs{});
    REQUIRE(free_art.success);

    MPC        bounded{art.template as<T>()};
    MPC        free{free_art.template as<T>()};
    const auto run = [](auto& law, double& peak) {
        ColVec<2> x{};
        peak = 0.0;
        for (size_t k = 0; k < 300; ++k) {
            const auto u = law.step(ColVec<1, T>{T{1}}, x.template as<T>());
            CHECK(law.last_status() == design::QPStatus::Success);
            x = plant_step(x, static_cast<double>(u(0)));
            peak = std::fmax(peak, x(1));
        }
        return x;
    };
    double peak_free = 0.0;
    double peak = 0.0;
    static_cast<void>(run(free, peak_free));
    const auto x = run(bounded, peak);
    CHECK(peak_free > 0.4);                                        // without the bound the speed runs past it
    CHECK(peak <= 0.3 + (std::is_same_v<T, float> ? 1e-5 : 1e-9)); // hard: plant == model, so the predicted cap holds
    CHECK(x(0) == doctest::Approx(1.0).epsilon(0.02));             // and the position still settles

    auto lim = bounded.limits();
    CHECK(lim.state.upper(0) == T(0.3));
    lim.state.upper(0) = T(0.15);
    REQUIRE(bounded.set_limits(lim));
    bounded.reset();
    static_cast<void>(run(bounded, peak));
    CHECK(peak <= 0.15 + (std::is_same_v<T, float> ? 1e-5 : 1e-9));
    lim.state.lower(0) = T(1.0); // above its upper end
    CHECK_FALSE(bounded.set_limits(lim));

    // Soft (the default): started past the bound, every tick solves and the speed comes back inside.
    const auto soft_art = design::mpc<20, 8>(double_integrator, move_weight(0.01), limits, design::StateLimits{.states = {1}, .lower = {-0.3}, .upper = {0.3}}, design::TrackedOutputs{});
    REQUIRE(soft_art.success);
    MPC    recovering{soft_art.template as<T>()};
    ColVec xs{0.0, 0.8};
    for (size_t k = 0; k < 100; ++k) {
        const auto u = recovering.step(ColVec<1, T>{T{0}}, xs.template as<T>());
        CHECK(recovering.last_status() == design::QPStatus::Success);
        xs = plant_step(xs, static_cast<double>(u(0)));
    }
    CHECK(std::fabs(xs(1)) <= 0.3);
}

TEST_CASE("StateLimits refusals: a state past the model's, one listed twice, a bad box, a negative softness") {
    static_assert(damp::is_same_v<decltype(design::StateLimits{.states = {0, 2}}), design::StateLimits<2>>); // NXC from the list
    CHECK(design::mpc<10, 4>(double_integrator, move_weight(0.05), design::StateLimits{.states = {2}}, design::TrackedOutputs{}).status == design::DesignStatus::invalid_parameter);
    CHECK(design::mpc<10, 4>(double_integrator, move_weight(0.05), design::StateLimits{.states = {1, 1}}, design::TrackedOutputs{}).status == design::DesignStatus::invalid_parameter);
    CHECK(design::mpc<10, 4>(double_integrator, move_weight(0.05), design::StateLimits{.states = {1}, .lower = {1.0}, .upper = {-1.0}}, design::TrackedOutputs{}).status == design::DesignStatus::invalid_parameter);
    CHECK(design::mpc<10, 4>(double_integrator, move_weight(0.05), design::StateLimits{.states = {1}, .lower_ecr = ColVec{-1.0}}, design::TrackedOutputs{}).status == design::DesignStatus::bad_weights);
}

TEST_CASE("OutputLimits on two of four outputs is the QP of every output listed, the others at ±max, in fewer rows") {
    //! Two double integrators, every state a controlled output, the positions tracked: the first
    //! axis's speed bounded softly, the second's position hard from above. Before OutputLimits the
    //! limits were MPCConstraints' y_lim on all four, the free ones at ±max.
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
    const design::InputLimits<2> limits{.lower = {-1.0, -1.0}, .upper = {1.0, 1.0}};
    design::OutputLimits         subset{.outputs = {1, 2}, .lower = {-0.3, -big}, .upper = {0.3, 0.4}};
    subset.upper_ecr(1) = 0.0;
    design::OutputLimits every{.outputs = {0, 1, 2, 3}, .lower = {-big, -0.3, -big, -big}, .upper = {big, 0.3, 0.4, big}};
    every.upper_ecr(2) = 0.0;

    const design::TrackedOutputs positions{{0, 2}}; // r on the positions, one per input
    const auto                   a = design::mpc<12, 4>(two_axes, w, limits, subset, positions);
    const auto                   b = design::mpc<12, 4>(two_axes, w, limits, every, positions);
    REQUIRE(a.success);
    REQUIRE(b.success);
    using A = std::remove_cvref_t<decltype(a.mpc)>;
    using B = std::remove_cvref_t<decltype(b.mpc)>;
    static_assert(A::NI == (2 * 8) + (2 * 12 * 2) && B::NI == (2 * 8) + (2 * 12 * 4)); // u's rows, then 2·NP rows per bounded output
    static_assert(std::remove_cvref_t<decltype(design::mpc<12, 4>(two_axes, w, limits, positions).mpc)>::NI == 2 * 8);
    CHECK(B::NI == 112); // the output rows of every output
    for (size_t i = 0; i < 4; ++i) {
        CHECK(a.mpc.scale_y(i) == b.mpc.scale_y(i));
    }
    CHECK(a.mpc.scale_y(1) == doctest::Approx(0.6)); // a listed output with both sides finite: its range
    for (size_t i = 0; i < A::NZ1; ++i) {
        for (size_t j = 0; j < A::NZ1; ++j) {
            CHECK(a.mpc.H(i, j) == doctest::Approx(b.mpc.H(i, j)).epsilon(1e-14));
        }
    }
    //! The rows are those of the listed outputs, step by step, in the order listed.
    for (size_t kb = 0; kb < 12; ++kb) {
        for (size_t j = 0; j < 2; ++j) {
            for (size_t c = 0; c < A::NZ1; ++c) {
                CHECK(a.mpc.A_con((2 * 8) + (kb * 2) + j, c) == b.mpc.A_con((2 * 8) + (kb * 4) + 1 + j, c));
            }
            CHECK(a.mpc.b_bound((2 * 8) + (kb * 2) + j) == b.mpc.b_bound((2 * 8) + (kb * 4) + 1 + j));
        }
    }

    //! The same commands, and the ones y_lim gave before the option.
    MPC          la{a};
    MPC          lb{b};
    ColVec<4>    x{};
    const ColVec r{1.0, 0.5}; // the speeds at the target's 0
    double       v_peak = 0.0;
    double       p_peak = 0.0;
    for (size_t k = 0; k < 60; ++k) {
        const auto ua = la.step(r, x);
        const auto ub = lb.step(r, x);
        REQUIRE(la.last_status() == design::QPStatus::Success);
        CHECK(ua(0) == doctest::Approx(ub(0)).epsilon(1e-12));
        CHECK(ua(1) == doctest::Approx(ub(1)).epsilon(1e-12));
        if (k == 3) {
            CHECK(ua(0) == doctest::Approx(0.00019054062607748445).epsilon(1e-8));
            CHECK(ua(1) == doctest::Approx(1.0).epsilon(1e-12));
        }
        if (k == 6) {
            CHECK(ua(1) == doctest::Approx(-0.35740042951566808).epsilon(1e-8));
        }
        if (k == 10) {
            CHECK(ua(1) == doctest::Approx(-0.71598701148296784).epsilon(1e-8));
        }
        x = (two_axes.A * x) + (two_axes.B * ua);
        v_peak = damp::max(v_peak, x(1));
        p_peak = damp::max(p_peak, x(2));
    }
    CHECK(v_peak <= 0.3 + 1e-3); // soft: about the bound
    CHECK(p_peak <= 0.4 + 1e-9); // hard: plant == model

    //! At run time the limits are the listed outputs', in the order listed.
    auto lim = la.limits();
    static_assert(std::is_same_v<decltype(lim.output), Bounds<2>>);
    CHECK(lim.output.upper(1) == 0.4);
    lim.output.upper(1) = 0.35;
    REQUIRE(la.set_limits(lim));
    CHECK(la.limits().output.upper(1) == 0.35);
}

TEST_CASE("OutputLimits refusals: an output past the controlled ones, one listed twice, a bad box, a negative softness") {
    CHECK(design::mpc<10, 4>(double_integrator, move_weight(0.05), design::OutputLimits{.outputs = {1}}, design::TrackedOutputs{}).status == design::DesignStatus::invalid_parameter);
    CHECK(design::mpc<10, 4>(double_integrator_full, design::MPCWeights<1, 2>{}, design::TrackedOutputs{{0}}, design::OutputLimits{.outputs = {1, 1}}).status == design::DesignStatus::invalid_parameter);
    CHECK(design::mpc<10, 4>(double_integrator, move_weight(0.05), design::OutputLimits{.outputs = {0}, .lower = {1.0}, .upper = {-1.0}}, design::TrackedOutputs{}).status == design::DesignStatus::invalid_parameter);
    CHECK(design::mpc<10, 4>(double_integrator, move_weight(0.05), design::OutputLimits{.outputs = {0}, .lower_ecr = {-1.0}}, design::TrackedOutputs{}).status == design::DesignStatus::bad_weights);
    CHECK(design::mpc<10, 4>(double_integrator_full, design::MPCWeights<1, 2>{}, design::TrackedOutputs{{0}}, design::OutputLimits{.outputs = {1}}).success); // the speed, output 1 of C = I
}

TEST_CASE("QPIterations sets the QP's budget; left out, the default 10·(NZ + 1 + NI + 1)") {
    const auto set = design::mpc<10, 4>(double_integrator, move_weight(0.05), design::QPIterations{40}, design::TrackedOutputs{});
    const auto left_out = design::mpc<10, 4>(double_integrator, move_weight(0.05), design::TrackedOutputs{});
    REQUIRE(set.success);
    REQUIRE(left_out.success);
    CHECK(set.mpc.max_qp_iterations == 40);
    CHECK(left_out.mpc.max_qp_iterations == 10 * (5 + 8 + 1));                                                                                                            // no Δu rows without an InputRateLimits, no output rows without an OutputLimits
    CHECK(design::mpc<10, 4>(double_integrator, move_weight(0.05), design::InputRateLimits<1>{}, design::TrackedOutputs{}).mpc.max_qp_iterations == 10 * (5 + 16 + 1));   // with the Δu rows
    CHECK(design::mpc<10, 4>(double_integrator, move_weight(0.05), design::OutputLimits{.outputs = {0}}, design::TrackedOutputs{}).mpc.max_qp_iterations == 340);         // with the output rows
    CHECK(design::mpc<10, 4>(double_integrator, move_weight(0.05), design::QPIterations{0}, design::TrackedOutputs{}).status == design::DesignStatus::invalid_parameter); // not one active-set step

    //! The runtime solves within it: a budget of one leaves a constrained tick MaxIterations.
    MPC  capped{design::mpc<10, 4>(double_integrator, move_weight(0.05), design::InputLimits{.lower = {-0.1}, .upper = {0.1}}, design::QPIterations{1}, design::TrackedOutputs{})};
    auto u = capped.step(1.0, ColVec<2>{});
    CHECK(capped.last_iterations() <= 1);
    CHECK(capped.last_status() == design::QPStatus::MaxIterations);
    CHECK(u(0) <= 0.1);
}

TEST_CASE("soft output constraint recovers from an infeasible start") {
    // Start above the output ceiling with limited authority: no admissible
    // input sequence satisfies the hard rows, so a hard configuration reports
    // Infeasible while the default (soft) configuration keeps solving.
    auto weights = move_weight(0.01);

    design::InputLimits<1> limits{.lower = {-1.0}, .upper = {1.0}};
    design::OutputLimits   ceiling{.outputs = {0}, .upper = {1.0}};

    const ColVec x0{2.0, 0.0}; // x₀ = 2 > y.upper
    const ColVec r{0.5};

    // Hard variant: first solve is infeasible and the fallback holds the input.
    {
        design::OutputLimits hard = ceiling;
        hard.upper_ecr(0) = 0.0;
        const auto art = design::mpc<20, 8>(double_integrator, weights, limits, hard, design::TrackedOutputs{});
        REQUIRE(art.success);
        MPC        controller{art};
        const auto u = controller.step(r, x0);
        CHECK(controller.last_status() == design::QPStatus::Infeasible);
        CHECK(u(0) == doctest::Approx(0.0)); // Δu falls back to zero
    }

    // Soft variant (default ECR = 1): solvable from the violated start, slack
    // active initially, and the loop still converges to the reference.
    {
        const auto art = design::mpc<20, 8>(double_integrator, weights, limits, ceiling, design::TrackedOutputs{});
        REQUIRE(art.success);
        MPC controller{art};

        ColVec<2>  x = x0;
        const auto u0 = controller.step(r, x);
        CHECK(controller.last_status() == design::QPStatus::Success);
        CHECK(controller.last_slack() > 0.0);
        x = plant_step(x, u0(0));

        for (size_t k = 0; k < 400; ++k) {
            const auto u = controller.step(r, x);
            CHECK(controller.last_status() == design::QPStatus::Success);
            x = plant_step(x, u(0));
        }
        CHECK(x(0) == doctest::Approx(0.5).epsilon(0.02));
        CHECK(controller.last_slack() == doctest::Approx(0.0).epsilon(1e-6));
    }
}

TEST_CASE("an infeasible tick holds the last applied input inside the u box") {
    design::InputLimits<1> hard{.lower = {-1.0}, .upper = {1.0}};
    const auto             art = design::mpc<20, 8>(double_integrator, move_weight(0.01), hard, design::OutputLimits{.outputs = {0}, .upper = {1.0}, .upper_ecr = {0.0}}, design::TrackedOutputs{});
    REQUIRE(art.success);
    MPC          controller{art};
    ColVec<2>    x{};
    const ColVec r{0.5};
    double       prev = 0.0;
    for (size_t k = 0; k < 20; ++k) {
        const auto u = controller.step(r, x);
        CHECK(u(0) >= -1.0);
        CHECK(u(0) <= 1.0);
        REQUIRE(controller.last_status() == design::QPStatus::Success);
        prev = u(0);
        x = plant_step(x, u(0));
    }
    x = ColVec{2.0, 0.0};
    const auto held = controller.step(r, x);
    CHECK(controller.last_status() == design::QPStatus::Infeasible);
    CHECK(held(0) == doctest::Approx(prev));
    CHECK(held(0) >= -1.0);
    CHECK(held(0) <= 1.0);
}

TEST_CASE("MaxIterations with hard output rows holds unless the iterate is primal feasible") {
    design::InputLimits<1>     limits{.lower = {-1.0}, .upper = {1.0}};
    const design::OutputLimits cap{.outputs = {0}, .upper = {1.02}, .upper_ecr = {0.0}};
    const auto                 art = design::mpc<30, 10>(double_integrator, move_weight(0.01), limits, cap, design::QPIterations{1}, design::TrackedOutputs{});
    REQUIRE(art.success);
    MPC          controller{art};
    ColVec<2>    x{};
    const ColVec r{1.0};
    const auto   u0 = controller.step(r, x);
    CHECK(u0(0) >= -1.0);
    CHECK(u0(0) <= 1.0);
    const auto st0 = controller.last_status();
    if (st0 == design::QPStatus::MaxIterations || st0 == design::QPStatus::Degenerate) {
        CHECK(u0(0) == doctest::Approx(0.0)); // no previous move
    }
    x = plant_step(x, u0(0));
    const auto u1 = controller.step(r, x);
    CHECK(u1(0) >= -1.0);
    CHECK(u1(0) <= 1.0);
    const auto st1 = controller.last_status();
    if (st1 != design::QPStatus::Success) {
        CHECK(u1(0) == doctest::Approx(u0(0)));
    }
}

TEST_CASE("hard output rows with ill-scaled u and y boxes stay inside both") {
    design::InputLimits<1> limits{.lower = {-1.0e6}, .upper = {1.0e6}};
    const auto             art = design::mpc<20, 8>(double_integrator, move_weight(0.01), limits, design::OutputLimits{.outputs = {0}, .upper = {1.02}, .upper_ecr = {0.0}}, design::TrackedOutputs{});
    REQUIRE(art.success);
    MPC          controller{art};
    ColVec<2>    x{};
    const ColVec r{1.0};
    for (size_t k = 0; k < 80; ++k) {
        const auto u = controller.step(r, x);
        CHECK(u(0) >= -1.0e6);
        CHECK(u(0) <= 1.0e6);
        x = plant_step(x, u(0));
        CHECK(x(0) <= 1.02 + 1e-6);
    }
}

TEST_CASE("measured-disturbance feedforward cancels a known load") {
    // The disturbance enters like the input (Bd = B). With the MD channel the
    // controller sees it in the prediction and cancels it exactly at steady
    // state; without it, the same load leaves a residual error.
    const Matrix<2, 1> Bd = double_integrator.B;

    auto weights = move_weight(0.05);

    design::InputLimits<1> limits{.lower = {-2.0}, .upper = {2.0}};

    const ColVec r{0.5};
    const double d_load = 0.4;

    // With MD feedforward:
    const auto art_md = design::mpc<15, 5>(double_integrator, weights, limits, design::MeasuredDisturbance{Bd}, design::TrackedOutputs{});
    REQUIRE(art_md.success);
    MPC       md_controller{art_md};
    ColVec<2> x_md{};
    double    u_md = 0.0;
    for (size_t k = 0; k < 400; ++k) {
        const auto u = md_controller.step(r, ColVec<1>{d_load}, x_md);
        u_md = u(0);
        x_md = plant_step(x_md, u_md + d_load);
    }

    // Baseline without MD, same load applied to the plant:
    const auto art_base = design::mpc<15, 5>(double_integrator, weights, limits, design::TrackedOutputs{});
    REQUIRE(art_base.success);
    MPC       base_controller{art_base};
    ColVec<2> x_base{};
    for (size_t k = 0; k < 400; ++k) {
        const auto u = base_controller.step(r, x_base);
        x_base = plant_step(x_base, u(0) + d_load);
    }

    // MD version: exact tracking, input settles at -d (net force zero).
    CHECK(x_md(0) == doctest::Approx(0.5).epsilon(1e-5));
    CHECK(u_md == doctest::Approx(-d_load).epsilon(1e-4));
    // And it beats the disturbance-blind baseline.
    CHECK(damp::abs(x_md(0) - 0.5) < damp::abs(x_base(0) - 0.5));
}

TEST_CASE("input target pulls the input when tracking is indifferent") {
    design::MPCWeights<1, 1> weights{};
    weights.Qy = Matrix<1, 1>{};          // no output cost:
    weights.Qy_terminal = Matrix<1, 1>{}; // pure input shaping
    weights.Ru = 1.0;
    weights.Rdu = 0.1;

    const auto art = design::mpc<5, 5>(double_integrator, weights, design::TrackedOutputs{});
    REQUIRE(art.success);

    MPC             controller{art};
    ColVec<2>       x{};
    const ColVec<1> r{};
    CHECK(controller.input_target()(0) == 0.0); // the operating point is the runtime's, zero until set
    controller.set_input_target(0.3);

    ColVec<1> u{};
    for (size_t k = 0; k < 100; ++k) {
        u = controller.step(r, x);
        x = plant_step(x, u(0));
    }
    CHECK(u(0) == doctest::Approx(0.3).epsilon(1e-6));

    // A new target takes effect without re-synthesis.
    controller.set_input_target(ColVec{-0.2});
    for (size_t k = 0; k < 100; ++k) {
        u = controller.step(r, x);
        x = plant_step(x, u(0));
    }
    CHECK(u(0) == doctest::Approx(-0.2).epsilon(1e-6));
}

TEST_CASE("automatic scaling is consistent with the unscaled problem") {
    // Finite ranges give su = sy = 20; compensating the weights by the scale
    // factors must reproduce the sentinel-range (unit-scale) trajectory.
    design::InputLimits<1>     wide{.lower = {-10.0}, .upper = {10.0}};
    const design::OutputLimits wide_y{.outputs = {0}, .lower = {-10.0}, .upper = {10.0}};

    auto weights_scaled = move_weight(0.05 * 400.0); // 0.05 · su²
    weights_scaled.Qy(0, 0) = 400.0;                 // 1 · sy²
    weights_scaled.Qy_terminal = weights_scaled.Qy;

    const auto weights_unit = move_weight(0.05);

    const auto art_scaled = design::mpc<12, 4>(double_integrator, weights_scaled, wide, wide_y, design::TrackedOutputs{});
    const auto art_unit = design::mpc<12, 4>(double_integrator, weights_unit, design::TrackedOutputs{});
    REQUIRE(art_scaled.success);
    REQUIRE(art_unit.success);
    CHECK(art_scaled.mpc.scale_u(0) == doctest::Approx(20.0));
    CHECK(art_scaled.mpc.scale_y(0) == doctest::Approx(20.0));

    MPC          scaled{art_scaled};
    MPC          unit{art_unit};
    ColVec<2>    xs{};
    ColVec<2>    xu{};
    const ColVec r{1.0};

    for (size_t k = 0; k < 100; ++k) {
        const auto us = scaled.step(r, xs);
        const auto uu = unit.step(r, xu);
        CHECK(us(0) == doctest::Approx(uu(0)).epsilon(1e-8));
        xs = plant_step(xs, us(0));
        xu = plant_step(xu, uu(0));
    }
}

TEST_CASE("scaling handles badly-conditioned mixed units") {
    // Milli-scale actuator driving a kilo-scale output: DC gain 1e6.
    constexpr StateSpace plant{
        .A = {{0.9}},
        .B = {{1000.0}},
        .C = {{100.0}},
        .Ts = 0.01,
    };

    design::InputLimits<1>     limits{.lower = {-1.0e-3}, .upper = {1.0e-3}};
    const design::OutputLimits range{.outputs = {0}, .lower = {0.0}, .upper = {1000.0}};

    const auto art = design::mpc<10, 4>(plant, design::MPCWeights<1, 1>{}, limits, range, design::TrackedOutputs{});
    REQUIRE(art.success);

    MPC          controller{art};
    ColVec<1>    x{};
    const ColVec r{500.0};

    for (size_t k = 0; k < 200; ++k) {
        const auto u = controller.step(r, x);
        CHECK(controller.last_status() == design::QPStatus::Success);
        CHECK(u(0) <= 1.0e-3 + 1e-12);
        CHECK(u(0) >= -1.0e-3 - 1e-12);
        x = ColVec<1>{0.9 * x(0) + 1000.0 * u(0)};
    }
    CHECK(100.0 * x(0) == doctest::Approx(500.0).epsilon(0.02));
}

TEST_CASE("mpc rejects invalid plants") {
    // Continuous-time (Ts = 0): must be discretized first.
    StateSpace<2, 1, 1> continuous = double_integrator;
    continuous.Ts = 0.0;
    CHECK(design::mpc<5, 5>(continuous, {}, design::TrackedOutputs{}).status == design::DesignStatus::invalid_sample_time);

    // Feedthrough (D != 0) is unsupported.
    StateSpace<2, 1, 1> feedthrough = double_integrator;
    feedthrough.D(0, 0) = 1.0;
    CHECK(!design::mpc<5, 5>(feedthrough, {}, design::TrackedOutputs{}).success);

    // Non-positive slack penalty cannot keep the Hessian PD.
    design::MPCWeights<1, 1> weights{};
    weights.ecr_weight = 0.0;
    CHECK(!design::mpc<5, 5>(double_integrator, weights, design::TrackedOutputs{}).success);
}

TEST_CASE("mpc is constexpr") {
    constexpr auto art = design::mpc<3, 3>(double_integrator, {}, design::TrackedOutputs{});
    static_assert(art.success);
    static_assert(art.mpc.max_qp_iterations > 0);
}

namespace {

/// Whether a law reports a state estimate: one with a filter.
template<typename L>
concept estimates_state = requires(const L& law) { law.state_estimate(); };

/// Whether a law reports a disturbance estimate: one whose filter estimates a disturbance.
template<typename L>
concept estimates_disturbance = requires(const L& law) { law.disturbance_estimate(); };

/// Whether a law's step takes these arguments.
template<typename L, typename... A>
concept steps_with = requires(L& law, const A&... a) { law.step(a...); };

/// Whether a law's control takes these arguments.
template<typename L, typename... A>
concept controls_with = requires(L& law, const A&... a) { law.control(a...); };

/// Whether a law's track takes these arguments.
template<typename L, typename... A>
concept tracks_with = requires(L& law, const A&... a) { law.track(a...); };

} // namespace

TEST_CASE("Law names an MPC member's type by its design") {
    static constexpr auto art = design::mpc<3, 3>(double_integrator, {}, design::TrackedOutputs{});
    struct Loop {
        Law<art, double> mpc; // the design's law
    };
    Loop         loop{};
    MPC          plain{art};
    const ColVec r{1.0};
    const ColVec x{0.0, 0.0};
    CHECK(loop.mpc.step(r, x)(0) == plain.step(r, x)(0));
    static_assert(!estimates_state<decltype(plain)> && !estimates_disturbance<decltype(plain)>); // no filter: the state is measured
    static_assert(std::is_same_v<decltype(plain.estimator), NoEstimator> && std::is_empty_v<NoEstimator>);

    //! With an OutputLimits, and offset-free with one: the law's QP carries the bounded outputs' rows.
    static constexpr design::OutputLimits cap{.outputs = {0}, .upper = {0.5}, .upper_ecr = {0.0}};
    static constexpr auto                 capped = design::mpc<3, 3>(double_integrator, {}, cap, design::DareCost{}, design::TrackedOutputs{});
    Law<capped, double>                   capped_law{};
    MPC                                   capped_plain{capped};
    static_assert(std::is_same_v<std::remove_cvref_t<decltype(capped_law.limits().output)>, Bounds<1>> && decltype(capped_plain)::NI1 == (2 * 3) + (2 * 3) + 1);
    CHECK(capped_law.step(r, x)(0) == capped_plain.step(r, x)(0));
    static constexpr StateSpace<2, 1, 1, double, 2, 1> measured{.A = double_integrator.A, .B = double_integrator.B, .C = double_integrator.C, .Ts = Ts};
    static constexpr auto                              of = design::mpc<3, 3>(measured, {}, design::KalmanNoise<2, 1, double, 1>{}, design::InputDisturbance{}, cap);
    Law<of, double>                                    of_law{};
    OFMPC                                              of_plain{of};
    CHECK(of_law.step(r, ColVec{0.0})(0) == doctest::Approx(of_plain.step(r, ColVec{0.0})(0)).epsilon(1e-12));
}

TEST_CASE("MPC through a filter without a disturbance: update(y) on the outputs, and Law runs it on constants") {
    //! design::mpc with a KalmanNoise and no disturbance spec: the filter reconstructs the state
    //! from the position alone and the QP is posed about the target of r. The plant starts away
    //! from the filter's estimate, so the filter has to find it, and then the loop settles on r.
    static constexpr auto kRes = design::mpc<15, 5>(double_integrator, design::MPCWeights<1, 1>{.Rdu = 0.05}, design::KalmanNoise{.process = diag(1e-4, 1e-4), .measurement = 1e-4}, design::TrackedOutputs{});
    static_assert(kRes.success);
    MPC law{kRes};
    using Filtered = decltype(law);
    static_assert(Filtered::NYM == 1 && Filtered::NDE == 0 && Filtered::NM == 0); // the reading, no disturbance state, nothing measured
    static_assert(estimates_state<Filtered> && !estimates_disturbance<Filtered>);
    static_assert(std::is_same_v<decltype(law.estimator), FixedGainObserver<2, 1, 1, double>>);
    static_assert(std::is_same_v<decltype(law.limits().disturbance), Bounds<0>>);
    static_assert(std::is_same_v<Law<kRes, double>::Runtime, Filtered>);
    static_assert(sizeof(Law<kRes, float>) < sizeof(Law<kRes, float>::Runtime)); // the filter's matrices are constants of the type

    Law<kRes, double> fixed;
    ColVec            x{0.4, 0.0};
    double            worst = 0.0;
    for (int k = 0; k < 300; ++k) {
        const ColVec<1> y{x(0)};
        const ColVec<1> u = law.step(ColVec{1.0}, y);
        worst = std::fmax(worst, std::fabs(u(0) - fixed.step(ColVec{1.0}, y)(0)));
        x = plant_step(x, u(0));
    }
    CHECK(x(0) == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(x(1) == doctest::Approx(0.0).scale(1.0).epsilon(1e-6));
    CHECK((law.state_estimate() - x).norm() < 1e-6);
    CHECK(worst < 1e-9);
}

TEST_CASE("MPC through a filter with neither a TrackedOutputs nor a disturbance spec is a regulator: step(y)") {
    //! The noise alone: the regulator of the measured state, through the filter. It takes no r,
    //! brings an unknown state to rest at 0, and its Law on constants gives the same commands.
    static constexpr auto kRes = design::mpc<15, 5>(double_integrator, design::MPCWeights<1, 1>{.Rdu = 0.05}, design::KalmanNoise{.process = diag(1e-4, 1e-4), .measurement = 1e-4});
    static_assert(kRes.success);
    MPC law{kRes};
    using Regulator = decltype(law);
    using V1 = ColVec<1>;
    static_assert(Regulator::NYM == 1 && Regulator::NDE == 0 && estimates_state<Regulator>);
    static_assert(controls_with<Regulator> && steps_with<Regulator, V1> && tracks_with<Regulator, V1, V1>); // control(), step(y), track(u, y)
    static_assert(!controls_with<Regulator, V1> && !steps_with<Regulator, V1, V1>);                         // no r
    static_assert(std::is_same_v<Law<kRes, double>::Runtime, Regulator>);

    Law<kRes, double> fixed;
    ColVec            x{0.4, 0.0};
    double            worst = 0.0;
    for (int k = 0; k < 300; ++k) {
        const V1 y{x(0)};
        const V1 u = law.step(y);
        worst = std::fmax(worst, std::fabs(u(0) - fixed.step(y)(0)));
        x = plant_step(x, u(0));
    }
    CHECK(x.norm() < 1e-6);
    CHECK((law.state_estimate() - x).norm() < 1e-6);
    CHECK(worst < 1e-9);
}

TEST_CASE("Law of a failed MPC design: a constexpr one does not compile, one made at start-up holds its input") {
    //! Law refuses a failed constexpr design with a static_assert whose error names its status, on
    //! the measured state as through a filter (checked by hand: a requires expression cannot see a
    //! static_assert). A design initialized at start-up is not known at compile time: its Law
    //! forms, valid() is false, and it holds its input, as an MPC built from it does.
    static const auto kStartUp = design::mpc<10, 4>(double_integrator, design::MPCWeights<1, 1>{.Qy = Matrix{{-1.0}}}, design::TrackedOutputs{});
    REQUIRE(kStartUp.status == design::DesignStatus::bad_weights);
    Law<kStartUp, double> law;
    CHECK_FALSE(law.valid());
    CHECK(law.step(ColVec{1.0}, ColVec{0.5, 0.0})(0) == 0.0);
}

TEST_CASE("Law<kRes, T> of an OFMPC on the design's constants takes set_limits as OFMPC does") {
    static constexpr StateSpace<2, 1, 1, double, 2, 1> measured{.A = double_integrator.A, .B = double_integrator.B, .C = double_integrator.C, .Ts = Ts};
    static constexpr design::OutputLimits              kSoftCap{.outputs = {0}, .upper = {0.8}};
    static constexpr auto                              kOnConstants = design::mpc<3, 3>(measured, {}, design::KalmanNoise<2, 1, double, 1>{}, design::InputDisturbance{}, kSoftCap);
    Law<kOnConstants, double>                          law{};
    OFMPC                                              plain{kOnConstants};
    auto                                               lim = law.limits();
    lim.input = {.lower = -0.4, .upper = 0.4};
    lim.output = {.upper = 0.5};
    lim.disturbance = {.lower = -0.5, .upper = 0.5};
    REQUIRE(law.set_limits(lim));
    REQUIRE(plain.set_limits(law.limits()));
    CHECK(plain.limits().disturbance.upper(0) == 0.5);

    ColVec<2> x{};
    for (size_t k = 0; k < 50; ++k) {
        const auto u = law.step(ColVec{0.6}, ColVec<1>{x(0)});
        REQUIRE(law.last_status() == design::QPStatus::Success);
        CHECK(plain.step(ColVec{0.6}, ColVec<1>{x(0)})(0) == doctest::Approx(u(0)).epsilon(1e-12));
        CHECK(std::fabs(u(0)) <= 0.4 + 1e-12);
        if (k == 10) {
            CHECK(u(0) == doctest::Approx(0.4).epsilon(1e-12)); // on the new limit, the design's being unbounded
        }
        if (k == 20) {
            CHECK(u(0) == doctest::Approx(-0.4).epsilon(1e-12));
        }
        x = plant_step(x, u(0));
    }
    CHECK(law.limits().output.upper(0) == 0.5);
}

TEST_CASE("an MPC with an InputRateLimits is the size every MPC was with Δu rows; without one it is smaller") {
    //! Recorded on x86-64 before InputRateLimits, when every QP carried the Δu rows (and with
    //! OFMPC's disturbance limit a public member: limits() and set_limits() read and write the same
    //! member, and a measured-state MPC holds none). With the option the law is that size; without
    //! it the 2·NC·NU rows, their bounds and their state map are gone.
    static constexpr StateSpace<2, 1, 1, double, 2, 1> measured{.A = double_integrator.A, .B = double_integrator.B, .C = double_integrator.C, .Ts = Ts};
    static constexpr design::InputLimits<1>            kHalf{.lower = {-0.5}, .upper = {0.5}};
    static constexpr design::InputRateLimits<1>        kMoves{.lower = {-0.2}, .upper = {0.2}};
    static constexpr design::OutputLimits              kCap{.outputs = {0}, .upper = {0.8}};
    static constexpr design::StateLimits               kSlow{.states = {1}, .lower = {-0.3}, .upper = {0.3}};
    static constexpr auto                              kMpc = design::mpc<10, 4>(double_integrator, {}, kHalf, kMoves, design::TrackedOutputs{});
    static constexpr auto                              kCapped = design::mpc<10, 4>(double_integrator, {}, kHalf, kMoves, kCap, kSlow, design::TrackedOutputs{});
    static constexpr auto                              kOf = design::mpc<10, 4>(measured, {}, kHalf, kMoves, design::KalmanNoise<2, 1, double, 1>{}, design::InputDisturbance{});
    static constexpr auto                              kOfCapped = design::mpc<10, 4>(measured, {}, kHalf, kMoves, kCap, design::KalmanNoise<2, 1, double, 1>{}, design::InputDisturbance{});
    static constexpr auto                              kFree = design::mpc<10, 4>(double_integrator, {}, kHalf, design::TrackedOutputs{});
    static constexpr auto                              kFreeCapped = design::mpc<10, 4>(double_integrator, {}, kHalf, kCap, kSlow, design::TrackedOutputs{});
    static constexpr auto                              kFreeOf = design::mpc<10, 4>(measured, {}, kHalf, design::KalmanNoise<2, 1, double, 1>{}, design::InputDisturbance{});
    static constexpr auto                              kFreeOfCapped = design::mpc<10, 4>(measured, {}, kHalf, kCap, design::KalmanNoise<2, 1, double, 1>{}, design::InputDisturbance{});
    static_assert(std::is_same_v<decltype(MPC{kMpc}.limits().disturbance), Bounds<0>>);
    static_assert(std::is_same_v<decltype(OFMPC{kOf}.limits().disturbance), Bounds<1>>);
    if constexpr (sizeof(void*) == 8) {
        CHECK(sizeof(MPC{kMpc.as<float>()}) == 1360);
        CHECK(sizeof(MPC{kCapped.as<float>()}) == 2840);
        CHECK(sizeof(OFMPC{kOf.as<float>()}) == 1592);
        CHECK(sizeof(OFMPC{kOfCapped.as<float>()}) == 2416);
        CHECK(sizeof(MPC{kMpc}) == 2600);
        CHECK(sizeof(MPC{kCapped}) == 5552);
        CHECK(sizeof(OFMPC{kOf}) == 3048);
        CHECK(sizeof(OFMPC{kOfCapped}) == 4688);
        CHECK(sizeof(Law<kMpc, float>) == 1360);
        CHECK(sizeof(Law<kCapped, float>) == 2840);
        CHECK(sizeof(Law<kOf, float>) == 1504);
        CHECK(sizeof(Law<kOfCapped, float>) == 2328);

        CHECK(sizeof(MPC{kFree.as<float>()}) == 1064);
        CHECK(sizeof(MPC{kFreeCapped.as<float>()}) == 2544);
        CHECK(sizeof(OFMPC{kFreeOf.as<float>()}) == 1264);
        CHECK(sizeof(OFMPC{kFreeOfCapped.as<float>()}) == 2088);
        CHECK(sizeof(MPC{kFree}) == 2008);
        CHECK(sizeof(MPC{kFreeCapped}) == 4968);
        CHECK(sizeof(OFMPC{kFreeOf}) == 2392);
        CHECK(sizeof(OFMPC{kFreeOfCapped}) == 4040);
        CHECK(sizeof(Law<kFree, float>) == 1064);
        CHECK(sizeof(Law<kFreeCapped, float>) == 2544);
        CHECK(sizeof(Law<kFreeOf, float>) == 1176);
        CHECK(sizeof(Law<kFreeOfCapped, float>) == 2000);
    }
    static_assert(sizeof(MPC{kFree}) < sizeof(MPC{kMpc}) && sizeof(OFMPC{kFreeOfCapped.as<float>()}) < sizeof(OFMPC{kOfCapped.as<float>()}));
}

TEST_CASE("InputRateLimits: the QP has its 2·NC·NU rows on Δu only with one, and a design without one plans as one left unbounded") {
    //! An unbounded InputRateLimits is the QP every design had before the option, rows that never
    //! bind: without one the QP is smaller and the commands are the same, on every kind of design.
    constexpr design::InputRateLimits<1> kUnbounded{};
    constexpr design::InputLimits<1>     kUnit{.lower = {-1.0}, .upper = {1.0}};
    const auto                           w = move_weight(0.05);
    const auto                           state = [](const ColVec<2>& x) { return x; };
    const auto                           position = [](const ColVec<2>& x) { return ColVec<1>{x(0)}; };
    const auto                           same = [](const auto& without, const auto& with, auto make, auto measure) {
        REQUIRE(without.success);
        REQUIRE(with.success);
        using A = std::remove_cvref_t<decltype(without.mpc)>;
        using B = std::remove_cvref_t<decltype(with.mpc)>;
        static_assert(A::NZR == 0 && B::NZR == A::NZ && B::NI == A::NI + (2 * A::NZ));
        static_assert(decltype(make(without))::NI1 == A::NI1 && decltype(make(with))::NI1 == B::NI1);
        CHECK(without.mpc.H == with.mpc.H);
        CHECK(without.mpc.F_x == with.mpc.F_x);
        auto      la = make(without);
        auto      lb = make(with);
        ColVec<2> x{};
        double    worst = 0.0;
        for (int k = 0; k < 80; ++k) {
            const auto ua = la.step(ColVec{1.0}, measure(x));
            const auto ub = lb.step(ColVec{1.0}, measure(x));
            CHECK(la.last_status() == lb.last_status());
            worst = damp::max(worst, damp::abs(ua(0) - ub(0)));
            x = plant_step(x, ua(0) + 0.05); // a load at the input
        }
        CHECK(worst <= 1e-12);
    };
    const auto mpc = [](const auto& res) { return MPC{res}; };
    const auto ofmpc = [](const auto& res) { return OFMPC{res}; };
    same(design::mpc<10, 4>(double_integrator, w, kUnit, design::TrackedOutputs{}), design::mpc<10, 4>(double_integrator, w, kUnit, kUnbounded, design::TrackedOutputs{}), mpc, state);
    same(design::mpc<10, 4>(double_integrator, w, design::TrackedOutputs{}), design::mpc<10, 4>(double_integrator, w, kUnbounded, design::TrackedOutputs{}), mpc, state);
    same(design::mpc<30, 10>(double_integrator, move_weight(0.01), soft_box, soft_y, design::TrackedOutputs{}), design::mpc<30, 10>(double_integrator, move_weight(0.01), soft_box, kUnbounded, soft_y, design::TrackedOutputs{}), mpc, state);
    constexpr design::StateLimits kSlow{.states = {1}, .lower = {-0.3}, .upper = {0.3}};
    same(design::mpc<20, 8>(double_integrator, w, kUnit, kSlow, design::TrackedOutputs{}), design::mpc<20, 8>(double_integrator, w, kUnit, kSlow, kUnbounded, design::TrackedOutputs{}), mpc, state);
    same(design::mpc<8, 8>(double_integrator, w, kUnit, design::DareCost{}, design::TrackedOutputs{}), design::mpc<8, 8>(double_integrator, w, kUnit, design::DareCost{}, kUnbounded, design::TrackedOutputs{}), mpc, state);
    constexpr design::OutputLimits kBox{.outputs = {0}, .lower = {-1.5}, .upper = {1.5}};
    same(design::mpc<8, 8>(double_integrator, w, kUnit, kBox, design::InvariantSet<128>{}, design::TrackedOutputs{}), design::mpc<8, 8>(double_integrator, w, kUnit, kBox, design::InvariantSet<128>{}, kUnbounded, design::TrackedOutputs{}), mpc, state);
    constexpr StateSpace<2, 1, 1, double, 2, 1>    measured{.A = double_integrator.A, .B = double_integrator.B, .C = double_integrator.C, .Ts = Ts};
    constexpr design::KalmanNoise<2, 1, double, 1> noise{.process = Matrix<2, 2>::identity() * 1e-3, .measurement = Matrix{{1e-4}}, .disturbance = Matrix{{1e-2}}};
    same(design::mpc<10, 4>(measured, w, kUnit, noise, design::InputDisturbance{}), design::mpc<10, 4>(measured, w, kUnit, noise, design::InputDisturbance{}, kUnbounded), ofmpc, position);
}

TEST_CASE("an MPC's move limit is in its limits() only with an InputRateLimits, set_limits moves it, and an InvariantSet keeps it") {
    constexpr design::InputLimits<1>     limits{.lower = {-0.5}, .upper = {0.5}};
    constexpr design::InputRateLimits<1> moves{.lower = {-0.2}, .upper = {0.2}};
    MPC                                  rated{design::mpc<20, 8>(double_integrator, move_weight(0.05), limits, moves, design::TrackedOutputs{})};
    MPC                                  plain{design::mpc<20, 8>(double_integrator, move_weight(0.05), limits, design::TrackedOutputs{})};
    static_assert(std::is_same_v<decltype(rated.limits().input_rate), Bounds<1>> && std::is_same_v<decltype(plain.limits().input_rate), Bounds<0>>);
    CHECK(plain.set_limits(plain.limits())); // the empty field: nothing to move

    auto lim = rated.limits();
    CHECK(lim.input_rate.upper(0) == 0.2);
    lim.input_rate = {.lower = -0.05, .upper = 0.05};
    REQUIRE(rated.set_limits(lim));
    ColVec<2> x{};
    double    prev = 0.0;
    for (int k = 0; k < 40; ++k) {
        const double u = rated.step(ColVec{1.0}, x)(0);
        REQUIRE(rated.last_status() == design::QPStatus::Success);
        CHECK(damp::abs(u - prev) <= 0.05 + 1e-12); // the new move limit, not the design's 0.2
        prev = u;
        x = plant_step(x, u);
    }
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const Bounds<1>& bad : {Bounds<1>{.lower = 0.1, .upper = -0.1}, Bounds<1>{.lower = nan}, Bounds<1>{.upper = nan}}) {
        auto refused = rated.limits();
        refused.input_rate = bad;
        CHECK_FALSE(rated.set_limits(refused));
        CHECK(rated.limits().input_rate.upper(0) == 0.05); // nothing changed
    }

    //! The terminal set was built for the move limit too: the limit in force is taken, another refused.
    constexpr design::OutputLimits position{.outputs = {0}, .lower = {-1.0}, .upper = {1.0}};
    MPC                            terminal{design::mpc<8, 8>(double_integrator, move_weight(0.05), design::InputLimits<1>{.lower = {-1.0}, .upper = {1.0}}, moves, position, design::InvariantSet<128>{}, design::TrackedOutputs{})};
    REQUIRE(terminal.valid());
    auto in_force = terminal.limits();
    CHECK(terminal.set_limits(in_force));
    in_force.input_rate.upper(0) = 0.1;
    CHECK_FALSE(terminal.set_limits(in_force));
    CHECK(terminal.limits().input_rate.upper(0) == 0.2);
}

TEST_CASE("Law<kRes, float> of a constexpr MPC with an InputRateLimits is the runtime's law") {
    static constexpr auto kRated = design::mpc<10, 4>(double_integrator, design::MPCWeights<1, 1>{.Rdu = {0.05}}, design::InputLimits<1>{.lower = {-0.5}, .upper = {0.5}}, design::InputRateLimits<1>{.lower = {-0.1}, .upper = {0.1}}, design::TrackedOutputs{});
    static_assert(kRated.success);
    static constinit Law<kRated, float> at_start{};
    Law<kRated, float>                  law{};
    MPC                                 held{kRated.as<float>()};
    static_assert(std::is_same_v<Law<kRated, float>::Runtime, decltype(held)> && decltype(held)::NI1 == (4 * 4) + 1);
    CHECK(at_start.limits().input_rate.upper(0) == 0.1f);
    ColVec<2, float> x{};
    float            prev = 0.0f;
    for (int k = 0; k < 40; ++k) {
        const auto u = law.step(ColVec{1.0f}, x);
        CHECK(u == held.step(ColVec{1.0f}, x));
        CHECK(damp::abs(u(0) - prev) <= 0.1f + 1e-6f);
        prev = u(0);
        x = ColVec<2, float>{plant_step(x.as<double>(), static_cast<double>(u(0))).as<float>()};
    }
}

TEST_CASE("artifact float conversion maps the unbounded sentinels safely") {
    constexpr auto art = design::mpc<3, 3>(double_integrator, {}, design::TrackedOutputs{});
    constexpr auto art_f = art.as<float>();
    static_assert(art_f.success);

    // Default constraints are unbounded: the double sentinel must arrive as
    // the float sentinel, not as an overflowed cast.
    static_assert(art_f.mpc.constraints.limits.input.upper(0) == std::numeric_limits<float>::max());
    static_assert(art_f.mpc.constraints.limits.input.lower(0) == -std::numeric_limits<float>::max());
    static_assert(art_f.mpc.b_bound(0) == std::numeric_limits<float>::max());
}

TEST_CASE("MPC runtime runs in float") {
    auto weights = move_weight(0.05);

    design::InputLimits<1> limits{.lower = {-0.5}, .upper = {0.5}};

    const auto art = design::mpc<10, 4>(double_integrator, weights, limits, design::TrackedOutputs{});
    REQUIRE(art.success);

    MPC              controller{art.as<float>()};
    ColVec<2, float> x{};
    const ColVec     r{1.0F};

    for (size_t k = 0; k < 150; ++k) {
        const auto u = controller.step(r, x);
        CHECK(u(0) <= 0.5F + 1e-6F);
        CHECK(u(0) >= -0.5F - 1e-6F);
        x = ColVec<2, float>{
            x(0) + (0.1F * x(1)) + (0.005F * u(0)),
            x(1) + (0.1F * u(0)),
        };
    }
    CHECK(x(0) == doctest::Approx(1.0).epsilon(0.05));
}

namespace {

// Solver policy that counts invocations and forwards to the default solver —
// exercises the pluggable-solver seam with a stateful policy.
struct CountingSolver {
    size_t calls = 0;

    template<size_t NV, size_t NI, typename T>
    [[nodiscard]] constexpr design::QPResult<NV, NI, T> operator()(
        const Matrix<NV, NV, T>& H,
        const ColVec<NV, T>&     f,
        const Matrix<NI, NV, T>& A,
        const ColVec<NI, T>&     b,
        size_t                   max_iterations
    ) {
        ++calls;
        return design::solve_qp(H, f, A, b, max_iterations);
    }
};

} // namespace

TEST_CASE("MPC accepts a custom QP solver policy") {
    auto weights = move_weight(0.05);

    const auto art = design::mpc<8, 4>(double_integrator, weights, design::TrackedOutputs{});
    REQUIRE(art.success);

    MPC reference{art};
    MPC custom{art, CountingSolver{}};

    ColVec       xr{0.5, 0.0};
    ColVec       xc{0.5, 0.0};
    const ColVec r{1.0};

    for (size_t k = 0; k < 20; ++k) {
        const auto ur = reference.step(r, xr);
        const auto uc = custom.step(r, xc);
        CHECK(uc(0) == doctest::Approx(ur(0)).epsilon(1e-12));
        xr = plant_step(xr, ur(0));
        xc = plant_step(xc, uc(0));
    }
    CHECK(custom.qp_solver().calls == 20);
}

TEST_CASE("horizon advisor reads the dominant time constant") {
    // First-order plant, pole at 0.9: tau = -Ts/ln(0.9) = 0.9491 s;
    // 4*tau = 3.796 s at Ts = 0.1 -> NP = 38, NC = 20% = 8.
    constexpr StateSpace plant{
        .A = {{0.9}},
        .B = {{1.0}},
        .C = {{1.0}},
        .Ts = 0.1,
    };

    constexpr auto h = design::suggest_mpc_horizon(plant);
    static_assert(h.success);
    CHECK(h.time_constant == doctest::Approx(0.94912).epsilon(1e-4));
    CHECK(h.prediction_horizon == 38);
    CHECK(h.control_horizon == 8);
}

TEST_CASE("horizon advisor covers the dead time, so its horizon is one design::mpc accepts") {
    constexpr StateSpace lag{.A = {{0.9}}, .B = {{0.1}}, .C = {{1.0}}, .Ts = Ts};
    constexpr auto       model = discretize(lag, Ts, InputDelay<5>{5.0 * Ts}).value();
    constexpr auto       h = design::suggest_mpc_horizon(model);
    static_assert(h.success);
    CHECK(h.prediction_horizon == 43); // 4τ = 3.80 s and the 0.5 s delay, at Ts = 0.1
    CHECK(design::mpc<h.prediction_horizon, h.control_horizon>(model, move_weight(0.05), design::TrackedOutputs{}).success);
    CHECK(design::suggest_mpc_horizon(lag).prediction_horizon == 38); // no delay: unchanged
}

TEST_CASE("horizon advisor falls back to the settling-time overload on integrating plants") {
    // A pure integrator chain has no strictly stable pole to read tau from.
    const auto from_poles = design::suggest_mpc_horizon(double_integrator);
    CHECK(!from_poles.success);

    // The explicit target works anywhere: 2 s at Ts = 0.1 -> NP = 20, NC = 4.
    const auto h = design::suggest_mpc_horizon(2.0, 0.1);
    REQUIRE(h.success);
    CHECK(h.prediction_horizon == 20);
    CHECK(h.control_horizon == 4);
}

TEST_CASE("to_ss of the LQR-equivalent design is the LQR gain on x and no memory of u_prev") {
    const Matrix Q{{1.0, 0.0}, {0.0, 0.1}};
    const Matrix R{{0.5}};
    const auto   lqr = design::discrete_lqr(double_integrator_full.A, double_integrator_full.B, Q, R);
    REQUIRE(lqr.success);

    design::MPCWeights<1, 2> weights{};
    weights.Qy = Q;
    weights.Qy_terminal = lqr.S;
    weights.Rdu = Matrix<1, 1>{};
    weights.Ru = R;

    const auto art = design::mpc<8, 8>(double_integrator_full, weights); // a free move every step: equals LQR
    REQUIRE(art.success);

    // State u_prev, inputs [r; x]: Δu = −K_x x − K_u u_prev with K_x the LQR gain and K_u = 1,
    // so u = u_prev + Δu = −K_lqr x and the block keeps nothing of u_prev (A_c = C_c = 0).
    const auto ss = art.mpc.to_ss();
    CHECK(ss.Ts == doctest::Approx(Ts));
    CHECK(ss.A(0, 0) == doctest::Approx(0.0).scale(1.0).epsilon(1e-8));
    CHECK(ss.C(0, 0) == doctest::Approx(0.0).scale(1.0).epsilon(1e-8));
    CHECK(ss.D(0, 2) == doctest::Approx(-lqr.K(0, 0)).epsilon(1e-8));
    CHECK(ss.D(0, 3) == doctest::Approx(-lqr.K(0, 1)).epsilon(1e-8));
}

TEST_CASE("to_ss realizes the unconstrained runtime, and closed around the plant is stable with unit DC gain") {
    const auto art = design::mpc<15, 5>(double_integrator, move_weight(0.05), design::TrackedOutputs{});
    REQUIRE(art.success);
    const auto ss = art.mpc.to_ss(); // state u_prev, in [r; x], out u

    MPC       runtime{art};
    ColVec<1> xc{};
    ColVec    xp{0.3, -0.1};
    for (int k = 0; k < 12; ++k) {
        const ColVec<1> r{(k < 6) ? 1.0 : -0.5};
        const ColVec<3> in{r[0], xp[0], xp[1]};
        const ColVec<1> u = ss.C * xc + ss.D * in;
        CHECK(u[0] == doctest::Approx(runtime.step(r, xp)[0]).epsilon(1e-9));
        xc = ss.A * xc + ss.B * in;
        xp = plant_step(xp, u[0]);
    }

    // The plant reports [0; x], so the law's inputs [r; x] close by unity positive feedback.
    const StateSpace plant_out{
        .A = double_integrator.A,
        .B = double_integrator.B,
        .C = Matrix{{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}},
        .Ts = Ts,
    };
    StateSpace<0, 3, 3> positive{.Ts = Ts};
    positive.D = -Matrix<3, 3>::identity();
    const auto loop = series(ss, plant_out);
    REQUIRE(loop.has_value());
    const auto cl = feedback(*loop, positive);
    REQUIRE(cl.has_value());

    // All closed-loop poles inside the unit circle...
    const auto eig = mat::compute_eigenvalues(cl->A);
    for (size_t i = 0; i < 3; ++i) {
        CHECK(damp::abs(eig.values[i]) < 1.0);
    }
    // ...and the r → y DC gain is exactly 1 (integral action).
    const auto dc = *eval_frf(*cl, damp::complex<double>{1.0, 0.0});
    CHECK(dc(1, 0).real() == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(dc(1, 0).imag() == doctest::Approx(0.0));
}

TEST_CASE("to_ss takes a measured disturbance between r and x, in the tick's order") {
    const auto art = design::mpc<15, 5>(double_integrator, move_weight(0.05), design::MeasuredDisturbance{double_integrator.B}, design::TrackedOutputs{});
    REQUIRE(art.success);
    const auto ss = art.mpc.to_ss(); // in [r; d; x]
    MPC        runtime{art};
    ColVec<1>  xc{};
    ColVec<2>  xp{};
    for (int k = 0; k < 12; ++k) {
        const ColVec    r{0.5};
        const ColVec<1> d{(k < 4) ? 0.0 : 0.4};
        const ColVec<4> in{r[0], d[0], xp[0], xp[1]};
        const ColVec<1> u = ss.C * xc + ss.D * in;
        CHECK(u[0] == doctest::Approx(runtime.step(r, d, xp)[0]).epsilon(1e-9));
        xc = ss.A * xc + ss.B * in;
        xp = plant_step(xp, u[0] + d[0]);
    }
}

TEST_CASE("MPC reset and bumpless seeding") {
    const auto art = design::mpc<5, 5>(double_integrator, {}, design::TrackedOutputs{});
    REQUIRE(art.success);

    MPC controller{art};
    controller.back_calculate(ColVec{0.3}); // the input another controller applied
    controller.predict();
    CHECK(controller.applied()(0) == doctest::Approx(0.3));

    controller.reset();
    CHECK(controller.applied()(0) == doctest::Approx(0.0));
    CHECK(controller.last_iterations() == 0);
    CHECK(controller.last_slack() == doctest::Approx(0.0));
}

TEST_CASE("MPC step is predict, update, then control; the no-input step holds the last command") {
    auto                   weights = move_weight(0.05);
    design::InputLimits<1> limits{.lower = {-1.0}, .upper = {1.0}};
    const auto             art = design::mpc<15, 5>(double_integrator, weights, limits, design::TrackedOutputs{});
    REQUIRE(art.success);

    // A downstream clamp tighter than the MPC's: the applied input is the base of the next move.
    MPC          fused{art};
    MPC          split{art};
    ColVec<2>    x{};
    ColVec<1>    u_applied{};
    const ColVec r{0.5};
    for (int k = 0; k < 40; ++k) {
        fused.back_calculate(u_applied);
        const ColVec<1> u_fused = fused.step(r, x);
        split.back_calculate(u_applied);
        split.predict();
        split.update(x);
        const ColVec<1> u_split = split.control(r);
        CHECK(u_fused(0) == doctest::Approx(u_split(0)).epsilon(1e-12));
        CHECK(fused.applied()(0) == u_applied(0));
        u_applied = ColVec<1>{damp::clamp(u_fused(0), -0.3, 0.3)};
        x = plant_step(x, u_applied(0));
    }

    // Without an applied input the step holds the command it returned last.
    MPC       implicit{art};
    MPC       explicit_u{art};
    ColVec<2> xi{};
    ColVec<2> xe{};
    for (int k = 0; k < 40; ++k) {
        const ColVec<1> ui = implicit.step(r, xi);
        explicit_u.back_calculate(explicit_u.last_command());
        const ColVec<1> ue = explicit_u.step(r, xe);
        CHECK(ui(0) == doctest::Approx(ue(0)).epsilon(1e-12));
        xi = plant_step(xi, ui(0));
        xe = plant_step(xe, ue(0));
    }
    CHECK(xi(0) == doctest::Approx(0.5).epsilon(1e-3));
}

namespace {

/// A damped mass, the position controlled, continuous.
constexpr StateSpace kDampedMass{.A = {{0.0, 1.0}, {-2.0, -0.6}}, .B = {{0.0}, {1.5}}, .C = {{1.0, 0.0}}};

/// The position at every @p every-th tick of a 3 s step to 1 under @p art's runtime, on the plant sampled at @p h.
template<typename Result>
damp::array<double, 30> step_response(const Result& art, double h, size_t every) {
    const auto              sampled = discretize(kDampedMass, h).value();
    MPC                     law{art};
    ColVec<2>               x{};
    damp::array<double, 30> y{};
    for (size_t k = 0; k < 30 * every; ++k) {
        if (k % every == 0) {
            y[k / every] = x[0];
        }
        x = sampled.A * x + sampled.B * law.step(ColVec{1.0}, x);
    }
    return y;
}

} // namespace

TEST_CASE("the continuous-plant design's QP is the Van Loan integral of its continuous cost") {
    // J = Σ_k ∫ q·(y − r)² dt over each period, + R_Δu·Δu²/Ts per move, + Ts·R_u·u² per move.
    constexpr double         h = 0.1;
    constexpr size_t         NP = 6;
    constexpr size_t         NC = 3;
    design::MPCWeights<1, 1> w{.Qy = 2.0, .Rdu = 0.3, .Ru = 0.5, .Qy_terminal = 2.0};
    const auto               art = design::mpc<NP, NC>(kDampedMass, h, w, design::TrackedOutputs{});
    REQUIRE(art.success);

    // Hand-built: the reference a constant state, the period's cost one Van Loan exponential.
    const Matrix A_aug{{0.0, 1.0, 0.0}, {-2.0, -0.6, 0.0}, {0.0, 0.0, 0.0}};
    const Matrix B_aug{{0.0}, {1.5}, {0.0}};
    const Matrix E{{1.0, 0.0, -1.0}}; // y − r on [x; r]
    const auto   period = design::discretize_lqr_cost(A_aug, B_aug, E.transpose() * (w.Qy * E), Matrix<1, 1>{}, h).value();
    const auto   sampled = discretize(kDampedMass, h).value();
    const ColVec x0{0.4, -0.2};
    const double u_prev = 0.3;
    const double r = 1.0;
    const auto   J = [&](const ColVec<NC>& moves) {
        double    cost = 0.0;
        ColVec<2> x = x0;
        double    u = u_prev;
        for (size_t k = 0; k < NP; ++k) {
            if (k < NC) {
                u += moves[k];
                cost += ((w.Rdu(0, 0) / h) * moves[k] * moves[k]) + (h * w.Ru(0, 0) * u * u);
            }
            const ColVec<3> p{x[0], x[1], r};
            cost += (p.transpose() * period.Q * p)(0, 0) + (2.0 * (p.transpose() * period.N)(0, 0) * u) + (period.R(0, 0) * u * u);
            x = sampled.A * x + sampled.B * ColVec<1>{u};
        }
        return cost;
    };

    // The QP at that point, unit scales and u_target = 0: ½ZᵀHZ + fᵀZ with f = F_x·x̃₀ − F_r·R′.
    const ColVec<3> xa{x0[0], x0[1], u_prev};
    ColVec<NP>      R{};
    for (size_t k = 0; k < NP; ++k) {
        R[k] = r;
    }
    const ColVec<NC> f = (art.mpc.F_x * xa) - (art.mpc.F_r * R);
    const double     scale = damp::abs(J(ColVec<NC>{}));
    for (size_t i = 0; i < NC; ++i) {
        ColVec<NC> ei{};
        ei[i] = 1.0;
        const double g = 0.5 * (J(ei) - J(ei * -1.0));
        CHECK(f[i] == doctest::Approx(g).epsilon(1e-9).scale(scale));
        for (size_t j = 0; j < NC; ++j) {
            ColVec<NC> ej{};
            ej[j] = 1.0;
            const double hij = J(ei + ej) - J(ei) - J(ej) + J(ColVec<NC>{});
            CHECK(art.mpc.H(i, j) == doctest::Approx(hij).epsilon(1e-9).scale(scale));
        }
    }
}

TEST_CASE("halving the period at fixed continuous weights keeps the closed loop") {
    // The same 2 s horizon and 0.4 s of moves at 0.1 s and 0.05 s. The weights are per second,
    // so the cost is the same functional of the response and the closed loop barely moves; the
    // same numbers taken per step (the discrete-plant form) weigh the output twice as often at
    // the finer period and give a different loop.
    const design::MPCWeights<1, 1> w{.Qy = 50.0, .Rdu = 0.5};
    const auto                     coarse = step_response(design::mpc<20, 4>(kDampedMass, 0.1, w, design::TrackedOutputs{}), 0.1, 1);
    const auto                     fine = step_response(design::mpc<40, 8>(kDampedMass, 0.05, w, design::TrackedOutputs{}), 0.05, 2);
    const auto                     coarse_d = step_response(design::mpc<20, 4>(discretize(kDampedMass, 0.1).value(), w, design::TrackedOutputs{}), 0.1, 1);
    const auto                     fine_d = step_response(design::mpc<40, 8>(discretize(kDampedMass, 0.05).value(), w, design::TrackedOutputs{}), 0.05, 2);
    double                         gap = 0.0;
    double                         gap_per_step = 0.0;
    for (size_t j = 0; j < 30; ++j) {
        gap = damp::max(gap, damp::abs(coarse[j] - fine[j]));
        gap_per_step = damp::max(gap_per_step, damp::abs(coarse_d[j] - fine_d[j]));
    }
    MESSAGE("largest position gap, 0.1 s against 0.05 s: continuous weights ", gap, ", per-step weights ", gap_per_step);
    CHECK(fine[29] == doctest::Approx(1.0).epsilon(1e-2)); // settling onto the step at 2.9 s
    CHECK(gap < 0.03);
    CHECK(gap < 0.5 * gap_per_step);
}

TEST_CASE("the continuous-plant terminal cost is the cost-to-go of the period cost: a free move every step gives one law at any horizon") {
    // With NC = NP and V_f the exact cost-to-go of the stage the QP sums, the first move is the
    // infinite-horizon law whatever NP is (Bellman), if and only if the two costs agree.
    const design::MPCWeights<1, 1> w{.Qy = 5.0, .Rdu = 0.2, .Ru = 0.1};
    const auto                     short_h = design::mpc<3, 3>(kDampedMass, 0.1, w, design::DareCost{}, design::TrackedOutputs{});
    const auto                     long_h = design::mpc<7, 7>(kDampedMass, 0.1, w, design::DareCost{}, design::TrackedOutputs{});
    REQUIRE(short_h.success);
    REQUIRE(long_h.success);
    const auto a = short_h.mpc.to_ss(); // the feedback from x and u_prev is exact with a terminal cost
    const auto b = long_h.mpc.to_ss();
    CHECK(a.A(0, 0) == doctest::Approx(b.A(0, 0)).epsilon(1e-8));
    CHECK(a.D(0, 1) == doctest::Approx(b.D(0, 1)).epsilon(1e-8));
    CHECK(a.D(0, 2) == doctest::Approx(b.D(0, 2)).epsilon(1e-8));
}

TEST_CASE("mpc rejects any nonzero feedthrough, whatever its units") {
    StateSpace sys{
        .A = {{1.0, 0.1}, {0.0, 1.0}},
        .B = {{0.005}, {0.1}},
        .C = {{1.0, 0.0}},
        .D = {{1e-13}}, // e.g. y in V, u in MW: real feedthrough, tiny number
        .Ts = 0.1,
    };
    CHECK_FALSE(design::mpc<5, 2>(sys, {}, design::TrackedOutputs{}).success);
    sys.D = Matrix<1, 1>{};
    CHECK(design::mpc<5, 2>(sys, {}, design::TrackedOutputs{}).success);
}

TEST_CASE("suggest_mpc_horizon does not read a defective integrator as a slow pole") {
    // Double integrator in a non-orthogonal basis with a large off-diagonal:
    // the computed eigenvalues sit O(√eps) off 1.
    const Matrix J{{1.0, 100.0}, {0.0, 1.0}};
    const Matrix Tm{{1.0, 0.3}, {0.7, 1.9}};
    const auto   Ti = Tm.inverse();
    REQUIRE(Ti.has_value());
    const StateSpace<2, 1, 1> sys{.A = Tm * J * *Ti, .B = {{0.0}, {1.0}}, .C = {{1.0, 0.0}}, .Ts = 0.1};
    CHECK_FALSE(design::suggest_mpc_horizon(sys).success);

    // A triple integrator's poles land O(eps^(1/3)) = 1e-6 off 1, outside a √eps band: it
    // suggested NP = 8686359.
    const Matrix J3{{1.0, 0.1, 0.005}, {0.0, 1.0, 0.1}, {0.0, 0.0, 1.0}};
    const Matrix T3{{1.0, 0.3, -0.2}, {0.7, 1.9, 0.4}, {0.1, -0.5, 1.3}};
    const auto   T3i = T3.inverse();
    REQUIRE(T3i.has_value());
    const StateSpace<3, 1, 1> triple{.A = T3 * J3 * *T3i, .B = {{0.0}, {0.0}, {1.0}}, .C = {{1.0, 0.0, 0.0}}, .Ts = 0.1};
    CHECK_FALSE(design::suggest_mpc_horizon(triple).success);
    const StateSpace slow{.A = {{0.99, 0.0}, {0.0, 0.5}}, .B = {{1.0}, {1.0}}, .C = {{1.0, 1.0}}, .Ts = 0.1};
    CHECK(design::suggest_mpc_horizon(slow).time_constant == doctest::Approx(-0.1 / damp::log(0.99)));

    // A genuinely stable pole still yields a suggestion.
    const StateSpace stable{.A = {{0.9}}, .B = {{1.0}}, .C = {{1.0}}, .Ts = 0.1};
    CHECK(design::suggest_mpc_horizon(stable).success);
    CHECK_FALSE(design::suggest_mpc_horizon(1e300, 1e-300).success); // NP not representable
}

TEST_CASE("a measured disturbance with a terminal cost settles on r: the target is the steady state under d") {
    //! The target ignored d, so V_f pulled the plan toward the d = 0 steady state: y = 0.969 for r = 1.
    const StateSpace plant{.A = {{0.9}}, .B = {{0.1}}, .C = {{1.0}}, .Ts = 0.1};
    const auto       art = design::mpc<10, 10>(plant, move_weight(0.1), design::DareCost{}, design::MeasuredDisturbance{Matrix{{0.1}}}, design::TrackedOutputs{});
    REQUIRE(art.success);
    MPC       law{art};
    double    x = 0.0;
    ColVec<1> u{};
    for (size_t k = 0; k < 400; ++k) {
        u = law.step(ColVec{1.0}, ColVec{2.0}, ColVec<1>{x});
        x = (0.9 * x) + (0.1 * u(0)) + (0.1 * 2.0);
    }
    CHECK(x == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(u(0) == doctest::Approx(-1.0).epsilon(1e-6));
}

TEST_CASE("as<float>() refuses a Hessian that does not factor in float") {
    //! It kept success: constinit compiled and every float tick returned NotPositiveDefinite.
    design::InputLimits<1> limits{.lower = {-1.0}, .upper = {1.0}};
    const auto             art = design::mpc<40, 10>(double_integrator, move_weight(1e-3), limits, design::TrackedOutputs{});
    REQUIRE(art.success);
    const auto art_f = art.as<float>();
    CHECK_FALSE(art_f.success);
    CHECK(art_f.status == design::DesignStatus::singular);
    CHECK_FALSE(art_f.mpc.success);
}

TEST_CASE("one hard output side by designated initializer leaves the other side soft") {
    //! y_lim_ecr was a Bounds: {.upper = {0.0}} set the lower softness to −max(), a hard lower
    //! bound and ±inf under as<float>().
    const design::InputLimits<1> limits{.lower = {-1.0}, .upper = {1.0}};
    const design::OutputLimits   band{
        .outputs = {0},
        .lower = {0.2},
        .upper = {1.0},
        .upper_ecr = {0.0}
    };
    const auto art = design::mpc<10, 3>(double_integrator, move_weight(0.1), limits, band, design::TrackedOutputs{});
    REQUIRE(art.success);
    MPC        controller{art};
    const auto u = controller.step(ColVec{0.5}, ColVec<2>{}); // y = 0 starts below the soft lower bound
    CHECK(controller.last_status() == design::QPStatus::Success);
    CHECK(u(0) > 0.0);
    constexpr size_t r_y_hi = 4 * 3;             // after the Δu and u rows of NC = 3 moves
    constexpr size_t r_y_lo = r_y_hi + (10 * 1); // after the NP = 10 upper output rows
    CHECK(art.mpc.A_con(r_y_hi, 3) == 0.0);      // the slack column: the upper side hard
    CHECK(art.mpc.A_con(r_y_lo, 3) == -1.0);     // the lower side soft
    CHECK(art.mpc.as<float>().A_con(r_y_lo, 3) == -1.0F);
}

TEST_CASE("design::mpc refuses NaN or inverted limits, negative softness and indefinite weights") {
    //! The design accepted what set_limits refuses (an input box {1, −1}: std::clamp with lo > hi on
    //! every tick), and output weights design::state_feedback refuses.
    design::InputLimits<1> inverted{.lower = {1.0}, .upper = {-1.0}};
    CHECK(design::mpc<6, 2>(double_integrator, move_weight(1.0), inverted, design::TrackedOutputs{}).status == design::DesignStatus::invalid_parameter);
    design::InputRateLimits<1> nan_du{};
    nan_du.upper(0) = std::numeric_limits<double>::quiet_NaN();
    CHECK(design::mpc<6, 2>(double_integrator, move_weight(1.0), nan_du, design::TrackedOutputs{}).status == design::DesignStatus::invalid_parameter);
    CHECK(design::mpc<6, 2>(double_integrator, move_weight(1.0), design::InputRateLimits{.lower = {0.2}, .upper = {-0.2}}, design::TrackedOutputs{}).status == design::DesignStatus::invalid_parameter);
    CHECK(design::mpc<6, 2>(double_integrator, move_weight(1.0), design::OutputLimits{.outputs = {0}, .upper_ecr = {-1.0}}, design::TrackedOutputs{}).status == design::DesignStatus::bad_weights);

    design::MPCWeights<1, 1> negative{};
    negative.Qy = -0.05;
    negative.Qy_terminal = -0.05;
    CHECK(design::mpc<6, 2>(double_integrator, negative, design::TrackedOutputs{}).status == design::DesignStatus::bad_weights);
    design::MPCWeights<1, 2> asymmetric{};
    asymmetric.Qy = Matrix{{1.0, 5.0}, {0.0, 1.0}};
    CHECK(design::mpc<6, 2>(double_integrator_full, asymmetric, design::TrackedOutputs{{0}}).status == design::DesignStatus::bad_weights);
    CHECK(design::mpc<6, 2>(double_integrator, move_weight(1.0), design::TrackedOutputs{}).success);
}

TEST_CASE("a horizon no output sees an input's first move in is refused") {
    //! It succeeded, and the law returned u = 0 at any error: no output in the horizon could
    //! pay for a move, and the move weight charged for one.
    const auto delayed = discretize(double_integrator, Ts, InputDelay{Ts}).value(); // [x; u(k−1)]
    const auto blind = design::mpc<1, 1>(delayed, move_weight(0.05), design::TrackedOutputs{});
    CHECK_FALSE(blind.success);
    CHECK(blind.status == design::DesignStatus::invalid_parameter);

    const auto sees = design::mpc<2, 1>(delayed, move_weight(0.05), design::TrackedOutputs{});
    REQUIRE(sees.success);
    MPC law{sees};
    CHECK(law.step(ColVec{1.0}, ColVec<3>{})(0) > 0.0);

    SUBCASE("a terminal cost weights the delay state, so it sees past the horizon") {
        const auto res = design::mpc<1, 1>(delayed, move_weight(0.05), design::DareCost{}, design::TrackedOutputs{});
        REQUIRE(res.success);
        MPC with_terminal{res};
        CHECK(with_terminal.step(ColVec{1.0}, ColVec<3>{})(0) > 0.0);
    }

    SUBCASE("the plant's relative degree, with no delay at all") {
        const StateSpace<2, 1, 1> euler{.A = double_integrator.A, .B = {{0.0}, {Ts}}, .C = double_integrator.C, .Ts = Ts}; // C·B = 0
        CHECK(design::mpc<1, 1>(euler, move_weight(0.05), design::TrackedOutputs{}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::mpc<2, 1>(euler, move_weight(0.05), design::TrackedOutputs{}).success);
    }

    SUBCASE("one input of two: the other reaching the output does not save it") {
        const StateSpace               two{.A = {{0.9, 0.1}, {0.0, 0.0}}, .B = {{0.1, 0.0}, {0.0, 1.0}}, .C = {{1.0, 0.0}}, .Ts = Ts}; // u₁ a step behind u₀
        const design::MPCWeights<2, 1> w{};
        CHECK(design::mpc<1, 1>(two, w, design::TrackedOutputs{}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::mpc<2, 1>(two, w, design::TrackedOutputs{}).success);
    }

    SUBCASE("relative to the plant's size: a small gain is not a blind input") {
        const StateSpace tiny{.A = double_integrator.A, .B = double_integrator.B * 1e-9, .C = double_integrator.C, .Ts = Ts};
        CHECK(design::mpc<1, 1>(tiny, move_weight(0.05), design::TrackedOutputs{}).success);
        const StateSpace tinier{.A = double_integrator.A, .B = double_integrator.B * 1e-11, .C = double_integrator.C, .Ts = Ts};
        CHECK(design::mpc<2, 1>(discretize(tinier, Ts, InputDelay{Ts}).value(), move_weight(0.05), design::TrackedOutputs{}).success);
    }

    SUBCASE("relative to the data: a delay state, float and the input's units do not refuse a horizon that sees") {
        //! The test compared C A^k B with max|C|·max|B|, and a delay state's 1 set max|B|.
        for (const float unit : {1.0F, 1e3F}) {
            CAPTURE(unit);
            constexpr float                 ts = 1e-3F;
            const StateSpace                plant{.A = {{1.0F, ts}, {0.0F, 1.0F}}, .B = {{0.5F * ts * ts * unit}, {ts * unit}}, .C = {{1.0F, 0.0F}}, .Ts = ts};
            const auto                      model = discretize(plant, ts, InputDelay{ts}).value();
            design::MPCWeights<1, 1, float> w{};
            w.Rdu = 1e-6F;
            CHECK(design::mpc<2, 1>(model, w, design::TrackedOutputs{}).success);
            CHECK(design::mpc<1, 1>(model, w, design::TrackedOutputs{}).status == design::DesignStatus::invalid_parameter);
        }
    }

    SUBCASE("in any basis of the same model: a blind horizon is refused, one that sees is not") {
        //! In a dense basis C A^k B of a blind horizon is roundoff, not 0, and was read as a path.
        const auto model = discretize(double_integrator, Ts, InputDelay<5>{5.0 * Ts}).value(); // a 5-sample delay, 7 states
        for (size_t seed = 1; seed <= 40; ++seed) {
            CAPTURE(seed);
            Matrix<7, 7> Q = Matrix<7, 7>::identity(); // a product of plane rotations: orthogonal and dense
            for (size_t i = 0; i < 7; ++i) {
                for (size_t j = i + 1; j < 7; ++j) {
                    const double a = (0.37 * static_cast<double>(seed * (i + 1))) + (0.11 * static_cast<double>(j));
                    Matrix<7, 7> G = Matrix<7, 7>::identity();
                    G(i, i) = std::cos(a);
                    G(j, j) = std::cos(a);
                    G(i, j) = -std::sin(a);
                    G(j, i) = std::sin(a);
                    Q = Q * G;
                }
            }
            const StateSpace rotated{.A = Q * model.A * Q.transpose(), .B = Q * model.B, .C = model.C * Q.transpose(), .Ts = Ts};
            CHECK(design::mpc<5, 1>(rotated, move_weight(0.05), design::TrackedOutputs{}).status == design::DesignStatus::invalid_parameter);
            CHECK(design::mpc<6, 1>(rotated, move_weight(0.05), design::TrackedOutputs{}).success);
        }
    }

    SUBCASE("an input an Ru weight steers to its target is not blind, while some input sees") {
        //! u₁ reaches y only a step later; at NP = 1 its Ru weight holds it at its target.
        const StateSpace         two{.A = {{0.9, 0.1}, {0.0, 0.0}}, .B = {{0.1, 0.0}, {0.0, 1.0}}, .C = {{1.0, 0.0}}, .Ts = Ts};
        design::MPCWeights<2, 1> w{};
        w.Ru = Matrix{{0.0, 0.0}, {0.0, 1.0}};
        const auto res = design::mpc<1, 1>(two, w, design::TrackedOutputs{});
        REQUIRE(res.success);
        MPC steered{res};
        steered.set_input_target({0.0, 5.0});
        CHECK(steered.step(ColVec{1.0}, ColVec<2>{})(1) > 0.0);

        auto only_ru = move_weight(0.05); // no input reaches the output: refused, Ru or not
        only_ru.Ru = 1.0;
        CHECK(design::mpc<1, 1>(delayed, only_ru, design::TrackedOutputs{}).status == design::DesignStatus::invalid_parameter);
    }

    SUBCASE("the offset-free design too") {
        const StateSpace<1, 1, 1, double, 1, 1>    lag{.A = {{0.9}}, .B = {{0.1}}, .C = {{1.0}}, .Ts = Ts};
        const auto                                 model = discretize(lag, Ts, InputDelay{Ts}).value();
        const design::KalmanNoise<1, 1, double, 1> noise{.process = 0.01, .measurement = 0.001, .disturbance = 0.01};
        const design::InputLimits<1>               none{};
        CHECK(design::mpc<1, 1>(model, move_weight(0.05), none, noise, design::InputDisturbance{}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::mpc<2, 1>(model, move_weight(0.05), none, noise, design::InputDisturbance{}).success);
    }
}

TEST_CASE_TEMPLATE("MPC with soft output bounds, started outside them, solves every tick and converges", T, float, double) {
    //! The float active set reported Infeasible on most such ticks, the input was held and the
    //! double integrator ran off to y = −356. Then it ended Degenerate on 17 to 35 of the 300
    //! ticks from x₀ = 1.5, which this test did not count, and on 295 from x₀ = 20: with the moves
    //! pinned by their limits a row is met through the slack alone, at a sine of 2·10⁻⁴ to the
    //! active normals (qp.hpp).
    const auto art = design::mpc<30, 10>(double_integrator, move_weight(0.01), soft_box, soft_moves, soft_y, design::TrackedOutputs{});
    REQUIRE(art.success);
    using Recorder = RecordingSolver<11, 101, T>;
    MPC controller{art.template as<T>(), Recorder{}};
    MPC reference{art}; // in double, from the same state and held input each tick

    double x0 = 0.0;
    size_t ticks = 0;
    SUBCASE("from x₀ = 1.5") {
        x0 = 1.5;
        ticks = 300;
    }
    SUBCASE("from x₀ = 20: the slack and its multiplier 40 times larger") {
        x0 = 20.0;
        ticks = 500;
    }
    CAPTURE(x0);

    ColVec<2>    x{x0, 0.0};
    size_t       unsolved = 0;
    size_t       reference_unsolved = 0;
    double       worst_move = 0.0;
    KktResiduals worst{};
    for (size_t k = 0; k < ticks; ++k) {
        const ColVec<2, T> x_read{x.template as<T>()}; // the state as the controller reads it
        const auto         u = controller.step(ColVec<1, T>{T{0}}, x_read);
        unsolved += (controller.last_status() != design::QPStatus::Success) ? 1U : 0U;
        const auto u_ref = reference.step(ColVec{0.0}, ColVec<2>{x_read.template as<double>()});
        reference_unsolved += (reference.last_status() != design::QPStatus::Success) ? 1U : 0U;
        worst_move = damp::max(worst_move, damp::abs(static_cast<double>(u(0)) - u_ref(0)));
        reference.back_calculate(ColVec<1>{static_cast<double>(u(0))});
        worst.take_worst(kkt_residuals(controller.qp_solver()));
        x = plant_step(x, static_cast<double>(u(0)));
    }
    const double e = static_cast<double>(eps<T>());
    CHECK(unsolved == 0);
    CHECK(reference_unsolved == 0);
    //! The move against the double controller's. Rounding H and f to ε moves the first move of
    //! an unconstrained plan by up to ε·Σⱼ|H⁻¹|₀ⱼ·(|fⱼ| + Σₖ|Hⱼₖ Zₖ|): the first row of H⁻¹ sums to
    //! 55 and the gradient's terms reach 600 here, 6.6·10⁴·ε in u (8·10⁻³ in float). Measured in
    //! float: 7.5·10⁻⁴ from 1.5, 1.5·10⁻³ from 20.
    CHECK(worst_move <= 65536.0 * e);
    //! The KKT conditions of every tick's solve. Rows: the solver's own tolerance, 64·(NV+1)·ε of
    //! a row's terms, and 16·ε on the rows with a multiplier (measured 0.4·ε). Stationarity:
    //! 4096·ε of the largest term (measured 430·ε): the multipliers reach 4·10⁶.
    CHECK(worst.primal <= ((64.0 * 12.0) + 16.0) * e);
    CHECK(worst.tight <= 16.0 * e);
    CHECK(worst.stationarity <= 4096.0 * e);
    CHECK(worst.lambda_min >= 0.0);
    CHECK(x(0) == doctest::Approx(0.0).scale(1.0).epsilon(1e-3));
    CHECK(x(1) == doctest::Approx(0.0).scale(1.0).epsilon(1e-3));
}

TEST_CASE("float MPC solves with a slack weight 10 times the default") {
    //! The margin the dependence floor of the active set leaves in single precision: the sine of
    //! a row met through the slack falls as 1/√ρ, 3.9·10⁻⁵ here against a floor of 1.1·10⁻⁵. At
    //! 100 times the default it is on the floor, and past that some ticks end Degenerate.
    auto weights = move_weight(0.01);
    weights.ecr_weight = 1.0e6;
    const auto art = design::mpc<30, 10>(double_integrator, weights, soft_box, soft_moves, soft_y, design::TrackedOutputs{});
    REQUIRE(art.success);
    for (const double x0 : {1.5, 20.0}) {
        CAPTURE(x0);
        MPC       controller{art.as<float>()};
        ColVec<2> x{x0, 0.0};
        size_t    unsolved = 0;
        for (size_t k = 0; k < 300; ++k) {
            const auto u = controller.step(ColVec{0.0F}, x.as<float>());
            unsolved += (controller.last_status() != design::QPStatus::Success) ? 1U : 0U;
            x = plant_step(x, static_cast<double>(u(0)));
        }
        CHECK(unsolved == 0);
    }
}

TEST_CASE_TEMPLATE("MPC rests on a soft output bound the reference lies beyond", T, float, double) {
    //! r = 1 above the soft ceiling y ≤ 0.5: the output settles on the ceiling plus the slack the
    //! cost buys, 30·0.5/(30 + ρ) = 1.5·10⁻⁴, with soft rows active on every tick from then on.
    //! On the Gram factor the float input chattered by 10⁻² there.
    const auto art = design::mpc<30, 10>(double_integrator, move_weight(0.01), soft_box, soft_moves, design::OutputLimits{.outputs = {0}, .upper = {0.5}}, design::TrackedOutputs{});
    REQUIRE(art.success);
    MPC controller{art.template as<T>()};
    MPC reference{art};

    ColVec<2> x{};
    size_t    unsolved = 0;
    double    worst_move = 0.0;
    double    y_low = 1.0;
    double    y_high = 0.0;
    double    u_tail = 0.0;
    for (size_t k = 0; k < 400; ++k) {
        const ColVec<2, T> x_read{x.template as<T>()};
        const auto         u = controller.step(ColVec<1, T>{T{1}}, x_read);
        unsolved += (controller.last_status() != design::QPStatus::Success) ? 1U : 0U;
        const auto u_ref = reference.step(ColVec{1.0}, ColVec<2>{x_read.template as<double>()});
        worst_move = damp::max(worst_move, damp::abs(static_cast<double>(u(0)) - u_ref(0)));
        reference.back_calculate(ColVec<1>{static_cast<double>(u(0))});
        x = plant_step(x, static_cast<double>(u(0)));
        if (k >= 300) {
            y_low = damp::min(y_low, x(0));
            y_high = damp::max(y_high, x(0));
            u_tail = damp::max(u_tail, damp::abs(static_cast<double>(u(0))));
        }
    }
    CHECK(unsolved == 0);
    CHECK(worst_move <= 65536.0 * static_cast<double>(eps<T>())); // as above; measured 8·10⁻⁴ in float
    CHECK(y_low > 0.5);
    CHECK(y_high < 0.5 + 2.0e-4);
    CHECK(static_cast<double>(controller.last_slack()) == doctest::Approx(1.5e-4).epsilon(0.05));
    CHECK(u_tail < 1e-3); // at rest on the bound: 2·10⁻⁴ in float
}

TEST_CASE_TEMPLATE("MPC under a hard output ceiling below the reference: a solved tick meets its rows", T, float, double) {
    //! r = 1 above a hard ceiling y ≤ 0.8: the output brakes onto the ceiling and stays, hard
    //! rows active on every tick. A QP on the edge of infeasible is a hard ceiling's lot (the
    //! plan rides it with no margin, and a rounded state can be past it); what a Success must
    //! mean is that the point meets every row. On the Gram factor 200 of these 300 float ticks
    //! returned Success with a row violated by up to 1.2·10⁵·ε of its terms.
    constexpr design::OutputLimits ceiling{.outputs = {0}, .upper = {0.8}, .upper_ecr = {0.0}};
    const auto                     run = [](const auto& art) {
        REQUIRE(art.success);
        using A = std::remove_cvref_t<decltype(art.mpc)>;
        using Recorder = RecordingSolver<A::NZ1, A::NI1, T>;
        constexpr size_t r_y = (2 * A::NZR) + (2 * A::NZ); // the upper output rows, after Δu's (with the move limit) and u's
        MPC              controller{art.template as<T>(), Recorder{}};

        ColVec<2>    x{};
        size_t       stalled = 0;    // MaxIterations or Degenerate: the solver gave up
        size_t       infeasible = 0; // the rounded state past the ceiling's reach
        size_t       on_ceiling = 0; // ticks with a hard output row active
        double       peak = 0.0;
        KktResiduals worst{};
        for (size_t k = 0; k < 300; ++k) {
            const ColVec<2, T> x_read{x.template as<T>()};
            const auto         u = controller.step(ColVec<1, T>{T{1}}, x_read);
            const auto         status = controller.last_status();
            infeasible += (status == design::QPStatus::Infeasible) ? 1U : 0U;
            stalled += (status == design::QPStatus::MaxIterations || status == design::QPStatus::Degenerate) ? 1U : 0U;
            if (status == design::QPStatus::Success) {
                worst.take_worst(kkt_residuals(controller.qp_solver()));
                bool active = false;
                for (size_t row = r_y; row < r_y + 20; ++row) { // the upper output rows
                    active = active || controller.qp_solver().result.lambda(row) > T{0};
                }
                on_ceiling += active ? 1U : 0U;
            }
            x = plant_step(x, static_cast<double>(u(0)));
            peak = damp::max(peak, x(0));
        }
        const double e = static_cast<double>(eps<T>());
        CHECK(stalled == 0);
        CHECK(infeasible <= 6); // 2 to 3 in float, none in double
        CHECK(on_ceiling > 250);
        CHECK(worst.primal <= ((64.0 * 10.0) + 16.0) * e); // the solver's row tolerance; measured 306·ε
        CHECK(worst.tight <= 16.0 * e);
        CHECK(worst.lambda_min >= 0.0);
        CHECK(peak <= 0.8 + 1e-6); // the plant is the model: over by the rounding of the state read, 10⁻⁷ in float
        CHECK(x(0) == doctest::Approx(0.8).epsilon(1e-5));
    };
    SUBCASE("the input free to move") {
        run(design::mpc<20, 8>(double_integrator, move_weight(0.01), soft_box, ceiling, design::TrackedOutputs{}));
    }
    SUBCASE("the moves limited") {
        run(design::mpc<20, 8>(double_integrator, move_weight(0.01), soft_box, soft_moves, ceiling, design::TrackedOutputs{}));
    }
}

TEST_CASE("back_calculate with a feedforward: the command is the law's share, and the loop settles on r") {
    //! The held input already carried u_ff, the user added it again: y = 1.092 for u_ff = 0.2.
    const StateSpace plant{.A = {{0.9}}, .B = {{0.1}}, .C = {{1.0}}, .Ts = 0.1};
    const auto       art = design::mpc<10, 3>(plant, move_weight(0.1), design::TrackedOutputs{});
    REQUIRE(art.success);
    MPC    law{art};
    double x = 0.0;
    double total = 0.0;
    for (size_t k = 0; k < 400; ++k) {
        const auto command = law.step(ColVec{1.0}, ColVec<1>{x});
        total = command(0) + 0.2;
        law.back_calculate(ColVec<1>{total}, ColVec{0.2});
        x = (0.9 * x) + (0.1 * total);
    }
    CHECK(x == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(total == doctest::Approx(1.0).epsilon(1e-6));
}

TEST_CASE("MPCWeights, InputLimits and the limits options take their sizes from the fields they set") {
    using design::InputLimits;
    using design::MPCWeights;
    using design::OutputLimits;
    //! One input and one controlled output unless a field says otherwise.
    constexpr MPCWeights siso{.Rdu = {1.6}};
    static_assert(damp::is_same_v<decltype(siso), const MPCWeights<1, 1>> && siso.Rdu(0, 0) == 1.6);
    constexpr MPCWeights two_outputs{.Qy = diag(20.0, 1.0), .Rdu = {0.1}};
    static_assert(damp::is_same_v<decltype(two_outputs), const MPCWeights<1, 2>> && two_outputs.Qy_terminal == two_outputs.Qy);
    constexpr MPCWeights two_inputs{.Rdu = diag(0.1, 0.2)};
    static_assert(damp::is_same_v<decltype(two_inputs), const MPCWeights<2, 1>>);

    constexpr InputLimits limits{.lower = {0.0}, .upper = {40.0}};
    static_assert(damp::is_same_v<decltype(limits), const InputLimits<1>> && limits.upper(0) == 40.0);
    constexpr InputLimits two_input_limits{.lower = ColVec{0.0, -1.0}, .upper = {40.0, 1.0}};
    static_assert(damp::is_same_v<decltype(two_input_limits), const InputLimits<2>>);
    //! An end left out is unbounded, as a Bounds' is; braces deduce nothing, a typed end gives the inputs and the scalar.
    constexpr InputLimits ceiling{.upper = {1.0}};
    static_assert(damp::is_same_v<decltype(ceiling), const InputLimits<1>> && ceiling.lower(0) == -damp::numeric_limits<double>::max());
    static_assert(damp::is_same_v<decltype(InputLimits{.lower = {-2.0}, .upper = {2.0}}), InputLimits<1>>);
    //! The move limit deduces the same way.
    constexpr design::InputRateLimits moves{.lower = {-0.2}, .upper = {0.2}};
    static_assert(damp::is_same_v<decltype(moves), const design::InputRateLimits<1>> && moves.upper(0) == 0.2);
    constexpr design::InputRateLimits up{.upper = {0.5}};
    static_assert(up.lower(0) == -damp::numeric_limits<double>::max() && up.upper(0) == 0.5);
    constexpr design::InputRateLimits two_moves{.lower = ColVec{-1.0f, -2.0f}, .upper = {1.0f, 2.0f}};
    static_assert(damp::is_same_v<decltype(two_moves), const design::InputRateLimits<2, float>>);
    constexpr design::InputRateLimits<2> spelled{.lower = {-1.0, -1.0}, .upper = {1.0, 1.0}};
    static_assert(spelled.upper(1) == 1.0);
    //! Inside a function template too, on g++ and clang, the spellings above (llvm/llvm-project#114870 is
    //! a field set to a constant declared outside the template: docs/known_limitations.md).
    static_assert(deduced_in_template<0>());

    //! The bounded outputs' count from the list; a side made hard with braces, which deduce nothing.
    constexpr OutputLimits pair{.outputs = {1, 2}, .lower = {-0.4, -0.26}, .upper = {0.4, 0.26}};
    static_assert(damp::is_same_v<decltype(pair), const OutputLimits<2>> && pair.outputs[1] == 2 && pair.lower_ecr(1) == 1.0);
    constexpr OutputLimits hard_cap{.outputs = {0}, .upper = {1.02}, .upper_ecr = {0.0}};
    static_assert(damp::is_same_v<decltype(hard_cap), const OutputLimits<1>> && hard_cap.upper_ecr(0) == 0.0 && hard_cap.lower_ecr(0) == 1.0);
    //! The terminal set's room is named, as its facet count is known only once the set is computed.
    static_assert(damp::is_same_v<decltype(design::InvariantSet<64>{.target_margin = 0.1}), design::InvariantSet<64, double>> && design::InvariantSet<64>{}.target_margin == 0.05);
    static_assert(design::detail::terminal_rows<design::InvariantSet<64>> == 64 && design::detail::terminal_rows<design::DareCost> == 0);
    static_assert(design::QPIterations{40}.budget == 40);
}

TEST_CASE("an MPC takes its reference and input target braced, one value per channel, or as a scalar") {
    design::MPCWeights<1, 2> weights{};
    weights.Ru = Matrix<1, 1>{0.1};
    const auto art = design::mpc<8, 4>(double_integrator_full, weights, design::TrackedOutputs{{0}});
    REQUIRE(art.success);
    const ColVec x{1.0, -0.5};

    //! One input: the target as its scalar or braced, as the ColVec it holds.
    MPC by_vector{art};
    MPC by_braces{art};
    MPC by_scalar{art};
    by_vector.set_input_target(ColVec{-0.2});
    by_braces.set_input_target({-0.2});
    by_scalar.set_input_target(-0.2);
    for (size_t k = 0; k < 4; ++k) {
        CHECK(by_braces.input_target_preview()[k] == ColVec{-0.2});
        CHECK(by_scalar.input_target_preview()[k] == ColVec{-0.2});
    }

    //! One held output: the reference as its scalar, or the ColVec it is.
    const double u_vector = by_vector.step(ColVec{0.3}, x)(0);
    CHECK(by_scalar.step(0.3, x)(0) == doctest::Approx(u_vector).epsilon(1e-12));

    //! Two inputs, both outputs held: the target and the reference braced, one value per channel.
    const StateSpace         two{.A = {{0.9, 0.1}, {0.0, 0.0}}, .B = {{0.1, 0.0}, {0.0, 1.0}}, .C = {{1.0, 0.0}, {0.0, 1.0}}, .Ts = Ts};
    design::MPCWeights<2, 2> w2{};
    w2.Ru = diag(0.1, 0.1);
    const auto art2 = design::mpc<2, 1>(two, w2, design::TrackedOutputs{});
    REQUIRE(art2.success);
    MPC two_inputs{art2};
    MPC two_braced{art2};
    MPC two_mixed{art2};
    MPC two_converted{art2};
    two_inputs.set_input_target({0.1, 0.2});
    CHECK(two_inputs.input_target_preview()[0] == ColVec{0.1, 0.2});
    const ColVec x2{0.4, -0.1};
    const double u_two = two_braced.step(ColVec{0.3, 0.0}, x2)(0);
    CHECK(two_mixed.step({0.3, 0.0}, x2)(0) == doctest::Approx(u_two).epsilon(1e-12));
    CHECK(two_converted.step({0.3, 0}, x2)(0) == doctest::Approx(u_two).epsilon(1e-12)); // any value that converts
}
