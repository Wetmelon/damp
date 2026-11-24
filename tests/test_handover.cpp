// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Bumpless handover: track(u, ...) in place of step(...) while another source drives the plant.
// Each tick it predicts and updates as step does, then leaves the memory that makes the next
// command start from u: the integrators preloaded, the filter and the delay states on u, an
// MPC's held input u. A Handover fades the mismatch of a law with nothing to preload.

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/law.hpp"
#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/controllers/oflqg.hpp"
#include "damp/controllers/state_feedback.hpp"
#include "damp/controllers/state_feedback_types.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/input_limits.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/toolbox/bounds.hpp"
#include "doctest.h"

using namespace damp;

namespace {

constexpr double kTs = 0.01;
constexpr double kFade = 0.2; // [s] the Handover's fade time

/// A mass on a spring and a damper: position, speed; force in; position out. At rest under a force u it sits at x = [u; 0].
constexpr StateSpace<2, 1, 1, double, 2, 1> plant() {
    return {
        .A = Matrix<2, 2>{{0.0, 1.0}, {-1.0, -0.5}},
        .B = Matrix<2, 1>{{0.0}, {1.0}},
        .C = Matrix<1, 2>{{1.0, 0.0}},
        .G = Matrix<2, 2>::identity(),
        .H = Matrix<1, 1>::identity(),
    };
}

/// The same mass pushed by two forces, the second at half the leverage: fewer integrators (one, on the position) than inputs.
constexpr StateSpace<2, 2, 1, double, 2, 1> two_forces() {
    return {
        .A = Matrix<2, 2>{{0.0, 1.0}, {-1.0, -0.5}},
        .B = Matrix<2, 2>{{0.0, 0.0}, {1.0, 0.5}},
        .C = Matrix<1, 2>{{1.0, 0.0}},
        .G = Matrix<2, 2>::identity(),
        .H = Matrix<1, 1>::identity(),
    };
}

constexpr Matrix<2, 2>                         kQ{{10.0, 0.0}, {0.0, 1.0}};
constexpr design::LQRWeights<2, 1>             kWeights{.Q = kQ, .R = 1.0};
constexpr design::KalmanNoise<2, 1>            kNoise{.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.1}}, .measurement = Matrix<1, 1>{{1e-4}}};
constexpr design::KalmanNoise<2, 1, double, 1> kWalk{.process = kNoise.process, .measurement = kNoise.measurement, .disturbance = Matrix<1, 1>{{1.0}}};
constexpr design::IntegralAction               kIntegral{.Q = 80.0};
constexpr design::Handover                     kHandover{.fade_time = kFade};

const auto kSampled = discretize(plant(), kTs).value();

/// One period of the plant under @p u, with a load @p w added to the force.
ColVec<2> advance(const ColVec<2>& x, double u, double w = 0.0) {
    return (kSampled.A * x) + (kSampled.B * ColVec<1>{u + w});
}

/// The command a copy of @p law would return now: what track left it set to.
template<typename Law, typename... Args>
ColVec<1> probe(const Law& law, const Args&... args) {
    Law copy = law;
    return copy.control(args...);
}

/// Whether @p M has a track overload taking @p A...: a regulator's take no r, a tracking law's do.
template<typename M, typename... A>
concept tracks = requires(M m, const A&... a) { m.track(a...); };

/// A NaN built from its bits at run time, past -ffinite-math-only.
double quiet_nan() {
    volatile std::uint64_t bits = 0x7ff8000000000000ULL;
    const std::uint64_t    b = bits;
    double                 out{};
    std::memcpy(&out, &b, sizeof(out));
    return out;
}

} // namespace

TEST_SUITE("Handover") {
    TEST_CASE("an LQI tracking a manual input hands over from it: its first command is u, and r moved at the handover moves it without a jump in the integral") {
        const auto res = design::state_feedback(plant(), kTs, kWeights, kIntegral);
        REQUIRE(res.success);
        constexpr double u = 0.7;
        const ColVec<1>  r{u}; // the operator holds the plant where r asks: nothing changes
        ColVec<2>        x{u, 0.0};
        LQI              law{res};
        for (int k = 0; k < 200; ++k) {
            law.track(ColVec<1>{u}, r, x);
            CHECK(probe(law, r)(0) == doctest::Approx(u).epsilon(1e-12)); // control at this state returns u
            x = advance(x, u);
        }
        const double xi = law.integral.xi(0);
        CHECK(law.step(r, x)(0) == doctest::Approx(u).epsilon(1e-12));

        //! r moved at the handover: the first command is still u (ξ integrates the error the track
        //! held, zero), and the next moves by the integrators' step of the new error, not a jump.
        LQI moved{res};
        x = ColVec<2>{u, 0.0};
        for (int k = 0; k < 50; ++k) {
            moved.track(ColVec<1>{u}, r, x);
            x = advance(x, u);
        }
        const ColVec<1> r_new{u + 0.3};
        const ColVec<1> u1 = moved.step(r_new, x);
        CHECK(u1(0) == doctest::Approx(u).epsilon(1e-12));
        CHECK(moved.integral.xi(0) == doctest::Approx(xi).epsilon(1e-12));
        x = advance(x, u1(0));
        const double    y = (kSampled.C * x)(0, 0);
        const double    xi_before = moved.integral.xi(0);
        const ColVec<1> u2 = moved.step(r_new, x);
        CHECK(moved.integral.xi(0) - xi_before == doctest::Approx(kTs * 0.3).epsilon(1e-9)); // Ts·(r_new − y), y still u
        CHECK(y == doctest::Approx(u).epsilon(1e-6));
        CHECK(damp::abs(u2(0) - u) < 0.1); // a move of the integrators' step, where the law untracked would jump
        LQI fresh{res};
        CHECK(damp::abs(fresh.step(r_new, x)(0) - u) > 0.5);
    }

    TEST_CASE("an LQI preloads its integrators at every tick of a transient, and an LQGI on its filter's estimate") {
        const auto lqi_res = design::state_feedback(plant(), kTs, kWeights, kIntegral);
        const auto lqgi_res = design::state_feedback(plant(), kTs, kWeights, kNoise, kIntegral);
        REQUIRE(lqi_res.success);
        REQUIRE(lqgi_res.success);
        LQI             lqi{lqi_res};
        LQGI            lqgi{lqgi_res};
        ColVec<2>       x{};
        const ColVec<1> r{1.0};
        for (int k = 0; k < 300; ++k) {
            const double    u = 0.5 + (0.3 * damp::sin(0.05 * static_cast<double>(k))); // an operator's hand, the plant moving under it
            const ColVec<1> y = kSampled.C * x;
            lqi.track(ColVec<1>{u}, r, x);
            lqgi.track(ColVec<1>{u}, r, y);
            CHECK(probe(lqi, r)(0) == doctest::Approx(u).epsilon(1e-12));
            CHECK(probe(lqgi, r)(0) == doctest::Approx(u).epsilon(1e-12));
            x = advance(x, u);
        }

        //! At a steady state the filter's estimate settles too, and the LQGI's first step is u.
        LQGI             settled{lqgi_res};
        constexpr double u = 0.6;
        x = ColVec<2>{u, 0.0};
        for (int k = 0; k < 3000; ++k) {
            settled.track(ColVec<1>{u}, ColVec<1>{u}, kSampled.C * x);
            x = advance(x, u);
        }
        CHECK(settled.step(ColVec<1>{u}, kSampled.C * x)(0) == doctest::Approx(u).epsilon(1e-12));
        CHECK(settled.applied()(0) == u); // the input the step's predict took: the tracked one
    }

    TEST_CASE("with fewer integrators than inputs the preload is the least-squares one, and a Handover fades what it leaves without the integrators reading it as saturation") {
        constexpr design::LQRWeights<2, 2> weights{.Q = kQ, .R = Matrix<2, 2>{{1.0, 0.0}, {0.0, 2.0}}};
        constexpr design::IntegralAction   unwinding{.Q = 80.0, .tracking_time = 0.05}; // back-calculation: an excess would move ξ at once
        const auto                         plain = design::state_feedback(two_forces(), kTs, weights, unwinding);
        const auto                         faded = design::state_feedback(two_forces(), kTs, weights, unwinding, kHandover);
        REQUIRE(plain.success);
        REQUIRE(faded.success);
        const auto      sampled = discretize(two_forces(), kTs).value();
        const ColVec<2> u{0.9, -0.4}; // not a command the law's one integrator reaches
        const ColVec<1> r{0.5};
        ColVec<2>       x{0.7, 0.0};
        LQI             ls{plain};
        LQI             fade{faded};
        for (int k = 0; k < 20; ++k) {
            ls.track(u, r, x);
            fade.track(u, r, x);
            x = (sampled.A * x) + (sampled.B * u);
        }
        //! The command nearest u: what is left of u is orthogonal to the integrator's column of −K.
        LQI             copy = ls;
        const ColVec<2> near = copy.control(r);
        const ColVec<2> rest = u - near;
        const auto      Ki = plain.integral_gain();
        CHECK(damp::abs((rest(0) * Ki(0, 0)) + (rest(1) * Ki(1, 0))) < 1e-12);
        CHECK(damp::abs(rest(0)) + damp::abs(rest(1)) > 0.1);

        //! With the Handover the first command is u, and the fade is a feedforward to the
        //! integrators: they integrate the error and nothing else, tick after tick.
        const ColVec<2> first = fade.step(r, x);
        CHECK(first(0) == doctest::Approx(u(0)).epsilon(1e-12));
        CHECK(first(1) == doctest::Approx(u(1)).epsilon(1e-12));
        ColVec<2> cmd = first;
        for (int k = 0; k < 30; ++k) {
            x = (sampled.A * x) + (sampled.B * cmd);
            const double xi = fade.integral.xi(0);
            const double e = r(0) - (sampled.C * x)(0, 0);
            const double held = fade.integral.e(0);
            cmd = fade.step(r, x);
            CHECK(fade.integral.xi(0) - xi == doctest::Approx(kTs * held).epsilon(1e-9));
            CHECK(fade.integral.e(0) == doctest::Approx(e).epsilon(1e-12));
        }
    }

    TEST_CASE("an LQR with a Handover starts from the tracked input and fades the mismatch at the given rate; without one it is a plain step") {
        const auto res = design::state_feedback(plant(), kTs, kWeights);
        const auto faded = design::state_feedback(plant(), kTs, kWeights, kHandover);
        REQUIRE(res.success);
        REQUIRE(faded.success);
        LQR       law{faded};
        LQR       plain{res};
        LQR       twin{res}; // never tracked: the law's own command at each state
        ColVec<2> x{};
        for (int k = 0; k < 100; ++k) {
            const double u = 0.8 * damp::cos(0.03 * static_cast<double>(k));
            law.track(ColVec<1>{u}, x);
            plain.track(ColVec<1>{u}, x);
            x = advance(x, u);
        }
        const double tracked = 0.8 * damp::cos(0.03 * 99);
        CHECK(plain.step(x)(0) == doctest::Approx(twin.step(x)(0)).epsilon(1e-14)); // nothing to preload: the plain step
        ColVec<1> cmd = law.step(x);
        CHECK(cmd(0) == doctest::Approx(tracked).epsilon(1e-12));
        const double decay = std::exp(-kTs / kFade);
        double       mismatch = cmd(0) + (res.K * x)(0, 0);
        REQUIRE(damp::abs(mismatch) > 0.1);
        for (int k = 0; k < 40; ++k) {
            x = advance(x, cmd(0));
            cmd = law.step(x);
            const double next = cmd(0) - twin.step(x)(0);
            CHECK(next == doctest::Approx(mismatch * decay).epsilon(1e-9));
            mismatch = next;
        }

        //! Never tracked, the Handover adds nothing; after reset() the next track starts over.
        LQR       idle{faded};
        LQR       own{res};
        ColVec<2> z{0.3, -0.1};
        for (int k = 0; k < 20; ++k) {
            const ColVec<1> a = idle.step(z);
            CHECK(a(0) == doctest::Approx(own.step(z)(0)).epsilon(1e-14));
            z = advance(z, a(0));
        }
        idle.track(ColVec<1>{2.0}, z);
        idle.reset();
        CHECK(idle.step(z)(0) == doctest::Approx(own.step(z)(0)).epsilon(1e-14));
    }

    TEST_CASE("a tracked input outside the input limits is one the law cannot emit: it hands over from u clamped") {
        constexpr design::InputLimits bound{.lower = {-1.0}, .upper = {1.0}};
        const auto                    lqi = design::state_feedback(plant(), kTs, kWeights, kIntegral, bound);
        const auto                    lqr = design::state_feedback(plant(), kTs, kWeights, bound, kHandover);
        REQUIRE(lqi.success);
        REQUIRE(lqr.success);
        constexpr double u = 1.5; // an operator past the law's limit, the plant at rest under it
        const ColVec<1>  r{u};
        ColVec<2>        x{u, 0.0};
        LQI              integral{lqi};
        LQR              fade{lqr};
        for (int k = 0; k < 50; ++k) {
            integral.track(ColVec<1>{u}, r, x);
            fade.track(ColVec<1>{u}, x);
            CHECK(probe(integral, r)(0) == doctest::Approx(1.0).epsilon(1e-12));
            x = advance(x, u);
        }
        const double xi = integral.integral.xi(0);
        CHECK(integral.step(r, x)(0) == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(integral.integral.xi(0) == doctest::Approx(xi).epsilon(1e-12)); // the rest of u is the operator's: no excess to hold or unwind
        CHECK(fade.step(x)(0) == doctest::Approx(1.0).epsilon(1e-12));
        for (int k = 0; k < 30; ++k) {
            x = advance(x, 1.0);
            CHECK(damp::abs(fade.step(x)(0)) <= 1.0); // the sum, clamped
        }
    }

    TEST_CASE("an OFLQG and an OFMPC tracking an operator against a load learn it into d-hat, and their first command is u") {
        constexpr double w = -0.25; // a load on the force
        constexpr double y_r = 0.6;
        constexpr double u = y_r - w; // the force that holds y at r against it
        const auto       of = design::state_feedback(plant(), kTs, kWeights, kWalk, design::InputDisturbance{});
        const auto       ofmpc = design::mpc<10, 3>(plant(), kTs, design::MPCWeights<1, 1>{}, kWalk, design::InputDisturbance{});
        REQUIRE(of.success);
        REQUIRE(ofmpc.success);
        OFLQG     lq{of};
        OFMPC     mpc{ofmpc};
        ColVec<2> x{y_r, 0.0};
        for (int k = 0; k < 3000; ++k) {
            const ColVec<1> y = kSampled.C * x;
            lq.track(ColVec<1>{u}, ColVec<1>{y_r}, y);
            mpc.track(ColVec<1>{u}, ColVec<1>{y_r}, y);
            x = advance(x, u, w);
        }
        CHECK(lq.disturbance_estimate()(0) == doctest::Approx(w).epsilon(1e-9));
        CHECK(mpc.disturbance_estimate()(0) == doctest::Approx(w).epsilon(1e-9));
        const ColVec<1> y = kSampled.C * x;
        CHECK(lq.step(ColVec<1>{y_r}, y)(0) == doctest::Approx(u).epsilon(1e-9));
        CHECK(mpc.step(ColVec<1>{y_r}, y)(0) == doctest::Approx(u).epsilon(1e-9));
        CHECK(mpc.applied()(0) == u);
    }

    TEST_CASE("an MPC holds the tracked input as u_prev: Rdu and the move box shape its first move away from it") {
        design::MPCWeights<1, 1> light{};
        design::MPCWeights<1, 1> heavy{};
        light.Rdu = Matrix<1, 1>{{1e-3}};
        heavy.Rdu = Matrix<1, 1>{{1.0}};
        constexpr design::InputRateLimits move{.lower = {-0.02}, .upper = {0.02}};
        const auto                        a = design::mpc<10, 3>(plant(), kTs, light, design::TrackedOutputs{});
        const auto                        b = design::mpc<10, 3>(plant(), kTs, heavy, design::TrackedOutputs{});
        const auto                        c = design::mpc<10, 3>(plant(), kTs, light, move, design::TrackedOutputs{});
        REQUIRE(a.success);
        REQUIRE(b.success);
        REQUIRE(c.success);
        MPC              fast{a};
        MPC              slow{b};
        MPC              boxed{c};
        MPC              untracked{a};
        constexpr double u = 0.5; // the operator pushes the mass off r = 0
        ColVec<2>        x{};
        for (int k = 0; k < 30; ++k) {
            fast.track(ColVec<1>{u}, ColVec<1>{0.0}, x);
            slow.track(ColVec<1>{u}, ColVec<1>{0.0}, x);
            boxed.track(ColVec<1>{u}, ColVec<1>{0.0}, x);
            x = advance(x, u);
        }
        const double first_fast = fast.step(ColVec<1>{0.0}, x)(0);
        const double first_slow = slow.step(ColVec<1>{0.0}, x)(0);
        const double first_boxed = boxed.step(ColVec<1>{0.0}, x)(0);
        CHECK(fast.applied()(0) == u); // the base of the first move
        CHECK(damp::abs(first_slow - u) < damp::abs(first_fast - u));
        CHECK(damp::abs(first_boxed - u) <= 0.02 + 1e-12);
        CHECK(damp::abs(untracked.step(ColVec<1>{0.0}, x)(0) - first_fast) > 0.1); // from u_prev = 0 instead
    }

    TEST_CASE("with an InputRateLimits the first command after a track moves at most the rate from the tracked u, on an LQ law as on an MPC") {
        //! The base of the move is the input applied: an LQR holds it alone, an LQG in its record,
        //! an LQI there beside its preload, an MPC as its held input. The operator holds the mass at
        //! rest off r = 0, where each regulator would move far from u on its own.
        constexpr double                     u = 0.5;
        constexpr double                     rate = 0.01;
        constexpr double                     slow = 1e-4; // the MPC's: its own first move is gentler than the LQR's
        constexpr design::InputRateLimits<1> moves{.lower = {-rate}, .upper = {rate}};
        constexpr design::InputRateLimits<1> slow_moves{.lower = {-slow}, .upper = {slow}};
        design::MPCWeights<1, 1>             light{};
        light.Rdu = Matrix<1, 1>{{1e-3}};
        const auto lqr = design::state_feedback(plant(), kTs, kWeights, moves);
        const auto lqg = design::state_feedback(plant(), kTs, kWeights, kNoise, moves);
        const auto mpc = design::mpc<10, 3>(plant(), kTs, light, slow_moves, design::TrackedOutputs{});
        REQUIRE(lqr.success);
        REQUIRE(lqg.success);
        REQUIRE(mpc.success);
        LQR       a{lqr};
        LQG       b{lqg};
        MPC       c{mpc};
        LQR       a_fresh{lqr}; // never tracked: the first move is from 0, as an MPC's
        MPC       c_fresh{mpc};
        ColVec<2> x{u, 0.0};
        for (int k = 0; k < 400; ++k) {
            a.track(ColVec<1>{u}, x);
            b.track(ColVec<1>{u}, kSampled.C * x);
            c.track(ColVec<1>{u}, ColVec<1>{0.0}, x);
            x = advance(x, u);
        }
        REQUIRE(-(lqr.K * x)(0, 0) < u - (2.0 * rate)); // the law's own command is far below u: the clamp binds
        CHECK(a.step(x)(0) == doctest::Approx(u - rate).epsilon(1e-12));
        CHECK(b.step(kSampled.C * x)(0) == doctest::Approx(u - rate).epsilon(1e-12));
        CHECK(c.step(ColVec<1>{0.0}, x)(0) == doctest::Approx(u - slow).epsilon(1e-12));
        CHECK(a_fresh.step(x)(0) == doctest::Approx(-rate).epsilon(1e-12));
        CHECK(c_fresh.step(ColVec<1>{0.0}, x)(0) == doctest::Approx(-slow).epsilon(1e-12));

        //! An LQI with fewer integrators than inputs: the preload leaves the command short of u, and
        //! the move is still from u, not from that command.
        constexpr design::LQRWeights<2, 2>   weights{.Q = kQ, .R = Matrix<2, 2>{{1.0, 0.0}, {0.0, 2.0}}};
        constexpr double                     tight = 1e-3;
        constexpr design::InputRateLimits<2> creep{.lower = {-tight, -tight}, .upper = {tight, tight}};
        const auto                           lqi = design::state_feedback(two_forces(), kTs, weights, kIntegral, creep);
        const auto                           free = design::state_feedback(two_forces(), kTs, weights, kIntegral); // the same preload, no clamp on the move
        REQUIRE(lqi.success);
        REQUIRE(free.success);
        const auto      sampled = discretize(two_forces(), kTs).value();
        const ColVec<2> v{0.9, -0.4}; // not a command the law's one integrator reaches
        const ColVec<1> r{0.5};
        ColVec<2>       z{0.7, 0.0};
        LQI             d{lqi};
        LQI             e{free};
        for (int k = 0; k < 20; ++k) {
            d.track(v, r, z);
            e.track(v, r, z);
            z = (sampled.A * z) + (sampled.B * v);
        }
        LQI             copy = e;
        const ColVec<2> near = copy.control(r);           // the command the preload reaches
        REQUIRE(damp::abs(near(0) - v(0)) > 2.0 * tight); // short of u by more than a move
        const ColVec<2> first = d.step(r, z);
        CHECK(damp::abs(first(0) - v(0)) <= tight + 1e-12);
        CHECK(damp::abs(first(1) - v(1)) <= tight + 1e-12);
    }

    TEST_CASE("an MPC regulator tracks as the LQR does, track(u, x) and its d and feedforward forms, and hands over as TrackedOutputs{}'s MPC at r = 0") {
        //! A design with no TrackedOutputs takes no r: each regulator step has its track twin, under
        //! the same requires clause. Its commands are those of the MPC that holds every output at
        //! r = 0, so after the same track the two hand over alike, both from u_prev = u.
        using Reg = std::remove_cvref_t<decltype(MPC{design::mpc<10, 3>(plant(), kTs, design::MPCWeights<1, 1>{})})>;
        using Every = std::remove_cvref_t<decltype(MPC{design::mpc<10, 3>(plant(), kTs, design::MPCWeights<1, 1>{}, design::TrackedOutputs{})})>;
        static_assert(tracks<Reg, ColVec<1>, ColVec<2>> && tracks<Reg, ColVec<1>, ColVec<2>, ColVec<1>, ColVec<2>>);
        static_assert(!tracks<Reg, ColVec<1>, ColVec<1>, ColVec<2>>); // no r to take
        static_assert(!tracks<Every, ColVec<1>, ColVec<2>>);          // it takes r
        static_assert(tracks<Every, ColVec<1>, ColVec<1>, ColVec<2>>);

        design::MPCWeights<1, 1> light{};
        light.Rdu = Matrix<1, 1>{{1e-3}};
        constexpr double                  u = 0.4;
        constexpr double                  w = 0.2;
        const design::MeasuredDisturbance load{.Bd = plant().B};
        const ColVec<1>                   r0{0.0};
        const ColVec<1>                   d{w};
        const ColVec<2>                   x_ref{0.3, 0.0};
        const ColVec<1>                   u_ref{0.3};
        const auto                        same = [](const ColVec<1>& a, const ColVec<1>& b) { CHECK(a(0) == doctest::Approx(b(0)).epsilon(1e-12)); };
        SUBCASE("track(u, x), and the feedforward form track(u, x_ref, u_ref, x)") {
            const auto reg = design::mpc<10, 3>(plant(), kTs, light);
            const auto every = design::mpc<10, 3>(plant(), kTs, light, design::TrackedOutputs{});
            REQUIRE(reg.success);
            REQUIRE(every.success);
            MPC       a{reg};
            MPC       b{every};
            MPC       af{reg};
            MPC       bf{every};
            MPC       untracked{reg};
            ColVec<2> x{u, 0.0};
            for (int k = 0; k < 30; ++k) {
                a.track(ColVec<1>{u}, x);
                b.track(ColVec<1>{u}, r0, x);
                af.track(ColVec<1>{u}, x_ref, u_ref, x);
                bf.track(ColVec<1>{u}, r0, x_ref, u_ref, x);
                x = advance(x, u);
            }
            const ColVec<1> first = a.step(x);
            CHECK(a.applied()(0) == u); // the base of the first move
            same(first, b.step(r0, x));
            same(af.step(x_ref, u_ref, x), bf.step(r0, x_ref, u_ref, x));
            CHECK(damp::abs(untracked.step(x)(0) - first(0)) > 0.1); // from u_prev = 0 instead
        }
        SUBCASE("under a measured d: track(u, d, x) and track(u, x_ref, u_ref, d, x)") {
            const auto reg = design::mpc<10, 3>(plant(), kTs, light, load);
            const auto every = design::mpc<10, 3>(plant(), kTs, light, design::TrackedOutputs{}, load);
            REQUIRE(reg.success);
            REQUIRE(every.success);
            using RegD = std::remove_cvref_t<decltype(MPC{reg})>;
            static_assert(tracks<RegD, ColVec<1>, ColVec<1>, ColVec<2>> && tracks<RegD, ColVec<1>, ColVec<2>, ColVec<1>, ColVec<1>, ColVec<2>>);
            static_assert(!tracks<RegD, ColVec<1>, ColVec<2>>); // d is not optional
            MPC       a{reg};
            MPC       b{every};
            MPC       af{reg};
            MPC       bf{every};
            ColVec<2> x{u + w, 0.0};
            for (int k = 0; k < 30; ++k) {
                a.track(ColVec<1>{u}, d, x);
                b.track(ColVec<1>{u}, r0, d, x);
                af.track(ColVec<1>{u}, x_ref, u_ref, d, x);
                bf.track(ColVec<1>{u}, r0, x_ref, u_ref, d, x);
                x = advance(x, u, w);
            }
            same(a.step(d, x), b.step(r0, d, x));
            CHECK(a.applied()(0) == u);
            same(af.step(x_ref, u_ref, d, x), bf.step(r0, x_ref, u_ref, d, x));
        }
    }

    TEST_CASE("the delay forms hold the tracked input in their delay states, the d forms feed d to the filter, the feedforward forms take the planned point") {
        constexpr double u = 0.4;
        SUBCASE("an LQI and an MPC whose input arrives late") {
            const auto lqi = design::state_feedback(plant(), kTs, kWeights, kIntegral, InputDelay<2>{0.015});
            const auto mpc = design::mpc<10, 3>(plant(), kTs, design::MPCWeights<1, 1>{}, InputDelay<2>{0.015}, design::TrackedOutputs{});
            REQUIRE(lqi.success);
            REQUIRE(mpc.success);
            LQI       a{lqi};
            MPC       b{mpc};
            ColVec<2> x{u, 0.0};
            for (int k = 0; k < 10; ++k) {
                a.track(ColVec<1>{u}, ColVec<1>{u}, x);
                b.track(ColVec<1>{u}, ColVec<1>{u}, x);
                x = advance(x, u);
            }
            CHECK(a.step(ColVec<1>{u}, x)(0) == doctest::Approx(u).epsilon(1e-12)); // the delay states at u: no command on its way to undo
            CHECK(b.step(ColVec<1>{u}, x)(0) == doctest::Approx(u).epsilon(1e-9));
        }
        SUBCASE("an LQI and an LQGI that read the load") {
            constexpr double                  w = 0.2;
            const design::MeasuredDisturbance load{.Bd = plant().B};
            const auto                        lqi = design::state_feedback(plant(), kTs, kWeights, kIntegral, load);
            const auto                        lqgi = design::state_feedback(plant(), kTs, kWeights, kNoise, kIntegral, load);
            REQUIRE(lqi.success);
            REQUIRE(lqgi.success);
            LQI             a{lqi};
            LQGI            b{lqgi};
            const ColVec<1> d{w};
            const ColVec<1> r{u + w};
            ColVec<2>       x{u + w, 0.0};
            for (int k = 0; k < 3000; ++k) {
                a.track(ColVec<1>{u}, r, d, x);
                b.track(ColVec<1>{u}, r, d, kSampled.C * x);
                CHECK(probe(a, r, d)(0) == doctest::Approx(u).epsilon(1e-12));
                x = advance(x, u, w);
            }
            CHECK(a.step(r, d, x)(0) == doctest::Approx(u).epsilon(1e-12));
            CHECK(b.step(r, d, kSampled.C * x)(0) == doctest::Approx(u).epsilon(1e-12)); // the filter predicted with [u; d]
        }
        SUBCASE("an LQI, an LQR with a Handover and an MPC about a planned point") {
            const auto      lqi = design::state_feedback(plant(), kTs, kWeights, kIntegral);
            const auto      lqr = design::state_feedback(plant(), kTs, kWeights, kHandover);
            const auto      mpc = design::mpc<10, 3>(plant(), kTs, design::MPCWeights<1, 1>{}, design::TrackedOutputs{});
            LQI             a{lqi};
            LQR             b{lqr};
            MPC             c{mpc};
            const ColVec<2> x_ref{0.3, 0.0};
            const ColVec<1> u_ref{0.3};
            ColVec<2>       x{0.1, 0.0};
            for (int k = 0; k < 20; ++k) {
                a.track(ColVec<1>{u}, ColVec<1>{0.0}, x_ref, u_ref, x);
                b.track(ColVec<1>{u}, x_ref, u_ref, x);
                c.track(ColVec<1>{u}, ColVec<1>{0.0}, x_ref, u_ref, x);
                CHECK(probe(a, ColVec<1>{0.0}, x_ref, u_ref)(0) == doctest::Approx(u).epsilon(1e-12));
                x = advance(x, u);
            }
            CHECK(b.step(x_ref, u_ref, x)(0) == doctest::Approx(u).epsilon(1e-12));
            static_cast<void>(c.step(ColVec<1>{0.0}, x_ref, u_ref, x));
            CHECK(c.applied()(0) == u);
        }
    }

    TEST_CASE("Law<kRes, float> tracks on constants as on members, and the same law with members agrees") {
        static constexpr auto kLqi = design::state_feedback(plant(), kTs, kWeights, kIntegral);
        static constexpr auto kLqrFaded = design::state_feedback(plant(), kTs, kWeights, kHandover);
        static constexpr auto kLqgFaded = design::state_feedback(plant(), kTs, kWeights, kNoise, kHandover);
        static constexpr auto kOfmpc = design::mpc<6, 2>(plant(), kTs, design::MPCWeights<1, 1>{}, kWalk, design::InputDisturbance{});
        static const auto     kAtInit = design::state_feedback(plant(), kTs, kWeights, kNoise, kIntegral);
        static_assert(detail::constant_design<kLqrFaded> && !detail::constant_design<kAtInit>);

        Law<kLqi, float>       lqi;
        Law<kLqrFaded, float>  lqr;
        Law<kLqgFaded, float>  lqg;
        Law<kOfmpc, float>     ofmpc;
        Law<kAtInit, float>    lqgi;
        LQI                    lqi_held{kLqi.as<float>()};
        LQR                    lqr_held{kLqrFaded.as<float>()};
        constexpr float        u = 0.4f;
        const ColVec<1, float> r{u};
        ColVec<2>              x{0.4, 0.0};
        for (int k = 0; k < 400; ++k) {
            const ColVec<2, float> xf = x.as<float>();
            const ColVec<1, float> y{xf(0)};
            lqi.track(ColVec<1, float>{u}, r, xf);
            lqi_held.track(ColVec<1, float>{u}, r, xf);
            lqr.track(ColVec<1, float>{u}, xf);
            lqr_held.track(ColVec<1, float>{u}, xf);
            lqg.track(ColVec<1, float>{u}, y);
            ofmpc.track(ColVec<1, float>{u}, r, y);
            lqgi.track(ColVec<1, float>{u}, r, y);
            x = advance(x, static_cast<double>(u));
        }
        const ColVec<2, float> xf = x.as<float>();
        const ColVec<1, float> y{xf(0)};
        constexpr float        tol = 64.0f * std::numeric_limits<float>::epsilon();
        const auto             near = [](float a, float b, float rel) { return damp::abs(a - b) <= rel * (1.0f + damp::max(damp::abs(a), damp::abs(b))); }; // as doctest::Approx
        CHECK(near(lqi.step(r, xf)(0), u, tol));
        CHECK(near(lqi_held.step(r, xf)(0), u, tol));
        CHECK(near(lqgi.step(r, y)(0), u, tol));
        CHECK(near(lqg.step(y)(0), u, tol));
        CHECK(near(ofmpc.step(r, y)(0), u, 1e-3f)); // d-hat learned in float
        CHECK(near(lqr.step(xf)(0), u, tol));
        CHECK(near(lqr_held.step(xf)(0), u, tol));
        for (int k = 0; k < 20; ++k) {
            CHECK(near(lqr.step(xf)(0), lqr_held.step(xf)(0), tol)); // the factor a constant of the type, or a member
        }
    }

    TEST_CASE("track adds nothing to a law without a Handover; the Handover holds the mismatch, its factor and a flag") {
        static constexpr auto kLqr = design::state_feedback(plant(), kTs, kWeights);
        static constexpr auto kLqi = design::state_feedback(plant(), kTs, kWeights, kIntegral);
        static constexpr auto kLqrFaded = design::state_feedback(plant(), kTs, kWeights, kHandover);
        static constexpr auto kLqiFaded = design::state_feedback(plant(), kTs, kWeights, kIntegral, kHandover);
        static_assert(std::is_same_v<decltype(LQR{kLqr}), LQR<2, 1, double, 2>>);
        static_assert(std::is_same_v<decltype(LQI{kLqi}), LQI<2, 1, 1, double, 1>>);
        static_assert(std::is_same_v<decltype(LQR{kLqrFaded}), LQR<2, 1, double, 2, 0, detail::HeldGains, 0, 0, 0, true>>);
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(kLqr)>, design::StateFeedbackResult<2, 1, double>>);
        static_assert(sizeof(kLqrFaded) == sizeof(kLqr) + sizeof(double)); // the fade time
        using Fade = detail::HandoverFade<1, float, true>;
        using ConstantFade = detail::HandoverFade<1, float, false>;
        static_assert(sizeof(ConstantFade) < sizeof(Fade)); // on a constexpr design's Law the factor is the type's
        if constexpr (sizeof(void*) == 8) {
            CHECK(sizeof(LQR{kLqr.as<float>()}) == 48); // test_lq_limits records these without track
            CHECK(sizeof(LQI{kLqi.as<float>()}) == 112);
            CHECK(sizeof(Law<kLqr, float>) == 20);
            CHECK(sizeof(Law<kLqi, float>) == 68);
            //! With a Handover: the mismatch and the factor (on constants the mismatch alone), the flag in the padding.
            CHECK(sizeof(LQR{kLqrFaded.as<float>()}) == 48 + (2 * sizeof(float)));
            CHECK(sizeof(LQI{kLqiFaded.as<float>()}) == 112 + (2 * sizeof(float)));
            CHECK(sizeof(Law<kLqrFaded, float>) == 20 + sizeof(float));
            CHECK(sizeof(Law<kLqiFaded, float>) == 68 + sizeof(float));
        }
    }

    TEST_CASE("a Handover's fade time is finite and positive, else the design is refused") {
        CHECK(design::state_feedback(plant(), kTs, kWeights, design::Handover{.fade_time = 0.0}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::state_feedback(plant(), kTs, kWeights, design::Handover{.fade_time = -0.1}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::state_feedback(plant(), kTs, kWeights, design::Handover{.fade_time = quiet_nan()}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::state_feedback(plant(), kTs, kWeights, design::Handover{.fade_time = std::numeric_limits<double>::infinity()}).status == design::DesignStatus::invalid_parameter);
        const auto ok = design::state_feedback(plant(), kTs, kWeights, kHandover);
        CHECK(ok.success);
        CHECK(ok.fade_time == kFade);
        CHECK(ok.K == design::state_feedback(plant(), kTs, kWeights).K); // the gain is the same
        CHECK(design::state_feedback(plant(), kTs, kWeights, design::Handover{.fade_time = 0.0}, InputDelay<1>{0.005}).status == design::DesignStatus::invalid_parameter);
    }
}
