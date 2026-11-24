// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// The limits of the LQ laws: the InputLimits clamp the command control returns, an InputRateLimits its
// move, the IntegralAction's limits the integral state. The oracle for each is the same design without limits, clamped outside the law and
// told with back_calculate: the two are one tick, on every runtime, so the integrators hold against
// the law's own clamp and a filter or the delay states move on the input that went out.

#include <cmath>
#include <cstddef>
#include <limits>

#include "damp/backend.hpp"
#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/controllers/oflqg.hpp"
#include "damp/controllers/state_feedback.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/input_limits.hpp"
#include "damp/design/qp.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/toolbox/bounds.hpp"
#include "doctest.h"

using namespace damp;

namespace {

constexpr double kTs = 0.01;
constexpr double kLimit = 1.05; // just above the input that holds r = 1 (k·r = 1)

/// A mass on a spring and a damper: position, speed; force in; position out. Noise on both states and the reading.
constexpr StateSpace<2, 1, 1, double, 2, 1> plant() {
    return {
        .A = Matrix<2, 2>{{0.0, 1.0}, {-1.0, -0.5}},
        .B = Matrix<2, 1>{{0.0}, {1.0}},
        .C = Matrix<1, 2>{{1.0, 0.0}},
        .G = Matrix<2, 2>::identity(),
        .H = Matrix<1, 1>::identity(),
    };
}

constexpr Matrix<2, 2>              Q{{10.0, 0.0}, {0.0, 1.0}};
constexpr Bounds<1>                 kBox{.lower = -kLimit, .upper = kLimit};
constexpr design::KalmanNoise<2, 1> kNoise{.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.1}}, .measurement = Matrix<1, 1>{{1e-4}}};
constexpr design::LQRWeights<2, 1>  kFree{.Q = Q, .R = 1.0};
constexpr design::InputLimits       kBound{.lower = kBox.lower, .upper = kBox.upper};
constexpr design::IntegralAction    kIntegral{.Q = 80.0}; // on kFree, with kBound or without, the LQI of each
constexpr double                    kStep = 0.002;        // a move limit that binds on a step of r: 0.2 a second
constexpr design::InputRateLimits   kRate{.lower = {-kStep}, .upper = {kStep}};

// The weights have no integrators to weight or limit, and no limits; the integral action carries its own.
template<typename C>
concept integral_fields = requires(C c) {
    c.Q_integral;
    c.integral;
};
template<typename C>
concept integral_action_fields = requires(C c) {
    c.Q;
    c.limits;
    c.tracking_time;
};
static_assert(!integral_fields<design::LQRWeights<2, 1>> && !integral_action_fields<design::LQRWeights<2, 1>>);
static_assert(integral_action_fields<design::IntegralAction<double>>);

/**
 * Run @p bounded (its design carries the InputLimits) and @p clamped (the same design without, clamped
 * outside and told with back_calculate), each on its own copy of @p model from rest.
 * @return the largest gap between their commands, relative to the limit; @p saturated counts the
 *         ticks the bounded law sat on its limit, @p y_end its output at the end
 */
template<size_t NXM, typename T, typename A, typename B, typename Feed, typename Control>
T gap(const StateSpace<NXM, 1, 1, T>& model, A& bounded, B& clamped, Feed feed, Control control, int& saturated, T& y_end) {
    const T        lim = static_cast<T>(kLimit);
    ColVec<NXM, T> za{};
    ColVec<NXM, T> zb{};
    T              worst{0};
    saturated = 0;
    for (int k = 0; k < 1500; ++k) {
        feed(bounded, za);
        feed(clamped, zb);
        const ColVec<1, T> ua = control(bounded);
        const ColVec<1, T> ub{damp::clamp(control(clamped)(0), -lim, lim)};
        clamped.back_calculate(ub);
        bounded.predict();
        clamped.predict();
        saturated += (damp::abs(ua(0)) >= lim) ? 1 : 0;
        worst = damp::max(worst, damp::abs(ua(0) - ub(0)));
        CHECK(damp::abs(ua(0)) <= lim);
        za = ColVec<NXM, T>((model.A * za) + (model.B * ua));
        zb = ColVec<NXM, T>((model.A * zb) + (model.B * ub));
    }
    y_end = (model.C * za)(0, 0);
    return worst / lim;
}

/**
 * Run @p rated (its design carries kRate) and @p outside (the same design without, its command's
 * move clamped outside and told with back_calculate), each on its own copy of @p model from rest,
 * toward r = 1.
 * @return the largest gap between their commands; @p bound counts the ticks the clamp bound, @p y_end
 *         the rated loop's output at the end
 */
template<size_t NXM, typename T, typename A, typename B, typename Feed>
T slew_gap(const StateSpace<NXM, 1, 1, T>& model, A& rated, B& outside, Feed feed, int& bound, T& y_end) {
    const T        step = static_cast<T>(kStep);
    ColVec<NXM, T> za{};
    ColVec<NXM, T> zb{};
    T              prev{0}; // the base of the first move: 0, as an MPC's held input
    T              last{0}; // the rated law's own
    T              worst{0};
    bound = 0;
    for (int k = 0; k < 3000; ++k) {
        feed(rated, za);
        feed(outside, zb);
        const ColVec<1, T> ua = rated.control(ColVec<1, T>{T{1}});
        const T            asked = outside.control(ColVec<1, T>{T{1}})(0);
        const ColVec<1, T> ub{damp::clamp(asked, prev - step, prev + step)};
        outside.back_calculate(ub);
        rated.predict();
        outside.predict();
        bound += (ub(0) != asked) ? 1 : 0;
        CHECK(damp::abs(ua(0) - last) <= step + (T{4} * eps<T>() * (damp::abs(last) + step))); // never more than the limit from its last command, to the rounding of the sum
        last = ua(0);
        worst = damp::max(worst, damp::abs(ua(0) - ub(0)));
        prev = ub(0);
        za = ColVec<NXM, T>((model.A * za) + (model.B * ua));
        zb = ColVec<NXM, T>((model.A * zb) + (model.B * ub));
    }
    y_end = (model.C * za)(0, 0);
    return worst;
}

/// A few hundred ticks of a closed loop keep a rounding difference near its size: well under this.
template<typename T>
constexpr T tol = T{1024} * eps<T>();

template<typename T>
StateSpace<2, 1, 1, T> sampled() {
    const auto d = discretize(plant(), kTs).value();
    return {.A = d.A.template as<T>(), .B = d.B.template as<T>(), .C = d.C.template as<T>(), .Ts = static_cast<T>(kTs)};
}

} // namespace

TEST_SUITE("LQ limits") {
    TEST_CASE_TEMPLATE("LQR: control returns the command clamped into its InputLimits, and no limit is no clamp", T, float, double) {
        const auto free = design::state_feedback(plant(), kTs, kFree);
        const auto bound = design::state_feedback(plant(), kTs, kFree, kBound);
        REQUIRE(free.success);
        REQUIRE(bound.success);
        CHECK((bound.K - free.K).norm() <= 4.0 * eps<double>() * free.K.norm()); // a limit is the runtime's: the gain is the same
        LQR                law{bound.template as<T>()};
        LQR                open{free.template as<T>()};
        const ColVec<2, T> far{T{100}, T{-40}};
        const ColVec<2, T> near{static_cast<T>(0.01), T{0}};
        law.update(far);
        open.update(far);
        const T u_far = open.control()(0);
        REQUIRE(damp::abs(u_far) > static_cast<T>(kLimit));
        CHECK(law.control()(0) == ((u_far > T{0}) ? static_cast<T>(kLimit) : -static_cast<T>(kLimit)));
        law.update(near);
        open.update(near);
        const T u_near = open.control()(0);
        CHECK(damp::abs(law.control()(0) - u_near) <= T{4} * eps<T>() * damp::abs(u_near)); // inside the box: the law itself

        // A rail that moves is set at run time.
        auto lim = law.limits();
        CHECK(lim.input.upper(0) == static_cast<T>(kLimit));
        lim.input = {.lower = T{-5}, .upper = T{5}};
        REQUIRE(law.set_limits(lim));
        CHECK(law.limits().input.upper(0) == T{5});
        law.update(far);
        CHECK(damp::abs(law.control()(0)) == T{5});
    }

    TEST_CASE_TEMPLATE("LQI with InputLimits is the law clamped outside with back_calculate: the integrator holds against its own clamp", T, float, double) {
        const auto bound = design::state_feedback(plant(), kTs, kFree, kBound, kIntegral);
        const auto free = design::state_feedback(plant(), kTs, kFree, kIntegral);
        REQUIRE(bound.success);
        LQI  bounded{bound.template as<T>()};
        LQI  clamped{free.template as<T>()};
        int  saturated = 0;
        T    y{0};
        auto feed = [](auto& law, const ColVec<2, T>& z) { law.update(z); };
        auto control = [](auto& law) { return law.control(ColVec<1, T>{T{1}}); };
        CHECK(gap(sampled<T>(), bounded, clamped, feed, control, saturated, y) <= tol<T>);
        CHECK(saturated > 10); // the step drives it onto the limit
        CHECK(static_cast<double>(y) == doctest::Approx(1.0).epsilon(1e-3));
    }

    TEST_CASE_TEMPLATE("LQI with an InputDelay: the delay states hold the clamped command", T, float, double) {
        const InputDelay<2> late{0.014};
        const auto          bound = design::state_feedback(plant(), kTs, kFree, kBound, late, kIntegral);
        const auto          free = design::state_feedback(plant(), kTs, kFree, late, kIntegral);
        REQUIRE(bound.success);
        const auto       d = discretize(plant(), kTs, late).value();
        const StateSpace model{.A = d.A.template as<T>(), .B = d.B.template as<T>(), .C = d.C.template as<T>(), .Ts = static_cast<T>(kTs)};
        LQI              bounded{bound.template as<T>()};
        LQI              clamped{free.template as<T>()};
        int              saturated = 0;
        T                y{0};
        auto             feed = [](auto& law, const ColVec<4, T>& z) { law.update(ColVec<2, T>{z(0), z(1)}); };
        auto             control = [](auto& law) { return law.control(ColVec<1, T>{T{1}}); };
        CHECK(gap(model, bounded, clamped, feed, control, saturated, y) <= tol<T>);
        CHECK(saturated > 10);
        CHECK(static_cast<double>(y) == doctest::Approx(1.0).epsilon(1e-3));
    }

    TEST_CASE_TEMPLATE("LQGI and OFLQG with InputLimits: the filter predicts with the clamped command", T, float, double) {
        auto feed = [](auto& law, const ColVec<2, T>& z) { law.update(ColVec<1, T>{z(0)}); };
        auto control = [](auto& law) { return law.control(ColVec<1, T>{T{1}}); };
        int  saturated = 0;
        T    y{0};

        const auto bound = design::state_feedback(plant(), kTs, kFree, kBound, kNoise, kIntegral);
        const auto free = design::state_feedback(plant(), kTs, kFree, kNoise, kIntegral);
        REQUIRE(bound.success);
        LQGI bounded{bound.template as<T>()};
        LQGI clamped{free.template as<T>()};
        CHECK(gap(sampled<T>(), bounded, clamped, feed, control, saturated, y) <= tol<T>);
        CHECK(saturated > 10);
        CHECK(static_cast<double>(y) == doctest::Approx(1.0).epsilon(1e-3));

        const design::KalmanNoise<2, 1, double, 1> noise{.process = kNoise.process, .measurement = kNoise.measurement, .disturbance = Matrix<1, 1>{{1.0}}};
        const auto                                 of_bound = design::state_feedback(plant(), kTs, kFree, kBound, noise, design::InputDisturbance{});
        const auto                                 of_free = design::state_feedback(plant(), kTs, kFree, noise, design::InputDisturbance{});
        REQUIRE(of_bound.success);
        OFLQG of_bounded{of_bound.template as<T>()};
        OFLQG of_clamped{of_free.template as<T>()};
        CHECK(gap(sampled<T>(), of_bounded, of_clamped, feed, control, saturated, y) <= tol<T>);
        CHECK(saturated > 0);
        CHECK(static_cast<double>(y) == doctest::Approx(1.0).epsilon(1e-3));

        const auto lqg_bound = design::state_feedback(plant(), kTs, kFree, kBound, kNoise);
        REQUIRE(lqg_bound.success);
        LQG regulator{lqg_bound.template as<T>()};
        regulator.update(ColVec<1, T>{T{100}});
        const T u = regulator.control()(0);
        CHECK(damp::abs(u) == static_cast<T>(kLimit));
        CHECK(regulator.last_command()(0) == u); // what predict moves the filter on

        // An LQG assembled from its halves keeps the regulator's limits.
        const auto parts = design::lqg_from_parts(design::kalman(plant(), kTs, kNoise), design::state_feedback(plant(), kTs, kFree, kBound));
        REQUIRE(parts.success);
        CHECK(parts.limits.input.upper(0) == kLimit);
    }

    TEST_CASE_TEMPLATE("the integral action's limits clamp the integral state, and move at run time", T, float, double) {
        auto ia = kIntegral;
        ia.limits = {.lower = -0.02, .upper = 0.05};
        const auto res = design::state_feedback(plant(), kTs, kFree, ia);
        REQUIRE(res.success);
        LQI        law{res.template as<T>()};
        const auto hold = [&](T r) { // the plant held at rest: a constant error
            for (int k = 0; k < 2000; ++k) {
                law.update(ColVec<2, T>{});
                static_cast<void>(law.control(ColVec<1, T>{r}));
                law.predict();
            }
            return law.integral.xi(0);
        };
        CHECK(hold(T{1}) == static_cast<T>(0.05));
        CHECK(hold(T{-1}) == static_cast<T>(-0.02));
        auto lim = law.limits();
        lim.integral.lower(0) = static_cast<T>(-0.5);
        REQUIRE(law.set_limits(lim));
        CHECK(hold(T{-1}) == static_cast<T>(-0.5));
    }

    TEST_CASE_TEMPLATE("the disturbance limit clamps OFLQG's disturbance estimate, and moves at run time", T, float, double) {
        const design::KalmanNoise<2, 1, double, 1> noise{.process = kNoise.process, .measurement = kNoise.measurement, .disturbance = Matrix<1, 1>{{1.0}}};
        const auto                                 res = design::state_feedback(plant(), kTs, kFree, noise, design::InputDisturbance{});
        REQUIRE(res.success);
        OFLQG        law{res.template as<T>()};
        const auto   model = sampled<T>();
        const T      load = static_cast<T>(0.5); // on the input, the plant held at r = 0 against it
        ColVec<2, T> x{};
        const auto   hold = [&]() {
            T worst{0};
            for (int k = 0; k < 3000; ++k) {
                law.update(ColVec<1, T>((model.C * x)));
                const ColVec<1, T> u = law.control(ColVec<1, T>{});
                law.predict();
                x = ColVec<2, T>((model.A * x) + (model.B * ColVec<1, T>(u + ColVec<1, T>{load})));
                worst = damp::max(worst, damp::abs(law.disturbance_estimate()(0)));
            }
            return worst;
        };
        auto lim = law.limits();
        lim.disturbance = {.lower = static_cast<T>(-0.2), .upper = static_cast<T>(0.2)};
        REQUIRE(law.set_limits(lim));
        CHECK(hold() <= static_cast<T>(0.2));
        CHECK(law.disturbance_estimate()(0) == static_cast<T>(0.2));
        lim.disturbance.upper(0) = T{1};
        REQUIRE(law.set_limits(lim));
        static_cast<void>(hold());
        CHECK(static_cast<double>(law.disturbance_estimate()(0)) == doctest::Approx(0.5).epsilon(1e-3));
        CHECK(static_cast<double>((model.C * x)(0, 0)) == doctest::Approx(0.0).epsilon(1e-3)); // offset-free again
    }

    TEST_CASE_TEMPLATE("a DemandLimit is the outer terms' demand clamped outside with back_calculate", T, float, double) {
        // Read as a cascade on the speed: position and the integral ask for a speed, at most 0.3. The
        // oracle is the same law with the clamp outside it, as a hand-written cascade limit would be.
        // A mass on a damper: no input holds it anywhere, so the demand is 0 at rest (with a spring the
        // outer terms would carry its force).
        constexpr StateSpace mass{.A = Matrix<2, 2>{{0.0, 1.0}, {0.0, -0.5}}, .B = Matrix<2, 1>{{0.0}, {1.0}}, .C = Matrix<1, 2>{{1.0, 0.0}}};
        const auto           res = design::state_feedback(mass, kTs, kFree, kIntegral);
        REQUIRE(res.success);
        LQI     bounded{res.template as<T>()};
        LQI     clamped{res.template as<T>()};
        const T lim = static_cast<T>(0.3);
        auto    limits = bounded.limits();
        limits.demand = DemandLimit{.state = 1, .outer = {0, 2}, .lower = -lim, .upper = lim}; // the position and the integrator
        REQUIRE(bounded.set_limits(limits));
        const auto&        K = clamped.K; // [p, v, ξ]
        const auto         d = discretize(mass, kTs).value();
        const auto         model = StateSpace{.A = d.A.template as<T>(), .B = d.B.template as<T>(), .C = d.C.template as<T>(), .Ts = static_cast<T>(kTs)};
        const ColVec<1, T> r{T{5}};
        ColVec<2, T>       xa{};
        ColVec<2, T>       xb{};
        T                  worst{0};
        T                  worst_free_demand{0};
        int                limited = 0;
        for (int k = 0; k < 3000; ++k) {
            bounded.update(xa);
            clamped.update(xb);
            const T ua = bounded.control(r)(0);
            T       ub = clamped.control(r)(0);
            const T outer = -((K(0, 0) * (xb(0) - T{0})) + (K(0, 2) * clamped.integral.xi(0)));
            const T asked = outer / K(0, 1);
            worst_free_demand = damp::max(worst_free_demand, damp::abs(asked));
            if (damp::abs(asked) > lim) {
                ub += (K(0, 1) * damp::clamp(asked, -lim, lim)) - outer;
                ++limited;
            }
            clamped.back_calculate(ColVec<1, T>{ub});
            bounded.predict();
            clamped.predict();
            worst = damp::max(worst, damp::abs(ua - ub));
            xa = ColVec<2, T>((model.A * xa) + (model.B * ColVec<1, T>{ua}));
            xb = ColVec<2, T>((model.A * xb) + (model.B * ColVec<1, T>{ub}));
        }
        CHECK(limited > 100);
        CHECK(worst_free_demand > lim); // the limit bound: without it the outer terms ask for more
        CHECK(worst <= T{64} * tol<T> * T{5});
        CHECK(static_cast<double>(xa(0)) == doctest::Approx(5.0).epsilon(1e-3)); // the integrator held, then settled
        // OFLQG: the same on the deviation from its target, d̂ standing in for the integral.
        const design::KalmanNoise<2, 1, double, 1> noise{.process = kNoise.process, .measurement = kNoise.measurement, .disturbance = Matrix<1, 1>{{1.0}}};
        constexpr StateSpace                       noisy{.A = mass.A, .B = mass.B, .C = mass.C, .G = Matrix<2, 2>::identity(), .H = Matrix<1, 1>::identity()};
        const auto                                 of = design::state_feedback(noisy, kTs, kFree, noise, design::InputDisturbance{});
        REQUIRE(of.success);
        OFLQG of_bounded{of.template as<T>()};
        OFLQG of_clamped{of.template as<T>()};
        auto  of_limits = of_bounded.limits();
        of_limits.demand = DemandLimit{.state = 1, .outer = {0}, .lower = -lim, .upper = lim};
        REQUIRE(of_bounded.set_limits(of_limits));
        const auto& Kf = of_clamped.K;
        xa = ColVec<2, T>{};
        xb = ColVec<2, T>{};
        worst = T{0};
        limited = 0;
        for (int k = 0; k < 3000; ++k) {
            of_bounded.update(ColVec<1, T>(model.C * xa));
            of_clamped.update(ColVec<1, T>(model.C * xb));
            const T ua = of_bounded.control(r)(0);
            T       ub = of_clamped.control(r)(0);
            const T outer = -(Kf(0, 0) * of_clamped.estimator.at_command()(0));
            const T asked = outer / Kf(0, 1);
            if (damp::abs(asked) > lim) {
                ub += (Kf(0, 1) * damp::clamp(asked, -lim, lim)) - outer;
                ++limited;
            }
            of_clamped.back_calculate(ColVec<1, T>{ub});
            of_bounded.predict();
            of_clamped.predict();
            worst = damp::max(worst, damp::abs(ua - ub));
            xa = ColVec<2, T>((model.A * xa) + (model.B * ColVec<1, T>{ua}));
            xb = ColVec<2, T>((model.A * xb) + (model.B * ColVec<1, T>{ub}));
        }
        CHECK(limited > 50);
        CHECK(worst <= T{64} * tol<T> * T{5});
        CHECK(static_cast<double>(xa(0)) == doctest::Approx(5.0).epsilon(1e-3));

        LQI off{res.template as<T>()}; // the demand limit off again: the plain law
        limits.demand = {};
        REQUIRE(bounded.set_limits(limits));
        bounded.reset();
        off.update(ColVec<2, T>{T{-3}, T{0}});
        bounded.update(ColVec<2, T>{T{-3}, T{0}});
        CHECK(bounded.control(r) == off.control(r));
    }

    TEST_CASE("a DemandLimit among the options: a constexpr design's Law holds it in no storage, and a bad one is refused") {
        static constexpr StateSpace  kMass{.A = Matrix<2, 2>{{0.0, 1.0}, {0.0, -0.5}}, .B = Matrix<2, 1>{{0.0}, {1.0}}, .C = Matrix<1, 2>{{1.0, 0.0}}};
        static constexpr DemandLimit kSpeed{.state = 1, .outer = {0, 2}, .lower = -0.3, .upper = 0.3}; // the position and the integrator ask for the speed
        static constexpr auto        kBounded = design::state_feedback(kMass, kTs, kFree, kSpeed, design::IntegralAction{.Q = 80.0});
        static constexpr auto        kPlain = design::state_feedback(kMass, kTs, kFree, kIntegral);
        static_assert(kBounded.success);
        static_assert(sizeof(Law<kBounded, float>) == sizeof(Law<kPlain, float>)); // a constant of the type
        Law<kBounded, float> constant{};
        LQI                  held{kBounded.as<float>()};
        CHECK(held.limits().demand.state == 1);
        CHECK(constant.limits().demand.same(held.limits().demand)); // limits() reports the constant
        for (int k = 0; k < 50; ++k) {
            const ColVec<2, float> x{-4.0f + (0.01f * static_cast<float>(k)), 0.0f};
            CHECK(constant.step(ColVec<1, float>{1.0f}, x) == held.step(ColVec<1, float>{1.0f}, x));
        }

        CHECK(design::state_feedback(kMass, kTs, kFree, DemandLimit{.state = 2, .outer = {0}}, design::IntegralAction{.Q = 80.0}).status == design::DesignStatus::invalid_parameter); // not a plant state
        CHECK(design::state_feedback(kMass, kTs, kFree, DemandLimit{.state = 1, .outer = {1}}, design::IntegralAction{.Q = 80.0}).status == design::DesignStatus::invalid_parameter); // asks of itself
        CHECK(design::state_feedback(kMass, kTs, kFree, DemandLimit{.state = 1, .outer = {0}, .lower = 1.0, .upper = -1.0}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::state_feedback(kMass, kTs, kFree, DemandLimit{.state = 1, .outer = {2}}).status == design::DesignStatus::invalid_parameter);                                       // the integrator, on a law without one
        CHECK(design::state_feedback(kMass, kTs, kFree, DemandLimit{.state = 1, .outer = {0, 0}}, design::IntegralAction{.Q = 80.0}).status == design::DesignStatus::invalid_parameter); // listed twice
        CHECK(design::state_feedback(kMass, kTs, kFree, DemandLimit{.outer = {0}}).status == design::DesignStatus::invalid_parameter);                                                   // no state
        const auto lqr = design::state_feedback(kMass, kTs, kFree, DemandLimit{.state = 1, .outer = {0}, .lower = -0.3, .upper = 0.3});
        REQUIRE(lqr.success);
        CHECK(LQR{lqr}.limits().demand.state == 1); // the plain regulator's result carries it too
    }

    TEST_CASE("a constexpr design's Law carries the limits, in float") {
        static constexpr auto kLqi = design::state_feedback(plant(), kTs, {.Q = Q, .R = 1.0}, design::InputLimits{.lower = {-1.05}, .upper = {1.05}}, design::IntegralAction{.Q = 80.0, .limits = {.lower = -0.05, .upper = 0.05}});
        static constexpr auto kLqr = design::state_feedback(plant(), kTs, {.Q = Q, .R = 1.0}, design::InputLimits{.lower = {-1.05}, .upper = {1.05}});
        static_assert(kLqi.success && kLqr.success);
        Law<kLqi, float> lqi;
        Law<kLqr, float> lqr;
        CHECK(lqi.limits().input.upper(0) == 1.05f);
        CHECK(lqi.limits().integral.lower(0) == -0.05f);
        CHECK(lqr.limits().input.lower(0) == -1.05f);
        lqr.update(ColVec<2, float>{100.0f, 0.0f});
        CHECK(lqr.control()(0) == -1.05f);

        // No limit given: unbounded in float too (an end at double's max is float's, not an overflow).
        static constexpr auto kOpen = design::state_feedback(plant(), kTs, {.Q = Q, .R = 1.0}, design::IntegralAction{.Q = 80.0});
        Law<kOpen, float>     open;
        CHECK(open.limits().input.upper(0) == damp::numeric_limits<float>::max());
        CHECK(open.limits().integral.lower(0) == -damp::numeric_limits<float>::max());
    }

    TEST_CASE("refusals: a lower limit above its upper one") {
        const Bounds<1> inverted{.lower = 1.0, .upper = -1.0};
        CHECK(design::state_feedback(plant(), kTs, {.Q = Q, .R = 1.0}, design::InputLimits{.lower = inverted.lower, .upper = inverted.upper}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::state_feedback(plant(), kTs, {.Q = Q, .R = 1.0}, design::IntegralAction{.Q = 80.0, .limits = inverted}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::state_feedback(plant(), kTs, {.Q = Q, .R = 1.0}, design::InputLimits{.lower = inverted.lower, .upper = inverted.upper}, kNoise).status == design::DesignStatus::invalid_parameter);
    }

    TEST_CASE("InputLimits is the LQ design and law LQRLimits gave, through the config too") {
        //! Recorded with LQRLimits{.u_lim = kBox} before InputLimits replaced it: the gains, and the
        //! LQI's commands on a step that drives it onto the limit.
        const auto lqr = design::state_feedback(plant(), kTs, kFree, kBound);
        REQUIRE(lqr.success);
        CHECK(lqr.K(0, 0) == doctest::Approx(2.2848741073255541).epsilon(1e-12));
        CHECK(lqr.K(0, 1) == doctest::Approx(1.9137912936610608).epsilon(1e-12));
        CHECK(lqr.limits.input.lower(0) == -kLimit);
        CHECK(lqr.limits.input.upper(0) == kLimit);

        const auto lqi = design::state_feedback(plant(), kTs, kFree, kBound, kIntegral);
        REQUIRE(lqi.success);
        CHECK(lqi.K(0, 0) == doctest::Approx(8.2223504663959446).epsilon(1e-12));
        CHECK(lqi.K(0, 1) == doctest::Approx(3.6893318033860898).epsilon(1e-12));
        CHECK(lqi.K(0, 2) == doctest::Approx(-8.7792408363330239).epsilon(1e-12));
        //! The commands at ticks 0, 5, 50, 150 and 299, then the output at the end.
        const auto run = [](const auto& res) {
            const auto             d = discretize(plant(), kTs).value();
            LQI                    law{res};
            ColVec<2>              z{};
            damp::array<double, 6> out{};
            size_t                 next = 0;
            for (int k = 0; k < 300; ++k) {
                law.update(z);
                const ColVec<1> u = law.control(ColVec<1>{1.0});
                law.predict();
                if (k == 0 || k == 5 || k == 50 || k == 150 || k == 299) {
                    out[next++] = u(0);
                }
                z = ColVec<2>((d.A * z) + (d.B * u));
            }
            out[next] = z(0);
            return out;
        };
        const auto got = run(lqi);
        CHECK(got[0] == 0.0);
        CHECK(got[1] == doctest::Approx(0.40694881271785843).epsilon(1e-9));
        CHECK(got[2] == doctest::Approx(1.0433463780799417).epsilon(1e-9));
        CHECK(got[3] == doctest::Approx(0.55416557619751838).epsilon(1e-9));
        CHECK(got[4] == doctest::Approx(0.92109412242830047).epsilon(1e-9));
        CHECK(got[5] == doctest::Approx(1.0220145809046597).epsilon(1e-9));
    }

    TEST_CASE("an LQR with an InputRateLimits is the same LQR followed by a slew clamp and back_calculate") {
        const auto rated = design::state_feedback(plant(), kTs, kFree, kBound, kRate);
        const auto plain = design::state_feedback(plant(), kTs, kFree, kBound);
        REQUIRE(rated.success);
        static_assert(std::is_same_v<decltype(rated.limits.input_rate), Bounds<1>> && std::is_same_v<decltype(plain.limits.input_rate), Bounds<0>>);
        CHECK((rated.K - plain.K).norm() <= 4.0 * eps<double>() * plain.K.norm()); // the gain is the same
        CHECK(rated.limits.input_rate.upper(0) == kStep);
        const auto d = discretize(plant(), kTs).value();
        //! @p pass: the share of the command a stage after the law lets through, reported to the
        //! rated law with back_calculate when @p report (else its base is its last command).
        const auto run = [&](bool report, double pass) {
            LQR       law{rated};
            LQR       outside{plain};
            ColVec<2> za{3.0, -1.0};
            ColVec<2> zb = za;
            double    prev = 0.0; // the base of the first move, as an MPC's held input
            double    worst = 0.0;
            int       bound = 0;
            for (int k = 0; k < 4000; ++k) {
                law.update(za);
                outside.update(zb);
                const double ua = law.control()(0);
                const double asked = outside.control()(0);
                const double ub = damp::clamp(asked, prev - kStep, prev + kStep); // the slew clamp, by hand
                bound += (ub != asked) ? 1 : 0;
                worst = damp::max(worst, damp::abs(ua - ub));
                CHECK(damp::abs(ua - prev) <= kStep + 1e-12); // never more than the limit from the input applied last tick
                prev = pass * ub;                             // what went out
                if (report) {
                    law.back_calculate(ColVec<1>{prev});
                }
                outside.back_calculate(ColVec<1>{prev});
                za = (d.A * za) + (d.B * ColVec<1>{prev});
                zb = (d.A * zb) + (d.B * ColVec<1>{prev});
            }
            CHECK(bound > 10); // the far start drives the move onto its limit
            CHECK(worst <= 1e-12);
            return za;
        };
        CHECK(run(false, 1.0).norm() <= 1e-6); // the regulator settles
        CHECK(run(true, 0.9).norm() <= 1e-6);
    }

    TEST_CASE("an LQ law clamps its move into the InputLimits as well; where the boxes do not meet, the InputLimits win") {
        constexpr double step = 0.1;
        LQR              law{design::state_feedback(plant(), kTs, kFree, kBound, design::InputRateLimits{.lower = {-step}, .upper = {step}})};
        law.update(ColVec<2>{-100.0, 0.0});
        double prev = 0.0;
        for (int k = 0; k < 20; ++k) { // up by a step a tick until the input box stops it
            const double u = law.control()(0);
            CHECK(u == doctest::Approx(damp::min(prev + step, kLimit)).epsilon(1e-15));
            prev = u;
        }
        CHECK(prev == kLimit);
        auto lim = law.limits();
        lim.input = {.lower = -0.5, .upper = 0.5}; // u_prev = 1.05 now lies outside the input box
        REQUIRE(law.set_limits(lim));
        CHECK(law.control()(0) == 0.5); // [0.95, 1.15] and [−0.5, 0.5] do not meet: the input box's nearest end
        law.update(ColVec<2>{100.0, 0.0});
        CHECK(law.control()(0) == doctest::Approx(0.5 - step).epsilon(1e-15)); // then a move from there

        //! After reset() the next move is from 0, as an MPC's.
        law.reset();
        CHECK(law.control()(0) == doctest::Approx(-step).epsilon(1e-15));
    }

    TEST_CASE_TEMPLATE("an LQ law's memory holds against its rate clamp as against one outside it: LQI, delayed LQI, LQGI, OFLQG", T, float, double) {
        //! The rated law and the free one clamped outside with back_calculate are one tick: the
        //! integrators compare the input applied with the one asked, a filter and the delay states
        //! move on the clamped command.
        const auto measured = [](auto& law, const ColVec<2, T>& z) { law.update(z); };
        const auto read = [](auto& law, const auto& z) { law.update(ColVec<1, T>{z(0)}); };
        int        bound = 0;
        T          y{0};

        LQI lqi{design::state_feedback(plant(), kTs, kFree, kRate, kIntegral).template as<T>()};
        LQI lqi_out{design::state_feedback(plant(), kTs, kFree, kIntegral).template as<T>()};
        CHECK(slew_gap(sampled<T>(), lqi, lqi_out, measured, bound, y) <= tol<T>);
        CHECK(bound > 10);
        CHECK(static_cast<double>(y) == doctest::Approx(1.0).epsilon(1e-3)); // the integrator held: no windup left to unwind

        const InputDelay<2> late{0.014};
        const auto          dl = discretize(plant(), kTs, late).value();
        const StateSpace    delayed{.A = dl.A.template as<T>(), .B = dl.B.template as<T>(), .C = dl.C.template as<T>(), .Ts = static_cast<T>(kTs)};
        LQI                 lqi_late{design::state_feedback(plant(), kTs, kFree, late, kRate, kIntegral).template as<T>()};
        LQI                 lqi_late_out{design::state_feedback(plant(), kTs, kFree, late, kIntegral).template as<T>()};
        const auto          plant_state = [](auto& law, const ColVec<4, T>& z) { law.update(ColVec<2, T>{z(0), z(1)}); };
        CHECK(slew_gap(delayed, lqi_late, lqi_late_out, plant_state, bound, y) <= tol<T>);
        CHECK(bound > 10);
        CHECK(static_cast<double>(y) == doctest::Approx(1.0).epsilon(1e-3));

        LQGI lqgi{design::state_feedback(plant(), kTs, kFree, kRate, kNoise, kIntegral).template as<T>()};
        LQGI lqgi_out{design::state_feedback(plant(), kTs, kFree, kNoise, kIntegral).template as<T>()};
        CHECK(slew_gap(sampled<T>(), lqgi, lqgi_out, read, bound, y) <= tol<T>);
        CHECK(bound > 10);
        CHECK(static_cast<double>(y) == doctest::Approx(1.0).epsilon(1e-3));

        const design::KalmanNoise<2, 1, double, 1> walk{.process = kNoise.process, .measurement = kNoise.measurement, .disturbance = Matrix<1, 1>{{1.0}}};
        OFLQG                                      of{design::state_feedback(plant(), kTs, kFree, kRate, walk, design::InputDisturbance{}).template as<T>()};
        OFLQG                                      of_out{design::state_feedback(plant(), kTs, kFree, walk, design::InputDisturbance{}).template as<T>()};
        CHECK(slew_gap(sampled<T>(), of, of_out, read, bound, y) <= tol<T>);
        CHECK(bound > 10);
        CHECK(static_cast<double>(y) == doctest::Approx(1.0).epsilon(1e-3));
    }

    TEST_CASE("an LQ law's move limit is in its limits() only with the option, and set_limits moves it or refuses a bad one") {
        LQR rated{design::state_feedback(plant(), kTs, kFree, kRate)};
        LQR plain{design::state_feedback(plant(), kTs, kFree)};
        static_assert(std::is_same_v<decltype(rated.limits().input_rate), Bounds<1>> && std::is_same_v<decltype(plain.limits().input_rate), Bounds<0>>);
        auto lim = rated.limits();
        CHECK(lim.input_rate.lower(0) == -kStep);
        lim.input_rate = {.lower = -0.25, .upper = 0.5};
        REQUIRE(rated.set_limits(lim));
        rated.update(ColVec<2>{-100.0, 0.0});
        CHECK(rated.control()(0) == 0.5);
        rated.update(ColVec<2>{100.0, 0.0});
        CHECK(rated.control()(0) == 0.25);
        const double nan = std::numeric_limits<double>::quiet_NaN();
        for (const Bounds<1>& bad : {Bounds<1>{.lower = 0.1, .upper = -0.1}, Bounds<1>{.lower = nan}, Bounds<1>{.upper = nan}}) {
            auto refused = rated.limits();
            refused.input_rate = bad;
            CHECK_FALSE(rated.set_limits(refused));
            CHECK(rated.limits().input_rate.upper(0) == 0.5); // nothing changed
        }
    }

    TEST_CASE("an LQ design refuses an InputRateLimits that is not a box (invalid_parameter) and ships no limits") {
        const double                               nan = std::numeric_limits<double>::quiet_NaN();
        constexpr design::InputRateLimits          inverted{.lower = {0.1}, .upper = {-0.1}};
        const design::KalmanNoise<2, 1, double, 1> walk{.process = kNoise.process, .measurement = kNoise.measurement, .disturbance = Matrix<1, 1>{{1.0}}};
        const auto                                 refused = [](const auto& res) { return !res.success && res.status == design::DesignStatus::invalid_parameter; };
        CHECK(refused(design::state_feedback(plant(), kTs, kFree, kBound, inverted)));                           // LQR
        CHECK(refused(design::state_feedback(plant(), kTs, kFree, inverted, kIntegral)));                        // LQI
        CHECK(refused(design::state_feedback(plant(), kTs, kFree, inverted, kNoise)));                           // LQG
        CHECK(refused(design::state_feedback(plant(), kTs, kFree, inverted, kNoise, kIntegral)));                // LQGI
        CHECK(refused(design::state_feedback(plant(), kTs, kFree, inverted, walk, design::InputDisturbance{}))); // OFLQG
        CHECK(refused(design::state_feedback(plant(), kTs, kFree, inverted, InputDelay<2>{0.014})));             // on the delayed plant
        CHECK(refused(design::state_feedback(discretize(plant(), kTs).value(), kFree, inverted)));               // a discrete plant
        CHECK(refused(design::state_feedback(plant(), kTs, kFree, design::InputRateLimits{.upper = {nan}})));
        CHECK(refused(design::state_feedback(plant(), kTs, kFree, design::InputRateLimits{.lower = {nan}})));
        const auto none = design::state_feedback(plant(), kTs, kFree, kBound, inverted);
        CHECK(none.limits.input.upper(0) == damp::numeric_limits<double>::max()); // a refused design ships no limits
        CHECK(none.limits.input_rate.upper(0) == damp::numeric_limits<double>::max());

        //! Unbounded, it changes no command: the law of the InputLimits alone.
        LQR spelled{design::state_feedback(plant(), kTs, kFree, kBound, design::InputRateLimits<1>{})};
        LQR alone{design::state_feedback(plant(), kTs, kFree, kBound)};
        for (const ColVec<2>& x : {ColVec<2>{3.0, -1.0}, ColVec<2>{0.01, 0.0}}) {
            spelled.update(x);
            alone.update(x);
            CHECK(spelled.control() == alone.control());
        }
    }

    TEST_CASE("a constexpr design's Law with an InputRateLimits is the runtime's law, and set_limits moves its move limit") {
        static constexpr auto             kRes = design::state_feedback(plant(), kTs, kFree, kBound, kRate, kIntegral);
        static constinit Law<kRes, float> at_start{};
        Law<kRes, float>                  constant{};
        LQI                               held{kRes.as<float>()};
        CHECK(at_start.limits().input_rate.upper(0) == static_cast<float>(kStep));
        for (int k = 0; k < 50; ++k) {
            const ColVec<2, float> x{-4.0f + (0.01f * static_cast<float>(k)), 0.0f};
            CHECK(constant.step(ColVec<1, float>{1.0f}, x) == held.step(ColVec<1, float>{1.0f}, x));
        }
        auto lim = constant.limits();
        lim.input_rate = {.lower = -0.01f, .upper = 0.01f};
        REQUIRE(constant.set_limits(lim));
        REQUIRE(held.set_limits(lim));
        for (int k = 0; k < 50; ++k) {
            const ColVec<2, float> x{-4.0f, 0.0f};
            const ColVec<1, float> u = constant.step(ColVec<1, float>{1.0f}, x);
            CHECK(u == held.step(ColVec<1, float>{1.0f}, x));
        }
    }

    TEST_CASE("a refused design ships no limits") {
        //! The refused result kept the inverted box, and a runtime built from it at run time took
        //! it with K = 0: control() returned a limit, not 0.
        const double    nan = std::numeric_limits<double>::quiet_NaN();
        const Bounds<1> inverted{.lower = 1.0, .upper = -1.0};
        const auto      lqr = design::state_feedback(plant(), kTs, {.Q = Q, .R = 1.0}, design::InputLimits{.lower = inverted.lower, .upper = inverted.upper});
        REQUIRE_FALSE(lqr.success);
        CHECK(lqr.limits.input.lower(0) == -damp::numeric_limits<double>::max());
        CHECK(lqr.limits.input.upper(0) == damp::numeric_limits<double>::max());
        LQR law{lqr};
        law.update(ColVec<2>{1.0, 1.0});
        CHECK(law.control()(0) == 0.0);
        CHECK_FALSE(law.valid());
        CHECK(LQR{design::state_feedback(plant(), kTs, {.Q = Q, .R = 1.0})}.valid());

        const auto lqi = design::state_feedback(plant(), kTs, {.Q = Q, .R = 1.0}, design::InputLimits{.lower = {nan}, .upper = {1.0}}, design::IntegralAction{.Q = 80.0, .limits = {.lower = -0.1, .upper = 0.1}});
        CHECK(lqi.status == design::DesignStatus::invalid_parameter);
        CHECK(lqi.limits.input.upper(0) == damp::numeric_limits<double>::max());
        CHECK(lqi.limits.integral.upper(0) == damp::numeric_limits<double>::max());

        //! An LQG assembled from a regulator whose box was edited by hand is refused as the design is.
        auto reg = design::state_feedback(plant(), kTs, kFree, kBound);
        REQUIRE(reg.success);
        reg.limits.input = inverted;
        const auto parts = design::lqg_from_parts(design::kalman(plant(), kTs, kNoise), reg);
        CHECK(parts.status == design::DesignStatus::invalid_parameter);
        CHECK(parts.limits.input.upper(0) == damp::numeric_limits<double>::max());
    }

    TEST_CASE_TEMPLATE("an integrator limit narrowed while the integrator is held pulls it in", T, float, double) {
        //! The integrators' limit was applied only when the integrator advanced: held against the input limits, it kept
        //! a value outside a limit set since.
        const auto res = design::state_feedback(plant(), kTs, kFree, kBound, kIntegral);
        REQUIRE(res.success);
        LQI        law{res.template as<T>()};
        const auto tick = [&] {
            law.update(ColVec<2, T>{});
            const T u = law.control(ColVec<1, T>{T{100}})(0); // far above what the input limits allow
            law.predict();
            return u;
        };
        law.integral.xi(0) = T{0.5};
        CHECK(tick() == static_cast<T>(kLimit));
        CHECK(law.integral.xi(0) == T{0.5}); // held against the clamp
        auto lim = law.limits();
        lim.integral.upper(0) = T{0.25};
        REQUIRE(law.set_limits(lim));
        CHECK(tick() == static_cast<T>(kLimit));
        CHECK(law.integral.xi(0) == T{0.25});
    }

    TEST_CASE_TEMPLATE("tracking_time: back-calculation takes min(Ts/T_t, 1) of the excess back through (du/dxi)+", T, float, double) {
        //! The integrators could only hold against a clamp; a tracking time unwinds them, as Kbc does on a PID.
        for (const double tt : {4.0 * kTs, kTs, 0.5 * kTs}) { // a share of 1/4, the whole excess, and one capped at it
            auto ia = kIntegral;
            ia.tracking_time = tt;
            const auto res = design::state_feedback(plant(), kTs, kFree, kBound, ia);
            REQUIRE(res.success);
            LQI     law{res.template as<T>()};
            const T du_dxi = -law.getK()(0, 2);
            CHECK(law.tracking_time == static_cast<T>(tt));
            CHECK(static_cast<double>(law.integral.unwind(0, 0) * du_dxi) == doctest::Approx(1.0).epsilon(1e-5));
            law.integral.xi(0) = T{3} * static_cast<T>(kLimit) / du_dxi; // asks three times the limit
            law.update(ColVec<2, T>{});
            const T u = law.control(ColVec<1, T>{T{5}})(0);
            REQUIRE(u == static_cast<T>(kLimit));
            const T xi = law.integral.xi(0);
            const T e = law.integral.e(0);
            const T asked = law.integral.asked(0);
            law.predict();
            const T share = damp::min(static_cast<T>(kTs / tt), T{1});
            const T expect = (xi + (static_cast<T>(kTs) * e)) + (share * ((u - asked) / du_dxi));
            CHECK(damp::abs(law.integral.xi(0) - expect) <= T{16} * eps<T>() * (damp::abs(xi) + damp::abs(expect) + static_cast<T>(kTs) * damp::abs(e)));
        }
    }

    TEST_CASE_TEMPLATE("tracking_time: the law's own clamp and one reported with back_calculate are one tick", T, float, double) {
        auto ia = kIntegral;
        ia.tracking_time = 0.05;
        const auto bound = design::state_feedback(plant(), kTs, kFree, kBound, ia);
        const auto free = design::state_feedback(plant(), kTs, kFree, ia);
        REQUIRE(bound.success);
        REQUIRE(free.success);
        LQI  bounded{bound.template as<T>()};
        LQI  clamped{free.template as<T>()};
        int  saturated = 0;
        T    y{0};
        auto feed = [](auto& law, const ColVec<2, T>& z) { law.update(z); };
        auto control = [](auto& law) { return law.control(ColVec<1, T>{T{1}}); };
        CHECK(gap(sampled<T>(), bounded, clamped, feed, control, saturated, y) <= tol<T>);
        CHECK(saturated > 10);
        CHECK(static_cast<double>(y) == doctest::Approx(1.0).epsilon(1e-3));
    }

    TEST_CASE("tracking_time 0 holds: noise that grazes the clamp leaves the law at r, where back-calculation takes it from one side") {
        //! 0.3.0 unwound the integrators by default, the whole excess each tick. Near its limit with
        //! noise on the reading, the law took every excursion that reached the clamp back from its
        //! integrator and none from below, and settled short of r by many times holding's offset.
        const auto model = sampled<double>();
        const auto offset = [&model](double tracking_time) {
            auto ia = kIntegral;
            ia.tracking_time = tracking_time;
            const auto res = design::state_feedback(plant(), kTs, kFree, kBound, ia);
            CHECK(res.success);
            LQI           law{res.as<float>()};
            ColVec<2>     x{};
            unsigned      seed = 12345U;
            double        sum = 0.0;
            constexpr int kTicks = 6000;
            constexpr int kTail = 3000;
            for (int k = 0; k < kTicks; ++k) {
                seed = (seed * 1664525U) + 1013904223U;                                               // a fixed sequence: the same on every host
                const double noise = 0.0346 * ((static_cast<double>(seed >> 8U) / 16777216.0) - 0.5); // uniform, σ ≈ 0.01
                law.update(ColVec<2, float>{static_cast<float>(x(0) + noise), static_cast<float>(x(1))});
                const float u = law.control(ColVec<1, float>{1.0f})(0);
                law.predict();
                x = (model.A * x) + (model.B * ColVec<1>{static_cast<double>(u)});
                if (k >= kTicks - kTail) {
                    sum += x(0) - 1.0;
                }
            }
            return sum / kTail;
        };
        const double held = offset(0.0);
        const double tracked = offset(kTs);
        CHECK(damp::abs(held) < 0.015);                    // about 0.006
        CHECK(damp::abs(tracked) > 4.0 * damp::abs(held)); // about 0.05
    }

    TEST_CASE("tracking_time: refused when negative or not finite, carried by LQGI and by a constexpr design's Law") {
        CHECK(design::state_feedback(plant(), kTs, kFree, kBound, kIntegral).tracking_time == 0.0); // the default: the integrators hold, as Kbc = 0 on a PID
        auto ia = kIntegral;
        ia.tracking_time = -0.1;
        CHECK(design::state_feedback(plant(), kTs, kFree, kBound, ia).status == design::DesignStatus::invalid_parameter);
        ia.tracking_time = std::numeric_limits<double>::quiet_NaN();
        CHECK(design::state_feedback(plant(), kTs, kFree, kBound, ia).status == design::DesignStatus::invalid_parameter);

        //! An LQGI design dropped the integrators' limits when they moved out of the integral action: the filter's path did not copy it.
        auto g = kIntegral;
        g.limits = {.lower = -0.3, .upper = 0.3};
        g.tracking_time = 0.05;
        const auto lqgi = design::state_feedback(plant(), kTs, kFree, kBound, kNoise, g);
        REQUIRE(lqgi.success);
        CHECK(lqgi.limits.integral.upper(0) == 0.3);
        CHECK(lqgi.tracking_time == 0.05);
        const LQGI law{lqgi.as<float>()};
        CHECK(law.limits().integral.lower(0) == -0.3f);
        CHECK(law.tracking_time == 0.05f);
        CHECK(law.integral.unwind(0, 0) != 0.0f);

        static constexpr design::IntegralAction kTracked{.Q = 80.0, .tracking_time = 0.05};
        static constexpr auto                   kRes = design::state_feedback(plant(), kTs, kFree, kBound, kTracked);
        static_assert(kRes.success);
        const Law<kRes, float> constant;
        const LQI              held{kRes.as<float>()};
        CHECK(constant.tracking_time == 0.05f);
        CHECK(constant.integral.unwind(0, 0) == held.integral.unwind(0, 0));
    }

    TEST_CASE("set_limits refuses a box whose lower end is above its upper one, or a demand limit the law cannot hold, and changes nothing") {
        //! The limits were public members, and nothing checked one set at run time: an inverted input limit
        //! clamped every command to its upper end.
        const auto res = design::state_feedback(plant(), kTs, kFree, kBound, kIntegral);
        REQUIRE(res.success);
        LQI             law{res};
        const auto      before = law.limits();
        const Bounds<1> inverted{.lower = 1.0, .upper = -1.0};
        const auto      refused = [&](const auto& lim) {
            const bool taken = law.set_limits(lim);
            const auto now = law.limits();
            return !taken && now.input.lower(0) == before.input.lower(0) && now.input.upper(0) == before.input.upper(0)
                && now.integral.lower(0) == before.integral.lower(0) && now.integral.upper(0) == before.integral.upper(0)
                && now.demand.same(before.demand);
        };
        auto lim = before;
        lim.input = inverted;
        lim.integral = {.lower = -0.1, .upper = 0.1}; // a good limit beside it is not taken either
        CHECK(refused(lim));
        lim = before;
        lim.integral = inverted;
        CHECK(refused(lim));
        lim = before;
        lim.demand = DemandLimit{.state = 2, .outer = {0}}; // the integrator: not a plant state
        CHECK(refused(lim));
        lim.demand = DemandLimit{.state = 1, .outer = {1}}; // asks of itself
        CHECK(refused(lim));
        lim.demand = DemandLimit{.state = 1, .outer = {3}}; // past the terms
        CHECK(refused(lim));
        lim.demand = DemandLimit{.state = 1, .outer = {0, 0}}; // listed twice
        CHECK(refused(lim));
        lim.demand = DemandLimit{.state = 1, .outer = {0}, .lower = inverted.lower(0), .upper = inverted.upper(0)};
        CHECK(refused(lim));
        law.update(ColVec<2>{100.0, 0.0});
        CHECK(law.control(ColVec<1>{0.0})(0) == -kLimit); // the design's box, still

        //! Every LQ law refuses one, and OFLQG an inverted disturbance limit too.
        const auto inverted_refused = [&](auto runtime) {
            auto       l = runtime.limits();
            const auto u = l.input;
            l.input = inverted;
            return !runtime.set_limits(l) && runtime.limits().input.lower(0) == u.lower(0) && runtime.limits().input.upper(0) == u.upper(0);
        };
        const design::KalmanNoise<2, 1, double, 1> walk{.process = kNoise.process, .measurement = kNoise.measurement, .disturbance = Matrix<1, 1>{{1.0}}};
        CHECK(inverted_refused(LQR{design::state_feedback(plant(), kTs, kFree, kBound)}));
        CHECK(inverted_refused(LQG{design::state_feedback(plant(), kTs, kFree, kBound, kNoise)}));
        CHECK(inverted_refused(LQGI{design::state_feedback(plant(), kTs, kFree, kBound, kNoise, kIntegral)}));
        OFLQG of{design::state_feedback(plant(), kTs, kFree, kBound, walk, design::InputDisturbance{})};
        CHECK(inverted_refused(of));
        auto d = of.limits();
        d.disturbance = inverted;
        CHECK_FALSE(of.set_limits(d));
        CHECK(of.limits().disturbance.upper(0) == damp::numeric_limits<double>::max());
        d.disturbance = {.lower = -0.2, .upper = 0.2};
        CHECK(of.set_limits(d));
        CHECK(of.limits().disturbance.upper(0) == 0.2);
    }

    TEST_CASE("a constexpr design's Law sets its limits as the runtime does, its demand limit a constant it reports and keeps") {
        static constexpr StateSpace  kMass{.A = Matrix<2, 2>{{0.0, 1.0}, {0.0, -0.5}}, .B = Matrix<2, 1>{{0.0}, {1.0}}, .C = Matrix<1, 2>{{1.0, 0.0}}};
        static constexpr DemandLimit kSpeed{.state = 1, .outer = {0, 2}, .lower = -0.3, .upper = 0.3};
        static constexpr auto        kRes = design::state_feedback(kMass, kTs, kFree, kBound, kSpeed, design::IntegralAction{.Q = 80.0});
        static_assert(kRes.success);
        static constinit Law<kRes, float> at_start{}; // set_limits in the constructor, constant-evaluated
        Law<kRes, float>                  constant{};
        LQI                               held{kRes.as<float>()};
        CHECK(at_start.limits().input.upper(0) == static_cast<float>(kLimit));
        CHECK(at_start.limits().integral.upper(0) == damp::numeric_limits<float>::max());
        auto lim = constant.limits();
        CHECK(lim.demand.same(held.limits().demand));
        CHECK(lim.demand.state == 1);
        lim.input = {.lower = -0.5f, .upper = 0.5f};
        lim.integral = {.lower = -0.2f, .upper = 0.2f};
        REQUIRE(constant.set_limits(lim));
        REQUIRE(held.set_limits(lim));
        CHECK(constant.limits().input.upper(0) == 0.5f);
        CHECK(constant.limits().integral.lower(0) == -0.2f);
        for (int k = 0; k < 50; ++k) {
            const ColVec<2, float> x{-4.0f + (0.01f * static_cast<float>(k)), 0.0f};
            const ColVec<1, float> u = constant.step(ColVec<1, float>{1.0f}, x);
            CHECK(u == held.step(ColVec<1, float>{1.0f}, x));
            CHECK(damp::abs(u(0)) <= 0.5f);
        }

        //! Another demand limit, or none, is refused on constants (nothing changed); the runtime with members takes it.
        auto other = lim;
        other.input = {.lower = -0.1f, .upper = 0.1f};
        other.demand = DemandLimit{.state = 1, .outer = {0}, .lower = -0.3f, .upper = 0.3f};
        CHECK_FALSE(constant.set_limits(other));
        other.demand = {};
        CHECK_FALSE(constant.set_limits(other));
        CHECK(constant.limits().input.upper(0) == 0.5f);
        CHECK(held.set_limits(other));
        CHECK(held.limits().demand.state == detail::DemandTerms<3, float>::none);
    }

    TEST_CASE("the limits take no storage of their own: each LQ law is the size it was with public limits") {
        //! Recorded on x86-64 with the limits public members: limits() and set_limits() read and
        //! write the same members, and the per-tick work is the same.
        static constexpr design::KalmanNoise<2, 1, double, 1> kWalk{.process = kNoise.process, .measurement = kNoise.measurement, .disturbance = Matrix<1, 1>{{1.0}}};
        static constexpr auto                                 kLqr = design::state_feedback(plant(), kTs, kFree);
        static constexpr auto                                 kLqi = design::state_feedback(plant(), kTs, kFree, kIntegral);
        static constexpr auto                                 kLqg = design::state_feedback(plant(), kTs, kFree, kNoise);
        static constexpr auto                                 kLqgi = design::state_feedback(plant(), kTs, kFree, kNoise, kIntegral);
        static constexpr auto                                 kOflqg = design::state_feedback(plant(), kTs, kFree, kWalk, design::InputDisturbance{});
        if constexpr (sizeof(void*) == 8) {
            CHECK(sizeof(LQR{kLqr.as<float>()}) == 48);
            CHECK(sizeof(LQI{kLqi.as<float>()}) == 112);
            CHECK(sizeof(LQG{kLqg.as<float>()}) == 128);
            CHECK(sizeof(LQGI{kLqgi.as<float>()}) == 176);
            CHECK(sizeof(OFLQG{kOflqg.as<float>()}) == 224);
            CHECK(sizeof(LQR{kLqr}) == 88);
            CHECK(sizeof(LQI{kLqi}) == 216);
            CHECK(sizeof(LQG{kLqg}) == 240);
            CHECK(sizeof(LQGI{kLqgi}) == 336);
            CHECK(sizeof(OFLQG{kOflqg}) == 424);
            CHECK(sizeof(Law<kLqr, float>) == 20);
            CHECK(sizeof(Law<kLqi, float>) == 68);
            CHECK(sizeof(Law<kLqg, float>) == 44);
            CHECK(sizeof(Law<kLqgi, float>) == 80);
            CHECK(sizeof(Law<kOflqg, float>) == 80);
        }
    }

    TEST_CASE("an InputRateLimits costs an LQ law its box, and the last command where it keeps no record; a law without one is the size it was") {
        //! The sizes without the option are the test above's, recorded before InputRateLimits. With
        //! it a law holds the box (NU ends a side), and an LQR, which keeps no record of its
        //! command, the base of its next move too; the others read it from the record they keep
        //! for their integrators, filter or delay states.
        static constexpr design::KalmanNoise<2, 1, double, 1> kWalk{.process = kNoise.process, .measurement = kNoise.measurement, .disturbance = Matrix<1, 1>{{1.0}}};
        static constexpr auto                                 kLqr = design::state_feedback(plant(), kTs, kFree, kRate);
        static constexpr auto                                 kLqi = design::state_feedback(plant(), kTs, kFree, kIntegral, kRate);
        static constexpr auto                                 kLqg = design::state_feedback(plant(), kTs, kFree, kNoise, kRate);
        static constexpr auto                                 kLqgi = design::state_feedback(plant(), kTs, kFree, kNoise, kIntegral, kRate);
        static constexpr auto                                 kOflqg = design::state_feedback(plant(), kTs, kFree, kWalk, design::InputDisturbance{}, kRate);
        static_assert(std::is_same_v<decltype(LQR{design::state_feedback(plant(), kTs, kFree)}), LQR<2, 1, double>>); // without: the type it was
        if constexpr (sizeof(void*) == 8) {
            CHECK(sizeof(LQR{kLqr.as<float>()}) == 48 + 16);
            CHECK(sizeof(LQI{kLqi.as<float>()}) == 112 + 8);
            CHECK(sizeof(LQG{kLqg.as<float>()}) == 128 + 8);
            CHECK(sizeof(LQGI{kLqgi.as<float>()}) == 176 + 8);
            CHECK(sizeof(OFLQG{kOflqg.as<float>()}) == 224 + 8);
            CHECK(sizeof(LQR{kLqr}) == 88 + 24);
            CHECK(sizeof(LQI{kLqi}) == 216 + 16);
            CHECK(sizeof(LQG{kLqg}) == 240 + 16);
            CHECK(sizeof(LQGI{kLqgi}) == 336 + 16);
            CHECK(sizeof(OFLQG{kOflqg}) == 424 + 16);
            CHECK(sizeof(Law<kLqr, float>) == 20 + 12);
            CHECK(sizeof(Law<kLqi, float>) == 68 + 8);
            CHECK(sizeof(Law<kLqg, float>) == 44 + 8);
            CHECK(sizeof(Law<kLqgi, float>) == 80 + 8);
            CHECK(sizeof(Law<kOflqg, float>) == 80 + 8);
        }
    }

    TEST_CASE("an LQ law's limits and an MPC's are read and set the same way, and both bind from the next tick") {
        constexpr StateSpace     double_integrator{.A = {{1.0, 0.1}, {0.0, 1.0}}, .B = {{0.005}, {0.1}}, .C = {{1.0, 0.0}}, .Ts = 0.1};
        design::MPCWeights<1, 1> w{};
        w.Rdu = 0.05;
        const auto lq = design::state_feedback(double_integrator, design::LQRWeights{.Q = Matrix<2, 2>::identity(), .R = {0.01}});
        const auto mpc = design::mpc<20, 8>(double_integrator, w); // a regulator, as the LQR is: control()
        REQUIRE(lq.success);
        REQUIRE(mpc.success);
        LQR             lqr{lq};
        MPC             law{mpc};
        const ColVec<2> far{-1.0, 0.0};
        lqr.update(far);
        law.update(far);
        REQUIRE(lqr.control()(0) > 0.1);
        REQUIRE(law.control()(0) > 0.1);

        auto lq_lim = lqr.limits();
        auto mpc_lim = law.limits();
        lq_lim.input = {.lower = -0.1, .upper = 0.1};
        mpc_lim.input = {.lower = -0.1, .upper = 0.1};
        REQUIRE(lqr.set_limits(lq_lim));
        REQUIRE(law.set_limits(mpc_lim));
        CHECK(lqr.control()(0) == 0.1);
        CHECK(law.control()(0) == doctest::Approx(0.1).epsilon(1e-12));
        CHECK(law.last_status() == design::QPStatus::Success);

        lq_lim.input = {.lower = 0.1, .upper = -0.1};
        mpc_lim.input = {.lower = 0.1, .upper = -0.1};
        CHECK_FALSE(lqr.set_limits(lq_lim));
        CHECK_FALSE(law.set_limits(mpc_lim));
        CHECK(lqr.limits().input.upper(0) == 0.1);
        CHECK(law.limits().input.upper(0) == 0.1);
    }
}
