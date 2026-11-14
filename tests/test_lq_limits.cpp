// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// The limits of the LQ laws: u_lim clamps the command control returns, integral_lim the integral
// state. The oracle for u_lim is the same design without limits, clamped outside the law and told
// with back_calculate: the two are one tick, on every runtime, so the integrators hold against the
// law's own clamp and a filter or the input registers move on the input that went out.

#include <cmath>
#include <cstddef>

#include "damp/controllers/law_detail.hpp"
#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/oflqg.hpp"
#include "damp/design/disturbance.hpp"
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

constexpr Matrix<2, 2>                Q{{10.0, 0.0}, {0.0, 1.0}};
constexpr Bounds<1>                   kBox{.lower = -kLimit, .upper = kLimit};
constexpr design::KalmanNoise<2, 1>   kNoise{.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.1}}, .measurement = Matrix<1, 1>{{1e-4}}};
constexpr design::LQRWeights<2, 1, 0> kFree{.Q = Q, .R = 1.0};
constexpr design::LQRWeights<2, 1, 0> kBound{.Q = Q, .R = 1.0, .u_lim = kBox};
constexpr design::LQRWeights<2, 1, 1> kFreeI{.Q = Q, .R = 1.0, .Q_integral = 80.0};
constexpr design::LQRWeights<2, 1, 1> kBoundI{.Q = Q, .R = 1.0, .Q_integral = 80.0, .u_lim = kBox};

/**
 * Run @p bounded (its design carries u_lim) and @p clamped (the same design without, clamped
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
    TEST_CASE_TEMPLATE("LQR: control returns the command clamped into u_lim, and no limit is no clamp", T, float, double) {
        const auto free = design::lqr(plant(), kTs, kFree);
        const auto bound = design::lqr(plant(), kTs, kBound);
        REQUIRE(free.success);
        REQUIRE(bound.success);
        CHECK((bound.K - free.K).norm() <= 4.0 * eps<double>() * free.K.norm()); // a limit is the runtime's: the gain is the same
        LQR                law{bound.template as<T>()};
        LQR                open{free.template as<T>()};
        const ColVec<2, T> far{{T{100}, T{-40}}};
        const ColVec<2, T> near{{static_cast<T>(0.01), T{0}}};
        law.update(far);
        open.update(far);
        const T u_far = open.control()(0);
        REQUIRE(damp::abs(u_far) > static_cast<T>(kLimit));
        CHECK(law.control()(0) == ((u_far > T{0}) ? static_cast<T>(kLimit) : -static_cast<T>(kLimit)));
        law.update(near);
        open.update(near);
        const T u_near = open.control()(0);
        CHECK(damp::abs(law.control()(0) - u_near) <= T{4} * eps<T>() * damp::abs(u_near)); // inside the box: the law itself

        // The limits are members: a rail that moves is set at run time.
        law.u_lim.upper(0) = T{5};
        law.u_lim.lower(0) = T{-5};
        law.update(far);
        CHECK(damp::abs(law.control()(0)) == T{5});
    }

    TEST_CASE_TEMPLATE("LQI with u_lim is the law clamped outside with back_calculate: the integrator holds against its own clamp", T, float, double) {
        const auto bound = design::lqi(plant(), kTs, kBoundI);
        const auto free = design::lqi(plant(), kTs, kFreeI);
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

    TEST_CASE_TEMPLATE("LQI with an InputDelay: the registers hold the clamped command", T, float, double) {
        const InputDelay<2> late{0.014};
        const auto          bound = design::lqi(plant(), kTs, late, kBoundI);
        const auto          free = design::lqi(plant(), kTs, late, kFreeI);
        REQUIRE(bound.success);
        const auto                   d = discretize(plant(), kTs, late).value();
        const StateSpace<4, 1, 1, T> model{.A = d.A.template as<T>(), .B = d.B.template as<T>(), .C = d.C.template as<T>(), .Ts = static_cast<T>(kTs)};
        LQI                          bounded{bound.template as<T>()};
        LQI                          clamped{free.template as<T>()};
        int                          saturated = 0;
        T                            y{0};
        auto                         feed = [](auto& law, const ColVec<4, T>& z) { law.update(ColVec<2, T>{{z(0), z(1)}}); };
        auto                         control = [](auto& law) { return law.control(ColVec<1, T>{T{1}}); };
        CHECK(gap(model, bounded, clamped, feed, control, saturated, y) <= tol<T>);
        CHECK(saturated > 10);
        CHECK(static_cast<double>(y) == doctest::Approx(1.0).epsilon(1e-3));
    }

    TEST_CASE_TEMPLATE("LQGI and OFLQG with u_lim: the filter predicts with the clamped command", T, float, double) {
        auto feed = [](auto& law, const ColVec<2, T>& z) { law.update(ColVec<1, T>{z(0)}); };
        auto control = [](auto& law) { return law.control(ColVec<1, T>{T{1}}); };
        int  saturated = 0;
        T    y{0};

        const auto bound = design::lqgi(plant(), kTs, kBoundI, kNoise);
        const auto free = design::lqgi(plant(), kTs, kFreeI, kNoise);
        REQUIRE(bound.success);
        LQGI bounded{bound.template as<T>()};
        LQGI clamped{free.template as<T>()};
        CHECK(gap(sampled<T>(), bounded, clamped, feed, control, saturated, y) <= tol<T>);
        CHECK(saturated > 10);
        CHECK(static_cast<double>(y) == doctest::Approx(1.0).epsilon(1e-3));

        const design::KalmanNoise<2, 1, double, 1> noise{.process = kNoise.process, .measurement = kNoise.measurement, .disturbance = Matrix<1, 1>{{1.0}}};
        const auto                                 of_bound = design::oflqg(plant(), kTs, kBound, noise, design::InputDisturbance{});
        const auto                                 of_free = design::oflqg(plant(), kTs, kFree, noise, design::InputDisturbance{});
        REQUIRE(of_bound.success);
        OFLQG of_bounded{of_bound.template as<T>()};
        OFLQG of_clamped{of_free.template as<T>()};
        CHECK(gap(sampled<T>(), of_bounded, of_clamped, feed, control, saturated, y) <= tol<T>);
        CHECK(saturated > 0);
        CHECK(static_cast<double>(y) == doctest::Approx(1.0).epsilon(1e-3));

        const auto lqg_bound = design::lqg(plant(), kTs, kBound, kNoise);
        REQUIRE(lqg_bound.success);
        LQG regulator{lqg_bound.template as<T>()};
        regulator.update(ColVec<1, T>{T{100}});
        const T u = regulator.control()(0);
        CHECK(damp::abs(u) == static_cast<T>(kLimit));
        CHECK(regulator.last_command()(0) == u); // what predict moves the filter on

        // An LQG assembled from its halves keeps the regulator's limits.
        const auto parts = design::lqg_from_parts(design::kalman(plant(), kTs, kNoise), design::lqr(plant(), kTs, kBound));
        REQUIRE(parts.success);
        CHECK(parts.u_lim.upper(0) == kLimit);
    }

    TEST_CASE_TEMPLATE("integral_lim clamps the integral state, and moves at run time", T, float, double) {
        design::LQRWeights<2, 1, 1> w = kFreeI;
        w.integral_lim = {.lower = -0.02, .upper = 0.05};
        const auto res = design::lqi(plant(), kTs, w);
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
        law.integral.lim.lower(0) = static_cast<T>(-0.5);
        CHECK(hold(T{-1}) == static_cast<T>(-0.5));
    }

    TEST_CASE("a constexpr design's Law carries the limits, in float") {
        static constexpr auto kLqi = design::lqi(plant(), kTs, {.Q = Q, .R = 1.0, .Q_integral = 80.0, .u_lim = {.lower = -1.05, .upper = 1.05}, .integral_lim = {.lower = -0.05, .upper = 0.05}});
        static constexpr auto kLqr = design::lqr(plant(), kTs, {.Q = Q, .R = 1.0, .u_lim = {.lower = -1.05, .upper = 1.05}});
        static_assert(kLqi.success && kLqr.success);
        Law<kLqi, float> lqi;
        Law<kLqr, float> lqr;
        CHECK(lqi.u_lim.upper(0) == 1.05f);
        CHECK(lqi.integral.lim.lower(0) == -0.05f);
        CHECK(lqr.u_lim.lower(0) == -1.05f);
        lqr.update(ColVec<2, float>{{100.0f, 0.0f}});
        CHECK(lqr.control()(0) == -1.05f);

        // No limit given: unbounded in float too (an end at double's max is float's, not an overflow).
        static constexpr auto kOpen = design::lqi(plant(), kTs, {.Q = Q, .R = 1.0, .Q_integral = 80.0});
        Law<kOpen, float>     open;
        CHECK(open.u_lim.upper(0) == damp::numeric_limits<float>::max());
        CHECK(open.integral.lim.lower(0) == -damp::numeric_limits<float>::max());
    }

    TEST_CASE("refusals: a lower limit above its upper one") {
        const Bounds<1> inverted{.lower = 1.0, .upper = -1.0};
        CHECK(design::lqr(plant(), kTs, {.Q = Q, .R = 1.0, .u_lim = inverted}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::lqi(plant(), kTs, {.Q = Q, .R = 1.0, .Q_integral = 80.0, .integral_lim = inverted}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::lqg(plant(), kTs, {.Q = Q, .R = 1.0, .u_lim = inverted}, kNoise).status == design::DesignStatus::invalid_parameter);
    }
}
