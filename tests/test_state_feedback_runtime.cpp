// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// StateFeedback, the one runtime of the LQ designs: it deduces from any of them, LQR, LQI, LQG, LQGI
// and OFLQG are the same type for their own kind and refuse another's, Law<kRes, T> names it, and
// each law is the type and the size it was as one of five structs.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <type_traits>
#include <utility>

#include "damp/backend.hpp"
#include "damp/controllers/detail/state_feedback_core.hpp"
#include "damp/controllers/law.hpp"
#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/oflqg.hpp"
#include "damp/controllers/state_feedback.hpp"
#include "damp/controllers/state_feedback_runtime.hpp"
#include "damp/controllers/state_feedback_types.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/input_limits.hpp"
#include "damp/design/riccati.hpp"
#include "damp/estimation/fixed_gain_observer.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "doctest.h"

using namespace damp;

namespace {

constexpr double kTs = 0.01;

/// A mass on a spring and a damper: position, speed; force in; position out. Noise on both states and the reading.
constexpr StateSpace<2, 1, 1, double, 2, 1> plant() {
    return {
        .A = Matrix{{0.0, 1.0}, {-1.0, -0.5}},
        .B = Matrix{{0.0}, {1.0}},
        .C = Matrix{{1.0, 0.0}},
        .G = Matrix<2, 2>::identity(),
        .H = Matrix<1, 1>::identity(),
    };
}

constexpr design::LQRWeights<2, 1>             kW{.Q = Matrix{{10.0, 0.0}, {0.0, 1.0}}, .R = 1.0};
constexpr design::KalmanNoise<2, 1>            kNoise{.process = Matrix{{0.01, 0.0}, {0.0, 0.1}}, .measurement = Matrix{{1e-4}}};
constexpr design::KalmanNoise<2, 1, double, 1> kWalk{.process = kNoise.process, .measurement = kNoise.measurement, .disturbance = Matrix{{1.0}}};
constexpr design::IntegralAction               kInt{.Q = 80.0};
constexpr design::InputRateLimits              kRate{.lower = {-0.002}, .upper = {0.002}};
constexpr design::Handover                     kFade{.fade_time = 0.1};
constexpr Matrix                               kLoad{{0.0}, {1.0}};

// One design of each kind, with each option that changes the law's type: tracked outputs, a
// delayed input, a late reading, a measured disturbance, a move limit, a handover's fade, a gain
// or a filter placed at poles.
inline constexpr auto kLqr = design::state_feedback(plant(), kTs, kW);
inline constexpr auto kLqrTrack = design::state_feedback(plant(), kTs, kW, design::TrackedOutputs{{0}});
inline constexpr auto kLqrDelayed = design::state_feedback(plant(), kTs, kW, InputDelay{kTs});
inline constexpr auto kLqrMeasured = design::state_feedback(plant(), kTs, kW, design::MeasuredDisturbance{kLoad});
inline constexpr auto kLqrTrackMeasured = design::state_feedback(plant(), kTs, kW, design::TrackedOutputs{{0}}, design::MeasuredDisturbance{kLoad});
inline constexpr auto kLqrRate = design::state_feedback(plant(), kTs, kW, kRate);
inline constexpr auto kLqrFaded = design::state_feedback(plant(), kTs, kW, kFade);
inline constexpr auto kLqrPlaced = design::state_feedback(plant(), kTs, design::Poles{-5.0, -6.0});
inline constexpr auto kLqi = design::state_feedback(plant(), kTs, kW, kInt);
inline constexpr auto kLqiDelayed = design::state_feedback(plant(), kTs, kW, kInt, InputDelay{kTs});
inline constexpr auto kLqiMeasured = design::state_feedback(plant(), kTs, kW, kInt, design::MeasuredDisturbance{kLoad});
inline constexpr auto kLqiRate = design::state_feedback(plant(), kTs, kW, kInt, kRate);
inline constexpr auto kLqiFaded = design::state_feedback(plant(), kTs, kW, kInt, kFade);
inline constexpr auto kLqg = design::state_feedback(plant(), kTs, kW, kNoise);
inline constexpr auto kLqgTrack = design::state_feedback(plant(), kTs, kW, kNoise, design::TrackedOutputs{{0}});
inline constexpr auto kLqgDelayed = design::state_feedback(plant(), kTs, kW, kNoise, InputDelay{kTs});
inline constexpr auto kLqgLate = design::state_feedback(plant(), kTs, kW, kNoise, OutputDelay{kTs});
inline constexpr auto kLqgMeasured = design::state_feedback(plant(), kTs, kW, kNoise, design::MeasuredDisturbance{kLoad});
inline constexpr auto kLqgRate = design::state_feedback(plant(), kTs, kW, kNoise, kRate);
inline constexpr auto kLqgFaded = design::state_feedback(plant(), kTs, kW, kNoise, kFade);
inline constexpr auto kLqgPlaced = design::state_feedback(plant(), kTs, design::Poles{-5.0, -6.0}, design::ObserverPoles{-20.0, -25.0});
inline constexpr auto kLqgi = design::state_feedback(plant(), kTs, kW, kNoise, kInt);
inline constexpr auto kLqgiDelayed = design::state_feedback(plant(), kTs, kW, kNoise, kInt, InputDelay{kTs});
inline constexpr auto kLqgiMeasured = design::state_feedback(plant(), kTs, kW, kNoise, kInt, design::MeasuredDisturbance{kLoad});
inline constexpr auto kLqgiRate = design::state_feedback(plant(), kTs, kW, kNoise, kInt, kRate);
inline constexpr auto kLqgiFaded = design::state_feedback(plant(), kTs, kW, kNoise, kInt, kFade);
inline constexpr auto kOflqg = design::state_feedback(plant(), kTs, kW, kWalk, design::InputDisturbance{});
inline constexpr auto kOflqgDelayed = design::state_feedback(plant(), kTs, kW, kWalk, design::InputDisturbance{}, InputDelay{kTs});
inline constexpr auto kOflqgLate = design::state_feedback(plant(), kTs, kW, kWalk, design::InputDisturbance{}, OutputDelay{kTs});
inline constexpr auto kOflqgMeasured = design::state_feedback(plant(), kTs, kW, kWalk, design::InputDisturbance{}, design::MeasuredDisturbance{kLoad});
inline constexpr auto kOflqgRate = design::state_feedback(plant(), kTs, kW, kWalk, design::InputDisturbance{}, kRate);
inline constexpr auto kOflqgFaded = design::state_feedback(plant(), kTs, kW, kWalk, design::InputDisturbance{}, kFade);

template<const auto& kRes>
using Res = damp::remove_cvref_t<decltype(kRes)>;

using PlantFilter = FixedGainObserver<2, 1, 1, double>;     // on [x]
using DelayFilter = FixedGainObserver<3, 1, 1, double>;     // on [x; v(k−1)], or [x; z(k−1)] for a late reading
using MeasuredFilter = FixedGainObserver<2, 2, 1, double>;  // on [x], d a known input after u
using LoadFilter = FixedGainObserver<3, 1, 1, double, 1>;   // on [x; d]
using LoadDelay = FixedGainObserver<4, 1, 1, double, 1>;    // on [x; v(k−1); d], or [x; z(k−1); d]
using LoadMeasured = FixedGainObserver<3, 2, 1, double, 1>; // on [x; d], the measured d after u
using Held = detail::HeldGains;

// Whether a name deduces a law from a design: the alias refuses another kind's at deduction.
template<typename R>
concept feedback_takes = requires { StateFeedback{std::declval<const R&>()}; };
template<typename R>
concept lqr_takes = requires { LQR{std::declval<const R&>()}; };
template<typename R>
concept lqi_takes = requires { LQI{std::declval<const R&>()}; };
template<typename R>
concept lqg_takes = requires { LQG{std::declval<const R&>()}; };
template<typename R>
concept lqgi_takes = requires { LQGI{std::declval<const R&>()}; };
template<typename R>
concept oflqg_takes = requires { OFLQG{std::declval<const R&>()}; };

constexpr int kLQR = 0;
constexpr int kLQI = 1;
constexpr int kLQG = 2;
constexpr int kLQGI = 3;
constexpr int kOFLQG = 4;

/// Exactly one of the five names takes @p R, the one @p kind names, and StateFeedback takes it.
template<typename R>
constexpr bool only(int kind) {
    const bool takes[] = {lqr_takes<R>, lqi_takes<R>, lqg_takes<R>, lqgi_takes<R>, oflqg_takes<R>};
    bool       ok = feedback_takes<R>;
    for (int k = 0; k < 5; ++k) {
        ok = ok && (takes[k] == (k == kind));
    }
    return ok;
}

// An alias named with another kind's filter, or LQI and LQGI with no integrators, does not form.
template<typename Est>
concept forms_lqg = requires { typename LQG<2, 1, 1, double, Est>; };
template<typename Est>
concept forms_lqgi = requires { typename LQGI<2, 1, 1, 1, double, Est>; };
template<typename Est>
concept forms_oflqg = requires { typename OFLQG<2, 1, 1, 1, double, Est>; };
template<size_t NI>
concept forms_lqi = requires { typename LQI<2, 1, NI, double, 1>; };
template<size_t NI>
concept forms_lqgi_with = requires { typename LQGI<2, 1, 1, NI, double, PlantFilter>; };

/// The measured-disturbance inputs a StateFeedback takes (NM).
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, typename Est, size_t NZ, typename Store, size_t NRU, size_t NM, size_t NDU, bool HO>
constexpr size_t measured_of(const StateFeedback<NX, NU, T, NI, NY, Est, NZ, Store, NRU, NM, NDU, HO>* /*law*/) {
    return NM;
}

/// A reading with some structure.
double reading(int k) {
    return (0.5 * std::sin(0.037 * k)) + (((k / 150) % 2 == 0) ? 0.3 : -0.2);
}

/// The tick's arguments at sample @p k: the measurement, r, a measured d and a planned trajectory.
template<typename L>
struct Arguments {
    static constexpr size_t NIn = L::Core::NIn;
    static constexpr size_t NR = L::Core::NR;
    static constexpr size_t NM = measured_of(static_cast<L*>(nullptr));

    ColVec<NIn, double> m{};
    ColVec<NR, double>  r{};
    ColVec<NM, double>  d{};
    ColVec<2, double>   x_ref{0.01, 0.0};
    ColVec<1, double>   u_ref{0.02};

    explicit Arguments(int k) {
        for (size_t i = 0; i < NIn; ++i) {
            m(i) = reading(k + (7 * static_cast<int>(i)));
        }
        for (size_t i = 0; i < NR; ++i) {
            r(i) = ((k / 200) % 2 == 0) ? 0.2 : -0.1;
        }
        for (size_t i = 0; i < NM; ++i) {
            d(i) = 0.3 * std::sin(0.011 * k);
        }
    }
};

/**
 * One tick through step, the feedforward form every third: track(u, …) while another source drives
 * the plant, then step(…) with a clamp after the law, told with back_calculate.
 */
template<typename L>
double step_tick(L& law, int k, bool tracking) {
    const Arguments<L>      a(k);
    const bool              ff = (k % 3) == 0;
    constexpr size_t        NR = Arguments<L>::NR;
    constexpr size_t        NM = Arguments<L>::NM;
    const ColVec<1, double> other{0.4 * std::sin(0.05 * k)};
    if (tracking) {
        if constexpr (NR == 0 && NM == 0) {
            ff ? law.track(other, a.x_ref, a.u_ref, a.m) : law.track(other, a.m);
        } else if constexpr (NR == 0) {
            ff ? law.track(other, a.x_ref, a.u_ref, a.d, a.m) : law.track(other, a.d, a.m);
        } else if constexpr (NM == 0) {
            ff ? law.track(other, a.r, a.x_ref, a.u_ref, a.m) : law.track(other, a.r, a.m);
        } else {
            ff ? law.track(other, a.r, a.x_ref, a.u_ref, a.d, a.m) : law.track(other, a.r, a.d, a.m);
        }
        return other(0);
    }
    ColVec<1, double> u{};
    if constexpr (NR == 0 && NM == 0) {
        u = ff ? law.step(a.x_ref, a.u_ref, a.m) : law.step(a.m);
    } else if constexpr (NR == 0) {
        u = ff ? law.step(a.x_ref, a.u_ref, a.d, a.m) : law.step(a.d, a.m);
    } else if constexpr (NM == 0) {
        u = ff ? law.step(a.r, a.x_ref, a.u_ref, a.m) : law.step(a.r, a.m);
    } else {
        u = ff ? law.step(a.r, a.x_ref, a.u_ref, a.d, a.m) : law.step(a.r, a.d, a.m);
    }
    const ColVec<1, double> applied{std::clamp(u(0), -5.0, 5.0)};
    law.back_calculate(applied);
    return applied(0);
}

/// One tick in the ISR's order, the feedforward form every third: update, control, back_calculate after a clamp, predict.
template<typename L>
double isr_tick(L& law, int k) {
    const Arguments<L> a(k);
    const bool         ff = (k % 3) == 0;
    constexpr size_t   NR = Arguments<L>::NR;
    constexpr size_t   NM = Arguments<L>::NM;
    law.update(a.m);
    ColVec<1, double> u{};
    if constexpr (NR == 0 && NM == 0) {
        u = ff ? law.control(a.x_ref, a.u_ref) : law.control();
    } else if constexpr (NR == 0) {
        u = ff ? law.control(a.x_ref, a.u_ref, a.d) : law.control(a.d);
    } else if constexpr (NM == 0) {
        u = ff ? law.control(a.r, a.x_ref, a.u_ref) : law.control(a.r);
    } else {
        u = ff ? law.control(a.r, a.x_ref, a.u_ref, a.d) : law.control(a.r, a.d);
    }
    const ColVec<1, double> applied{std::clamp(u(0), -5.0, 5.0)};
    law.back_calculate(applied);
    law.predict();
    return applied(0);
}

/// The law with members and the Law on the design's constants send the same commands, to rounding, in both orders.
template<const auto& kRes>
void same_commands() {
    double worst = 0.0;
    double scale = 0.0;
    {
        StateFeedback     held{kRes};
        Law<kRes, double> constant;
        for (int k = 0; k < 600; ++k) {
            const bool   tracking = (k % 200) < 40; // another source drives the plant, then the law takes over
            const double a = step_tick(held, k, tracking);
            const double b = step_tick(constant, k, tracking);
            worst = std::max(worst, std::abs(a - b));
            scale = std::max(scale, std::abs(a));
        }
    }
    {
        StateFeedback     held{kRes};
        Law<kRes, double> constant;
        for (int k = 0; k < 600; ++k) {
            const double a = isr_tick(held, k);
            const double b = isr_tick(constant, k);
            worst = std::max(worst, std::abs(a - b));
            scale = std::max(scale, std::abs(a));
        }
    }
    CAPTURE(worst);
    CAPTURE(scale);
    CHECK(scale > 0.0);
    CHECK(worst <= 1e-11 * std::max(scale, 1.0));
}

/// The sizes the law had as one of the five structs (x86-64): with members in double and float, and the Law on constants in float.
template<const auto& kRes>
void sizes_are(size_t held, size_t held_float, size_t law_float) {
    CHECK(sizeof(decltype(StateFeedback{kRes})) == held);
    CHECK(sizeof(decltype(StateFeedback{kRes.template as<float>()})) == held_float);
    CHECK(sizeof(Law<kRes, float>) == law_float);
}

} // namespace

// Each design deduces StateFeedback, its kind's name the same type, the one that name spelled out was.
static_assert(std::is_same_v<decltype(StateFeedback{kLqr}), decltype(LQR{kLqr})> && std::is_same_v<decltype(LQR{kLqr}), LQR<2, 1, double>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqrTrack}), decltype(LQR{kLqrTrack})> && std::is_same_v<decltype(LQR{kLqrTrack}), LQR<2, 1, double, 2, 1>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqrDelayed}), decltype(LQR{kLqrDelayed})> && std::is_same_v<decltype(LQR{kLqrDelayed}), LQR<2, 1, double, 2, 0, Held, 1>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqrMeasured}), decltype(LQR{kLqrMeasured})> && std::is_same_v<decltype(LQR{kLqrMeasured}), LQR<2, 1, double, 2, 0, Held, 0, 1>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqrTrackMeasured}), decltype(LQR{kLqrTrackMeasured})> && std::is_same_v<decltype(LQR{kLqrTrackMeasured}), LQR<2, 1, double, 2, 1, Held, 0, 1>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqrRate}), decltype(LQR{kLqrRate})> && std::is_same_v<decltype(LQR{kLqrRate}), LQR<2, 1, double, 2, 0, Held, 0, 0, 1>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqrFaded}), decltype(LQR{kLqrFaded})> && std::is_same_v<decltype(LQR{kLqrFaded}), LQR<2, 1, double, 2, 0, Held, 0, 0, 0, true>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqrPlaced}), decltype(LQR{kLqrPlaced})> && std::is_same_v<decltype(LQR{kLqrPlaced}), LQR<2, 1, double>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqi}), decltype(LQI{kLqi})> && std::is_same_v<decltype(LQI{kLqi}), LQI<2, 1, 1, double, 1>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqiDelayed}), decltype(LQI{kLqiDelayed})> && std::is_same_v<decltype(LQI{kLqiDelayed}), LQI<2, 1, 1, double, 1, Held, 1>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqiMeasured}), decltype(LQI{kLqiMeasured})> && std::is_same_v<decltype(LQI{kLqiMeasured}), LQI<2, 1, 1, double, 1, Held, 0, 1>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqiRate}), decltype(LQI{kLqiRate})> && std::is_same_v<decltype(LQI{kLqiRate}), LQI<2, 1, 1, double, 1, Held, 0, 0, 1>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqiFaded}), decltype(LQI{kLqiFaded})> && std::is_same_v<decltype(LQI{kLqiFaded}), LQI<2, 1, 1, double, 1, Held, 0, 0, 0, true>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqg}), decltype(LQG{kLqg})> && std::is_same_v<decltype(LQG{kLqg}), LQG<2, 1, 1, double, PlantFilter>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqgTrack}), decltype(LQG{kLqgTrack})> && std::is_same_v<decltype(LQG{kLqgTrack}), LQG<2, 1, 1, double, PlantFilter, 1>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqgDelayed}), decltype(LQG{kLqgDelayed})> && std::is_same_v<decltype(LQG{kLqgDelayed}), LQG<2, 1, 1, double, DelayFilter, 0, Held, 1>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqgLate}), decltype(LQG{kLqgLate})> && std::is_same_v<decltype(LQG{kLqgLate}), LQG<2, 1, 1, double, DelayFilter>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqgMeasured}), decltype(LQG{kLqgMeasured})> && std::is_same_v<decltype(LQG{kLqgMeasured}), LQG<2, 1, 1, double, MeasuredFilter, 0, Held, 0, 1>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqgRate}), decltype(LQG{kLqgRate})> && std::is_same_v<decltype(LQG{kLqgRate}), LQG<2, 1, 1, double, PlantFilter, 0, Held, 0, 0, 1>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqgFaded}), decltype(LQG{kLqgFaded})> && std::is_same_v<decltype(LQG{kLqgFaded}), LQG<2, 1, 1, double, PlantFilter, 0, Held, 0, 0, 0, true>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqgPlaced}), decltype(LQG{kLqgPlaced})> && std::is_same_v<decltype(LQG{kLqgPlaced}), LQG<2, 1, 1, double, PlantFilter>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqgi}), decltype(LQGI{kLqgi})> && std::is_same_v<decltype(LQGI{kLqgi}), LQGI<2, 1, 1, 1, double, PlantFilter>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqgiDelayed}), decltype(LQGI{kLqgiDelayed})> && std::is_same_v<decltype(LQGI{kLqgiDelayed}), LQGI<2, 1, 1, 1, double, DelayFilter, Held, 1>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqgiMeasured}), decltype(LQGI{kLqgiMeasured})> && std::is_same_v<decltype(LQGI{kLqgiMeasured}), LQGI<2, 1, 1, 1, double, MeasuredFilter, Held, 0, 1>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqgiRate}), decltype(LQGI{kLqgiRate})> && std::is_same_v<decltype(LQGI{kLqgiRate}), LQGI<2, 1, 1, 1, double, PlantFilter, Held, 0, 0, 1>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqgiFaded}), decltype(LQGI{kLqgiFaded})> && std::is_same_v<decltype(LQGI{kLqgiFaded}), LQGI<2, 1, 1, 1, double, PlantFilter, Held, 0, 0, 0, true>>);
static_assert(std::is_same_v<decltype(StateFeedback{kOflqg}), decltype(OFLQG{kOflqg})> && std::is_same_v<decltype(OFLQG{kOflqg}), OFLQG<2, 1, 1, 1, double, LoadFilter>>);
static_assert(std::is_same_v<decltype(StateFeedback{kOflqgDelayed}), decltype(OFLQG{kOflqgDelayed})> && std::is_same_v<decltype(OFLQG{kOflqgDelayed}), OFLQG<2, 1, 1, 1, double, LoadDelay, Held, 1>>);
static_assert(std::is_same_v<decltype(StateFeedback{kOflqgLate}), decltype(OFLQG{kOflqgLate})> && std::is_same_v<decltype(OFLQG{kOflqgLate}), OFLQG<2, 1, 1, 1, double, LoadDelay>>);
static_assert(std::is_same_v<decltype(StateFeedback{kOflqgMeasured}), decltype(OFLQG{kOflqgMeasured})> && std::is_same_v<decltype(OFLQG{kOflqgMeasured}), OFLQG<2, 1, 1, 1, double, LoadMeasured, Held, 0, 1>>);
static_assert(std::is_same_v<decltype(StateFeedback{kOflqgRate}), decltype(OFLQG{kOflqgRate})> && std::is_same_v<decltype(OFLQG{kOflqgRate}), OFLQG<2, 1, 1, 1, double, LoadFilter, Held, 0, 0, 1>>);
static_assert(std::is_same_v<decltype(StateFeedback{kOflqgFaded}), decltype(OFLQG{kOflqgFaded})> && std::is_same_v<decltype(OFLQG{kOflqgFaded}), OFLQG<2, 1, 1, 1, double, LoadFilter, Held, 0, 0, 0, true>>);

// The scalar is the result's, and a gain K of its own is an LQR on the measured state.
static_assert(std::is_same_v<decltype(StateFeedback{kLqr.as<float>()}), LQR<2, 1, float>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqg.as<float>()}), LQG<2, 1, 1, float, FixedGainObserver<2, 1, 1, float>>>);
static_assert(std::is_same_v<decltype(StateFeedback{kLqr.K}), LQR<2, 1, double>> && std::is_same_v<decltype(LQR{kLqr.K}), LQR<2, 1, double>>);
static_assert(std::is_same_v<LQR<2, 1, double>, StateFeedback<2, 1, double>>);
static_assert(std::is_same_v<LQI<2, 1, 1, double, 1>, StateFeedback<2, 1, double, 1, 1>>);

// Each name takes its own kind's design and refuses every other's; StateFeedback takes them all.
static_assert(only<Res<kLqr>>(kLQR) && only<Res<kLqrTrack>>(kLQR) && only<Res<kLqrDelayed>>(kLQR) && only<Res<kLqrMeasured>>(kLQR));
static_assert(only<Res<kLqrTrackMeasured>>(kLQR) && only<Res<kLqrRate>>(kLQR) && only<Res<kLqrFaded>>(kLQR) && only<Res<kLqrPlaced>>(kLQR));
static_assert(only<Res<kLqi>>(kLQI) && only<Res<kLqiDelayed>>(kLQI) && only<Res<kLqiMeasured>>(kLQI) && only<Res<kLqiRate>>(kLQI) && only<Res<kLqiFaded>>(kLQI));
static_assert(only<Res<kLqg>>(kLQG) && only<Res<kLqgTrack>>(kLQG) && only<Res<kLqgDelayed>>(kLQG) && only<Res<kLqgLate>>(kLQG));
static_assert(only<Res<kLqgMeasured>>(kLQG) && only<Res<kLqgRate>>(kLQG) && only<Res<kLqgFaded>>(kLQG) && only<Res<kLqgPlaced>>(kLQG));
static_assert(only<Res<kLqgi>>(kLQGI) && only<Res<kLqgiDelayed>>(kLQGI) && only<Res<kLqgiMeasured>>(kLQGI) && only<Res<kLqgiRate>>(kLQGI) && only<Res<kLqgiFaded>>(kLQGI));
static_assert(only<Res<kOflqg>>(kOFLQG) && only<Res<kOflqgDelayed>>(kOFLQG) && only<Res<kOflqgLate>>(kOFLQG));
static_assert(only<Res<kOflqgMeasured>>(kOFLQG) && only<Res<kOflqgRate>>(kOFLQG) && only<Res<kOflqgFaded>>(kOFLQG));
static_assert(only<Matrix<1, 2>>(kLQR));                                                        // a gain K of its own
static_assert(!feedback_takes<design::LQRResult<2, 1>> && !lqr_takes<design::LQRResult<2, 1>>); // the bare gain runs on no runtime

// Named with its arguments, an alias does not form for another kind's filter, or without integrators.
static_assert(forms_lqg<PlantFilter> && forms_lqg<DelayFilter> && !forms_lqg<NoEstimator> && !forms_lqg<LoadFilter>);
static_assert(forms_lqgi<PlantFilter> && !forms_lqgi<NoEstimator> && !forms_lqgi<LoadFilter>);
static_assert(forms_oflqg<LoadFilter> && !forms_oflqg<PlantFilter> && !forms_oflqg<NoEstimator>);
static_assert(forms_lqi<1> && !forms_lqi<0> && forms_lqgi_with<1> && !forms_lqgi_with<0>);

// Law<kRes, T> runs the same type: its Runtime is the deduced law, each kind's name.
static_assert(std::is_same_v<Law<kLqr, float>::Runtime, LQR<2, 1, float>>);
static_assert(std::is_same_v<Law<kLqrTrack, double>::Runtime, decltype(LQR{kLqrTrack})>);
static_assert(std::is_same_v<Law<kLqrDelayed, double>::Runtime, decltype(LQR{kLqrDelayed})>);
static_assert(std::is_same_v<Law<kLqrMeasured, double>::Runtime, decltype(LQR{kLqrMeasured})>);
static_assert(std::is_same_v<Law<kLqrFaded, double>::Runtime, decltype(LQR{kLqrFaded})>);
static_assert(std::is_same_v<Law<kLqi, double>::Runtime, decltype(LQI{kLqi})>);
static_assert(std::is_same_v<Law<kLqiRate, double>::Runtime, decltype(LQI{kLqiRate})>);
static_assert(std::is_same_v<Law<kLqg, double>::Runtime, decltype(LQG{kLqg})>);
static_assert(std::is_same_v<Law<kLqgLate, double>::Runtime, decltype(LQG{kLqgLate})>);
static_assert(std::is_same_v<Law<kLqgPlaced, double>::Runtime, decltype(LQG{kLqgPlaced})>);
static_assert(std::is_same_v<Law<kLqgi, double>::Runtime, decltype(LQGI{kLqgi})>);
static_assert(std::is_same_v<Law<kLqgiMeasured, double>::Runtime, decltype(LQGI{kLqgiMeasured})>);
static_assert(std::is_same_v<Law<kOflqg, double>::Runtime, decltype(OFLQG{kOflqg})>);
static_assert(std::is_same_v<Law<kOflqgDelayed, double>::Runtime, decltype(OFLQG{kOflqgDelayed})>);
static_assert(std::is_base_of_v<StateFeedback<2, 1, float, 0, 2, NoEstimator, 0, detail::ConstantGains<kLqr, float>>, Law<kLqr, float>>); // on the design's constants

// What the law has follows from its design: a data member it does not carry is an empty type, a function does not compile.
template<typename L>
concept has_state_estimate = requires(const L& l) { l.state_estimate(); };
template<typename L>
concept has_disturbance_estimate = requires(const L& l) { l.disturbance_estimate(); };
template<typename L>
concept has_last_command = requires(const L& l) {
    l.last_command();
    l.applied();
};
template<typename L>
concept has_integrators = requires(L l) { l.integral.xi; };
template<typename L>
concept has_target = requires(L l) { l.target.target_d; };
template<typename L>
inline constexpr bool has_estimator = !std::is_empty_v<decltype(L{}.estimator)>;
template<typename L>
inline constexpr bool has_tracking_time = !std::is_empty_v<decltype(L{}.tracking_time)>;

// The limits it reports: an empty box for each it does not have.
template<typename L, size_t NI, size_t ND, size_t NDU>
inline constexpr bool limits_sized = std::is_same_v<decltype(L{}.limits()), design::StateFeedbackLimits<2, 1, double, NI, ND, NDU>>;

template<typename C, typename... A>
concept controls_with = requires(C c, const A&... a) { c.control(a...); };
template<typename C, typename... A>
concept steps_with = requires(C c, const A&... a) { c.step(a...); };
template<typename C, typename... A>
concept tracks_with = requires(C c, const A&... a) { c.track(a...); };

using V1 = ColVec<1, double>; // u, u_ref, y, r, d
using X2 = ColVec<2, double>; // x, x_ref

using Lqr = decltype(LQR{kLqr});
using LqrTrack = decltype(LQR{kLqrTrack});
using LqrDelayed = decltype(LQR{kLqrDelayed});
using LqrMeasured = decltype(LQR{kLqrMeasured});
using LqrRate = decltype(LQR{kLqrRate});
using Lqi = decltype(LQI{kLqi});
using LqiMeasured = decltype(LQI{kLqiMeasured});
using Lqg = decltype(LQG{kLqg});
using LqgTrack = decltype(LQG{kLqgTrack});
using LqgMeasured = decltype(LQG{kLqgMeasured});
using Lqgi = decltype(LQGI{kLqgi});
using Oflqg = decltype(OFLQG{kOflqg});
using OflqgMeasured = decltype(OFLQG{kOflqgMeasured});

static_assert(!has_state_estimate<Lqr> && !has_last_command<Lqr> && !has_integrators<Lqr> && !has_target<Lqr> && !has_estimator<Lqr> && !has_tracking_time<Lqr>);
static_assert(has_target<LqrTrack> && has_target<LqrMeasured> && !has_integrators<LqrTrack> && !has_last_command<LqrTrack>);
static_assert(has_last_command<LqrDelayed> && !has_state_estimate<LqrDelayed>); // it keeps the inputs still on their way
static_assert(has_integrators<Lqi> && has_tracking_time<Lqi> && has_last_command<Lqi> && !has_state_estimate<Lqi> && !has_target<Lqi> && has_target<LqiMeasured>);
static_assert(has_state_estimate<Lqg> && has_estimator<Lqg> && !has_disturbance_estimate<Lqg> && !has_integrators<Lqg> && !has_target<Lqg> && has_target<LqgTrack>);
static_assert(has_state_estimate<Lqgi> && has_integrators<Lqgi> && !has_disturbance_estimate<Lqgi> && !has_target<Lqgi>);
static_assert(has_state_estimate<Oflqg> && has_disturbance_estimate<Oflqg> && has_target<Oflqg> && !has_integrators<Oflqg> && !has_tracking_time<Oflqg>);
static_assert(limits_sized<Lqr, 0, 0, 0> && limits_sized<LqrRate, 0, 0, 1> && limits_sized<Lqi, 1, 0, 0> && limits_sized<Lqgi, 1, 0, 0> && limits_sized<Oflqg, 0, 1, 0>);

// The tick each law takes: no r on a regulator, d on a measured disturbance's, and a track beside each step.
static_assert(controls_with<Lqr> && steps_with<Lqr, X2> && tracks_with<Lqr, V1, X2> && !controls_with<Lqr, V1>);
static_assert(controls_with<LqrMeasured, V1> && steps_with<LqrMeasured, V1, X2> && tracks_with<LqrMeasured, V1, V1, X2> && !controls_with<LqrMeasured>);
static_assert(controls_with<LqrMeasured, X2, V1, V1> && steps_with<LqrMeasured, X2, V1, V1, X2> && tracks_with<LqrMeasured, V1, X2, V1, V1, X2>);
static_assert(controls_with<Lqi, V1> && steps_with<Lqi, V1, X2> && tracks_with<Lqi, V1, V1, X2> && !controls_with<Lqi>);
static_assert(controls_with<LqiMeasured, V1, V1> && steps_with<LqiMeasured, V1, V1, X2> && tracks_with<LqiMeasured, V1, V1, V1, X2> && !controls_with<LqiMeasured, V1>);
static_assert(controls_with<LqgMeasured, V1> && steps_with<LqgMeasured, V1, V1> && !steps_with<LqgMeasured, V1, X2>);
static_assert(controls_with<OflqgMeasured, V1, V1> && controls_with<OflqgMeasured, V1, X2, V1, V1> && !controls_with<OflqgMeasured, V1>);

TEST_SUITE("StateFeedback") {
    TEST_CASE("each LQ law is the size it was as one of five structs") {
        //! Recorded on x86-64 on the base this was built on, where LQR, LQI, LQG, LQGI and OFLQG were
        //! five structs over the same core: the type is the core's, so nothing is added.
        if constexpr (sizeof(void*) == 8) {
            sizes_are<kLqr>(88, 48, 20);
            sizes_are<kLqrTrack>(120, 64, 20);
            sizes_are<kLqrDelayed>(136, 80, 40);
            sizes_are<kLqrMeasured>(120, 64, 20);
            sizes_are<kLqrTrackMeasured>(136, 72, 20);
            sizes_are<kLqrRate>(112, 64, 32);
            sizes_are<kLqrFaded>(104, 56, 24);
            sizes_are<kLqrPlaced>(88, 48, 20);
            sizes_are<kLqi>(216, 112, 68);
            sizes_are<kLqiDelayed>(240, 128, 72);
            sizes_are<kLqiMeasured>(248, 128, 68);
            sizes_are<kLqiRate>(232, 120, 76);
            sizes_are<kLqiFaded>(232, 120, 72);
            sizes_are<kLqg>(240, 128, 44);
            sizes_are<kLqgTrack>(272, 144, 44);
            sizes_are<kLqgDelayed>(320, 168, 48);
            sizes_are<kLqgLate>(312, 168, 48);
            sizes_are<kLqgMeasured>(320, 168, 56);
            sizes_are<kLqgRate>(256, 136, 52);
            sizes_are<kLqgFaded>(264, 136, 52);
            sizes_are<kLqgPlaced>(240, 128, 44);
            sizes_are<kLqgi>(336, 176, 80);
            sizes_are<kLqgiDelayed>(416, 216, 84);
            sizes_are<kLqgiMeasured>(416, 216, 92);
            sizes_are<kLqgiRate>(352, 184, 88);
            sizes_are<kLqgiFaded>(360, 184, 88);
            sizes_are<kOflqg>(424, 224, 80);
            sizes_are<kOflqgDelayed>(552, 288, 92);
            sizes_are<kOflqgLate>(544, 280, 92);
            sizes_are<kOflqgMeasured>(512, 264, 96);
            sizes_are<kOflqgRate>(440, 232, 88);
            sizes_are<kOflqgFaded>(448, 232, 88);
        }
    }

    TEST_CASE("StateFeedback law{res} of each kind sends the commands of its Law on constants, through step, track and the ISR's order") {
        same_commands<kLqr>();
        same_commands<kLqrTrack>();
        same_commands<kLqrDelayed>();
        same_commands<kLqrMeasured>();
        same_commands<kLqrTrackMeasured>();
        same_commands<kLqrRate>();
        same_commands<kLqrFaded>();
        same_commands<kLqrPlaced>();
        same_commands<kLqi>();
        same_commands<kLqiDelayed>();
        same_commands<kLqiMeasured>();
        same_commands<kLqiRate>();
        same_commands<kLqiFaded>();
        same_commands<kLqg>();
        same_commands<kLqgTrack>();
        same_commands<kLqgDelayed>();
        same_commands<kLqgLate>();
        same_commands<kLqgMeasured>();
        same_commands<kLqgRate>();
        same_commands<kLqgFaded>();
        same_commands<kLqgPlaced>();
        same_commands<kLqgi>();
        same_commands<kLqgiDelayed>();
        same_commands<kLqgiMeasured>();
        same_commands<kLqgiRate>();
        same_commands<kLqgiFaded>();
        same_commands<kOflqg>();
        same_commands<kOflqgDelayed>();
        same_commands<kOflqgLate>();
        same_commands<kOflqgMeasured>();
        same_commands<kOflqgRate>();
        same_commands<kOflqgFaded>();
    }

    TEST_CASE("a Law of a design at init holds a StateFeedback, and takes a redesign of the same kind") {
        static const auto    at_init = design::state_feedback(plant(), kTs, kW, kNoise);
        Law<at_init, double> law;
        static_assert(std::is_base_of_v<decltype(LQG{kLqg}), Law<at_init, double>>);
        law = LQG{design::state_feedback(plant(), kTs, {.Q = kW.Q * 10.0, .R = 1.0}, kNoise)};
        CHECK(law.valid());
        CHECK(law.K(0, 0) != doctest::Approx(kLqg.K(0, 0)));
    }
}
