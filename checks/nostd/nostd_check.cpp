// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Design and tick with no C++ standard headers: tools/contract_check.sh nostd builds this with
// -nostdinc++ and runs it. In a subfolder so the checks/ asm glob does not pick it up. ETL_NO_STL
// comes from etl_profile.h beside this file, the setup README.md gives.
#define DAMP_BACKEND_ETL
#define DAMP_MATH_BACKEND_FREESTANDING

#include "damp/control.hpp"
#include "damp/motor/foc.hpp"
#include "damp/motor/modulation.hpp"
#include "damp/motor/spm.hpp"

using namespace damp;

namespace {

// The optional this configuration provides: a non-POD value, empty and engaged, in constant evaluation.
struct NonPod {
    int v{7};
};
constexpr optional<NonPod> kEmpty{};
constexpr optional<NonPod> kNull{nullopt};
constexpr optional<NonPod> kFull{NonPod{3}};
constexpr optional<NonPod> kInPlace{in_place, NonPod{4}}; // how the solves build their result
static_assert(!kEmpty.has_value() && !kNull && kFull.has_value() && kFull->v == 3 && (*kFull).v == 3);
static_assert(kInPlace.has_value() && kInPlace->v == 4);
static_assert(kEmpty.value_or(NonPod{9}).v == 9 && kFull.value_or(NonPod{9}).v == 3);
static_assert([] {
    optional<NonPod> o{NonPod{1}};
    o.reset();
    const bool cleared = !o.has_value();
    o.emplace(NonPod{5});
    return cleared && o.value().v == 5;
}());

constexpr auto sys = StateSpace{.A = Matrix<2, 2>{{0.0, 1.0}, {9.81, -0.1}}, .B = Matrix<2, 1>{{0.0}, {1.0}}, .C = Matrix<1, 2>{{1.0, 0.0}}};

constexpr auto kLqr = design::lqr(sys, 0.01, {.Q = Matrix<2, 2>::identity() * 10.0, .R = 1.0});
constexpr auto kEso = design::eso<2>(40.0, 20.0, 1e-3);
constexpr auto kPid = design::pid(2.0, 5.0, 0.1).discretize(0.01);
constexpr StateSpace<2, 1, 1, double, 1, 1> noisy{.A = sys.A, .B = sys.B, .C = sys.C, .G = Matrix<2, 1>{{0.0}, {1.0}}, .H = Matrix<1, 1>{{1.0}}};
constexpr auto kLqg = design::lqg(noisy, 0.01, {.Q = Matrix<2, 2>::identity() * 10.0, .R = 1.0}, design::KalmanNoise<1, 1>{.process = 0.5, .measurement = 1e-4});
static_assert(kLqr.success && kEso.success && kLqg.success);

static constinit LQR           lqr{kLqr.as<float>()};
static constinit ESO           eso{kEso.as<float>()};
static constinit PIDController pid{kPid.as<float>()};
static constinit LQG           lqg{kLqg.as<float>()};
static constinit Law<kLqr, float> lqr_on_constants; // the same law, its gains constants of the type

} // namespace

int main() {
    const float u = lqr.step(ColVec<2, float>{0.1F, 0.0F})(0, 0);
    const float uc = lqr_on_constants.step(ColVec<2, float>{0.1F, 0.0F})(0, 0);
    (void)eso.step(0.1F, u);
    const float v = pid.control(1.0F, 0.0F);
    const float w = lqg.step(ColVec<1, float>{0.1F})(0, 0);
    const bool  ok = damp::isfinite(u) && (u < 0.0F) && damp::isfinite(eso.disturbance()) && damp::isfinite(v) && (v > 0.0F) && damp::isfinite(w) && (damp::abs(uc - u) <= 1e-5F * (1.0F + damp::abs(u)));
    return ok ? 0 : 1;
}
