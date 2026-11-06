// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// One construction path for the non-LQ runtimes: a design result in any scalar converts
// implicitly (`Law law = res;`), the runtime type deduces from the result, and
// `static constinit` deploys a double design without an explicit .as<float>().

#include <numbers>
#include <type_traits>

#include "damp/controllers/adrc.hpp"
#include "damp/controllers/esc.hpp"
#include "damp/controllers/lead_lag.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/controllers/pr.hpp"
#include "damp/controllers/repetitive.hpp"
#include "damp/controllers/smc.hpp"
#include "damp/controllers/stsmc.hpp"
#include "damp/design/pid_design.hpp"
#include "damp/estimation/excitation/chirp.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/estimation/luenberger.hpp"
#include "damp/filters/biquad.hpp"
#include "damp/filters/first_order_section.hpp"
#include "damp/filters/iir_design.hpp"
#include "damp/filters/lowpass.hpp"
#include "damp/matrix/core.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

constexpr double kTs = 1e-3;

// Double designs, deployed in float below with no .as<float>().
constexpr auto kPid = design::pid(2.0, 5.0, 0.1).discretize(kTs);
constexpr auto kAdrc = design::adrc<2>(10.0, 50.0, 1.0, kTs);
constexpr auto kSmc = design::smc(10.0, 8.0, 1.0, kTs);
constexpr auto kStsmc = design::stsmc(4.0, kTs, 10.0);
constexpr auto kPr = design::pr({.Kp = 0.0, .Ki = 10.0, .w0 = 2.0 * std::numbers::pi * 50.0, .wc = 6.0}, kTs);
constexpr auto kLead = design::lead(0.6, 20.0, kTs);
constexpr auto kEsc = design::esc({.dither_amplitude = 0.1, .dither_freq_hz = 50.0, .gain = 20.0}, kTs);
constexpr auto kRep = design::repetitive(1000.0, 50.0);
constexpr auto kLowPass = design::lowpass<2>(50.0, kTs);
constexpr auto kNotch = design::notch(50.0, 2.0, kTs);
constexpr auto kHighPass = design::highpass_1st(5.0, kTs);
constexpr auto kChirp = design::chirp(design::ChirpConfig<double>{.amplitude = 1.0, .f_start_hz = 1.0, .f_end_hz = 10.0, .duration_s = 1.0}, kTs);
constexpr auto kObserver = design::luenberger(
    StateSpace<2, 1, 1>{
        .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
        .B = Matrix<2, 1>{{0.005}, {0.1}},
        .C = Matrix<1, 2>{{1.0, 0.0}},
        .Ts = 0.1,
    },
    {0.3, 0.4}
);

static constinit PIDController<float>               pid = kPid;
static constinit ContinuousPID<float>               cpid = design::pid(2.0, 5.0, 0.1);
static constinit ADRCController<2, float>           adrc = kAdrc;
static constinit SMCController<float>               smc = kSmc;
static constinit STSMCController<float>             stsmc = kStsmc;
static constinit PRController<float>                pr = kPr;
static constinit LeadLagController<float>           lead = kLead;
static constinit ExtremumSeekingController<float>   esc = kEsc;
static constinit RepetitiveController<64, float, 0> rep = kRep;
static constinit LowPass<2, float>                  lowpass = kLowPass;
static constinit Biquad<float>                      notch = kNotch;
static constinit FirstOrderSection<float>           highpass = kHighPass;
static constinit Chirp<float>                       chirp = kChirp;
static constinit Luenberger<2, 1, 1, float>         observer = kObserver;

} // namespace

TEST_SUITE("assign from design") {
    TEST_CASE("a double design converts into a float runtime under constinit") {
        CHECK(pid.Kp == doctest::Approx(2.0f));
        CHECK(cpid.Ki == doctest::Approx(5.0f));
        CHECK(observer.valid());
        CHECK(static_cast<double>(lowpass(1.0f)) == doctest::Approx(kLowPass.b[0]).epsilon(1e-6));
        CHECK(static_cast<double>(highpass(1.0f)) == doctest::Approx(kHighPass.b0).epsilon(1e-6));
        CHECK(static_cast<double>(notch(1.0f)) == doctest::Approx(kNotch.b0).epsilon(1e-6));
        const Chirp chirp_d{kChirp}; // deduces Chirp<double>
        static_assert(std::is_same_v<std::remove_const_t<decltype(chirp_d)>, Chirp<double>>);
        CHECK(static_cast<double>(chirp.step(0.25f)) == doctest::Approx(chirp_d.step(0.25)).epsilon(1e-5));
        (void)adrc;
        (void)smc;
        (void)stsmc;
        (void)pr;
        (void)lead;
        (void)esc;
        (void)rep;
    }

    TEST_CASE("the runtime type deduces from the result's scalar") {
        PIDController             pd{kPid};
        ContinuousPID             cd{design::pid(1.0, 1.0, 0.0)};
        ADRCController            ad{kAdrc};
        SMCController             sd{kSmc};
        STSMCController           td{kStsmc};
        PRController              prd{kPr};
        LeadLagController         ld{kLead};
        ExtremumSeekingController ed{kEsc};
        LowPass                   lp{kLowPass};
        Biquad                    bq{kNotch};
        FirstOrderSection         hp{kHighPass};
        Luenberger                ob{kObserver};
        static_assert(std::is_same_v<decltype(pd), PIDController<double>>);
        static_assert(std::is_same_v<decltype(cd), ContinuousPID<double>>);
        static_assert(std::is_same_v<decltype(ad), ADRCController<2, double>>);
        static_assert(std::is_same_v<decltype(sd), SMCController<double>>);
        static_assert(std::is_same_v<decltype(td), STSMCController<double>>);
        static_assert(std::is_same_v<decltype(prd), PRController<double>>);
        static_assert(std::is_same_v<decltype(ld), LeadLagController<double>>);
        static_assert(std::is_same_v<decltype(ed), ExtremumSeekingController<double>>);
        static_assert(std::is_same_v<decltype(lp), LowPass<2, double>>);
        static_assert(std::is_same_v<decltype(bq), Biquad<double>>);
        static_assert(std::is_same_v<decltype(hp), FirstOrderSection<double>>);
        static_assert(std::is_same_v<decltype(ob), Luenberger<2, 1, 1, double>>);
        CHECK(pd.Kp == doctest::Approx(2.0));
    }

    TEST_CASE("if (result) reads a design result's success") {
        constexpr StateSpace<1, 1, 1, double, 1, 1> plant{
            .A = Matrix<1, 1>{{0.95}},
            .B = Matrix<1, 1>{{1.0}},
            .C = Matrix<1, 1>{{1.0}},
            .G = Matrix<1, 1>{{1.0}},
            .H = Matrix<1, 1>{{1.0}},
            .Ts = 1.0,
        };

        const auto lqr = design::lqr(plant, {.Q = Matrix<1, 1>{{1.0}}, .R = Matrix<1, 1>{{0.1}}});
        const auto lqr_bad = design::lqr(plant, {.Q = Matrix<1, 1>{{-1.0}}, .R = Matrix<1, 1>{{0.1}}}); // Q not PSD
        REQUIRE(lqr.success);
        REQUIRE_FALSE(lqr_bad.success);
        CHECK(static_cast<bool>(lqr));
        CHECK_FALSE(static_cast<bool>(lqr_bad));

        const auto kalman = design::kalman(plant, {.process = Matrix<1, 1>{{0.1}}, .measurement = Matrix<1, 1>{{0.25}}});
        REQUIRE(kalman.success);
        CHECK(static_cast<bool>(kalman));
        CHECK_FALSE(static_cast<bool>(decltype(kalman){}));

        constexpr auto pid = design::ziegler_nichols(10.0, 2.0);
        static_assert(pid); // contextual conversion under constant evaluation
        static_assert(!design::ziegler_nichols(10.0, 0.0));

        const auto mpc = design::mpc<4, 2>(plant, design::MPCWeights<1, 1>{}, {});
        REQUIRE(mpc.success);
        CHECK(static_cast<bool>(mpc));
        CHECK_FALSE(static_cast<bool>(decltype(mpc){}));

        int designed = 0;
        if (lqr) {
            ++designed;
        }
        if (kalman) {
            ++designed;
        }
        if (pid) {
            ++designed;
        }
        if (mpc) {
            ++designed;
        }
        if (lqr_bad) {
            ++designed;
        }
        CHECK(designed == 4);
    }

    TEST_CASE("a refused design is not a constant-evaluated runtime") {
        // A refused cutoff: the conversion compiles at run time (zero taps) but not under constant
        // evaluation, where the constructor calls refuse_failed_design(). That half cannot be a
        // static_assert (the constructor exists; only its evaluation is ill-formed):
        // static constinit LowPass<1, float> refused{design::lowpass<1>(0.0, kTs)}; // ill-formed
        static_assert(!design::lowpass<1>(0.0, kTs).success);
        LowPass<1, float> lp = design::lowpass<1>(0.0, kTs);
        CHECK(lp(1.0f) == 0.0f);
    }
}
