// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file test_contract.cpp
 * @brief The SISO runtime contract (#142), checked for every runtime at compile time
 *
 * Every SISO law holds its limits behind limits()/set_limits(), hands over with track(u, r, y),
 * reports its command with last_command() and whether it runs with valid(), and is built from its
 * design's result (concepts.hpp: SisoLaw). A law on the whole contract is listed under "Full
 * contract". What a law still lacks is listed under "Known gaps" as a static_assert that it lacks
 * it. A PR that closes a gap breaks that line, so it moves the law up rather than leave the list
 * stale, and a new law can't ship off the contract without saying so here.
 */

#include "damp/concepts.hpp"
#include "damp/controllers/adrc.hpp"
#include "damp/controllers/composition.hpp"
#include "damp/controllers/esc.hpp"
#include "damp/controllers/lead_lag.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/controllers/pr.hpp"
#include "damp/controllers/repetitive.hpp"
#include "damp/controllers/smc.hpp"
#include "damp/controllers/smith_predictor.hpp"
#include "damp/controllers/stsmc.hpp"
#include "damp/design/pid_design.hpp"
#include "doctest.h"

using namespace damp;

namespace {

using PI = PIController<float>;
using CascadePI = Cascade<PI, PI, float>;
using SwitchedPI = SwitchedController<PI, PI, PI, float>;

/** @brief A cascade's handover: track(u, r, y_outer, y_inner), as its step takes both readings. */
template<typename C, typename T>
concept CascadeTracking = requires(C c, T u, T r, T y_outer, T y_inner) { c.track(u, r, y_outer, y_inner); };

/** @brief An optimizer's handover: track(u), its operating point started from the input applied, as its step takes no reference. */
template<typename C, typename T>
concept OptimizerTracking = requires(C c, T u) { c.track(u); };

using ESC = ExtremumSeekingController<float>;
/** @brief The handover of a law with an error-form step(e): track(u, e). */
template<typename C, typename T>
concept ErrorTracking = requires(C c, T u, T e) { c.track(u, e); };

} // namespace

// ---------------------------------------------------------------------------------------------
// Full contract
// ---------------------------------------------------------------------------------------------

static_assert(SisoLaw<PIDController<float>, float>);
static_assert(SisoLaw<PIDController<double>, double>);
static_assert(SisoLaw<PIController<float>, float>);
static_assert(SisoLaw<PDController<float>, float>);
static_assert(SisoLaw<PController<float>, float>);
static_assert(SisoLaw<SMCController<float>, float>);
static_assert(SisoLaw<SMCController<double>, double>);
static_assert(SisoLaw<STSMCController<float>, float>);
static_assert(SisoLaw<ADRCController<1, float>, float>);
static_assert(SisoLaw<ADRCController<2, float>, float>);
static_assert(SisoLaw<SmithPredictor<16, float>, float>);
static_assert(SisoLaw<SmithPredictor<16, float, SMCController<float>>, float>);
static_assert(SisoLaw<SwitchedPI, float>);

// Cascade: not a SisoController (step(r, y_outer, y_inner)); its contract is limits() (the inner law's), track(u, r, y_outer, y_inner), last_command() and valid().
static_assert(HoldsLimits<CascadePI> && CascadeTracking<CascadePI, float> && ReportsCommand<CascadePI, float> && ReportsValidity<CascadePI>);

// ExtremumSeekingController: an optimizer (step(J, valid)), not a SisoController; its contract is limits(), track(u), last_command() and valid().
static_assert(HoldsLimits<ESC> && OptimizerTracking<ESC, float> && ReportsCommand<ESC, float> && ReportsValidity<ESC>);
static_assert(!SisoController<ESC, float>);
static_assert(SisoLaw<PRController<float>, float> && ErrorTracking<PRController<float>, float>);
static_assert(SisoLaw<PRController<double, 3>, double>); // a bank, design::Harmonics
static_assert(SisoLaw<LeadLagController<float>, float> && ErrorTracking<LeadLagController<float>, float>);
static_assert(SisoLaw<LeadLagController<double, 2>, double>); // design::lead_lag's series

// RepetitiveController: a plug-in add-on (step(e)), clamped with the law it adds to, so it holds no limits of its own.
static_assert(ErrorTracking<RepetitiveController<64, float>, float>);
static_assert(ReportsCommand<RepetitiveController<64, float>, float> && ReportsValidity<RepetitiveController<64, float>>);

// ---------------------------------------------------------------------------------------------
// Built from a design: every runtime that is a law (not a composite) takes its design's result.
// ---------------------------------------------------------------------------------------------

static_assert(damp::constructible_from<PIDController<float>, design::PIDResult<float>>);
static_assert(damp::constructible_from<PRController<float>, design::PRResult<float>>);
static_assert(damp::constructible_from<PRController<float, 3>, design::PRResult<float, 3>>);
static_assert(damp::constructible_from<LeadLagController<float>, design::LeadLagResult<float>>);
static_assert(damp::constructible_from<LeadLagController<float, 2>, design::LeadLagResult<float, 2>>);
static_assert(damp::constructible_from<RepetitiveController<64, float>, design::RepetitiveResult<float, 0>>);
static_assert(damp::constructible_from<SMCController<float>, design::SMCResult<float>>);
static_assert(damp::constructible_from<STSMCController<float>, design::STSMCResult<float>>);
static_assert(damp::constructible_from<ADRCController<2, float>, design::ADRCResult<2, float>>);
static_assert(damp::constructible_from<ExtremumSeekingController<float>, design::ESCResult<float>>);
static_assert(damp::constructible_from<SmithPredictor<16, float>, design::SmithPredictorResult<float, 16>>);

// ---------------------------------------------------------------------------------------------
// A constexpr design's Law holds its gains as constants of the type: it is smaller than its
// Runtime, which holds them as members, and where the Runtime retunes with set_gains the Law has
// none. Every family's Law has gains to drop, so each is strictly smaller; the repetitive and Smith
// runtimes have no set_gains to lose (a repetitive period is a buffer size; a Smith predictor
// retunes through its primary, which on the Law is the PID's constant-gains runtime).
// ---------------------------------------------------------------------------------------------

namespace {

inline constexpr double kPi = 3.14159265358979323846;

inline constexpr auto kPid = design::pid(design::Parallel{.Kp = 2.0, .Ki = 5.0, .Kd = 0.1}, 0.01, design::DerivativeFilter{.Tf = 0.02});
inline constexpr auto kPr = design::pr({.Kp = 0.5, .Ki = 200.0, .w0 = 2.0 * kPi * 50.0, .wc = 5.0}, 1e-4, design::Harmonics{{1, 5, 7}});
inline constexpr auto kLeadLag = design::lead_lag({.phase_max = 0.6, .wc = 20.0, .dc_gain_boost = 10.0}, 1e-3);
inline constexpr auto kRepetitive = design::repetitive_binomial<2>({.w0 = 2.0 * kPi * 50.0, .gain = 0.5}, 1e-4);
inline constexpr auto kSmc = design::smc({.lambda = 10.0, .k = 5.0, .b0 = 1.0}, 0.01);
inline constexpr auto kStsmc = design::stsmc({.disturbance_bound = 2.0, .lambda = 5.0}, 1e-3);
inline constexpr auto kAdrc = design::adrc<2>({.wc = 10.0, .wo = 50.0, .b0 = 2.0}, 1e-3);
inline constexpr auto kEsc = design::esc({.dither_amplitude = 0.05, .dither_omega = 2.0 * kPi * 40.0, .gain = 60.0, .lp_omega = 2.0 * kPi * 8.0}, 1e-3);
inline constexpr auto kSmith = design::smith_predictor(design::FOPDT{.K = 1.0, .tau = 2.0}, 0.1, InputDelay<12>{1.05}, design::SIMC{.tau_c = 1.0});

/** @brief A law whose gains set_gains changes at run time: it takes back the gains it reports. */
template<typename L>
concept Retunable = requires(L& l) { l.set_gains(l.gains()); };

/** @brief The Law of @p Res in float is smaller than its Runtime, and retunes only where its Runtime does not: none of it is held. */
template<const auto& Res>
inline constexpr bool gains_in_type = (sizeof(Law<Res, float>) < sizeof(typename Law<Res, float>::Runtime)) && !Retunable<Law<Res, float>>;

} // namespace

static_assert(gains_in_type<kPid> && Retunable<Law<kPid, float>::Runtime>);
static_assert(gains_in_type<kPr> && Retunable<Law<kPr, float>::Runtime>);
static_assert(gains_in_type<kLeadLag> && Retunable<Law<kLeadLag, float>::Runtime>);
static_assert(gains_in_type<kSmc> && Retunable<Law<kSmc, float>::Runtime>);
static_assert(gains_in_type<kStsmc> && Retunable<Law<kStsmc, float>::Runtime>);
static_assert(gains_in_type<kAdrc> && Retunable<Law<kAdrc, float>::Runtime>);
static_assert(gains_in_type<kEsc> && Retunable<Law<kEsc, float>::Runtime>);
static_assert(gains_in_type<kRepetitive> && !Retunable<Law<kRepetitive, float>::Runtime>);
static_assert(gains_in_type<kSmith> && !Retunable<Law<kSmith, float>::Runtime>);
static_assert(!Retunable<damp::remove_cvref_t<decltype(damp::declval<Law<kSmith, float>&>().primary())>>); // the primary's gains are the type's too
static_assert(Retunable<damp::remove_cvref_t<decltype(damp::declval<Law<kSmith, float>::Runtime&>().primary())>>);

// The Law of a constexpr design is still on the contract the Runtime is on.
static_assert(SisoLaw<Law<kPid, float>, float> && SisoLaw<Law<kPr, float>, float> && SisoLaw<Law<kLeadLag, float>, float>);
static_assert(SisoLaw<Law<kSmc, float>, float> && SisoLaw<Law<kStsmc, float>, float> && SisoLaw<Law<kAdrc, float>, float>);
static_assert(SisoLaw<Law<kSmith, float>, float>);
static_assert(HoldsLimits<Law<kEsc, float>> && OptimizerTracking<Law<kEsc, float>, float> && ReportsCommand<Law<kEsc, float>, float> && ReportsValidity<Law<kEsc, float>>);
static_assert(ErrorTracking<Law<kRepetitive, float>, float> && ReportsCommand<Law<kRepetitive, float>, float> && ReportsValidity<Law<kRepetitive, float>>);

// ---------------------------------------------------------------------------------------------
// Known gaps: none. A runtime that misses part of the contract is listed here, as a
// static_assert(!…) of what it lacks, until a PR closes it and moves the runtime up.
// ---------------------------------------------------------------------------------------------

TEST_SUITE("Runtime contract") {
    TEST_CASE("every SISO runtime is checked against the contract at compile time") {
        CHECK(SisoLaw<PIDController<float>, float>); // the static_asserts above are the test
    }
}
