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
static_assert(damp::constructible_from<SmithPredictor<16, float>, design::SmithPredictorResult<float>>);

// ---------------------------------------------------------------------------------------------
// Known gaps (#142): each PR deletes the lines it closes and moves the law up.
// ---------------------------------------------------------------------------------------------

// SMCController: a public u_lim, no track, no last_command.
static_assert(SisoController<SMCController<float>, float> && ReportsValidity<SMCController<float>>);
static_assert(!HoldsLimits<SMCController<float>>);
static_assert(!SisoTracking<SMCController<float>, float>);
static_assert(!ReportsCommand<SMCController<float>, float>);

// STSMCController: a public u_lim, no track.
static_assert(SisoController<STSMCController<float>, float> && ReportsCommand<STSMCController<float>, float> && ReportsValidity<STSMCController<float>>);
static_assert(!HoldsLimits<STSMCController<float>>);
static_assert(!SisoTracking<STSMCController<float>, float>);

// ADRCController: a public u_lim, no track.
static_assert(SisoController<ADRCController<2, float>, float> && ReportsCommand<ADRCController<2, float>, float> && ReportsValidity<ADRCController<2, float>>);
static_assert(!HoldsLimits<ADRCController<2, float>>);
static_assert(!SisoTracking<ADRCController<2, float>, float>);

// SmithPredictor: no limits() (its primary's), no track, no last_command.
static_assert(SisoController<SmithPredictor<16, float>, float> && ReportsValidity<SmithPredictor<16, float>>);
static_assert(!HoldsLimits<SmithPredictor<16, float>>);
static_assert(!SisoTracking<SmithPredictor<16, float>, float>);
static_assert(!ReportsCommand<SmithPredictor<16, float>, float>);

// SwitchedController: no limits(), no track, no valid().
static_assert(SisoController<SwitchedPI, float> && ReportsCommand<SwitchedPI, float>);
static_assert(!HoldsLimits<SwitchedPI>);
static_assert(!SisoTracking<SwitchedPI, float>);
static_assert(!ReportsValidity<SwitchedPI>);

// Cascade: not a SisoController (step(r, y_outer, y_inner)); its contract is limits(), track(u, r, y_outer, y_inner), last_command() and valid().
static_assert(ReportsCommand<CascadePI, float>);
static_assert(!HoldsLimits<CascadePI>);
static_assert(!CascadeTracking<CascadePI, float>);
static_assert(!ReportsValidity<CascadePI>);

// ExtremumSeekingController: an optimizer (step(J, valid)), not a SisoController; its contract is limits(), last_command() and valid().
static_assert(ReportsValidity<ExtremumSeekingController<float>>);
static_assert(!HoldsLimits<ExtremumSeekingController<float>>);
static_assert(!ReportsCommand<ExtremumSeekingController<float>, float>);

TEST_SUITE("Runtime contract") {
    TEST_CASE("every SISO runtime is checked against the contract at compile time") {
        CHECK(SisoLaw<PIDController<float>, float>); // the static_asserts above are the test
    }
}
