// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>

#include "damp/backend.hpp"
#include "damp/controllers/esc.hpp"
#include "damp/design/status.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {
// Run the ESC loop on a static scalar objective J(u); return the converged û.
template<typename Obj>
double seek(const design::ESCResult<double>& art, Obj J, int steps) {
    ExtremumSeekingController<double> esc(art);
    double                            u = esc.input();
    for (int k = 0; k < steps; ++k) {
        u = esc.step(J(u));
    }
    return esc.operating_point();
}

/** @brief The angular frequency [rad/s] of @p hz: the cases were written in Hz, as the spec was. */
template<typename T>
constexpr T w(T hz) {
    return T{2} * std::numbers::pi_v<T> * hz;
}

/** @brief A NaN that reaches the guard at run time, from bits through a volatile (tests/test_nan_guards.cpp). */
double runtime_nan() {
    volatile std::uint64_t bits = 0x7FF8000000000000ULL;
    return std::bit_cast<double>(static_cast<std::uint64_t>(bits));
}
} // namespace

TEST_SUITE("esc") {
    TEST_CASE("synthesize computes coefficients and validates") {
        const auto r = design::esc({.dither_amplitude = 0.1, .dither_omega = w(50.0), .gain = 20.0, .type = ExtremumType::Maximize}, 1e-3);
        REQUIRE(r.success);
        CHECK(r.direction == doctest::Approx(1.0));
        CHECK(r.dither_omega == doctest::Approx(2.0 * std::numbers::pi * 50.0));
        CHECK(r.hp_alpha > 0.0);
        CHECK(r.hp_alpha <= 1.0);

        // Minimize flips the integration direction.
        CHECK(design::esc({.dither_amplitude = 0.1, .dither_omega = w(50.0), .gain = 20.0, .type = ExtremumType::Minimize}, 1e-3).direction == doctest::Approx(-1.0));

        // Bad specs rejected.
        CHECK_FALSE(design::esc({.dither_amplitude = 0.0, .dither_omega = w(50.0), .gain = 20.0}, 1e-3).success);  // a <= 0
        CHECK_FALSE(design::esc({.dither_amplitude = 0.1, .dither_omega = 0.0, .gain = 20.0}, 1e-3).success);      // freq <= 0
        CHECK_FALSE(design::esc({.dither_amplitude = 0.1, .dither_omega = w(50.0), .gain = 0.0}, 1e-3).success);   // gain <= 0
        CHECK_FALSE(design::esc({.dither_amplitude = 0.1, .dither_omega = w(600.0), .gain = 20.0}, 1e-3).success); // dither >= Nyquist (500 Hz)
    }

    TEST_CASE("climbs to the maximum of an unknown quadratic map") {
        // J(u) = 5 - 0.5*(u - 2.3)^2, peak at u* = 2.3 (unknown to the ESC).
        const double u_star = 2.3;
        const auto   art = design::esc({.dither_amplitude = 0.05, .dither_omega = w(40.0), .gain = 60.0, .type = ExtremumType::Maximize, .lp_omega = w(8.0)}, 1e-3);
        REQUIRE(art.success);
        const double u_opt = seek(art, [&](double u) { return 5.0 - 0.5 * (u - u_star) * (u - u_star); }, 12000);
        CHECK(u_opt == doctest::Approx(u_star).epsilon(0.03)); // within ~3% (plus residual dither)
    }

    TEST_CASE("descends to the minimum when minimizing") {
        const double u_star = -1.5;
        const auto   art = design::esc({.dither_amplitude = 0.05, .dither_omega = w(40.0), .gain = 60.0, .type = ExtremumType::Minimize, .lp_omega = w(8.0), .u_init = 1.0}, 1e-3);
        REQUIRE(art.success);
        const double u_opt = seek(art, [&](double u) { return 2.0 + 0.5 * (u - u_star) * (u - u_star); }, 12000);
        CHECK(u_opt == doctest::Approx(u_star).epsilon(0.05));
    }

    TEST_CASE("tracks a slowly drifting optimum") {
        // u*(t) ramps; the ESC should follow it (online optimization, not one-shot).
        const auto art = design::esc({.dither_amplitude = 0.05, .dither_omega = w(40.0), .gain = 80.0, .type = ExtremumType::Maximize, .lp_omega = w(8.0)}, 1e-3);
        REQUIRE(art.success);
        ExtremumSeekingController<double> esc(art);
        double                            u = esc.input();
        double                            u_star = 0.0;
        for (int k = 0; k < 30000; ++k) {
            u_star = 1.0 + 0.5 * (k * 1e-3); // drifts from 1.0 upward
            u = esc.step(-(u - u_star) * (u - u_star));
        }
        CHECK(esc.operating_point() == doctest::Approx(u_star).epsilon(0.05)); // tracking, not lagging far
    }

    TEST_CASE("freezes the integrator on a degraded measurement") {
        const auto                        art = design::esc({.dither_amplitude = 0.05, .dither_omega = w(40.0), .gain = 60.0}, 1e-3);
        ExtremumSeekingController<double> esc(art);
        // Converge near the peak first.
        double u = esc.input();
        for (int k = 0; k < 8000; ++k) {
            u = esc.step(5.0 - 0.5 * (u - 2.0) * (u - 2.0));
        }
        const double held = esc.operating_point();
        // Now feed garbage with measurement_valid=false: û must not move.
        for (int k = 0; k < 2000; ++k) {
            u = esc.step(1e6 * std::sin(0.01 * k), false);
        }
        CHECK(esc.operating_point() == doctest::Approx(held)); // frozen exactly
    }

    TEST_CASE("MPPT wrapper clamps the operating point to its band") {
        // Maximize a "power curve" peaking outside the allowed band -> clamp holds.
        const auto art = design::esc_mppt({.dither_amplitude = 0.02, .dither_omega = w(100.0), .gain = 200.0, .u_init = 0.5}, 1e-4, design::InputLimits{.lower = {0.2}, .upper = {0.8}});
        REQUIRE(art.success);
        ExtremumSeekingController<double> esc(art);
        double                            u = esc.input();
        for (int k = 0; k < 20000; ++k) {
            u = esc.step(-(u - 1.5) * (u - 1.5)); // optimum at 1.5, above the 0.8 cap
        }
        CHECK(esc.operating_point() <= 0.8 + 1e-6); // pinned at the upper limit
        CHECK(esc.operating_point() > 0.7);         // and pushed up against it
    }

    TEST_CASE("invalid design is inert; reset restores the initial point") {
        ExtremumSeekingController<double> bad(design::esc({.dither_amplitude = 0.0, .dither_omega = w(50.0), .gain = 1.0}, 1e-3)); // invalid
        CHECK_FALSE(bad.valid());
        CHECK(bad.step(1.0) == doctest::Approx(0.0));

        ExtremumSeekingController<double> esc(design::esc({.dither_amplitude = 0.05, .dither_omega = w(40.0), .gain = 60.0, .type = ExtremumType::Maximize, .lp_omega = w(8.0), .u_init = 1.0}, 1e-3));
        for (int k = 0; k < 1000; ++k) {
            (void)esc.step(5.0 - 0.5 * (esc.input() - 2.0) * (esc.input() - 2.0));
        }
        esc.reset();
        CHECK(esc.operating_point() == doctest::Approx(1.0)); // back to u_init
        CHECK(esc.last_command() == esc.input());
    }

    TEST_CASE("float specialization seeks the maximum") {
        const auto                       art = design::esc({.dither_amplitude = 0.05f, .dither_omega = w(40.0f), .gain = 60.0f, .type = ExtremumType::Maximize, .lp_omega = w(8.0f)}, 1e-3f);
        ExtremumSeekingController<float> esc(art.as<float>());
        float                            u = esc.input();
        for (int k = 0; k < 12000; ++k) {
            u = esc.step(5.0f - 0.5f * (u - 1.7f) * (u - 1.7f));
        }
        CHECK(esc.operating_point() == doctest::Approx(1.7f).epsilon(0.05));
    }

    TEST_CASE("ESC is constexpr-evaluable") {
        constexpr double u_opt = []() consteval {
            auto                              art = design::esc({.dither_amplitude = 0.05, .dither_omega = w(40.0), .gain = 60.0, .type = ExtremumType::Maximize, .lp_omega = w(8.0)}, 1e-3);
            ExtremumSeekingController<double> esc(art);
            double                            u = esc.input();
            for (int k = 0; k < 8000; ++k) {
                u = esc.step(5.0 - 0.5 * (u - 1.0) * (u - 1.0));
            }
            return esc.operating_point();
        }();
        static_assert(u_opt > 0.9 && u_opt < 1.1, "ESC must converge at compile time");
        CHECK(u_opt == doctest::Approx(1.0).epsilon(0.05));
    }

    TEST_CASE("dither phasor keeps phase and amplitude over a long run") {
        // Frozen estimator: only the dither advances, so input() − û = a·sin(kωTs).
        constexpr double                  a = 0.1;
        const auto                        art = design::esc({.dither_amplitude = a, .dither_omega = w(50.0), .gain = 20.0, .type = ExtremumType::Maximize}, 1e-3);
        ExtremumSeekingController<double> esc(art);
        constexpr int                     N = 1000000;
        for (int k = 1; k <= N; ++k) {
            (void)esc.step(0.0, false);
        }
        const double phase = std::fmod(2.0 * std::numbers::pi * 50.0 * 1e-3 * N, 2.0 * std::numbers::pi);
        CHECK(esc.input() - esc.operating_point() == doctest::Approx(a * std::sin(phase)).epsilon(1e-9).scale(a));

        ExtremumSeekingController<float> escf(art.as<float>());
        float                            peak = 0.0f;
        for (int k = 1; k <= N; ++k) {
            (void)escf.step(0.0f, false);
            if (k > N - 100) { // last few cycles
                peak = std::max(peak, std::abs(escf.input() - escf.operating_point()));
            }
        }
        CHECK(peak == doctest::Approx(0.1f).epsilon(1e-4));
    }

    TEST_CASE("ESC: a flagged-bad objective does not leak into the estimate") {
        const auto art = design::esc({.dither_amplitude = 0.1, .dither_omega = w(5.0), .gain = 2.0}, 1e-3);
        REQUIRE(art.success);
        ExtremumSeekingController<double> a(art);
        ExtremumSeekingController<double> b(art);
        auto                              J = [](double u) { return -((u - 1.0) * (u - 1.0)); };
        double                            ua = a.input();
        double                            ub = b.input();
        for (int k = 0; k < 2000; ++k) {
            const bool fault = (k >= 500 && k < 510);
            ua = a.step(fault ? 1e6 : J(ua), !fault); // sensor spike, flagged
            ub = b.step(J(ub), !fault);               // same flag, sane value
        }
        CHECK(a.operating_point() == doctest::Approx(b.operating_point()).epsilon(1e-12));
    }

    TEST_CASE("ESC refuses NaN and negative corners, naming the cause") {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        // A NaN Ts used to pass, and the runtime's dither rotation cos(NaN) emitted NaN forever.
        const auto bad_ts = design::esc({.dither_amplitude = 0.1, .dither_omega = w(50.0), .gain = 20.0}, nan);
        CHECK_FALSE(bad_ts.success);
        CHECK(bad_ts.status == design::DesignStatus::invalid_sample_time);
        ExtremumSeekingController<double> inert(bad_ts);
        CHECK_FALSE(inert.valid());
        CHECK(design::esc({.dither_amplitude = 0.1, .dither_omega = w(50.0), .gain = 20.0}, 0.0).status == design::DesignStatus::invalid_sample_time);

        CHECK(design::esc({.dither_amplitude = nan, .dither_omega = w(50.0), .gain = 20.0}, 1e-3).status == design::DesignStatus::invalid_parameter);
        CHECK(design::esc({.dither_amplitude = 0.1, .dither_omega = w(50.0), .gain = nan}, 1e-3).status == design::DesignStatus::invalid_parameter);
        // A negative corner used to mean "no low-pass" or "the default high-pass".
        CHECK(design::esc({.dither_amplitude = 0.1, .dither_omega = w(50.0), .gain = 20.0, .lp_omega = -w(8.0)}, 1e-3).status == design::DesignStatus::invalid_parameter);
        CHECK(design::esc({.dither_amplitude = 0.1, .dither_omega = w(50.0), .gain = 20.0, .hp_omega = -w(8.0)}, 1e-3).status == design::DesignStatus::invalid_parameter);
        CHECK(design::esc({.dither_amplitude = 0.1, .dither_omega = w(600.0), .gain = 20.0}, 1e-3).status == design::DesignStatus::invalid_parameter);
    }

    TEST_CASE("InputLimits is a box or the design is refused; it starts û inside and moves at run time") {
        //! lower ≥ upper, or a NaN end, switched the clamp off without a word.
        const double nan = std::numeric_limits<double>::quiet_NaN();
        CHECK(design::esc({.dither_amplitude = 0.1, .dither_omega = w(50.0), .gain = 20.0}, 1e-3, design::InputLimits{.lower = {1.0}, .upper = {-1.0}}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::esc({.dither_amplitude = 0.1, .dither_omega = w(50.0), .gain = 20.0}, 1e-3, design::InputLimits{.lower = {nan}, .upper = {1.0}}).status == design::DesignStatus::invalid_parameter);

        const auto art = design::esc_mppt({.dither_amplitude = 0.02, .dither_omega = w(100.0), .gain = 200.0, .u_init = 2.0}, 1e-4, design::InputLimits{.lower = {0.2}, .upper = {0.8}});
        REQUIRE(art.success);
        ExtremumSeekingController<double> esc(art);
        CHECK(esc.operating_point() == 0.8); // u_init outside the box starts on it
        auto lim = esc.limits();
        lim.input.upper(0) = 0.5; // a rail that moves
        REQUIRE(esc.set_limits(lim));
        CHECK(esc.operating_point() == 0.5); // û moves inside at once
        double u = esc.input();
        for (int k = 0; k < 2000; ++k) {
            u = esc.step(-(u - 1.5) * (u - 1.5));
            CHECK(u <= 0.5);
        }
        CHECK(esc.operating_point() == 0.5);
        esc.reset();
        CHECK(esc.operating_point() == 0.5); // u_init into the box as it stands
    }

    TEST_CASE("the operating point's own limit keeps the dither unclipped; set_limits refuses a limit that is not a box") {
        const double a = 0.02;
        const auto   art = design::esc_mppt({.dither_amplitude = a, .dither_omega = w(100.0), .gain = 200.0, .u_init = 0.5}, 1e-4, design::InputLimits{.lower = {0.2}, .upper = {0.8}});
        REQUIRE(art.success);
        ExtremumSeekingController<double> esc(art);
        auto                              lim = esc.limits();
        lim.operating_point = {.lower = 0.2 + a, .upper = 0.8 - a}; // a dither's amplitude inside
        REQUIRE(esc.set_limits(lim));
        double u = esc.input();
        double lo = 1.0;
        double hi = 0.0;
        for (int k = 0; k < 20000; ++k) {
            u = esc.step(-(u - 1.5) * (u - 1.5)); // optimum above the band
            if (k > 15000) {
                lo = std::min(lo, u);
                hi = std::max(hi, u);
            }
        }
        CHECK(esc.operating_point() == doctest::Approx(0.8 - a));
        CHECK(hi - lo == doctest::Approx(2.0 * a).epsilon(1e-3)); // the whole dither, never clipped

        const double nan = runtime_nan();
        for (int end = 0; end < 4; ++end) {
            auto    bad = esc.limits();
            double* ends[] = {&bad.input.lower(0), &bad.input.upper(0), &bad.operating_point.lower(0), &bad.operating_point.upper(0)};
            *ends[end] = nan;
            CHECK_FALSE(esc.set_limits(bad));
            CHECK(esc.limits().operating_point.upper(0) == 0.8 - a);
            CHECK(esc.limits().input.upper(0) == 0.8);
        }
        auto inverted = esc.limits();
        inverted.operating_point = {.lower = 0.7, .upper = 0.3};
        CHECK_FALSE(esc.set_limits(inverted));
        CHECK(esc.operating_point() == doctest::Approx(0.8 - a)); // nothing changed
    }

    TEST_CASE("track(u) starts the operating point from the input applied, with no stale gradient") {
        //! An operator drove the plant at 1.2 while the law was parked; without a handover it resumed
        //! at its old û and its old gradient.
        const auto art = design::esc({.dither_amplitude = 0.05, .dither_omega = w(40.0), .gain = 60.0, .lp_omega = w(8.0)}, 1e-3);
        REQUIRE(art.success);
        ExtremumSeekingController<double> esc(art);
        double                            u = esc.input();
        for (int k = 0; k < 3000; ++k) {
            u = esc.step(5.0 - 0.5 * (u - 2.0) * (u - 2.0));
        }
        REQUIRE(esc.gradient() != 0.0);
        for (int k = 0; k < 100; ++k) {
            esc.track(1.2);
        }
        CHECK(esc.operating_point() == 1.2);
        CHECK(esc.last_command() == 1.2);
        CHECK(esc.gradient() == 0.0);
        const double first = esc.step(5.0 - 0.5 * (1.2 - 2.0) * (1.2 - 2.0)); // the objective at the input applied
        CHECK(esc.operating_point() == 1.2);                                  // the high-pass seeds: no gradient yet
        CHECK(damp::abs(first - 1.2) <= 0.05 + 1e-12);                        // u plus the dither
        u = first;
        for (int k = 0; k < 12000; ++k) {
            u = esc.step(5.0 - 0.5 * (u - 2.0) * (u - 2.0));
        }
        CHECK(esc.operating_point() == doctest::Approx(2.0).epsilon(0.03)); // and seeks from there
    }

    TEST_CASE("back_calculate keeps the operating point where a clamp after the law left the input") {
        //! With no way to report a clamp downstream, û climbed past it until the whole dither sat
        //! above the clamp and the objective stopped answering.
        const auto art = design::esc_mppt({.dither_amplitude = 0.02, .dither_omega = w(100.0), .gain = 200.0, .u_init = 0.5}, 1e-4);
        REQUIRE(art.success);
        ExtremumSeekingController<double> told(art);
        ExtremumSeekingController<double> untold(art);
        const auto                        objective = [](double u) { return -(u - 1.5) * (u - 1.5); }; // optimum above the clamp
        const double                      rail = 0.8;
        double                            ua = told.input();
        double                            ub = untold.input();
        for (int k = 0; k < 20000; ++k) {
            ua = damp::min(told.step(objective(ua)), rail);
            told.back_calculate(ua);
            ub = damp::min(untold.step(objective(ub)), rail);
        }
        CHECK(untold.operating_point() > rail);
        CHECK(told.operating_point() <= rail);
        CHECK(told.operating_point() >= rail - (2.0 * 0.02) - 1e-9); // within the dither's swing of it

        //! An input applied as commanded changes nothing: the two are one run.
        ExtremumSeekingController<double> plain(art);
        ExtremumSeekingController<double> echoed(art);
        double                            up = plain.input();
        double                            ue = echoed.input();
        for (int k = 0; k < 2000; ++k) {
            up = plain.step(-(up - 0.6) * (up - 0.6));
            ue = echoed.step(-(ue - 0.6) * (ue - 0.6));
            echoed.back_calculate(ue);
        }
        CHECK(echoed.operating_point() == plain.operating_point());
    }

    TEST_CASE("Law of an extremum-seeking design runs the same law") {
        static constexpr auto             kMppt = design::esc_mppt({.dither_amplitude = 0.02, .dither_omega = w(100.0), .gain = 200.0, .u_init = 0.5}, 1e-4, design::InputLimits{.lower = {0.2}, .upper = {0.8}});
        Law<kMppt, double>                law;
        ExtremumSeekingController<double> held{kMppt};
        double                            ul = law.input();
        double                            uh = held.input();
        for (int k = 0; k < 500; ++k) {
            ul = law.step(-(ul - 0.6) * (ul - 0.6));
            uh = held.step(-(uh - 0.6) * (uh - 0.6));
        }
        CHECK(ul == uh);
        CHECK(law.limits().input.upper(0) == 0.8);
    }
}
