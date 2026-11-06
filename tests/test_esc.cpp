// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

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
        u = esc.control(J(u));
    }
    return esc.operating_point();
}
} // namespace

TEST_SUITE("esc") {
    TEST_CASE("synthesize computes coefficients and validates") {
        const auto r = design::esc({.dither_amplitude = 0.1, .dither_freq_hz = 50.0, .gain = 20.0, .type = ExtremumType::Maximize}, 1e-3);
        REQUIRE(r.success);
        CHECK(r.direction == doctest::Approx(1.0));
        CHECK(r.dither_omega == doctest::Approx(2.0 * std::numbers::pi * 50.0));
        CHECK(r.hp_alpha > 0.0);
        CHECK(r.hp_alpha <= 1.0);

        // Minimize flips the integration direction.
        CHECK(design::esc({.dither_amplitude = 0.1, .dither_freq_hz = 50.0, .gain = 20.0, .type = ExtremumType::Minimize}, 1e-3).direction == doctest::Approx(-1.0));

        // Bad specs rejected.
        CHECK_FALSE(design::esc({.dither_amplitude = 0.0, .dither_freq_hz = 50.0, .gain = 20.0}, 1e-3).success);  // a <= 0
        CHECK_FALSE(design::esc({.dither_amplitude = 0.1, .dither_freq_hz = 0.0, .gain = 20.0}, 1e-3).success);   // freq <= 0
        CHECK_FALSE(design::esc({.dither_amplitude = 0.1, .dither_freq_hz = 50.0, .gain = 0.0}, 1e-3).success);   // gain <= 0
        CHECK_FALSE(design::esc({.dither_amplitude = 0.1, .dither_freq_hz = 600.0, .gain = 20.0}, 1e-3).success); // dither >= Nyquist (500 Hz)
    }

    TEST_CASE("climbs to the maximum of an unknown quadratic map") {
        // J(u) = 5 - 0.5*(u - 2.3)^2, peak at u* = 2.3 (unknown to the ESC).
        const double u_star = 2.3;
        const auto   art = design::esc({.dither_amplitude = 0.05, .dither_freq_hz = 40.0, .gain = 60.0, .type = ExtremumType::Maximize, .lp_cutoff_hz = 8.0}, 1e-3);
        REQUIRE(art.success);
        const double u_opt = seek(art, [&](double u) { return 5.0 - 0.5 * (u - u_star) * (u - u_star); }, 12000);
        CHECK(u_opt == doctest::Approx(u_star).epsilon(0.03)); // within ~3% (plus residual dither)
    }

    TEST_CASE("descends to the minimum when minimizing") {
        const double u_star = -1.5;
        const auto   art = design::esc({.dither_amplitude = 0.05, .dither_freq_hz = 40.0, .gain = 60.0, .type = ExtremumType::Minimize, .lp_cutoff_hz = 8.0, .u_init = 1.0}, 1e-3);
        REQUIRE(art.success);
        const double u_opt = seek(art, [&](double u) { return 2.0 + 0.5 * (u - u_star) * (u - u_star); }, 12000);
        CHECK(u_opt == doctest::Approx(u_star).epsilon(0.05));
    }

    TEST_CASE("tracks a slowly drifting optimum") {
        // u*(t) ramps; the ESC should follow it (online optimization, not one-shot).
        const auto art = design::esc({.dither_amplitude = 0.05, .dither_freq_hz = 40.0, .gain = 80.0, .type = ExtremumType::Maximize, .lp_cutoff_hz = 8.0}, 1e-3);
        REQUIRE(art.success);
        ExtremumSeekingController<double> esc(art);
        double                            u = esc.input();
        double                            u_star = 0.0;
        for (int k = 0; k < 30000; ++k) {
            u_star = 1.0 + 0.5 * (k * 1e-3); // drifts from 1.0 upward
            u = esc.control(-(u - u_star) * (u - u_star));
        }
        CHECK(esc.operating_point() == doctest::Approx(u_star).epsilon(0.05)); // tracking, not lagging far
    }

    TEST_CASE("freezes the integrator on a degraded measurement") {
        const auto                        art = design::esc({.dither_amplitude = 0.05, .dither_freq_hz = 40.0, .gain = 60.0}, 1e-3);
        ExtremumSeekingController<double> esc(art);
        // Converge near the peak first.
        double u = esc.input();
        for (int k = 0; k < 8000; ++k) {
            u = esc.control(5.0 - 0.5 * (u - 2.0) * (u - 2.0));
        }
        const double held = esc.operating_point();
        // Now feed garbage with measurement_valid=false: û must not move.
        for (int k = 0; k < 2000; ++k) {
            u = esc.control(1e6 * std::sin(0.01 * k), false);
        }
        CHECK(esc.operating_point() == doctest::Approx(held)); // frozen exactly
    }

    TEST_CASE("MPPT wrapper clamps the operating point to its band") {
        // Maximize a "power curve" peaking outside the allowed band -> clamp holds.
        const auto art = design::esc_mppt({.dither_amplitude = 0.02, .dither_freq_hz = 100.0, .gain = 200.0, .u_init = 0.5, .u_lim = {
                                                                                                                                .lower = 0.2,
                                                                                                                                .upper = 0.8,
                                                                                                                            }},
                                          1e-4);
        REQUIRE(art.success);
        ExtremumSeekingController<double> esc(art);
        double                            u = esc.input();
        for (int k = 0; k < 20000; ++k) {
            u = esc.control(-(u - 1.5) * (u - 1.5)); // optimum at 1.5, above the 0.8 cap
        }
        CHECK(esc.operating_point() <= 0.8 + 1e-6); // pinned at the upper limit
        CHECK(esc.operating_point() > 0.7);         // and pushed up against it
    }

    TEST_CASE("invalid design is inert; reset restores the initial point") {
        ExtremumSeekingController<double> bad(design::esc({.dither_amplitude = 0.0, .dither_freq_hz = 50.0, .gain = 1.0}, 1e-3)); // invalid
        CHECK_FALSE(bad.valid());
        CHECK(bad.control(1.0) == doctest::Approx(0.0));

        ExtremumSeekingController<double> esc(design::esc({.dither_amplitude = 0.05, .dither_freq_hz = 40.0, .gain = 60.0, .type = ExtremumType::Maximize, .lp_cutoff_hz = 8.0, .u_init = 1.0}, 1e-3));
        for (int k = 0; k < 1000; ++k) {
            (void)esc.control(5.0 - 0.5 * (esc.input() - 2.0) * (esc.input() - 2.0));
        }
        esc.reset();
        CHECK(esc.operating_point() == doctest::Approx(1.0)); // back to u_init
    }

    TEST_CASE("float specialization seeks the maximum") {
        const auto                       art = design::esc({.dither_amplitude = 0.05f, .dither_freq_hz = 40.0f, .gain = 60.0f, .type = ExtremumType::Maximize, .lp_cutoff_hz = 8.0f}, 1e-3f);
        ExtremumSeekingController<float> esc(art.as<float>());
        float                            u = esc.input();
        for (int k = 0; k < 12000; ++k) {
            u = esc.control(5.0f - 0.5f * (u - 1.7f) * (u - 1.7f));
        }
        CHECK(esc.operating_point() == doctest::Approx(1.7f).epsilon(0.05));
    }

    TEST_CASE("ESC is constexpr-evaluable") {
        constexpr double u_opt = []() consteval {
            auto                              art = design::esc({.dither_amplitude = 0.05, .dither_freq_hz = 40.0, .gain = 60.0, .type = ExtremumType::Maximize, .lp_cutoff_hz = 8.0}, 1e-3);
            ExtremumSeekingController<double> esc(art);
            double                            u = esc.input();
            for (int k = 0; k < 8000; ++k) {
                u = esc.control(5.0 - 0.5 * (u - 1.0) * (u - 1.0));
            }
            return esc.operating_point();
        }();
        static_assert(u_opt > 0.9 && u_opt < 1.1, "ESC must converge at compile time");
        CHECK(u_opt == doctest::Approx(1.0).epsilon(0.05));
    }

    TEST_CASE("dither phasor keeps phase and amplitude over a long run") {
        // Frozen estimator: only the dither advances, so input() − û = a·sin(kωTs).
        constexpr double                  a = 0.1;
        const auto                        art = design::esc({.dither_amplitude = a, .dither_freq_hz = 50.0, .gain = 20.0, .type = ExtremumType::Maximize}, 1e-3);
        ExtremumSeekingController<double> esc(art);
        constexpr int                     N = 1000000;
        for (int k = 1; k <= N; ++k) {
            (void)esc.control(0.0, false);
        }
        const double phase = std::fmod(2.0 * std::numbers::pi * 50.0 * 1e-3 * N, 2.0 * std::numbers::pi);
        CHECK(esc.input() - esc.operating_point() == doctest::Approx(a * std::sin(phase)).epsilon(1e-9).scale(a));

        ExtremumSeekingController<float> escf(art.as<float>());
        float                            peak = 0.0f;
        for (int k = 1; k <= N; ++k) {
            (void)escf.control(0.0f, false);
            if (k > N - 100) { // last few cycles
                peak = std::max(peak, std::abs(escf.input() - escf.operating_point()));
            }
        }
        CHECK(peak == doctest::Approx(0.1f).epsilon(1e-4));
    }

    TEST_CASE("ESC: a flagged-bad objective does not leak into the estimate") {
        const auto art = design::esc({.dither_amplitude = 0.1, .dither_freq_hz = 5.0, .gain = 2.0}, 1e-3);
        REQUIRE(art.success);
        ExtremumSeekingController<double> a(art);
        ExtremumSeekingController<double> b(art);
        auto                              J = [](double u) { return -((u - 1.0) * (u - 1.0)); };
        double                            ua = a.input();
        double                            ub = b.input();
        for (int k = 0; k < 2000; ++k) {
            const bool fault = (k >= 500 && k < 510);
            ua = a.control(fault ? 1e6 : J(ua), !fault); // sensor spike, flagged
            ub = b.control(J(ub), !fault);               // same flag, sane value
        }
        CHECK(a.operating_point() == doctest::Approx(b.operating_point()).epsilon(1e-12));
    }

    TEST_CASE("ESC refuses NaN and negative corners, naming the cause") {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        // A NaN Ts used to pass, and the runtime's dither rotation cos(NaN) emitted NaN forever.
        const auto bad_ts = design::esc({.dither_amplitude = 0.1, .dither_freq_hz = 50.0, .gain = 20.0}, nan);
        CHECK_FALSE(bad_ts.success);
        CHECK(bad_ts.status == design::DesignStatus::invalid_sample_time);
        ExtremumSeekingController<double> inert(bad_ts);
        CHECK_FALSE(inert.valid());

        CHECK(design::esc({.dither_amplitude = nan, .dither_freq_hz = 50.0, .gain = 20.0}, 1e-3).status == design::DesignStatus::invalid_parameter);
        CHECK(design::esc({.dither_amplitude = 0.1, .dither_freq_hz = 50.0, .gain = nan}, 1e-3).status == design::DesignStatus::invalid_parameter);
        // A negative corner used to mean "no low-pass" or "the default high-pass".
        CHECK(design::esc({.dither_amplitude = 0.1, .dither_freq_hz = 50.0, .gain = 20.0, .lp_cutoff_hz = -8.0}, 1e-3).status == design::DesignStatus::invalid_parameter);
        CHECK(design::esc({.dither_amplitude = 0.1, .dither_freq_hz = 50.0, .gain = 20.0, .hp_cutoff_hz = -8.0}, 1e-3).status == design::DesignStatus::invalid_parameter);
        CHECK(design::esc({.dither_amplitude = 0.1, .dither_freq_hz = 600.0, .gain = 20.0}, 1e-3).status == design::DesignStatus::invalid_parameter);
    }
}
