// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <numbers>

#include "damp/controllers/esc.hpp"

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
    return esc.estimate();
}
} // namespace

TEST_SUITE("esc") {
    TEST_CASE("synthesize computes coefficients and validates") {
        const auto r = design::esc(0.1, 50.0, 20.0, 1e-3, ExtremumType::Maximize);
        REQUIRE(r.success);
        CHECK(r.config.direction == doctest::Approx(1.0));
        CHECK(r.config.dither_omega == doctest::Approx(2.0 * std::numbers::pi * 50.0));
        CHECK(r.config.hp_alpha > 0.0);
        CHECK(r.config.hp_alpha <= 1.0);

        // Minimize flips the integration direction.
        CHECK(design::esc(0.1, 50.0, 20.0, 1e-3, ExtremumType::Minimize).config.direction == doctest::Approx(-1.0));

        // Bad specs rejected.
        CHECK_FALSE(design::esc(0.0, 50.0, 20.0, 1e-3).success);  // a <= 0
        CHECK_FALSE(design::esc(0.1, 0.0, 20.0, 1e-3).success);   // freq <= 0
        CHECK_FALSE(design::esc(0.1, 50.0, 0.0, 1e-3).success);   // gain <= 0
        CHECK_FALSE(design::esc(0.1, 600.0, 20.0, 1e-3).success); // dither >= Nyquist (500 Hz)
    }

    TEST_CASE("climbs to the maximum of an unknown quadratic map") {
        // J(u) = 5 - 0.5*(u - 2.3)^2, peak at u* = 2.3 (unknown to the ESC).
        const double u_star = 2.3;
        const auto   art = design::esc(0.05, 40.0, 60.0, 1e-3, ExtremumType::Maximize, 0.0, 8.0, 0.0);
        REQUIRE(art.success);
        const double u_opt = seek(art, [&](double u) { return 5.0 - 0.5 * (u - u_star) * (u - u_star); }, 12000);
        CHECK(u_opt == doctest::Approx(u_star).epsilon(0.03)); // within ~3% (plus residual dither)
    }

    TEST_CASE("descends to the minimum when minimizing") {
        const double u_star = -1.5;
        const auto   art = design::esc(0.05, 40.0, 60.0, 1e-3, ExtremumType::Minimize, 0.0, 8.0, 1.0);
        REQUIRE(art.success);
        const double u_opt = seek(art, [&](double u) { return 2.0 + 0.5 * (u - u_star) * (u - u_star); }, 12000);
        CHECK(u_opt == doctest::Approx(u_star).epsilon(0.05));
    }

    TEST_CASE("tracks a slowly drifting optimum") {
        // u*(t) ramps; the ESC should follow it (online optimization, not one-shot).
        const auto art = design::esc(0.05, 40.0, 80.0, 1e-3, ExtremumType::Maximize, 0.0, 8.0, 0.0);
        REQUIRE(art.success);
        ExtremumSeekingController<double> esc(art);
        double                            u = esc.input();
        double                            u_star = 0.0;
        for (int k = 0; k < 30000; ++k) {
            u_star = 1.0 + 0.5 * (k * 1e-3); // drifts from 1.0 upward
            u = esc.step(-(u - u_star) * (u - u_star));
        }
        CHECK(esc.estimate() == doctest::Approx(u_star).epsilon(0.05)); // tracking, not lagging far
    }

    TEST_CASE("freezes the integrator on a degraded measurement") {
        const auto                        art = design::esc(0.05, 40.0, 60.0, 1e-3);
        ExtremumSeekingController<double> esc(art);
        // Converge near the peak first.
        double u = esc.input();
        for (int k = 0; k < 8000; ++k) {
            u = esc.step(5.0 - 0.5 * (u - 2.0) * (u - 2.0));
        }
        const double held = esc.estimate();
        // Now feed garbage with measurement_valid=false: û must not move.
        for (int k = 0; k < 2000; ++k) {
            u = esc.step(1e6 * std::sin(0.01 * k), false);
        }
        CHECK(esc.estimate() == doctest::Approx(held)); // frozen exactly
    }

    TEST_CASE("MPPT wrapper clamps the operating point to its band") {
        // Maximize a "power curve" peaking outside the allowed band -> clamp holds.
        const auto art = design::esc_mppt(0.02, 100.0, 200.0, 1e-4, 0.5, 0.2, 0.8);
        REQUIRE(art.success);
        ExtremumSeekingController<double> esc(art);
        double                            u = esc.input();
        for (int k = 0; k < 20000; ++k) {
            u = esc.step(-(u - 1.5) * (u - 1.5)); // optimum at 1.5, above the 0.8 cap
        }
        CHECK(esc.estimate() <= 0.8 + 1e-6); // pinned at the upper limit
        CHECK(esc.estimate() > 0.7);         // and pushed up against it
    }

    TEST_CASE("invalid design is inert; reset restores the initial point") {
        ExtremumSeekingController<double> bad(design::esc(0.0, 50.0, 1.0, 1e-3)); // invalid
        CHECK_FALSE(bad.valid());
        CHECK(bad.step(1.0) == doctest::Approx(0.0));

        ExtremumSeekingController<double> esc(design::esc(0.05, 40.0, 60.0, 1e-3, ExtremumType::Maximize, 0.0, 8.0, 1.0));
        for (int k = 0; k < 1000; ++k) {
            (void)esc.step(5.0 - 0.5 * (esc.input() - 2.0) * (esc.input() - 2.0));
        }
        esc.reset();
        CHECK(esc.estimate() == doctest::Approx(1.0)); // back to u_init
    }

    TEST_CASE("float specialization seeks the maximum") {
        const auto                       art = design::esc(0.05f, 40.0f, 60.0f, 1e-3f, ExtremumType::Maximize, 0.0f, 8.0f, 0.0f);
        ExtremumSeekingController<float> esc(art.as<float>());
        float                            u = esc.input();
        for (int k = 0; k < 12000; ++k) {
            u = esc.step(5.0f - 0.5f * (u - 1.7f) * (u - 1.7f));
        }
        CHECK(esc.estimate() == doctest::Approx(1.7f).epsilon(0.05));
    }

    TEST_CASE("ESC is constexpr-evaluable") {
        constexpr double u_opt = []() consteval {
            auto                              art = design::esc(0.05, 40.0, 60.0, 1e-3, ExtremumType::Maximize, 0.0, 8.0, 0.0);
            ExtremumSeekingController<double> esc(art);
            double                            u = esc.input();
            for (int k = 0; k < 8000; ++k) {
                u = esc.step(5.0 - 0.5 * (u - 1.0) * (u - 1.0));
            }
            return esc.estimate();
        }();
        static_assert(u_opt > 0.9 && u_opt < 1.1, "ESC must converge at compile time");
        CHECK(u_opt == doctest::Approx(1.0).epsilon(0.05));
    }
}
