// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>

#include "damp/estimation/commissioning_frontends.hpp"
#include "damp/estimation/experiment_safety.hpp"
#include "damp/estimation/frequency_response.hpp"
#include "damp/estimation/successive_compensator.hpp"
#include "damp/math/math.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

TEST_SUITE("Experiment safety") {
    TEST_CASE("protect clamps and slew-limits") {
        ExperimentSafetyConfig<float> cfg{
            .command_bounds = Bounds<1, float>{-1.0f, 1.0f},
            .max_slew_rate = 10.0f, // units/s
            .timeout_s = 0.5f,
        };
        ExperimentSafety<2, float> saf{cfg};
        damp::array<float, 2>      state{1.0f, 2.0f};
        saf.arm(state);
        CHECK(saf.armed());

        const float dt = 0.01f;
        // Request +10 but slew allows only 0.1 per step
        float u = saf.protect(10.0f, dt);
        CHECK(u == doctest::Approx(0.1f).epsilon(1e-5));
        u = saf.protect(10.0f, dt);
        CHECK(u == doctest::Approx(0.2f).epsilon(1e-5));
        // Clamp
        for (int i = 0; i < 200; ++i) {
            u = saf.protect(10.0f, dt);
        }
        CHECK(u == doctest::Approx(1.0f).epsilon(1e-5));

        // Timeout
        bool expired = false;
        for (int i = 0; i < 60; ++i) {
            expired = saf.tick_timeout(dt);
        }
        CHECK(expired);
        CHECK(saf.timed_out());

        damp::array<float, 2> dirty{9.0f, 9.0f};
        saf.working() = dirty;
        saf.rollback(dirty);
        CHECK(dirty[0] == doctest::Approx(1.0f));
        CHECK(dirty[1] == doctest::Approx(2.0f));
        CHECK_FALSE(saf.armed());
    }

    TEST_CASE("commit ends experiment without restoring") {
        ExperimentSafety<1, double> saf{ExperimentSafetyConfig<double>{
            .command_bounds = Bounds<1, double>{-2.0, 2.0},
        }};
        damp::array<double, 1>      s{3.0};
        saf.arm(s);
        saf.working()[0] = 7.0;
        saf.commit();
        CHECK(saf.committed());
        CHECK_FALSE(saf.armed());
        CHECK(saf.working()[0] == doctest::Approx(7.0));
    }
}

TEST_SUITE("Commissioning front-ends") {
    TEST_CASE("resonance_autotune_from_frf places notch on peak") {
        damp::array<FrfPoint<double>, 5> table{};
        table[0] = {10.0, 0.5, 0.0, 1.0, true};
        table[1] = {20.0, 1.0, 0.0, 1.0, true};
        table[2] = {30.0, 3.0, 0.0, 1.0, true};
        table[3] = {40.0, 1.0, 0.0, 1.0, true};
        table[4] = {50.0, 0.5, 0.0, 1.0, true};
        ModeExtractorConfig<double> mcfg{};
        mcfg.require_valid_zeta = false;
        SuccessiveCompensatorConfig<double> cfg{.Ts = 0.001, .notch_Q = 8.0};
        const auto                          res = design::resonance_autotune_from_frf<5, 2>(table, cfg, mcfg);
        REQUIRE(res.success);
        CHECK(res.count >= 1);
        CHECK(res.features[0].frequency_hz == doctest::Approx(30.0));
        CHECK(res.kinds[0] == CompensatorKind::Notch);
    }

    TEST_CASE("frf_pid_autotune returns finite gains") {
        constexpr std::size_t             Nf = 16;
        damp::array<FrfPoint<double>, Nf> table{};
        // Integrator-like plant G = 5 / (j w)
        for (std::size_t i = 0; i < Nf; ++i) {
            const double f = 0.5 + static_cast<double>(i) * 1.0;
            const double w = 2.0 * damp::numbers::pi_v<double> * f;
            table[i] = {f, 5.0 / w, -damp::numbers::pi_v<double> / 2.0, 1.0, true};
        }
        const auto art = design::frf_pid_autotune(
            table,
            design::FrfPidAutotuneConfig<double>{
                .target_crossover_hz = 2.0,
                .target_phase_margin_deg = 60.0,
                .Ts = 0.01,
                .type = design::PIDType::PI,
            }
        );
        CHECK(art.success);
        CHECK(damp::isfinite(art.pid.Kp));
        CHECK(damp::isfinite(art.pid.Ki));
        CHECK(art.pid.Kp > 0.0);
    }
}

TEST_SUITE("Successive compensators FRF path") {
    TEST_CASE("valley maps to bandpass section") {
        SuccessiveCompensatorCommissioner<2, double> cm{SuccessiveCompensatorConfig<double>{
            .Ts = 0.001,
            .valley_map = CompensatorKind::Bandpass,
        }};
        ResonantMode<double>                         v{};
        v.kind = FrfFeatureKind::Valley;
        v.frequency_hz = 50.0;
        v.omega_n = 2.0 * damp::numbers::pi_v<double> * 50.0;
        v.zeta = 0.05;
        v.peak_gain = 0.1;
        v.valid = true;
        REQUIRE(cm.assign_feature(v));
        CHECK(cm.kind(0) == CompensatorKind::Bandpass);
        CHECK(cm.count() == 1);
    }
}
