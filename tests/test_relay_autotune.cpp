// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>
#include <numbers>

#include "damp/backend.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/design/pid_design.hpp"
#include "damp/estimation/relay_autotune.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

/// Simulate a discrete LTI plant under closed-loop relay autotuning.
template<size_t NX, typename T>
RelayAutotuneOutput<T> run_autotune(
    const StateSpace<NX, 1, 1, T, NX, 1>& plant_d,
    RelayAutotuner<T>&                    tuner,
    int                                   max_steps
) {
    ColVec<NX, T>          x{};
    RelayAutotuneOutput<T> out{};
    for (int i = 0; i < max_steps; ++i) {
        const T y = (plant_d.C * x)[0];
        out = tuner.step(y);
        const ColVec<1, T> u{out.u};
        x = plant_d.A * x + plant_d.B * u;
        if (out.status == RelayAutotuneStatus::Done) {
            break;
        }
        if (out.status == RelayAutotuneStatus::Failed) {
            break;
        }
    }
    return out;
}

} // namespace

TEST_SUITE("Relay Autotune Design") {
    TEST_CASE("config validation accepts sensible defaults") {
        constexpr design::RelayAutotuneConfig<double> cfg{};
        static_assert(cfg.valid());
        constexpr auto res = design::relay_autotune(cfg);
        static_assert(res.success);
    }

    TEST_CASE("config validation rejects non-positive amplitude") {
        constexpr design::RelayAutotuneConfig<double> cfg{.amplitude = 0.0};
        static_assert(!cfg.valid());
        constexpr auto res = design::relay_autotune(cfg);
        static_assert(!res.success);
    }

    TEST_CASE("config validation rejects negative hysteresis") {
        constexpr design::RelayAutotuneConfig<double> cfg{.hysteresis = -0.1};
        static_assert(!cfg.valid());
    }

    TEST_CASE("config validation rejects measure_cycles < 2") {
        constexpr design::RelayAutotuneConfig<double> cfg{.measure_cycles = 1};
        static_assert(!cfg.valid());
    }

    TEST_CASE("config validation rejects hysteresis >= amplitude") {
        // ε ≥ d cannot sustain a limit cycle (Kᵤ = 4d/(π√(a²−ε²)) ill-defined).
        constexpr design::RelayAutotuneConfig<double> cfg{.amplitude = 1.0, .hysteresis = 1.0};
        static_assert(!cfg.valid());
    }

    TEST_CASE(".as<float>() round-trips config values") {
        constexpr design::RelayAutotuneConfig<double> cfg{
            .amplitude = 2.5,
            .hysteresis = 0.05,
            .setpoint = 1.0,
            .u_bias = 0.25,
            .warmup_cycles = 3,
            .measure_cycles = 8,
            .period_tolerance = 0.1,
            .max_duration = 20.0,
        };
        constexpr auto res = design::relay_autotune(cfg);
        constexpr auto res_f = res.template as<float>();
        static_assert(res_f.success);
        static_assert(res_f.config.amplitude == 2.5f);
        static_assert(res_f.config.hysteresis == 0.05f);
        static_assert(res_f.config.warmup_cycles == 3);
        static_assert(res_f.config.measure_cycles == 8);
    }
} // TEST_SUITE

TEST_SUITE("Relay Autotune Runtime") {
    TEST_CASE("invalid design starts in Failed state and stays there") {
        constexpr design::RelayAutotuneConfig<float> cfg{.amplitude = 0.0f};
        constexpr auto                               res = design::relay_autotune(cfg);

        RelayAutotuner<float> tuner(res, 0.01f);
        CHECK(tuner.status() == RelayAutotuneStatus::Failed);

        const auto out = tuner.step(0.0f);
        CHECK(out.status == RelayAutotuneStatus::Failed);
    }

    TEST_CASE("non-positive Ts marks the design Failed") {
        constexpr design::RelayAutotuneConfig<float> cfg{};
        constexpr auto                               res = design::relay_autotune(cfg);

        RelayAutotuner<float> tuner(res, 0.0f);
        CHECK(tuner.status() == RelayAutotuneStatus::Failed);
    }

    TEST_CASE("idle transitions to warmup on first step") {
        constexpr design::RelayAutotuneConfig<float> cfg{
            .amplitude = 1.0f,
            .warmup_cycles = 2,
        };
        constexpr auto        res = design::relay_autotune(cfg);
        RelayAutotuner<float> tuner(res, 0.01f);

        REQUIRE(tuner.status() == RelayAutotuneStatus::Idle);
        const auto out = tuner.step(0.0f);
        CHECK(out.status == RelayAutotuneStatus::Warmup);
    }

    TEST_CASE("warmup_cycles = 0 skips straight to measuring") {
        constexpr design::RelayAutotuneConfig<float> cfg{.warmup_cycles = 0};
        constexpr auto                               res = design::relay_autotune(cfg);
        RelayAutotuner<float>                        tuner(res, 0.01f);

        const auto out = tuner.step(0.0f);
        CHECK(out.status == RelayAutotuneStatus::Measuring);
    }

    TEST_CASE("timeout marks the experiment Failed") {
        // Plant that never crosses the hysteresis band: constant y = 0,
        // setpoint = 0 with sizable hysteresis → no zero crossings.
        constexpr design::RelayAutotuneConfig<float> cfg{
            .amplitude = 1.0f,
            .hysteresis = 10.0f,   // huge band → never triggers
            .max_duration = 0.05f, // 5 ticks at Ts=0.01
        };
        constexpr auto        res = design::relay_autotune(cfg);
        RelayAutotuner<float> tuner(res, 0.01f);

        RelayAutotuneOutput<float> out{};
        for (int i = 0; i < 100; ++i) {
            out = tuner.step(0.0f);
            if (out.status == RelayAutotuneStatus::Failed) {
                break;
            }
        }
        CHECK(out.status == RelayAutotuneStatus::Failed);
    }

    TEST_CASE("recovers Ku and Tu for a three-pole plant 1/(s+1)^3") {
        // Continuous plant: 1/(s+1)^3 in controllable-canonical form.
        //   x_dot = A x + B u,  y = C x
        // Phase crosses -180° at omega_u = tan(60°) = sqrt(3), so:
        //   Tu = 2*pi/sqrt(3) ≈ 3.628 s
        //   |G(j*omega_u)| = 1/8  →  Ku = 8
        constexpr double Ts = 0.05;

        const StateSpace<3, 1, 1, double, 3, 1> plant{
            .A = {{0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {-1.0, -3.0, -3.0}},
            .B = {{0.0}, {0.0}, {1.0}},
            .C = {{1.0, 0.0, 0.0}},
            .Ts = 0.0,
        };
        const auto plant_d = *discretize(plant, Ts, DiscretizationMethod::ZOH);

        constexpr design::RelayAutotuneConfig<double> cfg{
            .amplitude = 1.0,
            .hysteresis = 0.02,
            .setpoint = 0.0,
            .warmup_cycles = 2,
            .measure_cycles = 6,
            .period_tolerance = 0.05,
            .max_duration = 60.0,
        };
        constexpr auto         res = design::relay_autotune(cfg);
        RelayAutotuner<double> tuner(res, Ts);

        const auto out = run_autotune(plant_d, tuner, static_cast<int>(cfg.max_duration / Ts) + 100);

        REQUIRE(out.status == RelayAutotuneStatus::Done);

        const double     Tu_expected = 2.0 * std::numbers::pi / std::numbers::sqrt3;
        constexpr double Ku_expected = 8.0;

        // Tu is a zero-crossing measurement and lands within a few percent.
        // Ku rests on the relay's first-harmonic describing function
        // Kᵤ ≈ 4d/(πa), which under-predicts the true ultimate gain by
        // 15–25% on low-order plants because the higher harmonics of the
        // square-wave relay output aren't fully filtered by 1/(s+1)³ — the
        // measured peak `a` exceeds the pure-fundamental amplitude that the
        // describing function assumes. This is inherent to symmetric-relay
        // autotuning; the asymmetric / biased-relay variant
        // recovers Ku more tightly and gives the static gain for AMIGO.
        CHECK(out.Tu == doctest::Approx(Tu_expected).epsilon(0.10));
        CHECK(out.Ku == doctest::Approx(Ku_expected).epsilon(0.25));

        // Feed (Ku, Tu) into Tyreus-Luyben (the modern drop-in over the
        // original Ziegler-Nichols formulas) and sanity-check the gains.
        const auto pid = design::tyreus_luyben(out.Ku, out.Tu, Ts);
        CHECK(pid.Kp > 0.0);
        CHECK(pid.Ki > 0.0);
        CHECK(pid.Kd > 0.0);
    }

    TEST_CASE("reset returns a finished tuner to Idle") {
        constexpr design::RelayAutotuneConfig<float> cfg{};
        constexpr auto                               res = design::relay_autotune(cfg);
        RelayAutotuner<float>                        tuner(res, 0.01f);

        (void)tuner.step(0.0f);
        REQUIRE(tuner.status() != RelayAutotuneStatus::Idle);

        tuner.reset();
        CHECK(tuner.status() == RelayAutotuneStatus::Idle);
        CHECK(tuner.Ku() == 0.0f);
        CHECK(tuner.Tu() == 0.0f);
    }
} // TEST_SUITE
