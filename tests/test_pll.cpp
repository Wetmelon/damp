// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <numbers>

#include "damp/backend.hpp"
#include "damp/filters/pll.hpp"
#include "damp/math/math.hpp"
#include "damp/math/transforms.hpp"
#include "damp/matrix/matrix.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;


namespace {

// Build a balanced three-phase set at the given amplitude/phase and Clarke it.
AlphaBeta<float> grid_ab(float amp, float wt) {
    const float            tp = 2.0f * std::numbers::pi_v<float> / 3.0f;
    const ColVec<3, float> abc = {amp * std::cos(wt), amp * std::cos(wt - tp), amp * std::cos(wt + tp)};
    return clarke_transform(abc);
}

} // namespace


TEST_SUITE("Single-Phase PLL") {
    TEST_CASE("Constructor initializes limits and outputs") {
        constexpr float Fnom = 50.0f;

        constexpr SinglePhasePLL<float> pll(Fnom);

        // Loop-filter output (frequency offset about nominal) is clamped to ±50%.
        static_assert(pll.loop_filter.u_max == Fnom * 0.5f);
        static_assert(pll.loop_filter.u_min == -Fnom * 0.5f);
        static_assert(pll.loop_filter.Kbc == 10.0f); // Kbc = Kp anti-windup constant
        static_assert(pll.frequency() == Fnom);
        static_assert(pll.phase() == 0.0f);

        CHECK(pll.loop_filter.Kp == doctest::Approx(10.0f));
        CHECK(pll.loop_filter.Ki == doctest::Approx(100.0f));
        CHECK(pll.frequency() == doctest::Approx(Fnom));
        CHECK(pll.phase() == doctest::Approx(0.0f));
    }

    TEST_CASE("Step keeps estimates finite and within configured bounds") {
        constexpr float Fnom = 50.0f;
        constexpr float Ts = 0.0001f;

        SinglePhasePLL<float> pll(Fnom);
        const float           two_pi = 2.0f * std::numbers::pi_v<float>;

        for (int i = 0; i < 5000; ++i) {
            const float t = static_cast<float>(i) * Ts;
            const float input = std::sin(two_pi * Fnom * t);
            pll.step(input, Ts);
        }

        CHECK(std::isfinite(pll.frequency()));
        CHECK(std::isfinite(pll.phase()));
        CHECK(pll.frequency() <= Fnom * 1.5f); // nominal + max offset
        CHECK(pll.frequency() >= Fnom * 0.5f); // nominal − max offset
        CHECK(pll.phase() >= 0.0f);
        CHECK(pll.phase() < (2.0f * std::numbers::pi_v<float>)+1e-4f); // wrapped to [0, 2π)
        CHECK(pll.frequency() == doctest::Approx(Fnom).epsilon(0.25f));
    }

    TEST_CASE("Tracks an off-nominal tone in the correct direction") {
        // Regression guard for the SOGI-mixer phase-detector sign: a faster input
        // must drive the estimate UP toward it, not down to the lower frequency
        // rail (the bug the +input·quadrature sign produced).
        constexpr float Fnom = 50.0f;
        constexpr float Fin = 55.0f; // above nominal
        constexpr float Ts = 0.001f;
        const float     two_pi = 2.0f * std::numbers::pi_v<float>;

        SinglePhasePLL<float> pll(Fnom);

        float phase = 0.0f;
        for (int i = 0; i < 30000; ++i) {
            phase += two_pi * Fin * Ts;
            pll.step(std::sin(phase), Ts);
        }

        CHECK(pll.frequency() > Fnom + 1.0f); // moved toward the input, not the rail
        CHECK(pll.frequency() == doctest::Approx(Fin).epsilon(0.1f));
    }

    TEST_CASE("Integrator leak defaults to zero (pure integrator)") {
        constexpr float Fnom = 50.0f;

        constexpr SinglePhasePLL<float> pll(Fnom);
        static_assert(pll.integrator_leak == 0.0f);
    }

    TEST_CASE("Integrator leak bleeds a parked offset back toward zero") {
        constexpr float Fnom = 50.0f;
        constexpr float Fin = 55.0f; // off-nominal drive
        constexpr float Ts = 0.001f;
        const float     two_pi = 2.0f * std::numbers::pi_v<float>;

        SinglePhasePLL<float> pll(Fnom); // gentle default gains (stable for 1-phase)

        // Phase 1 (no leak): an off-nominal tone drives the estimator away from
        // nominal, parking a frequency offset in the integrator.
        float phase = 0.0f;
        for (int i = 0; i < 15000; ++i) {
            phase += two_pi * Fin * Ts;
            pll.step(std::sin(phase), Ts);
        }
        const float freq_charged = pll.frequency();
        CHECK(freq_charged > Fnom + 0.5f); // offset parked, and in the correct direction

        // Phase 2: enable the leak and remove excitation. phase_error → 0, so the
        // only dynamics are the leak, which pulls the estimate back to nominal.
        pll.integrator_leak = 50.0f; // [1/s] → time constant 20 ms
        for (int i = 0; i < 15000; ++i) {
            pll.step(0.0f, Ts);
        }
        const float freq_bled = pll.frequency();

        CHECK(damp::abs(freq_bled - Fnom) < damp::abs(freq_charged - Fnom));
        CHECK(freq_bled == doctest::Approx(Fnom).epsilon(0.05f));
    }

    TEST_CASE("Reset restores nominal frequency and zero phase") {
        constexpr float Fnom = 60.0f;
        constexpr float Ts = 0.0001f;

        SinglePhasePLL<float> pll(Fnom);

        // Drive the estimator away from initial state.
        for (int i = 0; i < 1000; ++i) {
            pll.step(0.5f, Ts);
        }

        pll.reset();

        CHECK(pll.frequency() == doctest::Approx(Fnom));
        CHECK(pll.phase() == doctest::Approx(0.0f));
    }
}

TEST_SUITE("Three-Phase PLL") {
    TEST_CASE("Locks frequency and phase on a balanced set") {
        constexpr float Fnom = 50.0f;
        constexpr float Ts = 0.0001f;
        constexpr float two_pi = 2.0f * std::numbers::pi_v<float>;
        constexpr float Fin = 55.0f; // off-nominal, within ±50% lock range

        ThreePhasePLL<float> pll(Fnom);
        pll.loop_filter.Kp = 40.0f;   // PI loop tuned for ~5 Hz natural frequency,
        pll.loop_filter.Ki = 1000.0f; // ζ≈0.6 (defaults are too slow to pull 5 Hz)

        float phase = 0.7f; // arbitrary input phase offset
        for (int i = 0; i < 20000; ++i) {
            phase += two_pi * Fin * Ts;
            const ColVec<3, float> abc = {
                std::cos(phase),
                std::cos(phase - (two_pi / 3.0f)),
                std::cos(phase + (two_pi / 3.0f)),
            };
            pll.step(abc, Ts);
        }

        CHECK(pll.frequency() == doctest::Approx(Fin).epsilon(0.01f));

        // Phase estimate tracks the input phase (compare on the circle so the
        // cos residual is wrap-invariant (phase may sit in [0, 2π)).
        CHECK(std::cos(pll.phase() - phase) == doctest::Approx(1.0f).epsilon(0.02f));
    }

    TEST_CASE("Reset restores nominal frequency and zero phase") {
        constexpr float Fnom = 60.0f;
        constexpr float Ts = 0.0001f;

        ThreePhasePLL<float> pll(Fnom);
        pll.step({1.0f, -0.5f, -0.5f}, Ts);
        pll.reset();

        CHECK(pll.frequency() == doctest::Approx(Fnom));
        CHECK(pll.phase() == doctest::Approx(0.0f));
    }
}

TEST_SUITE("DSOGI-PLL") {

    TEST_CASE("Instantaneous sequence calculator: balanced input is pure positive") {
        // v = (1, 0); its 90° lag qv = (0, -1) for a positive-sequence rotation.
        const AlphaBeta<float> v = {1.0f, 0.0f};
        const AlphaBeta<float> qv = {0.0f, -1.0f};

        const auto pos = positive_sequence_ab(v, qv);
        const auto neg = negative_sequence_ab(v, qv);

        CHECK(pos.alpha == doctest::Approx(1.0f));
        CHECK(pos.beta == doctest::Approx(0.0f));
        CHECK(neg.abs() == doctest::Approx(0.0f).epsilon(1e-6f));
    }

    TEST_CASE("Locks frequency and phase on a balanced 50 Hz grid") {
        const float f0 = 50.0f;
        const float Ts = 1.0f / 10000.0f; // 10 kHz
        const float w = 2.0f * std::numbers::pi_v<float> * f0;

        DsogiPll<float> pll(f0);

        // Run for ~0.5 s to settle.
        float wt = 0.0f;
        for (int i = 0; i < 5000; ++i) {
            pll.step(grid_ab(1.0f, wt), Ts);
            wt += w * Ts;
        }

        CHECK(pll.frequency() == doctest::Approx(f0).epsilon(0.01f));
        // At lock the positive-sequence dq aligns: d → amplitude, q → 0.
        CHECK(pll.positive_dq().d == doctest::Approx(1.0f).epsilon(0.02f));
        CHECK(pll.positive_dq().q == doctest::Approx(0.0f).epsilon(0.02f));
        // Balanced input ⇒ negligible negative sequence.
        CHECK(pll.negative_sequence().abs() == doctest::Approx(0.0f).epsilon(0.02f));
    }

    TEST_CASE("Extracts positive sequence from an unbalanced grid") {
        const float f0 = 50.0f;
        const float Ts = 1.0f / 10000.0f;
        const float w = 2.0f * std::numbers::pi_v<float> * f0;

        DsogiPll<float> pll(f0);

        // Unbalanced: inject a negative-sequence component by scaling one phase.
        const float tp = 2.0f * std::numbers::pi_v<float> / 3.0f;
        float       wt = 0.0f;
        for (int i = 0; i < 6000; ++i) {
            const ColVec<3, float> abc = {std::cos(wt), 0.6f * std::cos(wt - tp), std::cos(wt + tp)};
            pll.step(clarke_transform(abc), Ts);
            wt += w * Ts;
        }

        // Still locks to the line frequency despite unbalance.
        CHECK(pll.frequency() == doctest::Approx(f0).epsilon(0.02f));
        // Negative sequence is non-trivial (unbalance is detected).
        CHECK(pll.negative_sequence().abs() > 0.05f);
    }

} // TEST_SUITE
