// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>

#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/math/transforms.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/motor/modulation.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using namespace damp::motor;

/**
 * @brief Tests for reference-frame transforms and power-electronics modulation
 */

TEST_SUITE("Transforms & Modulation") {
    TEST_CASE("Clarke transform") {
        // Test balanced three-phase
        const ColVec<3, float> abc = {1.0f, -0.5f, -0.5f};

        const auto [alpha, beta] = clarke_transform(abc);

        // For balanced system: α = (2a - b - c)/3, β = (b - c)/√3
        CHECK(alpha == doctest::Approx(1.0f));
        CHECK(beta == doctest::Approx(0.0f).epsilon(1e-6f));
    }

    TEST_CASE("Inverse Clarke transform") {
        // Test round-trip
        const ColVec<3, float> abc_orig = {1.0f, -0.5f, -0.5f};

        const auto ab = clarke_transform(abc_orig);
        const auto abc = inverse_clarke_transform(ab);

        CHECK(abc[0] == doctest::Approx(abc_orig[0]).epsilon(1e-6f));
        CHECK(abc[1] == doctest::Approx(abc_orig[1]).epsilon(1e-6f));
        CHECK(abc[2] == doctest::Approx(abc_orig[2]).epsilon(1e-6f));
    }

    TEST_CASE("Park transform") {
        // Test with θ = 0 (should be identity)
        const AlphaBeta<float> ab = {.alpha = 1.0f, .beta = 0.5f};
        const float            theta = 0.0f;

        const auto [d, q] = park_transform(ab, theta);

        CHECK(d == doctest::Approx(ab.alpha));
        CHECK(q == doctest::Approx(ab.beta));
    }

    TEST_CASE("Park transform with rotation") {
        const AlphaBeta<float> ab = {.alpha = 1.0f, .beta = 0.0f};
        const float            theta = std::numbers::pi_v<float> / 4.0f; // 45°

        const auto [d, q] = park_transform(ab, theta);

        const float expected_d = ab.alpha * std::cos(theta);
        const float expected_q = -ab.alpha * std::sin(theta);

        CHECK(d == doctest::Approx(expected_d).epsilon(1e-6f));
        CHECK(q == doctest::Approx(expected_q).epsilon(1e-6f));
    }

    TEST_CASE("Park d-axis convention: space vector at theta maps to ( |v|, 0 )") {
        // θ is the electrical angle of the d-axis: a pure αβ vector at angle θ
        // lands entirely on d (q leads d by +90° electrical).
        const float            theta = std::numbers::pi_v<float> / 3.0f;
        const float            mag = 2.0f;
        const AlphaBeta<float> ab = {
            .alpha = mag * std::cos(theta),
            .beta = mag * std::sin(theta),
        };
        const auto [d, q] = park_transform(ab, theta);
        CHECK(d == doctest::Approx(mag).epsilon(1e-6f));
        CHECK(q == doctest::Approx(0.0f).epsilon(1e-6f));
        // q-lead: rotate the space vector +90° → pure +q.
        const AlphaBeta<float> ab_q = {
            .alpha = mag * std::cos(theta + std::numbers::pi_v<float> / 2.0f),
            .beta = mag * std::sin(theta + std::numbers::pi_v<float> / 2.0f),
        };
        const auto dq_q = park_transform(ab_q, theta);
        CHECK(dq_q.d == doctest::Approx(0.0f).epsilon(1e-6f));
        CHECK(dq_q.q == doctest::Approx(mag).epsilon(1e-6f));
    }

    TEST_CASE("Inverse Park transform") {
        // Test round-trip
        const DirectQuadrature<float> dq_orig = {.d = 0.8f, .q = 0.3f};
        const float                   theta = std::numbers::pi_v<float> / 6.0f; // 30°

        const auto ab = inverse_park_transform(dq_orig, theta);
        const auto [d, q] = park_transform(ab, theta);

        CHECK(d == doctest::Approx(dq_orig.d).epsilon(1e-6f));
        CHECK(q == doctest::Approx(dq_orig.q).epsilon(1e-6f));
    }

    TEST_CASE("Clarke-Park combined transform") {
        // Test three-phase to dq
        const ColVec<3, float> abc = {
            std::cos(0.0f),
            std::cos(2.0f * std::numbers::pi_v<float> / 3.0f),
            std::cos(4.0f * std::numbers::pi_v<float> / 3.0f),
        };
        const float theta = 0.0f;

        const auto [d, q] = clarke_park_transform(abc, theta);

        // At θ = 0, d should be the amplitude, q should be 0
        CHECK(d == doctest::Approx(1.0f).epsilon(1e-6f)); // Clarke transform normalizes to 1.0
        CHECK(q == doctest::Approx(0.0f).epsilon(1e-6f));
    }

    TEST_CASE("Inverse Park-Clarke combined transform") {
        // Test dq to three-phase round-trip
        const DirectQuadrature<float> dq = {.d = 1.0f, .q = 0.5f};
        const float                   theta = std::numbers::pi_v<float> / 4.0f;

        const auto abc = inverse_park_clarke_transform(dq, theta);
        const auto [d2, q2] = clarke_park_transform(abc, theta);

        CHECK(d2 == doctest::Approx(dq.d).epsilon(1e-6f));
        CHECK(q2 == doctest::Approx(dq.q).epsilon(1e-6f));
    }

    TEST_CASE("SVM duty cycles") {
        // Test zero voltage
        const auto svm = svm_duty_cycles<float>({.alpha = 0.0f, .beta = 0.0f}, 100.0f);

        CHECK(svm.duties[0] == doctest::Approx(0.5f));
        CHECK(svm.duties[1] == doctest::Approx(0.5f));
        CHECK(svm.duties[2] == doctest::Approx(0.5f));
        CHECK_FALSE(svm.is_clipped);

        // Test maximum linear voltage (peak phase = Vdc/√3 with SVPWM injection)
        const float v_max = 100.0f / std::numbers::sqrt3_v<float>;
        const auto  svm_max = svm_duty_cycles<float>({.alpha = v_max, .beta = 0.0f}, 100.0f);

        CHECK(svm_max.duties[0] >= 0.0f);
        CHECK(svm_max.duties[0] <= 1.0f);
        CHECK(svm_max.duties[1] >= 0.0f);
        CHECK(svm_max.duties[1] <= 1.0f);
        CHECK(svm_max.duties[2] >= 0.0f);
        CHECK(svm_max.duties[2] <= 1.0f);
        CHECK_FALSE(svm_max.is_clipped); // exactly on the inscribed circle, not clipped
    }

    TEST_CASE("SVM Vdc <= 0 returns mid-rail duties (finite, no Inf)") {
        const AlphaBeta<float> v = {.alpha = 10.0f, .beta = -5.0f};
        for (float vdc : {0.0f, -48.0f}) {
            const auto svm = svm_duty_cycles<float>(v, vdc);
            CHECK(svm.duties[0] == doctest::Approx(0.5f));
            CHECK(svm.duties[1] == doctest::Approx(0.5f));
            CHECK(svm.duties[2] == doctest::Approx(0.5f));
            CHECK(svm.is_clipped);
            CHECK(std::isfinite(svm.duties[0]));
            CHECK(std::isfinite(svm.duties[1]));
            CHECK(std::isfinite(svm.duties[2]));
        }
    }

    TEST_CASE("modulation_duty_cycles(Svpwm) matches svm_duty_cycles") {
        const float            Vdc = 100.0f;
        const AlphaBeta<float> cmds[] = {
            {.alpha = 0.0f, .beta = 0.0f},
            {.alpha = 20.0f, .beta = -10.0f},
            {.alpha = Vdc / std::numbers::sqrt3_v<float>, .beta = 0.0f},
            {.alpha = -15.0f, .beta = 40.0f},
        };
        for (const auto& v : cmds) {
            const auto a = svm_duty_cycles<float>(v, Vdc);
            const auto b = modulation_duty_cycles<float>(v, Vdc, PwmScheme::Svpwm);
            for (size_t i = 0; i < 3; ++i) {
                CHECK(a.duties[i] == doctest::Approx(b.duties[i]).epsilon(1e-6f));
            }
            CHECK(a.is_clipped == b.is_clipped);
        }
    }

    TEST_CASE("SPWM has zero common-mode; line-line matches SVPWM in linear range") {
        const float            Vdc = 48.0f;
        const float            Vpk = 15.0f; // well inside SPWM limit Vdc/2
        const AlphaBeta<float> v = {.alpha = Vpk, .beta = 0.0f};

        const auto spwm = modulation_duty_cycles<float>(v, Vdc, PwmScheme::Spwm);
        const auto svpwm = modulation_duty_cycles<float>(v, Vdc, PwmScheme::Svpwm);

        // SPWM: d = 0.5 + v_phase/Vdc (no ZS) — phase A peak.
        CHECK(spwm.duties[0] == doctest::Approx(0.5f + Vpk / Vdc).epsilon(1e-5f));
        CHECK_FALSE(spwm.is_clipped);
        CHECK_FALSE(svpwm.is_clipped);

        // Line-to-line (duty differences) are ZS-invariant for continuous schemes.
        for (size_t i = 0; i < 3; ++i) {
            const size_t j = (i + 1) % 3;
            CHECK((spwm.duties[i] - spwm.duties[j]) == doctest::Approx(svpwm.duties[i] - svpwm.duties[j]).epsilon(1e-5f));
        }
    }

    TEST_CASE("THIPWM 1/6 injection and line-line identity vs SVPWM") {
        const float Vdc = 100.0f;
        const float Vpk = 40.0f; // linear for both THIPWM and SVPWM
        // θ = 0 → v0 = Vpk/6
        {
            const AlphaBeta<float> v = {.alpha = Vpk, .beta = 0.0f};
            const float            v0 = thipwm_zero_sequence(v);
            CHECK(v0 == doctest::Approx(Vpk / 6.0f).epsilon(1e-5f));

            const auto thi = modulation_duty_cycles<float>(v, Vdc, PwmScheme::Thipwm);
            const auto svm = modulation_duty_cycles<float>(v, Vdc, PwmScheme::Svpwm);
            CHECK_FALSE(thi.is_clipped);
            for (size_t i = 0; i < 3; ++i) {
                const size_t j = (i + 1) % 3;
                CHECK((thi.duties[i] - thi.duties[j]) == doctest::Approx(svm.duties[i] - svm.duties[j]).epsilon(1e-5f));
            }
        }
        // θ = 90° → cos(3θ) = cos(270°) = 0
        {
            const AlphaBeta<float> v = {.alpha = 0.0f, .beta = Vpk};
            CHECK(thipwm_zero_sequence(v) == doctest::Approx(0.0f).epsilon(1e-5f));
        }
        // θ = 60° → cos(180°) = -1 → v0 = -Vpk/6
        {
            const float            c60 = 0.5f;
            const float            s60 = std::numbers::sqrt3_v<float> / 2.0f;
            const AlphaBeta<float> v = {.alpha = Vpk * c60, .beta = Vpk * s60};
            CHECK(thipwm_zero_sequence(v) == doctest::Approx(-Vpk / 6.0f).epsilon(1e-5f));
        }
    }

    TEST_CASE("DPWMMAX clamps max duty to 1; DPWMMIN clamps min duty to 0") {
        const float Vdc = 100.0f;
        // Sweep electrical angle at moderate MI (linear hexagon interior).
        const float Vpk = 30.0f;
        const int   N = 24;
        for (int k = 0; k < N; ++k) {
            const float            th = 2.0f * std::numbers::pi_v<float> * static_cast<float>(k) / static_cast<float>(N);
            const AlphaBeta<float> v = {.alpha = Vpk * std::cos(th), .beta = Vpk * std::sin(th)};

            const auto dmax = modulation_duty_cycles<float>(v, Vdc, PwmScheme::DpwmMax);
            const auto dmin = modulation_duty_cycles<float>(v, Vdc, PwmScheme::DpwmMin);

            const float mx = std::max({dmax.duties[0], dmax.duties[1], dmax.duties[2]});
            const float mn = std::min({dmin.duties[0], dmin.duties[1], dmin.duties[2]});
            CHECK(mx == doctest::Approx(1.0f).epsilon(1e-5f));
            CHECK(mn == doctest::Approx(0.0f).epsilon(1e-5f));
            CHECK_FALSE(dmax.is_clipped);
            CHECK_FALSE(dmin.is_clipped);

            // Line-line preserved vs SVPWM (no overmod).
            const auto svm = svm_duty_cycles<float>(v, Vdc);
            for (size_t i = 0; i < 3; ++i) {
                const size_t j = (i + 1) % 3;
                CHECK((dmax.duties[i] - dmax.duties[j]) == doctest::Approx(svm.duties[i] - svm.duties[j]).epsilon(1e-5f));
                CHECK((dmin.duties[i] - dmin.duties[j]) == doctest::Approx(svm.duties[i] - svm.duties[j]).epsilon(1e-5f));
            }
        }
    }

    TEST_CASE("DPWM1 clamps either max to 1 or min to 0 each sample") {
        const float Vdc = 100.0f;
        const float Vpk = 35.0f;
        const int   N = 36;
        int         high_clamps = 0;
        int         low_clamps = 0;
        for (int k = 0; k < N; ++k) {
            const float            th = 2.0f * std::numbers::pi_v<float> * static_cast<float>(k) / static_cast<float>(N);
            const AlphaBeta<float> v = {.alpha = Vpk * std::cos(th), .beta = Vpk * std::sin(th)};
            const auto             d = modulation_duty_cycles<float>(v, Vdc, PwmScheme::Dpwm1);

            const float mx = std::max({d.duties[0], d.duties[1], d.duties[2]});
            const float mn = std::min({d.duties[0], d.duties[1], d.duties[2]});
            const bool  hi = std::abs(mx - 1.0f) < 1e-4f;
            const bool  lo = std::abs(mn - 0.0f) < 1e-4f;
            const bool  rail = hi || lo;
            CHECK(rail);
            if (hi) {
                ++high_clamps;
            }
            if (lo) {
                ++low_clamps;
            }
            CHECK_FALSE(d.is_clipped);
        }
        // Both polarities appear over a full fundamental (60° windows).
        CHECK(high_clamps > 0);
        CHECK(low_clamps > 0);
    }

    TEST_CASE("DPWM0/2/3 always rail-clamp one phase in linear range") {
        const float     Vdc = 80.0f;
        const float     Vpk = 25.0f;
        const PwmScheme schemes[] = {PwmScheme::Dpwm0, PwmScheme::Dpwm2, PwmScheme::Dpwm3};
        const int       N = 24;
        for (const auto scheme : schemes) {
            for (int k = 0; k < N; ++k) {
                const float            th = 2.0f * std::numbers::pi_v<float> * static_cast<float>(k) / static_cast<float>(N);
                const AlphaBeta<float> v = {.alpha = Vpk * std::cos(th), .beta = Vpk * std::sin(th)};
                const auto             d = modulation_duty_cycles<float>(v, Vdc, scheme);
                const float            mx = std::max({d.duties[0], d.duties[1], d.duties[2]});
                const float            mn = std::min({d.duties[0], d.duties[1], d.duties[2]});
                const bool             rail = (std::abs(mx - 1.0f) < 1e-4f) || (std::abs(mn - 0.0f) < 1e-4f);
                CHECK(rail);
                CHECK_FALSE(d.is_clipped);
            }
        }
    }

    TEST_CASE("DPWM1 peak angle clamps phase A high") {
        // θ = 0: phase A is the positive peak → DPWM1 high-clamps A.
        const float            Vdc = 100.0f;
        const float            Vpk = 40.0f;
        const AlphaBeta<float> v = {.alpha = Vpk, .beta = 0.0f};
        const auto             d = modulation_duty_cycles<float>(v, Vdc, PwmScheme::Dpwm1);
        CHECK(d.duties[0] == doctest::Approx(1.0f).epsilon(1e-5f));
        CHECK(d.duties[0] >= d.duties[1]);
        CHECK(d.duties[0] >= d.duties[2]);
    }

    TEST_CASE("modulation_index and linear/six-step voltage helpers") {
        const float Vdc = 100.0f;
        CHECK(linear_modulation_voltage(Vdc) == doctest::Approx(Vdc / std::numbers::sqrt3_v<float>).epsilon(1e-6f));
        CHECK(six_step_fundamental_voltage(Vdc) == doctest::Approx(2.0f * Vdc / std::numbers::pi_v<float>).epsilon(1e-6f));

        const float            Vlin = linear_modulation_voltage(Vdc);
        const AlphaBeta<float> v_lin = {.alpha = Vlin, .beta = 0.0f};
        CHECK(modulation_index(v_lin, Vdc) == doctest::Approx(1.0f).epsilon(1e-5f));
        CHECK(modulation_index({.alpha = 0.0f, .beta = 0.0f}, Vdc) == doctest::Approx(0.0f));
        CHECK(modulation_index(v_lin, 0.0f) == doctest::Approx(0.0f));

        // Inside circle: not clipped; well outside: clipped.
        CHECK_FALSE(svm_duty_cycles<float>({.alpha = 0.5f * Vlin, .beta = 0.0f}, Vdc).is_clipped);
        CHECK(svm_duty_cycles<float>({.alpha = 1.2f * Vlin, .beta = 0.0f}, Vdc).is_clipped);
    }

    TEST_CASE("six_step_duty_cycles: 180° square per phase, 120° apart") {
        // θ = 0 → A high, B low, C low (cos(0)>0, cos(-120)<0, cos(+120)<0)
        {
            const auto d = six_step_duty_cycles<float>({.alpha = 1.0f, .beta = 0.0f});
            CHECK(d.duties[0] == doctest::Approx(1.0f));
            CHECK(d.duties[1] == doctest::Approx(0.0f));
            CHECK(d.duties[2] == doctest::Approx(0.0f));
            CHECK(d.is_clipped);
        }
        // θ = 60° → A: cos(60)>0, B: cos(-60)>0, C: cos(180)<0 → {1,1,0}
        {
            const float c60 = 0.5f;
            const float s60 = std::numbers::sqrt3_v<float> / 2.0f;
            const auto  d = six_step_duty_cycles<float>({.alpha = c60, .beta = s60});
            CHECK(d.duties[0] == doctest::Approx(1.0f));
            CHECK(d.duties[1] == doctest::Approx(1.0f));
            CHECK(d.duties[2] == doctest::Approx(0.0f));
        }
        // Zero vector → mid-rail, still flagged non-linear
        {
            const auto d = six_step_duty_cycles<float>({.alpha = 0.0f, .beta = 0.0f});
            CHECK(d.duties[0] == doctest::Approx(0.5f));
            CHECK(d.is_clipped);
        }
    }

    TEST_CASE("Modulator wrapper forwards scheme") {
        Modulator<float>       mod{.scheme = PwmScheme::Thipwm};
        const float            Vdc = 50.0f;
        const AlphaBeta<float> v = {.alpha = 10.0f, .beta = -5.0f};
        const auto             a = mod.duties(v, Vdc);
        const auto             b = modulation_duty_cycles<float>(v, Vdc, PwmScheme::Thipwm);
        for (size_t i = 0; i < 3; ++i) {
            CHECK(a.duties[i] == doctest::Approx(b.duties[i]));
        }
    }

    TEST_CASE("constexpr modulation maps are compile-time evaluable") {
        constexpr AlphaBeta<double> v{.alpha = 10.0, .beta = 0.0};
        constexpr double            Vdc = 48.0;
        constexpr auto              svm = svm_duty_cycles(v, Vdc);
        constexpr auto              thi = modulation_duty_cycles(v, Vdc, PwmScheme::Thipwm);
        constexpr auto              dp1 = modulation_duty_cycles(v, Vdc, PwmScheme::Dpwm1);
        constexpr auto              mi = modulation_index(v, Vdc);
        static_assert(svm.duties[0] > 0.5);
        static_assert(thi.duties[0] > 0.5);
        static_assert(dp1.duties[0] == 1.0);
        static_assert(mi > 0.0);
        CHECK(svm.duties[0] > 0.5);
        CHECK(dp1.duties[0] == doctest::Approx(1.0));
        CHECK(mi > 0.0);
    }

    TEST_CASE("Zero sequence: common-mode DC lands in zero, not αβ") {
        // Identical DC bias on all three phases is pure common-mode.
        const float            d = 0.7f;
        const ColVec<3, float> abc = {d, d, d};
        const auto             ab = clarke_transform<float>(abc);

        CHECK(ab.alpha == doctest::Approx(0.0f).epsilon(1e-6f));
        CHECK(ab.beta == doctest::Approx(0.0f).epsilon(1e-6f));
        CHECK(zero_sequence<float>(abc) == doctest::Approx(d).epsilon(1e-6f));
    }

    TEST_CASE("Zero sequence: per-phase DC leaks into αβ") {
        // Offset on phase a only: 2d/3 into alpha, d/3 into zero.
        const float            d = 0.9f;
        const ColVec<3, float> abc = {d, 0.0f, 0.0f};
        const auto             ab = clarke_transform<float>(abc);

        CHECK(ab.alpha == doctest::Approx(2.0f * d / 3.0f).epsilon(1e-6f));
        CHECK(ab.beta == doctest::Approx(0.0f).epsilon(1e-6f));
        CHECK(zero_sequence<float>(abc) == doctest::Approx(d / 3.0f).epsilon(1e-6f));
    }

    TEST_CASE("Clarke round-trip abc → αβ + zero → abc") {
        const ColVec<3, float> abc = {1.3f, -0.4f, 0.2f};
        const auto             ab = clarke_transform<float>(abc);
        const auto             z = zero_sequence<float>(abc);
        const auto             rt = inverse_clarke_transform(ab, z);

        for (size_t i = 0; i < 3; ++i) {
            CHECK(rt[i] == doctest::Approx(abc[i]).epsilon(1e-6f));
        }
    }

    TEST_CASE("Symmetrical components: balanced set is pure positive sequence") {
        using Cplx = damp::complex<float>;
        // Balanced positive-sequence phasors: a, a·e^{-j120}, a·e^{-j240}
        const Cplx a = {1.0f, 0.0f};
        const Cplx b = {-0.5f, -std::numbers::sqrt3_v<float> / 2.0f}; // 1∠-120°
        const Cplx c = {-0.5f, std::numbers::sqrt3_v<float> / 2.0f};  // 1∠-240°

        const auto seq = symmetrical_components<float>({a, b, c});

        CHECK(seq.zero.abs() == doctest::Approx(0.0f).epsilon(1e-6f));
        CHECK(seq.negative.abs() == doctest::Approx(0.0f).epsilon(1e-6f));
        CHECK(seq.positive.real() == doctest::Approx(1.0f).epsilon(1e-6f));
        CHECK(seq.positive.imag() == doctest::Approx(0.0f).epsilon(1e-6f));
    }

    TEST_CASE("Symmetrical components: round-trip 012 → abc → 012") {
        using Cplx = damp::complex<float>;
        const ColVec<3, Cplx> abc = {Cplx{1.0f, 0.2f}, Cplx{-0.4f, -0.9f}, Cplx{0.1f, 0.7f}};

        const auto seq = symmetrical_components<float>(abc);
        const auto rt = inverse_symmetrical_components<float>(seq);

        for (size_t i = 0; i < 3; ++i) {
            CHECK(rt[i].real() == doctest::Approx(abc[i].real()).epsilon(1e-6f));
            CHECK(rt[i].imag() == doctest::Approx(abc[i].imag()).epsilon(1e-6f));
        }
    }

    TEST_CASE("Power-invariant Clarke: magnitude scaling and round-trip") {
        using damp::Convention;
        const ColVec<3, float> abc = {1.0f, -0.4f, 0.3f};

        const auto amp = clarke_transform<float>(abc); // amplitude-invariant (default)
        const auto pwr = clarke_transform<float, Convention::PowerInvariant>(abc);

        // Power-invariant αβ is √(3/2)× the amplitude-invariant αβ.
        const float ratio = std::sqrt(1.5f);
        CHECK(pwr.alpha == doctest::Approx(amp.alpha * ratio).epsilon(1e-6f));
        CHECK(pwr.beta == doctest::Approx(amp.beta * ratio).epsilon(1e-6f));

        // Round-trip through the power-invariant inverse with zero (exact).
        const auto z = zero_sequence<float, Convention::PowerInvariant>(abc);
        const auto rt = inverse_clarke_transform(pwr, z);
        for (size_t i = 0; i < 3; ++i) {
            CHECK(rt[i] == doctest::Approx(abc[i]).epsilon(1e-6f));
        }
    }

    TEST_CASE("Fused Park-Clarke round-trip with scalar zero passthrough") {
        const ColVec<3, float> abc = {1.3f, -0.4f, 0.2f};
        const float            theta = 0.9f;

        const auto dq = clarke_park_transform<float>(abc, theta);
        const auto z = zero_sequence<float>(abc);
        const auto rt = inverse_park_clarke_transform(dq, theta, z);
        for (size_t i = 0; i < 3; ++i) {
            CHECK(rt[i] == doctest::Approx(abc[i]).epsilon(1e-5f));
        }

        // Zero is rotation-invariant: same scalar before/after Park.
        CHECK(z == doctest::Approx(zero_sequence<float>(rt)).epsilon(1e-5f));
        // Three-wire inverse (zero = 0) drops residual common-mode.
        const auto rt0 = inverse_park_clarke_transform(dq, theta);
        CHECK(damp::abs(rt0[0] - abc[0]) > 0.1f);
    }

    TEST_CASE("Instantaneous power: balanced load, αβ and dq agree") {
        // Balanced 3φ voltage, unity-PF current (in phase), amplitude 1 each.
        const float            tp = 2.0f * std::numbers::pi_v<float> / 3.0f;
        const float            th = 0.37f;
        const ColVec<3, float> v = {std::cos(th), std::cos(th - tp), std::cos(th + tp)};
        const ColVec<3, float> i = v; // unity power factor, equal amplitude

        const auto p_ab = instantaneous_power<float>(clarke_transform<float>(v), clarke_transform<float>(i));

        // P = 3/2 · V·I (peak) = 1.5 for unit peak, balanced, unity PF; Q = 0.
        CHECK(p_ab.p == doctest::Approx(1.5f).epsilon(1e-5f));
        CHECK(p_ab.q == doctest::Approx(0.0f).epsilon(1e-5f));

        // dq form agrees (frame-invariant scalars).
        const auto p_dq = instantaneous_power<float>(clarke_park_transform<float>(v, th), clarke_park_transform<float>(i, th));
        CHECK(p_dq.p == doctest::Approx(p_ab.p).epsilon(1e-5f));
        CHECK(p_dq.q == doctest::Approx(p_ab.q).epsilon(1e-5f));
    }

    TEST_CASE("Instantaneous power: convention scaling differs by 3/2") {
        using damp::Convention;
        const ColVec<3, float> v = {1.0f, -0.5f, -0.5f};
        const ColVec<3, float> i = {0.8f, -0.2f, -0.6f};

        const auto p_amp = instantaneous_power<float>(clarke_transform<float>(v), clarke_transform<float>(i));
        const auto vp = clarke_transform<float, Convention::PowerInvariant>(v);
        const auto ip = clarke_transform<float, Convention::PowerInvariant>(i);
        const auto p_pwr = instantaneous_power<float, Convention::PowerInvariant>(vp, ip);

        // Both report the same real watts.
        CHECK(p_pwr.p == doctest::Approx(p_amp.p).epsilon(1e-5f));
        CHECK(p_pwr.q == doctest::Approx(p_amp.q).epsilon(1e-5f));
    }

    TEST_CASE("InstantaneousPower: apparent power, angle, power factor") {
        // v on d-axis, current at 45° (equal active/reactive) → φ = -45°.
        const auto s = instantaneous_power<float>(DirectQuadrature<float>{1.0f, 0.0f}, DirectQuadrature<float>{1.0f, 1.0f});

        CHECK(s.p == doctest::Approx(1.5f));
        CHECK(s.q == doctest::Approx(-1.5f));
        CHECK(s.abs() == doctest::Approx(1.5f * std::numbers::sqrt2_v<float>));
        CHECK(s.arg() == doctest::Approx(-std::numbers::pi_v<float> / 4.0f));
        CHECK(s.power_factor() == doctest::Approx(1.0f / std::numbers::sqrt2_v<float>));

        // Degenerate: zero power → zero pf, no division by zero.
        CHECK(InstantaneousPower<float>{}.power_factor() == doctest::Approx(0.0f));
    }

    TEST_CASE("Complex power via αβ conj/product matches instantaneous_power") {
        const AlphaBeta<float> v = {0.9f, -0.3f};
        const AlphaBeta<float> i = {0.4f, 0.7f};

        // S = 3/2 · V · conj(I), expressed directly with the complex operators.
        const AlphaBeta<float> s = (v * i.conj()) * 1.5f;
        const auto             p = instantaneous_power<float>(v, i);

        CHECK(s.alpha == doctest::Approx(p.p));
        CHECK(s.beta == doctest::Approx(p.q));

        // conj() and complex operator* sanity: j · conj(j) = j · (−j) = 1.
        const AlphaBeta<float> j = {0.0f, 1.0f};
        const AlphaBeta<float> one = j * j.conj();
        CHECK(one.alpha == doctest::Approx(1.0f));
        CHECK(one.beta == doctest::Approx(0.0f));
    }

    TEST_CASE("Modulator default SVPWM matches svm_duty_cycles; Spwm differs") {
        using damp::motor::Modulator;
        using damp::motor::PwmScheme;
        using damp::motor::svm_duty_cycles;

        const AlphaBeta<float> Vab{10.0f, 5.0f};
        const float            Vdc = 48.0f;

        Modulator<float> mod{}; // default Svpwm
        const auto       out_sv = mod.duties(Vab, Vdc);
        const auto       ref = svm_duty_cycles(Vab, Vdc);
        for (std::size_t i = 0; i < 3; ++i) {
            CHECK(out_sv.duties[i] == doctest::Approx(ref.duties[i]));
            CHECK(out_sv.duties[i] >= 0.0f);
            CHECK(out_sv.duties[i] <= 1.0f);
        }

        mod.scheme = PwmScheme::Spwm;
        const auto out_sp = mod.duties(Vab, Vdc);
        const bool differ = (out_sv.duties[0] != out_sp.duties[0]) || (out_sv.duties[1] != out_sp.duties[1])
                         || (out_sv.duties[2] != out_sp.duties[2]);
        CHECK(differ);
    }

} // TEST_SUITE
