// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "damp/backend.hpp"
#include "damp/estimation/frequency_response.hpp"
#include "damp/estimation/successive_compensator.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

TEST_SUITE("Second-order step fit") {
    TEST_CASE("recovers underdamped ring from analytic samples") {
        // Free response of underdamped mode (impulse residual shape):
        // e(t) = A e^{-ζ ωn t} cos(ωd t)
        constexpr double       fn = 20.0;
        constexpr double       wn = 2.0 * damp::numbers::pi_v<double> * fn;
        constexpr double       zeta = 0.08;
        constexpr double       wd = wn * std::sqrt(1.0 - (zeta * zeta));
        constexpr double       Ts = 0.001;
        constexpr std::size_t  N = 800;
        damp::array<double, N> y{};
        const double           y_ss = 1.0;
        for (std::size_t k = 0; k < N; ++k) {
            const double t = static_cast<double>(k) * Ts;
            const double env = std::exp(-zeta * wn * t);
            y[k] = y_ss + (0.4 * env * std::cos(wd * t));
        }
        const auto mode = design::fit_second_order_step(y, Ts, 0.01);
        REQUIRE(mode.valid);
        CHECK(mode.frequency_hz == doctest::Approx(fn).epsilon(0.08));
        CHECK(mode.zeta == doctest::Approx(zeta).epsilon(0.35));
    }
}

TEST_SUITE("Successive notch commissioner") {
    TEST_CASE("step capture recovers resonance and installs a deep notch") {
        constexpr double Ts = 0.001;
        constexpr double fn = 25.0;
        constexpr double wn = 2.0 * damp::numbers::pi_v<double> * fn;
        constexpr double zet = 0.05;

        // Plant: wn² / (s² + 2ζwn s + wn²)
        TransferFunction<1, 3, double> tf{
            .num = {wn * wn},
            .den = {wn * wn, 2.0 * zet * wn, 1.0},
        };
        const auto sys_d = *discretize(*tf.to_state_space(), Ts, DiscretizationMethod::ZOH);

        SuccessiveCompensatorCommissioner<2, double> cm{SuccessiveCompensatorConfig<double>{
            .Ts = Ts,
            .notch_Q = 10.0,
            .min_peak = 0.02,
            .min_frequency_hz = 5.0,
            .max_zeta = 0.4,
        }};

        constexpr std::size_t  N = 1000;
        damp::array<double, N> y{};
        ColVec<2, double>      x{};
        for (std::size_t k = 0; k < N; ++k) {
            const double u = 1.0;
            const auto   yv = sys_d.C * x + sys_d.D * ColVec<1, double>{u};
            y[k] = yv(0, 0);
            x = sys_d.A * x + sys_d.B * ColVec<1, double>{u};
        }

        REQUIRE(cm.assign_from_step_capture(y));
        CHECK(cm.count() == 1);
        CHECK(cm.mode(0).frequency_hz == doctest::Approx(fn).epsilon(0.12));

        // Notch loaded into cascade: drive sine at the *fitted* frequency; after
        // settle the peak should be small (RBJ null at f0).
        const double f_notch = cm.mode(0).frequency_hz;
        double       peak = 0.0;
        for (int k = 0; k < 4000; ++k) {
            const double t = static_cast<double>(k) * Ts;
            const double u = std::sin(2.0 * damp::numbers::pi_v<double> * f_notch * t);
            const double v = cm.compensate(u);
            if (k > 2000) {
                peak = std::max(peak, std::abs(v));
            }
        }
        CHECK(peak < 0.2);
    }

    TEST_CASE("assign_mode fills cascade; DC gain near unity") {
        SuccessiveCompensatorCommissioner<3, float> cm{SuccessiveCompensatorConfig<float>{
            .Ts = 0.001f,
            .notch_Q = 5.0f,
        }};
        ResonantMode<float>                         m{};
        m.frequency_hz = 40.0f;
        m.omega_n = 2.0f * damp::numbers::pi_v<float> * 40.0f;
        m.zeta = 0.05f;
        m.peak_gain = 1.0f;
        m.valid = true;
        CHECK(cm.assign_mode(m));
        CHECK(cm.count() == 1);
        // RBJ notch DC gain is ~1 (small discrete warping).
        CHECK(cm.compensate(1.0f) == doctest::Approx(1.0f).epsilon(0.05));
    }
}
