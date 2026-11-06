// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>

#include "damp/estimation/rls.hpp"
#include "damp/estimation/state_variable_filter.hpp"
#include "damp/matrix/colvec.hpp"
#include "doctest.h"

using namespace damp;
using namespace damp::estimation;

namespace {
constexpr double pi = std::numbers::pi_v<double>;
}

TEST_SUITE("StateVariableFilter") {
    TEST_CASE("design succeeds and rejects nonsense") {
        static_assert(state_variable_filter(2.0 * pi * 20.0, 1e-3).success);
        CHECK_FALSE(state_variable_filter(0.0, 1e-3).success);
        CHECK_FALSE(state_variable_filter(10.0, 0.0).success);
    }

    TEST_CASE("its derivatives are the filtered derivatives, with no differencing") {
        // Filter sin(ωt) and, separately, its true derivative ω·cos(ωt): the first filter's
        // derivative outputs must match the second filter's value outputs.
        constexpr double               Ts = 1e-4;
        constexpr auto                 svf = state_variable_filter(2.0 * pi * 50.0, Ts);
        StateVariableFilter<2, double> f{svf};
        const double                   w = 2.0 * pi * 3.0;
        double                         worst_d = 0.0;
        double                         worst_dd = 0.0;
        for (int k = 0; k < 20000; ++k) {
            const double t = k * Ts;
            f.step(ColVec<2, double>{std::sin(w * t), w * std::cos(w * t)});
            if (t > 0.5) {
                worst_d = std::max(worst_d, std::abs(f[0].dx - f[1].x));
                worst_dd = std::max(worst_dd, std::abs(f[0].ddx - f[1].dx));
            }
        }
        CHECK(worst_d < 0.01 * w);      // derivative amplitude is w
        CHECK(worst_dd < 0.01 * w * w); // second derivative amplitude is w²
    }

    TEST_CASE("identifies a belt's stiffness and damping from quantized encoder counts") {
        // F = k δ + c δ̇, with δ read through 4.8 nm encoder counts and F with a little noise.
        // The regression runs in µm so the RLS prior (p₀ = 10⁶) does not outweigh the data.
        constexpr double               k_true = 2.0e4; // [N/m]
        constexpr double               c_true = 15.0;  // [N·s/m]
        constexpr double               um = 1e6;       // [µm/m]
        constexpr double               Ts = 1e-3;
        constexpr double               count = 4.8e-9; // [m]
        constexpr auto                 svf = state_variable_filter(2.0 * pi * 30.0, Ts);
        StateVariableFilter<2, double> f{svf};
        RlsVector<2, double>           rls{rls_vector(RlsConfig<double, 2>{.lambda = 1.0, .p0 = 1.0e6})};
        std::uint32_t                  seed = 7u;
        for (int k = 0; k < 20000; ++k) {
            const double t = k * Ts;
            // Stretch from a few overlapping moves: 1–10 µm, 1–9 Hz.
            const double delta = (5e-6 * std::sin(2.0 * pi * 1.3 * t)) + (3e-6 * std::sin(2.0 * pi * 4.1 * t)) + (1e-6 * std::sin(2.0 * pi * 8.7 * t));
            const double ddelta = (5e-6 * 2.0 * pi * 1.3 * std::cos(2.0 * pi * 1.3 * t)) + (3e-6 * 2.0 * pi * 4.1 * std::cos(2.0 * pi * 4.1 * t))
                                + (1e-6 * 2.0 * pi * 8.7 * std::cos(2.0 * pi * 8.7 * t));
            seed = (seed * 1664525u) + 1013904223u;
            const double noise = 0.005 * ((static_cast<double>(seed >> 8) / 8388608.0) - 1.0); // ±5 mN
            const double force = (k_true * delta) + (c_true * ddelta) + noise;
            const double delta_counts = std::round(delta / count) * count;
            f.step(ColVec<2, double>{delta_counts * um, force});
            if (k > 200) { // past the filter's start-up transient
                (void)rls.update(ColVec<2, double>{f[0].x, f[0].dx}, f[1].x);
            }
        }
        CHECK(rls.state().theta[0] * um == doctest::Approx(k_true).epsilon(0.01)); // [N/µm] → [N/m]
        CHECK(rls.state().theta[1] * um == doctest::Approx(c_true).epsilon(0.05));
    }

    TEST_CASE("identifies a mass and its damping from position counts alone (uses the second derivative)") {
        // F = m ẍ + b ẋ; position read through 50 nm counts, float like a target. In mm the
        // regressors are ~10²–10⁵ and the parameters ~10⁻³, so a prior p₀ = 1 is both loose enough
        // for the data and small enough that φᵀPφ stays inside float's range.
        constexpr double              m_true = 0.35; // [kg]
        constexpr double              b_true = 1.5;  // [N·s/m]
        constexpr double              Ts = 1e-3;
        constexpr double              count = 50e-9; // [m]
        constexpr double              mm = 1e3;
        StateVariableFilter<2, float> f{state_variable_filter(2.0 * pi * 30.0, Ts)};
        RlsVector<2, float>           rls{rls_vector(RlsConfig<float, 2>{.lambda = 1.0f, .p0 = 1.0f})};
        for (int k = 0; k < 20000; ++k) {
            const double t = k * Ts;
            const double w1 = 2.0 * pi * 2.0;
            const double w2 = 2.0 * pi * 7.0;
            const double x = (2e-3 * std::sin(w1 * t)) + (2e-4 * std::sin(w2 * t));
            const double v = (2e-3 * w1 * std::cos(w1 * t)) + (2e-4 * w2 * std::cos(w2 * t));
            const double a = -(2e-3 * w1 * w1 * std::sin(w1 * t)) - (2e-4 * w2 * w2 * std::sin(w2 * t));
            const double force = (m_true * a) + (b_true * v);
            const double x_counts = std::round(x / count) * count;
            f.step(ColVec<2, float>{static_cast<float>(x_counts * mm), static_cast<float>(force)});
            if (k > 200) {
                (void)rls.update(ColVec<2, float>{f[0].ddx, f[0].dx}, f[1].x);
            }
        }
        CHECK(static_cast<double>(rls.state().theta[0]) * mm == doctest::Approx(m_true).epsilon(0.02));
        CHECK(static_cast<double>(rls.state().theta[1]) * mm == doctest::Approx(b_true).epsilon(0.1));
    }
}
