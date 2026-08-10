// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <numbers>
#include <type_traits>

#include "damp/filters/differentiator.hpp"
#include "damp/math/math.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {
constexpr double pi = 3.14159265358979323846;

// Quantize x to the nearest multiple of q (an encoder-count model).
double quantize(double x, double q) {
    return std::round(x / q) * q;
}
} // namespace

TEST_SUITE("differentiator") {
    TEST_CASE("recovers the derivative of a clean sinusoid") {
        // f = sin(w t), f_dot = w cos(w t), |f_ddot| <= w^2 = L.
        const double                      w = 2.0 * pi * 2.0; // 2 Hz
        const double                      dt = 1e-3;
        const double                      L = w * w;
        RobustExactDifferentiator<double> red(L, dt);

        double worst = 0.0;
        for (int k = 0; k < 4000; ++k) {
            const double t = k * dt;
            red.update(std::sin(w * t));
            if (k > 1000) { // after finite-time convergence
                worst = std::max(worst, std::abs(red.derivative() - (w * std::cos(w * t))));
            }
        }
        CHECK(worst < 0.05 * w); // derivative tracked to a few % of its amplitude
    }

    TEST_CASE("value() denoises the signal (z0 -> f)") {
        const double                      w = 2.0 * pi * 1.0;
        const double                      dt = 1e-3;
        RobustExactDifferentiator<double> red(w * w, dt);
        for (int k = 0; k < 3000; ++k) {
            red.update(std::sin(w * (k * dt)));
        }
        const double t = 3000 * dt;
        CHECK(red.value() == doctest::Approx(std::sin(w * t)).epsilon(0.02));
    }

    TEST_CASE("beats finite-difference on a quantized (encoder) signal") {
        // Slow move so quantization dominates: ~1 count/sample.
        const double w = 2.0 * pi * 0.5;    // 0.5 Hz
        const double A = 0.05;              // rad amplitude
        const double q = 2.0 * pi / 8192.0; // 2048-PPR quadrature encoder
        const double dt = 1e-3;
        const double L = A * w * w;

        RobustExactDifferentiator<double> red(5.0 * L, dt); // generous accel bound

        double sum_fd = 0.0;
        double sum_red = 0.0;
        double prev = quantize(0.0, q);
        int    n = 0;
        for (int k = 0; k < 6000; ++k) {
            const double t = k * dt;
            const double pos = quantize(A * std::sin(w * t), q);
            const double true_vel = A * w * std::cos(w * t);

            const double fd_vel = (pos - prev) / dt; // raw finite difference
            prev = pos;
            const double red_vel = red.update(pos);

            if (k > 2000) {
                sum_fd += (fd_vel - true_vel) * (fd_vel - true_vel);
                sum_red += (red_vel - true_vel) * (red_vel - true_vel);
                ++n;
            }
        }
        const double rms_fd = std::sqrt(sum_fd / n);
        const double rms_red = std::sqrt(sum_red / n);
        CHECK(rms_red < rms_fd);        // the differentiator is cleaner...
        CHECK(rms_red < 0.25 * rms_fd); // ...by a wide margin
    }

    TEST_CASE("invalid config is inert; reset re-seeds") {
        RobustExactDifferentiator<double> bad(-1.0, 1e-3);
        CHECK_FALSE(bad.valid());
        CHECK(bad.update(1.0) == doctest::Approx(0.0));

        RobustExactDifferentiator<double> red(10.0, 1e-3);
        red.reset(0.5, 2.0);
        CHECK(red.value() == doctest::Approx(0.5));
        CHECK(red.derivative() == doctest::Approx(2.0));
    }

    TEST_CASE("float specialization") {
        RobustExactDifferentiator<float> red(50.0f, 1e-3f);
        for (int k = 0; k < 2000; ++k) {
            red.update(std::sin(2.0f * std::numbers::pi_v<float> * 1.0f * (static_cast<float>(k) * 1e-3f)));
        }
        CHECK(std::isfinite(red.derivative()));
    }

    TEST_CASE("differentiator is constexpr-evaluable") {
        constexpr double err = []() consteval {
            const double                      w = 6.2831853;
            const double                      dt = 1e-3;
            RobustExactDifferentiator<double> red(w * w, dt);
            double                            last_t = 0.0;
            for (int k = 0; k < 2000; ++k) {
                last_t = k * dt;
                red.update(damp::sin(w * last_t));
            }
            return damp::abs(red.derivative() - (w * damp::cos(w * last_t)));
        }();
        static_assert(err < 1.0, "Levant differentiator must converge at compile time");
        CHECK(err < 1.0);
    }
}
