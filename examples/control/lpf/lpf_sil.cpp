// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file lpf_sil.cpp
 * @brief Low-pass filter — host step-response demo (calls lpf_filter.hpp)
 */

#include "damp/filters/filters.hpp"
#include "fmt/base.h"
#include "fmt/core.h"
#include "lpf_filter.hpp"

using namespace damp;
using namespace damp::examples_lpf;

static constinit LowPass<1, float> lpf{coeffs1};
static constinit LowPass<2, float> lpf2{{coeffs2.b0, coeffs2.b1, coeffs2.b2}, {coeffs2.a1, coeffs2.a2}};
static constinit LowPass<1, float> lpf1{fc, Ts}; // convenience ctor

int main() {
    fmt::print("===== LPF step response (host) =====\n\n");
    fmt::print("fc={:.0f} Hz, Ts={:.0f} us, 2nd-order zeta={:.3f}\n", fc, Ts * 1e6f, zeta);
    (void)lpf1;
    (void)lpf2;

    lpf.reset();

    for (int i = 0; i < 101; ++i) {
        const float output = filter_period(lpf, 1.0f);
        if (i % 10 == 0) {
            fmt::print("Step {}: Output = {:.4f}\n", i, output);
        }
    }

    return 0;
}
