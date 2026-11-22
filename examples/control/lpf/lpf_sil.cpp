// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file lpf_sil.cpp
 * @brief Low-pass filter — host step-response demo (calls lpf_filter.hpp)
 */

#include "damp/control.hpp"
#include "fmt/base.h"
#include "fmt/core.h"
#include "lpf_filter.hpp"

using namespace damp;
using namespace damp::examples_lpf;

static constinit LowPass<1, float> lpf{coeffs1};
static constinit LowPass<2, float> lpf2{coeffs2};
static constinit LowPass<1, float> lpf1{fc, Ts}; // convenience ctor

int main() {
    fmt::print("===== LPF step response (host) =====\n\n");
    fmt::print("fc={:.0f} Hz, Ts={:.0f} us, first order (the sketch's filter)\n", fc, Ts * 1e6f);
    (void)lpf1;
    (void)lpf2;

    lpf.reset();

    bool  ok = true;
    float prev = 0.0f;
    for (int i = 0; i < 101; ++i) {
        const float output = filter_period(lpf, 1.0f);
        if (i % 10 == 0) {
            fmt::print("Step {}: Output = {:.4f}\n", i, output);
        }
        if (output < prev || output > 1.0f) {
            fmt::print("FAIL: step {} at {:.6f}: a first-order step response rises to 1 and never past it\n", i, output);
            ok = false;
        }
        prev = output;
    }
    if (!(prev > 0.995f)) {
        fmt::print("FAIL: {:.4f} after 6 time constants: not unit DC gain\n", prev);
        ok = false;
    }
    return ok ? 0 : 1;
}
