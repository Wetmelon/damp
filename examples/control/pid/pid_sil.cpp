// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file pid_sil.cpp
 * @brief PI — thin host smoke (calls pid_controller.hpp)
 *
 * No ODE plant in this demo. Flashable path: pid_sketch.cpp.
 * For LTI SIL against the same PI API, see sim::simulate_lti_siso in tests.
 */

#include "damp/control.hpp"
#include "fmt/base.h"
#include "fmt/core.h"
#include "pid_controller.hpp"

using namespace damp;
using namespace damp::examples_pid;

static constinit PIController<float> controller = pi_design.as<float>();

int main() {
    fmt::print("===== PID / PI deploy (host smoke) =====\n\n");
    fmt::print("Ts={:.0f} Hz, Kp={:.1f}, Ki={:.1f} (unbounded)\n", 1.0f / Ts, 4.0, 8.0);
    fmt::print("First tick expect u ≈ Kp·e = 4 (y=0, r=1); Ki grows the command after that.\n\n");

    controller.reset();
    const float r = 1.0f;
    bool        ok = true;
    for (int i = 0; i < 10; ++i) {
        const float y = 0.0f; // mock process at rest
        const float u = control_period(controller, r, y);
        if (i == 0 || i == 9) {
            fmt::print("tick {}: u={:.4f}\n", i, u);
        }
        // e = 1 held: Kp·e plus the integral, Ki·Ts·e more each tick from the first
        const float expect = 4.0f + (8.0f * Ts * static_cast<float>(i + 1));
        if (damp::abs(u - expect) > 1e-4f) {
            fmt::print("FAIL: tick {} sent {:.6f}, not Kp·e + (k + 1)·Ki·Ts·e = {:.6f}\n", i, u, expect);
            ok = false;
        }
    }
    return ok ? 0 : 1;
}
