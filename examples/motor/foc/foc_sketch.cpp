// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file foc_sketch.cpp
 * @brief FOC current-loop — thin smoke sketch (host teaching demo)
 *
 * Prefer foc_switching for product DiD. This sketch only smokes FOController::current_controller.
 * Full PI vs I-P plots: foc_sil.cpp.
 */

#include "damp/math/transforms.hpp"
#include "damp/motor/foc.hpp"
#include "foc_controller.hpp"

using namespace damp;
using namespace damp::examples_foc;
using namespace damp::motor;

#if defined(ARDUINO)
#include <Arduino.h>
#else
static void delayMicroseconds(int) {}
#define setup setup_
#define loop loop_
#endif

// Create Field Oriented (current) Controller
static constinit FOController<T> foc = {Ldq, Rs, lambda_pm, bw};

void setup() {
    foc.enable();
}

void loop() {
    const DirectQuadrature<T> Idq_ref{T{0}, Iq_ref};
    const DirectQuadrature<T> Idq{T{0}, T{0}};

    (void)foc.current_controller(Idq_ref, Idq, omega_e, Ts, Vmax);
    delayMicroseconds(static_cast<int>(Ts * 1e6));
}

#if !defined(ARDUINO)
int main() {
    setup();
    for (int i = 0; i < 10; ++i) {
        loop();
    }
    return 0;
}
#endif
