// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file lpf_sketch.cpp
 * @brief Low-pass filter — flashable sketch
 *
 * Host: examples toolchain. Target: copy + lpf_filter.hpp; replace mock ADC.
 */

#include "damp/filters/filters.hpp"
#include "lpf_filter.hpp"

using namespace damp;
using namespace damp::examples_lpf;

#if defined(ARDUINO)
#include <Arduino.h>
#else
static void pinMode(int, int) {}
static void delay(int) {}

static constexpr int OUTPUT = 1;
#define setup setup_
#define loop loop_
#endif

static constinit LowPass<1, float> lpf{coeffs1};

static float read_sample() {
    return 0.0f;
}

static void write_filtered(float /*y*/) {}

void setup() {
    pinMode(5, OUTPUT);
    lpf.reset();
}

void loop() {
    const float x = read_sample();
    const float y = filter_period(lpf, x);
    write_filtered(y);
    delay(static_cast<int>(Ts * 1000.0f));
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
