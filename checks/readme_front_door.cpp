// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Keep in sync with the "Plant → law → tick" snippet in README.md. The ISR lines sit in main.

#include "damp/control.hpp"
using namespace damp;

constexpr double g = 9.81, L = 1.0, m = 1.0, b_damp = 0.1, Ts = 0.01;

constexpr auto sys = StateSpace{
    .A = Matrix<2, 2>{{0.0, 1.0}, {g / L, -b_damp / (m * L * L)}},
    .B = Matrix<2, 1>{{0.0}, {1.0 / (m * L * L)}},
    .C = Matrix<1, 2>{{1.0, 0.0}},
};
constexpr auto result = design::lqr(sys, Ts, {.Q = Matrix<2, 2>::identity() * 10.0, .R = 1.0});

constinit Law<result, float> pendulum; // the design's LQR, in float

int main() {
    // ISR: measurement x at this sample, command u
    ColVec<2, float> x{0.1f, 0.0f};
    ColVec<1, float> u = pendulum.step(x);
    (void)u;
}
