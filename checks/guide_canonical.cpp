// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Keep in sync with the Canonical path snippet in docs/guide.md.

#include "damp/backend.hpp"
#include "damp/control.hpp"
using namespace damp;

namespace {
constexpr double Ts = 0.01; // [sec] Controller sample time

constexpr double m = 1.0; // kg
constexpr double k = 1.0; // N/m
constexpr double c = 0.5; // N·s/m

// x = [position, speed], u = force, y = position. G(s) = 1 / (m s² + c s + k)
constexpr auto sys = StateSpace<2, 1, 1>{
    .A = Matrix<2, 2>{{0.0, 1.0}, {-k / m, -c / m}},
    .B = Matrix<2, 1>{{0.0}, {1.0 / m}},
    .C = Matrix<1, 2>{{1.0, 0.0}},
};
} // namespace

// Configure tuning parameters (control weights) for an LQI
constexpr auto weights = design::LQRWeights<2, 1, 1>{
    .Q = Matrix<2, 2>{
        {10.0, 0.0},
        {0.0, 1.0},
    },
    .R = 1.0,
    .Q_integral = 80.0,
};

// Discretize the plant and solve the DARE
constexpr auto res = design::lqi(sys, Ts, weights);
static_assert(res.success);

// Convert the design result to a single-precision LQI controller instance
constinit LQI law{res.as<float>()};

// Fake function representing some PWM/DAC/comms
static void write_force(float /*u*/) {}

int main() {
    ColVec<2, float> x{}; // position, speed from sensors

    const float r = 1.0f; // Target 
    const float u = law.step(r, x)[0];         // Tick the controller and get the first output
    const float u_sat = clamp(u, -1.2f, 1.2f); // Saturate the control output
    law.back_calculate(u_sat);                 // Make the controller aware of the saturation

    // Apply the saturated control
    write_force(u_sat);
}
