// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Keep in sync with the Canonical path snippet in docs/guide.md.

#include "damp/backend.hpp"
#include "damp/control.hpp"
using namespace damp;

// Fake function representing some PWM/DAC/comms
static void write_force(float /*u*/) {}

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
constexpr auto weights = design::LQRConfig<2, 1, 1>{
    .Q = Matrix<2, 2>{
        {10.0, 0.0},
        {0.0, 1.0},
    },
    .R = 1.0,
    .Q_integral = 80.0,
};

// Discretize the plant and solve the DARE
constexpr auto res = design::lqi(sys, Ts, weights);

// The design's LQI in single precision, its gains constants of the type; a failed design does not compile
constinit Law<res, float> law;

int main() {
    ColVec<2, float> x{}; // position, speed from sensors

    // Update the state (measurements)
    law.update(x);

    const float r = 1.0f;
    const float u = law.control(r)[0];         // Tick the controller and get the first output
    const float u_sat = clamp(u, -1.2f, 1.2f); // Saturate the control output
    law.back_calculate(u_sat);                 // Make the controller aware of the saturation
    law.predict();                             // Run the integrator (internally clamped)

    // Apply the saturated control
    write_force(u_sat);
}
