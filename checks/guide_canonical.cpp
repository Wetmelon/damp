// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Keep in sync with the Canonical path snippet in docs/guide.md.

#include "damp/control.hpp"
using namespace damp;

namespace {
constexpr double m = 1.0; // kg
constexpr double k = 1.0; // N/m
constexpr double c = 0.5; // N·s/m
// G(s) = 1 / (m s² + c s + k). Companion form is [position, speed] when m = 1.
constexpr auto G = TransferFunction<1, 3>{.num = {1.0}, .den = {m, c, k}};
constexpr auto sys = *G.to_state_space();
}

constexpr auto model = *discretize(sys, 0.01);
constexpr auto res = design::lqi(model, {.Q = Matrix<2, 2>{{10.0, 0.0}, {0.0, 1.0}}, .R = 1.0, .Q_integral = 80.0});
static_assert(res.success);
constinit LQI law{res.as<float>()};

static void write_force(float /*u*/) {} // PWM / DAC on the board

int main() {
    ColVec<2, float> x{}; // measured state from sensors
    const float      r = 1.0f;
    law.update(x);
    const float u = law.control(ColVec<1, float>{{r}})[0];
    const float u_sat = clamp(u, -1.2f, 1.2f);
    law.back_calculate(ColVec<1, float>{{u_sat}});
    law.predict();
    write_force(u_sat);
}
