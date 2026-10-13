// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file heater_oflqg_sil.cpp
 * @brief Heater with dead time — host SIL for the predictor-form OFLQG: the flashable tick on a
 *        hotend with the true transport delay, a colder room than assumed, and a fan at 240 s
 */

#include <cstddef>
#include <vector>

#include "damp/math/math.hpp"
#include "fmt/core.h"
#include "heater_oflqg_controller.hpp"

using namespace damp;
using namespace damp::examples_heater_oflqg;

namespace {

constexpr double kDt = 0.1;         // [s] plant integration step
constexpr int    kSub = 20;         // plant steps per control period
constexpr int    kLag = 50;         // plant steps of dead time (kDeadTime / kDt)
constexpr double kRoom = 22.0;      // [°C] colder than the controller assumes
constexpr double kSetpoint = 210.0; // [°C]
constexpr double kFanLoss = 5.0;    // [W] unmeasured, from kFanOn
constexpr double kFanOn = 240.0;    // [s]
constexpr double kEnd = 480.0;      // [s]

/// τ Ṫ = −(T − T_room) + K (P(t − θ) − loss), exact over each plant step.
struct Hotend {
    double              temperature{kRoom};
    std::vector<double> pipe = std::vector<double>(kLag, 0.0); ///< power on its way to the thermistor
    std::size_t         head{0};

    void step(double power, double loss) {
        const double arrived = pipe[head];
        pipe[head] = power;
        head = (head + 1) % pipe.size();
        const double target = kRoom + (kGain * (arrived - loss));
        const double a = damp::exp(-kDt / kTau);
        temperature = target + (a * (temperature - target));
    }
};

} // namespace

int main() {
    fmt::print("===== Hotend, 5 s dead time: offset-free LQG clamped at the actuator, two periods predicted (host SIL) =====\n\n");
    Heater    law{kLaw.as<float>()};
    Hotend    plant;
    double    overshoot = 0.0;
    double    worst_last_minute = 0.0;
    double    p_max = 0.0;
    const int periods = static_cast<int>(kEnd / kTs);
    for (int k = 0; k < periods; ++k) {
        const double t = k * kTs;
        const double p = static_cast<double>(heater_period(law, static_cast<float>(kSetpoint), static_cast<float>(plant.temperature)));
        p_max = damp::max(p_max, p);
        if (t < kFanOn) {
            overshoot = damp::max(overshoot, plant.temperature - kSetpoint);
        }
        if (t >= kEnd - 60.0) {
            worst_last_minute = damp::max(worst_last_minute, damp::abs(plant.temperature - kSetpoint));
        }
        for (int s = 0; s < kSub; ++s) {
            plant.step(p, (t >= kFanOn) ? kFanLoss : 0.0);
        }
    }
    fmt::print("heat-up overshoot {:.2f} °C   last-minute |e| under the fan {:.3f} °C   peak power {:.1f} W\n", overshoot, worst_last_minute, p_max);

    bool ok = true;
    if (!(overshoot < 0.5)) {
        fmt::print("FAIL: the heat-up overshoots (the delay is not being predicted)\n");
        ok = false;
    }
    if (!(worst_last_minute < 0.1)) {
        fmt::print("FAIL: an offset remains under the fan\n");
        ok = false;
    }
    if (!(p_max <= kPmax)) {
        fmt::print("FAIL: power above Pmax (the clamp)\n");
        ok = false;
    }
    return ok ? 0 : 1;
}
