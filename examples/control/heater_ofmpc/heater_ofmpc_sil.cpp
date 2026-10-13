// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file heater_ofmpc_sil.cpp
 * @brief Heater with dead time — host SIL: the delay-aware MPC vs the same MPC without the delay
 *
 * Both laws run through heater_period, the tick the sketch flashes. A second run steps the
 * setpoint down from 210 to 190 °C at 360 s with the MPC told about it 60 s ahead (preview)
 * and not told (held setpoint). The plant is integrated
 * finely with its own transport delay; the ambient is 3 °C colder than the controller assumes,
 * and a part-cooling fan switches on halfway (an unmeasured 5 W heat loss).
 */

#include <cstddef>
#include <vector>

#include "damp/backend.hpp"
#include "damp/math/math.hpp"
#include "damp/simulation/plot_plotly.hpp"
#include "fmt/base.h"
#include "fmt/core.h"
#include "heater_ofmpc_controller.hpp"
#include "plotlypp/figure.hpp"
#include "plotlypp/layout/layout.hpp"
#include "plotlypp/traces/scatter.hpp"

using namespace damp;
using namespace damp::examples_heater_ofmpc;
using plotlypp::Figure;
using plotlypp::Layout;
using plotlypp::Scatter;

namespace {

constexpr double kDt = 0.1;           // [s] plant integration step
constexpr int    kSub = 20;           // plant steps per control period (kTs / kDt)
constexpr int    kPlotEvery = 2;      // plant steps per plotted sample: 0.2 s, ten per control period
constexpr int    kLag = 50;           // plant steps of dead time (kDeadTime / kDt)
constexpr double kAmbientTrue = 22.0; // [°C] the room, colder than kAmbient
constexpr double kSetpoint = 210.0;   // [°C]
constexpr double kFanLoss = 5.0;      // [W] part-cooling fan, unmeasured
constexpr double kFanOn = 240.0;      // [s]
constexpr double kEnd = 480.0;        // [s]
constexpr double kStepDownAt = 360.0; // [s] the profile run: setpoint 210 → 190 °C
constexpr int    kPeriods = static_cast<int>(kEnd / kTs);

/// The hotend: τ Ṫ = −(T − T_amb) + K (P(t − θ) − fan loss), exact over each plant step.
struct Hotend {
    double              temperature{kAmbientTrue};
    std::vector<double> pipe = std::vector<double>(kLag, 0.0); ///< heater power on its way to the thermistor
    std::size_t         head{0};

    void step(double power, double fan_loss) {
        const double delayed = pipe[head];
        pipe[head] = power;
        head = (head + 1) % pipe.size();
        const double target = kAmbientTrue + (kGain * (delayed - fan_loss));
        const double a = damp::exp(-kDt / kTau);
        temperature = target + (a * (temperature - target));
    }
};

struct Trace {
    std::vector<double> temperature; ///< one per control period, for the metrics
    std::vector<double> power;
    std::vector<double> t_fine; ///< every kPlotEvery plant steps, for the plot
    std::vector<double> temperature_fine;
    std::vector<double> power_fine; ///< the held power at those instants
    double              p_min{kPmax};
    double              p_max{0.0};
};

/// Close the loop with one law: the same tick the sketch runs, on its own copy of the plant.
template<typename Law>
Trace run(Law& law) {
    Hotend plant;
    Trace  tr;
    for (int k = 0; k < kPeriods; ++k) {
        const double t = k * kTs;
        const double p = static_cast<double>(heater_period(law, static_cast<float>(kSetpoint), static_cast<float>(plant.temperature)));
        tr.p_min = damp::min(tr.p_min, p);
        tr.p_max = damp::max(tr.p_max, p);
        tr.temperature.push_back(plant.temperature);
        tr.power.push_back(p);
        for (int s = 0; s < kSub; ++s) {
            if (s % kPlotEvery == 0) {
                tr.t_fine.push_back(t + (s * kDt));
                tr.temperature_fine.push_back(plant.temperature);
                tr.power_fine.push_back(p);
            }
            plant.step(p, (t >= kFanOn) ? kFanLoss : 0.0);
        }
    }
    return tr;
}

/// The profile run's setpoint at time t.
double profile_setpoint(double t) {
    return (t < kStepDownAt) ? kSetpoint : (kSetpoint - 20.0);
}

/// Close the loop on the setpoint profile, the MPC told about it over its horizon or not.
template<typename Mpc>
Trace run_profile(Mpc& mpc, bool preview) {
    Hotend plant;
    Trace  tr;
    for (int k = 0; k < kPeriods; ++k) {
        const double  t = k * kTs;
        SetpointAhead ahead{};
        for (std::size_t j = 0; j < kHorizon; ++j) {
            ahead[j] = static_cast<float>(preview ? profile_setpoint(t + ((static_cast<double>(j) + 1.0) * kTs)) : profile_setpoint(t));
        }
        const double p = static_cast<double>(heater_period(mpc, ahead, static_cast<float>(plant.temperature)));
        tr.p_min = damp::min(tr.p_min, p);
        tr.p_max = damp::max(tr.p_max, p);
        tr.temperature.push_back(plant.temperature);
        tr.power.push_back(p);
        for (int s = 0; s < kSub; ++s) {
            if (s % kPlotEvery == 0) {
                tr.t_fine.push_back(t + (s * kDt));
                tr.temperature_fine.push_back(plant.temperature);
                tr.power_fine.push_back(p);
            }
            plant.step(p, (t >= kFanOn) ? kFanLoss : 0.0);
        }
    }
    return tr;
}

/// ∫|T − setpoint| dt over [t0, t1) of the profile run [°C·s].
double profile_error(const Trace& tr, double t0, double t1) {
    double sum = 0.0;
    for (int k = static_cast<int>(t0 / kTs); k < static_cast<int>(t1 / kTs); ++k) {
        sum += damp::abs(tr.temperature[static_cast<std::size_t>(k)] - profile_setpoint(k * kTs)) * kTs;
    }
    return sum;
}

/// Largest excursion above the setpoint from period k0 on.
double overshoot(const Trace& tr, int k0, int k1) {
    double worst = 0.0;
    for (int k = k0; k < k1; ++k) {
        worst = damp::max(worst, tr.temperature[static_cast<std::size_t>(k)] - kSetpoint);
    }
    return worst;
}

/// Largest |error| over periods k0 … k1−1.
double worst_error(const Trace& tr, int k0, int k1) {
    double worst = 0.0;
    for (int k = k0; k < k1; ++k) {
        worst = damp::max(worst, damp::abs(tr.temperature[static_cast<std::size_t>(k)] - kSetpoint));
    }
    return worst;
}

void add_line(Figure& fig, const std::vector<double>& t, const std::vector<double>& y, const char* name, std::size_t color, int panel) {
    auto tr = Scatter().x(t).y(y).mode({Scatter::Mode::Lines}).name(name).legend(plot::panel_legend_id(panel));
    tr.line([color, panel](auto& ln) {
        ln.color(plot::panel_color(color));
        if (panel == 2) {
            ln.shape(Scatter::Line::Shape::Hv); // the power is held over each period
        }
    });
    if (panel == 2) {
        tr.xaxis("x2").yaxis("y2");
    }
    fig.addTrace(damp::move(tr));
}

void write_preview_plot(const Trace& with, const Trace& without) {
    const std::vector<double>& t = with.t_fine;
    std::vector<double>        sp;
    for (double ti : t) {
        sp.push_back(profile_setpoint(ti));
    }
    Figure fig;
    add_line(fig, t, sp, "setpoint profile", 2, 1);
    add_line(fig, t, with.temperature_fine, "T, setpoint known 60 s ahead", 0, 1);
    add_line(fig, t, without.temperature_fine, "T, setpoint held", 1, 1);
    add_line(fig, t, with.power_fine, "P, known ahead", 0, 2);
    add_line(fig, t, without.power_fine, "P, held", 1, 2);
    const auto xd = plot::panel_x_domain();
    auto       layout = Layout()
                      .title([](auto& tt) { tt.text("Hotend, 210 → 190 °C at 360 s: offset-free MPC with and without the setpoint known over its horizon"); })
                      .xaxis(Layout::Xaxis().domain(xd).showticklabels(false).anchor("y"))
                      .yaxis(Layout::Yaxis().title([](auto& tt) { tt.text("T [°C]"); }).domain({0.45, 1.0}))
                      .xaxis(2, Layout::Xaxis().domain(xd).matches("x").title([](auto& tt) { tt.text("Time (s)"); }).anchor("y2"))
                      .yaxis(2, Layout::Yaxis().title([](auto& tt) { tt.text("P [W]"); }).domain({0.0, 0.40}))
                      .legend(plot::panel_legend(1.0))
                      .legend(2, plot::panel_legend(0.40))
                      .height(800);
    fig.setLayout(damp::move(layout));
    plot::write_html(fig, "plots/control/heater_ofmpc_preview.html");
}

void write_plot(const Trace& aware, const Trace& naive) {
    const std::vector<double>& t = aware.t_fine;
    const std::vector<double>  sp(t.size(), kSetpoint);
    Figure                     fig;
    add_line(fig, t, sp, "setpoint", 2, 1);
    add_line(fig, t, aware.temperature_fine, "T, MPC on the delayed model", 0, 1);
    add_line(fig, t, naive.temperature_fine, "T, MPC ignoring the delay", 1, 1);
    add_line(fig, t, aware.power_fine, "P, delayed model", 0, 2);
    add_line(fig, t, naive.power_fine, "P, ignoring the delay", 1, 2);
    const auto xd = plot::panel_x_domain();
    auto       layout = Layout()
                      .title([](auto& tt) { tt.text("Hotend, 5 s dead time: offset-free MPC with and without the delay in its model"); })
                      .xaxis(Layout::Xaxis().domain(xd).showticklabels(false).anchor("y"))
                      .yaxis(Layout::Yaxis().title([](auto& tt) { tt.text("T [°C]"); }).domain({0.45, 1.0}))
                      .xaxis(2, Layout::Xaxis().domain(xd).matches("x").title([](auto& tt) { tt.text("Time (s)"); }).anchor("y2"))
                      .yaxis(2, Layout::Yaxis().title([](auto& tt) { tt.text("P [W]"); }).domain({0.0, 0.40}))
                      .legend(plot::panel_legend(1.0))
                      .legend(2, plot::panel_legend(0.40))
                      .height(800);
    fig.setLayout(damp::move(layout));
    plot::write_html(fig, "plots/control/heater_ofmpc.html");
}

} // namespace

int main() {
    fmt::print("===== Hotend, offset-free MPC through a 2.5-period dead time (host SIL) =====\n\n");
    constexpr auto naive_law = design_mpc<0.0>(); // the same design, blind to the dead time
    static_assert(naive_law.success);
    OffsetFreeMPC aware_mpc{kLaw.as<float>()};
    OffsetFreeMPC naive_mpc{naive_law.as<float>()};
    const Trace   aware = run(aware_mpc);
    const Trace   naive = run(naive_mpc);

    const int fan = static_cast<int>(kFanOn / kTs);
    const int last_minute = kPeriods - static_cast<int>(60.0 / kTs);
    fmt::print("                  heat-up overshoot   worst |e| after fan   last-minute |e|   power range\n");
    const auto row = [&](const char* name, const Trace& tr) {
        fmt::print(
            "{:<16}  {:14.2f} °C  {:17.2f} °C  {:13.3f} °C   {:4.1f} … {:4.1f} W\n",
            name,
            overshoot(tr, 0, fan),
            worst_error(tr, fan, kPeriods),
            worst_error(tr, last_minute, kPeriods),
            tr.p_min,
            tr.p_max
        );
    };
    row("delayed model", aware);
    row("ignoring delay", naive);
    write_plot(aware, naive);
    fmt::print("\nPlot written to plots/control/heater_ofmpc.html\n");

    // The profile run: the same delay-aware MPC, told about the step down or not.
    OffsetFreeMPC preview_mpc{kLaw.as<float>()};
    OffsetFreeMPC held_mpc{kLaw.as<float>()};
    const Trace   with = run_profile(preview_mpc, true);
    const Trace   without = run_profile(held_mpc, false);
    const double  err_with = profile_error(with, kStepDownAt - 60.0, kEnd);
    const double  err_without = profile_error(without, kStepDownAt - 60.0, kEnd);
    fmt::print("\nSetpoint 210 → 190 °C at 360 s, ∫|e| dt over 300 … 480 s:  known 60 s ahead {:6.1f} °C·s   held {:6.1f} °C·s\n", err_with, err_without);
    write_preview_plot(with, without);
    fmt::print("Plot written to plots/control/heater_ofmpc_preview.html\n");

    bool       ok = true;
    const auto check = [&ok](bool cond, const char* what) {
        if (!cond) {
            fmt::print("FAIL: {}\n", what);
            ok = false;
        }
    };
    check(worst_error(aware, last_minute, kPeriods) < 0.1, "the delay-aware MPC holds the setpoint without offset under the fan");
    check(aware.p_min >= 0.0 && aware.p_max <= kPmax, "heater power stays within 0 … Pmax");
    check(naive.p_min >= 0.0 && naive.p_max <= kPmax, "heater power stays within 0 … Pmax (no-delay law)");
    check(overshoot(naive, 0, fan) > (3.0 * overshoot(aware, 0, fan)) + 1.0, "ignoring the delay overshoots clearly more");
    check(err_with < 0.3 * err_without, "the setpoint known ahead cuts the error over the step down to under a third");
    check(profile_error(with, kEnd - 40.0, kEnd) / 40.0 < 0.5, "and lands on the new setpoint");
    return ok ? 0 : 1;
}
