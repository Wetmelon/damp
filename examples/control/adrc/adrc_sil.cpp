// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file adrc_sil.cpp
 * @brief ADRC vs well-tuned PI-D on a rigid inertia (host SIL + plot)
 *
 * Same deploy objects as adrc_sketch.cpp. Two copies of the plant: load step
 * then a 3× inertia jump — the PI-D is honest for the nameplate, not a strawman.
 */

#include <vector>

#include "adrc_controller.hpp"
#include "damp/backend.hpp"
#include "damp/controllers/adrc.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/math/math.hpp"
#include "damp/simulation/plot_plotly.hpp"
#include "damp/toolbox/conditioning.hpp"
#include "fmt/base.h"
#include "fmt/core.h"
#include "plotlypp/figure.hpp"
#include "plotlypp/layout/layout.hpp"
#include "plotlypp/traces/scatter.hpp"

using namespace damp;
using namespace damp::examples_adrc;
using plotlypp::Figure;
using plotlypp::Layout;
using plotlypp::Scatter;

namespace {

struct Axis {
    double theta{0.0};
    double omega{0.0};
    double J{J_nom};

    void step(double u, double tau_L, double dt) {
        const double fric = (B_visc * omega) + (tau_c * soft_sign(omega, omega_c));
        const double alpha = (u - fric - tau_L) / J;
        omega += alpha * dt;
        theta += omega * dt;
    }
};

[[nodiscard]] double reference_at(double t) {
    return (t < 3.6) ? 1.0 : 0.0;
}

[[nodiscard]] double load_at(double t) {
    return (t >= 1.2) ? 0.45 : 0.0;
}

[[nodiscard]] double inertia_at(double t) {
    return (t >= 2.4) ? (3.0 * J_nom) : J_nom;
}

} // namespace

int main() {
    fmt::print("===== ADRC vs PI-D (rigid inertia SIL) =====\n\n");
    fmt::print("J_nom={:.3f} kg·m²  wc={:.0f} rad/s  wo={:.0f} rad/s  b0={:.1f}\n", J_nom, wc, wo, b0);
    fmt::print("PI-D: Kp={:.3f}  Ki={:.3f}  Kd={:.3f}  (c=0, Tf=1/(10 wc))\n", pid_d.Kp, pid_d.Ki, pid_d.Kd);
    fmt::print("Scenario: step r=1 → load 0.45 N·m @ 1.2 s → J×3 @ 2.4 s → r=0 @ 3.6 s\n\n");

    ADRCController<2, double> adrc{adrc_d, Ts};
    PIDController<double>     pid{pid_d.discretize(Ts)};

    Axis plant_adrc{};
    Axis plant_pid{};

    const double T_end = 5.2;
    const int    N = static_cast<int>(T_end / Ts);

    std::vector<double> t, r, y_a, y_p, u_a, u_p, fhat;
    t.reserve(static_cast<std::size_t>(N));
    r.reserve(static_cast<std::size_t>(N));
    y_a.reserve(static_cast<std::size_t>(N));
    y_p.reserve(static_cast<std::size_t>(N));
    u_a.reserve(static_cast<std::size_t>(N));
    u_p.reserve(static_cast<std::size_t>(N));
    fhat.reserve(static_cast<std::size_t>(N));

    double e2_adrc = 0.0;
    double e2_pid = 0.0;
    int    n_post = 0;

    for (int k = 0; k < N; ++k) {
        const double tk = static_cast<double>(k) * Ts;
        const double rk = reference_at(tk);
        const double tau_L = load_at(tk);
        plant_adrc.J = inertia_at(tk);
        plant_pid.J = inertia_at(tk);

        double ua = control_period_adrc(adrc, rk, plant_adrc.theta);
        double up = control_period_pid(pid, rk, plant_pid.theta);
        const double ua_sat = damp::clamp(ua, -u_max, u_max);
        const double up_sat = damp::clamp(up, -u_max, u_max);
        adrc.back_calculate(ua, ua_sat);
        pid.back_calculate(up, up_sat);

        plant_adrc.step(ua_sat, tau_L, Ts);
        plant_pid.step(up_sat, tau_L, Ts);

        t.push_back(tk);
        r.push_back(rk);
        y_a.push_back(plant_adrc.theta);
        y_p.push_back(plant_pid.theta);
        u_a.push_back(ua_sat);
        u_p.push_back(up_sat);
        fhat.push_back(adrc.observer()[2]);

        if (tk >= 1.2) {
            e2_adrc += (plant_adrc.theta - rk) * (plant_adrc.theta - rk);
            e2_pid += (plant_pid.theta - rk) * (plant_pid.theta - rk);
            ++n_post;
        }
    }

    const double rms_a = damp::sqrt(e2_adrc / static_cast<double>(n_post));
    const double rms_p = damp::sqrt(e2_pid / static_cast<double>(n_post));
    fmt::print("RMS |e| after first load (t≥1.2 s):  ADRC={:.4f}   PI-D={:.4f}\n", rms_a, rms_p);
    fmt::print(
        "Final: ADRC θ={:.4f}   PI-D θ={:.4f}   (r=0)\n",
        plant_adrc.theta,
        plant_pid.theta
    );

    Figure fig;
    fig.addTrace(
        Scatter()
            .x(t)
            .y(r)
            .mode({Scatter::Mode::Lines})
            .name("r")
            .legend(plot::panel_legend_id(1))
            .line([](auto& ln) { ln.color(plot::panel_color(2)).dash("dash"); })
    );
    fig.addTrace(
        Scatter()
            .x(t)
            .y(y_a)
            .mode({Scatter::Mode::Lines})
            .name("θ ADRC")
            .legend(plot::panel_legend_id(1))
            .line([](auto& ln) { ln.color(plot::panel_color(0)); })
    );
    fig.addTrace(
        Scatter()
            .x(t)
            .y(y_p)
            .mode({Scatter::Mode::Lines})
            .name("θ PI-D")
            .legend(plot::panel_legend_id(1))
            .line([](auto& ln) { ln.color(plot::panel_color(1)); })
    );
    fig.addTrace(
        Scatter()
            .x(t)
            .y(u_a)
            .mode({Scatter::Mode::Lines})
            .name("u ADRC")
            .xaxis("x2")
            .yaxis("y2")
            .legend(plot::panel_legend_id(2))
            .line([](auto& ln) { ln.color(plot::panel_color(0)); })
    );
    fig.addTrace(
        Scatter()
            .x(t)
            .y(u_p)
            .mode({Scatter::Mode::Lines})
            .name("u PI-D")
            .xaxis("x2")
            .yaxis("y2")
            .legend(plot::panel_legend_id(2))
            .line([](auto& ln) { ln.color(plot::panel_color(1)); })
    );
    fig.addTrace(
        Scatter()
            .x(t)
            .y(fhat)
            .mode({Scatter::Mode::Lines})
            .name("f̂")
            .xaxis("x3")
            .yaxis("y3")
            .legend(plot::panel_legend_id(3))
            .line([](auto& ln) { ln.color(plot::panel_color(0)); })
    );

    const auto xd = plot::panel_x_domain();
    auto       layout = Layout()
                      .title([](auto& tt) { tt.text("ADRC vs PI-D — load step, then J×3"); })
                      .xaxis(Layout::Xaxis().domain(xd).showticklabels(false).anchor("y"))
                      .yaxis(Layout::Yaxis().title([](auto& tt) { tt.text("θ [rad]"); }).domain({0.70, 1.0}))
                      .xaxis(2, Layout::Xaxis().domain(xd).matches("x").showticklabels(false).anchor("y2"))
                      .yaxis(2, Layout::Yaxis().title([](auto& tt) { tt.text("u [N·m]"); }).domain({0.37, 0.63}))
                      .xaxis(3, Layout::Xaxis().domain(xd).matches("x").title([](auto& tt) { tt.text("Time (s)"); }).anchor("y3"))
                      .yaxis(3, Layout::Yaxis().title([](auto& tt) { tt.text("f̂ [rad/s²]"); }).domain({0.0, 0.30}))
                      .legend(plot::panel_legend(1.0))
                      .legend(2, plot::panel_legend(0.63))
                      .legend(3, plot::panel_legend(0.30))
                      .height(900);
    fig.setLayout(damp::move(layout));
    plot::write_html(fig, "plots/control/adrc_vs_pid.html");
    fmt::print("\nPlot written to plots/control/adrc_vs_pid.html\n");
    return 0;
}
