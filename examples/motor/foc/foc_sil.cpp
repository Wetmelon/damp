// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file foc_sil.cpp
 * @brief Host teaching plot: FOC current-loop PI vs I-P (calls foc_controller.hpp)
 *
 * Prefer motor/foc_switching for the product DiD path. This file is a host-only
 * contrast plot: electrical dq plant, no PWM edges.
 *
 * Plot: examples/plots/motor/foc_current_loop.html, foc_id_iq.html
 */

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "damp/backend.hpp"
#include "damp/math/transforms.hpp"
#include "damp/motor/foc.hpp"
#include "damp/motor/spm.hpp"
#include "damp/simulation/plot_plotly.hpp"
#include "fmt/base.h"
#include "fmt/core.h"
#include "foc_controller.hpp"
#include "plotlypp/figure.hpp"
#include "plotlypp/layout/layout.hpp"
#include "plotlypp/traces/scatter.hpp"

using namespace damp;
using namespace damp::motor;
using namespace damp::examples_foc;

namespace {

struct RunLog {
    std::vector<double> ts, iq_ref, iq, id, vmag; // sampled at the plant sub-tick rate
    double              rise_ms = 0.0;            // 98% rise after the step
    double              overshoot_pct = 0.0;      // peak iq overshoot in the step window
    double              peak_vmag = 0.0;          // peak |Vdq| in the step window
    double              final_iq = 0.0;           // iq at end (post-disturbance recovery)
};

// Run the full scenario with proportional setpoint weight b (1 = PI, 0 = I-P).
RunLog run(T b) {
    FOController<T> foc = {Ldq, Rs, lambda_pm, bw, b};
    foc.enable();

    DirectQuadrature<T> Idq{0.0, 0.0};
    RunLog              log;

    double     peak_iq = 0.0;
    bool       rise_found = false;
    const T    dt = Ts / substeps;
    const auto n_ticks = static_cast<int>(std::lround(t_end / Ts));

    for (int k = 0; k <= n_ticks; ++k) {
        const T                   t = k * Ts;
        const T                   iq_ref = (t >= t_step) ? Iq_ref : T{0};
        const DirectQuadrature<T> Idq_ref{T{0}, iq_ref};

        const auto cmd = foc.current_controller(Idq_ref, Idq, omega_e, Ts);
        const T    vd = cmd.Vdq.d;
        const T    vq = cmd.Vdq.q + ((t >= t_dist) ? vq_dist : T{0});
        const T    vmag = cmd.Vdq.abs();

        for (int s = 0; s < substeps; ++s) {
            const T did = (vd - (Rs * Idq.d) + (omega_e * Ldq.q * Idq.q)) / Ldq.d;
            const T diq = (vq - (Rs * Idq.q) - (omega_e * Ldq.d * Idq.d) - (omega_e * lambda_pm)) / Ldq.q;
            Idq.d += did * dt;
            Idq.q += diq * dt;

            log.ts.push_back((t + ((s + 1) * dt)) * 1e3); // [ms]
            log.iq_ref.push_back(iq_ref);
            log.id.push_back(Idq.d);
            log.iq.push_back(Idq.q);
            log.vmag.push_back(vmag);
        }

        // Step-response metrics, measured in the window before the disturbance.
        if (t >= t_step && t < t_dist) {
            if (!rise_found && Idq.q >= 0.98 * Iq_ref) {
                log.rise_ms = (t - t_step) * 1e3;
                rise_found = true;
            }
            peak_iq = std::max(peak_iq, Idq.q);
            log.peak_vmag = std::max(log.peak_vmag, vmag);
        }
    }

    log.overshoot_pct = 100.0 * (peak_iq - Iq_ref) / Iq_ref;
    log.final_iq = Idq.q;
    return log;
}

} // namespace

int main() {
    const T Kt = damp::motor::spm::torque_constant_from_flux(pole_pairs, lambda_pm);
    const T Km = damp::motor::spm::motor_constant(Kt, Rs);

    fmt::print("===== PMSM FOC current loop: PI vs I-P =====\n");
    fmt::print("Kt = {:.4f} Nm/A   Km = {:.3f} Nm/sqrt(W)   Vmax = {:.2f} V   f_bw = {:.0f} Hz   fsw = {:.0f} Hz\n\n", Kt, Km, Vmax, bw / (2.0 * numbers::pi_v<T>), fsw);

    const RunLog pi = run(T{1}); // standard PI (P on error)
    const RunLog ip = run(T{0}); // I-P        (P on measurement)

    fmt::print("  structure |  98%% rise [ms] | overshoot [%%] | peak |Vdq| [V] | final iq [A]\n");
    fmt::print("  ----------+----------------+---------------+----------------+-------------\n");
    fmt::print("  PI  (b=1) | {:14.3f} | {:13.1f} | {:14.3f} | {:11.3f}\n", pi.rise_ms, pi.overshoot_pct, pi.peak_vmag, pi.final_iq);
    fmt::print("  I-P (b=0) | {:14.3f} | {:13.1f} | {:14.3f} | {:11.3f}\n", ip.rise_ms, ip.overshoot_pct, ip.peak_vmag, ip.final_iq);
    fmt::print("\n  Vmax = {:.2f} V — PI peaks higher into the voltage circle; both reject the\n", Vmax);
    fmt::print("  {:.0f} V disturbance to the same steady iq (identical poles / rejection).\n", vq_dist);

    // ----- Plot: PI vs I-P overlay, currents + |Vdq| ------------------------
    using namespace plotlypp;
    const std::vector<double> vmax_line(pi.ts.size(), Vmax);

    auto event_lines = [&](const char* yref) {
        return std::vector<Layout::Shape>{
            Layout::Shape().type(Layout::Shape::Type::Line).x0(t_step * 1e3).x1(t_step * 1e3).xref("x").y0(0).y1(1).yref(yref).line(Layout::Shape::Line().dash("dot").color("green")),
            Layout::Shape().type(Layout::Shape::Type::Line).x0(t_dist * 1e3).x1(t_dist * 1e3).xref("x").y0(0).y1(1).yref(yref).line(Layout::Shape::Line().dash("dot").color("red"))
        };
    };

    Figure fig;
    // Panel 1: q-axis current — reference vs PI vs I-P (color sequence restarts per panel).
    fig.addTrace(Scatter().x(pi.ts).y(pi.iq_ref).mode({Scatter::Mode::Lines}).name("iq_ref").legend("legend").line([](auto& ln) { ln.color("gray").dash("dash"); }));
    fig.addTrace(Scatter().x(pi.ts).y(pi.iq).mode({Scatter::Mode::Lines}).name("iq — PI (b=1)").legend("legend").line([](auto& ln) { ln.color(plot::panel_color(0)); }));
    fig.addTrace(Scatter().x(ip.ts).y(ip.iq).mode({Scatter::Mode::Lines}).name("iq — I-P (b=0)").legend("legend").line([](auto& ln) { ln.color(plot::panel_color(1)); }));

    // Panel 2: voltage command magnitude vs the circle limit.
    fig.addTrace(Scatter().x(pi.ts).y(pi.vmag).mode({Scatter::Mode::Lines}).name("|Vdq| — PI").xaxis("x2").yaxis("y2").legend("legend2").line([](auto& ln) { ln.color(plot::panel_color(0)); }));
    fig.addTrace(Scatter().x(ip.ts).y(ip.vmag).mode({Scatter::Mode::Lines}).name("|Vdq| — I-P").xaxis("x2").yaxis("y2").legend("legend2").line([](auto& ln) { ln.color(plot::panel_color(1)); }));
    fig.addTrace(Scatter().x(pi.ts).y(vmax_line).mode({Scatter::Mode::Lines}).name("Vmax").xaxis("x2").yaxis("y2").legend("legend2").line([](auto& ln) { ln.color("black").dash("dot"); }));

    auto shapes = event_lines("y");
    for (auto& sh : event_lines("y2")) {
        shapes.push_back(sh);
    }

    const auto xd = plot::panel_x_domain();
    auto       layout = Layout()
                      .title([](auto& t) { t.text("PMSM FOC Current Loop — PI vs I-P (step + voltage disturbance)"); })
                      .xaxis(Layout::Xaxis().domain(xd).showticklabels(false).anchor("y"))
                      .yaxis(Layout::Yaxis().title([](auto& t) { t.text("q-axis current (A)"); }).domain({0.55, 1.0}))
                      .xaxis(2, Layout::Xaxis().domain(xd).matches("x").title([](auto& t) { t.text("Time (ms)"); }).anchor("y2"))
                      .yaxis(2, Layout::Yaxis().title([](auto& t) { t.text("|Vdq| (V)"); }).domain({0.0, 0.45}))
                      .legend(plot::panel_legend(1.0))
                      .legend(2, plot::panel_legend(0.45))
                      .shapes(shapes)
                      .height(800);
    fig.setLayout(damp::move(layout));
    plot::write_html(fig, "plots/motor/foc_current_loop.html");
    fmt::print("\n  Plot written to plots/motor/foc_current_loop.html\n");

    // id–iq loci: SPMSM holds id≈0, so the path is nearly a vertical line in iq.
    {
        double id_peak = 0.0;
        double iq_lo = 0.0;
        double iq_hi = 0.0;
        for (const auto* log : {&pi, &ip}) {
            for (double v : log->id) {
                id_peak = std::max(id_peak, std::abs(v));
            }
            for (double v : log->iq) {
                iq_lo = std::min(iq_lo, v);
                iq_hi = std::max(iq_hi, v);
            }
        }
        const double id_span = std::max(id_peak * 1.2, 0.5);
        const double iq_pad = std::max(0.08 * (iq_hi - iq_lo), 0.5);

        auto fig_xy = plot::plot_xy(
            {
                {.name = "PI (b=1)", .x = pi.id, .y = pi.iq, .lines = true, .markers = false, .color = std::string(plot::panel_color(0))},
                {.name = "I-P (b=0)", .x = ip.id, .y = ip.iq, .lines = true, .markers = false, .color = std::string(plot::panel_color(1))},
            },
            "PMSM FOC — id–iq loci (PI vs I-P; id≈0 for SPMSM Id*=0)",
            "id [A]",
            "iq [A]",
            plot::XyPlotOpts{
                .equal_aspect = false,
                .x_lo = -id_span,
                .x_hi = id_span,
                .y_lo = iq_lo - iq_pad,
                .y_hi = iq_hi + iq_pad,
            }
        );
        plot::write_html(fig_xy, "plots/motor/foc_id_iq.html");
        fmt::print("  Plot written to plots/motor/foc_id_iq.html\n");
    }
    return 0;
}
