// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file encoder_velocity_sil.cpp
 * @brief Encoder / tach velocity estimation — host comparison + plots
 *
 * Uses encoder_velocity_estimator.hpp for PLL / Levant ticks. Compares raw FD,
 * LPF'd FD, PLL observer, and Levant differentiator on quantized sinusoids.
 */

#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

#include "damp/backend.hpp"
#include "damp/filters/differentiator.hpp"
#include "damp/simulation/plot_plotly.hpp"
#include "encoder_velocity_estimator.hpp"
#include "fmt/base.h"
#include "fmt/core.h"
#include "fmt/format.h"
#include "plotlypp/figure.hpp"
#include "plotlypp/layout/layout.hpp"
#include "plotlypp/traces/scatter.hpp"

using namespace damp;
using namespace damp::examples_encoder_velocity;

namespace {

constexpr double two_pi = 2.0 * std::numbers::pi;

struct Scenario {
    const char* label;
    double      counts_per_rev;
    double      bandwidth;   // PLL bandwidth [rad/s]
    double      amplitude;   // motion amplitude [turns]
    double      freq;        // motion frequency [Hz]
    double      lpf_tau;     // one-pole LPF time constant [s]
    double      zoom_half_s; // half-width of the position-zoom panel [s]
    const char* plot_file;
};

std::vector<double> decimate(const std::vector<double>& v, int d) {
    std::vector<double> out;
    for (size_t i = 0; i < v.size(); i += static_cast<size_t>(d)) {
        out.push_back(v[i]);
    }
    return out;
}

void run_scenario(const Scenario& s) {
    const double w = two_pi * s.freq;
    const double q = 1.0 / s.counts_per_rev;
    auto         true_pos = [&](double t) { return s.amplitude * std::sin(w * t); };
    auto         true_vel = [&](double t) { return s.amplitude * w * std::cos(w * t); };
    auto         encoder = [&](double t) { return std::round(true_pos(t) / q) * q; };

    const int steps = static_cast<int>(4.0 / s.freq / dt);
    const int settle = steps / 8;

    RobustExactDifferentiator<double> red(2.0 * s.amplitude * w * w, dt);
    const double                      lpf_alpha = dt / (s.lpf_tau + dt);
    double                            lpf_vel = 0.0;
    PllObserver                       pll(s.bandwidth);
    double                            pos_prev = encoder(0.0);

    std::vector<double> ts, p_true, p_enc, p_pll, p_red;
    std::vector<double> v_true, v_lpf, v_pll, v_red, e_lpf, e_pll, e_red;

    double sum_fd = 0, sum_lpf = 0, sum_pll = 0, sum_red = 0;
    double rev_fd = 0, rev_lpf = 0, rev_pll = 0, rev_red = 0;
    int    n = 0, n_rev = 0;

    for (int k = 0; k < steps; ++k) {
        const double t = k * dt;
        const double pos = encoder(t);
        const double vt = true_vel(t);

        const double fd = (pos - pos_prev) / dt;
        pos_prev = pos;
        lpf_vel += lpf_alpha * (fd - lpf_vel);
        const double vp = estimate_period_pll(pll, pos);
        const double rd = estimate_period_levant(red, pos);

        ts.push_back(t);
        p_true.push_back(true_pos(t));
        p_enc.push_back(pos);
        p_pll.push_back(pll.pos);
        p_red.push_back(red.value());
        v_true.push_back(vt);
        v_lpf.push_back(lpf_vel);
        v_pll.push_back(vp);
        v_red.push_back(rd);
        e_lpf.push_back(lpf_vel - vt);
        e_pll.push_back(vp - vt);
        e_red.push_back(rd - vt);

        if (k >= settle) {
            sum_fd += (fd - vt) * (fd - vt);
            sum_lpf += (lpf_vel - vt) * (lpf_vel - vt);
            sum_pll += (vp - vt) * (vp - vt);
            sum_red += (rd - vt) * (rd - vt);
            ++n;
            if (std::abs(vt) < 0.15 * s.amplitude * w) {
                rev_fd += (fd - vt) * (fd - vt);
                rev_lpf += (lpf_vel - vt) * (lpf_vel - vt);
                rev_pll += (vp - vt) * (vp - vt);
                rev_red += (rd - vt) * (rd - vt);
                ++n_rev;
            }
        }
    }

    auto rms = [](double s2, int cnt) { return std::sqrt(s2 / cnt); };
    fmt::print("--- {} ---\n", s.label);
    fmt::print(
        "  {:.0f} counts/rev ({:.2f} deg/count), bw={:.0f}, motion {:.2f} turns @ {:.2f} Hz, ~{:.2f} counts/sample peak\n",
        s.counts_per_rev,
        q * 360.0,
        s.bandwidth,
        s.amplitude,
        s.freq,
        (s.amplitude * w * dt) / q
    );
    fmt::print("  {:<26} {:>16} {:>18}\n", "estimator", "RMS err [turns/s]", "near-reversal");
    fmt::print("  {:<26} {:>16.4f} {:>18.4f}\n", "raw finite difference", rms(sum_fd, n), rms(rev_fd, n_rev));
    fmt::print("  {:<26} {:>16.4f} {:>18.4f}\n", "LPF'd diff", rms(sum_lpf, n), rms(rev_lpf, n_rev));
    fmt::print("  {:<26} {:>16.4f} {:>18.4f}\n", "PLL observer", rms(sum_pll, n), rms(rev_pll, n_rev));
    fmt::print("  {:<26} {:>16.4f} {:>18.4f}\n\n", "Levant differentiator", rms(sum_red, n), rms(rev_red, n_rev));

    using namespace plotlypp;
    const int rev_idx = static_cast<int>((1.0 / (4.0 * s.freq)) / dt);
    const int half = static_cast<int>(s.zoom_half_s / dt);
    const int z0 = rev_idx - half;
    const int z1 = rev_idx + half;
    auto      slice = [&](const std::vector<double>& v) {
        return std::vector<double>(v.begin() + z0, v.begin() + z1);
    };
    const auto tz = slice(ts);
    const int  D = 1 + steps / 2500;
    const auto td = decimate(ts, D);

    Figure fig;
    fig.addTrace(Scatter().x(tz).y(slice(p_true)).mode({Scatter::Mode::Lines}).name("true position").legend("legend").line([](auto& ln) {
        ln.color(plot::panel_color(0));
    }));
    fig.addTrace(Scatter().x(tz).y(slice(p_enc)).mode({Scatter::Mode::Lines}).name("encoder (quantized)").legend("legend").line([](auto& ln) {
        ln.color(plot::panel_color(1)).shape(Scatter::Line::Shape::Hv);
    }));
    fig.addTrace(Scatter().x(tz).y(slice(p_pll)).mode({Scatter::Mode::Lines}).name("PLL estimate").legend("legend").line([](auto& ln) {
        ln.color(plot::panel_color(2));
    }));
    fig.addTrace(Scatter().x(tz).y(slice(p_red)).mode({Scatter::Mode::Lines}).name("Levant value()").legend("legend").line([](auto& ln) {
        ln.color(plot::panel_color(3));
    }));

    fig.addTrace(Scatter().x(td).y(decimate(v_true, D)).mode({Scatter::Mode::Lines}).name("true velocity").xaxis("x2").yaxis("y2").legend("legend2").line([](auto& ln) {
        ln.color(plot::panel_color(0)).dash("dash");
    }));
    fig.addTrace(Scatter().x(td).y(decimate(v_lpf, D)).mode({Scatter::Mode::Lines}).name("LPF'd diff").xaxis("x2").yaxis("y2").legend("legend2").line([](auto& ln) {
        ln.color(plot::panel_color(1));
    }));
    fig.addTrace(Scatter().x(td).y(decimate(v_pll, D)).mode({Scatter::Mode::Lines}).name("PLL observer").xaxis("x2").yaxis("y2").legend("legend2").line([](auto& ln) {
        ln.color(plot::panel_color(2));
    }));
    fig.addTrace(Scatter().x(td).y(decimate(v_red, D)).mode({Scatter::Mode::Lines}).name("Levant").xaxis("x2").yaxis("y2").legend("legend2").line([](auto& ln) {
        ln.color(plot::panel_color(3));
    }));

    fig.addTrace(Scatter().x(td).y(decimate(e_lpf, D)).mode({Scatter::Mode::Lines}).name("LPF'd err").xaxis("x3").yaxis("y3").legend("legend3").line([](auto& ln) {
        ln.color(plot::panel_color(0));
    }));
    fig.addTrace(Scatter().x(td).y(decimate(e_pll, D)).mode({Scatter::Mode::Lines}).name("PLL err").xaxis("x3").yaxis("y3").legend("legend3").line([](auto& ln) {
        ln.color(plot::panel_color(1));
    }));
    fig.addTrace(Scatter().x(td).y(decimate(e_red, D)).mode({Scatter::Mode::Lines}).name("Levant err").xaxis("x3").yaxis("y3").legend("legend3").line([](auto& ln) {
        ln.color(plot::panel_color(2));
    }));

    const auto xd = plot::panel_x_domain();
    auto       layout = Layout()
                      .title([&](auto& t) { t.text(fmt::format("Encoder Velocity Estimation — {}", s.label)); })
                      .xaxis(Layout::Xaxis().domain(xd).title([](auto& t) { t.text("Time (s) — zoom on a reversal"); }).anchor("y"))
                      .yaxis(Layout::Yaxis().title([](auto& t) { t.text("Position (turns)"); }).domain({0.70, 1.0}))
                      .xaxis(2, Layout::Xaxis().domain(xd).showticklabels(false).anchor("y2"))
                      .yaxis(2, Layout::Yaxis().title([](auto& t) { t.text("Velocity (turns/s)"); }).domain({0.37, 0.63}))
                      .xaxis(3, Layout::Xaxis().domain(xd).matches("x2").title([](auto& t) { t.text("Time (s)"); }).anchor("y3"))
                      .yaxis(3, Layout::Yaxis().title([](auto& t) { t.text("Velocity error (turns/s)"); }).domain({0.0, 0.30}))
                      .legend(plot::panel_legend(1.0))
                      .legend(2, plot::panel_legend(0.63))
                      .legend(3, plot::panel_legend(0.30))
                      .height(1100);
    fig.setLayout(damp::move(layout));
    plot::write_html(fig, s.plot_file);
    fmt::print("  Plot written to {}\n\n", s.plot_file);
}

} // namespace

int main() {
    fmt::print("===== Servo Encoder Velocity Estimation =====\n");
    fmt::print("8 kHz loop. PLL = critically-damped tracking observer.\n\n");

    run_scenario({"14-bit absolute, bw=1000", 16384.0, 1000.0, 0.25, 1.0, 3.0e-3, 0.02, "plots/estimation/encoder_velocity_absolute.html"});
    run_scenario({"6-state hall, 7 pole pairs (42 cpr), bw=100", 42.0, 100.0, 1.0, 0.5, 8.0e-3, 0.15, "plots/estimation/encoder_velocity_hall.html"});

    fmt::print("Across both a fine 14-bit encoder and a coarse 42-count hall, the PLL\n");
    fmt::print("observer and the Levant differentiator stay neck-and-neck (within ~15%;\n");
    fmt::print("which one edges ahead flips with the regime), and both leave raw\n");
    fmt::print("differencing and the LPF well behind. The Levant block gets there\n");
    fmt::print("model-free, with no bandwidth to tune -- just an acceleration bound and\n");
    fmt::print("the standard gains -- and it's the matched rate estimator for the\n");
    fmt::print("super-twisting controller.\n");
    return 0;
}
