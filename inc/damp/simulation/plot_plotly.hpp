// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file plot_plotly.hpp
 * @brief Plotly visualization helpers (host-only)
 *
 * @defgroup plot_plotly Plotly++ Visualization
 * @brief Host Plotly HTML figures for Bode, sims, and traces (plotlypp + json)
 *
 * Host-only. Builds `plotlypp::Figure` from damp analysis/simulation types.
 * Requires plotlypp and nlohmann-json on the include path. Intentionally not
 * exported by @c workbench.hpp so the host umbrella stays free of those
 * third-party includes; include this header directly when you need figures.
 * Never pull into @c control.hpp or target firmware.
 *
 * plotlypp is an unmodified third-party pin (git submodule). Do not patch
 * `libs/plotlypp`. Wet-specific HTML shell policy (full-window layout, no
 * phantom horizontal scrollbar, skip re-emitting plotly.min.js) lives in
 * @ref write_html — prefer that over `Figure::writeHtml`.
 *
 * Usage:
 * @code
 *   auto sim = simulate(...);
 *   auto fig = plot::plot_simulation(sim, "Pendulum Response");
 *   plot::write_html(fig, "plots/pendulum.html");
 * @endcode
 */

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <limits>
#include <plotlypp/figure.hpp>
#include <plotlypp/layout/layout.hpp>
#include <plotlypp/traces/scatter.hpp>
#include <sstream>
#include <string>
#include <vector>

#include "damp/analysis/analysis.hpp"
#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/simulation/simulate.hpp"

namespace damp {
namespace plot {

/**
 * @brief Apply damp's full-window HTML shell to plotlypp-generated markup
 *
 * plotlypp's default shell uses `100vw` and re-`newPlot`s on resize. We rewrite
 * the shell so figures always fill the browser viewport (width and height),
 * ignoring any fixed layout.width / layout.height from the C++ side — those
 * only shrink or pin the plot in Edge/Chrome after a hard refresh on a small
 * window. Submodule stays unmodified.
 */
[[nodiscard]] inline std::string apply_wet_plotly_shell(std::string html) {
    // Avoid 100vw (includes scrollbar gutter → phantom horizontal scrollbar).
    const std::string old_css = "#plot { width: 100vw; height: 100vh; }";
    const std::string new_css =
        "html, body { margin: 0; width: 100%; height: 100%; overflow: hidden; }\n"
        "                    #plot { width: 100%; height: 100%; }";
    if (const auto p = html.find(old_css); p != std::string::npos) {
        html.replace(p, old_css.size(), new_css);
    }

    // Always size to the live viewport — do not floor on layout.height.
    const std::string marker = "const plotDiv = document.getElementById('plot');";
    const auto        mpos = html.find(marker);
    if (mpos == std::string::npos) {
        return html;
    }
    const auto script_end = html.find("</script>", mpos);
    if (script_end == std::string::npos) {
        return html;
    }

    const std::string wet_js =
        "const plotDiv = document.getElementById('plot');\n"
        "                    const fitViewport = () => {\n"
        "                        layout.autosize = true;\n"
        "                        delete layout.width;\n"
        "                        layout.height = window.innerHeight;\n"
        "                    };\n"
        "                    fitViewport();\n"
        "                    const config = {responsive: true};\n"
        "                    Plotly.newPlot(plotDiv, data, layout, config);\n"
        "                    let resizing = false;\n"
        "                    window.addEventListener('resize', () => {\n"
        "                        if (!resizing) {\n"
        "                            resizing = true;\n"
        "                            window.requestAnimationFrame(() => {\n"
        "                                fitViewport();\n"
        "                                Plotly.relayout(plotDiv, {height: layout.height, autosize: true});\n"
        "                                Plotly.Plots.resize(plotDiv);\n"
        "                                resizing = false;\n"
        "                            });\n"
        "                        }\n"
        "                    });\n"
        "                ";
    html.replace(mpos, script_end - mpos, wet_js);
    return html;
}

/**
 * @brief Walk up from @p dir until a directory named "plots" is found
 *
 * Used so HTML under plots/power/, plots/motor/, … share one plots/js/plotly.min.js
 * instead of a multi‑MB copy beside every domain folder.
 */
[[nodiscard]] inline std::filesystem::path plots_root_for(std::filesystem::path dir) {
    namespace fs = std::filesystem;
    if (dir.empty()) {
        return fs::path{"plots"};
    }
    fs::path cur = dir;
    for (int i = 0; i < 8; ++i) {
        if (cur.filename() == "plots") {
            return cur;
        }
        const fs::path parent = cur.parent_path();
        if (parent.empty() || parent == cur) {
            break;
        }
        cur = parent;
    }
    return dir; // not under plots/ — keep assets next to the HTML
}

/**
 * @brief Write a plotlypp figure to HTML with damp's shell (prefer over writeHtml)
 *
 * @param fig                  Figure to write
 * @param path                 Output .html path (e.g. plots/power/buck_sil.html)
 * @param include_js_resources If true, emit shared plotly.min.js only when missing
 *                             under the plots/ tree (or next to the HTML if outside)
 */
inline void write_html(
    const plotlypp::Figure&      fig,
    const std::filesystem::path& path,
    bool                         include_js_resources = true
) {
    namespace fs = std::filesystem;
    const fs::path html_dir = path.parent_path().empty() ? fs::path{"."} : path.parent_path();
    fs::create_directories(html_dir);

    const fs::path plots_root = plots_root_for(html_dir);
    const fs::path js_path = plots_root / "js" / "plotly.min.js";
    const bool     need_shared_js = include_js_resources && !fs::exists(js_path);

    // plotlypp always writes js/ next to the HTML. Seed the shared plots/js once
    // via a short-lived sibling HTML, then write the real figure without JS.
    if (need_shared_js) {
        fs::create_directories(plots_root / "js");
        const fs::path seed = plots_root / "_wet_plotly_js_seed.html";
        fig.writeHtml(seed, true);
        std::error_code ec;
        fs::remove(seed, ec);
        // plotlypp may have written plots/js/plotly.min.js already when seed was
        // at plots/_seed.html; if seed lived under plots/power/, move js up.
        const fs::path local_js = html_dir / "js" / "plotly.min.js";
        if (!fs::exists(js_path) && fs::exists(local_js)) {
            fs::create_directories(js_path.parent_path());
            fs::rename(local_js, js_path, ec);
            fs::remove(html_dir / "js", ec);
        }
    }

    fig.writeHtml(path, false);

    // Read via rdbuf: GCC 14 -Wnull-dereference false
    // positives inside libstdc++ when the iterator form is inlined here.
    std::ifstream      in(path, std::ios::binary);
    std::ostringstream oss;
    oss << in.rdbuf();
    std::string html = std::move(oss).str();
    in.close();

    // Point script src at the shared js (e.g. ../js/plotly.min.js from plots/power/)
    if (include_js_resources) {
        std::error_code ec;
        fs::path        rel = fs::relative(js_path, html_dir, ec);
        if (ec || rel.empty()) {
            rel = fs::path{"js"} / "plotly.min.js";
        }
        std::string       rel_s = rel.generic_string();
        const std::string old_src = "js/plotly.min.js";
        for (std::size_t p = 0; (p = html.find(old_src, p)) != std::string::npos;) {
            // only replace path-like occurrences (src="js/plotly.min.js")
            html.replace(p, old_src.size(), rel_s);
            p += rel_s.size();
        }
    }

    html = apply_wet_plotly_shell(std::move(html));
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    out << html;
}

/**
 * @brief Convert a ColVec<N,T> to std::vector<double> for plotlypp
 */
template<size_t N, typename T>
std::vector<double> to_std_vector(const ColVec<N, T>& v) {
    std::vector<double> out(N);
    for (size_t i = 0; i < N; ++i) {
        out[i] = static_cast<double>(v(i, 0));
    }
    return out;
}

/**
 * @brief Extract the i-th element from each vector entry into a std::vector<double>
 */
template<size_t N, typename T>
std::vector<double> extract_channel(const std::vector<ColVec<N, T>>& history, size_t channel) {
    std::vector<double> out;
    out.reserve(history.size());
    for (const auto& v : history) {
        out.push_back(static_cast<double>(v(channel, 0)));
    }
    return out;
}

/**
 * @brief Convert std::vector\<T\> to std::vector<double>
 */
template<typename T>
std::vector<double> to_double_vector(const std::vector<T>& v) {
    std::vector<double> out(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        out[i] = static_cast<double>(v[i]);
    }
    return out;
}

// ---------------------------------------------------------------------------
// Default multi-panel format (stacked shared-x figures)
//
// This is the house style for examples and host tooling. Prefer these helpers
// over hand-rolled Layout::Legend / colorway picks.
//
//   1. Linked x: secondary axes use .matches("x") (do not match axes that
//      intentionally use different domains, e.g. a zoom inset).
//   2. Per-panel legends: traces → panel_legend_id(row); layout →
//      .legend(i, panel_legend(y_top_of_band)).
//   3. Per-panel colors: panel_color(i) restarts at 0 each subplot (Plotly's
//      default colorway is global across the figure otherwise).
//   4. X domain: .domain(panel_x_domain()) on every stacked x-axis. Default
//      panel_x_right = 1.0 (full width); lower it only if a figure needs a
//      wider right gutter for long legend labels.
//   5. Legend x = panel_x_right + 0.01 (just past the plot column).
//
// HTML shell: use plot::write_html so damp's full-window
// responsive shell applies without patching the plotlypp submodule.
// ---------------------------------------------------------------------------

/// Right edge of the plot column in paper coords (default: full width).
inline constexpr double panel_x_right = 1.0;

/**
 * @brief Horizontal domain for stacked multi-panel x-axes
 *
 * Pass to every x-axis: `.domain(plot::panel_x_domain())`.
 */
[[nodiscard]] inline std::vector<double> panel_x_domain() {
    return {0.0, panel_x_right};
}

/**
 * @brief Legend id for subplot row (1-based): "legend", "legend2", …
 *
 * Pass to Scatter::legend(...). Layout legends are configured with
 * Layout::legend(index, ...) where index 1 is the default "legend".
 */
[[nodiscard]] inline std::string panel_legend_id(int row) {
    if (row <= 1) {
        return "legend";
    }
    return "legend" + std::to_string(row);
}

/**
 * @brief Per-panel legend box just right of a stacked subplot band
 *
 * Pair with @ref panel_x_domain on every x-axis (see multi-panel format above).
 *
 * @param y_top  Paper y of the top of the panel's domain (e.g. 1.0, 0.63, 0.30)
 */
[[nodiscard]] inline plotlypp::Layout::Legend panel_legend(double y_top) {
    using Lg = plotlypp::Layout::Legend;
    return Lg()
        .x(panel_x_right + 0.01)
        .y(y_top)
        .xanchor(Lg::Xanchor::Left)
        .yanchor(Lg::Yanchor::Top);
}

/**
 * @brief Plotly default colorway entry for index i within a single subplot
 *
 * Restart the sequence at 0 for each panel so the first series is always the
 * same blue, not the next color after the previous panel's last trace.
 *
 * @param index 0-based series index *within the panel* (wraps the palette)
 */
[[nodiscard]] inline const char* panel_color(std::size_t index) {
    // Plotly.js default qualitative colorway (layout.colorway).
    static constexpr const char* kColors[] = {
        "#636efa",
        "#EF553B",
        "#00cc96",
        "#ab63fa",
        "#FFA15A",
        "#19d3f3",
        "#FF6692",
        "#B6E880",
        "#FF97FF",
        "#FECB52",
    };
    constexpr std::size_t n = sizeof(kColors) / sizeof(kColors[0]);
    return kColors[index % n];
}

/**
 * @brief Plot simulation results with subplots for states, outputs, and inputs
 *
 * Creates a 3-row subplot: states on top, outputs in middle, control inputs on
 * bottom. X-axes match (linked zoom); each row has its own legend.
 *
 * @param sim   SimulationResult from simulate() or simulate_state_feedback()
 * @param title Plot title
 * @return plotlypp::Figure ready for .show() or .writeHtml()
 */
template<size_t NX, size_t NU, size_t NY, typename T>
plotlypp::Figure plot_simulation(
    const sim::SimulationResult<NX, NU, NY, T>& sim,
    const std::string&                          title = "Simulation"
) {
    using namespace plotlypp;

    auto   t = to_double_vector(sim.t);
    Figure fig;

    // States subplot (row 1) — color index restarts per panel
    for (size_t i = 0; i < NX; ++i) {
        auto trace = Scatter()
                         .x(t)
                         .y(extract_channel(sim.x, i))
                         .mode({Scatter::Mode::Lines})
                         .name("x" + std::to_string(i))
                         .legend(panel_legend_id(1))
                         .line([i](auto& ln) { ln.color(panel_color(i)); });
        fig.addTrace(damp::move(trace));
    }

    // Outputs subplot (row 2)
    for (size_t i = 0; i < NY; ++i) {
        auto trace = Scatter()
                         .x(t)
                         .y(extract_channel(sim.y, i))
                         .mode({Scatter::Mode::Lines})
                         .name("y" + std::to_string(i))
                         .xaxis("x2")
                         .yaxis("y2")
                         .legend(panel_legend_id(2))
                         .line([i](auto& ln) { ln.color(panel_color(i)); });
        fig.addTrace(damp::move(trace));
    }

    // Inputs subplot (row 3)
    for (size_t i = 0; i < NU; ++i) {
        auto trace = Scatter()
                         .x(t)
                         .y(extract_channel(sim.u, i))
                         .mode({Scatter::Mode::Lines})
                         .name("u" + std::to_string(i))
                         .xaxis("x3")
                         .yaxis("y3")
                         .legend(panel_legend_id(3))
                         .line([i](auto& ln) { ln.color(panel_color(i)); });
        fig.addTrace(damp::move(trace));
    }

    // Stacked domains (top → bottom) with gaps; legends sit at each band top.
    const auto xd = panel_x_domain();
    auto       layout = Layout()
                      .title([&](auto& tt) { tt.text(title); })
                      .xaxis(Layout::Xaxis().domain(xd).showticklabels(false).anchor("y"))
                      .yaxis(Layout::Yaxis().title([](auto& tt) { tt.text("States"); }).domain({0.70, 1.0}))
                      .xaxis(2, Layout::Xaxis().domain(xd).matches("x").showticklabels(false).anchor("y2"))
                      .yaxis(2, Layout::Yaxis().title([](auto& tt) { tt.text("Outputs"); }).domain({0.37, 0.63}))
                      .xaxis(3, Layout::Xaxis().domain(xd).matches("x").title([](auto& tt) { tt.text("Time (s)"); }).anchor("y3"))
                      .yaxis(3, Layout::Yaxis().title([](auto& tt) { tt.text("Inputs"); }).domain({0.0, 0.30}))
                      .legend(panel_legend(1.0))
                      .legend(2, panel_legend(0.63))
                      .legend(3, panel_legend(0.30))
                      .height(900);

    fig.setLayout(damp::move(layout));
    return fig;
}

/**
 * @brief Plot Bode magnitude and phase as subplots
 *
 * Magnitude on top, phase on bottom; log-frequency x-axes match; each panel
 * has its own legend.
 *
 * @param bode  BodeResult from analysis::bode()
 * @param title Plot title
 * @return plotlypp::Figure
 */
template<typename T>
plotlypp::Figure plot_bode(
    const analysis::BodeResult<T>& bode,
    const std::string&             title = "Bode Plot"
) {
    using namespace plotlypp;

    std::vector<double> omega, mag_db, phase_deg;
    omega.reserve(bode.points.size());
    mag_db.reserve(bode.points.size());
    phase_deg.reserve(bode.points.size());

    for (const auto& pt : bode.points) {
        omega.push_back(static_cast<double>(pt.omega));
        mag_db.push_back(static_cast<double>(pt.magnitude_db));
        phase_deg.push_back(static_cast<double>(pt.phase_deg));
    }

    auto mag_trace = Scatter()
                         .x(omega)
                         .y(mag_db)
                         .mode({Scatter::Mode::Lines})
                         .name("Magnitude")
                         .legend(panel_legend_id(1))
                         .line([](auto& ln) { ln.color(panel_color(0)); });

    auto phase_trace = Scatter()
                           .x(omega)
                           .y(phase_deg)
                           .mode({Scatter::Mode::Lines})
                           .name("Phase")
                           .xaxis("x2")
                           .yaxis("y2")
                           .legend(panel_legend_id(2))
                           .line([](auto& ln) { ln.color(panel_color(0)); });

    const auto xd = panel_x_domain();
    auto       layout = Layout()
                      .title([&](auto& tt) { tt.text(title); })
                      .xaxis(Layout::Xaxis().domain(xd).type(Layout::Xaxis::Type::Log).showticklabels(false).anchor("y"))
                      .yaxis(Layout::Yaxis().title([](auto& tt) { tt.text("Magnitude (dB)"); }).domain({0.55, 1.0}))
                      .xaxis(
                          2,
                          Layout::Xaxis()
                              .domain(xd)
                              .type(Layout::Xaxis::Type::Log)
                              .matches("x")
                              .title([](auto& tt) { tt.text("Frequency (rad/s)"); })
                              .anchor("y2")
                      )
                      .yaxis(2, Layout::Yaxis().title([](auto& tt) { tt.text("Phase (deg)"); }).domain({0.0, 0.45}))
                      .legend(panel_legend(1.0))
                      .legend(2, panel_legend(0.45))
                      .height(700);

    return Figure()
        .addTrace(damp::move(mag_trace))
        .addTrace(damp::move(phase_trace))
        .setLayout(damp::move(layout));
}

/**
 * @brief Simple line plot of time vs value
 *
 * @param time   Time vector
 * @param values Value vector
 * @param title  Plot title
 * @return plotlypp::Figure
 */
template<typename T>
plotlypp::Figure plot_line(
    const std::vector<T>& time,
    const std::vector<T>& values,
    const std::string&    title = "Plot"
) {
    using namespace plotlypp;

    auto trace = Scatter()
                     .x(to_double_vector(time))
                     .y(to_double_vector(values))
                     .mode({Scatter::Mode::Lines});

    auto layout = Layout()
                      .title([&](auto& t) { t.text(title); })
                      .xaxis(Layout::Xaxis().title([](auto& t) { t.text("Time (s)"); }))
                      .yaxis(Layout::Yaxis().title([](auto& t) { t.text("Value"); }));

    return Figure().addTrace(damp::move(trace)).setLayout(damp::move(layout));
}

/**
 * @brief One series for a planar (x, y) scatter / path plot
 *
 * Prefer this over two independent 1-D time series when the story is a path,
 * locus, or phase portrait (toolhead XY, motor A/B space, id–iq MTPA, …).
 */
struct XySeries {
    std::string         name;
    std::vector<double> x;
    std::vector<double> y;
    bool                lines{true};    ///< connect points in sample order
    bool                markers{false}; ///< mark samples (paths often lines-only)
    std::string         color{};        ///< empty → plotly default colorway
};

/**
 * @brief Optional axis framing for @ref plot_xy
 *
 * `equal_aspect` locks pixel scale (good for workspace circles / MTPA loci).
 * It also forces linked zoom — turn it off when free box-zoom is more useful
 * (e.g. an id–iq step locus that is nearly a vertical line at id ≈ 0).
 *
 * Finite `x_lo`/`x_hi` or `y_lo`/`y_hi` pin that axis; leave non-finite for auto.
 */
struct XyPlotOpts {
    bool   equal_aspect{false};
    double x_lo{std::numeric_limits<double>::quiet_NaN()};
    double x_hi{std::numeric_limits<double>::quiet_NaN()};
    double y_lo{std::numeric_limits<double>::quiet_NaN()};
    double y_hi{std::numeric_limits<double>::quiet_NaN()};
};

/**
 * @brief Multi-series planar scatter (true 2-D path / locus)
 *
 * Unlike @ref plot_line (one signal vs time), both axes are data. Use for
 * workspace paths, actuator-space paths, and phase portraits.
 *
 * @param series  One or more (x, y) series of equal length per series
 * @param title   Figure title
 * @param xlabel  X-axis label
 * @param ylabel  Y-axis label
 * @param opts    Aspect lock and optional fixed axis ranges
 * @return plotlypp::Figure
 */
inline plotlypp::Figure plot_xy(
    const std::vector<XySeries>& series,
    const std::string&           title = "XY",
    const std::string&           xlabel = "x",
    const std::string&           ylabel = "y",
    XyPlotOpts                   opts = {}
) {
    using namespace plotlypp;
    Figure fig;
    for (std::size_t i = 0; i < series.size(); ++i) {
        const auto& s = series[i];
        const bool  use_lines = s.lines;
        const bool  use_markers = s.markers || !s.lines;
        Scatter     tr;
        tr.x(s.x).y(s.y).name(s.name).legend("legend");
        if (use_lines && use_markers) {
            tr.mode({Scatter::Mode::Lines, Scatter::Mode::Markers});
        } else if (use_markers) {
            tr.mode({Scatter::Mode::Markers});
        } else {
            tr.mode({Scatter::Mode::Lines});
        }
        const std::string col = s.color.empty() ? std::string(panel_color(i)) : s.color;
        tr.line([&col](auto& ln) { ln.color(col).width(2.0); });
        if (use_markers) {
            tr.marker([&col](auto& mk) { mk.color(col); });
        }
        fig.addTrace(damp::move(tr));
    }

    auto xax = Layout::Xaxis().title([&](auto& t) { t.text(xlabel); });
    auto yax = Layout::Yaxis().title([&](auto& t) { t.text(ylabel); }).zeroline(true);
    xax.zeroline(true);
    if (opts.equal_aspect) {
        // Locks scale *and* couples zoom; use only when a circle must stay round.
        yax.scaleanchor("x").scaleratio(1.0);
    }
    if (damp::isfinite(opts.x_lo) && damp::isfinite(opts.x_hi) && opts.x_hi > opts.x_lo) {
        xax.range({opts.x_lo, opts.x_hi});
    }
    if (damp::isfinite(opts.y_lo) && damp::isfinite(opts.y_hi) && opts.y_hi > opts.y_lo) {
        yax.range({opts.y_lo, opts.y_hi});
    }
    fig.setLayout(Layout().title([&](auto& t) { t.text(title); }).xaxis(damp::move(xax)).yaxis(damp::move(yax)));
    return fig;
}

/// Convenience overload — equal-aspect flag only.
inline plotlypp::Figure plot_xy(
    const std::vector<XySeries>& series,
    const std::string&           title,
    const std::string&           xlabel,
    const std::string&           ylabel,
    bool                         equal_aspect
) {
    return plot_xy(series, title, xlabel, ylabel, XyPlotOpts{.equal_aspect = equal_aspect});
}

/**
 * @brief Plot step response data
 *
 * @param time_values Pair of {time_vector, response_vector} from step_response()
 * @param title Plot title
 */
template<typename T>
plotlypp::Figure plot_step(
    const damp::pair<std::vector<T>, std::vector<T>>& time_values,
    const std::string&                                title = "Step Response"
) {
    return plot_line(time_values.first, time_values.second, title);
}

// ============================================================================
// MATLAB®-style plot helpers (thin plotlypp wrappers over analysis:: results)
// ============================================================================

namespace detail {

/**
 * @brief Build a time-response figure with one line per (output, input) pair
 *
 * Shared implementation for stepplot() and impulseplot(): a `TimeResponse`
 * stores `y[k](i, j)` (output i from a canonical input on channel j), so this
 * emits NY·NU traces. Single-input systems are labelled `y<i>`; multi-input
 * systems `y<i> <- u<j>`.
 *
 * @tparam NY Number of outputs
 * @tparam NU Number of inputs
 * @param resp  Per-channel time response
 * @param title Figure title
 * @return plotlypp::Figure ready for .show() or .writeHtml()
 */
template<size_t NY, size_t NU, typename T>
plotlypp::Figure time_response_figure(const analysis::TimeResponse<NY, NU, T>& resp, const std::string& title) {
    using namespace plotlypp;

    const auto t = to_double_vector(resp.t);
    Figure     fig;
    for (size_t j = 0; j < NU; ++j) {
        for (size_t i = 0; i < NY; ++i) {
            std::vector<double> y;
            y.reserve(resp.y.size());
            for (const auto& yk : resp.y) {
                y.push_back(static_cast<double>(yk(i, j)));
            }
            const std::string name = (NU > 1) ? ("y" + std::to_string(i) + " <- u" + std::to_string(j)) : ("y" + std::to_string(i));
            fig.addTrace(Scatter().x(t).y(y).mode({Scatter::Mode::Lines}).name(name));
        }
    }
    fig.setLayout(
        Layout()
            .title([&](auto& tt) { tt.text(title); })
            .xaxis(Layout::Xaxis().title([](auto& tt) { tt.text("Time (s)"); }))
            .yaxis(Layout::Yaxis().title([](auto& tt) { tt.text("Amplitude"); }))
    );
    return fig;
}

/**
 * @brief Build a markers-only scatter of complex points (real vs imaginary)
 *
 * Used by pzplot() to draw poles and zeros on the complex plane.
 *
 * @param pts  Complex values to plot
 * @param name Trace name (legend label)
 * @param sym  Marker symbol (e.g. X for poles, CircleOpen for zeros)
 * @return plotlypp::Scatter trace
 */
template<typename T>
plotlypp::Scatter complex_scatter(const std::vector<damp::complex<T>>& pts, const std::string& name, plotlypp::Scatter::Marker::Symbol sym) {
    using namespace plotlypp;

    std::vector<double> re, im;
    re.reserve(pts.size());
    im.reserve(pts.size());
    for (const auto& p : pts) {
        re.push_back(static_cast<double>(p.real()));
        im.push_back(static_cast<double>(p.imag()));
    }
    return Scatter()
        .x(re)
        .y(im)
        .mode({Scatter::Mode::Markers})
        .name(name)
        .marker([sym](auto& m) { m.symbol(sym).size(10.0); });
}

} // namespace detail

/**
 * @brief Plot a step response, one trace per input/output pair
 *
 * MATLAB® equivalent: `stepplot(sys)`. Pair with `analysis::step`.
 *
 * @tparam NY Number of outputs
 * @tparam NU Number of inputs
 * @param resp  Step response from analysis::step()
 * @param title Figure title
 * @return plotlypp::Figure
 */
template<size_t NY, size_t NU, typename T>
plotlypp::Figure stepplot(const analysis::TimeResponse<NY, NU, T>& resp, const std::string& title = "Step Response") {
    return detail::time_response_figure(resp, title);
}

/**
 * @brief Plot an impulse response, one trace per input/output pair
 *
 * MATLAB® equivalent: `impulseplot(sys)`. Pair with `analysis::impulse`.
 *
 * @tparam NY Number of outputs
 * @tparam NU Number of inputs
 * @param resp  Impulse response from analysis::impulse()
 * @param title Figure title
 * @return plotlypp::Figure
 */
template<size_t NY, size_t NU, typename T>
plotlypp::Figure impulseplot(const analysis::TimeResponse<NY, NU, T>& resp, const std::string& title = "Impulse Response") {
    return detail::time_response_figure(resp, title);
}

/**
 * @brief Plot a forced (lsim) simulation, one trace per output
 *
 * MATLAB® equivalent: `lsimplot(sys, u, t)`. Pair with `analysis::lsim`.
 *
 * @tparam NX Number of states
 * @tparam NY Number of outputs
 * @param resp  Forced response from analysis::lsim()
 * @param title Figure title
 * @return plotlypp::Figure
 */
template<size_t NX, size_t NY, typename T>
plotlypp::Figure lsimplot(const analysis::LsimResult<NX, NY, T>& resp, const std::string& title = "Linear Simulation") {
    using namespace plotlypp;

    const auto t = to_double_vector(resp.t);
    Figure     fig;
    for (size_t i = 0; i < NY; ++i) {
        fig.addTrace(
            Scatter()
                .x(t)
                .y(extract_channel(resp.y, i))
                .mode({Scatter::Mode::Lines})
                .name("y" + std::to_string(i))
        );
    }
    fig.setLayout(
        Layout()
            .title([&](auto& tt) { tt.text(title); })
            .xaxis(Layout::Xaxis().title([](auto& tt) { tt.text("Time (s)"); }))
            .yaxis(Layout::Yaxis().title([](auto& tt) { tt.text("Amplitude"); }))
    );
    return fig;
}

/**
 * @brief Plot magnitude and phase Bode subplots
 *
 * MATLAB® equivalent: `bodeplot(sys)`. Thin alias of plot_bode().
 *
 * @param bode  Frequency response from analysis::bode()
 * @param title Figure title
 * @return plotlypp::Figure
 */
template<typename T>
plotlypp::Figure bodeplot(const analysis::BodeResult<T>& bode, const std::string& title = "Bode Plot") {
    return plot_bode(bode, title);
}

/**
 * @brief Plot a magnitude-only Bode diagram (log frequency, dB magnitude)
 *
 * MATLAB® equivalent: `bodemag(sys)`.
 *
 * @param bode  Frequency response from analysis::bode()
 * @param title Figure title
 * @return plotlypp::Figure
 */
template<typename T>
plotlypp::Figure bodemag(const analysis::BodeResult<T>& bode, const std::string& title = "Bode Magnitude") {
    using namespace plotlypp;

    std::vector<double> omega, mag_db;
    omega.reserve(bode.points.size());
    mag_db.reserve(bode.points.size());
    for (const auto& pt : bode.points) {
        omega.push_back(static_cast<double>(pt.omega));
        mag_db.push_back(static_cast<double>(pt.magnitude_db));
    }
    return Figure()
        .addTrace(Scatter().x(omega).y(mag_db).mode({Scatter::Mode::Lines}).name("Magnitude"))
        .setLayout(Layout().title([&](auto& tt) { tt.text(title); }).xaxis(Layout::Xaxis().title([](auto& tt) { tt.text("Frequency (rad/s)"); }).type(Layout::Xaxis::Type::Log)).yaxis(Layout::Yaxis().title([](auto& tt) { tt.text("Magnitude (dB)"); })));
}

/**
 * @brief Plot a Nyquist locus with the -1 critical point marked
 *
 * MATLAB® equivalent: `nyquistplot(sys)`. Pair with `analysis::nyquist`.
 *
 * @param nyq   Nyquist response from analysis::nyquist()
 * @param title Figure title
 * @return plotlypp::Figure
 */
template<typename T>
plotlypp::Figure nyquistplot(const analysis::NyquistResult<T>& nyq, const std::string& title = "Nyquist Plot") {
    using namespace plotlypp;

    std::vector<double> re, im;
    re.reserve(nyq.points.size());
    im.reserve(nyq.points.size());
    for (const auto& p : nyq.points) {
        re.push_back(static_cast<double>(p.real));
        im.push_back(static_cast<double>(p.imag));
    }
    Figure fig;
    fig.addTrace(Scatter().x(re).y(im).mode({Scatter::Mode::Lines}).name("L(jw)"));
    fig.addTrace(Scatter().x(std::vector<double>{-1.0}).y(std::vector<double>{0.0}).mode({Scatter::Mode::Markers}).name("-1").marker([](auto& m) { m.symbol(Scatter::Marker::Symbol::X).size(10.0); }));
    fig.setLayout(Layout().title([&](auto& tt) { tt.text(title); }).xaxis(Layout::Xaxis().title([](auto& tt) { tt.text("Real"); })).yaxis(Layout::Yaxis().title([](auto& tt) { tt.text("Imag"); })));
    return fig;
}

/**
 * @brief Plot a pole-zero map on the complex plane (poles as ×, zeros as ○)
 *
 * MATLAB® equivalent: `pzplot(sys)`. Pair with `analysis::pzmap`.
 *
 * @param pz    Pole-zero data from analysis::pzmap()
 * @param title Figure title
 * @return plotlypp::Figure
 */
template<typename T>
plotlypp::Figure pzplot(const analysis::PoleZeroMap<T>& pz, const std::string& title = "Pole-Zero Map") {
    using namespace plotlypp;

    Figure fig;
    fig.addTrace(detail::complex_scatter(pz.poles, "poles", Scatter::Marker::Symbol::X));
    fig.addTrace(detail::complex_scatter(pz.zeros, "zeros", Scatter::Marker::Symbol::CircleOpen));
    fig.setLayout(Layout().title([&](auto& tt) { tt.text(title); }).xaxis(Layout::Xaxis().title([](auto& tt) { tt.text("Real"); })).yaxis(Layout::Yaxis().title([](auto& tt) { tt.text("Imag"); })));
    return fig;
}

namespace detail {

/**
 * @brief Sample an M-circle (constant closed-loop |T|) in the G-plane and map to Nichols axes
 *
 * For @f$ |G/(1+G)| = M @f$ with @f$ M \ne 1 @f$, the locus is a circle of center
 * @f$ -M^2/(1-M^2) @f$ and radius @f$ |M/(1-M^2)| @f$. Samples outside a
 * magnitude window are dropped so the grid stays readable.
 */
inline void add_nichols_m_circle(
    plotlypp::Figure&  fig,
    double             M_lin,
    const std::string& name,
    double             mag_db_min = -40.0,
    double             mag_db_max = 40.0,
    size_t             n_samples = 180
) {
    using namespace plotlypp;
    if (M_lin <= 0.0 || damp::abs(M_lin - 1.0) < 1e-9) {
        return; // M = 1 is the vertical line Re(G) = -1/2 (omitted for simplicity)
    }
    const double M2 = M_lin * M_lin;
    const double den = 1.0 - M2;
    const double cx = -M2 / den;
    const double r = damp::abs(M_lin / den);

    std::vector<double> phase_deg;
    std::vector<double> mag_db;
    phase_deg.reserve(n_samples);
    mag_db.reserve(n_samples);

    constexpr double pi = 3.14159265358979323846;
    for (size_t i = 0; i < n_samples; ++i) {
        const double th = 2.0 * pi * static_cast<double>(i) / static_cast<double>(n_samples);
        const auto [s, c] = damp::sincos(th);
        const double re = cx + r * c;
        const double im = r * s;
        const double mag = damp::sqrt(re * re + im * im);
        if (mag <= 0.0) {
            continue;
        }
        const double mdb = 20.0 * damp::log10(mag);
        if (mdb < mag_db_min || mdb > mag_db_max) {
            continue;
        }
        phase_deg.push_back(damp::atan2(im, re) * 180.0 / pi);
        mag_db.push_back(mdb);
    }
    if (phase_deg.empty()) {
        return;
    }
    fig.addTrace(
        Scatter()
            .x(phase_deg)
            .y(mag_db)
            .mode({Scatter::Mode::Lines})
            .name(name)
            .line([](auto& ln) { ln.width(1.0).dash("dot"); })
    );
}

} // namespace detail

/**
 * @brief Plot a Nichols chart (open-loop phase vs magnitude) with M-circle grid
 *
 * Locus from @ref analysis::nichols / Bode data; overlays constant closed-loop
 * magnitude (M) contours and marks the critical point @f$(-180^\circ,\,0\,\mathrm{dB})@f$.
 *
 * MATLAB® equivalent: `nicholsplot(sys)`. Pair with `analysis::nichols`.
 *
 * @param nich  Nichols response from analysis::nichols()
 * @param title Figure title
 * @return plotlypp::Figure
 *
 * @note Compare with MATLAB®'s nicholsplot.
 */
template<typename T>
plotlypp::Figure nicholsplot(const analysis::NicholsResult<T>& nich, const std::string& title = "Nichols Chart") {
    using namespace plotlypp;

    std::vector<double> phase_deg;
    std::vector<double> mag_db;
    phase_deg.reserve(nich.points.size());
    mag_db.reserve(nich.points.size());
    for (const auto& p : nich.points) {
        phase_deg.push_back(static_cast<double>(p.phase_deg));
        mag_db.push_back(static_cast<double>(p.magnitude_db));
    }

    Figure fig;
    // M-grid first so the locus draws on top
    const double m_db[] = {12.0, 6.0, 3.0, 1.0, -1.0, -3.0, -6.0, -12.0};
    for (double mdb : m_db) {
        const double M = damp::pow(10.0, mdb / 20.0);
        detail::add_nichols_m_circle(fig, M, "M = " + std::to_string(static_cast<int>(mdb)) + " dB");
    }

    fig.addTrace(Scatter().x(phase_deg).y(mag_db).mode({Scatter::Mode::Lines}).name("L(jw)"));
    fig.addTrace(
        Scatter()
            .x(std::vector<double>{-180.0})
            .y(std::vector<double>{0.0})
            .mode({Scatter::Mode::Markers})
            .name("critical (-180, 0 dB)")
            .marker([](auto& m) { m.symbol(Scatter::Marker::Symbol::X).size(10.0); })
    );
    fig.setLayout(
        Layout()
            .title([&](auto& tt) { tt.text(title); })
            .xaxis(Layout::Xaxis().title([](auto& tt) { tt.text("Open-loop phase (deg)"); }))
            .yaxis(Layout::Yaxis().title([](auto& tt) { tt.text("Open-loop magnitude (dB)"); }))
    );
    return fig;
}

/**
 * @brief Plot singular-value frequency response (log frequency, dB)
 *
 * One trace per singular value @f$\sigma_i@f$ (largest first). Pair with
 * `analysis::sigma`.
 *
 * MATLAB® equivalent: `sigmaplot(sys)` / `sigma(sys)`.
 *
 * @tparam NS Number of singular values (= min(NY, NU))
 * @param sig   SigmaResult from analysis::sigma()
 * @param title Figure title
 * @return plotlypp::Figure
 *
 * @note Compare with MATLAB®'s sigmaplot / sigma.
 */
template<size_t NS, typename T>
plotlypp::Figure sigmaplot(const analysis::SigmaResult<NS, T>& sig, const std::string& title = "Singular Values") {
    using namespace plotlypp;

    std::vector<double> omega;
    omega.reserve(sig.points.size());
    for (const auto& pt : sig.points) {
        omega.push_back(static_cast<double>(pt.omega));
    }

    Figure fig;
    for (size_t k = 0; k < NS; ++k) {
        std::vector<double> sdb;
        sdb.reserve(sig.points.size());
        for (const auto& pt : sig.points) {
            sdb.push_back(static_cast<double>(pt.sigma_db[k]));
        }
        const std::string name = (NS == 1) ? "sigma" : ("sigma" + std::to_string(k + 1));
        fig.addTrace(Scatter().x(omega).y(sdb).mode({Scatter::Mode::Lines}).name(name));
    }
    fig.setLayout(
        Layout()
            .title([&](auto& tt) { tt.text(title); })
            .xaxis(Layout::Xaxis().title([](auto& tt) { tt.text("Frequency (rad/s)"); }).type(Layout::Xaxis::Type::Log))
            .yaxis(Layout::Yaxis().title([](auto& tt) { tt.text("Singular value (dB)"); }))
    );
    return fig;
}

/**
 * @brief Plot a root locus (closed-loop poles vs gain) on the complex plane
 *
 * One trace per pole branch; open-loop poles (×) and zeros (○) when present.
 * Pair with `analysis::rlocus`.
 *
 * MATLAB® equivalent: `rlocusplot(sys)`.
 *
 * @param rl    Root-locus data from analysis::rlocus()
 * @param title Figure title
 * @return plotlypp::Figure
 *
 * @note Compare with MATLAB®'s rlocusplot.
 */
template<typename T>
plotlypp::Figure rlocusplot(const analysis::RootLocusResult<T>& rl, const std::string& title = "Root Locus") {
    using namespace plotlypp;

    Figure fig;

    // Branch traces: poles[gain_idx][branch] — transpose to one path per branch
    size_t n_branches = 0;
    for (const auto& pk : rl.poles) {
        n_branches = damp::max(n_branches, pk.size());
    }
    for (size_t b = 0; b < n_branches; ++b) {
        std::vector<double> re;
        std::vector<double> im;
        re.reserve(rl.poles.size());
        im.reserve(rl.poles.size());
        for (const auto& pk : rl.poles) {
            if (b < pk.size()) {
                re.push_back(static_cast<double>(pk[b].real()));
                im.push_back(static_cast<double>(pk[b].imag()));
            }
        }
        fig.addTrace(
            Scatter()
                .x(re)
                .y(im)
                .mode({Scatter::Mode::Lines})
                .name("branch " + std::to_string(b + 1))
        );
    }

    if (!rl.open_loop_poles.empty()) {
        fig.addTrace(detail::complex_scatter(rl.open_loop_poles, "open-loop poles", Scatter::Marker::Symbol::X));
    }
    if (!rl.open_loop_zeros.empty()) {
        fig.addTrace(detail::complex_scatter(rl.open_loop_zeros, "open-loop zeros", Scatter::Marker::Symbol::CircleOpen));
    }

    fig.setLayout(
        Layout()
            .title([&](auto& tt) { tt.text(title); })
            .xaxis(Layout::Xaxis().title([](auto& tt) { tt.text("Real"); }))
            .yaxis(Layout::Yaxis().title([](auto& tt) { tt.text("Imag"); }))
    );
    return fig;
}

} // namespace plot
} // namespace damp
