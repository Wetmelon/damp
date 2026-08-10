// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file plot_check.hpp
 * @brief Test-only Plotly dumps for human inspection (not pass/fail)
 *
 * ## Plotting policy (tests vs examples vs library)
 *
 * | Audience | Role | Where | How |
 * |----------|------|-------|-----|
 * | **Unit tests** | Eyeball what `CHECK`s already proved | `tests/plots/` HTML (gitignored) | This header (`plotcheck::xy`) |
 * | **Examples** | Teach Design → SIL / story | `examples/plots/` HTML (gitignored) | Prefer `damp::plot::*`; always `plot::write_html` |
 * | **Library** | damp result types → figure | host-only | `damp/simulation/plot_plotly.hpp` (`damp::plot`) |
 *
 * Rules:
 * - **CHECK / static_assert own pass/fail.** Plots are side-output only. Never
 *   make an HTML file the oracle or a required CI artifact.
 * - **Tests always write under `tests/plots/`** — not the process cwd — so
 *   multi-config / CI paths stay predictable.
 * - **Filename = scenario** (`hybrid_buck_dcm_step.html`), not the suite name.
 * - Prefer this helper over raw plotlypp in tests. For multi-row sim figures,
 *   build with `damp::plot::plot_simulation` then write into `tests/plots/`.
 * - **Default multi-panel format** (`damp::plot` in `plot_plotly.hpp`):
 *   secondary x-axes use `.matches("x")` + `.domain(panel_x_domain())`; each
 *   panel has its own legend via `panel_legend_id(row)` +
 *   `Layout::legend(i, panel_legend(y))` (legend at `panel_x_right + 0.01`,
 *   default `panel_x_right = 1.0`); colors restart per panel with
 *   `panel_color(i)`. Do not hand-roll `Layout::Legend{x: 1.02, …}` or match
 *   axes that intentionally use different domains (e.g. a zoom inset).
 *   Write via `plot::write_html` (damp shell; do not call `Figure::writeHtml`).
 * - Plotting is **host-only** (plotlypp + `<vector>`). Never pull it into
 *   `control.hpp` or target firmware.
 *
 * Phase 1 hybrid sim: unit tests assert event times / mode sequence / bounds;
 * optional `plotcheck` dumps switch state and waveforms; the buck *example*
 * is the full narrative + nice multi-panel figure.
 */

#include <filesystem>
#include <string>
#include <vector>

#include "damp/simulation/plot_plotly.hpp"

namespace plotcheck {

struct Series {
    std::string         name;
    std::vector<double> x;
    std::vector<double> y;
    bool                markers = false; //!< true → markers-only; false → lines-only
};

/// Write a multi-series planar (x, y) scatter to tests/plots/<file>.
/// Prefer this for paths / loci / phase portraits — not two independent 1-D lines.
inline void xy(
    const std::string&         file,
    const std::string&         title,
    const std::string&         xlabel,
    const std::string&         ylabel,
    const std::vector<Series>& series,
    bool                       equal_aspect = false
) {
    std::vector<damp::plot::XySeries> xs;
    xs.reserve(series.size());
    for (const auto& s : series) {
        xs.push_back({
            .name = s.name,
            .x = s.x,
            .y = s.y,
            .lines = !s.markers,
            .markers = s.markers,
            .color = {},
        });
    }
    auto fig = damp::plot::plot_xy(xs, title, xlabel, ylabel, equal_aspect);
    std::filesystem::create_directories("tests/plots");
    damp::plot::write_html(fig, "tests/plots/" + file);
}

} // namespace plotcheck
