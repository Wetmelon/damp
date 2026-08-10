// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file examples/animate_2d.hpp
 * @brief Host-only multi-panel 2D animation (canvas) with play / scrub / theme.
 *
 * Companion to `animate_3d.hpp` for planar teaching plots. Supports mixed path
 * (2D parametric) and time-series panels in an arbitrary row×col grid.
 *
 * Static full curves stay on screen; the playhead draws the trail, time cursor,
 * and current samples. Not part of the embeddable library.
 */

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <string>
#include <utility>
#include <vector>

#include "fmt/core.h"
#include "fmt/format.h"
#include "nlohmann/json.hpp"

namespace damp::examples_plot2d {

using json = nlohmann::json;

[[nodiscard]] inline std::string json_js_literal(const json& j) {
    std::string s = j.dump();
    size_t      pos = 0;
    while ((pos = s.find('<', pos)) != std::string::npos) {
        s.replace(pos, 1, "\\u003c");
        pos += 6;
    }
    return s;
}

// ---------------------------------------------------------------------------
// General multi-panel grid
// ---------------------------------------------------------------------------

/**
 * @brief One time-series channel in a Time panel
 */
struct Series2d {
    std::string         name;
    std::string         color;
    std::vector<double> y; //!< length N (same as times)
};

/// Extra static path curve overlaid on a Path panel (e.g. cubic vs quintic map).
struct PathOverlay2d {
    std::string         name;
    std::string         color;
    std::vector<double> x;
    std::vector<double> y;
};

/**
 * @brief One cell in the animated grid (path or time series)
 */
struct Panel2d {
    enum class Kind { Path,
                      Time };

    Kind        kind{Kind::Time};
    std::string xlabel;
    std::string ylabel;

    // Path (Kind::Path): static curve (path_x, path_y) — any length; optional
    // per-frame playhead (path_mark_x/y, length N). If marks empty, playhead uses
    // path_x[i]/path_y[i] (curve must be length N). path_overlays = more curves.
    std::vector<double>        path_x;
    std::vector<double>        path_y;
    std::vector<double>        path_mark_x; ///< playhead x per frame (optional)
    std::vector<double>        path_mark_y; ///< playhead y per frame (optional)
    std::string                path_color{"#7c3aed"};
    std::vector<PathOverlay2d> path_overlays;
    /// If true (default), force a square data domain so XY paths keep aspect.
    /// Set false for non-spatial maps (e.g. cam s vs m) so the plot fills the
    /// panel width and lines up with time panels in the same column.
    bool path_equal_aspect{true};

    // Time (Kind::Time): one or more y(t) series against shared `times`
    std::vector<Series2d> series;
    bool                  y_auto{true}; //!< auto y-extent from data (+pad)
    double                y_lo{-1};
    double                y_hi{1};
    std::vector<double>   vlines; //!< vertical markers in data-x (time)
};

/**
 * @brief Write a multi-panel 2D animation HTML page
 *
 * Features: Play/Pause, Space toggle, scrubber, theme (system/light/dark),
 * grid + tick labels + axis titles, **draggable purple cursors** on every panel
 * (vertical time cursor on Time plots; crosshair on Path plots — both scrub the
 * shared playhead), path ghost/trail/marker (optional equal aspect), time faint curves /
 * trail / markers + vlines, column headers under the control bar, optional row
 * headers.
 *
 * @param path         e.g. plots/foo.html
 * @param title        page title
 * @param times        sample times [s], length N
 * @param dt_s         wall-clock step between frames for real-time play
 * @param nrows        grid rows
 * @param ncols        grid columns
 * @param panels       row-major panels, size nrows*ncols
 * @param col_headers  optional column titles (length ncols)
 * @param row_headers  optional row titles drawn on the left (length nrows)
 */
inline void write_animated_grid(
    const std::filesystem::path&    path,
    const std::string&              title,
    const std::vector<double>&      times,
    double                          dt_s,
    int                             nrows,
    int                             ncols,
    const std::vector<Panel2d>&     panels,
    const std::vector<std::string>& col_headers = {},
    const std::vector<std::string>& row_headers = {}
) {
    namespace fs = std::filesystem;
    if (times.empty() || nrows < 1 || ncols < 1) {
        return;
    }
    const size_t expect = static_cast<size_t>(nrows) * static_cast<size_t>(ncols);
    if (panels.size() < expect) {
        return;
    }
    if (!path.parent_path().empty()) {
        fs::create_directories(path.parent_path());
    }

    json jpanels = json::array();
    for (size_t i = 0; i < expect; ++i) {
        const auto& p = panels[i];
        json        jp;
        jp["kind"] = (p.kind == Panel2d::Kind::Path) ? "path" : "time";
        jp["xlabel"] = p.xlabel;
        jp["ylabel"] = p.ylabel;
        if (p.kind == Panel2d::Kind::Path) {
            jp["path_x"] = p.path_x;
            jp["path_y"] = p.path_y;
            jp["path_mark_x"] = p.path_mark_x;
            jp["path_mark_y"] = p.path_mark_y;
            jp["path_color"] = p.path_color.empty() ? "#7c3aed" : p.path_color;
            jp["equal_aspect"] = p.path_equal_aspect;
            json jov = json::array();
            for (const auto& ov : p.path_overlays) {
                jov.push_back({
                    {"name", ov.name},
                    {"color", ov.color.empty() ? "#1f77b4" : ov.color},
                    {"x", ov.x},
                    {"y", ov.y},
                });
            }
            jp["path_overlays"] = jov;
        } else {
            json jser = json::array();
            for (const auto& s : p.series) {
                jser.push_back({
                    {"name", s.name},
                    {"color", s.color.empty() ? "#1f77b4" : s.color},
                    {"y", s.y},
                });
            }
            jp["series"] = jser;
            jp["y_auto"] = p.y_auto;
            jp["y_lo"] = p.y_lo;
            jp["y_hi"] = p.y_hi;
            jp["vlines"] = p.vlines;
        }
        jpanels.push_back(jp);
    }

    json jcol = json::array();
    for (int c = 0; c < ncols; ++c) {
        if (c < static_cast<int>(col_headers.size())) {
            jcol.push_back(col_headers[static_cast<size_t>(c)]);
        } else {
            jcol.push_back("");
        }
    }
    json jrow = json::array();
    for (int r = 0; r < nrows; ++r) {
        if (r < static_cast<int>(row_headers.size())) {
            jrow.push_back(row_headers[static_cast<size_t>(r)]);
        } else {
            jrow.push_back("");
        }
    }

    const std::string html = std::string(R"(<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>)") + title + R"(</title>
  <style>
    :root, [data-theme="light"] {
      --bg: #f4f4f5; --fg: #18181b; --muted: #52525b;
      --btn-bg: #fff; --btn-border: #d4d4d8; --btn-hover: #f4f4f5;
      --panel: #fafafa; --grid: #d4d4d8; --grid-major: #a1a1aa;
      --axis: #71717a; --ghost: #a1a1aa;
      --cursor: #a855f7; --trail: #7c3aed;
    }
    [data-theme="dark"] {
      --bg: #111113; --fg: #e4e4e7; --muted: #a1a1aa;
      --btn-bg: #27272a; --btn-border: #52525b; --btn-hover: #3f3f46;
      --panel: #1c1c1f; --grid: #2e2e33; --grid-major: #52525b;
      --axis: #71717a; --ghost: #52525b;
      --cursor: #c084fc; --trail: #a78bfa;
    }
    html, body {
      margin: 0; height: 100%; overflow: hidden;
      background: var(--bg); color: var(--fg);
      font-family: system-ui, sans-serif;
    }
    #bar {
      position: absolute; z-index: 10; top: 0; left: 0; right: 0;
      display: flex; flex-wrap: wrap; gap: 10px; align-items: center;
      padding: 8px 12px;
      min-height: 40px;
      box-sizing: border-box;
      background: color-mix(in srgb, var(--bg) 92%, transparent);
      backdrop-filter: blur(4px);
      pointer-events: none;
    }
    #bar > * { pointer-events: auto; }
    #bar button {
      cursor: pointer; padding: 6px 14px; border-radius: 4px;
      border: 1px solid var(--btn-border); background: var(--btn-bg); color: var(--fg);
    }
    #bar button:hover { background: var(--btn-hover); }
    #tlabel { font-size: 13px; min-width: 7em; }
    #hint { font-size: 12px; color: var(--muted); flex: 1 1 auto; }
    #scrubWrap {
      position: absolute; z-index: 10; left: 12px; right: 12px; bottom: 8px;
      height: 28px; display: flex; align-items: center; gap: 10px;
      pointer-events: none;
    }
    #scrubWrap > * { pointer-events: auto; }
    #scrub { flex: 1; height: 28px; cursor: pointer; }
    #c { display: block; width: 100%; height: 100%; touch-action: none; }
  </style>
</head>
<body>
  <div id="bar">
    <button type="button" id="btnPlay">Play</button>
    <button type="button" id="btnPause">Pause</button>
    <button type="button" id="btnTheme" title="Theme: system / light / dark">Theme: system</button>
    <span id="tlabel">t = 0.00 s</span>
    <span id="hint">Space = play/pause · drag purple cursor on any panel · scrubber</span>
  </div>
  <canvas id="c"></canvas>
  <div id="scrubWrap">
    <input type="range" id="scrub" min="0" max="0" value="0" step="1" />
  </div>
  <script>
    const times = )" + json_js_literal(json(times))
                           + R"(;
    const dtSec = )" + std::to_string(dt_s)
                           + R"(;
    const nRows = )" + std::to_string(nrows)
                           + R"(;
    const nCols = )" + std::to_string(ncols)
                           + R"(;
    const panels = )" + json_js_literal(jpanels)
                           + R"(;
    const colHeaders = )" + json_js_literal(jcol)
                           + R"(;
    const rowHeaders = )" + json_js_literal(jrow)
                           + R"(;
    const nFrames = times.length;

    const canvas = document.getElementById('c');
    const ctx = canvas.getContext('2d');
    const scrub = document.getElementById('scrub');
    const btnTheme = document.getElementById('btnTheme');
    scrub.max = Math.max(0, nFrames - 1);

    const THEME_KEY = 'damp-anim-2d-theme';
    let themePref = localStorage.getItem(THEME_KEY) || 'system';

    function resolvedTheme() {
      if (themePref === 'light' || themePref === 'dark') return themePref;
      return window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
    }
    function applyTheme() {
      document.documentElement.setAttribute('data-theme', resolvedTheme());
      btnTheme.textContent = 'Theme: ' + themePref;
      readCss();
      draw();
    }
    btnTheme.onclick = function () {
      switch (themePref) {
        case 'system': themePref = 'light'; break;
        case 'light':  themePref = 'dark'; break;
        default:       themePref = 'system'; break;
      }
      localStorage.setItem(THEME_KEY, themePref);
      applyTheme();
    };
    window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', function () {
      if (themePref === 'system') applyTheme();
    });

    let frameIndex = 0;
    let playing = false;
    let playOriginWall = 0;
    let playOriginFrame = 0;
    let css = null;

    function readCss() {
      const s = getComputedStyle(document.documentElement);
      css = {
        bg: s.getPropertyValue('--bg').trim(),
        fg: s.getPropertyValue('--fg').trim(),
        muted: s.getPropertyValue('--muted').trim(),
        panel: s.getPropertyValue('--panel').trim(),
        grid: s.getPropertyValue('--grid').trim(),
        gridMajor: s.getPropertyValue('--grid-major').trim(),
        axis: s.getPropertyValue('--axis').trim(),
        ghost: s.getPropertyValue('--ghost').trim(),
        cursor: s.getPropertyValue('--cursor').trim(),
        trail: s.getPropertyValue('--trail').trim(),
      };
    }

    function niceStep(span, target) {
      if (!(span > 0) || !isFinite(span)) return 1;
      const raw = span / Math.max(target, 1);
      const pow = Math.pow(10, Math.floor(Math.log10(raw)));
      const n = raw / pow;
      let m = 1;
      if (n > 5) m = 10;
      else if (n > 2) m = 5;
      else if (n > 1) m = 2;
      return m * pow;
    }
    function ticks(lo, hi, target) {
      const step = niceStep(hi - lo, target);
      const start = Math.ceil((lo - 1e-12 * step) / step) * step;
      const out = [];
      for (let v = start; v <= hi + 1e-9 * step; v += step) {
        const t = Math.abs(v) < step * 1e-10 ? 0 : v;
        out.push(t);
        if (out.length > 40) break;
      }
      return { values: out, step: step };
    }
    function fmtTick(v, step) {
      const a = Math.abs(step);
      let dig = 0;
      if (a > 0 && a < 1) dig = Math.min(4, Math.ceil(-Math.log10(a)));
      else if (a >= 1 && a < 10) dig = 1;
      return Number(v).toFixed(dig);
    }

    function mapX(v, lo, hi, x0, x1) {
      return x0 + (x1 - x0) * (v - lo) / (hi - lo);
    }
    function mapY(v, lo, hi, y0, y1) {
      return y1 - (y1 - y0) * (v - lo) / (hi - lo);
    }

    function drawGrid(px0, py0, px1, py1, xlo, xhi, ylo, yhi, opts) {
      const o = opts || {};
      const nx = o.nx || 5;
      const ny = o.ny || 5;
      const showX = o.showX !== false;
      const showY = o.showY !== false;
      const xt = ticks(xlo, xhi, nx);
      const yt = ticks(ylo, yhi, ny);

      ctx.save();
      ctx.beginPath();
      ctx.rect(px0, py0, px1 - px0, py1 - py0);
      ctx.clip();
      for (const v of xt.values) {
        const X = mapX(v, xlo, xhi, px0, px1);
        const major = Math.abs(v) < 1e-12 * Math.max(xt.step, 1);
        ctx.strokeStyle = major ? css.gridMajor : css.grid;
        ctx.lineWidth = major ? 1.25 : 1;
        ctx.beginPath();
        ctx.moveTo(X, py0);
        ctx.lineTo(X, py1);
        ctx.stroke();
      }
      for (const v of yt.values) {
        const Y = mapY(v, ylo, yhi, py0, py1);
        const major = Math.abs(v) < 1e-12 * Math.max(yt.step, 1);
        ctx.strokeStyle = major ? css.gridMajor : css.grid;
        ctx.lineWidth = major ? 1.25 : 1;
        ctx.beginPath();
        ctx.moveTo(px0, Y);
        ctx.lineTo(px1, Y);
        ctx.stroke();
      }
      ctx.restore();

      ctx.strokeStyle = css.axis;
      ctx.lineWidth = 1.25;
      ctx.strokeRect(px0 + 0.5, py0 + 0.5, px1 - px0 - 1, py1 - py0 - 1);

      ctx.fillStyle = css.muted;
      ctx.font = '9px system-ui, sans-serif';
      if (showX) {
        ctx.textAlign = 'center';
        ctx.textBaseline = 'top';
        for (const v of xt.values) {
          const X = mapX(v, xlo, xhi, px0, px1);
          if (X < px0 - 2 || X > px1 + 2) continue;
          ctx.fillText(fmtTick(v, xt.step), X, py1 + 3);
        }
      }
      if (showY) {
        ctx.textAlign = 'right';
        ctx.textBaseline = 'middle';
        for (const v of yt.values) {
          const Y = mapY(v, ylo, yhi, py0, py1);
          if (Y < py0 - 2 || Y > py1 + 2) continue;
          ctx.fillText(fmtTick(v, yt.step), px0 - 4, Y);
        }
      }
    }

    function tText(i) {
      const t = (times && times.length > i) ? times[i] : (i * dtSec);
      return 't = ' + Number(t).toFixed(2) + ' s';
    }
    function setChrome(i) {
      document.getElementById('tlabel').textContent = tText(i);
      if (document.activeElement !== scrub) scrub.value = String(i);
    }

    // Raw min/max (no pad) — used for time axis so left/right plot edges are
    // exactly t₀ and t_end (dragging the purple cursor tracks the pointer).
    function extentRaw(arr) {
      let lo = Infinity, hi = -Infinity;
      for (let i = 0; i < arr.length; i++) {
        const v = arr[i];
        if (v < lo) lo = v;
        if (v > hi) hi = v;
      }
      if (!isFinite(lo)) { lo = 0; hi = 1; }
      if (hi - lo < 1e-12) { lo -= 1; hi += 1; }
      return [lo, hi];
    }
    function extent(arr) {
      const [lo, hi] = extentRaw(arr);
      const pad = 0.08 * (hi - lo);
      return [lo - pad, hi + pad];
    }
    function extent2(xs, ys) {
      const [x0, x1] = extent(xs);
      const [y0, y1] = extent(ys);
      const cx = 0.5 * (x0 + x1), cy = 0.5 * (y0 + y1);
      let hx = 0.5 * (x1 - x0), hy = 0.5 * (y1 - y0);
      const h = Math.max(hx, hy);
      return [cx - h, cx + h, cy - h, cy + h];
    }
    function extentSeriesY(series) {
      let lo = Infinity, hi = -Infinity;
      for (const s of series) {
        for (let i = 0; i < s.y.length; i++) {
          const v = s.y[i];
          if (v < lo) lo = v;
          if (v > hi) hi = v;
        }
      }
      if (!isFinite(lo)) { lo = 0; hi = 1; }
      if (hi - lo < 1e-12) { lo -= 1; hi += 1; }
      const pad = 0.08 * (hi - lo);
      return [lo - pad, hi + pad];
    }

    // Precompute axis ranges per panel.
    const panelBox = panels.map((p) => {
      if (p.kind === 'path') {
        let xs = (p.path_x || []).concat(p.path_mark_x || []);
        let ys = (p.path_y || []).concat(p.path_mark_y || []);
        for (const ov of (p.path_overlays || [])) {
          xs = xs.concat(ov.x || []);
          ys = ys.concat(ov.y || []);
        }
        // equal_aspect (default true): square data domain for planar XY.
        // false (cam s vs m): tight x like time plots (left = m_min, right = m_max);
        // y still padded so curves are not flush with top/bottom.
        const box = (p.equal_aspect === false)
          ? (() => { const [x0, x1] = extentRaw(xs); const [y0, y1] = extent(ys); return [x0, x1, y0, y1]; })()
          : extent2(xs, ys);
        return { type: 'path', box };
      }
      let ylo, yhi;
      if (p.y_auto) {
        [ylo, yhi] = extentSeriesY(p.series || []);
      } else {
        ylo = p.y_lo;
        yhi = p.y_hi;
        if (!(yhi > ylo)) { ylo = -1; yhi = 1; }
      }
      return { type: 'time', y: [ylo, yhi] };
    });
    // Tight time range: plot left edge = first sample, right edge = last sample.
    const tRange = (() => {
      if (!times.length) return [0, 1];
      const lo = times[0];
      const hi = times[times.length - 1];
      return (hi > lo) ? [lo, hi] : [lo, lo + 1];
    })();

    // Break the stroke when a segment jumps across most of the axis (e.g. leader
    // wrap 0.99 → 0.01 would otherwise draw a false chord across the plot).
    function strokePolyline(xs, ys, i0, i1, xlo, xhi, ylo, yhi, px0, px1, py0, py1) {
      if (i1 <= i0) return;
      const n = Math.min(xs.length, ys.length);
      if (n === 0) return;
      const a = Math.max(0, Math.min(i0, n - 1));
      const b = Math.max(0, Math.min(i1, n - 1));
      if (b < a) return;
      const xSpan = Math.max(Math.abs(xhi - xlo), 1e-12);
      const ySpan = Math.max(Math.abs(yhi - ylo), 1e-12);
      const jumpX = 0.45 * xSpan;
      const jumpY = 0.45 * ySpan;
      ctx.beginPath();
      let pen = false;
      let prevX = 0, prevY = 0;
      for (let i = a; i <= b; i++) {
        const X = mapX(xs[i], xlo, xhi, px0, px1);
        const Y = mapY(ys[i], ylo, yhi, py0, py1);
        if (!pen) {
          ctx.moveTo(X, Y);
          pen = true;
        } else if (Math.abs(xs[i] - xs[i - 1]) > jumpX || Math.abs(ys[i] - ys[i - 1]) > jumpY) {
          ctx.moveTo(X, Y);
        } else {
          ctx.lineTo(X, Y);
        }
        prevX = X; prevY = Y;
      }
      ctx.stroke();
    }

    function drawPanelBg(px0, py0, pw, ph) {
      ctx.fillStyle = css.panel;
      ctx.strokeStyle = css.grid;
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.rect(px0, py0, pw, ph);
      ctx.fill();
      ctx.stroke();
    }

    function drawText(str, x, y, opts) {
      if (!str) return;
      ctx.fillStyle = (opts && opts.color) || css.fg;
      ctx.font = (opts && opts.font) || '12px system-ui, sans-serif';
      ctx.textAlign = (opts && opts.align) || 'left';
      ctx.textBaseline = (opts && opts.baseline) || 'middle';
      ctx.fillText(str, x, y);
    }

    function hasColHeaders() {
      return colHeaders.some((h) => h && String(h).length > 0);
    }
    function hasRowHeaders() {
      return rowHeaders.some((h) => h && String(h).length > 0);
    }

    // Plot box padding inside each cell (must match drawFrame).
    const PLOT_PAD = { L: 32, R: 6, T: 6, B: 22 };

    let cachedCssW = 0, cachedCssH = 0;

    function ensureCanvasSize() {
      const dpr = Math.min(window.devicePixelRatio || 1, 2);
      const W = canvas.clientWidth || window.innerWidth || 1;
      const H = canvas.clientHeight || window.innerHeight || 1;
      if (W !== cachedCssW || H !== cachedCssH || canvas.width !== Math.floor(W * dpr)) {
        cachedCssW = W;
        cachedCssH = H;
        canvas.width = Math.floor(W * dpr);
        canvas.height = Math.floor(H * dpr);
      }
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      return { W, H, dpr };
    }

    function layout() {
      const { W, H } = ensureCanvasSize();

      const chromeH = 44;
      const colHeaderH = hasColHeaders() ? 22 : 0;
      const rowLabelW = hasRowHeaders() ? 28 : 0;
      const left = 16 + rowLabelW, right = 12, top = chromeH + colHeaderH, bottom = 48;
      const hgap = 12, vgap = 10;
      const innerW = W - left - right;
      const innerH = H - top - bottom;
      const colw = (innerW - Math.max(0, nCols - 1) * hgap) / nCols;
      const rowh = (innerH - Math.max(0, nRows - 1) * vgap) / nRows;

      const cells = [];
      for (let r = 0; r < nRows; r++) {
        for (let c = 0; c < nCols; c++) {
          const x = left + c * (colw + hgap);
          const y = top + r * (rowh + vgap);
          cells.push({ r, c, x, y, w: colw, h: rowh, idx: r * nCols + c });
        }
      }
      return { W, H, left, top, chromeH, colHeaderH, rowLabelW, cells, colw, rowh };
    }

    function plotBox(cell) {
      return {
        px0: cell.x + PLOT_PAD.L,
        py0: cell.y + PLOT_PAD.T,
        px1: cell.x + cell.w - PLOT_PAD.R,
        py1: cell.y + cell.h - PLOT_PAD.B,
      };
    }

    // Hit targets for scrubbing (rebuilt each draw; scrub recomputes geometry live).
    const hitPlots = [];

    function drawAxisTitles(xlabel, ylabel, px0, py0, px1, py1, panel) {
      drawText(xlabel, (px0 + px1) / 2, panel.y + panel.h - 2, {
        align: 'center', color: css.muted, font: '10px system-ui', baseline: 'bottom',
      });
      ctx.save();
      ctx.translate(panel.x + 10, (py0 + py1) / 2);
      ctx.rotate(-Math.PI / 2);
      drawText(ylabel, 0, 0, { align: 'center', color: css.muted, font: '10px system-ui' });
      ctx.restore();
    }

    function drawFrame(i) {
      if (!css) readCss();
      const L = layout();
      ctx.fillStyle = css.bg;
      ctx.fillRect(0, 0, L.W, L.H);

      // Column headers under control bar.
      if (L.colHeaderH > 0) {
        const headerY = L.chromeH + L.colHeaderH * 0.5;
        for (let c = 0; c < nCols; c++) {
          const cell = L.cells[c]; // r=0
          drawText(colHeaders[c], cell.x + cell.w / 2, headerY, {
            align: 'center', font: '13px system-ui,sans-serif', baseline: 'middle',
          });
        }
      }
      // Row headers on the left.
      if (L.rowLabelW > 0) {
        for (let r = 0; r < nRows; r++) {
          const cell = L.cells[r * nCols];
          ctx.save();
          ctx.translate(12, cell.y + cell.h / 2);
          ctx.rotate(-Math.PI / 2);
          drawText(rowHeaders[r], 0, 0, {
            align: 'center', color: css.muted, font: '12px system-ui,sans-serif',
          });
          ctx.restore();
        }
      }

      hitPlots.length = 0;

      for (const cell of L.cells) {
        const p = panels[cell.idx];
        const box = panelBox[cell.idx];
        const { px0, py0, px1, py1 } = plotBox(cell);
        drawPanelBg(cell.x, cell.y, cell.w, cell.h);

        if (p.kind === 'path') {
          const [xlo, xhi, ylo, yhi] = box.box;
          drawGrid(px0, py0, px1, py1, xlo, xhi, ylo, yhi, { nx: 5, ny: 5 });

          ctx.save();
          ctx.beginPath();
          ctx.rect(px0, py0, px1 - px0, py1 - py0);
          ctx.clip();

          const xs = p.path_x || [];
          const ys = p.path_y || [];
          const mx = p.path_mark_x || [];
          const my = p.path_mark_y || [];
          const hasMarks = mx.length > 0 && my.length > 0;
          const nCurve = Math.min(xs.length, ys.length);
          const nMark = hasMarks
            ? Math.min(mx.length, my.length, nFrames)
            : Math.min(xs.length, ys.length, nFrames);
          const lastCurve = Math.max(0, nCurve - 1);
          const ii = Math.min(i, Math.max(0, nMark - 1));

          // Static map / path (full curve). Jump-break avoids wrap chords.
          if (nCurve > 0) {
            const curveCol = p.path_color || css.trail;
            if (hasMarks) {
              // Full colored map (cam table); playhead is separate.
              ctx.strokeStyle = curveCol;
              ctx.lineWidth = 2.25;
              strokePolyline(xs, ys, 0, lastCurve, xlo, xhi, ylo, yhi, px0, px1, py0, py1);
            } else {
              // Legacy path-as-time-series: ghost full curve + colored trail to ii.
              ctx.strokeStyle = css.ghost;
              ctx.lineWidth = 1.75;
              strokePolyline(xs, ys, 0, lastCurve, xlo, xhi, ylo, yhi, px0, px1, py0, py1);
              if (nMark > 0) {
                ctx.strokeStyle = curveCol;
                ctx.lineWidth = 2.5;
                strokePolyline(xs, ys, 0, ii, xlo, xhi, ylo, yhi, px0, px1, py0, py1);
              }
            }
          }
          // Additional static overlays (e.g. cubic vs quintic on the same axes).
          for (const ov of (p.path_overlays || [])) {
            const ox = ov.x || [];
            const oy = ov.y || [];
            const on = Math.min(ox.length, oy.length);
            if (on < 2) continue;
            ctx.strokeStyle = ov.color || css.trail;
            ctx.lineWidth = 2.0;
            ctx.setLineDash([6, 4]);
            strokePolyline(ox, oy, 0, on - 1, xlo, xhi, ylo, yhi, px0, px1, py0, py1);
            ctx.setLineDash([]);
          }

          if (nMark > 0) {
            const hx = hasMarks ? mx[ii] : xs[ii];
            const hy = hasMarks ? my[ii] : ys[ii];
            const X = mapX(hx, xlo, xhi, px0, px1);
            const Y = mapY(hy, ylo, yhi, py0, py1);
            ctx.strokeStyle = css.cursor;
            ctx.lineWidth = 1.25;
            ctx.setLineDash([4, 3]);
            ctx.beginPath();
            ctx.moveTo(X, py0);
            ctx.lineTo(X, py1);
            ctx.moveTo(px0, Y);
            ctx.lineTo(px1, Y);
            ctx.stroke();
            ctx.setLineDash([]);
            ctx.fillStyle = css.cursor;
            ctx.beginPath();
            ctx.arc(X, Y, 5, 0, Math.PI * 2);
            ctx.fill();
          }
          ctx.restore();

          drawAxisTitles(p.xlabel || '', p.ylabel || '', px0, py0, px1, py1, cell);
          hitPlots.push({
            kind: 'path',
            idx: cell.idx,
            r: cell.r,
            c: cell.c,
            px0, px1, py0, py1,
            xlo, xhi, ylo, yhi,
          });
        } else {
          const [tlo, thi] = tRange;
          const [ylo, yhi] = box.y;
          const tNow = times[Math.min(i, nFrames - 1)];
          drawGrid(px0, py0, px1, py1, tlo, thi, ylo, yhi, { nx: 4, ny: 4 });

          ctx.save();
          ctx.beginPath();
          ctx.rect(px0, py0, px1 - px0, py1 - py0);
          ctx.clip();

          const series = p.series || [];
          let yCursor = null;
          for (let j = 0; j < series.length; j++) {
            const s = series[j];
            const col = s.color || '#1f77b4';
            const ys = s.y || [];
            const last = Math.min(nFrames, ys.length) - 1;
            if (last < 0) continue;
            const ii = Math.min(i, last);

            ctx.strokeStyle = col;
            ctx.lineWidth = 1.75;
            ctx.globalAlpha = 0.35;
            strokePolyline(times, ys, 0, last, tlo, thi, ylo, yhi, px0, px1, py0, py1);
            ctx.globalAlpha = 1;
            ctx.lineWidth = 2.25;
            strokePolyline(times, ys, 0, ii, tlo, thi, ylo, yhi, px0, px1, py0, py1);

            const X = mapX(times[ii], tlo, thi, px0, px1);
            const Y = mapY(ys[ii], ylo, yhi, py0, py1);
            if (j === 0) {
              yCursor = Y;
            }
            ctx.fillStyle = col;
            ctx.beginPath();
            ctx.arc(X, Y, 4, 0, Math.PI * 2);
            ctx.fill();
          }

          // Fixed vertical markers (e.g. T_sync)
          const vlines = p.vlines || [];
          for (const tv of vlines) {
            const xv = mapX(tv, tlo, thi, px0, px1);
            ctx.strokeStyle = css.muted;
            ctx.lineWidth = 1;
            ctx.setLineDash([3, 3]);
            ctx.beginPath();
            ctx.moveTo(xv, py0);
            ctx.lineTo(xv, py1);
            ctx.stroke();
            ctx.setLineDash([]);
          }

          // Time cursor at the active sample (pointer→frame mapping is exact
          // for uniform grids, so this sits under the mouse while dragging).
          const xc = mapX(tNow, tlo, thi, px0, px1);
          ctx.strokeStyle = css.cursor;
          ctx.lineWidth = 1.5;
          ctx.setLineDash([4, 3]);
          ctx.beginPath();
          ctx.moveTo(xc, py0);
          ctx.lineTo(xc, py1);
          if (yCursor !== null) {
            ctx.moveTo(px0, yCursor);
            ctx.lineTo(px1, yCursor);
          }
          ctx.stroke();
          ctx.setLineDash([]);
          // Grab handle on the vertical cursor (makes drag discoverable)
          ctx.fillStyle = css.cursor;
          ctx.beginPath();
          ctx.moveTo(xc, py0);
          ctx.lineTo(xc - 5, py0 - 7);
          ctx.lineTo(xc + 5, py0 - 7);
          ctx.closePath();
          ctx.fill();
          ctx.restore();

          // Per-panel series legend (top-right inside the plot box)
          if (series.length > 0) {
            const sw = 12, gap = 4, rowH = 13;
            const pad = 4;
            let maxTw = 0;
            ctx.font = '10px system-ui, sans-serif';
            for (const s of series) {
              const nm = s.name || '';
              maxTw = Math.max(maxTw, ctx.measureText(nm).width);
            }
            const boxW = pad + sw + 4 + maxTw + pad;
            const boxH = pad + series.length * rowH + pad - 2;
            let bx = px1 - boxW - 2;
            let by = py0 + 2;
            // Keep legend inside the plot when the panel is tight
            if (bx < px0 + 2) bx = px0 + 2;
            ctx.fillStyle = css.panel;
            ctx.globalAlpha = 0.88;
            ctx.fillRect(bx, by, boxW, boxH);
            ctx.globalAlpha = 1;
            ctx.strokeStyle = css.grid;
            ctx.lineWidth = 1;
            ctx.strokeRect(bx + 0.5, by + 0.5, boxW - 1, boxH - 1);
            for (let j = 0; j < series.length; j++) {
              const s = series[j];
              const col = s.color || '#1f77b4';
              const yy = by + pad + j * rowH + 5;
              ctx.strokeStyle = col;
              ctx.lineWidth = 2.5;
              ctx.beginPath();
              ctx.moveTo(bx + pad, yy);
              ctx.lineTo(bx + pad + sw, yy);
              ctx.stroke();
              drawText(s.name || '', bx + pad + sw + 4, yy, {
                align: 'left', color: css.fg, font: '10px system-ui,sans-serif', baseline: 'middle',
              });
            }
          }

          drawAxisTitles(p.xlabel || 't [s]', p.ylabel || '', px0, py0, px1, py1, cell);
          hitPlots.push({
            kind: 'time',
            idx: cell.idx,
            r: cell.r,
            c: cell.c,
            px0, px1, py0, py1,
          });
        }
      }
    }

    let draggingCursor = false;
    let dragPanelIdx = -1; // panel that started the drag (geometry refreshed each move)

    function showFrame(i) {
      if (nFrames === 0) return;
      i = ((i % nFrames) + nFrames) % nFrames;
      frameIndex = i;
      setChrome(i);
      drawFrame(i);
    }
    function draw() { showFrame(frameIndex); }

    function pause() { playing = false; }
    function play() {
      if (playing || nFrames === 0) return;
      playing = true;
      playOriginWall = performance.now();
      playOriginFrame = frameIndex;
    }
    function togglePlayPause() { if (playing) pause(); else play(); }

    // Map event position into CSS-pixel canvas space (matches setTransform drawing).
    function canvasPos(e) {
      const r = canvas.getBoundingClientRect();
      const cssW = canvas.clientWidth || r.width || 1;
      const cssH = canvas.clientHeight || r.height || 1;
      return {
        x: (e.clientX - r.left) * (cssW / Math.max(r.width, 1e-9)),
        y: (e.clientY - r.top) * (cssH / Math.max(r.height, 1e-9)),
      };
    }

    // Left edge of plot box → frame 0, right edge → last frame (uniform grids).
    function frameFromTimeBox(px0, px1, cx) {
      if (nFrames <= 1) return 0;
      const x = Math.max(px0, Math.min(px1, cx));
      const u = (x - px0) / Math.max(px1 - px0, 1e-9);
      return Math.round(u * (nFrames - 1));
    }

    // Path scrub: map pointer → data (prefer x = leader for cam maps), pick the
    // frame whose *playhead* is closest. Continuity bias keeps multi-cycle
    // leaders from flipping between identical wraps.
    function frameFromPathBox(idx, px0, py0, px1, py1, xlo, xhi, ylo, yhi, cx, cy) {
      const p = panels[idx];
      const mx = p.path_mark_x || [];
      const my = p.path_mark_y || [];
      const hasMarks = mx.length > 0 && my.length > 0;
      const xs = hasMarks ? mx : (p.path_x || []);
      const ys = hasMarks ? my : (p.path_y || []);
      const n = Math.min(xs.length, ys.length, nFrames);
      if (n === 0) return 0;
      // Pointer → data coordinates
      const u = (Math.max(px0, Math.min(px1, cx)) - px0) / Math.max(px1 - px0, 1e-9);
      const v = 1 - (Math.max(py0, Math.min(py1, cy)) - py0) / Math.max(py1 - py0, 1e-9);
      const xT = xlo + u * (xhi - xlo);
      const yT = ylo + v * (yhi - ylo);
      const xSpan = Math.max(Math.abs(xhi - xlo), 1e-12);
      const ySpan = Math.max(Math.abs(yhi - ylo), 1e-12);
      let best = 0, bestScore = Infinity;
      for (let k = 0; k < n; k++) {
        const dx = (xs[k] - xT) / xSpan;
        const dy = (ys[k] - yT) / ySpan;
        // Weight leader (x) more — cam maps are functions of θ.
        const cont = 0.002 * (Math.abs(k - frameIndex) / Math.max(n, 1));
        const score = 4 * dx * dx + dy * dy + cont;
        if (score < bestScore) { bestScore = score; best = k; }
      }
      return best;
    }

    function hitPanel(cx, cy) {
      for (const p of hitPlots) {
        if (cx >= p.px0 && cx <= p.px1 && cy >= p.py0 && cy <= p.py1) return p;
      }
      return null;
    }

    function nearCursor(panel, cx, cy) {
      if (panel.kind === 'time') {
        const [tlo, thi] = tRange;
        const xc = mapX(times[frameIndex], tlo, thi, panel.px0, panel.px1);
        // Full-height grab strip so time cursor is easy to catch
        return Math.abs(cx - xc) <= 12;
      }
      const p = panels[panel.idx];
      const mx = p.path_mark_x || [];
      const my = p.path_mark_y || [];
      const hasMarks = mx.length > 0 && my.length > 0;
      const xs = hasMarks ? mx : (p.path_x || []);
      const ys = hasMarks ? my : (p.path_y || []);
      if (xs.length === 0) return false;
      const ii = Math.min(frameIndex, xs.length - 1, ys.length - 1);
      const X = mapX(xs[ii], panel.xlo, panel.xhi, panel.px0, panel.px1);
      const Y = mapY(ys[ii], panel.ylo, panel.yhi, panel.py0, panel.py1);
      // Vertical strip + nearby point (path scrub is mostly horizontal in m)
      return Math.abs(cx - X) <= 14 || Math.hypot(cx - X, cy - Y) <= 16;
    }

    // Fresh geometry for the panel under scrub (recomputed from current layout).
    function livePanelGeom(idx) {
      const L = layout();
      const cell = L.cells[idx];
      if (!cell) return null;
      const box = plotBox(cell);
      const p = panels[idx];
      const kind = (p.kind === 'path') ? 'path' : 'time';
      let xlo, xhi, ylo, yhi;
      if (kind === 'path' && panelBox[idx] && panelBox[idx].box) {
        [xlo, xhi, ylo, yhi] = panelBox[idx].box;
      }
      return { kind, idx, ...box, xlo, xhi, ylo, yhi };
    }

    function scrubAt(e) {
      const { x, y } = canvasPos(e);
      let idx = dragPanelIdx;
      if (idx < 0) {
        const hit = hitPanel(x, y);
        if (!hit) return false;
        idx = hit.idx;
      }
      const geom = livePanelGeom(idx);
      if (!geom) return false;
      pause();
      let frame;
      if (geom.kind === 'time') {
        // Pointer x within THIS plot box → frame. Dragging outside clamps to ends.
        frame = frameFromTimeBox(geom.px0, geom.px1, x);
      } else {
        frame = frameFromPathBox(
          geom.idx, geom.px0, geom.py0, geom.px1, geom.py1,
          geom.xlo, geom.xhi, geom.ylo, geom.yhi, x, y
        );
      }
      showFrame(frame);
      return true;
    }

    canvas.addEventListener('pointerdown', (e) => {
      if (e.button !== 0) return;
      const { x, y } = canvasPos(e);
      // Need hitPlots from a draw; ensure one exists
      if (hitPlots.length === 0) draw();
      const panel = hitPanel(x, y);
      if (!panel) return;
      draggingCursor = true;
      dragPanelIdx = panel.idx;
      canvas.setPointerCapture(e.pointerId);
      canvas.style.cursor = panel.kind === 'time' ? 'ew-resize' : 'move';
      scrubAt(e);
      e.preventDefault();
    });
    canvas.addEventListener('pointermove', (e) => {
      if (draggingCursor) {
        scrubAt(e);
        return;
      }
      const { x, y } = canvasPos(e);
      const panel = hitPanel(x, y);
      if (panel && nearCursor(panel, x, y)) {
        canvas.style.cursor = panel.kind === 'time' ? 'ew-resize' : 'move';
      } else if (panel) {
        canvas.style.cursor = 'crosshair';
      } else {
        canvas.style.cursor = '';
      }
    });
    function endDrag(e) {
      if (!draggingCursor) return;
      draggingCursor = false;
      dragPanelIdx = -1;
      try { canvas.releasePointerCapture(e.pointerId); } catch (_) {}
      canvas.style.cursor = '';
      draw();
    }
    canvas.addEventListener('pointerup', endDrag);
    canvas.addEventListener('pointercancel', endDrag);

    document.getElementById('btnPlay').onclick = play;
    document.getElementById('btnPause').onclick = pause;
    window.addEventListener('keydown', (e) => {
      if (e.code !== 'Space' && e.key !== ' ') return;
      if (e.target && (e.target.tagName === 'INPUT' || e.target.tagName === 'TEXTAREA')) return;
      e.preventDefault();
      togglePlayPause();
    });
    scrub.addEventListener('input', () => {
      pause();
      showFrame(Number(scrub.value) || 0);
    });
    window.addEventListener('resize', () => draw());

    applyTheme();
    showFrame(0);

    function loop(now) {
      if (playing && nFrames > 0) {
        const elapsedSec = (now - playOriginWall) * 0.001;
        let i = playOriginFrame + Math.floor(elapsedSec / dtSec);
        if (i >= nFrames) {
          i = i % nFrames;
          playOriginWall = now;
          playOriginFrame = i;
        }
        if (i !== frameIndex) showFrame(i);
      }
      requestAnimationFrame(loop);
    }
    requestAnimationFrame(loop);
  </script>
</body>
</html>
)";

    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    out << html;
}

// ---------------------------------------------------------------------------
// Kinematic-maps convenience wrapper (path + p/v/a per machine column)
// ---------------------------------------------------------------------------

/**
 * @brief One machine column: planar path + per-component p/v/a vs time
 *
 * Layout of the written figure:
 *
 *   columns = machines
 *   rows    = path (2D) | position | velocity | acceleration
 */
struct KinematicColumn {
    std::string                      label;
    std::vector<double>              path_x;
    std::vector<double>              path_y;
    std::string                      path_xlabel{"x"};
    std::string                      path_ylabel{"y"};
    std::vector<std::vector<double>> pos; //!< pos[component][sample]
    std::vector<std::vector<double>> vel;
    std::vector<std::vector<double>> acc;
    std::vector<std::string>         colors; //!< per component (CSS/hex)
};

/**
 * @brief Write animated kinematic-map HTML via @ref write_animated_grid
 */
inline void write_animated_kinematic_maps(
    const std::filesystem::path&        path,
    const std::string&                  title,
    const std::vector<double>&          times,
    double                              dt_s,
    const std::vector<KinematicColumn>& cols
) {
    if (times.empty() || cols.empty()) {
        return;
    }

    const int                      nrows = 4;
    const int                      ncols = static_cast<int>(cols.size());
    std::vector<Panel2d>           panels(static_cast<size_t>(nrows * ncols));
    std::vector<std::string>       col_headers(static_cast<size_t>(ncols));
    const std::vector<std::string> row_headers{
        "path (2D)",
        "position",
        "velocity",
        "acceleration",
    };

    auto default_colors = [](size_t n) {
        static const char*       palette[] = {"#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"};
        std::vector<std::string> out;
        out.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            out.emplace_back(palette[i % 5]);
        }
        return out;
    };

    auto make_time_panel = [&](const KinematicColumn&                  c,
                               const std::vector<std::vector<double>>& series,
                               const char*                             ylabel) {
        Panel2d p;
        p.kind = Panel2d::Kind::Time;
        p.xlabel = "t [s]";
        p.ylabel = ylabel;
        p.y_auto = false;
        p.y_lo = -1.18;
        p.y_hi = 1.18;
        const auto colors = c.colors.empty() ? default_colors(series.size()) : c.colors;
        for (size_t j = 0; j < series.size(); ++j) {
            Series2d s;
            s.name = fmt::format("c{}", j);
            s.color = (j < colors.size()) ? colors[j] : "#1f77b4";
            s.y = series[j];
            p.series.push_back(std::move(s));
        }
        return p;
    };

    for (int c = 0; c < ncols; ++c) {
        const auto& col = cols[static_cast<size_t>(c)];
        col_headers[static_cast<size_t>(c)] = col.label;

        Panel2d path_panel;
        path_panel.kind = Panel2d::Kind::Path;
        path_panel.xlabel = col.path_xlabel;
        path_panel.ylabel = col.path_ylabel;
        path_panel.path_x = col.path_x;
        path_panel.path_y = col.path_y;
        path_panel.path_color = "#7c3aed";
        panels[static_cast<size_t>(0 * ncols + c)] = std::move(path_panel);

        panels[static_cast<size_t>(1 * ncols + c)] = make_time_panel(col, col.pos, "position (norm.)");
        panels[static_cast<size_t>(2 * ncols + c)] = make_time_panel(col, col.vel, "velocity (norm.)");
        panels[static_cast<size_t>(3 * ncols + c)] = make_time_panel(col, col.acc, "acceleration (norm.)");
    }

    write_animated_grid(
        path, title, times, dt_s, nrows, ncols, panels, col_headers, row_headers
    );
}

} // namespace damp::examples_plot2d
