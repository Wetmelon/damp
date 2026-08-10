// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file examples/animate_3d.hpp
 * @brief Host-only 3D animation writer using **three.js** + OrbitControls.
 *
 * Plotly scatter3d (WebGL) is fine for static 3D, but continuous restyle while
 * orbiting fights the camera (snap / freeze). three.js updates line buffers
 * independently of OrbitControls, so play + pan/rotate work together.
 *
 * Not part of the embeddable library. CDN three.js (no submodule).
 */

#include <algorithm>
#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "damp/kinematics/pose.hpp"
#include "damp/math/geometry.hpp"
#include "damp/matrix/colvec.hpp"
#include "fmt/core.h"
#include "nlohmann/json.hpp"

namespace damp::examples_plot3d {

using json = nlohmann::json;

struct Polyline3 {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
};

struct TraceStyle {
    std::string name;
    std::string color{"#1f77b4"}; //!< CSS/hex color for the line
    std::string mode{"lines"};    //!< "lines", "lines+markers", "markers" (markers → points)
    int         line_width{3};    //!< three.js linewidth is often ignored by WebGL; kept for API
    int         marker_size{4};
    /// "solid" or "dot" / "dash" / "longdash" → dashed LineDashedMaterial
    std::string dash{"solid"};
    /// If true, geometry is drawn once and never updated per frame.
    bool static_geometry{false};
};

struct SceneBox {
    std::array<double, 2> x{-1.0, 1.0};
    std::array<double, 2> y{-1.0, 1.0};
    std::array<double, 2> z{-1.0, 1.0};
    double                eye_x{1.55};
    double                eye_y{1.35};
    double                eye_z{1.15};
    double                center_x{0.0};
    double                center_y{0.0};
    double                center_z{0.0};
};

/**
 * @brief Axis-aligned scene + orbit camera from world-space bounds (Z-up)
 *
 * Pads the AABB, puts the OrbitControls target at the box centre, and places
 * the eye on a diagonal so the whole volume fits in a ~50° FOV. Prefer this over
 * hand-picked numbers when paths leave the origin (INS demos, long trajectories).
 */
[[nodiscard]] inline SceneBox fit_scene_box(
    double xmin,
    double xmax,
    double ymin,
    double ymax,
    double zmin,
    double zmax,
    double pad_frac = 0.15
) {
    if (xmax < xmin) {
        std::swap(xmax, xmin);
    }
    if (ymax < ymin) {
        std::swap(ymax, ymin);
    }
    if (zmax < zmin) {
        std::swap(zmax, zmin);
    }
    const double cx = 0.5 * (xmin + xmax);
    const double cy = 0.5 * (ymin + ymax);
    const double cz = 0.5 * (zmin + zmax);
    const double sx = xmax - xmin;
    const double sy = ymax - ymin;
    const double sz = zmax - zmin;
    const double span = std::max({sx, sy, sz, 1.0});
    const double pad = pad_frac * span + 0.5;

    SceneBox b;
    b.x = {xmin - pad, xmax + pad};
    b.y = {ymin - pad, ymax + pad};
    b.z = {zmin - pad, zmax + pad};
    b.center_x = cx;
    b.center_y = cy;
    b.center_z = cz;
    // Distance ~ span / tan(fov/2); 1.7× is comfortable for orbit start.
    const double dist = 1.7 * span;
    b.eye_x = cx + 0.55 * dist;
    b.eye_y = cy - 0.70 * dist;
    b.eye_z = cz + 0.45 * dist;
    return b;
}

/** Ground rectangle on z = 0 covering an XY range (closed loop). */
[[nodiscard]] inline Polyline3 ground_rect(double xmin, double xmax, double ymin, double ymax) {
    return Polyline3{
        {xmin, xmax, xmax, xmin, xmin},
        {ymin, ymin, ymax, ymax, ymin},
        {0.0, 0.0, 0.0, 0.0, 0.0},
    };
}

[[nodiscard]] inline std::string json_js_literal(const json& j) {
    std::string s = j.dump();
    size_t      pos = 0;
    while ((pos = s.find('<', pos)) != std::string::npos) {
        s.replace(pos, 1, "\\u003c");
        pos += 6;
    }
    return s;
}

[[nodiscard]] inline bool is_dashed(const std::string& dash) {
    return dash != "solid" && !dash.empty();
}

/**
 * @param path    e.g. plots/foo.html
 * @param styles  one style per trace
 * @param frames  frames[t][trace] = polyline
 * @param times_s slider / label times
 * @param dt_s    sim step used for wall-clock real-time play
 * @param title   page title
 * @param box     axis-ish extents + initial camera (Z-up scene)
 */
inline void write_animated_scatter3d(
    const std::filesystem::path&               path,
    const std::vector<TraceStyle>&             styles,
    const std::vector<std::vector<Polyline3>>& frames,
    const std::vector<double>&                 times_s,
    double                                     dt_s,
    const std::string&                         title,
    const SceneBox&                            box = {}
) {
    namespace fs = std::filesystem;
    if (frames.empty() || styles.empty() || frames[0].size() != styles.size()) {
        return;
    }
    if (!path.parent_path().empty()) {
        fs::create_directories(path.parent_path());
    }

    const size_t n_tr = styles.size();
    const size_t n_fr = frames.size();

    // Per-trace metadata + initial geometry; per-frame dynamic xyz only.
    json                traces_meta = json::array();
    json                dyn_idx = json::array();
    std::vector<size_t> dyn_map;
    for (size_t t = 0; t < n_tr; ++t) {
        const auto& st = styles[t];
        const auto& pl = frames[0][t];
        traces_meta.push_back({
            {"name", st.name},
            {"color", st.color},
            {"dashed", is_dashed(st.dash)},
            {"static", st.static_geometry},
            {"markers", st.mode.find("markers") != std::string::npos && st.mode.find("lines") == std::string::npos},
            {"x", pl.x},
            {"y", pl.y},
            {"z", pl.z},
        });
        if (!st.static_geometry) {
            dyn_idx.push_back(t);
            dyn_map.push_back(t);
        }
    }

    json frames_xyz = json::array();
    for (size_t i = 0; i < n_fr; ++i) {
        json one = json::array();
        for (size_t t : dyn_map) {
            const auto& pl = frames[i][t];
            one.push_back({{"x", pl.x}, {"y", pl.y}, {"z", pl.z}});
        }
        frames_xyz.push_back(std::move(one));
    }

    json scene = {
        {"x", json::array({box.x[0], box.x[1]})},
        {"y", json::array({box.y[0], box.y[1]})},
        {"z", json::array({box.z[0], box.z[1]})},
        {"eye", {{"x", box.eye_x}, {"y", box.eye_y}, {"z", box.eye_z}}},
        {"center", {{"x", box.center_x}, {"y", box.center_y}, {"z", box.center_z}}},
    };

    std::string html = std::string(R"(<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>)") + title + R"(</title>
  <style>
    :root, [data-theme="light"] {
      --bg: #f4f4f5;
      --fg: #18181b;
      --muted: #52525b;
      --btn-bg: #fff;
      --btn-border: #d4d4d8;
      --btn-hover: #f4f4f5;
      --legend-bg: rgba(255, 255, 255, 0.88);
      --clear: 0xf4f4f5;
      --grid-a: 0xaaaaaa;
      --grid-b: 0xcccccc;
    }
    [data-theme="dark"] {
      --bg: #111113;
      --fg: #e4e4e7;
      --muted: #a1a1aa;
      --btn-bg: #27272a;
      --btn-border: #52525b;
      --btn-hover: #3f3f46;
      --legend-bg: rgba(0, 0, 0, 0.62);
      --clear: 0x111113;
      --grid-a: 0x333333;
      --grid-b: 0x222222;
    }
    html, body {
      margin: 0; height: 100%; overflow: hidden;
      background: var(--bg); color: var(--fg);
      font-family: system-ui, sans-serif;
    }
    #bar {
      position: absolute; z-index: 10; top: 8px; left: 12px; right: 12px;
      display: flex; flex-wrap: wrap; gap: 10px; align-items: center;
      pointer-events: none;
    }
    #bar > * { pointer-events: auto; }
    #bar button {
      cursor: pointer; padding: 6px 14px; border-radius: 4px;
      border: 1px solid var(--btn-border); background: var(--btn-bg); color: var(--fg);
    }
    #bar button:hover { background: var(--btn-hover); }
    #tlabel { font-size: 13px; min-width: 7em; }
    #hint { font-size: 12px; color: var(--muted); }
    #scrubWrap {
      position: absolute; z-index: 10; left: 12px; right: 12px; bottom: 8px;
      height: 28px; display: flex; align-items: center; gap: 10px;
      pointer-events: none;
    }
    #scrubWrap > * { pointer-events: auto; }
    #scrub { flex: 1; height: 28px; cursor: pointer; }
    #c { display: block; width: 100%; height: 100%; }
    #legend {
      position: absolute; z-index: 10; top: 48px; right: 12px;
      background: var(--legend-bg); padding: 8px 10px; border-radius: 6px;
      font-size: 12px; max-width: 220px; pointer-events: none;
      border: 1px solid var(--btn-border);
    }
    #legend div { margin: 2px 0; display: flex; align-items: center; gap: 6px; }
    #legend i { width: 18px; height: 0; border-top-width: 3px; border-top-style: solid; display: inline-block; }
  </style>
  <script type="importmap">
  {
    "imports": {
      "three": "https://unpkg.com/three@0.160.0/build/three.module.js",
      "three/addons/": "https://unpkg.com/three@0.160.0/examples/jsm/"
    }
  }
  </script>
</head>
<body>
  <div id="bar">
    <button type="button" id="btnPlay">Play</button>
    <button type="button" id="btnPause">Pause</button>
    <button type="button" id="btnTheme" title="Theme: system / light / dark">Theme: system</button>
    <span id="tlabel">t = 0.00 s</span>
    <span id="hint">Space = play/pause · left-drag orbit · right-drag pan · scroll zoom</span>
  </div>
  <div id="legend"></div>
  <canvas id="c"></canvas>
  <div id="scrubWrap">
    <input type="range" id="scrub" min="0" max="0" value="0" step="1" />
  </div>
  <script type="module">
    import * as THREE from 'three';
    import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

    const traces = )" + json_js_literal(traces_meta)
                     + R"(;
    const framesXyz = )"
                     + json_js_literal(frames_xyz) + R"(;
    const dynIdx = )" + json_js_literal(dyn_idx)
                     + R"(;
    const times = )" + json_js_literal(json(times_s))
                     + R"(;
    const dtSec = )" + std::to_string(dt_s)
                     + R"(;
    const sceneBox = )"
                     + json_js_literal(scene) + R"(;
    const nFrames = framesXyz.length;

    const canvas = document.getElementById('c');
    const scrub = document.getElementById('scrub');
    const btnTheme = document.getElementById('btnTheme');
    scrub.max = Math.max(0, nFrames - 1);

    // Theme: system (default) → light → dark → system …
    const THEME_KEY = 'damp-anim-3d-theme';
    let themePref = localStorage.getItem(THEME_KEY) || 'system'; // system | light | dark
    let gridHelper = null;

    function resolvedTheme() {
      if (themePref === 'light' || themePref === 'dark') return themePref;
      return window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
    }

    function themeLabel() {
      return 'Theme: ' + themePref;
    }

    function cssClearColor() {
      // Match --clear hex tokens above (0xrrggbb as numbers in CSS vars is awkward;
      // hard-code the same pair used in :root / [data-theme=dark]).
      return resolvedTheme() === 'dark' ? 0x111113 : 0xf4f4f5;
    }

    function applyTheme() {
      const r = resolvedTheme();
      document.documentElement.setAttribute('data-theme', r);
      btnTheme.textContent = themeLabel();
      if (renderer) {
        renderer.setClearColor(cssClearColor(), 1);
      }
      if (gridHelper) {
        scene3.remove(gridHelper);
        gridHelper.geometry.dispose();
        gridHelper.material.dispose?.();
        // GridHelper has materials array
        if (Array.isArray(gridHelper.material)) {
          gridHelper.material.forEach((m) => m.dispose && m.dispose());
        }
        gridHelper = makeGrid();
        scene3.add(gridHelper);
      }
    }

    // Z-up scene (ENU / NED-friendly), not three.js default Y-up.
    const renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
    renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));

    const scene3 = new THREE.Scene();
    const spanX = Math.abs(sceneBox.x[1] - sceneBox.x[0]);
    const spanY = Math.abs(sceneBox.y[1] - sceneBox.y[0]);
    const spanZ = Math.abs(sceneBox.z[1] - sceneBox.z[0]);
    const span = Math.max(spanX, spanY, spanZ, 1);
    // Near/far scale with the scene so long INS paths stay pickable (not clipped,
    // not stuck with a near plane that eats the vehicle when zoomed in).
    const camera = new THREE.PerspectiveCamera(50, 1, Math.max(span * 1e-4, 0.01), span * 80);
    camera.up.set(0, 0, 1);
    camera.position.set(sceneBox.eye.x, sceneBox.eye.y, sceneBox.eye.z);

    const controls = new OrbitControls(camera, canvas);
    controls.target.set(sceneBox.center.x, sceneBox.center.y, sceneBox.center.z);
    controls.maxDistance = span * 40;
    controls.enableDamping = true;
    controls.dampingFactor = 0.08;
    controls.screenSpacePanning = true;
    controls.update();

    // Soft light so lines stay visible
    scene3.add(new THREE.AmbientLight(0xffffff, 0.85));

    function makeGrid() {
      const gx0 = sceneBox.x[0], gx1 = sceneBox.x[1];
      const gy0 = sceneBox.y[0], gy1 = sceneBox.y[1];
      const dark = resolvedTheme() === 'dark';
      const grid = new THREE.GridHelper(
        Math.max(gx1 - gx0, gy1 - gy0), 10,
        dark ? 0x444444 : 0x999999,
        dark ? 0x2a2a2a : 0xcccccc
      );
      grid.rotation.x = Math.PI / 2; // lie in XY (Z-up)
      grid.position.z = sceneBox.z[0];
      return grid;
    }

    applyTheme();

    function flatPositions(x, y, z) {
      const n = Math.min(x.length, y.length, z.length);
      const arr = new Float32Array(n * 3);
      for (let i = 0; i < n; i++) {
        arr[i * 3] = x[i];
        arr[i * 3 + 1] = y[i];
        arr[i * 3 + 2] = z[i];
      }
      return arr;
    }

    function makeObject(tr) {
      const pos = flatPositions(tr.x, tr.y, tr.z);
      const geo = new THREE.BufferGeometry();
      geo.setAttribute('position', new THREE.BufferAttribute(pos, 3));
      if (tr.markers) {
        const mat = new THREE.PointsMaterial({ color: tr.color, size: 0.06, sizeAttenuation: true });
        return new THREE.Points(geo, mat);
      }
      const mat = tr.dashed
        ? new THREE.LineDashedMaterial({ color: tr.color, dashSize: 0.06, gapSize: 0.04, linewidth: 2 })
        : new THREE.LineBasicMaterial({ color: tr.color });
      const line = new THREE.Line(geo, mat);
      if (tr.dashed) {
        line.computeLineDistances();
      }
      return line;
    }

    function setPositions(obj, x, y, z) {
      const pos = flatPositions(x, y, z);
      const attr = obj.geometry.getAttribute('position');
      if (!attr || attr.count !== pos.length / 3) {
        obj.geometry.setAttribute('position', new THREE.BufferAttribute(pos, 3));
      } else {
        attr.array.set(pos);
        attr.needsUpdate = true;
      }
      obj.geometry.computeBoundingSphere();
      if (obj.material && obj.material.isLineDashedMaterial) {
        obj.computeLineDistances();
      }
    }

    const objects = traces.map(makeObject);
    objects.forEach((o) => scene3.add(o));

    // Legend
    const leg = document.getElementById('legend');
    traces.forEach((tr) => {
      const row = document.createElement('div');
      const i = document.createElement('i');
      i.style.borderTopColor = tr.color;
      i.style.borderTopStyle = tr.dashed ? 'dashed' : 'solid';
      row.appendChild(i);
      row.appendChild(document.createTextNode(tr.name));
      leg.appendChild(row);
    });

    gridHelper = makeGrid();
    scene3.add(gridHelper);

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

    function tText(i) {
      const t = (times && times.length > i) ? times[i] : (i * dtSec);
      return 't = ' + Number(t).toFixed(2) + ' s';
    }

    function setChrome(i) {
      document.getElementById('tlabel').textContent = tText(i);
      if (document.activeElement !== scrub) {
        scrub.value = String(i);
      }
    }

    function showFrame(i) {
      if (nFrames === 0) return;
      i = ((i % nFrames) + nFrames) % nFrames;
      frameIndex = i;
      setChrome(i);
      const fr = framesXyz[i];
      for (let k = 0; k < dynIdx.length; k++) {
        const ti = dynIdx[k];
        setPositions(objects[ti], fr[k].x, fr[k].y, fr[k].z);
      }
    }

    function pause() {
      playing = false;
    }

    function play() {
      if (playing || nFrames === 0) return;
      playing = true;
      playOriginWall = performance.now();
      playOriginFrame = frameIndex;
    }

    function togglePlayPause() {
      if (playing) pause();
      else play();
    }

    function resize() {
      const w = canvas.clientWidth || window.innerWidth;
      const h = canvas.clientHeight || window.innerHeight;
      renderer.setSize(w, h, false);
      camera.aspect = w / Math.max(h, 1);
      camera.updateProjectionMatrix();
    }
    window.addEventListener('resize', resize);
    resize();

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
        if (i !== frameIndex) {
          showFrame(i);
        }
      }
      controls.update();
      renderer.render(scene3, camera);
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

/** Body-axis triad at a pose (body → parent via pose). */
[[nodiscard]] inline std::array<Polyline3, 3>
body_axes(const damp::Pose<double>& pose, double len) {
    std::array<Polyline3, 3> ax{};
    const damp::Vec3<double> p0 = static_cast<const damp::Vec3<double>&>(pose.translation);
    for (int a = 0; a < 3; ++a) {
        damp::Vec3<double> e{};
        e[static_cast<size_t>(a)] = len;
        const damp::Vec3<double> p1 = pose.transform_point(e);
        ax[static_cast<size_t>(a)] = Polyline3{{p0[0], p1[0]}, {p0[1], p1[1]}, {p0[2], p1[2]}};
    }
    return ax;
}

[[nodiscard]] inline std::array<Polyline3, 3>
body_axes_at(const damp::Quaternion<double>& q, const damp::Vec3<double>& origin, double len) {
    damp::Pose<double> pose;
    pose.orientation = q;
    pose.translation = damp::Translation3<double>(origin);
    return body_axes(pose, len);
}

[[nodiscard]] inline std::array<Polyline3, 3> world_axes(double len) {
    return {
        Polyline3{{0.0, len}, {0.0, 0.0}, {0.0, 0.0}},
        Polyline3{{0.0, 0.0}, {0.0, len}, {0.0, 0.0}},
        Polyline3{{0.0, 0.0}, {0.0, 0.0}, {0.0, len}},
    };
}

} // namespace damp::examples_plot3d
