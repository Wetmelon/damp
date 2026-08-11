// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file damp/toolbox/lookup.hpp
 * @brief Breakpoint lookup tables with interpolation and online adaptation.
 *
 * Fixed-size, `constexpr`, allocation-free interpolation tables for the
 * non-affine "raw → engineering" mappings that a AffineCal can't capture:
 * sensor linearization, gain-scheduling maps, fan / efficiency / torque-speed
 * curves, motor (id,iq) maps, hydraulic 3-parameter efficiency maps. 1-D
 * (linear / cubic / nearest), 2-D (bilinear), and 3-D (trilinear) over
 * monotonic breakpoints; optional adaptive cell updates for online maps.
 *
 * Breakpoints and grids are stored as the library's own ColVec / Matrix
 * types (3-D values as a flat damp::array) rather than raw heap arrays, so a
 * table composes with the linear-algebra core.
 *
 * @see toolbox/scaling.hpp for the affine / polynomial conversions.
 * @see roadmap #30 (cubic, 2-D extrapolate, adaptive tables, Lut3D)
 */

#include <cstddef>

#include "damp/backend.hpp"         // damp::clamp, damp::min, damp::max
#include "damp/matrix/matrix.hpp"   // Matrix, ColVec
#include "damp/toolbox/scaling.hpp" // lerp, inverse_lerp

namespace damp {

/// Out-of-range behaviour for a @ref Lut1D / @ref Lut2D / @ref Lut3D query beyond its breakpoints.
enum class Extrapolation {
    Clamp, ///< Hold the nearest endpoint value.
    Linear ///< Continue the slope of the nearest end segment.
};

/// In-range interpolant for @ref Lut1D (2-D stays bilinear; use @ref SplineSurface for C¹ grids).
enum class Interpolation {
    Linear, ///< Piecewise-linear between breakpoints.
    Cubic,  ///< Non-uniform Catmull-Rom through the nodes (C¹, local; uses @ref catmull_1d).
};

/// How @ref AdaptiveLut1D / @ref AdaptiveLut2D treat samples outside the breakpoint span.
enum class AdaptOutOfRange {
    Ignore,      ///< Drop the sample (no cell update).
    NearestCell, ///< Update the nearest endpoint cell.
};

/**
 * @brief Index of the interpolation segment containing @p x.
 *
 * For strictly increasing breakpoints @p xs, returns the largest `i` with
 * `xs[i] <= x`, clamped to `[0, N-2]` so `i` and `i+1` always bracket a valid
 * segment. O(log N) binary search.
 *
 * @tparam N Number of breakpoints (must be ≥ 2).
 */
template<size_t N, typename T>
[[nodiscard]] constexpr size_t lut_segment(const ColVec<N, T>& xs, T x) {
    static_assert(N >= 2, "lut_segment needs at least two breakpoints");
    if (x <= xs[1]) {
        return 0;
    }
    if (x >= xs[N - 2]) {
        return N - 2;
    }
    size_t lo = 1;
    size_t hi = N - 2;
    // lo + 1 < hi: for N == 2, hi = 0 < lo, and the size_t
    // subtraction would underflow into a phantom in-loop index (GCC -Warray-bounds).
    while (lo + 1 < hi) {
        const size_t mid = lo + ((hi - lo) / 2);
        if (xs[mid] <= x) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return lo;
}

// catmull_1d is defined below (after Lut2D); Lut1D::Cubic calls it — declare first.
template<size_t N, typename T, typename Yf>
[[nodiscard]] constexpr T catmull_1d(const ColVec<N, T>& xs, Yf y, T q);

/**
 * @brief 1-D interpolating lookup table over monotonic breakpoints.
 *
 * Stores `N` strictly-increasing breakpoints @ref xs and their values @ref ys
 * as ColVec. `operator()` does linear interpolation; @ref nearest does
 * nearest-neighbour. Queries outside `[xs[0], xs[N-1]]` follow the @ref oob
 * policy (clamp or linear extrapolation).
 *
 * @code
 * // NTC: resistance (kΩ) -> temperature (°C), falling curve.
 * constexpr Lut1D<3, float> ntc{{32.6f, 10.0f, 3.6f}, {0.0f, 25.0f, 50.0f}};
 * float t = ntc(measured_kohm);
 * @endcode
 *
 * @tparam N Number of breakpoints (must be ≥ 1).
 * @tparam T Scalar type (default: float)
 */
template<size_t N, typename T = float>
struct Lut1D {
    static_assert(N >= 1, "Lut1D needs at least one breakpoint");

    ColVec<N, T>  xs{};                          ///< Strictly increasing breakpoints.
    ColVec<N, T>  ys{};                          ///< Value at each breakpoint.
    Extrapolation oob{Extrapolation::Clamp};     ///< Out-of-range policy.
    Interpolation interp{Interpolation::Linear}; ///< In-range interpolant.

    /// Interpolated value at @p x (mode @ref interp; extrapolation per @ref oob).
    [[nodiscard]] constexpr T operator()(T x) const {
        if constexpr (N == 1) {
            return ys[0];
        } else {
            if (x <= xs[0] && oob == Extrapolation::Clamp) {
                return ys[0];
            }
            if (x >= xs[N - 1] && oob == Extrapolation::Clamp) {
                return ys[N - 1];
            }
            if (interp == Interpolation::Cubic) {
                // Catmull-Rom is clamp-only in catmull_1d; Linear oob uses end-segment slope.
                if (oob == Extrapolation::Linear && x < xs[0]) {
                    const size_t i = 0;
                    return lerp(ys[i], ys[i + 1], inverse_lerp(xs[i], xs[i + 1], x));
                }
                if (oob == Extrapolation::Linear && x > xs[N - 1]) {
                    const size_t i = N - 2;
                    return lerp(ys[i], ys[i + 1], inverse_lerp(xs[i], xs[i + 1], x));
                }
                return catmull_1d(xs, [&](size_t i) { return ys[i]; }, x);
            }
            const size_t i = lut_segment(xs, x);
            return lerp(ys[i], ys[i + 1], inverse_lerp(xs[i], xs[i + 1], x));
        }
    }

    /// Nearest-neighbour value at @p x (always clamps to the breakpoint range).
    [[nodiscard]] constexpr T nearest(T x) const {
        if constexpr (N == 1) {
            return ys[0];
        } else {
            if (x <= xs[0]) {
                return ys[0];
            }
            if (x >= xs[N - 1]) {
                return ys[N - 1];
            }
            const size_t i = lut_segment(xs, x);
            return (x - xs[i] <= xs[i + 1] - x) ? ys[i] : ys[i + 1];
        }
    }

    /// Index of the nearest breakpoint (for adaptive cell updates).
    [[nodiscard]] constexpr size_t nearest_index(T x) const {
        if constexpr (N == 1) {
            return 0;
        } else {
            if (x <= xs[0]) {
                return 0;
            }
            if (x >= xs[N - 1]) {
                return N - 1;
            }
            const size_t i = lut_segment(xs, x);
            return (x - xs[i] <= xs[i + 1] - x) ? i : i + 1;
        }
    }
};

/**
 * @brief 2-D bilinear interpolating lookup table over a regular grid.
 *
 * Row breakpoints @ref rows (length `R`) and column breakpoints @ref cols
 * (length `C`) index a grid of values @ref z (`z(r, c)`), stored as a
 * Matrix. `operator()(r, c)` bilinearly interpolates. Queries outside the
 * grid follow @ref oob per axis (clamp or linear end-segment extrapolation) —
 * the same policy as @ref Lut1D / @ref Lut3D.
 *
 * @code
 * // Efficiency map vs (speed, torque):
 * Lut2D<2, 2, float> eff{{0.0f, 1.0f}, {0.0f, 1.0f},
 *                        {{0.80f, 0.85f}, {0.88f, 0.92f}}};
 * float e = eff(speed, torque);
 * @endcode
 *
 * @tparam R Number of row breakpoints (must be ≥ 1).
 * @tparam C Number of column breakpoints (must be ≥ 1).
 * @tparam T Scalar type (default: float)
 *
 * @see Lut1D, Lut3D for 1-D / 3-D partners with the same @ref Extrapolation policy
 */
template<size_t R, size_t C, typename T = float>
struct Lut2D {
    static_assert(R >= 1 && C >= 1, "Lut2D needs at least one breakpoint per axis");

    ColVec<R, T>    rows{};                    ///< Strictly increasing row breakpoints.
    ColVec<C, T>    cols{};                    ///< Strictly increasing column breakpoints.
    Matrix<R, C, T> z{};                       ///< Grid values, z(row, col).
    Extrapolation   oob{Extrapolation::Clamp}; ///< Out-of-range policy (applied per axis).

    /// Bilinearly interpolated value at (@p r, @p c); out-of-grid per @ref oob.
    [[nodiscard]] constexpr T operator()(T r, T c) const {
        const T rc = (oob == Extrapolation::Clamp) ? damp::clamp(r, rows[0], rows[R - 1]) : r;
        const T cc = (oob == Extrapolation::Clamp) ? damp::clamp(c, cols[0], cols[C - 1]) : c;

        const size_t i = (R == 1) ? 0 : lut_segment(rows, damp::clamp(rc, rows[0], rows[R - 1]));
        const size_t j = (C == 1) ? 0 : lut_segment(cols, damp::clamp(cc, cols[0], cols[C - 1]));

        // inverse_lerp may leave [0,1] when oob == Linear and the query is outside.
        const T tr = (R == 1) ? T{0} : inverse_lerp(rows[i], rows[i + 1], rc);
        const T tc = (C == 1) ? T{0} : inverse_lerp(cols[j], cols[j + 1], cc);

        const size_t i1 = (R == 1) ? 0 : i + 1;
        const size_t j1 = (C == 1) ? 0 : j + 1;

        const T top = lerp(z(i, j), z(i, j1), tc);
        const T bot = lerp(z(i1, j), z(i1, j1), tc);
        return lerp(top, bot, tr);
    }
};

/**
 * @brief 3-D trilinear interpolating lookup table over a regular grid.
 *
 * Axis breakpoints @ref xs / @ref ys / @ref zs (lengths `NX` / `NY` / `NZ`) index
 * a dense value tensor @ref v with layout
 * @f$ v(i,j,k) \mapsto \mathrm{index}\, i\cdot(N_Y N_Z) + j\cdot N_Z + k @f$
 * (x slowest, z fastest). `operator()(x, y, z)` trilinearly interpolates using
 * the same @ref inverse_lerp / @ref lerp path as @ref Lut2D. Queries outside the
 * grid follow @ref oob per axis (clamp or linear end-segment extrapolation).
 *
 * Typical use: a three-parameter map such as hydraulic pump volumetric
 * efficiency @f$ \eta_v = f(V_g, \omega, p) @f$. Axis roles are caller-defined;
 * a hydraulic-shaped layout that matches the literature is:
 * - @ref xs = displacement fraction @f$ V_g @f$ (or swash fraction) ∈ [0, 1]
 * - @ref ys = shaft speed @f$ \omega @f$ [rad/s] (or rpm if the grid is in rpm)
 * - @ref zs = differential pressure @f$ p @f$ [Pa] (or bar if the grid is in bar)
 * - @ref v  = @f$ \eta_v @f$ (or another scalar map value) at each grid node
 *
 * Keep units consistent between the breakpoints and the query arguments; the
 * table does not convert units.
 *
 * @code
 * // η_v vs (V_g ∈ [0,1], ω rad/s, p bar) — hydraulic-shaped axes.
 * constexpr Lut3D<2, 2, 2, float> eta{
 *     {0.0f, 1.0f},       // xs: V_g
 *     {0.0f, 200.0f},     // ys: ω [rad/s]
 *     {0.0f, 300.0f},     // zs: p [bar]
 *     // v(i,j,k): flat, k (pressure) fastest
 *     {0.90f, 0.88f,  // i=0,j=0, k=0..1
 *      0.92f, 0.90f,  // i=0,j=1
 *      0.93f, 0.91f,  // i=1,j=0
 *      0.95f, 0.93f}  // i=1,j=1
 * };
 * float e = eta(0.5f, 100.0f, 150.0f); // V_g, ω [rad/s], p [bar]
 * @endcode
 *
 * @tparam NX Number of x-axis breakpoints (must be ≥ 1).
 * @tparam NY Number of y-axis breakpoints (must be ≥ 1).
 * @tparam NZ Number of z-axis breakpoints (must be ≥ 1).
 * @tparam T  Scalar type (default: float)
 *
 * @see Lut2D for the 2-D bilinear analog
 * @see roadmap #30
 */
template<size_t NX, size_t NY, size_t NZ, typename T = float>
struct Lut3D {
    static_assert(NX >= 1 && NY >= 1 && NZ >= 1, "Lut3D needs at least one breakpoint per axis");

    ColVec<NX, T>                xs{};                      ///< Strictly increasing x breakpoints.
    ColVec<NY, T>                ys{};                      ///< Strictly increasing y breakpoints.
    ColVec<NZ, T>                zs{};                      ///< Strictly increasing z breakpoints.
    damp::array<T, NX * NY * NZ> v{};                       ///< Dense values; see @ref at for (i,j,k).
    Extrapolation                oob{Extrapolation::Clamp}; ///< Out-of-range policy (applied per axis).

    /// Flat index for grid node (@p i, @p j, @p k); x slowest, z fastest.
    [[nodiscard]] static constexpr size_t index(size_t i, size_t j, size_t k) {
        return (i * (NY * NZ)) + (j * NZ) + k;
    }

    /// Mutable value at grid node (@p i, @p j, @p k).
    [[nodiscard]] constexpr T& at(size_t i, size_t j, size_t k) { return v[index(i, j, k)]; }

    /// Const value at grid node (@p i, @p j, @p k).
    [[nodiscard]] constexpr const T& at(size_t i, size_t j, size_t k) const {
        return v[index(i, j, k)];
    }

    /// Trilinearly interpolated value at (@p x, @p y, @p z); out-of-grid per @ref oob.
    [[nodiscard]] constexpr T operator()(T x, T y, T z) const {
        const T xc = (oob == Extrapolation::Clamp) ? damp::clamp(x, xs[0], xs[NX - 1]) : x;
        const T yc = (oob == Extrapolation::Clamp) ? damp::clamp(y, ys[0], ys[NY - 1]) : y;
        const T zc = (oob == Extrapolation::Clamp) ? damp::clamp(z, zs[0], zs[NZ - 1]) : z;

        const size_t i = (NX == 1) ? 0 : lut_segment(xs, damp::clamp(xc, xs[0], xs[NX - 1]));
        const size_t j = (NY == 1) ? 0 : lut_segment(ys, damp::clamp(yc, ys[0], ys[NY - 1]));
        const size_t k = (NZ == 1) ? 0 : lut_segment(zs, damp::clamp(zc, zs[0], zs[NZ - 1]));

        // inverse_lerp may leave [0,1] when oob == Linear and the query is outside.
        const T tx = (NX == 1) ? T{0} : inverse_lerp(xs[i], xs[i + 1], xc);
        const T ty = (NY == 1) ? T{0} : inverse_lerp(ys[j], ys[j + 1], yc);
        const T tz = (NZ == 1) ? T{0} : inverse_lerp(zs[k], zs[k + 1], zc);

        const size_t i1 = (NX == 1) ? 0 : i + 1;
        const size_t j1 = (NY == 1) ? 0 : j + 1;
        const size_t k1 = (NZ == 1) ? 0 : k + 1;

        const T c00 = lerp(at(i, j, k), at(i, j, k1), tz);
        const T c01 = lerp(at(i, j1, k), at(i, j1, k1), tz);
        const T c10 = lerp(at(i1, j, k), at(i1, j, k1), tz);
        const T c11 = lerp(at(i1, j1, k), at(i1, j1, k1), tz);
        const T c0 = lerp(c00, c01, ty);
        const T c1 = lerp(c10, c11, ty);
        return lerp(c0, c1, tx);
    }
};

/**
 * @brief 1-D non-uniform Catmull-Rom (cubic Hermite) evaluation.
 *
 * Interpolates the values @p y (indexed `y(i)`) at strictly-increasing
 * breakpoints @p xs, querying at @p q (clamped to the breakpoint range). Uses
 * local finite-difference tangents scaled for non-uniform spacing, so the curve
 * passes through every node (C¹, interpolating) with no linear solve — unlike a
 * natural cubic spline. Two breakpoints degrade to linear; one returns the value.
 * Used by @ref Lut1D with @ref Interpolation::Cubic.
 *
 * @tparam N  Number of breakpoints.
 * @tparam Yf Callable `size_t -> T` giving the value at breakpoint i.
 */
template<size_t N, typename T, typename Yf>
[[nodiscard]] constexpr T catmull_1d(const ColVec<N, T>& xs, Yf y, T q) {
    if constexpr (N == 1) {
        return y(0);
    } else {
        const T      qc = damp::clamp(q, xs[0], xs[N - 1]);
        const size_t i = lut_segment(xs, qc); // i in [0, N-2]
        const T      h = xs[i + 1] - xs[i];
        const T      t = (h > T{0}) ? (qc - xs[i]) / h : T{0};

        const T yi = y(i);
        const T yi1 = y(i + 1);
        // Catmull-Rom tangents (one-sided at the ends), scaled for spacing.
        // Zero-width segments (non-strict breakpoints) yield zero tangent rather
        // than Inf/NaN — breakpoints should still be strictly increasing.
        T mi{};
        if (i == 0) {
            mi = (h > T{0}) ? (yi1 - yi) / h : T{0};
        } else {
            const T span = xs[i + 1] - xs[i - 1];
            mi = (span > T{0}) ? (yi1 - y(i - 1)) / span : T{0};
        }
        T mi1{};
        if (i + 2 >= N) {
            mi1 = (h > T{0}) ? (yi1 - yi) / h : T{0};
        } else {
            const T span = xs[i + 2] - xs[i];
            mi1 = (span > T{0}) ? (y(i + 2) - yi) / span : T{0};
        }

        const T t2 = t * t;
        const T t3 = t2 * t;
        const T h00 = (T{2} * t3) - (T{3} * t2) + T{1};
        const T h10 = t3 - (T{2} * t2) + t;
        const T h01 = (T{-2} * t3) + (T{3} * t2);
        const T h11 = t3 - t2;
        return (h00 * yi) + (h10 * h * mi) + (h01 * yi1) + (h11 * h * mi1);
    }
}

/**
 * @brief 2-D interpolating surface with smooth (Catmull-Rom spline) blending.
 *
 * Same shape and query interface as @ref Lut2D — strictly-increasing @ref rows /
 * @ref cols breakpoints and a @ref z grid, `operator()(r, c)` — but interpolates
 * with a separable non-uniform Catmull-Rom spline instead of bilinear, giving a
 * C¹ surface that passes through every grid point. A drop-in surface anywhere a
 * callable `(r, c) -> T` is wanted (e.g. @ref damp::io::steer_map) when the coarse
 * facets of bilinear interpolation are undesirable. Queries clamp to the grid.
 *
 * Local (4-point-per-axis) and allocation-free — no spline pre-solve — so it
 * stays constexpr and on-target friendly. Interior overshoot is possible (as with
 * any Catmull-Rom); clamp downstream if the output must stay bounded.
 */
template<size_t R, size_t C, typename T = float>
struct SplineSurface {
    static_assert(R >= 1 && C >= 1, "SplineSurface needs at least one breakpoint per axis");

    ColVec<R, T>    rows{}; ///< Strictly increasing row breakpoints.
    ColVec<C, T>    cols{}; ///< Strictly increasing column breakpoints.
    Matrix<R, C, T> z{};    ///< Grid values, z(row, col).

    /// Catmull-Rom interpolated value at (@p r, @p c), clamped to the grid edges.
    [[nodiscard]] constexpr T operator()(T r, T c) const {
        return catmull_1d(
            rows, [&](size_t i) { return catmull_1d(cols, [&](size_t j) { return z(i, j); }, c); }, r
        );
    }
};

/**
 * @brief Adaptive 1-D lookup: fixed breakpoints, online-updated cell values.
 *
 * Breakpoints @ref xs stay fixed. Each sample `(x, measured)` updates the
 * nearest cell (@ref Lut1D::nearest_index). Two update laws:
 * - Running mean (`forgetting == false`): @f$ y \leftarrow y + (m-y)/n @f$ with
 *   per-cell sample count @ref n.
 * - Forgetting (`forgetting == true`): @f$ y \leftarrow (1-\alpha)y + \alpha m @f$
 *   with adaptation gain @ref alpha ∈ (0, 1].
 *
 * @ref enabled freezes updates; @ref reset restores @ref ys_init and clears counts;
 * @ref locked skips a cell when true. Out-of-range samples follow @ref adapt_oob.
 *
 * Query uses the same interpolant/extrapolation as @ref Lut1D via @ref table().
 *
 * @tparam N Number of breakpoints.
 * @tparam T Scalar type.
 */
template<size_t N, typename T = float>
struct AdaptiveLut1D {
    static_assert(N >= 1, "AdaptiveLut1D needs at least one breakpoint");

    ColVec<N, T>    xs{};
    ColVec<N, T>    ys{};
    ColVec<N, T>    ys_init{}; ///< Snapshot for @ref reset
    ColVec<N, T>    n{};       ///< Sample counts (running-mean mode)
    ColVec<N, T>    locked{};  ///< nonzero ⇒ cell does not adapt
    Extrapolation   oob{Extrapolation::Clamp};
    Interpolation   interp{Interpolation::Linear};
    AdaptOutOfRange adapt_oob{AdaptOutOfRange::Ignore};
    T               alpha{T{1}}; ///< Forgetting gain ∈ (0, 1]; ignored if !forgetting
    bool            forgetting{false};
    bool            enabled{true};

    /// Build from a static table (initial ys copied for reset).
    [[nodiscard]] static constexpr AdaptiveLut1D from_lut(const Lut1D<N, T>& lut) {
        return AdaptiveLut1D{
            .xs = lut.xs,
            .ys = lut.ys,
            .ys_init = lut.ys,
            .oob = lut.oob,
            .interp = lut.interp,
        };
    }

    [[nodiscard]] constexpr Lut1D<N, T> table() const {
        return Lut1D<N, T>{.xs = xs, .ys = ys, .oob = oob, .interp = interp};
    }

    [[nodiscard]] constexpr T operator()(T x) const { return table()(x); }

    constexpr void set_enabled(bool on) { enabled = on; }
    constexpr void reset() {
        ys = ys_init;
        n = ColVec<N, T>{};
    }

    [[nodiscard]] constexpr bool is_locked(size_t i) const { return locked[i] != T{0}; }
    constexpr void               set_locked(size_t i, bool on) { locked[i] = on ? T{1} : T{0}; }

    /// Incorporate one measurement; no-op if disabled or cell locked / ignored OOR.
    constexpr void observe(T x, T measured) {
        if (!enabled) {
            return;
        }
        if constexpr (N >= 1) {
            if ((x < xs[0] || x > xs[N - 1]) && adapt_oob == AdaptOutOfRange::Ignore) {
                return;
            }
            const size_t i = Lut1D<N, T>{.xs = xs, .ys = ys}.nearest_index(x);
            if (is_locked(i)) {
                return;
            }
            if (forgetting) {
                const T a = damp::clamp(alpha, T{0}, T{1});
                ys[i] = ((T{1} - a) * ys[i]) + (a * measured);
            } else {
                n[i] += T{1};
                ys[i] += (measured - ys[i]) / n[i];
            }
        }
    }
};

/**
 * @brief Adaptive 2-D lookup: fixed grid, online-updated cells (nearest grid node).
 *
 * Same adaptation knobs as @ref AdaptiveLut1D; query is bilinear via @ref table().
 *
 * @tparam R Row breakpoints.
 * @tparam C Column breakpoints.
 * @tparam T Scalar type.
 */
template<size_t R, size_t C, typename T = float>
struct AdaptiveLut2D {
    static_assert(R >= 1 && C >= 1, "AdaptiveLut2D needs at least one breakpoint per axis");

    ColVec<R, T>    rows{};
    ColVec<C, T>    cols{};
    Matrix<R, C, T> z{};
    Matrix<R, C, T> z_init{};
    Matrix<R, C, T> n{}; ///< counts
    // lock packed as 0/1 in a matrix of T to stay POD-friendly
    Matrix<R, C, T> lock{}; ///< nonzero ⇒ locked
    Extrapolation   oob{Extrapolation::Clamp};
    AdaptOutOfRange adapt_oob{AdaptOutOfRange::Ignore};
    T               alpha{T{1}};
    bool            forgetting{false};
    bool            enabled{true};

    [[nodiscard]] static constexpr AdaptiveLut2D from_lut(const Lut2D<R, C, T>& lut) {
        return AdaptiveLut2D{
            .rows = lut.rows,
            .cols = lut.cols,
            .z = lut.z,
            .z_init = lut.z,
            .oob = lut.oob,
        };
    }

    [[nodiscard]] constexpr Lut2D<R, C, T> table() const {
        return Lut2D<R, C, T>{.rows = rows, .cols = cols, .z = z, .oob = oob};
    }

    [[nodiscard]] constexpr T operator()(T r, T c) const { return table()(r, c); }

    constexpr void set_enabled(bool on) { enabled = on; }
    constexpr void reset() {
        z = z_init;
        n = Matrix<R, C, T>{};
    }

    [[nodiscard]] constexpr bool is_locked(size_t i, size_t j) const { return lock(i, j) != T{0}; }
    constexpr void               set_locked(size_t i, size_t j, bool on) { lock(i, j) = on ? T{1} : T{0}; }

    constexpr void observe(T r, T c, T measured) {
        if (!enabled) {
            return;
        }
        if ((r < rows[0] || r > rows[R - 1] || c < cols[0] || c > cols[C - 1])
            && adapt_oob == AdaptOutOfRange::Ignore) {
            return;
        }
        const size_t i = Lut1D<R, T>{.xs = rows, .ys = ColVec<R, T>{}}.nearest_index(r);
        const size_t j = Lut1D<C, T>{.xs = cols, .ys = ColVec<C, T>{}}.nearest_index(c);
        if (is_locked(i, j)) {
            return;
        }
        if (forgetting) {
            const T a = damp::clamp(alpha, T{0}, T{1});
            z(i, j) = ((T{1} - a) * z(i, j)) + (a * measured);
        } else {
            n(i, j) += T{1};
            z(i, j) += (measured - z(i, j)) / n(i, j);
        }
    }
};

} // namespace damp
