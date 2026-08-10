// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>

#include "damp/matrix/matrix.hpp"
#include "damp/toolbox/lookup.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

// Breakpoint lookup tables: 1-D linear/nearest + 2-D bilinear + 3-D trilinear.

TEST_SUITE("Lookup tables") {
    TEST_CASE("Lut1D hits breakpoints and interpolates between") {
        constexpr Lut1D<3, double> lut{{0.0, 10.0, 20.0}, {0.0, 100.0, 0.0}};
        // Exact breakpoints.
        CHECK(lut(0.0) == doctest::Approx(0.0));
        CHECK(lut(10.0) == doctest::Approx(100.0));
        CHECK(lut(20.0) == doctest::Approx(0.0));
        // Linear interpolation within segments.
        CHECK(lut(5.0) == doctest::Approx(50.0));
        CHECK(lut(15.0) == doctest::Approx(50.0));
    }

    TEST_CASE("Lut1D out-of-range: clamp vs linear extrapolation") {
        constexpr Lut1D<2, double> clamp_lut{{0.0, 10.0}, {0.0, 100.0}, Extrapolation::Clamp};
        CHECK(clamp_lut(-5.0) == doctest::Approx(0.0));   // held at first value
        CHECK(clamp_lut(15.0) == doctest::Approx(100.0)); // held at last value

        constexpr Lut1D<2, double> lin_lut{{0.0, 10.0}, {0.0, 100.0}, Extrapolation::Linear};
        CHECK(lin_lut(-5.0) == doctest::Approx(-50.0)); // slope continues
        CHECK(lin_lut(15.0) == doctest::Approx(150.0));
    }

    TEST_CASE("Lut1D nearest-neighbour") {
        constexpr Lut1D<3, double> lut{{0.0, 10.0, 20.0}, {1.0, 2.0, 3.0}};
        CHECK(lut.nearest(2.0) == doctest::Approx(1.0));  // closer to 0
        CHECK(lut.nearest(8.0) == doctest::Approx(2.0));  // closer to 10
        CHECK(lut.nearest(16.0) == doctest::Approx(3.0)); // closer to 20
    }

    TEST_CASE("Lut2D bilinear over a unit grid") {
        constexpr Lut2D<2, 2, double> lut{
            {0.0, 1.0}, // rows
            {0.0, 1.0}, // cols
            {{0.0, 10.0}, {20.0, 30.0}}
        };
        // Corners are the grid values.
        CHECK(lut(0.0, 0.0) == doctest::Approx(0.0));
        CHECK(lut(0.0, 1.0) == doctest::Approx(10.0));
        CHECK(lut(1.0, 0.0) == doctest::Approx(20.0));
        CHECK(lut(1.0, 1.0) == doctest::Approx(30.0));
        // Center is the average of all four.
        CHECK(lut(0.5, 0.5) == doctest::Approx(15.0));
        // Edge midpoints.
        CHECK(lut(0.0, 0.5) == doctest::Approx(5.0));
        CHECK(lut(0.5, 0.0) == doctest::Approx(10.0));
        // Out of range clamps to the grid edge.
        CHECK(lut(-1.0, -1.0) == doctest::Approx(0.0));
        CHECK(lut(2.0, 2.0) == doctest::Approx(30.0));
    }

    TEST_CASE("SplineSurface passes through nodes and stays smooth off them") {
        // Non-uniform breakpoints, a curved grid so the spline differs from bilinear.
        constexpr SplineSurface<3, 3, double> s{
            {0.0, 0.3, 1.0}, // rows (non-uniform)
            {0.0, 0.5, 1.0}, // cols
            {{0.0, 1.0, 0.0}, {1.0, 2.0, 1.0}, {0.0, 1.0, 0.0}}
        };

        // Interpolating property: exact at every grid node.
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                CHECK(s(s.rows[i], s.cols[j]) == doctest::Approx(s.z(i, j)));
            }
        }

        // Out of range clamps to the corner node.
        CHECK(s(-1.0, -1.0) == doctest::Approx(0.0));
        CHECK(s(2.0, 2.0) == doctest::Approx(0.0));

        // Between nodes the spline overshoots the bilinear value toward the
        // curvature (center row bulges up): sampled point exceeds plain average.
        const double bilinear_center = (0.0 + 1.0 + 1.0 + 2.0) / 4.0; // rows 0..1, cols 0..1 corner avg
        CHECK(s(0.15, 0.25) > bilinear_center - 1.0);                 // sane, finite, near the hump
        CHECK(s(0.3, 0.5) == doctest::Approx(2.0));                   // the peak node itself
    }

    TEST_CASE("catmull_1d reduces to linear for two breakpoints") {
        constexpr ColVec<2, double> xs{0.0, 2.0};
        auto                        y = [](size_t i) { return i == 0 ? 1.0 : 5.0; };
        CHECK(catmull_1d(xs, y, 0.0) == doctest::Approx(1.0));
        CHECK(catmull_1d(xs, y, 2.0) == doctest::Approx(5.0));
        CHECK(catmull_1d(xs, y, 1.0) == doctest::Approx(3.0)); // midpoint = linear
    }

    TEST_CASE("Lut1D cubic mode hits nodes") {
        constexpr Lut1D<4, double> lut{
            {0.0, 1.0, 2.0, 3.0},
            {0.0, 1.0, 0.0, 1.0},
            Extrapolation::Clamp,
            Interpolation::Cubic
        };
        CHECK(lut(0.0) == doctest::Approx(0.0));
        CHECK(lut(1.0) == doctest::Approx(1.0));
        CHECK(lut(2.0) == doctest::Approx(0.0));
        CHECK(lut(3.0) == doctest::Approx(1.0));
        CHECK(lut(0.5) > 0.0);
        CHECK(lut(0.5) < 1.0);
    }

    TEST_CASE("Lut2D linear extrapolation continues edge slope") {
        constexpr Lut2D<2, 2, double> clamp_lut{
            {0.0, 1.0},
            {0.0, 1.0},
            {{0.0, 10.0}, {20.0, 30.0}},
            Extrapolation::Clamp
        };
        CHECK(clamp_lut(-1.0, 0.0) == doctest::Approx(0.0));

        constexpr Lut2D<2, 2, double> lin_lut{
            {0.0, 1.0},
            {0.0, 1.0},
            {{0.0, 10.0}, {20.0, 30.0}},
            Extrapolation::Linear
        };
        // Along row=0, col past 1: top edge 0→10, slope 10 → at col=2 expect 20
        CHECK(lin_lut(0.0, 2.0) == doctest::Approx(20.0));
    }

    TEST_CASE("Lut3D constant table and exact grid nodes") {
        // v(i,j,k) = 7 everywhere.
        constexpr Lut3D<2, 2, 2, double> const_lut{
            {0.0, 1.0},
            {0.0, 1.0},
            {0.0, 1.0},
            {7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0, 7.0}
        };
        CHECK(const_lut(0.0, 0.0, 0.0) == doctest::Approx(7.0));
        CHECK(const_lut(0.3, 0.7, 0.5) == doctest::Approx(7.0));
        CHECK(const_lut(1.0, 1.0, 1.0) == doctest::Approx(7.0));

        // Separable plane f = x + y + z on the unit cube (values = i+j+k at nodes).
        // Flat layout: index = i*(2*2) + j*2 + k
        constexpr Lut3D<2, 2, 2, double> plane{
            {0.0, 1.0},
            {0.0, 1.0},
            {0.0, 1.0},
            {// i=0,j=0,k=0..1
             0.0, 1.0,
             // i=0,j=1,k=0..1
             1.0, 2.0,
             // i=1,j=0,k=0..1
             1.0, 2.0,
             // i=1,j=1,k=0..1
             2.0, 3.0
            }
        };
        // Exact at all eight nodes.
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                for (size_t k = 0; k < 2; ++k) {
                    const double x = static_cast<double>(i);
                    const double y = static_cast<double>(j);
                    const double z = static_cast<double>(k);
                    CHECK(plane(x, y, z) == doctest::Approx(x + y + z));
                    CHECK(plane.at(i, j, k) == doctest::Approx(x + y + z));
                }
            }
        }
        // Interior: trilinear of a linear field is exact.
        CHECK(plane(0.5, 0.5, 0.5) == doctest::Approx(1.5));
        CHECK(plane(0.25, 0.0, 0.75) == doctest::Approx(1.0));
    }

    TEST_CASE("Lut3D out-of-range: clamp vs linear extrapolation") {
        // f = x on [0,1]^3 (independent of y,z): nodes store the x breakpoint value.
        constexpr Lut3D<2, 2, 2, double> clamp_lut{
            {0.0, 1.0},
            {0.0, 1.0},
            {0.0, 1.0},
            {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0},
            Extrapolation::Clamp
        };
        CHECK(clamp_lut(-1.0, 0.5, 0.5) == doctest::Approx(0.0));
        CHECK(clamp_lut(2.0, 0.5, 0.5) == doctest::Approx(1.0));
        CHECK(clamp_lut(0.5, -1.0, 2.0) == doctest::Approx(0.5));

        constexpr Lut3D<2, 2, 2, double> lin_lut{
            {0.0, 1.0},
            {0.0, 1.0},
            {0.0, 1.0},
            {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0},
            Extrapolation::Linear
        };
        // Along x past 1: slope 1 → at x=2 expect 2.
        CHECK(lin_lut(2.0, 0.0, 0.0) == doctest::Approx(2.0));
        CHECK(lin_lut(-1.0, 0.0, 0.0) == doctest::Approx(-1.0));
    }

    TEST_CASE("Lut3D float constexpr construction smoke") {
        constexpr Lut3D<2, 2, 2, float> lut{
            {0.0f, 1.0f},
            {0.0f, 1.0f},
            {0.0f, 1.0f},
            {0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f}
        };
        static_assert(lut.at(0, 0, 0) == 0.0f);
        static_assert(lut.at(1, 1, 1) == 7.0f);
        constexpr float mid = lut(0.5f, 0.5f, 0.5f);
        CHECK(mid == doctest::Approx(3.5f));
    }

    TEST_CASE("Lut3D hydraulic-shaped axes (Vg, omega, p) trilinear mid-cell") {
        // η_v = f(V_g, ω [rad/s], p [bar]): flat layout k (pressure) fastest.
        // Nodes store η_v; mid-cell average of the eight corners.
        constexpr Lut3D<2, 2, 2, double> eta{
            {0.0, 1.0},   // V_g
            {0.0, 200.0}, // ω [rad/s]
            {0.0, 300.0}, // p [bar]
            {0.90, 0.88, 0.92, 0.90, 0.93, 0.91, 0.95, 0.93},
            Extrapolation::Clamp
        };
        // Exact at low-displacement, low-speed, low-pressure corner.
        CHECK(eta(0.0, 0.0, 0.0) == doctest::Approx(0.90));
        // Mid-cell (0.5, 100, 150): average of all eight nodes.
        const double mid_avg = (0.90 + 0.88 + 0.92 + 0.90 + 0.93 + 0.91 + 0.95 + 0.93) / 8.0;
        CHECK(eta(0.5, 100.0, 150.0) == doctest::Approx(mid_avg));
        // OOB clamp holds the high-pressure face at V_g=1, ω=200.
        CHECK(eta(1.0, 200.0, 500.0) == doctest::Approx(0.93));
    }

    TEST_CASE("Lut1D single breakpoint and cubic linear-oob stay finite") {
        constexpr Lut1D<1, double> one{{3.0}, {42.0}};
        CHECK(one(0.0) == doctest::Approx(42.0));
        CHECK(one.nearest(99.0) == doctest::Approx(42.0));

        constexpr Lut1D<3, double> cubic{
            {0.0, 1.0, 2.0},
            {0.0, 1.0, 0.0},
            Extrapolation::Linear,
            Interpolation::Cubic
        };
        // Linear oob continues the first/last segment slope of the cubic table.
        CHECK(cubic(-1.0) == doctest::Approx(-1.0)); // first segment 0→1 slope 1
        CHECK(cubic(3.0) == doctest::Approx(-1.0));  // last segment 1→0 slope -1
        CHECK(std::isfinite(cubic(1.5)));
    }

    TEST_CASE("AdaptiveLut1D running mean converges and lock freezes") {
        Lut1D<3, double> base{{0.0, 1.0, 2.0}, {0.0, 0.0, 0.0}};
        auto             a = AdaptiveLut1D<3, double>::from_lut(base);
        for (int k = 0; k < 20; ++k) {
            a.observe(1.0, 10.0);
        }
        CHECK(a(1.0) == doctest::Approx(10.0).epsilon(1e-9));

        a.set_locked(1, true);
        a.observe(1.0, 0.0);
        CHECK(a(1.0) == doctest::Approx(10.0).epsilon(1e-9)); // locked

        a.set_locked(1, false);
        a.set_enabled(false);
        a.observe(1.0, 0.0);
        CHECK(a(1.0) == doctest::Approx(10.0).epsilon(1e-9)); // disabled

        a.reset();
        CHECK(a(1.0) == doctest::Approx(0.0));
    }

    TEST_CASE("AdaptiveLut1D forgetting pulls toward new samples") {
        auto a = AdaptiveLut1D<2, double>::from_lut(Lut1D<2, double>{{0.0, 1.0}, {0.0, 0.0}});
        a.forgetting = true;
        a.alpha = 0.5;
        a.observe(0.0, 8.0);
        CHECK(a.ys[0] == doctest::Approx(4.0)); // 0.5*0 + 0.5*8
        a.observe(0.0, 8.0);
        CHECK(a.ys[0] == doctest::Approx(6.0));
    }
}
