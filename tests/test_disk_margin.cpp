// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>

#include "damp/analysis/analysis.hpp"
#include "damp/analysis/disk_margin.hpp"
#include "damp/matlab.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {
// MATLAB®'s diskmargin example: L = 25 / (s³ + 10s² + 10s + 10).
StateSpace<3, 1, 1> make_loop(double k = 25.0) {
    return *TransferFunction<1, 4>{{k}, {1.0, 10.0, 10.0, 10.0}}.to_state_space();
}
} // namespace

TEST_SUITE("Disk margin") {
    TEST_CASE("balanced disk (σ = 0) matches the reference: α, gain interval, phase, frequency") {
        //! Reference: SciPy, dense grid + Brent refinement of max |S − 1/2|; MATLAB® reports
        //! DiskMargin 0.4581, GainMargin [0.6273 1.5942], PhaseMargin ±25.8017, Frequency 1.9550.
        const auto dm = analysis::disk_margin(make_loop());
        REQUIRE(dm.stable);
        REQUIRE(dm.success);
        CHECK(dm.alpha == doctest::Approx(0.4580925477415598).epsilon(1e-9));
        CHECK(dm.frequency == doctest::Approx(1.9550269042822883).epsilon(1e-5));
        CHECK(dm.gain_low == doctest::Approx(0.6272780305506032).epsilon(1e-9));
        CHECK(dm.gain_high == doctest::Approx(1.5941894204747362).epsilon(1e-9));
        CHECK_FALSE(dm.gain_high_unbounded);
        CHECK(dm.phase_margin_deg == doctest::Approx(25.801709378449843).epsilon(1e-8));
        //! σ = 0: φ = 2·atan(α/2), and the gain interval is symmetric in dB.
        CHECK(dm.phase_margin_deg == doctest::Approx(2.0 * std::atan(dm.alpha / 2.0) * 180.0 / std::numbers::pi));
        CHECK(dm.gain_low * dm.gain_high == doctest::Approx(1.0));
    }

    TEST_CASE("skewed disks: σ = 1 is the classical 1/‖S‖∞ disk, σ = −0.5 leans toward T") {
        const auto s1 = analysis::disk_margin(make_loop(), 1.0);
        CHECK(s1.alpha == doctest::Approx(0.40214586422401916).epsilon(1e-9));
        CHECK(s1.gain_low == doctest::Approx(0.7131925611415784).epsilon(1e-9));
        CHECK(s1.gain_low == doctest::Approx(1.0 / (1.0 + s1.alpha)));
        CHECK(s1.gain_high == doctest::Approx(1.6726487953488132).epsilon(1e-9));
        CHECK(s1.phase_margin_deg == doctest::Approx(23.199416378802606).epsilon(1e-8));

        const auto sm = matlab::diskmargin(make_loop(), -0.5);
        CHECK(sm.skew == -0.5);
        CHECK(sm.alpha == doctest::Approx(0.47822106443826895).epsilon(1e-9));
        CHECK(sm.gain_low == doctest::Approx(0.5728472913176694).epsilon(1e-9));
        CHECK(sm.gain_high == doctest::Approx(1.5431585266290901).epsilon(1e-9));
        CHECK(sm.phase_margin_deg == doctest::Approx(27.0822395726329).epsilon(1e-8));
    }

    TEST_CASE("an integrator loop tolerates any gain increase and 90° of phase; an unstable loop has no margin") {
        //! L = 2/s: S − T = (s − 2)/(s + 2), all-pass, so α = 2 exactly.
        const StateSpace<1, 1, 1> integrator{.A = Matrix<1, 1>{{0.0}}, .B = Matrix<1, 1>{{1.0}}, .C = Matrix<1, 1>{{2.0}}};
        const auto                dm = analysis::disk_margin(integrator);
        REQUIRE(dm.success);
        CHECK(dm.alpha == doctest::Approx(2.0));
        CHECK(dm.gain_low == doctest::Approx(0.0).scale(1.0));
        CHECK((dm.gain_high_unbounded || dm.gain_high > 1e6));
        CHECK(dm.phase_margin_deg == doctest::Approx(90.0));

        //! k = 250: s³ + 10s² + 10s + 260 fails Routh (10·10 < 260).
        const auto unstable = analysis::disk_margin(make_loop(250.0));
        CHECK_FALSE(unstable.stable);
        CHECK_FALSE(unstable.success);
        CHECK(unstable.alpha == 0.0);
    }

    TEST_CASE("discrete loop: α = 2 / max|S − T| over the unit circle, and float") {
        const auto loop_d = *discretize(make_loop(), 0.05);
        const auto grid = analysis::linspace(1e-3, 0.999 * std::numbers::pi / 0.05, 4000);
        const auto dm = analysis::disk_margin(loop_d);
        REQUIRE(dm.success);
        CHECK(dm.frequency > 0.0);
        CHECK(dm.frequency < std::numbers::pi / 0.05);
        double peak = 0.0;
        for (const auto& p : analysis::nyquist(loop_d, grid).points) {
            const std::complex<double> L(static_cast<double>(p.real), static_cast<double>(p.imag));
            peak = std::max(peak, std::abs((1.0 - L) / (1.0 + L)));
        }
        CHECK(dm.alpha <= 2.0 / peak); // the exact peak is at least the grid's
        CHECK(dm.alpha == doctest::Approx(2.0 / peak).epsilon(1e-3));
        CHECK(dm.alpha < 0.4580925477415598); // the ZOH sample-and-hold costs margin

        const auto loop_f = make_loop().as<float>();
        const auto dm_f = analysis::disk_margin(loop_f);
        REQUIRE(dm_f.success);
        CHECK(dm_f.alpha == doctest::Approx(0.4580925477415598).epsilon(1e-4));
    }

    TEST_CASE("a lightly damped flexible mode between grid samples sets the margin") {
        //! L = 100k / (s(s + 1)(s² + 20ζs + 100)), k = 3.2ζ√101·10. Reference: SciPy, dense grid
        //! across the mode's half-power band + Brent on |S − 1/2|. A search of logspace(10⁻², 10², 200)
        //! with golden-section refinement finds only the crossover peak: α = 1.3214 (ζ = 10⁻³) and
        //! 1.9873 (ζ = 10⁻⁵).
        struct Case {
            double num;
            double den1;
            double den2;
            double alpha;
            double frequency;
        };
        constexpr Case cases[] = {
            {.num = 32.15960198758685, .den1 = 1.02, .den2 = 100.02, .alpha = 0.4803575799211788, .frequency = 9.982024166498777},
            {.num = 0.32159601987586856, .den1 = 1.0002, .den2 = 100.0002, .alpha = 0.4831772237376791, .frequency = 9.99982104761693},
        };
        for (const auto& c : cases) {
            CAPTURE(c.den1);
            const auto loop = TransferFunction<1, 5>{.num = {c.num}, .den = {1.0, c.den1, c.den2, 100.0, 0.0}}.to_state_space();
            REQUIRE(loop.has_value());
            const auto dm = analysis::disk_margin(*loop);
            REQUIRE(dm.success);
            CHECK(dm.alpha == doctest::Approx(c.alpha).epsilon(1e-10));
            CHECK(dm.frequency == doctest::Approx(c.frequency).epsilon(1e-8));
            CHECK(dm.gain_low * dm.gain_high == doctest::Approx(1.0));
            CHECK(dm.phase_margin_deg == doctest::Approx(2.0 * std::atan(dm.alpha / 2.0) * 180.0 / std::numbers::pi));
        }
    }

    TEST_CASE("loop-at-a-time: a diagonal loop gives each channel's SISO disk margin") {
        const auto                loop1 = make_loop();
        const StateSpace<4, 2, 2> loop{
            .A = blkdiag(loop1.A, Matrix<1, 1>{{0.0}}),
            .B = blkdiag(loop1.B, Matrix<1, 1>{{1.0}}),
            .C = blkdiag(loop1.C, Matrix<1, 1>{{2.0}}),
        };
        const auto lat = analysis::loop_at_a_time_disk_margin(loop);
        REQUIRE(lat.stable);
        REQUIRE(lat.success);
        CHECK(lat.channel[0].alpha == doctest::Approx(0.4580925477415598).epsilon(1e-9));
        CHECK(lat.channel[0].frequency == doctest::Approx(1.9550269042822883).epsilon(1e-5));
        CHECK(lat.channel[1].alpha == doctest::Approx(2.0));
        CHECK(lat.channel[1].phase_margin_deg == doctest::Approx(90.0));
    }

    TEST_CASE("loop-at-a-time: spinning satellite, each loop broken with the other closed") {
        //! L = C(sI − A)⁻¹BK, A = [0 a; −a 0], C = [1 a; −a 1], a = 10. With K = I, T = C/(s + 1),
        //! so S₁₁ = S₂₂ = s/(s + 1) and |S_ii − 1/2| = 1/2: α = 2, any gain increase and ±90° in
        //! either loop alone. With K = diag(1, 2), S₁₁(0) = −50 gives α₁ = 2/101, and
        //! S₂₂ = s/(s + 2) gives α₂ = 2.
        const auto satellite = [](double k2) {
            return StateSpace<2, 2, 2>{
                .A = Matrix<2, 2>{{0.0, 10.0}, {-10.0, 0.0}},
                .B = Matrix<2, 2>{{1.0, 0.0}, {0.0, k2}},
                .C = Matrix<2, 2>{{1.0, 10.0}, {-10.0, 1.0}},
            };
        };
        const auto lat = analysis::loop_at_a_time_disk_margin(satellite(1.0));
        REQUIRE(lat.success);
        for (const auto& ch : lat.channel) {
            CHECK(ch.alpha == doctest::Approx(2.0).epsilon(1e-12));
            CHECK((ch.gain_high_unbounded || ch.gain_high > 1e6));
            CHECK(ch.phase_margin_deg == doctest::Approx(90.0));
        }

        const auto skewed = matlab::diskmargin(satellite(2.0), 0.0);
        REQUIRE(skewed.success);
        CHECK(skewed.channel[0].alpha == doctest::Approx(2.0 / 101.0).epsilon(1e-12));
        CHECK(skewed.channel[0].frequency == 0.0);
        CHECK(skewed.channel[1].alpha == doctest::Approx(2.0).epsilon(1e-12));
    }

    TEST_CASE("loop-at-a-time: an unstable MIMO closed loop has no margins") {
        //! Positive feedback on the satellite: A + C has eigenvalues 1 ± 20j.
        const StateSpace<2, 2, 2> loop{
            .A = Matrix<2, 2>{{0.0, 10.0}, {-10.0, 0.0}},
            .B = Matrix<2, 2>{{-1.0, 0.0}, {0.0, -1.0}},
            .C = Matrix<2, 2>{{1.0, 10.0}, {-10.0, 1.0}},
        };
        const auto lat = analysis::loop_at_a_time_disk_margin(loop);
        CHECK_FALSE(lat.stable);
        CHECK_FALSE(lat.success);
        for (const auto& ch : lat.channel) {
            CHECK(ch.alpha == 0.0);
            CHECK_FALSE(ch.success);
        }
    }
}
