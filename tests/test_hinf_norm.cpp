// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>
#include <limits>

#include "damp/analysis/norms.hpp"
#include "damp/matlab.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

//! References not in closed form: a 400 001-point log or linear grid of σ̄ plus a 4001-point
//! grid across each pole's half-power band, then bounded Brent (SciPy) on a unit-width variable
//! around each local maximum, three times with the bracket shrunk 10³. The grid reads low on
//! these peaks by up to 2·10⁻⁴; the refined values agree with this file's results to the
//! realization's conditioning (≈ 10⁻¹² at ζ = 10⁻⁵).

namespace {
constexpr double kZeta[] = {1e-1, 1e-2, 1e-3, 1e-4, 1e-5};

// ωₙ² / (s² + 2ζωₙ s + ωₙ²) in controllable canonical form, ωₙ = 10.
StateSpace<2, 1, 1> resonance(double zeta, double d = 0.0) {
    return StateSpace<2, 1, 1>{
        .A = Matrix<2, 2>{{0.0, 1.0}, {-100.0, -20.0 * zeta}},
        .B = Matrix<2, 1>{{0.0}, {1.0}},
        .C = Matrix<1, 2>{{100.0, 0.0}},
        .D = Matrix<1, 1>{{d}},
    };
}

constexpr StateSpace<2, 1, 1> kResonance{
    .A = Matrix<2, 2>{{0.0, 1.0}, {-100.0, -0.2}},
    .B = Matrix<2, 1>{{0.0}, {1.0}},
    .C = Matrix<1, 2>{{100.0, 0.0}},
};
constexpr auto kResonancePeak = analysis::norm_hinf(kResonance);
static_assert(kResonancePeak.success);
static_assert(kResonancePeak.gain > 50.0025 && kResonancePeak.gain < 50.0026); // 1/(2ζ√(1−ζ²)), ζ = 0.01
} // namespace

TEST_SUITE("H-infinity norm") {
    TEST_CASE("second order: peak 1/(2ζ√(1−ζ²)) at ωₙ√(1−2ζ²) for ζ from 10⁻¹ to 10⁻⁵, canonical and modal") {
        for (const double zeta : kZeta) {
            CAPTURE(zeta);
            const double peak = 1.0 / (2.0 * zeta * std::sqrt(1.0 - (zeta * zeta)));
            const double w_peak = 10.0 * std::sqrt(1.0 - (2.0 * zeta * zeta));
            //! The peak is flat to (Δω/ζωₙ)², so the frequency is good to about √tol of the half-power width.
            const double w_tol = 1e-5 * zeta * 10.0;

            const auto canonical = analysis::norm_hinf(resonance(zeta));
            REQUIRE(canonical.success);
            CHECK(canonical.gain == doctest::Approx(peak).epsilon(1e-12));
            CHECK(std::abs(canonical.frequency - w_peak) <= w_tol);

            //! Real modal form of the same mode: −ζωₙ ± jω_d.
            const double              wd = 10.0 * std::sqrt(1.0 - (zeta * zeta));
            const StateSpace<2, 1, 1> modal{
                .A = Matrix<2, 2>{{-10.0 * zeta, wd}, {-wd, -10.0 * zeta}},
                .B = Matrix<2, 1>{{0.0}, {1.0}},
                .C = Matrix<1, 2>{{100.0 / wd, 0.0}},
            };
            const auto m = analysis::norm_hinf(modal);
            REQUIRE(m.success);
            CHECK(m.gain == doctest::Approx(peak).epsilon(1e-12));
            CHECK(std::abs(m.frequency - w_peak) <= w_tol);
        }
    }

    TEST_CASE("fourth order: (50s + 100) / ((s² + 0.6s + 1)(s² + 20ζs + 100)) from its transfer function") {
        //! Denominators expanded in double (numpy.polymul) and the SciPy reference peaks.
        constexpr double den[5][5] = {
            {1.0, 2.6, 102.2, 62.0, 100.0},
            {1.0, 0.8, 101.12, 60.2, 100.0},
            {1.0, 0.62, 101.012, 60.02, 100.0},
            {1.0, 0.602, 101.0012, 60.002, 100.0},
            {1.0, 0.6002, 101.00012, 60.0002, 100.0},
        };
        constexpr double gain[5] = {1.9374870304026381, 2.5710886936298367, 25.7055120626421, 257.05458331507543, 2570.5457794220392};
        constexpr double freq[5] = {0.9258054186853084, 9.997944391616468, 9.9999794506173, 9.999999794503573, 9.999999997944897};
        for (size_t k = 0; k < 5; ++k) {
            CAPTURE(kZeta[k]);
            const auto sys = TransferFunction<2, 5>{.num = {50.0, 100.0}, .den = {den[k][0], den[k][1], den[k][2], den[k][3], den[k][4]}}.to_state_space();
            REQUIRE(sys.has_value());
            const auto r = analysis::norm_hinf(*sys);
            REQUIRE(r.success);
            CHECK(r.gain == doctest::Approx(gain[k]).epsilon(1e-10));
            CHECK(r.frequency == doctest::Approx(freq[k]).epsilon(1e-5 * kZeta[k]));
        }
    }

    TEST_CASE("discrete resonance near Nyquist: poles r·e^{±j0.95π}, r = e^{−ζ·0.95π/√(1−ζ²)}") {
        //! A = r·[cos θ, −sin θ; sin θ, cos θ] (entries printed from Python), B = e₂, C = e₁, Ts = 1 ms.
        constexpr double a[5][2] = {
            {-0.7317305598970453, 0.1158947351976532},
            {-0.9586447618558084, 0.15183441407658071},
            {-0.9847449648667347, 0.15596828012282307},
            {-0.9873936077062482, 0.15638778393657898},
            {-0.9876588633478834, 0.1564297963029223},
        };
        constexpr double gain[5] = {1.3570799031327254, 16.749827265794362, 167.53118656903132, 1675.315157193318, 16753.151901083627};
        constexpr double freq[5] = {3141.592653589793, 2987.3508678793755, 2984.5411427638132, 2984.5133021455404, 2984.513023696363};
        for (size_t k = 0; k < 5; ++k) {
            CAPTURE(kZeta[k]);
            const StateSpace<2, 1, 1> sys{
                .A = Matrix<2, 2>{{a[k][0], -a[k][1]}, {a[k][1], a[k][0]}},
                .B = Matrix<2, 1>{{0.0}, {1.0}},
                .C = Matrix<1, 2>{{1.0, 0.0}},
                .Ts = 1e-3,
            };
            const auto r = analysis::norm_hinf(sys);
            REQUIRE(r.success);
            CHECK(r.gain == doctest::Approx(gain[k]).epsilon(1e-10));
            CHECK(r.frequency == doctest::Approx(freq[k]).epsilon(1e-5 * kZeta[k]));
        }
        //! ζ = 0.1 peaks at Nyquist: θ = π exactly, π/Ts.
        CHECK(analysis::norm_hinf(StateSpace<2, 1, 1>{
                                      .A = Matrix<2, 2>{{a[0][0], -a[0][1]}, {a[0][1], a[0][0]}},
                                      .B = Matrix<2, 1>{{0.0}, {1.0}},
                                      .C = Matrix<1, 2>{{1.0, 0.0}},
                                      .Ts = 1e-3,
                                  })
                  .frequency
              == doctest::Approx(3141.592653589793).epsilon(1e-12));
    }

    TEST_CASE("MIMO 3×2 with D ≠ 0: two lightly damped modes (ζ = 10⁻³ at 3 rad/s, 10⁻⁴ at 8 rad/s)") {
        const StateSpace<4, 2, 3> sys{
            .A = Matrix<4, 4>{
                {-0.003, 2.999998499999625, 0.0, 0.0},
                {-2.999998499999625, -0.003, 0.0, 0.0},
                {0.0, 0.0, -0.0008, 7.999999959999999},
                {0.0, 0.0, -7.999999959999999, -0.0008},
            },
            .B = Matrix<4, 2>{{1.0, 0.0}, {0.5, 1.0}, {0.0, 1.0}, {1.0, -0.5}},
            .C = Matrix<3, 4>{{1.0, 0.0, 0.2, 0.0}, {0.0, 1.0, 1.0, 0.0}, {0.3, 0.0, 0.0, 1.0}},
            .D = Matrix<3, 2>{{0.1, 0.0}, {0.0, -0.2}, {0.05, 0.05}},
        };
        const auto r = analysis::norm_hinf(sys);
        REQUIRE(r.success);
        CHECK(r.gain == doctest::Approx(1338.944093012454).epsilon(1e-10));
        CHECK(r.frequency == doctest::Approx(7.999999973350349).epsilon(1e-9));
    }

    TEST_CASE("D ≠ 0: a finite peak moved by the feedthrough, and a peak only approached as ω → ∞") {
        //! 0.5 + resonance (ζ = 10⁻³): the feedthrough adds in quadrature at ωₙ and shifts the peak.
        const auto r = analysis::norm_hinf(resonance(1e-3, 0.5));
        REQUIRE(r.success);
        CHECK(r.gain == doctest::Approx(500.00124999843763).epsilon(1e-10));
        CHECK(r.gain > 1.0 / (2e-3 * std::sqrt(1.0 - 1e-6))); // above the bare resonance
        CHECK(r.frequency == doctest::Approx(9.999980000068012).epsilon(1e-8));

        //! (s + 0.1)/(s + 1): |G(jω)| < 1 at every finite ω, → 1 = |D| as ω → ∞.
        const StateSpace<1, 1, 1> high_pass{.A = Matrix<1, 1>{{-1.0}}, .B = Matrix<1, 1>{{1.0}}, .C = Matrix<1, 1>{{-0.9}}, .D = Matrix<1, 1>{{1.0}}};
        const auto                hp = analysis::norm_hinf(high_pass);
        REQUIRE(hp.success);
        CHECK(hp.gain == 1.0);
        CHECK(hp.frequency == std::numeric_limits<double>::max());
    }

    TEST_CASE("an unstable or marginally stable G is refused, as norm_h2 refuses it") {
        const StateSpace<1, 1, 1> unstable{.A = Matrix<1, 1>{{0.5}}, .B = Matrix<1, 1>{{1.0}}, .C = Matrix<1, 1>{{1.0}}};
        const StateSpace<1, 1, 1> integrator{.A = Matrix<1, 1>{{0.0}}, .B = Matrix<1, 1>{{1.0}}, .C = Matrix<1, 1>{{1.0}}};
        const StateSpace<1, 1, 1> discrete{.A = Matrix<1, 1>{{1.02}}, .B = Matrix<1, 1>{{1.0}}, .C = Matrix<1, 1>{{1.0}}, .Ts = 0.1};
        const StateSpace<1, 1, 1> on_circle{.A = Matrix<1, 1>{{-1.0}}, .B = Matrix<1, 1>{{1.0}}, .C = Matrix<1, 1>{{1.0}}, .Ts = 0.1};
        CHECK_FALSE(analysis::norm_hinf(unstable).success);
        CHECK_FALSE(analysis::norm_hinf(integrator).success);
        CHECK_FALSE(analysis::norm_hinf(discrete).success);
        CHECK_FALSE(analysis::norm_hinf(on_circle).success);
        CHECK_FALSE(analysis::norm_h2(unstable).has_value());
        CHECK_FALSE(analysis::norm_h2(discrete).has_value());
        CHECK_FALSE(analysis::norm_hinf(resonance(1e-2), 0.0).success); // a tolerance must be positive
    }

    TEST_CASE("a transfer function that is identically zero has norm 0") {
        //! The driven mode is not observed and the observed mode is not driven.
        const StateSpace<2, 1, 1> hidden{
            .A = Matrix<2, 2>{{-1.0, 0.0}, {0.0, -2.0}},
            .B = Matrix<2, 1>{{1.0}, {0.0}},
            .C = Matrix<1, 2>{{0.0, 1.0}},
        };
        const auto r = analysis::norm_hinf(hidden);
        CHECK(r.success);
        CHECK(r.gain == 0.0);
    }

    TEST_CASE("float, and the MATLAB® alias with a tolerance") {
        for (const double zeta : {1e-1, 1e-2, 1e-3}) {
            CAPTURE(zeta);
            const auto r = analysis::norm_hinf(resonance(zeta).as<float>());
            REQUIRE(r.success);
            CHECK(static_cast<double>(r.gain) == doctest::Approx(1.0 / (2.0 * zeta * std::sqrt(1.0 - (zeta * zeta)))).epsilon(2e-6));
        }
        const auto loose = matlab::hinfnorm(resonance(1e-2), 1e-3);
        REQUIRE(loose.success);
        const double peak = 1.0 / (2e-2 * std::sqrt(1.0 - 1e-4));
        CHECK(loose.gain <= peak);
        CHECK(loose.gain >= peak / (1.0 + 2e-3));
    }
}
