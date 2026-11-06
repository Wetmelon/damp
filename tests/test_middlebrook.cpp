// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file test_middlebrook.cpp
 * @brief Tests for Middlebrook impedance stability analysis
 */
#include <vector>

#include "damp/analysis/analysis.hpp"
#include "damp/backend.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using namespace damp::analysis;

// ============================================================================
// Simple RC source impedance: Z_s(s) = R + 1/(sC)
//
// As an admittance system: Y_s(s) = sC / (sCR + 1)
//   State-space (admittance): x = i_C
//     dx/dt = -1/(RC) * x + 1/(RC) * v
//     i_out = x
//   Actually, use transfer function form:
//     Y_s = sC / (sRC + 1) = (C/RC) / (s + 1/RC) = (1/R) * s/(s + 1/RC)
//   But that has a zero at origin, let's use an impedance model directly.
//
// For testing, use simple 1st-order systems where we know the answer.
// ============================================================================

TEST_CASE("Middlebrook: impedance - pure resistance") {
    // Z(s) = R (constant), represented as D-only system with no states
    // StateSpace: A=[], B=[], C=[], D=[R]
    // We can't have NX=0, so use a trivially decoupled 1st-order system:
    //   dx/dt = -1 * x + 0 * u
    //   y = 0*x + R*u
    // This gives G(s) = 0*1/(s+1) + R = R for all s.

    constexpr double R = 10.0;
    constexpr auto   A = Matrix<1, 1, double>{{-1.0}};
    constexpr auto   B = Matrix<1, 1, double>{{0.0}};
    constexpr auto   C = Matrix<1, 1, double>{{0.0}};
    constexpr auto   D = Matrix<1, 1, double>{{R}};
    constexpr auto   sys = StateSpace<1, 1, 1>{A, B, C, D};

    auto omega = geomspace(1.0, 10000.0, 50);
    auto result = impedance(sys, omega);

    REQUIRE(result.points.size() == 50);

    for (const auto& pt : result.points) {
        CHECK(pt.magnitude == doctest::Approx(R).epsilon(1e-10));
        CHECK(pt.magnitude_db == doctest::Approx(20.0).epsilon(1e-6));
        CHECK(pt.phase_deg == doctest::Approx(0.0).epsilon(1e-6));
    }
}

TEST_CASE("Middlebrook: impedance - RL impedance") {
    // Z(s) = R + sL = L(s + R/L)
    // As a state-space impedance model (V = Z * I):
    //   Trick: we can represent Z(s) = R + sL as a non-proper TF,
    //   but SS requires proper TFs. Instead, represent as:
    //   Z(s) = R + sL via SS with derivative. We'll just check at DC (s=0).
    //
    // Actually, let's just use a first-order Z(s) = R/(1 + s*tau) (RL divider)
    // Z(s) = R / (1 + sL/R) = R * (R/L) / (s + R/L)
    //
    // SS form:
    //   dx/dt = -(R/L)*x + (R/L)*u  => scalar a = -R/L, b = R/L
    //   y     = R*x                   => actually let's derive properly
    //
    // Actually Z(s) = R/(1+sL/R) = (R^2/L)/(s + R/L)
    // So: A = [-R/L], B = [1], C = [R^2/L], D = [0]
    // Check: G(s) = C(sI-A)^{-1}B + D = (R^2/L)/(s+R/L) * 1 = R/(1+sL/R) ✓

    constexpr double R = 5.0;
    constexpr double L = 0.001;    // 1 mH
    constexpr double pole = R / L; // 5000

    constexpr auto A = Matrix<1, 1, double>{{-pole}};
    constexpr auto B = Matrix<1, 1, double>{{1.0}};
    constexpr auto C = Matrix<1, 1, double>{{R * pole}};
    constexpr auto D = Matrix<1, 1, double>{{0.0}};
    constexpr auto sys = StateSpace<1, 1, 1>{A, B, C, D};

    auto omega = geomspace(1.0, 1000000.0, 200);
    auto result = impedance(sys, omega);

    // At DC (ω→0): Z ≈ R = 5 Ω
    CHECK(result.points[0].magnitude == doctest::Approx(R).epsilon(1e-3));

    // At high frequency (ω >> pole): Z → 0 (short circuit)
    CHECK(result.points.back().magnitude < 0.1);

    // At ω = pole: |Z| = R/√2
    // Find the point closest to ω = pole
    double r_sqrt2 = R / damp::sqrt(2.0);
    for (const auto& pt : result.points) {
        if (damp::abs(pt.omega - pole) / pole < 0.1) {
            CHECK(pt.magnitude == doctest::Approx(r_sqrt2).epsilon(0.1));
            break;
        }
    }
}

TEST_CASE("Middlebrook: stable system - Zs << ZL") {
    // Source: low impedance, Z_s(s) = 0.1 / (s/1000 + 1)
    //   DC value = 0.1 Ω, pole at 1000 rad/s
    //   A_s = [-1000], B_s = [1], C_s = [100], D_s = [0]
    //   G(s) = 100/(s+1000) = 0.1/(1+s/1000) ✓

    constexpr auto A_s = Matrix<1, 1, double>{{-1000.0}};
    constexpr auto B_s = Matrix<1, 1, double>{{1.0}};
    constexpr auto C_s = Matrix<1, 1, double>{{100.0}};
    constexpr auto D_s = Matrix<1, 1, double>{{0.0}};
    constexpr auto Zs_sys = StateSpace<1, 1, 1>{A_s, B_s, C_s, D_s};

    // Load: high impedance, Z_L(s) = 100 / (s/500 + 1)
    //   DC value = 100 Ω, pole at 500 rad/s
    //   A_L = [-500], B_L = [1], C_L = [50000], D_L = [0]
    //   G(s) = 50000/(s+500) = 100/(1+s/500) ✓

    constexpr auto A_L = Matrix<1, 1, double>{{-500.0}};
    constexpr auto B_L = Matrix<1, 1, double>{{1.0}};
    constexpr auto C_L = Matrix<1, 1, double>{{50000.0}};
    constexpr auto D_L = Matrix<1, 1, double>{{0.0}};
    constexpr auto ZL_sys = StateSpace<1, 1, 1>{A_L, B_L, C_L, D_L};

    auto omega = geomspace(0.1, 100000.0, 200);

    // Middlebrook from impedance data
    auto Zs_data = impedance(Zs_sys, omega);
    auto ZL_data = impedance(ZL_sys, omega);
    auto result = middlebrook(Zs_data, ZL_data);

    // Sufficient condition should hold: |Zs| < |ZL| everywhere
    CHECK(result.is_stable_sufficient());

    // At DC: Tm = 0.1/100 = 0.001 => -60 dB
    CHECK(result.minor_loop_gain.points[0].magnitude_db == doctest::Approx(-60.0).epsilon(0.5));

    // Worst case margin should be > 1 (since Zs << ZL)
    auto [ratio, freq] = result.worst_case_margin();
    CHECK(ratio > 1.0);
}

TEST_CASE("Middlebrook: unstable system - Zs > ZL at some frequencies") {
    // Source: high impedance at mid frequencies
    // Z_s = 100 / (s/1000 + 1)  => DC = 100 Ω, drops above 1 krad/s
    constexpr auto A_s = Matrix<1, 1, double>{{-1000.0}};
    constexpr auto B_s = Matrix<1, 1, double>{{1.0}};
    constexpr auto C_s = Matrix<1, 1, double>{{100000.0}};
    constexpr auto D_s = Matrix<1, 1, double>{{0.0}};
    constexpr auto Zs_sys = StateSpace<1, 1, 1>{A_s, B_s, C_s, D_s};

    // Load: low impedance at low frequencies
    // Z_L = 1 / (s/10000 + 1)  => DC = 1 Ω
    constexpr auto A_L = Matrix<1, 1, double>{{-10000.0}};
    constexpr auto B_L = Matrix<1, 1, double>{{1.0}};
    constexpr auto C_L = Matrix<1, 1, double>{{10000.0}};
    constexpr auto D_L = Matrix<1, 1, double>{{0.0}};
    constexpr auto ZL_sys = StateSpace<1, 1, 1>{A_L, B_L, C_L, D_L};

    auto omega = geomspace(0.1, 100000.0, 200);

    auto Zs_data = impedance(Zs_sys, omega);
    auto ZL_data = impedance(ZL_sys, omega);
    auto result = middlebrook(Zs_data, ZL_data);

    // Sufficient condition should NOT hold: Zs > ZL at low frequencies
    CHECK_FALSE(result.is_stable_sufficient());

    // At DC: Tm = 100/1 = 100 => 40 dB
    CHECK(result.minor_loop_gain.points[0].magnitude_db == doctest::Approx(40.0).epsilon(0.5));

    // Worst case margin should be < 1
    auto [ratio, freq] = result.worst_case_margin();
    CHECK(ratio < 1.0);
}

TEST_CASE("Middlebrook: admittance models through impedance_from_admittance") {
    // Source admittance: Y_s(s) = (1/R_s) * s/(s + p_s) with R_s = 0.1, p_s = 1000
    // => Z_s = R_s * (s + p_s)/s  ... this has a zero at origin making Z_s → ∞ at DC
    //
    // Let's use proper admittances instead:
    // Source admittance: Y_s(s) = 1/R_s * p_s/(s + p_s) = (p_s/R_s)/(s + p_s)
    //   => Z_s(s) = R_s * (s + p_s)/p_s = R_s * (1 + s/p_s)
    //   At DC: Z_s = R_s. At high freq: Z_s → ∞
    //
    // Wait, that makes Z_s improper. Let's just use simple resistive admittances:
    // Y_s(s) = (1/R_s) / (s*tau_s + 1), R_s = 10, tau_s = 1/1000
    //   G(s) = (1/R_s) * (1/tau_s) / (s + 1/tau_s) = 100 / (s + 1000)
    //   => Z_s = R_s * (1 + s/1000)... no, Z = 1/Y = R_s*(s+1000)/1000... that's not right.
    //
    // Let me think again. Y(s) = 100/(s+1000). At DC: Y(0) = 0.1, so Z(0) = 10. ✓
    // At high freq: Y→0, Z→∞ (capacitor-like). That's fine for testing.

    constexpr double R_s = 10.0; // DC impedance
    constexpr double p_s = 1000.0;
    constexpr auto   A_ys = Matrix<1, 1, double>{{-p_s}};
    constexpr auto   B_ys = Matrix<1, 1, double>{{1.0}};
    constexpr auto   C_ys = Matrix<1, 1, double>{{p_s / R_s}}; // = 100
    constexpr auto   D_ys = Matrix<1, 1, double>{{0.0}};
    constexpr auto   Ys_sys = StateSpace<1, 1, 1>{A_ys, B_ys, C_ys, D_ys};

    // Load admittance: Y_L = 1/(s*tau_L + 1) / R_L, R_L = 1000
    // G(s) = (p_L/R_L)/(s+p_L) = 1/(s+500)
    constexpr double R_L = 1000.0;
    constexpr double p_L = 500.0;
    constexpr auto   A_yL = Matrix<1, 1, double>{{-p_L}};
    constexpr auto   B_yL = Matrix<1, 1, double>{{1.0}};
    constexpr auto   C_yL = Matrix<1, 1, double>{{p_L / R_L}}; // = 0.5
    constexpr auto   D_yL = Matrix<1, 1, double>{{0.0}};
    constexpr auto   YL_sys = StateSpace<1, 1, 1>{A_yL, B_yL, C_yL, D_yL};

    auto omega = geomspace(0.1, 100000.0, 100);

    auto result = middlebrook(impedance_from_admittance(Ys_sys, omega), impedance_from_admittance(YL_sys, omega));

    // At DC: Z_s = R_s = 10, Z_L = R_L = 1000
    // T_m = Z_s/Z_L = 0.01 => -40 dB
    CHECK(result.minor_loop_gain.points[0].magnitude_db == doctest::Approx(-40.0).epsilon(0.5));

    // Source impedance at DC
    CHECK(result.source_impedance.points[0].magnitude == doctest::Approx(R_s).epsilon(0.1));

    // Load impedance at DC
    CHECK(result.load_impedance.points[0].magnitude == doctest::Approx(R_L).epsilon(1.0));

    // System should be stable (sufficient condition)
    CHECK(result.is_stable_sufficient());
}

TEST_CASE("impedance: discrete plant evaluates on unit circle z=e^{jωTs}") {
    // Constant Z = R as a discrete system (Ts > 0). Magnitude must match continuous
    // D-only path; using jω on a discrete A would not be the intended evaluation.
    constexpr double R = 7.0;
    constexpr double Ts = 0.001;
    constexpr auto   A = Matrix<1, 1, double>{{-0.5}};
    constexpr auto   B = Matrix<1, 1, double>{{0.0}};
    constexpr auto   C = Matrix<1, 1, double>{{0.0}};
    constexpr auto   D = Matrix<1, 1, double>{{R}};
    auto             sys = StateSpace<1, 1, 1>{A, B, C, D};
    sys.Ts = Ts;
    REQUIRE(sys.is_discrete());

    auto omega = geomspace(1.0, 100.0, 20);
    auto result = impedance(sys, omega);
    REQUIRE(result.points.size() == 20);
    for (const auto& pt : result.points) {
        CHECK(pt.magnitude == doctest::Approx(R).epsilon(1e-10));
        CHECK(pt.phase_deg == doctest::Approx(0.0).epsilon(1e-6));
    }
}

TEST_CASE("Middlebrook: gain and phase margin delegation") {
    // Create a system where we can compute known margins
    // Source: Z_s = 1 (unit resistance, constant)
    constexpr auto A_s = Matrix<1, 1, double>{{-1.0}};
    constexpr auto B_s = Matrix<1, 1, double>{{0.0}};
    constexpr auto C_s = Matrix<1, 1, double>{{0.0}};
    constexpr auto D_s = Matrix<1, 1, double>{{1.0}};
    constexpr auto Zs_sys = StateSpace<1, 1, 1>{A_s, B_s, C_s, D_s};

    // Load: Z_L = 10 (constant, 10x bigger)
    constexpr auto A_L = Matrix<1, 1, double>{{-1.0}};
    constexpr auto B_L = Matrix<1, 1, double>{{0.0}};
    constexpr auto C_L = Matrix<1, 1, double>{{0.0}};
    constexpr auto D_L = Matrix<1, 1, double>{{10.0}};
    constexpr auto ZL_sys = StateSpace<1, 1, 1>{A_L, B_L, C_L, D_L};

    auto omega = geomspace(0.1, 100000.0, 100);
    auto Zs_data = impedance(Zs_sys, omega);
    auto ZL_data = impedance(ZL_sys, omega);
    auto result = middlebrook(Zs_data, ZL_data);

    // Tm = 1/10 = 0.1 = -20 dB at all frequencies, phase = 0°
    // With constant magnitude below 0 dB, there's no 0 dB crossing
    // so phase_margin() should return nullopt
    auto pm = result.phase_margin();
    CHECK_FALSE(pm.has_value());

    // Gain margin: no -180° crossing either
    auto gm = result.gain_margin();
    CHECK_FALSE(gm.has_value());

    // Sufficient condition
    CHECK(result.is_stable_sufficient());
}

TEST_CASE("Middlebrook from impedance data pairs samples by frequency") {
    analysis::ImpedanceResult<double> zs, zl;
    for (const double w : {1.0, 2.0, 3.0}) {
        zs.points.push_back({w, {1.0, 0.0}, 1.0, 0.0, 0.0});
    }
    // Load sweep lost its w = 2 sample (contour pole skipped).
    zl.points.push_back({1.0, {10.0, 0.0}, 10.0, 20.0, 0.0});
    zl.points.push_back({3.0, {4.0, 0.0}, 4.0, 12.0, 0.0});
    const auto mb = analysis::middlebrook(zs, zl);
    REQUIRE(mb.minor_loop_gain.points.size() == 2);
    CHECK(mb.minor_loop_gain.points[1].omega == doctest::Approx(3.0));
    CHECK(mb.minor_loop_gain.points[1].magnitude == doctest::Approx(0.25));
}
