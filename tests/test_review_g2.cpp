// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// Regression tests for the design / systems review (group g2). Most cases pin
// scale invariance: rescaling units (time, signal magnitude) must not change
// a structural answer or turn a valid design into a failure.

#include <cmath>
#include <cstddef>

#include "damp/backend.hpp"
#include "damp/design/linearization.hpp"
#include "damp/design/lyapunov.hpp"
#include "damp/design/minreal.hpp"
#include "damp/design/model_reduction.hpp"
#include "damp/design/pid_design.hpp"
#include "damp/design/qp.hpp"
#include "damp/design/stability.hpp"
#include "damp/math/complex.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/systems/zpk.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

TEST_SUITE("Review g2") {
    // ------------------------------------------------------------------ stability

    TEST_CASE("controllability_gramian rejects an unstable A") {
        // A X + X Aᵀ + B Bᵀ = 0 has the unique solution X = −0.5 for A = 1, which
        // is not a Gramian (the integral diverges). It must not be reported.
        const Matrix<1, 1> A{{1.0}};
        const Matrix<1, 1> B{{1.0}};
        CHECK_FALSE(stability::controllability_gramian(A, B).has_value());
        CHECK_FALSE(stability::observability_gramian(A, B).has_value());
        const Matrix<1, 1> Ad{{2.0}};
        CHECK_FALSE(stability::controllability_gramian(Ad, B, true).has_value());
    }

    TEST_CASE("is_controllable / is_observable are invariant to signal scaling") {
        const Matrix<2, 2> A{{-1.0, 0.0}, {0.0, -2.0}};
        const Matrix<2, 1> B{{1e-11}, {1e-11}};
        const Matrix<1, 2> C{{1e-11, 1e-11}};
        CHECK(stability::is_controllable(A, B));
        CHECK(stability::is_observable(A, C));
    }

    TEST_CASE("is_controllable is invariant to the time unit") {
        // Same plant in seconds and in milliseconds⁻¹: A scales by 1e3, so the
        // Krylov columns grow as 1, 1e3, 1e6, 1e9.
        Matrix<4, 4> A{};
        for (size_t i = 0; i < 4; ++i) {
            A(i, i) = -1e3 * static_cast<double>(i + 1);
        }
        const Matrix<4, 1> B{{1.0}, {1.0}, {1.0}, {1.0}};
        CHECK(stability::is_controllable(A, B));
        const StateSpace<4, 1, 1> sys{.A = A, .B = B, .C = B.transpose()};
        CHECK(design::minreal(sys).order == 4);
    }

    // ------------------------------------------------------------------ linearization

    TEST_CASE("linearize is accurate at a large operating point") {
        const auto f = [](const ColVec<1>& x, const ColVec<1>& u) {
            return ColVec<1>{(x(0) * x(0)) + (u(0) * u(0))};
        };
        const auto lin = design::linearize<1, 1>(f, ColVec<1>{1e10}, ColVec<1>{-3e10});
        CHECK(lin.A(0, 0) == doctest::Approx(2e10).epsilon(1e-6));
        CHECK(lin.B(0, 0) == doctest::Approx(-6e10).epsilon(1e-6));
    }

    // ------------------------------------------------------------------ PID rules

    TEST_CASE("AMIGO PI matches Åström & Hägglund (2006) §7.6") {
        constexpr double Ks = 1.0;
        constexpr double Ku = 2.0;
        constexpr double Tu = 1.0;
        const double     wu = 2.0 * damp::numbers::pi_v<double> / Tu;
        const double     T = std::sqrt((Ks * Ks * Ku * Ku) - 1.0) / wu;
        const double     L = (damp::numbers::pi_v<double> - std::atan(wu * T)) / wu;
        const double     K = (0.15 / Ks) + ((0.35 - (L * T / ((L + T) * (L + T)))) * T / (Ks * L));
        const double     Ti = (0.35 * L) + (13.0 * L * T * T / ((T * T) + (12.0 * L * T) + (7.0 * L * L)));
        const auto       pi = design::amigo_kappa_tau(Ks, Ku, Tu);
        CHECK(pi.Kp == doctest::Approx(K).epsilon(1e-12));
        CHECK(pi.Ki == doctest::Approx(K / Ti).epsilon(1e-12));
        CHECK(pi.Kd == 0.0);
    }

    TEST_CASE("Cohen–Coon PD matches the classic table") {
        // Kc = (τ/(K L))(5/4 + r/6), Td = L (6 − 2r)/(22 + 3r), r = L/τ
        const double K = 2.0;
        const double L = 1.0;
        const double tau = 4.0;
        const double r = L / tau;
        const double Kp = (tau / (K * L)) * (1.25 + (r / 6.0));
        const double Td = L * (6.0 - (2.0 * r)) / (22.0 + (3.0 * r));
        const auto   pd = design::cohen_coon(K, L, tau, 0.01, design::PIDType::PD);
        CHECK(pd.Kp == doctest::Approx(Kp).epsilon(1e-12));
        CHECK(pd.Kd == doctest::Approx(Kp * Td).epsilon(1e-12));
    }

    TEST_CASE("phase_margin_from_damping_ratio is continuous at ζ = 1") {
        const double below = design::phase_margin_from_damping_ratio(0.999999);
        const double at = design::phase_margin_from_damping_ratio(1.0);
        CHECK(at == doctest::Approx(below).epsilon(1e-4));
        // Exact second-order value: atan(2 / √(√5 − 2)) ≈ 76.35°
        CHECK(at == doctest::Approx(76.345).epsilon(1e-4));
        CHECK(design::phase_margin_from_damping_ratio(2.0) > at);
    }

    TEST_CASE("pid_from_bandwidth rejects a non-positive bandwidth") {
        for (const auto type : {design::PIDType::PI, design::PIDType::PD, design::PIDType::PID}) {
            const auto r = design::pid_from_bandwidth(0.0, 60.0, 0.01, type);
            CHECK(r.Kp == 0.0);
            CHECK(r.Ki == 0.0);
            CHECK(r.Kd == 0.0);
        }
    }

    TEST_CASE("ziegler_nichols_step rejects a non-positive time constant") {
        const auto r = design::ziegler_nichols_step(1.0, 0.5, -2.0, 0.01);
        CHECK(r.Kp == 0.0);
        CHECK(r.Ki == 0.0);
        CHECK(r.Kd == 0.0);
    }

    // ------------------------------------------------------------------ discretization

    TEST_CASE("ZOH is accurate for a nearly-singular A") {
        // B_d = (e^{a Ts} − 1)/a ≈ Ts; the A⁻¹(e^{ATs} − I) form cancels.
        const StateSpace<1, 1, 1> sys{.A = {{-1e-12}}, .B = {{1.0}}, .C = {{1.0}}};
        const auto                d = discretize(sys, 1e-3, DiscretizationMethod::ZOH);
        REQUIRE(d.has_value());
        CHECK(d->B(0, 0) == doctest::Approx(1e-3).epsilon(1e-12));
    }

    // ------------------------------------------------------------------ TF / ZPK

    TEST_CASE("TF realization accepts a uniformly scaled polynomial") {
        // 1e-13/(1e-13 s + 1e-13) = 1/(s + 1)
        const TransferFunction<1, 2> tf{.num = {1e-13}, .den = {1e-13, 1e-13}};
        const auto                   ss = tf.to_state_space();
        REQUIRE(ss.has_value());
        CHECK(ss->A(0, 0) == doctest::Approx(-1.0));
        CHECK((ss->C(0, 0) * ss->B(0, 0)) == doctest::Approx(1.0));
        const auto z = to_zpk(tf);
        REQUIRE(z.success);
        CHECK(z.zpk.gain == doctest::Approx(1.0));
    }

    TEST_CASE("SS→ZPK keeps a finite zero of a small-gain plant") {
        // G = 1e-14 (s + 2)/((s + 1)(s + 3))
        const StateSpace<2, 1, 1> sys{
            .A = {{0.0, 1.0}, {-3.0, -4.0}},
            .B = {{0.0}, {1.0}},
            .C = {{2e-14, 1e-14}},
        };
        const auto z = to_zpk(sys);
        REQUIRE(z.success);
        REQUIRE(z.n_zeros == 1);
        CHECK(z.zeros[0].real() == doctest::Approx(-2.0));
        CHECK(z.gain == doctest::Approx(1e-14));
    }

    TEST_CASE("ZPK dcgain is defined for a slow but nonzero pole") {
        // H = 1/(s + 1e-13): DC gain 1e13
        const ZPK<0, 1> sys{.poles = {damp::complex<double>{-1e-13, 0.0}}, .gain = 1.0};
        const auto      dc = sys.dcgain();
        REQUIRE(dc.has_value());
        CHECK(*dc == doctest::Approx(1e13));
    }

    TEST_CASE("ZPK cancel_matching is relative to the root magnitude") {
        // A right-half-plane zero must not cancel a left-half-plane pole just
        // because both are near the origin.
        const ZPK<1, 1> near_origin{
            .zeros = {damp::complex<double>{1e-13, 0.0}},
            .poles = {damp::complex<double>{-1e-13, 0.0}},
        };
        CHECK(near_origin.cancel_matching().n_cancelled == 0);
        // Identical to 12 digits at 1 MHz-scale: a numerical common factor.
        const ZPK<1, 1> fast{
            .zeros = {damp::complex<double>{-1e6 * (1.0 + 1e-12), 0.0}},
            .poles = {damp::complex<double>{-1e6, 0.0}},
        };
        CHECK(fast.cancel_matching().n_cancelled == 1);
    }

    // ------------------------------------------------------------------ interconnections

    TEST_CASE("series does not route measurement noise into process-noise columns") {
        // u → a → b. a's measurement noise v₁ is not a process-noise input w₁.
        const StateSpace<1, 1, 1, double, 1, 2> a{
            .A = {{-1.0}},
            .B = {{1.0}},
            .C = {{1.0}},
            .G = {{2.0}},
            .H = {{0.3, 0.7}},
        };
        const StateSpace<1, 1, 1, double, 1, 1> b{
            .A = {{-2.0}},
            .B = {{1.0}},
            .C = {{1.0}},
            .G = {{4.0}},
            .H = {{0.5}},
        };
        const auto ser = series(a, b);
        REQUIRE(ser.has_value());
        CHECK(ser->G(0, 0) == 2.0);
        CHECK(ser->G(1, 0) == 0.0); // w₁ only drives x₁
        CHECK(ser->G(1, 1) == 4.0);

        const auto fb = feedback(a, b);
        REQUIRE(fb.has_value());
        CHECK(fb->G(1, 0) == 0.0);
        CHECK(fb->H(0, 0) == doctest::Approx(0.3)); // D₁ = 0: y = C₁x₁ + H₁v₁
        CHECK(fb->H(0, 2) == 0.0);                  // v₂ reaches y only through x₁
    }

    TEST_CASE("sample-time match is relative") {
        const StateSpace<1, 1, 1> a{.A = {{0.5}}, .B = {{1.0}}, .C = {{1.0}}, .Ts = 1e-13};
        const StateSpace<1, 1, 1> b{.A = {{0.5}}, .B = {{1.0}}, .C = {{1.0}}, .Ts = 2e-13};
        CHECK_FALSE(series(a, b).has_value());
    }

    TEST_CASE("mix accepts a high but finite algebraic-loop gain") {
        // v = u + y with D = 1 − 1e-10: 1 − D is small but not lost to roundoff.
        const StateSpace<1, 1, 1> P{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}, .D = {{1.0 - 1e-13}}};
        const auto                m = mix(P, 0, 1.0, damp::array<SsFeed<double>, 1>{SsFeed<double>{0, 1.0}});
        CHECK(m.has_value());
    }

    // ------------------------------------------------------------------ QP

    TEST_CASE("solve_qp is invariant to scaling the cost") {
        // MATLAB quadprog example; scaling H and f by s leaves the minimizer.
        for (const double s : {1e-14, 1.0, 1e14}) {
            const Matrix<2, 2> H{{s, -s}, {-s, 2.0 * s}};
            const ColVec<2>    f{-2.0 * s, -6.0 * s};
            const Matrix<3, 2> A{{1.0, 1.0}, {-1.0, 2.0}, {2.0, 1.0}};
            const ColVec<3>    b{2.0, 2.0, 3.0};
            const auto         res = design::solve_qp(H, f, A, b);
            CAPTURE(s);
            REQUIRE(res.success);
            CHECK(res.x(0) == doctest::Approx(2.0 / 3.0).epsilon(1e-9));
            CHECK(res.x(1) == doctest::Approx(4.0 / 3.0).epsilon(1e-9));
        }
    }

    TEST_CASE("solve_qp is invariant to scaling the constraint rows") {
        const Matrix<2, 2> H{{1.0, -1.0}, {-1.0, 2.0}};
        const ColVec<2>    f{-2.0, -6.0};
        for (const double s : {1e-9, 1.0, 1e9}) {
            const Matrix<3, 2> A{{s, s}, {-s, 2.0 * s}, {2.0 * s, s}};
            const ColVec<3>    b{2.0 * s, 2.0 * s, 3.0 * s};
            const auto         res = design::solve_qp(H, f, A, b);
            CAPTURE(s);
            REQUIRE(res.success);
            CHECK(res.x(0) == doctest::Approx(2.0 / 3.0).epsilon(1e-9));
            CHECK(res.x(1) == doctest::Approx(4.0 / 3.0).epsilon(1e-9));
        }
    }

    TEST_CASE("interior-point QP is invariant to scaling the cost") {
        for (const double s : {1e-6, 1.0, 1e9}) {
            const Matrix<2, 2> H{{s, -s}, {-s, 2.0 * s}};
            const ColVec<2>    f{-2.0 * s, -6.0 * s};
            const Matrix<3, 2> A{{1.0, 1.0}, {-1.0, 2.0}, {2.0, 1.0}};
            const ColVec<3>    b{2.0, 2.0, 3.0};
            const auto         res = design::solve_qp_interior_point(H, f, A, b);
            CAPTURE(s);
            REQUIRE(res.success);
            CHECK(res.x(0) == doctest::Approx(2.0 / 3.0).epsilon(1e-6));
            CHECK(res.x(1) == doctest::Approx(4.0 / 3.0).epsilon(1e-6));
        }
    }
}
