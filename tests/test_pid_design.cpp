// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <array>
#include <cmath>
#include <numbers>

#include "damp/backend.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/design/pid_design.hpp"
#include "damp/design/status.hpp"
#include "damp/math/complex.hpp"
#include "damp/matlab.hpp"
#include "damp/matrix/core.hpp"
#include "damp/systems/realization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

/// s-plane pole whose ZOH image at Ts is the real z-plane pole z: s = ln(z)/Ts.
constexpr damp::complex<double> s_pole(double z, double Ts) {
    return {damp::log(z) / Ts, 0.0};
}

/// Closed-loop characteristic polynomial (highest power first) of b₀/(a₂ s² + a₁ s + a₀) under
/// Kp + Kᵢ/s + Kd s/(1 + Tf s): s(1 + Tf s)(a₂ s² + a₁ s + a₀) + b₀(Kp s(1 + Tf s) + Kᵢ(1 + Tf s) + Kd s²).
std::array<double, 5> pidf_closed_loop(const TransferFunction<1, 3>& plant, const design::PIDGains<double>& g) {
    const double b0 = plant.num[0];
    const double a2 = plant.den[0];
    const double a1 = plant.den[1];
    const double a0 = plant.den[2];
    return {
        a2 * g.Tf,
        a2 + (a1 * g.Tf),
        a1 + (a0 * g.Tf) + (b0 * ((g.Kp * g.Tf) + g.Kd)),
        a0 + (b0 * (g.Kp + (g.Ki * g.Tf))),
        b0 * g.Ki,
    };
}

/// Monic (s² + 2ζω s + ω²)(s + ωᵢ)(s + ω_f).
std::array<double, 5> placed_poly(double w, double zeta, double wi, double wf) {
    const double c1 = (2.0 * zeta * w) + wi;
    const double c2 = (w * w) + (2.0 * zeta * w * wi);
    const double c3 = w * w * wi;
    return {1.0, c1 + wf, c2 + (c1 * wf), c3 + (c2 * wf), c3 * wf};
}

/// The closed loop of the delivered gains (filter included) has the placed poles and the filter pole at −10ω.
void check_placed(const TransferFunction<1, 3>& plant, const design::PIDGains<double>& g, double w, double zeta, double wi) {
    const auto cl = pidf_closed_loop(plant, g);
    const auto want = placed_poly(w, zeta, wi, 10.0 * w);
    REQUIRE(cl[0] > 0.0);
    for (size_t i = 0; i < 5; ++i) {
        CHECK(cl[i] / cl[0] == doctest::Approx(want[i]).epsilon(1e-10));
    }
}

} // namespace

/**
 * @brief Tests for modelless PID design methods (pid_design.hpp)
 *
 * Validates Ziegler-Nichols, Cohen-Coon, SIMC, Lambda, pidtune,
 * Tyreus-Luyben, and pole placement PID tuning methods.
 */

TEST_SUITE("PID Design - Ziegler-Nichols Ultimate Gain") {
    TEST_CASE("ZN PID tuning from Ku and Tu") {
        constexpr double Ku = 10.0;
        constexpr double Tu = 2.0;

        constexpr auto result = design::ziegler_nichols(Ku, Tu);
        static_assert(result.Kp == 0.6 * Ku);             // 6.0
        static_assert(result.Ki == result.Kp / 1.0);      // Kp / (Tu/2) = 6.0
        static_assert(result.Kd == result.Kp * Tu / 8.0); // 6.0 * 0.25 = 1.5

        CHECK(result.Kp == doctest::Approx(6.0));
        CHECK(result.Ki == doctest::Approx(6.0));
        CHECK(result.Kd == doctest::Approx(1.5));
    }

    TEST_CASE("ZN PI tuning") {
        constexpr double Ku = 10.0;
        constexpr double Tu = 2.0;

        constexpr auto result = design::ziegler_nichols(Ku, Tu, PIDMode::PI);
        CHECK(result.Kp == doctest::Approx(4.5)); // 0.45 * 10
        CHECK(result.Ki == doctest::Approx(4.5 / (2.0 / 1.2)));
        CHECK(result.Kd == doctest::Approx(0.0));
    }

    TEST_CASE("ZN P-only tuning") {
        constexpr auto result = design::ziegler_nichols(10.0, 2.0, PIDMode::P);
        CHECK(result.Kp == doctest::Approx(5.0));
        CHECK(result.Ki == doctest::Approx(0.0));
        CHECK(result.Kd == doctest::Approx(0.0));
    }

    TEST_CASE("Runtime ZN matches compile-time") {
        auto ct = design::ziegler_nichols(10.0, 2.0, PIDMode::PID);
        auto rt = design::ziegler_nichols(10.0, 2.0, PIDMode::PID);
        CHECK(ct.Kp == doctest::Approx(rt.Kp));
        CHECK(ct.Ki == doctest::Approx(rt.Ki));
        CHECK(ct.Kd == doctest::Approx(rt.Kd));
    }

    TEST_CASE("ZN Tu=0 returns finite zeros for PI/PID") {
        const auto pi = design::ziegler_nichols(10.0, 0.0, PIDMode::PI);
        CHECK(std::isfinite(pi.Kp));
        CHECK(std::isfinite(pi.Ki));
        CHECK(std::isfinite(pi.Kd));
        CHECK(pi.Kp == doctest::Approx(0.0));
        CHECK(pi.Ki == doctest::Approx(0.0));
        CHECK(pi.Kd == doctest::Approx(0.0));

        const auto pid = design::ziegler_nichols(10.0, 0.0, PIDMode::PID);
        CHECK(std::isfinite(pid.Kp));
        CHECK(pid.Kp == doctest::Approx(0.0));
        CHECK(pid.Ki == doctest::Approx(0.0));
        CHECK(pid.Kd == doctest::Approx(0.0));
    }

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
}

TEST_SUITE("PID Design - Ziegler-Nichols Step Response") {
    TEST_CASE("ZN step response PID") {
        constexpr double K = 1.0;
        constexpr double L = 0.5;
        constexpr double tau = 2.0;

        constexpr auto result = design::ziegler_nichols_step(K, L, tau);
        // Kp = 1.2 * tau/(K*L) = 1.2 * 4.0 = 4.8
        CHECK(result.Kp == doctest::Approx(4.8));
        // Ki = Kp / (2*L) = 4.8 / 1.0 = 4.8
        CHECK(result.Ki == doctest::Approx(4.8));
        // Kd = Kp * 0.5 * L = 4.8 * 0.25 = 1.2
        CHECK(result.Kd == doctest::Approx(1.2));
    }

    TEST_CASE("ZN step K=0 or L=0 returns finite zeros") {
        const auto k0 = design::ziegler_nichols_step(0.0, 0.5, 2.0);
        CHECK(std::isfinite(k0.Kp));
        CHECK(std::isfinite(k0.Ki));
        CHECK(std::isfinite(k0.Kd));
        CHECK(k0.Kp == doctest::Approx(0.0));
        CHECK(k0.Ki == doctest::Approx(0.0));
        CHECK(k0.Kd == doctest::Approx(0.0));

        const auto l0 = design::ziegler_nichols_step(1.0, 0.0, 2.0);
        CHECK(std::isfinite(l0.Kp));
        CHECK(std::isfinite(l0.Ki));
        CHECK(std::isfinite(l0.Kd));
        CHECK(l0.Kp == doctest::Approx(0.0));
        CHECK(l0.Ki == doctest::Approx(0.0));
        CHECK(l0.Kd == doctest::Approx(0.0));
    }

    TEST_CASE("ziegler_nichols_step refuses PD, which the reaction-curve table does not define") {
        CHECK(design::ziegler_nichols_step(1.0, 0.5, 2.0, PIDMode::PD).status == design::DesignStatus::invalid_parameter);
    }

    TEST_CASE("ziegler_nichols_step rejects a non-positive time constant") {
        const auto r = design::ziegler_nichols_step(1.0, 0.5, -2.0);
        CHECK(r.Kp == 0.0);
        CHECK(r.Ki == 0.0);
        CHECK(r.Kd == 0.0);
    }
}

TEST_SUITE("PID Design - Tyreus-Luyben") {
    TEST_CASE("TL PID tuning") {
        constexpr double Ku = 10.0;
        constexpr double Tu = 2.0;

        constexpr auto result = design::tyreus_luyben(Ku, Tu);
        CHECK(result.Kp == doctest::Approx(Ku / 2.2));
        CHECK(result.Ki == doctest::Approx(result.Kp / (2.2 * Tu)));
        CHECK(result.Kd == doctest::Approx(result.Kp * Tu / 6.3));
    }

    TEST_CASE("TL PI tuning") {
        constexpr auto result = design::tyreus_luyben(10.0, 2.0, PIDMode::PI);
        CHECK(result.Kp == doctest::Approx(10.0 / 3.2));
    }

    TEST_CASE("TL refuses the P and PD modes it does not define") {
        CHECK(design::tyreus_luyben(10.0, 2.0, PIDMode::P).status == design::DesignStatus::invalid_parameter);
        CHECK(design::tyreus_luyben(10.0, 2.0, PIDMode::PD).status == design::DesignStatus::invalid_parameter);
    }

    TEST_CASE("TL more conservative than ZN") {
        constexpr auto zn = design::ziegler_nichols(10.0, 2.0);
        constexpr auto tl = design::tyreus_luyben(10.0, 2.0);
        // Tyreus-Luyben should have lower Kp (more conservative)
        CHECK(tl.Kp < zn.Kp);
    }
}

TEST_SUITE("PID Design - Cohen-Coon") {
    TEST_CASE("Cohen-Coon PID for FOPDT") {
        constexpr double K = 1.0;
        constexpr double L = 1.0;
        constexpr double tau = 4.0;

        constexpr auto result = design::cohen_coon(K, L, tau);
        // r = L/tau = 0.25
        // a = tau/(K*L) = 4.0
        // Kp = a * (4/3 + r/4) = 4.0 * (1.333 + 0.0625) = 5.583
        double r = L / tau;
        double a = tau / (K * L);
        double expected_Kp = a * ((4.0 / 3.0) + (r / 4.0));
        CHECK(result.Kp == doctest::Approx(expected_Kp).epsilon(1e-10));
        CHECK(result.Ki > 0.0);
        CHECK(result.Kd > 0.0);
    }

    TEST_CASE("Runtime Cohen-Coon matches compile-time") {
        auto ct = design::cohen_coon(1.0, 1.0, 4.0);
        auto rt = design::cohen_coon(1.0, 1.0, 4.0);
        CHECK(ct.Kp == doctest::Approx(rt.Kp).epsilon(1e-12));
        CHECK(ct.Ki == doctest::Approx(rt.Ki).epsilon(1e-12));
        CHECK(ct.Kd == doctest::Approx(rt.Kd).epsilon(1e-12));
    }

    TEST_CASE("Cohen-Coon K=0 or L=0 returns finite zeros") {
        const auto k0 = design::cohen_coon(0.0, 1.0, 4.0);
        CHECK(std::isfinite(k0.Kp));
        CHECK(k0.Kp == doctest::Approx(0.0));
        CHECK(k0.Ki == doctest::Approx(0.0));
        CHECK(k0.Kd == doctest::Approx(0.0));

        const auto l0 = design::cohen_coon(1.0, 0.0, 4.0);
        CHECK(std::isfinite(l0.Kp));
        CHECK(l0.Kp == doctest::Approx(0.0));
    }

    TEST_CASE("Cohen–Coon PD matches the classic table") {
        // Kc = (τ/(K L))(5/4 + r/6), Td = L (6 − 2r)/(22 + 3r), r = L/τ
        const double K = 2.0;
        const double L = 1.0;
        const double tau = 4.0;
        const double r = L / tau;
        const double Kp = (tau / (K * L)) * (1.25 + (r / 6.0));
        const double Td = L * (6.0 - (2.0 * r)) / (22.0 + (3.0 * r));
        const auto   pd = design::cohen_coon(K, L, tau, PIDMode::PD);
        CHECK(pd.Kp == doctest::Approx(Kp).epsilon(1e-12));
        CHECK(pd.Kd == doctest::Approx(Kp * Td).epsilon(1e-12));
    }

    TEST_CASE("Cohen–Coon PD refuses L/τ ≥ 3, where the table's Td is not positive") {
        CHECK(design::cohen_coon(1.0, 4.0, 1.0, PIDMode::PD).status == design::DesignStatus::invalid_parameter);
        CHECK(design::cohen_coon(1.0, 3.0, 1.0, PIDMode::PD).status == design::DesignStatus::invalid_parameter);
        CHECK(design::cohen_coon(1.0, 4.0, 1.0, PIDMode::PID).success); // the PID row stays positive
    }
}

TEST_SUITE("PID Design - SIMC") {
    TEST_CASE("SIMC PI tuning") {
        constexpr double K = 2.0;
        constexpr double L = 0.5;
        constexpr double tau = 3.0;
        constexpr double tau_c = 1.0;

        constexpr auto result = design::simc(K, L, tau, tau_c, PIDMode::PI);
        // Kp = tau / (K * (tau_c + L)) = 3.0 / (2.0 * 1.5) = 1.0
        CHECK(result.Kp == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(result.Ki > 0.0);
        CHECK(result.Kd == doctest::Approx(0.0));
    }

    TEST_CASE("SIMC K=0 returns finite zeros") {
        const auto z = design::simc(0.0, 0.5, 3.0, 1.0, PIDMode::PI);
        CHECK(std::isfinite(z.Kp));
        CHECK(std::isfinite(z.Ki));
        CHECK(z.Kp == doctest::Approx(0.0));
        CHECK(z.Ki == doctest::Approx(0.0));
    }
}

TEST_SUITE("PID Design - Lambda Tuning") {
    TEST_CASE("Lambda PI tuning") {
        constexpr double K = 1.0;
        constexpr double L = 0.5;
        constexpr double tau = 2.0;
        constexpr double lambda = 2.0; // Desired closed-loop time constant

        constexpr auto result = design::lambda_tuning(K, L, tau, lambda);
        // Kp = tau / (K * (lambda + L)) = 2.0 / (1.0 * 2.5) = 0.8
        CHECK(result.Kp == doctest::Approx(0.8).epsilon(1e-12));
        // Ki = Kp / tau = 0.8 / 2.0 = 0.4
        CHECK(result.Ki == doctest::Approx(0.4).epsilon(1e-12));
        CHECK(result.Kd == doctest::Approx(0.0));
    }

    TEST_CASE("Lambda K=0 or tau=0 returns finite zeros") {
        const auto k0 = design::lambda_tuning(0.0, 0.5, 2.0, 2.0);
        CHECK(std::isfinite(k0.Kp));
        CHECK(k0.Kp == doctest::Approx(0.0));
        CHECK(k0.Ki == doctest::Approx(0.0));

        const auto tau0 = design::lambda_tuning(1.0, 0.5, 0.0, 2.0);
        CHECK(std::isfinite(tau0.Kp));
        CHECK(tau0.Kp == doctest::Approx(0.0));
        CHECK(tau0.Ki == doctest::Approx(0.0));
    }
}

TEST_SUITE("PID Design - Margin helpers") {
    TEST_CASE("phase_margin_from_damping_ratio is continuous at ζ = 1") {
        const double below = design::phase_margin_from_damping_ratio(0.999999);
        const double at = design::phase_margin_from_damping_ratio(1.0);
        CHECK(at == doctest::Approx(below).epsilon(1e-4));
        // Exact second-order value: atan(2 / √(√5 − 2)) ≈ 76.35°
        CHECK(at == doctest::Approx(76.345).epsilon(1e-4));
        CHECK(design::phase_margin_from_damping_ratio(2.0) > at);
    }
}

TEST_SUITE("PID Design - Plant pidtune") {
    TEST_CASE("PI at ωc on 1/(s+1)") {
        constexpr StateSpace<1, 1, 1> sys{
            .A = Matrix<1, 1>{{-1.0}},
            .B = Matrix<1, 1>{{1.0}},
            .C = Matrix<1, 1>{{1.0}},
            .D = Matrix<1, 1>{{0.0}},
        };
        const auto pi = design::pidtune(sys, {.wc = 2.0, .phase_margin_deg = 60.0, .mode = PIDMode::PI});
        REQUIRE(pi.success);
        CHECK(pi.Kp > 0.0);
        CHECK(pi.Ki > 0.0);
        CHECK(pi.Kd == doctest::Approx(0.0));
        CHECK(pi.b == doctest::Approx(1.0));
        CHECK(pi.c == doctest::Approx(1.0));
    }

    TEST_CASE("PIDTuneSpec copies Tf and 2-DOF weights") {
        constexpr StateSpace<1, 1, 1> sys{
            .A = Matrix<1, 1>{{-1.0}},
            .B = Matrix<1, 1>{{1.0}},
            .C = Matrix<1, 1>{{1.0}},
            .D = Matrix<1, 1>{{0.0}},
        };
        design::PIDTuneSpec<double> spec{};
        spec.wc = 2.0;
        spec.phase_margin_deg = 70.0;
        spec.mode = PIDMode::PI;
        spec.Tf = 1e-3;
        spec.b = 0.0;
        spec.c = 0.0;
        const auto r = design::pidtune(sys, spec);
        REQUIRE(r.success);
        CHECK(r.Tf == doctest::Approx(1e-3));
        CHECK(r.b == doctest::Approx(0.0));
        CHECK(r.c == doctest::Approx(0.0));
        CHECK(r.Ki > 0.0);
    }

    TEST_CASE("1-DOF pidtune is unchanged when b=c=1") {
        constexpr StateSpace<1, 1, 1> sys{
            .A = Matrix<1, 1>{{-1.0}},
            .B = Matrix<1, 1>{{1.0}},
            .C = Matrix<1, 1>{{1.0}},
            .D = Matrix<1, 1>{{0.0}},
        };
        const auto                  a = design::pidtune(sys, {.wc = 2.0, .phase_margin_deg = 60.0, .mode = PIDMode::PI});
        design::PIDTuneSpec<double> spec{};
        spec.wc = 2.0;
        spec.phase_margin_deg = 60.0;
        spec.mode = PIDMode::PI;
        spec.b = 1.0;
        spec.c = 1.0;
        const auto b = design::pidtune(sys, spec);
        REQUIRE(a.success);
        REQUIRE(b.success);
        CHECK(a.Kp == doctest::Approx(b.Kp));
        CHECK(a.Ki == doctest::Approx(b.Ki));
        CHECK(a.Kd == doctest::Approx(b.Kd));
    }

    TEST_CASE("I-P pidtune matches 1-DOF T −3 dB") {
        constexpr StateSpace<1, 1, 1> sys{
            .A = Matrix<1, 1>{{-1.0}},
            .B = Matrix<1, 1>{{1.0}},
            .C = Matrix<1, 1>{{1.0}},
            .D = Matrix<1, 1>{{0.0}},
        };
        constexpr double            wc = 2.0;
        design::PIDTuneSpec<double> one{};
        one.wc = wc;
        one.phase_margin_deg = 60.0;
        one.mode = PIDMode::PI;
        one.b = 1.0;
        one.c = 1.0;
        design::PIDTuneSpec<double> ip = one;
        ip.b = 0.0;
        ip.c = 0.0;
        const auto r1 = design::pidtune(sys, one);
        const auto r0 = design::pidtune(sys, ip);
        REQUIRE(r1.success);
        REQUIRE(r0.success);
        CHECK(r0.Kp > r1.Kp);
        CHECK(r0.Ki > r1.Ki);

        const auto t_mag = [&](const design::PIDGains<double>& pid, double w) {
            const auto Gopt = eval_frf(sys, complex<double>{0.0, w});
            REQUIRE(Gopt.has_value());
            const auto P = (*Gopt)(0, 0);
            const auto L = pid.eval_fb(w) * P;
            const auto den = complex<double>{1.0, 0.0} + L;
            return abs((pid.eval_ff(w) * P) / den);
        };
        const auto t_bw = [&](const design::PIDGains<double>& pid) {
            const double dc = t_mag(pid, wc * 1e-3);
            const double thresh = dc / std::sqrt(2.0);
            double       lo = wc * 1e-3;
            double       hi = wc * 1e3;
            for (int k = 0; k < 20; ++k) {
                const double mid = std::sqrt(lo * hi);
                if (t_mag(pid, mid) >= thresh) {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            return std::sqrt(lo * hi);
        };
        const double bw1 = t_bw(r1);
        const double bw0 = t_bw(r0);
        CHECK(bw0 == doctest::Approx(bw1).epsilon(0.08));
        CHECK(bw1 == doctest::Approx(wc).epsilon(0.35));
    }

    TEST_CASE("PI-D pidtune on LC does not sign-flip") {
        constexpr double          L = 100e-6;
        constexpr double          C = 100e-6;
        constexpr double          R = 4.0;
        const StateSpace<2, 1, 1> sys{
            .A = Matrix<2, 2>{{0.0, -1.0 / L}, {1.0 / C, -1.0 / (R * C)}},
            .B = Matrix<2, 1>{{1.0 / L}, {0.0}},
            .C = Matrix<1, 2>{{0.0, 1.0}},
            .D = Matrix<1, 1>{{0.0}},
        };
        constexpr double            wc = 27650.638614364343;
        design::PIDTuneSpec<double> spec{};
        spec.wc = wc;
        spec.mode = PIDMode::PID;
        spec.b = 1.0;
        spec.c = 0.0;
        for (int k = 600; k <= 955; k += 5) {
            spec.phase_margin_deg = static_cast<double>(k) / 10.0;
            const auto tuned = design::pidtune(sys, spec);
            REQUIRE(tuned.success);
            CHECK(tuned.Kp > 0.0);
        }
        // The direct fit needs Kp < 0 here: refused rather than returned.
        spec.phase_margin_deg = 96.0;
        const auto past = design::pidtune(sys, spec);
        CHECK_FALSE(past.success);
        CHECK(past.status == design::DesignStatus::not_placeable);
    }

    TEST_CASE("pidtune refuses a law that cannot realize the crossover with positive gains") {
        using design::DesignStatus;
        // 1/(s(s+1)) at ωc = 1 needs phase lead for 60°: a PI would need Ki < 0.
        constexpr StateSpace<2, 1, 1> type1{.A = {{0.0, 1.0}, {0.0, -1.0}}, .B = {{0.0}, {1.0}}, .C = {{1.0, 0.0}}};
        CHECK(design::pidtune(type1, {.wc = 1.0, .mode = PIDMode::PI}).status == DesignStatus::not_placeable);
        // 1/(s+1)² at ωc = 2: the PI fit had Ki < 0 and an unstable closed loop.
        constexpr StateSpace<2, 1, 1> lag2{.A = {{-1.0, 1.0}, {0.0, -1.0}}, .B = {{0.0}, {1.0}}, .C = {{1.0, 0.0}}};
        CHECK(design::pidtune(lag2, {.wc = 2.0, .mode = PIDMode::PI}).status == DesignStatus::not_placeable);
        // 1/(s+1)³ at ωc = 3: every law needed Kp < 0.
        constexpr StateSpace<3, 1, 1> lag3{.A = {{-1.0, 1.0, 0.0}, {0.0, -1.0, 1.0}, {0.0, 0.0, -1.0}}, .B = {{0.0}, {0.0}, {1.0}}, .C = {{1.0, 0.0, 0.0}}};
        CHECK(design::pidtune(lag3, {.wc = 3.0, .mode = PIDMode::PI}).status == DesignStatus::not_placeable);
        CHECK(design::pidtune(lag3, {.wc = 3.0, .mode = PIDMode::PID}).status == DesignStatus::not_placeable);
        // 1/(s+1): a PD at ωc = 5 would need Kd < 0, a PID at ωc = 0.1 Kp < 0.
        constexpr StateSpace<1, 1, 1> lag1{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}};
        CHECK(design::pidtune(lag1, {.wc = 5.0, .mode = PIDMode::PD}).status == DesignStatus::not_placeable);
        CHECK(design::pidtune(lag1, {.wc = 0.1, .mode = PIDMode::PID}).status == DesignStatus::not_placeable);
    }

    TEST_CASE("pidtune refuses an unstable closed loop (P on 1/(s+1)³ above the ultimate frequency)") {
        // Kp = 1/|P(j3)| = 10^1.5 > 8, the ultimate gain of 1/(s+1)³.
        constexpr StateSpace<3, 1, 1> lag3{.A = {{-1.0, 1.0, 0.0}, {0.0, -1.0, 1.0}, {0.0, 0.0, -1.0}}, .B = {{0.0}, {0.0}, {1.0}}, .C = {{1.0, 0.0, 0.0}}};
        CHECK(design::pidtune(lag3, {.wc = 3.0, .mode = PIDMode::P}).status == design::DesignStatus::not_stabilizable);
        CHECK(design::pidtune(lag3, {.wc = 0.5, .mode = PIDMode::P}).success);
    }

    TEST_CASE("pidtune's certificate does not depend on the time unit") {
        using design::DesignStatus;
        // 1/(s+1)³ with time scaled by 1000 has every closed-loop root scaled by 1000. The Routh
        // test compared the leading coefficient with c[n] ~ ωc³ and refused the fast copy.
        constexpr StateSpace<3, 1, 1> lag3{.A = {{-1.0, 1.0, 0.0}, {0.0, -1.0, 1.0}, {0.0, 0.0, -1.0}}, .B = {{0.0}, {0.0}, {1.0}}, .C = {{1.0, 0.0, 0.0}}};
        constexpr double              a = 1000.0;
        constexpr StateSpace<3, 1, 1> fast{.A = {{-a, a, 0.0}, {0.0, -a, a}, {0.0, 0.0, -a}}, .B = {{0.0}, {0.0}, {a}}, .C = {{1.0, 0.0, 0.0}}};
        const auto                    slow_g = design::pidtune(lag3, {.wc = 0.5, .Tf = 0.2});
        const auto                    fast_g = design::pidtune(fast, {.wc = 500.0, .Tf = 2e-4});
        REQUIRE(slow_g.success);
        REQUIRE(fast_g.success);
        CHECK(fast_g.Kp == doctest::Approx(slow_g.Kp).epsilon(1e-9));
        CHECK(fast_g.Ki == doctest::Approx(a * slow_g.Ki).epsilon(1e-9));
        CHECK(fast_g.Kd == doctest::Approx(slow_g.Kd / a).epsilon(1e-9));
        // Past the ultimate gain (Kp = (1 + ωc²)^{3/2} > 8) both copies are refused.
        CHECK(design::pidtune(lag3, {.wc = 1.75, .mode = PIDMode::P}).status == DesignStatus::not_stabilizable);
        CHECK(design::pidtune(fast, {.wc = 1750.0, .mode = PIDMode::P}).status == DesignStatus::not_stabilizable);

        // A float PI at ωc = 1000 rad/s on 1/(s+1) was refused above ωc ≈ 200 rad/s.
        constexpr StateSpace<1, 1, 1, float> lag1f{.A = {{-1.0f}}, .B = {{1.0f}}, .C = {{1.0f}}};
        const auto                           f = design::pidtune(lag1f, design::PIDTuneSpec<float>{.wc = 1000.0f, .mode = PIDMode::PI});
        REQUIRE(f.success);
        CHECK(static_cast<double>(f.Kp) == doctest::Approx(-0.5 + (std::sqrt(3.0) / 2.0 * 1000.0)).epsilon(1e-5));
        // An R–L current loop 1000/(s + 100) at 1 kHz, in float.
        constexpr StateSpace<1, 1, 1, float> rl{.A = {{-100.0f}}, .B = {{1000.0f}}, .C = {{1.0f}}};
        CHECK(design::pidtune(rl, design::PIDTuneSpec<float>{.wc = 6283.0f, .mode = PIDMode::PI}).success);
    }

    TEST_CASE("pidtune on a reverse-acting plant returns the negated gains") {
        using design::DesignStatus;
        // −P needs −C for the same loop C·P: Kp, Ki, Kd ≤ 0, and Kbc keeps Ki's sign.
        constexpr StateSpace<1, 1, 1> lag1{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}};
        constexpr StateSpace<1, 1, 1> rev1{.A = lag1.A, .B = lag1.B, .C = {{-1.0}}};
        constexpr StateSpace<2, 1, 1> osc{.A = {{0.0, 1.0}, {-1.0, -0.5}}, .B = {{0.0}, {1.0}}, .C = {{1.0, 0.0}}};
        constexpr StateSpace<2, 1, 1> rev_osc{.A = osc.A, .B = osc.B, .C = {{-1.0, 0.0}}};
        const auto                    check_negated = [](const auto& direct, const auto& reverse, const design::PIDTuneSpec<double>& spec) {
            const auto g = design::pidtune(direct, spec);
            const auto r = design::pidtune(reverse, spec);
            REQUIRE(g.success);
            REQUIRE(r.success);
            CHECK(r.Kp == doctest::Approx(-g.Kp).epsilon(1e-12));
            CHECK(r.Ki == doctest::Approx(-g.Ki).epsilon(1e-12));
            CHECK(r.Kd == doctest::Approx(-g.Kd).epsilon(1e-12));
            CHECK(r.Kbc == doctest::Approx(-g.Kbc).epsilon(1e-12));
            CHECK(r.Kp < 0.0);
            CHECK(r.discretize(0.01).success);
        };
        check_negated(lag1, rev1, {.wc = 2.0, .mode = PIDMode::P});
        check_negated(lag1, rev1, {.wc = 2.0, .mode = PIDMode::PI});
        check_negated(lag1, rev1, {.wc = 2.0, .mode = PIDMode::PID});
        check_negated(lag1, rev1, {.wc = 2.0, .mode = PIDMode::PI, .b = 0.0}); // the 2-DOF crossover search
        check_negated(osc, rev_osc, {.wc = 2.0, .mode = PIDMode::PD});         // the plant lags 162° at ωc = 2: the PD leads
        check_negated(osc, rev_osc, {.wc = 2.0, .mode = PIDMode::PID});

        // A refusal mirrors too: 1/(s+1)³ at ωc = 3 is past the ultimate gain (P) and needs lead (PI).
        constexpr StateSpace<3, 1, 1> lag3{.A = {{-1.0, 1.0, 0.0}, {0.0, -1.0, 1.0}, {0.0, 0.0, -1.0}}, .B = {{0.0}, {0.0}, {1.0}}, .C = {{1.0, 0.0, 0.0}}};
        constexpr StateSpace<3, 1, 1> rev3{.A = lag3.A, .B = lag3.B, .C = {{-1.0, 0.0, 0.0}}};
        CHECK(design::pidtune(rev3, {.wc = 3.0, .mode = PIDMode::P}).status == DesignStatus::not_stabilizable);
        CHECK(design::pidtune(rev3, {.wc = 3.0, .mode = PIDMode::PI}).status == DesignStatus::not_placeable);
    }

    TEST_CASE("pidtune P meets the crossover: |Kp P(jωc)| = 1") {
        constexpr StateSpace<1, 1, 1> lag1{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}};
        const auto                    p = design::pidtune(lag1, {.wc = 2.0, .mode = PIDMode::P});
        REQUIRE(p.success);
        CHECK(p.Kp == doctest::Approx(std::sqrt(5.0)));
    }

    TEST_CASE("pidtune meets |L| = 1 and the phase margin with a derivative filter") {
        constexpr StateSpace<2, 1, 1> type1{.A = {{0.0, 1.0}, {0.0, -1.0}}, .B = {{0.0}, {1.0}}, .C = {{1.0, 0.0}}};
        constexpr StateSpace<3, 1, 1> lag3{.A = {{-1.0, 1.0, 0.0}, {0.0, -1.0, 1.0}, {0.0, 0.0, -1.0}}, .B = {{0.0}, {0.0}, {1.0}}, .C = {{1.0, 0.0, 0.0}}};
        for (const double Tf : {0.0, 0.1, 0.3}) {
            for (const auto mode : {PIDMode::PD, PIDMode::PID}) {
                const auto g = design::pidtune(type1, {.wc = 1.0, .mode = mode, .Tf = Tf});
                REQUIRE(g.success);
                const auto L = g.eval_fb(1.0) * (*eval_frf(type1, complex<double>{0.0, 1.0}))(0, 0);
                CHECK(abs(L) == doctest::Approx(1.0).epsilon(1e-12));
                CHECK(arg(L) * 180.0 / std::numbers::pi == doctest::Approx(-120.0).epsilon(1e-12));
            }
            const auto g = design::pidtune(lag3, {.wc = 1.0, .Tf = Tf});
            REQUIRE(g.success);
            CHECK(g.Ki == doctest::Approx(g.Kp * g.Kp / (4.0 * g.Kd))); // Ti = 4 Td
            const auto L = g.eval_fb(1.0) * (*eval_frf(lag3, complex<double>{0.0, 1.0}))(0, 0);
            CHECK(abs(L) == doctest::Approx(1.0).epsilon(1e-12));
            CHECK(arg(L) * 180.0 / std::numbers::pi == doctest::Approx(-120.0).epsilon(1e-12));
        }
    }

    TEST_CASE("pidtune refuses a non-finite or out-of-range spec") {
        using design::DesignStatus;
        constexpr StateSpace<1, 1, 1> lag1{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}};
        CHECK(design::pidtune(lag1, {.wc = 1.0, .phase_margin_deg = std::nan("")}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::pidtune(lag1, {.wc = 1.0, .phase_margin_deg = 0.0}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::pidtune(lag1, {.wc = 1.0, .phase_margin_deg = 180.0}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::pidtune(lag1, {.wc = 1.0, .Tf = -0.1}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::pidtune(lag1, {.wc = std::nan("")}).status == design::DesignStatus::invalid_parameter);
        constexpr StateSpace<1, 1, 1> sampled{.A = {{0.9}}, .B = {{1.0}}, .C = {{1.0}}, .Ts = 0.01};
        CHECK(design::pidtune(sampled, {.wc = 1.0}).status == DesignStatus::invalid_sample_time);
    }

    TEST_CASE("pidtune and pid_pole_placement accept a plant with noise channels") {
        constexpr StateSpace<2, 1, 1, double, 2, 1> noisy{
            .A = {{0.0, 1.0}, {0.0, -1.0}},
            .B = {{0.0}, {1.0}},
            .C = {{1.0, 0.0}},
            .G = Matrix<2, 2>::identity(),
            .H = Matrix<1, 1>::identity(),
        };
        constexpr StateSpace<2, 1, 1> plain{.A = noisy.A, .B = noisy.B, .C = noisy.C};
        const auto                    a = design::pidtune(noisy, {.wc = 1.0});
        REQUIRE(a.success);
        CHECK(a.Kp == design::pidtune(plain, {.wc = 1.0}).Kp);
        const auto p = design::pid_pole_placement(noisy, {.omega = 3.0});
        REQUIRE(p.success);
        CHECK(p.Kp == design::pid_pole_placement(plain, {.omega = 3.0}).Kp);
    }

    TEST_CASE("PIDGains to_tf_ff equals to_tf when b=c=1") {
        constexpr auto r = design::PIDGains{.Kp = 2.0, .Ki = 3.0, .Kd = 0.5, .Tf = 0.1};
        const auto     fb = r.to_tf();
        const auto     ff = r.to_tf_ff();
        for (size_t i = 0; i < 3; ++i) {
            CHECK(ff.num[i] == doctest::Approx(fb.num[i]));
            CHECK(ff.den[i] == doctest::Approx(fb.den[i]));
        }
        CHECK(abs(r.eval_ff(4.0) - r.eval_fb(4.0)) < 1e-12);
    }

    TEST_CASE("matlab::pidtune(sys, wc) is 60° PID") {
        constexpr StateSpace<2, 1, 1> sys{
            .A = Matrix<2, 2>{{0.0, 1.0}, {-1.0, -0.5}},
            .B = Matrix<2, 1>{{0.0}, {1.0}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .D = Matrix<1, 1>{{0.0}},
        };
        const auto a = matlab::pidtune(sys, 1.0);
        const auto b = design::pidtune(sys, {.wc = 1.0, .phase_margin_deg = 60.0, .mode = PIDMode::PID});
        REQUIRE(a.success);
        REQUIRE(b.success);
        CHECK(a.Kp == doctest::Approx(b.Kp));
        CHECK(a.Ki == doctest::Approx(b.Ki));
        CHECK(a.Kd == doctest::Approx(b.Kd));
    }

    TEST_CASE("PIDGains::to_ss realizes Kp + Ki/s + Kd s/(1 + Tf s)") {
        constexpr auto pi = design::pid(2.0, 4.0, 0.0).to_ss();
        static_assert(pi.has_value());
        static_assert(pi->A(0, 0) == 0.0);
        static_assert(pi->B(0, 0) == 1.0);
        static_assert(pi->C(0, 0) == 4.0);
        static_assert(pi->D(0, 0) == 2.0);
        constexpr auto L = *series(*pi, StateSpace<1, 1, 1>{
                                            .A = {{-1.0}},
                                            .B = {{1.0}},
                                            .C = {{1.0}},
                                        });
        CHECK(L.A.rows() == 3);

        // The filtered derivative used to be dropped: C(jω) must match to_tf at every ω.
        design::PIDGains<double> gains = design::pid(2.0, 4.0, 0.5);
        gains.Tf = 0.05;
        const auto C = gains.to_ss();
        REQUIRE(C.has_value());
        for (const double w : {0.1, 1.0, 10.0, 100.0}) {
            const auto ss = (*eval_frf(*C, damp::complex<double>{0.0, w}))(0, 0);
            const auto tf = gains.eval_fb(w);
            CHECK(ss.real() == doctest::Approx(tf.real()).epsilon(1e-12));
            CHECK(ss.imag() == doctest::Approx(tf.imag()).epsilon(1e-12));
        }
        // An unfiltered derivative is not proper.
        CHECK_FALSE(design::pid(2.0, 4.0, 0.5).to_ss().has_value());
    }
}

TEST_SUITE("PID Design - Pole Placement") {
    TEST_CASE("PI pole placement for first-order plant") {
        // Plant: G(s) = 1/(s+1), K=1, tau=1
        // Desired s-plane poles whose images are z=0.5 and z=0.3 (well inside unit circle)
        constexpr double K = 1.0;
        constexpr double tau = 1.0;
        constexpr double Ts = 0.1;
        constexpr double p1 = 0.5;
        constexpr double p2 = 0.3;

        constexpr auto result = design::pid_pole_placement(TransferFunction<1, 2>{.num = {K}, .den = {tau, 1.0}}, damp::array<damp::complex<double>, 2>{s_pole(p1, Ts), s_pole(p2, Ts)}, Ts);
        CHECK(result.Kp > 0.0);
        CHECK(result.Ki > 0.0);
        CHECK(result.Kd == doctest::Approx(0.0)); // PI only

        // Verify: closed-loop poles should be at p1, p2
        // Discretize plant: a = exp(-Ts/tau), b = K*(1-a)
        double a = std::exp(-Ts / tau);
        double b = K * (1.0 - a);

        // Backward-Euler PI: (z−a)(z−1) + b (Kp (z−1) + Ki Ts z)
        // = z² + (b Kp + b Ki Ts − a − 1) z + (a − b Kp)
        double c1 = (b * result.Kp) + (b * result.Ki * Ts) - a - 1.0;
        double c0 = a - (b * result.Kp);

        // Should equal desired: z² - (p1+p2)*z + p1*p2
        CHECK(c1 == doctest::Approx(-(p1 + p2)).epsilon(1e-10));
        CHECK(c0 == doctest::Approx(p1 * p2).epsilon(1e-10));
    }

    TEST_CASE("PID pole placement for first-order plant (3 poles)") {
        constexpr double K = 2.0;
        constexpr double tau = 0.5;
        constexpr double Ts = 0.01;
        constexpr double p1 = 0.8;
        constexpr double p2 = 0.7;
        constexpr double p3 = 0.6;

        constexpr auto result = design::pid_pole_placement(TransferFunction<1, 2>{.num = {K}, .den = {tau, 1.0}}, damp::array<damp::complex<double>, 3>{s_pole(p1, Ts), s_pole(p2, Ts), s_pole(p3, Ts)}, Ts);
        CHECK(result.Kp != 0.0);
        CHECK(result.Ki != 0.0);
        CHECK(result.Kd != 0.0);

        // Deploy C(z) = Kp + Ki Ts z/(z−1) + (Kd/Ts)(z−1)/z
        const double a = std::exp(-Ts / tau);
        const double b = K * (1.0 - a);
        const double c2 = -(1.0 + a) + (b * result.Kp) + (b * result.Ki * Ts) + (b * result.Kd / Ts);
        const double c1 = a - (b * result.Kp) - (2.0 * b * result.Kd / Ts);
        const double c0 = b * result.Kd / Ts;
        CHECK(c2 == doctest::Approx(-(p1 + p2 + p3)).epsilon(1e-10));
        CHECK(c1 == doctest::Approx((p1 * p2) + (p1 * p3) + (p2 * p3)).epsilon(1e-10));
        CHECK(c0 == doctest::Approx(-(p1 * p2 * p3)).epsilon(1e-10));
    }

    TEST_CASE("discrete placement takes s-plane poles, as place(sys, poles, Ts) does") {
        // s = −5 twice at Ts = 0.01 is z = e^(−0.05): a stable pair. Read as z-plane poles,
        // {−5, −5} gave an unstable loop with success.
        constexpr double Ts = 0.01;
        const auto       g = design::pid_pole_placement(TransferFunction<1, 2>{.num = {1.0}, .den = {1.0, 1.0}}, damp::array<damp::complex<double>, 2>{-5.0, -5.0}, Ts);
        REQUIRE(g.success);
        const double a = std::exp(-Ts);
        const double b = 1.0 - a;
        const double z = std::exp(-5.0 * Ts);
        CHECK(a - (b * g.Kp) == doctest::Approx(z * z).epsilon(1e-12));
        CHECK((1.0 + a) - (b * g.Kp) - (b * g.Ki * Ts) == doctest::Approx(2.0 * z).epsilon(1e-12));
    }

    TEST_CASE("discrete placement places a complex pair and refuses an unpaired complex pole") {
        constexpr double Ts = 0.01;
        constexpr auto   plant = TransferFunction<1, 2>{.num = {1.0}, .den = {1.0, 1.0}};
        const auto       g = design::pid_pole_placement(plant, damp::array<damp::complex<double>, 2>{complex<double>{-5.0, 3.0}, complex<double>{-5.0, -3.0}}, Ts);
        REQUIRE(g.success);
        // z₁z₂ = e^(−10 Ts), z₁ + z₂ = 2 e^(−5 Ts) cos(3 Ts)
        const double a = std::exp(-Ts);
        const double b = 1.0 - a;
        CHECK(a - (b * g.Kp) == doctest::Approx(std::exp(-10.0 * Ts)).epsilon(1e-12));
        CHECK((1.0 + a) - (b * g.Kp) - (b * g.Ki * Ts) == doctest::Approx(2.0 * std::exp(-5.0 * Ts) * std::cos(3.0 * Ts)).epsilon(1e-12));
        const auto lone = design::pid_pole_placement(plant, damp::array<damp::complex<double>, 2>{complex<double>{-5.0, 3.0}, complex<double>{-5.0, 3.0}}, Ts);
        CHECK(lone.status == design::DesignStatus::invalid_parameter);
        CHECK(design::pid_pole_placement(plant, damp::array<damp::complex<double>, 2>{std::nan(""), -5.0}, Ts).status == design::DesignStatus::invalid_parameter);
    }

    TEST_CASE("discrete placement refuses Re(s) ≥ 0, so a braced list of z-plane poles is refused") {
        // Read as s, {0.95, 0.96} is z = e^{0.95 Ts} > 1: an unstable loop that came back with success.
        constexpr auto plant = TransferFunction<1, 2>{.num = {1.0}, .den = {1.0, 1.0}};
        CHECK(design::pid_pole_placement<2>(plant, {0.95, 0.96}, 0.01).status == design::DesignStatus::invalid_parameter);
        CHECK(design::pid_pole_placement<3>(plant, {0.8, 0.7, 0.6}, 0.01).status == design::DesignStatus::invalid_parameter);
        CHECK(design::pid_pole_placement<2>(plant, {0.0, 0.0}, 0.01).status == design::DesignStatus::invalid_parameter);
        const auto axis = design::pid_pole_placement<2>(plant, {complex<double>{0.0, 3.0}, complex<double>{0.0, -3.0}}, 0.01);
        CHECK(axis.status == design::DesignStatus::invalid_parameter);
    }

    TEST_CASE("discrete placement keeps the ZOH gain accurate in float for light friction") {
        // 1/(1e-4 s + b_visc) at Ts = 1e-4: x = b_visc Ts/1e-4 is tiny and 1 − e^(−x) cancels.
        const damp::array<damp::complex<float>, 2>  pf{damp::complex<float>{-1053.6f, 0.0f}, damp::complex<float>{-1053.6f, 0.0f}};
        const damp::array<damp::complex<double>, 2> pd{-1053.6, -1053.6};
        for (const double bv : {1e-5, 1e-6, 1e-7}) {
            const auto f = design::pid_pole_placement(TransferFunction<1, 2, float>{.num = {1.0f}, .den = {1e-4f, static_cast<float>(bv)}}, pf, 1e-4f);
            const auto d = design::pid_pole_placement(TransferFunction<1, 2>{.num = {1.0}, .den = {1e-4, bv}}, pd, 1e-4);
            REQUIRE(f.success);
            CHECK(static_cast<double>(f.Kp) == doctest::Approx(d.Kp).epsilon(1e-5));
            CHECK(static_cast<double>(f.Ki) == doctest::Approx(d.Ki).epsilon(1e-4));
        }
    }

    TEST_CASE("Pole placement K=0 or Ts=0 returns finite zeros") {
        const auto k0 = design::pid_pole_placement(TransferFunction<1, 2>{.num = {0.0}, .den = {1.0, 1.0}}, damp::array<damp::complex<double>, 2>{-5.0, -10.0}, 0.1);
        CHECK(std::isfinite(k0.Kp));
        CHECK(std::isfinite(k0.Ki));
        CHECK(k0.Kp == doctest::Approx(0.0));
        CHECK(k0.Ki == doctest::Approx(0.0));

        const auto ts0 = design::pid_pole_placement(TransferFunction<1, 2>{.num = {1.0}, .den = {1.0, 1.0}}, damp::array<damp::complex<double>, 2>{-5.0, -10.0}, 0.0);
        CHECK(std::isfinite(ts0.Kp));
        CHECK(ts0.Kp == doctest::Approx(0.0));
        CHECK(ts0.Ki == doctest::Approx(0.0));
    }
}

TEST_SUITE("PID Design - Back-calculation seed") {
    // Plant 1/(s + 1) at Ts = 0.1: a = e^(−0.1), b = 1 − a. Poles 0.95 and 0.96 have
    // p₁p₂ > a, so Kp = (a − p₁p₂)/b < 0.
    constexpr double kTs = 0.1;
    constexpr auto   kPlant = TransferFunction<1, 2>{.num = {1.0}, .den = {1.0, 1.0}};

    TEST_CASE("a placement with Kp < 0 seeds Kbc > 0, the sign that opposes windup") {
        constexpr auto g = design::pid_pole_placement(kPlant, damp::array<damp::complex<double>, 2>{s_pole(0.95, kTs), s_pole(0.96, kTs)}, kTs);
        static_assert(g.success);
        static_assert(g.Kp < 0.0 && g.Ki > 0.0);
        CHECK(g.Kbc > 0.0);
        CHECK(g.Kbc == doctest::Approx(-g.Kp)); // |Kp| ≥ Ki·Ts here: T_t = |T_i|
        CHECK(g.Kbc / g.Ki > 0.0);              // tracking time T_t = Kbc/Ki
        static_assert(g.discretize(kTs).success);
    }

    TEST_CASE("the Kp < 0 placement's integrator does not wind up under a clamp") {
        auto gains = design::pid_pole_placement(kPlant, damp::array<damp::complex<double>, 2>{s_pole(0.95, kTs), s_pole(0.96, kTs)}, kTs);
        gains.u_lim = {.lower = -1.0, .upper = 1.0};
        PIController<double> law{gains.discretize(kTs)};
        PIController<double> clamp_only{[&] {
            auto g = gains;
            g.Kbc = 0.0; // output clamp only: the integrator runs away
            return g.discretize(kTs);
        }()};

        // r = 10 needs u = 10 on a unit-gain plant; the clamp holds u at 1 and y settles at 1.
        const double a = std::exp(-kTs);
        const double b = 1.0 - a;
        double       y = 0.0;
        double       y_open = 0.0;
        for (int k = 0; k < 400; ++k) {
            y = (a * y) + (b * law.control(10.0, y));
            y_open = (a * y_open) + (b * clamp_only.control(10.0, y_open));
        }
        // Saturated equilibrium: ΔI = Kᵢ·Ts·e + Kᵢ·Ts·(u_sat − u_unsat)/Kbc = 0 gives
        // u_unsat − u_sat = Kbc·e; the stored I (after the tick's correction) is then
        // u_sat + (Kbc − Kp − Kᵢ·Ts)·e.
        const double e = 10.0 - y;
        CHECK(y == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(law.integral == doctest::Approx(1.0 + ((gains.Kbc - gains.Kp - (gains.Ki * kTs)) * e)).epsilon(1e-6));
        CHECK(clamp_only.integral > 10.0 * law.integral);

        // Back to an attainable reference: the anti-windup law leaves saturation within a
        // few T_t; the clamp-only law is still pinned at the limit.
        int  ticks = 0;
        bool left = false;
        for (; ticks < 50 && !left; ++ticks) {
            const double u = law.control(0.5, y);
            left = u < 1.0;
            y = (a * y) + (b * u);
            y_open = (a * y_open) + (b * clamp_only.control(0.5, y_open));
        }
        CHECK(left);
        CHECK(ticks < 20);
        CHECK(clamp_only.control(0.5, y_open) == doctest::Approx(1.0));
    }

    TEST_CASE("an I-only placement (Kp = 0) seeds Kbc = Ki·Ts, a one-sample tracking time") {
        // s = −½ twice: z₁z₂ = e^(−kTs) = a ⇒ Kp = 0
        const auto g = design::pid_pole_placement(kPlant, damp::array<damp::complex<double>, 2>{-0.5, -0.5}, kTs);
        REQUIRE(g.success);
        CHECK(g.Kp == doctest::Approx(0.0));
        CHECK(g.Ki > 0.0);
        CHECK(g.Kbc == doctest::Approx(g.Ki * kTs));
        CHECK(g.Kbc / g.Ki == doctest::Approx(kTs));
    }

    TEST_CASE("a negative-gain plant seeds Kbc with Ki's sign") {
        constexpr auto g = design::pid_pole_placement(TransferFunction<1, 2>{.num = {-1.0}, .den = {1.0, 1.0}}, damp::array<damp::complex<double>, 2>{s_pole(0.5, kTs), s_pole(0.3, kTs)}, kTs);
        static_assert(g.Ki < 0.0);
        CHECK(g.Kbc < 0.0);
        CHECK(g.Kbc == doctest::Approx(g.Kp)); // Kp < 0 too: Kbc = −|Kp|
        static_assert(g.discretize(kTs).success);
    }

    TEST_CASE("pidtune seeds Kbc for a tracking time of 1/ωc") {
        constexpr StateSpace<1, 1, 1> sys{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}};
        for (const double wc : {1.5, 2.0, 8.0}) { // below ωc = tan 30° the plant lags too little: C would need Kp < 0
            for (const auto mode : {PIDMode::PI, PIDMode::PID}) {
                const auto g = design::pidtune(sys, {.wc = wc, .mode = mode});
                REQUIRE(g.success);
                REQUIRE(g.Ki > 0.0);
                CHECK(g.Kbc == doctest::Approx(g.Ki / wc));
                CHECK(g.Kbc / g.Ki == doctest::Approx(1.0 / wc));
            }
        }
    }

    TEST_CASE("the continuous placements keep Kbc = |Kp|") {
        constexpr auto pi = design::pid_pole_placement(TransferFunction<1, 2>{.num = {1.0}, .den = {1.2e-3, 4.0e-4}}, {.omega = 1000.0});
        CHECK(pi.Kbc == doctest::Approx(pi.Kp));
        constexpr auto slow = design::current_loop_pi(200e-6, 0.5, {.omega = 10.0}); // R/L above ω: Kp = ωL
        CHECK(slow.Kbc == doctest::Approx(10.0 * 200e-6));
        CHECK(slow.Kbc / slow.Ki == doctest::Approx(200e-6 / 0.5)); // T_t = T_i = L/R
    }
}

TEST_SUITE("PID Design - Double-integrator pole placement") {
    TEST_CASE("PD places a critically-damped pair on 1/(J s^2)") {
        constexpr double J = 0.05;
        constexpr double w = 8.0;
        constexpr auto   plant = TransferFunction<1, 3>{.num = {1.0}, .den = {J, 0.0, 0.0}};
        constexpr auto   r = design::pid_pole_placement(plant, {.omega = w});
        static_assert(r.Ki == 0.0);
        // J Tf s³ + J s² + (Kp Tf + Kd) s + Kp = J Tf (s + ω)²(s + 10ω)
        CHECK(r.Tf == doctest::Approx(1.0 / (12.0 * w)));
        CHECK(r.Kp == doctest::Approx(J * r.Tf * 10.0 * w * w * w));
        CHECK(r.Ki == doctest::Approx(0.0));
        CHECK(r.c == doctest::Approx(0.0));
        CHECK(r.Kbc == doctest::Approx(r.Kp));
        check_placed(plant, r, w, 1.0, 0.0);
    }

    TEST_CASE("PID places (s^2 + 2ζω s + ω^2)(s + ωi)") {
        constexpr double J = 0.05;
        constexpr double w = 8.0;
        constexpr double z = 1.0;
        constexpr double wi = w / 8.0;
        constexpr auto   plant = TransferFunction<1, 3>{.num = {1.0}, .den = {J, 0.0, 0.0}};
        constexpr auto   r = design::pid_pole_placement(plant, {.omega = w, .zeta = z, .omega_i = wi});
        static_assert(r.success);
        // The delivered loop, filter included: the requested poles and the filter pole at −10ω.
        check_placed(plant, r, w, z, wi);
        const auto under = design::pid_pole_placement(plant, {.omega = w, .zeta = 0.7, .omega_i = wi});
        check_placed(plant, under, w, 0.7, wi);
    }

    TEST_CASE("bad J / bandwidth / zeta / omega_i are refused with zero gains") {
        CHECK(design::pid_pole_placement(TransferFunction<1, 3>{.num = {1.0}, .den = {0.0, 0.0, 0.0}}, {.omega = 8.0}).Kp == 0.0);
        CHECK(design::pid_pole_placement(TransferFunction<1, 3>{.num = {1.0}, .den = {0.05, 0.0, 0.0}}, {.omega = -1.0}).Kp == 0.0);
        CHECK(design::pid_pole_placement(TransferFunction<1, 3>{.num = {1.0}, .den = {0.05, 0.0, 0.0}}, {.omega = 8.0, .zeta = 0.0}).Kp == 0.0);
        CHECK(design::pid_pole_placement(TransferFunction<1, 3>{.num = {1.0}, .den = {0.05, 0.0, 0.0}}, {.omega = 8.0, .zeta = 1.0, .omega_i = -0.1}).Kp == 0.0);
    }
}

TEST_SUITE("PID Design - Type Conversion") {
    TEST_CASE("Design result converts to PIDController via discretize") {
        constexpr double Ts = 0.01;
        constexpr auto   result = design::ziegler_nichols(10.0, 2.0);
        PIDController    controller(result.discretize(Ts).as<float>());
        CHECK(controller.Kp == doctest::Approx(6.0f).epsilon(1e-4));
        // Discrete Ki = continuous Ki * Ts
        CHECK(controller.Ki == doctest::Approx(static_cast<float>(6.0 * Ts)).epsilon(1e-4));
        CHECK(controller.Kd == doctest::Approx(1.5f).epsilon(1e-4));
    }

    TEST_CASE("Runtime design feeds into PIDController via discretize") {
        constexpr double Ts = 0.01;
        auto             result = design::simc(1.0, 0.5, 2.0, 1.0);
        PIDController    controller(result.discretize(Ts));
        double           u = controller.control(1.0, 0.0);
        CHECK(u != 0.0); // Non-zero output for non-zero error
    }

    TEST_CASE("PIDGains::to_tf yields Kp + Ki/s + Kd*s/(1+Tf*s)") {
        // Ideal PID (Tf = 0): C(s) = (Kd s^2 + Kp s + Ki)/s, highest power first.
        constexpr auto tf = design::pid(2.0, 3.0, 0.5).to_tf();
        CHECK(tf.num[0] == doctest::Approx(0.5)); // Kd
        CHECK(tf.num[1] == doctest::Approx(2.0)); // Kp
        CHECK(tf.num[2] == doctest::Approx(3.0)); // Ki
        CHECK(tf.den[0] == doctest::Approx(0.0)); // Tf
        CHECK(tf.den[1] == doctest::Approx(1.0));
        CHECK(tf.den[2] == doctest::Approx(0.0));

        // Filtered (Tf > 0): num = {Kp*Tf + Kd, Kp + Ki*Tf, Ki}, den = {Tf, 1, 0}.
        constexpr auto tff = design::PIDGains{.Kp = 2.0, .Ki = 3.0, .Kd = 0.5, .Tf = 0.1}.to_tf();
        CHECK(tff.den[0] == doctest::Approx(0.1));
        CHECK(tff.num[1] == doctest::Approx(2.0 + (3.0 * 0.1)));
        CHECK(tff.num[0] == doctest::Approx((2.0 * 0.1) + 0.5));
    }
}

TEST_SUITE("PID Design - Margin helpers from time-domain targets") {
    TEST_CASE("Damping ratio from overshoot percent") {
        constexpr double zeta_10 = design::damping_ratio_from_overshoot_percent(10.0);
        CHECK(zeta_10 == doctest::Approx(0.591).epsilon(1e-3));

        constexpr double zeta_0 = design::damping_ratio_from_overshoot_percent(0.0);
        CHECK(zeta_0 == doctest::Approx(1.0).epsilon(1e-12));
    }

    TEST_CASE("Phase margin estimate from damping ratio") {
        constexpr double pm = design::phase_margin_from_damping_ratio(0.7);
        CHECK(pm > 60.0);
        CHECK(pm < 70.0);
    }
}

TEST_SUITE("PID Design - PI Pole Placement (first-order plant)") {

    // Plant 1/(a1 s + a0) under PI closes to a1 s² + (a0+Kp) s + Ki, i.e.
    //   s² + 2ζω s + ω²  with  ω = √(Ki/a1),  ζ = (a0+Kp)/(2 a1 ω).
    TEST_CASE("places the closed-loop pole pair at omega_bw") {
        constexpr double a1 = 1.2e-3; // e.g. inductance L or inertia J
        constexpr double a0 = 4.0e-4; // e.g. resistance R or viscous b
        const double     omega = 2.0 * std::numbers::pi * 200.0;

        const auto g = design::pid_pole_placement(TransferFunction<1, 2>{.num = {1.0}, .den = {a1, a0}}, {.omega = omega});

        CHECK(g.Kp == doctest::Approx((2.0 * omega * a1) - a0));
        CHECK(g.Ki == doctest::Approx(a1 * omega * omega));
        CHECK(g.Kbc == doctest::Approx(g.Kp)); // anti-windup tracking = Kp

        const double omega_n = std::sqrt(g.Ki / a1);
        const double zeta = (a0 + g.Kp) / (2.0 * a1 * omega_n);
        CHECK(omega_n == doctest::Approx(omega)); // bandwidth as requested
        CHECK(zeta == doctest::Approx(1.0));      // critically damped (default)
    }

    TEST_CASE("honors a non-unity damping ratio") {
        constexpr double a1 = 5.0e-4, a0 = 0.0, omega = 300.0, zeta = 0.7;
        const auto       g = design::pid_pole_placement(TransferFunction<1, 2>{.num = {1.0}, .den = {a1, a0}}, {.omega = omega, .zeta = zeta});
        const double     omega_n = std::sqrt(g.Ki / a1);
        const double     z = (a0 + g.Kp) / (2.0 * a1 * omega_n);
        CHECK(omega_n == doctest::Approx(omega));
        CHECK(z == doctest::Approx(zeta));
    }

    // (FOController::tune's use of this kernel is checked in test_foc.cpp,
    //  where foc.hpp is in scope.)

    // TF overload: b0/(a1 s + a0) reduces to monic by dividing den by b0.
    TEST_CASE("transfer-function overload matches the scalar form") {
        constexpr double a1 = 1.0e-3, a0 = 0.8, omega = 1500.0;

        // Unit numerator → identical to the scalar (a1,a0) call.
        constexpr TransferFunction<1, 2, double> plant_unit{.num = {1.0}, .den = {a1, a0}};
        const auto                               from_tf = design::pid_pole_placement(plant_unit, {.omega = omega});
        const auto                               from_scalar = design::pid_pole_placement(TransferFunction<1, 2>{.num = {1.0}, .den = {a1, a0}}, {.omega = omega});
        CHECK(from_tf.Kp == doctest::Approx(from_scalar.Kp));
        CHECK(from_tf.Ki == doctest::Approx(from_scalar.Ki));

        // Numerator gain b0 ≠ 1 must divide through: b0/(a1 s+a0) ≡ 1/((a1/b0)s+(a0/b0)).
        constexpr double                         b0 = 2.5;
        constexpr TransferFunction<1, 2, double> plant_gain{.num = {b0}, .den = {a1, a0}};
        const auto                               from_gain = design::pid_pole_placement(plant_gain, {.omega = omega});
        const auto                               expected = design::pid_pole_placement(TransferFunction<1, 2>{.num = {1.0}, .den = {a1 / b0, a0 / b0}}, {.omega = omega});
        CHECK(from_gain.Kp == doctest::Approx(expected.Kp));
        CHECK(from_gain.Ki == doctest::Approx(expected.Ki));
    }

    // Degenerate plant (zero input gain, e.g. a default-constructed TF) must not
    // divide by zero — it returns an inert all-zero PI.
    TEST_CASE("transfer-function overload guards a zero numerator gain") {
        constexpr TransferFunction<1, 2, double> degenerate{}; // num={0}, den={0,0}
        const auto                               g = design::pid_pole_placement(degenerate, {.omega = 1000.0});
        CHECK(g.Kp == doctest::Approx(0.0));
        CHECK(g.Ki == doctest::Approx(0.0));
    }

    TEST_CASE("pid_pole_placement refuses a pair the plant pole already outruns (Kp ≤ 0)") {
        // 1/(L s + R): Kp = 2ζωL − R ≤ 0 when the requested pair is slower than R/L.
        constexpr double L = 1.0e-3, R = 2.0;
        const auto       slow = design::pid_pole_placement(TransferFunction<1, 2>{.num = {1.0}, .den = {L, R}}, {.omega = 500.0});
        CHECK_FALSE(slow.success);
        CHECK(slow.status == design::DesignStatus::not_placeable);
        CHECK(slow.Kp == 0.0);
        CHECK(design::pid_pole_placement(TransferFunction<1, 2>{.num = {1.0}, .den = {L, R}}, {.omega = 1001.0}).success);

        // A PI places two poles: an integrator pole on a first-order plant is refused.
        const auto third = design::pid_pole_placement(TransferFunction<1, 2>{.num = {1.0}, .den = {L, R}}, {.omega = 5000.0, .omega_i = 100.0});
        CHECK(third.status == design::DesignStatus::invalid_parameter);
    }

    TEST_CASE("current_loop_pi cancels an R/L pole faster than ω, and is continuous at ω = R/L") {
        // 1/(L s + R), R/L = 2000 rad/s. Below it: Kp = ωL, Ki = ωR (the PI zero at −R/L).
        constexpr double L = 1.0e-3, R = 2.0;
        constexpr auto   slow = design::current_loop_pi(L, R, {.omega = 500.0});
        static_assert(slow.success);
        CHECK(slow.Kp == doctest::Approx(500.0 * L));
        CHECK(slow.Ki == doctest::Approx(500.0 * R));

        // Above it, the critically damped pair; at ω = R/L both rules give Kp = R, Ki = R²/L.
        const auto below = design::current_loop_pi(L, R, {.omega = 1999.0});
        const auto above = design::current_loop_pi(L, R, {.omega = 2001.0});
        CHECK(below.Kp == doctest::Approx(R).epsilon(1e-3));
        CHECK(above.Kp == doctest::Approx(R).epsilon(1e-3));
        CHECK(below.Ki == doctest::Approx(R * R / L).epsilon(1e-3));
        CHECK(above.Ki == doctest::Approx(R * R / L).epsilon(1e-3));
        CHECK(above.Kp == doctest::Approx((2.0 * 2001.0 * L) - R));

        // A pair with Kp ≤ 0 above the pole (ζ ≤ ½) cancels it too; a bad input is still refused.
        CHECK(design::current_loop_pi(L, R, {.omega = 3000.0, .zeta = 0.3}).Kp == doctest::Approx(3000.0 * L));
        CHECK(design::current_loop_pi(L, R, {.omega = -1.0}).status == design::DesignStatus::invalid_parameter);
    }

    TEST_CASE("second-order plant with damping: Kd and Kp net out the plant's own terms") {
        // 1/(J s² + b s): PD closes to J s² + (b + Kd) s + Kp.
        constexpr double J = 0.05, b = 0.1, w = 8.0;
        constexpr auto   plant = TransferFunction<1, 3>{.num = {1.0}, .den = {J, b, 0.0}};
        constexpr auto   r = design::pid_pole_placement(plant, {.omega = w});
        static_assert(r.success);
        check_placed(plant, r, w, 1.0, 0.0);

        // More plant damping than the pair asks for needs Kd < 0: refused.
        const auto over = design::pid_pole_placement(TransferFunction<1, 3>{.num = {1.0}, .den = {J, 10.0, 0.0}}, {.omega = w});
        CHECK(over.status == design::DesignStatus::not_placeable);
    }

    TEST_CASE("discrete placement on an integrator plant uses b = b₀Ts/a₁") {
        // 1/s sampled at Ts: G(z) = Ts/(z − 1). PI poles at the images z = {0.5, 0.6}.
        constexpr double Ts = 0.01;
        constexpr auto   r = design::pid_pole_placement(TransferFunction<1, 2>{.num = {1.0}, .den = {1.0, 0.0}}, damp::array<damp::complex<double>, 2>{s_pole(0.5, Ts), s_pole(0.6, Ts)}, Ts);
        static_assert(r.success);
        // (z − 1)² + Ts (Kp (z − 1) + Ki Ts z) = z² − 1.1 z + 0.3
        CHECK(1.0 - (Ts * r.Kp) == doctest::Approx(0.3));
        CHECK((Ts * r.Kp) + (Ts * Ts * r.Ki) - 2.0 == doctest::Approx(-1.1));
    }

    TEST_CASE("state-space overloads give the transfer-function overloads' gains") {
        // First order: 2.5/(1e-3 s + 0.8) as ẋ = −800 x + u, y = 2500 x.
        constexpr double                   b0 = 2.5, a1 = 1.0e-3, a0 = 0.8, omega = 1500.0;
        constexpr StateSpace<1, 1, 1>      first{.A = {{-a0 / a1}}, .B = {{1.0}}, .C = {{b0 / a1}}};
        constexpr TransferFunction<1, 2>   first_tf{.num = {b0}, .den = {a1, a0}};
        constexpr TransferFunction<1, 2>   first_monic{.num = {b0 / a1}, .den = {1.0, a0 / a1}}; // ss2tf of `first`
        constexpr design::PIDGains<double> pi_ss = design::pid_pole_placement(first, {.omega = omega});
        static_assert(pi_ss.success);
        const auto pi_monic = design::pid_pole_placement(first_monic, {.omega = omega});
        CHECK(pi_ss.Kp == pi_monic.Kp); // the SS overload forwards ss2tf's TF unchanged
        CHECK(pi_ss.Ki == pi_monic.Ki);
        const auto pi_tf = design::pid_pole_placement(first_tf, {.omega = omega});
        CHECK(pi_ss.Kp == doctest::Approx(pi_tf.Kp).epsilon(1e-12));
        CHECK(pi_ss.Ki == doctest::Approx(pi_tf.Ki).epsilon(1e-12));
        CHECK(pi_ss.Kbc == doctest::Approx(pi_tf.Kbc).epsilon(1e-12));

        // Discrete counterpart on the same continuous plant.
        constexpr double Ts = 1.0e-4;
        const auto       dpi_ss = design::pid_pole_placement(first, damp::array<damp::complex<double>, 2>{-1000.0, -2000.0}, Ts);
        const auto       dpi_tf = design::pid_pole_placement(first_tf, damp::array<damp::complex<double>, 2>{-1000.0, -2000.0}, Ts);
        REQUIRE(dpi_ss.success);
        CHECK(dpi_ss.Kp == doctest::Approx(dpi_tf.Kp).epsilon(1e-12));
        CHECK(dpi_ss.Ki == doctest::Approx(dpi_tf.Ki).epsilon(1e-12));
        const auto dpid_ss = design::pid_pole_placement(first, damp::array<damp::complex<double>, 3>{-1000.0, -2000.0, -3000.0}, Ts);
        const auto dpid_tf = design::pid_pole_placement(first_tf, damp::array<damp::complex<double>, 3>{-1000.0, -2000.0, -3000.0}, Ts);
        REQUIRE(dpid_ss.success);
        CHECK(dpid_ss.Kd == doctest::Approx(dpid_tf.Kd).epsilon(1e-12));

        // Second order: mass-spring-damper 1/(m s² + b s + k) with states (x, v).
        constexpr double                 m = 0.05, b = 0.1, k = 2.0, w = 40.0;
        constexpr StateSpace<2, 1, 1>    msd{.A = {{0.0, 1.0}, {-k / m, -b / m}}, .B = {{0.0}, {1.0 / m}}, .C = {{1.0, 0.0}}};
        constexpr TransferFunction<1, 3> msd_tf{.num = {1.0}, .den = {m, b, k}};
        const auto                       pid_ss = design::pid_pole_placement(msd, {.omega = w, .omega_i = w / 8.0});
        const auto                       pid_tf = design::pid_pole_placement(msd_tf, {.omega = w, .omega_i = w / 8.0});
        REQUIRE(pid_ss.success);
        CHECK(pid_ss.Kp == doctest::Approx(pid_tf.Kp).epsilon(1e-12));
        CHECK(pid_ss.Ki == doctest::Approx(pid_tf.Ki).epsilon(1e-12));
        CHECK(pid_ss.Kd == doctest::Approx(pid_tf.Kd).epsilon(1e-12));
        CHECK(pid_ss.Tf == doctest::Approx(pid_tf.Tf).epsilon(1e-12));

        // The same plant in mixed coordinates is still all-pole: C·B cancels to rounding residue.
        const auto mixed = similarity_transform(msd, Matrix<2, 2>{{1000.0, 0.2}, {0.3, 1.7}});
        REQUIRE(mixed.success);
        const auto pid_mixed = design::pid_pole_placement(mixed.sys, {.omega = w, .omega_i = w / 8.0});
        REQUIRE(pid_mixed.success);
        CHECK(pid_mixed.Kp == doctest::Approx(pid_tf.Kp).epsilon(1e-9));
        CHECK(pid_mixed.Kd == doctest::Approx(pid_tf.Kd).epsilon(1e-9));
    }

    TEST_CASE("state-space overloads refuse a finite zero, D ≠ 0 and a discrete plant") {
        using design::DesignStatus;
        constexpr double m = 0.05, b = 0.1, k = 2.0;
        // C = [1, 0.01]: y = x + 0.01 v, so G(s) = (0.01 s + 1)/(m s² + b s + k) has a zero at −100.
        constexpr StateSpace<2, 1, 1> zero{.A = {{0.0, 1.0}, {-k / m, -b / m}}, .B = {{0.0}, {1.0 / m}}, .C = {{1.0, 0.01}}};
        const auto                    with_zero = design::pid_pole_placement(zero, {.omega = 40.0});
        CHECK_FALSE(with_zero.success);
        CHECK(with_zero.status == design::DesignStatus::invalid_parameter);
        CHECK(with_zero.Kp == 0.0);

        constexpr StateSpace<1, 1, 1> feedthrough{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}, .D = {{0.5}}};
        CHECK(design::pid_pole_placement(feedthrough, {.omega = 10.0}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::pid_pole_placement(feedthrough, damp::array<damp::complex<double>, 2>{-50.0, -90.0}, 0.01).status == design::DesignStatus::invalid_parameter);

        constexpr StateSpace<1, 1, 1> no_gain{.A = {{-1.0}}, .B = {{1.0}}, .C = {{0.0}}};
        CHECK(design::pid_pole_placement(no_gain, {.omega = 10.0}).status == design::DesignStatus::invalid_parameter);

        constexpr StateSpace<1, 1, 1> sampled{.A = {{0.9}}, .B = {{1.0}}, .C = {{1.0}}, .Ts = 0.01};
        CHECK(design::pid_pole_placement(sampled, {.omega = 10.0}).status == DesignStatus::invalid_sample_time);
        CHECK(design::pid_pole_placement(sampled, damp::array<damp::complex<double>, 2>{-50.0, -90.0}, 0.01).status == DesignStatus::invalid_sample_time);
    }

    TEST_CASE("amigo_kappa_tau matches FOPDT formula golden") {
        // Ks=1, Ku=2, Tu=1 → τ,L from the ultimate point, then the AMIGO PI rule
        // (Åström & Hägglund, Advanced PID Control, 2006, §7.6)
        constexpr double Ks = 1.0;
        constexpr double Ku = 2.0;
        constexpr double Tu = 1.0;
        const double     wu = 2.0 * damp::numbers::pi_v<double> / Tu;
        const double     prod = Ks * Ku;
        const double     tau = std::sqrt((prod * prod) - 1.0) / wu;
        const double     L = (damp::numbers::pi_v<double> - std::atan(wu * tau)) / wu;
        const double     Kp_exp = (0.15 / Ks) + ((0.35 - (L * tau / ((L + tau) * (L + tau)))) * tau / (Ks * L));
        const double     Ti_exp = (0.35 * L) + (13.0 * L * tau * tau / ((tau * tau) + (12.0 * L * tau) + (7.0 * L * L)));
        const double     Ki_exp = Kp_exp / Ti_exp;

        const auto pid = design::amigo_kappa_tau(Ks, Ku, Tu);
        CHECK(pid.Kp == doctest::Approx(Kp_exp).epsilon(1e-12));
        CHECK(pid.Ki == doctest::Approx(Ki_exp).epsilon(1e-12));
        CHECK(pid.Kd == doctest::Approx(0.0));

        // Non-physical / no real FOPDT → zero gains
        CHECK(design::amigo_kappa_tau(0.0, Ku, Tu).Kp == doctest::Approx(0.0));
        CHECK(design::amigo_kappa_tau(Ks, 1.0, Tu).Kp == doctest::Approx(0.0)); // Ks*Ku <= 1
    }
}

TEST_SUITE("PID Design - success / status") {
    TEST_CASE("tuning rules refuse non-physical inputs through status") {
        using design::DesignStatus;
        const auto zn = design::ziegler_nichols(10.0, 0.0, PIDMode::PI);
        CHECK_FALSE(zn.success);
        CHECK(zn.status == design::DesignStatus::invalid_parameter);
        CHECK(zn.Kp == doctest::Approx(0.0));
        CHECK(design::cohen_coon(0.0, 1.0, 4.0).status == design::DesignStatus::invalid_parameter);
        CHECK(design::ziegler_nichols(10.0, 2.0).success);

        // Non-finite inputs, Ku ≤ 0, tau ≤ 0 and L < 0 are refused, never returned as gains.
        const double nan = std::nan("");
        CHECK(design::ziegler_nichols(nan, 1.0).status == design::DesignStatus::invalid_parameter);
        CHECK(design::ziegler_nichols(0.0, 1.0).status == design::DesignStatus::invalid_parameter);
        CHECK(design::ziegler_nichols(0.0, 1.0, PIDMode::P).status == design::DesignStatus::invalid_parameter);
        CHECK(design::tyreus_luyben(0.0, 1.0).status == design::DesignStatus::invalid_parameter);
        CHECK(design::cohen_coon(nan, 1.0, 1.0).status == design::DesignStatus::invalid_parameter);
        CHECK(design::ziegler_nichols_step(1.0, nan, 1.0).status == design::DesignStatus::invalid_parameter);
        CHECK(design::amigo_kappa_tau(1.0, nan, 1.0).status == design::DesignStatus::invalid_parameter);
        CHECK(design::simc(1.0, 1.0, 0.0, 1.0, PIDMode::P).status == design::DesignStatus::invalid_parameter);
        CHECK(design::simc(1.0, 0.1, -1.0, 0.1, PIDMode::P).status == design::DesignStatus::invalid_parameter);
        CHECK(design::simc(1.0, -0.05, 1.0, 0.1, PIDMode::PID).status == design::DesignStatus::invalid_parameter);
        CHECK(design::simc(1.0, 0.0, 1.0, 0.1, PIDMode::PI).success); // L = 0 (delay-free plant) is valid
        CHECK(design::lambda_tuning(1.0, -0.5, 1.0, 1.0).status == design::DesignStatus::invalid_parameter);
        CHECK(design::lambda_tuning(1.0, 0.5, 1.0, 0.0).status == design::DesignStatus::invalid_parameter);
        CHECK(design::lambda_tuning(1.0, 0.5, 1.0, nan).status == design::DesignStatus::invalid_parameter);

        // discretize keeps a failure and refuses a non-positive sample time.
        CHECK(zn.discretize(0.01).status == design::DesignStatus::invalid_parameter);
        const auto bad_ts = design::pid(1.0, 2.0, 0.0).discretize(0.0);
        CHECK_FALSE(bad_ts.success);
        CHECK(bad_ts.status == DesignStatus::invalid_sample_time);
    }

    TEST_CASE("pidtune reports a refused crossover") {
        const StateSpace<1, 1, 1> sys{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}};
        const auto                r = design::pidtune(sys, design::PIDTuneSpec<double>{.wc = 0.0});
        CHECK_FALSE(r.success);
        CHECK(r.status == design::DesignStatus::invalid_parameter);
    }

    TEST_CASE("current_loop_pi returns the discrete gains with their status") {
        constexpr design::DiscretePIDGains<double> des = design::current_loop_pi(250e-6, 0.0, {.omega = 2e4}).discretize(1e-5);
        static_assert(des.success);
        static constinit PIController<float> pi = des.as<float>();
        CHECK(pi.Kp > 0.0f);
        CHECK(design::current_loop_pi(250e-6, 0.0, {.omega = 2e4}).discretize(0.0).status == design::DesignStatus::invalid_sample_time);
    }
}
