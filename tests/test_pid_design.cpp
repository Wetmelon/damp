// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/law.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/design/pid_design.hpp"
#include "damp/design/status.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matlab.hpp"
#include "damp/matrix/core.hpp"
#include "damp/systems/realization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using design::DesignStatus;

namespace {

/// The sample period of the designs whose rule is continuous: their gains do not depend on it.
constexpr double kTs = 0.01;

/// A quiet NaN from its bits, read through a volatile so the compiler cannot fold it.
double runtime_nan() {
    volatile std::uint64_t bits = 0x7ff8000000000000ull;
    return std::bit_cast<double>(static_cast<std::uint64_t>(bits));
}

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

/// LoopShape on @p sys at kTs: the result's gains.
template<typename Sys>
design::PIDResult<double> loop_shape(const Sys& sys, const design::LoopShape<double>& rule) {
    return design::pid(sys, kTs, rule);
}

/// The 1-DOF loop LoopShape fits on 1/(s + 1).
constexpr StateSpace kLag1{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}};

/// 1/(s + 1)³: ultimate frequency √3, ultimate gain 8.
constexpr StateSpace kLag3{.A = {{-1.0, 1.0, 0.0}, {0.0, -1.0, 1.0}, {0.0, 0.0, -1.0}}, .B = {{0.0}, {0.0}, {1.0}}, .C = {{1.0, 0.0, 0.0}}};

} // namespace

/**
 * @brief Tests for design::pid (pid_design.hpp): the rules on their models, the gains forms, the
 *        options, and the deployed result
 */

TEST_SUITE("PID Design - SIMC") {
    TEST_CASE("SIMC PI tuning") {
        constexpr auto result = design::pid(design::FOPDT{.K = 2.0, .L = 0.5, .tau = 3.0}, kTs, design::SIMC{.tau_c = 1.0, .mode = PIDMode::PI});
        static_assert(result.success);
        // Kp = tau / (K * (tau_c + L)) = 3.0 / (2.0 * 1.5) = 1.0
        CHECK(result.gains.Kp == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(result.gains.Ki > 0.0);
        CHECK(result.gains.Kd == doctest::Approx(0.0));
        CHECK(result.Ts == kTs);
    }

    TEST_CASE("SIMC is a PI by default, and a PID adds the dead-time derivative Kd = Kp·L/2") {
        constexpr design::FOPDT plant{.K = 2.0, .L = 0.5, .tau = 3.0};
        constexpr auto          pi = design::pid(plant, kTs, design::SIMC{.tau_c = 1.0});
        constexpr auto          pid = design::pid(plant, kTs, design::SIMC{.tau_c = 1.0, .mode = PIDMode::PID});
        CHECK(pi.gains.Kd == 0.0);
        CHECK(pi.gains.Ki == doctest::Approx(1.0 / 3.0)); // Ti = min(τ, 4(τc + L)) = 3
        CHECK(pid.gains.Kp == doctest::Approx(pi.gains.Kp));
        CHECK(pid.gains.Kd == doctest::Approx(pi.gains.Kp * 0.25));
    }

    TEST_CASE("SIMC refuses a non-physical model with zero gains") {
        const auto z = design::pid(design::FOPDT{.K = 0.0, .L = 0.5, .tau = 3.0}, kTs, design::SIMC{.tau_c = 1.0});
        CHECK_FALSE(z.success);
        CHECK(z.status == DesignStatus::invalid_parameter);
        CHECK(z.gains.Kp == 0.0);
        CHECK(z.gains.Ki == 0.0);
        CHECK(design::pid(design::FOPDT{.K = 1.0, .L = 1.0, .tau = 0.0}, kTs, design::SIMC{.tau_c = 1.0, .mode = PIDMode::P}).status == DesignStatus::invalid_parameter);
        CHECK(design::pid(design::FOPDT{.K = 1.0, .L = 0.1, .tau = -1.0}, kTs, design::SIMC{.tau_c = 0.1, .mode = PIDMode::P}).status == DesignStatus::invalid_parameter);
        CHECK(design::pid(design::FOPDT{.K = 1.0, .L = -0.05, .tau = 1.0}, kTs, design::SIMC{.tau_c = 0.1, .mode = PIDMode::PID}).status == DesignStatus::invalid_parameter);
        CHECK(design::pid(design::FOPDT{.K = 1.0, .L = 0.5, .tau = 1.0}, kTs, design::SIMC{.tau_c = -0.5}).status == DesignStatus::invalid_parameter); // τc + L = 0
        CHECK(design::pid(design::FOPDT{.K = 1.0, .L = 0.0, .tau = 1.0}, kTs, design::SIMC{.tau_c = 0.1}).success);                                    // L = 0 (delay-free plant) is valid
        CHECK(design::pid(design::FOPDT{.K = -1.0, .L = 0.5, .tau = 2.0}, kTs, design::SIMC{.tau_c = 1.0}).gains.Kp < 0.0);                            // reverse-acting
    }

    TEST_CASE("a SIMC design deploys at run time") {
        auto          result = design::pid(design::FOPDT{.K = 1.0, .L = 0.5, .tau = 2.0}, kTs, design::SIMC{.tau_c = 1.0});
        PIDController controller(result);
        double        u = controller.step(1.0, 0.0);
        CHECK(u != 0.0); // Non-zero output for non-zero error
    }
}

TEST_SUITE("PID Design - AMIGO") {
    TEST_CASE("AMIGO PI on an ultimate point matches Åström & Hägglund (2006) §7.6") {
        // Ks=1, Ku=2, Tu=1 → τ,L from the ultimate point, then the AMIGO PI rule
        constexpr double Ks = 1.0;
        constexpr double Ku = 2.0;
        constexpr double Tu = 1.0;
        const double     wu = 2.0 * damp::numbers::pi_v<double> / Tu;
        const double     T = std::sqrt((Ks * Ks * Ku * Ku) - 1.0) / wu;
        const double     L = (damp::numbers::pi_v<double> - std::atan(wu * T)) / wu;
        const double     K = (0.15 / Ks) + ((0.35 - (L * T / ((L + T) * (L + T)))) * T / (Ks * L));
        const double     Ti = (0.35 * L) + (13.0 * L * T * T / ((T * T) + (12.0 * L * T) + (7.0 * L * L)));
        const auto       pi = design::pid(design::Ultimate{.Ku = Ku, .Tu = Tu, .K = Ks}, kTs, design::AMIGO{});
        REQUIRE(pi.success);
        CHECK(pi.gains.Kp == doctest::Approx(K).epsilon(1e-12));
        CHECK(pi.gains.Ki == doctest::Approx(K / Ti).epsilon(1e-12));
        CHECK(pi.gains.Kd == 0.0);
    }

    TEST_CASE("AMIGO on an FOPDT gives the PI of the ultimate point with the same model") {
        // The FOPDT recovered from (Ks, Ku, Tu) has that ultimate point: both routes are one PI.
        constexpr double Ks = 1.5;
        constexpr double Ku = 3.0;
        constexpr double Tu = 2.0;
        const double     wu = 2.0 * std::numbers::pi / Tu;
        const double     tau = std::sqrt((Ks * Ks * Ku * Ku) - 1.0) / wu;
        const double     L = (std::numbers::pi - std::atan(wu * tau)) / wu;
        const auto       from_relay = design::pid(design::Ultimate{.Ku = Ku, .Tu = Tu, .K = Ks}, kTs, design::AMIGO{});
        const auto       from_model = design::pid(design::FOPDT{.K = Ks, .L = L, .tau = tau}, kTs, design::AMIGO{});
        REQUIRE(from_relay.success);
        REQUIRE(from_model.success);
        CHECK(from_model.gains.Kp == doctest::Approx(from_relay.gains.Kp).epsilon(1e-12));
        CHECK(from_model.gains.Ki == doctest::Approx(from_relay.gains.Ki).epsilon(1e-12));
        // The formula on the model directly.
        const double LpT = L + tau;
        const double Kp = (0.15 / Ks) + ((0.35 - (L * tau / (LpT * LpT))) * tau / (Ks * L));
        const double Ti = (0.35 * L) + (13.0 * L * tau * tau / ((tau * tau) + (12.0 * L * tau) + (7.0 * L * L)));
        CHECK(from_model.gains.Kp == doctest::Approx(Kp).epsilon(1e-12));
        CHECK(from_model.gains.Ki == doctest::Approx(Kp / Ti).epsilon(1e-12));
    }

    TEST_CASE("AMIGO on a reverse-acting FOPDT negates the gains") {
        const auto direct = design::pid(design::FOPDT{.K = 2.0, .L = 1.0, .tau = 5.0}, kTs, design::AMIGO{});
        const auto reverse = design::pid(design::FOPDT{.K = -2.0, .L = 1.0, .tau = 5.0}, kTs, design::AMIGO{});
        REQUIRE(direct.success);
        REQUIRE(reverse.success);
        CHECK(reverse.gains.Kp == doctest::Approx(-direct.gains.Kp).epsilon(1e-12));
        CHECK(reverse.gains.Ki == doctest::Approx(-direct.gains.Ki).epsilon(1e-12));
    }

    TEST_CASE("AMIGO refuses what has no AMIGO PI") {
        CHECK(design::pid(design::Ultimate{.Ku = 2.0, .Tu = 1.0}, kTs, design::AMIGO{}).status == DesignStatus::invalid_parameter);           // K not known
        CHECK(design::pid(design::Ultimate{.Ku = 1.0, .Tu = 1.0, .K = 1.0}, kTs, design::AMIGO{}).status == DesignStatus::invalid_parameter); // K·Ku ≤ 1
        CHECK(design::pid(design::Ultimate{.Ku = runtime_nan(), .Tu = 1.0, .K = 1.0}, kTs, design::AMIGO{}).status == DesignStatus::invalid_parameter);
        CHECK(design::pid(design::FOPDT{.K = 1.0, .L = 0.0, .tau = 1.0}, kTs, design::AMIGO{}).status == DesignStatus::invalid_parameter); // the rule divides by L
        CHECK(design::pid(design::FOPDT{.K = 0.0, .L = 1.0, .tau = 1.0}, kTs, design::AMIGO{}).gains.Kp == 0.0);
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

TEST_SUITE("PID Design - LoopShape") {
    TEST_CASE("PI at ωc on 1/(s+1)") {
        const auto pi = design::pid(kLag1, kTs, design::LoopShape{.wc = 2.0, .phase_margin = 60.0, .mode = PIDMode::PI});
        REQUIRE(pi.success);
        CHECK(pi.gains.Kp > 0.0);
        CHECK(pi.gains.Ki > 0.0);
        CHECK(pi.gains.Kd == doctest::Approx(0.0));
        CHECK(pi.gains.b == doctest::Approx(1.0));
        CHECK(pi.gains.c == doctest::Approx(1.0));
    }

    TEST_CASE("LoopShape takes the derivative filter and the 2-DOF weights from the options") {
        const auto r = design::pid(kLag1, kTs, design::LoopShape{.wc = 2.0, .phase_margin = 70.0, .mode = PIDMode::PI}, design::DerivativeFilter{.Tf = 1e-3}, design::SetpointWeights{.b = 0.0, .c = 0.0});
        REQUIRE(r.success);
        CHECK(r.gains.Tf == doctest::Approx(1e-3));
        CHECK(r.gains.b == doctest::Approx(0.0));
        CHECK(r.gains.c == doctest::Approx(0.0));
        CHECK(r.gains.Ki > 0.0);
    }

    TEST_CASE("1-DOF LoopShape is unchanged by SetpointWeights{1, 1}") {
        const design::LoopShape rule{.wc = 2.0, .phase_margin = 60.0, .mode = PIDMode::PI};
        const auto              a = design::pid(kLag1, kTs, rule);
        const auto              b = design::pid(kLag1, kTs, rule, design::SetpointWeights{.b = 1.0, .c = 1.0});
        REQUIRE(a.success);
        REQUIRE(b.success);
        CHECK(a.gains.Kp == doctest::Approx(b.gains.Kp));
        CHECK(a.gains.Ki == doctest::Approx(b.gains.Ki));
        CHECK(a.gains.Kd == doctest::Approx(b.gains.Kd));
    }

    TEST_CASE("I-P LoopShape matches 1-DOF T −3 dB") {
        constexpr StateSpace    sys = kLag1;
        constexpr double        wc = 2.0;
        const design::LoopShape rule{.wc = wc, .phase_margin = 60.0, .mode = PIDMode::PI};
        const auto              r1 = design::pid(sys, kTs, rule).gains;
        const auto              r0_res = design::pid(sys, kTs, rule, design::SetpointWeights{.b = 0.0, .c = 0.0});
        REQUIRE(r0_res.success);
        const auto r0 = r0_res.gains;
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

    TEST_CASE("PI-D LoopShape on LC does not sign-flip") {
        constexpr double L = 100e-6;
        constexpr double C = 100e-6;
        constexpr double R = 4.0;
        const StateSpace sys{
            .A = Matrix<2, 2>{{0.0, -1.0 / L}, {1.0 / C, -1.0 / (R * C)}},
            .B = Matrix<2, 1>{{1.0 / L}, {0.0}},
            .C = Matrix{{0.0, 1.0}},
            .D = Matrix{{0.0}},
        };
        constexpr double wc = 27650.638614364343;
        for (int k = 600; k <= 955; k += 5) {
            const auto tuned = design::pid(sys, 1e-5, design::LoopShape{.wc = wc, .phase_margin = static_cast<double>(k) / 10.0}, design::SetpointWeights{.b = 1.0, .c = 0.0});
            REQUIRE(tuned.success);
            CHECK(tuned.gains.Kp > 0.0);
        }
        // The direct fit needs Kp < 0 here: refused rather than returned.
        const auto past = design::pid(sys, 1e-5, design::LoopShape{.wc = wc, .phase_margin = 96.0}, design::SetpointWeights{.b = 1.0, .c = 0.0});
        CHECK_FALSE(past.success);
        CHECK(past.status == DesignStatus::not_placeable);
    }

    TEST_CASE("LoopShape refuses a law that cannot realize the crossover with positive gains") {
        // 1/(s(s+1)) at ωc = 1 needs phase lead for 60°: a PI would need Ki < 0.
        constexpr StateSpace type1{.A = {{0.0, 1.0}, {0.0, -1.0}}, .B = {{0.0}, {1.0}}, .C = {{1.0, 0.0}}};
        CHECK(loop_shape(type1, {.wc = 1.0, .mode = PIDMode::PI}).status == DesignStatus::not_placeable);
        // 1/(s+1)² at ωc = 2: the PI fit had Ki < 0 and an unstable closed loop.
        constexpr StateSpace lag2{.A = {{-1.0, 1.0}, {0.0, -1.0}}, .B = {{0.0}, {1.0}}, .C = {{1.0, 0.0}}};
        CHECK(loop_shape(lag2, {.wc = 2.0, .mode = PIDMode::PI}).status == DesignStatus::not_placeable);
        // 1/(s+1)³ at ωc = 3: every law needed Kp < 0.
        CHECK(loop_shape(kLag3, {.wc = 3.0, .mode = PIDMode::PI}).status == DesignStatus::not_placeable);
        CHECK(loop_shape(kLag3, {.wc = 3.0, .mode = PIDMode::PID}).status == DesignStatus::not_placeable);
        // 1/(s+1): a PD at ωc = 5 would need Kd < 0, a PID at ωc = 0.1 Kp < 0.
        CHECK(loop_shape(kLag1, {.wc = 5.0, .mode = PIDMode::PD}).status == DesignStatus::not_placeable);
        CHECK(loop_shape(kLag1, {.wc = 0.1, .mode = PIDMode::PID}).status == DesignStatus::not_placeable);
    }

    TEST_CASE("LoopShape refuses an unstable closed loop (P on 1/(s+1)³ above the ultimate frequency)") {
        // Kp = 1/|P(j3)| = 10^1.5 > 8, the ultimate gain of 1/(s+1)³.
        CHECK(loop_shape(kLag3, {.wc = 3.0, .mode = PIDMode::P}).status == DesignStatus::not_stabilizable);
        CHECK(loop_shape(kLag3, {.wc = 0.5, .mode = PIDMode::P}).success);
    }

    TEST_CASE("LoopShape's certificate does not depend on the time unit") {
        // 1/(s+1)³ with time scaled by 1000 has every closed-loop root scaled by 1000. The Routh
        // test compared the leading coefficient with c[n] ~ ωc³ and refused the fast copy.
        constexpr double     a = 1000.0;
        constexpr StateSpace fast{.A = {{-a, a, 0.0}, {0.0, -a, a}, {0.0, 0.0, -a}}, .B = {{0.0}, {0.0}, {a}}, .C = {{1.0, 0.0, 0.0}}};
        const auto           slow_g = design::pid(kLag3, kTs, design::LoopShape{.wc = 0.5}, design::DerivativeFilter{.Tf = 0.2});
        const auto           fast_g = design::pid(fast, kTs, design::LoopShape{.wc = 500.0}, design::DerivativeFilter{.Tf = 2e-4});
        REQUIRE(slow_g.success);
        REQUIRE(fast_g.success);
        CHECK(fast_g.gains.Kp == doctest::Approx(slow_g.gains.Kp).epsilon(1e-9));
        CHECK(fast_g.gains.Ki == doctest::Approx(a * slow_g.gains.Ki).epsilon(1e-9));
        CHECK(fast_g.gains.Kd == doctest::Approx(slow_g.gains.Kd / a).epsilon(1e-9));
        // Past the ultimate gain (Kp = (1 + ωc²)^{3/2} > 8) both copies are refused.
        CHECK(loop_shape(kLag3, {.wc = 1.75, .mode = PIDMode::P}).status == DesignStatus::not_stabilizable);
        CHECK(loop_shape(fast, {.wc = 1750.0, .mode = PIDMode::P}).status == DesignStatus::not_stabilizable);

        // A float PI at ωc = 1000 rad/s on 1/(s+1) was refused above ωc ≈ 200 rad/s.
        constexpr StateSpace lag1f{.A = {{-1.0f}}, .B = {{1.0f}}, .C = {{1.0f}}};
        const auto           f = design::pid(lag1f, 1e-4f, design::LoopShape{.wc = 1000.0f, .mode = PIDMode::PI});
        REQUIRE(f.success);
        CHECK(static_cast<double>(f.gains.Kp) == doctest::Approx(-0.5 + (std::sqrt(3.0) / 2.0 * 1000.0)).epsilon(1e-5));
        // An R–L current loop 1000/(s + 100) at 1 kHz, in float.
        constexpr StateSpace rl{.A = {{-100.0f}}, .B = {{1000.0f}}, .C = {{1.0f}}};
        CHECK(design::pid(rl, 1e-5f, design::LoopShape{.wc = 6283.0f, .mode = PIDMode::PI}).success);
    }

    TEST_CASE("LoopShape on a reverse-acting plant returns the negated gains") {
        // −P needs −C for the same loop C·P: Kp, Ki, Kd ≤ 0.
        constexpr StateSpace<1, 1, 1> rev1{.A = kLag1.A, .B = kLag1.B, .C = {{-1.0}}};
        constexpr StateSpace<2, 1, 1> osc{.A = {{0.0, 1.0}, {-1.0, -0.5}}, .B = {{0.0}, {1.0}}, .C = {{1.0, 0.0}}};
        constexpr StateSpace<2, 1, 1> rev_osc{.A = osc.A, .B = osc.B, .C = {{-1.0, 0.0}}};
        const auto                    check_negated = [](const auto& direct, const auto& reverse, const design::LoopShape<double>& rule, double b) {
            const auto g = design::pid(direct, kTs, rule, design::SetpointWeights{.b = b});
            const auto r = design::pid(reverse, kTs, rule, design::SetpointWeights{.b = b});
            REQUIRE(g.success);
            REQUIRE(r.success);
            CHECK(r.gains.Kp == doctest::Approx(-g.gains.Kp).epsilon(1e-12));
            CHECK(r.gains.Ki == doctest::Approx(-g.gains.Ki).epsilon(1e-12));
            CHECK(r.gains.Kd == doctest::Approx(-g.gains.Kd).epsilon(1e-12));
            CHECK(r.gains.Kp < 0.0);
            CHECK(PIDController<double>{r}.valid());
        };
        check_negated(kLag1, rev1, {.wc = 2.0, .mode = PIDMode::P}, 1.0);
        check_negated(kLag1, rev1, {.wc = 2.0, .mode = PIDMode::PI}, 1.0);
        check_negated(kLag1, rev1, {.wc = 2.0, .mode = PIDMode::PID}, 1.0);
        check_negated(kLag1, rev1, {.wc = 2.0, .mode = PIDMode::PI}, 0.0);  // the 2-DOF crossover search
        check_negated(osc, rev_osc, {.wc = 2.0, .mode = PIDMode::PD}, 1.0); // the plant lags 162° at ωc = 2: the PD leads
        check_negated(osc, rev_osc, {.wc = 2.0, .mode = PIDMode::PID}, 1.0);

        // A refusal mirrors too: 1/(s+1)³ at ωc = 3 is past the ultimate gain (P) and needs lead (PI).
        constexpr StateSpace<3, 1, 1> rev3{.A = kLag3.A, .B = kLag3.B, .C = {{-1.0, 0.0, 0.0}}};
        CHECK(loop_shape(rev3, {.wc = 3.0, .mode = PIDMode::P}).status == DesignStatus::not_stabilizable);
        CHECK(loop_shape(rev3, {.wc = 3.0, .mode = PIDMode::PI}).status == DesignStatus::not_placeable);
    }

    TEST_CASE("LoopShape P meets the crossover: |Kp P(jωc)| = 1") {
        const auto p = loop_shape(kLag1, {.wc = 2.0, .mode = PIDMode::P});
        REQUIRE(p.success);
        CHECK(p.gains.Kp == doctest::Approx(std::sqrt(5.0)));
    }

    TEST_CASE("LoopShape meets |L| = 1 and the phase margin with a derivative filter") {
        constexpr StateSpace type1{.A = {{0.0, 1.0}, {0.0, -1.0}}, .B = {{0.0}, {1.0}}, .C = {{1.0, 0.0}}};
        for (const double Tf : {0.0, 0.1, 0.3}) {
            for (const auto mode : {PIDMode::PD, PIDMode::PID}) {
                const auto res = design::pid(type1, kTs, design::LoopShape{.wc = 1.0, .mode = mode}, design::DerivativeFilter{.Tf = Tf});
                REQUIRE(res.success);
                const auto L = res.gains.eval_fb(1.0) * (*eval_frf(type1, complex<double>{0.0, 1.0}))(0, 0);
                CHECK(abs(L) == doctest::Approx(1.0).epsilon(1e-12));
                CHECK(arg(L) * 180.0 / std::numbers::pi == doctest::Approx(-120.0).epsilon(1e-12));
            }
            const auto res = design::pid(kLag3, kTs, design::LoopShape{.wc = 1.0}, design::DerivativeFilter{.Tf = Tf});
            REQUIRE(res.success);
            const auto& g = res.gains;
            CHECK(g.Ki == doctest::Approx(g.Kp * g.Kp / (4.0 * g.Kd))); // Ti = 4 Td
            const auto L = g.eval_fb(1.0) * (*eval_frf(kLag3, complex<double>{0.0, 1.0}))(0, 0);
            CHECK(abs(L) == doctest::Approx(1.0).epsilon(1e-12));
            CHECK(arg(L) * 180.0 / std::numbers::pi == doctest::Approx(-120.0).epsilon(1e-12));
        }
    }

    TEST_CASE("LoopShape refuses a non-finite or out-of-range spec") {
        const double nan = runtime_nan();
        CHECK(loop_shape(kLag1, {.wc = 1.0, .phase_margin = nan}).status == DesignStatus::invalid_parameter);
        CHECK(loop_shape(kLag1, {.wc = 1.0, .phase_margin = 0.0}).status == DesignStatus::invalid_parameter);
        CHECK(loop_shape(kLag1, {.wc = 1.0, .phase_margin = 180.0}).status == DesignStatus::invalid_parameter);
        CHECK(design::pid(kLag1, kTs, design::LoopShape{.wc = 1.0}, design::DerivativeFilter{.Tf = -0.1}).status == DesignStatus::invalid_parameter);
        CHECK(loop_shape(kLag1, {.wc = nan}).status == DesignStatus::invalid_parameter);
        CHECK(loop_shape(kLag1, {.wc = 0.0}).status == DesignStatus::invalid_parameter);
        constexpr StateSpace sampled{.A = {{0.9}}, .B = {{1.0}}, .C = {{1.0}}, .Ts = 0.01};
        CHECK(loop_shape(sampled, {.wc = 1.0}).status == DesignStatus::invalid_sample_time);
    }

    TEST_CASE("LoopShape and PolePlacement accept a plant with noise channels") {
        constexpr StateSpace noisy{
            .A = {{0.0, 1.0}, {0.0, -1.0}},
            .B = {{0.0}, {1.0}},
            .C = {{1.0, 0.0}},
            .G = Matrix<2, 2>::identity(),
            .H = Matrix<1, 1>::identity(),
        };
        constexpr StateSpace plain{.A = noisy.A, .B = noisy.B, .C = noisy.C};
        const auto           a = loop_shape(noisy, {.wc = 1.0});
        REQUIRE(a.success);
        CHECK(a.gains.Kp == doctest::Approx(loop_shape(plain, {.wc = 1.0}).gains.Kp).epsilon(1e-14));
        const auto p = design::pid(noisy, kTs, design::PolePlacement{.omega = 3.0});
        REQUIRE(p.success);
        CHECK(p.gains.Kp == doctest::Approx(design::pid(plain, kTs, design::PolePlacement{.omega = 3.0}).gains.Kp).epsilon(1e-14));
    }

    TEST_CASE("LoopShape on a transfer function fits its realization") {
        constexpr TransferFunction<1, 3> tf{.num = {1.0}, .den = {1.0, 0.5, 1.0}};
        constexpr StateSpace             ss{.A = {{0.0, 1.0}, {-1.0, -0.5}}, .B = {{0.0}, {1.0}}, .C = {{1.0, 0.0}}};
        const auto                       from_tf = loop_shape(tf, {.wc = 1.0});
        const auto                       from_ss = loop_shape(ss, {.wc = 1.0});
        REQUIRE(from_tf.success);
        CHECK(from_tf.gains.Kp == doctest::Approx(from_ss.gains.Kp).epsilon(1e-12));
        CHECK(from_tf.gains.Ki == doctest::Approx(from_ss.gains.Ki).epsilon(1e-12));
        CHECK(from_tf.gains.Kd == doctest::Approx(from_ss.gains.Kd).epsilon(1e-12));
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

    TEST_CASE("matlab::pidtune(sys, wc) is a 60° PID, continuous, and c2d deploys it") {
        constexpr StateSpace sys{
            .A = Matrix{{0.0, 1.0}, {-1.0, -0.5}},
            .B = Matrix{{0.0}, {1.0}},
            .C = Matrix{{1.0, 0.0}},
            .D = Matrix{{0.0}},
        };
        const auto a = matlab::pidtune(sys, 1.0);
        const auto b = design::pid(sys, kTs, design::LoopShape{.wc = 1.0, .phase_margin = 60.0, .mode = PIDMode::PID});
        REQUIRE(a.has_value());
        REQUIRE(b.success);
        CHECK(a->Kp == doctest::Approx(b.gains.Kp));
        CHECK(a->Ki == doctest::Approx(b.gains.Ki));
        CHECK(a->Kd == doctest::Approx(b.gains.Kd));
        CHECK(a->tracking_time == 0.0);
        const auto deployed = matlab::c2d(*a, kTs);
        REQUIRE(deployed.success);
        CHECK(deployed.gains.Kp == doctest::Approx(b.gains.Kp).epsilon(1e-14));
        CHECK(deployed.Ts == kTs);
        CHECK_FALSE(matlab::pidtune(kLag3, 3.0).has_value()); // the PID needs Kp < 0 at ωc = 3: refused
    }

    TEST_CASE("PIDGains::to_ss realizes Kp + Ki/s + Kd s/(1 + Tf s)") {
        constexpr auto pi = design::PIDGains{.Kp = 2.0, .Ki = 4.0, .Kd = 0.0}.to_ss();
        static_assert(pi.has_value());
        static_assert(pi->A(0, 0) == 0.0);
        static_assert(pi->B(0, 0) == 1.0);
        static_assert(pi->C(0, 0) == 4.0);
        static_assert(pi->D(0, 0) == 2.0);
        constexpr auto L = *series(*pi, StateSpace{
                                            .A = {{-1.0}},
                                            .B = {{1.0}},
                                            .C = {{1.0}},
                                        });
        CHECK(L.A.rows() == 3);

        // The filtered derivative used to be dropped: C(jω) must match to_tf at every ω.
        design::PIDGains<double> gains{.Kp = 2.0, .Ki = 4.0, .Kd = 0.5};
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
        CHECK_FALSE(design::PIDGains{.Kp = 2.0, .Ki = 4.0, .Kd = 0.5}.to_ss().has_value());
    }
}

TEST_SUITE("PID Design - LoopShape on a relay test's ultimate point") {
    TEST_CASE("LoopShape on Ultimate is the fit at that point of the plant it came from") {
        // 1/(s + 1)³ crosses −180° at ω = √3 with |P| = 1/8: Ku = 8, Tu = 2π/√3. The plant fit at
        // ωc = √3 and the ultimate point's are one law.
        const double wu = std::sqrt(3.0);
        const double Tu = 2.0 * std::numbers::pi / wu;
        for (const double pm : {30.0, 45.0, 60.0}) {
            for (const double Tf : {0.0, 0.05}) {
                const auto from_relay = design::pid(design::Ultimate{.Ku = 8.0, .Tu = Tu}, kTs, design::LoopShape{.phase_margin = pm}, design::DerivativeFilter{.Tf = Tf});
                const auto from_plant = design::pid(kLag3, kTs, design::LoopShape{.wc = wu, .phase_margin = pm}, design::DerivativeFilter{.Tf = Tf});
                REQUIRE(from_relay.success);
                REQUIRE(from_plant.success);
                CHECK(from_relay.gains.Kp == doctest::Approx(from_plant.gains.Kp).epsilon(1e-9));
                CHECK(from_relay.gains.Ki == doctest::Approx(from_plant.gains.Ki).epsilon(1e-9));
                CHECK(from_relay.gains.Kd == doctest::Approx(from_plant.gains.Kd).epsilon(1e-9));
                CHECK(from_relay.gains.Tf == Tf);
            }
        }
    }

    TEST_CASE("LoopShape on Ultimate is Åström and Hägglund's phase-margin PID: Kp = Ku·cos φm, Ti = 4 Td") {
        constexpr double Ku = 5.0;
        constexpr double Tu = 0.8;
        const double     wu = 2.0 * std::numbers::pi / Tu;
        const double     pm = 45.0 * std::numbers::pi / 180.0;
        const auto       res = design::pid(design::Ultimate{.Ku = Ku, .Tu = Tu}, kTs, design::LoopShape{.phase_margin = 45.0});
        REQUIRE(res.success);
        CHECK(res.gains.Kp == doctest::Approx(Ku * std::cos(pm)).epsilon(1e-12));
        // ω Td − 1/(ω Ti) = tan φm with Ti = 4 Td
        const double Td = res.gains.Kd / res.gains.Kp;
        const double Ti = res.gains.Kp / res.gains.Ki;
        CHECK(Ti == doctest::Approx(4.0 * Td).epsilon(1e-12));
        CHECK((wu * Td) - (1.0 / (wu * Ti)) == doctest::Approx(std::tan(pm)).epsilon(1e-12));
    }

    TEST_CASE("LoopShape on Ultimate refuses what the point cannot give") {
        const design::Ultimate relay{.Ku = 5.0, .Tu = 0.8};
        CHECK(design::pid(relay, kTs, design::LoopShape{.wc = 1.0}).status == DesignStatus::invalid_parameter); // the crossover is the point's
        CHECK(design::pid(relay, kTs, design::LoopShape{.mode = PIDMode::PI}).status == DesignStatus::not_placeable);
        CHECK(design::pid(relay, kTs, design::LoopShape{.mode = PIDMode::P}).status == DesignStatus::not_stabilizable);
        CHECK(design::pid(relay, kTs, design::LoopShape{.mode = PIDMode::PD}).success);
        CHECK(design::pid(design::Ultimate{.Ku = 0.0, .Tu = 0.8}, kTs, design::LoopShape{}).status == DesignStatus::invalid_parameter);
        CHECK(design::pid(design::Ultimate{.Ku = 5.0, .Tu = runtime_nan()}, kTs, design::LoopShape{}).status == DesignStatus::invalid_parameter);
        // The weights are the law's, and do not move the crossover: one point has no tracking map.
        const auto ip = design::pid(relay, kTs, design::LoopShape{}, design::SetpointWeights{.b = 0.0, .c = 0.0});
        const auto one = design::pid(relay, kTs, design::LoopShape{});
        REQUIRE(ip.success);
        CHECK(ip.gains.Kp == doctest::Approx(one.gains.Kp).epsilon(1e-14));
        CHECK(ip.gains.b == 0.0);
    }
}

TEST_SUITE("PID Design - Poles sampled at Ts") {
    TEST_CASE("PI pole placement for first-order plant") {
        // Plant: G(s) = 1/(s+1), K=1, tau=1
        // Desired s-plane poles whose images are z=0.5 and z=0.3 (well inside unit circle)
        constexpr double K = 1.0;
        constexpr double tau = 1.0;
        constexpr double Ts = 0.1;
        constexpr double p1 = 0.5;
        constexpr double p2 = 0.3;

        constexpr auto res = design::pid(TransferFunction<1, 2>{.num = {K}, .den = {tau, 1.0}}, Ts, design::Poles{s_pole(p1, Ts), s_pole(p2, Ts)});
        static_assert(res.success);
        constexpr auto result = res.gains;
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

        constexpr auto result = design::pid(TransferFunction<1, 2>{.num = {K}, .den = {tau, 1.0}}, Ts, design::Poles{s_pole(p1, Ts), s_pole(p2, Ts), s_pole(p3, Ts)}).gains;
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

    TEST_CASE("Poles are s-plane poles, as place(sys, poles, Ts) takes them") {
        // s = −5 twice at Ts = 0.01 is z = e^(−0.05): a stable pair. Read as z-plane poles,
        // {−5, −5} gave an unstable loop with success.
        constexpr double Ts = 0.01;
        const auto       res = design::pid(TransferFunction<1, 2>{.num = {1.0}, .den = {1.0, 1.0}}, Ts, design::Poles{-5.0, -5.0});
        REQUIRE(res.success);
        const auto&  g = res.gains;
        const double a = std::exp(-Ts);
        const double b = 1.0 - a;
        const double z = std::exp(-5.0 * Ts);
        CHECK(a - (b * g.Kp) == doctest::Approx(z * z).epsilon(1e-12));
        CHECK((1.0 + a) - (b * g.Kp) - (b * g.Ki * Ts) == doctest::Approx(2.0 * z).epsilon(1e-12));
    }

    TEST_CASE("Poles places a complex pair and refuses an unpaired complex pole") {
        constexpr double Ts = 0.01;
        constexpr auto   plant = TransferFunction<1, 2>{.num = {1.0}, .den = {1.0, 1.0}};
        const auto       res = design::pid(plant, Ts, design::Poles{complex<double>{-5.0, 3.0}, complex<double>{-5.0, -3.0}});
        REQUIRE(res.success);
        const auto& g = res.gains;
        // z₁z₂ = e^(−10 Ts), z₁ + z₂ = 2 e^(−5 Ts) cos(3 Ts)
        const double a = std::exp(-Ts);
        const double b = 1.0 - a;
        CHECK(a - (b * g.Kp) == doctest::Approx(std::exp(-10.0 * Ts)).epsilon(1e-12));
        CHECK((1.0 + a) - (b * g.Kp) - (b * g.Ki * Ts) == doctest::Approx(2.0 * std::exp(-5.0 * Ts) * std::cos(3.0 * Ts)).epsilon(1e-12));
        const auto lone = design::pid(plant, Ts, design::Poles{complex<double>{-5.0, 3.0}, complex<double>{-5.0, 3.0}});
        CHECK(lone.status == DesignStatus::invalid_parameter);
        CHECK(design::pid(plant, Ts, design::Poles{runtime_nan(), -5.0}).status == DesignStatus::invalid_parameter);
    }

    TEST_CASE("Poles refuses Re(s) ≥ 0, so a braced list of z-plane poles is refused") {
        // Read as s, {0.95, 0.96} is z = e^{0.95 Ts} > 1: an unstable loop that came back with success.
        constexpr auto plant = TransferFunction<1, 2>{.num = {1.0}, .den = {1.0, 1.0}};
        CHECK(design::pid(plant, 0.01, design::Poles{0.95, 0.96}).status == DesignStatus::invalid_parameter);
        CHECK(design::pid(plant, 0.01, design::Poles{0.8, 0.7, 0.6}).status == DesignStatus::invalid_parameter);
        CHECK(design::pid(plant, 0.01, design::Poles{0.0, 0.0}).status == DesignStatus::invalid_parameter);
        const auto axis = design::pid(plant, 0.01, design::Poles{complex<double>{0.0, 3.0}, complex<double>{0.0, -3.0}});
        CHECK(axis.status == DesignStatus::invalid_parameter);
    }

    TEST_CASE("Poles keeps the ZOH gain accurate in float for light friction") {
        // 1/(1e-4 s + b_visc) at Ts = 1e-4: x = b_visc Ts/1e-4 is tiny and 1 − e^(−x) cancels.
        for (const double bv : {1e-5, 1e-6, 1e-7}) {
            const auto f = design::pid(TransferFunction<1, 2, float>{.num = {1.0f}, .den = {1e-4f, static_cast<float>(bv)}}, 1e-4f, design::Poles{-1053.6f, -1053.6f});
            const auto d = design::pid(TransferFunction<1, 2>{.num = {1.0}, .den = {1e-4, bv}}, 1e-4, design::Poles{-1053.6, -1053.6});
            REQUIRE(f.success);
            CHECK(static_cast<double>(f.gains.Kp) == doctest::Approx(d.gains.Kp).epsilon(1e-5));
            CHECK(static_cast<double>(f.gains.Ki) == doctest::Approx(d.gains.Ki).epsilon(1e-4));
        }
    }

    TEST_CASE("Poles refuses K = 0 or Ts = 0 with zero gains") {
        const auto k0 = design::pid(TransferFunction<1, 2>{.num = {0.0}, .den = {1.0, 1.0}}, 0.1, design::Poles{-5.0, -10.0});
        CHECK(k0.status == DesignStatus::invalid_parameter);
        CHECK(k0.gains.Kp == 0.0);
        CHECK(k0.gains.Ki == 0.0);

        const auto ts0 = design::pid(TransferFunction<1, 2>{.num = {1.0}, .den = {1.0, 1.0}}, 0.0, design::Poles{-5.0, -10.0});
        CHECK(ts0.status == DesignStatus::invalid_sample_time);
        CHECK(ts0.gains.Kp == 0.0);
        CHECK(ts0.gains.Ki == 0.0);
    }

    TEST_CASE("Poles on an integrator plant uses b = b₀Ts/a₁") {
        // 1/s sampled at Ts: G(z) = Ts/(z − 1). PI poles at the images z = {0.5, 0.6}.
        constexpr double Ts = 0.01;
        constexpr auto   r = design::pid(TransferFunction<1, 2>{.num = {1.0}, .den = {1.0, 0.0}}, Ts, design::Poles{s_pole(0.5, Ts), s_pole(0.6, Ts)});
        static_assert(r.success);
        // (z − 1)² + Ts (Kp (z − 1) + Ki Ts z) = z² − 1.1 z + 0.3
        CHECK(1.0 - (Ts * r.gains.Kp) == doctest::Approx(0.3));
        CHECK((Ts * r.gains.Kp) + (Ts * Ts * r.gains.Ki) - 2.0 == doctest::Approx(-1.1));
    }
}

TEST_SUITE("PID Design - No seeded tracking time") {
    // Plant 1/(s + 1) at Ts = 0.1: a = e^(−0.1), b = 1 − a. Poles 0.95 and 0.96 have
    // p₁p₂ > a, so Kp = (a − p₁p₂)/b < 0.
    constexpr double kPlaceTs = 0.1;
    constexpr auto   kPlant = TransferFunction<1, 2>{.num = {1.0}, .den = {1.0, 1.0}};

    TEST_CASE("no rule seeds a tracking time: the integral holds against a clamp, as on the LQ laws") {
        constexpr auto placed = design::pid(kPlant, kPlaceTs, design::Poles{s_pole(0.95, kPlaceTs), s_pole(0.96, kPlaceTs)});
        static_assert(placed.success);
        static_assert(placed.gains.Kp < 0.0 && placed.gains.Ki > 0.0);
        CHECK(placed.gains.tracking_time == 0.0);
        CHECK(design::pid(kPlant, kTs, design::PolePlacement{.omega = 3.0}).gains.tracking_time == 0.0);
        CHECK(design::pid(TransferFunction<1, 3>{.num = {1.0}, .den = {0.05, 0.0, 0.0}}, kTs, design::PolePlacement{.omega = 8.0, .omega_i = 1.0}).gains.tracking_time == 0.0);
        for (const double wc : {1.5, 2.0, 8.0}) { // below ωc = tan 30° the plant lags too little: C would need Kp < 0
            for (const auto mode : {PIDMode::PI, PIDMode::PID}) {
                const auto g = loop_shape(kLag1, {.wc = wc, .mode = mode});
                REQUIRE(g.success);
                REQUIRE(g.gains.Ki > 0.0);
                CHECK(g.gains.tracking_time == 0.0);
            }
        }
        CHECK(design::pid(design::FOPDT{.K = 1.0, .L = 0.5, .tau = 2.0}, kTs, design::SIMC{.tau_c = 1.0}).gains.tracking_time == 0.0);
        CHECK(design::pid(design::FOPDT{.K = 1.0, .L = 0.5, .tau = 2.0}, kTs, design::AMIGO{}).gains.tracking_time == 0.0);
    }

    TEST_CASE("a Kp < 0 placement with a tracking time does not wind up under a clamp") {
        const auto                gains = design::pid(kPlant, kPlaceTs, design::Poles{s_pole(0.95, kPlaceTs), s_pole(0.96, kPlaceTs)}).gains;
        const double              tt = damp::max(damp::abs(gains.Kp / gains.Ki), kPlaceTs); // T_t = |T_i|, at least a sample
        const design::InputLimits lim{.lower = {-1.0}, .upper = {1.0}};
        PIController<double>      law{design::pid(kPlant, kPlaceTs, design::Poles{s_pole(0.95, kPlaceTs), s_pole(0.96, kPlaceTs)}, lim, design::AntiWindup{.tracking_time = tt})};
        PIController<double>      holding{design::pid(kPlant, kPlaceTs, design::Poles{s_pole(0.95, kPlaceTs), s_pole(0.96, kPlaceTs)}, lim)}; // the integral holds against the clamp
        REQUIRE(law.valid());
        REQUIRE(holding.valid());
        CHECK(law.gains().tracking_time == tt);

        // r = 10 needs u = 10 on a unit-gain plant; the clamp holds u at 1 and y settles at 1.
        const double a = std::exp(-kPlaceTs);
        const double b = 1.0 - a;
        double       y = 0.0;
        double       y_held = 0.0;
        double       at_limit = 0.0; // the held law's integral on the tick its output first met the limit
        bool         met = false;
        for (int k = 0; k < 400; ++k) {
            y = (a * y) + (b * law.step(10.0, y));
            const double u_held = holding.step(10.0, y_held);
            if (!met && u_held >= 1.0) {
                met = true;
                at_limit = holding.integral();
            }
            y_held = (a * y_held) + (b * u_held);
        }
        // Saturated equilibrium: ΔI = Kᵢ·Ts·e + (Ts/T_t)·(u_sat − u_unsat) = 0 gives
        // u_unsat − u_sat = Kᵢ·T_t·e; the stored I (after the tick's correction) is then
        // u_sat + (Kᵢ·T_t − Kp − Kᵢ·Ts)·e.
        const double e = 10.0 - y;
        CHECK(y == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(law.integral() == doctest::Approx(1.0 + (((gains.Ki * tt) - gains.Kp - (gains.Ki * kPlaceTs)) * e)).epsilon(1e-6));
        // With T_t = 0 the integral stops where the output met the limit and holds there: the
        // ticks of positive error that follow add nothing.
        REQUIRE(met);
        CHECK(y_held == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(holding.integral() == doctest::Approx(at_limit).epsilon(1e-12));
        CHECK(holding.integral() < 1.0 - (gains.Kp * 10.0)); // under u_sat − Kp·e at the first error

        // Back to an attainable reference: the back-calculating law leaves saturation within a
        // few T_t.
        int  ticks = 0;
        bool left = false;
        for (; ticks < 50 && !left; ++ticks) {
            const double u = law.step(0.5, y);
            left = u < 1.0;
            y = (a * y) + (b * u);
        }
        CHECK(left);
        CHECK(ticks < 20);
    }

    TEST_CASE("an I-only placement (Kp = 0) and a negative-gain plant deploy") {
        // s = −½ twice: z₁z₂ = e^(−kPlaceTs) = a ⇒ Kp = 0
        const auto g = design::pid(kPlant, kPlaceTs, design::Poles{-0.5, -0.5});
        REQUIRE(g.success);
        CHECK(g.gains.Kp == doctest::Approx(0.0));
        CHECK(g.gains.Ki > 0.0);
        constexpr auto neg = design::pid(TransferFunction<1, 2>{.num = {-1.0}, .den = {1.0, 1.0}}, kPlaceTs, design::Poles{s_pole(0.5, kPlaceTs), s_pole(0.3, kPlaceTs)});
        static_assert(neg.success);
        static_assert(neg.gains.Ki < 0.0);
    }
}

TEST_SUITE("PID Design - Double-integrator pole placement") {
    TEST_CASE("PD places a critically-damped pair on 1/(J s^2)") {
        constexpr double J = 0.05;
        constexpr double w = 8.0;
        constexpr auto   plant = TransferFunction<1, 3>{.num = {1.0}, .den = {J, 0.0, 0.0}};
        constexpr auto   r = design::pid(plant, kTs, design::PolePlacement{.omega = w}).gains;
        static_assert(r.Ki == 0.0);
        // J Tf s³ + J s² + (Kp Tf + Kd) s + Kp = J Tf (s + ω)²(s + 10ω)
        CHECK(r.Tf == doctest::Approx(1.0 / (12.0 * w)));
        CHECK(r.Kp == doctest::Approx(J * r.Tf * 10.0 * w * w * w));
        CHECK(r.Ki == doctest::Approx(0.0));
        CHECK(r.c == doctest::Approx(0.0));
        CHECK(r.tracking_time == 0.0);
        check_placed(plant, r, w, 1.0, 0.0);
    }

    TEST_CASE("PID places (s^2 + 2ζω s + ω^2)(s + ωi)") {
        constexpr double J = 0.05;
        constexpr double w = 8.0;
        constexpr double z = 1.0;
        constexpr double wi = w / 8.0;
        constexpr auto   plant = TransferFunction<1, 3>{.num = {1.0}, .den = {J, 0.0, 0.0}};
        constexpr auto   r = design::pid(plant, kTs, design::PolePlacement{.omega = w, .zeta = z, .omega_i = wi});
        static_assert(r.success);
        // The delivered loop, filter included: the requested poles and the filter pole at −10ω.
        check_placed(plant, r.gains, w, z, wi);
        const auto under = design::pid(plant, kTs, design::PolePlacement{.omega = w, .zeta = 0.7, .omega_i = wi});
        check_placed(plant, under.gains, w, 0.7, wi);
    }

    TEST_CASE("SetpointWeights given set both weights of a placement, which leave the poles") {
        constexpr auto plant = TransferFunction<1, 3>{.num = {1.0}, .den = {0.05, 0.0, 0.0}};
        const auto     on_measurement = design::pid(plant, kTs, design::PolePlacement{.omega = 8.0, .omega_i = 1.0});
        const auto     weighted = design::pid(plant, kTs, design::PolePlacement{.omega = 8.0, .omega_i = 1.0}, design::SetpointWeights{.b = 0.5});
        REQUIRE(weighted.success);
        CHECK(on_measurement.gains.c == 0.0);
        CHECK(weighted.gains.b == 0.5);
        CHECK(weighted.gains.c == 1.0); // the option's own c
        CHECK(weighted.gains.Kp == doctest::Approx(on_measurement.gains.Kp).epsilon(1e-14));
        check_placed(plant, weighted.gains, 8.0, 1.0, 1.0);
    }

    TEST_CASE("bad J / bandwidth / zeta / omega_i are refused with zero gains") {
        const auto placed = [](const TransferFunction<1, 3>& plant, const design::PolePlacement<double>& rule) { return design::pid(plant, kTs, rule).gains.Kp; };
        CHECK(placed(TransferFunction<1, 3>{.num = {1.0}, .den = {0.0, 0.0, 0.0}}, {.omega = 8.0}) == 0.0);
        CHECK(placed(TransferFunction<1, 3>{.num = {1.0}, .den = {0.05, 0.0, 0.0}}, {.omega = -1.0}) == 0.0);
        CHECK(placed(TransferFunction<1, 3>{.num = {1.0}, .den = {0.05, 0.0, 0.0}}, {.omega = 8.0, .zeta = 0.0}) == 0.0);
        CHECK(placed(TransferFunction<1, 3>{.num = {1.0}, .den = {0.05, 0.0, 0.0}}, {.omega = 8.0, .zeta = 1.0, .omega_i = -0.1}) == 0.0);
    }
}

TEST_SUITE("PID Design - PI Pole Placement (first-order plant)") {

    // Plant 1/(a1 s + a0) under PI closes to a1 s² + (a0+Kp) s + Ki, i.e.
    //   s² + 2ζω s + ω²  with  ω = √(Ki/a1),  ζ = (a0+Kp)/(2 a1 ω).
    TEST_CASE("places the closed-loop pole pair at omega_bw") {
        constexpr double a1 = 1.2e-3; // e.g. inductance L or inertia J
        constexpr double a0 = 4.0e-4; // e.g. resistance R or viscous b
        const double     omega = 2.0 * std::numbers::pi * 200.0;

        const auto g = design::pid(TransferFunction<1, 2>{.num = {1.0}, .den = {a1, a0}}, 1e-4, design::PolePlacement{.omega = omega}).gains;

        CHECK(g.Kp == doctest::Approx((2.0 * omega * a1) - a0));
        CHECK(g.Ki == doctest::Approx(a1 * omega * omega));
        CHECK(g.tracking_time == 0.0);

        const double omega_n = std::sqrt(g.Ki / a1);
        const double zeta = (a0 + g.Kp) / (2.0 * a1 * omega_n);
        CHECK(omega_n == doctest::Approx(omega)); // bandwidth as requested
        CHECK(zeta == doctest::Approx(1.0));      // critically damped (default)
    }

    TEST_CASE("honors a non-unity damping ratio") {
        constexpr double a1 = 5.0e-4, a0 = 0.0, omega = 300.0, zeta = 0.7;
        const auto       g = design::pid(TransferFunction<1, 2>{.num = {1.0}, .den = {a1, a0}}, 1e-4, design::PolePlacement{.omega = omega, .zeta = zeta}).gains;
        const double     omega_n = std::sqrt(g.Ki / a1);
        const double     z = (a0 + g.Kp) / (2.0 * a1 * omega_n);
        CHECK(omega_n == doctest::Approx(omega));
        CHECK(z == doctest::Approx(zeta));
    }

    // TF form: b0/(a1 s + a0) reduces to monic by dividing den by b0.
    TEST_CASE("a numerator gain b0 ≠ 1 divides through") {
        constexpr double a1 = 1.0e-3, a0 = 0.8, omega = 1500.0;
        constexpr double b0 = 2.5;
        const auto       from_gain = design::pid(TransferFunction<1, 2>{.num = {b0}, .den = {a1, a0}}, 1e-4, design::PolePlacement{.omega = omega}).gains;
        const auto       expected = design::pid(TransferFunction<1, 2>{.num = {1.0}, .den = {a1 / b0, a0 / b0}}, 1e-4, design::PolePlacement{.omega = omega}).gains;
        CHECK(from_gain.Kp == doctest::Approx(expected.Kp));
        CHECK(from_gain.Ki == doctest::Approx(expected.Ki));
    }

    // Degenerate plant (zero input gain, e.g. a default-constructed TF) must not
    // divide by zero — it is refused with zero gains.
    TEST_CASE("PolePlacement guards a zero numerator gain") {
        constexpr TransferFunction<1, 2, double> degenerate{}; // num={0}, den={0,0}
        const auto                               g = design::pid(degenerate, kTs, design::PolePlacement{.omega = 1000.0});
        CHECK(g.status == DesignStatus::invalid_parameter);
        CHECK(g.gains.Kp == doctest::Approx(0.0));
        CHECK(g.gains.Ki == doctest::Approx(0.0));
    }

    TEST_CASE("PolePlacement refuses a pair the plant pole already outruns (Kp ≤ 0)") {
        // 1/(L s + R): Kp = 2ζωL − R ≤ 0 when the requested pair is slower than R/L.
        constexpr double L = 1.0e-3, R = 2.0;
        const auto       slow = design::pid(TransferFunction<1, 2>{.num = {1.0}, .den = {L, R}}, 1e-5, design::PolePlacement{.omega = 500.0});
        CHECK_FALSE(slow.success);
        CHECK(slow.status == DesignStatus::not_placeable);
        CHECK(slow.gains.Kp == 0.0);
        CHECK(design::pid(TransferFunction<1, 2>{.num = {1.0}, .den = {L, R}}, 1e-5, design::PolePlacement{.omega = 1001.0}).success);

        // A PI places two poles: an integrator pole on a first-order plant is refused.
        const auto third = design::pid(TransferFunction<1, 2>{.num = {1.0}, .den = {L, R}}, 1e-5, design::PolePlacement{.omega = 5000.0, .omega_i = 100.0});
        CHECK(third.status == DesignStatus::invalid_parameter);
    }

    TEST_CASE("second-order plant with damping: Kd and Kp net out the plant's own terms") {
        // 1/(J s² + b s): PD closes to J s² + (b + Kd) s + Kp.
        constexpr double J = 0.05, b = 0.1, w = 8.0;
        constexpr auto   plant = TransferFunction<1, 3>{.num = {1.0}, .den = {J, b, 0.0}};
        constexpr auto   r = design::pid(plant, kTs, design::PolePlacement{.omega = w});
        static_assert(r.success);
        check_placed(plant, r.gains, w, 1.0, 0.0);

        // More plant damping than the pair asks for needs Kd < 0: refused.
        const auto over = design::pid(TransferFunction<1, 3>{.num = {1.0}, .den = {J, 10.0, 0.0}}, kTs, design::PolePlacement{.omega = w});
        CHECK(over.status == DesignStatus::not_placeable);
    }

    TEST_CASE("state-space plants give the transfer-function plants' gains") {
        // First order: 2.5/(1e-3 s + 0.8) as ẋ = −800 x + u, y = 2500 x.
        constexpr double                 b0 = 2.5, a1 = 1.0e-3, a0 = 0.8, omega = 1500.0;
        constexpr double                 Ts = 1.0e-4;
        constexpr StateSpace             first{.A = {{-a0 / a1}}, .B = {{1.0}}, .C = {{b0 / a1}}};
        constexpr TransferFunction<1, 2> first_tf{.num = {b0}, .den = {a1, a0}};
        constexpr TransferFunction<1, 2> first_monic{.num = {b0 / a1}, .den = {1.0, a0 / a1}}; // ss2tf of `first`
        constexpr auto                   pi_ss = design::pid(first, Ts, design::PolePlacement{.omega = omega});
        static_assert(pi_ss.success);
        const auto pi_monic = design::pid(first_monic, Ts, design::PolePlacement{.omega = omega});
        CHECK(pi_ss.gains.Kp == doctest::Approx(pi_monic.gains.Kp).epsilon(1e-14)); // the SS plant forwards ss2tf's TF unchanged (a constexpr and a run-time copy: a few ulps)
        CHECK(pi_ss.gains.Ki == doctest::Approx(pi_monic.gains.Ki).epsilon(1e-14));
        const auto pi_tf = design::pid(first_tf, Ts, design::PolePlacement{.omega = omega});
        CHECK(pi_ss.gains.Kp == doctest::Approx(pi_tf.gains.Kp).epsilon(1e-12));
        CHECK(pi_ss.gains.Ki == doctest::Approx(pi_tf.gains.Ki).epsilon(1e-12));

        // Sampled counterpart on the same continuous plant.
        const auto dpi_ss = design::pid(first, Ts, design::Poles{-1000.0, -2000.0});
        const auto dpi_tf = design::pid(first_tf, Ts, design::Poles{-1000.0, -2000.0});
        REQUIRE(dpi_ss.success);
        CHECK(dpi_ss.gains.Kp == doctest::Approx(dpi_tf.gains.Kp).epsilon(1e-12));
        CHECK(dpi_ss.gains.Ki == doctest::Approx(dpi_tf.gains.Ki).epsilon(1e-12));
        const auto dpid_ss = design::pid(first, Ts, design::Poles{-1000.0, -2000.0, -3000.0});
        const auto dpid_tf = design::pid(first_tf, Ts, design::Poles{-1000.0, -2000.0, -3000.0});
        REQUIRE(dpid_ss.success);
        CHECK(dpid_ss.gains.Kd == doctest::Approx(dpid_tf.gains.Kd).epsilon(1e-12));

        // Second order: mass-spring-damper 1/(m s² + b s + k) with states (x, v).
        constexpr double                 m = 0.05, b = 0.1, k = 2.0, w = 40.0;
        constexpr StateSpace             msd{.A = {{0.0, 1.0}, {-k / m, -b / m}}, .B = {{0.0}, {1.0 / m}}, .C = {{1.0, 0.0}}};
        constexpr TransferFunction<1, 3> msd_tf{.num = {1.0}, .den = {m, b, k}};
        const auto                       pid_ss = design::pid(msd, kTs, design::PolePlacement{.omega = w, .omega_i = w / 8.0});
        const auto                       pid_tf = design::pid(msd_tf, kTs, design::PolePlacement{.omega = w, .omega_i = w / 8.0});
        REQUIRE(pid_ss.success);
        CHECK(pid_ss.gains.Kp == doctest::Approx(pid_tf.gains.Kp).epsilon(1e-12));
        CHECK(pid_ss.gains.Ki == doctest::Approx(pid_tf.gains.Ki).epsilon(1e-12));
        CHECK(pid_ss.gains.Kd == doctest::Approx(pid_tf.gains.Kd).epsilon(1e-12));
        CHECK(pid_ss.gains.Tf == doctest::Approx(pid_tf.gains.Tf).epsilon(1e-12));

        // The same plant in mixed coordinates is still all-pole: C·B cancels to rounding residue.
        const auto mixed = similarity_transform(msd, Matrix{{1000.0, 0.2}, {0.3, 1.7}});
        REQUIRE(mixed.success);
        const auto pid_mixed = design::pid(mixed.sys, kTs, design::PolePlacement{.omega = w, .omega_i = w / 8.0});
        REQUIRE(pid_mixed.success);
        CHECK(pid_mixed.gains.Kp == doctest::Approx(pid_tf.gains.Kp).epsilon(1e-9));
        CHECK(pid_mixed.gains.Kd == doctest::Approx(pid_tf.gains.Kd).epsilon(1e-9));
    }

    TEST_CASE("state-space plants with a finite zero, D ≠ 0 or a sample time are refused") {
        constexpr double m = 0.05, b = 0.1, k = 2.0;
        // C = [1, 0.01]: y = x + 0.01 v, so G(s) = (0.01 s + 1)/(m s² + b s + k) has a zero at −100.
        constexpr StateSpace zero{.A = {{0.0, 1.0}, {-k / m, -b / m}}, .B = {{0.0}, {1.0 / m}}, .C = {{1.0, 0.01}}};
        const auto           with_zero = design::pid(zero, kTs, design::PolePlacement{.omega = 40.0});
        CHECK_FALSE(with_zero.success);
        CHECK(with_zero.status == DesignStatus::invalid_parameter);
        CHECK(with_zero.gains.Kp == 0.0);

        constexpr StateSpace feedthrough{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}, .D = {{0.5}}};
        CHECK(design::pid(feedthrough, kTs, design::PolePlacement{.omega = 10.0}).status == DesignStatus::invalid_parameter);
        CHECK(design::pid(feedthrough, 0.01, design::Poles{-50.0, -90.0}).status == DesignStatus::invalid_parameter);

        constexpr StateSpace no_gain{.A = {{-1.0}}, .B = {{1.0}}, .C = {{0.0}}};
        CHECK(design::pid(no_gain, kTs, design::PolePlacement{.omega = 10.0}).status == DesignStatus::invalid_parameter);

        constexpr StateSpace sampled{.A = {{0.9}}, .B = {{1.0}}, .C = {{1.0}}, .Ts = 0.01};
        CHECK(design::pid(sampled, kTs, design::PolePlacement{.omega = 10.0}).status == DesignStatus::invalid_sample_time);
        CHECK(design::pid(sampled, 0.01, design::Poles{-50.0, -90.0}).status == DesignStatus::invalid_sample_time);
    }
}

TEST_SUITE("PID Design - Gains you have") {
    TEST_CASE("Standard is Parallel with Ki = Kp/Ti and Kd = Kp·Td") {
        constexpr auto standard = design::pid(design::Standard{.Kp = 4.0, .Ti = 0.5, .Td = 0.0025}, 1e-3);
        constexpr auto parallel = design::pid(design::Parallel{.Kp = 4.0, .Ki = 8.0, .Kd = 0.01}, 1e-3);
        static_assert(standard.success);
        CHECK(standard.gains.Kp == parallel.gains.Kp);
        CHECK(standard.gains.Ki == doctest::Approx(parallel.gains.Ki).epsilon(1e-15));
        CHECK(standard.gains.Kd == doctest::Approx(parallel.gains.Kd).epsilon(1e-15));
        // Ti = 0 is no integral; a negative time is refused.
        CHECK(design::pid(design::Standard{.Kp = 4.0, .Ti = 0.0, .Td = 0.01}, 1e-3).gains.Ki == 0.0);
        CHECK(design::pid(design::Standard{.Kp = 4.0, .Ti = -0.5}, 1e-3).status == DesignStatus::invalid_parameter);
        CHECK(design::pid(design::Standard{.Kp = 4.0, .Ti = 0.5, .Td = -0.01}, 1e-3).status == DesignStatus::invalid_parameter);
        CHECK(design::pid(design::Standard{.Kp = runtime_nan(), .Ti = 0.5}, 1e-3).status == DesignStatus::invalid_parameter);
    }

    TEST_CASE("Parallel deploys the gains as given, and refuses a non-finite one") {
        constexpr auto res = design::pid(design::Parallel{.Kp = 2.0, .Ki = 5.0, .Kd = 0.1}, 0.01);
        static_assert(res.success);
        CHECK(res.gains.Kp == 2.0);
        CHECK(res.gains.Ki == 5.0);
        CHECK(res.gains.Kd == 0.1);
        CHECK(res.gains.b == 1.0);
        CHECK(res.gains.c == 1.0);
        CHECK(res.gains.Tf == 0.0);
        CHECK(res.gains.tracking_time == 0.0);
        CHECK(res.Ts == 0.01);
        CHECK(design::pid(design::Parallel{.Kp = runtime_nan()}, 0.01).status == DesignStatus::invalid_parameter);
        // Ki·Ts past the scalar's range: the runtime would refuse it, so the design does.
        CHECK(design::pid(design::Parallel{.Kp = 1.0, .Ki = 1e300}, 1e10).status == DesignStatus::invalid_parameter);
        // A negative gain is a reverse-acting law, not a refusal.
        CHECK(design::pid(design::Parallel{.Kp = -1.0, .Ki = 2.0}, 0.01).success);
    }

    TEST_CASE("Ts is the design's, refused unless finite and positive") {
        for (const double ts : {0.0, -0.01, runtime_nan()}) {
            CHECK(design::pid(design::Parallel{.Kp = 1.0}, ts).status == DesignStatus::invalid_sample_time);
            CHECK(design::pid(design::FOPDT{.K = 1.0, .L = 0.5, .tau = 2.0}, ts, design::SIMC{.tau_c = 1.0}).status == DesignStatus::invalid_sample_time);
            CHECK(design::pid(kLag1, ts, design::LoopShape{.wc = 2.0}).status == DesignStatus::invalid_sample_time);
        }
    }

    TEST_CASE("Design result converts to PIDController and ticks at Ts") {
        constexpr double Ts = 0.01;
        constexpr auto   result = design::pid(design::Parallel{.Kp = 6.0, .Ki = 6.0, .Kd = 1.5}, Ts);
        PIDController    controller(result.as<float>());
        CHECK(controller.gains().Kp == doctest::Approx(6.0f).epsilon(1e-4));
        CHECK(controller.gains().Ki == doctest::Approx(6.0f).epsilon(1e-4));
        CHECK(controller.gains().Kd == doctest::Approx(1.5f).epsilon(1e-4));
        // The tick's integral increment is Ki·Ts per unit error
        (void)controller.step(1.0f, 0.0f);
        CHECK(controller.integral() == doctest::Approx(static_cast<float>(6.0 * Ts)).epsilon(1e-4));
    }

    TEST_CASE("PIDGains::to_tf yields Kp + Ki/s + Kd*s/(1+Tf*s)") {
        // Ideal PID (Tf = 0): C(s) = (Kd s^2 + Kp s + Ki)/s, highest power first.
        constexpr auto tf = design::pid(design::Parallel{.Kp = 2.0, .Ki = 3.0, .Kd = 0.5}, kTs).gains.to_tf();
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

namespace {

/// A design deployed as a Law, and one refused at compile time.
inline constexpr auto kTemp = design::pid(design::FOPDT{.K = 2.5, .L = 4.0, .tau = 30.0}, 0.5, design::SIMC{.tau_c = 4.0}, design::InputLimits{.lower = {0.0}, .upper = {100.0}});
inline constexpr auto kRefused = design::pid(design::FOPDT{.K = 0.0, .L = 4.0, .tau = 30.0}, 0.5, design::SIMC{.tau_c = 4.0});
inline constexpr auto kNoTs = design::pid(design::Parallel{.Kp = 1.0}, 0.0);

} // namespace

TEST_SUITE("PID Design - Options and deploy") {
    TEST_CASE("each option reaches the result") {
        constexpr auto res = design::pid(
            design::Parallel{.Kp = 2.0, .Ki = 5.0, .Kd = 0.1},
            0.01,
            design::AntiWindup{.tracking_time = 0.2, .limits = {.lower = -3.0, .upper = 4.0}},
            design::SetpointWeights{.b = 0.4, .c = 0.0},
            design::InputLimits{.lower = {-12.0}, .upper = {12.0}},
            design::DerivativeFilter{.Tf = 0.02}
        );
        static_assert(res.success);
        CHECK(res.gains.Tf == 0.02);
        CHECK(res.gains.b == 0.4);
        CHECK(res.gains.c == 0.0);
        CHECK(res.gains.tracking_time == 0.2);
        CHECK(res.limits.input.lower(0) == -12.0);
        CHECK(res.limits.input.upper(0) == 12.0);
        CHECK(res.limits.integral.lower(0) == -3.0);
        CHECK(res.limits.integral.upper(0) == 4.0);
        // An end left out is unbounded, and no option leaves everything unbounded.
        constexpr auto upper_only = design::pid(design::Parallel{.Kp = 1.0}, 0.01, design::InputLimits{.upper = {5.0}});
        CHECK(upper_only.limits.input.lower(0) == -std::numeric_limits<double>::max());
        CHECK(upper_only.limits.input.upper(0) == 5.0);
        constexpr auto plain = design::pid(design::Parallel{.Kp = 1.0}, 0.01);
        CHECK(plain.limits.input.upper(0) == std::numeric_limits<double>::max());
        CHECK(plain.limits.integral.upper(0) == std::numeric_limits<double>::max());
        // The rules take them too.
        const auto simc = design::pid(design::FOPDT{.K = 1.0, .L = 0.5, .tau = 2.0}, kTs, design::SIMC{.tau_c = 1.0, .mode = PIDMode::PID}, design::DerivativeFilter{.Tf = 0.05}, design::AntiWindup{.tracking_time = 0.5}, design::InputLimits{.lower = {0.0}, .upper = {1.0}});
        REQUIRE(simc.success);
        CHECK(simc.gains.Tf == 0.05);
        CHECK(simc.gains.tracking_time == 0.5);
        CHECK(simc.limits.input.lower(0) == 0.0);
    }

    TEST_CASE("an option's number outside its domain is refused") {
        const double nan = runtime_nan();
        const auto   refused = [](const auto&... opts) { return design::pid(design::Parallel{.Kp = 1.0, .Ki = 1.0}, 0.01, opts...).status; };
        CHECK(refused(design::DerivativeFilter{.Tf = -0.01}) == DesignStatus::invalid_parameter);
        CHECK(refused(design::DerivativeFilter{.Tf = nan}) == DesignStatus::invalid_parameter);
        CHECK(refused(design::SetpointWeights{.b = nan}) == DesignStatus::invalid_parameter);
        CHECK(refused(design::SetpointWeights{.c = nan}) == DesignStatus::invalid_parameter);
        CHECK(refused(design::AntiWindup{.tracking_time = -1.0}) == DesignStatus::invalid_parameter);
        CHECK(refused(design::AntiWindup{.tracking_time = nan}) == DesignStatus::invalid_parameter);
        CHECK(refused(design::AntiWindup{.limits = {.lower = 0.5, .upper = -0.5}}) == DesignStatus::invalid_parameter);
        CHECK(refused(design::InputLimits{.lower = {1.0}, .upper = {-1.0}}) == DesignStatus::invalid_parameter);
        CHECK(refused(design::InputLimits{.lower = {nan}, .upper = {1.0}}) == DesignStatus::invalid_parameter);
        CHECK(refused(design::InputLimits{.lower = {0.0}, .upper = {0.0}}) == DesignStatus::ok); // pinned, not refused
        CHECK(refused(design::AntiWindup{.tracking_time = 0.5}) == DesignStatus::ok);
    }

    TEST_CASE("a constexpr design deploys as Law<kRes, float> and as a PIDController") {
        static_assert(kTemp.success);
        static constinit Law<kTemp, float> heater;
        static_assert(std::is_same_v<Law<kTemp, float>::Runtime, PIDController<float, PIDMode::PI>>); // SIMC's PI: the narrowest mode
        static constinit PIController<float> same = kTemp.as<float>();
        CHECK(heater.limits().input.upper(0) == 100.0f);
        CHECK(heater.sample_time() == 0.5f);
        for (int k = 0; k < 5; ++k) {
            const float u = heater.step(50.0f, 20.0f);
            CHECK(u == doctest::Approx(same.step(50.0f, 20.0f)).epsilon(1e-6));
            CHECK(u <= 100.0f);
        }
    }

    TEST_CASE("a failed design is refused at compile time") {
        static_assert(!kRefused);
        static_assert(kRefused.status == DesignStatus::invalid_parameter);
        static_assert(!kNoTs);
        static_assert(kNoTs.status == DesignStatus::invalid_sample_time);
        // static constinit Law<kRefused, float> law; // does not compile: design_refused::invalid_parameter()
        PIDController<double> at_run_time{kRefused};
        CHECK_FALSE(at_run_time.valid());
    }
}
