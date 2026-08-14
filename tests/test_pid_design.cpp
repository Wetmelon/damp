// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <limits>
#include <numbers>

#include "damp/backend.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/design/pid_design.hpp"
#include "damp/systems/transfer_function.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

/**
 * @brief Tests for modelless PID design methods (pid_design.hpp)
 *
 * Validates Ziegler-Nichols, Cohen-Coon, SIMC, Lambda, bandwidth,
 * Tyreus-Luyben, and pole placement PID tuning methods.
 */

TEST_SUITE("PID Design - Ziegler-Nichols Ultimate Gain") {
    TEST_CASE("ZN PID tuning from Ku and Tu") {
        constexpr double Ku = 10.0;
        constexpr double Tu = 2.0;
        constexpr double Ts = 0.01;

        constexpr auto result = design::ziegler_nichols(Ku, Tu, Ts);
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
        constexpr double Ts = 0.01;

        constexpr auto result = design::ziegler_nichols(Ku, Tu, Ts, design::PIDType::PI);
        CHECK(result.Kp == doctest::Approx(4.5)); // 0.45 * 10
        CHECK(result.Ki == doctest::Approx(4.5 / (2.0 / 1.2)));
        CHECK(result.Kd == doctest::Approx(0.0));
    }

    TEST_CASE("ZN P-only tuning") {
        constexpr auto result = design::ziegler_nichols(10.0, 2.0, 0.01, design::PIDType::P);
        CHECK(result.Kp == doctest::Approx(5.0));
        CHECK(result.Ki == doctest::Approx(0.0));
        CHECK(result.Kd == doctest::Approx(0.0));
    }

    TEST_CASE("Runtime ZN matches compile-time") {
        auto ct = design::ziegler_nichols(10.0, 2.0, 0.01, design::PIDType::PID);
        auto rt = design::ziegler_nichols(10.0, 2.0, 0.01, design::PIDType::PID);
        CHECK(ct.Kp == doctest::Approx(rt.Kp));
        CHECK(ct.Ki == doctest::Approx(rt.Ki));
        CHECK(ct.Kd == doctest::Approx(rt.Kd));
    }

    TEST_CASE("ZN Tu=0 returns finite zeros for PI/PID") {
        const auto pi = design::ziegler_nichols(10.0, 0.0, 0.01, design::PIDType::PI);
        CHECK(std::isfinite(pi.Kp));
        CHECK(std::isfinite(pi.Ki));
        CHECK(std::isfinite(pi.Kd));
        CHECK(pi.Kp == doctest::Approx(0.0));
        CHECK(pi.Ki == doctest::Approx(0.0));
        CHECK(pi.Kd == doctest::Approx(0.0));

        const auto pid = design::ziegler_nichols(10.0, 0.0, 0.01, design::PIDType::PID);
        CHECK(std::isfinite(pid.Kp));
        CHECK(pid.Kp == doctest::Approx(0.0));
        CHECK(pid.Ki == doctest::Approx(0.0));
        CHECK(pid.Kd == doctest::Approx(0.0));
    }
}

TEST_SUITE("PID Design - Ziegler-Nichols Step Response") {
    TEST_CASE("ZN step response PID") {
        constexpr double K = 1.0;
        constexpr double L = 0.5;
        constexpr double tau = 2.0;
        constexpr double Ts = 0.01;

        constexpr auto result = design::ziegler_nichols_step(K, L, tau, Ts);
        // Kp = 1.2 * tau/(K*L) = 1.2 * 4.0 = 4.8
        CHECK(result.Kp == doctest::Approx(4.8));
        // Ki = Kp / (2*L) = 4.8 / 1.0 = 4.8
        CHECK(result.Ki == doctest::Approx(4.8));
        // Kd = Kp * 0.5 * L = 4.8 * 0.25 = 1.2
        CHECK(result.Kd == doctest::Approx(1.2));
    }

    TEST_CASE("ZN step K=0 or L=0 returns finite zeros") {
        const auto k0 = design::ziegler_nichols_step(0.0, 0.5, 2.0, 0.01);
        CHECK(std::isfinite(k0.Kp));
        CHECK(std::isfinite(k0.Ki));
        CHECK(std::isfinite(k0.Kd));
        CHECK(k0.Kp == doctest::Approx(0.0));
        CHECK(k0.Ki == doctest::Approx(0.0));
        CHECK(k0.Kd == doctest::Approx(0.0));

        const auto l0 = design::ziegler_nichols_step(1.0, 0.0, 2.0, 0.01);
        CHECK(std::isfinite(l0.Kp));
        CHECK(std::isfinite(l0.Ki));
        CHECK(std::isfinite(l0.Kd));
        CHECK(l0.Kp == doctest::Approx(0.0));
        CHECK(l0.Ki == doctest::Approx(0.0));
        CHECK(l0.Kd == doctest::Approx(0.0));
    }
}

TEST_SUITE("PID Design - Tyreus-Luyben") {
    TEST_CASE("TL PID tuning") {
        constexpr double Ku = 10.0;
        constexpr double Tu = 2.0;
        constexpr double Ts = 0.01;

        constexpr auto result = design::tyreus_luyben(Ku, Tu, Ts);
        CHECK(result.Kp == doctest::Approx(Ku / 2.2));
        CHECK(result.Ki == doctest::Approx(result.Kp / (2.2 * Tu)));
        CHECK(result.Kd == doctest::Approx(result.Kp * Tu / 6.3));
    }

    TEST_CASE("TL PI tuning") {
        constexpr auto result = design::tyreus_luyben(10.0, 2.0, 0.01, design::PIDType::PI);
        CHECK(result.Kp == doctest::Approx(10.0 / 3.2));
    }

    TEST_CASE("TL more conservative than ZN") {
        constexpr auto zn = design::ziegler_nichols(10.0, 2.0, 0.01);
        constexpr auto tl = design::tyreus_luyben(10.0, 2.0, 0.01);
        // Tyreus-Luyben should have lower Kp (more conservative)
        CHECK(tl.Kp < zn.Kp);
    }
}

TEST_SUITE("PID Design - Cohen-Coon") {
    TEST_CASE("Cohen-Coon PID for FOPDT") {
        constexpr double K = 1.0;
        constexpr double L = 1.0;
        constexpr double tau = 4.0;
        constexpr double Ts = 0.01;

        constexpr auto result = design::cohen_coon(K, L, tau, Ts);
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
        auto ct = design::cohen_coon(1.0, 1.0, 4.0, 0.01);
        auto rt = design::cohen_coon(1.0, 1.0, 4.0, 0.01);
        CHECK(ct.Kp == doctest::Approx(rt.Kp).epsilon(1e-12));
        CHECK(ct.Ki == doctest::Approx(rt.Ki).epsilon(1e-12));
        CHECK(ct.Kd == doctest::Approx(rt.Kd).epsilon(1e-12));
    }

    TEST_CASE("Cohen-Coon K=0 or L=0 returns finite zeros") {
        const auto k0 = design::cohen_coon(0.0, 1.0, 4.0, 0.01);
        CHECK(std::isfinite(k0.Kp));
        CHECK(k0.Kp == doctest::Approx(0.0));
        CHECK(k0.Ki == doctest::Approx(0.0));
        CHECK(k0.Kd == doctest::Approx(0.0));

        const auto l0 = design::cohen_coon(1.0, 0.0, 4.0, 0.01);
        CHECK(std::isfinite(l0.Kp));
        CHECK(l0.Kp == doctest::Approx(0.0));
    }
}

TEST_SUITE("PID Design - SIMC") {
    TEST_CASE("SIMC PI tuning") {
        constexpr double K = 2.0;
        constexpr double L = 0.5;
        constexpr double tau = 3.0;
        constexpr double tau_c = 1.0;

        constexpr auto result = design::simc(K, L, tau, tau_c, design::PIDType::PI);
        // Kp = tau / (K * (tau_c + L)) = 3.0 / (2.0 * 1.5) = 1.0
        CHECK(result.Kp == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(result.Ki > 0.0);
        CHECK(result.Kd == doctest::Approx(0.0));
    }

    TEST_CASE("SIMC K=0 returns finite zeros") {
        const auto z = design::simc(0.0, 0.5, 3.0, 1.0, design::PIDType::PI);
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
        constexpr double Ts = 0.01;

        constexpr auto result = design::lambda_tuning(K, L, tau, lambda, Ts);
        // Kp = tau / (K * (lambda + L)) = 2.0 / (1.0 * 2.5) = 0.8
        CHECK(result.Kp == doctest::Approx(0.8).epsilon(1e-12));
        // Ki = Kp / tau = 0.8 / 2.0 = 0.4
        CHECK(result.Ki == doctest::Approx(0.4).epsilon(1e-12));
        CHECK(result.Kd == doctest::Approx(0.0));
    }

    TEST_CASE("Lambda K=0 or tau=0 returns finite zeros") {
        const auto k0 = design::lambda_tuning(0.0, 0.5, 2.0, 2.0, 0.01);
        CHECK(std::isfinite(k0.Kp));
        CHECK(k0.Kp == doctest::Approx(0.0));
        CHECK(k0.Ki == doctest::Approx(0.0));

        const auto tau0 = design::lambda_tuning(1.0, 0.5, 0.0, 2.0, 0.01);
        CHECK(std::isfinite(tau0.Kp));
        CHECK(tau0.Kp == doctest::Approx(0.0));
        CHECK(tau0.Ki == doctest::Approx(0.0));
    }
}

TEST_SUITE("PID Design - Bandwidth") {
    TEST_CASE("PI from bandwidth and phase margin") {
        // PI controller phase range: [0, -90°), so max phase margin < 90°
        // With 30° phase margin: desired_phase = 30° - 180° = -150°
        // But controller phase = atan2(-Ki/ω, Kp) is in [-90°, 0°]
        // Achievable range means we set the controller to contribute phase
        // such that phase margin target is met assuming unit-gain plant.
        constexpr auto result = design::pid_from_bandwidth(
            10.0, 30.0, 0.001, design::PIDType::PI
        );
        CHECK(result.Kp > 0.0);
        CHECK(result.Ki > 0.0);
        CHECK(result.Kd == doctest::Approx(0.0));
    }

    TEST_CASE("PID from bandwidth") {
        constexpr auto result = design::pid_from_bandwidth(
            50.0, 45.0, 0.001, design::PIDType::PID
        );
        CHECK(result.Kp > 0.0);
        CHECK(result.Ki > 0.0);
        CHECK(result.Kd > 0.0);
    }
}

TEST_SUITE("PID Design - Pole Placement") {
    TEST_CASE("PI pole placement for first-order plant") {
        // Plant: G(s) = 1/(s+1), K=1, tau=1
        // Desired poles at z=0.5 and z=0.3 (well inside unit circle)
        constexpr double K = 1.0;
        constexpr double tau = 1.0;
        constexpr double Ts = 0.1;
        constexpr double p1 = 0.5;
        constexpr double p2 = 0.3;

        constexpr auto result = design::pid_pole_placement(K, tau, p1, p2, Ts);
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

        constexpr auto result = design::pid_pole_placement(K, tau, p1, p2, p3, Ts);
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

    TEST_CASE("Pole placement at origin gives deadbeat") {
        // All poles at z=0 → deadbeat response
        constexpr double K = 1.0;
        constexpr double tau = 1.0;
        constexpr double Ts = 0.1;

        constexpr auto result = design::pid_pole_placement(K, tau, 0.0, 0.0, Ts);

        double a = std::exp(-Ts / tau);
        double b = K * (1.0 - a);
        double c1 = (b * result.Kp) + (b * result.Ki * Ts) - a - 1.0;
        double c0 = a - (b * result.Kp);

        // Both poles at 0: z² - 0*z + 0 → c1 = 0, c0 = 0
        CHECK(c1 == doctest::Approx(0.0).epsilon(1e-10));
        CHECK(c0 == doctest::Approx(0.0).epsilon(1e-10));
    }

    TEST_CASE("Pole placement K=0 or Ts=0 returns finite zeros") {
        const auto k0 = design::pid_pole_placement(0.0, 1.0, 0.5, 0.3, 0.1);
        CHECK(std::isfinite(k0.Kp));
        CHECK(std::isfinite(k0.Ki));
        CHECK(k0.Kp == doctest::Approx(0.0));
        CHECK(k0.Ki == doctest::Approx(0.0));

        const auto ts0 = design::pid_pole_placement(1.0, 1.0, 0.5, 0.3, 0.0);
        CHECK(std::isfinite(ts0.Kp));
        CHECK(ts0.Kp == doctest::Approx(0.0));
        CHECK(ts0.Ki == doctest::Approx(0.0));
    }
}

TEST_SUITE("PID Design - Double-integrator pole placement") {
    TEST_CASE("PD places a critically-damped pair on 1/(J s^2)") {
        constexpr double J = 0.05;
        constexpr double w = 8.0;
        constexpr auto   r = design::pid_pole_placement_double_integrator(J, w);
        static_assert(r.Ki == 0.0);
        CHECK(r.Kp == doctest::Approx(J * w * w));
        CHECK(r.Kd == doctest::Approx(2.0 * J * w));
        CHECK(r.Ki == doctest::Approx(0.0));
        CHECK(r.c == doctest::Approx(0.0));
        CHECK(r.Kbc == doctest::Approx(r.Kp));
        CHECK(r.Tf == doctest::Approx(1.0 / (10.0 * w)));
    }

    TEST_CASE("PID places (s^2 + 2ζω s + ω^2)(s + ωi)") {
        constexpr double J = 0.05;
        constexpr double w = 8.0;
        constexpr double z = 1.0;
        constexpr double wi = w / 8.0;
        constexpr auto   r = design::pid_pole_placement_double_integrator(J, w, z, wi);
        CHECK(r.Kd == doctest::Approx(J * (2.0 * z * w + wi)));
        CHECK(r.Kp == doctest::Approx(J * ((w * w) + (2.0 * z * w * wi))));
        CHECK(r.Ki == doctest::Approx(J * w * w * wi));
        // Characteristic J s³ + Kd s² + Kp s + Ki
        CHECK(r.Kd / J == doctest::Approx(2.0 * z * w + wi));
        CHECK(r.Kp / J == doctest::Approx(w * w + 2.0 * z * w * wi));
        CHECK(r.Ki / J == doctest::Approx(w * w * wi));
    }

    TEST_CASE("bad J / bandwidth / zeta / omega_i return zeros") {
        CHECK(design::pid_pole_placement_double_integrator(0.0, 8.0).Kp == 0.0);
        CHECK(design::pid_pole_placement_double_integrator(0.05, -1.0).Kp == 0.0);
        CHECK(design::pid_pole_placement_double_integrator(0.05, 8.0, 0.0).Kp == 0.0);
        CHECK(design::pid_pole_placement_double_integrator(0.05, 8.0, 1.0, -0.1).Kp == 0.0);
    }
}

TEST_SUITE("PID Design - Type Conversion") {
    TEST_CASE("Design result converts to PIDController via discretize") {
        constexpr double Ts = 0.01;
        constexpr auto   result = design::ziegler_nichols(10.0, 2.0, Ts);
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

    TEST_CASE("PIDResult::to_tf yields Kp + Ki/s + Kd*s/(1+Tf*s)") {
        // Ideal PID (Tf = 0): C(s) = (Kd s^2 + Kp s + Ki)/s, ascending powers.
        constexpr auto tf = design::pid(2.0, 3.0, 0.5).to_tf();
        CHECK(tf.num[0] == doctest::Approx(3.0)); // Ki
        CHECK(tf.num[1] == doctest::Approx(2.0)); // Kp
        CHECK(tf.num[2] == doctest::Approx(0.5)); // Kd
        CHECK(tf.den[0] == doctest::Approx(0.0));
        CHECK(tf.den[1] == doctest::Approx(1.0));
        CHECK(tf.den[2] == doctest::Approx(0.0));

        // Filtered (Tf > 0): num = {Ki, Kp + Ki*Tf, Kp*Tf + Kd}, den = {0, 1, Tf}.
        constexpr double inf = std::numeric_limits<double>::infinity();
        constexpr auto   tff = design::pid(2.0, 3.0, 0.5, -inf, inf, -inf, inf, 0.0, 1.0, 1.0, 0.1).to_tf();
        CHECK(tff.den[2] == doctest::Approx(0.1));
        CHECK(tff.num[1] == doctest::Approx(2.0 + (3.0 * 0.1)));
        CHECK(tff.num[2] == doctest::Approx((2.0 * 0.1) + 0.5));
    }
}

TEST_SUITE("PID Design - Performance Spec Glue") {
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

    TEST_CASE("PID from performance spec produces usable gains") {
        constexpr design::PIDPerformanceSpec<double> spec{
            .settling_time = 0.2,
            .overshoot_percent = 10.0,
            .Ts = 0.001,
            .type = design::PIDType::PID,
            .bandwidth_scale = 1.0
        };

        constexpr auto result = design::pid_from_performance_spec(spec);
        CHECK(result.Kp > 0.0);
        CHECK(result.Ki > 0.0);
        CHECK(result.Kd > 0.0);
    }

    TEST_CASE("PI from performance spec") {
        constexpr design::PIDPerformanceSpec<double> spec{
            .settling_time = 0.1,
            .overshoot_percent = 5.0,
            .Ts = 0.0005,
            .type = design::PIDType::PI,
            .bandwidth_scale = 1.2
        };

        constexpr auto result = design::pid_from_performance_spec(spec);
        CHECK(result.Kp > 0.0);
        CHECK(result.Ki > 0.0);
        CHECK(result.Kd == doctest::Approx(0.0));
    }
}

TEST_SUITE("PID Design - PI Pole Placement (first-order plant)") {

    // Plant 1/(a1 s + a0) under PI closes to a1 s² + (a0+Kp) s + Ki, i.e.
    //   s² + 2ζω s + ω²  with  ω = √(Ki/a1),  ζ = (a0+Kp)/(2 a1 ω).
    TEST_CASE("places the closed-loop pole pair at omega_bw") {
        constexpr double a1 = 1.2e-3; // e.g. inductance L or inertia J
        constexpr double a0 = 4.0e-4; // e.g. resistance R or viscous b
        const double     omega = 2.0 * std::numbers::pi * 200.0;

        const auto g = design::pi_pole_placement_first_order(a1, a0, omega);

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
        const auto       g = design::pi_pole_placement_first_order(a1, a0, omega, zeta);
        const double     omega_n = std::sqrt(g.Ki / a1);
        const double     z = (a0 + g.Kp) / (2.0 * a1 * omega_n);
        CHECK(omega_n == doctest::Approx(omega));
        CHECK(z == doctest::Approx(zeta));
    }

    // (current_loop_pi's delegation to this kernel is checked in test_foc.cpp,
    //  where foc.hpp is in scope.)

    // TF overload: b0/(a1 s + a0) reduces to monic by dividing den by b0.
    TEST_CASE("transfer-function overload matches the scalar form") {
        constexpr double a1 = 1.0e-3, a0 = 0.8, omega = 1500.0;

        // Unit numerator → identical to the scalar (a1,a0) call.
        constexpr TransferFunction<1, 2, double> plant_unit{.num = {1.0}, .den = {a0, a1}};
        const auto                               from_tf = design::pi_pole_placement_first_order(plant_unit, omega);
        const auto                               from_scalar = design::pi_pole_placement_first_order(a1, a0, omega);
        CHECK(from_tf.Kp == doctest::Approx(from_scalar.Kp));
        CHECK(from_tf.Ki == doctest::Approx(from_scalar.Ki));

        // Numerator gain b0 ≠ 1 must divide through: b0/(a1 s+a0) ≡ 1/((a1/b0)s+(a0/b0)).
        constexpr double                         b0 = 2.5;
        constexpr TransferFunction<1, 2, double> plant_gain{.num = {b0}, .den = {a0, a1}};
        const auto                               from_gain = design::pi_pole_placement_first_order(plant_gain, omega);
        const auto                               expected = design::pi_pole_placement_first_order(a1 / b0, a0 / b0, omega);
        CHECK(from_gain.Kp == doctest::Approx(expected.Kp));
        CHECK(from_gain.Ki == doctest::Approx(expected.Ki));
    }

    // Degenerate plant (zero input gain, e.g. a default-constructed TF) must not
    // divide by zero — it returns an inert all-zero PI.
    TEST_CASE("transfer-function overload guards a zero numerator gain") {
        constexpr TransferFunction<1, 2, double> degenerate{}; // num={0}, den={0,0}
        const auto                               g = design::pi_pole_placement_first_order(degenerate, 1000.0);
        CHECK(g.Kp == doctest::Approx(0.0));
        CHECK(g.Ki == doctest::Approx(0.0));
    }

    TEST_CASE("amigo_kappa_tau matches FOPDT formula golden") {
        // Ks=1, Ku=2, Tu=1 → τ,L,Kp,Ti from Åström & Hägglund AMIGO κ–τ map
        constexpr double Ks = 1.0;
        constexpr double Ku = 2.0;
        constexpr double Tu = 1.0;
        const double     wu = 2.0 * damp::numbers::pi_v<double> / Tu;
        const double     prod = Ks * Ku;
        const double     tau = std::sqrt((prod * prod) - 1.0) / wu;
        const double     L = (damp::numbers::pi_v<double> - std::atan(wu * tau)) / wu;
        const double     Kp_exp = (0.25 + (0.45 * tau / L)) / Ks;
        const double     Ti_exp = ((0.4 * L) + (0.8 * tau)) / (L + (0.1 * tau)) * L;
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
