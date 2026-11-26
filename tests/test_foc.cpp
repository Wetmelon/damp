// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <limits>
#include <numbers>

#include "damp/controllers/pid.hpp"
#include "damp/design/pid_design.hpp"
#include "damp/math/transforms.hpp"
#include "damp/motor/foc.hpp"
#include "damp/motor/spm.hpp"
#include "damp/systems/transfer_function.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

TEST_SUITE("FOC Controller") {
    constexpr float Ts = 1.0f / 8000.0f;

    TEST_CASE("Bandwidth tuning: Kp = ωL, Ki = ωR at every bandwidth") {
        const float Ld = 200e-6f;
        const float Lq = 350e-6f; // salient: Lq != Ld -> distinct gains
        const float R = 0.5f;
        const float omega_bw = 2.0f * std::numbers::pi_v<float> * 1000.0f; // 1 kHz, above R/Ld = 2500 rad/s

        damp::motor::FOController<float> foc({.d = Ld, .q = Lq}, R, 0.0f);
        CHECK(foc.tune(omega_bw));

        // Kp = ω·L, Ki = ω·R per axis: the PI zero at −Ki/Kp = −R/L cancels the plant pole.
        CHECK(foc.dctrl.gains().Kp == doctest::Approx(omega_bw * Ld));
        CHECK(foc.dctrl.gains().Ki == doctest::Approx(omega_bw * R));
        CHECK(foc.dctrl.gains().tracking_time == 0.0f); // the integral holds against the voltage circle
        CHECK(foc.qctrl.gains().Kp == doctest::Approx(omega_bw * Lq));
        CHECK(foc.qctrl.gains().Ki == doctest::Approx(omega_bw * R));
        CHECK(foc.dctrl.gains().Ki / foc.dctrl.gains().Kp == doctest::Approx(R / Ld));
        CHECK(foc.qctrl.gains().Ki / foc.qctrl.gains().Kp == doctest::Approx(R / Lq));

        // q-axis has larger inductance -> larger Kp, the same Ki
        CHECK(foc.qctrl.gains().Kp > foc.dctrl.gains().Kp);

        // Below the R/L pole the formula is the same.
        CHECK(foc.tune(10.0f));
        CHECK(foc.dctrl.gains().Kp == doctest::Approx(10.0f * Ld));
        CHECK(foc.qctrl.gains().Ki == doctest::Approx(10.0f * R));

        // An out-of-range input is refused and the gains kept.
        const float Kp_d = foc.dctrl.gains().Kp;
        CHECK_FALSE(foc.tune(-1.0f));
        CHECK(foc.dctrl.gains().Kp == Kp_d);
    }

    TEST_CASE("the tuned current loop is ω/(s + ω): a step reaches 1 − e^(−1) at t = 1/ω") {
        // L di/dt = v − R i under the PI on the measured period, sub-stepped: the closed loop is first order.
        constexpr double                  L = 200e-6;
        constexpr double                  R = 0.5;
        constexpr double                  omega = 6000.0; // above R/L = 2500 rad/s: the old placement put a pair here
        constexpr double                  dt = 1e-7;
        damp::motor::FOController<double> foc({.d = L, .q = L}, R, 0.0, omega);
        double                            iq = 0.0;
        const int                         n = static_cast<int>(std::lround(1.0 / (omega * dt)));
        for (int k = 0; k < n; ++k) {
            const auto cmd = foc.current_controller({.d = 0.0, .q = 1.0}, {.d = 0.0, .q = iq}, 0.0, dt);
            iq += dt * (cmd.Vdq.q - (R * iq)) / L;
        }
        CHECK(iq == doctest::Approx(1.0 - std::exp(-1.0)).epsilon(2e-3));
    }

    TEST_CASE("the damped pair of the old tuning is design::pid's PolePlacement on 1/(L s + R)") {
        constexpr float L = 200e-6f;
        constexpr float R = 0.5f;
        constexpr float omega = 6000.0f;
        constexpr auto  pair = damp::design::pid(damp::TransferFunction<1, 2, float>{.num = {1.0F}, .den = {L, R}}, Ts, damp::design::PolePlacement{.omega = omega, .zeta = 0.8f}, damp::design::SetpointWeights{.b = 0.0f});
        static_assert(pair.success);
        CHECK(pair.gains.Kp == doctest::Approx((2.0f * 0.8f * omega * L) - R));
        CHECK(pair.gains.Ki == doctest::Approx(L * omega * omega));
        CHECK(pair.gains.b == 0.0f);
        damp::motor::FOController<float> foc({.d = L, .q = L}, R, 0.0f);
        CHECK(foc.dctrl.set_gains(pair.gains)); // a placed pair runs on the FOC's PIs
        CHECK(foc.dctrl.gains().Kp == pair.gains.Kp);
    }

    TEST_CASE("Voltage-circle limit and anti-windup") {
        const float                      Vmax = 10.0f;
        damp::motor::FOController<float> foc({.d = 200e-6f, .q = 200e-6f}, 0.5f, 0.05f);
        foc.tune(2.0f * std::numbers::pi_v<float> * 800.0f);

        damp::DirectQuadrature<float> Idq_ref = {.d = 0.0f, .q = 50.0f}; // unreachable -> saturates
        damp::DirectQuadrature<float> Idq = {.d = 0.0f, .q = 0.0f};

        // Output magnitude must stay on/inside the voltage circle throughout.
        for (int k = 0; k < 200; ++k) {
            const auto  cmd = foc.current_controller(Idq_ref, Idq, 0.0f, Ts, Vmax);
            const float Vmag = std::sqrt((cmd.Vdq.d * cmd.Vdq.d) + (cmd.Vdq.q * cmd.Vdq.q));
            CHECK(Vmag <= doctest::Approx(Vmax).epsilon(1e-4f));
            CHECK(cmd.is_saturated); // unreachable ref -> always saturated here
        }

        // Anti-windup: the integrator must stop growing once saturated. Without
        // back-calculation it would ramp by Ki*e*Ts every step; with it, the
        // value over a further 200 saturated steps is essentially unchanged.
        const float i_settled = foc.qctrl.integral();
        for (int k = 0; k < 200; ++k) {
            (void)foc.current_controller(Idq_ref, Idq, 0.0f, Ts, Vmax);
        }
        CHECK(std::abs(foc.qctrl.integral() - i_settled) * foc.qctrl.gains().Ki < 1.0f);
    }

    TEST_CASE("current_controller returns a clamped voltage command") {
        damp::motor::FOController<float> foc({.d = 200e-6f, .q = 200e-6f}, 0.5f, 0.05f);
        foc.tune(2.0f * std::numbers::pi_v<float> * 800.0f);

        // Modest reference, ample voltage limit -> no saturation.
        const auto ok = foc.current_controller({.d = 0.0f, .q = 1.0f}, {}, 0.0f, Ts, 48.0f);
        CHECK_FALSE(ok.is_saturated);
        CHECK(ok.v_excess < 1.0f);

        // Huge reference on a tiny voltage limit -> command is clamped and flagged.
        foc.reset();
        damp::DirectQuadrature<float> big_ref = {.d = 0.0f, .q = 80.0f};
        damp::motor::DqCommand<float> st;
        for (int k = 0; k < 50; ++k) {
            st = foc.current_controller(big_ref, {}, 0.0f, Ts, 2.0f);
        }
        CHECK(st.is_saturated);
        CHECK(st.v_excess > 1.0f);
        CHECK(st.Vdq.abs() <= doctest::Approx(2.0f).epsilon(1e-4f));
    }

    TEST_CASE("Cross-axis decoupling feedforward is always applied") {
        // omega != 0, Id != 0 -> the q-axis sees an omega*Ld*Id decoupling term
        // even with PI gains and plant-inversion FF off.
        const float                      omega = 3000.0f;
        damp::motor::FOController<float> foc({.d = 200e-6f, .q = 350e-6f}, 0.5f, 0.1f);

        damp::DirectQuadrature<float> Idq = {.d = 5.0f, .q = 0.0f};
        const auto                    cmd = foc.current_controller(Idq, Idq, omega, Ts); // zero error -> PI contributes 0

        // Vq = omega*Ld*Id + omega*lambda ; Vd = -omega*Lq*Iq = 0
        CHECK(cmd.Vdq.q == doctest::Approx((omega * 200e-6f * 5.0f) + (omega * 0.1f)));
        CHECK(cmd.Vdq.d == doctest::Approx(-(omega * 350e-6f * 0.0f)));
    }

    TEST_CASE("plant_inversion_ff toggles the R*I + L*dI/dt term") {
        damp::motor::FOController<float> foc({.d = 200e-6f, .q = 200e-6f}, 0.5f, 0.0f);

        // Isolate the feedforward term: zero the auto-tuned PI gains so only the
        // model-inversion FF contributes to the command.
        foc.dctrl = {};
        foc.qctrl = {};

        damp::DirectQuadrature<float> Idq_ref = {.d = 0.0f, .q = 2.0f};
        damp::DirectQuadrature<float> Idq = {.d = 0.0f, .q = 0.0f};

        // Off (default): no FF and PI gains zeroed above -> zero command.
        CHECK(foc.current_controller(Idq_ref, Idq, 0.0f, Ts).Vdq.q == doctest::Approx(0.0f));

        // On: deadbeat L*dI/dt + R*I appears (R*Iq = 0 since Iq = 0).
        foc.plant_inversion_ff = true;
        const float expected_q = 200e-6f * ((2.0f - 0.0f) / Ts); // L * dIq/dt
        CHECK(foc.current_controller(Idq_ref, Idq, 0.0f, Ts).Vdq.q == doctest::Approx(expected_q));

        // R·I_meas term with nonzero measured current (L term zero when I_ref == I).
        Idq.q = 2.0f;
        CHECK(foc.current_controller(Idq_ref, Idq, 0.0f, Ts).Vdq.q == doctest::Approx(0.5f * 2.0f));

        // Ts ≤ 0: skip L·ΔI/Ts, still apply R·I (no /0).
        Idq.q = 0.0f;
        CHECK(foc.current_controller(Idq_ref, Idq, 0.0f, 0.0f).Vdq.q == doctest::Approx(0.0f));
        Idq.q = 1.0f;
        CHECK(foc.current_controller(Idq_ref, Idq, 0.0f, 0.0f).Vdq.q == doctest::Approx(0.5f));
    }

    TEST_CASE("PolePlacement on an R–L axis returns a clean 2-DOF design") {
        const double L = 250e-6;
        const double R = 0.35;
        const double wbw = 1500.0;
        const auto   placed = damp::design::pid(damp::TransferFunction<1, 2>{.num = {1.0}, .den = {L, R}}, 1e-4, damp::design::PolePlacement{.omega = wbw}); // zeta=1, b=1 defaults
        REQUIRE(placed.success);
        const auto& res = placed.gains;

        CHECK(res.Kp == doctest::Approx((2.0 * wbw * L) - R));
        CHECK(res.Ki == doctest::Approx(L * wbw * wbw));
        CHECK(res.Kd == doctest::Approx(0.0));
        CHECK(res.tracking_time == 0.0); // no seed: the integral holds against a clamp
        CHECK(res.b == doctest::Approx(1.0));

        // Regression: output/integrator limits must stay unbounded. A positional
        // brace-init slip {Kp, Ki, Kd, Ts} once landed Ts in the output-limit slot, which
        // would silently clamp any PIDController built straight from the result. The limits
        // are the options', unbounded unless given.
        constexpr double inf = std::numeric_limits<double>::max();
        CHECK(placed.limits.input.lower(0) == -inf);
        CHECK(placed.limits.input.upper(0) == inf);
        CHECK(placed.limits.integral.lower(0) == -inf);
        CHECK(placed.limits.integral.upper(0) == inf);

        // The setpoint weight propagates for the I-P structure.
        CHECK(damp::design::pid(damp::TransferFunction<1, 2>{.num = {1.0}, .den = {L, R}}, 1e-4, damp::design::PolePlacement{.omega = wbw}, damp::design::SetpointWeights{.b = 0.0}).gains.b == doctest::Approx(0.0));
    }

    TEST_CASE("tune() setpoint weight selects PI vs I-P") {
        const float                      wbw = 2.0f * std::numbers::pi_v<float> * 800.0f;
        damp::motor::FOController<float> foc({.d = 200e-6f, .q = 200e-6f}, 0.5f, 0.1f);

        foc.tune(wbw); // default b = 1 -> standard PI on both axes
        CHECK(foc.dctrl.gains().b == doctest::Approx(1.0f));
        CHECK(foc.qctrl.gains().b == doctest::Approx(1.0f));

        foc.tune(wbw, 0.0f); // b = 0 -> I-P
        CHECK(foc.dctrl.gains().b == doctest::Approx(0.0f));
        CHECK(foc.qctrl.gains().b == doctest::Approx(0.0f));
    }

    TEST_CASE("I-P removes the proportional step-kick") {
        // omega_e = 0, lambda = 0 -> all feedforward terms vanish, isolating the PI/I-P.
        const float wbw = 2.0f * std::numbers::pi_v<float> * 800.0f;

        auto make_foc = [&](float b) {
            damp::motor::FOController<float> foc({.d = 200e-6f, .q = 200e-6f}, 0.5f, 0.0f);
            foc.tune(wbw, b);
            return foc;
        };

        damp::DirectQuadrature<float> ref = {.d = 0.0f, .q = 1.0f};
        damp::DirectQuadrature<float> meas = {.d = 0.0f, .q = 0.0f};

        auto pi = make_foc(1.0f);
        auto ip = make_foc(0.0f);

        const float vq_pi = pi.current_controller(ref, meas, 0.0f, Ts).Vdq.q;
        const float vq_ip = ip.current_controller(ref, meas, 0.0f, Ts).Vdq.q;

        // Backward-Euler: on the first tick the integrator is seeded by e*Ts, so both
        // structures share the same integral term (Ki*e*Ts). The structural difference
        // is purely the proportional path: PI feeds Kp*(r - y), I-P feeds Kp*(0 - y),
        // so the step-kick PI carries over I-P is exactly Kp*r.
        CHECK((vq_pi - vq_ip) == doctest::Approx(pi.qctrl.gains().Kp * 1.0f)); // proportional step-kick
        CHECK(vq_ip == doctest::Approx(ip.qctrl.gains().Ki * Ts));             // I-P: only the integral acts
    }

    TEST_CASE("back-EMF saturation does not charge the current integrator") {
        damp::motor::FOController<float> foc({.d = 200e-6f, .q = 200e-6f}, 0.5f, 0.1f);
        foc.tune(1000.0f);
        damp::DirectQuadrature<float> idq{};
        const auto                    cmd = foc.current_controller(idq, idq, 8000.0f, Ts, 10.0f);
        CHECK(cmd.is_saturated);
        CHECK(cmd.Vdq.abs() == doctest::Approx(10.0f).epsilon(1e-4f));
        CHECK(std::abs(foc.qctrl.integral()) < 1.0f);
        CHECK(std::abs(foc.dctrl.integral()) < 1.0f);
    }

    TEST_CASE("non-positive Vmax forces zero voltage") {
        damp::motor::FOController<float> foc({.d = 200e-6f, .q = 200e-6f}, 0.5f, 0.0f);
        damp::DirectQuadrature<float>    ref{.d = 0.0f, .q = 5.0f};
        const auto                       cmd = foc.current_controller(ref, {}, 0.0f, Ts, 0.0f);
        CHECK(cmd.is_saturated);
        CHECK(cmd.Vdq.d == doctest::Approx(0.0f));
        CHECK(cmd.Vdq.q == doctest::Approx(0.0f));
    }

    TEST_CASE("track hands the current loops over from the applied voltage, feedforward included") {
        damp::motor::FOController<float> foc({.d = 200e-6f, .q = 300e-6f}, 0.5f, 0.01f);
        foc.tune(2.0f * std::numbers::pi_v<float> * 800.0f);
        const damp::DirectQuadrature<float> ref{.d = 0.0f, .q = 4.0f};
        const damp::DirectQuadrature<float> meas{.d = 0.3f, .q = 2.5f};
        const damp::DirectQuadrature<float> applied{.d = -1.5f, .q = 3.0f}; // what another source drove
        constexpr float                     omega_e = 900.0f;
        for (int k = 0; k < 5; ++k) {
            foc.track(applied, ref, meas, omega_e, Ts);
        }
        const auto cmd = foc.current_controller(ref, meas, omega_e, Ts);
        CHECK(cmd.Vdq.d == doctest::Approx(applied.d).epsilon(1e-4f));
        CHECK(cmd.Vdq.q == doctest::Approx(applied.q).epsilon(1e-4f));
    }

    TEST_CASE("SPM datasheet constant conversions") {
        using namespace damp::motor;
        using damp::motor::spm::flux_from_Kv;
        using damp::motor::spm::flux_from_torque_constant;
        using damp::motor::spm::iq_from_torque;
        using damp::motor::spm::motor_constant;
        using damp::motor::spm::torque_constant_from_flux;
        using damp::motor::spm::torque_constant_from_Kv;

        const double p = 7.0;
        const double lambda = 0.00295;

        const double Kt = torque_constant_from_flux(p, lambda);
        CHECK(Kt == doctest::Approx(1.5 * p * lambda));
        CHECK(flux_from_torque_constant(p, Kt) == doctest::Approx(lambda));
        CHECK(iq_from_torque(2.0, p, lambda) == doctest::Approx(2.0 / Kt));

        const double Kv = 270.0;
        CHECK(torque_constant_from_Kv(Kv) * Kv == doctest::Approx(8.2699).epsilon(1e-3));
        CHECK(flux_from_Kv(p, Kv) == doctest::Approx(flux_from_torque_constant(p, torque_constant_from_Kv(Kv))));
        CHECK(motor_constant(Kt, 0.039) == doctest::Approx(Kt / std::sqrt(0.039)));

        CHECK(voltage_circle_radius(24.0) == doctest::Approx(24.0 / std::sqrt(3.0)));
        CHECK(voltage_circle_radius(24.0, 0.9) == doctest::Approx(0.9 * 24.0 / std::sqrt(3.0)));

        // Unloaded SPMSM (Id = Iq = 0): ω_base = Vmax / λ
        const double Vmax = voltage_circle_radius(24.0);
        CHECK(Vmax / lambda == doctest::Approx(24.0 / std::sqrt(3.0) / lambda));
    }
}
