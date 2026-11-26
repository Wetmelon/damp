// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/adrc.hpp"
#include "damp/design/status.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/zpk.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

/// One forward-Euler step of the continuous ESO with the β gains: the oracle for the small-ω_o·Ts
/// comparison.
template<size_t NX>
struct EulerADRC {
    design::ADRCResult<NX, double> r{};
    damp::array<double, NX + 1>    z{};
    double                         u_prev{0.0};

    double step(double ref, double y) {
        const double                e = z[0] - y;
        damp::array<double, NX + 1> dz{};
        for (size_t i = 0; i < NX; ++i) {
            dz[i] = z[i + 1] - (r.beta[i] * e);
        }
        dz[NX - 1] += r.gains.b0 * u_prev;
        dz[NX] = -(r.beta[NX] * e);
        for (size_t i = 0; i <= NX; ++i) {
            z[i] += r.Ts * dz[i];
        }
        const double u = ((r.Kp * (ref - z[0])) - (r.Kd * z[1]) - z[NX]) / r.gains.b0;
        u_prev = u;
        return u;
    }
};

/// Double integrator ÿ = b·u + d, advanced exactly over one period of held u.
struct DoubleIntegrator {
    double b{1.0};
    double d{0.0};
    double y{0.0};
    double v{0.0};

    void step(double u, double Ts) {
        const double a = (b * u) + d;
        y += (v * Ts) + (0.5 * a * Ts * Ts);
        v += a * Ts;
    }
};

/// Largest |u| and largest |u_new − u_euler| over @p N ticks of a 2nd-order ADRC, the discrete
/// ESO and the Euler oracle each in closed loop with its own double integrator; a load step at N/2.
struct Gap {
    double max_u{0.0};
    double max_du{0.0};
};

Gap euler_gap(double wc, double wo, double Ts, size_t N) {
    const auto                res = design::adrc<2>({.wc = wc, .wo = wo, .b0 = 1.0}, Ts);
    ADRCController<2, double> zoh{res};
    EulerADRC<2>              euler{.r = res};
    DoubleIntegrator          p_zoh{};
    DoubleIntegrator          p_euler{};
    Gap                       g{};
    for (size_t k = 0; k < N; ++k) {
        p_zoh.d = p_euler.d = (k < N / 2) ? 0.0 : 0.5;
        const double u_zoh = zoh.step(1.0, p_zoh.y);
        const double u_euler = euler.step(1.0, p_euler.y);
        g.max_u = damp::max(g.max_u, damp::abs(u_euler));
        g.max_du = damp::max(g.max_du, damp::abs(u_zoh - u_euler));
        p_zoh.step(u_zoh, Ts);
        p_euler.step(u_euler, Ts);
    }
    return g;
}

/// Coefficients of (z − λ)ⁿ, highest power first: what the ESO's characteristic polynomial should be.
template<size_t N>
damp::array<double, N + 1> repeated_root_poly(double lambda) {
    damp::array<double, N + 1> c{};
    c[0] = 1.0;
    for (size_t k = 0; k < N; ++k) { // multiply by (z − λ)
        for (size_t i = k + 1; i > 0; --i) {
            c[i] = c[i] - (lambda * c[i - 1]);
        }
    }
    return c;
}

} // namespace

TEST_SUITE("Active Disturbance Rejection Control (ADRC)") {
    TEST_CASE("1st-order ADRC gain computation") {
        constexpr size_t NX = 1;

        double wc = 5.0;  // Controller bandwidth
        double wo = 20.0; // FixedGainObserver bandwidth
        double b0 = 1.0;  // Plant gain

        auto adrc_result = design::adrc<NX>({.wc = wc, .wo = wo, .b0 = b0}, 1e-3);
        CHECK(adrc_result.success);

        // Expected ESO gains for 1st-order system with poles at -wo
        damp::array<double, NX + 1> expected_beta = {
            2 * wo,
            wo * wo,
        };

        CHECK(adrc_result.beta.size() == expected_beta.size());
        for (size_t i = 0; i < expected_beta.size(); ++i) {
            CHECK(adrc_result.beta[i] == doctest::Approx(expected_beta[i]).epsilon(1e-6));
        }

        // Unscaled PD: b0 is applied only in the runtime law (u0 − f̂)/b0
        double expected_Kp = wc;
        double expected_Kd = 0.0;

        CHECK(adrc_result.Kp == doctest::Approx(expected_Kp).epsilon(1e-6));
        CHECK(adrc_result.Kd == doctest::Approx(expected_Kd).epsilon(1e-6));
    }

    TEST_CASE("2nd-order ADRC gain computation") {
        constexpr size_t NX = 2;

        double wc = 5.0;  // Controller bandwidth
        double wo = 20.0; // FixedGainObserver bandwidth
        double b0 = 1.0;  // Plant gain

        auto adrc_result = design::adrc<NX>({.wc = wc, .wo = wo, .b0 = b0}, 1e-3);
        CHECK(adrc_result.success);

        // Expected ESO gains for 2nd-order system with poles at -wo
        damp::array<double, NX + 1> expected_beta = {
            3 * wo,
            3 * wo * wo,
            wo * wo * wo,
        };

        CHECK(adrc_result.beta.size() == expected_beta.size());
        for (size_t i = 0; i < expected_beta.size(); ++i) {
            CHECK(adrc_result.beta[i] == doctest::Approx(expected_beta[i]).epsilon(1e-6));
        }

        // Unscaled PD: b0 is applied only in the runtime law (u0 − f̂)/b0
        double expected_Kp = wc * wc;
        double expected_Kd = 2 * wc;

        CHECK(adrc_result.Kp == doctest::Approx(expected_Kp).epsilon(1e-6));
        CHECK(adrc_result.Kd == doctest::Approx(expected_Kd).epsilon(1e-6));
    }

    TEST_CASE("invalid b0 / bandwidths report success=false and zero gains") {
        auto bad_b0 = design::adrc<1>({.wc = 5.0, .wo = 20.0, .b0 = 0.0}, 1e-3);
        CHECK_FALSE(bad_b0.success);
        CHECK(bad_b0.Kp == 0.0);
        CHECK(bad_b0.gains.b0 == 0.0);

        auto bad_wc = design::adrc<1>({.wc = -1.0, .wo = 20.0, .b0 = 1.0}, 1e-3);
        CHECK_FALSE(bad_wc.success);

        auto bad_wo = design::adrc<2>({.wc = 5.0, .wo = 0.0, .b0 = 1.0}, 1e-3);
        CHECK_FALSE(bad_wo.success);

        // Runtime from a failed design must not produce Inf/NaN.
        ADRCController<1, double> ctrl(bad_b0);
        CHECK_FALSE(ctrl.valid());
        CHECK(ctrl.step(1.0, 0.0) == 0.0);
    }

    TEST_CASE("converting ctor preserves gains and ESO state") {
        ADRCController<1, double> d(design::adrc<1>({.wc = 20.0, .wo = 100.0, .b0 = 1.0}, 1e-3));
        (void)d.step(1.0, 0.2);
        (void)d.step(1.0, 0.4);

        ADRCController<1, float> f(d);
        CHECK(f.valid());
        // Next step matches within float precision (state + gains carried over).
        CHECK(static_cast<double>(f.step(1.0f, 0.5f)) == doctest::Approx(d.step(1.0, 0.5)).epsilon(1e-5));
    }

    TEST_CASE("ADRCResult::as propagates success") {
        constexpr auto ok = design::adrc<1>({.wc = 5.0, .wo = 20.0, .b0 = 1.0}, 1e-3);
        static_assert(ok.success);
        constexpr auto ok_f = ok.as<float>();
        static_assert(ok_f.success);
        CHECK(ok_f.Kp == doctest::Approx(static_cast<float>(ok.Kp)));

        constexpr auto bad = design::adrc<1>({.wc = 5.0, .wo = 20.0, .b0 = -1.0}, 1e-3);
        static_assert(!bad.success);
        static_assert(!bad.as<float>().success);
    }

    TEST_CASE("1st-order runtime: tracks setpoint and rejects a constant disturbance") {
        const double Ts = 1e-3;
        auto         ctrl = ADRCController<1, double>(design::adrc<1>({.wc = 20.0, .wo = 100.0, .b0 = 1.0}, Ts));

        // Plant: ẏ = −a·y + b·u + d. The −a·y term and the constant d are both
        // unmodeled, lumped into the ESO's total-disturbance estimate.
        const double a = 2.0, b = 1.0, d = 0.5, r = 1.0;
        double       y = 0.0;
        for (int k = 0; k < 20000; ++k) {
            const double u = ctrl.step(r, y);
            y += (-a * y + b * u + d) * Ts;
        }
        CHECK(y == doctest::Approx(r).epsilon(0.01)); // zero steady-state error
    }

    TEST_CASE("2nd-order runtime: double integrator tracks with disturbance rejection") {
        const double Ts = 1e-3;
        auto         ctrl = ADRCController<2, double>(design::adrc<2>({.wc = 10.0, .wo = 50.0, .b0 = 1.0}, Ts));

        // Plant: ÿ = b·u + d (double integrator + constant load disturbance).
        const double b = 1.0, d = 0.3, r = 1.0;
        double       y = 0.0, v = 0.0;
        for (int k = 0; k < 30000; ++k) {
            const double u = ctrl.step(r, y);
            v += (b * u + d) * Ts;
            y += v * Ts;
        }
        CHECK(y == doctest::Approx(r).epsilon(0.01));   // converged to setpoint
        CHECK(v == doctest::Approx(0.0).epsilon(1e-3)); // and settled (no residual velocity)
    }

    TEST_CASE("2nd-order runtime: b0 != 1 scales the first command as 1/b0") {
        const double Ts = 1e-3;
        const double wc = 10.0;
        const double wo = 50.0;
        auto         c1 = ADRCController<2, double>(design::adrc<2>({.wc = wc, .wo = wo, .b0 = 1.0}, Ts));
        auto         c2 = ADRCController<2, double>(design::adrc<2>({.wc = wc, .wo = wo, .b0 = 2.0}, Ts));
        const double u1 = c1.step(1.0, 0.0);
        const double u2 = c2.step(1.0, 0.0);
        // Fresh ESO: z = 0, so u = Kp r / b0 = wc² r / b0 (no extra 1/b0).
        CHECK(u1 == doctest::Approx(wc * wc).epsilon(1e-12));
        CHECK(u2 == doctest::Approx(u1 / 2.0).epsilon(1e-12));

        double       y = 0.0, v = 0.0;
        const double b = 2.0, d = 0.3, r = 1.0;
        auto         ctrl = ADRCController<2, double>(design::adrc<2>({.wc = wc, .wo = wo, .b0 = b}, Ts));
        for (int k = 0; k < 30000; ++k) {
            const double u = ctrl.step(r, y);
            v += (b * u + d) * Ts;
            y += v * Ts;
        }
        CHECK(y == doctest::Approx(r).epsilon(0.01));
        CHECK(v == doctest::Approx(0.0).epsilon(1e-3));
    }

    TEST_CASE("reset clears the observer state") {
        const double Ts = 1e-3;
        auto         ctrl = ADRCController<1, double>(design::adrc<1>({.wc = 20.0, .wo = 100.0, .b0 = 1.0}, Ts));
        for (int k = 0; k < 100; ++k) {
            (void)ctrl.step(1.0, 0.0);
        }
        ctrl.reset();
        // After reset the ESO state is zero, so the first command equals the
        // fresh-controller command for the same inputs.
        ADRCController<1, double> fresh(design::adrc<1>({.wc = 20.0, .wo = 100.0, .b0 = 1.0}, Ts));
        CHECK(ctrl.step(1.0, 0.0) == doctest::Approx(fresh.step(1.0, 0.0)));
    }

    TEST_CASE("the runtime ticks at the design's Ts; back_calculate overrides u_prev") {
        const double Ts = 1e-3;
        auto         res = design::adrc<1>({.wc = 20.0, .wo = 100.0, .b0 = 1.0}, Ts);
        REQUIRE(res.success);
        CHECK(res.Ts == doctest::Approx(Ts));
        ADRCController<1, double> ctrl{res}; // Ts comes from the result

        const double u3 = ctrl.step(1.0, 0.0);
        ctrl.reset();
        const double u2 = ctrl.step(1.0, 0.0);
        CHECK(u2 == doctest::Approx(u3));

        // Warm the ESO with a nonzero measurement so u_prev enters the update.
        ctrl.reset();
        for (int k = 0; k < 20; ++k) {
            (void)ctrl.step(1.0, 0.5);
        }
        const double u_free = ctrl.step(1.0, 0.5);
        ctrl.back_calculate(0.0); // plant saturated to 0
        // Same inputs; ESO next tick uses u_prev=0 instead of u_free
        CHECK(ctrl.last_command() == u_free); // the command, as on the other laws
        CHECK(ctrl.applied() == 0.0);         // and what the ESO takes as applied
        const double u_after = ctrl.step(1.0, 0.5);
        CHECK(u_after != doctest::Approx(u_free));
        CHECK(ctrl.last_command() == doctest::Approx(u_after));
        CHECK(ctrl.applied() == ctrl.last_command());
    }

    TEST_CASE("Kp Kd beta do not absorb 1/b0 — only the tick does") {
        constexpr auto a = design::adrc<2>({.wc = 8.0, .wo = 40.0, .b0 = 1.0}, 1e-3);
        constexpr auto b = design::adrc<2>({.wc = 8.0, .wo = 40.0, .b0 = 4.0}, 1e-3);
        static_assert(a.success && b.success);
        CHECK(a.Kp == doctest::Approx(b.Kp));
        CHECK(a.Kd == doctest::Approx(b.Kd));
        CHECK(a.beta[0] == doctest::Approx(b.beta[0]));
        CHECK(a.beta[2] == doctest::Approx(b.beta[2]));
        CHECK(a.gains.b0 == doctest::Approx(1.0));
        CHECK(b.gains.b0 == doctest::Approx(4.0));

        ADRCController<1, double> c1(design::adrc<1>({.wc = 12.0, .wo = 60.0, .b0 = 0.5}, 1e-3));
        CHECK(c1.step(1.0, 0.0) == doctest::Approx(12.0 / 0.5).epsilon(1e-12));
    }

    TEST_CASE("1st-order: y→r, u→(a r − d)/b, ESO matches f = −b0 u at rest") {
        const double Ts = 5e-4;
        const double wc = 25.0, wo = 125.0;
        const double a = 3.0, b = 2.5, d = -0.8, r = 1.25;
        auto         ctrl = ADRCController<1, double>(design::adrc<1>({.wc = wc, .wo = wo, .b0 = b}, Ts));

        double    y = 0.0;
        double    u = 0.0;
        const int N = 20000; // 10 s
        for (int k = 0; k < N; ++k) {
            u = ctrl.step(r, y);
            y += (-a * y + b * u + d) * Ts;
        }
        const double u_ss = (a * r - d) / b;
        CHECK(y == doctest::Approx(r).epsilon(0.005));
        CHECK(u == doctest::Approx(u_ss).epsilon(0.02));
        CHECK(ctrl.observer()[0] == doctest::Approx(y).epsilon(0.01));
        // ẏ=0 ⇒ f = −b0 u. ESO f̂ is z[1] for NX=1.
        CHECK(ctrl.observer()[1] == doctest::Approx(-b * u).epsilon(0.05));
    }

    TEST_CASE("2nd-order: y→r, v→0, u→−d/b, ESO velocity and f̂") {
        const double Ts = 5e-4;
        const double wc = 15.0, wo = 75.0;
        const double b = 3.0, d = 0.6, r = -0.75;
        auto         ctrl = ADRCController<2, double>(design::adrc<2>({.wc = wc, .wo = wo, .b0 = b}, Ts));

        double y = 0.0, v = 0.0, u = 0.0;
        for (int k = 0; k < 30000; ++k) {
            u = ctrl.step(r, y);
            v += (b * u + d) * Ts;
            y += v * Ts;
        }
        CHECK(y == doctest::Approx(r).epsilon(0.01));
        CHECK(v == doctest::Approx(0.0).epsilon(2e-3));
        CHECK(u == doctest::Approx(-d / b).epsilon(0.03));
        CHECK(ctrl.observer()[0] == doctest::Approx(y).epsilon(0.015));
        CHECK(ctrl.observer()[1] == doctest::Approx(v).epsilon(0.02));
        CHECK(ctrl.observer()[2] == doctest::Approx(-b * u).epsilon(0.08));
    }

    TEST_CASE("settling scales with wc (critically-damped 2nd-order law)") {
        const double Ts = 5e-4;
        const double r = 1.0;
        const auto   settle_time = [&](double wc) {
            const double wo = 5.0 * wc;
            auto         ctrl = ADRCController<2, double>(design::adrc<2>({.wc = wc, .wo = wo, .b0 = 1.0}, Ts));
            double       y = 0.0, v = 0.0;
            const int    hold = static_cast<int>(0.05 / Ts); // 50 ms inside the band
            int          inside = 0;
            const int    N = static_cast<int>(3.0 / Ts);
            for (int k = 0; k < N; ++k) {
                const double u = ctrl.step(r, y);
                v += u * Ts;
                y += v * Ts;
                if (damp::abs(y - r) < 0.05) {
                    ++inside;
                    if (inside >= hold) {
                        return static_cast<double>(k - hold) * Ts;
                    }
                } else {
                    inside = 0;
                }
            }
            return 1e9;
        };

        const double t10 = settle_time(10.0);
        const double t20 = settle_time(20.0);
        REQUIRE(t10 < 1.0);
        REQUIRE(t20 < 1.0);
        // Ideal 2% settle is ~4/wc; doubling wc should cut settle time clearly.
        CHECK(t20 < 0.7 * t10);
        CHECK(t10 < 8.0 / 10.0); // 10/wc would be sloppy; 8/wc still allows ESO transient
    }

    TEST_CASE("faster ESO takes a smaller hit from a load step") {
        // First-order plant: a load step hits ẏ immediately, so peak |y−r|
        // is a clean readout of ESO bandwidth.
        const double Ts = 5e-4;
        const double wc = 10.0, r = 1.0, a = 1.0;
        const auto   peak_after_load = [&](double wo) {
            auto      ctrl = ADRCController<1, double>(design::adrc<1>({.wc = wc, .wo = wo, .b0 = 1.0}, Ts));
            double    y = 0.0, d = 0.0;
            const int warm = static_cast<int>(1.2 / Ts);
            for (int k = 0; k < warm; ++k) {
                const double u = ctrl.step(r, y);
                y += (-a * y + u + d) * Ts;
            }
            REQUIRE(damp::abs(y - r) < 0.03);
            d = 4.0;
            double    peak = 0.0;
            const int N = static_cast<int>(1.0 / Ts);
            for (int k = 0; k < N; ++k) {
                const double u = ctrl.step(r, y);
                y += (-a * y + u + d) * Ts;
                peak = damp::max(peak, damp::abs(y - r));
            }
            return peak;
        };

        const double p_slow = peak_after_load(20.0); // wo = 2 wc
        const double p_fast = peak_after_load(80.0); // wo = 8 wc
        REQUIRE(p_slow > 0.03);
        REQUIRE(p_fast > 0.0);
        CHECK(p_fast < 0.75 * p_slow);
    }

    TEST_CASE("b0 mismatch of 2x still tracks and rejects (robustness)") {
        const double Ts = 5e-4;
        // Designed for b0=1; plant is twice as sensitive.
        auto         ctrl = ADRCController<2, double>(design::adrc<2>({.wc = 12.0, .wo = 60.0, .b0 = 1.0}, Ts));
        const double b_plant = 2.0, d = 0.4, r = 1.0;
        double       y = 0.0, v = 0.0;
        for (int k = 0; k < 25000; ++k) {
            const double u = ctrl.step(r, y);
            v += (b_plant * u + d) * Ts;
            y += v * Ts;
        }
        CHECK(y == doctest::Approx(r).epsilon(0.03));
        CHECK(v == doctest::Approx(0.0).epsilon(5e-3));
    }

    TEST_CASE("1st-order ADRC on a lagged plant (unmodeled pole)") {
        // Design model: ẏ = f + b0 u. Truth has a 8 ms input lag.
        const double Ts = 2e-4;
        const double tau = 0.008;
        auto         ctrl = ADRCController<1, double>(design::adrc<1>({.wc = 20.0, .wo = 100.0, .b0 = 1.0}, Ts));
        double       y = 0.0, x = 0.0;
        const double r = 1.0, a = 1.5, d = 0.3;
        for (int k = 0; k < 25000; ++k) {
            const double u = ctrl.step(r, y);
            x += ((u - x) / tau) * Ts;
            y += (-a * y + x + d) * Ts;
        }
        CHECK(y == doctest::Approx(r).epsilon(0.03));
    }

    TEST_CASE("NX=2 on a first-order plant: extra ESO state stays quiet") {
        const double Ts = 5e-4;
        auto         ctrl = ADRCController<2, double>(design::adrc<2>({.wc = 15.0, .wo = 70.0, .b0 = 1.0}, Ts));
        double       y = 0.0;
        const double r = 0.8, a = 2.0, d = 0.2;
        for (int k = 0; k < 20000; ++k) {
            const double u = ctrl.step(r, y);
            y += (-a * y + u + d) * Ts;
        }
        CHECK(y == doctest::Approx(r).epsilon(0.02));
        CHECK(ctrl.observer()[1] == doctest::Approx(0.0).epsilon(0.05)); // ẏ ≈ 0
    }

    TEST_CASE("saturation + back_calculate: ESO sees applied u, then recovers") {
        const double Ts = 5e-4;
        // ẏ = −y + u, r=1 ⇒ u_ss=1. Clamp at 0.25 so the loop saturates hard.
        const double umax = 0.25;
        auto         with_aw = ADRCController<1, double>(design::adrc<1>({.wc = 20.0, .wo = 80.0, .b0 = 1.0}, Ts));
        auto         no_aw = ADRCController<1, double>(design::adrc<1>({.wc = 20.0, .wo = 80.0, .b0 = 1.0}, Ts));
        double       y_aw = 0.0, y_na = 0.0;
        const int    sat_steps = static_cast<int>(0.8 / Ts);
        for (int k = 0; k < sat_steps; ++k) {
            const double u_aw = with_aw.step(1.0, y_aw);
            const double u_na = no_aw.step(1.0, y_na);
            const double applied_aw = damp::clamp(u_aw, -umax, umax);
            const double applied_na = damp::clamp(u_na, -umax, umax);
            with_aw.back_calculate(applied_aw);
            // no_aw deliberately not told about the clamp
            y_aw += (-y_aw + applied_aw) * Ts;
            y_na += (-y_na + applied_na) * Ts;
        }
        // During sat the plant cannot reach r=1. Drop the reference to a feasible 0.15.
        const double r2 = 0.15;
        for (int k = 0; k < static_cast<int>(1.5 / Ts); ++k) {
            const double u_aw = with_aw.step(r2, y_aw);
            const double u_na = no_aw.step(r2, y_na);
            const double applied_aw = damp::clamp(u_aw, -umax, umax);
            const double applied_na = damp::clamp(u_na, -umax, umax);
            with_aw.back_calculate(applied_aw);
            y_aw += (-y_aw + applied_aw) * Ts;
            y_na += (-y_na + applied_na) * Ts;
        }
        CHECK(y_aw == doctest::Approx(r2).epsilon(0.03));
        // Anti-windup should finish at least as close as the lied-to ESO.
        CHECK(damp::abs(y_aw - r2) <= damp::abs(y_na - r2) + 0.02);
    }

    TEST_CASE("float deploy of a 2nd-order design tracks the same plant") {
        const float    Ts = 1e-3f;
        constexpr auto art = design::adrc<2>({.wc = 12.0, .wo = 60.0, .b0 = 2.0}, 1e-3);
        static_assert(art.success);
        ADRCController<2, float> ctrl{art.as<float>()};
        float                    y = 0.0f, v = 0.0f;
        const float              b = 2.0f, d = 0.25f, r = 1.0f;
        for (int k = 0; k < 20000; ++k) {
            const float u = ctrl.step(r, y);
            v += (b * u + d) * Ts;
            y += v * Ts;
        }
        CHECK(static_cast<double>(y) == doctest::Approx(1.0).epsilon(0.03));
        CHECK(static_cast<double>(v) == doctest::Approx(0.0).epsilon(8e-3));
    }

    TEST_CASE("a Ts of 0, once the continuous design, is refused and leaves the runtime invalid") {
        // The ESO has no discrete map at Ts = 0: the runtime used to report valid() and emit
        // Kp·r/b0 whatever y was, an open-loop command. The design refuses it now.
        CHECK(design::adrc<2>({.wc = 10.0, .wo = 50.0, .b0 = 1.0}, 0.0).status == design::DesignStatus::invalid_sample_time);
        auto ctrl = ADRCController<2, double>(design::adrc<2>({.wc = 10.0, .wo = 50.0, .b0 = 1.0}, 0.0));
        CHECK_FALSE(ctrl.valid());
        for (const double y : {0.0, 0.5, 5.0, -3.0}) {
            CHECK(ctrl.step(1.0, y) == 0.0);
        }
        CHECK(ctrl.observer()[0] == 0.0);
    }

    TEST_CASE("a negative or non-finite Ts is refused") {
        CHECK(design::adrc<1>({.wc = 10.0, .wo = 40.0, .b0 = 1.0}, -1e-3).status == design::DesignStatus::invalid_sample_time);
        CHECK(design::adrc<2>({.wc = 10.0, .wo = 40.0, .b0 = 1.0}, std::numeric_limits<double>::quiet_NaN()).status == design::DesignStatus::invalid_sample_time);
        CHECK(design::adrc<1>({.wc = -10.0, .wo = 40.0, .b0 = 1.0}, 1e-3).status == design::DesignStatus::invalid_parameter);
    }

    TEST_CASE("non-finite design knobs fail closed") {
        CHECK_FALSE(design::adrc<1>({.wc = std::numeric_limits<double>::infinity(), .wo = 10.0, .b0 = 1.0}, 1e-3).success);
        CHECK_FALSE(design::adrc<2>({.wc = 10.0, .wo = std::numeric_limits<double>::quiet_NaN(), .b0 = 1.0}, 1e-3).success);
        CHECK_FALSE(design::adrc<1>({.wc = 10.0, .wo = 10.0, .b0 = -2.0}, 1e-3).success);
    }

    TEST_CASE("1st-order to_tf is (Kp β2 + (Kp β1 + β2)s) / (b0 s (s + β1 + Kp))") {
        constexpr auto r = design::adrc<1>({.wc = 8.0, .wo = 40.0, .b0 = 2.0}, 1e-3);
        static_assert(r.success);
        const auto   tf = r.to_tf();
        const double Kp = r.Kp;
        const double b0 = r.gains.b0;
        const double b1 = r.beta[0];
        const double b2 = r.beta[1];
        // Leverrier returns a monic denominator, highest power first: s² + (β1 + Kp)s + 0.
        CHECK(tf.den[0] == doctest::Approx(1.0));
        CHECK(std::abs(tf.den[2]) < 1e-12);
        CHECK(tf.den[1] == doctest::Approx(b1 + Kp));
        CHECK(tf.num[2] == doctest::Approx((Kp * b2) / b0));
        CHECK(tf.num[1] == doctest::Approx(((Kp * b1) + b2) / b0));
    }

    TEST_CASE("to_ss matches to_tf FRF; 2-DOF y-channel is −C") {
        constexpr auto              r = design::adrc<2>({.wc = 12.0, .wo = 60.0, .b0 = 1.5}, 1e-3);
        const auto                  C = r.to_ss();
        const auto                  tf = r.to_tf();
        const auto                  full = r.to_ss_2dof();
        const damp::complex<double> s{0.0, 7.0};
        const auto                  Gc = *eval_frf(C, s);
        const auto                  Gy = *eval_frf(full, s);
        CHECK(Gc(0, 0).real() == doctest::Approx(-Gy(0, 1).real()).epsilon(1e-9));
        CHECK(Gc(0, 0).imag() == doctest::Approx(-Gy(0, 1).imag()).epsilon(1e-9));
        // Highest power first: entry 3 − k multiplies sᵏ.
        damp::complex<double> num{tf.num[3], 0.0};
        damp::complex<double> den{tf.den[3], 0.0};
        damp::complex<double> sk{1.0, 0.0};
        for (size_t k = 1; k < 4; ++k) {
            sk = sk * s;
            num = num + (tf.num[3 - k] * sk);
            den = den + (tf.den[3 - k] * sk);
        }
        const auto H = num / den;
        CHECK(Gc(0, 0).real() == doctest::Approx(H.real()).epsilon(1e-8));
        CHECK(Gc(0, 0).imag() == doctest::Approx(H.imag()).epsilon(1e-8));
    }

    TEST_CASE("failed ADRC to_ss / to_tf are zero") {
        constexpr auto bad = design::adrc<1>({.wc = 0.0, .wo = 10.0, .b0 = 1.0}, 1e-3);
        CHECK_FALSE(bad.success);
        const auto ss = bad.to_ss();
        CHECK(ss.A(0, 0) == 0.0);
        CHECK(ss.B(0, 0) == 0.0);
        CHECK(ss.C(0, 0) == 0.0);
        CHECK(ss.D(0, 0) == 0.0);
    }
    TEST_CASE("small ω_o·Ts: the discrete ESO matches the forward-Euler ESO it replaces") {
        // The two are discretizations of one continuous ESO: per step they differ by O((ω_o·Ts)²),
        // over an observer time constant (1/(ω_o·Ts) steps) by O(ω_o·Ts). So the command gap,
        // relative to the largest command, scales with ω_o·Ts: 0.55·ω_o·Ts at both rates here,
        // bounded by ω_o·Ts.
        for (const double woTs : {1e-3, 1e-2}) {
            const double Ts = 1e-4;
            const double wo = woTs / Ts;
            const Gap    g = euler_gap(wo / 5.0, wo, Ts, static_cast<size_t>(40.0 / (wo * Ts)));
            CAPTURE(woTs);
            CAPTURE(g.max_du / g.max_u);
            CHECK(g.max_du <= woTs * g.max_u);
        }
    }

    TEST_CASE("large ω_o·Ts: Euler diverges, the discrete ESO stays stable") {
        // Euler's ESO poles are 1 − ω_o·Ts (repeated): outside the unit circle for ω_o·Ts > 2.
        // The discrete ESO's are e^{−ω_o·Ts} for any ω_o·Ts.
        const double              Ts = 1e-3;
        const double              wo = 2.5 / Ts; // ω_o·Ts = 2.5
        const double              wc = 20.0;
        const auto                res = design::adrc<2>({.wc = wc, .wo = wo, .b0 = 1.0}, Ts);
        ADRCController<2, double> zoh{res};
        EulerADRC<2>              euler{.r = res};
        DoubleIntegrator          p_zoh{.d = 0.3};
        DoubleIntegrator          p_euler{.d = 0.3};
        int                       diverged_at = -1; // first tick |u| passes 1e6 (then Euler stops)
        for (int k = 0; k < 2000; ++k) {
            const double u_zoh = zoh.step(1.0, p_zoh.y);
            p_zoh.step(u_zoh, Ts);
            if (diverged_at < 0) {
                const double u_euler = euler.step(1.0, p_euler.y);
                if (damp::abs(u_euler) > 1e6) {
                    diverged_at = k;
                }
                p_euler.step(u_euler, Ts);
            }
        }
        CHECK(diverged_at > 0); // the Euler ESO blew up (here within 5 ticks)
        CHECK(diverged_at < 100);
        CHECK(p_zoh.y == doctest::Approx(1.0).epsilon(1e-3)); // the discrete one tracks
        CHECK(p_zoh.v == doctest::Approx(0.0).epsilon(1e-3));
        CHECK(zoh.observer()[2] == doctest::Approx(0.3).epsilon(1e-3)); // f̂ = d at rest (b0 = b)
    }

    TEST_CASE("the discrete ESO's poles are all at e^{−ω_o·Ts}") {
        for (const double woTs : {1e-3, 0.1, 1.0, 3.0}) {
            const double Ts = 1e-3;
            const double lambda = std::exp(-woTs);
            CAPTURE(woTs);
            const auto r1 = design::adrc<1>({.wc = 10.0, .wo = woTs / Ts, .b0 = 2.0}, Ts);
            const auto r2 = design::adrc<2>({.wc = 10.0, .wo = woTs / Ts, .b0 = 2.0}, Ts);
            const auto c1 = ::damp::detail::leverrier(r1.as_eso().error_dynamics()).den;
            const auto c2 = ::damp::detail::leverrier(r2.as_eso().error_dynamics()).den;
            const auto e1 = repeated_root_poly<2>(lambda);
            const auto e2 = repeated_root_poly<3>(lambda);
            // Coefficients of order 1–3: rounding of a few products, 1e-12 absolute.
            for (size_t i = 0; i < 3; ++i) {
                CHECK(c1[i] == doctest::Approx(e1[i]).epsilon(1e-12).scale(1.0));
            }
            for (size_t i = 0; i < 4; ++i) {
                CHECK(c2[i] == doctest::Approx(e2[i]).epsilon(1e-12).scale(1.0));
            }
        }
        // The design carries the matrices; a Ts of 0, once the continuous design that left them zero, is refused.
        constexpr auto rc = design::adrc<2>({.wc = 10.0, .wo = 50.0, .b0 = 1.0}, 0.0);
        static_assert(!rc.success && rc.status == design::DesignStatus::invalid_sample_time);
        constexpr auto rd = design::adrc<2>({.wc = 10.0, .wo = 50.0, .b0 = 1.0}, 1e-3);
        static_assert(rd.L(0) > 0.0 && rd.L(0) < 1.0);
    }

    TEST_CASE("reset(z0) removes the first-command kick from ẑ = 0") {
        const auto                res = design::adrc<2>({.wc = 10.0, .wo = 50.0, .b0 = 1.0}, 1e-3);
        ADRCController<2, double> cold{res};
        ADRCController<2, double> warm{res};
        warm.reset(ColVec{1.0, 0.0, 0.0});
        CHECK(damp::abs(cold.step(1.0, 1.0)) > 100.0);
        CHECK(warm.step(1.0, 1.0) == doctest::Approx(0.0).scale(1.0));
    }

    TEST_CASE_TEMPLATE("InputLimits alone is the law clamped outside and told with back_calculate", T, float, double) {
        //! The law had no limit of its own: without a back_calculate each tick its ESO took the
        //! unclamped command as applied and put the difference in the disturbance estimate.
        constexpr double     Ts = 5e-4;
        constexpr double     umax = 0.25;
        const auto           res = design::adrc<1>({.wc = 20.0, .wo = 80.0, .b0 = 1.0}, Ts, design::InputLimits{.lower = {-umax}, .upper = {umax}});
        const auto           open = design::adrc<1>({.wc = 20.0, .wo = 80.0, .b0 = 1.0}, Ts);
        ADRCController<1, T> bounded{res.template as<T>()};
        ADRCController<1, T> outside{open.template as<T>()};
        // Two inlined copies of the law, which fast-math may round differently: a few ulps of its
        // largest term, Kp·(r − ẑ₁)/b0 with |r| ≤ 1, not of u (clang 23 reaches 9 of them in float).
        const T tol = T{64} * eps<T>() * static_cast<T>(res.Kp / res.gains.b0);
        T       ya{0};
        T       yb{0};
        int     saturated = 0;
        for (int k = 0; k < 4600; ++k) {
            const T r = (k < 1600) ? T{1} : static_cast<T>(0.15); // ẏ = −y + u: r = 1 is out of reach, 0.15 is not
            const T ua = bounded.step(r, ya);
            const T ub = damp::clamp(outside.step(r, yb), static_cast<T>(-umax), static_cast<T>(umax));
            outside.back_calculate(ub);
            saturated += (damp::abs(ua) >= static_cast<T>(umax)) ? 1 : 0;
            CHECK(damp::abs(ua) <= static_cast<T>(umax));
            CHECK(damp::abs(ua - ub) <= tol);
            CHECK(bounded.applied() == ua);
            ya += static_cast<T>(Ts) * (ua - ya);
            yb += static_cast<T>(Ts) * (ub - yb);
        }
        CHECK(saturated > 100);
        CHECK(static_cast<double>(ya) == doctest::Approx(0.15).epsilon(0.03));
    }

    TEST_CASE("limits that are not a box are refused; a constexpr design's Law carries the box, set_limits moves it") {
        CHECK(design::adrc<1>({.wc = 20.0, .wo = 80.0, .b0 = 1.0}, 5e-4, design::InputLimits{.lower = {1.0}, .upper = {-1.0}}).status == design::DesignStatus::invalid_parameter);
        auto res = design::adrc<1>({.wc = 20.0, .wo = 80.0, .b0 = 1.0}, 5e-4);
        res.limits.input = {.lower = 1.0, .upper = -1.0};
        ADRCController<1, double> inert{res};
        CHECK_FALSE(inert.valid());
        CHECK(inert.step(1.0, 0.0) == 0.0);

        static constexpr auto kRes = design::adrc<1>({.wc = 20.0, .wo = 80.0, .b0 = 1.0}, 5e-4, design::InputLimits{.lower = {-0.25}, .upper = {0.25}});
        Law<kRes, float>      law;
        CHECK(law.limits().input.upper(0) == 0.25F);
        CHECK(law.step(1.0F, 0.0F) == 0.25F);
        auto lim = law.limits();
        lim.input.upper(0) = 0.125F; // a rail that moves
        REQUIRE(law.set_limits(lim));
        CHECK(law.step(1.0F, 0.0F) == 0.125F);
        lim.input.lower(0) = 1.0F; // lower above upper
        CHECK_FALSE(law.set_limits(lim));
        CHECK(law.limits().input.lower(0) == -0.25F);
        CHECK(law.gains().wo == 80.0F);
    }

    TEST_CASE("track runs the ESO on the applied input: the first step is the law's own command on a current estimate") {
        //! Without a handover the ESO ran on whatever the law last commanded, and its f̂ took the
        //! other source's input for a disturbance; the first step after it kicked by that error.
        const double              Ts = 1e-3;
        const double              b = 1.0;
        const auto                res = design::adrc<2>({.wc = 10.0, .wo = 50.0, .b0 = b}, Ts);
        ADRCController<2, double> law{res};
        DoubleIntegrator          plant{.b = b, .d = -0.4};
        const double              u_hold = 0.4; // an operator holds the plant at rest against the load
        for (int k = 0; k < 4000; ++k) {
            law.track(u_hold, 1.0, plant.y);
            plant.step(u_hold, Ts);
        }
        CHECK(law.last_command() == u_hold);
        CHECK(law.applied() == u_hold);
        CHECK(law.observer()[2] == doctest::Approx(-b * u_hold).epsilon(1e-6)); // f̂ = −b₀u at rest under u
        const double u_law = law.step(0.0, plant.y);                            // at r = y, ẏ = 0: the law's own command is u
        CHECK(u_law == doctest::Approx(u_hold).epsilon(1e-6));
    }

    TEST_CASE("set_gains retunes the ESO and the PD at the design's Ts; a new b₀ keeps the command") {
        const double              Ts = 1e-3;
        const auto                res = design::adrc<1>({.wc = 20.0, .wo = 80.0, .b0 = 1.0}, Ts);
        ADRCController<1, double> law{res};
        double                    y = 0.0;
        for (int k = 0; k < 3000; ++k) {
            const double u = law.step(0.5, y);
            y += Ts * (u - y + 0.2); // ẏ = −y + u + d
        }
        auto g = law.gains();
        CHECK(g.wc == 20.0);
        g.b0 = 2.0;
        g.wo = 120.0;
        REQUIRE(law.set_gains(g));
        CHECK(law.gains().b0 == 2.0);
        const double u_prev = law.last_command();
        const double u_next = law.step(0.5, y);
        CHECK(u_next == doctest::Approx(u_prev).epsilon(1e-3)); // f̂ moved by −Δb₀·u: no jump at a steady estimate
        auto bad = g;
        bad.wo = -1.0;
        CHECK_FALSE(law.set_gains(bad));
        CHECK(law.gains().wo == 120.0);
        CHECK(law.sample_time() == Ts);
    }
}

namespace {
inline constexpr auto kAxis = design::adrc<2>({.wc = 10.0, .wo = 50.0, .b0 = 2.0}, 1e-3, design::InputLimits{.lower = {-1.0}, .upper = {1.0}});

/** @brief A law whose gains set_gains changes at run time. */
template<typename L>
concept Retunable = requires(L& l) { l.set_gains(l.gains()); };

/** @brief @p law on the double integrator ÿ = 2u − 0.3 for @p steps ticks of 1 ms, its references a square wave: the commands it returns. */
template<typename T, typename L>
void run_axis(L& law, DoubleIntegrator& plant, int steps, double* commands) {
    for (int k = 0; k < steps; ++k) {
        const T r = ((k / 1000) % 2 == 0) ? static_cast<T>(0.5) : static_cast<T>(-0.2);
        const T u = law.step(r, static_cast<T>(plant.y));
        commands[k] = static_cast<double>(u);
        plant.step(static_cast<double>(u), 1e-3);
    }
}
} // namespace

TEST_SUITE("ADRC Law") {
    TEST_CASE("Law<kRes, T> of a constexpr design::adrc: its PD gains and b₀ constants of the type, beside its ESO's maps") {
        static_assert(std::is_same_v<Law<kAxis, float>::Runtime, ADRCController<2, float>>);
        //! K_p, K_d, b₀ and T_s are the type's, as the ESO's Φ, Γ and L are: the Law stores the ESO's
        //! estimate, the input it holds as applied, the limits and the last command, and it has no set_gains.
        static_assert(sizeof(Law<kAxis, float>) < sizeof(ADRCController<2, float>));
        static_assert(!Retunable<Law<kAxis, float>>);
        static_assert(Retunable<ADRCController<2, float>>);

        //! In double the same commands to a few ulps of the operands, closed loop.
        static constinit Law<kAxis, double>          law;
        static constinit Law<kAxis, double>::Runtime same{kAxis};
        DoubleIntegrator                             plant_a{.b = 2.0, .d = -0.3};
        DoubleIntegrator                             plant_b{.b = 2.0, .d = -0.3};
        static double                                ua[5000];
        static double                                ub[5000];
        run_axis<double>(law, plant_a, 5000, ua);
        run_axis<double>(same, plant_b, 5000, ub);
        double worst = 0.0;
        for (int k = 0; k < 5000; ++k) {
            worst = std::max(worst, std::abs(ua[k] - ub[k]));
        }
        CAPTURE(worst);
        CHECK(worst < 1e-10);
        CHECK(plant_a.y == doctest::Approx(0.5).epsilon(0.02)); // the loop holds the last reference against the load
        CHECK(law.gains().b0 == 2.0);
        CHECK(law.sample_time() == 1e-3);

        //! In float the commands agree as far as the float ESO at 1 kHz resolves them: its third state
        //! is ÿ, ∝ 1/T_s² of the readings' rounding, and the constant maps sum in another order. The
        //! same holds of the ESO's constant maps before the PD gains moved, so it is not this change.
        static constinit Law<kAxis, float>          law_f;
        static constinit Law<kAxis, float>::Runtime same_f{kAxis};
        DoubleIntegrator                            plant_f{.b = 2.0, .d = -0.3};
        DoubleIntegrator                            plant_g{.b = 2.0, .d = -0.3};
        run_axis<float>(law_f, plant_f, 5000, ua);
        run_axis<float>(same_f, plant_g, 5000, ub);
        double worst_f = 0.0;
        for (int k = 0; k < 5000; ++k) {
            worst_f = std::max(worst_f, std::abs(ua[k] - ub[k]));
        }
        CAPTURE(worst_f);
        CHECK(worst_f < 2e-3);

        //! Its limits stay settable at run time.
        auto lim = law_f.limits();
        lim.input.upper(0) = 0.25f;
        REQUIRE(law_f.set_limits(lim));
        CHECK(law_f.step(10.0f, 0.0f) == 0.25f);
        // law_f = ADRCController<2, float>{kAxis}; // ill-formed: nothing replaces the type's gains (Law::operator=)
        //! The Law's state carries over to a Runtime, which retunes.
        Law<kAxis, float>::Runtime copied{law_f};
        CHECK(copied.valid());
        CHECK(copied.applied() == law_f.applied());
        CHECK(copied.observer()[2] == law_f.observer()[2]);
        CHECK(copied.limits().input.upper(0) == 0.25f);
        for (int k = 0; k < 50; ++k) {
            CHECK(copied.step(0.1f, 0.05f) == doctest::Approx(law_f.step(0.1f, 0.05f)).epsilon(1e-5));
        }
        auto g = copied.gains();
        g.wc = 12.0f;
        CHECK(copied.set_gains(g));
    }
}
