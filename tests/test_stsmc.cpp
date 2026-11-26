// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <limits>

#include "damp/backend.hpp"
#include "damp/controllers/stsmc.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

// Closed-loop sim of the super-twisting algorithm on the canonical relative-
// degree-1 plant  ṡ = u + d(t),  integrated with explicit Euler at step Ts.
// Returns {worst |s| over the final `tail` fraction, worst |Δu| between
// consecutive samples over that tail}. The second number is the chattering
// metric: a *continuous* controller keeps it small; a sign() controller jumps
// by ~2k every crossing.
template<typename Ctrl, typename Dist>
damp::pair<double, double> run(Ctrl& c, Dist d, double Ts, int steps, double s0, double tail = 0.25) {
    double    s = s0;
    double    u_prev = 0.0;
    double    worst_s = 0.0;
    double    worst_du = 0.0;
    const int tail_start = static_cast<int>(steps * (1.0 - tail));
    for (int k = 0; k < steps; ++k) {
        const double u = c.step_sliding(s);
        if (k >= tail_start) {
            worst_s = std::max(worst_s, std::abs(s));
            worst_du = std::max(worst_du, std::abs(u - u_prev));
        }
        u_prev = u;
        s += Ts * (u + d(k * Ts)); // plant: ṡ = u + d
    }
    return {worst_s, worst_du};
}

} // namespace

TEST_SUITE("stsmc") {
    TEST_CASE("stsmc computes Levant/Moreno gains and validates") {
        SUBCASE("valid: k1 = 1.5*sqrt(L), k2 = 1.1*L") {
            const auto r = design::stsmc({.disturbance_bound = 4.0}, 1e-3);
            CHECK(r.success);
            CHECK(r.gains.k1 == doctest::Approx(1.5 * std::sqrt(4.0)));
            CHECK(r.gains.k2 == doctest::Approx(1.1 * 4.0));
        }
        SUBCASE("gain_margin scales both gains") {
            const auto r = design::stsmc({.disturbance_bound = 4.0, .gain_margin = 2.0}, 1e-3);
            CHECK(r.gains.k1 == doctest::Approx(2.0 * 1.5 * std::sqrt(4.0)));
            CHECK(r.gains.k2 == doctest::Approx(2.0 * 1.1 * 4.0));
        }
        SUBCASE("rejects bad specs") {
            CHECK_FALSE(design::stsmc({.disturbance_bound = 0.0}, 1e-3).success);                     // L <= 0
            CHECK_FALSE(design::stsmc({.disturbance_bound = 4.0, .k_lin = -1.0}, 1e-3).success);      // kl < 0
            CHECK_FALSE(design::stsmc({.disturbance_bound = 4.0, .gain_margin = 0.5}, 1e-3).success); // margin < 1
        }
        SUBCASE("the direct-gain form validates") {
            CHECK(design::stsmc({.k1 = 2.0, .k2 = 2.2}, 1e-3).success);
            CHECK_FALSE(design::stsmc({.k1 = -1.0, .k2 = 2.2}, 1e-3).success);
            CHECK_FALSE(design::stsmc({.k1 = 2.0, .k2 = 0.0}, 1e-3).success);
        }
    }

    TEST_CASE("classic super-twisting rejects a Lipschitz disturbance, continuously") {
        // d(t) = sin(2t) -> |d_dot| <= 2 = L. 1 kHz loop.
        const double Ts = 1e-3;
        const double L = 2.0;
        const auto   art = design::stsmc({.disturbance_bound = L}, 1e-3);
        REQUIRE(art.success);
        STSMCController<double> c(art);

        auto d = [](double t) { return std::sin(2.0 * t); };
        const auto [worst_s, worst_du] = run(c, d, Ts, 8000, 1.0); // 8 s

        CHECK(worst_s < 1e-2);  // s driven into a tight band despite the disturbance
        CHECK(worst_du < 0.05); // control is continuous (no sign()-style jumps)
        // The integral state tracks -d: at the end of an 8 s run it is near -sin(2*8).
        CHECK(c.integral_state() == doctest::Approx(-std::sin(2.0 * 8.0)).epsilon(0.1));
    }

    TEST_CASE("first-order SMC chatters where super-twisting does not (the whole point)") {
        // Same plant/disturbance; compare the steady-state control jumpiness.
        const double Ts = 1e-3;
        const double L = 2.0;
        auto         d = [](double t) { return std::sin(2.0 * t); };

        STSMCController<double> st(design::stsmc({.disturbance_bound = L}, 1e-3));
        const auto [st_s, st_du] = run(st, d, Ts, 8000, 1.0);

        // A bare sign()-based controller of comparable authority: u = -k*sign(s).
        struct SignSMC {
            double k;
            double step_sliding(double s) { return -k * static_cast<double>(damp::sgn(s)); }
        } sign_smc{3.0};
        const auto [sg_s, sg_du] = run(sign_smc, d, Ts, 8000, 1.0);

        // Super-twisting's steady-state control steps are far smaller than the
        // sign controller's ~2k bang-bang jumps.
        CHECK(st_du < 0.1);
        CHECK(sg_du > 1.0);
        CHECK(st_du < sg_du);
    }

    TEST_CASE("generalized STA (kl > 0) also converges") {
        const double Ts = 1e-3;
        const double L = 2.0;
        const auto   art = design::stsmc({.disturbance_bound = L, .k_lin = 1.5}, 1e-3);
        REQUIRE(art.success);
        REQUIRE(art.gains.k_lin == doctest::Approx(1.5));
        STSMCController<double> c(art);

        auto d = [](double t) { return 0.5 * std::sin(3.0 * t); };
        const auto [worst_s, worst_du] = run(c, d, Ts, 8000, 1.0);
        CHECK(worst_s < 1e-2);
        CHECK(worst_du < 0.1);
    }

    TEST_CASE("boundary layer keeps the canonical loop bounded and continuous") {
        const double Ts = 1e-3;
        const auto   art = design::stsmc({.disturbance_bound = 2.0, .phi = 0.05}, 1e-3);
        REQUIRE(art.success);
        STSMCController<double> c(art);

        auto d = [](double t) { return std::sin(2.0 * t); };
        const auto [worst_s, worst_du] = run(c, d, Ts, 8000, 1.0);
        CHECK(worst_s < 0.25);  // softened sign -> wider band than true sign (which hit <1e-2)
        CHECK(worst_du < 0.05); // still continuous
    }

    TEST_CASE("(r,y) overload builds sigma = -(lambda*e + e_dot) then applies the STA") {
        // The (r, y) form is a surface-builder over control(s). The first call
        // seeds e_prev = e (no derivative kick), so e_dot = 0 and s = lambda*e;
        // the next call exercises the backward difference.
        const double            Ts = 1e-3;
        const double            lambda = 5.0;
        STSMCController<double> via_ry(design::stsmc({.disturbance_bound = 2.0, .lambda = lambda}, 1e-3));
        STSMCController<double> via_s(design::stsmc({.disturbance_bound = 2.0, .lambda = lambda}, 1e-3));

        const double e1 = 1.0;         // r=1, y=0
        const double s1 = lambda * e1; // first call: e_dot = 0 (seeded)
        CHECK(via_ry.step(1.0, 0.0) == doctest::Approx(via_s.step_sliding(-s1)));

        const double e2 = 1.5 - 0.2;                        // r=1.5, y=0.2
        const double s2 = (lambda * e2) + ((e2 - e1) / Ts); // e_dot from backward difference
        CHECK(via_ry.step(1.5, 0.2) == doctest::Approx(via_s.step_sliding(-s2)));
    }

    TEST_CASE("cross-precision converting ctor preserves gains and integral state") {
        STSMCController<double> d(design::stsmc({.disturbance_bound = 2.0, .lambda = 5.0}, 1e-3));
        (void)d.step_sliding(1.0);
        (void)d.step_sliding(0.8);   // advance the integral state v
        STSMCController<float> f(d); // converting ctor (private members via friend)
        CHECK(f.valid());
        CHECK(f.integral_state() == doctest::Approx(static_cast<float>(d.integral_state())));
        // Continues from the copied state: next command matches within float precision.
        CHECK(f.step_sliding(0.6f) == doctest::Approx(static_cast<float>(d.step_sliding(0.6))).epsilon(1e-4));
    }

    TEST_CASE("a Ts that is not finite and positive is refused, and leaves the runtime invalid: both entry points emit 0") {
        //! Ts = 0 used to be the continuous design, for analysis; the runtime refused it. Every Ts
        //! that is not finite and positive is refused by the design now.
        CHECK(design::stsmc({.disturbance_bound = 2.0, .lambda = 5.0}, 0.0).status == design::DesignStatus::invalid_sample_time);
        STSMCController<double> c(design::stsmc({.disturbance_bound = 2.0, .lambda = 5.0}, 0.0));
        CHECK_FALSE(c.valid());
        CHECK(c.step(1.0, 0.0) == 0.0);
        CHECK(c.step_sliding(-5.0) == 0.0);
        CHECK(c.integral_state() == 0.0);
        CHECK(design::stsmc({.disturbance_bound = 2.0, .lambda = 5.0}, -1e-3).status == design::DesignStatus::invalid_sample_time);
        CHECK(design::stsmc({.k1 = 2.0, .k2 = 2.2}, -1e-3).status == design::DesignStatus::invalid_sample_time);
        CHECK(design::stsmc({.k1 = 2.0, .k2 = 2.2}, std::numeric_limits<double>::infinity()).status == design::DesignStatus::invalid_sample_time);
    }

    TEST_CASE("invalid design is inert; reset clears state") {
        STSMCController<double> bad(design::stsmc({.disturbance_bound = -1.0}, 1e-3));
        CHECK_FALSE(bad.valid());
        CHECK(bad.step_sliding(5.0) == doctest::Approx(0.0));

        STSMCController<double> c(design::stsmc({.disturbance_bound = 2.0}, 1e-3));
        (void)c.step_sliding(1.0);
        (void)c.step_sliding(1.0);
        CHECK(c.integral_state() != 0.0);
        c.reset();
        CHECK(c.integral_state() == doctest::Approx(0.0));
    }

    TEST_CASE("float deployment via as<float>()") {
        const auto             art = design::stsmc({.disturbance_bound = 2.0}, 1e-3);
        STSMCController<float> c(art.as<float>());
        const float            u = c.step_sliding(1.0f);
        CHECK(std::isfinite(u));
        CHECK(u < 0.0f); // s > 0 -> u pushes negative
    }

    TEST_CASE("super-twisting controller is constexpr-evaluable") {
        constexpr double final_s = []() consteval {
            auto                    art = design::stsmc({.disturbance_bound = 2.0}, 1e-3);
            STSMCController<double> c(art);
            double                  s = 1.0;
            for (int k = 0; k < 4000; ++k) {
                const double u = c.step_sliding(s);
                s += 1e-3 * (u + 0.5 * damp::sin(2.0 * k * 1e-3));
            }
            return s;
        }();
        static_assert(final_s < 0.05 && final_s > -0.05, "STA must converge at compile time");
        CHECK(final_s == doctest::Approx(0.0).epsilon(0.05));
    }

    TEST_CASE("STSMC (r, y) regulates a unit-gain double integrator under a Lipschitz load") {
        const double            Ts = 1e-3;
        STSMCController<double> ctrl(design::stsmc({.disturbance_bound = 2.0, .lambda = 5.0}, Ts));
        double                  y = 0.0;
        double                  v = 0.0;
        for (int k = 0; k < 10000; ++k) {
            const double u = ctrl.step(1.0, y);
            v += Ts * (u + (0.5 * std::sin(k * Ts)));
            y += Ts * v;
        }
        CHECK(std::abs(y - 1.0) < 0.02);
    }

    TEST_CASE("the smoothed sign reaches ±1, so |s| settles within about φ; NaN and λ < 0 are refused") {
        // ṡ = u + L·t with |ḋ| = L = 2 at the design bound k₂ = 1.1·L. With s/(|s| + φ) the
        // integral's gain reached L only at |s| = 10φ; smooth_sign is exact outside |s| < φ.
        constexpr double phi = 0.05;
        constexpr double Ts = 1e-4;
        const auto       art = design::stsmc({.disturbance_bound = 2.0, .phi = phi}, Ts);
        REQUIRE(art.success);
        STSMCController<double> c(art);
        double                  s = 0.0;
        double                  worst = 0.0;
        for (int k = 0; k < 200000; ++k) {
            const double t = static_cast<double>(k) * Ts;
            s += Ts * (c.step_sliding(s) + (2.0 * t));
            if (k >= 150000) {
                worst = std::max(worst, std::abs(s));
            }
        }
        CHECK(worst < phi);

        const double nan = std::numeric_limits<double>::quiet_NaN();
        CHECK(design::stsmc({.disturbance_bound = nan}, 1e-3).status == design::DesignStatus::invalid_parameter);
        CHECK(design::stsmc({.disturbance_bound = 2.0, .gain_margin = nan}, 1e-3).status == design::DesignStatus::invalid_parameter);
        CHECK(design::stsmc({.disturbance_bound = 2.0, .lambda = -5.0}, 1e-3).status == design::DesignStatus::invalid_parameter);
        CHECK(design::stsmc({.disturbance_bound = 2.0, .lambda = nan}, 1e-3).status == design::DesignStatus::invalid_parameter);
        CHECK(design::stsmc({.k1 = 2.0, .k2 = 2.2, .lambda = nan}, 1e-3).status == design::DesignStatus::invalid_parameter);
        CHECK(design::stsmc({.k1 = 2.0, .k2 = 2.2, .lambda = -1.0}, 1e-3).status == design::DesignStatus::invalid_parameter);
    }

    TEST_CASE_TEMPLATE("InputLimits: the integral state holds against the clamp, the law's own or one reported", T, float, double) {
        //! v ran on while the clamped command could not bring s to 0, and had to be run back
        //! before the law answered again.
        constexpr double   Ts = 1e-3;
        const auto         res = design::stsmc({.disturbance_bound = 2.0}, Ts, design::InputLimits{.lower = {-1.0}, .upper = {1.0}});
        const auto         open = design::stsmc({.disturbance_bound = 2.0}, Ts);
        STSMCController<T> bounded{res.template as<T>()};
        STSMCController<T> outside{open.template as<T>()};
        T                  sa{0.5};
        T                  sb{0.5};
        T                  v_on_rail{0};
        for (int k = 0; k < 14000; ++k) {
            const T d = (k < 2000) ? T{3} : T{0.5}; // ṡ = u + d: a load past the ±1 authority, then within it
            const T ua = bounded.step_sliding(sa);
            const T ub = damp::clamp(outside.step_sliding(sb), T{-1}, T{1});
            outside.back_calculate(ub);
            CHECK(damp::abs(ua) <= T{1});
            CHECK(damp::abs(ua - ub) <= T{64} * eps<T>()); // one tick, whichever stage clamps
            sa += static_cast<T>(Ts) * (ua + d);
            sb += static_cast<T>(Ts) * (ub + d);
            if (k == 1999) {
                v_on_rail = bounded.integral_state();
            }
        }
        CHECK(damp::abs(v_on_rail) <= static_cast<T>(0.01)); // held, where it would have run to −4.4
        CHECK(damp::abs(sa) < static_cast<T>(0.05));
        CHECK(static_cast<double>(bounded.integral_state()) == doctest::Approx(-0.5).epsilon(0.1)); // v settles on −d
        CHECK(bounded.last_command() == bounded.limits().input.saturate(bounded.last_command()));
    }

    TEST_CASE("limits that are not a box are refused by the design, and leave a runtime built by hand invalid") {
        CHECK(design::stsmc({.disturbance_bound = 2.0}, 1e-3, design::InputLimits{.lower = {1.0}, .upper = {-1.0}}).status == design::DesignStatus::invalid_parameter);
        auto res = design::stsmc({.disturbance_bound = 2.0}, 1e-3);
        res.limits.input = {.lower = 1.0, .upper = -1.0};
        STSMCController<double> inert{res};
        CHECK_FALSE(inert.valid());
        CHECK(inert.step_sliding(1.0) == 0.0);
        CHECK(res.as<float>().limits.input.lower(0) == 1.0F); // the box travels with the result
    }

    TEST_CASE("track preloads v: the next step at the same r and y, or s, returns u") {
        //! Without a handover the law started from its last v and the first command jumped.
        const auto              res = design::stsmc({.disturbance_bound = 2.0, .lambda = 5.0}, 1e-3, design::InputLimits{.lower = {-3.0}, .upper = {3.0}});
        STSMCController<double> law{res};
        for (int k = 0; k < 50; ++k) {
            (void)law.step(1.0, 0.01 * static_cast<double>(k));
        }
        double y = 0.4;
        for (int k = 0; k < 30; ++k) {
            y = 0.4 + (0.001 * static_cast<double>(k));
            law.track(1.7, 0.8, y); // an operator holds the plant
        }
        CHECK(law.last_command() == 1.7);
        CHECK(law.step(0.8, y) == doctest::Approx(1.7).epsilon(1e-12));
        law.track(9.0, 0.8, y); // past the limit: lands on it
        CHECK(law.step(0.8, y) == doctest::Approx(3.0).epsilon(1e-12));

        STSMCController<double> canon{design::stsmc({.disturbance_bound = 2.0}, 1e-3)};
        canon.track_sliding(-0.6, 0.3);
        CHECK(canon.step_sliding(0.3) == doctest::Approx(-0.6).epsilon(1e-12));
    }

    TEST_CASE("set_gains keeps v and retunes; set_limits moves the clamp; a bad value is refused") {
        STSMCController<double> law{design::stsmc({.disturbance_bound = 2.0}, 1e-3)};
        for (int k = 0; k < 100; ++k) {
            (void)law.step_sliding(0.5);
        }
        const double v = law.integral_state();
        auto         g = law.gains();
        g.k1 = 2.0 * g.k1;
        REQUIRE(law.set_gains(g));
        CHECK(law.integral_state() == v);
        CHECK(law.gains().k1 == g.k1);
        auto bad = g;
        bad.k2 = 0.0;
        CHECK_FALSE(law.set_gains(bad));
        CHECK(law.gains().k2 == g.k2);
        auto lim = law.limits();
        lim.input = {.lower = -0.1, .upper = 0.1};
        REQUIRE(law.set_limits(lim));
        CHECK(damp::abs(law.step_sliding(5.0)) <= 0.1);
        lim.input = {.lower = 0.1, .upper = -0.1};
        CHECK_FALSE(law.set_limits(lim));
        CHECK(law.limits().input.upper(0) == 0.1);
    }

    TEST_CASE("a constexpr super-twisting design's Law: its gains constants of the type") {
        static constexpr auto  kSt = design::stsmc({.disturbance_bound = 2.0, .lambda = 5.0}, 1e-3);
        Law<kSt, float>        law;
        STSMCController<float> held{kSt};
        static_assert(sizeof(Law<kSt, float>) < sizeof(STSMCController<float>));
        for (int k = 0; k < 20; ++k) {
            const float y = 0.01f * static_cast<float>(k);
            CHECK(law.step(1.0f, y) == doctest::Approx(held.step(1.0f, y)).epsilon(1e-5));
        }
        law.track(0.5f, 1.0f, 0.3f);
        CHECK(law.step(1.0f, 0.3f) == doctest::Approx(0.5f).epsilon(1e-5));
    }
}
