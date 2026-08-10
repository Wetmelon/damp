// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>

#include "damp/backend.hpp"
#include "damp/controllers/stsmc.hpp"
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
        const double u = c.control_sliding(s, Ts);
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
            const auto r = design::stsmc(4.0);
            CHECK(r.success);
            CHECK(r.k1 == doctest::Approx(1.5 * std::sqrt(4.0)));
            CHECK(r.k2 == doctest::Approx(1.1 * 4.0));
        }
        SUBCASE("gain_margin scales both gains") {
            const auto r = design::stsmc(4.0, 0.0, 0.0, 0.0, 2.0);
            CHECK(r.k1 == doctest::Approx(2.0 * 1.5 * std::sqrt(4.0)));
            CHECK(r.k2 == doctest::Approx(2.0 * 1.1 * 4.0));
        }
        SUBCASE("rejects bad specs") {
            CHECK_FALSE(design::stsmc(0.0).success);                     // L <= 0
            CHECK_FALSE(design::stsmc(4.0, 0.0, -1.0).success);          // kl < 0
            CHECK_FALSE(design::stsmc(4.0, 0.0, 0.0, 0.0, 0.5).success); // margin < 1
        }
        SUBCASE("direct-gain factory validates") {
            CHECK(design::stsmc_gains(2.0, 2.2).success);
            CHECK_FALSE(design::stsmc_gains(-1.0, 2.2).success);
            CHECK_FALSE(design::stsmc_gains(2.0, 0.0).success);
        }
    }

    TEST_CASE("classic super-twisting rejects a Lipschitz disturbance, continuously") {
        // d(t) = sin(2t) -> |d_dot| <= 2 = L. 1 kHz loop.
        const double Ts = 1e-3;
        const double L = 2.0;
        const auto   art = design::stsmc(L);
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

        STSMCController<double> st(design::stsmc(L));
        const auto [st_s, st_du] = run(st, d, Ts, 8000, 1.0);

        // A bare sign()-based controller of comparable authority: u = -k*sign(s).
        struct SignSMC {
            double k;
            double control_sliding(double s, double /*Ts*/) { return -k * static_cast<double>(damp::sgn(s)); }
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
        const auto   art = design::stsmc(L, 0.0, 1.5);
        REQUIRE(art.success);
        REQUIRE(art.k_lin == doctest::Approx(1.5));
        STSMCController<double> c(art);

        auto d = [](double t) { return 0.5 * std::sin(3.0 * t); };
        const auto [worst_s, worst_du] = run(c, d, Ts, 8000, 1.0);
        CHECK(worst_s < 1e-2);
        CHECK(worst_du < 0.1);
    }

    TEST_CASE("boundary layer keeps the canonical loop bounded and continuous") {
        const double Ts = 1e-3;
        const auto   art = design::stsmc(2.0, 0.0, 0.0, 0.05);
        REQUIRE(art.success);
        STSMCController<double> c(art);

        auto d = [](double t) { return std::sin(2.0 * t); };
        const auto [worst_s, worst_du] = run(c, d, Ts, 8000, 1.0);
        CHECK(worst_s < 0.25);  // softened sign -> wider band than true sign (which hit <1e-2)
        CHECK(worst_du < 0.05); // still continuous
    }

    TEST_CASE("(r,y) overload builds s = lambda*e + e_dot then applies the STA") {
        // The (r, y) form is a surface-builder over control(s). The first call
        // seeds e_prev = e (no derivative kick), so e_dot = 0 and s = lambda*e;
        // the next call exercises the backward difference.
        const double            Ts = 1e-3;
        const double            lambda = 5.0;
        STSMCController<double> via_ry(design::stsmc(2.0, lambda));
        STSMCController<double> via_s(design::stsmc(2.0, lambda));

        const double e1 = 1.0;         // r=1, y=0
        const double s1 = lambda * e1; // first call: e_dot = 0 (seeded)
        CHECK(via_ry.control(1.0, 0.0, Ts) == doctest::Approx(via_s.control_sliding(s1, Ts)));

        const double e2 = 1.5 - 0.2;                        // r=1.5, y=0.2
        const double s2 = (lambda * e2) + ((e2 - e1) / Ts); // e_dot from backward difference
        CHECK(via_ry.control(1.5, 0.2, Ts) == doctest::Approx(via_s.control_sliding(s2, Ts)));
    }

    TEST_CASE("cross-precision converting ctor preserves gains and integral state") {
        STSMCController<double> d(design::stsmc(2.0, 5.0));
        (void)d.control_sliding(1.0, 1e-3);
        (void)d.control_sliding(0.8, 1e-3); // advance the integral state v
        STSMCController<float> f(d);        // converting ctor (private members via friend)
        CHECK(f.valid());
        CHECK(f.integral_state() == doctest::Approx(static_cast<float>(d.integral_state())));
        // Continues from the copied state: next command matches within float precision.
        CHECK(f.control_sliding(0.6f, 1e-3f) == doctest::Approx(static_cast<float>(d.control_sliding(0.6, 1e-3))).epsilon(1e-4));
    }

    TEST_CASE("non-positive Ts holds: no divide in (r,y), no backward integration in (s)") {
        STSMCController<double> c(design::stsmc(2.0, 5.0));
        CHECK(c.control(1.0, 0.0, 0.0) == doctest::Approx(0.0)); // (r,y) Ts=0 -> inert, no NaN
        (void)c.control_sliding(0.5, 1e-3);
        const double v = c.integral_state();
        (void)c.control_sliding(0.5, 0.0); // (s) Ts=0 -> integral frozen
        CHECK(c.integral_state() == doctest::Approx(v));
    }

    TEST_CASE("invalid design is inert; reset clears state") {
        STSMCController<double> bad(design::stsmc(-1.0));
        CHECK_FALSE(bad.valid());
        CHECK(bad.control_sliding(5.0, 1e-3) == doctest::Approx(0.0));

        STSMCController<double> c(design::stsmc(2.0));
        (void)c.control_sliding(1.0, 1e-3);
        (void)c.control_sliding(1.0, 1e-3);
        CHECK(c.integral_state() != 0.0);
        c.reset();
        CHECK(c.integral_state() == doctest::Approx(0.0));
    }

    TEST_CASE("float deployment via as<float>()") {
        const auto             art = design::stsmc(2.0);
        STSMCController<float> c(art.as<float>());
        const float            u = c.control_sliding(1.0f, 1e-3f);
        CHECK(std::isfinite(u));
        CHECK(u < 0.0f); // s > 0 -> u pushes negative
    }

    TEST_CASE("super-twisting controller is constexpr-evaluable") {
        constexpr double final_s = []() consteval {
            auto                    art = design::stsmc(2.0);
            STSMCController<double> c(art);
            double                  s = 1.0;
            for (int k = 0; k < 4000; ++k) {
                const double u = c.control_sliding(s, 1e-3);
                s += 1e-3 * (u + 0.5 * damp::sin(2.0 * k * 1e-3));
            }
            return s;
        }();
        static_assert(final_s < 0.05 && final_s > -0.05, "STA must converge at compile time");
        CHECK(final_s == doctest::Approx(0.0).epsilon(0.05));
    }
}
