// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>
#include <numbers>

#include "damp/backend.hpp"
#include "damp/controllers/repetitive.hpp"
#include "damp/math/math.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {
constexpr double pi = std::numbers::pi;
} // namespace

TEST_SUITE("repetitive") {
    TEST_CASE("repetitive computes the period and validates the spec") {
        SUBCASE("valid: fs=1000, f0=50 -> N=20") {
            const auto r = design::repetitive(1000.0, 50.0, 1.0, 0.98, 1);
            CHECK(r.success);
            CHECK(r.config.period == 20);
            CHECK(r.config.gain == doctest::Approx(1.0));
            CHECK(r.config.q_filter == doctest::Approx(0.98));
            CHECK(r.config.lead == 1);
        }
        SUBCASE("rounds to nearest period") {
            CHECK(design::repetitive(1000.0, 49.0).config.period == 20); // 1000/49 = 20.4
            CHECK(design::repetitive(1000.0, 47.0).config.period == 21); // 1000/47 = 21.3
        }
        SUBCASE("rejects bad specs") {
            CHECK_FALSE(design::repetitive(1000.0, 1000.0).success);         // f0 >= fs
            CHECK_FALSE(design::repetitive(1000.0, -50.0).success);          // f0 < 0
            CHECK_FALSE(design::repetitive(1000.0, 50.0, 3.0).success);      // gain > 2
            CHECK_FALSE(design::repetitive(1000.0, 50.0, 1.0, 1.5).success); // q > 1
        }
    }

    TEST_CASE("drives a periodic tracking error to zero (internal-model principle)") {
        // Plant y[k] = u[k-1] (unit transport delay); repetitive controller alone.
        // With a perfect 1-sample model (lead m = 1, Q = 1, k_rc = 1) the periodic
        // tracking error must vanish after one period.
        constexpr size_t N = 20;
        const auto       design = design::repetitive(1000.0, 50.0, 1.0, 1.0, 1);
        REQUIRE(design.success);
        REQUIRE(design.config.period == N);

        RepetitiveController<64, double> rc(design);

        auto reference = [](size_t k) { return std::sin(2.0 * pi * static_cast<double>(k) / N); };

        double       u_prev = 0.0;
        double       worst_first_period = 0.0;
        double       worst_last_period = 0.0;
        const size_t periods = 8;
        for (size_t k = 0; k < periods * N; ++k) {
            const double y = u_prev;           // y[k] = u[k-1]
            const double e = reference(k) - y; // tracking error
            const double u = rc.step(e);
            u_prev = u;

            if (k < N) {
                worst_first_period = std::max(worst_first_period, std::abs(e));
            }
            if (k >= (periods - 1) * N) {
                worst_last_period = std::max(worst_last_period, std::abs(e));
            }
        }
        // The first period still has error (model not yet learned); the last is ~0.
        CHECK(worst_first_period > 0.1);
        CHECK(worst_last_period < 1e-9);
    }

    TEST_CASE("rejects a periodic disturbance") {
        // Plant y[k] = u[k-1] + d[k], periodic disturbance d, reference r = 0.
        constexpr size_t N = 24;
        const auto       design = design::repetitive(1200.0, 50.0, 1.0, 1.0, 1);
        REQUIRE(design.success);
        REQUIRE(design.config.period == N);

        RepetitiveController<64, double> rc(design);
        auto                             disturbance = [](size_t k) {
            return 0.5 * std::sin(2.0 * pi * static_cast<double>(k) / N)
                 + 0.2 * std::sin(2.0 * pi * 3.0 * static_cast<double>(k) / N); // fundamental + 3rd
        };

        double       u_prev = 0.0;
        double       worst_last_period = 0.0;
        const size_t periods = 10;
        for (size_t k = 0; k < periods * N; ++k) {
            const double y = u_prev + disturbance(k);
            const double e = 0.0 - y;
            const double u = rc.step(e);
            u_prev = u;
            if (k >= (periods - 1) * N) {
                worst_last_period = std::max(worst_last_period, std::abs(y));
            }
        }
        CHECK(worst_last_period < 1e-9); // both the fundamental and the 3rd are rejected
    }

    TEST_CASE("Q < 1 stays bounded (robustness roll-off) and reset clears state") {
        const auto                       design = design::repetitive(1000.0, 50.0, 1.0, 0.9, 1);
        RepetitiveController<64, double> rc(design);
        // Feed a unit-impulse error; with Q<1 the internal model decays, output bounded.
        double maxout = std::abs(rc.step(1.0));
        for (int k = 0; k < 2000; ++k) {
            maxout = std::max(maxout, std::abs(rc.step(0.0)));
        }
        CHECK(maxout < 2.0); // bounded (no marginal-stability blow-up)
        CHECK(maxout > 0.0); // and it actually produced corrections

        rc.reset();
        // After reset, a zero error gives zero output (empty internal model).
        CHECK(rc.step(0.0) == doctest::Approx(0.0));
    }

    TEST_CASE("invalid design yields a pass-through (zero correction)") {
        RepetitiveController<64, double> rc(design::repetitive(1000.0, 1000.0)); // invalid
        CHECK_FALSE(rc.valid());
        CHECK(rc.step(1.0) == doctest::Approx(0.0));
    }

    TEST_CASE("repetitive controller is constexpr-evaluable") {
        constexpr bool ok = []() consteval {
            auto                             d = design::repetitive(1000.0, 100.0, 1.0, 1.0, 0);
            RepetitiveController<32, double> rc(d);
            // Period N=10. Feed a constant error of 1: the first period emits no
            // correction (empty internal model), then the model from one period
            // ago surfaces, so the correction at step N (the 11th) equals 1.
            double last = 0.0;
            for (size_t k = 0; k < 11; ++k) {
                last = rc.step(1.0);
            }
            return d.success && d.config.period == 10 && damp::abs(last - 1.0) < 1e-9;
        }();
        static_assert(ok, "repetitive controller must work at compile time");
        CHECK(ok);
    }

    TEST_CASE("cross-precision converting ctor preserves the learned model") {
        const auto                       d = design::repetitive(1000.0, 50.0, 1.0, 0.99, 1);
        RepetitiveController<64, double> rc(d);

        // Learn from a periodic error for several periods.
        for (size_t k = 0; k < 100; ++k) {
            (void)rc.step(std::sin(2.0 * pi * static_cast<double>(k) / 20.0));
        }

        RepetitiveController<64, float> rf(rc); // double -> float converting ctor
        CHECK(rf.valid() == rc.valid());
        CHECK(rf.config().period == rc.config().period);
        CHECK(rf.config().gain == doctest::Approx(static_cast<float>(rc.config().gain)));

        // The learned buffer carried over: the next correction matches within
        // float precision (a freshly-built rf would instead emit a stale value).
        const double e_next = std::sin(2.0 * pi * 100.0 / 20.0);
        const double u_d = rc.step(e_next);
        const float  u_f = rf.step(static_cast<float>(e_next));
        CHECK(static_cast<double>(u_f) == doctest::Approx(u_d).epsilon(1e-4));
    }
}

// ---------------------------------------------------------------------------
// Zero-phase FIR Q-filter (binomial robustness filter)
// ---------------------------------------------------------------------------

namespace {
// Closed loop with a repetitive plug-in around a plant with a pure transport
// delay of `delay` samples and gain `g`: y[k] = g·u[k-delay] + d[k]. Reference 0,
// so e = -y. The repetitive lead is 0, so the delay is *uncompensated* — the high
// harmonics see a wrong phase. Returns the worst |y| over the run; a marginal
// (scalar Q=1) loop blows up there, a rolled-off (FIR Q) loop stays bounded.
template<typename Rc>
double closed_loop_peak(Rc& rc, int delay, double g, int steps) {
    double                 y = 0.0;
    double                 worst = 0.0;
    damp::array<double, 8> ud{}; // delay line of applied inputs
    size_t                 head = 0;
    for (int k = 0; k < steps; ++k) {
        const double d = (k == 0) ? 1.0 : 0.0; // impulse kick to excite the loop
        const double e = -y;
        const double u = rc.step(e);
        ud[head] = u;
        const double u_delayed = ud[(head + ud.size() - static_cast<size_t>(delay)) % ud.size()];
        head = (head + 1) % ud.size();
        y = g * u_delayed + d;
        worst = std::max(worst, std::abs(y));
    }
    return worst;
}
} // namespace

TEST_SUITE("repetitive FIR Q") {
    TEST_CASE("binomial Q taps are correct and unity-DC-gain") {
        SUBCASE("M=1 -> [1,2,1]/4") {
            const auto r = design::repetitive_binomial<1>(1000.0, 50.0);
            REQUIRE(r.success);
            CHECK(r.config.q_half == 1);
            CHECK(r.config.q_filter == doctest::Approx(0.5));    // 2/4
            CHECK(r.config.q_side[0] == doctest::Approx(0.25));  // 1/4
            CHECK(r.config.q_dc_gain() == doctest::Approx(1.0)); // 0.5 + 2*0.25
        }
        SUBCASE("M=2 -> [1,4,6,4,1]/16") {
            const auto r = design::repetitive_binomial<2>(1000.0, 50.0);
            REQUIRE(r.success);
            CHECK(r.config.q_filter == doctest::Approx(6.0 / 16.0));
            CHECK(r.config.q_side[0] == doctest::Approx(4.0 / 16.0));
            CHECK(r.config.q_side[1] == doctest::Approx(1.0 / 16.0));
            CHECK(r.config.q_dc_gain() == doctest::Approx(1.0));
        }
    }

    TEST_CASE("FIR-Q repetitive still rejects a multi-harmonic disturbance") {
        // Same delay-1 plant as the scalar tests, with a lead of 1 to match it.
        // The binomial Q has near-unity gain on the low harmonics (strong
        // rejection) but deliberately rolls off the higher ones (robustness), so
        // rejection is strong but not perfect — unlike the scalar Q=1 case.
        constexpr size_t N = 24;
        const auto       art = design::repetitive_binomial<1>(1200.0, 50.0, 1.0, 1);
        REQUIRE(art.success);
        REQUIRE(art.config.period == N);
        RepetitiveController<64, double, 1> rc(art);
        REQUIRE(rc.valid());

        auto disturbance = [&](size_t k) {
            return 0.5 * std::sin(2.0 * pi * static_cast<double>(k) / N)
                 + 0.2 * std::sin(2.0 * pi * 3.0 * static_cast<double>(k) / N);
        };
        double       u_prev = 0.0;
        double       worst_last = 0.0;
        const size_t periods = 40;
        for (size_t k = 0; k < periods * N; ++k) {
            const double y = u_prev + disturbance(k);
            const double u = rc.step(0.0 - y);
            u_prev = u;
            if (k >= (periods - 1) * N) {
                worst_last = std::max(worst_last, std::abs(y));
            }
        }
        CHECK(worst_last < 0.05); // both harmonics strongly attenuated (0.7 -> ~0.036, ~19x)
    }

    TEST_CASE("FIR Q stabilizes a loop where scalar Q=1 diverges") {
        // Uncompensated 1-sample plant delay (lead m=0): scalar Q=1 has loop gain
        // |1 − k_rc·z^-1| → up to 1.5 near Nyquist, so the marginal model amplifies
        // the high harmonics without bound; the binomial M=2 roll-off pulls the
        // stability term |Q − k_rc·z^-1| under 1 everywhere (peak ≈ 0.56).
        const int    delay = 1;
        const double g = 1.0;
        const int    steps = 4000;

        RepetitiveController<64, double, 0> scalar_q(design::repetitive(1000.0, 50.0, 0.5, 1.0, 0));
        RepetitiveController<64, double, 2> fir_q(design::repetitive_binomial<2>(1000.0, 50.0, 0.5, 0));
        REQUIRE(scalar_q.valid());
        REQUIRE(fir_q.valid());

        const double scalar_peak = closed_loop_peak(scalar_q, delay, g, steps);
        const double fir_peak = closed_loop_peak(fir_q, delay, g, steps);

        CHECK(scalar_peak > 10.0); // Q=1 blows up on the high harmonics
        CHECK(fir_peak < 5.0);     // the FIR roll-off keeps it bounded
        CHECK(fir_peak < scalar_peak);
    }

    TEST_CASE("as<float>() carries the FIR taps") {
        const auto art = design::repetitive_binomial<2>(1000.0, 50.0);
        const auto af = art.template as<float>();
        CHECK(af.success);
        CHECK(af.config.q_half == 2);
        CHECK(af.config.q_filter == doctest::Approx(6.0f / 16.0f));
        CHECK(af.config.q_side[1] == doctest::Approx(1.0f / 16.0f));
        RepetitiveController<64, float, 2> rc(af);
        CHECK(rc.valid());
    }

    TEST_CASE("FIR-Q repetitive is constexpr-evaluable") {
        constexpr bool ok = []() consteval {
            auto                                art = design::repetitive_binomial<1>(1000.0, 100.0, 1.0, 0);
            RepetitiveController<32, double, 1> rc(art);
            double                              last = 0.0;
            for (int k = 0; k < 40; ++k) {
                last = rc.step(1.0);
            }
            return art.success && rc.valid() && damp::abs(last) < 1e6; // bounded
        }();
        static_assert(ok, "FIR-Q repetitive must work at compile time");
        CHECK(ok);
    }
}
