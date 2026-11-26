// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/repetitive.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {
constexpr double pi = std::numbers::pi;

/** @brief [rad/s] The angular frequency of @p f [Hz]. */
constexpr double hz(double f) {
    return 2.0 * pi * f;
}
} // namespace

TEST_SUITE("repetitive") {
    TEST_CASE("repetitive computes the period and validates the spec") {
        SUBCASE("valid: 50 Hz at 1 kHz -> N=20") {
            const auto r = design::repetitive({.w0 = hz(50.0), .gain = 1.0, .q = 0.98, .lead = 1}, 1.0 / 1000.0);
            CHECK(r.success);
            CHECK(r.config.period == 20);
            CHECK(r.config.gain == doctest::Approx(1.0));
            CHECK(r.config.q_filter == doctest::Approx(0.98));
            CHECK(r.config.lead == 1);
        }
        SUBCASE("rounds to nearest period") {
            CHECK(design::repetitive({.w0 = hz(49.0)}, 1.0 / 1000.0).config.period == 20); // 1000/49 = 20.4
            CHECK(design::repetitive({.w0 = hz(47.0)}, 1.0 / 1000.0).config.period == 21); // 1000/47 = 21.3
        }
        SUBCASE("rejects bad specs") {
            CHECK_FALSE(design::repetitive({.w0 = hz(1000.0)}, 1.0 / 1000.0).success);                      // ω₀ at the sample rate: N = 1
            CHECK_FALSE(design::repetitive({.w0 = hz(-50.0)}, 1.0 / 1000.0).success);                       // ω₀ < 0
            CHECK_FALSE(design::repetitive({.w0 = hz(50.0), .gain = 3.0}, 1.0 / 1000.0).success);           // gain > 2
            CHECK_FALSE(design::repetitive({.w0 = hz(50.0), .gain = 1.0, .q = 1.5}, 1.0 / 1000.0).success); // q > 1
        }
    }

    TEST_CASE("drives a periodic tracking error to zero (internal-model principle)") {
        // Plant y[k] = u[k-1] (unit transport delay); repetitive controller alone.
        // With a perfect 1-sample model (lead m = 1, Q = 1, k_rc = 1) the periodic
        // tracking error must vanish after one period.
        constexpr size_t N = 20;
        const auto       design = design::repetitive({.w0 = hz(50.0), .gain = 1.0, .q = 1.0, .lead = 1}, 1.0 / 1000.0);
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
        const auto       design = design::repetitive({.w0 = hz(50.0), .gain = 1.0, .q = 1.0, .lead = 1}, 1.0 / 1200.0);
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
        const auto                       design = design::repetitive({.w0 = hz(50.0), .gain = 1.0, .q = 0.9, .lead = 1}, 1.0 / 1000.0);
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
        RepetitiveController<64, double> rc(design::repetitive({.w0 = hz(1000.0)}, 1.0 / 1000.0)); // invalid
        CHECK_FALSE(rc.valid());
        CHECK(rc.step(1.0) == doctest::Approx(0.0));
    }

    TEST_CASE("repetitive controller is constexpr-evaluable") {
        constexpr bool ok = []() consteval {
            auto                             d = design::repetitive({.w0 = hz(100.0), .gain = 1.0, .q = 1.0, .lead = 0}, 1.0 / 1000.0);
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
        const auto                       d = design::repetitive({.w0 = hz(50.0), .gain = 1.0, .q = 0.99, .lead = 1}, 1.0 / 1000.0);
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

    TEST_CASE("repetitive design rejects a fundamental whose period overflows size_t") {
        CHECK_FALSE(design::repetitive({.w0 = hz(1e-300)}, 1.0 / 1e300).success);
        CHECK(design::repetitive({.w0 = hz(50.0)}, 1.0 / 10000.0).config.period == 200);
        CHECK(design::repetitive_binomial<2>({.w0 = static_cast<float>(hz(50.0))}, 1e-4f).success); // float taps still sum to ≤ 1
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
            const auto r = design::repetitive_binomial<1>({.w0 = hz(50.0)}, 1.0 / 1000.0);
            REQUIRE(r.success);
            CHECK(r.config.q_half == 1);
            CHECK(r.config.q_filter == doctest::Approx(0.5));    // 2/4
            CHECK(r.config.q_side[0] == doctest::Approx(0.25));  // 1/4
            CHECK(r.config.q_dc_gain() == doctest::Approx(1.0)); // 0.5 + 2*0.25
        }
        SUBCASE("M=2 -> [1,4,6,4,1]/16") {
            const auto r = design::repetitive_binomial<2>({.w0 = hz(50.0)}, 1.0 / 1000.0);
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
        const auto       art = design::repetitive_binomial<1>({.w0 = hz(50.0), .gain = 1.0, .lead = 1}, 1.0 / 1200.0);
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
        // stability term |Q·(1 − k_rc·z^-1)| under 1 everywhere.
        const int    delay = 1;
        const double g = 1.0;
        const int    steps = 4000;

        RepetitiveController<64, double, 0> scalar_q(design::repetitive({.w0 = hz(50.0), .gain = 0.5, .q = 1.0, .lead = 0}, 1.0 / 1000.0));
        RepetitiveController<64, double, 2> fir_q(design::repetitive_binomial<2>({.w0 = hz(50.0), .gain = 0.5, .lead = 0}, 1.0 / 1000.0));
        REQUIRE(scalar_q.valid());
        REQUIRE(fir_q.valid());

        const double scalar_peak = closed_loop_peak(scalar_q, delay, g, steps);
        const double fir_peak = closed_loop_peak(fir_q, delay, g, steps);

        CHECK(scalar_peak > 10.0); // Q=1 blows up on the high harmonics
        CHECK(fir_peak < 5.0);     // the FIR roll-off keeps it bounded
        CHECK(fir_peak < scalar_peak);
    }

    TEST_CASE("as<float>() carries the FIR taps") {
        const auto art = design::repetitive_binomial<2>({.w0 = hz(50.0)}, 1.0 / 1000.0);
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
            auto                                art = design::repetitive_binomial<1>({.w0 = hz(100.0), .gain = 1.0, .lead = 0}, 1.0 / 1000.0);
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

    TEST_CASE("Q on the output tap stabilizes a plant resonance the memory-only tap diverges on") {
        // y[k] = −0.36·y[k−2] + 1.36·u[k−2]: DC gain 1, a peak of 2.1 at fs/4, where
        // |1 − k_rc·G| > 1. With Q on both the memory and the output, the loop condition is
        // |Q·(1 − k_rc·G)| < 1 (max 0.77 here); with Q on the memory only it was
        // |Q − k_rc·G| < 1 (max 1.34), and |e| grew past 1e21 in 40000 ticks.
        RepetitiveController<256, double, 2> rc(design::repetitive_binomial<2>({.w0 = hz(50.0), .gain = 0.5}, 1.0 / 10000.0));
        REQUIRE(rc.valid());
        double y1 = 0.0;
        double y2 = 0.0;
        double u1 = 0.0;
        double u2 = 0.0;
        double worst = 0.0;
        for (int k = 0; k < 40000; ++k) {
            const double y = (-0.36 * y2) + (1.36 * u2);
            const double e = std::sin(2.0 * std::numbers::pi * 50.0 * static_cast<double>(k) / 10000.0) - y;
            const double u = rc.step(e);
            y2 = y1;
            y1 = y;
            u2 = u1;
            u1 = u;
            if (k >= 40000 - 200) {
                worst = std::max(worst, std::abs(e));
            }
        }
        CHECK(worst < 1e-2);
    }

    TEST_CASE("the output tap's FIR window must lie in the past; NaN rates are refused") {
        // N = 20, M = 2: lead + M must stay below N.
        CHECK(design::repetitive_binomial<2>({.w0 = hz(50.0), .gain = 1.0, .lead = 17}, 1.0 / 1000.0).success);
        const auto late = design::repetitive_binomial<2>({.w0 = hz(50.0), .gain = 1.0, .lead = 18}, 1.0 / 1000.0);
        CHECK_FALSE(late.success);
        CHECK(late.status == design::DesignStatus::invalid_parameter);

        const double nan = std::numeric_limits<double>::quiet_NaN();
        CHECK(design::repetitive({.w0 = nan}, 1e-3).status == design::DesignStatus::invalid_parameter);
        CHECK(design::repetitive({.w0 = hz(50.0)}, nan).status == design::DesignStatus::invalid_sample_time);
        CHECK_FALSE(design::repetitive({.w0 = hz(50.0), .gain = nan}, 1e-3).success);
        CHECK(design::detail::repetitive_period(hz(50.0), nan) == 0);
        CHECK(design::detail::repetitive_period(nan, 1e-3) == 0);
    }

    TEST_CASE("back_calculate: the memory does not learn from samples the actuator could not follow") {
        //! With Q = 1 the memory integrates the error period over period. Behind a clamp the
        //! error at the clipped samples never closes, so the memory grew for as long as it lasted.
        constexpr size_t N = 20;
        const auto       design = design::repetitive({.w0 = hz(50.0), .gain = 1.0, .q = 1.0, .lead = 1}, 1.0 / 1000.0);
        REQUIRE(design.success);
        RepetitiveController<64, double> told(design);
        RepetitiveController<64, double> untold(design);
        const double                     rail = 0.5;
        double                           u_t = 0.0; // y[k] = u[k−1], the clamped command
        double                           u_u = 0.0;
        double                           worst_t = 0.0;
        double                           worst_u = 0.0;
        double                           told_peak = 0.0;
        double                           untold_peak = 0.0;
        for (size_t k = 0; k < 60 * N; ++k) {
            const double amplitude = (k < 50 * N) ? 1.0 : 0.4; // past the rail, then within it
            const double r = amplitude * std::sin(2.0 * pi * static_cast<double>(k) / static_cast<double>(N));
            const double e_t = r - u_t;
            const double e_u = r - u_u;
            const double c_t = told.step(e_t);
            u_t = damp::clamp(c_t, -rail, rail);
            told.back_calculate(u_t);
            const double c_u = untold.step(e_u);
            u_u = damp::clamp(c_u, -rail, rail);
            if (k >= 49 * N && k < 50 * N) {
                told_peak = damp::max(told_peak, damp::abs(c_t));
                untold_peak = damp::max(untold_peak, damp::abs(c_u));
            }
            if (k >= 59 * N) {
                worst_t = damp::max(worst_t, damp::abs(e_t));
                worst_u = damp::max(worst_u, damp::abs(e_u));
            }
        }
        CHECK(told_peak <= 1.0 + 1e-12); // what the first period learned, held
        CHECK(untold_peak > 10.0);       // half the reference more, every period
        CHECK(worst_t < 1e-9);           // learning again once the reference is within reach
        CHECK(worst_u > 0.05);           // still on the rail at the peaks
        CHECK(told.last_command() == doctest::Approx(u_t));
    }

    TEST_CASE("back_calculate with no lead: a clamp a ring ago does not stop this period's learning") {
        const auto design = design::repetitive({.w0 = hz(250.0), .gain = 1.0, .q = 1.0, .lead = 0}, 1.0 / 1000.0); // N = 4: the ring is the period
        REQUIRE(design.success);
        REQUIRE(design.config.period == 4);
        RepetitiveController<4, double> rc(design);
        static_cast<void>(rc.step(1.0));
        rc.back_calculate(-1.0); // cut off from above: tick 0's error is not learned
        for (int k = 1; k < 4; ++k) {
            static_cast<void>(rc.step(0.0));
        }
        static_cast<void>(rc.step(1.0)); // tick 4, sent as asked: learned
        for (int k = 5; k < 8; ++k) {
            static_cast<void>(rc.step(0.0));
        }
        CHECK(rc.step(0.0) == doctest::Approx(1.0)); // was 0: tick 4 read tick 0's flag
    }
}

namespace {
constexpr auto kRep = design::repetitive({.w0 = hz(50.0), .gain = 1.0, .q = 0.98, .lead = 1}, 1e-3); // N = 20
constexpr auto kBinomial = design::repetitive_binomial<2>({.w0 = hz(50.0), .gain = 0.5}, 1e-4);      // N = 200, M = 2
} // namespace

TEST_SUITE("repetitive: the runtime contract") {
    TEST_CASE("Ts is refused unless finite and positive; the period comes from Ts and ω₀") {
        // The design took the sample rate and the fundamental in Hz; it takes Ts and ω₀ [rad/s], as PR does.
        for (const double Ts : {0.0, -1e-3, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
            CAPTURE(Ts);
            CHECK(design::repetitive({.w0 = hz(50.0)}, Ts).status == design::DesignStatus::invalid_sample_time);
            CHECK(design::repetitive_binomial<1>({.w0 = hz(50.0)}, Ts).status == design::DesignStatus::invalid_sample_time);
        }
        CHECK(design::repetitive({.w0 = hz(50.0)}, 1e-4).config.period == 200);
        CHECK(design::repetitive({.w0 = hz(60.0)}, 1.0 / 12000.0).config.period == 200);
        // The spec's q scales the binomial filter, so q = 1 is the binomial filter itself.
        const auto scaled = design::repetitive_binomial<1>({.w0 = hz(50.0), .q = 0.9}, 1e-3);
        REQUIRE(scaled.success);
        CHECK(scaled.config.q_dc_gain() == doctest::Approx(0.9));
        CHECK(scaled.config.q_side[0] == doctest::Approx(0.9 * 0.25));
    }

    TEST_CASE("track: the internal model is held and kept in phase, the error not learned, u recorded") {
        // Two laws learn the same periodic error, then one steps with e = 0 while the other tracks a
        // nonzero error: the tracked one stored w[k] = Q·w[k−N] as the stepped one did, and did not
        // learn the error it was not acting on, so the two run on as one.
        RepetitiveController<64, double> stepped(kRep);
        RepetitiveController<64, double> tracked(kRep);
        for (size_t k = 0; k < 100; ++k) {
            const double e = std::sin(2.0 * pi * static_cast<double>(k) / 20.0);
            CHECK(stepped.step(e) == doctest::Approx(tracked.step(e)).epsilon(1e-14));
        }
        for (size_t k = 0; k < 13; ++k) { // not a whole period: the phase must move with the ticks
            static_cast<void>(stepped.step(0.0));
            tracked.track(0.25, 5.0);
            CHECK(tracked.last_command() == 0.25);
        }
        tracked.back_calculate(0.25); // the correction track recorded: no clamp to report
        for (size_t k = 0; k < 60; ++k) {
            const double e = std::sin(2.0 * pi * static_cast<double>(k + 113) / 20.0);
            CHECK(tracked.step(e) == doctest::Approx(stepped.step(e)).epsilon(1e-14)); // a few ulps: two inlined copies of one tick
        }
    }

    TEST_CASE("Law<kRes, T> of a repetitive design sizes the buffer to the design's period") {
        static constinit Law<kRep, float>      rep;
        static constinit Law<kBinomial, float> binomial;
        static_assert(std::is_same_v<Law<kRep, float>::Runtime, RepetitiveController<20, float>>);          // N + 0
        static_assert(std::is_same_v<Law<kBinomial, float>::Runtime, RepetitiveController<202, float, 2>>); // N + M
        static_assert(sizeof(Law<kRep, float>) < sizeof(RepetitiveController<64, float>));
        RepetitiveController<64, float> wide{kRep};
        CHECK(rep.valid());
        CHECK(binomial.valid());
        for (size_t k = 0; k < 200; ++k) {
            const float e = static_cast<float>(std::sin(2.0 * pi * static_cast<double>(k) / 20.0));
            CHECK(rep.step(e) == doctest::Approx(wide.step(e)).epsilon(1e-6));
        }
    }
}

TEST_SUITE("Repetitive Law") {
    TEST_CASE("Law<kRes, T> of a constexpr repetitive design: the period, gain, taps and lead constants of the type") {
        //! The tuning is the type's: the Law stores only the internal model, its clamp flags, the ring's
        //! index and the last error and command, so it is smaller than its runtime of the same buffer.
        //! Neither has a set_gains: a new period is a new buffer.
        static_assert(sizeof(Law<kBinomial, float>) < sizeof(Law<kBinomial, float>::Runtime));
        static_assert(sizeof(Law<kRep, float>) < sizeof(Law<kRep, float>::Runtime));

        static constinit Law<kBinomial, float>          law;
        static constinit Law<kBinomial, float>::Runtime same{kBinomial};
        CHECK(law.config().period == 200);
        CHECK(law.config().q_half == 2);
        CHECK(law.config().gain == 0.5f);
        const auto e_at = [](size_t k) { return static_cast<float>(std::sin(2.0 * pi * static_cast<double>(k) / 200.0) + (0.3 * std::sin(2.0 * pi * static_cast<double>(k) / 40.0))); };
        for (size_t k = 0; k < 2000; ++k) {
            const float e = e_at(k);
            const float u = law.step(e);
            const float v = same.step(e);
            CHECK(u == doctest::Approx(v).epsilon(1e-5)); // two inlined copies of the tick may round differently under -ffast-math, and the model carries it a period on
            if (k % 97 == 0) {
                law.back_calculate(0.5f * u); // a clamp after it, as the law it adds to reports
                same.back_calculate(0.5f * u);
            }
        }
        //! It has no limits of its own: a clamp is the law's it adds to, reported with back_calculate (above).
        // law = Law<kBinomial, float>::Runtime{kBinomial}; // ill-formed: nothing replaces the type's tuning (Law::operator=)
        //! The Law's learned model carries over to a Runtime.
        Law<kBinomial, float>::Runtime copied{law};
        CHECK(copied.valid());
        CHECK(copied.config().period == 200);
        for (size_t k = 2000; k < 2400; ++k) {
            const float e = e_at(k);
            const auto  a = copied.step(e);
            const auto  b = law.step(e);
            CHECK(a == doctest::Approx(b).epsilon(1e-5)); // two inlined copies of the tick may round differently under -ffast-math, and the law's memory carries it on
        }
    }
}
