// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)


#include "damp/backend.hpp"
#include "damp/concepts.hpp"
#include "damp/controllers/composition.hpp"
#include "damp/controllers/pid.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

// Minimal SisoController for switch tests.
template<typename T>
struct GainController {
    T                         k{T{1}};
    T                         last_r{T{0}};
    T                         last_y{T{0}};
    int                       resets{0};
    int                       reports{0};
    [[nodiscard]] constexpr T step(T r, T y) {
        last_r = r;
        last_y = y;
        return k * (r - y);
    }
    constexpr void back_calculate(T /*u_sat*/, T /*u_ff*/ = T{0}) { ++reports; }
    constexpr void reset() {
        last_r = T{0};
        last_y = T{0};
        ++resets;
    }
};

} // namespace

TEST_SUITE("Cascade and SwitchedController") {
    TEST_CASE("Cascade outer then inner") {
        GainController<double>                                          outer{2.0};
        GainController<double>                                          inner{3.0};
        Cascade<GainController<double>, GainController<double>, double> cas{outer, inner};
        // r=1, y_o=0 → u_o = 2; y_i=0 → u = 3*2 = 6
        CHECK(cas.step(1.0, 0.0, 0.0) == doctest::Approx(6.0));
        cas.reset();
        // Two measurements: not a SisoController (that contract is one y).
        static_assert(!SisoController<decltype(cas), double>);
    }

    TEST_CASE("SwitchedController modes") {
        GainController<float>                                                                          n{1.0f};
        GainController<float>                                                                          e{5.0f};
        GainController<float>                                                                          b{0.1f};
        SwitchedController<GainController<float>, GainController<float>, GainController<float>, float> sw{
            n, e, b
        };
        CHECK(sw.step(1.0f, 0.0f) == doctest::Approx(1.0f));
        sw.set_mode(SwitchMode::Experiment);
        CHECK(sw.step(1.0f, 0.0f) == doctest::Approx(5.0f));
        sw.set_mode(SwitchMode::Backup);
        CHECK(sw.step(1.0f, 0.0f) == doctest::Approx(0.1f));
        sw.inject(99.0f); // ignored outside Experiment
        CHECK(sw.last_command() == doctest::Approx(0.1f));
        sw.set_mode(SwitchMode::Experiment);
        sw.inject(7.0f);
        CHECK(sw.last_command() == doctest::Approx(7.0f));
    }

    TEST_CASE("SwitchedController set_mode resets newly selected path") {
        GainController<double>                                                                             n{1.0};
        GainController<double>                                                                             e{2.0};
        GainController<double>                                                                             b{3.0};
        SwitchedController<GainController<double>, GainController<double>, GainController<double>, double> sw{
            n, e, b
        };
        (void)sw.step(1.0, 0.25); // warm normal (last_y = 0.25)
        CHECK(sw.normal().last_y == doctest::Approx(0.25));

        const int e_before = sw.experiment().resets;
        sw.set_mode(SwitchMode::Experiment);
        CHECK(sw.experiment().resets == e_before + 1);
        CHECK(sw.experiment().last_r == doctest::Approx(0.0));
        CHECK(sw.experiment().last_y == doctest::Approx(0.0));

        const int n_before = sw.normal().resets;
        sw.set_mode(SwitchMode::Normal);
        CHECK(sw.normal().resets == n_before + 1);
        CHECK(sw.normal().last_y == doctest::Approx(0.0));

        sw.back_calculate(0.5); // to the active law only
        CHECK(sw.normal().reports == 1);
        CHECK(sw.experiment().reports == 0);
    }

    TEST_CASE("Cascade holds the outer integral while the inner command is clamped") {
        const auto                                                 outer_gains = design::PIDGains<double>{.Kp = 1.0, .Ki = 10.0}.discretize(0.01);
        const auto                                                 inner_gains = design::PIDGains<double>{.Kp = 2.0}.discretize(0.01, {.input = {.lower = -1.0, .upper = 1.0}});
        Cascade<PIController<double>, PController<double>, double> cas{PIController<double>{outer_gains}, PController<double>{inner_gains}};
        for (int k = 0; k < 100; ++k) {
            CHECK(cas.step(5.0, 0.0, 0.0) == 1.0); // the inner loop sits on its limit
        }
        CHECK(cas.outer().integral() == 0.0); // was 100·Ki·Ts·5 = 50: it ran on while the inner loop could not follow

        // The same through a clamp after the cascade, reported with back_calculate.
        const auto                                                 open_inner = design::PIDGains<double>{.Kp = 2.0}.discretize(0.01);
        Cascade<PIController<double>, PController<double>, double> told{PIController<double>{outer_gains}, PController<double>{open_inner}};
        for (int k = 0; k < 100; ++k) {
            told.back_calculate(damp::min(told.step(5.0, 0.0, 0.0), 1.0));
        }
        CHECK(told.outer().integral() == 0.0);

        // Once the inner loop follows again, the outer integral moves.
        static_cast<void>(cas.step(0.4, 0.0, 0.5)); // r_inner = 0.44: inner 2·(0.44 − 0.5) inside ±1
        CHECK(cas.outer().integral() > 0.0);
    }

    TEST_CASE("SwitchedController starts a PID switched in from the command the plant last received") {
        const auto gains = design::PIDGains<double>{.Kp = 1.0, .Ki = 10.0}.discretize(0.01);
        using PI = PIController<double>;
        SwitchedController<PI, PI, PI, double> sw{PI{gains}, PI{gains}, PI{gains}};
        double                                 u = 0.0;
        for (int k = 0; k < 50; ++k) {
            u = sw.step(1.0, 0.0); // the integral builds: u = 1 + 0.1·50
        }
        CHECK(u == doctest::Approx(6.0));
        sw.set_mode(SwitchMode::Backup);
        CHECK(sw.step(0.5, 0.5) == doctest::Approx(u)); // at zero error: was 0, the backup's reset state

        sw.back_calculate(4.0); // the plant received 4
        sw.set_mode(SwitchMode::Normal);
        CHECK(sw.step(0.5, 0.5) == doctest::Approx(4.0));
    }

    TEST_CASE("SwitchedController hands over exactly at r ≠ y with a setpoint weight b ≠ 1 (#87)") {
        //! The old track(u, y) preloaded the integral as if r = y: the first command of a PID switched
        //! in was u + Kp·(b·r − y) − Kp·(b − 1)·y, 2.4 off here. track(u, r, y) preloads at this r.
        const auto gains = design::PIDGains<double>{.Kp = 2.0, .Ki = 10.0, .Kd = 0.05, .b = 0.4, .c = 0.0, .Tf = 0.01}.discretize(0.01, {.input = {.lower = -20.0, .upper = 20.0}});
        using Law = PIDController<double>;
        SwitchedController<Law, Law, Law, double> sw{Law{gains}, Law{gains}, Law{gains}};
        double                                    u = 0.0;
        for (int k = 0; k < 50; ++k) {
            u = sw.step(3.0, 0.02 * static_cast<double>(k));
        }
        sw.set_mode(SwitchMode::Backup);
        CHECK(sw.step(2.5, 1.3) == doctest::Approx(u).epsilon(1e-12)); // r and y both moved at the switch
        CHECK(sw.backup().last_command() == doctest::Approx(u).epsilon(1e-12));
        const double next = sw.step(2.5, 1.3);
        CHECK(damp::abs(next - u) < 0.2); // then the law's own: one tick of the integral, no jump
    }
}
