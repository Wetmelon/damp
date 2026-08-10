// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>

#include "damp/controllers/composition.hpp"

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
    [[nodiscard]] constexpr T control(T r, T y) {
        last_r = r;
        last_y = y;
        return k * (r - y);
    }
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
        CHECK(cas.control(1.0, 0.0, 0.0) == doctest::Approx(6.0));
        cas.reset();
        static_assert(SisoController<decltype(cas), double>);
    }

    TEST_CASE("Cascade 2-arg control uses y for both loops") {
        GainController<double>                                          outer{2.0};
        GainController<double>                                          inner{3.0};
        Cascade<GainController<double>, GainController<double>, double> cas{outer, inner};
        // Same y for outer and inner: u_o = 2*(1-0)=2; u = 3*(2-0)=6
        CHECK(cas.control(1.0, 0.0) == doctest::Approx(6.0));
    }

    TEST_CASE("SwitchedController modes") {
        GainController<float>                                                                          n{1.0f};
        GainController<float>                                                                          e{5.0f};
        GainController<float>                                                                          b{0.1f};
        SwitchedController<GainController<float>, GainController<float>, GainController<float>, float> sw{
            n, e, b
        };
        CHECK(sw.control(1.0f, 0.0f) == doctest::Approx(1.0f));
        sw.set_mode(SwitchMode::Experiment);
        CHECK(sw.control(1.0f, 0.0f) == doctest::Approx(5.0f));
        sw.set_mode(SwitchMode::Backup);
        CHECK(sw.control(1.0f, 0.0f) == doctest::Approx(0.1f));
        sw.inject(99.0f); // ignored outside Experiment
        CHECK(sw.last_u() == doctest::Approx(0.1f));
        sw.set_mode(SwitchMode::Experiment);
        sw.inject(7.0f);
        CHECK(sw.last_u() == doctest::Approx(7.0f));
    }

    TEST_CASE("SwitchedController set_mode resets newly selected path") {
        GainController<double>                                                                             n{1.0};
        GainController<double>                                                                             e{2.0};
        GainController<double>                                                                             b{3.0};
        SwitchedController<GainController<double>, GainController<double>, GainController<double>, double> sw{
            n, e, b
        };
        (void)sw.control(1.0, 0.25); // warm normal (last_y = 0.25)
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
    }
}
