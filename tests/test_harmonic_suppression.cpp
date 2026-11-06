// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

#include "damp/backend.hpp"
#include "damp/controllers/harmonic_suppression.hpp"
#include "damp/controllers/pr.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {
constexpr double pi = 3.14159265358979323846;
} // namespace

TEST_SUITE("harmonic_suppression") {
    TEST_CASE("synthesize validates the spec") {
        const double                 w0 = 2.0 * pi * 50.0; // 50 Hz fundamental
        const double                 Ts = 1.0 / 6000.0;    // 6 kHz loop
        const damp::array<size_t, 4> harmonics{1, 5, 7, 11};

        SUBCASE("valid spec succeeds and tunes one resonator per harmonic") {
            const auto d = design::harmonic_suppressor({.Kp = 0.5, .Ki = 200.0, .w0 = w0, .wc = 5.0}, Ts, harmonics);
            CHECK(d.success);
            // Fundamental carries Kp; harmonics carry Kp = 0.
            CHECK(d.gains[0].Kp == doctest::Approx(0.5));
            CHECK(d.gains[1].Kp == doctest::Approx(0.0));
            CHECK(d.gains[3].w0 == doctest::Approx(w0 * 11.0));
        }
        SUBCASE("rejects non-positive Ts / frequency") {
            CHECK_FALSE(design::harmonic_suppressor({.Kp = 0.5, .Ki = 200.0, .w0 = w0, .wc = 5.0}, 0.0, harmonics).success);
            CHECK_FALSE(design::harmonic_suppressor({.Kp = 0.5, .Ki = 200.0, .w0 = -w0, .wc = 5.0}, Ts, harmonics).success);
        }
        SUBCASE("rejects a 0th harmonic and any harmonic at/above Nyquist") {
            CHECK_FALSE(
                design::harmonic_suppressor({.Kp = 0.5, .Ki = 200.0, .w0 = w0, .wc = 5.0}, Ts, damp::array<size_t, 2>{1, 0}).success
            );
            // 61st harmonic of 50 Hz = 3050 Hz > Nyquist (3 kHz).
            CHECK_FALSE(
                design::harmonic_suppressor({.Kp = 0.5, .Ki = 200.0, .w0 = w0, .wc = 5.0}, Ts, damp::array<size_t, 1>{61}).success
            );
        }
    }

    TEST_CASE("tracks a multi-harmonic reference to zero steady-state error") {
        // Mirror the single-PR tracking test: unit-gain plant y[k] = u[k], the bank
        // closes the loop on the error. Each PR resonator has (near-)infinite gain at
        // its harmonic, so a reference built from exactly those harmonics is tracked
        // with vanishing steady-state error.
        const double                 fs = 10000.0;
        const double                 Ts = 1.0 / fs;
        const double                 w0 = 2.0 * pi * 50.0;
        const damp::array<size_t, 2> harmonics{5, 7};

        const auto d = design::harmonic_suppressor({.Kp = 0.0, .Ki = 500.0, .w0 = w0, .wc = 10.0}, Ts, harmonics);
        REQUIRE(d.success);
        HarmonicSuppressor<2, double> sup(d);

        auto reference = [&](size_t k) {
            const double t = static_cast<double>(k) * Ts;
            return 1.0 * std::sin(2.0 * pi * 250.0 * t)  // 5th
                 + 0.7 * std::sin(2.0 * pi * 350.0 * t); // 7th
        };

        double       y = 0.0;
        double       error_sum = 0.0;
        const size_t n_steps = 20000; // 2 s at 10 kHz
        for (size_t k = 0; k < n_steps; ++k) {
            const double e = reference(k) - y;
            y = sup.control(e); // unit-gain plant: y[k] = u[k]
            if (k > 15000) {    // settled window
                error_sum += std::abs(e);
            }
        }
        CHECK(error_sum / 5000.0 < 0.05); // both targeted harmonics tracked
    }

    TEST_CASE("invalid design is inert (valid() false, finite output)") {
        const double                  Ts = 1.0 / 6000.0;
        const auto                    d = design::harmonic_suppressor({.Kp = 0.5, .Ki = 200.0, .w0 = -1.0, .wc = 5.0}, Ts, damp::array<size_t, 1>{5});
        HarmonicSuppressor<1, double> sup(d);
        CHECK_FALSE(sup.valid());
    }

    TEST_CASE("float specialization compiles and runs") {
        const float                  Ts = 1.0f / 6000.0f;
        const float                  w0 = 2.0f * 3.14159265f * 50.0f;
        const damp::array<size_t, 2> harmonics{1, 5};
        const auto                   d = design::harmonic_suppressor({.Kp = 0.5f, .Ki = 50.0f, .w0 = w0, .wc = 5.0f}, Ts, harmonics);
        REQUIRE(d.success);
        HarmonicSuppressor<2, float> sup(d);
        const float                  u = sup.control(1.0f);
        CHECK(std::isfinite(u));
    }

    TEST_CASE("as<U>() converts the design result precision") {
        const double                 w0 = 2.0 * pi * 50.0;
        const double                 Ts = 1.0 / 6000.0;
        const damp::array<size_t, 2> harmonics{1, 5};
        const auto                   d = design::harmonic_suppressor({.Kp = 0.5, .Ki = 200.0, .w0 = w0, .wc = 5.0}, Ts, harmonics);
        const auto                   df = d.template as<float>();
        CHECK(df.success);
        CHECK(df.gains[0].Kp == doctest::Approx(0.5f));
    }

    TEST_CASE("harmonic suppressor is constexpr-evaluable") {
        constexpr bool ok = []() consteval {
            const double                  w0 = 2.0 * 3.14159265358979323846 * 50.0;
            const double                  Ts = 1.0 / 6000.0;
            const damp::array<size_t, 2>  harmonics{1, 5};
            const auto                    d = design::harmonic_suppressor({.Kp = 0.5, .Ki = 100.0, .w0 = w0, .wc = 5.0}, Ts, harmonics);
            HarmonicSuppressor<2, double> sup(d);
            double                        u = 0.0;
            for (int k = 0; k < 5; ++k) {
                u = sup.control(1.0);
            }
            return d.success && damp::abs(u) < 1e6; // finite, bounded
        }();
        static_assert(ok, "harmonic suppressor must work at compile time");
        CHECK(ok);
    }

    TEST_CASE("Kp is applied once even if every resonator carries Kp") {
        const double                 w0 = 2.0 * pi * 50.0;
        const double                 Ts = 1.0 / 10000.0;
        const damp::array<size_t, 3> harmonics{1, 3, 5};
        const auto                   once = design::harmonic_suppressor({.Kp = 1.0, .Ki = 100.0, .w0 = w0, .wc = 10.0}, Ts, harmonics);
        auto                         stuffed = once;
        for (size_t i = 0; i < stuffed.gains.size(); ++i) {
            stuffed.gains[i].Kp = 1.0;
        }
        HarmonicSuppressor<3, double> a(once);
        HarmonicSuppressor<3, double> b(stuffed);
        CHECK(a.Kp() == doctest::Approx(1.0));
        CHECK(b.Kp() == doctest::Approx(1.0));
        CHECK(a.control(1.0) == doctest::Approx(b.control(1.0)));
    }

    TEST_CASE("runtime control(r,y), reset, and as<U>() preserve the bank") {
        const double                 w0 = 2.0 * pi * 50.0;
        const damp::array<size_t, 3> harmonics{1, 3, 5};
        const auto                   results = design::harmonic_suppressor({.Kp = 1.0, .Ki = 100.0, .w0 = w0, .wc = 10.0}, 1.0 / 10000.0, harmonics);

        HarmonicSuppressor<3, double> ctrl(results);
        CHECK(ctrl.control(1.0) != 0.0);

        HarmonicSuppressor<3, double> a(results);
        HarmonicSuppressor<3, double> b(results);
        CHECK(a.control(0.5, 0.1) == doctest::Approx(b.control(0.4)));

        ctrl.reset();
        CHECK(ctrl.control(1.0) != 0.0);

        HarmonicSuppressor<3, double> d(results);
        (void)d.control(0.5);
        const auto f = d.as<float>();
        CHECK(f.Kp() == doctest::Approx(1.0f));
        CHECK(f.resonator(0).resonant.last_output() == doctest::Approx(static_cast<float>(d.resonator(0).resonant.last_output())));
    }

    TEST_CASE("to_ss / to_tf DC gain is the fundamental Kp") {
        const double                 w0 = 2.0 * pi * 50.0;
        const double                 Ts = 1.0 / 6000.0;
        const damp::array<size_t, 2> harmonics{1, 5};
        const auto                   d = design::harmonic_suppressor({.Kp = 0.5, .Ki = 200.0, .w0 = w0, .wc = 5.0}, Ts, harmonics);
        REQUIRE(d.success);
        const auto ss = d.to_ss();
        const auto tf = d.to_tf();
        CHECK(ss.D(0, 0) == doctest::Approx(0.5));
        // Highest power first: DC gain is the constant-term ratio; den's s³ term is Σ 2wc.
        CHECK(tf.num[4] / tf.den[4] == doctest::Approx(0.5).epsilon(1e-8));
        CHECK(tf.den[1] == doctest::Approx(2.0 * 2.0 * 5.0).epsilon(1e-12));
        const auto G0 = *eval_frf(ss, damp::complex<double>{0.0, 0.0});
        CHECK(G0(0, 0).real() == doctest::Approx(0.5).epsilon(1e-8));
    }

    TEST_CASE("a bank without the fundamental retunes each harmonic; a failed resonator fails the bank") {
        const double                  w0 = 2.0 * pi * 50.0;
        const double                  Ts = 1.0 / 10000.0;
        const damp::array<size_t, 2>  harmonics{5, 7};
        HarmonicSuppressor<2, double> bank(design::harmonic_suppressor({.Kp = 0.0, .Ki = 100.0, .w0 = w0, .wc = 5.0}, Ts, harmonics));
        REQUIRE(bank.valid());
        bank.set_fundamental(2.0 * pi * 51.0);
        CHECK(bank.resonator(0).w0 == doctest::Approx(2.0 * pi * 255.0));
        CHECK(bank.resonator(1).w0 == doctest::Approx(2.0 * pi * 357.0));
        // A fundamental that would push the 7th past Nyquist (5 kHz) leaves the whole bank alone.
        bank.set_fundamental(2.0 * pi * 800.0);
        CHECK(bank.resonator(0).w0 == doctest::Approx(2.0 * pi * 255.0));

        const double inf = std::numeric_limits<double>::infinity();
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const auto   bad_ki = design::harmonic_suppressor({.Kp = 0.5, .Ki = inf, .w0 = w0, .wc = 5.0}, Ts, harmonics);
        CHECK_FALSE(bad_ki.success);
        CHECK(bad_ki.status == design::DesignStatus::invalid_parameter);
        CHECK_FALSE(design::harmonic_suppressor({.Kp = nan, .Ki = 100.0, .w0 = w0, .wc = 5.0}, Ts, harmonics).success);
        CHECK(design::harmonic_suppressor({.Kp = 0.5, .Ki = 100.0, .w0 = w0, .wc = 5.0}, nan, harmonics).status == design::DesignStatus::invalid_sample_time);

        // A hand-edited result whose resonator failed leaves the runtime invalid.
        auto edited = design::harmonic_suppressor({.Kp = 0.5, .Ki = 100.0, .w0 = w0, .wc = 5.0}, Ts, harmonics);
        edited.gains[1] = design::PRResult<double>{};
        HarmonicSuppressor<2, double> inert(edited);
        CHECK_FALSE(inert.valid());
        CHECK(inert.control(1.0) == 0.0);
    }
}
