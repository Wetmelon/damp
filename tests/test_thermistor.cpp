// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "damp/backend.hpp"
#include "damp/toolbox/thermistor.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

// NTC thermistor linearization: the design:: coefficient fits (Beta model and
// three-point Steinhart-Hart) plus the runtime Thermistor evaluation.

TEST_SUITE("Thermistor") {
    TEST_CASE("Beta model reads reference point and slope sign") {
        // 10 kΩ NTC, β = 3950, 25 °C (298.15 K) reference.
        const Thermistor<double> ntc{design::beta(10000.0, 298.15, 3950.0)};
        // At the reference resistance it must read the reference temperature.
        CHECK(ntc.celsius(10000.0) == doctest::Approx(25.0).epsilon(1e-6));
        // NTC: higher resistance -> colder; lower resistance -> hotter.
        CHECK(ntc.celsius(30000.0) < 25.0);
        CHECK(ntc.celsius(3000.0) > 25.0);
    }

    TEST_CASE("Steinhart-Hart fit reproduces its calibration points") {
        // Three (R, T) points generated from a known Beta model: take three
        // resistances and read their true temperature off the Beta curve. The S-H
        // fit must pass through all three points exactly.
        const Thermistor<double> truth{design::beta(10000.0, 298.15, 3950.0)};
        const double             r1 = 30000.0;
        const double             r2 = 10000.0;
        const double             r3 = 3000.0;
        const double             t1 = truth.kelvin(r1);
        const double             t2 = truth.kelvin(r2);
        const double             t3 = truth.kelvin(r3);

        const Thermistor<double> sh{
            design::steinhart_hart<double>({r1, t1}, {r2, t2}, {r3, t3})
        };
        CHECK(sh.kelvin(r1) == doctest::Approx(t1).epsilon(1e-9));
        CHECK(sh.kelvin(r2) == doctest::Approx(t2).epsilon(1e-9));
        CHECK(sh.kelvin(r3) == doctest::Approx(t3).epsilon(1e-9));
    }

    TEST_CASE("direct coefficient construction matches design::beta") {
        // The (a, b, c) constructor is equivalent to feeding the same coefficients
        // through design::beta (which is just the c = 0 special case).
        const auto               coeffs = design::beta(10000.0, 298.15, 3950.0);
        const Thermistor<double> direct{coeffs.a, coeffs.b, coeffs.c};
        const Thermistor<double> viacoeffs{coeffs};
        CHECK(direct.kelvin(15000.0) == doctest::Approx(viacoeffs.kelvin(15000.0)));
    }
}
