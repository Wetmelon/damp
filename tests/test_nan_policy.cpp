// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// NaN/Inf policy: from finite inputs the library never produces NaN or ±Inf, and a
// runtime object built from a failed (or default) design emits finite, safe output.
// Each case below is a path that used to break the rule.

#include <cstddef>
#include <limits>

#include "damp/controllers/adrc.hpp"
#include "damp/controllers/esc.hpp"
#include "damp/controllers/harmonic_suppression.hpp"
#include "damp/controllers/lead_lag.hpp"
#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/controllers/offset_free_mpc.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/controllers/pr.hpp"
#include "damp/controllers/repetitive.hpp"
#include "damp/controllers/smc.hpp"
#include "damp/controllers/smith_predictor.hpp"
#include "damp/controllers/stsmc.hpp"
#include "damp/estimation/ukf.hpp"
#include "damp/filters/complementary.hpp"
#include "damp/filters/moving_average.hpp"
#include "damp/math/math.hpp"
#include "damp/math/trig.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/toolbox/conditioning.hpp"
#include "damp/toolbox/scaling.hpp"
#include "damp/toolbox/thermistor.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

// Finite edge inputs: zero, unit, tiny, huge, both signs.
constexpr float kEdges[] = {0.0f, 1.0f, -1.0f, 1e-30f, -1e-30f, 1e30f, -1e30f};

template<size_t N, typename T>
bool all_finite(const ColVec<N, T>& v) {
    for (size_t i = 0; i < N; ++i) {
        if (!damp::isfinite(v[i])) {
            return false;
        }
    }
    return true;
}

} // namespace

TEST_SUITE("NaN/Inf policy") {
    TEST_CASE("runtime exp / pow saturate instead of overflowing to Inf") {
        float x = 200.0f; // not a constant expression: the backend path
        CHECK(damp::isfinite(damp::exp(x)));
        CHECK(damp::exp(x) == std::numeric_limits<float>::max());
        double xd = 1000.0;
        CHECK(damp::isfinite(damp::exp(xd)));
        CHECK(damp::isfinite(damp::pow(10.0f, x)));
        CHECK(damp::isfinite(damp::pow(1e-30f, -3)));
        CHECK(damp::pow(-1e30f, 3) == -std::numeric_limits<float>::max());
        CHECK(damp::exp(1.0f) == doctest::Approx(2.718281828f)); // unchanged in range
    }

    TEST_CASE("fast sin / cos stay finite and defined for huge angles") {
        // Past 2^30 half-turns the period index no longer fits an int (UB before).
        for (const double a : {1e10, -1e10, 1e300, -1e300}) {
            const auto [s, c] = detail::fast_sincos(a);
            CHECK(damp::isfinite(s));
            CHECK(damp::isfinite(c));
            CHECK(damp::isfinite(detail::fast_sin(a)));
        }
        for (const float a : {1e10f, -1e10f, 1e38f}) {
            CHECK(damp::isfinite(detail::fast_cos(a)));
        }
        CHECK(detail::fast_sin(1.0) == doctest::Approx(0.8414709848)); // unchanged in range
    }

    TEST_CASE("scalar complementary filter with zero tau and zero dt follows the measurement") {
        Complementary<float> f{0.0f};
        (void)f(1.0f, 0.0f, 0.01f);
        CHECK(f(2.0f, 5.0f, 0.0f) == 2.0f);
        Complementary<float> g{0.5f};
        (void)g(1.0f, 0.0f, 0.01f);
        CHECK(g(2.0f, 5.0f, 0.0f) == 1.0f); // dt = 0: hold
    }

    TEST_CASE("thermistor past its hot end reads the largest finite temperature") {
        const Thermistor<float> ntc{design::beta(10000.0f, 298.15f, 3950.0f)};
        for (const float r : {0.0f, 1e-6f, 0.01f, 1e-30f}) {
            const float t = ntc.kelvin(r);
            CHECK(damp::isfinite(t));
            CHECK(t > 1000.0f); // shorted sensor: hot, so protection trips
        }
        CHECK(ntc.kelvin(10000.0f) == doctest::Approx(298.15f));
        const Thermistor<float> blank{};
        CHECK(damp::isfinite(blank.kelvin(1000.0f)));
    }

    TEST_CASE("affine calibration without a slope has a finite inverse and fit") {
        const AffineCal<float> flat{0.0f, 3.0f};
        CHECK(flat.invert(5.0f) == 0.0f);
        const auto degenerate = two_point_cal(2.0f, 10.0f, 2.0f, 20.0f);
        CHECK(degenerate.gain == 0.0f);
        CHECK(degenerate.apply(7.0f) == 10.0f);
        CHECK(two_point_cal(0.0f, 0.0f, 2.0f, 4.0f).gain == 2.0f);
    }

    TEST_CASE("scaled dead band of the whole span passes nothing") {
        CHECK(scaled_deadband(1.5f, 1.0f) == 0.0f);
        CHECK(scaled_deadband(-3.0f, 2.0f) == 0.0f);
        CHECK(scaled_deadband(1.0f, 0.5f) == doctest::Approx(1.0f));
    }

    TEST_CASE("comb notch window never casts an out-of-range length") {
        CHECK(comb_notch_window(1e6f, 1e-30f) == 0);
        CHECK(comb_notch_window(1000.0f, 50.0f) == 20);
    }

    TEST_CASE("UKF with alpha = 0 falls back to the plain unscented transform") {
        UnscentedParams<float> p{};
        p.alpha = 0.0f;
        UnscentedKalmanFilter<2, 1, 1, float> ukf{ColVec<2, float>{1.0f, 2.0f}, Matrix<2, 2, float>::identity(), Matrix<2, 2, float>::identity() * 0.01f, p};
        const auto                            f = [](const ColVec<2, float>& x, const ColVec<1, float>&) { return x; };
        const auto                            h = [](const ColVec<2, float>& x, const ColVec<1, float>&) { return ColVec<1, float>{x[0]}; };
        CHECK(ukf.predict(f));
        CHECK(ukf.update(h, ColVec<1, float>{1.5f}, Matrix<1, 1, float>{{0.1f}}));
        CHECK(all_finite(ukf.state()));
        CHECK(ukf.state()[0] == doctest::Approx(1.45f).epsilon(0.05));
    }

    TEST_CASE("MPC from a failed design holds its input instead of dividing by zero scales") {
        MPC<2, 1, 1, 10, 3, 0, float> mpc{design::MPCArtifacts<2, 1, 1, 10, 3, 0, float>{}};
        for (const float r : kEdges) {
            const auto u = mpc.control(ColVec<1, float>{r}, ColVec<2, float>{r, -r});
            CHECK(all_finite(u));
            CHECK(u[0] == 0.0f);
        }
        OffsetFreeMPC<2, 1, 1, 10, 3, float> ofmpc{};
        CHECK(all_finite(ofmpc.control(ColVec<1, float>{1.0f}, ColVec<1, float>{0.5f})));
    }

    TEST_CASE("controllers from failed designs emit finite output on edge inputs") {
        PIDController<float, PIDMode::PID> pid{design::DiscretePIDResult<float>{}};
        PRController<float>                pr{design::PRResult<float>{}};
        LeadLagController<float>           lead{design::LeadLagResult<float>{}};
        ADRCController<2, float>           adrc{design::ADRCResult<2, float>{}};
        SMCController<float>               smc{design::SMCResult<float>{}};
        STSMCController<float>             stsmc{design::STSMCResult<float>{}};
        ExtremumSeekingController<float>   esc{design::ESCResult<float>{}};
        RepetitiveController<16, float>    rep{design::RepetitiveResult<float>{}};
        SmithPredictor<8, float>           smith{design::SmithPredictorResult<float>{}};
        HarmonicSuppressor<3, float>       harm{design::HarmonicSuppressorResult<3, float>{}};
        LQI<2, 1, 1, float>                lqi{design::LQIResult<2, 1, 1, float>{}};
        LQGI<2, 1, 1, float, 2, 1>         lqgi{design::LQGIResult<2, 1, 1, float, 2, 1>{}};
        LQG<2, 1, 1, float, 2, 1>          lqg{design::LQGResult<2, 1, 1, float, 2, 1>{}};
        for (const float r : kEdges) {
            for (const float y : kEdges) {
                CHECK(damp::isfinite(pid.control(r, y)));
                CHECK(damp::isfinite(pr.control(r, y)));
                CHECK(damp::isfinite(lead.control(r, y)));
                CHECK(damp::isfinite(adrc.control(r, y)));
                CHECK(damp::isfinite(smc.control(r, y)));
                CHECK(damp::isfinite(stsmc.control(r, y)));
                CHECK(damp::isfinite(esc.step(y)));
                CHECK(damp::isfinite(rep.step(r - y)));
                CHECK(damp::isfinite(smith.control(r, y)));
                CHECK(damp::isfinite(harm.control(r, y)));
                const ColVec<1, float> rv{r};
                const ColVec<1, float> yv{y};
                CHECK(all_finite(lqi.control(rv, yv, ColVec<2, float>{y, y})));
                CHECK(all_finite(lqgi.step(rv, yv)));
                CHECK(all_finite(lqg.step(yv)));
            }
        }
    }
}
