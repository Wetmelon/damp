// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>

#include "damp/concepts.hpp"
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
#include "damp/estimation/ekf.hpp"
#include "damp/estimation/excitation.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/estimation/luenberger.hpp"
#include "damp/estimation/mhe.hpp"
#include "damp/estimation/parameter_estimation.hpp"
#include "damp/estimation/relay_autotune.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

// ============================================================================
// Conformance: every runtime pinned to its taxonomy role. API drift on any of
// these becomes a compile error here, not a surprise in user code.
// ============================================================================

// --- SisoController: u = control(r, y) --------------------------------------
static_assert(SisoController<PIDController<double, PIDMode::PID>, double>);
static_assert(SisoController<PIDController<double, PIDMode::PI>, double>);
static_assert(SisoController<PIDController<float, PIDMode::P>, float>);
static_assert(SisoController<ADRCController<2, float>, float>);
static_assert(SisoController<SMCController<double>, double>);
static_assert(SisoController<STSMCController<float>, float>);
static_assert(SisoController<LeadLagController<double>, double>);
static_assert(SisoController<PRController<double>, double>);
static_assert(SisoController<MultiPRController<3, double>, double>);
static_assert(SisoController<HarmonicSuppressor<3, double>, double>);
static_assert(SisoController<SmithPredictor<16, double>, double>); // default PI primary
static_assert(SisoController<SmithPredictor<16, float, PIController<float>>, float>);

// Not SISO controllers in this vocabulary (deliberate):
static_assert(!SisoController<ExtremumSeekingController<double>, double>); // an optimizer: step(J)
static_assert(!SisoController<RepetitiveController<64, double>, double>);  // a plug-in add-on: step(e)
static_assert(!SisoController<RelayAutotuner<double>, double>);            // an experiment: step(y) + status

// --- OutputFeedbackController: u = control(r, y), self-contained tick -------
static_assert(OutputFeedbackController<OffsetFreeMPC<2, 1, 1, 10, 4, double>, 1, 1, double>);
// LQGI::control(r, y) matches the syntax only — its Kalman filter is caller-sequenced
// (see the @warning on it). LQGI::step(r, y) is the self-contained tick, and
// feedback/commit the saturation-aware split. Pinned here so the trap stays visible.
static_assert(OutputFeedbackController<LQGI<2, 1, 1, double, 2, 1>, 1, 1, double>);
static_assert(!OutputFeedbackController<LQG<2, 1, 1, double, 2, 1>, 1, 1, double>);     // law+estimator pair
static_assert(!OutputFeedbackController<MPC<2, 1, 1, 10, 4, 0, double>, 1, 1, double>); // needs x, not y

// --- StateFeedbackController: u = control(r, x), reference-first ------------
static_assert(StateFeedbackController<StateFeedback<2, 1, double>, 1, 2, 2, double>);    // r = x_ref
static_assert(StateFeedbackController<LQR<2, 1, double>, 1, 2, 2, double>);              // LQR is an alias of StateFeedback
static_assert(StateFeedbackController<MPC<2, 1, 1, 10, 4, 0, double>, 1, 1, 2, double>); // r = y_ref
static_assert(!StateFeedbackController<LQI<2, 1, 1, double>, 1, 1, 2, double>);          // tick is (r, y, x)

// --- StateEstimator: estimate(y, u) + state() --------------------------------
static_assert(StateEstimator<KalmanFilter<2, 1, 1, double, 2, 1>, 2, 1, 1, double>);
static_assert(StateEstimator<SteadyStateKalmanFilter<2, 1, 1, double, 2, 1>, 2, 1, 1, double>);
static_assert(StateEstimator<Luenberger<2, 1, 1, double>, 2, 1, 1, double>);
static_assert(StateEstimator<MHE<1, 1, 1, 4, double>, 1, 1, 1, double>);
static_assert(!StateEstimator<ExtendedKalmanFilter<2, 1, 1, double>, 2, 1, 1, double>); // model-callback API

// --- SignalSource: step() + done() -------------------------------------------
static_assert(SignalSource<Chirp<double>, double>);
static_assert(SignalSource<PRBS<double>, double>);
static_assert(SignalSource<StepTrain<double>, double>);
static_assert(SignalSource<Ramp<double>, double>);
static_assert(SignalSource<MultiSine<3, double>, double>);
static_assert(SignalSource<SteppedSine<4, double>, double>);
static_assert(SignalSource<Impulse<double>, double>);
static_assert(SignalSource<Step<double>, double>);

// --- ParameterEstimator: parameters() + valid() + confidence() ---------------
static_assert(ParameterEstimator<FirstOrderPlantEstimator<double>, 2, double>);
static_assert(ParameterEstimator<FirstOrderPlantEstimator<float>, 2, float>);

// ============================================================================
// Runtime equivalence: the normalized forms are exactly the originals.
// ============================================================================

TEST_CASE("stored-Ts 2-arg control equals the explicit 3-arg form") {
    constexpr double Ts = 1e-3;

    SUBCASE("PID (discrete vs ContinuousPID at fixed Ts)") {
        design::PIDResult<double> gains{};
        gains.Kp = 2.0;
        gains.Ki = 10.0;
        gains.Kd = 0.05;

        PIDController<double, PIDMode::PID> disc{gains.discretize(Ts)};
        ContinuousPID<double>               cont{gains};
        for (int k = 0; k < 50; ++k) {
            const double r = 1.0;
            const double y = 0.1 * k * Ts;
            CHECK(disc.control(r, y) == doctest::Approx(cont.control(r, y, Ts)).epsilon(1e-15));
        }
    }

    SUBCASE("ADRC") {
        const auto art = design::adrc<1>(4.0, 20.0, 1.0);
        REQUIRE(art.success);

        ADRCController<1, double> two{art, Ts};
        ADRCController<1, double> three{art};
        for (int k = 0; k < 50; ++k) {
            const double y = 0.05 * k * Ts;
            CHECK(two.control(1.0, y) == doctest::Approx(three.control(1.0, y, Ts)).epsilon(1e-15));
        }
    }

    SUBCASE("SMC") {
        const auto gains = design::smc(20.0, 5.0, 1.0);

        SMCController<double> two{gains, Ts};
        SMCController<double> three{gains};
        for (int k = 0; k < 50; ++k) {
            const double y = 0.02 * k;
            CHECK(two.control(1.0, y) == doctest::Approx(three.control(1.0, y, Ts)).epsilon(1e-15));
        }
    }

    SUBCASE("super-twisting") {
        const auto art = design::stsmc(2.0, 5.0);
        REQUIRE(art.success);

        STSMCController<double> two{art, Ts};
        STSMCController<double> three{art};
        for (int k = 0; k < 50; ++k) {
            const double y = 0.02 * k;
            CHECK(two.control(1.0, y) == doctest::Approx(three.control(1.0, y, Ts)).epsilon(1e-15));
        }
    }
}

TEST_CASE("fused estimate(y, u) equals the manual predict/update sequence") {
    constexpr double Ts = 0.1;

    StateSpace<2, 1, 1, double, 2, 1> sys{};
    sys.A = Matrix<2, 2>{{1.0, Ts}, {0.0, 1.0}};
    sys.B = Matrix<2, 1>{{0.5 * Ts * Ts}, {Ts}};
    sys.C = Matrix<1, 2>{{1.0, 0.0}};
    sys.Ts = Ts;

    SUBCASE("KalmanFilter") {
        const auto kd = design::kalman(sys, Matrix<2, 2>::identity(), Matrix<1, 1>{{0.01}});
        REQUIRE(kd.success);

        KalmanFilter<2, 1, 1, double, 2, 1> fused{kd};
        KalmanFilter<2, 1, 1, double, 2, 1> manual{kd};
        for (int k = 0; k < 30; ++k) {
            const ColVec<1> y{0.1 * k};
            const ColVec<1> u{0.5};
            const auto      xf = fused.estimate(y, u);
            manual.predict(u);
            (void)manual.update(y, u);
            CHECK(xf(0) == doctest::Approx(manual.state()(0)).epsilon(1e-15));
            CHECK(xf(1) == doctest::Approx(manual.state()(1)).epsilon(1e-15));
        }
    }

    SUBCASE("Luenberger") {
        const Matrix<2, 1>          L{{0.4}, {0.8}};
        Luenberger<2, 1, 1, double> fused{sys, L};
        Luenberger<2, 1, 1, double> manual{sys, L};
        for (int k = 0; k < 30; ++k) {
            const ColVec<1> y{0.1 * k};
            const ColVec<1> u{0.5};
            const auto      xf = fused.estimate(y, u);
            manual.step(y, u);
            CHECK(xf(0) == doctest::Approx(manual.state()(0)).epsilon(1e-15));
            CHECK(xf(1) == doctest::Approx(manual.state()(1)).epsilon(1e-15));
        }
    }

    SUBCASE("MHE") {
        StateSpace<1, 1, 1> tank{};
        tank.A = Matrix<1, 1>{{0.9}};
        tank.B = Matrix<1, 1>{{0.1}};
        tank.C = Matrix<1, 1>{{1.0}};
        tank.Ts = Ts;
        const auto art = design::mhe<4>(tank, Matrix<1, 1>{{1e-4}}, Matrix<1, 1>{{1e-3}});
        REQUIRE(art.success);

        MHE fused{art};
        MHE manual{art};
        for (int k = 0; k < 20; ++k) {
            const ColVec<1> y{0.3 + (0.01 * k)};
            const ColVec<1> u{0.1};
            const auto      xf = fused.estimate(y, u);
            const auto      xm = manual.update(y, u);
            CHECK(xf(0) == doctest::Approx(xm(0)).epsilon(1e-15));
            CHECK(fused.state()(0) == doctest::Approx(xf(0)));
        }
    }
}
