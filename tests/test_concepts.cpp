// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)


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
#include "damp/controllers/oflqg.hpp"
#include "damp/controllers/ofmpc.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/controllers/pr.hpp"
#include "damp/controllers/repetitive.hpp"
#include "damp/controllers/smc.hpp"
#include "damp/controllers/smith_predictor.hpp"
#include "damp/controllers/stsmc.hpp"
#include "damp/design/qp.hpp"
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
static_assert(SisoController<PIDController<double>, double>);
static_assert(SisoController<PIController<double>, double>);
static_assert(SisoController<PController<float>, float>);
static_assert(SisoController<ADRCController<2, float>, float>);
static_assert(SisoController<SMCController<double>, double>);
static_assert(SisoController<STSMCController<float>, float>);
static_assert(SisoController<LeadLagController<double>, double>);
static_assert(SisoController<PRController<double>, double>);
static_assert(SisoController<HarmonicSuppressor<3, double>, double>);
static_assert(SisoController<SmithPredictor<16, double>, double>); // default PI primary
static_assert(SisoController<SmithPredictor<16, float, PIController<float>>, float>);

// Not SISO controllers in this vocabulary (deliberate):
static_assert(!SisoController<ExtremumSeekingController<double>, double>); // an optimizer: control(J, valid)
static_assert(!SisoController<RepetitiveController<64, double>, double>);  // a plug-in add-on: control(e)
static_assert(!SisoController<RelayAutotuner<double>, double>);            // an experiment: step(y) + status

// --- FeedbackRuntime: lifecycle only (update, predict, back_calculate, reset) ---
// The seven LQ and MPC runtimes on a two-state, one-input, one-output plant.
using PlantFilter = Luenberger<2, 1, 1, double>;                                                                // on [x]
using LoadFilter = Luenberger<3, 1, 1, double, 1>;                                                              // on [x; d], one disturbance state
using LqrLaw = LQR<2, 1, double>;                                                                               // design::lqr
using LqrTrackLaw = LQR<2, 1, double, 2, 1>;                                                                    // design::lqr, one tracked output
using LqiLaw = LQI<2, 1, 1, double, 1>;                                                                         // design::lqi, one integrated output
using LqgLaw = LQG<2, 1, 1, double, PlantFilter>;                                                               // design::lqg
using LqgTrackLaw = LQG<2, 1, 1, double, PlantFilter, 1>;                                                       // design::lqg, one tracked output
using LqgiLaw = LQGI<2, 1, 1, 1, double, PlantFilter>;                                                          // design::lqgi
using OflqgLaw = OFLQG<2, 1, 1, 1, double, LoadFilter>;                                                         // design::oflqg, one controlled output
using MpcLaw = MPC<2, 1, 1, 10, 4, 0, double>;                                                                  // design::mpc
using OfmpcLaw = OFMPC<2, 1, 1, 10, 4, double, design::WarmStartActiveSetSolver<5, 37, double>, 0, LoadFilter>; // design::ofmpc

static_assert(FeedbackRuntime<LqrLaw, 1, 2, double>);   // m = x
static_assert(FeedbackRuntime<LqiLaw, 1, 2, double>);   // m = x
static_assert(FeedbackRuntime<LqgLaw, 1, 1, double>);   // m = y
static_assert(FeedbackRuntime<LqgiLaw, 1, 1, double>);  // m = y
static_assert(FeedbackRuntime<OflqgLaw, 1, 1, double>); // m = y
static_assert(FeedbackRuntime<MpcLaw, 1, 2, double>);   // m = x
static_assert(FeedbackRuntime<OfmpcLaw, 1, 1, double>); // m = y
static_assert(!FeedbackRuntime<LqrLaw, 1, 1, double>);  // needs x, not y
static_assert(!FeedbackRuntime<LqgLaw, 1, 2, double>);  // needs y, not x
static_assert(!FeedbackRuntime<MpcLaw, 1, 1, double>);  // needs x, not y

// FeedbackLaw: control(r), step(r, m). Regulators have no r.
static_assert(FeedbackLaw<LqiLaw, 1, 1, 2, double>);
static_assert(FeedbackLaw<LqgiLaw, 1, 1, 1, double>);
static_assert(FeedbackLaw<OflqgLaw, 1, 1, 1, double>);
static_assert(FeedbackLaw<MpcLaw, 1, 1, 2, double>);
static_assert(FeedbackLaw<OfmpcLaw, 1, 1, 1, double>);
static_assert(FeedbackLaw<LqrTrackLaw, 1, 1, 2, double>);
static_assert(FeedbackLaw<LqgTrackLaw, 1, 1, 1, double>);
static_assert(!FeedbackLaw<LqrLaw, 1, 1, 2, double>); // regulator: control(), step(x)
static_assert(!FeedbackLaw<LqgLaw, 1, 1, 1, double>);

// The references each law's control and step take. Every law has the basic call control(r),
// step(r, m), r the outputs' reference; LQR and LQG without tracked outputs have none, control()
// and step(m). The LQ laws add the feedforward form with (x_ref, u_ref), MPC the measured
// disturbance d when its design has one, and then only that form. No default arguments.
template<typename C, typename... A>
concept ControlsWith = requires(C c, const A&... a) { c.control(a...); };
template<typename C, typename... A>
concept StepsWith = requires(C c, const A&... a) { c.step(a...); };

using MpcMdLaw = MPC<2, 1, 1, 10, 4, 1, double>; // design::mpc with design::MeasuredDisturbance, one channel

using X2 = ColVec<2, double>; // x, x_ref
using V1 = ColVec<1, double>; // u_ref, y, r, d
using D0 = ColVec<0, double>; // a measured disturbance of no channels

// The basic call: control(r), step(r, m); the regulator control(), step(m).
static_assert(ControlsWith<LqrLaw> && StepsWith<LqrLaw, X2>);                   // control(), step(x)
static_assert(ControlsWith<LqgLaw> && StepsWith<LqgLaw, V1>);                   // control(), step(y)
static_assert(ControlsWith<LqrTrackLaw, V1> && StepsWith<LqrTrackLaw, V1, X2>); // control(r), step(r, x)
static_assert(ControlsWith<LqgTrackLaw, V1> && StepsWith<LqgTrackLaw, V1, V1>); // control(r), step(r, y)
static_assert(!ControlsWith<LqrLaw, X2> && !StepsWith<LqrLaw, X2, X2>);         // a state setpoint is the feedforward form's x_ref
static_assert(!ControlsWith<LqrTrackLaw> && !StepsWith<LqrTrackLaw, X2>);       // a tracking design always takes r
static_assert(ControlsWith<LqiLaw, V1> && StepsWith<LqiLaw, V1, X2>);           // control(r), step(r, x)
static_assert(ControlsWith<LqgiLaw, V1> && StepsWith<LqgiLaw, V1, V1>);         // control(r), step(r, y)
static_assert(ControlsWith<OflqgLaw, V1> && StepsWith<OflqgLaw, V1, V1>);       // control(r), step(r, y)
static_assert(ControlsWith<MpcLaw, V1> && StepsWith<MpcLaw, V1, X2>);           // control(r), step(r, x)
static_assert(ControlsWith<OfmpcLaw, V1> && StepsWith<OfmpcLaw, V1, V1>);       // control(r), step(r, y)

// The LQ feedforward form from a planned trajectory (x_ref, u_ref).
static_assert(ControlsWith<LqrLaw, X2, V1> && StepsWith<LqrLaw, X2, V1, X2>);                   // control(x_ref, u_ref), step(x_ref, u_ref, x)
static_assert(ControlsWith<LqiLaw, V1, X2, V1> && StepsWith<LqiLaw, V1, X2, V1, X2>);           // control(r, x_ref, u_ref), step(r, x_ref, u_ref, x)
static_assert(ControlsWith<LqgLaw, X2, V1> && StepsWith<LqgLaw, X2, V1, V1>);                   // control(x_ref, u_ref), step(x_ref, u_ref, y)
static_assert(ControlsWith<LqrTrackLaw, V1, X2, V1> && StepsWith<LqrTrackLaw, V1, X2, V1, X2>); // control(r, x_ref, u_ref), step(r, x_ref, u_ref, x)
static_assert(ControlsWith<LqgTrackLaw, V1, X2, V1> && StepsWith<LqgTrackLaw, V1, X2, V1, V1>); // control(r, x_ref, u_ref), step(r, x_ref, u_ref, y)
static_assert(ControlsWith<LqgiLaw, V1, X2, V1> && StepsWith<LqgiLaw, V1, X2, V1, V1>);         // control(r, x_ref, u_ref), step(r, x_ref, u_ref, y)
static_assert(ControlsWith<OflqgLaw, V1, X2, V1> && StepsWith<OflqgLaw, V1, X2, V1, V1>);       // control(r, x_ref, u_ref), step(r, x_ref, u_ref, y)
static_assert(!ControlsWith<LqiLaw, V1, X2> && !StepsWith<LqiLaw, V1, X2, X2>);                 // x_ref goes with u_ref, never alone

// MPC: the d argument follows the design. Without a measured disturbance there is none to pass,
// not even {}; with one, it is always passed.
template<typename C>
concept StepsWithEmptyD = requires(C c, const V1& r, const X2& x) { c.step(r, {}, x); };
template<typename C>
concept ControlsWithEmptyD = requires(C c, const V1& r) { c.control(r, {}); };

static_assert(!StepsWithEmptyD<MpcLaw> && !ControlsWithEmptyD<MpcLaw>);
static_assert(!ControlsWith<MpcLaw, V1, D0> && !StepsWith<MpcLaw, V1, D0, X2>);
static_assert(StepsWithEmptyD<MpcMdLaw> && ControlsWithEmptyD<MpcMdLaw>);
static_assert(ControlsWith<MpcMdLaw, V1, V1> && StepsWith<MpcMdLaw, V1, V1, X2>); // control(r, d), step(r, d, x)
static_assert(!ControlsWith<MpcMdLaw, V1> && !StepsWith<MpcMdLaw, V1, X2>);

// --- StateEstimator: step(y, u) + state() ------------------------------------
static_assert(StateEstimator<KalmanFilter<2, 1, 1, double, 2, 1>, 2, 1, 1, double>); // recursive
static_assert(StateEstimator<Luenberger<2, 1, 1, double>, 2, 1, 1, double>);         // fixed gain
static_assert(StateEstimator<Luenberger<3, 1, 1, double, 1>, 3, 1, 1, double>);      // fixed gain, anchored
static_assert(StateEstimator<MHE<1, 1, 1, 4, double>, 1, 1, 1, double>);
static_assert(StateEstimator<ReducedLuenberger<3, 1, double>, 3, 1, 1, double>); // reduced order
struct ConstantVelocityModel {                                                   // ẋ = [v, 0]
    [[nodiscard]] ColVec<2> operator()(const ColVec<2>& x, const ColVec<1>& /*u*/) const { return ColVec<2>{x[1], 0.0}; }
};
static_assert(!StateEstimator<ExtendedKalmanFilter<ConstantVelocityModel, double>, 2, 1, 1, double>); // update(h, y, R): no step

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

TEST_CASE("fixed-rate PID equals ContinuousPID ticked at the same Ts") {
    constexpr double Ts = 1e-3;

    SUBCASE("PID (discrete vs ContinuousPID at fixed Ts)") {
        design::PIDGains<double> gains{};
        gains.Kp = 2.0;
        gains.Ki = 10.0;
        gains.Kd = 0.05;

        PIDController<double> disc{gains.discretize(Ts)};
        ContinuousPID<double> cont{gains};
        for (int k = 0; k < 50; ++k) {
            const double r = 1.0;
            const double y = 0.1 * k * Ts;
            CHECK(disc.control(r, y) == doctest::Approx(cont.control(r, y, Ts)).epsilon(1e-15));
        }
    }
}

TEST_CASE("fused step(y, u) equals the manual predict/update sequence") {
    constexpr double Ts = 0.1;

    StateSpace<2, 1, 1, double, 2, 1> sys{};
    sys.A = Matrix<2, 2>{{1.0, Ts}, {0.0, 1.0}};
    sys.B = Matrix<2, 1>{{0.5 * Ts * Ts}, {Ts}};
    sys.C = Matrix<1, 2>{{1.0, 0.0}};
    sys.Ts = Ts;

    SUBCASE("KalmanFilter") {
        const auto kd = design::kalman(sys, {.process = Matrix<2, 2>::identity(), .measurement = 0.01});
        REQUIRE(kd.success);

        KalmanFilter<2, 1, 1, double, 2, 1> fused{kd};
        KalmanFilter<2, 1, 1, double, 2, 1> manual{kd};
        for (int k = 0; k < 30; ++k) {
            const ColVec<1> y{0.1 * k};
            const ColVec<1> u{0.5};
            const auto      xf = fused.step(y, u);
            manual.predict(u);
            (void)manual.update(y);
            CHECK(xf(0) == doctest::Approx(manual.state()(0)).epsilon(1e-15));
            CHECK(xf(1) == doctest::Approx(manual.state()(1)).epsilon(1e-15));
        }
    }

    SUBCASE("Luenberger") {
        const auto ld = design::luenberger(sys, {0.5, 0.6});
        REQUIRE(ld.success);
        Luenberger<2, 1, 1, double> fused{ld};
        Luenberger<2, 1, 1, double> manual{ld};
        for (int k = 0; k < 30; ++k) {
            const ColVec<1> y{0.1 * k};
            const ColVec<1> u{0.5};
            const auto      xf = fused.step(y, u);
            manual.predict(u);
            manual.update(y, u);
            CHECK(xf(0) == manual.state()(0));
            CHECK(xf(1) == manual.state()(1));
            CHECK(xf(0) == fused.state()(0));
        }
    }

    SUBCASE("MHE") {
        StateSpace<1, 1, 1, double, 1, 1> tank{};
        tank.A = Matrix<1, 1>{{0.9}};
        tank.B = Matrix<1, 1>{{0.1}};
        tank.C = Matrix<1, 1>{{1.0}};
        tank.Ts = Ts;
        const auto art = design::mhe<4>(tank, {.process = 1e-4, .measurement = 1e-3});
        REQUIRE(art.success);

        MHE est{art};
        for (int k = 0; k < 20; ++k) {
            const ColVec<1> y{0.3 + (0.01 * k)};
            const ColVec<1> u{0.1};
            const auto      xf = est.step(y, u); // no predict/update split: step returns the window's last state
            CHECK(est.state()(0) == doctest::Approx(xf(0)));
        }
    }
}

namespace {

/// Tick two copies of @p law in closed loop on @p sys from x = [1; 0], one through @p basic and one through @p full,
/// alternating step with update then control; the commands are equal at every tick (to rounding:
/// under -ffast-math the folded zeros of the basic call may reassociate the sums).
template<typename Law, typename Plant, typename Measure, typename Basic, typename Full>
void check_basic_equals_full(const Law& law, const Plant& sys, Measure measure, Basic basic, Full full) {
    Law       a = law;
    Law       b = law;
    ColVec<2> xa{1.0, 0.0};
    ColVec<2> xb{1.0, 0.0};
    for (int k = 0; k < 40; ++k) {
        const bool      split = (k % 2) == 1;
        const ColVec<1> ua = basic(a, measure(xa), split);
        const ColVec<1> ub = full(b, measure(xb), split);
        CHECK(ua(0) == doctest::Approx(ub(0)).epsilon(1e-12));
        xa = sys.A * xa + sys.B * ua;
        xb = sys.A * xb + sys.B * ub;
    }
    CHECK(xa(0) != 1.0); // the loop moved
}

} // namespace

TEST_CASE("each LQ law's basic step(r, m) and control(r) equal the feedforward form with x_ref = 0 and u_ref = 0") {
    constexpr double                  Ts = 0.1;
    StateSpace<2, 1, 1, double, 2, 1> sys{};
    sys.A = Matrix<2, 2>{{1.0, Ts}, {0.0, 1.0}};
    sys.B = Matrix<2, 1>{{0.5 * Ts * Ts}, {Ts}};
    sys.C = Matrix<1, 2>{{1.0, 0.0}};
    sys.Ts = Ts;

    const Matrix<2, 2>        Q = Matrix<2, 2>::identity();
    const Matrix<1, 1>        R{{0.1}};
    const Matrix<1, 1>        Qi{{10.0}};
    const design::KalmanNoise noise{.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, .measurement = Matrix<1, 1>{{0.1}}};
    const design::KalmanNoise load_noise{Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.01}}, Matrix<1, 1>{{0.1}}, Matrix<1, 1>{{1.0}}};

    const auto      state = [](const ColVec<2>& x) { return x; };
    const auto      output = [&](const ColVec<2>& x) { return ColVec<1>{(sys.C * x)(0, 0)}; };
    const ColVec<1> r{0.5}; // the controlled outputs' reference
    const ColVec<2> x0{};
    const ColVec<1> u0{};

    // The LQR and LQG regulators: control(), step(m) against control(0, 0), step(0, 0, m).
    const auto basic_reg = [&](auto& law, const auto& m, bool split) {
        if (split) {
            law.predict();
            law.update(m);
            return law.control();
        }
        return law.step(m);
    };
    const auto full_reg = [&](auto& law, const auto& m, bool split) {
        if (split) {
            law.predict();
            law.update(m);
            return law.control(x0, u0);
        }
        return law.step(x0, u0, m);
    };
    // The laws with a reference: control(r), step(r, m) against control(r, 0, 0), step(r, 0, 0, m).
    const auto basic_r = [&](auto& law, const auto& m, bool split) {
        if (split) {
            law.predict();
            law.update(m);
            return law.control(r);
        }
        return law.step(r, m);
    };
    const auto full_r = [&](auto& law, const auto& m, bool split) {
        if (split) {
            law.predict();
            law.update(m);
            return law.control(r, x0, u0);
        }
        return law.step(r, x0, u0, m);
    };

    SUBCASE("LQR") {
        const auto res = design::lqr(sys, {.Q = Q, .R = R});
        REQUIRE(res.success);
        check_basic_equals_full(LQR{res}, sys, state, basic_reg, full_reg);
    }
    SUBCASE("LQG") {
        const auto res = design::lqg(sys, {.Q = Q, .R = R}, noise);
        REQUIRE(res.success);
        check_basic_equals_full(LQG{res}, sys, output, basic_reg, full_reg);
    }
    SUBCASE("LQR with a tracked output") {
        const auto res = design::lqr(sys, {.Q = Q, .R = R}, design::TrackedOutputs{});
        REQUIRE(res.success);
        check_basic_equals_full(LQR{res}, sys, state, basic_r, full_r);
    }
    SUBCASE("LQG with a tracked output") {
        const auto res = design::lqg(sys, {.Q = Q, .R = R}, noise, design::TrackedOutputs{});
        REQUIRE(res.success);
        check_basic_equals_full(LQG{res}, sys, output, basic_r, full_r);
    }
    SUBCASE("LQI") {
        const auto res = design::lqi(sys, {.Q = Q, .R = R, .Q_integral = Qi});
        REQUIRE(res.success);
        check_basic_equals_full(LQI{res}, sys, state, basic_r, full_r);
    }
    SUBCASE("LQGI") {
        const auto res = design::lqgi(sys, {.Q = Q, .R = R, .Q_integral = Qi}, noise);
        REQUIRE(res.success);
        check_basic_equals_full(LQGI{res}, sys, output, basic_r, full_r);
    }
    SUBCASE("OFLQG") {
        const auto res = design::oflqg(sys, {.Q = Q, .R = R}, load_noise, design::InputDisturbance{});
        REQUIRE(res.success);
        check_basic_equals_full(OFLQG{res}, sys, output, basic_r, full_r);
    }
}
