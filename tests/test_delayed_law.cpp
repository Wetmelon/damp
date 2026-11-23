// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// The LQ designs with an InputDelay: the continuous plant, its weights and its noise, the delay,
// and a law that holds or carries the inputs still on their way. For the measured-state laws the
// oracle is the continuous cost itself: the true plant run on the delayed, held input, its cost
// summed exactly over sub-steps the delay's switch falls on, against the Riccati form the design
// predicts. A filter's law on the exact model with no noise is the measured law, tick for tick.

#include <cmath>
#include <cstddef>
#include <type_traits>
#include <vector>

#include "damp/backend.hpp"
#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/controllers/oflqg.hpp"
#include "damp/controllers/state_feedback.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/status.hpp"
#include "damp/design/steady_state_target.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "doctest.h"

using namespace damp;

namespace {

/// An unstable second-order plant (a pole at 1.8 rad/s): position, speed; force in; position out.
constexpr StateSpace<2, 1, 1> plant() {
    StateSpace<2, 1, 1> s{};
    s.A(0, 1) = 1.0;
    s.A(1, 0) = 4.0;
    s.A(1, 1) = -0.4;
    s.B(1, 0) = 1.0;
    s.C(0, 0) = 1.0;
    return s;
}

constexpr double kTs = 0.1;
constexpr int    kSub = 10;  // sub-steps per period
constexpr int    kLate = 14; // the delay in sub-steps: 0.14 s, one period and 0.4 of the next
constexpr double kTheta = 0.14;

/**
 * The continuous cost of the closed loop from x0 and ξ0, the registers empty: the true plant on
 * the input @p law sent kLate sub-steps earlier, ξ ramping on the error sampled at each period's
 * start, ∫ (xᵀQx + vᵀRv + ξᵀ·q_integral·ξ) dt summed exactly per sub-step (Van Loan on [x; ξ]
 * with [v; e] held), q_integral the IntegralAction's Q (unread for NI = 0). The law reads the
 * plant state and holds its delay states itself.
 */
template<size_t NI, typename L>
double continuous_cost(L& law, const design::LQRWeights<2, 1>& w, double q_integral, const ColVec<2>& x0, double xi0, int periods) {
    constexpr size_t      NE = 2 + NI;
    constexpr size_t      NH = 1 + NI;
    const auto            p = plant();
    const double          h = kTs / kSub;
    StateSpace<NE, NH, 1> ext{};
    ext.A.template block<2, 2>(0, 0) = p.A;
    ext.B(1, 0) = 1.0;
    Matrix<NE, NE> Q{};
    Q.template block<2, 2>(0, 0) = w.Q;
    Matrix<NH, NH> R{};
    R(0, 0) = w.R(0, 0);
    if constexpr (NI > 0) {
        ext.B(2, 1) = -1.0; // ξ̇ = −e, e = the position sampled at the period's start
        Q(2, 2) = q_integral;
    }
    const auto cost = design::discretize_lqr_cost(ext.A, ext.B, Q, R, h).value();
    const auto step = discretize(ext, h).value();

    ColVec<NE> s{};
    s(0) = x0(0);
    s(1) = x0(1);
    if constexpr (NI > 0) {
        s(2) = xi0;
        law.integral.xi(0) = xi0;
    }
    std::vector<double> sent;
    double              J = 0.0;
    for (int k = 0; k < periods; ++k) {
        law.update(ColVec<2>{s(0), s(1)});
        if constexpr (NI > 0) {
            sent.push_back(law.control(ColVec<1>{0.0})(0));
        } else {
            sent.push_back(law.control()(0));
        }
        law.back_calculate(ColVec<1>{sent.back()});
        law.predict();
        const double e = s(0);
        for (int i = 0; i < kSub; ++i) {
            const int  n = (k * kSub) + i;
            ColVec<NH> held{};
            held(0) = (n >= kLate) ? sent[static_cast<size_t>((n - kLate) / kSub)] : 0.0;
            if constexpr (NI > 0) {
                held(1) = e;
            }
            J += (s.t() * cost.Q * s)(0, 0) + (2.0 * (s.t() * cost.N * held)(0, 0)) + (held.t() * cost.R * held)(0, 0);
            s = ColVec<NE>((step.A * s) + (step.B * held));
        }
    }
    return J;
}

} // namespace

TEST_SUITE("Delayed laws") {
    TEST_CASE("LQR with an InputDelay: the closed loop's continuous cost is the one the design predicts") {
        const design::LQRWeights<2, 1> w{.Q = Matrix<2, 2>{{10.0, 0.0}, {0.0, 1.0}}, .R = Matrix<1, 1>{{0.1}}};
        const auto                     res = design::state_feedback(plant(), kTs, w, InputDelay<2>{kTheta});
        REQUIRE(res.success);
        REQUIRE(res.is_stable());
        LQR             law{res};
        const ColVec<2> x0{0.3, -0.2};
        const double    J = continuous_cost<0>(law, w, 0.0, x0, 0.0, 600);
        const ColVec<4> z0{0.3, -0.2, 0.0, 0.0}; // the registers empty
        const double    predicted = (z0.t() * res.S * z0)(0, 0);
        CHECK(J == doctest::Approx(predicted).epsilon(1e-8));
    }

    TEST_CASE("LQI with an InputDelay: the closed loop's continuous cost, integrator included, is the one the design predicts") {
        const design::LQRWeights<2, 1> w{.Q = Matrix<2, 2>{{10.0, 0.0}, {0.0, 1.0}}, .R = Matrix<1, 1>{{0.1}}};
        const design::IntegralAction   ia{.Q = 5.0};
        const auto                     res = design::state_feedback(plant(), kTs, w, InputDelay<2>{kTheta}, ia);
        REQUIRE(res.success);
        REQUIRE(res.is_stable());
        LQI             law{res};
        const ColVec<2> x0{0.3, -0.2};
        const double    J = continuous_cost<1>(law, w, ia.Q, x0, 0.1, 600);
        const ColVec<5> q0{0.3, -0.2, 0.0, 0.0, 0.1}; // [x; registers empty; ξ]
        const double    predicted = (q0.t() * res.S * q0)(0, 0);
        CHECK(J == doctest::Approx(predicted).epsilon(1e-8));
    }

    TEST_CASE("an InputDelay of zero is the design without one, the registers unread") {
        const design::LQRWeights<2, 1> w{.Q = Matrix<2, 2>{{10.0, 0.0}, {0.0, 1.0}}, .R = Matrix<1, 1>{{0.1}}};
        const auto                     late = design::state_feedback(plant(), kTs, w, InputDelay<2>{0.0}, design::IntegralAction{.Q = 5.0});
        const auto                     now = design::state_feedback(plant(), kTs, w, design::IntegralAction{.Q = 5.0});
        REQUIRE(late.success);
        REQUIRE(now.success);
        const double scale = now.K.norm();
        CHECK(damp::abs(late.K(0, 0) - now.K(0, 0)) <= 1e-9 * scale);
        CHECK(damp::abs(late.K(0, 1) - now.K(0, 1)) <= 1e-9 * scale);
        CHECK(damp::abs(late.K(0, 4) - now.K(0, 2)) <= 1e-9 * scale); // the integrator
        CHECK(damp::abs(late.K(0, 2)) <= 1e-9 * scale);
        CHECK(damp::abs(late.K(0, 3)) <= 1e-9 * scale);

        const auto lqr_late = design::state_feedback(plant(), kTs, w, InputDelay<1>{0.0});
        const auto lqr_now = design::state_feedback(plant(), kTs, w);
        REQUIRE(lqr_late.success);
        CHECK((lqr_late.K.template block<1, 2>(0, 0) - lqr_now.K).norm() <= 1e-9 * lqr_now.K.norm());
    }

    TEST_CASE("the law holds what back_calculate recorded, newest first, and a reference point holds u* in them") {
        const design::LQRWeights<2, 1> w{.Q = Matrix<2, 2>{{10.0, 0.0}, {0.0, 1.0}}, .R = Matrix<1, 1>{{0.1}}};
        const auto                     res = design::state_feedback(plant(), kTs, w, InputDelay<3>{0.25}, design::IntegralAction{.Q = 5.0});
        REQUIRE(res.success);
        LQI             law{res};
        const auto&     K = law.getK();
        const ColVec<2> x{0.1, 0.0};
        // Three ticks, each command clamped to a value of its own: the registers hold those.
        const double held[3] = {0.5, -0.25, 0.125};
        for (const double v : held) {
            law.update(x);
            static_cast<void>(law.control(ColVec<1>{0.0}));
            law.back_calculate(ColVec<1>{v});
            law.predict();
        }
        law.update(x);
        const double xi = law.integral.xi(0);
        const double u = law.control(ColVec<1>{0.0})(0);
        const double expected = -((K(0, 0) * 0.1) + (K(0, 2) * held[2]) + (K(0, 3) * held[1]) + (K(0, 4) * held[0]) + (K(0, 5) * xi));
        CHECK(u == doctest::Approx(expected).epsilon(1e-12));

        // About a reference point (x*, u*), registers that hold u* are at the reference.
        LQI             fresh{res};
        const ColVec<2> x_ref{0.2, 0.0};
        const ColVec<1> u_ref{0.7};
        for (int k = 0; k < 3; ++k) {
            fresh.update(x_ref);
            static_cast<void>(fresh.control(ColVec<1>{0.2}, x_ref, u_ref));
            fresh.back_calculate(u_ref);
            fresh.predict();
        }
        fresh.integral.xi(0) = 0.0;
        fresh.update(x_ref);
        CHECK(fresh.control(ColVec<1>{0.2}, x_ref, u_ref)(0) == doctest::Approx(0.7).epsilon(1e-12));
    }

    TEST_CASE("a constexpr design with an InputDelay deploys in float on constants") {
        static constexpr auto kRes = design::state_feedback(plant(), 0.01, {.Q = Matrix<2, 2>{{10.0, 0.0}, {0.0, 1.0}}, .R = Matrix<1, 1>{{0.1}}}, InputDelay<3>{0.022}, design::IntegralAction{.Q = 5.0});
        static_assert(kRes.success);
        Law<kRes, float> fast;
        LQI              ref{kRes};
        for (int k = 0; k < 20; ++k) {
            const double p = 0.1 * std::cos(0.3 * k);
            fast.update(ColVec<2, float>{static_cast<float>(p), 0.0f});
            ref.update(ColVec<2>{p, 0.0});
            const float  uf = fast.control(ColVec<1, float>{0.0f})(0);
            const double ud = ref.control(ColVec<1>{0.0})(0);
            CHECK(static_cast<double>(uf) == doctest::Approx(ud).epsilon(1e-4));
            fast.back_calculate(ColVec<1, float>{uf});
            ref.back_calculate(ColVec<1>{ud});
            fast.predict();
            ref.predict();
        }
    }

    TEST_CASE("refusals: a discrete plant, a delay past the registers, an input cost that is not positive definite") {
        const design::LQRWeights<2, 1> w{.Q = Matrix<2, 2>::identity(), .R = Matrix<1, 1>{{0.1}}};
        const auto                     discrete = discretize(plant(), kTs).value();
        CHECK(design::state_feedback(discrete, kTs, w, InputDelay<2>{kTheta}).status == design::DesignStatus::invalid_sample_time);
        CHECK(design::state_feedback(plant(), kTs, w, InputDelay<1>{kTheta}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::state_feedback(plant(), kTs, {.Q = w.Q, .R = Matrix<1, 1>{{0.0}}}, InputDelay<2>{kTheta}).status == design::DesignStatus::bad_weights);
    }
}

namespace {

const design::KalmanNoise<2, 1> kNoise{.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.1}}, .measurement = Matrix<1, 1>{{1e-4}}};
const design::LQRWeights<2, 1>  kW{.Q = Matrix<2, 2>{{10.0, 0.0}, {0.0, 1.0}}, .R = Matrix<1, 1>{{0.1}}};
const design::IntegralAction    kIntegral{.Q = 5.0};

/// The plant with its noise channels: process noise on both states (G = I), one on the reading.
constexpr StateSpace<2, 1, 1, double, 2, 1> noisy() {
    const auto p = plant();
    return {.A = p.A, .B = p.B, .C = p.C, .D = p.D, .G = Matrix<2, 2>::identity(), .H = Matrix<1, 1>::identity()};
}

/// The delayed plant at the samples, exactly: z = [x; v(k−1); v(k−2)].
StateSpace<4, 1, 1, double, 2, 1> delayed_model() { return discretize(noisy(), kTs, InputDelay<2>{kTheta}).value(); }

/**
 * Close @p measured (fed the plant state) and @p filtered (fed y) on the exact delayed model, each
 * on its own copy, both from rest at @p r, no noise: the filter's estimate is then the state, so
 * the commands agree. @return the largest gap between them over the run, relative to the largest command
 */
template<typename M, typename F, typename Control>
double gap(M& measured, F& filtered, Control control, ColVec<1>& y_end) {
    const auto model = delayed_model();
    ColVec<4>  zm{};
    ColVec<4>  zf{};
    double     worst = 0.0;
    double     scale = 0.0;
    for (int k = 0; k < 400; ++k) {
        measured.update(ColVec<2>{zm(0), zm(1)});
        filtered.update(ColVec<1>(model.C * zf));
        const ColVec<1> um = control(measured);
        const ColVec<1> uf = control(filtered);
        measured.back_calculate(um);
        filtered.back_calculate(uf);
        measured.predict();
        filtered.predict();
        worst = std::fmax(worst, std::fabs(um(0) - uf(0)));
        scale = std::fmax(scale, std::fabs(um(0)));
        zm = ColVec<4>((model.A * zm) + (model.B * um));
        zf = ColVec<4>((model.A * zf) + (model.B * uf));
    }
    y_end = ColVec<1>(model.C * zf);
    return worst / scale;
}

} // namespace

TEST_SUITE("Delayed laws with a filter") {
    TEST_CASE("LQG with an InputDelay: the delayed LQR's gain, a filter on the delayed model that carries the registers exactly") {
        const auto lqr = design::state_feedback(noisy(), kTs, kW, InputDelay<2>{kTheta});
        const auto res = design::state_feedback(noisy(), kTs, kW, kNoise, InputDelay<2>{kTheta});
        REQUIRE(res.success);
        CHECK((res.K - lqr.K).norm() == 0.0);
        const auto model = delayed_model();
        CHECK((res.kalman.sys.A - model.A).norm() == 0.0);
        CHECK((res.kalman.sys.B - model.B).norm() == 0.0);
        for (size_t i = 2; i < 4; ++i) { // the registers hold the inputs applied: no error, no reading moves them
            CHECK(damp::abs(res.kalman.P(i, i)) <= 1e-10 * res.kalman.P.norm());
            CHECK(damp::abs(res.kalman.M(i, 0)) <= 1e-10 * res.kalman.M.norm());
        }

        // With no delay the plant's block is the filter without one: the noise is the plant's.
        const auto now = design::kalman(noisy(), kTs, kNoise);
        const auto late = design::kalman(noisy(), kTs, kNoise, InputDelay<1>{0.0});
        REQUIRE(now.success);
        REQUIRE(late.success);
        CHECK((Matrix<2, 2>(late.P.block<2, 2>(0, 0)) - now.P).norm() <= 1e-9 * now.P.norm());
        CHECK((Matrix<2, 1>(late.M.block<2, 1>(0, 0)) - now.M).norm() <= 1e-9 * now.M.norm());
    }

    TEST_CASE("to_ss of an LQG with an InputDelay is the runtime's tick, on the delayed model") {
        const auto res = design::state_feedback(noisy(), kTs, kW, kNoise, InputDelay<2>{kTheta});
        REQUIRE(res.success);
        const auto ctrl = res.to_ss();
        static_assert(damp::remove_cvref_t<decltype(ctrl.A)>::rows() == 4);
        const auto model = delayed_model();
        LQG        law{res};
        ColVec<4>  z{0.1, 0.0, 0.0, 0.0};
        ColVec<4>  zc = z;
        ColVec<4>  xc{};
        double     worst = 0.0;
        double     scale = 0.0;
        for (int k = 0; k < 200; ++k) {
            law.update(model.C * z);
            const double u = law.control()(0);
            law.back_calculate(ColVec<1>{u});
            law.predict();
            const ColVec<1> yc(model.C * zc);
            const double    uc = ((ctrl.C * xc) + (ctrl.D * yc))(0, 0);
            xc = (ctrl.A * xc) + (ctrl.B * yc);
            worst = std::fmax(worst, std::fabs(u - uc));
            scale = std::fmax(scale, std::fabs(u));
            z = (model.A * z) + (model.B * ColVec<1>{u});
            zc = (model.A * zc) + (model.B * ColVec<1>{uc});
        }
        CHECK(worst <= 1e-9 * scale);
        CHECK(z.norm() < 1e-4);
    }

    TEST_CASE("to_ss of an OFLQG with an InputDelay is the runtime's tick at r = 0, d̂ in the gain") {
        const design::KalmanNoise<2, 1, double, 1> noise{.process = kNoise.process, .measurement = kNoise.measurement, .disturbance = Matrix<1, 1>{{1.0}}};
        const auto                                 res = design::state_feedback(noisy(), kTs, kW, noise, design::InputDisturbance{}, InputDelay<2>{kTheta});
        REQUIRE(res.success);
        const auto ctrl = res.to_ss();
        static_assert(damp::remove_cvref_t<decltype(ctrl.A)>::rows() == 5); // [x; v(k−1); v(k−2); d]
        const auto model = delayed_model();
        OFLQG      law{res};
        ColVec<4>  z{0.1, 0.0, 0.0, 0.0};
        ColVec<4>  zc = z;
        ColVec<5>  xc{};
        double     worst = 0.0;
        double     scale = 0.0;
        for (int k = 0; k < 200; ++k) {
            law.update(model.C * z);
            const double u = law.control(ColVec<1>{})(0);
            law.predict();
            const ColVec<1> yc(model.C * zc);
            const double    uc = ((ctrl.C * xc) + (ctrl.D * yc))(0, 0);
            xc = (ctrl.A * xc) + (ctrl.B * yc);
            worst = std::fmax(worst, std::fabs(u - uc));
            scale = std::fmax(scale, std::fabs(u));
            z = (model.A * z) + (model.B * ColVec<1>{u});
            zc = (model.A * zc) + (model.B * ColVec<1>{uc});
        }
        CHECK(worst <= 1e-9 * scale);
        CHECK(z.norm() < 1e-4);

        // Closed on a plant a period later than the design's: the loop the block makes is the runtime's there too.
        const auto later = discretize(noisy(), kTs, InputDelay<3>{kTheta + kTs}).value();
        OFLQG      late_law{res};
        ColVec<5>  w{0.1, 0.0, 0.0, 0.0, 0.0};
        ColVec<5>  wc = w;
        xc = ColVec<5>{};
        worst = 0.0;
        for (int k = 0; k < 100; ++k) {
            late_law.update(later.C * w);
            const double u = late_law.control(ColVec<1>{})(0);
            late_law.predict();
            const ColVec<1> yc(later.C * wc);
            const double    uc = ((ctrl.C * xc) + (ctrl.D * yc))(0, 0);
            xc = (ctrl.A * xc) + (ctrl.B * yc);
            worst = std::fmax(worst, std::fabs(u - uc));
            w = (later.A * w) + (later.B * ColVec<1>{u});
            wc = (later.A * wc) + (later.B * ColVec<1>{uc});
        }
        CHECK(worst <= 1e-9 * scale);
    }

    TEST_CASE("LQG and LQR with an InputDelay hold r alike, the target's registers at u_t") {
        const auto lqr = design::state_feedback(noisy(), kTs, kW, design::TrackedOutputs{}, InputDelay<2>{kTheta});
        const auto lqg = design::state_feedback(noisy(), kTs, kW, kNoise, design::TrackedOutputs{}, InputDelay<2>{kTheta});
        REQUIRE(lqr.success);
        REQUIRE(lqg.success);
        LQR       measured{lqr};
        LQG       filtered{lqg};
        ColVec<1> y{};
        CHECK(gap(measured, filtered, [](auto& law) { return law.control(ColVec<1>{0.5}); }, y) <= 1e-9);
        CHECK(y(0) == doctest::Approx(0.5).epsilon(1e-6));

        const auto target = design::steady_state_target(noisy(), kTs, InputDelay<2>{kTheta});
        REQUIRE(target.success);
        CHECK(target.target_r(2, 0) == doctest::Approx(target.target_r(4, 0)).epsilon(1e-12)); // [x_t; v(k−1); v(k−2); u_t]
        CHECK(target.target_r(3, 0) == doctest::Approx(target.target_r(4, 0)).epsilon(1e-12));
    }

    TEST_CASE("LQGI and LQI with an InputDelay integrate alike") {
        const auto lqi = design::state_feedback(noisy(), kTs, kW, InputDelay<2>{kTheta}, kIntegral);
        const auto lqgi = design::state_feedback(noisy(), kTs, kW, kNoise, InputDelay<2>{kTheta}, kIntegral);
        REQUIRE(lqi.success);
        REQUIRE(lqgi.success);
        CHECK((lqgi.K - lqi.K).norm() == 0.0);
        LQI       measured{lqi};
        LQGI      filtered{lqgi};
        ColVec<1> y{};
        CHECK(gap(measured, filtered, [](auto& law) { return law.control(ColVec<1>{0.3}); }, y) <= 1e-9);
        CHECK(y(0) == doctest::Approx(0.3).epsilon(1e-6));
    }

    TEST_CASE("OFLQG with an InputDelay holds r against a constant disturbance, through the delay or on the plant at once") {
        const design::KalmanNoise<2, 1, double, 1> noise{.process = kNoise.process, .measurement = kNoise.measurement, .disturbance = Matrix<1, 1>{{1.0}}};
        const auto                                 model = delayed_model();
        const auto                                 run = [&](const auto& res, const ColVec<4>& on_plant, double at_input) {
            OFLQG     law{res};
            ColVec<4> z{};
            for (int k = 0; k < 600; ++k) {
                law.update(ColVec<1>(model.C * z));
                const ColVec<1> u = law.control(ColVec<1>{0.2});
                law.back_calculate(u);
                law.predict();
                z = ColVec<4>((model.A * z) + (model.B * ColVec<1>{u(0) + at_input}) + (on_plant * 0.4));
            }
            return (model.C * z)(0, 0);
        };

        // An InputDisturbance adds to u, so it reaches the plant through the delay.
        const auto at_input = design::state_feedback(noisy(), kTs, kW, noise, design::InputDisturbance{}, InputDelay<2>{kTheta});
        REQUIRE(at_input.success);
        CHECK(run(at_input, ColVec<4>{}, 0.4) == doctest::Approx(0.2).epsilon(1e-6));

        // A load on the plant itself acts at once: its channel is the plant's, the registers untouched.
        const auto on_plant = design::state_feedback(noisy(), kTs, kW, noise, design::DisturbanceModel{.Bd = noisy().B}, InputDelay<2>{kTheta});
        REQUIRE(on_plant.success);
        const auto gamma = discretize(plant(), kTs).value().B;
        CHECK(run(on_plant, ColVec<4>{gamma(0, 0), gamma(1, 0), 0.0, 0.0}, 0.0) == doctest::Approx(0.2).epsilon(1e-6));
    }

    TEST_CASE("a constexpr LQG with an InputDelay deploys in float on constants") {
        static constexpr auto kRes = design::state_feedback(
            noisy(), 0.01, {.Q = Matrix<2, 2>{{10.0, 0.0}, {0.0, 1.0}}, .R = Matrix<1, 1>{{0.1}}},
            design::KalmanNoise{.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.1}}, .measurement = Matrix<1, 1>{{1e-4}}},
            InputDelay<3>{0.022}
        );
        static_assert(kRes.success);
        Law<kRes, float> fast;
        LQG              ref{kRes};
        for (int k = 0; k < 20; ++k) {
            const double y = 0.1 * std::cos(0.3 * k);
            fast.update(ColVec<1, float>{static_cast<float>(y)});
            ref.update(ColVec<1>{y});
            const float  uf = fast.control()(0);
            const double ud = ref.control()(0);
            CHECK(static_cast<double>(uf) == doctest::Approx(ud).epsilon(1e-4));
            fast.back_calculate(ColVec<1, float>{uf});
            ref.back_calculate(ColVec<1>{ud});
            fast.predict();
            ref.predict();
        }
    }

    TEST_CASE("refusals with a filter: a discrete plant, a delay past the registers") {
        const auto discrete = discretize(noisy(), kTs).value();
        CHECK(design::kalman(discrete, kTs, kNoise, InputDelay<2>{kTheta}).status == design::DesignStatus::invalid_sample_time);
        CHECK(design::state_feedback(noisy(), kTs, kW, kNoise, InputDelay<1>{kTheta}).status == design::DesignStatus::invalid_parameter);
        CHECK(design::state_feedback(noisy(), kTs, kW, kNoise, InputDelay<1>{kTheta}, kIntegral).status == design::DesignStatus::invalid_parameter);
    }
}

namespace {
template<typename R>
concept has_integral_gain = requires(const R& r) { r.integral_gain(); };
template<typename R>
concept has_delay_gain = requires(const R& r) { r.delay_gain(); };
} // namespace

TEST_CASE("state_gain, delay_gain and integral_gain are K's columns on x, the delay states and xi") {
    constexpr design::LQRWeights     w{.Q = Matrix<2, 2>{{10.0, 0.0}, {0.0, 1.0}}, .R = {0.1}};
    constexpr design::IntegralAction ia{.Q = 5.0};
    constexpr auto                   res = design::state_feedback(plant(), kTs, w, InputDelay<2>{kTheta}, ia); // K on [x(2); v(k-1), v(k-2); xi]
    static_assert(res.success);
    static_assert(res.state_gain()(0, 0) == res.K(0, 0) && res.state_gain()(0, 1) == res.K(0, 1));
    static_assert(res.delay_gain()(0, 0) == res.K(0, 2) && res.delay_gain()(0, 1) == res.K(0, 3));
    static_assert(res.integral_gain()(0, 0) == res.K(0, 4));

    //! A part the law does not have is not there: no delay states without a delay, no integrators on an LQR.
    constexpr auto lqr = design::state_feedback(plant(), kTs, design::LQRWeights{.Q = Matrix<2, 2>{{10.0, 0.0}, {0.0, 1.0}}, .R = {0.1}});
    static_assert(lqr.success && lqr.state_gain() == lqr.K);
    static_assert(!has_integral_gain<decltype(lqr)> && !has_delay_gain<decltype(lqr)>);
    constexpr auto undelayed = design::state_feedback(plant(), kTs, w, ia);
    static_assert(undelayed.success && !has_delay_gain<decltype(undelayed)> && undelayed.integral_gain()(0, 0) == undelayed.K(0, 2));
    CHECK(res.integral_gain()(0, 0) != 0.0);
}

TEST_CASE("the designs take their features as trailing options, in any order") {
    constexpr design::LQRWeights lq{.Q = Matrix<2, 2>{{10.0, 0.0}, {0.0, 1.0}}, .R = {0.1}};
    constexpr InputDelay<2>      late{kTheta};

    //! The integrators listed and the delay, either way round, are the same design.
    static_assert(design::state_feedback(plant(), kTs, lq, design::IntegralAction{.Q = 5.0, .outputs = {0}}, late).K == design::state_feedback(plant(), kTs, lq, late, design::IntegralAction{.Q = 5.0, .outputs = {0}}).K);
    static_assert(design::state_feedback(plant(), kTs, lq, late, design::IntegralAction{.Q = 5.0}).K == design::state_feedback(plant(), kTs, lq, late, design::IntegralAction{.Q = 5.0, .outputs = {0}}).K); // every output, by default

    //! The disturbance and the delay, of the filter, its target and the offset-free law.
    constexpr design::KalmanNoise<2, 1, double, 1> noise{.process = Matrix<2, 2>{{1e-3, 0.0}, {0.0, 1e-2}}, .measurement = 1e-3, .disturbance = 1e-2};
    static_assert(design::kalman(plant(), kTs, noise, design::InputDisturbance{}, late).M == design::kalman(plant(), kTs, noise, late, design::InputDisturbance{}).M);
    static_assert(design::steady_state_target(plant(), kTs, design::InputDisturbance{}, late).target_d == design::steady_state_target(plant(), kTs, design::InputDisturbance{}, late).target_d);
    constexpr auto ab = design::state_feedback(plant(), kTs, lq, noise, design::InputDisturbance{}, late);
    constexpr auto ba = design::state_feedback(plant(), kTs, lq, noise, late, design::InputDisturbance{});
    static_assert(ab.success && ab.K == ba.K && ab.kalman.M == ba.kalman.M);

    //! One result type, whatever the options: the regulator's has no integrators and no filter.
    static_assert(decltype(design::state_feedback(plant(), kTs, lq))::NK == 2 && decltype(design::state_feedback(plant(), kTs, lq, late))::NK == 4); // K on x, then on x and its 2 delay states
    static_assert(design::state_feedback(plant(), kTs, lq, late).delay_gain()(0, 0) == design::state_feedback(plant(), kTs, lq, late).K(0, 2));

    //! The MPC's limits on an output and a state, its terminal cost and its QP budget, in any order; bare braces are the default weights.
    constexpr auto                 sysd = discretize(plant(), kTs).value();
    constexpr design::MPCWeights   weights{.Qy = {20.0}, .Rdu = {0.1}};
    constexpr design::InputLimits  limits{.u_lim = {.lower = -1.0, .upper = 1.0}};
    constexpr design::OutputLimits cap{.outputs = {0}, .lim = {.upper = 1.0}, .upper_ecr = {0.0}};
    constexpr design::StateLimits  slow{.states = {1}, .lim = {.lower = -0.5, .upper = 0.5}};
    constexpr auto                 bounded = design::mpc<6, 6>(sysd, weights, limits, cap, slow, design::DareCost{}, design::QPIterations{30});
    constexpr auto                 reordered = design::mpc<6, 6>(sysd, weights, limits, design::QPIterations{30}, design::DareCost{}, slow, cap);
    static_assert(std::is_same_v<decltype(bounded), decltype(reordered)> && bounded.success && reordered.mpc.H == bounded.mpc.H && reordered.mpc.A_con == bounded.mpc.A_con && reordered.mpc.b_bound == bounded.mpc.b_bound);
    static_assert(reordered.mpc.max_qp_iterations == 30 && reordered.mpc.P == bounded.mpc.P && decltype(bounded.mpc)::NI == (4 * 6) + (2 * 6) + (2 * 6));
    static_assert(design::mpc<4, 2>(sysd, {}).success && design::mpc<4, 2>(plant(), kTs, {}).success); // on a discrete plant, and on the plant at Ts
}
