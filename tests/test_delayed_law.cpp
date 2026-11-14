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
#include <vector>

#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/oflqg.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/steady_state_target.hpp"
#include "damp/estimation/kalman.hpp"
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
 * start, ∫ (xᵀQx + vᵀRv + ξᵀQ_integral ξ) dt summed exactly per sub-step (Van Loan on [x; ξ]
 * with [v; e] held). The law reads the plant state and holds its registers itself.
 */
template<size_t NI, typename L>
double continuous_cost(L& law, const design::LQRConfig<2, 1, NI>& w, const ColVec<2>& x0, double xi0, int periods) {
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
        Q(2, 2) = w.Q_integral(0, 0);
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
        const design::LQRConfig<2, 1, 0> w{.Q = Matrix<2, 2>{{10.0, 0.0}, {0.0, 1.0}}, .R = Matrix<1, 1>{{0.1}}};
        const auto                       res = design::lqr(plant(), kTs, InputDelay<2>{kTheta}, w);
        REQUIRE(res.success);
        REQUIRE(res.is_stable());
        LQR             law{res};
        const ColVec<2> x0{{0.3, -0.2}};
        const double    J = continuous_cost<0>(law, w, x0, 0.0, 600);
        const ColVec<4> z0{{0.3, -0.2, 0.0, 0.0}}; // the registers empty
        const double    predicted = (z0.t() * res.S * z0)(0, 0);
        CHECK(J == doctest::Approx(predicted).epsilon(1e-8));
    }

    TEST_CASE("LQI with an InputDelay: the closed loop's continuous cost, integrator included, is the one the design predicts") {
        const design::LQRConfig<2, 1, 1> w{.Q = Matrix<2, 2>{{10.0, 0.0}, {0.0, 1.0}}, .R = Matrix<1, 1>{{0.1}}, .Q_integral = Matrix<1, 1>{{5.0}}};
        const auto                       res = design::lqi(plant(), kTs, InputDelay<2>{kTheta}, w);
        REQUIRE(res.success);
        REQUIRE(res.is_stable());
        LQI             law{res};
        const ColVec<2> x0{{0.3, -0.2}};
        const double    J = continuous_cost<1>(law, w, x0, 0.1, 600);
        const ColVec<5> q0{{0.3, -0.2, 0.0, 0.0, 0.1}}; // [x; registers empty; ξ]
        const double    predicted = (q0.t() * res.S * q0)(0, 0);
        CHECK(J == doctest::Approx(predicted).epsilon(1e-8));
    }

    TEST_CASE("an InputDelay of zero is the design without one, the registers unread") {
        const design::LQRConfig<2, 1, 1> w{.Q = Matrix<2, 2>{{10.0, 0.0}, {0.0, 1.0}}, .R = Matrix<1, 1>{{0.1}}, .Q_integral = Matrix<1, 1>{{5.0}}};
        const auto                       late = design::lqi(plant(), kTs, InputDelay<2>{0.0}, w);
        const auto                       now = design::lqi(plant(), kTs, w);
        REQUIRE(late.success);
        REQUIRE(now.success);
        const double scale = now.K.norm();
        CHECK(damp::abs(late.K(0, 0) - now.K(0, 0)) <= 1e-9 * scale);
        CHECK(damp::abs(late.K(0, 1) - now.K(0, 1)) <= 1e-9 * scale);
        CHECK(damp::abs(late.K(0, 4) - now.K(0, 2)) <= 1e-9 * scale); // the integrator
        CHECK(damp::abs(late.K(0, 2)) <= 1e-9 * scale);
        CHECK(damp::abs(late.K(0, 3)) <= 1e-9 * scale);

        const design::LQRConfig<2, 1, 0> wr{.Q = w.Q, .R = w.R};
        const auto                       lqr_late = design::lqr(plant(), kTs, InputDelay<1>{0.0}, wr);
        const auto                       lqr_now = design::lqr(plant(), kTs, wr);
        REQUIRE(lqr_late.success);
        CHECK((lqr_late.K.template block<1, 2>(0, 0) - lqr_now.K).norm() <= 1e-9 * lqr_now.K.norm());
    }

    TEST_CASE("the law holds what back_calculate recorded, newest first, and a reference point holds u* in them") {
        const design::LQRConfig<2, 1, 1> w{.Q = Matrix<2, 2>{{10.0, 0.0}, {0.0, 1.0}}, .R = Matrix<1, 1>{{0.1}}, .Q_integral = Matrix<1, 1>{{5.0}}};
        const auto                       res = design::lqi(plant(), kTs, InputDelay<3>{0.25}, w);
        REQUIRE(res.success);
        LQI             law{res};
        const auto&     K = law.getK();
        const ColVec<2> x{{0.1, 0.0}};
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
        const ColVec<2> x_ref{{0.2, 0.0}};
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
        static constexpr auto kRes = design::lqi(plant(), 0.01, InputDelay<3>{0.022}, {.Q = Matrix<2, 2>{{10.0, 0.0}, {0.0, 1.0}}, .R = Matrix<1, 1>{{0.1}}, .Q_integral = Matrix<1, 1>{{5.0}}});
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
        const design::LQRConfig<2, 1, 0> w{.Q = Matrix<2, 2>::identity(), .R = Matrix<1, 1>{{0.1}}};
        const auto                       discrete = discretize(plant(), kTs).value();
        CHECK(design::lqr(discrete, kTs, InputDelay<2>{kTheta}, w).status == design::DesignStatus::invalid_sample_time);
        CHECK(design::lqr(plant(), kTs, InputDelay<1>{kTheta}, w).status == design::DesignStatus::invalid_parameter);
        CHECK(design::lqr(plant(), kTs, InputDelay<2>{kTheta}, {.Q = w.Q, .R = Matrix<1, 1>{{0.0}}}).status == design::DesignStatus::bad_weights);
    }
}

namespace {

const design::KalmanNoise<2, 1>  kNoise{.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.1}}, .measurement = Matrix<1, 1>{{1e-4}}};
const design::LQRConfig<2, 1, 0> kW{.Q = Matrix<2, 2>{{10.0, 0.0}, {0.0, 1.0}}, .R = Matrix<1, 1>{{0.1}}};
const design::LQRConfig<2, 1, 1> kWi{.Q = kW.Q, .R = kW.R, .Q_integral = Matrix<1, 1>{{5.0}}};

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
        const auto lqr = design::lqr(noisy(), kTs, InputDelay<2>{kTheta}, kW);
        const auto res = design::lqg(noisy(), kTs, InputDelay<2>{kTheta}, kW, kNoise);
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
        const auto late = design::kalman(noisy(), kTs, InputDelay<1>{0.0}, kNoise);
        REQUIRE(now.success);
        REQUIRE(late.success);
        CHECK((Matrix<2, 2>(late.P.block<2, 2>(0, 0)) - now.P).norm() <= 1e-9 * now.P.norm());
        CHECK((Matrix<2, 1>(late.M.block<2, 1>(0, 0)) - now.M).norm() <= 1e-9 * now.M.norm());
    }

    TEST_CASE("to_ss of an LQG with an InputDelay is the runtime's tick, on the delayed model") {
        const auto res = design::lqg(noisy(), kTs, InputDelay<2>{kTheta}, kW, kNoise);
        REQUIRE(res.success);
        const auto ctrl = res.to_ss();
        static_assert(damp::remove_cvref_t<decltype(ctrl.A)>::rows() == 4);
        const auto model = delayed_model();
        LQG        law{res};
        ColVec<4>  z{{0.1, 0.0, 0.0, 0.0}};
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

    TEST_CASE("LQG and LQR with an InputDelay hold r alike, the target's registers at u_t") {
        const auto lqr = design::lqr(noisy(), kTs, InputDelay<2>{kTheta}, kW, design::TrackedOutputs{});
        const auto lqg = design::lqg(noisy(), kTs, InputDelay<2>{kTheta}, kW, kNoise, design::TrackedOutputs{});
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
        const auto lqi = design::lqi(noisy(), kTs, InputDelay<2>{kTheta}, kWi);
        const auto lqgi = design::lqgi(noisy(), kTs, InputDelay<2>{kTheta}, kWi, kNoise);
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
        const auto at_input = design::oflqg(noisy(), kTs, InputDelay<2>{kTheta}, kW, noise, design::InputDisturbance{});
        REQUIRE(at_input.success);
        CHECK(run(at_input, ColVec<4>{}, 0.4) == doctest::Approx(0.2).epsilon(1e-6));

        // A load on the plant itself acts at once: its channel is the plant's, the registers untouched.
        const auto on_plant = design::oflqg(noisy(), kTs, InputDelay<2>{kTheta}, kW, noise, design::DisturbanceModel<2, 1, 1>{.Bd = noisy().B});
        REQUIRE(on_plant.success);
        const auto gamma = discretize(plant(), kTs).value().B;
        CHECK(run(on_plant, ColVec<4>{{gamma(0, 0), gamma(1, 0), 0.0, 0.0}}, 0.0) == doctest::Approx(0.2).epsilon(1e-6));
    }

    TEST_CASE("a constexpr LQG with an InputDelay deploys in float on constants") {
        static constexpr auto kRes = design::lqg(
            noisy(), 0.01, InputDelay<3>{0.022}, {.Q = Matrix<2, 2>{{10.0, 0.0}, {0.0, 1.0}}, .R = Matrix<1, 1>{{0.1}}},
            {.process = Matrix<2, 2>{{0.01, 0.0}, {0.0, 0.1}}, .measurement = Matrix<1, 1>{{1e-4}}}
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
        CHECK(design::kalman(discrete, kTs, InputDelay<2>{kTheta}, kNoise).status == design::DesignStatus::invalid_sample_time);
        CHECK(design::lqg(noisy(), kTs, InputDelay<1>{kTheta}, kW, kNoise).status == design::DesignStatus::invalid_parameter);
        CHECK(design::lqgi(noisy(), kTs, InputDelay<1>{kTheta}, kWi, kNoise).status == design::DesignStatus::invalid_parameter);
    }
}
