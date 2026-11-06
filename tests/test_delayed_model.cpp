// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file test_delayed_model.cpp
 * @brief Laws on a plant with an input delay of at most one sample, and the applied input they
 *        predict with.
 *
 * The delay is part of the plant model: discretize(plant, Ts, InputDelay{θ}) appends one register
 * state v(k−1), checked against the block's closed-form response. A law designed on that model
 * (per-step weights, none on the register) is, for one whole sample, the delay-free law on the
 * state predicted one sample ahead, K₀·[Φ Γ], with the delay-free closed-loop pole; its Kalman
 * filter's register state is the input back_calculate recorded; offset-free LQG and MPC on the
 * model hold the true delayed plant at r under a load. Without a delay: the Van Loan noise of the
 * offset-free LQG's disturbance augment, and back_calculate's anti-windup against recorded runs.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <type_traits>
#include <utility>

#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/controllers/oflqg.hpp"
#include "damp/controllers/ofmpc.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/riccati.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "doctest.h"

using namespace damp;

namespace {

/// The heater of examples/control/heater_ofmpc: K = 6 °C/W, τ = 30 s, sampled at Ts = 2 s here.
constexpr double kA = -1.0 / 30.0; // [1/s] the block's pole
constexpr double kB = 6.0 / 30.0;  // [°C/(W·s)] the heater's gain on it
constexpr double kTs = 2.0;
constexpr double kFrac = 1.4; // [s] a delay within one sample; θ = kTs is one whole sample

constexpr StateSpace<1, 1, 1, double, 1, 1> heater_block() {
    return StateSpace<1, 1, 1, double, 1, 1>{.A = {{kA}}, .B = {{kB}}, .C = {{1.0}}}; // G = H = 1: the noise on the block and the reading
}

/// The heater sampled at kTs without a delay: the model of the delay-blind designs.
constexpr auto kPlain = discretize(heater_block(), kTs);

/**
 * The heater with an input delay θ ≤ Ts as the true plant, from the block's closed-form response
 * (no discretize): within a period it holds v(k−1) for θ, then v(k), so
 * x(k+1) = e^{a·Ts}·x + Γ₀·v(k) + Γ₁·v(k−1), Γ₀ = (e^{a(Ts−θ)} − 1)·b/a, Γ₁ = e^{a(Ts−θ)}·(e^{aθ} − 1)·b/a,
 * and the register holds v(k−1). State [x; v(k−1)].
 */
StateSpace<2, 1, 1> delayed_heater(double theta) {
    const double        tail = std::exp(kA * (kTs - theta));
    StateSpace<2, 1, 1> m{};
    m.A = Matrix<2, 2>{{std::exp(kA * kTs), tail * (std::exp(kA * theta) - 1.0) * kB / kA}, {0.0, 0.0}};
    m.B = Matrix<2, 1>{{(tail - 1.0) * kB / kA}, {1.0}};
    m.C = Matrix<1, 2>{{1.0, 0.0}};
    m.Ts = kTs;
    return m;
}

/// Any discrete model as the true plant, driven by u + w.
template<size_t N>
struct ModelPlant {
    StateSpace<N, 1, 1> m{};
    ColVec<N>           x{};
    double              y() const { return (m.C * x)(0, 0); }
    void                step(double u, double w) { x = m.A * x + m.B * ColVec<1>{u + w}; }
};

/// Largest |a − b| over two matrices of one shape.
template<size_t R, size_t C>
double max_diff(const Matrix<R, C>& a, const Matrix<R, C>& b) {
    double worst = 0.0;
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            worst = std::max(worst, std::abs(a(i, j) - b(i, j)));
        }
    }
    return worst;
}

/// The delay-free cost per step, and the same on the delayed model: the temperature weighted, the register not.
constexpr design::LQRWeights<1, 1> kLq{.Q = 1.0, .R = 0.02};
constexpr design::LQRWeights<2, 1> kLqModel{.Q = {{1.0, 0.0}, {0.0, 0.0}}, .R = 0.02};

/// The offset-free designs' noise per step: model error on the block, the reading, the heat loss's walk.
constexpr design::KalmanNoise<1, 1, double, 1> kNoise{.process = 0.02, .measurement = 0.01, .disturbance = 0.1};

/// The scalar DARE's gain in closed form: P = q + φ²P − (φγP)²/(r + γ²P) is γ²P² + (r(1 − φ²) − qγ²)P − qr = 0, K₀ = φγP/(r + γ²P).
double scalar_lqr_gain(double phi, double gam, double q, double r) {
    const double b = (r * (1.0 - (phi * phi))) - (q * gam * gam);
    const double P = (-b + std::sqrt((b * b) + (4.0 * gam * gam * q * r))) / (2.0 * gam * gam);
    return phi * gam * P / (r + (gam * gam * P));
}

/// A setpoint run on the true plant: the last output, and the largest excursion past r before the load.
struct SetpointRun {
    double y_end = 0.0;
    double overshoot = 0.0;
};

/// Run @p law to @p r from rest on @p truth through the 0 … 40 W clamp, a 2 W heat loss from halfway.
template<typename Law>
SetpointRun setpoint_run(Law& law, const StateSpace<2, 1, 1>& truth, double r) {
    ModelPlant<2> plant{truth};
    SetpointRun   out{};
    for (int k = 0; k < 400; ++k) {
        const double u = std::clamp(law.step(ColVec<1>{r}, ColVec<1>{plant.y()})[0], 0.0, 40.0);
        law.back_calculate(ColVec<1>{u});
        plant.step(u, (k >= 200) ? -2.0 : 0.0);
        if (k < 200) {
            out.overshoot = std::max(out.overshoot, plant.y() - r);
        }
    }
    out.y_end = plant.y();
    return out;
}

/// Lightly damped mass-spring (ω = 2 rad/s, ζ = 0.1): position, speed; force in; position out; noise on both states.
constexpr StateSpace<2, 1, 1, double, 2, 1> mass_spring() {
    StateSpace<2, 1, 1, double, 2, 1> s{};
    s.A(0, 1) = 1.0;
    s.A(1, 0) = -4.0;
    s.A(1, 1) = -0.4;
    s.B(1, 0) = 1.0;
    s.C(0, 0) = 1.0;
    return s; // G = I, H = I by default
}

constexpr double   kTsSpring = 0.01; // the mass-spring's sampling time [s]
const auto         kPlantD = discretize(mass_spring(), kTsSpring).value();
const Matrix<2, 2> kQ{{100.0, 0.0}, {0.0, 1.0}};
const Matrix<1, 1> kR{{0.01}};
const Matrix<1, 1> kQint{{50.0 / (kTsSpring * kTsSpring)}}; // on ∫e
const Matrix<2, 2> kQn{{1e-4, 0.0}, {0.0, 1e-2}};
const Matrix<1, 1> kRn{{1e-4}};

} // namespace

TEST_SUITE("Input delay: a model of the plant") {
    TEST_CASE("InputDelay: one register state, for no delay, a fraction of a sample or a whole sample") {
        static_assert(std::is_same_v<decltype(InputDelay{kFrac}), InputDelay<double>>); // the scalar from the seconds
        constexpr auto none = discretize(heater_block(), kTs, InputDelay{0.0});
        constexpr auto frac = discretize(heater_block(), kTs, InputDelay{kFrac});
        constexpr auto whole = discretize(heater_block(), kTs, InputDelay{kTs});
        static_assert(none.has_value() && frac.has_value() && whole.has_value());
        static_assert(std::is_same_v<decltype(none), decltype(whole)>); // NX + NU states for any delay up to a sample

        // against the closed-form response, the register loading the command in each
        for (const double theta : {0.0, kFrac, kTs}) {
            CAPTURE(theta);
            const auto got = *discretize(heater_block(), kTs, InputDelay{theta});
            const auto ref = delayed_heater(theta);
            CHECK(got.A(0, 0) == doctest::Approx(ref.A(0, 0)).epsilon(1e-12));
            CHECK(got.A(0, 1) == doctest::Approx(ref.A(0, 1)).epsilon(1e-12));
            CHECK(got.B(0, 0) == doctest::Approx(ref.B(0, 0)).epsilon(1e-12));
            CHECK(got.A(1, 0) == 0.0);
            CHECK(got.A(1, 1) == 0.0);
            CHECK(got.B(1, 0) == 1.0);
            CHECK(got.C(0, 0) == 1.0);
            CHECK(got.C(0, 1) == 0.0);
            // the heat of one period is the same, split at θ between the two commands
            CHECK(got.A(0, 1) + got.B(0, 0) == doctest::Approx(kPlain->B(0, 0)).epsilon(1e-12));
        }

        // within a few eps of Ts the delay is the whole sample; beyond it, refused
        constexpr double eps = std::numeric_limits<double>::epsilon();
        const auto       snapped = discretize(heater_block(), kTs, InputDelay{kTs * (1.0 + (4.0 * eps))});
        REQUIRE(snapped.has_value());
        CHECK(snapped->B(0, 0) == 0.0);
        for (const double theta : {kTs * (1.0 + 1e-9), 1.5 * kTs, std::numeric_limits<double>::infinity(), 1e300}) {
            CAPTURE(theta);
            CHECK_FALSE(discretize(heater_block(), kTs, InputDelay{theta}).has_value());
        }
    }

    TEST_CASE("measured state, one whole sample: the law on the model is the delay-free law on the predicted state") {
        // x(k+1) = Φx + Γv(k−1) is known when v(k) is chosen, and v(k) first acts on x(k+2): the
        // cost from k is a constant plus the delay-free cost of the predicted state, so the law on
        // [x; v(k−1)] is K₀·[Φ Γ], and its poles are the delay-free pole and one at z = 0.
        const double phi = std::exp(kA * kTs);
        const double gam = (phi - 1.0) * kB / kA;
        const double k0 = scalar_lqr_gain(phi, gam, kLq.Q(0, 0), kLq.R(0, 0));
        const double pole = phi - (gam * k0); // the delay-free closed loop

        constexpr auto model = *discretize(heater_block(), kTs, InputDelay{kTs});
        constexpr auto res = design::lqr(model, kLqModel);
        static_assert(res.success);
        CHECK(res.K(0, 0) == doctest::Approx(k0 * phi).epsilon(1e-10));
        CHECK(res.K(0, 1) == doctest::Approx(k0 * gam).epsilon(1e-10));
        const bool first_small = damp::abs(res.e[0]) < damp::abs(res.e[1]);
        const auto at_zero = first_small ? res.e[0] : res.e[1];
        const auto moving = first_small ? res.e[1] : res.e[0];
        CHECK(damp::abs(at_zero) <= 1e-10);
        CHECK(moving.real() == doctest::Approx(pole).epsilon(1e-10));
        CHECK(damp::abs(moving.imag()) <= 1e-12);

        // On the plant with the true delay, the law measures [x; u(k−1)], the last input applied:
        // the first period runs open loop on the input sent before (none), then the command meets
        // the state it was computed for and the temperature decays at the delay-free pole.
        LQR           law{res};
        ModelPlant<2> plant{delayed_heater(kTs), ColVec<2>{10.0, 0.0}};
        double        u_prev = 0.0;
        double        prev = plant.y();
        for (int k = 0; k < 20; ++k) {
            law.update(ColVec<2>{plant.y(), u_prev});
            const double u = law.control()[0];
            plant.step(u, 0.0);
            u_prev = u;
            const double expect = (k == 0) ? phi * prev : pole * prev;
            CHECK(plant.y() == doctest::Approx(expect).epsilon(1e-9));
            prev = plant.y();
        }
        CHECK(std::abs(prev) < 1e-3 * 10.0);

        // The delay-blind gain K₀ on the same plant: x(k+1) = Φx(k) − ΓK₀x(k−1), poles of modulus
        // √(ΓK₀) (complex here: Φ² < 4ΓK₀), a ring that decays by √(ΓK₀) per step instead of the pole above.
        REQUIRE(phi * phi < 4.0 * gam * k0);
        const auto blind = design::lqr(*kPlain, kLq);
        REQUIRE(blind.success);
        CHECK(blind.K(0, 0) == doctest::Approx(k0).epsilon(1e-10));
        const double  rho = std::sqrt(gam * k0);
        const double  turn = std::acos(phi / (2.0 * rho)); // the poles' angle per step
        LQR           blind_law{blind};
        ModelPlant<2> blind_plant{delayed_heater(kTs), ColVec<2>{10.0, 0.0}};
        double        x_km1 = 0.0; // x(k−1): the input before the first is none
        double        x_k = blind_plant.y();
        int           sign_changes = 0;
        double        last_turn = 0.0; // the largest |x| over the last full turn
        for (int k = 0; k < 20; ++k) {
            blind_plant.step(blind_law.step(ColVec<1>{blind_plant.y()})[0], 0.0);
            const double expect = (phi * x_k) - ((k == 0) ? 0.0 : gam * k0 * x_km1);
            CHECK(blind_plant.y() == doctest::Approx(expect).epsilon(1e-9));
            sign_changes += ((blind_plant.y() > 0.0) != (x_k > 0.0)) ? 1 : 0;
            x_km1 = x_k;
            x_k = blind_plant.y();
            if (static_cast<double>(19 - k) * turn < 2.0 * std::numbers::pi) {
                last_turn = std::max(last_turn, std::abs(x_k));
            }
        }
        CHECK(sign_changes >= static_cast<int>(20.0 * turn / std::numbers::pi)); // the phase turns 20·turn, a sign change each π
        CHECK(last_turn > 0.5 * 10.0 * std::pow(rho, 20));                       // the envelope from x(0) = 10, where the law above is at 10·pole²⁰
        MESSAGE("delay-blind: ", sign_changes, " sign changes in 20 steps, |x| ", last_turn, " over the last turn; on the model: ", std::abs(prev));

        // the same on a two-state plant: dlqr's gain on the plant, times [Φ Γ]
        const auto lq0 = design::discrete_lqr(kPlantD.A, kPlantD.B, kQ, kR);
        REQUIRE(lq0.success);
        Matrix<3, 3> Q3{};
        Q3.block<2, 2>(0, 0) = kQ; // on the plant states, not the register
        const auto spring = design::lqr(discretize(kPlantD, kTsSpring, InputDelay{kTsSpring}).value(), {.Q = Q3, .R = kR});
        REQUIRE(spring.success);
        Matrix<1, 3> K_pred{};
        K_pred.block<1, 2>(0, 0) = lq0.K * kPlantD.A;
        K_pred(0, 2) = (lq0.K * kPlantD.B)(0, 0);
        CHECK((spring.K - K_pred).norm() < 1e-9 * K_pred.norm());
    }

    TEST_CASE("the Kalman filter's register state is the input back_calculate recorded, else the last command") {
        // The register is a command already sent: no noise reaches it, so the filter's gain on it
        // is zero and predict loads it with the input the plant received.
        constexpr auto model = *discretize(heater_block(), kTs, InputDelay{kFrac});
        constexpr auto res = design::lqg(model, kLqModel, design::KalmanNoise<1, 1>{.process = 0.02, .measurement = 0.01});
        static_assert(res.success);
        CHECK(damp::abs(res.kalman.M(1, 0)) <= 1e-12 * damp::abs(res.kalman.M(0, 0)));

        for (const bool record : {true, false}) {
            CAPTURE(record);
            LQG           law{res};
            ModelPlant<2> plant{delayed_heater(kFrac), ColVec<2>{-40.0, 0.0}}; // 40 °C under the setpoint
            double        sent = 0.0;                                          // the input the plant received, or the command, at k − 1
            int           clamped = 0;
            for (int k = 0; k < 30; ++k) {
                const double u = law.step(ColVec<1>{plant.y()})[0]; // predict with k − 1's input, update, control
                CHECK(law.estimator.state_estimate()(1) == doctest::Approx(sent).epsilon(1e-12).scale(10.0));
                const double applied = std::clamp(u, -10.0, 10.0);
                clamped += (applied != u) ? 1 : 0;
                if (record) {
                    law.back_calculate(ColVec<1>{applied});
                }
                plant.step(applied, 0.0);
                sent = record ? applied : u;
            }
            CHECK(clamped >= 3);
        }
    }

    TEST_CASE("offset-free LQG and MPC on the model hold the true delayed plant at r under a load; the delay-blind designs overshoot") {
        constexpr design::MPCWeights<1, 1>     w{.Rdu = 10.0}; // per step: Qy = 1
        constexpr design::MPCConstraints<1, 1> c{.u_lim = {
                                                     .lower = 0.0,
                                                     .upper = 40.0,
                                                 }};
        constexpr design::InputDisturbance     load{};
        const auto                             lqg_blind = design::oflqg(*kPlain, kLq, kNoise, load);
        const auto                             mpc_blind = design::ofmpc<30, 4>(*kPlain, w, c, kNoise, load);
        REQUIRE(lqg_blind.success);
        REQUIRE(mpc_blind.success);

        for (const double theta : {kFrac, kTs}) {
            CAPTURE(theta);
            const auto model = discretize(heater_block(), kTs, InputDelay{theta});
            REQUIRE(model.has_value());
            const auto lqg = design::oflqg(*model, kLqModel, kNoise, load);
            const auto mpc = design::ofmpc<30, 4>(*model, w, c, kNoise, load);
            REQUIRE(lqg.success);
            REQUIRE(mpc.success);
            using Est = std::remove_cvref_t<decltype(OFMPC{mpc}.estimator())>;
            static_assert(Est::NXE == 2 && std::is_same_v<std::remove_cvref_t<decltype(std::declval<const Est&>().state())>, ColVec<3>>); // [x; v(k−1)], then d

            const auto truth = delayed_heater(theta);
            OFLQG      lqg_law{lqg};
            OFLQG      lqg_blind_law{lqg_blind};
            OFMPC      mpc_law{mpc};
            OFMPC      mpc_blind_law{mpc_blind};
            const auto a = setpoint_run(lqg_law, truth, 100.0);
            const auto b = setpoint_run(lqg_blind_law, truth, 100.0);
            const auto m = setpoint_run(mpc_law, truth, 185.0);
            const auto n = setpoint_run(mpc_blind_law, truth, 185.0);
            MESSAGE("θ = ", theta, " s overshoot: LQG ", a.overshoot, " (blind ", b.overshoot, "), MPC ", m.overshoot, " (blind ", n.overshoot, ")");
            CHECK(std::abs(a.y_end - 100.0) < 1e-3); // offset-free under the 2 W loss
            CHECK(std::abs(m.y_end - 185.0) < 1e-3);
            CHECK(a.overshoot < b.overshoot);
            CHECK(m.overshoot < n.overshoot);
        }
    }

    TEST_CASE("offset-free LQG from a continuous plant: the Van Loan noise of the disturbance augment, the measurement variance as given") {
        constexpr double                               q = 0.01;  // process density on the block temperature
        constexpr double                               qd = 0.05; // the heat loss's random walk, a density
        constexpr design::KalmanNoise<1, 1, double, 1> noise{.process = q, .measurement = 0.01, .disturbance = qd};
        constexpr double                               a = kA;
        constexpr double                               b = kB;

        // the regulator is lqr(plant, Ts, w), and the noise the closed form of
        // ∫ e^{Aη} Qc e^{Aᵀη} dη on the input-disturbance augment [a b; 0 0]
        const auto res = design::oflqg(heater_block(), kTs, kLq, noise, design::InputDisturbance{});
        REQUIRE(res.success);
        CHECK(max_diff(res.K, design::lqr(heater_block(), kTs, kLq).K) == 0.0);
        const double e1 = std::exp(a * kTs) - 1.0;
        const double e2 = std::exp(2.0 * a * kTs) - 1.0;
        const double g = b / a;
        const double xx = (q * e2 / (2.0 * a)) + (qd * g * g * ((e2 / (2.0 * a)) - (2.0 * e1 / a) + kTs));
        const double xd = qd * g * ((e1 / a) - kTs);
        CHECK(res.kalman.Q(0, 0) == doctest::Approx(xx).epsilon(1e-12));
        CHECK(res.kalman.Q(0, 1) == doctest::Approx(xd).epsilon(1e-12));
        CHECK(res.kalman.Q(1, 0) == doctest::Approx(xd).epsilon(1e-12));
        CHECK(res.kalman.Q(1, 1) == doctest::Approx(qd * kTs).epsilon(1e-12));
        CHECK(res.kalman.R(0, 0) == 0.01); // a variance per reading: no Ts scaling

        // offset-free MPC from the continuous plant: the same estimator noise
        const auto ofmpc = design::ofmpc<30, 4>(heater_block(), kTs, {}, {}, noise, design::InputDisturbance{});
        REQUIRE(ofmpc.success);
        CHECK(max_diff(ofmpc.kalman.Q, res.kalman.Q) == 0.0);
        CHECK(ofmpc.kalman.R(0, 0) == 0.01);
    }
}

namespace {

/// Ticks of the recorded runs, and their commands from the AppliedInput step these replace.
constexpr int kAt[] = {5, 40, 90, 130, 199, 260, 299};

double setpoint_at(int k) {
    return ((k / 120) % 2 == 0) ? 0.6 : -0.4;
}

/// Compare a run's commands at kAt with the recorded ones, relative to the largest recorded.
struct Recorded {
    const double* golden;
    size_t        next = 0;
    double        worst = 0.0;
    void          at(int k, double u) {
        if (next < 7 && kAt[next] == k) {
            worst = std::max(worst, std::abs(u - golden[next]));
            ++next;
        }
    }
    [[nodiscard]] double relative() const {
        double scale = 0.0;
        for (size_t i = 0; i < 7; ++i) {
            scale = std::max(scale, std::abs(golden[i]));
        }
        return worst / scale;
    }
};

} // namespace

TEST_SUITE("back_calculate: the applied input through the law's own last command") {
    // Commands recorded from the AppliedInput step, step(…, {.u = u_sat, .u_ff = u_ff}), which
    // back_calculate(u_sat, u_ff) at the end of each tick replaces.
    const design::KalmanNoise<2, 1> kNoiseSpring{.process = kQn, .measurement = kRn};

    TEST_CASE("LQGI with a clamp; without back_calculate the step predicts with the last command") {
        constexpr double clamped[] = {33.210771066535429, 21.698602436461023, 0.30528446320147395, -24.341242954954595, -4.0558405401675799, 67.446982029701559, 40.217762773741939};
        constexpr double unclamped[] = {65.106839026183209, 6.5445551085095417, 2.4076864909603728, -15.855437674242751, -1.7288642913350714, -58.284596307229833, 1.4312379783334563};
        const auto       res = design::lqgi(kPlantD, {.Q = kQ, .R = kR, .Q_integral = kQint}, kNoiseSpring);
        for (const bool clamp : {true, false}) {
            LQGI      law{res};
            ColVec<2> x{};
            Recorded  rec{clamp ? clamped : unclamped};
            for (int k = 0; k < 300; ++k) {
                const double u = law.step(ColVec<1>{setpoint_at(k)}, ColVec<1>{(kPlantD.C * x)(0, 0)})[0];
                double       in = u;
                if (clamp) {
                    in = std::clamp(u, -0.8, 0.8);
                    law.back_calculate(ColVec<1>{in});
                }
                x = kPlantD.A * x + kPlantD.B * ColVec<1>{in};
                rec.at(k, u);
            }
            CHECK(rec.next == 7);
            CHECK(rec.relative() < 1e-12);
        }
    }

    TEST_CASE("a feedforward inside the limits is not saturation: the integrator runs as without one") {
        auto      with_ff = LQI{design::lqi(kPlantD, {.Q = kQ, .R = kR, .Q_integral = kQint})};
        auto      plain = with_ff;
        ColVec<2> x{};
        for (int k = 0; k < 200; ++k) {
            const double u = with_ff.step(ColVec<1>{setpoint_at(k)}, x)[0];
            static_cast<void>(plain.step(ColVec<1>{setpoint_at(k)}, x));
            const double ff = 0.3 * std::sin(0.05 * k);
            with_ff.back_calculate(ColVec<1>{u + ff}, ColVec<1>{ff}); // no clamp: u_sat = u_cmd + u_ff
            x = kPlantD.A * x + kPlantD.B * ColVec<1>{u + ff};
        }
        CHECK(with_ff.integral.xi[0] == plain.integral.xi[0]);

        // The same input reported without its feedforward reads as a clamp and holds an integrator.
        auto      misread = LQI{design::lqi(kPlantD, {.Q = kQ, .R = kR, .Q_integral = kQint})};
        auto      reference = misread;
        ColVec<2> z{};
        bool      held = false;
        for (int k = 0; k < 200; ++k) {
            const double u = misread.step(ColVec<1>{setpoint_at(k)}, z)[0];
            static_cast<void>(reference.step(ColVec<1>{setpoint_at(k)}, z));
            const double ff = 0.3 * std::sin(0.05 * k);
            misread.back_calculate(ColVec<1>{u + ff});
            z = kPlantD.A * z + kPlantD.B * ColVec<1>{u + ff};
            held = held || (misread.integral.xi[0] != reference.integral.xi[0]);
        }
        CHECK(held);
    }
}
