// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/estimation/parameter_estimation.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

// Deterministic PRBS-like excitation (reproducible tests). The hold length
// matters: identification needs energy at the plant's own timescale, so each
// level is held for `hold` samples — the same reason the PRBS generator has a
// clock divider. A fast alternating input is filtered away by a slow plant and
// leaves the pole nearly unobservable.
constexpr damp::array<double, 12> pattern{
    1.0,
    -1.0,
    1.0,
    1.0,
    -1.0,
    -1.0,
    1.0,
    -1.0,
    -1.0,
    1.0,
    1.0,
    -1.0,
};
constexpr double excitation(size_t k, size_t hold) {
    return pattern[(k / hold) % pattern.size()];
}

} // namespace

TEST_CASE("first-order gain and time constant are recovered from noiseless data") {
    // Plant: tau*y' + y = K*u with K = 2, tau = 0.5 s, sampled at 10 ms.
    constexpr double Ts = 0.01;
    constexpr double K = 2.0;
    constexpr double tau = 0.5;
    const double     a = damp::exp(-Ts / tau);
    const double     b = K * (1.0 - a);

    design::FirstOrderPlantEstimatorConfig<double> cfg{};
    cfg.Ts = Ts;
    cfg.initial_covariance = 1e4; // weak prior: don't bias the near-unit pole
    FirstOrderPlantEstimator estimator{cfg};

    double y = 0.0;
    double u_prev = 0.0;
    double confidence_early = 0.0;
    for (size_t k = 0; k < 1200; ++k) {
        y = (a * y) + (b * u_prev);
        const double u = excitation(k, 25); // hold ~tau/2: energy where the pole lives
        (void)estimator.update(u, y);
        u_prev = u;
        if (k == 10) {
            confidence_early = estimator.confidence();
        }
    }

    REQUIRE(estimator.valid());
    CHECK(estimator.gain() == doctest::Approx(K).epsilon(0.01));
    CHECK(estimator.time_constant() == doctest::Approx(tau).epsilon(0.01));
    CHECK(estimator.parameters()(0) == doctest::Approx(estimator.gain()));
    CHECK(estimator.parameters()(1) == doctest::Approx(estimator.time_constant()));
    CHECK(estimator.confidence() > confidence_early); // excitation shrinks the covariance
}

TEST_CASE("recovery holds within 10% under measurement noise") {
    constexpr double Ts = 0.01;
    constexpr double K = 2.0;
    constexpr double tau = 0.5;
    const double     a = damp::exp(-Ts / tau);
    const double     b = K * (1.0 - a);

    design::FirstOrderPlantEstimatorConfig<double> cfg{};
    cfg.Ts = Ts;
    cfg.initial_covariance = 1e4;
    FirstOrderPlantEstimator estimator{cfg};

    double y = 0.0;
    double u_prev = 0.0;
    for (size_t k = 0; k < 3000; ++k) {
        y = (a * y) + (b * u_prev);
        const double u = excitation(k, 25);
        const double y_meas = y + (0.02 * pattern[(k + 5) % pattern.size()]); // deterministic "sensor noise"
        (void)estimator.update(u, y_meas);
        u_prev = u;
    }

    REQUIRE(estimator.valid());
    CHECK(estimator.gain() == doctest::Approx(K).epsilon(0.10));
    CHECK(estimator.time_constant() == doctest::Approx(tau).epsilon(0.10));
}

TEST_CASE("without excitation the estimator honestly refuses to vouch") {
    design::FirstOrderPlantEstimatorConfig<double> cfg{};
    cfg.Ts = 0.01;
    FirstOrderPlantEstimator estimator{cfg};

    CHECK(!estimator.valid()); // no data at all
    CHECK(estimator.confidence() == doctest::Approx(0.0));
    // Physical readouts are gated: zeros when !valid (no /0 or log of garbage)
    CHECK(estimator.gain() == doctest::Approx(0.0));
    CHECK(estimator.time_constant() == doctest::Approx(0.0));

    for (size_t k = 0; k < 200; ++k) {
        (void)estimator.update(0.0, 0.0); // flat-lined plant: zero information
    }
    CHECK(!estimator.valid());
    CHECK(estimator.confidence() < 0.5);
    CHECK(estimator.gain() == doctest::Approx(0.0));
}

TEST_CASE("first-order FO plant maps to mechanical J, b on a servo axis") {
    // Same FirstOrderPlantEstimator: J*omega' + b*omega = tau
    //   K = 1/b,  tau_mech = J/b  →  b = 1/K,  J = tau_mech/K  (example recipe).
    constexpr double Ts = 0.01;
    constexpr double J = 0.01;
    constexpr double b_damp = 0.05;
    const double     a = damp::exp(-b_damp * Ts / J); // tau_mech = J/b = 0.2 s
    const double     b_d = (1.0 / b_damp) * (1.0 - a);

    design::FirstOrderPlantEstimatorConfig<double> cfg{};
    cfg.Ts = Ts;
    cfg.initial_covariance = 1e4;
    FirstOrderPlantEstimator estimator{cfg};

    double omega = 0.0;
    double torque_prev = 0.0;
    for (size_t k = 0; k < 1200; ++k) {
        omega = (a * omega) + (b_d * torque_prev);
        const double torque = 0.5 * excitation(k, 10);
        (void)estimator.update(torque, omega);
        torque_prev = torque;
    }

    REQUIRE(estimator.valid());
    REQUIRE(estimator.gain() > 0.0);
    const double b_est = 1.0 / estimator.gain();
    const double J_est = estimator.time_constant() / estimator.gain();
    CHECK(J_est == doctest::Approx(J).epsilon(0.01));
    CHECK(b_est == doctest::Approx(b_damp).epsilon(0.01));
}

TEST_CASE("drift monitor fires on a plant change and stays quiet otherwise") {
    ParameterDriftMonitor<double> monitor{0.1, 0.01, 3.0};

    // Baseline: small residuals — no drift while the model explains the data.
    bool fired_during_baseline = false;
    for (size_t k = 0; k < 400; ++k) {
        fired_during_baseline = monitor.update(0.01 * pattern[k % pattern.size()]) || fired_during_baseline;
    }
    CHECK(!fired_during_baseline);
    CHECK(monitor.drift_score() == doctest::Approx(1.0).epsilon(0.2));

    // The plant moves: residuals jump by an order of magnitude.
    bool fired_after_change = false;
    for (size_t k = 0; k < 100; ++k) {
        fired_after_change = monitor.update(0.15 * pattern[k % pattern.size()]) || fired_after_change;
    }
    CHECK(fired_after_change);

    monitor.reset();
    CHECK(!monitor.drift_detected());
}

TEST_CASE("adaptation gate blocks invalid, unconfident, and wild updates") {
    AdaptationGate<2, double> gate{0.5, ColVec<2>{0.5, 0.1}};

    const ColVec<2> nominal{2.0, 0.5};

    CHECK(!gate.evaluate(false, 0.9, nominal).allow_update); // estimator does not vouch
    CHECK(!gate.evaluate(true, 0.3, nominal).allow_update);  // confidence below the floor
    CHECK(!gate.has_accepted());

    const auto first = gate.evaluate(true, 0.9, nominal); // first acceptance sets the baseline
    CHECK(first.allow_update);
    CHECK(gate.accepted_parameters()(0) == doctest::Approx(2.0));

    CHECK(!gate.evaluate(true, 0.9, ColVec<2>{3.0, 0.5}).allow_update); // parameter 0 jumps by 1.0 > 0.5
    CHECK(gate.accepted_parameters()(0) == doctest::Approx(2.0));       // baseline unchanged on a block

    const auto small_step = gate.evaluate(true, 0.9, ColVec<2>{2.3, 0.45});
    CHECK(small_step.allow_update);
    CHECK(gate.accepted_parameters()(0) == doctest::Approx(2.3));
}

TEST_CASE("grey-box estimation is constexpr") {
    constexpr double recovered_gain = [] {
        constexpr double Ts = 0.01;
        constexpr double a = 0.98;       // tau ~ 0.495 s
        constexpr double b = 2.0 * 0.02; // K = 2

        design::FirstOrderPlantEstimatorConfig<double> cfg{};
        cfg.Ts = Ts;
        cfg.initial_covariance = 1e4;
        FirstOrderPlantEstimator<double> estimator{cfg};

        double y = 0.0;
        double u_prev = 0.0;
        for (size_t k = 0; k < 600; ++k) {
            y = (a * y) + (b * u_prev);
            const double u = excitation(k, 25);
            (void)estimator.update(u, y);
            u_prev = u;
        }
        return estimator.gain();
    }();
    static_assert(recovered_gain > 1.9 && recovered_gain < 2.1);
}
