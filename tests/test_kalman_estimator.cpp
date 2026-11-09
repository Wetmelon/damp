// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file test_kalman_estimator.cpp
 * @brief The shared core of the laws: the options by type, the disturbance specs, design::kalman
 *        with a disturbance spec, the FixedGainObserver running its gain, and SteadyStateTarget.
 *
 * The disturbance-augmented filter and the target are checked against the offset-free design they
 * compose into (design::oflqg, OFLQG): the same gains, maps and commands. A delayed plant is a
 * model, discretize(plant, Ts, InputDelay{θ}), designed on like any other: its register holds a
 * command already sent, so the filter puts no gain on it and
 * its plant part is the undelayed filter, and the target is the undelayed steady state with the
 * register holding the input that reaches the plant. The law's filter, run in deviations from the
 * target it is anchored on, is checked against the textbook offset-free law on the absolute state,
 * u = u_t − K(x̂ − x_t): the same estimates and commands.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/lqg.hpp"
#include "damp/controllers/oflqg.hpp"
#include "damp/controllers/state_feedback_core.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/options.hpp"
#include "damp/design/status.hpp"
#include "damp/design/steady_state_target.hpp"
#include "damp/estimation/fixed_gain_observer.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "doctest.h"
#include "register_form.hpp"

using namespace damp;

namespace {

constexpr double kTs = 0.01;

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

/// The mass-spring with a position and a speed sensor.
constexpr StateSpace<2, 1, 2, double, 2, 2> two_sensors() {
    StateSpace<2, 1, 2, double, 2, 2> s{};
    s.A = mass_spring().A;
    s.B = mass_spring().B;
    s.C(0, 0) = 1.0;
    s.C(1, 1) = 1.0;
    return s;
}

/// The discrete mass-spring with its process noise on the states (G = I), as the offset-free designs model it.
StateSpace<2, 1, 1, double, 2, 1> plant_d() {
    auto p = discretize(mass_spring(), kTs).value();
    p.G = Matrix<2, 2>::identity();
    return p;
}

const auto         kPlantD = plant_d();
const Matrix<2, 2> kQp{{1e-6, 0.0}, {0.0, 1e-4}};
const Matrix<1, 1> kRn{{1e-6}};
const Matrix<1, 1> kQd{{1e-3}};
const Matrix<2, 2> kQ{{1e3, 0.0}, {0.0, 1.0}};
const Matrix<3, 3> kQ3{{1e3, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 0.0}}; // kQ on [x; v(k−1)], the register unweighted
const Matrix<1, 1> kR{{1e-2}};

template<typename A, typename B>
double gap(const A& a, const B& b) {
    return (a - b).norm() / std::max(1.0, b.norm());
}

/// Gap relative to the reference, whatever its size.
template<typename A, typename B>
double rel(const A& a, const B& b) {
    return (a - b).norm() / b.norm();
}

/// Leading N entries of a vector.
template<size_t N, size_t M>
ColVec<N> lead(const ColVec<M>& v) {
    ColVec<N> out{};
    for (size_t i = 0; i < N; ++i) {
        out(i) = v(i);
    }
    return out;
}

/// The rows of @p m on the plant of a delayed model [x; v(k−1)]: m, then NU zero rows for the register.
template<size_t NU, size_t NX, size_t NC>
Matrix<NX + NU, NC> padded(const Matrix<NX, NC>& m) {
    Matrix<NX + NU, NC> out{};
    out.template block<NX, NC>(0, 0) = m;
    return out;
}

/**
 * The filter on a delayed model [x; v(k−1)] against the plain filter on the undelayed plant. The
 * register is loaded with the command and gets no process noise: it is known exactly, so its
 * covariance and gain are zero and the plant's part solves the undelayed Riccati equation.
 */
template<size_t NX, size_t NXD, size_t NU, size_t NY, size_t NW, size_t NV>
void check_register_known(
    const StateSpace<NXD, NU, NY, double, NW, NV>& delayed,
    const StateSpace<NX, NU, NY, double, NW, NV>&  undelayed,
    const design::KalmanNoise<NW, NV>&             n
) {
    static_assert(NXD == NX + NU, "The delayed model is the plant, then one register per input.");
    const auto filt = design::kalman(delayed, n);
    const auto ref = design::kalman(undelayed, n);
    static_assert(std::is_same_v<std::remove_cvref_t<decltype(filt)>, design::KalmanResult<NXD, NU, NY, double, NW, NV>>);
    REQUIRE(filt.success);
    REQUIRE(ref.success);
    const Matrix<NX, NY>  L_x = filt.M.template block<NX, NY>(0, 0).to_matrix();
    const Matrix<NU, NY>  L_v = filt.M.template block<NU, NY>(NX, 0).to_matrix();
    const Matrix<NX, NX>  P_x = filt.P.template block<NX, NX>(0, 0).to_matrix();
    const Matrix<NU, NXD> P_v = filt.P.template block<NU, NXD>(NX, 0).to_matrix();
    MESSAGE("L_x ", rel(L_x, ref.M), ", L_v ", L_v.norm() / ref.M.norm(), ", P_x ", rel(P_x, ref.P), ", P_v ", P_v.norm() / ref.P.norm());
    CHECK(rel(L_x, ref.M) < 1e-12);
    CHECK(L_v.norm() < 1e-12 * ref.M.norm());
    CHECK(rel(P_x, ref.P) < 1e-12);
    CHECK(P_v.norm() < 1e-12 * ref.P.norm());
}

} // namespace

TEST_SUITE("Design options by type") {
    using damp::detail::option_count;
    using damp::detail::option_t;
    using design::detail::is_disturbance;
    using design::detail::is_kalman_noise;
    using Noise = design::KalmanNoise<2, 1, double, 1>;

    TEST_CASE("each kind found in any order, its default when absent, its sizes in its type") {
        static_assert(option_count<is_disturbance, Noise, design::OutputDisturbance<2>> == 1);
        static_assert(option_count<is_kalman_noise, design::OutputDisturbance<2>, Noise> == 1);
        static_assert(option_count<is_disturbance, Noise> == 0);
        static_assert(std::is_same_v<option_t<is_disturbance, design::NoDisturbance, Noise, design::OutputDisturbance<2>>, design::OutputDisturbance<2>>);
        static_assert(std::is_same_v<option_t<is_disturbance, design::NoDisturbance, Noise>, design::NoDisturbance>);

        constexpr auto spec = damp::detail::option<is_disturbance>(design::NoDisturbance{}, Noise{}, design::OutputDisturbance<2>{{0, 3}});
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(spec)>, design::OutputDisturbance<2>>);
        static_assert(spec.outputs[1] == 3);
        constexpr auto none = damp::detail::option<is_disturbance>(design::NoDisturbance{}, Noise{});
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(none)>, design::NoDisturbance>);
    }

    TEST_CASE("an unknown type or a kind given twice is refused; a delay is the plant's model, not an option") {
        using TargetSet = design::detail::target_options;
        static_assert(TargetSet::known<> && TargetSet::unique<>);
        static_assert(TargetSet::known<design::OutputDisturbance<>> && TargetSet::known<design::DisturbanceModel<2, 1, 1>>);
        static_assert(!TargetSet::known<design::OutputDisturbance<>, double>);
        static_assert(!TargetSet::unique<design::InputDisturbance<>, design::NoDisturbance>);                                       // two disturbance specs
        static_assert(TargetSet::known<Noise, design::InputDisturbance<>> && TargetSet::unique<Noise, design::InputDisturbance<>>); // the target takes the noise and ignores it
        static_assert(!TargetSet::unique<Noise, Noise>);
        //! discretize(plant, Ts, InputDelay{θ}) is the delayed plant: the target does not take the delay.
        //! design::kalman takes the disturbance as an argument, not an option.
        static_assert(!TargetSet::known<InputDelay<>>);
    }
}

TEST_SUITE("Disturbance specs") {
    TEST_CASE("sizes: one per input, per output, per column of G, or the model's") {
        using design::detail::disturbance_traits;
        static_assert(disturbance_traits<design::NoDisturbance, 2, 3, 4>::states == 0);
        static_assert(disturbance_traits<design::InputDisturbance<>, 2, 3, 4>::states == 2);
        static_assert(disturbance_traits<design::OutputDisturbance<2>, 2, 3, 4>::states == 3);
        static_assert(disturbance_traits<design::OutputDisturbance<2>, 2, 3, 4>::listed == 2);
        static_assert(disturbance_traits<design::PlantDisturbance<>, 2, 3, 4>::states == 4);
        static_assert(disturbance_traits<design::DisturbanceModel<5, 3, 1>, 2, 3, 4>::states == 1);
        static_assert(design::detail::disturbance_states<2, 3, 4, design::KalmanNoise<4, 3, double, 3>, design::OutputDisturbance<>> == 3);
        static_assert(design::detail::disturbance_states<2, 3, 4, design::KalmanNoise<4, 3>> == 0);
    }

    TEST_CASE("the controlled outputs: the listed ones, checked, or all") {
        constexpr auto all = design::detail::selected_outputs<3>(damp::array<size_t, 0>{});
        static_assert(all && (*all)[0] == 0 && (*all)[2] == 2);
        constexpr auto two = design::detail::selected_outputs<3>(damp::array<size_t, 2>{{2, 0}});
        static_assert(two && (*two)[0] == 2 && (*two)[1] == 0);
        static_assert(!design::detail::selected_outputs<3>(damp::array<size_t, 1>{{3}}));    // out of range
        static_assert(!design::detail::selected_outputs<3>(damp::array<size_t, 2>{{1, 1}})); // twice
    }

    TEST_CASE("KalmanNoise carries the disturbance's random walk only with a disturbance model") {
        constexpr design::KalmanNoise plain{Matrix<2, 2>::identity(), Matrix<1, 1>{{0.1}}};
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(plain)>, design::KalmanNoise<2, 1, double, 0>>);
        constexpr design::KalmanNoise with_d{Matrix<2, 2>::identity(), Matrix<1, 1>{{0.1}}, Matrix<1, 1>{{0.5}}};
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(with_d)>, design::KalmanNoise<2, 1, double, 1>>);
        static_assert(with_d.disturbance(0, 0) == 0.5);
    }
}

TEST_SUITE("design::kalman with a disturbance spec") {
    TEST_CASE("no options is the plain filter; on a delayed plant the register is known: no gain on it, the plant's part the undelayed filter") {
        const design::KalmanNoise<2, 1> n{.process = kQp, .measurement = kRn};
        const auto                      plain = design::kalman(kPlantD, n);
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(plain)>, design::KalmanResult<2, 1, 1, double, 2, 1>>);
        REQUIRE(plain.success);
        CHECK((design::kalman(mass_spring(), kTs, n).M - design::lqg(mass_spring(), kTs, {.Q = kQ, .R = kR}, n).kalman.M).norm() == 0.0);

        //! a whole sample on the discrete plant: x⁺ = A x + B v(k−1)
        check_register_known<2>(discretize(kPlantD, kTs, InputDelay{kTs}).value(), kPlantD, n);

        //! a fraction of a period from the continuous plant: x⁺ = Φx + Γ₀v(k) + Γ₁v(k−1), the process noise on the plant's states
        auto part = discretize(mass_spring(), kTs, InputDelay{0.4 * kTs}).value();
        part.G = padded<1>(Matrix<2, 2>::identity());
        check_register_known<2>(part, kPlantD, n);
    }

    TEST_CASE("an input disturbance is the offset-free filter and target of oflqg") {
        const auto old = design::oflqg(kPlantD, {.Q = kQ, .R = kR}, design::KalmanNoise{kQp, kRn, kQd}, design::InputDisturbance{});
        REQUIRE(old.success);
        const auto filt = design::kalman(kPlantD, design::KalmanNoise{.process = kQp, .measurement = kRn, .disturbance = kQd}, design::InputDisturbance{});
        const auto tg = design::steady_state_target(kPlantD, design::InputDisturbance{});
        REQUIRE(filt.success);
        REQUIRE(tg.success);
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(filt)>, design::KalmanResult<3, 1, 1, double, 3, 1, 1>>);
        CHECK(gap(filt.M, old.kalman.M) < 1e-12);
        CHECK(gap(filt.sys.A, old.kalman.sys.A) == 0.0);
        CHECK(gap(tg.target_d, old.target.target_d) < 1e-12);
        CHECK(gap(tg.target_r, old.target.target_r) < 1e-12);
    }

    TEST_CASE("on a delayed plant: the composed design's filter and target, the undelayed steady state with the register at the input") {
        // The mass-spring holds x = r with u = 4r − d under a load d at its input (ẍ = −4x − 0.4ẋ + u + d).
        // ZOH keeps that steady state and a delay moves none: the register holds what reaches the plant.
        const design::KalmanNoise<2, 1, double, 1> nn{kQp, kRn, kQd};

        //! a whole sample; the load on the plant, undelayed: B_d = [B; 0], so the register stays exact
        const auto                              whole = discretize(kPlantD, kTs, InputDelay{kTs}).value();
        const design::DisturbanceModel<3, 1, 1> load{.Bd = padded<1>(kPlantD.B)};
        const auto                              old = design::oflqg(whole, {.Q = kQ3, .R = kR}, nn, load);
        REQUIRE(old.success);
        const auto filt = design::kalman(whole, nn, load);
        const auto tg = design::steady_state_target(whole, load);
        REQUIRE(filt.success);
        REQUIRE(tg.success);
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(filt)>, design::KalmanResult<4, 1, 1, double, 3, 1, 1>>); // [x; v(k−1); d]
        CHECK(gap(filt.M, old.kalman.M) < 1e-12);
        CHECK(gap(tg.target_d, old.target.target_d) < 1e-12);
        CHECK(gap(tg.target_r, old.target.target_r) < 1e-12);
        CHECK(std::abs(filt.M(2, 0)) < 1e-12 * filt.M.norm());                       // a command already sent: no gain on the register
        CHECK(gap(tg.target_r, Matrix<4, 1>{{1.0}, {0.0}, {4.0}, {4.0}}) < 1e-12);   // [x_t; v_t; u_t] per unit r
        CHECK(gap(tg.target_d, Matrix<4, 1>{{0.0}, {0.0}, {-1.0}, {-1.0}}) < 1e-12); // u_t = −d, and the register holds it

        //! a fraction of a period from the continuous plant; d adds to the command, so it enters the register too: B_d = B = [Γ₀; I]
        auto part = discretize(mass_spring(), kTs, InputDelay{0.4 * kTs}).value();
        part.G = padded<1>(Matrix<2, 2>::identity());
        const auto old_c = design::oflqg(part, {.Q = kQ3, .R = kR}, nn, design::InputDisturbance{});
        REQUIRE(old_c.success);
        const auto filt_c = design::kalman(part, nn, design::InputDisturbance{});
        const auto tg_c = design::steady_state_target(part, design::InputDisturbance{});
        REQUIRE(filt_c.success);
        REQUIRE(tg_c.success);
        CHECK(gap(filt_c.M, old_c.kalman.M) < 1e-12);
        CHECK(gap(tg_c.target_d, old_c.target.target_d) < 1e-12);
        CHECK(gap(tg_c.target_r, old_c.target.target_r) < 1e-12);
        CHECK(gap(tg_c.target_r, Matrix<4, 1>{{1.0}, {0.0}, {4.0}, {4.0}}) < 1e-12);
        CHECK(gap(tg_c.target_d, Matrix<4, 1>{{0.0}, {0.0}, {0.0}, {-1.0}}) < 1e-12); // the register holds u_t + d = 4r
    }

    TEST_CASE("a non-finite DisturbanceModel::Bd is refused, on a continuous plant and on a discrete one") {
        //! It returned success with a NaN target_d: the target has no filter to refuse it.
        const design::DisturbanceModel<2, 1, 1> bad{.Bd = Matrix<2, 1>{{0.0}, {std::numeric_limits<double>::quiet_NaN()}}};
        const auto                              continuous = design::steady_state_target(mass_spring(), kTs, bad);
        CHECK_FALSE(continuous.success);
        CHECK(continuous.status == design::DesignStatus::bad_weights);
        CHECK(design::steady_state_target(plant_d(), bad).status == design::DesignStatus::bad_weights);
        CHECK(design::steady_state_target(mass_spring(), kTs, design::DisturbanceModel<2, 1, 1>{.Bd = mass_spring().B}).success);
    }

    TEST_CASE("a discrete plant with a non-finite entry is refused, as a continuous one is") {
        //! The discrete plant skipped discretize's check: a NaN A gave success with a finite target
        //! (NZ > NU), a NaN B came back singular.
        const double nan = std::numeric_limits<double>::quiet_NaN();
        auto         nan_a = plant_d();
        nan_a.A(0, 0) = nan;
        CHECK(design::steady_state_target(nan_a, design::InputDisturbance{}).status == design::DesignStatus::bad_weights);
        auto nan_b = plant_d();
        nan_b.B(1, 0) = nan;
        CHECK(design::steady_state_target(nan_b).status == design::DesignStatus::bad_weights);
        auto nan_c = discretize(two_sensors(), kTs).value();
        nan_c.C(1, 1) = nan;
        CHECK(design::steady_state_target(nan_c).status == design::DesignStatus::bad_weights);
        CHECK(design::steady_state_target(plant_d()).success);
    }

    TEST_CASE("more outputs than inputs: the target is the equilibrium closest to r") {
        //! The least squares on the stacked [A − I, B; C, D] traded the equilibrium rows against
        //! the outputs: x_t = [1.0089, 1.9911], u_t = 1.9429 for r = [1, 2], O(Ts) off an
        //! equilibrium, and the loop settled where the weights and Ts put it.
        const StateSpace<2, 1, 2> chain{.A = {{-1.0, 1.0}, {0.0, -1.0}}, .B = {{0.0}, {1.0}}, .C = {{1.0, 0.0}, {0.0, 1.0}}};
        const auto                tg = design::steady_state_target(chain, 0.1);
        REQUIRE(tg.success);
        const ColVec<2> r{1.0, 2.0};
        const auto      s = tg.target_r * r;
        CHECK(s(0) == doctest::Approx(1.5).epsilon(1e-9)); // x₁ = x₂ = u at every steady state: the mean of r
        CHECK(s(1) == doctest::Approx(1.5).epsilon(1e-9));
        CHECK(s(2) == doctest::Approx(1.5).epsilon(1e-9));
    }

    TEST_CASE("an explicit model and the listed outputs") {
        auto sys = discretize(two_sensors(), kTs).value();
        sys.G = Matrix<2, 2>::identity(); // the process noise on the states
        const auto old = design::oflqg(sys, {.Q = kQ, .R = kR}, design::KalmanNoise{kQp, Matrix<2, 2>{{1e-6, 0.0}, {0.0, 1e-4}}, kQd}, design::DisturbanceModel<2, 2, 1, 1>{.Bd = sys.B, .outputs = {{0}}});
        REQUIRE(old.success);
        const design::DisturbanceModel<2, 2, 1, 1> spec{.Bd = sys.B, .outputs = {{0}}};
        const auto                                 filt = design::kalman(sys, design::KalmanNoise{.process = kQp, .measurement = Matrix<2, 2>{{1e-6, 0.0}, {0.0, 1e-4}}, .disturbance = kQd}, spec);
        const auto                                 tg = design::steady_state_target(sys, spec);
        REQUIRE(filt.success);
        REQUIRE(tg.success);
        CHECK(gap(filt.M, old.kalman.M) < 1e-12);
        CHECK(gap(tg.target_d, old.target.target_d) < 1e-12);
        CHECK(gap(tg.target_r, old.target.target_r) < 1e-12);

        //! a force through the plant's G that the sensor also reads through H: B_d = G, C_d = H
        StateSpace<2, 1, 1, double, 1, 1> g{.A = mass_spring().A, .B = mass_spring().B, .C = mass_spring().C};
        g.G = mass_spring().B;
        g.H = Matrix<1, 1>{{0.5}};
        const auto gd = discretize(g, kTs).value();
        const auto filt_g = design::kalman(gd, design::KalmanNoise{.process = Matrix<1, 1>{{1e-4}}, .measurement = Matrix<1, 1>{{1e-6}}, .disturbance = kQd}, design::PlantDisturbance{});
        REQUIRE(filt_g.success);
        CHECK(gap(filt_g.sys.A.block<2, 1>(0, 2).to_matrix(), gd.G) == 0.0);
        CHECK(gap(filt_g.sys.C.block<1, 1>(0, 2).to_matrix(), gd.H) == 0.0);
    }

    TEST_CASE("refusals: an undetectable augment, an output out of range, a continuous plant") {
        //! an output bias on an integrator: rank [A − I, 0; C, I] < NX + ND
        const StateSpace<1, 1, 1, double, 1, 1> integrator{.A = Matrix<1, 1>{{1.0}}, .B = Matrix<1, 1>{{kTs}}, .C = Matrix<1, 1>{{1.0}}, .Ts = kTs};
        const auto                              undetectable = design::kalman(integrator, design::KalmanNoise{.process = Matrix<1, 1>{{1e-4}}, .measurement = Matrix<1, 1>{{1e-6}}, .disturbance = kQd}, design::OutputDisturbance{});
        CHECK(undetectable.status == design::DesignStatus::not_detectable);
        CHECK(design::kalman(integrator, design::KalmanNoise{.process = Matrix<1, 1>{{1e-4}}, .measurement = Matrix<1, 1>{{1e-6}}, .disturbance = kQd}, design::InputDisturbance{}).success);

        CHECK(design::steady_state_target(kPlantD, design::InputDisturbance<1>{{1}}).status == design::DesignStatus::bad_weights);
        CHECK(design::kalman(mass_spring(), design::KalmanNoise{.process = kQp, .measurement = kRn, .disturbance = kQd}, design::InputDisturbance{}).status == design::DesignStatus::invalid_sample_time);
    }
}

TEST_SUITE("FixedGainObserver and SteadyStateTarget") {
    TEST_CASE("the disturbance spec sets the noise's disturbance count, so the noise is a braced list") {
        const auto braced = design::kalman(kPlantD, {.process = kQp, .measurement = kRn, .disturbance = kQd}, design::InputDisturbance{});
        const auto typed = design::kalman(kPlantD, design::KalmanNoise{.process = kQp, .measurement = kRn, .disturbance = kQd}, design::InputDisturbance{});
        REQUIRE(braced.success);
        CHECK((braced.M - typed.M).norm() == 0.0);
    }

    TEST_CASE("no options: the plain filter, with the input it holds") {
        const auto        res = design::kalman(kPlantD, {.process = kQp, .measurement = kRn});
        FixedGainObserver est{res};
        static_assert(sizeof(KalmanFilter<2, 1, 1, double, 2, 1>) > sizeof(est)); // no P, Q, R, G or H
        ColVec<2> x{0.1, 0.0};
        ColVec<2> ref{};
        ColVec<1> u_held{};
        for (int k = 0; k < 100; ++k) {
            const ColVec<1> u{{0.2 * std::sin(0.1 * k)}};
            const ColVec<1> y{(kPlantD.C * x)(0, 0)};
            est.update(y);
            ref = ref + res.M * (y - (kPlantD.C * ref + kPlantD.D * u_held)); // the textbook update, the feedthrough on the held input
            est.predict(u);
            ref = kPlantD.A * ref + kPlantD.B * u;
            u_held = u;
            x = kPlantD.A * x + kPlantD.B * u;
        }
        CHECK((est.state_estimate() - ref).norm() <= 1e-14); // not 0: -ffast-math may order the two sides differently
        CHECK((est.at_command() - ref).norm() <= 1e-14);
    }

    TEST_CASE("the deviation observer matches the absolute observer, and anchor shifts the deviation exactly") {
        const auto res = design::kalman(
            kPlantD, design::KalmanNoise{.process = kQp, .measurement = kRn, .disturbance = kQd}, design::InputDisturbance{}
        );
        REQUIRE(res.success);
        const SteadyStateTarget target{design::steady_state_target(kPlantD, design::InputDisturbance{})};

        FixedGainObserver<3, 1, 1, double, 1> dev{res};
        FixedGainObserver<3, 1, 1, double>    absf{res}; // the design on [x; d] as an absolute observer
        static_assert(std::is_same_v<decltype(FixedGainObserver{res}), decltype(dev)>);
        ColVec<2> x{0.05, 0.0};
        ColVec<1> u{{0.1}};
        for (int k = 0; k < 40; ++k) {
            const ColVec<1> y{{(kPlantD.C * x)(0, 0)}};
            dev.predict(u);
            absf.predict(u);
            if (k == 16) {
                dev.update(y, u); // an explicit feedthrough input
                absf.update(y, u);
            } else {
                dev.update(y);
                absf.update(y, u);
            }
            ColVec<1> d_abs{};
            d_abs(0) = absf.state()(2);
            CHECK((dev.state_estimate() - lead<2>(absf.state())).norm() < 1e-9);
            CHECK((dev.disturbance_estimate() - d_abs).norm() < 1e-9);
            CHECK((dev.gain() - absf.gain()).norm() == 0.0);
            if (k == 15) {
                const auto xhat = dev.state_estimate();
                const auto dhat = dev.disturbance_estimate();
                const auto p = target.equilibrium(ColVec<1>{{0.2}}, dhat);
                dev.anchor(p.x, p.u);
                CHECK((dev.state_estimate() - xhat).norm() < 1e-12);
                CHECK((dev.disturbance_estimate() - dhat).norm() < 1e-12);
                CHECK((dev.anchor_state() - p.x).norm() == 0.0); // the anchor moved; the estimate did not
            }
            x = kPlantD.A * x + kPlantD.B * u;
            u = ColVec<1>{{0.1 * std::sin(0.2 * k)}};
        }
    }

    TEST_CASE("the estimator, the target and the gain run the offset-free law's commands, on a delayed plant") {
        const design::KalmanNoise<2, 1, double, 1> nn{kQp, kRn, kQd};
        const auto                                 whole = discretize(kPlantD, kTs, InputDelay{kTs}).value(); // [x; v(k−1)]
        const design::DisturbanceModel<3, 1, 1>    load{.Bd = padded<1>(kPlantD.B)};                          // the load on the plant, undelayed
        const auto                                 old = design::oflqg(whole, {.Q = kQ3, .R = kR}, nn, load);
        REQUIRE(old.success);
        OFLQG             law{old};
        FixedGainObserver est{design::kalman(whole, nn, load)};
        SteadyStateTarget target{design::steady_state_target(whole, load)};
        static_assert(decltype(est)::NXE == 3);
        const Matrix<1, 3> K = old.K; // on the law's model [x; v(k−1)]

        //! the plant with the true delay: the command lands one period late
        ColVec<2>        x{};
        double           sent = 0.0;
        ColVec<1>        u_new{};
        double           worst = 0.0;
        double           scale = 0.0;
        constexpr int    n = 1500;
        constexpr double kForce = 0.5;
        for (int k = 0; k < n; ++k) {
            const ColVec<1> r{{((k / 300) % 2 == 0) ? 0.3 : -0.2}};
            const ColVec<1> y{{(kPlantD.C * x)(0, 0)}};
            const double    u_old = law.step(r, y)[0];

            est.predict(u_new);
            est.update(y);
            target.anchor(est, r);
            u_new = est.anchor_input() - (K * lead<3>(est.at_command()));

            worst = std::max(worst, std::abs(u_old - u_new[0]));
            scale = std::max(scale, std::abs(u_old));
            const double force = (k > n / 3) ? kForce : 0.0;
            x = kPlantD.A * x + kPlantD.B * ColVec<1>{{sent + force}};
            sent = u_old;
        }
        MESSAGE("|Δu| / max |u| = ", worst / scale, ", |d̂ − load| = ", std::abs(est.disturbance_estimate()[0] - kForce));
        CHECK(worst <= 1e-12 * scale);
        CHECK(std::abs(x[0] - 0.3) < 1e-3); // offset-free under the load
        CHECK((est.disturbance_estimate() - law.disturbance_estimate()).norm() < 1e-9);
        CHECK(std::abs(est.disturbance_estimate()[0] - kForce) < 1e-3); // B_d = B: d̂ is the force
    }

    TEST_CASE("the target from r alone, and with a disturbance about a reference point") {
        const SteadyStateTarget target{design::steady_state_target(kPlantD)};
        const auto              p = target.equilibrium(ColVec<1>{{0.3}});
        CHECK(std::abs((kPlantD.C * p.x)(0, 0) - 0.3) < 1e-12);
        CHECK((kPlantD.A * p.x + kPlantD.B * p.u - p.x).norm() < 1e-12); // a steady state

        const SteadyStateTarget with_d{design::steady_state_target(kPlantD, design::InputDisturbance{})};
        const auto              q = with_d.equilibrium(ColVec<1>{{0.3}}, ColVec<1>{{0.1}});
        CHECK((kPlantD.A * q.x + kPlantD.B * (q.u + ColVec<1>{{0.1}}) - q.x).norm() < 1e-12); // the input cancels the load
        CHECK(std::abs((kPlantD.C * q.x)(0, 0) - 0.3) < 1e-12);
    }
}

namespace {

/// Two coupled masses with four readings: both positions read locally (rows 0, 1) and by the
/// motor encoders (rows 2, 3, with a little of the speed). Force on each mass.
constexpr StateSpace<4, 2, 4, double, 4, 4> four_readings() {
    StateSpace<4, 2, 4, double, 4, 4> s{};
    s.A(0, 1) = 1.0;
    s.A(1, 0) = -5.0;
    s.A(1, 1) = -0.4;
    s.A(1, 2) = 4.0;
    s.A(2, 3) = 1.0;
    s.A(3, 0) = 4.0;
    s.A(3, 2) = -6.0;
    s.A(3, 3) = -0.3;
    s.B(1, 0) = 1.0;
    s.B(3, 1) = 1.0;
    s.C(0, 0) = 1.0;
    s.C(1, 2) = 1.0;
    s.C(2, 0) = 1.0;
    s.C(2, 1) = 0.1;
    s.C(3, 2) = 1.0;
    s.C(3, 3) = 0.05;
    return s; // G = I, H = I by default
}

/**
 * The textbook offset-free law (test oracle): the steady-state gain on the absolute [x; d], the
 * target for r under d̂ about the reference point (x*, u*), and u = u_t + u* − K·(x̂ − x_t − x*).
 * OFLQG runs the same filter in deviations from that target instead, carrying the
 * reference point's motion between anchors; in exact arithmetic the two are one.
 */
template<size_t NXK, size_t NU, size_t NY, size_t NW, size_t NV, size_t ND, size_t NZ>
struct TextbookLaw {
    static constexpr size_t NX = NXK - ND;

    FixedGainObserver<NXK, NU, NY, double>    kf;
    SteadyStateTarget<NX, NU, ND, NZ, double> target;
    Matrix<NU, NX>                            K;

    TextbookLaw(const design::KalmanResult<NXK, NU, NY, double, NW, NV, ND>& kalman, const design::SteadyStateTargetResult<NX, NU, ND, NZ, double>& tg, const Matrix<NU, NX>& K_)
        : kf(kalman), target(tg), K(K_) {}

    /// One tick: the filter on the input the plant held over the last period (and holds now), then the target and the gain.
    ColVec<NU> step(const ColVec<NU>& u_applied, const ColVec<NZ>& r, const ColVec<NX>& x_ref, const ColVec<NU>& u_ref, const ColVec<NY>& y) {
        kf.predict(u_applied);
        kf.update(y, u_applied);
        const auto p = target.equilibrium(r, disturbance_estimate());
        return ColVec<NU>(p.u + u_ref - (K * ColVec<NX>(state_estimate() - p.x - x_ref)));
    }

    [[nodiscard]] ColVec<NX> state_estimate() const { return lead<NX>(kf.state()); }

    [[nodiscard]] ColVec<ND> disturbance_estimate() const {
        ColVec<ND> d{};
        for (size_t i = 0; i < ND; ++i) {
            d(i) = kf.state()(NX + i);
        }
        return d;
    }
};

/// The textbook law of an offset-free design result: its filter, target and gain.
template<size_t NX, size_t NU, size_t NY, size_t NXK, size_t NW, size_t NV, size_t ND, size_t NZ>
TextbookLaw<NXK, NU, NY, NW, NV, ND, NZ> textbook_law(const design::StateFeedbackResult<NX, NU, double, 0, NY, design::KalmanResult<NXK, NU, NY, double, NW, NV, ND>, NZ>& art) {
    static_assert(NXK == NX + ND, "The filter's model is the law's, then the disturbance.");
    return TextbookLaw<NXK, NU, NY, NW, NV, ND, NZ>{art.kalman, art.target, art.K};
}

/// The worst gap over a run, relative to the largest value seen.
struct Gap {
    double worst = 0.0;
    double scale = 0.0;
    template<typename A, typename B>
    void add(const A& a, const B& b) {
        worst = std::max(worst, (a - b).norm());
        scale = std::max(scale, b.norm());
    }
    [[nodiscard]] double relative() const { return worst / std::max(scale, 1e-300); }
};

/// The true plant: a model run exactly, with a load at its inputs.
template<size_t N>
struct Truth {
    StateSpace<N, 2, 4, double, 4, 4> m{};
    ColVec<N>                         x{};
    [[nodiscard]] ColVec<4>           y(const ColVec<2>& u) const { return (m.C * x) + (m.D * u); }
    void                              step(const ColVec<2>& u, const ColVec<2>& load) { x = (m.A * x) + (m.B * ColVec<2>(u + load)); }
};

const Matrix<4, 4> kQn4 = Matrix<4, 4>::identity() * 1e-4;
const Matrix<4, 4> kRn4{{1e-6, 0.0, 0.0, 0.0}, {0.0, 1e-6, 0.0, 0.0}, {0.0, 0.0, 1e-5, 0.0}, {0.0, 0.0, 0.0, 1e-5}};
const Matrix<2, 2> kQd2 = Matrix<2, 2>::identity() * 1e-2;
const Matrix<4, 4> kQ4{{1e3, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 1e3, 0.0}, {0.0, 0.0, 0.0, 0.0}};
const Matrix<6, 6> kQ6 = [] {
    Matrix<6, 6> q{}; // kQ4 on [x; v(k−1)], the registers unweighted
    q.block<4, 4>(0, 0) = kQ4;
    return q;
}();
const Matrix<2, 2> kR2 = Matrix<2, 2>::identity() * 1e-3;

} // namespace

TEST_SUITE("FixedGainObserver: the anchored law against the textbook form") {
    TEST_CASE("an input delay and a load at both inputs: the textbook law's estimates and commands") {
        //! The corexy shape: the commands one period late, a load at both inputs; the model is the plant with its registers
        auto plant = discretize(four_readings(), kTs).value();
        plant.G = Matrix<4, 4>::identity();
        const auto                                 whole = discretize(plant, kTs, InputDelay{kTs}).value(); // [x; v(k−1)]
        const design::DisturbanceModel<6, 4, 2, 2> load_model{.Bd = padded<2>(plant.B), .outputs = {{0, 1}}};
        const auto                                 art = design::oflqg(whole, {.Q = kQ6, .R = kR2}, design::KalmanNoise{kQn4, kRn4, kQd2}, load_model);
        REQUIRE(art.success);
        OFLQG law{art};
        static_assert(decltype(law.estimator)::NXE == 6 && decltype(law.disturbance_estimate())::rows() == 2);
        auto oracle = textbook_law(art);

        //! the truth: the register form, built apart from discretize; the load enters with the command
        Truth<6>        t{test::register_form<1>(four_readings(), kTs, kTs)};
        const ColVec<2> load{0.3, -0.2};
        ColVec<2>       u{};
        Gap             du;
        Gap             dx;
        Gap             dd;
        constexpr int   n = 600;
        for (int k = 0; k < n; ++k) {
            const ColVec<2> r{((k / 150) % 2 == 0) ? 0.2 : -0.1, 0.05 * std::sin(0.02 * k)};
            const ColVec<4> y = t.y(u);
            const ColVec<2> u_law = law.step(r, y);
            const ColVec<2> u_book = oracle.step(u, r, ColVec<6>{}, ColVec<2>{}, y);

            du.add(u_law, u_book);
            dx.add(law.state_estimate(), oracle.state_estimate());
            dd.add(law.disturbance_estimate(), oracle.disturbance_estimate());
            u = u_law;
            t.step(u, (k > n / 3) ? load : ColVec<2>{});
        }
        MESSAGE("commands ", du.relative(), ", estimate (registers included) ", dx.relative(), ", disturbance ", dd.relative(), ", |d̂ − load| ", (law.disturbance_estimate() - load).norm());
        CHECK(du.relative() < 1e-12);
        CHECK(dx.relative() < 1e-12);
        CHECK(dd.relative() < 1e-12);
        CHECK((law.disturbance_estimate() - load).norm() < 1e-3); // B_d = Γ on the plant: d̂ is the load
        CHECK(t.x.norm() < 1.0);                                  // the loop held
    }

    TEST_CASE("feedthrough, an input disturbance, no input delay, about a moving reference point") {
        //! The anchor drifts with the reference point, which is no trajectory of the model
        auto pc = four_readings();
        pc.D(2, 0) = 0.2;
        pc.D(0, 1) = 0.1;
        auto plant = discretize(pc, kTs).value();
        plant.G = Matrix<4, 4>::identity();
        const auto art = design::oflqg(plant, {.Q = kQ4, .R = kR2}, design::KalmanNoise{kQn4, kRn4, kQd2}, design::InputDisturbance<2>{{0, 1}});
        REQUIRE(art.success);
        OFLQG law{art};
        auto  oracle = textbook_law(art);

        Truth<4>      t{plant};
        ColVec<2>     u{};
        Gap           du;
        Gap           dx;
        constexpr int n = 400;
        for (int k = 0; k < n; ++k) {
            const double    s = 0.1 * std::sin(0.5 * k * kTs);
            const double    c = 0.05 * std::cos(0.5 * k * kTs);
            const ColVec<4> x_ref{s, c, -s, -c};
            const ColVec<2> u_ref{s, -2.0 * s};
            const ColVec<4> y = t.y(u);
            const ColVec<2> u_law = law.step({}, x_ref, u_ref, y);
            const ColVec<2> u_book = oracle.step(u, ColVec<2>{}, x_ref, u_ref, y);

            du.add(u_law, u_book);
            dx.add(law.state_estimate(), oracle.state_estimate());
            u = u_law;
            t.step(u, (k > n / 3) ? ColVec<2>{0.2, 0.1} : ColVec<2>{});
        }
        MESSAGE("commands ", du.relative(), ", estimate ", dx.relative());
        CHECK(du.relative() < 1e-12);
        CHECK(dx.relative() < 1e-12);
        CHECK(t.x.norm() < 1.0); // the loop held
    }

    TEST_CASE("no disturbance, a fraction of a period of input delay: the register known, the plant's estimate the undelayed filter's") {
        constexpr double theta = 0.6 * kTs;
        auto             whole = discretize(four_readings(), kTs, InputDelay{theta}).value(); // [x; v(k−1)], the hold split at θ
        whole.G = padded<2>(Matrix<4, 4>::identity());                                        // the process noise on the plant's states
        const auto art = design::lqg(whole, {.Q = kQ6, .R = kR2}, design::KalmanNoise{.process = kQn4, .measurement = kRn4});
        REQUIRE(art.success);
        LQG law{art};

        //! The oracle: the plain filter on the undelayed plant, fed what the plant receives over each
        //! period, Γ₀u(k) + Γ₁u(k−1), with Γ₀, Γ₁ from the register form (built apart from discretize)
        auto plant = discretize(four_readings(), kTs).value();
        plant.G = Matrix<4, 4>::identity();
        const auto plain = design::kalman(plant, {.process = kQn4, .measurement = kRn4});
        REQUIRE(plain.success);
        const auto         truth = test::register_form<1>(four_readings(), kTs, theta);
        const Matrix<4, 2> gamma0 = truth.B.block<4, 2>(0, 0);
        const Matrix<4, 2> gamma1 = truth.A.block<4, 2>(0, 4);

        Truth<6>      t{truth, ColVec<6>{0.1, 0.0, -0.1}};
        ColVec<4>     x_book{};
        ColVec<2>     u{};        // u(k−1)
        ColVec<2>     u_before{}; // u(k−2)
        Gap           dx;
        Gap           dv;
        constexpr int n = 300;
        for (int k = 0; k < n; ++k) {
            const ColVec<4> y = t.y(u);
            const ColVec<2> u_law = law.step(y);
            x_book = (plain.sys.A * x_book) + (gamma0 * u) + (gamma1 * u_before);
            x_book = x_book + (plain.M * (y - (plain.sys.C * x_book)));

            const ColVec<6> est = law.state_estimate();
            dx.add(lead<4>(est), x_book);
            dv.add(ColVec<2>{est(4), est(5)}, u); // the register: the command sent last period
            u_before = u;
            u = u_law;
            t.step(u, ColVec<2>{});
        }
        MESSAGE("plant estimate ", dx.relative(), ", register ", dv.relative());
        CHECK(dx.relative() < 1e-12);
        CHECK(dv.relative() < 1e-12);
        CHECK(std::abs(t.x[0]) < 1e-3); // regulated
    }
}
