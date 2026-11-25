// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// design::observer: the placed filter, design::kalman's counterpart, with the filter designs'
// options, on the Kalman design's result and runtime.

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/law.hpp"
#include "damp/controllers/state_feedback.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/pole_placement.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/fixed_gain_observer.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/estimation/observer.hpp"
#include "damp/math/complex.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

using Cplx = damp::complex<double>;

constexpr double kTs = 0.01;

/// A double integrator: position read, the input a force.
constexpr StateSpace kCart{.A = Matrix{{0.0, 1.0}, {0.0, 0.0}}, .B = Matrix{{0.0}, {1.0}}, .C = Matrix{{1.0, 0.0}}};

/// The same plant sampled at kTs.
constexpr auto kCartD = *discretize(kCart, kTs);

/// Two lags: a bias on the reading is observable here (on the double integrator it is one mode with the position).
constexpr StateSpace kLags{.A = Matrix{{-1.0, 0.5}, {0.0, -2.0}}, .B = Matrix{{0.0}, {1.0}}, .C = Matrix{{1.0, 0.0}}};

/// A three-state, two-input, two-output plant.
constexpr StateSpace kTwoInput{
    .A = Matrix{{0.0, 1.0, 0.0}, {0.0, -0.5, 1.0}, {0.0, 0.0, -2.0}},
    .B = Matrix{{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}},
    .C = Matrix{{1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}},
};

constexpr design::LQRWeights<2, 1> kWeights{.Q = Matrix{{10.0, 0.0}, {0.0, 1.0}}, .R = {0.1}};

/// Two designs of the same result type, cv aside.
template<typename A, typename B>
inline constexpr bool same_type = std::is_same_v<std::remove_cvref_t<A>, std::remove_cvref_t<B>>;

/// A quiet NaN from its bits, read through a volatile so the compiler cannot fold it.
double runtime_nan() {
    volatile std::uint64_t bits = 0x7ff8000000000000ull;
    return std::bit_cast<double>(static_cast<std::uint64_t>(bits));
}

/// -Inf from its bits, read through a volatile.
double runtime_minus_inf() {
    volatile std::uint64_t bits = 0xfff0000000000000ull;
    return std::bit_cast<double>(static_cast<std::uint64_t>(bits));
}

/// Every pole of @p want is an eigenvalue of @p A, to @p tol, each claimed once.
template<size_t N>
bool has_spectrum(const Matrix<N, N>& A, const std::array<Cplx, N>& want, double tol) {
    const auto ev = mat::compute_eigenvalues(A);
    if (!ev.success) {
        return false;
    }
    std::array<bool, N> claimed{};
    for (size_t i = 0; i < N; ++i) {
        size_t best = N;
        double best_d = 0.0;
        for (size_t j = 0; j < N; ++j) {
            const double d = (Cplx(ev.values[j]) - want[i]).abs();
            if (!claimed[j] && (best == N || d < best_d)) {
                best = j;
                best_d = d;
            }
        }
        if (best == N || !(best_d <= tol)) {
            return false;
        }
        claimed[best] = true;
    }
    return true;
}

/// ‖Π (A − pᵢ I)‖ relative to ‖A‖^N: zero when the poles are A's characteristic roots with their
/// multiplicities (Cayley–Hamilton), however defective; real or conjugate-paired poles.
template<size_t N>
double characteristic_residual(const Matrix<N, N>& A, const std::array<Cplx, N>& poles) {
    Matrix<N, N> re = Matrix<N, N>::identity(); // the real part of the running product
    Matrix<N, N> im{};                          // and its imaginary part
    for (size_t i = 0; i < N; ++i) {
        const Matrix<N, N> shifted = A - (Matrix<N, N>::identity() * poles[i].real());
        const Matrix<N, N> nre = (re * shifted) + (im * poles[i].imag());
        const Matrix<N, N> nim = (im * shifted) - (re * poles[i].imag());
        re = nre;
        im = nim;
    }
    const double scale = std::pow(1.0 + A.norm(), static_cast<double>(N));
    return (re.norm() + im.norm()) / scale;
}

/// The poles mapped to z = e^{s·Ts}.
template<size_t N>
std::array<Cplx, N> to_z(const std::array<Cplx, N>& s, double Ts) {
    std::array<Cplx, N> z{};
    for (size_t i = 0; i < N; ++i) {
        z[i] = design::detail::pole_s_to_z(s[i], Ts);
    }
    return z;
}

/// The error dynamics of the current-form filter a design runs: (I − M·C)·A on its model.
template<typename Filter>
auto error_dynamics(const Filter& kf) {
    constexpr size_t NXK = decltype(kf.sys.A)::rows();
    return Matrix<NXK, NXK>((Matrix<NXK, NXK>::identity() - (kf.M * kf.sys.C)) * kf.sys.A);
}

/// The poles in a function template: bare braces deduce their count there as at namespace scope.
template<typename T>
constexpr auto observer_in_template(T wo) {
    return design::observer(kCart, kTs, {-wo, -1.2 * wo, -1.4 * wo}, design::InputDisturbance{});
}

/// A constexpr placed filter on [x; d], deployed on its constants.
constexpr auto kObserved = design::observer(kCart, kTs, {-30.0, -35.0, -40.0}, design::InputDisturbance{});

} // namespace

TEST_SUITE("design::observer: the placed filter") {
    TEST_CASE("it is state_feedback's placed filter, and without options place_observer's gain, in the Kalman design's type") {
        //! With a disturbance spec: the filter of state_feedback with ObserverPoles, the same model and type.
        const auto res = design::observer(kCart, kTs, {-30.0, -35.0, -40.0}, design::InputDisturbance{});
        REQUIRE(res.success);
        const auto law = design::state_feedback(kCart, kTs, design::Poles{-8.0, -12.0}, design::ObserverPoles{-30.0, -35.0, -40.0}, design::InputDisturbance{});
        REQUIRE(law.success);
        static_assert(same_type<decltype(res), decltype(law.kalman)>);
        CHECK((res.M - law.kalman.M).norm() <= 1e-10 * law.kalman.M.norm()); // M reaches 250: the placement's rounding grows with it
        CHECK((res.L - law.kalman.L).norm() <= 1e-10 * law.kalman.L.norm());
        CHECK((res.sys.A - law.kalman.sys.A).norm() <= 1e-15 * res.sys.A.norm());
        CHECK((res.sys.C - law.kalman.sys.C).norm() <= 1e-15 * res.sys.C.norm());
        const auto kf = design::kalman(kCart, kTs, design::KalmanNoise{.process = Matrix<2, 2>::identity(), .measurement = 0.01, .disturbance = 1.0}, design::InputDisturbance{});
        static_assert(same_type<decltype(res), decltype(kf)>); // design::kalman's result, on its model
        CHECK((res.sys.A - kf.sys.A).norm() <= 1e-15 * res.sys.A.norm());
        CHECK((res.sys.G - kf.sys.G).norm() == 0.0); // G = I, written
        CHECK(res.P.norm() == 0.0);                  // no noise: no covariance
        CHECK(res.Q.norm() == 0.0);

        //! Without options: design::place_observer's M and L on the ZOH model.
        const auto plain = design::observer(kCart, kTs, {-30.0, -35.0});
        const auto bare = design::place_observer(kCart, {-30.0, -35.0}, kTs);
        REQUIRE(plain.success);
        REQUIRE(bare.success);
        CHECK((plain.M - bare.M).norm() <= 1e-12 * bare.M.norm());
        CHECK((plain.L - bare.L).norm() <= 1e-12 * bare.L.norm());
        CHECK((plain.sys.A - bare.sys.A).norm() <= 1e-15 * bare.sys.A.norm());
        CHECK((plain.sys.B - bare.sys.B).norm() <= 1e-15 * bare.sys.B.norm());
        static_assert(same_type<decltype(plain), decltype(design::kalman(kCart, kTs, design::KalmanNoise{.process = Matrix<2, 2>::identity(), .measurement = 0.01}))>);

        //! A discrete plant: z-plane poles, place_observer's gain on it.
        const auto discrete = design::observer(kCartD, {0.6, 0.65});
        const auto discrete_bare = design::place_observer(kCartD, {0.6, 0.65});
        REQUIRE(discrete.success);
        CHECK((discrete.M - discrete_bare.M).norm() <= 1e-12 * discrete_bare.M.norm());
        CHECK(discrete.sys.Ts == kTs);
    }

    TEST_CASE("the error poles are where they are placed: s-plane on the sampled form, z-plane on a discrete plant") {
        const auto res = design::observer(kCart, kTs, {-30.0, -35.0, -40.0}, design::InputDisturbance{});
        REQUIRE(res.success);
        CHECK(has_spectrum(error_dynamics(res), to_z(std::array<Cplx, 3>{Cplx{-30.0, 0.0}, Cplx{-35.0, 0.0}, Cplx{-40.0, 0.0}}, kTs), 1e-10));
        CHECK((res.L - (res.sys.A * res.M)).norm() <= 1e-15 * res.L.norm() + 1e-15); // L = A·M

        //! A conjugate pair.
        const auto pair = design::observer(kCart, kTs, {damp::complex{-30.0, 20.0}, damp::complex{-30.0, -20.0}});
        REQUIRE(pair.success);
        CHECK(has_spectrum(error_dynamics(pair), to_z(std::array<Cplx, 2>{Cplx{-30.0, 20.0}, Cplx{-30.0, -20.0}}, kTs), 1e-12));

        //! z-plane on a discrete plant, with a bias on the reading.
        const auto lags_d = *discretize(kLags, kTs);
        const auto biased = design::observer(lags_d, {0.6, 0.65, 0.7}, design::OutputDisturbance{});
        REQUIRE(biased.success);
        CHECK(has_spectrum(error_dynamics(biased), {Cplx{0.6, 0.0}, Cplx{0.65, 0.0}, Cplx{0.7, 0.0}}, 1e-8)); // M reaches 475 on modes near z = 1: the placement's rounding

        //! Two readings, KNV: the multi-output placement.
        const auto two = design::observer(kTwoInput, kTs, {-20.0, -25.0, -30.0});
        REQUIRE(two.success);
        CHECK(has_spectrum(error_dynamics(two), to_z(std::array<Cplx, 3>{Cplx{-20.0, 0.0}, Cplx{-25.0, 0.0}, Cplx{-30.0, 0.0}}, kTs), 1e-12));
    }

    TEST_CASE("bare braces and ObserverPoles give the same design, a list in another scalar converted") {
        const auto braces = design::observer(kCart, kTs, {-30.0, -35.0, -40.0}, design::InputDisturbance{});
        const auto spelled = design::observer(kCart, kTs, design::ObserverPoles{-30.0, -35.0, -40.0}, design::InputDisturbance{});
        const auto single = design::observer(kCart, kTs, design::ObserverPoles{-30.0f, -35.0f, -40.0f}, design::InputDisturbance{});
        static_assert(same_type<decltype(braces), decltype(spelled)>);
        static_assert(same_type<decltype(braces), decltype(single)>); // the plant's scalar
        CHECK((braces.M - spelled.M).norm() <= 1e-12 * braces.M.norm());
        CHECK((braces.M - single.M).norm() <= 1e-10 * braces.M.norm()); // -35.0f and -40.0f are exact, so the same poles
        const auto discrete = design::observer(kCartD, {0.6, 0.65});
        const auto discrete_spelled = design::observer(kCartD, design::ObserverPoles{0.6, 0.65});
        static_assert(same_type<decltype(discrete), decltype(discrete_spelled)>);
        CHECK((discrete.M - discrete_spelled.M).norm() <= 1e-12 * discrete.M.norm());

        //! In a function template, at compile time.
        constexpr double wo = 30.0;
        constexpr auto   in_template = observer_in_template(wo);
        static_assert(in_template.success);
        const auto at_run_time = design::observer(kCart, kTs, {-wo, -1.2 * wo, -1.4 * wo}, design::InputDisturbance{});
        CHECK((in_template.M - at_run_time.M).norm() <= 1e-10 * at_run_time.M.norm());
    }

    TEST_CASE("delays: the delay states' error poles at 0, M with no rows on the input's") {
        const auto zp = to_z(std::array<Cplx, 2>{Cplx{-30.0, 0.0}, Cplx{-35.0, 0.0}}, kTs);
        //! An InputDelay: its delay states hold the inputs sent, which the readings do not correct.
        const auto in = design::observer(kCart, kTs, {-30.0, -35.0}, InputDelay<2>{0.015});
        REQUIRE(in.success);
        CHECK(characteristic_residual(error_dynamics(in), {zp[0], zp[1], Cplx{}, Cplx{}}) <= 1e-12);
        CHECK(in.M(2, 0) == 0.0);
        CHECK(in.M(3, 0) == 0.0);
        const auto in_law = design::state_feedback(kCart, kTs, kWeights, design::ObserverPoles{-30.0, -35.0}, InputDelay<2>{0.015});
        static_assert(same_type<decltype(in), decltype(in_law.kalman)>);
        CHECK((in.M - in_law.kalman.M).norm() <= 1e-12 * in.M.norm());

        //! With an InputDisturbance d rides through the delay states: they take d's correction.
        const auto zd = to_z(std::array<Cplx, 3>{Cplx{-30.0, 0.0}, Cplx{-35.0, 0.0}, Cplx{-40.0, 0.0}}, kTs);
        const auto of = design::observer(kCart, kTs, {-30.0, -35.0, -40.0}, design::InputDisturbance{}, InputDelay<2>{0.015});
        REQUIRE(of.success);
        CHECK(characteristic_residual(error_dynamics(of), {zd[0], zd[1], zd[2], Cplx{}, Cplx{}}) <= 1e-12);
        CHECK(of.M(2, 0) == doctest::Approx(of.M(4, 0)).epsilon(1e-12));
        CHECK(of.M(3, 0) == doctest::Approx(of.M(4, 0)).epsilon(1e-12));

        //! A late reading: its delay state is placed at 0 with the rest, on design::kalman's model.
        const auto late = design::observer(kCart, kTs, {-30.0, -35.0}, OutputDelay{0.005});
        REQUIRE(late.success);
        CHECK(characteristic_residual(error_dynamics(late), {zp[0], zp[1], Cplx{}}) <= 1e-12);
        const auto late_kf = design::kalman(kCart, kTs, design::KalmanNoise{.process = Matrix<2, 2>::identity(), .measurement = 0.01}, OutputDelay{0.005});
        static_assert(same_type<decltype(late), decltype(late_kf)>);
        CHECK((late.sys.A - late_kf.sys.A).norm() <= 1e-15 * late.sys.A.norm());

        //! A bias on a late reading: d reaches the reading as it was when it was taken, as design::kalman's model has it.
        const auto late_bias = design::observer(kLags, kTs, {-30.0, -35.0, -40.0}, design::OutputDisturbance{}, OutputDelay{0.005});
        const auto late_bias_kf = design::kalman(kLags, kTs, design::KalmanNoise{.process = Matrix<2, 2>::identity(), .measurement = 0.01, .disturbance = 1.0}, design::OutputDisturbance{}, OutputDelay{0.005});
        REQUIRE(late_bias.success);
        REQUIRE(late_bias_kf.success);
        static_assert(same_type<decltype(late_bias), decltype(late_bias_kf)>);
        CHECK((late_bias.sys.A - late_bias_kf.sys.A).norm() <= 1e-15 * late_bias.sys.A.norm());
        CHECK((late_bias.sys.C - late_bias_kf.sys.C).norm() <= 1e-15 * late_bias.sys.C.norm());
        const auto zb = to_z(std::array<Cplx, 3>{Cplx{-30.0, 0.0}, Cplx{-35.0, 0.0}, Cplx{-40.0, 0.0}}, kTs);
        CHECK(characteristic_residual(error_dynamics(late_bias), {zb[0], zb[1], zb[2], Cplx{}}) <= 1e-12);
        //! A law's placed filter is on the same model.
        const auto late_bias_law = design::state_feedback(kLags, kTs, kWeights, design::ObserverPoles{-30.0, -35.0, -40.0}, design::OutputDisturbance{}, OutputDelay{0.005});
        REQUIRE(late_bias_law.success);
        CHECK((late_bias_law.kalman.sys.A - late_bias_kf.sys.A).norm() <= 1e-15 * late_bias.sys.A.norm());
        CHECK((late_bias_law.kalman.M - late_bias.M).norm() <= 1e-12 * late_bias.M.norm());
    }

    TEST_CASE("a measured disturbance: known inputs after u, the gain the same, the estimate exact under it") {
        constexpr design::MeasuredDisturbance measured{Matrix{{0.0}, {1.0}}};
        const auto                            res = design::observer(kCart, kTs, {-30.0, -35.0}, measured);
        const auto                            plain = design::observer(kCart, kTs, {-30.0, -35.0});
        REQUIRE(res.success);
        static_assert(decltype(res.sys.B)::cols() == 2); // [u; d]
        CHECK((res.M - plain.M).norm() <= 1e-12 * plain.M.norm());
        const auto both = *discretize(StateSpace{.A = kCart.A, .B = Matrix{{0.0, 0.0}, {1.0, 1.0}}, .C = kCart.C}, kTs);
        CHECK((res.sys.B - both.B).norm() <= 1e-15);
        CHECK(res.sys.D(0, 1) == 0.0);
        const auto kf = design::kalman(kCart, kTs, design::KalmanNoise{.process = Matrix<2, 2>::identity(), .measurement = 0.01}, measured);
        static_assert(same_type<decltype(res), decltype(kf)>);
        CHECK((res.sys.B - kf.sys.B).norm() <= 1e-15);

        //! The observer predicts with [u; d], the d the plant held over the period: its error decays at the poles.
        FixedGainObserver obs{res};
        static_assert(std::is_same_v<decltype(obs), FixedGainObserver<2, 2, 1, double, 0>>);
        ColVec    x{0.2, -0.1};
        ColVec<2> held{};
        double    worst = 0.0;
        for (int k = 0; k < 300; ++k) {
            const auto xhat = obs.step(kCart.C * x, held);
            if (k >= 200) {
                worst = std::max(worst, (xhat - x).norm() / (1.0 + x.norm()));
            }
            const double t = static_cast<double>(k);
            const ColVec ud{std::sin(0.05 * t), 0.7 + (0.3 * std::cos(0.02 * t))};
            x = (both.A * x) + (both.B * ud);
            held = ud;
        }
        CHECK(worst <= 1e-9);
    }

    TEST_CASE("it deploys on FixedGainObserver and on Law<kRes, float> of a constexpr design") {
        static_assert(kObserved.success);
        FixedGainObserver obs{kObserved};
        static_assert(std::is_same_v<decltype(obs), FixedGainObserver<3, 1, 1, double, 1>>);
        static_assert(std::is_same_v<Law<kObserved, float>::Runtime, FixedGainObserver<3, 1, 1, float, 1>>);

        //! A load at the input the observer is not told of, the command taking most of it back so the state stays near 0.
        static constinit Law<kObserved, float> law;
        FixedGainObserver                      held{kObserved.as<float>()};
        ColVec                                 x{0.0, 0.0};
        ColVec                                 x_read{0.0, 0.0};
        constexpr double                       load = 0.5;
        double                                 apart = 0.0;
        double                                 u_held = 0.0; // the input the plant held over the last period
        for (int k = 0; k < 400; ++k) {
            const ColVec y{static_cast<float>(x[0])};
            const ColVec uf{static_cast<float>(u_held)};
            static_cast<void>(obs.step(ColVec{x[0]}, ColVec{u_held}));
            const auto a = law.step(y, uf);
            const auto b = held.step(y, uf);
            apart = std::max(apart, static_cast<double>((a - b).norm()) / (1.0 + static_cast<double>(b.norm())));
            x_read = x;
            u_held = -load + (0.3 * std::sin(0.03 * static_cast<double>(k)));
            x = (kCartD.A * x) + (kCartD.B * ColVec{u_held + load});
        }
        CHECK(apart <= 1e-4); // the same arithmetic in float, the zeros of the constants skipped
        CHECK(obs.disturbance_estimate()[0] == doctest::Approx(load).epsilon(1e-9));
        CHECK(obs.state_estimate()[0] == doctest::Approx(x_read[0]).epsilon(1e-9));
        CHECK(law.disturbance_estimate()[0] == doctest::Approx(load).epsilon(1e-2)); // float: the reading's rounding over Ts²/2
    }

    TEST_CASE("refusals at run time: a pole that is not finite, an unobservable mode, a plant in the other domain") {
        //! Not finite: NaN, and −∞ in the s-plane, which z = e^{s·Ts} would map to a finite 0.
        for (const double bad : {runtime_nan(), runtime_minus_inf()}) {
            CHECK(design::observer(kCart, kTs, {-30.0, bad}).status == design::DesignStatus::invalid_parameter);
            CHECK(design::observer(kCart, kTs, {bad, -35.0, -40.0}, design::InputDisturbance{}).status == design::DesignStatus::invalid_parameter);
            CHECK(design::observer(kCartD, {0.6, bad}).status == design::DesignStatus::invalid_parameter);
            CHECK(design::observer(kCart, kTs, design::ObserverPoles{bad, -35.0}, InputDelay{0.004}).status == design::DesignStatus::invalid_parameter);
        }
        //! A mode the reading does not see.
        constexpr StateSpace split{.A = Matrix{{-1.0, 0.0}, {0.0, -2.0}}, .B = Matrix{{1.0}, {1.0}}, .C = Matrix{{1.0, 0.0}}};
        CHECK(design::observer(split, kTs, {-30.0, -35.0}).status == design::DesignStatus::not_placeable);
        //! A bias on the double integrator's reading is one mode with the position.
        CHECK(design::observer(kCart, kTs, {-30.0, -35.0, -40.0}, design::OutputDisturbance{}).status == design::DesignStatus::not_placeable);
        //! The discrete form takes a discrete plant, the sampled one a continuous plant and a usable Ts.
        CHECK(design::observer(kCart, {0.6, 0.65}).status == design::DesignStatus::invalid_sample_time);
        CHECK(design::observer(kCartD, kTs, {-30.0, -35.0}).status == design::DesignStatus::invalid_sample_time);
        CHECK(design::observer(kCart, 0.0, {-30.0, -35.0}).status == design::DesignStatus::invalid_sample_time);
        //! A failed design deploys as an open-loop predictor, valid() false.
        const FixedGainObserver refused{design::observer(split, kTs, {-30.0, -35.0})};
        CHECK_FALSE(refused.valid());

        //! At compile time (each one a static_assert):
        //!   design::observer(kCart, kTs, {-30.0, -35.0}, design::InputDisturbance{});        // observer_pole_count<3, 2>
        //!   design::observer(kCart, kTs, design::ObserverPoles{-30.0, -35.0, -40.0});      // observer_pole_count<2, 3>
        //!   design::observer(kCart, kTs, {-30.0, -35.0}, design::KalmanNoise{...});         // a KalmanNoise is design::kalman's
        //!   design::observer(kCart, kTs, {-30.0, -35.0}, design::ObserverPoles{...});       // the poles are the argument, not an option
        //!   design::observer(kCartD, {0.6, 0.65}, OutputDelay{0.005});                     // a delay is sampled with the plant
        //!   design::observer(kCart, kTs, {-30.0, -35.0}, design::IntegralAction{.Q = 1.0}); // an option of a type this design does not take
        //!   design::observer(kCart, kTs, {-30.0, -35.0, -40.0}, design::DisturbanceModel{...}); // not on the plant's states and outputs
    }
}
