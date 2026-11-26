// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// design::state_feedback with Poles in the weights' slot and ObserverPoles in place of a
// KalmanNoise: the placed gain and the placed filter, with every option, on the LQ runtimes.

#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/adrc.hpp"
#include "damp/controllers/detail/lq_cost.hpp"
#include "damp/controllers/law.hpp"
#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/controllers/oflqg.hpp"
#include "damp/controllers/state_feedback.hpp"
#include "damp/controllers/state_feedback_types.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/pole_placement.hpp"
#include "damp/design/status.hpp"
#include "damp/design/steady_state_target.hpp"
#include "damp/estimation/fixed_gain_observer.hpp"
#include "damp/estimation/kalman.hpp"
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

/// A three-state, two-input plant: the multi-input (KNV) path.
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

/// The poles in a function template: deduced there as at namespace scope (clang's llvm#114870 is about designated fields).
template<typename T>
constexpr auto placed_in_template(T wc) {
    return design::state_feedback(kCart, kTs, design::Poles{-wc, -1.5 * wc}, design::ObserverPoles{-4.0 * wc, -5.0 * wc});
}

constexpr auto kPlacedDiscrete = design::state_feedback(kCartD, design::Poles{0.9, 0.95});
constexpr auto kPlacedObserver = design::state_feedback(kCart, kTs, design::Poles{-8.0, -12.0}, design::ObserverPoles{-30.0, -35.0, -40.0}, design::InputDisturbance{});

} // namespace

TEST_SUITE("State feedback placed: Poles and ObserverPoles") {
    TEST_CASE("Poles deduce their count and scalar from bare braces, real entries converting, in a template too") {
        static_assert(std::is_same_v<decltype(design::Poles{-5.0, -5.0}), design::Poles<2, double>>);
        static_assert(std::is_same_v<decltype(design::Poles{Cplx{-1.0, 2.0}, Cplx{-1.0, -2.0}}), design::Poles<2, double>>);
        static_assert(std::is_same_v<decltype(design::Poles{-3.0, damp::complex{-1.0, 2.0}, damp::complex{-1.0, -2.0}}), design::Poles<3, double>>);
        static_assert(std::is_same_v<decltype(design::Poles{-5.0f}), design::Poles<1, float>>);
        static_assert(std::is_same_v<decltype(design::ObserverPoles{-40.0, -40.0, -40.0}), design::ObserverPoles<3, double>>);
        constexpr design::Poles p{-5.0, damp::complex{-1.0, 2.0}};
        static_assert(p.poles[0].real() == -5.0 && p.poles[0].imag() == 0.0 && p.poles[1].imag() == 2.0);
        constexpr auto in_template = placed_in_template(8.0);
        static_assert(in_template.success);
        CHECK(in_template.K(0, 0) == doctest::Approx(design::place(kCart, {-8.0, -12.0}, kTs).K(0, 0)).epsilon(1e-12));
    }

    TEST_CASE("a discrete plant takes z-plane poles: eig(A − BK) there, K design::place's") {
        static_assert(kPlacedDiscrete.success && kPlacedDiscrete.Ts == kTs);
        const auto gain = design::place(kCartD, {0.9, 0.95});
        REQUIRE(gain.success);
        for (size_t j = 0; j < 2; ++j) {
            CHECK(std::abs(kPlacedDiscrete.K(0, j) - gain.K(0, j)) <= 1e-12 * gain.K.norm());
        }
        const Matrix<2, 2> Acl = kCartD.A - (kCartD.B * kPlacedDiscrete.K);
        CHECK(has_spectrum(Acl, {Cplx{0.9, 0.0}, Cplx{0.95, 0.0}}, 1e-12));
        CHECK(has_spectrum(Acl, {kPlacedDiscrete.e[0], kPlacedDiscrete.e[1]}, 1e-12)); // e is the closed loop's
        CHECK(kPlacedDiscrete.S.norm() == 0.0);                                        // no Riccati equation
    }

    TEST_CASE("the sampled form takes s-plane poles in rad/s, z = e^{s·Ts}: design::place(sys, poles, Ts)'s gain") {
        const std::array<Cplx, 2> real{Cplx{-8.0, 0.0}, Cplx{-12.0, 0.0}};
        const std::array<Cplx, 2> pair{Cplx{-6.0, 8.0}, Cplx{-6.0, -8.0}};
        const auto                res = design::state_feedback(kCart, kTs, design::Poles{-8.0, -12.0});
        const auto                res_pair = design::state_feedback(kCart, kTs, design::Poles{damp::complex{-6.0, 8.0}, damp::complex{-6.0, -8.0}});
        REQUIRE(res.success);
        REQUIRE(res_pair.success);
        const auto gain = design::place(kCart, {-8.0, -12.0}, kTs);
        const auto gain_pair = design::place(kCart, {pair[0], pair[1]}, kTs);
        for (size_t j = 0; j < 2; ++j) {
            CHECK(std::abs(res.K(0, j) - gain.K(0, j)) <= 1e-12 * gain.K.norm());
            CHECK(std::abs(res_pair.K(0, j) - gain_pair.K(0, j)) <= 1e-12 * gain_pair.K.norm());
        }
        CHECK(has_spectrum((kCartD.A - (kCartD.B * res.K)), to_z(real, kTs), 1e-12));
        CHECK(has_spectrum((kCartD.A - (kCartD.B * res_pair.K)), to_z(pair, kTs), 1e-12));

        //! Several inputs: KNV, design::place's gain on the same pair.
        const auto two = design::state_feedback(kTwoInput, kTs, design::Poles{-2.0, -3.0, -4.0});
        const auto two_gain = design::place(kTwoInput, {-2.0, -3.0, -4.0}, kTs);
        REQUIRE(two.success);
        CHECK((two.K - two_gain.K).norm() <= 1e-12 * two_gain.K.norm());
        const auto two_d = *discretize(kTwoInput, kTs);
        CHECK(has_spectrum((two_d.A - (two_d.B * two.K)), to_z(std::array<Cplx, 3>{Cplx{-2.0, 0.0}, Cplx{-3.0, 0.0}, Cplx{-4.0, 0.0}}, kTs), 1e-12));
    }

    TEST_CASE("with TrackedOutputs the target is the LQ design's: it is the model's, not the gain's") {
        const auto placed = design::state_feedback(kCart, kTs, design::Poles{-8.0, -12.0}, design::TrackedOutputs{{0}});
        const auto lq = design::state_feedback(kCart, kTs, kWeights, design::TrackedOutputs{{0}});
        REQUIRE(placed.success);
        REQUIRE(lq.success);
        static_assert(same_type<decltype(placed), decltype(lq)>);
        CHECK((placed.target.target_r - lq.target.target_r).norm() <= 1e-12);

        //! The law holds r: at the target the command is the target's input.
        LQR law{placed};
        law.update(ColVec{0.5, 0.0});
        CHECK(law.control(ColVec{0.5})[0] == doctest::Approx(0.0).epsilon(1e-12));
    }

    TEST_CASE("with an IntegralAction the poles cover the integrators: eig of the integral model at them") {
        //! Its .Q is the LQ design's cost on the integrators: it sizes them, and a placed design reads no cost.
        const auto res = design::state_feedback(kCart, kTs, design::Poles{-8.0, -12.0, -4.0}, design::IntegralAction{.Q = 80.0});
        REQUIRE(res.success);
        static_assert(same_type<decltype(res), decltype(design::state_feedback(kCart, kTs, kWeights, design::IntegralAction{.Q = 80.0}))>);
        const auto im = design::detail::integral_model(kCartD.A, kCartD.B, kCartD.C, kCartD.D, Matrix{{1.0}}, kTs);
        const auto Acl = (im.A - (im.B * res.K));
        const auto z = to_z(std::array<Cplx, 3>{Cplx{-8.0, 0.0}, Cplx{-12.0, 0.0}, Cplx{-4.0, 0.0}}, kTs);
        CHECK(has_spectrum(Acl, z, 1e-12));
        const auto gain = design::place(StateSpace{.A = im.A, .B = im.B, .C = Matrix<1, 3>{}, .Ts = kTs}, {z[0], z[1], z[2]});
        CHECK((res.K - gain.K).norm() <= 1e-12 * gain.K.norm());
        CHECK(res.integral.step == kTs);

        //! The LQI runs it, and holds the output at r.
        LQI        law{res};
        StateSpace plant = kCartD;
        ColVec<2>  x{};
        for (int k = 0; k < 2000; ++k) {
            const auto u = law.step(ColVec{1.0}, x);
            x = (plant.A * x) + (plant.B * (u + ColVec{0.3})); // a load the integrator takes up
        }
        CHECK(x[0] == doctest::Approx(1.0).epsilon(1e-9));
    }

    TEST_CASE("with an InputDelay the delay states are deadbeat: their poles at 0, on top of the list") {
        const auto res = design::state_feedback(kCart, kTs, design::Poles{-8.0, -12.0}, InputDelay<2>{0.015});
        REQUIRE(res.success);
        static_assert(same_type<decltype(res), decltype(design::state_feedback(kCart, kTs, kWeights, InputDelay<2>{0.015}))>);
        const auto                model = *discretize(kCart, kTs, InputDelay<2>{0.015});
        const auto                Acl = (model.A - (model.B * res.K));
        const auto                z = to_z(std::array<Cplx, 2>{Cplx{-8.0, 0.0}, Cplx{-12.0, 0.0}}, kTs);
        const std::array<Cplx, 4> want{z[0], z[1], Cplx{}, Cplx{}};
        CHECK(characteristic_residual(Acl, want) <= 1e-12);
        CHECK(has_spectrum(Acl, want, 1e-6)); // a pole at 0 twice, one input: a Jordan chain, its eigenvalues ε^{1/2} apart
        const auto gain = design::place(StateSpace{.A = model.A, .B = model.B, .C = Matrix<1, 4>{}, .Ts = kTs}, {z[0], z[1], Cplx{}, Cplx{}});
        CHECK((res.K - gain.K).norm() <= 1e-12 * gain.K.norm());

        //! An IntegralAction with the delay: the list is the plant's and the integrator's.
        const auto lqi = design::state_feedback(kCart, kTs, design::Poles{-8.0, -12.0, -4.0}, design::IntegralAction{.Q = 1.0}, InputDelay{0.004});
        REQUIRE(lqi.success);
        const auto delayed = *discretize(kCart, kTs, InputDelay{0.004});
        const auto im = design::detail::integral_model(delayed.A, delayed.B, delayed.C, delayed.D, Matrix{{1.0}}, kTs);
        const auto zi = to_z(std::array<Cplx, 3>{Cplx{-8.0, 0.0}, Cplx{-12.0, 0.0}, Cplx{-4.0, 0.0}}, kTs);
        CHECK(characteristic_residual((im.A - (im.B * lqi.K)), {zi[0], zi[1], zi[2], Cplx{}}) <= 1e-12);
    }

    TEST_CASE("ObserverPoles place the filter's error dynamics (I − MC)A, without a disturbance spec and with one") {
        //! LQG-shaped: the model is the plant's, M design::place_observer's.
        const auto lqg = design::state_feedback(kCart, kTs, kWeights, design::ObserverPoles{-30.0, -35.0});
        REQUIRE(lqg.success);
        const auto lqg_kalman = design::state_feedback(kCart, kTs, kWeights, design::KalmanNoise{.process = Matrix<2, 2>::identity(), .measurement = 0.01});
        static_assert(same_type<decltype(lqg), decltype(lqg_kalman)>); // the same filter type, so the same runtime
        const auto z = to_z(std::array<Cplx, 2>{Cplx{-30.0, 0.0}, Cplx{-35.0, 0.0}}, kTs);
        CHECK(has_spectrum(error_dynamics(lqg.kalman), z, 1e-12));
        const auto observer = design::place_observer(kCart, {-30.0, -35.0}, kTs);
        REQUIRE(observer.success);
        CHECK((lqg.kalman.M - observer.M).norm() <= 1e-12 * observer.M.norm());
        CHECK((lqg.kalman.L - (lqg.kalman.sys.A * lqg.kalman.M)).norm() <= 1e-15 * lqg.kalman.L.norm() + 1e-15);
        CHECK((lqg.kalman.sys.A - lqg_kalman.kalman.sys.A).norm() == 0.0); // the Kalman design's model
        CHECK((lqg.kalman.sys.G - lqg_kalman.kalman.sys.G).norm() == 0.0);

        //! OFLQG-shaped: the poles are the plant's and the disturbance's.
        static_assert(kPlacedObserver.success);
        const auto oflqg_kalman = design::state_feedback(
            kCart, kTs, design::Poles{-8.0, -12.0}, design::KalmanNoise{.process = Matrix<2, 2>::identity(), .measurement = 0.01, .disturbance = 1.0}, design::InputDisturbance{}
        );
        static_assert(same_type<decltype(kPlacedObserver), decltype(oflqg_kalman)>);
        CHECK((kPlacedObserver.kalman.sys.A - oflqg_kalman.kalman.sys.A).norm() == 0.0);
        CHECK((kPlacedObserver.kalman.sys.G - oflqg_kalman.kalman.sys.G).norm() == 0.0);
        const auto zd = to_z(std::array<Cplx, 3>{Cplx{-30.0, 0.0}, Cplx{-35.0, 0.0}, Cplx{-40.0, 0.0}}, kTs);
        CHECK(has_spectrum(error_dynamics(kPlacedObserver.kalman), zd, 1e-10));                  // M reaches 250: the placement's rounding
        CHECK((kPlacedObserver.target.target_d - oflqg_kalman.target.target_d).norm() <= 1e-12); // the target is the model's

        //! A bias on the reading of a discrete plant: z-plane poles. (On the double integrator a bias and
        //! the position are one mode at z = 1, which no filter tells apart.)
        constexpr StateSpace lag{.A = Matrix{{0.9, 0.05}, {0.0, 0.8}}, .B = Matrix{{0.0}, {0.1}}, .C = Matrix{{1.0, 0.0}}, .Ts = kTs};
        const auto           biased = design::state_feedback(lag, kWeights, design::ObserverPoles{0.6, 0.65, 0.7}, design::OutputDisturbance{});
        REQUIRE(biased.success);
        CHECK(has_spectrum(error_dynamics(biased.kalman), {Cplx{0.6, 0.0}, Cplx{0.65, 0.0}, Cplx{0.7, 0.0}}, 1e-12));
    }

    TEST_CASE("ObserverPoles with delays: the delay states' error poles at 0, the input's not corrected") {
        const auto zp = to_z(std::array<Cplx, 2>{Cplx{-30.0, 0.0}, Cplx{-35.0, 0.0}}, kTs);
        //! An InputDelay: its delay states hold the inputs sent, and M has no rows on them.
        const auto in = design::state_feedback(kCart, kTs, design::Poles{-8.0, -12.0}, design::ObserverPoles{-30.0, -35.0}, InputDelay<2>{0.015});
        REQUIRE(in.success);
        static_assert(same_type<decltype(in), decltype(design::state_feedback(kCart, kTs, kWeights, design::KalmanNoise{.process = Matrix<2, 2>::identity(), .measurement = 0.01}, InputDelay<2>{0.015}))>);
        CHECK(characteristic_residual(error_dynamics(in.kalman), {zp[0], zp[1], Cplx{}, Cplx{}}) <= 1e-12);
        CHECK(in.kalman.M(2, 0) == 0.0);
        CHECK(in.kalman.M(3, 0) == 0.0);

        //! With an InputDisturbance d rides through the delay states: they take d's correction.
        const auto zd = to_z(std::array<Cplx, 3>{Cplx{-30.0, 0.0}, Cplx{-35.0, 0.0}, Cplx{-40.0, 0.0}}, kTs);
        const auto of = design::state_feedback(kCart, kTs, design::Poles{-8.0, -12.0}, design::ObserverPoles{-30.0, -35.0, -40.0}, design::InputDisturbance{}, InputDelay<2>{0.015});
        REQUIRE(of.success);
        CHECK(characteristic_residual(error_dynamics(of.kalman), {zd[0], zd[1], zd[2], Cplx{}, Cplx{}}) <= 1e-12);
        CHECK(of.kalman.M(2, 0) == doctest::Approx(of.kalman.M(4, 0)).epsilon(1e-12));
        CHECK(of.kalman.M(3, 0) == doctest::Approx(of.kalman.M(4, 0)).epsilon(1e-12));

        //! Two readings, three delay states per input: more poles at 0 than readings, which the core leaves out.
        const auto two = design::state_feedback(kTwoInput, kTs, design::LQRWeights<3, 2>{.Q = Matrix<3, 3>::identity(), .R = Matrix<2, 2>::identity()}, design::ObserverPoles{-20.0, -25.0, -30.0}, InputDelay<3>{0.025});
        REQUIRE(two.success);
        const auto          zt = to_z(std::array<Cplx, 3>{Cplx{-20.0, 0.0}, Cplx{-25.0, 0.0}, Cplx{-30.0, 0.0}}, kTs);
        std::array<Cplx, 9> want{};
        for (size_t i = 0; i < 3; ++i) {
            want[i] = zt[i];
        }
        CHECK(characteristic_residual(error_dynamics(two.kalman), want) <= 1e-12);

        //! A late reading: its delay state is placed at 0 with the rest.
        const auto late = design::state_feedback(kCart, kTs, kWeights, design::ObserverPoles{-30.0, -35.0}, OutputDelay{0.005});
        REQUIRE(late.success);
        CHECK(characteristic_residual(error_dynamics(late.kalman), {zp[0], zp[1], Cplx{}}) <= 1e-12);
        CHECK(has_spectrum(error_dynamics(late.kalman), {zp[0], zp[1], Cplx{}}, 1e-9));
    }

    TEST_CASE("the mixed combinations: LQRWeights with ObserverPoles, Poles with a KalmanNoise") {
        const design::KalmanNoise noise{.process = Matrix<2, 2>::identity(), .measurement = 0.01};
        const auto                lq_kalman = design::state_feedback(kCart, kTs, kWeights, noise);
        const auto                lq_placed = design::state_feedback(kCart, kTs, kWeights, design::ObserverPoles{-30.0, -35.0});
        const auto                placed_kalman = design::state_feedback(kCart, kTs, design::Poles{-8.0, -12.0}, noise);
        REQUIRE(lq_kalman.success);
        REQUIRE(lq_placed.success);
        REQUIRE(placed_kalman.success);
        CHECK((lq_placed.K - lq_kalman.K).norm() == 0.0);                   // the LQ gain either way
        CHECK((placed_kalman.kalman.M - lq_kalman.kalman.M).norm() == 0.0); // the Kalman gain either way
        CHECK((placed_kalman.K - design::place(kCart, {-8.0, -12.0}, kTs).K).norm() <= 1e-12 * placed_kalman.K.norm());
        static_assert(same_type<decltype(lq_placed), decltype(placed_kalman)>);

        //! An LQGI: the poles cover the integrator, the filter is placed.
        const auto lqgi = design::state_feedback(kCart, kTs, design::Poles{-8.0, -12.0, -4.0}, design::ObserverPoles{-30.0, -35.0}, design::IntegralAction{.Q = 1.0});
        REQUIRE(lqgi.success);
        static_assert(std::remove_cvref_t<decltype(lqgi)>::kIntegral && std::remove_cvref_t<decltype(lqgi)>::kEstimated);

        //! A Handover, like every other option: the LQ design's type with it.
        const auto handed = design::state_feedback(kCart, kTs, design::Poles{-8.0, -12.0}, design::ObserverPoles{-30.0, -35.0}, design::Handover{.fade_time = 0.1});
        REQUIRE(handed.success);
        static_assert(same_type<decltype(handed), decltype(design::state_feedback(kCart, kTs, kWeights, noise, design::Handover{.fade_time = 0.1}))>);
        CHECK((handed.K - placed_kalman.K).norm() == 0.0);
    }

    TEST_CASE("ADRC is an OFLQG of the chain of integrators: its ESO is ObserverPoles at −wo, its PD gain not ZOH-placed") {
        //! ADRC's observer is the ZOH-discretized ESO in current form, (I − L·C)·Φ, every pole at
        //! e^{−wo·Ts}; ObserverPoles at −wo on the chain with an InputDisturbance (d in u's units, f = b0·d)
        //! are the same filter. Its PD gains are the continuous design's (wc, or wc² and 2·wc), so its
        //! discrete loop has eig(Φ − Γ·K) at 1 − wc·Ts (NX = 1), not e^{−wc·Ts}: Poles at −wc are the
        //! ZOH-exact gain, and Poles at ADRC's discrete loop poles on the discrete chain are ADRC.
        constexpr double wc = 8.0;
        constexpr double wo = 40.0;
        constexpr double b0 = 2.5;
        constexpr double Ts = 1e-3;
        const double     lambda = std::exp(-wo * Ts);

        //! NX = 1
        {
            constexpr StateSpace chain{.A = Matrix{{0.0}}, .B = Matrix<1, 1>{{b0}}, .C = Matrix{{1.0}}};
            const auto           adrc_res = design::adrc<1>(wc, wo, b0, Ts);
            REQUIRE(adrc_res.success);
            const auto sampled = design::state_feedback(chain, Ts, design::Poles{-wc}, design::ObserverPoles{-wo, -wo}, design::InputDisturbance{});
            REQUIRE(sampled.success);
            static_assert(std::is_same_v<decltype(OFLQG{sampled}), OFLQG<1, 1, 1, 1, double, FixedGainObserver<2, 1, 1, double, 1>>>);
            //! The filter is the ESO: M on [y; d] is L on [y; f], d = f / b0.
            CHECK(sampled.kalman.M(0, 0) == doctest::Approx(adrc_res.L(0)).epsilon(1e-12));
            CHECK(sampled.kalman.M(1, 0) * b0 == doctest::Approx(adrc_res.L(1)).epsilon(1e-12));
            //! The gain is the ZOH-exact one, (1 − e^{−wc·Ts})/(b0·Ts); ADRC's is wc/b0.
            CHECK(sampled.K(0, 0) == doctest::Approx((1.0 - std::exp(-wc * Ts)) / (b0 * Ts)).epsilon(1e-12));
            CHECK(adrc_res.Kp / b0 == doctest::Approx(wc / b0).epsilon(1e-15));

            //! Poles at ADRC's discrete loop pole 1 − wc·Ts on the discrete chain: ADRC's commands.
            const auto chain_d = *discretize(chain, Ts);
            const auto equal = design::state_feedback(chain_d, design::Poles{1.0 - (wc * Ts)}, design::ObserverPoles{lambda, lambda}, design::InputDisturbance{});
            REQUIRE(equal.success);
            CHECK(equal.K(0, 0) == doctest::Approx(adrc_res.Kp / b0).epsilon(1e-12));
            ADRCController<1, double> adrc{adrc_res};
            OFLQG                     law{equal};
            OFLQG                     zoh{sampled};
            double                    x = 0.0;
            double                    worst = 0.0;
            double                    apart = 0.0;
            for (int k = 0; k < 3000; ++k) {
                const double y = x;
                const double r = (k < 1500) ? 1.0 : -0.5;
                const double u_adrc = adrc.step(r, y);
                const double u_law = law.step(ColVec<1>{r}, ColVec<1>{y})[0];
                const double u_zoh = zoh.step(ColVec<1>{r}, ColVec<1>{y})[0];
                worst = std::max(worst, std::abs(u_law - u_adrc) / (1.0 + std::abs(u_adrc)));
                apart = std::max(apart, std::abs(u_zoh - u_adrc) / (1.0 + std::abs(u_adrc)));
                x += b0 * Ts * (u_adrc + ((k > 800) ? 0.4 : 0.0)); // the plant, a load on its input
            }
            CHECK(worst <= 1e-9);
            CHECK(apart > 1e-4); // the ZOH-placed gain is a different law
        }

        //! NX = 2
        {
            constexpr StateSpace chain{.A = Matrix{{0.0, 1.0}, {0.0, 0.0}}, .B = Matrix<2, 1>{{0.0}, {b0}}, .C = Matrix{{1.0, 0.0}}};
            const auto           adrc_res = design::adrc<2>(wc, wo, b0, Ts);
            REQUIRE(adrc_res.success);
            const auto sampled = design::state_feedback(chain, Ts, design::Poles{-wc, -wc}, design::ObserverPoles{-wo, -wo, -wo}, design::InputDisturbance{});
            REQUIRE(sampled.success);
            for (size_t i = 0; i < 2; ++i) {
                CHECK(sampled.kalman.M(i, 0) == doctest::Approx(adrc_res.L(i)).epsilon(1e-9));
            }
            CHECK(sampled.kalman.M(2, 0) * b0 == doctest::Approx(adrc_res.L(2)).epsilon(1e-9));

            //! ADRC's discrete loop: eig(Φ − Γ·[wc², 2·wc]/b0) on the ZOH chain.
            const auto         chain_d = *discretize(chain, Ts);
            const Matrix<1, 2> K_adrc{{adrc_res.Kp / b0, adrc_res.Kd / b0}};
            const auto         loop = mat::compute_eigenvalues((chain_d.A - (chain_d.B * K_adrc)));
            REQUIRE(loop.success);
            const Cplx z0 = loop.values[0];
            const Cplx z1 = loop.values[1];
            const auto equal = design::state_feedback(chain_d, design::Poles{z0, z1}, design::ObserverPoles{lambda, lambda, lambda}, design::InputDisturbance{});
            REQUIRE(equal.success);
            CHECK((equal.K - K_adrc).norm() <= 1e-9 * K_adrc.norm());
            ADRCController<2, double> adrc{adrc_res};
            OFLQG                     law{equal};
            OFLQG                     zoh{sampled};
            ColVec<2>                 x{};
            double                    worst = 0.0;
            double                    apart = 0.0;
            for (int k = 0; k < 3000; ++k) {
                const double y = x[0];
                const double r = (k < 1500) ? 1.0 : -0.5;
                const double u_adrc = adrc.step(r, y);
                const double u_law = law.step(ColVec<1>{r}, ColVec<1>{y})[0];
                const double u_zoh = zoh.step(ColVec<1>{r}, ColVec<1>{y})[0];
                worst = std::max(worst, std::abs(u_law - u_adrc) / (1.0 + std::abs(u_adrc)));
                apart = std::max(apart, std::abs(u_zoh - u_adrc) / (1.0 + std::abs(u_adrc)));
                x = (chain_d.A * x) + (chain_d.B * ColVec<1>{u_adrc + ((k > 800) ? 0.4 : 0.0)});
            }
            CHECK(worst <= 1e-9);
            CHECK(apart > 1e-4);
        }
    }

    TEST_CASE("Law<kRes, float> of a constexpr placed design: the LQ runtime's commands, sizes unchanged") {
        static constinit Law<kPlacedDiscrete, float> law;
        LQR<2, 1, float>                             held{kPlacedDiscrete.as<float>()};
        const ColVec                                 x{0.3f, -0.1f};
        CHECK(law.step(x)[0] == doctest::Approx(held.step(x)[0]).epsilon(1e-6));
        static_assert(std::is_same_v<Law<kPlacedDiscrete, float>::Runtime, LQR<2, 1, float>>);

        static constinit Law<kPlacedObserver, float> of;
        Law<kPlacedObserver, float>::Runtime         of_held{kPlacedObserver.as<float>()};
        const ColVec                                 r{1.0f};
        const ColVec                                 y{0.2f};
        CHECK(of.step(r, y)[0] == doctest::Approx(of_held.step(r, y)[0]).epsilon(1e-5));

        //! The same types as the LQ designs of the same shape, so the same runtimes and sizes.
        constexpr auto lq = design::state_feedback(kCartD, kWeights);
        constexpr auto lq_of = design::state_feedback(kCart, kTs, kWeights, design::KalmanNoise{.process = Matrix<2, 2>::identity(), .measurement = 0.01, .disturbance = 1.0}, design::InputDisturbance{});
        static_assert(same_type<decltype(kPlacedDiscrete), decltype(lq)>);
        static_assert(same_type<decltype(kPlacedObserver), decltype(lq_of)>);
        static_assert(std::is_same_v<Law<kPlacedObserver, float>::Runtime, decltype(OFLQG{lq_of.as<float>()})>);
    }

    TEST_CASE("design::mpc takes ObserverPoles as state_feedback does: the same placed filter, the same runtime") {
        constexpr design::MPCWeights w{.Qy = {{10.0}}, .Rdu = {0.1}};
        //! Without a disturbance spec: the LQG-shaped filter of the same poles.
        const auto regulator = design::mpc<50, 10>(kCart, kTs, w, design::ObserverPoles{-30.0, -35.0}, design::DareCost{});
        REQUIRE(regulator.success);
        const auto lq_placed = design::state_feedback(kCart, kTs, kWeights, design::ObserverPoles{-30.0, -35.0});
        CHECK((regulator.kalman.M - lq_placed.kalman.M).norm() <= 1e-12 * lq_placed.kalman.M.norm());

        //! Offset-free: the filter on [x; d] at the plant's and the disturbance's poles, the OFLQG's. M
        //! reaches 250, so the placement's rounding is relative to its conditioning, not to eps.
        const auto of = design::mpc<50, 10>(kCart, kTs, w, design::ObserverPoles{-30.0, -35.0, -40.0}, design::InputDisturbance{}, design::DareCost{});
        REQUIRE(of.success);
        const auto oflqg = design::state_feedback(kCart, kTs, design::Poles{-8.0, -12.0}, design::ObserverPoles{-30.0, -35.0, -40.0}, design::InputDisturbance{});
        CHECK((of.kalman.M - oflqg.kalman.M).norm() <= 1e-10 * oflqg.kalman.M.norm());
        CHECK((of.kalman.L - oflqg.kalman.L).norm() <= 1e-10 * oflqg.kalman.L.norm());
        CHECK(has_spectrum(error_dynamics(of.kalman), to_z(std::array<Cplx, 3>{Cplx{-30.0, 0.0}, Cplx{-35.0, 0.0}, Cplx{-40.0, 0.0}}, kTs), 1e-10));
        //! The Kalman design's type, so the same runtime.
        const auto of_kalman = design::mpc<50, 10>(
            kCart, kTs, w, design::KalmanNoise{.process = Matrix<2, 2>::identity(), .measurement = 0.01, .disturbance = 1.0}, design::InputDisturbance{}, design::DareCost{}
        );
        static_assert(same_type<decltype(of), decltype(of_kalman)>);

        //! In closed loop it holds r under a load at the input it is not told of, and estimates the load.
        OFMPC            law{of};
        ColVec           x{0.0, 0.0};
        constexpr double load = 0.5;
        for (size_t k = 0; k < 3000; ++k) { // 30 s: the move weight makes the loop slow
            const double u = law.step(ColVec{1.0}, ColVec{x(0)})(0);
            x = (kCartD.A * x) + (kCartD.B * ColVec{u + load});
        }
        CHECK(x(0) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(law.estimator.disturbance_estimate()(0) == doctest::Approx(load).epsilon(1e-6));

        //! At compile time (each one a static_assert):
        //!   design::mpc<50, 10>(kCart, kTs, w, design::ObserverPoles{-30.0, -35.0}, design::InputDisturbance{}); // observer_pole_count<3, 2>
        //!   design::mpc<50, 10>(kCart, kTs, w, design::ObserverPoles{-30.0, -35.0}, design::KalmanNoise{...});  // one filter
    }

    TEST_CASE("refusals: an unreachable mode, a repeat placement does not take, a pole that is not finite") {
        //! A mode the input does not reach, and one the reading does not see: not_placeable.
        constexpr StateSpace split{.A = Matrix{{-1.0, 0.0}, {0.0, -2.0}}, .B = Matrix{{1.0}, {0.0}}, .C = Matrix{{1.0, 0.0}}};
        CHECK(design::state_feedback(split, kTs, design::Poles{-3.0, -4.0}).status == design::DesignStatus::not_placeable);
        CHECK(design::state_feedback(kCart, kTs, kWeights, design::ObserverPoles{-30.0, -35.0}).success);
        CHECK(design::state_feedback(split, kTs, design::Poles{-3.0, -2.0}, design::ObserverPoles{-30.0, -35.0}).status == design::DesignStatus::not_placeable);

        //! A pole repeated more often than there are inputs: KNV places independent eigenvectors.
        CHECK(design::state_feedback(kTwoInput, kTs, design::Poles{-2.0, -2.0, -2.0}).status == design::DesignStatus::not_placeable);
        //! So is a deadbeat chain of two delay states per input on two inputs.
        CHECK(design::state_feedback(kTwoInput, kTs, design::Poles{-2.0, -3.0, -4.0}, InputDelay<2>{0.015}).status == design::DesignStatus::not_placeable);
        CHECK(design::state_feedback(kTwoInput, kTs, design::Poles{-2.0, -3.0, -4.0}, InputDelay{0.004}).success); // one each: NU poles at 0

        //! Not finite: NaN, and −∞ in the s-plane, which z = e^{s·Ts} would map to a finite 0.
        for (const double bad : {runtime_nan(), runtime_minus_inf()}) {
            CHECK(design::state_feedback(kCart, kTs, design::Poles{bad, -12.0}).status == design::DesignStatus::invalid_parameter);
            CHECK(design::state_feedback(kCartD, design::Poles{bad, 0.9}).status == design::DesignStatus::invalid_parameter);
            CHECK(design::state_feedback(kCart, kTs, kWeights, design::ObserverPoles{-30.0, bad}).status == design::DesignStatus::invalid_parameter);
            CHECK(design::state_feedback(kCart, kTs, design::Poles{-8.0, -12.0}, design::ObserverPoles{bad, -35.0, -40.0}, design::InputDisturbance{}).status == design::DesignStatus::invalid_parameter);
        }
        //! At compile time (each one a static_assert naming its counts or the option given twice):
        //!   design::state_feedback(kCart, kTs, design::Poles{-8.0});                                       // controller_pole_count<2, 1>
        //!   design::state_feedback(kCart, kTs, design::Poles{-8.0, -12.0}, design::IntegralAction{.Q = 1.0}); // controller_pole_count<3, 2>
        //!   design::state_feedback(kCart, kTs, kWeights, design::ObserverPoles{-30.0, -35.0}, design::InputDisturbance{}); // observer_pole_count<3, 2>
        //!   design::state_feedback(kCart, kTs, kWeights, design::ObserverPoles{-30.0, -35.0}, design::KalmanNoise{...}); // an option given twice
        //!   design::state_feedback(kCart, kTs, kWeights, design::Poles{-8.0, -12.0});                       // the weights given twice
    }
}
