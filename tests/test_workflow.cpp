// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>
#include <limits>
#include <numbers>

#include "damp/backend.hpp"
#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/pr.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/math/complex.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

// Unity positive feedback, v = r + w: closes a loop whose output w is the law's own input.
template<size_t N>
StateSpace<0, N, N> unity_positive(double Ts) {
    StateSpace<0, N, N> s{.Ts = Ts};
    s.D = -Matrix<N, N>::identity();
    return s;
}

// Every expected pole has its own closed-loop pole within tol.
template<size_t NC, size_t NE>
bool contains_poles(const ColVec<NC, complex<double>>& closed_loop, const ColVec<NE, complex<double>>& expected, double tol) {
    damp::array<bool, NC> used{};
    for (size_t i = 0; i < NE; ++i) {
        size_t best = NC;
        double nearest = std::numeric_limits<double>::max();
        for (size_t j = 0; j < NC; ++j) {
            const double d = (closed_loop[j] - expected[i]).abs();
            if (!used[j] && d < nearest) {
                nearest = d;
                best = j;
            }
        }
        if (best == NC || nearest > tol) {
            return false;
        }
        used[best] = true;
    }
    return true;
}

constexpr StateSpace<2, 1, 1, double, 2, 1> make_plant(double Ts) {
    return StateSpace<2, 1, 1, double, 2, 1>{
        .A = Matrix<2, 2>{{1.0, Ts}, {0.0, 1.0}},
        .B = Matrix<2, 1>{{0.5 * Ts * Ts}, {Ts}},
        .C = Matrix<1, 2>{{1.0, 0.0}},
        .D = Matrix<1, 1>::zeros(),
        .G = Matrix<2, 2>::identity(),
        .H = Matrix<1, 1>::identity(),
        .Ts = Ts,
    };
}

constexpr Matrix<2, 2> kQx = Matrix<2, 2>::identity(); // on [x]
constexpr Matrix<1, 1> kQInt{{10.0}};                  // on ∫e
constexpr Matrix<1, 1> kR{{0.1}};
constexpr Matrix<2, 2> kQKf{{0.01, 0.0}, {0.0, 0.01}};
constexpr Matrix<1, 1> kRKf{{0.1}};

} // namespace

TEST_SUITE("Workflow: design, result, runtime") {
    TEST_CASE("lqg: the runtime from the result ticks, and to_ss closed around the plant has the regulator and estimator poles") {
        const auto sys = make_plant(0.1);
        const auto res = design::lqg(sys, {.Q = Matrix<2, 2>::identity(), .R = kR}, design::KalmanNoise{.process = kQKf, .measurement = kRKf});
        REQUIRE(res.success);

        LQG runtime{res.as<float>()};
        CHECK(runtime.step(ColVec<1, float>{1.0f})(0, 0) != doctest::Approx(0.0f));

        // Separation: the loop's poles are those of A − BK and of the estimator error A(I − LC).
        const auto loop = series(res.to_ss(), sys); // y → u → y
        REQUIRE(loop.has_value());
        const auto cl = feedback(*loop, unity_positive<1>(sys.Ts));
        REQUIRE(cl.has_value());
        const auto                 observer = mat::compute_eigenvalues(sys.A * (Matrix<2, 2>::identity() - (res.kalman.M * sys.C))).values;
        ColVec<4, complex<double>> expected{};
        for (size_t i = 0; i < 2; ++i) {
            expected[i] = res.e[i];
            expected[2 + i] = observer[i];
        }
        CHECK(contains_poles(mat::compute_eigenvalues(cl->A).values, expected, 1e-9));
    }

    TEST_CASE("lqg + PR: both runtimes from their results, the PR on the tracking error") {
        const auto sys = make_plant(0.001);
        const auto lqg = design::lqg(sys, {.Q = Matrix<2, 2>::identity(), .R = kR}, design::KalmanNoise{.process = kQKf, .measurement = kRKf});
        const auto pr = design::pr({.Kp = 0.0, .Ki = 10.0, .w0 = 2.0 * std::numbers::pi, .wc = 5.0}, 0.001);
        REQUIRE(lqg.success);
        REQUIRE(pr.success);
        CHECK(pr.as<float>().Ki == doctest::Approx(10.0f));

        LQG                    regulator{lqg.as<float>()};
        PRController<float>    internal_model{pr.as<float>()};
        const float            r = 1.0f;
        const ColVec<1, float> y{0.0f};
        regulator.update(y);
        ColVec<1, float> u = regulator.control();
        u(0, 0) += internal_model.control(r - y(0, 0));
        regulator.back_calculate(u);
        regulator.predict();
        CHECK(u(0, 0) != doctest::Approx(0.0f));
    }

    TEST_CASE("lqi: the runtime from the result ticks") {
        StateSpace<2, 1, 1> sys{
            .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
            .B = Matrix<2, 1>{{0.005}, {0.1}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .D = Matrix<1, 1>::zeros(),
            .Ts = 0.1
        };
        const auto res = design::lqi(sys, {.Q = kQx, .R = kR, .Q_integral = kQInt});
        REQUIRE(res.success);

        auto                   runtime = LQI{res.as<float>()};
        const ColVec<2, float> x{1.0f, 0.0f};
        const ColVec<1, float> y{1.0f};
        const ColVec<1, float> r{0.0f};
        CHECK(runtime.step(r, x)(0, 0) != doctest::Approx(0.0f));
    }

    TEST_CASE("lqgi: the runtime from the result ticks, and to_ss closed around the plant has the servo and estimator poles") {
        const auto sys = make_plant(0.1);
        const auto res = design::lqgi(sys, {.Q = kQx, .R = kR, .Q_integral = kQInt}, design::KalmanNoise{.process = kQKf, .measurement = kRKf});
        REQUIRE(res.success);

        LQGI                   runtime{res.as<float>()};
        const ColVec<1, float> y{0.0f};
        const ColVec<1, float> r{1.0f};

        // u[k] = -K·[x̂[k]; xi[k]] uses the *current* integral state (the DARE augmentation
        // convention), so the first tick from a zero estimate and zero xi commands nothing;
        // xi picks up the error and the second tick drives the plant.
        CHECK(runtime.step(r, y)(0, 0) == doctest::Approx(0.0f));
        CHECK(runtime.step(r, y)(0, 0) != doctest::Approx(0.0f));

        // The plant reports [0; y] so the law's inputs [r; y] close by unity positive feedback.
        const StateSpace<2, 1, 2> plant_out{
            .A = sys.A,
            .B = sys.B,
            .C = Matrix<2, 2>{{0.0, 0.0}, {1.0, 0.0}},
            .Ts = sys.Ts,
        };
        const auto loop = series(res.to_ss(), plant_out);
        REQUIRE(loop.has_value());
        const auto cl = feedback(*loop, unity_positive<2>(sys.Ts));
        REQUIRE(cl.has_value());
        const auto                 observer = mat::compute_eigenvalues(sys.A * (Matrix<2, 2>::identity() - (res.kalman.M * sys.C))).values;
        ColVec<5, complex<double>> expected{};
        for (size_t i = 0; i < 3; ++i) {
            expected[i] = res.e[i];
        }
        for (size_t i = 0; i < 2; ++i) {
            expected[3 + i] = observer[i];
        }
        CHECK(contains_poles(mat::compute_eigenvalues(cl->A).values, expected, 1e-9));
    }

    TEST_CASE("lqi: the runtime drives the steady-state tracking error to zero") {
        StateSpace<2, 1, 1> sys{
            .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
            .B = Matrix<2, 1>{{0.005}, {0.1}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .D = Matrix<1, 1>::zeros(),
            .Ts = 0.1
        };
        const auto res = design::lqi(sys, {.Q = kQx, .R = kR, .Q_integral = kQInt});
        REQUIRE(res.success);

        // Closed-loop simulation against the (linear, discrete) design plant.
        // The runtime advances the integrator by xi += Ts·(r − y), the step the
        // design's cost was taken with (IntegralDesign::step), and y reaches r.
        const auto Af = sys.A.as<float>();
        const auto Bf = sys.B.as<float>();
        const auto Cf = sys.C.as<float>();

        auto             runtime = LQI{res.as<float>()};
        ColVec<2, float> x{0.0f, 0.0f};
        const float      r = 1.0f;
        float            y = 0.0f;

        for (int k = 0; k < 5000; ++k) {
            y = (Cf * x)(0, 0);
            const auto u = runtime.step(ColVec<1, float>{r}, x);
            x = Af * x + Bf * u;
        }

        CHECK(y == doctest::Approx(1.0f).epsilon(0.01));
    }

    TEST_CASE("lqi: design at compile time, constinit runtime") {
        static constexpr auto res = [] {
            StateSpace<2, 1, 1> sys{
                .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
                .B = Matrix<2, 1>{{0.005}, {0.1}},
                .C = Matrix<1, 2>{{1.0, 0.0}},
                .D = Matrix<1, 1>::zeros(),
                .Ts = 0.1
            };
            return design::lqi(sys, {.Q = kQx, .R = kR, .Q_integral = kQInt});
        }();
        static constinit LQI law{res.as<float>()}; // a failed design does not compile

        static_assert(res.success, "LQI must synthesize at compile time");
        CHECK(law.K(0, 0) == doctest::Approx(static_cast<float>(res.K(0, 0))));
    }

    TEST_CASE("lqi / lqgi integrate a subset of the outputs: the output list reaches to_ss and the runtime") {
        // Position and velocity sensors on a double integrator; integrate the position only.
        constexpr StateSpace<2, 1, 2, double, 2, 2> sys{
            .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
            .B = Matrix<2, 1>{{0.005}, {0.1}},
            .C = Matrix<2, 2>::identity(),
            .G = Matrix<2, 2>::identity(),
            .H = Matrix<2, 2>::identity(),
            .Ts = 0.1,
        };

        constexpr auto by_index = design::lqi(sys, {.Q = kQx, .R = kR, .Q_integral = kQInt}, design::IntegralAction<1>{{0}});
        static_assert(by_index.success);
        // The law integrates r − y₀ only: in [r; y; x], the integrator row reads [Ts, −Ts, 0, 0, 0] (xi += Ts·e).
        const auto servo = by_index.to_ss();
        CHECK(servo.B(0, 0) == doctest::Approx(sys.Ts));
        CHECK(servo.B(0, 1) == doctest::Approx(-sys.Ts));
        CHECK(servo.B(0, 2) == doctest::Approx(0.0));

        // The runtime takes one reference and both measurements, and holds y₀ at r.
        auto             lqi = LQI{by_index.as<float>()};
        ColVec<2, float> x{};
        for (int k = 0; k < 600; ++k) {
            const auto u = lqi.step(ColVec<1, float>{1.0F}, x);
            x = sys.A.as<float>() * x + sys.B.as<float>() * u;
        }
        CHECK(x(0) == doctest::Approx(1.0F).epsilon(1e-3));

        // LQGI: the same argument, an estimator underneath.
        const Matrix<2, 2> R_kf{{0.1, 0.0}, {0.0, 0.1}};
        const auto         lqgi_index = design::lqgi(sys, {.Q = kQx, .R = kR, .Q_integral = kQInt}, design::KalmanNoise{.process = kQKf, .measurement = R_kf}, design::IntegralAction<1>{{0}});
        REQUIRE(lqgi_index.success);
        CHECK(lqgi_index.K(0, 2) == doctest::Approx(by_index.K(0, 2)));
        // In [r; y], the integrator row (state NX) reads −Ts·y₀ only.
        const auto lqgi_ss = lqgi_index.to_ss();
        CHECK(lqgi_ss.B(2, 1) == doctest::Approx(-sys.Ts));
        CHECK(lqgi_ss.B(2, 2) == doctest::Approx(0.0));
        LQGI      lqgi{lqgi_index.as<float>()};
        ColVec<2> xp{};
        for (int k = 0; k < 600; ++k) {
            const auto u = lqgi.step(ColVec<1, float>{1.0F}, (sys.C * xp).as<float>());
            xp = sys.A * xp + sys.B * u.as<double>();
        }
        CHECK(xp(0) == doctest::Approx(1.0).epsilon(0.02));

        // An index past the last output is an inconsistent design input, in either design.
        const auto lqi_bad = design::lqi(sys, {.Q = kQx, .R = kR, .Q_integral = kQInt}, design::IntegralAction<1>{{2}});
        const auto lqgi_bad = design::lqgi(sys, {.Q = kQx, .R = kR, .Q_integral = kQInt}, design::KalmanNoise{.process = kQKf, .measurement = R_kf}, design::IntegralAction<1>{{2}});
        CHECK_FALSE(lqi_bad.success);
        CHECK(lqi_bad.status == design::DesignStatus::bad_weights);
        CHECK_FALSE(lqgi_bad.success);
        CHECK(lqgi_bad.status == design::DesignStatus::bad_weights);
    }
}
