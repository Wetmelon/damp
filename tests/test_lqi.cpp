// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/status.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {
// Discrete double integrator at Ts = 0.1 s, position measured. NX=2, NU=1, NY=1.
constexpr StateSpace<2, 1, 1, double, 1, 1> make_plant() {
    return StateSpace<2, 1, 1, double, 1, 1>{
        .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
        .B = Matrix<2, 1>{{0.005}, {0.1}},
        .C = Matrix<1, 2>{{1.0, 0.0}}, // measure position
        .Ts = 0.1,
    };
}

constexpr auto design_lqi() {
    return design::lqi(make_plant(), {.Q = Matrix<2, 2>{{10.0, 0.0}, {0.0, 1.0}}, .R = 1.0, .Q_integral = 500.0});
}
} // namespace

TEST_SUITE("LQI") {
    TEST_CASE("LQIResult::as<U>() conversion") {
        constexpr auto lqi_d = design::lqi(StateSpace<1, 1, 1, double, 1, 1>{.A = Matrix<1, 1>{{1.0}}, .B = Matrix<1, 1>{{1.0}}, .C = Matrix<1, 1>{{1.0}}, .Ts = 1.0}, {.Q = 1.0, .R = 1.0, .Q_integral = 1.0});

        constexpr auto lqi_f = lqi_d.as<float>();

        static_assert(lqi_f.success);
        CHECK(lqi_f.success);
        CHECK(lqi_f.K(0, 0) != 0.0f);
        CHECK(lqi_f.K(0, 1) != 0.0f);
    }

    TEST_CASE("own integrator matches manual augmented state") {
        constexpr auto r = design_lqi();
        auto           a = LQI{r}; // uses internal integrator
        auto           b = LQI{r}; // fed a manual augmented vector

        ColVec<1>       xi{{0.0}};
        const ColVec<2> states[] = {{{0.5}, {0.0}}, {{0.4}, {-0.1}}, {{0.2}, {0.05}}};
        const double    refs[] = {1.0, 1.0, 1.0};
        for (int k = 0; k < 3; ++k) {
            const ColVec<1> y{{states[k][0]}};
            const ColVec<1> rf{{refs[k]}};

            const ColVec<3> x_aug{{states[k][0]}, {states[k][1]}, {xi[0]}};
            const ColVec<1> ub = -(b.K * x_aug);        // manual
            const ColVec<1> ua = a.step(rf, states[k]); // owns xi: integrates last tick's error, then the law
            CHECK(ua[0] == doctest::Approx(ub[0]).epsilon(1e-12));
            CHECK(a.integral.xi[0] == doctest::Approx(xi[0]).epsilon(1e-12));

            xi[0] += make_plant().Ts * (refs[k] - y[0]); // mirror the controller's update (xi += Ts·e), which the next step applies
        }
    }

    TEST_CASE("closed loop drives output to the reference with zero steady-state error") {
        const auto sys = make_plant();
        const auto r = design_lqi();
        REQUIRE(r.success);

        auto            ctrl = LQI{r};
        ColVec<2>       x{{0.0}, {0.0}};
        const ColVec<1> ref{{1.0}};
        for (int k = 0; k < 500; ++k) {
            const ColVec<1> y{{x[0]}};
            const ColVec<1> u = ctrl.step(ref, x);
            x = sys.A * x + sys.B * u;
        }
        CHECK(x[0] == doctest::Approx(1.0).epsilon(1e-3)); // tracks position
        CHECK(x[1] == doctest::Approx(0.0).epsilon(1e-3)); // at rest

        ctrl.reset();
        CHECK(ctrl.integral.xi[0] == doctest::Approx(0.0));
    }

    TEST_CASE("LQIResult::to_ss realizes the runtime, and closed around the plant has the design's poles") {
        const auto sys = make_plant();
        const auto res = design_lqi();
        REQUIRE(res.success);

        const auto ss = res.to_ss(); // state xi, in [r; x], out u
        CHECK(ss.Ts == doctest::Approx(sys.Ts));
        auto      runtime = LQI{res};
        ColVec<1> xc{};
        ColVec<2> xp{{0.3}, {-0.1}};
        for (int k = 0; k < 8; ++k) {
            const ColVec<1> r{{1.0}};
            const ColVec<3> in{r[0], xp[0], xp[1]};
            const ColVec<1> u = ss.C * xc + ss.D * in;
            CHECK(u[0] == doctest::Approx(runtime.step(r, xp)[0]).epsilon(1e-12));
            xc = ss.A * xc + ss.B * in;
            xp = sys.A * xp + sys.B * u;
        }

        // The plant reports [0; x] so the law's inputs [r; x] close by unity positive feedback.
        const StateSpace<2, 1, 3> plant_out{
            .A = sys.A,
            .B = sys.B,
            .C = Matrix<3, 2>{{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}},
            .Ts = sys.Ts,
        };
        StateSpace<0, 3, 3> positive{.Ts = sys.Ts};
        positive.D = -Matrix<3, 3>::identity();
        const auto loop = series(ss, plant_out);
        REQUIRE(loop.has_value());
        const auto cl = feedback(*loop, positive);
        REQUIRE(cl.has_value());

        const auto poles = mat::compute_eigenvalues(cl->A).values; // [xi; x]
        for (size_t i = 0; i < 3; ++i) {
            double nearest = std::numeric_limits<double>::max();
            for (size_t j = 0; j < 3; ++j) {
                nearest = damp::min(nearest, (poles[j] - res.e[i]).abs());
            }
            CHECK(nearest < 1e-9);
        }

        // Integral action: unit DC gain from r to y.
        const auto x_dc = mat::solve(Matrix<3, 3>::identity() - cl->A, cl->B);
        REQUIRE(x_dc.has_value());
        const auto dc = cl->C * (*x_dc) + cl->D;
        CHECK(dc(1, 0) == doctest::Approx(1.0).epsilon(1e-9));
    }

    TEST_CASE("integrating a subset of outputs matches a hand-built single-output model") {
        // Same double integrator, two sensors (position, velocity); integrate position only.
        constexpr StateSpace<2, 1, 2, double> sys2{
            .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
            .B = Matrix<2, 1>{{0.005}, {0.1}},
            .C = Matrix<2, 2>::identity(),
            .Ts = 0.1
        };
        constexpr design::LQRConfig<2, 1, 1> w{.Q = Matrix<2, 2>{{10.0, 0.0}, {0.0, 1.0}}, .R = 1.0, .Q_integral = 500.0}; // design_lqi's
        constexpr auto                       picked = design::lqi(sys2, w, design::IntegralAction<1>{{0}});
        constexpr auto                       by_hand = design_lqi(); // C = [1 0]
        static_assert(picked.success);
        for (size_t j = 0; j < 3; ++j) {
            CHECK(picked.K(0, j) == doctest::Approx(by_hand.K(0, j)).epsilon(1e-12));
        }
        CHECK(picked.integral.Ci(0, 0) == 1.0);
        CHECK(picked.integral.Ci(0, 1) == 0.0);

        // The selector form takes the reference for the integrated output and the measured state (y = x here, C = I).
        auto            ctrl = LQI{picked.as<float>()};
        ColVec<2>       x{};
        const ColVec<1> r{{1.0}};
        for (int k = 0; k < 500; ++k) {
            const ColVec<2, float> y{static_cast<float>(x[0]), static_cast<float>(x[1])};
            const ColVec<1, float> u = ctrl.step(r.as<float>(), y);
            x = sys2.A * x + sys2.B * u.as<double>();
        }
        CHECK(x[0] == doctest::Approx(1.0).epsilon(1e-3));

        // An index past the last output is an inconsistent design input.
        constexpr auto bad = design::lqi(sys2, w, design::IntegralAction<1>{{2}});
        static_assert(!bad.success && bad.status == design::DesignStatus::bad_weights);
        // Q_integral left out of the braces: the integrator is a mode at z = 1 the cost does not
        // see (the DARE would give Kᵢ = 0 and a closed-loop pole on the unit circle).
        constexpr auto unweighted = design::lqi(sys2, {.Q = w.Q, .R = w.R});
        static_assert(!unweighted.success && unweighted.status == design::DesignStatus::bad_weights);
        // Integrating both outputs of a single-input plant cannot be held apart: not stabilizable.
        const auto both = design::lqi(sys2, {.Q = Matrix<2, 2>::identity(), .R = 1.0, .Q_integral = Matrix<2, 2>::identity()}, design::IntegralAction<2>{{0, 1}});
        CHECK_FALSE(both.success);
        auto idle = LQI{both};
        idle.update(ColVec<2>{});
        CHECK(idle.control(ColVec<2>{{1.0}, {1.0}})[0] == 0.0); // failed design: zero output
    }

    TEST_CASE("control / predict hold the integrator while the actuator is clamped") {
        const auto sys = make_plant();
        const auto r = design_lqi();
        REQUIRE(r.success);
        const double u_max = 0.2;
        const auto   clamp_u = [&](const ColVec<1>& u) { return ColVec<1>{{damp::clamp(u[0], -u_max, u_max)}}; };

        const auto run = [&](bool aware) {
            auto      ctrl = LQI{r};
            ColVec<2> x{};
            double    peak = 0.0;
            double    xi_peak = 0.0;
            for (int k = 0; k < 2000; ++k) {
                const ColVec<1> y{{x[0]}};
                ctrl.update(x);
                const ColVec<1> u_req = ctrl.control(ColVec<1>{{1.0}});
                const ColVec<1> u = clamp_u(u_req);
                ctrl.back_calculate(aware ? u : u_req);
                ctrl.predict(); // the unaware loop reports the unclamped command
                x = sys.A * x + sys.B * u;
                peak = damp::max(peak, x[0]);
                xi_peak = damp::max(xi_peak, damp::abs(ctrl.integral.xi[0]));
            }
            return damp::array<double, 3>{peak, xi_peak, x[0]};
        };
        const auto naive = run(false);
        const auto aware = run(true);
        CHECK(aware[0] < naive[0]);                            // less overshoot
        CHECK(aware[1] < naive[1]);                            // the integral stays smaller
        CHECK(aware[2] == doctest::Approx(1.0).epsilon(1e-3)); // still offset-free

        // Feedforward before the clamp is not mistaken for saturation.
        auto            a = LQI{r};
        auto            b = LQI{r};
        const ColVec<1> y{{0.1}};
        const ColVec<2> x{{0.1}, {0.0}};
        const ColVec<1> u_ff{{0.05}};
        a.update(x);
        const ColVec<1> ua = a.control(ColVec<1>{{0.2}});
        a.back_calculate(ua);
        a.predict();
        b.update(x);
        const ColVec<1> ub = b.control(ColVec<1>{{0.2}});
        b.back_calculate(ub + u_ff, u_ff);
        b.predict();
        CHECK(a.integral.xi[0] != 0.0);
        CHECK(b.integral.xi[0] == doctest::Approx(a.integral.xi[0]));
    }

    TEST_CASE("asymmetric clamp in both directions, then recovery to r") {
        const auto sys = make_plant();
        const auto r = design_lqi();
        REQUIRE(r.success);
        auto       ctrl = LQI{r};
        ColVec<2>  x{};
        const auto clamp_u = [](double u) { return damp::clamp(u, -0.15, 0.4); };
        bool       hit_hi = false;
        bool       hit_lo = false;
        for (int k = 0; k < 400; ++k) {
            ctrl.update(x);
            const double u = clamp_u(ctrl.control(ColVec<1>{{3.0}})[0]);
            hit_hi = hit_hi || (u == 0.4);
            ctrl.back_calculate(ColVec<1>{{u}});
            ctrl.predict();
            x = sys.A * x + sys.B * ColVec<1>{{u}};
        }
        CHECK(hit_hi);
        const double xi_hi = damp::abs(ctrl.integral.xi[0]);
        for (int k = 0; k < 400; ++k) {
            ctrl.update(x);
            const double u = clamp_u(ctrl.control(ColVec<1>{{-3.0}})[0]);
            hit_lo = hit_lo || (u == -0.15);
            ctrl.back_calculate(ColVec<1>{{u}});
            ctrl.predict();
            x = sys.A * x + sys.B * ColVec<1>{{u}};
        }
        CHECK(hit_lo);
        for (int k = 0; k < 2500; ++k) {
            const double u = ctrl.step(ColVec<1>{{1.0}}, x)[0];
            x = sys.A * x + sys.B * ColVec<1>{{u}};
        }
        CHECK(x[0] == doctest::Approx(1.0).epsilon(1e-3));
        CHECK(damp::abs(ctrl.integral.xi[0]) < xi_hi);
    }

    TEST_CASE("LQI deduces sizes and scalar from its design") {
        static constexpr auto res = design_lqi();
        static constinit auto law = LQI{res.as<float>()}; // sizes and scalar deduced from the result
        static_assert(std::is_same_v<decltype(law), LQI<2, 1, 1, float, 1>>);
        // A double design does not narrow to the float default.
        static_assert(std::is_same_v<decltype(LQI{res}), LQI<2, 1, 1, double, 1>>);
        CHECK(law.K(0, 0) == doctest::Approx(static_cast<float>(res.K(0, 0))));
    }
}

TEST_CASE("LQI weights: the cross term N reaches the augmented DARE") {
    const auto         sys = make_plant();
    const Matrix<2, 2> Qx{{10.0, 0.0}, {0.0, 1.0}};
    const Matrix<1, 1> Qi{{5.0}};
    const Matrix<1, 1> R{{1.0}};
    const Matrix<2, 1> N{{0.1}, {0.2}};
    Matrix<3, 3>       Q{}; // on [x; xi]: Q and Q_integral on the diagonal
    Q.block<2, 2>(0, 0) = Qx;
    Q(2, 2) = Qi(0, 0);
    Matrix<3, 1> N_aug{};
    N_aug.block<2, 1>(0, 0) = N;
    Matrix<3, 3> A_aug{};
    A_aug.block<2, 2>(0, 0) = sys.A;
    A_aug.block<1, 2>(2, 0) = -sys.Ts * sys.C; // xi += Ts·e
    A_aug(2, 2) = 1.0;
    Matrix<3, 1> B_aug{};
    B_aug.block<2, 1>(0, 0) = sys.B;

    const auto by_hand = design::discrete_lqr(A_aug, B_aug, Q, R, N_aug);
    const auto got = design::lqi(sys, {.Q = Qx, .R = R, .N = N, .Q_integral = Qi});
    REQUIRE(by_hand.success);
    REQUIRE(got.success);
    CHECK((got.K - by_hand.K).norm() == 0.0);
    CHECK((got.K - design::lqi(sys, {.Q = Qx, .R = R, .Q_integral = Qi}).K).norm() > 0.0); // N is not dropped
}

TEST_SUITE("LQI status") {
    TEST_CASE("a continuous plant is refused; Ts and status ride on the result") {
        auto continuous = make_plant();
        continuous.Ts = 0.0;
        const auto refused = design::lqi(continuous, {.Q = Matrix<2, 2>::identity(), .R = 1.0, .Q_integral = 1.0});
        CHECK_FALSE(refused.success);
        CHECK(refused.status == design::DesignStatus::invalid_sample_time);

        constexpr auto ok = design_lqi();
        static_assert(ok.success && ok.status == design::DesignStatus::ok && ok.Ts == 0.1);
        CHECK(ok.as<float>().status == design::DesignStatus::ok);
        CHECK(ok.as<float>().Ts == 0.1f);

        constexpr auto bad_index = design::lqi(make_plant(), {.Q = Matrix<2, 2>::identity(), .R = 1.0, .Q_integral = 1.0}, design::IntegralAction<1>{{1}});
        static_assert(!bad_index.success && bad_index.status == design::DesignStatus::bad_weights);
    }
}

namespace {
// Continuous double integrator, position measured: the plant of the sampled-data LQI tests.
constexpr StateSpace<2, 1, 1> make_continuous_plant() {
    return StateSpace<2, 1, 1>{
        .A = Matrix<2, 2>{{0.0, 1.0}, {0.0, 0.0}},
        .B = Matrix<2, 1>{{0.0}, {1.0}},
        .C = Matrix<1, 2>{{1.0, 0.0}},
    };
}

// Continuous cost on [position; velocity], u and ∫e dt.
constexpr design::LQRConfig<2, 1, 1> continuous_weights() {
    return {.Q = Matrix<2, 2>::diagonal({10.0, 1.0}), .R = 1.0, .Q_integral = 100.0};
}

// Closed-loop poles mapped to the s-plane, s = log(z)/Ts.
template<size_t N>
damp::array<std::complex<double>, N> s_plane(const damp::array<damp::complex<double>, N>& z, double Ts) {
    damp::array<std::complex<double>, N> s{};
    for (size_t i = 0; i < N; ++i) {
        s[i] = std::log(std::complex<double>(z[i].real(), z[i].imag())) / Ts;
    }
    return s;
}

// Largest relative distance from a pole of a to the nearest pole of b.
template<size_t N>
double pole_mismatch(const damp::array<std::complex<double>, N>& a, const damp::array<std::complex<double>, N>& b) {
    double worst = 0.0;
    for (const auto& sa : a) {
        double nearest = std::numeric_limits<double>::infinity();
        for (const auto& sb : b) {
            nearest = std::min(nearest, std::abs(sa - sb) / std::abs(sa));
        }
        worst = std::max(worst, nearest);
    }
    return worst;
}
} // namespace

TEST_SUITE("LQI sampled-data (continuous plant, Ts)") {
    TEST_CASE("the discrete weights are the cost integral over one sample of the runtime's loop") {
        //! Between samples: the plant on the held u, xi ramping on the held error −y_k (r = 0).
        constexpr auto plant = make_continuous_plant();
        const auto     w = design::LQRCost<3, 1>{.Q = Matrix<3, 3>{{10.0, 0.5, 0.2}, {0.5, 1.0, 0.0}, {0.2, 0.0, 100.0}}, .R = 2.0, .N = Matrix<3, 1>{{0.3}, {0.1}, {0.4}}};
        const double   Ts = 0.05;
        const auto     cost = design::detail::sampled_integral_cost(plant.A, plant.B, plant.C, plant.D, Matrix<1, 1>{{1.0}}, w, Ts).value(); // with x–∫e cross terms

        const double p0 = 0.7, v0 = -1.3, xi0 = 0.4, u = 0.9;
        const auto   integrand = [&](double t) {
            const ColVec<3> z{p0 + v0 * t + 0.5 * u * t * t, v0 + u * t, xi0 - p0 * t};
            return (z.t() * w.Q * z)(0, 0) + 2.0 * (z.t() * w.N)(0, 0) * u + w.R(0, 0) * u * u;
        };
        const int n = 2000; // Simpson: the integrand is a quartic in t
        double    sum = integrand(0.0) + integrand(Ts);
        for (int k = 1; k < n; ++k) {
            sum += ((k % 2) != 0 ? 4.0 : 2.0) * integrand(Ts * k / n);
        }
        const double    exact = sum * Ts / (3.0 * n);
        const ColVec<3> z0{p0, v0, xi0};
        const double    discrete = (z0.t() * cost.Q * z0)(0, 0) + 2.0 * (z0.t() * cost.N)(0, 0) * u + cost.R(0, 0) * u * u;
        CHECK(discrete == doctest::Approx(exact).epsilon(1e-12));
    }

    TEST_CASE("the same continuous weights give the same continuous closed loop at 1 kHz and 10 kHz") {
        constexpr auto plant = make_continuous_plant();
        const auto     fast = design::lqi(plant, 1e-4, continuous_weights());
        const auto     slow = design::lqi(plant, 1e-3, continuous_weights());
        REQUIRE(fast.success);
        REQUIRE(slow.success);
        CHECK(fast.integral.step == 1e-4);
        CHECK(slow.integral.step == 1e-3);
        CHECK(fast.Ts == 1e-4);

        const auto   s_fast = s_plane(fast.e, 1e-4);
        const auto   s_slow = s_plane(slow.e, 1e-3);
        const double mismatch = pole_mismatch(s_slow, s_fast);
        for (size_t i = 0; i < 3; ++i) {
            MESSAGE("s at 10 kHz: " << s_fast[i] << ", at 1 kHz: " << s_slow[i]);
        }
        MESSAGE("Ts form, 1 kHz vs 10 kHz: worst relative pole gap " << mismatch);
        CHECK(mismatch < 1e-5);
        const auto   coarse = design::lqi(plant, 1e-2, continuous_weights()); // 100 Hz, still ~30× the loop bandwidth
        const double coarse_mismatch = pole_mismatch(s_plane(coarse.e, 1e-2), s_fast);
        MESSAGE("Ts form, 100 Hz vs 10 kHz: worst relative pole gap " << coarse_mismatch);
        CHECK(coarse_mismatch < 1e-4); // the gap shrinks as Ts²

        //! The discrete-plant form integrates with xi += Ts·e too, so the same numbers give the same loop at each rate.
        const auto discrete_plant = [&](double Ts) {
            const auto zoh = discretize(plant, Ts, DiscretizationMethod::ZOH);
            return design::lqi(*zoh, continuous_weights());
        };
        const auto   d_fast = discrete_plant(1e-4);
        const auto   d_slow = discrete_plant(1e-3);
        const double d_mismatch = pole_mismatch(s_plane(d_slow.e, 1e-3), s_plane(d_fast.e, 1e-4));
        MESSAGE("discrete-plant form, 1 kHz vs 10 kHz: worst relative pole gap " << d_mismatch);
        CHECK(d_mismatch < 1e-5);
    }

    TEST_CASE("it deploys constinit in float and holds the output at r against an input disturbance") {
        static constexpr auto res = design::lqi(make_continuous_plant(), 1e-3, continuous_weights());
        static constinit auto law = LQI{res.as<float>()};
        static_assert(res.success && res.integral.step == 1e-3);
        CHECK(law.integral.step == 1e-3f);

        //! SIL: the plant itself, sampled exactly (ZOH), with a constant input disturbance.
        const auto      plant_d = *discretize(make_continuous_plant(), 1e-3, DiscretizationMethod::ZOH);
        const ColVec<1> r{1.0};
        ColVec<2>       x{};
        for (int k = 0; k < 20000; ++k) {
            const ColVec<1> y = plant_d.C * x;
            const ColVec<1> u = law.step(r.as<float>(), x.as<float>()).as<double>();
            x = plant_d.A * x + plant_d.B * (u + ColVec<1>{0.5});
        }
        CHECK(x[0] == doctest::Approx(1.0).epsilon(1e-4));
        CHECK(damp::abs(x[1]) < 1e-4);
        //! xi = ∫e dt, so the integral term carries the disturbance: Ki·xi = −0.5 − Kx·x.
        CHECK(static_cast<double>(law.integral.xi[0]) * res.K(0, 2) == doctest::Approx(0.5 - res.K(0, 0)).epsilon(1e-3));
    }

    TEST_CASE("Ts form: selectors, and the rows that refuse") {
        constexpr auto plant = make_continuous_plant();
        constexpr auto by_index = design::lqi(plant, 1e-3, continuous_weights(), design::IntegralAction<1>{{0}});
        constexpr auto every = design::lqi(plant, 1e-3, continuous_weights()); // Ci = I on the one output
        static_assert(by_index.success && every.success && by_index.K == every.K);
        static_assert(design::lqi(plant, 1e-3, continuous_weights(), design::IntegralAction<1>{{1}}).status == design::DesignStatus::bad_weights);
        CHECK(design::lqi(plant.as<float>(), 1e-3f, {.Q = continuous_weights().Q.as<float>(), .R = 1.0f, .Q_integral = continuous_weights().Q_integral.as<float>()}).success);

        //! The discrete-plant form integrates with forward Euler at the plant's Ts, as MATLAB®'s lqi does.
        CHECK(design_lqi().integral.step == make_plant().Ts);

        auto discrete = plant;
        discrete.Ts = 1e-3;
        CHECK(design::lqi(discrete, 1e-3, continuous_weights()).status == design::DesignStatus::invalid_sample_time);
        for (const double bad : {0.0, -1e-3, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
            const auto refused = design::lqi(plant, bad, continuous_weights());
            CHECK_FALSE(refused.success);
            CHECK(refused.status == design::DesignStatus::invalid_sample_time);
        }
    }
}

namespace {
// ẋ = −x + u, y = x + D·u: a first-order plant with feedthrough, continuous and sampled (ZOH) at 0.1 s.
constexpr StateSpace<1, 1, 1, double, 1, 1> feedthrough_plant(double D, double Ts) {
    const double a = (Ts > 0.0) ? damp::exp(-Ts) : -1.0;
    const double b = (Ts > 0.0) ? 1.0 - a : 1.0;
    return {.A = Matrix<1, 1>{{a}}, .B = Matrix<1, 1>{{b}}, .C = Matrix<1, 1>{{1.0}}, .D = Matrix<1, 1>{{D}}, .Ts = Ts};
}
} // namespace

TEST_SUITE("LQI on a plant with feedthrough") {
    TEST_CASE("the design and the runtime integrate the same output: the deployed loop has the design's poles and y holds r") {
        for (const double D : {0.5, -2.0}) {
            const auto plant = feedthrough_plant(D, 0.1);
            const auto w = design::LQRConfig<1, 1, 1>{.Q = 1.0, .R = 1.0, .Q_integral = 10.0};
            for (const bool sampled : {false, true}) {
                const auto res = sampled ? design::lqi(feedthrough_plant(D, 0.0), 0.1, w) : design::lqi(plant, w);
                REQUIRE(res.success);
                REQUIRE(res.is_stable());
                // Plant + to_ss at r = 0: [x; ξ]⁺ with u = C_c ξ + D_c·[r; x].
                const auto   c = res.to_ss();
                Matrix<2, 2> Acl{};
                Acl(0, 0) = plant.A(0, 0) + (plant.B(0, 0) * c.D(0, 1));
                Acl(0, 1) = plant.B(0, 0) * c.C(0, 0);
                Acl(1, 0) = c.B(0, 1);
                Acl(1, 1) = c.A(0, 0);
                const auto ev = mat::compute_eigenvalues(Acl);
                REQUIRE(ev.success);
                for (const auto& p : ev.values) {
                    const bool matched = (p - res.e[0]).abs() < 1e-9 || (p - res.e[1]).abs() < 1e-9;
                    CHECK(matched);
                }
                // The runtime around the plant: the output y = x + D·u, not just C·x, reaches r.
                auto   law = LQI{res};
                double x = 0.0;
                double y = 0.0;
                for (int k = 0; k < 2000; ++k) {
                    const double u = law.step(ColVec<1>{1.0}, ColVec<1>{x})[0];
                    y = x + (D * u);
                    x = (plant.A(0, 0) * x) + (plant.B(0, 0) * u);
                }
                CHECK(y == doctest::Approx(1.0).epsilon(1e-9));
            }
        }
    }

    TEST_CASE("a failed design is not reported stable") {
        // C = I on a double integrator, integrating the velocity: no stabilizing gain holds it.
        const StateSpace<2, 1, 2> plant{.A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}}, .B = Matrix<2, 1>{{0.005}, {0.1}}, .C = Matrix<2, 2>::identity(), .Ts = 0.1};
        const auto                res = design::lqi(plant, {.Q = Matrix<2, 2>::identity(), .R = 1.0, .Q_integral = 1.0}, design::IntegralAction<1>{{1}});
        CHECK_FALSE(res.success);
        CHECK_FALSE(res.is_stable());
    }
}
