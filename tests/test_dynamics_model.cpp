// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>
#include <random>

#include "damp/backend.hpp"
#include "damp/concepts.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/design/linearization.hpp"
#include "damp/estimation/ekf.hpp"
#include "damp/estimation/ukf.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "damp/simulation/integrator.hpp"
#include "damp/simulation/simulate.hpp"
#include "damp/simulation/solver.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

// One nonlinear model, written once as a Dynamics ẋ = f(x, u) and a
// Measurement y = h(x, u), drives everything that takes a nonlinear model: the
// SIL plant (sim::simulate_sampled), design::linearize feeding design::lqr,
// and the ExtendedKalmanFilter and UnscentedKalmanFilter estimating from a
// noisy h. The model carries analytic Jacobians; wrapping it in WithoutJacobian
// hides them, and the central differences must agree.

namespace {

// Damped pendulum, θ from the upright: θ̈ = (g/l)·sin θ − c·θ̇ + b·u.
struct Pendulum {
    double g_over_l{9.81}; ///< [1/s²]
    double c{0.5};         ///< [1/s] viscous damping
    double b{1.0};         ///< [rad/(s²·N·m)] torque gain

    [[nodiscard]] ColVec<2> operator()(const ColVec<2>& x, const ColVec<1>& u) const {
        return {x[1], (g_over_l * damp::sin(x[0])) - (c * x[1]) + (b * u[0])};
    }
    [[nodiscard]] DynamicsJacobian<2, 1> jacobian(const ColVec<2>& x, const ColVec<1>& /*u*/) const {
        return {.A = {{0.0, 1.0}, {g_over_l * damp::cos(x[0]), -c}}, .B = {{0.0}, {b}}};
    }
};

// Horizontal position of the bob, y = l·sin θ.
struct BobPosition {
    double length{1.0}; ///< [m]

    [[nodiscard]] ColVec<1> operator()(const ColVec<2>& x, const ColVec<1>& /*u*/) const {
        return ColVec<1>{length * damp::sin(x[0])};
    }
    [[nodiscard]] MeasurementJacobian<2, 1, 1> jacobian(const ColVec<2>& x, const ColVec<1>& /*u*/) const {
        return {.C = {{length * damp::cos(x[0]), 0.0}}};
    }
};

static_assert(AnalyticDynamics<Pendulum, 2, 1, double>);
static_assert(AnalyticMeasurement<BobPosition, 2, 1, 1, double>);

// A time-varying f(t, x, u) is not a Dynamics (the simulators still take it).
constexpr auto time_varying = [](double t, const ColVec<2>& x, const ColVec<1>& /*u*/) { return ColVec<2>{x[1], t}; };
static_assert(!Dynamics<decltype(time_varying), 2, 1, double>);

// The same f or h, without its Jacobian.
template<typename Model>
struct WithoutJacobian {
    Model model{};

    [[nodiscard]] auto operator()(const ColVec<2>& x, const ColVec<1>& u) const { return model(x, u); }
};

template<size_t R, size_t C>
[[nodiscard]] double max_abs_diff(const Matrix<R, C>& a, const Matrix<R, C>& b) {
    double m = 0.0;
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            m = damp::max(m, damp::abs(a(i, j) - b(i, j)));
        }
    }
    return m;
}

} // namespace

TEST_SUITE("Dynamics model") {
    TEST_CASE("analytic and central-difference Jacobians agree") {
        const Pendulum                     pendulum{};
        const BobPosition                  sensor{};
        const WithoutJacobian<Pendulum>    f_fd{pendulum};
        const WithoutJacobian<BobPosition> h_fd{sensor};
        static_assert(Dynamics<decltype(f_fd), 2, 1, double> && !AnalyticDynamics<decltype(f_fd), 2, 1, double>);
        static_assert(Measurement<decltype(h_fd), 2, 1, 1, double> && !AnalyticMeasurement<decltype(h_fd), 2, 1, 1, double>);

        const ColVec<2> points[] = {{0.0, 0.0}, {0.3, -0.5}, {2.0, 1.5}, {-3.0, 4.0}};
        for (const ColVec<2>& x : points) {
            const ColVec<1> u{0.4};
            const auto      exact = design::linearize<2, 1, 1>(pendulum, sensor, x, u);
            const auto      fd = design::linearize<2, 1, 1>(f_fd, h_fd, x, u);

            // A central difference with the default step ∛eps carries a relative
            // error of order ∛eps² against the Jacobian's own scale.
            const double scale = exact.A.norm() + exact.B.norm() + exact.C.norm();
            const double tol = 10.0 * damp::cbrt(eps<double>()) * damp::cbrt(eps<double>()) * scale;
            CHECK(max_abs_diff(exact.A, fd.A) <= tol);
            CHECK(max_abs_diff(exact.B, fd.B) <= tol);
            CHECK(max_abs_diff(exact.C, fd.C) <= tol);
            CHECK(max_abs_diff(exact.D, fd.D) <= tol);
            CHECK(max_abs_diff(exact.y_op, fd.y_op) == 0.0);
        }
    }

    TEST_CASE("one pendulum model drives the SIL plant, linearize + LQR, the EKF and the UKF") {
        const Pendulum    pendulum{};
        const BobPosition sensor{};
        const double      Ts = 0.01;

        // Design: linearize at the upright (sin θ ≈ θ), then a digital LQR.
        const auto lin = design::linearize<2, 1, 1>(pendulum, sensor, ColVec<2>{}, ColVec<1>{});
        CHECK(lin.A(1, 0) == doctest::Approx(pendulum.g_over_l));
        CHECK(lin.C(0, 0) == doctest::Approx(sensor.length));
        const auto res = design::lqr(lin.to_state_space(), Ts, {.Q = Matrix<2, 2>::diagonal({10.0, 1.0}), .R = Matrix<1, 1>{{0.1}}});
        REQUIRE(res.success);
        LQR law{res};

        // Estimators on the same model; the bob position is read with 1 cm noise.
        const double       sigma = 0.01;
        const Matrix<2, 2> Q = Matrix<2, 2>::diagonal({1e-8, 1e-6});
        const Matrix<2, 2> P0 = Matrix<2, 2>::diagonal({0.1, 0.1});

        ExtendedKalmanFilter<Pendulum, double>                  ekf{pendulum, {.Ts = Ts, .Q = Q, .P0 = P0}};
        ExtendedKalmanFilter<WithoutJacobian<Pendulum>, double> ekf_fd{WithoutJacobian<Pendulum>{pendulum}, {.Ts = Ts, .Q = Q, .P0 = P0}}; // central differences
        UnscentedKalmanFilter<Pendulum, double>                 ukf{pendulum, {.Ts = Ts, .Q = Q, .P0 = P0}};
        const Matrix<1, 1>                                      R{sigma * sigma};
        const WithoutJacobian<BobPosition>                      h_fd{sensor};

        std::mt19937                     rng(20261027);
        std::normal_distribution<double> noise(0.0, sigma);
        ColVec<1>                        u_held{};
        double                           ekf_sq[2] = {0.0, 0.0};
        double                           ukf_sq[2] = {0.0, 0.0};
        double                           fd_gap = 0.0;
        int                              n = 0;

        // Each tick: measure, update, control on the EKF estimate, predict.
        const auto controller = [&](double t, const ColVec<2>& x_true) -> ColVec<1> {
            const ColVec<1> y{sensor(x_true, u_held)[0] + noise(rng)};
            REQUIRE(ekf.update(sensor, y, R));
            REQUIRE(ekf_fd.update(h_fd, y, R));
            REQUIRE(ukf.update(sensor, y, R));
            if (t >= 1.0) {
                for (size_t i = 0; i < 2; ++i) {
                    ekf_sq[i] += (ekf.state()[i] - x_true[i]) * (ekf.state()[i] - x_true[i]);
                    ukf_sq[i] += (ukf.state()[i] - x_true[i]) * (ukf.state()[i] - x_true[i]);
                    fd_gap = damp::max(fd_gap, damp::abs(ekf_fd.state()[i] - ekf.state()[i]));
                }
                ++n;
            }
            u_held = law.step(ekf.state());
            ekf.predict(u_held);
            ekf_fd.predict(u_held);
            REQUIRE(ukf.predict(u_held));
            return u_held;
        };

        // Truth: the same Pendulum object, integrated with RK4 at Ts/10.
        const auto           full_state = [](const ColVec<2>& x) { return x; };
        sim::RK4<2>          rk4;
        sim::FixedStepSolver solver(rk4, Ts / 10.0);
        const auto           run = sim::simulate_sampled<2, 1, 2>(pendulum, full_state, controller, solver, Ts, ColVec<2>{0.3, 0.0}, {0.0, 5.0}, 10);

        // The LQR, fed the EKF estimate, brings the pendulum from 0.3 rad to the upright.
        CHECK(damp::abs(run.x.back()[0]) < 0.01);
        CHECK(damp::abs(run.x.back()[1]) < 0.05);

        // Both filters track the truth, the unmeasured rate included.
        REQUIRE(n > 0);
        const double count = static_cast<double>(n);
        CHECK(damp::sqrt(ekf_sq[0] / count) < sigma);
        CHECK(damp::sqrt(ukf_sq[0] / count) < sigma);
        CHECK(damp::sqrt(ekf_sq[1] / count) < 5.0 * sigma);
        CHECK(damp::sqrt(ukf_sq[1] / count) < 5.0 * sigma);

        // Central differences change the EKF only at the level of their own error.
        CHECK(fd_gap <= 1e3 * damp::cbrt(eps<double>()) * damp::cbrt(eps<double>()) * sigma);
    }
}
