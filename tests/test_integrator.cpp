// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>

#include "damp/matrix/matrix.hpp"
#include "damp/simulation/integrator.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

// Exercises the fixed-step integrators in damp::sim, with emphasis on the
// implicit LTI paths (BackwardEuler / BDF2 / Trapezoidal), which solve
// (I − αhA)·x_next = rhs each step. These overloads are templated and were
// previously never instantiated, so this suite both compiles and numerically
// validates them.
TEST_SUITE("integrator") {
    // -----------------------------------------------------------------------
    // AdaptiveStepIntegrator concept: only embedded pairs opt in.
    // -----------------------------------------------------------------------
    static_assert(sim::AdaptiveStepIntegrator<sim::DP45<1, double>, 1, double>);
    static_assert(sim::AdaptiveStepIntegrator<sim::RK23<1, double>, 1, double>);
    static_assert(sim::AdaptiveStepIntegrator<sim::TRBDF2<1, double>, 1, double>);
    static_assert(!sim::AdaptiveStepIntegrator<sim::RK4<1, double>, 1, double>);
    static_assert(!sim::AdaptiveStepIntegrator<sim::ForwardEuler<1, double>, 1, double>);
    static_assert(!sim::AdaptiveStepIntegrator<sim::BDF2<1, double>, 1, double>);
    static_assert(sim::DP45<1, double>::error_order == 5);
    static_assert(sim::RK23<1, double>::error_order == 3);
    static_assert(sim::TRBDF2<1, double>::error_order == 3);

    TEST_CASE("implicit LTI one-step matches the closed form (dx/dt = -x)") {
        // Scalar decay: A = -1, B = 0. Exact one-step factors are known.
        const Matrix<1, 1, double> A{{-1.0}};
        const Matrix<1, 1, double> B{{0.0}};
        const ColVec<1, double>    x0{1.0};
        const ColVec<1, double>    u{0.0};
        const double               h = 0.1;

        sim::BackwardEuler<1, double> be;
        // (1 + h)·x_next = x  ->  x_next = 1 / (1 + h).
        CHECK(be.evolve(A, B, x0, u, h).x[0] == doctest::Approx(1.0 / (1.0 + h)));

        sim::Trapezoidal<1, double> tr;
        // (1 + h/2)·x_next = (1 - h/2)·x.
        CHECK(tr.evolve(A, B, x0, u, h).x[0] == doctest::Approx((1.0 - (0.5 * h)) / (1.0 + 0.5 * h)));

        sim::BDF2<1, double> bdf;
        // First BDF2 step falls back to Backward Euler.
        CHECK(bdf.evolve(A, B, x0, u, h).x[0] == doctest::Approx(1.0 / (1.0 + h)));

        // TR-BDF2 one step (Hosea–Shampine, γ = 2−√2) vs closed form.
        const double           g = 2.0 - std::sqrt(2.0);
        const double           gh2 = g * h * 0.5;
        const double           xg = (1.0 - gh2) / (1.0 + gh2); // stage-1 trapezoidal factor
        const double           c_g = 1.0 / (g * (2.0 - g));
        const double           c_n = -((1.0 - g) * (1.0 - g)) / (g * (2.0 - g));
        const double           x_trbdf2 = (c_g * xg + c_n) / (1.0 + gh2); // A = -1, B = 0
        sim::TRBDF2<1, double> trbdf2;
        CHECK(trbdf2.evolve(A, B, x0, u, h).x[0] == doctest::Approx(x_trbdf2));
    }

    TEST_CASE("implicit LTI integrators converge to e^{-t}") {
        const Matrix<1, 1, double> A{{-1.0}};
        const Matrix<1, 1, double> B{{0.0}};
        const ColVec<1, double>    u{0.0};
        const double               h = 0.001;
        const int                  N = 1000; // integrate to t = 1
        const double               exact = std::exp(-1.0);

        SUBCASE("Backward Euler (1st order)") {
            sim::BackwardEuler<1, double> be;
            ColVec<1, double>             x{1.0};
            for (int k = 0; k < N; ++k) {
                x = be.evolve(A, B, x, u, h).x;
            }
            CHECK(x[0] == doctest::Approx(exact).epsilon(1e-3));
        }
        SUBCASE("Trapezoidal (2nd order)") {
            sim::Trapezoidal<1, double> tr;
            ColVec<1, double>           x{1.0};
            for (int k = 0; k < N; ++k) {
                x = tr.evolve(A, B, x, u, h).x;
            }
            CHECK(x[0] == doctest::Approx(exact).epsilon(1e-5));
        }
        SUBCASE("BDF2 (2nd order, multistep)") {
            sim::BDF2<1, double> bdf;
            ColVec<1, double>    x{1.0};
            for (int k = 0; k < N; ++k) {
                x = bdf.evolve(A, B, x, u, h).x;
            }
            CHECK(x[0] == doctest::Approx(exact).epsilon(1e-4));
        }
        SUBCASE("TR-BDF2 (2nd order, one-step L-stable)") {
            sim::TRBDF2<1, double> trbdf2;
            ColVec<1, double>      x{1.0};
            for (int k = 0; k < N; ++k) {
                x = trbdf2.evolve(A, B, x, u, h).x;
            }
            CHECK(x[0] == doctest::Approx(exact).epsilon(1e-4));
        }
    }

    TEST_CASE("TR-BDF2 damps a stiff mode without trapezoidal ringing") {
        // dx/dt = -1000 x, h = 0.1 ⇒ |λh| = 100. Trapezoidal |R| ≈ 1 (rings);
        // TR-BDF2 is L-stable so |x| shrinks well below 1 after one step.
        const Matrix<1, 1, double> A{{-1000.0}};
        const Matrix<1, 1, double> B{{0.0}};
        const ColVec<1, double>    x0{1.0};
        const ColVec<1, double>    u{0.0};
        const double               h = 0.1;

        sim::Trapezoidal<1, double> tr;
        const double                x_tr = tr.evolve(A, B, x0, u, h).x[0];
        CHECK(std::abs(x_tr) > 0.9); // A-stable but |R(∞)| = 1 → near-unit oscillation

        sim::TRBDF2<1, double> trbdf2;
        const auto             step = trbdf2.evolve(A, B, x0, u, h);
        CHECK(std::abs(step.x[0]) < 0.1); // L-stable composite damps the stiff mode
        CHECK(step.error.norm() > 0.0);   // embedded estimate is populated
    }

    TEST_CASE("TR-BDF2 embedded error estimate scales as O(h^3)") {
        const Matrix<1, 1, double> A{{-1.0}};
        const Matrix<1, 1, double> B{{0.0}};
        const ColVec<1, double>    x{1.0};
        const ColVec<1, double>    u{0.0};

        sim::TRBDF2<1, double> trbdf2;
        const double           e1 = trbdf2.evolve(A, B, x, u, 0.1).error.norm();
        const double           e2 = trbdf2.evolve(A, B, x, u, 0.05).error.norm();

        CHECK(e1 > 0.0);
        CHECK(e2 < e1 / 4.0); // ~O(h^3) → factor ~8 on halving; allow slack
        const double true_err = std::abs(trbdf2.evolve(A, B, x, u, 0.1).x[0] - std::exp(-0.1));
        // Estimate and true local error should be the same order of magnitude.
        CHECK(e1 < 10.0 * true_err);
        CHECK(true_err < 10.0 * e1);
    }

    TEST_CASE("implicit LTI integrators track a forced steady state (x_ss = u)") {
        // dx/dt = -x + u  ->  steady state x = u.
        const Matrix<1, 1, double> A{{-1.0}};
        const Matrix<1, 1, double> B{{1.0}};
        const ColVec<1, double>    u{2.0};
        const double               h = 0.01;

        sim::BackwardEuler<1, double> be;
        ColVec<1, double>             x{0.0};
        for (int k = 0; k < 3000; ++k) {
            x = be.evolve(A, B, x, u, h).x;
        }
        CHECK(x[0] == doctest::Approx(2.0).epsilon(1e-3));
    }

    TEST_CASE("implicit 2-state solve agrees with the exact (matrix-exponential) integrator") {
        // Damped oscillator: stable, non-diagonal A exercises the 2x2 linear solve.
        const Matrix<2, 2, double> A{{0.0, 1.0}, {-4.0, -0.5}};
        const Matrix<2, 1, double> B{{0.0}, {0.0}};
        const ColVec<1, double>    u{0.0};
        const ColVec<2, double>    x0{1.0, 0.0};
        const double               h = 0.0005;
        const int                  N = 2000; // t = 1

        // Reference: exact discrete update via the matrix exponential.
        sim::Exact<2, double> ex;
        ColVec<2, double>     xref = x0;
        for (int k = 0; k < N; ++k) {
            xref = ex.evolve(A, B, xref, u, h).x;
        }

        sim::BackwardEuler<2, double> be;
        ColVec<2, double>             xbe = x0;
        for (int k = 0; k < N; ++k) {
            xbe = be.evolve(A, B, xbe, u, h).x;
        }
        CHECK(xbe[0] == doctest::Approx(xref[0]).epsilon(2e-2));
        CHECK(xbe[1] == doctest::Approx(xref[1]).epsilon(2e-2));

        sim::Trapezoidal<2, double> tr;
        ColVec<2, double>           xtr = x0;
        for (int k = 0; k < N; ++k) {
            xtr = tr.evolve(A, B, xtr, u, h).x;
        }
        CHECK(xtr[0] == doctest::Approx(xref[0]).epsilon(1e-3));
        CHECK(xtr[1] == doctest::Approx(xref[1]).epsilon(1e-3));
    }

    TEST_CASE("implicit LTI integrator is constexpr-evaluable") {
        constexpr double y = []() consteval {
            Matrix<1, 1, double>          A{{-1.0}};
            Matrix<1, 1, double>          B{{0.0}};
            ColVec<1, double>             x{1.0};
            ColVec<1, double>             u{0.0};
            sim::BackwardEuler<1, double> be;
            for (int k = 0; k < 10; ++k) {
                x = be.evolve(A, B, x, u, 0.1).x;
            }
            return x[0];
        }();
        // 10 Backward-Euler steps of h = 0.1: x = 1 / 1.1^10.
        static_assert(y > 0.38 && y < 0.39, "implicit integrator must work at compile time");
        CHECK(y == doctest::Approx(1.0 / std::pow(1.1, 10)));
    }

    // -----------------------------------------------------------------------
    // Explicit integrators: ForwardEuler, Heun, RK3, RK23, DP45, SymplecticEuler
    // All tested on scalar decay dx/dt = -x, x(0) = 1, exact solution e^{-1}.
    // -----------------------------------------------------------------------

    TEST_CASE("ForwardEuler LTI converges to e^{-1}") {
        const Matrix<1, 1, double> A{{-1.0}};
        const Matrix<1, 1, double> B{{0.0}};
        const ColVec<1, double>    u{0.0};
        const double               h = 0.001;
        const int                  N = 1000;
        const double               exact = std::exp(-1.0);

        sim::ForwardEuler<1, double> fe;
        ColVec<1, double>            x{1.0};
        for (int k = 0; k < N; ++k) {
            x = fe.evolve(A, B, x, u, h).x;
        }
        // Forward Euler is 1st-order; with h = 0.001 expect ~1e-3 relative error.
        CHECK(x[0] == doctest::Approx(exact).epsilon(1e-2));
    }

    TEST_CASE("ForwardEuler nonlinear converges to e^{-1}") {
        // f(t, x) = -x, same decay, nonlinear overload.
        auto f = [](double, const ColVec<1, double>& x) {
            return ColVec<1, double>{-x[0]};
        };
        const double h = 0.001;
        const int    N = 1000;
        const double exact = std::exp(-1.0);

        sim::ForwardEuler<1, double> fe;
        ColVec<1, double>            x{1.0};
        for (int k = 0; k < N; ++k) {
            x = fe.evolve(f, x, k * h, h).x;
        }
        CHECK(x[0] == doctest::Approx(exact).epsilon(1e-2));
    }

    TEST_CASE("Heun LTI converges to e^{-1}") {
        const Matrix<1, 1, double> A{{-1.0}};
        const Matrix<1, 1, double> B{{0.0}};
        const ColVec<1, double>    u{0.0};
        const double               h = 0.01;
        const int                  N = 100;
        const double               exact = std::exp(-1.0);

        sim::Heun<1, double> heun;
        ColVec<1, double>    x{1.0};
        for (int k = 0; k < N; ++k) {
            x = heun.evolve(A, B, x, u, h).x;
        }
        // Heun is 2nd-order; with h = 0.01 expect ~1e-4 error.
        CHECK(x[0] == doctest::Approx(exact).epsilon(1e-4));
    }

    TEST_CASE("Heun nonlinear: dx/dt = cos(t) -> sin(t)") {
        // Integrating dx/dt = cos(t) from t=0 to t=1 gives x(1) = sin(1).
        auto f = [](double t, const ColVec<1, double>&) {
            return ColVec<1, double>{std::cos(t)};
        };
        const double h = 0.001;
        const int    N = 1000;

        sim::Heun<1, double> heun;
        ColVec<1, double>    x{0.0};
        for (int k = 0; k < N; ++k) {
            x = heun.evolve(f, x, k * h, h).x;
        }
        CHECK(x[0] == doctest::Approx(std::sin(1.0)).epsilon(1e-5));
    }

    TEST_CASE("RK3 LTI converges to e^{-1}") {
        const Matrix<1, 1, double> A{{-1.0}};
        const Matrix<1, 1, double> B{{0.0}};
        const ColVec<1, double>    u{0.0};
        const double               h = 0.01;
        const int                  N = 100;
        const double               exact = std::exp(-1.0);

        sim::RK3<1, double> rk3;
        ColVec<1, double>   x{1.0};
        for (int k = 0; k < N; ++k) {
            x = rk3.evolve(A, B, x, u, h).x;
        }
        // RK3 is 3rd-order; with h = 0.01 expect ~1e-6 error.
        CHECK(x[0] == doctest::Approx(exact).epsilon(1e-6));
    }

    TEST_CASE("RK3 nonlinear: dx/dt = -x") {
        auto f = [](double, const ColVec<1, double>& x) {
            return ColVec<1, double>{-x[0]};
        };
        const double h = 0.01;
        const int    N = 100;
        const double exact = std::exp(-1.0);

        sim::RK3<1, double> rk3;
        ColVec<1, double>   x{1.0};
        for (int k = 0; k < N; ++k) {
            x = rk3.evolve(f, x, k * h, h).x;
        }
        CHECK(x[0] == doctest::Approx(exact).epsilon(1e-6));
    }

    TEST_CASE("RK23 LTI converges to e^{-1} and populates error estimate") {
        const Matrix<1, 1, double> A{{-1.0}};
        const Matrix<1, 1, double> B{{0.0}};
        const ColVec<1, double>    u{0.0};
        const double               h = 0.01;
        const int                  N = 100;
        const double               exact = std::exp(-1.0);

        sim::RK23<1, double> rk23;
        ColVec<1, double>    x{1.0};
        double               last_error = 0.0;
        for (int k = 0; k < N; ++k) {
            auto res = rk23.evolve(A, B, x, u, h);
            x = res.x;
            last_error = res.error.norm();
        }
        // RK23 returns the 3rd-order solution; with h = 0.01 expect ~1e-6.
        CHECK(x[0] == doctest::Approx(exact).epsilon(1e-6));
        // Error estimate must be non-negative (it is a norm).
        CHECK(last_error >= 0.0);
    }

    TEST_CASE("RK23 nonlinear: dx/dt = -x") {
        auto f = [](double, const ColVec<1, double>& x) {
            return ColVec<1, double>{-x[0]};
        };
        const double h = 0.01;
        const int    N = 100;
        const double exact = std::exp(-1.0);

        sim::RK23<1, double> rk23;
        ColVec<1, double>    x{1.0};
        for (int k = 0; k < N; ++k) {
            x = rk23.evolve(f, x, k * h, h).x;
        }
        CHECK(x[0] == doctest::Approx(exact).epsilon(1e-6));
    }

    TEST_CASE("DP45 LTI converges to e^{-1} with high accuracy") {
        const Matrix<1, 1, double> A{{-1.0}};
        const Matrix<1, 1, double> B{{0.0}};
        const ColVec<1, double>    u{0.0};
        const double               h = 0.1;
        const int                  N = 10;
        const double               exact = std::exp(-1.0);

        sim::DP45<1, double> dp45;
        ColVec<1, double>    x{1.0};
        for (int k = 0; k < N; ++k) {
            x = dp45.evolve(A, B, x, u, h).x;
        }
        // DP45 advances with the 5th-order solution; with h = 0.1 (10 steps to t = 1) expect ~1e-9.
        CHECK(x[0] == doctest::Approx(exact).epsilon(1e-9));
    }

    TEST_CASE("DP45 nonlinear: dx/dt = cos(t) -> sin(t)") {
        auto f = [](double t, const ColVec<1, double>&) {
            return ColVec<1, double>{std::cos(t)};
        };
        const double h = 0.1;
        const int    N = 10;

        sim::DP45<1, double> dp45;
        ColVec<1, double>    x{0.0};
        for (int k = 0; k < N; ++k) {
            x = dp45.evolve(f, x, k * h, h).x;
        }
        CHECK(x[0] == doctest::Approx(std::sin(1.0)).epsilon(1e-9));
    }

    TEST_CASE("DP45 embedded error estimate tracks the true local error") {
        // dx/dt = -x from x = 1: local truncation error of the 5th-order solution over
        // one step is ~h⁶/6! ~ 1e-9 at h = 0.1; the embedded 4(5) difference estimates
        // the 4th-order member's error (~h⁵), so it must be small, positive, and within
        // a couple of orders of magnitude of that scale — and shrink ~32x when h halves.
        const Matrix<1, 1, double> A{{-1.0}};
        const Matrix<1, 1, double> B{{0.0}};
        const ColVec<1, double>    u{0.0};

        sim::DP45<1, double>    dp45;
        const ColVec<1, double> x{1.0};

        const double e1 = dp45.evolve(A, B, x, u, 0.1).error.norm();
        const double e2 = dp45.evolve(A, B, x, u, 0.05).error.norm();

        CHECK(e1 > 0.0);
        CHECK(e1 < 1e-6);      // sane magnitude for h = 0.1
        CHECK(e2 < e1 / 16.0); // ~O(h^5) scaling on halving
        const double true_err = std::abs(dp45.evolve(A, B, x, u, 0.1).x[0] - std::exp(-0.1));
        CHECK(true_err < e1); // 5th-order solution beats the 4th-order error estimate
    }

    TEST_CASE("SymplecticEuler advances a harmonic oscillator without secular energy growth") {
        // ÿ + y = 0 as first-order system x = [q, v], A = [[0,1],[-1,0]].
        // Forward Euler energy grows; symplectic Euler keeps it O(1) bounded.
        const Matrix<2, 2, double> A{{0.0, 1.0}, {-1.0, 0.0}};
        const Matrix<2, 1, double> B{{0.0}, {0.0}};
        const ColVec<1, double>    u{0.0};
        const double               h = 0.01;
        const int                  N = 10000; // t = 100 (many periods)

        auto energy = [](const ColVec<2, double>& x) {
            return 0.5 * (x[0] * x[0] + x[1] * x[1]);
        };

        sim::SymplecticEuler<2, double> se;
        ColVec<2, double>               xse{1.0, 0.0};
        const double                    E0 = energy(xse);
        double                          E_max_se = E0;
        double                          E_min_se = E0;
        for (int k = 0; k < N; ++k) {
            xse = se.evolve(A, B, xse, u, h).x;
            const double E = energy(xse);
            E_max_se = std::max(E_max_se, E);
            E_min_se = std::min(E_min_se, E);
        }
        // Energy oscillates but stays within a few percent of E0 over long horizon.
        CHECK(E_max_se < E0 * 1.05);
        CHECK(E_min_se > E0 * 0.95);

        sim::ForwardEuler<2, double> fe;
        ColVec<2, double>            xfe{1.0, 0.0};
        for (int k = 0; k < N; ++k) {
            xfe = fe.evolve(A, B, xfe, u, h).x;
        }
        // Forward Euler has secular energy growth on the same horizon.
        CHECK(energy(xfe) > E0 * 1.5);
    }

    TEST_CASE("explicit integrators are constexpr-evaluable (ForwardEuler)") {
        constexpr double y = []() consteval {
            Matrix<1, 1, double>         A{{-1.0}};
            Matrix<1, 1, double>         B{{0.0}};
            ColVec<1, double>            x{1.0};
            ColVec<1, double>            u{0.0};
            sim::ForwardEuler<1, double> fe;
            for (int k = 0; k < 10; ++k) {
                x = fe.evolve(A, B, x, u, 0.1).x;
            }
            return x[0];
        }();
        // 10 Forward-Euler steps of h = 0.1: x = (1 - 0.1)^10 = 0.9^10.
        static_assert(y > 0.34 && y < 0.36, "ForwardEuler must work at compile time");
        CHECK(y == doctest::Approx(std::pow(0.9, 10)));
    }
}
