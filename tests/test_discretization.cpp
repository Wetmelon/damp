// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <numbers>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/design/stability.hpp"
#include "damp/design/status.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

/**
 * @brief Tests for discretization methods (ZOH, Tustin)
 *
 * These tests verify the accuracy and correctness of continuous-to-discrete
 * system transformations.
 */

TEST_SUITE("Control Design: Discretization") {
    // Test 1: Simple 1st-order system discretization via ZOH
    // System: dx/dt = -x + u (time constant = 1 second)
    // Analytical solution: x[k+1] = exp(-Ts)*x[k] + (1-exp(-Ts))*u[k]
    TEST_CASE("ZOH Discretization: 1st-order system") {
        double     Ts = 0.1; // 100ms sampling time
        StateSpace sys{
            Matrix<1, 1>{{-1.0}}, // A: dx/dt = -x
            Matrix<1, 1>{{1.0}},  // B: input gain
            Matrix<1, 1>{{1.0}}   // C: output is state
        };

        // Use new unified discretize() function
        StateSpace sys_d = *discretize(sys, Ts, DiscretizationMethod::ZOH);

        // Verify A_d ≈ exp(-Ts)
        double A_d_expected = std::exp(-Ts);
        double A_d_actual = sys_d.A(0, 0);
        CHECK(doctest::Approx(A_d_actual).epsilon(1e-4) == A_d_expected);

        // Verify B_d ≈ 1 - exp(-Ts)
        double B_d_expected = 1.0 - std::exp(-Ts);
        double B_d_actual = sys_d.B(0, 0);
        CHECK(doctest::Approx(B_d_actual).epsilon(1e-4) == B_d_expected);
    }

    // Test 2: Tustin Discretization
    TEST_CASE("Tustin Discretization: 1st-order system") {
        double     Ts = 0.1;
        StateSpace sys{
            Matrix<1, 1>{{-1.0}},
            Matrix<1, 1>{{1.0}},
            Matrix<1, 1>{{1.0}}
        };

        // Use new unified discretize() function
        StateSpace sys_d = *discretize(sys, Ts, DiscretizationMethod::Tustin);

        // For Tustin with A=-1, Ts=0.1:
        // M = (I - A*Ts/2)^{-1} = (1 + 0.05)^{-1} = 1/1.05
        // A_d = M * (I + A*Ts/2) = (1/1.05) * 0.95 = (2-Ts)/(2+Ts)
        // B_d = Ts * M * B = 0.1/1.05
        // C_d = C * M = 1/1.05
        // D_d = D + (Ts/2) * C_d * B = 0 + 0.05/1.05
        double expected_A_d = (2.0 - Ts) / (2.0 + Ts);       // 0.904762...
        double expected_B_d = Ts / (1.0 + Ts / 2.0);         // 0.095238...
        double expected_C_d = 1.0 / (1.0 + Ts / 2.0);        // 0.952381...
        double expected_D_d = (Ts / 2.0) / (1.0 + Ts / 2.0); // 0.047619...

        CHECK(sys_d.A(0, 0) == doctest::Approx(expected_A_d).epsilon(1e-12));
        CHECK(sys_d.B(0, 0) == doctest::Approx(expected_B_d).epsilon(1e-12));
        CHECK(sys_d.C(0, 0) == doctest::Approx(expected_C_d).epsilon(1e-12));
        CHECK(sys_d.D(0, 0) == doctest::Approx(expected_D_d).epsilon(1e-12));

        // Discrete pole should be inside unit circle (stable system stays stable)
        CHECK(std::abs(sys_d.A(0, 0)) < 1.0);
    }

    // Test 3: ZOH vs Tustin comparison using unified API
    TEST_CASE("ZOH and Tustin produce different discretizations") {
        double     Ts = 0.5;
        StateSpace sys{
            Matrix<1, 1>{{-2.0}},
            Matrix<1, 1>{{1.0}},
            Matrix<1, 1>{{1.0}}
        };

        StateSpace sys_d_zoh = *discretize(sys, Ts, DiscretizationMethod::ZOH);
        StateSpace sys_d_tustin = *discretize(sys, Ts, DiscretizationMethod::Tustin);

        // They should produce different results
        CHECK_NE(sys_d_zoh.A(0, 0), doctest::Approx(sys_d_tustin.A(0, 0)));
    }

    // Test 4: Default method is ZOH
    TEST_CASE("eval_frf returns nullopt at a continuous pole") {
        StateSpace<1, 1, 1> sys{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}};
        // s = -1 is a pole of 1/(s+1)
        CHECK_FALSE(eval_frf(sys, damp::complex<double>{-1.0, 0.0}).has_value());
        // s = j is fine
        CHECK(eval_frf(sys, damp::complex<double>{0.0, 1.0}).has_value());
    }

    TEST_CASE("Default discretization method is ZOH") {
        double     Ts = 0.1;
        StateSpace sys{
            Matrix<1, 1>{{-1.0}},
            Matrix<1, 1>{{1.0}},
            Matrix<1, 1>{{1.0}}
        };

        StateSpace sys_d_default = *discretize(sys, Ts);
        StateSpace sys_d_zoh = *discretize(sys, Ts, DiscretizationMethod::ZOH);

        CHECK(doctest::Approx(sys_d_default.A(0, 0)) == sys_d_zoh.A(0, 0));
        CHECK(doctest::Approx(sys_d_default.B(0, 0)) == sys_d_zoh.B(0, 0));
    }

    // ZOH maps process-noise input G with the same ∫e^{Aτ}G dτ as B
    TEST_CASE("ZOH G_d matches B_d integral for identical columns") {
        double                            Ts = 0.1;
        StateSpace<1, 1, 1, double, 1, 0> sys{
            .A = {{-1.0}},
            .B = {{1.0}},
            .C = {{1.0}},
            .D = {{0.0}},
            .G = {{1.0}}, // same column as B
        };

        auto sys_d = *discretize(sys, Ts, DiscretizationMethod::ZOH);

        // G_d must equal B_d, not the Forward-Euler G*Ts
        CHECK(sys_d.G(0, 0) == doctest::Approx(sys_d.B(0, 0)).epsilon(1e-12));
        CHECK(sys_d.G(0, 0) == doctest::Approx(1.0 - std::exp(-Ts)).epsilon(1e-10));
        CHECK(std::abs(sys_d.G(0, 0) - Ts) > 1e-6); // not FE
    }

    TEST_CASE("ForwardEuler: Ad = I + A Ts, Bd = B Ts") {
        const double Ts = 0.1;
        StateSpace   sys{
            .A = Matrix<1, 1>{{-1.0}},
            .B = Matrix<1, 1>{{1.0}},
            .C = Matrix<1, 1>{{1.0}},
        };
        const auto sys_d = *discretize(sys, Ts, DiscretizationMethod::ForwardEuler);
        CHECK(sys_d.A(0, 0) == doctest::Approx(1.0 - Ts)); // 0.9
        CHECK(sys_d.B(0, 0) == doctest::Approx(Ts));
        CHECK(sys_d.Ts == doctest::Approx(Ts));
    }

    TEST_CASE("discretize refuses Ts <= 0 and a discrete model at another rate, and passes one at the same rate") {
        const StateSpace<1, 1, 1> cont{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}};
        CHECK(discretize(cont, 0.0).status == design::DesignStatus::invalid_sample_time);
        CHECK(discretize(cont, -0.01).status == design::DesignStatus::invalid_sample_time);
        static_assert(!discretize(StateSpace<1, 1, 1>{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}}, -1.0));

        const StateSpace<1, 1, 1> disc{.A = {{-0.5}}, .B = {{0.1}}, .C = {{1.0}}, .Ts = 0.01};
        const auto                other_rate = discretize(disc, 0.001);
        CHECK_FALSE(other_rate);
        CHECK(other_rate.status == design::DesignStatus::invalid_sample_time);
        const auto same_rate = discretize(disc, 0.01, DiscretizationMethod::Tustin);
        REQUIRE(same_rate);
        CHECK(same_rate->Ts == 0.01);
        CHECK(same_rate->A(0, 0) == -0.5);
        CHECK(same_rate->B(0, 0) == 0.1);

        //! Ts < 0 (MATLAB's unspecified rate) or NaN is not a continuous model.
        StateSpace<1, 1, 1> unspecified = cont;
        unspecified.Ts = -1.0;
        CHECK(discretize(unspecified, 0.1).status == design::DesignStatus::invalid_sample_time);
        unspecified.Ts = std::numeric_limits<double>::quiet_NaN();
        CHECK(discretize(unspecified, 0.1).status == design::DesignStatus::invalid_sample_time);
    }

    TEST_CASE("discretize refuses a model with a non-finite entry, for every method") {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double inf = std::numeric_limits<double>::infinity();
        for (const auto method : {DiscretizationMethod::ForwardEuler, DiscretizationMethod::ZOH, DiscretizationMethod::Tustin}) {
            for (const double bad : {nan, inf}) {
                StateSpace<2, 1, 1, double, 1> in_A{.A = {{0.0, 1.0}, {-4.0, -1.0}}, .B = {{0.0}, {1.0}}, .C = {{1.0, 0.0}}};
                in_A.A(1, 1) = bad;
                StateSpace<2, 1, 1, double, 1> in_B = in_A;
                in_B.A(1, 1) = -1.0;
                in_B.B(1, 0) = bad;
                StateSpace<2, 1, 1, double, 1> in_G = in_B;
                in_G.B(1, 0) = 1.0;
                in_G.G(0, 0) = bad;
                CHECK(discretize(in_A, 0.01, method).status == design::DesignStatus::bad_weights);
                CHECK(discretize(in_B, 0.01, method).status == design::DesignStatus::bad_weights);
                CHECK(discretize(in_G, 0.01, method).status == design::DesignStatus::bad_weights);
                CHECK(discretize(in_A, 0.01, InputDelay{0.004}).status == design::DesignStatus::bad_weights);
            }
        }
    }

    TEST_CASE("discretize refuses a model that overflows within one period") {
        const StateSpace<1, 1, 1> fast{.A = {{1000.0}}, .B = {{1.0}}, .C = {{1.0}}};
        CHECK(discretize(fast, 10.0).status == design::DesignStatus::not_converged);
        const StateSpace<2, 1, 1> coupled{.A = {{1000.0, 1.0}, {0.0, 1000.0}}, .B = {{0.0}, {1.0}}, .C = {{1.0, 0.0}}};
        CHECK(discretize(coupled, 10.0).status == design::DesignStatus::not_converged);
        CHECK(discretize(coupled, 1e-3)); // the same plant at a rate it fits
    }

    TEST_CASE("ZOH of a double integrator uses the singular-A integral") {
        // A is nilpotent, so A⁻¹(e^{ATs}−I) is unavailable. Exact ZOH:
        // Ad = [1, Ts; 0, 1], Bd = [Ts²/2; Ts].
        const double Ts = 0.1;
        StateSpace   sys{
            .A = Matrix<2, 2>{{0.0, 1.0}, {0.0, 0.0}},
            .B = Matrix<2, 1>{{0.0}, {1.0}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
        };
        const auto sys_d = discretize(sys, Ts, DiscretizationMethod::ZOH);
        REQUIRE(sys_d);
        CHECK(sys_d->A(0, 0) == doctest::Approx(1.0));
        CHECK(sys_d->A(0, 1) == doctest::Approx(Ts));
        CHECK(sys_d->A(1, 0) == doctest::Approx(0.0));
        CHECK(sys_d->A(1, 1) == doctest::Approx(1.0));
        CHECK(sys_d->B(0, 0) == doctest::Approx(0.5 * Ts * Ts));
        CHECK(sys_d->B(1, 0) == doctest::Approx(Ts));
    }

    TEST_CASE("ZOH is accurate for a nearly-singular A") {
        // B_d = (e^{a Ts} − 1)/a ≈ Ts; the A⁻¹(e^{ATs} − I) form cancels.
        const StateSpace<1, 1, 1> sys{.A = {{-1e-12}}, .B = {{1.0}}, .C = {{1.0}}};
        const auto                d = discretize(sys, 1e-3, DiscretizationMethod::ZOH);
        REQUIRE(d.has_value());
        CHECK(d->B(0, 0) == doctest::Approx(1e-3).epsilon(1e-12));
    }
}

TEST_SUITE("DiscretizationResult") {
    using design::DesignStatus;

    // A lag with a noise input on each side, so G and H ride through every path.
    constexpr StateSpace<2, 1, 1, double, 1, 1> make_lag() {
        return StateSpace<2, 1, 1, double, 1, 1>{
            .A = Matrix<2, 2>{{-1.0, 0.5}, {0.0, -3.0}},
            .B = Matrix<2, 1>{{0.0}, {2.0}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .D = Matrix<1, 1>{{0.25}},
            .G = Matrix<2, 1>{{1.0}, {0.5}},
            .H = Matrix<1, 1>{{0.1}},
        };
    }

    // Same model to rounding: -ffast-math may round an inlined copy of the method differently.
    template<size_t R, size_t C>
    void check_close(const Matrix<R, C>& a, const Matrix<R, C>& b) {
        for (size_t i = 0; i < R; ++i) {
            for (size_t j = 0; j < C; ++j) {
                CHECK(a(i, j) == doctest::Approx(b(i, j)).epsilon(1e-14).scale(1.0));
            }
        }
    }

    template<size_t NX, size_t NU, size_t NY, size_t NW, size_t NV>
    void check_model_close(const StateSpace<NX, NU, NY, double, NW, NV>& a, const StateSpace<NX, NU, NY, double, NW, NV>& b) {
        check_close(a.A, b.A);
        check_close(a.B, b.B);
        check_close(a.C, b.C);
        check_close(a.D, b.D);
        check_close(a.G, b.G);
        check_close(a.H, b.H);
        CHECK(a.Ts == b.Ts);
    }

    TEST_CASE("a success reads like an optional and holds the same matrices as the method") {
        constexpr auto sys = make_lag();
        constexpr auto zoh_c = discretize(sys, 0.1);
        static_assert(zoh_c.success && zoh_c.status == DesignStatus::ok && zoh_c.has_value());
        static_assert(static_cast<bool>(zoh_c) && !!zoh_c);
        static_assert(!std::is_convertible_v<decltype(zoh_c), bool>); // contextual only, like std::optional
        static_assert(zoh_c->Ts == 0.1 && (*zoh_c).A == zoh_c.sys.A && zoh_c.value().B == zoh_c.sys.B);

        //! Every accessor is the held model, and that model is the method's.
        const auto zoh = discretize(sys, 0.1);
        CHECK(&*zoh == &zoh.sys);
        CHECK(&zoh.value() == &zoh.sys);
        CHECK(zoh.operator->() == &zoh.sys);
        const auto direct = detail::discretize_zoh_impl(sys, 0.1).value();
        check_model_close(*zoh, direct);
        check_model_close(*zoh_c, direct);

        const auto euler = discretize(sys, 0.1, DiscretizationMethod::ForwardEuler);
        REQUIRE(euler);
        check_model_close(*euler, detail::discretize_forward_euler_impl(sys, 0.1));

        //! Tustin: the closed forms of a first-order lag (as in the test above).
        const StateSpace<1, 1, 1> lag{.A = {{-1.0}}, .B = {{1.0}}, .C = {{1.0}}};
        const auto                tustin = discretize(lag, 0.1, DiscretizationMethod::Tustin);
        REQUIRE(tustin.has_value());
        CHECK(tustin->A(0, 0) == doctest::Approx(1.9 / 2.1).epsilon(1e-14));
        CHECK(tustin->D(0, 0) == doctest::Approx(0.05 / 1.05).epsilon(1e-14));

        //! A whole-sample delay: the register model built on the plain ZOH model.
        const auto late = discretize(sys, 0.1, InputDelay{0.1});
        REQUIRE(late);
        CHECK(late.status == DesignStatus::ok);
        check_close(Matrix<2, 2>(late->A.block<2, 2>(0, 0)), direct.A);
        check_close(Matrix<2, 1>(late->A.block<2, 1>(0, 2)), direct.B);
        check_close(Matrix<2, 1>(late->G.block<2, 1>(0, 0)), direct.G);

        //! Mutable access writes through to the held model; as<U>() keeps the status.
        auto copy = zoh;
        copy->Ts = 0.2;
        (*copy).D(0, 0) = 1.0;
        CHECK(copy.sys.Ts == 0.2);
        CHECK(copy.value().D(0, 0) == 1.0);
        constexpr auto zoh_f = zoh_c.as<float>();
        static_assert(zoh_f.success && zoh_f.status == DesignStatus::ok);
        CHECK(zoh_f->A(0, 0) == doctest::Approx(direct.A(0, 0)).epsilon(1e-6));
    }

    TEST_CASE("each refusal names its cause") {
        constexpr auto sys = make_lag();
        const double   nan = std::numeric_limits<double>::quiet_NaN();
        const double   inf = std::numeric_limits<double>::infinity();

        //! Plain discretize.
        const StateSpace<1, 1, 1> pole_at_2_over_Ts{.A = {{4.0}}, .B = {{1.0}}, .C = {{1.0}}}; // I − A·Ts/2 = 0 at Ts = 0.5
        const auto                singular = discretize(pole_at_2_over_Ts, 0.5, DiscretizationMethod::Tustin);
        CHECK_FALSE(singular);
        CHECK_FALSE(singular.success);
        CHECK(singular.status == DesignStatus::singular);
        CHECK(singular->A(0, 0) == 0.0); // the zero model
        CHECK(discretize(sys, 0.01, DiscretizationMethod::Tustin, std::numbers::pi / 0.01).status == DesignStatus::invalid_parameter);
        CHECK(discretize(sys, 0.01, DiscretizationMethod::Tustin, -1.0).status == DesignStatus::invalid_parameter);
        CHECK(discretize(sys, 0.01, DiscretizationMethod::Tustin, nan).status == DesignStatus::invalid_parameter);
        CHECK(discretize(sys, 0.01, static_cast<DiscretizationMethod>(7)).status == DesignStatus::invalid_parameter);
        CHECK(discretize(sys, nan).status == DesignStatus::invalid_sample_time);
        CHECK(discretize(sys, inf, DiscretizationMethod::Tustin).status == DesignStatus::invalid_sample_time);
        const StateSpace<2, 1, 1> nan_plant{.A = {{0.0, 1.0}, {-4.0, nan}}, .B = {{0.0}, {1.0}}, .C = {{1.0, 0.0}}};
        CHECK(discretize(nan_plant, 0.01).status == DesignStatus::bad_weights); // expm used to return I: Ad = I, Bd = 0
        CHECK(discretize(nan_plant, 0.01, InputDelay{0.004}).status == DesignStatus::bad_weights);

        //! Input delay.
        const auto fast = *discretize(sys, 0.1);
        CHECK(discretize(sys, 0.0, InputDelay{0.0}).status == DesignStatus::invalid_sample_time);
        CHECK(discretize(sys, -0.1, InputDelay{0.0}).status == DesignStatus::invalid_sample_time);
        CHECK(discretize(sys, nan, InputDelay{0.0}).status == DesignStatus::invalid_sample_time);
        CHECK(discretize(sys, 0.1, InputDelay{-0.01}).status == DesignStatus::invalid_parameter);
        CHECK(discretize(sys, 0.1, InputDelay{nan}).status == DesignStatus::invalid_parameter);
        CHECK(discretize(sys, 0.1, InputDelay{0.15}).status == DesignStatus::invalid_parameter); // longer than one sample
        CHECK(discretize(sys, 0.1, InputDelay{0.04}, DiscretizationMethod::Tustin).status == DesignStatus::invalid_parameter);
        CHECK(discretize(fast, 0.1, InputDelay{0.04}).status == DesignStatus::invalid_sample_time); // a fraction needs the continuous plant
        CHECK(discretize(fast, 0.2, InputDelay{0.2}).status == DesignStatus::invalid_sample_time);  // a discrete model at another rate
        CHECK(discretize(fast, 0.1, InputDelay{0.1}).status == DesignStatus::ok);
        const auto late_singular = discretize(pole_at_2_over_Ts, 0.5, InputDelay{0.5}, DiscretizationMethod::Tustin);
        CHECK_FALSE(late_singular.has_value());
        CHECK(late_singular.status == DesignStatus::singular); // the plain step's cause

        //! Resample.
        CHECK(resample(sys, 0.1).status == DesignStatus::invalid_sample_time); // continuous
        CHECK(resample(fast, 0.0).status == DesignStatus::invalid_sample_time);
        CHECK(resample(fast, nan).status == DesignStatus::invalid_sample_time);
        const StateSpace<1, 1, 1> alternating{.A = {{-0.5}}, .B = {{1.0}}, .C = {{1.0}}, .Ts = 0.1};
        CHECK(resample(alternating, 0.15).status == DesignStatus::not_converged); // z < 0: no real logarithm
        CHECK(resample(alternating, 0.2).status == DesignStatus::ok);
    }

    TEST_CASE("value() of a success is a constant expression; a refusal is inspected through status") {
        //! value() of a failed result calls the non-constexpr trap, so only a success reaches it here.
        constexpr auto ok = discretize(make_lag(), 0.1);
        constexpr auto A = ok.value().A;
        static_assert(A(1, 1) > 0.0 && A(1, 1) < 1.0);
        constexpr auto refused = discretize(make_lag(), 0.1, InputDelay{0.2});
        static_assert(!refused && refused.status == DesignStatus::invalid_parameter);
    }
}

TEST_SUITE("Tustin pre-warping") {
    // Second-order plant with ω₀ = 50 rad/s sampled at 50 Hz: ω₀Ts = 1, so the plain
    // bilinear map evaluates the continuous response at (2/Ts)·tan(0.5) = 54.6 rad/s instead
    // of 50 and misses the resonance. Pre-warped at ω₀ the two responses agree there exactly.
    TEST_CASE("the discrete response matches the continuous one at the pre-warp frequency") {
        constexpr double    w0 = 50.0;
        constexpr double    zeta = 0.1;
        constexpr double    Ts = 0.02;
        StateSpace<2, 1, 1> sys{};
        sys.A(0, 1) = 1.0;
        sys.A(1, 0) = -w0 * w0;
        sys.A(1, 1) = -2.0 * zeta * w0;
        sys.B(1, 0) = w0 * w0;
        sys.C(0, 0) = 1.0;

        const auto plain = discretize(sys, Ts, DiscretizationMethod::Tustin);
        const auto warped = discretize(sys, Ts, DiscretizationMethod::Tustin, w0);
        REQUIRE(plain.has_value());
        REQUIRE(warped.has_value());
        CHECK(warped->Ts == Ts);

        const auto Gc = eval_frf(sys, damp::complex<double>{0.0, w0});
        const auto z = damp::complex<double>{std::cos(w0 * Ts), std::sin(w0 * Ts)};
        const auto Gp = eval_frf(*plain, z);
        const auto Gw = eval_frf(*warped, z);
        REQUIRE(Gc.has_value());
        REQUIRE(Gp.has_value());
        REQUIRE(Gw.has_value());
        const double ref = (*Gc)(0, 0).abs();
        CHECK(((*Gw)(0, 0) - (*Gc)(0, 0)).abs() <= 1e-9 * ref);
        CHECK(((*Gp)(0, 0) - (*Gc)(0, 0)).abs() > 1e-2 * ref); // the warping the option removes

        // Both maps keep the DC gain (z = 1 ↔ s = 0) and stability.
        const auto dc_w = eval_frf(*warped, damp::complex<double>{1.0, 0.0});
        REQUIRE(dc_w.has_value());
        CHECK((*dc_w)(0, 0).real() == doctest::Approx(1.0).epsilon(1e-9));
        CHECK(stability::is_stable(*warped));
    }

    TEST_CASE("a pre-warp frequency at or past Nyquist is refused") {
        StateSpace<1, 1, 1> sys{};
        sys.A(0, 0) = -1.0;
        sys.B(0, 0) = 1.0;
        sys.C(0, 0) = 1.0;
        const double Ts = 0.01;
        CHECK_FALSE(discretize(sys, Ts, DiscretizationMethod::Tustin, std::numbers::pi / Ts).has_value());
        CHECK(discretize(sys, Ts, DiscretizationMethod::Tustin, 0.9 * std::numbers::pi / Ts).has_value());
        CHECK(discretize(sys, Ts, DiscretizationMethod::ZOH, 1e9).has_value()); // ignored by other methods
    }
}

namespace {
// Damped oscillator (complex poles) with a process-noise input, continuous.
constexpr StateSpace<2, 1, 1, double, 1, 1> make_resample_plant() {
    return StateSpace<2, 1, 1, double, 1, 1>{
        .A = Matrix<2, 2>{{0.0, 1.0}, {-4.0, -0.8}},
        .B = Matrix<2, 1>{{0.0}, {2.0}},
        .C = Matrix<1, 2>{{1.0, 0.5}},
        .D = Matrix<1, 1>{{0.1}},
        .G = Matrix<2, 1>{{1.0}, {0.0}},
        .H = Matrix<1, 1>{{1.0}},
    };
}

template<size_t NX, size_t NU, size_t NY, size_t NW, size_t NV>
void check_same_model(const StateSpace<NX, NU, NY, double, NW, NV>& a, const StateSpace<NX, NU, NY, double, NW, NV>& b, double tol) {
    CHECK(a.Ts == doctest::Approx(b.Ts));
    for (size_t i = 0; i < NX; ++i) {
        for (size_t j = 0; j < NX; ++j) {
            CHECK(a.A(i, j) == doctest::Approx(b.A(i, j)).epsilon(tol).scale(1.0));
        }
        for (size_t j = 0; j < NU; ++j) {
            CHECK(a.B(i, j) == doctest::Approx(b.B(i, j)).epsilon(tol).scale(1.0));
        }
        for (size_t j = 0; j < NW; ++j) {
            CHECK(a.G(i, j) == doctest::Approx(b.G(i, j)).epsilon(tol).scale(1.0));
        }
    }
    CHECK(a.C == b.C);
    CHECK(a.D == b.D);
}
} // namespace

TEST_SUITE("Resampling (d2d / upsample)") {
    TEST_CASE("resample: a ZOH model at 0.1 s is the ZOH model at the new rate") {
        constexpr auto plant = make_resample_plant();
        constexpr auto fast = discretize(plant, 0.1);
        static_assert(fast.has_value());

        //! Integer ratio (matrix power) and fractional ratios (through the logarithm), both ways.
        for (const double Ts_new : {0.3, 0.25, 0.04, 0.1}) {
            const auto resampled = resample(*fast, Ts_new);
            REQUIRE(resampled.has_value());
            check_same_model(*resampled, *discretize(plant, Ts_new), 1e-9);
        }

        //! Round trip.
        const auto slow = resample(*fast, 0.25);
        REQUIRE(slow.has_value());
        const auto back = resample(*slow, 0.1);
        REQUIRE(back.has_value());
        check_same_model(*back, *fast, 1e-9);
    }

    TEST_CASE("resample is constexpr, in double and float") {
        constexpr auto fast = *discretize(make_resample_plant(), 0.1);
        constexpr auto slow = resample(fast, 0.25);
        static_assert(slow.has_value() && slow->Ts == 0.25);
        constexpr auto direct = *discretize(make_resample_plant(), 0.25);
        static_assert(damp::abs(slow->A(1, 0) - direct.A(1, 0)) < 1e-10);

        constexpr auto fast_f = fast.as<float>();
        constexpr auto slow_f = resample(fast_f, 0.25f);
        static_assert(slow_f.has_value());
        CHECK(slow_f->A(1, 0) == doctest::Approx(direct.A(1, 0)).epsilon(1e-4));
    }

    TEST_CASE_TEMPLATE("resample keeps a graded thermal model at a few eps", T, float, double) {
        // States in units 10⁴ apart; the logarithm of the ZOH block [[Ad, Bd], [0, I]] took its
        // square roots on the unbalanced matrix, one Schur form each: 2878 eps (float), 2·10⁴ (double).
        const StateSpace<3, 1, 3, T> sys{
            .A = Matrix<3, 3, T>{{T(-0.5), T(2e-3), T(0)}, {T(300), T(-8), T(1e-2)}, {T(0), T(5e3), T(-40)}},
            .B = Matrix<3, 1, T>{{T(10)}, {T(0)}, {T(0)}},
            .C = Matrix<3, 3, T>::identity(),
            .D = Matrix<3, 1, T>{},
        };
        const auto fast = discretize(sys, T(1e-4));
        const auto direct = discretize(sys, T(2.5e-4));
        REQUIRE(fast.has_value());
        REQUIRE(direct.has_value());
        const auto slow = resample(*fast, T(2.5e-4));
        REQUIRE(slow.has_value());
        Matrix<3, 4, T> got{};
        Matrix<3, 4, T> want{};
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                got(i, j) = slow->A(i, j);
                want(i, j) = direct->A(i, j);
            }
            got(i, 3) = slow->B(i, 0);
            want(i, 3) = direct->B(i, 0);
        }
        CHECK(Matrix<3, 4, T>(got - want).norm() <= T(32) * eps<T>() * want.norm());
    }

    TEST_CASE("resample refuses a continuous model, a bad rate, and a mode with no continuous parent") {
        constexpr auto plant = make_resample_plant();
        CHECK_FALSE(resample(plant, 0.1).has_value()); // continuous: use discretize
        const auto fast = *discretize(plant, 0.1);
        CHECK_FALSE(resample(fast, 0.0).has_value());
        CHECK_FALSE(resample(fast, -0.1).has_value());
        CHECK_FALSE(resample(fast, std::numeric_limits<double>::quiet_NaN()).has_value());

        //! A pole at z = −0.5 alternates sign every sample: no real logarithm, so no fractional
        //! resample; an integer ratio is still a plain power.
        const StateSpace<1, 1, 1, double> alternating{
            .A = Matrix<1, 1>{{-0.5}},
            .B = Matrix<1, 1>{{1.0}},
            .C = Matrix<1, 1>{{1.0}},
            .Ts = 0.1
        };
        CHECK_FALSE(resample(alternating, 0.15).has_value());
        const auto twice = resample(alternating, 0.2);
        REQUIRE(twice.has_value());
        CHECK(twice->A(0, 0) == doctest::Approx(0.25));
        CHECK(twice->B(0, 0) == doctest::Approx(0.5)); // A·B + B
    }

    TEST_CASE("upsample<L>: sys(zᴸ), the impulse response with L − 1 zeros between samples") {
        constexpr auto sys = *discretize(make_resample_plant(), 0.1);
        constexpr auto up = upsample<3>(sys);
        static_assert(up.Ts == sys.Ts / 3.0);
        static_assert(up.A(0, 4) == sys.A(0, 0) && up.A(2, 0) == 1.0 && up.C(0, 4) == sys.C(0, 0));

        ColVec<2> x{};
        ColVec<6> xu{};
        for (size_t k = 0; k < 30; ++k) {
            const double u = (k == 0) ? 1.0 : 0.0;
            const double y = (sys.C * x)(0, 0) + (sys.D(0, 0) * u);
            x = sys.A * x + sys.B * ColVec<1>{{u}};
            for (size_t r = 0; r < 3; ++r) {
                const double uu = (k == 0 && r == 0) ? 1.0 : 0.0;
                const double yu = (up.C * xu)(0, 0) + (up.D(0, 0) * uu);
                xu = up.A * xu + up.B * ColVec<1>{{uu}};
                CHECK(yu == doctest::Approx(r == 0 ? y : 0.0).epsilon(1e-12).scale(1.0));
            }
        }

        //! L = 1 is the model itself; the noise input rides with B.
        constexpr auto same = upsample<1>(sys);
        static_assert(same.A == sys.A && same.B == sys.B && same.G == sys.G && same.Ts == sys.Ts);
        static_assert(up.G(0, 0) == sys.G(0, 0) && up.G(2, 0) == 0.0);
    }

    TEST_CASE("upsample puts no noise on the chained copies when NW == NX·L") {
        constexpr StateSpace<1, 1, 1, double, 2> sys{.A = {{0.5}}, .B = {{1.0}}, .C = {{1.0}}, .G = {{0.3, 0.7}}, .Ts = 0.1};
        constexpr auto                           up = upsample<2>(sys);
        static_assert(up.G(0, 0) == 0.3 && up.G(0, 1) == 0.7);
        static_assert(up.G(1, 0) == 0.0 && up.G(1, 1) == 0.0);
    }
}

namespace {
// H(e^{jω}) of a z-domain TF, highest power first (Horner; oracle: std::complex).
template<size_t M>
std::complex<double> eval_z(const TransferFunction<M, M, double>& H, double w) {
    const std::complex<double> z = std::polar(1.0, w);
    std::complex<double>       num{};
    std::complex<double>       den{};
    for (size_t k = 0; k < M; ++k) {
        num = (num * z) + H.num[k];
        den = (den * z) + H.den[k];
    }
    return num / den;
}

// Exact DC group delay of B(z)/A(z): τ = Σp·aₚ/Σaₚ − Σp·bₚ/Σbₚ, entry k multiplying z^p with p = M − 1 − k.
template<size_t M>
double dc_group_delay(const TransferFunction<M, M, double>& H) {
    double sb = 0.0;
    double kb = 0.0;
    double sa = 0.0;
    double ka = 0.0;
    for (size_t k = 0; k < M; ++k) {
        const auto p = static_cast<double>(M - 1 - k);
        sb += H.num[k];
        kb += p * H.num[k];
        sa += H.den[k];
        ka += p * H.den[k];
    }
    return (ka / sa) - (kb / sb);
}
} // namespace

TEST_SUITE("Fractional delay (thiran)") {
    TEST_CASE("first and second order match Thiran's closed forms") {
        //! N = 1: H = (a + z⁻¹)/(1 + a z⁻¹) = (a z + 1)/(z + a), a = (1 − D)/(1 + D).
        constexpr auto H1 = fractional_delay<1>(0.05, 0.1); // D = 0.5
        static_assert(H1.has_value());
        CHECK(H1->num[0] == doctest::Approx(1.0 / 3.0));
        CHECK(H1->num[1] == doctest::Approx(1.0));
        CHECK(H1->den[0] == doctest::Approx(1.0));
        CHECK(H1->den[1] == doctest::Approx(1.0 / 3.0));

        //! N = 2: a₁ = −2(D − 2)/(D + 1), a₂ = (D − 1)(D − 2)/((D + 1)(D + 2)).
        const double D = 1.7;
        const auto   H2 = fractional_delay<2>(D * 1e-3, 1e-3);
        REQUIRE(H2.has_value());
        const double a1 = -2.0 * (D - 2.0) / (D + 1.0);
        const double a2 = (D - 1.0) * (D - 2.0) / ((D + 1.0) * (D + 2.0));
        //! den = {1, a₁, a₂} (MATLAB® thiran's den), num its reverse.
        CHECK(H2->den[0] == doctest::Approx(1.0));
        CHECK(H2->den[1] == doctest::Approx(a1).epsilon(1e-12));
        CHECK(H2->den[2] == doctest::Approx(a2).epsilon(1e-12));
        CHECK(H2->num[1] == doctest::Approx(a1).epsilon(1e-12));
        CHECK(H2->num[0] == doctest::Approx(a2).epsilon(1e-12));
    }

    TEST_CASE("all-pass with group delay D at DC, stable, for orders 1 to 5") {
        const double Ts = 0.01;
        const auto   check_order = [&]<size_t N>(double D) {
            const auto H = fractional_delay<N>(D * Ts, Ts);
            REQUIRE(H.has_value());
            CHECK(dc_group_delay(*H) == doctest::Approx(D).epsilon(1e-10));
            for (const double w : {0.01, 0.5, 1.5, 3.0}) {
                CHECK(std::abs(eval_z(*H, w)) == doctest::Approx(1.0).epsilon(1e-12));
            }
            //! Low-frequency phase ≈ −D·ω.
            CHECK(-std::arg(eval_z(*H, 1e-3)) / 1e-3 == doctest::Approx(D).epsilon(1e-6));
            //! Poles inside the unit circle: the companion realization is stable.
            const auto sys = H->to_state_space();
            REQUIRE(sys.has_value());
            auto sampled = *sys; // the TF carries no sample time; the realization is discrete
            sampled.Ts = 1.0;
            CHECK(stability::is_stable(sampled));
        };
        check_order.operator()<1>(0.3);
        check_order.operator()<2>(1.25);
        check_order.operator()<3>(2.6);
        check_order.operator()<4>(3.9);
        check_order.operator()<5>(4.5);
    }

    TEST_CASE("a whole delay is the pure delay; D ≤ N − 1 and bad inputs are refused") {
        constexpr auto H = fractional_delay<3>(0.3, 0.1); // D = 3 = N: 1/z³
        static_assert(H.has_value());
        CHECK(H->num[3] == doctest::Approx(1.0));
        CHECK(H->num[1] == doctest::Approx(0.0).scale(1.0));
        CHECK(H->num[0] == doctest::Approx(0.0).scale(1.0));
        CHECK(H->den[0] == doctest::Approx(1.0));
        CHECK(H->den[3] == doctest::Approx(0.0).scale(1.0));

        CHECK_FALSE(fractional_delay<3>(0.15, 0.1).has_value()); // D = 1.5 ≤ 2: unstable
        CHECK_FALSE(fractional_delay<1>(-0.2, 0.1).has_value()); // D = −2 ≤ 0
        CHECK_FALSE(fractional_delay<2>(0.15, 0.0).has_value());
        CHECK_FALSE(fractional_delay<2>(std::numeric_limits<double>::quiet_NaN(), 0.1).has_value());

        //! Float.
        constexpr auto Hf = fractional_delay<2>(0.0125f, 0.01f); // D = 1.25
        static_assert(Hf.has_value() && Hf->den[0] == 1.0f);
        CHECK(Hf->den[2] == doctest::Approx((0.25 * -0.75) / (2.25 * 3.25)).epsilon(1e-5));
    }
}
