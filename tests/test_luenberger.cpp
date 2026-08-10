// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>

#include "damp/backend.hpp"
#include "damp/estimation/luenberger.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using namespace damp::mat;

TEST_SUITE("Luenberger") {
    TEST_CASE("luenberger places error poles at requested (real) locations") {
        constexpr Matrix<2, 2> A{{1.0, 0.1}, {0.0, 1.0}};
        constexpr Matrix<1, 2> C{{1.0, 0.0}};

        const auto result = design::luenberger<2, 1>(A, C, ColVec<2>{0.3, 0.4});
        REQUIRE(result.success);

        // Error dynamics matrix A - L*C must have eigenvalues at {0.3, 0.4}.
        const Matrix<2, 2> A_err = A - result.L * C;
        const auto         ev = compute_eigenvalues<2, double>(A_err).values;

        const double a = ev[0].real();
        const double b = ev[1].real();
        const bool   matches = (a == doctest::Approx(0.3).epsilon(1e-9) && b == doctest::Approx(0.4).epsilon(1e-9))
                          || (a == doctest::Approx(0.4).epsilon(1e-9) && b == doctest::Approx(0.3).epsilon(1e-9));
        CHECK(matches);
        CHECK(ev[0].imag() == doctest::Approx(0.0).epsilon(1e-9));
        CHECK(ev[1].imag() == doctest::Approx(0.0).epsilon(1e-9));
        CHECK(result.is_stable());
    }

    TEST_CASE("luenberger handles complex-conjugate pole pairs") {
        constexpr Matrix<2, 2> A{{1.0, 0.1}, {0.0, 1.0}};
        constexpr Matrix<1, 2> C{{1.0, 0.0}};

        const ColVec<2, damp::complex<double>> poles{
            damp::complex<double>{0.5, 0.2},
            damp::complex<double>{0.5, -0.2}
        };
        const auto result = design::luenberger<2, 1>(A, C, poles);
        REQUIRE(result.success);

        // Gain must be real (conjugate pairs cancel imaginary parts).
        const Matrix<2, 2> A_err = A - result.L * C;
        const auto         ev = compute_eigenvalues<2, double>(A_err).values;
        // Both eigenvalues share real part 0.5 and |imag| 0.2.
        CHECK(ev[0].real() == doctest::Approx(0.5).epsilon(1e-9));
        CHECK(ev[1].real() == doctest::Approx(0.5).epsilon(1e-9));
        CHECK(std::abs(ev[0].imag()) == doctest::Approx(0.2).epsilon(1e-9));
    }

    TEST_CASE("multi-output observer places error poles (NY > 1)") {
        // 3 states, 2 measured outputs — the dual place() path that Ackermann could not do.
        constexpr Matrix<3, 3> A{{1.0, 0.1, 0.0}, {0.0, 1.0, 0.1}, {0.0, 0.0, 1.0}};
        constexpr Matrix<2, 3> C{{1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}};

        const auto result = design::luenberger<3, 2>(A, C, ColVec<3>{0.2, 0.3, 0.4});
        REQUIRE(result.success);

        const Matrix<3, 3> A_err = A - result.L * C;
        const auto         ev = compute_eigenvalues<3, double>(A_err).values;
        // Match the requested spectrum {0.2, 0.3, 0.4} as a set.
        double sum = 0.0;
        double prod = 1.0;
        for (size_t i = 0; i < 3; ++i) {
            CHECK(ev[i].imag() == doctest::Approx(0.0).epsilon(1e-9));
            sum += ev[i].real();
            prod *= ev[i].real();
        }
        CHECK(sum == doctest::Approx(0.9).epsilon(1e-9));    // 0.2 + 0.3 + 0.4
        CHECK(prod == doctest::Approx(0.024).epsilon(1e-9)); // 0.2 * 0.3 * 0.4
        CHECK(result.is_stable());
    }

    TEST_CASE("unobservable system reports failure") {
        // C sees only state 0; state 1 is decoupled and unobservable.
        constexpr Matrix<2, 2> A{{1.0, 0.0}, {0.0, 1.0}};
        constexpr Matrix<1, 2> C{{1.0, 0.0}};

        const auto result = design::luenberger<2, 1>(A, C, ColVec<2>{0.3, 0.4});
        CHECK_FALSE(result.success);
    }

    TEST_CASE("observer estimate converges to the true state") {
        constexpr Matrix<2, 2>    A{{1.0, 0.1}, {0.0, 1.0}};
        constexpr Matrix<2, 1>    B{{0.005}, {0.1}};
        constexpr Matrix<1, 2>    C{{1.0, 0.0}};
        const StateSpace<2, 1, 1> sys{
            .A = A,
            .B = B,
            .C = C,
            .D = Matrix<1, 1>::zeros(),
            .Ts = 0.1
        };

        const auto result = design::luenberger(sys, ColVec<2>{0.2, 0.3});
        REQUIRE(result.success);

        Luenberger<2, 1, 1, double> observer(sys, result);
        observer.reset(ColVec<2>{0.0, 0.0});

        ColVec<2>       x{1.0, -0.5}; // unknown true state
        const ColVec<1> u{0.1};

        for (int k = 0; k < 60; ++k) {
            const ColVec<1> y = C * x;
            observer.step(y, u);
            x = A * x + B * u;
        }

        const auto xhat = observer.state();
        CHECK(xhat(0, 0) == doctest::Approx(x(0, 0)).epsilon(1e-4));
        CHECK(xhat(1, 0) == doctest::Approx(x(1, 0)).epsilon(1e-4));
    }

    TEST_CASE("design is constexpr-evaluable") {
        constexpr auto result = [] {
            constexpr Matrix<2, 2> A{{1.0, 0.1}, {0.0, 1.0}};
            constexpr Matrix<1, 2> C{{1.0, 0.0}};
            return design::luenberger<2, 1>(A, C, ColVec<2>{0.3, 0.4});
        }();
        static_assert(result.success, "observer design must evaluate at compile time");

        constexpr auto result_f = result.as<float>();
        static_assert(result_f.success);
        CHECK(result_f.success);
    }

    TEST_CASE("set_state clamps the estimate between steps") {
        constexpr Matrix<2, 2> A{{1.0, 0.1}, {0.0, 1.0}};
        constexpr Matrix<2, 1> B{{0.0}, {0.1}};
        constexpr Matrix<1, 2> C{{1.0, 0.0}};
        constexpr Matrix<1, 1> D{{0.0}};
        const auto             result = design::luenberger<2, 1>(A, C, ColVec<2>{0.3, 0.4});
        REQUIRE(result.success);

        const StateSpace<2, 1, 1>   sys{.A = A, .B = B, .C = C, .D = D};
        Luenberger<2, 1, 1, double> obs(sys, result);

        obs.step(ColVec<1>{1.0});
        // Force a non-physical estimate, then clamp it.
        obs.set_state(ColVec<2>{99.0, -99.0});
        CHECK(obs.state()(0, 0) == doctest::Approx(99.0));

        obs.set_state(1, 0.0); // clamp the velocity component
        CHECK(obs.state()(1, 0) == doctest::Approx(0.0));
    }
}

TEST_SUITE("Reduced-Order Luenberger") {
    TEST_CASE("encoder case: estimate velocity from measured position") {
        constexpr double          Ts = 0.01;
        constexpr Matrix<2, 2>    A{{1.0, Ts}, {0.0, 1.0}}; // [position; velocity]
        constexpr Matrix<2, 1>    B{{0.0}, {Ts}};           // acceleration input
        constexpr Matrix<1, 2>    C{{1.0, 0.0}};            // measure position (encoder)
        const StateSpace<2, 1, 1> sys{
            .A = A,
            .B = B,
            .C = C,
            .D = Matrix<1, 1>::zeros(),
            .Ts = Ts
        };

        // One unmeasured state (velocity) -> one error pole.
        const auto result = design::reduced_luenberger(sys, ColVec<1>{0.5});
        REQUIRE(result.success);
        CHECK(result.F(0, 0) == doctest::Approx(0.5).epsilon(1e-9)); // scalar error pole

        ReducedLuenberger<2, 1, double> observer(result);
        observer.reset();

        ColVec<2>       x{1.0, 2.0}; // true: position 1, velocity 2 (velocity unknown)
        const ColVec<1> u{0.1};
        for (int k = 0; k < 50; ++k) {
            observer.step(C * x, u);
            x = A * x + B * u;
        }
        observer.step(C * x, u); // reconstruct estimate for the current state

        const auto xhat = observer.state();
        CHECK(xhat(0, 0) == doctest::Approx(x(0, 0)).epsilon(1e-9)); // position read exactly from y
        CHECK(xhat(1, 0) == doctest::Approx(x(1, 0)).epsilon(1e-3)); // velocity converged
    }

    TEST_CASE("two unmeasured states: error poles placed and estimate converges") {
        constexpr double          Ts = 0.01;
        constexpr Matrix<3, 3>    A{{1.0, Ts, 0.0}, {0.0, 1.0, Ts}, {0.0, 0.0, 1.0}};
        constexpr Matrix<3, 1>    B{{0.0}, {0.0}, {Ts}};
        constexpr Matrix<1, 3>    C{{1.0, 0.0, 0.0}};
        const StateSpace<3, 1, 1> sys{
            .A = A,
            .B = B,
            .C = C,
            .D = Matrix<1, 1>::zeros(),
            .Ts = Ts
        };

        const auto result = design::reduced_luenberger(sys, ColVec<2>{0.4, 0.5});
        REQUIRE(result.success);

        // Error dynamics F must have eigenvalues at {0.4, 0.5}.
        const auto   ev = compute_eigenvalues<2, double>(result.F).values;
        const double a = ev[0].real();
        const double b = ev[1].real();
        const bool   matches = (a == doctest::Approx(0.4).epsilon(1e-9) && b == doctest::Approx(0.5).epsilon(1e-9))
                          || (a == doctest::Approx(0.5).epsilon(1e-9) && b == doctest::Approx(0.4).epsilon(1e-9));
        CHECK(matches);

        ReducedLuenberger<3, 1, double> observer(result);
        observer.reset();

        ColVec<3>       x{0.5, -1.0, 2.0};
        const ColVec<1> u{0.2};
        for (int k = 0; k < 80; ++k) {
            observer.step(C * x, u);
            x = A * x + B * u;
        }
        observer.step(C * x, u);

        const auto xhat = observer.state();
        CHECK(xhat(0, 0) == doctest::Approx(x(0, 0)).epsilon(1e-9));
        CHECK(xhat(1, 0) == doctest::Approx(x(1, 0)).epsilon(1e-3));
        CHECK(xhat(2, 0) == doctest::Approx(x(2, 0)).epsilon(1e-3));
    }

    TEST_CASE("reduced-order design is constexpr-evaluable") {
        constexpr auto result = [] {
            constexpr Matrix<2, 2> A{{1.0, 0.01}, {0.0, 1.0}};
            constexpr Matrix<2, 1> B{{0.0}, {0.01}};
            constexpr Matrix<1, 2> C{{1.0, 0.0}};
            return design::reduced_luenberger<2, 1, 1>(A, B, C, ColVec<1>{0.5});
        }();
        static_assert(result.success, "reduced-order observer must evaluate at compile time");

        constexpr auto result_f = result.as<float>();
        static_assert(result_f.success);
        CHECK(result_f.success);
    }

    TEST_CASE("set_internal_state overwrites the persistent recursion state") {
        constexpr double          Ts = 0.01;
        constexpr Matrix<2, 2>    A{{1.0, Ts}, {0.0, 1.0}};
        constexpr Matrix<2, 1>    B{{0.0}, {Ts}};
        constexpr Matrix<1, 2>    C{{1.0, 0.0}};
        const StateSpace<2, 1, 1> sys{.A = A, .B = B, .C = C, .D = Matrix<1, 1>::zeros(), .Ts = Ts};

        const auto result = design::reduced_luenberger(sys, ColVec<1>{0.5});
        REQUIRE(result.success);

        ReducedLuenberger<2, 1, double> obs(result);
        obs.step(ColVec<1>{1.0});

        // The persistent state is z (NX-1 = 1 element here), not the
        // reconstructed full state. Writing z and stepping with a steady
        // measurement reconstructs x = Tinv*[y; z + L*y].
        obs.set_internal_state(ColVec<1>{0.0});
        CHECK(obs.internal_state()(0, 0) == doctest::Approx(0.0));

        obs.set_internal_state(0, 3.0);
        CHECK(obs.internal_state()(0, 0) == doctest::Approx(3.0));
    }
}
