// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <limits>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/fixed_gain_observer.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using namespace damp::mat;

namespace {

// A discrete plant (Ts = 0.1) with the given A and C, for the observer designs; B is unused by them.
template<size_t NX, size_t NY>
constexpr StateSpace<NX, 1, NY> observer_plant(const Matrix<NX, NX>& A, const Matrix<NY, NX>& C) {
    return StateSpace<NX, 1, NY>{.A = A, .C = C, .Ts = 0.1};
}

// The current form's error dynamics: e[k|k] = (I − M C) A e[k−1|k−1].
template<size_t NX, size_t NY>
Matrix<NX, NX> current_error(const Matrix<NX, NX>& A, const Matrix<NY, NX>& C, const Matrix<NX, NY>& M) {
    return (Matrix<NX, NX>::identity() - (M * C)) * A;
}

// Two real eigenvalues match {p, q} as a set.
bool real_pair_is(const damp::array<damp::complex<double>, 2>& ev, double p, double q) {
    const double a = ev[0].real();
    const double b = ev[1].real();
    const bool   imag_zero = damp::abs(ev[0].imag()) < 1e-9 && damp::abs(ev[1].imag()) < 1e-9;
    return imag_zero
        && ((a == doctest::Approx(p).epsilon(1e-9) && b == doctest::Approx(q).epsilon(1e-9))
            || (a == doctest::Approx(q).epsilon(1e-9) && b == doctest::Approx(p).epsilon(1e-9)));
}

} // namespace

TEST_SUITE("FixedGainObserver") {
    TEST_CASE("place_observer places error poles at requested (real) locations") {
        constexpr Matrix<2, 2> A{{1.0, 0.1}, {0.0, 1.0}};
        constexpr Matrix<1, 2> C{{1.0, 0.0}};

        const auto result = design::place_observer(observer_plant<2, 1>(A, C), {0.3, 0.4});
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

        // The current form runs M with A M = L: (I − M C) A = A⁻¹ (A − L C) A has the same poles.
        CHECK((A * result.M - result.L).norm() < 1e-12);
        const auto ev_c = compute_eigenvalues<2, double>(current_error(A, C, result.M)).values;
        CHECK(real_pair_is({ev_c[0], ev_c[1]}, 0.3, 0.4));
    }

    TEST_CASE("place_observer handles complex-conjugate pole pairs") {
        constexpr Matrix<2, 2> A{{1.0, 0.1}, {0.0, 1.0}};
        constexpr Matrix<1, 2> C{{1.0, 0.0}};

        const damp::array<damp::complex<double>, 2> poles{damp::complex{0.5, 0.2}, damp::complex{0.5, -0.2}};
        const auto                                  result = design::place_observer(observer_plant<2, 1>(A, C), poles);
        REQUIRE(result.success);

        // Gain must be real (conjugate pairs cancel imaginary parts).
        const Matrix<2, 2> A_err = A - result.L * C;
        const auto         ev = compute_eigenvalues<2, double>(A_err).values;
        // Both eigenvalues share real part 0.5 and |imag| 0.2.
        CHECK(ev[0].real() == doctest::Approx(0.5).epsilon(1e-9));
        CHECK(ev[1].real() == doctest::Approx(0.5).epsilon(1e-9));
        CHECK(std::abs(ev[0].imag()) == doctest::Approx(0.2).epsilon(1e-9));

        const auto ev_c = compute_eigenvalues<2, double>(current_error(A, C, result.M)).values;
        CHECK(ev_c[0].real() == doctest::Approx(0.5).epsilon(1e-9));
        CHECK(ev_c[1].real() == doctest::Approx(0.5).epsilon(1e-9));
        CHECK(std::abs(ev_c[0].imag()) == doctest::Approx(0.2).epsilon(1e-9));
    }

    TEST_CASE("multi-output observer places error poles (NY > 1)") {
        // 3 states, 2 measured outputs — the dual place() path that Ackermann could not do.
        constexpr Matrix<3, 3> A{{1.0, 0.1, 0.0}, {0.0, 1.0, 0.1}, {0.0, 0.0, 1.0}};
        constexpr Matrix<2, 3> C{{1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}};

        const auto result = design::place_observer(observer_plant<3, 2>(A, C), {0.2, 0.3, 0.4});
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

        // The current form: the same spectrum for (I − M C) A.
        const auto ev_c = compute_eigenvalues<3, double>(current_error(A, C, result.M)).values;
        double     sum_c = 0.0;
        double     prod_c = 1.0;
        for (size_t i = 0; i < 3; ++i) {
            CHECK(ev_c[i].imag() == doctest::Approx(0.0).epsilon(1e-9));
            sum_c += ev_c[i].real();
            prod_c *= ev_c[i].real();
        }
        CHECK(sum_c == doctest::Approx(0.9).epsilon(1e-9));
        CHECK(prod_c == doctest::Approx(0.024).epsilon(1e-9));
    }

    TEST_CASE("multi-output observer places a repeated pole beside a complex pair") {
        // Dual of a 4-state, 2-input placement: the complex KNV path refused the repeat at −2.
        constexpr Matrix<4, 4>    A{{0.0, 0.0, 0.0, 1.0}, {1.0, 0.0, 0.0, 2.0}, {0.0, 1.0, 0.0, 3.0}, {0.0, 0.0, 1.0, 4.0}};
        constexpr Matrix<2, 4>    C{{0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 1.0}};
        const StateSpace<4, 1, 2> sys{.A = A, .C = C};
        const auto                result = design::place_observer(sys, {damp::complex<double>{-2.0, 0.0}, damp::complex<double>{-2.0, 0.0}, damp::complex<double>{-1.0, 1.0}, damp::complex<double>{-1.0, -1.0}});
        REQUIRE(result.success);
        const auto ev = compute_eigenvalues<4, double>(A - (result.L * C)).values;
        int        at_minus_two = 0;
        for (size_t i = 0; i < 4; ++i) {
            if (damp::abs(ev[i] - damp::complex<double>{-2.0, 0.0}) < 1e-6) {
                ++at_minus_two;
            }
        }
        CHECK(at_minus_two == 2);
    }

    TEST_CASE("single-output observer places nearly repeated poles") {
        // The dual of a single-input cluster: the gain is unique, so the cluster at 0.5 is placed.
        constexpr Matrix<3, 3> A{{1.0, 0.1, 0.0}, {0.0, 1.0, 0.1}, {0.0, 0.0, 1.0}};
        constexpr Matrix<1, 3> C{{1.0, 0.0, 0.0}};
        const auto             result = design::place_observer(observer_plant<3, 1>(A, C), {0.5, 0.5 * (1.0 + 1e-9), 0.5 * (1.0 + 2e-9)});
        REQUIRE(result.success);
        const Matrix<3, 3> A_err = A - (result.L * C);
        CHECK(A_err.trace() == doctest::Approx(1.5).epsilon(1e-8)); // (z − 0.5)³
        CHECK(det(A_err) == doctest::Approx(0.125).epsilon(1e-8));
    }

    TEST_CASE("unobservable system reports failure") {
        // C sees only state 0; state 1 is decoupled and unobservable.
        constexpr Matrix<2, 2> A{{1.0, 0.0}, {0.0, 1.0}};
        constexpr Matrix<1, 2> C{{1.0, 0.0}};

        const auto result = design::place_observer(observer_plant<2, 1>(A, C), {0.3, 0.4});
        CHECK_FALSE(result.success);
        CHECK(result.status == design::DesignStatus::not_placeable);
    }

    TEST_CASE("result carries the plant and the achieved poles; the runtime deduces from it") {
        constexpr StateSpace<2, 1, 1> sys{
            .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
            .B = Matrix<2, 1>{{0.005}, {0.1}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .Ts = 0.1,
        };
        constexpr auto result = design::place_observer(sys, {0.3, 0.4});
        static_assert(result.success);
        static_assert(result.sys.Ts == 0.1);

        // e is eig(A − LC) from the gain, sorted here for the comparison.
        const double e0 = damp::min(result.e[0].real(), result.e[1].real());
        const double e1 = damp::max(result.e[0].real(), result.e[1].real());
        CHECK(e0 == doctest::Approx(0.3).epsilon(1e-9));
        CHECK(e1 == doctest::Approx(0.4).epsilon(1e-9));

        static constinit FixedGainObserver observer{result.as<float>()}; // no plant passed again
        static_assert(std::is_same_v<decltype(observer), FixedGainObserver<2, 1, 1, float>>);
        CHECK(observer.valid());
        CHECK(observer.gain()(1, 0) == doctest::Approx(static_cast<float>(result.M(1, 0)))); // the runtime holds M
        CHECK((sys.A * result.M - result.L).norm() < 1e-12);
    }

    TEST_CASE("continuous plant with a trailing Ts: ZOH model, z = e^{s Ts} poles") {
        constexpr StateSpace<2, 1, 1> plant{
            .A = Matrix<2, 2>{{0.0, 1.0}, {0.0, 0.0}},
            .B = Matrix<2, 1>{{0.0}, {1.0}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
        };
        constexpr double Ts = 0.01;
        const auto       sampled = design::place_observer(plant, {-20.0, -30.0}, Ts);
        REQUIRE(sampled.success);
        CHECK(sampled.sys.Ts == doctest::Approx(Ts));

        const auto plant_d = *discretize(plant, Ts, DiscretizationMethod::ZOH);
        const auto direct = design::place_observer(plant_d, {std::exp(-20.0 * Ts), std::exp(-30.0 * Ts)});
        REQUIRE(direct.success);
        CHECK(sampled.L(0, 0) == doctest::Approx(direct.L(0, 0)).epsilon(1e-9));
        CHECK(sampled.L(1, 0) == doctest::Approx(direct.L(1, 0)).epsilon(1e-9));
        CHECK(sampled.M(1, 0) == doctest::Approx(direct.M(1, 0)).epsilon(1e-9));

        // A continuous design has no current form: M = 0.
        const auto continuous = design::place_observer(plant, {-20.0, -30.0});
        REQUIRE(continuous.success);
        CHECK(continuous.M.norm() == 0.0);
        CHECK(continuous.L.norm() > 0.0);

        // A discrete plant is already sampled: a trailing Ts is refused, not applied twice.
        const auto twice = design::place_observer(plant_d, {-20.0, -30.0}, Ts);
        CHECK_FALSE(twice.success);
        CHECK(twice.status == design::DesignStatus::invalid_sample_time);
        CHECK(design::reduced_order_observer(plant_d, {-20.0}, Ts).status == design::DesignStatus::invalid_sample_time);
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

        const auto result = design::place_observer(sys, {0.2, 0.3});
        REQUIRE(result.success);

        FixedGainObserver<2, 1, 1, double> observer(result);
        observer.reset(ColVec<2>{0.0, 0.0});

        ColVec<2>       x{1.0, -0.5}; // unknown true state
        const ColVec<1> u{0.1};

        for (int k = 0; k < 60; ++k) {
            x = A * x + B * u; // the input over the last period moved the plant
            const ColVec<1> y = C * x;
            observer.step(y, u); // predict with u, then correct with the reading
        }

        const auto xhat = observer.state();
        CHECK(xhat(0, 0) == doctest::Approx(x(0, 0)).epsilon(1e-4));
        CHECK(xhat(1, 0) == doctest::Approx(x(1, 0)).epsilon(1e-4));
    }

    TEST_CASE("a singular A (a delay register): M placed on (A, C A), the register's mode pinned at z = 0") {
        // x⁺ = 0.9 x + 0.1 v, v⁺ = u: v is the command sent last period, so A has a zero row.
        constexpr Matrix<2, 2>    A{{0.9, 0.1}, {0.0, 0.0}};
        constexpr Matrix<2, 1>    B{{0.0}, {1.0}};
        constexpr Matrix<1, 2>    C{{1.0, 0.0}};
        const StateSpace<2, 1, 1> sys{.A = A, .B = B, .C = C, .Ts = 0.1};

        const auto result = design::place_observer(sys, {0.0, 0.3});
        REQUIRE(result.success);
        CHECK((A * result.M - result.L).norm() < 1e-12); // L = A M
        CHECK(real_pair_is(result.e, 0.0, 0.3));
        const auto ev_c = compute_eigenvalues<2, double>(current_error(A, C, result.M)).values;
        CHECK(real_pair_is({ev_c[0], ev_c[1]}, 0.0, 0.3));
        const auto ev_p = compute_eigenvalues<2, double>(A - result.L * C).values;
        CHECK(real_pair_is({ev_p[0], ev_p[1]}, 0.0, 0.3));

        // The estimate converges.
        FixedGainObserver observer{result};
        ColVec<2>         x{0.4, -0.2};
        for (int k = 0; k < 40; ++k) {
            const ColVec<1> u{{0.2 * damp::sin(0.3 * k)}};
            x = A * x + B * u;
            observer.step(C * x, u);
        }
        CHECK((observer.state() - x).norm() < 1e-12);

        // The register's mode is unobservable through C A: a pole anywhere else is not placeable.
        const auto pinned = design::place_observer(sys, {0.2, 0.3});
        CHECK_FALSE(pinned.success);
        CHECK(pinned.status == design::DesignStatus::not_placeable);
        CHECK(pinned.M.norm() == 0.0);
        CHECK(pinned.L.norm() == 0.0);
    }

    TEST_CASE("a well-conditioned A takes the solve A M = L: e, the runtime's poles, at the placed ones") {
        constexpr Matrix<2, 2> A{{1.0, 0.1}, {0.0, 1.0}};
        constexpr Matrix<1, 2> C{{1.0, 0.0}};
        CHECK(design::detail::solvable_for_innovation_gain(A));

        const auto result = design::place_observer(observer_plant<2, 1>(A, C), {0.3, 0.4});
        REQUIRE(result.success);
        const double kappa = one_norm(A) * one_norm(*solve(A, Matrix<2, 2>::identity()));
        const double tol = 256.0 * 2.0 * kappa * damp::eps<double>();
        const double lo = damp::min(result.e[0].real(), result.e[1].real());
        const double hi = damp::max(result.e[0].real(), result.e[1].real());
        CHECK(damp::abs(lo - 0.3) <= tol);
        CHECK(damp::abs(hi - 0.4) <= tol);
        CHECK(damp::abs(result.e[0].imag()) <= tol);
        CHECK(damp::abs(result.e[1].imag()) <= tol);
        CHECK((A * result.M - result.L).norm() <= tol * result.L.norm()); // L = A M to rounding: the solve
    }

    TEST_CASE("an ill-conditioned A (a mode much faster than the sample): M placed on (A, C A), never solved") {
        // κ₁(A) = 0.9 / 1e-12 = 9e11 > 1/√ε ≈ 6.7e7: A M = L would lose M's digits along the fast
        // mode, so M is placed directly on (A, C A).
        constexpr Matrix<2, 2>    A{{0.9, 0.0}, {0.0, 1e-12}};
        constexpr Matrix<2, 1>    B{{1.0}, {1.0}};
        constexpr Matrix<1, 2>    C{{1.0, 1.0}}; // sees both states
        const StateSpace<2, 1, 1> sys{.A = A, .B = B, .C = C, .Ts = 0.1};
        CHECK_FALSE(design::detail::solvable_for_innovation_gain(A));

        //! a pole near 0 for the fast mode: M bounded, e at the request
        const auto fast = design::place_observer(sys, {0.5, 0.0});
        REQUIRE(fast.success);
        CHECK((A * fast.M - fast.L).norm() == 0.0); // L = A M exactly: the direct placement
        CHECK(fast.M.norm() < 2.0 / C.norm());      // comparable to 1/‖C‖: M = [4/9, 5/9]
        CHECK(fast.M(0, 0) == doctest::Approx(4.0 / 9.0).epsilon(1e-12));
        CHECK(fast.M(1, 0) == doctest::Approx(5.0 / 9.0).epsilon(1e-12));
        CHECK(real_pair_is(fast.e, 0.5, 0.0));
        const auto ev_c = compute_eigenvalues<2, double>(current_error(A, C, fast.M)).values;
        CHECK(real_pair_is({ev_c[0], ev_c[1]}, 0.5, 0.0)); // e is the runtime's eig((I − M C) A)

        // The float runtime: the estimate converges with no amplified rounding along the fast mode.
        FixedGainObserver observer{fast.as<float>()};
        ColVec<2>         x{0.4, -0.2};
        double            scale = x.norm();
        for (int k = 0; k < 40; ++k) {
            const ColVec<1> u{{0.2 * damp::sin(0.3 * k)}};
            x = A * x + B * u;
            scale = damp::max(scale, x.norm());
            observer.step((C * x).as<float>(), u.as<float>());
        }
        CHECK((observer.state().as<double>() - x).norm() <= 64.0 * static_cast<double>(damp::eps<float>()) * scale);
    }

    TEST_CASE("a slow pole requested for a mode the plant forgets: the gain it demands, or not_placeable") {
        // A = diag(0.9, λ), C = [1 1]: det((I − M C) A) = 0.9 λ (1 − C M) must equal 0.5 · 0.4, so
        // C M ≈ −0.2 / (0.9 λ). No bounded gain moves the fast mode's current-form pole away from 0.
        // The KNV placement on (A, C A) starts refusing between λ = 1e-12 and 1e-13 (its rank test
        // on the eigenvectors is 64·NX·ε relative); λ = 1e-12 sits at that edge, so the two cases
        // below are taken clear of it on either side.
        constexpr Matrix<1, 2> C{{1.0, 1.0}};

        //! λ = 1e-10 (κ₁(A) = 9e9: the direct placement): the pair (A, C A) is still observable in
        //! double, so the request is honoured with a gain of order 1/λ.
        constexpr Matrix<2, 2> A{{0.9, 0.0}, {0.0, 1e-10}};
        CHECK_FALSE(design::detail::solvable_for_innovation_gain(A));
        const auto slow = design::place_observer(observer_plant<2, 1>(A, C), {0.5, 0.4});
        REQUIRE(slow.success);
        CHECK(slow.M.norm() > 0.1 / 1e-10);
        const double kappa = one_norm(A) * one_norm(*solve(A, Matrix<2, 2>::identity()));
        const double tol = 64.0 * kappa * damp::eps<double>(); // the placement's rounding on the nearly unobservable mode
        const double lo = damp::min(slow.e[0].real(), slow.e[1].real());
        const double hi = damp::max(slow.e[0].real(), slow.e[1].real());
        CHECK(damp::abs(lo - 0.4) <= tol);
        CHECK(damp::abs(hi - 0.5) <= tol);

        //! λ = e^{−62}, below ε against the slow mode: the fast mode is numerically unobservable
        //! through C A, so the placement refuses. A pole near 0 for it is placed with M bounded.
        const Matrix<2, 2> A62{{0.9, 0.0}, {0.0, damp::exp(-62.0)}};
        const auto         refused = design::place_observer(observer_plant<2, 1>(A62, C), {0.5, 0.4});
        CHECK_FALSE(refused.success);
        CHECK(refused.status == design::DesignStatus::not_placeable);
        CHECK(refused.M.norm() == 0.0);
        CHECK(refused.L.norm() == 0.0);
        const auto placed = design::place_observer(observer_plant<2, 1>(A62, C), {0.5, 0.0});
        REQUIRE(placed.success);
        CHECK(placed.M.norm() < 2.0 / C.norm());
        CHECK(real_pair_is(placed.e, 0.5, 0.0));
    }

    TEST_CASE("the current form on design::kalman's gain: error poles of (I − M C) A, and the estimate converges") {
        StateSpace<2, 1, 1, double, 2, 1> sys{};
        sys.A = {{1.0, 0.1}, {0.0, 1.0}};
        sys.B = {{0.005}, {0.1}};
        sys.C = {{1.0, 0.0}};
        sys.Ts = 0.1;
        const auto kf = design::kalman(sys, {.process = Matrix<2, 2>::identity() * 1e-3, .measurement = 1e-2});
        REQUIRE(kf.success);

        FixedGainObserver obs{kf};
        static_assert(std::is_same_v<decltype(obs), FixedGainObserver<2, 1, 1, double>>);
        CHECK((obs.gain() - kf.M).norm() == 0.0);
        const auto ev = compute_eigenvalues<2, double>(current_error(sys.A, sys.C, kf.M)).values;
        CHECK(ev[0].abs() < 1.0);
        CHECK(ev[1].abs() < 1.0);

        ColVec<2> x{1.0, -0.5};
        for (int k = 0; k < 400; ++k) {
            const ColVec<1> u{{0.1}};
            x = sys.A * x + sys.B * u;
            obs.step(sys.C * x, u);
        }
        CHECK((obs.state() - x).norm() < 1e-9);
    }

    TEST_CASE("from the model and an innovation gain; a continuous design is inert at run time") {
        constexpr StateSpace<2, 1, 1> sys{
            .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
            .B = Matrix<2, 1>{{0.005}, {0.1}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .Ts = 0.1,
        };
        constexpr auto result = design::place_observer(sys, {0.3, 0.4});
        static_assert(result.success);

        //! the matrices: the same observer as the design's
        FixedGainObserver from_design{result};
        FixedGainObserver from_matrices{sys, result.M};
        static_assert(std::is_same_v<decltype(from_matrices), FixedGainObserver<2, 1, 1, double>>);
        CHECK(from_matrices.valid());
        for (int k = 0; k < 10; ++k) {
            const ColVec<1> y{{0.1 * k}};
            const ColVec<1> u{{0.5}};
            CHECK((from_design.step(y, u) - from_matrices.step(y, u)).norm() == 0.0);
        }

        //! a continuous design or model is not a sampled observer: L = 0, valid() false
        const auto continuous = design::place_observer(StateSpace<2, 1, 1>{.A = {{0.0, 1.0}, {0.0, 0.0}}, .B = {{0.0}, {1.0}}, .C = {{1.0, 0.0}}}, {-20.0, -30.0});
        REQUIRE(continuous.success);
        FixedGainObserver inert{continuous};
        CHECK_FALSE(inert.valid());
        CHECK(inert.gain().norm() == 0.0);
        FixedGainObserver inert_m{StateSpace<2, 1, 1>{.A = sys.A, .B = sys.B, .C = sys.C}, result.M};
        CHECK_FALSE(inert_m.valid());

        //! the same observer in another scalar
        const FixedGainObserver<2, 1, 1, float> as_float{from_design};
        CHECK(as_float.state()(0) == doctest::Approx(static_cast<float>(from_design.state()(0))));
    }

    TEST_CASE("design is constexpr-evaluable") {
        constexpr auto result = [] {
            constexpr Matrix<2, 2> A{{1.0, 0.1}, {0.0, 1.0}};
            constexpr Matrix<1, 2> C{{1.0, 0.0}};
            return design::place_observer(observer_plant<2, 1>(A, C), {0.3, 0.4});
        }();
        static_assert(result.success, "observer design must evaluate at compile time");

        constexpr auto result_f = result.as<float>();
        static_assert(result_f.success);
        CHECK(result_f.success);
    }

    TEST_CASE("set_state clamps the estimate between steps") {
        constexpr Matrix<2, 2>    A{{1.0, 0.1}, {0.0, 1.0}};
        constexpr Matrix<2, 1>    B{{0.0}, {0.1}};
        constexpr Matrix<1, 2>    C{{1.0, 0.0}};
        constexpr Matrix<1, 1>    D{{0.0}};
        const StateSpace<2, 1, 1> sys{.A = A, .B = B, .C = C, .D = D, .Ts = 0.1};
        const auto                result = design::place_observer(sys, {0.3, 0.4});
        REQUIRE(result.success);

        FixedGainObserver<2, 1, 1, double> obs(result);

        obs.step(ColVec<1>{1.0}, {});
        // Force a non-physical estimate, then clamp it.
        obs.set_state(ColVec<2>{99.0, -99.0});
        CHECK(obs.state()(0, 0) == doctest::Approx(99.0));

        obs.set_state(1, 0.0); // clamp the velocity component
        CHECK(obs.state()(1, 0) == doctest::Approx(0.0));
    }
}

TEST_SUITE("Reduced-Order FixedGainObserver") {
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
        const auto result = design::reduced_order_observer(sys, {0.5});
        REQUIRE(result.success);
        CHECK(result.F(0, 0) == doctest::Approx(0.5).epsilon(1e-9)); // scalar error pole

        ReducedOrderObserver<2, 1, double> observer(result);
        CHECK(observer.valid());
        observer.reset();

        // A continuous design (Ts = 0) is not a sampled observer: at run time it is inert.
        const auto continuous = design::reduced_order_observer(StateSpace<2, 1, 1>{.A = {{0.0, 1.0}, {0.0, 0.0}}, .B = {{0.0}, {1.0}}, .C = C}, {-20.0});
        REQUIRE(continuous.success);
        ReducedOrderObserver<2, 1, double> inert(continuous);
        CHECK_FALSE(inert.valid());
        (void)inert.step(ColVec<1>{1.0}, ColVec<1>{0.5});
        CHECK(inert.state()(0, 0) == 0.0);
        CHECK(inert.internal_state()(0, 0) == 0.0);

        ColVec<2>       x{1.0, 2.0}; // true: position 1, velocity 2 (velocity unknown)
        const ColVec<1> u{0.1};
        for (int k = 0; k < 50; ++k) {
            x = A * x + B * u;
            (void)observer.step(C * x, u); // u: the input applied over the period that just ended
        }

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

        const auto result = design::reduced_order_observer(sys, {0.4, 0.5});
        REQUIRE(result.success);

        // Error dynamics F must have eigenvalues at {0.4, 0.5}.
        const auto   ev = compute_eigenvalues<2, double>(result.F).values;
        const double a = ev[0].real();
        const double b = ev[1].real();
        const bool   matches = (a == doctest::Approx(0.4).epsilon(1e-9) && b == doctest::Approx(0.5).epsilon(1e-9))
                            || (a == doctest::Approx(0.5).epsilon(1e-9) && b == doctest::Approx(0.4).epsilon(1e-9));
        CHECK(matches);

        ReducedOrderObserver<3, 1, double> observer(result);
        observer.reset();

        ColVec<3>       x{0.5, -1.0, 2.0};
        const ColVec<1> u{0.2};
        for (int k = 0; k < 80; ++k) {
            x = A * x + B * u;
            (void)observer.step(C * x, u);
        }

        const auto xhat = observer.state();
        CHECK(xhat(0, 0) == doctest::Approx(x(0, 0)).epsilon(1e-9));
        CHECK(xhat(1, 0) == doctest::Approx(x(1, 0)).epsilon(1e-3));
        CHECK(xhat(2, 0) == doctest::Approx(x(2, 0)).epsilon(1e-3));
    }

    TEST_CASE("step(y, u) takes the input of the period that just ended, as every estimator") {
        //! A time-varying input: the input of the last period drives z, and step returns the
        //! estimate, so u[k] = −K x̂[k] can follow the call.
        constexpr double                   Ts = 0.01;
        constexpr Matrix<2, 2>             A{{1.0, Ts}, {0.0, 1.0}};
        constexpr Matrix<2, 1>             B{{0.5 * Ts * Ts}, {Ts}};
        constexpr Matrix<1, 2>             C{{1.0, 0.0}};
        const StateSpace<2, 1, 1>          sys{.A = A, .B = B, .C = C, .Ts = Ts};
        ReducedOrderObserver<2, 1, double> observer(design::reduced_order_observer(sys, {0.5}));
        REQUIRE(observer.valid());

        ColVec<2> x{1.0, 2.0};
        double    worst = 0.0;
        for (int k = 0; k < 400; ++k) {
            const ColVec<1> u{10.0 * damp::sin(0.05 * static_cast<double>(k))};
            x = A * x + B * u;
            const ColVec<2> xhat = observer.step(C * x, u);
            if (k >= 100) {
                worst = damp::max(worst, damp::abs(xhat(1, 0) - x(1, 0)));
            }
        }
        CHECK(worst < 1e-9); // was 7.5e-3: z advanced with the next period's input
    }

    TEST_CASE("refused: a feedthrough, a non-finite model entry") {
        //! The measured state is read as y: with D ≠ 0 the estimate would be off by D u.
        constexpr Matrix<2, 2> A{{1.0, 0.1}, {0.0, 1.0}};
        constexpr Matrix<2, 1> B{{0.005}, {0.1}};
        constexpr Matrix<1, 2> C{{1.0, 0.0}};
        const auto             with_d = design::reduced_order_observer(StateSpace<2, 1, 1>{.A = A, .B = B, .C = C, .D = {{0.5}}, .Ts = 0.1}, {0.5});
        CHECK_FALSE(with_d.success);
        CHECK(with_d.status == design::DesignStatus::bad_weights);

        //! Placement reads only (A, C); a NaN in B or D would ship to the runtime.
        StateSpace<2, 1, 1> nan_b{.A = A, .B = B, .C = C, .Ts = 0.1};
        nan_b.B(1, 0) = std::numeric_limits<double>::quiet_NaN();
        CHECK(design::reduced_order_observer(nan_b, {0.5}).status == design::DesignStatus::bad_weights);
        CHECK(design::place_observer(nan_b, {0.3, 0.4}).status == design::DesignStatus::bad_weights);
        StateSpace<2, 1, 1> nan_d{.A = A, .B = B, .C = C, .Ts = 0.1};
        nan_d.D(0, 0) = std::numeric_limits<double>::quiet_NaN();
        CHECK(design::place_observer(nan_d, {0.3, 0.4}).status == design::DesignStatus::bad_weights);

        //! A refused design's poles are zero-filled, inside the unit circle: not stable.
        const auto unobservable = design::place_observer(observer_plant<2, 1>(Matrix<2, 2>{{0.5, 0.0}, {0.0, 1.5}}, C), {0.1, 0.2});
        REQUIRE_FALSE(unobservable.success);
        CHECK_FALSE(unobservable.is_stable());
        const auto unobservable_reduced = design::reduced_order_observer(observer_plant<2, 1>(Matrix<2, 2>{{0.5, 0.0}, {0.0, 1.5}}, C), {0.1});
        REQUIRE_FALSE(unobservable_reduced.success);
        CHECK_FALSE(unobservable_reduced.is_stable());
    }

    TEST_CASE("reduced-order design is constexpr-evaluable") {
        constexpr auto result = [] {
            constexpr Matrix<2, 2> A{{1.0, 0.01}, {0.0, 1.0}};
            constexpr Matrix<2, 1> B{{0.0}, {0.01}};
            constexpr Matrix<1, 2> C{{1.0, 0.0}};
            return design::reduced_order_observer(StateSpace<2, 1, 1>{.A = A, .B = B, .C = C, .Ts = 0.1}, {0.5});
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

        const auto result = design::reduced_order_observer(sys, {0.5});
        REQUIRE(result.success);

        ReducedOrderObserver<2, 1, double> obs(result);
        (void)obs.step(ColVec<1>{1.0}, {});

        // The persistent state is z (NX-1 = 1 element here), not the
        // reconstructed full state. Writing z and stepping with a steady
        // measurement reconstructs x = Tinv*[y; z + L*y].
        obs.set_internal_state(ColVec<1>{0.0});
        CHECK(obs.internal_state()(0, 0) == doctest::Approx(0.0));

        obs.set_internal_state(0, 3.0);
        CHECK(obs.internal_state()(0, 0) == doctest::Approx(3.0));
    }
}
