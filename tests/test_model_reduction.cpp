// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/design/model_reduction.hpp"
#include "damp/design/stability.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using namespace damp::design;

namespace {

/// Max-abs entry-wise difference.
template<size_t R, size_t C>
double max_abs_diff(const Matrix<R, C>& A, const Matrix<R, C>& B) {
    double m = 0.0;
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            m = damp::max(m, std::abs(A(i, j) - B(i, j)));
        }
    }
    return m;
}

/// Continuous DC gain G(0) = D − C A^{-1} B via solve.
template<size_t NX, size_t NU, size_t NY>
Matrix<NY, NU> cont_dc(const StateSpace<NX, NU, NY>& sys) {
    const auto X = mat::solve(-sys.A, sys.B);
    REQUIRE(X.has_value());
    return sys.C * (*X) + sys.D;
}

/// Frobenius-ish residual of continuous Lyapunov A X + X Aᵀ + Q.
template<size_t N>
double lyap_resid(const Matrix<N, N>& A, const Matrix<N, N>& X, const Matrix<N, N>& Q) {
    const Matrix<N, N> R = A * X + X * A.transpose() + Q;
    double             m = 0.0;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            m = damp::max(m, std::abs(R(i, j)));
        }
    }
    return m;
}

/// Two-time-scale SISO plant: strong slow mode + weak fast mode.
StateSpace<2, 1, 1> two_scale_plant() {
    // ẋ = [−1 0; 0 −20] x + [1; 0.05] u,  y = [1 0.05] x
    // Dominant input-output path is the slow pole; fast mode is weakly coupled.
    return StateSpace<2, 1, 1>{
        .A = Matrix<2, 2>{{-1.0, 0.0}, {0.0, -20.0}},
        .B = Matrix<2, 1>{{1.0}, {0.05}},
        .C = Matrix<1, 2>{{1.0, 0.05}},
        .D = Matrix<1, 1>{{0.0}},
        .Ts = 0.0
    };
}

/// DC motor, x = [i, ω], u = V, y = [i, ω]: R = 1 Ω, L = 1 mH, Kt = Ke = 0.1, J = 1e-3, b = 1e-4.
constexpr StateSpace<2, 1, 2> dc_motor{
    .A = Matrix<2, 2>{{-1000.0, -100.0}, {100.0, -0.1}},
    .B = Matrix<2, 1>{{1000.0}, {0.0}},
    .C = Matrix<2, 2>{{1.0, 0.0}, {0.0, 1.0}},
};

/// Diagonal plant with poles −1, −100, −2 and unit B, C.
constexpr StateSpace<3, 1, 1> three_poles{
    .A = Matrix<3, 3>{{-1.0, 0.0, 0.0}, {0.0, -100.0, 0.0}, {0.0, 0.0, -2.0}},
    .B = Matrix<3, 1>{{1.0}, {1.0}, {1.0}},
    .C = Matrix<1, 3>{{1.0, 1.0, 1.0}},
};

/// True when eliminate_states accepts the index list for an NX-state model.
template<size_t NX, size_t... Elim>
constexpr bool can_eliminate = requires(const StateSpace<NX, 1, 1>& sys) { eliminate_states<Elim...>(sys); };

static_assert(can_eliminate<2, 0>);
static_assert(can_eliminate<3, 2, 0>);
static_assert(!can_eliminate<2>);       // nothing eliminated
static_assert(!can_eliminate<2, 0, 1>); // every state eliminated
static_assert(!can_eliminate<2, 2>);    // index out of range
static_assert(!can_eliminate<3, 1, 1>); // duplicate index

/// Stable fully coupled 2nd-order plant (both modes well excited).
StateSpace<2, 1, 1> coupled_2nd_order() {
    return StateSpace<2, 1, 1>{
        .A = Matrix<2, 2>{{-2.0, 1.0}, {0.0, -3.0}},
        .B = Matrix<2, 1>{{1.0}, {1.0}},
        .C = Matrix<1, 2>{{1.0, 0.0}},
        .D = Matrix<1, 1>{{0.0}},
        .Ts = 0.0
    };
}

} // namespace

TEST_SUITE("model reduction") {

    TEST_CASE("balanced_realization: Hankel SV descending and gramians balanced") {
        const auto sys = coupled_2nd_order();
        const auto bal = balanced_realization(sys);
        REQUIRE(bal.success);

        CHECK(bal.sigma(0) >= bal.sigma(1));
        CHECK(bal.sigma(0) > 0.0);
        CHECK(bal.sigma(1) > 0.0);

        // In balanced coords Wc ≈ Wo ≈ diag(σ)
        const auto Wc = stability::controllability_gramian(bal.balanced.A, bal.balanced.B, false);
        const auto Wo = stability::observability_gramian(bal.balanced.A, bal.balanced.C, false);
        REQUIRE(Wc.has_value());
        REQUIRE(Wo.has_value());

        Matrix<2, 2> Sigma{};
        Sigma(0, 0) = bal.sigma(0);
        Sigma(1, 1) = bal.sigma(1);

        CHECK(max_abs_diff(*Wc, Sigma) < 1e-8);
        CHECK(max_abs_diff(*Wo, Sigma) < 1e-8);
        CHECK(max_abs_diff(*Wc, *Wo) < 1e-8);

        // Lyapunov residual still holds for the balanced system
        CHECK(lyap_resid(bal.balanced.A, *Wc, bal.balanced.B * bal.balanced.B.transpose()) < 1e-8);
    }

    TEST_CASE("balanced_realization: preserves transfer function (DC + poles)") {
        const auto sys = coupled_2nd_order();
        const auto bal = balanced_realization(sys);
        REQUIRE(bal.success);

        // DC gain invariant under similarity
        const auto dc0 = cont_dc(sys);
        const auto dc1 = cont_dc(bal.balanced);
        CHECK(std::abs(dc0(0, 0) - dc1(0, 0)) < 1e-10);

        // Eigenvalues of A preserved (similarity)
        const auto e0 = mat::compute_eigenvalues(sys.A);
        const auto e1 = mat::compute_eigenvalues(bal.balanced.A);
        REQUIRE(e0.success);
        REQUIRE(e1.success);

        // Match poles as a multiset (sort by real part)
        damp::array<double, 2> r0{e0.values[0].real(), e0.values[1].real()};
        damp::array<double, 2> r1{e1.values[0].real(), e1.values[1].real()};
        if (r0[0] > r0[1]) {
            damp::swap(r0[0], r0[1]);
        }
        if (r1[0] > r1[1]) {
            damp::swap(r1[0], r1[1]);
        }
        CHECK(r0[0] == doctest::Approx(r1[0]).epsilon(1e-10));
        CHECK(r0[1] == doctest::Approx(r1[1]).epsilon(1e-10));

        // Transform is invertible: T * T^{-1} ≈ I via reconstruct from factors
        // x = T x_b ⇒ check that C_bal = C T matches bal.balanced.C (already by construction)
        CHECK(max_abs_diff(bal.balanced.D, sys.D) < 1e-15);
    }

    TEST_CASE("hankel_singular_values matches balanced_realization sigma") {
        const auto sys = coupled_2nd_order();
        const auto bal = balanced_realization(sys);
        const auto hsv = hankel_singular_values(sys);
        REQUIRE(bal.success);
        REQUIRE(hsv.has_value());
        CHECK((*hsv)(0) == doctest::Approx(bal.sigma(0)));
        CHECK((*hsv)(1) == doctest::Approx(bal.sigma(1)));
    }

    TEST_CASE("balanced_reduction drops the weakly coupled mode") {
        const auto sys = two_scale_plant();
        const auto bal = balanced_realization(sys);
        REQUIRE(bal.success);

        // Second HSV should be much smaller than the first
        CHECK(bal.sigma(1) < 0.1 * bal.sigma(0));

        const auto red = balanced_reduction<1>(sys);
        REQUIRE(red.success);
        CHECK(red.error_bound == doctest::Approx(2.0 * bal.sigma(1)));
        CHECK(red.sigma_kept(0) == doctest::Approx(bal.sigma(0)));

        // Reduced model is stable (1st-order lag, A < 0)
        CHECK(red.reduced.A(0, 0) < 0.0);

        // DC gain of truncation is close to full-order DC (not exact — that is MatchDC's job)
        const auto dc_full = cont_dc(sys);
        const auto dc_red = cont_dc(red.reduced);
        // Full DC ≈ 1/1 + 0.05*0.05/20 = 1 + 0.000125 = 1.000125
        // Truncation should still be order-1 accurate
        CHECK(std::abs(dc_red(0, 0) - dc_full(0, 0)) < 0.05);
    }

    TEST_CASE("balanced_reduction MatchDC preserves DC gain exactly") {
        const auto sys = two_scale_plant();
        const auto red = balanced_reduction<1>(sys, ModelReductionMethod::MatchDC);
        REQUIRE(red.success);

        const auto dc_full = cont_dc(sys);
        const auto dc_red = cont_dc(red.reduced);
        CHECK(std::abs(dc_red(0, 0) - dc_full(0, 0)) < 1e-10);

        // Reduced pole stays stable
        CHECK(red.reduced.A(0, 0) < 0.0);
    }

    TEST_CASE("balanced_reduction of a 4th-order cascade keeps dominant dynamics") {
        // Four real poles: −1, −2, −50, −80 with decaying input/output coupling
        StateSpace<4, 1, 1> sys{
            .A = Matrix<4, 4>{
                {-1.0, 0.0, 0.0, 0.0},
                {0.0, -2.0, 0.0, 0.0},
                {0.0, 0.0, -50.0, 0.0},
                {0.0, 0.0, 0.0, -80.0}
            },
            .B = Matrix<4, 1>{{1.0}, {0.8}, {0.02}, {0.01}},
            .C = Matrix<1, 4>{{1.0, 0.8, 0.02, 0.01}},
            .D = Matrix<1, 1>{{0.0}},
            .Ts = 0.0
        };

        const auto hsv = hankel_singular_values(sys);
        REQUIRE(hsv.has_value());
        // Descending
        for (size_t i = 1; i < 4; ++i) {
            CHECK((*hsv)(i - 1) >= (*hsv)(i));
        }
        // Fast weakly-coupled modes → small tail
        CHECK((*hsv)(2) + (*hsv)(3) < 0.2 * ((*hsv)(0) + (*hsv)(1)));

        const auto red = balanced_reduction<2>(sys);
        REQUIRE(red.success);
        CHECK(red.error_bound == doctest::Approx(2.0 * ((*hsv)(2) + (*hsv)(3))));

        // Reduced DC should stay near full DC
        const auto dc_full = cont_dc(sys);
        const auto dc_red = cont_dc(red.reduced);
        CHECK(std::abs(dc_red(0, 0) - dc_full(0, 0)) / std::abs(dc_full(0, 0)) < 0.05);

        // MatchDC residualization for order 2 is even tighter on DC
        const auto md = balanced_reduction<2>(sys, ModelReductionMethod::MatchDC);
        REQUIRE(md.success);
        const auto dc_md = cont_dc(md.reduced);
        CHECK(std::abs(dc_md(0, 0) - dc_full(0, 0)) < 1e-9);
    }

    TEST_CASE("discrete balanced_realization succeeds on a Schur-stable plant") {
        // Discrete double-integrator-ish stable plant
        StateSpace<2, 1, 1> sys{
            .A = Matrix<2, 2>{{0.9, 0.1}, {0.0, 0.8}},
            .B = Matrix<2, 1>{{0.05}, {0.1}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .D = Matrix<1, 1>{{0.0}},
            .Ts = 0.01
        };

        const auto bal = balanced_realization(sys);
        REQUIRE(bal.success);
        CHECK(bal.sigma(0) >= bal.sigma(1));
        CHECK(bal.balanced.is_discrete());
        CHECK(bal.balanced.Ts == doctest::Approx(0.01));

        const auto red = balanced_reduction<1>(sys);
        REQUIRE(red.success);
        CHECK(red.reduced.is_discrete());
    }

    TEST_CASE("balanced_realization fails cleanly on an unstable plant") {
        StateSpace<1, 1, 1> sys{
            .A = Matrix<1, 1>{{1.0}}, // unstable continuous pole
            .B = Matrix<1, 1>{{1.0}},
            .C = Matrix<1, 1>{{1.0}},
            .Ts = 0.0
        };
        const auto bal = balanced_realization(sys);
        CHECK_FALSE(bal.success);
        CHECK_FALSE(hankel_singular_values(sys).has_value());
    }


    TEST_CASE("controllability_gramian rejects an unstable A") {
        // A X + X Aᵀ + B Bᵀ = 0 has the unique solution X = −0.5 for A = 1, which
        // is not a Gramian (the integral diverges). It must not be reported.
        const Matrix<1, 1> A{{1.0}};
        const Matrix<1, 1> B{{1.0}};
        CHECK_FALSE(stability::controllability_gramian(A, B).has_value());
        CHECK_FALSE(stability::observability_gramian(A, B).has_value());
        const Matrix<1, 1> Ad{{2.0}};
        CHECK_FALSE(stability::controllability_gramian(Ad, B, true).has_value());
    }
} // TEST_SUITE

TEST_SUITE("eliminate_states") {

    TEST_CASE("residualizing the motor current gives the closed-form speed model") {
        // i = (u − Ke ω)/R  ⇒  J ω̇ = Kt (u − Ke ω)/R − b ω
        constexpr auto red = eliminate_states<0>(dc_motor);
        static_assert(red.has_value());
        CHECK(red->A(0, 0) == doctest::Approx(-10.1).epsilon(1e-12));
        CHECK(red->B(0, 0) == doctest::Approx(100.0).epsilon(1e-12));
        CHECK(red->C(0, 0) == doctest::Approx(-0.1).epsilon(1e-12)); // y_i = (u − Ke ω)/R
        CHECK(red->C(1, 0) == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(red->D(0, 0) == doctest::Approx(1.0).epsilon(1e-12));
        CHECK(red->D(1, 0) == doctest::Approx(0.0));
        CHECK(max_abs_diff(cont_dc(*red), cont_dc(dc_motor)) < 1e-12);
    }

    TEST_CASE("Truncate drops the state with its rows and columns") {
        const auto red = eliminate_states<0>(dc_motor, ModelReductionMethod::Truncate);
        REQUIRE(red.has_value());
        CHECK(red->A(0, 0) == -0.1);
        CHECK(red->B(0, 0) == 0.0);
        CHECK(red->C(0, 0) == 0.0);
        CHECK(red->C(1, 0) == 1.0);
        CHECK(max_abs_diff(red->D, dc_motor.D) == 0.0);
    }

    TEST_CASE("kept states stay in their original order") {
        const auto mid = eliminate_states<1>(three_poles);
        REQUIRE(mid.has_value());
        CHECK(mid->A(0, 0) == -1.0);
        CHECK(mid->A(1, 1) == -2.0);
        CHECK(mid->D(0, 0) == doctest::Approx(0.01));
        CHECK(cont_dc(*mid)(0, 0) == doctest::Approx(cont_dc(three_poles)(0, 0)));

        const auto outer = eliminate_states<2, 0>(three_poles);
        REQUIRE(outer.has_value());
        CHECK(outer->A(0, 0) == -100.0);
        CHECK(outer->D(0, 0) == doctest::Approx(1.5));
    }

    TEST_CASE("discrete residualization holds the eliminated state at steady state") {
        // x₀⁺ = x₀ ⇒ x₀ = (1 − 0.5)⁻¹ u: full DC 1/0.5 + 1/0.1 = 12
        const StateSpace<2, 1, 1> sys{
            .A = Matrix<2, 2>{{0.5, 0.0}, {0.0, 0.9}},
            .B = Matrix<2, 1>{{1.0}, {1.0}},
            .C = Matrix<1, 2>{{1.0, 1.0}},
            .Ts = 0.1,
        };
        const auto red = eliminate_states<0>(sys);
        REQUIRE(red.has_value());
        CHECK(red->Ts == 0.1);
        CHECK(red->A(0, 0) == doctest::Approx(0.9));
        CHECK(red->D(0, 0) == doctest::Approx(2.0));
        CHECK(red->C(0, 0) / (1.0 - red->A(0, 0)) * red->B(0, 0) + red->D(0, 0) == doctest::Approx(12.0));
    }

    TEST_CASE("an integrator among the eliminated states cannot be residualized") {
        const StateSpace<2, 1, 1> sys{
            .A = Matrix<2, 2>{{-1.0, 1.0}, {0.0, 0.0}},
            .B = Matrix<2, 1>{{0.0}, {1.0}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
        };
        CHECK_FALSE(eliminate_states<1>(sys).has_value());
        CHECK(eliminate_states<1>(sys, ModelReductionMethod::Truncate).has_value());
    }
} // TEST_SUITE

TEST_SUITE("balanced reduction of a model with negligible Hankel singular values") {
    TEST_CASE_TEMPLATE("square-root truncation keeps the leading σ; the realization reports the singular tail", T, float, double) {
        // The suite's 4th-order cascade: σ = 0.646, 0.0138, 3.7e-6, 2.5e-8 (scipy). In float the
        // trailing σ sit below NX·64·ε·σ₁, which once refused hankel_singular_values and every
        // balanced_reduction.
        const StateSpace<4, 1, 1, T> sys = StateSpace<4, 1, 1>{
            .A = Matrix<4, 4>{{-1.0, 0.0, 0.0, 0.0}, {0.0, -2.0, 0.0, 0.0}, {0.0, 0.0, -50.0, 0.0}, {0.0, 0.0, 0.0, -80.0}},
            .B = Matrix<4, 1>{{1.0}, {0.8}, {0.02}, {0.01}},
            .C = Matrix<1, 4>{{1.0, 0.8, 0.02, 0.01}},
        }
                                               .template as<T>();
        const auto hsv = hankel_singular_values(sys);
        REQUIRE(hsv.has_value());
        CHECK((*hsv)(0) == doctest::Approx(0.646246).epsilon(1e-4));
        CHECK((*hsv)(1) == doctest::Approx(0.0137552).epsilon(1e-3));

        const auto red = balanced_reduction<2>(sys);
        REQUIRE(red.success);
        CHECK(red.status == DesignStatus::ok);
        const auto dc_full = cont_dc(sys.template as<double>());
        const auto dc_red = cont_dc(red.reduced.template as<double>());
        CHECK(std::abs(dc_red(0, 0) - dc_full(0, 0)) / std::abs(dc_full(0, 0)) < 0.05);

        if constexpr (std::is_same_v<T, float>) {
            const auto bal = balanced_realization(sys);
            CHECK_FALSE(bal.success);
            CHECK(bal.status == DesignStatus::singular);
            CHECK(balanced_reduction<2>(sys, ModelReductionMethod::MatchDC).status == DesignStatus::singular);
        }
    }

    TEST_CASE("a non-finite model is refused with bad_weights") {
        StateSpace<2, 1, 1> sys{
            .A = Matrix<2, 2>{{-1.0, 0.0}, {0.0, -2.0}},
            .B = Matrix<2, 1>{{1.0}, {1.0}},
            .C = Matrix<1, 2>{{1.0, std::numeric_limits<double>::quiet_NaN()}},
        };
        CHECK(balanced_realization(sys).status == DesignStatus::bad_weights);
        CHECK(balanced_reduction<1>(sys).status == DesignStatus::bad_weights);
        CHECK_FALSE(hankel_singular_values(sys).has_value());
    }
} // TEST_SUITE
