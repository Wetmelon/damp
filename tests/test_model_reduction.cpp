// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>

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

    TEST_CASE("balreal: Hankel SV descending and gramians balanced") {
        const auto sys = coupled_2nd_order();
        const auto bal = balreal(sys);
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

    TEST_CASE("balreal: preserves transfer function (DC + poles)") {
        const auto sys = coupled_2nd_order();
        const auto bal = balreal(sys);
        REQUIRE(bal.success);

        // DC gain invariant under similarity
        const auto dc0 = cont_dc(sys);
        const auto dc1 = cont_dc(bal.balanced);
        CHECK(std::abs(dc0(0, 0) - dc1(0, 0)) < 1e-10);

        // Eigenvalues of A preserved (similarity)
        const auto e0 = mat::compute_eigenvalues(sys.A);
        const auto e1 = mat::compute_eigenvalues(bal.balanced.A);
        REQUIRE(e0.converged);
        REQUIRE(e1.converged);

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

    TEST_CASE("hankelsv matches balreal sigma") {
        const auto sys = coupled_2nd_order();
        const auto bal = balreal(sys);
        const auto hsv = hankelsv(sys);
        REQUIRE(bal.success);
        REQUIRE(hsv.has_value());
        CHECK((*hsv)(0) == doctest::Approx(bal.sigma(0)));
        CHECK((*hsv)(1) == doctest::Approx(bal.sigma(1)));
    }

    TEST_CASE("balred drops the weakly coupled mode") {
        const auto sys = two_scale_plant();
        const auto bal = balreal(sys);
        REQUIRE(bal.success);

        // Second HSV should be much smaller than the first
        CHECK(bal.sigma(1) < 0.1 * bal.sigma(0));

        const auto red = balred<1>(sys);
        REQUIRE(red.success);
        CHECK(red.error_bound == doctest::Approx(2.0 * bal.sigma(1)));
        CHECK(red.sigma_kept(0) == doctest::Approx(bal.sigma(0)));

        // Reduced model is stable (1st-order lag, A < 0)
        CHECK(red.reduced.A(0, 0) < 0.0);

        // DC gain of truncation is close to full-order DC (not exact — that is modred's job)
        const auto dc_full = cont_dc(sys);
        const auto dc_red = cont_dc(red.reduced);
        // Full DC ≈ 1/1 + 0.05*0.05/20 = 1 + 0.000125 = 1.000125
        // Truncation should still be order-1 accurate
        CHECK(std::abs(dc_red(0, 0) - dc_full(0, 0)) < 0.05);
    }

    TEST_CASE("modred MatchDC preserves DC gain exactly") {
        const auto sys = two_scale_plant();
        const auto red = modred<1>(sys, ModelReductionMethod::MatchDC);
        REQUIRE(red.success);

        const auto dc_full = cont_dc(sys);
        const auto dc_red = cont_dc(red.reduced);
        CHECK(std::abs(dc_red(0, 0) - dc_full(0, 0)) < 1e-10);

        // Reduced pole stays stable
        CHECK(red.reduced.A(0, 0) < 0.0);
    }

    TEST_CASE("modred Truncate matches balred") {
        const auto sys = two_scale_plant();
        const auto a = balred<1>(sys);
        const auto b = modred<1>(sys, ModelReductionMethod::Truncate);
        REQUIRE(a.success);
        REQUIRE(b.success);
        CHECK(max_abs_diff(a.reduced.A, b.reduced.A) < 1e-12);
        CHECK(max_abs_diff(a.reduced.B, b.reduced.B) < 1e-12);
        CHECK(max_abs_diff(a.reduced.C, b.reduced.C) < 1e-12);
        CHECK(a.error_bound == doctest::Approx(b.error_bound));
    }

    TEST_CASE("balred of a 4th-order cascade keeps dominant dynamics") {
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

        const auto hsv = hankelsv(sys);
        REQUIRE(hsv.has_value());
        // Descending
        for (size_t i = 1; i < 4; ++i) {
            CHECK((*hsv)(i - 1) >= (*hsv)(i));
        }
        // Fast weakly-coupled modes → small tail
        CHECK((*hsv)(2) + (*hsv)(3) < 0.2 * ((*hsv)(0) + (*hsv)(1)));

        const auto red = balred<2>(sys);
        REQUIRE(red.success);
        CHECK(red.error_bound == doctest::Approx(2.0 * ((*hsv)(2) + (*hsv)(3))));

        // Reduced DC should stay near full DC
        const auto dc_full = cont_dc(sys);
        const auto dc_red = cont_dc(red.reduced);
        CHECK(std::abs(dc_red(0, 0) - dc_full(0, 0)) / std::abs(dc_full(0, 0)) < 0.05);

        // MatchDC residualization for order 2 is even tighter on DC
        const auto md = modred<2>(sys, ModelReductionMethod::MatchDC);
        REQUIRE(md.success);
        const auto dc_md = cont_dc(md.reduced);
        CHECK(std::abs(dc_md(0, 0) - dc_full(0, 0)) < 1e-9);
    }

    TEST_CASE("discrete balreal succeeds on a Schur-stable plant") {
        // Discrete double-integrator-ish stable plant
        StateSpace<2, 1, 1> sys{
            .A = Matrix<2, 2>{{0.9, 0.1}, {0.0, 0.8}},
            .B = Matrix<2, 1>{{0.05}, {0.1}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .D = Matrix<1, 1>{{0.0}},
            .Ts = 0.01
        };

        const auto bal = balreal(sys);
        REQUIRE(bal.success);
        CHECK(bal.sigma(0) >= bal.sigma(1));
        CHECK(bal.balanced.is_discrete());
        CHECK(bal.balanced.Ts == doctest::Approx(0.01));

        const auto red = balred<1>(sys);
        REQUIRE(red.success);
        CHECK(red.reduced.is_discrete());
    }

    TEST_CASE("balreal fails cleanly on an unstable plant") {
        StateSpace<1, 1, 1> sys{
            .A = Matrix<1, 1>{{1.0}}, // unstable continuous pole
            .B = Matrix<1, 1>{{1.0}},
            .C = Matrix<1, 1>{{1.0}},
            .Ts = 0.0
        };
        const auto bal = balreal(sys);
        CHECK_FALSE(bal.success);
        CHECK_FALSE(hankelsv(sys).has_value());
    }

} // TEST_SUITE
