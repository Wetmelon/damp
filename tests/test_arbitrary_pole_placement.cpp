// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>

#include "damp/backend.hpp"
#include "damp/design/pole_placement.hpp"
#include "damp/math/complex.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/svd.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using namespace damp::design;
using damp::complex;

namespace {

// Rigorously verify that A − B·K has exactly the requested Jordan structure.
//
// The Jordan structure at an eigenvalue λ is fixed by the rank sequence
// rₖ = rank((Acl − λI)ᵏ): the kernel dimension dₖ = N − rₖ must equal
// Σ min(pᵢ, k) over the block orders pᵢ. Matching this for k = 1…max(pᵢ)
// pins down every block size and the algebraic multiplicity — and unlike
// eigenvalue extraction it stays well-conditioned for defective spectra.
//
// Singular values of the kᵗʰ power are thresholded against σ_max(Acl − λI)ᵏ:
// a relative-to-self rank (MATLAB® semantics) would call the numerically-zero
// high powers "full rank", so the operator scale must come from the base matrix.
template<size_t N, size_t NB>
bool has_jordan_structure(
    const Matrix<N, N, double>&           Acl,
    const damp::array<JordanBlock<>, NB>& blocks,
    double                                rank_tol = 1e-7
) {
    using Cplx = complex<double>;
    const Matrix<N, N, Cplx> Ac = Acl.template as<Cplx>();
    // Operator scale for the rank threshold: tie it to the closed-loop matrix,
    // not to (Acl − λI), which can be the zero matrix for a semisimple λ.
    const double          acl_scale = mat::svd(Ac).singular_values[0];
    damp::array<bool, NB> done{};

    for (size_t bi = 0; bi < NB; ++bi) {
        if (done[bi]) {
            continue;
        }
        const Cplx lam = blocks[bi].eigenvalue;

        // Gather all mini-block orders for this eigenvalue.
        damp::array<size_t, NB> orders{};
        size_t                  g = 0;
        size_t                  max_p = 0;
        for (size_t bj = 0; bj < NB; ++bj) {
            if (!done[bj] && damp::abs(blocks[bj].eigenvalue - lam) < 1e-7) {
                orders[g] = blocks[bj].size;
                max_p = std::max(max_p, blocks[bj].size);
                ++g;
                done[bj] = true;
            }
        }

        // Shifted pencil (Acl − λI); its powers scale at most as (‖Acl‖+|λ|)ᵏ.
        Matrix<N, N, Cplx> shifted = Ac;
        for (size_t i = 0; i < N; ++i) {
            shifted(i, i) -= lam;
        }
        const double base = acl_scale + damp::abs(lam);

        // Kernel dimension of each power, thresholded at rank_tol·baseᵏ.
        Matrix<N, N, Cplx> power = Matrix<N, N, Cplx>::identity();
        for (size_t k = 1; k <= max_p; ++k) {
            power = power * shifted;
            const auto   s = mat::svd(power);
            const double thresh = rank_tol * std::pow(base, static_cast<double>(k));
            size_t       rank = 0;
            for (size_t i = 0; i < N; ++i) {
                if (s.singular_values[i] > thresh) {
                    ++rank;
                }
            }
            size_t expected = 0;
            for (size_t t = 0; t < g; ++t) {
                expected += std::min(orders[t], k);
            }
            if ((N - rank) != expected) {
                return false;
            }
        }
    }
    return true;
}

// A continuous plant with the given A and B (outputs = states).
template<size_t NX, size_t NU>
constexpr StateSpace<NX, NU, NX> plant(const Matrix<NX, NX, double>& A, const Matrix<NX, NU, double>& B) {
    return StateSpace<NX, NU, NX>{.A = A, .B = B, .C = Matrix<NX, NX, double>::identity()};
}

template<size_t R, size_t C>
double max_abs(const Matrix<R, C, double>& Mx) {
    double m = 0.0;
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            m = std::max(m, std::abs(Mx(i, j)));
        }
    }
    return m;
}

} // namespace

TEST_SUITE("place_jordan") {
    TEST_CASE("SISO distinct poles match the unique Ackermann gain") {
        // SISO pole placement has a unique solution, so place_jordan must agree
        // with the independent Ackermann formula to machine precision.
        const Matrix<3, 3, double>          A = {{0, 1, 0}, {0, 0, 1}, {-6, -11, -6}};
        const Matrix<3, 1, double>          B = {{0}, {0}, {1}};
        const damp::array<JordanBlock<>, 3> blocks = {{{{-1, 0}, 1}, {{-2, 0}, 1}, {{-3, 0}, 1}}};

        const auto K = place_jordan(plant(A, B), blocks);
        REQUIRE(K.success);

        // Ackermann by hand on the companion form: s³ + 6s² + 11s + 6 is already the open loop.
        CHECK(max_abs(K.K) < 1e-9);
        const auto Kknv = place(plant(A, B), {-1.0, -2.0, -3.0});
        REQUIRE(Kknv.success);
        CHECK(max_abs(K.K - Kknv.K) < 1e-9);
        CHECK(has_jordan_structure(A - (B * K.K), blocks));
    }

    TEST_CASE("MIMO distinct poles") {
        const Matrix<4, 4, double>          A = {{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {1, 2, 3, 4}};
        const Matrix<4, 2, double>          B = {{0, 0}, {1, 0}, {0, 0}, {0, 1}};
        const damp::array<JordanBlock<>, 4> blocks = {{{{-1, 0}, 1}, {{-2, 0}, 1}, {{-3, 0}, 1}, {{-4, 0}, 1}}};

        const auto K = place_jordan(plant(A, B), blocks);
        REQUIRE(K.success);
        CHECK(has_jordan_structure(A - (B * K.K), blocks));
    }

    TEST_CASE("Repeated semisimple eigenvalue (two 1×1 blocks)") {
        // -2 with two independent eigenvectors: A − B·K must be exactly −2·I.
        const Matrix<2, 2, double>          A = {{0, 1}, {0, 0}};
        const Matrix<2, 2, double>          B = Matrix<2, 2, double>::identity();
        const damp::array<JordanBlock<>, 2> blocks = {{{{-2, 0}, 1}, {{-2, 0}, 1}}};

        const auto K = place_jordan(plant(A, B), blocks);
        REQUIRE(K.success);
        const Matrix<2, 2, double> Acl = A - (B * K.K);
        CHECK(has_jordan_structure(Acl, blocks));
        // Semisimple ⇒ minimal polynomial degree 1 ⇒ (Acl + 2I) = 0.
        CHECK(max_abs(Acl + (2.0 * Matrix<2, 2, double>::identity())) < 1e-10);
    }

    TEST_CASE("Defective real block: deadbeat with a single 3×3 Jordan block") {
        // All poles at 0 in one chain ⇒ nilpotent of index exactly 3.
        const Matrix<3, 3, double>          A = {{0, 1, 0}, {0, 0, 1}, {-6, -11, -6}};
        const Matrix<3, 1, double>          B = {{0}, {0}, {1}};
        const damp::array<JordanBlock<>, 1> blocks = {{{{0, 0}, 3}}};

        const auto K = place_jordan(plant(A, B), blocks);
        REQUIRE(K.success);
        const Matrix<3, 3, double> Acl = A - (B * K.K);
        CHECK(has_jordan_structure(Acl, blocks));
        // Independent confirmation: nilpotency index exactly 3.
        const Matrix<3, 3, double> A2 = Acl * Acl;
        CHECK(max_abs(A2) > 0.1);        // index > 2
        CHECK(max_abs(A2 * Acl) < 1e-9); // (A−BK)³ = 0
    }

    TEST_CASE("single-input place() assigns a repeated pole as the place_jordan chain") {
        // With one input a double pole has one eigenvector: place() assigns it as one
        // 2×2 Jordan block, the same gain place_jordan gives for that structure.
        const Matrix<2, 2, double>          A = {{0, 1}, {0, 0}};
        const Matrix<2, 1, double>          B = {{0}, {1}};
        const damp::array<JordanBlock<>, 1> blocks = {{{{-3, 0}, 2}}};

        const auto Kj = place_jordan(plant(A, B), blocks);
        REQUIRE(Kj.success);
        const Matrix<2, 2, double> Acl = A - (B * Kj.K);
        CHECK(has_jordan_structure(Acl, blocks)); // one 2×2 block, not two 1×1

        const auto Kp = place(plant(A, B), {-3.0, -3.0});
        REQUIRE(Kp.success);
        CHECK(max_abs(Kp.K - Kj.K) < 1e-9);
    }

    TEST_CASE("Complex conjugate pair yields a real gain") {
        const Matrix<2, 2, double>          A = {{0, 1}, {0, 0}};
        const Matrix<2, 1, double>          B = {{0}, {1}};
        const damp::array<JordanBlock<>, 2> blocks = {{{{-1, 2}, 1}, {{-1, -2}, 1}}};

        const auto K = place_jordan(plant(A, B), blocks);
        REQUIRE(K.success);
        CHECK(std::isfinite(K.K(0, 0)));
        CHECK(has_jordan_structure(A - (B * K.K), blocks));
    }

    TEST_CASE("Defective complex pair: two 2×2 Jordan blocks at -1±1j") {
        // Four-integrator chain, controllability indices {2,2}.
        const Matrix<4, 4, double>          A = {{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {0, 0, 0, 0}};
        const Matrix<4, 2, double>          B = {{0, 0}, {1, 0}, {0, 0}, {0, 1}};
        const damp::array<JordanBlock<>, 2> blocks = {{{{-1, 1}, 2}, {{-1, -1}, 2}}};

        const auto K = place_jordan(plant(A, B), blocks);
        REQUIRE(K.success);
        const Matrix<4, 4, double> Acl = A - (B * K.K);
        CHECK(has_jordan_structure(Acl, blocks));

        // Independent confirmation via the real quadratic factor
        // Q = Acl² − 2σ·Acl + |λ|²·I, which is nilpotent of index 2: scale-
        // relative ‖Q²‖ ≪ ‖Q‖² while Q itself is not nilpotent of index 1.
        const auto                 I4 = Matrix<4, 4, double>::identity();
        const Matrix<4, 4, double> Q = (Acl * Acl) - (-2.0 * Acl) + (2.0 * I4);
        const double               nq = Q.norm();
        CHECK((Q * Q).norm() < 1e-10 * nq * nq); // (defective) index 2
        CHECK(nq > 0.1);                         // not already zero
    }

    TEST_CASE("Mixed structure: a defective real block, a simple pole, a complex pair") {
        // n = 5: a 2×2 Jordan block at -2, a 1×1 at -5, and a simple pair ±j.
        const Matrix<5, 5, double> A = {
            {0, 1, 0, 0, 0},
            {0, 0, 1, 0, 0},
            {0, 0, 0, 1, 0},
            {0, 0, 0, 0, 1},
            {0, 0, 0, 0, 0},
        };
        const Matrix<5, 2, double>          B = {{0, 0}, {1, 0}, {0, 0}, {0, 0}, {0, 1}};
        const damp::array<JordanBlock<>, 4> blocks = {{{{-2, 0}, 2}, {{-5, 0}, 1}, {{0, 1}, 1}, {{0, -1}, 1}}};

        const auto K = place_jordan(plant(A, B), blocks);
        REQUIRE(K.success);
        CHECK(has_jordan_structure(A - (B * K.K), blocks));
    }

    TEST_CASE("Inadmissible requests are refused") {
        const Matrix<3, 3, double> A = {{0, 1, 0}, {0, 0, 1}, {-6, -11, -6}};
        const Matrix<3, 1, double> B = {{0}, {0}, {1}};

        // Block sizes do not sum to NX (2 ≠ 3).
        const damp::array<JordanBlock<>, 1> too_few = {{{{-1, 0}, 2}}};
        CHECK_FALSE(place_jordan(plant(A, B), too_few).success);

        // Single-input system: an eigenvalue cannot have two independent
        // eigenvectors (g = 2 > NU = 1), so a repeated semisimple pole fails.
        const damp::array<JordanBlock<>, 3> two_chains = {{{{-1, 0}, 1}, {{-1, 0}, 1}, {{-2, 0}, 1}}};
        CHECK_FALSE(place_jordan(plant(A, B), two_chains).success);
    }

    TEST_CASE("Conjugate blocks with different orders are refused (a real gain mirrors λ onto λ̄)") {
        const Matrix<4, 4, double>          A = {{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {1, 2, 3, 4}};
        const Matrix<4, 2, double>          B = {{0, 0}, {1, 0}, {0, 0}, {0, 1}};
        const damp::array<JordanBlock<>, 3> mismatched = {{{{-1, 1}, 2}, {{-1, -1}, 1}, {{-1, -1}, 1}}};
        const auto                          r = place_jordan(plant(A, B), mismatched);
        CHECK_FALSE(r.success);
        CHECK(r.status == DesignStatus::not_placeable);
        const damp::array<JordanBlock<>, 2> matched = {{{{-1, 1}, 2}, {{-1, -1}, 2}}};
        CHECK(place_jordan(plant(A, B), matched).success);
    }

    TEST_CASE("Slow distinct eigenvalues are not merged, and a non-finite one is refused") {
        const Matrix<2, 2, double>          A = {{0, 1}, {0, 0}};
        const Matrix<2, 1, double>          B = {{0}, {1}};
        const damp::array<JordanBlock<>, 2> slow = {{{{-1e-8, 0}, 1}, {{-2e-8, 0}, 1}}};
        CHECK(place_jordan(plant(A, B), slow).success);
        const damp::array<JordanBlock<>, 2> bad = {{{{std::nan(""), 0}, 1}, {{-2, 0}, 1}}};
        CHECK(place_jordan(plant(A, B), bad).status == DesignStatus::invalid_parameter);
    }
}

TEST_SUITE("place_jordan optimization") {
    // A 4-state, 2-input system where the canonical parameter is far from optimal
    // (its gain is ~49), leaving plenty of room for the optimizer to improve.
    const Matrix<4, 4, double>          A = {{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {1, 2, 3, 4}};
    const Matrix<4, 2, double>          B = {{0, 0}, {1, 0}, {0, 0}, {0, 1}};
    const damp::array<JordanBlock<>, 4> poles = {{{{-1, 0}, 1}, {{-2, 0}, 1}, {{-3, 0}, 1}, {{-4, 0}, 1}}};

    TEST_CASE("Robust (Method 1) preserves the spectrum and improves conditioning/gain") {
        const auto canonical = place_jordan(plant(A, B), poles);
        REQUIRE(canonical.success);
        const auto r = place_jordan(plant(A, B), poles, {.objective = JordanObjective::ConditionNumber, .alpha = 1.0});
        REQUIRE(r.success);

        // The optimizer only moves within exact-placement solutions.
        CHECK(has_jordan_structure(A - (B * r.K), poles));
        CHECK(r.cond_fro < 30.0);               // well-conditioned eigenvectors
        CHECK(r.gain_fro < canonical.K.norm()); // and less gain than canonical
        CHECK(r.converged);
    }

    TEST_CASE("Minimum gain (alpha=0) drives the feedback norm below the robust solution") {
        const auto robust = place_jordan(plant(A, B), poles, {.objective = JordanObjective::ConditionNumber, .alpha = 1.0});
        const auto mingain = place_jordan(plant(A, B), poles, {.objective = JordanObjective::ConditionNumber, .alpha = 0.0});
        REQUIRE(robust.success);
        REQUIRE(mingain.success);

        CHECK(has_jordan_structure(A - (B * mingain.K), poles));
        CHECK(mingain.gain_fro <= robust.gain_fro + 1e-6);
    }

    TEST_CASE("Method 2 minimizes departure from normality") {
        const auto m1 = place_jordan(plant(A, B), poles, {.objective = JordanObjective::ConditionNumber, .alpha = 1.0});
        const auto m2 = place_jordan(plant(A, B), poles, {.objective = JordanObjective::DepartureFromNormality, .alpha = 1.0});
        REQUIRE(m1.success);
        REQUIRE(m2.success);

        CHECK(has_jordan_structure(A - (B * m2.K), poles));
        CHECK(m2.departure_fro >= 0.0);
        CHECK(m2.departure_fro <= m1.departure_fro + 1e-6); // M2 targets this measure
    }

    TEST_CASE("Optimization preserves a defective Jordan structure") {
        // Two 2×2 Jordan blocks at -1±1j on the integrator chain.
        const Matrix<4, 4, double>          Achain = {{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {0, 0, 0, 0}};
        const Matrix<4, 2, double>          Bchain = {{0, 0}, {1, 0}, {0, 0}, {0, 1}};
        const damp::array<JordanBlock<>, 2> defective = {{{{-1, 1}, 2}, {{-1, -1}, 2}}};

        const auto r = place_jordan(plant(Achain, Bchain), defective, {.objective = JordanObjective::ConditionNumber, .alpha = 0.5});
        REQUIRE(r.success);
        CHECK(has_jordan_structure(Achain - (Bchain * r.K), defective));
    }

    TEST_CASE("Inadmissible request is refused") {
        const damp::array<JordanBlock<>, 1> too_few = {{{{-1, 0}, 2}}}; // 2 ≠ 4
        CHECK_FALSE(place_jordan(plant(A, B), too_few, {.objective = JordanObjective::ConditionNumber}).success);
    }

    TEST_CASE("Runs at compile time (constexpr)") {
        // The analytic-gradient BFGS fits the constexpr op budget for small
        // systems, so a robust gain can be baked in at compile time.
        constexpr Matrix<2, 2, double>          Ac = {{0, 1}, {0, 0}};
        constexpr Matrix<2, 1, double>          Bc = {{0}, {1}};
        constexpr damp::array<JordanBlock<>, 2> pc = {{{{-1, 0}, 1}, {{-2, 0}, 1}}};
        constexpr auto                          r = place_jordan(plant(Ac, Bc), pc, {.objective = JordanObjective::ConditionNumber});
        static_assert(r.success);
        static_assert(r.converged);
        static_assert(r.cond_fro < 10.0);
    }

    TEST_CASE("JordanPlaceResult.as<float>() preserves metrics") {
        const auto r = place_jordan(plant(A, B), poles, {.objective = JordanObjective::ConditionNumber, .alpha = 1.0});
        REQUIRE(r.success);
        const auto rf = r.as<float>();
        CHECK(rf.converged == r.converged);
        CHECK(rf.iterations == r.iterations);
        CHECK(static_cast<double>(rf.cond_fro) == doctest::Approx(r.cond_fro).epsilon(1e-6));
        CHECK(static_cast<double>(rf.gain_fro) == doctest::Approx(r.gain_fro).epsilon(1e-6));
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                CHECK(static_cast<double>(rf.K(i, j)) == doctest::Approx(r.K(i, j)).epsilon(1e-6));
            }
        }
    }
}
