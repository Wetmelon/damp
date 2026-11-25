// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

#include "damp/backend.hpp"
#include "damp/design/qp.hpp"
#include "damp/math/math.hpp"
#include "damp/matlab.hpp"
#include "damp/matrix/matrix.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

// MATLAB® quadprog documentation example:
//   H = [1 -1; -1 2], f = [-2; -6], A = [1 1; -1 2; 2 1], b = [2; 2; 3]
//   -> x = [2/3; 4/3], fval = -74/9, constraints 1 and 2 active.
struct QuadprogExample {
    Matrix<2, 2> H{{1.0, -1.0}, {-1.0, 2.0}};
    ColVec<2>    f{-2.0, -6.0};
    Matrix<3, 2> A{{1.0, 1.0}, {-1.0, 2.0}, {2.0, 1.0}};
    ColVec<3>    b{2.0, 2.0, 3.0};
};

/// The minimizer of a QP of 4 variables and 8 rows by enumeration: the oracle of the active set.
struct EnumeratedQP {
    bool        feasible{false};
    long double objective{0.0L};
};

/// Gaussian elimination with partial pivoting on the leading n × n of M; false on a pivot at roundoff (dependent rows).
bool solve_in_place(long double M[8][8], long double rhs[8], size_t n) {
    for (size_t c = 0; c < n; ++c) {
        size_t      pivot = c;
        long double column = 0.0L;
        for (size_t r = 0; r < n; ++r) {
            column = std::max(column, std::fabs(M[r][c]));
            if (r > c && std::fabs(M[r][c]) > std::fabs(M[pivot][c])) {
                pivot = r;
            }
        }
        if (!(std::fabs(M[pivot][c]) > 1e-13L * column)) {
            return false;
        }
        for (size_t j = 0; j < n; ++j) {
            std::swap(M[pivot][j], M[c][j]);
        }
        std::swap(rhs[pivot], rhs[c]);
        for (size_t r = c + 1; r < n; ++r) {
            const long double m = M[r][c] / M[c][c];
            for (size_t j = c; j < n; ++j) {
                M[r][j] -= m * M[c][j];
            }
            rhs[r] -= m * rhs[c];
        }
    }
    for (size_t i = n; i-- > 0;) {
        for (size_t j = i + 1; j < n; ++j) {
            rhs[i] -= M[i][j] * rhs[j];
        }
        rhs[i] /= M[i][i];
    }
    return true;
}

/// H's 1-norm condition number, ‖H‖₁·‖H⁻¹‖₁, in long double (the largest finite value for a singular H).
long double condition_1(const Matrix<4, 4, float>& H) {
    long double h_norm = 0.0L;
    long double inv_norm = 0.0L;
    for (size_t c = 0; c < 4; ++c) {
        long double M[8][8]{};
        long double e[8]{};
        long double column = 0.0L;
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                M[i][j] = static_cast<long double>(H(i, j));
            }
            column += std::fabs(static_cast<long double>(H(i, c)));
        }
        e[c] = 1.0L;
        if (!solve_in_place(M, e, 4)) {
            return std::numeric_limits<long double>::max();
        }
        long double inv_column = 0.0L;
        for (size_t i = 0; i < 4; ++i) {
            inv_column += std::fabs(e[i]);
        }
        h_norm = std::max(h_norm, column);
        inv_norm = std::max(inv_norm, inv_column);
    }
    return h_norm * inv_norm;
}

/// Every set of at most 4 active rows: its KKT system solved in long double, kept when the point
/// meets all the rows and its multipliers are nonnegative. A strictly convex QP has one such point.
EnumeratedQP enumerate_qp(const Matrix<4, 4, float>& H, const ColVec<4, float>& f, const Matrix<8, 4, float>& A, const ColVec<8, float>& b) {
    EnumeratedQP best{};
    for (unsigned mask = 0; mask < 256U; ++mask) {
        size_t rows[4] = {};
        size_t count = 0;
        for (size_t i = 0; i < 8 && count <= 4; ++i) {
            if ((mask & (1U << i)) != 0U) {
                if (count < 4) {
                    rows[count] = i;
                }
                ++count;
            }
        }
        if (count > 4) {
            continue;
        }
        long double M[8][8] = {};
        long double rhs[8] = {};
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                M[i][j] = static_cast<long double>(H(i, j));
            }
            rhs[i] = -static_cast<long double>(f(i));
        }
        for (size_t a = 0; a < count; ++a) {
            for (size_t j = 0; j < 4; ++j) {
                M[4 + a][j] = static_cast<long double>(A(rows[a], j));
                M[j][4 + a] = static_cast<long double>(A(rows[a], j));
            }
            rhs[4 + a] = static_cast<long double>(b(rows[a]));
        }
        if (!solve_in_place(M, rhs, 4 + count)) {
            continue;
        }
        bool kkt = true;
        for (size_t a = 0; a < count; ++a) {
            kkt = kkt && rhs[4 + a] >= -1e-12L * (1.0L + std::fabs(rhs[4 + a]));
        }
        long double objective = 0.0L;
        for (size_t i = 0; i < 8 && kkt; ++i) {
            long double ax = -static_cast<long double>(b(i));
            long double terms = std::fabs(static_cast<long double>(b(i)));
            for (size_t j = 0; j < 4; ++j) {
                ax += static_cast<long double>(A(i, j)) * rhs[j];
                terms += std::fabs(static_cast<long double>(A(i, j)) * rhs[j]);
            }
            kkt = kkt && ax <= 1e-12L * terms;
        }
        if (!kkt) {
            continue;
        }
        for (size_t i = 0; i < 4; ++i) {
            long double hx = 0.0L;
            for (size_t j = 0; j < 4; ++j) {
                hx += static_cast<long double>(H(i, j)) * rhs[j];
            }
            objective += rhs[i] * ((0.5L * hx) + static_cast<long double>(f(i)));
        }
        if (!best.feasible || objective < best.objective) {
            best = EnumeratedQP{true, objective};
        }
    }
    return best;
}

} // namespace

TEST_CASE("solve_qp reproduces the MATLAB® quadprog example") {
    const QuadprogExample qp;

    const auto res = design::solve_qp(qp.H, qp.f, qp.A, qp.b);
    REQUIRE(res.success);
    CHECK(res.status == design::QPStatus::Success);

    CHECK(res.x(0) == doctest::Approx(2.0 / 3.0).epsilon(1e-9));
    CHECK(res.x(1) == doctest::Approx(4.0 / 3.0).epsilon(1e-9));
    CHECK(res.objective == doctest::Approx(-74.0 / 9.0).epsilon(1e-9));

    // Multipliers (hand KKT solve): lambda = [28/9, 4/9, 0].
    CHECK(res.lambda(0) == doctest::Approx(28.0 / 9.0).epsilon(1e-9));
    CHECK(res.lambda(1) == doctest::Approx(4.0 / 9.0).epsilon(1e-9));
    CHECK(res.lambda(2) == doctest::Approx(0.0));
}

TEST_CASE("solve_qp satisfies the KKT conditions") {
    const QuadprogExample qp;
    const auto            res = design::solve_qp(qp.H, qp.f, qp.A, qp.b);
    REQUIRE(res.success);

    // Primal feasibility + dual feasibility + complementary slackness.
    for (size_t i = 0; i < 3; ++i) {
        double slack = -qp.b(i);
        for (size_t c = 0; c < 2; ++c) {
            slack += qp.A(i, c) * res.x(c);
        }
        CHECK(slack <= 1e-9);          // A x <= b
        CHECK(res.lambda(i) >= -1e-9); // lambda >= 0
        CHECK(res.lambda(i) * slack == doctest::Approx(0.0).epsilon(1e-9));
    }

    // Stationarity: H x + f + A' lambda = 0.
    for (size_t r = 0; r < 2; ++r) {
        double g = qp.f(r);
        for (size_t c = 0; c < 2; ++c) {
            g += qp.H(r, c) * res.x(c);
        }
        for (size_t i = 0; i < 3; ++i) {
            g += qp.A(i, r) * res.lambda(i);
        }
        CHECK(g == doctest::Approx(0.0).epsilon(1e-9));
    }
}

TEST_CASE("solve_qp is constexpr") {
    constexpr QuadprogExample qp;
    constexpr auto            res = design::solve_qp(qp.H, qp.f, qp.A, qp.b);
    static_assert(res.success);
    static_assert(res.x(0) - (2.0 / 3.0) < 1e-9 && res.x(0) - (2.0 / 3.0) > -1e-9);
    static_assert(res.x(1) - (4.0 / 3.0) < 1e-9 && res.x(1) - (4.0 / 3.0) > -1e-9);
}

TEST_CASE("solve_qp returns the unconstrained minimum when nothing binds") {
    const QuadprogExample qp;

    // All bounds at the unbounded sentinel: rows are skipped entirely.
    ColVec<3> b_open{};
    for (size_t i = 0; i < 3; ++i) {
        b_open(i) = design::unbounded_bound<double>();
    }
    const auto res = design::solve_qp(qp.H, qp.f, qp.A, b_open);
    REQUIRE(res.success);
    CHECK(res.iterations == 0);

    // x = -H^{-1} f = [10; 8]
    CHECK(res.x(0) == doctest::Approx(10.0));
    CHECK(res.x(1) == doctest::Approx(8.0));
    CHECK(res.lambda(0) == doctest::Approx(0.0));
}

TEST_CASE("solve_qp clamps against simple box constraints") {
    // min 1/2||x - (2,2)||^2 s.t. x <= (1,1): solution x = (1,1), lambda = (1,1).
    const auto   H = Matrix<2, 2>::identity();
    const ColVec f{-2.0, -2.0};
    const Matrix A{{1.0, 0.0}, {0.0, 1.0}};
    const ColVec b{1.0, 1.0};

    const auto res = design::solve_qp(H, f, A, b);
    REQUIRE(res.success);
    CHECK(res.x(0) == doctest::Approx(1.0));
    CHECK(res.x(1) == doctest::Approx(1.0));
    CHECK(res.lambda(0) == doctest::Approx(1.0));
    CHECK(res.lambda(1) == doctest::Approx(1.0));
    CHECK(res.objective == doctest::Approx(-3.0));
}

TEST_CASE("solve_qp detects infeasibility") {
    // x <= -1 and -x <= -1 (i.e. x >= 1) cannot both hold.
    const Matrix H{{1.0}};
    const ColVec f{0.0};
    const Matrix A{{1.0}, {-1.0}};
    const ColVec b{-1.0, -1.0};

    const auto res = design::solve_qp(H, f, A, b);
    CHECK(!res.success);
    CHECK(res.status == design::QPStatus::Infeasible);
}

TEST_CASE("solve_qp rejects an indefinite Hessian") {
    const Matrix H{{-1.0}};
    const ColVec f{1.0};
    const Matrix A{{1.0}};
    const ColVec b{1.0};

    const auto res = design::solve_qp(H, f, A, b);
    CHECK(!res.success);
    CHECK(res.status == design::QPStatus::NotPositiveDefinite);
}

TEST_CASE("solve_qp works in float") {
    const QuadprogExample qp;
    const auto            res = design::solve_qp(
        qp.H.as<float>(), ColVec<2, float>{qp.f.as<float>()},
        qp.A.as<float>(), ColVec<3, float>{qp.b.as<float>()}
    );
    REQUIRE(res.success);
    CHECK(res.x(0) == doctest::Approx(2.0F / 3.0F).epsilon(1e-5));
    CHECK(res.x(1) == doctest::Approx(4.0F / 3.0F).epsilon(1e-5));
}

TEST_CASE("float solve_qp holds a row on a small variable beside a large one") {
    // min ½‖x‖² − [1000, 0.02]·x s.t. x₂ ≤ 0.01. A slack tolerance of 64·(NV+1)·ε·‖x‖∞ (0.02 in
    // float) accepted the unconstrained x₂ = 0.02 as feasible and returned Success at twice the limit.
    const auto res = design::solve_qp(
        Matrix<2, 2, float>::identity(), ColVec{-1000.0F, -0.02F},
        Matrix{{0.0F, 1.0F}}, ColVec{0.01F}
    );
    REQUIRE(res.success);
    CHECK(res.x(0) == doctest::Approx(1000.0F).epsilon(1e-6));
    CHECK(res.x(1) <= 0.01F * (1.0F + 1e-5F));
}

TEST_CASE_TEMPLATE("solve_qp meets a row through a heavily weighted slack", T, float, double) {
    //! min ½x² − 2x + ρε² s.t. x ≤ ½, x + ε ≥ 1, ε ≥ 0, with a soft-constrained MPC's ρ = 10⁵.
    //! Once x ≤ ½ is active, x + ε ≥ 1 moves ε alone: in the metric of H⁻¹ its normal meets the
    //! active one at a sine of 1/√(2ρ) = 2.2·10⁻³. The Gram factor of the active normals tested
    //! the sine's square, 5·10⁻⁶, against a float floor of 4.6·10⁻⁵ and returned Degenerate.
    constexpr T           rho = static_cast<T>(1.0e5);
    const Matrix<2, 2, T> H{{T{1}, T{0}}, {T{0}, T{2} * rho}};
    const ColVec<2, T>    f{T{-2}, T{0}};
    const Matrix<3, 2, T> A{{T{1}, T{0}}, {T{-1}, T{-1}}, {T{0}, T{-1}}};
    const ColVec<3, T>    b{T{0.5}, T{-1}, T{0}};

    const auto res = design::solve_qp(H, f, A, b);
    REQUIRE(res.status == design::QPStatus::Success);
    const double row_tol = 64.0 * static_cast<double>(eps<T>()); // each active row to the roundoff of its terms
    CHECK(static_cast<double>(res.x(0)) == doctest::Approx(0.5).epsilon(row_tol));
    CHECK(static_cast<double>(res.x(1)) == doctest::Approx(0.5).epsilon(row_tol));
    // Stationarity gives λ = [ρ + 3/2, ρ, 0], known to ε over the sine.
    const double lambda_tol = row_tol * damp::sqrt(2.0e5);
    CHECK(static_cast<double>(res.lambda(0)) == doctest::Approx(1.0e5 + 1.5).epsilon(lambda_tol));
    CHECK(static_cast<double>(res.lambda(1)) == doctest::Approx(1.0e5).epsilon(lambda_tol));
    CHECK(res.lambda(2) == T{0});
}

TEST_CASE_TEMPLATE("the factors of the active normals hold them through adds and drops", T, float, double) {
    //! Q·[R; 0] is the set of active normals and Q is orthogonal after every change: the solver
    //! reads the dependence test, the minimizer on the active face and its multipliers from them.
    constexpr size_t   NV = 5;
    unsigned long long s = 99;
    const auto         next = [&s] { // uniform in [−1, 1)
        s = (s * 6364136223846793005ULL) + 1442695040888963407ULL;
        return static_cast<T>((static_cast<double>((s >> 11U) & ((1ULL << 53U) - 1)) / static_cast<double>(1ULL << 52U)) - 1.0);
    };
    const auto random_vector = [&next] {
        ColVec<NV, T> v{};
        for (size_t i = 0; i < NV; ++i) {
            v(i) = next();
        }
        return v;
    };
    const T tol = T{16} * static_cast<T>(NV) * eps<T>(); // of vectors of length about 1

    design::detail::ActiveNormals<NV, NV, T> normals{};
    damp::array<ColVec<NV, T>, NV>           held{}; // the active normals, in the factor's order
    const auto                               check = [&] {
        for (size_t i = 0; i < NV; ++i) {
            for (size_t j = 0; j < NV; ++j) {
                T dot{0};
                for (size_t r = 0; r < NV; ++r) {
                    dot += normals.Q(r, i) * normals.Q(r, j);
                }
                CHECK(damp::abs(dot - ((i == j) ? T{1} : T{0})) <= tol);
            }
        }
        for (size_t a = 0; a < normals.k; ++a) {
            for (size_t r = 0; r < NV; ++r) {
                T entry{0};
                for (size_t j = 0; j <= a; ++j) {
                    entry += normals.Q(r, j) * normals.R(j, a);
                }
                CHECK(damp::abs(entry - held[a](r)) <= tol * T{4});
            }
        }
    };
    const auto add = [&](const ColVec<NV, T>& v) {
        held[normals.k] = v;
        normals.append(normals.rotate(v));
    };
    const auto drop = [&](size_t a, ColVec<NV, T>& candidate) {
        for (size_t i = a; i + 1 < normals.k; ++i) {
            held[i] = held[i + 1];
        }
        normals.remove(a, candidate);
    };

    for (size_t i = 0; i < 4; ++i) {
        add(random_vector());
    }
    check();
    const ColVec<NV, T> v = random_vector();
    ColVec<NV, T>       c = normals.rotate(v); // a candidate turns with Q through a drop
    drop(1, c);
    check();
    const ColVec<NV, T> c_fresh = normals.rotate(v);
    for (size_t i = 0; i < NV; ++i) {
        CHECK(damp::abs(c(i) - c_fresh(i)) <= tol * T{4});
    }
    add(v);
    drop(0, c);
    drop(normals.k - 1, c);
    add(random_vector());
    check();
    REQUIRE(normals.k == 3);

    // The point of the face ñₐᵀy = βₐ nearest y₀ meets every active row, and is y₀ in the free
    // directions: y − y₀ has no part outside the span of the active normals.
    const ColVec<NV, T> y0 = random_vector();
    ColVec<NV, T>       beta{};
    for (size_t a = 0; a < normals.k; ++a) {
        beta(a) = next();
    }
    const ColVec<NV, T> y = normals.nearest(y0, beta);
    for (size_t a = 0; a < normals.k; ++a) {
        T row{0};
        for (size_t r = 0; r < NV; ++r) {
            row += held[a](r) * y(r);
        }
        CHECK(damp::abs(row - beta(a)) <= tol * T{8});
    }
    for (size_t j = normals.k; j < NV; ++j) {
        T free{0};
        for (size_t r = 0; r < NV; ++r) {
            free += normals.Q(r, j) * (y(r) - y0(r));
        }
        CHECK(damp::abs(free) <= tol * T{8});
    }
}

TEST_CASE("float solve_qp agrees with double on random QPs") {
    //! 4000 QPs of 4 variables and 8 rows, the float solve against the double solve of the same
    //! float data. On the Gram factor 8 ended Degenerate, 2 reported Infeasible on a feasible
    //! problem, and one Success was 0.27 off in x with an objective 82 % above the optimum.
    unsigned long long s = 777;
    const auto         next = [&s] { // uniform in [−1, 1)
        s = (s * 6364136223846793005ULL) + 1442695040888963407ULL;
        return (static_cast<double>((s >> 11U) & ((1ULL << 53U) - 1)) / static_cast<double>(1ULL << 52U)) - 1.0;
    };
    size_t degenerate = 0;
    size_t disagree = 0;
    size_t solved = 0;
    double worst_x = 0.0;
    double worst_objective = 0.0;
    for (size_t t = 0; t < 4000; ++t) {
        Matrix<4, 4, float> L{};
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j <= i; ++j) {
                L(i, j) = static_cast<float>(next());
            }
            L(i, i) = 0.2F + damp::abs(L(i, i));
        }
        const Matrix<4, 4, float> H = L * L.transpose();
        ColVec<4, float>          f{};
        Matrix<8, 4, float>       A{};
        ColVec<8, float>          b{};
        for (size_t i = 0; i < 4; ++i) {
            f(i) = static_cast<float>(3.0 * next());
        }
        for (size_t i = 0; i < 8; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                A(i, j) = static_cast<float>(next());
            }
            b(i) = static_cast<float>(next() + 0.6); // feasible on 97 % of them
        }
        const auto single = design::solve_qp(H, f, A, b);
        const auto exact = design::solve_qp(H.as<double>(), ColVec<4>{f.as<double>()}, A.as<double>(), ColVec<8>{b.as<double>()});
        degenerate += (single.status == design::QPStatus::Degenerate) ? 1U : 0U;
        disagree += (single.status != exact.status) ? 1U : 0U;
        if (!single.success || !exact.success) {
            continue;
        }
        ++solved;
        double dx = 0.0;
        double scale = 1.0;
        double objective = 0.0; // of the float point, on the same data in double
        for (size_t i = 0; i < 4; ++i) {
            dx = damp::max(dx, damp::abs(static_cast<double>(single.x(i)) - exact.x(i)));
            scale = damp::max(scale, damp::abs(exact.x(i)));
            double hx = 0.0;
            for (size_t j = 0; j < 4; ++j) {
                hx += static_cast<double>(H(i, j)) * static_cast<double>(single.x(j));
            }
            objective += static_cast<double>(single.x(i)) * ((0.5 * hx) + static_cast<double>(f(i)));
        }
        worst_x = damp::max(worst_x, dx / scale);
        worst_objective = damp::max(worst_objective, damp::abs(objective - exact.objective) / (1.0 + damp::abs(exact.objective)));
    }
    CHECK(degenerate == 0);
    CHECK(disagree <= 4); // a row within the float tolerance of infeasible may read either way
    CHECK(solved > 3800);
    CHECK(worst_x < 1e-2);         // measured 1.1·10⁻³: a vertex is known to ε over the sine its rows meet at
    CHECK(worst_objective < 1e-3); // measured 3·10⁻⁵
}

TEST_CASE_TEMPLATE("solve_qp with exactly dependent rows agrees with an enumeration of the active sets", T, float, double) {
    //! Rows 4 … 7 are sums or differences of rows 0 … 3 (dyadic entries: exact), their bounds
    //! consistent, tighter or contradictory, so every problem holds dependent normals and the
    //! infeasible ones need a Farkas certificate. The factor of H is graded by 1 and by 10: the
    //! roundoff of a transformed normal grows with it, past the dependence floor in double.
    //! On the Gram factor, of 3000 at each grading: 18 and 42 feasible problems Infeasible in
    //! float and 13 and 105 in double, 29 and 752 Degenerate in float, and 4 and 259 float
    //! Success points off the optimum.
    unsigned long long s = 4242;
    const auto         next = [&s] { // uniform in [−1, 1)
        s = (s * 6364136223846793005ULL) + 1442695040888963407ULL;
        return (static_cast<double>((s >> 11U) & ((1ULL << 53U) - 1)) / static_cast<double>(1ULL << 52U)) - 1.0;
    };
    const auto dyadic = [](double v) { return static_cast<float>(std::round(v * 64.0) / 64.0); };
    for (const double grade : {1.0, 0.1}) {
        CAPTURE(grade);
        size_t off_optimum = 0;
        size_t false_infeasible = 0;
        size_t missed_infeasible = 0;
        size_t unsolved = 0;
        size_t infeasible = 0;
        for (size_t t = 0; t < 3000; ++t) {
            Matrix<4, 4, float> L{};
            for (size_t i = 0; i < 4; ++i) {
                for (size_t j = 0; j < i; ++j) {
                    L(i, j) = static_cast<float>(next());
                }
                L(i, i) = static_cast<float>((0.2 + std::abs(next())) * std::pow(grade, static_cast<double>(i) / 3.0));
            }
            const Matrix<4, 4, float> H = L * L.transpose();
            ColVec<4, float>          f{};
            Matrix<8, 4, float>       A{};
            ColVec<8, float>          b{};
            for (size_t i = 0; i < 4; ++i) {
                f(i) = static_cast<float>(3.0 * next());
            }
            for (size_t i = 0; i < 4; ++i) {
                for (size_t j = 0; j < 4; ++j) {
                    A(i, j) = dyadic(next());
                }
                b(i) = dyadic(next() + 0.6);
            }
            for (size_t i = 4; i < 8; ++i) {
                const size_t first = i - 4;
                const size_t second = (i - 3) % 4;
                const float  sign = (t % 2 == 0) ? 1.0F : -1.0F;
                for (size_t j = 0; j < 4; ++j) {
                    A(i, j) = A(first, j) + (sign * A(second, j));
                }
                b(i) = b(first) + (sign * b(second)) + ((t % 3 == 0) ? 0.0F : dyadic(0.5 * next()));
            }
            const EnumeratedQP exact = enumerate_qp(H, f, A, b);
            const auto         res = design::solve_qp(H.template as<T>(), ColVec<4, T>{f.template as<T>()}, A.template as<T>(), ColVec<8, T>{b.template as<T>()});
            infeasible += exact.feasible ? 0U : 1U;
            if (res.status == design::QPStatus::Infeasible) {
                false_infeasible += exact.feasible ? 1U : 0U;
            } else if (res.status != design::QPStatus::Success) {
                ++unsolved;
            } else if (!exact.feasible) {
                ++missed_infeasible;
            } else {
                long double objective = 0.0L; // of the returned point, on the data as given
                for (size_t i = 0; i < 4; ++i) {
                    long double hx = 0.0L;
                    for (size_t j = 0; j < 4; ++j) {
                        hx += static_cast<long double>(H(i, j)) * static_cast<long double>(res.x(j));
                    }
                    objective += static_cast<long double>(res.x(i)) * ((0.5L * hx) + static_cast<long double>(f(i)));
                }
                // Graded by 10, an H can be past what float solves: cond₁(H)·ε near 1 leaves its factor
                // no digit in the weak direction, and the objective off by whatever the rounding gives
                // (one draw with cond₁ 5·10⁶ lands 1.8·10⁵ ε off on clang 23, under 10⁴ ε on g++). Such
                // a draw, 282 of the 3000 in float at grade 0.1 and none in double, is not graded on its
                // objective; the others reach 1.4·10³ ε in float and 5.5·10³ ε in double.
                const bool        gradable = condition_1(H) * static_cast<long double>(eps<T>()) < 0.1L;
                const long double gap = std::fabs(objective - exact.objective) / (1.0L + std::fabs(exact.objective));
                off_optimum += (gradable && gap > 1.0e5L * static_cast<long double>(eps<T>())) ? 1U : 0U;
            }
        }
        CHECK(infeasible > 100);
        CHECK(off_optimum == 0);
        CHECK(false_infeasible == 0);
        CHECK(missed_infeasible == 0);
        CHECK(unsolved == 0);
    }
}

TEST_CASE("equality-constrained solve_qp finds the minimum-norm allocation") {
    // min 1/2||x||^2 s.t. x1 + x2 + x3 = 1: x = (1/3, 1/3, 1/3), equality
    // multiplier -1/3 (exercises the negative-multiplier lower twin).
    const auto      H = Matrix<3, 3>::identity();
    const ColVec<3> f{};
    const Matrix    A{{1.0, 0.0, 0.0}};
    const ColVec<1> b{design::unbounded_bound<double>()}; // no inequality binds
    const Matrix    Aeq{{1.0, 1.0, 1.0}};
    const ColVec    beq{1.0};

    const auto res = design::solve_qp(H, f, A, b, Aeq, beq);
    REQUIRE(res.success);
    for (size_t v = 0; v < 3; ++v) {
        CHECK(res.x(v) == doctest::Approx(1.0 / 3.0).epsilon(1e-9));
    }
    // Equality multiplier = upper twin - lower twin = -1/3.
    const double lambda_eq = res.lambda(1) - res.lambda(2);
    CHECK(lambda_eq == doctest::Approx(-1.0 / 3.0).epsilon(1e-9));
}

TEST_CASE("equality constraint composes with active inequalities") {
    // The quadprog example restricted to the line x1 = x2: along the line the
    // cost decreases until inequality rows 1 and 3 clamp it at x = (1, 1).
    const QuadprogExample qp;
    const Matrix          Aeq{{1.0, -1.0}};
    const ColVec          beq{0.0};

    const auto res = design::solve_qp(qp.H, qp.f, qp.A, qp.b, Aeq, beq);
    REQUIRE(res.success);
    CHECK(res.x(0) == doctest::Approx(1.0).epsilon(1e-9));
    CHECK(res.x(1) == doctest::Approx(1.0).epsilon(1e-9));
}

TEST_CASE("inconsistent equalities are reported infeasible") {
    const Matrix    H{{1.0, 0.0}, {0.0, 1.0}};
    const ColVec<2> f{};
    const Matrix    A{{1.0, 0.0}};
    const ColVec<1> b{design::unbounded_bound<double>()};
    const Matrix    Aeq{{1.0, 1.0}, {1.0, 1.0}};
    const ColVec    beq{1.0, 2.0}; // x1+x2 = 1 and = 2

    const auto res = design::solve_qp(H, f, A, b, Aeq, beq);
    CHECK(!res.success);
    CHECK(res.status == design::QPStatus::Infeasible);
}

TEST_CASE("equality-constrained solve_qp is constexpr") {
    constexpr auto res = [] {
        constexpr auto      H = Matrix<2, 2>::identity();
        constexpr ColVec    f{-1.0, 0.0};
        constexpr Matrix    A{{1.0, 0.0}};
        constexpr ColVec<1> b{design::unbounded_bound<double>()};
        constexpr Matrix    Aeq{{1.0, 1.0}};
        constexpr ColVec    beq{1.0};
        return design::solve_qp(H, f, A, b, Aeq, beq);
    }();
    static_assert(res.success);
    // min 1/2||x||^2 - x1 s.t. x1+x2=1 -> x = (1, 0).
    static_assert(res.x(0) - 1.0 < 1e-9 && res.x(0) - 1.0 > -1e-9);
    static_assert(res.x(1) < 1e-9 && res.x(1) > -1e-9);
}

TEST_CASE("matlab::quadprog equality overload returns the same minimizer") {
    const auto      H = Matrix<3, 3>::identity();
    const ColVec<3> f{};
    const Matrix    A{{1.0, 0.0, 0.0}};
    const ColVec<1> b{design::unbounded_bound<double>()};
    const Matrix    Aeq{{1.0, 1.0, 1.0}};
    const ColVec    beq{1.0};

    const auto x = *matlab::quadprog(H, f, A, b, Aeq, beq);
    for (size_t v = 0; v < 3; ++v) {
        CHECK(x(v) == doctest::Approx(1.0 / 3.0).epsilon(1e-9));
    }
}

TEST_CASE("matlab::quadprog alias returns the same minimizer") {
    const QuadprogExample qp;
    const auto            x = *matlab::quadprog(qp.H, qp.f, qp.A, qp.b);
    const auto            ref = design::solve_qp(qp.H, qp.f, qp.A, qp.b).x;
    CHECK(x(0) == doctest::Approx(ref(0)));
    CHECK(x(1) == doctest::Approx(ref(1)));
}

TEST_CASE("solve_qp is invariant to scaling the cost") {
    // MATLAB quadprog example; scaling H and f by s leaves the minimizer.
    for (const double s : {1e-14, 1.0, 1e14}) {
        const Matrix<2, 2> H{{s, -s}, {-s, 2.0 * s}};
        const ColVec<2>    f{-2.0 * s, -6.0 * s};
        const Matrix       A{{1.0, 1.0}, {-1.0, 2.0}, {2.0, 1.0}};
        const ColVec       b{2.0, 2.0, 3.0};
        const auto         res = design::solve_qp(H, f, A, b);
        CAPTURE(s);
        REQUIRE(res.success);
        CHECK(res.x(0) == doctest::Approx(2.0 / 3.0).epsilon(1e-9));
        CHECK(res.x(1) == doctest::Approx(4.0 / 3.0).epsilon(1e-9));
    }
}

TEST_CASE("solve_qp is invariant to scaling the constraint rows") {
    const Matrix H{{1.0, -1.0}, {-1.0, 2.0}};
    const ColVec f{-2.0, -6.0};
    for (const double s : {1e-9, 1.0, 1e9}) {
        const Matrix<3, 2> A{{s, s}, {-s, 2.0 * s}, {2.0 * s, s}};
        const ColVec<3>    b{2.0 * s, 2.0 * s, 3.0 * s};
        const auto         res = design::solve_qp(H, f, A, b);
        CAPTURE(s);
        REQUIRE(res.success);
        CHECK(res.x(0) == doctest::Approx(2.0 / 3.0).epsilon(1e-9));
        CHECK(res.x(1) == doctest::Approx(4.0 / 3.0).epsilon(1e-9));
    }
}

TEST_CASE("solve_qp proves infeasibility once NV constraints are active, and reports no false Degenerate") {
    //! With NV normals active a new one lies in their span, but roundoff left κ above the
    //! dependence threshold: the full step then hit "cannot activate more than NV" and returned
    //! Degenerate (118 of these 4000 random QPs), which costs solve_miqp its infeasibility proofs.
    unsigned long long s = 777;
    const auto         next = [&s] { // uniform in [−1, 1)
        s = (s * 6364136223846793005ULL) + 1442695040888963407ULL;
        return (static_cast<double>((s >> 11U) & ((1ULL << 53U) - 1)) / static_cast<double>(1ULL << 52U)) - 1.0;
    };
    size_t degenerate = 0;
    size_t infeasible = 0;
    size_t disagree = 0;
    for (size_t t = 0; t < 4000; ++t) {
        Matrix<4, 4> L{};
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j <= i; ++j) {
                L(i, j) = next();
            }
            L(i, i) = 0.2 + damp::abs(L(i, i));
        }
        const Matrix<4, 4> H = L * L.transpose();
        ColVec<4>          f{};
        Matrix<8, 4>       A{};
        ColVec<8>          b{};
        for (size_t i = 0; i < 4; ++i) {
            f(i) = 3.0 * next();
        }
        for (size_t i = 0; i < 8; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                A(i, j) = next();
            }
            b(i) = next() - 0.3;
        }
        const auto res = design::solve_qp(H, f, A, b);
        degenerate += (res.status == design::QPStatus::Degenerate) ? 1U : 0U;
        if (res.status == design::QPStatus::Infeasible) {
            ++infeasible;
            disagree += design::solve_qp_interior_point(H, f, A, b).success ? 1U : 0U;
        }
    }
    CHECK(degenerate == 0);
    CHECK(infeasible > 1000);
    CHECK(disagree == 0);
}

TEST_CASE("solve_qp and the interior point refuse non-finite data instead of returning success") {
    //! A NaN f gave x = NaN with Success, and a NaN bᵢ dropped its row silently (x = [3, 3]).
    const double       nan = std::numeric_limits<double>::quiet_NaN();
    const Matrix<2, 2> H = Matrix<2, 2>::identity();
    const ColVec       f{-3.0, -3.0};
    const Matrix       A{{1.0, 1.0}};
    CHECK(design::solve_qp(H, ColVec<2>{nan, -3.0}, A, ColVec{2.0}).status == design::QPStatus::NonFinite);
    CHECK(design::solve_qp(H, f, A, ColVec<1>{nan}).status == design::QPStatus::NonFinite);
    CHECK(design::solve_qp(H, f, Matrix<1, 2>{{nan, 1.0}}, ColVec{2.0}).status == design::QPStatus::NonFinite);
    CHECK(design::solve_qp_interior_point(H, f, A, ColVec<1>{nan}).status == design::QPStatus::NonFinite);
    CHECK(design::solve_qp_interior_point(H, ColVec<2>{nan, -3.0}, A, ColVec{2.0}).status == design::QPStatus::NonFinite);
    const auto ok = design::solve_qp(H, f, A, ColVec{2.0});
    REQUIRE(ok.success);
    CHECK(ok.x(0) == doctest::Approx(1.0));
}
