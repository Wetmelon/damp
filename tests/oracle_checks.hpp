// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file oracle_checks.hpp
 * @brief The oracle conformance suite's case loader, tolerances and per-family checks;
 *        the test_oracle_*.cpp files run them.
 */

#include <array>
#include <cstddef>
#include <cstdint>
#include <tuple>
#include <type_traits>
#include <utility>

#include "damp/backend.hpp"
#include "damp/design/lyapunov.hpp"
#include "damp/design/pole_placement.hpp"
#include "damp/design/riccati.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/decomposition.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/functions.hpp"
#include "damp/matrix/generalized_eigen.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/solve.hpp"
#include "damp/matrix/svd.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/systems/zeros.hpp"
#include "damp/systems/zpk.hpp"
#include "oracle_cases.hpp"
#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

// Oracle conformance: expm, care and dare, the factorizations (LU, QR, Cholesky,
// LDLᵀ, SVD), eigenvalues and the Schur basis, logm, sqrtm, lyap / dlyap, ZOH and
// Tustin discretization, the SISO conversions (ss2tf, ss2zpk, tf2zpk, zpk2tf,
// tf2ss) and pole placement, against scipy cases (or exact rational arithmetic on
// the same doubles) that tools/gen_oracle.py writes into oracle_cases.hpp (run the
// script by hand; make never calls Python). Fixed-seed property tests close the
// file: QR, SVD, LU / Cholesky solves, Lyapunov residuals and ZOH resampling
// through the logarithm (the d2c step) on random matrices from a SplitMix64 stream.
//
// Passing does not mean matching scipy's digits (G2.4). The primary checks are
// properties a correct answer has whatever the reference says:
//   expm       relative Frobenius error ≤ C_EXPM · n · ε · max(1, κ), where κ =
//              scipy.linalg.expm_cond is the relative condition number of exp at
//              A. A method with backward error of order n·ε (the design target of
//              Al-Mohy and Higham's scaling and squaring) has forward error of
//              order κ·n·ε. expm has no residual, so this is its primary check.
//   care/dare  relative Riccati residual ‖res‖ / Σ‖terms‖ ≤ C_RES · n · ε ·
//              max(1, res_floor), X symmetric and PSD to C_RES · n · ε, and the
//              closed loop A − BK strictly stable (Re λ < 0 / |λ| < 1) by Damp's own
//              eigenvalues. res_floor is the residual of the scipy X in exact
//              arithmetic, in units of n·ε: rounding X moves the residual through the
//              closed-loop Lyapunov operator, ε‖A_cl‖‖X‖, which can dwarf Σ‖terms‖
//              (care_cheap_control_3: fast poles near −1/√R, res_floor 87). No answer
//              in double does better than the correctly rounded one; elsewhere it is
//              below 1 and the bound is C_RES · n · ε.
// The scipy X is a secondary regression check: relative error ≤ C_GOLD · n · ε ·
// max(1, κ), κ a first-order condition estimate of X; skipped when that bound
// reaches 1 (no digit of X is determined). The scipy closed-loop margin (LAPACK
// eigenvalues) is the oracle for Damp's eigenvalue routine, to the square root
// of that bound, since a nearly defective pole moves as √ of the perturbation.
//
// The constants leave about 3× headroom over the worst case observed with GCC
// 14 and -ffast-math: expm 14 (a 60 rad rotation, through ~7 squarings),
// residual 10 (care_random_5x2 in float), gold 0.2.
//
// The later families compare only what is unique: reconstructions and residuals
// (backward errors, in the units of the bound the algorithm meets), orthogonality,
// |det|, singular values, eigenvalues matched to the reference, the triangular
// factors fixed by a sign convention (R with diag ≥ 0, Cholesky L, unpivoted LDLᵀ),
// polynomial coefficients, Markov parameters for tf2ss (any realization) and
// closed-loop eigenvalues for place (K is not unique for several inputs). A root of
// multiplicity m (a Jordan block of size m) moves as the m-th root of the
// perturbation, so roots are checked as |Δz|ᵐ against the stored sensitivity.
// Their constants leave about 3× headroom too: C_FACTOR over 3.1 (orthogonality of
// the SVD's V in the property runs), C_FWD over 4.5 (sqrtm_jordan_3), C_ROOT over
// 1.3 (a zero of ss2tf_biproper_3 in double). ZOH is an expm and shares C_EXPM.
//
// Everything is measured in the case's natural units x₀ = D⁻¹x (D = 1 unless the
// case is graded), so a case that only mixes units must be solved as well as
// the benign problem underneath. Metrics are evaluated in double for the double
// and the float design alike (float: on the float-rounded data); ε is the design
// type's epsilon. Every case runs in double. Float runs unless float_ok = false,
// which the generator sets when κ·ε_float ≥ 1, or when float data sits within 4× of
// a refusal floor the library documents (cholesky's pivot floor, qr_decompose's
// rank-loss floor, lyap / dlyap's singularity band): there a float run tests the
// refusal, not the answer.
//
// @see "Functions of Matrices" (Higham, 2008), §3.1 (condition) and §10.3 (expm)
// @see "Numerical condition of the algebraic Riccati equation" (Byers, 1985)
// @see "Accuracy and Stability of Numerical Algorithms" (Higham, 2002), §9, §10, §11, §16, §19

namespace oracle_checks {

inline constexpr double C_EXPM = 50.0;  ///< expm forward-error constant
inline constexpr double C_RES = 32.0;   ///< Riccati residual / symmetry / definiteness constant
inline constexpr double C_GOLD = 100.0; ///< secondary scipy-X constant (on top of the κ estimate)

template<typename T>
constexpr double eps_of() {
    return static_cast<double>(eps<T>());
}

template<typename T>
constexpr const char* type_name() {
    return std::is_same_v<T, float> ? "float" : "double";
}

/// Row-major case data → Matrix<R, C, T> (rounded to T).
template<size_t R, size_t C, typename T>
constexpr Matrix<R, C, T> load(const std::array<double, R * C>& a) {
    Matrix<R, C, T> M;
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            M(i, j) = static_cast<T>(a[(i * C) + j]);
        }
    }
    return M;
}

/// diag(l) · M · diag(r)
template<size_t R, size_t C>
constexpr Matrix<R, C> scale(const std::array<double, R>& l, Matrix<R, C> M, const std::array<double, C>& r) {
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            M(i, j) *= l[i] * r[j];
        }
    }
    return M;
}

template<size_t N>
constexpr std::array<double, N> recip(const std::array<double, N>& d) {
    std::array<double, N> out{};
    for (size_t i = 0; i < N; ++i) {
        out[i] = 1.0 / d[i];
    }
    return out;
}

template<size_t N>
constexpr std::array<double, N> ones() {
    std::array<double, N> out{};
    out.fill(1.0);
    return out;
}

template<size_t R, size_t C>
constexpr double rel_diff(const Matrix<R, C>& X, const Matrix<R, C>& Ref) {
    const double err = Matrix<R, C>(X - Ref).norm();
    const double den = Ref.norm();
    return den > 0.0 ? err / den : err;
}

/// m scaled by a volatile 1: an input the optimizer cannot constant-fold.
template<typename M>
M opaque(const M& m) {
    volatile int one = 1;
    return m * static_cast<typename M::value_type>(one);
}

template<typename Tuple, typename F>
void for_each_case(F&& f) {
    [&]<size_t... I>(std::index_sequence<I...>) {
        (f.template operator()<std::tuple_element_t<I, Tuple>>(), ...);
    }(std::make_index_sequence<std::tuple_size_v<Tuple>>{});
}

// Cases this suite caught the library failing. Library numerics are out of scope
// for the suite's patch, so they run under doctest::may_fail() below instead of
// in the gating test cases. Specialize known_bug<C> = true for a case found
// failing, with the symptom, and remove it once its bug is fixed. None is open:
// the Newton polish and residual certificate (CARE and DARE), the PBH test on
// the CARE certificate and the bit-pattern finiteness check on the DARE iterate
// closed every case this list held.
template<typename C>
inline constexpr bool known_bug = false;

// ---- expm ----

template<typename C>
constexpr double expm_bound(double e) {
    return C_EXPM * static_cast<double>(C::n) * e * damp::max(1.0, C::kappa);
}

/// ‖E₀ − E₀_ref‖_F / ‖E₀_ref‖_F in natural units E₀ = D⁻¹ E D.
template<typename C, typename T>
constexpr double expm_error(const Matrix<C::n, C::n, T>& E) {
    const auto di = recip(C::units);
    const auto E0 = scale(di, Matrix<C::n, C::n>(E), C::units);
    const auto R0 = scale(di, load<C::n, C::n, double>(C::expA), C::units);
    return rel_diff(E0, R0);
}

template<typename C, typename T>
void check_expm() {
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const double err = expm_error<C, T>(mat::expm(load<C::n, C::n, T>(C::A)).value());
    INFO("relative error ", err, ", κ ", C::kappa);
    CHECK(err <= expm_bound<C>(eps_of<T>()));
}

/// C1.2: the same expm at compile time (abstract IEEE) and at run time (-ffast-math).
template<typename C, typename T>
void check_expm_constexpr() {
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    constexpr auto E_ct = mat::expm(load<C::n, C::n, T>(C::A)).value();
    static_assert(expm_error<C, T>(E_ct) <= expm_bound<C>(eps_of<T>()));
    const auto E_rt = mat::expm(opaque(load<C::n, C::n, T>(C::A))).value();
    CHECK(rel_diff(Matrix<C::n, C::n>(E_rt), Matrix<C::n, C::n>(E_ct)) <= expm_bound<C>(eps_of<T>()));
}

// ---- Riccati ----

/// Primary and secondary metrics of a Riccati design, in natural units, in double.
struct Report {
    double residual{1.0};  ///< ‖res‖_F / Σ‖terms‖_F
    double asymmetry{1.0}; ///< ‖X − Xᵀ‖_F / ‖X‖_F
    double min_eig{-1.0};  ///< λ_min(X) / ‖X‖_F
    double margin{2.0};    ///< max Re λ (care) or max |λ| (dare) of A − BK
    bool   eig_converged{false};
    double gold{1.0}; ///< ‖X − X_scipy‖_F / ‖X_scipy‖_F
};

template<typename C>
constexpr double residual_bound(double e) {
    return C_RES * static_cast<double>(C::nx) * e;
}

/// residual_bound, raised to the residual the correctly rounded scipy X attains (res_floor).
template<typename C>
constexpr double residual_floor_bound(double e) {
    return residual_bound<C>(e) * damp::max(1.0, C::res_floor);
}

template<typename C>
constexpr double gold_bound(double e) {
    return C_GOLD * static_cast<double>(C::nx) * e * damp::max(1.0, C::kappa);
}

template<size_t NX>
constexpr double max_margin(const ColVec<NX, damp::complex<double>>& ev, bool discrete) {
    double m = discrete ? ev[0].abs() : ev[0].real();
    for (size_t i = 1; i < NX; ++i) {
        m = damp::max(m, discrete ? ev[i].abs() : ev[i].real());
    }
    return m;
}

template<bool Discrete, typename C, typename T>
constexpr Report riccati_report(const Matrix<C::nx, C::nx, T>& Xt) {
    constexpr size_t NX = C::nx;
    constexpr size_t NU = C::nu;
    using MX = Matrix<NX, NX>;
    const auto d = C::units;
    const auto di = recip(d);
    // Natural units: A₀ = D⁻¹AD, B₀ = D⁻¹B, Q₀ = DQD, N₀ = DN, X₀ = DXD.
    const MX             A = scale(di, MX(load<NX, NX, T>(C::A)), d);
    const Matrix<NX, NU> B = scale(di, Matrix<NX, NU>(load<NX, NU, T>(C::B)), ones<NU>());
    const MX             Q = scale(d, MX(load<NX, NX, T>(C::Q)), d);
    const Matrix<NU, NU> R(load<NU, NU, T>(C::R));
    const Matrix<NX, NU> N = scale(d, Matrix<NX, NU>(load<NX, NU, T>(C::N)), ones<NU>());
    const MX             X = scale(d, MX(Xt), d);

    // K = (R + BᵀXB)⁻¹(BᵀXA + Nᵀ) or R⁻¹(BᵀX + Nᵀ); the terms of the residual.
    const auto K = Discrete ? mat::solve(Matrix<NU, NU>(R + B.t() * X * B), Matrix<NU, NX>(B.t() * X * A + N.t()))
                            : mat::solve(R, Matrix<NU, NX>(B.t() * X + N.t()));
    if (!K) {
        return Report{};
    }
    const MX           Acl = A - B * K.value();
    damp::array<MX, 4> terms{};
    if constexpr (Discrete) {
        terms = {MX(A.t() * X * A), X * -1.0, MX((A.t() * X * B + N) * K.value()) * -1.0, Q};
    } else {
        terms = {MX(A.t() * X), MX(X * A), MX((X * B + N) * K.value()) * -1.0, Q};
    }
    MX     res = MX::zeros();
    double den = 0.0;
    for (const MX& t : terms) {
        res = res + t;
        den += t.norm();
    }

    Report     rep;
    const auto ex = mat::compute_eigenvalues(MX((X + X.t()) * 0.5));
    const auto ec = mat::compute_eigenvalues(Acl);
    rep.residual = res.norm() / den;
    rep.asymmetry = MX(X - X.t()).norm() / X.norm();
    rep.eig_converged = ex.success && ec.success;
    double lmin = ex.values[0].real();
    for (size_t i = 1; i < NX; ++i) {
        lmin = damp::min(lmin, ex.values[i].real());
    }
    rep.min_eig = lmin / X.norm();
    rep.margin = max_margin<NX>(ec.values, Discrete);
    if constexpr (requires { C::X; }) {
        rep.gold = rel_diff(X, scale(d, load<NX, NX, double>(C::X), d));
    }
    return rep;
}

template<bool Discrete, typename C, typename T>
constexpr auto run_riccati(const Matrix<C::nx, C::nx, T>& A) {
    const auto B = load<C::nx, C::nu, T>(C::B);
    const auto Q = load<C::nx, C::nx, T>(C::Q);
    const auto R = load<C::nu, C::nu, T>(C::R);
    const auto N = load<C::nx, C::nu, T>(C::N);
    if constexpr (Discrete) {
        return dare(A, B, Q, R, N);
    } else {
        return care(A, B, Q, R, N);
    }
}

template<bool Discrete, typename C, typename T>
void check_riccati() {
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const auto X = run_riccati<Discrete, C, T>(load<C::nx, C::nx, T>(C::A));
    CHECK(X.has_value());
    if (!X) {
        return;
    }
    const Report r = riccati_report<Discrete, C, T>(X.value());
    const double e = eps_of<T>();
    INFO("residual ", r.residual, ", asymmetry ", r.asymmetry, ", λmin(X)/‖X‖ ", r.min_eig);
    INFO("margin ", r.margin, " (scipy ", C::margin, "), gold error ", r.gold, ", κ ", C::kappa);
    CHECK(r.eig_converged);
    CHECK(r.residual <= residual_floor_bound<C>(e));
    CHECK(r.asymmetry <= residual_bound<C>(e));
    CHECK(r.min_eig >= -residual_bound<C>(e));
    CHECK(r.margin < (Discrete ? 1.0 : 0.0));
    const double gb = gold_bound<C>(e);
    if (gb < 1.0) {
        CHECK(r.gold <= gb);
        CHECK(damp::abs(r.margin - C::margin) <= damp::sqrt(gb) * damp::max(1.0, damp::abs(C::margin)));
    }
}

template<bool Discrete, typename C, typename T>
void check_riccati_fails() {
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    CHECK_FALSE(run_riccati<Discrete, C, T>(load<C::nx, C::nx, T>(C::A)).has_value());
}

/// C1.2: the design at compile time certifies, and the run-time design agrees with it.
template<bool Discrete, typename C, typename T>
void check_riccati_constexpr() {
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    constexpr auto X_ct = run_riccati<Discrete, C, T>(load<C::nx, C::nx, T>(C::A));
    static_assert(X_ct.has_value());
    constexpr Report r = riccati_report<Discrete, C, T>(X_ct.value());
    static_assert(r.residual <= residual_floor_bound<C>(eps_of<T>()));
    static_assert(r.margin < (Discrete ? 1.0 : 0.0));
    const auto X_rt = run_riccati<Discrete, C, T>(opaque(load<C::nx, C::nx, T>(C::A)));
    CHECK(X_rt.has_value());
    if (X_rt) {
        using MX = Matrix<C::nx, C::nx>;
        CHECK(rel_diff(MX(X_rt.value()), MX(X_ct.value())) <= gold_bound<C>(eps_of<T>()));
    }
}

template<bool Discrete, typename Tuple, bool Bugs>
void check_riccati_cases() {
    for_each_case<Tuple>([]<typename C>() {
        if constexpr (known_bug<C> == Bugs) {
            check_riccati<Discrete, C, double>();
            if constexpr (C::float_ok) {
                check_riccati<Discrete, C, float>();
            }
        }
    });
}

template<bool Discrete, typename Tuple, bool Bugs>
void check_riccati_fail_cases() {
    for_each_case<Tuple>([]<typename C>() {
        if constexpr (known_bug<C> == Bugs) {
            check_riccati_fails<Discrete, C, double>();
            check_riccati_fails<Discrete, C, float>();
        }
    });
}

// ---- Factorizations, eigenvalues, matrix functions, Lyapunov, LTI conversions ----

/// err ≤ limit · unit; the report gives err / unit, the margin in the check's own units.
inline void check_within(const char* what, double err, double unit, double limit) {
    INFO(what, ": ", err, " = ", (unit > 0.0 ? err / unit : err), " units (limit ", limit, ")");
    CHECK(err <= limit * unit);
}

/// Run a check on every case of a tuple in double, and in float unless the case is double only.
template<typename Tuple, typename F>
void run_cases(F&& f) {
    for_each_case<Tuple>([&]<typename C>() {
        f.template operator()<C, double>();
        if constexpr (C::float_ok) {
            f.template operator()<C, float>();
        }
    });
}

/// The same matrix in double (a copy when it already is).
template<size_t R, size_t C, typename T>
Matrix<R, C> dbl(const Matrix<R, C, T>& M) {
    if constexpr (std::is_same_v<T, double>) {
        return M;
    } else {
        return Matrix<R, C>(M);
    }
}

template<typename T>
damp::complex<double> dbl(const damp::complex<T>& z) {
    return damp::complex<double>{static_cast<double>(z.real()), static_cast<double>(z.imag())};
}

template<size_t R, size_t C>
Matrix<R, C> abs_of(const Matrix<R, C>& M) {
    Matrix<R, C> out;
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            out(i, j) = damp::abs(M(i, j));
        }
    }
    return out;
}

/// xᵐ for a small whole m (a root multiplicity or Jordan block size).
inline double ipow(double x, double m) {
    double p = 1.0;
    for (int k = 0; k < static_cast<int>(m); ++k) {
        p *= x;
    }
    return p;
}

/// Greedy match of computed roots to the first @p count reference roots, each checked as
/// |Δz|ᵐ ≤ limit · ne · scaleᵢ (scaleᵢ in units of zᵐ: an m-fold root moves as the m-th root
/// of the perturbation).
template<size_t N>
void check_roots(
    const char*                                 what,
    const std::array<damp::complex<double>, N>& got,
    size_t                                      count,
    const std::array<double, N>&                re,
    const std::array<double, N>&                im,
    const std::array<double, N>&                scale,
    const std::array<double, N>&                mult,
    double                                      ne,
    double                                      limit
) {
    std::array<bool, N> used{};
    for (size_t i = 0; i < count; ++i) {
        const damp::complex<double> ref{re[i], im[i]};
        size_t                      best = N;
        double                      dist = 0.0;
        for (size_t j = 0; j < count; ++j) {
            const double d = (got[j] - ref).abs();
            if (!used[j] && (best == N || d < dist)) {
                best = j;
                dist = d;
            }
        }
        if (best == N) {
            return;
        }
        used[best] = true;
        INFO(what, " ", i, ": reference (", re[i], ", ", im[i], "), got (", got[best].real(), ", ", got[best].imag(), ")");
        check_within(what, ipow(dist, mult[i]), ne * scale[i], limit);
    }
}

template<size_t N>
std::array<double, N> filled(double v) {
    std::array<double, N> a{};
    a.fill(v);
    return a;
}

template<size_t R, size_t C>
bool is_symmetric_data(const std::array<double, R * C>& a) {
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            if (a[(i * C) + j] != a[(j * C) + i]) {
                return false;
            }
        }
    }
    return true;
}

inline constexpr double C_FACTOR = 12.0; ///< backward error, orthogonality and residual constant
inline constexpr double C_FWD = 16.0;    ///< forward error against the reference, on top of its condition number
inline constexpr double C_ROOT = 4.0;    ///< eigenvalue, root and pole constant, on top of its sensitivity

// ---- LU, QR, Cholesky, LDLᵀ, SVD ----

template<typename C, typename T>
void check_lu() {
    constexpr size_t N = C::n;
    using MD = Matrix<N, N>;
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const double ne = static_cast<double>(N) * eps_of<T>();
    const auto   A = load<N, N, T>(C::A);
    const auto   f = mat::lu_decomposition(A);
    CHECK(f.has_value());
    if (!f) {
        return;
    }
    const auto& [Lt, Ut, piv] = f.value();
    const MD L = dbl(Lt);
    const MD Ad = dbl(A);
    MD       U = MD::zeros(); // the upper triangle, which lu_solve reads
    MD       PA;
    bool     shape = true;
    for (size_t i = 0; i < N; ++i) {
        shape = shape && piv[i] < N && L(i, i) == 1.0;
        for (size_t j = 0; j < N; ++j) {
            PA(i, j) = Ad(piv[i] < N ? piv[i] : 0, j);
            U(i, j) = j >= i ? static_cast<double>(Ut(i, j)) : 0.0;
            shape = shape && (j <= i || L(i, j) == 0.0) && (j == i || piv[i] != piv[j]);
        }
    }
    CHECK(shape);
    // Partial pivoting is backward stable relative to |L||U| (Higham 2002, Thm 9.3).
    const double lu_size = MD(abs_of(L) * abs_of(U)).norm();
    check_within("lu ‖PA − LU‖ / ‖|L||U|‖", MD(PA - L * U).norm(), ne * lu_size, C_FACTOR);
    // det and x move by ‖A⁻¹‖ times that backward error.
    const double cond = C::inv_norm * lu_size;
    double       det = 1.0;
    for (size_t i = 0; i < N; ++i) {
        det *= damp::abs(U(i, i));
    }
    check_within("lu |det A| relative error", damp::abs(det - C::absdet) / C::absdet, ne * cond, C_FWD);

    const auto x = mat::lu_solve(A, load<N, 1, T>(C::b));
    CHECK(x.has_value());
    if (!x) {
        return;
    }
    const Matrix<N, 1> xd = dbl(x.value());
    const Matrix<N, 1> b = dbl(load<N, 1, T>(C::b));
    check_within("lu_solve ‖b − Ax‖ / (‖|L||U|‖‖x‖)", Matrix<N, 1>(b - Ad * xd).norm(), ne * lu_size * xd.norm(), C_FACTOR);
    check_within("lu_solve forward error", rel_diff(xd, load<N, 1, double>(C::x)), ne * cond, C_FWD);
}

template<typename C, typename T>
void check_qr() {
    constexpr size_t R = C::rows;
    constexpr size_t K = C::cols;
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const double ne = static_cast<double>(R) * eps_of<T>();
    const auto   A = load<R, K, T>(C::A);
    const auto   qr = mat::qr_decompose(A);
    CHECK(qr.full_rank());
    const Matrix<R, K> Q = dbl(qr.Q);
    const Matrix<K, K> Rm = dbl(qr.R);
    const auto         di = recip(C::units);
    check_within("qr ‖QᵀQ − I‖", Matrix<K, K>(Q.t() * Q - Matrix<K, K>::identity()).norm(), ne, C_FACTOR);
    // Householder QR is backward stable column by column, ‖Δaⱼ‖ ≤ γ̃_{mn}‖aⱼ‖ (Higham 2002,
    // Thm 19.4): unit columns, m·n·ε per column.
    const auto res = scale(ones<R>(), Matrix<R, K>(Q * Rm - dbl(A)), di);
    check_within("qr ‖(QR − A)D⁻¹‖, unit columns", res.norm(), ne * static_cast<double>(K) * damp::sqrt(static_cast<double>(K)), C_FACTOR);
    bool triangular = true;
    for (size_t i = 0; i < K; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            triangular = triangular && (i == j ? Rm(i, i) >= 0.0 : Rm(i, j) == 0.0);
        }
    }
    CHECK(triangular);
    const double err = rel_diff(scale(ones<K>(), Rm, di), scale(ones<K>(), load<K, K, double>(C::R), di));
    check_within("qr R forward error, unit columns", err, ne * damp::max(1.0, C::kappa), C_FWD);
}

template<typename C, typename T>
void check_cholesky() {
    constexpr size_t N = C::n;
    using MD = Matrix<N, N>;
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const double ne = static_cast<double>(N) * eps_of<T>();
    const auto   A = load<N, N, T>(C::A);
    const auto   Lo = mat::cholesky(A);
    CHECK(Lo.has_value());
    if (!Lo) {
        return;
    }
    const MD   L = dbl(Lo.value());
    const MD   Ad = dbl(A);
    const auto d = C::units;
    const auto di = recip(d);
    bool       triangular = true;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i; j < N; ++j) {
            triangular = triangular && (i == j ? L(i, i) > 0.0 : L(i, j) == 0.0);
        }
    }
    CHECK(triangular);
    // |ΔA| ≤ γ|L||Lᵀ| and |L||Lᵀ| ≤ √(aᵢᵢaⱼⱼ): the error is small in natural units D⁻¹AD⁻¹.
    check_within("cholesky ‖D⁻¹(LLᵀ − A)D⁻¹‖", scale(di, MD(L * L.t() - Ad), di).norm(), ne, C_FACTOR);
    const double err = rel_diff(scale(di, L, ones<N>()), scale(di, load<N, N, double>(C::L), ones<N>()));
    check_within("cholesky L forward error, natural units", err, ne * damp::max(1.0, C::kappa), C_FWD);

    const auto x = mat::cholesky_solve(A, load<N, 1, T>(C::b));
    CHECK(x.has_value());
    if (!x) {
        return;
    }
    const Matrix<N, 1> xd = dbl(x.value());
    const Matrix<N, 1> b = dbl(load<N, 1, T>(C::b));
    const double       den = (scale(di, Ad, di).norm() * scale(d, xd, ones<1>()).norm()) + scale(di, b, ones<1>()).norm();
    check_within("cholesky_solve residual, natural units", scale(di, Matrix<N, 1>(Ad * xd - b), ones<1>()).norm(), ne * den, C_FACTOR);
    const auto xr = load<N, 1, double>(C::x);
    check_within("cholesky_solve forward error, natural units", rel_diff(scale(d, xd, ones<1>()), scale(d, xr, ones<1>())), ne * damp::max(1.0, C::kappa), C_FWD);
}

template<typename C, typename T>
void check_cholesky_fails() {
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const auto A = load<C::n, C::n, T>(C::A);
    CHECK_FALSE(mat::cholesky(A).has_value());
    Matrix<C::n, 1, T> b{};
    b(0, 0) = T{1};
    CHECK_FALSE(mat::cholesky_solve(A, b).has_value());
}

template<typename C, typename T>
void check_ldlt() {
    constexpr size_t N = C::n;
    using MD = Matrix<N, N>;
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const double ne = static_cast<double>(N) * eps_of<T>();
    const auto   A = load<N, N, T>(C::A);
    const auto   f = mat::ldlt(A);
    CHECK(f.has_value());
    if (!f) {
        return;
    }
    const auto d = C::units;
    const auto di = recip(d);
    // Natural units A₀ = D⁻¹AD⁻¹: L₀ = D⁻¹LD, d₀ = d / D².
    const MD A0 = scale(di, dbl(A), di);
    const MD L0 = scale(di, dbl(f->L), d);
    MD       D0 = MD::zeros();
    size_t   neg = 0;
    bool     unit_lower = true;
    for (size_t i = 0; i < N; ++i) {
        D0(i, i) = static_cast<double>(f->d[i]) * di[i] * di[i];
        neg += (f->d[i] < 0) ? 1U : 0U;
        for (size_t j = i; j < N; ++j) {
            unit_lower = unit_lower && f->L(i, j) == (i == j ? T{1} : T{0});
        }
    }
    CHECK(unit_lower);
    CHECK(neg == C::neg);
    // Unpivoted LDLᵀ is backward stable relative to |L||D||Lᵀ| (Higham 2002, §11.1).
    const double growth = MD(abs_of(L0) * abs_of(D0) * abs_of(L0).t()).norm();
    check_within("ldlt ‖LDLᵀ − A‖ / ‖|L||D||Lᵀ|‖", MD(L0 * D0 * L0.t() - A0).norm(), ne * growth, C_FACTOR);
    const double cond = damp::max(1.0, C::kappa) * growth / A0.norm();
    const MD     Lref0 = scale(di, load<N, N, double>(C::L), d);
    check_within("ldlt L forward error", rel_diff(L0, Lref0), ne * cond, C_FWD);
    for (size_t i = 0; i < N; ++i) {
        const double ref = C::d[i] * di[i] * di[i];
        check_within("ldlt d forward error", damp::abs(D0(i, i) - ref) / damp::abs(ref), ne * cond, C_FWD);
    }

    const auto x = mat::ldlt_solve(A, load<N, 1, T>(C::b));
    CHECK(x.has_value());
    if (!x) {
        return;
    }
    const Matrix<N, 1> xd = dbl(x.value());
    const Matrix<N, 1> b = dbl(load<N, 1, T>(C::b));
    const double       den = (growth * scale(d, xd, ones<1>()).norm()) + scale(di, b, ones<1>()).norm();
    check_within("ldlt_solve residual, natural units", scale(di, Matrix<N, 1>(dbl(A) * xd - b), ones<1>()).norm(), ne * den, C_FACTOR);
    const auto xr = load<N, 1, double>(C::x);
    check_within("ldlt_solve forward error, natural units", rel_diff(scale(d, xd, ones<1>()), scale(d, xr, ones<1>())), ne * cond, C_FWD);
}

template<typename C, typename T>
void check_ldlt_fails() {
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    CHECK_FALSE(mat::ldlt(load<C::n, C::n, T>(C::A)).has_value());
}

template<typename C, typename T>
void check_svd() {
    constexpr size_t M = C::rows;
    constexpr size_t N = C::cols;
    constexpr size_t K = M < N ? M : N;
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const double ne = static_cast<double>(M > N ? M : N) * eps_of<T>();
    const auto   A = load<M, N, T>(C::A);
    const auto   r = mat::svd(A);
    CHECK(r.success);
    const Matrix<M, M> U = dbl(r.singular_U);
    const Matrix<N, N> V = dbl(r.singular_V);
    Matrix<M, N>       S = Matrix<M, N>::zeros();
    bool               sorted = true;
    for (size_t k = 0; k < K; ++k) {
        S(k, k) = static_cast<double>(r.singular_values[k]);
        sorted = sorted && S(k, k) >= 0.0 && (k == 0 || S(k, k) <= S(k - 1, k - 1));
    }
    CHECK(sorted);
    check_within("svd ‖UᵀU − I‖", Matrix<M, M>(U.t() * U - Matrix<M, M>::identity()).norm(), ne, C_FACTOR);
    check_within("svd ‖VᵀV − I‖", Matrix<N, N>(V.t() * V - Matrix<N, N>::identity()).norm(), ne, C_FACTOR);
    const Matrix<M, N> Ad = dbl(A);
    check_within("svd ‖UΣVᵀ − A‖ / ‖A‖", Matrix<M, N>(U * S * V.t() - Ad).norm(), ne * Ad.norm(), C_FACTOR);
    // Weyl: a backward error E moves every singular value by at most ‖E‖₂.
    for (size_t k = 0; k < K; ++k) {
        check_within("svd σ error / σ_max", damp::abs(S(k, k) - C::s[k]), ne * C::s[0], C_FACTOR);
    }
}

// ---- Eigenvalues and the real Schur basis ----

template<typename C, typename T>
void check_eig() {
    constexpr size_t N = C::n;
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const double ne = static_cast<double>(N) * eps_of<T>();
    const auto   r = mat::compute_eigenvalues(load<N, N, T>(C::A));
    CHECK(r.success);
    const auto                           A0 = scale(recip(C::units), dbl(load<N, N, T>(C::A)), C::units);
    const double                         nA = A0.norm();
    std::array<damp::complex<double>, N> got{};
    std::array<double, N>                sc{};
    damp::complex<double>                sum{0.0, 0.0};
    for (size_t i = 0; i < N; ++i) {
        got[i] = dbl(r.values[i]);
        sum = sum + got[i];
        sc[i] = C::ev_cond[i] * ipow(nA, C::ev_order[i]);
    }
    check_roots("eigenvalue |Δλ|ᵏ / (cond·‖A‖ᵏ)", got, N, C::ev_re, C::ev_im, sc, C::ev_order, ne, C_ROOT);
    // The sum of a cluster is as well conditioned as the trace, whatever the Jordan structure.
    check_within("eigenvalue Σλ − tr A", (sum - damp::complex<double>{A0.trace(), 0.0}).abs(), ne * nA, C_FACTOR);
    for (size_t i = 0; i < N; ++i) {
        if (got[i].imag() != 0.0) {
            double nearest = -1.0;
            for (size_t j = 0; j < N; ++j) {
                const double dj = (got[j] - damp::conj(got[i])).abs();
                nearest = (j != i && (nearest < 0.0 || dj < nearest)) ? dj : nearest;
            }
            check_within("eigenvalue conjugate partner", nearest, ne * nA, C_FACTOR);
        }
    }
}

/// compute_eigenvalues' schur_vectors: the real Schur basis Z of the balanced matrix A_b = D⁻¹AD,
/// mapped back as x = D·z with unit columns. Undo the map, then Zᵀ·A_b·Z must be quasi upper
/// triangular with Z orthogonal, its 1×1 and 2×2 diagonal blocks carrying the eigenvalues.
template<typename C, typename T>
void check_schur() {
    constexpr size_t N = C::n;
    using MD = Matrix<N, N>;
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const double ne = static_cast<double>(N) * eps_of<T>();
    const auto   A = load<N, N, T>(C::A);
    const auto   r = mat::compute_eigenvalues(A);
    const auto   e = mat::detail::balance_exponents(A);
    const MD     Ab = dbl(mat::detail::balanced(A, e));
    MD           Z;
    for (size_t i = 0; i < N; ++i) {
        double nrm = 0.0;
        for (size_t j = 0; j < N; ++j) {
            Z(j, i) = mat::detail::scale_by_pow2(static_cast<double>(r.schur_vectors(j, i)), -e[j]);
            nrm += Z(j, i) * Z(j, i);
        }
        for (size_t j = 0; j < N; ++j) {
            Z(j, i) /= damp::sqrt(nrm);
        }
    }
    const double nA = Ab.norm();
    check_within("schur ‖ZᵀZ − I‖", MD(Z.t() * Z - MD::identity()).norm(), ne, C_FACTOR);
    const MD T_ = Z.t() * Ab * Z;
    double   below = 0.0;
    for (size_t i = 0; i < N; ++i) {
        const damp::complex<double> lam = dbl(r.values[i]);
        const bool                  pair = lam.imag() != 0.0 && i + 1 < N;
        for (size_t row = i + (pair ? 2 : 1); row < N; ++row) {
            below += T_(row, i) * T_(row, i);
        }
        if (pair) {
            for (size_t row = i + 2; row < N; ++row) {
                below += T_(row, i + 1) * T_(row, i + 1);
            }
            const double tr = T_(i, i) + T_(i + 1, i + 1);
            const double det = (T_(i, i) * T_(i + 1, i + 1)) - (T_(i, i + 1) * T_(i + 1, i));
            check_within("schur 2×2 block trace − 2 Re λ", damp::abs(tr - (2.0 * lam.real())), ne * nA, C_FACTOR);
            check_within("schur 2×2 block det − |λ|²", damp::abs(det - (lam.abs() * lam.abs())), ne * nA * nA, C_FACTOR);
            ++i;
        } else {
            check_within("schur diagonal − λ", damp::abs(T_(i, i) - lam.real()), ne * nA, C_FACTOR);
        }
    }
    check_within("schur part below the quasi-triangle", damp::sqrt(below), ne * nA, C_FACTOR);
}

// ---- logm, sqrtm ----

template<typename C, typename T>
void check_logm() {
    constexpr size_t N = C::n;
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const double ne = static_cast<double>(N) * eps_of<T>();
    const auto   X = mat::log(load<N, N, T>(C::A));
    CHECK(X.has_value());
    if (!X) {
        return;
    }
    const auto di = recip(C::units);
    const auto X0 = scale(di, dbl(X.value()), C::units);
    const auto R0 = scale(di, load<N, N, double>(C::X), C::units);
    check_within("logm forward error", rel_diff(X0, R0), ne * damp::max(1.0, C::kappa), C_FWD);
}

template<typename C, typename T>
void check_sqrtm() {
    constexpr size_t N = C::n;
    using MD = Matrix<N, N>;
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const double ne = static_cast<double>(N) * eps_of<T>();
    const auto   X = mat::sqrt(load<N, N, T>(C::A));
    CHECK(X.has_value());
    if (!X) {
        return;
    }
    const auto di = recip(C::units);
    const MD   A0 = scale(di, dbl(load<N, N, T>(C::A)), C::units);
    const MD   X0 = scale(di, dbl(X.value()), C::units);
    // The Schur method is backward stable (Higham 2008, §6.2): the residual reaches α·ε with
    // α = ‖X‖²/‖A‖. The bound keeps the κ factor of a forward-stable method, α·κ·ε.
    check_within("sqrtm ‖X² − A‖ / ‖A‖", MD(X0 * X0 - A0).norm() / A0.norm(), ne * C::alpha * damp::max(1.0, C::kappa), C_FACTOR);
    const auto R0 = scale(di, load<N, N, double>(C::X), C::units);
    check_within("sqrtm forward error", rel_diff(X0, R0), ne * damp::max(1.0, C::kappa), C_FWD);
}

// ---- Lyapunov ----

template<bool Discrete, typename C, typename T>
auto run_lyapunov() {
    const auto A = load<C::n, C::n, T>(C::A);
    const auto Q = load<C::n, C::n, T>(C::Q);
    if constexpr (Discrete) {
        return dlyap(A, Q);
    } else {
        return lyap(A, Q);
    }
}

template<bool Discrete, typename C, typename T>
void check_lyapunov() {
    constexpr size_t N = C::n;
    using MD = Matrix<N, N>;
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const double ne = static_cast<double>(N) * eps_of<T>();
    const auto   X = run_lyapunov<Discrete, C, T>();
    CHECK(X.has_value());
    if (!X) {
        return;
    }
    // Natural units: A₀ = D⁻¹AD, Q₀ = D⁻¹QD⁻¹, X₀ = D⁻¹XD⁻¹.
    const auto d = C::units;
    const auto di = recip(d);
    const MD   A0 = scale(di, dbl(load<N, N, T>(C::A)), d);
    const MD   Q0 = scale(di, dbl(load<N, N, T>(C::Q)), di);
    const MD   X0 = scale(di, dbl(X.value()), di);
    // Normwise backward error of a Sylvester-type equation (Higham 2002, §16.2).
    MD     res;
    double den = 0.0;
    if constexpr (Discrete) {
        res = A0 * X0 * A0.t() - X0 + Q0;
        den = (A0.norm() * A0.norm() * X0.norm()) + X0.norm() + Q0.norm();
    } else {
        res = A0 * X0 + X0 * A0.t() + Q0;
        den = (2.0 * A0.norm() * X0.norm()) + Q0.norm();
    }
    check_within(Discrete ? "dlyap residual" : "lyap residual", res.norm(), ne * den, C_FACTOR);
    if (is_symmetric_data<N, N>(C::Q)) {
        check_within("lyapunov ‖X − Xᵀ‖ / ‖X‖", MD(X0 - X0.t()).norm(), ne * X0.norm(), C_FACTOR);
    }
    const MD R0 = scale(di, load<N, N, double>(C::X), di);
    check_within(Discrete ? "dlyap forward error" : "lyap forward error", rel_diff(X0, R0), ne * damp::max(1.0, C::kappa), C_FWD);
}

template<bool Discrete, typename C, typename T>
void check_lyapunov_fails() {
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    CHECK_FALSE(run_lyapunov<Discrete, C, T>().has_value());
}

// ---- discretize ----

template<typename C, typename T>
void check_zoh() {
    constexpr size_t NX = C::nx;
    constexpr size_t NU = C::nu;
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    StateSpace<NX, NU, 1, T> sys{};
    sys.A = load<NX, NX, T>(C::A);
    sys.B = load<NX, NU, T>(C::B);
    const auto sd = discretize(sys, static_cast<T>(C::Ts), DiscretizationMethod::ZOH);
    CHECK(sd.has_value());
    if (!sd) {
        return;
    }
    CHECK(sd->Ts == static_cast<T>(C::Ts));
    // Natural units: A_d₀ = D⁻¹A_dD, B_d₀ = D⁻¹B_d; compare the [A_d B_d] block of exp(M·Ts).
    const auto          di = recip(C::units);
    const auto          Ad0 = scale(di, dbl(sd->A), C::units);
    const auto          Bd0 = scale(di, dbl(sd->B), ones<NU>());
    const auto          Ar0 = scale(di, load<NX, NX, double>(C::Ad), C::units);
    const auto          Br0 = scale(di, load<NX, NU, double>(C::Bd), ones<NU>());
    Matrix<NX, NX + NU> got;
    Matrix<NX, NX + NU> ref;
    for (size_t i = 0; i < NX; ++i) {
        for (size_t j = 0; j < NX + NU; ++j) {
            got(i, j) = j < NX ? Ad0(i, j) : Bd0(i, j - NX);
            ref(i, j) = j < NX ? Ar0(i, j) : Br0(i, j - NX);
        }
    }
    const double ne = static_cast<double>(NX + NU) * eps_of<T>();
    check_within("zoh [Ad Bd] forward error", rel_diff(got, ref), ne * damp::max(1.0, C::kappa), C_EXPM);
}

template<typename C, typename T>
StateSpace<C::nx, C::nu, C::ny, T> tustin_plant() {
    StateSpace<C::nx, C::nu, C::ny, T> sys{};
    sys.A = load<C::nx, C::nx, T>(C::A);
    sys.B = load<C::nx, C::nu, T>(C::B);
    sys.C = load<C::ny, C::nx, T>(C::C);
    sys.D = load<C::ny, C::nu, T>(C::D);
    return sys;
}

template<typename C, typename T>
void check_tustin() {
    constexpr size_t NX = C::nx;
    constexpr size_t NU = C::nu;
    constexpr size_t NY = C::ny;
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const auto sd = discretize(tustin_plant<C, T>(), static_cast<T>(C::Ts), DiscretizationMethod::Tustin);
    CHECK(sd.has_value());
    if (!sd) {
        return;
    }
    // Natural units: x = Dx₀, so A₀ = D⁻¹AD, B₀ = D⁻¹B, C₀ = CD; compare [[A_d B_d], [C_d D_d]].
    const auto               d = C::units;
    const auto               di = recip(d);
    const auto               Ad0 = scale(di, dbl(sd->A), d);
    const auto               Bd0 = scale(di, dbl(sd->B), ones<NU>());
    const auto               Cd0 = scale(ones<NY>(), dbl(sd->C), d);
    const auto               Dd0 = dbl(sd->D);
    const auto               Ar0 = scale(di, load<NX, NX, double>(C::Ad), d);
    const auto               Br0 = scale(di, load<NX, NU, double>(C::Bd), ones<NU>());
    const auto               Cr0 = scale(ones<NY>(), load<NY, NX, double>(C::Cd), d);
    const auto               Dr0 = load<NY, NU, double>(C::Dd);
    Matrix<NX + NY, NX + NU> got;
    Matrix<NX + NY, NX + NU> ref;
    for (size_t i = 0; i < NX + NY; ++i) {
        for (size_t j = 0; j < NX + NU; ++j) {
            const bool top = i < NX;
            const bool left = j < NX;
            got(i, j) = top ? (left ? Ad0(i, j) : Bd0(i, j - NX)) : (left ? Cd0(i - NX, j) : Dd0(i - NX, j - NX));
            ref(i, j) = top ? (left ? Ar0(i, j) : Br0(i, j - NX)) : (left ? Cr0(i - NX, j) : Dr0(i - NX, j - NX));
        }
    }
    const double ne = static_cast<double>(NX + NU) * eps_of<T>();
    check_within("tustin forward error", rel_diff(got, ref), ne * damp::max(1.0, C::kappa), C_FWD);
}

template<typename C, typename T>
void check_tustin_fails() {
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    CHECK_FALSE(discretize(tustin_plant<C, T>(), static_cast<T>(C::Ts), DiscretizationMethod::Tustin).has_value());
}

// ---- SISO conversions ----

template<typename C, typename T>
StateSpace<C::nx, 1, 1, T> siso_plant() {
    StateSpace<C::nx, 1, 1, T> sys{};
    sys.A = load<C::nx, C::nx, T>(C::A);
    sys.B = load<C::nx, 1, T>(C::B);
    sys.C = load<1, C::nx, T>(C::C);
    sys.D(0, 0) = static_cast<T>(C::D);
    return sys;
}

template<typename C, typename T>
void check_ss2tf() {
    constexpr size_t NX = C::nx;
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const double ne = static_cast<double>(NX) * eps_of<T>();
    const auto   tf = to_transfer_function(siso_plant<C, T>());
    for (size_t k = 0; k <= NX; ++k) {
        INFO("coefficient of s^", NX - k);
        check_within("ss2tf den coefficient", damp::abs(static_cast<double>(tf.den[k]) - C::den[k]), ne * C::den_scale[k], C_FACTOR);
        check_within("ss2tf num coefficient", damp::abs(static_cast<double>(tf.num[k]) - C::num[k]), ne * C::num_scale[k], C_FACTOR);
    }
}

template<typename C, typename T>
void check_ss2zpk() {
    constexpr size_t NX = C::nx;
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const double ne = static_cast<double>(NX) * eps_of<T>();
    const auto   z = to_zpk(siso_plant<C, T>());
    CHECK(z.success);
    CHECK(z.n_zeros == C::n_zeros);
    std::array<damp::complex<double>, NX> poles{};
    std::array<damp::complex<double>, NX> zeros{};
    std::array<double, NX>                sc{};
    for (size_t i = 0; i < NX; ++i) {
        poles[i] = dbl(z.poles[i]);
        zeros[i] = dbl(z.zeros[i]);
        sc[i] = C::p_cond[i] * ipow(C::norm_a, C::p_order[i]);
    }
    check_roots("ss2zpk pole", poles, NX, C::p_re, C::p_im, sc, C::p_order, ne, C_ROOT);
    if (z.n_zeros == C::n_zeros) {
        check_roots("ss2zpk zero", zeros, C::n_zeros, C::z_re, C::z_im, C::z_scale, filled<NX>(1.0), ne, C_ROOT);
    }
    const double gscale = C::D != 0.0 ? damp::abs(C::D) : C::num_scale[NX - C::n_zeros];
    check_within("ss2zpk gain", damp::abs(static_cast<double>(z.gain) - C::gain), ne * gscale, C_FACTOR);
}

template<typename C, typename T>
TransferFunction<C::nnum, C::nden, T> load_tf() {
    TransferFunction<C::nnum, C::nden, T> tf{};
    for (size_t i = 0; i < C::nnum; ++i) {
        tf.num[i] = static_cast<T>(C::num[i]);
    }
    for (size_t i = 0; i < C::nden; ++i) {
        tf.den[i] = static_cast<T>(C::den[i]);
    }
    return tf;
}

template<typename C, typename T>
void check_tf2zpk() {
    constexpr size_t NZ = C::nnum - 1;
    constexpr size_t NP = C::nden - 1;
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const double ne = static_cast<double>(NP) * eps_of<T>();
    const auto   r = to_zpk(load_tf<C, T>());
    CHECK(r.success);
    std::array<damp::complex<double>, NZ> zeros{};
    std::array<damp::complex<double>, NP> poles{};
    for (size_t i = 0; i < NZ; ++i) {
        zeros[i] = dbl(r.zpk.zeros[i]);
    }
    for (size_t i = 0; i < NP; ++i) {
        poles[i] = dbl(r.zpk.poles[i]);
    }
    check_roots("tf2zpk zero |Δz|ᵐ", zeros, NZ, C::z_re, C::z_im, C::z_scale, C::z_mult, ne, C_ROOT);
    check_roots("tf2zpk pole |Δp|ᵐ", poles, NP, C::p_re, C::p_im, C::p_scale, C::p_mult, ne, C_ROOT);
    check_within("tf2zpk gain", damp::abs(static_cast<double>(r.zpk.gain) - C::gain), ne * damp::abs(C::gain), C_FACTOR);
}

template<typename C, typename T>
void check_zpk2tf() {
    constexpr size_t NZ = C::nz;
    constexpr size_t NP = C::np;
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const double   ne = static_cast<double>(NP) * eps_of<T>();
    ZPK<NZ, NP, T> m{};
    for (size_t i = 0; i < NZ; ++i) {
        m.zeros[i] = damp::complex<T>{static_cast<T>(C::z_re[i]), static_cast<T>(C::z_im[i])};
    }
    for (size_t i = 0; i < NP; ++i) {
        m.poles[i] = damp::complex<T>{static_cast<T>(C::p_re[i]), static_cast<T>(C::p_im[i])};
    }
    m.gain = static_cast<T>(C::gain);
    const auto tf = m.to_transfer_function();
    for (size_t k = 0; k <= NZ; ++k) {
        INFO("numerator coefficient ", k);
        check_within("zpk2tf num coefficient", damp::abs(static_cast<double>(tf.num[k]) - C::num[k]), ne * C::num_scale[k], C_FACTOR);
    }
    for (size_t k = 0; k <= NP; ++k) {
        INFO("denominator coefficient ", k);
        check_within("zpk2tf den coefficient", damp::abs(static_cast<double>(tf.den[k]) - C::den[k]), ne * C::den_scale[k], C_FACTOR);
    }
}

/// tf2ss: the realization is not unique, its Markov parameters h₀ = D, hₖ = C·Aᵏ⁻¹·B are.
template<typename C, typename T>
void check_tf2ss() {
    constexpr size_t NX = C::nden - 1;
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const auto ss = load_tf<C, T>().to_state_space();
    CHECK(ss.has_value());
    if (!ss) {
        return;
    }
    const Matrix<NX, NX> A = dbl(ss->A);
    const Matrix<NX, 1>  B = dbl(ss->B);
    const Matrix<1, NX>  Cm = dbl(ss->C);
    Matrix<NX, 1>        AkB = B;
    Matrix<NX, 1>        AkB_abs = abs_of(B);
    for (size_t k = 0; k < 2 * NX + 1; ++k) {
        // Each product rounds once per factor: |Δhₖ| ≲ (k + 1)·ε·|C||A|ᵏ⁻¹|B|.
        const double h = k == 0 ? static_cast<double>(ss->D(0, 0)) : (Cm * AkB)(0, 0);
        const double sc = k == 0 ? damp::abs(h) : (abs_of(Cm) * AkB_abs)(0, 0);
        INFO("Markov parameter h", k);
        const double ne = static_cast<double>(NX * (k + 1)) * eps_of<T>();
        check_within("tf2ss Markov parameter", damp::abs(h - C::markov[k]), ne * sc, C_FACTOR);
        if (k > 0) {
            AkB = A * AkB;
            AkB_abs = abs_of(A) * AkB_abs;
        }
    }
}

// ---- Pole placement ----

template<typename C, typename T>
void check_place() {
    constexpr size_t NX = C::nx;
    constexpr size_t NU = C::nu;
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const double             ne = static_cast<double>(NX) * eps_of<T>();
    StateSpace<NX, NU, 1, T> sys{};
    sys.A = load<NX, NX, T>(C::A);
    sys.B = load<NX, NU, T>(C::B);
    sys.Ts = static_cast<T>(C::Ts);
    damp::array<damp::complex<T>, NX> poles{};
    for (size_t i = 0; i < NX; ++i) {
        poles[i] = damp::complex<T>{static_cast<T>(C::p_re[i]), static_cast<T>(C::p_im[i])};
    }
    const auto r = design::place(sys, poles);
    CHECK(r.success);
    if (!r.success) {
        return;
    }
    const Matrix<NX, NX> A = dbl(sys.A);
    const Matrix<NX, NU> B = dbl(sys.B);
    const Matrix<NU, NX> K = dbl(r.K);
    const auto           ev = mat::compute_eigenvalues(Matrix<NX, NX>(A - B * K));
    CHECK(ev.success);
    // Bauer–Fike on A − BK with the backward error of K: ε·(‖A‖ + ‖B‖‖K‖).
    const double                          size = A.norm() + (B.norm() * K.norm());
    std::array<damp::complex<double>, NX> got{};
    std::array<double, NX>                sc{};
    for (size_t i = 0; i < NX; ++i) {
        got[i] = ev.values[i];
        sc[i] = C::p_cond[i] * ipow(size, C::p_order[i]);
    }
    check_roots("place closed-loop pole", got, NX, C::p_re, C::p_im, sc, C::p_order, ne, C_ROOT);
}

// ---- Generalized eigenvalues (QZ) and invariant zeros ----

/// σ_min of the complex matrix X + iY, from its real form [X −Y; Y X] (every singular value twice).
template<size_t R, size_t K>
double smallest_singular_value(const Matrix<R, K>& X, const Matrix<R, K>& Y) {
    Matrix<2 * R, 2 * K> M{};
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < K; ++j) {
            M(i, j) = X(i, j);
            M(i, K + j) = -Y(i, j);
            M(R + i, j) = Y(i, j);
            M(R + i, K + j) = X(i, j);
        }
    }
    const auto s = mat::svd(M);
    return s.singular_values[(R < K ? 2 * R : 2 * K) - 1];
}

/// Finite generalized eigenvalues: the count, each within C_ROOT of its condition, and the backward
/// error σ_min(A − λB)/(‖A‖ + |λ|‖B‖) of each.
template<typename C, typename T>
void check_qz() {
    constexpr size_t N = C::n;
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const double ne = static_cast<double>(N) * eps_of<T>();
    const auto   A = load<N, N, T>(C::A);
    const auto   B = load<N, N, T>(C::B);
    const auto   r = mat::generalized_eigenvalues(A, B);
    CHECK(r.success);
    CHECK(r.n_finite == C::n_finite);
    std::array<damp::complex<double>, N> got{};
    for (size_t i = 0; i < N; ++i) {
        got[i] = dbl(r.values[i]);
    }
    if (r.n_finite == C::n_finite) {
        check_roots("qz eigenvalue |Δλ| / cond", got, C::n_finite, C::ev_re, C::ev_im, C::ev_cond, filled<N>(1.0), ne, C_ROOT);
    }
    const Matrix<N, N> Ad = dbl(A);
    const Matrix<N, N> Bd = dbl(B);
    for (size_t i = 0; i < r.n_finite; ++i) {
        const double smin = smallest_singular_value(Matrix<N, N>(Ad - (got[i].real() * Bd)), Matrix<N, N>(Bd * -got[i].imag()));
        check_within("qz backward error σ_min(A − λB) / (‖A‖ + |λ|‖B‖)", smin / (Ad.norm() + (got[i].abs() * Bd.norm())), ne, C_FACTOR);
    }
}

template<typename C, typename T>
StateSpace<C::nx, C::nu, C::ny, T> mimo_plant() {
    StateSpace<C::nx, C::nu, C::ny, T> sys{};
    sys.A = load<C::nx, C::nx, T>(C::A);
    sys.B = load<C::nx, C::nu, T>(C::B);
    sys.C = load<C::ny, C::nx, T>(C::C);
    sys.D = load<C::ny, C::nu, T>(C::D);
    return sys;
}

/// Invariant zeros: the count, each within C_ROOT of its condition, and the backward error
/// σ_min(S₀ − z₀N₀)/‖S₀‖ of each in natural units (z₀ = z/α, S₀ the natural system matrix).
template<typename C, typename T>
void check_tzero() {
    constexpr size_t NX = C::nx;
    constexpr size_t NU = C::nu;
    constexpr size_t NY = C::ny;
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const double ne = static_cast<double>(NX + (NU > NY ? NU : NY)) * eps_of<T>();
    const auto   sys = mimo_plant<C, T>();
    const auto   r = transmission_zeros(sys);
    CHECK(r.success);
    CHECK(r.n_zeros == C::n_zeros);
    std::array<damp::complex<double>, NX> got{};
    for (size_t i = 0; i < NX; ++i) {
        got[i] = dbl(r.zeros[i]);
    }
    if (r.n_zeros == C::n_zeros) {
        check_roots("tzero zero |Δz| / cond", got, C::n_zeros, C::z_re, C::z_im, C::z_scale, filled<NX>(1.0), ne, C_ROOT);
    }
    const auto               ix = recip(C::units_x);
    const auto               iu = recip(C::units_u);
    const auto               iy = recip(C::units_y);
    Matrix<NX + NY, NX + NU> S0{};
    Matrix<NX + NY, NX + NU> N0{};
    S0.template block<NX, NX>(0, 0) = scale(ix, dbl(sys.A), C::units_x) * (1.0 / C::alpha);
    S0.template block<NX, NU>(0, NX) = scale(ix, dbl(sys.B), iu) * (1.0 / C::alpha);
    S0.template block<NY, NX>(NX, 0) = scale(iy, dbl(sys.C), C::units_x);
    S0.template block<NY, NU>(NX, NX) = scale(iy, dbl(sys.D), iu);
    for (size_t i = 0; i < NX; ++i) {
        N0(i, i) = 1.0;
    }
    for (size_t i = 0; i < r.n_zeros; ++i) {
        const damp::complex<double> z0 = got[i] * (1.0 / C::alpha);
        const double                smin = smallest_singular_value(Matrix<NX + NY, NX + NU>(S0 - (z0.real() * N0)), Matrix<NX + NY, NX + NU>(N0 * -z0.imag()));
        check_within("tzero backward error σ_min(S₀ − z₀N₀) / ‖S₀‖", smin / S0.norm(), ne, C_FACTOR);
    }
}

/// transmission_zeros on the SISO conversion cases: the zeros of ss2zpk's exact reference.
template<typename C, typename T>
void check_ss2tf_tzero() {
    constexpr size_t NX = C::nx;
    INFO(doctest::String(C::name), " (", doctest::String(type_name<T>()), ")");
    const double ne = static_cast<double>(NX) * eps_of<T>();
    const auto   r = transmission_zeros(siso_plant<C, T>());
    CHECK(r.success);
    CHECK(r.n_zeros == C::n_zeros);
    std::array<damp::complex<double>, NX> got{};
    for (size_t i = 0; i < NX; ++i) {
        got[i] = dbl(r.zeros[i]);
    }
    if (r.n_zeros == C::n_zeros) {
        check_roots("tzero SISO zero", got, C::n_zeros, C::z_re, C::z_im, C::z_scale, filled<NX>(1.0), ne, C_ROOT);
    }
}

// ---- Fixed-seed properties ----

/// SplitMix64 (Steele, Lea & Flood, 2014): the same stream on every platform for a given seed.
struct SplitMix64 {
    std::uint64_t state{0};

    constexpr std::uint64_t next() {
        state += 0x9E3779B97F4A7C15ULL;
        std::uint64_t z = state;
        z = (z ^ (z >> 30U)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27U)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31U);
    }

    /// Uniform on [−1, 1).
    constexpr double uniform() { return (static_cast<double>(next() >> 11U) * 0x1.0p-52) - 1.0; }
};

inline constexpr std::uint64_t kSeeds = 16;

template<size_t R, size_t C, typename T>
Matrix<R, C, T> random_matrix(SplitMix64& g) {
    Matrix<R, C, T> M;
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            M(i, j) = static_cast<T>(g.uniform());
        }
    }
    return M;
}

template<size_t R, size_t K, typename T>
void qr_property() {
    const double ne = static_cast<double>(R) * eps_of<T>();
    for (std::uint64_t s = 1; s <= kSeeds; ++s) {
        INFO(R, "×", K, " ", doctest::String(type_name<T>()), " seed ", s);
        SplitMix64         g{s};
        const auto         A = random_matrix<R, K, T>(g);
        const auto         qr = mat::qr_decompose(A);
        const Matrix<R, K> Q = dbl(qr.Q);
        const Matrix<K, K> Rm = dbl(qr.R);
        const Matrix<R, K> Ad = dbl(A);
        check_within("property qr ‖QR − A‖ / ‖A‖", Matrix<R, K>(Q * Rm - Ad).norm(), ne * Ad.norm(), C_FACTOR);
        check_within("property qr ‖QᵀQ − I‖", Matrix<K, K>(Q.t() * Q - Matrix<K, K>::identity()).norm(), ne, C_FACTOR);
        bool triangular = true;
        for (size_t i = 0; i < K; ++i) {
            for (size_t j = 0; j <= i; ++j) {
                triangular = triangular && (i == j ? Rm(i, i) >= 0.0 : Rm(i, j) == 0.0);
            }
        }
        CHECK(triangular);
    }
}

template<size_t M, size_t N, typename T>
void svd_property() {
    constexpr size_t K = M < N ? M : N;
    const double     ne = static_cast<double>(M > N ? M : N) * eps_of<T>();
    for (std::uint64_t s = 1; s <= kSeeds; ++s) {
        INFO(M, "×", N, " ", doctest::String(type_name<T>()), " seed ", s);
        SplitMix64         g{s};
        const auto         A = random_matrix<M, N, T>(g);
        const auto         r = mat::svd(A);
        const Matrix<M, M> U = dbl(r.singular_U);
        const Matrix<N, N> V = dbl(r.singular_V);
        Matrix<M, N>       S = Matrix<M, N>::zeros();
        for (size_t k = 0; k < K; ++k) {
            S(k, k) = static_cast<double>(r.singular_values[k]);
        }
        const Matrix<M, N> Ad = dbl(A);
        CHECK(r.success);
        check_within("property svd ‖UΣVᵀ − A‖ / ‖A‖", Matrix<M, N>(U * S * V.t() - Ad).norm(), ne * Ad.norm(), C_FACTOR);
        check_within("property svd ‖UᵀU − I‖", Matrix<M, M>(U.t() * U - Matrix<M, M>::identity()).norm(), ne, C_FACTOR);
        check_within("property svd ‖VᵀV − I‖", Matrix<N, N>(V.t() * V - Matrix<N, N>::identity()).norm(), ne, C_FACTOR);
    }
}

template<size_t N, typename T>
void solve_property() {
    using MD = Matrix<N, N>;
    const double ne = static_cast<double>(N) * eps_of<T>();
    for (std::uint64_t s = 1; s <= kSeeds; ++s) {
        INFO(N, "×", N, " ", doctest::String(type_name<T>()), " seed ", s);
        SplitMix64 g{s};
        const auto A = random_matrix<N, N, T>(g);
        const auto b = random_matrix<N, 1, T>(g);
        const auto f = mat::lu_decomposition(A);
        const auto x = mat::lu_solve(A, b);
        CHECK(f.has_value());
        CHECK(x.has_value());
        if (f && x) {
            const double       lu_size = MD(abs_of(dbl(std::get<0>(*f))) * abs_of(dbl(std::get<1>(*f)))).norm();
            const Matrix<N, 1> xd = dbl(*x);
            check_within("property lu_solve ‖b − Ax‖ / (‖|L||U|‖‖x‖)", Matrix<N, 1>(dbl(b) - dbl(A) * xd).norm(), ne * lu_size * xd.norm(), C_FACTOR);
        }
        // An SPD matrix: G Gᵀ/N + I/10.
        const auto Gm = random_matrix<N, N, T>(g);
        const auto S = Matrix<N, N, T>(Gm * Gm.t() * static_cast<T>(1.0 / static_cast<double>(N)) + Matrix<N, N, T>::identity() * static_cast<T>(0.1));
        const auto y = mat::cholesky_solve(S, b);
        CHECK(y.has_value());
        if (y) {
            const Matrix<N, 1> yd = dbl(*y);
            const MD           Sd = dbl(S);
            const double       den = (Sd.norm() * yd.norm()) + dbl(b).norm();
            check_within("property cholesky_solve ‖b − Sy‖ / (‖S‖‖y‖ + ‖b‖)", Matrix<N, 1>(dbl(b) - Sd * yd).norm(), ne * den, C_FACTOR);
        }
    }
}

template<size_t N, typename T>
void lyapunov_property() {
    using MD = Matrix<N, N>;
    const double ne = static_cast<double>(N) * eps_of<T>();
    for (std::uint64_t s = 1; s <= kSeeds; ++s) {
        INFO(N, "×", N, " ", doctest::String(type_name<T>()), " seed ", s);
        SplitMix64 g{s};
        const auto M = random_matrix<N, N, T>(g);
        const auto Gm = random_matrix<N, 2, T>(g);
        const auto Q = Matrix<N, N, T>(Gm * Gm.t());
        // Hurwitz: every |λ(M)| ≤ ‖M‖_F. Schur stable: ρ(0.9·M/‖M‖_F) ≤ 0.9.
        const auto Ac = Matrix<N, N, T>(M - Matrix<N, N, T>::identity() * (M.norm() + static_cast<T>(0.1)));
        const auto Ad = Matrix<N, N, T>(M * (static_cast<T>(0.9) / M.norm()));
        const auto Xc = lyap(Ac, Q);
        const auto Xd = dlyap(Ad, Q);
        CHECK(Xc.has_value());
        CHECK(Xd.has_value());
        const MD Qd = dbl(Q);
        if (Xc) {
            const MD A = dbl(Ac);
            const MD X = dbl(*Xc);
            const MD res = A * X + X * A.t() + Qd;
            check_within("property lyap residual", res.norm(), ne * ((2.0 * A.norm() * X.norm()) + Qd.norm()), C_FACTOR);
        }
        if (Xd) {
            const MD A = dbl(Ad);
            const MD X = dbl(*Xd);
            const MD res = A * X * A.t() - X + Qd;
            check_within("property dlyap residual", res.norm(), ne * ((A.norm() * A.norm() * X.norm()) + X.norm() + Qd.norm()), C_FACTOR);
        }
    }
}

/// d2c ∘ c2d: resampling a ZOH model to a fractional ratio takes the principal logarithm of
/// the discrete block (the d2c step) and must land on the direct discretization.
template<size_t NX, size_t NU, typename T>
void resample_property() {
    const double ne = static_cast<double>(NX + NU) * eps_of<T>();
    for (std::uint64_t s = 1; s <= kSeeds; ++s) {
        INFO(NX, " states, ", NU, " inputs, ", doctest::String(type_name<T>()), " seed ", s);
        SplitMix64               g{s};
        StateSpace<NX, NU, 1, T> sys{};
        const auto               M = random_matrix<NX, NX, T>(g);
        const auto               Bm = random_matrix<NX, NU, T>(g);
        sys.A = M * (T{1} / M.norm());
        sys.B = Bm * (T{1} / Bm.norm());
        const auto sd = discretize(sys, T{1}, DiscretizationMethod::ZOH);
        const auto fast = discretize(sys, static_cast<T>(0.7), DiscretizationMethod::ZOH);
        CHECK(sd.has_value());
        CHECK(fast.has_value());
        if (!sd || !fast) {
            continue;
        }
        const auto rs = resample(*sd, static_cast<T>(0.7));
        CHECK(rs.has_value());
        if (!rs) {
            continue;
        }
        Matrix<NX, NX + NU> got;
        Matrix<NX, NX + NU> ref;
        for (size_t i = 0; i < NX; ++i) {
            for (size_t j = 0; j < NX + NU; ++j) {
                got(i, j) = static_cast<double>(j < NX ? rs->A(i, j) : rs->B(i, j - NX));
                ref(i, j) = static_cast<double>(j < NX ? fast->A(i, j) : fast->B(i, j - NX));
            }
        }
        check_within("property resample through log vs direct ZOH", rel_diff(got, ref), ne, C_FWD);
    }
}

} // namespace oracle_checks
