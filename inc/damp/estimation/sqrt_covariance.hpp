// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file sqrt_covariance.hpp
 * @brief The square-root covariance kernel of the recursive Kalman filters: a filter carries the
 *        lower-triangular factor S of its covariance, P = S Sᵀ, and never the covariance itself
 *
 * A covariance stored as P holds its eigenvalues to eps·‖P‖. Once cond(P) passes about 1/eps the
 * small ones are rounding noise before any update runs, whatever formula the update uses. A prior
 * much wider than a reading's noise does that on the first tick: on a state plus a sensor bias,
 * from P₀/R ≈ 10⁵ in float (10⁹ in double), the Joseph-form filter's P went indefinite within a
 * few ticks (the third at 10⁶) and all but a dozen of 20 000 readings were refused. The factor
 * has the square root of that condition number, so the same arithmetic holds twice the range:
 * the float filter on S runs that problem to P₀/R = 10¹² with the output error of the double one.
 *
 * | Step | Function | Method |
 * | --- | --- | --- |
 * | P₀, Q, R → factor | semidefinite_factor(), process_noise_factor(), reading_noise_factor() | Cholesky that takes a zero pivot; P₀ and Q once at construction |
 * | A P Aᵀ + G Q Gᵀ | sqrt_time_update() | Householder triangularization of [A S, G Q½], in S's own storage |
 * | (I − K C) P | sqrt_measurement_update() | Carlson's triangular update, one whitened scalar reading at a time |
 *
 * Every function is in place on fixed-size storage, with loops bounded by the template sizes, and
 * `constexpr`.
 *
 * @see "Discrete square root filtering: a survey of current techniques" (Kaminski, Bryson &
 *      Schmidt, 1971), IEEE TAC 16(6), https://doi.org/10.1109/TAC.1971.1099816
 * @see "Fast triangular formulation of the square root filter" (Carlson, 1973), AIAA Journal
 *      11(9), https://doi.org/10.2514/3.6907
 * @see "Factorization Methods for Discrete Sequential Estimation" (Bierman, 1977), ch. IV
 *      (Householder orthogonal transformations) and ch. V (sequential square root data processing)
 * @see "Kalman Filtering: Theory and Practice Using MATLAB" (Grewal & Andrews, 2015), §7.4–7.5
 */

#include "damp/backend.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/pattern_matrix.hpp"

namespace damp::detail {

/**
 * @brief The covariance P = S Sᵀ of a lower-triangular factor S, exactly symmetric
 *
 * Forming P rounds it to eps·‖P‖: read it for display and checks, and keep computing on S.
 *
 * @param S Lower-triangular factor (the strict upper triangle is not read)
 * @return P = S Sᵀ
 */
template<size_t N, typename T>
[[nodiscard]] constexpr Matrix<N, N, T> covariance_of_factor(const Matrix<N, N, T>& S) {
    Matrix<N, N, T> P{detail::unfilled};
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            T s{0};
            for (size_t k = 0; k <= j; ++k) {
                s += S(i, k) * S(j, k);
            }
            P(i, j) = s;
            P(j, i) = s;
        }
    }
    return P;
}

/**
 * @brief The lower-triangular factor S of a positive semidefinite P, P = S Sᵀ, with a zero column
 *        where a pivot vanishes
 *
 * The Cholesky recursion in state order, for a covariance that may be singular: all zeros, a state
 * known exactly (a delay register, whose row and column are zero), or states that are exact
 * combinations of earlier ones. mat::cholesky() refuses those, since a solve cannot use them; a
 * factor can.
 *
 * Pivot j is @f$ d_j = p_{jj} - \sum_{k<j} s_{jk}^2 @f$, and its rounding is
 * @f$ (j+1)\,\varepsilon\,(|p_{jj}| + \sum_{k<j} s_{jk}^2) @f$. A pivot above that is a column
 * like any other, however small against the rest of P: the test is relative to the state's own
 * variance, so it does not depend on the states' units. A pivot within that band, or negative by
 * no more than @f$ 64 N \varepsilon \max_i |p_{ii}| @f$ (the roundoff of whatever computed P), is
 * a zero: its column is zero, and the entries it leaves out have to vanish with it,
 * @f$ r_{ij}^2 \le 4\,\tau_j\,(|p_{ii}| + \tau) @f$ for the band @f$ \tau_j @f$ of the pivot. A
 * pivot more negative than that, a left-out entry larger than that, or an entry that is not
 * finite means P is not a covariance, and it is refused.
 *
 * There is no pivoting: the factor has to be triangular in the state's own order. Each pivot's
 * test is therefore against the states before it only, as in mat::cholesky().
 *
 * @param P Symmetric positive semidefinite matrix (only the lower triangle is read)
 * @param S The factor, lower triangular with a non-negative diagonal; written only on success
 * @return false if P has a non-finite entry or is indefinite by more than roundoff: S is unchanged
 *
 * @see "Matrix Computations" (Golub & Van Loan, 2013), §4.2
 * @see "Analysis of the Cholesky Decomposition of a Semi-definite Matrix" (Higham, 1990), in
 *      Reliable Numerical Computation: the factorization with complete pivoting, an order a
 *      triangular factor of the state cannot take
 */
template<size_t N, typename T>
[[nodiscard]] constexpr bool semidefinite_factor(const Matrix<N, N, T>& P, Matrix<N, N, T>& S) {
    T top{0};
    for (size_t i = 0; i < N; ++i) {
        if (!damp::isfinite(P(i, i))) {
            return false;
        }
        top = damp::max(top, damp::abs(P(i, i)));
    }
    const T         floor = static_cast<T>(N) * T{64} * eps<T>() * top; // roundoff of the computation that produced P
    Matrix<N, N, T> L{};
    for (size_t j = 0; j < N; ++j) {
        T d = P(j, j);
        T mag = damp::abs(P(j, j));
        for (size_t k = 0; k < j; ++k) {
            d -= L(j, k) * L(j, k);
            mag += L(j, k) * L(j, k);
        }
        if (!damp::isfinite(d)) {
            return false; // a NaN or ±Inf below the diagonal reaches its row's pivot
        }
        const T    own = static_cast<T>(j + 1) * eps<T>() * mag; // rounding of the pivot's own terms
        const T    band = damp::max(own, floor);
        const bool zero = !(d > own);
        if (zero && (d < -band)) {
            return false; // indefinite
        }
        const T ljj = zero ? T{0} : damp::sqrt(d);
        L(j, j) = ljj;
        for (size_t i = j + 1; i < N; ++i) {
            T r = P(i, j);
            for (size_t k = 0; k < j; ++k) {
                r -= L(i, k) * L(j, k);
            }
            if (!zero) {
                L(i, j) = r / ljj;
            } else if (!damp::isfinite(r) || ((r * r) > (T{4} * band * (damp::abs(P(i, i)) + floor)))) {
                return false; // a zero variance with a covariance: indefinite
            }
        }
    }
    S = L;
    return true;
}

/**
 * @brief The factor a filter starts from: of @p P when it is a covariance; otherwise refused under
 *        constant evaluation, and the identity at run time
 *
 * The same rule as a failed design handed to a runtime: it does not compile where it is
 * constant-evaluated, and at run time the filter starts from P = I.
 */
template<size_t N, typename T>
[[nodiscard]] constexpr Matrix<N, N, T> covariance_factor_or_identity(const Matrix<N, N, T>& P) {
    Matrix<N, N, T> S{};
    if (semidefinite_factor(P, S)) {
        return S;
    }
    if (damp::is_constant_evaluated()) {
        refuse_failed_design(design::DesignStatus::bad_weights); // constinit / constexpr: a P₀ that is not a covariance does not compile
    }
    return Matrix<N, N, T>::identity();
}

/**
 * @brief Apply row @p i's reflection of triangularize() to the rows below it: row k ← row k − s·v,
 *        with column i's sign turned when @p flip
 *
 * A row the reflection does not move (its product with v is exactly zero: a state known exactly,
 * a delay register) is left as it is.
 */
template<size_t N, size_t M, typename T>
constexpr void reflect_rows_below(Matrix<N, N, T>& W, Matrix<N, M, T>& E, size_t i, size_t width, T v0, T beta, bool flip) {
    for (size_t k = i + 1; k < N; ++k) {
        T dot = W(k, i) * v0;
        for (size_t j = i + 1; j < N; ++j) {
            dot += W(k, j) * W(i, j);
        }
        for (size_t j = 0; j < width; ++j) {
            dot += E(k, j) * E(i, j);
        }
        if (dot == T{0}) {
            if (flip) {
                W(k, i) = -W(k, i);
            }
            continue;
        }
        const T s = dot * beta;
        const T wki = W(k, i) - (s * v0);
        W(k, i) = flip ? -wki : wki;
        for (size_t j = i + 1; j < N; ++j) {
            W(k, j) -= s * W(i, j);
        }
        for (size_t j = 0; j < width; ++j) {
            E(k, j) -= s * E(i, j);
        }
    }
}

/**
 * @brief [W E] Θ = [S 0] for an orthogonal Θ: W becomes the lower-triangular S with
 *        S Sᵀ = W Wᵀ + E Eᵀ, by Householder reflections, in place
 *
 * Row i's reflection acts on columns i… of W and on E, and maps the row onto its diagonal entry;
 * it is then applied to the rows below. The vector is the row itself with its first entry
 * replaced by @f$ x_0 + \operatorname{sign}(x_0)\,\mu @f$, @f$ \mu = \lVert\text{row}\rVert @f$,
 * so nothing cancels and nothing is copied, and the column's sign is chosen so the diagonal is
 * @f$ +\mu @f$.
 *
 * Householder and not modified Gram–Schmidt, for three reasons. It runs in the storage of W and E
 * with no second array, where Gram–Schmidt keeps the orthogonalized rows beside the factor. A row
 * that is zero (a state known exactly, a delay register) is no reflection; Gram–Schmidt divides
 * by its norm. And its error is a backward error of a few eps on each row of [W E] relative to
 * that row's own norm, so a state in small units keeps its digits beside one in large units.
 *
 * Squares of the entries are summed, as P's own entries would be, so the factor's range is the
 * square root of the scalar's: a row whose squares underflow (below 1e-19 in float) is taken as
 * zero.
 *
 * Cost: one square root per row and one division per row but the last, and about ⅔N³ + M N²
 * multiply-adds, or N³ for a lower-trapezoidal E with M = N.
 *
 * @tparam Trapezoid E(i, j) is zero for j > i and is not read there: the factor
 *                   process_noise_factor() returns for NW ≥ NX
 * @param W N × N, read in full; on return lower triangular with a non-negative diagonal
 * @param E N × M; overwritten with the reflections' workspace
 *
 * @see "Matrix Computations" (Golub & Van Loan, 2013), §5.1.2–5.1.3 and §5.2.2
 * @see "Accuracy and Stability of Numerical Algorithms" (Higham, 2002), §19.3 (error analysis of
 *      Householder computations)
 */
template<bool Trapezoid, size_t N, size_t M, typename T>
constexpr void triangularize(Matrix<N, N, T>& W, Matrix<N, M, T>& E) {
    for (size_t i = 0; i < N; ++i) {
        const size_t width = (Trapezoid && (i + 1 < M)) ? (i + 1) : M; // E's columns row i reaches
        T            tail{0};
        for (size_t j = i + 1; j < N; ++j) {
            tail += W(i, j) * W(i, j);
        }
        for (size_t j = 0; j < width; ++j) {
            tail += E(i, j) * E(i, j);
        }
        const T    x0 = W(i, i);
        const bool flip = !(x0 < T{0}); // the plain reflection leaves −sign(x₀)·μ on the diagonal
        if (!(tail >= damp::numeric_limits<T>::min())) {
            // Nothing to annihilate: the row is its diagonal entry already.
            if (!flip) {
                for (size_t k = i; k < N; ++k) {
                    W(k, i) = -W(k, i);
                }
            }
        } else {
            const T mu = damp::sqrt((x0 * x0) + tail);
            if (i + 1 < N) {
                const T v0 = flip ? (x0 + mu) : (x0 - mu);
                reflect_rows_below(W, E, i, width, v0, T{1} / (mu * damp::abs(v0)), flip); // β = 2 / vᵀv
            }
            W(i, i) = mu;
        }
        for (size_t j = i + 1; j < N; ++j) {
            W(i, j) = T{0};
        }
    }
}

/** @brief Columns of the process-noise factor of an NX-state model with NW noise inputs: NW when there are fewer inputs than states, else NX. */
template<size_t NX, size_t NW>
inline constexpr size_t noise_factor_columns = damp::min(NW, NX);

/** @brief Whether process_noise_factor() returns a lower-triangular factor (NW ≥ NX) or the dense G Q½ (NW < NX). */
template<size_t NX, size_t NW>
inline constexpr bool noise_factor_is_triangular = (NW >= NX);

/**
 * @brief A factor Γ of the process noise on the state, Γ Γᵀ = G Q Gᵀ, for sqrt_time_update()
 *
 * Q is factored by semidefinite_factor(), so a singular Q is taken as it is: a noise input with no
 * variance, or a model whose G has zero rows for its delay registers. With fewer noise inputs
 * than states Γ is G Q½, NX × NW and dense. Otherwise it is triangularized once to NX × NX lower
 * triangular, which is all the time update reads and costs it a third less.
 *
 * Computed once, where Q is set, never in a tick. A Q that is not positive semidefinite or not
 * finite is refused under constant evaluation; at run time Γ is zero and the filter adds no
 * process noise.
 *
 * @param G Process-noise input matrix (NX × NW)
 * @param Q Process-noise covariance on the inputs of G (NW × NW, positive semidefinite)
 * @return Γ: NX × noise_factor_columns<NX, NW>
 */
template<size_t NX, size_t NW, typename T>
[[nodiscard]] constexpr Matrix<NX, noise_factor_columns<NX, NW>, T> process_noise_factor(const Matrix<NX, NW, T>& G, const Matrix<NW, NW, T>& Q) {
    Matrix<NW, NW, T> Lq{};
    const bool        finite_G = mat::all_of(G, [](T v) { return damp::isfinite(v); });
    if (!finite_G || !semidefinite_factor(Q, Lq)) {
        if (damp::is_constant_evaluated()) {
            refuse_failed_design(design::DesignStatus::bad_weights); // constinit / constexpr: a Q that is not a covariance does not compile
        }
        return {};
    }
    Matrix<NX, NW, T> GL = G * Lq;
    if constexpr (noise_factor_is_triangular<NX, NW>) {
        Matrix<NX, NX, T> W{};
        triangularize<false>(W, GL);
        return W;
    } else {
        return GL;
    }
}

/** @brief The columns c ≥ @p J of row @p R of pattern @p P that are not zero: the terms of that row times column J of a lower-triangular matrix. */
template<const auto& P, size_t R, size_t J>
inline constexpr auto nonzero_columns_from = [] {
    using Pattern = damp::remove_cvref_t<decltype(P)>;
    struct Columns {
        size_t                               n{0};
        damp::array<size_t, Pattern::cols()> col{};
    };
    Columns out{};
    for (size_t c = J; c < Pattern::cols(); ++c) {
        if (P(R, c) != PatternEntry::zero) {
            out.col[out.n] = c;
            ++out.n;
        }
    }
    return out;
}();

/** @brief Row @p R of @p F times column @p J of the lower-triangular @p S: a term for each nonzero of the row at or below the diagonal of that column, and nothing else. */
template<size_t R, size_t J, const auto& P, size_t N, typename T>
[[nodiscard, gnu::always_inline]] constexpr T pattern_row_times_factor_column(const PatternMatrix<P, T>& F, const Matrix<N, N, T>& S) {
    return [&]<size_t... I>(damp::index_sequence<I...>) __attribute__((always_inline)) -> T {
        if constexpr (sizeof...(I) == 0) {
            return T{0};
        } else {
            return (... + F.template times<R, nonzero_columns_from<P, R, J>.col[I]>(S(nonzero_columns_from<P, R, J>.col[I], J)));
        }
    }(damp::make_index_sequence<nonzero_columns_from<P, R, J>.n>{});
}

/**
 * @brief @p S ← A S in place for a lower-triangular S: column j of the product reads column j of S
 *        only, so one column of scratch is all it takes
 *
 * The result is dense; triangularize() makes it a factor again.
 */
template<size_t N, typename T>
constexpr void times_factor_in_place(const Matrix<N, N, T>& A, Matrix<N, N, T>& S) {
    for (size_t j = 0; j < N; ++j) {
        T column[N]; // written before it is read
        for (size_t i = 0; i < N; ++i) {
            T s{0};
            for (size_t k = j; k < N; ++k) {
                s += A(i, k) * S(k, j);
            }
            column[i] = s;
        }
        for (size_t i = 0; i < N; ++i) {
            S(i, j) = column[i];
        }
    }
}

/** @brief @p S ← A S in place for an A whose zeros and ones are known at compile time: no term for a zero, no multiply for a one. */
template<const auto& P, size_t N, typename T>
[[gnu::always_inline]] constexpr void times_factor_in_place(const PatternMatrix<P, T>& A, Matrix<N, N, T>& S) {
    static_assert(P.rows() == N && P.cols() == N, "A S takes a square A of S's size.");
    [&]<size_t... J>(damp::index_sequence<J...>) __attribute__((always_inline)) {
        const auto column = [&]<size_t C>() __attribute__((always_inline)) {
            T col[N]; // written before it is read
            [&]<size_t... I>(damp::index_sequence<I...>) __attribute__((always_inline)) {
                ((col[I] = pattern_row_times_factor_column<I, C>(A, S)), ...);
                ((S(I, C) = col[I]), ...);
            }(damp::make_index_sequence<N>{});
        };
        (column.template operator()<J>(), ...);
    }(damp::make_index_sequence<N>{});
}

/**
 * @brief Time update of the factor: S ← the lower-triangular factor of A S Sᵀ Aᵀ + Γ Γᵀ
 *
 * The pre-array [A S, Γ] is triangularized by orthogonal reflections (triangularize()), so the
 * sum A P Aᵀ + G Q Gᵀ is never formed and cannot lose a small eigenvalue to a large one. A S is
 * written over S a column at a time, and the reflections run on a copy of Γ: the stack is NX
 * scalars and that copy, 0.7 kB at twelve states in float on a Cortex-M4F against the 1.8 kB of
 * the covariance form's products.
 *
 * Cost for NX states, Γ lower triangular (NW ≥ NX): ½NX³ multiply-adds for A S and NX³ for the
 * triangularization, with NX square roots and divisions.
 *
 * @tparam Triangular Whether Γ is lower triangular (noise_factor_is_triangular)
 * @param S Factor of the covariance, lower triangular; updated in place
 * @param A State transition matrix, dense or on a pattern (its zeros then cost nothing in A S)
 * @param Gamma Factor of the process noise on the state, from process_noise_factor()
 *
 * @see "Discrete square root filtering: a survey of current techniques" (Kaminski, Bryson &
 *      Schmidt, 1971), IEEE TAC 16(6)
 */
template<bool Triangular, size_t NX, size_t NQ, typename T, typename AMatrix>
constexpr void sqrt_time_update(Matrix<NX, NX, T>& S, const AMatrix& A, const Matrix<NX, NQ, T>& Gamma) {
    times_factor_in_place(A, S);
    Matrix<NX, NQ, T> E{detail::unfilled}; // the reflections overwrite it: a copy, of the part they read
    for (size_t i = 0; i < NX; ++i) {
        const size_t width = (Triangular && (i + 1 < NQ)) ? (i + 1) : NQ;
        for (size_t j = 0; j < width; ++j) {
            E(i, j) = Gamma(i, j);
        }
    }
    triangularize<Triangular>(S, E);
}

/** @brief @p F = C S for a lower-triangular S and a dense C. */
template<size_t NX, size_t NY, typename T>
constexpr void reading_times_factor(const Matrix<NY, NX, T>& C, const Matrix<NX, NX, T>& S, Matrix<NY, NX, T>& F) {
    for (size_t m = 0; m < NY; ++m) {
        for (size_t k = 0; k < NX; ++k) {
            T s{0};
            for (size_t i = k; i < NX; ++i) {
                s += C(m, i) * S(i, k);
            }
            F(m, k) = s;
        }
    }
}

/** @brief @p F = C S for a lower-triangular S and a C whose zeros and ones are known at compile time. */
template<const auto& P, size_t NX, size_t NY, typename T>
[[gnu::always_inline]] constexpr void reading_times_factor(const PatternMatrix<P, T>& C, const Matrix<NX, NX, T>& S, Matrix<NY, NX, T>& F) {
    static_assert(PatternMatrix<P, T>::NR == NY && PatternMatrix<P, T>::NC == NX, "C is NY × NX.");
    for_index_pairs<NY, NX>([&]<size_t M, size_t K>() __attribute__((always_inline)) { F(M, K) = pattern_row_times_factor_column<M, K>(C, S); });
}

/**
 * @brief A reading's noise covariance R factored for sqrt_measurement_update(): R = L Lᵀ, and
 *        whether R is a covariance at all
 *
 * L is semidefinite_factor()'s, so R may be singular: a row of the reading with no noise of its
 * own, a perfect measurement, has a zero on L's diagonal. An R that is not positive semidefinite
 * or not finite has no factor: @c valid is false, and an update with it is refused.
 */
template<size_t NY, typename T>
struct [[nodiscard]] ReadingNoiseFactor {
    Matrix<NY, NY, T> L{};          ///< Lower-triangular factor of R; a zero column for a row with no noise of its own
    bool              valid{false}; ///< R is positive semidefinite and finite

    /** @brief The same factor in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr ReadingNoiseFactor<NY, U> as() const {
        return {L.template as<U>(), valid};
    }
};

/** @brief The factor of a reading's noise covariance R (positive semidefinite; only its lower triangle is read), or an invalid one when R is indefinite or not finite. */
template<size_t NY, typename T>
[[nodiscard]] constexpr ReadingNoiseFactor<NY, T> reading_noise_factor(const Matrix<NY, NY, T>& R) {
    ReadingNoiseFactor<NY, T> out{};
    out.valid = semidefinite_factor(R, out.L);
    return out;
}

/**
 * @brief The size of the terms row @p m of C S is summed from, ‖|C_m|·|S|‖: what a perfect row's
 *        ‖f‖ is measured against, so that "zero" is the rounding of those terms and not a literal
 */
template<size_t NX, typename T, typename CMatrix>
[[nodiscard]] constexpr T reading_row_reach(const CMatrix& C, const Matrix<NX, NX, T>& S, size_t m) {
    T sum{0};
    for (size_t k = 0; k < NX; ++k) {
        T column{0};
        for (size_t i = k; i < NX; ++i) {
            column += damp::abs(C(m, i)) * damp::abs(S(i, k));
        }
        sum += column * column;
    }
    return damp::sqrt(sum);
}

/**
 * @brief The scale of row @p m of a reading that has no noise of its own: 1/√(eps·‖f‖²), its noise
 *        taken at the rounding of its own prior variance ‖f‖² = c̃ P c̃ᵀ
 *
 * A perfect measurement would set the variance of what it reads to zero. It is given the rounding
 * of that variance instead, eps times the prior, and is then a reading like any other, with
 * α = 1 + 1/eps: Carlson's update needs no special case. The estimate takes the reading to
 * within eps of the innovation, and the factor keeps √eps of the prior standard deviation there,
 * well above the rounding of its own entries, so a later reading still finds a variance to
 * divide by and not that rounding.
 *
 * @param terms The size of the terms f was summed from (reading_row_reach(), and the rows
 *              subtracted from it)
 * @return false if ‖f‖ is not finite, or is zero to the rounding of those terms, within
 *         64·(NX + NY)·eps of their size: no noise and no prior variance, C P Cᵀ + R singular
 */
template<size_t NX, size_t NY, typename T>
[[nodiscard]] constexpr bool perfect_row_scale(const Matrix<NY, NX, T>& F, size_t m, T terms, T& scale) {
    T norm2{0};
    for (size_t k = 0; k < NX; ++k) {
        norm2 += F(m, k) * F(m, k);
    }
    const T norm = damp::sqrt(norm2);
    const T rounding = T{64} * static_cast<T>(NX + NY) * eps<T>() * terms;
    if (!damp::isfinite(norm2) || !(norm > rounding)) {
        return false;
    }
    scale = T{1} / (damp::sqrt(eps<T>()) * norm);
    return true;
}

/**
 * @brief Whiten a reading in place with the factor L of its noise: each row of F and of ν less the
 *        rows before it through L, then scaled to unit noise
 *
 * The whitened reading has unit, uncorrelated noise, so its rows can be taken one at a time. A
 * row with a positive pivot is divided by it. A row with a zero pivot has no noise left once the
 * earlier rows are subtracted, a perfect measurement, and is scaled by perfect_row_scale().
 * Only then are the sizes of the rows' terms followed, for that row's test.
 *
 * @param F  C S on entry, C̃ S on return
 * @param nu The innovation on entry, whitened on return
 * @param L  Factor of the reading's noise covariance (ReadingNoiseFactor::L)
 * @param C  The reading's output matrix, read for a perfect row's test only
 * @param S  The factor F came from, read for a perfect row's test only
 * @return false if a perfect row has no prior variance either, or if a row's 1 + ‖f‖² is not
 *         finite: a NaN or ±Inf in C, or an innovation variance past the scalar's range. F and ν
 *         are scratch; nothing of the filter is written.
 */
template<size_t NX, size_t NY, typename T, typename CMatrix>
[[nodiscard]] constexpr bool whiten_reading(Matrix<NY, NX, T>& F, ColVec<NY, T>& nu, const Matrix<NY, NY, T>& L, const CMatrix& C, const Matrix<NX, NX, T>& S) {
    size_t perfect_end = 0; // one past the last row with no noise of its own
    for (size_t m = 0; m < NY; ++m) {
        if (!(L(m, m) > T{0})) {
            perfect_end = m + 1;
        }
    }
    T reach[NY] = {}; // the size of each whitened row's terms, kept up to the last perfect row
    for (size_t m = 0; m < NY; ++m) {
        const bool followed = m < perfect_end;
        T          terms = followed ? reading_row_reach(C, S, m) : T{0};
        T          v = nu[m];
        for (size_t j = 0; j < m; ++j) {
            const T l = L(m, j);
            v -= l * nu[j];
            for (size_t k = 0; k < NX; ++k) {
                F(m, k) -= l * F(j, k);
            }
            if (followed) {
                terms += damp::abs(l) * reach[j];
            }
        }
        T scale{0};
        if (L(m, m) > T{0}) {
            scale = T{1} / L(m, m);
        } else if (!perfect_row_scale(F, m, terms, scale)) {
            return false; // nothing to learn from it: no noise and no prior variance
        }
        nu[m] = v * scale;
        if (followed) {
            reach[m] = terms * scale;
        }
        T alpha{1};
        for (size_t k = 0; k < NX; ++k) {
            F(m, k) *= scale;
            alpha += F(m, k) * F(m, k);
        }
        if (!damp::isfinite(alpha)) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Carlson's triangular update for row @p m of a whitened reading: S ← S T with
 *        T Tᵀ = I − f fᵀ/α, f = Sᵀc̃ᵀ the row of @p F and α = 1 + fᵀf
 *
 * With @f$ \alpha_k = 1 + \sum_{j \ge k} f_j^2 @f$, taken from the last column to the first,
 * @f[
 *   s^+_{ik} = \sqrt{\tfrac{\alpha_{k+1}}{\alpha_k}}\, s_{ik}
 *              - \frac{f_k}{\sqrt{\alpha_{k+1}\alpha_k}} \sum_{j>k} s_{ij} f_j ,
 * @f]
 * and the sum, run to the first column, is row i of P c̃ᵀ: the gain is that over α₀. Each column
 * is one plane rotation of [√α, f; ·, S] with the gain column left unscaled, so the update is an
 * orthogonal triangularization at three multiplies an entry. S stays lower triangular, and its
 * diagonal is only ever scaled by a factor in (0, 1]: no entry of the factor is the difference of
 * two covariances.
 *
 * The same T is applied to the rows of F after m, which are the later readings' f on the updated
 * factor, and their innovations move with the estimate.
 *
 * Columns past the last nonzero of f are left alone: a reading of the leading states only costs
 * its own columns.
 */
template<size_t NX, size_t NY, typename T>
constexpr void carlson_update(Matrix<NX, NX, T>& S, ColVec<NX, T>& x, Matrix<NY, NX, T>& F, ColVec<NY, T>& nu, size_t m) {
    T      b[NX]; // √(α_{k+1}/α_k)
    T      c[NX]; // f_k/√(α_{k+1} α_k)
    T      alpha{1};
    T      root{1};
    T      inv_root{1};
    size_t top = 0; // one past the last column f moves
    for (size_t k = NX; k-- > 0;) {
        const T f = F(m, k);
        if (f == T{0}) {
            b[k] = T{1}; // an exact zero of f is the identity on its column: no square root
            c[k] = T{0};
            continue;
        }
        if (top == 0) {
            top = k + 1;
        }
        const T alpha_k = alpha + (f * f);
        const T root_k = damp::sqrt(alpha_k);
        const T inv_root_k = T{1} / root_k;
        b[k] = root * inv_root_k;
        c[k] = (f * inv_root_k) * inv_root;
        alpha = alpha_k;
        root = root_k;
        inv_root = inv_root_k;
    }
    const T gain = nu[m] * (inv_root * inv_root); // ν̃/α₀
    for (size_t i = 0; i < NX; ++i) {
        T e{0};
        for (size_t k = ((i + 1 < top) ? (i + 1) : top); k-- > 0;) {
            const T old = S(i, k);
            S(i, k) = (b[k] * old) - (c[k] * e);
            e += old * F(m, k);
        }
        x[i] += e * gain;
    }
    for (size_t r = m + 1; r < NY; ++r) {
        T e{0};
        for (size_t k = top; k-- > 0;) {
            const T old = F(r, k);
            F(r, k) = (b[k] * old) - (c[k] * e);
            e += old * F(m, k);
        }
        nu[r] -= e * gain;
    }
}

/**
 * @brief Measurement update of the factor and the estimate with a reading of innovation ν, output
 *        matrix C and noise covariance R = L Lᵀ
 *
 * The reading is whitened with R's factor, ν̃ = L⁻¹ν and C̃ = L⁻¹C, so a full R is handled, and
 * its rows are then taken one at a time by Carlson's update (carlson_update()). F = C̃ S is formed
 * once: each row's update is applied to the rows after it, so no row is formed twice.
 *
 * R may be singular. A row with no noise of its own, a perfect measurement (R = 0, or a
 * zero-variance row), is taken at the rounding of its own prior variance (perfect_row_scale()):
 * the estimate then reads C x̂ = y to within eps of the innovation, and what the row reads keeps
 * eps times its prior variance in the factor.
 *
 * Chosen against the two other factored updates, in float, on priors up to P₀/R = 10¹² and
 * (NX, NY) = (2, 1), (6, 3) and (12, 4). All three hold the covariance to the same accuracy: a
 * median error of 2·10⁻⁵ of the covariance itself at twelve states, and of 0.2 to 0.5 of it on
 * the first ticks of the widest prior. The array update, a triangularization of the whole
 * pre-array [[L, C S], [0, S]], takes 1.4 times the tick's run time and an (NX + NY)² array:
 * 1.1 kB of stack against 0.48 kB at (12, 4) on a Cortex-M4F. Bierman's U-D update takes no
 * square root and about the same stack; with its time update the tick runs 1.6 times faster at two
 * states and 1.35 times slower at twelve. It carries U and D, not a factor S, and that time
 * update orthogonalizes a copy of the whole NX × (NX + NW) pre-array, 1.4 kB against 0.7 kB at
 * twelve states.
 *
 * Every refusal is decided before anything is written, so a refused reading leaves @p S and
 * @p x bit for bit: an R that is not positive semidefinite or not finite (@p noise is not
 * valid), a NaN or ±Inf in C, an innovation variance past the range of T, or a perfect row whose
 * prior variance is zero as well, to the rounding of the terms it is summed from: a perfect
 * reading of a state known exactly. That last one is the singular C P Cᵀ + R the covariance
 * form refused; with any noise, or any prior variance, C P Cᵀ + R is positive definite by
 * construction. Perfect rows that read the same thing twice are both taken, the second against
 * the rounding the first left.
 *
 * Cost for NX states and NY readings: about 2·NY·NX² multiply-adds, NY·NX square roots and
 * divisions, and NY × NX + 2·NX + NY scalars of stack: 0.5 kB for the whole update at (12, 4) in
 * float on a Cortex-M4F, against the 4.6 kB of the Joseph form. A reading with a perfect row
 * costs ½·NX² more multiply-adds for each row up to it.
 *
 * @param S     Factor of the covariance, lower triangular; updated in place
 * @param x     Estimate; corrected in place
 * @param nu    Innovation y − ŷ (NY rows)
 * @param C     Output matrix of the reading, dense or on a pattern
 * @param noise The reading's noise covariance, factored (reading_noise_factor())
 * @return false if the reading cannot be used: S and x are unchanged
 *
 * @see "Fast triangular formulation of the square root filter" (Carlson, 1973), AIAA Journal 11(9)
 * @see "Factorization Methods for Discrete Sequential Estimation" (Bierman, 1977), ch. V
 */
template<size_t NX, size_t NY, typename T, typename CMatrix>
[[nodiscard]] constexpr bool sqrt_measurement_update(Matrix<NX, NX, T>& S, ColVec<NX, T>& x, const ColVec<NY, T>& nu, const CMatrix& C, const ReadingNoiseFactor<NY, T>& noise) {
    if (!noise.valid) {
        return false; // R is not a covariance
    }
    Matrix<NY, NX, T> F{detail::unfilled};
    ColVec<NY, T>     whitened = nu;
    reading_times_factor(C, S, F);
    if (!whiten_reading(F, whitened, noise.L, C, S)) {
        return false;
    }
    for (size_t m = 0; m < NY; ++m) {
        carlson_update(S, x, F, whitened, m);
    }
    return true;
}

} // namespace damp::detail
