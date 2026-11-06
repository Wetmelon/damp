// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file zeros.hpp
 * @brief Invariant zeros of a state-space model, any number of inputs and outputs:
 *        transmission_zeros (MATLAB®'s tzero)
 *
 * The zeros are the finite λ at which the system matrix
 * @f[
 *   P(\lambda) = \begin{bmatrix} A - \lambda I & B \\ C & D \end{bmatrix}
 * @f]
 * loses rank below its normal rank. The pencil is reduced by orthogonal transformations to a
 * square regular pencil with the same finite zeros (Emami-Naeini & Van Dooren), whose
 * generalized eigenvalues the QZ algorithm then computes (mat::generalized_eigenvalues). No
 * transfer-function coefficients are formed and D is never inverted.
 *
 * @code
 * // G(s) = [1/(s+1) 2/(s+2); 1/(s+3) 1/(s+4)]: det G(s)·∏(s+k) = −(s² + 5s + 2)
 * constexpr StateSpace<4, 2, 2> sys{
 *     .A = Matrix<4, 4>{{-1.0, 0.0, 0.0, 0.0}, {0.0, -2.0, 0.0, 0.0}, {0.0, 0.0, -3.0, 0.0}, {0.0, 0.0, 0.0, -4.0}},
 *     .B = Matrix<4, 2>{{1.0, 0.0}, {0.0, 1.0}, {1.0, 0.0}, {0.0, 1.0}},
 *     .C = Matrix<2, 4>{{1.0, 2.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 1.0}},
 * };
 * constexpr auto z = transmission_zeros(sys);
 * static_assert(z.success && z.n_zeros == 2); // z.zeros[0 .. 2) = (−5 ± √17)/2
 * @endcode
 */

#include <cstddef>
#include <limits>

#include "damp/backend.hpp"
#include "damp/design/status.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/core.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/generalized_eigen.hpp"
#include "state_space.hpp"

namespace damp {

/**
 * @brief The finite invariant zeros of a state-space model and their count
 *
 * @c zeros[0..n_zeros) hold them; a complex pair occupies two consecutive entries. At most NX
 * zeros are finite.
 *
 * @tparam NX States of the model
 * @tparam T  Floating-point scalar type
 */
template<size_t NX, typename T = double>
struct [[nodiscard]] TransmissionZerosResult {
    damp::array<damp::complex<T>, NX> zeros{};                                     ///< Finite zeros, the first n_zeros entries
    size_t                            n_zeros{0};                                  ///< Number of finite zeros
    bool                              success{false};                              ///< status == DesignStatus::ok
    design::DesignStatus              status{design::DesignStatus::not_converged}; ///< ok, bad_weights (a non-finite entry) or not_converged (QZ)

    /** @brief True when the zeros were found (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The same zeros in the scalar @p U. */
    template<typename U>
    [[nodiscard]] constexpr TransmissionZerosResult<NX, U> as() const {
        TransmissionZerosResult<NX, U> out{};
        for (size_t i = 0; i < NX; ++i) {
            out.zeros[i] = damp::complex<U>{static_cast<U>(zeros[i].real()), static_cast<U>(zeros[i].imag())};
        }
        out.n_zeros = n_zeros;
        out.success = success;
        out.status = status;
        return out;
    }
};

namespace detail {

/** @brief e with |x| = 2^e·f, f ∈ [1, 2); 0 for x = 0 or non-finite x. */
template<typename T>
[[nodiscard]] constexpr int binary_exponent(T x) {
    T v = damp::abs(x);
    if (!(v > T{0}) || !damp::isfinite(v)) {
        return 0;
    }
    constexpr int cap = std::numeric_limits<T>::max_exponent - std::numeric_limits<T>::min_exponent + std::numeric_limits<T>::digits;
    int           e = 0;
    while (v >= T{2} && e < cap) {
        v *= static_cast<T>(0.5);
        ++e;
    }
    while (v < T{1} && e > -cap) {
        v *= T{2};
        --e;
    }
    return e;
}

/**
 * @brief The system matrix [A B; C D] of a model with n states, m inputs and p outputs, held in
 *        the leading (n + p) × (n + m) block of a K × K buffer (zero outside it)
 */
template<size_t K, typename T>
struct SystemPencil {
    Matrix<K, K, T> S{};              ///< [A B; C D]
    size_t          n{0};             ///< States
    size_t          m{0};             ///< Inputs
    size_t          p{0};             ///< Outputs
    int             time_exponent{0}; ///< The pencil's zeros times 2^time_exponent are the model's
};

/** @brief A Householder reflector H = I − β·u·uᵀ acting on the indices [first, last). */
template<size_t K, typename T>
struct IndexReflector {
    damp::array<T, K> u{};      ///< Householder vector on [first, last), scaled to its largest entry
    T                 beta{0};  ///< 2/(uᵀu), or 0 for the identity
    size_t            first{0}; ///< First index acted on
    size_t            last{0};  ///< One past the last index acted on
};

/** @brief The reflector on [first, last) with H·x = ∓‖x‖·e_target, x(i) the entry at index i. */
template<size_t K, typename T, typename Entry>
[[nodiscard]] constexpr IndexReflector<K, T> reflector_onto(size_t first, size_t last, size_t target, Entry x) {
    IndexReflector<K, T> h{};
    h.first = first;
    h.last = last;
    T scale = T{0};
    for (size_t i = first; i < last; ++i) {
        scale = damp::max(scale, damp::abs(x(i)));
    }
    if (!(scale > T{0})) {
        return h;
    }
    T nrm_sq = T{0};
    for (size_t i = first; i < last; ++i) {
        h.u[i] = x(i) / scale;
        nrm_sq += h.u[i] * h.u[i];
    }
    const T nrm = damp::sqrt(nrm_sq);
    const T ut = h.u[target];
    h.u[target] = ut + damp::copysign(nrm, ut);
    h.beta = T{1} / (nrm * (nrm + damp::abs(ut)));
    return h;
}

/** @brief Rows [h.first, h.last) of M on the columns [0, c1) ← H·(those rows). */
template<size_t K, typename T>
constexpr void reflect_pencil_rows(Matrix<K, K, T>& M, const IndexReflector<K, T>& h, size_t c1) {
    for (size_t j = 0; j < c1; ++j) {
        T w = T{0};
        for (size_t i = h.first; i < h.last; ++i) {
            w += h.u[i] * M(i, j);
        }
        w *= h.beta;
        for (size_t i = h.first; i < h.last; ++i) {
            M(i, j) -= w * h.u[i];
        }
    }
}

/** @brief Columns [h.first, h.last) of M on the rows [0, r1) ← (those columns)·H. */
template<size_t K, typename T>
constexpr void reflect_pencil_cols(Matrix<K, K, T>& M, const IndexReflector<K, T>& h, size_t r1) {
    for (size_t i = 0; i < r1; ++i) {
        T w = T{0};
        for (size_t j = h.first; j < h.last; ++j) {
            w += M(i, j) * h.u[j];
        }
        w *= h.beta;
        for (size_t j = h.first; j < h.last; ++j) {
            M(i, j) -= w * h.u[j];
        }
    }
}

/**
 * @brief The system matrix of @p sys scaled by powers of two (exact): states balanced as A is,
 *        time by the size of A, inputs and outputs to unit columns of [B; D] and rows of [C D]
 *
 * With Dₓ the balancing diagonal of A, 2^a nearest ‖Dₓ⁻¹ADₓ‖_F from below, and Dᵤ, D_y the input
 * and output scalings, the pencil is [Dₓ⁻¹ADₓ/2^a, Dₓ⁻¹BDᵤ/2^a; D_yCDₓ, D_yDDᵤ]. Its zeros are the
 * model's divided by 2^a, so the rank decisions on it do not depend on the units of time, of
 * the inputs or of the outputs, and depend on the state units only as far as balancing A
 * removes them.
 */
template<size_t K, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr SystemPencil<K, T> scaled_system_pencil(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    SystemPencil<K, T> out{.n = NX, .m = NU, .p = NY};
    const auto         e = mat::detail::balance_exponents(sys.A);
    const auto         A = mat::detail::balanced(sys.A, e);
    const int          a = binary_exponent(A.norm());
    out.time_exponent = a;
    auto& S = out.S;
    for (size_t i = 0; i < NX; ++i) {
        for (size_t j = 0; j < NX; ++j) {
            S(i, j) = mat::detail::scale_by_pow2(A(i, j), -a);
        }
        for (size_t j = 0; j < NU; ++j) {
            S(i, NX + j) = mat::detail::scale_by_pow2(sys.B(i, j), -e[i] - a);
        }
    }
    for (size_t k = 0; k < NY; ++k) {
        for (size_t j = 0; j < NX; ++j) {
            S(NX + k, j) = mat::detail::scale_by_pow2(sys.C(k, j), e[j]);
        }
        for (size_t j = 0; j < NU; ++j) {
            S(NX + k, NX + j) = sys.D(k, j);
        }
    }
    for (size_t j = NX; j < NX + NU; ++j) {
        const int u = binary_exponent(scaled_norm2<T>(NX + NY, [&](size_t i) { return damp::abs(S(i, j)); }));
        for (size_t i = 0; i < NX + NY; ++i) {
            S(i, j) = mat::detail::scale_by_pow2(S(i, j), -u);
        }
    }
    for (size_t i = NX; i < NX + NY; ++i) {
        const int y = binary_exponent(scaled_norm2<T>(NX + NU, [&](size_t j) { return damp::abs(S(i, j)); }));
        for (size_t j = 0; j < NX + NU; ++j) {
            S(i, j) = mat::detail::scale_by_pow2(S(i, j), -y);
        }
    }
    return out;
}

/**
 * @brief Row-compress D by column-pivoted Householder reflections on the outputs; returns its rank
 *
 * U·D = [D̄; 0] with D̄ of full row rank τ, U applied to all of [C D]. A pivot column whose norm
 * is at most @p tol ends the compression, and what remains below the first τ rows is set to zero.
 */
template<size_t K, typename T>
[[nodiscard]] constexpr size_t compress_feedthrough_rows(SystemPencil<K, T>& P, T tol) {
    auto&                S = P.S;
    const size_t         n = P.n;
    const size_t         rows = P.n + P.p;
    damp::array<bool, K> used{};
    size_t               tau = 0;
    while (tau < P.p && tau < P.m) {
        const size_t r0 = n + tau;
        size_t       best = P.m;
        T            best_norm{0};
        for (size_t c = 0; c < P.m; ++c) {
            const T nc = used[c] ? T{0} : scaled_norm2<T>(rows - r0, [&](size_t k) { return damp::abs(S(r0 + k, n + c)); });
            if (nc > best_norm) {
                best_norm = nc;
                best = c;
            }
        }
        if (best == P.m || !(best_norm > tol)) {
            break;
        }
        used[best] = true;
        const auto h = reflector_onto<K, T>(r0, rows, r0, [&](size_t i) { return S(i, n + best); });
        reflect_pencil_rows(S, h, n + P.m);
        for (size_t i = r0 + 1; i < rows; ++i) {
            S(i, n + best) = T{0};
        }
        ++tau;
    }
    for (size_t i = n + tau; i < rows; ++i) {
        for (size_t j = n; j < n + P.m; ++j) {
            S(i, j) = T{0}; // below the rank tolerance
        }
    }
    return tau;
}

/**
 * @brief Compress the rows [r0, r1) of S on the columns [0, e) into its trailing columns, one row
 *        at a time, largest remaining row first, while that row's norm exceeds @p tol; returns
 *        how many columns it took
 *
 * Each step reflects the columns [0, e) so the chosen row ends with a single entry at column
 * e − 1, then shortens e. @p also_apply(H) then applies the same reflector elsewhere: to the
 * rows [0, e) for a change of state coordinates (which keeps the λI block), or to the columns
 * of the λ-part.
 */
template<size_t K, typename T, typename AlsoApply>
constexpr size_t compress_into_trailing_columns(
    SystemPencil<K, T>& P,
    size_t              r0,
    size_t              r1,
    size_t              e,
    size_t              steps,
    T                   tol,
    AlsoApply           also_apply
) {
    auto&                S = P.S;
    damp::array<bool, K> done{};
    size_t               taken = 0;
    while (taken < steps && taken < e) {
        const size_t width = e - taken;
        size_t       best = r1;
        T            best_norm{0};
        for (size_t r = r0; r < r1; ++r) {
            const T nr = done[r] ? T{0} : scaled_norm2<T>(width, [&](size_t j) { return damp::abs(S(r, j)); });
            if (nr > best_norm) {
                best_norm = nr;
                best = r;
            }
        }
        if (best == r1 || !(best_norm > tol)) {
            break;
        }
        done[best] = true;
        const auto h = reflector_onto<K, T>(0, width, width - 1, [&](size_t j) { return S(best, j); });
        reflect_pencil_cols(S, h, P.n + P.p);
        also_apply(h);
        for (size_t j = 0; j + 1 < width; ++j) {
            S(best, j) = T{0};
        }
        ++taken;
    }
    return taken;
}

/**
 * @brief Reduce the system to one with the same finite zeros whose D has full row rank
 *
 * Each pass row-compresses D to [D̄; 0] (rank τ). If rows are left over, their C part C̃ is
 * column-compressed by a change of state coordinates to [0 Ĉ], Ĉ of full column rank ρ. Those
 * ρ states are then fixed by the left-over outputs, and the pencil loses them: their state
 * equations become outputs, so the new model has n − ρ states and ρ + τ outputs, and the
 * left-over rows are dropped. In the stored system matrix that is deleting the ρ state columns
 * and the left-over rows. With ρ = 0 the left-over rows are zero and are dropped alone. Each
 * pass with ρ > 0 removes a state, so at most n + 1 passes run.
 *
 * @see "Computation of zeros of linear multivariable systems" (Emami-Naeini & Van Dooren, 1982),
 *      Automatica 18(4), §3 (the procedure reduce)
 */
template<size_t K, typename T>
constexpr void reduce_to_full_row_rank_feedthrough(SystemPencil<K, T>& P, T tol) {
    auto& S = P.S;
    for (size_t pass = 0; pass <= K; ++pass) {
        const size_t n = P.n;
        const size_t rows = P.n + P.p;
        const size_t cols = P.n + P.m;
        const size_t tau = compress_feedthrough_rows(P, tol);
        if (tau == P.p) {
            return;
        }
        const size_t rho = compress_into_trailing_columns(P, n + tau, rows, n, P.p - tau, tol, [&](const IndexReflector<K, T>& h) {
            reflect_pencil_rows(S, h, cols); // H·S·H: a change of state coordinates
        });
        for (size_t i = n + tau; i < rows; ++i) {
            for (size_t j = 0; j < n - rho; ++j) {
                S(i, j) = T{0}; // below the rank tolerance
            }
        }
        for (size_t i = 0; i < n + tau; ++i) {
            for (size_t j = n; j < cols; ++j) {
                S(i, j - rho) = S(i, j);
            }
        }
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = (i < n + tau) ? cols - rho : 0; j < cols; ++j) {
                S(i, j) = T{0}; // outside the new (n − ρ + ρ + τ) × (n − ρ + m) block
            }
        }
        P.n = n - rho;
        P.p = rho + tau;
        if (rho == 0) {
            return;
        }
    }
}

/** @brief The dual system (Aᵀ, Cᵀ, Bᵀ, Dᵀ): the transposed system matrix, inputs and outputs swapped. */
template<size_t K, typename T>
constexpr void transpose_system(SystemPencil<K, T>& P) {
    const size_t extent = P.n + damp::max(P.m, P.p);
    for (size_t i = 0; i < extent; ++i) {
        for (size_t j = 0; j < i; ++j) {
            const T t = P.S(i, j);
            P.S(i, j) = P.S(j, i);
            P.S(j, i) = t;
        }
    }
    const size_t m = P.m;
    P.m = P.p;
    P.p = m;
}

} // namespace detail

/**
 * @brief The finite invariant zeros of a state-space model with any number of inputs and outputs
 *
 * The zeros of the system matrix P(λ) = [A − λI, B; C, D]: the finite λ at which it loses rank
 * below its normal rank. For a minimal realization these are the transmission zeros. A
 * non-minimal one adds its decoupling zeros: an unobservable mode of a square or tall model, an
 * uncontrollable mode of a square or wide one. Continuous or discrete alike (the zeros are in s
 * or z).
 *
 * The model is first scaled by powers of two (exact): states balanced as A is, time by the size
 * of A, inputs and outputs to unit size. Then the reduction of Emami-Naeini and Van Dooren: the
 * procedure reduce (detail::reduce_to_full_row_rank_feedthrough) leaves a model with the same
 * finite zeros whose D has full row rank, the same on the dual model leaves a square D of full
 * rank, and an orthogonal compression of the columns of [C D] onto D leaves the square regular
 * pencil A_f − λE_f, n_f ≤ NX, whose generalized eigenvalues (QZ) are the zeros. Each rank
 * decision takes as zero what is at most (NX + NY)(NX + NU)·ε·‖S‖_F, S the scaled system
 * matrix, so a D at roundoff level is treated as zero and its zeros as infinite, and a small D
 * above it gives the large finite zeros it implies. Only orthogonal transformations are used:
 * the zeros are those of a system within that tolerance of the scaled model.
 *
 * A zero of a non-square model exists only for special data (a zero common to every output, a
 * decoupling zero); a generic non-square model has none.
 *
 * Statuses: DesignStatus::bad_weights for a non-finite entry; DesignStatus::not_converged when
 * the QZ iteration does not converge.
 *
 * @note Compare with MATLAB®'s tzero(sys); matlab::tzero is the alias.
 * @see "Computation of zeros of linear multivariable systems" (Emami-Naeini & Van Dooren, 1982),
 *      Automatica 18(4)
 * @see "The computation of Kronecker's canonical form of a singular pencil" (Van Dooren, 1979),
 *      Linear Algebra Appl. 27
 * @see "An Algorithm for Generalized Matrix Eigenvalue Problems" (Moler & Stewart, 1973), SIAM J.
 *      Numer. Anal. 10(2)
 *
 * @param sys Model, continuous or discrete, any NU and NY
 * @return TransmissionZerosResult: zeros[0 .. n_zeros), success and status
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr TransmissionZerosResult<NX, T> transmission_zeros(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    using Result = TransmissionZerosResult<NX, T>;
    const auto finite_entry = [](T v) { return damp::isfinite(v); };
    if (!mat::all_of(sys.A, finite_entry) || !mat::all_of(sys.B, finite_entry) || !mat::all_of(sys.C, finite_entry)
        || !mat::all_of(sys.D, finite_entry)) {
        return detail::refused<Result>(design::DesignStatus::bad_weights);
    }
    Result out{};
    if constexpr (NX == 0) {
        return detail::set_status(out, design::DesignStatus::ok);
    } else {
        constexpr size_t K = NX + ((NU > NY) ? NU : NY);
        auto             P = detail::scaled_system_pencil<K>(sys);
        const T          tol = static_cast<T>((NX + NY) * (NX + NU)) * eps<T>() * P.S.norm();
        for (size_t round = 0; round < (2 * K) + 2; ++round) {
            detail::reduce_to_full_row_rank_feedthrough(P, tol);
            if (P.p == P.m) {
                break;
            }
            detail::transpose_system(P); // D is wide: reduce the dual
        }
        if (P.p != P.m) {
            return detail::refused<Result>(design::DesignStatus::not_converged);
        }
        //! Compress [C D] onto its last p columns: [C D]·V = [0 D̃], D̃ square and of full rank.
        Matrix<K, K, T> E{};
        for (size_t i = 0; i < P.n; ++i) {
            E(i, i) = T{1};
        }
        detail::compress_into_trailing_columns(P, P.n, P.n + P.p, P.n + P.m, P.p, T{0}, [&](const detail::IndexReflector<K, T>& h) {
            detail::reflect_pencil_cols(E, h, P.n);
        });
        Matrix<NX, NX, T> Af{};
        Matrix<NX, NX, T> Ef{};
        for (size_t i = 0; i < P.n; ++i) {
            for (size_t j = 0; j < P.n; ++j) {
                Af(i, j) = P.S(i, j);
                Ef(i, j) = E(i, j);
            }
        }
        damp::array<damp::complex<T>, NX> lambda{};
        damp::array<bool, NX>             finite{};
        if (!mat::detail::pencil_eigenvalues(Af, Ef, P.n, lambda, finite)) {
            return detail::refused<Result>(design::DesignStatus::not_converged);
        }
        for (size_t i = 0; i < P.n; ++i) {
            if (finite[i]) {
                out.zeros[out.n_zeros] = damp::complex<T>{
                    mat::detail::scale_by_pow2(lambda[i].real(), P.time_exponent),
                    mat::detail::scale_by_pow2(lambda[i].imag(), P.time_exponent),
                };
                ++out.n_zeros;
            }
        }
        return detail::set_status(out, design::DesignStatus::ok);
    }
}

} // namespace damp
