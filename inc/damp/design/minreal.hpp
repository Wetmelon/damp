// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file minreal.hpp
 * @brief Minimal realization (cancel uncontrollable / unobservable modes)
 *
 * Removes dynamics that do not appear in the input–output map by projecting
 * onto the controllable-and-observable (Kalman) subspace. Continuous and
 * discrete SISO/MIMO @c StateSpace models are supported. This is structural
 * pole–zero cancellation / Kalman reduction — not balanced truncation
 * (balreal / balred / modred).
 *
 * Algorithm (fixed-size, stack-only):
 * 1. Controllable subspace of @f$ (A, B) @f$ by the block-Arnoldi staircase
 *    (stability::detail::controllable_subspace): each stair adds the new directions
 *    of @f$ A @f$ applied to the previous one, tested against the relative tolerance.
 * 2. Orthogonal change of basis that places that subspace first; form the
 *    controllable subsystem @f$ (A_c, B_c, C_c) @f$.
 * 3. Observable subspace of @f$ (A_c, C_c) @f$: the same staircase on the dual
 *    @f$ (A_c^\top, C_c^\top) @f$.
 * 4. Project @f$ (A,B,C,D) @f$ onto that CO subspace. The I/O map is preserved;
 *    eliminated modes do not contribute to @f$ C(sI-A)^{-1}B + D @f$.
 *
 * The rank of the Krylov matrices @f$ [B, AB, \dots] @f$ would answer 1 and 3 in exact
 * arithmetic, but their columns align with the dominant mode and the matrices lose
 * rank numerically from about ten states, or sooner for modes spread over decades,
 * which dropped genuine modes. The staircase tests each stair on its own.
 *
 * The reduced model is returned both as a typed extract (when the order is
 * known at compile time) and embedded in the original state dimension (leading
 * principal block = minimal realization; trailing states zero).
 *
 * Example: drop an uncontrollable mode
 * @code
 * #include "damp/design/minreal.hpp"
 * using namespace damp;
 *
 * // A = diag(-1, -2), only the first state is driven; both appear in C
 * constexpr StateSpace<2, 1, 1> sys{
 *     .A = Matrix<2, 2>{{-1.0, 0.0}, {0.0, -2.0}},
 *     .B = Matrix<2, 1>{{1.0}, {0.0}},
 *     .C = Matrix<1, 2>{{1.0, 1.0}},
 * };
 * constexpr auto mr = design::minreal(sys);
 * static_assert(mr.success);
 * static_assert(mr.order == 1);           // second mode cancelled
 * constexpr auto red = mr.extract<1>();   // StateSpace<1,1,1>
 * static_assert(red.has_value());
 * @endcode
 *
 * @note Compare with MATLAB®'s minreal(sys).
 * @see balreal(), balred(), modred() for Hankel-based order reduction
 * @see stability::controllability_matrix(), stability::observability_matrix()
 * @see Kailath, "Linear Systems" (1980), §6.4 (Kalman decomposition)
 * @see Chen, "Linear System Theory and Design" (3rd ed.), §6.4–6.5
 */

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/design/stability.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "damp/matrix/svd.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace design {

/**
 * @brief Minimal realization result
 *
 * @c system embeds the reduced dynamics in the original state dimension: the
 * leading @c order × @c order block of @c A (and matching @c B / @c C
 * rows/columns) is a minimal realization; trailing dynamics are zero.
 * Use @ref extract to obtain a typed @c StateSpace of size @c order.
 *
 * @note Compare with MATLAB®'s minreal(sys) output (reduced @c ss).
 * @see minreal(), minimal_realization()
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
struct MinimalRealizationResult {
    StateSpace<NX, NU, NY, T, NW, NV> system{};       ///< Reduced model (leading @c order states)
    size_t                            order{0};       ///< Dimension of controllable ∩ observable subspace
    size_t                            eliminated{0};  ///< @c NX − @c order
    bool                              success{false}; ///< true if SVDs converged and projection succeeded

    template<typename U>
    [[nodiscard]] constexpr MinimalRealizationResult<NX, NU, NY, U, NW, NV> as() const {
        return MinimalRealizationResult<NX, NU, NY, U, NW, NV>{
            system.template as<U>(),
            order,
            eliminated,
            success
        };
    }

    /**
     * @brief Typed reduced model of compile-time order @p NR
     *
     * @tparam NR Expected reduced order (must match @c order at runtime)
     * @return Reduced @c StateSpace, or @c damp::nullopt if @c !success or @c order != NR
     */
    template<size_t NR>
    [[nodiscard]] constexpr damp::optional<StateSpace<NR, NU, NY, T, NW, NV>> extract() const {
        static_assert(NR > 0 && NR <= NX, "minreal extract: require 0 < NR <= NX");
        if (!success || order != NR) {
            return damp::nullopt;
        }
        if constexpr (NR == NX) {
            return system;
        } else {
            StateSpace<NR, NU, NY, T, NW, NV> out{};
            for (size_t i = 0; i < NR; ++i) {
                for (size_t j = 0; j < NR; ++j) {
                    out.A(i, j) = system.A(i, j);
                }
                for (size_t j = 0; j < NU; ++j) {
                    out.B(i, j) = system.B(i, j);
                }
            }
            for (size_t i = 0; i < NY; ++i) {
                for (size_t j = 0; j < NR; ++j) {
                    out.C(i, j) = system.C(i, j);
                }
                for (size_t j = 0; j < NU; ++j) {
                    out.D(i, j) = system.D(i, j);
                }
            }
            if constexpr (NW > 0) {
                for (size_t i = 0; i < NR; ++i) {
                    for (size_t j = 0; j < NW; ++j) {
                        out.G(i, j) = system.G(i, j);
                    }
                }
            }
            out.H = system.H;
            out.Ts = system.Ts;
            return out;
        }
    }
};

namespace detail {

/// Y = M(:, 0:ncols)' * X  (left-multiply by basis columns transpose).
template<size_t N, size_t C, typename T>
[[nodiscard]] constexpr Matrix<N, C, T> mul_basis_cols_t(
    const Matrix<N, N, T>& basis,
    size_t                 ncols,
    const Matrix<N, C, T>& X
) {
    Matrix<N, C, T> Y{};
    const size_t    nc = (ncols < N) ? ncols : N;
    for (size_t j = 0; j < C; ++j) {
        for (size_t i = 0; i < nc; ++i) {
            T sum = T{0};
            for (size_t k = 0; k < N; ++k) {
                sum += basis(k, i) * X(k, j);
            }
            Y(i, j) = sum;
        }
    }
    return Y;
}

/// S = T(:,0:n)' * A * T(:,0:n) stored in the leading n×n of an N×N matrix.
template<size_t N, typename T>
[[nodiscard]] constexpr Matrix<N, N, T> restrict_A(
    const Matrix<N, N, T>& A,
    const Matrix<N, N, T>& Tbasis,
    size_t                 n
) {
    Matrix<N, N, T> out{};
    if (n == 0) {
        return out;
    }
    // W = A * T(:,0:n)
    Matrix<N, N, T> W{};
    for (size_t j = 0; j < n; ++j) {
        for (size_t i = 0; i < N; ++i) {
            T sum = T{0};
            for (size_t k = 0; k < N; ++k) {
                sum += A(i, k) * Tbasis(k, j);
            }
            W(i, j) = sum;
        }
    }
    // out(i,j) = Tcol_i · Wcol_j
    for (size_t j = 0; j < n; ++j) {
        for (size_t i = 0; i < n; ++i) {
            T sum = T{0};
            for (size_t k = 0; k < N; ++k) {
                sum += Tbasis(k, i) * W(k, j);
            }
            out(i, j) = sum;
        }
    }
    return out;
}

/// C_r = C * T(:,0:n) stored in the leading columns of NY×N.
template<size_t NY, size_t N, typename T>
[[nodiscard]] constexpr Matrix<NY, N, T> restrict_C(
    const Matrix<NY, N, T>& C,
    const Matrix<N, N, T>&  Tbasis,
    size_t                  n
) {
    Matrix<NY, N, T> out{};
    for (size_t j = 0; j < n; ++j) {
        for (size_t i = 0; i < NY; ++i) {
            T sum = T{0};
            for (size_t k = 0; k < N; ++k) {
                sum += C(i, k) * Tbasis(k, j);
            }
            out(i, j) = sum;
        }
    }
    return out;
}

} // namespace detail

/**
 * @brief Minimal realization — cancel uncontrollable and unobservable modes
 *
 * Projects @p sys onto the controllable-and-observable subspace via SVD of the
 * controllability and observability matrices (orthogonal range test). Continuous
 * and discrete systems are treated the same (structural property of (A,B,C)).
 *
 * MIMO and SISO are both supported. Numerical rank uses a relative singular-value
 * threshold @p tol · σ_max (same convention as @c mat::rank_from_svd). Modes with
 * singular values at or below that floor are treated as cancelled.
 *
 * @note This is not balanced truncation. Use balred / modred to drop
 *       weakly coupled but still controllable and observable modes.
 *
 * @note Compare with MATLAB®'s minreal(sys, tol).
 *
 * @see minimal_realization(), balreal(), balred(), modred()
 * @see stability::is_controllable(), stability::is_observable()
 * @see Kailath, "Linear Systems" (1980), §6.4
 *
 * @tparam NX Number of states
 * @tparam NU Number of inputs
 * @tparam NY Number of outputs
 * @tparam T  Scalar type (default double)
 * @tparam NW Process-noise columns (passed through the same basis)
 * @tparam NV Measurement-noise columns (unchanged)
 * @param sys State-space model (continuous if Ts == 0, else discrete)
 * @param tol Relative tolerance on each stair's singular values (default √eps, as MATLAB®'s minreal)
 * @return MinimalRealizationResult with reduced model, order, success flag
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr MinimalRealizationResult<NX, NU, NY, T, NW, NV>
minreal(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T tol = damp::sqrt(eps<T>())) {
    if constexpr (NX == 0) {
        return MinimalRealizationResult<NX, NU, NY, T, NW, NV>{
            .system = sys,
            .order = 0,
            .eliminated = 0,
            .success = true,
        };
    }

    // --- Controllable subspace (staircase) -------------------------------------
    const auto cs = stability::detail::controllable_subspace(sys.A, sys.B, tol);
    if (!cs.converged) {
        return MinimalRealizationResult<NX, NU, NY, T, NW, NV>{};
    }
    const size_t             nc = cs.dim;
    const Matrix<NX, NX, T>& Uc = cs.basis; // leading nc columns: orthonormal basis

    if (nc == 0) {
        // No controllable dynamics: I/O is pure feedthrough.
        return MinimalRealizationResult<NX, NU, NY, T, NW, NV>{
            .system = StateSpace<NX, NU, NY, T, NW, NV>{
                .D = sys.D,
                .H = sys.H,
                .Ts = sys.Ts,
            },
            .order = 0,
            .eliminated = NX,
            .success = true,
        };
    }

    // --- Observable subspace of the controllable subsystem ------------------
    // (Ac, Cc) = (Ucᵀ A Uc, C Uc) in the leading nc coordinates; the dual staircase on
    // (Acᵀ, Ccᵀ) stays inside them because the padding is zero.
    const Matrix<NX, NX, T> Ac_pad = detail::restrict_A(sys.A, Uc, nc);
    const Matrix<NY, NX, T> Cc_pad = detail::restrict_C(sys.C, Uc, nc);
    const auto              os = stability::detail::controllable_subspace(
        Matrix<NX, NX, T>(Ac_pad.transpose()), Matrix<NX, NY, T>(Cc_pad.transpose()), tol
    );
    if (!os.converged) {
        return MinimalRealizationResult<NX, NU, NY, T, NW, NV>{};
    }
    const size_t             no = os.dim;
    const Matrix<NX, NX, T>& Vo = os.basis;

    if (no == 0) {
        return MinimalRealizationResult<NX, NU, NY, T, NW, NV>{
            .system = StateSpace<NX, NU, NY, T, NW, NV>{
                .D = sys.D,
                .H = sys.H,
                .Ts = sys.Ts,
            },
            .order = 0,
            .eliminated = NX,
            .success = true,
        };
    }

    // Already minimal: return the original realization unchanged (identity basis).
    if (nc == NX && no == NX) {
        return MinimalRealizationResult<NX, NU, NY, T, NW, NV>{
            .system = sys,
            .order = NX,
            .eliminated = 0,
            .success = true,
        };
    }

    // CO basis in original coordinates: T = Uc * Vo(:, 0:no)
    // (Vo's significant right vectors live in the controllable coordinates.)
    Matrix<NX, NX, T> Tbasis{};
    for (size_t j = 0; j < no; ++j) {
        for (size_t i = 0; i < NX; ++i) {
            T sum = T{0};
            for (size_t k = 0; k < NX; ++k) {
                sum += Uc(i, k) * Vo(k, j);
            }
            Tbasis(i, j) = sum;
        }
    }

    // Project: Ar = T' A T, Br = T' B, Cr = C T (leading no block)
    if constexpr (NW > 0) {
        return MinimalRealizationResult<NX, NU, NY, T, NW, NV>{
            .system = StateSpace<NX, NU, NY, T, NW, NV>{
                .A = detail::restrict_A(sys.A, Tbasis, no),
                .B = detail::mul_basis_cols_t(Tbasis, no, sys.B),
                .C = detail::restrict_C(sys.C, Tbasis, no),
                .D = sys.D,
                .G = detail::mul_basis_cols_t(Tbasis, no, sys.G),
                .H = sys.H,
                .Ts = sys.Ts,
            },
            .order = no,
            .eliminated = NX - no,
            .success = true,
        };
    } else {
        return MinimalRealizationResult<NX, NU, NY, T, NW, NV>{
            .system = StateSpace<NX, NU, NY, T, NW, NV>{
                .A = detail::restrict_A(sys.A, Tbasis, no),
                .B = detail::mul_basis_cols_t(Tbasis, no, sys.B),
                .C = detail::restrict_C(sys.C, Tbasis, no),
                .D = sys.D,
                .H = sys.H,
                .Ts = sys.Ts,
            },
            .order = no,
            .eliminated = NX - no,
            .success = true,
        };
    }
}

/**
 * @brief Descriptive alias for @ref minreal
 * @see minreal()
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr MinimalRealizationResult<NX, NU, NY, T, NW, NV>
minimal_realization(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T tol = damp::sqrt(eps<T>())) {
    return minreal(sys, tol);
}

} // namespace design
} // namespace damp
