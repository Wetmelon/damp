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
 * (balanced_reduction / eliminate_states).
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
 * @see balanced_realization(), balanced_reduction() for Hankel-based order reduction; eliminate_states() to eliminate chosen states
 * @see stability::controllability_matrix(), stability::observability_matrix()
 * @see Kailath, "Linear Systems" (1980), §6.4 (Kalman decomposition)
 * @see Chen, "Linear System Theory and Design" (3rd ed.), §6.4–6.5
 */


#include "damp/backend.hpp"
#include "damp/design/stability.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "damp/systems/realization.hpp"
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
 * @see minreal()
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
struct [[nodiscard]] MinimalRealizationResult {
    StateSpace<NX, NU, NY, T, NW, NV> system{};                            ///< Reduced model (leading @c order states)
    size_t                            order{0};                            ///< Dimension of controllable ∩ observable subspace
    size_t                            eliminated{0};                       ///< @c NX − @c order
    bool                              success{false};                      ///< status == DesignStatus::ok
    DesignStatus                      status{DesignStatus::not_converged}; ///< ok, or why the reduction is refused

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The same reduced model in the scalar @p U. */
    template<typename U>
    [[nodiscard]] constexpr MinimalRealizationResult<NX, NU, NY, U, NW, NV> as() const {
        return MinimalRealizationResult<NX, NU, NY, U, NW, NV>{
            system.template as<U>(),
            order,
            eliminated,
            success,
            status
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

/** @brief Y = M(:, 0:ncols)' * X  (left-multiply by basis columns transpose). */
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

/** @brief S = T(:,0:n)' * A * T(:,0:n) stored in the leading n×n of an N×N matrix. */
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

/**
 * @brief [B G], each block divided by its Frobenius norm (when nonzero): the inputs the
 *        controllable staircase runs on. Scaling a block leaves its range unchanged.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr Matrix<NX, NU + NW, T> input_blocks(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    if constexpr (NW == 0) {
        return sys.B;
    } else {
        const T                b_norm = sys.B.norm();
        const T                g_norm = sys.G.norm();
        const T                b_scale = (b_norm > T{0}) ? T{1} / b_norm : T{1};
        const T                g_scale = (g_norm > T{0}) ? T{1} / g_norm : T{1};
        Matrix<NX, NU + NW, T> out{};
        for (size_t r = 0; r < NX; ++r) {
            for (size_t c = 0; c < NU; ++c) {
                out(r, c) = sys.B(r, c) * b_scale;
            }
            for (size_t c = 0; c < NW; ++c) {
                out(r, NU + c) = sys.G(r, c) * g_scale;
            }
        }
        return out;
    }
}

/** @brief C_r = C * T(:,0:n) stored in the leading columns of NY×N. */
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
 * Projects @p sys onto its controllable-and-observable subspace, found by a block-Arnoldi
 * staircase (stability::detail::controllable_subspace). Stair 0 is range(B); each later
 * stair is A, scaled by 1/‖A‖∞, applied to the previous stair's new directions, then
 * projected off the basis found so far by classical Gram–Schmidt run twice. An SVD of
 * each stair sets its rank: the left singular vectors whose singular values exceed
 * @p tol (tol·σ_max(B) for stair 0) extend the orthonormal basis, and the staircase stops
 * at a stair with none. The same staircase on the dual @f$ (A_c^\top, C_c^\top) @f$ of the
 * controllable part gives its observable subspace. Continuous and discrete systems are
 * treated the same (a structural property of (A, B, C)).
 *
 * The states are balanced first (balance_states: D⁻¹AD, D⁻¹B, CD with D = diag(2ᵉ) over
 * [A B; C 0], exact in floating point), so the verdict does not depend on the units of each
 * state: a pressure in Pa next to a position in m couples through entries 10¹³ apart, and
 * unbalanced that coupling reads as a cancellation. The staircase on (A, B) runs on
 * [B G] when the model has noise inputs, so a state only the process noise drives (the
 * disturbance state of an offset-free augmentation) is kept; each block is scaled to unit
 * norm first, so a large noise gain cannot hide a small input. MIMO and SISO are
 * both supported. A mode reached only through couplings at or below the tolerance,
 * relative to the balanced plant's fastest rate ‖A‖∞, is treated as cancelled.
 *
 * Refused with DesignStatus::bad_weights when an entry of the model is not finite, and
 * with DesignStatus::not_converged when an SVD does not converge.
 *
 * @note This is not balanced truncation. Use balanced_reduction to drop weakly coupled
 *       but still controllable and observable modes, or eliminate_states to eliminate
 *       chosen fast states.
 *
 * @note Compare with MATLAB®'s minreal(sys, tol).
 *
 * @see balanced_realization(), balanced_reduction(), eliminate_states()
 * @see stability::is_controllable(), stability::is_observable()
 * @see "The Generalized Eigenstructure Problem in Linear System Theory" (Van Dooren, 1981), IEEE TAC 26(1)
 *
 * @tparam NX Number of states
 * @tparam NU Number of inputs
 * @tparam NY Number of outputs
 * @tparam T  Scalar type (default double)
 * @tparam NW Process-noise columns (inputs to the staircase, then passed through the same basis)
 * @tparam NV Measurement-noise columns (unchanged)
 * @param sys State-space model (continuous if Ts == 0, else discrete)
 * @param tol Relative tolerance on each stair's singular values (default √eps, as MATLAB®'s minreal)
 * @return MinimalRealizationResult with reduced model, order, success and status
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr MinimalRealizationResult<NX, NU, NY, T, NW, NV>
minreal(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T tol = damp::sqrt(eps<T>())) {
    using Result = MinimalRealizationResult<NX, NU, NY, T, NW, NV>;
    if constexpr (NX == 0) {
        Result out{.system = sys};
        return ::damp::detail::set_status(out, DesignStatus::ok);
    }
    const auto finite = [](T v) { return damp::isfinite(v); };
    if (!mat::all_of(sys.A, finite) || !mat::all_of(sys.B, finite) || !mat::all_of(sys.C, finite) || !mat::all_of(sys.D, finite)
        || !mat::all_of(sys.G, finite) || !mat::all_of(sys.H, finite)) {
        return ::damp::detail::refused<Result>(DesignStatus::bad_weights);
    }

    // --- Balance the states: the staircases and the projection run on D⁻¹AD, D⁻¹B, CD, D⁻¹G --
    const StateSpace<NX, NU, NY, T, NW, NV> bal = balance_states(sys).sys;

    // --- Controllable subspace (staircase) of (A, [B G]) ------------------------
    const auto cs = stability::detail::controllable_subspace(bal.A, detail::input_blocks(bal), tol);
    if (!cs.converged) {
        return ::damp::detail::refused<Result>(DesignStatus::not_converged);
    }
    const size_t             nc = cs.dim;
    const Matrix<NX, NX, T>& Uc = cs.basis; // leading nc columns: orthonormal basis

    if (nc == 0) {
        // No controllable dynamics: I/O is pure feedthrough.
        Result out{.system = StateSpace<NX, NU, NY, T, NW, NV>{.D = sys.D, .H = sys.H, .Ts = sys.Ts}, .eliminated = NX};
        return ::damp::detail::set_status(out, DesignStatus::ok);
    }

    // --- Observable subspace of the controllable subsystem ------------------
    // (Ac, Cc) = (Ucᵀ A Uc, C Uc) in the leading nc coordinates; the dual staircase on
    // (Acᵀ, Ccᵀ) stays inside them because the padding is zero.
    const Matrix<NX, NX, T> Ac_pad = detail::restrict_A(bal.A, Uc, nc);
    const Matrix<NY, NX, T> Cc_pad = detail::restrict_C(bal.C, Uc, nc);
    const auto              os = stability::detail::controllable_subspace(
        Matrix<NX, NX, T>(Ac_pad.transpose()), Matrix<NX, NY, T>(Cc_pad.transpose()), tol
    );
    if (!os.converged) {
        return ::damp::detail::refused<Result>(DesignStatus::not_converged);
    }
    const size_t             no = os.dim;
    const Matrix<NX, NX, T>& Vo = os.basis;

    if (no == 0) {
        Result out{.system = StateSpace<NX, NU, NY, T, NW, NV>{.D = sys.D, .H = sys.H, .Ts = sys.Ts}, .eliminated = NX};
        return ::damp::detail::set_status(out, DesignStatus::ok);
    }

    // Already minimal: return the original realization unchanged (identity basis).
    if (nc == NX && no == NX) {
        Result out{.system = sys, .order = NX};
        return ::damp::detail::set_status(out, DesignStatus::ok);
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

    // Project the balanced model: Ar = T' A T, Br = T' B, Cr = C T (leading no block)
    Result out{
        .system = StateSpace<NX, NU, NY, T, NW, NV>{
            .A = detail::restrict_A(bal.A, Tbasis, no),
            .B = detail::mul_basis_cols_t(Tbasis, no, bal.B),
            .C = detail::restrict_C(bal.C, Tbasis, no),
            .D = sys.D,
            .H = sys.H,
            .Ts = sys.Ts,
        },
        .order = no,
        .eliminated = NX - no,
    };
    if constexpr (NW > 0) {
        out.system.G = detail::mul_basis_cols_t(Tbasis, no, bal.G);
    }
    return ::damp::detail::set_status(out, DesignStatus::ok);
}

} // namespace design
} // namespace damp
