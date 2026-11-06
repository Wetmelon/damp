// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file model_reduction.hpp
 * @brief Balanced realization and model-order reduction (Moore / Glover)
 *
 * Square-root balanced realization from the controllability and observability
 * Gramians, balanced reduction (`balanced_reduction`, by truncation or DC-matched
 * residualization), and elimination of chosen states in the model's own
 * coordinates (`eliminate_states`). All fixed-size and stack-only. The MATLAB®
 * spellings balreal / hankelsv / balred / modred are in matlab::.
 *
 * Example: balance a second-order lag and truncate the weak mode
 * @code
 * #include "damp/design/model_reduction.hpp"
 * using namespace damp;
 *
 * constexpr StateSpace<2, 1, 1> sys{
 *     .A = Matrix<2, 2>{{-1.0, 0.0}, {0.0, -10.0}},
 *     .B = Matrix<2, 1>{{1.0}, {1.0}},
 *     .C = Matrix<1, 2>{{1.0, 0.1}},
 * };
 * constexpr auto bal = design::balanced_realization(sys);
 * static_assert(bal.success);
 * constexpr auto red = design::balanced_reduction<1>(sys);
 * static_assert(red.success);
 * // a-priori bound: ‖G − G_r‖_∞ ≤ red.error_bound = 2·Σ σ_discarded
 * @endcode
 *
 * @see Moore, "Principal Component Analysis in Linear Systems," IEEE TAC 26(1),
 *      1981, https://doi.org/10.1109/TAC.1981.1102568
 * @see Glover, "All optimal Hankel-norm approximations…," Int. J. Control 39(6),
 *      1984
 * @see Antoulas, Approximation of Large-Scale Dynamical Systems, SIAM, 2005
 */

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/design/lyapunov.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace design {

/**
 * @brief How balanced_reduction and eliminate_states eliminate states
 *
 * Truncate drops states outright. MatchDC residualizes them (singular
 * perturbation) so the reduced model matches the full-order DC gain.
 */
enum class ModelReductionMethod {
    Truncate, ///< Drop states; D unchanged
    MatchDC   ///< Residualize (ẋ₂ = 0 or x₂⁺ = x₂); exact DC gain
};

/**
 * @brief Balanced realization result (Moore square-root method)
 *
 * In balanced coordinates the Gramians equal @f$ W_c = W_o = \mathrm{diag}(\sigma) @f$
 * where @f$ \sigma @f$ are the Hankel singular values (descending). The original
 * state maps as @f$ x = T x_b @f$.
 *
 * @note Compare with MATLAB®'s `[sysb, g] = balreal(sys)`.
 * @see balanced_realization(), hankel_singular_values(), balanced_reduction(), eliminate_states()
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
struct [[nodiscard]] BalancedRealizationResult {
    StateSpace<NX, NU, NY, T, NW, NV> balanced{};                          ///< Balanced state-space model
    ColVec<NX, T>                     sigma{};                             ///< Hankel singular values (descending)
    Matrix<NX, NX, T>                 transform{};                         ///< Balancing transform: x = T x_b
    bool                              success{false};                      ///< status == DesignStatus::ok
    DesignStatus                      status{DesignStatus::not_converged}; ///< ok, or why the realization is refused

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The same realization in the scalar @p U. */
    template<typename U>
    [[nodiscard]] constexpr BalancedRealizationResult<NX, NU, NY, U, NW, NV> as() const {
        return BalancedRealizationResult<NX, NU, NY, U, NW, NV>{
            balanced.template as<U>(),
            sigma.template as<U>(),
            transform.template as<U>(),
            success,
            status
        };
    }
};

/**
 * @brief Reduced-order model from balanced truncation or residualization
 *
 * @tparam NR Reduced state dimension
 *
 * @note Compare with MATLAB®'s `balred` output.
 * @see balanced_reduction(), BalancedRealizationResult
 */
template<size_t NR, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
struct [[nodiscard]] BalancedReductionResult {
    StateSpace<NR, NU, NY, T, NW, NV> reduced{};                           ///< Reduced-order model
    ColVec<NR, T>                     sigma_kept{};                        ///< Hankel SV of retained states
    T                                 error_bound{T{0}};                   ///< 2 Σ σ_discarded (truncation a-priori ‖·‖∞ bound)
    bool                              success{false};                      ///< status == DesignStatus::ok
    DesignStatus                      status{DesignStatus::not_converged}; ///< ok, or why the reduction is refused

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The same reduced model in the scalar @p U. */
    template<typename U>
    [[nodiscard]] constexpr BalancedReductionResult<NR, NU, NY, U, NW, NV> as() const {
        return BalancedReductionResult<NR, NU, NY, U, NW, NV>{
            reduced.template as<U>(),
            sigma_kept.template as<U>(),
            static_cast<U>(error_bound),
            success,
            status
        };
    }
};

namespace detail {

/** @brief Copy a fixed-size block out of a parent matrix into an owning Matrix. */
template<size_t BR, size_t BC, size_t R, size_t C, typename T>
[[nodiscard]] constexpr Matrix<BR, BC, T>
take_block(const Matrix<R, C, T>& M, size_t r0, size_t c0) {
    Matrix<BR, BC, T> out{};
    if constexpr (BR > 0 && BC > 0) {
        for (size_t i = 0; i < BR; ++i) {
            for (size_t j = 0; j < BC; ++j) {
                out(i, j) = M(r0 + i, c0 + j);
            }
        }
    }
    return out;
}

/**
 * @brief Truncate a partitioned state-space to the leading NR states.
 */
template<
    size_t NR,
    size_t NX,
    size_t NU,
    size_t NY,
    size_t NW,
    size_t NV,
    typename T>
[[nodiscard]] constexpr StateSpace<NR, NU, NY, T, NW, NV>
truncate_states(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    static_assert(NR > 0 && NR <= NX, "truncate_states: 0 < NR <= NX");
    if constexpr (NW > 0) {
        return StateSpace<NR, NU, NY, T, NW, NV>{
            .A = take_block<NR, NR>(sys.A, 0, 0),
            .B = take_block<NR, NU>(sys.B, 0, 0),
            .C = take_block<NY, NR>(sys.C, 0, 0),
            .D = sys.D,
            .G = take_block<NR, NW>(sys.G, 0, 0),
            .H = sys.H,
            .Ts = sys.Ts,
        };
    } else {
        return StateSpace<NR, NU, NY, T, NW, NV>{
            .A = take_block<NR, NR>(sys.A, 0, 0),
            .B = take_block<NR, NU>(sys.B, 0, 0),
            .C = take_block<NY, NR>(sys.C, 0, 0),
            .D = sys.D,
            .H = sys.H,
            .Ts = sys.Ts,
        };
    }
}

/**
 * @brief Singular-perturbation residualization of the trailing NX−NR states.
 *
 * Continuous: set @f$ \dot x_2 = 0 @f$ and eliminate @f$ x_2 @f$.
 * Discrete: set @f$ x_2^+ = x_2 @f$ (steady state) and eliminate @f$ x_2 @f$.
 * Returns nullopt if the eliminated dynamics block is singular.
 */
template<
    size_t NR,
    size_t NX,
    size_t NU,
    size_t NY,
    size_t NW,
    size_t NV,
    typename T>
[[nodiscard]] constexpr damp::optional<StateSpace<NR, NU, NY, T, NW, NV>>
residualize_states(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    static_assert(NR > 0 && NR <= NX, "residualize_states: 0 < NR <= NX");

    if constexpr (NR == NX) {
        return sys;
    } else {
        constexpr size_t NE = NX - NR; // eliminated

        const auto A11 = take_block<NR, NR>(sys.A, 0, 0);
        const auto A12 = take_block<NR, NE>(sys.A, 0, NR);
        const auto A21 = take_block<NE, NR>(sys.A, NR, 0);
        const auto A22 = take_block<NE, NE>(sys.A, NR, NR);
        const auto B1 = take_block<NR, NU>(sys.B, 0, 0);
        const auto B2 = take_block<NE, NU>(sys.B, NR, 0);
        const auto C1 = take_block<NY, NR>(sys.C, 0, 0);
        const auto C2 = take_block<NY, NE>(sys.C, 0, NR);

        // M^{-1} [A21 B2 G2]  with M = A22 (cont., via −M^{-1}) or (I−A22) (disc.)
        Matrix<NE, NE, T> M{};
        T                 sign = T{1}; // discrete: +; continuous applied as − below
        if (sys.is_continuous()) {
            M = A22;
            sign = T{-1};
        } else {
            M = Matrix<NE, NE, T>::identity() - A22;
            sign = T{1};
        }

        const auto Minv_A21 = mat::solve(M, A21);
        if (!Minv_A21) {
            return damp::nullopt;
        }
        const auto Minv_B2 = mat::solve(M, B2);
        if (!Minv_B2) {
            return damp::nullopt;
        }

        const auto Ar = A11 + (sign * (A12 * (*Minv_A21)));
        const auto Br = B1 + (sign * (A12 * (*Minv_B2)));
        const auto Cr = C1 + (sign * (C2 * (*Minv_A21)));
        const auto Dr = sys.D + (sign * (C2 * (*Minv_B2)));

        if constexpr (NW > 0) {
            const auto G1 = take_block<NR, NW>(sys.G, 0, 0);
            const auto G2 = take_block<NE, NW>(sys.G, NR, 0);
            const auto Minv_G2 = mat::solve(M, G2);
            if (!Minv_G2) {
                return damp::nullopt;
            }
            return StateSpace<NR, NU, NY, T, NW, NV>{
                .A = Ar,
                .B = Br,
                .C = Cr,
                .D = Dr,
                .G = G1 + (sign * (A12 * (*Minv_G2))),
                .H = sys.H,
                .Ts = sys.Ts,
            };
        } else {
            return StateSpace<NR, NU, NY, T, NW, NV>{
                .A = Ar,
                .B = Br,
                .C = Cr,
                .D = Dr,
                .H = sys.H,
                .Ts = sys.Ts,
            };
        }
    }
}

/** @brief True when Elim lists at least one and fewer than NX distinct state indices, each below NX. */
template<size_t NX, size_t... Elim>
[[nodiscard]] constexpr bool valid_elimination() {
    if constexpr (sizeof...(Elim) == 0 || sizeof...(Elim) >= NX) {
        return false;
    } else {
        constexpr damp::array<size_t, sizeof...(Elim)> elim{Elim...};
        for (size_t i = 0; i < elim.size(); ++i) {
            if (elim[i] >= NX) {
                return false;
            }
            for (size_t j = 0; j < i; ++j) {
                if (elim[j] == elim[i]) {
                    return false;
                }
            }
        }
        return true;
    }
}

/** @brief State order with the kept states first (ascending) and the eliminated ones last. */
template<size_t NX, size_t... Elim>
[[nodiscard]] constexpr damp::array<size_t, NX> kept_then_eliminated() {
    constexpr damp::array<size_t, sizeof...(Elim)> elim{Elim...};
    damp::array<size_t, NX>                        order{};
    size_t                                         n = 0;
    for (size_t i = 0; i < NX; ++i) {
        if (!((Elim == i) || ...)) {
            order[n] = i;
            ++n;
        }
    }
    for (const size_t e : elim) {
        order[n] = e;
        ++n;
    }
    return order;
}

/** @brief The same system with state i of the result being state order[i] of @p sys. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StateSpace<NX, NU, NY, T, NW, NV>
permute_states(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const damp::array<size_t, NX>& order) {
    Matrix<NX, NX, T> P{}; // x = P x', exact in floating point
    for (size_t i = 0; i < NX; ++i) {
        P(order[i], i) = T{1};
    }
    auto out = sys;
    out.A = P.transpose() * sys.A * P;
    out.B = P.transpose() * sys.B;
    out.C = sys.C * P;
    if constexpr (NW > 0) {
        out.G = P.transpose() * sys.G;
    }
    return out;
}

/**
 * @brief Square-root factors and Hankel singular values shared by balanced_realization /
 *        hankel_singular_values / balanced_reduction.
 */
template<size_t NX, size_t NU, size_t NY, typename T>
struct SquareRootBalance {
    Matrix<NX, NX, T> Lc{};                                ///< Wc = Lc Lcᵀ (lower triangular, Hammarling)
    Matrix<NX, NX, T> Lo{};                                ///< Wo = Lo Loᵀ
    Matrix<NX, NX, T> U{};                                 ///< left singular vectors of Loᵀ Lc
    Matrix<NX, NX, T> V{};                                 ///< right singular vectors of Loᵀ Lc
    ColVec<NX, T>     sigma{};                             ///< Hankel singular values
    size_t            rank{0};                             ///< leading σᵢ above NX·64·ε·σ₁: the states Σ^{-1/2} can scale
    DesignStatus      status{DesignStatus::not_converged}; ///< ok when the factors and the SVD succeeded
};

template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr SquareRootBalance<NX, NU, NY, T>
square_root_balance_factors(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    SquareRootBalance<NX, NU, NY, T> out{};
    static_assert(NX > 0, "square_root_balance_factors requires NX > 0");
    const auto finite = [](T v) { return damp::isfinite(v); };
    if (!mat::all_of(sys.A, finite) || !mat::all_of(sys.B, finite) || !mat::all_of(sys.C, finite) || !mat::all_of(sys.D, finite)) {
        out.status = DesignStatus::bad_weights;
        return out;
    }

    // Hammarling factors Wc = Rcᵀ Rc, Wo = Roᵀ Ro, computed without forming either Gramian:
    // a Cholesky factorization of a computed Gramian loses its small eigenvalues, which are
    // the Hankel singular values this method exists to resolve.
    const bool discrete = sys.is_discrete();
    const auto Rc = discrete ? discrete_lyapunov_cholesky(sys.A, sys.B) : lyapunov_cholesky(sys.A, sys.B);
    const auto At = sys.A.transpose();
    const auto Ct = sys.C.transpose();
    const auto Ro = discrete ? discrete_lyapunov_cholesky(At, Ct) : lyapunov_cholesky(At, Ct);
    if (!Rc || !Ro) {
        return out;
    }
    out.Lc = Rc->transpose();
    out.Lo = Ro->transpose();

    // Loᵀ Lc = U Σ Vᵀ  →  σ_i are the Hankel singular values
    const Matrix<NX, NX, T> product = out.Lo.transpose() * out.Lc;
    const auto              svd_r = mat::svd(product);
    if (!svd_r.success) {
        return out;
    }
    out.U = svd_r.singular_U;
    out.V = svd_r.singular_V;
    for (size_t i = 0; i < NX; ++i) {
        out.sigma(i) = static_cast<T>(svd_r.singular_values[i]);
    }

    // The σᵢ are accurate to roundoff of σ₁ even when small (Hammarling factors, Jacobi SVD), but
    // Σ^{-1/2} in a balancing transform amplifies the noise in the singular vectors of a σᵢ lost
    // to roundoff: rank counts the leading states a transform can scale.
    const T tol = static_cast<T>(NX) * T{64} * eps<T>();
    const T smax = out.sigma(0);
    while (out.rank < NX && smax > T{0} && out.sigma(out.rank) > tol * smax) {
        ++out.rank;
    }
    out.status = DesignStatus::ok;
    return out;
}

/**
 * @brief The leading NR columns of the balancing transform and rows of its inverse:
 *        T = Lc V₁ Σ₁^{-1/2}, Tinv = Σ₁^{-1/2} U₁ᵀ Loᵀ over the NR largest σ
 *
 * NR = NX is the full balancing transform; NR < NX is the square-root truncation, which needs
 * only σ_NR above roundoff.
 */
template<size_t NR, size_t NX, typename T>
struct BalancingTransforms {
    Matrix<NX, NR, T> Tmat{}; ///< x = Tmat x_b
    Matrix<NR, NX, T> Tinv{}; ///< x_b = Tinv x
};

/** @copybrief BalancingTransforms */
template<size_t NR, size_t NX, size_t NU, size_t NY, typename T>
[[nodiscard]] constexpr BalancingTransforms<NR, NX, T> make_balancing_transforms(const SquareRootBalance<NX, NU, NY, T>& f) {
    BalancingTransforms<NR, NX, T> out{};
    const Matrix<NX, NX, T>        LcV = f.Lc * f.V;
    const Matrix<NX, NX, T>        LoU = f.Lo * f.U;
    for (size_t j = 0; j < NR; ++j) {
        const T si = T{1} / damp::sqrt(f.sigma(j));
        for (size_t i = 0; i < NX; ++i) {
            out.Tmat(i, j) = LcV(i, j) * si;
            out.Tinv(j, i) = LoU(i, j) * si;
        }
    }
    return out;
}

/** @brief The model in the coordinates x = T x_b of the leading NR balanced states. */
template<size_t NR, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StateSpace<NR, NU, NY, T, NW, NV>
project_balanced(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const BalancingTransforms<NR, NX, T>& tf) {
    StateSpace<NR, NU, NY, T, NW, NV> out{
        .A = tf.Tinv * sys.A * tf.Tmat,
        .B = tf.Tinv * sys.B,
        .C = sys.C * tf.Tmat,
        .D = sys.D,
        .H = sys.H,
        .Ts = sys.Ts,
    };
    if constexpr (NW > 0) {
        out.G = tf.Tinv * sys.G;
    }
    return out;
}

template<size_t NX, typename T>
[[nodiscard]] constexpr T discarded_error_bound(const ColVec<NX, T>& sigma, size_t keep) {
    T sum = T{0};
    for (size_t i = keep; i < NX; ++i) {
        sum += sigma(i);
    }
    return T{2} * sum;
}

} // namespace detail

/**
 * @brief Balanced realization via the square-root (Moore / Laub) method
 *
 * Computes the Cholesky factors of the controllability/observability Gramians,
 * @f$ W_c = L_c L_c^\top @f$, @f$ W_o = L_o L_o^\top @f$, directly by Hammarling's method
 * (lyapunov_cholesky; the Gramians themselves are never formed), the SVD
 * @f$ L_o^\top L_c = U \Sigma V^\top @f$, and the balancing transform
 * @f$ T = L_c V \Sigma^{-1/2} @f$ so that in coordinates @f$ x = T x_b @f$
 * both Gramians equal @f$ \mathrm{diag}(\sigma) @f$.
 *
 * Continuous or discrete is selected from @p sys.Ts. Requires a stable,
 * minimal realization (positive-definite Gramians). Refused with DesignStatus::bad_weights
 * for a non-finite model, DesignStatus::not_converged when the Gramian factors (an unstable
 * model) or the SVD fail, and DesignStatus::singular when a σᵢ is lost to roundoff of σ₁
 * (a non-minimal model: balanced_reduction() by truncation still reduces it).
 *
 * @note Compare with MATLAB®'s `[sysb, g] = balreal(sys)`.
 * @see hankel_singular_values(), balanced_reduction(), eliminate_states(), minreal()
 * @see Moore (1981); Laub et al., "Computation of system balancing
 *      transformations," IEEE TAC 32(2), 1987
 *
 * @param sys Stable state-space system (continuous if Ts == 0, else discrete)
 * @return Balanced model, Hankel singular values, transform T, success and status
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr BalancedRealizationResult<NX, NU, NY, T, NW, NV>
balanced_realization(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    using Result = BalancedRealizationResult<NX, NU, NY, T, NW, NV>;
    const auto factors = detail::square_root_balance_factors(sys);
    if (factors.status != DesignStatus::ok) {
        return ::damp::detail::refused<Result>(factors.status);
    }
    if (factors.rank < NX) {
        return ::damp::detail::refused<Result>(DesignStatus::singular);
    }
    const auto tf = detail::make_balancing_transforms<NX>(factors);
    Result     out{.balanced = detail::project_balanced(sys, tf), .sigma = factors.sigma, .transform = tf.Tmat};
    return ::damp::detail::set_status(out, DesignStatus::ok);
}

/**
 * @brief Hankel singular values of a stable state-space system
 *
 * Reuses the square-root Gramian path of balanced_realization without forming the
 * balanced realization. Values are non-increasing and non-negative, small ones
 * included: a non-minimal model has σᵢ at roundoff of σ₁.
 *
 * @note Compare with MATLAB®'s `hsv = hankelsv(sys)` / `hsvd(sys)`.
 * @see balanced_realization(), balanced_reduction()
 *
 * @param sys Stable state-space system
 * @return σ vector, or damp::nullopt for a non-finite model or when the Gramian factors
 *         (an unstable model) or the SVD fail
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr damp::optional<ColVec<NX, T>>
hankel_singular_values(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    const auto factors = detail::square_root_balance_factors(sys);
    if (factors.status != DesignStatus::ok) {
        return damp::nullopt;
    }
    return factors.sigma;
}

/**
 * @brief Balanced reduction to @p NR states, by truncation or DC-matched residualization
 *
 * Balances @p sys, which orders the states by descending Hankel singular value,
 * then eliminates the trailing @f$ n - \mathrm{NR} @f$ states:
 * - @ref ModelReductionMethod::Truncate (default) drops them. The a-priori error bound is
 *   @f$ \lVert G - G_r \rVert_\infty \le 2 \sum_{k=\mathrm{NR}+1}^{n} \sigma_k @f$
 *   (stored in @c error_bound). Only the leading NR singular vectors are used (square-root
 *   truncation), so σ_NR must be above roundoff of σ₁ and the discarded σ need not be: a
 *   non-minimal model reduces.
 * - @ref ModelReductionMethod::MatchDC residualizes them (singular perturbation), so
 *   the reduced DC gain equals the full-order one. It needs the full balanced
 *   realization, so every σᵢ above roundoff (balanced_realization()).
 *
 * Balanced states mix the original ones. To eliminate chosen states and keep the
 * others physical, use eliminate_states().
 *
 * Refused as balanced_realization(), plus DesignStatus::singular when σ_NR is lost to
 * roundoff or MatchDC meets a singular A₂₂ (I − A₂₂ discrete).
 *
 * @tparam NR Reduced order (1 … NX)
 * @note Compare with MATLAB®'s `balred(sys, NR, balredOptions('StateElimMethod', …))`,
 *       whose default is MatchDC.
 * @see balanced_realization(), eliminate_states(), hankel_singular_values()
 * @see "Singular perturbation approximation of balanced systems" (Liu & Anderson, 1989)
 * @see "Truncated balanced realization of a stable non-minimal state-space system" (Tombs & Postlethwaite, 1987)
 *
 * @param sys    Stable full-order system
 * @param method Truncate (default) or MatchDC
 * @return Reduced model, kept σ, error bound (informative for MatchDC), success and status
 */
template<
    size_t NR,
    size_t NX,
    size_t NU,
    size_t NY,
    size_t NW = 0,
    size_t NV = 0,
    typename T = double>
[[nodiscard]] constexpr BalancedReductionResult<NR, NU, NY, T, NW, NV> balanced_reduction(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    ModelReductionMethod                     method = ModelReductionMethod::Truncate
) {
    static_assert(NR > 0 && NR <= NX, "balanced_reduction: require 0 < NR <= NX");
    using Result = BalancedReductionResult<NR, NU, NY, T, NW, NV>;

    const auto factors = detail::square_root_balance_factors(sys);
    if (factors.status != DesignStatus::ok) {
        return ::damp::detail::refused<Result>(factors.status);
    }
    const size_t needed = (method == ModelReductionMethod::Truncate) ? NR : NX;
    if (factors.rank < needed) {
        return ::damp::detail::refused<Result>(DesignStatus::singular);
    }

    Result out{};
    for (size_t i = 0; i < NR; ++i) {
        out.sigma_kept(i) = factors.sigma(i);
    }
    out.error_bound = detail::discarded_error_bound(factors.sigma, NR);

    if (method == ModelReductionMethod::Truncate) {
        out.reduced = detail::project_balanced(sys, detail::make_balancing_transforms<NR>(factors));
        return ::damp::detail::set_status(out, DesignStatus::ok);
    }
    const auto red = detail::residualize_states<NR>(detail::project_balanced(sys, detail::make_balancing_transforms<NX>(factors)));
    if (!red) {
        return ::damp::detail::set_status(out, DesignStatus::singular);
    }
    out.reduced = *red;
    return ::damp::detail::set_status(out, DesignStatus::ok);
}

/**
 * @brief Eliminate the listed states, by DC-matched residualization or truncation
 *
 * Works in the coordinates of @p sys, so the kept states stay physical (currents,
 * pressures, positions) and the reduced model still matches the sensors, weights and
 * observers written for them. The reduced states are the kept ones in their original
 * order.
 *
 * With @f$ x_1 @f$ kept and @f$ x_2 @f$ eliminated:
 * - @ref ModelReductionMethod::MatchDC (default) sets @f$ \dot x_2 = 0 @f$ and
 *   substitutes the solution for @f$ x_2 @f$, so the DC gain equals the full-order one:
 *   @f[ A_r = A_{11} - A_{12} A_{22}^{-1} A_{21},\quad B_r = B_1 - A_{12} A_{22}^{-1} B_2, @f]
 *   @f[ C_r = C_1 - C_2 A_{22}^{-1} A_{21},\quad D_r = D - C_2 A_{22}^{-1} B_2. @f]
 *   A discrete model sets @f$ x_2^+ = x_2 @f$, which replaces @f$ -A_{22}^{-1} @f$
 *   with @f$ (I - A_{22})^{-1} @f$.
 * - @ref ModelReductionMethod::Truncate drops @f$ x_2 @f$ with its rows and columns.
 *
 * The approximation holds when the eliminated states are stable and much faster than
 * the kept ones. That is the caller's judgment and is not checked.
 *
 * @code
 * // x = [i, ω]: drop the electrical pole of a DC motor, keep the speed
 * constexpr auto speed = design::eliminate_states<0>(motor);
 * static_assert(speed.has_value());
 * @endcode
 *
 * @tparam Elim Indices of the states to eliminate: distinct, each below NX, at least one state kept
 * @note Compare with MATLAB®'s `modred(sys, elim, 'MatchDC'|'Truncate')`; MATLAB® indices are 1-based.
 * @see balanced_reduction() to reduce by Hankel singular value instead
 * @see "Singular Perturbation Methods in Control: Analysis and Design" (Kokotović, Khalil & O'Reilly, 1986)
 *
 * @param sys    Full-order system
 * @param method MatchDC (default) or Truncate
 * @return Reduced model, or nullopt when MatchDC meets a singular A₂₂ (I − A₂₂ discrete),
 *         i.e. an integrator among the eliminated states
 */
template<size_t... Elim, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
    requires(detail::valid_elimination<NX, Elim...>())
[[nodiscard]] constexpr damp::optional<StateSpace<NX - sizeof...(Elim), NU, NY, T, NW, NV>> eliminate_states(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    ModelReductionMethod                     method = ModelReductionMethod::MatchDC
) {
    constexpr size_t NR = NX - sizeof...(Elim);
    const auto       reordered = detail::permute_states(sys, detail::kept_then_eliminated<NX, Elim...>());
    if (method == ModelReductionMethod::Truncate) {
        return detail::truncate_states<NR>(reordered);
    }
    return detail::residualize_states<NR>(reordered);
}

} // namespace design
} // namespace damp
