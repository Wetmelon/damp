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
 * Gramians, balanced truncation (`balred`), and singular-perturbation
 * residualization (`modred` with DC matching). All fixed-size and stack-only.
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
 * constexpr auto bal = design::balreal(sys);
 * static_assert(bal.success);
 * constexpr auto red = design::balred<1>(sys);
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
#include "damp/design/stability.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace design {

/**
 * @brief Method for eliminating states in modred
 *
 * Truncate drops states outright. MatchDC residualizes them (singular
 * perturbation) so the reduced model matches the full-order DC gain.
 */
enum class ModelReductionMethod {
    Truncate, ///< Drop states; D unchanged (balanced truncation)
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
 * @see balreal(), hankelsv(), balred(), modred()
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
struct BalancedRealizationResult {
    StateSpace<NX, NU, NY, T, NW, NV> balanced{};     ///< Balanced state-space model
    ColVec<NX, T>                     sigma{};        ///< Hankel singular values (descending)
    Matrix<NX, NX, T>                 transform{};    ///< Balancing transform: @f$ x = T x_b @f$
    bool                              success{false}; ///< true if Gramians, Cholesky, SVD all OK

    template<typename U>
    [[nodiscard]] constexpr BalancedRealizationResult<NX, NU, NY, U, NW, NV> as() const {
        return BalancedRealizationResult<NX, NU, NY, U, NW, NV>{
            balanced.template as<U>(),
            sigma.template as<U>(),
            transform.template as<U>(),
            success
        };
    }
};

/**
 * @brief Reduced-order model from balanced truncation or residualization
 *
 * @tparam NR Reduced state dimension
 *
 * @note Compare with MATLAB®'s `balred` / `modred` output.
 * @see balred(), modred(), BalancedRealizationResult
 */
template<size_t NR, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
struct BalancedReductionResult {
    StateSpace<NR, NU, NY, T, NW, NV> reduced{};         ///< Reduced-order model
    ColVec<NR, T>                     sigma_kept{};      ///< Hankel SV of retained states
    T                                 error_bound{T{0}}; ///< @f$ 2\sum \sigma_{\mathrm{discarded}} @f$ (truncation a-priori ‖·‖_∞ bound)
    bool                              success{false};    ///< true if balancing and reduction succeeded

    template<typename U>
    [[nodiscard]] constexpr BalancedReductionResult<NR, NU, NY, U, NW, NV> as() const {
        return BalancedReductionResult<NR, NU, NY, U, NW, NV>{
            reduced.template as<U>(),
            sigma_kept.template as<U>(),
            static_cast<U>(error_bound),
            success
        };
    }
};

namespace detail {

/// Symmetrize (X + Xᵀ)/2 to clean Lyapunov residual asymmetry before Cholesky.
template<size_t N, typename T>
[[nodiscard]] constexpr Matrix<N, N, T> symmetrize(const Matrix<N, N, T>& X) {
    return static_cast<T>(0.5) * (X + X.transpose());
}

/// Copy a fixed-size block out of a parent matrix into an owning Matrix.
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

/**
 * @brief Square-root factors and Hankel singular values shared by balreal / hankelsv.
 */
template<size_t NX, size_t NU, size_t NY, typename T>
struct SquareRootBalance {
    Matrix<NX, NX, T> Lc{};    ///< Wc = Lc Lcᵀ (lower Cholesky)
    Matrix<NX, NX, T> Lo{};    ///< Wo = Lo Loᵀ
    Matrix<NX, NX, T> U{};     ///< left singular vectors of Loᵀ Lc
    Matrix<NX, NX, T> V{};     ///< right singular vectors of Loᵀ Lc
    ColVec<NX, T>     sigma{}; ///< Hankel singular values
    bool              ok{false};
};

template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr SquareRootBalance<NX, NU, NY, T>
square_root_balance_factors(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    SquareRootBalance<NX, NU, NY, T> out{};
    static_assert(NX > 0, "square_root_balance_factors requires NX > 0");

    const bool discrete = sys.is_discrete();
    const auto Wc_opt = stability::controllability_gramian(sys.A, sys.B, discrete);
    const auto Wo_opt = stability::observability_gramian(sys.A, sys.C, discrete);
    if (!Wc_opt || !Wo_opt) {
        return out;
    }

    const auto Lc_opt = mat::cholesky(symmetrize(*Wc_opt));
    const auto Lo_opt = mat::cholesky(symmetrize(*Wo_opt));
    if (!Lc_opt || !Lo_opt) {
        return out;
    }
    out.Lc = *Lc_opt;
    out.Lo = *Lo_opt;

    // Loᵀ Lc = U Σ Vᵀ  →  σ_i are the Hankel singular values
    const Matrix<NX, NX, T> product = out.Lo.transpose() * out.Lc;
    const auto              svd_r = mat::svd(product);
    if (!svd_r.converged) {
        return out;
    }
    out.U = svd_r.singular_U;
    out.V = svd_r.singular_V;
    for (size_t i = 0; i < NX; ++i) {
        out.sigma(i) = static_cast<T>(svd_r.singular_values[i]);
    }

    // Reject non-positive singular values (non-minimal / numerically singular).
    const T tol = default_tol<T>();
    const T smax = out.sigma(0);
    if (!(smax > T{0})) {
        return out;
    }
    for (size_t i = 0; i < NX; ++i) {
        if (!(out.sigma(i) > tol * smax)) {
            return out;
        }
    }

    out.ok = true;
    return out;
}

/// Balancing transforms: T = Lc V Σ^{-1/2}, Tinv = Σ^{-1/2} Uᵀ Loᵀ.
template<size_t NX, typename T>
struct BalancingTransforms {
    Matrix<NX, NX, T> Tmat{};
    Matrix<NX, NX, T> Tinv{};
};

template<size_t NX, typename T>
[[nodiscard]] constexpr BalancingTransforms<NX, T> make_balancing_transforms(
    const Matrix<NX, NX, T>& Lc,
    const Matrix<NX, NX, T>& Lo,
    const Matrix<NX, NX, T>& U,
    const Matrix<NX, NX, T>& V,
    const ColVec<NX, T>&     sigma
) {
    Matrix<NX, NX, T> Si{};
    for (size_t i = 0; i < NX; ++i) {
        Si(i, i) = T{1} / damp::sqrt(sigma(i));
    }
    return BalancingTransforms<NX, T>{
        Lc * V * Si,
        Si * U.transpose() * Lo.transpose()
    };
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
 * Computes controllability/observability Gramians, their Cholesky factors
 * @f$ W_c = L_c L_c^\top @f$, @f$ W_o = L_o L_o^\top @f$, the SVD
 * @f$ L_o^\top L_c = U \Sigma V^\top @f$, and the balancing transform
 * @f$ T = L_c V \Sigma^{-1/2} @f$ so that in coordinates @f$ x = T x_b @f$
 * both Gramians equal @f$ \mathrm{diag}(\sigma) @f$.
 *
 * Continuous or discrete is selected from @p sys.Ts. Requires a stable,
 * minimal realization (positive-definite Gramians).
 *
 * @note Compare with MATLAB®'s `[sysb, g] = balreal(sys)`.
 * @see hankelsv(), balred(), modred(), minreal()
 * @see Moore (1981); Laub et al., "Computation of system balancing
 *      transformations," IEEE TAC 32(2), 1987
 *
 * @param sys Stable state-space system (continuous if Ts == 0, else discrete)
 * @return Balanced model, Hankel singular values, transform T, success flag
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr BalancedRealizationResult<NX, NU, NY, T, NW, NV>
balreal(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    const auto factors = detail::square_root_balance_factors(sys);
    if (!factors.ok) {
        return BalancedRealizationResult<NX, NU, NY, T, NW, NV>{};
    }

    const auto tf = detail::make_balancing_transforms(
        factors.Lc, factors.Lo, factors.U, factors.V, factors.sigma
    );

    if constexpr (NW > 0) {
        return BalancedRealizationResult<NX, NU, NY, T, NW, NV>{
            .balanced = StateSpace<NX, NU, NY, T, NW, NV>{
                .A = tf.Tinv * sys.A * tf.Tmat,
                .B = tf.Tinv * sys.B,
                .C = sys.C * tf.Tmat,
                .D = sys.D,
                .G = tf.Tinv * sys.G,
                .H = sys.H,
                .Ts = sys.Ts,
            },
            .sigma = factors.sigma,
            .transform = tf.Tmat,
            .success = true,
        };
    } else {
        return BalancedRealizationResult<NX, NU, NY, T, NW, NV>{
            .balanced = StateSpace<NX, NU, NY, T, NW, NV>{
                .A = tf.Tinv * sys.A * tf.Tmat,
                .B = tf.Tinv * sys.B,
                .C = sys.C * tf.Tmat,
                .D = sys.D,
                .H = sys.H,
                .Ts = sys.Ts,
            },
            .sigma = factors.sigma,
            .transform = tf.Tmat,
            .success = true,
        };
    }
}

/**
 * @brief Hankel singular values of a stable state-space system
 *
 * Reuses the square-root Gramian path of balreal without forming the
 * balanced realization. Values are non-increasing and non-negative.
 *
 * @note Compare with MATLAB®'s `hsv = hankelsv(sys)` / `hsvd(sys)`.
 * @see balreal(), balred()
 *
 * @param sys Stable state-space system
 * @return σ vector, or damp::nullopt if Gramians/Cholesky/SVD fail
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr damp::optional<ColVec<NX, T>>
hankelsv(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    const auto factors = detail::square_root_balance_factors(sys);
    if (!factors.ok) {
        return damp::nullopt;
    }
    return factors.sigma;
}

/**
 * @brief Balanced truncation to @p NR states
 *
 * Balances @p sys then drops the states with the smallest Hankel singular
 * values. The a-priori error bound on the result is
 * @f$ \lVert G - G_r \rVert_\infty \le 2 \sum_{k=\mathrm{NR}+1}^{n} \sigma_k @f$
 * (stored in @c error_bound).
 *
 * @tparam NR Reduced order (1 … NX)
 * @note Compare with MATLAB®'s `balred(sys, NR)` (truncation).
 * @see balreal(), modred(), hankelsv()
 *
 * @param sys Stable full-order system
 * @return Reduced model, kept σ, error bound, success
 */
template<
    size_t NR,
    size_t NX,
    size_t NU,
    size_t NY,
    size_t NW = 0,
    size_t NV = 0,
    typename T = double>
[[nodiscard]] constexpr BalancedReductionResult<NR, NU, NY, T, NW, NV>
balred(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    static_assert(NR > 0 && NR <= NX, "balred: require 0 < NR <= NX");

    const auto bal = balreal(sys);
    if (!bal.success) {
        return BalancedReductionResult<NR, NU, NY, T, NW, NV>{};
    }

    ColVec<NR, T> sigma_kept{};
    for (size_t i = 0; i < NR; ++i) {
        sigma_kept(i) = bal.sigma(i);
    }
    return BalancedReductionResult<NR, NU, NY, T, NW, NV>{
        .reduced = detail::truncate_states<NR>(bal.balanced),
        .sigma_kept = sigma_kept,
        .error_bound = detail::discarded_error_bound(bal.sigma, NR),
        .success = true,
    };
}

/**
 * @brief Model reduction by truncation or DC-matched residualization
 *
 * Balances @p sys (ordering states by descending Hankel singular value) then
 * eliminates the trailing @f$ n - \mathrm{NR} @f$ states:
 * - @ref ModelReductionMethod::Truncate — same as balred
 * - @ref ModelReductionMethod::MatchDC — singular perturbation so DC gain of
 *   the reduced model matches the full-order system
 *
 * @tparam NR Number of states to keep
 * @note Compare with MATLAB®'s `modred(sys, elim, 'MatchDC'|'Truncate')`
 *       (elim = states NR…n−1 after balancing).
 * @see balred(), balreal()
 *
 * @param sys    Stable full-order system
 * @param method Truncate or MatchDC (default MatchDC)
 * @return Reduced model; @c error_bound is the truncation bound (informative for MatchDC)
 */
template<
    size_t NR,
    size_t NX,
    size_t NU,
    size_t NY,
    size_t NW = 0,
    size_t NV = 0,
    typename T = double>
[[nodiscard]] constexpr BalancedReductionResult<NR, NU, NY, T, NW, NV>
modred(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    ModelReductionMethod                     method = ModelReductionMethod::MatchDC
) {
    static_assert(NR > 0 && NR <= NX, "modred: require 0 < NR <= NX");

    const auto bal = balreal(sys);
    if (!bal.success) {
        return BalancedReductionResult<NR, NU, NY, T, NW, NV>{};
    }

    ColVec<NR, T> sigma_kept{};
    for (size_t i = 0; i < NR; ++i) {
        sigma_kept(i) = bal.sigma(i);
    }
    const T error_bound = detail::discarded_error_bound(bal.sigma, NR);

    if (method == ModelReductionMethod::Truncate) {
        return BalancedReductionResult<NR, NU, NY, T, NW, NV>{
            .reduced = detail::truncate_states<NR>(bal.balanced),
            .sigma_kept = sigma_kept,
            .error_bound = error_bound,
            .success = true,
        };
    }

    const auto red = detail::residualize_states<NR>(bal.balanced);
    if (!red) {
        return BalancedReductionResult<NR, NU, NY, T, NW, NV>{
            .sigma_kept = sigma_kept,
            .error_bound = error_bound,
        };
    }
    return BalancedReductionResult<NR, NU, NY, T, NW, NV>{
        .reduced = *red,
        .sigma_kept = sigma_kept,
        .error_bound = error_bound,
        .success = true,
    };
}

// ---------------------------------------------------------------------------
// Descriptive aliases (AGENTS.md primary-name style)
// ---------------------------------------------------------------------------

/**
 * @brief Descriptive alias for balreal
 * @see balreal()
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr auto
balanced_realization(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    return balreal(sys);
}

/**
 * @brief Descriptive alias for balred
 * @see balred()
 */
template<
    size_t NR,
    size_t NX,
    size_t NU,
    size_t NY,
    size_t NW = 0,
    size_t NV = 0,
    typename T = double>
[[nodiscard]] constexpr auto
balanced_truncation(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    return balred<NR>(sys);
}

/**
 * @brief Descriptive alias for @ref hankelsv
 * @see hankelsv()
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr auto
hankel_singular_values(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    return hankelsv(sys);
}

} // namespace design
} // namespace damp
