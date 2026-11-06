// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file norms.hpp
 * @brief DC gain and H2 / H∞ norms for state-space systems
 */

#include <cstddef>
#include <limits>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/design/stability.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/concat.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/functions.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "damp/matrix/solve.hpp"
#include "damp/matrix/svd.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace analysis {

/**
 * @brief Compute DC gain of a continuous-time system
 *
 * DC gain = C * (-A)^{-1} * B + D (for continuous-time systems)
 * DC gain = C * (I - A)^{-1} * B + D (for discrete-time systems)
 *
 * @return DC gain matrix, or nullopt if A is singular at the evaluation point
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr damp::optional<Matrix<NY, NU, T>> dcgain(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys
) noexcept {
    Matrix<NX, NX, T> M;
    if (sys.is_continuous()) {
        // G(0) = C * A^{-1} * B + D, but we want -A^{-1} since G(s) = C(sI-A)^{-1}B + D at s=0
        M = -sys.A;
    } else {
        // G(1) = C * (I - A)^{-1} * B + D
        M = Matrix<NX, NX, T>::identity() - sys.A;
    }
    // M⁻¹·B by solving M·X = B (more accurate than forming the inverse).
    auto X = mat::solve(M, sys.B);
    if (!X) {
        return damp::nullopt;
    }
    return sys.C * (*X) + sys.D;
}

/**
 * @brief H2 norm of a state-space system.
 *
 * @f$ \lVert G \rVert_2^2 = \mathrm{tr}(C W_c C^\top) @f$ (continuous) or
 * @f$ \mathrm{tr}(C W_c C^\top + D D^\top) @f$ (discrete), where @f$ W_c @f$ is the
 * controllability Gramian. The continuous H2 norm is finite only for a strictly
 * proper system (@f$ D = 0 @f$); a nonzero @f$ D @f$ makes it unbounded and this
 * returns damp::nullopt (the project builds under @c -ffinite-math-only, so we
 * signal "infinite" with nullopt rather than a @c +inf sentinel).
 *
 * @note Compare with MATLAB®'s @c norm(sys,2). Requires @f$ A @f$ stable.
 * @see norm_hinf(), stability::controllability_gramian().
 *
 * @return ‖G‖₂, or damp::nullopt if the Gramian solve fails or
 *         the continuous norm is infinite (D != 0).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr damp::optional<T> norm_h2(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys
) {
    const bool discrete = sys.is_discrete();

    // Nonzero feedthrough: infinite in continuous time, finite (DᵀD term) in discrete.
    T d_sumsq = T{0};
    for (size_t i = 0; i < NY; ++i) {
        for (size_t j = 0; j < NU; ++j) {
            d_sumsq += sys.D(i, j) * sys.D(i, j);
        }
    }
    if (!discrete && d_sumsq > T{0}) {
        return damp::nullopt; // ‖G‖₂ = ∞ for a non-strictly-proper continuous system.
    }

    const auto Wc = stability::controllability_gramian(sys.A, sys.B, discrete);
    if (!Wc) {
        return damp::nullopt;
    }

    // trace(C Wc Cᵀ) = Σ over outputs of (row_i(C) · Wc · row_i(C)ᵀ).
    const Matrix<NY, NX, T> CW = sys.C * Wc.value();
    T                       sumsq = T{0};
    for (size_t i = 0; i < NY; ++i) {
        for (size_t k = 0; k < NX; ++k) {
            sumsq += CW(i, k) * sys.C(i, k);
        }
    }
    if (discrete) {
        sumsq += d_sumsq;
    }
    return damp::sqrt(sumsq);
}

/**
 * @brief H∞ norm of a stable system and the frequency of its peak
 *
 * The true norm lies in [gain, (1 + 2·tol)·gain], and σ̄(G) equals @c gain at @c frequency.
 */
template<typename T = double>
struct [[nodiscard]] HinfNormResult {
    T    gain{0};        ///< ‖G‖∞ to the relative tolerance (a lower bound attained at frequency)
    T    frequency{0};   ///< [rad/s] Where σ̄(G) peaks; the largest finite T when the peak is σ̄(D), approached as ω → ∞ (MATLAB®'s Inf)
    bool success{false}; ///< false: A not stable, an eigenvalue or Cholesky solve failed, or the iteration cap was reached
};

namespace detail {

/** @brief σ̄(M): |m| for a 1×1 matrix, the largest singular value otherwise. */
template<size_t NR, size_t NC, typename E>
[[nodiscard]] constexpr scalar_type_t<E> largest_singular_value(const Matrix<NR, NC, E>& M) {
    if constexpr (NR == 1 && NC == 1) {
        return damp::abs(M(0, 0));
    } else {
        return mat::svd(M).singular_values[0];
    }
}

/** @brief Peak of σ̄ over a continuous frequency axis, before mapping back to the caller's domain. */
template<typename T>
struct HinfPeak {
    T    gain{0};            ///< Best lower bound found
    T    omega{0};           ///< Continuous frequency where gain was evaluated
    bool at_infinity{false}; ///< gain is σ̄(D), approached as ω → ∞
    bool converged{false};   ///< The Hamiltonian test certified gain to the tolerance
};

/**
 * @brief Hamiltonian whose imaginary eigenvalues jω are the frequencies where γ is a singular value of G(jω)
 *
 * With @f$ R = \gamma^2 I - D^\top D @f$, positive definite for @f$ \gamma > \bar\sigma(D) @f$,
 * @f[
 *   H(\gamma) = \begin{bmatrix} A + B R^{-1} D^\top C & B R^{-1} B^\top \\
 *               -C^\top (I + D R^{-1} D^\top) C & -(A + B R^{-1} D^\top C)^\top \end{bmatrix}.
 * @f]
 * R⁻¹ is applied by two Cholesky solves; nullopt when R is not numerically positive definite.
 *
 * @see "A Bisection Method for Computing the H∞ Norm of a Transfer Matrix and Related Problems"
 *      (Boyd, Balakrishnan & Kabamba, 1989), Mathematics of Control, Signals and Systems 2(3), Theorem 1
 */
template<size_t NX, size_t NU, size_t NY, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<2 * NX, 2 * NX, T>> hinf_hamiltonian(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NY, NX, T>& C,
    const Matrix<NY, NU, T>& D,
    T                        gamma
) {
    const Matrix<NU, NU, T> R = (Matrix<NU, NU, T>::identity() * (gamma * gamma)) - (D.transpose() * D);
    const Matrix<NU, NX, T> DtC = D.transpose() * C;
    const auto              RDtC = mat::cholesky_solve(R, DtC);
    const auto              RBt = mat::cholesky_solve(R, B.transpose());
    if (!RDtC || !RBt) {
        return damp::nullopt;
    }
    const Matrix<NX, NX, T> Ah = A + (B * (*RDtC));
    const Matrix<NX, NX, T> Gh = B * (*RBt);
    const Matrix<NX, NX, T> Qh = (C.transpose() * C) + (DtC.transpose() * (*RDtC));
    return vcat(hcat(Ah, Gh), hcat(-Qh, -Ah.transpose()));
}

/** @brief G ≡ D: every Markov parameter C Aᵏ B, k < NX, is exactly zero (Cayley–Hamilton). */
template<size_t NX, size_t NU, size_t NY, typename T>
[[nodiscard]] constexpr bool zero_markov_parameters(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B, const Matrix<NY, NX, T>& C) {
    Matrix<NX, NU, T> AkB = B;
    for (size_t k = 0; k < NX; ++k) {
        const Matrix<NY, NU, T> M = C * AkB;
        for (size_t i = 0; i < NY; ++i) {
            for (size_t j = 0; j < NU; ++j) {
                if (M(i, j) != T{0}) {
                    return false;
                }
            }
        }
        AkB = A * AkB;
    }
    return true;
}

/**
 * @brief Bruinsma–Steinbuch peak search on a stable continuous realization (A, B, C, D)
 *
 * @p sigma_at(ω) returns σ̄ at the continuous frequency ω and @p sigma_infinity its limit
 * σ̄(D); both may come from an equivalent system (the discrete original of a bilinear map).
 * @p seeds are the pole magnitudes.
 *
 * Each pass forms H(γ) at γ = (1 + 2·tol)·gain. With no eigenvalue within
 * √ε·‖H‖₁ of the imaginary axis, γ bounds the norm and the search stops. Otherwise the
 * sorted |Im λ| are the frequencies where a singular value of G crosses γ, σ̄ is evaluated
 * at the midpoint of each pair of neighbours, and the largest value is the new gain. The
 * axis test is wide: an eigenvalue taken for a crossing only adds a midpoint, and a
 * midpoint inside an interval where σ̄ > γ still lies inside it. So a pass whose midpoints
 * do not beat γ also certifies the bound (in exact arithmetic), and stops the search at the
 * numerical floor where the computed crossings no longer separate.
 */
template<size_t NX, size_t NU, size_t NY, typename T, typename SigmaAt>
[[nodiscard]] constexpr HinfPeak<T> hinf_peak(
    const Matrix<NX, NX, T>&  A,
    const Matrix<NX, NU, T>&  B,
    const Matrix<NY, NX, T>&  C,
    const Matrix<NY, NU, T>&  D,
    const damp::array<T, NX>& seeds,
    const SigmaAt&            sigma_at,
    T                         sigma_infinity,
    T                         tol
) {
    constexpr size_t max_iterations = 30;
    HinfPeak<T>      peak{.gain = sigma_at(T{0})};
    for (size_t i = 0; i < NX; ++i) {
        const T s = sigma_at(seeds[i]);
        if (s > peak.gain) {
            peak.gain = s;
            peak.omega = seeds[i];
        }
    }
    if (sigma_infinity > peak.gain) {
        peak.gain = sigma_infinity;
        peak.at_infinity = true;
    }
    if (!(peak.gain > T{0})) {
        peak.converged = zero_markov_parameters(A, B, C); // zero at every seed: G ≡ 0, or no start
        return peak;
    }

    for (size_t it = 0; it < max_iterations; ++it) {
        const T    gamma = (T{1} + (T{2} * tol)) * peak.gain;
        const auto H = hinf_hamiltonian(A, B, C, D, gamma);
        if (!H) {
            return peak;
        }
        const auto eig = mat::compute_eigenvalues(*H);
        if (!eig.success) {
            return peak;
        }

        //! Crossing frequencies, one per conjugate pair, sorted ascending.
        const T                axis_tol = damp::sqrt(eps<T>()) * mat::one_norm(*H);
        damp::array<T, 2 * NX> crossing{};
        size_t                 count = 0;
        for (size_t i = 0; i < 2 * NX; ++i) {
            const auto& lambda = eig.values[i];
            if (lambda.imag() < T{0} || damp::abs(lambda.real()) > axis_tol) {
                continue;
            }
            size_t k = count++;
            for (; k > 0 && crossing[k - 1] > lambda.imag(); --k) {
                crossing[k] = crossing[k - 1];
            }
            crossing[k] = lambda.imag();
        }
        if (count == 0) {
            peak.converged = true;
            return peak;
        }

        T            best = T{0};
        T            best_omega = T{0};
        const size_t n_mid = (count > 1) ? count - 1 : 1;
        for (size_t k = 0; k < n_mid; ++k) {
            const T w = (count > 1) ? (crossing[k] + crossing[k + 1]) / T{2} : crossing[0];
            const T s = sigma_at(w);
            if (s > best) {
                best = s;
                best_omega = w;
            }
        }
        if (best > peak.gain) {
            peak.gain = best;
            peak.omega = best_omega;
            peak.at_infinity = false;
        }
        if (!(best > gamma)) {
            peak.converged = true;
            return peak;
        }
    }
    return peak;
}

} // namespace detail

/**
 * @brief H∞ norm of a stable state-space system, sup σ̄(G) over the imaginary axis or the unit circle, and its frequency
 *
 * Computed to a relative tolerance by the Bruinsma–Steinbuch algorithm: a lower bound
 * @f$ \gamma_{lb} @f$ from DC, ω → ∞ and the pole magnitudes, then the imaginary-axis
 * eigenvalues of the Hamiltonian @f$ H(\gamma) @f$ at @f$ \gamma = (1 + 2\,\mathrm{tol})\gamma_{lb} @f$.
 * Those eigenvalues are the frequencies where a singular value of @f$ G(j\omega) @f$ equals
 * @f$ \gamma @f$; σ̄ at the midpoints between them raises the bound, and the iteration
 * converges quadratically. It stops when @f$ H(\gamma) @f$ has no imaginary eigenvalue, so
 * @f$ \gamma_{lb} \le \lVert G \rVert_\infty \le (1 + 2\,\mathrm{tol})\gamma_{lb} @f$, with a
 * cap of 30 passes.
 *
 * A discrete system is mapped to the continuous one with the same values on the boundary,
 * @f$ G_c(s) = G\!\left(\frac{1+s}{1-s}\right) @f$: @f$ A_c = (A+I)^{-1}(A-I) @f$,
 * @f$ B_c = \sqrt2 (A+I)^{-1}B @f$, @f$ C_c = \sqrt2\,C(A+I)^{-1} @f$,
 * @f$ D_c = D - C(A+I)^{-1}B @f$, by solves with @f$ A + I @f$ (nonsingular for a stable A).
 * The unit circle @f$ z = e^{j\theta} @f$ maps to @f$ s = j\tan(\theta/2) @f$, so the norms
 * agree, and σ̄ is evaluated on the discrete G at @f$ \theta = 2\arctan\omega_c @f$.
 *
 * An unstable or marginally stable A (Re λ ≥ 0, or |λ| ≥ 1 for a discrete system) has no
 * H∞ norm: @c success = false, as norm_h2() returns nullopt. So does a @p tol that is not
 * positive.
 *
 * @note Compare with MATLAB®'s [ninf, fpeak] = hinfnorm(sys, tol) / norm(sys, Inf);
 *       matlab::hinfnorm is the alias.
 * @see "A Fast Algorithm to Compute the H∞-Norm of a Transfer Function Matrix"
 *      (Bruinsma & Steinbuch, 1990), Systems & Control Letters 14(4)
 * @see "A Regularity Result for the Singular Values of a Transfer Matrix and a Quadratically
 *      Convergent Algorithm for Computing its L∞-Norm" (Boyd & Balakrishnan, 1990),
 *      Systems & Control Letters 15(1)
 * @see norm_h2(), sigma()
 *
 * @param sys State-space system (continuous or discrete)
 * @param tol Relative tolerance on the norm, > 0 (default 1e-12, 1e-6 for float)
 * @return Peak gain, its frequency [rad/s] and whether it converged
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr HinfNormResult<T> norm_hinf(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    std::type_identity_t<T>                  tol = default_tol<T>()
) {
    using Cplx = damp::complex<T>;
    HinfNormResult<T> out{};
    if (!(tol > T{0})) {
        return out;
    }
    if constexpr (NX == 0) {
        out.gain = detail::largest_singular_value(sys.D);
        out.success = true;
        return out;
    } else {
        const bool discrete = sys.is_discrete();
        const auto eig = mat::compute_eigenvalues(sys.A);
        if (!eig.success) {
            return out;
        }
        for (size_t i = 0; i < NX; ++i) {
            const bool stable = discrete ? (damp::abs(eig.values[i]) < T{1}) : (eig.values[i].real() < T{0});
            if (!stable) {
                return out;
            }
        }
        const auto sigma_at_s = [&](Cplx s_or_z) -> T {
            const auto G = eval_frf(sys, s_or_z);
            return G ? detail::largest_singular_value(*G) : T{0}; // nullopt only on a pole, which a stable A keeps off the boundary
        };

        if (!discrete) {
            damp::array<T, NX> seeds{};
            for (size_t i = 0; i < NX; ++i) {
                seeds[i] = eig.values[i].abs();
            }
            const auto peak = detail::hinf_peak(
                sys.A, sys.B, sys.C, sys.D, seeds, [&](T w) { return sigma_at_s(Cplx{T{0}, w}); },
                detail::largest_singular_value(sys.D), tol
            );
            out.gain = peak.gain;
            out.frequency = peak.at_infinity ? std::numeric_limits<T>::max() : peak.omega;
            out.success = peak.converged;
            return out;
        }

        //! Bilinear map to the continuous equivalent, by solves with A + I.
        const Matrix<NX, NX, T> I = Matrix<NX, NX, T>::identity();
        const Matrix<NX, NX, T> ApI = sys.A + I;
        const auto              Ac = mat::solve(ApI, Matrix<NX, NX, T>(sys.A - I));
        const auto              X = mat::solve(ApI, sys.B);
        const auto              Yt = mat::solve(Matrix<NX, NX, T>(ApI.transpose()), Matrix<NX, NY, T>(sys.C.transpose()));
        if (!Ac || !X || !Yt) {
            return out;
        }
        const T                 root2 = damp::sqrt(T{2});
        const Matrix<NX, NU, T> Bc = (*X) * root2;
        const Matrix<NY, NX, T> Cc = Yt->transpose() * root2;
        const Matrix<NY, NU, T> Dc = sys.D - (sys.C * (*X));
        damp::array<T, NX>      seeds{};
        const Cplx              one{T{1}, T{0}};
        for (size_t i = 0; i < NX; ++i) {
            seeds[i] = ((eig.values[i] - one) / (eig.values[i] + one)).abs(); // |s| of the mapped pole
        }
        const auto theta_of = [](T w) { return T{2} * damp::atan(w); };
        const auto peak = detail::hinf_peak(
            *Ac, Bc, Cc, Dc, seeds,
            [&](T w) {
            const auto [s, c] = damp::sincos(theta_of(w));
            return sigma_at_s(Cplx{c, s});
        },
            sigma_at_s(Cplx{T{-1}, T{0}}), tol
        );
        out.gain = peak.gain;
        out.frequency = (peak.at_infinity ? damp::numbers::pi_v<T> : theta_of(peak.omega)) / sys.Ts;
        out.success = peak.converged;
        return out;
    }
}

} // namespace analysis
} // namespace damp
