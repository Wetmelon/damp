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
#include <numbers>
#include <vector>

#include "damp/backend.hpp"
#include "damp/design/stability.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/matrix.hpp"
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
 * @brief H∞ norm of a state-space system: sup_ω σ̄(G(jω)).
 *
 * The peak of the largest singular value of the frequency response over all
 * frequencies (continuous: @f$ s=j\omega @f$; discrete: @f$ z=e^{j\theta} @f$ over
 * the unit circle). Estimated by sampling a log-spaced grid seeded with each
 * pole's natural frequency (peaks sit near lightly-damped poles), taking the max.
 *
 * @note Compare with MATLAB®'s @c norm(sys,Inf) / @c hinfnorm. Requires @f$ A @f$
 *       stable (otherwise the H∞ norm is unbounded).
 * @note ponytail: frequency-sweep estimate, O(n_points) FRF evals. Pole-seeded so
 *       it catches resonant peaks; a needle-sharp peak between samples can still be
 *       under-reported. Swap to the exact Hamiltonian/Boyd–Balakrishnan bisection
 *       (imaginary-axis eigenvalue test, reusing the eigen machinery) if you need a
 *       guaranteed bound for robustness certification. Raise @p n_points meanwhile.
 *
 * @see norm_h2(), eval_frf().
 *
 * @param sys      State-space system (continuous or discrete), assumed stable.
 * @param n_points Number of log-spaced frequency samples (default 1024).
 * @return Estimated ‖G‖∞.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] T norm_hinf(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    size_t                                   n_points = 1024
) {
    using Cplx = damp::complex<T>;
    const bool discrete = sys.is_discrete();

    // σ̄(G) at a single evaluation point on the stability boundary.
    // Contour poles (eval_frf nullopt) mean ‖G‖∞ is unbounded; under finite-math
    // builds we report max() rather than fabricating +inf.
    const auto sigma_max = [&](Cplx s) -> T {
        const auto G_opt = eval_frf(sys, s);
        if (!G_opt) {
            return std::numeric_limits<T>::max();
        }
        const Matrix<NY, NU, Cplx>& G = *G_opt;
        if constexpr (NY == 1 && NU == 1) {
            return G(0, 0).abs(); // SISO: largest singular value is just |G|.
        } else {
            return mat::svd(G).singular_values[0];
        }
    };

    T peak = T{0};

    if (discrete) {
        // Sweep the upper unit circle z = e^{jθ}, θ ∈ (0, π].
        const T pi = std::numbers::pi_v<T>;
        for (size_t k = 1; k <= n_points; ++k) {
            const T theta = pi * static_cast<T>(k) / static_cast<T>(n_points);
            const auto [s, c] = damp::sincos(theta);
            peak = damp::max(peak, sigma_max(Cplx(c, s)));
        }
        return peak;
    }

    // Continuous: log-spaced ω grid spanning the pole frequencies, plus the pole
    // natural frequencies themselves as explicit samples (peaks live near them).
    const auto eig = mat::compute_eigenvalues(sys.A);
    T          w_min = std::numeric_limits<T>::max();
    T          w_max = T{0};
    for (size_t i = 0; i < NX; ++i) {
        const T w = eig.values[i].abs();
        if (w > T{0}) {
            w_min = damp::min(w_min, w);
            w_max = damp::max(w_max, w);
        }
        if (w > T{0}) {
            peak = damp::max(peak, sigma_max(Cplx(T{0}, w))); // seed at each |pole|.
        }
    }
    if (w_max == T{0}) { // pure integrator(s): pick a default decade window.
        w_min = static_cast<T>(1e-3);
        w_max = static_cast<T>(1e3);
    }
    const T lo = damp::log10(w_min) - T{3};
    const T hi = damp::log10(w_max) + T{3};
    for (size_t k = 0; k <= n_points; ++k) {
        const T ex = lo + (hi - lo) * static_cast<T>(k) / static_cast<T>(n_points);
        const T w = damp::pow(T{10}, ex);
        peak = damp::max(peak, sigma_max(Cplx(T{0}, w)));
    }
    return peak;
}

} // namespace analysis
} // namespace damp
