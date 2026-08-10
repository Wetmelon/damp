// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file poles.hpp
 * @brief Open-loop poles, damping/natural frequency, pole-zero maps, root locus
 */

#include <cstddef>
#include <limits>
#include <vector>

#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"


namespace damp {
namespace analysis {

/**
 * @brief Compute open-loop poles (eigenvalues of A matrix)
 *
 * Returns the eigenvalue vector from the Francis QR eigensolver when it
 * converges. If the solver does not converge, returns a zero vector (same
 * conservative posture as @ref is_stable_continuous, which treats
 * non-convergence as "not proven stable").
 *
 * @param A State matrix
 * @return Vector of pole locations as complex numbers, or zeros if eigen failed
 */
template<size_t NX, typename T = double>
[[nodiscard]] constexpr ColVec<NX, damp::complex<T>> poles(const Matrix<NX, NX, T>& A) {
    auto eigen = mat::compute_eigenvalues(A);
    if (!eigen.converged) {
        return ColVec<NX, damp::complex<T>>{};
    }
    return eigen.values;
}

/**
 * @brief Check continuous-time stability
 *
 * A continuous system is stable if all eigenvalues have Re(λ) < 0.
 *
 * @param A State matrix
 * @return true if all eigenvalues in the left half-plane
 */
template<size_t NX, typename T = double>
[[nodiscard]] constexpr bool is_stable_continuous(const Matrix<NX, NX, T>& A) {
    auto eigen = mat::compute_eigenvalues(A);
    if (!eigen.converged) {
        return false;
    }
    for (size_t i = 0; i < NX; ++i) {
        if (eigen.values[i].real() >= T{0}) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Natural frequency and damping ratio for each pole
 */
template<typename T = double>
struct PoleInfo {
    damp::complex<T> location{};
    T                natural_freq{};  ///< ωn = |pole| (rad/s)
    T                damping_ratio{}; ///< ζ = -Re(pole)/|pole|
    T                time_constant{}; ///< τ = -1/Re(pole) (seconds, for stable poles)
};

/**
 * @brief Compute natural frequency and damping for each pole
 *
 * @param A State matrix
 * @return Array of PoleInfo with ωn, ζ, τ for each pole
 */
template<size_t NX, typename T = double>
[[nodiscard]] constexpr damp::array<PoleInfo<T>, NX>
damp(const Matrix<NX, NX, T>& A) {
    auto                         eigen = mat::compute_eigenvalues(A);
    damp::array<PoleInfo<T>, NX> info{};
    for (size_t i = 0; i < NX; ++i) {
        auto p = eigen.values[i];
        T    wn = damp::abs(p);
        T    zeta = (wn > T{0}) ? -p.real() / wn : T{0};
        T    tau = (p.real() < T{0}) ? T{-1} / p.real() : std::numeric_limits<T>::max();
        info[i] = PoleInfo<T>{p, wn, zeta, tau};
    }
    return info;
}

// ============================================================================
// Pole-Zero Maps
// ============================================================================

/**
 * @brief Poles and zeros of a system, for pole-zero plotting
 *
 * Stored as runtime vectors of complex values (a pole-zero map feeds a scatter
 * plot, and the zero count is data-dependent). MATLAB® equivalent: the data
 * returned by `pzmap`.
 */
template<typename T = double>
struct PoleZeroMap {
    std::vector<damp::complex<T>> poles; ///< Pole locations
    std::vector<damp::complex<T>> zeros; ///< Zero locations
};

/**
 * @brief Roots of a polynomial given in ascending powers (MATLAB® `roots`, reversed order)
 *
 * For coefficients c[0] + c[1]·x + … + c[N-1]·x^{N-1}, returns the N-1 roots as
 * the eigenvalues of the companion matrix. The highest-order coefficient
 * c[N-1] must be nonzero (no trailing-zero padding).
 */
template<size_t N, typename T>
[[nodiscard]] std::vector<damp::complex<T>> poly_roots(const damp::array<T, N>& c) {
    std::vector<damp::complex<T>> out;
    if constexpr (N >= 2) {
        constexpr size_t M = N - 1;
        Matrix<M, M, T>  companion{};
        for (size_t i = 0; i + 1 < M; ++i) {
            companion(i, i + 1) = T{1};
        }
        for (size_t j = 0; j < M; ++j) {
            companion(M - 1, j) = -c[j] / c[N - 1];
        }
        const auto eig = mat::compute_eigenvalues(companion);
        out.reserve(M);
        for (size_t i = 0; i < M; ++i) {
            out.push_back(eig.values[i]);
        }
    }
    return out;
}

/**
 * @brief Pole-zero map of a SISO transfer function (MATLAB® `pzmap(tf)`)
 *
 * Poles are the roots of the denominator, zeros the roots of the numerator.
 */
template<size_t Nnum, size_t Nden, typename T>
[[nodiscard]] PoleZeroMap<T> pzmap(const TransferFunction<Nnum, Nden, T>& tf) {
    return {poly_roots(tf.den), poly_roots(tf.num)};
}

/**
 * @brief Pole map of a state matrix (MATLAB® `pzmap(sys)`, poles only)
 *
 * Returns the eigenvalues of A as poles. Transmission zeros are not computed
 * (they require a generalized/QZ eigensolver, not yet available).
 */
template<size_t NX, typename T>
[[nodiscard]] PoleZeroMap<T> pzmap(const Matrix<NX, NX, T>& A) {
    PoleZeroMap<T> r;
    const auto     eig = mat::compute_eigenvalues(A);
    r.poles.reserve(NX);
    for (size_t i = 0; i < NX; ++i) {
        r.poles.push_back(eig.values[i]);
    }
    return r;
}

/**
 * @brief Pole map of a state-space system (MATLAB® `pzmap(sys)`, poles only)
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] PoleZeroMap<T> pzmap(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    return pzmap(sys.A);
}

// ============================================================================
// Root locus (SISO closed-loop poles vs gain)
// ============================================================================

/**
 * @brief Root-locus data: closed-loop poles along a gain grid
 *
 * For unity negative feedback with scalar gain @f$ k @f$, the closed-loop
 * state matrix of a SISO plant @f$(A,B,C,D)@f$ is
 * @f[
 *   A_{\mathrm{cl}}(k) = A - \frac{k}{1 + k D}\, B C
 * @f]
 * (the usual @f$ A - kBC @f$ when @f$ D = 0 @f$). Each sample stores the
 * eigenvalues of @f$ A_{\mathrm{cl}}(k) @f$, ordered for continuity along the
 * gain grid by nearest-neighbor matching to the previous sample.
 *
 * @note Compare with MATLAB®'s rlocus / rlocusplot (minimal gain-grid form;
 *       no adaptive branch tracking beyond nearest-neighbor reordering).
 */
template<typename T = double>
struct RootLocusResult {
    std::vector<T>                             gains;           ///< Gain samples k
    std::vector<std::vector<damp::complex<T>>> poles;           ///< poles[i] = closed-loop poles at gains[i]
    std::vector<damp::complex<T>>              open_loop_poles; ///< Eigenvalues of A (k = 0)
    std::vector<damp::complex<T>>              open_loop_zeros; ///< Transmission zeros when available (TF path)
};

namespace detail {

/**
 * @brief Reorder @p curr so each entry is nearest (in C) to the matching @p prev entry
 *
 * Greedy sequential matching; good enough for smooth gain grids on modest NX.
 */
template<typename T>
[[nodiscard]] std::vector<damp::complex<T>> match_pole_order(
    const std::vector<damp::complex<T>>& prev,
    std::vector<damp::complex<T>>        curr
) {
    if (prev.size() != curr.size() || prev.empty()) {
        return curr;
    }
    std::vector<damp::complex<T>> ordered(curr.size());
    std::vector<bool>             used(curr.size(), false);
    for (size_t i = 0; i < prev.size(); ++i) {
        size_t best = 0;
        T      best_d = std::numeric_limits<T>::max();
        for (size_t j = 0; j < curr.size(); ++j) {
            if (used[j]) {
                continue;
            }
            const T d = damp::abs(curr[j] - prev[i]);
            if (d < best_d) {
                best_d = d;
                best = j;
            }
        }
        used[best] = true;
        ordered[i] = curr[best];
    }
    return ordered;
}

} // namespace detail

/**
 * @brief Root locus of a SISO state-space plant over an explicit gain grid
 *
 * Samples closed-loop poles of unity negative feedback @f$ u = -k y @f$ at each
 * gain in @p gains. Points where @f$ |1 + k D| @f$ is near zero (infinite
 * closed-loop pole from direct feedthrough) are skipped.
 *
 * @param sys   SISO plant (continuous or discrete; same A_cl formula)
 * @param gains Gain samples (typically non-negative, increasing)
 * @return RootLocusResult with per-gain pole vectors and open-loop poles
 *
 * @note Compare with MATLAB®'s r = rlocus(sys, k).
 */
template<size_t NX, size_t NW, size_t NV, typename T>
[[nodiscard]] RootLocusResult<T> rlocus(
    const StateSpace<NX, 1, 1, T, NW, NV>& sys,
    const std::vector<T>&                  gains
) {
    RootLocusResult<T> result;
    result.gains.reserve(gains.size());
    result.poles.reserve(gains.size());

    // Open-loop poles
    {
        const auto eig = mat::compute_eigenvalues(sys.A);
        result.open_loop_poles.reserve(NX);
        for (size_t i = 0; i < NX; ++i) {
            result.open_loop_poles.push_back(eig.values[i]);
        }
    }

    const T d = sys.D(0, 0);
    const T eps = damp::default_tol<T>() * (T{1} + damp::abs(d));

    std::vector<damp::complex<T>> prev;
    for (const T k : gains) {
        const T den = T{1} + k * d;
        if (damp::abs(den) < eps) {
            continue; // direct-feedthrough singularity
        }
        const T scale = k / den;
        // A_cl = A - scale * B * C  (outer product B C is NX×NX)
        Matrix<NX, NX, T> A_cl = sys.A;
        for (size_t i = 0; i < NX; ++i) {
            for (size_t j = 0; j < NX; ++j) {
                A_cl(i, j) = A_cl(i, j) - scale * sys.B(i, 0) * sys.C(0, j);
            }
        }

        const auto eig = mat::compute_eigenvalues(A_cl);
        if (!eig.converged) {
            continue;
        }

        std::vector<damp::complex<T>> poles_k;
        poles_k.reserve(NX);
        for (size_t i = 0; i < NX; ++i) {
            poles_k.push_back(eig.values[i]);
        }
        if (!prev.empty()) {
            poles_k = detail::match_pole_order(prev, damp::move(poles_k));
        }
        prev = poles_k;
        result.gains.push_back(k);
        result.poles.push_back(damp::move(poles_k));
    }

    return result;
}

/**
 * @brief Root locus of a SISO transfer function over an explicit gain grid
 *
 * Realizes the TF in companion form and delegates to the state-space path.
 * Open-loop zeros are filled from the numerator roots. Returns an empty
 * result (zeros still filled when possible) if companion realization fails
 * (zero leading denominator coefficient).
 *
 * @note Compare with MATLAB®'s r = rlocus(sys, k).
 */
template<size_t Nnum, size_t Nden, typename T>
[[nodiscard]] RootLocusResult<T> rlocus(
    const TransferFunction<Nnum, Nden, T>& tf,
    const std::vector<T>&                  gains
) {
    RootLocusResult<T> result;
    result.open_loop_zeros = poly_roots(tf.num);
    const auto ss = tf.to_state_space();
    if (!ss) {
        return result;
    }
    auto ss_rl = rlocus(*ss, gains);
    result.gains = damp::move(ss_rl.gains);
    result.poles = damp::move(ss_rl.poles);
    result.open_loop_poles = damp::move(ss_rl.open_loop_poles);
    return result;
}

} // namespace analysis
} // namespace damp
