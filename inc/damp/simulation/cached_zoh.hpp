// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file cached_zoh.hpp
 * @brief Cached ZOH maps for continuous LTI plants (expm once, step many times)
 *
 * Exact recomputes the augmented matrix exponential on every call. When the
 * sample period @c h is fixed and the input is held (ZOH) — the usual case for
 * dense plot samples between control updates — form Ad = e^{A h} and
 * Bd = ∫₀ʰ e^{Aτ} B dτ once, then advance with Discrete:
 *
 *     x⁺ = Ad x + Bd u
 *
 * Identical to one Exact step of size @c h; far cheaper when @c h does not change.
 *
 * @code
 *   StateSpace sys{.A = ..., .B = ..., .C = ..., .Ts = 0};
 *   auto zoh = CachedZoh<2, 1>::from(sys, 1e-4);  // expm once
 *   ColVec<2> x = x0;
 *   ColVec<1> u{1.0};  // held
 *   for (int k = 0; k < 1000; ++k) {
 *       // sample (t, x, C*x + D*u) for plotting
 *       x = zoh.step(x, u);
 *   }
 * @endcode
 *
 * @see Exact, Discrete, discretize(..., DiscretizationMethod::ZOH), analysis::lsim
 * @see "Feedback Control of Dynamic Systems" (Franklin et al., 2015), §8.3
 */

#include <cstddef>
#include <vector>

#include "damp/backend.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/functions.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#include "integrator.hpp"
#include "simulate.hpp"

namespace damp::sim {

/**
 * @brief Cached discrete ZOH maps Ad, Bd for step size @c h
 *
 * Built from continuous (A, B) via the same augmented expm as Exact:
 *
 *     [ Ad  Bd ]   [ A  B ] h
 *     [  0   I ] = e
 *
 * @tparam NX Number of states
 * @tparam NU Number of inputs
 * @tparam T  Scalar type
 */
template<size_t NX, size_t NU, typename T = double>
struct CachedZoh {
    Matrix<NX, NX, T> Ad{}; ///< Discrete A: e^{A h}
    Matrix<NX, NU, T> Bd{}; ///< Discrete B: ∫₀ʰ e^{A τ} B dτ
    T                 h{T{0}};

    /**
     * @brief Build Ad, Bd from continuous A, B and step @p step
     *
     * @param A     Continuous system matrix (NX × NX)
     * @param B     Continuous input matrix (NX × NU)
     * @param step  [s] ZOH sample period (must be > 0)
     */
    [[nodiscard]] static CachedZoh from(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B, T step) {
        CachedZoh z{};
        z.h = step;
        if (!(step > T{0})) {
            z.Ad = Matrix<NX, NX, T>::identity();
            z.Bd = Matrix<NX, NU, T>::zeros();
            return z;
        }

        // Same augmented expm as Exact::evolve — matches exact_lti_step for this h.
        Matrix<NX + NU, NX + NU, T> M = Matrix<NX + NU, NX + NU, T>::zeros();
        M.template block<NX, NX>(0, 0) = A;
        M.template block<NX, NU>(0, NX) = B;
        const auto expM = mat::expm(M * step);
        z.Ad = expM.template block<NX, NX>(0, 0);
        z.Bd = expM.template block<NX, NU>(0, NX);
        return z;
    }

    /**
     * @brief Build from a continuous StateSpace (uses sys.A, sys.B only)
     */
    template<size_t NY, size_t NW, size_t NV>
    [[nodiscard]] static CachedZoh from(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T step) {
        return from(sys.A, sys.B, step);
    }

    /**
     * @brief One discrete ZOH step: x⁺ = Ad x + Bd u
     *
     * @note Uses @ref Discrete (no integration — pure discrete map).
     */
    [[nodiscard]] ColVec<NX, T> step(const ColVec<NX, T>& x, const ColVec<NU, T>& u) const {
        return Discrete<NX, T>{}.evolve(Ad, Bd, x, u).x;
    }

    /**
     * @brief Discrete StateSpace for this step size (C, D from continuous plant)
     *
     * @param C  Output matrix
     * @param D  Feedthrough matrix
     * @return StateSpace with Ts = h
     */
    template<size_t NY>
    [[nodiscard]] StateSpace<NX, NU, NY, T> as_discrete_ss(
        const Matrix<NY, NX, T>& C,
        const Matrix<NY, NU, T>& D = {}
    ) const {
        return StateSpace<NX, NU, NY, T>{
            .A = Ad,
            .B = Bd,
            .C = C,
            .D = D,
            .Ts = h,
        };
    }
};

/**
 * @brief Simulate a continuous LTI plant with fixed ZOH step and held input
 *
 * Forms @ref CachedZoh once, then records @p n_steps + 1 samples (initial plus
 * each post-step state). Input @p u is held constant for the whole run — the
 * dense-plotting case when the command does not change between samples.
 *
 * @param sys     Continuous plant (Ts == 0)
 * @param x0      Initial state
 * @param u       Held input
 * @param h       [s] Sample / plot period
 * @param n_steps Number of discrete steps
 * @param t0      [s] Start time for the time vector (default 0)
 * @return SimulationResult with t, x, y, u histories
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] SimulationResult<NX, NU, NY, T> simulate_cached_zoh(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const ColVec<NX, T>&                     x0,
    const ColVec<NU, T>&                     u,
    T                                        h,
    size_t                                   n_steps,
    T                                        t0 = T{0}
) {
    const auto                      zoh = CachedZoh<NX, NU, T>::from(sys, h);
    SimulationResult<NX, NU, NY, T> sim;
    sim.t.reserve(n_steps + 1);
    sim.x.reserve(n_steps + 1);
    sim.y.reserve(n_steps + 1);
    sim.u.reserve(n_steps + 1);

    ColVec<NX, T> x = x0;
    T             t = t0;
    for (size_t k = 0; k <= n_steps; ++k) {
        sim.t.push_back(t);
        sim.x.push_back(x);
        sim.y.push_back(sys.C * x + sys.D * u);
        sim.u.push_back(u);
        if (k < n_steps) {
            x = zoh.step(x, u);
            t += h;
        }
    }
    return sim;
}

/**
 * @brief Simulate continuous LTI over @p t_span with fixed sample period @p h
 *
 * Same as the n_steps overload; step count is derived from the span. The final
 * sample is at or just before @c tf (partial last step is not taken — use
 * Exact / @ref exact_lti_step for a remainder shorter than @p h).
 *
 * @param sys    Continuous plant
 * @param x0     Initial state
 * @param u      Held input
 * @param t_span {t0, tf}
 * @param h      [s] Sample period
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] SimulationResult<NX, NU, NY, T> simulate_cached_zoh(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const ColVec<NX, T>&                     x0,
    const ColVec<NU, T>&                     u,
    const damp::pair<T, T>&                  t_span,
    T                                        h
) {
    const T span = t_span.second - t_span.first;
    if (!(h > T{0}) || !(span > T{0})) {
        SimulationResult<NX, NU, NY, T> sim;
        sim.t.push_back(t_span.first);
        sim.x.push_back(x0);
        sim.y.push_back(sys.C * x0 + sys.D * u);
        sim.u.push_back(u);
        return sim;
    }
    const size_t n_steps = static_cast<size_t>(span / h + static_cast<T>(0.5)); // nearest count of full steps
    return simulate_cached_zoh(sys, x0, u, h, n_steps, t_span.first);
}

/**
 * @brief Multi-rate ZOH: held control samples, dense plant samples for plotting
 *
 * For each control interval @p k, holds @p u_ctrl[k] and advances the plant with
 * the cached map @p n_fine times (plot stride @c h = control_period / n_fine).
 * Records every fine sample. @p u_ctrl.size() control holds are applied.
 *
 * @param sys      Continuous plant
 * @param x0       Initial state
 * @param u_ctrl   One input vector per control period (held for n_fine fine steps)
 * @param h_fine   [s] Fine sample period (must match CachedZoh construction)
 * @param n_fine   Fine steps per control hold ( ≥ 1 )
 * @param t0       [s] Start time
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] SimulationResult<NX, NU, NY, T> simulate_cached_zoh_multirate(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const ColVec<NX, T>&                     x0,
    const std::vector<ColVec<NU, T>>&        u_ctrl,
    T                                        h_fine,
    size_t                                   n_fine,
    T                                        t0 = T{0}
) {
    if (n_fine == 0) {
        n_fine = 1;
    }
    const auto                      zoh = CachedZoh<NX, NU, T>::from(sys, h_fine);
    const size_t                    n_total = u_ctrl.size() * n_fine;
    SimulationResult<NX, NU, NY, T> sim;
    sim.t.reserve(n_total + 1);
    sim.x.reserve(n_total + 1);
    sim.y.reserve(n_total + 1);
    sim.u.reserve(n_total + 1);

    ColVec<NX, T> x = x0;
    T             t = t0;
    ColVec<NU, T> u0 = u_ctrl.empty() ? ColVec<NU, T>{} : u_ctrl.front();
    sim.t.push_back(t);
    sim.x.push_back(x);
    sim.y.push_back(sys.C * x + sys.D * u0);
    sim.u.push_back(u0);

    for (size_t k = 0; k < u_ctrl.size(); ++k) {
        const ColVec<NU, T>& u = u_ctrl[k];
        for (size_t m = 0; m < n_fine; ++m) {
            x = zoh.step(x, u);
            t += h_fine;
            sim.t.push_back(t);
            sim.x.push_back(x);
            sim.y.push_back(sys.C * x + sys.D * u);
            sim.u.push_back(u);
        }
    }
    return sim;
}

} // namespace damp::sim
