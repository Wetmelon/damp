// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file mhe.hpp
 * @brief Moving Horizon Estimation: constrained sliding-window state estimation
 *        (`design::mhe` + `damp::MHE`, roadmap #15).
 *
 * The estimation counterpart of MPC: instead of the Kalman filter's recursive
 * update, each tick solves a weighted least-squares problem over the last
 * NH + 1 window states @f$ z = [x_0; \dots; x_{NH}] @f$,
 * @f[
 *   J = \lVert x_0 - \bar x \rVert^2_{P_\infty^{-1}}
 *     + \sum_{k=0}^{NH-1} \lVert x_{k+1} - A x_k - B u_k \rVert^2_{Q^{-1}}
 *     + \sum_{k=0}^{NH} \lVert y_k - C x_k \rVert^2_{R^{-1}},
 * @f]
 * subject to box constraints on every window state. The dynamics residuals are
 * the process noises, the measurement residuals the sensor noises, and the
 * arrival term anchors the window start to the prior mean @f$ \bar x @f$ with
 * the steady-state Kalman covariance @f$ P_\infty @f$ as its weight — on a
 * linear, unconstrained problem the window endpoint therefore reproduces the
 * steady-state Kalman filter estimate.
 *
 * Relative to a Kalman filter, MHE can enforce state/input constraints on the
 * estimate (e.g. non-negative level, pressure, or state of charge via
 * `x_min = 0`). Use MHE when those constraints matter; keep the Kalman filter
 * when they do not.
 *
 * The QP is fixed-size and allocation-free (Hessian and constraint rows are
 * constant; only the gradient changes with data), solved by the warm-started
 * active-set solver — at steady state the active constraint set barely changes
 * between ticks, so solves are a few iterations. Approximations vs. the exact
 * (full-information) estimator, both standard practice: the arrival weight is
 * held at the steady-state @f$ P_\infty^{-1} @f$ rather than propagated, and
 * the arrival mean is updated with the previous solve's smoothed estimate.
 *
 * Feedthrough is not supported (D = 0), matching the MPC. Bounds use the
 * `numeric_limits` sentinels (never infinity; -ffinite-math-only).
 *
 * @note Compare with MATLAB®'s nlmhe / nlmheMultistage restricted to the linear
 *       case, with mhestate's bounds as MHEConstraints.
 * @see estimation/kalman.hpp for the unconstrained recursive equivalent,
 *      design/qp.hpp for the solver, controllers/mpc.hpp for the control-side
 *      mirror of this formulation
 * @see Rao, Rawlings & Mayne, "Constrained State Estimation for Nonlinear
 *      Discrete-Time Systems: Stability and Moving Horizon Approximations,"
 *      IEEE TAC 48(2), 2003, https://doi.org/10.1109/TAC.2003.812777
 * @see Rawlings, Mayne & Diehl, "Model Predictive Control," 2nd ed., 2017, ch. 4
 *
 * @code
 * // Level estimator that can never report a negative level:
 * constexpr design::MHEConstraints<1> limits{.x_min = ColVec<1>{0.0}};
 * constexpr auto art = design::mhe<8>(tank, Q, R, limits);
 * static_assert(art.success);
 * MHE estimator{art.as<float>()};
 * // each tick: x_hat = estimator.update(y, u_previous);
 * @endcode
 */

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/design/qp.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/block.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/solve.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace design {

/**
 * @struct MHEConstraints
 * @brief Box constraints on the estimated states (applied to every window state)
 *
 * Defaults are unbounded (numeric_limits sentinels) — with no finite bounds the
 * estimator is an unconstrained least-squares smoother matching the Kalman
 * filter at the window endpoint.
 */
template<size_t NX, typename T = double>
struct MHEConstraints {
    ColVec<NX, T> x_min = detail::filled_vector<NX>(std::numeric_limits<T>::lowest()); ///< State lower bound
    ColVec<NX, T> x_max = detail::filled_vector<NX>(std::numeric_limits<T>::max());    ///< State upper bound

    template<typename U>
    [[nodiscard]] constexpr MHEConstraints<NX, U> as() const {
        return MHEConstraints<NX, U>{
            detail::convert_bounds<U>(x_min),
            detail::convert_bounds<U>(x_max),
        };
    }
};

/**
 * @struct MHEArtifacts
 * @brief Window QP data produced by mhe(), consumed by damp::MHE
 *
 * The Hessian and constraint rows are constant; the runtime rebuilds only the
 * gradient from its measurement/input ring buffers each tick using the stored
 * weight maps. Use .as<float>() for embedded deployment (bound sentinels
 * convert safely).
 */
template<size_t NX, size_t NU, size_t NY, size_t NH, typename T = double>
struct MHEArtifacts {
    static constexpr size_t NV = NX * (NH + 1); ///< Decision variables (window states)
    static constexpr size_t NIB = 2 * NV;       ///< Box-constraint rows

    Matrix<NV, NV, T>     H{};                 ///< Window QP Hessian (positive definite)
    Matrix<NIB, NV, T>    A_con{};             ///< Box rows: [I; −I]
    ColVec<NIB, T>        b_bound{};           ///< [x_max…; −x_min…] (constant; sentinel rows disabled)
    Matrix<NX, NX, T>     M_a{};               ///< 2·P∞⁻¹ (arrival gradient map)
    Matrix<NX, NY, T>     M_cy{};              ///< 2·CᵀR⁻¹ (measurement gradient map)
    Matrix<NX, NU, T>     M_qu{};              ///< 2·Q⁻¹B (process gradient map, next state)
    Matrix<NX, NU, T>     M_aqu{};             ///< 2·AᵀQ⁻¹B (process gradient map, current state)
    MHEConstraints<NX, T> constraints{};       ///< Original limits
    size_t                max_qp_iterations{}; ///< Per-tick QP iteration budget
    bool                  success{false};      ///< true if the synthesis validated

    template<typename U>
    [[nodiscard]] constexpr MHEArtifacts<NX, NU, NY, NH, U> as() const {
        return MHEArtifacts<NX, NU, NY, NH, U>{
            H.template as<U>(),
            A_con.template as<U>(),
            detail::convert_bounds<U>(b_bound),
            M_a.template as<U>(),
            M_cy.template as<U>(),
            M_qu.template as<U>(),
            M_aqu.template as<U>(),
            constraints.template as<U>(),
            max_qp_iterations,
            success,
        };
    }
};

/**
 * @brief Synthesize a constrained moving-horizon estimator
 *
 * Builds the constant window Hessian (block tridiagonal from the dynamics and
 * measurement weights), the box-constraint rows, and the gradient maps. The
 * arrival weight is the inverse of the steady-state Kalman covariance from
 * design::kalman on the same (Q, R) — so the unconstrained estimator agrees
 * with the steady-state Kalman filter. Fails (success = false) if the plant is
 * not discrete, has feedthrough (D ≠ 0), any weight inversion or the arrival
 * Riccati solve fails, or the Hessian is not positive definite.
 *
 * Q_process is the trust in the model (smaller = dynamics enforced harder);
 * R_measurement the trust in the sensor — the same covariances a Kalman filter
 * design uses, and a good design reuses exactly those numbers.
 *
 * @tparam NH  Estimation horizon (window covers NH+1 states / measurements)
 * @param sys                Discrete-time plant (Ts > 0, D = 0)
 * @param Q_process          Process-noise covariance (NX × NX, positive definite)
 * @param R_measurement      Measurement-noise covariance (NY × NY, positive definite)
 * @param constraints        State box limits (default: unbounded)
 * @param max_qp_iterations  Per-tick QP budget (0 → 10·(NV + NIB))
 * @return MHEArtifacts for damp::MHE
 */
template<size_t NH, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV_SS>
[[nodiscard]] constexpr MHEArtifacts<NX, NU, NY, NH, T> mhe(
    const StateSpace<NX, NU, NY, T, NW, NV_SS>& sys,
    const Matrix<NX, NX, T>&                    Q_process,
    const Matrix<NY, NY, T>&                    R_measurement,
    const MHEConstraints<NX, T>&                constraints = {},
    size_t                                      max_qp_iterations = 0
) {
    static_assert(NH >= 1, "The estimation window must cover at least two states.");

    constexpr size_t NV = NX * (NH + 1);
    constexpr size_t NIB = 2 * NV;

    MHEArtifacts<NX, NU, NY, NH, T> art{};
    art.constraints = constraints;
    art.max_qp_iterations = (max_qp_iterations == 0) ? 10 * (NV + NIB) : max_qp_iterations;

    if (!(sys.Ts > T{0})) {
        return art; // discrete plants only
    }
    for (size_t r = 0; r < NY; ++r) {
        for (size_t c = 0; c < NU; ++c) {
            if (damp::abs(sys.D(r, c)) > damp::default_tol<T>()) {
                return art; // feedthrough unsupported
            }
        }
    }

    // Weights: W_q = Q⁻¹, W_r = R⁻¹, W_a = P∞⁻¹ (steady-state Kalman covariance).
    const auto Wq_opt = mat::solve(Q_process, Matrix<NX, NX, T>::identity());
    const auto Wr_opt = mat::solve(R_measurement, Matrix<NY, NY, T>::identity());
    if (!Wq_opt || !Wr_opt) {
        return art;
    }
    const Matrix<NX, NX, T> Wq = Wq_opt.value();
    const Matrix<NY, NY, T> Wr = Wr_opt.value();

    StateSpace<NX, NU, NY, T, NX, NY> sys_noise{};
    sys_noise.A = sys.A;
    sys_noise.B = sys.B;
    sys_noise.C = sys.C;
    sys_noise.Ts = sys.Ts;
    const auto kf = kalman(sys_noise, Q_process, R_measurement);
    if (!kf.success) {
        return art;
    }
    const auto Wa_opt = mat::solve(kf.P, Matrix<NX, NX, T>::identity());
    if (!Wa_opt) {
        return art;
    }
    const Matrix<NX, NX, T> Wa = Wa_opt.value();

    // Constant Hessian (×2 for the ½zᵀHz convention): block tridiagonal.
    const Matrix<NX, NX, T> AtWqA = sys.A.transpose() * (Wq * sys.A);
    const Matrix<NX, NX, T> AtWq = sys.A.transpose() * Wq;
    const Matrix<NX, NX, T> CtWrC = sys.C.transpose() * (Wr * sys.C);
    for (size_t k = 0; k <= NH; ++k) {
        Matrix<NX, NX, T> diag = CtWrC;
        if (k == 0) {
            diag += Wa;
        }
        if (k < NH) {
            diag += AtWqA;
        }
        if (k > 0) {
            diag += Wq;
        }
        art.H.template block<NX, NX>(k * NX, k * NX) = T{2} * diag;
        if (k < NH) {
            art.H.template block<NX, NX>(k * NX, (k + 1) * NX) = T{-2} * AtWq;
            art.H.template block<NX, NX>((k + 1) * NX, k * NX) = T{-2} * AtWq.transpose();
        }
    }

    // Gradient maps for the runtime.
    art.M_a = T{2} * Wa;
    art.M_cy = T{2} * (sys.C.transpose() * Wr);
    art.M_qu = T{2} * (Wq * sys.B);
    art.M_aqu = T{2} * (sys.A.transpose() * (Wq * sys.B));

    // Box rows: x_k ≤ x_max and −x_k ≤ −x_min for every window state.
    for (size_t v = 0; v < NV; ++v) {
        const size_t state = v % NX;
        art.A_con(v, v) = T{1};
        art.A_con(NV + v, v) = T{-1};
        art.b_bound(v) = constraints.x_max(state);
        const T lo = constraints.x_min(state);
        art.b_bound(NV + v) = (lo <= -unbounded_bound<T>() / T{4}) ? unbounded_bound<T>() : -lo;
    }

    if (!mat::cholesky(art.H)) {
        return art;
    }

    art.success = true;
    return art;
}

} // namespace design

/**
 * @brief Runtime moving-horizon estimator (fixed per-tick iteration budget)
 *
 * Feed one measurement + previously applied input per tick; the estimator
 * slides its window, rebuilds the QP gradient from its ring buffers, solves
 * under the iteration budget, and returns the constrained estimate of the
 * current state. The first call primes the whole window with the first sample,
 * so there is no warm-up branch — early estimates simply assume the system was
 * at rest at the first measurement. Seed set_initial_state() beforehand when a
 * better prior is available.
 *
 * The smoothed window states are available via window_state() — index 0 is the
 * oldest, NH the current estimate.
 *
 * @see design::mhe, KalmanFilter for the unconstrained equivalent
 */
template<size_t NX, size_t NU, size_t NY, size_t NH, typename T = float, typename Solver = design::WarmStartActiveSetSolver<NX*(NH + 1), 2 * NX*(NH + 1), T>>
    requires std::is_floating_point_v<T>
class MHE {
public:
    using Artifacts = design::MHEArtifacts<NX, NU, NY, NH, T>;
    static constexpr size_t NV = Artifacts::NV;
    static constexpr size_t NIB = Artifacts::NIB;

    constexpr MHE() = default;
    /// Failed artifacts: valid_ = false; update holds last estimate and skips the QP.
    constexpr explicit MHE(const Artifacts& artifacts)
        : art(artifacts), valid_(artifacts.success) {}

    /// Seed the arrival-state prior (call before the first update)
    constexpr void set_initial_state(const ColVec<NX, T>& x0) { x_bar = x0; }

    /// true if constructed from a successful design::mhe result
    [[nodiscard]] constexpr bool valid() const { return valid_; }

    /**
     * @brief Slide the window and solve for the current state estimate
     *
     * @param y  Measurement taken now
     * @param u  Input applied over the previous interval (drove the plant to now)
     * @return Constrained estimate of the current state
     */
    constexpr ColVec<NX, T> update(const ColVec<NY, T>& y, const ColVec<NU, T>& u = ColVec<NU, T>{}) {
        if (!valid_) {
            return state();
        }
        if (!primed) {
            for (size_t k = 0; k <= NH; ++k) {
                ys[k] = y;
            }
            for (size_t k = 0; k < NH; ++k) {
                us[k] = u;
            }
            primed = true;
        } else {
            for (size_t k = 0; k < NH; ++k) {
                ys[k] = ys[k + 1];
            }
            ys[NH] = y;
            for (size_t k = 0; k + 1 < NH; ++k) {
                us[k] = us[k + 1];
            }
            us[NH - 1] = u;
            // The window start advanced by one: its prior mean is the previous
            // solve's (smoothed) estimate of that state.
            for (size_t i = 0; i < NX; ++i) {
                x_bar(i) = z(NX + i);
            }
        }

        // Gradient from the ring buffers (Hessian and bounds are constant).
        ColVec<NV, T> f{};
        for (size_t k = 0; k <= NH; ++k) {
            const auto gy = art.M_cy * ys[k];
            for (size_t i = 0; i < NX; ++i) {
                f(k * NX + i) -= gy(i);
            }
        }
        for (size_t k = 0; k < NH; ++k) {
            const auto gu_next = art.M_qu * us[k];
            const auto gu_curr = art.M_aqu * us[k];
            for (size_t i = 0; i < NX; ++i) {
                f(k * NX + i) += gu_curr(i);
                f((k + 1) * NX + i) -= gu_next(i);
            }
        }
        const auto ga = art.M_a * x_bar;
        for (size_t i = 0; i < NX; ++i) {
            f(i) -= ga(i);
        }

        const auto qp = qp_solver_(art.H, f, art.A_con, art.b_bound, art.max_qp_iterations);
        last_status_ = qp.status;
        last_iterations_ = qp.iterations;
        if (qp.status == design::QPStatus::Success || qp.status == design::QPStatus::MaxIterations) {
            z = qp.x;
        }

        ColVec<NX, T> estimate{};
        for (size_t i = 0; i < NX; ++i) {
            estimate(i) = z(NH * NX + i);
        }
        return estimate;
    }

    /// update() under the StateEstimator concept name: slide, solve, return the estimate.
    constexpr ColVec<NX, T> estimate(const ColVec<NY, T>& y, const ColVec<NU, T>& u = ColVec<NU, T>{}) {
        return update(y, u);
    }

    /// Current state estimate under the StateEstimator concept name.
    [[nodiscard]] constexpr ColVec<NX, T> state() const { return window_state(NH); }

    /// Smoothed window state k (0 = oldest, NH = current)
    [[nodiscard]] constexpr ColVec<NX, T> window_state(size_t k) const {
        ColVec<NX, T> x{};
        const size_t  kk = (k <= NH) ? k : NH;
        for (size_t i = 0; i < NX; ++i) {
            x(i) = z(kk * NX + i);
        }
        return x;
    }

    /// Current state estimate (last window state of the most recent solve)
    [[nodiscard]] constexpr ColVec<NX, T> state_estimate() const { return window_state(NH); }

    constexpr void reset() {
        primed = false;
        z = ColVec<NV, T>{};
        x_bar = ColVec<NX, T>{};
        last_status_ = design::QPStatus::Success;
        last_iterations_ = 0;
    }

    [[nodiscard]] constexpr design::QPStatus last_status() const { return last_status_; }
    [[nodiscard]] constexpr size_t           last_iterations() const { return last_iterations_; }
    [[nodiscard]] constexpr Solver&          qp_solver() { return qp_solver_; }

private:
    Artifacts                          art{};
    Solver                             qp_solver_{};
    damp::array<ColVec<NY, T>, NH + 1> ys{};
    damp::array<ColVec<NU, T>, NH>     us{};
    ColVec<NV, T>                      z{};
    ColVec<NX, T>                      x_bar{};
    design::QPStatus                   last_status_{design::QPStatus::Success};
    size_t                             last_iterations_{0};
    bool                               valid_{false};
    bool                               primed{false};
};

/// Deduce the runtime from its artifacts: MHE estimator{art};
template<size_t NX, size_t NU, size_t NY, size_t NH, typename T>
MHE(const design::MHEArtifacts<NX, NU, NY, NH, T>&) -> MHE<NX, NU, NY, NH, T>;

} // namespace damp
