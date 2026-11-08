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
 * Relative to a Kalman filter, MHE can enforce state constraints on the
 * estimate (e.g. non-negative level, pressure, or state of charge via
 * a lower bound of 0). Use MHE when those constraints matter; keep the Kalman filter
 * when they do not.
 *
 * The QP is fixed-size and allocation-free (Hessian and constraint rows are
 * constant; only the gradient changes with data), solved by the warm-started
 * active-set solver — at steady state the active constraint set barely changes
 * between ticks, so solves are a few iterations. Approximation vs. the exact
 * (full-information) estimator, standard practice: the arrival weight is held
 * at the steady-state a-priori @f$ P_\infty^{-1} @f$ rather than propagated.
 * The arrival mean is the matching a-priori prediction, advanced each tick by
 * one steady-state Kalman step over the sample leaving the window,
 * @f$ \bar x \leftarrow A\bigl(\bar x + M(y_0 - C\bar x)\bigr) + B u_0 @f$ with the
 * innovation gain M (the filtering prior of Rao et al.). Using the previous solve's smoothed
 * estimate instead would count the measurements still inside the window twice.
 *
 * The noise enters through the plant's G and H, as in design::kalman: the window weights
 * @f$ G Q G^\top @f$ and @f$ H R H^\top @f$.
 *
 * Feedthrough is not supported (D = 0), matching the MPC. Bounds use the
 * `numeric_limits` sentinels (never infinity; -ffinite-math-only).
 *
 * @note Compare with MATLAB®'s nlmhe / nlmheMultistage restricted to the linear
 *       case, with mhestate's bounds as MHEConstraints.
 * @see estimation/kalman.hpp for the unconstrained recursive equivalent,
 *      design/qp.hpp for the solver, controllers/mpc.hpp for the control-side
 *      mirror of this formulation
 * @see "Constrained State Estimation for Nonlinear Discrete-Time Systems:
 *      Stability and Moving Horizon Approximations" (Rao, Rawlings & Mayne, 2003)
 * @see "Model Predictive Control: Theory, Computation, and Design"
 *      (Rawlings, Mayne & Diehl, 2017), ch. 4
 *
 * @code
 * // Level estimator that can never report a negative level:
 * constexpr design::MHEConstraints<1> limits{.x_lim = {.lower = 0.0, .upper = std::numeric_limits<double>::max()}};
 * constexpr auto res = design::mhe<8>(tank, {.process = Q, .measurement = R}, limits);
 * static constinit MHE estimator{res.as<float>()}; // a failed design does not compile
 * // each tick: x_hat = estimator.step(y, u_previous);
 * @endcode
 */


#include "damp/backend.hpp"
#include "damp/design/qp.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/block.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "damp/matrix/solve.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/toolbox/bounds.hpp"

namespace damp {

namespace design {

/**
 * @struct MHEConstraints
 * @brief Box constraints on the estimated states (applied to every window state) and the QP budget
 *
 * Defaults are unbounded (numeric_limits sentinels) — with no finite bounds the
 * estimator is an unconstrained least-squares smoother matching the Kalman
 * filter at the window endpoint.
 */
template<size_t NX, typename T = double>
struct MHEConstraints {
    Bounds<NX, T> x_lim{};              ///< State bounds (default unbounded); no NaN end, lower ≤ upper
    size_t        max_qp_iterations{0}; ///< Per-tick QP budget (0 → 10·(NV + NIB))

    template<typename U>
    [[nodiscard]] constexpr MHEConstraints<NX, U> as() const {
        return MHEConstraints<NX, U>{x_lim.template as<U>(), max_qp_iterations};
    }
};

/**
 * @struct MHEResult
 * @brief Window QP data produced by mhe(), consumed by damp::MHE
 *
 * The Hessian and constraint rows are constant; the runtime rebuilds only the
 * gradient from its measurement/input ring buffers each tick using the stored
 * weight maps. Use .as<float>() for embedded deployment (bound sentinels
 * convert safely); it refuses (singular) a Hessian that does not factor in float.
 */
template<size_t NX, size_t NU, size_t NY, size_t NH, typename T = double>
struct [[nodiscard]] MHEResult {
    static constexpr size_t NV = NX * (NH + 1); ///< Decision variables (window states)
    static constexpr size_t NIB = 2 * NV;       ///< Box-constraint rows

    Matrix<NV, NV, T>     H{};                                 ///< Window QP Hessian (positive definite)
    Matrix<NIB, NV, T>    A_con{};                             ///< Box rows: [I; −I]
    ColVec<NIB, T>        b_bound{};                           ///< [x_lim.upper…; −x_lim.lower…] (constant; sentinel rows disabled)
    Matrix<NX, NX, T>     M_a{};                               ///< 2·P∞⁻¹ (arrival gradient map)
    Matrix<NX, NY, T>     M_cy{};                              ///< 2·CᵀR⁻¹ (measurement gradient map)
    Matrix<NX, NU, T>     M_qu{};                              ///< 2·Q⁻¹B (process gradient map, next state)
    Matrix<NX, NU, T>     M_aqu{};                             ///< 2·AᵀQ⁻¹B (process gradient map, current state)
    Matrix<NX, NX, T>     Ap_x{};                              ///< A(I − M C), M the innovation gain (arrival-prior update, prior term)
    Matrix<NX, NY, T>     Ap_y{};                              ///< A·M, the Kalman predictor gain (arrival-prior update, measurement term)
    Matrix<NX, NU, T>     Ap_u{};                              ///< B (arrival-prior update, input term)
    MHEConstraints<NX, T> constraints{};                       ///< Original limits
    size_t                max_qp_iterations{};                 ///< Per-tick QP iteration budget
    bool                  success{false};                      ///< status == DesignStatus::ok
    DesignStatus          status{DesignStatus::not_converged}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr MHEResult<NX, NU, NY, NH, U> as() const {
        MHEResult<NX, NU, NY, NH, U> out{
            H.template as<U>(),
            A_con.template as<U>(),
            detail::convert_bounds<U>(b_bound),
            M_a.template as<U>(),
            M_cy.template as<U>(),
            M_qu.template as<U>(),
            M_aqu.template as<U>(),
            Ap_x.template as<U>(),
            Ap_y.template as<U>(),
            Ap_u.template as<U>(),
            constraints.template as<U>(),
            max_qp_iterations,
            success,
            status,
        };
        //! The design certified chol(H) in its own scalar. Noise on fewer channels than states
        //! weights the noise-free directions by 1/√ε_T, past 1/ε_U in condition for a float U,
        //! and every tick would fail with NotPositiveDefinite.
        if (out.success && !mat::cholesky(out.H)) {
            out.success = false;
            out.status = DesignStatus::singular;
        }
        return out;
    }
};

namespace detail {

/**
 * @brief The window's covariance W + δ I, δ = √ε·max diag(W): invertible when W is singular
 *        (noise on fewer channels than states, deterministic delay registers).
 */
template<size_t N, typename T>
[[nodiscard]] constexpr Matrix<N, N, T> window_covariance(const Matrix<N, N, T>& W) {
    T w_max{0};
    for (size_t i = 0; i < N; ++i) {
        w_max = damp::max(w_max, damp::abs(W(i, i)));
    }
    return W + (Matrix<N, N, T>::identity() * (damp::sqrt(eps<T>()) * w_max));
}

/** @brief State limits the window accepts: no NaN end (∓max() for none) and no lower end above its upper one. */
template<size_t N, typename T>
[[nodiscard]] constexpr bool mhe_limits_valid(const Bounds<N, T>& b) {
    for (size_t i = 0; i < N; ++i) {
        if (!damp::isfinite(b.lower(i)) || !damp::isfinite(b.upper(i)) || !(b.lower(i) <= b.upper(i))) {
            return false;
        }
    }
    return true;
}

} // namespace detail

/**
 * @brief Synthesize a constrained moving-horizon estimator
 *
 * Builds the constant window Hessian (block tridiagonal from the dynamics and
 * measurement weights), the box-constraint rows, and the gradient maps. The noise
 * enters through the plant's G and H, as in design::kalman: the window weights are the
 * inverses of @f$ G Q G^\top + \delta_q I @f$ and @f$ H R H^\top + \delta_r I @f$, each
 * @f$ \delta = \sqrt{\varepsilon}\, \max_i W_{ii} @f$ of its own matrix. G Q Gᵀ is singular
 * whenever the noise drives fewer channels than states (a G with fewer columns than rows, or
 * the delay registers of discretize(plant, Ts, InputDelay{θ}), which get no noise); δ keeps
 * the window solvable. The arrival weight is the inverse of the steady-state Kalman
 * covariance on the same two covariances, so the unconstrained estimator agrees with that
 * Kalman filter.
 *
 * Status: invalid_sample_time if the plant is not discrete; bad_weights on a non-finite model
 * or noise entry, a feedthrough D ≠ 0 or a zero covariance; invalid_parameter on a state limit
 * with a NaN end or lower > upper; design::kalman's status if the arrival filter fails;
 * singular if P∞ or the window Hessian is not positive definite.
 *
 * Q is the trust in the model (smaller = dynamics enforced harder); R the trust in the
 * sensor — the same covariances a Kalman filter design uses, and a good design reuses exactly
 * those numbers.
 *
 * @tparam NH  Estimation horizon (window covers NH+1 states / measurements)
 * @param sys          Discrete-time plant (Ts > 0, D = 0); G and H carry the noise
 * @param noise        Process-noise covariance Q (NW × NW) and measurement-noise covariance R
 *                     (NV × NV), positive semidefinite
 * @param constraints  State box limits (default: unbounded) and the per-tick QP budget
 * @return MHEResult for damp::MHE
 */
template<size_t NH, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV_SS>
[[nodiscard]] constexpr MHEResult<NX, NU, NY, NH, T> mhe(
    const StateSpace<NX, NU, NY, T, NW, NV_SS>& sys,
    const KalmanNoise<NW, NV_SS, T>&            noise,
    const MHEConstraints<NX, T>&                constraints = {}
) {
    static_assert(NH >= 1, "The estimation window must cover at least two states.");

    constexpr size_t NV = NX * (NH + 1);
    constexpr size_t NIB = 2 * NV;

    MHEResult<NX, NU, NY, NH, T> art{};
    art.constraints = constraints;
    art.max_qp_iterations = (constraints.max_qp_iterations == 0) ? 10 * (NV + NIB) : constraints.max_qp_iterations;

    if (!damp::finite_positive(sys.Ts)) {
        return ::damp::detail::set_status(art, DesignStatus::invalid_sample_time); // discrete plants only
    }
    if (!detail::finite_filter_model(sys, noise.process, noise.measurement)) {
        return ::damp::detail::set_status(art, DesignStatus::bad_weights);
    }
    if (!detail::mhe_limits_valid(constraints.x_lim)) {
        return ::damp::detail::set_status(art, DesignStatus::invalid_parameter); // as design::mpc's limits
    }
    if (!mat::all_of(sys.D, [](T v) { return v == T{0}; })) {
        return ::damp::detail::set_status(art, DesignStatus::bad_weights); // feedthrough unsupported
    }

    // Weights: W_q = Q⁻¹, W_r = R⁻¹, W_a = P∞⁻¹ (steady-state Kalman covariance), on the noise
    // as it reaches the states and the outputs.
    const Matrix<NX, NX, T> Q_process = detail::window_covariance(Matrix<NX, NX, T>(sys.G * noise.process * sys.G.t()));
    const Matrix<NY, NY, T> R_measurement = detail::window_covariance(Matrix<NY, NY, T>(sys.H * noise.measurement * sys.H.t()));
    const auto              Wq_opt = mat::solve(Q_process, Matrix<NX, NX, T>::identity());
    const auto              Wr_opt = mat::solve(R_measurement, Matrix<NY, NY, T>::identity());
    if (!Wq_opt || !Wr_opt) {
        return ::damp::detail::set_status(art, DesignStatus::bad_weights); // a zero or indefinite covariance
    }
    const Matrix<NX, NX, T> Wq = Wq_opt.value();
    const Matrix<NY, NY, T> Wr = Wr_opt.value();

    const StateSpace<NX, NU, NY, T, NX, NY> sys_noise{.A = sys.A, .B = sys.B, .C = sys.C, .Ts = sys.Ts};
    const auto                              kf = kalman(sys_noise, {.process = Q_process, .measurement = R_measurement});
    if (!kf.success) {
        return ::damp::detail::set_status(art, kf.status);
    }
    const auto Wa_opt = mat::solve(kf.P, Matrix<NX, NX, T>::identity());
    if (!Wa_opt) {
        return ::damp::detail::set_status(art, DesignStatus::singular);
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

    // Arrival prior: one steady-state Kalman step over the sample leaving the window.
    art.Ap_x = sys.A * (Matrix<NX, NX, T>::identity() - (kf.M * sys.C));
    art.Ap_y = sys.A * kf.M;
    art.Ap_u = sys.B;

    // Box rows: x_k ≤ x_lim.upper and −x_k ≤ −x_lim.lower for every window state.
    for (size_t v = 0; v < NV; ++v) {
        const size_t state = v % NX;
        art.A_con(v, v) = T{1};
        art.A_con(NV + v, v) = T{-1};
        art.b_bound(v) = constraints.x_lim.upper(state);
        const T lo = constraints.x_lim.lower(state);
        art.b_bound(NV + v) = (lo <= -unbounded_bound<T>() / T{4}) ? unbounded_bound<T>() : -lo;
    }

    if (!mat::cholesky(art.H)) {
        return ::damp::detail::set_status(art, DesignStatus::singular);
    }
    return ::damp::detail::set_status(art, DesignStatus::ok);
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
    requires damp::is_floating_point_v<T>
class MHE {
public:
    using Result = design::MHEResult<NX, NU, NY, NH, T>;
    static constexpr size_t NV = Result::NV;
    static constexpr size_t NIB = Result::NIB;

    constexpr MHE() = default;
    /**
     * @brief From design::mhe. A failed design does not compile when constant-evaluated; at run
     *        time valid() is false and step() holds the last estimate, skipping the QP.
     */
    constexpr MHE(const Result& result) // NOLINT(google-explicit-constructor): assign-from-design
        : art(result), valid_(result.success) {
        if (!result.success && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(result.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    /** @brief Seed the arrival-state prior (call before the first step) */
    constexpr void set_initial_state(const ColVec<NX, T>& x0) { x_bar = x0; }

    /** @brief true if constructed from a successful design::mhe result */
    [[nodiscard]] constexpr bool valid() const { return valid_; }

    /**
     * @brief One tick: slide the window and solve for the current state estimate
     *
     * @param y  Measurement taken now
     * @param u  Input applied over the previous interval (drove the plant to now)
     * @return Constrained estimate of the current state
     */
    constexpr ColVec<NX, T> step(const ColVec<NY, T>& y, const ColVec<NU, T>& u) {
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
            // The window start advances by one: carry its a-priori mean through
            // one steady-state Kalman step over the departing sample (y₀, u₀).
            x_bar = (art.Ap_x * x_bar) + (art.Ap_y * ys[0]) + (art.Ap_u * us[0]);
            for (size_t k = 0; k < NH; ++k) {
                ys[k] = ys[k + 1];
            }
            ys[NH] = y;
            for (size_t k = 0; k + 1 < NH; ++k) {
                us[k] = us[k + 1];
            }
            us[NH - 1] = u;
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
        if (qp.status == design::QPStatus::Success) {
            z = qp.x;
        } else if (qp.status == design::QPStatus::MaxIterations) {
            // The dual iterate is not primal feasible: project it onto the box, this QP's whole
            // feasible set.
            for (size_t v = 0; v < NV; ++v) {
                z(v) = damp::clamp(qp.x(v), art.constraints.x_lim.lower(v % NX), art.constraints.x_lim.upper(v % NX));
            }
        }

        ColVec<NX, T> estimate{};
        for (size_t i = 0; i < NX; ++i) {
            estimate(i) = z(NH * NX + i);
        }
        return estimate;
    }

    /** @brief Current state estimate under the StateEstimator concept name. */
    [[nodiscard]] constexpr ColVec<NX, T> state() const { return window_state(NH); }

    /** @brief Smoothed window state k (0 = oldest, NH = current) */
    [[nodiscard]] constexpr ColVec<NX, T> window_state(size_t k) const {
        ColVec<NX, T> x{};
        const size_t  kk = (k <= NH) ? k : NH;
        for (size_t i = 0; i < NX; ++i) {
            x(i) = z(kk * NX + i);
        }
        return x;
    }

    /** @brief Current state estimate (last window state of the most recent solve) */
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
    Result                             art{};
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

/** @brief Deduce the runtime from its design: MHE estimator{res}; */
template<size_t NX, size_t NU, size_t NY, size_t NH, typename T>
MHE(const design::MHEResult<NX, NU, NY, NH, T>&) -> MHE<NX, NU, NY, NH, T>;

} // namespace damp
