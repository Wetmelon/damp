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
 *     + \sum_{k=0}^{NH} \lVert y_k - C x_k \rVert^2_{R^{-1}}
 *     + \rho_\varepsilon \varepsilon^2,
 * @f]
 * subject to the limits a StateLimits sets on the states it lists, at every window state. The
 * dynamics residuals are the process noises, the measurement residuals the sensor noises, and
 * the arrival term anchors the window start to the prior mean @f$ \bar x @f$ with the
 * steady-state Kalman covariance @f$ P_\infty @f$ as its weight — on a linear, unconstrained
 * problem the window endpoint therefore reproduces the steady-state Kalman filter estimate.
 *
 * Relative to a Kalman filter, MHE can enforce state constraints on the
 * estimate (e.g. non-negative level, pressure, or state of charge via
 * a lower bound of 0). Use MHE when those constraints matter; keep the Kalman filter
 * when they do not.
 *
 * Soft limits, as design::mpc's. A StateLimits is soft by default: one slack
 * @f$ \varepsilon \ge 0 @f$ relaxes every soft row, @f$ x_i / s_i \le \bar x_i / s_i + V_i \varepsilon @f$
 * with @f$ s_i @f$ the state's range (1 where an end is unbounded) and @f$ V_i @f$ the row's
 * ECR coefficient (StateLimits::lower_ecr, upper_ecr), at the cost
 * @f$ \rho_\varepsilon \varepsilon^2 @f$, @f$ \rho_\varepsilon = 10^5 @f$ (MPCWeights::ecr_weight's
 * default). The residuals of J are in standard deviations of their noise, so a violation of
 * @f$ s_i V_i / 100 @f$ costs about as much as one reading 3σ off: a limit that holds against
 * noise but gives way to a run of readings that disagree with it (a sensor offset, a limit set
 * wrong). A coefficient
 * of 0 makes a row hard, and the estimate then never leaves it; a larger one lets it go more
 * readily (a violation v costs @f$ \rho_\varepsilon (v / (s_i V_i))^2 @f$).
 *
 * The QP is fixed-size and allocation-free (Hessian and constraint rows are
 * constant; only the gradient changes with data), solved by the warm-started
 * active-set solver — at steady state the active constraint set barely changes
 * between ticks, so solves are a few iterations. It carries rows for the listed states alone,
 * and with no StateLimits none at all: the solve is the unconstrained minimum. Approximation vs.
 * the exact (full-information) estimator, standard practice: the arrival weight is held
 * at the steady-state a-priori @f$ P_\infty^{-1} @f$ rather than propagated.
 * The arrival mean is the matching a-priori prediction, advanced each tick by
 * one steady-state Kalman step over the sample leaving the window,
 * @f$ \bar x \leftarrow A\bigl(\bar x + M(y_0 - C\bar x)\bigr) + B u_0 @f$ with the
 * innovation gain M (the filtering prior of Rao et al.). Using the previous solve's smoothed
 * estimate instead would count the measurements still inside the window twice.
 *
 * The noise enters through the plant's G and H, as in design::kalman: the window weights
 * @f$ G Q G^\top @f$ and @f$ H R H^\top @f$. A discrete plant takes Q per step,
 * design::mhe<NH>(sys, noise, …); a continuous one, design::mhe<NH>(plant, Ts, noise, …), takes
 * it as a density, discretized exactly over the period (Van Loan), as design::kalman's sampled
 * form does. R is the variance of one reading on either.
 *
 * Feedthrough is not supported (D = 0), matching the MPC. Bounds use the
 * `numeric_limits` sentinels (never infinity; -ffinite-math-only).
 *
 * @note Compare with MATLAB®'s nlmhe / nlmheMultistage restricted to the linear
 *       case, with mhestate's bounds as a hard StateLimits.
 * @see estimation/kalman.hpp for the unconstrained recursive equivalent,
 *      design/qp.hpp for the solver, controllers/mpc.hpp for the control-side
 *      mirror of this formulation
 * @see "Constrained State Estimation for Nonlinear Discrete-Time Systems:
 *      Stability and Moving Horizon Approximations" (Rao, Rawlings & Mayne, 2003)
 * @see "Model Predictive Control: Theory, Computation, and Design"
 *      (Rawlings, Mayne & Diehl, 2017), ch. 4
 * @see "Predictive Control with Constraints" (Maciejowski, 2002), §3.4: the soft constraints
 *
 * @code
 * // Level estimator that can never report a negative level (a hard lower end), its QP capped at 20 iterations a tick:
 * constexpr design::StateLimits level{.states = {0}, .lower = {0.0}, .lower_ecr = {0.0}};
 * constexpr auto res = design::mhe<8>(tank, {.process = Q, .measurement = R}, level, design::QPIterations{20});
 * static constinit MHE estimator{res.as<float>()}; // a failed design does not compile
 * // each tick: x_hat = estimator.step(y, u_previous);
 * // From the continuous model at Ts = 0.1, the process noise a density:
 * constexpr auto sampled = design::mhe<8>(tank_c, 0.1, {.process = Qc, .measurement = R}, level);
 * @endcode
 */


#include "damp/backend.hpp"
#include "damp/controllers/detail/law_detail.hpp"
#include "damp/controllers/law.hpp" // IWYU pragma: export
#include "damp/design/input_limits.hpp"
#include "damp/design/options.hpp"
#include "damp/design/qp.hpp"           // IWYU pragma: export
#include "damp/design/state_limits.hpp" // IWYU pragma: export
#include "damp/design/status.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/block.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "damp/matrix/solve.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/toolbox/bounds.hpp" // IWYU pragma: export

namespace damp {

namespace design {

/**
 * @brief The limits an MHE takes at run time: the bounds on the states its design's StateLimits
 *        lists, in physical units and in the order of StateLimits::states. One left out is
 *        unbounded, and with no StateLimits state is empty.
 *
 * The softness (ECR), the weights and the QP budget stay as designed. `auto l = mhe.limits();
 * l.state.lower(0) = …; mhe.set_limits(l);` moves one, as on an MPC (MPCLimits::state). A design's
 * result holds the limits its runtime starts with in the same struct (MHEResult::limits).
 *
 * @tparam NXC States with bounds (the design's StateLimits; 0: none, state empty)
 */
template<size_t NXC = 0, typename T = double>
struct MHELimits {
    Bounds<NXC, T> state{}; ///< Bounds on the states the design's StateLimits lists

    /** @brief The limits in another scalar (an unbounded end stays one). */
    template<typename U>
    [[nodiscard]] constexpr MHELimits<NXC, U> as() const {
        return {state.template as<U>()};
    }
};

/**
 * @struct MHEResult
 * @brief Window QP data produced by mhe(), consumed by damp::MHE
 *
 * The Hessian and constraint rows are constant; the runtime rebuilds only the
 * gradient from its measurement/input ring buffers each tick using the stored
 * weight maps, and set_limits rewrites the bound vector. The QP is over the window states and,
 * with a StateLimits, the slack ε of its soft rows. Use .as<float>() for embedded deployment
 * (bound sentinels convert safely); it refuses (singular) a Hessian that does not factor in float.
 *
 * @tparam NH  Estimation horizon (the window holds NH + 1 states)
 * @tparam NXC States with bounds (the design's StateLimits; 0: none, no rows)
 */
template<size_t NX, size_t NU, size_t NY, size_t NH, typename T = double, size_t NXC = 0>
struct [[nodiscard]] MHEResult {
    /** @brief The design in another scalar: as<U>()'s type. */
    template<typename U>
    using as_t = MHEResult<NX, NU, NY, NH, U, NXC>;

    using scalar = T;                                 ///< The scalar of its matrices
    static constexpr size_t nx = NX;                  ///< States
    static constexpr size_t nu = NU;                  ///< Inputs
    static constexpr size_t ny = NY;                  ///< Outputs
    static constexpr size_t nh = NH;                  ///< Estimation horizon
    static constexpr size_t nxc = NXC;                ///< Bounded states
    static constexpr size_t NV = NX * (NH + 1);       ///< Window states
    static constexpr size_t NS = (NXC > 0) ? 1 : 0;   ///< Slack ε of the soft rows: one with a StateLimits, none without
    static constexpr size_t NZ = NV + NS;             ///< QP decision variables: the window states, then ε
    static constexpr size_t NIB = 2 * NXC * (NH + 1); ///< Box rows: the listed states' upper rows, then their lower ones, at every window state
    static constexpr size_t NI = NIB + NS;            ///< QP constraint rows: the box rows, then ε ≥ 0

    Matrix<NZ, NZ, T>        H{};                                           ///< Window QP Hessian (positive definite), 2ρ_ε on ε
    Matrix<NI, NZ, T>        A_con{};                                       ///< Box rows on the scaled listed states, −V on ε, then −ε ≤ 0
    ColVec<NI, T>            b_bound{};                                     ///< The scaled bounds [upper…; −lower…; 0] (sentinel rows disabled)
    Matrix<NX, NX, T>        M_a{};                                         ///< 2·P∞⁻¹ (arrival gradient map)
    Matrix<NX, NY, T>        M_cy{};                                        ///< 2·CᵀR⁻¹ (measurement gradient map)
    Matrix<NX, NU, T>        M_qu{};                                        ///< 2·Q⁻¹B (process gradient map, next state)
    Matrix<NX, NU, T>        M_aqu{};                                       ///< 2·AᵀQ⁻¹B (process gradient map, current state)
    Matrix<NX, NX, T>        Ap_x{};                                        ///< A(I − M C), M the innovation gain (arrival-prior update, prior term)
    Matrix<NX, NY, T>        Ap_y{};                                        ///< A·M, the Kalman predictor gain (arrival-prior update, measurement term)
    Matrix<NX, NU, T>        Ap_u{};                                        ///< B (arrival-prior update, input term)
    MHELimits<NXC, T>        limits{};                                      ///< The limits in force: the design's, or the last set_limits
    damp::array<size_t, NXC> states{};                                      ///< The bounded states (StateLimits::states)
    ColVec<NXC, T>           lower_ecr = Matrix<NXC, 1, T>::constant(T{1}); ///< Softness of the limits.state.lower rows (0: hard)
    ColVec<NXC, T>           upper_ecr = Matrix<NXC, 1, T>::constant(T{1}); ///< Softness of the limits.state.upper rows (0: hard)
    ColVec<NXC, T>           scale_x{};                                     ///< Scale factors of the bounded states (their range, else 1)
    size_t                   max_qp_iterations{};                           ///< Per-tick QP iteration budget
    bool                     success{false};                                ///< status == DesignStatus::ok
    DesignStatus             status{DesignStatus::not_converged};           ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The design in another scalar; refused (singular) when H does not factor in it. */
    template<typename U>
    [[nodiscard]] constexpr as_t<U> as() const {
        as_t<U> out{
            .H = H.template as<U>(),
            .A_con = A_con.template as<U>(),
            .b_bound = detail::convert_bounds<U>(b_bound),
            .M_a = M_a.template as<U>(),
            .M_cy = M_cy.template as<U>(),
            .M_qu = M_qu.template as<U>(),
            .M_aqu = M_aqu.template as<U>(),
            .Ap_x = Ap_x.template as<U>(),
            .Ap_y = Ap_y.template as<U>(),
            .Ap_u = Ap_u.template as<U>(),
            .limits = limits.template as<U>(),
            .states = states,
            .lower_ecr = ColVec<NXC, U>{lower_ecr.template as<U>()},
            .upper_ecr = ColVec<NXC, U>{upper_ecr.template as<U>()},
            .scale_x = ColVec<NXC, U>{scale_x.template as<U>()},
            .max_qp_iterations = max_qp_iterations,
            .success = success,
            .status = status,
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

} // namespace design

namespace detail {

/**
 * @brief Write the box rows of the window QP's bound vector from @p lim with the design's scales:
 *        the listed states' upper rows, then their lower ones, at every window state. The design
 *        writes the first, set_limits each one after.
 */
template<size_t NX, size_t NU, size_t NY, size_t NH, typename T, size_t NXC>
constexpr void mhe_box_bounds(design::MHEResult<NX, NU, NY, NH, T, NXC>& art, const design::MHELimits<NXC, T>& lim) {
    constexpr size_t lower_rows = NXC * (NH + 1);
    for (size_t k = 0; k <= NH; ++k) {
        for (size_t j = 0; j < NXC; ++j) {
            const T s = art.scale_x(j);
            art.b_bound((k * NXC) + j) = design::detail::scaled_bound(lim.state.upper(j), s);
            art.b_bound(lower_rows + (k * NXC) + j) = design::detail::scaled_bound(-lim.state.lower(j), s);
        }
    }
}

} // namespace detail

namespace design {

namespace detail {

using ::damp::detail::mhe_box_bounds;

/** @brief The options design::mhe takes after the noise. */
using mhe_options = option_set<kind::StateLimits, kind::QPIterations>;

/** @brief The slack penalty ρ_ε of an MHE's soft rows: 10⁵, MPCWeights::ecr_weight's default. */
template<typename T>
inline constexpr T mhe_ecr_weight = static_cast<T>(1.0e5);

/** @brief The listed states of a design::mhe call, their limits and the softness of their rows: none without a StateLimits. */
template<size_t NXC, typename T>
struct MHEStateSet {
    MHELimits<NXC, T>        limits{};                                      ///< Their bounds (StateLimits::lower, upper)
    damp::array<size_t, NXC> states{};                                      ///< The bounded states (StateLimits::states)
    ColVec<NXC, T>           lower_ecr = Matrix<NXC, 1, T>::constant(T{1}); ///< Softness of the lower rows
    ColVec<NXC, T>           upper_ecr = Matrix<NXC, 1, T>::constant(T{1}); ///< Softness of the upper rows
};

/** @brief The StateLimits among @p opts in the scalar @p T, as the design holds it: none without one. */
template<typename T, typename... Opts>
[[nodiscard]] constexpr MHEStateSet<listed_by<kind::StateLimits, Opts...>, T> mhe_state_set(const Opts&... opts) {
    constexpr size_t    NXC = listed_by<kind::StateLimits, Opts...>;
    MHEStateSet<NXC, T> out{};
    if constexpr (NXC > 0) {
        const auto limits = the<kind::StateLimits>(opts...);
        copy_listed(out.states, limits.states);
        out.limits.state = box_of<T>(limits.lower, limits.upper);
        out.lower_ecr = ColVec<NXC, T>{limits.lower_ecr.template as<T>()};
        out.upper_ecr = ColVec<NXC, T>{limits.upper_ecr.template as<T>()};
    } else {
        (static_cast<void>(opts), ...);
    }
    return out;
}

/** @brief The QP budget among @p opts: QPIterations::budget, else 10·(NZ + NI), the window QP's variables and rows. */
template<typename Result, typename... Opts>
[[nodiscard]] constexpr size_t mhe_budget_of(const Opts&... opts) {
    return option<kind::QPIterations>(QPIterations{10 * (Result::NZ + Result::NI)}, opts...).budget;
}

/**
 * @brief The window's covariance W + δ I, δ = √ε·max diag(W): invertible when W is singular
 *        (noise on fewer channels than states, deterministic delay states).
 */
template<size_t N, typename T>
[[nodiscard]] constexpr Matrix<N, N, T> window_covariance(const Matrix<N, N, T>& W) {
    T w_max{0};
    for (size_t i = 0; i < N; ++i) {
        w_max = damp::max(w_max, damp::abs(W(i, i)));
    }
    return W + (Matrix<N, N, T>::identity() * (damp::sqrt(eps<T>()) * w_max));
}

/** @brief The listed states are the model's, each once, inside a box (invalid_parameter), and their softness finite and ≥ 0 (bad_weights); else ok. */
template<size_t NX, size_t NXC, typename T>
[[nodiscard]] constexpr DesignStatus mhe_limits_status(const MHEStateSet<NXC, T>& set) {
    if (!indices_valid<NX>(set.states) || !set.limits.state.valid()) {
        return DesignStatus::invalid_parameter; // an index past the model's or listed twice, or a NaN, infinite or inverted end
    }
    const auto softness = [](T v) { return damp::isfinite(v) && v >= T{0}; };
    if (!mat::all_of(set.lower_ecr, softness) || !mat::all_of(set.upper_ecr, softness)) {
        return DesignStatus::bad_weights; // as design::mpc's ECR
    }
    return DesignStatus::ok;
}

/** @brief The constant window Hessian (×2 for the ½zᵀHz convention): block tridiagonal in the window states, 2ρ_ε on the slack. */
template<size_t NX, size_t NU, size_t NY, size_t NH, typename T, size_t NXC>
constexpr void mhe_hessian(
    MHEResult<NX, NU, NY, NH, T, NXC>& art,
    const Matrix<NX, NX, T>&           A,
    const Matrix<NY, NX, T>&           C,
    const Matrix<NX, NX, T>&           Wq,
    const Matrix<NY, NY, T>&           Wr,
    const Matrix<NX, NX, T>&           Wa
) {
    const Matrix<NX, NX, T> AtWqA = A.transpose() * (Wq * A);
    const Matrix<NX, NX, T> AtWq = A.transpose() * Wq;
    const Matrix<NX, NX, T> CtWrC = C.transpose() * (Wr * C);
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
    if constexpr (NXC > 0) {
        constexpr size_t slack = MHEResult<NX, NU, NY, NH, T, NXC>::NV;
        art.H(slack, slack) = T{2} * mhe_ecr_weight<T>;
    }
}

/**
 * @brief The rows on the listed states, scaled by their range, at every window state: x_i/s_i ≤ upper/s_i + V⁺ε,
 *        −x_i/s_i ≤ −lower/s_i + V⁻ε, then −ε ≤ 0. None without a StateLimits.
 */
template<size_t NX, size_t NU, size_t NY, size_t NH, typename T, size_t NXC>
constexpr void mhe_state_rows(MHEResult<NX, NU, NY, NH, T, NXC>& art) {
    if constexpr (NXC > 0) {
        using Art = MHEResult<NX, NU, NY, NH, T, NXC>;
        constexpr size_t lower_rows = NXC * (NH + 1);
        for (size_t j = 0; j < NXC; ++j) {
            art.scale_x(j) = scale_from_range(art.limits.state.lower(j), art.limits.state.upper(j));
        }
        for (size_t k = 0; k <= NH; ++k) {
            for (size_t j = 0; j < NXC; ++j) {
                const size_t col = (k * NX) + art.states[j];
                const size_t hi = (k * NXC) + j;
                art.A_con(hi, col) = T{1} / art.scale_x(j);
                art.A_con(lower_rows + hi, col) = T{-1} / art.scale_x(j);
                art.A_con(hi, Art::NV) = -art.upper_ecr(j);
                art.A_con(lower_rows + hi, Art::NV) = -art.lower_ecr(j);
            }
        }
        art.A_con(Art::NIB, Art::NV) = T{-1};
        art.b_bound(Art::NIB) = T{0};
        mhe_box_bounds(art, art.limits);
    }
}

/**
 * @brief design::mhe on the discrete model @p sys, its noise inputs G and H given (noise_model):
 *        the window QP, its weights, its rows on the listed states and the arrival filter.
 */
template<size_t NH, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NVN, size_t NXC>
[[nodiscard]] constexpr MHEResult<NX, NU, NY, NH, T, NXC> mhe_design(
    const StateSpace<NX, NU, NY, T, NW, NVN>& sys,
    const KalmanNoise<NW, NVN, T>&            noise,
    const MHEStateSet<NXC, T>&                set,
    size_t                                    budget
) {
    static_assert(NH >= 1, "The estimation window must cover at least two states.");
    MHEResult<NX, NU, NY, NH, T, NXC> art{};
    art.limits = set.limits;
    art.states = set.states;
    art.lower_ecr = set.lower_ecr;
    art.upper_ecr = set.upper_ecr;
    art.max_qp_iterations = budget;

    if (!damp::finite_positive(sys.Ts)) {
        return set_status(art, DesignStatus::invalid_sample_time); // a discrete plant here; a continuous one takes mhe(plant, Ts, …)
    }
    if (!finite_filter_model(sys, noise.process, noise.measurement)) {
        return set_status(art, DesignStatus::bad_weights);
    }
    if (const auto st = mhe_limits_status<NX>(set); st != DesignStatus::ok) {
        return set_status(art, st);
    }
    if (art.max_qp_iterations == 0) {
        return set_status(art, DesignStatus::invalid_parameter); // QPIterations{0}: not one active-set step, as design::mpc
    }
    if (!mat::all_of(sys.D, [](T v) { return v == T{0}; })) {
        return set_status(art, DesignStatus::bad_weights); // feedthrough unsupported
    }

    // Weights: W_q = Q⁻¹, W_r = R⁻¹, W_a = P∞⁻¹ (steady-state Kalman covariance), on the noise
    // as it reaches the states and the outputs.
    const Matrix<NX, NX, T> Q_process = window_covariance(Matrix<NX, NX, T>(sys.G * noise.process * sys.G.t()));
    const Matrix<NY, NY, T> R_measurement = window_covariance(Matrix<NY, NY, T>(sys.H * noise.measurement * sys.H.t()));
    const auto              Wq = mat::solve(Q_process, Matrix<NX, NX, T>::identity());
    const auto              Wr = mat::solve(R_measurement, Matrix<NY, NY, T>::identity());
    if (!Wq || !Wr) {
        return set_status(art, DesignStatus::bad_weights); // a zero or indefinite covariance
    }

    const StateSpace<NX, NU, NY, T, NX, NY> sys_noise{.A = sys.A, .B = sys.B, .C = sys.C, .Ts = sys.Ts};
    const auto                              kf = kalman(sys_noise, KalmanNoise{.process = Q_process, .measurement = R_measurement});
    if (!kf.success) {
        return set_status(art, kf.status);
    }
    const auto Wa = mat::solve(kf.P, Matrix<NX, NX, T>::identity());
    if (!Wa) {
        return set_status(art, DesignStatus::singular);
    }

    mhe_hessian(art, sys.A, sys.C, *Wq, *Wr, *Wa);

    // Gradient maps for the runtime.
    art.M_a = T{2} * *Wa;
    art.M_cy = T{2} * (sys.C.transpose() * *Wr);
    art.M_qu = T{2} * (*Wq * sys.B);
    art.M_aqu = T{2} * (sys.A.transpose() * (*Wq * sys.B));

    // Arrival prior: one steady-state Kalman step over the sample leaving the window.
    art.Ap_x = sys.A * (Matrix<NX, NX, T>::identity() - (kf.M * sys.C));
    art.Ap_y = sys.A * kf.M;
    art.Ap_u = sys.B;

    mhe_state_rows(art);

    if (!mat::cholesky(art.H)) {
        return set_status(art, DesignStatus::singular);
    }
    return set_status(art, DesignStatus::ok);
}

/**
 * @brief design::mhe on the continuous plant @p sys sampled at @p Ts: its ZOH model, the process
 *        noise density discretized over the period (Van Loan) on G = I, the measurement noise as
 *        given, as design::kalman's sampled form reads them.
 */
template<size_t NH, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NVN, size_t NXC>
[[nodiscard]] constexpr MHEResult<NX, NU, NY, NH, T, NXC> mhe_sampled(
    const StateSpace<NX, NU, NY, T, NW, NVN>& sys,
    T                                         Ts,
    const KalmanNoise<NW, NVN, T>&            noise,
    const MHEStateSet<NXC, T>&                set,
    size_t                                    budget
) {
    using Result = MHEResult<NX, NU, NY, NH, T, NXC>;
    const auto m = sampled_model(sys, Ts, Domain::sampled);
    if (m.status != DesignStatus::ok) {
        return refused<Result>(m.status); // a discrete plant, a bad Ts, a non-finite entry or an overflow
    }
    const auto Qd = discretize_process_noise(sys.A, sys.G, noise.process, Ts);
    if (!Qd) {
        return refused<Result>(DesignStatus::bad_weights); // a non-finite noise entry
    }
    const StateSpace<NX, NU, NY, T, NX, NVN> sys_d{.A = m.law.A, .B = m.law.B, .C = sys.C, .D = sys.D, .G = Matrix<NX, NX, T>::identity(), .H = sys.H, .Ts = Ts};
    return mhe_design<NH>(sys_d, KalmanNoise<NX, NVN, T>{.process = *Qd, .measurement = noise.measurement}, set, budget);
}

/** @brief The result of design::mhe<NH> on a plant of NX states, NU inputs and NY outputs, with @p Opts. */
template<size_t NX, size_t NU, size_t NY, size_t NH, typename T, typename... Opts>
using mhe_result_t = MHEResult<NX, NU, NY, NH, T, listed_by<kind::StateLimits, Opts...>>;

} // namespace detail

/**
 * @brief Synthesize a constrained moving-horizon estimator on a discrete plant
 *
 * Builds the constant window Hessian (block tridiagonal from the dynamics and
 * measurement weights), the rows on the listed states, and the gradient maps. The noise
 * enters through the plant's G and H, as in design::kalman: the window weights are the
 * inverses of @f$ G Q G^\top + \delta_q I @f$ and @f$ H R H^\top + \delta_r I @f$, each
 * @f$ \delta = \sqrt{\varepsilon}\, \max_i W_{ii} @f$ of its own matrix. G Q Gᵀ is singular
 * whenever the noise drives fewer channels than states (a G with fewer columns than rows, or
 * the delay states of discretize(plant, Ts, InputDelay{θ}), which get no noise); δ keeps
 * the window solvable. The arrival weight is the inverse of the steady-state Kalman
 * covariance on the same two covariances, so the unconstrained estimator agrees with that
 * Kalman filter.
 *
 * | Design | Runtime |
 * | --- | --- |
 * | mhe<NH>(plant, [Ts,] noise) | MHE: the unconstrained window, which agrees with the steady-state Kalman filter |
 * | … plus StateLimits{.states, .lower, .upper} | the same, the states it lists bounded at every window state, soft by default (`.lower_ecr = {0.0}` for a hard lower end) |
 * | … plus QPIterations{n} | the same, its QP capped at n iterations a tick |
 *
 * The plant and the noise are positional, so braces size the noise from the plant:
 * `mhe<4>(tank, {.process = 1e-4, .measurement = 1e-3})`. The rest trails as options, at most
 * one of each kind, in any order: the limits on the states (none: no rows) and the QP budget
 * (none: 10·(NZ + NI), the window QP's variables and rows).
 *
 * Status: invalid_sample_time if the plant is not discrete (a continuous one takes
 * mhe(plant, Ts, …)); bad_weights on a non-finite model or noise entry, a feedthrough D ≠ 0, a
 * zero covariance or a negative or non-finite ECR; invalid_parameter on a StateLimits index past
 * the plant's states or listed twice, an end NaN or ±Inf, lower > upper, or a QPIterations{0};
 * design::kalman's status if the arrival filter fails; singular if P∞ or the window Hessian is
 * not positive definite.
 *
 * Q is the trust in the model (smaller = dynamics enforced harder); R the trust in the
 * sensor — the same covariances a Kalman filter design uses, and a good design reuses exactly
 * those numbers.
 *
 * @see StateLimits, QPIterations, damp::MHE, design::kalman
 *
 * @tparam NH     Estimation horizon (window covers NH+1 states / measurements)
 * @param plant   Discrete-time plant (Ts > 0, D = 0); G and H carry the noise
 * @param noise   Process-noise covariance Q per step (NW × NW) and measurement-noise variance R
 *                of one reading (NV × NV), positive semidefinite
 * @param opts    In any order: StateLimits{.states, .lower, .upper}, the limits on the states it
 *                lists (soft by default); QPIterations{n}, the per-tick QP budget
 * @return MHEResult: deploy with `MHE estimator{res}`
 */
template<size_t NH, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename... Opts>
[[nodiscard]] constexpr detail::mhe_result_t<NX, NU, NY, NH, T, Opts...> mhe(
    const StateSpace<NX, NU, NY, T, NW, NV>&      plant,
    const detail::plant_noise<NX, NY, T, NW, NV>& noise,
    const Opts&... opts
) {
    detail::mhe_options::check<Opts...>();
    using Result = detail::mhe_result_t<NX, NU, NY, NH, T, Opts...>;
    return detail::mhe_design<NH>(detail::noise_model(plant), noise, detail::mhe_state_set<T>(opts...), detail::mhe_budget_of<Result>(opts...));
}

/**
 * @brief Synthesize a constrained moving-horizon estimator from the continuous plant: the window
 *        on its ZOH model at @p Ts
 *
 * design::mhe(discretize(plant, Ts), …)'s window, the options as there, with the noise read as
 * design::kalman(plant, Ts, …) reads it: the process noise is a power spectral density, the
 * intensity of the white noise in ẋ = Ax + Bu + Gw, discretized exactly over the period (Van
 * Loan, discretize_process_noise), so the same density means the same at any rate; the
 * measurement noise is the variance of one reading, with no Ts scaling. The unconstrained window
 * then agrees with design::kalman(plant, Ts, KalmanNoise{…})'s filter.
 *
 * Refused with DesignStatus::invalid_sample_time when @p plant is already discrete or @p Ts is
 * not finite and positive, and as the discrete form is otherwise.
 *
 * @see "Computing Integrals Involving the Matrix Exponential" (Van Loan, 1978), IEEE TAC 23(3)
 *
 * @tparam NH     Estimation horizon (window covers NH+1 states / measurements)
 * @param plant   Continuous-time plant (Ts = 0, D = 0); G and H carry the noise
 * @param Ts      Sample period [s]
 * @param noise   Process-noise density Q (NW × NW) and measurement-noise variance R of one
 *                reading (NV × NV), positive semidefinite
 * @param opts    As the discrete form's
 * @return MHEResult: deploy with `MHE estimator{res}`
 */
template<size_t NH, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename... Opts>
[[nodiscard]] constexpr detail::mhe_result_t<NX, NU, NY, NH, T, Opts...> mhe(
    const StateSpace<NX, NU, NY, T, NW, NV>&      plant,
    damp::type_identity_t<T>                      Ts,
    const detail::plant_noise<NX, NY, T, NW, NV>& noise,
    const Opts&... opts
) {
    detail::mhe_options::check<Opts...>();
    using Result = detail::mhe_result_t<NX, NU, NY, NH, T, Opts...>;
    return detail::mhe_sampled<NH>(detail::noise_model(plant), Ts, noise, detail::mhe_state_set<T>(opts...), detail::mhe_budget_of<Result>(opts...));
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
 * oldest, NH the current estimate. limits() and set_limits() read and move the bounds on the
 * listed states at run time (design::MHELimits), and last_slack() reports how far the last
 * solve relaxed the soft ones.
 *
 * @tparam NXC    States with bounds (the design's StateLimits; 0: none, no rows)
 * @tparam Solver The QP solver (default: the warm-started active set of the design's QP)
 *
 * @see design::mhe, KalmanFilter for the unconstrained equivalent
 */
template<
    size_t NX, size_t NU, size_t NY, size_t NH, typename T = float, size_t NXC = 0,
    typename Solver = design::WarmStartActiveSetSolver<design::MHEResult<NX, NU, NY, NH, T, NXC>::NZ, design::MHEResult<NX, NU, NY, NH, T, NXC>::NI, T>>
    requires damp::is_floating_point_v<T>
class MHE {
public:
    using Result = design::MHEResult<NX, NU, NY, NH, T, NXC>; ///< The design it runs
    using Limits = design::MHELimits<NXC, T>;                 ///< What limits() returns and set_limits() takes
    static constexpr size_t NV = Result::NV;                  ///< Window states
    static constexpr size_t NZ = Result::NZ;                  ///< QP decision variables
    static constexpr size_t NI = Result::NI;                  ///< QP constraint rows

    constexpr MHE() = default;
    /**
     * @brief From design::mhe. A failed design does not compile when constant-evaluated; at run
     *        time valid() is false and step() holds the last estimate, skipping the QP.
     */
    constexpr MHE(const Result& result) // NOLINT(google-explicit-constructor): assign-from-design
        : art(result), valid_(result.success) {
        if (!result.success && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(result.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
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
        slide(y, u);
        const auto qp = qp_solver_(art.H, gradient(), art.A_con, art.b_bound, art.max_qp_iterations);
        last_status_ = qp.status;
        last_iterations_ = qp.iterations;
        if (qp.status == design::QPStatus::Success || qp.status == design::QPStatus::MaxIterations) {
            for (size_t v = 0; v < NV; ++v) {
                z(v) = qp.x(v);
            }
            if constexpr (NXC > 0) {
                last_slack_ = qp.x(NV);
                hold_hard_limits();
            }
        }
        return state();
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

    /**
     * @brief Move the bounds on the listed states (design::MHELimits::state, in the order of the
     *        design's StateLimits::states)
     *
     * The QP keeps its Hessian, constraint rows, scaling and softness; the bound vector takes the
     * new limits from the next step. Refused (false, limits unchanged) when an end is NaN or ±Inf
     * or a lower end is above its upper one.
     */
    [[nodiscard]] constexpr bool set_limits(const Limits& lim) {
        if (!lim.state.valid()) {
            return false;
        }
        detail::mhe_box_bounds(art, lim);
        art.limits = lim;
        return true;
    }

    /** @brief The limits in force (the design's, or the last set_limits). */
    [[nodiscard]] constexpr const Limits& limits() const { return art.limits; }

    constexpr void reset() {
        primed = false;
        z = ColVec<NV, T>{};
        x_bar = ColVec<NX, T>{};
        last_status_ = design::QPStatus::Success;
        last_iterations_ = 0;
        last_slack_ = T{0};
    }

    /** @brief The status of the last QP solve. */
    [[nodiscard]] constexpr design::QPStatus last_status() const { return last_status_; }
    /** @brief The iterations of the last QP solve. */
    [[nodiscard]] constexpr size_t last_iterations() const { return last_iterations_; }
    /** @brief Soft-limit slack ε from the last solve: 0 when no soft limit was relaxed (and always without a StateLimits). */
    [[nodiscard]] constexpr T last_slack() const { return last_slack_; }
    /** @brief Access the solver policy (stateful solvers expose warm-start caches here) */
    [[nodiscard]] constexpr Solver& qp_solver() { return qp_solver_; }

private:
    /** @brief Take the new sample into the window: the first primes it, each later one advances its start. */
    constexpr void slide(const ColVec<NY, T>& y, const ColVec<NU, T>& u) {
        if (!primed) {
            for (size_t k = 0; k <= NH; ++k) {
                ys[k] = y;
            }
            for (size_t k = 0; k < NH; ++k) {
                us[k] = u;
            }
            primed = true;
            return;
        }
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

    /** @brief The QP gradient from the ring buffers (the Hessian and the rows are constant; none on ε). */
    [[nodiscard]] constexpr ColVec<NZ, T> gradient() const {
        ColVec<NZ, T> f{};
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
        return f;
    }

    /**
     * @brief Project the window onto its hard ends (ECR 0), the part of the QP's feasible set a
     *        slack cannot widen: a solution leaves it by rounding on an active bound (−1e-17 for a
     *        bound of 0 when multiply-adds are fused), a dual iterate that ran out of iterations by
     *        more. A soft end holds for some ε, so the estimate is left as the QP gave it.
     */
    constexpr void hold_hard_limits() {
        for (size_t k = 0; k <= NH; ++k) {
            for (size_t j = 0; j < NXC; ++j) {
                T& x = z((k * NX) + art.states[j]);
                if (art.lower_ecr(j) == T{0}) {
                    x = damp::max(x, art.limits.state.lower(j));
                }
                if (art.upper_ecr(j) == T{0}) {
                    x = damp::min(x, art.limits.state.upper(j));
                }
            }
        }
    }

    Result                             art{};
    Solver                             qp_solver_{};
    damp::array<ColVec<NY, T>, NH + 1> ys{};
    damp::array<ColVec<NU, T>, NH>     us{};
    ColVec<NV, T>                      z{};
    ColVec<NX, T>                      x_bar{};
    design::QPStatus                   last_status_{design::QPStatus::Success};
    size_t                             last_iterations_{0};
    T                                  last_slack_{0};
    bool                               valid_{false};
    bool                               primed{false};
};

/** @brief Deduce the runtime from its design: MHE estimator{res}; */
template<size_t NX, size_t NU, size_t NY, size_t NH, typename T, size_t NXC>
MHE(const design::MHEResult<NX, NU, NY, NH, T, NXC>&) -> MHE<NX, NU, NY, NH, T, NXC>;

namespace detail {

/** @brief The runtime a design result runs as, in the scalar @p T: none for other results. */
template<typename R, typename T>
struct mhe_runtime {};

/** @brief design::mhe runs as MHE, its QP data in members (set_limits rewrites it). */
template<size_t NX, size_t NU, size_t NY, size_t NH, typename U, size_t NXC, typename T>
struct mhe_runtime<design::MHEResult<NX, NU, NY, NH, U, NXC>, T> {
    using type = MHE<NX, NU, NY, NH, T, NXC>; ///< The runtime
};

/** @brief A design::mhe design. */
template<const auto& Res, typename T>
concept mhe_design = requires { typename mhe_runtime<design_t<Res>, T>::type; };

} // namespace detail

/** @brief Law<kRes, T> of a design::mhe: its MHE, default-constructed from the design. A failed `constexpr` design does not compile, and the error names its status. */
template<const auto& Res, typename T>
    requires detail::mhe_design<Res, T>
struct Law<Res, T> : detail::mhe_runtime<detail::design_t<Res>, T>::type {
    using Runtime = detail::mhe_runtime<detail::design_t<Res>, T>::type; ///< The runtime the design runs as
    static_assert(detail::member_law_runs<Res>());                       // a failed constexpr design does not compile: the error names its status

    /** @brief The design's estimator. */
    constexpr Law() : Runtime{Res.template as<T>()} {}

    using Runtime::operator=;
};

} // namespace damp
