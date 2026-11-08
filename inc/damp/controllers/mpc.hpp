// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file mpc.hpp
 * @brief Constrained linear MPC: condensed-QP synthesis (`design::mpc`) and the fixed-budget
 *        runtime (`damp::MPC`), on the measured state or through a Kalman estimator, offset-free
 *        with a disturbance model, with an opt-in DARE terminal cost and Gilbert–Tan positively
 *        invariant terminal set.
 *
 * Finite-horizon constrained control of a discrete-time plant
 * @f$ x_{k+1} = A x_k + B u_k + B_d d_k,\; y_k = C x_k @f$ (roadmap #14), where
 * @f$ u @f$ are the manipulated inputs and @f$ d @f$ optional disturbances: measured ones
 * (feedforward inputs the controller knows but does not command) or the estimate of an
 * unmeasured one (offset-free control, see design::mpc).
 * The move increments @f$ \Delta u @f$ are the decision variables (the
 * "velocity form"): the state is augmented as @f$ \tilde x = [x;\, u_{prev};\, d] @f$
 * so the held input and the disturbance are carried by the prediction
 * model and the controller has built-in integral action. The cost over
 * prediction horizon NP / control horizon NC is
 * @f[
 *   J = \sum_{k=1}^{NP} \lVert y_k - r \rVert^2_{Q_y}
 *     + \sum_{k=0}^{NC-1} \lVert \Delta u_k \rVert^2_{R_{\Delta u}}
 *     + \sum_{k=0}^{NC-1} \lVert u_k - u_{target} \rVert^2_{R_u}
 *     + \rho_\varepsilon\, \varepsilon^2,
 * @f]
 * with @f$ Q_y @f$ replaced by @f$ Q_{y,N} @f$ on the final step (a terminal
 * output weight; set it to the DARE solution with @f$ C = I @f$ for the
 * guaranteed-stability short-horizon setup). Opt-in @ref damp::design::MPCTerminal "MPCTerminal::dare_cost"
 * replaces that last output weight with @f$ V_f = z_N^\top P z_N @f$, the cost-to-go of this
 * stage cost itself: a DARE on @f$ z = [x - x_t;\; u'_{prev} - u'_t] @f$ with the move as
 * input, so it is a control Lyapunov function for the Δu form. @ref damp::design::MPCTerminal "MPCTerminal::invariant_set"
 * adds the maximal positively invariant set of that terminal law, augmented with the constant
 * target (MPC for tracking), as hard inequalities on @f$ z_N @f$: one set valid for every
 * target that keeps @c terminal_target_margin of each half-range clear of its bounds, the
 * target moving each row's offset at run time. Default
 * @ref damp::design::MPCTerminal "MPCTerminal::none" is the plain Qy / Qy_terminal cost.
 * Subject to box constraints on @f$ u @f$, @f$ \Delta u @f$ (hard), and @f$ y @f$
 * (soft by default) over the horizon.
 *
 * That cost, with weights per step, is the discrete-plant form's. The continuous-plant form
 * design::mpc(plant, Ts, …) takes the weights per second and integrates the output error over
 * each period instead of sampling it,
 * @f[
 *   J = \sum_{k=0}^{NP-1} \int_{kT_s}^{(k+1)T_s} \lVert y(t) - r_k \rVert^2_{Q_y} dt
 *     + \sum_{k=0}^{NC-1} T_s \Bigl( \bigl\lVert \tfrac{\Delta u_k}{T_s} \bigr\rVert^2_{R_{\Delta u}}
 *     + \lVert u_k - u_{target} \rVert^2_{R_u} \Bigr) + \rho_\varepsilon\, \varepsilon^2,
 * @f]
 * with @f$ y(t) = C x(t) + D u_k @f$ under the held input, @f$ R_{\Delta u} @f$ on the rate
 * @f$ \Delta u/T_s @f$ (the Riemann sum of @f$ \int \dot u^\top R_{\Delta u} \dot u\, dt @f$)
 * and @f$ Q_{y,N} @f$ on the last period. Each period's integral is a quadratic form in
 * @f$ [x_k;\, u_k;\, d;\, r_k] @f$, cross terms included, computed exactly by one Van Loan
 * exponential, as the sampled-data LQ
 * designs discretize their cost. The weights then mean the same at any @f$ T_s @f$: halving
 * the period, with the horizons in steps doubled, leaves the closed loop nearly unchanged.
 *
 * The runtime poses the QP about the steady-state target (x_t, u_t) of the reference: the
 * terminal ingredients are on the deviation from it, and the gradient is a sum of small terms
 * however far the plant sits from zero. Without terminal ingredients this is the QP in
 * absolute coordinates, exactly in exact arithmetic.
 *
 * Soft constraints (ECR). A single slack @f$ \varepsilon \ge 0 @f$ relaxes
 * the output rows: @f$ y \le y_{max} + V\varepsilon @f$ with per-channel
 * softness coefficients V (`MPCConstraints::y_lower_ecr` / `y_upper_ecr`, default 1)
 * and penalty @f$ \rho_\varepsilon @f$ (`MPCWeights::ecr_weight`). V = 0 makes
 * a row hard. Input and rate constraints are always hard (actuator physics).
 * This is MATLAB®'s equal-concern-for-relaxation scheme: one slack, so an
 * over-tight output limit degrades gracefully instead of making the QP
 * infeasible.
 *
 * Scaling. Per-channel scale factors are derived automatically from the
 * constraint ranges (max − min where both are finite, else 1) and applied
 * internally, so the QP conditioning does not depend on physical units and the
 * weights act on dimensionless signals — MATLAB®'s ScaleFactor semantics.
 * Disturbances are not rescaled (unit scale).
 *
 * Weighting a state that is not an output: add it as an extra
 * (unmeasured) output row in C and weight that row in Qy — the toolbox-standard
 * recipe. On the measured state (no KalmanNoise) the measurement is the state x and the
 * plant's C still names the controlled outputs: C = [1 0] on a position–velocity plant tracks
 * the position with a scalar r and Qy. There is deliberately no separate state weight
 * (MATLAB®'s linear MPC has none either).
 *
 * r already is the tracking reference (NY). There is no LQ-style state-reference tick
 * `step(r, x_ref, u_ref, y)`: the QP has no ‖x − x_ref‖ term (it does not force x → x_ref except
 * as a consequence of a feasible y → r), a state reference at one instant does not give the
 * reference at the NP prediction steps the plan needs.
 *
 * Synthesis condenses the problem: predictions @f$ Y = \Phi \tilde x_0 +
 * \Gamma Z @f$ over the horizon give a dense strictly convex QP in
 * @f$ [Z; \varepsilon] @f$, @f$ Z = [\Delta u_0; \dots; \Delta u_{NC-1}] @f$,
 * solved each tick by the Goldfarb–Idnani active-set solver (design/qp.hpp)
 * under a hard iteration budget. Everything is fixed-size and allocation-free.
 *
 * Feedthrough is not supported: the controlled outputs must have @f$ D = 0 @f$ (synthesis
 * fails otherwise). Offset-free tracking under plant/model mismatch or unmeasured
 * disturbances needs a disturbance model and an estimator (design::ofmpc, ofmpc.hpp); the velocity form alone rejects only what its integral action sees
 * through the measured state.
 *
 * Unbounded constraint entries use the `numeric_limits` sentinels
 * (`-max()`/`max()`, never infinity — the library builds with
 * -ffinite-math-only); sentinel rows cost nothing in the solver.
 *
 * @note Compare with MATLAB®'s mpc(plant, Ts, p, m, weights) + mpcmove, with
 *       Weights.OV/MV/MVRate/ECR, ManipulatedVariables Targets, measured
 *       disturbances via setmpcsignals, and the built-in estimator (setindist, setoutdist).
 * @see design/qp.hpp for the QP solver
 * @see design/invariant_set.hpp for the Gilbert–Tan terminal set
 * @see estimation/kalman.hpp for the estimator's design, estimation/luenberger.hpp for the observer
 *      that runs it, and design/steady_state_target.hpp for the target
 * @see García, Prett & Morari, "Model Predictive Control: Theory and Practice —
 *      A Survey," Automatica 25(3), 1989, https://doi.org/10.1016/0005-1098(89)90002-2
 * @see Maciejowski, "Predictive Control with Constraints," Prentice Hall, 2002
 *      (condensed Δu formulation, ch. 2–3; soft constraints §3.4)
 * @see Rawlings, Mayne & Diehl, "Model Predictive Control: Theory, Computation,
 *      and Design," 2nd ed., Nob Hill, 2017
 * @see "Linear Systems with State and Control Constraints: The Theory and
 *      Application of Maximal Output Admissible Sets" (Gilbert & Tan, 1991)
 * @see "MPC for Tracking Piecewise Constant References for Constrained Linear Systems"
 *      (Limón, Alvarado, Alamo & Camacho, 2008), Automatica 44(9): the target-augmented set
 * @see "Computing Integrals Involving the Matrix Exponential" (Van Loan, 1978), IEEE TAC 23(3):
 *      the continuous-plant form's period integrals
 *
 * @code
 * // Double integrator, Ts = 0.1, the state measured: track the position with |u| ≤ 0.5, |Δu| ≤ 0.2.
 * constexpr StateSpace<2, 1, 1> sys{
 *     .A = {{1.0, 0.1}, {0.0, 1.0}},
 *     .B = {{0.005}, {0.1}},
 *     .C = {{1.0, 0.0}}, // the controlled output: the position
 *     .Ts = 0.1,
 * };
 * constexpr design::MPCConstraints<1, 1> limits{.u_lim = {.lower = -0.5, .upper = 0.5}, .du_lim = {.lower = -0.2, .upper = 0.2}};
 * constexpr auto art = design::mpc<15, 5>(sys, {}, limits);
 * static constinit MPC controller{art.as<float>()}; // a failed design does not compile
 * // each tick: u = controller.step(ColVec<1, float>{r}, x); // the state x measured
 * //            controller.back_calculate(u_applied);         // when a downstream stage clamps it
 *
 * // From the position reading alone, offset-free under an unmeasured load at the input:
 * constexpr StateSpace<2, 1, 1, double, 2, 1> plant{.A = sys.A, .B = sys.B, .C = sys.C, .Ts = 0.1};
 * constexpr auto of = design::ofmpc<15, 5>(plant, {}, {.u_lim = {.lower = -0.5, .upper = 0.5}}, design::KalmanNoise{Qn, Rn, Qd},
 *                                          design::InputDisturbance{}); // ofmpc.hpp
 * static constinit OFMPC offset_free{of.as<float>()};
 * // each tick: u = offset_free.step(ColVec<1, float>{r}, ColVec<1, float>{y});
 * @endcode
 */

#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/law_detail.hpp"
#include "damp/controllers/state_feedback_core.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/invariant_set.hpp"
#include "damp/design/options.hpp"
#include "damp/design/qp.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/status.hpp"
#include "damp/design/steady_state_target.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/estimation/luenberger.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/block.hpp"
#include "damp/matrix/decomposition.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/functions.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/svd.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/toolbox/bounds.hpp"

namespace damp {

namespace design {

namespace detail {

// filled_vector / convert_bounds (sentinel machinery) live in design/qp.hpp's
// detail — shared with the MHE and any other bound-carrying design function.

/** @brief Scale factor from a constraint range: (max − min) when both finite, else 1. */
template<typename T>
[[nodiscard]] constexpr T scale_from_range(T lo, T hi) {
    constexpr T thresh = std::numeric_limits<T>::max() / T{4};
    if (lo <= -thresh || lo >= thresh || hi >= thresh || hi <= -thresh) {
        return T{1};
    }
    const T span = hi - lo;
    return (span > T{0}) ? span : T{1};
}

/** @brief Divide a bound by a scale factor, leaving unbounded sentinels untouched. */
template<typename T>
[[nodiscard]] constexpr T scaled_bound(T bound, T scale) {
    constexpr T thresh = std::numeric_limits<T>::max() / T{4};
    if (bound >= thresh || bound <= -thresh) {
        return bound;
    }
    return bound / scale;
}

} // namespace detail

/**
 * @brief Opt-in terminal cost / set for `design::mpc` (default: none)
 *
 * `dare_cost` replaces the last-step Qy_terminal with V_f = z_Nᵀ P z_N, z = [x − x_t; u′_prev − u′_t],
 * the cost-to-go of the MPC's own Δu stage cost (needs R_Δu + R_u positive definite).
 * `invariant_set` adds that cost and the maximal positively invariant set of the terminal
 * law Δu′ = −K·z, augmented with the constant target, valid for targets inside the boxes
 * tightened by MPCConstraints::terminal_target_margin. That set has more facets than a
 * zero-target one (a 2-state double integrator ~48): raise NT_max to fit, and design it at
 * initialization, since constant-evaluated it exceeds GCC's default -fconstexpr-ops-limit.
 * dare_cost alone stays cheap enough for constinit.
 *
 * Both rest on shifting the last plan and appending the terminal move, which needs a free move
 * at every step: dare_cost's stability argument holds with NC = NP, and invariant_set is
 * refused (invalid_parameter) with NC < NP. The set's rows are hard and the target is the
 * steady state of r[0], so the QP is feasible only from states that reach the set within NP
 * steps for that target: a reference step outside that region makes every tick Infeasible, and
 * the runtime holds its input. Size NP for the reference range.
 */
enum class MPCTerminal : std::uint8_t {
    none = 0,      ///< Today's Qy / Qy_terminal behaviour, no extra set
    dare_cost,     ///< V_f: the Δu stage cost's cost-to-go on z, no extra set
    invariant_set, ///< That cost + the target-augmented invariant set (implies dare_cost; needs NC = NP)
};

/**
 * @struct MPCWeights
 * @brief Cost weights for mpc()
 *
 * Weights act on the internally scaled (dimensionless) signals, so they are
 * unit-free once constraint ranges are given — MATLAB®'s ScaleFactor semantics.
 *
 * Their time base is the plant's. On a discrete plant they are per prediction step: Qy on each
 * sampled output, RΔu on each move, Ru on each move's input. From a continuous plant
 * (design::mpc(plant, Ts, …)) they are per second: Qy and Ru weight the output error and the
 * input deviation integrated over time, RΔu the input rate Δu/Ts, so a change of Ts keeps
 * their meaning (see design::mpc).
 * All weight matrices are symmetric positive semidefinite (design::mpc refuses others with
 * bad_weights, as design::lqr does); the QP Hessian must come out positive definite, which holds
 * whenever RΔu ≻ 0 or Ru ≻ 0.
 *
 * @note Mirrors the MATLAB® mpc object's Weights.OV / Weights.MVRate /
 *       Weights.MV (+ MV Targets) / Weights.ECR.
 */
template<size_t NU, size_t NY, typename T = double>
struct MPCWeights {
    Matrix<NY, NY, T> Qy = Matrix<NY, NY, T>::identity();  ///< Output tracking weight, steps 1..NP−1 (per second from a continuous plant)
    Matrix<NU, NU, T> Rdu = Matrix<NU, NU, T>::identity(); ///< Move weight on Δu (on the rate Δu/Ts, per second, from a continuous plant)
    Matrix<NU, NU, T> Ru{};                                ///< Input-deviation weight on (u − u_target) (default 0; per second from a continuous plant)
    Matrix<NY, NY, T> Qy_terminal = Qy;                    ///< Output weight on the final step NP, or the last period (default: Qy)
    ColVec<NU, T>     u_target{};                          ///< Preferred input operating point (physical units)
    T                 ecr_weight{static_cast<T>(1.0e5)};   ///< Soft-constraint slack penalty ρ_ε

    /** @brief The weights in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr MPCWeights<NU, NY, U> as() const {
        return MPCWeights<NU, NY, U>{
            Qy.template as<U>(),
            Rdu.template as<U>(),
            Ru.template as<U>(),
            Qy_terminal.template as<U>(),
            ColVec<NU, U>{u_target.template as<U>()},
            static_cast<U>(ecr_weight),
        };
    }
};

/**
 * @struct MPCConstraints
 * @brief Box constraints for mpc()
 *
 * Bound defaults are unbounded (numeric_limits sentinels). u and Δu bounds are
 * hard and apply over the control horizon; y bounds apply over the whole
 * prediction horizon and are soft by default (ECR coefficient 1): y_lower_ecr(i) is the
 * softness of the y_lim.lower(i) rows and y_upper_ecr(i) that of the y_lim.upper(i) rows. Set
 * one to 0 to make that side of the channel hard. A larger ECR relaxes that bound more readily
 * relative to the others. design::mpc refuses a NaN limit or a lower limit above its upper one
 * (invalid_parameter), and a negative or non-finite ECR (bad_weights).
 *
 * @note Mirrors the MATLAB® mpc object's ManipulatedVariables Min/Max/RateMin/
 *       RateMax and OutputVariables Min/Max (+ MinECR/MaxECR).
 */
template<size_t NU, size_t NY, typename T = double>
struct MPCConstraints {
    Bounds<NU, T> u_lim{};                                      ///< Input bounds (hard)
    Bounds<NU, T> du_lim{};                                     ///< Move bounds on Δu (hard)
    Bounds<NY, T> y_lim{};                                      ///< Output bounds (soft through y_lower_ecr, y_upper_ecr)
    ColVec<NY, T> y_lower_ecr{detail::filled_vector<NY>(T{1})}; ///< Softness of the y_lim.lower rows, ≥ 0 (0 = hard; MATLAB's MinECR)
    ColVec<NY, T> y_upper_ecr{detail::filled_vector<NY>(T{1})}; ///< Softness of the y_lim.upper rows, ≥ 0 (0 = hard; MATLAB's MaxECR)
    MPCTerminal   terminal{MPCTerminal::none};                  ///< Opt-in DARE cost / Gilbert–Tan set
    T             terminal_target_margin{static_cast<T>(0.05)}; ///< invariant_set: half-range fraction a tracked target keeps clear of each bound
    size_t        max_qp_iterations{0};                         ///< Per-tick QP iteration budget (0: 10·(NZ+1 + NI+1))

    /** @brief The constraints in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr MPCConstraints<NU, NY, U> as() const {
        return MPCConstraints<NU, NY, U>{
            u_lim.template as<U>(),
            du_lim.template as<U>(),
            y_lim.template as<U>(),
            ColVec<NY, U>{y_lower_ecr.template as<U>()},
            ColVec<NY, U>{y_upper_ecr.template as<U>()},
            terminal,
            static_cast<U>(terminal_target_margin),
            max_qp_iterations,
        };
    }
};

/**
 * @brief The limits an MPC takes at run time (MATLAB®'s online constraints): u, Δu and the
 *        controlled outputs, in physical units. A field left out is unbounded.
 *
 * The softness (ECR), the terminal ingredients and the QP budget stay in MPCConstraints: they
 * are part of the design. `auto l = mpc.limits(); l.u_lim.upper = …; mpc.set_limits(l);` changes one.
 */
template<size_t NU, size_t NY, typename T = double>
struct MPCLimits {
    Bounds<NU, T> u_lim{};  ///< Input bounds
    Bounds<NU, T> du_lim{}; ///< Move bounds on Δu
    Bounds<NY, T> y_lim{};  ///< Output bounds
};

/**
 * @struct MPCArtifacts
 * @brief Condensed-QP data produced by mpc() (MPCResult::mpc), consumed by damp::MPC
 *
 * Holds the dense QP in [Z; ε] (move increments + soft-constraint slack):
 * Hessian H, gradient maps (f = F_x·x̃0 − F_r·r' − F_ut·u'_target with
 * x̃0 = [x; u'_prev; d]), the stacked constraint matrix A_con with its
 * state-dependent bounds (b = b_bound + B_map·x̃0), the per-channel scale
 * factors, and the physical clamping limits. Primes denote internally scaled
 * quantities. Use .as<float>() for embedded deployment (bound sentinels
 * convert safely).
 *
 * @tparam NX       States of the prediction model
 * @tparam NY       Controlled outputs
 * @tparam ND       Disturbance channels: the measured disturbances, or the estimated ones
 */
template<size_t NX, size_t NU, size_t NY, size_t NP, size_t NC, size_t ND = 0, typename T = double, size_t NT_max = 0>
struct [[nodiscard]] MPCArtifacts {
    static constexpr size_t NXA = NX + NU + ND;                     ///< Augmented state size [x; u'_prev; d]
    static constexpr size_t NYS = NP * NY;                          ///< Stacked outputs over the prediction horizon
    static constexpr size_t NZ = NC * NU;                           ///< Move-increment decision variables
    static constexpr size_t NZ1 = NZ + 1;                           ///< QP decision size (moves + slack ε)
    static constexpr size_t NI = (4 * NZ) + (2 * NP * NY) + NT_max; ///< Box-constraint rows (+ terminal capacity)
    static constexpr size_t NI1 = NI + 1;                           ///< QP constraint rows (+ ε ≥ 0)
    static constexpr size_t NI_box = (4 * NZ) + (2 * NP * NY);      ///< Rows before the terminal block

    Matrix<NZ1, NZ1, T>         H{};                                 ///< QP Hessian (positive definite)
    Matrix<NZ, NXA, T>          F_x{};                               ///< Gradient state map
    Matrix<NZ, NYS, T>          F_r{};                               ///< Gradient reference map, one scaled reference per prediction step (preview)
    Matrix<NZ, NZ, T>           F_ut{};                              ///< Gradient input-target map, one scaled target per move
    Matrix<NI1, NZ1, T>         A_con{};                             ///< Constraint rows: A_con·[Z;ε] ≤ b_bound + B_map·x̃0
    ColVec<NI1, T>              b_bound{};                           ///< Constant bound term (sentinel rows disabled)
    Matrix<NI1, NXA, T>         B_map{};                             ///< State-dependent bound term
    Matrix<NY, NX, T>           C{};                                 ///< Output map of the prediction model (the target's output)
    Matrix<NY, ND, T>           Cd{};                                ///< The disturbance channels into the outputs (C_{d,z}; zero for measured ones)
    ColVec<NU, T>               scale_u{};                           ///< Input scale factors
    ColVec<NY, T>               scale_y{};                           ///< Output scale factors
    ColVec<NU, T>               u_target{};                          ///< Initial input target (physical units)
    MPCConstraints<NU, NY, T>   constraints{};                       ///< Original limits (runtime clamp fallback)
    size_t                      max_qp_iterations{};                 ///< Per-tick QP iteration budget
    Matrix<NX + NU, NX + NU, T> P{};                                 ///< Terminal cost on z = [x − x_t; u′_prev − u′_t] (zero when unused)
    Matrix<NU, NX + NU, T>      K_terminal{};                        ///< Terminal move law Δu′ = −K·z
    Matrix<NT_max, NX + NU, T>  G_terminal{};                        ///< Terminal-set rows on z; first n_terminal used
    Matrix<NT_max, NY + NU, T>  G_terminal_target{};                 ///< Their target terms: G·z + G_target·[y_t; u_t] ≤ h
    ColVec<NT_max, T>           h_terminal{};                        ///< Terminal-set offsets
    size_t                      n_terminal{0};                       ///< Live terminal-set rows
    MPCTerminal                 terminal{MPCTerminal::none};         ///< Requested terminal mode
    DesignStatus                status{DesignStatus::not_converged}; ///< ok, or why the design failed
    bool                        success{false};                      ///< status == DesignStatus::ok
    T                           Ts{0};                               ///< Sample time of the law [s]

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The artifacts in another scalar; refused (singular) when H does not factor in it. */
    template<typename U>
    [[nodiscard]] constexpr MPCArtifacts<NX, NU, NY, NP, NC, ND, U, NT_max> as() const {
        MPCArtifacts<NX, NU, NY, NP, NC, ND, U, NT_max> out{
            .H = H.template as<U>(),
            .F_x = F_x.template as<U>(),
            .F_r = F_r.template as<U>(),
            .F_ut = F_ut.template as<U>(),
            .A_con = A_con.template as<U>(),
            .b_bound = detail::convert_bounds<U>(b_bound),
            .B_map = B_map.template as<U>(),
            .C = C.template as<U>(),
            .Cd = Cd.template as<U>(),
            .scale_u = ColVec<NU, U>{scale_u.template as<U>()},
            .scale_y = ColVec<NY, U>{scale_y.template as<U>()},
            .u_target = ColVec<NU, U>{u_target.template as<U>()},
            .constraints = constraints.template as<U>(),
            .max_qp_iterations = max_qp_iterations,
            .P = P.template as<U>(),
            .K_terminal = K_terminal.template as<U>(),
            .G_terminal = G_terminal.template as<U>(),
            .G_terminal_target = G_terminal_target.template as<U>(),
            .h_terminal = detail::convert_bounds<U>(h_terminal),
            .n_terminal = n_terminal,
            .terminal = terminal,
            .status = status,
            .success = success,
            .Ts = static_cast<U>(Ts),
        };
        //! The design certified chol(H) in its own scalar; a Hessian near 1/ε_U in condition can
        //! fail in U, and the runtime would then hold its input on every tick.
        if (out.success && !mat::cholesky(out.H)) {
            out.success = false;
            out.status = DesignStatus::singular;
        }
        return out;
    }

    /**
     * @brief The unconstrained receding-horizon law as a discrete controller block: state u_prev,
     *        inputs [r; d; x], output u
     *
     * While no constraint is active the QP's minimizer is linear in its data,
     * @f$ Z^* = -H^{-1}(F_x \tilde x_0 - F_r R' - F_{ut} U'_t) @f$, so the first move is the gain
     * @f$ \Delta u = -K_x x - K_u u_{prev} - K_d d + K_r r @f$ in physical units (r held over the
     * horizon, @f$ K_r @f$ the sum of the preview gains). The runtime applies
     * @f$ u = u_{prev} + \Delta u @f$ and holds it as the next @f$ u_{prev} @f$:
     * @f[
     *   A_c = I - K_u,\quad B_c = [K_r\;\; -K_d\;\; -K_x],\quad C_c = A_c,\quad D_c = B_c .
     * @f]
     * The inputs are the tick's arguments in order: r for the NY controlled outputs, the ND
     * disturbance channels (none when ND = 0), and the state x of the law's model the QP is
     * solved from. Close it around a plant that reports [0; 0; x] (the reference and the
     * disturbance entering through the first rows) for the loop's poles, margins and DC gain,
     * as an LQ result's to_ss() is closed. The input target enters the law as a constant offset,
     * which a linear block does not carry: the block is the law at u_target = 0.
     *
     * Exact for MPCTerminal::none. A terminal cost is posed about the target of the reference,
     * which these artifacts do not hold: the block keeps the feedback from x and u_prev (the
     * loop's poles and margins) but leaves the terminal cost's pull to that target out of its r
     * and d gains. A design whose Hessian does not factor gives the zero block.
     *
     * @note Compare with MATLAB®'s ss(mpcobj): the controller while no constraint is active.
     *
     * @return StateSpace with NU states (u_prev), NY + ND + NX inputs ([r; d; x]), NU outputs (u)
     */
    [[nodiscard]] constexpr StateSpace<NU, NY + ND + NX, NU, T> to_ss() const {
        constexpr size_t                         NIN = NY + ND + NX;
        StateSpace<NU, NIN, NU, T>               out{.Ts = Ts};
        const Matrix<NZ, NZ, T>                  Hz = H.template block<NZ, NZ>(0, 0).to_matrix();
        const damp::optional<Matrix<NZ, NXA, T>> Gx = mat::cholesky_solve(Hz, F_x);
        const damp::optional<Matrix<NZ, NYS, T>> Gr = mat::cholesky_solve(Hz, F_r);
        if (!Gx || !Gr) {
            return out;
        }
        //! Back to physical units: Δu = S_u·Δu′, x̃ = [x; S_u⁻¹u_prev; d], r′ = S_y⁻¹r, r′ held at every step.
        Matrix<NU, NIN, T> Bc{};
        Matrix<NU, NU, T>  Ac = Matrix<NU, NU, T>::identity();
        for (size_t i = 0; i < NU; ++i) {
            const T su = scale_u(i);
            for (size_t j = 0; j < NY; ++j) {
                T held{0};
                for (size_t kb = 0; kb < NP; ++kb) {
                    held += (*Gr)(i, (kb * NY) + j);
                }
                Bc(i, j) = (su * held) / scale_y(j);
            }
            for (size_t j = 0; j < ND; ++j) {
                Bc(i, NY + j) = -su * (*Gx)(i, NX + NU + j);
            }
            for (size_t j = 0; j < NX; ++j) {
                Bc(i, NY + ND + j) = -su * (*Gx)(i, j);
            }
            for (size_t j = 0; j < NU; ++j) {
                Ac(i, j) -= (su * (*Gx)(i, NX + j)) / scale_u(j);
            }
        }
        out.A = Ac;
        out.B = Bc;
        out.C = Ac;
        out.D = Bc;
        return out;
    }
};

namespace detail {

/** @brief Limits an MPC accepts: every end not NaN (∓max() for none) and no lower end above its upper one. */
template<size_t N, typename T>
[[nodiscard]] constexpr bool mpc_limits_valid(const Bounds<N, T>& b) {
    for (size_t j = 0; j < N; ++j) {
        if (!damp::isfinite(b.lower(j)) || !damp::isfinite(b.upper(j)) || !(b.lower(j) <= b.upper(j))) {
            return false;
        }
    }
    return true;
}

/** @brief A weight matrix design::mpc accepts: finite, symmetric and positive semidefinite (as design::lqr's). */
template<size_t N, typename T>
[[nodiscard]] constexpr bool mpc_weight_valid(const Matrix<N, N, T>& W) {
    return mat::is_symmetric_or_hermitian(W) && ::damp::detail::is_psd(W);
}

/** @brief The weights and softness design::mpc accepts: PSD weights, a finite positive ECR weight, finite ECR coefficients ≥ 0. */
template<size_t NU, size_t NY, typename T>
[[nodiscard]] constexpr bool mpc_weights_valid(const MPCWeights<NU, NY, T>& w, const MPCConstraints<NU, NY, T>& c) {
    const auto softness = [](T v) { return damp::isfinite(v) && v >= T{0}; };
    return mpc_weight_valid(w.Qy) && mpc_weight_valid(w.Qy_terminal) && mpc_weight_valid(w.Rdu) && mpc_weight_valid(w.Ru)
        && damp::isfinite(w.ecr_weight) && w.ecr_weight > T{0} && mat::all_of(w.u_target, [](T v) { return damp::isfinite(v); })
        && mat::all_of(c.y_lower_ecr, softness) && mat::all_of(c.y_upper_ecr, softness);
}

/**
 * @brief Write the box rows of the QP's bound vector from @p lim with the design's scales: Δu and
 *        u over the control horizon, the controlled outputs over the prediction horizon.
 */
template<size_t NX, size_t NU, size_t NY, size_t NP, size_t NC, size_t ND, typename T, size_t NT_max>
constexpr void mpc_box_bounds(MPCArtifacts<NX, NU, NY, NP, NC, ND, T, NT_max>& art, const MPCLimits<NU, NY, T>& lim) {
    constexpr size_t NZ = NC * NU;
    constexpr size_t NYS = NP * NY;
    constexpr size_t r_du_hi = 0;
    constexpr size_t r_du_lo = NZ;
    constexpr size_t r_u_hi = 2 * NZ;
    constexpr size_t r_u_lo = 3 * NZ;
    constexpr size_t r_y_hi = 4 * NZ;
    constexpr size_t r_y_lo = (4 * NZ) + NYS;
    for (size_t i = 0; i < NC; ++i) {
        for (size_t j = 0; j < NU; ++j) {
            const T su = art.scale_u(j);
            art.b_bound(r_du_hi + (i * NU) + j) = scaled_bound(lim.du_lim.upper(j), su);
            art.b_bound(r_du_lo + (i * NU) + j) = scaled_bound(-lim.du_lim.lower(j), su);
            art.b_bound(r_u_hi + (i * NU) + j) = scaled_bound(lim.u_lim.upper(j), su);
            art.b_bound(r_u_lo + (i * NU) + j) = scaled_bound(-lim.u_lim.lower(j), su);
        }
    }
    for (size_t kb = 0; kb < NP; ++kb) {
        for (size_t j = 0; j < NY; ++j) {
            const T sy = art.scale_y(j);
            art.b_bound(r_y_hi + (kb * NY) + j) = scaled_bound(lim.y_lim.upper(j), sy);
            art.b_bound(r_y_lo + (kb * NY) + j) = scaled_bound(-lim.y_lim.lower(j), sy);
        }
    }
}

template<size_t NP, size_t NC, size_t NXA, size_t NX, size_t NU, size_t NY, typename T>
struct MpcPred {
    static constexpr size_t NYS = NP * NY;
    static constexpr size_t NZ = NC * NU;
    Matrix<NYS, NXA, T>     Phi{};
    Matrix<NYS, NZ, T>      Gamma{};
    Matrix<NX, NXA, T>      Phi_x{};
    Matrix<NX, NZ, T>       Gamma_x{};
};

template<size_t NP, size_t NC, size_t NXA, size_t NX, size_t NU, size_t NY, typename T>
[[nodiscard]] constexpr MpcPred<NP, NC, NXA, NX, NU, NY, T> mpc_predictions(
    const Matrix<NXA, NXA, T>& Aa,
    const Matrix<NXA, NU, T>&  Ba,
    const Matrix<NY, NXA, T>&  Ca
) {
    MpcPred<NP, NC, NXA, NX, NU, NY, T> pred{};
    damp::array<Matrix<NY, NU, T>, NP>  CAB{};
    damp::array<Matrix<NX, NU, T>, NP>  EAB{};
    Matrix<NXA, NXA, T>                 Apow = Matrix<NXA, NXA, T>::identity();
    for (size_t m = 0; m < NP; ++m) {
        const auto ABa = Apow * Ba;
        CAB[m] = Ca * ABa;
        EAB[m] = ABa.template block<NX, NU>(0, 0);
        Apow = Apow * Aa;
        pred.Phi.template block<NY, NXA>(m * NY, 0) = Ca * Apow;
    }
    pred.Phi_x = Apow.template block<NX, NXA>(0, 0);
    for (size_t k = 1; k <= NP; ++k) {
        const size_t jmax = (k < NC) ? k : NC;
        for (size_t j = 0; j < jmax; ++j) {
            pred.Gamma.template block<NY, NU>((k - 1) * NY, j * NU) = CAB[k - 1 - j];
        }
    }
    for (size_t j = 0; j < NC; ++j) {
        pred.Gamma_x.template block<NX, NU>(0, j * NU) = EAB[NP - 1 - j];
    }
    return pred;
}

/** @brief Terminal dynamics on z = [e; v], e = x − x_t, v = u′_prev − u′_t, with the move Δu′ as input. */
template<size_t NX, size_t NU, typename T>
struct MpcTerminalModel {
    Matrix<NX + NU, NX + NU, T> A{}; ///< [A B′; 0 I]
    Matrix<NX + NU, NU, T>      B{}; ///< [B′; I]
};

template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr MpcTerminalModel<NX, NU, T> mpc_terminal_model(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& Bs) {
    MpcTerminalModel<NX, NU, T> m{};
    m.A.template block<NX, NX>(0, 0) = A;
    m.A.template block<NX, NU>(0, NX) = Bs;
    m.A.template block<NU, NU>(NX, NX) = Matrix<NU, NU, T>::identity();
    m.B.template block<NX, NU>(0, 0) = Bs;
    m.B.template block<NU, NU>(NX, 0) = Matrix<NU, NU, T>::identity();
    return m;
}

/**
 * @brief The output part of one stage's cost in the deviation from the target: a quadratic form
 *        eᵀ·Q_x·e + 2·eᵀ·N_xu·u′ + u′ᵀ·R_uu·u′ on the state e and the input u′ applied over it.
 */
template<size_t NX, size_t NU, typename T>
struct MpcStageCost {
    Matrix<NX, NX, T> Q_x{};  ///< On the state
    Matrix<NX, NU, T> N_xu{}; ///< State–input cross term
    Matrix<NU, NU, T> R_uu{}; ///< On the input
};

/**
 * @brief The terminal cost-to-go of the MPC's own stage cost: an LQR on z = [e; v] with input Δu′
 *
 * Past the horizon the stage cost is ℓ(e, u′) + Δᵀ·R_Δu·Δ + (v + Δ)ᵀ·R_u·(v + Δ), the applied
 * input being u′ = v + Δ, with ℓ the output cost of @p stage: (C′e)ᵀQ_y(C′e) on a discrete plant,
 * the integral over the period on a continuous one. So Q = [Q_x N_xu; N_xuᵀ R_uu + R_u], cross
 * term N = [N_xu; R_uu + R_u], R = R_uu + R_Δu + R_u. P is then the exact infinite-horizon
 * cost-to-go, so V_f = zᵀPz is a control Lyapunov function for the stage cost, and K the
 * terminal move law Δu′ = −K·z.
 */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr DesignStatus mpc_dare_lqr(
    const Matrix<NX, NX, T>&       A,
    const Matrix<NX, NU, T>&       Bs,
    const MpcStageCost<NX, NU, T>& stage,
    const Matrix<NU, NU, T>&       Rdu,
    const Matrix<NU, NU, T>&       Ru,
    Matrix<NX + NU, NX + NU, T>&   P,
    Matrix<NU, NX + NU, T>&        K
) {
    constexpr size_t NZD = NX + NU;
    if (!mat::cholesky(Matrix<NU, NU, T>(Rdu + Ru))) {
        return DesignStatus::bad_weights; // R_Δu + R_u must be positive definite
    }
    const Matrix<NU, NU, T> Rv = stage.R_uu + Ru;
    const auto              m = mpc_terminal_model(A, Bs);
    Matrix<NZD, NZD, T>     Q{};
    Q.template block<NX, NX>(0, 0) = stage.Q_x;
    Q.template block<NX, NU>(0, NX) = stage.N_xu;
    Q.template block<NU, NX>(NX, 0) = stage.N_xu.transpose();
    Q.template block<NU, NU>(NX, NX) = Rv;
    Matrix<NZD, NU, T> N{};
    N.template block<NX, NU>(0, 0) = stage.N_xu;
    N.template block<NU, NU>(NX, 0) = Rv;
    const auto lqr = discrete_lqr(m.A, m.B, Q, Matrix<NU, NU, T>(Rv + Rdu), N);
    if (lqr.status != DesignStatus::ok) {
        return lqr.status;
    }
    if (!lqr.is_stable()) {
        return DesignStatus::not_stabilizable;
    }
    P = lqr.S;
    K = lqr.K;
    return DesignStatus::ok;
}

template<size_t NG, size_t NX, typename T>
constexpr void mpc_push_box_row(
    Matrix<NG, NX, T>&   G0,
    ColVec<NG, T>&       h0,
    size_t&              n,
    const ColVec<NX, T>& a,
    T                    beta
) {
    if (bound_is_unbounded(beta) || n >= NG) {
        return;
    }
    for (size_t j = 0; j < NX; ++j) {
        G0(n, j) = a(j);
    }
    h0(n) = beta;
    ++n;
}

/**
 * @brief Rows of the target-augmented terminal set on w = [z; q], z = [e; v], q = [y_t; u_t]:
 * y = y_t + C·e, u = u_t + S_u·(v − K·z), Δu = −S_u·K·z inside their boxes, and the target
 * q inside the boxes tightened by @p margin of each half-range (the Gilbert–Tan condition
 * that makes the set finitely determined with the constant target's unit eigenvalue).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NA = NX + NU + NY + NU, size_t NG = (4 * NY) + (6 * NU)>
constexpr void mpc_terminal_rows(
    const Matrix<NY, NX, T>&         C,
    const Matrix<NU, NX + NU, T>&    K,
    const ColVec<NU, T>&             scale_u,
    const ColVec<NY, T>&             scale_y,
    const MPCConstraints<NU, NY, T>& c,
    Matrix<NG, NA, T>&               G0,
    ColVec<NG, T>&                   h0
) {
    constexpr size_t NZD = NX + NU;
    h0 = filled_vector<NG>(unbounded_bound<T>());
    size_t     n = 0;
    const T    half = static_cast<T>(0.5) * c.terminal_target_margin;
    const auto push_pair = [&](const ColVec<NA, T>& a, T hi, T lo) {
        mpc_push_box_row(G0, h0, n, a, hi);
        mpc_push_box_row(G0, h0, n, ColVec<NA, T>(a * T{-1}), -lo);
    };
    for (size_t i = 0; i < NY; ++i) {
        ColVec<NA, T> a{};
        for (size_t j = 0; j < NX; ++j) {
            a(j) = C(i, j);
        }
        a(NZD + i) = T{1};
        push_pair(a, c.y_lim.upper(i), c.y_lim.lower(i));
        ColVec<NA, T> t{};
        t(NZD + i) = T{1};
        const T m = half * scale_y(i);
        push_pair(t, c.y_lim.upper(i) - m, c.y_lim.lower(i) + m);
    }
    for (size_t i = 0; i < NU; ++i) {
        ColVec<NA, T> a{};
        ColVec<NA, T> da{};
        for (size_t j = 0; j < NZD; ++j) {
            a(j) = -scale_u(i) * K(i, j);
            da(j) = -scale_u(i) * K(i, j);
        }
        a(NX + i) += scale_u(i);
        a(NZD + NY + i) = T{1};
        push_pair(a, c.u_lim.upper(i), c.u_lim.lower(i));
        push_pair(da, c.du_lim.upper(i), c.du_lim.lower(i));
        ColVec<NA, T> t{};
        t(NZD + NY + i) = T{1};
        const T m = half * scale_u(i);
        push_pair(t, c.u_lim.upper(i) - m, c.u_lim.lower(i) + m);
    }
}

template<size_t NP, size_t NC, size_t ND, size_t NT_max, size_t NX, size_t NU, size_t NY, typename T>
constexpr void mpc_stack_terminal_rows(
    MPCArtifacts<NX, NU, NY, NP, NC, ND, T, NT_max>& art,
    const Matrix<NX + NU, NX + NU + ND, T>&          Phi_z,
    const Matrix<NX + NU, NC * NU, T>&               Gamma_z
) {
    constexpr size_t NZ = NC * NU;
    constexpr size_t NXA = NX + NU + ND;
    constexpr size_t r_term = (4 * NZ) + (2 * NP * NY);
    if constexpr (NT_max > 0) {
        art.A_con.template block<NT_max, NZ>(r_term, 0) = art.G_terminal * Gamma_z;
        art.B_map.template block<NT_max, NXA>(r_term, 0) = -(art.G_terminal * Phi_z);
        for (size_t i = 0; i < NT_max; ++i) {
            art.b_bound(r_term + i) = (i < art.n_terminal) ? art.h_terminal(i) : unbounded_bound<T>();
        }
    }
}

/**
 * @brief The continuous plant a design from it integrates its stage cost over: the dynamics, the
 *        controlled outputs, the disturbance channel and the period. Off on a discrete plant.
 */
template<size_t NXP, size_t NU, size_t NY, size_t ND, typename T>
struct MpcContinuousCost {
    Matrix<NXP, NXP, T> A{};       ///< Continuous state matrix
    Matrix<NXP, NU, T>  B{};       ///< Continuous input matrix
    Matrix<NXP, ND, T>  Bd{};      ///< Continuous disturbance channel
    Matrix<NY, NXP, T>  C{};       ///< Controlled outputs
    Matrix<NY, NU, T>   D{};       ///< Their feedthrough
    T                   Ts{0};     ///< Period [s]
    bool                on{false}; ///< Whether the design is from the continuous plant
};

/**
 * @brief The output cost of one period, ∫ (y′(t) − r′)ᵀ·Q_y·(y′(t) − r′) dt over [0, Ts), as a
 *        quadratic form on s = [x; u′; d; r′] at the period's start (Van Loan)
 *
 * Over the period the input is held, u = S_u·u′; d and r′ are constant, y′ = S_y⁻¹·(C·x + D·u).
 * The reference and the disturbance are carried as constant states of the plant, so one Van Loan
 * exponential gives every cross term exactly: the state's
 * motion within the period under the held input, the disturbance, and the reference.
 *
 * @see "Computing Integrals Involving the Matrix Exponential" (Van Loan, 1978), IEEE TAC 23(3)
 */
template<size_t NXP, size_t NU, size_t NY, size_t ND, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NXP + NU + ND + NY, NXP + NU + ND + NY, T>> mpc_period_cost(
    const MpcContinuousCost<NXP, NU, NY, ND, T>& plant,
    const Matrix<NY, NY, T>&                     Qy,
    const ColVec<NU, T>&                         scale_u,
    const ColVec<NY, T>&                         scale_y
) {
    constexpr size_t    NPS = NXP + ND + NY;     // p = [x; d; r′], d and r′ held
    constexpr size_t    NW = NPS + NU;           // w = [p; u], the Van Loan order
    constexpr size_t    NS = NXP + NU + ND + NY; // s = [x; u′; d; r′]
    Matrix<NPS, NPS, T> Ap{};
    Ap.template block<NXP, NXP>(0, 0) = plant.A;
    Ap.template block<NXP, ND>(0, NXP) = plant.Bd;
    Matrix<NPS, NU, T> Bp{};
    Bp.template block<NXP, NU>(0, 0) = plant.B;
    //! The scaled output error e = E·p + F·u: E = [S_y⁻¹C, 0, −I], F = S_y⁻¹D
    Matrix<NY, NPS, T> E{};
    Matrix<NY, NU, T>  F{};
    for (size_t i = 0; i < NY; ++i) {
        for (size_t j = 0; j < NXP; ++j) {
            E(i, j) = plant.C(i, j) / scale_y(i);
        }
        E(i, NXP + ND + i) = T{-1};
        for (size_t j = 0; j < NU; ++j) {
            F(i, j) = plant.D(i, j) / scale_y(i);
        }
    }
    const Matrix<NPS, NY, T>  EtQ = E.transpose() * Qy;
    const LQRCost<NPS, NU, T> integrand{.Q = EtQ * E, .R = F.transpose() * (Qy * F), .N = EtQ * F};
    Matrix<NW, NW, T>         Wp{};
    const auto                c_opt = discretize_lqr_cost(Ap, Bp, integrand.Q, integrand.R, plant.Ts, integrand.N);
    if (!c_opt) {
        return damp::nullopt;
    }
    const auto& c = *c_opt;
    Wp.template block<NPS, NPS>(0, 0) = c.Q;
    Wp.template block<NPS, NU>(0, NPS) = c.N;
    Wp.template block<NU, NPS>(NPS, 0) = c.N.transpose();
    Wp.template block<NU, NU>(NPS, NPS) = c.R;
    //! w = P·s: the QP's order [x; u′; d; r′], the input scaled, u = S_u·u′
    Matrix<NW, NS, T> P{};
    for (size_t j = 0; j < NXP; ++j) {
        P(j, j) = T{1};
    }
    for (size_t j = 0; j < ND; ++j) {
        P(NXP + j, NXP + NU + j) = T{1};
    }
    for (size_t j = 0; j < NY; ++j) {
        P(NXP + ND + j, NXP + NU + ND + j) = T{1};
    }
    for (size_t j = 0; j < NU; ++j) {
        P(NPS + j, NXP + j) = scale_u(j);
    }
    const Matrix<NS, NS, T> W = P.transpose() * (Wp * P);
    return (W + W.transpose()) * static_cast<T>(0.5);
}

/** @brief The output part of the stage cost on the law's state and the scaled input (r′ = 0, d = 0), from a period cost on [x; u′; d; r′]. */
template<size_t NX, size_t NU, size_t NS, typename T>
[[nodiscard]] constexpr MpcStageCost<NX, NU, T> mpc_stage_of(const Matrix<NS, NS, T>& W) {
    return MpcStageCost<NX, NU, T>{
        .Q_x = W.template block<NX, NX>(0, 0),
        .N_xu = W.template block<NX, NU>(0, NX),
        .R_uu = W.template block<NU, NU>(NX, NX),
    };
}

/**
 * @brief The output cost of the horizon in the condensed QP: on the sampled outputs y′₁ … y′_NP
 *        (a discrete plant), or on the period integrals of a continuous one.
 *
 * Sets ½·the Hessian @p Hm and ½·the state gradient map @p Fx of the output cost, and art.F_r.
 * The continuous cost of period k is s_kᵀ·W·s_k on s_k = [q_k; r′_k], q_k = [x_k; u′_k; d] with
 * u′_k the input held over the period (the input rows of x̃_{k+1}); the last period is on
 * Qy_terminal unless a terminal cost replaces it.
 *
 * @return false if a period's exponential fails (a non-finite weight)
 */
template<size_t NP, size_t NC, size_t NXP, size_t NX, size_t NU, size_t NY, size_t ND, size_t NT_max, typename T, size_t NXA = NX + NU + ND, size_t NZ = NC * NU>
[[nodiscard]] constexpr bool mpc_output_cost(
    MPCArtifacts<NX, NU, NY, NP, NC, ND, T, NT_max>& art,
    const MpcPred<NP, NC, NXA, NX, NU, NY, T>&       pred,
    const Matrix<NXA, NXA, T>&                       Aa,
    const Matrix<NXA, NU, T>&                        Ba,
    const MPCWeights<NU, NY, T>&                     weights,
    const MpcContinuousCost<NXP, NU, NY, ND, T>&     plant,
    bool                                             dare_on,
    Matrix<NZ, NZ, T>&                               Hm,
    Matrix<NZ, NXA, T>&                              Fx
) {
    if (!plant.on) {
        const Matrix<NY, NY, T> Q_last = dare_on ? Matrix<NY, NY, T>{} : weights.Qy_terminal;
        Matrix<NP * NY, NZ, T>  QG{};
        Matrix<NP * NY, NXA, T> QPhi{};
        for (size_t kb = 0; kb < NP; ++kb) {
            const Matrix<NY, NY, T>& Qk = (kb == NP - 1) ? Q_last : weights.Qy;
            QG.template block<NY, NZ>(kb * NY, 0) = Qk * pred.Gamma.template block<NY, NZ>(kb * NY, 0).to_matrix();
            QPhi.template block<NY, NXA>(kb * NY, 0) = Qk * pred.Phi.template block<NY, NXA>(kb * NY, 0).to_matrix();
        }
        Hm = pred.Gamma.transpose() * QG;
        Fx = pred.Gamma.transpose() * QPhi;
        art.F_r = T{2} * QG.transpose();
        return true;
    }
    static_assert(NX == NXP, "The prediction model is the plant.");
    const auto W_opt = mpc_period_cost(plant, weights.Qy, art.scale_u, art.scale_y);
    const auto W_last_opt = dare_on ? W_opt : mpc_period_cost(plant, weights.Qy_terminal, art.scale_u, art.scale_y);
    if (!W_opt || !W_last_opt) {
        return false;
    }
    const auto&         W = *W_opt;
    const auto&         W_last = *W_last_opt;
    Matrix<NXA, NXA, T> Tx = Matrix<NXA, NXA, T>::identity(); // x̃_k = Tx·x̃₀ + TZ·Z
    Matrix<NXA, NZ, T>  TZ{};
    for (size_t k = 0; k < NP; ++k) {
        const Matrix<NXA, NXA, T> Tx1 = Aa * Tx;
        Matrix<NXA, NZ, T>        TZ1 = Aa * TZ;
        if (k < NC) {
            TZ1.template block<NXA, NU>(0, k * NU) = TZ1.template block<NXA, NU>(0, k * NU).to_matrix() + Ba;
        }
        //! q_k = Mq·x̃₀ + MZ·Z: the state and disturbance rows of x̃_k, the input rows of x̃_{k+1}
        Matrix<NXA, NXA, T> Mq = Tx;
        Matrix<NXA, NZ, T>  MZ = TZ;
        Mq.template block<NU, NXA>(NX, 0) = Tx1.template block<NU, NXA>(NX, 0);
        MZ.template block<NU, NZ>(NX, 0) = TZ1.template block<NU, NZ>(NX, 0);
        const auto&              Wk = (k == NP - 1) ? W_last : W;
        const Matrix<NZ, NXA, T> MZtW = MZ.transpose() * Wk.template block<NXA, NXA>(0, 0).to_matrix();
        Hm = Hm + (MZtW * MZ);
        Fx = Fx + (MZtW * Mq);
        art.F_r.template block<NZ, NY>(0, k * NY) = (MZ.transpose() * Wk.template block<NXA, NY>(0, NXA).to_matrix()) * T{-2};
        Tx = Tx1;
        TZ = TZ1;
    }
    return true;
}

/**
 * @brief The first step k at which each input's move reaches an output, or @p max_steps if it has
 *        not by then
 *
 * Output y's response to input j's move, k steps on, is (C A^k B)(y, j). Roundoff leaves a true zero
 * of it within about n·(k+1)·ε of (|C| |A|^k |B|)(y, j), so a response above 16 times that is a path.
 * The two scale together with the units of x, u and y, and a structural zero is exactly 0 in both,
 * so neither the units nor a dense basis of the same model changes the answer. Each input's column
 * is renormalized every step, so neither overflows on an unstable A.
 */
template<size_t NX, size_t NU, size_t NY, typename T>
[[nodiscard]] constexpr damp::array<size_t, NU> first_reaching_steps(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B, const Matrix<NY, NX, T>& C, size_t max_steps) {
    damp::array<size_t, NU> first{};
    for (size_t j = 0; j < NU; ++j) {
        ColVec<NX, T> v{}; // A^k B's column j, then scaled
        ColVec<NX, T> w{}; // |A|^k |B|'s, scaled the same
        for (size_t i = 0; i < NX; ++i) {
            v(i) = B(i, j);
            w(i) = damp::abs(B(i, j));
        }
        first[j] = max_steps;
        for (size_t k = 0; k < max_steps && first[j] == max_steps; ++k) {
            const T tol = T{16} * static_cast<T>(NX * (k + 1)) * eps<T>();
            for (size_t y = 0; y < NY; ++y) {
                T response{0};
                T bound{0};
                for (size_t i = 0; i < NX; ++i) {
                    response += C(y, i) * v(i);
                    bound += damp::abs(C(y, i)) * w(i);
                }
                first[j] = (damp::abs(response) > tol * bound) ? k : first[j];
            }
            ColVec<NX, T> v_next{};
            ColVec<NX, T> w_next{};
            T             scale{0};
            for (size_t r = 0; r < NX; ++r) {
                for (size_t c = 0; c < NX; ++c) {
                    v_next(r) += A(r, c) * v(c);
                    w_next(r) += damp::abs(A(r, c)) * w(c);
                }
                scale = damp::max(scale, w_next(r));
            }
            if (!(scale > T{0})) {
                break; // A^k B is exactly 0 from here on
            }
            v = v_next / scale;
            w = w_next / scale;
        }
    }
    return first;
}

/**
 * @brief Whether the first move is used: every input with no Ru weight reaches a predicted output
 *        within NP steps, and some input does
 *
 * Only the first move is applied, so an input whose first move no output sees is never used: a
 * horizon no longer than the input delay, or than the plant's relative degree, leaves it at rest.
 * An Ru weight still steers its input to the input target, as long as some input reaches an output.
 */
template<size_t NP, size_t NXA, size_t NU, size_t NY, typename T>
[[nodiscard]] constexpr bool mpc_first_moves_reach_outputs(
    const Matrix<NXA, NXA, T>& A,
    const Matrix<NXA, NU, T>&  B,
    const Matrix<NY, NXA, T>&  C,
    const Matrix<NU, NU, T>&   Ru
) {
    const auto first = first_reaching_steps(A, B, C, NP);
    bool       any = false;
    for (size_t j = 0; j < NU; ++j) {
        const bool reaches = first[j] < NP;
        any = any || reaches;
        if (!reaches && !(Ru(j, j) > T{0})) {
            return false;
        }
    }
    return any;
}

template<size_t NP, size_t NC, size_t ND, size_t NT_max = 0, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NXP>
[[nodiscard]] constexpr MPCArtifacts<NX, NU, NY, NP, NC, ND, T, NT_max> mpc_impl(
    const StateSpace<NX, NU, NY, T, NW, NV>&     sys,
    const Matrix<NX, ND, T>&                     Bd,
    const Matrix<NY, ND, T>&                     Cd,
    const MPCWeights<NU, NY, T>&                 weights,
    const MPCConstraints<NU, NY, T>&             constraints,
    const MpcContinuousCost<NXP, NU, NY, ND, T>& plant
) {
    static_assert(NP >= 1, "Prediction horizon must be at least one step.");
    static_assert(NC >= 1 && NC <= NP, "Control horizon must satisfy 1 <= NC <= NP.");

    constexpr size_t NXA = NX + NU + ND;
    constexpr size_t NZ = NC * NU;
    constexpr size_t NZ1 = NZ + 1;
    constexpr size_t NYS = NP * NY;
    constexpr size_t NI_box = (4 * NZ) + (2 * NYS);
    constexpr size_t NI = NI_box + NT_max;

    using Art = MPCArtifacts<NX, NU, NY, NP, NC, ND, T, NT_max>;
    Art art{};
    art.constraints = constraints;
    art.C = sys.C;
    art.Cd = Cd;
    art.Ts = sys.Ts;
    art.u_target = weights.u_target;
    art.terminal = constraints.terminal;
    art.max_qp_iterations = (constraints.max_qp_iterations == 0) ? 10 * (NZ1 + NI + 1) : constraints.max_qp_iterations;
    for (size_t i = 0; i < NT_max; ++i) {
        art.h_terminal(i) = unbounded_bound<T>();
    }

    if (constraints.terminal == MPCTerminal::invariant_set && NT_max == 0) {
        return ::damp::detail::set_status(art, DesignStatus::bad_weights);
    }
    if (!(sys.Ts > T{0}) || !damp::isfinite(sys.Ts)) {
        return ::damp::detail::set_status(art, DesignStatus::invalid_sample_time);
    }
    for (size_t r = 0; r < NY; ++r) {
        for (size_t c = 0; c < NU; ++c) {
            if (sys.D(r, c) != T{0}) {
                return ::damp::detail::set_status(art, DesignStatus::bad_weights);
            }
        }
    }
    if (!mpc_weights_valid(weights, constraints)) {
        return ::damp::detail::set_status(art, DesignStatus::bad_weights); // as design::lqr: indefinite, asymmetric or non-finite
    }
    if (!mpc_limits_valid(constraints.u_lim) || !mpc_limits_valid(constraints.du_lim) || !mpc_limits_valid(constraints.y_lim)) {
        return ::damp::detail::set_status(art, DesignStatus::invalid_parameter); // as set_limits: a NaN or inverted limit
    }
    if (constraints.terminal == MPCTerminal::invariant_set && NC < NP) {
        return ::damp::detail::set_status(art, DesignStatus::invalid_parameter); // the shifted plan needs a free move at every step
    }

    for (size_t j = 0; j < NU; ++j) {
        art.scale_u(j) = detail::scale_from_range(constraints.u_lim.lower(j), constraints.u_lim.upper(j));
    }
    for (size_t i = 0; i < NY; ++i) {
        art.scale_y(i) = detail::scale_from_range(constraints.y_lim.lower(i), constraints.y_lim.upper(i));
    }

    Matrix<NX, NU, T> Bs{};
    for (size_t i = 0; i < NX; ++i) {
        for (size_t j = 0; j < NU; ++j) {
            Bs(i, j) = sys.B(i, j) * art.scale_u(j);
        }
    }
    Matrix<NY, NX, T> Cs{};
    for (size_t i = 0; i < NY; ++i) {
        for (size_t j = 0; j < NX; ++j) {
            Cs(i, j) = sys.C(i, j) / art.scale_y(i);
        }
    }

    Matrix<NXA, NXA, T> Aa{};
    Aa.template block<NX, NX>(0, 0) = sys.A;
    Aa.template block<NX, NU>(0, NX) = Bs;
    Aa.template block<NU, NU>(NX, NX) = Matrix<NU, NU, T>::identity();
    if constexpr (ND > 0) {
        Aa.template block<NX, ND>(0, NX + NU) = Bd;
        Aa.template block<ND, ND>(NX + NU, NX + NU) = Matrix<ND, ND, T>::identity();
    }
    Matrix<NXA, NU, T> Ba{};
    Ba.template block<NX, NU>(0, 0) = Bs;
    Ba.template block<NU, NU>(NX, 0) = Matrix<NU, NU, T>::identity();
    Matrix<NY, NXA, T> Ca{};
    Ca.template block<NY, NX>(0, 0) = Cs;
    if constexpr (ND > 0) {
        for (size_t i = 0; i < NY; ++i) { // the outputs carry an estimated output bias, so the y limits see it
            for (size_t j = 0; j < ND; ++j) {
                Ca(i, NX + NU + j) = Cd(i, j) / art.scale_y(i);
            }
        }
    }

    const auto pred = mpc_predictions<NP, NC, NXA, NX, NU, NY, T>(Aa, Ba, Ca);
    //! A terminal cost weights the state at NP, delay registers included, so it sees past the horizon.
    if (constraints.terminal == MPCTerminal::none && !mpc_first_moves_reach_outputs<NP>(Aa, Ba, Ca, weights.Ru)) {
        return ::damp::detail::set_status(art, DesignStatus::invalid_parameter);
    }

    //! On a continuous plant the move and input weights are per second: R_Δu on the rate Δu/Ts, R_u over the period
    const bool              dare_on = constraints.terminal != MPCTerminal::none;
    const Matrix<NU, NU, T> Rdu = plant.on ? Matrix<NU, NU, T>(weights.Rdu * (T{1} / plant.Ts)) : weights.Rdu;
    const Matrix<NU, NU, T> Ru = plant.on ? Matrix<NU, NU, T>(weights.Ru * plant.Ts) : weights.Ru;
    Matrix<NZ, NZ, T>       Hm{};
    Matrix<NZ, NXA, T>      Fx{};
    if (!mpc_output_cost(art, pred, Aa, Ba, weights, plant, dare_on, Hm, Fx)) {
        return ::damp::detail::set_status(art, DesignStatus::bad_weights); // the plant was sampled: a non-finite weight
    }
    for (size_t i = 0; i < NC; ++i) {
        for (size_t j = 0; j < NC; ++j) {
            const size_t m = (i > j) ? i : j;
            Hm.template block<NU, NU>(i * NU, j * NU) += static_cast<T>(NC - m) * Ru;
        }
        Hm.template block<NU, NU>(i * NU, i * NU) += Rdu;
    }
    art.H.template block<NZ, NZ>(0, 0) = T{2} * Hm;
    art.H(NZ, NZ) = T{2} * weights.ecr_weight;

    for (size_t i = 0; i < NC; ++i) {
        Fx.template block<NU, NU>(i * NU, NX) += static_cast<T>(NC - i) * Ru;
        for (size_t j = i; j < NC; ++j) {
            art.F_ut.template block<NU, NU>(i * NU, j * NU) = T{2} * Ru;
        }
    }
    art.F_x = T{2} * Fx;

    //! z_N = [x_N − x_t; u′_{N−1} − u′_t] = Φ_z·x̃ + Γ_z·Z: the input rows select u′_prev and sum the moves.
    Matrix<NX + NU, NXA, T> Phi_z{};
    Matrix<NX + NU, NZ, T>  Gamma_z{};
    Phi_z.template block<NX, NXA>(0, 0) = pred.Phi_x;
    Phi_z.template block<NU, NU>(NX, NX) = Matrix<NU, NU, T>::identity();
    Gamma_z.template block<NX, NZ>(0, 0) = pred.Gamma_x;
    for (size_t j = 0; j < NC; ++j) {
        Gamma_z.template block<NU, NU>(NX, j * NU) = Matrix<NU, NU, T>::identity();
    }
    if (dare_on) {
        const auto stage = plant.on ? mpc_stage_of<NX, NU>(*mpc_period_cost(plant, weights.Qy, art.scale_u, art.scale_y)) // formed above
                                    : MpcStageCost<NX, NU, T>{.Q_x = Cs.transpose() * (weights.Qy * Cs)};
        const auto st = mpc_dare_lqr(sys.A, Bs, stage, Rdu, Ru, art.P, art.K_terminal);
        if (st != DesignStatus::ok) {
            return ::damp::detail::set_status(art, st);
        }
        const auto HP = Gamma_z.transpose() * (art.P * Gamma_z);
        const auto FP = Gamma_z.transpose() * (art.P * Phi_z);
        art.H.template block<NZ, NZ>(0, 0) = art.H.template block<NZ, NZ>(0, 0).to_matrix() + (T{2} * HP);
        art.F_x = art.F_x + (T{2} * FP);
    }

    constexpr size_t r_du_hi = 0;
    constexpr size_t r_du_lo = NZ;
    constexpr size_t r_u_hi = 2 * NZ;
    constexpr size_t r_u_lo = 3 * NZ;
    constexpr size_t r_y_hi = 4 * NZ;
    constexpr size_t r_y_lo = (4 * NZ) + NYS;

    const auto I_u = Matrix<NU, NU, T>::identity();
    for (size_t i = 0; i < NC; ++i) {
        art.A_con.template block<NU, NU>(r_du_hi + (i * NU), i * NU) = I_u;
        art.A_con.template block<NU, NU>(r_du_lo + (i * NU), i * NU) = -I_u;
        for (size_t j = 0; j <= i; ++j) {
            art.A_con.template block<NU, NU>(r_u_hi + (i * NU), j * NU) = I_u;
            art.A_con.template block<NU, NU>(r_u_lo + (i * NU), j * NU) = -I_u;
        }
        art.B_map.template block<NU, NU>(r_u_hi + (i * NU), NX) = -I_u;
        art.B_map.template block<NU, NU>(r_u_lo + (i * NU), NX) = I_u;
    }
    art.A_con.template block<NYS, NZ>(r_y_hi, 0) = pred.Gamma;
    art.A_con.template block<NYS, NZ>(r_y_lo, 0) = -pred.Gamma;
    art.B_map.template block<NYS, NXA>(r_y_hi, 0) = -pred.Phi;
    art.B_map.template block<NYS, NXA>(r_y_lo, 0) = pred.Phi;
    for (size_t kb = 0; kb < NP; ++kb) {
        for (size_t j = 0; j < NY; ++j) {
            art.A_con(r_y_hi + (kb * NY) + j, NZ) = -constraints.y_upper_ecr(j);
            art.A_con(r_y_lo + (kb * NY) + j, NZ) = -constraints.y_lower_ecr(j);
        }
    }
    detail::mpc_box_bounds(art, MPCLimits<NU, NY, T>{constraints.u_lim, constraints.du_lim, constraints.y_lim});

    if constexpr (NT_max > 0) {
        if (constraints.terminal == MPCTerminal::invariant_set) {
            if (!(constraints.terminal_target_margin > T{0}) || !(constraints.terminal_target_margin < T{1})) {
                return ::damp::detail::set_status(art, DesignStatus::invalid_parameter);
            }
            //! The set of w = [z; q] under w⁺ = blkdiag(A_z − B_z·K, I)·w: valid for every target
            //! q = [y_t; u_t] in the tightened boxes, the target entering each row's offset.
            constexpr size_t  NZD = NX + NU;
            constexpr size_t  NA = NZD + NY + NU;
            constexpr size_t  NG = (4 * NY) + (6 * NU);
            const auto        tm = mpc_terminal_model(sys.A, Bs);
            Matrix<NA, NA, T> A_w{};
            A_w.template block<NZD, NZD>(0, 0) = tm.A - (tm.B * art.K_terminal);
            A_w.template block<NY + NU, NY + NU>(NZD, NZD) = Matrix<NY + NU, NY + NU, T>::identity();
            Matrix<NG, NA, T> G0{};
            ColVec<NG, T>     h0{};
            mpc_terminal_rows<NX, NU, NY, T>(sys.C, art.K_terminal, art.scale_u, art.scale_y, constraints, G0, h0);
            const auto inv = maximal_invariant_set<NT_max>(A_w, G0, h0);
            if (!inv.success) {
                return ::damp::detail::set_status(art, inv.status);
            }
            //! Keep the rows that constrain z; the rest bound the target alone (a design assumption).
            size_t kept = 0;
            for (size_t r = 0; r < inv.n_rows; ++r) {
                T zsum{0};
                for (size_t j = 0; j < NZD; ++j) {
                    zsum += damp::abs(inv.G(r, j));
                }
                if (!(zsum > T{0})) {
                    continue;
                }
                for (size_t j = 0; j < NZD; ++j) {
                    art.G_terminal(kept, j) = inv.G(r, j);
                }
                for (size_t j = 0; j < NY + NU; ++j) {
                    art.G_terminal_target(kept, j) = inv.G(r, NZD + j);
                }
                art.h_terminal(kept) = inv.h(r);
                ++kept;
            }
            art.n_terminal = kept;
        }
    }
    mpc_stack_terminal_rows(art, Phi_z, Gamma_z);
    art.A_con(NI, NZ) = T{-1};
    art.b_bound(NI) = T{0};

    if (!mat::cholesky(art.H)) {
        return ::damp::detail::set_status(art, DesignStatus::singular);
    }
    return ::damp::detail::set_status(art, DesignStatus::ok);
}

} // namespace detail

/**
 * @brief A measured disturbance: a known input the controller does not command, entering the
 *        state update through @c Bd (a load torque, a feed rate, the ambient temperature).
 *
 * Given to design::mpc on the measured state, it adds the channel to the prediction model, and
 * the runtime takes the disturbance each tick (MPC::step(r, d, x)). It is held over the
 * horizon and not rescaled, and the steady-state target is the one under it, so a terminal
 * cost or set holds about the steady state d moves the plant to. B_d is in the plant's time domain: continuous on a continuous
 * plant (discretized with it by ZOH), per step on a discrete one.
 *
 * @note Compare with MATLAB®'s measured-disturbance channels (setmpcsignals 'MD').
 * @tparam NX Plant states; NM Measured-disturbance inputs
 */
template<size_t NX, size_t NM, typename T = double>
struct MeasuredDisturbance {
    Matrix<NX, NM, T> Bd{}; ///< How the measured disturbance enters the state update
};

/** @brief Sizes from the matrix: `MeasuredDisturbance{Bd}`. */
template<size_t NX, size_t NM, typename T>
MeasuredDisturbance(Matrix<NX, NM, T>) -> MeasuredDisturbance<NX, NM, T>;

/**
 * @struct MPCResult
 * @brief An MPC design: the condensed QP, the estimator when there is one, and the steady-state
 *        target the runtime poses the QP about; consumed by damp::MPC
 *
 * @tparam Law       The QP on the law's model, MPCArtifacts
 * @tparam Estimator design::kalman's result with the design's options, or NoEstimator when the
 *                   measurement is the state
 * @tparam Target    SteadyStateTargetResult on the estimator's model (the law's without one)
 */
template<typename Law, typename Estimator, typename Target>
struct [[nodiscard]] MPCResult {
    Law                             mpc{};                               ///< The Δu QP on the law's model
    [[no_unique_address]] Estimator kalman{};                            ///< The filter's design (nothing when the measurement is the state)
    Target                          target{};                            ///< The maps from the reference and the disturbance to (x_t, u_t)
    bool                            success{false};                      ///< status == DesignStatus::ok
    DesignStatus                    status{DesignStatus::not_converged}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The design in another scalar; refused when the QP's Hessian does not factor in it. */
    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        using LawU = decltype(mpc.template as<U>());
        using EstimatorU = decltype(kalman.template as<U>());
        using TargetU = decltype(target.template as<U>());
        const LawU         law = mpc.template as<U>();
        const DesignStatus st = (status == DesignStatus::ok) ? law.status : status;
        return MPCResult<LawU, EstimatorU, TargetU>{law, kalman.template as<U>(), target.template as<U>(), st == DesignStatus::ok, st};
    }
};

namespace detail {

/** @brief Not the limits of an MPC. */
template<typename>
struct is_mpc_constraints : std::false_type {};

/** @brief MPCConstraints. */
template<size_t NU, size_t NY, typename T>
struct is_mpc_constraints<MPCConstraints<NU, NY, T>> : std::true_type {};

/** @brief Not a measured disturbance. */
template<typename>
struct is_measured_disturbance : std::false_type {};

/** @brief MeasuredDisturbance. */
template<size_t NX, size_t NM, typename T>
struct is_measured_disturbance<MeasuredDisturbance<NX, NM, T>> : std::true_type {};

/** @brief The options design::mpc takes after the weights. */
using mpc_options = ::damp::detail::option_set<
    is_mpc_constraints, is_kalman_noise, is_disturbance, is_measured_disturbance>;

/** @brief Whether a design::mpc call has an estimator: a KalmanNoise among its options. */
template<typename... Opts>
inline constexpr bool mpc_estimated = ::damp::detail::option_count<is_kalman_noise, Opts...> > 0;

/** @brief The controlled outputs of a design::mpc call: those its disturbance spec lists, else every output. */
template<size_t NU, size_t NY, size_t NW, typename... Opts>
inline constexpr size_t mpc_outputs = (disturbance_traits<disturbance_option_t<Opts...>, NU, NY, NW>::listed == 0)
                                        ? NY
                                        : disturbance_traits<disturbance_option_t<Opts...>, NU, NY, NW>::listed;

/** @brief The first status that is not ok among @p statuses, else ok. */
template<typename... S>
[[nodiscard]] constexpr DesignStatus first_failure(S... statuses) {
    DesignStatus out = DesignStatus::ok;
    ((out = (out == DesignStatus::ok) ? statuses : out), ...);
    return out;
}

/**
 * @brief The nominal steady-state target maps: (x_t, u_t) for a reference r and a measured
 *        disturbance d (@p Bd), the weighted least-squares steady state of the controlled
 *        outputs among the model's equilibria (detail::equilibrium_target, with @p Q the output
 *        weight in physical units)
 *
 * Always an equilibrium of the model under d, so the QP posed about it is the absolute one and
 * the terminal ingredients hold about the steady state the measured d moves the plant to.
 */
template<size_t NZ, size_t ND, size_t NXE, size_t NU, size_t NY, typename T, size_t NWE, size_t NVE, size_t NXL>
[[nodiscard]] constexpr SteadyStateTargetResult<NXE, NU, ND, NZ, T> nominal_target(
    const StateSpace<NXE, NU, NY, T, NWE, NVE>& model,
    const Matrix<NZ, NXL, T>&                   Cz,
    const Matrix<NZ, NU, T>&                    Dz,
    const Matrix<NZ, NZ, T>&                    Q,
    const Matrix<NXE, ND, T>&                   Bd
) {
    static_assert(NXL <= NXE, "The controlled outputs read the leading states of the model.");
    constexpr size_t   NS = NXE + NU;
    Matrix<NXE, NS, T> E{};
    E.template block<NXE, NXE>(0, 0) = model.A - Matrix<NXE, NXE, T>::identity();
    E.template block<NXE, NU>(0, NXE) = model.B;
    Matrix<NZ, NS, T> S{};
    S.template block<NZ, NXL>(0, 0) = Cz;
    S.template block<NZ, NU>(0, NXE) = Dz;
    return equilibrium_target<NXE, NU>(E, S, Bd, Matrix<NZ, ND, T>{}, Q);
}

/**
 * @brief The output weight of the cost in physical units: S_y⁻¹·Q_y·S_y⁻¹ with the QP's output
 *        scales of a design that succeeded (@p art), else the identity.
 */
template<size_t NZ, typename T, typename Art>
[[nodiscard]] constexpr Matrix<NZ, NZ, T> physical_output_weight(const Matrix<NZ, NZ, T>& Qy, const Art& art) {
    if (!art.success) {
        return Matrix<NZ, NZ, T>::identity(); // a refused design has no scales; its target is not used
    }
    Matrix<NZ, NZ, T> Q{};
    for (size_t i = 0; i < NZ; ++i) {
        for (size_t j = 0; j < NZ; ++j) {
            Q(i, j) = Qy(i, j) / (art.scale_y(i) * art.scale_y(j));
        }
    }
    return Q;
}

/**
 * @brief The continuous stage cost of a design from @p sys at @p Ts: the plant, the controlled
 *        outputs @p rows of C and D and the disturbance's continuous channel @p Bd_c; off on a
 *        discrete plant (@p Discrete).
 */
template<bool Discrete, size_t NZ, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t ND>
[[nodiscard]] constexpr MpcContinuousCost<NX, NU, NZ, ND, T> mpc_continuous_cost(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const damp::array<size_t, NZ>&           rows,
    const Matrix<NX, ND, T>&                 Bd_c,
    T                                        Ts
) {
    MpcContinuousCost<NX, NU, NZ, ND, T> out{};
    if constexpr (!Discrete) {
        out = MpcContinuousCost<NX, NU, NZ, ND, T>{.A = sys.A, .B = sys.B, .Bd = Bd_c, .Ts = Ts, .on = true};
        for (size_t i = 0; i < NZ; ++i) {
            for (size_t j = 0; j < NX; ++j) {
                out.C(i, j) = sys.C(rows[i], j);
            }
            for (size_t j = 0; j < NU; ++j) {
                out.D(i, j) = sys.D(rows[i], j);
            }
        }
    } else {
        static_cast<void>(sys);
        static_cast<void>(rows);
        static_cast<void>(Bd_c);
        static_cast<void>(Ts);
    }
    return out;
}

/**
 * @brief The MPC on the measured state: the measurement is x, the controlled outputs C·x (D = 0),
 *        the QP on the plant's model, the
 *        measured-disturbance channel, and the nominal target.
 */
template<bool Discrete, size_t NP, size_t NC, size_t NT_max, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NXD, size_t NM>
[[nodiscard]] constexpr auto mpc_measured(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        Ts,
    const MPCWeights<NU, NY, T>&             weights,
    const MPCConstraints<NU, NY, T>&         constraints,
    const MeasuredDisturbance<NXD, NM, T>&   measured
) {
    static_assert(NM == 0 || NXD == NX, "The measured disturbance enters the plant's states: MeasuredDisturbance<NX, NM>.");
    constexpr size_t NXL = NX;
    using Law = MPCArtifacts<NXL, NU, NY, NP, NC, NM, T, NT_max>;
    using Target = SteadyStateTargetResult<NXL, NU, NM, NY, T>;
    using Result = MPCResult<Law, NoEstimator, Target>;
    const auto m = sampled_model(sys, Ts, Discrete);
    if (m.status != DesignStatus::ok) {
        return ::damp::detail::refused<Result>(m.status);
    }
    Matrix<NXL, NM, T> Bd{};
    if constexpr (NM > 0) {
        if constexpr (Discrete) {
            Bd.template block<NX, NM>(0, 0) = measured.Bd;
        } else {
            const auto d = ::damp::discretize(StateSpace<NX, NM, 1, T>{.A = sys.A, .B = measured.Bd}, Ts);
            if (!d) {
                return ::damp::detail::refused<Result>(d.status); // a non-finite B_d (bad_weights)
            }
            Bd.template block<NX, NM>(0, 0) = d->B;
        }
    } else {
        static_cast<void>(measured);
    }
    Result     out{};
    const auto cost = mpc_continuous_cost<Discrete>(sys, *selected_outputs<NY>(damp::array<size_t, 0>{}), measured.Bd, Ts);
    out.mpc = mpc_impl<NP, NC, NM, NT_max>(m.law, Bd, Matrix<NY, NM, T>{}, weights, constraints, cost);
    out.target = nominal_target(m.law, m.law.C, m.law.D, physical_output_weight(weights.Qy, out.mpc), Bd);
    return ::damp::detail::set_status(out, out.mpc.status);
}

/**
 * @brief The MPC through an estimator: the filter of design::kalman with the same options, the
 *        target (design::steady_state_target with a disturbance, the nominal one without), and
 *        the QP on the law's model with the controlled outputs and the disturbance channel.
 */
template<
    bool Discrete, size_t NP, size_t NC, size_t NT_max, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NZ, size_t NWN, size_t NVN,
    size_t NDN, typename Spec>
[[nodiscard]] constexpr auto mpc_estimated_design(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        Ts,
    const MPCWeights<NU, NZ, T>&             weights,
    const MPCConstraints<NU, NZ, T>&         constraints,
    const KalmanNoise<NWN, NVN, T, NDN>&     noise,
    const Spec&                              spec
) {
    constexpr size_t ND = disturbance_traits<Spec, NU, NY, NW>::states;
    static_assert(NWN == NW && NVN == NV, "KalmanNoise is sized to the plant's noise inputs: KalmanNoise<NW, NV> for the columns of G and H.");
    static_assert(NDN == ND, "KalmanNoise carries the random walk of the spec's ND disturbances: KalmanNoise<NW, NV, T, ND>.");
    constexpr size_t NXL = NX;
    constexpr size_t NXE = NXL;
    const auto       kf = kalman_design<Discrete>(sys, Ts, noise, spec);
    using Target = SteadyStateTargetResult<NXE, NU, ND, NZ, T>;
    using Law = MPCArtifacts<NXL, NU, NZ, NP, NC, ND, T, NT_max>;
    using Result = MPCResult<Law, std::remove_cvref_t<decltype(kf)>, Target>;
    Result     out{.kalman = kf};
    const auto m = sampled_model(sys, Ts, Discrete);
    const auto rows = selected_outputs<NY>(listed_outputs(spec));
    if (m.status != DesignStatus::ok || !rows) {
        return ::damp::detail::set_status(out, first_failure(kf.status, m.status, DesignStatus::bad_weights)); // a delay refused, or a controlled output out of range or listed twice
    }
    //! The QP predicts and constrains the controlled outputs z, the listed rows of the law's model
    StateSpace<NXL, NU, NZ, T> law{.A = m.law.A, .B = m.law.B, .Ts = m.law.Ts};
    for (size_t i = 0; i < NZ; ++i) {
        for (size_t j = 0; j < NXL; ++j) {
            law.C(i, j) = m.law.C((*rows)[i], j);
        }
        for (size_t j = 0; j < NU; ++j) {
            law.D(i, j) = m.law.D((*rows)[i], j);
        }
    }
    const auto         channels = disturbance_channels(spec, sys, m.law, Ts);
    const auto         cost = mpc_continuous_cost<Discrete>(sys, *rows, channels.Bd_c, Ts);
    Matrix<NXL, ND, T> Bd{};
    if constexpr (ND > 0) {
        //! Offset-free: the target under the estimate, which rides the QP's disturbance channel
        Bd = channels.Bd.template block<NXL, ND>(0, 0);
        Matrix<NZ, ND, T> Cd{}; // the controlled rows of C_d: the target holds z + C_{d,z}·d̂ at r
        for (size_t i = 0; i < NZ; ++i) {
            for (size_t j = 0; j < ND; ++j) {
                Cd(i, j) = channels.Cd((*rows)[i], j);
            }
        }
        out.mpc = mpc_impl<NP, NC, ND, NT_max>(law, Bd, Cd, weights, constraints, cost);
        out.target = target_design<Discrete>(sys, Ts, spec, physical_output_weight(weights.Qy, out.mpc));
    } else {
        out.mpc = mpc_impl<NP, NC, ND, NT_max>(law, Bd, Matrix<NZ, ND, T>{}, weights, constraints, cost);
        out.target = nominal_target(m.law, law.C, law.D, physical_output_weight(weights.Qy, out.mpc), Bd);
    }
    return ::damp::detail::set_status(out, first_failure(kf.status, out.target.status, out.mpc.status));
}

/** @brief design::mpc in the plant's time domain (@p Discrete) or sampled at @p Ts: the options sorted by kind. */
template<bool Discrete, size_t NP, size_t NC, size_t NT_max, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NZ, typename... Opts>
[[nodiscard]] constexpr auto mpc_design(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T Ts, const MPCWeights<NU, NZ, T>& weights, const Opts&... opts) {
    mpc_options::check<Opts...>();
    using Limits = ::damp::detail::option_t<is_mpc_constraints, MPCConstraints<NU, NZ, T>, Opts...>;
    static_assert(std::is_same_v<Limits, MPCConstraints<NU, NZ, T>>, "The limits are on the inputs and the controlled outputs: MPCConstraints<NU, NZ>.");
    const auto constraints = ::damp::detail::option<is_mpc_constraints>(MPCConstraints<NU, NZ, T>{}, opts...);
    if constexpr (mpc_estimated<Opts...>) {
        static_assert(
            ::damp::detail::option_count<is_measured_disturbance, Opts...> == 0,
            "A measured disturbance with an estimator is not supported: the filter does not take it as an input."
        );
        return mpc_estimated_design<Discrete, NP, NC, NT_max>(
            sys, Ts, weights, constraints, ::damp::detail::option<is_kalman_noise>(KalmanNoise<NW, NV, T>{}, opts...),
            ::damp::detail::option<is_disturbance>(NoDisturbance{}, opts...)
        );
    } else {
        static_assert(::damp::detail::option_count<is_disturbance, Opts...> == 0, "A disturbance spec needs KalmanNoise: the disturbance is estimated.");
        return mpc_measured<Discrete, NP, NC, NT_max>(
            sys, Ts, weights, constraints, ::damp::detail::option<is_measured_disturbance>(MeasuredDisturbance<NX, 0, T>{}, opts...)
        );
    }
}

} // namespace detail

/**
 * @brief Synthesize a constrained linear MPC (condensed dense QP, Δu form) on a discrete plant
 *
 * Builds the prediction matrices for the augmented system
 * @f[
 *   \tilde A = \begin{bmatrix} A & B' & B_d \\ 0 & I & 0 \\ 0 & 0 & I \end{bmatrix},\quad
 *   \tilde B = \begin{bmatrix} B' \\ I \\ 0 \end{bmatrix},\quad
 *   \tilde C = [\,C' \;\; 0 \;\; 0\,]
 * @f]
 * (primes: internally scaled B/C), condenses the horizon cost into
 * @f$ \tfrac{1}{2} \zeta^\top H \zeta + f^\top \zeta @f$ over
 * @f$ \zeta = [Z; \varepsilon] @f$ with
 * @f$ f = [F_x \tilde x_0 - F_r r' - F_{ut} u'_{target};\; 0] @f$, and stacks
 * the hard u/Δu and soft y box constraints as
 * @f$ A_{con} \zeta \le b_{bound} + B_{map}\tilde x_0 @f$.
 * Beyond NC−1 the input is held (Δu = 0), which the augmented model carries automatically.
 *
 * | Design | Runtime |
 * | --- | --- |
 * | mpc<NP, NC>(plant, [Ts,] weights, constraints) | MPC: the state measured |
 * | mpc<NP, NC>(plant, [Ts,] weights, constraints, MeasuredDisturbance{Bd}) | MPC with a known disturbance input, given each tick |
 * | ofmpc<NP, NC>(plant, [Ts,] weights, constraints, noise, disturbance) | OFMPC: a Kalman filter on [x; d], offset-free (below) |
 *
 * The constraints are the box limits on u, Δu and the controlled outputs, the terminal
 * ingredients and the QP budget (`{}`: unbounded, none).
 *
 * Measured state. design::mpc's runtime takes the plant state x as its measurement
 * (update(x)), and the plant's C still defines the controlled outputs y = C·x that r, Qy and
 * the y bounds refer to (NY = rows of C), as the LQI on a measured state integrates r − C·x.
 * For design::ofmpc, C and D are the measured outputs, the filter reconstructs the state from
 * them, and the QP is on the controlled outputs, the rows of C the disturbance spec lists (all
 * of them without one).
 *
 * The target. The runtime poses the QP about the steady state (x_t, u_t) of the reference:
 * under the estimated disturbance with a disturbance spec (below), otherwise the nominal one,
 * the steady state among the model's equilibria closest to r in the Qy-weighted norm (an
 * output without cost does not pull it). An equilibrium leaves the QP the absolute one; the
 * terminal cost and set are on the deviation from it.
 *
 * Input delay. A delay is a model of the plant: design on discretize(plant, Ts, InputDelay<N>{θ}),
 * whose register states the QP's model carries like any other. The first outputs of the horizon
 * are ones no decision can change, so a longer delay wants a longer NP.
 *
 * Horizon. Without a terminal cost, each input's first move has to reach a predicted output
 * within NP steps, or the law would never use that input: a shorter horizon, NP = 1 on a model
 * with a one-sample delay for one, is refused (invalid_parameter). An input with an Ru weight is
 * the exception, since Ru steers it to its input target, as long as some input reaches an output.
 * The test is relative to the data: the units of x, u and y and the model's basis do not change it.
 * suggest_mpc_horizon(sys) counts the dead time into the horizon it suggests.
 *
 * ## Offset-free control
 *
 * design::ofmpc is the offset-free MPC of Muske and
 * Badgwell and of Pannocchia and Rawlings. The plant model is augmented with ND integrating
 * disturbances,
 * @f[
 *   x_{k+1} = A x_k + B u_k + B_d d_k,\qquad d_{k+1} = d_k,\qquad
 *   y_k = C x_k + C_d d_k + D u_k ,
 * @f]
 * where the spec chooses where d enters: InputDisturbance (@f$ B_d = B,\ C_d = D @f$, one per
 * input: a load or an actuator bias, the usual choice, detectable on an integrating plant),
 * OutputDisturbance (@f$ B_d = 0,\ C_d = I @f$, a bias on each reading, MATLAB®'s default),
 * PlantDisturbance (@f$ B_d = G,\ C_d = H @f$) or an explicit DisturbanceModel{.Bd, .Cd} (a
 * force on a load that an accelerometer also feels).
 *
 * The Kalman filter estimates @f$ [\hat x;\, \hat d] @f$ from the readings (design::kalman with
 * the same options; @c noise.disturbance is the random walk of d, how fast the estimate may
 * move against the measurement noise). The estimate rides the QP's disturbance channels
 * @f$ B_d @f$ and @f$ C_{d,z} @f$, so the predictions and the output limits carry it. Each
 * tick the steady-state target holds the controlled outputs @f$ z = C_z x + C_{d,z} \hat d @f$
 * at r under the estimate,
 * @f[
 *   \begin{bmatrix} A-I & B \\ C_z & 0 \end{bmatrix}
 *   \begin{bmatrix} x_t \\ u_t \end{bmatrix}
 *   =
 *   \begin{bmatrix} -B_d\hat d \\ r - C_{d,z}\hat d \end{bmatrix},
 * @f]
 * solved once at design for unit d and r (design::steady_state_target), and the QP is posed in
 * deviations from (x_t, u_t). At a steady state of the closed loop the innovation is zero, so
 * the readings match the model with @f$ \hat d @f$ and z = r whatever the true disturbance or
 * model error: no offset under constant disturbances and plant/model mismatch.
 *
 * The spec's @c outputs lists the rows of the measured outputs held offset-free, z (all of them
 * when it lists none); the weights, the y limits and the reference are on z. Offset-free
 * tracking of every controlled output needs NZ ≤ NU and a steady state that reaches any r;
 * with NZ = NU the target is unique or the design is refused (DesignStatus::singular), and with
 * more controlled outputs than inputs it is the least-squares compromise. The augment must be
 * detectable, rank [A − I, B_d; C, C_d] = NX + ND (so ND ≤ NY), or the design is refused with
 * DesignStatus::not_detectable: an output disturbance on a plant with a pole at z = 1 is not,
 * which is why MATLAB® falls back to input disturbances on integrating plants.
 *
 * The estimator holds deviations from the target (Luenberger with the fixed gain of
 * design::kalman, anchored by SteadyStateTarget each tick), so the QP's gradient stays a sum of
 * small terms however far the plant sits from zero (a temperature at 200 °C, a position in
 * metres at micrometre precision).
 *
 * Fails (success = false) if the plant is not discrete (DesignStatus::invalid_sample_time; use
 * the overload with the sample time), the controlled outputs have feedthrough, the slack
 * penalty is non-positive, the Hessian is not positive definite, or the filter or the target is
 * refused.
 *
 * @note Compare with MATLAB®'s mpc(plant, Ts, p, m) — p = NP, m = NC — with its built-in
 *       estimator (setindist / setoutdist, getEstimator) and mpcmove(mpcobj, xc, ym, r).
 * @see MPCWeights, MPCConstraints, damp::MPC, design::kalman, design::steady_state_target
 * @see Maciejowski, "Predictive Control with Constraints," 2002, ch. 2–3
 * @see "Disturbance modeling for offset-free linear model predictive control" (Muske &
 *      Badgwell, 2002), Journal of Process Control 12(5)
 * @see "Disturbance Models for Offset-Free Model-Predictive Control" (Pannocchia & Rawlings,
 *      2003), AIChE Journal 49(2)
 * @see "Offset-free reference tracking with model predictive control" (Maeder, Borrelli &
 *      Morari, 2009), Automatica 45(9)
 *
 * @tparam NP     Prediction horizon (steps)
 * @tparam NC     Control horizon (moves, NC ≤ NP)
 * @tparam NT_max Terminal-set row capacity (default 0; required > 0 for invariant_set)
 * @param sys      Discrete-time plant (Ts > 0); C the controlled outputs, D = 0
 * @param weights  Cost weights on the controlled outputs, per prediction step (default: Qy = I, RΔu = I)
 * @param constraints The limits ({}: unbounded)
 * @return MPCResult: deploy with `MPC law{res}`
 */
template<size_t NP, size_t NC, size_t NT_max = 0, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr auto mpc(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const MPCWeights<NU, NY, T>&             weights,
    const MPCConstraints<NU, NY, T>&         constraints
) {
    return detail::mpc_design<true, NP, NC, NT_max>(sys, T{0}, weights, constraints);
}

/** @brief design::mpc on a discrete plant with a measured disturbance input @p md (Bd), given each tick: runs as MPC. */
template<size_t NP, size_t NC, size_t NT_max = 0, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NM>
[[nodiscard]] constexpr auto mpc(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const MPCWeights<NU, NY, T>&             weights,
    const MPCConstraints<NU, NY, T>&         constraints,
    const MeasuredDisturbance<NX, NM, T>&    md
) {
    return detail::mpc_design<true, NP, NC, NT_max>(sys, T{0}, weights, constraints, md);
}

/**
 * @brief Synthesize a constrained linear MPC from the continuous plant and its continuous cost: the QP on its ZOH model at @p Ts
 *
 * The QP on design::mpc(discretize(plant, Ts), …)'s model, the options as there, with the cost
 * in continuous time: the weights are per second, like those of the sampled-data LQ designs
 * (lqr(plant, Ts, …)). Qy weights the output error integrated over each period, exactly under
 * the held input (Van Loan, with every cross term between the state, the input, the
 * disturbance and the reference), Qy_terminal the last period's, Ru the input's deviation
 * integrated over the period, and RΔu the input rate Δu/Ts: a move costs Δuᵀ·RΔu·Δu / Ts. A
 * terminal cost (MPCTerminal) is the cost-to-go of that period cost. The same weights at a
 * shorter Ts, with the horizons in steps rescaled to the same seconds, give nearly the same
 * loop. The y bounds hold at the samples. With KalmanNoise the estimator is
 * continuous in, discrete out (design::kalman(plant, Ts, noise, …)): the process and
 * disturbance noise are densities, discretized exactly over the period (Van Loan), and the
 * measurement noise is the variance of one reading.
 *
 * Refused with DesignStatus::invalid_sample_time when @p plant is already discrete or @p Ts is
 * not finite and positive, and as the discrete form is otherwise.
 *
 * @tparam NP     Prediction horizon (steps)
 * @tparam NC     Control horizon (moves, NC ≤ NP)
 * @tparam NT_max Terminal-set row capacity (default 0; required > 0 for invariant_set)
 * @see "Computing Integrals Involving the Matrix Exponential" (Van Loan, 1978), IEEE TAC 23(3)
 * @see lqr(plant, Ts, …) for the same discretization of an LQ cost
 *
 * @param plant    Continuous-time plant (Ts = 0)
 * @param Ts       Control period [s]
 * @param weights  Cost weights on the controlled outputs, per second
 * @param constraints The limits ({}: unbounded)
 * @return MPCResult: deploy with `MPC law{res}`; success = false if refused or the plant does not discretize
 */
template<size_t NP, size_t NC, size_t NT_max = 0, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr auto mpc(
    const StateSpace<NX, NU, NY, T, NW, NV>& plant,
    std::type_identity_t<T>                  Ts,
    const MPCWeights<NU, NY, T>&             weights,
    const MPCConstraints<NU, NY, T>&         constraints
) {
    return detail::mpc_design<false, NP, NC, NT_max>(plant, Ts, weights, constraints);
}

/** @brief The sampled-data design::mpc with a measured disturbance input @p md (Bd continuous): runs as MPC. */
template<size_t NP, size_t NC, size_t NT_max = 0, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NM>
[[nodiscard]] constexpr auto mpc(
    const StateSpace<NX, NU, NY, T, NW, NV>& plant,
    std::type_identity_t<T>                  Ts,
    const MPCWeights<NU, NY, T>&             weights,
    const MPCConstraints<NU, NY, T>&         constraints,
    const MeasuredDisturbance<NX, NM, T>&    md
) {
    return detail::mpc_design<false, NP, NC, NT_max>(plant, Ts, weights, constraints, md);
}


/**
 * @struct MPCHorizonSuggestion
 * @brief Advisory NP/NC horizon values from suggest_mpc_horizon()
 *
 * NP and NC are template parameters of mpc, so they cannot be
 * chosen at runtime — read the suggestion once at design time and plug the
 * numbers back in as template arguments. When the plant model is constexpr the
 * whole round trip is compile-time:
 *
 * @code
 * constexpr auto h = design::suggest_mpc_horizon(sys);
 * static_assert(h.success);
 * constexpr auto art = design::mpc<h.prediction_horizon, h.control_horizon>(sys, {}, {});
 * @endcode
 */
template<typename T = double>
struct MPCHorizonSuggestion {
    size_t prediction_horizon{0}; ///< Suggested NP
    size_t control_horizon{0};    ///< Suggested NC
    T      time_constant{};       ///< Dominant open-loop time constant τ [s] (pole-based overload)
    T      settling_time{};       ///< Settling-time target the horizon covers [s]
    bool   success{false};        ///< false: no usable dominant pole — use the settling-time overload
};

/**
 * @brief Suggest MPC horizons from an explicit settling-time target
 *
 * The prediction horizon covers the target settling time
 * (@f$ NP = \lceil t_s / T_s \rceil @f$, at least 2) and the control horizon is
 * the standard 20% of it (at least 2) — moves beyond that have little effect
 * on the first move but cost QP size. A very large NP suggests the controller
 * rate is faster than the plant needs; consider a slower MPC tick.
 *
 * @note Compare with the MPC Toolbox™ guidelines for choosing p and m.
 *
 * @param settling_time  Closed-loop settling-time target [s]
 * @param Ts             Controller sample time [s]
 * @return Advisory horizon values (success = false if inputs are non-positive)
 */
template<typename T>
[[nodiscard]] constexpr MPCHorizonSuggestion<T> suggest_mpc_horizon(T settling_time, T Ts) {
    MPCHorizonSuggestion<T> s{};
    if (!(settling_time > T{0}) || !(Ts > T{0})) {
        return s;
    }
    const T steps = settling_time / Ts;
    // NP is a template argument: anything near size_t range is not a usable
    // suggestion, and the cast below would be UB past it.
    if (!(steps < static_cast<T>(std::numeric_limits<size_t>::max() / 2))) {
        return s;
    }
    auto np = static_cast<size_t>(steps);
    if (static_cast<T>(np) < steps) {
        ++np;
    }
    s.prediction_horizon = (np < 2) ? 2 : np;
    const size_t nc = (s.prediction_horizon + 4) / 5; // ceil(NP/5) = 20%
    s.control_horizon = (nc < 2) ? 2 : nc;
    s.settling_time = settling_time;
    s.success = true;
    return s;
}

/**
 * @brief Suggest MPC horizons from the plant's dominant time constant
 *
 * Finds the slowest strictly stable discrete pole of A, converts it to a time
 * constant @f$ \tau = -T_s / \ln|\lambda| @f$, and sizes the horizon to cover
 * @f$ k \cdot \tau @f$ (default k = 4, the ~2% settling criterion) plus the dead time: the
 * steps before every input's move reaches an output, an input delay's registers included, which
 * design::mpc needs the horizon to cover.
 *
 * Fails (success = false) when the plant has no strictly stable pole to read a
 * timescale from — pure integrator chains (e.g. the double integrator) and
 * unstable plants have no intrinsic settling time. Use the explicit
 * settling-time overload there with your closed-loop target.
 *
 * @param sys                    Discrete-time plant (Ts > 0)
 * @param settle_time_constants  Settling target in time constants (default 4 ≈ 2%)
 * @return Advisory horizon values
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr MPCHorizonSuggestion<T> suggest_mpc_horizon(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        settle_time_constants = T{4}
) {
    MPCHorizonSuggestion<T> s{};
    if (!(sys.Ts > T{0})) {
        return s;
    }

    const auto eig = mat::compute_eigenvalues(sys.A);
    //! An m-fold defective unit eigenvalue (an integrator chain) is computed only to
    //! O((eps·‖A‖)^(1/m)): a triple integrator's poles land 1e-6 off 1 in double, and reading
    //! them as a slow stable pole suggests an absurd NP. Their mean is accurate, so the m poles
    //! nearest 1 form an integrator cluster when they lie within 4·NX·(eps·‖A‖)^(1/m) of 1 and
    //! their mean within 4·NX·√(eps·‖A‖); the largest such m is skipped.
    const T                 rel = eps<T>() * sys.A.norm();
    damp::array<size_t, NX> order{};
    damp::array<T, NX>      dist{};
    for (size_t i = 0; i < NX; ++i) {
        order[i] = i;
        dist[i] = damp::abs(eig.values(i) - T{1});
    }
    for (size_t i = 1; i < NX; ++i) {
        for (size_t j = i; j > 0 && dist[order[j]] < dist[order[j - 1]]; --j) {
            damp::swap(order[j], order[j - 1]);
        }
    }
    size_t cluster = 0;
    if (rel > T{0}) {
        const T band_scale = T{4} * static_cast<T>(NX);
        for (size_t m = NX; m > 0 && cluster == 0; --m) {
            damp::complex<T> mean{};
            for (size_t k = 0; k < m; ++k) {
                mean += eig.values(order[k]);
            }
            mean = mean / static_cast<T>(m);
            const T spread = band_scale * damp::exp(damp::log(rel) / static_cast<T>(m));
            if (dist[order[m - 1]] <= spread && damp::abs(mean - T{1}) <= band_scale * damp::sqrt(rel)) {
                cluster = m;
            }
        }
    }
    damp::array<bool, NX> integrator{};
    for (size_t k = 0; k < cluster; ++k) {
        integrator[order[k]] = true;
    }

    T tau_max = T{0};
    for (size_t i = 0; i < NX; ++i) {
        const T mag = damp::abs(eig.values(i));
        if (integrator[i] || !(mag > T{0}) || !(mag < T{1})) {
            continue; // integrators, unstable and deadbeat poles contribute no timescale
        }
        const T tau = -sys.Ts / damp::log(mag);
        tau_max = (tau > tau_max) ? tau : tau_max;
    }
    if (!(tau_max > T{0})) {
        return s; // no strictly stable pole — supply a settling-time target instead
    }

    size_t dead = 0; // steps before every input's move reaches an output; one that never does adds none
    for (const size_t k : detail::first_reaching_steps(sys.A, sys.B, sys.C, NX)) {
        dead = (k < NX && k > dead) ? k : dead;
    }
    s = suggest_mpc_horizon((settle_time_constants * tau_max) + (static_cast<T>(dead) * sys.Ts), sys.Ts);
    s.time_constant = tau_max;
    return s;
}

} // namespace design

namespace detail {

/** @brief Whether @p Design is a filter design of these sizes and scalar, with any noise inputs. */
template<typename Design, size_t NXK, size_t NU, size_t NY, typename T, size_t ND>
inline constexpr bool is_kalman_design = false;

/** @brief A design::kalman result on NXK states, NU inputs, NY readings and ND disturbance states. */
template<size_t NXK, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t ND>
inline constexpr bool is_kalman_design<design::KalmanResult<NXK, NU, NY, T, NW, NV, ND>, NXK, NU, NY, T, ND> = true;

/** @brief The estimator a runtime MPC holds: none (the measurement is the state). */
template<typename Estimator>
struct mpc_estimator_traits {
    static constexpr bool   estimated = false; ///< Whether there is a filter
    static constexpr size_t states = 0;        ///< NXK: the filter's states
    static constexpr size_t outputs = 0;       ///< NY: the readings it takes
    static constexpr size_t disturbances = 0;  ///< ND: its disturbance states

    /** @brief The design results it is built from: none. */
    template<typename Design>
    static constexpr bool built_from = std::is_same_v<Design, NoEstimator>;
};

/** @brief A Luenberger observer: its sizes, and the design::kalman results it runs (any noise inputs). */
template<size_t NXK, size_t NU, size_t NY, typename T, size_t ND, typename Store>
struct mpc_estimator_traits<Luenberger<NXK, NU, NY, T, ND, Store>> {
    static constexpr bool   estimated = true;  ///< Whether there is a filter
    static constexpr size_t states = NXK;      ///< NXK: the filter's states
    static constexpr size_t outputs = NY;      ///< NY: the readings it takes
    static constexpr size_t disturbances = ND; ///< ND: its disturbance states

    /** @brief The design results it is built from: design::kalman's, of its sizes and scalar. */
    template<typename Design>
    static constexpr bool built_from = is_kalman_design<Design, NXK, NU, NY, T, ND>;
};

/** @brief No measured state: the estimator holds the state. */
struct NoMeasuredState {};

/**
 * @brief Runtime constrained MPC controller (fixed per-tick iteration budget)
 *
 * Each tick solves the condensed QP from design::mpc() at the state the command will meet and
 * applies the first move: u = u_prev + Δu₀. The applied command is always clamped to the u box,
 * and to the Δu box whenever the held input lies inside the u box, so even a budget-exhausted or
 * failed solve never emits an out-of-range actuator command (the QP's intermediate iterates are
 * dual, not primal, feasible). A MaxIterations or Degenerate iterate is applied when the
 * outputs are soft (the default). If any output ECR is 0 or a terminal set is present, the
 * iterate is applied only when it is primal feasible; otherwise the tick holds the last
 * applied input, as Infeasible does. Holding the input does not make the predicted outputs or
 * a terminal set feasible. The soft-constraint
 * slack actually used is reported by last_slack(). @p r is the reference for the NY controlled
 * outputs (held, or a Preview), not a state x_ref: the cost is ‖y − r‖²_Qy plus move /
 * input-target terms, with no ‖x − x_ref‖ term (see the file comment).
 *
 * One tick is @ref step(): @ref predict holds the input the plant received over [k−1, k)
 * (the base of the next Δu) and runs the estimator's time update, @ref update takes the
 * measurement m_k, and @ref control solves for u_k. The three are public, for the low-latency
 * order update → control → write the output → back_calculate → predict (the trailing predict is
 * the next sample's time update). When a downstream stage can
 * clamp or override the command, @ref back_calculate records the input that went out, so the
 * next move and the estimator start from it; without it predict holds the last command.
 *
 * - No estimator (Estimator = NoEstimator): the measurement is the plant state x; the reference
 *   is on the controlled outputs C·x.
 * - Luenberger (the fixed gain of design::kalman): the measurement is the plant's outputs; the
 *   observer reconstructs the state.
 *   With a disturbance model the estimator holds deviations from the steady-state target,
 *   re-anchored each tick on the reference now (SteadyStateTarget), and the law is offset-free.
 *   Without one the QP is posed about the nominal target of the reference.
 *
 * @tparam NX        Plant states
 * @tparam NY        Controlled outputs: the reference, the weights and the y limits (the rows of C on the measured state)
 * @tparam ND        Measured-disturbance inputs (0: none)
 * @tparam Solver    QP solver policy (default: the warm-started Goldfarb–Idnani active set)
 * @tparam NT_max    Terminal-set row capacity (default 0)
 * @tparam Estimator Luenberger, or NoEstimator for the measured state; it holds the
 *                   disturbance when there is one
 *
 * @note Compare with MATLAB®'s mpcmove(mpcobj, xc, ym, r, v).
 * @see design::mpc
 */
template<size_t NX, size_t NU, size_t NY, size_t NP, size_t NC, size_t ND, typename T, typename Solver, size_t NT_max, typename Estimator>
    requires std::is_floating_point_v<T>
class MPCCore {
    using Traits = detail::mpc_estimator_traits<Estimator>;

public:
    static constexpr bool   estimated = Traits::estimated;                  ///< Whether a Kalman filter reconstructs the state
    static constexpr size_t NDE = Traits::disturbances;                     ///< Estimated disturbance states (0: not offset-free)
    static constexpr size_t NXM = NX;                                       ///< States of the prediction model
    static constexpr size_t NYM = estimated ? Traits::outputs : NX;         ///< Measured outputs (the state without an estimator)
    static constexpr size_t NXE = estimated ? (Traits::states - NDE) : NXM; ///< States of the estimator's model before the disturbance
    static_assert(!estimated || ND == 0, "A measured disturbance with an estimator is not supported.");
    static_assert(std::is_same_v<Estimator, NoEstimator> || estimated, "The estimator is a Luenberger observer, or NoEstimator for the measured state.");

    using Artifacts = design::MPCArtifacts<NXM, NU, NY, NP, NC, ND + NDE, T, NT_max>;
    using TargetResult = design::SteadyStateTargetResult<NXE, NU, ND + NDE, NY, T>;
    /** @brief A design it is built from: the QP and the target, and @p Kalman, the estimator's design (NoEstimator without one). */
    template<typename Kalman>
    using ResultWith = design::MPCResult<Artifacts, Kalman, TargetResult>;
    static constexpr size_t NXA = Artifacts::NXA;
    static constexpr size_t NZ = Artifacts::NZ;
    static constexpr size_t NZ1 = Artifacts::NZ1;
    static constexpr size_t NI1 = Artifacts::NI1;
    static constexpr size_t NYS = Artifacts::NYS;
    /** @brief The output reference over the prediction horizon, one value per step (MATLAB's `ref` with preview). */
    using Preview = damp::array<ColVec<NY, T>, NP>;
    /** @brief The input target over the control horizon, one value per move (MATLAB's `mv.target` with preview). */
    using InputPreview = damp::array<ColVec<NU, T>, NC>;

    constexpr MPCCore() = default;

    /** @brief From the design. A failed design is a compile error when constant-evaluated (`static constinit MPC law{art};`). */
    template<typename Kalman>
        requires(Traits::template built_from<Kalman>)
    constexpr explicit MPCCore(const ResultWith<Kalman>& result)
        : art(result.mpc), estimator_(make_estimator(result.kalman)), target_(target_of(result.target)), u_target_(held_target(result.mpc.u_target)) {
        if (!result.success && std::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(result.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    /** @brief Construct with a specific solver instance (for stateful solver policies) */
    template<typename Kalman>
        requires(Traits::template built_from<Kalman>)
    constexpr MPCCore(const ResultWith<Kalman>& result, const Solver& solver)
        : art(result.mpc), qp_solver_(solver), estimator_(make_estimator(result.kalman)), target_(target_of(result.target)), u_target_(held_target(result.mpc.u_target)) {
        if (!result.success && std::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(result.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    /**
     * @brief The output reference over the prediction horizon: one value held, or a preview
     *
     * Built implicitly from a ColVec<NY> (held over the NP steps) or from a Preview (the
     * reference at each of steps 1 … NP), so every tick takes either. From a continuous plant
     * r[k] is the reference over period k, [k·Ts, (k+1)·Ts) (the cost integrates each period),
     * so r[0] covers the period starting now. With a preview the plan
     * sees a setpoint change before it happens and starts moving ahead of it; a trajectory
     * generator hands in its next NP points each tick (repeat the last known value to fill).
     * Output references, not state references.
     */
    struct Reference {
        Preview r{}; ///< The output reference at steps 1 … NP (periods 0 … NP − 1 from a continuous plant)

        constexpr Reference(const ColVec<NY, T>& held_value) : r(held(held_value)) {} // NOLINT(google-explicit-constructor)
        constexpr Reference(const Preview& preview) : r(preview) {}                   // NOLINT(google-explicit-constructor)

        /** @brief A SISO reference from its scalar, held over the NP steps (as ColVec<1>'s). */
        constexpr Reference(T held_value) // NOLINT(google-explicit-constructor)
            requires(NY == 1)
            : r(held(ColVec<NY, T>{held_value})) {}
    };

    /**
     * @brief Time update with the input the plant received over the last period (back_calculate's,
     *        else the last command): the base of the next Δu and the estimator's predict.
     */
    constexpr void predict() {
        applied_.take();
        if constexpr (estimated) {
            estimator_.predict(applied_.applied);
        }
    }

    /**
     * @brief Measurement update with @p y: the plant's outputs through the filter, or the plant
     *        state itself without an estimator.
     */
    constexpr void update(const ColVec<NYM, T>& y) {
        if constexpr (estimated) {
            estimator_.update(y);
        } else {
            x_ = y;
        }
    }

    /**
     * @brief The control move at the current state: u = u_prev + Δu₀
     *
     * The law alone: it solves from the held input and records the command, but does not hold
     * it; @ref predict does, with the input the plant actually received. The QP is posed about
     * the steady-state target of the reference now, r[0]: under the estimated disturbance when
     * the design is offset-free (re-anchoring the estimator on it), the nominal target
     * otherwise. With d periods of delay the preview the plan reads is r shifted by d.
     *
     * @param r  Controlled-output reference (NY held, or a Preview); not a state x_ref
     * @param d  Measured disturbance, held over the horizon (ND > 0; nothing otherwise)
     * @return Input command u (clamped to the configured u and Δu boxes)
     */
    [[nodiscard]] constexpr ColVec<NU, T> control(const Reference& r, const ColVec<ND, T>& d) {
        const Preview& ahead = r.r;
        if constexpr (NDE > 0) {
            static_cast<void>(d); // no measured disturbance alongside an estimated one
            target_.anchor(estimator_, r.r[0]);
            const auto x_dev = ::damp::detail::leading<NXM>(estimator_.at_command());
            const auto x_a = ::damp::detail::leading<NXM>(estimator_.anchor_state());
            const auto d_a = estimator_.anchor_disturbance();
            return record(control_impl(ahead, ColVec<NXM, T>(x_dev + x_a), d_a, ColVec<NDE, T>{}, x_a, estimator_.anchor_input(), x_dev));
        } else {
            const ColVec<NXM, T> x = law_state();
            const auto           t = target_.equilibrium(r.r[0], d); // the steady state under the measured d
            const auto           x_t = ::damp::detail::leading<NXM>(t.x);
            return record(control_impl(ahead, x, d, ColVec<ND, T>{}, x_t, t.u, ColVec<NXM, T>(x - x_t)));
        }
    }

    /**
     * @brief Record the input the plant received, @p u_sat (after any downstream clamp or
     *        override, with any feedforward @p u_ff added outside the law): the next @ref predict
     *        holds it as the base of the next Δu and gives it to the estimator.
     *
     * The next move plans the total input from u_sat and returns the law's share, the planned
     * total less u_ff: adding u_ff back outside the law applies the plan, and the total lies in
     * the u box. Also seeds the held input for a bumpless takeover from another controller:
     * record the input that controller applied, and the next step moves from it.
     */
    constexpr void back_calculate(const ColVec<NU, T>& u_sat, const ColVec<NU, T>& u_ff = ColVec<NU, T>{}) { applied_.record(u_sat, u_ff); }

    /**
     * @brief One tick: predict(), update(y), then the move.
     * @param r Controlled-output reference (NY held, or a Preview); not a state x_ref
     * @param y The measurement: the plant's outputs with an estimator, the state without
     * @return Input command u (clamped to the configured u and Δu boxes)
     */
    [[nodiscard]] constexpr ColVec<NU, T> step(const Reference& r, const ColVec<ND, T>& d, const ColVec<NYM, T>& y) {
        predict();
        update(y);
        return control(r, d);
    }

    /** @brief Reset the controller memory (the held input returns to zero); a state estimate is left in place. */
    constexpr void reset() {
        applied_ = ::damp::detail::AppliedRecord<NU, T>{};
        if constexpr (estimated) {
            estimator_.reset();
        }
        last_status_ = design::QPStatus::Success;
        last_iterations_ = 0;
        last_slack_ = T{0};
        last_cost_ = T{0};
    }

    /** @brief Override the preferred input operating point (physical units) */
    constexpr void set_input_target(const ColVec<NU, T>& u_target) { u_target_ = held_target(u_target); }

    /** @brief Input target per move over the control horizon (MATLAB's mv.target with preview); weighted by Ru. */
    constexpr void set_input_target(const InputPreview& u_target) { u_target_ = u_target; }

    /**
     * @brief Change the limits on u, Δu and the controlled outputs (MATLAB®'s online constraints)
     *
     * The QP keeps its Hessian, constraint matrix, scaling and softness; the bound vector and the
     * fallback clamp take the new limits from the next control. Refused (false, limits unchanged)
     * when a limit is NaN, a lower limit is above its upper one, or the design has a terminal
     * invariant set, which was built for the design's limits.
     */
    constexpr bool set_limits(const design::MPCLimits<NU, NY, T>& lim) {
        if (art.terminal == design::MPCTerminal::invariant_set) {
            return false;
        }
        if (!design::detail::mpc_limits_valid(lim.u_lim) || !design::detail::mpc_limits_valid(lim.du_lim) || !design::detail::mpc_limits_valid(lim.y_lim)) {
            return false;
        }
        design::detail::mpc_box_bounds(art, lim);
        art.constraints.u_lim = lim.u_lim;
        art.constraints.du_lim = lim.du_lim;
        art.constraints.y_lim = lim.y_lim;
        return true;
    }

    /** @brief The limits in force (the design's, or the last set_limits). */
    [[nodiscard]] constexpr design::MPCLimits<NU, NY, T> limits() const {
        const auto& c = art.constraints;
        return {c.u_lim, c.du_lim, c.y_lim};
    }

    /** @brief Access the solver policy (stateful solvers expose warm-start caches here) */
    [[nodiscard]] constexpr Solver&       qp_solver() { return qp_solver_; }
    [[nodiscard]] constexpr const Solver& qp_solver() const { return qp_solver_; }

    /** @brief The estimator (the observer and its anchor). */
    [[nodiscard]] constexpr const Estimator& estimator() const
        requires(estimated)
    {
        return estimator_;
    }

    /** @brief The input the last predict held: the base of the next Δu (back_calculate's, else the command before it). */
    [[nodiscard]] constexpr const ColVec<NU, T>& applied() const { return applied_.applied; }
    /** @brief The command the last control returned: the input a predict without back_calculate holds. */
    [[nodiscard]] constexpr const ColVec<NU, T>& last_command() const { return applied_.command; }
    /** @brief The input target of the first move (physical units). */
    [[nodiscard]] constexpr const ColVec<NU, T>& input_target() const { return u_target_[0]; }
    /** @brief The input target of every move over the control horizon. */
    [[nodiscard]] constexpr const InputPreview& input_target_preview() const { return u_target_; }
    /** @brief The status of the last QP solve. */
    [[nodiscard]] constexpr design::QPStatus last_status() const { return last_status_; }
    /** @brief The iterations of the last QP solve. */
    [[nodiscard]] constexpr size_t last_iterations() const { return last_iterations_; }
    /** @brief Soft-constraint slack ε from the last solve (0 = no output bound was relaxed) */
    [[nodiscard]] constexpr T last_slack() const { return last_slack_; }
    /** @brief QP objective ½ζᵀHζ + fᵀζ from the last solve (for a Lyapunov / monotonicity check) */
    [[nodiscard]] constexpr T last_cost() const { return last_cost_; }

    /** @brief The law's model state: the estimate (absolute; with the register of a delay), or the last measurement. */
    [[nodiscard]] constexpr ColVec<NX, T> state_estimate() const {
        if constexpr (estimated) {
            return ::damp::detail::leading<NX>(estimator_.state_estimate());
        } else {
            return x_;
        }
    }

    /** @brief The estimated disturbance d̂ (what the integrating model has absorbed). */
    [[nodiscard]] constexpr ColVec<NDE, T> disturbance_estimate() const
        requires(NDE > 0)
    {
        return estimator_.disturbance_estimate();
    }

    /**
     * @brief NP copies of one reference: the held-setpoint case of a preview. Built by copy-construction
     * (no element assignment), so it is a constant expression for a constinit controller.
     */
    [[nodiscard]] static constexpr Preview held(const ColVec<NY, T>& r) {
        return [&]<size_t... I>(damp::index_sequence<I...>) { return Preview{{(static_cast<void>(I), r)...}}; }(damp::make_index_sequence<NP>{});
    }

    /** @brief NC copies of one input target. */
    [[nodiscard]] static constexpr InputPreview held_target(const ColVec<NU, T>& u) {
        return [&]<size_t... I>(damp::index_sequence<I...>) { return InputPreview{{(static_cast<void>(I), u)...}}; }(damp::make_index_sequence<NC>{});
    }

private:
    using Target = SteadyStateTarget<NXE, NU, ND + NDE, NY, T>;
    using MeasuredState = std::conditional_t<estimated, detail::NoMeasuredState, ColVec<NX, T>>;

    /** @brief The estimator from its design; nothing on the measured state. */
    template<typename Kalman>
    [[nodiscard]] static constexpr Estimator make_estimator(const Kalman& design) {
        if constexpr (requires { Estimator(design); }) {
            return Estimator(design);
        } else {
            static_cast<void>(design); // nothing to build: no filter, or its matrices are constants of its type
            return Estimator{};
        }
    }

    /** @brief The runtime target from its design, refused only by this constructor. */
    [[nodiscard]] static constexpr Target target_of(const TargetResult& design) {
        Target t{};
        t.target_d = design.target_d;
        t.target_r = design.target_r;
        return t;
    }

    /** @brief The law's state the command decided now acts on (absolute): the estimate, or the measurement. */
    [[nodiscard]] constexpr ColVec<NXM, T> law_state() const {
        if constexpr (estimated) {
            return ::damp::detail::leading<NXM>(estimator_.at_command());
        } else {
            return x_;
        }
    }

    [[nodiscard]] constexpr ColVec<NU, T> record(const ColVec<NU, T>& u) {
        applied_.command = u;
        return u;
    }

    /**
     * @brief The QP about the target (x_t, u_t): @p x the absolute state for the bounds,
     *        @p x_dev = x − x_t for the gradient, @p dq the disturbance channels (absolute) and
     *        @p dq_dev their deviation from the target's (zero: the target holds under the
     *        estimated or measured disturbance).
     */
    constexpr ColVec<NU, T> control_impl(
        const Preview&             r,
        const ColVec<NXM, T>&      x,
        const ColVec<ND + NDE, T>& dq,
        const ColVec<ND + NDE, T>& dq_dev,
        const ColVec<NXM, T>&      x_t,
        const ColVec<NU, T>&       u_t,
        const ColVec<NXM, T>&      x_dev
    ) {
        constexpr size_t NDQ = ND + NDE;
        // A failed (or default) design has zero scale factors: hold the last input.
        if (!art.success) {
            return ColVec<NU, T>(applied_.applied - applied_.feedforward);
        }
        constexpr T sentinel = design::unbounded_bound<T>();
        constexpr T skip_threshold = sentinel / T{4};

        // x̃0 = [x; u'_prev; d] for the bounds, and its deviation from the target for the gradient.
        ColVec<NXA, T> xa{};
        ColVec<NXA, T> xa_dev{};
        for (size_t i = 0; i < NXM; ++i) {
            xa(i) = x(i);
            xa_dev(i) = x_dev(i);
        }
        for (size_t j = 0; j < NU; ++j) {
            xa(NXM + j) = applied_.applied(j) / art.scale_u(j);
            xa_dev(NXM + j) = (applied_.applied(j) - u_t(j)) / art.scale_u(j);
        }
        for (size_t j = 0; j < NDQ; ++j) {
            xa(NXM + NU + j) = dq(j);
            xa_dev(NXM + NU + j) = dq_dev(j);
        }

        // Scaled deviations: the reference at every prediction step from C x_t, the input target
        // at every move from u_t. The target is an equilibrium of the prediction model, so its own
        // free response contributes nothing and this f is the absolute form's in exact
        // arithmetic, without its difference of two large terms in float. The bounds keep the
        // absolute x̃0.
        ColVec<NY, T> y_t = art.C * x_t; // the target's output, with any estimated output bias
        if constexpr (NDQ > 0) {
            y_t = y_t + (art.Cd * dq);
        }
        ColVec<NYS, T> r_dev{};
        for (size_t kb = 0; kb < NP; ++kb) {
            for (size_t i = 0; i < NY; ++i) {
                r_dev((kb * NY) + i) = (r[kb](i) - y_t(i)) / art.scale_y(i);
            }
        }
        ColVec<NZ, T> ut_dev{};
        for (size_t jm = 0; jm < NC; ++jm) {
            for (size_t j = 0; j < NU; ++j) {
                ut_dev((jm * NU) + j) = (u_target_[jm](j) - u_t(j)) / art.scale_u(j);
            }
        }
        const ColVec<NZ, T> f_z = art.F_x * xa_dev - art.F_r * r_dev - art.F_ut * ut_dev;
        ColVec<NZ1, T>      f{};
        for (size_t i = 0; i < NZ; ++i) {
            f(i) = f_z(i);
        }

        const auto       shift = art.B_map * xa;
        const auto       shift_dev = art.B_map * xa_dev;
        ColVec<NI1, T>   b{};
        constexpr size_t r_term = (4 * NZ) + (2 * NP * NY);
        for (size_t i = 0; i < NI1; ++i) {
            const T sh = (i >= r_term && i < r_term + NT_max) ? shift_dev(i) : shift(i);
            b(i) = (art.b_bound(i) >= skip_threshold) ? sentinel : art.b_bound(i) + sh;
        }
        if constexpr (NT_max > 0) {
            //! Terminal rows are G·z_N + G_target·q ≤ h: the target q = [y_t; u_t] moves each offset.
            ColVec<NY + NU, T> q{};
            for (size_t i = 0; i < NY; ++i) {
                q(i) = y_t(i);
            }
            for (size_t j = 0; j < NU; ++j) {
                q(NY + j) = u_t(j);
            }
            const ColVec<NT_max, T> tq = art.G_terminal_target * q;
            for (size_t row = 0; row < art.n_terminal; ++row) {
                b(r_term + row) -= tq(row);
            }
        }

        // The default WarmStartActiveSetSolver caches chol(H) across ticks and
        // seeds each solve with the previous active set; the artifacts' QP
        // matrices are constant by construction, so its cache contract holds.
        const auto qp = qp_solver_(art.H, f, art.A_con, b, art.max_qp_iterations);
        last_status_ = qp.status;
        last_iterations_ = qp.iterations;
        last_slack_ = qp.x(NZ);
        last_cost_ = qp.objective;

        bool take = qp.status == design::QPStatus::Success;
        if (!take && (qp.status == design::QPStatus::MaxIterations || qp.status == design::QPStatus::Degenerate)) {
            take = !hard_output_or_terminal() || qp_primal_feasible(qp.x, b);
        }
        ColVec<NU, T> du{};
        if (take) {
            for (size_t j = 0; j < NU; ++j) {
                du(j) = qp.x(j) * art.scale_u(j); // back to physical units
            }
        }

        // Safety clamp: guarantees box feasibility of the applied command. The plan is the total
        // input; a feedforward recorded with it is added outside the law, so the command is the
        // law's share.
        du = art.constraints.du_lim.saturate(du);
        return ColVec<NU, T>(art.constraints.u_lim.saturate(applied_.applied + du) - applied_.feedforward);
    }

    [[nodiscard]] constexpr bool hard_output_or_terminal() const {
        if (art.n_terminal > 0) {
            return true;
        }
        for (size_t i = 0; i < NY; ++i) {
            if (art.constraints.y_lower_ecr(i) == T{0} || art.constraints.y_upper_ecr(i) == T{0}) {
                return true;
            }
        }
        return false;
    }

    /** @brief True when A_con z ≤ b on every enabled (non-sentinel) row. */
    [[nodiscard]] constexpr bool qp_primal_feasible(const ColVec<Artifacts::NZ1, T>& z, const ColVec<Artifacts::NI1, T>& b) const {
        constexpr T skip = design::unbounded_bound<T>() / T{4};
        constexpr T c = T{32};
        const T     e = eps<T>();
        for (size_t i = 0; i < Artifacts::NI1; ++i) {
            if (!(b(i) < skip)) {
                continue;
            }
            T az = T{0};
            T s = damp::abs(b(i));
            for (size_t j = 0; j < Artifacts::NZ1; ++j) {
                const T term = art.A_con(i, j) * z(j);
                az += term;
                s += damp::abs(term);
            }
            if (az > b(i) + (c * e * s)) {
                return false;
            }
        }
        return true;
    }

    Artifacts                            art{};
    Solver                               qp_solver_{};
    [[no_unique_address]] Estimator      estimator_{}; ///< The filter and its anchor (no storage on the measured state)
    [[no_unique_address]] MeasuredState  x_{};         ///< The last measured state (no storage with an estimator)
    Target                               target_{};    ///< The maps to the steady-state target
    ::damp::detail::AppliedRecord<NU, T> applied_{};   ///< The last command, and the held input (the base of Δu)
    InputPreview                         u_target_{};
    design::QPStatus                     last_status_{design::QPStatus::Success};
    size_t                               last_iterations_{0};
    T                                    last_slack_{0};
    T                                    last_cost_{0};
};

} // namespace detail

/**
 * @ingroup discrete_controllers
 * @brief Constrained linear MPC on the measured state (design::mpc)
 *
 * The tick is update(x) then control(r), or step(r, x) for both; @p r is the controlled outputs'
 * reference (held, or a Preview over the horizon). A design with a measured disturbance takes it
 * as well, control(r, d) and step(r, d, x); one without has no d argument. back_calculate(u_sat)
 * records what went out, predict() holds it as the base of the next move. A failed design does
 * not compile when constant-evaluated.
 *
 * @code
 * static constinit MPC law{kRes.as<float>()};
 * law.update(x);
 * u = law.control(r);
 * @endcode
 *
 * @note Compare with MATLAB®'s mpcmove(mpcobj, xc, ym, r, v) on a full-state model.
 */
template<
    size_t NX,
    size_t NU,
    size_t NY,
    size_t NP,
    size_t NC,
    size_t ND = 0,
    typename T = float,
    typename Solver = design::WarmStartActiveSetSolver<(NC * NU) + 1, (4 * NC * NU) + (2 * NP * NY) + 1, T>,
    size_t NT_max = 0>
class MPC : private detail::MPCCore<NX, NU, NY, NP, NC, ND, T, Solver, NT_max, NoEstimator> {
    using Core = detail::MPCCore<NX, NU, NY, NP, NC, ND, T, Solver, NT_max, NoEstimator>;

public:
    using Reference = Core::Reference;                     ///< r: held, or a Preview
    using Preview = Core::Preview;                         ///< The reference at steps 1 … NP
    using InputPreview = Core::InputPreview;               ///< The input target per move
    using Result = Core::template ResultWith<NoEstimator>; ///< The design it is built from
    static constexpr size_t NZ1 = Core::NZ1;               ///< QP decision size (moves + slack)
    static constexpr size_t NI1 = Core::NI1;               ///< QP constraint rows

    constexpr MPC() = default;

    /** @brief From design::mpc. A failed design is a compile error when constant-evaluated. */
    constexpr explicit MPC(const Result& result) : Core(result) {}

    /** @brief From design::mpc with a specific solver instance (for stateful solver policies). */
    constexpr MPC(const Result& result, const Solver& solver) : Core(result, solver) {}

    /** @brief Hold the measured state @p x for the next control. */
    constexpr void update(const ColVec<NX, T>& x) { Core::update(x); }

    /** @brief The move at the state of the last update: u = u_prev + Δu₀, clamped to the u and Δu boxes. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const Reference& r)
        requires(ND == 0)
    {
        return Core::control(r, ColVec<0, T>{});
    }

    /** @brief The move at the state of the last update under the measured disturbance @p d, held over the horizon. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const Reference& r, const ColVec<ND, T>& d)
        requires(ND > 0)
    {
        return Core::control(r, d);
    }

    /** @brief One tick: predict(), update(@p x), then control(@p r). */
    [[nodiscard]] constexpr ColVec<NU, T> step(const Reference& r, const ColVec<NX, T>& x)
        requires(ND == 0)
    {
        return Core::step(r, ColVec<0, T>{}, x);
    }

    /** @brief One tick: predict(), update(@p x), then control(@p r, @p d). */
    [[nodiscard]] constexpr ColVec<NU, T> step(const Reference& r, const ColVec<ND, T>& d, const ColVec<NX, T>& x)
        requires(ND > 0)
    {
        return Core::step(r, d, x);
    }

    using Core::applied;
    using Core::back_calculate;
    using Core::input_target;
    using Core::input_target_preview;
    using Core::last_command;
    using Core::last_cost;
    using Core::last_iterations;
    using Core::last_slack;
    using Core::last_status;
    using Core::limits;
    using Core::predict;
    using Core::qp_solver;
    using Core::reset;
    using Core::set_input_target;
    using Core::set_limits;
    using Core::state_estimate;
};

/** @brief MPC law{res}: the sizes, scalar and terminal capacity of design::mpc's result. */
template<size_t NXM, size_t NU, size_t NY, size_t NP, size_t NC, size_t ND, typename T, size_t NT>
MPC(const design::MPCResult<design::MPCArtifacts<NXM, NU, NY, NP, NC, ND, T, NT>, NoEstimator, design::SteadyStateTargetResult<NXM, NU, ND, NY, T>>&)
    -> MPC<NXM, NU, NY, NP, NC, ND, T, design::WarmStartActiveSetSolver<(NC * NU) + 1, (4 * NC * NU) + (2 * NP * NY) + NT + 1, T>, NT>;

/** @brief MPC law{res, solver}: as MPC law{res}, with that solver. */
template<size_t NXM, size_t NU, size_t NY, size_t NP, size_t NC, size_t ND, typename T, size_t NT, typename Solver>
MPC(const design::MPCResult<design::MPCArtifacts<NXM, NU, NY, NP, NC, ND, T, NT>, NoEstimator, design::SteadyStateTargetResult<NXM, NU, ND, NY, T>>&, const Solver&)
    -> MPC<NXM, NU, NY, NP, NC, ND, T, Solver, NT>;

namespace detail {

/** @brief The runtime an MPC design runs as, in the scalar @p T: none for other results (ofmpc.hpp adds OFMPC). */
template<typename R, typename T>
struct mpc_runtime {};

/** @brief design::mpc runs as MPC. */
template<size_t NXM, size_t NU, size_t NY, size_t NP, size_t NC, size_t ND, typename U, size_t NT, typename T>
struct mpc_runtime<design::MPCResult<design::MPCArtifacts<NXM, NU, NY, NP, NC, ND, U, NT>, NoEstimator, design::SteadyStateTargetResult<NXM, NU, ND, NY, U>>, T> {
    using type = MPC<NXM, NU, NY, NP, NC, ND, T, design::WarmStartActiveSetSolver<(NC * NU) + 1, (4 * NC * NU) + (2 * NP * NY) + NT + 1, T>, NT>; ///< MPC
};

} // namespace detail

namespace detail {

/** @brief An MPC design. */
template<const auto& Res, typename T>
concept mpc_design = requires { typename mpc_runtime<std::remove_cvref_t<decltype(Res)>, T>::type; };

/** @brief An MPC design with a filter: OFMPC, whose filter can run on the design's constants. */
template<const auto& Res, typename T>
concept filtered_mpc_design = mpc_design<Res, T> && mpc_runtime<std::remove_cvref_t<decltype(Res)>, T>::filtered;

} // namespace detail

/** @brief Law<kRes, T> of a `constexpr` OFMPC design: the same law, its filter on the design's matrices as constants of the type. The QP's data stays in members: set_limits rewrites it. */
template<const auto& Res, typename T>
    requires(detail::filtered_mpc_design<Res, T> && detail::constant_design<Res>)
struct Law<Res, T> : detail::mpc_runtime<std::remove_cvref_t<decltype(Res)>, T>::template on_constants<Res> {
    using Base = detail::mpc_runtime<std::remove_cvref_t<decltype(Res)>, T>::template on_constants<Res>; ///< The runtime on constants
    using Runtime = detail::mpc_runtime<std::remove_cvref_t<decltype(Res)>, T>::type;                    ///< The runtime with members: the type of a member a redesign assigns
    static_assert(detail::constant_design_runs(Res.success, Res.status));                                // a failed design does not compile: the error names its status

    /** @brief The design's law. */
    constexpr Law() : Base{Res.template as<T>()} {}
};

/** @brief Law<kRes, T> of an MPC design: MPC or OFMPC, default-constructed from the design. */
template<const auto& Res, typename T>
    requires detail::mpc_design<Res, T>
struct Law<Res, T> : detail::mpc_runtime<std::remove_cvref_t<decltype(Res)>, T>::type {
    using Runtime = detail::mpc_runtime<std::remove_cvref_t<decltype(Res)>, T>::type; ///< The runtime the design runs as

    /** @brief The design's law. */
    constexpr Law() : Runtime{Res.template as<T>()} {}

    using Runtime::operator=;
};

} // namespace damp
