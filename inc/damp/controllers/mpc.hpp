// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file mpc.hpp
 * @brief Constrained linear MPC: condensed-QP synthesis (`design::mpc`) and the fixed-budget
 *        runtimes, `damp::MPC` on the measured state and `damp::OFMPC` through a Kalman
 *        estimator, offset-free with a disturbance model, with an opt-in DARE terminal cost and Gilbert–Tan positively
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
 *   J = \sum_{k=1}^{NP} \lVert y_k - y_{r,k} \rVert^2_{Q_y}
 *     + \sum_{k=0}^{NC-1} \lVert \Delta u_k \rVert^2_{R_{\Delta u}}
 *     + \sum_{k=0}^{NC-1} \lVert u_k - u_{t,k} \rVert^2_{R_u}
 *     + \rho_\varepsilon\, \varepsilon^2,
 * @f]
 * on every output of the plant. The reference r is on the outputs the steady-state target holds
 * (every output, or those a TrackedOutputs or a disturbance spec lists): @f$ y_{r,k} @f$ is r_k
 * itself when the target holds every output, else the target's output for r_k, the listed outputs
 * at r_k and the others where that steady state puts them. @f$ u_{t,k} @f$ is the target's input
 * for r_k, the steady input that holds it, so Ru leaves no offset; an input target set at run time
 * (set_input_target) takes its place. With @f$ Q_y @f$ replaced by @f$ Q_{y,N} @f$ on the final step (a terminal
 * output weight; set it to the DARE solution with @f$ C = I @f$ for the
 * guaranteed-stability short-horizon setup). The option @ref damp::design::DareCost "DareCost"
 * replaces that last output weight with @f$ V_f = z_N^\top P z_N @f$, the cost-to-go of this
 * stage cost itself: a DARE on @f$ z = [x - x_t;\; u'_{prev} - u'_t] @f$ with the move as
 * input, so it is a control Lyapunov function for the Δu form. @ref damp::design::InvariantSet "InvariantSet"
 * adds the maximal positively invariant set of that terminal law, augmented with the constant
 * target (MPC for tracking), as hard inequalities on @f$ z_N @f$: one set valid for every
 * target that keeps @c target_margin of each half-range clear of its bounds, the
 * target moving each row's offset at run time. Without either the cost is the plain
 * Qy / Qy_terminal one. Subject to box constraints on @f$ u @f$ and @f$ \Delta u @f$ (hard, an
 * @ref damp::design::InputLimits "InputLimits"), and on the outputs an
 * @ref damp::design::OutputLimits "OutputLimits" lists (soft by default) over the horizon.
 *
 * That cost, with weights per step, is the discrete-plant form's. The continuous-plant form
 * design::mpc(plant, Ts, …) takes the weights per second and integrates the output error over
 * each period instead of sampling it,
 * @f[
 *   J = \sum_{k=0}^{NP-1} \int_{kT_s}^{(k+1)T_s} \lVert y(t) - y_{r,k} \rVert^2_{Q_y} dt
 *     + \sum_{k=0}^{NC-1} T_s \Bigl( \bigl\lVert \tfrac{\Delta u_k}{T_s} \bigr\rVert^2_{R_{\Delta u}}
 *     + \lVert u_k - u_{t,k} \rVert^2_{R_u} \Bigr) + \rho_\varepsilon\, \varepsilon^2,
 * @f]
 * with @f$ y(t) = C x(t) + D u_k @f$ under the held input, @f$ R_{\Delta u} @f$ on the rate
 * @f$ \Delta u/T_s @f$ (the Riemann sum of @f$ \int \dot u^\top R_{\Delta u} \dot u\, dt @f$)
 * and @f$ Q_{y,N} @f$ on the last period. Each period's integral is a quadratic form in
 * @f$ [x_k;\, u_k;\, d;\, y_{r,k}] @f$, cross terms included, computed exactly by one Van Loan
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
 * softness coefficients V (`OutputLimits::lower_ecr` / `upper_ecr`, default 1)
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
 * plant's C still names the outputs the QP weights: C = [1 0] on a position–velocity plant tracks
 * the position with a scalar r and Qy, and C = I with TrackedOutputs{{0}} weights both states
 * while r holds the position alone. There is deliberately no separate state weight
 * (MATLAB®'s linear MPC has none either).
 *
 * r is the reference of the outputs the target holds. There is no LQ-style state-reference tick
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
 * Feedthrough is not supported: the outputs the QP weighs or bounds must have @f$ D = 0 @f$
 * (synthesis fails otherwise). Offset-free tracking under plant/model mismatch or unmeasured
 * disturbances needs a disturbance model and an estimator (design::mpc with a KalmanNoise and a disturbance spec); the velocity form alone rejects only what its integral action sees
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
 * @see estimation/kalman.hpp for the estimator's design, estimation/fixed_gain_observer.hpp for the observer
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
 *     .C = {{1.0, 0.0}}, // the output the QP weights and r holds: the position
 *     .Ts = 0.1,
 * };
 * constexpr design::InputLimits limits{.u_lim = {.lower = -0.5, .upper = 0.5}, .du_lim = {.lower = -0.2, .upper = 0.2}};
 * constexpr auto art = design::mpc<15, 5>(sys, {}, limits);
 * // The position kept under 1.02, hard, and a DARE terminal cost:
 * constexpr auto capped = design::mpc<15, 15>(sys, {}, limits, design::OutputLimits{.outputs = {0}, .lim = {.upper = 1.02}, .upper_ecr = {0.0}},
 *                                             design::DareCost{});
 * static constinit MPC controller{art.as<float>()}; // a failed design does not compile
 * // each tick: u = controller.step(ColVec<1, float>{r}, x); // the state x measured
 * //            controller.back_calculate(u_applied);         // when a downstream stage clamps it
 *
 * // From the position reading alone, offset-free under an unmeasured load at the input:
 * constexpr StateSpace<2, 1, 1, double, 2, 1> plant{.A = sys.A, .B = sys.B, .C = sys.C, .Ts = 0.1};
 * constexpr auto of = design::mpc<15, 5>(plant, {}, design::InputLimits{.u_lim = {.lower = -0.5, .upper = 0.5}},
 *                                          design::KalmanNoise{Qn, Rn, Qd}, design::InputDisturbance{});
 * static constinit OFMPC offset_free{of.as<float>()};
 * // each tick: u = offset_free.step(ColVec<1, float>{r}, ColVec<1, float>{y});
 * @endcode
 */

#include "damp/backend.hpp"
#include "damp/controllers/law.hpp" // IWYU pragma: export
#include "damp/controllers/state_feedback_types.hpp"
#include "damp/design/disturbance.hpp"  // IWYU pragma: export
#include "damp/design/input_limits.hpp" // IWYU pragma: export
#include "damp/design/invariant_set.hpp"
#include "damp/design/options.hpp"
#include "damp/design/qp.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/status.hpp"
#include "damp/design/steady_state_target.hpp" // IWYU pragma: export
#include "damp/estimation/fixed_gain_observer.hpp"
#include "damp/estimation/kalman.hpp" // IWYU pragma: export
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/block.hpp"
#include "damp/matrix/decomposition.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/toolbox/bounds.hpp" // IWYU pragma: export

namespace damp {

namespace design {

namespace detail {

// filled_vector / convert_bounds (sentinel machinery) live in design/qp.hpp's
// detail — shared with the MHE and any other bound-carrying design function.

/** @brief Scale factor from a constraint range: (max − min) when both finite, else 1. */
template<typename T>
[[nodiscard]] constexpr T scale_from_range(T lo, T hi) {
    constexpr T thresh = damp::numeric_limits<T>::max() / T{4};
    if (lo <= -thresh || lo >= thresh || hi >= thresh || hi <= -thresh) {
        return T{1};
    }
    const T span = hi - lo;
    return (span > T{0}) ? span : T{1};
}

/** @brief Divide a bound by a scale factor, leaving unbounded sentinels untouched. */
template<typename T>
[[nodiscard]] constexpr T scaled_bound(T bound, T scale) {
    constexpr T thresh = damp::numeric_limits<T>::max() / T{4};
    if (bound >= thresh || bound <= -thresh) {
        return bound;
    }
    return bound / scale;
}

/** @brief The terminal ingredient a design::mpc call was given, as its result records it (internal: the options are DareCost and InvariantSet). */
enum class MPCTerminal : uint8_t {
    none = 0,      ///< No terminal option: the plain Qy / Qy_terminal cost, no extra set
    dare_cost,     ///< DareCost: V_f, the Δu stage cost's cost-to-go on z, no extra set
    invariant_set, ///< InvariantSet: that cost and the target-augmented invariant set
};

} // namespace detail

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
 * bad_weights, as design::state_feedback does); the QP Hessian must come out positive definite, which holds
 * whenever RΔu ≻ 0 or Ru ≻ 0.
 *
 * Named, the weights take their sizes from their fields, NU and NY 1 unless a field says
 * otherwise: `MPCWeights{.Rdu = {1.6}}` is MPCWeights<1, 1>, `MPCWeights{.Qy = diag(20.0, 1.0)}`
 * MPCWeights<1, 2>, and InputLimits the same: `InputLimits{.u_lim = {.lower = 0.0, .upper = 40.0}}`.
 * Qy and Qy_terminal are on every output of the plant, the rows of C; a zero row and column
 * leave an output out of the cost.
 *
 * Ru pulls each input toward the steady-state target's input u_t, the steady input that holds
 * the reference (under the estimated or measured disturbance), so it damps the input without an
 * offset. The runtime's set_input_target replaces u_t with an operating point of its own, and
 * clear_input_target returns to u_t; neither is part of the cost's design.
 *
 * @note Mirrors the MATLAB® mpc object's Weights.OV / Weights.MVRate /
 *       Weights.MV / Weights.ECR.
 */
template<size_t NU = 1, size_t NY = 1, typename T = double>
struct MPCWeights {
    Matrix<NY, NY, T> Qy = Matrix<NY, NY, T>::identity();  ///< Output tracking weight, steps 1..NP−1 (per second from a continuous plant)
    Matrix<NU, NU, T> Rdu = Matrix<NU, NU, T>::identity(); ///< Move weight on Δu (on the rate Δu/Ts, per second, from a continuous plant)
    Matrix<NU, NU, T> Ru{};                                ///< Input-deviation weight on u − u_t, the target's input or the runtime's input target (default 0; per second from a continuous plant)
    Matrix<NY, NY, T> Qy_terminal = Qy;                    ///< Output weight on the final step NP, or the last period (default: Qy)
    T                 ecr_weight{static_cast<T>(1.0e5)};   ///< Soft-constraint slack penalty ρ_ε

    /** @brief The weights in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr MPCWeights<NU, NY, U> as() const {
        return MPCWeights<NU, NY, U>{
            Qy.template as<U>(),
            Rdu.template as<U>(),
            Ru.template as<U>(),
            Qy_terminal.template as<U>(),
            static_cast<U>(ecr_weight),
        };
    }
};

/**
 * @struct OutputLimits
 * @brief Limits on an MPC's outputs: among design::mpc's options
 *
 * lim bounds the NYC outputs listed in outputs over the whole prediction horizon, and the QP
 * carries rows for those alone. The indices are the plant's outputs, the rows of C the QP
 * predicts and Qy weights, as a TrackedOutputs' and a disturbance spec's are.
 * `OutputLimits{.outputs = {1, 2}, .lim = {.lower = {-0.4, -0.26}, .upper = {0.4, 0.26}}}` bounds
 * outputs 1 and 2 and leaves the others free; outputs is an array, so the limits take
 * NYC from its list. Each listed output's scale factor is its range (max − min where both are
 * finite, else 1).
 *
 * Soft by default (ECR coefficient 1): lower_ecr(i) is the softness of the lim.lower(i) rows and
 * upper_ecr(i) that of the lim.upper(i) rows. Set one to 0 to make that side hard. A larger ECR
 * relaxes that bound more readily relative to the others. A listed output left unbounded has its
 * rows all the same, for set_limits to move at run time (MPCLimits::y_lim, in the order of
 * outputs). design::mpc refuses an index past the plant's outputs or listed twice, a NaN limit
 * or a lower limit above its upper one (invalid_parameter), and a negative or non-finite ECR
 * (bad_weights).
 *
 * @note Mirrors the MATLAB® mpc object's OutputVariables Min/Max (+ MinECR/MaxECR).
 */
template<size_t NYC, typename T = double>
struct OutputLimits {
    static_assert(NYC > 0, "OutputLimits lists the outputs it bounds: .outputs = {i, ...}, at least one; leave it out for none.");

    size_t         outputs[NYC]{};                              ///< Indices of the plant's outputs bounded, each once (an array, so bare braces give NYC)
    Bounds<NYC, T> lim{};                                       ///< Their bounds (soft through lower_ecr, upper_ecr)
    ColVec<NYC, T> lower_ecr{detail::filled_vector<NYC>(T{1})}; ///< Softness of the lim.lower rows, ≥ 0 (0 = hard; MATLAB's MinECR)
    ColVec<NYC, T> upper_ecr{detail::filled_vector<NYC>(T{1})}; ///< Softness of the lim.upper rows, ≥ 0 (0 = hard; MATLAB's MaxECR)
};

/**
 * @struct DareCost
 * @brief A terminal cost for design::mpc, among its options: at most one terminal ingredient
 *
 * Replaces the last step's Qy_terminal with V_f = z_Nᵀ P z_N, z = [x − x_t; u′_prev − u′_t], the
 * cost-to-go of the MPC's own Δu stage cost: a DARE on z with the move as input (needs
 * R_Δu + R_u positive definite), so V_f is a control Lyapunov function of the Δu form. Its
 * stability argument rests on shifting the last plan and appending the terminal move, which needs
 * a free move at every step: it holds with NC = NP. Cheap enough for constinit. `DareCost{}`.
 */
struct DareCost {};

/**
 * @struct InvariantSet
 * @brief A terminal cost and set for design::mpc, among its options: DareCost's cost and the
 *        maximal positively invariant set of its terminal law (at most one terminal ingredient)
 *
 * The set of the terminal law Δu′ = −K·z, augmented with the constant target (MPC for tracking),
 * as hard rows on z_N: one set valid for every target that keeps target_margin of each
 * half-range clear of its bounds (the InputLimits and the OutputLimits), the target moving each
 * row's offset at run time. `InvariantSet{.target_margin = 0.05}`. That set has more facets than
 * a zero-target one (a 2-state double integrator ~48): give it room, mpc<NP, NC, NT_max>, and
 * design it at initialization, since constant-evaluated it exceeds GCC's default
 * -fconstexpr-ops-limit.
 *
 * The set's rows are hard and the target is the steady state of r[0], so the QP is feasible only
 * from states that reach the set within NP steps for that target: a reference step outside that
 * region makes every tick Infeasible, and the runtime holds its input. Size NP for the reference
 * range. design::mpc refuses NC < NP (the shifted plan needs a free move at every step), a margin
 * outside (0, 1) and a StateLimits, whose rows the set does not cover (invalid_parameter), and
 * NT_max = 0 (bad_weights). The runtime's set_limits is refused, the set being built for the
 * design's limits.
 */
template<typename T = double>
struct InvariantSet {
    T target_margin{static_cast<T>(0.05)}; ///< The fraction of each half-range a tracked target keeps clear of its bounds, in (0, 1)
};

/**
 * @struct QPIterations
 * @brief The QP's iteration budget per tick, among design::mpc's options (left out:
 *        10·(NZ + 1 + NI + 1), NZ the moves and NI the constraint rows)
 *
 * `QPIterations{40}` caps each tick's active-set solve at 40 iterations, which bounds the tick's
 * worst-case time. A solve that runs out reports MaxIterations through last_status(), and the
 * runtime applies or holds as damp::MPC describes. design::mpc refuses a budget of 0
 * (invalid_parameter).
 */
struct QPIterations {
    size_t budget{0}; ///< Active-set iterations per tick, ≥ 1
};

/**
 * @struct StateLimits
 * @brief Limits on states an MPC does not track: among design::mpc's options
 *
 * lim bounds the NXC states of the prediction model listed in states (a delayed model's delay
 * states come after the plant's states), over the whole prediction horizon, soft by default like
 * OutputLimits: a lean, a speed or a temperature the law must keep within without
 * tracking it. `StateLimits{.states = {2}, .lim = {.lower = -0.26, .upper = 0.26}}` bounds state 2;
 * states is an array, so the limits take NXC from its list. lower_ecr(i) and upper_ecr(i) are the
 * softness of each side (0 = hard). design::mpc refuses a state index past the model's or listed
 * twice, a NaN limit or a lower limit above its upper one, and state limits with an
 * InvariantSet, whose rows do not cover the states (invalid_parameter); a negative
 * or non-finite ECR is refused with bad_weights. The runtime moves lim with set_limits (MPCLimits::x_lim).
 */
template<size_t NXC, typename T = double>
struct StateLimits {
    static_assert(NXC > 0, "StateLimits lists the states it bounds: .states = {i, ...}, at least one.");

    size_t         states[NXC]{};                               ///< Indices of the prediction model's states bounded, each once (an array, so bare braces give NXC)
    Bounds<NXC, T> lim{};                                       ///< Their bounds (soft through lower_ecr, upper_ecr)
    ColVec<NXC, T> lower_ecr{detail::filled_vector<NXC>(T{1})}; ///< Softness of the lim.lower rows, ≥ 0 (0 = hard)
    ColVec<NXC, T> upper_ecr{detail::filled_vector<NXC>(T{1})}; ///< Softness of the lim.upper rows, ≥ 0 (0 = hard)
};

namespace detail {

/**
 * @brief The options design::mpc's QP reads (InputLimits, OutputLimits, StateLimits, the terminal
 *        ingredient, QPIterations) as the design reads them: NYC bounded outputs and
 *        NXC bounded states (0: none).
 */
template<size_t NU, size_t NYC = 0, typename T = double, size_t NXC = 0>
struct MPCConstraintSet {
    Bounds<NU, T>            u_lim{};                                       ///< Input bounds (hard)
    Bounds<NU, T>            du_lim{};                                      ///< Move bounds on Δu (hard)
    damp::array<size_t, NYC> y_outputs{};                                   ///< The bounded outputs (OutputLimits::outputs)
    Bounds<NYC, T>           y_lim{};                                       ///< Their bounds (OutputLimits::lim)
    ColVec<NYC, T>           y_lower_ecr{detail::filled_vector<NYC>(T{1})}; ///< Softness of the y_lim.lower rows
    ColVec<NYC, T>           y_upper_ecr{detail::filled_vector<NYC>(T{1})}; ///< Softness of the y_lim.upper rows
    damp::array<size_t, NXC> x_states{};                                    ///< The bounded states (StateLimits::states)
    Bounds<NXC, T>           x_lim{};                                       ///< Their bounds (StateLimits::lim)
    ColVec<NXC, T>           x_lower_ecr{detail::filled_vector<NXC>(T{1})}; ///< Softness of the x_lim.lower rows
    ColVec<NXC, T>           x_upper_ecr{detail::filled_vector<NXC>(T{1})}; ///< Softness of the x_lim.upper rows
    MPCTerminal              terminal{MPCTerminal::none};                   ///< The terminal option given (DareCost, InvariantSet; none)
    T                        terminal_target_margin{static_cast<T>(0.05)};  ///< InvariantSet::target_margin
    damp::optional<size_t>   max_qp_iterations{};                           ///< QPIterations::budget (none: the default)

    /** @brief The constraints in another scalar. */
    template<typename U>
    [[nodiscard]] constexpr MPCConstraintSet<NU, NYC, U, NXC> as() const {
        return MPCConstraintSet<NU, NYC, U, NXC>{
            u_lim.template as<U>(),
            du_lim.template as<U>(),
            y_outputs,
            y_lim.template as<U>(),
            ColVec<NYC, U>{y_lower_ecr.template as<U>()},
            ColVec<NYC, U>{y_upper_ecr.template as<U>()},
            x_states,
            x_lim.template as<U>(),
            ColVec<NXC, U>{x_lower_ecr.template as<U>()},
            ColVec<NXC, U>{x_upper_ecr.template as<U>()},
            terminal,
            static_cast<U>(terminal_target_margin),
            max_qp_iterations,
        };
    }
};

/** @brief An OutputLimits: a kind of design::mpc option. */
template<typename O>
struct is_output_limits : damp::false_type {};

/** @brief OutputLimits of any length. */
template<size_t NYC, typename T>
struct is_output_limits<OutputLimits<NYC, T>> : damp::true_type {};

/** @brief A StateLimits: a kind of design::mpc option. */
template<typename O>
struct is_state_limits : damp::false_type {};

/** @brief StateLimits of any length. */
template<size_t NXC, typename T>
struct is_state_limits<StateLimits<NXC, T>> : damp::true_type {};

/** @brief A terminal ingredient: a kind of design::mpc option, DareCost or InvariantSet. */
template<typename O>
struct is_mpc_terminal : damp::false_type {};

/** @brief DareCost. */
template<>
struct is_mpc_terminal<DareCost> : damp::true_type {};

/** @brief InvariantSet. */
template<typename T>
struct is_mpc_terminal<InvariantSet<T>> : damp::true_type {};

/** @brief A QPIterations: a kind of design::mpc option. */
template<typename O>
struct is_qp_iterations : damp::false_type {};

/** @brief QPIterations. */
template<>
struct is_qp_iterations<QPIterations> : damp::true_type {};

/** @brief How many entries a limits option lists: none without one. */
template<typename O>
inline constexpr size_t limits_listed = 0;

/** @brief An OutputLimits lists NYC outputs. */
template<size_t NYC, typename T>
inline constexpr size_t limits_listed<OutputLimits<NYC, T>> = NYC;

/** @brief A StateLimits lists NXC states. */
template<size_t NXC, typename T>
inline constexpr size_t limits_listed<StateLimits<NXC, T>> = NXC;

/** @brief Copy a limits option's index list into the set's array. */
template<size_t N>
constexpr void copy_listed(damp::array<size_t, N>& to, const size_t (&from)[N]) {
    for (size_t i = 0; i < N; ++i) {
        to[i] = from[i];
    }
}

/**
 * @brief The options among @p opts that the QP reads, as the design's set on the plant's NU inputs
 *        in its scalar @p T: the InputLimits (none: unbounded), the OutputLimits (none: no bounded
 *        output), the StateLimits (none: no bounded state), the terminal ingredient (none: the
 *        plain cost) and the QP budget (none: the default).
 */
template<size_t NU, typename T, typename... Opts>
[[nodiscard]] constexpr auto constraint_set(const Opts&... opts) {
    constexpr size_t                  NYC = limits_listed<::damp::detail::option_t<is_output_limits, void, Opts...>>;
    constexpr size_t                  NXC = limits_listed<::damp::detail::option_t<is_state_limits, void, Opts...>>;
    const auto                        input = input_limits_of<NU, T>(opts...);
    MPCConstraintSet<NU, NYC, T, NXC> out{.u_lim = input.u_lim, .du_lim = input.du_lim};
    if constexpr (NYC > 0) {
        const auto limits = ::damp::detail::option<is_output_limits>(0, opts...);
        copy_listed(out.y_outputs, limits.outputs);
        out.y_lim = limits.lim.template as<T>();
        out.y_lower_ecr = ColVec<NYC, T>{limits.lower_ecr.template as<T>()};
        out.y_upper_ecr = ColVec<NYC, T>{limits.upper_ecr.template as<T>()};
    }
    if constexpr (NXC > 0) {
        const auto limits = ::damp::detail::option<is_state_limits>(0, opts...);
        copy_listed(out.x_states, limits.states);
        out.x_lim = limits.lim.template as<T>();
        out.x_lower_ecr = ColVec<NXC, T>{limits.lower_ecr.template as<T>()};
        out.x_upper_ecr = ColVec<NXC, T>{limits.upper_ecr.template as<T>()};
    }
    using Terminal = ::damp::detail::option_t<is_mpc_terminal, void, Opts...>;
    if constexpr (damp::is_same_v<Terminal, DareCost>) {
        out.terminal = MPCTerminal::dare_cost;
    } else if constexpr (!damp::is_same_v<Terminal, void>) {
        out.terminal = MPCTerminal::invariant_set;
        out.terminal_target_margin = static_cast<T>(::damp::detail::option<is_mpc_terminal>(0, opts...).target_margin);
    }
    if constexpr (::damp::detail::option_count<is_qp_iterations, Opts...> > 0) {
        out.max_qp_iterations = ::damp::detail::option<is_qp_iterations>(0, opts...).budget;
    }
    return out;
}

} // namespace detail

/**
 * @brief The limits an MPC takes at run time (MATLAB®'s online constraints): u, Δu, the
 *        outputs its design's OutputLimits lists and the states its StateLimits lists,
 *        in physical units. A field left out is unbounded.
 *
 * The softness (ECR), the terminal ingredient and the QP budget stay as designed. `auto l =
 * mpc.limits(); l.u_lim.upper = …; mpc.set_limits(l);` changes one. y_lim(i) and x_lim(i) are
 * the bounds of the i-th output and state their option lists.
 */
template<size_t NU, size_t NYC = 0, typename T = double, size_t NXC = 0>
struct MPCLimits {
    Bounds<NU, T>  u_lim{};  ///< Input bounds
    Bounds<NU, T>  du_lim{}; ///< Move bounds on Δu
    Bounds<NYC, T> y_lim{};  ///< Bounds on the outputs the design's OutputLimits lists
    Bounds<NXC, T> x_lim{};  ///< Bounds on the states the design's StateLimits lists
};

/**
 * @struct MPCArtifacts
 * @brief Condensed-QP data produced by mpc() (MPCResult::mpc), consumed by damp::MPC
 *
 * Holds the dense QP in [Z; ε] (move increments + soft-constraint slack):
 * Hessian H, gradient maps (f = F_x·x̃0 − F_r·y'_r − F_ut·u'_t with
 * x̃0 = [x; u'_prev; d], y'_r the reference on the outputs at each step and u'_t the input
 * target at each move), the stacked constraint matrix A_con with its
 * state-dependent bounds (b = b_bound + B_map·x̃0), the per-channel scale
 * factors, and the physical clamping limits. Primes denote internally scaled
 * quantities. Use .as<float>() for embedded deployment (bound sentinels
 * convert safely).
 *
 * @tparam NX       States of the prediction model
 * @tparam NY       Outputs the QP predicts and weighs: every row of the plant's C
 * @tparam ND       Disturbance channels: the measured disturbances, or the estimated ones
 * @tparam NXC      States with bounds (the design's StateLimits; 0: none)
 * @tparam NYC      Outputs with bounds (the design's OutputLimits; 0: none)
 */
template<size_t NX, size_t NU, size_t NY, size_t NP, size_t NC, size_t ND = 0, typename T = double, size_t NT_max = 0, size_t NXC = 0, size_t NYC = 0>
struct [[nodiscard]] MPCArtifacts {
    static constexpr size_t NXA = NX + NU + ND;                         ///< Augmented state size [x; u'_prev; d]
    static constexpr size_t NYS = NP * NY;                              ///< Stacked outputs over the prediction horizon (the reference)
    static constexpr size_t NYCS = NP * NYC;                            ///< Stacked bounded outputs over the prediction horizon
    static constexpr size_t NXS = NP * NXC;                             ///< Stacked bounded states over the prediction horizon
    static constexpr size_t NZ = NC * NU;                               ///< Move-increment decision variables
    static constexpr size_t NZ1 = NZ + 1;                               ///< QP decision size (moves + slack ε)
    static constexpr size_t NI_box = (4 * NZ) + (2 * NYCS) + (2 * NXS); ///< Rows before the terminal block: Δu, u, the bounded y, x
    static constexpr size_t NI = NI_box + NT_max;                       ///< Box-constraint rows (+ terminal capacity)
    static constexpr size_t NI1 = NI + 1;                               ///< QP constraint rows (+ ε ≥ 0)

    Matrix<NZ1, NZ1, T>                       H{};                                 ///< QP Hessian (positive definite)
    Matrix<NZ, NXA, T>                        F_x{};                               ///< Gradient state map
    Matrix<NZ, NYS, T>                        F_r{};                               ///< Gradient reference map, one scaled reference per prediction step (preview)
    Matrix<NZ, NZ, T>                         F_ut{};                              ///< Gradient input-target map, one scaled target per move
    Matrix<NI1, NZ1, T>                       A_con{};                             ///< Constraint rows: A_con·[Z;ε] ≤ b_bound + B_map·x̃0
    ColVec<NI1, T>                            b_bound{};                           ///< Constant bound term (sentinel rows disabled)
    Matrix<NI1, NXA, T>                       B_map{};                             ///< State-dependent bound term
    Matrix<NY, NX, T>                         C{};                                 ///< Output map of the prediction model (the target's output)
    Matrix<NY, ND, T>                         Cd{};                                ///< The disturbance channels into the outputs (C_{d,z}; zero for measured ones)
    ColVec<NU, T>                             scale_u{};                           ///< Input scale factors
    ColVec<NY, T>                             scale_y{};                           ///< Output scale factors (1 on an output the OutputLimits does not list)
    ColVec<NXC, T>                            scale_x{};                           ///< Scale factors of the bounded states
    detail::MPCConstraintSet<NU, NYC, T, NXC> constraints{};                       ///< Original limits (runtime clamp fallback)
    size_t                                    max_qp_iterations{};                 ///< Per-tick QP iteration budget
    Matrix<NX + NU, NX + NU, T>               P{};                                 ///< Terminal cost on z = [x − x_t; u′_prev − u′_t] (zero when unused)
    Matrix<NU, NX + NU, T>                    K_terminal{};                        ///< Terminal move law Δu′ = −K·z
    Matrix<NT_max, NX + NU, T>                G_terminal{};                        ///< Terminal-set rows on z; first n_terminal used
    Matrix<NT_max, NY + NU, T>                G_terminal_target{};                 ///< Their target terms: G·z + G_target·[y_t; u_t] ≤ h
    ColVec<NT_max, T>                         h_terminal{};                        ///< Terminal-set offsets
    size_t                                    n_terminal{0};                       ///< Live terminal-set rows
    detail::MPCTerminal                       terminal{detail::MPCTerminal::none}; ///< The terminal option it was given (internal: none, DareCost or InvariantSet)
    DesignStatus                              status{DesignStatus::not_converged}; ///< ok, or why the design failed
    bool                                      success{false};                      ///< status == DesignStatus::ok
    T                                         Ts{0};                               ///< Sample time of the law [s]

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The artifacts in another scalar; refused (singular) when H does not factor in it. */
    template<typename U>
    [[nodiscard]] constexpr MPCArtifacts<NX, NU, NY, NP, NC, ND, U, NT_max, NXC, NYC> as() const {
        MPCArtifacts<NX, NU, NY, NP, NC, ND, U, NT_max, NXC, NYC> out{
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
            .scale_x = ColVec<NXC, U>{scale_x.template as<U>()},
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
     * The inputs are the tick's arguments in order: r, the QP's reference on the NY outputs (the
     * tick's r when the target holds every output), the ND disturbance channels (none when
     * ND = 0), and the state x of the law's model the QP is solved from. Close it around a plant
     * that reports [0; 0; x] (the reference and the disturbance entering through the first rows)
     * for the loop's poles, margins and DC gain, as an LQ result's to_ss() is closed. The block is
     * the law with its input target at zero: Ru's default anchor is the steady-state target's
     * input, which these artifacts do not hold (it does not enter with Ru = 0).
     *
     * Exact without a terminal option and without Ru. A terminal cost and Ru are posed about the
     * target of the reference, which these artifacts do not hold: the block keeps the feedback from
     * x and u_prev (the loop's poles and margins) but leaves their pull to that target out of its
     * r and d gains. A design whose Hessian does not factor gives the zero block.
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

/** @brief A weight matrix design::mpc accepts: finite, symmetric and positive semidefinite (as design::state_feedback's). */
template<size_t N, typename T>
[[nodiscard]] constexpr bool mpc_weight_valid(const Matrix<N, N, T>& W) {
    return mat::is_symmetric_or_hermitian(W) && ::damp::detail::is_psd(W);
}

/** @brief The weights and softness design::mpc accepts: PSD weights, a finite positive ECR weight, finite ECR coefficients ≥ 0. */
template<size_t NU, size_t NY, typename T, size_t NYC, size_t NXC>
[[nodiscard]] constexpr bool mpc_weights_valid(const MPCWeights<NU, NY, T>& w, const MPCConstraintSet<NU, NYC, T, NXC>& c) {
    const auto softness = [](T v) { return damp::isfinite(v) && v >= T{0}; };
    return mpc_weight_valid(w.Qy) && mpc_weight_valid(w.Qy_terminal) && mpc_weight_valid(w.Rdu) && mpc_weight_valid(w.Ru)
        && damp::isfinite(w.ecr_weight) && w.ecr_weight > T{0} && mat::all_of(c.y_lower_ecr, softness) && mat::all_of(c.y_upper_ecr, softness)
        && mat::all_of(c.x_lower_ecr, softness) && mat::all_of(c.x_upper_ecr, softness);
}

/**
 * @brief Whether the QP weighs output @p r (a nonzero row of Qy or Qy_terminal) or bounds it
 *        (OutputLimits lists it). A NaN or ±Inf entry is nonzero, tested by its bits.
 */
template<size_t NU, size_t NY, typename T, size_t NYC, size_t NXC>
[[nodiscard]] constexpr bool mpc_output_used(const MPCWeights<NU, NY, T>& w, const MPCConstraintSet<NU, NYC, T, NXC>& c, size_t r) {
    const auto nonzero = [](T v) { return !damp::isfinite(v) || v != T{0}; };
    bool       used = false;
    for (size_t j = 0; j < NY; ++j) {
        used = used || nonzero(w.Qy(r, j)) || nonzero(w.Qy_terminal(r, j));
    }
    for (size_t k = 0; k < NYC; ++k) {
        used = used || c.y_outputs[k] == r;
    }
    return used;
}

/** @brief Each index of @p list is below @p N and listed once. */
template<size_t N, size_t M>
[[nodiscard]] constexpr bool mpc_indices_valid(const damp::array<size_t, M>& list) {
    for (size_t i = 0; i < M; ++i) {
        if (list[i] >= N) {
            return false;
        }
        for (size_t j = 0; j < i; ++j) {
            if (list[j] == list[i]) {
                return false;
            }
        }
    }
    return true;
}

/** @brief The bounded states are the model's and the bounded outputs the plant's, each listed once, inside valid boxes. */
template<size_t NX, size_t NY, size_t NU, size_t NYC, typename T, size_t NXC>
[[nodiscard]] constexpr bool mpc_signal_bounds_valid(const MPCConstraintSet<NU, NYC, T, NXC>& c) {
    return mpc_indices_valid<NX>(c.x_states) && c.x_lim.valid() && mpc_indices_valid<NY>(c.y_outputs) && c.y_lim.valid();
}

/**
 * @brief Write the box rows of the QP's bound vector from @p lim with the design's scales: Δu and
 *        u over the control horizon, the bounded outputs and states over the
 *        prediction horizon.
 */
template<size_t NX, size_t NU, size_t NY, size_t NP, size_t NC, size_t ND, typename T, size_t NT_max, size_t NXC, size_t NYC>
constexpr void mpc_box_bounds(MPCArtifacts<NX, NU, NY, NP, NC, ND, T, NT_max, NXC, NYC>& art, const MPCLimits<NU, NYC, T, NXC>& lim) {
    constexpr size_t NZ = NC * NU;
    constexpr size_t NYCS = NP * NYC;
    constexpr size_t r_du_hi = 0;
    constexpr size_t r_du_lo = NZ;
    constexpr size_t r_u_hi = 2 * NZ;
    constexpr size_t r_u_lo = 3 * NZ;
    constexpr size_t r_y_hi = 4 * NZ;
    constexpr size_t r_y_lo = (4 * NZ) + NYCS;
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
        for (size_t j = 0; j < NYC; ++j) {
            const T sy = art.scale_y(art.constraints.y_outputs[j]);
            art.b_bound(r_y_hi + (kb * NYC) + j) = scaled_bound(lim.y_lim.upper(j), sy);
            art.b_bound(r_y_lo + (kb * NYC) + j) = scaled_bound(-lim.y_lim.lower(j), sy);
        }
    }
    constexpr size_t NXS = NP * NXC;
    constexpr size_t r_x_hi = (4 * NZ) + (2 * NYCS);
    constexpr size_t r_x_lo = r_x_hi + NXS;
    for (size_t kb = 0; kb < NP; ++kb) {
        for (size_t j = 0; j < NXC; ++j) {
            const T sx = art.scale_x(j);
            art.b_bound(r_x_hi + (kb * NXC) + j) = scaled_bound(lim.x_lim.upper(j), sx);
            art.b_bound(r_x_lo + (kb * NXC) + j) = scaled_bound(-lim.x_lim.lower(j), sx);
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
 * that makes the set finitely determined with the constant target's unit eigenvalue). The output
 * rows are those of the bounded outputs (OutputLimits), in the order it lists them.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NYC, size_t NXC, size_t NA = NX + NU + NY + NU, size_t NG = (4 * NYC) + (6 * NU)>
constexpr void mpc_terminal_rows(
    const Matrix<NY, NX, T>&                 C,
    const Matrix<NU, NX + NU, T>&            K,
    const ColVec<NU, T>&                     scale_u,
    const ColVec<NY, T>&                     scale_y,
    const MPCConstraintSet<NU, NYC, T, NXC>& c,
    Matrix<NG, NA, T>&                       G0,
    ColVec<NG, T>&                           h0
) {
    constexpr size_t NZD = NX + NU;
    h0 = filled_vector<NG>(unbounded_bound<T>());
    size_t     n = 0;
    const T    half = static_cast<T>(0.5) * c.terminal_target_margin;
    const auto push_pair = [&](const ColVec<NA, T>& a, T hi, T lo) {
        mpc_push_box_row(G0, h0, n, a, hi);
        mpc_push_box_row(G0, h0, n, ColVec<NA, T>(a * T{-1}), -lo);
    };
    for (size_t k = 0; k < NYC; ++k) {
        const size_t  i = c.y_outputs[k];
        ColVec<NA, T> a{};
        for (size_t j = 0; j < NX; ++j) {
            a(j) = C(i, j);
        }
        a(NZD + i) = T{1};
        push_pair(a, c.y_lim.upper(k), c.y_lim.lower(k));
        ColVec<NA, T> t{};
        t(NZD + i) = T{1};
        const T m = half * scale_y(i);
        push_pair(t, c.y_lim.upper(k) - m, c.y_lim.lower(k) + m);
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

template<size_t NP, size_t NC, size_t ND, size_t NT_max, size_t NX, size_t NU, size_t NY, typename T, size_t NXC, size_t NYC>
constexpr void mpc_stack_terminal_rows(
    MPCArtifacts<NX, NU, NY, NP, NC, ND, T, NT_max, NXC, NYC>& art,
    const Matrix<NX + NU, NX + NU + ND, T>&                    Phi_z,
    const Matrix<NX + NU, NC * NU, T>&                         Gamma_z
) {
    constexpr size_t NZ = NC * NU;
    constexpr size_t NXA = NX + NU + ND;
    constexpr size_t r_term = MPCArtifacts<NX, NU, NY, NP, NC, ND, T, NT_max, NXC, NYC>::NI_box;
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
 *        outputs, the disturbance channel and the period. Off on a discrete plant.
 */
template<size_t NXP, size_t NU, size_t NY, size_t ND, typename T>
struct MpcContinuousCost {
    Matrix<NXP, NXP, T> A{};       ///< Continuous state matrix
    Matrix<NXP, NU, T>  B{};       ///< Continuous input matrix
    Matrix<NXP, ND, T>  Bd{};      ///< Continuous disturbance channel
    Matrix<NY, NXP, T>  C{};       ///< The outputs the cost is on
    Matrix<NY, NU, T>   D{};       ///< Their feedthrough
    T                   Ts{0};     ///< Period [s]
    bool                on{false}; ///< Whether the design is from the continuous plant

    damp::array<::damp::detail::RegisterDelay<T>, NU> inputs{}; ///< Each input's delay, d·Ts + τ (none: zero)
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
    //! The scaled output error e = E·p + Ds·u: E = [S_y⁻¹C, 0, −I], Ds = S_y⁻¹D
    Matrix<NY, NPS, T> E{};
    Matrix<NY, NU, T>  Ds{};
    for (size_t i = 0; i < NY; ++i) {
        for (size_t j = 0; j < NXP; ++j) {
            E(i, j) = plant.C(i, j) / scale_y(i);
        }
        E(i, NXP + ND + i) = T{-1};
        for (size_t j = 0; j < NU; ++j) {
            Ds(i, j) = plant.D(i, j) / scale_y(i);
        }
    }
    const Matrix<NPS, NY, T>  EtQ = E.transpose() * Qy;
    const LQRCost<NPS, NU, T> integrand{.Q = EtQ * E, .R = Ds.transpose() * (Qy * Ds), .N = EtQ * Ds};
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

/**
 * @brief mpc_period_cost on a plant whose input arrives late: the output cost of one period as a
 *        quadratic form on s = [z; u′; d; r′], z = [x; v(k−1); …; v(k−N)] the delayed model's state
 *
 * Within the period input j holds v_j(k − d_j − 1) until its switch at τ_j, then v_j(k − d_j), the
 * command u′ itself when d_j = 0 (discretize with an InputDelay). The switches cut the period
 * into stretches over which every input holds one value. Over each, Van Loan gives the exact
 * cost of [x; d; r′] from the stretch's start under the held value, and the exponential carries
 * the start to the next stretch. The start and the held values are linear maps of s, so the
 * stretches sum to a quadratic form on it, as delayed_cost does for the LQ designs.
 */
template<size_t NX, size_t NXP, size_t NU, size_t NY, size_t ND, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NX + NU + ND + NY, NX + NU + ND + NY, T>> mpc_delayed_period_cost(
    const MpcContinuousCost<NXP, NU, NY, ND, T>& plant,
    const Matrix<NY, NY, T>&                     Qy,
    const ColVec<NU, T>&                         scale_u,
    const ColVec<NY, T>&                         scale_y
) {
    constexpr size_t    NPS = NXP + ND + NY;    // p = [x; d; r′], d and r′ held
    constexpr size_t    NW = NPS + NU;          // w = [p; v], v the input the plant holds
    constexpr size_t    NS = NX + NU + ND + NY; // s = [z; u′; d; r′]
    Matrix<NPS, NPS, T> Ap{};
    Ap.template block<NXP, NXP>(0, 0) = plant.A;
    Ap.template block<NXP, ND>(0, NXP) = plant.Bd;
    Matrix<NPS, NU, T> Bp{};
    Bp.template block<NXP, NU>(0, 0) = plant.B;
    Matrix<NY, NPS, T> E{};
    Matrix<NY, NU, T>  Ds{};
    for (size_t i = 0; i < NY; ++i) {
        for (size_t j = 0; j < NXP; ++j) {
            E(i, j) = plant.C(i, j) / scale_y(i);
        }
        E(i, NXP + ND + i) = T{-1};
        for (size_t j = 0; j < NU; ++j) {
            Ds(i, j) = plant.D(i, j) / scale_y(i);
        }
    }
    const Matrix<NPS, NY, T>  EtQ = E.transpose() * Qy;
    const Matrix<NPS, NPS, T> Qp = EtQ * E;
    const Matrix<NU, NU, T>   Rp = Ds.transpose() * (Qy * Ds);
    const Matrix<NPS, NU, T>  Np = EtQ * Ds;
    Matrix<NW, NW, T>         flow{}; // d/dt [p; v] = [Ap·p + Bp·v; 0]
    flow.template block<NPS, NPS>(0, 0) = Ap;
    flow.template block<NPS, NU>(0, NPS) = Bp;

    //! The stretches: 0, each switch once in time order, then Ts.
    damp::array<T, NU + 2> edge{};
    size_t                 cuts = 1;
    for (size_t j = 0; j < NU; ++j) {
        const T tau = plant.inputs[j].fraction;
        bool    seen = !(tau > T{0});
        for (size_t c = 1; c < cuts; ++c) {
            seen = seen || (edge[c] == tau);
        }
        if (!seen) {
            size_t c = cuts++;
            for (; c > 1 && edge[c - 1] > tau; --c) {
                edge[c] = edge[c - 1];
            }
            edge[c] = tau;
        }
    }
    edge[cuts] = plant.Ts;

    //! p at the period's start as a map of s: the plant's states of z, d and r′
    Matrix<NPS, NS, T> start{};
    for (size_t i = 0; i < NXP; ++i) {
        start(i, i) = T{1};
    }
    for (size_t j = 0; j < ND + NY; ++j) {
        start(NXP + j, NX + NU + j) = T{1};
    }
    Matrix<NS, NS, T> total{};
    for (size_t c = 0; c < cuts; ++c) {
        const T a = edge[c];
        const T h = edge[c + 1] - a;
        if (!(h > T{0})) {
            continue;
        }
        Matrix<NW, NS, T> M{}; // [p; v] at the stretch's start
        M.template block<NPS, NS>(0, 0) = start;
        for (size_t j = 0; j < NU; ++j) {
            const bool   before = plant.inputs[j].fraction > T{0} && a < plant.inputs[j].fraction;
            const size_t age = plant.inputs[j].whole + (before ? 1 : 0);
            if (age == 0) {
                M(NPS + j, NX + j) = scale_u(j); // the command itself, u = S_u·u′
            } else {
                M(NPS + j, NXP + ((age - 1) * NU) + j) = T{1}; // a delay state, in physical units
            }
        }
        const auto c_opt = discretize_lqr_cost(Ap, Bp, Qp, Rp, h, Np);
        const auto F = mat::expm(Matrix<NW, NW, T>(flow * h));
        if (!c_opt || !F) {
            return damp::nullopt;
        }
        Matrix<NW, NW, T> Wp{};
        Wp.template block<NPS, NPS>(0, 0) = c_opt->Q;
        Wp.template block<NPS, NU>(0, NPS) = c_opt->N;
        Wp.template block<NU, NPS>(NPS, 0) = c_opt->N.transpose();
        Wp.template block<NU, NU>(NPS, NPS) = c_opt->R;
        total = total + (M.transpose() * (Wp * M));
        start = Matrix<NPS, NW, T>(F->template block<NPS, NW>(0, 0)) * M; // p at the next stretch's start
    }
    return (total + total.transpose()) * static_cast<T>(0.5);
}

/** @brief The output cost of one period on the law's model of @p NX states: the plant's own (mpc_period_cost), or the delayed one when the model carries delay states. */
template<size_t NX, size_t NXP, size_t NU, size_t NY, size_t ND, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NX + NU + ND + NY, NX + NU + ND + NY, T>> mpc_model_period_cost(
    const MpcContinuousCost<NXP, NU, NY, ND, T>& plant,
    const Matrix<NY, NY, T>&                     Qy,
    const ColVec<NU, T>&                         scale_u,
    const ColVec<NY, T>&                         scale_y
) {
    if constexpr (NX == NXP) {
        return mpc_period_cost(plant, Qy, scale_u, scale_y);
    } else {
        static_assert(NX > NXP && (NX - NXP) % NU == 0, "The law's model is the plant, then the delay states of its inputs.");
        return mpc_delayed_period_cost<NX>(plant, Qy, scale_u, scale_y);
    }
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
template<size_t NP, size_t NC, size_t NXP, size_t NX, size_t NU, size_t NY, size_t ND, size_t NT_max, typename T, size_t NXC, size_t NYC, size_t NXA = NX + NU + ND, size_t NZ = NC * NU>
[[nodiscard]] constexpr bool mpc_output_cost(
    MPCArtifacts<NX, NU, NY, NP, NC, ND, T, NT_max, NXC, NYC>& art,
    const MpcPred<NP, NC, NXA, NX, NU, NY, T>&                 pred,
    const Matrix<NXA, NXA, T>&                                 Aa,
    const Matrix<NXA, NU, T>&                                  Ba,
    const MPCWeights<NU, NY, T>&                               weights,
    const MpcContinuousCost<NXP, NU, NY, ND, T>&               plant,
    bool                                                       dare_on,
    Matrix<NZ, NZ, T>&                                         Hm,
    Matrix<NZ, NXA, T>&                                        Fx
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
    const auto W_opt = mpc_model_period_cost<NX>(plant, weights.Qy, art.scale_u, art.scale_y);
    const auto W_last_opt = dare_on ? W_opt : mpc_model_period_cost<NX>(plant, weights.Qy_terminal, art.scale_u, art.scale_y);
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

template<size_t NP, size_t NC, size_t ND, size_t NT_max = 0, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NXP, size_t NYC, size_t NXC>
[[nodiscard]] constexpr MPCArtifacts<NX, NU, NY, NP, NC, ND, T, NT_max, NXC, NYC> mpc_impl(
    const StateSpace<NX, NU, NY, T, NW, NV>&     sys,
    const Matrix<NX, ND, T>&                     Bd,
    const Matrix<NY, ND, T>&                     Cd,
    const MPCWeights<NU, NY, T>&                 weights,
    const MPCConstraintSet<NU, NYC, T, NXC>&     constraints,
    const MpcContinuousCost<NXP, NU, NY, ND, T>& plant
) {
    static_assert(NP >= 1, "Prediction horizon must be at least one step.");
    static_assert(NC >= 1 && NC <= NP, "Control horizon must satisfy 1 <= NC <= NP.");

    constexpr size_t NXA = NX + NU + ND;
    constexpr size_t NZ = NC * NU;
    constexpr size_t NZ1 = NZ + 1;
    constexpr size_t NYCS = NP * NYC;
    constexpr size_t NXS = NP * NXC;

    using Art = MPCArtifacts<NX, NU, NY, NP, NC, ND, T, NT_max, NXC, NYC>;
    constexpr size_t NI = Art::NI;
    Art              art{};
    art.constraints = constraints;
    art.C = sys.C;
    art.Cd = Cd;
    art.Ts = sys.Ts;
    art.terminal = constraints.terminal;
    art.max_qp_iterations = constraints.max_qp_iterations.value_or(10 * (NZ1 + NI + 1));
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
            if ((!damp::isfinite(sys.D(r, c)) || sys.D(r, c) != T{0}) && mpc_output_used(weights, constraints, r)) {
                return ::damp::detail::set_status(art, DesignStatus::bad_weights); // the prediction is y = C·x: no feedthrough (a NaN one either) on an output the QP weighs or bounds
            }
        }
    }
    if (!mpc_weights_valid(weights, constraints)) {
        return ::damp::detail::set_status(art, DesignStatus::bad_weights); // as design::state_feedback: indefinite, asymmetric or non-finite
    }
    if (!constraints.u_lim.valid() || !constraints.du_lim.valid()) {
        return ::damp::detail::set_status(art, DesignStatus::invalid_parameter); // as set_limits: a NaN or inverted limit
    }
    if (!mpc_signal_bounds_valid<NX, NY>(constraints) || (NXC > 0 && constraints.terminal == MPCTerminal::invariant_set)) {
        return ::damp::detail::set_status(art, DesignStatus::invalid_parameter); // an index past the model's or listed twice, a bad box, or a set that does not cover the states
    }
    if (constraints.terminal == MPCTerminal::invariant_set && NC < NP) {
        return ::damp::detail::set_status(art, DesignStatus::invalid_parameter); // the shifted plan needs a free move at every step
    }
    if (art.max_qp_iterations == 0) {
        return ::damp::detail::set_status(art, DesignStatus::invalid_parameter); // QPIterations{0}: not one active-set step
    }

    for (size_t j = 0; j < NU; ++j) {
        art.scale_u(j) = detail::scale_from_range(constraints.u_lim.lower(j), constraints.u_lim.upper(j));
    }
    for (size_t i = 0; i < NY; ++i) {
        art.scale_y(i) = T{1}; // an output no OutputLimits lists has no range
    }
    for (size_t k = 0; k < NYC; ++k) {
        art.scale_y(constraints.y_outputs[k]) = detail::scale_from_range(constraints.y_lim.lower(k), constraints.y_lim.upper(k));
    }
    for (size_t i = 0; i < NXC; ++i) {
        art.scale_x(i) = detail::scale_from_range(constraints.x_lim.lower(i), constraints.x_lim.upper(i));
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
        for (size_t i = 0; i < NY; ++i) { // the outputs carry an estimated output bias, so the output limits see it
            for (size_t j = 0; j < ND; ++j) {
                Ca(i, NX + NU + j) = Cd(i, j) / art.scale_y(i);
            }
        }
    }

    const auto pred = mpc_predictions<NP, NC, NXA, NX, NU, NY, T>(Aa, Ba, Ca);
    //! A terminal cost weights the state at NP, delay states included, so it sees past the horizon.
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
        const auto stage = plant.on ? mpc_stage_of<NX, NU>(*mpc_model_period_cost<NX>(plant, weights.Qy, art.scale_u, art.scale_y)) // formed above
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
    constexpr size_t r_y_lo = (4 * NZ) + NYCS;

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
    //! The bounded outputs over the horizon: their rows of the output prediction, step by step.
    for (size_t kb = 0; kb < NP; ++kb) {
        for (size_t j = 0; j < NYC; ++j) {
            const size_t hi = r_y_hi + (kb * NYC) + j;
            const size_t lo = r_y_lo + (kb * NYC) + j;
            const size_t row = (kb * NY) + constraints.y_outputs[j];
            for (size_t c = 0; c < NZ; ++c) {
                art.A_con(hi, c) = pred.Gamma(row, c);
                art.A_con(lo, c) = -pred.Gamma(row, c);
            }
            for (size_t c = 0; c < NXA; ++c) {
                art.B_map(hi, c) = -pred.Phi(row, c);
                art.B_map(lo, c) = pred.Phi(row, c);
            }
            art.A_con(hi, NZ) = -constraints.y_upper_ecr(j);
            art.A_con(lo, NZ) = -constraints.y_lower_ecr(j);
        }
    }
    if constexpr (NXC > 0) {
        //! The bounded states over the horizon, scaled: the state rows of the same prediction, x = S_x·Ã^k·x̃ + ….
        constexpr size_t    r_x_hi = (4 * NZ) + (2 * NYCS);
        constexpr size_t    r_x_lo = r_x_hi + NXS;
        Matrix<NXC, NXA, T> Sx{};
        for (size_t i = 0; i < NXC; ++i) {
            Sx(i, constraints.x_states[i]) = T{1} / art.scale_x(i);
        }
        const auto xpred = mpc_predictions<NP, NC, NXA, NX, NU, NXC, T>(Aa, Ba, Sx);
        art.A_con.template block<NXS, NZ>(r_x_hi, 0) = xpred.Gamma;
        art.A_con.template block<NXS, NZ>(r_x_lo, 0) = -xpred.Gamma;
        art.B_map.template block<NXS, NXA>(r_x_hi, 0) = -xpred.Phi;
        art.B_map.template block<NXS, NXA>(r_x_lo, 0) = xpred.Phi;
        for (size_t kb = 0; kb < NP; ++kb) {
            for (size_t j = 0; j < NXC; ++j) {
                art.A_con(r_x_hi + (kb * NXC) + j, NZ) = -constraints.x_upper_ecr(j);
                art.A_con(r_x_lo + (kb * NXC) + j, NZ) = -constraints.x_lower_ecr(j);
            }
        }
    }
    detail::mpc_box_bounds(art, MPCLimits<NU, NYC, T, NXC>{constraints.u_lim, constraints.du_lim, constraints.y_lim, constraints.x_lim});

    if constexpr (NT_max > 0) {
        if (constraints.terminal == MPCTerminal::invariant_set) {
            if (!damp::finite_positive(constraints.terminal_target_margin) || !(constraints.terminal_target_margin < T{1})) {
                return ::damp::detail::set_status(art, DesignStatus::invalid_parameter);
            }
            //! The set of w = [z; q] under w⁺ = blkdiag(A_z − B_z·K, I)·w: valid for every target
            //! q = [y_t; u_t] in the tightened boxes, the target entering each row's offset.
            constexpr size_t  NZD = NX + NU;
            constexpr size_t  NA = NZD + NY + NU;
            constexpr size_t  NG = (4 * NYC) + (6 * NU);
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
 * @tparam Target    SteadyStateTargetResult on the estimator's model (the law's without one), on
 *                   the outputs it holds
 * @tparam NRU       Delay states the runtime holds per input: a measured-state MPC on a plant whose
 *                   input arrives late (InputDelay), NRU·NU states after the plant's (0: none; a
 *                   filter carries its own)
 * @tparam NH        Outputs the target holds, as a TrackedOutputs or the disturbance spec lists
 *                   them (0: every output, r the QP's reference on each)
 */
template<typename Law, typename Estimator, typename Target, size_t NRU = 0, size_t NH = 0>
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
        return MPCResult<LawU, EstimatorU, TargetU, NRU, NH>{law, kalman.template as<U>(), target.template as<U>(), st == DesignStatus::ok, st};
    }
};

namespace detail {

/** @brief Not a measured disturbance. */
template<typename>
struct is_measured_disturbance : damp::false_type {};

/** @brief MeasuredDisturbance. */
template<size_t NX, size_t NM, typename T>
struct is_measured_disturbance<MeasuredDisturbance<NX, NM, T>> : damp::true_type {};

/** @brief The options design::mpc takes after the weights. */
using mpc_options = ::damp::detail::option_set<
    is_input_limits, is_kalman_noise, is_disturbance, is_measured_disturbance, is_tracked_option, is_output_limits, is_state_limits, is_mpc_terminal,
    is_qp_iterations, is_input_delay, is_output_delay>;

/**
 * @brief The compile-time check of a design::mpc call's trailing options: a kind it does not take
 *        or one given twice is refused, with a message naming the spelling for the limits, the
 *        held outputs, the terminal ingredient and the QP budget.
 */
template<typename... Opts>
constexpr void check_mpc_options() {
    static_assert(::damp::detail::option_count<is_input_limits, Opts...> <= 1, "One InputLimits bounds every input: InputLimits{.u_lim = {...}, .du_lim = {...}}.");
    static_assert(::damp::detail::option_count<is_output_limits, Opts...> <= 1, "One OutputLimits lists every bounded output: OutputLimits{.outputs = {i, ...}, .lim = {...}}.");
    static_assert(::damp::detail::option_count<is_state_limits, Opts...> <= 1, "One StateLimits lists every bounded state: StateLimits{.states = {i, ...}, .lim = {...}}.");
    static_assert(::damp::detail::option_count<is_tracked_option, Opts...> <= 1, "One TrackedOutputs lists every output the target holds at r: TrackedOutputs{{i, ...}}.");
    static_assert(
        ::damp::detail::option_count<is_tracked_option, Opts...> == 0 || !(is_disturbance_model<damp::remove_cvref_t<Opts>> || ...),
        "A disturbance spec holds the outputs it lists at r: give them as its .outputs, InputDisturbance{{i, ...}} or "
        "DisturbanceModel{.Bd = ..., .outputs = {i, ...}}, not as TrackedOutputs."
    );
    static_assert(
        ::damp::detail::option_count<is_mpc_terminal, Opts...> <= 1,
        "One terminal ingredient: DareCost{} for the cost, or InvariantSet{.target_margin = ...} for the cost and the set, not both."
    );
    static_assert(::damp::detail::option_count<is_qp_iterations, Opts...> <= 1, "One QP budget: QPIterations{n}.");
    static_assert(!(damp::is_integral_v<damp::remove_cvref_t<Opts>> || ...), "The QP's iteration budget is QPIterations{n}; the terminal set's room is mpc<NP, NC, NT_max>.");
    mpc_options::check<Opts...>();
}

/** @brief Whether a design::mpc call has an estimator: a KalmanNoise among its options. */
template<typename... Opts>
inline constexpr bool mpc_estimated = ::damp::detail::option_count<is_kalman_noise, Opts...> > 0;

/** @brief The first status that is not ok among @p statuses, else ok. */
template<typename... S>
[[nodiscard]] constexpr DesignStatus first_failure(S... statuses) {
    DesignStatus out = DesignStatus::ok;
    ((out = (out == DesignStatus::ok) ? statuses : out), ...);
    return out;
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
 * @brief The weight of the target's least squares (more held outputs than inputs): the rows and
 *        columns of the cost's output weight @p W on the outputs @p spec holds, every output when it
 *        lists none. An index out of range or listed twice gives the identity: the target refuses
 *        the list.
 */
template<size_t NU, size_t NW, typename Spec, size_t NY, typename T>
[[nodiscard]] constexpr target_weight<Spec, NU, NY, NW, T> held_weight(const Matrix<NY, NY, T>& W, const Spec& spec) {
    constexpr size_t NZ = target_outputs<Spec, NU, NY, NW>;
    const auto       rows = selected_outputs<NY>(listed_outputs(spec));
    if (!rows) {
        return Matrix<NZ, NZ, T>::identity();
    }
    Matrix<NZ, NZ, T> Q{};
    for (size_t i = 0; i < NZ; ++i) {
        for (size_t j = 0; j < NZ; ++j) {
            Q(i, j) = W((*rows)[i], (*rows)[j]);
        }
    }
    return Q;
}

/**
 * @brief The spec a measured-state MPC's target is designed on: the outputs @p held lists, under the
 *        measured disturbance as a disturbance of the steady state that no reading feels
 *        (DisturbanceModel{.Bd, .outputs}, C_d = 0), or the TrackedOutputs alone without one.
 */
template<size_t NXD, size_t NM, typename T, size_t NS>
[[nodiscard]] constexpr auto measured_target_spec(const MeasuredDisturbance<NXD, NM, T>& measured, const TrackedOutputs<NS>& held) {
    if constexpr (NM == 0) {
        static_cast<void>(measured);
        return held;
    } else {
        DisturbanceModel<NXD, 0, NM, NS, T> spec{.Bd = measured.Bd};
        if constexpr (NS > 0) {
            for (size_t i = 0; i < NS; ++i) {
                spec.outputs[i] = held.outputs[i];
            }
        } else {
            static_cast<void>(held);
        }
        return spec;
    }
}

/** @brief The spec a filtered MPC's target is designed on: the disturbance model, which holds the outputs it lists, else the outputs @p held lists (TrackedOutputs, no disturbance). */
template<typename Spec, size_t NS>
[[nodiscard]] constexpr auto target_spec_of(const Spec& spec, const TrackedOutputs<NS>& held) {
    if constexpr (is_disturbance_model<Spec>) {
        static_cast<void>(held);
        return spec;
    } else {
        static_cast<void>(spec);
        return held;
    }
}

/**
 * @brief The continuous stage cost of a design from @p sys at @p Ts: the plant, its outputs and the
 *        disturbance's continuous channel @p Bd_c; off on a discrete plant (@p Discrete).
 */
template<bool Discrete, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t ND>
[[nodiscard]] constexpr MpcContinuousCost<NX, NU, NY, ND, T> mpc_continuous_cost(
    const StateSpace<NX, NU, NY, T, NW, NV>&                 sys,
    const Matrix<NX, ND, T>&                                 Bd_c,
    T                                                        Ts,
    const damp::array<::damp::detail::RegisterDelay<T>, NU>& inputs = {}
) {
    MpcContinuousCost<NX, NU, NY, ND, T> out{};
    if constexpr (!Discrete) {
        out = MpcContinuousCost<NX, NU, NY, ND, T>{.A = sys.A, .B = sys.B, .Bd = Bd_c, .C = sys.C, .D = sys.D, .Ts = Ts, .on = true, .inputs = inputs};
    } else {
        static_cast<void>(sys);
        static_cast<void>(Bd_c);
        static_cast<void>(Ts);
        static_cast<void>(inputs);
    }
    return out;
}

/** @brief Each input's delay as discretize split it, d·Ts + τ: none without an InputDelay. */
template<size_t NU, typename T, typename Delay>
[[nodiscard]] constexpr damp::array<::damp::detail::RegisterDelay<T>, NU> input_delays_of(const Delay& delay, T Ts) {
    if constexpr (is_input_delay<Delay>::value) {
        const auto split = ::damp::detail::channel_register_delays<NU>(delay.seconds, Ts, delay_registers<Delay>);
        return split ? *split : damp::array<::damp::detail::RegisterDelay<T>, NU>{}; // sampled_model refused a delay discretize does not take
    } else {
        static_cast<void>(delay);
        static_cast<void>(Ts);
        return {};
    }
}

/**
 * @brief The MPC on the measured state: the measurement is x, the QP on the plant's model and its
 *        outputs C·x with the measured-disturbance channel, and the target of
 *        design::steady_state_target on the outputs @p held lists, under the measured d.
 */
template<
    bool Discrete, size_t NP, size_t NC, size_t NT_max, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NYC, size_t NXC, size_t NXD, size_t NM,
    size_t NS, typename Delay = NoDelay>
[[nodiscard]] constexpr auto mpc_measured(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        Ts,
    const MPCWeights<NU, NY, T>&             weights,
    const MPCConstraintSet<NU, NYC, T, NXC>& constraints,
    const MeasuredDisturbance<NXD, NM, T>&   measured,
    const TrackedOutputs<NS>&                held,
    const Delay&                             delay = Delay{}
) {
    static_assert(NM == 0 || NXD == NX, "The measured disturbance enters the plant's states: MeasuredDisturbance<NX, NM>.");
    constexpr size_t NXL = NX + (delay_registers<Delay> * NU); // the plant, then the delay states
    const auto       spec = measured_target_spec(measured, held);
    using Law = MPCArtifacts<NXL, NU, NY, NP, NC, NM, T, NT_max, NXC, NYC>;
    using Target = SteadyStateTargetResult<NXL, NU, NM, target_outputs<damp::remove_cvref_t<decltype(spec)>, NU, NY, NW>, T>;
    using Result = MPCResult<Law, NoEstimator, Target, delay_registers<Delay>, NS>; // the runtime holds the delay states
    const auto m = sampled_model(sys, Ts, Discrete, delay);
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
    }
    Result     out{};
    const auto cost = mpc_continuous_cost<Discrete>(sys, measured.Bd, Ts, input_delays_of<NU>(delay, Ts));
    out.mpc = mpc_impl<NP, NC, NM, NT_max>(m.law, Bd, Matrix<NY, NM, T>{}, weights, constraints, cost);
    out.target = target_design<Discrete>(sys, Ts, spec, held_weight<NU, NW>(physical_output_weight(weights.Qy, out.mpc), spec), delay);
    return ::damp::detail::set_status(out, first_failure(out.mpc.status, out.target.status));
}

/**
 * @brief The MPC through an estimator: the filter of design::kalman with the same options, the
 *        QP on the law's model and every output with the disturbance channel, and the target of
 *        design::steady_state_target: under the disturbance spec, holding the outputs it lists,
 *        else on the outputs @p held lists.
 */
template<
    bool Discrete, size_t NP, size_t NC, size_t NT_max, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NWN, size_t NVN, size_t NDN,
    typename Spec, size_t NS, size_t NYC, size_t NXC, typename Delay = NoDelay>
[[nodiscard]] constexpr auto mpc_estimated_design(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        Ts,
    const MPCWeights<NU, NY, T>&             weights,
    const MPCConstraintSet<NU, NYC, T, NXC>& constraints,
    const KalmanNoise<NWN, NVN, T, NDN>&     noise,
    const Spec&                              spec,
    const TrackedOutputs<NS>&                held,
    const Delay&                             delay = Delay{}
) {
    constexpr size_t ND = disturbance_traits<Spec, NU, NY, NW>::states;
    static_assert(NWN == NW && NVN == NV, "KalmanNoise is sized to the plant's noise inputs: KalmanNoise<NW, NV> for the columns of G and H.");
    static_assert(NDN == ND, "KalmanNoise carries the random walk of the spec's ND disturbances: KalmanNoise<NW, NV, T, ND>.");
    constexpr size_t NXL = NX + (delay_registers<Delay> * NU); // the QP's: the plant, then the input's delay states
    constexpr size_t NXE = NX + delay_states<Delay, NU, NY>;   // the filter's: then the late readings' delay states
    const auto&      in = input_delay_of(delay);
    const auto       kf = kalman_design<Discrete>(sys, Ts, noise, spec, delay);
    const auto       aim = target_spec_of(spec, held);
    using Aim = damp::remove_cvref_t<decltype(aim)>;
    using Target = SteadyStateTargetResult<NXE, NU, ND, target_outputs<Aim, NU, NY, NW>, T>;
    using Law = MPCArtifacts<NXL, NU, NY, NP, NC, ND, T, NT_max, NXC, NYC>;
    using Result = MPCResult<Law, damp::remove_cvref_t<decltype(kf)>, Target, 0, disturbance_traits<Aim, NU, NY, NW>::listed>;
    Result     out{.kalman = kf};
    const auto m = sampled_model(sys, Ts, Discrete, in);
    if (m.status != DesignStatus::ok) {
        return ::damp::detail::set_status(out, first_failure(kf.status, m.status)); // a delay refused
    }
    //! The QP predicts, weighs and bounds every output of the law's model, the estimate riding B_d and C_d
    const auto         channels = disturbance_channels(spec, sys, m.law, Ts);
    const auto         cost = mpc_continuous_cost<Discrete>(sys, channels.Bd_c, Ts, input_delays_of<NU>(in, Ts));
    Matrix<NXL, ND, T> Bd{};
    if constexpr (ND > 0) {
        Bd = channels.Bd.template block<NXL, ND>(0, 0);
    }
    out.mpc = mpc_impl<NP, NC, ND, NT_max>(m.law, Bd, channels.Cd, weights, constraints, cost);
    //! On the filter's model: under d̂ with a disturbance spec; a late reading's delay states hold the steady output
    out.target = target_design<Discrete>(sys, Ts, aim, held_weight<NU, NW>(physical_output_weight(weights.Qy, out.mpc), aim), delay);
    return ::damp::detail::set_status(out, first_failure(kf.status, out.target.status, out.mpc.status));
}

/** @brief design::mpc in the plant's time domain (@p Discrete) or sampled at @p Ts, once the options pass check_mpc_options: the options sorted by kind. */
template<bool Discrete, size_t NP, size_t NC, size_t NT_max, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename... Opts>
[[nodiscard]] constexpr auto mpc_design(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T Ts, const MPCWeights<NU, NY, T>& weights, const Opts&... opts) {
    const auto constraints = constraint_set<NU, T>(opts...);
    const auto held = ::damp::detail::option<is_tracked_option>(TrackedOutputs<>{}, opts...);
    static_assert(
        !Discrete || ::damp::detail::option_count<is_input_delay, Opts...> + ::damp::detail::option_count<is_output_delay, Opts...> == 0,
        "A delay is sampled with the plant: give the continuous plant and Ts, mpc<NP, NC>(plant, Ts, weights, ..., InputDelay<N>{theta}, OutputDelay<M>{theta})."
    );
    if constexpr (mpc_estimated<Opts...>) {
        static_assert(
            ::damp::detail::option_count<is_measured_disturbance, Opts...> == 0,
            "A measured disturbance with an estimator is not supported: the filter does not take it as an input."
        );
        return mpc_estimated_design<Discrete, NP, NC, NT_max>(
            noise_model(sys), Ts, weights, constraints,
            ::damp::detail::option<is_kalman_noise>(KalmanNoise<noise_inputs<NX, NW>, noise_inputs<NY, NV>, T>{}, opts...),
            ::damp::detail::option<is_disturbance>(NoDisturbance{}, opts...), held, delays_of(opts...)
        );
    } else {
        static_assert(::damp::detail::option_count<is_disturbance, Opts...> == 0, "A disturbance spec needs KalmanNoise: the disturbance is estimated.");
        static_assert(::damp::detail::option_count<is_output_delay, Opts...> == 0, "An OutputDelay is a late reading: it needs a KalmanNoise, a filter to carry it.");
        return mpc_measured<Discrete, NP, NC, NT_max>(
            sys, Ts, weights, constraints, ::damp::detail::option<is_measured_disturbance>(MeasuredDisturbance<NX, 0, T>{}, opts...), held,
            ::damp::detail::option<is_input_delay>(NoDelay{}, opts...)
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
 * @f$ f = [F_x \tilde x_0 - F_r y'_r - F_{ut} u'_t;\; 0] @f$, and stacks
 * the hard u/Δu and soft y box constraints as
 * @f$ A_{con} \zeta \le b_{bound} + B_{map}\tilde x_0 @f$.
 * Beyond NC−1 the input is held (Δu = 0), which the augmented model carries automatically.
 *
 * | Design | Runtime |
 * | --- | --- |
 * | mpc<NP, NC>(plant, [Ts,] weights) | MPC: the state measured |
 * | mpc<NP, NC>(plant, [Ts,] weights, MeasuredDisturbance{Bd}) | MPC with a known disturbance input, given each tick |
 * | mpc<NP, NC>(plant, [Ts,] weights, noise, disturbance) | OFMPC: a Kalman filter on [x; d], offset-free (below) |
 * | … plus InputLimits{.u_lim, .du_lim} | the same, with hard limits on u and Δu |
 * | … plus OutputLimits{.outputs, .lim} | the same, with limits on the outputs it lists |
 * | … plus StateLimits{.states, .lim} | the same, with limits on states it does not track |
 * | … plus TrackedOutputs{{i, …}}, without a disturbance spec | the same, r on the outputs it lists, held through the target |
 * | … plus DareCost{} or mpc<NP, NC, NT_max> and InvariantSet{.target_margin} | the same, with a terminal cost, or that cost and a terminal set |
 * | … plus QPIterations{n} | the same, its QP capped at n iterations a tick |
 *
 * The plant and the weights are positional, so braces size the weights from the plant. The rest
 * trails as options, at most one of each kind, in any order: the hard limits on u and Δu (none:
 * unbounded), the limits on the outputs and on the states, the outputs r holds (none: every
 * output), the terminal ingredient (none: the plain cost) and the QP budget (none: the default,
 * 10·(NZ + 1 + NI + 1)).
 *
 * Outputs. The QP predicts, weighs and bounds every output of the plant, the rows of C: Qy is
 * NY × NY and an OutputLimits index is a row of C. On the measured state design::mpc's runtime
 * takes the plant state x as its measurement (update(x)), and C names those outputs, y = C·x, as
 * the LQI on a measured state integrates r − C·x; a state the cost weighs is a row of C. With a
 * KalmanNoise, C and D are the measured outputs, which the filter reconstructs the state from.
 *
 * The target. The runtime poses the QP about the steady state (x_t, u_t) of the reference, the
 * target of design::steady_state_target with the design's options, as the LQ laws do. It holds at
 * r the outputs a TrackedOutputs lists, or with a disturbance spec those its outputs lists, under
 * the estimated disturbance; every output when neither lists any, and under a measured
 * disturbance on the measured state. r is on those outputs, in their order: TrackedOutputs{{0}}
 * on a plant with four outputs takes a scalar r. With NZ = NU held outputs the steady state is
 * unique or the design is refused (DesignStatus::singular), with fewer the minimum-norm steady
 * state that holds every r (refused when none does), and with more the equilibrium closest to r
 * in the Qy-weighted norm (an output without cost does not pull it). The QP's reference on each
 * output at each step is r itself when the target holds every output, else that steady state's
 * output for r at that step, so an output r does not list settles where the target puts it,
 * the lean at which gravity balances a load. Ru weighs the input's deviation from the target's
 * (set_input_target replaces it), so it leaves no offset. An equilibrium leaves the QP the
 * absolute one; the terminal cost and set are on the deviation from it. A TrackedOutputs index
 * past the plant's outputs or listed twice is refused (DesignStatus::bad_weights), as a disturbance
 * spec's is, and a TrackedOutputs with a disturbance spec does not compile: the spec's outputs are
 * the ones it holds.
 *
 * Input delay. On the sampled form, mpc<NP, NC>(plant, Ts, weights, …, InputDelay<N>{θ})
 * designs on the plant whose input arrives θ late, the model discretize(plant, Ts, InputDelay<N>{θ}):
 * the plant, then the delay states v(k−1) … v(k−N), the commands still on their way. The weights
 * and the noise are the plant's, as without the delay; the output cost is integrated exactly over
 * each period, the plant holding the older command until each input's switch. A measured-state MPC
 * takes the plant's state, update(x), and holds the delay states itself from the input
 * back_calculate recorded (else its command); a filter carries them, as design::kalman with the
 * same delay does. A delay on a discrete plant does not compile: it is sampled with the plant. A
 * discrete model with the delay built in is designed on like any plant, its delay states then part
 * of the measured state and its noise inputs named. With a KalmanNoise, OutputDelay<M>{θ} too: the
 * reading arrives θ late, and the filter and the target are design::kalman's and
 * steady_state_target's with the same delays, the late readings' delay states the filter's alone;
 * the QP's model and its weights are as without it. The first outputs of the horizon are ones no
 * decision can change, so a longer delay wants a longer NP.
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
 * design::mpc with a KalmanNoise and a disturbance spec is the offset-free MPC of Muske and
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
 * @f$ B_d @f$ and @f$ C_d @f$, so the predictions and the output limits carry it. Each
 * tick the steady-state target holds the outputs the spec lists, @f$ z = C_z x + C_{d,z} \hat d @f$,
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
 * when it lists none), as in the OFLQG; r is on z, while the weights and an OutputLimits are on
 * every output, so DisturbanceModel{.Bd = …, .outputs = {0}} on a plant with four readings
 * weighs all four and holds the first at a scalar r, the others at the target's. Offset-free
 * tracking of every listed output needs NZ ≤ NU and a steady state that reaches any r;
 * with NZ = NU the target is unique or the design is refused (DesignStatus::singular), and with
 * more listed outputs than inputs it is the least-squares compromise. The augment must be
 * detectable, rank [A − I, B_d; C, C_d] = NX + ND (so ND ≤ NY), or the design is refused with
 * DesignStatus::not_detectable: an output disturbance on a plant with a pole at z = 1 is not,
 * which is why MATLAB® falls back to input disturbances on integrating plants.
 *
 * The estimator holds deviations from the target (FixedGainObserver with the fixed gain of
 * design::kalman, anchored by SteadyStateTarget each tick), so the QP's gradient stays a sum of
 * small terms however far the plant sits from zero (a temperature at 200 °C, a position in
 * metres at micrometre precision).
 *
 * Fails (success = false) if the plant is not discrete (DesignStatus::invalid_sample_time; use
 * the overload with the sample time), an output the QP weighs or bounds has feedthrough, the slack
 * penalty is non-positive, the Hessian is not positive definite, or the filter or the target is
 * refused.
 *
 * @note Compare with MATLAB®'s mpc(plant, Ts, p, m) — p = NP, m = NC — with its built-in
 *       estimator (setindist / setoutdist, getEstimator) and mpcmove(mpcobj, xc, ym, r).
 * @see MPCWeights, InputLimits, OutputLimits, StateLimits, TrackedOutputs, DareCost, InvariantSet, QPIterations,
 *      damp::MPC, design::kalman, design::steady_state_target
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
 * @tparam NT_max Terminal-set row capacity (default 0; required > 0 for an InvariantSet)
 * @param sys      Discrete-time plant (Ts > 0); C its outputs, D = 0 on those the QP weighs or bounds
 * @param weights  Cost weights on every output of the plant, per prediction step (default: Qy = I, RΔu = I)
 * @param opts     In any order: InputLimits{.u_lim, .du_lim}, the hard limits on u and Δu;
 *                 MeasuredDisturbance<NX, NM>{.Bd}, a disturbance input measured each tick
 *                 (control(r, d)); or a KalmanNoise, the MPC then on a Kalman filter's estimate from the
 *                 measured outputs C and D, offset-free with a disturbance spec (runs as OFMPC: the
 *                 filter on [x; d], the QP posed about the steady-state target under d̂, which holds
 *                 the rows the spec lists, all of them without a list); TrackedOutputs{{i, …}}
 *                 without a disturbance spec, the outputs r holds through the target;
 *                 OutputLimits{.outputs, .lim}, limits on the outputs it lists;
 *                 StateLimits{.states, .lim}, limits on states the QP does not track; DareCost{} or
 *                 InvariantSet{.target_margin}, a terminal ingredient; and QPIterations{n}, the QP's budget
 * @return MPCResult: deploy with `MPC law{res}`
 */
template<size_t NP, size_t NC, size_t NT_max = 0, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename... Opts>
[[nodiscard]] constexpr auto mpc(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const damp::type_identity_t<MPCWeights<NU, NY, T>>& weights, const Opts&... opts) {
    detail::check_mpc_options<Opts...>();
    return detail::mpc_design<true, NP, NC, NT_max>(sys, T{0}, weights, opts...);
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
 * terminal cost (DareCost, InvariantSet) is the cost-to-go of that period cost. The same weights at a
 * shorter Ts, with the horizons in steps rescaled to the same seconds, give nearly the same
 * loop. An OutputLimits holds at the samples. With KalmanNoise the estimator is
 * continuous in, discrete out (design::kalman(plant, Ts, noise, …)): the process and
 * disturbance noise are densities, discretized exactly over the period (Van Loan), and the
 * measurement noise is the variance of one reading.
 *
 * Refused with DesignStatus::invalid_sample_time when @p plant is already discrete or @p Ts is
 * not finite and positive, and as the discrete form is otherwise.
 *
 * @tparam NP     Prediction horizon (steps)
 * @tparam NC     Control horizon (moves, NC ≤ NP)
 * @tparam NT_max Terminal-set row capacity (default 0; required > 0 for an InvariantSet)
 * @see "Computing Integrals Involving the Matrix Exponential" (Van Loan, 1978), IEEE TAC 23(3)
 * @see lqr(plant, Ts, …) for the same discretization of an LQ cost
 *
 * @param plant    Continuous-time plant (Ts = 0)
 * @param Ts       Control period [s]
 * @param weights  Cost weights on every output of the plant, per second
 * @param opts     As the discrete form's, a MeasuredDisturbance's B_d continuous; InputDelay<N>{θ}, the
 *                 plant's input θ late, and OutputDelay<M>{θ} with a KalmanNoise, the reading θ late
 *                 (Input delay, above)
 * @return MPCResult: deploy with `MPC law{res}`; success = false if refused or the plant does not discretize
 */
template<size_t NP, size_t NC, size_t NT_max = 0, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename... Opts>
[[nodiscard]] constexpr auto mpc(
    const StateSpace<NX, NU, NY, T, NW, NV>&            plant,
    damp::type_identity_t<T>                            Ts,
    const damp::type_identity_t<MPCWeights<NU, NY, T>>& weights,
    const Opts&... opts
) {
    detail::check_mpc_options<Opts...>();
    return detail::mpc_design<false, NP, NC, NT_max>(plant, Ts, weights, opts...);
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
 * constexpr auto art = design::mpc<h.prediction_horizon, h.control_horizon>(sys, {});
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
    if (!damp::finite_positive(settling_time) || !damp::finite_positive(Ts)) {
        return s;
    }
    const T steps = settling_time / Ts;
    // NP is a template argument: anything near size_t range is not a usable
    // suggestion, and the cast below would be UB past it.
    if (!(steps < static_cast<T>(damp::numeric_limits<size_t>::max() / 2))) {
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
 * steps before every input's move reaches an output, an input delay's delay states included, which
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
    if (!damp::finite_positive(sys.Ts)) {
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
        dist[i] = damp::abs(eig.values[i] - T{1});
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
                mean += eig.values[order[k]];
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
        const T mag = damp::abs(eig.values[i]);
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
    static constexpr bool built_from = damp::is_same_v<Design, NoEstimator>;
};

/** @brief A FixedGainObserver: its sizes, and the design::kalman results it runs (any noise inputs). */
template<size_t NXK, size_t NU, size_t NY, typename T, size_t ND, typename Store>
struct mpc_estimator_traits<FixedGainObserver<NXK, NU, NY, T, ND, Store>> {
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
 * dual, not primal, feasible).
 *
 * A solve that does not converge. last_status() is the only report of it: the tick returns a
 * command either way. A MaxIterations or Degenerate iterate is applied when the outputs are soft
 * (the default). It is the active set's minimizer over the limits it had taken in when it
 * stopped: the plan of a QP with the others left out, whose first move the boxes then clamp.
 * Holding the input instead is the worse fallback on a loop that has to move. On a double
 * integrator brought back from outside its soft limits (NP = 30, NC = 10) with a budget of 12
 * iterations where a tick needs up to 17, the applied iterates settle the loop at 1.4 times the
 * converged Σy², and holding leaves the input where the last solved tick put it: the output
 * runs away. If any output or state ECR is 0 or a terminal set is present, the iterate is applied only
 * when it is primal feasible; otherwise the tick holds the last applied input, as Infeasible
 * does. Holding the input does not make the predicted outputs or a terminal set feasible, and
 * a hard limit the plan rides at full input can read Infeasible in single precision, the
 * rounded state lying past what the input can still reach (2 of 300 ticks braking onto a hard
 * ceiling).
 *
 * The soft-constraint
 * slack actually used is reported by last_slack(). @p r is the reference for the outputs the
 * target holds (held, or a Preview), not a state x_ref: the cost is ‖y − y_r‖²_Qy on every output,
 * y_r the target's output for r when the target holds fewer than all of them, plus move and
 * input-target terms, with no ‖x − x_ref‖ term (see the file comment). The input target is the
 * target's input for r at each move, unless set_input_target set one.
 *
 * One tick is @ref step(): @ref predict holds the input the plant received over [k−1, k)
 * (the base of the next Δu) and runs the estimator's time update, @ref update takes the
 * measurement m_k, and @ref control solves for u_k. The three are public, for the low-latency
 * order update → control → write the output → back_calculate → predict (the trailing predict is
 * the next sample's time update). When a downstream stage can
 * clamp or override the command, @ref back_calculate records the input that went out, so the
 * next move and the estimator start from it; without it predict holds the last command.
 *
 * - No estimator (Estimator = NoEstimator): the measurement is the plant state x; the outputs
 *   are C·x.
 * - FixedGainObserver (the fixed gain of design::kalman): the measurement is the plant's outputs; the
 *   observer reconstructs the state.
 *   With a disturbance model the estimator holds deviations from the steady-state target,
 *   re-anchored each tick on the reference now (SteadyStateTarget), and the law is offset-free.
 *   update clamps d̂ into disturbance_lim (unbounded by default), so a load past it (a wall the
 *   plant is held against) does not wind d̂ and the target up.
 *   Without one the QP is posed about the target of the reference without a disturbance.
 *
 * @tparam NX        Plant states
 * @tparam NY        Outputs: the QP's weights and predictions, every row of the plant's C
 * @tparam ND        Measured-disturbance inputs (0: none)
 * @tparam Solver    QP solver policy (default: the warm-started Goldfarb–Idnani active set)
 * @tparam NT_max    Terminal-set row capacity (default 0)
 * @tparam Estimator FixedGainObserver, or NoEstimator for the measured state; it holds the
 *                   disturbance when there is one
 * @tparam NXC       States with bounds (the design's StateLimits; 0: none)
 * @tparam NYC       Outputs with bounds (the design's OutputLimits; 0: none)
 * @tparam NRU       Delay states the law holds per input (a measured state with an InputDelay; 0: none)
 * @tparam NH        Outputs the target holds, listed by a TrackedOutputs or the disturbance spec
 *                   (0: every output, r the QP's reference on each)
 *
 * @note Compare with MATLAB®'s mpcmove(mpcobj, xc, ym, r, v).
 * @see design::mpc
 */
template<
    size_t NX, size_t NU, size_t NY, size_t NP, size_t NC, size_t ND, typename T, typename Solver, size_t NT_max, typename Estimator, size_t NXC = 0, size_t NYC = 0,
    size_t NRU = 0, size_t NH = 0>
    requires damp::is_floating_point_v<T>
class MPCCore {
    using Traits = detail::mpc_estimator_traits<Estimator>;

public:
    static constexpr bool   estimated = Traits::estimated;                       ///< Whether a Kalman filter reconstructs the state
    static constexpr size_t NDE = Traits::disturbances;                          ///< Estimated disturbance states (0: not offset-free)
    static constexpr size_t NXM = NX;                                            ///< States of the prediction model
    static constexpr size_t NYM = estimated ? Traits::outputs : NX - (NRU * NU); ///< Measured outputs (the plant's state without an estimator; the law holds its delay states)
    static constexpr size_t NXE = estimated ? (Traits::states - NDE) : NXM;      ///< States of the estimator's model before the disturbance
    static_assert(!estimated || ND == 0, "A measured disturbance with an estimator is not supported.");
    static_assert(damp::is_same_v<Estimator, NoEstimator> || estimated, "The estimator is a FixedGainObserver, or NoEstimator for the measured state.");
    static_assert(NRU == 0 || !estimated, "A filter carries the delay states itself.");

    static constexpr size_t NREF = (NH == 0) ? NY : NH; ///< References r: one per output the target holds

    using Artifacts = design::MPCArtifacts<NXM, NU, NY, NP, NC, ND + NDE, T, NT_max, NXC, NYC>;
    using TargetResult = design::SteadyStateTargetResult<NXE, NU, ND + NDE, NREF, T>;
    using DisturbanceLim = detail::maybe<(NDE > 0), Bounds<NDE, T>, detail::OffSlot::disturbance_lim>;

    /** @brief Limits on the disturbance estimate d̂, in d's units: update clamps it, as OFLQG's (NDE > 0; default unbounded). */
    [[no_unique_address]] DisturbanceLim disturbance_lim{};
    /** @brief The filter and its anchor, as the LQ laws' estimator: restart it with estimator.reset(x0, d0) (no storage on the measured state). */
    [[no_unique_address]] Estimator estimator{};
    /** @brief A design it is built from: the QP and the target, and @p Kalman, the estimator's design (NoEstimator without one). */
    template<typename Kalman>
    using ResultWith = design::MPCResult<Artifacts, Kalman, TargetResult, NRU, NH>;
    static constexpr size_t NXA = Artifacts::NXA;
    static constexpr size_t NZ = Artifacts::NZ;
    static constexpr size_t NZ1 = Artifacts::NZ1;
    static constexpr size_t NI1 = Artifacts::NI1;
    static constexpr size_t NYS = Artifacts::NYS;
    /** @brief The reference over the prediction horizon, one value per step for the outputs the target holds (MATLAB's `ref` with preview). */
    using Preview = damp::array<ColVec<NREF, T>, NP>;

    /** @brief The input target over the control horizon, one value per move (MATLAB's `mv.target` with preview). */
    using InputPreview = damp::array<ColVec<NU, T>, NC>;

    constexpr MPCCore() = default;

    /** @brief From the design. A failed design is a compile error when constant-evaluated (`static constinit MPC law{art};`). */
    template<typename Kalman>
        requires(Traits::template built_from<Kalman>)
    constexpr explicit MPCCore(const ResultWith<Kalman>& result)
        : estimator(make_estimator(result.kalman)), art(artifacts_of(result)), target_(target_of(result.target)), y_r_(held_map(art, target_)) {
        if (!result.success && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(result.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    /** @brief Construct with a specific solver instance (for stateful solver policies) */
    template<typename Kalman>
        requires(Traits::template built_from<Kalman>)
    constexpr MPCCore(const ResultWith<Kalman>& result, const Solver& solver)
        : estimator(make_estimator(result.kalman)), art(artifacts_of(result)), qp_solver_(solver), target_(target_of(result.target)), y_r_(held_map(art, target_)) {
        if (!result.success && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(result.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    /**
     * @brief The reference over the prediction horizon: one value held, or a preview
     *
     * On the outputs the target holds, in their order: every output, or those a TrackedOutputs or
     * the disturbance spec lists. Built implicitly from a ColVec<NREF> (held over the NP steps),
     * from one value per output in braces (`control({0.3, 0.0})`), from the scalar of a single
     * one, or from a Preview (the
     * reference at each of steps 1 … NP), so every tick takes any of them. From a continuous plant
     * r[k] is the reference over period k, [k·Ts, (k+1)·Ts) (the cost integrates each period),
     * so r[0] covers the period starting now. With a preview the plan
     * sees a setpoint change before it happens and starts moving ahead of it; a trajectory
     * generator hands in its next NP points each tick (repeat the last known value to fill).
     * Output references, not state references.
     */
    struct Reference {
        Preview r{}; ///< The output reference at steps 1 … NP (periods 0 … NP − 1 from a continuous plant)

        constexpr Reference(const ColVec<NREF, T>& held_value) : r(held(held_value)) {} // NOLINT(google-explicit-constructor)
        constexpr Reference(const Preview& preview) : r(preview) {}                     // NOLINT(google-explicit-constructor)

        /** @brief A reference on one output from its scalar, held over the NP steps (as ColVec<1>'s). */
        constexpr Reference(T held_value) // NOLINT(google-explicit-constructor)
            requires(NREF == 1)
            : r(held(ColVec<NREF, T>{held_value})) {}

        /** @brief One value per output, braced at the call (`control({0.3, 0.0})`), held over the NP steps. */
        template<typename... V>
            requires(NREF > 1 && sizeof...(V) == NREF && (damp::convertible_to<V, T> && ...))
        constexpr Reference(V... held_values) // NOLINT(google-explicit-constructor)
            : r(held(ColVec<NREF, T>{static_cast<T>(held_values)...})) {}
    };

    /**
     * @brief An input target over the control horizon, weighted by Ru in place of the target's
     *        input: held (a ColVec<NU>, one value per input in braces, or the scalar of a single
     *        input) or per move (an InputPreview, MATLAB's mv.target with preview).
     */

    struct InputTarget {
        InputPreview u{}; ///< The target at moves 0 … NC − 1

        constexpr InputTarget(const ColVec<NU, T>& held_value) : u(held_target(held_value)) {} // NOLINT(google-explicit-constructor)
        constexpr InputTarget(const InputPreview& preview) : u(preview) {}                     // NOLINT(google-explicit-constructor)

        /** @brief A single input's target from its scalar, held over the NC moves. */
        constexpr InputTarget(T held_value) // NOLINT(google-explicit-constructor)
            requires(NU == 1)
            : u(held_target(ColVec<NU, T>{held_value})) {}

        /** @brief One value per input, braced at the call (`set_input_target({0.1, 0.2})`), held over the NC moves. */
        template<typename... V>
            requires(NU > 1 && sizeof...(V) == NU && (damp::convertible_to<V, T> && ...))
        constexpr InputTarget(V... held_values) // NOLINT(google-explicit-constructor)
            : u(held_target(ColVec<NU, T>{static_cast<T>(held_values)...})) {}
    };

    /**
     * @brief Time update with the input the plant received over the last period (back_calculate's,
     *        else the last command): the base of the next Δu and the estimator's predict.
     */
    constexpr void predict() {
        applied_.take();
        if constexpr (estimated) {
            estimator.predict(applied_.applied);
        } else if constexpr (NRU > 0) {
            //! The delay states [v(k−1); …; v(k−NRU)] move a period on: the input the plant just received is the newest.
            for (size_t b = NRU - 1; b > 0; --b) {
                for (size_t j = 0; j < NU; ++j) {
                    x_(NYM + (b * NU) + j) = x_(NYM + ((b - 1) * NU) + j);
                }
            }
            for (size_t j = 0; j < NU; ++j) {
                x_(NYM + j) = applied_.applied(j);
            }
        }
    }

    /**
     * @brief Measurement update with @p y: the plant's outputs through the filter, or the plant
     *        state itself without an estimator.
     */
    constexpr void update(const ColVec<NYM, T>& y) {
        if constexpr (estimated) {
            estimator.update(y);
            if constexpr (NDE > 0) {
                estimator.limit_disturbance(disturbance_lim);
            }
        } else if constexpr (NRU == 0) {
            x_ = y;
        } else {
            for (size_t i = 0; i < NYM; ++i) {
                x_(i) = y(i); // the plant's state; the delay states are the law's own
            }
        }
    }

    /**
     * @brief The control move at the current state: u = u_prev + Δu₀
     *
     * The law alone: it solves from the held input and records the command, but does not hold
     * it; @ref predict does, with the input the plant actually received. The QP is posed about
     * the steady-state target of the reference now, r[0]: under the estimated disturbance when
     * the design is offset-free (re-anchoring the estimator on it), under the measured one
     * otherwise. Each step's reference on the outputs is r there, or the target's output for it
     * when the target holds fewer than every output, and each move's input target the target's
     * input for r there, unless set_input_target set one. With d periods of delay the preview the
     * plan reads is r shifted by d.
     *
     * @param r  Reference on the outputs the target holds (held, or a Preview); not a state x_ref
     * @param d  Measured disturbance, held over the horizon (ND > 0; nothing otherwise)
     * @return Input command u (clamped to the configured u and Δu boxes)
     */
    [[nodiscard]] constexpr ColVec<NU, T> control(const Reference& r, const ColVec<ND, T>& d) {
        const Preview& ahead = r.r;
        if constexpr (NDE > 0) {
            static_cast<void>(d); // no measured disturbance alongside an estimated one
            target_.anchor(estimator, r.r[0]);
            const auto x_dev = ::damp::detail::leading<NXM>(estimator.at_command());
            const auto x_a = ::damp::detail::leading<NXM>(estimator.anchor_state());
            const auto d_a = estimator.anchor_disturbance();
            return record(control_impl(ahead, ColVec<NXM, T>(x_dev + x_a), d_a, ColVec<NDE, T>{}, x_a, estimator.anchor_input(), x_dev));
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
     * @param r Reference on the outputs the target holds (held, or a Preview); not a state x_ref
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
            estimator.reset();
        } else if constexpr (NRU > 0) {
            x_ = MeasuredState{}; // no input on its way, as the held input
        }
        last_status_ = design::QPStatus::Success;
        last_iterations_ = 0;
        last_slack_ = T{0};
        last_cost_ = T{0};
    }

    /**
     * @brief Set the input target Ru pulls toward, in place of the target's input (physical units):
     *        held, or per move over the control horizon (InputTarget), until clear_input_target.
     *
     * An operating point of the application's (a preferred valve opening, the cheaper of two
     * actuators): Ru then trades it against the outputs, and an input target off the steady
     * input that holds r leaves an offset.
     */
    constexpr void set_input_target(const InputTarget& u_target) {
        u_target_ = u_target.u;
        follow_target_ = false;
    }

    /** @brief Undo set_input_target: from the next control Ru pulls toward the target's input again, the default. */
    constexpr void clear_input_target() { follow_target_ = true; }


    /**
     * @brief Change the limits on u, Δu and the outputs and states the design bounds (its OutputLimits and StateLimits; MATLAB®'s online constraints)
     *
     * The QP keeps its Hessian, constraint matrix, scaling and softness; the bound vector and the
     * fallback clamp take the new limits from the next control. Refused (false, limits unchanged)
     * when a limit is NaN, a lower limit is above its upper one, or the design has a terminal
     * invariant set, which was built for the design's limits.
     */
    [[nodiscard]] constexpr bool set_limits(const design::MPCLimits<NU, NYC, T, NXC>& lim) {
        if (art.terminal == design::detail::MPCTerminal::invariant_set) {
            return false;
        }
        if (!lim.u_lim.valid() || !lim.du_lim.valid() || !lim.y_lim.valid() || !lim.x_lim.valid()) {
            return false;
        }
        design::detail::mpc_box_bounds(art, lim);
        art.constraints.u_lim = lim.u_lim;
        art.constraints.du_lim = lim.du_lim;
        art.constraints.y_lim = lim.y_lim;
        art.constraints.x_lim = lim.x_lim;
        return true;
    }

    /** @brief The limits in force (the design's, or the last set_limits). */
    [[nodiscard]] constexpr design::MPCLimits<NU, NYC, T, NXC> limits() const {
        const auto& c = art.constraints;
        return {c.u_lim, c.du_lim, c.y_lim, c.x_lim};
    }

    /** @brief Access the solver policy (stateful solvers expose warm-start caches here) */
    [[nodiscard]] constexpr Solver&       qp_solver() { return qp_solver_; }
    [[nodiscard]] constexpr const Solver& qp_solver() const { return qp_solver_; }

    /** @brief The input the last predict held: the base of the next Δu (back_calculate's, else the command before it). */
    [[nodiscard]] constexpr const ColVec<NU, T>& applied() const { return applied_.applied; }
    /** @brief The command the last control returned: the input a predict without back_calculate holds. */
    [[nodiscard]] constexpr const ColVec<NU, T>& last_command() const { return applied_.command; }
    /** @brief The input target of the first move (physical units): the one set_input_target set, else the target's input the last control used (zero before it). */
    [[nodiscard]] constexpr const ColVec<NU, T>& input_target() const { return u_target_[0]; }
    /** @brief The input target of every move over the control horizon, as input_target(). */
    [[nodiscard]] constexpr const InputPreview& input_target_preview() const { return u_target_; }
    /** @brief The status of the last QP solve. */
    [[nodiscard]] constexpr design::QPStatus last_status() const { return last_status_; }
    /** @brief False when built from a failed design: control holds the last input. */
    [[nodiscard]] constexpr bool valid() const { return art.success; }
    /** @brief The iterations of the last QP solve. */
    [[nodiscard]] constexpr size_t last_iterations() const { return last_iterations_; }
    /** @brief Soft-constraint slack ε from the last solve (0 = no output bound was relaxed) */
    [[nodiscard]] constexpr T last_slack() const { return last_slack_; }
    /** @brief QP objective ½ζᵀHζ + fᵀζ from the last solve (for a Lyapunov / monotonicity check) */
    [[nodiscard]] constexpr T last_cost() const { return last_cost_; }

    /** @brief The law's model state: the estimate (absolute; with the register of a delay), or the last measurement. */
    [[nodiscard]] constexpr ColVec<NX, T> state_estimate() const {
        if constexpr (estimated) {
            return ::damp::detail::leading<NX>(estimator.state_estimate());
        } else {
            return x_;
        }
    }

    /** @brief The estimated disturbance d̂ (what the integrating model has absorbed). */
    [[nodiscard]] constexpr ColVec<NDE, T> disturbance_estimate() const
        requires(NDE > 0)
    {
        return estimator.disturbance_estimate();
    }

    /**
     * @brief NP copies of one reference: the held-setpoint case of a preview. Built by copy-construction
     * (no element assignment), so it is a constant expression for a constinit controller.
     */
    [[nodiscard]] static constexpr Preview held(const ColVec<NREF, T>& r) {
        return [&]<size_t... I>(damp::index_sequence<I...>) { return Preview{{(static_cast<void>(I), r)...}}; }(damp::make_index_sequence<NP>{});
    }

    /** @brief NC copies of one input target. */
    [[nodiscard]] static constexpr InputPreview held_target(const ColVec<NU, T>& u) {
        return [&]<size_t... I>(damp::index_sequence<I...>) { return InputPreview{{(static_cast<void>(I), u)...}}; }(damp::make_index_sequence<NC>{});
    }

private:
    using Target = SteadyStateTarget<NXE, NU, ND + NDE, NREF, T>;
    using MeasuredState = damp::conditional_t<estimated, detail::NoMeasuredState, ColVec<NX, T>>;
    using HeldMap = detail::maybe<(NH > 0), Matrix<NY, NREF, T>, detail::OffSlot::held_outputs>;

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

    /** @brief The QP of a design, refused with the design when its filter or target was: valid() is false and control holds its input. */
    template<typename Kalman>
    [[nodiscard]] static constexpr Artifacts artifacts_of(const ResultWith<Kalman>& result) {
        Artifacts a = result.mpc;
        if (!result.success && a.success) {
            a.success = false;
            a.status = result.status;
        }
        return a;
    }

    /** @brief The runtime target from its design, refused only by this constructor. */

    [[nodiscard]] static constexpr Target target_of(const TargetResult& design) {
        Target t{};
        t.target_d = design.target_d;
        t.target_r = design.target_r;
        return t;
    }

    /**
     * @brief The target's outputs per unit reference, C·Nₓ (NH > 0): a step of the preview moves
     *        the QP's reference on every output by it, the held outputs by the step of r.
     */
    [[nodiscard]] static constexpr HeldMap held_map(const Artifacts& a, const Target& t) {
        if constexpr (NH > 0) {
            Matrix<NY, NREF, T> m{};
            for (size_t i = 0; i < NY; ++i) {
                for (size_t j = 0; j < NREF; ++j) {
                    T sum{0};
                    for (size_t k = 0; k < NXM; ++k) {
                        sum += a.C(i, k) * t.target_r(k, j);
                    }
                    m(i, j) = sum;
                }
            }
            return m;
        } else {
            static_cast<void>(a);
            static_cast<void>(t);
            return HeldMap{};
        }
    }


    /** @brief The law's state the command decided now acts on (absolute): the estimate, or the measurement. */
    [[nodiscard]] constexpr ColVec<NXM, T> law_state() const {
        if constexpr (estimated) {
            return ::damp::detail::leading<NXM>(estimator.at_command());
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
            if constexpr (NH > 0) {
                //! The target's output for r[kb]: y_t moved by C·Nₓ·(r[kb] − r[0]), the held outputs at r[kb], the others where it puts them
                const ColVec<NY, T> dy = y_r_ * ColVec<NREF, T>(r[kb] - r[0]);
                for (size_t i = 0; i < NY; ++i) {
                    r_dev((kb * NY) + i) = dy(i) / art.scale_y(i);
                }
            } else {
                for (size_t i = 0; i < NY; ++i) {
                    r_dev((kb * NY) + i) = (r[kb](i) - y_t(i)) / art.scale_y(i);
                }
            }
        }
        //! Ru's anchor at each move: the target's input for r there (u_t moved by Nᵤ·(r[jm] − r[0])), or the input target set
        ColVec<NZ, T> ut_dev{};
        for (size_t jm = 0; jm < NC; ++jm) {
            for (size_t j = 0; j < NU; ++j) {
                if (follow_target_) {
                    T du{0};
                    for (size_t c = 0; c < NREF; ++c) {
                        du += target_.target_r(NXE + j, c) * (r[jm](c) - r[0](c));
                    }
                    u_target_[jm](j) = u_t(j) + du;
                    ut_dev((jm * NU) + j) = du / art.scale_u(j);
                } else {
                    ut_dev((jm * NU) + j) = (u_target_[jm](j) - u_t(j)) / art.scale_u(j);
                }
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
        constexpr size_t r_term = Artifacts::NI_box;
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
        for (size_t i = 0; i < NYC; ++i) {
            if (art.constraints.y_lower_ecr(i) == T{0} || art.constraints.y_upper_ecr(i) == T{0}) {
                return true;
            }
        }
        for (size_t i = 0; i < NXC; ++i) {
            if (art.constraints.x_lower_ecr(i) == T{0} || art.constraints.x_upper_ecr(i) == T{0}) {
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
    [[no_unique_address]] MeasuredState  x_{};                 ///< The last measured state (no storage with an estimator)
    Target                               target_{};            ///< The maps to the steady-state target
    [[no_unique_address]] HeldMap        y_r_{};               ///< The target's outputs per unit reference (NH > 0)
    ::damp::detail::AppliedRecord<NU, T> applied_{};           ///< The last command, and the held input (the base of Δu)
    InputPreview                         u_target_{};          ///< The input target of each move: the target's input the last control used, or the one set
    bool                                 follow_target_{true}; ///< Whether Ru pulls toward the target's input (no set_input_target since clear_input_target)
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
 * The tick is update(x) then control(r), or step(r, x) for both; @p r is the reference of the
 * outputs the target holds (held, or a Preview over the horizon). A design with a measured disturbance takes it
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
    typename Solver = design::WarmStartActiveSetSolver<(NC * NU) + 1, (4 * NC * NU) + 1, T>,
    size_t NT_max = 0,
    size_t NXC = 0,
    size_t NYC = 0,
    size_t NRU = 0,
    size_t NH = 0>
class MPC : private detail::MPCCore<NX, NU, NY, NP, NC, ND, T, Solver, NT_max, NoEstimator, NXC, NYC, NRU, NH> {
    using Core = detail::MPCCore<NX, NU, NY, NP, NC, ND, T, Solver, NT_max, NoEstimator, NXC, NYC, NRU, NH>;

public:
    using Reference = Core::Reference;                     ///< r: held, or a Preview
    using Preview = Core::Preview;                         ///< The reference at steps 1 … NP
    using InputPreview = Core::InputPreview;               ///< The input target per move
    using InputTarget = Core::InputTarget;                 ///< The input target: held, braced, scalar or per move
    using Result = Core::template ResultWith<NoEstimator>; ///< The design it is built from
    static constexpr size_t NZ1 = Core::NZ1;               ///< QP decision size (moves + slack)
    static constexpr size_t NI1 = Core::NI1;               ///< QP constraint rows

    constexpr MPC() = default;

    /** @brief From design::mpc. A failed design is a compile error when constant-evaluated. */
    constexpr explicit MPC(const Result& result) : Core(result) {}

    /** @brief From design::mpc with a specific solver instance (for stateful solver policies). */
    constexpr MPC(const Result& result, const Solver& solver) : Core(result, solver) {}

    /** @brief Hold the measured state @p x for the next control: the plant's (a design with an InputDelay holds the delay states itself). */
    constexpr void update(const ColVec<Core::NYM, T>& x) { Core::update(x); }

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
    [[nodiscard]] constexpr ColVec<NU, T> step(const Reference& r, const ColVec<Core::NYM, T>& x)
        requires(ND == 0)
    {
        return Core::step(r, ColVec<0, T>{}, x);
    }

    /** @brief One tick: predict(), update(@p x), then control(@p r, @p d). */
    [[nodiscard]] constexpr ColVec<NU, T> step(const Reference& r, const ColVec<ND, T>& d, const ColVec<Core::NYM, T>& x)
        requires(ND > 0)
    {
        return Core::step(r, d, x);
    }

    using Core::applied;
    using Core::back_calculate;
    using Core::clear_input_target;
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
    using Core::valid;
};

/** @brief MPC law{res}: the sizes, scalar, terminal capacity, bounded states and held outputs of design::mpc's result. */
template<size_t NXM, size_t NU, size_t NY, size_t NP, size_t NC, size_t ND, typename T, size_t NT, size_t NXC, size_t NYC, size_t NRU, size_t NZT, size_t NH>
MPC(const design::MPCResult<design::MPCArtifacts<NXM, NU, NY, NP, NC, ND, T, NT, NXC, NYC>, NoEstimator, design::SteadyStateTargetResult<NXM, NU, ND, NZT, T>, NRU, NH>&)
    -> MPC<NXM, NU, NY, NP, NC, ND, T, design::WarmStartActiveSetSolver<(NC * NU) + 1, design::MPCArtifacts<NXM, NU, NY, NP, NC, ND, T, NT, NXC, NYC>::NI1, T>, NT, NXC, NYC, NRU, NH>;

/** @brief MPC law{res, solver}: as MPC law{res}, with that solver. */
template<size_t NXM, size_t NU, size_t NY, size_t NP, size_t NC, size_t ND, typename T, size_t NT, size_t NXC, size_t NYC, size_t NRU, size_t NZT, size_t NH, typename Solver>
MPC(const design::MPCResult<design::MPCArtifacts<NXM, NU, NY, NP, NC, ND, T, NT, NXC, NYC>, NoEstimator, design::SteadyStateTargetResult<NXM, NU, ND, NZT, T>, NRU, NH>&, const Solver&)
    -> MPC<NXM, NU, NY, NP, NC, ND, T, Solver, NT, NXC, NYC, NRU, NH>;

namespace detail {

/** @brief The runtime an MPC design runs as, in the scalar @p T: none for other results (OFMPC below adds its own). */
template<typename R, typename T>
struct mpc_runtime {};

/** @brief design::mpc runs as MPC. */
template<size_t NXM, size_t NU, size_t NY, size_t NP, size_t NC, size_t ND, typename U, size_t NT, size_t NXC, size_t NYC, size_t NRU, size_t NZT, size_t NH, typename T>
struct mpc_runtime<design::MPCResult<design::MPCArtifacts<NXM, NU, NY, NP, NC, ND, U, NT, NXC, NYC>, NoEstimator, design::SteadyStateTargetResult<NXM, NU, ND, NZT, U>, NRU, NH>, T> {
    using type = MPC<NXM, NU, NY, NP, NC, ND, T, design::WarmStartActiveSetSolver<(NC * NU) + 1, design::MPCArtifacts<NXM, NU, NY, NP, NC, ND, T, NT, NXC, NYC>::NI1, T>, NT, NXC, NYC, NRU, NH>; ///< MPC
};

} // namespace detail

namespace detail {

/** @brief An MPC design. */
template<const auto& Res, typename T>
concept mpc_design = requires { typename mpc_runtime<damp::remove_cvref_t<decltype(Res)>, T>::type; };

/** @brief An MPC design with a filter: OFMPC, whose filter can run on the design's constants. */
template<const auto& Res, typename T>
concept filtered_mpc_design = mpc_design<Res, T> && mpc_runtime<damp::remove_cvref_t<decltype(Res)>, T>::filtered;

} // namespace detail

/** @brief Law<kRes, T> of a `constexpr` OFMPC design: the same law, its filter on the design's matrices as constants of the type. The QP's data stays in members: set_limits rewrites it. */
template<const auto& Res, typename T>
    requires(detail::filtered_mpc_design<Res, T> && detail::constant_design<Res>)
struct Law<Res, T> : detail::mpc_runtime<damp::remove_cvref_t<decltype(Res)>, T>::template on_constants<Res> {
    using Base = detail::mpc_runtime<damp::remove_cvref_t<decltype(Res)>, T>::template on_constants<Res>; ///< The runtime on constants
    using Runtime = detail::mpc_runtime<damp::remove_cvref_t<decltype(Res)>, T>::type;                    ///< The runtime with members: the type of a member a redesign assigns
    static_assert(detail::constant_design_runs(Res.success, Res.status));                                 // a failed design does not compile: the error names its status

    /** @brief The design's law. */
    constexpr Law() : Base{Res.template as<T>()} {}
};

/** @brief Law<kRes, T> of an MPC design: MPC or OFMPC, default-constructed from the design. */
template<const auto& Res, typename T>
    requires detail::mpc_design<Res, T>
struct Law<Res, T> : detail::mpc_runtime<damp::remove_cvref_t<decltype(Res)>, T>::type {
    using Runtime = detail::mpc_runtime<damp::remove_cvref_t<decltype(Res)>, T>::type; ///< The runtime the design runs as

    /** @brief The design's law. */
    constexpr Law() : Runtime{Res.template as<T>()} {}

    using Runtime::operator=;
};

} // namespace damp

// The offset-free MPC: design::mpc with a noise and a disturbance spec runs as OFMPC.
namespace damp {

/**
 * @ingroup discrete_controllers
 * @brief Offset-free constrained linear MPC: a Kalman filter on [x; d] and the QP about the
 *        steady-state target under d̂ (design::mpc with a KalmanNoise and a disturbance spec)
 *
 * The estimator is a FixedGainObserver with design::kalman's fixed gain, holding deviations
 * from the target. The tick is update(y) then control(r), or step(r, y) for both; @p r is the
 * reference of the outputs the target holds (held, or a Preview). back_calculate(u_sat) records what went
 * out, predict() moves the observer with it. update() clamps d̂ into disturbance_lim, a member
 * (unbounded by default), in d's units, as OFLQG's. A failed design does not compile when
 * constant-evaluated.
 *
 * @note Compare with MATLAB®'s mpcmove(mpcobj, xc, ym, r) with its default disturbance model.
 */
template<size_t NX, size_t NU, size_t NY, size_t NP, size_t NC, typename T, typename Solver, size_t NT_max, typename Estimator, size_t NXC = 0, size_t NYC = 0, size_t NH = 0>
class OFMPC : private detail::MPCCore<NX, NU, NY, NP, NC, 0, T, Solver, NT_max, Estimator, NXC, NYC, 0, NH> {
    using Core = detail::MPCCore<NX, NU, NY, NP, NC, 0, T, Solver, NT_max, Estimator, NXC, NYC, 0, NH>;

public:
    using Reference = Core::Reference;       ///< r: held, or a Preview
    using Preview = Core::Preview;           ///< The reference at steps 1 … NP
    using InputPreview = Core::InputPreview; ///< The input target per move
    using InputTarget = Core::InputTarget;   ///< The input target: held, braced, scalar or per move

    static constexpr size_t NYM = Core::NYM; ///< Measured outputs
    static constexpr size_t NDE = Core::NDE; ///< Estimated disturbance states
    static constexpr size_t NZ1 = Core::NZ1; ///< QP decision size (moves + slack)
    static constexpr size_t NI1 = Core::NI1; ///< QP constraint rows
    static_assert(NDE > 0, "OFMPC runs design::mpc with a KalmanNoise and a disturbance spec: the filter estimates a disturbance.");

    constexpr OFMPC() = default;

    /** @brief From design::mpc with a KalmanNoise and a disturbance spec (its filter with any noise inputs). A failed design is a compile error when constant-evaluated. */
    template<typename Kalman>
    constexpr explicit OFMPC(const design::MPCResult<typename Core::Artifacts, Kalman, typename Core::TargetResult, 0, NH>& result) : Core(result) {}

    /** @brief From that design with a specific solver instance (for stateful solver policies). */
    template<typename Kalman>
    constexpr OFMPC(const design::MPCResult<typename Core::Artifacts, Kalman, typename Core::TargetResult, 0, NH>& result, const Solver& solver) : Core(result, solver) {}

    /** @brief The filter's measurement update with @p y. */
    constexpr void update(const ColVec<NYM, T>& y) { Core::update(y); }

    /** @brief The move about the target for @p r under d̂: u = u_prev + Δu₀, clamped to the u and Δu boxes. */
    [[nodiscard]] constexpr ColVec<NU, T> control(const Reference& r) { return Core::control(r, ColVec<0, T>{}); }

    /** @brief One tick: predict(), update(@p y), then control(@p r). */
    [[nodiscard]] constexpr ColVec<NU, T> step(const Reference& r, const ColVec<NYM, T>& y) { return Core::step(r, ColVec<0, T>{}, y); }

    using Core::applied;
    using Core::back_calculate;
    using Core::clear_input_target;
    using Core::disturbance_estimate;
    using Core::disturbance_lim;
    using Core::estimator;
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
    using Core::valid;
};

/** @brief OFMPC law{res}: the sizes, scalar, filter, target and held outputs of an offset-free design::mpc result. */
template<size_t NXM, size_t NU, size_t NY, size_t NP, size_t NC, size_t ND, typename T, size_t NT, size_t NXC, size_t NYC, size_t NXK, size_t NYM, size_t NW, size_t NV, size_t NXE, size_t NZT, size_t NH>
OFMPC(const design::MPCResult<design::MPCArtifacts<NXM, NU, NY, NP, NC, ND, T, NT, NXC, NYC>, design::KalmanResult<NXK, NU, NYM, T, NW, NV, ND>, design::SteadyStateTargetResult<NXE, NU, ND, NZT, T>, 0, NH>&)
    -> OFMPC<NXM, NU, NY, NP, NC, T, design::WarmStartActiveSetSolver<(NC * NU) + 1, design::MPCArtifacts<NXM, NU, NY, NP, NC, ND, T, NT, NXC, NYC>::NI1, T>, NT, FixedGainObserver<NXK, NU, NYM, T, ND>, NXC, NYC, NH>;

/** @brief OFMPC law{res, solver}: as OFMPC law{res}, with that solver. */
template<size_t NXM, size_t NU, size_t NY, size_t NP, size_t NC, size_t ND, typename T, size_t NT, size_t NXC, size_t NYC, size_t NXK, size_t NYM, size_t NW, size_t NV, size_t NXE, size_t NZT, size_t NH, typename Solver>
OFMPC(const design::MPCResult<design::MPCArtifacts<NXM, NU, NY, NP, NC, ND, T, NT, NXC, NYC>, design::KalmanResult<NXK, NU, NYM, T, NW, NV, ND>, design::SteadyStateTargetResult<NXE, NU, ND, NZT, T>, 0, NH>&, const Solver&) -> OFMPC<NXM, NU, NY, NP, NC, T, Solver, NT, FixedGainObserver<NXK, NU, NYM, T, ND>, NXC, NYC, NH>;

namespace detail {

/** @brief An offset-free design::mpc runs as OFMPC. */
template<size_t NXM, size_t NU, size_t NY, size_t NP, size_t NC, size_t ND, typename U, size_t NT, size_t NXC, size_t NYC, size_t NXK, size_t NYM, size_t NW, size_t NV, size_t NXE, size_t NZT, size_t NH, typename T>
struct mpc_runtime<design::MPCResult<design::MPCArtifacts<NXM, NU, NY, NP, NC, ND, U, NT, NXC, NYC>, design::KalmanResult<NXK, NU, NYM, U, NW, NV, ND>, design::SteadyStateTargetResult<NXE, NU, ND, NZT, U>, 0, NH>, T> {
    using type = OFMPC<NXM, NU, NY, NP, NC, T, design::WarmStartActiveSetSolver<(NC * NU) + 1, design::MPCArtifacts<NXM, NU, NY, NP, NC, ND, T, NT, NXC, NYC>::NI1, T>, NT, FixedGainObserver<NXK, NU, NYM, T, ND>, NXC, NYC, NH>; ///< OFMPC
    static constexpr bool filtered = true;                                                                                                                                                                                         ///< Its filter can run on the design's constants
    template<const auto& Res>
    using on_constants = OFMPC<NXM, NU, NY, NP, NC, T, design::WarmStartActiveSetSolver<(NC * NU) + 1, design::MPCArtifacts<NXM, NU, NY, NP, NC, ND, T, NT, NXC, NYC>::NI1, T>, NT, FixedGainObserver<NXK, NU, NYM, T, ND, ConstantObserver<Res, T, NXK - ND>>, NXC, NYC, NH>; ///< OFMPC with its filter on the design's constants
};

} // namespace detail

} // namespace damp
