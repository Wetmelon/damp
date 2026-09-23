// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file mpc.hpp
 * @brief Constrained linear MPC: condensed-QP synthesis (`design::state_mpc`)
 *        and the fixed-budget runtime controller (`damp::MPC`).
 *
 * Finite-horizon constrained control of a discrete-time plant
 * @f$ x_{k+1} = A x_k + B u_k + B_d d_k,\; y_k = C x_k @f$ (roadmap #14), where
 * @f$ u @f$ are the manipulated inputs and @f$ d @f$ optional measured
 * disturbances (feedforward inputs the controller knows but does not command).
 * The move increments @f$ \Delta u @f$ are the decision variables (the
 * "velocity form"): the state is augmented as @f$ \tilde x = [x;\, u_{prev};\, d] @f$
 * so the held input and the measured disturbance are carried by the prediction
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
 * guaranteed-stability short-horizon setup). Subject to box constraints on
 * @f$ u @f$, @f$ \Delta u @f$ (hard), and @f$ y @f$ (soft by default) over the
 * horizon.
 *
 * Soft constraints (ECR). A single slack @f$ \varepsilon \ge 0 @f$ relaxes
 * the output rows: @f$ y \le y_{max} + V\varepsilon @f$ with per-channel
 * softness coefficients V (`MPCConstraints::y_max_ecr`/`y_min_ecr`, default 1)
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
 * Measured disturbances are not rescaled (unit scale).
 *
 * Weighting a state that is not an output: add it as an extra
 * (unmeasured) output row in C and weight that row in Qy — the toolbox-standard
 * recipe. E.g. to damp velocity x₁ of a position-output plant, use
 * C = [1 0; 0 1], NY = 2, r = [r_pos; 0], Qy = diag(q_pos, q_vel). There is
 * deliberately no separate state weight (MATLAB®'s linear MPC has none either).
 *
 * r already is the tracking reference (NY). There is no LQR-style
 * `control(x_ref, x)` overload: the QP has no ‖x − x_ref‖ term (it does not
 * force x → 0 except as a consequence of a feasible y → r), and that
 * two-argument signature would collide with `control(r, x)` whenever NY = NX.
 *
 * Synthesis condenses the problem: predictions @f$ Y = \Phi \tilde x_0 +
 * \Gamma Z @f$ over the horizon give a dense strictly convex QP in
 * @f$ [Z; \varepsilon] @f$, @f$ Z = [\Delta u_0; \dots; \Delta u_{NC-1}] @f$,
 * solved each tick by the Goldfarb–Idnani active-set solver (design/qp.hpp)
 * under a hard iteration budget. Everything is fixed-size and allocation-free;
 * per decisions.md D16 this header ships behind `workbench.hpp` (include it
 * directly for on-target builds — nothing here allocates).
 *
 * Feedthrough is not supported: plants must have @f$ D = 0 @f$ (synthesis
 * fails otherwise). Offset-free tracking under plant/model mismatch or
 * *unmeasured* disturbances additionally needs the disturbance-model
 * augmentation planned in roadmap #14 — the velocity form alone rejects only
 * what its integral action sees through the measured state.
 *
 * Unbounded constraint entries use the `numeric_limits` sentinels
 * (`lowest()`/`max()`, never infinity — the library builds with
 * -ffinite-math-only); sentinel rows cost nothing in the solver.
 *
 * @note Compare with MATLAB®'s mpc(plant, Ts, p, m, weights) + mpcmove, with
 *       Weights.OV/MV/MVRate/ECR, ManipulatedVariables Targets, and measured
 *       disturbances via setmpcsignals.
 * @see design/qp.hpp for the QP solver
 * @see García, Prett & Morari, "Model Predictive Control: Theory and Practice —
 *      A Survey," Automatica 25(3), 1989, https://doi.org/10.1016/0005-1098(89)90002-2
 * @see Maciejowski, "Predictive Control with Constraints," Prentice Hall, 2002
 *      (condensed Δu formulation, ch. 2–3; soft constraints §3.4)
 * @see Rawlings, Mayne & Diehl, "Model Predictive Control: Theory, Computation,
 *      and Design," 2nd ed., Nob Hill, 2017
 *
 * @code
 * // Double integrator, Ts = 0.1: track r with |u| ≤ 0.5, |Δu| ≤ 0.2.
 * constexpr StateSpace<2, 1, 1> sys{
 *     .A = {{1.0, 0.1}, {0.0, 1.0}},
 *     .B = {{0.005}, {0.1}},
 *     .C = {{1.0, 0.0}},
 *     .Ts = 0.1,
 * };
 * constexpr design::MPCWeights<1, 1> weights{};   // Qy = I, RΔu = I
 * constexpr design::MPCConstraints<1, 1> limits{
 *     .u_min = ColVec<1>{-0.5}, .u_max = ColVec<1>{0.5},
 *     .du_min = ColVec<1>{-0.2}, .du_max = ColVec<1>{0.2},
 * };
 * constexpr auto art = design::state_mpc<15, 5>(sys, weights, limits);
 * static_assert(art.success);
 * MPC controller{art.as<float>()};
 * // each tick: u = controller.control(r, x);
 * @endcode
 */

#include <cstddef>
#include <limits>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/design/qp.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/block.hpp"
#include "damp/matrix/decomposition.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace design {

namespace detail {

// filled_vector / convert_bounds (sentinel machinery) live in design/qp.hpp's
// detail — shared with the MHE and any other bound-carrying design function.

/// Scale factor from a constraint range: (max − min) when both finite, else 1.
template<typename T>
[[nodiscard]] constexpr T scale_from_range(T lo, T hi) {
    constexpr T thresh = std::numeric_limits<T>::max() / T{4};
    if (lo <= -thresh || lo >= thresh || hi >= thresh || hi <= -thresh) {
        return T{1};
    }
    const T span = hi - lo;
    return (span > T{0}) ? span : T{1};
}

/// Divide a bound by a scale factor, leaving unbounded sentinels untouched.
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
 * @struct MPCWeights
 * @brief Cost weights for mpc()
 *
 * Weights act on the internally scaled (dimensionless) signals, so they are
 * unit-free once constraint ranges are given — MATLAB®'s ScaleFactor semantics.
 * All weight matrices are symmetric positive semidefinite; the QP Hessian must
 * come out positive definite, which holds whenever RΔu ≻ 0 or Ru ≻ 0.
 *
 * @note Mirrors the MATLAB® mpc object's Weights.OV / Weights.MVRate /
 *       Weights.MV (+ MV Targets) / Weights.ECR.
 */
template<size_t NU, size_t NY, typename T = double>
struct MPCWeights {
    Matrix<NY, NY, T> Qy = Matrix<NY, NY, T>::identity();  ///< Output tracking weight, steps 1..NP−1
    Matrix<NU, NU, T> Rdu = Matrix<NU, NU, T>::identity(); ///< Move-increment weight on Δu
    Matrix<NU, NU, T> Ru{};                                ///< Input-deviation weight on (u − u_target) (default 0)
    Matrix<NY, NY, T> Qy_terminal = Qy;                    ///< Output weight on the final step NP (default: Qy)
    ColVec<NU, T>     u_target{};                          ///< Preferred input operating point (physical units)
    T                 ecr_weight{static_cast<T>(1.0e5)};   ///< Soft-constraint slack penalty ρ_ε

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
 * prediction horizon and are soft by default (ECR coefficient 1) — set a
 * channel's ECR to 0 to make it hard. Larger ECR = that channel's bound is
 * relaxed more readily relative to the others.
 *
 * @note Mirrors the MATLAB® mpc object's ManipulatedVariables Min/Max/RateMin/
 *       RateMax and OutputVariables Min/Max (+ MinECR/MaxECR).
 */
template<size_t NU, size_t NY, typename T = double>
struct MPCConstraints {
    ColVec<NU, T> u_min = detail::filled_vector<NU>(std::numeric_limits<T>::lowest());  ///< Input lower bound
    ColVec<NU, T> u_max = detail::filled_vector<NU>(std::numeric_limits<T>::max());     ///< Input upper bound
    ColVec<NU, T> du_min = detail::filled_vector<NU>(std::numeric_limits<T>::lowest()); ///< Move lower bound
    ColVec<NU, T> du_max = detail::filled_vector<NU>(std::numeric_limits<T>::max());    ///< Move upper bound
    ColVec<NY, T> y_min = detail::filled_vector<NY>(std::numeric_limits<T>::lowest());  ///< Output lower bound
    ColVec<NY, T> y_max = detail::filled_vector<NY>(std::numeric_limits<T>::max());     ///< Output upper bound
    ColVec<NY, T> y_min_ecr = detail::filled_vector<NY>(T{1});                          ///< Softness of y_min rows (0 = hard)
    ColVec<NY, T> y_max_ecr = detail::filled_vector<NY>(T{1});                          ///< Softness of y_max rows (0 = hard)

    template<typename U>
    [[nodiscard]] constexpr MPCConstraints<NU, NY, U> as() const {
        return MPCConstraints<NU, NY, U>{
            detail::convert_bounds<U>(u_min),
            detail::convert_bounds<U>(u_max),
            detail::convert_bounds<U>(du_min),
            detail::convert_bounds<U>(du_max),
            detail::convert_bounds<U>(y_min),
            detail::convert_bounds<U>(y_max),
            ColVec<NY, U>{y_min_ecr.template as<U>()},
            ColVec<NY, U>{y_max_ecr.template as<U>()},
        };
    }
};

/**
 * @struct MPCArtifacts
 * @brief Condensed-QP data produced by mpc(), consumed by damp::MPC
 *
 * Holds the dense QP in [Z; ε] (move increments + soft-constraint slack):
 * Hessian H, gradient maps (f = F_x·x̃0 − F_r·r' − F_ut·u'_target with
 * x̃0 = [x; u'_prev; d]), the stacked constraint matrix A_con with its
 * state-dependent bounds (b = b_bound + B_map·x̃0), the per-channel scale
 * factors, and the physical clamping limits. Primes denote internally scaled
 * quantities. Use .as<float>() for embedded deployment (bound sentinels
 * convert safely).
 */
template<size_t NX, size_t NU, size_t NY, size_t NP, size_t NC, size_t ND = 0, typename T = double>
struct MPCArtifacts {
    static constexpr size_t NXA = NX + NU + ND;            ///< Augmented state size [x; u'_prev; d]
    static constexpr size_t NZ = NC * NU;                  ///< Move-increment decision variables
    static constexpr size_t NZ1 = NZ + 1;                  ///< QP decision size (moves + slack ε)
    static constexpr size_t NI = (4 * NZ) + (2 * NP * NY); ///< Box-constraint rows
    static constexpr size_t NI1 = NI + 1;                  ///< QP constraint rows (+ ε ≥ 0)

    Matrix<NZ1, NZ1, T>       H{};                 ///< QP Hessian (positive definite)
    Matrix<NZ, NXA, T>        F_x{};               ///< Gradient state map
    Matrix<NZ, NY, T>         F_r{};               ///< Gradient reference map (scaled reference)
    Matrix<NZ, NU, T>         F_ut{};              ///< Gradient input-target map (scaled target)
    Matrix<NI1, NZ1, T>       A_con{};             ///< Constraint rows: A_con·[Z;ε] ≤ b_bound + B_map·x̃0
    ColVec<NI1, T>            b_bound{};           ///< Constant bound term (sentinel rows disabled)
    Matrix<NI1, NXA, T>       B_map{};             ///< State-dependent bound term
    ColVec<NU, T>             scale_u{};           ///< Input scale factors
    ColVec<NY, T>             scale_y{};           ///< Output scale factors
    ColVec<NU, T>             u_target{};          ///< Initial input target (physical units)
    MPCConstraints<NU, NY, T> constraints{};       ///< Original limits (runtime clamp fallback)
    size_t                    max_qp_iterations{}; ///< Per-tick QP iteration budget
    bool                      success{false};      ///< true if the synthesis validated

    template<typename U>
    [[nodiscard]] constexpr MPCArtifacts<NX, NU, NY, NP, NC, ND, U> as() const {
        return MPCArtifacts<NX, NU, NY, NP, NC, ND, U>{
            H.template as<U>(),
            F_x.template as<U>(),
            F_r.template as<U>(),
            F_ut.template as<U>(),
            A_con.template as<U>(),
            detail::convert_bounds<U>(b_bound),
            B_map.template as<U>(),
            ColVec<NU, U>{scale_u.template as<U>()},
            ColVec<NY, U>{scale_y.template as<U>()},
            ColVec<NU, U>{u_target.template as<U>()},
            constraints.template as<U>(),
            max_qp_iterations,
            success,
        };
    }
};

namespace detail {

template<size_t NP, size_t NC, size_t ND, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr MPCArtifacts<NX, NU, NY, NP, NC, ND, T> mpc_impl(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NX, ND, T>&                 Bd,
    const MPCWeights<NU, NY, T>&             weights,
    const MPCConstraints<NU, NY, T>&         constraints,
    size_t                                   max_qp_iterations
) {
    static_assert(NP >= 1, "Prediction horizon must be at least one step.");
    static_assert(NC >= 1 && NC <= NP, "Control horizon must satisfy 1 <= NC <= NP.");

    constexpr size_t NXA = NX + NU + ND;
    constexpr size_t NZ = NC * NU;
    constexpr size_t NZ1 = NZ + 1;
    constexpr size_t NYS = NP * NY;
    constexpr size_t NI = (4 * NZ) + (2 * NYS);

    MPCArtifacts<NX, NU, NY, NP, NC, ND, T> art{};
    art.constraints = constraints;
    art.u_target = weights.u_target;
    art.max_qp_iterations = (max_qp_iterations == 0) ? 10 * (NZ1 + NI + 1) : max_qp_iterations;

    if (!(sys.Ts > T{0})) {
        return art; // discrete plants only — discretize first (c2d)
    }

    for (size_t r = 0; r < NY; ++r) {
        for (size_t c = 0; c < NU; ++c) {
            if (damp::abs(sys.D(r, c)) > damp::default_tol<T>()) {
                return art; // feedthrough unsupported
            }
        }
    }

    if (!(weights.ecr_weight > T{0})) {
        return art; // slack penalty must be positive to keep the Hessian PD
    }

    // Automatic per-channel scaling from constraint ranges (MATLAB® ScaleFactor).
    for (size_t j = 0; j < NU; ++j) {
        art.scale_u(j) = detail::scale_from_range(constraints.u_min(j), constraints.u_max(j));
    }
    for (size_t i = 0; i < NY; ++i) {
        art.scale_y(i) = detail::scale_from_range(constraints.y_min(i), constraints.y_max(i));
    }

    // Scaled model: u = Su·u', y = Sy·y'  →  B' = B·Su, C' = Sy⁻¹·C.
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

    // Augmented Δu-form system over x̃ = [x; u'_prev; d].
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

    // Prediction: Y = Phi x̃0 + Gamma Z, with CAB[m] = C̃ Ã^m B̃.
    Matrix<NYS, NXA, T>                Phi{};
    Matrix<NYS, NZ, T>                 Gamma{};
    damp::array<Matrix<NY, NU, T>, NP> CAB{};
    Matrix<NXA, NXA, T>                Apow = Matrix<NXA, NXA, T>::identity();
    for (size_t m = 0; m < NP; ++m) {
        CAB[m] = Ca * (Apow * Ba);
        Apow = Apow * Aa;
        Phi.template block<NY, NXA>(m * NY, 0) = Ca * Apow; // C̃ Ã^{m+1}
    }
    for (size_t k = 1; k <= NP; ++k) {
        const size_t jmax = (k < NC) ? k : NC;
        for (size_t j = 0; j < jmax; ++j) {
            Gamma.template block<NY, NU>((k - 1) * NY, j * NU) = CAB[k - 1 - j];
        }
    }

    // Q̄-weighted stacks (Qy on steps 1..NP−1, Qy_terminal on step NP).
    Matrix<NYS, NZ, T>  QG{};
    Matrix<NYS, NXA, T> QPhi{};
    Matrix<NYS, NY, T>  Qstack{};
    for (size_t kb = 0; kb < NP; ++kb) {
        const Matrix<NY, NY, T>& Qk = (kb == NP - 1) ? weights.Qy_terminal : weights.Qy;
        QG.template block<NY, NZ>(kb * NY, 0) = Qk * Gamma.template block<NY, NZ>(kb * NY, 0).to_matrix();
        QPhi.template block<NY, NXA>(kb * NY, 0) = Qk * Phi.template block<NY, NXA>(kb * NY, 0).to_matrix();
        Qstack.template block<NY, NY>(kb * NY, 0) = Qk;
    }

    // Hessian: 2(ΓᵀQ̄Γ + R̄Δu + TᵀR̄uT) plus the slack corner 2ρ_ε. With T the
    // block lower-triangular ones map (Z → u sequence),
    // (TᵀR̄uT)_{ij} = (NC − max(i,j))·Ru.
    Matrix<NZ, NZ, T> Hm = Gamma.transpose() * QG;
    for (size_t i = 0; i < NC; ++i) {
        for (size_t j = 0; j < NC; ++j) {
            const size_t m = (i > j) ? i : j;
            Hm.template block<NU, NU>(i * NU, j * NU) += static_cast<T>(NC - m) * weights.Ru;
        }
        Hm.template block<NU, NU>(i * NU, i * NU) += weights.Rdu;
    }
    art.H.template block<NZ, NZ>(0, 0) = T{2} * Hm;
    art.H(NZ, NZ) = T{2} * weights.ecr_weight;

    // Gradient maps: f = F_x x̃0 − F_r r' − F_ut u'_target. The TᵀR̄uS term hits
    // the u'_prev columns of x̃0 (S stacks [0 I 0] over the control horizon);
    // the target map F_ut shares the same (NC−i)·Ru block-column structure.
    Matrix<NZ, NXA, T> Fx = Gamma.transpose() * QPhi;
    for (size_t i = 0; i < NC; ++i) {
        Fx.template block<NU, NU>(i * NU, NX) += static_cast<T>(NC - i) * weights.Ru;
        art.F_ut.template block<NU, NU>(i * NU, 0) = T{2} * static_cast<T>(NC - i) * weights.Ru;
    }
    art.F_x = T{2} * Fx;
    art.F_r = T{2} * (Gamma.transpose() * Qstack);

    // Constraint stack: A_con [Z; ε] ≤ b_bound + B_map x̃0.
    // Rows: [0,NZ) Δu' ≤ du'_max · [NZ,2NZ) −Δu' ≤ −du'_min ·
    //       [2NZ,3NZ) u' ≤ u'_max · [3NZ,4NZ) −u' ≤ −u'_min ·
    //       [4NZ,4NZ+NYS) y' − V·ε ≤ y'_max · [4NZ+NYS,NI) −y' − V·ε ≤ −y'_min ·
    //       row NI: −ε ≤ 0.
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
        for (size_t j = 0; j < NU; ++j) {
            const T su = art.scale_u(j);
            art.b_bound(r_du_hi + (i * NU) + j) = detail::scaled_bound(constraints.du_max(j), su);
            art.b_bound(r_du_lo + (i * NU) + j) = detail::scaled_bound(-constraints.du_min(j), su);
            art.b_bound(r_u_hi + (i * NU) + j) = detail::scaled_bound(constraints.u_max(j), su);
            art.b_bound(r_u_lo + (i * NU) + j) = detail::scaled_bound(-constraints.u_min(j), su);
        }
    }
    art.A_con.template block<NYS, NZ>(r_y_hi, 0) = Gamma;
    art.A_con.template block<NYS, NZ>(r_y_lo, 0) = -Gamma;
    art.B_map.template block<NYS, NXA>(r_y_hi, 0) = -Phi;
    art.B_map.template block<NYS, NXA>(r_y_lo, 0) = Phi;
    for (size_t kb = 0; kb < NP; ++kb) {
        for (size_t j = 0; j < NY; ++j) {
            const T sy = art.scale_y(j);
            art.b_bound(r_y_hi + (kb * NY) + j) = detail::scaled_bound(constraints.y_max(j), sy);
            art.b_bound(r_y_lo + (kb * NY) + j) = detail::scaled_bound(-constraints.y_min(j), sy);
            art.A_con(r_y_hi + (kb * NY) + j, NZ) = -constraints.y_max_ecr(j);
            art.A_con(r_y_lo + (kb * NY) + j, NZ) = -constraints.y_min_ecr(j);
        }
    }
    // Slack non-negativity: −ε ≤ 0.
    art.A_con(NI, NZ) = T{-1};
    art.b_bound(NI) = T{0};

    // The runtime QP requires H ≻ 0; validate once here.
    if (!mat::cholesky(art.H)) {
        return art;
    }

    art.success = true;
    return art;
}

} // namespace detail

/**
 * @brief Synthesize a constrained linear MPC (condensed dense QP, Δu form)
 *
 * Builds the prediction matrices for the augmented system
 * @f[
 *   \tilde A = \begin{bmatrix} A & B' \\ 0 & I \end{bmatrix},\quad
 *   \tilde B = \begin{bmatrix} B' \\ I \end{bmatrix},\quad
 *   \tilde C = [\,C' \;\; 0\,]
 * @f]
 * (primes: internally scaled B/C), condenses the horizon cost into
 * @f$ \tfrac{1}{2} \zeta^\top H \zeta + f^\top \zeta @f$ over
 * @f$ \zeta = [Z; \varepsilon] @f$ with
 * @f$ f = [F_x \tilde x_0 - F_r r' - F_{ut} u'_{target};\; 0] @f$, and stacks
 * the hard u/Δu and soft y box constraints as
 * @f$ A_{con} \zeta \le b_{bound} + B_{map}\tilde x_0 @f$.
 * Fails (success = false) if the plant is not discrete, has feedthrough
 * (D ≠ 0), the slack penalty is non-positive, or the resulting Hessian is not
 * positive definite.
 *
 * Beyond NC−1 the input is held (Δu = 0), which the augmented model carries
 * automatically.
 *
 * @note Compare with MATLAB®'s mpc(plant, Ts, p, m) — p = NP, m = NC.
 * @see MPCWeights, MPCConstraints, damp::MPC
 * @see Maciejowski, "Predictive Control with Constraints," 2002, ch. 2–3
 *
 * @tparam NP  Prediction horizon (steps)
 * @tparam NC  Control horizon (moves, NC ≤ NP; default NP)
 * @param sys                Discrete-time plant (Ts > 0, D = 0)
 * @param weights            Cost weights (default: Qy = I, RΔu = I)
 * @param constraints        Box limits (default: unbounded)
 * @param max_qp_iterations  Per-tick QP budget (0 → 10·(NZ+1 + NI+1))
 * @return MPCArtifacts for damp::MPC
 */
template<size_t NP, size_t NC = NP, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr MPCArtifacts<NX, NU, NY, NP, NC, 0, T> state_mpc(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const MPCWeights<NU, NY, T>&             weights = {},
    const MPCConstraints<NU, NY, T>&         constraints = {},
    size_t                                   max_qp_iterations = 0
) {
    return detail::mpc_impl<NP, NC, 0>(sys, Matrix<NX, 0, T>{}, weights, constraints, max_qp_iterations);
}

/**
 * @brief Synthesize a constrained linear MPC with measured-disturbance inputs
 *
 * Same as the plain overload, with the plant extended to
 * @f$ x_{k+1} = A x_k + B u_k + B_d d_k @f$ where @f$ d @f$ is a measured
 * disturbance (known feedforward input — load torque, feed rate, ambient
 * temperature) supplied to the runtime each tick via
 * @ref damp::MPC::control(const ColVec<NY,T>&, const ColVec<NX,T>&, const ColVec<ND,T>&).
 * The disturbance is held constant over the prediction horizon (the same
 * simplification as the reference; preview is a planned follow-on) and is not
 * rescaled.
 *
 * @note Compare with MATLAB®'s measured-disturbance channels (setmpcsignals 'MD').
 *
 * @tparam NP  Prediction horizon (steps)
 * @tparam NC  Control horizon (moves, NC ≤ NP; default NP)
 * @param sys  Discrete-time plant for the manipulated inputs (Ts > 0, D = 0)
 * @param Bd   Measured-disturbance input matrix (NX × ND)
 * @param weights            Cost weights (default: Qy = I, RΔu = I)
 * @param constraints        Box limits (default: unbounded)
 * @param max_qp_iterations  Per-tick QP budget (0 → 10·(NZ+1 + NI+1))
 * @return MPCArtifacts for damp::MPC
 */
template<size_t NP, size_t NC = NP, size_t NX, size_t NU, size_t NY, size_t ND, size_t NW, size_t NV, typename T>
[[nodiscard]] constexpr MPCArtifacts<NX, NU, NY, NP, NC, ND, T> state_mpc(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NX, ND, T>&                 Bd,
    const MPCWeights<NU, NY, T>&             weights = {},
    const MPCConstraints<NU, NY, T>&         constraints = {},
    size_t                                   max_qp_iterations = 0
) {
    return detail::mpc_impl<NP, NC, ND>(sys, Bd, weights, constraints, max_qp_iterations);
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
 * constexpr auto art = design::state_mpc<h.prediction_horizon, h.control_horizon>(sys);
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
    auto    np = static_cast<size_t>(steps);
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
 * @f$ k \cdot \tau @f$ (default k = 4, the ~2% settling criterion).
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

    const auto  eig = mat::compute_eigenvalues(sys.A);
    constexpr T stable_edge = T{1} - static_cast<T>(1e-9); // |λ| beyond this has no finite settling time
    constexpr T origin_tol = static_cast<T>(1e-12);        // deadbeat poles contribute no timescale

    T tau_max = T{0};
    for (size_t i = 0; i < NX; ++i) {
        const T mag = damp::abs(eig.values(i));
        if (mag <= origin_tol || mag >= stable_edge) {
            continue;
        }
        const T tau = -sys.Ts / damp::log(mag);
        tau_max = (tau > tau_max) ? tau : tau_max;
    }
    if (!(tau_max > T{0})) {
        return s; // no strictly stable pole — supply a settling-time target instead
    }

    s = suggest_mpc_horizon(settle_time_constants * tau_max, sys.Ts);
    s.time_constant = tau_max;
    return s;
}

/**
 * @struct MPCAnalysisModels
 * @brief LTI models of the unconstrained MPC loop, for margin/robustness analysis
 *
 * Inside the constraint-inactive region the receding-horizon law is the static
 * gain @f$ \Delta u_0 = -K [x; u_{prev}] + K_r\, r @f$ (physical units), so the
 * closed loop is an ordinary LTI system — run it through the existing
 * frequency-domain tooling (analysis::bode, matlab::margin/bandwidth,
 * analysis::poles) to evaluate a candidate weight set, exactly as the LQG/LQI
 * analysis models are used. This mirrors design::build_lqi_analysis_models.
 *
 * Constraints are ignored by construction: margins read from these models
 * describe the loop while no constraint is active.
 */
template<size_t NX, size_t NU, size_t NY, typename T = double>
struct MPCAnalysisModels {
    Matrix<NU, NX + NU, T>         K{};           ///< Unconstrained feedback gain on [x; u_prev]
    Matrix<NU, NY, T>              Kr{};          ///< Unconstrained reference gain
    StateSpace<NX + NU, NY, NY, T> closed_loop{}; ///< r → y closed loop over [x; u_prev]
    bool                           success{false};
};

/**
 * @brief Build the unconstrained-MPC LTI analysis models
 *
 * Solves the unconstrained condensed QP symbolically
 * (@f$ Z^* = -H^{-1}(F_x \tilde x - F_r r') @f$), extracts the first-move gain,
 * converts it back to physical units, and closes the loop over the augmented
 * state @f$ [x; u_{prev}] @f$. Measured-disturbance channels are held at zero.
 *
 * @note The MPC analog of MATLAB®'s ss(mpcobj) — valid while constraints are
 *       inactive.
 * @see MPCAnalysisModels, design::state_mpc
 *
 * @param sys        The plant the controller was synthesized for
 * @param artifacts  Result of design::state_mpc
 * @return Gains + closed-loop StateSpace (success = false if H is singular)
 */
template<size_t NP, size_t NC, size_t NX, size_t NU, size_t NY, size_t ND, size_t NW, size_t NV, typename T>
[[nodiscard]] constexpr MPCAnalysisModels<NX, NU, NY, T> build_mpc_analysis_models(
    const StateSpace<NX, NU, NY, T, NW, NV>&       sys,
    const MPCArtifacts<NX, NU, NY, NP, NC, ND, T>& artifacts
) {
    constexpr size_t NZ = MPCArtifacts<NX, NU, NY, NP, NC, ND, T>::NZ;

    if (!artifacts.success) {
        return MPCAnalysisModels<NX, NU, NY, T>{};
    }

    // Unconstrained first move in the scaled frame: Δu'0 = −[H⁻¹F_x]₀ x̃ + [H⁻¹F_r]₀ r'.
    const Matrix<NZ, NZ, T> Hz = artifacts.H.template block<NZ, NZ>(0, 0).to_matrix();
    const auto              Gx_opt = mat::solve(Hz, artifacts.F_x);
    const auto              Gr_opt = mat::solve(Hz, artifacts.F_r);
    if (!Gx_opt || !Gr_opt) {
        return MPCAnalysisModels<NX, NU, NY, T>{};
    }

    // Back to physical units: Δu = Su·Δu'0, x̃ = [x; Su⁻¹u_prev; d], r' = Sy⁻¹r.
    Matrix<NU, NX + NU, T> K{};
    Matrix<NU, NY, T>      Kr{};
    for (size_t i = 0; i < NU; ++i) {
        const T su = artifacts.scale_u(i);
        for (size_t j = 0; j < NX; ++j) {
            K(i, j) = su * Gx_opt.value()(i, j);
        }
        for (size_t j = 0; j < NU; ++j) {
            K(i, NX + j) = (su * Gx_opt.value()(i, NX + j)) / artifacts.scale_u(j);
        }
        for (size_t j = 0; j < NY; ++j) {
            Kr(i, j) = (su * Gr_opt.value()(i, j)) / artifacts.scale_y(j);
        }
    }

    // Closed loop over [x; u]: u⁺ = u + Δu, Δu = −K[x;u] + Kr·r.
    Matrix<NX + NU, NX + NU, T> A_open{};
    A_open.template block<NX, NX>(0, 0) = sys.A;
    A_open.template block<NX, NU>(0, NX) = sys.B;
    A_open.template block<NU, NU>(NX, NX) = Matrix<NU, NU, T>::identity();

    Matrix<NX + NU, NU, T> B_du{};
    B_du.template block<NX, NU>(0, 0) = sys.B;
    B_du.template block<NU, NU>(NX, 0) = Matrix<NU, NU, T>::identity();

    Matrix<NY, NX + NU, T> C_cl{};
    C_cl.template block<NY, NX>(0, 0) = sys.C;

    return MPCAnalysisModels<NX, NU, NY, T>{
        .K = K,
        .Kr = Kr,
        .closed_loop = StateSpace<NX + NU, NY, NY, T>{
            .A = A_open - (B_du * K),
            .B = B_du * Kr,
            .C = C_cl,
            .Ts = sys.Ts,
        },
        .success = true,
    };
}

} // namespace design

/**
 * @brief Runtime constrained MPC controller (fixed per-tick iteration budget)
 *
 * Each tick solves the condensed QP from design::state_mpc() at the
 * measured state and applies the first move: u = u_prev + Δu₀. The applied
 * command is always clamped to the u and Δu boxes, so even a budget-exhausted
 * or failed solve never emits an out-of-range actuator command (the QP's
 * intermediate iterates are dual, not primal, feasible). The soft-constraint
 * slack actually used is reported by last_slack().
 *
 * Full state feedback; pair with an observer/Kalman filter when x is not
 * measured. @p r is the output reference (NY), not a state x_ref — the cost
 * is ‖y − r‖²_Qy plus move / input-target terms, with no ‖x − x_ref‖ term.
 *
 * @note Compare with MATLAB®'s mpcmove(mpcobj, xc, ym, r, v).
 * @see design::state_mpc
 */
template<size_t NX, size_t NU, size_t NY, size_t NP, size_t NC, size_t ND = 0, typename T = float, typename Solver = design::WarmStartActiveSetSolver<(NC * NU) + 1, (4 * NC * NU) + (2 * NP * NY) + 1, T>>
    requires std::is_floating_point_v<T>
class MPC {
public:
    using Artifacts = design::MPCArtifacts<NX, NU, NY, NP, NC, ND, T>;
    static constexpr size_t NXA = Artifacts::NXA;
    static constexpr size_t NZ = Artifacts::NZ;
    static constexpr size_t NZ1 = Artifacts::NZ1;
    static constexpr size_t NI1 = Artifacts::NI1;

    constexpr MPC() = default;
    constexpr explicit MPC(const Artifacts& artifacts) : art(artifacts), u_target_(artifacts.u_target) {}
    /// Construct with a specific solver instance (for stateful solver policies)
    constexpr MPC(const Artifacts& artifacts, const Solver& solver)
        : art(artifacts), qp_solver_(solver), u_target_(artifacts.u_target) {}

    /**
     * @brief Compute the control move for the current tick (no measured disturbance)
     *
     * @p r is the output tracking reference (held constant over the horizon), not a
     * state reference. The condensed QP has no ‖x − x_ref‖ term, so this is not
     * u = −K(x − x_ref) and there is no LQR-style control(x_ref, x) overload
     * (that signature would collide with this one whenever NY = NX). To weight
     * or track a state that is not an output, add an extra C row and a matching
     * channel of r / Qy (see the file comment).
     *
     * @param r  Output reference r (NY), held over the horizon; not a state x_ref
     * @param x  Measured/estimated plant state
     * @return Input command u (clamped to the configured u and Δu boxes)
     */
    constexpr ColVec<NU, T> control(const ColVec<NY, T>& r, const ColVec<NX, T>& x) {
        return control_impl(r, x, ColVec<ND, T>{});
    }

    /**
     * @brief Compute the control move with a measured disturbance
     *
     * Same reference meaning as control(const ColVec<NY,T>&, const ColVec<NX,T>&):
     * @p r is the NY output reference, not x_ref.
     *
     * @param r  Output reference r (NY), held over the horizon; not a state x_ref
     * @param x  Measured/estimated plant state
     * @param d  Measured disturbance (held constant over the horizon)
     * @return Input command u (clamped to the configured u and Δu boxes)
     */
    constexpr ColVec<NU, T> control(const ColVec<NY, T>& r, const ColVec<NX, T>& x, const ColVec<ND, T>& d)
        requires(ND > 0)
    {
        return control_impl(r, x, d);
    }

    /// Reset controller memory (held input returns to zero)
    constexpr void reset() {
        u_prev = ColVec<NU, T>{};
        last_status_ = design::QPStatus::Success;
        last_iterations_ = 0;
        last_slack_ = T{0};
    }

    /// Seed the held input, e.g. for bumpless takeover from another controller
    constexpr void set_previous_control(const ColVec<NU, T>& u) { u_prev = u; }

    /// Override the preferred input operating point (physical units)
    constexpr void set_input_target(const ColVec<NU, T>& u_target) { u_target_ = u_target; }

    /// Access the solver policy (stateful solvers expose warm-start caches here)
    [[nodiscard]] constexpr Solver&       qp_solver() { return qp_solver_; }
    [[nodiscard]] constexpr const Solver& qp_solver() const { return qp_solver_; }

    [[nodiscard]] constexpr const ColVec<NU, T>& previous_control() const { return u_prev; }
    [[nodiscard]] constexpr const ColVec<NU, T>& input_target() const { return u_target_; }
    [[nodiscard]] constexpr design::QPStatus     last_status() const { return last_status_; }
    [[nodiscard]] constexpr size_t               last_iterations() const { return last_iterations_; }
    /// Soft-constraint slack ε from the last solve (0 = no output bound was relaxed)
    [[nodiscard]] constexpr T last_slack() const { return last_slack_; }

private:
    constexpr ColVec<NU, T> control_impl(const ColVec<NY, T>& r, const ColVec<NX, T>& x, const ColVec<ND, T>& d) {
        constexpr T sentinel = design::unbounded_bound<T>();
        constexpr T skip_threshold = sentinel / T{4};

        ColVec<NXA, T> xa{};
        for (size_t i = 0; i < NX; ++i) {
            xa(i) = x(i);
        }
        for (size_t j = 0; j < NU; ++j) {
            xa(NX + j) = u_prev(j) / art.scale_u(j);
        }
        for (size_t j = 0; j < ND; ++j) {
            xa(NX + NU + j) = d(j);
        }

        ColVec<NY, T> r_scaled{};
        for (size_t i = 0; i < NY; ++i) {
            r_scaled(i) = r(i) / art.scale_y(i);
        }
        ColVec<NU, T> ut_scaled{};
        for (size_t j = 0; j < NU; ++j) {
            ut_scaled(j) = u_target_(j) / art.scale_u(j);
        }

        const ColVec<NZ, T> f_z = art.F_x * xa - art.F_r * r_scaled - art.F_ut * ut_scaled;
        ColVec<NZ1, T>      f{};
        for (size_t i = 0; i < NZ; ++i) {
            f(i) = f_z(i);
        }

        const auto     shift = art.B_map * xa;
        ColVec<NI1, T> b{};
        for (size_t i = 0; i < NI1; ++i) {
            b(i) = (art.b_bound(i) >= skip_threshold) ? sentinel : art.b_bound(i) + shift(i);
        }

        // The default WarmStartActiveSetSolver caches chol(H) across ticks and
        // seeds each solve with the previous active set; the artifacts' QP
        // matrices are constant by construction, so its cache contract holds.
        const auto qp = qp_solver_(art.H, f, art.A_con, b, art.max_qp_iterations);
        last_status_ = qp.status;
        last_iterations_ = qp.iterations;
        last_slack_ = qp.x(NZ);

        ColVec<NU, T> du{};
        if (qp.status == design::QPStatus::Success || qp.status == design::QPStatus::MaxIterations) {
            for (size_t j = 0; j < NU; ++j) {
                du(j) = qp.x(j) * art.scale_u(j); // back to physical units
            }
        }

        // Safety clamp: guarantees box feasibility of the applied command.
        ColVec<NU, T> u{};
        for (size_t j = 0; j < NU; ++j) {
            du(j) = damp::clamp(du(j), art.constraints.du_min(j), art.constraints.du_max(j));
            u(j) = damp::clamp(u_prev(j) + du(j), art.constraints.u_min(j), art.constraints.u_max(j));
        }
        u_prev = u;
        return u;
    }

    Artifacts        art{};
    Solver           qp_solver_{};
    ColVec<NU, T>    u_prev{};
    ColVec<NU, T>    u_target_{};
    design::QPStatus last_status_{design::QPStatus::Success};
    size_t           last_iterations_{0};
    T                last_slack_{0};
};

/// Deduce the runtime from its artifacts: MPC controller{art};
template<size_t NX, size_t NU, size_t NY, size_t NP, size_t NC, size_t ND, typename T>
MPC(const design::MPCArtifacts<NX, NU, NY, NP, NC, ND, T>&) -> MPC<NX, NU, NY, NP, NC, ND, T>;

/// Deduce the runtime from artifacts + a solver instance: MPC controller{art, solver};
template<size_t NX, size_t NU, size_t NY, size_t NP, size_t NC, size_t ND, typename T, typename Solver>
MPC(const design::MPCArtifacts<NX, NU, NY, NP, NC, ND, T>&, const Solver&) -> MPC<NX, NU, NY, NP, NC, ND, T, Solver>;

} // namespace damp
