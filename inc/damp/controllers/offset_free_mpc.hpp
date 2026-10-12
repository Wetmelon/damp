// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file offset_free_mpc.hpp
 * @brief Offset-free constrained MPC: the damp::MPC controller bundled with a
 *        disturbance-augmented Kalman estimator (`design::mpc`
 *        + `damp::OffsetFreeMPC`).
 *
 * Closes the two usability gaps of the bare MPC runtime in one feature, the
 * same way LQG bundles LQR with a Kalman filter (controllers/lqg.hpp):
 *
 * 1. Output feedback — the plant state no longer needs to be measured; the
 *    estimator reconstructs it from the measured outputs.
 * 2. Offset-free tracking — the plant model is augmented with an
 *    integrating disturbance @f$ d @f$:
 *    @f[
 *      x_{k+1} = A x_k + B u_k + B_d d_k, \qquad d_{k+1} = d_k, \qquad
 *      y_k = C x_k + C_d d_k + D u_k,
 *    @f]
 *    and the Kalman filter estimates @f$ [\hat x;\, \hat d] @f$. The estimate
 *    @f$ \hat d @f$ is fed into the MPC's measured-disturbance channel
 *    (@f$ B_d @f$), so predictions are disturbance-corrected and the closed loop
 *    holds the controlled outputs @f$ z = C_z x + D_z u @f$ at @f$ r @f$ with no
 *    steady-state error under constant disturbances and plant/model mismatch —
 *    the cases the bare velocity form cannot reject.
 *
 * The disturbance model is the one design::discrete_oflqg takes, OffsetFreeModel:
 * - where the disturbance enters (@f$ B_d @f$): at the input, or where it physically
 *   acts, e.g. a force on a load behind a compliant transmission;
 * - which sensors feel it directly (@f$ C_d @f$), e.g. an accelerometer on that load;
 * - which outputs the MPC tracks and constrains (@f$ C_z, D_z @f$), separately from the
 *   measured outputs @f$ C @f$: estimate with every sensor, hold only the ones that matter.
 *
 * Without a model, design::mpc uses MATLAB®'s input-disturbance (`setindist`)
 * model with one integrator per manipulated input (@f$ B_d = B,\ C_d = 0 @f$)
 * and tracks every measured output (@f$ C_z = C @f$). It matches the physics of
 * servo/drive load disturbances and stays detectable on integrating plants — an
 * output-integrating disturbance on a plant with poles at z = 1 (e.g. a double
 * integrator) is undetectable, which is why MATLAB® silently falls back to input
 * models there.
 *
 * Detectability (Pannocchia & Rawlings). The augmented model is detectable iff
 * @f$ (A, C) @f$ is and
 * @f$ \operatorname{rank}\begin{bmatrix} A - I & B_d \\ C & C_d \end{bmatrix} = N_X + N_D @f$,
 * which needs at least as many sensors as disturbances: ND ≤ NY is enforced at
 * compile time, the rank at synthesis time, and an undetectable @f$ (A, C) @f$
 * fails the filter's DARE (success = false either way). Offset-free tracking of
 * every controlled output further needs NZ ≤ NU and a steady state that reaches
 * any r; with more controlled outputs than inputs the MPC settles on its
 * weighted least-squares compromise.
 *
 * Delays: design on the model from discretize(plant, Ts, InputDelay<D>{…}) and/or
 * OutputDelay<D>{}, with the OffsetFreeModel and Q_process still sized to the plant:
 * they are the leading block, and the delay states get no disturbance, no process
 * noise and no controlled output.
 *
 * @note Compare with MATLAB®'s mpc default estimator: getindist/setindist +
 *       getEstimator/setEstimator + mpcstate; mpcmove(mpcobj, xc, ym, r).
 * @see controllers/mpc.hpp for the underlying controller,
 *      estimation/kalman.hpp for the estimator, controllers/lqg.hpp for the
 *      bundling pattern this mirrors, controllers/offset_free_lqg.hpp for the
 *      unconstrained infinite-horizon form
 * @see Pannocchia & Rawlings, "Disturbance Models for Offset-Free
 *      Model-Predictive Control," AIChE Journal 49(2), 2003,
 *      https://doi.org/10.1002/aic.690490213
 * @see Muske & Badgwell, "Disturbance modeling for offset-free linear model
 *      predictive control," Journal of Process Control 12(5), 2002,
 *      https://doi.org/10.1016/S0959-1524(01)00051-8
 *
 * @code
 * constexpr auto art = design::mpc<15, 5>(sys, weights, limits); // input disturbance, z = y
 * static_assert(art.success);
 * OffsetFreeMPC controller{art.as<float>()};
 * // each tick, from the output measurement alone:
 * //   u = controller.control(r, y);
 *
 * // A load force sensed by an accelerometer; hold the load position only.
 * design::OffsetFreeModel<4, 1, 2, 1, 1> load{.Bd = G_load, .Cd = accel_coupling, .Cz = load_position};
 * constexpr auto art2 = design::mpc<40, 8>(sys, load, weights_z, limits_z);
 * @endcode
 */

#include <cstddef>

#include "damp/controllers/mpc.hpp"
#include "damp/controllers/offset_free_lqg.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/block.hpp"
#include "damp/matrix/functions.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace design {

/**
 * @struct OffsetFreeMPCArtifacts
 * @brief Combined MPC + disturbance-augmented Kalman design, consumed by damp::OffsetFreeMPC
 * @tparam ND Disturbance count (default: one per input); NZ Controlled outputs (default: the measured ones)
 */
template<size_t NX, size_t NU, size_t NY, size_t NP, size_t NC, typename T = double, size_t ND = NU, size_t NZ = NY>
struct OffsetFreeMPCArtifacts {
    static constexpr size_t NXK = NX + ND; ///< Estimator state size [x; d]

    MPCArtifacts<NX, NU, NZ, NP, NC, ND, T> mpc{};    ///< MPC on the controlled outputs, disturbance channel B_d
    KalmanResult<NXK, NU, NY, T, NXK, NY>   kalman{}; ///< Steady-state filter on the augmented model
    bool                                    success{false};

    template<typename U>
    [[nodiscard]] constexpr OffsetFreeMPCArtifacts<NX, NU, NY, NP, NC, U, ND, NZ> as() const {
        return OffsetFreeMPCArtifacts<NX, NU, NY, NP, NC, U, ND, NZ>{
            mpc.template as<U>(),
            kalman.template as<U>(),
            success,
        };
    }
};

/**
 * @brief Synthesize an offset-free constrained MPC with a chosen disturbance model
 *
 * Runs design::state_mpc on the controlled outputs (C_z, D_z) with the
 * measured-disturbance channel wired to B_d, and design::kalman on the
 * disturbance-augmented model
 * @f[
 *   \tilde A = \begin{bmatrix} A & B_d \\ 0 & I \end{bmatrix},\quad
 *   \tilde B = \begin{bmatrix} B \\ 0 \end{bmatrix},\quad
 *   \tilde C = [\,C \;\; C_d\,],\quad \tilde D = D.
 * @f]
 * Fails (success = false) if the Pannocchia–Rawlings rank condition fails
 * (rank [A−I, B_d; C, C_d] < NX + ND), the MPC synthesis fails (which includes
 * D_z ≠ 0), or the filter DARE does not converge (which includes an undetectable
 * (A, C)). Weights, constraints and the reference act on z, not y.
 *
 * Q_disturbance sets how fast the disturbance estimate moves: larger values
 * relative to R_measurement give faster rejection of load steps at the cost of
 * more measurement noise feeding through to the input.
 *
 * @tparam NP  Prediction horizon (steps)
 * @tparam NC  Control horizon (moves, NC ≤ NP; default NP)
 * @param sys             Discrete-time plant (Ts > 0); C, D are the measured outputs
 * @param model           Disturbance entry Bd, its direct sensor coupling Cd, controlled outputs Cz, Dz (Dz = 0)
 * @param weights         MPC cost weights on z (default: Qy = I, RΔu = I)
 * @param constraints     MPC box limits on u, Δu and z (default: unbounded)
 * @param Q_process       Process-noise covariance on x (default I); plant-sized like @p model on a delayed model
 * @param R_measurement   Measurement-noise covariance (NY × NY, default I)
 * @param Q_disturbance   Random-walk covariance on d (ND × ND, default I)
 * @param max_qp_iterations  Per-tick QP budget (0 → default)
 * @return OffsetFreeMPCArtifacts for damp::OffsetFreeMPC
 */
template<
    size_t NP, size_t NC = NP, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NXM, size_t ND,
    size_t NZ, size_t NXQ = NXM>
[[nodiscard]] constexpr OffsetFreeMPCArtifacts<NX, NU, NY, NP, NC, T, ND, NZ> mpc(
    const StateSpace<NX, NU, NY, T, NW, NV>&       sys,
    const OffsetFreeModel<NXM, NU, NY, ND, NZ, T>& model,
    const MPCWeights<NU, NZ, T>&                   weights = {},
    const MPCConstraints<NU, NZ, T>&               constraints = {},
    const Matrix<NXQ, NXQ, T>&                     Q_process = Matrix<NXQ, NXQ, T>::identity(),
    const Matrix<NY, NY, T>&                       R_measurement = Matrix<NY, NY, T>::identity(),
    const Matrix<ND, ND, T>&                       Q_disturbance = Matrix<ND, ND, T>::identity(),
    size_t                                         max_qp_iterations = 0
) {
    static_assert(ND >= 1, "An offset-free model estimates at least one disturbance.");
    static_assert(ND <= NY, "Estimating ND disturbances needs at least ND independent measurements (ND <= NY).");
    static_assert(NXM <= NX && NXQ <= NX, "The model and Q_process cover the plant: the leading states of sys.");

    constexpr size_t NXK = NX + ND;

    OffsetFreeMPCArtifacts<NX, NU, NY, NP, NC, T, ND, NZ> art{};

    // Delay states, if any, follow the plant's: no disturbance, no controlled output.
    const auto Bd = damp::detail::pad_leading<NX, ND>(model.Bd);
    const auto Cz = damp::detail::pad_leading<NZ, NX>(model.Cz);

    // Pannocchia–Rawlings detectability condition for the augmented model:
    // rank [A−I, Bd; C, Cd] must be NX + ND.
    Matrix<NX + NY, NX + ND, T> rank_test{};
    rank_test.template block<NX, NX>(0, 0) = sys.A - Matrix<NX, NX, T>::identity();
    rank_test.template block<NX, ND>(0, NX) = Bd;
    rank_test.template block<NY, NX>(NX, 0) = sys.C;
    rank_test.template block<NY, ND>(NX, NX) = model.Cd;
    if (mat::rank(rank_test) != NX + ND) {
        return art;
    }

    // Controller: predicts and constrains z; the estimated disturbance rides the
    // measured-disturbance channel.
    const StateSpace<NX, NU, NZ, T> sys_z{.A = sys.A, .B = sys.B, .C = Cz, .D = model.Dz, .Ts = sys.Ts};
    art.mpc = detail::mpc_impl<NP, NC, ND>(sys_z, Bd, weights, constraints, max_qp_iterations);

    // Estimator: steady-state Kalman filter on the augmented model (G = H = I).
    StateSpace<NXK, NU, NY, T, NXK, NY> sys_aug{};
    sys_aug.A.template block<NX, NX>(0, 0) = sys.A;
    sys_aug.A.template block<NX, ND>(0, NX) = Bd;
    sys_aug.A.template block<ND, ND>(NX, NX) = Matrix<ND, ND, T>::identity();
    sys_aug.B.template block<NX, NU>(0, 0) = sys.B;
    sys_aug.C.template block<NY, NX>(0, 0) = sys.C;
    sys_aug.C.template block<NY, ND>(0, NX) = model.Cd;
    sys_aug.D = sys.D;
    sys_aug.Ts = (sys.Ts > T{0}) ? sys.Ts : T{1}; // discrete DARE requires Ts > 0

    Matrix<NXK, NXK, T> Q_aug{};
    Q_aug.template block<NXQ, NXQ>(0, 0) = Q_process; // delay states, if any: no process noise
    Q_aug.template block<ND, ND>(NX, NX) = Q_disturbance;

    art.kalman = kalman(sys_aug, Q_aug, R_measurement);

    art.success = art.mpc.success && art.kalman.success;
    return art;
}

/**
 * @brief Synthesize an offset-free constrained MPC, input disturbance, every measured output tracked
 *
 * The default form of the model overload: one integrating disturbance per input,
 * entering with it (B_d = B) but not the sensors (C_d = 0), and z = y (C_z = C,
 * D_z = D). Plants with feedthrough (D ≠ 0) fail; NU ≤ NY is enforced at compile time.
 *
 * @tparam NP  Prediction horizon (steps)
 * @tparam NC  Control horizon (moves, NC ≤ NP; default NP)
 * @param sys             Discrete-time plant (Ts > 0, D = 0)
 * @param weights         MPC cost weights (default: Qy = I, RΔu = I)
 * @param constraints     MPC box limits (default: unbounded)
 * @param Q_process       State process-noise covariance (default I). On a model with delay
 *                        states (discretize with InputDelay / OutputDelay) give it plant-sized:
 *                        the delay states that follow get none.
 * @param R_measurement   Measurement-noise covariance (default I)
 * @param Q_disturbance   Disturbance drift covariance (default I)
 * @param max_qp_iterations  Per-tick QP budget (0 → default)
 * @return OffsetFreeMPCArtifacts for damp::OffsetFreeMPC
 */
template<size_t NP, size_t NC = NP, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NXQ = NX>
[[nodiscard]] constexpr OffsetFreeMPCArtifacts<NX, NU, NY, NP, NC, T> mpc(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const MPCWeights<NU, NY, T>&             weights = {},
    const MPCConstraints<NU, NY, T>&         constraints = {},
    const Matrix<NXQ, NXQ, T>&               Q_process = Matrix<NXQ, NXQ, T>::identity(),
    const Matrix<NY, NY, T>&                 R_measurement = Matrix<NY, NY, T>::identity(),
    const Matrix<NU, NU, T>&                 Q_disturbance = Matrix<NU, NU, T>::identity(),
    size_t                                   max_qp_iterations = 0
) {
    static_assert(NU <= NY, "Estimating NU input disturbances needs at least NU independent measurements (NU <= NY).");
    const OffsetFreeModel<NX, NU, NY, NU, NY, T> model{.Bd = sys.B, .Cz = sys.C, .Dz = sys.D};
    return mpc<NP, NC>(sys, model, weights, constraints, Q_process, R_measurement, Q_disturbance, max_qp_iterations);
}

} // namespace design

/**
 * @brief Runtime offset-free MPC: constrained MPC + disturbance-augmented Kalman filter
 *
 * Mirrors the LQG bundling pattern: the controller and estimator are plain
 * members, `predict`/`update` forward to the filter for caller-sequenced use,
 * and `control(r, y)` runs the whole tick from the output measurement alone —
 * predict with the previously applied input, update with y, then solve the MPC
 * at the estimated state with the estimated disturbance on the feedforward
 * channel. @p r is the reference for the NZ controlled outputs (the measured
 * outputs in the default form); there is no x_ref overload (the estimator
 * reconstructs x, and control(x_ref) would collide with control(r) whenever
 * NX = NZ).
 *
 * @note Compare with MATLAB®'s mpcmove(mpcobj, xc, ym, r) with the default
 *       built-in estimator.
 * @see design::mpc
 */
template<size_t NX, size_t NU, size_t NY, size_t NP, size_t NC, typename T = float, size_t ND = NU, size_t NZ = NY>
struct OffsetFreeMPC {
    static constexpr bool   output_feedback_tick = true;
    static constexpr size_t NXK = NX + ND;
    using Artifacts = design::OffsetFreeMPCArtifacts<NX, NU, NY, NP, NC, T, ND, NZ>;

    MPC<NX, NU, NZ, NP, NC, ND, T>        mpc{};
    KalmanFilter<NXK, NU, NY, T, NXK, NY> kf{};

    constexpr OffsetFreeMPC() = default;
    constexpr explicit OffsetFreeMPC(const Artifacts& artifacts) : mpc(artifacts.mpc), kf(artifacts.kalman) {}

    /// Advance the estimator with the applied input (caller-sequenced use)
    constexpr void predict(const ColVec<NU, T>& u) { kf.predict(u); }

    /// Correct the estimator with a measurement (caller-sequenced use)
    constexpr bool update(const ColVec<NY, T>& y, const ColVec<NU, T>& u = ColVec<NU, T>{}) { return kf.update(y, u); }

    /// MPC move at the current estimate (call predict/update first). @p r is the NZ controlled-output reference, not x_ref.
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NZ, T>& r) {
        return mpc.control(r, state_estimate(), disturbance_estimate());
    }

    /**
     * @brief Full tick from the output measurement alone
     *
     * Runs predict (with the previously applied input), update (with y), then
     * the MPC solve. Equivalent to the caller-sequenced predict/update/control
     * cycle.
     *
     * @param r  Controlled-output reference (NZ); not a state x_ref
     * @param y  Measured plant output (NY)
     * @return Input command u (clamped to the configured boxes)
     */
    constexpr ColVec<NU, T> control(const ColVec<NZ, T>& r, const ColVec<NY, T>& y) {
        const ColVec<NU, T> u = feedback(r, y);
        commit(u);
        return u;
    }

    /**
     * @brief First half of a saturation-aware tick: estimator update, then this tick's move
     *
     * Predict with the previously committed input, update with @p y, and solve the MPC at
     * the estimate. The move is not held until commit, so a downstream clamp or override
     * reaches both the next Δu and the estimator's next predict.
     *
     * @param r  Controlled-output reference (NZ); not a state x_ref
     * @param y  Measured plant output (NY)
     * @return Input command u (clamped to the configured boxes)
     */
    [[nodiscard]] constexpr ColVec<NU, T> feedback(const ColVec<NZ, T>& r, const ColVec<NY, T>& y) {
        kf.predict(mpc.previous_control());
        kf.update(y, mpc.previous_control());
        return mpc.feedback(r, state_estimate(), disturbance_estimate());
    }

    /// Second half of the tick: the input the plant actually received (post-clamp).
    constexpr void commit(const ColVec<NU, T>& u_applied) { mpc.commit(u_applied); }

    /// Estimated plant state x̂
    [[nodiscard]] constexpr ColVec<NX, T> state_estimate() const {
        ColVec<NX, T> x{};
        for (size_t i = 0; i < NX; ++i) {
            x(i) = kf.state()(i);
        }
        return x;
    }

    /// Estimated disturbance d̂ (what the integrating model has absorbed)
    [[nodiscard]] constexpr ColVec<ND, T> disturbance_estimate() const {
        ColVec<ND, T> d{};
        for (size_t j = 0; j < ND; ++j) {
            d(j) = kf.state()(NX + j);
        }
        return d;
    }

    /// Reset the controller memory; the state estimate is left in place
    constexpr void reset() { mpc.reset(); }
};

/// Deduce the runtime from its artifacts: OffsetFreeMPC controller{art};
template<size_t NX, size_t NU, size_t NY, size_t NP, size_t NC, typename T, size_t ND, size_t NZ>
OffsetFreeMPC(const design::OffsetFreeMPCArtifacts<NX, NU, NY, NP, NC, T, ND, NZ>&) -> OffsetFreeMPC<NX, NU, NY, NP, NC, T, ND, NZ>;

} // namespace damp
