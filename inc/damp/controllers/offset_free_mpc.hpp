// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file mpc.hpp
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
 *    *integrating input disturbance* @f$ w @f$:
 *    @f[
 *      x_{k+1} = A x_k + B (u_k + w_k), \qquad w_{k+1} = w_k, \qquad y_k = C x_k,
 *    @f]
 *    and the Kalman filter estimates @f$ [\hat x;\, \hat w] @f$. The estimate
 *    @f$ \hat w @f$ is fed into the MPC's measured-disturbance channel
 *    (@f$ B_d = B @f$), so predictions are disturbance-corrected and the
 *    closed loop settles with zero steady-state error under constant input
 *    disturbances and plant/model gain mismatch — the cases the bare velocity
 *    form cannot reject.
 *
 * Disturbance model choice. This is MATLAB®'s input-disturbance
 * (`setindist`) model with one integrator per manipulated input. The input
 * model is the first cut (rather than MATLAB®'s output-integrator default)
 * because it matches the physics of servo/drive load disturbances, reuses the
 * MPC's measured-disturbance channel unchanged, and stays detectable on
 * integrating plants — an output-integrating disturbance on a plant with poles
 * at z = 1 (e.g. a double integrator) is undetectable, which is why MATLAB®
 * silently falls back to input models there. The output-integrating model is a
 * planned follow-on. Estimating NU disturbances needs at least NU independent
 * measurements: NU ≤ NY is enforced at compile time, and the
 * Pannocchia–Rawlings rank condition
 * @f$ \operatorname{rank}\begin{bmatrix} A - I & B \\ C & 0 \end{bmatrix} = N_X + N_U @f$
 * is checked at synthesis time (success = false when it fails).
 *
 * @note Compare with MATLAB®'s mpc default estimator: getindist/setindist +
 *       getEstimator/setEstimator + mpcstate; mpcmove(mpcobj, xc, ym, r).
 * @see controllers/mpc.hpp for the underlying controller,
 *      estimation/kalman.hpp for the estimator, controllers/lqg.hpp for the
 *      bundling pattern this mirrors
 * @see Pannocchia & Rawlings, "Disturbance Models for Offset-Free
 *      Model-Predictive Control," AIChE Journal 49(2), 2003,
 *      https://doi.org/10.1002/aic.690490213
 * @see Muske & Badgwell, "Disturbance modeling for offset-free linear model
 *      predictive control," Journal of Process Control 12(5), 2002,
 *      https://doi.org/10.1016/S0959-1524(01)00051-8
 *
 * @code
 * constexpr auto art = design::mpc<15, 5>(sys, weights, limits);
 * static_assert(art.success);
 * OffsetFreeMPC controller{art.as<float>()};
 * // each tick, from the output measurement alone:
 * //   u = controller.control(r, y);
 * @endcode
 */

#include <cstddef>

#include "damp/controllers/mpc.hpp"
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
 */
template<size_t NX, size_t NU, size_t NY, size_t NP, size_t NC, typename T = double>
struct OffsetFreeMPCArtifacts {
    static constexpr size_t NXK = NX + NU; ///< Estimator state size [x; w]

    MPCArtifacts<NX, NU, NY, NP, NC, NU, T> mpc{};    ///< MPC with the disturbance channel B_d = B
    KalmanResult<NXK, NU, NY, T, NXK, NY>   kalman{}; ///< Steady-state filter on the augmented model
    bool                                    success{false};

    template<typename U>
    [[nodiscard]] constexpr OffsetFreeMPCArtifacts<NX, NU, NY, NP, NC, U> as() const {
        return OffsetFreeMPCArtifacts<NX, NU, NY, NP, NC, U>{
            mpc.template as<U>(),
            kalman.template as<U>(),
            success,
        };
    }
};

/**
 * @brief Synthesize an offset-free constrained MPC (controller + estimator)
 *
 * Runs design::state_mpc with the measured-disturbance channel wired to
 * the input matrix (B_d = B), and design::kalman on the input-disturbance-
 * augmented model
 * @f[
 *   \tilde A = \begin{bmatrix} A & B \\ 0 & I \end{bmatrix},\quad
 *   \tilde B = \begin{bmatrix} B \\ 0 \end{bmatrix},\quad
 *   \tilde C = [\,C \;\; 0\,].
 * @f]
 * Fails (success = false) if the underlying MPC synthesis fails, the
 * augmented-model detectability rank condition fails, or the filter DARE does
 * not converge.
 *
 * Q_disturbance sets how fast the disturbance estimate moves: larger values
 * relative to R_measurement give faster rejection of load steps at the cost of
 * more measurement noise feeding through to the input.
 *
 * @tparam NP  Prediction horizon (steps)
 * @tparam NC  Control horizon (moves, NC ≤ NP; default NP)
 * @param sys             Discrete-time plant (Ts > 0, D = 0)
 * @param weights         MPC cost weights (default: Qy = I, RΔu = I)
 * @param constraints     MPC box limits (default: unbounded)
 * @param Q_process       State process-noise covariance (default I)
 * @param Q_disturbance   Disturbance drift covariance (default I)
 * @param R_measurement   Measurement-noise covariance (default I)
 * @param max_qp_iterations  Per-tick QP budget (0 → default)
 * @return OffsetFreeMPCArtifacts for damp::OffsetFreeMPC
 */
template<size_t NP, size_t NC = NP, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr OffsetFreeMPCArtifacts<NX, NU, NY, NP, NC, T> mpc(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const MPCWeights<NU, NY, T>&             weights = {},
    const MPCConstraints<NU, NY, T>&         constraints = {},
    const Matrix<NX, NX, T>&                 Q_process = Matrix<NX, NX, T>::identity(),
    const Matrix<NU, NU, T>&                 Q_disturbance = Matrix<NU, NU, T>::identity(),
    const Matrix<NY, NY, T>&                 R_measurement = Matrix<NY, NY, T>::identity(),
    size_t                                   max_qp_iterations = 0
) {
    static_assert(NU <= NY, "Estimating NU input disturbances needs at least NU independent measurements (NU <= NY).");

    constexpr size_t NXK = NX + NU;

    OffsetFreeMPCArtifacts<NX, NU, NY, NP, NC, T> art{};

    // Pannocchia–Rawlings detectability condition for the augmented model:
    // rank [A−I, B; C, 0] must be NX + NU.
    Matrix<NX + NY, NX + NU, T> rank_test{};
    rank_test.template block<NX, NX>(0, 0) = sys.A - Matrix<NX, NX, T>::identity();
    rank_test.template block<NX, NU>(0, NX) = sys.B;
    rank_test.template block<NY, NX>(NX, 0) = sys.C;
    if (mat::rank(rank_test) != NX + NU) {
        return art;
    }

    // Controller: the estimated disturbance rides the measured-disturbance channel.
    art.mpc = detail::mpc_impl<NP, NC, NU>(sys, sys.B, weights, constraints, max_qp_iterations);

    // Estimator: steady-state Kalman filter on the augmented model (G = H = I).
    StateSpace<NXK, NU, NY, T, NXK, NY> sys_aug{};
    sys_aug.A.template block<NX, NX>(0, 0) = sys.A;
    sys_aug.A.template block<NX, NU>(0, NX) = sys.B;
    sys_aug.A.template block<NU, NU>(NX, NX) = Matrix<NU, NU, T>::identity();
    sys_aug.B.template block<NX, NU>(0, 0) = sys.B;
    sys_aug.C.template block<NY, NX>(0, 0) = sys.C;
    sys_aug.Ts = (sys.Ts > T{0}) ? sys.Ts : T{1}; // discrete DARE requires Ts > 0

    Matrix<NXK, NXK, T> Q_aug{};
    Q_aug.template block<NX, NX>(0, 0) = Q_process;
    Q_aug.template block<NU, NU>(NX, NX) = Q_disturbance;

    art.kalman = kalman(sys_aug, Q_aug, R_measurement);

    art.success = art.mpc.success && art.kalman.success;
    return art;
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
 * channel.
 *
 * @note Compare with MATLAB®'s mpcmove(mpcobj, xc, ym, r) with the default
 *       built-in estimator.
 * @see design::mpc
 */
template<size_t NX, size_t NU, size_t NY, size_t NP, size_t NC, typename T = float>
struct OffsetFreeMPC {
    static constexpr size_t NXK = NX + NU;
    using Artifacts = design::OffsetFreeMPCArtifacts<NX, NU, NY, NP, NC, T>;

    MPC<NX, NU, NY, NP, NC, NU, T>        mpc{};
    KalmanFilter<NXK, NU, NY, T, NXK, NY> kf{};

    constexpr OffsetFreeMPC() = default;
    constexpr explicit OffsetFreeMPC(const Artifacts& artifacts) : mpc(artifacts.mpc), kf(artifacts.kalman) {}

    /// Advance the estimator with the applied input (caller-sequenced use)
    constexpr void predict(const ColVec<NU, T>& u) { kf.predict(u); }

    /// Correct the estimator with a measurement (caller-sequenced use)
    constexpr bool update(const ColVec<NY, T>& y, const ColVec<NU, T>& u = ColVec<NU, T>{}) { return kf.update(y, u); }

    /// MPC move at the current estimate (call predict/update first)
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NY, T>& r) {
        return mpc.control(r, state_estimate(), disturbance_estimate());
    }

    /**
     * @brief Full tick from the output measurement alone
     *
     * Runs predict (with the previously applied input), update (with y), then
     * the MPC solve. Equivalent to the caller-sequenced predict/update/control
     * cycle.
     *
     * @param r  Output reference
     * @param y  Measured plant output
     * @return Input command u (clamped to the configured boxes)
     */
    constexpr ColVec<NU, T> control(const ColVec<NY, T>& r, const ColVec<NY, T>& y) {
        kf.predict(mpc.previous_control());
        kf.update(y, mpc.previous_control());
        return control(r);
    }

    /// Estimated plant state x̂
    [[nodiscard]] constexpr ColVec<NX, T> state_estimate() const {
        ColVec<NX, T> x{};
        for (size_t i = 0; i < NX; ++i) {
            x(i) = kf.state()(i);
        }
        return x;
    }

    /// Estimated input disturbance ŵ (what the integrating model has absorbed)
    [[nodiscard]] constexpr ColVec<NU, T> disturbance_estimate() const {
        ColVec<NU, T> w{};
        for (size_t j = 0; j < NU; ++j) {
            w(j) = kf.state()(NX + j);
        }
        return w;
    }

    /// Reset the controller memory; the state estimate is left in place
    constexpr void reset() { mpc.reset(); }
};

/// Deduce the runtime from its artifacts: OffsetFreeMPC controller{art};
template<size_t NX, size_t NU, size_t NY, size_t NP, size_t NC, typename T>
OffsetFreeMPC(const design::OffsetFreeMPCArtifacts<NX, NU, NY, NP, NC, T>&) -> OffsetFreeMPC<NX, NU, NY, NP, NC, T>;

} // namespace damp
