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
 * OutputDelay<D>{}, with the OffsetFreeModel and the process noise still sized to the plant:
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
 * static constinit OffsetFreeMPC controller{art.as<float>()}; // a failed design does not compile
 * // each tick, from the output measurement alone:
 * //   u = controller.control(r, y);
 *
 * // A load force sensed by an accelerometer; hold the load position only.
 * design::OffsetFreeModel<4, 1, 2, 1, 1> load{.Bd = G_load, .Cd = accel_coupling, .Cz = load_position};
 * constexpr auto art2 = design::mpc<40, 8>(sys, load, weights_z, limits_z);
 * @endcode
 */

#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/delay_predictor.hpp"
#include "damp/controllers/mpc.hpp"
#include "damp/controllers/offset_free_lqg.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/block.hpp"
#include "damp/matrix/functions.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/svd.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace design {

/**
 * @struct OffsetFreeMPCArtifacts
 * @brief Combined MPC + disturbance-augmented Kalman design, consumed by damp::OffsetFreeMPC
 * @tparam ND Disturbance count (default: one per input); NZ Controlled outputs (default: the measured ones)
 * @tparam D  Whole periods of input delay the runtime predicts (design::mpc<NP, NC, Periods>); the runtime deduces it from here
 */
template<size_t NX, size_t NU, size_t NY, size_t NP, size_t NC, typename T = double, size_t ND = NU, size_t NZ = NY, size_t D = 0>
struct OffsetFreeMPCArtifacts {
    static constexpr size_t NXK = NX + ND; ///< Estimator state size [x; d]

    MPCArtifacts<NX, NU, NZ, NP, NC, ND, T> mpc{};      ///< MPC on the controlled outputs, disturbance channel B_d
    KalmanResult<NXK, NU, NY, T, NXK, NY>   kalman{};   ///< Steady-state filter on the augmented model
    Matrix<NX + NU, ND, T>                  target_d{}; ///< [x_t; u_t] per unit estimated disturbance
    Matrix<NX + NU, NZ, T>                  target_r{}; ///< [x_t; u_t] per unit controlled-output reference
    bool                                    success{false};

    template<typename U>
    [[nodiscard]] constexpr OffsetFreeMPCArtifacts<NX, NU, NY, NP, NC, U, ND, NZ, D> as() const {
        return OffsetFreeMPCArtifacts<NX, NU, NY, NP, NC, U, ND, NZ, D>{
            mpc.template as<U>(),
            kalman.template as<U>(),
            target_d.template as<U>(),
            target_r.template as<U>(),
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
 * The disturbance noise sets how fast the disturbance estimate moves: larger values
 * relative to the measurement noise give faster rejection of load steps at the cost of
 * more measurement noise feeding through to the input.
 *
 * @tparam NP       Prediction horizon (steps)
 * @tparam NC       Control horizon (moves, NC ≤ NP)
 * @tparam Periods  Input delay the runtime predicts, in whole sampling periods, as a floating value
 *                  (default 0.0): a discrete model carries no fraction of a period, so pass the
 *                  continuous plant and Ts for one. The runtime deduces it from the result.
 * @param sys             Discrete-time plant (Ts > 0); C, D are the measured outputs
 * @param model           Disturbance entry Bd, its direct sensor coupling Cd, controlled outputs Cz, Dz (Dz = 0)
 * @param weights         MPC cost weights on z (default: Qy = I, RΔu = I)
 * @param constraints     MPC box limits on u, Δu and z (default: unbounded)
 * @param noise           Process (plant-sized like @p model on a delayed model), measurement and
 *                        disturbance covariances (default I)
 * @param max_qp_iterations  Per-tick QP budget (0 → default)
 * @return OffsetFreeMPCArtifacts for damp::OffsetFreeMPC
 */
template<
    size_t NP, size_t NC, auto Periods = 0.0, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NXM,
    size_t ND, size_t NZ, size_t NXQ = NXM>
[[nodiscard]] constexpr OffsetFreeMPCArtifacts<NX, NU, NY, NP, NC, T, ND, NZ, ::damp::detail::whole_periods(Periods)> mpc(
    const StateSpace<NX, NU, NY, T, NW, NV>&       sys,
    const OffsetFreeModel<NXM, NU, NY, ND, NZ, T>& model,
    const MPCWeights<NU, NZ, T>&                   weights = {},
    const MPCConstraints<NU, NZ, T>&               constraints = {},
    const OffsetFreeNoise<NXQ, NY, ND, T>&         noise = {},
    size_t                                         max_qp_iterations = 0
) {
    static_assert(std::is_floating_point_v<std::remove_cvref_t<decltype(Periods)>>, "design::mpc<NP, NC, Periods>: give the input delay in periods as a floating value (2.0, or kDeadTime / kTs).");
    static_assert(Periods >= 0, "design::mpc: the input delay is not negative.");
    static_assert(
        !::damp::detail::has_fractional_period(Periods),
        "design::mpc: a discrete model is delayed by whole periods; pass the continuous plant and Ts for a fraction of a period."
    );
    constexpr size_t D = ::damp::detail::whole_periods(Periods);
    static_assert(ND >= 1, "An offset-free model estimates at least one disturbance.");
    static_assert(ND <= NY, "Estimating ND disturbances needs at least ND independent measurements (ND <= NY).");
    static_assert(NXM <= NX && NXQ <= NX, "The model and the process noise cover the plant: the leading states of sys.");

    constexpr size_t NXK = NX + ND;

    OffsetFreeMPCArtifacts<NX, NU, NY, NP, NC, T, ND, NZ, D> art{};

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

    // Steady-state target maps: [A−I, B; C_z, D_z][x_t; u_t] = [−B_d d; r], solved once for
    // unit d and r (least squares when NZ ≠ NU). The runtime regulates about this target, so
    // the QP gradient is formed from small deviations however far the plant sits from zero.
    Matrix<NX + NZ, NX + NU, T> M{};
    M.template block<NX, NX>(0, 0) = sys.A - Matrix<NX, NX, T>::identity();
    M.template block<NX, NU>(0, NX) = sys.B;
    M.template block<NZ, NX>(NX, 0) = Cz;
    M.template block<NZ, NU>(NX, NX) = model.Dz;
    if constexpr (NZ == NU) {
        if (mat::rank(M) != NX + NU) {
            return art; // no unique steady state for these outputs (as discrete_oflqg refuses)
        }
    }
    const Matrix<NX + NU, NX + NZ, T> Mp = mat::pseudo_inverse(M);
    Matrix<NX + NZ, ND, T>            rhs_d{};
    rhs_d.template block<NX, ND>(0, 0) = -Bd;
    Matrix<NX + NZ, NZ, T> rhs_r{};
    rhs_r.template block<NZ, NZ>(NX, 0) = Matrix<NZ, NZ, T>::identity();
    art.target_d = Mp * rhs_d;
    art.target_r = Mp * rhs_r;

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
    Q_aug.template block<NXQ, NXQ>(0, 0) = noise.process; // delay states, if any: no process noise
    Q_aug.template block<ND, ND>(NX, NX) = noise.disturbance;

    art.kalman = kalman(sys_aug, Q_aug, noise.measurement);

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
 * @tparam NP       Prediction horizon (steps)
 * @tparam NC       Control horizon (moves, NC ≤ NP)
 * @tparam Periods  Input delay the runtime predicts, in whole periods, as a floating value (default 0.0)
 * @param sys             Discrete-time plant (Ts > 0, D = 0)
 * @param weights         MPC cost weights (default: Qy = I, RΔu = I)
 * @param constraints     MPC box limits (default: unbounded)
 * @param noise           Process, measurement and disturbance covariances (default I). On a
 *                        model with delay states (discretize with InputDelay / OutputDelay)
 *                        give the process noise plant-sized: the delay states get none.
 * @param max_qp_iterations  Per-tick QP budget (0 → default)
 * @return OffsetFreeMPCArtifacts for damp::OffsetFreeMPC
 */
template<size_t NP, size_t NC, auto Periods = 0.0, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NXQ = NX>
[[nodiscard]] constexpr OffsetFreeMPCArtifacts<NX, NU, NY, NP, NC, T, NU, NY, ::damp::detail::whole_periods(Periods)> mpc(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const MPCWeights<NU, NY, T>&             weights = {},
    const MPCConstraints<NU, NY, T>&         constraints = {},
    const OffsetFreeNoise<NXQ, NY, NU, T>&   noise = {},
    size_t                                   max_qp_iterations = 0
) {
    static_assert(NU <= NY, "Estimating NU input disturbances needs at least NU independent measurements (NU <= NY).");
    const OffsetFreeModel<NX, NU, NY, NU, NY, T> model{.Bd = sys.B, .Cz = sys.C, .Dz = sys.D};
    return mpc<NP, NC, Periods>(sys, model, weights, constraints, noise, max_qp_iterations);
}

/**
 * @brief Synthesize an offset-free constrained MPC from the continuous plant, its input delay
 *        in sampling periods
 *
 * Write the delay as @f$ d + f @f$ periods, @f$ 0 \le f < 1 @f$. The d whole periods are
 * predicted at run time (InputDelayPredictor); a fraction f puts one register state in the
 * model (discretize with PredictedInputDelay), and a whole delay needs none. The process
 * noise is sized to the plant.
 *
 * @code
 * constexpr auto art = design::mpc<30, 10, kDeadTime / kTs>(  // 5 s at Ts = 2 s: 2.5 periods
 *     plant, kTs, {.Rdu = {{10.0}}}, {.u_min = {0.0}, .u_max = {40.0}});
 * static constinit OffsetFreeMPC law{art.as<float>()};       // sizes and the delay from art
 * @endcode
 *
 * @tparam NP       Prediction horizon (steps)
 * @tparam NC       Control horizon (moves, NC ≤ NP)
 * @tparam Periods  Input delay in sampling periods, as a floating value (default 0.0: none)
 * @param plant           Continuous-time plant (or a discrete one at Ts, for a whole delay)
 * @param sampling_time   Control period Ts [s]
 * @return OffsetFreeMPCArtifacts on NX states, NX + NU with a fraction of a period;
 *         success = false if the plant does not discretize.
 */
template<size_t NP, size_t NC, auto Periods = 0.0, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr OffsetFreeMPCArtifacts<
    NX + (::damp::detail::has_fractional_period(Periods) ? NU : 0), NU, NY, NP, NC, T, NU, NY, ::damp::detail::whole_periods(Periods)>
mpc(const StateSpace<NX, NU, NY, T, NW, NV>& plant,
    T                                        sampling_time,
    const MPCWeights<NU, NY, T>&             weights = {},
    const MPCConstraints<NU, NY, T>&         constraints = {},
    const OffsetFreeNoise<NX, NY, NU, T>&    noise = {},
    size_t                                   max_qp_iterations = 0) {
    static_assert(std::is_floating_point_v<std::remove_cvref_t<decltype(Periods)>>, "design::mpc<NP, NC, Periods>: give the input delay in periods as a floating value (2.0, or kDeadTime / kTs).");
    static_assert(Periods >= 0, "design::mpc: the input delay is not negative.");
    constexpr size_t D = ::damp::detail::whole_periods(Periods);
    constexpr double whole = static_cast<double>(D);
    if constexpr (::damp::detail::has_fractional_period(Periods)) {
        const auto model = discretize(plant, sampling_time, PredictedInputDelay<D, T>{static_cast<T>(Periods) * sampling_time});
        if (!model) {
            return {};
        }
        return mpc<NP, NC, whole>(*model, weights, constraints, noise, max_qp_iterations);
    } else {
        const auto model = discretize(plant, sampling_time);
        if (!model) {
            return {};
        }
        return mpc<NP, NC, whole>(*model, weights, constraints, noise, max_qp_iterations);
    }
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
 * The estimator holds deviations from an anchor, the steady-state target (x_t, u_t) for the
 * current reference under the estimated disturbance, re-anchored each tick: measurements
 * enter as y − C x_a − C_d d_a − D u_a and inputs as u − u_a, so the state the filter stores
 * and the gradient the QP forms stay small however far the plant sits from zero (a
 * temperature at 200 °C, a position in metres at micrometre precision). Equal to the
 * absolute-coordinate filter in exact arithmetic; in float it is the difference between a
 * held input that is constant and one that quantizes. The same device as OFLQG's error
 * coordinates, anchored on the target instead of a supplied trajectory.
 *
 * @tparam D Whole periods of input delay handled by prediction (InputDelayPredictor): design
 *           with the delay in periods (design::mpc<NP, NC, kDeadTime / kTs>(plant, Ts, …));
 *           OffsetFreeMPC law{art} deduces it.
 *
 * @note Compare with MATLAB®'s mpcmove(mpcobj, xc, ym, r) with the default
 *       built-in estimator.
 * @see design::mpc
 */
template<size_t NX, size_t NU, size_t NY, size_t NP, size_t NC, typename T = float, size_t ND = NU, size_t NZ = NY, size_t D = 0>
struct OffsetFreeMPC {
    static constexpr bool   output_feedback_tick = true;
    static constexpr size_t NXK = NX + ND;
    using Artifacts = design::OffsetFreeMPCArtifacts<NX, NU, NY, NP, NC, T, ND, NZ, D>;

    /// The controlled-output reference over the prediction horizon (MATLAB's `ref` with preview).
    using Preview = typename MPC<NX, NU, NZ, NP, NC, ND, T>::Preview;

    MPC<NX, NU, NZ, NP, NC, ND, T>                   mpc{};
    SteadyStateKalmanFilter<NXK, NU, NY, T, NXK, NY> kf{};       ///< Fixed-L estimator on [x̂ − x_a; d̂ − d_a]: deviations from the anchor
    Matrix<NX + NU, ND, T>                           target_d{}; ///< [x_t; u_t] per unit d̂
    Matrix<NX + NU, NZ, T>                           target_r{}; ///< [x_t; u_t] per unit r
    InputDelayPredictor<NXK, NU, (D > 0) ? D : 1, T> delay{};    ///< the commands in flight (D > 0)
    ColVec<NX, T>                                    x_a{};      ///< anchor: the steady-state target the estimate is relative to
    ColVec<NU, T>                                    u_a{};      ///< its input
    ColVec<ND, T>                                    d_a{};      ///< and the disturbance it holds under

    constexpr OffsetFreeMPC() = default;
    /// From the design. A failed design is a compile error when constant-evaluated
    /// (`static constinit OffsetFreeMPC law{art};`), so no separate static_assert is needed.
    constexpr explicit OffsetFreeMPC(const Artifacts& artifacts)
        : mpc(artifacts.mpc), kf(artifacts.kalman), target_d(artifacts.target_d), target_r(artifacts.target_r), delay(artifacts.kalman.sys) {
        if (!artifacts.success && std::is_constant_evaluated()) {
            detail::refuse_failed_design();
        }
    }

    /// Advance the estimator with the applied input (caller-sequenced use; absolute units).
    constexpr void predict(const ColVec<NU, T>& u) { kf.predict(ColVec<NU, T>(u - u_a)); }

    /// Correct the estimator with a measurement (caller-sequenced use; absolute units).
    constexpr bool update(const ColVec<NY, T>& y, const ColVec<NU, T>& u = ColVec<NU, T>{}) {
        return kf.update(measured_deviation(y), ColVec<NU, T>(u - u_a));
    }

    /// MPC move at the current estimate (call predict/update first). @p r is the NZ controlled-output reference, not x_ref.
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NZ, T>& r) { return control(mpc.held(r)); }

    /// @copydoc control(const ColVec<NZ,T>&)
    /// With the reference known over the horizon.
    [[nodiscard]] constexpr ColVec<NU, T> control(const Preview& r) {
        const ColVec<NU, T> u = feedback_at_estimate(r);
        mpc.commit(u);
        return u;
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
    constexpr ColVec<NU, T> control(const ColVec<NZ, T>& r, const ColVec<NY, T>& y) { return control(mpc.held(r), y); }

    /// Single-input single-output tick: scalar reference and measurement, scalar input.
    constexpr T control(T r, T y)
        requires(NU == 1 && NY == 1 && NZ == 1)
    {
        return control(ColVec<1, T>{r}, ColVec<1, T>{y})[0];
    }

    /**
     * @brief Full tick with the reference known over the horizon (preview)
     *
     * @p r holds the controlled-output reference at each of the NP prediction steps, so the
     * plan sees a setpoint change before it happens: a trajectory generator hands in its next
     * NP points each tick, repeating the last known value to fill. The single-reference tick
     * is this one with NP copies of r.
     */
    constexpr ColVec<NU, T> control(const Preview& r, const ColVec<NY, T>& y) {
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
    [[nodiscard]] constexpr ColVec<NU, T> feedback(const ColVec<NZ, T>& r, const ColVec<NY, T>& y) { return feedback(mpc.held(r), y); }

    /// Single-input single-output first half (see feedback).
    [[nodiscard]] constexpr T feedback(T r, T y)
        requires(NU == 1 && NY == 1 && NZ == 1)
    {
        return feedback(ColVec<1, T>{r}, ColVec<1, T>{y})[0];
    }

    /// @copydoc feedback(const ColVec<NZ,T>&, const ColVec<NY,T>&)
    /// With the reference known over the horizon (see control(const Preview&, y)).
    [[nodiscard]] constexpr ColVec<NU, T> feedback(const Preview& r, const ColVec<NY, T>& y) {
        if constexpr (D > 0) {
            predict(delay.applied_last()); // the plant moved on the command that reached it
            update(y, delay.applied_now());
        } else {
            predict(mpc.previous_control());
            update(y, mpc.previous_control());
        }
        return feedback_at_estimate(r);
    }

    /// Steady-state pair (x_t, u_t) that holds z = r under the estimated disturbance.
    constexpr void equilibrium(const ColVec<NZ, T>& r, const ColVec<ND, T>& d, ColVec<NX, T>& xt, ColVec<NU, T>& ut) const {
        const ColVec<NX + NU, T> s = (target_d * d) + (target_r * r);
        for (size_t i = 0; i < NX; ++i) {
            xt(i) = s(i);
        }
        for (size_t i = 0; i < NU; ++i) {
            ut(i) = s(NX + i);
        }
    }

    /// Second half of the tick: the input the plant actually received (post-clamp).
    constexpr void commit(const ColVec<NU, T>& u_applied) {
        mpc.commit(u_applied);
        if constexpr (D > 0) {
            delay.push(u_applied);
        }
    }

    /// Single-input second half: the input the plant actually received (see commit).
    constexpr void commit(T u_applied)
        requires(NU == 1)
    {
        commit(ColVec<1, T>{u_applied});
    }

    /// Estimated plant state x̂ (absolute: the anchor plus the deviation the filter holds)
    [[nodiscard]] constexpr ColVec<NX, T> state_estimate() const {
        ColVec<NX, T> x{};
        for (size_t i = 0; i < NX; ++i) {
            x(i) = kf.state()(i) + x_a(i);
        }
        return x;
    }

    /// Estimated disturbance d̂ (what the integrating model has absorbed)
    [[nodiscard]] constexpr ColVec<ND, T> disturbance_estimate() const {
        ColVec<ND, T> d{};
        for (size_t j = 0; j < ND; ++j) {
            d(j) = kf.state()(NX + j) + d_a(j);
        }
        return d;
    }

    /// Reset the controller memory; the state estimate and its anchor are left in place
    constexpr void reset() {
        mpc.reset();
        delay.reset();
    }

private:
    /// y − C x_a − C_d d_a − D u_a: the measurement the deviation state explains.
    [[nodiscard]] constexpr ColVec<NY, T> measured_deviation(const ColVec<NY, T>& y) const {
        const auto&    m = kf.model();
        ColVec<NXK, T> anchor{};
        for (size_t i = 0; i < NX; ++i) {
            anchor(i) = x_a(i);
        }
        for (size_t j = 0; j < ND; ++j) {
            anchor(NX + j) = d_a(j);
        }
        return ColVec<NY, T>(y - (m.C * anchor) - (m.D * u_a));
    }

    /// Move the anchor to the target for @p r under the current d̂: the deviations shift by the
    /// anchor change, exactly (no model involved between two anchors at the same instant), and
    /// the disturbance deviation returns to zero.
    constexpr void anchor_on_target(const ColVec<NZ, T>& r) {
        const ColVec<ND, T> d_hat = disturbance_estimate();
        ColVec<NX, T>       x_t{};
        ColVec<NU, T>       u_t{};
        equilibrium(r, d_hat, x_t, u_t);
        ColVec<NXK, T> s = kf.state();
        for (size_t i = 0; i < NX; ++i) {
            s(i) += x_a(i) - x_t(i);
        }
        for (size_t j = 0; j < ND; ++j) {
            s(NX + j) = T{0};
        }
        kf.set_state(s);
        x_a = x_t;
        u_a = u_t;
        d_a = d_hat;
    }

    /// The deviation the command decided now will meet: the filter's, predicted through the
    /// commands in flight (relative to u_a; the anchor is the equilibrium for d_a = d̂) when D > 0.
    [[nodiscard]] constexpr ColVec<NX, T> deviation_at_command() const {
        ColVec<NXK, T> e = kf.state();
        if constexpr (D > 0) {
            e = delay.predict(e, u_a);
        }
        ColVec<NX, T> out{};
        for (size_t i = 0; i < NX; ++i) {
            out(i) = e(i);
        }
        return out;
    }

    /// Re-anchor on the target for the reference now (r[0]), then the MPC move from the
    /// deviation with the whole preview (see MPC::feedback_deviation). With D periods
    /// predicted the plan starts at k + D, so the preview it sees is r shifted by D, the last
    /// value held: the command decided now meets the plant D periods from now.
    [[nodiscard]] constexpr ColVec<NU, T> feedback_at_estimate(const Preview& r) {
        anchor_on_target(r[0]);
        if constexpr (D > 0) {
            Preview ahead{};
            for (size_t j = 0; j < NP; ++j) {
                ahead[j] = r[(j + D < NP) ? (j + D) : (NP - 1)];
            }
            return mpc.feedback_deviation(ahead, deviation_at_command(), d_a, x_a, u_a);
        } else {
            return mpc.feedback_deviation(r, deviation_at_command(), d_a, x_a, u_a);
        }
    }
};

/// Deduce the runtime from its artifacts: OffsetFreeMPC controller{art};
template<size_t NX, size_t NU, size_t NY, size_t NP, size_t NC, typename T, size_t ND, size_t NZ, size_t D>
OffsetFreeMPC(const design::OffsetFreeMPCArtifacts<NX, NU, NY, NP, NC, T, ND, NZ, D>&) -> OffsetFreeMPC<NX, NU, NY, NP, NC, T, ND, NZ, D>;

} // namespace damp
