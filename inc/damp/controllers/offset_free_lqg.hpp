// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file offset_free_lqg.hpp
 * @brief Offset-free LQG: LQR + disturbance-augmented Kalman + steady-state target
 *
 * Packs the unconstrained infinite-horizon form of OffsetFreeMPC: discrete LQR on the
 * plant, a Kalman filter on a disturbance augment
 * @f[
 *   x_{k+1} = A x_k + B u_k + B_d d_k,\qquad d_{k+1} = d_k,\qquad
 *   y_k = C x_k + C_d d_k + D u_k,
 * @f]
 * and a steady-state target that holds the controlled outputs @f$ z = C_z x + D_z u @f$ at
 * @f$ r @f$ under the estimated disturbance:
 * @f[
 *   u = u_t - K(\hat x - x_t),\qquad
 *   \begin{bmatrix} A-I & B \\ C_z & D_z \end{bmatrix}
 *   \begin{bmatrix} x_t \\ u_t \end{bmatrix}
 *   =
 *   \begin{bmatrix} -B_d\hat d \\ r \end{bmatrix}.
 * @f]
 *
 * Three choices make it fit a real plant, all in OffsetFreeModel:
 * - where the disturbance enters (@f$ B_d @f$): at the input (the default, Pannocchia &
 *   Rawlings' input disturbance) or where it physically acts, e.g. friction on a load
 *   behind a compliant transmission;
 * - which sensors feel it directly (@f$ C_d @f$), e.g. an accelerometer on that load;
 * - which outputs are held offset-free (@f$ C_z @f$), separately from the measured
 *   outputs @f$ C @f$: estimate with every sensor, hold only the ones that matter.
 *
 * The target matrix is constant, so its solution is precomputed at design time: the
 * runtime target is two small matrix–vector products, not a linear solve.
 *
 * Trajectory tracking: pass the reference state and input (x*, u*) each tick and the
 * controlled outputs follow @f$ C_z x^* + D_z u^* @f$. The filter then runs on the error
 * @f$ \tilde x = x - x^* @f$: measurements enter as @f$ y - C x^* - D u^* @f$ and the input
 * as @f$ u - u^* @f$, so its state holds small numbers however far the plant travels. When
 * the reference moves, the error is re-anchored by @f$ A x^*_k + B u^*_k - x^*_{k+1} @f$,
 * which makes it exactly the absolute-coordinate filter: the reference need not be a
 * trajectory of the model.
 *
 * Delays: design on the model from discretize(plant, Ts, InputDelay<D>{…}) and/or
 * OutputDelay<D>{}, with the OffsetFreeModel and the weights still sized to the plant. The
 * tracking form then takes the plant-state reference x* alone and keeps the reference for
 * the delay states (the past u*, the outputs x* produced) itself.
 *
 * @note Compare with MATLAB®'s lqg + setindist input-disturbance estimator,
 *       without a QP. Constrained / preview use OffsetFreeMPC.
 * @see controllers/lqg.hpp, controllers/offset_free_mpc.hpp
 * @see Pannocchia & Rawlings, "Disturbance Models for Offset-Free
 *      Model-Predictive Control," AIChE Journal 49(2), 2003
 *
 * @code
 * constexpr auto res = design::discrete_oflqg(sys_d, Q, R, {.process = Qn, .measurement = Rn});  // discrete plant; input disturbance, z = y
 * static constinit OFLQG ctrl{res};  // a failed design does not compile
 * // each tick: auto u = ctrl.control(r, y);
 * @endcode
 */

#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/delay_predictor.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/matrix/block.hpp"
#include "damp/matrix/functions.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/solve.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "lqr.hpp"

namespace damp {

namespace design {

/**
 * @brief Where the disturbance acts, what sees it, and which outputs are held offset-free.
 * @tparam ND Disturbance count; NZ Controlled-output count (must equal NU).
 */
template<size_t NX, size_t NU, size_t NY, size_t ND, size_t NZ, typename T = double>
struct OffsetFreeModel {
    Matrix<NX, ND, T> Bd{}; ///< how the disturbance enters the discrete state update (e.g. G of the discretized plant)
    Matrix<NY, ND, T> Cd{}; ///< how it enters the measurements directly (e.g. an accelerometer)
    Matrix<NZ, NX, T> Cz{}; ///< controlled outputs, held at r offset-free
    Matrix<NZ, NU, T> Dz{}; ///< controlled-output feedthrough
};

/**
 * @brief The offset-free estimator's noise: on the plant states, the sensors, and the disturbance.
 * @tparam NX Plant states the process noise covers: the leading block on a model with delay states.
 */
template<size_t NX, size_t NY, size_t ND, typename T = double>
struct OffsetFreeNoise {
    Matrix<NX, NX, T> process = Matrix<NX, NX, T>::identity();     ///< process-noise covariance on x: model error
    Matrix<NY, NY, T> measurement = Matrix<NY, NY, T>::identity(); ///< measurement-noise covariance
    Matrix<ND, ND, T> disturbance = Matrix<ND, ND, T>::identity(); ///< random-walk covariance on d: how fast the disturbance may move
};

/// Sizes from the matrices: `OffsetFreeNoise{Qp, Rn, Qd}` (plant-sized Qp on a delayed model).
template<size_t NX, size_t NY, size_t ND, typename T>
OffsetFreeNoise(Matrix<NX, NX, T>, Matrix<NY, NY, T>, Matrix<ND, ND, T>) -> OffsetFreeNoise<NX, NY, ND, T>;

/**
 * @struct OFLQGResult
 * @brief Offset-free LQG design: LQR on the plant, Kalman on [x; d], precomputed target
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t ND = NU, size_t NZ = NY, size_t D = 0>
struct OFLQGResult {
    static constexpr size_t NXK = NX + ND; ///< Estimator state size [x; d]

    LQRResult<NX, NU, T>                  lqr{};      ///< Regulator on the original plant
    KalmanResult<NXK, NU, NY, T, NXK, NY> kalman{};   ///< Filter on the disturbance augment
    Matrix<NX + NU, ND, T>                target_d{}; ///< [x_t; u_t] per unit estimated disturbance
    Matrix<NX + NU, NZ, T>                target_r{}; ///< [x_t; u_t] per unit controlled-output reference
    bool                                  success{false};

    template<typename U>
    [[nodiscard]] constexpr OFLQGResult<NX, NU, NY, U, ND, NZ, D> as() const {
        return OFLQGResult<NX, NU, NY, U, ND, NZ, D>{
            lqr.template as<U>(),
            kalman.template as<U>(),
            target_d.template as<U>(),
            target_r.template as<U>(),
            success,
        };
    }
};

namespace detail {

/// Steady-state target: [A−I, B; Cz, Dz][x_t; u_t] = [−Bd d; r], solved once for unit d and r.
template<size_t NX, size_t NU, size_t ND, size_t NZ, typename T>
[[nodiscard]] constexpr bool oflqg_target(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, ND, T>& Bd,
    const Matrix<NZ, NX, T>& Cz,
    const Matrix<NZ, NU, T>& Dz,
    Matrix<NX + NU, ND, T>&  target_d,
    Matrix<NX + NU, NZ, T>&  target_r
) {
    Matrix<NX + NZ, NX + NU, T> M{};
    M.template block<NX, NX>(0, 0) = A - Matrix<NX, NX, T>::identity();
    M.template block<NX, NU>(0, NX) = B;
    M.template block<NZ, NX>(NX, 0) = Cz;
    M.template block<NZ, NU>(NX, NX) = Dz;
    Matrix<NX + NZ, ND, T> rhs_d{};
    rhs_d.template block<NX, ND>(0, 0) = -Bd;
    Matrix<NX + NZ, NZ, T> rhs_r{};
    rhs_r.template block<NZ, NZ>(NX, 0) = Matrix<NZ, NZ, T>::identity();
    const auto td = mat::solve(M, rhs_d);
    const auto tr = mat::solve(M, rhs_r);
    if (!td || !tr) {
        return false;
    }
    target_d = td.value();
    target_r = tr.value();
    return true;
}

} // namespace detail

/**
 * @brief Discrete offset-free LQG with a chosen disturbance model and controlled outputs
 *
 * LQR is designed on @p sys (Q on x only). The filter is designed on the disturbance
 * augment defined by @p model. Fails if the steady-state target is not unique
 * ([A−I, B; Cz, Dz] singular), or if either DARE fails (which includes an augment the
 * sensors cannot detect).
 *
 * On a model with delay states (discretize with InputDelay / OutputDelay), size @p model
 * and the weights to the plant (NXM states): they are the leading block, and the delay
 * states get no cost, no process noise and no disturbance. The design is otherwise the
 * one for the undelayed plant, swapped onto the delayed model.
 *
 * @tparam Periods        Input delay the runtime predicts, in whole sampling periods, as a floating value
 *                        (default 0.0); a fraction of a period needs the continuous plant and Ts.
 *                        OFLQG law{res} deduces it
 * @param sys             Discrete-time plant (Ts > 0); C, D are the measured outputs
 * @param model           Disturbance entry Bd, its direct sensor coupling Cd, controlled outputs Cz, Dz
 * @param Q_lqr           State cost (NXM × NXM)
 * @param R_lqr           Input cost (NU × NU)
 * @param noise           Process (NXM × NXM), measurement and disturbance covariances (default I)
 * @param N               LQR cross term (default 0)
 */
template<auto Periods = 0.0, size_t NX, size_t NU, size_t NY, size_t ND, size_t NZ, size_t NXM, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr OFLQGResult<NX, NU, NY, T, ND, NZ, ::damp::detail::whole_periods(Periods)> discrete_oflqg(
    const StateSpace<NX, NU, NY, T, NW, NV>&       sys,
    const OffsetFreeModel<NXM, NU, NY, ND, NZ, T>& model,
    const Matrix<NXM, NXM, T>&                     Q_lqr,
    const Matrix<NU, NU, T>&                       R_lqr,
    const OffsetFreeNoise<NXM, NY, ND, T>&         noise = {},
    const Matrix<NXM, NU, T>&                      N = Matrix<NXM, NU, T>{}
) {
    static_assert(std::is_floating_point_v<std::remove_cvref_t<decltype(Periods)>>, "design::discrete_oflqg<Periods>: give the input delay in periods as a floating value (2.0, or kDeadTime / kTs).");
    static_assert(Periods >= 0, "design::discrete_oflqg: the input delay is not negative.");
    static_assert(
        !::damp::detail::has_fractional_period(Periods),
        "design::discrete_oflqg: a discrete model is delayed by whole periods; pass the continuous plant and Ts for a fraction of a period."
    );
    static_assert(NZ == NU, "OFLQG holds as many outputs offset-free as it has inputs (NZ = NU).");
    static_assert(NXM <= NX, "The model and weights cover the plant: the leading states of sys.");
    constexpr size_t                      D = ::damp::detail::whole_periods(Periods);
    constexpr size_t                      NXK = NX + ND;
    OFLQGResult<NX, NU, NY, T, ND, NZ, D> art{};
    // Delay states, if any, follow the plant's: no disturbance, no controlled output, no cost.
    const auto Bd = damp::detail::pad_leading<NX, ND>(model.Bd);
    const auto Cz = damp::detail::pad_leading<NZ, NX>(model.Cz);
    const auto Qx = damp::detail::pad_leading<NX, NX>(Q_lqr);
    const auto Nx = damp::detail::pad_leading<NX, NU>(N);
    if (!detail::oflqg_target(sys.A, sys.B, Bd, Cz, model.Dz, art.target_d, art.target_r)) {
        return art;
    }
    art.lqr = discrete_lqr(sys.A, sys.B, Qx, R_lqr, Nx);

    StateSpace<NXK, NU, NY, T, NXK, NY> sys_aug{};
    sys_aug.A.template block<NX, NX>(0, 0) = sys.A;
    sys_aug.A.template block<NX, ND>(0, NX) = Bd;
    sys_aug.A.template block<ND, ND>(NX, NX) = Matrix<ND, ND, T>::identity();
    sys_aug.B.template block<NX, NU>(0, 0) = sys.B;
    sys_aug.C.template block<NY, NX>(0, 0) = sys.C;
    sys_aug.C.template block<NY, ND>(0, NX) = model.Cd;
    sys_aug.D = sys.D;
    sys_aug.Ts = (sys.Ts > T{0}) ? sys.Ts : T{1};

    Matrix<NXK, NXK, T> Q_aug{};
    Q_aug.template block<NXM, NXM>(0, 0) = noise.process;
    Q_aug.template block<ND, ND>(NX, NX) = noise.disturbance;
    art.kalman = kalman(sys_aug, Q_aug, noise.measurement);
    art.kalman.sys.Ts = sys.Ts;
    art.success = art.lqr.success && art.kalman.success;
    return art;
}

/**
 * @brief Discrete offset-free LQG, input disturbance, every measured output held at r
 *
 * The classic square form: d enters with the input (Bd = B) but not the sensors (Cd = 0), and
 * z = y (Cz = C, Dz = D).
 * Requires NU = NY. On a model with delay states, @p Q_lqr and the process noise may be sized
 * to the plant (NXQ < NX): the delay states that follow get no cost and no process noise.
 */
template<auto Periods = 0.0, size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0, size_t NXQ = NX>
[[nodiscard]] constexpr OFLQGResult<NX, NU, NY, T, NU, NY, ::damp::detail::whole_periods(Periods)> discrete_oflqg(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NXQ, NXQ, T>&               Q_lqr,
    const Matrix<NU, NU, T>&                 R_lqr,
    const OffsetFreeNoise<NXQ, NY, NU, T>&   noise = {},
    const Matrix<NXQ, NU, T>&                N = Matrix<NXQ, NU, T>{}
) {
    static_assert(NU == NY, "The square OFLQG form needs NU = NY; pass an OffsetFreeModel to pick NU controlled outputs.");
    static_assert(NXQ <= NX, "The weights cover the plant: the leading states of sys.");
    const OffsetFreeModel<NX, NU, NY, NU, NY, T> model{.Bd = sys.B, .Cz = sys.C, .Dz = sys.D};
    const OffsetFreeNoise<NX, NY, NU, T>         padded{
                .process = damp::detail::pad_leading<NX, NX>(noise.process),
                .measurement = noise.measurement,
                .disturbance = noise.disturbance,
    };
    return discrete_oflqg<Periods>(sys, model, damp::detail::pad_leading<NX, NX>(Q_lqr), R_lqr, padded, damp::detail::pad_leading<NX, NU>(N));
}

/**
 * @brief Offset-free LQG from the continuous plant, its input delay in sampling periods
 *
 * Write the delay as @f$ d + f @f$ periods, @f$ 0 \le f < 1 @f$. The d whole periods are
 * predicted at run time (InputDelayPredictor); a fraction f puts one register state in the
 * model (discretize with PredictedInputDelay), and a whole delay needs none. @p Q_lqr, the
 * process noise and @p N are sized to the plant; the register state gets no cost and no
 * process noise.
 *
 * @code
 * constexpr auto res = design::discrete_oflqg<kDeadTime / kTs>(plant, kTs, Q, R, {.process = Qn, .measurement = Rn});
 * static constinit OFLQG law{res.as<float>()};  // sizes and the delay from res
 * @endcode
 *
 * @tparam Periods  Input delay in sampling periods, as a floating value (default 0.0: none)
 * @return OFLQGResult on NX states, NX + NU with a fraction of a period; success = false if
 *         the plant does not discretize.
 */
template<auto Periods = 0.0, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr OFLQGResult<
    NX + (::damp::detail::has_fractional_period(Periods) ? NU : 0), NU, NY, T, NU, NY, ::damp::detail::whole_periods(Periods)>
discrete_oflqg(
    const StateSpace<NX, NU, NY, T, NW, NV>& plant,
    T                                        sampling_time,
    const Matrix<NX, NX, T>&                 Q_lqr,
    const Matrix<NU, NU, T>&                 R_lqr,
    const OffsetFreeNoise<NX, NY, NU, T>&    noise = {},
    const Matrix<NX, NU, T>&                 N = Matrix<NX, NU, T>{}
) {
    static_assert(std::is_floating_point_v<std::remove_cvref_t<decltype(Periods)>>, "design::discrete_oflqg<Periods>: give the input delay in periods as a floating value (2.0, or kDeadTime / kTs).");
    static_assert(Periods >= 0, "design::discrete_oflqg: the input delay is not negative.");
    constexpr size_t D = ::damp::detail::whole_periods(Periods);
    constexpr double whole = static_cast<double>(D);
    if constexpr (::damp::detail::has_fractional_period(Periods)) {
        const auto model = discretize(plant, sampling_time, PredictedInputDelay<D, T>{static_cast<T>(Periods) * sampling_time});
        if (!model) {
            return {};
        }
        return discrete_oflqg<whole>(*model, Q_lqr, R_lqr, noise, N);
    } else {
        const auto model = discretize(plant, sampling_time);
        if (!model) {
            return {};
        }
        return discrete_oflqg<whole>(*model, Q_lqr, R_lqr, noise, N);
    }
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Offset-free LQG: u = u_t − K(x̂ − x_t) with d̂ in the steady-state target
 *
 * control(r, y) is a self-contained tick (OutputFeedbackController); control(y, x*, u*)
 * tracks a reference trajectory, with x* on the plant states only when the model carries
 * delay states. Prefer feedback/commit when the applied input can saturate: commit takes
 * the input actually applied.
 *
 * @tparam D Whole samples of input delay handled by prediction (InputDelayPredictor): the
 *           estimator runs on the command the plant actually received and the law on the
 *           state the new command will meet. Design with the delay in periods
 *           (discrete_oflqg<kDeadTime / kTs>(plant, Ts, …)); OFLQG law{res} deduces it. Setpoint form only; the tracking form needs D = 0
 *           (put the delay in the model with InputDelay instead).
 */
template<size_t NX, size_t NU, size_t NY, typename T = float, size_t ND = NU, size_t NZ = NY, size_t D = 0>
struct OFLQG {
    static constexpr bool   output_feedback_tick = true;
    static constexpr size_t NXK = NX + ND;
    using Result = design::OFLQGResult<NX, NU, NY, T, ND, NZ, D>;

    StateFeedback<NX, NU, T>                         lqr{};
    SteadyStateKalmanFilter<NXK, NU, NY, T, NXK, NY> kf{};
    Matrix<NX + NU, ND, T>                           target_d{};
    Matrix<NX + NU, NZ, T>                           target_r{};
    ColVec<NU, T>                                    u_prev{}; ///< last applied input, relative to its reference
    InputDelayPredictor<NXK, NU, (D > 0) ? D : 1, T> delay{};  ///< the commands in flight (D > 0)

    constexpr OFLQG() = default;

    template<typename U>
    constexpr OFLQG(const design::OFLQGResult<NX, NU, NY, U, ND, NZ, D>& result) // NOLINT
        : lqr(result.lqr),
          kf(result.kalman),
          target_d(result.target_d.template as<T>()),
          target_r(result.target_r.template as<T>()),
          delay(result.kalman.sys) {
        if (!result.success && std::is_constant_evaluated()) {
            detail::refuse_failed_design(); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    constexpr void predict(const ColVec<NU, T>& u = ColVec<NU, T>{}) { kf.predict(u); }
    constexpr bool update(const ColVec<NY, T>& y, const ColVec<NU, T>& u = ColVec<NU, T>{}) {
        return kf.update(y, u);
    }

    /// Estimated plant state (first NX of the augmented filter), relative to the current reference
    [[nodiscard]] constexpr ColVec<NX, T> state_estimate() const {
        ColVec<NX, T> x{};
        for (size_t i = 0; i < NX; ++i) {
            x(i) = kf.state()(i);
        }
        return x;
    }

    /// Estimated disturbance d̂
    [[nodiscard]] constexpr ColVec<ND, T> disturbance_estimate() const {
        ColVec<ND, T> d{};
        for (size_t j = 0; j < ND; ++j) {
            d(j) = kf.state()(NX + j);
        }
        return d;
    }

    /// Steady-state pair (x_t, u_t) that holds z = r under disturbance @p d.
    constexpr void equilibrium(const ColVec<NZ, T>& r, const ColVec<ND, T>& d, ColVec<NX, T>& xt, ColVec<NU, T>& ut) const {
        const ColVec<NX + NU, T> s = (target_d * d) + (target_r * r);
        for (size_t i = 0; i < NX; ++i) {
            xt(i) = s(i);
        }
        for (size_t i = 0; i < NU; ++i) {
            ut(i) = s(NX + i);
        }
    }

    /// The plant state the command decided now will meet: the estimate, predicted through the
    /// commands in flight when D > 0.
    [[nodiscard]] constexpr ColVec<NX, T> state_at_command() const {
        if constexpr (D > 0) {
            const ColVec<NXK, T> xp = delay.predict(kf.state());
            ColVec<NX, T>        x{};
            for (size_t i = 0; i < NX; ++i) {
                x(i) = xp(i);
            }
            return x;
        } else {
            return state_estimate();
        }
    }

    /// Regulator at the current estimate: u = u_t − K(x̂ − x_t), x̂ predicted to when the command lands
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NZ, T>& r) const {
        ColVec<NX, T> xt{};
        ColVec<NU, T> ut{};
        equilibrium(r, disturbance_estimate(), xt, ut);
        return ut + lqr.control(xt, state_at_command());
    }

    /// Setpoint tick, first half: update on @p y, return the input that holds z at @p r.
    [[nodiscard]] constexpr ColVec<NU, T> feedback(const ColVec<NZ, T>& r, const ColVec<NY, T>& y) {
        move_reference(ColVec<NX, T>{}, ColVec<NU, T>{});
        if constexpr (D > 0) {
            kf.update(y, delay.applied_now()); // feedthrough sees the command that reached the plant
        } else {
            kf.update(y, u_prev);
        }
        return control(r);
    }

    /// Tracking tick, first half: update on @p y, return the input that keeps z on the
    /// reference trajectory (@p x_ref, @p u_ref).
    [[nodiscard]] constexpr ColVec<NU, T> feedback(const ColVec<NY, T>& y, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref) {
        return track(y, x_ref, u_ref);
    }

    /**
     * @brief Tracking tick on a model with delay states: the reference covers the plant only.
     *
     * @p x_ref holds the NR plant states; the delay states that follow them (discretize with
     * InputDelay / OutputDelay) are filled in from the reference history, the past u* and
     * the outputs the reference produced, so the caller never keeps that history.
     */
    template<size_t NR>
        requires(NR < NX)
    [[nodiscard]] constexpr ColVec<NU, T> feedback(const ColVec<NY, T>& y, const ColVec<NR, T>& x_ref, const ColVec<NU, T>& u_ref) {
        const auto& m = kf.model();
        return track(y, detail::extend_reference<NX>(m.A, m.B, x_ref_, u_ref_, x_ref), u_ref);
    }

    /// Second half of either tick: the input actually applied (after any clamp).
    constexpr void commit(const ColVec<NU, T>& u_applied) {
        u_prev = u_applied - u_ref_;
        if constexpr (D > 0) {
            kf.predict(delay.applied_now()); // the plant moves on the command that reached it
            delay.push(u_prev);
        } else {
            kf.predict(u_prev);
        }
    }

    /// Self-contained tick from the measurement (and the controlled-output reference r)
    constexpr ColVec<NU, T> control(const ColVec<NZ, T>& r, const ColVec<NY, T>& y) {
        const auto u = feedback(r, y);
        commit(u);
        return u;
    }

    /// Single-input single-output setpoint tick: scalar reference and measurement, scalar input.
    constexpr T control(T r, T y)
        requires(NU == 1 && NY == 1 && NZ == 1)
    {
        return control(ColVec<1, T>{r}, ColVec<1, T>{y})[0];
    }

    /// Single-input single-output setpoint tick, first half (see feedback).
    [[nodiscard]] constexpr T feedback(T r, T y)
        requires(NU == 1 && NY == 1 && NZ == 1)
    {
        return feedback(ColVec<1, T>{r}, ColVec<1, T>{y})[0];
    }

    /// Single-input second half: the input actually applied (see commit).
    constexpr void commit(T u_applied)
        requires(NU == 1)
    {
        commit(ColVec<1, T>{u_applied});
    }

    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NZ, T>& r, const ColVec<NY, T>& y) {
        return control(r, y);
    }

    /// Self-contained tracking tick
    constexpr ColVec<NU, T> control(const ColVec<NY, T>& y, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref) {
        const auto u = feedback(y, x_ref, u_ref);
        commit(u);
        return u;
    }

    /// Self-contained tracking tick, reference on the plant states only (see feedback).
    template<size_t NR>
        requires(NR < NX)
    constexpr ColVec<NU, T> control(const ColVec<NY, T>& y, const ColVec<NR, T>& x_ref, const ColVec<NU, T>& u_ref) {
        const auto u = feedback(y, x_ref, u_ref);
        commit(u);
        return u;
    }

    /// Clear the applied-input memory and the commands in flight (the estimator state, and the reference it is relative to, are left untouched).
    constexpr void reset() {
        u_prev = ColVec<NU, T>{};
        delay.reset();
    }

private:
    ColVec<NX, T> x_ref_{}; ///< reference the error state is currently relative to
    ColVec<NU, T> u_ref_{};

    /// Tracking update on the full reference: measurement in error coordinates, then the law.
    [[nodiscard]] constexpr ColVec<NU, T> track(const ColVec<NY, T>& y, const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref) {
        static_assert(D == 0, "OFLQG tracking with a predicted delay is not supported: put the delay in the model (InputDelay) instead.");
        move_reference(x_ref, u_ref);
        const auto&   m = kf.model();
        ColVec<NY, T> y_err = y - (m.D * u_ref);
        for (size_t i = 0; i < NY; ++i) {
            for (size_t j = 0; j < NX; ++j) {
                y_err(i) -= m.C(i, j) * x_ref(j);
            }
        }
        kf.update(y_err, u_prev);
        return u_ref + control(ColVec<NZ, T>{});
    }

    /// Re-anchor the predicted error on a new reference: x̃ += A x*ₖ + B u*ₖ − x*ₖ₊₁.
    constexpr void move_reference(const ColVec<NX, T>& x_ref, const ColVec<NU, T>& u_ref) {
        const auto&    m = kf.model();
        ColVec<NXK, T> xk = kf.state();
        for (size_t i = 0; i < NX; ++i) {
            T shift = -x_ref(i);
            for (size_t j = 0; j < NX; ++j) {
                shift += m.A(i, j) * x_ref_(j);
            }
            for (size_t j = 0; j < NU; ++j) {
                shift += m.B(i, j) * u_ref_(j);
            }
            xk(i) += shift;
        }
        kf.set_state(xk);
        x_ref_ = x_ref;
        u_ref_ = u_ref;
    }
};

template<size_t NX, size_t NU, size_t NY, typename T, size_t ND, size_t NZ, size_t D>
OFLQG(const design::OFLQGResult<NX, NU, NY, T, ND, NZ, D>&) -> OFLQG<NX, NU, NY, T, ND, NZ, D>;

} // namespace damp
