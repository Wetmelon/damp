// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file offset_free_lqg.hpp
 * @brief Offset-free LQG: LQR + input-disturbance Kalman + DC target
 *
 * Packs the unconstrained infinite-horizon form of OffsetFreeMPC: discrete
 * LQR on the original plant, a Kalman filter on the input-disturbance
 * augment
 * @f[
 *   x_{k+1} = A x_k + B(u_k + w_k),\qquad w_{k+1} = w_k,\qquad y_k = C x_k,
 * @f]
 * and a DC solve that uses @f$ \hat w @f$ so @f$ y \to r @f$ under a constant
 * load. @f$ K @f$ stays @f$ N_U \times N_X @f$; the bias never enters the
 * regulator. Runtime:
 * @f[
 *   u = u_t - K(\hat x - x_t),\qquad
 *   \begin{bmatrix} A-I & B \\ C & D \end{bmatrix}
 *   \begin{bmatrix} x_t \\ u_t \end{bmatrix}
 *   =
 *   \begin{bmatrix} -B\hat w \\ r \end{bmatrix}.
 * @f]
 *
 * @note Compare with MATLAB®'s lqg + setindist input-disturbance estimator,
 *       without a QP. Constrained / preview use OffsetFreeMPC.
 * @see controllers/lqg.hpp, controllers/offset_free_mpc.hpp
 * @see Pannocchia & Rawlings, "Disturbance Models for Offset-Free
 *      Model-Predictive Control," AIChE Journal 49(2), 2003
 *
 * @code
 * constexpr auto res = design::discrete_oflqg(sys, Q, R, Qn, Rn);
 * static_assert(res.success);
 * OFLQG ctrl{res};
 * auto u = ctrl.control(r, y);
 * @endcode
 */

#include <cstddef>

#include "damp/estimation/kalman.hpp"
#include "damp/matrix/block.hpp"
#include "damp/matrix/functions.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/solve.hpp"
#include "damp/systems/state_space.hpp"
#include "lqr.hpp"

namespace damp {

namespace design {

/**
 * @struct OFLQGResult
 * @brief Offset-free LQG design: LQR on the plant, Kalman on [x; w]
 */
template<size_t NX, size_t NU, size_t NY, typename T = double>
struct OFLQGResult {
    static constexpr size_t NXK = NX + NU; ///< Estimator state size [x; w]

    LQRResult<NX, NU, T>                  lqr{};    ///< Regulator on the original plant
    KalmanResult<NXK, NU, NY, T, NXK, NY> kalman{}; ///< Filter on the input-disturbance augment
    StateSpace<NX, NU, NY, T>             sys{};    ///< Discrete plant used for the DC target
    bool                                  success{false};

    template<typename U>
    [[nodiscard]] constexpr OFLQGResult<NX, NU, NY, U> as() const {
        return OFLQGResult<NX, NU, NY, U>{
            lqr.template as<U>(),
            kalman.template as<U>(),
            sys.template as<U>(),
            success,
        };
    }
};

/**
 * @brief Discrete offset-free LQG (LQR + bias-augmented Kalman)
 *
 * LQR is designed on @p sys (Q on x only). The filter is designed on the
 * input-disturbance augment used by design::mpc (OffsetFreeMPC). Fails if
 * the Pannocchia–Rawlings rank test fails or either DARE fails.
 *
 * Square plants (NU = NY) are required so the DC target (x_t, u_t) is unique.
 *
 * @param sys             Discrete-time plant (Ts > 0)
 * @param Q_lqr           State cost (NX × NX)
 * @param R_lqr           Input cost (NU × NU)
 * @param Q_process       Process-noise covariance on x (NX × NX)
 * @param R_measurement   Measurement-noise covariance (NY × NY)
 * @param Q_disturbance   Random-walk covariance on w (NU × NU, default I)
 * @param N               LQR cross term (default 0)
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr OFLQGResult<NX, NU, NY, T> discrete_oflqg(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NX, NX, T>&                 Q_lqr,
    const Matrix<NU, NU, T>&                 R_lqr,
    const Matrix<NX, NX, T>&                 Q_process,
    const Matrix<NY, NY, T>&                 R_measurement,
    const Matrix<NU, NU, T>&                 Q_disturbance = Matrix<NU, NU, T>::identity(),
    const Matrix<NX, NU, T>&                 N = Matrix<NX, NU, T>{}
) {
    static_assert(NU == NY, "OFLQG DC target needs a square plant (NU = NY).");
    constexpr size_t           NXK = NX + NU;
    OFLQGResult<NX, NU, NY, T> art{};
    art.sys.A = sys.A;
    art.sys.B = sys.B;
    art.sys.C = sys.C;
    art.sys.D = sys.D;
    art.sys.Ts = sys.Ts;

    Matrix<NX + NY, NX + NU, T> rank_test{};
    rank_test.template block<NX, NX>(0, 0) = sys.A - Matrix<NX, NX, T>::identity();
    rank_test.template block<NX, NU>(0, NX) = sys.B;
    rank_test.template block<NY, NX>(NX, 0) = sys.C;
    rank_test.template block<NY, NU>(NX, NX) = sys.D;
    if (mat::rank(rank_test) != NX + NU) {
        return art;
    }

    art.lqr = discrete_lqr(sys.A, sys.B, Q_lqr, R_lqr, N);

    StateSpace<NXK, NU, NY, T, NXK, NY> sys_aug{};
    sys_aug.A.template block<NX, NX>(0, 0) = sys.A;
    sys_aug.A.template block<NX, NU>(0, NX) = sys.B;
    sys_aug.A.template block<NU, NU>(NX, NX) = Matrix<NU, NU, T>::identity();
    sys_aug.B.template block<NX, NU>(0, 0) = sys.B;
    sys_aug.C.template block<NY, NX>(0, 0) = sys.C;
    sys_aug.Ts = (sys.Ts > T{0}) ? sys.Ts : T{1};

    Matrix<NXK, NXK, T> Q_aug{};
    Q_aug.template block<NX, NX>(0, 0) = Q_process;
    Q_aug.template block<NU, NU>(NX, NX) = Q_disturbance;
    art.kalman = kalman(sys_aug, Q_aug, R_measurement);
    art.kalman.sys.Ts = sys.Ts;
    art.success = art.lqr.success && art.kalman.success;
    return art;
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Offset-free LQG: u = u_t − K(x̂ − x_t) with ŵ in the DC target
 *
 * control(r, y) is a self-contained tick (OutputFeedbackController). Prefer
 * feedback/commit when the applied input can saturate.
 */
template<size_t NX, size_t NU, size_t NY, typename T = float>
struct OFLQG {
    static constexpr bool   output_feedback_tick = true;
    static constexpr size_t NXK = NX + NU;
    using Result = design::OFLQGResult<NX, NU, NY, T>;

    StateFeedback<NX, NU, T>                         lqr{};
    SteadyStateKalmanFilter<NXK, NU, NY, T, NXK, NY> kf{};
    StateSpace<NX, NU, NY, T>                        sys{};
    ColVec<NU, T>                                    u_prev{};

    constexpr OFLQG() = default;

    template<typename U>
    constexpr OFLQG(const design::OFLQGResult<NX, NU, NY, U>& result) // NOLINT
        : lqr(result.lqr), kf(result.kalman), sys(result.sys.template as<T>()) {}

    constexpr void predict(const ColVec<NU, T>& u = ColVec<NU, T>{}) { kf.predict(u); }
    constexpr bool update(const ColVec<NY, T>& y, const ColVec<NU, T>& u = ColVec<NU, T>{}) {
        return kf.update(y, u);
    }

    /// Estimated plant state (first NX of the augmented filter)
    [[nodiscard]] constexpr ColVec<NX, T> state_estimate() const {
        ColVec<NX, T> x{};
        for (size_t i = 0; i < NX; ++i) {
            x(i) = kf.state()(i);
        }
        return x;
    }

    /// Estimated input disturbance ŵ
    [[nodiscard]] constexpr ColVec<NU, T> disturbance_estimate() const {
        ColVec<NU, T> w{};
        for (size_t j = 0; j < NU; ++j) {
            w(j) = kf.state()(NX + j);
        }
        return w;
    }

    /**
     * @brief Discrete DC pair so y → r under the estimated load ŵ
     *
     * Solves @f$ (A-I)x_t + B u_t = -B w @f$, @f$ C x_t + D u_t = r @f$.
     */
    [[nodiscard]] constexpr bool equilibrium(
        const ColVec<NY, T>& r,
        const ColVec<NU, T>& w,
        ColVec<NX, T>&       xt,
        ColVec<NU, T>&       ut
    ) const {
        constexpr size_t N = NX + NY;
        Matrix<N, N, T>  M{};
        ColVec<N, T>     b{};
        M.template block<NX, NX>(0, 0) = sys.A - Matrix<NX, NX, T>::identity();
        M.template block<NX, NU>(0, NX) = sys.B;
        M.template block<NY, NX>(NX, 0) = sys.C;
        M.template block<NY, NU>(NX, NX) = sys.D;
        const ColVec<NX, T> rhs_x = -(sys.B * w);
        for (size_t i = 0; i < NX; ++i) {
            b(i) = rhs_x(i);
        }
        for (size_t i = 0; i < NY; ++i) {
            b(NX + i) = r(i);
        }
        const auto sol = mat::solve(M, b);
        if (!sol) {
            return false;
        }
        xt.template segment<NX>(0) = (*sol).template block<NX, 1>(0, 0);
        ut.template segment<NU>(0) = (*sol).template block<NU, 1>(NX, 0);
        return true;
    }

    /// Regulator at the current estimate: u = u_t − K(x̂ − x_t)
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NY, T>& r) const {
        ColVec<NX, T> xt{};
        ColVec<NU, T> ut{};
        if (!equilibrium(r, disturbance_estimate(), xt, ut)) {
            // Same matrix the rank test factored. Hold the last applied input
            // rather than dropping r and regulating with −K x̂.
            return u_prev;
        }
        return ut + lqr.control(xt, state_estimate());
    }

    [[nodiscard]] constexpr ColVec<NU, T> feedback(const ColVec<NY, T>& r, const ColVec<NY, T>& y) {
        kf.update(y, u_prev);
        return control(r);
    }

    constexpr void commit(const ColVec<NU, T>& u_applied) {
        kf.predict(u_applied);
        u_prev = u_applied;
    }

    /// Self-contained tick from the output measurement (and r)
    constexpr ColVec<NU, T> control(const ColVec<NY, T>& r, const ColVec<NY, T>& y) {
        const auto u = feedback(r, y);
        commit(u);
        return u;
    }

    [[nodiscard]] constexpr ColVec<NU, T> step(const ColVec<NY, T>& r, const ColVec<NY, T>& y) {
        return control(r, y);
    }

    /// Clear the applied-input memory (the estimator state is left untouched).
    constexpr void reset() { u_prev = ColVec<NU, T>{}; }
};

template<size_t NX, size_t NU, size_t NY, typename T>
OFLQG(const design::OFLQGResult<NX, NU, NY, T>&) -> OFLQG<NX, NU, NY, T>;

} // namespace damp
