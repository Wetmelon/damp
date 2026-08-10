// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file kalman.hpp
 * @brief Steady-state Kalman filter design and runtime estimator
 *
 * @code
 * using namespace damp;
 * constexpr auto res = design::kalman(sys, Q, R);
 * static_assert(res.success);
 * KalmanFilter filt(res);
 * filt.predict(u);
 * filt.update(y, u);
 * auto xhat = filt.state();
 * @endcode
 */

#include <cstddef>

#include "damp/design/riccati.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/solve.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace design {

/**
 * @struct KalmanResult
 * @brief Steady-state Kalman filter design result
 *
 * Contains the optimal estimator gain L and steady-state error covariance P
 * that minimize E[‖x − x̂‖²] for the given noise statistics.
 * Use .as<float>() to convert for embedded deployment.
 *
 * @see "Optimal State Estimation" (Simon, 2006), §5.3
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
struct KalmanResult {
    StateSpace<NX, NU, NY, T, NW, NV> sys{};          ///< System model
    Matrix<NW, NW, T>                 Q{};            ///< Process noise covariance
    Matrix<NV, NV, T>                 R{};            ///< Measurement noise covariance
    Matrix<NX, NY, T>                 L{};            ///< Kalman gain (steady-state)
    Matrix<NX, NX, T>                 P{};            ///< Error covariance (steady-state)
    bool                              success{false}; ///< Indicates filter design success

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return KalmanResult<NX, NU, NY, U, NW, NV>{
            sys.template as<U>(),
            Q.template as<U>(),
            R.template as<U>(),
            L.template as<U>(),
            P.template as<U>(),
            success
        };
    }
};

/**
 * @brief Steady-state Kalman filter design
 *
 * Computes the optimal steady-state estimator gain L for the discrete system:
 *
 *     x[k+1] = Ax[k] + Bu[k] + Gw[k]
 *     y[k]   = Cx[k] + Du[k] + Hv[k]
 *
 * where w ~ N(0, Q) and v ~ N(0, R). The gain minimizes the steady-state
 * error covariance P = E[(x − x̂)(x − x̂)ᵀ] by solving the filter DARE:
 *
 *     P = APAᵀ + GQGᵀ − APCᵀ(CPCᵀ + HRHᵀ)⁻¹CPAᵀ
 *
 * then computing L = PCᵀ(CPCᵀ + HRHᵀ)⁻¹.
 *
 * Taxonomy. This is the design product (L, P, plant). The underlying
 * Riccati solve is the control-form design::dare on the dual:
 *
 *     P = dare(Aᵀ, Cᵀ, GQGᵀ, HRHᵀ)
 *
 * Use design::dare with that dual only when P alone is required; use this
 * function when the estimator gain L (or LQG’s estimator half) is the goal.
 * Requires a discrete plant (Ts > 0).
 *
 * @note Compare with MATLAB®'s kalman(sys, Q, R) or dlqe(A, G, C, Q, R).
 *
 * @see design::discrete_lqe — dual-of-LQR spelling (same implementation)
 * @see design::dare — control-form DARE solver; filter form via dual (Aᵀ, Cᵀ, …)
 * @see design::discrete_lqg — pairs this L with LQR’s K
 * @see "Optimal State Estimation" (Simon, 2006), §5.3–5.4
 *
 * @param sys  Discrete-time state-space system (Ts > 0)
 * @param Q    Process noise covariance (NW × NW, positive semidefinite)
 * @param R    Measurement noise covariance (NV × NV, positive definite)
 * @return KalmanResult containing steady-state gain L and error covariance P
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr KalmanResult<NX, NU, NY, T, NW, NV> kalman(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NW, NW, T>&                 Q,
    const Matrix<NV, NV, T>&                 R
) {
    KalmanResult<NX, NU, NY, T, NW, NV> result{sys, Q, R};

    // Filter DARE is discrete-time only (same gate as design::mhe).
    // Continuous plants would produce a wrong L with success=true if we proceeded.
    if (!(sys.Ts > T{0})) {
        return result;
    }

    // Compute effective noise covariances accounting for G and H. The result
    // dimensions (NX×NX and NY×NY) follow from the matrix ops, so let them deduce.
    const auto Q_eff = sys.G * Q * sys.G.t();
    const auto R_eff = sys.H * R * sys.H.t();

    // Fast path: R_eff ≈ 0 with square, invertible C → L = C⁻¹ (solve C L = I).
    // The inverse itself is the gain deliverable; use LU solve (no bare .inverse()).
    const T r_eps = std::is_same_v<T, float> ? static_cast<T>(1e-6) : static_cast<T>(1e-10);
    if (R_eff.norm() < r_eps) {
        if constexpr (NY == NX) {
            const auto L_opt = mat::lu_solve(sys.C, Matrix<NX, NX, T>::identity());
            if (L_opt) {
                result.P = static_cast<T>(0.5) * (Q_eff + Q_eff.transpose()); // keep P symmetric
                result.L = L_opt.value();
                result.success = true;
                return result;
            }
        }
    }

    // Solve filter DARE: P = A*P*A' + Q_eff - A*P*C'*(C*P*C' + R_eff)^{-1}*C*P*A'
    // dare() handles R ≥ 0 (falls back to RDE iteration when R is singular)
    const auto dare_opt = dare(sys.A.transpose(), sys.C.transpose(), Q_eff, R_eff);
    if (!dare_opt) {
        return result;
    }
    // Scrub DARE round-off asymmetry before gain solve (S = CPCᵀ + R must stay SPD).
    {
        const auto& P_raw = dare_opt.value();
        result.P = static_cast<T>(0.5) * (P_raw + P_raw.transpose());
    }

    // Compute Kalman gain: L = PCᵀS⁻¹ → solve S Lᵀ = CP (Cholesky, LU if S semi-def)
    const Matrix<NY, NY, T> S = sys.C * result.P * sys.C.t() + R_eff;
    const Matrix<NY, NX, T> CP = sys.C * result.P;
    auto                    L_opt = mat::cholesky_solve(S, CP);
    if (!L_opt) {
        L_opt = mat::lu_solve(S, CP);
    }
    if (!L_opt) {
        return result; // singular innovation covariance S
    }
    result.L = L_opt.value().transpose();

    result.success = true;
    return result;
}

/**
 * @brief Discrete linear-quadratic estimator (LQE) — dual-of-LQR spelling of @ref kalman
 *
 * Same steady-state Kalman design as kalman: filter DARE on the dual of
 * @ref discrete_lqr. Prefer @ref kalman in new DiD code; use this name when pairing
 * with @ref discrete_lqr for textbook LQG duality.
 *
 * @note Compare with MATLAB®'s [L,P,Z,E] = dlqe(A, G, C, Q, R) (sys form stores G/H on the plant).
 * @see kalman, discrete_lqr
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr KalmanResult<NX, NU, NY, T, NW, NV> discrete_lqe(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NW, NW, T>&                 Q,
    const Matrix<NV, NV, T>&                 R
) {
    return kalman(sys, Q, R);
}

/**
 * @brief Discrete LQE from matrices (MATLAB® @c dlqe argument shape)
 *
 * Builds a discrete plant with process noise input @p G, identity measurement
 * noise map @f$ H = I @f$, zero @p B/@p D (unused for the gain), and @p Ts = 1
 * (marker that the model is discrete), then calls @ref kalman.
 *
 * @warning The returned @c KalmanResult::sys is a gain-design carrier only
 *          (B/D zeroed, Ts marked discrete). Deploy KalmanFilter from a real
 *          plant + this L/P, or use the discrete_lqe(sys, Q, R) overload.
 *
 * @param A  State matrix (NX × NX)
 * @param G  Process-noise input matrix (NX × NW)
 * @param C  Output matrix (NY × NX)
 * @param Q  Process noise covariance (NW × NW)
 * @param R  Measurement noise covariance (NY × NY)
 * @return KalmanResult with steady-state L and P (not closed-loop poles like MATLAB E)
 *
 * @note Compare with MATLAB®'s [L,P,Z,E] = dlqe(A, G, C, Q, R) — Damp returns L and P
 *       on @ref KalmanResult (not Z/E).
 * @see discrete_lqe(sys, Q, R), kalman, discrete_lqr
 */
template<size_t NX, size_t NW, size_t NY, typename T = double>
[[nodiscard]] constexpr KalmanResult<NX, 1, NY, T, NW, NY> discrete_lqe(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NW, T>& G,
    const Matrix<NY, NX, T>& C,
    const Matrix<NW, NW, T>& Q,
    const Matrix<NY, NY, T>& R
) {
    StateSpace<NX, 1, NY, T, NW, NY> sys{};
    sys.A = A;
    sys.G = G;
    sys.C = C;
    sys.H = Matrix<NY, NY, T>::identity();
    sys.Ts = T{1}; // discrete marker; matrices are already discrete
    return kalman(sys, Q, R);
}

} // namespace design

/**
 * @brief Runtime Kalman filter for embedded systems
 *
 * Implements the standard predict/update cycle:
 *
 *     Predict:  x̂[k|k−1] = Ax̂[k−1] + Bu[k−1]
 *               P[k|k−1]  = AP[k−1]Aᵀ + GQGᵀ
 *
 *     Update:   K = P[k|k−1]Cᵀ(CP[k|k−1]Cᵀ + HRHᵀ)⁻¹
 *               x̂[k|k] = x̂[k|k−1] + K(y[k] − Cx̂[k|k−1])
 *               P[k|k]  = (I − KC)P(I − KC)ᵀ + KHRHᵀKᵀ   [Joseph form]
 *
 * The Joseph form update is used for numerical stability (maintains P symmetry
 * and positive-definiteness even with finite-precision arithmetic).
 *
 * @see "Optimal State Estimation" (Simon, 2006), §5.1, §7.6.2 (Joseph form)
 *
 * @tparam NX Number of states
 * @tparam NU Number of inputs
 * @tparam NY Number of outputs
 * @tparam NW Number of process noise inputs
 * @tparam NV Number of measurement noise inputs
 * @tparam T  Scalar type (default: float — embedded runtime)
 */
template<size_t NX, size_t NU, size_t NY, typename T = float, size_t NW = 0, size_t NV = 0>
struct KalmanFilter {
    constexpr KalmanFilter() = default;

    constexpr KalmanFilter(
        const StateSpace<NX, NU, NY, T, NW, NV>& sys_,
        const Matrix<NW, NW, T>&                 Q_,
        const Matrix<NV, NV, T>&                 R_,
        const ColVec<NX, T>&                     x0 = ColVec<NX, T>{},
        const Matrix<NX, NX, T>&                 P0 = Matrix<NX, NX, T>::identity()
    ) : sys(sys_), x(x0), P(P0), Q(Q_), R(R_) {}

    // Construct directly from a design result (any scalar). The result's
    // `sys` already carries A/B/C/D/G/H. Steady-state P seeds P0 (identity if
    // the design did not converge).
    template<typename U>
    constexpr KalmanFilter(const design::KalmanResult<NX, NU, NY, U, NW, NV>& result) // NOLINT
        : sys(result.sys.template as<T>()),
          x(ColVec<NX, T>{}),
          P(result.success ? result.P.template as<T>() : Matrix<NX, NX, T>::identity()),
          Q(result.Q.template as<T>()),
          R(result.R.template as<T>()) {}

    // Type conversion constructor
    template<typename U>
    constexpr KalmanFilter(const KalmanFilter<NX, NU, NY, U, NW, NV>& other) // NOLINT
        : sys(other.model()),
          x(other.state()),
          P(other.covariance()),
          Q(other.process_noise_covariance()),
          R(other.measurement_noise_covariance()),
          innov(other.innovation()) {}

    // Predict: x[k+1|k] = A*x[k|k] + B*u[k], P[k+1|k] = APAᵀ + GQGᵀ
    constexpr void predict(const ColVec<NU, T>& u = ColVec<NU, T>{}) {
        x = sys.A * x + sys.B * u;
        P = quadratic_form(sys.A, P) + quadratic_form(sys.G, Q);
    }

    // Measurement update: returns false if innovation covariance is singular
    constexpr bool update(const ColVec<NY, T>& y, const ColVec<NU, T>& u = ColVec<NU, T>{}) {
        const auto y_pred = sys.C * x + sys.D * u;
        innov = y - y_pred;

        // S = CPCᵀ + HRHᵀ
        const auto S = quadratic_form(sys.C, P) + quadratic_form(sys.H, R);

        // K = PCᵀS⁻¹ → solve S Kᵀ = C P (Cholesky; S must be SPD). Fail closed if singular.
        auto K_opt = mat::cholesky_solve(S, sys.C * P);
        if (!K_opt) {
            K_opt = mat::lu_solve(S, sys.C * P);
        }
        if (!K_opt) {
            return false;
        }

        const Matrix<NX, NY, T> K = K_opt.value().transpose();
        x = x + K * innov;

        // Joseph form via quadratic_form: exact P symmetry under finite precision
        //     P = (I − KC) P (I − KC)ᵀ + (KH) R (KH)ᵀ
        const auto I_KC = Matrix<NX, NX, T>::identity() - (K * sys.C);
        const auto KH = K * sys.H;
        P = quadratic_form(I_KC, P) + quadratic_form(KH, R);

        return true;
    }

    // Measurement update with an explicit measurement model (C, D, R), for a sensor
    // whose output equation differs from the filter's nominal sys.C/sys.D — i.e.
    // multi-sensor fusion where heterogeneous measurements (e.g. an angle encoder and
    // a load accelerometer) share one filter. Same dimensions as the nominal update.
    // Returns false if the innovation covariance is singular.
    constexpr bool update(
        const ColVec<NY, T>&     y,
        const Matrix<NY, NX, T>& C_meas,
        const Matrix<NY, NU, T>& D_meas,
        const Matrix<NV, NV, T>& R_meas,
        const ColVec<NU, T>&     u = ColVec<NU, T>{}
    ) {
        const auto y_pred = C_meas * x + D_meas * u;
        innov = y - y_pred;

        const auto S = quadratic_form(C_meas, P) + quadratic_form(sys.H, R_meas);

        auto K_opt = mat::cholesky_solve(S, C_meas * P);
        if (!K_opt) {
            K_opt = mat::lu_solve(S, C_meas * P);
        }
        if (!K_opt) {
            return false;
        }

        const Matrix<NX, NY, T> K = K_opt.value().transpose();
        x = x + K * innov;

        const auto I_KC = Matrix<NX, NX, T>::identity() - K * C_meas;
        const auto KH = K * sys.H;
        P = quadratic_form(I_KC, P) + quadratic_form(KH, R_meas);

        return true;
    }

    /**
     * @brief Fused per-tick estimate: predict with u, correct with y, return the state.
     *
     * The StateEstimator concept shape — one call = one complete filter tick,
     * so generic code cannot skip the time update. The explicit predict()/
     * update() split stays for multirate use (predict every tick, update when
     * a measurement arrives). A singular innovation covariance skips the
     * correction (the prediction stands), matching update()'s false return.
     *
     * @param y Measurement taken now
     * @param u Input applied over the previous interval
     * @return The updated state estimate
     */
    constexpr const ColVec<NX, T>& estimate(const ColVec<NY, T>& y, const ColVec<NU, T>& u = ColVec<NU, T>{}) {
        predict(u);
        (void)update(y, u);
        return x;
    }

    // Accessors
    [[nodiscard]] constexpr const auto& model() const { return sys; }
    [[nodiscard]] constexpr const auto& innovation() const { return innov; }
    [[nodiscard]] constexpr const auto& state() const { return x; }
    [[nodiscard]] constexpr const auto& covariance() const { return P; }
    [[nodiscard]] constexpr const auto& process_noise_covariance() const { return Q; }
    [[nodiscard]] constexpr const auto& measurement_noise_covariance() const { return R; }

    // Mutators. set_state lets the caller enforce physical constraints on the
    // estimate after an update() — clamp a state to its actuator/sensor range,
    // wrap an angle, or zero a state that has gone non-physical — before the next
    // predict()/update(). Clamping the estimate keeps the filter's state vector
    // meaningful; the next predict() then propagates from the constrained value.
    constexpr void set_state(const ColVec<NX, T>& x_new) { x = x_new; }
    constexpr void set_state(size_t i, T value) { x[i] = value; }
    constexpr void set_covariance(const Matrix<NX, NX, T>& P_new) { P = P_new; }

    // Re-initialize the estimate and clear the innovation (model and noise
    // covariances are kept). Restores the filter to a known starting belief.
    constexpr void reset(
        const ColVec<NX, T>&     x0 = ColVec<NX, T>{},
        const Matrix<NX, NX, T>& P0 = Matrix<NX, NX, T>::identity()
    ) {
        x = x0;
        P = P0;
        innov = ColVec<NY, T>{};
    }

private:
    StateSpace<NX, NU, NY, T, NW, NV> sys{}; //< System Model (discrete)

    ColVec<NX, T>     x{};     //< System state (predicted)
    Matrix<NX, NX, T> P{};     //< State estimate error covariance
    Matrix<NW, NW, T> Q{};     //< Process noise covariance
    Matrix<NV, NV, T> R{};     //< Measurement noise covariance
    ColVec<NY, T>     innov{}; //< Residual between measurements and estimated measurements
};

/**
 * @brief Steady-state (fixed-gain) Kalman estimator for LQG-class designs
 *
 * Uses the designed gain L from design::kalman and never updates P online:
 *
 *     Predict:  x̂ ← A x̂ + B u
 *     Update:   x̂ ← x̂ + L (y − C x̂ − D u)
 *
 * Same structure as a Luenberger observer; L is the steady-state Kalman gain
 * from the filter DARE (classical LQG separation). Prefer this for LQG/LQGI
 * deploy. Use KalmanFilter when the gain must stay time-varying (Joseph
 * covariance recursion, multi-sensor R/C changes).
 *
 * @see design::kalman
 * @see KalmanFilter — full recursive filter
 * @see "Optimal Control" (Anderson & Moore, 1990), §8
 */
template<size_t NX, size_t NU, size_t NY, typename T = float, size_t NW = 0, size_t NV = 0>
struct SteadyStateKalmanFilter {
    constexpr SteadyStateKalmanFilter() = default;

    constexpr SteadyStateKalmanFilter(
        const StateSpace<NX, NU, NY, T, NW, NV>& sys_,
        const Matrix<NX, NY, T>&                 L_,
        const ColVec<NX, T>&                     x0 = ColVec<NX, T>{},
        const Matrix<NX, NX, T>&                 P_ss_ = Matrix<NX, NX, T>::identity()
    )
        : sys(sys_), L(L_), P_ss(P_ss_), x(x0) {}

    /// From design (any scalar): plant + steady-state L and P∞ (zero gain if design failed).
    template<typename U>
    constexpr SteadyStateKalmanFilter(const design::KalmanResult<NX, NU, NY, U, NW, NV>& result)
        : sys(result.sys.template as<T>()),
          L(result.success ? result.L.template as<T>() : Matrix<NX, NY, T>{}),
          P_ss(result.success ? result.P.template as<T>() : Matrix<NX, NX, T>::identity()),
          x(ColVec<NX, T>{}) {}

    template<typename U>
    constexpr SteadyStateKalmanFilter(const SteadyStateKalmanFilter<NX, NU, NY, U, NW, NV>& other)
        : sys(other.model()),
          L(other.gain().template as<T>()),
          P_ss(other.covariance().template as<T>()),
          x(other.state().template as<T>()),
          innov(other.innovation().template as<T>()) {}

    /// Time update: x̂ ← A x̂ + B u (no covariance recursion).
    constexpr void predict(const ColVec<NU, T>& u = ColVec<NU, T>{}) { x = sys.A * x + sys.B * u; }

    /// Measurement update with fixed L. Always succeeds (no Cholesky).
    constexpr bool update(const ColVec<NY, T>& y, const ColVec<NU, T>& u = ColVec<NU, T>{}) {
        innov = y - (sys.C * x + sys.D * u);
        x = x + L * innov;
        return true;
    }

    /// Fused tick: predict then update; returns the state (StateEstimator shape).
    constexpr const ColVec<NX, T>& estimate(const ColVec<NY, T>& y, const ColVec<NU, T>& u = ColVec<NU, T>{}) {
        predict(u);
        (void)update(y, u);
        return x;
    }

    [[nodiscard]] constexpr const auto& model() const { return sys; }
    [[nodiscard]] constexpr const auto& gain() const { return L; }
    /// Steady-state design covariance P∞ (not updated online).
    [[nodiscard]] constexpr const auto& covariance() const { return P_ss; }
    [[nodiscard]] constexpr const auto& innovation() const { return innov; }
    [[nodiscard]] constexpr const auto& state() const { return x; }

    constexpr void set_state(const ColVec<NX, T>& x_new) { x = x_new; }
    constexpr void set_state(size_t i, T value) { x[i] = value; }
    constexpr void set_gain(const Matrix<NX, NY, T>& L_new) { L = L_new; }

    constexpr void reset(const ColVec<NX, T>& x0 = ColVec<NX, T>{}) {
        x = x0;
        innov = ColVec<NY, T>{};
    }

private:
    StateSpace<NX, NU, NY, T, NW, NV> sys{};
    Matrix<NX, NY, T>                 L{};
    Matrix<NX, NX, T>                 P_ss{}; ///< Design P∞ only (not recursed)
    ColVec<NX, T>                     x{};
    ColVec<NY, T>                     innov{};
};

} // namespace damp
