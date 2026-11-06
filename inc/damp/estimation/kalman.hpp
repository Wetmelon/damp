// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file kalman.hpp
 * @brief Steady-state Kalman filter design (design::kalman) and the recursive KalmanFilter
 *
 * design::kalman solves the filter DARE for the steady-state covariance P∞ and gains: the
 * innovation gain M (MATLAB®'s Mx) and the predictor gain L = A M (MATLAB®'s L). For a fixed
 * (steady-state) gain, run the design on Luenberger (luenberger.hpp), the fixed-gain observer the
 * LQ and MPC laws hold: `Luenberger obs{design::kalman(sys, noise)}`. KalmanFilter propagates its
 * covariance each tick, for a reading whose R or C changes at run time.
 *
 * @code
 * using namespace damp;
 * // discrete plant: per-step process covariance, per-reading measurement variance
 * constexpr auto res = design::kalman(sys, {.process = Qd, .measurement = sigma2});
 * // continuous plant: process-noise density (Van Loan'd), the same per-reading variance
 * constexpr auto res_c = design::kalman(sys_c, Ts, {.process = Qc, .measurement = sigma2});
 * // a filter on [x; d]: noise.disturbance is the random walk of d
 * // (the disturbance spec sets the count, so the noise is a plain braced list)
 * constexpr auto res_d = design::kalman(sys, {.process = Qd, .measurement = sigma2, .disturbance = Qdd}, design::InputDisturbance{});
 * static constinit Luenberger obs{res};     // the fixed gain M (luenberger.hpp); a failed design does not compile
 * static constinit KalmanFilter filt{res};  // the recursive filter, P∞ its initial P
 * // each tick: auto xhat = obs.step(y, u);  // or obs.predict(u); obs.update(y, u); for a split tick
 * @endcode
 */


#include "damp/backend.hpp"
#include "damp/controllers/law_detail.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/solve.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace design {

/**
 * @brief The estimator's noise: process noise on the plant's noise inputs (through G),
 *        measurement noise on the sensors (through H)
 *
 * The meaning of @c process follows the plant's time domain, as LQRWeights does. On a
 * continuous plant (the `kalman(sys, Ts, noise)` form) it is a power spectral density, the
 * intensity of white noise in ẋ = Ax + Bu + Gw; the design integrates it over one sample
 * (Van Loan, @ref discretize_process_noise), so the same density gives the same filter
 * behaviour at any rate. On a discrete plant it is the covariance added per step.
 *
 * @c measurement is always the variance of each reading as the firmware gets it, with no Ts
 * scaling on either plant: log the sensor at rest and take the variance, or from a datasheet
 * σ² = ASD² × B with B the noise bandwidth of the configured sensor filter. MATLAB®'s kalmd
 * converts a density with R/Ts, which assumes an integrate-and-dump sampler whose noise
 * bandwidth tracks Ts; damp does not (matlab::kalmd keeps that convention).
 *
 * @code
 * design::kalman(continuous_sys, Ts, {.process = Qc, .measurement = sigma2});
 * design::kalman(discrete_sys, {.process = Qd, .measurement = sigma2});
 * @endcode
 *
 * With a disturbance model (a filter on [x; d]), @c disturbance is the random walk of d: how
 * fast the disturbance may move, in the plant's time domain like @c process (a density on a
 * continuous plant, a covariance per step on a discrete one). Without one it has no entries.
 *
 * @tparam NW Process-noise inputs (columns of G); NV Measurement-noise inputs (columns of H).
 * @tparam ND Disturbance states (0: no disturbance model).
 */
template<size_t NW, size_t NV, typename T = double, size_t ND = 0>
struct KalmanNoise {
    /**
     * @brief Process noise on the NW inputs of G (NW × NW, positive semidefinite)
     *
     * On a continuous plant, a power spectral density: the intensity of the white noise in
     * ẋ = Ax + Bu + Gw. The design integrates it over one sample (Van Loan), so the same
     * density means the same at any rate. On a discrete plant, the covariance added per step.
     */
    Matrix<NW, NW, T> process = Matrix<NW, NW, T>::identity();
    /**
     * @brief Variance of each reading as the firmware gets it (NV × NV, positive definite)
     *
     * The same on either plant, with no Ts scaling. Log the sensor at rest and take the
     * variance, or from a datasheet σ² = ASD² × B with B the noise bandwidth of the configured
     * sensor filter. MATLAB®'s kalmd R/Ts assumes an integrate-and-dump sampler; damp does not.
     */
    Matrix<NV, NV, T> measurement = Matrix<NV, NV, T>::identity();
    /**
     * @brief Random walk of the disturbance d (ND × ND, positive definite; a d with no walk is
     *        refused with DesignStatus::bad_weights), used only with a disturbance model: a
     *        density on a continuous plant, a covariance per step on a discrete one.
     */
    Matrix<ND, ND, T> disturbance = Matrix<ND, ND, T>::identity();
};

/** @brief Sizes from the matrices: `KalmanNoise{Qn, Rn}`. */
template<size_t NW, size_t NV, typename T>
KalmanNoise(Matrix<NW, NW, T>, Matrix<NV, NV, T>) -> KalmanNoise<NW, NV, T>;

/** @brief Sizes from the matrices: `KalmanNoise{Qn, Rn, Qd}`, with the disturbance's random walk. */
template<size_t NW, size_t NV, typename T, size_t ND>
KalmanNoise(Matrix<NW, NW, T>, Matrix<NV, NV, T>, Matrix<ND, ND, T>) -> KalmanNoise<NW, NV, T, ND>;

/**
 * @struct KalmanResult
 * @brief Steady-state Kalman filter design result
 *
 * The steady-state gains and error covariance that minimize E[‖x − x̂‖²] for the given noise
 * statistics, under MATLAB®'s names:
 *
 * | Member | What it is | MATLAB® |
 * | --- | --- | --- |
 * | M | innovation gain P Cᵀ (C P Cᵀ + H R Hᵀ)⁻¹: x̂[k\|k] = x̂[k\|k−1] + M (y − C x̂[k\|k−1] − D u) | kalman's Mx, dlqe's M |
 * | L | predictor gain A M: x̂[k+1\|k] = A x̂[k\|k−1] + B u + L (y − C x̂[k\|k−1] − D u) | kalman's L |
 * | P | a priori error covariance P[k\|k−1], the filter DARE's solution | kalman's P |
 *
 * Luenberger runs M in the current form (`Luenberger obs{res}`); KalmanFilter starts from P.
 * Use .as<float>() to convert for embedded deployment.
 *
 * With a disturbance (design::kalman(sys[, Ts], noise, disturbance)), @c sys is the filter's model:
 * the plant, then the ND disturbance states; @c M, @c L and @c P are the steady-state gains and
 * covariance on that whole model.
 *
 * @tparam ND Disturbance states, last in the filter's state (0: none)
 *
 * @see "Optimal State Estimation" (Simon, 2006), §5.3
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0, size_t ND = 0>
struct [[nodiscard]] KalmanResult {
    StateSpace<NX, NU, NY, T, NW, NV> sys{};                               ///< System model
    Matrix<NW, NW, T>                 Q{};                                 ///< Process noise covariance
    Matrix<NV, NV, T>                 R{};                                 ///< Measurement noise covariance
    Matrix<NX, NY, T>                 L{};                                 ///< Predictor gain A M (MATLAB®'s L): x̂[k+1|k] = A x̂ + B u + L (y − C x̂ − D u)
    Matrix<NX, NY, T>                 M{};                                 ///< Innovation gain P Cᵀ (C P Cᵀ + H R Hᵀ)⁻¹ (MATLAB®'s Mx): the current-form gain
    Matrix<NX, NX, T>                 P{};                                 ///< A priori error covariance (steady-state)
    bool                              success{false};                      ///< status == DesignStatus::ok
    DesignStatus                      status{DesignStatus::not_converged}; ///< ok, or why the design failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return KalmanResult<NX, NU, NY, U, NW, NV, ND>{
            sys.template as<U>(),
            Q.template as<U>(),
            R.template as<U>(),
            L.template as<U>(),
            M.template as<U>(),
            P.template as<U>(),
            success,
            status,
        };
    }
};

namespace detail {

/** @brief True when every entry of the model (A, B, C, D, G, H) is finite. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr bool finite_model(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    const auto finite = [](T v) { return damp::isfinite(v); };
    return mat::all_of(sys.A, finite) && mat::all_of(sys.B, finite) && mat::all_of(sys.C, finite) && mat::all_of(sys.D, finite)
        && mat::all_of(sys.G, finite) && mat::all_of(sys.H, finite);
}

/** @brief True when every entry of the filter's model (A, B, C, D, G, H) and noise (Q, R) is finite. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr bool finite_filter_model(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const Matrix<NW, NW, T>& Q, const Matrix<NV, NV, T>& R) {
    const auto finite = [](T v) { return damp::isfinite(v); };
    return finite_model(sys) && mat::all_of(Q, finite) && mat::all_of(R, finite);
}

/** @brief The filter DARE on a discrete model: the plain steady-state Kalman design (see kalman). */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr KalmanResult<NX, NU, NY, T, NW, NV> kalman_plain(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const KalmanNoise<NW, NV, T>&            noise
) {
    const auto&                         Q = noise.process;
    const auto&                         R = noise.measurement;
    KalmanResult<NX, NU, NY, T, NW, NV> result{sys, Q, R};

    // Filter DARE is discrete-time only (same gate as design::mhe).
    // Continuous plants would produce a wrong gain with success=true if we proceeded.
    if (!damp::finite_positive(sys.Ts)) {
        return damp::detail::set_status(result, DesignStatus::invalid_sample_time);
    }
    // The DARE never reads B or D, so a non-finite one would ship in result.sys with success.
    if (!finite_filter_model(sys, Q, R)) {
        return damp::detail::set_status(result, DesignStatus::bad_weights);
    }

    // Compute effective noise covariances accounting for G and H. The result
    // dimensions (NX×NX and NY×NY) follow from the matrix ops, so let them deduce.
    const auto Q_eff = sys.G * Q * sys.G.t();
    const auto R_eff = sys.H * R * sys.H.t();

    // Solve filter DARE: P = A*P*A' + Q_eff - A*P*C'*(C*P*C' + R_eff)^{-1}*C*P*A'
    // No deadbeat shortcut on a small absolute R: a precise sensor is not a
    // perfect measurement, and M = C⁻¹ would ship with success == true.
    // The DARE handles R ≥ 0 (falls back to RDE iteration when R is singular).
    const auto sol = damp::detail::solve_dare(
        Matrix<NX, NX, T>(sys.A.transpose()), Matrix<NX, NY, T>(sys.C.transpose()), Matrix<NX, NX, T>(Q_eff), Matrix<NY, NY, T>(R_eff), Matrix<NX, NY, T>{}, DareMethod::Auto
    );
    if (sol.status != DesignStatus::ok) {
        return damp::detail::set_status(result, damp::detail::dual_status(sol.status));
    }
    // Scrub DARE round-off asymmetry before gain solve (S = CPCᵀ + R must stay SPD).
    result.P = static_cast<T>(0.5) * (sol.X + sol.X.transpose());

    // Innovation gain: M = PCᵀS⁻¹ → solve S Mᵀ = CP (Cholesky, LU if S semi-def)
    const Matrix<NY, NY, T> S = sys.C * result.P * sys.C.t() + R_eff;
    const auto              M_opt = mat::solve(S, Matrix<NY, NX, T>(sys.C * result.P));
    if (!M_opt) {
        return damp::detail::set_status(result, DesignStatus::singular); // singular innovation covariance S
    }
    result.M = M_opt.value().transpose();
    result.L = sys.A * result.M; // the predictor gain: x̂[k+1|k] = A x̂[k|k]
    return damp::detail::set_status(result, DesignStatus::ok);
}

} // namespace detail

/**
 * @brief Discretize continuous process noise over one sample (Van Loan)
 *
 * For @f$ \dot x = A x + G w @f$ with white @f$ w @f$ of intensity (power spectral density)
 * @f$ Q_c @f$, the noise a sampled model accumulates over one period has covariance
 * @f[
 *   Q_d = \int_0^{T_s} e^{A\tau} G Q_c G^\top e^{A^\top \tau}\, d\tau ,
 * @f]
 * which enters the discrete model through @f$ G_d = I @f$. One matrix exponential gives it
 * exactly (Van Loan, 1978):
 * @f[
 *   \exp\!\left( \begin{bmatrix} -A & G Q_c G^\top \\ 0 & A^\top \end{bmatrix} T_s \right)
 *     = \begin{bmatrix} M_{11} & M_{12} \\ 0 & M_{22} \end{bmatrix}, \quad
 *   Q_d = M_{22}^\top M_{12}.
 * @f]
 * For a short period @f$ Q_d \approx G Q_c G^\top T_s @f$; the exact form keeps the
 * correlation the dynamics build between states within a sample (a double integrator's
 * position picks up @f$ T_s^3/3 @f$ of the velocity noise). The exponential is taken over a
 * fraction of the period short enough that the product does not cancel on a fast stable
 * mode, then doubled up to @f$ T_s @f$ (detail::van_loan_integral()).
 *
 * @note This is the process-noise step of MATLAB®'s kalmd(sys, Qn, Rn, Ts).
 * @see "Computing Integrals Involving the Matrix Exponential" (Van Loan, 1978),
 *      IEEE TAC 23(3), https://doi.org/10.1109/TAC.1978.1101743
 * @see "Optimal State Estimation" (Simon, 2006), §8.1
 * @see discretize_lqr_cost — the same construction for an LQR cost integral
 *
 * @param A   Continuous state matrix (NX × NX)
 * @param G   Process-noise input matrix (NX × NW)
 * @param Qc  Process-noise intensity (NW × NW, positive semidefinite)
 * @param Ts  Sampling time [s], finite and positive
 * @return Discrete process-noise covariance Q_d (NX × NX, symmetric), or damp::nullopt on a
 *         non-finite entry or result
 */
template<size_t NX, size_t NW, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>> discretize_process_noise(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NW, T>& G,
    const Matrix<NW, NW, T>& Qc,
    T                        Ts
) {
    //! Q_d = ∫₀^Ts e^{Aτ} G Q_c Gᵀ e^{Aᵀτ} dτ: the Van Loan integral with F = Aᵀ
    const Matrix<NX, NX, T> W = G * Qc * G.t();
    return ::damp::detail::van_loan_integral(Matrix<NX, NX, T>(A.t()), W, Ts);
}

namespace detail {

/** @brief The sampled-data Kalman design on a continuous plant (see kalman(sys, Ts, noise)). */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr KalmanResult<NX, NU, NY, T, NX, NV> kalman_sampled(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        Ts,
    const KalmanNoise<NW, NV, T>&            noise
) {
    using Result = KalmanResult<NX, NU, NY, T, NX, NV>;
    if (sys.Ts != T{0} || !damp::finite_positive(Ts)) {
        return damp::detail::refused<Result>(DesignStatus::invalid_sample_time); // one discretization, at a valid Ts
    }
    const auto plant_d = ::damp::discretize(StateSpace<NX, NU, NY, T>{.A = sys.A, .B = sys.B, .C = sys.C, .D = sys.D}, Ts);
    if (!plant_d) {
        return damp::detail::refused<Result>(plant_d.status); // a non-finite A or B, or an overflow
    }
    const auto Qd = discretize_process_noise(sys.A, sys.G, noise.process, Ts);
    if (!Qd) {
        return damp::detail::refused<Result>(DesignStatus::bad_weights); // a non-finite G or noise entry
    }
    const StateSpace<NX, NU, NY, T, NX, NV> sys_d{
        .A = plant_d->A,
        .B = plant_d->B,
        .C = sys.C,
        .D = sys.D,
        .G = Matrix<NX, NX, T>::identity(),
        .H = sys.H,
        .Ts = Ts,
    };
    return kalman_plain(sys_d, KalmanNoise<NX, NV, T>{.process = *Qd, .measurement = noise.measurement});
}

/** @brief A plain design result as the result of a design with options: its ND. */
template<size_t ND, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr KalmanResult<NX, NU, NY, T, NW, NV, ND> with_disturbance_states(const KalmanResult<NX, NU, NY, T, NW, NV>& result) {
    return KalmanResult<NX, NU, NY, T, NW, NV, ND>{result.sys, result.Q, result.R, result.L, result.M, result.P, result.success, result.status};
}

/**
 * @brief The Kalman design, plain or on [x; d]: the filter's model (the plant, then the disturbance),
 *        in the plant's own time domain (@p Discrete) or sampled at @p Ts.
 */
template<bool Discrete, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t ND, typename Spec>
[[nodiscard]] constexpr auto kalman_design(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T Ts, const KalmanNoise<NW, NV, T, ND>& noise, const Spec& spec) {
    constexpr size_t NXE = NX;
    constexpr size_t NXK = NXE + ND;
    if constexpr (ND == 0) {
        static_cast<void>(spec);
        if constexpr (Discrete) {
            static_cast<void>(Ts);
            return kalman_plain(sys, noise);
        } else {
            return kalman_sampled(sys, Ts, noise);
        }
    } else {
        static_assert(ND <= NY, "Estimating ND disturbances needs at least ND independent readings (ND <= NY).");
        constexpr size_t NWK = Discrete ? (NW + ND) : NXK;
        using Result = KalmanResult<NXK, NU, NY, T, NWK, NV, ND>;
        const auto m = sampled_model(sys, Ts, Discrete);
        if (m.status != DesignStatus::ok) {
            return ::damp::detail::refused<Result>(m.status);
        }
        if (!finite_filter_model(m.law, noise.process, noise.measurement)) {
            return ::damp::detail::refused<Result>(DesignStatus::bad_weights); // before the detectability test reads C
        }
        const auto ch = disturbance_channels(spec, sys, m.law, Ts);
        if (!augment_detectable(m.law.A, m.law.C, ch.Bd, ch.Cd)) {
            return ::damp::detail::refused<Result>(DesignStatus::not_detectable);
        }
        auto                aug = disturbance_augment<NWK>(m.law, ch.Bd, ch.Cd, Discrete ? sys.Ts : Ts);
        Matrix<NWK, NWK, T> Q{};
        if constexpr (Discrete) {
            //! The process noise through G, the disturbance's random walk on d
            aug.G.template block<NXE, NW>(0, 0) = m.law.G;
            aug.G.template block<ND, ND>(NXE, NW) = Matrix<ND, ND, T>::identity();
            Q.template block<NW, NW>(0, 0) = noise.process;
            Q.template block<ND, ND>(NW, NW) = noise.disturbance;
        } else {
            //! One period of ẋ = Ax + B_d d + Gw, ḋ = w_d (Van Loan): the walk of d moves x within
            //! the period.
            constexpr size_t            NA = NX + ND;
            Matrix<NA, NA, T>           Ac{};
            Matrix<NA, NW + ND, T>      Gc{};
            Matrix<NW + ND, NW + ND, T> Qc{};
            Ac.template block<NX, NX>(0, 0) = sys.A;
            Ac.template block<NX, ND>(0, NX) = ch.Bd_c;
            Gc.template block<NX, NW>(0, 0) = sys.G;
            Gc.template block<ND, ND>(NX, NW) = Matrix<ND, ND, T>::identity();
            Qc.template block<NW, NW>(0, 0) = noise.process;
            Qc.template block<ND, ND>(NW, NW) = noise.disturbance;
            const auto Qd = discretize_process_noise(Ac, Gc, Qc, Ts);
            if (!Qd) {
                return ::damp::detail::refused<Result>(DesignStatus::bad_weights); // a non-finite noise or channel entry
            }
            aug.G = Matrix<NXK, NXK, T>::identity();
            for (size_t i = 0; i < NA; ++i) {
                for (size_t j = 0; j < NA; ++j) {
                    Q((i < NX) ? i : (i - NX + NXE), (j < NX) ? j : (j - NX + NXE)) = (*Qd)(i, j);
                }
            }
        }
        return with_disturbance_states<ND>(kalman_plain(aug, KalmanNoise<NWK, NV, T>{.process = Q, .measurement = noise.measurement}));
    }
}

/** @brief True when the noise's ND is the disturbance's and an explicit model is on this plant. */
template<typename Disturbance, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t ND>
consteval bool kalman_disturbance_ok() {
    if constexpr (disturbance_traits<Disturbance, NU, NY, NW>::states != ND) {
        return false;
    } else if constexpr (requires(const Disturbance& d) { d.Bd; }) {
        constexpr size_t NZ = disturbance_traits<Disturbance, NU, NY, NW>::listed;
        return damp::is_same_v<Disturbance, DisturbanceModel<NX, NY, ND, NZ, T>>;
    } else {
        return true;
    }
}

/** @brief The message for a disturbance that does not fit the noise or the plant. Not a success path. */
template<typename Disturbance, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t ND>
consteval void reject_kalman_disturbance() {
    constexpr bool sized = disturbance_traits<Disturbance, NU, NY, NW>::states == ND;
    if constexpr (!sized) {
        static_assert(sized, "KalmanNoise carries the random walk of the disturbance's ND states: KalmanNoise<NW, NV, T, ND>.");
    } else if constexpr (requires(const Disturbance& d) { d.Bd; }) {
        constexpr size_t NZ = disturbance_traits<Disturbance, NU, NY, NW>::listed;
        static_assert(
            damp::is_same_v<Disturbance, DisturbanceModel<NX, NY, ND, NZ, T>>,
            "The DisturbanceModel is on the plant's states and outputs, in the plant's scalar."
        );
    }
}

/**
 * @brief kalman(sys[, Ts], noise, disturbance), once the disturbance fits. A mistake stops here,
 *        before the filter DARE.
 */
template<bool Discrete, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t ND, typename Disturbance>
[[nodiscard]] constexpr auto kalman_with_disturbance(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        Ts,
    const KalmanNoise<NW, NV, T, ND>&        noise,
    const Disturbance&                       disturbance
) {
    if constexpr (kalman_disturbance_ok<Disturbance, NX, NU, NY, T, NW, ND>()) {
        return kalman_design<Discrete>(sys, Ts, noise, disturbance);
    } else {
        reject_kalman_disturbance<Disturbance, NX, NU, NY, T, NW, ND>();
        return KalmanResult<NX, NU, NY, T>{};
    }
}

} // namespace detail

/**
 * @brief Steady-state Kalman filter design
 *
 * Computes the optimal steady-state estimator gains for the discrete system:
 *
 *     x[k+1] = Ax[k] + Bu[k] + Gw[k]
 *     y[k]   = Cx[k] + Du[k] + Hv[k]
 *
 * where w ~ N(0, Q) and v ~ N(0, R), Q = noise.process and R = noise.measurement (see
 * KalmanNoise). The gain minimizes the steady-state
 * error covariance P = E[(x − x̂)(x − x̂)ᵀ] by solving the filter DARE:
 *
 *     P = APAᵀ + GQGᵀ − APCᵀ(CPCᵀ + HRHᵀ)⁻¹CPAᵀ
 *
 * then computing the innovation gain M = PCᵀ(CPCᵀ + HRHᵀ)⁻¹ (MATLAB®'s Mx, the current-form
 * gain) and the predictor gain L = A M (MATLAB®'s L).
 *
 * Taxonomy. This is the design product (L, M, P, plant). The underlying
 * Riccati solve is the control-form dare on the dual:
 *
 *     P = dare(Aᵀ, Cᵀ, GQGᵀ, HRHᵀ)
 *
 * Use dare with that dual only when P alone is required; use this
 * function when the estimator gains (or LQG’s estimator half) are the goal.
 *
 * The result runs on Luenberger, the fixed-gain observer (`Luenberger obs{res}`: M held, no
 * per-tick covariance), or on KalmanFilter, which starts from P∞ and propagates P each tick.
 *
 * Requires a discrete plant (Ts > 0): there is no continuous-time (Kalman–Bucy) design, so a
 * continuous plant is refused with DesignStatus::invalid_sample_time. For a continuous plant
 * with a noise density, pass the sample time instead: kalman(sys, Ts, noise) discretizes
 * both. A dual DARE that cannot be stabilized means (A, C) is not detectable and is reported
 * as DesignStatus::not_detectable. Process noise that leaves a mode on the unit circle
 * unexcited (GQGᵀ does not reach it, e.g. a constant bias with no noise on it) has no
 * stabilizing filter and is refused with DesignStatus::bad_weights. Without a disturbance
 * argument this is the plain filter on x.
 *
 * @note Compare with MATLAB®'s [kest, L, P, Mx] = kalman(sys, Q, R) or [M, P] = dlqe(A, G, C, Q, R)
 *       — the latter is matlab::dlqe in matlab.hpp.
 *
 * @see dare — control-form DARE solver; filter form via dual (Aᵀ, Cᵀ, …)
 * @see design::lqg — pairs this M with LQR’s K
 * @see "Optimal State Estimation" (Simon, 2006), §5.3–5.4
 *
 * @param sys    Discrete-time state-space system (Ts > 0)
 * @param noise  {.process, .measurement}: per-step process covariance (NW × NW, positive
 *               semidefinite) and per-reading measurement variance (NV × NV, positive definite)
 * @return KalmanResult containing the steady-state gains L and M, error covariance P and status
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr auto kalman(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const KalmanNoise<NW, NV, T>&            noise
) {
    return detail::kalman_design<true>(sys, T{0}, noise, NoDisturbance{});
}

/**
 * @brief Steady-state Kalman filter on [x; d]
 *
 * As kalman(sys, noise), with ND integrating disturbances after the plant state. @p disturbance
 * says where d enters: InputDisturbance (B_d = B, C_d = D), OutputDisturbance (a bias on each
 * reading, B_d = 0, C_d = I), PlantDisturbance (the plant's G and H) or DisturbanceModel
 * (explicit B_d and C_d). The process noise enters through G as without one, and
 * @c noise.disturbance is the random walk of d (G on [x; d] is [G 0; 0 I]). The augment must be
 * detectable, rank [A − I, B_d; C, C_d] = NX + ND (so ND ≤ NY), or the design is refused with
 * DesignStatus::not_detectable; a @c noise.disturbance that is not positive definite is refused
 * with DesignStatus::bad_weights. A delayed plant is the model discretize(plant, Ts, InputDelay{θ})
 * returns.
 *
 * @see "Disturbance Models for Offset-Free Model-Predictive Control" (Pannocchia & Rawlings,
 *      2003), AIChE Journal 49(2)
 * @see "Disturbance modeling for offset-free linear model predictive control" (Muske &
 *      Badgwell, 2002), Journal of Process Control 12(5)
 *
 * @param sys          Discrete-time state-space system (Ts > 0)
 * @param noise        {.process, .measurement, .disturbance}: the per-step process covariance
 *                     (positive semidefinite), the random walk of d (positive definite: each d
 *                     is a mode at z = 1 that only its walk moves, so a zero walk is refused
 *                     with DesignStatus::bad_weights), and the per-reading measurement variance
 *                     (positive definite)
 * @param disturbance  Where d enters: InputDisturbance, OutputDisturbance, PlantDisturbance or
 *                     DisturbanceModel
 * @return KalmanResult on [x; d]: steady-state gains L and M, error covariance P and status
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Disturbance, size_t ND = detail::disturbance_traits<Disturbance, NU, NY, NW>::states>
    requires(detail::is_disturbance<Disturbance>::value && !damp::is_same_v<Disturbance, NoDisturbance>)
[[nodiscard]] constexpr auto kalman(
    const StateSpace<NX, NU, NY, T, NW, NV>&                 sys,
    const damp::type_identity_t<KalmanNoise<NW, NV, T, ND>>& noise,
    const Disturbance&                                       disturbance
) {
    return detail::kalman_with_disturbance<true>(sys, T{0}, noise, disturbance);
}

/**
 * @brief Discrete steady-state Kalman filter for a continuous plant with a process-noise density
 *
 * The sampled-data Kalman design: the plant is discretized by ZOH at @p Ts and the process
 * noise density by @ref discretize_process_noise (Van Loan), then the discrete filter is
 * designed as by kalman(sys, noise). The result is the discrete model and its noise: @c sys
 * has @c Ts = Ts and @c G = I (Q_d is a full NX × NX covariance, so the process noise has NX
 * inputs), @c Q = Q_d and @c R = noise.measurement. It deploys as any other KalmanResult.
 *
 * @c noise.process is a density, like the weights of lqr(sys, Ts, w) are a cost integral, so
 * the same number means the same at any rate. @c noise.measurement is the variance of each
 * reading and is used as-is: log the sensor at rest and take the variance, or from a
 * datasheet σ² = ASD² × B with B the noise bandwidth of the configured sensor filter.
 * MATLAB®'s kalmd R/Ts assumes an integrate-and-dump sampler; damp does not.
 *
 * Dispatch, as for lqr(sys, Ts, w):
 *
 * | `sys.Ts` | `Ts` argument | Result |
 * | --- | --- | --- |
 * | 0 | > 0, finite | ZOH + Van Loan noise, then the filter DARE; `sys.Ts = Ts` |
 * | > 0 | any | refused: invalid_sample_time (no second discretization) |
 * | any | ≤ 0 or non-finite | refused: invalid_sample_time |
 *
 * @note Compare with MATLAB®'s kalmd(sys, Qn, Rn, Ts), which passes Rn/Ts as the
 *       measurement variance; matlab::kalmd keeps that convention and forwards here.
 * @see "Optimal State Estimation" (Simon, 2006), §8.1 (continuous and discrete white noise)
 *
 * @param sys    Continuous plant (sys.Ts = 0; uses A, B, C, D, G, H)
 * @param Ts     Sampling time of the filter [s]
 * @param noise  {.process, .measurement}: process-noise density on the NW inputs of G
 *               (positive semidefinite) and per-reading measurement variance on the NV
 *               inputs of H (positive definite, no Ts scaling)
 * @return KalmanResult on the discretized plant (NW = NX), with L, M, P and status
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr auto kalman(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    damp::type_identity_t<T>                 Ts,
    const KalmanNoise<NW, NV, T>&            noise
) {
    return detail::kalman_design<false>(sys, Ts, noise, NoDisturbance{});
}

/**
 * @brief Sampled-data Kalman filter on [x; d] for a continuous plant
 *
 * As kalman(sys, Ts, noise), with a disturbance argument. The process noise on [x; d] is one
 * period of ẋ = Ax + Bu + B_d d + Gw, ḋ = w_d with the densities @c noise.process on w and
 * @c noise.disturbance on w_d (Van Loan, cross terms included: the walk of d moves x within the
 * period), so G on [x; d] is I.
 *
 * @see "Disturbance Models for Offset-Free Model-Predictive Control" (Pannocchia & Rawlings,
 *      2003), AIChE Journal 49(2)
 * @see "Disturbance modeling for offset-free linear model predictive control" (Muske &
 *      Badgwell, 2002), Journal of Process Control 12(5)
 *
 * @param sys          Continuous plant (sys.Ts = 0)
 * @param Ts           Sampling time of the filter [s]
 * @param noise        {.process, .measurement, .disturbance}: process-noise density, per-reading
 *                     measurement variance (no Ts scaling), and the density of d's random walk
 * @param disturbance  Where d enters: InputDisturbance, OutputDisturbance, PlantDisturbance or
 *                     DisturbanceModel
 * @return KalmanResult on the discretized augment [x; d] (G = I), with L, M, P and status
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Disturbance, size_t ND = detail::disturbance_traits<Disturbance, NU, NY, NW>::states>
    requires(detail::is_disturbance<Disturbance>::value && !damp::is_same_v<Disturbance, NoDisturbance>)
[[nodiscard]] constexpr auto kalman(
    const StateSpace<NX, NU, NY, T, NW, NV>&                 sys,
    damp::type_identity_t<T>                                 Ts,
    const damp::type_identity_t<KalmanNoise<NW, NV, T, ND>>& noise,
    const Disturbance&                                       disturbance
) {
    return detail::kalman_with_disturbance<false>(sys, Ts, noise, disturbance);
}

} // namespace design

/**
 * @brief Configuration of a KalmanFilter built from a model rather than from design::kalman
 *
 * An aggregate; name the fields you set. The model is the constructor's first argument, so
 * `KalmanFilter kf{sys, {.Q = Q, .R = R}}` deduces the filter from it. R is the covariance on
 * the inputs of H: the model's own reading (update(y)) has covariance H R Hᵀ.
 *
 * @tparam NW Process-noise inputs of the model (columns of G)
 * @tparam NV Measurement-noise inputs (columns of H)
 */
template<size_t NX, size_t NU, size_t NY, typename T = float, size_t NW = NX, size_t NV = NY>
struct KalmanFilterConfig {
    Matrix<NW, NW, T> Q{Matrix<NW, NW, T>::identity()};  ///< Process-noise covariance on the inputs of G, per step
    Matrix<NV, NV, T> R{Matrix<NV, NV, T>::identity()};  ///< Measurement-noise covariance on the inputs of H, per reading
    ColVec<NX, T>     x0{};                              ///< Initial estimate x̂₀
    Matrix<NX, NX, T> P0{Matrix<NX, NX, T>::identity()}; ///< Initial covariance P₀
};

/**
 * @brief The recursive Kalman filter: the covariance propagates each tick, and the gain is
 *        formed from it in Joseph form
 *
 * The linear member of the recursive family beside the ExtendedKalmanFilter, the
 * UnscentedKalmanFilter and the ErrorStateKalmanFilter, with the same shape: the filter holds
 * its model and configuration, predict(u) takes the input applied now, and each update takes
 * one reading. Built from design::kalman (P∞ seeds P) or from the model and a
 * KalmanFilterConfig. predict(u) propagates P ← A P Aᵀ + G Q Gᵀ; update(y) corrects with the
 * model's own C, D and H R Hᵀ; update(y, C, R) and update(y, C, D, R) correct with a reading of
 * any size NY2 whose C, D or R differ from the model's: a second sensor, a channel whose
 * variance or geometry moves, a reading that covers part of the outputs. The feedthrough takes
 * the input of the last predict, the input the plant holds when the reading is taken.
 *
 * Each update solves S = C P Cᵀ + R by Cholesky (S is symmetric positive definite), forms
 * K = P Cᵀ S⁻¹ and updates P in Joseph form. It returns false, and leaves the estimate and P
 * unchanged, when S is not positive definite.
 *
 * @code
 * KalmanFilter kf{sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0}};
 * kf.predict(u);
 * (void)kf.update(y);               // the model's own C, D and R
 * (void)kf.update(theta, C_enc, r); // a second sensor, y = C x
 * @endcode
 *
 * For a fixed (steady-state) gain, run the design on Luenberger: `Luenberger obs{design::kalman(sys, noise)}`
 * — about 2–13 % of the recursive filter's multiply-accumulates per tick (2.4 % at 10 states, 6 outputs),
 * no per-tick solve, no failure path. Luenberger is also the estimator inside LQG, LQGI, OFLQG and
 * OFMPC, with the disturbance anchor the offset-free laws use.
 *
 * A design on [x; d] (design::kalman with a disturbance) runs here as the absolute filter on the
 * whole augmented state. An input delay (at most one sample) is in the model
 * discretize(plant, Ts, InputDelay{θ}) returns: its register v(k−1) is a state like any other.
 *
 * @tparam NX Filter states
 * @tparam NW Process-noise inputs of the filter's model (columns of G)
 * @tparam NV Measurement-noise inputs (columns of H)
 *
 * @see "Optimal State Estimation" (Simon, 2006), §5.1, §7.6.2 (Joseph form)
 * @see "Optimal Control" (Anderson & Moore, 1990), §8
 */
template<size_t NX, size_t NU, size_t NY, typename T = float, size_t NW = NX, size_t NV = NY>
struct KalmanFilter {
    constexpr KalmanFilter() = default;

    /**
     * @brief From design::kalman (any scalar): P∞ seeds P, and the design's Q and R are the noise.
     *
     * A design with a disturbance model runs as the absolute filter on [x; d]. A failed design
     * does not compile when constant-evaluated; at run time its P is the identity.
     */
    template<typename U, size_t ND>
    constexpr KalmanFilter(const design::KalmanResult<NX, NU, NY, U, NW, NV, ND>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : sys(result.sys.template as<T>()), P_(result.success ? result.P.template as<T>() : Matrix<NX, NX, T>::identity()), Q_(result.Q.template as<T>()), R_(result.R.template as<T>()) {
        if (!result.success && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(result.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    /**
     * @brief From a discrete model, its noise covariances and the initial estimate.
     *
     * A continuous model (Ts not positive and finite) does not compile when constant-evaluated.
     */
    constexpr KalmanFilter(const StateSpace<NX, NU, NY, T, NW, NV>& model, const KalmanFilterConfig<NX, NU, NY, T, NW, NV>& config)
        : sys(model), x(config.x0), P_(config.P0), Q_(config.Q), R_(config.R) {
        if (!damp::finite_positive(model.Ts) && damp::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(design::DesignStatus::invalid_sample_time); // constant-evaluated: A is not a transition matrix
        }
    }

    /** @brief The same filter in another scalar. */
    template<typename U>
    constexpr KalmanFilter(const KalmanFilter<NX, NU, NY, U, NW, NV>& other) // NOLINT(google-explicit-constructor): a scalar conversion
        : sys(other.sys.template as<T>()),
          x(other.x.template as<T>()),
          u_prev_(other.u_prev_.template as<T>()),
          P_(other.P_.template as<T>()),
          Q_(other.Q_.template as<T>()),
          R_(other.R_.template as<T>()) {}

    /** @brief Time update with the input applied now: x̂ ← A x̂ + B u, P ← A P Aᵀ + G Q Gᵀ. The updates' feedthrough takes u until the next predict. */
    constexpr void predict(const ColVec<NU, T>& u) {
        u_prev_ = u;
        x = (sys.A * x) + (sys.B * u);
        P_ = quadratic_form(sys.A, P_) + quadratic_form(sys.G, Q_);
    }

    /** @brief Measurement update with the model's own C, D and H R Hᵀ. False if S is not positive definite. */
    [[nodiscard]] constexpr bool update(const ColVec<NY, T>& y) {
        return correct(ColVec<NY, T>(y - ((sys.C * x) + (sys.D * u_prev_))), sys.C, quadratic_form(sys.H, R_));
    }

    /**
     * @brief Measurement update with a linear reading y = C x + v of any size, v ~ N(0, R)
     *
     * @param y Reading (NY2 rows)
     * @param C Output matrix of the reading on the filter's state
     * @param R Covariance of the reading's noise
     * @return false if S = C P Cᵀ + R is not positive definite: the estimate and P are unchanged
     */
    template<size_t NY2>
    [[nodiscard]] constexpr bool update(const ColVec<NY2, T>& y, const Matrix<NY2, NX, T>& C, const damp::type_identity_t<Matrix<NY2, NY2, T>>& R) {
        return correct(ColVec<NY2, T>(y - (C * x)), C, R);
    }

    /**
     * @brief Measurement update with a linear reading y = C x + D u + v of any size, u the input
     *        of the last predict
     *
     * @param y Reading (NY2 rows)
     * @param C Output matrix of the reading on the filter's state
     * @param D Feedthrough of the held input
     * @param R Covariance of the reading's noise
     * @return false if S = C P Cᵀ + R is not positive definite: the estimate and P are unchanged
     */
    template<size_t NY2>
    [[nodiscard]] constexpr bool update(
        const ColVec<NY2, T>&                             y,
        const Matrix<NY2, NX, T>&                         C,
        const Matrix<NY2, NU, T>&                         D,
        const damp::type_identity_t<Matrix<NY2, NY2, T>>& R
    ) {
        return correct(ColVec<NY2, T>(y - ((C * x) + (D * u_prev_))), C, R);
    }

    /**
     * @brief Fused tick: predict with @p u, correct with @p y from the model's own C, D and R,
     *        return the state.
     *
     * One call is one complete filter tick (StateEstimator). An S that is not positive definite
     * skips the correction and the prediction stands.
     *
     * @param y Measurement taken now
     * @param u Input applied over the previous interval
     * @return The updated state
     */
    constexpr const ColVec<NX, T>& step(const ColVec<NY, T>& y, const ColVec<NU, T>& u) {
        predict(u);
        (void)update(y);
        return x;
    }

    /** @brief The state estimate x̂. */
    [[nodiscard]] constexpr const ColVec<NX, T>& state() const { return x; }

    /** @brief The running error covariance P. */
    [[nodiscard]] constexpr const Matrix<NX, NX, T>& covariance() const { return P_; }

    /** @brief The filter's model. */
    [[nodiscard]] constexpr const StateSpace<NX, NU, NY, T, NW, NV>& model() const { return sys; }

    /** @brief Process-noise covariance Q. */
    [[nodiscard]] constexpr const Matrix<NW, NW, T>& process_noise_covariance() const { return Q_; }

    /** @brief Measurement-noise covariance R of the model's own reading, on the inputs of H. */
    [[nodiscard]] constexpr const Matrix<NV, NV, T>& measurement_noise_covariance() const { return R_; }

    /** @brief Overwrite the estimate between updates (a clamp to a physical bound, an angle wrap); P and the held input are kept. */
    constexpr void set_state(const ColVec<NX, T>& x_new) { x = x_new; }

    /** @brief Forget the applied input; the estimate and its covariance are left in place. */
    constexpr void reset() { u_prev_ = ColVec<NU, T>{}; }

    /** @brief Replace the estimate and its covariance, and forget the applied input. Q and R stay. */
    constexpr void reset(const ColVec<NX, T>& x0, const Matrix<NX, NX, T>& P0) {
        x = x0;
        P_ = P0;
        u_prev_ = ColVec<NU, T>{};
    }

private:
    template<size_t, size_t, size_t, typename, size_t, size_t>
    friend struct KalmanFilter;

    /**
     * @brief Joseph update with innovation ν, output matrix C and reading covariance R
     *
     * S = C P Cᵀ + R is symmetric positive definite, so K = P Cᵀ S⁻¹ comes from S Kᵀ = C P by
     * Cholesky. P ← (I − K C) P (I − K C)ᵀ + K R Kᵀ keeps P symmetric and positive semidefinite
     * for any K.
     */
    template<size_t NYS>
    [[nodiscard]] constexpr bool correct(const ColVec<NYS, T>& nu, const Matrix<NYS, NX, T>& C, const Matrix<NYS, NYS, T>& R) {
        const Matrix<NYS, NYS, T> S = quadratic_form(C, P_) + R;
        const auto                K_t = mat::cholesky_solve(S, Matrix<NYS, NX, T>(C * P_));
        if (!K_t) {
            return false; // S not positive definite: fail closed
        }
        const Matrix<NX, NYS, T> K = K_t.value().transpose();
        x += K * nu;
        P_ = quadratic_form(Matrix<NX, NX, T>(Matrix<NX, NX, T>::identity() - (K * C)), P_) + quadratic_form(K, R);
        return true;
    }

    StateSpace<NX, NU, NY, T, NW, NV> sys{};     ///< The filter's model
    ColVec<NX, T>                     x{};       ///< Estimate
    ColVec<NU, T>                     u_prev_{}; ///< The input the plant holds, for the feedthrough
    Matrix<NX, NX, T>                 P_{};      ///< Error covariance
    Matrix<NW, NW, T>                 Q_{};      ///< Process noise covariance
    Matrix<NV, NV, T>                 R_{};      ///< Measurement noise covariance of the model's own reading
};

/** @brief Deduce the filter from its design, keeping the result's scalar: KalmanFilter kf{res.as<float>()}. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t ND>
KalmanFilter(const design::KalmanResult<NX, NU, NY, T, NW, NV, ND>&) -> KalmanFilter<NX, NU, NY, T, NW, NV>;

/** @brief Deduce the filter from its model: KalmanFilter kf{sys, {.Q = Q, .R = R}}. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
KalmanFilter(const StateSpace<NX, NU, NY, T, NW, NV>&, const KalmanFilterConfig<NX, NU, NY, T, NW, NV>&) -> KalmanFilter<NX, NU, NY, T, NW, NV>;

} // namespace damp
