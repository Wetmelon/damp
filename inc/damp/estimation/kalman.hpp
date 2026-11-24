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
 * (steady-state) gain, run the design on FixedGainObserver (fixed_gain_observer.hpp), the fixed-gain observer the
 * LQ and MPC laws hold: `FixedGainObserver obs{design::kalman(sys, noise)}`. KalmanFilter propagates its
 * covariance each tick, for a reading whose R or C changes at run time.
 *
 * @code
 * using namespace damp;
 * // discrete plant: per-step process covariance, per-reading measurement variance
 * constexpr auto res = design::kalman(sys, design::KalmanNoise{.process = Qd, .measurement = sigma2});
 * // continuous plant: process-noise density (Van Loan'd), the same per-reading variance
 * constexpr auto res_c = design::kalman(sys_c, Ts, design::KalmanNoise{.process = Qc, .measurement = sigma2});
 * // a filter on [x; d]: noise.disturbance is the random walk of d
 * // (the disturbance spec sets the count, so the noise is a plain braced list)
 * constexpr auto res_d = design::kalman(sys, design::KalmanNoise{.process = Qd, .measurement = sigma2, .disturbance = Qdd}, design::InputDisturbance{});
 * static constinit FixedGainObserver obs{res};     // the fixed gain M (fixed_gain_observer.hpp); a failed design does not compile
 * static constinit KalmanFilter filt{res};  // the recursive filter, P∞ its initial P
 * // each tick: auto xhat = obs.step(y, u);  // or obs.predict(u); obs.update(y, u); for a split tick
 * @endcode
 */

#include "damp/backend.hpp"
#include "damp/controllers/detail/law_detail.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/options.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/sqrt_covariance.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/const_matrix.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/pattern_matrix.hpp"
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
 * A design takes the noise among its trailing options, and the noise deduces its sizes from its
 * fields: a matrix gives its own, a scalar is 1 × 1 (`KalmanNoise{.process = diag(1e-6, 1e-2),
 * .measurement = 0.01}` is KalmanNoise<2, 1>). The design refuses one not sized to its plant.
 *
 * @code
 * design::kalman(continuous_sys, Ts, design::KalmanNoise{.process = Qc, .measurement = sigma2});
 * design::kalman(discrete_sys, design::KalmanNoise{.process = Qd, .measurement = sigma2});
 * @endcode
 *
 * With a disturbance model (a filter on [x; d]), @c disturbance is the random walk of d: how
 * fast the disturbance may move, in the plant's time domain like @c process (a density on a
 * continuous plant, a covariance per step on a discrete one). Without one it has no entries.
 *
 * @tparam NW Process-noise inputs (columns of G); NV Measurement-noise inputs (columns of H). A plant
 *            that declares none (`StateSpace<NX, NU, NY>`, NW = NV = 0) has noise on every state and
 *            reading: the design reads it with G = I and H = I, so the noise is NX × NX and NY × NY.
 *            A model discretized with an InputDelay or an OutputDelay counts its delay states among
 *            its states, so they take noise too: the continuous plant and the delay as an option
 *            keep them exact.
 * @tparam ND Disturbance states (0: no disturbance model).
 */
template<size_t NW, size_t NV, typename T = double, size_t ND = 0>
struct KalmanNoise {
    using option_kind = kind::KalmanNoise; ///< Its kind of design option
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

// Sizes from the fields, in order: `KalmanNoise{.process = diag(…), .measurement = 0.01}`. A matrix
// gives its size, a scalar (one floating type throughout) is 1 × 1, and .disturbance is optional. The
// scalar guides take a floating type only, so two matrices of one type deduce as matrices.

/** @brief Sizes from the fields: process and measurement matrices. */
template<size_t NW, size_t NV, typename T>
KalmanNoise(Matrix<NW, NW, T>, Matrix<NV, NV, T>) -> KalmanNoise<NW, NV, T>;
/** @brief Sizes from the fields: a process matrix, one reading. */
template<size_t NW, typename T>
    requires damp::is_floating_point_v<T>
KalmanNoise(Matrix<NW, NW, T>, T) -> KalmanNoise<NW, 1, T>;
/** @brief Sizes from the fields: one process input, a measurement matrix. */
template<size_t NV, typename T>
    requires damp::is_floating_point_v<T>
KalmanNoise(T, Matrix<NV, NV, T>) -> KalmanNoise<1, NV, T>;
/** @brief Sizes from the fields: one process input, one reading. */
template<typename T>
    requires damp::is_floating_point_v<T>
KalmanNoise(T, T) -> KalmanNoise<1, 1, T>;
/** @brief Sizes from the fields, with the disturbance's walk: three matrices. */
template<size_t NW, size_t NV, size_t ND, typename T>
KalmanNoise(Matrix<NW, NW, T>, Matrix<NV, NV, T>, Matrix<ND, ND, T>) -> KalmanNoise<NW, NV, T, ND>;
/** @brief Sizes from the fields, with the disturbance's walk: one disturbance. */
template<size_t NW, size_t NV, typename T>
    requires damp::is_floating_point_v<T>
KalmanNoise(Matrix<NW, NW, T>, Matrix<NV, NV, T>, T) -> KalmanNoise<NW, NV, T, 1>;
/** @brief Sizes from the fields, with the disturbance's walk: one reading. */
template<size_t NW, size_t ND, typename T>
    requires damp::is_floating_point_v<T>
KalmanNoise(Matrix<NW, NW, T>, T, Matrix<ND, ND, T>) -> KalmanNoise<NW, 1, T, ND>;
/** @brief Sizes from the fields, with the disturbance's walk: one reading, one disturbance. */
template<size_t NW, typename T>
    requires damp::is_floating_point_v<T>
KalmanNoise(Matrix<NW, NW, T>, T, T) -> KalmanNoise<NW, 1, T, 1>;
/** @brief Sizes from the fields, with the disturbance's walk: one process input. */
template<size_t NV, size_t ND, typename T>
    requires damp::is_floating_point_v<T>
KalmanNoise(T, Matrix<NV, NV, T>, Matrix<ND, ND, T>) -> KalmanNoise<1, NV, T, ND>;
/** @brief Sizes from the fields, with the disturbance's walk: one process input, one disturbance. */
template<size_t NV, typename T>
    requires damp::is_floating_point_v<T>
KalmanNoise(T, Matrix<NV, NV, T>, T) -> KalmanNoise<1, NV, T, 1>;
/** @brief Sizes from the fields, with the disturbance's walk: one process input, one reading. */
template<size_t ND, typename T>
    requires damp::is_floating_point_v<T>
KalmanNoise(T, T, Matrix<ND, ND, T>) -> KalmanNoise<1, 1, T, ND>;
/** @brief Sizes from the fields, with the disturbance's walk: scalars throughout. */
template<typename T>
    requires damp::is_floating_point_v<T>
KalmanNoise(T, T, T) -> KalmanNoise<1, 1, T, 1>;

namespace detail {

/** @brief The noise inputs a Kalman design reads: @p N (one per state, or per reading) when the plant declares none (@p NW = 0). */
template<size_t N, size_t NW>
inline constexpr size_t noise_inputs = (NW == 0) ? N : NW;

/** @brief KalmanNoise sized to the plant: NX × NX process and NY × NY measurement noise when it declares no noise inputs. */
template<size_t NX, size_t NY, typename T, size_t NW, size_t NV, size_t ND = 0>
using plant_noise = KalmanNoise<noise_inputs<NX, NW>, noise_inputs<NY, NV>, T, ND>;

/**
 * @brief The plant as a Kalman design reads it: with no process-noise inputs (NW = 0), noise on
 *        every state (G = I); with no measurement-noise inputs (NV = 0), on every reading (H = I).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr decltype(auto) noise_model(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    if constexpr (NW > 0 && NV > 0) {
        return (sys); // the plant itself
    } else {
        StateSpace<NX, NU, NY, T, noise_inputs<NX, NW>, noise_inputs<NY, NV>> out{.A = sys.A, .B = sys.B, .C = sys.C, .D = sys.D, .Ts = sys.Ts};
        if constexpr (NW > 0) {
            out.G = sys.G;
        }
        if constexpr (NV > 0) {
            out.H = sys.H;
        }
        return out;
    }
}

} // namespace detail

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
 * FixedGainObserver runs M in the current form (`FixedGainObserver obs{res}`); KalmanFilter starts from P.
 * Use .as<float>() to convert for embedded deployment.
 *
 * With a disturbance (design::kalman(sys[, Ts], noise, disturbance)), @c sys is the filter's model:
 * the plant, then the ND disturbance states; @c M, @c L and @c P are the steady-state gains and
 * covariance on that whole model.
 *
 * A state-feedback design with ObserverPoles in place of its KalmanNoise carries its filter in
 * this shape too, on the same model: @c M and @c L placed, @c Q, @c R and @c P zero (no noise).
 *
 * @tparam ND Disturbance states, last in the filter's state (0: none)
 *
 * @see "Optimal State Estimation" (Simon, 2006), §5.3
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0, size_t ND = 0>
struct [[nodiscard]] KalmanResult {
    /** @brief The design in another scalar: as<U>()'s type. */
    template<typename U>
    using as_t = KalmanResult<NX, NU, NY, U, NW, NV, ND>;

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
        return set_status(result, DesignStatus::invalid_sample_time);
    }
    // The DARE never reads B or D, so a non-finite one would ship in result.sys with success.
    if (!finite_filter_model(sys, Q, R)) {
        return set_status(result, DesignStatus::bad_weights);
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
        return set_status(result, dual_status(sol.status));
    }
    // Scrub DARE round-off asymmetry before gain solve (S = CPCᵀ + R must stay SPD).
    result.P = static_cast<T>(0.5) * (sol.X + sol.X.transpose());

    // Innovation gain: M = PCᵀS⁻¹ → solve S Mᵀ = CP (Cholesky, LU if S semi-def)
    const Matrix<NY, NY, T> S = sys.C * result.P * sys.C.t() + R_eff;
    const auto              M_opt = mat::solve(S, Matrix<NY, NX, T>(sys.C * result.P));
    if (!M_opt) {
        return set_status(result, DesignStatus::singular); // singular innovation covariance S
    }
    result.M = M_opt.value().transpose();
    result.L = sys.A * result.M; // the predictor gain: x̂[k+1|k] = A x̂[k|k]
    return set_status(result, DesignStatus::ok);
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

/**
 * @brief The process noise over one period on the model of a filter whose readings arrive late:
 *        on the continuous states c (the plant's, and d's with a disturbance) at the period's end,
 *        and on each late reading as its first delay state loads it (Van Loan)
 *
 * A late reading's first delay state loads R_i·c(s_i), s_i = Ts − τ_i into the period for a
 * fraction τ_i of a sample (0 for whole samples: c(k) itself, no new noise). With
 * n(t) = ∫₀ᵗ e^{A(t−η)} G w(η) dη and Q(t) its covariance, cov(n(Ts)) = Q(Ts),
 * cov(n(Ts), R_i n(s_i)) = e^{A(Ts−s_i)} Q(s_i) R_iᵀ and cov(R_i n(s_i), R_j n(s_j)) =
 * R_i e^{A(s_i−m)} Q(m) e^{Aᵀ(s_j−m)} R_jᵀ, m = min(s_i, s_j). The input's delay states and the
 * older readings take none. @p row maps each state of c to the model's, @p first the first
 * reading's delay state.
 * @return nullopt on a non-finite entry or result
 */
template<size_t NQ, size_t NC, size_t NG, size_t NY, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NQ, NQ, T>> late_reading_noise(
    const Matrix<NC, NC, T>&                              A,
    const Matrix<NC, NG, T>&                              G,
    const Matrix<NG, NG, T>&                              Qc,
    const Matrix<NY, NC, T>&                              R,
    const damp::array<::damp::detail::DelaySplit<T>, NY>& outputs,
    T                                                     Ts,
    const damp::array<size_t, NC>&                        row,
    size_t                                                first
) {
    const auto Q_T = discretize_process_noise(A, G, Qc, Ts);
    if (!Q_T) {
        return damp::nullopt;
    }
    Matrix<NQ, NQ, T> Q{};
    for (size_t a = 0; a < NC; ++a) {
        for (size_t b = 0; b < NC; ++b) {
            Q(row[a], row[b]) = (*Q_T)(a, b);
        }
    }
    damp::array<T, NY> s{};
    for (size_t i = 0; i < NY; ++i) {
        s[i] = (outputs[i].fraction > T{0}) ? Ts - outputs[i].fraction : T{0};
    }
    const auto held = [&](T t) { return discretize_process_noise(A, G, Qc, t); };
    for (size_t i = 0; i < NY; ++i) {
        if (!(s[i] > T{0})) {
            continue; // whole samples: the reading loads c(k), which holds no new noise
        }
        const auto Qs = held(s[i]);
        const auto E = mat::expm(Matrix<NC, NC, T>(A * (Ts - s[i])));
        if (!Qs || !E) {
            return damp::nullopt;
        }
        const Matrix<NC, NC, T> EQ = *E * *Qs;
        for (size_t a = 0; a < NC; ++a) {
            T v{0};
            for (size_t c = 0; c < NC; ++c) {
                v += EQ(a, c) * R(i, c);
            }
            Q(row[a], first + i) = v;
            Q(first + i, row[a]) = v;
        }
        for (size_t j = 0; j < NY; ++j) {
            if (!(s[j] > T{0})) {
                continue;
            }
            const T    m = damp::min(s[i], s[j]);
            const auto Qm = held(m);
            const auto Ei = mat::expm(Matrix<NC, NC, T>(A * (s[i] - m)));
            const auto Ej = mat::expm(Matrix<NC, NC, T>(A * (s[j] - m)));
            if (!Qm || !Ei || !Ej) {
                return damp::nullopt;
            }
            const Matrix<NC, NC, T> M = *Ei * *Qm * Ej->t();
            T                       v{0};
            for (size_t a = 0; a < NC; ++a) {
                for (size_t b = 0; b < NC; ++b) {
                    v += R(i, a) * M(a, b) * R(j, b);
                }
            }
            Q(first + i, first + j) = v;
        }
    }
    return Q;
}

/**
 * @brief The sampled-data Kalman design on a continuous plant (see kalman(sys, Ts, noise)), on
 *        discretize(sys, Ts) or, with @p delay, on discretize(sys, Ts, delay): the filter carries
 *        the delay states.
 *
 * - No delay: the plant's ZOH model and Q_d over the period (discretize_process_noise), G = I.
 * - An InputDelay: the process noise acts on the plant, not through the delay: Q_d is the plant's,
 *   the delay states take none, and G = [I; 0]. The delay states hold the inputs applied, which
 *   predict shifts in, so their covariance is zero and the readings do not move them.
 * - A late reading (PlantDelays with an OutputDelay), and its input too: the process noise acts
 *   over the period on the plant's state, and on each late reading up to when it was taken
 *   (late_reading_noise); the input's delay states and the older readings take none. The
 *   measurement noise is the reading's as it arrives.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Delay = NoDelay>
[[nodiscard]] constexpr auto kalman_sampled(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T Ts, const KalmanNoise<NW, NV, T>& noise, const Delay& delay = Delay{}) {
    constexpr size_t NXE = NX + delay_states<Delay, NU, NY>;
    if constexpr (damp::is_same_v<Delay, NoDelay>) {
        using Result = KalmanResult<NX, NU, NY, T, NX, NV>;
        if (!damp::isfinite(sys.Ts) || sys.Ts != T{0} || !damp::finite_positive(Ts)) {
            return refused<Result>(DesignStatus::invalid_sample_time); // one discretization, at a valid Ts
        }
        const auto plant_d = ::damp::discretize(StateSpace{.A = sys.A, .B = sys.B, .C = sys.C, .D = sys.D}, Ts);
        if (!plant_d) {
            return refused<Result>(plant_d.status); // a non-finite A or B, or an overflow
        }
        const auto Qd = discretize_process_noise(sys.A, sys.G, noise.process, Ts);
        if (!Qd) {
            return refused<Result>(DesignStatus::bad_weights); // a non-finite G or noise entry
        }
        const StateSpace sys_d{.A = plant_d->A, .B = plant_d->B, .C = sys.C, .D = sys.D, .G = Matrix<NX, NX, T>::identity(), .H = sys.H, .Ts = Ts};
        return kalman_plain(sys_d, KalmanNoise<NX, NV, T>{.process = *Qd, .measurement = noise.measurement});
    } else if constexpr (OptionOf<Delay, kind::InputDelay>) {
        using Result = KalmanResult<NXE, NU, NY, T, NX, NV>;
        const auto m = sampled_model(sys, Ts, Domain::sampled, delay);
        if (m.status != DesignStatus::ok) {
            return refused<Result>(m.status);
        }
        const auto Qd = discretize_process_noise(sys.A, sys.G, noise.process, Ts);
        if (!Qd) {
            return refused<Result>(DesignStatus::bad_weights); // a non-finite noise entry
        }
        StateSpace<NXE, NU, NY, T, NX, NV> model{.A = m.law.A, .B = m.law.B, .C = m.law.C, .D = m.law.D, .H = sys.H, .Ts = Ts};
        model.G = Matrix<NXE, NX, T>{};
        model.G.template block<NX, NX>(0, 0) = Matrix<NX, NX, T>::identity();
        return kalman_plain(model, KalmanNoise<NX, NV, T>{.process = *Qd, .measurement = noise.measurement});
    } else {
        using Result = KalmanResult<NXE, NU, NY, T, NXE, NV>;
        const auto m = sampled_model(sys, Ts, Domain::sampled, delay);
        const auto outputs = ::damp::detail::channel_delay_splits<NY>(delay.output.seconds, Ts, output_delay_states<Delay>);
        if (m.status != DesignStatus::ok || !outputs) {
            return refused<Result>((m.status != DesignStatus::ok) ? m.status : DesignStatus::invalid_parameter);
        }
        damp::array<size_t, NX> row{};
        for (size_t i = 0; i < NX; ++i) {
            row[i] = i;
        }
        const auto Q = late_reading_noise<NXE>(sys.A, sys.G, noise.process, sys.C, *outputs, Ts, row, NX + (input_delay_states<Delay> * NU));
        if (!Q) {
            return refused<Result>(DesignStatus::bad_weights); // a non-finite noise entry
        }
        StateSpace<NXE, NU, NY, T, NXE, NV> model{.A = m.law.A, .B = m.law.B, .C = m.law.C, .D = m.law.D, .G = Matrix<NXE, NXE, T>::identity(), .H = sys.H, .Ts = Ts};
        return kalman_plain(model, KalmanNoise<NXE, NV, T>{.process = *Q, .measurement = noise.measurement});
    }
}

/** @brief A plain design result as the result of a design with options: its ND. */
template<size_t ND, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr KalmanResult<NX, NU, NY, T, NW, NV, ND> with_disturbance_states(const KalmanResult<NX, NU, NY, T, NW, NV>& result) {
    return KalmanResult<NX, NU, NY, T, NW, NV, ND>{result.sys, result.Q, result.R, result.L, result.M, result.P, result.success, result.status};
}

/** @brief A filter design with the channel @p Bm of NM > 0 measured inputs after u (with_measured_inputs): B and D widened, the gain kept. */
template<size_t NXE, size_t NM, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t ND>
[[nodiscard]] constexpr KalmanResult<NX, NU + NM, NY, T, NW, NV, ND> widened(const KalmanResult<NX, NU, NY, T, NW, NV, ND>& result, const damp::optional<Matrix<NXE, NM, T>>& Bm) {
    KalmanResult<NX, NU + NM, NY, T, NW, NV, ND> out{.Q = result.Q, .R = result.R, .L = result.L, .M = result.M, .P = result.P, .success = result.success, .status = result.status};
    out.sys = StateSpace<NX, NU + NM, NY, T, NW, NV>{.A = result.sys.A, .C = result.sys.C, .G = result.sys.G, .H = result.sys.H, .Ts = result.sys.Ts};
    out.sys.B.template block<NX, NU>(0, 0) = result.sys.B;
    out.sys.D.template block<NY, NU>(0, 0) = result.sys.D;
    if (Bm) {
        out.sys.B.template block<NXE, NM>(0, NU) = *Bm;
    } else if (out.success) {
        out = refused<decltype(out)>(DesignStatus::bad_weights); // a non-finite B_d
    }
    return out;
}

/**
 * @brief A filter design that takes a measured disturbance as NM known inputs after u: its model's
 *        B is [B, B_m] and D is [D, 0], B_m the disturbance's channel @p Bm on the model before
 *        the estimated disturbance (measured_channel; zero on d's rows, which do not move).
 *
 * The filter DARE reads neither B nor D, so M, L and P are the design's without it: d is not
 * noise. A channel that is not finite (nullopt) refuses the design with DesignStatus::bad_weights.
 * With none (NM = 0) the design is returned as it is.
 */
template<size_t NXE, size_t NM, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t ND>
[[nodiscard]] constexpr KalmanResult<NX, NU + NM, NY, T, NW, NV, ND> with_measured_inputs(
    const KalmanResult<NX, NU, NY, T, NW, NV, ND>& result,
    const damp::optional<Matrix<NXE, NM, T>>&      Bm
) {
    static_assert(NXE + ND == NX, "The measured channel is on the filter's model before its disturbance states.");
    if constexpr (NM == 0) {
        static_cast<void>(Bm);
        return result;
    } else {
        return widened(result, Bm);
    }
}

/**
 * @brief The Kalman design, plain or on [x; d]: the filter's model (the plant, its input delay states
 *        when @p delay is an InputDelay, then the disturbance), in the plant's own time domain
 *        (@p D plant) or sampled at @p Ts.
 */
template<Domain D, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t ND, typename Spec, typename Delay = NoDelay>
[[nodiscard]] constexpr auto kalman_design(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        Ts,
    const KalmanNoise<NW, NV, T, ND>&        noise,
    const Spec&                              spec,
    const Delay&                             delay = Delay{}
) {
    constexpr size_t NXE = NX + delay_states<Delay, NU, NY>;
    constexpr size_t NXK = NXE + ND;
    if constexpr (ND == 0) {
        static_cast<void>(spec);
        if constexpr (D == Domain::plant) {
            static_cast<void>(Ts);
            return kalman_plain(sys, noise);
        } else {
            return kalman_sampled(sys, Ts, noise, delay);
        }
    } else {
        static_assert(ND <= NY, "Estimating ND disturbances needs at least ND independent readings (ND <= NY).");
        constexpr size_t NWK = (D == Domain::plant) ? (NW + ND) : NXK;
        using Result = KalmanResult<NXK, NU, NY, T, NWK, NV, ND>;
        const auto m = sampled_model(sys, Ts, D, delay);
        if (m.status != DesignStatus::ok) {
            return refused<Result>(m.status);
        }
        if (!finite_filter_model(m.law, noise.process, noise.measurement)) {
            return refused<Result>(DesignStatus::bad_weights); // before the detectability test reads C
        }
        const auto given = disturbance_channels(spec, sys, m.law, Ts);
        const auto ch = on_late_readings(spec, given, sys, Ts, delay); // d on a late reading as it was when taken
        if (!augment_detectable(m.law.A, m.law.C, ch.Bd, ch.Cd)) {
            return refused<Result>(DesignStatus::not_detectable);
        }
        auto                aug = disturbance_augment<NWK>(m.law, ch.Bd, ch.Cd, (D == Domain::plant) ? sys.Ts : Ts);
        Matrix<NWK, NWK, T> Q{};
        if constexpr (D == Domain::plant) {
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
                return refused<Result>(DesignStatus::bad_weights); // a non-finite noise or channel entry
            }
            aug.G = Matrix<NXK, NXK, T>::identity();
            for (size_t i = 0; i < NA; ++i) {
                for (size_t j = 0; j < NA; ++j) {
                    Q((i < NX) ? i : (i - NX + NXE), (j < NX) ? j : (j - NX + NXE)) = (*Qd)(i, j);
                }
            }
            if constexpr (output_delay_states<Delay> > 0) {
                //! A late reading also holds the noise up to when it was taken: [x; d] from it, the
                //! readings on x, and on d through C_d where d is held (an InputDisturbance's rides u).
                const auto outputs = ::damp::detail::channel_delay_splits<NY>(delay.output.seconds, Ts, output_delay_states<Delay>);
                if (!outputs) {
                    return refused<Result>(DesignStatus::invalid_parameter);
                }
                Matrix<NY, NA, T> R{};
                R.template block<NY, NX>(0, 0) = sys.C;
                if constexpr (held_on_late_readings<Spec, Delay>) {
                    R.template block<NY, ND>(0, NX) = given.Cd;
                }
                damp::array<size_t, NA> row{};
                for (size_t i = 0; i < NA; ++i) {
                    row[i] = (i < NX) ? i : (i - NX + NXE);
                }
                const auto Ql = late_reading_noise<NXK>(Ac, Gc, Qc, R, *outputs, Ts, row, NX + (input_delay_states<Delay> * NU));
                if (!Ql) {
                    return refused<Result>(DesignStatus::bad_weights);
                }
                Q = *Ql;
            }
        }
        return with_disturbance_states<ND>(kalman_plain(aug, KalmanNoise<NWK, NV, T>{.process = Q, .measurement = noise.measurement}));
    }
}

/** @brief The options kalman takes: its noise (required), a disturbance spec, a MeasuredDisturbance, an InputDelay, an OutputDelay. */
using kalman_options = option_set<kind::InputDelay, kind::OutputDelay, kind::Disturbance, kind::MeasuredDisturbance, kind::KalmanNoise>;

/** @brief A filter design @p kf on @p sys with the measured disturbance @p measured as known inputs after u (with_measured_inputs); itself without one (NM = 0). */
template<Domain D, size_t NXK, size_t NU, size_t NY, typename T, size_t NWK, size_t NVK, size_t ND, size_t NX, size_t NYP, size_t NW, size_t NV, size_t NM, typename Delay>
[[nodiscard]] constexpr auto measured_filter(
    const KalmanResult<NXK, NU, NY, T, NWK, NVK, ND>& kf,
    const StateSpace<NX, NU, NYP, T, NW, NV>&         sys,
    T                                                 Ts,
    const MeasuredDisturbance<NX, NM, T>&             measured,
    const Delay&                                      delay
) {
    return with_measured_inputs(kf, measured_channel<D>(sys, Ts, measured, delay));
}

/** @brief The disturbance states of the disturbance spec among @p Opts (0 without one). */
template<size_t NX, size_t NU, size_t NY, size_t NW, typename... Opts>
inline constexpr size_t disturbance_states_in = disturbance_states<option_t<kind::Disturbance, NoDisturbance, Opts...>, NU, NY, noise_inputs<NX, NW>>;

/** @brief The KalmanNoise a filter on this plant with the options @p Opts takes. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename... Opts>
using noise_for = plant_noise<NX, NY, T, NW, NV, disturbance_states_in<NX, NU, NY, NW, Opts...>>;

/** @brief The KalmanNoise among @p Opts, checked against the size the plant needs. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename... Opts>
[[nodiscard]] constexpr noise_for<NX, NU, NY, T, NW, NV, Opts...> noise_of(const Opts&... opts) {
    using Noise = option_t<kind::KalmanNoise, noise_for<NX, NU, NY, T, NW, NV, Opts...>, Opts...>;
    constexpr bool sized = damp::is_same_v<Noise, noise_for<NX, NU, NY, T, NW, NV, Opts...>>;
    static_assert(
        sized,
        "The KalmanNoise is sized to the plant: .process NX x NX (or the columns of G when the plant names them), .measurement NY x NY, "
        "and with a disturbance spec .disturbance ND x ND (its states), in the plant's scalar."
    );
    if constexpr (sized) {
        return option<kind::KalmanNoise>(noise_for<NX, NU, NY, T, NW, NV, Opts...>{}, opts...);
    } else {
        (static_cast<void>(opts), ...);
        return noise_for<NX, NU, NY, T, NW, NV, Opts...>{}; // the message says why: nothing further
    }
}

/** @brief design::kalman with its trailing options, in the plant's own domain (@p D plant) or at @p Ts. */
template<Domain D, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename... Opts>
[[nodiscard]] constexpr auto kalman_with(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T Ts, const Opts&... opts) {
    kalman_options::check<Opts...>();
    static_assert(option_count<kind::KalmanNoise, Opts...> == 1, "kalman needs its noise among the options: kalman(plant, Ts, KalmanNoise{.process = ..., .measurement = ...}).");
    constexpr bool kModelFits = check_delays<D, Opts...>() && check_disturbance<NX, NY, T, Opts...>();
    const auto     noise = noise_of<NX, NU, NY, T, NW, NV>(opts...);
    const auto     delay = delays_of(opts...);
    const auto     measured = measured_of<NX, T>(opts...);
    using Disturbance = option_t<kind::Disturbance, NoDisturbance, Opts...>;
    if constexpr (!kModelFits) {
        static_cast<void>(noise);
        static_cast<void>(delay);
        static_cast<void>(measured);
        return KalmanResult<NX, NU, NY, T>{}; // the message says why: nothing further
    } else if constexpr (DisturbanceSpec<Disturbance>) {
        const auto kf = kalman_design<D>(noise_model(sys), Ts, noise, the<kind::Disturbance>(opts...), delay);
        return measured_filter<D>(kf, sys, Ts, measured, delay);
    } else {
        return measured_filter<D>(kalman_design<D>(noise_model(sys), Ts, noise, NoDisturbance{}, delay), sys, Ts, measured, delay);
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
 * The result runs on FixedGainObserver, the fixed-gain observer (`FixedGainObserver obs{res}`: M held, no
 * per-tick covariance), or on KalmanFilter, which starts from P∞ and propagates P each tick.
 *
 * Requires a discrete plant (Ts > 0): there is no continuous-time (Kalman–Bucy) design, so a
 * continuous plant is refused with DesignStatus::invalid_sample_time. For a continuous plant
 * with a noise density, pass the sample time instead: kalman(sys, Ts, noise) discretizes
 * both. A dual DARE that cannot be stabilized means (A, C) is not detectable and is reported
 * as DesignStatus::not_detectable. Process noise that leaves a mode on the unit circle
 * unexcited (GQGᵀ does not reach it, e.g. a constant bias with no noise on it) has no
 * stabilizing filter and is refused with DesignStatus::bad_weights.
 *
 * @p opts trail in any order: a disturbance spec puts the filter on [x; d], ND integrating
 * disturbances after the plant state; without one it is the plain filter on x. A
 * MeasuredDisturbance{.Bd} adds known inputs after u, B on the model [B, B_d] and D [D, 0], the
 * gain unchanged (see kalman(sys, Ts, opts…)). The spec says
 * where d enters: InputDisturbance (B_d = B, C_d = D), OutputDisturbance (a bias on each reading,
 * B_d = 0, C_d = I), PlantDisturbance (the plant's G and H) or DisturbanceModel (explicit B_d and
 * C_d). The process noise enters through G as without one, and @c noise.disturbance is the random
 * walk of d (G on [x; d] is [G 0; 0 I]). The augment must be detectable, rank [A − I, B_d; C, C_d]
 * = NX + ND (so ND ≤ NY), or the design is refused with DesignStatus::not_detectable; a
 * @c noise.disturbance that is not positive definite is refused with DesignStatus::bad_weights.
 * A delayed plant, its input late or its readings, is the model discretize(plant, Ts,
 * InputDelay{θ}, OutputDelay{θ}) returns, or kalman(plant, Ts, noise, InputDelay<N>{θ},
 * OutputDelay<M>{θ}) from the continuous plant; a delay among this form's options is refused at
 * compile time.
 *
 * @note Compare with MATLAB®'s [kest, L, P, Mx] = kalman(sys, Q, R) or [M, P] = dlqe(A, G, C, Q, R)
 *       — the latter is matlab::dlqe in matlab.hpp.
 *
 * @see dare — control-form DARE solver; filter form via dual (Aᵀ, Cᵀ, …)
 * @see design::state_feedback with a KalmanNoise — pairs this M with LQR’s K
 * @see "Optimal State Estimation" (Simon, 2006), §5.3–5.4
 *
 * @param sys    Discrete-time state-space system (Ts > 0)
 * @see "Disturbance Models for Offset-Free Model-Predictive Control" (Pannocchia & Rawlings,
 *      2003), AIChE Journal 49(2)
 * @see "Disturbance modeling for offset-free linear model predictive control" (Muske &
 *      Badgwell, 2002), Journal of Process Control 12(5)
 * @param opts   The KalmanNoise{.process, .measurement[, .disturbance]}, required: per-step process
 *               covariance (NW × NW, positive semidefinite), per-reading measurement variance
 *               (NV × NV, positive definite) and, with a disturbance spec, the random walk of d
 *               (positive definite); and a disturbance spec, or none
 * @return KalmanResult, on [x; d] with a disturbance spec: the gains L and M, error covariance P and status
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename... Opts>
[[nodiscard]] constexpr auto kalman(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const Opts&... opts) {
    return detail::kalman_with<detail::Domain::plant>(sys, T{0}, opts...);
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
 * @p opts trail in any order, at most one of each kind:
 * - A disturbance spec, as for kalman(sys, noise, opts…). The process noise on [x; d] is one
 *   period of ẋ = Ax + Bu + B_d d + Gw, ḋ = w_d with the densities @c noise.process on w and
 *   @c noise.disturbance on w_d (Van Loan, cross terms included: the walk of d moves x within
 *   the period), so G on [x; d] is I.
 * - MeasuredDisturbance{.Bd}: a disturbance read each tick, a known input after u. The model's B
 *   is [B, B_d] (B_d continuous, ZOH'd with the plant; none on the input's delay states; on a late
 *   reading's what the plant moved under it until the reading was taken) and D is [D, 0]; the
 *   noise and so M, L and P are the design's without it, d being no noise. Predict with [u; d],
 *   the d the plant held over the period. Beside a disturbance spec the readings no longer carry
 *   the part of a load d accounts for, so d̂ estimates only the rest.
 * - InputDelay<N>{θ}: the filter is on the model discretize(sys, Ts, delay), the plant then the
 *   delay states v(k−1) … v(k−N) of the inputs still on their way. The process noise acts on
 *   the plant, not through the delay: Q_d is the plant's over the period (Van Loan) and the delay
 *   states take none, so G on the model is [I; 0]. The filter carries them exactly: predict
 *   shifts the applied input in, their covariance is zero and no reading moves them. A
 *   disturbance enters where its spec says on that model: an InputDisturbance adds to u, so it
 *   reaches the plant through the delay states; a load on the plant itself, which acts at once,
 *   is DisturbanceModel{.Bd = sys.B, .Cd = sys.D} or PlantDisturbance. A delay discretize
 *   refuses (negative, not finite or longer than N samples) is refused with its status.
 * - OutputDelay<M>{θ}: each reading arrives θ late, one θ for every output or a ColVec of one
 *   per output. The filter is on discretize(sys, Ts, [InputDelay,] OutputDelay): the plant, the
 *   input's delay states, then M delay states per output holding the readings already taken. Its
 *   y is the reading that arrives now, so update(y) takes the late reading as it is, and the
 *   estimate is the current state. The process noise is the plant's: a reading taken s into the
 *   period holds the noise up to then, C·Q(s)·Cᵀ, and shares e^{A(Ts−s)}·Q(s)·Cᵀ with the state
 *   at the period's end (Van Loan, exact; a reading taken at a sample holds none), so G on the
 *   model is I. The measurement noise is on the reading as it arrives. An InputDisturbance adds
 *   to u, which the model carries into the late readings; another spec's d is a state of the
 *   plant held over the period, so each reading holds d as it was when taken: what the plant moved
 *   under it until then and C_d·d, and in the noise its walk up to then. The model is the delayed
 *   discretization of [A, B_d; 0, 0] read through [C, C_d], d last.
 *
 * @code
 * constexpr auto kf = design::kalman(plant, 0.01, design::KalmanNoise{.process = Qc, .measurement = sigma2}, InputDelay<3>{0.022});
 * constexpr auto camera = design::kalman(plant, 0.01, design::KalmanNoise{.process = Qc, .measurement = sigma2}, OutputDelay<2>{0.015});
 * @endcode
 * @param opts   The KalmanNoise{.process, .measurement[, .disturbance]}, required: process-noise
 *               density on the NW inputs of G (positive semidefinite), per-reading measurement
 *               variance on the NV inputs of H (positive definite, no Ts scaling) and, with a
 *               disturbance spec, the density of d's random walk; a disturbance spec, an
 *               InputDelay and an OutputDelay, each or none
 * @return KalmanResult on the discretized plant, its delay states and d (G = I, or [I; 0] with an InputDelay alone), with L, M, P and status
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename... Opts>
[[nodiscard]] constexpr auto kalman(const StateSpace<NX, NU, NY, T, NW, NV>& sys, damp::type_identity_t<T> Ts, const Opts&... opts) {
    return detail::kalman_with<detail::Domain::sampled>(sys, Ts, opts...);
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

namespace detail {

/** @brief Where a KalmanFilter keeps its model and noise: in the object, any value at run time. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
struct HeldKalman {
    static constexpr bool constant = false;                                    ///< The matrices are members
    using Model = StateSpace<NX, NU, NY, T, NW, NV>;                           ///< A, B, C, D, G, H and Ts
    using ProcessNoise = Matrix<NW, NW, T>;                                    ///< Q
    using ReadingNoise = Matrix<NV, NV, T>;                                    ///< R of the model's own reading
    using ProcessFactor = Matrix<NX, detail::noise_factor_columns<NX, NW>, T>; ///< Γ, Γ Γᵀ = G Q Gᵀ (detail::process_noise_factor)
    using ReadingFactor = detail::ReadingNoiseFactor<NY, T>;                   ///< H R Hᵀ factored, and whether it is a covariance
};

/**
 * @brief A `constexpr` design::kalman's model and noise as constants of the type: the store of the
 *        KalmanFilter that Law<kRes, T>::Recursive names
 *
 * A and C are also patterns holding their values, so the products A S and C S have no term for an
 * exact zero (a delay state's) and no multiply for an exact one. The factors of G Q Gᵀ and
 * H R Hᵀ are constants, and the covariance starts at the design's steady state, factored in the
 * design's own scalar.
 *
 * @tparam Res A design::kalman result with static storage, `constexpr`
 * @tparam T   The runtime scalar
 */
template<const auto& Res, typename T>
struct ConstantKalman {
    static constexpr bool constant = true;                                                                          ///< The matrices are constants of the type
    static constexpr auto design = Res.template as<T>();                                                            ///< The design, in the runtime's scalar
    static_assert(constant_design_runs(design.success, design.status, design.sys.Ts));                              // a failed or continuous design does not compile
    static constexpr auto kA = design.sys.A;                                                                        ///< A
    static constexpr auto kB = design.sys.B;                                                                        ///< B
    static constexpr auto kC = design.sys.C;                                                                        ///< C
    static constexpr auto kD = design.sys.D;                                                                        ///< D
    static constexpr auto kGamma = process_noise_factor(Res.sys.G, Res.Q).template as<T>();                         ///< Γ, Γ Γᵀ = G Q Gᵀ: the process noise on the state
    static constexpr auto kReadingFactor = reading_noise_factor(quadratic_form(Res.sys.H, Res.R)).template as<T>(); ///< H R Hᵀ factored, the noise of the model's own reading
    static constexpr auto model = design.sys;                                                                       ///< What model() reads
    static constexpr auto process_noise = design.Q;                                                                 ///< What process_noise_covariance() reads
    static constexpr auto reading_noise = design.R;                                                                 ///< What measurement_noise_covariance() reads
    static constexpr auto initial_factor = covariance_factor_or_identity(Res.P).template as<T>();                   ///< The factor at construction: of the design's steady-state P

    /** @brief A, B, C and D, by the names the filter reads them under. */
    struct Model {
        static constexpr ConstMatrix<kA> A{}; ///< A
        static constexpr ConstMatrix<kB> B{}; ///< B
        static constexpr ConstMatrix<kC> C{}; ///< C
        static constexpr ConstMatrix<kD> D{}; ///< D
    };
    /** @brief Q is process_noise: no member. */
    struct ProcessNoise {};
    /** @brief R is reading_noise: no member. */
    struct ReadingNoise {};
    /** @brief Γ is kGamma: no member. */
    struct ProcessFactor {};
    /** @brief The reading's factor is kReadingFactor: no member. */
    struct ReadingFactor {};
};

} // namespace detail

/**
 * @brief The recursive Kalman filter: the covariance propagates each tick as its triangular
 *        factor, and the gain is formed from it
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
 * Square-root form. The filter carries the lower-triangular factor S of the covariance,
 * P = S Sᵀ, and never P itself (estimation/sqrt_covariance.hpp). predict triangularizes
 * [A S, G Q½] by Householder reflections; an update whitens the reading with the factor of
 * its R and applies Carlson's triangular update one whitened row at a time. A
 * covariance stored as P loses its small eigenvalues once cond(P) passes about 1/eps: from a
 * prior 10⁵ times a reading's variance, the float filter on P went indefinite within a few
 * ticks and refused nearly every later reading (10⁹ in double). The factor holds that range squared.
 * covariance() forms P = S Sᵀ on request; covariance_factor() is S.
 *
 * R may be singular: a reading, or a row of one, with no noise of its own is a perfect
 * measurement, and is taken at the rounding of its own prior variance, so the estimate reads
 * C x̂ = y to within eps of the innovation and the factor keeps eps times the prior variance
 * there. An update returns false, and leaves the estimate and the factor bit for bit, when the
 * reading cannot be used: an R that is not positive semidefinite or not finite, a NaN or ±Inf
 * in C, an innovation variance past the scalar's range, or a perfect reading of what is already
 * known exactly (no noise and no prior variance: C P Cᵀ + R singular).
 *
 * P₀ and Q are factored once, at construction and in reset(x0, P0), and may be singular: zeros,
 * a state known exactly, a noise input with no variance. One that is not positive semidefinite
 * or not finite does not compile where the filter is constant-evaluated; at run time P₀ becomes
 * the identity, as a failed design's, and such a Q adds no process noise.
 *
 * @code
 * KalmanFilter kf{sys, {.Q = Q, .R = R, .x0 = x0, .P0 = P0}};
 * kf.predict(u);
 * (void)kf.update(y);               // the model's own C, D and R
 * (void)kf.update(theta, C_enc, r); // a second sensor, y = C x
 * @endcode
 *
 * For a fixed (steady-state) gain, run the design on FixedGainObserver: `FixedGainObserver obs{design::kalman(sys, noise)}`
 * — about 6 % of the recursive filter's multiply-accumulates per tick at 10 states and 6 outputs,
 * no square root, no failure path. FixedGainObserver is also the estimator inside LQG, LQGI, OFLQG and
 * an MPC through a filter, with the disturbance anchor the offset-free laws use.
 *
 * A design on [x; d] (design::kalman with a disturbance) runs here as the absolute filter on the
 * whole augmented state. An input delay is in the model discretize(plant, Ts, InputDelay<N>{θ})
 * returns, which design::kalman(plant, Ts, noise, InputDelay<N>{θ}) designs on from the continuous
 * plant: its delay states v(k−1) … v(k−N) are states like any other, with zero rows in the factor.
 *
 * On a `constexpr` design::kalman, `Law<kRes, float>::Recursive` is this filter with the design's
 * model and noise as constants of the type (@p Store): no RAM for them, and the products A S and
 * C S have no term for an exact zero, a delay state's. The factor starts at the design's
 * steady state. A reading of another C, D or R (update(y, C, R)) runs dense, as with members.
 *
 * @tparam NX Filter states
 * @tparam NW Process-noise inputs of the filter's model (columns of G)
 * @tparam NV Measurement-noise inputs (columns of H)
 *
 * @see "Optimal State Estimation" (Simon, 2006), §5.1
 * @see "Discrete square root filtering: a survey of current techniques" (Kaminski, Bryson &
 *      Schmidt, 1971), IEEE TAC 16(6)
 * @see "Fast triangular formulation of the square root filter" (Carlson, 1973), AIAA Journal 11(9)
 */
template<size_t NX, size_t NU, size_t NY, typename T = float, size_t NW = NX, size_t NV = NY, typename Store = detail::HeldKalman<NX, NU, NY, T, NW, NV>>
struct KalmanFilter {
    constexpr KalmanFilter() = default;

    /**
     * @brief From design::kalman (any scalar): P∞ seeds P, and the design's Q and R are the noise.
     *
     * A design with a disturbance model runs as the absolute filter on [x; d]. A failed design
     * does not compile when constant-evaluated; at run time its P is the identity. P∞, Q and R
     * are factored in the design's scalar, then converted.
     */
    template<typename U, size_t ND>
        requires(!Store::constant)
    constexpr KalmanFilter(const design::KalmanResult<NX, NU, NY, U, NW, NV, ND>& result) // NOLINT(google-explicit-constructor): assign-from-design
        : sys(model_of(result)),
          S_(detail::covariance_factor_or_identity(result.success ? result.P : Matrix<NX, NX, U>::identity()).template as<T>()),
          Q_(result.Q.template as<T>()),
          R_(result.R.template as<T>()),
          Gamma_(detail::process_noise_factor(result.sys.G, result.Q).template as<T>()),
          Lr_(detail::reading_noise_factor(quadratic_form(result.sys.H, result.R)).template as<T>()) {}

    /**
     * @brief From a discrete model, its noise covariances and the initial estimate.
     *
     * A continuous model (Ts not positive and finite) does not compile when constant-evaluated,
     * nor does a P0 or a Q that is not positive semidefinite.
     */
    constexpr KalmanFilter(const StateSpace<NX, NU, NY, T, NW, NV>& model, const KalmanFilterConfig<NX, NU, NY, T, NW, NV>& config)
        requires(!Store::constant)
        : sys(model),
          x(config.x0),
          S_(detail::covariance_factor_or_identity(config.P0)),
          Q_(config.Q),
          R_(config.R),
          Gamma_(detail::process_noise_factor(model.G, config.Q)),
          Lr_(detail::reading_noise_factor(quadratic_form(model.H, config.R))) {
        if (!damp::finite_positive(model.Ts) && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(design::DesignStatus::invalid_sample_time); // constant-evaluated: A is not a transition matrix
        }
    }

    /** @brief The same filter in another scalar: the factors are converted, not refactored. */
    template<typename U>
        requires(!Store::constant)
    constexpr KalmanFilter(const KalmanFilter<NX, NU, NY, U, NW, NV>& other) // NOLINT(google-explicit-constructor): a scalar conversion
        : sys(other.sys.template as<T>()),
          x(other.x.template as<T>()),
          u_prev_(other.u_prev_.template as<T>()),
          S_(other.S_.template as<T>()),
          Q_(other.Q_.template as<T>()),
          R_(other.R_.template as<T>()),
          Gamma_(other.Gamma_.template as<T>()),
          Lr_(other.Lr_.template as<T>()) {}

    /** @brief Time update with the input applied now: x̂ ← A x̂ + B u, P ← A P Aᵀ + G Q Gᵀ on its factor. The updates' feedthrough takes u until the next predict. */
    constexpr void predict(const ColVec<NU, T>& u) {
        u_prev_ = u;
        x = (sys.A * x) + (sys.B * u);
        if constexpr (Store::constant) {
            detail::sqrt_time_update<kTriangularNoise>(S_, detail::as_pattern_matrix<Store::kA>, Store::kGamma);
        } else {
            detail::sqrt_time_update<kTriangularNoise>(S_, sys.A, Gamma_);
        }
    }

    /** @brief Measurement update with the model's own C, D and H R Hᵀ. False if H R Hᵀ is not positive semidefinite, C P Cᵀ + H R Hᵀ is singular, or C P Cᵀ is past the scalar's range. */
    [[nodiscard]] constexpr bool update(const ColVec<NY, T>& y) {
        const ColVec<NY, T> nu = y - ((sys.C * x) + (sys.D * u_prev_));
        if constexpr (Store::constant) {
            return record(detail::sqrt_measurement_update(S_, x, nu, detail::as_pattern_matrix<Store::kC>, Store::kReadingFactor));
        } else {
            return record(detail::sqrt_measurement_update(S_, x, nu, sys.C, Lr_));
        }
    }

    /**
     * @brief Measurement update with a linear reading y = C x + v of any size, v ~ N(0, R)
     *
     * @param y Reading (NY2 rows)
     * @param C Output matrix of the reading on the filter's state
     * @param R Covariance of the reading's noise, positive semidefinite: a zero-variance row is a perfect measurement
     * @return false if R is not positive semidefinite, C P Cᵀ + R is singular or the reading is not finite: the estimate and the covariance are unchanged
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
     * @param R Covariance of the reading's noise, positive semidefinite: a zero-variance row is a perfect measurement
     * @return false if R is not positive semidefinite, C P Cᵀ + R is singular or the reading is not finite: the estimate and the covariance are unchanged
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
     * @brief Measurement update with a linear reading y = C x + v whose C has zeros and ones known
     *        at compile time: no term for a zero and no multiply for a one in C S
     *
     * @param y Reading (C's rows)
     * @param C Output matrix of the reading on the filter's state, on its pattern
     * @param R Covariance of the reading's noise, positive semidefinite: a zero-variance row is a perfect measurement
     * @return false if R is not positive semidefinite, C P Cᵀ + R is singular or the reading is not finite: the estimate and the covariance are unchanged
     */
    template<const auto& Pattern, size_t NY2>
    [[nodiscard]] constexpr bool update(const ColVec<NY2, T>& y, const PatternMatrix<Pattern, T>& C, const damp::type_identity_t<Matrix<NY2, NY2, T>>& R) {
        return correct(ColVec<NY2, T>(y - (C * x)), C, R);
    }

    /** @brief Measurement update with a linear reading y = C x + D u + v, C on its pattern: update(y, C, R) with the held input's feedthrough. */
    template<const auto& Pattern, size_t NY2>
    [[nodiscard]] constexpr bool update(
        const ColVec<NY2, T>&                             y,
        const PatternMatrix<Pattern, T>&                  C,
        const Matrix<NY2, NU, T>&                         D,
        const damp::type_identity_t<Matrix<NY2, NY2, T>>& R
    ) {
        return correct(ColVec<NY2, T>(y - ((C * x) + (D * u_prev_))), C, R);
    }

    /**
     * @brief Fused tick: predict with @p u, correct with @p y from the model's own C, D and R,
     *        return the state.
     *
     * One call is one complete filter tick (StateEstimator). A reading update(y) refuses skips
     * the correction and the prediction stands; last_status() says so.
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

    /** @brief What the last measurement update did: taken, or refused (update returned false, step() went on with the prediction). */
    [[nodiscard]] constexpr UpdateStatus last_status() const { return last_status_; }

    /** @brief The running error covariance P = S Sᵀ, formed from the factor on each call: exactly symmetric, rounded to eps·‖P‖. */
    [[nodiscard]] constexpr Matrix<NX, NX, T> covariance() const { return detail::covariance_of_factor(S_); }

    /** @brief The factor S the filter carries, P = S Sᵀ: lower triangular with a non-negative diagonal, and a zero row for a state known exactly. */
    [[nodiscard]] constexpr const Matrix<NX, NX, T>& covariance_factor() const { return S_; }

    /** @brief The filter's model. */
    [[nodiscard]] constexpr const StateSpace<NX, NU, NY, T, NW, NV>& model() const {
        if constexpr (Store::constant) {
            return Store::model;
        } else {
            return sys;
        }
    }

    /** @brief Process-noise covariance Q. */
    [[nodiscard]] constexpr const Matrix<NW, NW, T>& process_noise_covariance() const {
        if constexpr (Store::constant) {
            return Store::process_noise;
        } else {
            return Q_;
        }
    }

    /** @brief Measurement-noise covariance R of the model's own reading, on the inputs of H. */
    [[nodiscard]] constexpr const Matrix<NV, NV, T>& measurement_noise_covariance() const {
        if constexpr (Store::constant) {
            return Store::reading_noise;
        } else {
            return R_;
        }
    }

    /** @brief Overwrite the estimate between updates (a clamp to a physical bound, an angle wrap); P and the held input are kept. */
    constexpr void set_state(const ColVec<NX, T>& x_new) { x = x_new; }

    /** @brief Forget the applied input and the last status; the estimate and its covariance are left in place. */
    constexpr void reset() {
        u_prev_ = ColVec<NU, T>{};
        last_status_ = UpdateStatus::ok;
    }

    /**
     * @brief Replace the estimate and its covariance, and forget the applied input. Q and R stay.
     *
     * @p P0 is factored here and may be singular. One that is not positive semidefinite or not
     * finite does not compile when constant-evaluated; at run time the covariance becomes the
     * identity.
     */
    constexpr void reset(const ColVec<NX, T>& x0, const Matrix<NX, NX, T>& P0) {
        x = x0;
        S_ = detail::covariance_factor_or_identity(P0);
        u_prev_ = ColVec<NU, T>{};
        last_status_ = UpdateStatus::ok;
    }

private:
    template<size_t, size_t, size_t, typename, size_t, size_t, typename>
    friend struct KalmanFilter;

    static constexpr bool kTriangularNoise = detail::noise_factor_is_triangular<NX, NW>; ///< Whether Γ is lower triangular (NW ≥ NX) or the dense G Q½

    /** @brief A design's model in the filter's scalar, or no constant expression when the design failed: the first member built, so the refusal carries the design's status. */
    template<typename Result>
    [[nodiscard]] static constexpr StateSpace<NX, NU, NY, T, NW, NV> model_of(const Result& result) {
        if (!result.success && damp::is_constant_evaluated()) {
            detail::refuse_failed_design(result.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
        return result.sys.template as<T>();
    }

    /** @brief The factor at construction: of the design's steady state on constants, else zeros until a constructor sets it. */
    [[nodiscard]] static constexpr Matrix<NX, NX, T> initial_factor() {
        if constexpr (Store::constant) {
            return Store::initial_factor;
        } else {
            return {};
        }
    }

    /** @brief The update with innovation ν, output matrix C (dense or on its pattern) and reading covariance R, factored here. */
    template<size_t NYS, typename CMatrix>
    [[nodiscard]] constexpr bool correct(const ColVec<NYS, T>& nu, const CMatrix& C, const Matrix<NYS, NYS, T>& R) {
        return record(detail::sqrt_measurement_update(S_, x, nu, C, detail::reading_noise_factor(R))); // an R that is not a covariance: fail closed
    }

    /** @brief Keep an update's outcome for last_status(), and pass it on. */
    constexpr bool record(bool taken) {
        last_status_ = taken ? UpdateStatus::ok : UpdateStatus::refused;
        return taken;
    }

    [[no_unique_address]] typename Store::Model         sys{};                          ///< The filter's model
    ColVec<NX, T>                                       x{};                            ///< Estimate
    ColVec<NU, T>                                       u_prev_{};                      ///< The input the plant holds, for the feedthrough
    Matrix<NX, NX, T>                                   S_{initial_factor()};           ///< Lower-triangular factor of the error covariance, P = S Sᵀ
    [[no_unique_address]] typename Store::ProcessNoise  Q_{};                           ///< Process noise covariance
    [[no_unique_address]] typename Store::ReadingNoise  R_{};                           ///< Measurement noise covariance of the model's own reading
    [[no_unique_address]] typename Store::ProcessFactor Gamma_{};                       ///< Γ, Γ Γᵀ = G Q Gᵀ, factored once
    [[no_unique_address]] typename Store::ReadingFactor Lr_{};                          ///< H R Hᵀ factored once; not valid when it is not a covariance, and update(y) then refuses
    UpdateStatus                                        last_status_{UpdateStatus::ok}; ///< What the last measurement update did
};

/** @brief Deduce the filter from its design, keeping the result's scalar: KalmanFilter kf{res.as<float>()}. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t ND>
KalmanFilter(const design::KalmanResult<NX, NU, NY, T, NW, NV, ND>&) -> KalmanFilter<NX, NU, NY, T, NW, NV>;

/** @brief Deduce the filter from its model: KalmanFilter kf{sys, {.Q = Q, .R = R}}. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
KalmanFilter(const StateSpace<NX, NU, NY, T, NW, NV>&, const KalmanFilterConfig<NX, NU, NY, T, NW, NV>&) -> KalmanFilter<NX, NU, NY, T, NW, NV>;

} // namespace damp
