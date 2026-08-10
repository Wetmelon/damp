// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file parameter_estimation.hpp
 * @brief Online grey-box parameter estimation: recover named physical
 *        parameters (gain, time constant, inertia, damping) from streaming
 *        input/output data, with plausibility and confidence gating.
 *
 * The library's system-identification stance is grey-box first: the model
 * structure is known physics with unknown parameters, and the estimator's job
 * is to report those parameters in physical units with an honest statement of
 * whether they can be trusted yet. This is the pattern the motor calibrator
 * (`motor/calibration.hpp` — R/L from a PRBS experiment) already uses,
 * generalized: a linear-in-parameters recursive regression
 * (`estimation/rls.hpp`) on the ZOH-discretized model,
 * followed by an exact algebraic extraction back to the physical parameters.
 *
 * The two grey-box mechanisms, and where each lives:
 *
 * 1. Recursive regression + physical extraction (this header): for
 *    parameters that enter a linear-in-parameters regression — plant gains,
 *    time constants, inertia/damping, winding R/L. Excite (see
 *    `estimation/excitation.hpp`), stream (u, y), read parameters when
 *    valid() && confidence() clear your gate.
 * 2. Parameter-as-state augmentation: for slowly-varying offsets observed
 *    through known dynamics — load torque (`motor/mechanical_estimator.hpp`),
 *    gyro bias (`estimation/eskf.hpp`), input disturbances
 *    (`controllers/mpc.hpp`). Augment the state, let the Kalman
 *    machinery estimate it; no new tooling needed.
 *
 * The adaptation-policy pieces (`ParameterDriftMonitor`, `AdaptationGate`)
 * close the loop from estimate to action: detect that the plant has moved,
 * and decide whether a proposed parameter update is safe to apply to whatever
 * consumes it (a retuned gain, `PRController::set_frequency`,
 * `Biquad::set_coefficients`, a controller re-synthesis). Gated updates are
 * the default posture; continuous adaptation is the gate held open.
 *
 * Everything here is fixed-size, allocation-free, and constexpr-capable.
 *
 * @note Compare with the MATLAB® Recursive Polynomial Model Estimator refit to
 *       physical parameterizations, and Simulink®'s inertia-estimation blocks.
 * @see estimation/rls.hpp for the regression engine,
 *      estimation/excitation.hpp for test signals, concepts.hpp for the
 *      ParameterEstimator concept this models
 * @see Ljung, "System Identification: Theory for the User," 2nd ed., 1999,
 *      ch. 11 (recursive estimation); Åström & Wittenmark, "Adaptive
 *      Control," 2nd ed., 1995, ch. 2 (real-time parameter estimation)
 *
 * @code
 * // First-order plant (e.g. velocity axis: tau*omega' + omega = K*tau_em):
 * design::FirstOrderPlantEstimatorConfig<float> cfg{.Ts = 1e-3f};
 * FirstOrderPlantEstimator estimator{cfg};
 * // each tick: estimator.update(u, y);
 * // when estimator.valid() && estimator.confidence() is high enough:
 * //   K = estimator.gain(); tau = estimator.time_constant();
 * // Mechanics map (example only — not a separate type): for J*w' + b*w = tau,
 * //   b = 1/K,  J = tau/K  (require K > 0).
 * @endcode
 *
 * @see examples/estimation/fo_plant_inertia/ for the J/b recipe
 */

#include <cstddef>
#include <limits>
#include <type_traits>

#include "damp/estimation/rls.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"

namespace damp {

namespace design {

/**
 * @struct FirstOrderPlantEstimatorConfig
 * @brief Configuration for the first-order grey-box estimator
 */
template<typename T = double>
struct FirstOrderPlantEstimatorConfig {
    T Ts{};                                    ///< Sample time [s] (required > 0)
    T forgetting{static_cast<T>(0.999)};       ///< RLS forgetting factor λ ∈ (0, 1]
    T initial_covariance{static_cast<T>(100)}; ///< RLS P₀ (large = weak prior)
    T confidence_scale{static_cast<T>(1e-3)};  ///< Covariance trace at which confidence() reads 0.5

    [[nodiscard]] constexpr bool valid() const {
        return Ts > T{0} && forgetting > T{0} && forgetting <= T{1} && initial_covariance > T{0} && confidence_scale > T{0};
    }
};

} // namespace design

/**
 * @brief Online estimator for a first-order plant's gain and time constant
 *
 * Identifies @f$ \tau \dot y + y = K u @f$ from streaming (u, y) samples.
 * Internally: recursive least squares on the exact ZOH-discrete form
 * @f$ y_k = a\, y_{k-1} + b\, u_{k-1} @f$, then the algebraic extraction
 * @f[
 *   \tau = -T_s / \ln a, \qquad K = b / (1 - a)
 * @f]
 * (the motor calibrator's R/L math, generalized to any first-order physics:
 * velocity loops, thermal paths, hydraulic flow, RC stages).
 *
 * valid() is the physical-plausibility gate: the RLS must have data and the
 * estimated pole must be strictly inside (0, 1) — a stable plant with a finite
 * time constant. confidence() maps the RLS covariance trace through
 * scale/(scale + trace): near 0 before excitation, toward 1 as persistent
 * excitation shrinks the covariance. It is a monotone gating heuristic, not a
 * probability.
 *
 * The plant must be excited to identify: a constant input at steady state
 * carries no information (the covariance stays large and confidence() stays
 * low — which is exactly the honest answer). Pair with a PRBS or chirp from
 * estimation/excitation.hpp during commissioning.
 *
 * @see design::FirstOrderPlantEstimatorConfig, ParameterEstimator (concepts.hpp)
 */
template<typename T = float>
    requires std::is_floating_point_v<T>
class FirstOrderPlantEstimator {
public:
    constexpr FirstOrderPlantEstimator() = default;

    constexpr explicit FirstOrderPlantEstimator(const design::FirstOrderPlantEstimatorConfig<T>& config)
        : rls_(estimation::RlsConfig<T>{config.forgetting, config.initial_covariance}),
          Ts_(config.Ts),
          confidence_scale_(config.confidence_scale),
          config_valid_(config.valid()) {}

    /**
     * @brief Feed one sample pair: the applied input and the resulting measurement
     *
     * @param u Input applied this tick
     * @param y Output measured this tick
     * @return false if the regression update was rejected (degenerate regressor)
     */
    constexpr bool update(T u, T y) {
        bool accepted = true;
        if (has_prev_) {
            const ColVec<2, T> phi{y_prev_, u_prev_};
            accepted = rls_.update(phi, y);
        }
        u_prev_ = u;
        y_prev_ = y;
        has_prev_ = true;
        return accepted;
    }

    /// Discrete pole estimate a (y_k = a·y_{k-1} + b·u_{k-1}); raw RLS theta, not PE-gated
    [[nodiscard]] constexpr T pole() const { return rls_.state().theta(0); }

    /**
     * @brief Steady-state gain K = b / (1 − a) [output units per input unit]
     *
     * Returns 0 when @ref valid is false (avoids /0 or non-physical K from a
     * pole at 1 or a non-initialized RLS). Callers must still gate on
     * @ref valid and @ref confidence before deploying K.
     */
    [[nodiscard]] constexpr T gain() const {
        if (!valid()) {
            return T{0};
        }
        const T a = rls_.state().theta(0);
        const T b = rls_.state().theta(1);
        return b / (T{1} - a);
    }

    /**
     * @brief Time constant τ = −Ts / ln(a) [s]
     *
     * Returns 0 when @ref valid is false (avoids log of non-positive a).
     */
    [[nodiscard]] constexpr T time_constant() const {
        if (!valid()) {
            return T{0};
        }
        return -Ts_ / damp::log(rls_.state().theta(0));
    }

    /// {K, τ} — the ParameterEstimator concept readout (physical units; zeros if !valid)
    [[nodiscard]] constexpr ColVec<2, T> parameters() const { return ColVec<2, T>{gain(), time_constant()}; }

    /**
     * @brief Physical-plausibility gate: initialized RLS and a stable pole in (0, 1)
     *
     * This is not a PE certificate. Flat or non-exciting (u, y) can still leave
     * @c valid true with a weak prior after a lucky first sample — check
     * @ref confidence (and optionally an external PE monitor) before trusting K/τ.
     */
    [[nodiscard]] constexpr bool valid() const {
        constexpr T edge = static_cast<T>(1e-6);
        const T     a = rls_.state().theta(0);
        return config_valid_ && rls_.valid() && rls_.state().initialized && a > edge && a < T{1} - edge;
    }

    /**
     * @brief Monotone excitation/convergence indicator in [0, 1): scale/(scale + covariance trace)
     *
     * Near 0 without PE (covariance stays large); rises only as the regressor
     * shrinks P. Heuristic gating score, not a calibrated probability.
     */
    [[nodiscard]] constexpr T confidence() const {
        if (!rls_.state().initialized) {
            return T{0};
        }
        const T trace = rls_.state().covariance(0, 0) + rls_.state().covariance(1, 1);
        return confidence_scale_ / (confidence_scale_ + trace);
    }

    /// Latest one-step prediction residual (drives ParameterDriftMonitor)
    [[nodiscard]] constexpr T residual() const { return rls_.state().residual; }

    constexpr void reset() {
        rls_.reset();
        u_prev_ = T{0};
        y_prev_ = T{0};
        has_prev_ = false;
    }

private:
    estimation::RlsVector<2, T> rls_{};
    T                           Ts_{};
    T                           confidence_scale_{static_cast<T>(1e-3)};
    T                           u_prev_{};
    T                           y_prev_{};
    bool                        has_prev_{false};
    bool                        config_valid_{false};
};

/**
 * @brief Residual-based drift detector: has the plant moved away from the model?
 *
 * Tracks the estimator's one-step prediction residual with two exponential
 * moving averages — a fast one (recent behaviour) and a slow one (the learned
 * baseline). drift_score() is their ratio: ≈1 while the model explains the
 * data, rising when the plant changes faster than the estimator has adapted.
 * drift_detected() fires when the score crosses the threshold after the
 * baseline has had time to form.
 *
 * Feed it the residual() of a parameter estimator (or a Kalman innovation)
 * every tick. Typical use: trigger a re-identification campaign or open an
 * AdaptationGate review when it fires.
 */
template<typename T = float>
    requires std::is_floating_point_v<T>
class ParameterDriftMonitor {
public:
    constexpr ParameterDriftMonitor() = default;

    /**
     * @param fast_alpha  EWMA weight of the recent-residual tracker (e.g. 0.05)
     * @param slow_alpha  EWMA weight of the baseline tracker (e.g. 0.002; must be < fast_alpha)
     * @param threshold   drift_score() level that flags drift (e.g. 3)
     */
    constexpr ParameterDriftMonitor(T fast_alpha, T slow_alpha, T threshold)
        : fast_alpha_(fast_alpha), slow_alpha_(slow_alpha), threshold_(threshold) {}

    /// Feed one residual sample; returns drift_detected() after the update
    constexpr bool update(T residual) {
        const T mag = damp::abs(residual);
        fast_ += fast_alpha_ * (mag - fast_);
        slow_ += slow_alpha_ * (mag - slow_);
        ++samples_;
        return drift_detected();
    }

    /// Recent-to-baseline residual ratio (≈1 = model still explains the data)
    [[nodiscard]] constexpr T drift_score() const {
        constexpr T floor_eps = std::numeric_limits<T>::epsilon() * T{100};
        return fast_ / ((slow_ > floor_eps) ? slow_ : floor_eps);
    }

    /// True once the baseline exists (≥ 2/slow_alpha samples) and the score crosses the threshold
    [[nodiscard]] constexpr bool drift_detected() const {
        const auto warmup = static_cast<size_t>(T{2} / slow_alpha_);
        return samples_ >= warmup && drift_score() > threshold_;
    }

    constexpr void reset() {
        fast_ = T{0};
        slow_ = T{0};
        samples_ = 0;
    }

private:
    T      fast_alpha_{static_cast<T>(0.05)};
    T      slow_alpha_{static_cast<T>(0.002)};
    T      threshold_{static_cast<T>(3)};
    T      fast_{};
    T      slow_{};
    size_t samples_{0};
};

/**
 * @brief Gated-adaptation policy: decide whether a parameter update is safe to apply
 *
 * The seam between estimation and action. A proposed parameter set is allowed
 * through only when the estimator vouches for it (valid), the confidence
 * clears the configured floor, and the step from the last accepted set is
 * within the per-parameter bound (no confident-but-wild jumps into the
 * consumer). Accepted sets are remembered as the next comparison baseline.
 *
 * Continuous adaptation is this gate with a zero confidence floor and
 * unbounded steps; gated (commissioning-style) adaptation is the default
 * posture per decisions.md O2.
 */
template<size_t NPARAM, typename T = float>
    requires std::is_floating_point_v<T>
class AdaptationGate {
public:
    /// Outcome of one evaluation
    struct Decision {
        T    confidence{};   ///< The confidence the decision was based on
        bool allow_update{}; ///< true: apply the proposed parameters
    };

    constexpr AdaptationGate() = default;

    /**
     * @param min_confidence  Confidence floor an update must clear (0 = always)
     * @param max_step        Per-parameter |Δ| bound vs. the last accepted set
     *                        (numeric_limits max sentinel = unbounded)
     */
    constexpr AdaptationGate(T min_confidence, const ColVec<NPARAM, T>& max_step)
        : min_confidence_(min_confidence), max_step_(max_step) {}

    /// Evaluate a proposed parameter set; remembers it as baseline when allowed
    constexpr Decision evaluate(bool estimator_valid, T confidence, const ColVec<NPARAM, T>& proposed) {
        Decision decision{confidence, false};
        if (!estimator_valid || confidence < min_confidence_) {
            return decision;
        }
        if (has_accepted_) {
            for (size_t i = 0; i < NPARAM; ++i) {
                if (damp::abs(proposed(i) - accepted_(i)) > max_step_(i)) {
                    return decision;
                }
            }
        }
        accepted_ = proposed;
        has_accepted_ = true;
        decision.allow_update = true;
        return decision;
    }

    /// Last parameter set that passed the gate
    [[nodiscard]] constexpr const ColVec<NPARAM, T>& accepted_parameters() const { return accepted_; }
    [[nodiscard]] constexpr bool                     has_accepted() const { return has_accepted_; }

    constexpr void reset() {
        accepted_ = ColVec<NPARAM, T>{};
        has_accepted_ = false;
    }

private:
    static constexpr ColVec<NPARAM, T> unbounded_steps() {
        ColVec<NPARAM, T> v{};
        for (size_t i = 0; i < NPARAM; ++i) {
            v(i) = std::numeric_limits<T>::max();
        }
        return v;
    }

    T                 min_confidence_{T{0}};
    ColVec<NPARAM, T> max_step_ = unbounded_steps();
    ColVec<NPARAM, T> accepted_{};
    bool              has_accepted_{false};
};

} // namespace damp
