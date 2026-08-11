// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file concepts.hpp
 * @brief The library-wide runtime taxonomy: C++20 concepts pinning the per-tick
 *        vocabulary every controller, estimator, and signal source follows.
 *
 * First-principles roles, by signal signature (reference r, measurement y,
 * command u, state x — reference-first argument order everywhere):
 *
 * | Role                     | Per-tick shape                    | Examples                          |
 * |--------------------------|-----------------------------------|-----------------------------------|
 * | SisoController           | `u = control(r, y)`               | PID, ADRC, SMC, lead-lag, PR      |
 * | OutputFeedbackController | `u = control(r, y)` (vector)      | OffsetFreeMPC, LQGI               |
 * | StateFeedbackController  | `u = control(r, x)`               | LQR (r = x_ref), MPC (r = y_ref)  |
 * | StateEstimator           | `x̂ = estimate(y, u)`, `state()`   | KalmanFilter, Luenberger, MHE       |
 * | SignalSource             | `u = step()`, `done()`            | Chirp, PRBS, StepTrain, Ramp, MultiSine, SteppedSine, Impulse, Step |
 * | ParameterEstimator       | `parameters()`, `valid()`,        | FirstOrderPlantEstimator          |
 * |                          | `confidence()`                    |                                   |
 *
 * Naming convention the concepts encode: controllers say `control`, estimators
 * say `update`/`estimate`, sources and experiments say `step`. `estimate(y, u)`
 * is the fused estimator tick (time update + measurement update + readout), so
 * generic code cannot silently skip the prediction step; the explicit
 * `predict`/`update` split remains available for multirate use.
 *
 * Deliberately outside the taxonomy (documented, pinned by the conformance
 * tests in tests/test_concepts.cpp):
 * - Experiments (RelayAutotuner, PhaseParameterCalibrator): controllers
 *   with a lifecycle and a product; they get a concept when a generic
 *   controller-switch consumer exists.
 * - EKF/UKF/ESKF: their per-tick model-callback API is an intentional,
 *   different shape.
 * - LQG/LQGI: law + estimator *pairs*. Prefer @c step for a self-contained
 *   tick (LQG::step(y), LQGI::step(r, y)). The split predict/update/control
 *   path remains for multirate use and post-saturation commit; LQGI::control(r, y)
 *   matches OutputFeedbackController syntax but does NOT advance its Kalman
 *   filter — see the note on that method.
 *
 * Conformance is enforced by static_asserts in tests/test_concepts.cpp: any
 * runtime API drift becomes a compile error there.
 */

#include <concepts>
#include <cstddef>
#include <type_traits>

#include "damp/matrix/matrix.hpp"

namespace damp {

/**
 * @brief Scalar output-feedback controller: u = control(r, y), fixed rate.
 *
 * One call is one complete tick. Sample time is stored (constructor parameter /
 * public member), never a control() argument — this is the common denominator
 * that dt-baked laws can satisfy, so it is the concept contract.
 *
 * A 3-arg variable-rate control(r, y, Ts) overload stays *outside* the concept
 * and is allowed only where the sample time enters the law linearly at runtime
 * (ContinuousPID: `I += Ki*e*Ts` with I the integral term, `deriv = Δ/Ts`), so
 * a measured dt is exact and free. Canonical fixed-rate controllers bake Ts at
 * design time (PIDController via discretize, Kalman F/Q, LQG, PR) and must NOT
 * expose a runtime Ts: that would be a lie unless the coeffs are re-discretized
 * every tick. Multirate = several blocks each storing their own fixed Ts, not
 * one block driven at two rates.
 */
template<typename C, typename T>
concept SisoController = std::is_floating_point_v<T> && requires(C c, T r, T y) {
    { c.control(r, y) } -> std::convertible_to<T>;
    c.reset();
};

/**
 * @brief Vector output-feedback controller: u = control(r, y), self-contained tick.
 *
 * Semantic contract beyond the requires-clause: one call performs the
 * complete tick, including any internal estimator. OffsetFreeMPC satisfies
 * both syntax and semantics; LQGI satisfies only the syntax (its filter is
 * caller-sequenced).
 */
template<typename C, size_t NU, size_t NY, typename T>
concept OutputFeedbackController = requires(C c, const ColVec<NY, T>& r, const ColVec<NY, T>& y) {
    { c.control(r, y) } -> std::convertible_to<ColVec<NU, T>>;
    c.reset();
};

/**
 * @brief State-feedback law: u = control(r, x), reference-first.
 *
 * NR is the reference dimension: a state reference for LQR (NR = NX), an
 * output reference for MPC (NR = NY). Pair with a StateEstimator to build an
 * output-feedback controller (the LQG / OffsetFreeMPC pattern).
 */
template<typename C, size_t NU, size_t NR, size_t NX, typename T>
concept StateFeedbackController = requires(C c, const ColVec<NR, T>& r, const ColVec<NX, T>& x) {
    { c.control(r, x) } -> std::convertible_to<ColVec<NU, T>>;
    c.reset();
};

/**
 * @brief State estimator: x̂ = estimate(y, u) — the fused per-tick form.
 *
 * estimate() runs the complete filter tick (time update with the applied input
 * u, measurement update with y) and returns the estimate; state() reads it
 * back without advancing.
 */
template<typename E, size_t NX, size_t NU, size_t NY, typename T>
concept StateEstimator = requires(E e, const ColVec<NY, T>& y, const ColVec<NU, T>& u) {
    { e.estimate(y, u) } -> std::convertible_to<ColVec<NX, T>>;
    { e.state() } -> std::convertible_to<ColVec<NX, T>>;
    e.reset();
};

/**
 * @brief Self-clocked signal source: u = step(), finished when done().
 *
 * The excitation generators' shape; also the natural fit for trajectory
 * playback and test stimuli.
 */
template<typename S, typename T>
concept SignalSource = std::is_floating_point_v<T> && requires(S s) {
    { s.step() } -> std::convertible_to<T>;
    { s.done() } -> std::convertible_to<bool>;
    s.reset();
};

/**
 * @brief Online grey-box parameter estimator: physical parameters + gating.
 *
 * parameters() reports the identified physical quantities (gains, time
 * constants, inertias — physical units, not regression coefficients);
 * valid() is the physical-plausibility acceptance check; confidence() a
 * monotone excitation/convergence indicator in [0, 1] for adaptation gating.
 */
template<typename E, size_t NPARAM, typename T>
concept ParameterEstimator = std::is_floating_point_v<T> && requires(E e) {
    { e.parameters() } -> std::convertible_to<ColVec<NPARAM, T>>;
    { e.valid() } -> std::convertible_to<bool>;
    { e.confidence() } -> std::convertible_to<T>;
    e.reset();
};

} // namespace damp
