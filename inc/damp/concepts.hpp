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
 * | SisoController           | `u = step(r, y)`, `back_calculate(u_sat, u_ff)` | PID, ADRC, SMC, lead-lag, PR      |
 * | SisoLaw                  | SisoController + `limits()`/`set_limits()`, `track(u, r, y)`, `last_command()`, `valid()` | PID (tests/test_contract.cpp lists the others' gaps) |
 * | FeedbackRuntime          | `update(m)`, `predict()`, `back_calculate(u)` | LQR, LQI, LQG, LQGI, OFLQG, MPC, OFMPC (lifecycle) |
 * | FeedbackLaw              | FeedbackRuntime + `control(r)`, `step(r, m)`, `track(u, r, m)` | LQI, LQGI, OFLQG, OFMPC, LQR/LQG/MPC with TrackedOutputs |
 * | StateEstimator           | `x̂ = step(y, u)`, `state()`       | FixedGainObserver (fixed gain), KalmanFilter (recursive), MHE, ESO |
 * | SignalSource             | `u = step()`, `done()`            | Chirp, PRBS, StepTrain, Ramp, MultiSine, SteppedSine, Impulse, Step |
 * | ParameterEstimator       | `parameters()`, `valid()`,        | FirstOrderPlantEstimator          |
 * |                          | `confidence()`                    |                                   |
 *
 * Naming convention the concepts encode: estimators say `predict` (time update)
 * and `update` (measurement update), an LQ or MPC law says `control` for its command, and `step`
 * combines one runtime's per-tick calls into one tick. `step(y, u)` is the fused
 * estimator tick (time update + measurement update + readout), so generic code
 * cannot silently skip the prediction step; the explicit `predict`/`update`
 * split remains available for multirate use. A SISO law (PID, PR, lead-lag, SMC, ADRC, …) has
 * no split: `step(r, y)` is its whole tick. Sources and experiments say `step`.
 *
 * Every LQ and MPC runtime (LQR, LQI, LQG, LQGI, OFLQG, MPC, OFMPC) has the same
 * verbs. `update(m)` takes the measurement at sample k (the state itself without
 * a filter). `control(…)` takes only the references and returns u_k.
 * `step(…, m)` is predict, then update, then control. The low-latency ISR order
 * is update → control → write → back_calculate → predict: the trailing predict
 * is the next tick's time update. For a measured-state LQR, predict is empty, so
 * `step(x)` equals update then control. `back_calculate(u_sat, u_ff)` records the
 * input the plant received; without it the next predict assumes the last command
 * went out unchanged.
 *
 * Deliberately outside the taxonomy (documented, pinned by the conformance
 * tests in tests/test_concepts.cpp):
 * - Experiments (RelayAutotuner, PhaseParameterCalibrator): controllers
 *   with a lifecycle and a product; they get a concept when a generic
 *   controller-switch consumer exists.
 * - EKF/UKF/ESKF: predict(u) and one update per reading (update(h, y, R),
 *   update_innovation(ν, H, R)), with no fused step(y, u): one filter fuses
 *   readings of several sizes, so there is no single y.
 *
 * Nonlinear models have one shape, written once: Dynamics ẋ = f(x, u) and
 * Measurement y = h(x, u). The same objects drive design::linearize, the
 * ExtendedKalmanFilter, the UnscentedKalmanFilter and the sim::simulate* plants.
 * A model may carry its own Jacobians (AnalyticDynamics, AnalyticMeasurement);
 * without them the consumers take central differences.
 *
 * Conformance is enforced by static_asserts in tests/test_concepts.cpp: any
 * runtime API drift becomes a compile error there.
 */

#include "damp/backend.hpp"
#include "damp/matrix/matrix.hpp"

namespace damp {

/**
 * @brief Scalar output-feedback controller: u = step(r, y), fixed rate.
 *
 * One call is one complete tick. Sample time comes from the design result the
 * runtime was built from, never a constructor or step() argument; this is the
 * common denominator that dt-baked laws can satisfy, so it is the concept contract.
 *
 * A 3-arg variable-rate step(r, y, dt) overload stays *outside* the concept
 * and is allowed only where the sample time enters the law in closed form at runtime
 * (PIDController::step(r, y, dt): `I += Ki*e*dt` with I the integral term, the filtered
 * derivative's coefficients from one reciprocal 1/(dt + Tf)), so a measured dt is exact and
 * cheap. Fixed-rate controllers bake Ts at design time (PIDController's step(r, y) at
 * design::pid's Ts, Kalman F/Q, LQG, PR) and must NOT expose a runtime Ts: that would be a lie
 * unless the coeffs are re-discretized every tick. Multirate = several blocks each storing their own fixed Ts, not
 * one block driven at two rates.
 *
 * `back_calculate(u_sat, u_ff)` reports the input the plant received after a clamp
 * downstream, and @p u_ff what was added to this law's command before it. Every law
 * has it, a no-op where nothing winds up (P, PD, SMC), so a composite (Cascade,
 * SwitchedController, SmithPredictor) passes a clamp on to the law it holds.
 */
template<typename C, typename T>
concept SisoController = damp::is_floating_point_v<T> && requires(C c, T r, T y, T u_sat, T u_ff) {
    { c.step(r, y) } -> damp::convertible_to<T>;
    c.back_calculate(u_sat, u_ff);
    c.reset();
};

/**
 * @brief A law whose run-time limits are read with limits() and set with set_limits(lim)
 *
 * set_limits refuses limits that are not a box (an end NaN or ±Inf, a lower end above its upper
 * one): it returns false and changes nothing. A public limit member can't refuse, which is why every
 * law holds its limits behind the pair.
 */
template<typename C>
concept HoldsLimits = requires(C c, const C& held) {
    held.limits();
    { c.set_limits(held.limits()) } -> damp::convertible_to<bool>;
};

/**
 * @brief A SISO law that hands over without a bump: track(u, r, y) in place of step(r, y) while
 *        another source drives the plant, u the input it applies
 *
 * No mode flag: the law tracks for as long as the caller calls track, and the next step(r, y) at
 * that r and y starts from u. A law with nothing to preload keeps its memory current.
 */
template<typename C, typename T>
concept SisoTracking = requires(C c, T u, T r, T y) { c.track(u, r, y); };

/** @brief A law that reports the command its last tick returned: last_command(). */
template<typename C, typename T>
concept ReportsCommand = requires(const C& c) {
    { c.last_command() } -> damp::convertible_to<T>;
};

/** @brief A law that says whether it runs: valid() is false when it was built from a design it refuses. */
template<typename C>
concept ReportsValidity = requires(const C& c) {
    { c.valid() } -> damp::convertible_to<bool>;
};

/**
 * @brief The whole SISO runtime contract: a SisoController that holds its limits, tracks, reports
 *        its command and says whether it runs
 *
 * tests/test_contract.cpp checks every SISO runtime against it, and lists each gap a law still has.
 */
template<typename C, typename T>
concept SisoLaw = SisoController<C, T> && HoldsLimits<C> && SisoTracking<C, T> && ReportsCommand<C, T> && ReportsValidity<C>;

/**
 * @brief Shared lifecycle of an LQ or MPC runtime: measurement, time update, applied input.
 *
 * This is the state-management surface, not the full law. control() and step() take different
 * references per design (none for a regulator, r for tracking / integral / MPC), so they are
 * FeedbackLaw. NM is the measurement: y with a filter, x without one.
 *
 * Discrete-time index: m is the reading at sample k. predict() uses the input the plant
 * held over [k−1, k). The low-latency order update → control → write → back_calculate →
 * predict is the same cycle with the time update at the end of the ISR.
 */
template<typename C, size_t NU, size_t NM, typename T>
concept FeedbackRuntime = requires(C c, const ColVec<NM, T>& m, const ColVec<NU, T>& u) {
    c.update(m);
    c.predict();
    c.back_calculate(u);
    c.reset();
};

/**
 * @brief An LQ or MPC law that holds outputs at r: control(r), step(r, m), and track(u, r, m) in
 *        place of step while another source drives the plant (u the input applied).
 *
 * LQI, LQGI, OFLQG, OFMPC, and LQR/LQG/MPC with TrackedOutputs. A regulator (LQR/LQG/MPC
 * with no tracked outputs) is FeedbackRuntime only: control(), step(m), track(u, m).
 */
template<typename C, size_t NU, size_t NR, size_t NM, typename T>
concept FeedbackLaw = FeedbackRuntime<C, NU, NM, T> && requires(C c, const ColVec<NU, T>& u, const ColVec<NR, T>& r, const ColVec<NM, T>& m) {
    { c.control(r) } -> damp::convertible_to<ColVec<NU, T>>;
    { c.step(r, m) } -> damp::convertible_to<ColVec<NU, T>>;
    c.track(u, r, m);
};

/**
 * @brief State estimator: x̂ = step(y, u) — the fused per-tick form.
 *
 * step() runs the complete filter tick (time update with the applied input
 * u, measurement update with y) and returns the estimate; state() reads it
 * back without advancing.
 */
template<typename E, size_t NX, size_t NU, size_t NY, typename T>
concept StateEstimator = requires(E e, const ColVec<NY, T>& y, const ColVec<NU, T>& u) {
    { e.step(y, u) } -> damp::convertible_to<ColVec<NX, T>>;
    { e.state() } -> damp::convertible_to<ColVec<NX, T>>;
    e.reset();
};

/**
 * @brief Self-clocked signal source: u = step(), finished when done().
 *
 * The excitation generators' shape; also the natural fit for trajectory
 * playback and test stimuli.
 */
template<typename S, typename T>
concept SignalSource = damp::is_floating_point_v<T> && requires(S s) {
    { s.step() } -> damp::convertible_to<T>;
    { s.done() } -> damp::convertible_to<bool>;
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
concept ParameterEstimator = damp::is_floating_point_v<T> && requires(E e) {
    { e.parameters() } -> damp::convertible_to<ColVec<NPARAM, T>>;
    { e.valid() } -> damp::convertible_to<bool>;
    { e.confidence() } -> damp::convertible_to<T>;
    e.reset();
};

/**
 * @brief Continuous-time nonlinear dynamics ẋ = f(x, u): a callable returning ẋ as a ColVec<NX, T>.
 *
 * A lambda, a function or a struct with operator(). The same object is the model of
 * design::linearize (A and B at an operating point), of ExtendedKalmanFilter and
 * UnscentedKalmanFilter (each holds it, reads NX and NU from its call operator, and integrates
 * it over its sample time with one RK4 step) and the
 * plant of the continuous sim::simulate* functions. Time does not enter: a time-varying
 * model carries t as a state or captures it.
 *
 * A model may also provide its Jacobians through a member jacobian(x, u) returning
 * DynamicsJacobian {A = ∂f/∂x, B = ∂f/∂u} (AnalyticDynamics). Every consumer uses them when
 * present and central differences otherwise.
 *
 * @code
 * // Damped pendulum about the upright, θ from vertical: θ̈ = (g/l)·sin θ − c·θ̇ + b·u
 * struct Pendulum {
 *     double g_over_l = 9.81, c = 0.5, b = 1.0;
 *     ColVec<2> operator()(const ColVec<2>& x, const ColVec<1>& u) const {
 *         return {x[1], (g_over_l * damp::sin(x[0])) - (c * x[1]) + (b * u[0])};
 *     }
 *     DynamicsJacobian<2, 1> jacobian(const ColVec<2>& x, const ColVec<1>&) const { // optional
 *         return {.A = {{0.0, 1.0}, {g_over_l * damp::cos(x[0]), -c}}, .B = {{0.0}, {b}}};
 *     }
 * };
 * static_assert(AnalyticDynamics<Pendulum, 2, 1, double>);
 * @endcode
 */
template<typename F, size_t NX, size_t NU, typename T>
concept Dynamics = requires(const F& f, const ColVec<NX, T>& x, const ColVec<NU, T>& u) {
    { f(x, u) } -> damp::convertible_to<ColVec<NX, T>>;
};

/**
 * @brief Nonlinear measurement y = h(x, u): a callable returning y as a ColVec<NY, T>.
 *
 * The sensor model of design::linearize (C and D) and of the EKF and UKF measurement updates.
 * A model may provide its Jacobians through a member jacobian(x, u) returning
 * MeasurementJacobian {C = ∂h/∂x, D = ∂h/∂u} (AnalyticMeasurement); otherwise the consumers
 * take central differences.
 */
template<typename H, size_t NX, size_t NU, size_t NY, typename T>
concept Measurement = requires(const H& h, const ColVec<NX, T>& x, const ColVec<NU, T>& u) {
    { h(x, u) } -> damp::convertible_to<ColVec<NY, T>>;
};

/** @brief Analytic Jacobians of a Dynamics model at (x, u), returned by its jacobian(x, u). */
template<size_t NX, size_t NU, typename T = double>
struct DynamicsJacobian {
    Matrix<NX, NX, T> A{}; ///< ∂f/∂x
    Matrix<NX, NU, T> B{}; ///< ∂f/∂u
};

/** @brief Analytic Jacobians of a Measurement model at (x, u), returned by its jacobian(x, u). */
template<size_t NX, size_t NU, size_t NY, typename T = double>
struct MeasurementJacobian {
    Matrix<NY, NX, T> C{}; ///< ∂h/∂x
    Matrix<NY, NU, T> D{}; ///< ∂h/∂u
};

/**
 * @brief Dynamics that provide their own Jacobians: f.jacobian(x, u) → DynamicsJacobian.
 *
 * A model whose jacobian(x, u) returns anything else is a compile error in the consumers,
 * not a silent fall-back to central differences.
 */
template<typename F, size_t NX, size_t NU, typename T>
concept AnalyticDynamics = Dynamics<F, NX, NU, T> && requires(const F& f, const ColVec<NX, T>& x, const ColVec<NU, T>& u) {
    { f.jacobian(x, u) } -> damp::convertible_to<DynamicsJacobian<NX, NU, T>>;
};

namespace detail {

/** @brief An R × C PatternMatrix of scalar T: what PatternMeasurement's jacobian(x, u) returns. */
template<typename M, size_t R, size_t C, typename T>
concept pattern_matrix_of = is_pattern_matrix<M> && damp::is_same_v<typename M::Scalar, T> && (M::NR == R) && (M::NC == C);

} // namespace detail

/**
 * @brief A Measurement whose jacobian(x, u) returns C = ∂h/∂x as an NY × NX PatternMatrix: zeros
 *        and ones known at compile time, the other entries computed at (x, u)
 *
 * The EKF's update then has no term for a zero and no multiply for a one. A consumer that needs
 * D = ∂h/∂u (design::linearize) takes it by central differences.
 */
template<typename H, size_t NX, size_t NU, size_t NY, typename T>
concept PatternMeasurement = Measurement<H, NX, NU, NY, T> && requires(const H& h, const ColVec<NX, T>& x, const ColVec<NU, T>& u) {
    requires detail::pattern_matrix_of<damp::remove_cvref_t<decltype(h.jacobian(x, u))>, NY, NX, T>;
};

/** @brief A Measurement that provides its own Jacobians: h.jacobian(x, u) → MeasurementJacobian. */
template<typename H, size_t NX, size_t NU, size_t NY, typename T>
concept AnalyticMeasurement = Measurement<H, NX, NU, NY, T> && requires(const H& h, const ColVec<NX, T>& x, const ColVec<NU, T>& u) {
    { h.jacobian(x, u) } -> damp::convertible_to<MeasurementJacobian<NX, NU, NY, T>>;
};

} // namespace damp
