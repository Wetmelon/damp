# Changelog

All notable changes to this project are documented here.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project uses [Semantic Versioning](https://semver.org/) for tags.

0.x releases may break API without a deprecation cycle (see [AGENTS.md](AGENTS.md)
and [known limitations](docs/known_limitations.md)).

## [Unreleased]

Core API wave: one runtime per LQ/MPC design, one set of tick verbs, delay in
the plant, and the matching numerics and heater examples.

This section describes the tree at HEAD. It replaces the unreleased delay-in-the-law
work that landed on main after 0.1.0 (`InputDelayPredictor`, `PredictedInputDelay`,
delay templates on OFLQG / MPC). v0.1.0 never shipped those names.

### Added

- `design::oflqg` / `OFLQG` and `design::ofmpc` / `OFMPC`, each in its own header,
  over one shared state-feedback core. `Law<kRes, T>` names a law from its design
  result for a struct member.
- `TrackedOutputs` so LQR and LQG can `control(r)` / `step(r, m)` through the
  model's steady-state target.
- `design::output_inverse` / `OutputInverse` and `design::state_feedforward`.
- `InputDelay<N>{θ}` on `discretize`: θ ≤ N·Ts in N registers, NX + N·NU states,
  exact under ZOH; `InputDelay{θ}` is N = 1. Every law and estimator takes the model
  as it is. Tustin pre-warp: `discretize(sys, Ts, Tustin, w0)`.
- `DiscretizationResult` from `discretize` / `c2d` / `d2d`. `lft`. Canonical
  forms in `systems/realization.hpp`. Transmission zeros (QZ) and `norm_hinf` /
  disk margins on the host.
- `matrix/concat.hpp`: size-deducing `diag` / `blkdiag` / `hcat` / `vcat`.
  `mat::ldlt`. Encoders under `damp/encoder/`.
- `lyapunov_cholesky` / `dlyapchol` (Hammarling). SOS IIR cascade.
  `sim::quantize`.
- Examples: `heater_ofmpc` and `heater_lqg` on a 5 s hotend, the dead time one
  sample of the 5 s control period.
- [docs/guide.md](docs/guide.md): task table from plant to law, plus MATLAB® /
  firmware / Python-control name maps (`c2d`, `pidtune`, `back_calculate`, ISR).
  Canonical example [`examples/control/lqi/`](examples/control/lqi/). CMake
  build-tree alias `wetmelon::damp`.
- `FeedbackLaw`: `control(r)` and `step(r, m)` on top of the `FeedbackRuntime`
  lifecycle (`update` / `predict` / `back_calculate`).
- Host contract check for `workbench.hpp` alone (`make workbench-check`). The
  guide's canonical snippet is compiled as `checks/guide_canonical.cpp`.
- The ETL backend takes any copy of ETL on the include path; `libs/etl` is only what
  this repo's tests use. README has the steps.
- `design::eso` / `ESO`: the Han–Gao extended-state observer on its own
  (ADRC is this observer plus the PD law that cancels f̂). The discrete current
  estimator in innovation form: `predict` advances with Φ, Γ, `update` corrects
  with L. Any finite b₀. `reset(z0)` on `ESO` and `ADRCController`.

### Changed

- One tick vocabulary: estimators `predict` / `update`, controllers `control`,
  `step` for a whole tick, `back_calculate(u_sat, u_ff = {})` for the applied
  input. Ts lives on the design; runtimes never take it.
- `LQRWeights` template order is NX, NU, NI (0 for a regulator), T:
  `LQRWeights<2, 1, 1>` is LQI on double.
- `design::lqr` follows the plant's time domain (CARE, or ZOH + Van Loan + DARE
  with a Ts). `PIDGains` is the PID result; one `pidtune` and one
  `pid_pole_placement`. `design::kalman` takes a disturbance spec.
  `Luenberger` is the fixed-gain observer; `KalmanFilter` is recursive only.
- Every design reports `success` / `DesignStatus` and `explicit operator bool()`.
  A refused `constinit` law names the cause. `ok` is that design's checks
  (weights, convergence, a closed-loop pole test where the design runs one), not
  a robustness certificate. The Riccati residual lives on the solver result.
- MPC and MHE ship in `control.hpp` (allocation-free). Weights are plant-sized.
  An infeasible MPC tick holds the last applied input (still
  inside the u box). A MaxIterations or Degenerate iterate is applied when the
  outputs are soft. With a hard output row or a terminal set it is applied only
  if primal feasible (row scale |b_i| + Σ|A_ij z_j|); otherwise the tick holds.
  Holding does not make predicted outputs or a terminal set feasible.
  A horizon too short for an input's first move to reach any output is refused
  (`invalid_parameter`) unless a terminal cost sees past it; it used to succeed and
  leave that input at rest.
- Numerics: Riccati certificates, SciPy oracles, Schur / Lyapunov / QR, finite
  inputs give finite outputs.
- Compilers: g++ 11+ and clang++ 19+. MSVC is not supported. Host CMake builds
  use `-mtune=native`. Every product sketch must fit clang's default constexpr
  budget (`tools/clang_budget_check.sh`).

### Removed

- `design/synthesis.hpp` (`*_bundle` wrappers). `StateFeedback` as a public
  runtime. `discrete_oflqg`, `offset_free_lqg.hpp`, `offset_free_mpc.hpp`.
- Delay in the design or runtime (`InputDelayPredictor`, `PredictedInputDelay`,
  `OutputDelay`, delay template parameters on OFLQG / MPC).
- `toolbox/encoder.hpp` (now `damp/encoder/`). `toolbox/thermal.hpp`.
- `estimation/commissioning_frontends.hpp`, `estimation/experiment_safety.hpp`,
  `estimation/successive_compensator.hpp`, `filters/highpass.hpp` (`HighPass`).

## [0.1.0] — 2026-08-09

First public release of Damp: header-only C++20 control design that ships in the
same tree as firmware (design-is-deploy: same types in SIL and on the MCU).

License: Boost Software License 1.0. Third-party host-only components are listed
in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

### Entry points

| Include | Role |
| ------- | ---- |
| `damp/control.hpp` | Embeddable core: matrix, systems, `design::`, core design-is-deploy controllers/estimators (incl. ESKF), filters, fixed-step integrators, toolbox |
| `damp/workbench.hpp` | Host: control plus Bode/margins, ODE solvers, closed-loop `simulate`, `matlab::` short names |

Drop `inc/` into an existing firmware tree, or use CMake (`find_package` / FetchContent).

### Library

- Matrix / math: fixed-size stack `Matrix`, views, decompositions, solve; pluggable math backend; geometry (quaternion / DCM)
- LTI systems: `StateSpace`, transfer functions, ZPK, discretization (ZOH / Tustin / Euler), interconnections
- Design (not PWM-rate): place, Riccati / DARE, PID rules, minreal / model reduction, stability, synthesis bundles → Result + `.as<U>()`
- Core design-is-deploy controllers: PID, PR, LQR / LQI / LQG / LQGI, lead-lag, ADRC, SMC / STSMC, Smith predictor, action governor
- Core design-is-deploy estimators: Kalman, EKF, UKF, ESKF + attitude fusion, Luenberger, DOB, RLS
- Filters: biquad family and common signal blocks; robust exact differentiator
- Simulation: fixed-step integrators on-target; host ODE solvers and closed-loop simulate via workbench
- Host analysis: Bode, Nyquist, margins, step/impulse/lsim, poles; MATLAB®-style aliases

Runtime controllers and estimators default to `float`; design synthesis defaults to `double`.

### Examples and validation

- Design-is-deploy folder layout under `examples/control/` and `examples/estimation/`
- Unit tests and embedded/freestanding contracts for the core surface

### Documentation

- [README.md](README.md)
- [REFERENCE.md](REFERENCE.md)
- [docs/known_limitations.md](docs/known_limitations.md)
- [examples/](examples/)

### Known limitations

See [docs/known_limitations.md](docs/known_limitations.md).

[0.1.0]: https://github.com/Wetmelon/damp/releases/tag/v0.1.0
