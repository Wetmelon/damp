# Changelog

All notable changes to this project are documented here.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project uses [Semantic Versioning](https://semver.org/) for tags.

0.x releases may break API without a deprecation cycle (see [AGENTS.md](AGENTS.md)
and [known limitations](docs/known_limitations.md)).

## [Unreleased]

### Added

- `design::discrete_oflqg` / `OFLQG`: offset-free LQG (discrete LQR, input-disturbance Kalman filter, DC target).
- `discretize` takes an input and output delay; every LQ / MPC / estimator design accepts plant-sized weights on delayed models.
- Tustin pre-warping: `discretize(sys, Ts, Tustin, w0)` (and `c2d`) matches the discrete frequency response exactly at ω₀.
- `mat::ldlt` / `mat::ldlt_solve`: LDLᵀ without pivoting for symmetric indefinite matrices (quasi-definite KKT included).
- `TransferFunction::to_observable_canonical_form` (and `to_controllable_canonical_form` as the named twin of `to_state_space`).
- LQI / LQGI integrate a chosen subset of outputs; saturation-aware `feedback` / `commit` across LQI, LQGI, MPC and OffsetFreeMPC.

### Changed

- `OffsetFreeMPC` takes `OffsetFreeModel` and its noise in OFLQG's order (`Q_process`, `R_measurement`, `Q_disturbance`); `mpc` likewise takes `R` before the input weight.
- Numerics: DARE convergence relative to the solution, Cholesky refuses pivots below the roundoff of their terms, `compute_eigenvalues` balances its argument, `is_stabilizable` equilibrates before PBH, `minreal` uses the staircase controllable subspace.
- NaN/Inf audit: finite inputs and failed designs give finite output across controllers, filters and toolbox.
- Host builds use `-mtune=native` instead of `-march=native` (g++ 14 MinGW misaligned AVX-512 stack stores).
- Prefer prvalue returns when assembling `StateSpace` / design results and geometry
  values (style; same behavior under optimization).
- Matrix: `ColVec`/`RowVec` `segment`, richer `RowView`/`ColView` assign and indexing;
  INS sparse covariance uses column/row views.
- Target smoke (`make targets`): full green PlatformIO board matrix, quiet `pio -s`,
  product flags `-O3 -ffast-math` (no LTO), framework warning demotion for Arduino cores.

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
