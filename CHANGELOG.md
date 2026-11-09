# Changelog

All notable changes to this project are documented here.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project uses [Semantic Versioning](https://semver.org/) for tags.

0.x releases may break API without a deprecation cycle (see [AGENTS.md](AGENTS.md)
and [known limitations](docs/known_limitations.md)).

## [0.2.0] — 2026-11-08

One runtime and one set of tick verbs for every LQ and MPC design, input and
output delay in the plant model, and a core that builds with no C++ standard
headers, with the numerics, examples and compiler workarounds that go with them.
The API breaks from 0.1.0 throughout: see Changed and Removed.

### Added

- [docs/known_issues.md](docs/known_issues.md): compiler bugs that change results. They are GCC 15+ folding `pow(0, y)` to 0 under `-ffast-math`, the MinGW AVX stack-alignment fault (GCC PR 54412), GCC 15.1/15.2 vectorizer wrong code, and what fused multiply-adds do to exact comparisons.
- `design::oflqg` / `OFLQG` and `design::ofmpc` / `OFMPC`, each in its own header,
  over one shared state-feedback core. `Law<kRes, T>` names a law from its design
  result for a struct member.
- `TrackedOutputs` so LQR and LQG can `control(r)` / `step(r, m)` through the
  model's steady-state target.
- `design::output_inverse` / `OutputInverse` and `design::state_feedforward`.
- `InputDelay<N>{θ}` on `discretize`: θ ≤ N·Ts in N registers, NX + N·NU states,
  exact under ZOH; `InputDelay{θ}` is N = 1. Every law and estimator takes the model
  as it is. Tustin pre-warp: `discretize(sys, Ts, Tustin, w0)`.
- `OutputDelay<N>{θ}` on `discretize`: a reading that arrives θ after it was taken, held in N
  output registers, NX + N·NY states, exact under ZOH for a fraction of a sample too. A filter on
  the model estimates the current state from the late readings, and MPC counts the dead time.
  `discretize(sys, Ts, InputDelay<N>{θu}, OutputDelay<M>{θy})` takes both, in that state order.
- One delay per channel: `InputDelay<N, ColVec<NU>>{{θ₁, θ₂}}` and `OutputDelay<N, ColVec<NY>>`,
  a delay for each input or output in the same N registers; the scalar form is still one delay for
  every channel.
- `Law<kRes, T>` of a `constexpr` design runs on `ConstMatrix`: the same ticks with the matrices
  constants of the type, so they take no RAM and a coefficient that is exactly 0 is no code. It covers
  LQR, LQI, LQG, LQGI, OFLQG, `FixedGainObserver`, `ReducedOrderObserver`, `ESO`, `ADRCController`,
  `StateFeedforward`, `OutputInverse`, and OFMPC's filter (its QP data stays in members). Nothing can
  replace them: a member a redesign assigns is `Law<kRes, T>::Runtime`. A design initialized at
  start-up keeps members. `Law` now also covers `design::kalman`, `design::place_observer`,
  `design::reduced_order_observer`, `design::eso` and `design::adrc`.
- `Law<kRes, T>::Recursive` of a `design::kalman`: the recursive `KalmanFilter` of the same model and
  noise. On a `constexpr` design its model and noise are constants of the type, so `A P Aᵀ` and the
  model's own update have no term for an exact zero and no multiply for an exact one: a 2-state plant
  with five delay registers takes 1,018 estimated cycles a tick on a Cortex-M4F against 8,266 with
  members, and 228 B against 568 B. `KalmanFilter` takes a trailing store parameter for it.
- `SparsityPattern` and `PatternMatrix` (`matrix/pattern_matrix.hpp`): a matrix whose zeros and ones are
  known at compile time and whose other entries are runtime values, for a Jacobian whose numbers change
  every tick and whose shape does not. A product has no term for a zero and no multiply for a one, and
  `quadratic_form_in_place` computes `F S Fᵀ` with no copy of S. `ErrorStateKalmanFilter` takes one in
  `predict(F)`, `update_innovation(ν, H, R)` and `reset_error_state(G)`, and the INS and the attitude
  filter use them throughout. On a Cortex-M4F the 15-state `F P Fᵀ` takes 3,906 estimated cycles and
  a 3-row update 8,331, against 37,618 and 51,463 for the dense ones.
- `KalmanFilter::update(y, C, R)` and `update(y, C, D, R)` take a `PatternMatrix` C, and the EKF's
  measurement model may return its C = ∂h/∂x as one from `jacobian(x, u)` (`PatternMeasurement`).
  On a Cortex-M4F, in float, a 7-state `KalmanFilter` reading one state takes 764 estimated cycles
  against 3,914 dense, and a 10-state EKF's 2-row update 2,922 against 15,734.
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
- The README's "Plant → law → tick" snippet is compiled as `checks/readme_front_door.cpp`, on the same Cortex-M7 build as the guide's.
- The ETL backend takes any copy of ETL on the include path; `libs/etl` is only what
  this repo's tests use. README has the steps.
- `design::eso` / `ESO`: the Han–Gao extended-state observer on its own
  (ADRC is this observer plus the PD law that cancels f̂). The discrete current
  estimator in innovation form: `predict` advances with Φ, Γ, `update` corrects
  with L. Any finite b₀. `reset(z0)` on `ESO` and `ADRCController`.

### Changed

- `mat::EigenResult::values` is a `damp::array<damp::complex<T>, N>`, as the pole lists,
  `GeneralizedEigenResult::values` and the transmission zeros already are; it was a `ColVec`.
  An eigenvalue list is a list, not a vector: index it with `values[i]`.
- The examples and the guide's canonical snippet deploy as `Law<kRes, float>`, on the design's constants: cart-pole, pendulum,
  LQI, ADRC and the heater OFMPC sketches and SILs, and the reaction-wheel SIL. Their tick functions take `Law<kRes, T>&`.
- README and the guide ("Where the design runs") show the three deploy forms side by side: `Law<kRes, float>` on
  constants, a compile-time design with its gains in RAM, and a design run on the MCU at init.
- The fixed-gain observer family is renamed for what it is: `Luenberger` → `FixedGainObserver`,
  `ReducedLuenberger` → `ReducedOrderObserver`, `design::luenberger` → `design::place_observer`
  (`LuenbergerResult` → `PlaceObserverResult`), `design::reduced_luenberger` →
  `design::reduced_order_observer` (`ReducedLuenbergerResult` → `ReducedOrderObserverResult`), and
  `estimation/luenberger.hpp` → `estimation/fixed_gain_observer.hpp`. `FixedGainObserver` runs a fixed
  gain from `design::kalman` or `design::place_observer`; `KalmanFilter` is the recursive filter.
- Matrix operations no longer zero-fill a result whose every entry they write, triangular solves
  build their result inside the `optional`, and an operator on a plain `Matrix` returns its result
  with no copy. That removes `memset` and `memcpy` calls from the run-time path. On a Cortex-M4F
  with full newlib a dense 15-state `F P Fᵀ` goes from 41,567 to 37,618 estimated cycles and the
  dense 3-row ESKF update from 61,809 to 51,463; with newlib-nano, whose `memset` works a byte at a
  time, that update goes from 146,340 to 65,501. A triangular solve checks its whole diagonal first.
- `Matrix` `+`, `-` and `+=`, `-=` keep constant indices at run time up to 64 entries, so a local
  past GCC's 16-iteration unroll limit (a `ColVec<18>`) still scalar-replaces and stays out of memory.
- One tick vocabulary: estimators `predict` / `update`, controllers `control`,
  `step` for a whole tick, `back_calculate(u_sat, u_ff = {})` for the applied
  input. Ts lives on the design; runtimes never take it.
- `LQRWeights` template order is NX, NU, NI (0 for a regulator), T:
  `LQRWeights<2, 1, 1>` is LQI on double.
- `design::lqr` follows the plant's time domain (CARE, or ZOH + Van Loan + DARE
  with a Ts). `PIDGains` is the PID result; one `pidtune` and one
  `pid_pole_placement`. `design::kalman` takes a disturbance spec.
  `FixedGainObserver` is the fixed-gain observer; `KalmanFilter` is recursive only.
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
  (`invalid_parameter`) unless a terminal cost sees past it or an `Ru` weight steers
  that input to its target; it used to succeed and leave that input at rest. The test
  is relative to the data, so units and the model's basis do not change it, and
  `suggest_mpc_horizon(sys)` counts the dead time into the horizon it suggests.
- Numerics: Riccati certificates, SciPy oracles, Schur / Lyapunov / QR, finite
  inputs give finite outputs.
- The embeddable core no longer needs any C++ standard header on the ETL backend:
  add `ETL_NO_STL` to `etl_profile.h`. Traits,
  limits and concepts go through `damp::` aliases, and `make nostd-check` builds and
  runs a constexpr design and its tick with `-nostdinc++`.
- With no `damp_profile.hpp` and no backend macro, Damp builds on the host defaults without a warning, so `g++ -std=c++20 -I inc` is quiet. The warning stays for a backend macro set with no profile (`-DDAMP_BACKEND_ETL`, `DAMP_MATH_BACKEND_*`), because every file in the program must agree on the backend.
- Compilers: g++ 11+ and clang++ 19+. MSVC is not supported. Host CMake builds
  use `-mtune=native`. Every product sketch must fit clang's default constexpr
  budget (`tools/clang_budget_check.sh`).

### Removed

- `ins_propagate_covariance` and `ins_F_mul`: `ErrorStateKalmanFilter::predict(F)` on the sparse
  transition does their job and costs less. It also carries a δx not yet injected through F, as the
  generic `predict` always did; after an injection δx is zero, so nothing changes there.
  `ins_error_jacobian` returns that `PatternMatrix`, and `ins_inject` the sparse reset Jacobian;
  `.to_matrix()` gives the dense one.
- `design/synthesis.hpp` (`*_bundle` wrappers). `StateFeedback` as a public
  runtime. `offset_free_mpc.hpp`: `OFMPC` replaces it.
- `toolbox/encoder.hpp` (now `damp/encoder/`). `toolbox/thermal.hpp`.
- `estimation/commissioning_frontends.hpp`, `estimation/experiment_safety.hpp`,
  `estimation/successive_compensator.hpp`, `filters/highpass.hpp` (`HighPass`).

### Fixed

- GCC 14 no longer crashes (an internal compiler error in `build_data_member_initialization`) on a `constexpr` constructor holding a local whose member's default initializer copies a `Matrix`, as an `RlsConfig` does through its `Bounds`. The copy constructor reached `unset_elements`, which is not `constexpr`; it now builds its unset array itself, with the same code.
- `UnscentedKalmanFilter::NSIG` is `NSP`. MinGW's `<signal.h>` defines `NSIG` as a macro, so GCC 11 on Windows could not compile `control.hpp`.
- Results no longer depend on the compiler rounding a cancellation to exactly 0. GCC 15 under `-ffast-math`, and any build that fuses multiply-adds (Cortex-M4F/M7, x86-64-v3), broke these:
  - `design::eso` decides "ω_o·T_s below eps of T" on the product. GCC 15 folded the old test `1 − e^{−ω_o T_s} > 0` into `ω_o T_s > 0`, so a design whose gain was exactly 0 reported success.
  - The implicit integrators (Backward Euler, BDF2, TR-BDF2, trapezoidal) refuse an I − hA that is singular to working precision. 1 − 0.1·10 is −5.6e-17 when fused, and the solve accepted that as a pivot.
  - `design::lowpass` and the other transfer-function stability checks judge float coefficients in double. Float's own recursion passed a pole at |z| = 1.011 when its multiply-adds were fused.
  - `MHE` clamps a converged QP solution onto its box, as it already did one that ran out of iterations. A bound of 0 could read −1e-17.
- `lyap` and `dlyap` return an exactly symmetric X when Q is symmetric. Before, the gap between the two triangles grew with the conditioning.
- `design::place_jordan` on a defective structure no longer fails under GCC 15.2 `-O3 -ffast-math` on x86-64, which vectorizes the chain vectors' open-coded complex multiply-accumulate wrongly. They are `Matrix` products now.
- GCC 15 warnings are gone: `-Wuninitialized` on the `Matrix` copy constructor's seed, which now has one return so it is built in place, and `-Wmaybe-uninitialized` on the pole-placement BFGS matrix.
- The PLLs treat an amplitude below the smallest normal value as no signal. A lost input could stick in the subnormals and hold the phase error at full scale, so `integrator_leak` never bled the offset back to nominal on a CPU that keeps subnormals. The tests missed it because linking with `-ffast-math` flushes them on x86.

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

[0.2.0]: https://github.com/Wetmelon/damp/releases/tag/v0.2.0
[0.1.0]: https://github.com/Wetmelon/damp/releases/tag/v0.1.0
