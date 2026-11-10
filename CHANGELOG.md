# Changelog

All notable changes to this project are documented here.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project uses [Semantic Versioning](https://semver.org/) for tags.

0.x releases may break API without a deprecation cycle (see [AGENTS.md](AGENTS.md)
and [known limitations](docs/known_limitations.md)).

## [Unreleased]

### Fixed

- Balancing takes a matrix's grading out in full again, and keeps it only where it pays, so `logm` of a graded matrix in float has the accuracy of the ungraded one (0.2.1 lost up to 16×, and failed a test on some hosts). A matrix of one scale is still left unscaled.

## [0.2.1] — 2026-11-10

Numerical robustness fixes across the matrix core, tested in float first, and Damp on the PlatformIO registry.

### Added

- `library.json`: Damp on the PlatformIO registry, `lib_deps = wetmelon/damp@^0.2.1`.

### Changed

- `mat::ldlt` pivots, P·A·Pᵀ = L·D·Lᴴ with 1×1 and 2×2 blocks, so `LDLT` gains `perm` and `e`. `ldlt_solve(A, B)` is unchanged.
- `QRDecomposition::full_rank()` is removed: it was not a rank test. Decide rank with `mat::rank`.
- An eigenvalue, Schur entry, generalized eigenvalue or singular value that overflows the type is a failure (`success = false`), not an Inf.
- For a matrix with isolated eigenvalues, `EigenResult::schur_vectors` includes the balancing permutation, and `values` can come in a different order.

### Fixed

- Balancing counts the diagonal in its row and column norms, so a matrix whose diagonal dominates is no longer over-scaled. This fixes `sqrtm_jordan_3`, which failed `make` on Intel and generic-tuned hosts. `balance_states` does the same.
- Eigenvalues and Schur forms balance with a permutation step: triangular, cascade and integrator-augmented matrices get exact eigenvalues, where a 4×4 with all-negative eigenvalues could return a positive one.
- Real Schur forms, `lyap`, `sqrt` and `log` hold from 1e-30 to 1e30 in float and from 1e-300 to 1e300 in double. A complex pair used to come back wrong past about 1e±19 in float and 1e±155 in double.
- A lower-triangular 2×2 Schur block is exchanged exactly instead of losing its small eigenvalue to cancellation.
- The QR sweep's reflector stays symmetric under `-ffast-math`, so the Schur form keeps the accuracy of a build without it.
- QZ (`generalized_eigenvalues`, `transmission_zeros`) solves its 2×2 blocks without cancellation: close eigenvalues keep their digits, and a real pair no longer comes out complex.
- QZ scales the pencil by powers of two, so ‖A‖/‖B‖ no longer has to stay within about 1e±18 in float.
- `svd` prescales by a power of two, so a matrix near the ends of the type's range no longer returns all-zero singular values as a success.
- `cholesky`, LU and the triangular solves refuse a NaN or Inf pivot under `-ffast-math`. `KalmanFilter::update` with a NaN in R returns false and leaves the state untouched.
- `ldlt` stays accurate on indefinite matrices (a float case went from 3.6e-3 relative error to 4.7e-9) and refuses NaN under `-ffast-math`.
- `state_feedforward` and the decoupling in `output_inverse` decide rank by singular values, so a rank-deficient B is refused whatever its column order.
- Norms pass a NaN or Inf entry through, so a `residual.norm() <= tol` check no longer passes on a NaN residual.

## [0.2.0] — 2026-11-09

One runtime and one set of tick verbs for every LQ and MPC design, input and output delay in the plant model, and a core that builds with no C++ standard headers. The API breaks from 0.1.0 throughout: see Changed and Removed.

### Added

- `design::oflqg` / `OFLQG` and `design::ofmpc` / `OFMPC`, output-feedback LQG and MPC, over one shared state-feedback core.
- `TrackedOutputs`: LQR and LQG can `control(r)` / `step(r, m)` through the model's steady-state target.
- `design::output_inverse` / `OutputInverse` and `design::state_feedforward`.
- `InputDelay<N>{θ}` and `OutputDelay<N>{θ}` on `discretize`: input and measurement delay in N registers, exact under ZOH, one delay per channel if needed. Tustin pre-warp: `discretize(sys, Ts, Tustin, w0)`.
- `Law<kRes, T>`: a `constexpr` design runs on its constants, so its matrices take no RAM and its exact zeros cost no code. It covers the LQ family, the observers, ESO, ADRC, `StateFeedforward`, `OutputInverse`, OFMPC's filter, and `design::kalman` through `Law<kRes, T>::Recursive`.
- `SparsityPattern` and `PatternMatrix`: zeros and ones fixed at compile time, other entries at run time, used by the ESKF, the INS and `KalmanFilter::update`. A 15-state `F P Fᵀ` on a Cortex-M4F drops from 37,618 to 3,906 estimated cycles.
- `DiscretizationResult`, `lft` and canonical forms; transmission zeros (QZ), `norm_hinf` and disk margins on the host.
- `diag` / `blkdiag` / `hcat` / `vcat`, `mat::ldlt`, `lyapunov_cholesky` / `dlyapchol`, an SOS IIR cascade, `sim::quantize`, and encoders under `damp/encoder/`.
- `design::eso` / `ESO`, the Han–Gao extended-state observer on its own, and `reset(z0)` on `ESO` and `ADRCController`.
- `FeedbackLaw`: `control(r)` and `step(r, m)` on top of the `FeedbackRuntime` lifecycle.
- Examples `heater_ofmpc` and `heater_lqg`; [docs/guide.md](docs/guide.md), from plant to law with MATLAB®, firmware and python-control names; [docs/known_issues.md](docs/known_issues.md), compiler bugs that change results.
- `make workbench-check`, and the README and guide snippets compiled as checks. The ETL backend takes any copy of ETL on the include path. CMake alias `wetmelon::damp`.

### Changed

- One tick vocabulary: estimators `predict` / `update`, controllers `control`, `step` for a whole tick, and `back_calculate(u_sat, u_ff = {})` for the applied input. Ts lives on the design.
- Every design reports `success` / `DesignStatus` and `explicit operator bool()`, a refused `constinit` law names the cause, and every `*Result` struct is `[[nodiscard]]`.
- `Luenberger` → `FixedGainObserver`, `ReducedLuenberger` → `ReducedOrderObserver`, `design::luenberger` → `design::place_observer`, with their results and header renamed to match.
- `mat::EigenResult::values` is a `damp::array<damp::complex<T>, N>`, like the other eigenvalue lists.
- `LQRWeights` takes NX, NU, NI, T (`LQRWeights<2, 1, 1>` is LQI); `design::lqr` follows the plant's time domain; `PIDGains` is the PID result; `design::kalman` takes a disturbance spec.
- MPC and MHE are allocation-free in `control.hpp`. An infeasible MPC tick holds the last applied input, and a horizon too short for an input to reach any output is refused unless a terminal cost or an `Ru` weight covers it.
- The core builds on ETL with no C++ standard headers (`ETL_NO_STL`), checked by `make nostd-check`. With no profile and no backend macro it builds on the host defaults without a warning.
- The examples and the guide deploy as `Law<kRes, float>`, and the README and the guide show the three deploy forms side by side.
- Matrix operations skip the zero-fill of results they write in full, taking `memset` and `memcpy` off the run-time path: a dense 3-row ESKF update drops from 61,809 to 51,463 estimated cycles on a Cortex-M4F.
- Compilers: GCC 11.3+ and Clang 19+; Clang on Windows targets MinGW's libstdc++. MSVC is not supported.
- Build: examples and `make target-smoke` link statically, `make` no longer needs `arm-none-eabi-g++`, and `make embedded-check` checks only Damp's own headers.

### Removed

- `CONFIG_BACKEND` and `CONFIG_COMPILER_PREFIX` in `tup.config`.
- `ins_propagate_covariance` and `ins_F_mul`: `ErrorStateKalmanFilter::predict(F)` on the sparse transition replaces them.
- `design/synthesis.hpp` (the `*_bundle` wrappers), `StateFeedback` as a public runtime, and `offset_free_mpc.hpp`, which `OFMPC` replaces.
- `toolbox/encoder.hpp` (now `damp/encoder/`), `toolbox/thermal.hpp`, `estimation/commissioning_frontends.hpp`, `estimation/experiment_safety.hpp`, `estimation/successive_compensator.hpp` and `filters/highpass.hpp`.

### Fixed

- GCC 11 and 14 no longer crash (internal compiler error) on a `constexpr` constructor that copies a `Matrix`.
- `UnscentedKalmanFilter::NSIG` is now `NSP`, since MinGW's `<signal.h>` defines `NSIG`.
- Results no longer depend on a cancellation rounding to exactly 0 under GCC 15 `-ffast-math` or fused multiply-adds: `design::eso`, the implicit integrators, the transfer-function stability checks and `MHE`.
- `lyap` and `dlyap` return an exactly symmetric X for a symmetric Q.
- `design::place_jordan` no longer fails under GCC 15.2's vectorizer, and GCC 15's warnings are gone.
- The PLLs treat a subnormal amplitude as no signal, so `integrator_leak` bleeds the offset back.
- `assoc_barrier` and `discretize` with a NaN `Ts` work on GCC 11, `design::lowpass<3>` / `<4>` keep unit DC gain in float, `math_backend.hpp` honours the profile when included first, and `dpwm_zero_sequence` is constexpr on Clang 19.
- Build: `make embedded-check` reports compile errors, `make tidy` / `iwyu` / `fix_includes` work with clang tools, CMake's `docs` target exists only at the top level, and the README's FetchContent snippet skips the submodules.

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

[Unreleased]: https://github.com/Wetmelon/damp/compare/v0.2.1...HEAD
[0.2.1]: https://github.com/Wetmelon/damp/releases/tag/v0.2.1
[0.2.0]: https://github.com/Wetmelon/damp/releases/tag/v0.2.0
[0.1.0]: https://github.com/Wetmelon/damp/releases/tag/v0.1.0
