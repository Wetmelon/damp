# Changelog

All notable changes to this project are documented here.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project uses [Semantic Versioning](https://semver.org/) for tags.

0.x releases may break API without a deprecation cycle (see [AGENTS.md](AGENTS.md)
and [known limitations](docs/known_limitations.md)).

## [Unreleased]

### Changed

- `LQRConfig`, `LQIConfig`, `MPCWeights` and `MPCConstraints` take their sizes from the fields a named one sets, with one input, integrator and controlled output unless a field says otherwise: `LQIConfig{.Q = diag(10.0, 1.0), .R = {1.0}, .Q_integral = {80.0}}` is `LQIConfig<2, 1, 1>`. A single weight is braced; the positional guides that kept `{.Q, .R, .demand_lim}` from deducing are gone.
- A Kalman design (`kalman`, `lqg`, `lqgi`, `oflqg`, `mpc` and `ofmpc` with noise, `mhe`) reads a plant that declares no noise inputs (`StateSpace<NX, NU, NY>`) as noise on every state and reading, G = I and H = I, so its `KalmanNoise` is NX × NX and NY × NY and the plant needs no `<…, T, NW, NV>`. Before, such a plant took only a `KalmanNoise<0, 0>`.
- `PRSpec::Tt`, `PRResult::Tt` and `PRController::Tt` are `tracking_time`, the name the LQI and LQGI runtimes give the same back-calculation tracking time.
- Every example SIL checks what it prints and returns nonzero when a result leaves its band, so `make` fails on a change in behavior, not only on a crash. Before, 2 of 18 did.

### Fixed

- The FOC SIL reported a 98% rise time one tick short (0 ms for the PI), and the LPF SIL called its first-order filter second-order Butterworth. The workflow SIL ran 2 s, before its 1 Hz PR mode settled; it runs 20 s and ends at r.

## [0.3.1] — 2026-11-17

The LQ integrators hold against a clamp by default again, as a PID's integral does.

### Changed

- The LQI and LQGI integrators hold against a clamp (conditional integration) when `tracking_time` is 0, the default, as a PID's integral does with `Kbc = 0`. A `tracking_time` above 0 unwinds them by back-calculation, and one up to Ts takes the whole excess each tick, 0.3.0's default: set `tracking_time = Ts` to keep 0.3.0's behaviour. Near its limit with noise on the reading, that default took each excursion that reached the clamp back from the integrators, from one side only, and an LQI settled short of r by about ten times the offset of holding; a slow integrator gave back for seconds what a short clamp took from it. Back-calculation still recovers from a long clamp without holding's overshoot.
- An LQ design whose (∂u/∂ξ)⁺ cannot be formed is refused only when it asks for back-calculation.

## [0.3.0] — 2026-11-17

Designs for a plant whose input arrives late, one limits and anti-windup contract for every law (with state limits on MPC and a bound on the offset-free laws' disturbance estimate), a Kalman filter and an MPC that hold in float, no fault at run time, and Damp after `Arduino.h`. The API breaks: `LQRWeights` is now `LQRConfig`, or `LQIConfig` with integral action, and a `SisoController` needs `back_calculate`.

### Added

- `MPCConstraints` takes `x_states` and `x_lim`, with their softness: an MPC or OFMPC bounds states it does not track over the horizon, soft or hard, and `set_limits` moves them. The type names how many, `MPCConstraints<NU, NY, T, NXC>`.
- `last_status()` on `KalmanFilter` and `ExtendedKalmanFilter`: whether the last measurement update was taken or refused (`UpdateStatus`), so the refusal `step()` drops is still there to read. MPC's `last_status()` reports its QP the same way.
- `OFMPC::disturbance_lim` clamps the disturbance estimate d̂ at each update, as on `OFLQG`, and `estimator()` has a mutable overload to restart the filter or set its estimate.
- `to_ss()` on an OFLQG design: the law as a discrete block at r = 0, its filter on [x; d] and the target's d̂ terms folded into the gain, to close on a plant other than the design's (one with another delay, say) for its poles.
- `LQRConfig` and `LQIConfig` take `demand_lim`, a `DemandLimit`: the LQ laws bound what the gain's outer terms ask of one state, read as a cascade (a balancing robot's position and integral asking for a lean), and the integrators unwind against it. Off by default, and free in a `constexpr` design's `Law`.
- `OFLQG::disturbance_lim` clamps the disturbance estimate d̂ at each update, so a load past it (a wall the plant is held against) does not wind d̂ and the target up, the way `integral_lim` bounds an LQI's integrator. Unbounded by default.
- `design::lqr`, `lqi`, `lqg`, `lqgi`, `oflqg`, `kalman` and `steady_state_target` take an `InputDelay` after Ts: the continuous plant, its weights and its noise, the cost and the noise discretized exactly on the delayed model. A measured-state law holds the registers itself, so `update` takes the plant state; a filter carries them with no noise on them.
- `LQRConfig`, the weights of an LQ design, takes `u_lim`, and `LQIConfig` also `integral_lim`: the LQR, LQI, LQG, LQGI and OFLQG runtimes clamp their command and their integral state to them, and unwind the integrators against that clamp with no `back_calculate`. Both are members of the runtime (`law.u_lim`, `law.integral.lim`), so a limit can move at run time.
- `Bounds::valid()`: no end NaN and no lower end above its upper one, the test `design::mpc` and `design::mhe` already applied.
- `u_lim` on `ADRCResult`, `STSMCResult` and `LeadLagResult`, and as a member of their runtimes and of `ExtremumSeekingController`: the command is clamped into it and the law's memory holds against that clamp.
- `back_calculate(u_sat, u_ff)` on `STSMCController`, `LeadLagController`, `ExtremumSeekingController`, `SmithPredictor`, `RepetitiveController` and `ClassicalDOB`: a clamp after the law is reported the same way to every law with memory.
- `back_calculate` on every SISO law and in the `SisoController` concept: P, PD and SMC take it as a no-op, and `HarmonicSuppressor`, `Cascade` and `SwitchedController` pass a clamp on to the laws they hold.
- `u_lim` on `PRResult`, `HarmonicSuppressorResult` and `SMCResult` and as a member of their runtimes: the command is clamped into it, and the resonators of a PR law or a harmonic bank hold against the clamp.
- `valid()` on the PID runtimes, the LQ laws, `MPC` and `OFMPC`, as on the other laws: false when built from a design they refuse, and the law sends 0 (an MPC holds its input).

### Changed

- `StateFeedbackResult` and the LQR, LQI, LQG, LQGI and OFLQG runtimes take a last template argument, NRU, the input registers of the law's model (0: none); K, S and e are on [x; v(k−1); …; v(k−NRU)] then the integrators.
- `damp::min`, `damp::max` and `damp::clamp` return by value on both backends, so a result bound to a reference cannot dangle on a temporary argument. The backends now agree: a tie gives `min` and `max` their first argument (a NaN too, without `-ffast-math`), and `clamp` with lo > hi gives hi, where `std::clamp` is undefined.
- `LQRWeights<NX, NU, NI>` is `LQRConfig<NX, NU>` for `lqr`, `lqg` and `oflqg`, and `LQIConfig<NX, NU, NI>` with `Q_integral` for `lqi` and `lqgi`: each has only its law's fields, the limits with the cost. A design called with braces, `{.Q = Q, .R = R}`, does not change.
- `mat::pseudo_inverse`, `mat::null_space` and `mat::two_norm` return an optional, empty when the SVD fails, where they returned NaN; `matlab::pinv` and `matlab::null` follow. `design::steady_state_target` with more outputs than inputs refuses a failed factorization, where it could report NaN maps as a success.
- `+` and `−` on matrices and views need equal shapes: a dimension of size 1 is no longer stretched to fit (NumPy's broadcast), so a missing transpose does not compile. An offset on every entry is `A + scalar`.
- A PID with an output limit and `Kbc = 0` holds its integral against the clamp (conditional integration), where it wound up, and `back_calculate` does the same for a clamp downstream. A `Kbc` still unwinds it by back-calculation.
- `PIDGains::discretize`, the PID runtimes' constructors and `design::esc` refuse a limit that is not a box, as the LQ and MPC designs do. In `ESCSpec`, `lower ≥ upper` no longer switches the clamp off.
- `ADRCController::last_command()` is the command the law returned, as on the other laws; the input its observer takes as applied is `applied()`.
- `Bounds::as` saturates a finite end past the narrower type's range.
- `state_feedback_core.hpp` holds only the runtime the LQ laws share. `LQRConfig`, `IntegralAction`, `NoEstimator` and `StateFeedbackResult` are in `state_feedback_types.hpp`, and the design behind them in `controllers/detail/`; `lqr.hpp`, `lqi.hpp`, `lqg.hpp`, `lqgi.hpp` and `oflqg.hpp` include all of it as before.
- `solve_qp` at its iteration budget: a solve that its last allowed iteration completes is `Success`, and a `MaxIterations` point is the minimizer on the constraints taken in so far.
- A float MPC riding a hard output limit at full input can read `Infeasible` on a tick whose rounded state is past the limit's reach, and hold its input. It returned `Success` there with the limit's row violated.
- `KalmanFilter` and `ExtendedKalmanFilter` carry the triangular factor S of their covariance, P = S Sᵀ: `covariance()` returns P by value, formed on each call, and `covariance_factor()` is S.
- A perfect measurement on `KalmanFilter` or `ExtendedKalmanFilter` (R = 0, or a zero-variance row of R) is taken at the rounding of its own prior variance: the estimate reads it to eps, and the covariance keeps eps times the prior there instead of zero.
- A P₀ or a Q of `KalmanFilter` or `ExtendedKalmanFilter` may be singular, and one that is not positive semidefinite or not finite does not compile where the filter is constant-evaluated; at run time P₀ becomes the identity and such a Q adds no process noise.
- The update of `KalmanFilter` and `ExtendedKalmanFilter` takes far less stack, about a ninth for a large filter; a `KalmanFilter` that holds its matrices is larger by the factors of its Q and R.
- Damp never faults at run time. A PID runtime built from gains it refuses is invalid where its constructor trapped, `DiscretizationResult::value()` of a refused result is the zero model, and a matrix initializer longer than the matrix drops the extra entries; each is still a compile error under constant evaluation.
- `Cascade` holds the outer law while the inner law's command is clamped, at its `u_lim` or by a clamp reported with `back_calculate`, so an outer integrator no longer winds up while the inner loop cannot follow.
- `SwitchedController` starts a law with a tracking mode, a PID, from the command the plant last received when it is switched in, where it started from its reset state.
- `PRController` with `Tt = 0` holds its resonant section against a clamp reported with `back_calculate`, where it ignored the report.
- The integrators' limits of the LQI and LQGI runtimes and their design results are `integral_lim`, the name `LQIConfig` gives them, where they were `integral.lim`.
- The LQI and LQGI integrators unwind against a clamp by back-calculation, where they held (conditional integration): a law recovers from a long clamp without the overshoot. `tracking_time` on `LQIConfig` and the runtimes slows it; the default takes the whole excess back each tick.

### Fixed

- Damp's headers compile after `Arduino.h`, the order an `.ino` sketch gets: Damp undefines Arduino's `abs` macro, as `<cstdlib>` does, and no longer names anything its `F`, `PI` or binary-constant (`B1`) macros rewrite.
- A float MPC with soft output limits solves its QP on every tick. Once the moves sat on their limits the active-set solver ended `Degenerate`, or called a feasible problem `Infeasible`, on up to every tick of a recovery from outside the limits; it now keeps an orthogonal factorization of the active constraints.
- `solve_qp` in float no longer returns `Success` with a constraint violated past its tolerance: the point is recomputed on the active constraints after each step, and every solved point meets its rows.
- PID back-calculation with a `Kbc` under half of `Ki·Ts` diverged, and overshot under `Ki·Ts`; it takes at most the whole excess per tick.
- An LQ design refused for its limits no longer carries them into a runtime, and an `integral_lim` narrowed at run time pulls the integrator in on the next tick.
- `PRController<float>` and `HarmonicSuppressor<float>` hold a resonance far below the sample rate, where the output used to lose its gain at ω₀ or run to Inf (f₀/fs = 5e-5). The resonant section, `PRController::resonant`, is a `PRResonator` on the output and its increment instead of a `Biquad`.
- `KalmanFilter` and `ExtendedKalmanFilter` in float keep their covariance under a prior far wider than a reading's noise: stored as P it went indefinite within a few ticks once P₀/R passed about 10⁵, and nearly every later update was refused.
- `PRController` with a tracking time `Tt` under the sample period overshot the realizable command on each unwind, and diverged under half of it; it takes at most the whole excess per tick.

## [0.2.2] — 2026-11-10

A fix for graded matrices, which 0.2.1's balancing change broke.

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

[0.3.1]: https://github.com/Wetmelon/damp/releases/tag/v0.3.1
[0.3.0]: https://github.com/Wetmelon/damp/releases/tag/v0.3.0
[0.2.2]: https://github.com/Wetmelon/damp/releases/tag/v0.2.2
[0.2.1]: https://github.com/Wetmelon/damp/releases/tag/v0.2.1
[0.2.0]: https://github.com/Wetmelon/damp/releases/tag/v0.2.0
[0.1.0]: https://github.com/Wetmelon/damp/releases/tag/v0.1.0
