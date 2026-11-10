# Changelog

All notable changes to this project are documented here.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project uses [Semantic Versioning](https://semver.org/) for tags.

0.x releases may break API without a deprecation cycle (see [AGENTS.md](AGENTS.md)
and [known limitations](docs/known_limitations.md)).

## [Unreleased]

### Fixed

- Balancing (`balance_exponents`, behind eigenvalues and Schur forms, zeros, stability, Lyapunov, the model inverse, `sqrt` and `log`) counts the diagonal in its row and column norms, as LAPACK's xGEBAL has since 3.5.0 (James, Langou & Lowery). Without it a matrix whose diagonal dominates was over-balanced. `sqrtm_jordan_3` took κ(D) = 16 for a 1.16× smaller norm, and mapping its root back by D cost up to 9× in accuracy. With GCC 14 on Intel or generic tuning, that was 21 units of forward error against a limit of 16, and `make` failed. It now gets D = I, at 5.2 units in its worst build, and the oracle drops its Clang 19 allowance. Graded matrices still balance. `balance_states` (`matlab::ssbal` / `prescale`) counts A's diagonal the same way, so a state whose own dynamics dominate its couplings is left unscaled. The Hamiltonian balance keeps Benner's zero diagonal.
- `cholesky` refuses a pivot that is not finite, tested on its bits (`damp::isfinite`) before its floor. Under `-ffast-math` GCC folds the floor test `!(d > floor)` into `d ≤ floor`, which a NaN fails, so a NaN anywhere in the lower triangle factored: `cholesky([[NaN, 0], [0, 1]])` returned a NaN factor and `cholesky_solve` a NaN solution, and `KalmanFilter::update(y, C, R)` with a NaN in R returned true and wrote NaN into x and P, where it documents false with both unchanged. The EKF, ESKF and UKF updates, the Joseph update and `solve_qp`'s `cholesky(H)` reach the same factorization. All of them now refuse, in float and double, and the Kalman update leaves x and P bit for bit. The IEEE build refused already; an Inf was refused by accident, through its infinite floor. The same bit test now guards LU's pivot, which a finite A near overflow turns into Inf − Inf = NaN (at `-O1` and `-O2`, `solve` and `inverse` returned NaN for `[[1, 0, 3e38], [1, 1, −3e38], [1, 2, −3e38]]` in float), and the diagonal of the triangular solves.
- `design::state_feedforward` and the decoupling matrix of `design::output_inverse` decide rank with `mat::rank` (the singular values) instead of the diagonal of an unpivoted QR, whose verdict depends on the column order. For B = [b₀, b₀ + 10⁻⁴·b₂, b₂], rank 2 with σ₃/σ₁ = 10⁻¹⁷, every |R(j, j)| cleared its column's floor and `state_feedforward` returned success with entries of B⁺ up to 2.25·10¹⁶ (1.5·10⁷ in float, where σ₃/σ₁ is 2·10⁻⁸); the same columns as [b₀, b₂, b₁] were refused. Both orders are now refused. `QRDecomposition::full_rank()` is removed: these were its only callers, and it was not a rank test. `qr_decompose` still stores a column lost to roundoff as R(j, j) = 0, so a solve on R refuses, and documents that this misses a dependency the column order hides.
- `Matrix::norm` and the other norms on `detail::scaled_norm2` (Householder vectors, balancing, the SVD's column norms, zeros, realizations) return the first magnitude that is not finite, tested on its bits, so a NaN or ±Inf entry gives a NaN or +Inf norm, as LAPACK's xNRM2 does since 3.10. `damp::max` drops a NaN against the running maximum, so `ColVec{NaN, 0}.norm()` and `{NaN, NaN}.norm()` were 0, and in float under `-ffast-math` `{1, NaN}` was finite too; a certificate `r.norm() <= tol` then passed on an all-NaN residual (`sqrtm`, the CARE basis skew check, `modal_form`). A complex entry is measured by its two parts, since `hypot` reads (NaN, 0) as 0. Under IEEE arithmetic those checks now refuse; under `-ffinite-math-only` a comparison with NaN is not reliable, and a check holds there only with `damp::isfinite` on the norm, as `sqrtm`'s has.
- `mat::ldlt` pivots: P·A·Pᵀ = L·D·Lᴴ with 1×1 and 2×2 blocks in D, chosen by the rook search of Ashcraft, Grimes & Lewis (bounded Bunch–Kaufman, α = (1 + √17)/8), so |L| ≤ 2.78 and the active entries grow at most 2.56× per step for every symmetric or Hermitian matrix. Unpivoted, the only guard was a floor for pivots that are pure roundoff, and growth went unreported: [[δ, 1], [1, 0]] (κ = 1) grew as 1/δ, so `ldlt_solve` returned x with a relative error of 3.6·10⁻³ in float at δ = 10⁻⁵ and 2.5·10⁻⁷ in double at δ = 10⁻¹⁰, and a 4×4 with κ = 4.6 lost 6.4·10⁻³ in float and 2.1·10⁻⁷ in double, all as successes. They now solve to 4.7·10⁻⁹ and 0, and 3.9·10⁻⁸ and 1.8·10⁻¹⁶, within an ulp of LU. On random 7×7 indefinite matrices the normwise backward error of the solve fell from 2600 units of n·ε in float (and 82 of 300 refused) to 0.12. Nonsingular matrices with a vanishing leading minor, such as [[0, 1], [1, 0]], now factor. A 2×2 block is solved scaled by its off-diagonal w, so neither its determinant nor |w|² is formed: a float block at 10⁻²⁵ or 10²⁵ solves to 10⁻⁷. `LDLT` gains `e`, the subdiagonal of D (nonzero on the first row of each 2×2 block), and `perm`, with (P·A·Pᵀ)(i, j) = A(perm[i], perm[j]); `d` is still D's diagonal. Under -ffast-math a NaN entry passed the pivot test (GCC folds `!(|d| > floor)`), and 14 of 30 placements of NaN or ±Inf in a 4×4 came back as factorizations; entries are now tested on their bit pattern, and an overflowing Schur complement is refused too. The doc no longer claims that a KKT matrix with E ⪰ 0 factors in any order: Vanderbei's quasi-definite result needs E ≻ 0.
- `generalized_eigenvalues` (QZ, and `transmission_zeros` through it) solves each converged 2×2 block for λ − σ, σ the diagonal ratio aᵢᵢ/bᵢᵢ nearer zero, so the discriminant is (half the gap between the two ratios)² plus a product, the shifted form LAPACK's xLAG2 credits to C. Van Loan. The unshifted quadratic det(A − λB) = 0 lost ε·λ²/gap to cancellation, and the sign of the cancelled discriminant decided real against complex. A rotation by θ with B = I came out as a real double eigenvalue at θ = 1e-4 in float and 1e-8 in double; at θ = 1e-6 the imaginary part was 4e-5 off, Day's matrix at h = 1e-10 lost its imaginary part, 216 of 40,000 symmetric 2×2 pencils with close eigenvalues came out complex in double (234 under x86-64-v3), and random orthogonal 8×8 pencils had a worst backward error of 4.9e-13 against LAPACK's 3.0e-16. Now the rotations' imaginary parts are right to a few ε of themselves, Day's matrix gets ±5e-11, no symmetric pencil comes out complex in float or double, the orthogonal pencils reach 4.0e-16, and over 200,000 random 2×2 blocks in each of six families (tiny bᵢᵢ, a rotation times a triangular B, close diagonals, graded) the worst backward error is 1.4ε in float and 1.5ε in double, from up to 990ε and 2.3·10⁷ε. The real/complex decision compares |p| with √(−q) and forms neither p² nor q, so a small block does not underflow in float.
- `generalized_eigenvalues` and `transmission_zeros` scale A and B by powers of two (exact) so their largest entries lie in [1, 2), and multiply each eigenvalue back by 2^(e_A − e_B). The double-shift sweep starts from a first column quadratic in A·B⁻¹, which overflowed or underflowed once ‖A‖/‖B‖ left [2⁻⁶⁰, 2⁶⁰] (about 1e±18) in float and [2⁻⁴⁹⁷, 2⁵⁰⁸] in double, and the pencil then failed; near the bottom of the range (A and B both 2⁻¹¹⁷ in float) it failed at a ratio of one. Over a grid of (2ᵃA, 2ᵇB) for a 5×5 pencil with normal entries and eigenvalues, 407 of 901 scalings failed in float and 32,133 of 63,365 in double; none fails now, and every one has the digits of the unscaled pencil. A finite eigenvalue too large for the type now fails the result instead of coming back +Inf as a success.
- `svd` (and `rank`, `pseudo_inverse`, `null_space`, `two_norm` with it) prescales A by an exact power of two, entry by entry, instead of multiplying by 1/max|a|. Past 1/min (8.5e37 in float, 4.5e307 in double) that reciprocal is subnormal, and -ffast-math flushed it to zero: A became 0 and every σ came back 0 as a success, so float diag(1e38, 1e37) and double diag(1e308, 1e307) had rank 0. For a subnormal max|a| it overflowed, and diag(3e-310, 1e-310) got σ = 0 the same way without flush-to-zero (constant-evaluated, it did not compile). A matrix that needs scaling down now takes the least power of two that keeps its column norms finite, so float diag(3e38, 1, 0.1) keeps σ₂ = 1 and σ₃ = 0.1 to a few ε. A singular value above the largest value of the type (float [[3e38, 1e38], [2e38, 1e38]], σ₁ = 3.86e38) fails the SVD instead of returning +Inf, or 0 under flush-to-zero, as a success.
- `standardize_schur_block` exchanges a lower-triangular 2×2 block [[a, 0], [c, d]] by the exact permutation [[0, −1], [1, 0]], as LAPACK's xLANV2 does, and returns (d, a). The real-root path formed a as d + (a − d), which cancels when |a| ≪ |d|: [[−1e-3, 0], [1, −1e3]] gave −1e-3 with a relative error of 2.4e-11 in double and 2.3e-2 in float, [[−1e-8, 0], [1, −1e4]] 8e-5 in double. Both are now exact, in the QR iteration and in `reorder_schur`'s re-standardization.

## [0.2.0] — 2026-11-09

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
- Compilers: g++ 11.3+ and clang++ 19+. MSVC is not supported. Host CMake builds
  use `-mtune=native`. Every product sketch must fit clang's default constexpr
  budget (`tools/clang_budget_check.sh`).
- Every `*Result` struct is `[[nodiscard]]`, so a discarded design, analysis or
  simulation result warns at the call site.
- `CONFIG_CXX=clang++` builds on Windows: it targets MinGW's libstdc++, since Damp does not support MSVC.
- Examples and `make target-smoke` link statically, as the tests do, so they run whichever libstdc++ DLL is on PATH.
- `make` no longer needs `arm-none-eabi-g++`: the Cortex-M7 asm dumps, which also compile the README and guide snippets, skip themselves without it.
- `make embedded-check` looks for heap container includes in Damp's own headers, not in what a standard library reaches.

### Removed

- `CONFIG_BACKEND` and `CONFIG_COMPILER_PREFIX` in `tup.config`, which did not build. The ETL backend has its own contract checks, and PlatformIO does cross builds.
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

- GCC 14 no longer crashes (an internal compiler error in `build_data_member_initialization`) on a `constexpr` constructor holding a local whose member's default initializer copies a `Matrix`, as an `RlsConfig` does through its `Bounds`. The copy constructor reached `unset_elements`, which is not `constexpr`; it now builds its unset array itself, with the same code. GCC 11 crashed on the same code for the copy constructor's `is_constant_evaluated` branch; there it copies the array whole.
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
- `damp/math/math_backend.hpp` honours the profile when it is the first Damp header included.
- `assoc_barrier` works on GCC 11 for x86, where float and double `sin` lost accuracy under `-ffast-math`.
- `design::lowpass<3>` and `<4>` in float keep unit DC gain under `-ffast-math` whatever the compiler's summation order.
- `discretize` refuses a NaN `Ts` under `-ffast-math` on GCC 11.
- `dpwm_zero_sequence` is constexpr on Clang 19.
- `make embedded-check` reports a compile error instead of passing, and the Clang constexpr budget check finds MinGW's libstdc++ whatever is first on PATH.
- `make tidy`, `make iwyu` and `make fix_includes` drop GCC's `-fconstexpr-ops-limit`, which clang tools reject.
- CMake's `docs` target exists only at the top level, so it cannot clash with a parent project's.
- The README's FetchContent snippet skips the submodules: `libs/` is only for this repo's own tests and examples.

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

[Unreleased]: https://github.com/Wetmelon/damp/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/Wetmelon/damp/releases/tag/v0.2.0
[0.1.0]: https://github.com/Wetmelon/damp/releases/tag/v0.1.0
