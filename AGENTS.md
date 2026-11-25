# AGENTS.md

Guidance for humans and LLMs working in this repository. Prefer **short rules +
links** over re-stating the full library map here.

## Project

Header-only C++20 control library (`damp`): stack-only fixed sizes, no core heap,
`constexpr`-friendly design, firmware-ready runtime ticks.

**Product promise:** compile-time / init-time design in the same tree as firmware
(variants as `constexpr` / template params), closed in SIL with the same runtime
objects.

The public tree is the design-is-deploy **core**. Motor FOC is opt-in under
`damp/motor/` (not in `control.hpp`); there are no power or motion packs in this
tree.

## Where to look

| Topic | Doc |
| ----- | --- |
| Public API table | [REFERENCE.md](REFERENCE.md) (regen: `make refs` or `python3 tools/gen_reference.py`) |
| Product front door | [README.md](README.md), [docs/guide.md](docs/guide.md) (task table + MATLAB® / firmware names), [CHANGELOG.md](CHANGELOG.md) |

## Build and test

[Tup](https://gittup.org/tup/) + Makefile. Compiler path: `tup.config` (see
`tup.config.default`).

```bash
make              # format + compiledb + build + run tests
make tests        # build tests, then run
make clang-tests  # build the three test runners with clang too, then run (about 70 s)
tup --quiet tests # compile tests only
```

**Final check: always run `make`.** Before reporting work as done, run plain
`make` and confirm it exits 0. It is the only command that covers everything:
formatting, the test runner, both backend runners, examples (built *and* run)
and the embedded/freestanding/workbench/no-C++-headers contract checks.
It doesn't regenerate `REFERENCE*.md`: a PR leaves them alone, so stacked PRs
don't conflict on them, and `make refs` regenerates them on `main` after merges.
`tup tests`, `test_runner`, or a hand-built subset are for iterating, not a
substitute. `make` builds with g++ (and checks that the product sketches fit
clang's constexpr budget); a change to deduction, templates or NaN handling
also wants `make clang-tests`, which builds the runners with clang
(`clang-tests/`, skipped without clang++ 19+ on PATH).

**Windows: run `make` from PowerShell or cmd, not Git Bash.** Under Git Bash,
Git's MSYS `sh` comes first on PATH. MSYS tools write `*-pipe-nt-*` files that
tup rejects ("was written to, but is not in .tup/db"). The coreutils also go
missing inside tup jobs (`basename: command not found`). From PowerShell, `sh`
resolves to `C:\Tools\binutils\bin\sh.exe` (busybox), which works.

**Worktrees go outside the repo** (e.g. `C:/workspace/damp-<name>`). tup treats
any nested directory containing `tup.config` as a variant and refuses to build
the main tree.

**Test runner** (path is OS-specific; tup writes under `tests/build/`):

```bash
# Unix / Git Bash
./tests/build/test_runner
# Windows (tup often emits .exe)
./tests/build/test_runner.exe
```

Filter suites/cases:

```bash
./tests/build/test_runner.exe -ts="Matrix Functions"
./tests/build/test_runner.exe -tc="Matrix exponential - diagonal"
```

Refresh `compile_commands.json` (no full rebuild):

```bash
tup --quiet compiledb
```

**Python tools** — Makefile picks `python3` on Unix and `py -3` on Windows
(`PYTHON` override: `make refs PYTHON=python3`). Prefer shebang
`#!/usr/bin/env python3` in scripts. On Windows without `python3` on PATH:
`py -3 tools/...` or `make refs`.

**Toolchains:** g++ 11.3+ (tested with 14) and clang++ 19+ (C++20); MSVC is not
supported. Flags live in `Tuprules.lua`. Format runs via `make` / `clang-format`.

**Golden data** for the oracle tests: `tools/gen_oracle.py` writes
`tests/oracle_cases.hpp` (SciPy; run by hand, `make` never calls Python for it).

**Tup paths** are relative to the Tupfile’s directory (e.g. `tests/build/`), not
the repo root.

## Architecture (essentials)

| Include | Role |
| ------- | ---- |
| `damp/control.hpp` | Embeddable core |
| `damp/workbench.hpp` | Host: analysis, sim, `matlab::` |

- `damp::design::` — heavy CST synthesis only (not PWM-rate).
- **Types:** `Matrix` / `ColVec` / `RowVec` / views; `StateSpace`; Result structs
  with `.as<U>()`, `success`, and `explicit operator bool()`; geometry `DCM` /
  `Quaternion` / `Euler` (`EulerZYX` = aerospace YPR, `EulerXYZ` = robotics RPY).
- **LTI shape:** design continuous SS → `discretize(..., ZOH)` → step matrices.
  Prefer ZOH over hand-rolled Euler when an LTI model exists. An input delay is
  `discretize(sys, Ts, InputDelay<N>{θ})`, θ ≤ N·Ts in N delay states of the
  plant (`InputDelay{θ}` is N = 1), and a late reading `OutputDelay<N>{θ}` the
  same way; a `ColVec` of seconds gives one delay per channel. A delay is a model
  of the plant, never a parameter of a law's tick. `state_feedback` (every LQ
  law), `kalman`, `steady_state_target` and `mpc` take the continuous plant and
  the delay as a trailing option,
  `fn(plant, Ts, config, …, InputDelay<N>{θ})`, and design on that model, the
  cost and the noise discretized exactly; a measured-state law holds the delay
  states (`update` takes the plant state), a filter carries them. `kalman`, and
  a `state_feedback` or `mpc` with a `KalmanNoise`, take `OutputDelay<M>{θ}` the
  same way: the filter carries the late readings' delay states, the gain and the
  QP read none of them. Many samples of dead time is a slower sample or
  `SmithPredictor`.
- **Design calls:** what every call needs is positional (the plant, `Ts` for a
  continuous one, the weights: `state_feedback`'s `LQRWeights`, or `Poles{…}` for the
  gain placed there), so braces size it from the plant; the rest is a
  trailing pack recognized by type (`design/options.hpp`: `KalmanNoise`, `ObserverPoles`,
  `InputDelay`, `OutputDelay`, `IntegralAction`, `TrackedOutputs`, `InputLimits`, `InputRateLimits`, `DemandLimit`,
  `OutputLimits`, `StateLimits`, `DareCost`, `InvariantSet`, `QPIterations`, a disturbance spec, `MeasuredDisturbance`, `Handover`, …),
  at most one of each kind, in any order, refused at compile time otherwise. What the gain does not depend on is an
  option, not a weight: the limits on u and on Δu are `InputLimits{.lower, .upper}` and
  `InputRateLimits{.lower, .upper}` in both families (an LQ law clamps its command and its move, an MPC's QP carries
  rows for each, the Δu rows only with an `InputRateLimits`), and an LQ design bounds what its outer terms ask of a
  state with a `DemandLimit{.state, .outer, .lower, .upper}`; an MPC's output and state limits, its terminal ingredient and its QP budget are
  `OutputLimits{.outputs, .lower, .upper}`, `StateLimits{.states, .lower, .upper}`, `DareCost` or `InvariantSet` and
  `QPIterations`, and an operating
  point is the runtime's (`set_input_target`; Ru is measured from the target's input otherwise).
  Both families hold r through one target, `steady_state_target`'s: on the outputs a disturbance
  spec's `.outputs` lists, or without one a `TrackedOutputs`, every output when the one given
  lists none, at most one per input; with neither the law is a regulator and takes no r, an MPC
  as an LQR. An MPC's QP weighs and bounds every plant output, the others toward the target's values. A
  delay is an option of every sampled design (`InputDelay<N>{θ}` with the
  continuous plant and `Ts`), the runtime holding or filtering its delay states;
  a model with the delay built in is the plant route, for a model one holds. The
  options pick the law: `state_feedback`'s runtime is `StateFeedback`, and five
  names name what the options produce, each an alias of it that refuses another
  kind's design: LQR, LQI with an `IntegralAction`, LQG with a `KalmanNoise` (or
  `ObserverPoles`, the filter placed), LQGI with both, OFLQG with a filter and a
  disturbance spec, its weights `LQRWeights` or `Poles` for any of them, the same result and runtime; `mpc` runs as one
  runtime, `MPC`, on the measured state or through a filter with a `KalmanNoise`, and `OFMPC` names the
  offset-free one, with a noise and a disturbance spec (an alias: `OFMPC law{res}` refuses a design without a
  disturbance); a `MeasuredDisturbance` adds to any of them, the target under d and a filter's known
  input after u. An option carries its own numbers, so the
  weights are the same with or without it: `IntegralAction{.Q, .outputs, .limits,
  .tracking_time}` the integrators' cost and limits, `KalmanNoise` the noise and
  the disturbance's walk, `InputDelay<N>{θ}` and `OutputDelay<M>{θ}` the delays. A weight an option's
  size comes from is a value (a scalar or a Matrix); `KalmanNoise` deduces its
  sizes from its fields. A limit option's ends are flat, `.lower` and `.upper`, each unbounded when left out, as a
  `Bounds`' are (vectors written flat, `{-0.5}` for one channel; a `DemandLimit`'s scalars); a box inside an option
  that is not a limit is `.limits` (`IntegralAction`). A list of indices is a C array (`DemandLimit::outer`,
  `TrackedOutputs::outputs`, `OutputLimits::outputs`, `StateLimits::states`, the disturbance specs' `outputs`), so its length
  is deduced from bare braces, which no class member's is; the type with none
  listed (NZ = 0) is a specialization without the array. Two overloads per design, `fn(sys, …, opts...)` and
  `fn(plant, Ts, …, opts...)`. The pack is the one spelling: no options or
  config struct. A named option in a pack is left out, combined, forwarded and
  deduced inside a template; a designated-initializer aggregate's fields would
  carry the options' sizes, and clang does not deduce one inside a template
  (llvm #114870). A new option is a new type in the pack, which names its kind,
  `using option_kind = design::kind::…;` (a new kind is a tag in `design/options.hpp`).
- **Ticks:** estimators `predict` / `update`; controllers `control`; `step` is
  predict then update then control (m_k in, u_k out). Low-latency ISR: update →
  control → write → `back_calculate` → predict. `back_calculate(u_sat, u_ff = {})`
  records the applied input. While another source drives the plant, every LQ law
  and MPC takes `track(u, …)` in place of `step(…)`, a twin of each overload with
  the applied input first (void, no mode flag): predict, update, then the memory
  that makes control here return u (u recorded as applied, an LQI's ξ preloaded
  through `unwind`, an MPC's u_prev = u, no QP solved). A law with nothing to
  preload takes `design::Handover{.fade_time}`, which fades the first command's
  mismatch; `track` and a law without it cost `step` nothing.
- **Limits:** a law's output limit is a `Bounds`, unbounded by default: an
  `InputLimits` among an LQ or MPC design's options, `u_lim` in the design input
  or on the result of the others, refused when it is not a box. The LQ laws and
  the MPCs read their run-time limits with `limits()` and set them with
  `set_limits(lim)`, which refuses one that is not a box (false, nothing
  changed): `StateFeedbackLimits` and `MPCLimits`, which start with `input` (then
  `input_rate`) and end with `disturbance`, and a design result holds the same
  struct as `limits`. The other runtimes hold `u_lim` as a public member. The
  tick returns a command inside its limits, and with nothing else set the law's
  memory holds against that clamp (conditional integration, an observer or
  internal model on the clamped command). A clamp after the law is reported with
  `back_calculate`, to the same effect. A move limit is an `InputRateLimits`
  (clamped on an LQ law, rows of an MPC's QP); state and output limits are
  `design::mpc` (`OutputLimits`, `StateLimits`), `design::mhe` (`x_lim`) or an
  `ActionGovernor` around any law. A new law with memory gets all of it.
- **Errors:** `optional` / `bool success` / `DesignStatus` / `valid()` /
  compile-time constraints. No exceptions, and no run-time faults anywhere: a
  refusal is a compile error only under constant evaluation; at run time it is
  reported through the result or `valid()`, and the application decides.

**`matrix.hpp` includes** (order is irrelevant): `block`, `colvec`, `concat`,
`const_matrix`, `core`, `decomposition`, `eigen`, `functions`, `generalized_eigen`,
`pattern_matrix`, `rowvec`, `solve`, `svd`, `views`.

## Contributor rules

### Design-is-deploy three tiers (laws only)

```text
design::fn(...)  →  Result { gains, success, .as<U>() }  →  Runtime.tick(...)
```

Plant builders and host sim are not this ladder.

### Naming

- Descriptive English primary names (`design::state_feedback`); MATLAB-style aliases
  (`matlab::dlqr`) call through.
- Results: `FooResult`. Runtime ticks: short `PascalCase` (`LQR`, `PIDController`).

### Code style

- `constexpr`-first for design; runtime = light mat-vec / scalar.
- **Solve, don’t invert** — `mat::solve` / `cholesky_solve` / `lu_solve`.
- Braces always on `if`/`for`/`while`. ~One-screen functions (~50 lines).
- Templated literals: `T{1}` / `T{0.5}`; non-dyadic decimals
  `static_cast<T>(0.98)`.
- A literal names no size: a `ColVec` is written flat, `ColVec{a, b}` and `ColVec{x}`, a `Matrix` by rows, `Matrix{{a}, {b}}`, and both deduce their size and scalar from the entries (integers give double). A literal gives every entry, its shape checked at compile time (nothing is zero-filled; a few entries are set on `Matrix<R, C>{}`). The size stays for entries of mixed types, a scalar other than the entries' (`Matrix<2, 2, float>{{1.0, 0.0}, {0.0, 1.0}}`), and where CTAD does not reach: a non-static member, a parameter, a return type, a template argument. Braces around each entry of a `ColVec` build the same vector and draw clang's `-Wbraced-scalar-init`.
- Tolerances: scale `eps<T>()` to the data (size and magnitude) so results don't
  depend on units; `default_tol<T>()` only for convergence of an
  already-relative quantity. Always `damp::` math and backed aliases in
  `inc/damp` and tests (tests may use `std::` math as an oracle only). That
  includes traits, limits and concepts (`damp::is_same_v`, `damp::numeric_limits`,
  `damp::floating_point`): `make nostd-check` builds the core with no C++ standard
  headers. New names go in both branches of `backend.hpp`. No structured bindings
  on `damp::tuple`; use `damp::get`.
- No virtuals, heap, exceptions, or global mutable state in core.
- No backwards-compat burden: clarity over API stability.
- **`detail` namespaces:** an internal lives in the `detail` of the namespace whose code calls it:
  `damp::detail` for the runtimes and systems, `design::detail` for the designs, `mat::detail`,
  `sim::detail` and so on. Inside `damp::design`, `detail::` means `design::detail`, so
  `::damp::detail::` spelled there says a name is on the wrong side. A helper both sides call
  stays in `damp::detail`, and `design::detail` names it with a using-declaration, as
  `status.hpp` does `refused` and `set_status`. No one-off `foo_detail` namespaces.
- **Headers:** a name users spell is in a public header; `detail/` holds `detail::`
  only. A public header exports (`// IWYU pragma: export`) what its users name to
  call and deploy it: `Law` (`controllers/law.hpp`, through the header that
  specializes it), its configs and results, and the option types only its designs
  take (`KalmanNoise`, the disturbance specs). Vocabulary types (`Matrix`,
  `StateSpace`, `Bounds`, `InputDelay`) are the user's own includes, or come with
  `control.hpp`, the front door, which exports every public header,
  `backend.hpp` and `math/math.hpp`. A part of an umbrella (`matrix/core.hpp`, …) says
  `// IWYU pragma: private, include "damp/matrix/matrix.hpp"`, so the cleaner asks
  for the umbrella. Every `detail/` header says `// IWYU pragma: private`, and a
  private header's exports stop at it: a public header exports the public ones
  its users name itself (`lqr.hpp` exports `law.hpp`, not only through
  `detail/lq_law.hpp`). The examples include `control.hpp` (`workbench.hpp` in a
  SIL that simulates), and only the headers outside it (`motor/`, `kinematics/`,
  `plot_plotly.hpp`) by name. A `detail/` include that `make iwyu` asks for in an
  example is a missing export.

### Docs (Doxygen)

- Public surfaces need `@brief` (feeds `REFERENCE.md`). `@file` once per header.
- **Math placement**
  - `@brief` / `@param` / `@return` / `@retval` / `@tparam` / `///<` one-liners:
    **Unicode** (α, ω, ≤, Kᵢ, L_g h, [0,1], …) — no `@f$...@f$` on those tags
    (keeps REFERENCE and IDE hovers readable).
  - Detailed description body: LaTeX fine — `@f$...@f$` / `@f[...@f]`.
- **Citations:** primary source on the algorithm with
  `@see "Title" (Author, Year), §…`. Doxygen owns the bibliography — no separate
  reference list in repo prose. Match `@see` style on nearby headers of the same
  family when adding algorithms.
- MATLAB note form: `@note Compare with MATLAB®'s ...`.
- House style: `make doxygen-style-check`. Prefer `///<` for members.
- **Algorithm detail lives in the Doxygen** of the function that implements it: the
  method, its source (`@see`), what it guards against, its limits. Not in the
  CHANGELOG.

### CHANGELOG and release notes

- CHANGELOG: one line per entry, one or two sentences: what changed and what it
  means for a user, at most one number. Mechanism, measurements and citations go
  in the Doxygen and the commit message.
- GitHub release notes are shorter still: a summary line, one line per
  user-visible change, and a link to the CHANGELOG.
- Credit the published method, never "as LAPACK does": that reads as copied code.

### Examples (design-is-deploy)

Example folders: deploy header + sketch + host SIL + derivation. Gold:
`examples/control/cart_pole/`. SIL must **call** the flashable tick — never
re-derive the cascade “for the sim.” A SIL checks what it prints: it returns
nonzero when a result leaves its band (`FAIL:` and the claim), so `make` fails
on a change in behavior, not only on a crash.

### Numerical / `-ffast-math`

Host tests and PIO target smoke use `-ffast-math`. Design `constexpr` stays
abstract IEEE; runtime hardens around fast-math (`damp::isfinite` is bit-pattern
based; no fragile cancellations; no ±Inf sentinels).

A guard that refuses a floating-point value tests it with `damp::isfinite`,
`damp::finite_positive` or `damp::finite_non_negative`, never with a comparison
alone (`!(x > 0)`, `x <= 0`, `!(lo <= hi)`, `x != x`, `std::isnan`): under
`-ffinite-math-only` a comparison on a NaN may come out either way (g++ `-O3`
folds `!(x > 0)` to `x <= 0`, and x86 reads `NaN == 0` as true). A value that
must be finite refuses +Inf too; unbounded is `numeric_limits::max()`. A test of
a guard builds its NaN from bits at run time (through a volatile), as
`tests/test_nan_guards.cpp` does.

The same code must also hold where multiply-adds fuse (Cortex-M4F/M7, x86-64-v3) and under GCC 15, which folds more of fast-math. Two rules follow. A decision never waits for a cancellation to round to exactly 0: test the operands, or a tolerance scaled to them. A test never expects two inlined copies of the same formula to agree bit for bit: compare to a few ulps of the operands. Check a change with `-march=x86-64-v3` on the host.

### Safety profile

Inspired by MISRA C++, AUTOSAR C++14, JSF AV C++ and the C++ Core Guidelines;
not a compliance target. Core = everything reachable from `damp/control.hpp`.
Host headers (`analysis/`, `simulation/`, `matlab.hpp`, `workbench.hpp`) are
exempt.

**Enforced by the build** — keep it that way: `-fno-exceptions -fno-rtti` in
`checks/`; the `Tuprules.lua` warning set (conversions, double promotion, casts,
fallthrough, shadowing); clang-tidy.

**Runtime ticks** (`control()`, `update()`, `step()`, …)

- No allocation, no exceptions, no recursion, no unbounded loops.
- Bounded worst-case time: every loop bound is a template size or a fixed cap.
- No undefined behavior: no signed overflow, no out-of-range float→int casts, no
  out-of-bounds indexing. Guard, or use unsigned/modular arithmetic.

**Design functions**

- Iterative solvers carry an iteration cap and report non-convergence through
  `success` / `optional`. Never return an unconverged value as a success.
- Certify outputs where a cheap check exists (residual, closed-loop stability).
- Tolerances are relative to the data (scaled `eps<T>()`); no absolute literals.

**Types and initialization**

- Every member has a default member initializer.
- `detail::unfilled` skips `Matrix`'s run-time zero-fill. Use it only for a result whose
  every entry is written before it is read; constant evaluation still starts from zeros.
- Single-argument constructors are `explicit`, except documented
  assign-from-design constructors.
- A runtime constructed from a design result refuses a failed one when
  constant-evaluated: `refuse_failed_design(result.status)` names the cause.
  `static constinit LQR law{res};` then needs no separate
  `static_assert(res.success)`, and `Law<kRes, T>` of a `constexpr` design
  refuses one wherever it is declared.
- No virtuals, no RTTI, no type erasure on the product path.

**Annotations — only where they carry meaning**

- `constexpr` wherever possible. `[[nodiscard]]` on value-returning functions and
  on Result structs (`struct [[nodiscard]] FooResult`).
- `noexcept` on user-declared move operations and `swap`. Elsewhere optional:
  `-fno-exceptions` is the guarantee.
- No `final` / `override` (there are no virtuals). No trap paths: nothing
  reachable at run time is `[[noreturn]]`.

### Scope

Controls / DSP lane only. Generic embedded containers →
[ETL](https://www.etlcpp.com) (`libs/etl`), not reimplemented here.

**Backends:** `damp_profile.hpp` maps `array`/`optional`/… to `std::` or `etl::`.
The one ETL-shaped type Damp owns is `damp::optional` under `ETL_NO_STL`, because
`etl::optional` of anything but a fundamental type is not constexpr there.
Host workbench may use vector/string; plotlypp for plots.

### What not to do

- Runtime polymorphism / type erasure on the product path
- Heap in embeddable core
- Exceptions for control flow
- Config macros instead of templates/args
- Inventing ETL-shaped utilities

## Third-party (host)

| Lib | Role |
| --- | ---- |
| fmt | tests / examples |
| nlohmann/json | tests / examples |
| plotlypp | plotting |
| etl | freestanding / embedded companion |
| expected.hpp | available; not core default yet |
| doctest | tests |
