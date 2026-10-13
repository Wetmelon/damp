# AGENTS.md

Guidance for humans and LLMs working in this repository. Prefer **short rules +
links** over re-stating the full library map here.

## Project

Header-only C++20 control library (`damp`): stack-only fixed sizes, no core heap,
`constexpr`-friendly design, firmware-ready runtime ticks.

**Product promise:** compile-time / init-time design in the same tree as firmware
(variants as `constexpr` / template params), closed in SIL with the same runtime
objects.

Public v0.1 is the design-is-deploy **core** (no vertical motor/power/motion packs
in this tree).

## Where to look

| Topic | Doc |
| ----- | --- |
| Public API table | [REFERENCE.md](REFERENCE.md) (regen: `make refs` or `python3 tools/gen_reference.py`) |
| Product front door | [README.md](README.md), [CHANGELOG.md](CHANGELOG.md) |

## Build and test

[Tup](https://gittup.org/tup/) + Makefile. Compiler path: `tup.config` (see
`tup.config.default`).

```bash
make              # format + compiledb + build + run tests
make tests        # build tests, then run
tup --quiet tests # compile tests only
```

**Final check: always run `make`.** Before reporting work as done, run plain
`make` and confirm it exits 0. It is the only command that covers everything:
formatting, the test runner, both backend runners, examples (built *and* run),
the embedded/freestanding contract checks, and the `REFERENCE.md` regen.
`tup tests`, `test_runner`, or a hand-built subset are for iterating, not a
substitute.

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

**Toolchains:** GCC 10+, Clang 12+, MSVC 2022+ (C++20). Flags live in
`Tuprules.lua`. Format runs via `make` / `clang-format`.

**Tup paths** are relative to the Tupfile’s directory (e.g. `tests/build/`), not
the repo root.

## Architecture (essentials)

| Include | Role |
| ------- | ---- |
| `damp/control.hpp` | Embeddable core |
| `damp/workbench.hpp` | Host: analysis, sim, `matlab::` |

- `damp::design::` — heavy CST synthesis only (not PWM-rate).
- **Types:** `Matrix` / `ColVec` / `RowVec` / views; `StateSpace`; Result structs
  with `.as<U>()` + `success`; geometry `DCM` / `Quaternion` / `Euler`
  (`EulerZYX` = aerospace YPR, `EulerXYZ` = robotics RPY).
- **LTI shape:** design continuous SS → `discretize(..., ZOH)` → step matrices.
  Prefer ZOH over hand-rolled Euler when an LTI model exists.
- **Errors:** `optional` / `bool success` / compile-time constraints. No
  exceptions on the tick path.

**`matrix.hpp` bottom-includes** (do not reorder): `block`, `cholesky`, `colvec`,
`matrix_functions`, `rowvec`, `views`.

## Contributor rules

### Design-is-deploy three tiers (laws only)

```text
design::fn(...)  →  Result { gains, success, .as<U>() }  →  Runtime.tick(...)
```

Plant builders and host sim are not this ladder.

### Naming

- Descriptive English primary names (`discrete_lqr`); MATLAB-style aliases
  (`dlqr`) call through.
- Results: `FooResult`. Runtime ticks: short `PascalCase` (`LQR`, `PIDController`).

### Code style

- `constexpr`-first for design; runtime = light mat-vec / scalar.
- **Solve, don’t invert** — `mat::solve` / `cholesky_solve` / `lu_solve`.
- Braces always on `if`/`for`/`while`. ~One-screen functions (~50 lines).
- Templated literals: `T{1}` / `T{0.5}`; non-dyadic decimals
  `static_cast<T>(0.98)`.
- Tolerances: scale `eps<T>()` to the data (size and magnitude) so results don't
  depend on units; `default_tol<T>()` only for convergence of an
  already-relative quantity. Always `damp::` math and backed aliases in
  `inc/damp` and tests (tests may use `std::` math as an oracle only).
- No virtuals, heap, exceptions, or global mutable state in core.
- No backwards-compat burden: clarity over API stability.

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

### Examples (design-is-deploy)

Example folders: deploy header + sketch + host SIL + derivation. Gold:
`examples/control/cart_pole/`. SIL must **call** the flashable tick — never
re-derive the cascade “for the sim.”

### Numerical / `-ffast-math`

Host tests and PIO target smoke use `-ffast-math`. Design `constexpr` stays
abstract IEEE; runtime hardens around fast-math (`damp::isfinite` is bit-pattern
based; no fragile cancellations; no ±Inf sentinels).

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
- Single-argument constructors are `explicit`, except documented
  assign-from-design constructors.
- A runtime constructed from a design result refuses a failed one when
  constant-evaluated: `if (!r.success && std::is_constant_evaluated())
  ::damp::detail::refuse_failed_design();`. `static constinit Law law{res};`
  then needs no separate `static_assert(res.success)`.
- No virtuals, no RTTI, no type erasure on the product path.

**Annotations — only where they carry meaning**

- `constexpr` wherever possible. `[[nodiscard]]` on value-returning functions and
  on Result structs (`struct [[nodiscard]] FooResult`).
- `noexcept` on user-declared move operations and `swap`. Elsewhere optional:
  `-fno-exceptions` is the guarantee.
- No `final` / `override` (there are no virtuals). `[[noreturn]]` on trap paths.

### Scope

Controls / DSP lane only. Generic embedded containers →
[ETL](https://www.etlcpp.com) (`libs/etl`), not reimplemented here.

**Backends:** `damp_profile.hpp` maps `array`/`optional`/… to `std::` or `etl::`.
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
