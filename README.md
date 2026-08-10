# Damp

Compile-time controls for embedded C++. Header-only, C++20.

v0.1.0 (alpha) · [CHANGELOG](CHANGELOG.md) · BSL-1.0 ([LICENSE](LICENSE))

After enough years writing firmware for control systems, one will see the same issues:

1. Everyone reimplements the same algorithms on every project (often incorrectly).
2. Model-based control (LQR, ESKF, …) is desirable, but the integration effort involved in shipping real model-based controls is huge
3. Every real controls project turns into a variant management problem.  Maintaining one source of truth across many variants, engineers, and tools over several years is difficult.

Damp addresses these by keeping design and deploy in one tree. Heavy synthesis runs at compile time, the MCU keeps thin runtime objects. Plant models and controllers can live in the same headers with hardware variants as `constexpr` parameters rather than generated code hand-copied or awkwardly wired into the build system.

Example: Inverted-pendulum LQR in the shape that lands on target:

```cpp
#include "damp/control.hpp"
using namespace damp;

constexpr double g = 9.81, L = 1.0, m = 1.0, b_damp = 0.1, Ts = 0.01;

constexpr auto sys = StateSpace{
    .A = Matrix<2, 2>{{0.0, 1.0}, {g / L, -b_damp / (m * L * L)}},
    .B = Matrix<2, 1>{{0.0}, {1.0 / (m * L * L)}},
    .C = Matrix<1, 2>{{1.0, 0.0}},
};
constexpr auto Q = Matrix<2, 2>::identity() * 10.0;
constexpr auto R = Matrix<1, 1>{{1.0}};

constexpr auto result = design::lqrd(sys.A, sys.B, Q, R, Ts);
static_assert(result.success);   // build fails if DARE did not converge

// float gain on the MCU; the Riccati solver does not ship in the .elf
constinit LQR<2, 1> pendulum_controller{result.as<float>()};

// ISR: u = -K x  →  a couple of multiply-accumulates for this plant
// ColVec<1, float> u = pendulum_controller.control(x);
```

That is the usual design-is-deploy path (PID, LQR, ESKF, and most other laws look the
same). Worked examples: [`examples/control/cart_pole/`](examples/control/cart_pole/),
[`examples/control/pid/`](examples/control/pid/),
[`examples/estimation/eskf/`](examples/estimation/eskf/).

| Habit | What it does |
| ----- | ------------ |
| `constexpr` on design | Plants, gains, and `design::…` results evaluate at compile time (or init). Pair with `static_assert(result.success)` so a failed design is a build failure, not a silent wrong gain. |
| `.as<float>()` | Design stays in `double` for numerics; deploy is `float` for the MCU. Convert once at the boundary instead of mixing precisions on every tick. |
| `constinit` on the runtime object | Controllers and estimators in static storage with constant initialization: state lands in the binary image (`.hex` / load image), not via a dynamic constructor during `__start`. No heap and no first-use setup in the interrupt. |
| `StateSpace` / `TransferFunction` / `ZPK` | Shared LTI types for plants and interconnections (series `*`, feedback `/`, `discretize`, …). Prefer these over ad-hoc `A`/`B` arrays local to one module. |
| `damp::` math (and aliases), not `std::` on target | `damp::sin` / `sqrt` / `exp` / … go through `MathBackend<T>` (constexpr-friendly at design time, swappable on target). Same for `damp::array`, `optional`, `clamp`, and the other backed aliases — raw `std::` skips that path and freestanding builds. |

Also: keep `design::…` out of the PWM ISR; prefer `mat::solve` over forming inverses;
runtime types default to `float`, design to `double`.

## Usage

Firmware only needs `inc/` on the include path. The default backend is the normal C++
standard library — no git submodules required for `damp/control.hpp`.

```bash
g++ -std=c++20 -I path/to/inc your_firmware.cpp
```

```cpp
#include "damp/control.hpp"      // MCU path
// #include "damp/workbench.hpp" // host only (Bode, ODE, simulate, matlab::)
```

| Include | When |
| ------- | ---- |
| `damp/control.hpp` | MCU / embeddable core — heap-free umbrella |
| `damp/workbench.hpp` | Host analysis and simulation (do not pull into firmware) |

Small utilities (timers, scaling, bounds, …) live under `toolbox/` and ride along with
`control.hpp`. Submodules under `libs/` (fmt, json, plotlypp, etl) are for host tests,
examples, Plotly, and the optional freestanding ETL backend
(`-DDAMP_BACKEND_ETL` + `-I libs/etl/include`) — see
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). Initialize them only if you need those:

```bash
git clone --recurse-submodules https://github.com/Wetmelon/damp.git
# or after a plain clone:
git submodule update --init --recursive
```

Optional CMake: `find_package` / FetchContent → link `wetmelon::damp`.

GCC 10+, Clang 12+, MSVC 2022+.

## Examples

Copy an example folder rather than starting from a blank `main`. This release ships
classical control and estimation demos.

| Want | Example | Headers |
| ---- | ------- | ------- |
| PID / LQR / cart-pole | [`control/cart_pole/`](examples/control/cart_pole/), [`pendulum/`](examples/control/pendulum/), [`pid/`](examples/control/pid/) | `controllers/lqr.hpp` or `pid.hpp`, `systems/state_space.hpp` |
| Attitude ESKF / INS | [`estimation/eskf/`](examples/estimation/eskf/), [`ins_navigator/`](examples/estimation/ins_navigator/) | `estimation/eskf.hpp`, `estimation/sensor_fusion.hpp` |

```text
examples/<domain>/<name>/
  <name>_controller.hpp   # (or _estimator / _filter / _deploy) — flashes
  <name>_sketch.cpp       # thin board loop
  <name>_sil.cpp          # host plant + plots; calls the same tick
  <name>_derivation.md
```

SIL should call the deploy tick — not reimplement the controller for the sim. Details:
[examples/README.md](examples/README.md).

## What’s in the tree

Embeddable core: stack matrices, `constexpr` math, LTI types (SS / TF / ZPK) and
discretize, design-is-deploy controllers and estimators (PID, LQR family, ADRC, SMC,
KF/EKF/UKF, ESKF, …), filters, fixed-step integrators, geometry, toolbox helpers.

Host (via `workbench.hpp`): Bode/margins, ODE solvers, closed-loop simulate, and a
set of MATLAB®-style short names in `matlab.hpp` (not full toolbox parity).

Full name dump: [REFERENCE.md](REFERENCE.md).

```text
inc/damp/
  control.hpp       # embeddable core
  workbench.hpp     # host extras
  math/  matrix/  systems/  controllers/  design/
  estimation/  filters/  toolbox/  simulation/  analysis/
tests/  examples/  cmake/  targets/  docs/
```

## Building this repo

Only if you are developing Damp itself — not required just to use the headers.

```bash
make            # format, build, run tests
make tests
make examples
```

[tup](https://gittup.org/tup/) under the Makefile. Copy `tup.config.default` →
`tup.config` if you need a local compiler path.

## Trademarks

MATLAB®, Simulink®, and related MathWorks product names are trademarks of The
MathWorks, Inc. Arduino® is a trademark of Arduino SA. PlatformIO® is a trademark
of PlatformIO Labs UG. Other names may be trademarks of their owners. Use here
does not imply endorsement. The C++ names `matlab::` / `matlab.hpp` are library
API identifiers, not MathWorks products.
