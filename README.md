# Damp

Compile-time controls for embedded C++. Header-only, C++20.

v0.1.0 (alpha) · [CHANGELOG](CHANGELOG.md) · [Guide](docs/guide.md) · BSL-1.0 ([LICENSE](LICENSE))

`design::lqr` (and the rest) run at compile time or init. The MCU keeps a small
runtime and calls `step` every tick. SIL uses that same object. Hardware variants
are `constexpr` parameters in the same headers.

## Add it

Point the compiler at `inc/`. Nothing to link for the core.

```bash
g++ -std=c++20 -I path/to/inc your_firmware.cpp
```

```cpp
#include "damp/control.hpp"      // MCU path
// #include "damp/workbench.hpp" // host only (Bode, ODE, simulate, matlab::)
```

CMake, FetchContent or `find_package`, same target:

```cmake
include(FetchContent)
FetchContent_Declare(damp
    GIT_REPOSITORY https://github.com/Wetmelon/damp.git
    GIT_TAG v0.1.0   # 0.x: pin a commit or tag; the API can move
)
FetchContent_MakeAvailable(damp)
target_link_libraries(your_fw PRIVATE wetmelon::damp)
```

GCC 11+ (tested with 14) and Clang 19+. MSVC is not supported. [Known issues](docs/known_issues.md) lists the compiler bugs that change results, among them GCC 15's `pow(0, y)` under `-ffast-math`.

## Using ETL instead of the standard library

Damp includes ETL as `<etl/...>`, so ETL can live anywhere: the `libs/etl` submodule is only what this repo's own tests use. Put ETL's `include/` folder on the include path, and add a `damp_profile.hpp` to a folder on the include path too:

```cpp
#pragma once
#define DAMP_BACKEND_ETL
#define DAMP_MATH_BACKEND_FREESTANDING // constexpr series math at run time; no libm
```

```bash
g++ -std=c++20 -I. -Ipath/to/damp/inc -isystem path/to/etl/include main.cpp
```

With CMake, `add_subdirectory` both trees and `target_link_libraries(your_fw PRIVATE wetmelon::damp etl::etl)`.

If the toolchain has no C++ standard headers, add an `etl_profile.h` to the same folder:

```cpp
#pragma once
#define ETL_NO_STL
#define ETL_NO_CPP_NAN_SUPPORT // quiet_NaN() from the NAN macro, never libm's nanf()
```

ETL reads this file itself, so these apply even if a file includes `<etl/...>` before Damp. In `damp_profile.hpp` they can arrive too late, and ETL then silently uses the standard library. On Clang 19 and 20, include Damp before ETL: with no C++ headers, `etl/tuple.h` needs two declarations that Damp provides.

Use the profile rather than `-DDAMP_BACKEND_ETL`: every file in the program must agree on the backend, and Damp warns when it sees a backend macro but no profile. Without a C++ header on the toolchain, `DAMP_MATH_BACKEND_DAMP` and the default math backend are out, because they need `<cmath>`: use `DAMP_MATH_BACKEND_FREESTANDING` or your own `DAMP_MATH_BACKEND_HEADER`. The host-only `damp/workbench.hpp` still needs the standard library. [Known limitations](docs/known_limitations.md) has what is tested.

## Plant → law → tick

```cpp
#include "damp/control.hpp"
using namespace damp;

constexpr double g = 9.81, L = 1.0, m = 1.0, b_damp = 0.1, Ts = 0.01;

constexpr auto sys = StateSpace{
    .A = Matrix<2, 2>{{0.0, 1.0}, {g / L, -b_damp / (m * L * L)}},
    .B = Matrix<2, 1>{{0.0}, {1.0 / (m * L * L)}},
    .C = Matrix<1, 2>{{1.0, 0.0}},
};
constexpr auto result = design::lqr(sys, Ts, {.Q = Matrix<2, 2>::identity() * 10.0, .R = 1.0});

constinit Law<result, float> pendulum; // the design's LQR, in float

// ISR: measurement x at this sample, command u
ColVec<2, float> x{0.1f, 0.0f};
ColVec<1, float> u = pendulum.step(x);
```

That is the design-is-deploy path. The deploy line chooses where the design runs:

| You write | The design runs | The gains live | A failed design | Tick cost |
| --------- | --------------- | -------------- | --------------- | --------- |
| `constinit Law<result, float> law;` | In the compiler | Constants of the type, in flash | Does not compile | Nonzero terms only |
| `constinit LQR law{result.as<float>()};` | In the compiler | Members, in RAM | Does not compile | Dense |
| `LQR law{design::lqr(plant, Ts, w).as<float>()};` | On the MCU, at init | Members, in RAM | The result's `success` is false, and the law keeps zero gains | Dense, plus the design at init |

`Law<result, float>` is the default when the plant and the weights are known at build time. The design cannot slip to run time unnoticed: a template argument has to be a constant expression, so a result computed at run time does not compile there. The middle row keeps a compile-time design replaceable by a redesign at run time; the last row is for a plant measured or identified at start-up. More in [docs/guide.md](docs/guide.md#where-the-design-runs).

PID, LQG, MPC, ESKF look the same: design result → runtime → `step`. Flagship (LQI + clamp + `back_calculate`): [`examples/control/lqi/`](examples/control/lqi/).

If you know MATLAB® (`c2d`, `lqr`, `pidtune`, `kalman`) or you are writing an ISR:
[docs/guide.md](docs/guide.md) maps those names onto Damp.

| Habit | What it does |
| ----- | ------------ |
| `constexpr` on design | Plants, gains, and `design::…` results evaluate at compile time (or init). A runtime built from a `constexpr` result refuses a failed one, so a failed design is a build failure, not a silent wrong gain. At init, check `result.success`. |
| `Law<result, float>` or `.as<float>()` | Design stays in `double` for numerics; deploy is `float` for the MCU. Convert once at the boundary instead of mixing precisions on every tick. |
| `constinit` on the runtime object | Controllers and estimators in static storage with constant initialization: state lands in the binary image (`.hex` / load image), not via a dynamic constructor during `__start`. No heap and no first-use setup in the interrupt. |
| `StateSpace` / `TransferFunction` / `ZPK` | Shared LTI types for plants and interconnections (series `*`, feedback `/`, `discretize`, …). Prefer these over ad-hoc `A`/`B` arrays local to one module. |
| `damp::` math (and aliases), not `std::` on target | `damp::sin` / `sqrt` / `exp` / … go through `MathBackend<T>` (constexpr-friendly at design time, swappable on target). Same for `damp::array`, `optional`, `clamp`, and the other backed aliases — raw `std::` skips that path and freestanding builds. |

Also: keep `design::…` out of the PWM ISR; prefer `mat::solve` over forming inverses;
runtime types default to `float`, design to `double`.

## Usage

| Include | When |
| ------- | ---- |
| `damp/control.hpp` | MCU / embeddable core — heap-free umbrella |
| `damp/workbench.hpp` | Host analysis and simulation (do not pull into firmware) |

The default backend is the normal C++ standard library — no git submodules for
`damp/control.hpp`. Small utilities (timers, scaling, bounds, …) live under
`toolbox/` and ride along with `control.hpp`. Submodules under `libs/` (fmt, json,
plotlypp, etl) are for host tests, examples and Plotly — see
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). Initialize them only if you need those:

```bash
git clone --recurse-submodules https://github.com/Wetmelon/damp.git
# or after a plain clone:
git submodule update --init --recursive
```

## Examples

Copy an example folder. Task table: [docs/guide.md](docs/guide.md).

| Want | Example | Headers |
| ---- | ------- | ------- |
| LQI + clamp | [`control/lqi/`](examples/control/lqi/) | `controllers/lqi.hpp` |
| PID / LQR / cart-pole | [`control/cart_pole/`](examples/control/cart_pole/), [`pendulum/`](examples/control/pendulum/), [`pid/`](examples/control/pid/) | `controllers/lqr.hpp` or `pid.hpp` |
| Attitude ESKF / INS | [`estimation/eskf/`](examples/estimation/eskf/), [`ins_navigator/`](examples/estimation/ins_navigator/) | `estimation/eskf.hpp` |

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

Embeddable core: stack matrices, `constexpr` math, LTI types (SS / TF / ZPK) and discretize (with input delay), design-is-deploy controllers and estimators (PID, LQR family, OFLQG, MPC, ADRC, ESO, SMC, fixed-gain and reduced-order observers, KF/EKF/UKF, ESKF, …), filters, encoders, fixed-step integrators, geometry, toolbox helpers.

Host (via `workbench.hpp`): Bode/margins, ODE solvers, closed-loop simulate, and a
set of MATLAB®-style short names in `matlab.hpp` (not full toolbox parity).

Full name dump: [REFERENCE.md](REFERENCE.md).

```text
inc/damp/
  control.hpp       # embeddable core
  workbench.hpp     # host extras
  math/  matrix/  systems/  controllers/  design/
  estimation/  filters/  encoder/  toolbox/  simulation/  analysis/
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
