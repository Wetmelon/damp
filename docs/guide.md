# Using Damp

Short paths from a plant to a tick. Exhaustive names: [REFERENCE.md](../REFERENCE.md).
Rules: [AGENTS.md](../AGENTS.md). Limits: [known_limitations.md](known_limitations.md).

The plant is the process under control (motor, heater, pendulum). Damp holds a
model of it: command $`u`$ in, measurement $`y`$ (and often state $`x`$) out —
`StateSpace`, a transfer function, or `Dynamics`.

## Canonical path

[`examples/control/lqi/`](../examples/control/lqi/) is the one to copy first: transfer
function → `to_state_space` → `discretize` → `design::lqi` → `LQI` → clamp →
`back_calculate` → SIL on the same tick.

```cpp
#include "damp/control.hpp"
using namespace damp;

namespace {
constexpr double m = 1.0; // kg
constexpr double k = 1.0; // N/m
constexpr double c = 0.5; // N·s/m
// G(s) = 1 / (m s² + c s + k). Companion form is [position, speed] when m = 1.
constexpr auto G = TransferFunction<1, 3>{.num = {1.0}, .den = {m, c, k}};
constexpr auto sys = *G.to_state_space();
}

constexpr auto model = *discretize(sys, 0.01);
constexpr auto res = design::lqi(model, {.Q = Matrix<2, 2>{{10.0, 0.0}, {0.0, 1.0}}, .R = 1.0, .Q_integral = 80.0});
static_assert(res.success);
constinit LQI law{res.as<float>()};

static void write_force(float /*u*/) {} // PWM / DAC on the board

int main() {
    ColVec<2, float> x{}; // measured state from sensors
    const float      r = 1.0f;
    law.update(x);
    const float u = law.control(ColVec<1, float>{{r}})[0];
    const float u_sat = clamp(u, -1.2f, 1.2f);
    law.back_calculate(ColVec<1, float>{{u_sat}});
    law.predict();
    write_force(u_sat);
}
```

The SIL in that folder runs this tick next to a naive loop that reports the unsaturated
request. The integrator winds on the naive path and overshoots; `back_calculate(u_sat)`
holds it while the clamp is against the move.

## Two headers

| Include | When |
| ------- | ---- |
| `damp/control.hpp` | Firmware and SIL ticks. No heap. |
| `damp/workbench.hpp` | Host only: Bode, margins, `simulate`, `matlab::c2d` / `lqr` / `pidtune` / `bode` |

`matlab::` is MATLAB-inspired spelling, not a drop-in Control System Toolbox.
Argument order, defaults, and return types follow the Damp primary. Design and
conversion aliases (`lqr`, `c2d`, `pidtune`, `place`, …) return that primary's
Result (`success`, `status`, `.as<float>()`). `ss` / `tf` / `zpk` return system
types. Bode and friends return analysis grids. New code uses `design::` and
`discretize`.

## I searched for…

| I typed | Use |
| ------- | --- |
| `ss`, `tf`, `zpk` | `StateSpace`, `TransferFunction`, `ZPK`; `matlab::ss` / `tf` / `zpk` on the host |
| `c2d`, `d2c`, `d2d` | `discretize` / `c2d` / `d2d`; `matlab::c2d` |
| `InputDelay`, dead time | `discretize(sys, Ts, InputDelay{θ})` with θ ≤ Ts. Longer: slower Ts or `SmithPredictor`. Not a field on the law. |
| `lqr`, `dlqr`, `lqrd` | `design::lqr(sys, …)` or `design::lqr(sys, Ts, …)`; `matlab::lqr` / `dlqr` / `lqrd` |
| `lqi`, `lqg`, `lqgtrack` | `design::lqi`, `design::lqg`, `design::lqgi` |
| `place`, `acker` | `design::place`; `matlab::place` / `acker` |
| `pidtune`, Ziegler–Nichols | `design::pidtune`, `design::ziegler_nichols_*`; then `.discretize(Ts)` → `PIDController` |
| `kalman`, `lqe`, `dlqe`, `kalmd` | `design::kalman` → `Luenberger` (fixed gain on the law) or `KalmanFilter` (P each tick) |
| observer, Luenberger | `Luenberger{design::kalman(...)}` |
| anti-windup, sat, `umax` | `back_calculate(u_sat)` after the driver clamp; PID `u_lim` on `PIDGains` |
| `mpcmove`, constraints | `design::mpc` / `ofmpc`, `MPC.step(r, y)` |
| `bode`, `nyquist`, `margin`, `rlocus` | `workbench.hpp`; `analysis::` or `matlab::bode` |
| `lsim`, `step`, `forced_response` | `sim::simulate*` on the host |
| ISR, float, no malloc | `control.hpp`, `.as<float>()`, `constinit`, `step` in the interrupt |

## Coming from MATLAB®

Keep the plant as `StateSpace` (or `matlab::ss` / `tf`). Discretize with `discretize` or `matlab::c2d`. Design with `matlab::lqr` if you want that spelling; deploy with `LQR law{res.as<float>()}` and `law.step(x)` — there is no `lsim` of a gain matrix on the MCU.

`[K,S,e] = lqr(sys,Q,R)` is `design::lqr(sys, {.Q, .R})` (CARE if continuous, DARE if discrete). Sampled-data (`lqrd`) is `design::lqr(sys, Ts, {.Q, .R})`.

Integral action is `design::lqi`; it builds the augmented plant. Actuator saturation on that law is `back_calculate`. Input and output limits over a horizon are `design::mpc`.

## Coming from firmware

One object, one tick. Design at compile time or init. Do not run Riccati in the PWM interrupt.

```cpp
#include "damp/control.hpp"
constinit LQI law{kRes.as<float>()}; // kRes from design::lqi
// ISR:
auto u = law.step(r, x);           // or update / control / back_calculate / predict
write_pwm(u[0]);
```

PID: `design::pidtune` or `design::pid(Kp, Ki, Kd)` → `.discretize(Ts)` → `PIDController`. Limits go on the gains (`u_lim`), not as extra ISR code, unless a downstream stage clamps — then `back_calculate`.

## Coming from Python Control / SciPy

`control.ss`, `control.lqr`, `control.c2d` map as in the MATLAB table. Frequency-domain helpers allocate (`std::vector`) and live in `workbench.hpp`. There is no Python runtime; you generate C++ and flash the tick.

## I want to…

| I want to | Start with | Example |
| --------- | ---------- | ------- |
| Regulate a fully measured plant | `design::lqr` → `LQR`, `step(x)` | [`examples/control/pendulum/`](../examples/control/pendulum/) |
| Hold listed outputs at r, trusting the model | `TrackedOutputs` + `design::lqr` / `lqg` | [`examples/control/cart_pole/`](../examples/control/cart_pole/) |
| Reject a constant disturbance / model offset | `design::lqi` / `lqgi`, or `design::oflqg` | [`examples/control/lqi/`](../examples/control/lqi/), [`heater_lqg/`](../examples/control/heater_lqg/) |
| Estimate the state from y | `design::lqg` or `design::kalman` → `KalmanFilter` / `Luenberger` | [`examples/estimation/eskf/`](../examples/estimation/eskf/) |
| Enforce input or output limits | `design::mpc` / `ofmpc`, `back_calculate` after a clamp | [`examples/control/heater_ofmpc/`](../examples/control/heater_ofmpc/) |
| A SISO loop with a PID | `design::pidtune` → `PIDController` | [`examples/control/pid/`](../examples/control/pid/) |
| A nonlinear plant | `Dynamics` / `Measurement`, `design::linearize`, EKF/UKF | [`examples/control/cart_pole/`](../examples/control/cart_pole/) SIL |
| Frequency response, Bode, margins | `workbench.hpp` | [`examples/control/workflow_end_to_end/`](../examples/control/workflow_end_to_end/) |

Copy an example folder. SIL calls the same tick as the sketch.

## LQR with TrackedOutputs, or LQI

Both can hold an output at r.

- **TrackedOutputs + LQR / LQG** uses the model's steady-state maps $`x^* = N_x r`$, $`u^* = N_u r`$. Offset-free only if the model is right and there is no unmatched load.
- **LQI / LQGI** adds an integrator on $`r - C_i y`$. It rejects a constant disturbance and a constant model error. Use this when the plant will not match the design model.
- **OFLQG / OFMPC** estimate a disturbance and hold r through a steady-state target. Same job as LQI when a disturbance model fits; MPC when the actuator or output limits bind.

## Tick

`step(r, m)` is predict, then update(m), then control. Low-latency ISR: update → control → write → `back_calculate` → predict. After a clamp, pass the applied input to `back_calculate(u_sat, u_ff)`.

## Add the library

Headers only: `-I path/to/inc` and `#include "damp/control.hpp"`. CMake: link `wetmelon::damp` (FetchContent or `find_package`). Host analysis: `damp/workbench.hpp`.
