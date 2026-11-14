# Using Damp

Short paths from a plant to a tick. Exhaustive names: [REFERENCE.md](../REFERENCE.md).
Rules: [AGENTS.md](../AGENTS.md). Limits: [known_limitations.md](known_limitations.md).
Compiler bugs: [known_issues.md](known_issues.md).

The plant is the process under control (motor, heater, pendulum). Damp holds a
model of it: command $`u`$ in, measurement $`y`$ (and often state $`x`$) out —
`StateSpace`, a transfer function, or `Dynamics`.

## Canonical path

[`examples/control/lqi/`](../examples/control/lqi/) is the one to copy first:
`StateSpace` in the coordinates you measure → `discretize` → `design::lqi` →
`LQI` → clamp → `back_calculate` → SIL on the same tick.

```cpp
#include "damp/backend.hpp"
#include "damp/control.hpp"
using namespace damp;

// Fake function representing some PWM/DAC/comms
static void write_force(float /*u*/) {}

namespace {
constexpr double Ts = 0.01; // [sec] Controller sample time

constexpr double m = 1.0; // kg
constexpr double k = 1.0; // N/m
constexpr double c = 0.5; // N·s/m

// x = [position, speed], u = force, y = position. G(s) = 1 / (m s² + c s + k)
constexpr auto sys = StateSpace<2, 1, 1>{
    .A = Matrix<2, 2>{{0.0, 1.0}, {-k / m, -c / m}},
    .B = Matrix<2, 1>{{0.0}, {1.0 / m}},
    .C = Matrix<1, 2>{{1.0, 0.0}},
};
} // namespace

// Configure tuning parameters (control weights) for an LQI
constexpr auto weights = design::LQRWeights<2, 1, 1>{
    .Q = Matrix<2, 2>{
        {10.0, 0.0},
        {0.0, 1.0},
    },
    .R = 1.0,
    .Q_integral = 80.0,
};

// Discretize the plant and solve the DARE
constexpr auto res = design::lqi(sys, Ts, weights);

// The design's LQI in single precision, its gains constants of the type; a failed design does not compile
constinit Law<res, float> law;

int main() {
    ColVec<2, float> x{}; // position, speed from sensors

    // Update the state (measurements)
    law.update(x);

    const float r = 1.0f;
    const float u = law.control(r)[0];         // Tick the controller and get the first output
    const float u_sat = clamp(u, -1.2f, 1.2f); // Saturate the control output
    law.back_calculate(u_sat);                 // Make the controller aware of the saturation
    law.predict();                             // Run the integrator (internally clamped)

    // Apply the saturated control
    write_force(u_sat);
}
```

The SIL in that folder runs this tick next to a naive loop that reports the unsaturated
request. The integrator winds on the naive path and overshoots; `back_calculate(u_sat)`
holds it while the clamp is against the move.

## Where the design runs

The design functions are the same in every case; the line that deploys the law decides where the design runs and where its numbers live.

| You write | The design runs | The gains live | A failed design | Tick cost |
| --------- | --------------- | -------------- | --------------- | --------- |
| `constinit Law<kRes, float> law;` | In the compiler | Constants of the type, in flash | Does not compile | Nonzero terms only |
| `constinit LQR law{kRes.as<float>()};` | In the compiler | Members, in RAM | Does not compile | Dense |
| `LQR law{design::lqr(plant, Ts, w).as<float>()};` | On the MCU, at init | Members, in RAM | The result's `success` is false, and the law keeps zero gains | Dense, plus the design at init |

- **Fixed at compile time, `Law<kRes, float>`.** The default when the plant and the weights are known at build time. `kRes` is a named `constexpr` result with static storage (namespace scope, or a `static` member or local). The matrices are constants of the type: they take no RAM, an entry that is exactly 0 is no code and a 1 is no multiply. Nothing can replace them at run time. The constants sit in flash, so on a part with flash wait states and no cache, run the tick from RAM or tightly coupled memory.
- **Designed at compile time, gains in RAM.** For a law a redesign replaces at run time: assign another result of the same shape. A struct member of that kind is `Law<kRes, float>::Runtime`. The design is still checked when the program is built, and every tick pays for the dense matrices.
- **Designed on the MCU at init.** For a plant measured or identified at start-up. The design's iterations are capped and it reports `success` and a `DesignStatus`; check them before the first tick. Keep `design::…` out of the ISR.

## Two headers

| Include              | When                                                                             |
| -------------------- | -------------------------------------------------------------------------------- |
| `damp/control.hpp`   | Firmware and SIL ticks. No heap.                                                 |
| `damp/workbench.hpp` | Host only: Bode, margins, `simulate`, `matlab::c2d` / `lqr` / `pidtune` / `bode` |

`matlab::` is MATLAB-inspired spelling, not a drop-in Control System Toolbox.
Argument order, defaults, and return types follow the Damp primary. Design and
conversion aliases (`lqr`, `c2d`, `pidtune`, `place`, …) return that primary's
Result (`success`, `status`, `.as<float>()`). `ss` / `tf` / `zpk` return system
types. Bode and friends return analysis grids. New code uses `design::` and
`discretize`.

## I searched for…

| I typed                               | Use                                                                                                              |
| ------------------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| `ss`, `tf`, `zpk`                     | `StateSpace`, `TransferFunction`, `ZPK`; `matlab::ss` / `tf` / `zpk` on the host                                 |
| `c2d`, `d2c`, `d2d`                   | `discretize` / `c2d` / `d2d`; `matlab::c2d`                                                                      |
| `InputDelay`, dead time               | `discretize(sys, Ts, InputDelay<N>{θ})` with θ ≤ N·Ts (`InputDelay{θ}`: N = 1). One delay per input: `InputDelay<N, ColVec<NU>>{{θ₁, θ₂}}`. Many samples: slower Ts or `SmithPredictor`. Not a field on the law. Every sampled LQ design and `kalman` take it in the design, `design::lqgi(plant, Ts, InputDelay<N>{θ}, weights, noise)`: the plant's weights and noise, the cost exact on the delayed model, and a law that holds or carries the registers itself. MPC is designed on the model. |
| `OutputDelay`, a late sensor          | `discretize(sys, Ts, OutputDelay<N>{θ})`: the reading that arrives now was taken θ ago, held in N registers. One delay per output: `OutputDelay<N, ColVec<NY>>{{θ₁, θ₂}}`. Both: `discretize(sys, Ts, InputDelay{θu}, OutputDelay{θy})`. A filter on the model estimates the current state. |
| A tick that costs only its nonzero terms | `static constinit Law<kRes, float> law;` on a `constexpr` design: the same ticks, the matrices constants of the type, no RAM for them and no code for an exact 0 (MPC keeps its QP data in members, and OFMPC's filter alone is on constants). A member a redesign assigns is `Law<kRes, float>::Runtime`. On a `design::kalman`, `Law<kRes, float>::Recursive` is the recursive `KalmanFilter` on the same constants. |
| A Jacobian with a fixed shape, in a filter | `SparsityPattern` + `PatternMatrix<kPattern, float>`: `eskf.predict(F)`, `eskf.update_innovation(nu, H, R)`, `kf.update(y, C, R)`, or an EKF measurement whose `jacobian(x, u)` returns one. The zeros and ones are known at compile time, so they cost nothing; the values stay runtime. |
| `lqr`, `dlqr`, `lqrd`                 | `design::lqr(sys, …)` or `design::lqr(sys, Ts, …)`; `matlab::lqr` / `dlqr` / `lqrd`                              |
| `lqi`, `lqg`, `lqgtrack`              | `design::lqi`, `design::lqg`, `design::lqgi`                                                                     |
| `place`, `acker`                      | `design::place`; `matlab::place` / `acker`                                                                       |
| `pidtune`, Ziegler–Nichols            | `design::pidtune`, `design::ziegler_nichols_*`; then `.discretize(Ts)` → `PIDController`                         |
| `kalman`, `lqe`, `dlqe`, `kalmd`      | `design::kalman` → `FixedGainObserver` (fixed gain on the law) or `KalmanFilter` (P each tick)                          |
| observer, FixedGainObserver                  | `FixedGainObserver{design::kalman(...)}`                                                                                |
| ESO, ADRC                             | `design::eso` → `ESO` (observer only); `design::adrc` → `ADRCController` (ESO + PD)                              |
| anti-windup, sat, `umax`              | `back_calculate(u_sat)` after the driver clamp; PID `u_lim` on `PIDGains`                                        |
| `mpcmove`, constraints                | `design::mpc` / `ofmpc`, `MPC.step(r, y)`                                                                        |
| `bode`, `nyquist`, `margin`, `rlocus` | `workbench.hpp`; `analysis::` or `matlab::bode`                                                                  |
| `lsim`, `step`, `forced_response`     | `sim::simulate*` on the host                                                                                     |
| ISR, float, no malloc                 | `control.hpp`, `.as<float>()`, `constinit`, `step` in the interrupt                                              |

## Coming from MATLAB®

Keep the plant as `StateSpace` (or `matlab::ss` / `tf`). Discretize with `discretize` or `matlab::c2d`. Design with `matlab::lqr` if you want that spelling; deploy with `LQR law{res.as<float>()}` and `law.step(x)` — there is no `lsim` of a gain matrix on the MCU.

`[K,S,e] = lqr(sys,Q,R)` is `design::lqr(sys, {.Q, .R})` (CARE if continuous, DARE if discrete). Sampled-data (`lqrd`) is `design::lqr(sys, Ts, {.Q, .R})`.

Integral action is `design::lqi`; it builds the augmented plant. Actuator saturation on that law is `back_calculate`. Input and output limits over a horizon are `design::mpc`.

## Coming from firmware

One object, one tick. Design at compile time or init. Do not run Riccati in the PWM interrupt.

```cpp
#include "damp/control.hpp"
constinit Law<kRes, float> law; // kRes from design::lqi
// ISR:
auto u = law.step(r, x);           // or update / control / back_calculate / predict
write_pwm(u[0]);
```

PID: `design::pidtune` or `design::pid(Kp, Ki, Kd)` → `.discretize(Ts)` → `PIDController`. Limits go on the gains (`u_lim`), not as extra ISR code, unless a downstream stage clamps — then `back_calculate`.

## Coming from Python Control / SciPy

`control.ss`, `control.lqr`, `control.c2d` map as in the MATLAB table. Frequency-domain helpers allocate (`std::vector`) and live in `workbench.hpp`. There is no Python runtime; you generate C++ and flash the tick.

## I want to…

| I want to                                    | Start with                                                        | Example                                                                                               |
| -------------------------------------------- | ----------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------- |
| Regulate a fully measured plant              | `design::lqr` → `LQR`, `step(x)`                                  | [`examples/control/pendulum/`](../examples/control/pendulum/)                                         |
| Hold listed outputs at r, trusting the model | `TrackedOutputs` + `design::lqr` / `lqg`                          | [`examples/control/cart_pole/`](../examples/control/cart_pole/)                                       |
| Reject a constant disturbance / model offset | `design::lqi` / `lqgi`, or `design::oflqg`                        | [`examples/control/lqi/`](../examples/control/lqi/), [`heater_lqg/`](../examples/control/heater_lqg/) |
| Estimate the state from y                    | `design::lqg` or `design::kalman` → `KalmanFilter` / `FixedGainObserver` | [`examples/estimation/eskf/`](../examples/estimation/eskf/)                                           |
| Estimate a lumped disturbance from y, u      | `design::eso` → `ESO` (no plant A, B, C)                          | [`inc/damp/estimation/eso.hpp`](../inc/damp/estimation/eso.hpp)                                       |
| A SISO loop with ADRC (ESO + PD)             | `design::adrc` → `ADRCController`                                 | [`examples/control/adrc/`](../examples/control/adrc/)                                                 |
| Enforce input or output limits               | `design::mpc` / `ofmpc`, `back_calculate` after a clamp           | [`examples/control/heater_ofmpc/`](../examples/control/heater_ofmpc/)                                 |
| A SISO loop with a PID                       | `design::pidtune` → `PIDController`                               | [`examples/control/pid/`](../examples/control/pid/)                                                   |
| A nonlinear plant                            | `Dynamics` / `Measurement`, `design::linearize`, EKF/UKF          | [`examples/control/cart_pole/`](../examples/control/cart_pole/) SIL                                   |
| Frequency response, Bode, margins            | `workbench.hpp`                                                   | [`examples/control/workflow_end_to_end/`](../examples/control/workflow_end_to_end/)                   |

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
