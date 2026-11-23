# Using Damp

Short paths from a plant to a tick. Exhaustive names: [REFERENCE.md](../REFERENCE.md).
Rules: [AGENTS.md](../AGENTS.md). Limits: [known_limitations.md](known_limitations.md).
Compiler bugs: [known_issues.md](known_issues.md).

The plant is the process under control (motor, heater, pendulum). Damp holds a
model of it: command $`u`$ in, measurement $`y`$ (and often state $`x`$) out —
`StateSpace`, a transfer function, or `Dynamics`.

## Canonical path

[`examples/control/lqi/`](../examples/control/lqi/) is the one to copy first:
`StateSpace` in the coordinates you measure → `discretize` → `design::state_feedback` with an `IntegralAction{.Q}` →
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
constexpr auto sys = StateSpace{
    .A = {{0.0, 1.0}, {-k / m, -c / m}},
    .B = {{0.0}, {1.0 / m}},
    .C = {{1.0, 0.0}},
};
} // namespace

// Configure tuning parameters (control weights) for an LQI
constexpr auto weights = design::LQRWeights{
    .Q = Matrix<2, 2>{
        {10.0, 0.0},
        {0.0, 1.0},
    },
    .R = {1.0},
};

// Discretize the plant and solve the DARE
constexpr auto res = design::state_feedback(sys, Ts, weights, design::IntegralAction{.Q = 80.0}); // on the position

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
| `LQR law{design::state_feedback(plant, Ts, w).as<float>()};` | On the MCU, at init | Members, in RAM | The result's `success` is false, and the law keeps zero gains | Dense, plus the design at init |

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
| `InputDelay`, dead time               | `discretize(sys, Ts, InputDelay<N>{θ})` with θ ≤ N·Ts (`InputDelay{θ}`: N = 1). One delay per input: `InputDelay<N, ColVec<NU>>{{θ₁, θ₂}}`. Many samples: slower Ts or `SmithPredictor`. Not a field on the law. Every sampled LQ design and `kalman` take it in the design, `design::state_feedback(plant, Ts, weights, noise, InputDelay<N>{θ}, design::IntegralAction{})`: the plant's weights and noise, the cost exact on the delayed model, and a law that holds or carries the registers itself. MPC is designed on the model. |
| `OutputDelay`, a late sensor          | `discretize(sys, Ts, OutputDelay<N>{θ})`: the reading that arrives now was taken θ ago, held in N registers. One delay per output: `OutputDelay<N, ColVec<NY>>{{θ₁, θ₂}}`. Both: `discretize(sys, Ts, InputDelay{θu}, OutputDelay{θy})`. A filter on the model estimates the current state. |
| A tick that costs only its nonzero terms | `static constinit Law<kRes, float> law;` on a `constexpr` design: the same ticks, the matrices constants of the type, no RAM for them and no code for an exact 0 (MPC keeps its QP data in members, and OFMPC's filter alone is on constants). A member a redesign assigns is `Law<kRes, float>::Runtime`. On a `design::kalman`, `Law<kRes, float>::Recursive` is the recursive `KalmanFilter` on the same constants. |
| A Jacobian with a fixed shape, in a filter | `SparsityPattern` + `PatternMatrix<kPattern, float>`: `eskf.predict(F)`, `eskf.update_innovation(nu, H, R)`, `kf.update(y, C, R)`, or an EKF measurement whose `jacobian(x, u)` returns one. The zeros and ones are known at compile time, so they cost nothing; the values stay runtime. |
| `lqr`, `dlqr`, `lqrd`                 | `design::state_feedback(sys, …)` or `design::state_feedback(sys, Ts, …)`; `matlab::lqr` / `dlqr` / `lqrd`                              |
| `lqi`, `lqg`, `lqgtrack`              | `design::state_feedback` with an `IntegralAction`, a `KalmanNoise`, or both                                      |
| `place`, `acker`                      | `design::place`; `matlab::place` / `acker`                                                                       |
| `pidtune`, Ziegler–Nichols            | `design::pidtune`, `design::ziegler_nichols_*`; then `.discretize(Ts)` → `PIDController`                         |
| `kalman`, `lqe`, `dlqe`, `kalmd`      | `design::kalman` → `FixedGainObserver` (fixed gain on the law) or `KalmanFilter` (P each tick)                          |
| observer, FixedGainObserver                  | `FixedGainObserver{design::kalman(...)}`                                                                                |
| ESO, ADRC                             | `design::eso` → `ESO` (observer only); `design::adrc` → `ADRCController` (ESO + PD)                              |
| anti-windup, sat, `umax`              | `u_lim` on the design or its result: the law clamps to it and holds its memory against that clamp, with nothing else to set. `back_calculate(u_sat)` after a clamp downstream of the law. See [Limits and anti-windup](#limits-and-anti-windup) |
| `mpcmove`, constraints                | `design::mpc` (`MPC`, or `OFMPC` with a `KalmanNoise` and a disturbance spec), `MPC.step(r, y)`                  |
| `bode`, `nyquist`, `margin`, `rlocus` | `workbench.hpp`; `analysis::` or `matlab::bode`                                                                  |
| `lsim`, `step`, `forced_response`     | `sim::simulate*` on the host                                                                                     |
| ISR, float, no malloc                 | `control.hpp`, `.as<float>()`, `constinit`, `step` in the interrupt                                              |

## Coming from MATLAB®

Keep the plant as `StateSpace` (or `matlab::ss` / `tf`). Discretize with `discretize` or `matlab::c2d`. Design with `matlab::lqr` if you want that spelling; deploy with `LQR law{res.as<float>()}` and `law.step(x)` — there is no `lsim` of a gain matrix on the MCU.

`[K,S,e] = lqr(sys,Q,R)` is `design::state_feedback(sys, {.Q, .R})` (CARE if continuous, DARE if discrete). Sampled-data (`lqrd`) is `design::state_feedback(sys, Ts, {.Q, .R})`.

Integral action is `design::state_feedback` with an `IntegralAction{.Q}`, its own cost on the integrators; it builds the augmented plant. Actuator saturation on that law is an `LQRLimits{.u_lim}` among its options, which the law clamps to and holds its integrators against (a `tracking_time` in the `IntegralAction` unwinds them instead), or `back_calculate` after a clamp downstream. Input and output limits over a horizon are `design::mpc`.

A design's options trail its call in any order, recognized by their type. For a caller who would name each one, or keep a design's inputs in one object, `damp/design/config.hpp` is a helper over them: `design::state_feedback(sys, Ts, design::LQConfig{.weights = weights, .integral = design::IntegralAction{.Q = 80.0}})` is the same design, and `KalmanConfig`, `MPCConfig` and `TargetConfig` do the same for `kalman`, `mpc` and `steady_state_target`. Every field holds a named object (`design::LQRWeights{…}`, not bare braces), which the plant cannot size inside a config.

## Coming from firmware

One object, one tick. Design at compile time or init. Do not run Riccati in the PWM interrupt.

```cpp
#include "damp/control.hpp"
constinit Law<kRes, float> law; // kRes from design::state_feedback with an IntegralAction
// ISR:
auto u = law.step(r, x);           // or update / control / back_calculate / predict
write_pwm(u[0]);
```

PID: `design::pidtune` or `design::pid(Kp, Ki, Kd)` → `.discretize(Ts)` → `PIDController`. Limits go on the gains (`u_lim`), not as extra ISR code, unless a downstream stage clamps — then `back_calculate`. The output limit is all the anti-windup a PID needs: the integral holds against it.

## Coming from Python Control / SciPy

`control.ss`, `control.lqr`, `control.c2d` map as in the MATLAB table. Frequency-domain helpers allocate (`std::vector`) and live in `workbench.hpp`. There is no Python runtime; you generate C++ and flash the tick.

## I want to…

| I want to                                    | Start with                                                        | Example                                                                                               |
| -------------------------------------------- | ----------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------- |
| Regulate a fully measured plant              | `design::state_feedback` → `LQR`, `step(x)`                                  | [`examples/control/pendulum/`](../examples/control/pendulum/)                                         |
| Hold listed outputs at r, trusting the model | `design::state_feedback` with `TrackedOutputs`                     | [`examples/control/cart_pole/`](../examples/control/cart_pole/)                                       |
| Reject a constant disturbance / model offset | `design::state_feedback` with an `IntegralAction`, or a `KalmanNoise` and a disturbance spec | [`examples/control/lqi/`](../examples/control/lqi/), [`heater_lqg/`](../examples/control/heater_lqg/) |
| Estimate the state from y                    | `design::state_feedback` with a `KalmanNoise` or `design::kalman` → `KalmanFilter` / `FixedGainObserver` | [`examples/estimation/eskf/`](../examples/estimation/eskf/)                                           |
| Estimate a lumped disturbance from y, u      | `design::eso` → `ESO` (no plant A, B, C)                          | [`inc/damp/estimation/eso.hpp`](../inc/damp/estimation/eso.hpp)                                       |
| A SISO loop with ADRC (ESO + PD)             | `design::adrc` → `ADRCController`                                 | [`examples/control/adrc/`](../examples/control/adrc/)                                                 |
| Enforce input or output limits               | `design::mpc`, `back_calculate` after a clamp                     | [`examples/control/heater_ofmpc/`](../examples/control/heater_ofmpc/)                                 |
| A SISO loop with a PID                       | `design::pidtune` → `PIDController`                               | [`examples/control/pid/`](../examples/control/pid/)                                                   |
| A nonlinear plant                            | `Dynamics` / `Measurement`, `design::linearize`, EKF/UKF          | [`examples/control/cart_pole/`](../examples/control/cart_pole/) SIL                                   |
| Frequency response, Bode, margins            | `workbench.hpp`                                                   | [`examples/control/workflow_end_to_end/`](../examples/control/workflow_end_to_end/)                   |

Copy an example folder. SIL calls the same tick as the sketch.

## LQR with TrackedOutputs, or LQI

Both can hold an output at r.

- **TrackedOutputs + LQR / LQG** uses the model's steady-state maps $`x^* = N_x r`$, $`u^* = N_u r`$. Offset-free only if the model is right and there is no unmatched load.
- **LQI / LQGI** adds an integrator on $`r - C_i y`$. It rejects a constant disturbance and a constant model error. Use this when the plant will not match the design model.
- **OFLQG / OFMPC** estimate a disturbance and hold r through a steady-state target. Same job as LQI when a disturbance model fits; MPC when the actuator or output limits bind.

## The laws at a glance

What each runtime measures, the references and limits it takes, and what it does with a clamp after it (`back_calculate`). Every runtime has `valid()`: one built from a design it refuses sends 0 (an MPC holds its input), and when constant-evaluated it does not compile. Nothing faults at run time.

### State-space laws

| Runtime | Measures | References | Limits | Clamp after the law |
|---|---|---|---|---|
| `LQR` | state x: `update(x)` | none (regulator); `x_ref`/`u_ref`; `r` on tracked outputs (`TrackedOutputs`) | `u_lim`, `demand_lim` | `back_calculate`; `predict` shifts the `InputDelay` registers |
| `LQG` | y, through a Kalman filter: `update(y)` | none (regulator); `x_ref`/`u_ref`; `r` on tracked outputs (`TrackedOutputs`) | `u_lim`, `demand_lim` | `back_calculate`; the filter runs on the applied input |
| `LQI` | state x: `update(x)` | `r` on the integrated outputs, plus optional `x_ref`/`u_ref` | `u_lim`, `demand_lim`, `integral_lim` | ξ holds against the `u_lim` clamp by itself, or unwinds over `tracking_time`; `back_calculate` |
| `LQGI` | y, through a Kalman filter: `update(y)` | `r` on the integrated outputs, plus optional `x_ref`/`u_ref` | `u_lim`, `demand_lim`, `integral_lim` | ξ holds against the `u_lim` clamp by itself, or unwinds over `tracking_time`; `back_calculate`; the filter runs on the applied input |
| `OFLQG` | y, through a Kalman filter on [x; d]: `update(y)` | `r` on the controlled outputs, plus optional `x_ref`/`u_ref` | `u_lim`, `demand_lim`, `disturbance_lim` | `back_calculate`; the filter runs on the applied input |
| `MPC` | state x: `update(x)` | `r` held or previewed; measured disturbance `d`; input target | `u_lim`, `du_lim`, `y_lim`, `x_lim` | `back_calculate`; `last_status()` |
| `OFMPC` | y, through a Kalman filter on [x; d]: `update(y)` | `r` held or previewed; input target | `u_lim`, `du_lim`, `y_lim`, `x_lim`, `disturbance_lim` | `back_calculate`; `last_status()`; the filter runs on the applied input |

On `MPC` and `OFMPC`, `y_lim` and `x_lim` (on the states listed in `x_states`) are soft by default, the terminal set is part of the design, and `set_limits()` moves `u_lim`, `du_lim`, `y_lim` and `x_lim` at run time.

### SISO laws

| Runtime | Tick | Reference | Limits | Clamp after the law |
|---|---|---|---|---|
| `PIDController` (PID, PI) | `control(r, y)` | `r`, with setpoint weights `b`, `c` | `u_lim`, `integral_lim` | `back_calculate` (`Kbc`; 0 means I holds), `track(u, y)`, `disable(u)` / `enable()` |
| `PDController` | `control(r, y)` | `r`, with setpoint weights `b`, `c` | `u_lim` | `back_calculate` (no-op); `disable(u)` / `enable()` |
| `PController` | `control(r, y)` | `r`, with setpoint weight `b` | `u_lim` | `back_calculate` (no-op); `disable(u)` / `enable()` |
| `ContinuousPID` | `control(r, y, Ts)` | `r`, with setpoint weights `b`, `c` | `u_lim`, `integral_lim` | `back_calculate`, `track(u, y, Ts)`, `disable(u)` / `enable()` |
| `LeadLagController` | `control(e)` or `control(r, y)` | error | `u_lim` | `back_calculate` |
| `PRController` | `control(e)` or `control(r, y)` | error | `u_lim` | resonant section holds (`tracking_time` = 0) or unwinds over `tracking_time`; `back_calculate`; `set_frequency` |
| `HarmonicSuppressor` | `control(e)` or `control(r, y)` | error | `u_lim` | resonators hold, fundamental first; `back_calculate`; `set_fundamental` |
| `RepetitiveController` | `control(e)` | error | none (an add-on: the clamp is on the sum) | `back_calculate` |
| `SmithPredictor` | `control(r, y)` | `r` | the primary law's `u_lim` | `back_calculate` (the model and the primary law) |
| `ADRCController` | `control(r, y)` | `r` | `u_lim` | the observer runs on the clamped command; `back_calculate` |
| `SMCController` | `control(r, y)` | `r` | `u_lim` (output also within ±k/b₀) | `back_calculate` (no-op) |
| `STSMCController` | `control(r, y)` or `control_sliding(s)` | `r`, or the sliding variable s | `u_lim` | v holds against the clamp; `back_calculate` |
| `ExtremumSeekingController` | `control(J, valid)` | none | `u_lim` on û and on û + dither | `back_calculate` |

### Composition and safety filter

| Runtime | Tick | Limits | Clamp after the law |
|---|---|---|---|
| `Cascade` | `control(r, y_outer, y_inner)` | each loop's own `u_lim` | `back_calculate`; the outer law holds while the inner command is clamped |
| `SwitchedController` | `set_mode`, `control(r, y)`, `inject(u)` | each law's own `u_lim` | `back_calculate` goes to the active law; a PID switched in starts from the last applied command, others from reset |
| `ActionGovernor` | `filter(u_des)` | A·u ≤ b (box, affine or CBF rows); `set_bounds(b)` each tick | if the projection fails it sends the last safe u |

## Limits and anti-windup

One contract for every law:

- **`u_lim`** is the output limit, a `Bounds` (`{.lower = -12.0, .upper = 12.0}`), unbounded by default. It is a field of the design input where the design takes a struct, and of the result otherwise. The runtime holds it as a public member, so a rail that moves is set per tick. A limit that is not a box (an end not finite, or lower above upper) is refused.
- **The tick returns a command inside `u_lim`**, and the law's memory holds against that clamp with nothing else to set: an integral does not run on behind it, and an observer or internal model takes the clamped command as the input applied.
- **`back_calculate(u_sat, u_ff = {})`** reports a clamp after the law: `u_sat` is what the plant received, `u_ff` what was added to the law's command before the clamp. The law's memory holds against it the same way.

| Law | Output limit | What holds against the clamp |
| --- | --- | --- |
| `PIDController`, `PIController`, `ContinuousPID` | `u_lim` in `PIDGains` | The integral: it gives back what the tick added, up to what the clamp cut off. With a `Kbc` it unwinds by back-calculation instead |
| `LQR`, `LQI`, `LQG`, `LQGI`, `OFLQG` | `u_lim` in an `LQRLimits`, and `lim` in the `IntegralAction` | The integrators: they hold (conditional integration), or with a `tracking_time` in the `IntegralAction` unwind by back-calculation; the filter and the input registers move on the clamped command; `OFLQG`'s d̂ stays in `disturbance_lim` |
| `MPC`, `OFMPC` | `u_lim`, `du_lim`, `y_lim`, `x_lim` in `MPCConstraints` | By construction: the limits are constraints of the QP; `OFMPC`'s d̂ stays in `disturbance_lim` |
| `ADRCController` | `u_lim` on `ADRCResult` | The ESO runs on the clamped command |
| `STSMCController` | `u_lim` on `STSMCResult` | The integral state |
| `LeadLagController` | `u_lim` on `LeadLagResult` | The recursion runs on the clamped output |
| `PRController` | `u_lim` on `PRResult` | The resonant section gives back what the tick added (`tracking_time` = 0), or unwinds over `tracking_time` |
| `HarmonicSuppressor` | `u_lim` on `HarmonicSuppressorResult` | The resonators, fundamental first, never more in all than the clamp cut |
| `SMCController` | `u_lim` on `SMCResult` | Nothing to hold: it keeps only the last error |
| `ExtremumSeekingController` | `u_lim` in `ESCSpec` | The operating point stays inside, and follows a clamp reported to it |
| `SmithPredictor` | its primary's | The primary's memory; the internal model runs on the input reported |
| `RepetitiveController` | none: it returns a correction | A sample the clamp cut off is not learned |
| `ClassicalDOB` | none: it returns the corrected command | The estimate takes the input reported |
| `Cascade` | each loop's own | The outer law holds at the last command the inner loop could follow |

`integral_lim` bounds the integral a law stores, so its unit is that law's: the I term in output units on a PID, ξ ≈ ∫e·dt on the LQ laws. `disturbance_lim` does the same for the d̂ of an `OFLQG` or an `OFMPC`, in d's units: a member of the runtime only, set on the law.

A limit on a state or an output needs a prediction: `design::mpc` (`y_lim` on the controlled outputs, `x_lim` on any state, soft or hard), `design::mhe` (`x_lim` on the estimate), or an `ActionGovernor` around any law. A `DemandLimit` among an LQ design's options is not one: it bounds what the gain's outer terms ask of a state, read as a cascade inside one gain (the lean a balancing robot's position and integral ask for). The integrators hold against it as against `u_lim`, or unwind with a `tracking_time`. It is a heuristic on the gain: nothing certifies the loop while the clamp binds. [demand_lim.md](demand_lim.md) derives it from first principles, with its limits and the published work nearest to it.

References follow the same split. Every law takes the output reference `r`. The state-feedback laws also take a reference point inside the law, `control(r, x_ref, u_ref)`, which the clamp covers. A feedforward added outside any law is `u_ff`, passed to `back_calculate` so it is not taken for saturation. MPC takes `r` or a `Preview` of it, and an input preference through `set_input_target`.

Not covered yet: `STSMCController`'s integral state has no limit of its own, and `SwitchedController` starts a law without a tracking mode (anything but a PID) from its reset state.

## Tick

`step(r, m)` is predict, then update(m), then control. Low-latency ISR: update → control → write → `back_calculate` → predict. After a clamp, pass the applied input to `back_calculate(u_sat, u_ff)`.

Readings, references and reported inputs are the application's to check, with `damp::isfinite` (on the bits, so it holds under `-ffast-math`): what to do when a sensor fails, stop, hold or switch to a backup, is the application's call. A NaN or ±Inf that reaches a law stays in its memory until that is reset: `reset()` clears the integrators, a filter restarts with `estimator.reset(x0)`, and `law = decltype(law){}` starts the whole law over.

## Add the library

Headers only: `-I path/to/inc` and `#include "damp/control.hpp"`. CMake: link `wetmelon::damp` (FetchContent or `find_package`). Host analysis: `damp/workbench.hpp`.
