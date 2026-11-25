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
    .Q = Matrix{
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
| `InputDelay`, dead time               | In the design, on the continuous plant with `Ts`: `design::state_feedback(plant, Ts, weights, …, InputDelay<N>{θ})`, and the same on `kalman`, `steady_state_target` and `mpc`, θ ≤ N·Ts (`InputDelay{θ}`: N = 1; one per input: `InputDelay<N, ColVec<NU>>{{θ₁, θ₂}}`). The plant's weights and noise, the cost exact on the delayed model, and a law that holds the delay states, the commands still on their way, itself (`update(x)` takes the plant's state) or carries them in its filter. Not a field on the law. Many samples: slower Ts or `SmithPredictor`. A discrete model with the delay built in, `discretize(sys, Ts, InputDelay<N>{θ})`, is designed on like any plant: its delay states are part of the measured state, and its noise names `G`. |
| `OutputDelay`, a late sensor          | In the filter's design, on the continuous plant with `Ts` and a `KalmanNoise`: `design::kalman(plant, Ts, noise, OutputDelay<N>{θ})`, and the same on `state_feedback` (LQG, LQGI, OFLQG), `steady_state_target` and `mpc` (an `MPC` through its filter, `OFMPC` among them). The reading that arrives now was taken θ ago; the filter holds the readings on their way in N delay states per output and estimates the current state, and the gain or the QP reads none of them. One delay per output: `OutputDelay<N, ColVec<NY>>{{θ₁, θ₂}}`; with an `InputDelay` too, both options. A discrete model with it built in, `discretize(sys, Ts, InputDelay{θu}, OutputDelay{θy})`, is designed on like any plant. |
| A tick that costs only its nonzero terms | `static constinit Law<kRes, float> law;` on a `constexpr` design: the same ticks, the matrices constants of the type, no RAM for them and no code for an exact 0 (MPC keeps its QP data in members, and only its filter, when it has one, is on constants). A member a redesign assigns is `Law<kRes, float>::Runtime`. On a `design::kalman`, `Law<kRes, float>::Recursive` is the recursive `KalmanFilter` on the same constants. |
| A Jacobian with a fixed shape, in a filter | `SparsityPattern` + `PatternMatrix<kPattern, float>`: `eskf.predict(F)`, `eskf.update_innovation(nu, H, R)`, `kf.update(y, C, R)`, or an EKF measurement whose `jacobian(x, u)` returns one. The zeros and ones are known at compile time, so they cost nothing; the values stay runtime. |
| `lqr`, `dlqr`, `lqrd`                 | `design::state_feedback(sys, …)` or `design::state_feedback(sys, Ts, …)` for a law; `matlab::lqr` / `dlqr` / `lqrd` (and `design::discrete_lqr`) for the bare gain |
| `lqi`, `lqg`, `lqgtrack`              | `design::state_feedback` with an `IntegralAction`, a `KalmanNoise`, or both                                      |
| `place`, `acker`                      | `design::place` for a gain (`matlab::place` / `acker`); `design::state_feedback(sys, [Ts,] design::Poles{…}, …)` for a law, the gain placed and every option of the LQ design (a target, limits, integrators, a filter, delays), on `StateFeedback`, the LQ runtime |
| `estim`, `reg` with placed gains      | `design::ObserverPoles{…}` in place of a `KalmanNoise` among `state_feedback`'s options: the law's filter placed, its error poles eig((I − MC)A) of the current form it runs; with `Poles` too for the regulator. `matlab::estim` / `reg` build the block from gains you hold |
| `pidtune`, Ziegler–Nichols            | `design::pidtune`, `design::ziegler_nichols_*`; then `.discretize(Ts)` → `PIDController`                         |
| `kalman`, `lqe`, `dlqe`, `kalmd`      | `design::kalman` → `FixedGainObserver` (fixed gain on the law) or `KalmanFilter` (P each tick)                          |
| observer, FixedGainObserver                  | `FixedGainObserver{design::kalman(...)}`                                                                                |
| ESO, ADRC                             | `design::eso` → `ESO` (observer only); `design::adrc` → `ADRCController` (ESO + PD)                              |
| anti-windup, sat, `umax`              | An `InputLimits` among an LQ or MPC design's options, `u_lim` on another design or its result: the law clamps to it and holds its memory against that clamp, with nothing else to set. `back_calculate(u_sat)` after a clamp downstream of the law. See [Limits and anti-windup](#limits-and-anti-windup) |
| bumpless transfer, manual mode, takeover | `law.track(u, …)` in place of `law.step(…)` on every LQ law and MPC while another source drives the plant, then `step` hands over from u; `design::Handover{.fade_time}` among `state_feedback`'s options for a law with nothing to preload. See [Bumpless handover](#bumpless-handover) |
| `mpcmove`, constraints                | `design::mpc` → `MPC`, on the state or through a filter with a `KalmanNoise`; `OFMPC`, the offset-free one, with a disturbance spec too. `MPC.step(r, m)` with a `TrackedOutputs{}` (every output) or a list; with neither a regulator, as `lqr` and `lqg`: `MPC.step(m)`, no r |
| `ManipulatedVariables` Min/Max, RateMin/RateMax | Min/Max → `design::InputLimits{.lower = {…}, .upper = {…}}`; RateMin/RateMax → `design::InputRateLimits{.lower = {…}, .upper = {…}}`, per sample. Among `mpc`'s options, hard rows on u and on Δu, the Δu rows only with an `InputRateLimits`; among `state_feedback`'s, a clamp on the command and on its move about the input applied last tick. Left out, unbounded |
| `OutputVariables` Min/Max, MinECR/MaxECR | `design::OutputLimits{.outputs = {i, …}, .lower = {…}, .upper = {…}}` among `mpc`'s options, the plant's outputs: rows for the ones it lists alone, soft by default (`.lower_ecr`, `.upper_ecr`; 0 is hard) |
| `Weights.OV`, a reference on some outputs | `MPCWeights::Qy` on every output of the plant; r on the outputs the target holds: a `TrackedOutputs{{i, …}}` among `mpc`'s options, or with a disturbance spec its `.outputs`. The others are weighed toward that steady state's values; `TrackedOutputs{}` holds every output, and with neither the MPC is a regulator with no r. At most one per input: with more references than inputs MATLAB® settles at a weighted least-squares compromise that holds none of them, and Damp asks which outputs r holds, so a plant with more outputs than inputs lists them (a longer list, or `TrackedOutputs{}`, does not compile) |
| `setterminal`, `Optimizer.MaxIterations` | `design::DareCost{}`, or `design::InvariantSet<NT_max>{.target_margin}` with room for NT_max rows; `design::QPIterations{n}` |
| `ManipulatedVariables` Target         | `Ru` weighs the input from the steady-state target's by default; `set_input_target(u)` on the runtime replaces it, `clear_input_target()` returns to it |
| `setmpcsignals` MD, a measured disturbance | `design::MeasuredDisturbance{.Bd = Bd}` among the options of `mpc`, `state_feedback` (every LQ law), `kalman` and `steady_state_target`, B_d continuous with `Ts`: the target is the steady state under d, `target_d` per unit [d̂; d], and the runtime takes it after the reference and its feedforward, `step(r, d, y)`. A filter takes d as a known input after u, so beside a disturbance spec d̂ estimates only the part it does not measure. A regulator holds every output at 0 under d, so on a plant with more outputs than inputs it does not compile: `TrackedOutputs{{i, …}}`, at most one per input, with r = 0 |
| `bode`, `nyquist`, `margin`, `rlocus` | `workbench.hpp`; `analysis::` or `matlab::bode`                                                                  |
| `lsim`, `step`, `forced_response`     | `sim::simulate*` on the host                                                                                     |
| ISR, float, no malloc                 | `control.hpp`, `.as<float>()`, `constinit`, `step` in the interrupt                                              |

## Coming from MATLAB®

Keep the plant as `StateSpace` (or `matlab::ss` / `tf`). Discretize with `discretize` or `matlab::c2d`. Design the law with `design::state_feedback` (`matlab::lqr` and `dlqr` give the bare gain, as MATLAB®'s do, which no runtime takes); deploy with `LQR law{res.as<float>()}` and `law.step(x)` — there is no `lsim` of a gain matrix on the MCU.

`[K,S,e] = lqr(sys,Q,R)` is `design::state_feedback(sys, {.Q, .R})` (CARE if continuous, DARE if discrete). Sampled-data (`lqrd`) is `design::state_feedback(sys, Ts, {.Q, .R})`.

Integral action is `design::state_feedback` with an `IntegralAction{.Q}`, its own cost on the integrators; it builds the augmented plant. Actuator saturation on that law is an `InputLimits{.lower, .upper}` among its options, and a slew limit an `InputRateLimits{.lower, .upper}`, which the law clamps to and holds its integrators against (a `tracking_time` in the `IntegralAction` unwinds them instead), or `back_calculate` after a clamp downstream. The same options bound an MPC's u and Δu, there over a horizon, with its output and state limits.

A design's options trail its call in any order, recognized by their type.

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
| Regulate a fully measured plant              | `design::state_feedback` → `StateFeedback` (an `LQR`), `step(x)`             | [`examples/control/pendulum/`](../examples/control/pendulum/)                                         |
| Hold listed outputs at r, trusting the model | `design::state_feedback` or `design::mpc` with `TrackedOutputs`    | [`examples/control/cart_pole/`](../examples/control/cart_pole/)                                       |
| Reject a constant disturbance / model offset | `design::state_feedback` with an `IntegralAction`, or a `KalmanNoise` and a disturbance spec | [`examples/control/lqi/`](../examples/control/lqi/), [`heater_lqg/`](../examples/control/heater_lqg/) |
| Feed a measured load forward                 | `design::state_feedback` or `design::mpc` with a `MeasuredDisturbance`, a disturbance spec beside it for the part it does not measure | [`tests/test_measured_disturbance.cpp`](../tests/test_measured_disturbance.cpp) |
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

- **TrackedOutputs + LQR / LQG / MPC** uses the model's steady-state maps $`x^* = N_x r`$, $`u^* = N_u r`$, those of `design::steady_state_target` with the same options. Offset-free only if the model is right and there is no unmatched load. An MPC weighs every output and pulls the ones r does not name toward $`C x^*`$.
- **LQI / LQGI** adds an integrator on $`r - C_i y`$. It rejects a constant disturbance and a constant model error. Use this when the plant will not match the design model.
- **OFLQG / OFMPC** estimate a disturbance and hold r through a steady-state target. Same job as LQI when a disturbance model fits; MPC when the actuator or output limits bind.

## The laws at a glance

What each runtime measures, the references and limits it takes, and what it does with a clamp after it (`back_calculate`). Every runtime has `valid()`: one built from a design it refuses sends 0 (an MPC holds its input), and when constant-evaluated it does not compile. Nothing faults at run time. A NaN or ±Inf parameter is refused like any other out of its domain, but a reading is not checked on the tick: a NaN measurement gives a NaN command, so check readings with `damp::isfinite` before the tick ([known limitations](known_limitations.md)).

### State-space laws

The five LQ rows are one runtime, `StateFeedback`, which deduces from any LQ design (`StateFeedback law{res}`, or `Law<kRes, T>`) and has the members its design carries. `LQR`, `LQG`, `LQI`, `LQGI` and `OFLQG` name what the options produce: each is an alias of it that refuses another kind's design, so `LQG law{res}` of an LQI's does not compile. The two MPC rows are one runtime too, `MPC`, which deduces from any `design::mpc` (on the measured state, or through a filter with a `KalmanNoise`); `OFMPC` is an alias of it that refuses a design whose filter estimates no disturbance.

| Runtime | Measures | References | Limits | Clamp after the law |
|---|---|---|---|---|
| `LQR` | state x: `update(x)` | none (regulator); `x_ref`/`u_ref`; `r` on tracked outputs (`TrackedOutputs`); measured disturbance `d` | `input`, `input_rate`, `demand` | `back_calculate`; `predict` shifts the `InputDelay`'s delay states |
| `LQG` | y, through a Kalman filter: `update(y)` | none (regulator); `x_ref`/`u_ref`; `r` on tracked outputs (`TrackedOutputs`); measured disturbance `d` | `input`, `input_rate`, `demand` | `back_calculate`; the filter runs on the applied input |
| `LQI` | state x: `update(x)` | `r` on the integrated outputs, plus optional `x_ref`/`u_ref`; measured disturbance `d` | `input`, `input_rate`, `integral`, `demand` | ξ holds against the input and move clamps by itself, or unwinds over `tracking_time`; `back_calculate` |
| `LQGI` | y, through a Kalman filter: `update(y)` | `r` on the integrated outputs, plus optional `x_ref`/`u_ref`; measured disturbance `d` | `input`, `input_rate`, `integral`, `demand` | ξ holds against the input and move clamps by itself, or unwinds over `tracking_time`; `back_calculate`; the filter runs on the applied input |
| `OFLQG` | y, through a Kalman filter on [x; d]: `update(y)` | `r` on the controlled outputs, plus optional `x_ref`/`u_ref`; measured disturbance `d` | `input`, `input_rate`, `demand`, `disturbance` | `back_calculate`; the filter runs on the applied input |
| `MPC` | state x: `update(x)`; or y, through a Kalman filter with a `KalmanNoise`: `update(y)` | `r` held or previewed, on the outputs the target holds (`TrackedOutputs`, else every output), plus optional `x_ref`/`u_ref`; measured disturbance `d`; input target | `input`, `input_rate`, `output`, `state` | `back_calculate`; `last_status()`; a filter runs on the applied input |
| `OFMPC` (the `MPC` of an offset-free design) | y, through a Kalman filter on [x; d]: `update(y)` | `r` held or previewed, on the outputs the disturbance spec lists (else every output), plus optional `x_ref`/`u_ref`; measured disturbance `d`; input target | `input`, `input_rate`, `output`, `state`, `disturbance` | `back_calculate`; `last_status()`; the filter runs on the applied input |

Every law in this table reads its run-time limits with `limits()` and sets them with `set_limits()`, the fields its Limits column lists: an LQ law's `StateFeedbackLimits` (`input`, `input_rate`, `integral`, `demand`, `disturbance`), an MPC's `MPCLimits` (`input`, `input_rate`, `output`, `state`, `disturbance`), a field the law does not have empty: `input_rate` without an `InputRateLimits`, `integral` without integrators, `disturbance` without a disturbance model. A design's result holds the limits its law starts with in the same struct, `res.limits`. `auto lim = mpc.limits(); lim.output.upper(0) = v_max; mpc.set_limits(lim);` On `MPC` and `OFMPC`, the QP weighs every output of the plant, `output` and `state` bound the outputs an `OutputLimits` option lists and the states a `StateLimits` lists, soft by default; the terminal ingredient (`DareCost`, `InvariantSet`) is part of the design, and with an `InvariantSet` `set_limits()` refuses a change of `input`, `input_rate` or `output`, which the set was built for.

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

- **The output limit** is a `Bounds` (`{.lower = -12.0, .upper = 12.0}`), unbounded by default: an `InputLimits{.lower = {-12.0}, .upper = {12.0}}` among an LQ or MPC design's options, and `u_lim`, a field of the design input where the design takes a struct, or of the result, on the others. A limit that is not a box (an end not finite, or lower above upper) is refused. A rail that moves is set per tick: the LQ laws and the MPCs read their run-time limits with `limits()` and set them with `set_limits(lim)`, which refuses one that is not a box and changes nothing (`auto l = law.limits(); l.input = {…}; law.set_limits(l);`); the other laws hold `u_lim` as a public member.
- **The tick returns a command inside its limits**, and the law's memory holds against that clamp with nothing else to set: an integral does not run on behind it, and an observer or internal model takes the clamped command as the input applied.
- **`back_calculate(u_sat, u_ff = {})`** reports a clamp after the law: `u_sat` is what the plant received, `u_ff` what was added to the law's command before the clamp. The law's memory holds against it the same way.

| Law | Output limit | What holds against the clamp |
| --- | --- | --- |
| `PIDController`, `PIController`, `ContinuousPID` | `u_lim` in `PIDGains` | The integral: it gives back what the tick added, up to what the clamp cut off. With a `Kbc` it unwinds by back-calculation instead |
| `LQR`, `LQI`, `LQG`, `LQGI`, `OFLQG` | An `InputLimits` and an `InputRateLimits`, `limits` in the `IntegralAction` and a `DemandLimit`; at run time `set_limits()` | The integrators: they hold (conditional integration), or with a `tracking_time` in the `IntegralAction` unwind by back-calculation; the filter and the delay states move on the clamped command; `OFLQG`'s d̂ stays in `disturbance_lim` |
| `MPC`, `OFMPC` | An `InputLimits`, an `InputRateLimits`, an `OutputLimits` and a `StateLimits`; at run time `set_limits()` | By construction: the limits are constraints of the QP; `OFMPC`'s d̂ stays in `disturbance` |
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

The integrators' limit bounds the integral a law stores, so its unit is that law's: the I term in output units on a PID (`integral_lim`), ξ ≈ ∫e·dt on the LQ laws (`IntegralAction::limits`, at run time `integral`). `disturbance` does the same for the d̂ of an `OFLQG` or an `OFMPC`, in d's units: the runtime's alone, set on the law with `set_limits()`. A move limit, an `InputRateLimits`, is per sample on both families; a `SlewLimiter`'s rates are per second.

A limit on a state or an output needs a prediction: `design::mpc` (an `OutputLimits` on any output, a `StateLimits` on any state, soft or hard), `design::mhe` (`x_lim` on the estimate), or an `ActionGovernor` around any law. A `DemandLimit` among an LQ design's options is not one: it bounds what the gain's outer terms ask of a state, read as a cascade inside one gain (the lean a balancing robot's position and integral ask for). The integrators hold against it as against the input limits, or unwind with a `tracking_time`. It is a heuristic on the gain: nothing certifies the loop while the clamp binds. [demand_lim.md](demand_lim.md) derives it from first principles, with its limits and the published work nearest to it.

References follow the same split. Every law takes the output reference `r`. The state-feedback laws and the MPCs also take a reference point inside the law, `control(r, x_ref, u_ref)` and `step(r, x_ref, u_ref, y)`, which the clamp covers; with a measured disturbance `control(r, x_ref, u_ref, d)` and `step(r, x_ref, u_ref, d, x)`, the reference and its feedforward first, then what was read (a regulator, which has no r, `control(d)` and `step(d, x)`). An MPC poses its QP about that point and carries it on the model over the horizon, as the LQ gain does, but its output and state limits stay on the absolute outputs and states: a limit that should ride the reference is moved at run time with `set_limits()`. A feedforward added outside any law is `u_ff`, passed to `back_calculate` so it is not taken for saturation. MPC takes `r` or a `Preview` of it, each step's target from it, and weighs the input from the target's: an input preference replaces that through `set_input_target`, and `clear_input_target` restores it.

Not covered yet: `STSMCController`'s integral state has no limit of its own, and `SwitchedController` starts a law without a tracking mode (anything but a PID) from its reset state.

## Bumpless handover

While another source drives the plant (an operator, another law), call the law's `track(u, …)` in place of `step(…)`: each `step` overload has a twin with step's arguments and the applied input first. It predicts and updates as `step` does, then sets the law's memory so that its command at this state, with these arguments, would be u. There is no mode flag: the law tracks for as long as `track` is called, and the first `step` after it starts from u.

```cpp
law.track(u_applied, r, y); // in place of step(r, y) while another source drives
u = law.step(r, y);         // handover: the first command starts from u_applied
```

| Law | What `track` leaves |
| --- | --- |
| `LQI`, `LQGI` | ξ preloaded through (∂u/∂ξ)⁺, back-calculation's map: exact with as many integrators as inputs, the least-squares ξ with fewer, then clamped into the integrators' limits (`IntegralAction::limits`) |
| `LQG`, `LQGI`, `OFLQG`, `OFMPC` | the filter predicts with u, and the measured d: an offset-free filter learns the load, so at a steady state the first command is u |
| any law with an `InputDelay` | the delay states hold u |
| any law with an `InputRateLimits` | the move's base is u: the first command moves at most the rate from u, on an LQ law as on an MPC |
| `MPC`, `OFMPC` | u_prev = u: `Rdu` and, with an `InputRateLimits`, the Δu box shape the first moves. No QP is solved while tracking, and the solver's warm start stays as the last solve left it |
| `LQR`, `LQG`, an `OFLQG` away from a steady state | nothing to preload: the first step is the law's own command. With `design::Handover{.fade_time}` among the design's options it is u, and the mismatch e = u − u_law decays by exp(−Ts/fade_time) each tick, within the input limits; the integrators take e as a feedforward |

A u outside the input limits is one the law cannot emit: an LQ law hands over from u clamped into them, while the filter and the delay states take the u the plant received. `track` adds nothing to `step`, and a law without a `Handover` is the type and size it was. A law on the low-latency tick calls `predict()` once after the last `track`, as `step` does at its start.

## Tick

`step(r, m)` is predict, then update(m), then control. Low-latency ISR: update → control → write → `back_calculate` → predict. After a clamp, pass the applied input to `back_calculate(u_sat, u_ff)`. While another source drives the plant, `track(u, r, m)` takes the place of `step(r, m)` ([Bumpless handover](#bumpless-handover)).

Readings, references and reported inputs are the application's to check, with `damp::isfinite` (on the bits, so it holds under `-ffast-math`): what to do when a sensor fails, stop, hold or switch to a backup, is the application's call. A NaN or ±Inf that reaches a law stays in its memory until that is reset: `reset()` clears the integrators, a filter restarts with `estimator.reset(x0)`, and `law = decltype(law){}` starts the whole law over.

## Add the library

Headers only: `-I path/to/inc` and `#include "damp/control.hpp"`. CMake: link `wetmelon::damp` (FetchContent or `find_package`). Host analysis: `damp/workbench.hpp`.
