# Heater with dead time — offset-free LQG, the dead time in one sample

The fewest lines that design and step an offset-free LQG on a plant with a transport delay. The physics, the fitted numbers and the delay handling are in [`../heater_ofmpc/`](../heater_ofmpc/heater_ofmpc_derivation.md); this folder is the same hotend under the unconstrained law that MPC reduces to. Numbers here match [`heater_lqg_controller.hpp`](heater_lqg_controller.hpp) and [`heater_lqg_sil.cpp`](heater_lqg_sil.cpp).

---

## The law

Take the horizon and the power limits away from the MPC and what is left is an offset-free LQG: an LQR on the same model with the same output weight, the same disturbance-augmented Kalman filter, and a steady-state target computed from the estimated heat loss $`\hat w`$ so that $`y \to r`$ under it (`design::oflqg(kModel, weights, noise, design::InputDisturbance{})`, run as `OFLQG`). It is the law the MPC applies whenever no constraint is active.

```cpp
constexpr auto kModel = *discretize(kPlant, kTs, InputDelay{kDeadTime}); // [T; P(k−1)]: the dead time is one period
constexpr auto kLaw = design::oflqg(kModel, {.Q = Q, .R = R}, noise, design::InputDisturbance{});
using Heater = Law<kLaw, float>; // the runtime, named by the result: OFLQG
```

`kPlant` is the continuous block. `discretize` samples it at `kTs` with the dead time in it: the power decided now reaches the block over the next period, so the model is the block temperature and one register state, the power of the last period. The design is on that discrete model, so its numbers are per period: `Q` and `R` weight each sampled temperature and power, and the process and disturbance noise are the covariances added per period (0.01 and 0.05 per second over 5 s). The measurement noise is per reading, the thermistor's variance as the firmware samples it. `Q` weights the block temperature, `R = 0.02` the power itself: an LQR has no rate weight, where the MPC's `Rdu` is on the power change, so the two are tuned to comparable aggressiveness rather than identical cost. The filter's process noise is on the block alone; the register state is a past power, known exactly.

The control period is the dead time, 5 s, so the power decided now reaches the block over the next period: one register state holds it in the model. The model holds at most one sample of delay; a longer dead time means sampling slower, which a 30 s thermal plant tolerates (six samples per time constant), or, for a single loop, `SmithPredictor`. θ can be any value up to the period, a value fitted at commissioning, say; the runtime type is the same for all of them.

The law knows nothing of $`0 \le P \le P_{\max}`$. The tick clamps its request at the actuator and, through `update(y)`, `control(r)`, then `back_calculate(u_sat)` and `predict()` with the clamped power, tells the estimator what was actually applied; without that the filter would believe in heat that never went in and the estimate of $`w`$ would wind up during the heat-up.

Tick in code: `heater_period(controller, setpoint, temperature)`, the same function in the sketch and the SIL. The sketch flashes a `constinit` controller designed at compile time.

---

## What the SIL shows

The hotend with its true 5 s transport delay, a room 3 °C colder than the controller assumes, and a 5 W part-cooling fan from 240 s. Cold start to 210 °C.

| Heat-up overshoot | Last-minute error under the fan | Peak power |
| --- | --- | --- |
| 0.01 °C | 0.000 °C | 40 W |

On this plant the clamped LQG matches the MPC of `heater_ofmpc/`: only full power during the heat-up ever binds, and a clamp with the estimator told about it handles that one constraint as well as the QP. The MPC earns its keep where constraints bind during normal operation, or where there are several of them.

The SIL fails on an overshoot above 0.5 °C, on an offset above 0.1 °C under the fan, or on power above 40 W.

---

## Forking it

- Put the fitted $`K`$, $`\tau`$, $`\theta`$ and $`T_s`$ in the nameplate and the model takes the dead time in seconds, `discretize(kPlant, kTs, InputDelay{kDeadTime})`, as long as it is at most `kTs`. The runtime type follows.
- A trajectory to track rather than a setpoint to hold is the `step({}, x_ref, u_ref, y)` form: `x_ref` is the model state along it, the block temperature and the last period's power $`[T; P(k-1)]`$, `u_ref` the power, and `r = {}` because the reference point carries the trajectory.
