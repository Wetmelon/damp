# Heater with dead time — offset-free LQG in the predictor form

The fewest lines that design and step an offset-free LQG on a plant with a transport delay. The physics, the fitted numbers and the delay handling are in [`../heater_ofmpc/`](../heater_ofmpc/heater_ofmpc_derivation.md); this folder is the same hotend under the unconstrained law that MPC reduces to. Numbers here match [`heater_oflqg_controller.hpp`](heater_oflqg_controller.hpp) and [`heater_oflqg_sil.cpp`](heater_oflqg_sil.cpp).

---

## The law

Take the horizon and the power limits away from the MPC and what is left is an offset-free LQG: an LQR on the same model with the same output weight, the same disturbance-augmented Kalman filter, and a steady-state target computed from the estimated heat loss $`\hat w`$ so that $`y \to r`$ under it (`design::discrete_oflqg`). It is the law the MPC applies whenever no constraint is active.

```cpp
constexpr auto kLaw = design::discrete_oflqg<kDeadTime / kTs>(kPlant, kTs, Q, R, noise); // 2.5 periods: the block and one register tap
using Heater = decltype(OFLQG{kLaw.as<float>()});                            // the runtime, deduced from the result
```

`Q` weights the block temperature, `R = 0.02` the power itself: an LQR has no rate weight, where the MPC's `Rdu` is on the power change, so the two are tuned to comparable aggressiveness rather than identical cost. The filter's process noise is on the block alone; the register state is a past power, known exactly.

The law knows nothing of $`0 \le P \le P_{\max}`$. The tick clamps its request at the actuator and, through `feedback` / `commit`, tells the estimator what was actually applied; without that the filter would believe in heat that never went in and the estimate of $`w`$ would wind up during the heat-up.

Tick in code: `heater_period(controller, setpoint, temperature)`, the same function in the sketch and the SIL. The sketch flashes a `constinit` controller designed at compile time.

---

## What the SIL shows

The hotend with its true 5 s transport delay, a room 3 °C colder than the controller assumes, and a 5 W part-cooling fan from 240 s. Cold start to 210 °C.

| Heat-up overshoot | Last-minute error under the fan | Peak power |
| --- | --- | --- |
| 0.08 °C | 0.000 °C | 40 W |

On this plant the clamped LQG matches the MPC of `heater_ofmpc/`: only full power during the heat-up ever binds, and a clamp with the estimator told about it handles that one constraint as well as the QP. The MPC earns its keep where constraints bind during normal operation, or where there are several of them.

The SIL fails on an overshoot above 0.5 °C, on an offset above 0.1 °C under the fan, or on power above 40 W.

---

## Forking it

- Put the fitted $`K`$, $`\tau`$, $`\theta`$ and $`T_s`$ in the nameplate and the design takes the dead time in periods, `kDeadTime / kTs`. The model size and the runtime type follow.
- A trajectory to track rather than a setpoint to hold is the `feedback(y, x_ref, u_ref)` form, which needs the delay in the model (`InputDelay<D>`, the register form) rather than predicted.
