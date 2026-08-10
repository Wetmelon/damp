# PID / PI deploy sketch — plant notes

Teaching notes for [`pid_controller.hpp`](pid_controller.hpp).

This file is a firmware-shaped Design Is Deploy sketch: design a PI once,
`discretize` + `.as<float>()`, then call `controller.control(r, y)` in `loop()`.
There is no ODE plant inside the example — the process variable is whatever
`read_process_variable()` returns (mock `0` on host; ADC on target).

For a host SIL of a linear plant against the same PI API, use
`sim::simulate_lti_siso` (see `tests/test_simulate.cpp`) with a plant you choose.

---

## What is in the sketch (nameplate)

| Symbol | Code | Value | Meaning |
| ------ | ---- | ----- | ------- |
| $`T_s`$ | `Ts` | $`0.01`$ s | Control period (100 Hz) |
| $`K_p`$ | `design::pid` | $`4.0`$ | Proportional gain |
| $`K_i`$ | | $`8.0`$ | Integral gain |
| $`K_d`$ | | $`0.0`$ | Derivative off → pure PI |
| $`b`$ | default | $`1.0`$ | Setpoint weight on $`K_p`$ (library default) |
| $`u_{\min},u_{\max}`$ | default | unbounded | Actuator limits (optional 4th/5th args) |

Call site — signature is $`(K_p, K_i, K_d [, u_{\min}, u_{\max}, \ldots])`$, not
$`(K_p, K_i, K_d, N, b)`$. Passing $`0, 1`$ after $`K_d`$ clamps the command to
$`[0,1]`$ (easy footgun). Omit limits for an unbounded teaching PI; $`b`$ stays 1.

```text
design::pid(4.0, 8.0, 0.0).discretize(Ts).as<float>()
→ constinit PIController<float>
```

I/O contract:

- Reference $`r`$ — `reference` (demo starts at $`1.0`$).
- Measurement $`y`$ — process variable (same units as $`r`$).
- Command $`u =`$ `controller.control(r, y)` — actuator units (duty, volts, …).

On hardware: replace mocks with ADC/PWM; keep the same controller object and
period (timer ISR or `delay` approximating $`T_s`$).

---

## Controller law (what the plant “sees”)

Ideal continuous PI (teaching form; runtime is the library’s discrete
`PIController` after Tustin-style discretization at $`T_s`$):

```math
U(s)
=
\Bigl(
K_p + \frac{K_i}{s}
\Bigr)
E(s)
,\qquad
e = r - y
\quad\text{(with setpoint weight } b \text{ on the } K_p \text{ path when } b \ne 1\text{)}
```

With the sketch numbers ($`K_d = 0`$, $`b = 1`$):

```math
C(s) = 4 + \frac{8}{s}
=
\frac{4s + 8}{s}
```

There is no plant $`P(s)`$ in the `.cpp`. Closed-loop behavior is entirely
determined by whatever process sits between `write_actuator(u)` and
`read_process_variable()`.

---

## Suggested plant for host SIL (not in this file)

If you attach a first-order process for teaching SIL (common hobby plant):

```math
P(s) = \frac{K}{\tau s + 1}
```

then the continuous open loop is $`C(s)P(s)`$ and you can tune $`K_p,K_i`$ against
rise time / overshoot offline before flashing the same `PIController`. Prefer
library maps when a known plant model exists (`design::pi_pole_placement_first_order`,
etc.) — this sketch intentionally uses fixed `design::pid(...)` gains to show the
deploy path with minimal includes.

---

## Design → deploy steps

1. Design in `double` on the host (`design::pid` / placement maps).
2. `.discretize(Ts)` then `.as<float>()` into `static constinit PIController<float>`.
3. `setup()`: `controller.reset()`.
4. `loop()` / ISR: $`y \leftarrow`$ sensor, $`u \leftarrow control(r,y)`$, write actuator.
5. Do not pull `damp/control.hpp` on a small MCU for this demo — lean includes only.

---

## Files

| File | Role |
| ---- | ---- |
| [`pid_derivation.md`](pid_derivation.md) | This note |
| [`pid_controller.hpp`](pid_controller.hpp) | Nameplate + `design::pid` + `control_period` |
| [`pid_sketch.cpp`](pid_sketch.cpp) | Flashable setup/loop smoke |
| [`pid_sil.cpp`](pid_sil.cpp) | Host smoke ticks (mock $`y=0`$, no ODE plant) |

Related fuller DiD with plants: [`cart_pole/`](cart_pole/),
[`../pendulum/`](../pendulum/),
[`../multirate_rig/`](../multirate_rig/).
