# ADRC vs PID — turning a flywheel to an angle

Read this before forking the sketch. Numbers here must match
[`adrc_controller.hpp`](adrc_controller.hpp) and [`adrc_sil.cpp`](adrc_sil.cpp).

**What you are looking at.** A motor turns a flywheel. You command an
angle. Two controllers try the same job, same speed of response, same
torque limit:

- **PID** (proportional–integral–derivative) — the usual industrial
  controller. Tuned carefully for the flywheel as written on the
  datasheet.
- **ADRC** (active disturbance rejection control) — the thing this
  folder is meant to teach. It barely models the machine. It watches
  the angle and estimates “everything else,” then cancels that.

The host sim (`adrc_sil.cpp`) runs **two copies** of the same flywheel,
one per controller, and plots them together. The flashable sketch only
ships ADRC.

---

## The idea in one paragraph

PID is told: “this is a 0.05 kg·m² inertia; make the angle loop as
fast as 8 rad/s.” That tune is good **until the real machine is not
that datasheet** — extra friction, someone hangs a weight on the
shaft, the load torque steps. PID’s gains are then the wrong size.

ADRC is told almost nothing: “torque mostly makes angular
acceleration, and the scale is about $`1/J`$.” Friction, load, and a
wrong inertia are lumped into one leftover $`f`$. A small observer
(the **ESO**, extended state observer) estimates $`f`$ from the
measured angle and subtracts it. You do not need a speed sensor.

---

## The flywheel (what the sim integrates)

Torque in, minus drag, minus a constant-ish rubbing torque, minus an
external load, equals inertia times angular acceleration:

```math
J(t)\,\ddot\theta = u - B\dot\theta - \tau_c\,\mathrm{soft\_sign}(\dot\theta,\omega_c) - \tau_L(t)
```

| What | Code | Value | Plain meaning |
| ---- | ---- | ----- | ------------- |
| $`J`$ (design) | `J_nom` | 0.05 kg·m² | Inertia you designed for (datasheet) |
| $`B`$ | `B_visc` | 0.03 N·m·s/rad | Speed-proportional drag (oil / bearings) |
| $`\tau_c`$ | `tau_c` | 0.12 N·m | Rubbing / Coulomb friction (opposes motion) |
| $`\omega_c`$ | `omega_c` | 0.04 rad/s | How sharply that rubbing turns on through zero speed |
| $`u_{\max}`$ | `u_max` | 2.5 N·m | Motor torque limit (both controllers) |
| $`T_s`$ | `Ts` | 0.001 s | Control tick (1000 times per second) |
| loop $`\omega_c`$ | `wc` | 8 rad/s | How snappy you asked the *angle* loop to be |
| observer $`\omega_o`$ | `wo` | 40 rad/s | How fast ADRC’s leftover-estimator runs (5× the loop) |
| $`b_0`$ | `b0` | $`1/J = 20`$ | “How much acceleration per newton-metre,” the only plant number ADRC gets |

`soft_sign` is a smooth sign function (library helper). It is just
“rubbing that does not chatter when the wheel stops.”

**What happens in time**

| Time | What we do to the machine |
| ---- | ------------------------- |
| 0 s | Ask for 1 radian |
| 1.2 s | Hang a 0.45 N·m load on the shaft |
| 2.4 s | Triple the inertia (payload on the wheel) |
| 3.6 s | Ask to go back to 0 |

So $`J(t) = J_{\mathrm{nom}}`$ until 2.4 s, then $`3J_{\mathrm{nom}}`$.
$`\tau_L = 0`$ until 1.2 s, then 0.45 N·m.

---

## The two tunings (same requested speed)

Both are aimed at the same loop speed `wc = 8` rad/s. That is the
fair part.

### ADRC (what you would flash)

`design::adrc<2>(wc, wo, b0)`:

- $`K_p = \omega_c^2`$, $`K_d = 2\omega_c`$ — a PD on the *estimated*
  angle and speed
- $`b_0 = 1/J_{\mathrm{nom}}`$ — “torque to acceleration”
- each tick:

```math
u = \bigl(K_p(r - \hat\theta) - K_d\hat\omega - \hat f\bigr)/b_0
```

$`\hat\theta`$, $`\hat\omega`$, $`\hat f`$ come from the observer, not
from extra sensors. If the motor saturates, call `back_calculate` so
the observer is told the torque that actually happened, not the
torque the law wanted.

Tick in code: `control_period_adrc`.

### PID (honest datasheet tune)

This is **not** a deliberately bad PID. It is
`design::pid_pole_placement_double_integrator(J, wc, 1, wc/8)`:
place a critically-damped PD pair at $`-\omega_c`$ on
$`G(s)=1/(J s^2)`$, plus a slower integrator pole at $`-\omega_c/8`$
so a constant load does not leave a standing error.

```math
(s^2 + 2\omega_c s + \omega_c^2)\,(s + \omega_c/8)
```

- $`K_p`$ — push toward the angle error
- $`K_d`$ — damp with speed (taken from the *measurement*, not the
  setpoint, so a step in $`r`$ does not spike torque)
- $`K_i`$ — slowly wind until a constant load is cancelled
- $`T_f = 1/(10\omega_c)`$ — a little filter on the D term

Then `.discretize(Ts)` and `control_period_pid`. Same torque limit and
the same “tell the integrator what torque actually went out” hook.

### After the payload

When $`J`$ becomes $`3J`$, every PID gain is three times too big
for the new plant. ADRC still only needs $`b_0`$ in the right
ballpark; the observer folds the mismatch into $`\hat f`$ and
cancels it.

---

## Files

| File | What it is |
| ---- | ---------- |
| [`adrc_controller.hpp`](adrc_controller.hpp) | Datasheet numbers, both designs, the ticks |
| [`adrc_sketch.cpp`](adrc_sketch.cpp) | Tiny program you could copy onto a board (ADRC only) |
| [`adrc_sil.cpp`](adrc_sil.cpp) | Host sim: two flywheels + plot |
| [`adrc_derivation.md`](adrc_derivation.md) | This note |

Plot after you run the sim: `examples/plots/control/adrc_vs_pid.html`
(angle, torque, and ADRC’s leftover estimate $`\hat f`$).
