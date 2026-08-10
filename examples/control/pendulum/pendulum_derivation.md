# Simple pendulum (upright LQR) — plant model

Teaching notes for [`pendulum_controller.hpp`](pendulum_controller.hpp) /
[`pendulum_sil.cpp`](pendulum_sil.cpp).
Read this before forking the sketch. Code numbers must match the equations here.

Design Is Deploy fence: the controller header linearizes about upright, designs
discrete LQR, and exposes `control_period`. SIL integrates the nonlinear plant
with the same `controller.control(x)` path.

Related full design-is-deploy path (4-state cart–pole): [`../cart_pole/`](../cart_pole/)
([`cart_pole_derivation.md`](../cart_pole/cart_pole_derivation.md)).

---

## Physical setup

| Symbol | Code | Value | Meaning |
| ------ | ---- | ----- | ------- |
| $`g`$ | `g` | $`9.81`$ m/s² | Gravity |
| $`\ell`$ | `L` | $`1.0`$ m | Pivot → bob length |
| $`m`$ | `m` | $`1.0`$ kg | Point mass at tip |
| $`b`$ | `b_damp` | $`0.1`$ N·m·s/rad | Viscous damping at the pivot |
| $`T_s`$ | `Ts` | $`0.01`$ s | Sample period (100 Hz control) |

States (column vector $`x`$):

```math
x =
\begin{bmatrix} \theta \\ \dot\theta \end{bmatrix}
=
\begin{bmatrix}
\text{angle from upright [rad]} \\
\text{angular rate [rad/s]}
\end{bmatrix}
```

Input $`u`$ — torque about the pivot [N·m].

Angle convention: $`\theta = 0`$ is upright (unstable equilibrium). Positive
$`\theta`$ is a small tip; gravity pulls the bob *away* from upright for small $`\theta`$.
This is the inverted (upright) convention, not hanging-down $`\theta = \pi`$.

Modeling assumptions:

- Point mass at distance $`\ell`$; no separate rod inertia about CoM.
- Single viscous term $`b\,\dot\theta`$ (pivot damping only).
- Actuator applies torque $`u`$ directly (no motor dynamics).

---

## Nonlinear dynamics (SIL plant)

Used in `pendulum_sil.cpp` as `pendulum_f`. From Newton–Euler about
the pivot (or Lagrange with generalized coordinate $`\theta`$):

```math
m \ell^2 \ddot\theta = m g \ell \sin\theta - b\,\dot\theta + u
```

```math
\ddot\theta
=
\frac{g}{\ell}\sin\theta
-
\frac{b}{m\ell^2}\,\dot\theta
+
\frac{1}{m\ell^2}\,u
```

First-order form $`\dot x = f(x,u)`$ with $`x = [\theta,\, \dot\theta]^\top`$:

```math
\begin{aligned}
\dot x_1 &= x_2 \\
\dot x_2 &= \frac{g}{\ell}\sin x_1 - \frac{b}{m\ell^2}\, x_2 + \frac{1}{m\ell^2}\,u
\end{aligned}
```

(Code: `theta_ddot = (g/L)*sin(theta) - (b_damp/(m*L*L))*theta_dot + (1/(m*L*L))*torque`.)

SIL integrates with RK4 at $`\Delta t = 0.001`$ s; control is the discrete law at
$`T_s`$ (via `simulate_state_feedback`).

Initial condition in the demo: $`x_0 = [0.5236,\; 0]^\top`$ (~$`30^\circ`$ tip).

---

## Linearization at upright ($`\theta = 0`$)

At the origin $`x = 0`$, $`u = 0`$: $`\sin\theta \approx \theta`$. Jacobian matches
`linearize_upright()`:

```math
\begin{aligned}
\ddot\theta
&=
\frac{g}{\ell}\,\theta
-
\frac{b}{m\ell^2}\,\dot\theta
+
\frac{1}{m\ell^2}\,u
\end{aligned}
```

In matrix form $`\dot x = A x + B u`$:

```math
A =
\begin{bmatrix}
0 & 1 \\
g/\ell & -b/(m\ell^2)
\end{bmatrix}
,\qquad
B =
\begin{bmatrix}
0 \\
1/(m\ell^2)
\end{bmatrix}
```

With the nameplate numbers:

```math
A =
\begin{bmatrix}
0 & 1 \\
9.81 & -0.1
\end{bmatrix}
,\qquad
B =
\begin{bmatrix}
0 \\
1
\end{bmatrix}
```

```math
\begin{aligned}
a_{21}
&=
g/\ell = 9.81
\\
a_{22}
&=
-b/(m\ell^2) = -0.1
\\
b_2
&=
1/(m\ell^2) = 1
\end{aligned}
```

Output for design (angle only): $`C = [1\ 0]`$, $`y = \theta`$. LQR uses full state
feedback, not output feedback.

Open-loop $`A`$ has a positive eigenvalue $`\approx +\sqrt{g/\ell}`$ — upright is unstable.

---

## Transfer function (linear SISO)

From the linearized ODE (Laplace, zero IC):

```math
P(s)
=
\frac{\Theta(s)}{U(s)}
=
\frac{1/(m\ell^2)}{s^2 + \bigl(b/(m\ell^2)\bigr)s - g/\ell}
=
\frac{1}{s^2 + 0.1\, s - 9.81}
```

(Note the sign of the constant term: $`-g/\ell`$ for the upright linearization.)

---

## Control design (what the controller header does)

1. Continuous plant $`(A,B)`$ as above.
2. Discrete LQR via `design::discrete_lqr_from_continuous(A, B, Q, R, Ts)`.
3. Cost weights:

| State / input | Weight | Intent |
| ------------- | ------ | ------ |
| $`\theta`$ | $`Q_{11} = 10`$ | Angle regulation |
| $`\dot\theta`$ | $`Q_{22} = 10`$ | Rate damping |
| $`u`$ | $`R = 1`$ | Control effort |

($`Q = 10\,I_2`$.)

4. Law: $`u = -K x`$ (`LQR::control` / `control_period`).
5. Host SIL keeps `double`; flash path: `lqr_d.as<float>()` into `LQR<2,1,float>`.

Host redesign in SIL: same plant, $`Q = 50\,I`$ — stiffer recovery; phase
portrait compares $`Q=10`$ vs $`Q=50`$. Also prints a redesign about
$`\theta_{\mathrm{op}} = 0.5`$ rad (not the flash path).

Tuning tip: raise $`Q_{11}`$ for tighter upright recovery; raise $`R`$ if torque saturates.

---

## Files

| File | Role |
| ---- | ---- |
| [`pendulum_derivation.md`](pendulum_derivation.md) | This derivation |
| [`pendulum_controller.hpp`](pendulum_controller.hpp) | Nameplate + design + tick |
| [`pendulum_sketch.cpp`](pendulum_sketch.cpp) | Float deploy smoke |
| [`pendulum_sil.cpp`](pendulum_sil.cpp) | Nonlinear SIL + plots + host redesign |
| [`../cart_pole/`](../cart_pole/) | 4-state cart–pole design-is-deploy path |

Host plots: `examples/plots/control/pendulum_sim.html`,
`pendulum_sim_high_q.html`, `pendulum_phase.html`.
