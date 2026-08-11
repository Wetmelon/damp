# Workflow end-to-end — nonlinear plant and linearization

Teaching notes for [`workflow_end_to_end_sil.cpp`](workflow_end_to_end_sil.cpp).

Design Is Deploy fence: linearize a nonlinear plant → discrete LQGI + PR →
frequency analysis → nonlinear SIL calling the same float runtime.

---

## Nameplate / rates

| Symbol | Code | Value | Meaning |
| ------ | ---- | ----- | ------- |
| $`T_s`$ | `kTs` | $`0.001`$ s | 1 kHz sample |
| $`x_{\mathrm{op}}`$ | `kXop` | $`[0,\,0]^\top`$ | Linearization point |
| $`u_{\mathrm{op}}`$ | `kUop` | $`0`$ | Input at equilibrium |
| $`r`$ | `reference` | $`0.5`$ | Closed-loop setpoint (SIL) |
| $`x_0`$ | | $`[0.25,\,0]^\top`$ | SIL initial state |
| horizon | | $`[0,\,2]`$ s | Simulation |

---

## Nonlinear plant

State $`x = [x_1,\, x_2]^\top`$, scalar input $`u`$, output $`y = x_1`$:

```math
\begin{aligned}
\dot x_1
&=
x_2
\\
\dot x_2
&=
-0.8\, x_2 - 2.0\,\sin(x_1) + 1.5\, u
\end{aligned}
```

(`plant_nonlinear` / `plant_output` in the example.) First-order mechanical
analogy: position $`x_1`$, rate $`x_2`$, viscous damping, nonlinear “gravity-like”
$`\sin`$ restoring, force gain 1.5.

---

## Design plant vs SIL truth

The same `plant_nonlinear` is used for finite-difference linearization and for
the nonlinear closed-loop SIL (no second plant). Design is at
$`(x_{\mathrm{op}}, u_{\mathrm{op}}) = (0, 0)`$ while SIL tracks $`r = 0.5`$ from
$`x_0 = [0.25, 0]^\top`$ — intentional op-point ≠ setpoint so the LQGI integral
(and PR) must supply the equilibrium bias
$`u_{\mathrm{eq}} = (2/1.5)\sin(0.5)`$ on the nonlinear plant.

## Linearization at the origin

At $`(x,u)=(0,0)`$: $`\sin(x_1)\approx x_1`$, so

```math
\begin{aligned}
\dot x_1
&=
x_2
\\
\dot x_2
&=
-2.0\, x_1 - 0.8\, x_2 + 1.5\, u
\end{aligned}
```

In matrix form $`\dot x = A x + B u`$, $`y = C x`$:

```math
A
=
\begin{bmatrix}
0 & 1 \\
-2 & -0.8
\end{bmatrix}
,\qquad
B
=
\begin{bmatrix}
0 \\
1.5
\end{bmatrix}
,\qquad
C
=
\begin{bmatrix}
1 & 0
\end{bmatrix}
,\qquad
D = 0
```

Code uses `design::linearize<2,1,1>(f,h,x_op,u_op)` (central finite differences); analytic
Jacobians match the above at the origin. Noise embedding for Kalman design:

```math
G = I_2
,\qquad
H = 1
```

Discrete plant: `discretize(sys_c, Ts, ZOH)`.

---

## Design artifacts (deploy path)

1. PID seed (optional side path): `pid_from_performance_spec` with
   $`t_s=0.20`$, overshoot 10%, PI, $`T_s`$.
2. LQGI: `design::lqgi_bundle(sys_d, Q_aug, R, Q_kf, R_kf)` with

```math
Q_{\mathrm{aug}}
=
\operatorname{diag}(20,\, 2,\, 80)
,\qquad
R = 0.25
```

```math
Q_{\mathrm{kf}}
=
\begin{bmatrix}
10^{-3} & 0 \\
0 & 10^{-2}
\end{bmatrix}
,\qquad
R_{\mathrm{kf}} = 5\cdot 10^{-3}
```

Success gated before runtime use. Float runtime from the bundle.

3. PR: `design::pr(0, 10, 2\pi, 6, Ts)` — resonant path on the tracking error;
   float `PRController`.

SIL control law (double plant, float controllers):

```math
u
=
u_{\mathrm{LQGI}}(r,\, y) + u_{\mathrm{PR}}(r - y)
```

Integrator: RK4 fixed step $`T_s`$. Host also computes discrete Bode / Nyquist on
the open-loop linearized plant.

---

## Files

| File | Role |
| ---- | ---- |
| [`workflow_end_to_end_derivation.md`](workflow_end_to_end_derivation.md) | This plant + design map |
| [`workflow_end_to_end_controller.hpp`](workflow_end_to_end_controller.hpp) | Nameplate plant, rates, ops |
| [`workflow_end_to_end_sketch.cpp`](workflow_end_to_end_sketch.cpp) | Float PR smoke (LQGI baked on host) |
| [`workflow_end_to_end_sil.cpp`](workflow_end_to_end_sil.cpp) | Linearize → LQGI+PR → nonlinear SIL |

No dedicated plot file in this example (console summary only).
