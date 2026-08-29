# PMSM FOC current loop (PI vs I-P)

Notes for [`foc_sil.cpp`](foc_sil.cpp) / [`foc_controller.hpp`](foc_controller.hpp).
Implementation: `damp/motor/foc.hpp` (`FOController`).

---

## Field-oriented control

FOC regulates stator current in the rotor $`dq`$ frame instead of in $`abc`$ or
$`\alpha\beta`$. With the electrical angle $`\theta_e`$ from the rotor:

1. Measure $`i_{abc}`$ → Clarke → Park($`\theta_e`$) → $`(i_d, i_q)`$.
2. Current regulators produce $`(v_d^*, v_q^*)`$.
3. Inverse Park → $`v_{\alpha\beta}`$ → PWM (SVPWM in a full drive; this SIL
   applies $`v_{dq}`$ directly to a dq plant).

For a surface PM machine ($`L_d = L_q`$), electromagnetic torque is

```math
T_e = \tfrac{3}{2}\,p\,\lambda\, i_q
```

so $`i_q`$ is the torque channel and $`i_d`$ is held at $`0`$ (MTPA). Cross-coupling
and back-EMF appear as known terms in the $`dq`$ voltage equations; FOC cancels
them with feedforward so each axis looks like a first-order $`R`$-$`L`$ plant to
its PI.

---

## Plant

Turnigy D5065 270KV surface PMSM. $`\omega_e`$ held fixed so only the current
loop is exercised.

| Symbol | Code | Value | Meaning |
| ------ | ---- | ----- | ------- |
| $`R_s`$ | `Rs` | $`0.039`$ Ω | Phase resistance |
| $`L_s`$ | `Ls` | $`16`$ µH | $`L_d = L_q`$ |
| $`K_t`$ | `Kt_spec` | $`0.031`$ Nm/A | Torque constant (amplitude) |
| $`p`$ | `pole_pairs` | $`7`$ | Pole pairs |
| $`\lambda`$ | `lambda_pm` | from $`K_t`$ | PM flux [Wb] |
| $`V_{\mathrm{dc}}`$ | `Vdc` | $`24`$ V | DC bus |
| $`\omega_e`$ | `omega_e` | $`600`$ rad/s | Fixed electrical speed |
| $`f_{\mathrm{sw}}`$ | `fsw` | $`20`$ kHz | Sample rate |
| $`T_s`$ | `Ts` | $`1/f_{\mathrm{sw}}`$ | Control period |
| $`\omega_{\mathrm{bw}}`$ | `bw` | $`2\pi f_{\mathrm{sw}}/13`$ | Current-loop bandwidth |

### Coupled $`dq`$ equations

```math
\begin{aligned}
L_d\,\dot i_d &= v_d - R_s\, i_d + \omega_e L_q\, i_q \\
L_q\,\dot i_q &= v_q - R_s\, i_q - \omega_e L_d\, i_d - \omega_e \lambda
\end{aligned}
```

### State-space form

Treat currents as the state, applied $`dq`$ voltages as the input, and measured
currents as the output:

| Symbol | Definition | Size |
| ------ | ---------- | ---- |
| $`x`$ | $`[i_d,\ i_q]^\top`$ | $`2\times 1`$ |
| $`u`$ | $`[v_d,\ v_q]^\top`$ | $`2\times 1`$ |
| $`y`$ | $`x`$ (currents measured) | $`2\times 1`$ |
| $`e`$ | affine back-EMF term $`[0,\ -\omega_e\lambda/L_q]^\top`$ | $`2\times 1`$ |

```math
\dot x = A x + B u + e,\qquad y = C x + D u
```

with

```math
A =
\begin{bmatrix}
-R_s/L_d & \omega_e L_q/L_d \\
-\omega_e L_d/L_q & -R_s/L_q
\end{bmatrix},
\quad
B =
\begin{bmatrix}
1/L_d & 0 \\
0 & 1/L_q
\end{bmatrix},
\quad
C = I_2,
\quad
D = 0
```

For this surface machine ($`L_d = L_q = L_s`$)

```math
A =
\begin{bmatrix}
-R_s/L_s & \omega_e \\
-\omega_e & -R_s/L_s
\end{bmatrix},
\quad
B = \frac{1}{L_s} I_2,
\quad
e =
\begin{bmatrix}
0 \\
-\omega_e\lambda/L_s
\end{bmatrix}.
```

SIL integrates $`\dot x`$ with forward Euler (`substeps = 100` per $`T_s`$). No
mechanical states — $`\omega_e`$ is an exogenous constant (frozen $`A,e`$ over the
run).

### Decoupled axis used for tuning

After ideal decoupling / back-EMF cancellation, each axis is an independent
SISO plant. For one axis, redefine

| Symbol | Definition |
| ------ | ---------- |
| $`x`$ | current $`i`$ ($`i_d`$ or $`i_q`$) |
| $`u`$ | regulator voltage $`v`$ (PI output for that axis) |
| $`y`$ | $`x`$ |

```math
\dot x = A x + B u,\qquad y = C x + D u
```

```math
A = -\frac{R}{L},\qquad B = \frac{1}{L},\qquad C = 1,\qquad D = 0
```

Transfer function $`G(s) = C(sI-A)^{-1}B + D = 1/(L s + R)`$. That is the plant
passed to `pi_pole_placement_first_order` / `FOController::tune`.

---

## Current controller (`FOController::current_controller`)

### Decoupling feedforward

```math
\begin{aligned}
v_d^{\mathrm{ff}} &= -\omega_e L_q\, i_q^* \\
v_q^{\mathrm{ff}} &= \omega_e (L_d\, i_d^* + \lambda)
\end{aligned}
```

(`decoupling_feedforward`, $`R = 0`$, on the **references**.) Applied voltages are

```math
u =
\begin{bmatrix} v_d^{\mathrm{ff}} \\ v_q^{\mathrm{ff}} \end{bmatrix}
+
\begin{bmatrix} u_d \\ u_q \end{bmatrix}
```

where $`u_d,u_q`$ are the PI outputs.

### PI with setpoint weight $`b`$

On each axis, with reference $`r`$ and measurement $`y`$:

```math
\begin{aligned}
u &= K_p\bigl(b\, r - y\bigr) + I \\
\dot I &= K_i (r - y)
\end{aligned}
```

| $`b`$ | Structure | Proportional term |
| ----- | --------- | ----------------- |
| $`1`$ | PI | on error $`r - y`$ |
| $`0`$ | I-P | on measurement $`-y`$ only |

Optional plant-inversion FF ($`R i + L \Delta i / T_s`$) exists in the header but
is off here.

### Voltage circle and anti-windup

The linear SVPWM limit on the $`dq`$ voltage vector is

```math
V_{\max} = \frac{V_{\mathrm{dc}}}{\sqrt{3}}.
```

In code that is `voltage_circle_radius(Vdc)`. If $`\|v_{dq}^*\| > V_{\max}`$, the
command is scaled onto the circle and both PI integrators are back-calculated
so they do not wind up against the rail.

---

## Gain placement

`FOController::tune` places each axis on $`G(s) = 1/(L s + R)`$ with bandwidth
$`\omega_{\mathrm{bw}}`$ and $`\zeta = 1`$:

```math
K_p = 2\zeta\omega_{\mathrm{bw}} L - R,\qquad
K_i = L\,\omega_{\mathrm{bw}}^2,\qquad
K_{bc} = K_p
```

Same poles for PI and I-P; only $`b`$ differs.
$`\omega_{\mathrm{bw}} \approx 2\pi f_{\mathrm{sw}}/13`$ keeps the continuous design
valid at the sample rate.

Identical poles ⇒ identical linear disturbance rejection. The step response
differs: PI applies a full $`K_p \Delta r`$ kick; I-P ramps only through the
integral.

---

## Scenario in `foc_sil`

1. $`i_d^* = 0`$, $`i_q^*`$ steps $`0 \to 8`$ A at $`t = 0.5`$ ms.
2. Unmodeled $`v_q`$ disturbance $`-1`$ V at $`t = 3`$ ms.
3. Two runs: $`b = 1`$ (PI) and $`b = 0`$ (I-P).

Plots: `examples/plots/motor/foc_current_loop.html`,
`examples/plots/motor/foc_id_iq.html`.

| File | Role |
| ---- | ---- |
| [`foc_controller.hpp`](foc_controller.hpp) | Nameplate / rates |
| [`foc_sil.cpp`](foc_sil.cpp) | PI vs I-P SIL |
| [`foc_sketch.cpp`](foc_sketch.cpp) | Minimal `current_controller` tick |
