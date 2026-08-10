# Reaction-wheel roll — plant model

Teaching notes for [`reaction_wheel_controller.hpp`](reaction_wheel_controller.hpp).

One linearized plant, many controller synthesis paths (cascade PI, ADRC/SMC/STSMC,
place / LQR / LQI, LQG / LQGI, offset-free MPC). Section 1 is firmware-shaped;
Section 2 only prints gains and runs one tick each — not closed-loop plant SIL.

---

## Physical setup

| Symbol | Code | Value | Meaning |
| ------ | ---- | ----- | ------- |
| $`\omega_0`$ | `w0` | $`1.0`$ rad/s | Open-loop unstable natural frequency |
| $`b_u`$ | `b_u` | $`1.0`$ rad/s² per N·m | Input gain: $`\dot\omega`$ per unit wheel torque |
| $`T_s`$ | `Ts` | $`0.001`$ s | Sample period (1 kHz control) |
| $`\|u\|_{\max}`$ | `u_max` | $`12`$ N·m | Wheel torque limit (cascade / MPC) |
| $`\|\omega^*\|_{\max}`$ | `w_cmd_max` | $`8`$ rad/s | Outer-loop rate command ceiling |

States:

```math
x =
\begin{bmatrix} \theta \\ \omega \end{bmatrix}
=
\begin{bmatrix}
\text{body roll [rad]} \\
\text{body roll rate [rad/s]}
\end{bmatrix}
```

Input $`u`$ — reaction-wheel torque [N·m] about the same body axis.

Output $`y = \theta`$ (roll). Rate $`\omega`$ is available for cascade / full-state
designs (gyro); LQG/LQGI estimate state from $`y`$ only.

Physical picture: body roll about upright with a single reaction wheel on that
axis. A wheel can only torque its spin axis — body-frame SISO by construction.
$`\omega_0`$ is a teaching stand-in for $`\sqrt{m g \ell / I}`$-style lean dynamics
(unstable upright); the code does not expand $`m,\ell,I`$ separately.

Sensing (not simulated here): treat $`\theta`$ as a fused roll estimate
(gyro integrated + accel tilt / complementary filter) and $`\omega`$ from the gyro.
Full 6-axis ESKF is out of scope.

---

## Continuous dynamics (linear plant)

As written in the file header and `plant` aggregate:

```math
\begin{aligned}
\dot\theta &= \omega \\
\dot\omega &= \omega_0^2\,\theta + b_u\, u
\end{aligned}
```

No viscous term on $`\omega`$ in this nameplate. In matrix form $`\dot x = A x + B u`$:

```math
A =
\begin{bmatrix}
0 & 1 \\
\omega_0^2 & 0
\end{bmatrix}
,\qquad
B =
\begin{bmatrix}
0 \\
b_u
\end{bmatrix}
```

With the nameplate numbers ($`\omega_0 = 1`$, $`b_u = 1`$):

```math
A =
\begin{bmatrix}
0 & 1 \\
1 & 0
\end{bmatrix}
,\qquad
B =
\begin{bmatrix}
0 \\
1
\end{bmatrix}
,\qquad
C =
\begin{bmatrix}
1 & 0
\end{bmatrix}
```

Open-loop eigenvalues of $`A`$ are $`\pm\omega_0`$ (one RHP pole) — upright is unstable.

Process / measurement noise channels for Kalman-based designs (identity, teaching
scale — not calibrated sensors):

```math
G = I_{2\times 2}
,\qquad
H = I_{1\times 1}
```

Discrete plant for LQ* / place / MPC:

```math
(A_d, B_d, \ldots) = \mathrm{ZOH}\bigl(A, B; T_s\bigr)
\quad\text{via}\quad
\texttt{discretize(plant, Ts, ZOH)}
```

---

## Transfer functions (SISO slices used by cascade)

Full linear plant:

```math
P(s)
=
\frac{\Theta(s)}{U(s)}
=
\frac{b_u}{s^2 - \omega_0^2}
=
\frac{1}{s^2 - 1}
```

Cascade successive-loop-closure approximations (when the inner rate loop is tight):

```math
\frac{\Omega(s)}{U(s)}
\approx
\frac{b_u}{s}
=
\frac{1}{(1/b_u)\, s}
\qquad
\text{(rate plant for inner PI)}
```

```math
\frac{\Theta(s)}{\Omega^*(s)}
\approx
\frac{1}{s}
\qquad
\text{(angle plant for outer PI when rate tracks } \omega^*\text{)}
```

Bandwidth separation in the demo: $`w_{\mathrm{rate}} = 25`$ rad/s,
$`w_{\mathrm{angle}} = 5`$ rad/s (~$`\times 5`$).

---

## Control design map (what Section 1 synthesizes)

All designs share the plant above. Weights for model-based linear laws:

| Weight | Value | Intent |
| ------ | ----- | ------ |
| $`Q`$ | $`\mathrm{diag}(40, 4)`$ | Angle heavily, rate lightly |
| $`R`$ | $`1`$ | Unit torque penalty |
| $`Q_{\mathrm{aug}}`$ | $`\mathrm{diag}(40, 4, 80)`$ | LQI/LQGI + integral on $`\xi`$ |
| $`Q_{\mathrm{kf}}`$ | $`\mathrm{diag}(10^{-4}, 10^{-3})`$ | Process noise (larger on $`\omega`$) |
| $`R_{\mathrm{kf}}`$ | $`10^{-3}`$ | Measurement noise on $`\theta`$ |

| Family | Library path | Notes |
| ------ | ------------ | ----- |
| Cascade PI | `pi_pole_placement_first_order` on $`b_u/s`$ then $`1/s`$ | Industry default on rate-capable axis |
| ADRC | `design::adrc<2>(wc=8, wo=40, b0=b_u)` | Treats $`\ddot\theta = f + b_0 u`$ |
| SMC / STSMC | `design::smc` / `design::stsmc` | Sliding on $`e,\dot e`$ (ST continuous $`u`$) |
| Place | `design::place(A,B,{-8,-12}, Ts)` | Discrete poles from continuous intent |
| LQR | `design::discrete_lqr(Ad,Bd,Q,R)` | Infinite-horizon discrete |
| LQI | `design::discrete_lqi(sysd, Q_aug, R)` | Integrator on $`r-y`$ |
| LQG / LQGI | `design::discrete_lqg` / `discrete_lqgi` | Separation + Kalman from $`y`$ |
| OF-MPC | `design::mpc<10,4>` | Synthesis only; no constraint SIL |

Cascade tick: $`\omega^* =`$ outer angle PI, $`u =`$ inner rate PI
(`cascade_period`). Full-state laws use $`u = -Kx`$ (plus integral / KF / MPC as appropriate).

---

## Design plant vs SIL truth

There is no nonlinear closed-loop plant SIL. Section 2 / `reaction_wheel_sil.cpp`
only prints gains and runs one tick of each law on a fixed $`x_0`$ — intentional
synthesis gallery, not gate-and-validate against an ODE. Cascade PI designs use
the successive-loop-closure approximations $`b_u/s`$ and $`1/s`$ (not the full
unstable $`b_u/(s^2-\omega_0^2)`$ plant); full-state laws use the linear $`plant`
aggregate above.

## Files

| File | Role |
| ---- | ---- |
| [`reaction_wheel_derivation.md`](reaction_wheel_derivation.md) | This plant + design map |
| [`reaction_wheel_controller.hpp`](reaction_wheel_controller.hpp) | Nameplate plant + all designs + `cascade_period` |
| [`reaction_wheel_sketch.cpp`](reaction_wheel_sketch.cpp) | Cascade float deploy smoke |
| [`reaction_wheel_sil.cpp`](reaction_wheel_sil.cpp) | Host: print gains + one tick each law |

No closed-loop plant plot in this example — gains and single-tick $`u`$ only.
