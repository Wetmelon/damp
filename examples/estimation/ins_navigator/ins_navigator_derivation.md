# 15-state INS navigator — design / deploy model

Teaching notes for [`ins_navigator_estimator.hpp`](ins_navigator_estimator.hpp) /
[`ins_navigator_sketch.cpp`](ins_navigator_sketch.cpp).
Same sensor densities as the host SIL [`../ins_eskf/`](../ins_eskf/)
(DiD parity). Design path: `design::ins_eskf_design` → `InsNavigator<float>`.

Canonical references: Solà (ESKF, 2017); Groves, *Principles of GNSS, Inertial, and Multisensor Integrated Navigation* (2nd ed.).

---

## Physical setup / nameplate

| Symbol | Code | Value | Meaning |
| ------ | ---- | ----- | ------- |
| $`\Delta t`$ | `kTs` | $`0.01`$ s | IMU period (100 Hz) |
| Frame | `kFrame` | ENU | East–North–Up local level |
| $`\sigma_g`$ | `kGyroNd` | $`0.003`$ rad/s/√Hz | Gyro noise density |
| $`\sigma_a`$ | `kAccelNd` | $`0.03`$ m/s²/√Hz | Accel noise density |
| $`\sigma_{bg}`$ | `kBgRw` | $`10^{-4}`$ rad/s^{3/2} | Gyro bias RW density |
| $`\sigma_{ba}`$ | `kBaRw` | $`0.001`$ m/s^{5/2} (≡ m/s³/√Hz) | Accel bias RW density |
| $`\sigma_p`$ | `kPosStd` | $`1.0`$ m | Default position aid 1-σ |
| $`\ell_b`$ | `kBaselineBody` | $`(1,0,0)`$ m | Dual-antenna baseline (body +x) |
| $`R_\psi`$ | sketch | $`(0.02)^2`$ | Heading aid variance (~1°) |

Design call:

```text
design::ins_eskf_design(kGyroNd, kAccelNd, kBgRw, kBaRw, kTs, kPosStd).as<float>()
```

Default $`P_0`$ 1-σ (library defaults, not overridden): attitude $`0.1`$ rad, velocity $`1`$ m/s, position $`10`$ m, $`b_g`$ $`0.01`$ rad/s, $`b_a`$ $`0.1`$ m/s².

---

## Nominal navigation state

`InsState` integrated by strapdown mechanization (`mechanize_step` inside `predict`):

```math
x =
\{
q,\;
v,\;
p,\;
b_g,\;
b_a
\}
```

| Quantity | Units | Meaning |
| -------- | ----- | ------- |
| $`q`$ | — | Body → nav quaternion |
| $`v`$ | m/s | Velocity in nav (ENU) |
| $`p`$ | m | Position in nav |
| $`b_g`$ | rad/s | Gyro bias (body) |
| $`b_a`$ | m/s² | Accel bias (body) |

Bias-corrected IMU and one-step mechanization (flat Earth, no Earth rate / transport rate):

```math
\begin{aligned}
\omega &= \omega_m - b_g
\\
a_b &= a_m - b_a
\\
a_n &= R(q)\,a_b + g_n
\\
v^+ &= v + a_n\,\Delta t
\\
p^+ &= p + v\,\Delta t + \tfrac12 a_n (\Delta t)^2
\\
q^+ &= q \otimes \Delta q(\omega\,\Delta t)
\end{aligned}
```

ENU gravity: $`g_n = (0,0,-g)`$ with $`g = 9.80665`$ (`kStandardGravity`).
At rest, body specific force $`a_b = -R(q)^\top g_n`$ (`specific_force_at_rest`); mock accel in the sketch uses that helper.

---

## Error state (15)

```math
\delta x =
\begin{bmatrix}
\delta\theta \\
\delta v \\
\delta p \\
\delta b_g \\
\delta b_a
\end{bmatrix}
\in \mathbb{R}^{15}
```

Index blocks (`ins_err`): $`\delta\theta`$ 0–2, $`\delta v`$ 3–5, $`\delta p`$ 6–8, $`\delta b_g`$ 9–11, $`\delta b_a`$ 12–14.

Right-multiplicative attitude error: $`q_{\mathrm{true}} = q \otimes \exp(\delta\theta)`$.

Continuous error model (Solà local ESKF; $`\omega = \omega_m - b_g`$, $`a_b = a_m - b_a`$, $`R = R(q)`$ body→nav):

```math
\begin{aligned}
\dot{\delta\theta}
&=
-[\omega]_\times\,\delta\theta
- \delta b_g
\\
\dot{\delta v}
&=
-R[a_b]_\times\,\delta\theta
- R\,\delta b_a
\\
\dot{\delta p}
&=
\delta v
\\
\dot{\delta b}_g
&=
0
\\
\dot{\delta b}_a
&=
0
\end{aligned}
```

Discretization: $`F_d = I + F_c\,\Delta t`$ (`ins_error_jacobian`). Leading blocks:

```math
\begin{aligned}
\delta\theta^+
&\approx
\bigl(I - [\omega]_\times\Delta t\bigr)\delta\theta
- \Delta t\,\delta b_g
\\
\delta v^+
&\approx
\delta v
- R[a_b]_\times\Delta t\,\delta\theta
- R\,\Delta t\,\delta b_a
\\
\delta p^+
&\approx
\delta p
+ \Delta t\,\delta v
\end{aligned}
```

Process noise input matrix $`G = I`$ (discrete-additive $`Q`$ from design).

---

## Process noise $`Q`$ (discrete)

Per axis, continuous–discrete integrated RW (`ins_eskf_design`):

```math
\begin{aligned}
Q_{\theta\theta}
&=
\sigma_g^2\,\Delta t + \sigma_{bg}^2\frac{\Delta t^3}{3}
\\
Q_{\theta,b_g}
&=
-\sigma_{bg}^2\frac{\Delta t^2}{2}
\\
Q_{b_g b_g}
&=
\sigma_{bg}^2\,\Delta t
\\
Q_{vv}
&=
\sigma_a^2\,\Delta t + \sigma_{ba}^2\frac{\Delta t^3}{3}
\\
Q_{v,b_a}
&=
-\sigma_{ba}^2\frac{\Delta t^2}{2}
\\
Q_{b_a b_a}
&=
\sigma_{ba}^2\,\Delta t
\\
Q_{pp}
&=
\sigma_a^2\frac{\Delta t^3}{3}
\\
Q_{pv}
&=
\sigma_a^2\frac{\Delta t^2}{2}
\end{aligned}
```

Accuracy note. These are closed-form per-axis integrals for white rate/accel
noise plus bias random walk (position from double-integrating white accel). They are
not a full van Loan $`Q_d`$ of the 15-state continuous $`F`$ (skew-$`\omega`$ /
$`R[a_b]_\times`$ coupling of process noise into other states is omitted). At typical
IMU rates this is the usual product approximation; large $`\|\omega\|\Delta t`$ or
aggressive specific-force dynamics need a denser $`Q`$ or smaller $`\Delta t`$.

---

## Measurement models (sparse aids)

Default position $`R`$ from design (isotropic):

```math
R_p = \sigma_p^2\, I_3
= I_3
\quad(\sigma_p = 1\,\mathrm{m})
```

Position (`update_position`): $`z = p_m - p`$, selector $`H`$ = $`I_3`$ on $`\delta p`$.

Heading (`update_heading`): dual-antenna / yaw aid. Forward direction from baseline $`\ell_b`$ mapped to nav; ENU heading

```math
\psi = \operatorname{atan2}(E, N)
```

Sketch variance $`R_\psi = (0.02)^2`$. Optional body baseline $`\ell_b = (1,0,0)`$ m.

Control period (sketch): `predict(imu, kTs)` every tick; if GPS mock true → `update_position`; if dual-antenna mock true → `update_heading`.

---

## Design → deploy

1. Bake $`Q,R,P_0`$ at compile time with `design::ins_eskf_design(...)`.
2. `static constinit InsNavigator<float> nav{kInsDesign, kFrame}`.
3. ISR / RTOS: `nav.predict` + sparse absolute aids; read `nav.state().p/v/q/b_g/b_a`.

No plant ODE is closed on target — the “plant” is the vehicle + IMU; the navigator is open-loop mechanization corrected by absolute measurements.

---

## Files

| File | Role |
| ---- | ---- |
| [`ins_navigator_derivation.md`](ins_navigator_derivation.md) | This derivation |
| [`ins_navigator_estimator.hpp`](ins_navigator_estimator.hpp) | Nameplate, design, `estimate_period` |
| [`ins_navigator_sketch.cpp`](ins_navigator_sketch.cpp) | Flashable sketch |
| [`ins_navigator_sil.cpp`](ins_navigator_sil.cpp) | Finite-tick host smoke |
| [`../ins_eskf/`](../ins_eskf/) | Host 3D SIL of the same design |
| `inc/damp/estimation/ins_eskf.hpp` | Design + `InsNavigator` |
| `inc/damp/estimation/ins_mechanization.hpp` | Strapdown step, frames, `InsState` |
