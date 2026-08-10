# IMU attitude open-loop — gyro integrate teaching model

Teaching notes for [`imu_pose_estimator.hpp`](imu_pose_estimator.hpp) /
[`imu_pose_sil.cpp`](imu_pose_sil.cpp).
Host-only. Orientation dead-reckoning from body rates — no filter, no position,
no accel/mag aiding. Drift under constant gyro bias is why MARG ESKF / dual-antenna
heading exist ([`../eskf/`](../eskf/)).

---

## Physical setup

| Symbol | Code | Value | Meaning |
| ------ | ---- | ----- | ------- |
| $`\Delta t`$ | `kTsSil` | $`1/60`$ s | Integration step |
| $`t_{\mathrm{end}}`$ | `t_end` | $`10`$ s | Scenario length |
| $`b_g`$ | `kBgZ` / `gyro_bias()` | $`(0, 0, 0.05)`$ rad/s | Constant yaw bias on the raw gyro |

True body rate schedule $`\omega(t)`$ [rad/s] (`omega_true`):

```math
\omega(t) =
\begin{cases}
(0,\; 0,\; 0.6) & t < 3~\mathrm{s} \quad\text{(yaw about body }z\text{)} \\
(0,\; 0.5,\; 0) & 3 \le t < 6~\mathrm{s} \quad\text{(pitch about body }y\text{)} \\
(0.4,\; 0,\; 0) & 6 \le t < 9~\mathrm{s} \quad\text{(roll about body }x\text{)} \\
(0,\; 0,\; 0) & t \ge 9~\mathrm{s}
\end{cases}
```

Initial orientation: $`q = I`$ for both ideal and open-loop paths.

---

## Frames

| Axes drawn | Role |
| ---------- | ---- |
| World / local-level (fixed, longdash) | Reference frame the body is drawn in (nav-like fixed axes; not an estimate) |
| Ideal body (solid) | Integrate true $`\omega(t)`$ |
| Open-loop body (dotted) | Integrate $`\omega(t) + b_g`$ as a raw gyro would |

Body axes are drawn at the origin (attitude only; no translation).

---

## Kinematic model

Unit quaternion $`q`$ maps body → world. Discrete step (`estimate_period` →
`Quaternion::integrate_body_rates`), first-order right-multiplicative:

```math
\begin{aligned}
\Delta q
&=
\bigl(1,\; \tfrac12\omega_x\Delta t,\; \tfrac12\omega_y\Delta t,\; \tfrac12\omega_z\Delta t\bigr)
\\
q^+
&=
\operatorname{normalize}\bigl(q \otimes \Delta q\bigr)
\end{aligned}
```

Ideal path: integrate true $`\omega(t)`$.

Open-loop biased path (not an estimator — pure open-loop sensor integration):

```math
q_{\mathrm{ol}}^+
=
q_{\mathrm{ol}}.\mathrm{integrate\_body\_rates}\bigl(\omega(t) + b_g,\, \Delta t\bigr)
```

with $`b_g = (0,0,0.05)`$ rad/s constant. The dotted triad peels away in yaw from the integrated bias.

Modeling assumptions:

- Perfect knowledge of true $`\omega`$ for the ideal path; no quantization or white noise.
- Bias is constant and only on the open-loop path.
- No accelerometer or magnetometer correction; no ESKF injection.
- Single rigid body at a fixed point (pose without translation).

---

## Why this is not “pose under bias”

The dotted axes are not a filtered estimate of orientation given biased sensors.
They are the trajectory you get if you integrate a bad gyro and never correct it.
Bounded attitude under real IMUs needs vector aids (accel / mag) or absolute
orientation / heading measurements — the ESKF path in `design::eskf_marg` /
`InsNavigator`.

---

## Outputs

Host plot: `plots/estimation/imu_pose_3d.html`.

---

## Files

| File | Role |
| ---- | ---- |
| [`imu_pose_derivation.md`](imu_pose_derivation.md) | This derivation |
| [`imu_pose_estimator.hpp`](imu_pose_estimator.hpp) | Nameplate, schedule, tick |
| [`imu_pose_sketch.cpp`](imu_pose_sketch.cpp) | Thin smoke |
| [`imu_pose_sil.cpp`](imu_pose_sil.cpp) | Host animation |
| `inc/damp/math/geometry.hpp` | `Quaternion::integrate_body_rates` |
