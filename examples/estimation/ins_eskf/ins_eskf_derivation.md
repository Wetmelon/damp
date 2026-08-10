# InsNavigator 3D SIL — plant / filter model

Teaching notes for [`ins_eskf_estimator.hpp`](ins_eskf_estimator.hpp) /
[`ins_eskf_sil.cpp`](ins_eskf_sil.cpp).
Host-only animation of the same design densities as
[`../ins_navigator/`](../ins_navigator/)
(green truth / orange free-run / cyan aided).

---

## Physical setup (truth schedule)

| Symbol | Code | Value | Meaning |
| ------ | ---- | ----- | ------- |
| Frame | `kFrame` | ENU | Local-level axes |
| $`\Delta t`$ | `kTs` / `dt` | $`0.01`$ s | IMU period (matches deploy) |
| $`t_{\mathrm{end}}`$ | `t_end` | $`18`$ s | Scenario length |
| Aid rate | `aid_every` | every 50 samples | Absolute aids at 2 Hz |
| $`b_a^{\mathrm{sensor}}`$ | `sensor_ba` | $`(0.05, 0, 0)`$ m/s² | Unmodeled body-x accel bias on IMU stream |

Truth kinematics (nav acceleration and body rate commands):

```math
a_n(t) =
\begin{cases}
(1.2,\; 0,\; 0) & t < 3~\mathrm{s} \\
(0,\; 0,\; 0) & \text{otherwise}
\end{cases}
,\qquad
\omega_b(t) =
\begin{cases}
(0,\; 0,\; 0.35) & 4 \le t < 7~\mathrm{s} \\
(0,\; 0,\; 0) & \text{otherwise}
\end{cases}
```

($`a_n`$ in m/s² East; $`\omega_b`$ yaw rate in rad/s.)

Initial state: $`q = I`$, $`v = 0`$, $`p = 0`$, biases zero on the truth and filter nominal until estimated.

---

## Synthetic IMU (truth-based)

Ideal body specific force from known nav acceleration and orientation:

```math
a_b = R(q)^\top \bigl(a_n - g_n\bigr)
,\qquad
\omega_m = \omega_b
```

with $`g_n = \mathrm{gravity\_nav}(\mathrm{ENU})`$. Biased stream shared by free-run and filter:

```math
a_m = a_b + b_a^{\mathrm{sensor}}
,\qquad
b_a^{\mathrm{sensor}} = (0.05,\; 0,\; 0)
```

Truth itself advances with ideal IMU via `mechanize_step` (perfect sensors, zero bias state).

---

## Three trajectories compared

| Trace | Model |
| ----- | ----- |
| Truth (green) | `mechanize_step` on ideal IMU from true $`a_n,\omega_b`$ |
| Free-run (orange) | Open-loop `mechanize_step` on biased IMU; $`b_a`$ not estimated → residual double-integrated accel |
| Aided (cyan) | `InsNavigator`: `predict(imu)` each step; every 50 steps `update_pose_heading(p_{\mathrm{truth}}, \psi_{\mathrm{truth}}, R_p, R_\psi)` |

Free-run position error grows roughly as $`\tfrac12 |b_a|\,t^2`$ while the bias is unobserved. Sparse position alone weakly observes $`\delta b_a`$; the design keeps velocity process noise tied to $`\sigma_a`$ / $`\sigma_{ba}`$ so the filter can attribute residual specific force to bias rather than white velocity noise.

---

## Filter design (same as deploy)

```text
design::ins_eskf_design(
  /*gyro_nd*/ 0.003, /*accel_nd*/ 0.03,
  /*bg_rw*/ 0.0001, /*ba_rw*/ 0.001,
  dt, /*pos_std*/ 1.0)
```

Full 15-state error dynamics, discrete $`Q`$ (per-axis integrated RW — not full
van Loan of $`F`$; see accuracy note there), default $`P_0`$, and position
$`R_p = I_3`$ are documented in
[`../ins_navigator/ins_navigator_derivation.md`](../ins_navigator/ins_navigator_derivation.md).

Heading measurement covariance in this SIL:

```math
R_\psi = (0.02)^2
```

Truth heading for aiding: `ins_heading(truth, frame)` (ENU $`\psi = \operatorname{atan2}(E,N)`$ of forward).

---

## Nominal / error models (summary)

Mechanization (nominal):

```math
\begin{aligned}
\omega &= \omega_m - b_g
\\
a_b &= a_m - b_a
\\
a_n &= R(q)\,a_b + g_n
\\
v^+ &= v + a_n\Delta t
\\
p^+ &= p + v\Delta t + \tfrac12 a_n(\Delta t)^2
\\
q^+ &= q \otimes \Delta q(\omega\Delta t)
\end{aligned}
```

Error state $`\delta x = [\delta\theta;\,\delta v;\,\delta p;\,\delta b_g;\,\delta b_a]`$:

```math
\begin{aligned}
\dot{\delta\theta}
&=
-[\omega]_\times\delta\theta - \delta b_g
\\
\dot{\delta v}
&=
-R[a_b]_\times\delta\theta - R\,\delta b_a
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

---

## Outputs

Host plot: `plots/estimation/ins_eskf_3d.html` (Play animation).
RGB triad = body axes at the true pose; cyan marker = aided vehicle.

---

## Files

| File | Role |
| ---- | ---- |
| [`ins_eskf_derivation.md`](ins_eskf_derivation.md) | This derivation |
| [`ins_eskf_estimator.hpp`](ins_eskf_estimator.hpp) | Nameplate / design (DiD) |
| [`ins_eskf_sketch.cpp`](ins_eskf_sketch.cpp) | Thin navigator smoke |
| [`ins_eskf_sil.cpp`](ins_eskf_sil.cpp) | Host SIL + Plotly animation |
| [`../ins_navigator/`](../ins_navigator/) | Deploy nameplate / flashable sketch |
| [`../ins_mechanization/`](../ins_mechanization/) | Open-loop mechanization-only demo |
