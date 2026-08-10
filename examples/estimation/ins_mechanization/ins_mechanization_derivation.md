# Strapdown INS mechanization — open-loop model

Teaching notes for [`ins_mechanization_estimator.hpp`](ins_mechanization_estimator.hpp) /
[`ins_mechanization_sil.cpp`](ins_mechanization_sil.cpp).
Host-only. No Kalman filter: IMU samples → $`p,v,q`$ by dead-reckoning
(`mechanize_step`). Bias peel-away motivates absolute aiding (see
[`../ins_eskf/`](../ins_eskf/)).

Canonical references: Solà (ESKF kinematics, 2017); Groves (strapdown mechanization).

---

## Physical setup

| Symbol | Code | Value | Meaning |
| ------ | ---- | ----- | ------- |
| Frame | `kFrame` | ENU | East–North–Up |
| $`\Delta t`$ | `kTsSil` | $`1/60`$ s | Integration step (~60 Hz animation) |
| $`t_{\mathrm{end}}`$ | `t_end` | $`12`$ s | Scenario length |
| $`b_a^{\mathrm{sensor}}`$ | `sensor_ba` | $`(0.05, 0, 0)`$ m/s² | Constant body-x accel bias on one stream |

Truth motion schedule (same shape as the ESKF SIL, shorter horizon):

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

Initial: $`q = I`$, $`v = p = 0`$, nominal biases zero.

---

## Frames and gravity

| `NavFrame` | Axes $`(x,y,z)`$ | Gravity $`g_n`$ |
| ---------- | ---------------- | --------------- |
| NED | North, East, Down | $`(0,0,+g)`$ |
| ENU (this demo) | East, North, Up | $`(0,0,-g)`$ |

$`g = 9.80665`$ m/s² (`kStandardGravity`). Orientation $`q`$ maps body → nav:
`q.rotate(v_body)` is in nav components.

---

## IMU measurement model (synthetic)

Body-frame sample:

```math
\begin{aligned}
\omega_m &= \omega_b
\\
a_b &= R(q)^\top \bigl(a_n - g_n\bigr)
\end{aligned}
```

so that $`a_n = R(q)\,a_b + g_n`$ recovers the commanded nav acceleration.
At rest ($`a_n = 0`$): $`a_b = -R^\top g_n`$ (`specific_force_at_rest`); with
$`q = I`$ in ENU, $`a_b \approx (0,0,+g)`$.

Biased free-run path adds a constant sensor offset not stored in `InsState.b_a`:

```math
a_m = a_b + (0.05,\; 0,\; 0)
```

---

## Mechanization step

State: `InsState` $`\{q, v, p, b_g, b_a\}`$. One sample
(`mechanize_step` / `estimate_period`):

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
p^+ &= p + v\,\Delta t + \tfrac12 a_n\,(\Delta t)^2
\\
q^+ &= q \otimes \Delta q(\omega\,\Delta t)
\end{aligned}
```

Attitude integrate is first-order right-multiplicative
(`Quaternion::integrate_body_rates`). Flat Earth: no Earth rate, transport rate, or WGS84 curvature.

Modeling assumptions (what this demo implements):

- Sensors fixed to the body (strapdown); rotation done in software.
- Constant gravity vector in local level; no Coriolis / Schuler terms.
- Open-loop only — no GPS, ZUPT, or heading aid.
- Ideal path uses perfect IMU; orange path uses accel bias with $`b_a = 0`$ in the state (uncompensated).

---

## Paths compared

| Trace | Integration |
| ----- | ----------- |
| Truth / ideal mechanization (green) | Ideal IMU; tracks truth by construction |
| Open-loop + accel bias (orange) | Same recipe on $`a_m = a_b + b_a^{\mathrm{sensor}}`$; position drifts |

Without estimating $`b_a`$, residual body-x specific force integrates twice into position.
Absolute aiding is the ESKF half (`InsNavigator`); this example stops at the open-loop recipe.

---

## Outputs

Host plot: `plots/estimation/ins_mechanization_3d.html`.
RGB triad = body axes at the true vehicle; ground square on $`z = 0`$.

---

## Files

| File | Role |
| ---- | ---- |
| [`ins_mechanization_derivation.md`](ins_mechanization_derivation.md) | This derivation |
| [`ins_mechanization_estimator.hpp`](ins_mechanization_estimator.hpp) | Nameplate, synthetic IMU, tick |
| [`ins_mechanization_sketch.cpp`](ins_mechanization_sketch.cpp) | Thin rest smoke |
| [`ins_mechanization_sil.cpp`](ins_mechanization_sil.cpp) | Host animation |
| `inc/damp/estimation/ins_mechanization.hpp` | `mechanize_step`, frames, `InsState` |
| [`../ins_eskf/`](../ins_eskf/) | Same motion with position + heading ESKF |
