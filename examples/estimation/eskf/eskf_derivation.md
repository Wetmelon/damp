# MARG attitude ESKF — design / deploy model

Teaching notes for [`eskf_estimator.hpp`](eskf_estimator.hpp) /
[`eskf_sketch.cpp`](eskf_sketch.cpp).
Design path: `design::eskf_marg` → `ESKFOrientationFilter<float, 6>`.

Canonical reference: Solà et al., *Quaternion kinematics for the error-state Kalman filter* (2017).

---

## Physical setup

| Symbol | Code | Value | Meaning |
| ------ | ---- | ----- | ------- |
| $`\Delta t`$ | `dt` | $`0.01`$ s | Sample period (100 Hz) |
| $`\sigma_g`$ | `kGyroNd` | $`0.003`$ rad/s/√Hz | Gyro noise density |
| $`\sigma_a`$ | `kAccelNd` | $`0.03`$ m/s²/√Hz | Accel noise density |
| $`\sigma_m`$ | `kMagNd` | $`0.3`$ /√Hz | Mag noise density (same units as mag meas) |
| $`\sigma_{bg}`$ | `kBgRw` | $`10^{-4}`$ rad/s^{3/2} | Gyro bias random-walk density |
| $`g`$ | default `g_nav` | $`9.81`$ m/s² | Gravity magnitude (ENU) |

Call site:

```text
design::eskf_marg(kGyroNd, kAccelNd, kMagNd, kBgRw, dt)
```

Default initial 1-σ (not overridden in the sketch): attitude $`0.1`$ rad, gyro bias $`0.01`$ rad/s.

Nav / body convention. Quaternion $`q`$ is body → world. Default gravity and field in the MARG `update`:

```math
g_n = \begin{bmatrix} 0 \\ 0 \\ -g \end{bmatrix}
\quad\text{(ENU)},
\qquad
m_n = \begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix}
```

Mock sensors at rest, $`q = I`$: specific force $`a \approx (0,0,+g)`$, gyro zero, mag $`m_b = m_n`$.

---

## Nominal state

```math
x = \{ q,\; b_g \}
```

| Quantity | Meaning |
| -------- | ------- |
| $`q`$ | Body → world unit quaternion |
| $`b_g`$ | Gyro bias [rad/s], body |

Predict (bias-corrected body rate, right-multiplicative integrate):

```math
\begin{aligned}
\omega &= \omega_m - b_g
\\
q &\leftarrow q \otimes \Delta q(\omega\,\Delta t)
\end{aligned}
```

with first-order $`\Delta q`$ as in `Quaternion::integrate_body_rates`.

---

## Error state (6)

```math
\delta x =
\begin{bmatrix}
\delta\theta \\
\delta b_g
\end{bmatrix}
\in \mathbb{R}^{6}
```

Right-multiplicative attitude error: $`q_{\mathrm{true}} = q \otimes \exp(\delta\theta)`$.

Continuous error dynamics (aligned with the attitude block of the INS ESKF):

```math
\begin{aligned}
\dot{\delta\theta}
&=
-[\omega]_\times\,\delta\theta
- \delta b_g
\\
\dot{\delta b}_g
&=
0
\end{aligned}
```

Discrete first-order map used in the filter:

```math
\delta\theta^+
\approx
\bigl(I - [\omega]_\times\Delta t\bigr)\delta\theta
- \Delta t\,\delta b_g
```

After the Kalman update, inject and reset:

```math
q \leftarrow q \otimes \exp(\delta\theta),
\qquad
b_g \leftarrow b_g + \delta b_g,
\qquad
\delta x \leftarrow 0
```

---

## Process noise $`Q`$ (discrete, continuous–discrete RW)

Per axis $`i`$ (blocks for $`\delta\theta`$ then $`\delta b_g`$; off-axis zero):

```math
\begin{aligned}
Q_{\theta\theta}
&=
\sigma_g^2\,\Delta t
+ \sigma_{bg}^2\,\frac{\Delta t^3}{3}
\\
Q_{\theta b}
&=
-\sigma_{bg}^2\,\frac{\Delta t^2}{2}
\\
Q_{bb}
&=
\sigma_{bg}^2\,\Delta t
\end{aligned}
```

(`design::detail::fill_attitude_bias_q` in `eskf.hpp`.)

Accuracy note. Per-axis first-order integrated RW (attitude + gyro bias), not a
full van Loan of the continuous attitude error $`F`$ (skew-$`\omega`$ process-noise
coupling omitted). Same class of approximation as the 15-state INS $`Q`$ in
[`../ins_navigator/ins_navigator_derivation.md`](../ins_navigator/ins_navigator_derivation.md).

---

## Measurement model

Specific force (tilt / pitch–roll), body frame:

```math
a_{\mathrm{pred}}
=
-R(q)^\top g_n
```

(matches `specific_force_at_rest`.)

Residual $`z_a = a_m - a_{\mathrm{pred}}`$; right-multiplicative Jacobian

```math
H_a = [a_{\mathrm{pred}}]_\times
```

(columns on $`\delta\theta`$; zeros on $`\delta b_g`$).

Magnetometer (yaw / heading), body field:

```math
m_{\mathrm{pred}}
=
R(q)^\top m_n
\qquad
H_m = [m_{\mathrm{pred}}]_\times
```

MARG stacks $`y = [a; m]`$ (6×1). Mag update is independent of the accel magnitude gate (default $`\gamma = 0`$: always update tilt).

Measurement noise from continuous densities:

```math
R_a = \frac{\sigma_a^2}{\Delta t}\, I_3
,\qquad
R_m = \frac{\sigma_m^2}{\Delta t}\, I_3
,\qquad
R = \operatorname{diag}(R_a, R_m)
```

Initial covariance (diagonal, design defaults):

```math
P_0(\delta\theta_i,\delta\theta_i) = (0.1)^2
,\qquad
P_0(\delta b_{g,i},\delta b_{g,i}) = (0.01)^2
```

---

## Design → deploy

1. `design::eskf_marg(σ_g, σ_a, σ_m, σ_bg, dt)` → `ESKFResult` with $`Q,R,P_0`$.
2. `static constinit ESKFOrientationFilter<float, 6> filt = kDesign.as<float>()`.
3. Loop: `estimate_period(filt, accel, gyro, mag)` → `filt.orientation().to_euler<ZYX>()`.

Observability: accel observes pitch/roll via gravity; mag observes yaw when $`m_n`$ is usable. Without mag, use `design::eskf_imu` / `NY=3`.

---

## Files

| File | Role |
| ---- | ---- |
| [`eskf_derivation.md`](eskf_derivation.md) | This derivation |
| [`eskf_estimator.hpp`](eskf_estimator.hpp) | Nameplate, design, `estimate_period` |
| [`eskf_sketch.cpp`](eskf_sketch.cpp) | Flashable sketch (mock IMU) |
| [`eskf_sil.cpp`](eskf_sil.cpp) | Finite-tick host smoke |
| `inc/damp/estimation/eskf.hpp` | Design $`Q/R/P_0`$ |
| `inc/damp/estimation/sensor_fusion.hpp` | `ESKFOrientationFilter` runtime |
