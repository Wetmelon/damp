# Encoder / tach velocity estimation — plant and estimators

Teaching notes for [`encoder_velocity_estimator.hpp`](encoder_velocity_estimator.hpp) /
[`encoder_velocity_sil.cpp`](encoder_velocity_sil.cpp).
Host comparison of four velocity estimators on a quantized position measurement
at an 8 kHz loop rate. Positions in turns (revolutions); velocities in turns/s.

---

## Physical setup

| Symbol | Code | Value | Meaning |
| ------ | ---- | ----- | ------- |
| $`\Delta t`$ | `dt` | $`1/8000`$ s | Sample period (8 kHz current-loop rate) |
| $`N_c`$ | `counts_per_rev` | scenario | Encoder counts per revolution |
| $`q`$ | `q` | $`1/N_c`$ | Quantization step [turns/count] |
| $`A`$ | `amplitude` | scenario | Motion amplitude [turns] |
| $`f`$ | `freq` | scenario | Motion frequency [Hz] |
| $`\omega`$ | `w` | $`2\pi f`$ | Angular frequency [rad/s] of the sinusoid |

True continuous motion (scalar axis in turns):

```math
\begin{aligned}
\theta(t)
&=
A\,\sin(\omega t)
\\
\dot\theta(t)
&=
A\,\omega\,\cos(\omega t)
\\
\ddot\theta(t)
&=
-A\,\omega^2\,\sin(\omega t)
\end{aligned}
```

Encoder (quantized position):

```math
\theta_m(t)
=
q \cdot \operatorname{round}\bigl(\theta(t)/q\bigr)
,\qquad
q = 1/N_c
```

Peak counts per sample scale as $`|\dot\theta|_{\max}\Delta t / q = A\omega\Delta t / q`$.
Near reversals, motion is a fraction of a count per sample — raw differencing is mostly quantization noise.

---

## Scenarios (as coded)

| Scenario | $`N_c`$ | PLL $`\mathrm{bw}`$ [rad/s] | $`A`$ [turns] | $`f`$ [Hz] | LPF $`\tau`$ [s] |
| -------- | ------- | --------------------------- | ------------- | ---------- | ---------------- |
| 14-bit absolute | 16384 | 1000 | 0.25 | 1.0 | $`3\times 10^{-3}`$ |
| 6-state hall × 7 pole pairs | 42 | 100 | 1.0 | 0.5 | $`8\times 10^{-3}`$ |

Metrics: RMS velocity error after settle, and the same restricted to near-reversal samples $`|\dot\theta| < 0.15\,A\omega`$.

---

## Estimators

### 1. Raw finite difference

```math
v_{\mathrm{fd}}[k]
=
\frac{\theta_m[k] - \theta_m[k-1]}{\Delta t}
```

Zero lag; high quantization noise at low speed.

### 2. One-pole LPF of the finite difference

```math
\begin{aligned}
\alpha
&=
\frac{\Delta t}{\tau + \Delta t}
\\
v_{\mathrm{lpf}}
&\leftarrow
v_{\mathrm{lpf}} + \alpha\bigl(v_{\mathrm{fd}} - v_{\mathrm{lpf}}\bigr)
\end{aligned}
```

### 3. Critically-damped PLL (`PllObserver` in the header)

```math
k_p = 2\,\mathrm{bw}
,\qquad
k_i = 0.25\,k_p^2
```

### 4. Levant robust exact differentiator

Library type `RobustExactDifferentiator`. Bound $`L = 2\,A\,\omega^2`$ in the SIL.

---

## Outputs

| Plot | Scenario |
| ---- | -------- |
| `plots/estimation/encoder_velocity_absolute.html` | 14-bit absolute |
| `plots/estimation/encoder_velocity_hall.html` | 42 cpr hall |

---

## Files

| File | Role |
| ---- | ---- |
| [`encoder_velocity_derivation.md`](encoder_velocity_derivation.md) | This derivation |
| [`encoder_velocity_estimator.hpp`](encoder_velocity_estimator.hpp) | Nameplate, PLL, tick helpers |
| [`encoder_velocity_sketch.cpp`](encoder_velocity_sketch.cpp) | Thin smoke |
| [`encoder_velocity_sil.cpp`](encoder_velocity_sil.cpp) | Host comparison + plots |
| `inc/damp/filters/differentiator.hpp` | `RobustExactDifferentiator` |
