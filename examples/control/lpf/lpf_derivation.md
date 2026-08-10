# Low-pass filter design — model note

Teaching notes for [`lpf_filter.hpp`](lpf_filter.hpp).

There is no plant ODE in this file — only filter models: continuous
prototypes discretized for runtime `LowPass` objects, then a step-response smoke
test on the first-order filter.

---

## Nameplate (as in the example)

| Quantity | Code | Value |
| -------- | ---- | ----- |
| Cutoff $`f_c`$ | `10.0f` | 10 Hz |
| Sample time $`T_s`$ | `0.001f` | 1 ms (1 kHz) |
| 2nd-order damping $`\zeta`$ | `0.707f` | Butterworth-ish |

Objects:

- `design::lowpass_1st(fc, Ts)` → `LowPass<1, float>`
- `design::lowpass_2nd(fc, Ts, zeta)` → `LowPass<2, float>` from $`(b_0,b_1,b_2),(a_1,a_2)`$
- Convenience ctor `LowPass<1,float>{fc, Ts}` (same first-order design)

---

## Continuous prototypes

First order — cutoff $`\omega_c = 2\pi f_c`$, time constant $`\tau = 1/\omega_c`$:

```math
H_1(s)
=
\frac{1}{\tau s + 1}
=
\frac{\omega_c}{s + \omega_c}
```

Second order (natural frequency $`\omega_c`$, damping $`\zeta`$):

```math
H_2(s)
=
\frac{\omega_c^2}{s^2 + 2\zeta\omega_c s + \omega_c^2}
```

With $`f_c = 10`$ Hz: $`\omega_c \approx 62.83`$ rad/s. $`\zeta \approx 0.707`$ gives a
maximally flat passband for the second-order low-pass.

---

## Discrete design (library path)

`design::lowpass_1st` / `lowpass_2nd` use Tustin (bilinear) discretization
$`s \leftarrow (2/T_s)\,(1 - z^{-1})/(1 + z^{-1})`$, with numerator taps forced so
unit DC gain holds by construction (robust under `-ffast-math` reassociation).

First-order difference equation shape:

```math
y[k]
=
b_0\, x[k] + b_1\, x[k-1] - a_1\, y[k-1]
```

with $`H(z=1) = 1`$ (DC pass-through). Second order adds $`b_2, a_2`$ taps.

Runtime: `lpf(1.0f)` steps a unit input; the example prints every 10th sample of
the step response.

---

## Files

| File | Role |
| ---- | ---- |
| [`lpf_derivation.md`](lpf_derivation.md) | This short filter model |
| [`lpf_filter.hpp`](lpf_filter.hpp) | Nameplate + `design::lowpass_*` coeffs + `filter_period` |
| [`lpf_sketch.cpp`](lpf_sketch.cpp) | Flashable filter smoke |
| [`lpf_sil.cpp`](lpf_sil.cpp) | Host step-response print (no plant ODE) |

No plant plot; console step samples only.
