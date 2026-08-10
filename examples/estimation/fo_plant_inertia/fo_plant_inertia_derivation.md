# First-order $`J,b`$ plant — grey-box estimator model

Teaching notes for [`fo_plant_inertia_estimator.hpp`](fo_plant_inertia_estimator.hpp) /
[`fo_plant_inertia_sil.cpp`](fo_plant_inertia_sil.cpp).
Online identification of mechanical inertia and viscous friction from streaming
$`(\tau, \omega)`$ using `FirstOrderPlantEstimator`, then an application map
to $`J`$ and $`b`$ (not a separate library type).

Canonical background: Ljung, *System Identification* (RLS); Åström & Wittenmark, *Adaptive Control*.

---

## Physical plant (true)

Single rotational axis:

```math
J\,\dot\omega + b\,\omega = \tau
```

| Symbol | Code | Value | Meaning |
| ------ | ---- | ----- | ------- |
| $`J`$ | `J_true` | $`0.01`$ kg·m² | Inertia |
| $`b`$ | `b_true` | $`0.05`$ N·m·s/rad | Viscous friction |
| $`T_s`$ | `Ts` | $`0.01`$ s | Sample period |
| $`\omega`$ | `omega` | — | Angular velocity [rad/s] |
| $`\tau`$ | `torque` | — | Applied torque [N·m] |

First-order form $`\tau_m \dot y + y = K u`$ with $`y = \omega`$, $`u = \tau`$:

```math
K = \frac{1}{b}
,\qquad
\tau_m = \frac{J}{b}
```

(so $`b = 1/K`$ and $`J = \tau_m / K`$ when $`K > 0`$).

---

## Discrete plant used in the SIL

Exact ZOH discretization of $`J\dot\omega + b\omega = \tau`$ (constant $`\tau`$ over $`T_s`$):

```math
\begin{aligned}
a
&=
\exp\bigl(-b\,T_s / J\bigr)
\\
b_d
&=
\frac{1}{b}\bigl(1 - a\bigr)
\\
\omega[k]
&=
a\,\omega[k-1] + b_d\,\tau[k-1]
\end{aligned}
```

---

## Estimator model (`FirstOrderPlantEstimator`)

Identifies the grey-box continuous parameters of $`\tau_m\,\dot y + y = K\, u`$
via RLS. Config in the header: `kCfg` (forgetting $`0.999`$, initial covariance $`10^4`$).

### Excitation

PRBS-like hold pattern (`kPattern`, hold 10 samples), $`\tau = 0.5 \times \mathrm{level}`$, 1200 samples.

## Mechanical map (`map_mechanical`)

```math
b_{\mathrm{est}} = 1/K,\qquad J_{\mathrm{est}} = \tau_m / K
```

---

## Files

| File | Role |
| ---- | ---- |
| [`fo_plant_inertia_derivation.md`](fo_plant_inertia_derivation.md) | This derivation |
| [`fo_plant_inertia_estimator.hpp`](fo_plant_inertia_estimator.hpp) | Nameplate, tick, map |
| [`fo_plant_inertia_sketch.cpp`](fo_plant_inertia_sketch.cpp) | Thin smoke |
| [`fo_plant_inertia_sil.cpp`](fo_plant_inertia_sil.cpp) | Host identification demo |
| `inc/damp/estimation/parameter_estimation.hpp` | `FirstOrderPlantEstimator` |
| `inc/damp/estimation/rls.hpp` | Underlying RLS |
