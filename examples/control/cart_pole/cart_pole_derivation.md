# Cart-pole — plant model

Teaching notes for [`cart_pole_controller.hpp`](cart_pole_controller.hpp) /
[`cart_pole_sil.cpp`](cart_pole_sil.cpp).
Read this before forking the sketch. Code numbers must match the equations here.

Related external reading (different $`I`$-based parameterization, continuous design):
[UMich CTMS — Inverted Pendulum](https://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=SystemModeling).
This folder is the design-is-deploy path (discrete LQR, nonlinear SIL, flashable sketch).

---

## Physical setup

| Symbol | Code | Value | Meaning |
| ------ | ---- | ----- | ------- |
| $`M`$ | `M` | $`1.0`$ kg | Cart mass |
| $`m`$ | `m_pole` | $`0.1`$ kg | Pole mass (point mass) |
| $`\ell`$ | `L` | $`0.5`$ m | Distance pivot → pole CoM |
| $`b`$ | `b_fric` | $`0.1`$ N·s/m | Viscous friction on the cart |
| $`g`$ | `g` | $`9.81`$ m/s² | Gravity |
| $`T_s`$ | `Ts` | $`0.01`$ s | Sample period (100 Hz control) |

States (column vector $`x`$):

```math
x =
\begin{bmatrix} p \\ \dot p \\ \theta \\ \dot\theta \end{bmatrix}
=
\begin{bmatrix}
\text{cart position [m]} \\
\text{cart velocity [m/s]} \\
\text{pole angle from upright [rad]} \\
\text{pole rate [rad/s]}
\end{bmatrix}
```

Input $`u = F`$ — horizontal force on the cart [N].

Angle convention: $`\theta = 0`$ is upright (unstable equilibrium). Positive $`\theta`$ is a small tip of the pole; gravity pulls the pole *away* from upright for small $`\theta`$.

Modeling assumptions (what this demo implements):

- Planar motion; no wheel slip model — force $`F`$ acts directly on the cart.
- Pole is a point mass at distance $`\ell`$ (no separate moment of inertia $`I`$ about CoM).
- Friction only on the cart ($`b \dot p`$); no damping on $`\dot\theta`$.
- No actuator dynamics, no encoder quantization in the plant.

The CTMS page keeps a rod inertia $`I`$ and a common denominator $`p = I(M+m)+M m \ell^2`$. This demo uses the lighter “point-mass pole” form below; both are standard teaching plants.

---

## Nonlinear dynamics (SIL plant)

Used in `cart_pole_sil.cpp` as `cart_pole_f`. Let $`s = \sin\theta`$, $`c = \cos\theta`$.

```math
D = M + m\, s^2
```

```math
\begin{aligned}
\ddot p
&=
\frac{
  F - b\,\dot p + m\, s\,(\ell\,\dot\theta^2 - g\, c)
}{D}
\\
\ddot\theta
&=
\frac{\ddot p\, c}{\ell}
+
\frac{(M+m)\, g\, s}{\ell\, D}
\end{aligned}
```

First-order form for integration: $`\dot x = f(x,u)`$ with $`\dot p = x_2`$, $`\dot\theta = x_4`$.

---

## Linearization at upright ($`\theta = 0`$)

At the origin $`x = 0`$, $`u = 0`$: $`s \approx \theta`$, $`c \approx 1`$, $`\dot\theta = 0`$, $`D \to M`$.

Jacobian of the nonlinear plant at the origin (matches `linearize_upright()`).
Gravity couples into cart acceleration ($`a_{23}`$); pole-rate row inherits that
through $`\ddot\theta \supset \ddot p\, c/\ell`$, which simplifies $`a_{43}`$ to
$`g/\ell`$:

```math
\begin{aligned}
\ddot p
&=
-\frac{b}{M}\,\dot p
-
\frac{m g}{M}\,\theta
+
\frac{1}{M}\, F
\\
\ddot\theta
&=
-\frac{b}{M \ell}\,\dot p
+
\frac{g}{\ell}\,\theta
+
\frac{1}{M \ell}\, F
\end{aligned}
```

In matrix form $`\dot x = A x + B u`$:

```math
A =
\begin{bmatrix}
0 & 1 & 0 & 0 \\
0 & -b/M & -m g/M & 0 \\
0 & 0 & 0 & 1 \\
0 & -b/(M\ell) & g/\ell & 0
\end{bmatrix}
,\qquad
B =
\begin{bmatrix}
0 \\
1/M \\
0 \\
1/(M\ell)
\end{bmatrix}
```

With the nameplate numbers:

```math
A =
\begin{bmatrix}
0 & 1 & 0 & 0 \\
0 & -0.1 & -0.981 & 0 \\
0 & 0 & 0 & 1 \\
0 & -0.2 & 19.62 & 0
\end{bmatrix}
,\qquad
B =
\begin{bmatrix}
0 \\ 1 \\ 0 \\ 2
\end{bmatrix}
```

```math
\begin{aligned}
a_{23}
&=
-m g / M
=
-0.981
\\
a_{43}
&=
g/\ell
=
9.81 / 0.5
=
19.62
\end{aligned}
```

Outputs (optional sensors for $`C`$ — design uses full state for LQR):

```math
C =
\begin{bmatrix}
1 & 0 & 0 & 0 \\
0 & 0 & 1 & 0
\end{bmatrix}
,\qquad
y =
\begin{bmatrix} p \\ \theta \end{bmatrix}
```

Open-loop: $`A`$ has a right-half-plane eigenvalue from the upright pole — the plant is open-loop unstable (as expected).

---

## Transfer functions

Cart position and pole angle are coupled through gravity ($`a_{23}`$), so the
SISO maps $`p(s)/F(s)`$ and $`\theta(s)/F(s)`$ are the corresponding rows of

```math
G(s) = C(sI - A)^{-1} B
```

with the $`A,B,C`$ above (fourth-order, one RHP pole from upright $`a_{43}=g/\ell`$).
LQR design uses the state-space form directly; expand $`G(s)`$ only if you need Bode
slices of a single output.

---

## Control design (what `cart_pole_controller.hpp` does)

1. Continuous plant $`(A,B)`$ as above.
2. Discrete LQR via `design::discrete_lqr_from_continuous(A, B, Q, R, Ts)` (ZOH-style design path in the library).
3. Cost weights (diagonal $`Q`$, $`R = 1`$):

| State | $`Q_{ii}`$ | Intent |
| ----- | ---------- | ------ |
| $`p`$ | 10 | Moderate cart regulation |
| $`\dot p`$ | 0.1 | Light velocity penalty |
| $`\theta`$ | 100 | Keep the pole upright tightly |
| $`\dot\theta`$ | 10 | Damp pole rate |
| $`F`$ | $`R = 1`$ | Unit control effort |

4. Law: $`u = -K x`$ (`LQR::control`).
5. Sketch uses `lqr_d.as<float>()`; SIL may keep double. Same $`K`$ source.

Tuning tip: raise $`Q_{33}`$ (angle) for stiffer upright recovery; raise $`R`$ if force saturates.

---

## Files

| File | Role |
| ---- | ---- |
| [`cart_pole_derivation.md`](cart_pole_derivation.md) | This derivation |
| [`cart_pole_controller.hpp`](cart_pole_controller.hpp) | Nameplate, $`A,B`$, discrete LQR, `control_period` |
| [`cart_pole_sketch.cpp`](cart_pole_sketch.cpp) | Flashable sketch |
| [`cart_pole_sil.cpp`](cart_pole_sil.cpp) | Nonlinear plant + closed-loop SIL + plot |

Host plot: `examples/plots/control/cart_pole_lqr.html`.
