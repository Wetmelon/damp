# Position LQI with a clamp

The shortest path that still shows design-is-deploy: a mass-spring-damper, `design::lqi`, an actuator that cannot do what the law asked, and `back_calculate` so the integrator sees the force that actually went out.

Numbers match [`lqi_controller.hpp`](lqi_controller.hpp) ($`m`$, $`k`$, $`c`$). The sketch and the SIL call [`control_period`](lqi_controller.hpp): update → control → clamp → `back_calculate` → predict.

## Plant

Mass, spring, damper. The state is what you measure: $`x = [\mathrm{position},\;\mathrm{speed}]`$, $`u`$ is force, $`y`$ is position.

```math
\dot x =
\begin{bmatrix} 0 & 1 \\ -k/m & -c/m \end{bmatrix} x
+ \begin{bmatrix} 0 \\ 1/m \end{bmatrix} u,\qquad
y = \begin{bmatrix} 1 & 0 \end{bmatrix} x
```

Numbers in [`lqi_controller.hpp`](lqi_controller.hpp): $`m = 1\,\mathrm{kg}`$, $`k = 1\,\mathrm{N/m}`$, $`c = 0.5\,\mathrm{N\cdot s/m}`$. Same input-output map as $`G(s) = 1/(m s^2 + c s + k)`$. `discretize(sys, 0.01)` is the 100 Hz model the design and the SIL share. The actuator is $`|u| \le 1.2`$ (steady force at $`r = 1`$ is $`k = 1`$).

## Why LQI

A step r = 1 saturates. Without `back_calculate`, the integrator keeps charging on the unsaturated request and overshoots. With `back_calculate(u_sat)`, predict holds the integrator while the clamp is against it, then the loop settles on r.

LQR with `TrackedOutputs` would hold r only if the model and the plant match and nothing loads the input. LQI is the path when they will not.

## Tick

Low-latency ISR, the same objects on the board and in SIL:

```cpp
law.update(x);
auto u = law.control(r);
auto u_sat = clamp(u, u_min, u_max);
law.back_calculate(u_sat);
law.predict();
write(u_sat);
```
