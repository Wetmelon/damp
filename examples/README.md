# Examples

Examples are the main entry point for this public tree: fork a design-is-deploy
demo under `control/` or `estimation/` and adjust nameplate / plant / rates.

Host demos and hybrid SILs. Basenames are unique across folders so the tup
build drops everything in `build/obj/` and `build/*.exe`. Consumers do not need
tup — open an example folder as a reference and build with your own tree;
`make examples` is for Damp maintainers.

## Start here

| Domain | Folder | Reference layout |
| ------ | ------ | ---------------- |
| Classical pendulum / PID / LQR | `control/` | [`control/cart_pole/`](control/cart_pole/) |
| IMU / INS / ESKF | `estimation/` | [`estimation/eskf/`](estimation/eskf/) |

Also useful: [`control/pid/`](control/pid/), [`control/pendulum/`](control/pendulum/),
[`control/adrc/`](control/adrc/) (ADRC vs PI-D SIL),
[`estimation/ins_navigator/`](estimation/ins_navigator/).

Motor, power, and full motion pack demos are not in this public tree.

## Design-is-deploy layout

Fence = file split between flashable deploy and host-only SIL:

```text
examples/<domain>/<name>/
  <name>_<role>.hpp       # deploy header — role by type
  <name>_sketch.cpp       # copy onto target
  <name>_sil.cpp          # host plant + plots
  <name>_derivation.md    # plant notes (optional)
```

| Role | When |
| ---- | ---- |
| `_controller` | LQR, PID, cascade |
| `_estimator` | ESKF, INS, pose / state estimator |
| `_filter` | LPF, biquad, sense conditioning |

Example: [`control/cart_pole/`](control/cart_pole/) (`cart_pole_controller.hpp`).
Full rules: [AGENTS.md](../AGENTS.md). Register new folders in `examples/Tupfile.lua`.

## Deploy-tick contract

SIL validates the same objects that flash:

1. Deploy header owns nameplate, design, `constinit` runtime objects, and a
   single period API (`control_period`, `estimate_period`, …).
2. Sketch includes the deploy header; thin `setup`/`loop` with mock I/O only.
3. SIL includes the same deploy header and calls that period API — never a second
   cascade or copy-pasted law for the sim.
4. Plant may differ from the design model; document both in `*_derivation.md` when they do.
5. Boundary is the sensor/actuator cast (double plant ↔ float control) when host
   numerics stay double.

Anti-pattern: `*_sil.cpp` re-derives `Kp`/`Ki` without calling the header tick.
Reference layouts: `control/cart_pole/`, `estimation/eskf/`.

## Layout

| Folder | Content |
| ------ | ------- |
| `control/` | Classic control example folders (`cart_pole/`, `pid/`, …) |
| `estimation/` | `eskf/`, `ins_navigator/`, `ins_eskf/`, … |
| *(root)* | Shared helpers (`animate_*.hpp`, `damp_profile.hpp`) + small utilities |

## Plots

HTML lives under `plots/<domain>/` (relative to `examples/` CWD):

```text
plots/
  js/            # shared plotly.min.js
  control/
  estimation/
```

Tup declares plot paths in `Tupfile.lua` `plot_outputs` — update that map when
you add or rename a default plot.

```bash
make examples          # build + run stale examples
make examples-force    # wipe stamps + plot trees, full re-run
```
