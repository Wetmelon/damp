# Hosted SIL plot gallery

Interactive Plotly HTML written by example SILs (maintainer runs: `make examples`).
Paths are relative to `examples/` CWD. Gitignored — regenerate locally; this
index lists the default plot names so you can jump from plot → example folder.

```bash
make examples          # build + run stale examples (writes plots/)
make examples-force    # wipe stamps + re-run
```

Open a file in a browser after generation. Shared Plotly JS: `js/plotly.min.js`.

## control/

| Plot | Example folder |
| ---- | -------------- |
| `cart_pole_lqr.html` | [`../control/cart_pole/`](../control/cart_pole/) |
| `pendulum_sim.html`, `pendulum_phase.html`, `pendulum_sim_high_q.html` | [`../control/pendulum/`](../control/pendulum/) |
| `input_shaper_resonance.html` | [`../control/input_shaper/`](../control/input_shaper/) |
| `multirate_rig.html` | [`../control/multirate_rig/`](../control/multirate_rig/) |


## estimation/

| Plot | Example folder |
| ---- | -------------- |
| `imu_pose_3d.html` | [`../estimation/imu_pose/`](../estimation/imu_pose/) |
| `ins_eskf_3d.html` | [`../estimation/ins_eskf/`](../estimation/ins_eskf/) |
| `ins_mechanization_3d.html` | [`../estimation/ins_mechanization/`](../estimation/ins_mechanization/) |
| `encoder_velocity_*.html` | [`../estimation/encoder_velocity/`](../estimation/encoder_velocity/) |

ESKF finite-tick smoke has no HTML plot; see [`../estimation/eskf/`](../estimation/eskf/).

