# Known limitations

v0.1. Breaking renames are allowed (no deprecation cycle). See [AGENTS.md](../AGENTS.md).

| Topic | Status |
| ----- | ------ |
| API stability | 0.x — breaking renames allowed |
| Example coverage | Teaching demos for the main paths; not one example per REFERENCE row |
| H∞ synthesis / mixsyn, μ | Not implemented; `norm_hinf` and SISO / loop-at-a-time disk margins ship |
| Motor / power / motion packs | Not in this public tree; the FOC slice (`foc.hpp`, `modulation.hpp`, `spm.hpp`) is the exception |
| No C++ standard headers | The core builds on ETL with `ETL_NO_STL`, math `FREESTANDING` or `HEADER` ([README](../README.md#using-etl-instead-of-the-standard-library)). `make nostd-check` emulates this with `-nostdinc++` on GCC; also checked by hand on Clang and arm-none-eabi. Not tried on a toolchain that truly ships none, such as avr-gcc. `damp/workbench.hpp` still needs the standard library |
| Delay in a design | At most one sample, in the plant: `discretize(sys, Ts, InputDelay{θ})`. Longer dead time: sample slower, or `SmithPredictor` |

Implementation notes:

- GCC: designs of about 10 states or more at compile time need `-fconstexpr-ops-limit=67108864` (2^26); smaller designs fit the default.
- Kalman/ESKF: thin-H Joseph path; some large-Δt Φ assumptions (documented in headers).
- Host analysis: some FRF helpers continuous-only unless routed via discrete `StateSpace`.
- Model inverse (`design::output_inverse`): square only, with a vector relative degree. A zero within √ε of the unit circle is refused (the −1 sampling zero of a ZOH-sampled double integrator).
- MPC: with a hard output row or a terminal set, a MaxIterations or Degenerate iterate that fails the primal row check is held (last applied input). That hold does not restore output or terminal-set feasibility. Soft-output ticks still apply the exhausted iterate, then the u box.
