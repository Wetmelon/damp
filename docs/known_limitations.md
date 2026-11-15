# Known limitations

0.x: breaking renames are allowed (no deprecation cycle). See [AGENTS.md](../AGENTS.md). Compiler bugs that change results are in [known_issues.md](known_issues.md).

| Topic | Status |
| ----- | ------ |
| API stability | 0.x — breaking renames allowed |
| Example coverage | Teaching demos for the main paths; not one example per REFERENCE row |
| H∞ synthesis / mixsyn, μ | Not implemented; `norm_hinf` and SISO / loop-at-a-time disk margins ship |
| Motor / power / motion packs | Not in this public tree; the FOC slice (`foc.hpp`, `modulation.hpp`, `spm.hpp`) is the exception |
| No C++ standard headers | The core builds on ETL with `ETL_NO_STL`, math `FREESTANDING` or `HEADER` ([README](../README.md#using-etl-instead-of-the-standard-library)). `make nostd-check` emulates this with `-nostdinc++` on GCC; also checked by hand on Clang 19 and 23 and arm-none-eabi. Not tried on a toolchain that truly ships none, such as avr-gcc. `damp/workbench.hpp` still needs the standard library |
| Delay in a design | In the plant, on the input and on the measurement: `discretize(sys, Ts, InputDelay<N>{θ}, OutputDelay<M>{θ})`, NX + N·NU + M·NY states, one delay per channel or one for all. The sampled LQ designs, `kalman` and `steady_state_target` take the continuous plant and the `InputDelay` in the design, `design::lqg(sys, Ts, InputDelay<N>{θ}, weights, noise)`; MPC and OFMPC are designed on the model. An `InputDisturbance` on a delayed plant reaches it through the registers, and its walk moves x only through them: exact for a delay of a sample or more, it leaves out the part of the last period after a shorter delay's switch. A `constexpr` design's `Law` pays only for the registers' nonzero terms; a runtime holding its matrices pays a dense square of the state per tick for a filter. Every channel gets the longest delay's registers. A fraction of a sample on the measurement is exact for the noise held over each period (discretize's convention for G), not for white noise of a density. Many samples of dead time: sample slower, or `SmithPredictor`. A delay that changes from tick to tick is not modelled |

Implementation notes:

- GCC: designs of about 10 states or more at compile time need `-fconstexpr-ops-limit=67108864` (2^26); smaller designs fit the default. Clang's default `-fconstexpr-steps` runs out sooner, at a few states (a 2-state plant with `InputDelay<5>` is 7): raise it, `-fconstexpr-steps=67108864`.
- GCC 11: `Mat2 A = {{…}}` does not deduce its scalar type; write `Mat2<double>`. GCC 11.2 and older cannot assign a `std::optional` at compile time, which some designs do.
- Clang 23 with `-ffast-math` on x86: float `asin` and `acos` reach 1.4e-6 near |x| = 1. Clang 19 and non-fast-math builds hold the documented 3.1e-7 and 4.8e-7.
- Kalman/ESKF: thin-H Joseph path; some large-Δt Φ assumptions (documented in headers).
- Host analysis: some FRF helpers continuous-only unless routed via discrete `StateSpace`.
- Model inverse (`design::output_inverse`): square only, with a vector relative degree. A zero within √ε of the unit circle is refused (the −1 sampling zero of a ZOH-sampled double integrator).
- MPC: with a hard output row or a terminal set, a MaxIterations or Degenerate iterate that fails the primal row check is held (last applied input). That hold does not restore output or terminal-set feasibility. Soft-output ticks still apply the exhausted iterate, then the u box.

## Stack

Matrices are passed and returned by value, so a routine's peak stack is a multiple of the matrices it works on. A `constexpr` design costs the target nothing; a design run at start-up on the target needs the stack below, and so does a recursive filter's tick. Worst-case call chain in float on Cortex-M4F (arm-none-eabi-g++ 14.2, `-O2`), one input and one output, C library frames not counted:

| Call | 2 states | 4 states | 8 states | 12 states |
| ---- | -------- | -------- | -------- | --------- |
| `discretize(sys, Ts)` | 1.7 KB | 4.2 KB | 12.2 KB | 24.9 KB |
| `design::lqr(sys, Ts, config)` | 5.6 KB | 15.1 KB | 47.2 KB | 97.7 KB |
| `design::kalman(sys, Ts, noise)` | 2.8 KB | 9.6 KB | 37.7 KB | 83.7 KB |
| `mat::expm(A)` | 0.6 KB | 2.4 KB | 8.1 KB | 17.8 KB |
| `mat::solve(A, b)` | 0.3 KB | 0.8 KB | 2.1 KB | 4.5 KB |
| `KalmanFilter::update(y)` | 0.2 KB | 0.8 KB | 1.8 KB | 3.7 KB |
| `KalmanFilter::predict(u)` | under 0.1 KB | 0.3 KB | 0.8 KB | 1.8 KB |

The sampled designs take most of theirs in the matrix exponential of the Van Loan block, 2·(NX + NU) rows for `lqr` and 2·NX for `kalman`. `KalmanFilter::update` holds about six NX×NX temporaries. A fixed-gain law or observer (`LQR`, `LQG`, `FixedGainObserver`) holds vectors only.
