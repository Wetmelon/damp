# Changelog

All notable changes to this project are documented here.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project uses [Semantic Versioning](https://semver.org/) for tags.

0.x releases may break API without a deprecation cycle (see [AGENTS.md](AGENTS.md)
and [known limitations](docs/known_limitations.md)).

## [Unreleased]

### Added

- `InputRateLimits{.lower, .upper}`, the limits on the move Δu per sample, is an option of `design::state_feedback` (every LQ law) as of `design::mpc`: an LQ law's command moves at most that far from the input applied last tick, inside its `InputLimits`, and its integrators, filter and delay states hold against that clamp as against one reported with `back_calculate`. A law without one is the type and size it was.
- `MeasuredDisturbance{.Bd}` is an option of `design::state_feedback` (every LQ law), `design::kalman` and `design::steady_state_target` as it is of `design::mpc`, and goes with a `KalmanNoise` and a disturbance spec in both families: the law's target is the steady state under the d it reads each tick (`steady_state_target`'s, its `target_d` per unit [d̂; d]), `step(r, d, y)` after the reference and its feedforward (`step(d, x)` on a regulator), and a filter takes d as a known input after u, so d̂ estimates only the part of a load it does not measure. It moves from `mpc.hpp` to `design/disturbance.hpp`; a design without it is the type and size it was.
- `MPC` and `OFMPC` take a feedforward from a planned trajectory, as the LQ laws do: `control(r, x_ref, u_ref)` and `step(r, x_ref, u_ref, y)`, or with a measured disturbance `control(r, x_ref, u_ref, d)` and `step(r, x_ref, u_ref, d, x)`, pose the QP about the target for r plus (x_ref, u_ref), the delay states at u_ref, so a moving reference is handed in on the plant's own readings instead of a frame shifted by hand. The output and state limits stay absolute; `control(r)` costs what it did.
- `design::mpc` and `design::steady_state_target` take a `TrackedOutputs{{i, …}}` without a disturbance spec, as `state_feedback` does: r is on the outputs it lists, held through the steady-state target, while the MPC's QP weighs every output, each of the others toward the value that steady state gives it. With a disturbance spec, its `.outputs` are the held ones and a `TrackedOutputs` does not compile, with the same message from all three designs.
- `design::mpc` takes an `InputDelay` on the sampled form, as `state_feedback`, `kalman` and `steady_state_target` do: `design::mpc<10, 5>(plant, Ts, weights, noise, design::InputDisturbance{}, InputDelay<3>{0.022})` designs on the plant whose input arrives late, with the plant's per-second weights and noise densities, the output cost integrated exactly over each period. A measured-state `MPC` takes the plant's state and holds the delay states itself, as the LQ laws do; an `OFMPC`'s filter carries them. The guide and the limitations lead with the option; designing on `discretize(sys, Ts, InputDelay<N>{θ})` remains for a model one holds.
- `OutputDelay<M>{θ}` is a trailing option of `design::kalman(plant, Ts, noise, …)`, of `design::state_feedback(plant, Ts, weights, noise, …)` (LQG, LQGI, OFLQG) and of `design::mpc<NP, NC>(plant, Ts, weights, noise, …)` (OFMPC), and of `steady_state_target`: the filter is on `discretize(plant, Ts, [InputDelay,] OutputDelay)`, the delay states of the readings already taken after the plant's, its process noise and a disturbance spec's d exact up to when each reading was taken, so `update(y)` takes the reading as it arrives and the estimate is the current state. The gain is the plant's and the QP's model the plant's with its input delay, and neither reads those delay states; the target holds them at the steady output. The guide and the limitations lead with the option. It goes with an `InputDelay`; on a discrete plant, or on a law without a `KalmanNoise`, it is refused at compile time.
- `StateSpace{.A = {{…}}, .B = {{…}}, .C = {{…}}, .G = Matrix<NX, NX>::identity(), .H = Matrix<NY, NY>::identity()}` deduces its noise inputs from the named G and H, with or without `.D` and `.Ts`, and so does `StateSpace{.A = A, .B = B, .C = C, .G = G, .H = H}` from named matrices: a plant with noise inputs no longer names its sizes.
- `state_gain()`, `delay_gain()` and `integral_gain()` on a state-feedback design's result, `state_gain()` on an `LQRResult`: K's columns on the plant's states, on the delay states of a plant whose input arrives late (the commands still on their way) and on the integrators, so a law that reads its gain needs no index into [x; v; ξ].

### Changed

- Breaking: `design::mpc` with no `TrackedOutputs` and no disturbance spec is a regulator, as `design::state_feedback` is: its `MPC` takes no r, `control()` and `step(x)` (with a measured disturbance `control(d)` and `step(d, x)`, about the steady state that holds every output at 0), and stores no map of r. An MPC that took r on every output, as 0.3's did by default, spells it `design::mpc<NP, NC>(plant, weights, TrackedOutputs{})`, with the same commands; `MPC<…, NH>` counts the outputs it holds, 0 for a regulator, every output by default. `design::steady_state_target` with no list, spec or measured disturbance holds nothing and no longer compiles: `TrackedOutputs{}`.
- Breaking for a plant with more outputs than inputs: a steady-state target holds at most one output per input, so a `TrackedOutputs` (`TrackedOutputs{}` too) or a disturbance spec's `.outputs` that holds more, and a regulator under a `MeasuredDisturbance`, no longer compile, where the target was a least-squares compromise that held none of them. Name the outputs r holds, `TrackedOutputs{{i, …}}` (r = 0 to regulate them) or the spec's `.outputs = {i, …}`.
- Every LQ law and MPC reads its run-time limits with `limits()` and sets them with `set_limits(lim)` (a breaking change): the LQ laws' public `u_lim`, `integral_lim`, `demand_lim` and `disturbance_lim` and the OFMPC's `disturbance_lim` move into `design::StateFeedbackLimits` (`input`, `input_rate`, `integral`, `demand`, `disturbance`) and `design::MPCLimits` (`input`, `input_rate`, `output`, `state`, `disturbance`, for 0.3's `u_lim`, `du_lim`, `y_lim` and `x_lim`), so `law.u_lim = {…}` becomes `auto l = law.limits(); l.input = {…}; law.set_limits(l);`, which refuses a NaN or infinite end, a lower end above its upper one or a demand limit the law cannot hold, and changes nothing. An LQ design's result holds the same struct as `limits` in place of its `u_lim`, `integral_lim` and `demand_lim` (`res.limits.input`), and so does an `LQRResult`. A `constexpr` design's `Law` refuses a demand limit other than its own and an MPC with an `InvariantSet` a change of its QP's limits, its `disturbance` still moving; the laws' sizes are unchanged.
- `design::mpc` loses its third template parameter, the terminal set's room, which the `InvariantSet` option now carries: `mpc<NP, NC, NT_max>(…)` with an invariant terminal set becomes `mpc<NP, NC>(…, InvariantSet<NT_max>{…})`. An `InvariantSet` with no room, with NC < NP or with a `StateLimits` no longer compiles, where the design refused it at run time.
- `design::mhe` takes its QP budget as a trailing `QPIterations{n}`, as `design::mpc` does, and `MHEConstraints` keeps only `x_lim`: `mhe<NH>(plant, noise, MHEConstraints{.x_lim = …, .max_qp_iterations = n})` becomes `mhe<NH>(plant, noise, MHEConstraints{.x_lim = …}, QPIterations{n})`, and `QPIterations{0}` is refused (`invalid_parameter`) where a 0 meant the default. `QPIterations` moves to `design/qp.hpp`, which `mpc.hpp` and `mhe.hpp` export.
- An offset-free MPC with an `InputDelay` carries the delay states per input in its type, `OFMPC<…, NXC, NYC, NRU, NH>` as `MPC`'s, so its `x_ref` is on the plant's states as the LQ laws' is; `OFMPC law{res}` and `Law<kRes, T>` deduce it.

- An MPC's weights and output limits are on every output of the plant, `MPCWeights<NU, NY>` sized from it and `OutputLimits` indexing the rows of C, and a disturbance spec's `.outputs` names the outputs r holds, as in the OFLQG, no longer the only outputs the QP sees. A design that listed outputs and weighed only those gives them their weights and the others zeros, the same design.
- Ru weighs the input's deviation from the steady-state target's input, the steady input that holds r, so an Ru weight no longer leaves an offset and an offset-free MPC runs as one `step(r, y)`. `set_input_target` replaces that anchor with an operating point, and `clear_input_target()` returns to it; before, Ru pulled toward zero until an input target was set.
- An MPC without a disturbance spec, on the measured state (a `MeasuredDisturbance` too) or through a filter, takes its target from `steady_state_target`'s design on the held outputs, as the LQ laws do. With as many held outputs as inputs a steady state that is not unique, or with fewer one that cannot hold every r, is refused (`singular`) where the nearest equilibrium was taken.
- `InputLimits{.lower, .upper}` is the one option for the limits on u and `InputRateLimits{.lower, .upper}` for those on Δu, among the trailing options of `design::state_feedback` and `design::mpc`, and they replace `LQRLimits` and `MPCConstraints`' `u_lim` and `du_lim`: the MPC's positional constraints argument is gone, `design::mpc<NP, NC>(plant, Ts, weights, InputLimits{.lower = {…}, .upper = {…}}, InputRateLimits{.lower = {…}, .upper = {…}}, …)`, its inputs and moves unbounded when they are left out.
- An MPC's QP carries rows on Δu only with an `InputRateLimits`: without one the QP, and the law that holds it, is 2·NC·NU rows smaller, the commands the same and `MPCLimits::input_rate` empty. A design that bounds Δu is the QP and the law it was. The default solver of an `MPC<…>` spelled out is sized for a QP without them.
- The limit options take flat `.lower` and `.upper` ends, each unbounded when left out, in place of a nested box (a breaking change for `DemandLimit`, whose `.lim` 0.3 shipped): `InputLimits{.lower = {-12.0}, .upper = {12.0}}`, `OutputLimits{.outputs = {1, 2}, .lower = {…}, .upper = {…}}`, `StateLimits{.states = {2}, .lower = {…}, .upper = {…}}` and `DemandLimit{.state = 2, .outer = {0, 1, 4}, .lower = -0.26, .upper = 0.26}`, the sizes deducing as before. `IntegralAction`'s box is `.limits`.
- An MPC's output limits, terminal ingredient and QP budget are trailing options, as its state limits are. `OutputLimits{.outputs = {1, 2}, .lower = {…}, .upper = {…}}` (with `.lower_ecr`, `.upper_ecr`) bounds only the outputs it lists, and the QP carries rows for those alone, so a design without one has no output rows; `DareCost{}` or `InvariantSet<64>{.target_margin = 0.05}`, at most one, replaces `terminal` and `terminal_target_margin`, and `MPCTerminal` is internal; `QPIterations{n}` replaces `max_qp_iterations`. The runtime's `output` limits are on the listed outputs.
- `MPCWeights` loses `u_target`: an input target is an operating point the runtime holds, `set_input_target(u)`.
- An MPC's limits on states it does not track are an option of their own, `StateLimits{.states = {1}, .lower = {-0.3}, .upper = {0.3}}`, which takes its count from its list, in place of `MPCConstraints`' `x_states`, `x_lim`, `x_lower_ecr` and `x_upper_ecr` (`.lower_ecr`, `.upper_ecr`). The runtime's `limits()` and `set_limits()` keep them, as `state`.
- A `DisturbanceModel` that leaves `.Cd` out names no sizes: `DisturbanceModel{.Bd = sys.B}` is `DisturbanceModel<NX, 0, ND>`, NY 0 meaning d reaches no reading, and the design reads C_d as zero on the plant's outputs.
- The disturbance specs take their count of held outputs from their list, as `TrackedOutputs` does: `InputDisturbance{{0, 2}}` is `InputDisturbance<2>`, and `DisturbanceModel{.Bd = Bd, .Cd = Cd, .outputs = {0}}` deduces every size. Their `outputs` is a C array; with none listed the spec has no `outputs` member.
- The library, the tests and the examples write `StateSpace{…}` with no template arguments wherever the initializer gives every size, 298 places; each was checked to deduce exactly the type it named before.
- `design::state_feedback(plant, [Ts,] weights, opts...)` replaces `lqr`, `lqi`, `lqg`, `lqgi` and `oflqg`: its options pick the law, an `IntegralAction` for LQI, a `KalmanNoise` for LQG, both for LQGI, a noise and a disturbance spec for OFLQG, `TrackedOutputs` for a tracking law. `design::mpc` takes the noise and the disturbance of what was `ofmpc`, and `design::kalman` its noise, the same way. Every design's features trail as a pack recognized by type, in any order, and an `InputDelay` moves there from before the weights: `state_feedback(plant, Ts, weights, KalmanNoise{.process = Qc, .measurement = 0.01}, InputDelay<3>{θ})`. A `KalmanNoise` deduces its sizes from its fields, a scalar for one channel. Integral action carries its own numbers, so the weights are the same with or without it: `IntegralAction{.Q = 80.0}` holds the integrators' cost, and `.outputs`, `.limits` and `.tracking_time` what `LQIConfig`'s `Q_integral`, `integral_lim` and `tracking_time` held; `LQIConfig` is internal. The cost is `LQRWeights{.Q, .R, .N}`, which `LQRConfig` was without its limits: those are options, as the gain does not depend on them, `InputLimits{.lower = {-12.0}, .upper = {12.0}}` and a `DemandLimit`. A `DemandLimit` lists its outer terms by index, the plant states then the integrators, `DemandLimit{.state = 2, .outer = {0, 1, 4}, .lower = …, .upper = …}`, and takes their count from its list; it replaces `.demand_lim`'s one flag per term and `IntegralAction`'s `demand_outer`, and a runtime takes one at run time through `set_limits()`, `l.demand = DemandLimit{…}` (`= {}` turns it off). `TrackedOutputs{{0, 2}}` takes its count from its list too. 44 design overloads become 8, and every LQ design returns a `StateFeedbackResult`; `matlab::` keeps the MATLAB® names.
- `state_feedback_core.hpp`, the runtime the LQ laws share and all of it `detail::`, moves to `controllers/detail/`.
- The header a design is in deploys its law: `state_feedback.hpp` includes the LQR, LQI, LQG, LQGI and OFLQG runtimes, and `ofmpc.hpp` folds into `mpc.hpp`, so `Law<kRes, T>` of an LQI or an OFMPC no longer needs a second header. `lqr.hpp` and the other runtime headers no longer bring the design.
- `Law<kRes, T>` is declared in `controllers/law.hpp`, which every header that specializes it exports, and the rest of `law_detail.hpp`, which `law.hpp` exports, moves to `controllers/detail/`. `lqg.hpp`, `lqgi.hpp` and `oflqg.hpp` export their law, configs and `KalmanNoise` (and `oflqg.hpp` the disturbance specs) as `lqr.hpp` and `lqi.hpp` export theirs; `mpc.hpp` exports `KalmanNoise` and the disturbance specs. `control.hpp` exports `backend.hpp` (`damp::clamp` and the other aliases), `math/math.hpp` (`damp::abs`, `damp::sqrt`, …) and the public headers it reached only through others (`complex.hpp`, `steady_state_target.hpp`, `invariant_set.hpp`, `model_result.hpp`, `dynamics.hpp`, `ins_mechanization.hpp`), the parts of `matrix.hpp` mark themselves private to it, and so does every `detail/` header (`law_detail.hpp` to `law.hpp`, `state_feedback_core.hpp` and `design/detail/config.hpp` among them) and the other internal ones (`constexpr_math.hpp`, `math_backend.hpp` and `trig.hpp` to `math.hpp`, `matrix_traits.hpp` to `matrix.hpp`, `joseph.hpp` to `eskf.hpp`), each exported by the header it names; `control.hpp` exports `split_real.hpp`. clang-include-cleaner no longer adds `law_detail.hpp`, `detail/lq_law.hpp`, `backend.hpp` or `matrix/core.hpp` to a file that includes `control.hpp`, or a part of `matrix.hpp` to one that includes a single law's header, and the examples include `control.hpp` (`workbench.hpp` in the three SILs that simulate) in place of the 2 to 20 headers each named.
- An MPC or OFMPC takes its reference and its input target braced, one value per channel (`control({0.3, 0.0})`, `set_input_target({0.1, 0.2})`), or as the scalar of a single channel. `set_input_target` takes an `InputTarget`, held or per move, so a braced `{x}` no longer reads as both.
- An OFMPC's filter is its public `estimator` member, as an OFLQG's is, no longer an `estimator()` accessor. `FixedGainObserver::reset(x0, d0)` restarts the estimate at x0 on the law's model and d0 on the disturbance, absolute as `state_estimate()` and `disturbance_estimate()` report them; before, an offset-free filter's `reset(x0)` took the deviation from its anchor.
- `LQRConfig`, `LQIConfig` and `MPCWeights` take their sizes from the fields a named one sets, with one input, integrator and controlled output unless a field says otherwise: `LQIConfig{.Q = diag(10.0, 1.0), .R = {1.0}, .Q_integral = {80.0}}` is `LQIConfig<2, 1, 1>`. A single weight is braced; the positional guides that kept `{.Q, .R, .demand_lim}` from deducing are gone.
- A Kalman design (`kalman`, `lqg`, `lqgi`, `oflqg`, `mpc` and `ofmpc` with noise, `mhe`) reads a plant that declares no noise inputs (`StateSpace<NX, NU, NY>`) as noise on every state and reading, G = I and H = I, so its `KalmanNoise` is NX × NX and NY × NY and the plant needs no `<…, T, NW, NV>`. Before, such a plant took only a `KalmanNoise<0, 0>`.
- `PRSpec::Tt`, `PRResult::Tt` and `PRController::Tt` are `tracking_time`, the name the LQI and LQGI runtimes give the same back-calculation tracking time.
- Every example SIL checks what it prints and returns nonzero when a result leaves its band, so `make` fails on a change in behavior, not only on a crash. Before, 2 of 18 did.
- `StateSpace::is_discrete()` and `is_continuous()` are both false for a `Ts` that is NaN, ±Inf or negative: such a model has no time domain, and every design, conversion and analysis that reads one refuses it (`invalid_sample_time`, or an empty result). Under `-ffast-math` a NaN `Ts` read as continuous.

### Fixed

- A `DemandLimit` whose box has an infinite or NaN end is refused (`invalid_parameter`), by the design and by `set_limits()`: its conversion to the law's scalar read an infinite end as unbounded before the box was tested.
- A `constexpr` design refused for its sample time now says why plain braces cause it. `fn(plant, {}, {…})` is the sampled form with Ts = 0, for a discrete plant as for a continuous one; the compile error's note at `design_refused::invalid_sample_time` names both right spellings, `fn(plant, weights, …)` for a discrete plant and `fn(plant, Ts, weights, …)` for a continuous one.
- Inside a template, clang deduces nothing from a designated field set to a literal or a constant ([llvm/llvm-project#114870](https://github.com/llvm/llvm-project/issues/114870)), so a test's `IntegralAction{.Q = 1e6}` in a function template did not compile there: the tests name the type, `IntegralAction<double>{.Q = 1e6}`, the design's message says to, and the [known limitations](docs/known_limitations.md) list the options it reaches. `make clang-tests` builds and runs the three test runners with clang as well, and the tests and examples build on clang without a warning, a `ColVec` written flat, `ColVec<2>{a, b}`, where `ColVec<2>{{a}, {b}}` drew clang's `-Wbraced-scalar-init`.
- Under `-ffast-math` a NaN passed the guards meant to refuse it, since g++ folds `!(x > 0)` into `x <= 0`, which a NaN fails: a cutoff, Q, sample time, gain, limit or config field in the filters, the controller and observer runtimes, the estimators, the PID and MPC helpers, the toolbox, the motor helpers and the host simulation. Each guard now tests finiteness by the bit pattern and refuses a NaN, or a +Inf where the value must be finite, as it documents; a PWM duty map parks mid-rail on a NaN bus voltage, and `analysis::arange` with a NaN step returns empty instead of looping forever. Inputs that had no check get one: `SmithPredictor`'s model pole, `SogiFll`'s FLL gain, the disk-margin skew, the adaptive ODE solver's options, the discrete PID runtimes' gains and others.
- An MPC or OFMPC built at run time from a design whose filter or target was refused reports `valid()` false and holds its input, as one with a refused QP does; it ran on the refused target before.
- The FOC SIL reported a 98% rise time one tick short (0 ms for the PI), and the LPF SIL called its first-order filter second-order Butterworth. The workflow SIL ran 2 s, before its 1 Hz PR mode settled; it runs 20 s and ends at r.

## [0.3.1] — 2026-11-17

The LQ integrators hold against a clamp by default again, as a PID's integral does.

### Changed

- The LQI and LQGI integrators hold against a clamp (conditional integration) when `tracking_time` is 0, the default, as a PID's integral does with `Kbc = 0`. A `tracking_time` above 0 unwinds them by back-calculation, and one up to Ts takes the whole excess each tick, 0.3.0's default: set `tracking_time = Ts` to keep 0.3.0's behaviour. Near its limit with noise on the reading, that default took each excursion that reached the clamp back from the integrators, from one side only, and an LQI settled short of r by about ten times the offset of holding; a slow integrator gave back for seconds what a short clamp took from it. Back-calculation still recovers from a long clamp without holding's overshoot.
- An LQ design whose (∂u/∂ξ)⁺ cannot be formed is refused only when it asks for back-calculation.

## [0.3.0] — 2026-11-17

Designs for a plant whose input arrives late, one limits and anti-windup contract for every law (with state limits on MPC and a bound on the offset-free laws' disturbance estimate), a Kalman filter and an MPC that hold in float, no fault at run time, and Damp after `Arduino.h`. The API breaks: `LQRWeights` is now `LQRConfig`, or `LQIConfig` with integral action, and a `SisoController` needs `back_calculate`.

### Added

- `MPCConstraints` takes `x_states` and `x_lim`, with their softness: an MPC or OFMPC bounds states it does not track over the horizon, soft or hard, and `set_limits` moves them. The type names how many, `MPCConstraints<NU, NY, T, NXC>`.
- `last_status()` on `KalmanFilter` and `ExtendedKalmanFilter`: whether the last measurement update was taken or refused (`UpdateStatus`), so the refusal `step()` drops is still there to read. MPC's `last_status()` reports its QP the same way.
- `OFMPC::disturbance_lim` clamps the disturbance estimate d̂ at each update, as on `OFLQG`, and `estimator()` has a mutable overload to restart the filter or set its estimate.
- `to_ss()` on an OFLQG design: the law as a discrete block at r = 0, its filter on [x; d] and the target's d̂ terms folded into the gain, to close on a plant other than the design's (one with another delay, say) for its poles.
- `LQRConfig` and `LQIConfig` take `demand_lim`, a `DemandLimit`: the LQ laws bound what the gain's outer terms ask of one state, read as a cascade (a balancing robot's position and integral asking for a lean), and the integrators unwind against it. Off by default, and free in a `constexpr` design's `Law`.
- `OFLQG::disturbance_lim` clamps the disturbance estimate d̂ at each update, so a load past it (a wall the plant is held against) does not wind d̂ and the target up, the way `integral_lim` bounds an LQI's integrator. Unbounded by default.
- `design::lqr`, `lqi`, `lqg`, `lqgi`, `oflqg`, `kalman` and `steady_state_target` take an `InputDelay` after Ts: the continuous plant, its weights and its noise, the cost and the noise discretized exactly on the delayed model. A measured-state law holds the registers itself, so `update` takes the plant state; a filter carries them with no noise on them.
- `LQRConfig`, the weights of an LQ design, takes `u_lim`, and `LQIConfig` also `integral_lim`: the LQR, LQI, LQG, LQGI and OFLQG runtimes clamp their command and their integral state to them, and unwind the integrators against that clamp with no `back_calculate`. Both are members of the runtime (`law.u_lim`, `law.integral.lim`), so a limit can move at run time.
- `Bounds::valid()`: no end NaN and no lower end above its upper one, the test `design::mpc` and `design::mhe` already applied.
- `u_lim` on `ADRCResult`, `STSMCResult` and `LeadLagResult`, and as a member of their runtimes and of `ExtremumSeekingController`: the command is clamped into it and the law's memory holds against that clamp.
- `back_calculate(u_sat, u_ff)` on `STSMCController`, `LeadLagController`, `ExtremumSeekingController`, `SmithPredictor`, `RepetitiveController` and `ClassicalDOB`: a clamp after the law is reported the same way to every law with memory.
- `back_calculate` on every SISO law and in the `SisoController` concept: P, PD and SMC take it as a no-op, and `HarmonicSuppressor`, `Cascade` and `SwitchedController` pass a clamp on to the laws they hold.
- `u_lim` on `PRResult`, `HarmonicSuppressorResult` and `SMCResult` and as a member of their runtimes: the command is clamped into it, and the resonators of a PR law or a harmonic bank hold against the clamp.
- `valid()` on the PID runtimes, the LQ laws, `MPC` and `OFMPC`, as on the other laws: false when built from a design they refuse, and the law sends 0 (an MPC holds its input).

### Changed

- `StateFeedbackResult` and the LQR, LQI, LQG, LQGI and OFLQG runtimes take a last template argument, NRU, the input registers of the law's model (0: none); K, S and e are on [x; v(k−1); …; v(k−NRU)] then the integrators.
- `damp::min`, `damp::max` and `damp::clamp` return by value on both backends, so a result bound to a reference cannot dangle on a temporary argument. The backends now agree: a tie gives `min` and `max` their first argument (a NaN too, without `-ffast-math`), and `clamp` with lo > hi gives hi, where `std::clamp` is undefined.
- `LQRWeights<NX, NU, NI>` is `LQRConfig<NX, NU>` for `lqr`, `lqg` and `oflqg`, and `LQIConfig<NX, NU, NI>` with `Q_integral` for `lqi` and `lqgi`: each has only its law's fields, the limits with the cost. A design called with braces, `{.Q = Q, .R = R}`, does not change.
- `mat::pseudo_inverse`, `mat::null_space` and `mat::two_norm` return an optional, empty when the SVD fails, where they returned NaN; `matlab::pinv` and `matlab::null` follow. `design::steady_state_target` with more outputs than inputs refuses a failed factorization, where it could report NaN maps as a success.
- `+` and `−` on matrices and views need equal shapes: a dimension of size 1 is no longer stretched to fit (NumPy's broadcast), so a missing transpose does not compile. An offset on every entry is `A + scalar`.
- A PID with an output limit and `Kbc = 0` holds its integral against the clamp (conditional integration), where it wound up, and `back_calculate` does the same for a clamp downstream. A `Kbc` still unwinds it by back-calculation.
- `PIDGains::discretize`, the PID runtimes' constructors and `design::esc` refuse a limit that is not a box, as the LQ and MPC designs do. In `ESCSpec`, `lower ≥ upper` no longer switches the clamp off.
- `ADRCController::last_command()` is the command the law returned, as on the other laws; the input its observer takes as applied is `applied()`.
- `Bounds::as` saturates a finite end past the narrower type's range.
- `state_feedback_core.hpp` holds only the runtime the LQ laws share. `LQRConfig`, `IntegralAction`, `NoEstimator` and `StateFeedbackResult` are in `state_feedback_types.hpp`, and the design behind them in `controllers/detail/`; `lqr.hpp`, `lqi.hpp`, `lqg.hpp`, `lqgi.hpp` and `oflqg.hpp` include all of it as before.
- `solve_qp` at its iteration budget: a solve that its last allowed iteration completes is `Success`, and a `MaxIterations` point is the minimizer on the constraints taken in so far.
- A float MPC riding a hard output limit at full input can read `Infeasible` on a tick whose rounded state is past the limit's reach, and hold its input. It returned `Success` there with the limit's row violated.
- `KalmanFilter` and `ExtendedKalmanFilter` carry the triangular factor S of their covariance, P = S Sᵀ: `covariance()` returns P by value, formed on each call, and `covariance_factor()` is S.
- A perfect measurement on `KalmanFilter` or `ExtendedKalmanFilter` (R = 0, or a zero-variance row of R) is taken at the rounding of its own prior variance: the estimate reads it to eps, and the covariance keeps eps times the prior there instead of zero.
- A P₀ or a Q of `KalmanFilter` or `ExtendedKalmanFilter` may be singular, and one that is not positive semidefinite or not finite does not compile where the filter is constant-evaluated; at run time P₀ becomes the identity and such a Q adds no process noise.
- The update of `KalmanFilter` and `ExtendedKalmanFilter` takes far less stack, about a ninth for a large filter; a `KalmanFilter` that holds its matrices is larger by the factors of its Q and R.
- Damp never faults at run time. A PID runtime built from gains it refuses is invalid where its constructor trapped, `DiscretizationResult::value()` of a refused result is the zero model, and a matrix initializer longer than the matrix drops the extra entries; each is still a compile error under constant evaluation.
- `Cascade` holds the outer law while the inner law's command is clamped, at its `u_lim` or by a clamp reported with `back_calculate`, so an outer integrator no longer winds up while the inner loop cannot follow.
- `SwitchedController` starts a law with a tracking mode, a PID, from the command the plant last received when it is switched in, where it started from its reset state.
- `PRController` with `Tt = 0` holds its resonant section against a clamp reported with `back_calculate`, where it ignored the report.
- The integrators' limits of the LQI and LQGI runtimes and their design results are `integral_lim`, the name `LQIConfig` gives them, where they were `integral.lim`.
- The LQI and LQGI integrators unwind against a clamp by back-calculation, where they held (conditional integration): a law recovers from a long clamp without the overshoot. `tracking_time` on `LQIConfig` and the runtimes slows it; the default takes the whole excess back each tick.

### Fixed

- Damp's headers compile after `Arduino.h`, the order an `.ino` sketch gets: Damp undefines Arduino's `abs` macro, as `<cstdlib>` does, and no longer names anything its `F`, `PI` or binary-constant (`B1`) macros rewrite.
- A float MPC with soft output limits solves its QP on every tick. Once the moves sat on their limits the active-set solver ended `Degenerate`, or called a feasible problem `Infeasible`, on up to every tick of a recovery from outside the limits; it now keeps an orthogonal factorization of the active constraints.
- `solve_qp` in float no longer returns `Success` with a constraint violated past its tolerance: the point is recomputed on the active constraints after each step, and every solved point meets its rows.
- PID back-calculation with a `Kbc` under half of `Ki·Ts` diverged, and overshot under `Ki·Ts`; it takes at most the whole excess per tick.
- An LQ design refused for its limits no longer carries them into a runtime, and an `integral_lim` narrowed at run time pulls the integrator in on the next tick.
- `PRController<float>` and `HarmonicSuppressor<float>` hold a resonance far below the sample rate, where the output used to lose its gain at ω₀ or run to Inf (f₀/fs = 5e-5). The resonant section, `PRController::resonant`, is a `PRResonator` on the output and its increment instead of a `Biquad`.
- `KalmanFilter` and `ExtendedKalmanFilter` in float keep their covariance under a prior far wider than a reading's noise: stored as P it went indefinite within a few ticks once P₀/R passed about 10⁵, and nearly every later update was refused.
- `PRController` with a tracking time `Tt` under the sample period overshot the realizable command on each unwind, and diverged under half of it; it takes at most the whole excess per tick.

## [0.2.2] — 2026-11-10

A fix for graded matrices, which 0.2.1's balancing change broke.

### Fixed

- Balancing takes a matrix's grading out in full again, and keeps it only where it pays, so `logm` of a graded matrix in float has the accuracy of the ungraded one (0.2.1 lost up to 16×, and failed a test on some hosts). A matrix of one scale is still left unscaled.

## [0.2.1] — 2026-11-10

Numerical robustness fixes across the matrix core, tested in float first, and Damp on the PlatformIO registry.

### Added

- `library.json`: Damp on the PlatformIO registry, `lib_deps = wetmelon/damp@^0.2.1`.

### Changed

- `mat::ldlt` pivots, P·A·Pᵀ = L·D·Lᴴ with 1×1 and 2×2 blocks, so `LDLT` gains `perm` and `e`. `ldlt_solve(A, B)` is unchanged.
- `QRDecomposition::full_rank()` is removed: it was not a rank test. Decide rank with `mat::rank`.
- An eigenvalue, Schur entry, generalized eigenvalue or singular value that overflows the type is a failure (`success = false`), not an Inf.
- For a matrix with isolated eigenvalues, `EigenResult::schur_vectors` includes the balancing permutation, and `values` can come in a different order.

### Fixed

- Balancing counts the diagonal in its row and column norms, so a matrix whose diagonal dominates is no longer over-scaled. This fixes `sqrtm_jordan_3`, which failed `make` on Intel and generic-tuned hosts. `balance_states` does the same.
- Eigenvalues and Schur forms balance with a permutation step: triangular, cascade and integrator-augmented matrices get exact eigenvalues, where a 4×4 with all-negative eigenvalues could return a positive one.
- Real Schur forms, `lyap`, `sqrt` and `log` hold from 1e-30 to 1e30 in float and from 1e-300 to 1e300 in double. A complex pair used to come back wrong past about 1e±19 in float and 1e±155 in double.
- A lower-triangular 2×2 Schur block is exchanged exactly instead of losing its small eigenvalue to cancellation.
- The QR sweep's reflector stays symmetric under `-ffast-math`, so the Schur form keeps the accuracy of a build without it.
- QZ (`generalized_eigenvalues`, `transmission_zeros`) solves its 2×2 blocks without cancellation: close eigenvalues keep their digits, and a real pair no longer comes out complex.
- QZ scales the pencil by powers of two, so ‖A‖/‖B‖ no longer has to stay within about 1e±18 in float.
- `svd` prescales by a power of two, so a matrix near the ends of the type's range no longer returns all-zero singular values as a success.
- `cholesky`, LU and the triangular solves refuse a NaN or Inf pivot under `-ffast-math`. `KalmanFilter::update` with a NaN in R returns false and leaves the state untouched.
- `ldlt` stays accurate on indefinite matrices (a float case went from 3.6e-3 relative error to 4.7e-9) and refuses NaN under `-ffast-math`.
- `state_feedforward` and the decoupling in `output_inverse` decide rank by singular values, so a rank-deficient B is refused whatever its column order.
- Norms pass a NaN or Inf entry through, so a `residual.norm() <= tol` check no longer passes on a NaN residual.

## [0.2.0] — 2026-11-09

One runtime and one set of tick verbs for every LQ and MPC design, input and output delay in the plant model, and a core that builds with no C++ standard headers. The API breaks from 0.1.0 throughout: see Changed and Removed.

### Added

- `design::oflqg` / `OFLQG` and `design::ofmpc` / `OFMPC`, output-feedback LQG and MPC, over one shared state-feedback core.
- `TrackedOutputs`: LQR and LQG can `control(r)` / `step(r, m)` through the model's steady-state target.
- `design::output_inverse` / `OutputInverse` and `design::state_feedforward`.
- `InputDelay<N>{θ}` and `OutputDelay<N>{θ}` on `discretize`: input and measurement delay in N registers, exact under ZOH, one delay per channel if needed. Tustin pre-warp: `discretize(sys, Ts, Tustin, w0)`.
- `Law<kRes, T>`: a `constexpr` design runs on its constants, so its matrices take no RAM and its exact zeros cost no code. It covers the LQ family, the observers, ESO, ADRC, `StateFeedforward`, `OutputInverse`, OFMPC's filter, and `design::kalman` through `Law<kRes, T>::Recursive`.
- `SparsityPattern` and `PatternMatrix`: zeros and ones fixed at compile time, other entries at run time, used by the ESKF, the INS and `KalmanFilter::update`. A 15-state `F P Fᵀ` on a Cortex-M4F drops from 37,618 to 3,906 estimated cycles.
- `DiscretizationResult`, `lft` and canonical forms; transmission zeros (QZ), `norm_hinf` and disk margins on the host.
- `diag` / `blkdiag` / `hcat` / `vcat`, `mat::ldlt`, `lyapunov_cholesky` / `dlyapchol`, an SOS IIR cascade, `sim::quantize`, and encoders under `damp/encoder/`.
- `design::eso` / `ESO`, the Han–Gao extended-state observer on its own, and `reset(z0)` on `ESO` and `ADRCController`.
- `FeedbackLaw`: `control(r)` and `step(r, m)` on top of the `FeedbackRuntime` lifecycle.
- Examples `heater_ofmpc` and `heater_lqg`; [docs/guide.md](docs/guide.md), from plant to law with MATLAB®, firmware and python-control names; [docs/known_issues.md](docs/known_issues.md), compiler bugs that change results.
- `make workbench-check`, and the README and guide snippets compiled as checks. The ETL backend takes any copy of ETL on the include path. CMake alias `wetmelon::damp`.

### Changed

- One tick vocabulary: estimators `predict` / `update`, controllers `control`, `step` for a whole tick, and `back_calculate(u_sat, u_ff = {})` for the applied input. Ts lives on the design.
- Every design reports `success` / `DesignStatus` and `explicit operator bool()`, a refused `constinit` law names the cause, and every `*Result` struct is `[[nodiscard]]`.
- `Luenberger` → `FixedGainObserver`, `ReducedLuenberger` → `ReducedOrderObserver`, `design::luenberger` → `design::place_observer`, with their results and header renamed to match.
- `mat::EigenResult::values` is a `damp::array<damp::complex<T>, N>`, like the other eigenvalue lists.
- `LQRWeights` takes NX, NU, NI, T (`LQRWeights<2, 1, 1>` is LQI); `design::lqr` follows the plant's time domain; `PIDGains` is the PID result; `design::kalman` takes a disturbance spec.
- MPC and MHE are allocation-free in `control.hpp`. An infeasible MPC tick holds the last applied input, and a horizon too short for an input to reach any output is refused unless a terminal cost or an `Ru` weight covers it.
- The core builds on ETL with no C++ standard headers (`ETL_NO_STL`), checked by `make nostd-check`. With no profile and no backend macro it builds on the host defaults without a warning.
- The examples and the guide deploy as `Law<kRes, float>`, and the README and the guide show the three deploy forms side by side.
- Matrix operations skip the zero-fill of results they write in full, taking `memset` and `memcpy` off the run-time path: a dense 3-row ESKF update drops from 61,809 to 51,463 estimated cycles on a Cortex-M4F.
- Compilers: GCC 11.3+ and Clang 19+; Clang on Windows targets MinGW's libstdc++. MSVC is not supported.
- Build: examples and `make target-smoke` link statically, `make` no longer needs `arm-none-eabi-g++`, and `make embedded-check` checks only Damp's own headers.

### Removed

- `CONFIG_BACKEND` and `CONFIG_COMPILER_PREFIX` in `tup.config`.
- `ins_propagate_covariance` and `ins_F_mul`: `ErrorStateKalmanFilter::predict(F)` on the sparse transition replaces them.
- `design/synthesis.hpp` (the `*_bundle` wrappers), `StateFeedback` as a public runtime, and `offset_free_mpc.hpp`, which `OFMPC` replaces.
- `toolbox/encoder.hpp` (now `damp/encoder/`), `toolbox/thermal.hpp`, `estimation/commissioning_frontends.hpp`, `estimation/experiment_safety.hpp`, `estimation/successive_compensator.hpp` and `filters/highpass.hpp`.

### Fixed

- GCC 11 and 14 no longer crash (internal compiler error) on a `constexpr` constructor that copies a `Matrix`.
- `UnscentedKalmanFilter::NSIG` is now `NSP`, since MinGW's `<signal.h>` defines `NSIG`.
- Results no longer depend on a cancellation rounding to exactly 0 under GCC 15 `-ffast-math` or fused multiply-adds: `design::eso`, the implicit integrators, the transfer-function stability checks and `MHE`.
- `lyap` and `dlyap` return an exactly symmetric X for a symmetric Q.
- `design::place_jordan` no longer fails under GCC 15.2's vectorizer, and GCC 15's warnings are gone.
- The PLLs treat a subnormal amplitude as no signal, so `integrator_leak` bleeds the offset back.
- `assoc_barrier` and `discretize` with a NaN `Ts` work on GCC 11, `design::lowpass<3>` / `<4>` keep unit DC gain in float, `math_backend.hpp` honours the profile when included first, and `dpwm_zero_sequence` is constexpr on Clang 19.
- Build: `make embedded-check` reports compile errors, `make tidy` / `iwyu` / `fix_includes` work with clang tools, CMake's `docs` target exists only at the top level, and the README's FetchContent snippet skips the submodules.

## [0.1.0] — 2026-08-09

First public release of Damp: header-only C++20 control design that ships in the
same tree as firmware (design-is-deploy: same types in SIL and on the MCU).

License: Boost Software License 1.0. Third-party host-only components are listed
in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

### Entry points

| Include | Role |
| ------- | ---- |
| `damp/control.hpp` | Embeddable core: matrix, systems, `design::`, core design-is-deploy controllers/estimators (incl. ESKF), filters, fixed-step integrators, toolbox |
| `damp/workbench.hpp` | Host: control plus Bode/margins, ODE solvers, closed-loop `simulate`, `matlab::` short names |

Drop `inc/` into an existing firmware tree, or use CMake (`find_package` / FetchContent).

### Library

- Matrix / math: fixed-size stack `Matrix`, views, decompositions, solve; pluggable math backend; geometry (quaternion / DCM)
- LTI systems: `StateSpace`, transfer functions, ZPK, discretization (ZOH / Tustin / Euler), interconnections
- Design (not PWM-rate): place, Riccati / DARE, PID rules, minreal / model reduction, stability, synthesis bundles → Result + `.as<U>()`
- Core design-is-deploy controllers: PID, PR, LQR / LQI / LQG / LQGI, lead-lag, ADRC, SMC / STSMC, Smith predictor, action governor
- Core design-is-deploy estimators: Kalman, EKF, UKF, ESKF + attitude fusion, Luenberger, DOB, RLS
- Filters: biquad family and common signal blocks; robust exact differentiator
- Simulation: fixed-step integrators on-target; host ODE solvers and closed-loop simulate via workbench
- Host analysis: Bode, Nyquist, margins, step/impulse/lsim, poles; MATLAB®-style aliases

Runtime controllers and estimators default to `float`; design synthesis defaults to `double`.

### Examples and validation

- Design-is-deploy folder layout under `examples/control/` and `examples/estimation/`
- Unit tests and embedded/freestanding contracts for the core surface

### Documentation

- [README.md](README.md)
- [REFERENCE.md](REFERENCE.md)
- [docs/known_limitations.md](docs/known_limitations.md)
- [examples/](examples/)

### Known limitations

See [docs/known_limitations.md](docs/known_limitations.md).

[0.3.1]: https://github.com/Wetmelon/damp/releases/tag/v0.3.1
[0.3.0]: https://github.com/Wetmelon/damp/releases/tag/v0.3.0
[0.2.2]: https://github.com/Wetmelon/damp/releases/tag/v0.2.2
[0.2.1]: https://github.com/Wetmelon/damp/releases/tag/v0.2.1
[0.2.0]: https://github.com/Wetmelon/damp/releases/tag/v0.2.0
[0.1.0]: https://github.com/Wetmelon/damp/releases/tag/v0.1.0
