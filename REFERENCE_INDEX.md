# API Reference — Alphabetical Index

Auto-generated from `@brief` doc comments in `inc/damp/`. Regenerate with `python tools/gen_reference.py`. Grouped-by-domain view: [REFERENCE.md](REFERENCE.md).

| Name | Kind | Domain | Description |
| ---- | ---- | ------ | ----------- |
| [`abs`](inc/damp/math/complex.hpp#L362) | function | Scalar math, complex & frames | Compute magnitude (absolute value) of a complex number (+1 more overload) |
| [`acker`](inc/damp/matlab.hpp#L861) | function | MATLAB®-style aliases (host) | Single-input pole placement (MATLAB®'s acker) |
| [`acos`](inc/damp/math/math.hpp#L149) | function | Scalar math, complex & frames | Arccosine ∈ [0, π]. Input is clamped to [−1, 1] in both paths |
| [`ActionGovernor`](inc/damp/controllers/action_governor.hpp#L235) | block | Runtime controllers | Runtime action governor — QP projection of u_des onto A u ≤ b |
| [`ActiveSetSolver`](inc/damp/design/qp.hpp#L966) | block | Design-time synthesis (not PWM-rate) | Default QP solver policy: the Goldfarb–Idnani active-set solve_qp() |
| [`AdaptationGate`](inc/damp/estimation/parameter_estimation.hpp#L338) | block | Observers & estimators | Gated-adaptation policy: decide whether a parameter update is safe to apply |
| [`adaptive_solve`](inc/damp/simulation/solver.hpp#L828) | function | Simulation / SIL harness (host) | Adaptive-step solve in one call |
| [`AdaptiveLut1D`](inc/damp/toolbox/lookup.hpp#L427) | block | Embedded helpers (controls-adjacent utilities) | Adaptive 1-D lookup: fixed breakpoints, online-updated cell values |
| [`AdaptiveLut2D`](inc/damp/toolbox/lookup.hpp#L500) | block | Embedded helpers (controls-adjacent utilities) | Adaptive 2-D lookup: fixed grid, online-updated cells (nearest grid node) |
| [`AdaptiveOptions`](inc/damp/simulation/solver.hpp#L349) | block | Simulation / SIL harness (host) | Options for adaptive-step ODE integration |
| [`AdaptiveStepIntegrator`](inc/damp/simulation/integrator.hpp#L133) | concept | Simulation / SIL harness (host) | Integrator that reports a genuine embedded local-error estimate |
| [`AdaptiveStepSolver`](inc/damp/simulation/solver.hpp#L404) | block | Simulation / SIL harness (host) | Adaptive-step ODE solver (+1 more overload) |
| [`AdaptOutOfRange`](inc/damp/toolbox/lookup.hpp#L47) | enum | Embedded helpers (controls-adjacent utilities) | How AdaptiveLut1D / AdaptiveLut2D treat samples outside the breakpoint span |
| [`AdmmSettings`](inc/damp/design/qp.hpp#L1067) | block | Design-time synthesis (not PWM-rate) | Tuning parameters for the ADMM QP solver |
| [`AdmmSolver`](inc/damp/design/qp.hpp#L1148) | block | Design-time synthesis (not PWM-rate) | ADMM (OSQP-style) QP solver policy — warm-started, fixed-cost iterations |
| [`adrc`](inc/damp/controllers/adrc.hpp#L207) | function | Design-time synthesis (not PWM-rate) | Active Disturbance Rejection Control design |
| [`ADRCController`](inc/damp/controllers/adrc.hpp#L413) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its design, keeping the result's scalar |
| [`ADRCController`](inc/damp/controllers/adrc.hpp#L258) | block | Runtime controllers | Active Disturbance Rejection Control (ADRC) |
| [`ADRCResult`](inc/damp/controllers/adrc.hpp#L50) | block | Design-time synthesis (not PWM-rate) | Active Disturbance Rejection Control design result |
| [`AffineCal`](inc/damp/toolbox/scaling.hpp#L81) | block | Embedded helpers (controls-adjacent utilities) | Affine sensor calibration `y = gain·x + offset` |
| [`all_of`](inc/damp/matrix/core.hpp#L1551) | function | Linear algebra | Whether p holds for every entry of A: the matrix counterpart of std::all_of |
| [`allmargin`](inc/damp/matlab.hpp#L1326) | function | MATLAB®-style aliases (host) | Gain, phase, and delay margins of a SISO loop over a frequency grid |
| [`AllMarginResult`](inc/damp/matlab.hpp#L1305) | block | MATLAB®-style aliases (host) | All classical margins including delay margin (superset of margin) |
| [`allpass`](inc/damp/filters/iir_design.hpp#L723) | function | Filters & signal conditioning | Second-order all-pass filter design (continuous-time) (+1 more overload) |
| [`AlphaBeta`](inc/damp/math/transforms.hpp#L163) | block | Scalar math, complex & frames | Alpha-beta (stationary-frame) component pair |
| [`amigo_kappa_tau`](inc/damp/design/pid_design.hpp#L257) | function | Design-time synthesis (not PWM-rate) | AMIGO PI from ultimate gain/period and static gain Kₛ |
| [`AnalogCrossMode`](inc/damp/toolbox/conditioning.hpp#L465) | enum | Embedded helpers (controls-adjacent utilities) | How two redundant analog channels should relate in engineering units |
| [`AnalogCrossResult`](inc/damp/toolbox/conditioning.hpp#L487) | block | Embedded helpers (controls-adjacent utilities) | Output of cross_check_analog |
| [`AnalogCrossStatus`](inc/damp/toolbox/conditioning.hpp#L473) | enum | Embedded helpers (controls-adjacent utilities) | Result of comparing two redundant analog channels |
| [`AnalogInput`](inc/damp/toolbox/io.hpp#L47) | block | Embedded helpers (controls-adjacent utilities) | A single analog input: range/fault check on the raw reading, then affine calibration to engineering units |
| [`AnalyticDynamics`](inc/damp/concepts.hpp#L238) | concept | Core, configuration & backend vocabulary | Dynamics that provide their own Jacobians: f.jacobian(x, u) → DynamicsJacobian |
| [`AnalyticMeasurement`](inc/damp/concepts.hpp#L264) | concept | Core, configuration & backend vocabulary | A Measurement that provides its own Jacobians: h.jacobian(x, u) → MeasurementJacobian |
| [`any_of`](inc/damp/matrix/core.hpp#L1575) | function | Linear algebra | Whether p holds for some entry of A: the matrix counterpart of std::any_of |
| [`apply_wet_plotly_shell`](inc/damp/simulation/plot_plotly.hpp#L66) | function | Simulation / SIL harness (host) | Apply damp's full-window HTML shell to plotlypp-generated markup |
| [`arange`](inc/damp/analysis/linspace.hpp#L115) | function | Frequency-domain analysis (host) | Half-open arithmetic range [start, stop) with step step (+2 more overloads) |
| [`arg`](inc/damp/math/complex.hpp#L373) | function | Scalar math, complex & frames | Compute argument (phase angle) of a complex number |
| [`array_as`](inc/damp/backend.hpp#L419) | function | Core, configuration & backend vocabulary | Element-wise static_cast of a damp::array to scalar type U (precision rebind) |
| [`asin`](inc/damp/math/math.hpp#L130) | function | Scalar math, complex & frames | Arcsine ∈ [−π/2, π/2]. Input is clamped to [−1, 1] in both paths so behavior matches at compile and run time (std::asin would return NaN for \|x\| > 1) |
| [`atan`](inc/damp/math/math.hpp#L116) | function | Scalar math, complex & frames | Single-argument arctangent ∈ (−π/2, π/2) |
| [`atan2`](inc/damp/math/math.hpp#L104) | function | Scalar math, complex & frames | Two-argument arctangent, atan2(y, x) ∈ [−π, π] |
| [`AxisInput`](inc/damp/toolbox/io.hpp#L148) | block | Embedded helpers (controls-adjacent utilities) | Operator-axis conditioning chain (joystick / RC stick → command) |
| [`AxisPosition`](inc/damp/encoder/axis_position.hpp#L46) | block | Embedded helpers (controls-adjacent utilities) | The position of an axis, rotary or linear: an exact count carrying its own scale |
| [`backward_substitute_transpose`](inc/damp/matrix/solve.hpp#L77) | function | Linear algebra | Backward substitution for Lᴴx = b (Lᵀx = b for real T) |
| [`BackwardEuler`](inc/damp/simulation/integrator.hpp#L361) | block | Simulation / SIL harness (host) | Backward Euler integrator |
| [`balance_states`](inc/damp/systems/realization.hpp#L83) | function | LTI systems (SS / TF / ZPK / discretize) | Balance a model by a diagonal state scaling, so A, B and C have comparable magnitudes |
| [`balanced_realization`](inc/damp/design/model_reduction.hpp#L483) | function | Design-time synthesis (not PWM-rate) | Balanced realization via the square-root (Moore / Laub) method |
| [`balanced_reduction`](inc/damp/design/model_reduction.hpp#L560) | function | Design-time synthesis (not PWM-rate) | Balanced reduction to NR states, by truncation or DC-matched residualization |
| [`BalancedRealizationResult`](inc/damp/design/model_reduction.hpp#L76) | block | Design-time synthesis (not PWM-rate) | Balanced realization result (Moore square-root method) |
| [`BalancedReductionResult`](inc/damp/design/model_reduction.hpp#L108) | block | Design-time synthesis (not PWM-rate) | Reduced-order model from balanced truncation or residualization |
| [`balreal`](inc/damp/matlab.hpp#L344) | function | MATLAB®-style aliases (host) | MATLAB® short alias for design::balanced_realization |
| [`balred`](inc/damp/matlab.hpp#L363) | function | MATLAB®-style aliases (host) | MATLAB® short alias for design::balanced_reduction to NR states |
| [`bandpass`](inc/damp/filters/iir_design.hpp#L614) | function | Filters & signal conditioning | Second-order band-pass filter design (continuous-time) (+1 more overload) |
| [`bandwidth`](inc/damp/matlab.hpp#L1414) | function | MATLAB®-style aliases (host) | -3 dB bandwidth of a SISO system over a frequency grid |
| [`BDF2`](inc/damp/simulation/integrator.hpp#L420) | block | Simulation / SIL harness (host) | Backward Differentiation Formula 2 (BDF2) integrator |
| [`Biquad`](inc/damp/filters/biquad.hpp#L536) | block | Filters & signal conditioning | Second-order IIR (biquad) section runtime (+1 more overload) |
| [`BiquadCascade`](inc/damp/filters/biquad.hpp#L626) | block | Filters & signal conditioning | Cascade of second-order sections (SOS) for higher-order IIR filters (+1 more overload) |
| [`bit_cast`](inc/damp/backend.hpp#L163) | function | Core, configuration & backend vocabulary | Reinterpret the bits of from as To, usable in constant evaluation (etl::bit_cast is not, without the STL) |
| [`bit_width`](inc/damp/backend.hpp#L169) | function | Core, configuration & backend vocabulary | Bits needed to represent value, as int like std::bit_width (etl::bit_width returns T) |
| [`BLINK`](inc/damp/toolbox/iec61131.hpp#L539) | block | Embedded helpers (controls-adjacent utilities) | BLINK (free-running square-wave / flasher) |
| [`blkdiag`](inc/damp/matrix/concat.hpp#L133) | function | Linear algebra | Block-diagonal assembly of MatrixLike arguments |
| [`Block`](inc/damp/matrix/block.hpp#L37) | block | Linear algebra | Block view (non-owning) into a parent matrix |
| [`bode`](inc/damp/analysis/frequency.hpp#L267) | function | Frequency-domain analysis (host) | Compute Bode plot data for a SISO state-space system (+2 more overloads) |
| [`bodemag`](inc/damp/simulation/plot_plotly.hpp#L798) | function | Simulation / SIL harness (host) | Plot a magnitude-only Bode diagram (log frequency, dB magnitude) |
| [`bodeplot`](inc/damp/simulation/plot_plotly.hpp#L784) | function | Simulation / SIL harness (host) | Plot magnitude and phase Bode subplots |
| [`BodeResult`](inc/damp/analysis/frequency.hpp#L59) | block | Frequency-domain analysis (host) | Bode plot data for a SISO system |
| [`Bounds`](inc/damp/toolbox/bounds.hpp#L52) | block | Embedded helpers (controls-adjacent utilities) | A per-channel closed-interval box constraint |
| [`butterworth_lowpass`](inc/damp/filters/iir_design.hpp#L781) | function | Filters & signal conditioning | Butterworth low-pass filter design (+1 more overload) |
| [`Button`](inc/damp/toolbox/io.hpp#L195) | block | Embedded helpers (controls-adjacent utilities) | Debounced momentary push-button with edge and long-press detection |
| [`c2d`](inc/damp/matlab.hpp#L426) | function | MATLAB®-style aliases (host) | MATLAB® interface function c2d to discretize a continuous-time state-space system (+1 more overload) |
| [`CachedZoh`](inc/damp/simulation/cached_zoh.hpp#L63) | block | Simulation / SIL harness (host) | Cached discrete ZOH maps Ad, Bd for step size h |
| [`canon`](inc/damp/matlab.hpp#L501) | function | MATLAB®-style aliases (host) | MATLAB® canon(sys) / canon(sys, 'modal'): real modal canonical form and its transform (+1 more overload) |
| [`canonical_phase_margin`](inc/damp/analysis/frequency.hpp#L150) | function | Frequency-domain analysis (host) | Normalize phase margin to (-180, 180] |
| [`care`](inc/damp/design/riccati.hpp#L1490) | function | Design-time synthesis (not PWM-rate) | Solve the Continuous-time Algebraic Riccati Equation (CARE) |
| [`carlson_update`](inc/damp/estimation/sqrt_covariance.hpp#L608) | function | Observers & estimators | Carlson's triangular update for row m of a whitened reading: S ← S T with T Tᵀ = I − f fᵀ/α, f = Sᵀc̃ᵀ the row of CS and α = 1 + fᵀf |
| [`Cascade`](inc/damp/controllers/composition.hpp#L62) | block | Runtime controllers | Series cascade of two SISO controllers: outer → inner reference |
| [`catmull_1d`](inc/damp/toolbox/lookup.hpp#L339) | function | Embedded helpers (controls-adjacent utilities) | 1-D non-uniform Catmull-Rom (cubic Hermite) evaluation |
| [`cbf_relative_degree_1`](inc/damp/controllers/action_governor.hpp#L203) | function | Design-time synthesis (not PWM-rate) | Build a relative-degree-1 CBF inequality row |
| [`CBFConstraint`](inc/damp/controllers/action_governor.hpp#L187) | block | Runtime controllers | One affine row from a relative-degree-1 CBF condition |
| [`cbrt`](inc/damp/math/math.hpp#L91) | function | Scalar math, complex & frames | Cube root (preserves sign for negative x) |
| [`ceil`](inc/damp/math/math.hpp#L380) | function | Scalar math, complex & frames | Ceiling — smallest integer ≥ x |
| [`chirp`](inc/damp/estimation/excitation/chirp.hpp#L120) | function | Observers & estimators | Build a chirp design payload from a configuration |
| [`Chirp`](inc/damp/estimation/excitation/chirp.hpp#L138) | block | Observers & estimators | Linear or logarithmic chirp runtime generator (+1 more overload) |
| [`ChirpConfig`](inc/damp/estimation/excitation/chirp.hpp#L38) | block | Observers & estimators | Configuration for a sine chirp excitation |
| [`ChirpMode`](inc/damp/estimation/excitation/chirp.hpp#L24) | enum | Observers & estimators | Chirp sweep law |
| [`ChirpResult`](inc/damp/estimation/excitation/chirp.hpp#L81) | block | Observers & estimators | Chirp design payload |
| [`cholesky`](inc/damp/matrix/decomposition.hpp#L129) | function | Linear algebra | Cholesky decomposition for positive-definite matrices |
| [`cholesky_solve`](inc/damp/matrix/solve.hpp#L163) | function | Linear algebra | Solve AX = B via Cholesky (A = LLᴴ) |
| [`clamp`](inc/damp/backend.hpp#L376) | function | Core, configuration & backend vocabulary | x clamped to [lo, hi], by value: min(max(lo, x), hi) |
| [`clarke_park_transform`](inc/damp/math/transforms.hpp#L403) | function | Scalar math, complex & frames | Fused Clarke-Park transform (abc → dq) |
| [`clarke_transform`](inc/damp/math/transforms.hpp#L276) | function | Scalar math, complex & frames | Clarke transform (abc → αβ) |
| [`classical_dob`](inc/damp/estimation/dob.hpp#L411) | function | Observers & estimators | Synthesize a classical disturbance observer from a nominal plant and Q-filter |
| [`ClassicalDOB`](inc/damp/estimation/dob.hpp#L446) | block | Observers & estimators | Classical Pn^-1·Q disturbance observer runtime (bolt-on compensator) (+1 more overload) |
| [`ClassicalDobResult`](inc/damp/estimation/dob.hpp#L346) | block | Observers & estimators | Design result for the classical Pn^-1·Q disturbance observer |
| [`classify_range`](inc/damp/toolbox/conditioning.hpp#L367) | function | Embedded helpers (controls-adjacent utilities) | Classify x against the valid span and the wire-fault limits of range |
| [`closed_loop_poles`](inc/damp/design/stability.hpp#L570) | function | Design-time synthesis (not PWM-rate) | Compute closed-loop poles (eigenvalues) with state feedback |
| [`cohen_coon`](inc/damp/design/pid_design.hpp#L304) | function | Design-time synthesis (not PWM-rate) | Cohen-Coon tuning from first-order-plus-dead-time model |
| [`ColVec`](inc/damp/matrix/colvec.hpp#L28) | block | Linear algebra | Concrete Column vector specialization of Matrix<N, 1, T> |
| [`ColView`](inc/damp/matrix/views.hpp#L331) | block | Linear algebra | Non-owning column view of a matrix |
| [`comb_notch_window`](inc/damp/filters/moving_average.hpp#L34) | function | Filters & signal conditioning | Window length for a moving-average comb that notches f_notch and all its harmonics: N = round(1 / (f_notch·Ts)) |
| [`combine_quadrature`](inc/damp/encoder/sincos.hpp#L119) | function | Embedded helpers (controls-adjacent utilities) | Absolute position from a quadrature edge counter and the fine phase |
| [`CommandProjectionResult`](inc/damp/controllers/action_governor.hpp#L89) | block | Design-time synthesis (not PWM-rate) | Result of an affine command projection |
| [`companion_form`](inc/damp/systems/realization.hpp#L614) | function | LTI systems (SS / TF / ZPK / discretize) | Companion form of a SISO model, MATLAB®'s canon(sys, 'companion'): B = e₁, −α in A's last column |
| [`Complementary`](inc/damp/filters/complementary.hpp#L32) | block | Filters & signal conditioning | Scalar (1-D) complementary filter — fuse a fast rate with a slow absolute |
| [`ComplementaryFilter`](inc/damp/estimation/sensor_fusion.hpp#L59) | block | Observers & estimators | Simple complementary filter for orientation estimation |
| [`complex`](inc/damp/math/complex.hpp#L38) | block | Scalar math, complex & frames | Constexpr complex number class for compile-time computations |
| [`compute_eigenvalues`](inc/damp/matrix/eigen.hpp#L1155) | function | Linear algebra | Compute the eigenvalues (and Schur vectors) of a real square matrix |
| [`constant_power_torque_limit`](inc/damp/motor/foc.hpp#L41) | function | Motor control pack (if present) | Constant-power torque ceiling Tₘₐₓ = P_rated / \|ω\| |
| [`ConstantInertiaFeedforward`](inc/damp/toolbox/actuator.hpp#L199) | block | Embedded helpers (controls-adjacent utilities) | Per-axis decoupled torque feedforward: `τ = J·a + b·v + τ_c·sign(v) + g` |
| [`ConstMatrix`](inc/damp/matrix/const_matrix.hpp#L87) | block | Linear algebra | A matrix that is a constant of its type: no storage, and a product with a vector that has no term for an entry that is exactly 0 |
| [`continuous_lpf_exact_step`](inc/damp/filters/lowpass.hpp#L36) | function | Filters & signal conditioning | Exact step of continuous first-order LPF ẏ = −ω_c (y − u) |
| [`continuous_lqr`](inc/damp/design/riccati.hpp#L1748) | function | Design-time synthesis (not PWM-rate) | Continuous-time Linear-Quadratic Regulator design |
| [`ContinuousPID`](inc/damp/controllers/pid.hpp#L1185) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its design, keeping the result's scalar |
| [`ContinuousPID`](inc/damp/controllers/pid.hpp#L1042) | block | Runtime controllers | Continuous-gain PID with per-tick sample time (variable-rate secondary form) |
| [`controllability_gramian`](inc/damp/design/stability.hpp#L381) | function | Design-time synthesis (not PWM-rate) | Continuous/discrete controllability Gramian W_c |
| [`controllability_matrix`](inc/damp/design/stability.hpp#L53) | function | Design-time synthesis (not PWM-rate) | Compute the controllability matrix [B, AB, A²B, ..., A^(N-1)B] |
| [`controllability_staircase`](inc/damp/systems/realization.hpp#L789) | function | LTI systems (SS / TF / ZPK / discretize) | Controllability staircase form: the controllable part split off by an orthogonal transform |
| [`controllable_canonical_form`](inc/damp/systems/realization.hpp#L516) | function | LTI systems (SS / TF / ZPK / discretize) | Controllable canonical form of a SISO model: ones on A's superdiagonal, −α in its last row, B = eₙ |
| [`Convention`](inc/damp/math/transforms.hpp#L68) | enum | Scalar math, complex & frames | Scaling convention for the Clarke/Park family |
| [`ConvertsExactly`](inc/damp/matrix/matrix_traits.hpp#L77) | concept | Linear algebra | Matrix elements of type From convert to To without rounding: the same type, or a real into the complex of its own precision |
| [`copysign`](inc/damp/math/math.hpp#L481) | function | Scalar math, complex & frames | Copy sign — magnitude of mag with the sign of sgn_src |
| [`cos`](inc/damp/matrix/functions.hpp#L1352) | function | Linear algebra | Matrix cosine via scaling and double-angle reconstruction |
| [`cos`](inc/damp/math/math.hpp#L168) | function | Scalar math, complex & frames | Cosine |
| [`cosh`](inc/damp/matrix/functions.hpp#L1386) | function | Linear algebra | Matrix hyperbolic cosine cosh(A) = (exp(A) + exp(−A))/2 |
| [`Counter`](inc/damp/toolbox/logic.hpp#L272) | block | Embedded helpers (controls-adjacent utilities) | Edge-counting up/down counter: increments on each rising edge of up, decrements on each rising edge of down. Returns the running count, saturating at the limits of T (no wrap / signed overflow) |
| [`covariance_factor_or_identity`](inc/damp/estimation/sqrt_covariance.hpp#L168) | function | Observers & estimators | The factor a filter starts from: of P when it is a covariance; otherwise refused under constant evaluation, and the identity at run time |
| [`covariance_of_factor`](inc/damp/estimation/sqrt_covariance.hpp#L69) | function | Observers & estimators | The covariance P = S Sᵀ of a lower-triangular factor S, exactly symmetric |
| [`cross_check_analog`](inc/damp/toolbox/conditioning.hpp#L541) | function | Embedded helpers (controls-adjacent utilities) | Cross-check two calibrated analog readings for redundant sensors |
| [`CTD`](inc/damp/toolbox/iec61131.hpp#L351) | block | Embedded helpers (controls-adjacent utilities) | CTD Counter (Count Down) |
| [`ctrb`](inc/damp/matlab.hpp#L387) | function | MATLAB®-style aliases (host) | MATLAB® short alias for controllability_matrix (+1 more overload) |
| [`ctrbf`](inc/damp/matlab.hpp#L542) | function | MATLAB®-style aliases (host) | MATLAB® ctrbf(A, B, C): controllability staircase form, its orthogonal transform and rank |
| [`CTU`](inc/damp/toolbox/iec61131.hpp#L311) | block | Embedded helpers (controls-adjacent utilities) | CTU Counter (Count Up) |
| [`CTUD`](inc/damp/toolbox/iec61131.hpp#L391) | block | Embedded helpers (controls-adjacent utilities) | CTUD Counter (Count Up Down) |
| [`current_loop_pi`](inc/damp/design/pid_design.hpp#L1302) | function | Design-time synthesis (not PWM-rate) | PI gains for an R–L current loop 1/(L s + R) at any bandwidth ω |
| [`d2d`](inc/damp/matlab.hpp#L459) | function | MATLAB®-style aliases (host) | MATLAB® d2d: resample a ZOH-discrete model to a new sample time |
| [`damp`](inc/damp/analysis/poles.hpp#L70) | function | Frequency-domain analysis (host) | Compute natural frequency and damping for each pole |
| [`damping_ratio_from_overshoot_percent`](inc/damp/design/pid_design.hpp#L862) | function | Design-time synthesis (not PWM-rate) | Map percent overshoot target to equivalent damping ratio |
| [`dare`](inc/damp/design/riccati.hpp#L1321) | function | Design-time synthesis (not PWM-rate) | Solve the Discrete Algebraic Riccati Equation (DARE) |
| [`DareCost`](inc/damp/controllers/mpc.hpp#L330) | block | Runtime controllers | A terminal cost for design::mpc, among its options: at most one terminal ingredient |
| [`DareMethod`](inc/damp/design/riccati.hpp#L1188) | enum | Design-time synthesis (not PWM-rate) | DARE solver selection for dare() |
| [`db2mag`](inc/damp/math/math.hpp#L549) | function | Scalar math, complex & frames | Decibels to magnitude, 10^(db/20) |
| [`dcgain`](inc/damp/analysis/norms.hpp#L42) | function | Frequency-domain analysis (host) | Compute DC gain of a continuous-time system |
| [`DCM`](inc/damp/math/geometry.hpp#L83) | block | Scalar math, complex & frames | Rotation matrix, active body→parent (SO(3) wrapper over Mat3) |
| [`deadband`](inc/damp/toolbox/conditioning.hpp#L47) | function | Embedded helpers (controls-adjacent utilities) | Dead zone over `[lower, upper]`, matching Simulink®'s Dead Zone block (+1 more overload) |
| [`Debounce`](inc/damp/toolbox/logic.hpp#L211) | block | Embedded helpers (controls-adjacent utilities) | Debounce: the output adopts in only after in differs from the current output continuously for stable_time. Rejects contact bounce and brief glitches. (Not an IEC block — the one everyone hand-rolls.) |
| [`decode_sincos`](inc/damp/encoder/sincos.hpp#L91) | function | Embedded helpers (controls-adjacent utilities) | Trim a (sin, cos) pair and decode its phase |
| [`decoupling_feedforward`](inc/damp/motor/foc.hpp#L69) | function | Motor control pack (if present) | Decoupling + back-EMF feedforward at R = 0 (reference currents) |
| [`default_tol`](inc/damp/matrix/matrix_traits.hpp#L122) | function | Linear algebra | Type-appropriate default tolerance for floating-point comparisons. float  ~7 decimal digits  → 1e-6 double ~15 decimal digits → 1e-12 |
| [`deg2rad`](inc/damp/math/math.hpp#L571) | function | Scalar math, complex & frames | Degrees to radians, deg·π/180 |
| [`Delay`](inc/damp/filters/delay.hpp#L26) | block | Filters & signal conditioning | Discrete-time delay buffer |
| [`DemandLimit`](inc/damp/toolbox/bounds.hpp#L218) | block | Embedded helpers (controls-adjacent utilities) | A limit on what a state-feedback law's outer terms ask of one of its states: a cascade's outer limit, inside one gain |
| [`DesignStatus`](inc/damp/design/status.hpp#L34) | enum | Design-time synthesis (not PWM-rate) | Why a design failed, one value per cause class (not per code branch) |
| [`det`](inc/damp/matrix/functions.hpp#L174) | function | Linear algebra | Matrix determinant det(A) |
| [`DFF`](inc/damp/toolbox/iec61131.hpp#L448) | block | Embedded helpers (controls-adjacent utilities) | D Flip-Flop (edge-triggered data latch) |
| [`diag`](inc/damp/matrix/concat.hpp#L53) | function | Linear algebra | Square diagonal matrix from scalars: diag(10.0, 0.1, 100.0, 10.0) (+4 more overloads) |
| [`Diagonal`](inc/damp/matrix/views.hpp#L30) | block | Linear algebra | Diagonal view of a square matrix |
| [`differential_drive`](inc/damp/toolbox/io.hpp#L262) | function | Embedded helpers (controls-adjacent utilities) | Differential (tank) drive mixer: throttle + turn → left/right |
| [`DirectQuadrature`](inc/damp/math/transforms.hpp#L83) | block | Scalar math, complex & frames | Direct-quadrature (rotor-frame) component pair |
| [`Discrete`](inc/damp/simulation/integrator.hpp#L151) | block | Simulation / SIL harness (host) | Discrete-time integrator (no integration, just one step) |
| [`discrete_lqr`](inc/damp/design/riccati.hpp#L1614) | function | Design-time synthesis (not PWM-rate) | Discrete-time Linear-Quadratic Regulator design (+1 more overload) |
| [`discrete_lyapunov_cholesky`](inc/damp/design/lyapunov.hpp#L518) | function | Design-time synthesis (not PWM-rate) | Cholesky factor of the discrete Lyapunov solution: A X Aᵀ − X + B Bᵀ = 0, X = RᵀR |
| [`DiscretePIDGains`](inc/damp/controllers/pid.hpp#L86) | block | Runtime controllers | Fixed-rate discrete PID coefficients (canonical deploy form) |
| [`DiscretizationMethod`](inc/damp/systems/discretization.hpp#L28) | enum | LTI systems (SS / TF / ZPK / discretize) | Discretization methods for continuous-time state-space systems |
| [`DiscretizationResult`](inc/damp/systems/discretization.hpp#L51) | block | LTI systems (SS / TF / ZPK / discretize) | A discretized or resampled model with its success / status, read like an optional |
| [`discretize`](inc/damp/systems/discretization.hpp#L331) | function | LTI systems (SS / TF / ZPK / discretize) | Discretize a continuous-time state-space system (+3 more overloads) |
| [`discretize_lqr_cost`](inc/damp/design/riccati.hpp#L1683) | function | Design-time synthesis (not PWM-rate) | Discretize a continuous LQR cost integral over one sample (Van Loan) |
| [`discretize_process_noise`](inc/damp/estimation/kalman.hpp#L362) | function | Observers & estimators | Discretize continuous process noise over one sample (Van Loan) |
| [`disk_margin`](inc/damp/analysis/disk_margin.hpp#L158) | function | Frequency-domain analysis (host) | Disk margin of a SISO loop L under unity negative feedback |
| [`diskmargin`](inc/damp/matlab.hpp#L1343) | function | MATLAB®-style aliases (host) | Disk margin of a SISO loop (+1 more overload) |
| [`DiskMarginResult`](inc/damp/analysis/disk_margin.hpp#L36) | block | Frequency-domain analysis (host) | Disk margin of a SISO loop: the disk size α and the gain and phase margins it guarantees |
| [`DisturbanceModel`](inc/damp/design/disturbance.hpp#L128) | block | Design-time synthesis (not PWM-rate) | An explicit disturbance model: d enters the plant's state through B_d and the readings through C_d (+1 more overload) |
| [`DisturbanceModelFor`](inc/damp/design/disturbance.hpp#L331) | concept | Design-time synthesis (not PWM-rate) | A disturbance spec on a plant of NX states and NY outputs in T: a DisturbanceModel of those sizes and scalar, or a spec that takes the plant's own |
| [`DisturbanceSpec`](inc/damp/design/disturbance.hpp#L315) | concept | Design-time synthesis (not PWM-rate) | A disturbance spec with states: InputDisturbance, OutputDisturbance, PlantDisturbance or DisturbanceModel (every spec but NoDisturbance) |
| [`DLATCH`](inc/damp/toolbox/iec61131.hpp#L479) | block | Embedded helpers (controls-adjacent utilities) | D Latch (level-sensitive / transparent latch) |
| [`dlqe`](inc/damp/matlab.hpp#L971) | function | MATLAB®-style aliases (host) | Discrete linear-quadratic estimator from matrices (MATLAB®'s dlqe) |
| [`dlqr`](inc/damp/matlab.hpp#L919) | function | MATLAB®-style aliases (host) | Discrete-time Linear-Quadratic Regulator design (MATLAB®'s dlqr) |
| [`dlyap`](inc/damp/design/lyapunov.hpp#L478) | function | Design-time synthesis (not PWM-rate) | Solve the discrete-time Lyapunov (Stein) equation A X Aᵀ − X + Q = 0 |
| [`dlyap`](inc/damp/matlab.hpp#L1434) | function | MATLAB®-style aliases (host) | MATLAB® alias for the discrete Lyapunov solve AXAᵀ−X+Q=0 |
| [`dlyapchol`](inc/damp/matlab.hpp#L1454) | function | MATLAB®-style aliases (host) | MATLAB® alias for the Cholesky factor of the discrete Lyapunov solution (X = RᵀR) |
| [`dob`](inc/damp/estimation/dob.hpp#L150) | function | Observers & estimators | Validate and package DOB configuration into a runtime-ready design result |
| [`DOB`](inc/damp/estimation/dob.hpp#L175) | block | Observers & estimators | Lightweight SISO disturbance observer runtime (+1 more overload) |
| [`DOBConfig`](inc/damp/estimation/dob.hpp#L30) | block | Observers & estimators | Configuration for a first-order disturbance observer |
| [`DP45`](inc/damp/simulation/integrator.hpp#L795) | block | Simulation / SIL harness (host) | Dormand-Prince 5(4) adaptive integrator — the ode45 pair |
| [`dpwm_zero_sequence`](inc/damp/motor/modulation.hpp#L267) | function | Motor control pack (if present) | DPWM zero-sequence for a chosen discontinuous scheme |
| [`DqCommand`](inc/damp/motor/foc.hpp#L83) | block | Motor control pack (if present) | Clamped dq voltage command from FOController::current_controller |
| [`DrivePair`](inc/damp/toolbox/io.hpp#L234) | block | Embedded helpers (controls-adjacent utilities) | Left/right actuator pair (differential/arcade drive output) |
| [`DsogiPll`](inc/damp/filters/pll.hpp#L313) | block | Filters & signal conditioning | Dual-SOGI three-phase positive-sequence PLL (DSOGI-PLL) |
| [`duties_from_phase_voltages`](inc/damp/motor/modulation.hpp#L362) | function | Motor control pack (if present) | Map phase voltages + zero-sequence to clamped half-bridge duties |
| [`Dynamics`](inc/damp/concepts.hpp#L200) | concept | Core, configuration & backend vocabulary | Continuous-time nonlinear dynamics ẋ = f(x, u): a callable returning ẋ as a ColVec<NX, T> |
| [`DynamicsJacobian`](inc/damp/concepts.hpp#L219) | block | Core, configuration & backend vocabulary | Analytic Jacobians of a Dynamics model at (x, u), returned by its jacobian(x, u) |
| [`eig`](inc/damp/matlab.hpp#L670) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the eigenvalues of a square matrix |
| [`EigenResult`](inc/damp/matrix/eigen.hpp#L40) | block | Linear algebra | Eigenvalue computation result |
| [`eliminate_states`](inc/damp/design/model_reduction.hpp#L632) | function | Design-time synthesis (not PWM-rate) | Eliminate the listed states, by DC-matched residualization or truncation |
| [`EncoderGate`](inc/damp/encoder/gate.hpp#L108) | block | Embedded helpers (controls-adjacent utilities) | Validates encoder readings before they reach the tracker |
| [`EncoderGateConfig`](inc/damp/encoder/gate.hpp#L42) | block | Embedded helpers (controls-adjacent utilities) | Configuration for EncoderGate. Positions in counts, as the tracker's |
| [`EncoderTracker`](inc/damp/encoder/tracker.hpp#L70) | block | Embedded helpers (controls-adjacent utilities) | Encoder position (exact) and speed (smoothed) from absolute counts |
| [`EncoderTrackerConfig`](inc/damp/encoder/tracker.hpp#L23) | block | Embedded helpers (controls-adjacent utilities) | Configuration for EncoderTracker |
| [`EncoderVerdict`](inc/damp/encoder/gate.hpp#L30) | enum | Embedded helpers (controls-adjacent utilities) | What EncoderGate made of a reading |
| [`enu_from_ned`](inc/damp/estimation/ins_mechanization.hpp#L115) | function | Observers & estimators | Map a NED vector into ENU (axis permute) |
| [`eps`](inc/damp/matrix/matrix_traits.hpp#L112) | function | Linear algebra | Machine epsilon of T's real scalar (MATLAB® eps). Scale it by problem size and magnitude to build a relative tolerance |
| [`ErrorStateJacobian`](inc/damp/estimation/eskf.hpp#L285) | block | Observers & estimators | ESKF prediction Jacobians F and G (nominal state updated outside) |
| [`ErrorStateKalmanFilter`](inc/damp/estimation/eskf.hpp#L338) | block | Observers & estimators | Runtime error-state KF: tracks δx and P, not the full nominal state (+1 more overload) |
| [`ErrorStateKalmanFilterConfig`](inc/damp/estimation/eskf.hpp#L309) | block | Observers & estimators | Configuration of an ErrorStateKalmanFilter: the process noise and the initial covariance |
| [`ErrorTolerance`](inc/damp/simulation/solver.hpp#L302) | block | Simulation / SIL harness (host) | Scalar-or-vector local-error tolerance for adaptive ODE control |
| [`esc`](inc/damp/controllers/esc.hpp#L156) | function | Design-time synthesis (not PWM-rate) | Synthesize an extremum-seeking controller |
| [`esc_mppt`](inc/damp/controllers/esc.hpp#L200) | function | Design-time synthesis (not PWM-rate) | MPPT-flavored ESC: maximize a power measurement by perturbing the operating point (e.g. converter duty or reference voltage) |
| [`ESCResult`](inc/damp/controllers/esc.hpp#L90) | block | Design-time synthesis (not PWM-rate) | Extremum-seeking design result: the discrete realization the runtime runs |
| [`ESCSpec`](inc/damp/controllers/esc.hpp#L70) | block | Runtime controllers | The extremum-seeking design input: dither, gain, filter corners and clamp; Ts is passed beside it |
| [`eskf_imu`](inc/damp/estimation/eskf.hpp#L217) | function | Observers & estimators | Design Q, R, P₀ for a 6-state IMU attitude ESKF (gyro + accel) |
| [`eskf_marg`](inc/damp/estimation/eskf.hpp#L256) | function | Observers & estimators | Design Q, R, P₀ for a 6-state MARG attitude ESKF (gyro + accel + mag) |
| [`ESKFOrientationFilter`](inc/damp/estimation/sensor_fusion.hpp#L381) | block | Observers & estimators | Turnkey 6-state attitude ESKF: owns q, b_g, and the error filter |
| [`ESKFPredictFn`](inc/damp/estimation/eskf.hpp#L297) | concept | Observers & estimators | Callable that returns ErrorStateJacobian for one predict step |
| [`ESKFResult`](inc/damp/estimation/eskf.hpp#L83) | block | Observers & estimators | Design payload for an ESKF: Q, R, P₀, and success |
| [`eso`](inc/damp/estimation/eso.hpp#L208) | function | Observers & estimators | Extended-state observer: poles at −ω_o, plant gain b₀, optional sample time |
| [`ESO`](inc/damp/estimation/eso.hpp#L311) | block | Observers & estimators | Extended-state observer (ESO) (+1 more overload) |
| [`ESOResult`](inc/damp/estimation/eso.hpp#L53) | block | Observers & estimators | Extended-state observer design: continuous β and, for Ts > 0, the discrete ESO |
| [`estim`](inc/damp/matlab.hpp#L732) | function | MATLAB®-style aliases (host) | Form the state estimator from a plant and an estimator gain |
| [`Euler`](inc/damp/math/geometry.hpp#L223) | block | Scalar math, complex & frames | Intrinsic Euler angles for sequence Order (default ZYX aerospace YPR) |
| [`EulerOrder`](inc/damp/math/geometry.hpp#L52) | enum | Scalar math, complex & frames | Intrinsic Euler angle sequence (axis order of successive principal rotations) |
| [`eval_frf`](inc/damp/systems/state_space.hpp#L286) | function | LTI systems (SS / TF / ZPK / discretize) | Evaluate frequency response of a state-space system |
| [`EventSimConfig`](inc/damp/simulation/hybrid.hpp#L67) | block | Simulation / SIL harness (host) | Configuration for event-driven hybrid / multi-rate style runs |
| [`every_output`](inc/damp/design/detail/config.hpp#L107) | function | Design-time synthesis (not PWM-rate) | The default of a list of N outputs: left out |
| [`Exact`](inc/damp/simulation/integrator.hpp#L185) | block | Simulation / SIL harness (host) | Exact integrator for LTI systems |
| [`exact_lti_step`](inc/damp/simulation/hybrid.hpp#L304) | function | Simulation / SIL harness (host) | Exact open-loop step of a single continuous LTI system |
| [`exp`](inc/damp/math/math.hpp#L231) | function | Scalar math, complex & frames | Exponential function |
| [`expm`](inc/damp/matrix/functions.hpp#L697) | function | Linear algebra | Matrix exponential by scaling and squaring (Al-Mohy & Higham 2009) |
| [`expm`](inc/damp/matlab.hpp#L682) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the matrix exponential (mat::expm) |
| [`expo`](inc/damp/toolbox/conditioning.hpp#L202) | function | Embedded helpers (controls-adjacent utilities) | Exponential response curve `y = (1−k)·x + k·x³` (RC "expo") |
| [`ExtendedKalmanFilter`](inc/damp/estimation/ekf.hpp#L129) | block | Observers & estimators | Extended Kalman filter: continuous nonlinear dynamics, sampled measurements (+1 more overload) |
| [`ExtendedKalmanFilterConfig`](inc/damp/estimation/ekf.hpp#L41) | block | Observers & estimators | Configuration of an ExtendedKalmanFilter: its sample time, the process noise and the initial estimate |
| [`extract_channel`](inc/damp/simulation/plot_plotly.hpp#L225) | function | Simulation / SIL harness (host) | Extract the i-th element from each vector entry into a std::vector<double> |
| [`extract_frf_features`](inc/damp/estimation/frequency_response.hpp#L808) | function | Observers & estimators | Extract peaks and valleys in one pass over the FRF table |
| [`extract_modes`](inc/damp/estimation/frequency_response.hpp#L783) | function | Observers & estimators | Extract resonant peaks (local \|G\| maxima) from an FRF table |
| [`extract_valleys`](inc/damp/estimation/frequency_response.hpp#L797) | function | Observers & estimators | Extract anti-resonance valleys (local \|G\| minima) from an FRF table |
| [`Extrapolation`](inc/damp/toolbox/lookup.hpp#L35) | enum | Embedded helpers (controls-adjacent utilities) | Out-of-range behaviour for a Lut1D / Lut2D / Lut3D query beyond its breakpoints |
| [`ExtremumSeekingController`](inc/damp/controllers/esc.hpp#L384) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its design, keeping the result's scalar |
| [`ExtremumSeekingController`](inc/damp/controllers/esc.hpp#L236) | block | Runtime controllers | Extremum-seeking controller runtime (model-free online optimizer) |
| [`ExtremumType`](inc/damp/controllers/esc.hpp#L48) | enum | Runtime controllers | Whether ESC climbs to a maximum or descends to a minimum of the objective |
| [`eye`](inc/damp/matlab.hpp#L610) | function | MATLAB®-style aliases (host) | Create an identity matrix of size n x n |
| [`F_TRIG`](inc/damp/toolbox/iec61131.hpp#L151) | block | Embedded helpers (controls-adjacent utilities) | F_TRIG (Falling Edge Trigger) |
| [`factorial`](inc/damp/trajectory/trajectory_types.hpp#L208) | function | Trajectory value types | k! (exact for modest k in float/double; used up through nonic BVP / jets) |
| [`falling_factorial`](inc/damp/trajectory/trajectory_types.hpp#L221) | function | Trajectory value types | Falling factorial i·(i−1)···(i−k+1) = i! / (i−k)! — the k-th derivative coefficient of tⁱ. Zero when i < k |
| [`FallingEdge`](inc/damp/toolbox/logic.hpp#L47) | block | Embedded helpers (controls-adjacent utilities) | Falling-edge detector: true on the tick x goes true → false |
| [`feedback`](inc/damp/systems/state_space.hpp#L444) | function | LTI systems (SS / TF / ZPK / discretize) | Negative feedback: y = sys1(u − sys2(y)) (+1 more overload) |
| [`feedback_read`](inc/damp/systems/state_space.hpp#L924) | function | LTI systems (SS / TF / ZPK / discretize) | Unity (or weighted) negative feedback around G_loop, read G_out |
| [`FeedbackLaw`](inc/damp/concepts.hpp#L122) | concept | Core, configuration & backend vocabulary | An LQ or MPC law that holds outputs at r: control(r), step(r, m), and track(u, r, m) in place of step while another source drives the plant (u the input applied) |
| [`FeedbackRuntime`](inc/damp/concepts.hpp#L107) | concept | Core, configuration & backend vocabulary | Shared lifecycle of an LQ or MPC runtime: measurement, time update, applied input |
| [`Filtered`](inc/damp/estimation/state_variable_filter.hpp#L65) | block | Observers & estimators | One signal after the state-variable filter: its filtered value and first two derivatives |
| [`FilterSpec`](inc/damp/controllers/state_feedback_types.hpp#L232) | concept | Runtime controllers | The spec of a law's filter: a KalmanNoise for the Kalman gain, or ObserverPoles for a placed one |
| [`finite_non_negative`](inc/damp/math/math.hpp#L523) | function | Scalar math, complex & frames | True if x is finite and non-negative (x ≥ 0) |
| [`finite_positive`](inc/damp/math/math.hpp#L513) | function | Scalar math, complex & frames | True if x is finite and strictly greater than zero |
| [`fir1`](inc/damp/filters/fir.hpp#L244) | function | Filters & signal conditioning | Window-method FIR design (normalized frequency, single cutoff) (+1 more overload) |
| [`fir1_hz`](inc/damp/filters/fir.hpp#L385) | function | Filters & signal conditioning | Window-method FIR design from a cutoff in Hz and the sample time (+1 more overload) |
| [`fir_window`](inc/damp/filters/fir.hpp#L421) | function | Filters & signal conditioning | Alias for fir1 (descriptive name) (+1 more overload) |
| [`FirDesignResult`](inc/damp/filters/fir.hpp#L65) | block | Filters & signal conditioning | Window-method FIR design result |
| [`fires_at`](inc/damp/simulation/multirate.hpp#L73) | function | Simulation / SIL harness (host) | True when base tick tick is a firing instant for divisor |
| [`FirFilter`](inc/damp/filters/fir.hpp#L463) | block | Filters & signal conditioning | Direct-form FIR filter runtime (tapped delay line + dot product) (+1 more overload) |
| [`FirstOrderCoeffs`](inc/damp/filters/iir_design.hpp#L30) | block | Filters & signal conditioning | DSP coefficients for first-order IIR filter |
| [`FirstOrderPlantEstimator`](inc/damp/estimation/parameter_estimation.hpp#L125) | block | Observers & estimators | Online estimator for a first-order plant's gain and time constant |
| [`FirstOrderPlantEstimatorConfig`](inc/damp/estimation/parameter_estimation.hpp#L83) | block | Observers & estimators | Configuration for the first-order grey-box estimator |
| [`FirstOrderSection`](inc/damp/filters/first_order_section.hpp#L43) | block | Filters & signal conditioning | First-order IIR section runtime, Direct Form I (+1 more overload) |
| [`FirType`](inc/damp/filters/fir.hpp#L36) | enum | Filters & signal conditioning | Ideal frequency-selective response for fir1 / fir_window |
| [`FirWindow`](inc/damp/filters/fir.hpp#L48) | enum | Filters & signal conditioning | Window applied to the ideal truncated impulse response |
| [`fixed_solve`](inc/damp/simulation/solver.hpp#L286) | function | Simulation / SIL harness (host) | Fixed-step solve in one call |
| [`FixedGainObserver`](inc/damp/estimation/fixed_gain_observer.hpp#L616) | block | Observers & estimators | The fixed-gain state observer: the innovation gain M held, the current form, and the optional disturbance anchor (+4 more overloads) |
| [`FixedStepSolver`](inc/damp/simulation/solver.hpp#L145) | block | Simulation / SIL harness (host) | Fixed-step ODE solver (+1 more overload) |
| [`floor`](inc/damp/math/math.hpp#L368) | function | Scalar math, complex & frames | Floor — largest integer ≤ x |
| [`flux_from_Kv`](inc/damp/motor/spm.hpp#L54) | function | Motor control pack (if present) | λ from Kᵥ via torque_constant_from_Kv |
| [`flux_from_torque_constant`](inc/damp/motor/spm.hpp#L36) | function | Motor control pack (if present) | λ from datasheet Kₜ (amplitude / peak-per-phase convention) |
| [`fmod`](inc/damp/math/math.hpp#L412) | function | Scalar math, complex & frames | Floating-point remainder, x − y·trunc(x/y) (sign of x), matching std::fmod's truncated-quotient convention |
| [`FOController`](inc/damp/motor/foc.hpp#L98) | block | Motor control pack (if present) | dq current regulator: PI + decoupling FF → clamped Vdq |
| [`forward_substitute`](inc/damp/matrix/solve.hpp#L53) | function | Linear algebra | Forward substitution for Lx = b |
| [`ForwardEuler`](inc/damp/simulation/integrator.hpp#L241) | block | Simulation / SIL harness (host) | Forward Euler integrator |
| [`fractional_delay`](inc/damp/systems/discretization.hpp#L513) | function | LTI systems (SS / TF / ZPK / discretize) | Fractional-delay all-pass filter of order N (Thiran): z⁻ᴰ for a real D = delay / Ts |
| [`FrequencyPoint`](inc/damp/analysis/frequency.hpp#L37) | block | Frequency-domain analysis (host) | Single-point frequency response result |
| [`FrequencyResponseEstimator`](inc/damp/estimation/frequency_response.hpp#L166) | block | Observers & estimators | Lock-in FRF estimator over a fixed table of NFreq bins |
| [`frf_s_or_z`](inc/damp/analysis/frequency.hpp#L1019) | function | Frequency-domain analysis (host) | FRF evaluation point: jω (continuous) or e^{jωTs} (discrete) |
| [`FrfFeatureKind`](inc/damp/estimation/frequency_response.hpp#L435) | enum | Observers & estimators | Peak (resonance) or valley (anti-resonance) on an FRF magnitude curve |
| [`FrfFeatureResult`](inc/damp/estimation/frequency_response.hpp#L518) | block | Observers & estimators | Combined peak + valley reduction of an FRF table |
| [`FrfMargins`](inc/damp/estimation/frequency_response.hpp#L553) | block | Observers & estimators | Open-loop margins read from a discrete FRF table (on-target) |
| [`FrfPoint`](inc/damp/estimation/frequency_response.hpp#L55) | block | Observers & estimators | One measured frequency-response point (on-target FRF table entry) |
| [`full_qr`](inc/damp/matrix/decomposition.hpp#L808) | function | Linear algebra | Full QR factorization via Householder reflections (real or complex T) |
| [`FullQR`](inc/damp/matrix/decomposition.hpp#L790) | block | Linear algebra | Result of a full (complete) QR factorization |
| [`gain_margin_unwrapped`](inc/damp/analysis/frequency.hpp#L211) | function | Frequency-domain analysis (host) | Find the gain margin using the unwrapped phase trajectory |
| [`GainSpec`](inc/damp/controllers/state_feedback_types.hpp#L228) | concept | Runtime controllers | The spec of a law's gain: LQRWeights for the LQ gain, or Poles for a placed one. The weights are an argument, never an option |
| [`gang_of_four`](inc/damp/analysis/sensitivity.hpp#L207) | function | Frequency-domain analysis (host) | Gang of four of a plant P and a controller C on the error, u = C(r − y) |
| [`gang_of_seven`](inc/damp/analysis/sensitivity.hpp#L150) | function | Frequency-domain analysis (host) | Gang of seven of a plant and a two-degree-of-freedom controller u = K[r; y] |
| [`GangOfFourResult`](inc/damp/analysis/sensitivity.hpp#L29) | block | Frequency-domain analysis (host) | Gang of four of a loop u = C(r − y): S, T, PS and CS, each a state-space model on the closed loop's states |
| [`GangOfSevenResult`](inc/damp/analysis/sensitivity.hpp#L45) | block | Frequency-domain analysis (host) | Gang of seven of a two-degree-of-freedom loop u = K[r; y]: the gang of four of its feedback path plus the reference responses |
| [`generalized_eigenvalues`](inc/damp/matrix/generalized_eigen.hpp#L622) | function | Linear algebra | Generalized eigenvalues of the real pencil A − λB: the λ with det(A − λB) = 0 |
| [`GeneralizedEigenResult`](inc/damp/matrix/generalized_eigen.hpp#L62) | block | Linear algebra | Generalized eigenvalues of A − λB: the finite ones, and how many there are |
| [`geomspace`](inc/damp/analysis/linspace.hpp#L68) | function | Frequency-domain analysis (host) | Geometric sequence from start to end (endpoint values, not exponents) (+1 more overload) |
| [`Goertzel`](inc/damp/filters/spectral.hpp#L51) | block | Filters & signal conditioning | Generalized Goertzel single-bin DFT — amplitude/phase at one frequency |
| [`gram`](inc/damp/matlab.hpp#L1510) | function | MATLAB®-style aliases (host) | MATLAB® alias for the controllability/observability Gramian of a system |
| [`gravity_nav`](inc/damp/estimation/ins_mechanization.hpp#L92) | function | Observers & estimators | Gravity vector in the local-level frame [m/s²] |
| [`h_pattern`](inc/damp/toolbox/io.hpp#L285) | function | Embedded helpers (controls-adjacent utilities) | H-pattern (dual-lever / tank) drive: one lever per track, no mixing |
| [`Handover`](inc/damp/controllers/state_feedback_types.hpp#L364) | block | Runtime controllers | Bumpless handover for a law with nothing to preload: the first command after track() is the input tracked, and the law's own takes over with the time constant fade_time. Among design::state_feedback's trailing options |
| [`hankel_singular_values`](inc/damp/design/model_reduction.hpp#L513) | function | Design-time synthesis (not PWM-rate) | Hankel singular values of a stable state-space system |
| [`hankelsv`](inc/damp/matlab.hpp#L353) | function | MATLAB®-style aliases (host) | MATLAB® short alias for design::hankel_singular_values |
| [`harmonic_suppressor`](inc/damp/controllers/harmonic_suppression.hpp#L129) | function | Design-time synthesis (not PWM-rate) | Synthesize a multi-resonant harmonic suppressor |
| [`HarmonicAnalyzer`](inc/damp/filters/spectral.hpp#L203) | block | Filters & signal conditioning | Harmonic analyzer — a Goertzel bank over a fundamental and K−1 harmonics |
| [`HarmonicSuppressor`](inc/damp/controllers/harmonic_suppression.hpp#L314) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its design, keeping the result's scalar |
| [`HarmonicSuppressor`](inc/damp/controllers/harmonic_suppression.hpp#L173) | block | Runtime controllers | Multi-resonant harmonic suppressor — a parallel bank of PR resonators |
| [`HarmonicSuppressorResult`](inc/damp/controllers/harmonic_suppression.hpp#L49) | block | Design-time synthesis (not PWM-rate) | Design result for a multi-resonant harmonic suppressor |
| [`hcat`](inc/damp/matrix/concat.hpp#L157) | function | Linear algebra | Horizontal concatenation (same row count) |
| [`heading_from_baseline_nav`](inc/damp/estimation/ins_eskf.hpp#L419) | function | Observers & estimators | Heading [rad] from a nav-frame baseline vector (dual-antenna difference) |
| [`Heun`](inc/damp/simulation/integrator.hpp#L878) | block | Simulation / SIL harness (host) | Heun's method (Improved Euler, RK2) integrator |
| [`highpass_1st`](inc/damp/filters/iir_design.hpp#L422) | function | Filters & signal conditioning | First-order high-pass filter design (continuous-time) (+1 more overload) |
| [`highpass_2nd`](inc/damp/filters/iir_design.hpp#L558) | function | Filters & signal conditioning | Second-order high-pass filter design (continuous-time) (+1 more overload) |
| [`highpass_sections`](inc/damp/filters/biquad.hpp#L498) | function | Filters & signal conditioning | Nth-order Butterworth high-pass at cutoff fc as second-order sections, any order |
| [`highshelf`](inc/damp/filters/iir_design.hpp#L910) | function | Filters & signal conditioning | High-shelf EQ filter: boost or cut everything above fc |
| [`hinfnorm`](inc/damp/matlab.hpp#L1536) | function | MATLAB®-style aliases (host) | MATLAB® alias for the H∞ system norm [ninf, fpeak] = hinfnorm(sys, tol) |
| [`HinfNormResult`](inc/damp/analysis/norms.hpp#L125) | block | Frequency-domain analysis (host) | H∞ norm of a stable system and the frequency of its peak |
| [`holonomic_drive`](inc/damp/toolbox/io.hpp#L386) | function | Embedded helpers (controls-adjacent utilities) | Four-wheel holonomic drive mixer (mecanum / omni) |
| [`HolonomicOutput`](inc/damp/toolbox/io.hpp#L355) | block | Embedded helpers (controls-adjacent utilities) | Four-wheel actuator set (mecanum/holonomic drive output), X-config |
| [`HybridEventAction`](inc/damp/simulation/hybrid.hpp#L104) | block | Simulation / SIL harness (host) | Action applied when an event fires (after Exact advance to the event time) |
| [`HybridSimulationResult`](inc/damp/simulation/hybrid.hpp#L79) | block | Simulation / SIL harness (host) | Result of a piecewise-LTI hybrid simulation |
| [`hypot`](inc/damp/math/math.hpp#L75) | function | Scalar math, complex & frames | Euclidean distance hypot(x, y) = √(x² + y²), without overflow |
| [`Hysteresis`](inc/damp/toolbox/conditioning.hpp#L282) | block | Embedded helpers (controls-adjacent utilities) | Hysteresis comparator (Schmitt trigger): bool output with separate on/off thresholds to reject chatter |
| [`impedance`](inc/damp/analysis/frequency.hpp#L1090) | function | Frequency-domain analysis (host) | Compute impedance frequency response from a SISO impedance transfer function |
| [`impedance_from_admittance`](inc/damp/analysis/frequency.hpp#L1045) | function | Frequency-domain analysis (host) | Compute impedance frequency response from a SISO admittance system |
| [`ImpedanceResult`](inc/damp/analysis/frequency.hpp#L915) | block | Frequency-domain analysis (host) | Result of impedance frequency response evaluation |
| [`impulse`](inc/damp/analysis/time_response.hpp#L198) | function | Frequency-domain analysis (host) | Impulse response of a (MIMO) state-space system (+1 more overload) |
| [`impulse`](inc/damp/estimation/excitation/impulse.hpp#L71) | function | Observers & estimators | Build an impulse design payload |
| [`Impulse`](inc/damp/estimation/excitation/impulse.hpp#L87) | block | Observers & estimators | Rectangular impulse (finite-width Dirac stand-in) |
| [`ImpulseConfig`](inc/damp/estimation/excitation/impulse.hpp#L30) | block | Observers & estimators | Configuration for a one-shot rectangular impulse (Dirac stand-in) |
| [`impulseplot`](inc/damp/simulation/plot_plotly.hpp#L735) | function | Simulation / SIL harness (host) | Plot an impulse response, one trace per input/output pair |
| [`ImpulseResult`](inc/damp/estimation/excitation/impulse.hpp#L50) | block | Observers & estimators | Impulse design payload |
| [`ImuSample`](inc/damp/estimation/ins_mechanization.hpp#L124) | block | Observers & estimators | IMU sample in the body frame |
| [`infinity_norm`](inc/damp/matrix/functions.hpp#L40) | function | Linear algebra | Infinity norm ‖A‖∞: maximum absolute row sum |
| [`initial`](inc/damp/analysis/time_response.hpp#L246) | function | Frequency-domain analysis (host) | Initial-condition (free) response of a (MIMO) state-space system |
| [`InputDelay`](inc/damp/systems/discretization.hpp#L560) | block | LTI systems (SS / TF / ZPK / discretize) | An input delay of θ = seconds, in up to N delay states (0 ≤ θ ≤ N·Ts): each delay state state holds the inputs already sent, v(k−1) … v(k−N). `discretize(sys, Ts, InputDelay<N>{θ})` returns that model; design on it like any other. `InputDelay{θ}` is N = 1 (+1 more overload) |
| [`InputDisturbance`](inc/damp/design/disturbance.hpp#L54) | block | Design-time synthesis (not PWM-rate) | An input disturbance: d adds to u, B_d = B and C_d = D, one disturbance per input (+1 more overload) |
| [`InputLimits`](inc/damp/design/input_limits.hpp#L56) | block | Design-time synthesis (not PWM-rate) | The limits on a law's command u, per input: among design::state_feedback's and design::mpc's options |
| [`InputRateLimits`](inc/damp/design/input_limits.hpp#L92) | block | Design-time synthesis (not PWM-rate) | The limits on a law's move Δu = u_k − u_{k−1} per period, per input: among design::state_feedback's and design::mpc's options |
| [`ins_aid_position`](inc/damp/estimation/ins_eskf.hpp#L386) | function | Observers & estimators | Predicted antenna / marker position with body lever-arm |
| [`ins_apply_correction`](inc/damp/estimation/ins_eskf.hpp#L372) | function | Observers & estimators | Inject filter correction into x and reset the error state |
| [`ins_error_jacobian`](inc/damp/estimation/ins_eskf.hpp#L280) | function | Observers & estimators | Discrete error-state transition for one IMU step, on its sparsity pattern |
| [`ins_eskf_design`](inc/damp/estimation/ins_eskf.hpp#L191) | function | Observers & estimators | Build Q / R / P0 for a 15-state INS ESKF from IMU noise densities |
| [`ins_heading`](inc/damp/estimation/ins_eskf.hpp#L404) | function | Observers & estimators | Horizontal heading [rad] of a body axis expressed in the nav frame |
| [`ins_inject`](inc/damp/estimation/ins_eskf.hpp#L323) | function | Observers & estimators | Inject δx into the nominal INS state (right-multiplicative q) |
| [`ins_predict`](inc/damp/estimation/ins_eskf.hpp#L352) | function | Observers & estimators | Mechanize nominal state and propagate the error covariance |
| [`ins_update_heading`](inc/damp/estimation/ins_eskf.hpp#L507) | function | Observers & estimators | Sparse dual-antenna / yaw heading update (scalar) |
| [`ins_update_orientation`](inc/damp/estimation/ins_eskf.hpp#L545) | function | Observers & estimators | Sparse orientation update from a measured body→nav quaternion |
| [`ins_update_pose`](inc/damp/estimation/ins_eskf.hpp#L569) | function | Observers & estimators | Generic pose alias: position then orientation (two sparse updates) |
| [`ins_update_pose_heading`](inc/damp/estimation/ins_eskf.hpp#L588) | function | Observers & estimators | Pose alias: position + dual-antenna heading (two sparse updates) |
| [`ins_update_position`](inc/damp/estimation/ins_eskf.hpp#L434) | function | Observers & estimators | Sparse position update: y = p + Rℓ |
| [`ins_update_velocity`](inc/damp/estimation/ins_eskf.hpp#L466) | function | Observers & estimators | Sparse velocity update: y = v at the IMU origin (no lever-arm) |
| [`ins_update_zupt`](inc/damp/estimation/ins_eskf.hpp#L485) | function | Observers & estimators | Zero-velocity update (ZUPT): velocity measurement of zero |
| [`InsNavigator`](inc/damp/estimation/ins_eskf.hpp#L622) | block | Observers & estimators | Runtime strapdown navigator: nominal InsState + 15-state ESKF (+2 more overloads) |
| [`InsState`](inc/damp/estimation/ins_mechanization.hpp#L138) | block | Observers & estimators | Nominal strapdown navigation state (external to the ESKF error state) |
| [`instantaneous_power`](inc/damp/math/transforms.hpp#L495) | function | Scalar math, complex & frames | Instantaneous active and reactive power from dq quantities (+1 more overload) |
| [`InstantaneousPower`](inc/damp/math/transforms.hpp#L456) | block | Scalar math, complex & frames | Instantaneous active and reactive power |
| [`IntegralAction`](inc/damp/controllers/state_feedback_types.hpp#L330) | block | Runtime controllers | Integral action: integrators on r − y, with the cost and the limits that go with them. Among a design's trailing options it makes the law an LQI (an LQGI with a KalmanNoise) |
| [`IntegralDesign`](inc/damp/controllers/state_feedback_types.hpp#L371) | block | Runtime controllers | The integral action of a design: the integrated outputs, their maps of the state and the input, the per-sample step, and what back-calculation moves the integrators by |
| [`IntegrationResult`](inc/damp/simulation/integrator.hpp#L61) | block | Simulation / SIL harness (host) | Result of an integration step |
| [`InteriorPointSolver`](inc/damp/design/qp.hpp#L1709) | block | Design-time synthesis (not PWM-rate) | Interior-point QP solver policy for damp::MPC |
| [`Interpolation`](inc/damp/toolbox/lookup.hpp#L41) | enum | Embedded helpers (controls-adjacent utilities) | In-range interpolant for Lut1D (2-D stays bilinear; use SplineSurface for C¹ grids) |
| [`InvariantSet`](inc/damp/controllers/mpc.hpp#L360) | block | Runtime controllers | A terminal cost and set for design::mpc, among its options: DareCost's cost and the maximal positively invariant set of its terminal law, with room for NT_max rows (at most one terminal ingredient) |
| [`InvariantSetResult`](inc/damp/design/invariant_set.hpp#L39) | block | Design-time synthesis (not PWM-rate) | Maximal positively invariant polytope G e ≤ h of x⁺ = A x |
| [`inverse`](inc/damp/matrix/solve.hpp#L316) | function | Linear algebra | Matrix inverse A⁻¹ via mat::solve(A, I) |
| [`inverse_clarke_transform`](inc/damp/math/transforms.hpp#L307) | function | Scalar math, complex & frames | Inverse Clarke transform (αβ → abc) |
| [`inverse_deadband`](inc/damp/toolbox/conditioning.hpp#L79) | function | Embedded helpers (controls-adjacent utilities) | Inverse dead zone: add an offset to overcome a physical dead zone (valve overlap, static friction, motor stiction), with independent negative/positive offsets (+1 more overload) |
| [`inverse_deadband_sine`](inc/damp/toolbox/conditioning.hpp#L162) | function | Embedded helpers (controls-adjacent utilities) | Sine-shaped inverse dead zone (full boost once \|x\|≥w) |
| [`inverse_deadband_soft`](inc/damp/toolbox/conditioning.hpp#L138) | function | Embedded helpers (controls-adjacent utilities) | Smooth inverse dead zone via soft-sign (tanh-like Coulomb boost) |
| [`inverse_lerp`](inc/damp/toolbox/scaling.hpp#L46) | function | Embedded helpers (controls-adjacent utilities) | Inverse of lerp: the fraction t such that `lerp(a, b, t) == x` |
| [`inverse_park_clarke_transform`](inc/damp/math/transforms.hpp#L422) | function | Scalar math, complex & frames | Fused inverse Park-Clarke transform (dq → abc) |
| [`inverse_park_transform`](inc/damp/math/transforms.hpp#L381) | function | Scalar math, complex & frames | Inverse Park transform (dq → αβ) |
| [`inverse_symmetrical_components`](inc/damp/math/transforms.hpp#L570) | function | Scalar math, complex & frames | Inverse symmetrical-component transform (012 → abc) |
| [`iq_from_torque`](inc/damp/motor/spm.hpp#L70) | function | Motor control pack (if present) | i_q = T_e / Kₜ for id = 0 |
| [`is_constant_evaluated`](inc/damp/backend.hpp#L156) | function | Core, configuration & backend vocabulary | True during constant evaluation; the builtin needs no header |
| [`is_controllable`](inc/damp/design/stability.hpp#L446) | function | Design-time synthesis (not PWM-rate) | Check if a system is controllable |
| [`is_fault`](inc/damp/toolbox/conditioning.hpp#L332) | function | Embedded helpers (controls-adjacent utilities) | True for a wire fault (FaultLow/FaultHigh) — i.e. not a real reading at all |
| [`is_matrix_element`](inc/damp/matrix/matrix_traits.hpp#L39) | block | Linear algebra | True if T may be a Matrix element type |
| [`is_observable`](inc/damp/design/stability.hpp#L460) | function | Design-time synthesis (not PWM-rate) | Check if a system is observable |
| [`is_stabilizable`](inc/damp/design/stability.hpp#L483) | function | Design-time synthesis (not PWM-rate) | Check if (A, B) of sys is a stabilizable pair, in the domain of sys |
| [`is_stable`](inc/damp/design/stability.hpp#L504) | function | Design-time synthesis (not PWM-rate) | Stability of a state-space model, in the domain of sys |
| [`is_symmetric_or_hermitian`](inc/damp/matrix/decomposition.hpp#L90) | function | Linear algebra | True if A equals its transpose (real) or conjugate transpose (complex) |
| [`is_usable`](inc/damp/toolbox/conditioning.hpp#L506) | function | Embedded helpers (controls-adjacent utilities) | True when a channel reading may still be used (in-span or saturated) |
| [`is_valid`](inc/damp/toolbox/conditioning.hpp#L327) | function | Embedded helpers (controls-adjacent utilities) | True only for the in-span status |
| [`isfinite`](inc/damp/math/math.hpp#L500) | function | Scalar math, complex & frames | Finiteness test — false for NaN and ±∞ |
| [`iso_c_drive`](inc/damp/toolbox/io.hpp#L304) | function | Embedded helpers (controls-adjacent utilities) | ISO-C (single-stick, coordinated / curvature) drive |
| [`isstable`](inc/damp/matlab.hpp#L1367) | function | MATLAB®-style aliases (host) | Stability predicate for a state-space model (+1 more overload) |
| [`Jet`](inc/damp/trajectory/trajectory_types.hpp#L97) | block | Trajectory value types | Truncated jet of a scalar motion sample: d[k] = s⁽ᵏ⁾ |
| [`JordanBlock`](inc/damp/design/pole_placement.hpp#L488) | block | Design-time synthesis (not PWM-rate) | One Jordan mini-block of a desired closed-loop spectrum |
| [`JordanObjective`](inc/damp/design/pole_placement.hpp#L1198) | enum | Design-time synthesis (not PWM-rate) | How place_jordan spends the free parameter of the Klein–Moore family (Schmid et al.) |
| [`JordanOptions`](inc/damp/design/pole_placement.hpp#L1211) | block | Design-time synthesis (not PWM-rate) | Options for place_jordan |
| [`JordanPlaceResult`](inc/damp/design/pole_placement.hpp#L1226) | block | Design-time synthesis (not PWM-rate) | place_jordan result: the PlaceResult fields plus the achieved robustness metrics |
| [`joseph_update`](inc/damp/estimation/joseph.hpp#L40) | function | Observers & estimators | P ← (I − K H) P (I − K H)ᵀ + K R Kᵀ with K = P Hᵀ S⁻¹ and S = H P Hᵀ + R, for an H whose zeros and ones are known at compile time |
| [`kalman`](inc/damp/estimation/kalman.hpp#L800) | function | Observers & estimators | Steady-state Kalman filter design (+1 more overload) |
| [`KalmanFilter`](inc/damp/estimation/kalman.hpp#L1030) | block | Observers & estimators | The recursive Kalman filter: the covariance propagates each tick as its triangular factor, and the gain is formed from it (+2 more overloads) |
| [`KalmanFilterConfig`](inc/damp/estimation/kalman.hpp#L895) | block | Observers & estimators | Configuration of a KalmanFilter built from a model rather than from design::kalman |
| [`KalmanNoise`](inc/damp/estimation/kalman.hpp#L90) | block | Observers & estimators | The estimator's noise: process noise on the plant's noise inputs (through G), measurement noise on the sensors (through H) (+12 more overloads) |
| [`KalmanResult`](inc/damp/estimation/kalman.hpp#L227) | block | Observers & estimators | Steady-state Kalman filter design result |
| [`kalmd`](inc/damp/matlab.hpp#L994) | function | MATLAB®-style aliases (host) | Discrete Kalman estimator for a continuous plant and continuous noise (MATLAB®'s kalmd) |
| [`lag`](inc/damp/controllers/lead_lag.hpp#L198) | function | Design-time synthesis (not PWM-rate) | Design a lag compensator from desired low-frequency gain boost |
| [`lambda_tuning`](inc/damp/design/pid_design.hpp#L433) | function | Design-time synthesis (not PWM-rate) | Lambda tuning for FOPDT model |
| [`Latch`](inc/damp/toolbox/logic.hpp#L65) | block | Embedded helpers (controls-adjacent utilities) | Set/reset latch. @tparam SetDominant which input wins when both are asserted (default: set-dominant, e.g. a trip overriding a clear) |
| [`Law`](inc/damp/design/model_inverse.hpp#L950) | block | Design-time synthesis (not PWM-rate) | Law<kRes, T> of a `constexpr` model inverse: the same runtime on the design's matrices as constants of the type (+1 more overload) |
| [`Law`](inc/damp/estimation/eso.hpp#L460) | block | Observers & estimators | Law<kRes, T> of a design::eso initialized at start-up: the ESO it runs as. A runtime of the same shape assigned replaces it (+5 more overloads) |
| [`Law`](inc/damp/controllers/law.hpp#L60) | block | Runtime controllers | The runtime a design result deploys as, named by the result: `Law<kRes, float> law;` (+6 more overloads) |
| [`LDLT`](inc/damp/matrix/decomposition.hpp#L176) | block | Linear algebra | LDLᴴ factors with symmetric pivoting: P·A·Pᵀ = L·D·Lᴴ, D block diagonal with 1×1 and 2×2 blocks |
| [`ldlt`](inc/damp/matrix/decomposition.hpp#L483) | function | Linear algebra | Symmetric indefinite LDLᴴ with rook pivoting: P·A·Pᵀ = L·D·Lᴴ, D block diagonal (1×1 and 2×2) |
| [`ldlt_solve`](inc/damp/matrix/solve.hpp#L194) | function | Linear algebra | Solve AX = B via pivoted LDLᴴ (P·A·Pᵀ = L·D·Lᴴ, symmetric indefinite) |
| [`lead`](inc/damp/controllers/lead_lag.hpp#L145) | function | Design-time synthesis (not PWM-rate) | Design a lead compensator from desired phase boost at a target frequency |
| [`lead_lag`](inc/damp/controllers/lead_lag.hpp#L259) | function | Design-time synthesis (not PWM-rate) | Design a lead-lag compensator (cascade of lead + lag sections) |
| [`lead_lag_direct`](inc/damp/controllers/lead_lag.hpp#L294) | function | Design-time synthesis (not PWM-rate) | Direct lead-lag specification from zero/pole locations |
| [`LeadLagController`](inc/damp/controllers/lead_lag.hpp#L464) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its design, keeping the result's scalar |
| [`LeadLagController`](inc/damp/controllers/lead_lag.hpp#L325) | block | Runtime controllers | Discrete lead-lag compensator |
| [`LeadLagResult`](inc/damp/controllers/lead_lag.hpp#L64) | block | Design-time synthesis (not PWM-rate) | Lead-lag compensator design result |
| [`LeadLagSeriesResult`](inc/damp/controllers/lead_lag.hpp#L225) | block | Design-time synthesis (not PWM-rate) | Cascaded lead+lag design result (2nd-order StateSpace + success) |
| [`lerp`](inc/damp/toolbox/scaling.hpp#L33) | function | Embedded helpers (controls-adjacent utilities) | Linear interpolation between a and b by fraction t |
| [`lft`](inc/damp/systems/state_space.hpp#L558) | function | LTI systems (SS / TF / ZPK / discretize) | Lower linear fractional transformation: K closes the loop from P's last outputs to its last inputs |
| [`LIMIT`](inc/damp/toolbox/iec61131.hpp#L590) | function | Embedded helpers (controls-adjacent utilities) | LIMIT (IEC 61131-3 selection function): clamp in to [mn, mx] |
| [`linear_modulation_voltage`](inc/damp/motor/modulation.hpp#L125) | function | Motor control pack (if present) | Peak phase voltage at the linear SVPWM hexagon limit |
| [`linear_screw`](inc/damp/toolbox/actuator.hpp#L168) | function | Embedded helpers (controls-adjacent utilities) | Build a ServoAxis for a linear axis driven by a leadscrew/belt |
| [`LinearizationResult`](inc/damp/design/linearization.hpp#L45) | block | Design-time synthesis (not PWM-rate) | Result of nonlinear operating-point linearization |
| [`linearize`](inc/damp/design/linearization.hpp#L117) | function | Design-time synthesis (not PWM-rate) | Linearize a Dynamics / Measurement model about an operating point (+1 more overload) |
| [`linmod`](inc/damp/matlab.hpp#L589) | function | MATLAB®-style aliases (host) | MATLAB®-style nonlinear linearization about an operating point |
| [`linspace`](inc/damp/analysis/linspace.hpp#L36) | function | Frequency-domain analysis (host) | num evenly spaced values from start to end inclusive (+1 more overload) |
| [`locate_zero_crossing_exact`](inc/damp/simulation/hybrid.hpp#L334) | function | Simulation / SIL harness (host) | Locate a state zero-crossing of g(x) on an Exact LTI segment |
| [`log`](inc/damp/matrix/functions.hpp#L1029) | function | Linear algebra | Principal matrix logarithm via inverse scaling and squaring |
| [`log`](inc/damp/math/math.hpp#L249) | function | Scalar math, complex & frames | Natural logarithm |
| [`log10`](inc/damp/math/math.hpp#L431) | function | Scalar math, complex & frames | Base-10 logarithm, log10(x) = ln(x) / ln(10) |
| [`logm`](inc/damp/matlab.hpp#L708) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the principal matrix logarithm (mat::log) |
| [`LogPolicy`](inc/damp/simulation/hybrid.hpp#L58) | enum | Simulation / SIL harness (host) | How often to record samples along a hybrid trajectory |
| [`logspace`](inc/damp/matlab.hpp#L1229) | function | MATLAB®-style aliases (host) | n points from 10ᵃ to 10ᵇ, logarithmically spaced (decade exponents, as MATLAB®'s) |
| [`loop_at_a_time_disk_margin`](inc/damp/analysis/disk_margin.hpp#L197) | function | Frequency-domain analysis (host) | Loop-at-a-time disk margins of a square MIMO loop L under unity negative feedback |
| [`loop_metrics`](inc/damp/analysis/frequency.hpp#L485) | function | Frequency-domain analysis (host) | One-call loop analysis: compute L/S/T response and return compact metrics (+1 more overload) |
| [`loop_response`](inc/damp/analysis/frequency.hpp#L590) | function | Frequency-domain analysis (host) | Compute open-loop L, sensitivity S, complementary sensitivity T, and Nyquist data (+1 more overload) |
| [`LoopDiskMarginResult`](inc/damp/analysis/disk_margin.hpp#L55) | block | Frequency-domain analysis (host) | Loop-at-a-time disk margins of a square MIMO loop: one SISO disk margin per channel |
| [`LoopResponseResult`](inc/damp/analysis/frequency.hpp#L420) | block | Frequency-domain analysis (host) | Open-loop and closed-loop frequency response package |
| [`LoopSummary`](inc/damp/analysis/frequency.hpp#L453) | block | Frequency-domain analysis (host) | Compact loop summary metrics for quick stability/robustness checks |
| [`LowerTriangle`](inc/damp/matrix/views.hpp#L164) | block | Linear algebra | Lower triangular view of a square matrix |
| [`lowpass`](inc/damp/filters/lowpass.hpp#L97) | function | Filters & signal conditioning | Discretize a continuous low-pass TF (Tustin) into direct-form coefficients (+1 more overload) |
| [`LowPass`](inc/damp/filters/lowpass.hpp#L223) | block | Filters & signal conditioning | Nth-order low-pass filter (direct form I) (+1 more overload) |
| [`lowpass_1st`](inc/damp/filters/iir_design.hpp#L365) | function | Filters & signal conditioning | First-order low-pass filter design (continuous-time) (+1 more overload) |
| [`lowpass_2nd`](inc/damp/filters/iir_design.hpp#L479) | function | Filters & signal conditioning | Second-order low-pass filter design (continuous-time) (+1 more overload) |
| [`lowpass_sections`](inc/damp/filters/biquad.hpp#L462) | function | Filters & signal conditioning | Nth-order Butterworth low-pass at cutoff fc as second-order sections, any order |
| [`LowPassResult`](inc/damp/filters/lowpass.hpp#L55) | block | Filters & signal conditioning | Direct-form coefficients of an Nth-order low-pass, with success / status |
| [`lowshelf`](inc/damp/filters/iir_design.hpp#L878) | function | Filters & signal conditioning | Low-shelf EQ filter: boost or cut everything below fc |
| [`lqg`](inc/damp/matlab.hpp#L1109) | function | MATLAB®-style aliases (host) | Linear-Quadratic-Gaussian regulator design |
| [`lqg_from_parts`](inc/damp/controllers/lqg.hpp#L56) | function | Design-time synthesis (not PWM-rate) | Assemble an LQG design from separately computed Kalman and LQR designs |
| [`lqgreg`](inc/damp/matlab.hpp#L1125) | function | MATLAB®-style aliases (host) | Combine separate Kalman filter and LQR designs into an LQG controller |
| [`lqgtrack`](inc/damp/matlab.hpp#L1146) | function | MATLAB®-style aliases (host) | Linear-Quadratic-Gaussian design with integral action for tracking (+2 more overloads) |
| [`lqi`](inc/damp/matlab.hpp#L1070) | function | MATLAB®-style aliases (host) | Linear-Quadratic Integral design for tracking (+2 more overloads) |
| [`lqr`](inc/damp/matlab.hpp#L904) | function | MATLAB®-style aliases (host) | Continuous-time LQR design (MATLAB®'s lqr) |
| [`lqr_gain`](inc/damp/design/riccati.hpp#L1358) | function | Design-time synthesis (not PWM-rate) | Optimal LQR state-feedback gain from a Riccati solution |
| [`LQRCost`](inc/damp/design/riccati.hpp#L1581) | block | Design-time synthesis (not PWM-rate) | A quadratic cost on (x, u): state cost Q, input cost R, cross term N |
| [`lqrd`](inc/damp/matlab.hpp#L934) | function | MATLAB®-style aliases (host) | Sampled-data LQR from a continuous plant: ZOH dynamics + Van Loan cost (MATLAB®'s lqrd) (+1 more overload) |
| [`LQRResult`](inc/damp/design/riccati.hpp#L1518) | block | Design-time synthesis (not PWM-rate) | The LQR gain: design::state_feedback on a plant, discrete_lqr and continuous_lqr on given matrices |
| [`LQRWeights`](inc/damp/controllers/state_feedback_types.hpp#L135) | block | Runtime controllers | The cost of an LQ design (design::state_feedback): what to penalize and how much |
| [`lqry`](inc/damp/matlab.hpp#L1020) | function | MATLAB®-style aliases (host) | Output-weighted continuous LQR (state cost Q = Cᵀ Q_y C) (+1 more overload) |
| [`lsim`](inc/damp/analysis/time_response.hpp#L284) | function | Frequency-domain analysis (host) | Forced time response of a (MIMO) state-space system to an input signal (+1 more overload) |
| [`LsimInfo`](inc/damp/analysis/time_response.hpp#L513) | block | Frequency-domain analysis (host) | Transient characteristics of an arbitrary response signal |
| [`lsiminfo`](inc/damp/analysis/time_response.hpp#L533) | function | Frequency-domain analysis (host) | Compute transient characteristics from an output/time signal |
| [`lsimplot`](inc/damp/simulation/plot_plotly.hpp#L751) | function | Simulation / SIL harness (host) | Plot a forced (lsim) simulation, one trace per output |
| [`LsimResult`](inc/damp/analysis/time_response.hpp#L66) | block | Frequency-domain analysis (host) | Result of a single-trajectory simulation: time, output, and state history |
| [`lu_decomposition`](inc/damp/matrix/decomposition.hpp#L534) | function | Linear algebra | LU decomposition with partial pivoting |
| [`lu_solve`](inc/damp/matrix/solve.hpp#L238) | function | Linear algebra | Solve AX = B via LU with partial pivoting |
| [`Lut1D`](inc/damp/toolbox/lookup.hpp#L108) | block | Embedded helpers (controls-adjacent utilities) | 1-D interpolating lookup table over monotonic breakpoints |
| [`Lut2D`](inc/damp/toolbox/lookup.hpp#L200) | block | Embedded helpers (controls-adjacent utilities) | 2-D bilinear interpolating lookup table over a regular grid |
| [`Lut3D`](inc/damp/toolbox/lookup.hpp#L274) | block | Embedded helpers (controls-adjacent utilities) | 3-D trilinear interpolating lookup table over a regular grid |
| [`lut_segment`](inc/damp/toolbox/lookup.hpp#L62) | function | Embedded helpers (controls-adjacent utilities) | Index of the interpolation segment containing x |
| [`lyap`](inc/damp/design/lyapunov.hpp#L452) | function | Design-time synthesis (not PWM-rate) | Solve the continuous-time Lyapunov equation A X + X Aᵀ + Q = 0 |
| [`lyap`](inc/damp/matlab.hpp#L1424) | function | MATLAB®-style aliases (host) | MATLAB® alias for the continuous Lyapunov solve AX+XAᵀ+Q=0 |
| [`lyapchol`](inc/damp/matlab.hpp#L1444) | function | MATLAB®-style aliases (host) | MATLAB® alias for the Cholesky factor of the continuous Lyapunov solution (X = RᵀR) |
| [`lyapunov_cholesky`](inc/damp/design/lyapunov.hpp#L500) | function | Design-time synthesis (not PWM-rate) | Cholesky factor of the continuous Lyapunov solution: A X + X Aᵀ + B Bᵀ = 0, X = RᵀR |
| [`MadgwickFilter`](inc/damp/estimation/sensor_fusion.hpp#L133) | block | Observers & estimators | Madgwick gradient-descent AHRS filter |
| [`mag2db`](inc/damp/math/math.hpp#L535) | function | Scalar math, complex & frames | Magnitude to decibels, 20·log10(mag) |
| [`MahonyFilter`](inc/damp/estimation/sensor_fusion.hpp#L175) | block | Observers & estimators | Mahony nonlinear complementary filter with PI correction |
| [`make1DOF`](inc/damp/matlab.hpp#L305) | function | MATLAB®-style aliases (host) | Force 1-DOF setpoint weights on a PID design result (b=c=1) |
| [`make2DOF`](inc/damp/matlab.hpp#L321) | function | MATLAB®-style aliases (host) | Apply 2-DOF setpoint weights on a PID design result |
| [`margin`](inc/damp/matlab.hpp#L1278) | function | MATLAB®-style aliases (host) | Gain and phase margins of a SISO loop over a frequency grid |
| [`MarginResult`](inc/damp/matlab.hpp#L1257) | block | MATLAB®-style aliases (host) | Gain/phase margins and their crossover frequencies |
| [`margins_from_frf`](inc/damp/estimation/frequency_response.hpp#L840) | function | Observers & estimators | Gain / phase margins from a fixed on-target FRF table |
| [`Matrix`](inc/damp/matrix/core.hpp#L296) | block | Linear algebra | Fixed-size, stack-allocated matrix for linear algebra operations |
| [`MatrixLike`](inc/damp/matrix/matrix_traits.hpp#L95) | concept | Linear algebra | Concept for any type that provides 2D matrix-like element access |
| [`MatrixLikeOf`](inc/damp/matrix/matrix_traits.hpp#L105) | concept | Linear algebra | Concept for a MatrixLike type with specific dimensions |
| [`max`](inc/damp/backend.hpp#L339) | function | Core, configuration & backend vocabulary | The larger of a and b, by value; a tie returns a, and so does a NaN without -ffast-math (+1 more overload) |
| [`max_iq`](inc/damp/motor/spm.hpp#L112) | function | Motor control pack (if present) | Max positive i_q for an SPM (i_d = 0) at electrical speed |
| [`max_torque_at_speed`](inc/damp/motor/spm.hpp#L176) | function | Motor control pack (if present) | SPM max-torque point at a mechanical speed (i_d = 0 only) |
| [`maximal_invariant_set`](inc/damp/design/invariant_set.hpp#L249) | function | Design-time synthesis (not PWM-rate) | Maximal positively invariant set of x⁺ = A x under G₀ x ≤ h₀ |
| [`MeasuredDisturbance`](inc/damp/design/disturbance.hpp#L182) | block | Design-time synthesis (not PWM-rate) | A measured disturbance: a known input the controller does not command, entering the plant's state through Bd (a load torque, a feed rate, the ambient temperature), read each tick (+1 more overload) |
| [`Measurement`](inc/damp/concepts.hpp#L213) | concept | Core, configuration & backend vocabulary | Nonlinear measurement y = h(x, u): a callable returning y as a ColVec<NY, T> |
| [`MeasurementJacobian`](inc/damp/concepts.hpp#L226) | block | Core, configuration & backend vocabulary | Analytic Jacobians of a Measurement model at (x, u), returned by its jacobian(x, u) |
| [`mechanize_step`](inc/damp/estimation/ins_mechanization.hpp#L209) | function | Observers & estimators | One dead-reckoning step from a raw IMU sample (subtracts the state's biases) |
| [`mechanize_step_from_corrected`](inc/damp/estimation/ins_mechanization.hpp#L174) | function | Observers & estimators | One dead-reckoning step from an IMU sample (strapdown mechanization) |
| [`MedianFilter`](inc/damp/filters/median.hpp#L33) | block | Filters & signal conditioning | Sliding-window median filter — nonlinear spike/outlier rejection |
| [`mhe`](inc/damp/estimation/mhe.hpp#L217) | function | Observers & estimators | Synthesize a constrained moving-horizon estimator |
| [`MHE`](inc/damp/estimation/mhe.hpp#L342) | block | Observers & estimators | Runtime moving-horizon estimator (fixed per-tick iteration budget) (+1 more overload) |
| [`MHEConstraints`](inc/damp/estimation/mhe.hpp#L97) | block | Observers & estimators | Box constraints on the estimated states, applied to every window state |
| [`MHEResult`](inc/damp/estimation/mhe.hpp#L116) | block | Observers & estimators | Window QP data produced by mhe(), consumed by damp::MHE |
| [`middlebrook`](inc/damp/analysis/frequency.hpp#L1133) | function | Frequency-domain analysis (host) | Middlebrook stability analysis of a cascaded source and load from their impedances |
| [`MiddlebrookResult`](inc/damp/analysis/frequency.hpp#L939) | block | Frequency-domain analysis (host) | Result of Middlebrook minor loop gain analysis |
| [`min`](inc/damp/backend.hpp#L333) | function | Core, configuration & backend vocabulary | The smaller of a and b, by value; a tie returns a, and so does a NaN without -ffast-math (+1 more overload) |
| [`MinimalRealizationResult`](inc/damp/design/minreal.hpp#L88) | block | Design-time synthesis (not PWM-rate) | Minimal realization result |
| [`minmax`](inc/damp/backend.hpp#L388) | function | Core, configuration & backend vocabulary | Ordered {min, max} pair returned by value (+1 more overload) |
| [`minreal`](inc/damp/design/minreal.hpp#L310) | function | Design-time synthesis (not PWM-rate) | Minimal realization — cancel uncontrollable and unobservable modes |
| [`minreal`](inc/damp/matlab.hpp#L335) | function | MATLAB®-style aliases (host) | MATLAB® short alias for design::minreal, with its default tolerance √ε |
| [`minreal_zpk`](inc/damp/systems/zpk.hpp#L827) | function | LTI systems (SS / TF / ZPK / discretize) | Cancel matching pole-zero pairs on a ZPK model |
| [`mix`](inc/damp/systems/state_space.hpp#L806) | function | LTI systems (SS / TF / ZPK / discretize) | Fold plant outputs into one input: v = g u + Σ kᵢ yᵢ (+1 more overload) |
| [`mix_outputs`](inc/damp/systems/state_space.hpp#L883) | function | LTI systems (SS / TF / ZPK / discretize) | Linear combination of outputs: y = Σ αⱼ yⱼ |
| [`modal_form`](inc/damp/systems/realization.hpp#L287) | function | LTI systems (SS / TF / ZPK / discretize) | Real modal canonical form: A block diagonal, one block per mode |
| [`ModeExtractorConfig`](inc/damp/estimation/frequency_response.hpp#L473) | block | Observers & estimators | Configuration for FRF peak / valley extraction |
| [`ModeExtractorResult`](inc/damp/estimation/frequency_response.hpp#L490) | block | Observers & estimators | Result of reducing an FRF table to at most MaxModes peaks (resonances) |
| [`ModelReductionMethod`](inc/damp/design/model_reduction.hpp#L60) | enum | Design-time synthesis (not PWM-rate) | How balanced_reduction and eliminate_states eliminate states |
| [`modred`](inc/damp/matlab.hpp#L376) | function | MATLAB®-style aliases (host) | MATLAB® short alias for design::eliminate_states (0-based indices) |
| [`modulation_duty_cycles`](inc/damp/motor/modulation.hpp#L401) | function | Motor control pack (if present) | Carrier-based VSI duty cycles from an αβ voltage command |
| [`modulation_index`](inc/damp/motor/modulation.hpp#L163) | function | Motor control pack (if present) | Modulation index relative to the linear SVPWM circle |
| [`modulation_zero_sequence`](inc/damp/motor/modulation.hpp#L328) | function | Motor control pack (if present) | Zero-sequence for any PwmScheme |
| [`Modulator`](inc/damp/motor/modulation.hpp#L498) | block | Motor control pack (if present) | Thin scheme-holding wrapper for FOC / deploy paths |
| [`motor_constant`](inc/damp/motor/spm.hpp#L62) | function | Motor control pack (if present) | Motor constant Kₘ = Kₜ / √R [Nm/√W] |
| [`MovingAverage`](inc/damp/filters/moving_average.hpp#L67) | block | Filters & signal conditioning | Moving-average (boxcar) filter — also a DC-preserving harmonic-notch comb |
| [`mpc`](inc/damp/controllers/mpc.hpp#L2143) | function | Design-time synthesis (not PWM-rate) | Synthesize a constrained linear MPC (condensed dense QP, Δu form) on a discrete plant (+1 more overload) |
| [`MPC`](inc/damp/controllers/mpc.hpp#L3572) | function | Design-time synthesis (not PWM-rate) | MPC law{res}: every size of design::mpc's result, its scalar and its filter: on the measured state, or through its filter's FixedGainObserver (OFMPC with a disturbance spec); with a MeasuredDisturbance the QP and the target have NDQ = ND̂ + ND disturbance channels (+1 more overload) |
| [`MPC`](inc/damp/controllers/mpc.hpp#L3303) | block | Runtime controllers | Constrained linear MPC (design::mpc): on the measured state, or on the outputs through a Kalman filter, offset-free when the filter estimates a disturbance (OFMPC) |
| [`MPCArtifacts`](inc/damp/controllers/mpc.hpp#L571) | block | Runtime controllers | Condensed-QP data produced by mpc() (MPCResult::mpc), consumed by damp::MPC |
| [`MPCHorizonSuggestion`](inc/damp/controllers/mpc.hpp#L2207) | block | Runtime controllers | Advisory NP/NC horizon values from suggest_mpc_horizon() |
| [`MPCLimits`](inc/damp/controllers/mpc.hpp#L417) | block | Runtime controllers | The limits an MPC takes at run time (MATLAB®'s online constraints): u, Δu, the outputs its design's OutputLimits lists, the states its StateLimits lists and the disturbance estimate d̂, in physical units. A field left out is unbounded, and one the law does not have is empty |
| [`MPCResult`](inc/damp/controllers/mpc.hpp#L1710) | block | Design-time synthesis (not PWM-rate) | An MPC design: the condensed QP, the estimator when there is one, and the steady-state target the runtime poses the QP about; consumed by damp::MPC |
| [`MPCWeights`](inc/damp/controllers/mpc.hpp#L265) | block | Runtime controllers | Cost weights for mpc() |
| [`mstogi`](inc/damp/filters/sogi.hpp#L153) | function | Filters & signal conditioning | Mixed Second/Third-Order Generalized Integrator (MSTOGI) (+1 more overload) |
| [`MSTOGI`](inc/damp/filters/sogi.hpp#L299) | block | Filters & signal conditioning | Runtime MSTOGI with exact resonator and forward-Euler washout |
| [`multi_sine`](inc/damp/estimation/excitation/multi_sine.hpp#L137) | function | Observers & estimators | Build a multi-sine design payload from a configuration |
| [`MultiRateConfig`](inc/damp/simulation/multirate.hpp#L81) | block | Simulation / SIL harness (host) | Configuration for multi-rate closed-loop runs |
| [`MultiSine`](inc/damp/estimation/excitation/multi_sine.hpp#L211) | block | Observers & estimators | Sum-of-tones multi-sine runtime generator (+1 more overload) |
| [`MultiSineConfig`](inc/damp/estimation/excitation/multi_sine.hpp#L60) | block | Observers & estimators | Configuration for fixed-component multi-sine excitation |
| [`MultiSineResult`](inc/damp/estimation/excitation/multi_sine.hpp#L102) | block | Observers & estimators | Multi-sine design payload |
| [`MUX`](inc/damp/toolbox/iec61131.hpp#L609) | function | Embedded helpers (controls-adjacent utilities) | MUX (IEC 61131-3 multiplexer): select input k of N (0-based) |
| [`NavFrame`](inc/damp/estimation/ins_mechanization.hpp#L76) | enum | Observers & estimators | Local-level navigation frame for strapdown mechanization |
| [`nearbyint`](inc/damp/math/math.hpp#L396) | function | Scalar math, complex & frames | Round to nearest integer; ties to even (IEEE default / `FE_TONEAREST`) |
| [`ned_from_enu`](inc/damp/estimation/ins_mechanization.hpp#L105) | function | Observers & estimators | Map an ENU vector into NED (axis permute) |
| [`negative_sequence_ab`](inc/damp/filters/pll.hpp#L274) | function | Filters & signal conditioning | Instantaneous negative-sequence αβ from a quadrature signal pair |
| [`nichols`](inc/damp/analysis/frequency.hpp#L723) | function | Frequency-domain analysis (host) | Build Nichols points from existing Bode data (+2 more overloads) |
| [`nicholsplot`](inc/damp/simulation/plot_plotly.hpp#L944) | function | Simulation / SIL harness (host) | Plot a Nichols chart (open-loop phase vs magnitude) with M-circle grid |
| [`NicholsPoint`](inc/damp/analysis/frequency.hpp#L694) | block | Frequency-domain analysis (host) | Single-point Nichols chart sample (open-loop phase vs magnitude) |
| [`NicholsResult`](inc/damp/analysis/frequency.hpp#L706) | block | Frequency-domain analysis (host) | Nichols chart locus across a frequency sweep |
| [`NlmsFilter`](inc/damp/filters/fir.hpp#L573) | block | Filters & signal conditioning | Normalized LMS adaptive FIR filter |
| [`NoDisturbance`](inc/damp/design/disturbance.hpp#L40) | block | Design-time synthesis (not PWM-rate) | No disturbance model (the default): the plain filter on x |
| [`NoEstimator`](inc/damp/controllers/state_feedback_types.hpp#L46) | block | Runtime controllers | The estimator of a law that measures its state: none, an empty type, so a `[[no_unique_address]]` member of it takes no storage |
| [`norm`](inc/damp/matlab.hpp#L1526) | function | MATLAB®-style aliases (host) | MATLAB® alias for the H2 system norm norm(sys,2) |
| [`norm_h2`](inc/damp/analysis/norms.hpp#L81) | function | Frequency-domain analysis (host) | H2 norm of a state-space system |
| [`norm_hinf`](inc/damp/analysis/norms.hpp#L343) | function | Frequency-domain analysis (host) | H∞ norm of a stable state-space system, sup σ̄(G) over the imaginary axis or the unit circle, and its frequency |
| [`notch`](inc/damp/filters/iir_design.hpp#L666) | function | Filters & signal conditioning | Second-order band-reject (notch) filter design (continuous-time) (+1 more overload) |
| [`null`](inc/damp/matlab.hpp#L655) | function | MATLAB®-style aliases (host) | MATLAB® short alias for an orthonormal null-space basis: an optional, empty when the SVD fails |
| [`null_space`](inc/damp/matrix/svd.hpp#L478) | function | Linear algebra | Orthonormal basis for the null space {x : A·x = 0} via SVD |
| [`NullSpace`](inc/damp/matrix/svd.hpp#L458) | block | Linear algebra | Orthonormal basis for the null space (kernel) of a matrix |
| [`nyquist`](inc/damp/analysis/frequency.hpp#L502) | function | Frequency-domain analysis (host) | Compute Nyquist data for a SISO state-space system (+1 more overload) |
| [`nyquistplot`](inc/damp/simulation/plot_plotly.hpp#L823) | function | Simulation / SIL harness (host) | Plot a Nyquist locus with the -1 critical point marked |
| [`NyquistPoint`](inc/damp/analysis/frequency.hpp#L374) | block | Frequency-domain analysis (host) | Single-point Nyquist response data |
| [`NyquistResult`](inc/damp/analysis/frequency.hpp#L386) | block | Frequency-domain analysis (host) | Nyquist response data across a frequency sweep |
| [`observability_gramian`](inc/damp/design/stability.hpp#L409) | function | Design-time synthesis (not PWM-rate) | Continuous/discrete observability Gramian W_o |
| [`observability_matrix`](inc/damp/design/stability.hpp#L83) | function | Design-time synthesis (not PWM-rate) | Compute the observability matrix [C; CA; CA²; ...; CA^(N-1)] |
| [`observability_staircase`](inc/damp/systems/realization.hpp#L848) | function | LTI systems (SS / TF / ZPK / discretize) | Observability staircase form: the observable part split off by an orthogonal transform |
| [`observable_canonical_form`](inc/damp/systems/realization.hpp#L566) | function | LTI systems (SS / TF / ZPK / discretize) | Observable canonical form of a SISO model: ones on A's subdiagonal, −α in its last column, C = eₙᵀ |
| [`ObserverPoles`](inc/damp/controllers/state_feedback_types.hpp#L224) | function | Design-time synthesis (not PWM-rate) | ObserverPoles{p, …}: the count from the list, the scalar from its first entry |
| [`ObserverPoles`](inc/damp/controllers/state_feedback_types.hpp#L214) | block | Runtime controllers | The error poles of a placed filter (design::state_feedback): a trailing option, in place of KalmanNoise |
| [`obsv`](inc/damp/matlab.hpp#L403) | function | MATLAB®-style aliases (host) | MATLAB® short alias for observability_matrix (+1 more overload) |
| [`obsvf`](inc/damp/matlab.hpp#L552) | function | MATLAB®-style aliases (host) | MATLAB® obsvf(A, B, C): observability staircase form, its orthogonal transform and rank |
| [`OffDelayTimer`](inc/damp/toolbox/logic.hpp#L135) | block | Embedded helpers (controls-adjacent utilities) | Off-delay timer: output goes true immediately when in is true and stays true until in has been false continuously for delay |
| [`OnDelayTimer`](inc/damp/toolbox/logic.hpp#L99) | block | Embedded helpers (controls-adjacent utilities) | On-delay timer: output goes true once in has been held true continuously for delay; drops immediately when in goes false |
| [`one_norm`](inc/damp/matrix/functions.hpp#L69) | function | Linear algebra | One-norm ‖A‖₁: maximum absolute column sum |
| [`ones`](inc/damp/matlab.hpp#L619) | function | MATLAB®-style aliases (host) | N×N matrix of ones |
| [`optional`](inc/damp/backend.hpp#L186) | block | Core, configuration & backend vocabulary | The optional of a toolchain with no C++ standard headers (ETL_NO_STL) |
| [`OptionOf`](inc/damp/design/options.hpp#L60) | concept | Design-time synthesis (not PWM-rate) | An option of the kind Kind: a type that names it, `using option_kind = Kind;` |
| [`output_inverse`](inc/damp/design/model_inverse.hpp#L683) | function | Design-time synthesis (not PWM-rate) | The stable inverse of a discrete model: the state and input of its trajectory through a reference y on the controlled outputs, from y ₖ … y ₖ₊NP₋₁ (+1 more overload) |
| [`OutputDelay`](inc/damp/systems/discretization.hpp#L589) | block | LTI systems (SS / TF / ZPK / discretize) | A measurement delay of θ = seconds, in up to N delay states (0 ≤ θ ≤ N·Ts): the reading that arrives at sample k is y(t_k − θ). `discretize(sys, Ts, OutputDelay<N>{θ})` returns that model; design on it like any other. `OutputDelay{θ}` is N = 1 (+1 more overload) |
| [`OutputDisturbance`](inc/damp/design/disturbance.hpp#L75) | block | Design-time synthesis (not PWM-rate) | An output disturbance: a bias on each reading, B_d = 0 and C_d = I, one per output (+1 more overload) |
| [`OutputInverse`](inc/damp/design/model_inverse.hpp#L861) | block | Design-time synthesis (not PWM-rate) | The stable inverse at run time: the reference point (x , u ) from the output reference over the preview, seven matrix–vector products and the zero dynamics' state (+1 more overload) |
| [`OutputInverseResult`](inc/damp/design/model_inverse.hpp#L144) | block | Design-time synthesis (not PWM-rate) | The stable inverse of a discrete model over a preview of NP output samples |
| [`OutputLimits`](inc/damp/controllers/mpc.hpp#L308) | block | Runtime controllers | Limits on an MPC's outputs: among design::mpc's options |
| [`pade`](inc/damp/matlab.hpp#L1396) | function | MATLAB®-style aliases (host) | Padé approximation of order N (1 or 2) of the pure delay e^{−sT} |
| [`pade_delay_1st`](inc/damp/filters/pade.hpp#L33) | function | Filters & signal conditioning | First-order Pade approximation of time delay (+1 more overload) |
| [`pade_delay_2nd`](inc/damp/filters/pade.hpp#L85) | function | Filters & signal conditioning | Second-order Pade approximation of time delay (+1 more overload) |
| [`panel_color`](inc/damp/simulation/plot_plotly.hpp#L316) | function | Simulation / SIL harness (host) | Plotly default colorway entry for index i within a single subplot |
| [`panel_legend`](inc/damp/simulation/plot_plotly.hpp#L299) | function | Simulation / SIL harness (host) | Per-panel legend box just right of a stacked subplot band |
| [`panel_legend_id`](inc/damp/simulation/plot_plotly.hpp#L285) | function | Simulation / SIL harness (host) | Legend id for subplot row (1-based): "legend", "legend2", … |
| [`panel_x_domain`](inc/damp/simulation/plot_plotly.hpp#L275) | function | Simulation / SIL harness (host) | Horizontal domain for stacked multi-panel x-axes |
| [`parallel`](inc/damp/systems/state_space.hpp#L377) | function | LTI systems (SS / TF / ZPK / discretize) | Parallel connection (shared input, summed outputs) |
| [`ParameterDriftMonitor`](inc/damp/estimation/parameter_estimation.hpp#L251) | block | Observers & estimators | Residual-based drift detector: has the plant moved away from the model? |
| [`ParameterEstimator`](inc/damp/concepts.hpp#L164) | concept | Core, configuration & backend vocabulary | Online grey-box parameter estimator: physical parameters + gating |
| [`park_transform`](inc/damp/math/transforms.hpp#L354) | function | Scalar math, complex & frames | Park transform (αβ → dq) |
| [`pattern_row_times_factor_column`](inc/damp/estimation/sqrt_covariance.hpp#L344) | function | Observers & estimators | Row R of F times column J of the lower-triangular S: a term for each nonzero of the row at or below the diagonal of that column, and nothing else |
| [`PatternEntry`](inc/damp/matrix/pattern_matrix.hpp#L22) | enum | Linear algebra | What an entry of a sparsity pattern is known to be at compile time |
| [`PatternMatrix`](inc/damp/matrix/pattern_matrix.hpp#L226) | block | Linear algebra | The runtime values of a sparsity pattern, and nothing else: an entry the pattern knows has no storage, and in a product it has no term (a 0) or no multiply (a 1) |
| [`PatternMeasurement`](inc/damp/concepts.hpp#L258) | concept | Core, configuration & backend vocabulary | A Measurement whose jacobian(x, u) returns C = ∂h/∂x as an NY × NX PatternMatrix: zeros and ones known at compile time, the other entries computed at (x, u) |
| [`peaking`](inc/damp/filters/iir_design.hpp#L853) | function | Filters & signal conditioning | Peaking (bell) EQ filter: boost or cut a band around f0 |
| [`perfect_row_scale`](inc/damp/estimation/sqrt_covariance.hpp#L509) | function | Observers & estimators | The scale of row m of a reading that has no noise of its own: 1/√(eps·‖f‖²), its noise taken at the rounding of its own prior variance ‖f‖² = c̃ P c̃ᵀ |
| [`Periodic`](inc/damp/toolbox/timing.hpp#L177) | block | Embedded helpers (controls-adjacent utilities) | Periodic trigger — fires once per elapsed period |
| [`phase_margin_from_damping_ratio`](inc/damp/design/pid_design.hpp#L892) | function | Design-time synthesis (not PWM-rate) | Approximate phase margin from damping ratio |
| [`phase_margin_unwrapped`](inc/damp/analysis/frequency.hpp#L169) | function | Frequency-domain analysis (host) | Find the phase margin using the unwrapped phase trajectory |
| [`PhaseUnwrapper`](inc/damp/encoder/unwrap.hpp#L40) | block | Embedded helpers (controls-adjacent utilities) | Counts whole periods from successive fine phases (no hardware counter) |
| [`pid`](inc/damp/controllers/pid.hpp#L380) | function | Design-time synthesis (not PWM-rate) | Continuous PID gains Kₚ, Kᵢ, K_d with every other knob at its default |
| [`pid`](inc/damp/matlab.hpp#L250) | function | MATLAB®-style aliases (host) | MATLAB®-style parallel-form continuous PID constructor |
| [`pid_pole_placement`](inc/damp/design/pid_design.hpp#L964) | function | Design-time synthesis (not PWM-rate) | PI gains that place the closed-loop poles of a first-order plant b₀/(a₁ s + a₀) (+4 more overloads) |
| [`PIDController`](inc/damp/controllers/pid.hpp#L1004) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its gains, keeping their scalar: PIDController pid{gains} is a PID; through an alias the mode comes from the alias: PIController pi{gains}; |
| [`PIDController`](inc/damp/controllers/pid.hpp#L549) | block | Runtime controllers | Fixed-rate discrete 2-DOF PID (canonical runtime) (+3 more overloads) |
| [`PIDGains`](inc/damp/controllers/pid.hpp#L150) | block | Runtime controllers | 2-DOF continuous-time PID gains, limits and setpoint weights |
| [`PIDMode`](inc/damp/controllers/pid.hpp#L393) | enum | Runtime controllers | Compile-time selection of the PID control-law structure |
| [`PIDPolePlacementSpec`](inc/damp/design/pid_design.hpp#L916) | block | Design-time synthesis (not PWM-rate) | Closed-loop targets for the continuous pid_pole_placement forms |
| [`PIDRuntimeMode`](inc/damp/controllers/pid.hpp#L418) | enum | Runtime controllers | Runtime operating mode for PIDController / ContinuousPID |
| [`pidstd`](inc/damp/matlab.hpp#L273) | function | MATLAB®-style aliases (host) | Standard-form continuous PID constructor (1-DOF) |
| [`pidstd2`](inc/damp/matlab.hpp#L292) | function | MATLAB®-style aliases (host) | Standard-form continuous 2-DOF PID constructor |
| [`pidtune`](inc/damp/design/pid_design.hpp#L811) | function | Design-time synthesis (not PWM-rate) | Plant-aware PID tune at a crossover |
| [`pidtune`](inc/damp/matlab.hpp#L1194) | function | MATLAB®-style aliases (host) | Plant-aware PID tune at crossover ωc (MATLAB® pidtune spelling): 60° PID |
| [`PIDTuneSpec`](inc/damp/design/pid_design.hpp#L453) | block | Design-time synthesis (not PWM-rate) | Knobs for plant-aware pidtune |
| [`PiecewiseLTI`](inc/damp/simulation/hybrid.hpp#L93) | block | Simulation / SIL harness (host) | Plant: NModes continuous LTI pieces (same state/input/output dimensions) |
| [`pinv`](inc/damp/matlab.hpp#L642) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the Moore–Penrose pseudoinverse: an optional, empty when the SVD fails |
| [`place`](inc/damp/design/pole_placement.hpp#L1162) | function | Design-time synthesis (not PWM-rate) | Robust pole placement (KNV): the gain K with eig(A − B·K) at poles, in the domain of sys (+1 more overload) |
| [`place`](inc/damp/matlab.hpp#L885) | function | MATLAB®-style aliases (host) | Robust multi-input pole placement (MATLAB®'s place) |
| [`place_jordan`](inc/damp/design/pole_placement.hpp#L1686) | function | Design-time synthesis (not PWM-rate) | Exact pole placement with an arbitrary Jordan structure (Schmid–Ntogramatzidis–Nguyen–Pandey / Klein–Moore parametric form) (+1 more overload) |
| [`place_observer`](inc/damp/estimation/fixed_gain_observer.hpp#L307) | function | Observers & estimators | Design a Luenberger observer by robust pole placement, in the domain of sys (+1 more overload) |
| [`PlaceObserverResult`](inc/damp/estimation/fixed_gain_observer.hpp#L117) | block | Observers & estimators | Luenberger observer design result: the plant, the gains L and M, and the achieved error poles |
| [`PlaceResult`](inc/damp/design/pole_placement.hpp#L1079) | block | Design-time synthesis (not PWM-rate) | Pole-placement design result: the gain K, the achieved closed-loop poles e, and Ts |
| [`PlantDisturbance`](inc/damp/design/disturbance.hpp#L99) | block | Design-time synthesis (not PWM-rate) | The plant's own disturbance channel: B_d = G and C_d = H, one disturbance per column of G (H needs as many columns) (+1 more overload) |
| [`plot_bode`](inc/damp/simulation/plot_plotly.hpp#L424) | function | Simulation / SIL harness (host) | Plot Bode magnitude and phase as subplots |
| [`plot_line`](inc/damp/simulation/plot_plotly.hpp#L493) | function | Simulation / SIL harness (host) | Simple line plot of time vs value |
| [`plot_simulation`](inc/damp/simulation/plot_plotly.hpp#L345) | function | Simulation / SIL harness (host) | Plot simulation results with subplots for states, outputs, and inputs |
| [`plot_step`](inc/damp/simulation/plot_plotly.hpp#L623) | function | Simulation / SIL harness (host) | Plot step response data |
| [`plot_xy`](inc/damp/simulation/plot_plotly.hpp#L558) | function | Simulation / SIL harness (host) | Multi-series planar scatter (true 2-D path / locus) (+1 more overload) |
| [`plots_root_for`](inc/damp/simulation/plot_plotly.hpp#L120) | function | Simulation / SIL harness (host) | Walk up from dir until a directory named "plots" is found |
| [`pole`](inc/damp/matlab.hpp#L1241) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the open-loop poles of a system |
| [`PoleInfo`](inc/damp/analysis/poles.hpp#L52) | block | Frequency-domain analysis (host) | Natural frequency and damping ratio for each pole |
| [`Poles`](inc/damp/controllers/state_feedback_types.hpp#L191) | function | Design-time synthesis (not PWM-rate) | Poles{p, …}: the count from the list, the scalar from its first entry |
| [`poles`](inc/damp/analysis/poles.hpp#L40) | function | Frequency-domain analysis (host) | Compute open-loop poles (eigenvalues of A matrix) |
| [`Poles`](inc/damp/controllers/state_feedback_types.hpp#L181) | block | Runtime controllers | The closed-loop poles of a placed design (design::state_feedback): the weights' slot, in place of LQRWeights |
| [`PoleZeroMap`](inc/damp/analysis/poles.hpp#L98) | block | Frequency-domain analysis (host) | Poles and zeros of a system, for pole-zero plotting |
| [`poly_horner`](inc/damp/toolbox/scaling.hpp#L136) | function | Embedded helpers (controls-adjacent utilities) | Evaluate a polynomial at x by Horner's method |
| [`poly_roots`](inc/damp/analysis/poles.hpp#L146) | function | Frequency-domain analysis (host) | Roots of a polynomial given highest power first (MATLAB® `roots`) |
| [`poly_roots`](inc/damp/systems/zpk.hpp#L122) | function | LTI systems (SS / TF / ZPK / discretize) | Compute roots of a real polynomial given highest power first |
| [`PolyRootsResult`](inc/damp/systems/zpk.hpp#L91) | block | LTI systems (SS / TF / ZPK / discretize) | Roots of a real polynomial given highest power first |
| [`Pose`](inc/damp/kinematics/pose.hpp#L61) | block | Kinematics / pose | Rigid-body pose: a translation and an orientation (unit quaternion) |
| [`positive_sequence_ab`](inc/damp/filters/pll.hpp#L253) | function | Filters & signal conditioning | Instantaneous positive-sequence αβ from a quadrature signal pair |
| [`pow`](inc/damp/matrix/functions.hpp#L1132) | function | Linear algebra | Integer matrix power via binary exponentiation (+1 more overload) |
| [`pow`](inc/damp/math/math.hpp#L279) | function | Scalar math, complex & frames | Power function with a real exponent, base^exponent (+1 more overload) |
| [`pr`](inc/damp/controllers/pr.hpp#L218) | function | Design-time synthesis (not PWM-rate) | Design a Proportional-Resonant controller |
| [`pr_harmonics`](inc/damp/controllers/pr.hpp#L256) | function | Design-time synthesis (not PWM-rate) | Design multiple-harmonic PR controller gains |
| [`prbs`](inc/damp/estimation/excitation/prbs.hpp#L163) | function | Observers & estimators | Build a PRBS design payload from a configuration |
| [`PRBS`](inc/damp/estimation/excitation/prbs.hpp#L186) | block | Observers & estimators | Maximal-length PRBS runtime generator (+1 more overload) |
| [`PRBSConfig`](inc/damp/estimation/excitation/prbs.hpp#L79) | block | Observers & estimators | Configuration for maximal-length pseudo-random binary excitation |
| [`PRBSResult`](inc/damp/estimation/excitation/prbs.hpp#L123) | block | Observers & estimators | PRBS design payload |
| [`PRController`](inc/damp/controllers/pr.hpp#L618) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its design, keeping the result's scalar |
| [`PRController`](inc/damp/controllers/pr.hpp#L436) | block | Runtime controllers | Discrete Proportional-Resonant Controller |
| [`prescale`](inc/damp/matlab.hpp#L482) | function | MATLAB®-style aliases (host) | MATLAB® prescale: diagonal state scaling that balances A, B and C |
| [`process_noise_factor`](inc/damp/estimation/sqrt_covariance.hpp#L305) | function | Observers & estimators | A factor Γ of the process noise on the state, Γ Γᵀ = G Q Gᵀ, for sqrt_time_update() |
| [`project_affine`](inc/damp/controllers/action_governor.hpp#L135) | function | Design-time synthesis (not PWM-rate) | Project u_des onto the polyhedron A u ≤ b (Euclidean) |
| [`project_box`](inc/damp/controllers/action_governor.hpp#L74) | function | Design-time synthesis (not PWM-rate) | Euclidean projection of u onto a Bounds box |
| [`PRResonator`](inc/damp/controllers/pr.hpp#L326) | block | Runtime controllers | Resonant section of PRController, run on its output y and the output's increment d |
| [`PRResult`](inc/damp/controllers/pr.hpp#L123) | block | Design-time synthesis (not PWM-rate) | Proportional-Resonant controller design result |
| [`PRSpec`](inc/damp/controllers/pr.hpp#L110) | block | Runtime controllers | The gains and frequencies of a proportional-resonant design; Ts is passed beside it |
| [`pseudo_inverse`](inc/damp/matrix/svd.hpp#L427) | function | Linear algebra | Moore–Penrose pseudoinverse A⁺ via SVD |
| [`PulseTimer`](inc/damp/toolbox/logic.hpp#L170) | block | Embedded helpers (controls-adjacent utilities) | Pulse timer (non-retriggerable): a rising edge of in emits a fixed |
| [`PwmScheme`](inc/damp/motor/modulation.hpp#L77) | enum | Motor control pack (if present) | Carrier-based three-phase VSI modulation scheme |
| [`pzmap`](inc/damp/analysis/poles.hpp#L169) | function | Frequency-domain analysis (host) | Pole-zero map of a SISO transfer function (MATLAB® `pzmap(tf)`) (+2 more overloads) |
| [`pzplot`](inc/damp/simulation/plot_plotly.hpp#L850) | function | Simulation / SIL harness (host) | Plot a pole-zero map on the complex plane (poles as ×, zeros as ○) |
| [`QPIterations`](inc/damp/design/qp.hpp#L114) | block | Design-time synthesis (not PWM-rate) | A QP's iteration budget per tick, among the options of design::mpc and design::mhe (left out: each design's default) |
| [`QPResult`](inc/damp/design/qp.hpp#L82) | block | Design-time synthesis (not PWM-rate) | Result of a dense QP solve |
| [`QPStatus`](inc/damp/design/qp.hpp#L49) | enum | Design-time synthesis (not PWM-rate) | Termination status of a QP solve |
| [`qr_decompose`](inc/damp/matrix/decomposition.hpp#L719) | function | Linear algebra | Thin QR via Householder reflections |
| [`QRDecomposition`](inc/damp/matrix/decomposition.hpp#L679) | block | Linear algebra | Thin QR factorization A = QR (Householder) |
| [`QuadMode`](inc/damp/encoder/quadrature.hpp#L45) | enum | Embedded helpers (controls-adjacent utilities) | Quadrature decode resolution (edges counted per A/B cycle) |
| [`quadprog`](inc/damp/matlab.hpp#L1465) | function | MATLAB®-style aliases (host) | MATLAB® alias for the dense inequality-constrained QP solve (+1 more overload) |
| [`quadratic_form`](inc/damp/matrix/core.hpp#L1380) | function | Linear algebra | Symmetric congruence (quadratic) form  S = M X Mᵀ (+1 more overload) |
| [`quadratic_form_in_place`](inc/damp/matrix/pattern_matrix.hpp#L352) | function | Linear algebra | S ← F S Fᵀ in place, for a symmetric S: every multiply the pattern can skip, skipped |
| [`QuadratureDecoder`](inc/damp/encoder/quadrature.hpp#L73) | block | Embedded helpers (controls-adjacent utilities) | Software A/B quadrature decoder with optional index |
| [`quantize`](inc/damp/simulation/quantize.hpp#L66) | function | Simulation / SIL harness (host) | ADC reading of x: x rounded to the nearest multiple of lsb, then clamped to [lo, hi] |
| [`Quaternion`](inc/damp/math/geometry.hpp#L359) | block | Scalar math, complex & frames | Unit quaternion rotation (w, x, y, z) (Hamilton product) |
| [`R_TRIG`](inc/damp/toolbox/iec61131.hpp#L116) | block | Embedded helpers (controls-adjacent utilities) | R_TRIG (Rising Edge Trigger) |
| [`rad2deg`](inc/damp/math/math.hpp#L560) | function | Scalar math, complex & frames | Radians to degrees, rad·180/π |
| [`ramp`](inc/damp/estimation/excitation/ramp.hpp#L99) | function | Observers & estimators | Build a ramp design payload from a configuration |
| [`Ramp`](inc/damp/estimation/excitation/ramp.hpp#L115) | block | Observers & estimators | Rate-limited ramp runtime generator (+1 more overload) |
| [`RampConfig`](inc/damp/estimation/excitation/ramp.hpp#L31) | block | Observers & estimators | Configuration for a slew-rate-limited ramp excitation |
| [`RampResult`](inc/damp/estimation/excitation/ramp.hpp#L63) | block | Observers & estimators | Ramp design payload |
| [`RangeMonitor`](inc/damp/toolbox/conditioning.hpp#L401) | block | Embedded helpers (controls-adjacent utilities) | Analog-input range/fault monitor (NAMUR NE43 pattern) |
| [`rank`](inc/damp/design/stability.hpp#L424) | function | Design-time synthesis (not PWM-rate) | Numerical rank of a matrix (forwards to mat::rank) (+1 more overload) |
| [`rank`](inc/damp/matrix/functions.hpp#L236) | function | Linear algebra | Numerical rank: the number of singular values above rel_tol·σₘₐₓ |
| [`rank_from_svd`](inc/damp/matrix/svd.hpp#L391) | function | Linear algebra | Numerical rank from a precomputed SVD result |
| [`reading_noise_factor`](inc/damp/estimation/sqrt_covariance.hpp#L469) | function | Observers & estimators | The factor of a reading's noise covariance R (positive semidefinite; only its lower triangle is read), or an invalid one when R is indefinite or not finite |
| [`reading_row_reach`](inc/damp/estimation/sqrt_covariance.hpp#L480) | function | Observers & estimators | The size of the terms row m of C S is summed from, ‖\|C_m\|·\|S\|‖: what a perfect row's ‖f‖ is measured against, so that "zero" is the rounding of those terms and not a literal |
| [`reading_times_factor`](inc/damp/estimation/sqrt_covariance.hpp#L428) | function | Observers & estimators | CS = C S for a lower-triangular S and a dense C (+1 more overload) |
| [`ReadingNoiseFactor`](inc/damp/estimation/sqrt_covariance.hpp#L456) | block | Observers & estimators | A reading's noise covariance R factored for sqrt_measurement_update(): R = L Lᵀ, and whether R is a covariance at all |
| [`reduced_order_observer`](inc/damp/estimation/fixed_gain_observer.hpp#L517) | function | Observers & estimators | Design a reduced-order (Gopinath) observer by pole placement, in the domain of sys (+1 more overload) |
| [`ReducedOrderObserver`](inc/damp/estimation/fixed_gain_observer.hpp#L1087) | block | Observers & estimators | Reduced-order (Gopinath) state observer (runtime) (+1 more overload) |
| [`ReducedOrderObserverResult`](inc/damp/estimation/fixed_gain_observer.hpp#L355) | block | Observers & estimators | Reduced-order (Gopinath) observer design result |
| [`RedundantAnalogPair`](inc/damp/toolbox/io.hpp#L91) | block | Embedded helpers (controls-adjacent utilities) | Dual-channel analog input with same-slope or opposite-slope cross-check |
| [`reflect_rows_below`](inc/damp/estimation/sqrt_covariance.hpp#L187) | function | Observers & estimators | Apply row i's reflection of triangularize() to the rows below it: row k ← row k − s·v, with column i's sign turned when flip |
| [`reg`](inc/damp/matlab.hpp#L763) | function | MATLAB®-style aliases (host) | Form the observer-based compensator y → u from a plant, K and L |
| [`relay_autotune`](inc/damp/estimation/relay_autotune.hpp#L226) | function | Observers & estimators | Build a validated relay-autotune design payload |
| [`RelayAutotuneConfig`](inc/damp/estimation/relay_autotune.hpp#L123) | block | Observers & estimators | Configuration for the relay-feedback autotuning experiment |
| [`RelayAutotuneOutput`](inc/damp/estimation/relay_autotune.hpp#L258) | block | Observers & estimators | Per-tick output of RelayAutotuner::step |
| [`RelayAutotuner`](inc/damp/estimation/relay_autotune.hpp#L277) | block | Observers & estimators | Runtime relay-feedback autotuner (+1 more overload) |
| [`RelayAutotuneResult`](inc/damp/estimation/relay_autotune.hpp#L185) | block | Observers & estimators | Relay-autotuner design payload |
| [`RelayAutotuneStatus`](inc/damp/estimation/relay_autotune.hpp#L242) | enum | Observers & estimators | Lifecycle state of a RelayAutotuner |
| [`repetitive`](inc/damp/controllers/repetitive.hpp#L211) | function | Design-time synthesis (not PWM-rate) | Synthesize a repetitive controller with a scalar robustness filter Q |
| [`repetitive_binomial`](inc/damp/controllers/repetitive.hpp#L244) | function | Design-time synthesis (not PWM-rate) | Synthesize a repetitive controller with a binomial zero-phase FIR Q |
| [`RepetitiveConfig`](inc/damp/controllers/repetitive.hpp#L74) | block | Runtime controllers | Repetitive-controller tuning + period (with optional zero-phase FIR Q) |
| [`RepetitiveController`](inc/damp/controllers/repetitive.hpp#L288) | block | Runtime controllers | Plug-in repetitive controller runtime (fixed-size internal model) |
| [`RepetitiveResult`](inc/damp/controllers/repetitive.hpp#L137) | block | Design-time synthesis (not PWM-rate) | Design result for the repetitive controller |
| [`resample`](inc/damp/systems/discretization.hpp#L402) | function | LTI systems (SS / TF / ZPK / discretize) | Resample a ZOH-discrete model to a new sample time |
| [`rescale`](inc/damp/toolbox/scaling.hpp#L67) | function | Embedded helpers (controls-adjacent utilities) | Affine map of x from the input range to the output range |
| [`ResonantMode`](inc/damp/estimation/frequency_response.hpp#L446) | block | Observers & estimators | One extracted FRF extremum (peak or valley) |
| [`RisingEdge`](inc/damp/toolbox/logic.hpp#L33) | block | Embedded helpers (controls-adjacent utilities) | Rising-edge detector: true on the tick x goes false → true |
| [`RK23`](inc/damp/simulation/integrator.hpp#L933) | block | Simulation / SIL harness (host) | Bogacki-Shampine 2(3) adaptive integrator |
| [`RK3`](inc/damp/simulation/integrator.hpp#L1000) | block | Simulation / SIL harness (host) | Classical 3rd-order Runge-Kutta (RK3) integrator |
| [`RK4`](inc/damp/simulation/integrator.hpp#L725) | block | Simulation / SIL harness (host) | Classical 4th-order Runge-Kutta (RK4) integrator |
| [`rlocus`](inc/damp/analysis/poles.hpp#L289) | function | Frequency-domain analysis (host) | Root locus of a SISO state-space plant over an explicit gain grid (+1 more overload) |
| [`rlocusplot`](inc/damp/simulation/plot_plotly.hpp#L1042) | function | Simulation / SIL harness (host) | Plot a root locus (closed-loop poles vs gain) on the complex plane |
| [`rls_vector`](inc/damp/estimation/rls.hpp#L135) | function | Observers & estimators | Build vector RLS design payload |
| [`RlsConfig`](inc/damp/estimation/rls.hpp#L62) | block | Observers & estimators | Common RLS configuration |
| [`RlsVector`](inc/damp/estimation/rls.hpp#L150) | block | Observers & estimators | Vector runtime RLS estimator (NP parameters) (+1 more overload) |
| [`RlsVectorResult`](inc/damp/estimation/rls.hpp#L106) | block | Observers & estimators | Vector RLS design payload for N parameters |
| [`RlsVectorState`](inc/damp/estimation/rls.hpp#L92) | block | Observers & estimators | Vector RLS runtime state for N parameters |
| [`RobustExactDifferentiator`](inc/damp/filters/differentiator.hpp#L62) | block | Filters & signal conditioning | First-order robust exact differentiator (super-twisting differentiator) |
| [`RootLocusResult`](inc/damp/analysis/poles.hpp#L230) | block | Frequency-domain analysis (host) | Root-locus data: closed-loop poles along a gain grid |
| [`rotary_gearbox`](inc/damp/toolbox/actuator.hpp#L151) | function | Embedded helpers (controls-adjacent utilities) | Build a ServoAxis for a rotary joint behind a gearbox |
| [`RowVec`](inc/damp/matrix/rowvec.hpp#L27) | block | Linear algebra | Row vector specialization of Matrix<1, N, T> |
| [`RowView`](inc/damp/matrix/views.hpp#L202) | block | Linear algebra | Non-owning row view of a matrix |
| [`RS`](inc/damp/toolbox/iec61131.hpp#L89) | block | Embedded helpers (controls-adjacent utilities) | RS Latch (Reset-Set Latch) |
| [`SamePrecision`](inc/damp/matrix/matrix_traits.hpp#L84) | concept | Linear algebra | Two matrix element types mix in one operator without rounding: equal, or a real and the complex of its own precision |
| [`scalar_of`](inc/damp/design/detail/config.hpp#L64) | block | Design-time synthesis (not PWM-rate) | The scalar of a weight's value; double left out (+2 more overloads) |
| [`scaled_deadband`](inc/damp/toolbox/conditioning.hpp#L183) | function | Embedded helpers (controls-adjacent utilities) | Center dead zone that rescales the surviving range back to full span |
| [`schroeder_multi_sine`](inc/damp/estimation/excitation/multi_sine.hpp#L191) | function | Observers & estimators | Build a multi-sine design payload with Schroeder phases applied |
| [`second_order_sections`](inc/damp/filters/biquad.hpp#L372) | function | Filters & signal conditioning | Factor a discrete transfer function in z into second-order sections, for BiquadCascade (+1 more overload) |
| [`SecondOrderCoeffs`](inc/damp/filters/iir_design.hpp#L71) | block | Filters & signal conditioning | DSP coefficients for second-order IIR filter |
| [`SecondOrderSectionsResult`](inc/damp/filters/biquad.hpp#L40) | block | Filters & signal conditioning | An IIR filter as a cascade of second-order sections: the result BiquadCascade runs |
| [`SEL`](inc/damp/toolbox/iec61131.hpp#L598) | function | Embedded helpers (controls-adjacent utilities) | SEL (IEC 61131-3 binary selection): g ? in1 : in0 |
| [`select`](inc/damp/systems/state_space.hpp#L765) | function | LTI systems (SS / TF / ZPK / discretize) | SISO channel slice y_iy / u_iu |
| [`semidefinite_factor`](inc/damp/estimation/sqrt_covariance.hpp#L116) | function | Observers & estimators | The lower-triangular factor S of a positive semidefinite P, P = S Sᵀ, with a zero column where a pivot vanishes |
| [`SequenceComponents`](inc/damp/math/transforms.hpp#L513) | block | Scalar math, complex & frames | Symmetrical (sequence) components of a three-phase phasor set |
| [`series`](inc/damp/systems/state_space.hpp#L321) | function | LTI systems (SS / TF / ZPK / discretize) | Series connection: sys2 follows sys1 (u → sys1 → sys2 → y) |
| [`ServoAxis`](inc/damp/toolbox/actuator.hpp#L105) | block | Embedded helpers (controls-adjacent utilities) | One servoactuator transmission: SI joint unit ⟷ drive (motor) units |
| [`ServoBank`](inc/damp/toolbox/actuator.hpp#L231) | block | Embedded helpers (controls-adjacent utilities) | A bank of ServoAxis transmissions: maps a synchronized multi-axis TrajectoryState array straight to drive commands in one call |
| [`ServoCommand`](inc/damp/toolbox/actuator.hpp#L82) | block | Embedded helpers (controls-adjacent utilities) | A drive-native servoactuator setpoint: position, velocity, torque |
| [`sgn`](inc/damp/math/math.hpp#L448) | function | Scalar math, complex & frames | Sign function — −1 if val < 0, 1 if val > 0, 0 if val == 0 |
| [`sigma`](inc/damp/analysis/frequency.hpp#L836) | function | Frequency-domain analysis (host) | Singular-value frequency response of a (possibly MIMO) state-space system (+1 more overload) |
| [`sigmaplot`](inc/damp/simulation/plot_plotly.hpp#L999) | function | Simulation / SIL harness (host) | Plot singular-value frequency response (log frequency, dB) |
| [`SigmaPoint`](inc/damp/analysis/frequency.hpp#L779) | block | Frequency-domain analysis (host) | Singular values of G(jω) at one frequency |
| [`SigmaResult`](inc/damp/analysis/frequency.hpp#L793) | block | Frequency-domain analysis (host) | Singular-value frequency response over a grid |
| [`SignalRange`](inc/damp/toolbox/conditioning.hpp#L350) | block | Embedded helpers (controls-adjacent utilities) | The measuring span and wire-fault limits of an analog input (NAMUR NE43 pattern) |
| [`SignalSource`](inc/damp/concepts.hpp#L149) | concept | Core, configuration & backend vocabulary | Self-clocked signal source: u = step(), finished when done() |
| [`SignalStatus`](inc/damp/toolbox/conditioning.hpp#L318) | enum | Embedded helpers (controls-adjacent utilities) | Classification of an analog input against its valid/fault bands |
| [`simc`](inc/damp/design/pid_design.hpp#L372) | function | Design-time synthesis (not PWM-rate) | SIMC (Skogestad Internal Model Control) tuning for FOPDT models |
| [`similarity_transform`](inc/damp/systems/realization.hpp#L372) | function | LTI systems (SS / TF / ZPK / discretize) | Change the state basis to z = T·x: (T·A·T⁻¹, T·B, C·T⁻¹, D), and T·G for the noise input |
| [`simulate`](inc/damp/simulation/simulate.hpp#L159) | function | Simulation / SIL harness (host) | Simulate a nonlinear plant with a controller in closed loop |
| [`simulate_cached_zoh`](inc/damp/simulation/cached_zoh.hpp#L157) | function | Simulation / SIL harness (host) | Simulate a continuous LTI plant with fixed ZOH step and held input (+1 more overload) |
| [`simulate_cached_zoh_multirate`](inc/damp/simulation/cached_zoh.hpp#L236) | function | Simulation / SIL harness (host) | Multi-rate ZOH: held control samples, dense plant samples for plotting |
| [`simulate_discrete`](inc/damp/simulation/simulate.hpp#L529) | function | Simulation / SIL harness (host) | Simulate a discrete-time system with a controller |
| [`simulate_discrete_nonlinear`](inc/damp/simulation/simulate.hpp#L478) | function | Simulation / SIL harness (host) | Simulate a discrete-time nonlinear plant with a controller |
| [`simulate_lti`](inc/damp/simulation/simulate.hpp#L359) | function | Simulation / SIL harness (host) | Simulate a continuous LTI system with a controller |
| [`simulate_lti_siso`](inc/damp/simulation/simulate.hpp#L421) | function | Simulation / SIL harness (host) | SIL: continuous SISO LTI plant under a sampled SISO controller with constant reference |
| [`simulate_multirate`](inc/damp/simulation/multirate.hpp#L142) | function | Simulation / SIL harness (host) | Multi-rate closed-loop simulation with MultiRateConfig (+1 more overload) |
| [`simulate_piecewise_lti`](inc/damp/simulation/hybrid.hpp#L185) | function | Simulation / SIL harness (host) | Simulate a piecewise-LTI plant with Exact integration between events |
| [`simulate_sampled`](inc/damp/simulation/simulate.hpp#L232) | function | Simulation / SIL harness (host) | Simulate a continuous plant under a discrete (sampled) controller — multi-rate |
| [`simulate_state_feedback`](inc/damp/simulation/simulate.hpp#L303) | function | Simulation / SIL harness (host) | Simulate a nonlinear plant with state-feedback controller |
| [`simulate_two_rate`](inc/damp/simulation/multirate.hpp#L266) | function | Simulation / SIL harness (host) | Two-rate cascade with MultiRateConfig (+1 more overload) |
| [`simulate_two_rate_lti`](inc/damp/simulation/multirate.hpp#L415) | function | Simulation / SIL harness (host) | Two-rate cascade on a continuous LTI plant (A, B, C, D) |
| [`SimulationResult`](inc/damp/simulation/simulate.hpp#L128) | block | Simulation / SIL harness (host) | Result of a closed-loop simulation |
| [`sin`](inc/damp/matrix/functions.hpp#L1335) | function | Linear algebra | Matrix sine via scaling and double-angle reconstruction |
| [`sin`](inc/damp/math/math.hpp#L181) | function | Scalar math, complex & frames | Sine |
| [`sincos`](inc/damp/matrix/functions.hpp#L1319) | function | Linear algebra | Compute sin(A) and cos(A) together via scaling and double-angle reconstruction |
| [`sincos`](inc/damp/math/math.hpp#L199) | function | Scalar math, complex & frames | Combined sine and cosine, {sin(x), cos(x)} |
| [`SinCosCalibration`](inc/damp/encoder/sincos.hpp#L65) | block | Embedded helpers (controls-adjacent utilities) | Per-channel gain/offset trims for the sin/cos front end (hardware calibration) |
| [`SinCosSample`](inc/damp/encoder/sincos.hpp#L77) | block | Embedded helpers (controls-adjacent utilities) | One decoded sin/cos sample: the fine phase and the signal amplitude |
| [`SinglePhasePLL`](inc/damp/filters/pll.hpp#L67) | block | Filters & signal conditioning | Single-Phase PLL |
| [`sinh`](inc/damp/matrix/functions.hpp#L1370) | function | Linear algebra | Matrix hyperbolic sine sinh(A) = (exp(A) − exp(−A))/2 |
| [`siso_ref`](inc/damp/simulation/simulate.hpp#L115) | function | Simulation / SIL harness (host) | Build a SisoReferenceAdapter for a SISO `control(r,y)` controller |
| [`SisoController`](inc/damp/concepts.hpp#L89) | concept | Core, configuration & backend vocabulary | Scalar output-feedback controller: u = control(r, y), fixed rate |
| [`SisoReferenceAdapter`](inc/damp/simulation/simulate.hpp#L96) | block | Simulation / SIL harness (host) | Adapt a SISO `control(r, y)` controller for callables that want `u = f(y)` or `u = f(t, y)` (the older simulate contracts) |
| [`six_step_duty_cycles`](inc/damp/motor/modulation.hpp#L467) | function | Motor control pack (if present) | Classical six-step (full-wave) duty pattern from the αβ angle |
| [`six_step_fundamental_voltage`](inc/damp/motor/modulation.hpp#L144) | function | Motor control pack (if present) | Six-step (square-wave) fundamental peak phase voltage |
| [`SlewLimiter`](inc/damp/toolbox/conditioning.hpp#L218) | block | Embedded helpers (controls-adjacent utilities) | Slew-rate limiter: bound how fast the output may follow the target |
| [`smc`](inc/damp/controllers/smc.hpp#L91) | function | Design-time synthesis (not PWM-rate) | Bundle hand-picked SMC parameters into an SMCResult |
| [`SMCController`](inc/damp/controllers/smc.hpp#L280) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its design, keeping the result's scalar |
| [`SMCController`](inc/damp/controllers/smc.hpp#L183) | block | Runtime controllers | First-order sliding-mode controller (SMC) for a SISO plant |
| [`SMCResult`](inc/damp/controllers/smc.hpp#L37) | block | Design-time synthesis (not PWM-rate) | Tuning parameters for a first-order sliding-mode controller |
| [`smith_predictor_from_fopdt`](inc/damp/controllers/smith_predictor.hpp#L254) | function | Design-time synthesis (not PWM-rate) | Design a Smith predictor from FOPDT parameters |
| [`smith_predictor_from_pid`](inc/damp/controllers/smith_predictor.hpp#L287) | function | Design-time synthesis (not PWM-rate) | Pack a Smith predictor from an existing continuous PID and FOPDT plant |
| [`SmithPredictor`](inc/damp/controllers/smith_predictor.hpp#L329) | block | Runtime controllers | Smith predictor runtime: SISO primary + FO model + pure delay |
| [`SmithPredictorResult`](inc/damp/controllers/smith_predictor.hpp#L97) | block | Design-time synthesis (not PWM-rate) | Smith predictor design result (discrete PID + FO model + delay samples) |
| [`smooth_sign`](inc/damp/math/math.hpp#L461) | function | Scalar math, complex & frames | Smoothed sign: an odd quintic through ±width, exactly ±1 outside it |
| [`soft_sign`](inc/damp/toolbox/conditioning.hpp#L105) | function | Embedded helpers (controls-adjacent utilities) | Smooth unit direction x / √(x²+ε²) ∈ (−1,1) |
| [`sogi`](inc/damp/filters/sogi.hpp#L69) | function | Filters & signal conditioning | Second-Order Generalized Integrator (SOGI) design (+1 more overload) |
| [`SOGI`](inc/damp/filters/sogi.hpp#L226) | block | Filters & signal conditioning | Runtime SOGI quadrature-signal generator (live-frequency resonator) |
| [`SogiFll`](inc/damp/filters/sogi.hpp#L396) | block | Filters & signal conditioning | SOGI with a Frequency-Locked Loop — self-tuning single-tone tracker |
| [`solve`](inc/damp/matrix/solve.hpp#L288) | function | Linear algebra | Solve AX = B — Cholesky when Hermitian, else LU (+2 more overloads) |
| [`solve_miqp`](inc/damp/design/qp.hpp#L1762) | function | Design-time synthesis (not PWM-rate) | Mixed-integer QP by branch and bound over the active-set relaxation |
| [`solve_qp`](inc/damp/design/qp.hpp#L873) | function | Design-time synthesis (not PWM-rate) | Solve a strictly convex inequality-constrained QP (dual active-set) (+1 more overload) |
| [`solve_qp_admm`](inc/damp/design/qp.hpp#L1437) | function | Design-time synthesis (not PWM-rate) | One-shot ADMM QP solve (cold start) |
| [`solve_qp_interior_point`](inc/damp/design/qp.hpp#L1481) | function | Design-time synthesis (not PWM-rate) | Interior-point QP solver (primal-dual path following) |
| [`SolveResult`](inc/damp/simulation/solver.hpp#L86) | block | Simulation / SIL harness (host) | Result of an ODE solve operation |
| [`SparsityPattern`](inc/damp/matrix/pattern_matrix.hpp#L53) | block | Linear algebra | Which entries of an R × C matrix are exactly 0, exactly 1, or known only at run time |
| [`specific_force_at_rest`](inc/damp/estimation/ins_mechanization.hpp#L226) | function | Observers & estimators | Specific force [m/s²] a stationary IMU measures for the given orientation |
| [`SplineSurface`](inc/damp/toolbox/lookup.hpp#L393) | block | Embedded helpers (controls-adjacent utilities) | 2-D interpolating surface with smooth (Catmull-Rom spline) blending |
| [`SplitReal`](inc/damp/math/split_real.hpp#L47) | block | Scalar math, complex & frames | A real number split into an exact integer part and a floating-point remainder |
| [`sqrt`](inc/damp/matrix/functions.hpp#L983) | function | Linear algebra | Principal matrix square root by the real Schur method |
| [`sqrt`](inc/damp/math/complex.hpp#L337) | function | Scalar math, complex & frames | Compute complex square root (constexpr) (+1 more overload) |
| [`sqrt_measurement_update`](inc/damp/estimation/sqrt_covariance.hpp#L704) | function | Observers & estimators | Measurement update of the factor and the estimate with a reading of innovation ν, output matrix C and noise covariance R = L Lᵀ |
| [`sqrt_time_update`](inc/damp/estimation/sqrt_covariance.hpp#L414) | function | Observers & estimators | Time update of the factor: S ← the lower-triangular factor of A S Sᵀ Aᵀ + Γ Γᵀ |
| [`sqrtm`](inc/damp/matlab.hpp#L695) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the principal matrix square root (mat::sqrt) |
| [`square_matrix`](inc/damp/design/detail/config.hpp#L90) | function | Design-time synthesis (not PWM-rate) | A weight's value as an N × N Matrix in T (its size checked by the caller), or fallback when it was left out |
| [`SquareWeight`](inc/damp/design/detail/config.hpp#L86) | concept | Design-time synthesis (not PWM-rate) | A weight's value that is an N × N matrix: a square Matrix of that size, or a floating-point scalar for N = 1 |
| [`SR`](inc/damp/toolbox/iec61131.hpp#L61) | block | Embedded helpers (controls-adjacent utilities) | SR Latch (Set-dominant Set-Reset Latch) |
| [`ss`](inc/damp/matlab.hpp#L122) | function | MATLAB®-style aliases (host) | MATLAB®-style state-space model constructor |
| [`ss2ss`](inc/damp/matlab.hpp#L532) | function | MATLAB®-style aliases (host) | MATLAB® ss2ss(sys, T): the model in the state z = T·x |
| [`ss2tf`](inc/damp/systems/zpk.hpp#L806) | function | LTI systems (SS / TF / ZPK / discretize) | MATLAB®-style alias for SISO to_transfer_function |
| [`ss2zpk`](inc/damp/matlab.hpp#L228) | function | MATLAB®-style aliases (host) | MATLAB® short alias for SISO to_zpk (state space → zero-pole-gain) |
| [`ss_gain`](inc/damp/systems/state_space.hpp#L752) | function | LTI systems (SS / TF / ZPK / discretize) | Static gain y = k u (no states) |
| [`ssbal`](inc/damp/matlab.hpp#L491) | function | MATLAB®-style aliases (host) | MATLAB® ssbal: diagonal state balancing of [A B; C 0] |
| [`SsFeed`](inc/damp/systems/state_space.hpp#L741) | block | LTI systems (SS / TF / ZPK / discretize) | One term of a linear output mix: k · y_iy |
| [`stability_margin`](inc/damp/design/stability.hpp#L532) | function | Design-time synthesis (not PWM-rate) | Distance from the poles of sys to its domain's stability boundary |
| [`StaircaseFormResult`](inc/damp/systems/realization.hpp#L643) | block | LTI systems (SS / TF / ZPK / discretize) | A model in staircase coordinates, the orthogonal transform that took it there, and its rank |
| [`state_feedback`](inc/damp/controllers/state_feedback.hpp#L366) | function | Design-time synthesis (not PWM-rate) | The linear-quadratic state-feedback law (DARE on a discrete plant, CARE on a continuous one): runs as StateFeedback, an LQR, LQI, LQG, LQGI or OFLQG as its options say (+3 more overloads) |
| [`state_feedforward`](inc/damp/design/model_inverse.hpp#L100) | function | Design-time synthesis (not PWM-rate) | The plant-inversion feedforward of a discrete model: the input that moves x⁺ = A·x + B·u along a state reference, u ₖ the least-squares solution of B·u = x ₖ₊₁ − A·x ₖ |
| [`state_variable_filter`](inc/damp/estimation/state_variable_filter.hpp#L97) | function | Observers & estimators | Design the filter F(s) = ω³/(s + ω)³ sampled at Ts |
| [`StateEstimator`](inc/damp/concepts.hpp#L136) | concept | Core, configuration & backend vocabulary | State estimator: x̂ = step(y, u) — the fused per-tick form |
| [`StateFeedback`](inc/damp/controllers/state_feedback_runtime.hpp#L345) | function | Design-time synthesis (not PWM-rate) | StateFeedback law{res}: the sizes, scalar, filter and target of design::state_feedback's result (+1 more overload) |
| [`StateFeedback`](inc/damp/controllers/state_feedback_runtime.hpp#L131) | block | Runtime controllers | The runtime of every LQ design, u = u − K·(x̂ − x ) on the measured state or a filter's estimate: LQR, LQI, LQG, LQGI and OFLQG name its kinds |
| [`StateFeedbackLimits`](inc/damp/design/input_limits.hpp#L116) | block | Design-time synthesis (not PWM-rate) | The limits an LQ law takes at run time: u, its move, the integrators ξ, what the outer terms ask of a state and the disturbance estimate d̂, in physical units. A field left out is unbounded (the demand limit off), and one the law does not have is empty |
| [`StateFeedbackResult`](inc/damp/controllers/state_feedback_types.hpp#L427) | block | Design-time synthesis (not PWM-rate) | A state-feedback design: the gain, and the parts the law's options add |
| [`StateFeedforward`](inc/damp/design/model_inverse.hpp#L818) | block | Design-time synthesis (not PWM-rate) | The plant-inversion feedforward at run time: u = B⁺·(x ₖ₊₁ − x ₖ) − B⁺(A − I)·x ₖ, two matrix–vector products (+1 more overload) |
| [`StateFeedforwardResult`](inc/damp/design/model_inverse.hpp#L60) | block | Design-time synthesis (not PWM-rate) | The plant-inversion feedforward of a discrete model: u ₖ = B⁺·(x ₖ₊₁ − x ₖ) − B⁺(A − I)·x ₖ |
| [`StateLimits`](inc/damp/controllers/mpc.hpp#L388) | block | Runtime controllers | Limits on states an MPC does not track: among design::mpc's options |
| [`StateScaling`](inc/damp/systems/realization.hpp#L39) | block | LTI systems (SS / TF / ZPK / discretize) | A model in balanced state coordinates, and the diagonal scaling that took it there |
| [`StateSpace`](inc/damp/systems/state_space.hpp#L157) | block | LTI systems (SS / TF / ZPK / discretize) | State-space representation for linear time-invariant systems (discrete or continuous) (+16 more overloads) |
| [`StateSpaceResult`](inc/damp/design/model_result.hpp#L58) | block | Design-time synthesis (not PWM-rate) | A designed state-space model with its success / status |
| [`StateSpaceZPKResult`](inc/damp/systems/zpk.hpp#L656) | block | LTI systems (SS / TF / ZPK / discretize) | SISO state-space → ZPK conversion result with runtime zero count |
| [`StateTransformResult`](inc/damp/systems/realization.hpp#L125) | block | LTI systems (SS / TF / ZPK / discretize) | A model in new state coordinates z = transform·x, and the transform that took it there |
| [`StateVariableFilter`](inc/damp/estimation/state_variable_filter.hpp#L127) | block | Observers & estimators | State-variable filter over NS signals: filtered values and derivatives, no differentiation |
| [`StateVariableFilterResult`](inc/damp/estimation/state_variable_filter.hpp#L76) | block | Observers & estimators | Discretized state-variable filter (shared by every channel) |
| [`steady_state_target`](inc/damp/design/steady_state_target.hpp#L492) | function | Design-time synthesis (not PWM-rate) | The steady-state target of a law on a discrete plant: the maps from the estimated disturbance and the reference to [x_t; u_t] (+1 more overload) |
| [`steady_state_vdq`](inc/damp/motor/spm.hpp#L78) | function | Motor control pack (if present) | Steady-state Vdq (rotor frame) |
| [`steady_state_voltage_magnitude`](inc/damp/motor/spm.hpp#L91) | function | Motor control pack (if present) | \|Vdq\| at a steady-state operating point |
| [`SteadyStateTarget`](inc/damp/design/steady_state_target.hpp#L541) | block | Design-time synthesis (not PWM-rate) | The steady-state target at run time: (x_t, u_t) from the reference and the estimated disturbance, two matrix–vector products, and the re-anchoring of a FixedGainObserver on it (+1 more overload) |
| [`SteadyStateTargetResult`](inc/damp/design/steady_state_target.hpp#L92) | block | Design-time synthesis (not PWM-rate) | The steady-state target maps: [x_t; u_t] per unit disturbance and per unit reference |
| [`steer_map`](inc/damp/toolbox/io.hpp#L344) | function | Embedded helpers (controls-adjacent utilities) | Curve-driven steering map — arbitrary `(throttle, turn)` geometry |
| [`steinhart_hart`](inc/damp/toolbox/thermistor.hpp#L127) | function | Embedded helpers (controls-adjacent utilities) | Fit the Steinhart-Hart coefficients from three calibration points |
| [`step`](inc/damp/analysis/time_response.hpp#L164) | function | Frequency-domain analysis (host) | Step response of a (MIMO) state-space system (+1 more overload) |
| [`step`](inc/damp/estimation/excitation/step.hpp#L98) | function | Observers & estimators | Build a single-step design payload |
| [`Step`](inc/damp/estimation/excitation/step.hpp#L115) | block | Observers & estimators | Single step with optional \|dy/dt\| settle detection (+1 more overload) |
| [`step_train`](inc/damp/estimation/excitation/step_train.hpp#L98) | function | Observers & estimators | Build a step-train design payload from a configuration |
| [`StepConfig`](inc/damp/estimation/excitation/step.hpp#L34) | block | Observers & estimators | Configuration for a single step with optional settle detection |
| [`StepInfo`](inc/damp/analysis/time_response.hpp#L381) | block | Frequency-domain analysis (host) | Step-response characteristics of a single output signal |
| [`stepinfo`](inc/damp/analysis/time_response.hpp#L402) | function | Frequency-domain analysis (host) | Compute step-response characteristics from an output/time signal (+1 more overload) |
| [`stepped_sine`](inc/damp/estimation/excitation/stepped_sine.hpp#L109) | function | Observers & estimators | Build a stepped-sine design payload from a configuration |
| [`SteppedSine`](inc/damp/estimation/excitation/stepped_sine.hpp#L134) | block | Observers & estimators | Stepped-sine excitation — one pure tone at a time across a frequency table (+1 more overload) |
| [`SteppedSineConfig`](inc/damp/estimation/excitation/stepped_sine.hpp#L32) | block | Observers & estimators | Configuration for single-tone stepped-sine excitation |
| [`SteppedSineResult`](inc/damp/estimation/excitation/stepped_sine.hpp#L72) | block | Observers & estimators | Stepped-sine design payload |
| [`stepplot`](inc/damp/simulation/plot_plotly.hpp#L719) | function | Simulation / SIL harness (host) | Plot a step response, one trace per input/output pair |
| [`StepResult`](inc/damp/estimation/excitation/step.hpp#L66) | block | Observers & estimators | Single-step design payload |
| [`StepTrain`](inc/damp/estimation/excitation/step_train.hpp#L114) | block | Observers & estimators | Alternating +/- step train runtime generator (+1 more overload) |
| [`StepTrainConfig`](inc/damp/estimation/excitation/step_train.hpp#L30) | block | Observers & estimators | Configuration for alternating +/- step excitation |
| [`StepTrainResult`](inc/damp/estimation/excitation/step_train.hpp#L62) | block | Observers & estimators | Step-train design payload |
| [`Stopwatch`](inc/damp/toolbox/timing.hpp#L98) | block | Embedded helpers (controls-adjacent utilities) | Free-running elapsed-time accumulator |
| [`stsmc`](inc/damp/controllers/stsmc.hpp#L101) | function | Design-time synthesis (not PWM-rate) | Synthesize super-twisting gains from a disturbance-derivative bound |
| [`stsmc_gains`](inc/damp/controllers/stsmc.hpp#L145) | function | Design-time synthesis (not PWM-rate) | Super-twisting controller from gains you specify directly |
| [`STSMCController`](inc/damp/controllers/stsmc.hpp#L363) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its design, keeping the result's scalar |
| [`STSMCController`](inc/damp/controllers/stsmc.hpp#L207) | block | Runtime controllers | Super-twisting controller (second-order sliding mode) |
| [`STSMCResult`](inc/damp/controllers/stsmc.hpp#L35) | block | Design-time synthesis (not PWM-rate) | Super-twisting (second-order sliding-mode) controller design result |
| [`subtract`](inc/damp/systems/state_space.hpp#L639) | function | LTI systems (SS / TF / ZPK / discretize) | Differencing connection: outputs y = y₁ − y₂ |
| [`suggest_mpc_horizon`](inc/damp/controllers/mpc.hpp#L2231) | function | Design-time synthesis (not PWM-rate) | Suggest MPC horizons from an explicit settling-time target (+1 more overload) |
| [`summarize_loop_response`](inc/damp/analysis/frequency.hpp#L465) | function | Frequency-domain analysis (host) | Summarize loop_response() results into one compact metrics struct |
| [`svd`](inc/damp/matrix/svd.hpp#L357) | function | Linear algebra | Full singular value decomposition A = U·Σ·Vᴴ (one-sided Jacobi) |
| [`svd`](inc/damp/matlab.hpp#L633) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the singular value decomposition |
| [`SVDResult`](inc/damp/matrix/svd.hpp#L48) | block | Linear algebra | Result of a full singular value decomposition A = U·Σ·Vᴴ |
| [`svm_duty_cycles`](inc/damp/motor/modulation.hpp#L444) | function | Motor control pack (if present) | Space-vector PWM duty cycles from an αβ voltage command |
| [`SvmDuties`](inc/damp/motor/modulation.hpp#L101) | block | Motor control pack (if present) | Result of a duty-map: half-bridge duties plus an over-modulation flag |
| [`svpwm_zero_sequence`](inc/damp/motor/modulation.hpp#L198) | function | Motor control pack (if present) | Min-max zero-sequence injection for space-vector PWM |
| [`SwitchedController`](inc/damp/controllers/composition.hpp#L173) | block | Runtime controllers | Switch between normal, experiment, and backup SISO laws |
| [`SwitchMode`](inc/damp/controllers/composition.hpp#L146) | enum | Runtime controllers | Which path owns the plant command |
| [`symmetrical_components`](inc/damp/math/transforms.hpp#L541) | function | Scalar math, complex & frames | Forward symmetrical-component (Fortescue) transform (abc → 012) |
| [`SymplecticEuler`](inc/damp/simulation/integrator.hpp#L296) | block | Simulation / SIL harness (host) | Semi-implicit (symplectic) Euler for mechanical systems |
| [`tan`](inc/damp/math/math.hpp#L214) | function | Scalar math, complex & frames | Tangent |
| [`TargetPoint`](inc/damp/design/steady_state_target.hpp#L525) | block | Design-time synthesis (not PWM-rate) | A state and the input at it: a steady state (SteadyStateTarget), or a point of a reference trajectory (OutputInverse) |
| [`tf`](inc/damp/matlab.hpp#L86) | function | MATLAB®-style aliases (host) | MATLAB®-style transfer function constructor (+1 more overload) |
| [`tf2sos`](inc/damp/matlab.hpp#L183) | function | MATLAB®-style aliases (host) | MATLAB®-style transfer function → second-order sections, from MATLAB®'s z⁻¹ vectors (+1 more overload) |
| [`tf2zpk`](inc/damp/matlab.hpp#L158) | function | MATLAB®-style aliases (host) | MATLAB® short alias for to_zpk (transfer function → zero-pole-gain) |
| [`TFF`](inc/damp/toolbox/iec61131.hpp#L504) | block | Embedded helpers (controls-adjacent utilities) | T Flip-Flop (toggle on rising edge) |
| [`Thermistor`](inc/damp/toolbox/thermistor.hpp#L182) | block | Embedded helpers (controls-adjacent utilities) | NTC thermistor linearization (resistance → temperature) |
| [`thermistor_beta`](inc/damp/toolbox/thermistor.hpp#L81) | function | Embedded helpers (controls-adjacent utilities) | Fit NTC coefficients from the Beta-parameter model |
| [`ThermistorResult`](inc/damp/toolbox/thermistor.hpp#L42) | block | Embedded helpers (controls-adjacent utilities) | Fitted NTC coefficients in Steinhart-Hart form |
| [`thipwm_zero_sequence`](inc/damp/motor/modulation.hpp#L227) | function | Motor control pack (if present) | Third-harmonic injection (THIPWM) zero-sequence — 1/6 of fundamental |
| [`thiran`](inc/damp/matlab.hpp#L573) | function | MATLAB®-style aliases (host) | MATLAB® thiran: fractional-delay all-pass of order N for a delay tau at Ts |
| [`ThreePhasePLL`](inc/damp/filters/pll.hpp#L167) | block | Filters & signal conditioning | Synchronous-reference-frame (SRF) PLL for balanced three-phase input |
| [`Timeout`](inc/damp/toolbox/timing.hpp#L133) | block | Embedded helpers (controls-adjacent utilities) | One-shot timeout |
| [`TimeResponse`](inc/damp/analysis/time_response.hpp#L51) | block | Frequency-domain analysis (host) | Multi-channel time-domain response sampled on a time grid |
| [`times_factor_in_place`](inc/damp/estimation/sqrt_covariance.hpp#L361) | function | Observers & estimators | S ← A S in place for a lower-triangular S: column j of the product reads column j of S only, so one column of scratch is all it takes (+1 more overload) |
| [`times_transpose`](inc/damp/matrix/pattern_matrix.hpp#L328) | function | Linear algebra | A dense matrix times Hᵀ for a sparse H: A Hᵀ |
| [`to_coeffs`](inc/damp/filters/iir_design.hpp#L131) | function | Filters & signal conditioning | Convert StateSpace system to first-order DSP coefficients (+4 more overloads) |
| [`to_double_vector`](inc/damp/simulation/plot_plotly.hpp#L238) | function | Simulation / SIL harness (host) | Convert std::vector\<T\> to std::vector<double> |
| [`to_std_vector`](inc/damp/simulation/plot_plotly.hpp#L213) | function | Simulation / SIL harness (host) | Convert a ColVec<N,T> to std::vector<double> for plotlypp |
| [`to_transfer_function`](inc/damp/systems/zpk.hpp#L776) | function | LTI systems (SS / TF / ZPK / discretize) | SISO state-space → transfer function, size NX+1 / NX+1 |
| [`to_zpk`](inc/damp/systems/zpk.hpp#L463) | function | LTI systems (SS / TF / ZPK / discretize) | Convert a SISO transfer function to zero-pole-gain form (+1 more overload) |
| [`TOF`](inc/damp/toolbox/iec61131.hpp#L221) | block | Embedded helpers (controls-adjacent utilities) | TOF Timer (Timer Off Delay) |
| [`Toggle`](inc/damp/toolbox/logic.hpp#L243) | block | Embedded helpers (controls-adjacent utilities) | Toggle (T flip-flop): output flips on each rising edge of in |
| [`TON`](inc/damp/toolbox/iec61131.hpp#L181) | block | Embedded helpers (controls-adjacent utilities) | TON Timer (Timer On Delay) |
| [`Tone`](inc/damp/estimation/excitation/multi_sine.hpp#L27) | block | Observers & estimators | One sinusoidal component in a multi-sine excitation |
| [`torque`](inc/damp/motor/spm.hpp#L102) | function | Motor control pack (if present) | Electromagnetic torque of an SPM (i_d = 0): T_e = Kₜ i_q |
| [`torque_constant_from_flux`](inc/damp/motor/spm.hpp#L28) | function | Motor control pack (if present) | Kₜ = 1½ p λ [Nm/A] (amplitude-invariant) |
| [`torque_constant_from_Kv`](inc/damp/motor/spm.hpp#L45) | function | Motor control pack (if present) | Kₜ from hobby Kᵥ [RPM/V] (peak line-to-line / bus-volt sense) |
| [`torque_speed_envelope`](inc/damp/motor/spm.hpp#L207) | function | Motor control pack (if present) | Fixed table of SPM torque–speed envelope samples (host / LUT bake) |
| [`TorqueSpeedPoint`](inc/damp/motor/spm.hpp#L148) | block | Motor control pack (if present) | One sample of the SPM torque–speed envelope (MTPA i_d=0 branch) |
| [`TP`](inc/damp/toolbox/iec61131.hpp#L259) | block | Embedded helpers (controls-adjacent utilities) | TP Timer (Timer Pulse) |
| [`TrackedOutputs`](inc/damp/design/steady_state_target.hpp#L68) | block | Design-time synthesis (not PWM-rate) | The outputs a law without integrators or a disturbance model holds at r through its steady-state target (design::state_feedback and design::mpc, with or without a KalmanNoise, and design::steady_state_target, the target alone): `TrackedOutputs{{0}}`, in that order (the references r follow it). Left out, the law is a regulator and takes no r (+1 more overload) |
| [`TrajectoryBoundary`](inc/damp/trajectory/trajectory_types.hpp#L180) | block | Trajectory value types | Boundary conditions at one endpoint of a polynomial trajectory: a position and its time derivatives through jerk |
| [`TrajectoryLimits`](inc/damp/trajectory/trajectory_types.hpp#L47) | block | Trajectory value types | Asymmetric kinematic limits for a trapezoidal or S-curve motion profile |
| [`TrajectoryState`](inc/damp/trajectory/trajectory_types.hpp#L78) | block | Trajectory value types | A point on a motion profile: commanded position, velocity, acceleration |
| [`TransferFunction`](inc/damp/systems/transfer_function.hpp#L115) | block | LTI systems (SS / TF / ZPK / discretize) | SISO polynomial transfer function G(s) = num(s)/den(s) |
| [`TransferFunctionResult`](inc/damp/design/model_result.hpp#L34) | block | Design-time synthesis (not PWM-rate) | A designed transfer function with its success / status |
| [`transform`](inc/damp/matrix/core.hpp#L1527) | function | Linear algebra | Apply f to every entry, B(i, j) = f(A(i, j)): the matrix counterpart of std::transform |
| [`Transform4`](inc/damp/math/geometry.hpp#L782) | block | Scalar math, complex & frames | 4×4 homogeneous transform (SE(3) interop / DH export) |
| [`Translation3`](inc/damp/kinematics/pose.hpp#L38) | block | Kinematics / pose | A 3-D translation — a thin Vec3 with domain-named conveniences |
| [`transmission_zeros`](inc/damp/systems/zeros.hpp#L406) | function | LTI systems (SS / TF / ZPK / discretize) | The finite invariant zeros of a state-space model with any number of inputs and outputs |
| [`TransmissionZerosResult`](inc/damp/systems/zeros.hpp#L56) | block | LTI systems (SS / TF / ZPK / discretize) | The finite invariant zeros of a state-space model and their count |
| [`TransposeView`](inc/damp/matrix/views.hpp#L455) | block | Linear algebra | Non-owning transpose view of a matrix (zero-copy) |
| [`Trapezoidal`](inc/damp/simulation/integrator.hpp#L667) | block | Simulation / SIL harness (host) | Trapezoidal (Tustin) integrator |
| [`TRBDF2`](inc/damp/simulation/integrator.hpp#L542) | block | Simulation / SIL harness (host) | TR-BDF2 composite integrator with an embedded error estimate (the ode23tb scheme) |
| [`triangularize`](inc/damp/estimation/sqrt_covariance.hpp#L247) | function | Observers & estimators | [W E] Θ = [S 0] for an orthogonal Θ: W becomes the lower-triangular S with S Sᵀ = W Wᵀ + E Eᵀ, by Householder reflections, in place |
| [`two_norm`](inc/damp/matrix/functions.hpp#L106) | function | Linear algebra | Spectral norm ‖A‖₂ = σₘₐₓ(A) |
| [`two_point_cal`](inc/damp/toolbox/scaling.hpp#L110) | function | Embedded helpers (controls-adjacent utilities) | Fit an AffineCal through two `(raw, engineering)` points |
| [`TwoRateSimulationResult`](inc/damp/simulation/multirate.hpp#L108) | block | Simulation / SIL harness (host) | Result of a two-rate cascade simulation |
| [`tyreus_luyben`](inc/damp/design/pid_design.hpp#L216) | function | Design-time synthesis (not PWM-rate) | Tyreus-Luyben tuning from ultimate gain and ultimate period |
| [`tzero`](inc/damp/matlab.hpp#L562) | function | MATLAB®-style aliases (host) | MATLAB® tzero(sys): the finite invariant zeros of a model with any number of inputs and outputs |
| [`unbounded_bound`](inc/damp/design/qp.hpp#L67) | function | Design-time synthesis (not PWM-rate) | Sentinel bound treated as "no constraint" on that row |
| [`UnscentedKalmanFilter`](inc/damp/estimation/ukf.hpp#L130) | block | Observers & estimators | Unscented (sigma-point) Kalman Filter: continuous nonlinear dynamics, sampled measurements (+1 more overload) |
| [`UnscentedKalmanFilterConfig`](inc/damp/estimation/ukf.hpp#L79) | block | Observers & estimators | Configuration of an UnscentedKalmanFilter: its sample time, the process noise, the initial estimate and the unscented transform's tuning |
| [`UnscentedParams`](inc/damp/estimation/ukf.hpp#L57) | block | Observers & estimators | Tuning parameters for the scaled unscented transform |
| [`Unset`](inc/damp/design/detail/config.hpp#L31) | block | Design-time synthesis (not PWM-rate) | A config field left out: that feature is off, or the weight takes its default |
| [`unwrap_phase_deg`](inc/damp/analysis/frequency.hpp#L127) | function | Frequency-domain analysis (host) | Unwrap phase data in degrees to avoid +/-180 discontinuities |
| [`UpdateStatus`](inc/damp/estimation/sqrt_covariance.hpp#L51) | enum | Observers & estimators | What a recursive filter's last measurement update did, its last_status(): taken, or refused with the estimate and the covariance unchanged |
| [`UpperTriangle`](inc/damp/matrix/views.hpp#L123) | block | Linear algebra | Upper triangular view of a square matrix |
| [`upsample`](inc/damp/systems/discretization.hpp#L462) | function | LTI systems (SS / TF / ZPK / discretize) | Upsample a discrete model by an integer factor: sys(z) → sys(zᴸ) at Ts / L |
| [`upsample`](inc/damp/matlab.hpp#L472) | function | MATLAB®-style aliases (host) | MATLAB® upsample: sys(z) → sys(zᴸ) at Ts / L |
| [`vcat`](inc/damp/matrix/concat.hpp#L180) | function | Linear algebra | Vertical concatenation (same column count) |
| [`voltage_circle_radius`](inc/damp/motor/foc.hpp#L57) | function | Motor control pack (if present) | SVPWM voltage-circle radius Vₘₐₓ = m · Vdc / √3 |
| [`WarmStartActiveSetSolver`](inc/damp/design/qp.hpp#L1002) | block | Design-time synthesis (not PWM-rate) | Warm-started active-set solver policy — the damp::MPC default |
| [`weight_size`](inc/damp/design/detail/config.hpp#L52) | function | Design-time synthesis (not PWM-rate) | The size a weight's value gives: its rows for a Matrix, 1 for a scalar, 0 left out |
| [`whiten_reading`](inc/damp/estimation/sqrt_covariance.hpp#L542) | function | Observers & estimators | Whiten a reading in place with the factor L of its noise: each row of CS and of ν less the rows before it through L, then scaled to unit noise |
| [`with_schroeder_phases`](inc/damp/estimation/excitation/multi_sine.hpp#L168) | function | Observers & estimators | Assign Schroeder low-crest-factor phases to a multi-sine tone table |
| [`wrap`](inc/damp/math/math.hpp#L593) | function | Scalar math, complex & frames | Wrap x into the half-open interval [min, max) (period max − min) |
| [`wrap_pi`](inc/damp/math/math.hpp#L636) | function | Scalar math, complex & frames | Wrap an angle to [−π, π) |
| [`wrap_two_pi`](inc/damp/math/math.hpp#L654) | function | Scalar math, complex & frames | Wrap an angle to [0, 2π) |
| [`wrapped_delta`](inc/damp/encoder/quadrature.hpp#L39) | function | Embedded helpers (controls-adjacent utilities) | Signed difference between two unsigned counter readings, wrap-safe |
| [`write_html`](inc/damp/simulation/plot_plotly.hpp#L147) | function | Simulation / SIL harness (host) | Write a plotlypp figure to HTML with damp's shell (prefer over writeHtml) |
| [`write_hybrid_trace_npy`](inc/damp/simulation/npy_export.hpp#L101) | function | Simulation / SIL harness (host) | Pack HybridSimulationResult into (N, 2+NX+NU+NY) float64: t, mode, x..., u..., y |
| [`write_npy_f64`](inc/damp/simulation/npy_export.hpp#L42) | function | Simulation / SIL harness (host) | Write a C-order float64 array as a .npy file |
| [`XyPlotOpts`](inc/damp/simulation/plot_plotly.hpp#L537) | block | Simulation / SIL harness (host) | Optional axis framing for plot_xy |
| [`XySeries`](inc/damp/simulation/plot_plotly.hpp#L519) | block | Simulation / SIL harness (host) | One series for a planar (x, y) scatter / path plot |
| [`zero_sequence`](inc/damp/math/transforms.hpp#L249) | function | Scalar math, complex & frames | Zero-sequence (common-mode) scalar of a three-phase set |
| [`ziegler_nichols`](inc/damp/design/pid_design.hpp#L114) | function | Design-time synthesis (not PWM-rate) | Ziegler-Nichols tuning from ultimate gain and ultimate period |
| [`ziegler_nichols_step`](inc/damp/design/pid_design.hpp#L165) | function | Design-time synthesis (not PWM-rate) | Ziegler-Nichols step response method (reaction curve) |
| [`zp2sos`](inc/damp/matlab.hpp#L219) | function | MATLAB®-style aliases (host) | MATLAB®-style zero-pole-gain → second-order sections, H(z) = k·Π(z − zᵢ)/Π(z − pⱼ) |
| [`ZPK`](inc/damp/systems/zpk.hpp#L192) | block | LTI systems (SS / TF / ZPK / discretize) | Zero-pole-gain (ZPK) representation of a SISO LTI system |
| [`zpk`](inc/damp/matlab.hpp#L145) | function | MATLAB®-style aliases (host) | MATLAB®-style zero-pole-gain model constructor |
| [`zpk2tf`](inc/damp/matlab.hpp#L167) | function | MATLAB®-style aliases (host) | MATLAB® short alias for ZPK::to_transfer_function |
| [`ZPKResult`](inc/damp/systems/zpk.hpp#L427) | block | LTI systems (SS / TF / ZPK / discretize) | Result of converting a transfer function to zero-pole-gain form |
