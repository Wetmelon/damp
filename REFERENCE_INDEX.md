# API Reference — Alphabetical Index

Auto-generated from `@brief` doc comments in `inc/damp/`. Regenerate with `python tools/gen_reference.py`. Grouped-by-domain view: [REFERENCE.md](REFERENCE.md).

| Name | Kind | Domain | Description |
| ---- | ---- | ------ | ----------- |
| [`abs`](inc/damp/math/complex.hpp#L358) | function | Scalar math, complex & frames | Compute magnitude (absolute value) of a complex number (+1 more overload) |
| [`acker`](inc/damp/matlab.hpp#L561) | function | MATLAB®-style aliases (host) | Pole placement for state-feedback control |
| [`ackermann`](inc/damp/design/pole_placement.hpp#L1582) | function | Design-time synthesis (not PWM-rate) | Single-input pole placement via Ackermann's formula |
| [`acos`](inc/damp/math/math.hpp#L150) | function | Scalar math, complex & frames | Arccosine ∈ [0, π]. Input is clamped to [−1, 1] in both paths |
| [`ActionGovernor`](inc/damp/controllers/action_governor.hpp#L331) | block | Runtime controllers | Runtime action governor — QP projection of u_des onto A u ≤ b |
| [`ActiveSetSolver`](inc/damp/design/qp.hpp#L618) | block | Design-time synthesis (not PWM-rate) | Default QP solver policy: the Goldfarb–Idnani active-set solve_qp() |
| [`AdaptationGate`](inc/damp/estimation/parameter_estimation.hpp#L315) | block | Observers & estimators | Gated-adaptation policy: decide whether a parameter update is safe to apply |
| [`adaptive_solve`](inc/damp/simulation/solver.hpp#L774) | function | Simulation / SIL harness (host) | Adaptive-step solve in one call |
| [`AdaptiveLut1D`](inc/damp/toolbox/lookup.hpp#L427) | block | Embedded helpers (controls-adjacent utilities) | Adaptive 1-D lookup: fixed breakpoints, online-updated cell values |
| [`AdaptiveLut2D`](inc/damp/toolbox/lookup.hpp#L502) | block | Embedded helpers (controls-adjacent utilities) | Adaptive 2-D lookup: fixed grid, online-updated cells (nearest grid node) |
| [`AdaptiveOptions`](inc/damp/simulation/solver.hpp#L306) | block | Simulation / SIL harness (host) | Options for adaptive-step ODE integration |
| [`AdaptiveStepIntegrator`](inc/damp/simulation/integrator.hpp#L65) | concept | Simulation / SIL harness (host) | Integrator that reports a genuine embedded local-error estimate |
| [`AdaptiveStepSolver`](inc/damp/simulation/solver.hpp#L358) | block | Simulation / SIL harness (host) | Adaptive-step ODE solver (+1 more overload) |
| [`AdaptOutOfRange`](inc/damp/toolbox/lookup.hpp#L48) | enum | Embedded helpers (controls-adjacent utilities) | How AdaptiveLut1D / AdaptiveLut2D treat samples outside the breakpoint span |
| [`AdmmSettings`](inc/damp/design/qp.hpp#L713) | block | Design-time synthesis (not PWM-rate) | Tuning parameters for the ADMM QP solver |
| [`AdmmSolver`](inc/damp/design/qp.hpp#L752) | block | Design-time synthesis (not PWM-rate) | ADMM (OSQP-style) QP solver policy — warm-started, fixed-cost iterations |
| [`adrc`](inc/damp/controllers/adrc.hpp#L157) | function | Design-time synthesis (not PWM-rate) | Active Disturbance Rejection Control design |
| [`ADRCController`](inc/damp/controllers/adrc.hpp#L219) | block | Runtime controllers | Active Disturbance Rejection Control (ADRC) |
| [`ADRCResult`](inc/damp/controllers/adrc.hpp#L36) | block | Design-time synthesis (not PWM-rate) | Active Disturbance Rejection Control design result |
| [`AffineCal`](inc/damp/toolbox/scaling.hpp#L83) | block | Embedded helpers (controls-adjacent utilities) | Affine sensor calibration `y = gain·x + offset` |
| [`allmargin`](inc/damp/matlab.hpp#L926) | function | MATLAB®-style aliases (host) | Gain, phase, and delay margins of a SISO loop over a frequency grid |
| [`AllMarginResult`](inc/damp/matlab.hpp#L905) | block | MATLAB®-style aliases (host) | All classical margins including delay margin (superset of margin) |
| [`allpass`](inc/damp/filters/iir_design.hpp#L860) | function | Filters & signal conditioning | Second-order all-pass filter (Tustin of the analog prototype) |
| [`allpass_2nd_continuous`](inc/damp/filters/iir_design.hpp#L434) | function | Filters & signal conditioning | Second-order all-pass filter design (continuous-time) |
| [`AlphaBeta`](inc/damp/math/transforms.hpp#L152) | block | Scalar math, complex & frames | Alpha-beta (stationary-frame) component pair |
| [`amigo_kappa_tau`](inc/damp/design/pid_design.hpp#L240) | function | Design-time synthesis (not PWM-rate) | AMIGO PI from ultimate gain/period and static gain Kₛ |
| [`AnalogCrossMode`](inc/damp/toolbox/conditioning.hpp#L498) | enum | Embedded helpers (controls-adjacent utilities) | How two redundant analog channels should relate in engineering units |
| [`AnalogCrossResult`](inc/damp/toolbox/conditioning.hpp#L520) | block | Embedded helpers (controls-adjacent utilities) | Output of cross_check_analog |
| [`AnalogCrossStatus`](inc/damp/toolbox/conditioning.hpp#L506) | enum | Embedded helpers (controls-adjacent utilities) | Result of comparing two redundant analog channels |
| [`AnalogInput`](inc/damp/toolbox/io.hpp#L46) | block | Embedded helpers (controls-adjacent utilities) | A single analog input: range/fault check on the raw reading, then affine calibration to engineering units |
| [`AngleUnwrapper`](inc/damp/toolbox/encoder.hpp#L190) | block | Embedded helpers (controls-adjacent utilities) | Wrap-safe accumulator: successive wrapped absolute readings → continuous position |
| [`apply_wet_plotly_shell`](inc/damp/simulation/plot_plotly.hpp#L64) | function | Simulation / SIL harness (host) | Apply damp's full-window HTML shell to plotlypp-generated markup |
| [`arange`](inc/damp/analysis/linspace.hpp#L134) | function | Frequency-domain analysis (host) | Half-open arithmetic range [start, stop) with step step (+2 more overloads) |
| [`arcade_drive`](inc/damp/toolbox/io.hpp#L361) | function | Embedded helpers (controls-adjacent utilities) | Arcade (single-stick) drive mixer: x (steer) + y (throttle) → left/right |
| [`arg`](inc/damp/math/complex.hpp#L369) | function | Scalar math, complex & frames | Compute argument (phase angle) of a complex number |
| [`asin`](inc/damp/math/math.hpp#L131) | function | Scalar math, complex & frames | Arcsine ∈ [−π/2, π/2]. Input is clamped to [−1, 1] in both paths so behavior matches at compile and run time (std::asin would return NaN for \|x\| > 1) |
| [`atan`](inc/damp/math/math.hpp#L117) | function | Scalar math, complex & frames | Single-argument arctangent ∈ (−π/2, π/2) |
| [`atan2`](inc/damp/math/math.hpp#L105) | function | Scalar math, complex & frames | Two-argument arctangent, atan2(y, x) ∈ [−π, π] |
| [`AxisInput`](inc/damp/toolbox/io.hpp#L191) | block | Embedded helpers (controls-adjacent utilities) | Operator-axis conditioning chain (joystick / RC stick → command) |
| [`backward_substitute_transpose`](inc/damp/matrix/solve.hpp#L59) | function | Linear algebra | Backward substitution for Lᵀx = b (real) / Lᴴx = b layout |
| [`BackwardEuler`](inc/damp/simulation/integrator.hpp#L295) | block | Simulation / SIL harness (host) | Backward Euler integrator |
| [`balanced_realization`](inc/damp/design/model_reduction.hpp#L569) | function | Design-time synthesis (not PWM-rate) | Descriptive alias for balreal |
| [`balanced_truncation`](inc/damp/design/model_reduction.hpp#L586) | function | Design-time synthesis (not PWM-rate) | Descriptive alias for balred |
| [`BalancedRealizationResult`](inc/damp/design/model_reduction.hpp#L74) | block | Design-time synthesis (not PWM-rate) | Balanced realization result (Moore square-root method) |
| [`BalancedReductionResult`](inc/damp/design/model_reduction.hpp#L100) | block | Design-time synthesis (not PWM-rate) | Reduced-order model from balanced truncation or residualization |
| [`balreal`](inc/damp/design/model_reduction.hpp#L383) | function | Design-time synthesis (not PWM-rate) | Balanced realization via the square-root (Moore / Laub) method |
| [`balred`](inc/damp/design/model_reduction.hpp#L471) | function | Design-time synthesis (not PWM-rate) | Balanced truncation to NR states |
| [`bandpass`](inc/damp/filters/iir_design.hpp#L811) | function | Filters & signal conditioning | Second-order band-pass filter (constant 0 dB peak gain) |
| [`bandpass_continuous`](inc/damp/filters/iir_design.hpp#L378) | function | Filters & signal conditioning | Second-order band-pass filter design (continuous-time) |
| [`bandwidth`](inc/damp/matlab.hpp#L1031) | function | MATLAB®-style aliases (host) | -3 dB bandwidth of a SISO system over a frequency grid |
| [`bandwidth_from_settling_time`](inc/damp/design/pid_design.hpp#L819) | function | Design-time synthesis (not PWM-rate) | Map settling-time and damping-ratio targets to a bandwidth estimate |
| [`BDF2`](inc/damp/simulation/integrator.hpp#L366) | block | Simulation / SIL harness (host) | Backward Differentiation Formula 2 (BDF2) integrator |
| [`beta`](inc/damp/toolbox/thermistor.hpp#L77) | function | Embedded helpers (controls-adjacent utilities) | Fit NTC coefficients from the Beta-parameter model |
| [`Biquad`](inc/damp/filters/biquad.hpp#L34) | block | Filters & signal conditioning | Second-order IIR (biquad) section runtime |
| [`BiquadCascade`](inc/damp/filters/biquad.hpp#L101) | block | Filters & signal conditioning | Cascade of second-order sections (SOS) for higher-order IIR filters |
| [`BLINK`](inc/damp/toolbox/iec61131.hpp#L550) | block | Embedded helpers (controls-adjacent utilities) | BLINK (free-running square-wave / flasher) |
| [`blkdiag`](inc/damp/matlab.hpp#L356) | function | MATLAB®-style aliases (host) | Block diagonal matrix construction |
| [`Block`](inc/damp/matrix/block.hpp#L38) | block | Linear algebra | Block view (non-owning) into a parent matrix |
| [`bode`](inc/damp/analysis/frequency.hpp#L250) | function | Frequency-domain analysis (host) | Compute Bode plot data for a SISO state-space system (+2 more overloads) |
| [`bode_discrete`](inc/damp/analysis/frequency.hpp#L358) | function | Frequency-domain analysis (host) | Compute Bode plot data for a discrete-time SISO state-space system |
| [`bodemag`](inc/damp/simulation/plot_plotly.hpp#L796) | function | Simulation / SIL harness (host) | Plot a magnitude-only Bode diagram (log frequency, dB magnitude) |
| [`bodeplot`](inc/damp/simulation/plot_plotly.hpp#L782) | function | Simulation / SIL harness (host) | Plot magnitude and phase Bode subplots |
| [`BodeResult`](inc/damp/analysis/frequency.hpp#L60) | block | Frequency-domain analysis (host) | Bode plot data for a SISO system |
| [`Bounds`](inc/damp/toolbox/bounds.hpp#L43) | block | Embedded helpers (controls-adjacent utilities) | A per-channel closed-interval box constraint |
| [`BoxCommandFilter`](inc/damp/controllers/action_governor.hpp#L119) | block | Runtime controllers | Runtime box command filter — closed-form clamp per tick |
| [`build_lqg_analysis_models`](inc/damp/design/synthesis.hpp#L211) | function | Design-time synthesis (not PWM-rate) | Build analysis models from an LQG design |
| [`build_lqgi_analysis_models`](inc/damp/design/synthesis.hpp#L297) | function | Design-time synthesis (not PWM-rate) | Build analysis models from an LQGI servo design |
| [`build_lqi_analysis_models`](inc/damp/design/synthesis.hpp#L266) | function | Design-time synthesis (not PWM-rate) | Build analysis models from an LQI servo design |
| [`build_mpc_analysis_models`](inc/damp/controllers/mpc.hpp#L695) | function | Design-time synthesis (not PWM-rate) | Build the unconstrained-MPC LTI analysis models |
| [`butterworth_lowpass`](inc/damp/filters/iir_design.hpp#L464) | function | Filters & signal conditioning | Butterworth low-pass filter design (+1 more overload) |
| [`Button`](inc/damp/toolbox/io.hpp#L238) | block | Embedded helpers (controls-adjacent utilities) | Debounced momentary push-button with edge and long-press detection |
| [`c2d`](inc/damp/matlab.hpp#L303) | function | MATLAB®-style aliases (host) | MATLAB® interface function c2d to discretize a continuous-time state-space system (+1 more overload) |
| [`CachedZoh`](inc/damp/simulation/cached_zoh.hpp#L62) | block | Simulation / SIL harness (host) | Cached discrete ZOH maps Ad, Bd for step size h |
| [`canonical_phase_margin`](inc/damp/analysis/frequency.hpp#L148) | function | Frequency-domain analysis (host) | Normalize phase margin to (-180, 180] |
| [`care`](inc/damp/design/riccati.hpp#L763) | function | Design-time synthesis (not PWM-rate) | Solve the Continuous-time Algebraic Riccati Equation (CARE) |
| [`Cascade`](inc/damp/controllers/composition.hpp#L52) | block | Runtime controllers | Series cascade of two SISO controllers: outer → inner reference |
| [`catmull_1d`](inc/damp/toolbox/lookup.hpp#L339) | function | Embedded helpers (controls-adjacent utilities) | 1-D non-uniform Catmull-Rom (cubic Hermite) evaluation |
| [`cauer_thermal_ss_ambient`](inc/damp/toolbox/thermal.hpp#L135) | function | Embedded helpers (controls-adjacent utilities) | Continuous Cauer RC ladder in absolute temperature with ambient input |
| [`cauer_thermal_ss_ambient_mimo`](inc/damp/toolbox/thermal.hpp#L162) | function | Embedded helpers (controls-adjacent utilities) | Continuous Cauer RC ladder, absolute temps, all nodes as outputs |
| [`cbf_relative_degree_1`](inc/damp/controllers/action_governor.hpp#L299) | function | Design-time synthesis (not PWM-rate) | Build a relative-degree-1 CBF inequality row |
| [`CBFConstraint`](inc/damp/controllers/action_governor.hpp#L283) | block | Runtime controllers | One affine row from a relative-degree-1 CBF condition |
| [`cbrt`](inc/damp/math/math.hpp#L92) | function | Scalar math, complex & frames | Cube root (preserves sign for negative x) |
| [`ceil`](inc/damp/math/math.hpp#L329) | function | Scalar math, complex & frames | Ceiling — smallest integer ≥ x |
| [`chirp`](inc/damp/estimation/excitation/chirp.hpp#L115) | function | Observers & estimators | Build a chirp design payload from a configuration |
| [`Chirp`](inc/damp/estimation/excitation/chirp.hpp#L132) | block | Observers & estimators | Linear or logarithmic chirp runtime generator |
| [`ChirpConfig`](inc/damp/estimation/excitation/chirp.hpp#L40) | block | Observers & estimators | Configuration for a sine chirp excitation |
| [`ChirpMode`](inc/damp/estimation/excitation/chirp.hpp#L26) | enum | Observers & estimators | Chirp sweep law |
| [`ChirpResult`](inc/damp/estimation/excitation/chirp.hpp#L83) | block | Observers & estimators | Chirp design payload |
| [`cholesky`](inc/damp/matrix/decomposition.hpp#L82) | function | Linear algebra | Cholesky decomposition for positive-definite matrices |
| [`cholesky_solve`](inc/damp/matrix/solve.hpp#L159) | function | Linear algebra | Solve AX = B via Cholesky (A = LLᴴ) |
| [`clarke_park_transform`](inc/damp/math/transforms.hpp#L380) | function | Scalar math, complex & frames | Fused Clarke-Park transform (abc → dq) |
| [`clarke_transform`](inc/damp/math/transforms.hpp#L253) | function | Scalar math, complex & frames | Clarke transform (abc → αβ) |
| [`classical_dob`](inc/damp/estimation/dob.hpp#L353) | function | Observers & estimators | Synthesize a classical disturbance observer from a nominal plant and Q-filter |
| [`ClassicalDOB`](inc/damp/estimation/dob.hpp#L385) | block | Observers & estimators | Classical Pn^-1·Q disturbance observer runtime (bolt-on compensator) |
| [`ClassicalDobResult`](inc/damp/estimation/dob.hpp#L296) | block | Observers & estimators | Design result for the classical Pn^-1·Q disturbance observer |
| [`classify_range`](inc/damp/toolbox/conditioning.hpp#L396) | function | Embedded helpers (controls-adjacent utilities) | Classify x against the four band edges `[fault_lo (valid_lo, valid_hi) fault_hi]` (assumed ordered, non-decreasing) |
| [`closed_loop_poles`](inc/damp/design/stability.hpp#L335) | function | Design-time synthesis (not PWM-rate) | Compute closed-loop poles (eigenvalues) with state feedback |
| [`cohen_coon`](inc/damp/design/pid_design.hpp#L283) | function | Design-time synthesis (not PWM-rate) | Cohen-Coon tuning from first-order-plus-dead-time model |
| [`ColVec`](inc/damp/matrix/colvec.hpp#L27) | block | Linear algebra | Concrete Column vector specialization of Matrix<N, 1, T> |
| [`ColView`](inc/damp/matrix/views.hpp#L289) | block | Linear algebra | Non-owning column view of a matrix |
| [`comb_notch_window`](inc/damp/filters/moving_average.hpp#L34) | function | Filters & signal conditioning | Window length for a moving-average comb that notches f_notch and all its harmonics: N = round(fs / f_notch) |
| [`CommandProjectionResult`](inc/damp/controllers/action_governor.hpp#L185) | block | Design-time synthesis (not PWM-rate) | Result of an affine command projection |
| [`compensate_coulomb_viscous`](inc/damp/toolbox/conditioning.hpp#L178) | function | Embedded helpers (controls-adjacent utilities) | Coulomb + viscous friction compensation (static inverse, massless) |
| [`compensate_stribeck`](inc/damp/toolbox/conditioning.hpp#L205) | function | Embedded helpers (controls-adjacent utilities) | Stribeck-style friction compensation (static + Coulomb + viscous) |
| [`compensator_from_feature`](inc/damp/estimation/successive_compensator.hpp#L179) | function | Observers & estimators | Map an FRF feature to RBJ biquad coefficients |
| [`CompensatorKind`](inc/damp/estimation/successive_compensator.hpp#L38) | enum | Observers & estimators | Which biquad family to load for a feature |
| [`Complementary`](inc/damp/filters/complementary.hpp#L32) | block | Filters & signal conditioning | Scalar (1-D) complementary filter — fuse a fast rate with a slow absolute |
| [`ComplementaryFilter`](inc/damp/estimation/sensor_fusion.hpp#L64) | block | Observers & estimators | Simple complementary filter for orientation estimation |
| [`complex`](inc/damp/math/complex.hpp#L38) | block | Scalar math, complex & frames | Constexpr complex number class for compile-time computations |
| [`compute_eigenvalues`](inc/damp/matrix/eigen.hpp#L382) | function | Linear algebra | Compute the eigenvalues (and Schur vectors) of a real square matrix |
| [`constant_power_torque_limit`](inc/damp/motor/foc.hpp#L41) | function | Motor control pack (if present) | Constant-power torque ceiling Tₘₐₓ = P_rated / \|ω\| |
| [`ConstantInertiaFeedforward`](inc/damp/toolbox/actuator.hpp#L197) | block | Embedded helpers (controls-adjacent utilities) | Per-axis decoupled torque feedforward: `τ = J·a + b·v + τ_c·sign(v) + g` |
| [`continuous_lpf_exact_step`](inc/damp/filters/lowpass.hpp#L37) | function | Filters & signal conditioning | Exact step of continuous first-order LPF ẏ = −ω_c (y − u) |
| [`continuous_lqr`](inc/damp/controllers/lqr.hpp#L382) | function | Design-time synthesis (not PWM-rate) | Continuous-time Linear-Quadratic Regulator design (+1 more overload) |
| [`ContinuousPID`](inc/damp/controllers/pid.hpp#L685) | block | Runtime controllers | Continuous-gain PID with per-tick sample time (variable-rate secondary form) |
| [`controllability_gramian`](inc/damp/design/stability.hpp#L118) | function | Design-time synthesis (not PWM-rate) | Continuous/discrete controllability Gramian W_c |
| [`controllability_matrix`](inc/damp/design/stability.hpp#L52) | function | Design-time synthesis (not PWM-rate) | Compute the controllability matrix [B, AB, A²B, ..., A^(N-1)B] |
| [`Convention`](inc/damp/math/transforms.hpp#L69) | enum | Scalar math, complex & frames | Scaling convention for the Clarke/Park family |
| [`copysign`](inc/damp/math/math.hpp#L406) | function | Scalar math, complex & frames | Copy sign — magnitude of mag with the sign of sgn_src |
| [`cos`](inc/damp/matrix/functions.hpp#L738) | function | Linear algebra | Matrix cosine via scaling and double-angle reconstruction |
| [`cos`](inc/damp/math/math.hpp#L169) | function | Scalar math, complex & frames | Cosine |
| [`cosh`](inc/damp/matrix/functions.hpp#L809) | function | Linear algebra | Matrix hyperbolic cosine cosh(A) = (exp(A) + exp(−A))/2 |
| [`Counter`](inc/damp/toolbox/logic.hpp#L273) | block | Embedded helpers (controls-adjacent utilities) | Edge-counting up/down counter: increments on each rising edge of up, decrements on each rising edge of down. Returns the running count |
| [`cross_check_analog`](inc/damp/toolbox/conditioning.hpp#L573) | function | Embedded helpers (controls-adjacent utilities) | Cross-check two calibrated analog readings for redundant sensors |
| [`CTD`](inc/damp/toolbox/iec61131.hpp#L364) | block | Embedded helpers (controls-adjacent utilities) | CTD Counter (Count Down) |
| [`ctrb`](inc/damp/matlab.hpp#L265) | function | MATLAB®-style aliases (host) | MATLAB® short alias for controllability_matrix (+1 more overload) |
| [`CTU`](inc/damp/toolbox/iec61131.hpp#L324) | block | Embedded helpers (controls-adjacent utilities) | CTU Counter (Count Up) |
| [`CTUD`](inc/damp/toolbox/iec61131.hpp#L404) | block | Embedded helpers (controls-adjacent utilities) | CTUD Counter (Count Up Down) |
| [`damp`](inc/damp/analysis/poles.hpp#L90) | function | Frequency-domain analysis (host) | Compute natural frequency and damping for each pole |
| [`damping_ratio_from_overshoot_percent`](inc/damp/design/pid_design.hpp#L758) | function | Design-time synthesis (not PWM-rate) | Map percent overshoot target to equivalent damping ratio |
| [`dare`](inc/damp/design/riccati.hpp#L634) | function | Design-time synthesis (not PWM-rate) | Solve the Discrete Algebraic Riccati Equation (DARE) |
| [`db2mag`](inc/damp/math/math.hpp#L470) | function | Scalar math, complex & frames | Decibels to magnitude, 10^(db/20) |
| [`dcgain`](inc/damp/analysis/norms.hpp#L40) | function | Frequency-domain analysis (host) | Compute DC gain of a continuous-time system |
| [`DCM`](inc/damp/math/geometry.hpp#L63) | block | Scalar math, complex & frames | Direction cosine matrix — 3×3 rotation (SO(3) wrapper over Mat3) |
| [`deadband`](inc/damp/toolbox/conditioning.hpp#L47) | function | Embedded helpers (controls-adjacent utilities) | Dead zone over `[lower, upper]`, matching Simulink®'s Dead Zone block (+1 more overload) |
| [`Debounce`](inc/damp/toolbox/logic.hpp#L213) | block | Embedded helpers (controls-adjacent utilities) | Debounce: the output adopts in only after in differs from the current output continuously for stable_time. Rejects contact bounce and brief glitches. (Not an IEC block — the one everyone hand-rolls.) |
| [`decoupling_feedforward`](inc/damp/motor/foc.hpp#L69) | function | Motor control pack (if present) | Decoupling + back-EMF feedforward at R = 0 (reference currents) |
| [`default_tol`](inc/damp/matrix/matrix_traits.hpp#L86) | function | Linear algebra | Type-appropriate default tolerance for floating-point comparisons. float  ~7 decimal digits  → 1e-6 double ~15 decimal digits → 1e-12 |
| [`deg2rad`](inc/damp/math/math.hpp#L492) | function | Scalar math, complex & frames | Degrees to radians, deg·π/180 |
| [`Delay`](inc/damp/filters/delay.hpp#L27) | block | Filters & signal conditioning | Discrete-time delay buffer |
| [`det`](inc/damp/matrix/functions.hpp#L182) | function | Linear algebra | Matrix determinant det(A) |
| [`DFF`](inc/damp/toolbox/iec61131.hpp#L461) | block | Embedded helpers (controls-adjacent utilities) | D Flip-Flop (edge-triggered data latch) |
| [`diag`](inc/damp/matlab.hpp#L379) | function | MATLAB®-style aliases (host) | Returns a square diagonal matrix from the given array (+1 more overload) |
| [`Diagonal`](inc/damp/matrix/views.hpp#L49) | block | Linear algebra | Diagonal view of a square matrix |
| [`differential_drive`](inc/damp/toolbox/io.hpp#L340) | function | Embedded helpers (controls-adjacent utilities) | Differential (tank) drive mixer: throttle + turn → left/right |
| [`DirectQuadrature`](inc/damp/math/transforms.hpp#L84) | block | Scalar math, complex & frames | Direct-quadrature (rotor-frame) component pair |
| [`Discrete`](inc/damp/simulation/integrator.hpp#L83) | block | Simulation / SIL harness (host) | Discrete-time integrator (no integration, just one step) |
| [`discrete_lqe`](inc/damp/estimation/kalman.hpp#L162) | function | Observers & estimators | Discrete linear-quadratic estimator (LQE) — dual-of-LQR spelling of kalman (+1 more overload) |
| [`discrete_lqg`](inc/damp/controllers/lqg.hpp#L125) | function | Design-time synthesis (not PWM-rate) | Discrete Linear-Quadratic-Gaussian regulator design |
| [`discrete_lqgi`](inc/damp/controllers/lqgi.hpp#L133) | function | Design-time synthesis (not PWM-rate) | Discrete LQG with integral action (LQI + Kalman) for output tracking |
| [`discrete_lqi`](inc/damp/controllers/lqi.hpp#L93) | function | Design-time synthesis (not PWM-rate) | Discrete Linear-Quadratic-Integral (LQI) design for output tracking |
| [`discrete_lqr`](inc/damp/controllers/lqr.hpp#L121) | function | Design-time synthesis (not PWM-rate) | Discrete-time Linear-Quadratic Regulator design |
| [`discrete_lqr_from_continuous`](inc/damp/controllers/lqr.hpp#L259) | function | Design-time synthesis (not PWM-rate) | Design discrete LQR from continuous-time system via discretization (+1 more overload) |
| [`DiscretePIDResult`](inc/damp/controllers/pid.hpp#L46) | block | Design-time synthesis (not PWM-rate) | Fixed-rate discrete PID coefficients (canonical deploy form) |
| [`DiscretizationMethod`](inc/damp/systems/discretization.hpp#L27) | enum | LTI systems (SS / TF / ZPK / discretize) | Discretization methods for continuous-time state-space systems |
| [`discretize`](inc/damp/systems/discretization.hpp#L268) | function | LTI systems (SS / TF / ZPK / discretize) | Discretize a continuous-time state-space system |
| [`discretize_lqr_cost`](inc/damp/controllers/lqr.hpp#L187) | function | Design-time synthesis (not PWM-rate) | Discretize a continuous LQR cost integral over one sample (Van Loan) |
| [`DLATCH`](inc/damp/toolbox/iec61131.hpp#L492) | block | Embedded helpers (controls-adjacent utilities) | D Latch (level-sensitive / transparent latch) |
| [`dlqr`](inc/damp/controllers/lqr.hpp#L304) | function | Design-time synthesis (not PWM-rate) | Discrete-time LQR design (MATLAB®-style short name) |
| [`dlqr`](inc/damp/matlab.hpp#L625) | function | MATLAB®-style aliases (host) | Discrete-time Linear-Quadratic Regulator design |
| [`dlyap`](inc/damp/design/lyapunov.hpp#L130) | function | Design-time synthesis (not PWM-rate) | Solve the discrete-time Lyapunov (Stein) equation A X Aᵀ − X + Q = 0 |
| [`dlyap`](inc/damp/matlab.hpp#L1051) | function | MATLAB®-style aliases (host) | MATLAB® alias for the discrete Lyapunov solve AXAᵀ−X+Q=0 |
| [`dob`](inc/damp/estimation/dob.hpp#L129) | function | Observers & estimators | Validate and package DOB configuration into a runtime-ready design result |
| [`DOB`](inc/damp/estimation/dob.hpp#L154) | block | Observers & estimators | Lightweight SISO disturbance observer runtime |
| [`DOBConfig`](inc/damp/estimation/dob.hpp#L32) | block | Observers & estimators | Configuration for a first-order disturbance observer |
| [`DP45`](inc/damp/simulation/integrator.hpp#L757) | block | Simulation / SIL harness (host) | Dormand-Prince 5(4) adaptive integrator — the ode45 pair |
| [`dpwm_zero_sequence`](inc/damp/motor/modulation.hpp#L265) | function | Motor control pack (if present) | DPWM zero-sequence for a chosen discontinuous scheme |
| [`DqCommand`](inc/damp/motor/foc.hpp#L83) | block | Motor control pack (if present) | Clamped dq voltage command from FOController::current_controller |
| [`DrivePair`](inc/damp/toolbox/io.hpp#L315) | block | Embedded helpers (controls-adjacent utilities) | Left/right actuator pair (differential/arcade drive output) |
| [`DsogiPll`](inc/damp/filters/pll.hpp#L282) | block | Filters & signal conditioning | Dual-SOGI three-phase positive-sequence PLL (DSOGI-PLL) |
| [`duties_from_phase_voltages`](inc/damp/motor/modulation.hpp#L357) | function | Motor control pack (if present) | Map phase voltages + zero-sequence to clamped half-bridge duties |
| [`eig`](inc/damp/matlab.hpp#L480) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the eigenvalues of a square matrix |
| [`EigenResult`](inc/damp/matrix/eigen.hpp#L34) | block | Linear algebra | Eigenvalue computation result |
| [`EKFMeasFn`](inc/damp/estimation/ekf.hpp#L82) | concept | Observers & estimators | Concept for EKF measurement functions |
| [`EKFStateFn`](inc/damp/estimation/ekf.hpp#L54) | concept | Observers & estimators | Concept for EKF state functions |
| [`enu_from_ned`](inc/damp/estimation/ins_mechanization.hpp#L116) | function | Observers & estimators | Map a NED vector into ENU (axis permute) |
| [`ErrorStateJacobian`](inc/damp/estimation/eskf.hpp#L247) | block | Observers & estimators | ESKF prediction Jacobians F and G (nominal state updated outside) |
| [`ErrorStateKalmanFilter`](inc/damp/estimation/eskf.hpp#L294) | block | Observers & estimators | Runtime error-state KF: tracks δx and P, not the full nominal state |
| [`ErrorTolerance`](inc/damp/simulation/solver.hpp#L259) | block | Simulation / SIL harness (host) | Scalar-or-vector local-error tolerance for adaptive ODE control |
| [`esc`](inc/damp/controllers/esc.hpp#L153) | function | Design-time synthesis (not PWM-rate) | Synthesize an extremum-seeking controller |
| [`esc_lpf_alpha`](inc/damp/controllers/esc.hpp#L131) | function | Design-time synthesis (not PWM-rate) | First-order discrete low-pass coefficient for corner wc [rad/s] at Ts |
| [`esc_mppt`](inc/damp/controllers/esc.hpp#L200) | function | Design-time synthesis (not PWM-rate) | MPPT-flavored ESC: maximize a power measurement by perturbing the operating point (e.g. converter duty or reference voltage) |
| [`ESCConfig`](inc/damp/controllers/esc.hpp#L59) | block | Runtime controllers | Extremum-seeking controller configuration (discrete realization) |
| [`ESCResult`](inc/damp/controllers/esc.hpp#L103) | block | Design-time synthesis (not PWM-rate) | Design result for the extremum-seeking controller |
| [`eskf_imu`](inc/damp/estimation/eskf.hpp#L179) | function | Observers & estimators | Design Q, R, P₀ for a 6-state IMU attitude ESKF (gyro + accel) |
| [`eskf_marg`](inc/damp/estimation/eskf.hpp#L218) | function | Observers & estimators | Design Q, R, P₀ for a 6-state MARG attitude ESKF (gyro + accel + mag) |
| [`ESKFMeasFn`](inc/damp/estimation/eskf.hpp#L270) | concept | Observers & estimators | Callable that returns a measurement linearization for one update |
| [`ESKFOrientationFilter`](inc/damp/estimation/sensor_fusion.hpp#L399) | block | Observers & estimators | Turnkey 6-state attitude ESKF: owns q, b_g, and the error filter |
| [`ESKFPredictFn`](inc/damp/estimation/eskf.hpp#L259) | concept | Observers & estimators | Callable that returns ErrorStateJacobian for one predict step |
| [`ESKFResult`](inc/damp/estimation/eskf.hpp#L50) | block | Observers & estimators | Design payload for an ESKF: Q, R, P₀, and success |
| [`estim`](inc/damp/matlab.hpp#L497) | function | MATLAB®-style aliases (host) | Form state estimator from system and estimator gain |
| [`Euler`](inc/damp/math/geometry.hpp#L205) | block | Scalar math, complex & frames | Intrinsic Euler angles for sequence Order (default ZYX aerospace YPR) |
| [`EulerOrder`](inc/damp/math/geometry.hpp#L39) | enum | Scalar math, complex & frames | Intrinsic Euler angle sequence (axis order of successive principal rotations) |
| [`eval_frf`](inc/damp/systems/state_space.hpp#L169) | function | LTI systems (SS / TF / ZPK / discretize) | Evaluate frequency response of a state-space system |
| [`EventSimConfig`](inc/damp/simulation/hybrid.hpp#L67) | block | Simulation / SIL harness (host) | Configuration for event-driven hybrid / multi-rate style runs |
| [`Exact`](inc/damp/simulation/integrator.hpp#L117) | block | Simulation / SIL harness (host) | Exact integrator for LTI systems |
| [`exact_lti_step`](inc/damp/simulation/hybrid.hpp#L315) | function | Simulation / SIL harness (host) | Exact open-loop step of a single continuous LTI system |
| [`exp`](inc/damp/math/math.hpp#L232) | function | Scalar math, complex & frames | Exponential function |
| [`ExperimentSafety`](inc/damp/estimation/experiment_safety.hpp#L62) | block | Observers & estimators | Experiment lifecycle + command protection for on-target commissioning |
| [`ExperimentSafetyConfig`](inc/damp/estimation/experiment_safety.hpp#L39) | block | Observers & estimators | Configuration for ExperimentSafety |
| [`expm`](inc/damp/matrix/functions.hpp#L323) | function | Linear algebra | Matrix exponential via scaling and squaring with Padé approximant of degree 13 |
| [`expo`](inc/damp/toolbox/conditioning.hpp#L252) | function | Embedded helpers (controls-adjacent utilities) | Exponential response curve `y = (1−k)·x + k·x³` (RC "expo") |
| [`ExtendedKalmanFilter`](inc/damp/estimation/ekf.hpp#L112) | block | Observers & estimators | Extended Kalman Filter for nonlinear discrete-time systems |
| [`extract_channel`](inc/damp/simulation/plot_plotly.hpp#L223) | function | Simulation / SIL harness (host) | Extract the i-th element from each vector entry into a std::vector<double> |
| [`extract_frf_features`](inc/damp/estimation/frequency_response.hpp#L740) | function | Observers & estimators | Extract peaks and valleys in one pass over the FRF table |
| [`extract_modes`](inc/damp/estimation/frequency_response.hpp#L715) | function | Observers & estimators | Extract resonant peaks (local \|G\| maxima) from an FRF table |
| [`extract_valleys`](inc/damp/estimation/frequency_response.hpp#L729) | function | Observers & estimators | Extract anti-resonance valleys (local \|G\| minima) from an FRF table |
| [`Extrapolation`](inc/damp/toolbox/lookup.hpp#L36) | enum | Embedded helpers (controls-adjacent utilities) | Out-of-range behaviour for a Lut1D / Lut2D / Lut3D query beyond its breakpoints |
| [`ExtremumSeekingController`](inc/damp/controllers/esc.hpp#L238) | block | Runtime controllers | Extremum-seeking controller runtime (model-free online optimizer) |
| [`ExtremumType`](inc/damp/controllers/esc.hpp#L47) | enum | Runtime controllers | Whether ESC climbs to a maximum or descends to a minimum of the objective |
| [`eye`](inc/damp/matlab.hpp#L420) | function | MATLAB®-style aliases (host) | Create an identity matrix of size n x n |
| [`F_TRIG`](inc/damp/toolbox/iec61131.hpp#L153) | block | Embedded helpers (controls-adjacent utilities) | F_TRIG (Falling Edge Trigger) |
| [`factorial`](inc/damp/trajectory/trajectory_types.hpp#L206) | function | Trajectory value types | k! (exact for modest k in float/double; used up through nonic BVP / jets) |
| [`falling_factorial`](inc/damp/trajectory/trajectory_types.hpp#L217) | function | Trajectory value types | Falling factorial i·(i−1)···(i−k+1) = i! / (i−k)! — the k-th derivative coefficient of tⁱ. Zero when i < k |
| [`FallingEdge`](inc/damp/toolbox/logic.hpp#L44) | block | Embedded helpers (controls-adjacent utilities) | Falling-edge detector: true on the tick x goes true → false |
| [`feedback`](inc/damp/systems/state_space.hpp#L314) | function | LTI systems (SS / TF / ZPK / discretize) | Negative feedback: y = sys1(u − sys2(y)) (+1 more overload) |
| [`feedback_read`](inc/damp/systems/state_space.hpp#L683) | function | LTI systems (SS / TF / ZPK / discretize) | Unity (or weighted) negative feedback around G_loop, read G_out |
| [`finish_extremum`](inc/damp/estimation/frequency_response.hpp#L561) | function | Observers & estimators | Half-power (peaks) or double-power (valleys) bandwidth → ζ |
| [`finite_non_negative`](inc/damp/math/math.hpp#L448) | function | Scalar math, complex & frames | True if x is finite and non-negative (x ≥ 0) |
| [`finite_positive`](inc/damp/math/math.hpp#L438) | function | Scalar math, complex & frames | True if x is finite and strictly greater than zero |
| [`fir1`](inc/damp/filters/fir.hpp#L218) | function | Filters & signal conditioning | Window-method FIR design (normalized frequency, single cutoff) (+1 more overload) |
| [`fir1_hz`](inc/damp/filters/fir.hpp#L372) | function | Filters & signal conditioning | Window-method FIR design from frequencies in Hz (+1 more overload) |
| [`fir_window`](inc/damp/filters/fir.hpp#L408) | function | Filters & signal conditioning | Alias for fir1 (descriptive name) (+1 more overload) |
| [`FirDesignResult`](inc/damp/filters/fir.hpp#L67) | block | Filters & signal conditioning | Window-method FIR design result |
| [`fires_at`](inc/damp/simulation/multirate.hpp#L73) | function | Simulation / SIL harness (host) | True when base tick tick is a firing instant for divisor |
| [`FirFilter`](inc/damp/filters/fir.hpp#L450) | block | Filters & signal conditioning | Direct-form FIR filter runtime (tapped delay line + dot product) |
| [`FirstOrderCoeffs`](inc/damp/filters/iir_design.hpp#L30) | block | Filters & signal conditioning | DSP coefficients for first-order IIR filter |
| [`FirstOrderPlantEstimator`](inc/damp/estimation/parameter_estimation.hpp#L126) | block | Observers & estimators | Online estimator for a first-order plant's gain and time constant |
| [`FirstOrderPlantEstimatorConfig`](inc/damp/estimation/parameter_estimation.hpp#L85) | block | Observers & estimators | Configuration for the first-order grey-box estimator |
| [`FirType`](inc/damp/filters/fir.hpp#L38) | enum | Filters & signal conditioning | Ideal frequency-selective response for fir1 / fir_window |
| [`FirWindow`](inc/damp/filters/fir.hpp#L50) | enum | Filters & signal conditioning | Window applied to the ideal truncated impulse response |
| [`fit_second_order_step`](inc/damp/estimation/successive_compensator.hpp#L85) | function | Observers & estimators | Fit one underdamped second-order mode from a step / ring-down capture |
| [`fixed_solve`](inc/damp/simulation/solver.hpp#L243) | function | Simulation / SIL harness (host) | Fixed-step solve in one call |
| [`FixedStepSolver`](inc/damp/simulation/solver.hpp#L118) | block | Simulation / SIL harness (host) | Fixed-step ODE solver (+1 more overload) |
| [`floor`](inc/damp/math/math.hpp#L317) | function | Scalar math, complex & frames | Floor — largest integer ≤ x |
| [`flux_from_Kv`](inc/damp/motor/spm.hpp#L55) | function | Motor control pack (if present) | λ from Kᵥ via torque_constant_from_Kv |
| [`flux_from_torque_constant`](inc/damp/motor/spm.hpp#L37) | function | Motor control pack (if present) | λ from datasheet Kₜ (amplitude / peak-per-phase convention) |
| [`fmod`](inc/damp/math/math.hpp#L361) | function | Scalar math, complex & frames | Floating-point remainder, x − y·trunc(x/y) (sign of x), matching std::fmod's truncated-quotient convention |
| [`FOController`](inc/damp/motor/foc.hpp#L98) | block | Motor control pack (if present) | dq current regulator: PI + decoupling FF → clamped Vdq |
| [`forward_substitute`](inc/damp/matrix/solve.hpp#L37) | function | Linear algebra | Forward substitution for Lx = b |
| [`ForwardEuler`](inc/damp/simulation/integrator.hpp#L169) | block | Simulation / SIL harness (host) | Forward Euler integrator |
| [`FrequencyPoint`](inc/damp/analysis/frequency.hpp#L38) | block | Frequency-domain analysis (host) | Single-point frequency response result |
| [`FrequencyResponseEstimator`](inc/damp/estimation/frequency_response.hpp#L106) | block | Observers & estimators | Lock-in FRF estimator over a fixed table of NFreq bins |
| [`frf_pid_autotune`](inc/damp/estimation/commissioning_frontends.hpp#L152) | function | Observers & estimators | PID gains from an open-loop FRF table + margin / bandwidth targets |
| [`frf_s_or_z`](inc/damp/analysis/frequency.hpp#L974) | function | Frequency-domain analysis (host) | FRF evaluation point: jω (continuous) or e^{jωTs} (discrete) |
| [`FrfFeatureKind`](inc/damp/estimation/frequency_response.hpp#L373) | enum | Observers & estimators | Peak (resonance) or valley (anti-resonance) on an FRF magnitude curve |
| [`FrfFeatureResult`](inc/damp/estimation/frequency_response.hpp#L453) | block | Observers & estimators | Combined peak + valley reduction of an FRF table |
| [`FrfMargins`](inc/damp/estimation/frequency_response.hpp#L485) | block | Observers & estimators | Open-loop margins read from a discrete FRF table (on-target) |
| [`FrfPidAutotuneConfig`](inc/damp/estimation/commissioning_frontends.hpp#L109) | block | Observers & estimators | Configuration for FRF-based PID gain selection |
| [`FrfPidAutotuneResult`](inc/damp/estimation/commissioning_frontends.hpp#L125) | block | Observers & estimators | FRF PID autotune hand-off (gains + success) |
| [`FrfPoint`](inc/damp/estimation/frequency_response.hpp#L60) | block | Observers & estimators | One measured frequency-response point (on-target FRF table entry) |
| [`frobenius_norm`](inc/damp/matrix/functions.hpp#L96) | function | Linear algebra | Frobenius norm ‖A‖F = √(Σᵢⱼ \|aᵢⱼ\|²) |
| [`full_qr`](inc/damp/matrix/decomposition.hpp#L291) | function | Linear algebra | Full QR factorization via Householder reflections (real or complex T) |
| [`FullQR`](inc/damp/matrix/decomposition.hpp#L273) | block | Linear algebra | Result of a full (complete) QR factorization |
| [`gain_margin_unwrapped`](inc/damp/analysis/frequency.hpp#L202) | function | Frequency-domain analysis (host) | Find gain margin using unwrapped phase trajectory |
| [`geomspace`](inc/damp/analysis/linspace.hpp#L88) | function | Frequency-domain analysis (host) | Geometric sequence from start to end (endpoint values, not exponents) |
| [`Goertzel`](inc/damp/filters/spectral.hpp#L46) | block | Filters & signal conditioning | Generalized Goertzel single-bin DFT — amplitude/phase at one frequency |
| [`gram`](inc/damp/matlab.hpp#L1098) | function | MATLAB®-style aliases (host) | MATLAB® alias for the controllability/observability Gramian of a system |
| [`gravity_nav`](inc/damp/estimation/ins_mechanization.hpp#L93) | function | Observers & estimators | Gravity vector in the local-level frame [m/s²] |
| [`h_pattern`](inc/damp/toolbox/io.hpp#L377) | function | Embedded helpers (controls-adjacent utilities) | H-pattern (dual-lever / tank) drive: one lever per track, no mixing |
| [`hankel_singular_values`](inc/damp/design/model_reduction.hpp#L596) | function | Design-time synthesis (not PWM-rate) | Descriptive alias for hankelsv |
| [`hankelsv`](inc/damp/design/model_reduction.hpp#L439) | function | Design-time synthesis (not PWM-rate) | Hankel singular values of a stable state-space system |
| [`harmonic_suppressor`](inc/damp/controllers/harmonic_suppression.hpp#L113) | function | Design-time synthesis (not PWM-rate) | Synthesize a multi-resonant harmonic suppressor |
| [`HarmonicAnalyzer`](inc/damp/filters/spectral.hpp#L173) | block | Filters & signal conditioning | Harmonic analyzer — a Goertzel bank over a fundamental and K−1 harmonics |
| [`HarmonicSuppressor`](inc/damp/controllers/harmonic_suppression.hpp#L154) | block | Runtime controllers | Multi-resonant harmonic suppressor — a parallel bank of PR resonators |
| [`HarmonicSuppressorResult`](inc/damp/controllers/harmonic_suppression.hpp#L48) | block | Design-time synthesis (not PWM-rate) | Design result for a multi-resonant harmonic suppressor |
| [`heading_from_baseline_nav`](inc/damp/estimation/ins_eskf.hpp#L536) | function | Observers & estimators | Heading [rad] from a nav-frame baseline vector (dual-antenna difference) |
| [`Heun`](inc/damp/simulation/integrator.hpp#L840) | block | Simulation / SIL harness (host) | Heun's method (Improved Euler, RK2) integrator |
| [`HighPass`](inc/damp/filters/highpass.hpp#L40) | block | Filters & signal conditioning | First-order high-pass (washout) filter runtime |
| [`highpass_1st`](inc/damp/filters/iir_design.hpp#L169) | function | Filters & signal conditioning | First-order high-pass filter design (Tustin / bilinear) (+1 more overload) |
| [`highpass_2nd`](inc/damp/filters/iir_design.hpp#L360) | function | Filters & signal conditioning | Second-order continuous high-pass at default Butterworth Q (1/√2) (+1 more overload) |
| [`highpass_2nd_continuous`](inc/damp/filters/iir_design.hpp#L338) | function | Filters & signal conditioning | Second-order high-pass filter design (continuous-time) |
| [`highshelf`](inc/damp/filters/iir_design.hpp#L935) | function | Filters & signal conditioning | High-shelf EQ filter: boost or cut everything above fc |
| [`hinfnorm`](inc/damp/matlab.hpp#L1120) | function | MATLAB®-style aliases (host) | MATLAB® alias for the H∞ system norm norm(sys,Inf) / hinfnorm(sys) |
| [`holonomic_drive`](inc/damp/toolbox/io.hpp#L495) | function | Embedded helpers (controls-adjacent utilities) | Four-wheel holonomic drive mixer (mecanum / omni) |
| [`HolonomicOutput`](inc/damp/toolbox/io.hpp#L463) | block | Embedded helpers (controls-adjacent utilities) | Four-wheel actuator set (mecanum/holonomic drive output), X-config |
| [`HybridEventAction`](inc/damp/simulation/hybrid.hpp#L102) | block | Simulation / SIL harness (host) | Action applied when an event fires (after Exact advance to the event time) |
| [`HybridSimulationResult`](inc/damp/simulation/hybrid.hpp#L80) | block | Simulation / SIL harness (host) | Result of a piecewise-LTI hybrid simulation |
| [`hypot`](inc/damp/math/math.hpp#L76) | function | Scalar math, complex & frames | Euclidean distance hypot(x, y) = √(x² + y²), without overflow |
| [`Hysteresis`](inc/damp/toolbox/conditioning.hpp#L332) | block | Embedded helpers (controls-adjacent utilities) | Hysteresis comparator (Schmitt trigger): bool output with separate on/off thresholds to reject chatter |
| [`impedance`](inc/damp/analysis/frequency.hpp#L999) | function | Frequency-domain analysis (host) | Compute impedance frequency response from a SISO admittance system |
| [`impedance_direct`](inc/damp/analysis/frequency.hpp#L1039) | function | Frequency-domain analysis (host) | Compute impedance frequency response from a SISO impedance transfer function |
| [`ImpedanceResult`](inc/damp/analysis/frequency.hpp#L872) | block | Frequency-domain analysis (host) | Result of impedance frequency response evaluation |
| [`impulse`](inc/damp/analysis/time_response.hpp#L185) | function | Frequency-domain analysis (host) | Impulse response of a (MIMO) state-space system (+1 more overload) |
| [`impulse`](inc/damp/estimation/excitation/impulse.hpp#L70) | function | Observers & estimators | Build an impulse design payload |
| [`Impulse`](inc/damp/estimation/excitation/impulse.hpp#L86) | block | Observers & estimators | Rectangular impulse (finite-width Dirac stand-in) |
| [`ImpulseConfig`](inc/damp/estimation/excitation/impulse.hpp#L32) | block | Observers & estimators | Configuration for a one-shot rectangular impulse (Dirac stand-in) |
| [`impulseplot`](inc/damp/simulation/plot_plotly.hpp#L733) | function | Simulation / SIL harness (host) | Plot an impulse response, one trace per input/output pair |
| [`ImpulseResult`](inc/damp/estimation/excitation/impulse.hpp#L52) | block | Observers & estimators | Impulse design payload |
| [`ImuSample`](inc/damp/estimation/ins_mechanization.hpp#L125) | block | Observers & estimators | IMU sample in the body frame |
| [`InductorCurrentPIResult`](inc/damp/design/pid_design.hpp#L1137) | block | Design-time synthesis (not PWM-rate) | Fixed-rate inductor current-loop PI (topology-agnostic) |
| [`infinity_norm`](inc/damp/matrix/functions.hpp#L35) | function | Linear algebra | Infinity norm ‖A‖∞: maximum absolute row sum |
| [`initial`](inc/damp/analysis/time_response.hpp#L219) | function | Frequency-domain analysis (host) | Initial-condition (free) response of a (MIMO) state-space system |
| [`ins_aid_position`](inc/damp/estimation/ins_eskf.hpp#L503) | function | Observers & estimators | Predicted antenna / marker position with body lever-arm |
| [`ins_apply_correction`](inc/damp/estimation/ins_eskf.hpp#L489) | function | Observers & estimators | Inject filter correction into x and reset the error state |
| [`ins_error_jacobian`](inc/damp/estimation/ins_eskf.hpp#L307) | function | Observers & estimators | Discrete error-state transition and noise input for one IMU step |
| [`ins_eskf_design`](inc/damp/estimation/ins_eskf.hpp#L219) | function | Observers & estimators | Build Q / R / P0 for a 15-state INS ESKF from IMU noise densities |
| [`ins_heading`](inc/damp/estimation/ins_eskf.hpp#L521) | function | Observers & estimators | Horizontal heading [rad] of a body axis expressed in the nav frame |
| [`ins_inject`](inc/damp/estimation/ins_eskf.hpp#L433) | function | Observers & estimators | Inject δx into the nominal INS state (right-multiplicative q) |
| [`ins_predict`](inc/damp/estimation/ins_eskf.hpp#L466) | function | Observers & estimators | Mechanize nominal state and propagate the error covariance |
| [`ins_propagate_covariance`](inc/damp/estimation/ins_eskf.hpp#L403) | function | Observers & estimators | P ← F P Fᵀ + Q using sparse F (exact for ins_error_jacobian structure) |
| [`ins_update_heading`](inc/damp/estimation/ins_eskf.hpp#L638) | function | Observers & estimators | Sparse dual-antenna / yaw heading update (scalar) |
| [`ins_update_orientation`](inc/damp/estimation/ins_eskf.hpp#L676) | function | Observers & estimators | Sparse orientation update from a measured body→nav quaternion |
| [`ins_update_pose`](inc/damp/estimation/ins_eskf.hpp#L701) | function | Observers & estimators | Generic pose alias: position then orientation (two sparse updates) |
| [`ins_update_pose_heading`](inc/damp/estimation/ins_eskf.hpp#L720) | function | Observers & estimators | Pose alias: position + dual-antenna heading (two sparse updates) |
| [`ins_update_position`](inc/damp/estimation/ins_eskf.hpp#L551) | function | Observers & estimators | Sparse position update: y = p + Rℓ (+1 more overload) |
| [`ins_update_velocity`](inc/damp/estimation/ins_eskf.hpp#L597) | function | Observers & estimators | Sparse velocity update: y = v at the IMU origin (no lever-arm) |
| [`ins_update_zupt`](inc/damp/estimation/ins_eskf.hpp#L616) | function | Observers & estimators | Zero-velocity update (ZUPT): velocity measurement of zero |
| [`InsEskfResult`](inc/damp/estimation/ins_eskf.hpp#L180) | block | Observers & estimators | Design result for the 15-state INS error-state filter |
| [`InsNavigator`](inc/damp/estimation/ins_eskf.hpp#L754) | block | Observers & estimators | Runtime strapdown navigator: nominal InsState + 15-state ESKF |
| [`InsState`](inc/damp/estimation/ins_mechanization.hpp#L139) | block | Observers & estimators | Nominal strapdown navigation state (external to the ESKF error state) |
| [`instantaneous_power`](inc/damp/math/transforms.hpp#L472) | function | Scalar math, complex & frames | Instantaneous active and reactive power from dq quantities (+1 more overload) |
| [`InstantaneousPower`](inc/damp/math/transforms.hpp#L433) | block | Scalar math, complex & frames | Instantaneous active and reactive power |
| [`IntegrationResult`](inc/damp/simulation/integrator.hpp#L31) | block | Simulation / SIL harness (host) | Result of an integration step |
| [`InteriorPointSolver`](inc/damp/design/qp.hpp#L1163) | block | Design-time synthesis (not PWM-rate) | Interior-point QP solver policy for damp::MPC |
| [`Interpolation`](inc/damp/toolbox/lookup.hpp#L42) | enum | Embedded helpers (controls-adjacent utilities) | In-range interpolant for Lut1D (2-D stays bilinear; use SplineSurface for C¹ grids) |
| [`inverse`](inc/damp/matrix/solve.hpp#L248) | function | Linear algebra | Matrix inverse A⁻¹ via mat::solve(A, I) |
| [`inverse_clarke_transform`](inc/damp/math/transforms.hpp#L284) | function | Scalar math, complex & frames | Inverse Clarke transform (αβ → abc) |
| [`inverse_deadband`](inc/damp/toolbox/conditioning.hpp#L79) | function | Embedded helpers (controls-adjacent utilities) | Inverse dead zone: add an offset to overcome a physical dead zone (valve overlap, static friction, motor stiction), with independent negative/positive offsets (+1 more overload) |
| [`inverse_deadband_sine`](inc/damp/toolbox/conditioning.hpp#L152) | function | Embedded helpers (controls-adjacent utilities) | Sine-shaped inverse dead zone (full boost once \|x\|≥w) |
| [`inverse_deadband_soft`](inc/damp/toolbox/conditioning.hpp#L128) | function | Embedded helpers (controls-adjacent utilities) | Smooth inverse dead zone via soft-sign (tanh-like Coulomb boost) |
| [`inverse_lerp`](inc/damp/toolbox/scaling.hpp#L48) | function | Embedded helpers (controls-adjacent utilities) | Inverse of lerp: the fraction t such that `lerp(a, b, t) == x` |
| [`inverse_park_clarke_transform`](inc/damp/math/transforms.hpp#L405) | function | Scalar math, complex & frames | Fused inverse Park-Clarke transform (dq → abc) |
| [`inverse_park_transform`](inc/damp/math/transforms.hpp#L358) | function | Scalar math, complex & frames | Inverse Park transform (dq → αβ) |
| [`inverse_symmetrical_components`](inc/damp/math/transforms.hpp#L551) | function | Scalar math, complex & frames | Inverse symmetrical-component transform (012 → abc) |
| [`iq_from_torque`](inc/damp/motor/spm.hpp#L71) | function | Motor control pack (if present) | i_q = T_e / Kₜ for id = 0 |
| [`is_closed_loop_stable_discrete`](inc/damp/design/stability.hpp#L245) | function | Design-time synthesis (not PWM-rate) | Check closed-loop stability for discrete system with state feedback |
| [`is_controllable`](inc/damp/design/stability.hpp#L176) | function | Design-time synthesis (not PWM-rate) | Check if a system is controllable |
| [`is_fault`](inc/damp/toolbox/conditioning.hpp#L382) | function | Embedded helpers (controls-adjacent utilities) | True for a wire fault (FaultLow/FaultHigh) — i.e. not a real reading at all |
| [`is_matrix_element`](inc/damp/matrix/matrix_traits.hpp#L41) | block | Linear algebra | True if T may be a Matrix element type |
| [`is_observable`](inc/damp/design/stability.hpp#L194) | function | Design-time synthesis (not PWM-rate) | Check if a system is observable |
| [`is_stabilizable`](inc/damp/design/riccati.hpp#L44) | function | Design-time synthesis (not PWM-rate) | Check if (A, B) is a stabilizable pair |
| [`is_stable_continuous`](inc/damp/analysis/poles.hpp#L58) | function | Frequency-domain analysis (host) | Check continuous-time stability |
| [`is_stable_discrete`](inc/damp/design/stability.hpp#L215) | function | Design-time synthesis (not PWM-rate) | Check if a discrete-time system matrix A is stable |
| [`is_usable`](inc/damp/toolbox/conditioning.hpp#L539) | function | Embedded helpers (controls-adjacent utilities) | True when a channel reading may still be used (in-span or saturated) |
| [`is_valid`](inc/damp/toolbox/conditioning.hpp#L377) | function | Embedded helpers (controls-adjacent utilities) | True only for the in-span status |
| [`isfinite`](inc/damp/math/math.hpp#L425) | function | Scalar math, complex & frames | Finiteness test — false for NaN and ±∞ |
| [`iso_c_drive`](inc/damp/toolbox/io.hpp#L412) | function | Embedded helpers (controls-adjacent utilities) | ISO-C (single-stick, coordinated / curvature) drive |
| [`iso_s_drive`](inc/damp/toolbox/io.hpp#L393) | function | Embedded helpers (controls-adjacent utilities) | ISO-S (single-stick, speed-summed) drive — arcade travel stick |
| [`isstable`](inc/damp/matlab.hpp#L945) | function | MATLAB®-style aliases (host) | Continuous-time stability predicate on a state matrix (+2 more overloads) |
| [`Jet`](inc/damp/trajectory/trajectory_types.hpp#L95) | block | Trajectory value types | Truncated jet of a scalar motion sample: d[k] = s⁽ᵏ⁾ |
| [`JordanBlock`](inc/damp/design/pole_placement.hpp#L831) | block | Design-time synthesis (not PWM-rate) | One Jordan mini-block of a desired closed-loop spectrum |
| [`JordanObjective`](inc/damp/design/pole_placement.hpp#L1165) | enum | Design-time synthesis (not PWM-rate) | Robustness objective for place_jordan_optimal (the paper's two methods) |
| [`kalman`](inc/damp/estimation/kalman.hpp#L103) | function | Observers & estimators | Steady-state Kalman filter design |
| [`KalmanFilter`](inc/damp/estimation/kalman.hpp#L236) | block | Observers & estimators | Runtime Kalman filter for embedded systems |
| [`KalmanResult`](inc/damp/estimation/kalman.hpp#L45) | block | Observers & estimators | Steady-state Kalman filter design result |
| [`lag`](inc/damp/controllers/lead_lag.hpp#L181) | function | Design-time synthesis (not PWM-rate) | Design a lag compensator from desired low-frequency gain boost |
| [`lambda_tuning`](inc/damp/design/pid_design.hpp#L401) | function | Design-time synthesis (not PWM-rate) | Lambda tuning for FOPDT model |
| [`Latch`](inc/damp/toolbox/logic.hpp#L62) | block | Embedded helpers (controls-adjacent utilities) | Set/reset latch. @tparam SetDominant which input wins when both are asserted (default: set-dominant, e.g. a trip overriding a clear) |
| [`lead`](inc/damp/controllers/lead_lag.hpp#L132) | function | Design-time synthesis (not PWM-rate) | Design a lead compensator from desired phase boost at a target frequency |
| [`lead_lag`](inc/damp/controllers/lead_lag.hpp#L234) | function | Design-time synthesis (not PWM-rate) | Design a lead-lag compensator (cascade of lead + lag sections) |
| [`lead_lag_direct`](inc/damp/controllers/lead_lag.hpp#L264) | function | Design-time synthesis (not PWM-rate) | Direct lead-lag specification from zero/pole locations |
| [`LeadLagController`](inc/damp/controllers/lead_lag.hpp#L281) | block | Runtime controllers | Discrete lead-lag compensator |
| [`LeadLagResult`](inc/damp/controllers/lead_lag.hpp#L59) | block | Design-time synthesis (not PWM-rate) | Lead-lag compensator design result |
| [`LeadLagSeriesResult`](inc/damp/controllers/lead_lag.hpp#L204) | block | Design-time synthesis (not PWM-rate) | Cascaded lead+lag design result (2nd-order StateSpace + success) |
| [`lerp`](inc/damp/toolbox/scaling.hpp#L35) | function | Embedded helpers (controls-adjacent utilities) | Linear interpolation between a and b by fraction t |
| [`LeverrierResult`](inc/damp/systems/zpk.hpp#L500) | block | LTI systems (SS / TF / ZPK / discretize) | Faddeev–LeVerrier characteristic polynomial and adjoint coefficient matrices |
| [`LIMIT`](inc/damp/toolbox/iec61131.hpp#L598) | function | Embedded helpers (controls-adjacent utilities) | LIMIT (IEC 61131-3 selection function): clamp in to [mn, mx] |
| [`linear_modulation_voltage`](inc/damp/motor/modulation.hpp#L127) | function | Motor control pack (if present) | Peak phase voltage at the linear SVPWM hexagon limit |
| [`linear_screw`](inc/damp/toolbox/actuator.hpp#L166) | function | Embedded helpers (controls-adjacent utilities) | Build a ServoAxis for a linear axis driven by a leadscrew/belt |
| [`LinearizationResult`](inc/damp/design/linearization.hpp#L45) | block | Design-time synthesis (not PWM-rate) | Result of nonlinear operating-point linearization |
| [`linearize`](inc/damp/design/linearization.hpp#L139) | function | Design-time synthesis (not PWM-rate) | Linearize nonlinear dynamics and output maps about an operating point (+1 more overload) |
| [`linmod`](inc/damp/matlab.hpp#L336) | function | MATLAB®-style aliases (host) | MATLAB®-style nonlinear linearization about an operating point |
| [`locate_zero_crossing_exact`](inc/damp/simulation/hybrid.hpp#L335) | function | Simulation / SIL harness (host) | Locate a state zero-crossing of g(x) on an Exact LTI segment |
| [`log`](inc/damp/matrix/functions.hpp#L511) | function | Linear algebra | Principal matrix logarithm via inverse scaling and squaring |
| [`log`](inc/damp/math/math.hpp#L249) | function | Scalar math, complex & frames | Natural logarithm |
| [`log10`](inc/damp/math/math.hpp#L380) | function | Scalar math, complex & frames | Base-10 logarithm, log10(x) = ln(x) / ln(10) |
| [`logm`](inc/damp/matrix/functions.hpp#L550) | function | Linear algebra | @brief MATLAB®-style alias for log |
| [`LogPolicy`](inc/damp/simulation/hybrid.hpp#L58) | enum | Simulation / SIL harness (host) | How often to record samples along a hybrid trajectory |
| [`loop_metrics`](inc/damp/analysis/frequency.hpp#L480) | function | Frequency-domain analysis (host) | One-call loop analysis: compute L/S/T response and return compact metrics (+1 more overload) |
| [`loop_response`](inc/damp/analysis/frequency.hpp#L558) | function | Frequency-domain analysis (host) | Compute open-loop L, sensitivity S, complementary sensitivity T, and Nyquist data (+1 more overload) |
| [`LoopResponseResult`](inc/damp/analysis/frequency.hpp#L415) | block | Frequency-domain analysis (host) | Open-loop and closed-loop frequency response package |
| [`LoopSummary`](inc/damp/analysis/frequency.hpp#L448) | block | Frequency-domain analysis (host) | Compact loop summary metrics for quick stability/robustness checks |
| [`LowerTriangle`](inc/damp/matrix/views.hpp#L157) | block | Linear algebra | Lower triangular view of a square matrix |
| [`LowPass`](inc/damp/filters/lowpass.hpp#L49) | block | Filters & signal conditioning | Nth-order low-pass filter |
| [`lowpass_1st`](inc/damp/filters/iir_design.hpp#L136) | function | Filters & signal conditioning | First-order low-pass filter design (+1 more overload) |
| [`lowpass_2nd`](inc/damp/filters/iir_design.hpp#L221) | function | Filters & signal conditioning | Second-order low-pass filter design (+1 more overload) |
| [`lowpass_2nd_continuous`](inc/damp/filters/iir_design.hpp#L272) | function | Filters & signal conditioning | Second-order low-pass filter design (continuous-time) |
| [`lowshelf`](inc/damp/filters/iir_design.hpp#L903) | function | Filters & signal conditioning | Low-shelf EQ filter: boost or cut everything below fc |
| [`lqg`](inc/damp/matlab.hpp#L744) | function | MATLAB®-style aliases (host) | Linear-Quadratic-Gaussian regulator design |
| [`LQG`](inc/damp/controllers/lqg.hpp#L181) | block | Runtime controllers | Linear-Quadratic-Gaussian (LQG) controller |
| [`lqg_bundle`](inc/damp/design/synthesis.hpp#L329) | function | Design-time synthesis (not PWM-rate) | Synthesize the full LQG artifact bundle in one call |
| [`lqg_from_parts`](inc/damp/controllers/lqg.hpp#L154) | function | Design-time synthesis (not PWM-rate) | Assemble an LQG design from separately computed Kalman and LQR results |
| [`lqg_pr_bundle`](inc/damp/design/synthesis.hpp#L356) | function | Design-time synthesis (not PWM-rate) | Synthesize a SISO LQG + PR design with internal-model compensation |
| [`LQGAnalysisModels`](inc/damp/design/synthesis.hpp#L39) | block | Design-time synthesis (not PWM-rate) | Analysis-oriented models produced from an LQG design |
| [`LQGArtifacts`](inc/damp/design/synthesis.hpp#L89) | block | Design-time synthesis (not PWM-rate) | Synthesis artifact bundle: design + analysis models + runtime bundle |
| [`LQGI`](inc/damp/controllers/lqgi.hpp#L167) | block | Runtime controllers | Linear-Quadratic-Gaussian-Integral (LQGI) controller |
| [`lqgi_bundle`](inc/damp/design/synthesis.hpp#L403) | function | Design-time synthesis (not PWM-rate) | Synthesize the full LQGI servo artifact bundle in one call |
| [`LQGIAnalysisModels`](inc/damp/design/synthesis.hpp#L152) | block | Design-time synthesis (not PWM-rate) | Analysis-oriented models produced from an LQGI design |
| [`LQGIArtifacts`](inc/damp/design/synthesis.hpp#L180) | block | Design-time synthesis (not PWM-rate) | Synthesis artifact bundle: LQGI servo design + analysis + ready-to-run controller |
| [`LQGIResult`](inc/damp/controllers/lqgi.hpp#L34) | block | Design-time synthesis (not PWM-rate) | LQGI design result |
| [`LQGPRArtifacts`](inc/damp/design/synthesis.hpp#L133) | block | Design-time synthesis (not PWM-rate) | Synthesis artifact bundle for SISO LQG + PR |
| [`LQGPRRuntimeBundle`](inc/damp/design/synthesis.hpp#L100) | block | Design-time synthesis (not PWM-rate) | Runtime bundle for SISO LQG + PR internal model compensation |
| [`lqgreg`](inc/damp/matlab.hpp#L760) | function | MATLAB®-style aliases (host) | Combine separate Kalman filter and LQR designs into an LQG controller |
| [`LQGResult`](inc/damp/controllers/lqg.hpp#L46) | block | Design-time synthesis (not PWM-rate) | LQG design result |
| [`LQGRuntimeBundle`](inc/damp/design/synthesis.hpp#L50) | block | Design-time synthesis (not PWM-rate) | Runtime bundle for LQG control |
| [`lqgtrack`](inc/damp/matlab.hpp#L772) | function | MATLAB®-style aliases (host) | Linear-Quadratic-Gaussian design with integral action for tracking |
| [`lqi`](inc/damp/matlab.hpp#L731) | function | MATLAB®-style aliases (host) | Linear-Quadratic Integral design for tracking |
| [`LQI`](inc/damp/controllers/lqi.hpp#L155) | block | Runtime controllers | Linear-Quadratic-Integral (LQI) controller |
| [`lqi_bundle`](inc/damp/design/synthesis.hpp#L385) | function | Design-time synthesis (not PWM-rate) | Synthesize the full LQI servo artifact bundle in one call |
| [`LQIAnalysisModels`](inc/damp/design/synthesis.hpp#L143) | block | Design-time synthesis (not PWM-rate) | Analysis-oriented models produced from an LQI design |
| [`LQIArtifacts`](inc/damp/design/synthesis.hpp#L165) | block | Design-time synthesis (not PWM-rate) | Synthesis artifact bundle: LQI servo design + analysis + ready-to-run controller |
| [`LQIResult`](inc/damp/controllers/lqi.hpp#L48) | block | Design-time synthesis (not PWM-rate) | LQI design result |
| [`lqr`](inc/damp/matlab.hpp#L610) | function | MATLAB®-style aliases (host) | Continuous-time LQR design (MATLAB®'s lqr) |
| [`lqr_gain`](inc/damp/design/riccati.hpp#L714) | function | Design-time synthesis (not PWM-rate) | Optimal LQR state-feedback gain from a Riccati solution |
| [`LQRCost`](inc/damp/controllers/lqr.hpp#L148) | block | Runtime controllers | Discretized LQR cost weights (Q, R, N) for a sampled-data problem |
| [`lqrd`](inc/damp/controllers/lqr.hpp#L323) | function | Design-time synthesis (not PWM-rate) | Sampled-data LQR from continuous plant (MATLAB®-style short name) (+1 more overload) |
| [`lqrd`](inc/damp/matlab.hpp#L640) | function | MATLAB®-style aliases (host) | Design discrete LQR from continuous-time system via discretization (+1 more overload) |
| [`LQRResult`](inc/damp/controllers/lqr.hpp#L56) | block | Design-time synthesis (not PWM-rate) | Linear-Quadratic Regulator design result |
| [`lqry`](inc/damp/matlab.hpp#L687) | function | MATLAB®-style aliases (host) | Output-weighted continuous LQR (state cost Q = Cᵀ Q_y C) (+1 more overload) |
| [`lsim`](inc/damp/analysis/time_response.hpp#L257) | function | Frequency-domain analysis (host) | Forced time response of a (MIMO) state-space system to an input signal (+1 more overload) |
| [`LsimInfo`](inc/damp/analysis/time_response.hpp#L484) | block | Frequency-domain analysis (host) | Transient characteristics of an arbitrary response signal |
| [`lsiminfo`](inc/damp/analysis/time_response.hpp#L503) | function | Frequency-domain analysis (host) | Compute transient characteristics from an output/time signal |
| [`lsimplot`](inc/damp/simulation/plot_plotly.hpp#L749) | function | Simulation / SIL harness (host) | Plot a forced (lsim) simulation, one trace per output |
| [`LsimResult`](inc/damp/analysis/time_response.hpp#L66) | block | Frequency-domain analysis (host) | Result of a single-trajectory simulation: time, output, and state history |
| [`lu_decomposition`](inc/damp/matrix/decomposition.hpp#L131) | function | Linear algebra | LU decomposition with partial pivoting |
| [`lu_solve`](inc/damp/matrix/solve.hpp#L187) | function | Linear algebra | Solve AX = B via LU with partial pivoting |
| [`luenberger`](inc/damp/estimation/luenberger.hpp#L138) | function | Observers & estimators | Design a Luenberger observer by robust pole placement (matrix form) (+3 more overloads) |
| [`Luenberger`](inc/damp/estimation/luenberger.hpp#L400) | block | Observers & estimators | Luenberger state observer (runtime) |
| [`LuenbergerResult`](inc/damp/estimation/luenberger.hpp#L92) | block | Observers & estimators | Luenberger observer design result |
| [`Lut1D`](inc/damp/toolbox/lookup.hpp#L108) | block | Embedded helpers (controls-adjacent utilities) | 1-D interpolating lookup table over monotonic breakpoints |
| [`Lut2D`](inc/damp/toolbox/lookup.hpp#L200) | block | Embedded helpers (controls-adjacent utilities) | 2-D bilinear interpolating lookup table over a regular grid |
| [`Lut3D`](inc/damp/toolbox/lookup.hpp#L274) | block | Embedded helpers (controls-adjacent utilities) | 3-D trilinear interpolating lookup table over a regular grid |
| [`lut_segment`](inc/damp/toolbox/lookup.hpp#L63) | function | Embedded helpers (controls-adjacent utilities) | Index of the interpolation segment containing x |
| [`lyap`](inc/damp/design/lyapunov.hpp#L105) | function | Design-time synthesis (not PWM-rate) | Solve the continuous-time Lyapunov equation A X + X Aᵀ + Q = 0 |
| [`lyap`](inc/damp/matlab.hpp#L1041) | function | MATLAB®-style aliases (host) | MATLAB® alias for the continuous Lyapunov solve AX+XAᵀ+Q=0 |
| [`MadgwickFilter`](inc/damp/estimation/sensor_fusion.hpp#L105) | block | Observers & estimators | Madgwick gradient-descent AHRS filter |
| [`mag2db`](inc/damp/math/math.hpp#L459) | function | Scalar math, complex & frames | Magnitude to decibels, 20·log10(mag) |
| [`MahonyFilter`](inc/damp/estimation/sensor_fusion.hpp#L163) | block | Observers & estimators | Mahony nonlinear complementary filter with PI correction |
| [`make1DOF`](inc/damp/matlab.hpp#L227) | function | MATLAB®-style aliases (host) | Force 1-DOF setpoint weights on a PID design result (b=c=1) |
| [`make2DOF`](inc/damp/matlab.hpp#L243) | function | MATLAB®-style aliases (host) | Apply 2-DOF setpoint weights on a PID design result |
| [`margin`](inc/damp/matlab.hpp#L878) | function | MATLAB®-style aliases (host) | Gain and phase margins of a SISO loop over a frequency grid |
| [`MarginResult`](inc/damp/matlab.hpp#L857) | block | MATLAB®-style aliases (host) | Gain/phase margins and their crossover frequencies |
| [`margins_from_frf`](inc/damp/estimation/frequency_response.hpp#L771) | function | Observers & estimators | Gain / phase margins from a fixed on-target FRF table |
| [`Matrix`](inc/damp/matrix/core.hpp#L141) | block | Linear algebra | Fixed-size, stack-allocated matrix for linear algebra operations |
| [`MatrixLike`](inc/damp/matrix/matrix_traits.hpp#L70) | concept | Linear algebra | Concept for any type that provides 2D matrix-like element access |
| [`MatrixLikeOf`](inc/damp/matrix/matrix_traits.hpp#L80) | concept | Linear algebra | Concept for a MatrixLike type with specific dimensions |
| [`max`](inc/damp/matrix/core.hpp#L970) | function | Linear algebra | Addition of two MatrixLike types (with broadcasting support) (+1 more overload) |
| [`max_iq`](inc/damp/motor/spm.hpp#L111) | function | Motor control pack (if present) | Max positive i_q for an SPM (i_d = 0) at electrical speed |
| [`max_torque_at_speed`](inc/damp/motor/spm.hpp#L172) | function | Motor control pack (if present) | SPM max-torque point at a mechanical speed (i_d = 0 only) |
| [`MeasJacobian`](inc/damp/estimation/ekf.hpp#L69) | block | Observers & estimators | Measurement prediction result from the user's observation function |
| [`mecanum_drive`](inc/damp/toolbox/io.hpp#L514) | function | Embedded helpers (controls-adjacent utilities) | Mecanum (H-layout, 45° rollers) mixer — alias for holonomic_drive |
| [`mechanize_step_from_corrected`](inc/damp/estimation/ins_mechanization.hpp#L173) | function | Observers & estimators | One dead-reckoning step from an IMU sample (strapdown mechanization) |
| [`MedianFilter`](inc/damp/filters/median.hpp#L34) | block | Filters & signal conditioning | Sliding-window median filter — nonlinear spike/outlier rejection |
| [`mhe`](inc/damp/estimation/mhe.hpp#L166) | function | Observers & estimators | Synthesize a constrained moving-horizon estimator |
| [`MHE`](inc/damp/estimation/mhe.hpp#L283) | block | Observers & estimators | Runtime moving-horizon estimator (fixed per-tick iteration budget) (+1 more overload) |
| [`MHEArtifacts`](inc/damp/estimation/mhe.hpp#L110) | block | Observers & estimators | Window QP data produced by mhe(), consumed by damp::MHE |
| [`MHEConstraints`](inc/damp/estimation/mhe.hpp#L87) | block | Observers & estimators | Box constraints on the estimated states (applied to every window state) |
| [`middlebrook`](inc/damp/analysis/frequency.hpp#L1083) | function | Frequency-domain analysis (host) | Middlebrook stability analysis for cascaded source-load systems (+1 more overload) |
| [`MiddlebrookResult`](inc/damp/analysis/frequency.hpp#L896) | block | Frequency-domain analysis (host) | Result of Middlebrook minor loop gain analysis |
| [`minimal_realization`](inc/damp/design/minreal.hpp#L428) | function | Design-time synthesis (not PWM-rate) | Descriptive alias for minreal |
| [`MinimalRealizationResult`](inc/damp/design/minreal.hpp#L84) | block | Design-time synthesis (not PWM-rate) | Minimal realization result |
| [`minmax`](inc/damp/backend.hpp#L164) | function | Core, configuration & backend vocabulary | Ordered {min, max} pair returned by value (+1 more overload) |
| [`minreal`](inc/damp/design/minreal.hpp#L284) | function | Design-time synthesis (not PWM-rate) | Minimal realization — cancel uncontrollable and unobservable modes |
| [`minreal`](inc/damp/matlab.hpp#L257) | function | MATLAB®-style aliases (host) | MATLAB® short alias for design::minreal |
| [`minreal_zpk`](inc/damp/systems/zpk.hpp#L867) | function | LTI systems (SS / TF / ZPK / discretize) | Cancel matching pole-zero pairs on a ZPK model |
| [`mix`](inc/damp/systems/state_space.hpp#L572) | function | LTI systems (SS / TF / ZPK / discretize) | Fold plant outputs into one input: $`v = g u + \sum k_i y_i`$ (+1 more overload) |
| [`mix_outputs`](inc/damp/systems/state_space.hpp#L640) | function | LTI systems (SS / TF / ZPK / discretize) | Linear combination of outputs: $`y = \sum \alpha_j y_j`$ |
| [`ModeExtractorConfig`](inc/damp/estimation/frequency_response.hpp#L411) | block | Observers & estimators | Configuration for FRF peak / valley extraction |
| [`ModeExtractorResult`](inc/damp/estimation/frequency_response.hpp#L428) | block | Observers & estimators | Result of reducing an FRF table to at most MaxModes peaks (resonances) |
| [`ModelReductionMethod`](inc/damp/design/model_reduction.hpp#L58) | enum | Design-time synthesis (not PWM-rate) | Method for eliminating states in modred |
| [`modred`](inc/damp/design/model_reduction.hpp#L518) | function | Design-time synthesis (not PWM-rate) | Model reduction by truncation or DC-matched residualization |
| [`modulation_duty_cycles`](inc/damp/motor/modulation.hpp#L396) | function | Motor control pack (if present) | Carrier-based VSI duty cycles from an αβ voltage command |
| [`modulation_index`](inc/damp/motor/modulation.hpp#L165) | function | Motor control pack (if present) | Modulation index relative to the linear SVPWM circle |
| [`modulation_zero_sequence`](inc/damp/motor/modulation.hpp#L323) | function | Motor control pack (if present) | Zero-sequence for any PwmScheme |
| [`Modulator`](inc/damp/motor/modulation.hpp#L493) | block | Motor control pack (if present) | Thin scheme-holding wrapper for FOC / deploy paths |
| [`motor_constant`](inc/damp/motor/spm.hpp#L63) | function | Motor control pack (if present) | Motor constant Kₘ = Kₜ / √R [Nm/√W] |
| [`MovingAverage`](inc/damp/filters/moving_average.hpp#L62) | block | Filters & signal conditioning | Moving-average (boxcar) filter — also a DC-preserving harmonic-notch comb |
| [`MPC`](inc/damp/controllers/mpc.hpp#L923) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its artifacts: MPC controller{art}; (+1 more overload) |
| [`mpc`](inc/damp/controllers/offset_free_mpc.hpp#L132) | function | Design-time synthesis (not PWM-rate) | Synthesize an offset-free constrained MPC (controller + estimator) |
| [`MPC`](inc/damp/controllers/mpc.hpp#L776) | block | Runtime controllers | Runtime constrained MPC controller (fixed per-tick iteration budget) |
| [`MPCAnalysisModels`](inc/damp/controllers/mpc.hpp#L671) | block | Runtime controllers | LTI models of the unconstrained MPC loop, for margin/robustness analysis |
| [`MPCArtifacts`](inc/damp/controllers/mpc.hpp#L240) | block | Runtime controllers | Condensed-QP data produced by mpc(), consumed by damp::MPC |
| [`MPCConstraints`](inc/damp/controllers/mpc.hpp#L202) | block | Runtime controllers | Box constraints for mpc() |
| [`MPCHorizonSuggestion`](inc/damp/controllers/mpc.hpp#L566) | block | Runtime controllers | Advisory NP/NC horizon values from suggest_mpc_horizon() |
| [`MPCWeights`](inc/damp/controllers/mpc.hpp#L167) | block | Runtime controllers | Cost weights for mpc() |
| [`mstogi`](inc/damp/filters/sogi.hpp#L117) | function | Filters & signal conditioning | Mixed Second/Third-Order Generalized Integrator (MSTOGI) (+1 more overload) |
| [`MSTOGI`](inc/damp/filters/sogi.hpp#L242) | block | Filters & signal conditioning | Runtime MSTOGI with exact resonator and forward-Euler washout |
| [`multi_sine`](inc/damp/estimation/excitation/multi_sine.hpp#L133) | function | Observers & estimators | Build a multi-sine design payload from a configuration |
| [`MultiPRController`](inc/damp/controllers/pr.hpp#L369) | block | Runtime controllers | Multi-harmonic PR Controller |
| [`MultiRateConfig`](inc/damp/simulation/multirate.hpp#L81) | block | Simulation / SIL harness (host) | Configuration for multi-rate closed-loop runs |
| [`MultiSine`](inc/damp/estimation/excitation/multi_sine.hpp#L200) | block | Observers & estimators | Sum-of-tones multi-sine runtime generator |
| [`MultiSineConfig`](inc/damp/estimation/excitation/multi_sine.hpp#L62) | block | Observers & estimators | Configuration for fixed-component multi-sine excitation |
| [`MultiSineResult`](inc/damp/estimation/excitation/multi_sine.hpp#L104) | block | Observers & estimators | Multi-sine design payload |
| [`MUX`](inc/damp/toolbox/iec61131.hpp#L617) | function | Embedded helpers (controls-adjacent utilities) | MUX (IEC 61131-3 multiplexer): select input k of N (0-based) |
| [`NavFrame`](inc/damp/estimation/ins_mechanization.hpp#L77) | enum | Observers & estimators | Local-level navigation frame for strapdown mechanization |
| [`nearbyint`](inc/damp/math/math.hpp#L345) | function | Scalar math, complex & frames | Round to nearest integer; ties to even (IEEE default / `FE_TONEAREST`) |
| [`ned_from_enu`](inc/damp/estimation/ins_mechanization.hpp#L106) | function | Observers & estimators | Map an ENU vector into NED (axis permute) |
| [`negative_sequence_ab`](inc/damp/filters/pll.hpp#L244) | function | Filters & signal conditioning | Instantaneous negative-sequence αβ from a quadrature signal pair |
| [`nichols`](inc/damp/analysis/frequency.hpp#L683) | function | Frequency-domain analysis (host) | Build Nichols points from existing Bode data (+2 more overloads) |
| [`nicholsplot`](inc/damp/simulation/plot_plotly.hpp#L936) | function | Simulation / SIL harness (host) | Plot a Nichols chart (open-loop phase vs magnitude) with M-circle grid |
| [`NicholsPoint`](inc/damp/analysis/frequency.hpp#L654) | block | Frequency-domain analysis (host) | Single-point Nichols chart sample (open-loop phase vs magnitude) |
| [`NicholsResult`](inc/damp/analysis/frequency.hpp#L666) | block | Frequency-domain analysis (host) | Nichols chart locus across a frequency sweep |
| [`NlmsFilter`](inc/damp/filters/fir.hpp#L546) | block | Filters & signal conditioning | Normalized LMS adaptive FIR filter |
| [`norm`](inc/damp/matlab.hpp#L1111) | function | MATLAB®-style aliases (host) | MATLAB® alias for the H2 system norm norm(sys,2) |
| [`norm_h2`](inc/damp/analysis/norms.hpp#L76) | function | Frequency-domain analysis (host) | H2 norm of a state-space system |
| [`norm_hinf`](inc/damp/analysis/norms.hpp#L134) | function | Frequency-domain analysis (host) | H∞ norm of a state-space system: sup_ω σ̄(G(jω)) |
| [`notch`](inc/damp/filters/iir_design.hpp#L789) | function | Filters & signal conditioning | Second-order band-reject (notch) filter |
| [`notch_continuous`](inc/damp/filters/iir_design.hpp#L406) | function | Filters & signal conditioning | Second-order band-reject (notch) filter design (continuous-time) |
| [`null`](inc/damp/matlab.hpp#L466) | function | MATLAB®-style aliases (host) | MATLAB® short alias for an orthonormal null-space basis |
| [`null_space`](inc/damp/matrix/svd.hpp#L350) | function | Linear algebra | Orthonormal basis for the null space {x : A·x = 0} via SVD |
| [`NullSpace`](inc/damp/matrix/svd.hpp#L331) | block | Linear algebra | Orthonormal basis for the null space (kernel) of a matrix |
| [`nyquist`](inc/damp/analysis/frequency.hpp#L497) | function | Frequency-domain analysis (host) | Compute Nyquist data for a SISO state-space system (+1 more overload) |
| [`nyquistplot`](inc/damp/simulation/plot_plotly.hpp#L821) | function | Simulation / SIL harness (host) | Plot a Nyquist locus with the -1 critical point marked |
| [`NyquistPoint`](inc/damp/analysis/frequency.hpp#L369) | block | Frequency-domain analysis (host) | Single-point Nyquist response data |
| [`NyquistResult`](inc/damp/analysis/frequency.hpp#L381) | block | Frequency-domain analysis (host) | Nyquist response data across a frequency sweep |
| [`observability_gramian`](inc/damp/design/stability.hpp#L146) | function | Design-time synthesis (not PWM-rate) | Continuous/discrete observability Gramian W_o |
| [`observability_matrix`](inc/damp/design/stability.hpp#L82) | function | Design-time synthesis (not PWM-rate) | Compute the observability matrix [C; CA; CA²; ...; CA^(N-1)] |
| [`obsv`](inc/damp/matlab.hpp#L281) | function | MATLAB®-style aliases (host) | MATLAB® short alias for observability_matrix (+1 more overload) |
| [`OffDelayTimer`](inc/damp/toolbox/logic.hpp#L133) | block | Embedded helpers (controls-adjacent utilities) | Off-delay timer: output goes true immediately when in is true and stays true until in has been false continuously for delay |
| [`OffsetFreeMPC`](inc/damp/controllers/offset_free_mpc.hpp#L262) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its artifacts: OffsetFreeMPC controller{art}; |
| [`OffsetFreeMPC`](inc/damp/controllers/offset_free_mpc.hpp#L198) | block | Runtime controllers | Runtime offset-free MPC: constrained MPC + disturbance-augmented Kalman filter |
| [`OffsetFreeMPCArtifacts`](inc/damp/controllers/offset_free_mpc.hpp#L84) | block | Runtime controllers | Combined MPC + disturbance-augmented Kalman design, consumed by damp::OffsetFreeMPC |
| [`omniwheel_drive`](inc/damp/toolbox/io.hpp#L520) | function | Embedded helpers (controls-adjacent utilities) | Omni-wheel (X-layout) mixer — alias for holonomic_drive |
| [`OnDelayTimer`](inc/damp/toolbox/logic.hpp#L96) | block | Embedded helpers (controls-adjacent utilities) | On-delay timer: output goes true once in has been held true continuously for delay; drops immediately when in goes false |
| [`one_norm`](inc/damp/matrix/functions.hpp#L64) | function | Linear algebra | One-norm ‖A‖₁: maximum absolute column sum |
| [`OptimalJordanPlacement`](inc/damp/design/pole_placement.hpp#L1172) | block | Design-time synthesis (not PWM-rate) | Result of optimized arbitrary pole placement (place_jordan_optimal) |
| [`OutputFeedbackController`](inc/damp/concepts.hpp#L87) | concept | Core, configuration & backend vocabulary | Vector output-feedback controller: u = control(r, y), self-contained tick |
| [`pade`](inc/damp/matlab.hpp#L987) | function | MATLAB®-style aliases (host) | First-order Padé approximation of pure delay e^{−sT} (+1 more overload) |
| [`pade2`](inc/damp/matlab.hpp#L1009) | function | MATLAB®-style aliases (host) | Second-order Padé approximation of pure delay e^{−sT} (+1 more overload) |
| [`pade_delay_1st`](inc/damp/filters/iir_design.hpp#L528) | function | Filters & signal conditioning | First-order Pade approximation of time delay (+1 more overload) |
| [`pade_delay_2nd`](inc/damp/filters/iir_design.hpp#L572) | function | Filters & signal conditioning | Second-order Pade approximation of time delay (+1 more overload) |
| [`panel_color`](inc/damp/simulation/plot_plotly.hpp#L314) | function | Simulation / SIL harness (host) | Plotly default colorway entry for index i within a single subplot |
| [`panel_legend`](inc/damp/simulation/plot_plotly.hpp#L297) | function | Simulation / SIL harness (host) | Per-panel legend box just right of a stacked subplot band |
| [`panel_legend_id`](inc/damp/simulation/plot_plotly.hpp#L283) | function | Simulation / SIL harness (host) | Legend id for subplot row (1-based): "legend", "legend2", … |
| [`panel_x_domain`](inc/damp/simulation/plot_plotly.hpp#L273) | function | Simulation / SIL harness (host) | Horizontal domain for stacked multi-panel x-axes |
| [`parallel`](inc/damp/systems/state_space.hpp#L256) | function | LTI systems (SS / TF / ZPK / discretize) | Parallel connection (shared input, summed outputs) |
| [`ParameterDriftMonitor`](inc/damp/estimation/parameter_estimation.hpp#L252) | block | Observers & estimators | Residual-based drift detector: has the plant moved away from the model? |
| [`ParameterEstimator`](inc/damp/concepts.hpp#L142) | concept | Core, configuration & backend vocabulary | Online grey-box parameter estimator: physical parameters + gating |
| [`park_transform`](inc/damp/math/transforms.hpp#L331) | function | Scalar math, complex & frames | Park transform (αβ → dq) |
| [`peaking`](inc/damp/filters/iir_design.hpp#L878) | function | Filters & signal conditioning | Peaking (bell) EQ filter: boost or cut a band around f0 |
| [`Periodic`](inc/damp/toolbox/timing.hpp#L125) | block | Embedded helpers (controls-adjacent utilities) | Periodic trigger — fires once per elapsed period |
| [`phase_margin_from_damping_ratio`](inc/damp/design/pid_design.hpp#L788) | function | Design-time synthesis (not PWM-rate) | Approximate phase margin from damping ratio |
| [`phase_margin_unwrapped`](inc/damp/analysis/frequency.hpp#L166) | function | Frequency-domain analysis (host) | Find phase margin using unwrapped phase trajectory |
| [`pi_pole_placement_first_order`](inc/damp/design/pid_design.hpp#L1005) | function | Design-time synthesis (not PWM-rate) | PI gains that place the closed-loop poles of a first-order plant (+1 more overload) |
| [`pid`](inc/damp/controllers/pid.hpp#L275) | function | Design-time synthesis (not PWM-rate) | 2-DOF continuous PID controller design |
| [`pid`](inc/damp/matlab.hpp#L168) | function | MATLAB®-style aliases (host) | MATLAB®-style parallel-form continuous PID constructor |
| [`pid_from_bandwidth`](inc/damp/design/pid_design.hpp#L435) | function | Design-time synthesis (not PWM-rate) | Design PID from desired bandwidth and phase margin |
| [`pid_from_performance_spec`](inc/damp/design/pid_design.hpp#L851) | function | Design-time synthesis (not PWM-rate) | Design PID directly from settling-time and overshoot targets |
| [`pid_pole_placement`](inc/damp/design/pid_design.hpp#L888) | function | Design-time synthesis (not PWM-rate) | Direct PID pole placement for a first-order-plus-dead-time model (+1 more overload) |
| [`pid_pole_placement_double_integrator`](inc/damp/design/pid_design.hpp#L1087) | function | Design-time synthesis (not PWM-rate) | PID (or PD) pole placement for a rigid inertia $`G(s) = 1/(J s^2)`$ |
| [`PIDController`](inc/damp/controllers/pid.hpp#L365) | block | Runtime controllers | Fixed-rate discrete 2-DOF PID (canonical runtime) (+2 more overloads) |
| [`PIDMode`](inc/damp/controllers/pid.hpp#L300) | enum | Runtime controllers | Compile-time selection of the PID control-law structure |
| [`PIDPerformanceSpec`](inc/damp/design/pid_design.hpp#L833) | block | Design-time synthesis (not PWM-rate) | Time-domain performance targets for quick PID synthesis |
| [`PIDResult`](inc/damp/controllers/pid.hpp#L100) | block | Design-time synthesis (not PWM-rate) | 2-DOF continuous-time PID controller design result |
| [`PIDRuntimeMode`](inc/damp/controllers/pid.hpp#L324) | enum | Runtime controllers | Runtime operating mode for PIDController / ContinuousPID |
| [`pidstd`](inc/damp/matlab.hpp#L196) | function | MATLAB®-style aliases (host) | Standard-form continuous PID constructor (1-DOF) |
| [`pidstd2`](inc/damp/matlab.hpp#L214) | function | MATLAB®-style aliases (host) | Standard-form continuous 2-DOF PID constructor |
| [`pidtune`](inc/damp/design/pid_design.hpp#L666) | function | Design-time synthesis (not PWM-rate) | Plant-aware PID tune at a crossover (+3 more overloads) |
| [`pidtune`](inc/damp/matlab.hpp#L792) | function | MATLAB®-style aliases (host) | Plant-aware PID tune (MATLAB® pidtune spelling) |
| [`PidTuneSpec`](inc/damp/design/pid_design.hpp#L509) | block | Design-time synthesis (not PWM-rate) | Knobs for plant-aware pidtune |
| [`PIDType`](inc/damp/design/pid_design.hpp#L50) | enum | Design-time synthesis (not PWM-rate) | PID controller type selection for tuning methods |
| [`PiecewiseLTI`](inc/damp/simulation/hybrid.hpp#L94) | block | Simulation / SIL harness (host) | Plant: NModes continuous LTI pieces (same state/input/output dimensions) |
| [`pinv`](inc/damp/matlab.hpp#L453) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the Moore–Penrose pseudoinverse |
| [`place`](inc/damp/design/pole_placement.hpp#L88) | function | Design-time synthesis (not PWM-rate) | Robust multi-input pole placement (Kautsky–Nichols–Van Dooren, real poles) (+5 more overloads) |
| [`place`](inc/damp/matlab.hpp#L589) | function | MATLAB®-style aliases (host) | Robust multi-input pole placement (MATLAB®'s place) |
| [`place_discrete`](inc/damp/design/pole_placement.hpp#L581) | function | Design-time synthesis (not PWM-rate) | Discrete plant + z-plane poles → discrete \(K\) (+3 more overloads) |
| [`place_jordan`](inc/damp/design/pole_placement.hpp#L1138) | function | Design-time synthesis (not PWM-rate) | Exact pole placement with an arbitrary Jordan structure (Schmid–Ntogramatzidis–Nguyen–Pandey / Klein–Moore parametric form) |
| [`place_jordan_optimal`](inc/damp/design/pole_placement.hpp#L1237) | function | Design-time synthesis (not PWM-rate) | Robust / minimum-gain arbitrary pole placement (Schmid et al., Methods 1–2) |
| [`place_observer`](inc/damp/design/pole_placement.hpp#L659) | function | Design-time synthesis (not PWM-rate) | Continuous observer gain: place eigenvalues of \(A - LC\) (s-plane poles) (+5 more overloads) |
| [`place_observer_discrete`](inc/damp/design/pole_placement.hpp#L769) | function | Design-time synthesis (not PWM-rate) | Discrete plant + z-plane poles → discrete observer gain \(L\) (+3 more overloads) |
| [`plot_bode`](inc/damp/simulation/plot_plotly.hpp#L422) | function | Simulation / SIL harness (host) | Plot Bode magnitude and phase as subplots |
| [`plot_line`](inc/damp/simulation/plot_plotly.hpp#L491) | function | Simulation / SIL harness (host) | Simple line plot of time vs value |
| [`plot_simulation`](inc/damp/simulation/plot_plotly.hpp#L343) | function | Simulation / SIL harness (host) | Plot simulation results with subplots for states, outputs, and inputs |
| [`plot_step`](inc/damp/simulation/plot_plotly.hpp#L621) | function | Simulation / SIL harness (host) | Plot step response data |
| [`plot_xy`](inc/damp/simulation/plot_plotly.hpp#L556) | function | Simulation / SIL harness (host) | Multi-series planar scatter (true 2-D path / locus) (+1 more overload) |
| [`plots_root_for`](inc/damp/simulation/plot_plotly.hpp#L118) | function | Simulation / SIL harness (host) | Walk up from dir until a directory named "plots" is found |
| [`pole`](inc/damp/matlab.hpp#L843) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the open-loop poles of a system |
| [`PoleInfo`](inc/damp/analysis/poles.hpp#L75) | block | Frequency-domain analysis (host) | Natural frequency and damping ratio for each pole |
| [`poles`](inc/damp/analysis/poles.hpp#L41) | function | Frequency-domain analysis (host) | Compute open-loop poles (eigenvalues of A matrix) |
| [`PoleZeroMap`](inc/damp/analysis/poles.hpp#L115) | block | Frequency-domain analysis (host) | Poles and zeros of a system, for pole-zero plotting |
| [`poly_horner`](inc/damp/toolbox/scaling.hpp#L134) | function | Embedded helpers (controls-adjacent utilities) | Evaluate a polynomial at x by Horner's method |
| [`poly_roots`](inc/damp/analysis/poles.hpp#L128) | function | Frequency-domain analysis (host) | Roots of a polynomial given in ascending powers (MATLAB® `roots`, reversed order) |
| [`poly_roots`](inc/damp/systems/zpk.hpp#L118) | function | LTI systems (SS / TF / ZPK / discretize) | Compute roots of an ascending-power real polynomial |
| [`PolyRootsResult`](inc/damp/systems/zpk.hpp#L89) | block | LTI systems (SS / TF / ZPK / discretize) | Roots of a real polynomial given in ascending powers |
| [`Pose`](inc/damp/kinematics/pose.hpp#L60) | block | Kinematics / pose | Rigid-body pose: a translation and an orientation (unit quaternion) |
| [`positive_sequence_ab`](inc/damp/filters/pll.hpp#L223) | function | Filters & signal conditioning | Instantaneous positive-sequence αβ from a quadrature signal pair |
| [`pow`](inc/damp/matrix/functions.hpp#L570) | function | Linear algebra | Integer matrix power via binary exponentiation (+1 more overload) |
| [`pow`](inc/damp/math/math.hpp#L268) | function | Scalar math, complex & frames | Power function, base^exponent (+1 more overload) |
| [`pr`](inc/damp/controllers/pr.hpp#L153) | function | Design-time synthesis (not PWM-rate) | Design a Proportional-Resonant controller |
| [`pr_harmonics`](inc/damp/controllers/pr.hpp#L185) | function | Design-time synthesis (not PWM-rate) | Design multiple-harmonic PR controller gains |
| [`prbs`](inc/damp/estimation/excitation/prbs.hpp#L158) | function | Observers & estimators | Build a PRBS design payload from a configuration |
| [`PRBS`](inc/damp/estimation/excitation/prbs.hpp#L181) | block | Observers & estimators | Maximal-length PRBS runtime generator |
| [`PRBSConfig`](inc/damp/estimation/excitation/prbs.hpp#L81) | block | Observers & estimators | Configuration for maximal-length pseudo-random binary excitation |
| [`PRBSResult`](inc/damp/estimation/excitation/prbs.hpp#L125) | block | Observers & estimators | PRBS design payload |
| [`PRController`](inc/damp/controllers/pr.hpp#L222) | block | Runtime controllers | Discrete Proportional-Resonant Controller |
| [`project_affine`](inc/damp/controllers/action_governor.hpp#L231) | function | Design-time synthesis (not PWM-rate) | Project u_des onto the polyhedron A u ≤ b (Euclidean) |
| [`project_box`](inc/damp/controllers/action_governor.hpp#L78) | function | Design-time synthesis (not PWM-rate) | Euclidean projection of u onto the axis-aligned box [umin, umax] (+1 more overload) |
| [`ProportionalCurrentMap`](inc/damp/toolbox/io.hpp#L139) | block | Embedded helpers (controls-adjacent utilities) | Affine map from a normalized command to a proportional current (mA) |
| [`PRResult`](inc/damp/controllers/pr.hpp#L72) | block | Design-time synthesis (not PWM-rate) | Proportional-Resonant controller design result |
| [`pseudo_inverse`](inc/damp/matrix/svd.hpp#L303) | function | Linear algebra | Moore–Penrose pseudoinverse A⁺ via SVD |
| [`PulseTimer`](inc/damp/toolbox/logic.hpp#L170) | block | Embedded helpers (controls-adjacent utilities) | Pulse timer (non-retriggerable): a rising edge of in emits a fixed |
| [`PwmScheme`](inc/damp/motor/modulation.hpp#L79) | enum | Motor control pack (if present) | Carrier-based three-phase VSI modulation scheme |
| [`pzmap`](inc/damp/analysis/poles.hpp#L154) | function | Frequency-domain analysis (host) | Pole-zero map of a SISO transfer function (MATLAB® `pzmap(tf)`) (+2 more overloads) |
| [`pzplot`](inc/damp/simulation/plot_plotly.hpp#L848) | function | Simulation / SIL harness (host) | Plot a pole-zero map on the complex plane (poles as ×, zeros as ○) |
| [`QPResult`](inc/damp/design/qp.hpp#L82) | block | Design-time synthesis (not PWM-rate) | Result of a dense QP solve |
| [`QPStatus`](inc/damp/design/qp.hpp#L51) | enum | Design-time synthesis (not PWM-rate) | Termination status of a QP solve |
| [`qr_decompose`](inc/damp/matrix/decomposition.hpp#L223) | function | Linear algebra | Thin QR via modified Gram–Schmidt |
| [`QRDecomposition`](inc/damp/matrix/decomposition.hpp#L191) | block | Linear algebra | Thin QR factorization A = QR (modified Gram–Schmidt) |
| [`QuadMode`](inc/damp/toolbox/encoder.hpp#L47) | enum | Embedded helpers (controls-adjacent utilities) | Quadrature decode resolution (edges counted per A/B cycle) |
| [`quadprog`](inc/damp/matlab.hpp#L1062) | function | MATLAB®-style aliases (host) | MATLAB® alias for the dense inequality-constrained QP solve (+1 more overload) |
| [`quadratic_form`](inc/damp/matrix/core.hpp#L1071) | function | Linear algebra | Symmetric congruence (quadratic) form  S = M X Mᵀ |
| [`QuadratureDecoder`](inc/damp/toolbox/encoder.hpp#L69) | block | Embedded helpers (controls-adjacent utilities) | Software A/B quadrature decoder with optional index |
| [`Quaternion`](inc/damp/math/geometry.hpp#L381) | block | Scalar math, complex & frames | Unit quaternion rotation (w, x, y, z) (Hamilton product) |
| [`R_TRIG`](inc/damp/toolbox/iec61131.hpp#L118) | block | Embedded helpers (controls-adjacent utilities) | R_TRIG (Rising Edge Trigger) |
| [`rad2deg`](inc/damp/math/math.hpp#L481) | function | Scalar math, complex & frames | Radians to degrees, rad·180/π |
| [`ramp`](inc/damp/estimation/excitation/ramp.hpp#L94) | function | Observers & estimators | Build a ramp design payload from a configuration |
| [`Ramp`](inc/damp/estimation/excitation/ramp.hpp#L110) | block | Observers & estimators | Rate-limited ramp runtime generator |
| [`RampConfig`](inc/damp/estimation/excitation/ramp.hpp#L32) | block | Observers & estimators | Configuration for a slew-rate-limited ramp excitation |
| [`RampResult`](inc/damp/estimation/excitation/ramp.hpp#L64) | block | Observers & estimators | Ramp design payload |
| [`RangeMonitor`](inc/damp/toolbox/conditioning.hpp#L430) | block | Embedded helpers (controls-adjacent utilities) | Analog-input range/fault monitor (NAMUR NE43 pattern) |
| [`rank`](inc/damp/design/stability.hpp#L163) | function | Design-time synthesis (not PWM-rate) | Compute rank of a matrix via Gaussian elimination with partial pivoting |
| [`rank`](inc/damp/matrix/functions.hpp#L259) | function | Linear algebra | Matrix rank via Gaussian elimination with partial pivoting |
| [`rank_from_svd`](inc/damp/matrix/svd.hpp#L271) | function | Linear algebra | Numerical rank from a precomputed SVD result |
| [`reduced_luenberger`](inc/damp/estimation/luenberger.hpp#L271) | function | Observers & estimators | Design a reduced-order (Gopinath) observer by pole placement (matrix form) (+3 more overloads) |
| [`ReducedLuenberger`](inc/damp/estimation/luenberger.hpp#L507) | block | Observers & estimators | Reduced-order (Gopinath) state observer (runtime) |
| [`ReducedLuenbergerResult`](inc/damp/estimation/luenberger.hpp#L217) | block | Observers & estimators | Reduced-order (Gopinath) observer design result |
| [`RedundantAnalogPair`](inc/damp/toolbox/io.hpp#L90) | block | Embedded helpers (controls-adjacent utilities) | Dual-channel analog input with same-slope or opposite-slope cross-check |
| [`reg`](inc/damp/matlab.hpp#L526) | function | MATLAB®-style aliases (host) | Form dynamic regulator from system, state-feedback gain, and estimator gain |
| [`relay_autotune`](inc/damp/estimation/relay_autotune.hpp#L213) | function | Observers & estimators | Build a validated relay-autotune design payload |
| [`RelayAutotuneConfig`](inc/damp/estimation/relay_autotune.hpp#L120) | block | Observers & estimators | Configuration for the relay-feedback autotuning experiment |
| [`RelayAutotuneOutput`](inc/damp/estimation/relay_autotune.hpp#L238) | block | Observers & estimators | Per-tick output of RelayAutotuner::step |
| [`RelayAutotuner`](inc/damp/estimation/relay_autotune.hpp#L257) | block | Observers & estimators | Runtime relay-feedback autotuner |
| [`RelayAutotuneResult`](inc/damp/estimation/relay_autotune.hpp#L180) | block | Observers & estimators | Relay-autotuner design payload |
| [`RelayAutotuneStatus`](inc/damp/estimation/relay_autotune.hpp#L222) | enum | Observers & estimators | Lifecycle state of a RelayAutotuner |
| [`repetitive`](inc/damp/controllers/repetitive.hpp#L171) | function | Design-time synthesis (not PWM-rate) | Synthesize a repetitive controller with a scalar robustness filter Q |
| [`repetitive_binomial`](inc/damp/controllers/repetitive.hpp#L223) | function | Design-time synthesis (not PWM-rate) | Synthesize a repetitive controller with a binomial zero-phase FIR Q |
| [`RepetitiveConfig`](inc/damp/controllers/repetitive.hpp#L70) | block | Runtime controllers | Repetitive-controller tuning + period (with optional zero-phase FIR Q) |
| [`RepetitiveController`](inc/damp/controllers/repetitive.hpp#L273) | block | Runtime controllers | Plug-in repetitive controller runtime (fixed-size internal model) |
| [`RepetitiveResult`](inc/damp/controllers/repetitive.hpp#L134) | block | Design-time synthesis (not PWM-rate) | Design result for the repetitive controller |
| [`rescale`](inc/damp/toolbox/scaling.hpp#L69) | function | Embedded helpers (controls-adjacent utilities) | Affine map of x from the input range to the output range |
| [`resonance_autotune_from_frf`](inc/damp/estimation/commissioning_frontends.hpp#L76) | function | Observers & estimators | Map FRF peaks and valleys to a fixed bank of biquad compensators |
| [`ResonanceAutotuneResult`](inc/damp/estimation/commissioning_frontends.hpp#L44) | block | Observers & estimators | Result of resonance / anti-resonance compensator design from an FRF |
| [`ResonantMode`](inc/damp/estimation/frequency_response.hpp#L384) | block | Observers & estimators | One extracted FRF extremum (peak or valley) |
| [`RisingEdge`](inc/damp/toolbox/logic.hpp#L30) | block | Embedded helpers (controls-adjacent utilities) | Rising-edge detector: true on the tick x goes false → true |
| [`RK23`](inc/damp/simulation/integrator.hpp#L902) | block | Simulation / SIL harness (host) | Bogacki-Shampine 2(3) adaptive integrator |
| [`RK3`](inc/damp/simulation/integrator.hpp#L969) | block | Simulation / SIL harness (host) | Classical 3rd-order Runge-Kutta (RK3) integrator |
| [`RK4`](inc/damp/simulation/integrator.hpp#L683) | block | Simulation / SIL harness (host) | Classical 4th-order Runge-Kutta (RK4) integrator |
| [`rlocus`](inc/damp/analysis/poles.hpp#L262) | function | Frequency-domain analysis (host) | Root locus of a SISO state-space plant over an explicit gain grid (+1 more overload) |
| [`rlocusplot`](inc/damp/simulation/plot_plotly.hpp#L1033) | function | Simulation / SIL harness (host) | Plot a root locus (closed-loop poles vs gain) on the complex plane |
| [`rls`](inc/damp/estimation/rls.hpp#L155) | function | Observers & estimators | Build scalar RLS design payload |
| [`Rls`](inc/damp/estimation/rls.hpp#L194) | block | Observers & estimators | Scalar runtime RLS estimator |
| [`rls_vector`](inc/damp/estimation/rls.hpp#L171) | function | Observers & estimators | Build vector RLS design payload |
| [`RlsConfig`](inc/damp/estimation/rls.hpp#L41) | block | Observers & estimators | Common RLS configuration |
| [`RlsResult`](inc/damp/estimation/rls.hpp#L83) | block | Observers & estimators | Scalar RLS design payload |
| [`RlsState`](inc/damp/estimation/rls.hpp#L69) | block | Observers & estimators | Scalar RLS runtime state |
| [`RlsVector`](inc/damp/estimation/rls.hpp#L289) | block | Observers & estimators | Vector runtime RLS estimator (NP parameters) |
| [`RlsVectorResult`](inc/damp/estimation/rls.hpp#L126) | block | Observers & estimators | Vector RLS design payload for N parameters |
| [`RlsVectorState`](inc/damp/estimation/rls.hpp#L112) | block | Observers & estimators | Vector RLS runtime state for N parameters |
| [`RobustExactDifferentiator`](inc/damp/filters/differentiator.hpp#L63) | block | Filters & signal conditioning | First-order robust exact differentiator (super-twisting differentiator) |
| [`RootLocusResult`](inc/damp/analysis/poles.hpp#L203) | block | Frequency-domain analysis (host) | Root-locus data: closed-loop poles along a gain grid |
| [`rotary_gearbox`](inc/damp/toolbox/actuator.hpp#L149) | function | Embedded helpers (controls-adjacent utilities) | Build a ServoAxis for a rotary joint behind a gearbox |
| [`rows`](inc/damp/estimation/ins_eskf.hpp#L348) | function | Observers & estimators | Sparse y = F x for the 15-state INS first-order structure (G = I path) |
| [`RowVec`](inc/damp/matrix/rowvec.hpp#L29) | block | Linear algebra | Row vector specialization of Matrix<1, N, T> |
| [`RowView`](inc/damp/matrix/views.hpp#L189) | block | Linear algebra | Non-owning row view of a matrix |
| [`RS`](inc/damp/toolbox/iec61131.hpp#L91) | block | Embedded helpers (controls-adjacent utilities) | RS Latch (Reset-Set Latch) |
| [`scaled_deadband`](inc/damp/toolbox/conditioning.hpp#L233) | function | Embedded helpers (controls-adjacent utilities) | Center dead zone that rescales the surviving range back to full span |
| [`schroeder_multi_sine`](inc/damp/estimation/excitation/multi_sine.hpp#L180) | function | Observers & estimators | Build a multi-sine design payload with Schroeder phases applied |
| [`SecondOrderCoeffs`](inc/damp/filters/iir_design.hpp#L65) | block | Filters & signal conditioning | DSP coefficients for second-order IIR filter |
| [`SEL`](inc/damp/toolbox/iec61131.hpp#L606) | function | Embedded helpers (controls-adjacent utilities) | SEL (IEC 61131-3 binary selection): g ? in1 : in0 |
| [`select`](inc/damp/systems/state_space.hpp#L536) | function | LTI systems (SS / TF / ZPK / discretize) | SISO channel slice $`y_{\mathrm{iy}} / u_{\mathrm{iu}}`$ |
| [`SequenceComponents`](inc/damp/math/transforms.hpp#L490) | block | Scalar math, complex & frames | Symmetrical (sequence) components of a three-phase phasor set |
| [`series`](inc/damp/systems/state_space.hpp#L199) | function | LTI systems (SS / TF / ZPK / discretize) | Series connection: sys2 follows sys1 (u → sys1 → sys2 → y) (+1 more overload) |
| [`ServoAxis`](inc/damp/toolbox/actuator.hpp#L107) | block | Embedded helpers (controls-adjacent utilities) | One servoactuator transmission: SI joint unit ⟷ drive (motor) units |
| [`ServoBank`](inc/damp/toolbox/actuator.hpp#L229) | block | Embedded helpers (controls-adjacent utilities) | A bank of ServoAxis transmissions: maps a synchronized multi-axis TrajectoryState array straight to drive commands in one call |
| [`ServoCommand`](inc/damp/toolbox/actuator.hpp#L84) | block | Embedded helpers (controls-adjacent utilities) | A drive-native servoactuator setpoint: position, velocity, torque |
| [`sgn`](inc/damp/math/math.hpp#L397) | function | Scalar math, complex & frames | Sign function — −1 if val < 0, 1 if val > 0, 0 if val == 0 |
| [`sigma`](inc/damp/analysis/frequency.hpp#L796) | function | Frequency-domain analysis (host) | Singular-value frequency response of a (possibly MIMO) state-space system (+1 more overload) |
| [`sigmaplot`](inc/damp/simulation/plot_plotly.hpp#L990) | function | Simulation / SIL harness (host) | Plot singular-value frequency response (log frequency, dB) |
| [`SigmaPoint`](inc/damp/analysis/frequency.hpp#L739) | block | Frequency-domain analysis (host) | Singular values of G(jω) at one frequency |
| [`SigmaResult`](inc/damp/analysis/frequency.hpp#L753) | block | Frequency-domain analysis (host) | Singular-value frequency response over a grid |
| [`SignalSource`](inc/damp/concepts.hpp#L127) | concept | Core, configuration & backend vocabulary | Self-clocked signal source: u = step(), finished when done() |
| [`SignalStatus`](inc/damp/toolbox/conditioning.hpp#L368) | enum | Embedded helpers (controls-adjacent utilities) | Classification of an analog input against its valid/fault bands |
| [`simc`](inc/damp/design/pid_design.hpp#L341) | function | Design-time synthesis (not PWM-rate) | SIMC (Skogestad Internal Model Control) tuning for FOPDT models |
| [`simulate`](inc/damp/simulation/simulate.hpp#L132) | function | Simulation / SIL harness (host) | Simulate a nonlinear plant with a controller in closed loop |
| [`simulate_cached_zoh`](inc/damp/simulation/cached_zoh.hpp#L148) | function | Simulation / SIL harness (host) | Simulate a continuous LTI plant with fixed ZOH step and held input (+1 more overload) |
| [`simulate_cached_zoh_multirate`](inc/damp/simulation/cached_zoh.hpp#L227) | function | Simulation / SIL harness (host) | Multi-rate ZOH: held control samples, dense plant samples for plotting |
| [`simulate_discrete`](inc/damp/simulation/simulate.hpp#L478) | function | Simulation / SIL harness (host) | Simulate a discrete-time system with a controller |
| [`simulate_discrete_nonlinear`](inc/damp/simulation/simulate.hpp#L429) | function | Simulation / SIL harness (host) | Simulate a discrete-time nonlinear plant with a controller |
| [`simulate_lti`](inc/damp/simulation/simulate.hpp#L319) | function | Simulation / SIL harness (host) | Simulate a continuous LTI system with a controller |
| [`simulate_lti_siso`](inc/damp/simulation/simulate.hpp#L377) | function | Simulation / SIL harness (host) | SIL: continuous SISO LTI plant under a sampled SISO controller with constant reference |
| [`simulate_multirate`](inc/damp/simulation/multirate.hpp#L128) | function | Simulation / SIL harness (host) | Multi-rate closed-loop simulation with MultiRateConfig (+1 more overload) |
| [`simulate_piecewise_lti`](inc/damp/simulation/hybrid.hpp#L182) | function | Simulation / SIL harness (host) | Simulate a piecewise-LTI plant with Exact integration between events |
| [`simulate_sampled`](inc/damp/simulation/simulate.hpp#L205) | function | Simulation / SIL harness (host) | Simulate a continuous plant under a discrete (sampled) controller — multi-rate |
| [`simulate_state_feedback`](inc/damp/simulation/simulate.hpp#L264) | function | Simulation / SIL harness (host) | Simulate a nonlinear plant with state-feedback controller |
| [`simulate_two_rate`](inc/damp/simulation/multirate.hpp#L244) | function | Simulation / SIL harness (host) | Two-rate cascade with MultiRateConfig (+1 more overload) |
| [`simulate_two_rate_lti`](inc/damp/simulation/multirate.hpp#L388) | function | Simulation / SIL harness (host) | Two-rate cascade on a continuous LTI plant (A, B, C) |
| [`SimulationResult`](inc/damp/simulation/simulate.hpp#L103) | block | Simulation / SIL harness (host) | Result of a closed-loop simulation |
| [`sin`](inc/damp/matrix/functions.hpp#L723) | function | Linear algebra | Matrix sine via scaling and double-angle reconstruction |
| [`sin`](inc/damp/math/math.hpp#L182) | function | Scalar math, complex & frames | Sine |
| [`sincos`](inc/damp/matrix/functions.hpp#L653) | function | Linear algebra | Compute sin(A) and cos(A) together via scaling and double-angle reconstruction |
| [`sincos`](inc/damp/math/math.hpp#L200) | function | Scalar math, complex & frames | Combined sine and cosine, {sin(x), cos(x)} |
| [`SinglePhasePLL`](inc/damp/filters/pll.hpp#L43) | block | Filters & signal conditioning | Single-Phase PLL |
| [`sinh`](inc/damp/matrix/functions.hpp#L793) | function | Linear algebra | Matrix hyperbolic sine sinh(A) = (exp(A) − exp(−A))/2 |
| [`siso_ref`](inc/damp/simulation/simulate.hpp#L90) | function | Simulation / SIL harness (host) | Build a SisoReferenceAdapter for a SISO `control(r,y)` controller |
| [`SisoController`](inc/damp/concepts.hpp#L74) | concept | Core, configuration & backend vocabulary | Scalar output-feedback controller: u = control(r, y), fixed rate |
| [`SisoReferenceAdapter`](inc/damp/simulation/simulate.hpp#L71) | block | Simulation / SIL harness (host) | Adapt a SISO `control(r, y)` controller for callables that want `u = f(y)` or `u = f(t, y)` (the older simulate contracts) |
| [`six_step_duty_cycles`](inc/damp/motor/modulation.hpp#L462) | function | Motor control pack (if present) | Classical six-step (full-wave) duty pattern from the αβ angle |
| [`six_step_fundamental_voltage`](inc/damp/motor/modulation.hpp#L146) | function | Motor control pack (if present) | Six-step (square-wave) fundamental peak phase voltage |
| [`SlewLimiter`](inc/damp/toolbox/conditioning.hpp#L268) | block | Embedded helpers (controls-adjacent utilities) | Slew-rate limiter: bound how fast the output may follow the target |
| [`smc`](inc/damp/controllers/smc.hpp#L79) | function | Design-time synthesis (not PWM-rate) | Bundle hand-picked SMC parameters into an SMCResult |
| [`SMCController`](inc/damp/controllers/smc.hpp#L156) | block | Runtime controllers | First-order sliding-mode controller (SMC) for a SISO plant |
| [`SMCResult`](inc/damp/controllers/smc.hpp#L33) | block | Design-time synthesis (not PWM-rate) | Tuning parameters for a first-order sliding-mode controller |
| [`smith_predictor_from_fopdt`](inc/damp/controllers/smith_predictor.hpp#L181) | function | Design-time synthesis (not PWM-rate) | Design a Smith predictor from FOPDT parameters |
| [`smith_predictor_from_pid`](inc/damp/controllers/smith_predictor.hpp#L247) | function | Design-time synthesis (not PWM-rate) | Pack a Smith predictor from an existing continuous PID and FOPDT plant |
| [`SmithPredictor`](inc/damp/controllers/smith_predictor.hpp#L307) | block | Runtime controllers | Smith predictor runtime: SISO primary + FO model + pure delay |
| [`SmithPredictorResult`](inc/damp/controllers/smith_predictor.hpp#L85) | block | Design-time synthesis (not PWM-rate) | Smith predictor design result (discrete PID + FO model + delay samples) |
| [`soft_sign`](inc/damp/toolbox/conditioning.hpp#L103) | function | Embedded helpers (controls-adjacent utilities) | Smooth unit direction x / √(x²+ε²) ∈ (−1,1) |
| [`sogi`](inc/damp/filters/sogi.hpp#L39) | function | Filters & signal conditioning | Second-Order Generalized Integrator (SOGI) design (+1 more overload) |
| [`SOGI`](inc/damp/filters/sogi.hpp#L177) | block | Filters & signal conditioning | Runtime SOGI wrapper around design::sogi(w0, alpha, Ts) |
| [`SogiFll`](inc/damp/filters/sogi.hpp#L330) | block | Filters & signal conditioning | SOGI with a Frequency-Locked Loop — self-tuning single-tone tracker |
| [`solve`](inc/damp/matrix/solve.hpp#L82) | function | Linear algebra | Solve lower-triangular system LX = B via forward substitution (+2 more overloads) |
| [`solve_miqp`](inc/damp/design/qp.hpp#L1215) | function | Design-time synthesis (not PWM-rate) | Mixed-integer QP by branch and bound over the active-set relaxation |
| [`solve_qp`](inc/damp/design/qp.hpp#L525) | function | Design-time synthesis (not PWM-rate) | Solve a strictly convex inequality-constrained QP (dual active-set) (+1 more overload) |
| [`solve_qp_admm`](inc/damp/design/qp.hpp#L916) | function | Design-time synthesis (not PWM-rate) | One-shot ADMM QP solve (cold start) |
| [`solve_qp_interior_point`](inc/damp/design/qp.hpp#L960) | function | Design-time synthesis (not PWM-rate) | Interior-point QP solver (primal-dual path following) |
| [`SolveResult`](inc/damp/simulation/solver.hpp#L59) | block | Simulation / SIL harness (host) | Result of an ODE solve operation |
| [`specific_force_at_rest`](inc/damp/estimation/ins_mechanization.hpp#L212) | function | Observers & estimators | Specific force [m/s²] a stationary IMU measures for the given orientation |
| [`SplineSurface`](inc/damp/toolbox/lookup.hpp#L393) | block | Embedded helpers (controls-adjacent utilities) | 2-D interpolating surface with smooth (Catmull-Rom spline) blending |
| [`sqrt`](inc/damp/matrix/functions.hpp#L446) | function | Linear algebra | Matrix square root via Denman–Beavers iteration |
| [`sqrt`](inc/damp/math/complex.hpp#L331) | function | Scalar math, complex & frames | Compute complex square root (constexpr) (+1 more overload) |
| [`sqrtm`](inc/damp/matrix/functions.hpp#L483) | function | Linear algebra | @brief MATLAB®-style alias for sqrt |
| [`SR`](inc/damp/toolbox/iec61131.hpp#L63) | block | Embedded helpers (controls-adjacent utilities) | SR Latch (Set-dominant Set-Reset Latch) |
| [`ss`](inc/damp/matlab.hpp#L98) | function | MATLAB®-style aliases (host) | MATLAB®-style state-space model constructor |
| [`ss2tf`](inc/damp/systems/zpk.hpp#L843) | function | LTI systems (SS / TF / ZPK / discretize) | MATLAB®-style alias for SISO to_transfer_function |
| [`ss2zpk`](inc/damp/systems/zpk.hpp#L813) | function | LTI systems (SS / TF / ZPK / discretize) | MATLAB®-style alias for SISO to_zpk (state-space) |
| [`ss_gain`](inc/damp/systems/state_space.hpp#L510) | function | LTI systems (SS / TF / ZPK / discretize) | Static gain $`y = k u`$ (no states) |
| [`ss_integrator`](inc/damp/systems/state_space.hpp#L520) | function | LTI systems (SS / TF / ZPK / discretize) | Integrator $`y = (k/s)\, u`$ |
| [`SsFeed`](inc/damp/systems/state_space.hpp#L499) | block | LTI systems (SS / TF / ZPK / discretize) | One term of a linear output mix: $`k \cdot y_{\mathrm{iy}}`$ |
| [`stability_margin_continuous`](inc/damp/design/stability.hpp#L268) | function | Design-time synthesis (not PWM-rate) | Compute stability margin for continuous system |
| [`stability_margin_discrete`](inc/damp/design/stability.hpp#L297) | function | Design-time synthesis (not PWM-rate) | Compute stability margin for discrete system |
| [`state_mpc`](inc/damp/controllers/mpc.hpp#L507) | function | Design-time synthesis (not PWM-rate) | Synthesize a constrained linear MPC (condensed dense QP, Δu form) (+1 more overload) |
| [`StateEstimator`](inc/damp/concepts.hpp#L114) | concept | Core, configuration & backend vocabulary | State estimator: x̂ = estimate(y, u) — the fused per-tick form |
| [`StateFeedback`](inc/damp/controllers/lqr.hpp#L454) | block | Runtime controllers | Runtime full-state feedback law u = −Kx |
| [`StateFeedbackController`](inc/damp/concepts.hpp#L101) | concept | Core, configuration & backend vocabulary | State-feedback law: u = control(r, x), reference-first |
| [`StateJacobian`](inc/damp/estimation/ekf.hpp#L41) | block | Observers & estimators | State prediction result from the user's dynamics function |
| [`StateSpace`](inc/damp/systems/state_space.hpp#L128) | block | LTI systems (SS / TF / ZPK / discretize) | State-space representation for linear time-invariant systems (discrete or continuous) |
| [`StateSpaceZPKResult`](inc/damp/systems/zpk.hpp#L570) | block | LTI systems (SS / TF / ZPK / discretize) | SISO state-space → ZPK conversion result with runtime zero count |
| [`steady_state_vdq`](inc/damp/motor/spm.hpp#L79) | function | Motor control pack (if present) | Steady-state Vdq (rotor frame) |
| [`steady_state_voltage_magnitude`](inc/damp/motor/spm.hpp#L92) | function | Motor control pack (if present) | \|Vdq\| at a steady-state operating point |
| [`SteadyStateKalmanFilter`](inc/damp/estimation/kalman.hpp#L413) | block | Observers & estimators | Steady-state (fixed-gain) Kalman estimator for LQG-class designs |
| [`steer_map`](inc/damp/toolbox/io.hpp#L452) | function | Embedded helpers (controls-adjacent utilities) | Curve-driven steering map — arbitrary `(throttle, turn)` geometry |
| [`steinhart_hart`](inc/damp/toolbox/thermistor.hpp#L118) | function | Embedded helpers (controls-adjacent utilities) | Fit the Steinhart-Hart coefficients from three calibration points |
| [`step`](inc/damp/analysis/time_response.hpp#L153) | function | Frequency-domain analysis (host) | Step response of a (MIMO) state-space system (+1 more overload) |
| [`step`](inc/damp/estimation/excitation/step.hpp#L91) | function | Observers & estimators | Build a single-step design payload |
| [`Step`](inc/damp/estimation/excitation/step.hpp#L108) | block | Observers & estimators | Single step with optional \|dy/dt\| settle detection |
| [`step_train`](inc/damp/estimation/excitation/step_train.hpp#L94) | function | Observers & estimators | Build a step-train design payload from a configuration |
| [`StepConfig`](inc/damp/estimation/excitation/step.hpp#L35) | block | Observers & estimators | Configuration for a single step with optional settle detection |
| [`StepInfo`](inc/damp/analysis/time_response.hpp#L354) | block | Frequency-domain analysis (host) | Step-response characteristics of a single output signal |
| [`stepinfo`](inc/damp/analysis/time_response.hpp#L374) | function | Frequency-domain analysis (host) | Compute step-response characteristics from an output/time signal (+1 more overload) |
| [`stepped_sine`](inc/damp/estimation/excitation/stepped_sine.hpp#L104) | function | Observers & estimators | Build a stepped-sine design payload from a configuration |
| [`SteppedSine`](inc/damp/estimation/excitation/stepped_sine.hpp#L125) | block | Observers & estimators | Stepped-sine excitation — one pure tone at a time across a frequency table |
| [`SteppedSineConfig`](inc/damp/estimation/excitation/stepped_sine.hpp#L33) | block | Observers & estimators | Configuration for single-tone stepped-sine excitation |
| [`SteppedSineResult`](inc/damp/estimation/excitation/stepped_sine.hpp#L73) | block | Observers & estimators | Stepped-sine design payload |
| [`stepplot`](inc/damp/simulation/plot_plotly.hpp#L717) | function | Simulation / SIL harness (host) | Plot a step response, one trace per input/output pair |
| [`StepResult`](inc/damp/estimation/excitation/step.hpp#L67) | block | Observers & estimators | Single-step design payload |
| [`StepTrain`](inc/damp/estimation/excitation/step_train.hpp#L110) | block | Observers & estimators | Alternating +/- step train runtime generator |
| [`StepTrainConfig`](inc/damp/estimation/excitation/step_train.hpp#L32) | block | Observers & estimators | Configuration for alternating +/- step excitation |
| [`StepTrainResult`](inc/damp/estimation/excitation/step_train.hpp#L64) | block | Observers & estimators | Step-train design payload |
| [`Stopwatch`](inc/damp/toolbox/timing.hpp#L43) | block | Embedded helpers (controls-adjacent utilities) | Free-running elapsed-time accumulator |
| [`stsmc`](inc/damp/controllers/stsmc.hpp#L87) | function | Design-time synthesis (not PWM-rate) | Synthesize super-twisting gains from a disturbance-derivative bound |
| [`stsmc_gains`](inc/damp/controllers/stsmc.hpp#L123) | function | Design-time synthesis (not PWM-rate) | Super-twisting controller from gains you specify directly |
| [`STSMCController`](inc/damp/controllers/stsmc.hpp#L179) | block | Runtime controllers | Super-twisting controller (second-order sliding mode) |
| [`STSMCResult`](inc/damp/controllers/stsmc.hpp#L30) | block | Design-time synthesis (not PWM-rate) | Super-twisting (second-order sliding-mode) controller design result |
| [`subtract`](inc/damp/systems/state_space.hpp#L397) | function | LTI systems (SS / TF / ZPK / discretize) | Differencing connection: outputs y = y₁ − y₂ |
| [`SuccessiveCompensatorCommissioner`](inc/damp/estimation/successive_compensator.hpp#L208) | block | Observers & estimators | Successive biquad compensator bank (notch and/or band-pass) |
| [`SuccessiveCompensatorConfig`](inc/damp/estimation/successive_compensator.hpp#L50) | block | Observers & estimators | Configuration for successive compensator commissioning |
| [`suggest_mpc_horizon`](inc/damp/controllers/mpc.hpp#L590) | function | Design-time synthesis (not PWM-rate) | Suggest MPC horizons from an explicit settling-time target (+1 more overload) |
| [`summarize_loop_response`](inc/damp/analysis/frequency.hpp#L460) | function | Frequency-domain analysis (host) | Summarize loop_response() results into one compact metrics struct |
| [`svd`](inc/damp/matrix/svd.hpp#L246) | function | Linear algebra | Full singular value decomposition A = U·Σ·Vᴴ (one-sided Jacobi) |
| [`svd`](inc/damp/matlab.hpp#L444) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the singular value decomposition |
| [`SVDResult`](inc/damp/matrix/svd.hpp#L228) | block | Linear algebra | Result of a full singular value decomposition A = U·Σ·Vᴴ |
| [`svm_duty_cycles`](inc/damp/motor/modulation.hpp#L439) | function | Motor control pack (if present) | Space-vector PWM duty cycles from an αβ voltage command |
| [`SvmDuties`](inc/damp/motor/modulation.hpp#L103) | block | Motor control pack (if present) | Result of a duty-map: half-bridge duties plus an over-modulation flag |
| [`svpwm_zero_sequence`](inc/damp/motor/modulation.hpp#L200) | function | Motor control pack (if present) | Min-max zero-sequence injection for space-vector PWM |
| [`Switch`](inc/damp/toolbox/io.hpp#L282) | block | Embedded helpers (controls-adjacent utilities) | Debounced maintained switch (toggle/selector contact) with change flag |
| [`SwitchedController`](inc/damp/controllers/composition.hpp#L117) | block | Runtime controllers | Bumpless-ish switch between normal, experiment, and backup SISO laws |
| [`SwitchMode`](inc/damp/controllers/composition.hpp#L94) | enum | Runtime controllers | Which path owns the plant command |
| [`symmetrical_components`](inc/damp/math/transforms.hpp#L520) | function | Scalar math, complex & frames | Forward symmetrical-component (Fortescue) transform (abc → 012) |
| [`SymplecticEuler`](inc/damp/simulation/integrator.hpp#L229) | block | Simulation / SIL harness (host) | Semi-implicit (symplectic) Euler for mechanical systems |
| [`Tachometer`](inc/damp/toolbox/encoder.hpp#L144) | block | Embedded helpers (controls-adjacent utilities) | Pulse-based speed (tachometer) with frequency/period crossover |
| [`tan`](inc/damp/math/math.hpp#L215) | function | Scalar math, complex & frames | Tangent |
| [`tf`](inc/damp/matlab.hpp#L62) | function | MATLAB®-style aliases (host) | MATLAB®-style transfer function constructor (+1 more overload) |
| [`tf2zpk`](inc/damp/systems/zpk.hpp#L471) | function | LTI systems (SS / TF / ZPK / discretize) | MATLAB®-style alias for to_zpk (transfer function) |
| [`TFF`](inc/damp/toolbox/iec61131.hpp#L517) | block | Embedded helpers (controls-adjacent utilities) | T Flip-Flop (toggle on rising edge) |
| [`thermal_kalman_plant`](inc/damp/toolbox/thermal.hpp#L195) | function | Embedded helpers (controls-adjacent utilities) | Prepare a discretized thermal plant for design::kalman |
| [`ThermalKalmanObserver`](inc/damp/toolbox/thermal.hpp#L248) | block | Embedded helpers (controls-adjacent utilities) | Runtime thermal-network Kalman observer |
| [`Thermistor`](inc/damp/toolbox/thermistor.hpp#L166) | block | Embedded helpers (controls-adjacent utilities) | NTC thermistor linearization (resistance → temperature) |
| [`ThermistorCoeffs`](inc/damp/toolbox/thermistor.hpp#L41) | block | Embedded helpers (controls-adjacent utilities) | Fitted NTC coefficients in Steinhart-Hart form |
| [`thipwm_zero_sequence`](inc/damp/motor/modulation.hpp#L226) | function | Motor control pack (if present) | Third-harmonic injection (THIPWM) zero-sequence — 1/6 of fundamental |
| [`ThreePhasePLL`](inc/damp/filters/pll.hpp#L141) | block | Filters & signal conditioning | Synchronous-reference-frame (SRF) PLL for balanced three-phase input |
| [`Timeout`](inc/damp/toolbox/timing.hpp#L78) | block | Embedded helpers (controls-adjacent utilities) | One-shot timeout |
| [`TimeResponse`](inc/damp/analysis/time_response.hpp#L51) | block | Frequency-domain analysis (host) | Multi-channel time-domain response sampled on a time grid |
| [`to_coeffs`](inc/damp/filters/iir_design.hpp#L617) | function | Filters & signal conditioning | Convert StateSpace system to first-order DSP coefficients (+3 more overloads) |
| [`to_double_vector`](inc/damp/simulation/plot_plotly.hpp#L236) | function | Simulation / SIL harness (host) | Convert std::vector\<T\> to std::vector<double> |
| [`to_std_vector`](inc/damp/simulation/plot_plotly.hpp#L211) | function | Simulation / SIL harness (host) | Convert a ColVec<N,T> to std::vector<double> for plotlypp |
| [`to_transfer_function`](inc/damp/systems/zpk.hpp#L826) | function | LTI systems (SS / TF / ZPK / discretize) | SISO state-space → transfer function (Leverrier), size NX+1 / NX+1 |
| [`to_zpk`](inc/damp/systems/zpk.hpp#L429) | function | LTI systems (SS / TF / ZPK / discretize) | Convert a SISO transfer function to zero-pole-gain form (+1 more overload) |
| [`TOF`](inc/damp/toolbox/iec61131.hpp#L226) | block | Embedded helpers (controls-adjacent utilities) | TOF Timer (Timer Off Delay) |
| [`Toggle`](inc/damp/toolbox/logic.hpp#L245) | block | Embedded helpers (controls-adjacent utilities) | Toggle (T flip-flop): output flips on each rising edge of in |
| [`TON`](inc/damp/toolbox/iec61131.hpp#L183) | block | Embedded helpers (controls-adjacent utilities) | TON Timer (Timer On Delay) |
| [`Tone`](inc/damp/estimation/excitation/multi_sine.hpp#L29) | block | Observers & estimators | One sinusoidal component in a multi-sine excitation |
| [`torque`](inc/damp/motor/spm.hpp#L103) | function | Motor control pack (if present) | Electromagnetic torque of an SPM (i_d = 0): T_e = Kₜ i_q |
| [`torque_constant_from_flux`](inc/damp/motor/spm.hpp#L29) | function | Motor control pack (if present) | Kₜ = 1½ p λ [Nm/A] (amplitude-invariant) |
| [`torque_constant_from_Kv`](inc/damp/motor/spm.hpp#L46) | function | Motor control pack (if present) | Kₜ from hobby Kᵥ [RPM/V] (peak line-to-line / bus-volt sense) |
| [`torque_speed_envelope`](inc/damp/motor/spm.hpp#L203) | function | Motor control pack (if present) | Fixed table of SPM torque–speed envelope samples (host / LUT bake) |
| [`TorqueSpeedPoint`](inc/damp/motor/spm.hpp#L144) | block | Motor control pack (if present) | One sample of the SPM torque–speed envelope (MTPA i_d=0 branch) |
| [`TP`](inc/damp/toolbox/iec61131.hpp#L271) | block | Embedded helpers (controls-adjacent utilities) | TP Timer (Timer Pulse) |
| [`TrajectoryBoundary`](inc/damp/trajectory/trajectory_types.hpp#L178) | block | Trajectory value types | Boundary conditions at one endpoint of a polynomial trajectory: a position and its time derivatives through jerk |
| [`TrajectoryLimits`](inc/damp/trajectory/trajectory_types.hpp#L45) | block | Trajectory value types | Asymmetric kinematic limits for a trapezoidal or S-curve motion profile |
| [`TrajectoryState`](inc/damp/trajectory/trajectory_types.hpp#L76) | block | Trajectory value types | A point on a motion profile: commanded position, velocity, acceleration |
| [`TransferFunction`](inc/damp/systems/transfer_function.hpp#L74) | block | LTI systems (SS / TF / ZPK / discretize) | SISO polynomial transfer function G(s) = num(s)/den(s) |
| [`Transform4`](inc/damp/math/geometry.hpp#L826) | block | Scalar math, complex & frames | 4×4 homogeneous transform (SE(3) interop / DH export) |
| [`Translation3`](inc/damp/kinematics/pose.hpp#L37) | block | Kinematics / pose | A 3-D translation — a thin Vec3 with domain-named conveniences |
| [`TransposeView`](inc/damp/matrix/views.hpp#L384) | block | Linear algebra | Non-owning transpose view of a matrix (zero-copy) |
| [`Trapezoidal`](inc/damp/simulation/integrator.hpp#L623) | block | Simulation / SIL harness (host) | Trapezoidal (Tustin) integrator |
| [`TRBDF2`](inc/damp/simulation/integrator.hpp#L497) | block | Simulation / SIL harness (host) | TR-BDF2 composite integrator — the stiff adaptive pair (ode23tb) |
| [`two_norm`](inc/damp/matrix/functions.hpp#L128) | function | Linear algebra | Spectral norm ‖A‖₂ = σₘₐₓ(A) |
| [`two_point_cal`](inc/damp/toolbox/scaling.hpp#L111) | function | Embedded helpers (controls-adjacent utilities) | Fit an AffineCal through two `(raw, engineering)` points |
| [`TwoRateSimulationResult`](inc/damp/simulation/multirate.hpp#L94) | block | Simulation / SIL harness (host) | Result of a two-rate cascade simulation |
| [`tyreus_luyben`](inc/damp/design/pid_design.hpp#L178) | function | Design-time synthesis (not PWM-rate) | Tyreus-Luyben tuning from ultimate gain and ultimate period |
| [`UKFMeasFn`](inc/damp/estimation/ukf.hpp#L66) | concept | Observers & estimators | Concept for UKF measurement functions |
| [`UKFStateFn`](inc/damp/estimation/ukf.hpp#L55) | concept | Observers & estimators | Concept for UKF state (process) functions |
| [`unbounded_bound`](inc/damp/design/qp.hpp#L67) | function | Design-time synthesis (not PWM-rate) | Sentinel bound treated as "no constraint" on that row |
| [`UnscentedKalmanFilter`](inc/damp/estimation/ukf.hpp#L129) | block | Observers & estimators | Unscented (sigma-point) Kalman Filter for nonlinear discrete-time systems |
| [`UnscentedParams`](inc/damp/estimation/ukf.hpp#L82) | block | Observers & estimators | Tuning parameters for the scaled unscented transform |
| [`unwrap_phase_deg`](inc/damp/analysis/frequency.hpp#L125) | function | Frequency-domain analysis (host) | Unwrap phase data in degrees to avoid +/-180 discontinuities |
| [`UpperTriangle`](inc/damp/matrix/views.hpp#L122) | block | Linear algebra | Upper triangular view of a square matrix |
| [`voltage_circle_radius`](inc/damp/motor/foc.hpp#L57) | function | Motor control pack (if present) | SVPWM voltage-circle radius Vₘₐₓ = m · Vdc / √3 |
| [`WarmStartActiveSetSolver`](inc/damp/design/qp.hpp#L654) | block | Design-time synthesis (not PWM-rate) | Warm-started active-set solver policy — the damp::MPC default |
| [`with_schroeder_phases`](inc/damp/estimation/excitation/multi_sine.hpp#L160) | function | Observers & estimators | Assign Schroeder low-crest-factor phases to a multi-sine tone table |
| [`wrap`](inc/damp/math/math.hpp#L509) | function | Scalar math, complex & frames | Wrap x into the half-open interval [min, max) (period max − min) |
| [`wrap_pi`](inc/damp/math/math.hpp#L544) | function | Scalar math, complex & frames | Wrap an angle to [−π, π) |
| [`wrap_two_pi`](inc/damp/math/math.hpp#L562) | function | Scalar math, complex & frames | Wrap an angle to [0, 2π) |
| [`wrapped_delta`](inc/damp/toolbox/encoder.hpp#L41) | function | Embedded helpers (controls-adjacent utilities) | Signed difference between two unsigned counter readings, wrap-safe |
| [`write_html`](inc/damp/simulation/plot_plotly.hpp#L145) | function | Simulation / SIL harness (host) | Write a plotlypp figure to HTML with damp's shell (prefer over writeHtml) |
| [`write_hybrid_trace_npy`](inc/damp/simulation/npy_export.hpp#L100) | function | Simulation / SIL harness (host) | Pack HybridSimulationResult into (N, 2+NX+NU+NY) float64: t, mode, x..., u..., y |
| [`write_npy_f64`](inc/damp/simulation/npy_export.hpp#L40) | function | Simulation / SIL harness (host) | Write a C-order float64 array as a .npy file |
| [`XyPlotOpts`](inc/damp/simulation/plot_plotly.hpp#L535) | block | Simulation / SIL harness (host) | Optional axis framing for plot_xy |
| [`XySeries`](inc/damp/simulation/plot_plotly.hpp#L517) | block | Simulation / SIL harness (host) | One series for a planar (x, y) scatter / path plot |
| [`zero_sequence`](inc/damp/math/transforms.hpp#L226) | function | Scalar math, complex & frames | Zero-sequence (common-mode) scalar of a three-phase set |
| [`ziegler_nichols`](inc/damp/design/pid_design.hpp#L78) | function | Design-time synthesis (not PWM-rate) | Ziegler-Nichols tuning from ultimate gain and ultimate period |
| [`ziegler_nichols_step`](inc/damp/design/pid_design.hpp#L127) | function | Design-time synthesis (not PWM-rate) | Ziegler-Nichols step response method (reaction curve) |
| [`ZPK`](inc/damp/systems/zpk.hpp#L177) | block | LTI systems (SS / TF / ZPK / discretize) | Zero-pole-gain (ZPK) representation of a SISO LTI system |
| [`zpk`](inc/damp/matlab.hpp#L121) | function | MATLAB®-style aliases (host) | MATLAB®-style zero-pole-gain model constructor |
| [`zpk2tf`](inc/damp/systems/zpk.hpp#L480) | function | LTI systems (SS / TF / ZPK / discretize) | MATLAB®-style alias for ZPK::to_transfer_function |
| [`ZPKResult`](inc/damp/systems/zpk.hpp#L399) | block | LTI systems (SS / TF / ZPK / discretize) | Result of converting a transfer function to zero-pole-gain form |
