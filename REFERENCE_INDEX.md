# API Reference — Alphabetical Index

Auto-generated from `@brief` doc comments in `inc/damp/`. Regenerate with `python tools/gen_reference.py`. Grouped-by-domain view: [REFERENCE.md](REFERENCE.md).

| Name | Kind | Domain | Description |
| ---- | ---- | ------ | ----------- |
| [`abs`](inc/damp/math/complex.hpp#L362) | function | Scalar math, complex & frames | Compute magnitude (absolute value) of a complex number (+1 more overload) |
| [`acker`](inc/damp/matlab.hpp#L562) | function | MATLAB®-style aliases (host) | Pole placement for state-feedback control |
| [`ackermann`](inc/damp/design/pole_placement.hpp#L1596) | function | Design-time synthesis (not PWM-rate) | Single-input pole placement via Ackermann's formula |
| [`acos`](inc/damp/math/math.hpp#L150) | function | Scalar math, complex & frames | Arccosine ∈ [0, π]. Input is clamped to [−1, 1] in both paths |
| [`ActionGovernor`](inc/damp/controllers/action_governor.hpp#L331) | block | Runtime controllers | Runtime action governor — QP projection of u_des onto A u ≤ b |
| [`ActiveSetSolver`](inc/damp/design/qp.hpp#L639) | block | Design-time synthesis (not PWM-rate) | Default QP solver policy: the Goldfarb–Idnani active-set solve_qp() |
| [`AdaptationGate`](inc/damp/estimation/parameter_estimation.hpp#L319) | block | Observers & estimators | Gated-adaptation policy: decide whether a parameter update is safe to apply |
| [`adaptive_solve`](inc/damp/simulation/solver.hpp#L802) | function | Simulation / SIL harness (host) | Adaptive-step solve in one call |
| [`AdaptiveLut1D`](inc/damp/toolbox/lookup.hpp#L428) | block | Embedded helpers (controls-adjacent utilities) | Adaptive 1-D lookup: fixed breakpoints, online-updated cell values |
| [`AdaptiveLut2D`](inc/damp/toolbox/lookup.hpp#L501) | block | Embedded helpers (controls-adjacent utilities) | Adaptive 2-D lookup: fixed grid, online-updated cells (nearest grid node) |
| [`AdaptiveOptions`](inc/damp/simulation/solver.hpp#L333) | block | Simulation / SIL harness (host) | Options for adaptive-step ODE integration |
| [`AdaptiveStepIntegrator`](inc/damp/simulation/integrator.hpp#L91) | concept | Simulation / SIL harness (host) | Integrator that reports a genuine embedded local-error estimate |
| [`AdaptiveStepSolver`](inc/damp/simulation/solver.hpp#L388) | block | Simulation / SIL harness (host) | Adaptive-step ODE solver (+1 more overload) |
| [`AdaptOutOfRange`](inc/damp/toolbox/lookup.hpp#L48) | enum | Embedded helpers (controls-adjacent utilities) | How AdaptiveLut1D / AdaptiveLut2D treat samples outside the breakpoint span |
| [`AdmmSettings`](inc/damp/design/qp.hpp#L734) | block | Design-time synthesis (not PWM-rate) | Tuning parameters for the ADMM QP solver |
| [`AdmmSolver`](inc/damp/design/qp.hpp#L773) | block | Design-time synthesis (not PWM-rate) | ADMM (OSQP-style) QP solver policy — warm-started, fixed-cost iterations |
| [`adrc`](inc/damp/controllers/adrc.hpp#L157) | function | Design-time synthesis (not PWM-rate) | Active Disturbance Rejection Control design |
| [`ADRCController`](inc/damp/controllers/adrc.hpp#L219) | block | Runtime controllers | Active Disturbance Rejection Control (ADRC) |
| [`ADRCResult`](inc/damp/controllers/adrc.hpp#L36) | block | Design-time synthesis (not PWM-rate) | Active Disturbance Rejection Control design result |
| [`AffineCal`](inc/damp/toolbox/scaling.hpp#L83) | block | Embedded helpers (controls-adjacent utilities) | Affine sensor calibration `y = gain·x + offset` |
| [`allmargin`](inc/damp/matlab.hpp#L952) | function | MATLAB®-style aliases (host) | Gain, phase, and delay margins of a SISO loop over a frequency grid |
| [`AllMarginResult`](inc/damp/matlab.hpp#L931) | block | MATLAB®-style aliases (host) | All classical margins including delay margin (superset of margin) |
| [`allpass`](inc/damp/filters/iir_design.hpp#L901) | function | Filters & signal conditioning | Second-order all-pass filter (Tustin of the analog prototype) |
| [`allpass_2nd_continuous`](inc/damp/filters/iir_design.hpp#L462) | function | Filters & signal conditioning | Second-order all-pass filter design (continuous-time) |
| [`AlphaBeta`](inc/damp/math/transforms.hpp#L152) | block | Scalar math, complex & frames | Alpha-beta (stationary-frame) component pair |
| [`amigo_kappa_tau`](inc/damp/design/pid_design.hpp#L238) | function | Design-time synthesis (not PWM-rate) | AMIGO PI from ultimate gain/period and static gain Kₛ |
| [`AnalogCrossMode`](inc/damp/toolbox/conditioning.hpp#L508) | enum | Embedded helpers (controls-adjacent utilities) | How two redundant analog channels should relate in engineering units |
| [`AnalogCrossResult`](inc/damp/toolbox/conditioning.hpp#L530) | block | Embedded helpers (controls-adjacent utilities) | Output of cross_check_analog |
| [`AnalogCrossStatus`](inc/damp/toolbox/conditioning.hpp#L516) | enum | Embedded helpers (controls-adjacent utilities) | Result of comparing two redundant analog channels |
| [`AnalogInput`](inc/damp/toolbox/io.hpp#L46) | block | Embedded helpers (controls-adjacent utilities) | A single analog input: range/fault check on the raw reading, then affine calibration to engineering units |
| [`AngleUnwrapper`](inc/damp/toolbox/encoder.hpp#L200) | block | Embedded helpers (controls-adjacent utilities) | Wrap-safe accumulator: successive wrapped absolute readings → continuous position |
| [`apply_wet_plotly_shell`](inc/damp/simulation/plot_plotly.hpp#L64) | function | Simulation / SIL harness (host) | Apply damp's full-window HTML shell to plotlypp-generated markup |
| [`arange`](inc/damp/analysis/linspace.hpp#L151) | function | Frequency-domain analysis (host) | Half-open arithmetic range [start, stop) with step step (+2 more overloads) |
| [`arcade_drive`](inc/damp/toolbox/io.hpp#L361) | function | Embedded helpers (controls-adjacent utilities) | Arcade (single-stick) drive mixer: x (steer) + y (throttle) → left/right |
| [`arg`](inc/damp/math/complex.hpp#L373) | function | Scalar math, complex & frames | Compute argument (phase angle) of a complex number |
| [`array_as`](inc/damp/backend.hpp#L195) | function | Core, configuration & backend vocabulary | Element-wise static_cast of a damp::array to scalar type U (precision rebind) |
| [`asin`](inc/damp/math/math.hpp#L131) | function | Scalar math, complex & frames | Arcsine ∈ [−π/2, π/2]. Input is clamped to [−1, 1] in both paths so behavior matches at compile and run time (std::asin would return NaN for \|x\| > 1) |
| [`atan`](inc/damp/math/math.hpp#L117) | function | Scalar math, complex & frames | Single-argument arctangent ∈ (−π/2, π/2) |
| [`atan2`](inc/damp/math/math.hpp#L105) | function | Scalar math, complex & frames | Two-argument arctangent, atan2(y, x) ∈ [−π, π] |
| [`AxisInput`](inc/damp/toolbox/io.hpp#L191) | block | Embedded helpers (controls-adjacent utilities) | Operator-axis conditioning chain (joystick / RC stick → command) |
| [`backward_substitute_transpose`](inc/damp/matrix/solve.hpp#L59) | function | Linear algebra | Backward substitution for Lᵀx = b (real) / Lᴴx = b layout |
| [`BackwardEuler`](inc/damp/simulation/integrator.hpp#L315) | block | Simulation / SIL harness (host) | Backward Euler integrator |
| [`balanced_realization`](inc/damp/design/model_reduction.hpp#L570) | function | Design-time synthesis (not PWM-rate) | Descriptive alias for balreal |
| [`balanced_truncation`](inc/damp/design/model_reduction.hpp#L587) | function | Design-time synthesis (not PWM-rate) | Descriptive alias for balred |
| [`BalancedRealizationResult`](inc/damp/design/model_reduction.hpp#L74) | block | Design-time synthesis (not PWM-rate) | Balanced realization result (Moore square-root method) |
| [`BalancedReductionResult`](inc/damp/design/model_reduction.hpp#L100) | block | Design-time synthesis (not PWM-rate) | Reduced-order model from balanced truncation or residualization |
| [`balreal`](inc/damp/design/model_reduction.hpp#L384) | function | Design-time synthesis (not PWM-rate) | Balanced realization via the square-root (Moore / Laub) method |
| [`balred`](inc/damp/design/model_reduction.hpp#L472) | function | Design-time synthesis (not PWM-rate) | Balanced truncation to NR states |
| [`bandpass`](inc/damp/filters/iir_design.hpp#L852) | function | Filters & signal conditioning | Second-order band-pass filter (constant 0 dB peak gain) |
| [`bandpass_continuous`](inc/damp/filters/iir_design.hpp#L406) | function | Filters & signal conditioning | Second-order band-pass filter design (continuous-time) |
| [`bandwidth`](inc/damp/matlab.hpp#L1057) | function | MATLAB®-style aliases (host) | -3 dB bandwidth of a SISO system over a frequency grid |
| [`bandwidth_from_settling_time`](inc/damp/design/pid_design.hpp#L869) | function | Design-time synthesis (not PWM-rate) | Map settling-time and damping-ratio targets to a bandwidth estimate |
| [`BDF2`](inc/damp/simulation/integrator.hpp#L374) | block | Simulation / SIL harness (host) | Backward Differentiation Formula 2 (BDF2) integrator |
| [`beta`](inc/damp/toolbox/thermistor.hpp#L77) | function | Embedded helpers (controls-adjacent utilities) | Fit NTC coefficients from the Beta-parameter model |
| [`Biquad`](inc/damp/filters/biquad.hpp#L34) | block | Filters & signal conditioning | Second-order IIR (biquad) section runtime |
| [`BiquadCascade`](inc/damp/filters/biquad.hpp#L101) | block | Filters & signal conditioning | Cascade of second-order sections (SOS) for higher-order IIR filters |
| [`BLINK`](inc/damp/toolbox/iec61131.hpp#L542) | block | Embedded helpers (controls-adjacent utilities) | BLINK (free-running square-wave / flasher) |
| [`blkdiag`](inc/damp/matlab.hpp#L361) | function | MATLAB®-style aliases (host) | Block diagonal matrix construction |
| [`Block`](inc/damp/matrix/block.hpp#L38) | block | Linear algebra | Block view (non-owning) into a parent matrix |
| [`bode`](inc/damp/analysis/frequency.hpp#L258) | function | Frequency-domain analysis (host) | Compute Bode plot data for a SISO state-space system (+2 more overloads) |
| [`bode_discrete`](inc/damp/analysis/frequency.hpp#L371) | function | Frequency-domain analysis (host) | Compute Bode plot data for a discrete-time SISO state-space system |
| [`bodemag`](inc/damp/simulation/plot_plotly.hpp#L796) | function | Simulation / SIL harness (host) | Plot a magnitude-only Bode diagram (log frequency, dB magnitude) |
| [`bodeplot`](inc/damp/simulation/plot_plotly.hpp#L782) | function | Simulation / SIL harness (host) | Plot magnitude and phase Bode subplots |
| [`BodeResult`](inc/damp/analysis/frequency.hpp#L60) | block | Frequency-domain analysis (host) | Bode plot data for a SISO system |
| [`Bounds`](inc/damp/toolbox/bounds.hpp#L43) | block | Embedded helpers (controls-adjacent utilities) | A per-channel closed-interval box constraint |
| [`BoxCommandFilter`](inc/damp/controllers/action_governor.hpp#L119) | block | Runtime controllers | Runtime box command filter — closed-form clamp per tick |
| [`build_lqg_analysis_models`](inc/damp/design/synthesis.hpp#L214) | function | Design-time synthesis (not PWM-rate) | Build analysis models from an LQG design |
| [`build_lqgi_analysis_models`](inc/damp/design/synthesis.hpp#L300) | function | Design-time synthesis (not PWM-rate) | Build analysis models from an LQGI servo design |
| [`build_lqi_analysis_models`](inc/damp/design/synthesis.hpp#L269) | function | Design-time synthesis (not PWM-rate) | Build analysis models from an LQI servo design |
| [`build_mpc_analysis_models`](inc/damp/controllers/mpc.hpp#L705) | function | Design-time synthesis (not PWM-rate) | Build the unconstrained-MPC LTI analysis models |
| [`butterworth_lowpass`](inc/damp/filters/iir_design.hpp#L492) | function | Filters & signal conditioning | Butterworth low-pass filter design (+1 more overload) |
| [`Button`](inc/damp/toolbox/io.hpp#L238) | block | Embedded helpers (controls-adjacent utilities) | Debounced momentary push-button with edge and long-press detection |
| [`c2d`](inc/damp/matlab.hpp#L308) | function | MATLAB®-style aliases (host) | MATLAB® interface function c2d to discretize a continuous-time state-space system (+1 more overload) |
| [`CachedZoh`](inc/damp/simulation/cached_zoh.hpp#L62) | block | Simulation / SIL harness (host) | Cached discrete ZOH maps Ad, Bd for step size h |
| [`canonical_phase_margin`](inc/damp/analysis/frequency.hpp#L150) | function | Frequency-domain analysis (host) | Normalize phase margin to (-180, 180] |
| [`care`](inc/damp/design/riccati.hpp#L1285) | function | Design-time synthesis (not PWM-rate) | Solve the Continuous-time Algebraic Riccati Equation (CARE) |
| [`Cascade`](inc/damp/controllers/composition.hpp#L52) | block | Runtime controllers | Series cascade of two SISO controllers: outer → inner reference |
| [`catmull_1d`](inc/damp/toolbox/lookup.hpp#L340) | function | Embedded helpers (controls-adjacent utilities) | 1-D non-uniform Catmull-Rom (cubic Hermite) evaluation |
| [`cauer_thermal_ss_ambient`](inc/damp/toolbox/thermal.hpp#L128) | function | Embedded helpers (controls-adjacent utilities) | Continuous Cauer RC ladder in absolute temperature with ambient input |
| [`cauer_thermal_ss_ambient_mimo`](inc/damp/toolbox/thermal.hpp#L155) | function | Embedded helpers (controls-adjacent utilities) | Continuous Cauer RC ladder, absolute temps, all nodes as outputs |
| [`cbf_relative_degree_1`](inc/damp/controllers/action_governor.hpp#L299) | function | Design-time synthesis (not PWM-rate) | Build a relative-degree-1 CBF inequality row |
| [`CBFConstraint`](inc/damp/controllers/action_governor.hpp#L283) | block | Runtime controllers | One affine row from a relative-degree-1 CBF condition |
| [`cbrt`](inc/damp/math/math.hpp#L92) | function | Scalar math, complex & frames | Cube root (preserves sign for negative x) |
| [`ceil`](inc/damp/math/math.hpp#L330) | function | Scalar math, complex & frames | Ceiling — smallest integer ≥ x |
| [`chirp`](inc/damp/estimation/excitation/chirp.hpp#L115) | function | Observers & estimators | Build a chirp design payload from a configuration |
| [`Chirp`](inc/damp/estimation/excitation/chirp.hpp#L132) | block | Observers & estimators | Linear or logarithmic chirp runtime generator |
| [`ChirpConfig`](inc/damp/estimation/excitation/chirp.hpp#L40) | block | Observers & estimators | Configuration for a sine chirp excitation |
| [`ChirpMode`](inc/damp/estimation/excitation/chirp.hpp#L26) | enum | Observers & estimators | Chirp sweep law |
| [`ChirpResult`](inc/damp/estimation/excitation/chirp.hpp#L83) | block | Observers & estimators | Chirp design payload |
| [`cholesky`](inc/damp/matrix/decomposition.hpp#L74) | function | Linear algebra | Cholesky decomposition for positive-definite matrices |
| [`cholesky_solve`](inc/damp/matrix/solve.hpp#L157) | function | Linear algebra | Solve AX = B via Cholesky (A = LLᴴ) |
| [`clarke_park_transform`](inc/damp/math/transforms.hpp#L380) | function | Scalar math, complex & frames | Fused Clarke-Park transform (abc → dq) |
| [`clarke_transform`](inc/damp/math/transforms.hpp#L253) | function | Scalar math, complex & frames | Clarke transform (abc → αβ) |
| [`classical_dob`](inc/damp/estimation/dob.hpp#L360) | function | Observers & estimators | Synthesize a classical disturbance observer from a nominal plant and Q-filter |
| [`ClassicalDOB`](inc/damp/estimation/dob.hpp#L392) | block | Observers & estimators | Classical Pn^-1·Q disturbance observer runtime (bolt-on compensator) |
| [`ClassicalDobResult`](inc/damp/estimation/dob.hpp#L303) | block | Observers & estimators | Design result for the classical Pn^-1·Q disturbance observer |
| [`classify_range`](inc/damp/toolbox/conditioning.hpp#L406) | function | Embedded helpers (controls-adjacent utilities) | Classify x against the four band edges `[fault_lo (valid_lo, valid_hi) fault_hi]` (assumed ordered, non-decreasing) |
| [`closed_loop_poles`](inc/damp/design/stability.hpp#L378) | function | Design-time synthesis (not PWM-rate) | Compute closed-loop poles (eigenvalues) with state feedback |
| [`cohen_coon`](inc/damp/design/pid_design.hpp#L285) | function | Design-time synthesis (not PWM-rate) | Cohen-Coon tuning from first-order-plus-dead-time model |
| [`ColVec`](inc/damp/matrix/colvec.hpp#L27) | block | Linear algebra | Concrete Column vector specialization of Matrix<N, 1, T> |
| [`ColView`](inc/damp/matrix/views.hpp#L295) | block | Linear algebra | Non-owning column view of a matrix |
| [`comb_notch_window`](inc/damp/filters/moving_average.hpp#L34) | function | Filters & signal conditioning | Window length for a moving-average comb that notches f_notch and all its harmonics: N = round(fs / f_notch) |
| [`CommandProjectionResult`](inc/damp/controllers/action_governor.hpp#L185) | block | Design-time synthesis (not PWM-rate) | Result of an affine command projection |
| [`compensate_coulomb_viscous`](inc/damp/toolbox/conditioning.hpp#L188) | function | Embedded helpers (controls-adjacent utilities) | Coulomb + viscous friction compensation (static inverse, massless) |
| [`compensate_stribeck`](inc/damp/toolbox/conditioning.hpp#L215) | function | Embedded helpers (controls-adjacent utilities) | Stribeck-style friction compensation (static + Coulomb + viscous) |
| [`compensator_from_feature`](inc/damp/estimation/successive_compensator.hpp#L179) | function | Observers & estimators | Map an FRF feature to RBJ biquad coefficients |
| [`CompensatorKind`](inc/damp/estimation/successive_compensator.hpp#L38) | enum | Observers & estimators | Which biquad family to load for a feature |
| [`Complementary`](inc/damp/filters/complementary.hpp#L32) | block | Filters & signal conditioning | Scalar (1-D) complementary filter — fuse a fast rate with a slow absolute |
| [`ComplementaryFilter`](inc/damp/estimation/sensor_fusion.hpp#L56) | block | Observers & estimators | Simple complementary filter for orientation estimation |
| [`complex`](inc/damp/math/complex.hpp#L38) | block | Scalar math, complex & frames | Constexpr complex number class for compile-time computations |
| [`compute_eigenvalues`](inc/damp/matrix/eigen.hpp#L421) | function | Linear algebra | Compute the eigenvalues (and Schur vectors) of a real square matrix |
| [`constant_power_torque_limit`](inc/damp/motor/foc.hpp#L41) | function | Motor control pack (if present) | Constant-power torque ceiling Tₘₐₓ = P_rated / \|ω\| |
| [`ConstantInertiaFeedforward`](inc/damp/toolbox/actuator.hpp#L197) | block | Embedded helpers (controls-adjacent utilities) | Per-axis decoupled torque feedforward: `τ = J·a + b·v + τ_c·sign(v) + g` |
| [`continuous_lpf_exact_step`](inc/damp/filters/lowpass.hpp#L37) | function | Filters & signal conditioning | Exact step of continuous first-order LPF ẏ = −ω_c (y − u) |
| [`continuous_lqr`](inc/damp/controllers/lqr.hpp#L390) | function | Design-time synthesis (not PWM-rate) | Continuous-time Linear-Quadratic Regulator design (+1 more overload) |
| [`ContinuousPID`](inc/damp/controllers/pid.hpp#L798) | block | Runtime controllers | Continuous-gain PID with per-tick sample time (variable-rate secondary form) |
| [`controllability_gramian`](inc/damp/design/stability.hpp#L161) | function | Design-time synthesis (not PWM-rate) | Continuous/discrete controllability Gramian W_c |
| [`controllability_matrix`](inc/damp/design/stability.hpp#L52) | function | Design-time synthesis (not PWM-rate) | Compute the controllability matrix [B, AB, A²B, ..., A^(N-1)B] |
| [`Convention`](inc/damp/math/transforms.hpp#L69) | enum | Scalar math, complex & frames | Scaling convention for the Clarke/Park family |
| [`copysign`](inc/damp/math/math.hpp#L407) | function | Scalar math, complex & frames | Copy sign — magnitude of mag with the sign of sgn_src |
| [`cos`](inc/damp/matrix/functions.hpp#L1062) | function | Linear algebra | Matrix cosine via scaling and double-angle reconstruction |
| [`cos`](inc/damp/math/math.hpp#L169) | function | Scalar math, complex & frames | Cosine |
| [`cosh`](inc/damp/matrix/functions.hpp#L1133) | function | Linear algebra | Matrix hyperbolic cosine cosh(A) = (exp(A) + exp(−A))/2 |
| [`Counter`](inc/damp/toolbox/logic.hpp#L273) | block | Embedded helpers (controls-adjacent utilities) | Edge-counting up/down counter: increments on each rising edge of up, decrements on each rising edge of down. Returns the running count, saturating at the limits of T (no wrap / signed overflow) |
| [`cross_check_analog`](inc/damp/toolbox/conditioning.hpp#L584) | function | Embedded helpers (controls-adjacent utilities) | Cross-check two calibrated analog readings for redundant sensors |
| [`CTD`](inc/damp/toolbox/iec61131.hpp#L354) | block | Embedded helpers (controls-adjacent utilities) | CTD Counter (Count Down) |
| [`ctrb`](inc/damp/matlab.hpp#L270) | function | MATLAB®-style aliases (host) | MATLAB® short alias for controllability_matrix (+1 more overload) |
| [`CTU`](inc/damp/toolbox/iec61131.hpp#L314) | block | Embedded helpers (controls-adjacent utilities) | CTU Counter (Count Up) |
| [`CTUD`](inc/damp/toolbox/iec61131.hpp#L394) | block | Embedded helpers (controls-adjacent utilities) | CTUD Counter (Count Up Down) |
| [`damp`](inc/damp/analysis/poles.hpp#L93) | function | Frequency-domain analysis (host) | Compute natural frequency and damping for each pole |
| [`damping_ratio_from_overshoot_percent`](inc/damp/design/pid_design.hpp#L814) | function | Design-time synthesis (not PWM-rate) | Map percent overshoot target to equivalent damping ratio |
| [`dare`](inc/damp/design/riccati.hpp#L1161) | function | Design-time synthesis (not PWM-rate) | Solve the Discrete Algebraic Riccati Equation (DARE) |
| [`db2mag`](inc/damp/math/math.hpp#L471) | function | Scalar math, complex & frames | Decibels to magnitude, 10^(db/20) |
| [`dcgain`](inc/damp/analysis/norms.hpp#L40) | function | Frequency-domain analysis (host) | Compute DC gain of a continuous-time system |
| [`DCM`](inc/damp/math/geometry.hpp#L63) | block | Scalar math, complex & frames | Direction cosine matrix — 3×3 rotation (SO(3) wrapper over Mat3) |
| [`deadband`](inc/damp/toolbox/conditioning.hpp#L47) | function | Embedded helpers (controls-adjacent utilities) | Dead zone over `[lower, upper]`, matching Simulink®'s Dead Zone block (+1 more overload) |
| [`Debounce`](inc/damp/toolbox/logic.hpp#L212) | block | Embedded helpers (controls-adjacent utilities) | Debounce: the output adopts in only after in differs from the current output continuously for stable_time. Rejects contact bounce and brief glitches. (Not an IEC block — the one everyone hand-rolls.) |
| [`decoupling_feedforward`](inc/damp/motor/foc.hpp#L69) | function | Motor control pack (if present) | Decoupling + back-EMF feedforward at R = 0 (reference currents) |
| [`default_tol`](inc/damp/matrix/matrix_traits.hpp#L94) | function | Linear algebra | Type-appropriate default tolerance for floating-point comparisons. float  ~7 decimal digits  → 1e-6 double ~15 decimal digits → 1e-12 |
| [`deg2rad`](inc/damp/math/math.hpp#L493) | function | Scalar math, complex & frames | Degrees to radians, deg·π/180 |
| [`Delay`](inc/damp/filters/delay.hpp#L27) | block | Filters & signal conditioning | Discrete-time delay buffer |
| [`det`](inc/damp/matrix/functions.hpp#L151) | function | Linear algebra | Matrix determinant det(A) |
| [`DFF`](inc/damp/toolbox/iec61131.hpp#L451) | block | Embedded helpers (controls-adjacent utilities) | D Flip-Flop (edge-triggered data latch) |
| [`diag`](inc/damp/matlab.hpp#L384) | function | MATLAB®-style aliases (host) | Returns a square diagonal matrix from the given array (+1 more overload) |
| [`Diagonal`](inc/damp/matrix/views.hpp#L49) | block | Linear algebra | Diagonal view of a square matrix |
| [`differential_drive`](inc/damp/toolbox/io.hpp#L340) | function | Embedded helpers (controls-adjacent utilities) | Differential (tank) drive mixer: throttle + turn → left/right |
| [`DirectQuadrature`](inc/damp/math/transforms.hpp#L84) | block | Scalar math, complex & frames | Direct-quadrature (rotor-frame) component pair |
| [`Discrete`](inc/damp/simulation/integrator.hpp#L109) | block | Simulation / SIL harness (host) | Discrete-time integrator (no integration, just one step) |
| [`discrete_lqe`](inc/damp/estimation/kalman.hpp#L158) | function | Observers & estimators | Discrete linear-quadratic estimator (LQE) — dual-of-LQR spelling of kalman (+1 more overload) |
| [`discrete_lqg`](inc/damp/controllers/lqg.hpp#L125) | function | Design-time synthesis (not PWM-rate) | Discrete Linear-Quadratic-Gaussian regulator design |
| [`discrete_lqgi`](inc/damp/controllers/lqgi.hpp#L133) | function | Design-time synthesis (not PWM-rate) | Discrete LQG with integral action (LQI + Kalman) for output tracking |
| [`discrete_lqi`](inc/damp/controllers/lqi.hpp#L93) | function | Design-time synthesis (not PWM-rate) | Discrete Linear-Quadratic-Integral (LQI) design for output tracking |
| [`discrete_lqr`](inc/damp/controllers/lqr.hpp#L121) | function | Design-time synthesis (not PWM-rate) | Discrete-time Linear-Quadratic Regulator design |
| [`discrete_lqr_from_continuous`](inc/damp/controllers/lqr.hpp#L259) | function | Design-time synthesis (not PWM-rate) | Design discrete LQR from continuous-time system via discretization (+1 more overload) |
| [`discrete_oflqg`](inc/damp/controllers/offset_free_lqg.hpp#L101) | function | Design-time synthesis (not PWM-rate) | Discrete offset-free LQG (LQR + bias-augmented Kalman) |
| [`DiscretePIDResult`](inc/damp/controllers/pid.hpp#L46) | block | Design-time synthesis (not PWM-rate) | Fixed-rate discrete PID coefficients (canonical deploy form) |
| [`DiscretizationMethod`](inc/damp/systems/discretization.hpp#L26) | enum | LTI systems (SS / TF / ZPK / discretize) | Discretization methods for continuous-time state-space systems |
| [`discretize`](inc/damp/systems/discretization.hpp#L207) | function | LTI systems (SS / TF / ZPK / discretize) | Discretize a continuous-time state-space system |
| [`discretize_lqr_cost`](inc/damp/controllers/lqr.hpp#L187) | function | Design-time synthesis (not PWM-rate) | Discretize a continuous LQR cost integral over one sample (Van Loan) |
| [`DLATCH`](inc/damp/toolbox/iec61131.hpp#L482) | block | Embedded helpers (controls-adjacent utilities) | D Latch (level-sensitive / transparent latch) |
| [`dlqr`](inc/damp/controllers/lqr.hpp#L312) | function | Design-time synthesis (not PWM-rate) | Discrete-time LQR design (MATLAB®-style short name) |
| [`dlqr`](inc/damp/matlab.hpp#L631) | function | MATLAB®-style aliases (host) | Discrete-time Linear-Quadratic Regulator design |
| [`dlyap`](inc/damp/design/lyapunov.hpp#L130) | function | Design-time synthesis (not PWM-rate) | Solve the discrete-time Lyapunov (Stein) equation A X Aᵀ − X + Q = 0 |
| [`dlyap`](inc/damp/matlab.hpp#L1077) | function | MATLAB®-style aliases (host) | MATLAB® alias for the discrete Lyapunov solve AXAᵀ−X+Q=0 |
| [`dob`](inc/damp/estimation/dob.hpp#L136) | function | Observers & estimators | Validate and package DOB configuration into a runtime-ready design result |
| [`DOB`](inc/damp/estimation/dob.hpp#L161) | block | Observers & estimators | Lightweight SISO disturbance observer runtime |
| [`DOBConfig`](inc/damp/estimation/dob.hpp#L32) | block | Observers & estimators | Configuration for a first-order disturbance observer |
| [`DP45`](inc/damp/simulation/integrator.hpp#L721) | block | Simulation / SIL harness (host) | Dormand-Prince 5(4) adaptive integrator — the ode45 pair |
| [`dpwm_zero_sequence`](inc/damp/motor/modulation.hpp#L269) | function | Motor control pack (if present) | DPWM zero-sequence for a chosen discontinuous scheme |
| [`DqCommand`](inc/damp/motor/foc.hpp#L83) | block | Motor control pack (if present) | Clamped dq voltage command from FOController::current_controller |
| [`DrivePair`](inc/damp/toolbox/io.hpp#L315) | block | Embedded helpers (controls-adjacent utilities) | Left/right actuator pair (differential/arcade drive output) |
| [`DsogiPll`](inc/damp/filters/pll.hpp#L282) | block | Filters & signal conditioning | Dual-SOGI three-phase positive-sequence PLL (DSOGI-PLL) |
| [`duties_from_phase_voltages`](inc/damp/motor/modulation.hpp#L361) | function | Motor control pack (if present) | Map phase voltages + zero-sequence to clamped half-bridge duties |
| [`eig`](inc/damp/matlab.hpp#L488) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the eigenvalues of a square matrix |
| [`EigenResult`](inc/damp/matrix/eigen.hpp#L35) | block | Linear algebra | Eigenvalue computation result |
| [`EKFMeasFn`](inc/damp/estimation/ekf.hpp#L81) | concept | Observers & estimators | Concept for EKF measurement functions |
| [`EKFStateFn`](inc/damp/estimation/ekf.hpp#L53) | concept | Observers & estimators | Concept for EKF state functions |
| [`enu_from_ned`](inc/damp/estimation/ins_mechanization.hpp#L116) | function | Observers & estimators | Map a NED vector into ENU (axis permute) |
| [`eps`](inc/damp/matrix/matrix_traits.hpp#L86) | function | Linear algebra | Machine epsilon of T's real scalar (MATLAB® eps). Scale it by problem size and magnitude to build a relative tolerance |
| [`ErrorStateJacobian`](inc/damp/estimation/eskf.hpp#L281) | block | Observers & estimators | ESKF prediction Jacobians F and G (nominal state updated outside) |
| [`ErrorStateKalmanFilter`](inc/damp/estimation/eskf.hpp#L328) | block | Observers & estimators | Runtime error-state KF: tracks δx and P, not the full nominal state |
| [`ErrorTolerance`](inc/damp/simulation/solver.hpp#L286) | block | Simulation / SIL harness (host) | Scalar-or-vector local-error tolerance for adaptive ODE control |
| [`esc`](inc/damp/controllers/esc.hpp#L153) | function | Design-time synthesis (not PWM-rate) | Synthesize an extremum-seeking controller |
| [`esc_lpf_alpha`](inc/damp/controllers/esc.hpp#L131) | function | Design-time synthesis (not PWM-rate) | First-order discrete low-pass coefficient for corner wc [rad/s] at Ts |
| [`esc_mppt`](inc/damp/controllers/esc.hpp#L200) | function | Design-time synthesis (not PWM-rate) | MPPT-flavored ESC: maximize a power measurement by perturbing the operating point (e.g. converter duty or reference voltage) |
| [`ESCConfig`](inc/damp/controllers/esc.hpp#L59) | block | Runtime controllers | Extremum-seeking controller configuration (discrete realization) |
| [`ESCResult`](inc/damp/controllers/esc.hpp#L103) | block | Design-time synthesis (not PWM-rate) | Design result for the extremum-seeking controller |
| [`eskf_imu`](inc/damp/estimation/eskf.hpp#L213) | function | Observers & estimators | Design Q, R, P₀ for a 6-state IMU attitude ESKF (gyro + accel) |
| [`eskf_marg`](inc/damp/estimation/eskf.hpp#L252) | function | Observers & estimators | Design Q, R, P₀ for a 6-state MARG attitude ESKF (gyro + accel + mag) |
| [`ESKFMeasFn`](inc/damp/estimation/eskf.hpp#L304) | concept | Observers & estimators | Callable that returns a measurement linearization for one update |
| [`ESKFOrientationFilter`](inc/damp/estimation/sensor_fusion.hpp#L355) | block | Observers & estimators | Turnkey 6-state attitude ESKF: owns q, b_g, and the error filter |
| [`ESKFPredictFn`](inc/damp/estimation/eskf.hpp#L293) | concept | Observers & estimators | Callable that returns ErrorStateJacobian for one predict step |
| [`ESKFResult`](inc/damp/estimation/eskf.hpp#L84) | block | Observers & estimators | Design payload for an ESKF: Q, R, P₀, and success |
| [`estim`](inc/damp/matlab.hpp#L505) | function | MATLAB®-style aliases (host) | Form state estimator from system and estimator gain |
| [`Euler`](inc/damp/math/geometry.hpp#L205) | block | Scalar math, complex & frames | Intrinsic Euler angles for sequence Order (default ZYX aerospace YPR) |
| [`EulerOrder`](inc/damp/math/geometry.hpp#L39) | enum | Scalar math, complex & frames | Intrinsic Euler angle sequence (axis order of successive principal rotations) |
| [`eval_frf`](inc/damp/systems/state_space.hpp#L169) | function | LTI systems (SS / TF / ZPK / discretize) | Evaluate frequency response of a state-space system |
| [`EventSimConfig`](inc/damp/simulation/hybrid.hpp#L68) | block | Simulation / SIL harness (host) | Configuration for event-driven hybrid / multi-rate style runs |
| [`Exact`](inc/damp/simulation/integrator.hpp#L143) | block | Simulation / SIL harness (host) | Exact integrator for LTI systems |
| [`exact_lti_step`](inc/damp/simulation/hybrid.hpp#L295) | function | Simulation / SIL harness (host) | Exact open-loop step of a single continuous LTI system |
| [`exp`](inc/damp/math/math.hpp#L232) | function | Scalar math, complex & frames | Exponential function |
| [`ExperimentSafety`](inc/damp/estimation/experiment_safety.hpp#L62) | block | Observers & estimators | Experiment lifecycle + command protection for on-target commissioning |
| [`ExperimentSafetyConfig`](inc/damp/estimation/experiment_safety.hpp#L39) | block | Observers & estimators | Configuration for ExperimentSafety |
| [`expm`](inc/damp/matrix/functions.hpp#L701) | function | Linear algebra | Matrix exponential by scaling and squaring (Al-Mohy & Higham 2009) |
| [`expo`](inc/damp/toolbox/conditioning.hpp#L262) | function | Embedded helpers (controls-adjacent utilities) | Exponential response curve `y = (1−k)·x + k·x³` (RC "expo") |
| [`ExtendedKalmanFilter`](inc/damp/estimation/ekf.hpp#L111) | block | Observers & estimators | Extended Kalman Filter for nonlinear discrete-time systems |
| [`extract_channel`](inc/damp/simulation/plot_plotly.hpp#L223) | function | Simulation / SIL harness (host) | Extract the i-th element from each vector entry into a std::vector<double> |
| [`extract_frf_features`](inc/damp/estimation/frequency_response.hpp#L740) | function | Observers & estimators | Extract peaks and valleys in one pass over the FRF table |
| [`extract_modes`](inc/damp/estimation/frequency_response.hpp#L715) | function | Observers & estimators | Extract resonant peaks (local \|G\| maxima) from an FRF table |
| [`extract_valleys`](inc/damp/estimation/frequency_response.hpp#L729) | function | Observers & estimators | Extract anti-resonance valleys (local \|G\| minima) from an FRF table |
| [`Extrapolation`](inc/damp/toolbox/lookup.hpp#L36) | enum | Embedded helpers (controls-adjacent utilities) | Out-of-range behaviour for a Lut1D / Lut2D / Lut3D query beyond its breakpoints |
| [`ExtremumSeekingController`](inc/damp/controllers/esc.hpp#L239) | block | Runtime controllers | Extremum-seeking controller runtime (model-free online optimizer) |
| [`ExtremumType`](inc/damp/controllers/esc.hpp#L47) | enum | Runtime controllers | Whether ESC climbs to a maximum or descends to a minimum of the objective |
| [`eye`](inc/damp/matlab.hpp#L425) | function | MATLAB®-style aliases (host) | Create an identity matrix of size n x n |
| [`F_TRIG`](inc/damp/toolbox/iec61131.hpp#L154) | block | Embedded helpers (controls-adjacent utilities) | F_TRIG (Falling Edge Trigger) |
| [`factorial`](inc/damp/trajectory/trajectory_types.hpp#L206) | function | Trajectory value types | k! (exact for modest k in float/double; used up through nonic BVP / jets) |
| [`falling_factorial`](inc/damp/trajectory/trajectory_types.hpp#L217) | function | Trajectory value types | Falling factorial i·(i−1)···(i−k+1) = i! / (i−k)! — the k-th derivative coefficient of tⁱ. Zero when i < k |
| [`FallingEdge`](inc/damp/toolbox/logic.hpp#L48) | block | Embedded helpers (controls-adjacent utilities) | Falling-edge detector: true on the tick x goes true → false |
| [`feedback`](inc/damp/systems/state_space.hpp#L323) | function | LTI systems (SS / TF / ZPK / discretize) | Negative feedback: y = sys1(u − sys2(y)) (+1 more overload) |
| [`feedback_read`](inc/damp/systems/state_space.hpp#L693) | function | LTI systems (SS / TF / ZPK / discretize) | Unity (or weighted) negative feedback around G_loop, read G_out |
| [`finish_extremum`](inc/damp/estimation/frequency_response.hpp#L561) | function | Observers & estimators | Half-power (peaks) or double-power (valleys) bandwidth → ζ |
| [`finite_non_negative`](inc/damp/math/math.hpp#L449) | function | Scalar math, complex & frames | True if x is finite and non-negative (x ≥ 0) |
| [`finite_positive`](inc/damp/math/math.hpp#L439) | function | Scalar math, complex & frames | True if x is finite and strictly greater than zero |
| [`fir1`](inc/damp/filters/fir.hpp#L241) | function | Filters & signal conditioning | Window-method FIR design (normalized frequency, single cutoff) (+1 more overload) |
| [`fir1_hz`](inc/damp/filters/fir.hpp#L381) | function | Filters & signal conditioning | Window-method FIR design from frequencies in Hz (+1 more overload) |
| [`fir_window`](inc/damp/filters/fir.hpp#L417) | function | Filters & signal conditioning | Alias for fir1 (descriptive name) (+1 more overload) |
| [`FirDesignResult`](inc/damp/filters/fir.hpp#L68) | block | Filters & signal conditioning | Window-method FIR design result |
| [`fires_at`](inc/damp/simulation/multirate.hpp#L73) | function | Simulation / SIL harness (host) | True when base tick tick is a firing instant for divisor |
| [`FirFilter`](inc/damp/filters/fir.hpp#L459) | block | Filters & signal conditioning | Direct-form FIR filter runtime (tapped delay line + dot product) |
| [`FirstOrderCoeffs`](inc/damp/filters/iir_design.hpp#L30) | block | Filters & signal conditioning | DSP coefficients for first-order IIR filter |
| [`FirstOrderPlantEstimator`](inc/damp/estimation/parameter_estimation.hpp#L126) | block | Observers & estimators | Online estimator for a first-order plant's gain and time constant |
| [`FirstOrderPlantEstimatorConfig`](inc/damp/estimation/parameter_estimation.hpp#L85) | block | Observers & estimators | Configuration for the first-order grey-box estimator |
| [`FirType`](inc/damp/filters/fir.hpp#L39) | enum | Filters & signal conditioning | Ideal frequency-selective response for fir1 / fir_window |
| [`FirWindow`](inc/damp/filters/fir.hpp#L51) | enum | Filters & signal conditioning | Window applied to the ideal truncated impulse response |
| [`fit_second_order_step`](inc/damp/estimation/successive_compensator.hpp#L85) | function | Observers & estimators | Fit one underdamped second-order mode from a step / ring-down capture |
| [`fixed_solve`](inc/damp/simulation/solver.hpp#L270) | function | Simulation / SIL harness (host) | Fixed-step solve in one call |
| [`FixedStepSolver`](inc/damp/simulation/solver.hpp#L145) | block | Simulation / SIL harness (host) | Fixed-step ODE solver (+1 more overload) |
| [`floor`](inc/damp/math/math.hpp#L318) | function | Scalar math, complex & frames | Floor — largest integer ≤ x |
| [`flux_from_Kv`](inc/damp/motor/spm.hpp#L55) | function | Motor control pack (if present) | λ from Kᵥ via torque_constant_from_Kv |
| [`flux_from_torque_constant`](inc/damp/motor/spm.hpp#L37) | function | Motor control pack (if present) | λ from datasheet Kₜ (amplitude / peak-per-phase convention) |
| [`fmod`](inc/damp/math/math.hpp#L362) | function | Scalar math, complex & frames | Floating-point remainder, x − y·trunc(x/y) (sign of x), matching std::fmod's truncated-quotient convention |
| [`FOController`](inc/damp/motor/foc.hpp#L98) | block | Motor control pack (if present) | dq current regulator: PI + decoupling FF → clamped Vdq |
| [`forward_substitute`](inc/damp/matrix/solve.hpp#L37) | function | Linear algebra | Forward substitution for Lx = b |
| [`ForwardEuler`](inc/damp/simulation/integrator.hpp#L195) | block | Simulation / SIL harness (host) | Forward Euler integrator |
| [`FrequencyPoint`](inc/damp/analysis/frequency.hpp#L38) | block | Frequency-domain analysis (host) | Single-point frequency response result |
| [`FrequencyResponseEstimator`](inc/damp/estimation/frequency_response.hpp#L106) | block | Observers & estimators | Lock-in FRF estimator over a fixed table of NFreq bins |
| [`frf_pid_autotune`](inc/damp/estimation/commissioning_frontends.hpp#L152) | function | Observers & estimators | PID gains from an open-loop FRF table + margin / bandwidth targets |
| [`frf_s_or_z`](inc/damp/analysis/frequency.hpp#L987) | function | Frequency-domain analysis (host) | FRF evaluation point: jω (continuous) or e^{jωTs} (discrete) |
| [`FrfFeatureKind`](inc/damp/estimation/frequency_response.hpp#L373) | enum | Observers & estimators | Peak (resonance) or valley (anti-resonance) on an FRF magnitude curve |
| [`FrfFeatureResult`](inc/damp/estimation/frequency_response.hpp#L453) | block | Observers & estimators | Combined peak + valley reduction of an FRF table |
| [`FrfMargins`](inc/damp/estimation/frequency_response.hpp#L485) | block | Observers & estimators | Open-loop margins read from a discrete FRF table (on-target) |
| [`FrfPidAutotuneConfig`](inc/damp/estimation/commissioning_frontends.hpp#L109) | block | Observers & estimators | Configuration for FRF-based PID gain selection |
| [`FrfPidAutotuneResult`](inc/damp/estimation/commissioning_frontends.hpp#L125) | block | Observers & estimators | FRF PID autotune hand-off (gains + success) |
| [`FrfPoint`](inc/damp/estimation/frequency_response.hpp#L60) | block | Observers & estimators | One measured frequency-response point (on-target FRF table entry) |
| [`frobenius_norm`](inc/damp/matrix/functions.hpp#L99) | function | Linear algebra | Frobenius norm ‖A‖F = √(Σᵢⱼ \|aᵢⱼ\|²) |
| [`full_qr`](inc/damp/matrix/decomposition.hpp#L310) | function | Linear algebra | Full QR factorization via Householder reflections (real or complex T) |
| [`FullQR`](inc/damp/matrix/decomposition.hpp#L292) | block | Linear algebra | Result of a full (complete) QR factorization |
| [`gain_margin_unwrapped`](inc/damp/analysis/frequency.hpp#L204) | function | Frequency-domain analysis (host) | Find gain margin using unwrapped phase trajectory |
| [`geomspace`](inc/damp/analysis/linspace.hpp#L105) | function | Frequency-domain analysis (host) | Geometric sequence from start to end (endpoint values, not exponents) |
| [`Goertzel`](inc/damp/filters/spectral.hpp#L46) | block | Filters & signal conditioning | Generalized Goertzel single-bin DFT — amplitude/phase at one frequency |
| [`gram`](inc/damp/matlab.hpp#L1133) | function | MATLAB®-style aliases (host) | MATLAB® alias for the controllability/observability Gramian of a system |
| [`gravity_nav`](inc/damp/estimation/ins_mechanization.hpp#L93) | function | Observers & estimators | Gravity vector in the local-level frame [m/s²] |
| [`h_pattern`](inc/damp/toolbox/io.hpp#L377) | function | Embedded helpers (controls-adjacent utilities) | H-pattern (dual-lever / tank) drive: one lever per track, no mixing |
| [`hankel_singular_values`](inc/damp/design/model_reduction.hpp#L597) | function | Design-time synthesis (not PWM-rate) | Descriptive alias for hankelsv |
| [`hankelsv`](inc/damp/design/model_reduction.hpp#L440) | function | Design-time synthesis (not PWM-rate) | Hankel singular values of a stable state-space system |
| [`harmonic_suppressor`](inc/damp/controllers/harmonic_suppression.hpp#L113) | function | Design-time synthesis (not PWM-rate) | Synthesize a multi-resonant harmonic suppressor |
| [`HarmonicAnalyzer`](inc/damp/filters/spectral.hpp#L183) | block | Filters & signal conditioning | Harmonic analyzer — a Goertzel bank over a fundamental and K−1 harmonics |
| [`HarmonicSuppressor`](inc/damp/controllers/harmonic_suppression.hpp#L154) | block | Runtime controllers | Multi-resonant harmonic suppressor — a parallel bank of PR resonators |
| [`HarmonicSuppressorResult`](inc/damp/controllers/harmonic_suppression.hpp#L48) | block | Design-time synthesis (not PWM-rate) | Design result for a multi-resonant harmonic suppressor |
| [`heading_from_baseline_nav`](inc/damp/estimation/ins_eskf.hpp#L510) | function | Observers & estimators | Heading [rad] from a nav-frame baseline vector (dual-antenna difference) |
| [`Heun`](inc/damp/simulation/integrator.hpp#L804) | block | Simulation / SIL harness (host) | Heun's method (Improved Euler, RK2) integrator |
| [`HighPass`](inc/damp/filters/highpass.hpp#L40) | block | Filters & signal conditioning | First-order high-pass (washout) filter runtime |
| [`highpass_1st`](inc/damp/filters/iir_design.hpp#L196) | function | Filters & signal conditioning | First-order high-pass filter design (Tustin / bilinear) (+1 more overload) |
| [`highpass_2nd`](inc/damp/filters/iir_design.hpp#L388) | function | Filters & signal conditioning | Second-order continuous high-pass at default Butterworth Q (1/√2) (+1 more overload) |
| [`highpass_2nd_continuous`](inc/damp/filters/iir_design.hpp#L366) | function | Filters & signal conditioning | Second-order high-pass filter design (continuous-time) |
| [`highshelf`](inc/damp/filters/iir_design.hpp#L976) | function | Filters & signal conditioning | High-shelf EQ filter: boost or cut everything above fc |
| [`hinfnorm`](inc/damp/matlab.hpp#L1155) | function | MATLAB®-style aliases (host) | MATLAB® alias for the H∞ system norm norm(sys,Inf) / hinfnorm(sys) |
| [`holonomic_drive`](inc/damp/toolbox/io.hpp#L495) | function | Embedded helpers (controls-adjacent utilities) | Four-wheel holonomic drive mixer (mecanum / omni) |
| [`HolonomicOutput`](inc/damp/toolbox/io.hpp#L463) | block | Embedded helpers (controls-adjacent utilities) | Four-wheel actuator set (mecanum/holonomic drive output), X-config |
| [`HybridEventAction`](inc/damp/simulation/hybrid.hpp#L102) | block | Simulation / SIL harness (host) | Action applied when an event fires (after Exact advance to the event time) |
| [`HybridSimulationResult`](inc/damp/simulation/hybrid.hpp#L80) | block | Simulation / SIL harness (host) | Result of a piecewise-LTI hybrid simulation |
| [`hypot`](inc/damp/math/math.hpp#L76) | function | Scalar math, complex & frames | Euclidean distance hypot(x, y) = √(x² + y²), without overflow |
| [`Hysteresis`](inc/damp/toolbox/conditioning.hpp#L342) | block | Embedded helpers (controls-adjacent utilities) | Hysteresis comparator (Schmitt trigger): bool output with separate on/off thresholds to reject chatter |
| [`impedance`](inc/damp/analysis/frequency.hpp#L1012) | function | Frequency-domain analysis (host) | Compute impedance frequency response from a SISO admittance system |
| [`impedance_direct`](inc/damp/analysis/frequency.hpp#L1052) | function | Frequency-domain analysis (host) | Compute impedance frequency response from a SISO impedance transfer function |
| [`ImpedanceResult`](inc/damp/analysis/frequency.hpp#L885) | block | Frequency-domain analysis (host) | Result of impedance frequency response evaluation |
| [`impulse`](inc/damp/analysis/time_response.hpp#L194) | function | Frequency-domain analysis (host) | Impulse response of a (MIMO) state-space system (+1 more overload) |
| [`impulse`](inc/damp/estimation/excitation/impulse.hpp#L70) | function | Observers & estimators | Build an impulse design payload |
| [`Impulse`](inc/damp/estimation/excitation/impulse.hpp#L86) | block | Observers & estimators | Rectangular impulse (finite-width Dirac stand-in) |
| [`ImpulseConfig`](inc/damp/estimation/excitation/impulse.hpp#L32) | block | Observers & estimators | Configuration for a one-shot rectangular impulse (Dirac stand-in) |
| [`impulseplot`](inc/damp/simulation/plot_plotly.hpp#L733) | function | Simulation / SIL harness (host) | Plot an impulse response, one trace per input/output pair |
| [`ImpulseResult`](inc/damp/estimation/excitation/impulse.hpp#L52) | block | Observers & estimators | Impulse design payload |
| [`ImuSample`](inc/damp/estimation/ins_mechanization.hpp#L125) | block | Observers & estimators | IMU sample in the body frame |
| [`inductor_current_pi`](inc/damp/design/pid_design.hpp#L1194) | function | Design-time synthesis (not PWM-rate) | Fixed-rate inductor current-loop PI (topology-agnostic) |
| [`InductorCurrentPIResult`](inc/damp/design/pid_design.hpp#L1165) | block | Design-time synthesis (not PWM-rate) | Result of inductor_current_pi(): discrete PI gains plus success flag |
| [`infinity_norm`](inc/damp/matrix/functions.hpp#L38) | function | Linear algebra | Infinity norm ‖A‖∞: maximum absolute row sum |
| [`initial`](inc/damp/analysis/time_response.hpp#L242) | function | Frequency-domain analysis (host) | Initial-condition (free) response of a (MIMO) state-space system |
| [`ins_aid_position`](inc/damp/estimation/ins_eskf.hpp#L477) | function | Observers & estimators | Predicted antenna / marker position with body lever-arm |
| [`ins_apply_correction`](inc/damp/estimation/ins_eskf.hpp#L463) | function | Observers & estimators | Inject filter correction into x and reset the error state |
| [`ins_error_jacobian`](inc/damp/estimation/ins_eskf.hpp#L281) | function | Observers & estimators | Discrete error-state transition and noise input for one IMU step |
| [`ins_eskf_design`](inc/damp/estimation/ins_eskf.hpp#L193) | function | Observers & estimators | Build Q / R / P0 for a 15-state INS ESKF from IMU noise densities |
| [`ins_heading`](inc/damp/estimation/ins_eskf.hpp#L495) | function | Observers & estimators | Horizontal heading [rad] of a body axis expressed in the nav frame |
| [`ins_inject`](inc/damp/estimation/ins_eskf.hpp#L411) | function | Observers & estimators | Inject δx into the nominal INS state (right-multiplicative q) |
| [`ins_predict`](inc/damp/estimation/ins_eskf.hpp#L440) | function | Observers & estimators | Mechanize nominal state and propagate the error covariance |
| [`ins_propagate_covariance`](inc/damp/estimation/ins_eskf.hpp#L377) | function | Observers & estimators | P ← F P Fᵀ + Q using sparse F (exact for ins_error_jacobian structure) |
| [`ins_update_heading`](inc/damp/estimation/ins_eskf.hpp#L609) | function | Observers & estimators | Sparse dual-antenna / yaw heading update (scalar) |
| [`ins_update_orientation`](inc/damp/estimation/ins_eskf.hpp#L647) | function | Observers & estimators | Sparse orientation update from a measured body→nav quaternion |
| [`ins_update_pose`](inc/damp/estimation/ins_eskf.hpp#L672) | function | Observers & estimators | Generic pose alias: position then orientation (two sparse updates) |
| [`ins_update_pose_heading`](inc/damp/estimation/ins_eskf.hpp#L691) | function | Observers & estimators | Pose alias: position + dual-antenna heading (two sparse updates) |
| [`ins_update_position`](inc/damp/estimation/ins_eskf.hpp#L525) | function | Observers & estimators | Sparse position update: y = p + Rℓ (+1 more overload) |
| [`ins_update_velocity`](inc/damp/estimation/ins_eskf.hpp#L568) | function | Observers & estimators | Sparse velocity update: y = v at the IMU origin (no lever-arm) |
| [`ins_update_zupt`](inc/damp/estimation/ins_eskf.hpp#L587) | function | Observers & estimators | Zero-velocity update (ZUPT): velocity measurement of zero |
| [`InsEskfResult`](inc/damp/estimation/ins_eskf.hpp#L150) | block | Observers & estimators | Design result for the 15-state INS error-state filter |
| [`InsNavigator`](inc/damp/estimation/ins_eskf.hpp#L725) | block | Observers & estimators | Runtime strapdown navigator: nominal InsState + 15-state ESKF |
| [`InsState`](inc/damp/estimation/ins_mechanization.hpp#L139) | block | Observers & estimators | Nominal strapdown navigation state (external to the ESKF error state) |
| [`instantaneous_power`](inc/damp/math/transforms.hpp#L472) | function | Scalar math, complex & frames | Instantaneous active and reactive power from dq quantities (+1 more overload) |
| [`InstantaneousPower`](inc/damp/math/transforms.hpp#L433) | block | Scalar math, complex & frames | Instantaneous active and reactive power |
| [`IntegrationResult`](inc/damp/simulation/integrator.hpp#L57) | block | Simulation / SIL harness (host) | Result of an integration step |
| [`InteriorPointSolver`](inc/damp/design/qp.hpp#L1193) | block | Design-time synthesis (not PWM-rate) | Interior-point QP solver policy for damp::MPC |
| [`Interpolation`](inc/damp/toolbox/lookup.hpp#L42) | enum | Embedded helpers (controls-adjacent utilities) | In-range interpolant for Lut1D (2-D stays bilinear; use SplineSurface for C¹ grids) |
| [`inverse`](inc/damp/matrix/solve.hpp#L246) | function | Linear algebra | Matrix inverse A⁻¹ via mat::solve(A, I) |
| [`inverse_clarke_transform`](inc/damp/math/transforms.hpp#L284) | function | Scalar math, complex & frames | Inverse Clarke transform (αβ → abc) |
| [`inverse_deadband`](inc/damp/toolbox/conditioning.hpp#L79) | function | Embedded helpers (controls-adjacent utilities) | Inverse dead zone: add an offset to overcome a physical dead zone (valve overlap, static friction, motor stiction), with independent negative/positive offsets (+1 more overload) |
| [`inverse_deadband_sine`](inc/damp/toolbox/conditioning.hpp#L162) | function | Embedded helpers (controls-adjacent utilities) | Sine-shaped inverse dead zone (full boost once \|x\|≥w) |
| [`inverse_deadband_soft`](inc/damp/toolbox/conditioning.hpp#L138) | function | Embedded helpers (controls-adjacent utilities) | Smooth inverse dead zone via soft-sign (tanh-like Coulomb boost) |
| [`inverse_lerp`](inc/damp/toolbox/scaling.hpp#L48) | function | Embedded helpers (controls-adjacent utilities) | Inverse of lerp: the fraction t such that `lerp(a, b, t) == x` |
| [`inverse_park_clarke_transform`](inc/damp/math/transforms.hpp#L405) | function | Scalar math, complex & frames | Fused inverse Park-Clarke transform (dq → abc) |
| [`inverse_park_transform`](inc/damp/math/transforms.hpp#L358) | function | Scalar math, complex & frames | Inverse Park transform (dq → αβ) |
| [`inverse_symmetrical_components`](inc/damp/math/transforms.hpp#L551) | function | Scalar math, complex & frames | Inverse symmetrical-component transform (012 → abc) |
| [`iq_from_torque`](inc/damp/motor/spm.hpp#L71) | function | Motor control pack (if present) | i_q = T_e / Kₜ for id = 0 |
| [`is_closed_loop_stable_discrete`](inc/damp/design/stability.hpp#L288) | function | Design-time synthesis (not PWM-rate) | Check closed-loop stability for discrete system with state feedback |
| [`is_controllable`](inc/damp/design/stability.hpp#L228) | function | Design-time synthesis (not PWM-rate) | Check if a system is controllable |
| [`is_fault`](inc/damp/toolbox/conditioning.hpp#L392) | function | Embedded helpers (controls-adjacent utilities) | True for a wire fault (FaultLow/FaultHigh) — i.e. not a real reading at all |
| [`is_matrix_element`](inc/damp/matrix/matrix_traits.hpp#L42) | block | Linear algebra | True if T may be a Matrix element type |
| [`is_observable`](inc/damp/design/stability.hpp#L242) | function | Design-time synthesis (not PWM-rate) | Check if a system is observable |
| [`is_stabilizable`](inc/damp/design/riccati.hpp#L45) | function | Design-time synthesis (not PWM-rate) | Check if (A, B) is a stabilizable pair |
| [`is_stable_continuous`](inc/damp/analysis/poles.hpp#L58) | function | Frequency-domain analysis (host) | Check continuous-time stability |
| [`is_stable_discrete`](inc/damp/design/stability.hpp#L258) | function | Design-time synthesis (not PWM-rate) | Check if a discrete-time system matrix A is stable |
| [`is_symmetric_or_hermitian`](inc/damp/matrix/decomposition.hpp#L36) | function | Linear algebra | True if A equals its transpose (real) or conjugate transpose (complex) |
| [`is_usable`](inc/damp/toolbox/conditioning.hpp#L549) | function | Embedded helpers (controls-adjacent utilities) | True when a channel reading may still be used (in-span or saturated) |
| [`is_valid`](inc/damp/toolbox/conditioning.hpp#L387) | function | Embedded helpers (controls-adjacent utilities) | True only for the in-span status |
| [`isfinite`](inc/damp/math/math.hpp#L426) | function | Scalar math, complex & frames | Finiteness test — false for NaN and ±∞ |
| [`iso_c_drive`](inc/damp/toolbox/io.hpp#L412) | function | Embedded helpers (controls-adjacent utilities) | ISO-C (single-stick, coordinated / curvature) drive |
| [`iso_s_drive`](inc/damp/toolbox/io.hpp#L393) | function | Embedded helpers (controls-adjacent utilities) | ISO-S (single-stick, speed-summed) drive — arcade travel stick |
| [`isstable`](inc/damp/matlab.hpp#L971) | function | MATLAB®-style aliases (host) | Continuous-time stability predicate on a state matrix (+2 more overloads) |
| [`Jet`](inc/damp/trajectory/trajectory_types.hpp#L95) | block | Trajectory value types | Truncated jet of a scalar motion sample: d[k] = s⁽ᵏ⁾ |
| [`JordanBlock`](inc/damp/design/pole_placement.hpp#L845) | block | Design-time synthesis (not PWM-rate) | One Jordan mini-block of a desired closed-loop spectrum |
| [`JordanObjective`](inc/damp/design/pole_placement.hpp#L1179) | enum | Design-time synthesis (not PWM-rate) | Robustness objective for place_jordan_optimal (the paper's two methods) |
| [`kalman`](inc/damp/estimation/kalman.hpp#L103) | function | Observers & estimators | Steady-state Kalman filter design |
| [`KalmanFilter`](inc/damp/estimation/kalman.hpp#L232) | block | Observers & estimators | Runtime Kalman filter for embedded systems |
| [`KalmanResult`](inc/damp/estimation/kalman.hpp#L45) | block | Observers & estimators | Steady-state Kalman filter design result |
| [`lag`](inc/damp/controllers/lead_lag.hpp#L181) | function | Design-time synthesis (not PWM-rate) | Design a lag compensator from desired low-frequency gain boost |
| [`lambda_tuning`](inc/damp/design/pid_design.hpp#L403) | function | Design-time synthesis (not PWM-rate) | Lambda tuning for FOPDT model |
| [`Latch`](inc/damp/toolbox/logic.hpp#L66) | block | Embedded helpers (controls-adjacent utilities) | Set/reset latch. @tparam SetDominant which input wins when both are asserted (default: set-dominant, e.g. a trip overriding a clear) |
| [`lead`](inc/damp/controllers/lead_lag.hpp#L132) | function | Design-time synthesis (not PWM-rate) | Design a lead compensator from desired phase boost at a target frequency |
| [`lead_lag`](inc/damp/controllers/lead_lag.hpp#L234) | function | Design-time synthesis (not PWM-rate) | Design a lead-lag compensator (cascade of lead + lag sections) |
| [`lead_lag_direct`](inc/damp/controllers/lead_lag.hpp#L264) | function | Design-time synthesis (not PWM-rate) | Direct lead-lag specification from zero/pole locations |
| [`LeadLagController`](inc/damp/controllers/lead_lag.hpp#L281) | block | Runtime controllers | Discrete lead-lag compensator |
| [`LeadLagResult`](inc/damp/controllers/lead_lag.hpp#L59) | block | Design-time synthesis (not PWM-rate) | Lead-lag compensator design result |
| [`LeadLagSeriesResult`](inc/damp/controllers/lead_lag.hpp#L204) | block | Design-time synthesis (not PWM-rate) | Cascaded lead+lag design result (2nd-order StateSpace + success) |
| [`lerp`](inc/damp/toolbox/scaling.hpp#L35) | function | Embedded helpers (controls-adjacent utilities) | Linear interpolation between a and b by fraction t |
| [`leverrier_num_scale`](inc/damp/systems/zpk.hpp#L562) | function | LTI systems (SS / TF / ZPK / discretize) | Magnitude scale of each strictly proper Leverrier numerator coefficient |
| [`LeverrierResult`](inc/damp/systems/zpk.hpp#L496) | block | LTI systems (SS / TF / ZPK / discretize) | Faddeev–LeVerrier characteristic polynomial and adjoint coefficient matrices |
| [`LIMIT`](inc/damp/toolbox/iec61131.hpp#L593) | function | Embedded helpers (controls-adjacent utilities) | LIMIT (IEC 61131-3 selection function): clamp in to [mn, mx] |
| [`linear_modulation_voltage`](inc/damp/motor/modulation.hpp#L127) | function | Motor control pack (if present) | Peak phase voltage at the linear SVPWM hexagon limit |
| [`linear_screw`](inc/damp/toolbox/actuator.hpp#L166) | function | Embedded helpers (controls-adjacent utilities) | Build a ServoAxis for a linear axis driven by a leadscrew/belt |
| [`LinearizationResult`](inc/damp/design/linearization.hpp#L44) | block | Design-time synthesis (not PWM-rate) | Result of nonlinear operating-point linearization |
| [`linearize`](inc/damp/design/linearization.hpp#L170) | function | Design-time synthesis (not PWM-rate) | Linearize nonlinear dynamics and output maps about an operating point (+1 more overload) |
| [`linmod`](inc/damp/matlab.hpp#L341) | function | MATLAB®-style aliases (host) | MATLAB®-style nonlinear linearization about an operating point |
| [`linspace`](inc/damp/analysis/linspace.hpp#L38) | function | Frequency-domain analysis (host) | num evenly spaced values from start to end inclusive (+1 more overload) |
| [`locate_zero_crossing_exact`](inc/damp/simulation/hybrid.hpp#L319) | function | Simulation / SIL harness (host) | Locate a state zero-crossing of g(x) on an Exact LTI segment |
| [`log`](inc/damp/matrix/functions.hpp#L823) | function | Linear algebra | Principal matrix logarithm via inverse scaling and squaring |
| [`log`](inc/damp/math/math.hpp#L249) | function | Scalar math, complex & frames | Natural logarithm |
| [`log10`](inc/damp/math/math.hpp#L381) | function | Scalar math, complex & frames | Base-10 logarithm, log10(x) = ln(x) / ln(10) |
| [`logm`](inc/damp/matrix/functions.hpp#L866) | function | Linear algebra | @brief MATLAB®-style alias for log |
| [`LogPolicy`](inc/damp/simulation/hybrid.hpp#L59) | enum | Simulation / SIL harness (host) | How often to record samples along a hybrid trajectory |
| [`logspace`](inc/damp/analysis/linspace.hpp#L69) | function | Frequency-domain analysis (host) | num logarithmically spaced values from start to end (endpoint values) (+1 more overload) |
| [`loop_metrics`](inc/damp/analysis/frequency.hpp#L493) | function | Frequency-domain analysis (host) | One-call loop analysis: compute L/S/T response and return compact metrics (+1 more overload) |
| [`loop_response`](inc/damp/analysis/frequency.hpp#L571) | function | Frequency-domain analysis (host) | Compute open-loop L, sensitivity S, complementary sensitivity T, and Nyquist data (+1 more overload) |
| [`LoopResponseResult`](inc/damp/analysis/frequency.hpp#L428) | block | Frequency-domain analysis (host) | Open-loop and closed-loop frequency response package |
| [`LoopSummary`](inc/damp/analysis/frequency.hpp#L461) | block | Frequency-domain analysis (host) | Compact loop summary metrics for quick stability/robustness checks |
| [`LowerTriangle`](inc/damp/matrix/views.hpp#L157) | block | Linear algebra | Lower triangular view of a square matrix |
| [`LowPass`](inc/damp/filters/lowpass.hpp#L49) | block | Filters & signal conditioning | Nth-order low-pass filter |
| [`lowpass_1st`](inc/damp/filters/iir_design.hpp#L161) | function | Filters & signal conditioning | First-order low-pass filter design (+1 more overload) |
| [`lowpass_2nd`](inc/damp/filters/iir_design.hpp#L249) | function | Filters & signal conditioning | Second-order low-pass filter design (+1 more overload) |
| [`lowpass_2nd_continuous`](inc/damp/filters/iir_design.hpp#L300) | function | Filters & signal conditioning | Second-order low-pass filter design (continuous-time) |
| [`lowshelf`](inc/damp/filters/iir_design.hpp#L944) | function | Filters & signal conditioning | Low-shelf EQ filter: boost or cut everything below fc |
| [`lqg`](inc/damp/matlab.hpp#L760) | function | MATLAB®-style aliases (host) | Linear-Quadratic-Gaussian regulator design |
| [`LQG`](inc/damp/controllers/lqg.hpp#L181) | block | Runtime controllers | Linear-Quadratic-Gaussian (LQG) controller |
| [`lqg_bundle`](inc/damp/design/synthesis.hpp#L332) | function | Design-time synthesis (not PWM-rate) | Synthesize the full LQG artifact bundle in one call |
| [`lqg_from_parts`](inc/damp/controllers/lqg.hpp#L154) | function | Design-time synthesis (not PWM-rate) | Assemble an LQG design from separately computed Kalman and LQR results |
| [`lqg_pr_bundle`](inc/damp/design/synthesis.hpp#L359) | function | Design-time synthesis (not PWM-rate) | Synthesize a SISO LQG + PR design with internal-model compensation |
| [`LQGAnalysisModels`](inc/damp/design/synthesis.hpp#L39) | block | Design-time synthesis (not PWM-rate) | Analysis-oriented models produced from an LQG design |
| [`LQGArtifacts`](inc/damp/design/synthesis.hpp#L89) | block | Design-time synthesis (not PWM-rate) | Synthesis artifact bundle: design + analysis models + runtime bundle |
| [`LQGI`](inc/damp/controllers/lqgi.hpp#L167) | block | Runtime controllers | Linear-Quadratic-Gaussian-Integral (LQGI) controller |
| [`lqgi_bundle`](inc/damp/design/synthesis.hpp#L406) | function | Design-time synthesis (not PWM-rate) | Synthesize the full LQGI servo artifact bundle in one call |
| [`LQGIAnalysisModels`](inc/damp/design/synthesis.hpp#L152) | block | Design-time synthesis (not PWM-rate) | Analysis-oriented models produced from an LQGI design |
| [`LQGIArtifacts`](inc/damp/design/synthesis.hpp#L180) | block | Design-time synthesis (not PWM-rate) | Synthesis artifact bundle: LQGI servo design + analysis + ready-to-run controller |
| [`LQGIResult`](inc/damp/controllers/lqgi.hpp#L34) | block | Design-time synthesis (not PWM-rate) | LQGI design result |
| [`LQGPRArtifacts`](inc/damp/design/synthesis.hpp#L133) | block | Design-time synthesis (not PWM-rate) | Synthesis artifact bundle for SISO LQG + PR |
| [`LQGPRRuntimeBundle`](inc/damp/design/synthesis.hpp#L100) | block | Design-time synthesis (not PWM-rate) | Runtime bundle for SISO LQG + PR internal model compensation |
| [`lqgreg`](inc/damp/matlab.hpp#L779) | function | MATLAB®-style aliases (host) | Combine separate Kalman filter and LQR designs into an LQG controller |
| [`LQGResult`](inc/damp/controllers/lqg.hpp#L46) | block | Design-time synthesis (not PWM-rate) | LQG design result |
| [`LQGRuntimeBundle`](inc/damp/design/synthesis.hpp#L50) | block | Design-time synthesis (not PWM-rate) | Runtime bundle for LQG control |
| [`lqgtrack`](inc/damp/matlab.hpp#L794) | function | MATLAB®-style aliases (host) | Linear-Quadratic-Gaussian design with integral action for tracking |
| [`lqi`](inc/damp/matlab.hpp#L741) | function | MATLAB®-style aliases (host) | Linear-Quadratic Integral design for tracking |
| [`LQI`](inc/damp/controllers/lqi.hpp#L155) | block | Runtime controllers | Linear-Quadratic-Integral (LQI) controller |
| [`lqi_bundle`](inc/damp/design/synthesis.hpp#L388) | function | Design-time synthesis (not PWM-rate) | Synthesize the full LQI servo artifact bundle in one call |
| [`LQIAnalysisModels`](inc/damp/design/synthesis.hpp#L143) | block | Design-time synthesis (not PWM-rate) | Analysis-oriented models produced from an LQI design |
| [`LQIArtifacts`](inc/damp/design/synthesis.hpp#L165) | block | Design-time synthesis (not PWM-rate) | Synthesis artifact bundle: LQI servo design + analysis + ready-to-run controller |
| [`LQIResult`](inc/damp/controllers/lqi.hpp#L48) | block | Design-time synthesis (not PWM-rate) | LQI design result |
| [`lqr`](inc/damp/matlab.hpp#L612) | function | MATLAB®-style aliases (host) | Continuous-time LQR design (MATLAB®'s lqr) |
| [`lqr_gain`](inc/damp/design/riccati.hpp#L1235) | function | Design-time synthesis (not PWM-rate) | Optimal LQR state-feedback gain from a Riccati solution |
| [`LQRCost`](inc/damp/controllers/lqr.hpp#L148) | block | Runtime controllers | Discretized LQR cost weights (Q, R, N) for a sampled-data problem |
| [`lqrd`](inc/damp/controllers/lqr.hpp#L331) | function | Design-time synthesis (not PWM-rate) | Sampled-data LQR from continuous plant (MATLAB®-style short name) (+1 more overload) |
| [`lqrd`](inc/damp/matlab.hpp#L646) | function | MATLAB®-style aliases (host) | Design discrete LQR from continuous-time system via discretization (+1 more overload) |
| [`LQRResult`](inc/damp/controllers/lqr.hpp#L56) | block | Design-time synthesis (not PWM-rate) | Linear-Quadratic Regulator design result |
| [`lqry`](inc/damp/matlab.hpp#L693) | function | MATLAB®-style aliases (host) | Output-weighted continuous LQR (state cost Q = Cᵀ Q_y C) (+1 more overload) |
| [`lsim`](inc/damp/analysis/time_response.hpp#L280) | function | Frequency-domain analysis (host) | Forced time response of a (MIMO) state-space system to an input signal (+1 more overload) |
| [`LsimInfo`](inc/damp/analysis/time_response.hpp#L509) | block | Frequency-domain analysis (host) | Transient characteristics of an arbitrary response signal |
| [`lsiminfo`](inc/damp/analysis/time_response.hpp#L529) | function | Frequency-domain analysis (host) | Compute transient characteristics from an output/time signal |
| [`lsimplot`](inc/damp/simulation/plot_plotly.hpp#L749) | function | Simulation / SIL harness (host) | Plot a forced (lsim) simulation, one trace per output |
| [`LsimResult`](inc/damp/analysis/time_response.hpp#L66) | block | Frequency-domain analysis (host) | Result of a single-trajectory simulation: time, output, and state history |
| [`lu_decomposition`](inc/damp/matrix/decomposition.hpp#L125) | function | Linear algebra | LU decomposition with partial pivoting |
| [`lu_solve`](inc/damp/matrix/solve.hpp#L185) | function | Linear algebra | Solve AX = B via LU with partial pivoting |
| [`luenberger`](inc/damp/estimation/luenberger.hpp#L138) | function | Observers & estimators | Design a Luenberger observer by robust pole placement (matrix form) (+3 more overloads) |
| [`Luenberger`](inc/damp/estimation/luenberger.hpp#L400) | block | Observers & estimators | Luenberger state observer (runtime) |
| [`LuenbergerResult`](inc/damp/estimation/luenberger.hpp#L92) | block | Observers & estimators | Luenberger observer design result |
| [`Lut1D`](inc/damp/toolbox/lookup.hpp#L109) | block | Embedded helpers (controls-adjacent utilities) | 1-D interpolating lookup table over monotonic breakpoints |
| [`Lut2D`](inc/damp/toolbox/lookup.hpp#L201) | block | Embedded helpers (controls-adjacent utilities) | 2-D bilinear interpolating lookup table over a regular grid |
| [`Lut3D`](inc/damp/toolbox/lookup.hpp#L275) | block | Embedded helpers (controls-adjacent utilities) | 3-D trilinear interpolating lookup table over a regular grid |
| [`lut_segment`](inc/damp/toolbox/lookup.hpp#L63) | function | Embedded helpers (controls-adjacent utilities) | Index of the interpolation segment containing x |
| [`lyap`](inc/damp/design/lyapunov.hpp#L105) | function | Design-time synthesis (not PWM-rate) | Solve the continuous-time Lyapunov equation A X + X Aᵀ + Q = 0 |
| [`lyap`](inc/damp/matlab.hpp#L1067) | function | MATLAB®-style aliases (host) | MATLAB® alias for the continuous Lyapunov solve AX+XAᵀ+Q=0 |
| [`MadgwickFilter`](inc/damp/estimation/sensor_fusion.hpp#L114) | block | Observers & estimators | Madgwick gradient-descent AHRS filter |
| [`mag2db`](inc/damp/math/math.hpp#L460) | function | Scalar math, complex & frames | Magnitude to decibels, 20·log10(mag) |
| [`MahonyFilter`](inc/damp/estimation/sensor_fusion.hpp#L155) | block | Observers & estimators | Mahony nonlinear complementary filter with PI correction |
| [`make1DOF`](inc/damp/matlab.hpp#L232) | function | MATLAB®-style aliases (host) | Force 1-DOF setpoint weights on a PID design result (b=c=1) |
| [`make2DOF`](inc/damp/matlab.hpp#L248) | function | MATLAB®-style aliases (host) | Apply 2-DOF setpoint weights on a PID design result |
| [`margin`](inc/damp/matlab.hpp#L904) | function | MATLAB®-style aliases (host) | Gain and phase margins of a SISO loop over a frequency grid |
| [`MarginResult`](inc/damp/matlab.hpp#L883) | block | MATLAB®-style aliases (host) | Gain/phase margins and their crossover frequencies |
| [`margins_from_frf`](inc/damp/estimation/frequency_response.hpp#L772) | function | Observers & estimators | Gain / phase margins from a fixed on-target FRF table |
| [`Matrix`](inc/damp/matrix/core.hpp#L159) | block | Linear algebra | Fixed-size, stack-allocated matrix for linear algebra operations |
| [`MatrixLike`](inc/damp/matrix/matrix_traits.hpp#L71) | concept | Linear algebra | Concept for any type that provides 2D matrix-like element access |
| [`MatrixLikeOf`](inc/damp/matrix/matrix_traits.hpp#L81) | concept | Linear algebra | Concept for a MatrixLike type with specific dimensions |
| [`max`](inc/damp/matrix/core.hpp#L1012) | function | Linear algebra | Addition of two MatrixLike types (with broadcasting support) (+1 more overload) |
| [`max_iq`](inc/damp/motor/spm.hpp#L111) | function | Motor control pack (if present) | Max positive i_q for an SPM (i_d = 0) at electrical speed |
| [`max_torque_at_speed`](inc/damp/motor/spm.hpp#L172) | function | Motor control pack (if present) | SPM max-torque point at a mechanical speed (i_d = 0 only) |
| [`MeasJacobian`](inc/damp/estimation/ekf.hpp#L68) | block | Observers & estimators | Measurement prediction result from the user's observation function |
| [`mecanum_drive`](inc/damp/toolbox/io.hpp#L514) | function | Embedded helpers (controls-adjacent utilities) | Mecanum (H-layout, 45° rollers) mixer — alias for holonomic_drive |
| [`mechanize_step`](inc/damp/estimation/ins_mechanization.hpp#L210) | function | Observers & estimators | One dead-reckoning step from a raw IMU sample (subtracts the state's biases) |
| [`mechanize_step_from_corrected`](inc/damp/estimation/ins_mechanization.hpp#L175) | function | Observers & estimators | One dead-reckoning step from an IMU sample (strapdown mechanization) |
| [`MedianFilter`](inc/damp/filters/median.hpp#L34) | block | Filters & signal conditioning | Sliding-window median filter — nonlinear spike/outlier rejection |
| [`mhe`](inc/damp/estimation/mhe.hpp#L176) | function | Observers & estimators | Synthesize a constrained moving-horizon estimator |
| [`MHE`](inc/damp/estimation/mhe.hpp#L298) | block | Observers & estimators | Runtime moving-horizon estimator (fixed per-tick iteration budget) (+1 more overload) |
| [`MHEArtifacts`](inc/damp/estimation/mhe.hpp#L114) | block | Observers & estimators | Window QP data produced by mhe(), consumed by damp::MHE |
| [`MHEConstraints`](inc/damp/estimation/mhe.hpp#L91) | block | Observers & estimators | Box constraints on the estimated states (applied to every window state) |
| [`middlebrook`](inc/damp/analysis/frequency.hpp#L1096) | function | Frequency-domain analysis (host) | Middlebrook stability analysis for cascaded source-load systems (+1 more overload) |
| [`MiddlebrookResult`](inc/damp/analysis/frequency.hpp#L909) | block | Frequency-domain analysis (host) | Result of Middlebrook minor loop gain analysis |
| [`minimal_realization`](inc/damp/design/minreal.hpp#L434) | function | Design-time synthesis (not PWM-rate) | Descriptive alias for minreal |
| [`MinimalRealizationResult`](inc/damp/design/minreal.hpp#L84) | block | Design-time synthesis (not PWM-rate) | Minimal realization result |
| [`minmax`](inc/damp/backend.hpp#L164) | function | Core, configuration & backend vocabulary | Ordered {min, max} pair returned by value (+1 more overload) |
| [`minreal`](inc/damp/design/minreal.hpp#L284) | function | Design-time synthesis (not PWM-rate) | Minimal realization — cancel uncontrollable and unobservable modes |
| [`minreal`](inc/damp/matlab.hpp#L262) | function | MATLAB®-style aliases (host) | MATLAB® short alias for design::minreal |
| [`minreal_zpk`](inc/damp/systems/zpk.hpp#L882) | function | LTI systems (SS / TF / ZPK / discretize) | Cancel matching pole-zero pairs on a ZPK model |
| [`mix`](inc/damp/systems/state_space.hpp#L580) | function | LTI systems (SS / TF / ZPK / discretize) | Fold plant outputs into one input: v = g u + Σ kᵢ yᵢ (+1 more overload) |
| [`mix_outputs`](inc/damp/systems/state_space.hpp#L650) | function | LTI systems (SS / TF / ZPK / discretize) | Linear combination of outputs: y = Σ αⱼ yⱼ |
| [`ModeExtractorConfig`](inc/damp/estimation/frequency_response.hpp#L411) | block | Observers & estimators | Configuration for FRF peak / valley extraction |
| [`ModeExtractorResult`](inc/damp/estimation/frequency_response.hpp#L428) | block | Observers & estimators | Result of reducing an FRF table to at most MaxModes peaks (resonances) |
| [`ModelReductionMethod`](inc/damp/design/model_reduction.hpp#L58) | enum | Design-time synthesis (not PWM-rate) | Method for eliminating states in modred |
| [`modred`](inc/damp/design/model_reduction.hpp#L519) | function | Design-time synthesis (not PWM-rate) | Model reduction by truncation or DC-matched residualization |
| [`modulation_duty_cycles`](inc/damp/motor/modulation.hpp#L400) | function | Motor control pack (if present) | Carrier-based VSI duty cycles from an αβ voltage command |
| [`modulation_index`](inc/damp/motor/modulation.hpp#L165) | function | Motor control pack (if present) | Modulation index relative to the linear SVPWM circle |
| [`modulation_zero_sequence`](inc/damp/motor/modulation.hpp#L327) | function | Motor control pack (if present) | Zero-sequence for any PwmScheme |
| [`Modulator`](inc/damp/motor/modulation.hpp#L497) | block | Motor control pack (if present) | Thin scheme-holding wrapper for FOC / deploy paths |
| [`motor_constant`](inc/damp/motor/spm.hpp#L63) | function | Motor control pack (if present) | Motor constant Kₘ = Kₜ / √R [Nm/√W] |
| [`MovingAverage`](inc/damp/filters/moving_average.hpp#L62) | block | Filters & signal conditioning | Moving-average (boxcar) filter — also a DC-preserving harmonic-notch comb |
| [`MPC`](inc/damp/controllers/mpc.hpp#L933) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its artifacts: MPC controller{art}; (+1 more overload) |
| [`mpc`](inc/damp/controllers/offset_free_mpc.hpp#L132) | function | Design-time synthesis (not PWM-rate) | Synthesize an offset-free constrained MPC (controller + estimator) |
| [`MPC`](inc/damp/controllers/mpc.hpp#L786) | block | Runtime controllers | Runtime constrained MPC controller (fixed per-tick iteration budget) |
| [`MPCAnalysisModels`](inc/damp/controllers/mpc.hpp#L681) | block | Runtime controllers | LTI models of the unconstrained MPC loop, for margin/robustness analysis |
| [`MPCArtifacts`](inc/damp/controllers/mpc.hpp#L241) | block | Runtime controllers | Condensed-QP data produced by mpc(), consumed by damp::MPC |
| [`MPCConstraints`](inc/damp/controllers/mpc.hpp#L203) | block | Runtime controllers | Box constraints for mpc() |
| [`MPCHorizonSuggestion`](inc/damp/controllers/mpc.hpp#L569) | block | Runtime controllers | Advisory NP/NC horizon values from suggest_mpc_horizon() |
| [`MPCWeights`](inc/damp/controllers/mpc.hpp#L168) | block | Runtime controllers | Cost weights for mpc() |
| [`mstogi`](inc/damp/filters/sogi.hpp#L120) | function | Filters & signal conditioning | Mixed Second/Third-Order Generalized Integrator (MSTOGI) (+1 more overload) |
| [`MSTOGI`](inc/damp/filters/sogi.hpp#L252) | block | Filters & signal conditioning | Runtime MSTOGI with exact resonator and forward-Euler washout |
| [`multi_sine`](inc/damp/estimation/excitation/multi_sine.hpp#L133) | function | Observers & estimators | Build a multi-sine design payload from a configuration |
| [`MultiPRController`](inc/damp/controllers/pr.hpp#L369) | block | Runtime controllers | Multi-harmonic PR Controller |
| [`MultiRateConfig`](inc/damp/simulation/multirate.hpp#L81) | block | Simulation / SIL harness (host) | Configuration for multi-rate closed-loop runs |
| [`MultiSine`](inc/damp/estimation/excitation/multi_sine.hpp#L200) | block | Observers & estimators | Sum-of-tones multi-sine runtime generator |
| [`MultiSineConfig`](inc/damp/estimation/excitation/multi_sine.hpp#L62) | block | Observers & estimators | Configuration for fixed-component multi-sine excitation |
| [`MultiSineResult`](inc/damp/estimation/excitation/multi_sine.hpp#L104) | block | Observers & estimators | Multi-sine design payload |
| [`MUX`](inc/damp/toolbox/iec61131.hpp#L612) | function | Embedded helpers (controls-adjacent utilities) | MUX (IEC 61131-3 multiplexer): select input k of N (0-based) |
| [`NavFrame`](inc/damp/estimation/ins_mechanization.hpp#L77) | enum | Observers & estimators | Local-level navigation frame for strapdown mechanization |
| [`nearbyint`](inc/damp/math/math.hpp#L346) | function | Scalar math, complex & frames | Round to nearest integer; ties to even (IEEE default / `FE_TONEAREST`) |
| [`ned_from_enu`](inc/damp/estimation/ins_mechanization.hpp#L106) | function | Observers & estimators | Map an ENU vector into NED (axis permute) |
| [`negative_sequence_ab`](inc/damp/filters/pll.hpp#L244) | function | Filters & signal conditioning | Instantaneous negative-sequence αβ from a quadrature signal pair |
| [`nichols`](inc/damp/analysis/frequency.hpp#L696) | function | Frequency-domain analysis (host) | Build Nichols points from existing Bode data (+2 more overloads) |
| [`nicholsplot`](inc/damp/simulation/plot_plotly.hpp#L936) | function | Simulation / SIL harness (host) | Plot a Nichols chart (open-loop phase vs magnitude) with M-circle grid |
| [`NicholsPoint`](inc/damp/analysis/frequency.hpp#L667) | block | Frequency-domain analysis (host) | Single-point Nichols chart sample (open-loop phase vs magnitude) |
| [`NicholsResult`](inc/damp/analysis/frequency.hpp#L679) | block | Frequency-domain analysis (host) | Nichols chart locus across a frequency sweep |
| [`NlmsFilter`](inc/damp/filters/fir.hpp#L557) | block | Filters & signal conditioning | Normalized LMS adaptive FIR filter |
| [`norm`](inc/damp/matlab.hpp#L1146) | function | MATLAB®-style aliases (host) | MATLAB® alias for the H2 system norm norm(sys,2) |
| [`norm_h2`](inc/damp/analysis/norms.hpp#L76) | function | Frequency-domain analysis (host) | H2 norm of a state-space system |
| [`norm_hinf`](inc/damp/analysis/norms.hpp#L134) | function | Frequency-domain analysis (host) | H∞ norm of a state-space system: sup_ω σ̄(G(jω)) |
| [`notch`](inc/damp/filters/iir_design.hpp#L830) | function | Filters & signal conditioning | Second-order band-reject (notch) filter |
| [`notch_continuous`](inc/damp/filters/iir_design.hpp#L434) | function | Filters & signal conditioning | Second-order band-reject (notch) filter design (continuous-time) |
| [`null`](inc/damp/matlab.hpp#L474) | function | MATLAB®-style aliases (host) | MATLAB® short alias for an orthonormal null-space basis |
| [`null_space`](inc/damp/matrix/svd.hpp#L381) | function | Linear algebra | Orthonormal basis for the null space {x : A·x = 0} via SVD |
| [`NullSpace`](inc/damp/matrix/svd.hpp#L362) | block | Linear algebra | Orthonormal basis for the null space (kernel) of a matrix |
| [`nyquist`](inc/damp/analysis/frequency.hpp#L510) | function | Frequency-domain analysis (host) | Compute Nyquist data for a SISO state-space system (+1 more overload) |
| [`nyquistplot`](inc/damp/simulation/plot_plotly.hpp#L821) | function | Simulation / SIL harness (host) | Plot a Nyquist locus with the -1 critical point marked |
| [`NyquistPoint`](inc/damp/analysis/frequency.hpp#L382) | block | Frequency-domain analysis (host) | Single-point Nyquist response data |
| [`NyquistResult`](inc/damp/analysis/frequency.hpp#L394) | block | Frequency-domain analysis (host) | Nyquist response data across a frequency sweep |
| [`observability_gramian`](inc/damp/design/stability.hpp#L189) | function | Design-time synthesis (not PWM-rate) | Continuous/discrete observability Gramian W_o |
| [`observability_matrix`](inc/damp/design/stability.hpp#L82) | function | Design-time synthesis (not PWM-rate) | Compute the observability matrix [C; CA; CA²; ...; CA^(N-1)] |
| [`obsv`](inc/damp/matlab.hpp#L286) | function | MATLAB®-style aliases (host) | MATLAB® short alias for observability_matrix (+1 more overload) |
| [`OffDelayTimer`](inc/damp/toolbox/logic.hpp#L136) | block | Embedded helpers (controls-adjacent utilities) | Off-delay timer: output goes true immediately when in is true and stays true until in has been false continuously for delay |
| [`OffsetFreeMPC`](inc/damp/controllers/offset_free_mpc.hpp#L262) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its artifacts: OffsetFreeMPC controller{art}; |
| [`OffsetFreeMPC`](inc/damp/controllers/offset_free_mpc.hpp#L198) | block | Runtime controllers | Runtime offset-free MPC: constrained MPC + disturbance-augmented Kalman filter |
| [`OffsetFreeMPCArtifacts`](inc/damp/controllers/offset_free_mpc.hpp#L84) | block | Runtime controllers | Combined MPC + disturbance-augmented Kalman design, consumed by damp::OffsetFreeMPC |
| [`OFLQG`](inc/damp/controllers/offset_free_lqg.hpp#L158) | block | Runtime controllers | Offset-free LQG: u = u_t − K(x̂ − x_t) with ŵ in the DC target |
| [`OFLQGResult`](inc/damp/controllers/offset_free_lqg.hpp#L62) | block | Design-time synthesis (not PWM-rate) | Offset-free LQG design: LQR on the plant, Kalman on [x; w] |
| [`omniwheel_drive`](inc/damp/toolbox/io.hpp#L520) | function | Embedded helpers (controls-adjacent utilities) | Omni-wheel (X-layout) mixer — alias for holonomic_drive |
| [`OnDelayTimer`](inc/damp/toolbox/logic.hpp#L100) | block | Embedded helpers (controls-adjacent utilities) | On-delay timer: output goes true once in has been held true continuously for delay; drops immediately when in goes false |
| [`one_norm`](inc/damp/matrix/functions.hpp#L67) | function | Linear algebra | One-norm ‖A‖₁: maximum absolute column sum |
| [`ones`](inc/damp/matlab.hpp#L434) | function | MATLAB®-style aliases (host) | N×N matrix of ones |
| [`OptimalJordanPlacement`](inc/damp/design/pole_placement.hpp#L1186) | block | Design-time synthesis (not PWM-rate) | Result of optimized arbitrary pole placement (place_jordan_optimal) |
| [`OutputFeedbackController`](inc/damp/concepts.hpp#L87) | concept | Core, configuration & backend vocabulary | Vector output-feedback controller: u = control(r, y), self-contained tick |
| [`pade`](inc/damp/matlab.hpp#L1013) | function | MATLAB®-style aliases (host) | First-order Padé approximation of pure delay e^{−sT} (+1 more overload) |
| [`pade2`](inc/damp/matlab.hpp#L1035) | function | MATLAB®-style aliases (host) | Second-order Padé approximation of pure delay e^{−sT} (+1 more overload) |
| [`pade_delay_1st`](inc/damp/filters/iir_design.hpp#L560) | function | Filters & signal conditioning | First-order Pade approximation of time delay (+1 more overload) |
| [`pade_delay_2nd`](inc/damp/filters/iir_design.hpp#L604) | function | Filters & signal conditioning | Second-order Pade approximation of time delay (+1 more overload) |
| [`panel_color`](inc/damp/simulation/plot_plotly.hpp#L314) | function | Simulation / SIL harness (host) | Plotly default colorway entry for index i within a single subplot |
| [`panel_legend`](inc/damp/simulation/plot_plotly.hpp#L297) | function | Simulation / SIL harness (host) | Per-panel legend box just right of a stacked subplot band |
| [`panel_legend_id`](inc/damp/simulation/plot_plotly.hpp#L283) | function | Simulation / SIL harness (host) | Legend id for subplot row (1-based): "legend", "legend2", … |
| [`panel_x_domain`](inc/damp/simulation/plot_plotly.hpp#L273) | function | Simulation / SIL harness (host) | Horizontal domain for stacked multi-panel x-axes |
| [`parallel`](inc/damp/systems/state_space.hpp#L260) | function | LTI systems (SS / TF / ZPK / discretize) | Parallel connection (shared input, summed outputs) |
| [`ParameterDriftMonitor`](inc/damp/estimation/parameter_estimation.hpp#L252) | block | Observers & estimators | Residual-based drift detector: has the plant moved away from the model? |
| [`ParameterEstimator`](inc/damp/concepts.hpp#L142) | concept | Core, configuration & backend vocabulary | Online grey-box parameter estimator: physical parameters + gating |
| [`park_transform`](inc/damp/math/transforms.hpp#L331) | function | Scalar math, complex & frames | Park transform (αβ → dq) |
| [`peaking`](inc/damp/filters/iir_design.hpp#L919) | function | Filters & signal conditioning | Peaking (bell) EQ filter: boost or cut a band around f0 |
| [`Periodic`](inc/damp/toolbox/timing.hpp#L167) | block | Embedded helpers (controls-adjacent utilities) | Periodic trigger — fires once per elapsed period |
| [`phase_margin_from_damping_ratio`](inc/damp/design/pid_design.hpp#L844) | function | Design-time synthesis (not PWM-rate) | Approximate phase margin from damping ratio |
| [`phase_margin_unwrapped`](inc/damp/analysis/frequency.hpp#L168) | function | Frequency-domain analysis (host) | Find phase margin using unwrapped phase trajectory |
| [`pi_pole_placement_first_order`](inc/damp/design/pid_design.hpp#L1050) | function | Design-time synthesis (not PWM-rate) | PI gains that place the closed-loop poles of a first-order plant (+1 more overload) |
| [`pid`](inc/damp/controllers/pid.hpp#L275) | function | Design-time synthesis (not PWM-rate) | 2-DOF continuous PID controller design |
| [`pid`](inc/damp/matlab.hpp#L173) | function | MATLAB®-style aliases (host) | MATLAB®-style parallel-form continuous PID constructor |
| [`pid_from_bandwidth`](inc/damp/design/pid_design.hpp#L437) | function | Design-time synthesis (not PWM-rate) | Design PID from desired bandwidth and phase margin |
| [`pid_from_performance_spec`](inc/damp/design/pid_design.hpp#L901) | function | Design-time synthesis (not PWM-rate) | Design PID directly from settling-time and overshoot targets |
| [`pid_pole_placement`](inc/damp/design/pid_design.hpp#L932) | function | Design-time synthesis (not PWM-rate) | Direct discrete PI pole placement for a first-order plant (+1 more overload) |
| [`pid_pole_placement_double_integrator`](inc/damp/design/pid_design.hpp#L1132) | function | Design-time synthesis (not PWM-rate) | PID (or PD) pole placement for a rigid inertia G(s) = 1/(J s²) |
| [`PIDController`](inc/damp/controllers/pid.hpp#L429) | block | Runtime controllers | Fixed-rate discrete 2-DOF PID (canonical runtime) (+3 more overloads) |
| [`PIDMode`](inc/damp/controllers/pid.hpp#L300) | enum | Runtime controllers | Compile-time selection of the PID control-law structure |
| [`PIDPerformanceSpec`](inc/damp/design/pid_design.hpp#L883) | block | Design-time synthesis (not PWM-rate) | Time-domain performance targets for quick PID synthesis |
| [`PIDResult`](inc/damp/controllers/pid.hpp#L100) | block | Design-time synthesis (not PWM-rate) | 2-DOF continuous-time PID controller design result |
| [`PIDRuntimeMode`](inc/damp/controllers/pid.hpp#L325) | enum | Runtime controllers | Runtime operating mode for PIDController / ContinuousPID |
| [`pidstd`](inc/damp/matlab.hpp#L201) | function | MATLAB®-style aliases (host) | Standard-form continuous PID constructor (1-DOF) |
| [`pidstd2`](inc/damp/matlab.hpp#L219) | function | MATLAB®-style aliases (host) | Standard-form continuous 2-DOF PID constructor |
| [`pidtune`](inc/damp/design/pid_design.hpp#L750) | function | Design-time synthesis (not PWM-rate) | Plant-aware PID tune at a crossover (+3 more overloads) |
| [`pidtune`](inc/damp/matlab.hpp#L817) | function | MATLAB®-style aliases (host) | Plant-aware PID tune (MATLAB® pidtune spelling) |
| [`PidTuneSpec`](inc/damp/design/pid_design.hpp#L515) | block | Design-time synthesis (not PWM-rate) | Knobs for plant-aware pidtune |
| [`PIDType`](inc/damp/design/pid_design.hpp#L48) | enum | Design-time synthesis (not PWM-rate) | PID controller type selection for tuning methods |
| [`PiecewiseLTI`](inc/damp/simulation/hybrid.hpp#L94) | block | Simulation / SIL harness (host) | Plant: NModes continuous LTI pieces (same state/input/output dimensions) |
| [`pinv`](inc/damp/matlab.hpp#L461) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the Moore–Penrose pseudoinverse |
| [`place`](inc/damp/design/pole_placement.hpp#L88) | function | Design-time synthesis (not PWM-rate) | Robust multi-input pole placement (Kautsky–Nichols–Van Dooren, real poles) (+5 more overloads) |
| [`place`](inc/damp/matlab.hpp#L590) | function | MATLAB®-style aliases (host) | Robust multi-input pole placement (MATLAB®'s place) |
| [`place_discrete`](inc/damp/design/pole_placement.hpp#L595) | function | Design-time synthesis (not PWM-rate) | Discrete plant + z-plane poles → discrete \(K\) (+3 more overloads) |
| [`place_jordan`](inc/damp/design/pole_placement.hpp#L1152) | function | Design-time synthesis (not PWM-rate) | Exact pole placement with an arbitrary Jordan structure (Schmid–Ntogramatzidis–Nguyen–Pandey / Klein–Moore parametric form) |
| [`place_jordan_optimal`](inc/damp/design/pole_placement.hpp#L1251) | function | Design-time synthesis (not PWM-rate) | Robust / minimum-gain arbitrary pole placement (Schmid et al., Methods 1–2) |
| [`place_observer`](inc/damp/design/pole_placement.hpp#L673) | function | Design-time synthesis (not PWM-rate) | Continuous observer gain: place eigenvalues of \(A - LC\) (s-plane poles) (+5 more overloads) |
| [`place_observer_discrete`](inc/damp/design/pole_placement.hpp#L783) | function | Design-time synthesis (not PWM-rate) | Discrete plant + z-plane poles → discrete observer gain \(L\) (+3 more overloads) |
| [`plot_bode`](inc/damp/simulation/plot_plotly.hpp#L422) | function | Simulation / SIL harness (host) | Plot Bode magnitude and phase as subplots |
| [`plot_line`](inc/damp/simulation/plot_plotly.hpp#L491) | function | Simulation / SIL harness (host) | Simple line plot of time vs value |
| [`plot_simulation`](inc/damp/simulation/plot_plotly.hpp#L343) | function | Simulation / SIL harness (host) | Plot simulation results with subplots for states, outputs, and inputs |
| [`plot_step`](inc/damp/simulation/plot_plotly.hpp#L621) | function | Simulation / SIL harness (host) | Plot step response data |
| [`plot_xy`](inc/damp/simulation/plot_plotly.hpp#L556) | function | Simulation / SIL harness (host) | Multi-series planar scatter (true 2-D path / locus) (+1 more overload) |
| [`plots_root_for`](inc/damp/simulation/plot_plotly.hpp#L118) | function | Simulation / SIL harness (host) | Walk up from dir until a directory named "plots" is found |
| [`pole`](inc/damp/matlab.hpp#L868) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the open-loop poles of a system |
| [`PoleInfo`](inc/damp/analysis/poles.hpp#L75) | block | Frequency-domain analysis (host) | Natural frequency and damping ratio for each pole |
| [`poles`](inc/damp/analysis/poles.hpp#L41) | function | Frequency-domain analysis (host) | Compute open-loop poles (eigenvalues of A matrix) |
| [`PoleZeroMap`](inc/damp/analysis/poles.hpp#L121) | block | Frequency-domain analysis (host) | Poles and zeros of a system, for pole-zero plotting |
| [`poly_horner`](inc/damp/toolbox/scaling.hpp#L134) | function | Embedded helpers (controls-adjacent utilities) | Evaluate a polynomial at x by Horner's method |
| [`poly_roots`](inc/damp/analysis/poles.hpp#L169) | function | Frequency-domain analysis (host) | Roots of a polynomial given in ascending powers (MATLAB® `roots`, reversed order) |
| [`poly_roots`](inc/damp/systems/zpk.hpp#L119) | function | LTI systems (SS / TF / ZPK / discretize) | Compute roots of an ascending-power real polynomial |
| [`PolyRootsResult`](inc/damp/systems/zpk.hpp#L91) | block | LTI systems (SS / TF / ZPK / discretize) | Roots of a real polynomial given in ascending powers |
| [`Pose`](inc/damp/kinematics/pose.hpp#L60) | block | Kinematics / pose | Rigid-body pose: a translation and an orientation (unit quaternion) |
| [`positive_sequence_ab`](inc/damp/filters/pll.hpp#L223) | function | Filters & signal conditioning | Instantaneous positive-sequence αβ from a quadrature signal pair |
| [`pow`](inc/damp/matrix/functions.hpp#L886) | function | Linear algebra | Integer matrix power via binary exponentiation (+1 more overload) |
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
| [`pseudo_inverse`](inc/damp/matrix/svd.hpp#L334) | function | Linear algebra | Moore–Penrose pseudoinverse A⁺ via SVD |
| [`PulseTimer`](inc/damp/toolbox/logic.hpp#L171) | block | Embedded helpers (controls-adjacent utilities) | Pulse timer (non-retriggerable): a rising edge of in emits a fixed |
| [`PwmScheme`](inc/damp/motor/modulation.hpp#L79) | enum | Motor control pack (if present) | Carrier-based three-phase VSI modulation scheme |
| [`pzmap`](inc/damp/analysis/poles.hpp#L188) | function | Frequency-domain analysis (host) | Pole-zero map of a SISO transfer function (MATLAB® `pzmap(tf)`) (+2 more overloads) |
| [`pzplot`](inc/damp/simulation/plot_plotly.hpp#L848) | function | Simulation / SIL harness (host) | Plot a pole-zero map on the complex plane (poles as ×, zeros as ○) |
| [`QPResult`](inc/damp/design/qp.hpp#L82) | block | Design-time synthesis (not PWM-rate) | Result of a dense QP solve |
| [`QPStatus`](inc/damp/design/qp.hpp#L51) | enum | Design-time synthesis (not PWM-rate) | Termination status of a QP solve |
| [`qr_decompose`](inc/damp/matrix/decomposition.hpp#L230) | function | Linear algebra | Thin QR via modified Gram–Schmidt |
| [`QRDecomposition`](inc/damp/matrix/decomposition.hpp#L194) | block | Linear algebra | Thin QR factorization A = QR (modified Gram–Schmidt) |
| [`QuadMode`](inc/damp/toolbox/encoder.hpp#L47) | enum | Embedded helpers (controls-adjacent utilities) | Quadrature decode resolution (edges counted per A/B cycle) |
| [`quadprog`](inc/damp/matlab.hpp#L1088) | function | MATLAB®-style aliases (host) | MATLAB® alias for the dense inequality-constrained QP solve (+1 more overload) |
| [`quadratic_form`](inc/damp/matrix/core.hpp#L1113) | function | Linear algebra | Symmetric congruence (quadratic) form  S = M X Mᵀ |
| [`QuadratureDecoder`](inc/damp/toolbox/encoder.hpp#L75) | block | Embedded helpers (controls-adjacent utilities) | Software A/B quadrature decoder with optional index |
| [`Quaternion`](inc/damp/math/geometry.hpp#L360) | block | Scalar math, complex & frames | Unit quaternion rotation (w, x, y, z) (Hamilton product) |
| [`R_TRIG`](inc/damp/toolbox/iec61131.hpp#L119) | block | Embedded helpers (controls-adjacent utilities) | R_TRIG (Rising Edge Trigger) |
| [`rad2deg`](inc/damp/math/math.hpp#L482) | function | Scalar math, complex & frames | Radians to degrees, rad·180/π |
| [`raise_loop_for_tracking`](inc/damp/design/pid_design.hpp#L663) | function | Design-time synthesis (not PWM-rate) | Raise the loop crossover so T −3 dB matches the 1-DOF target |
| [`ramp`](inc/damp/estimation/excitation/ramp.hpp#L94) | function | Observers & estimators | Build a ramp design payload from a configuration |
| [`Ramp`](inc/damp/estimation/excitation/ramp.hpp#L110) | block | Observers & estimators | Rate-limited ramp runtime generator |
| [`RampConfig`](inc/damp/estimation/excitation/ramp.hpp#L32) | block | Observers & estimators | Configuration for a slew-rate-limited ramp excitation |
| [`RampResult`](inc/damp/estimation/excitation/ramp.hpp#L64) | block | Observers & estimators | Ramp design payload |
| [`RangeMonitor`](inc/damp/toolbox/conditioning.hpp#L440) | block | Embedded helpers (controls-adjacent utilities) | Analog-input range/fault monitor (NAMUR NE43 pattern) |
| [`rank`](inc/damp/design/stability.hpp#L206) | function | Design-time synthesis (not PWM-rate) | Numerical rank of a matrix (forwards to mat::rank) (+1 more overload) |
| [`rank`](inc/damp/matrix/functions.hpp#L233) | function | Linear algebra | Numerical rank: the number of singular values above rel_tol·σₘₐₓ |
| [`rank_from_svd`](inc/damp/matrix/svd.hpp#L302) | function | Linear algebra | Numerical rank from a precomputed SVD result |
| [`reduced_luenberger`](inc/damp/estimation/luenberger.hpp#L271) | function | Observers & estimators | Design a reduced-order (Gopinath) observer by pole placement (matrix form) (+3 more overloads) |
| [`ReducedLuenberger`](inc/damp/estimation/luenberger.hpp#L507) | block | Observers & estimators | Reduced-order (Gopinath) state observer (runtime) |
| [`ReducedLuenbergerResult`](inc/damp/estimation/luenberger.hpp#L217) | block | Observers & estimators | Reduced-order (Gopinath) observer design result |
| [`RedundantAnalogPair`](inc/damp/toolbox/io.hpp#L90) | block | Embedded helpers (controls-adjacent utilities) | Dual-channel analog input with same-slope or opposite-slope cross-check |
| [`reg`](inc/damp/matlab.hpp#L538) | function | MATLAB®-style aliases (host) | Form the observer-based compensator y → u from a plant, K and L |
| [`relay_autotune`](inc/damp/estimation/relay_autotune.hpp#L214) | function | Observers & estimators | Build a validated relay-autotune design payload |
| [`RelayAutotuneConfig`](inc/damp/estimation/relay_autotune.hpp#L121) | block | Observers & estimators | Configuration for the relay-feedback autotuning experiment |
| [`RelayAutotuneOutput`](inc/damp/estimation/relay_autotune.hpp#L239) | block | Observers & estimators | Per-tick output of RelayAutotuner::step |
| [`RelayAutotuner`](inc/damp/estimation/relay_autotune.hpp#L258) | block | Observers & estimators | Runtime relay-feedback autotuner |
| [`RelayAutotuneResult`](inc/damp/estimation/relay_autotune.hpp#L181) | block | Observers & estimators | Relay-autotuner design payload |
| [`RelayAutotuneStatus`](inc/damp/estimation/relay_autotune.hpp#L223) | enum | Observers & estimators | Lifecycle state of a RelayAutotuner |
| [`repetitive`](inc/damp/controllers/repetitive.hpp#L199) | function | Design-time synthesis (not PWM-rate) | Synthesize a repetitive controller with a scalar robustness filter Q |
| [`repetitive_binomial`](inc/damp/controllers/repetitive.hpp#L235) | function | Design-time synthesis (not PWM-rate) | Synthesize a repetitive controller with a binomial zero-phase FIR Q |
| [`RepetitiveConfig`](inc/damp/controllers/repetitive.hpp#L71) | block | Runtime controllers | Repetitive-controller tuning + period (with optional zero-phase FIR Q) |
| [`RepetitiveController`](inc/damp/controllers/repetitive.hpp#L280) | block | Runtime controllers | Plug-in repetitive controller runtime (fixed-size internal model) |
| [`RepetitiveResult`](inc/damp/controllers/repetitive.hpp#L137) | block | Design-time synthesis (not PWM-rate) | Design result for the repetitive controller |
| [`rescale`](inc/damp/toolbox/scaling.hpp#L69) | function | Embedded helpers (controls-adjacent utilities) | Affine map of x from the input range to the output range |
| [`resonance_autotune_from_frf`](inc/damp/estimation/commissioning_frontends.hpp#L76) | function | Observers & estimators | Map FRF peaks and valleys to a fixed bank of biquad compensators |
| [`ResonanceAutotuneResult`](inc/damp/estimation/commissioning_frontends.hpp#L44) | block | Observers & estimators | Result of resonance / anti-resonance compensator design from an FRF |
| [`ResonantMode`](inc/damp/estimation/frequency_response.hpp#L384) | block | Observers & estimators | One extracted FRF extremum (peak or valley) |
| [`RisingEdge`](inc/damp/toolbox/logic.hpp#L34) | block | Embedded helpers (controls-adjacent utilities) | Rising-edge detector: true on the tick x goes false → true |
| [`RK23`](inc/damp/simulation/integrator.hpp#L859) | block | Simulation / SIL harness (host) | Bogacki-Shampine 2(3) adaptive integrator |
| [`RK3`](inc/damp/simulation/integrator.hpp#L926) | block | Simulation / SIL harness (host) | Classical 3rd-order Runge-Kutta (RK3) integrator |
| [`RK4`](inc/damp/simulation/integrator.hpp#L651) | block | Simulation / SIL harness (host) | Classical 4th-order Runge-Kutta (RK4) integrator |
| [`rlocus`](inc/damp/analysis/poles.hpp#L296) | function | Frequency-domain analysis (host) | Root locus of a SISO state-space plant over an explicit gain grid (+1 more overload) |
| [`rlocusplot`](inc/damp/simulation/plot_plotly.hpp#L1033) | function | Simulation / SIL harness (host) | Plot a root locus (closed-loop poles vs gain) on the complex plane |
| [`rls`](inc/damp/estimation/rls.hpp#L174) | function | Observers & estimators | Build scalar RLS design payload |
| [`Rls`](inc/damp/estimation/rls.hpp#L213) | block | Observers & estimators | Scalar runtime RLS estimator |
| [`rls_vector`](inc/damp/estimation/rls.hpp#L190) | function | Observers & estimators | Build vector RLS design payload |
| [`RlsConfig`](inc/damp/estimation/rls.hpp#L60) | block | Observers & estimators | Common RLS configuration |
| [`RlsResult`](inc/damp/estimation/rls.hpp#L102) | block | Observers & estimators | Scalar RLS design payload |
| [`RlsState`](inc/damp/estimation/rls.hpp#L88) | block | Observers & estimators | Scalar RLS runtime state |
| [`RlsVector`](inc/damp/estimation/rls.hpp#L313) | block | Observers & estimators | Vector runtime RLS estimator (NP parameters) |
| [`RlsVectorResult`](inc/damp/estimation/rls.hpp#L145) | block | Observers & estimators | Vector RLS design payload for N parameters |
| [`RlsVectorState`](inc/damp/estimation/rls.hpp#L131) | block | Observers & estimators | Vector RLS runtime state for N parameters |
| [`RobustExactDifferentiator`](inc/damp/filters/differentiator.hpp#L63) | block | Filters & signal conditioning | First-order robust exact differentiator (super-twisting differentiator) |
| [`RootLocusResult`](inc/damp/analysis/poles.hpp#L237) | block | Frequency-domain analysis (host) | Root-locus data: closed-loop poles along a gain grid |
| [`rotary_gearbox`](inc/damp/toolbox/actuator.hpp#L149) | function | Embedded helpers (controls-adjacent utilities) | Build a ServoAxis for a rotary joint behind a gearbox |
| [`rows`](inc/damp/estimation/ins_eskf.hpp#L322) | function | Observers & estimators | Sparse y = F x for the 15-state INS first-order structure (G = I path) |
| [`RowVec`](inc/damp/matrix/rowvec.hpp#L29) | block | Linear algebra | Row vector specialization of Matrix<1, N, T> |
| [`RowView`](inc/damp/matrix/views.hpp#L189) | block | Linear algebra | Non-owning row view of a matrix |
| [`RS`](inc/damp/toolbox/iec61131.hpp#L92) | block | Embedded helpers (controls-adjacent utilities) | RS Latch (Reset-Set Latch) |
| [`scaled_deadband`](inc/damp/toolbox/conditioning.hpp#L243) | function | Embedded helpers (controls-adjacent utilities) | Center dead zone that rescales the surviving range back to full span |
| [`schroeder_multi_sine`](inc/damp/estimation/excitation/multi_sine.hpp#L180) | function | Observers & estimators | Build a multi-sine design payload with Schroeder phases applied |
| [`SecondOrderCoeffs`](inc/damp/filters/iir_design.hpp#L65) | block | Filters & signal conditioning | DSP coefficients for second-order IIR filter |
| [`SEL`](inc/damp/toolbox/iec61131.hpp#L601) | function | Embedded helpers (controls-adjacent utilities) | SEL (IEC 61131-3 binary selection): g ? in1 : in0 |
| [`select`](inc/damp/systems/state_space.hpp#L544) | function | LTI systems (SS / TF / ZPK / discretize) | SISO channel slice y_iy / u_iu |
| [`SequenceComponents`](inc/damp/math/transforms.hpp#L490) | block | Scalar math, complex & frames | Symmetrical (sequence) components of a three-phase phasor set |
| [`series`](inc/damp/systems/state_space.hpp#L204) | function | LTI systems (SS / TF / ZPK / discretize) | Series connection: sys2 follows sys1 (u → sys1 → sys2 → y) (+1 more overload) |
| [`ServoAxis`](inc/damp/toolbox/actuator.hpp#L107) | block | Embedded helpers (controls-adjacent utilities) | One servoactuator transmission: SI joint unit ⟷ drive (motor) units |
| [`ServoBank`](inc/damp/toolbox/actuator.hpp#L229) | block | Embedded helpers (controls-adjacent utilities) | A bank of ServoAxis transmissions: maps a synchronized multi-axis TrajectoryState array straight to drive commands in one call |
| [`ServoCommand`](inc/damp/toolbox/actuator.hpp#L84) | block | Embedded helpers (controls-adjacent utilities) | A drive-native servoactuator setpoint: position, velocity, torque |
| [`sgn`](inc/damp/math/math.hpp#L398) | function | Scalar math, complex & frames | Sign function — −1 if val < 0, 1 if val > 0, 0 if val == 0 |
| [`sigma`](inc/damp/analysis/frequency.hpp#L809) | function | Frequency-domain analysis (host) | Singular-value frequency response of a (possibly MIMO) state-space system (+1 more overload) |
| [`sigmaplot`](inc/damp/simulation/plot_plotly.hpp#L990) | function | Simulation / SIL harness (host) | Plot singular-value frequency response (log frequency, dB) |
| [`SigmaPoint`](inc/damp/analysis/frequency.hpp#L752) | block | Frequency-domain analysis (host) | Singular values of G(jω) at one frequency |
| [`SigmaResult`](inc/damp/analysis/frequency.hpp#L766) | block | Frequency-domain analysis (host) | Singular-value frequency response over a grid |
| [`SignalSource`](inc/damp/concepts.hpp#L127) | concept | Core, configuration & backend vocabulary | Self-clocked signal source: u = step(), finished when done() |
| [`SignalStatus`](inc/damp/toolbox/conditioning.hpp#L378) | enum | Embedded helpers (controls-adjacent utilities) | Classification of an analog input against its valid/fault bands |
| [`simc`](inc/damp/design/pid_design.hpp#L343) | function | Design-time synthesis (not PWM-rate) | SIMC (Skogestad Internal Model Control) tuning for FOPDT models |
| [`simulate`](inc/damp/simulation/simulate.hpp#L133) | function | Simulation / SIL harness (host) | Simulate a nonlinear plant with a controller in closed loop |
| [`simulate_cached_zoh`](inc/damp/simulation/cached_zoh.hpp#L148) | function | Simulation / SIL harness (host) | Simulate a continuous LTI plant with fixed ZOH step and held input (+1 more overload) |
| [`simulate_cached_zoh_multirate`](inc/damp/simulation/cached_zoh.hpp#L227) | function | Simulation / SIL harness (host) | Multi-rate ZOH: held control samples, dense plant samples for plotting |
| [`simulate_discrete`](inc/damp/simulation/simulate.hpp#L489) | function | Simulation / SIL harness (host) | Simulate a discrete-time system with a controller |
| [`simulate_discrete_nonlinear`](inc/damp/simulation/simulate.hpp#L442) | function | Simulation / SIL harness (host) | Simulate a discrete-time nonlinear plant with a controller |
| [`simulate_lti`](inc/damp/simulation/simulate.hpp#L323) | function | Simulation / SIL harness (host) | Simulate a continuous LTI system with a controller |
| [`simulate_lti_siso`](inc/damp/simulation/simulate.hpp#L385) | function | Simulation / SIL harness (host) | SIL: continuous SISO LTI plant under a sampled SISO controller with constant reference |
| [`simulate_multirate`](inc/damp/simulation/multirate.hpp#L128) | function | Simulation / SIL harness (host) | Multi-rate closed-loop simulation with MultiRateConfig (+1 more overload) |
| [`simulate_piecewise_lti`](inc/damp/simulation/hybrid.hpp#L183) | function | Simulation / SIL harness (host) | Simulate a piecewise-LTI plant with Exact integration between events |
| [`simulate_sampled`](inc/damp/simulation/simulate.hpp#L206) | function | Simulation / SIL harness (host) | Simulate a continuous plant under a discrete (sampled) controller — multi-rate |
| [`simulate_state_feedback`](inc/damp/simulation/simulate.hpp#L268) | function | Simulation / SIL harness (host) | Simulate a nonlinear plant with state-feedback controller |
| [`simulate_two_rate`](inc/damp/simulation/multirate.hpp#L245) | function | Simulation / SIL harness (host) | Two-rate cascade with MultiRateConfig (+1 more overload) |
| [`simulate_two_rate_lti`](inc/damp/simulation/multirate.hpp#L390) | function | Simulation / SIL harness (host) | Two-rate cascade on a continuous LTI plant (A, B, C) |
| [`SimulationResult`](inc/damp/simulation/simulate.hpp#L104) | block | Simulation / SIL harness (host) | Result of a closed-loop simulation |
| [`sin`](inc/damp/matrix/functions.hpp#L1047) | function | Linear algebra | Matrix sine via scaling and double-angle reconstruction |
| [`sin`](inc/damp/math/math.hpp#L182) | function | Scalar math, complex & frames | Sine |
| [`sincos`](inc/damp/matrix/functions.hpp#L977) | function | Linear algebra | Compute sin(A) and cos(A) together via scaling and double-angle reconstruction |
| [`sincos`](inc/damp/math/math.hpp#L200) | function | Scalar math, complex & frames | Combined sine and cosine, {sin(x), cos(x)} |
| [`SinglePhasePLL`](inc/damp/filters/pll.hpp#L48) | block | Filters & signal conditioning | Single-Phase PLL |
| [`sinh`](inc/damp/matrix/functions.hpp#L1117) | function | Linear algebra | Matrix hyperbolic sine sinh(A) = (exp(A) − exp(−A))/2 |
| [`siso_ref`](inc/damp/simulation/simulate.hpp#L91) | function | Simulation / SIL harness (host) | Build a SisoReferenceAdapter for a SISO `control(r,y)` controller |
| [`SisoController`](inc/damp/concepts.hpp#L74) | concept | Core, configuration & backend vocabulary | Scalar output-feedback controller: u = control(r, y), fixed rate |
| [`SisoReferenceAdapter`](inc/damp/simulation/simulate.hpp#L72) | block | Simulation / SIL harness (host) | Adapt a SISO `control(r, y)` controller for callables that want `u = f(y)` or `u = f(t, y)` (the older simulate contracts) |
| [`six_step_duty_cycles`](inc/damp/motor/modulation.hpp#L466) | function | Motor control pack (if present) | Classical six-step (full-wave) duty pattern from the αβ angle |
| [`six_step_fundamental_voltage`](inc/damp/motor/modulation.hpp#L146) | function | Motor control pack (if present) | Six-step (square-wave) fundamental peak phase voltage |
| [`SlewLimiter`](inc/damp/toolbox/conditioning.hpp#L278) | block | Embedded helpers (controls-adjacent utilities) | Slew-rate limiter: bound how fast the output may follow the target |
| [`smc`](inc/damp/controllers/smc.hpp#L79) | function | Design-time synthesis (not PWM-rate) | Bundle hand-picked SMC parameters into an SMCResult |
| [`SMCController`](inc/damp/controllers/smc.hpp#L160) | block | Runtime controllers | First-order sliding-mode controller (SMC) for a SISO plant |
| [`SMCResult`](inc/damp/controllers/smc.hpp#L33) | block | Design-time synthesis (not PWM-rate) | Tuning parameters for a first-order sliding-mode controller |
| [`smith_predictor_from_fopdt`](inc/damp/controllers/smith_predictor.hpp#L181) | function | Design-time synthesis (not PWM-rate) | Design a Smith predictor from FOPDT parameters |
| [`smith_predictor_from_pid`](inc/damp/controllers/smith_predictor.hpp#L247) | function | Design-time synthesis (not PWM-rate) | Pack a Smith predictor from an existing continuous PID and FOPDT plant |
| [`SmithPredictor`](inc/damp/controllers/smith_predictor.hpp#L307) | block | Runtime controllers | Smith predictor runtime: SISO primary + FO model + pure delay |
| [`SmithPredictorResult`](inc/damp/controllers/smith_predictor.hpp#L85) | block | Design-time synthesis (not PWM-rate) | Smith predictor design result (discrete PID + FO model + delay samples) |
| [`soft_sign`](inc/damp/toolbox/conditioning.hpp#L105) | function | Embedded helpers (controls-adjacent utilities) | Smooth unit direction x / √(x²+ε²) ∈ (−1,1) |
| [`sogi`](inc/damp/filters/sogi.hpp#L39) | function | Filters & signal conditioning | Second-Order Generalized Integrator (SOGI) design (+1 more overload) |
| [`SOGI`](inc/damp/filters/sogi.hpp#L187) | block | Filters & signal conditioning | Runtime SOGI quadrature-signal generator (live-frequency resonator) |
| [`SogiFll`](inc/damp/filters/sogi.hpp#L340) | block | Filters & signal conditioning | SOGI with a Frequency-Locked Loop — self-tuning single-tone tracker |
| [`solve`](inc/damp/matrix/solve.hpp#L82) | function | Linear algebra | Solve lower-triangular system LX = B via forward substitution (+2 more overloads) |
| [`solve_miqp`](inc/damp/design/qp.hpp#L1245) | function | Design-time synthesis (not PWM-rate) | Mixed-integer QP by branch and bound over the active-set relaxation |
| [`solve_qp`](inc/damp/design/qp.hpp#L546) | function | Design-time synthesis (not PWM-rate) | Solve a strictly convex inequality-constrained QP (dual active-set) (+1 more overload) |
| [`solve_qp_admm`](inc/damp/design/qp.hpp#L937) | function | Design-time synthesis (not PWM-rate) | One-shot ADMM QP solve (cold start) |
| [`solve_qp_interior_point`](inc/damp/design/qp.hpp#L981) | function | Design-time synthesis (not PWM-rate) | Interior-point QP solver (primal-dual path following) |
| [`SolveResult`](inc/damp/simulation/solver.hpp#L86) | block | Simulation / SIL harness (host) | Result of an ODE solve operation |
| [`specific_force_at_rest`](inc/damp/estimation/ins_mechanization.hpp#L227) | function | Observers & estimators | Specific force [m/s²] a stationary IMU measures for the given orientation |
| [`SplineSurface`](inc/damp/toolbox/lookup.hpp#L394) | block | Embedded helpers (controls-adjacent utilities) | 2-D interpolating surface with smooth (Catmull-Rom spline) blending |
| [`sqrt`](inc/damp/matrix/functions.hpp#L749) | function | Linear algebra | Matrix square root via Denman–Beavers iteration |
| [`sqrt`](inc/damp/math/complex.hpp#L337) | function | Scalar math, complex & frames | Compute complex square root (constexpr) (+1 more overload) |
| [`sqrtm`](inc/damp/matrix/functions.hpp#L793) | function | Linear algebra | @brief MATLAB®-style alias for sqrt |
| [`SR`](inc/damp/toolbox/iec61131.hpp#L64) | block | Embedded helpers (controls-adjacent utilities) | SR Latch (Set-dominant Set-Reset Latch) |
| [`ss`](inc/damp/matlab.hpp#L103) | function | MATLAB®-style aliases (host) | MATLAB®-style state-space model constructor |
| [`ss2tf`](inc/damp/systems/zpk.hpp#L858) | function | LTI systems (SS / TF / ZPK / discretize) | MATLAB®-style alias for SISO to_transfer_function |
| [`ss2zpk`](inc/damp/systems/zpk.hpp#L828) | function | LTI systems (SS / TF / ZPK / discretize) | MATLAB®-style alias for SISO to_zpk (state-space) |
| [`ss_gain`](inc/damp/systems/state_space.hpp#L518) | function | LTI systems (SS / TF / ZPK / discretize) | Static gain y = k u (no states) |
| [`ss_integrator`](inc/damp/systems/state_space.hpp#L528) | function | LTI systems (SS / TF / ZPK / discretize) | Integrator y = (k/s) u |
| [`SsFeed`](inc/damp/systems/state_space.hpp#L507) | block | LTI systems (SS / TF / ZPK / discretize) | One term of a linear output mix: k · y_iy |
| [`stability_margin_continuous`](inc/damp/design/stability.hpp#L311) | function | Design-time synthesis (not PWM-rate) | Compute stability margin for continuous system |
| [`stability_margin_discrete`](inc/damp/design/stability.hpp#L340) | function | Design-time synthesis (not PWM-rate) | Compute stability margin for discrete system |
| [`state_mpc`](inc/damp/controllers/mpc.hpp#L510) | function | Design-time synthesis (not PWM-rate) | Synthesize a constrained linear MPC (condensed dense QP, Δu form) (+1 more overload) |
| [`StateEstimator`](inc/damp/concepts.hpp#L114) | concept | Core, configuration & backend vocabulary | State estimator: x̂ = estimate(y, u) — the fused per-tick form |
| [`StateFeedback`](inc/damp/controllers/lqr.hpp#L445) | block | Runtime controllers | Runtime full-state feedback law u = −Kx |
| [`StateFeedbackController`](inc/damp/concepts.hpp#L101) | concept | Core, configuration & backend vocabulary | State-feedback law: u = control(r, x), reference-first |
| [`StateJacobian`](inc/damp/estimation/ekf.hpp#L40) | block | Observers & estimators | State prediction result from the user's dynamics function |
| [`StateSpace`](inc/damp/systems/state_space.hpp#L128) | block | LTI systems (SS / TF / ZPK / discretize) | State-space representation for linear time-invariant systems (discrete or continuous) |
| [`StateSpaceZPKResult`](inc/damp/systems/zpk.hpp#L592) | block | LTI systems (SS / TF / ZPK / discretize) | SISO state-space → ZPK conversion result with runtime zero count |
| [`steady_state_vdq`](inc/damp/motor/spm.hpp#L79) | function | Motor control pack (if present) | Steady-state Vdq (rotor frame) |
| [`steady_state_voltage_magnitude`](inc/damp/motor/spm.hpp#L92) | function | Motor control pack (if present) | \|Vdq\| at a steady-state operating point |
| [`SteadyStateKalmanFilter`](inc/damp/estimation/kalman.hpp#L383) | block | Observers & estimators | Steady-state (fixed-gain) Kalman estimator for LQG-class designs |
| [`steer_map`](inc/damp/toolbox/io.hpp#L452) | function | Embedded helpers (controls-adjacent utilities) | Curve-driven steering map — arbitrary `(throttle, turn)` geometry |
| [`steinhart_hart`](inc/damp/toolbox/thermistor.hpp#L122) | function | Embedded helpers (controls-adjacent utilities) | Fit the Steinhart-Hart coefficients from three calibration points |
| [`step`](inc/damp/analysis/time_response.hpp#L160) | function | Frequency-domain analysis (host) | Step response of a (MIMO) state-space system (+1 more overload) |
| [`step`](inc/damp/estimation/excitation/step.hpp#L91) | function | Observers & estimators | Build a single-step design payload |
| [`Step`](inc/damp/estimation/excitation/step.hpp#L108) | block | Observers & estimators | Single step with optional \|dy/dt\| settle detection |
| [`step_train`](inc/damp/estimation/excitation/step_train.hpp#L94) | function | Observers & estimators | Build a step-train design payload from a configuration |
| [`StepConfig`](inc/damp/estimation/excitation/step.hpp#L35) | block | Observers & estimators | Configuration for a single step with optional settle detection |
| [`StepInfo`](inc/damp/analysis/time_response.hpp#L377) | block | Frequency-domain analysis (host) | Step-response characteristics of a single output signal |
| [`stepinfo`](inc/damp/analysis/time_response.hpp#L398) | function | Frequency-domain analysis (host) | Compute step-response characteristics from an output/time signal (+1 more overload) |
| [`stepped_sine`](inc/damp/estimation/excitation/stepped_sine.hpp#L104) | function | Observers & estimators | Build a stepped-sine design payload from a configuration |
| [`SteppedSine`](inc/damp/estimation/excitation/stepped_sine.hpp#L125) | block | Observers & estimators | Stepped-sine excitation — one pure tone at a time across a frequency table |
| [`SteppedSineConfig`](inc/damp/estimation/excitation/stepped_sine.hpp#L33) | block | Observers & estimators | Configuration for single-tone stepped-sine excitation |
| [`SteppedSineResult`](inc/damp/estimation/excitation/stepped_sine.hpp#L73) | block | Observers & estimators | Stepped-sine design payload |
| [`stepplot`](inc/damp/simulation/plot_plotly.hpp#L717) | function | Simulation / SIL harness (host) | Plot a step response, one trace per input/output pair |
| [`StepResult`](inc/damp/estimation/excitation/step.hpp#L67) | block | Observers & estimators | Single-step design payload |
| [`StepTrain`](inc/damp/estimation/excitation/step_train.hpp#L110) | block | Observers & estimators | Alternating +/- step train runtime generator |
| [`StepTrainConfig`](inc/damp/estimation/excitation/step_train.hpp#L32) | block | Observers & estimators | Configuration for alternating +/- step excitation |
| [`StepTrainResult`](inc/damp/estimation/excitation/step_train.hpp#L64) | block | Observers & estimators | Step-train design payload |
| [`Stopwatch`](inc/damp/toolbox/timing.hpp#L90) | block | Embedded helpers (controls-adjacent utilities) | Free-running elapsed-time accumulator |
| [`stsmc`](inc/damp/controllers/stsmc.hpp#L87) | function | Design-time synthesis (not PWM-rate) | Synthesize super-twisting gains from a disturbance-derivative bound |
| [`stsmc_gains`](inc/damp/controllers/stsmc.hpp#L123) | function | Design-time synthesis (not PWM-rate) | Super-twisting controller from gains you specify directly |
| [`STSMCController`](inc/damp/controllers/stsmc.hpp#L180) | block | Runtime controllers | Super-twisting controller (second-order sliding mode) |
| [`STSMCResult`](inc/damp/controllers/stsmc.hpp#L30) | block | Design-time synthesis (not PWM-rate) | Super-twisting (second-order sliding-mode) controller design result |
| [`subtract`](inc/damp/systems/state_space.hpp#L405) | function | LTI systems (SS / TF / ZPK / discretize) | Differencing connection: outputs y = y₁ − y₂ |
| [`SuccessiveCompensatorCommissioner`](inc/damp/estimation/successive_compensator.hpp#L208) | block | Observers & estimators | Successive biquad compensator bank (notch and/or band-pass) |
| [`SuccessiveCompensatorConfig`](inc/damp/estimation/successive_compensator.hpp#L50) | block | Observers & estimators | Configuration for successive compensator commissioning |
| [`suggest_mpc_horizon`](inc/damp/controllers/mpc.hpp#L593) | function | Design-time synthesis (not PWM-rate) | Suggest MPC horizons from an explicit settling-time target (+1 more overload) |
| [`summarize_loop_response`](inc/damp/analysis/frequency.hpp#L473) | function | Frequency-domain analysis (host) | Summarize loop_response() results into one compact metrics struct |
| [`svd`](inc/damp/matrix/svd.hpp#L277) | function | Linear algebra | Full singular value decomposition A = U·Σ·Vᴴ (one-sided Jacobi) |
| [`svd`](inc/damp/matlab.hpp#L452) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the singular value decomposition |
| [`SVDResult`](inc/damp/matrix/svd.hpp#L259) | block | Linear algebra | Result of a full singular value decomposition A = U·Σ·Vᴴ |
| [`svm_duty_cycles`](inc/damp/motor/modulation.hpp#L443) | function | Motor control pack (if present) | Space-vector PWM duty cycles from an αβ voltage command |
| [`SvmDuties`](inc/damp/motor/modulation.hpp#L103) | block | Motor control pack (if present) | Result of a duty-map: half-bridge duties plus an over-modulation flag |
| [`svpwm_zero_sequence`](inc/damp/motor/modulation.hpp#L200) | function | Motor control pack (if present) | Min-max zero-sequence injection for space-vector PWM |
| [`Switch`](inc/damp/toolbox/io.hpp#L282) | block | Embedded helpers (controls-adjacent utilities) | Debounced maintained switch (toggle/selector contact) with change flag |
| [`SwitchedController`](inc/damp/controllers/composition.hpp#L117) | block | Runtime controllers | Bumpless-ish switch between normal, experiment, and backup SISO laws |
| [`SwitchMode`](inc/damp/controllers/composition.hpp#L94) | enum | Runtime controllers | Which path owns the plant command |
| [`symmetrical_components`](inc/damp/math/transforms.hpp#L520) | function | Scalar math, complex & frames | Forward symmetrical-component (Fortescue) transform (abc → 012) |
| [`SymplecticEuler`](inc/damp/simulation/integrator.hpp#L250) | block | Simulation / SIL harness (host) | Semi-implicit (symplectic) Euler for mechanical systems |
| [`Tachometer`](inc/damp/toolbox/encoder.hpp#L153) | block | Embedded helpers (controls-adjacent utilities) | Pulse-based speed (tachometer) with frequency/period crossover |
| [`tan`](inc/damp/math/math.hpp#L215) | function | Scalar math, complex & frames | Tangent |
| [`tf`](inc/damp/matlab.hpp#L65) | function | MATLAB®-style aliases (host) | MATLAB®-style transfer function constructor (+1 more overload) |
| [`tf2zpk`](inc/damp/systems/zpk.hpp#L467) | function | LTI systems (SS / TF / ZPK / discretize) | MATLAB®-style alias for to_zpk (transfer function) |
| [`TFF`](inc/damp/toolbox/iec61131.hpp#L507) | block | Embedded helpers (controls-adjacent utilities) | T Flip-Flop (toggle on rising edge) |
| [`thermal_kalman_plant`](inc/damp/toolbox/thermal.hpp#L188) | function | Embedded helpers (controls-adjacent utilities) | Prepare a discretized thermal plant for design::kalman |
| [`ThermalKalmanObserver`](inc/damp/toolbox/thermal.hpp#L241) | block | Embedded helpers (controls-adjacent utilities) | Runtime thermal-network Kalman observer |
| [`Thermistor`](inc/damp/toolbox/thermistor.hpp#L177) | block | Embedded helpers (controls-adjacent utilities) | NTC thermistor linearization (resistance → temperature) |
| [`ThermistorCoeffs`](inc/damp/toolbox/thermistor.hpp#L41) | block | Embedded helpers (controls-adjacent utilities) | Fitted NTC coefficients in Steinhart-Hart form |
| [`thipwm_zero_sequence`](inc/damp/motor/modulation.hpp#L229) | function | Motor control pack (if present) | Third-harmonic injection (THIPWM) zero-sequence — 1/6 of fundamental |
| [`ThreePhasePLL`](inc/damp/filters/pll.hpp#L141) | block | Filters & signal conditioning | Synchronous-reference-frame (SRF) PLL for balanced three-phase input |
| [`Timeout`](inc/damp/toolbox/timing.hpp#L125) | block | Embedded helpers (controls-adjacent utilities) | One-shot timeout |
| [`TimeResponse`](inc/damp/analysis/time_response.hpp#L51) | block | Frequency-domain analysis (host) | Multi-channel time-domain response sampled on a time grid |
| [`to_coeffs`](inc/damp/filters/iir_design.hpp#L649) | function | Filters & signal conditioning | Convert StateSpace system to first-order DSP coefficients (+4 more overloads) |
| [`to_double_vector`](inc/damp/simulation/plot_plotly.hpp#L236) | function | Simulation / SIL harness (host) | Convert std::vector\<T\> to std::vector<double> |
| [`to_std_vector`](inc/damp/simulation/plot_plotly.hpp#L211) | function | Simulation / SIL harness (host) | Convert a ColVec<N,T> to std::vector<double> for plotlypp |
| [`to_transfer_function`](inc/damp/systems/zpk.hpp#L841) | function | LTI systems (SS / TF / ZPK / discretize) | SISO state-space → transfer function (Leverrier), size NX+1 / NX+1 |
| [`to_zpk`](inc/damp/systems/zpk.hpp#L432) | function | LTI systems (SS / TF / ZPK / discretize) | Convert a SISO transfer function to zero-pole-gain form (+1 more overload) |
| [`TOF`](inc/damp/toolbox/iec61131.hpp#L224) | block | Embedded helpers (controls-adjacent utilities) | TOF Timer (Timer Off Delay) |
| [`Toggle`](inc/damp/toolbox/logic.hpp#L244) | block | Embedded helpers (controls-adjacent utilities) | Toggle (T flip-flop): output flips on each rising edge of in |
| [`TON`](inc/damp/toolbox/iec61131.hpp#L184) | block | Embedded helpers (controls-adjacent utilities) | TON Timer (Timer On Delay) |
| [`Tone`](inc/damp/estimation/excitation/multi_sine.hpp#L29) | block | Observers & estimators | One sinusoidal component in a multi-sine excitation |
| [`torque`](inc/damp/motor/spm.hpp#L103) | function | Motor control pack (if present) | Electromagnetic torque of an SPM (i_d = 0): T_e = Kₜ i_q |
| [`torque_constant_from_flux`](inc/damp/motor/spm.hpp#L29) | function | Motor control pack (if present) | Kₜ = 1½ p λ [Nm/A] (amplitude-invariant) |
| [`torque_constant_from_Kv`](inc/damp/motor/spm.hpp#L46) | function | Motor control pack (if present) | Kₜ from hobby Kᵥ [RPM/V] (peak line-to-line / bus-volt sense) |
| [`torque_speed_envelope`](inc/damp/motor/spm.hpp#L203) | function | Motor control pack (if present) | Fixed table of SPM torque–speed envelope samples (host / LUT bake) |
| [`TorqueSpeedPoint`](inc/damp/motor/spm.hpp#L144) | block | Motor control pack (if present) | One sample of the SPM torque–speed envelope (MTPA i_d=0 branch) |
| [`TP`](inc/damp/toolbox/iec61131.hpp#L262) | block | Embedded helpers (controls-adjacent utilities) | TP Timer (Timer Pulse) |
| [`TrajectoryBoundary`](inc/damp/trajectory/trajectory_types.hpp#L178) | block | Trajectory value types | Boundary conditions at one endpoint of a polynomial trajectory: a position and its time derivatives through jerk |
| [`TrajectoryLimits`](inc/damp/trajectory/trajectory_types.hpp#L45) | block | Trajectory value types | Asymmetric kinematic limits for a trapezoidal or S-curve motion profile |
| [`TrajectoryState`](inc/damp/trajectory/trajectory_types.hpp#L76) | block | Trajectory value types | A point on a motion profile: commanded position, velocity, acceleration |
| [`TransferFunction`](inc/damp/systems/transfer_function.hpp#L74) | block | LTI systems (SS / TF / ZPK / discretize) | SISO polynomial transfer function G(s) = num(s)/den(s) |
| [`Transform4`](inc/damp/math/geometry.hpp#L802) | block | Scalar math, complex & frames | 4×4 homogeneous transform (SE(3) interop / DH export) |
| [`Translation3`](inc/damp/kinematics/pose.hpp#L37) | block | Kinematics / pose | A 3-D translation — a thin Vec3 with domain-named conveniences |
| [`TransposeView`](inc/damp/matrix/views.hpp#L396) | block | Linear algebra | Non-owning transpose view of a matrix (zero-copy) |
| [`Trapezoidal`](inc/damp/simulation/integrator.hpp#L594) | block | Simulation / SIL harness (host) | Trapezoidal (Tustin) integrator |
| [`TRBDF2`](inc/damp/simulation/integrator.hpp#L485) | block | Simulation / SIL harness (host) | TR-BDF2 composite integrator — the stiff adaptive pair (ode23tb) |
| [`two_norm`](inc/damp/matrix/functions.hpp#L132) | function | Linear algebra | Spectral norm ‖A‖₂ = σₘₐₓ(A) |
| [`two_point_cal`](inc/damp/toolbox/scaling.hpp#L111) | function | Embedded helpers (controls-adjacent utilities) | Fit an AffineCal through two `(raw, engineering)` points |
| [`TwoRateSimulationResult`](inc/damp/simulation/multirate.hpp#L94) | block | Simulation / SIL harness (host) | Result of a two-rate cascade simulation |
| [`tyreus_luyben`](inc/damp/design/pid_design.hpp#L178) | function | Design-time synthesis (not PWM-rate) | Tyreus-Luyben tuning from ultimate gain and ultimate period |
| [`UKFMeasFn`](inc/damp/estimation/ukf.hpp#L64) | concept | Observers & estimators | Concept for UKF measurement functions |
| [`UKFStateFn`](inc/damp/estimation/ukf.hpp#L53) | concept | Observers & estimators | Concept for UKF state (process) functions |
| [`unbounded_bound`](inc/damp/design/qp.hpp#L67) | function | Design-time synthesis (not PWM-rate) | Sentinel bound treated as "no constraint" on that row |
| [`UnscentedKalmanFilter`](inc/damp/estimation/ukf.hpp#L127) | block | Observers & estimators | Unscented (sigma-point) Kalman Filter for nonlinear discrete-time systems |
| [`UnscentedParams`](inc/damp/estimation/ukf.hpp#L80) | block | Observers & estimators | Tuning parameters for the scaled unscented transform |
| [`unwrap_phase_deg`](inc/damp/analysis/frequency.hpp#L127) | function | Frequency-domain analysis (host) | Unwrap phase data in degrees to avoid +/-180 discontinuities |
| [`UpperTriangle`](inc/damp/matrix/views.hpp#L122) | block | Linear algebra | Upper triangular view of a square matrix |
| [`voltage_circle_radius`](inc/damp/motor/foc.hpp#L57) | function | Motor control pack (if present) | SVPWM voltage-circle radius Vₘₐₓ = m · Vdc / √3 |
| [`WarmStartActiveSetSolver`](inc/damp/design/qp.hpp#L675) | block | Design-time synthesis (not PWM-rate) | Warm-started active-set solver policy — the damp::MPC default |
| [`with_schroeder_phases`](inc/damp/estimation/excitation/multi_sine.hpp#L160) | function | Observers & estimators | Assign Schroeder low-crest-factor phases to a multi-sine tone table |
| [`wrap`](inc/damp/math/math.hpp#L510) | function | Scalar math, complex & frames | Wrap x into the half-open interval [min, max) (period max − min) |
| [`wrap_pi`](inc/damp/math/math.hpp#L545) | function | Scalar math, complex & frames | Wrap an angle to [−π, π) |
| [`wrap_two_pi`](inc/damp/math/math.hpp#L563) | function | Scalar math, complex & frames | Wrap an angle to [0, 2π) |
| [`wrapped_delta`](inc/damp/toolbox/encoder.hpp#L41) | function | Embedded helpers (controls-adjacent utilities) | Signed difference between two unsigned counter readings, wrap-safe |
| [`write_html`](inc/damp/simulation/plot_plotly.hpp#L145) | function | Simulation / SIL harness (host) | Write a plotlypp figure to HTML with damp's shell (prefer over writeHtml) |
| [`write_hybrid_trace_npy`](inc/damp/simulation/npy_export.hpp#L100) | function | Simulation / SIL harness (host) | Pack HybridSimulationResult into (N, 2+NX+NU+NY) float64: t, mode, x..., u..., y |
| [`write_npy_f64`](inc/damp/simulation/npy_export.hpp#L41) | function | Simulation / SIL harness (host) | Write a C-order float64 array as a .npy file |
| [`XyPlotOpts`](inc/damp/simulation/plot_plotly.hpp#L535) | block | Simulation / SIL harness (host) | Optional axis framing for plot_xy |
| [`XySeries`](inc/damp/simulation/plot_plotly.hpp#L517) | block | Simulation / SIL harness (host) | One series for a planar (x, y) scatter / path plot |
| [`zero_sequence`](inc/damp/math/transforms.hpp#L226) | function | Scalar math, complex & frames | Zero-sequence (common-mode) scalar of a three-phase set |
| [`ziegler_nichols`](inc/damp/design/pid_design.hpp#L76) | function | Design-time synthesis (not PWM-rate) | Ziegler-Nichols tuning from ultimate gain and ultimate period |
| [`ziegler_nichols_step`](inc/damp/design/pid_design.hpp#L127) | function | Design-time synthesis (not PWM-rate) | Ziegler-Nichols step response method (reaction curve) |
| [`ZPK`](inc/damp/systems/zpk.hpp#L175) | block | LTI systems (SS / TF / ZPK / discretize) | Zero-pole-gain (ZPK) representation of a SISO LTI system |
| [`zpk`](inc/damp/matlab.hpp#L126) | function | MATLAB®-style aliases (host) | MATLAB®-style zero-pole-gain model constructor |
| [`zpk2tf`](inc/damp/systems/zpk.hpp#L476) | function | LTI systems (SS / TF / ZPK / discretize) | MATLAB®-style alias for ZPK::to_transfer_function |
| [`ZPKResult`](inc/damp/systems/zpk.hpp#L402) | block | LTI systems (SS / TF / ZPK / discretize) | Result of converting a transfer function to zero-pole-gain form |
