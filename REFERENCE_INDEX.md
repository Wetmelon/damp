# API Reference — Alphabetical Index

Auto-generated from `@brief` doc comments in `inc/damp/`. Regenerate with `python tools/gen_reference.py`. Grouped-by-domain view: [REFERENCE.md](REFERENCE.md).

| Name | Kind | Domain | Description |
| ---- | ---- | ------ | ----------- |
| [`abs`](inc/damp/math/complex.hpp#L362) | function | Scalar math, complex & frames | Compute magnitude (absolute value) of a complex number (+1 more overload) |
| [`accumulate_poly`](inc/damp/systems/transfer_function.hpp#L38) | function | LTI systems (SS / TF / ZPK / discretize) | dst += scale·src, aligned on the constant term (the last entry: coefficients are highest power first) |
| [`acker`](inc/damp/matlab.hpp#L849) | function | MATLAB®-style aliases (host) | Single-input pole placement (MATLAB®'s acker) |
| [`acos`](inc/damp/math/math.hpp#L152) | function | Scalar math, complex & frames | Arccosine ∈ [0, π]. Input is clamped to [−1, 1] in both paths |
| [`ActionGovernor`](inc/damp/controllers/action_governor.hpp#L237) | block | Runtime controllers | Runtime action governor — QP projection of u_des onto A u ≤ b |
| [`ActiveSetSolver`](inc/damp/design/qp.hpp#L703) | block | Design-time synthesis (not PWM-rate) | Default QP solver policy: the Goldfarb–Idnani active-set solve_qp() |
| [`AdaptationGate`](inc/damp/estimation/parameter_estimation.hpp#L339) | block | Observers & estimators | Gated-adaptation policy: decide whether a parameter update is safe to apply |
| [`adaptive_solve`](inc/damp/simulation/solver.hpp#L814) | function | Simulation / SIL harness (host) | Adaptive-step solve in one call |
| [`AdaptiveLut1D`](inc/damp/toolbox/lookup.hpp#L428) | block | Embedded helpers (controls-adjacent utilities) | Adaptive 1-D lookup: fixed breakpoints, online-updated cell values |
| [`AdaptiveLut2D`](inc/damp/toolbox/lookup.hpp#L501) | block | Embedded helpers (controls-adjacent utilities) | Adaptive 2-D lookup: fixed grid, online-updated cells (nearest grid node) |
| [`AdaptiveOptions`](inc/damp/simulation/solver.hpp#L349) | block | Simulation / SIL harness (host) | Options for adaptive-step ODE integration |
| [`AdaptiveStepIntegrator`](inc/damp/simulation/integrator.hpp#L115) | concept | Simulation / SIL harness (host) | Integrator that reports a genuine embedded local-error estimate |
| [`AdaptiveStepSolver`](inc/damp/simulation/solver.hpp#L404) | block | Simulation / SIL harness (host) | Adaptive-step ODE solver (+1 more overload) |
| [`AdaptOutOfRange`](inc/damp/toolbox/lookup.hpp#L48) | enum | Embedded helpers (controls-adjacent utilities) | How AdaptiveLut1D / AdaptiveLut2D treat samples outside the breakpoint span |
| [`AdmmSettings`](inc/damp/design/qp.hpp#L804) | block | Design-time synthesis (not PWM-rate) | Tuning parameters for the ADMM QP solver |
| [`AdmmSolver`](inc/damp/design/qp.hpp#L885) | block | Design-time synthesis (not PWM-rate) | ADMM (OSQP-style) QP solver policy — warm-started, fixed-cost iterations |
| [`adrc`](inc/damp/controllers/adrc.hpp#L206) | function | Design-time synthesis (not PWM-rate) | Active Disturbance Rejection Control design |
| [`ADRCController`](inc/damp/controllers/adrc.hpp#L384) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its design, keeping the result's scalar |
| [`ADRCController`](inc/damp/controllers/adrc.hpp#L257) | block | Runtime controllers | Active Disturbance Rejection Control (ADRC) |
| [`ADRCResult`](inc/damp/controllers/adrc.hpp#L51) | block | Design-time synthesis (not PWM-rate) | Active Disturbance Rejection Control design result |
| [`AffineCal`](inc/damp/toolbox/scaling.hpp#L83) | block | Embedded helpers (controls-adjacent utilities) | Affine sensor calibration `y = gain·x + offset` |
| [`all_of`](inc/damp/matrix/core.hpp#L1453) | function | Linear algebra | Whether p holds for every entry of A: the matrix counterpart of std::all_of |
| [`allmargin`](inc/damp/matlab.hpp#L1314) | function | MATLAB®-style aliases (host) | Gain, phase, and delay margins of a SISO loop over a frequency grid |
| [`AllMarginResult`](inc/damp/matlab.hpp#L1293) | block | MATLAB®-style aliases (host) | All classical margins including delay margin (superset of margin) |
| [`allpass`](inc/damp/filters/iir_design.hpp#L718) | function | Filters & signal conditioning | Second-order all-pass filter design (continuous-time) (+1 more overload) |
| [`AlphaBeta`](inc/damp/math/transforms.hpp#L164) | block | Scalar math, complex & frames | Alpha-beta (stationary-frame) component pair |
| [`amigo_kappa_tau`](inc/damp/design/pid_design.hpp#L253) | function | Design-time synthesis (not PWM-rate) | AMIGO PI from ultimate gain/period and static gain Kₛ |
| [`AnalogCrossMode`](inc/damp/toolbox/conditioning.hpp#L467) | enum | Embedded helpers (controls-adjacent utilities) | How two redundant analog channels should relate in engineering units |
| [`AnalogCrossResult`](inc/damp/toolbox/conditioning.hpp#L489) | block | Embedded helpers (controls-adjacent utilities) | Output of cross_check_analog |
| [`AnalogCrossStatus`](inc/damp/toolbox/conditioning.hpp#L475) | enum | Embedded helpers (controls-adjacent utilities) | Result of comparing two redundant analog channels |
| [`AnalogInput`](inc/damp/toolbox/io.hpp#L47) | block | Embedded helpers (controls-adjacent utilities) | A single analog input: range/fault check on the raw reading, then affine calibration to engineering units |
| [`AnalyticDynamics`](inc/damp/concepts.hpp#L233) | concept | Core, configuration & backend vocabulary | Dynamics that provide their own Jacobians: f.jacobian(x, u) → DynamicsJacobian |
| [`AnalyticMeasurement`](inc/damp/concepts.hpp#L239) | concept | Core, configuration & backend vocabulary | A Measurement that provides its own Jacobians: h.jacobian(x, u) → MeasurementJacobian |
| [`any_of`](inc/damp/matrix/core.hpp#L1477) | function | Linear algebra | Whether p holds for some entry of A: the matrix counterpart of std::any_of |
| [`apply_wet_plotly_shell`](inc/damp/simulation/plot_plotly.hpp#L66) | function | Simulation / SIL harness (host) | Apply damp's full-window HTML shell to plotlypp-generated markup |
| [`arange`](inc/damp/analysis/linspace.hpp#L115) | function | Frequency-domain analysis (host) | Half-open arithmetic range [start, stop) with step step (+2 more overloads) |
| [`arg`](inc/damp/math/complex.hpp#L373) | function | Scalar math, complex & frames | Compute argument (phase angle) of a complex number |
| [`array_as`](inc/damp/backend.hpp#L195) | function | Core, configuration & backend vocabulary | Element-wise static_cast of a damp::array to scalar type U (precision rebind) |
| [`asin`](inc/damp/math/math.hpp#L133) | function | Scalar math, complex & frames | Arcsine ∈ [−π/2, π/2]. Input is clamped to [−1, 1] in both paths so behavior matches at compile and run time (std::asin would return NaN for \|x\| > 1) |
| [`atan`](inc/damp/math/math.hpp#L119) | function | Scalar math, complex & frames | Single-argument arctangent ∈ (−π/2, π/2) |
| [`atan2`](inc/damp/math/math.hpp#L107) | function | Scalar math, complex & frames | Two-argument arctangent, atan2(y, x) ∈ [−π, π] |
| [`AxisInput`](inc/damp/toolbox/io.hpp#L148) | block | Embedded helpers (controls-adjacent utilities) | Operator-axis conditioning chain (joystick / RC stick → command) |
| [`AxisPosition`](inc/damp/encoder/axis_position.hpp#L47) | block | Embedded helpers (controls-adjacent utilities) | The position of an axis, rotary or linear: an exact count carrying its own scale |
| [`backward_substitute_transpose`](inc/damp/matrix/solve.hpp#L61) | function | Linear algebra | Backward substitution for Lᴴx = b (Lᵀx = b for real T) |
| [`BackwardEuler`](inc/damp/simulation/integrator.hpp#L343) | block | Simulation / SIL harness (host) | Backward Euler integrator |
| [`balance_states`](inc/damp/systems/realization.hpp#L82) | function | LTI systems (SS / TF / ZPK / discretize) | Balance a model by a diagonal state scaling, so A, B and C have comparable magnitudes |
| [`balanced_realization`](inc/damp/design/model_reduction.hpp#L478) | function | Design-time synthesis (not PWM-rate) | Balanced realization via the square-root (Moore / Laub) method |
| [`balanced_reduction`](inc/damp/design/model_reduction.hpp#L555) | function | Design-time synthesis (not PWM-rate) | Balanced reduction to NR states, by truncation or DC-matched residualization |
| [`BalancedRealizationResult`](inc/damp/design/model_reduction.hpp#L77) | block | Design-time synthesis (not PWM-rate) | Balanced realization result (Moore square-root method) |
| [`BalancedReductionResult`](inc/damp/design/model_reduction.hpp#L109) | block | Design-time synthesis (not PWM-rate) | Reduced-order model from balanced truncation or residualization |
| [`balreal`](inc/damp/matlab.hpp#L343) | function | MATLAB®-style aliases (host) | MATLAB® short alias for design::balanced_realization |
| [`balred`](inc/damp/matlab.hpp#L362) | function | MATLAB®-style aliases (host) | MATLAB® short alias for design::balanced_reduction to NR states |
| [`bandpass`](inc/damp/filters/iir_design.hpp#L609) | function | Filters & signal conditioning | Second-order band-pass filter design (continuous-time) (+1 more overload) |
| [`bandwidth`](inc/damp/matlab.hpp#L1402) | function | MATLAB®-style aliases (host) | -3 dB bandwidth of a SISO system over a frequency grid |
| [`BDF2`](inc/damp/simulation/integrator.hpp#L403) | block | Simulation / SIL harness (host) | Backward Differentiation Formula 2 (BDF2) integrator |
| [`Biquad`](inc/damp/filters/biquad.hpp#L538) | block | Filters & signal conditioning | Second-order IIR (biquad) section runtime (+1 more overload) |
| [`BiquadCascade`](inc/damp/filters/biquad.hpp#L628) | block | Filters & signal conditioning | Cascade of second-order sections (SOS) for higher-order IIR filters (+1 more overload) |
| [`BLINK`](inc/damp/toolbox/iec61131.hpp#L542) | block | Embedded helpers (controls-adjacent utilities) | BLINK (free-running square-wave / flasher) |
| [`blkdiag`](inc/damp/matrix/concat.hpp#L135) | function | Linear algebra | Block-diagonal assembly of MatrixLike arguments |
| [`Block`](inc/damp/matrix/block.hpp#L38) | block | Linear algebra | Block view (non-owning) into a parent matrix |
| [`bode`](inc/damp/analysis/frequency.hpp#L267) | function | Frequency-domain analysis (host) | Compute Bode plot data for a SISO state-space system (+2 more overloads) |
| [`bodemag`](inc/damp/simulation/plot_plotly.hpp#L798) | function | Simulation / SIL harness (host) | Plot a magnitude-only Bode diagram (log frequency, dB magnitude) |
| [`bodeplot`](inc/damp/simulation/plot_plotly.hpp#L784) | function | Simulation / SIL harness (host) | Plot magnitude and phase Bode subplots |
| [`BodeResult`](inc/damp/analysis/frequency.hpp#L59) | block | Frequency-domain analysis (host) | Bode plot data for a SISO system |
| [`Bounds`](inc/damp/toolbox/bounds.hpp#L52) | block | Embedded helpers (controls-adjacent utilities) | A per-channel closed-interval box constraint |
| [`butterworth_lowpass`](inc/damp/filters/iir_design.hpp#L776) | function | Filters & signal conditioning | Butterworth low-pass filter design (+1 more overload) |
| [`Button`](inc/damp/toolbox/io.hpp#L195) | block | Embedded helpers (controls-adjacent utilities) | Debounced momentary push-button with edge and long-press detection |
| [`c2d`](inc/damp/matlab.hpp#L425) | function | MATLAB®-style aliases (host) | MATLAB® interface function c2d to discretize a continuous-time state-space system (+1 more overload) |
| [`CachedZoh`](inc/damp/simulation/cached_zoh.hpp#L63) | block | Simulation / SIL harness (host) | Cached discrete ZOH maps Ad, Bd for step size h |
| [`canon`](inc/damp/matlab.hpp#L500) | function | MATLAB®-style aliases (host) | MATLAB® canon(sys) / canon(sys, 'modal'): real modal canonical form and its transform (+1 more overload) |
| [`canonical_phase_margin`](inc/damp/analysis/frequency.hpp#L150) | function | Frequency-domain analysis (host) | Normalize phase margin to (-180, 180] |
| [`care`](inc/damp/design/riccati.hpp#L1491) | function | Design-time synthesis (not PWM-rate) | Solve the Continuous-time Algebraic Riccati Equation (CARE) |
| [`Cascade`](inc/damp/controllers/composition.hpp#L50) | block | Runtime controllers | Series cascade of two SISO controllers: outer → inner reference |
| [`catmull_1d`](inc/damp/toolbox/lookup.hpp#L340) | function | Embedded helpers (controls-adjacent utilities) | 1-D non-uniform Catmull-Rom (cubic Hermite) evaluation |
| [`cbf_relative_degree_1`](inc/damp/controllers/action_governor.hpp#L205) | function | Design-time synthesis (not PWM-rate) | Build a relative-degree-1 CBF inequality row |
| [`CBFConstraint`](inc/damp/controllers/action_governor.hpp#L189) | block | Runtime controllers | One affine row from a relative-degree-1 CBF condition |
| [`cbrt`](inc/damp/math/math.hpp#L94) | function | Scalar math, complex & frames | Cube root (preserves sign for negative x) |
| [`ceil`](inc/damp/math/math.hpp#L383) | function | Scalar math, complex & frames | Ceiling — smallest integer ≥ x |
| [`chirp`](inc/damp/estimation/excitation/chirp.hpp#L123) | function | Observers & estimators | Build a chirp design payload from a configuration |
| [`Chirp`](inc/damp/estimation/excitation/chirp.hpp#L141) | block | Observers & estimators | Linear or logarithmic chirp runtime generator (+1 more overload) |
| [`ChirpConfig`](inc/damp/estimation/excitation/chirp.hpp#L41) | block | Observers & estimators | Configuration for a sine chirp excitation |
| [`ChirpMode`](inc/damp/estimation/excitation/chirp.hpp#L27) | enum | Observers & estimators | Chirp sweep law |
| [`ChirpResult`](inc/damp/estimation/excitation/chirp.hpp#L84) | block | Observers & estimators | Chirp design payload |
| [`cholesky`](inc/damp/matrix/decomposition.hpp#L124) | function | Linear algebra | Cholesky decomposition for positive-definite matrices |
| [`cholesky_solve`](inc/damp/matrix/solve.hpp#L159) | function | Linear algebra | Solve AX = B via Cholesky (A = LLᴴ) |
| [`clarke_park_transform`](inc/damp/math/transforms.hpp#L404) | function | Scalar math, complex & frames | Fused Clarke-Park transform (abc → dq) |
| [`clarke_transform`](inc/damp/math/transforms.hpp#L277) | function | Scalar math, complex & frames | Clarke transform (abc → αβ) |
| [`classical_dob`](inc/damp/estimation/dob.hpp#L404) | function | Observers & estimators | Synthesize a classical disturbance observer from a nominal plant and Q-filter |
| [`ClassicalDOB`](inc/damp/estimation/dob.hpp#L439) | block | Observers & estimators | Classical Pn^-1·Q disturbance observer runtime (bolt-on compensator) (+1 more overload) |
| [`ClassicalDobResult`](inc/damp/estimation/dob.hpp#L339) | block | Observers & estimators | Design result for the classical Pn^-1·Q disturbance observer |
| [`classify_range`](inc/damp/toolbox/conditioning.hpp#L369) | function | Embedded helpers (controls-adjacent utilities) | Classify x against the valid span and the wire-fault limits of range |
| [`closed_loop_poles`](inc/damp/design/stability.hpp#L564) | function | Design-time synthesis (not PWM-rate) | Compute closed-loop poles (eigenvalues) with state feedback |
| [`closed_loop_stable`](inc/damp/design/pid_design.hpp#L684) | function | Design-time synthesis (not PWM-rate) | Closed-loop certificate: den_P·den_C + num_P·num_C is Hurwitz for C_fb = Kp + Kᵢ/s + Kd s/(1 + Tf s) |
| [`cohen_coon`](inc/damp/design/pid_design.hpp#L300) | function | Design-time synthesis (not PWM-rate) | Cohen-Coon tuning from first-order-plus-dead-time model |
| [`ColVec`](inc/damp/matrix/colvec.hpp#L30) | block | Linear algebra | Concrete Column vector specialization of Matrix<N, 1, T> |
| [`ColView`](inc/damp/matrix/views.hpp#L332) | block | Linear algebra | Non-owning column view of a matrix |
| [`comb_notch_window`](inc/damp/filters/moving_average.hpp#L34) | function | Filters & signal conditioning | Window length for a moving-average comb that notches f_notch and all its harmonics: N = round(1 / (f_notch·Ts)) |
| [`combine_quadrature`](inc/damp/encoder/sincos.hpp#L120) | function | Embedded helpers (controls-adjacent utilities) | Absolute position from a quadrature edge counter and the fine phase |
| [`CommandProjectionResult`](inc/damp/controllers/action_governor.hpp#L91) | block | Design-time synthesis (not PWM-rate) | Result of an affine command projection |
| [`companion_form`](inc/damp/systems/realization.hpp#L617) | function | LTI systems (SS / TF / ZPK / discretize) | Companion form of a SISO model, MATLAB®'s canon(sys, 'companion'): B = e₁, −α in A's last column |
| [`Complementary`](inc/damp/filters/complementary.hpp#L32) | block | Filters & signal conditioning | Scalar (1-D) complementary filter — fuse a fast rate with a slow absolute |
| [`ComplementaryFilter`](inc/damp/estimation/sensor_fusion.hpp#L61) | block | Observers & estimators | Simple complementary filter for orientation estimation |
| [`complex`](inc/damp/math/complex.hpp#L38) | block | Scalar math, complex & frames | Constexpr complex number class for compile-time computations |
| [`compute_eigenvalues`](inc/damp/matrix/eigen.hpp#L763) | function | Linear algebra | Compute the eigenvalues (and Schur vectors) of a real square matrix |
| [`constant_power_torque_limit`](inc/damp/motor/foc.hpp#L44) | function | Motor control pack (if present) | Constant-power torque ceiling Tₘₐₓ = P_rated / \|ω\| |
| [`ConstantInertiaFeedforward`](inc/damp/toolbox/actuator.hpp#L201) | block | Embedded helpers (controls-adjacent utilities) | Per-axis decoupled torque feedforward: `τ = J·a + b·v + τ_c·sign(v) + g` |
| [`ConstMatrix`](inc/damp/matrix/const_matrix.hpp#L88) | block | Linear algebra | A matrix that is a constant of its type: no storage, and a product with a vector that has no term for an entry that is exactly 0 |
| [`continuous_lpf_exact_step`](inc/damp/filters/lowpass.hpp#L39) | function | Filters & signal conditioning | Exact step of continuous first-order LPF ẏ = −ω_c (y − u) |
| [`continuous_lqr`](inc/damp/design/riccati.hpp#L1744) | function | Design-time synthesis (not PWM-rate) | Continuous-time Linear-Quadratic Regulator design |
| [`ContinuousPID`](inc/damp/controllers/pid.hpp#L1071) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its design, keeping the result's scalar |
| [`ContinuousPID`](inc/damp/controllers/pid.hpp#L941) | block | Runtime controllers | Continuous-gain PID with per-tick sample time (variable-rate secondary form) |
| [`controllability_gramian`](inc/damp/design/stability.hpp#L382) | function | Design-time synthesis (not PWM-rate) | Continuous/discrete controllability Gramian W_c |
| [`controllability_matrix`](inc/damp/design/stability.hpp#L55) | function | Design-time synthesis (not PWM-rate) | Compute the controllability matrix [B, AB, A²B, ..., A^(N-1)B] |
| [`controllability_staircase`](inc/damp/systems/realization.hpp#L792) | function | LTI systems (SS / TF / ZPK / discretize) | Controllability staircase form: the controllable part split off by an orthogonal transform |
| [`controllable_canonical_form`](inc/damp/systems/realization.hpp#L519) | function | LTI systems (SS / TF / ZPK / discretize) | Controllable canonical form of a SISO model: ones on A's superdiagonal, −α in its last row, B = eₙ |
| [`Convention`](inc/damp/math/transforms.hpp#L69) | enum | Scalar math, complex & frames | Scaling convention for the Clarke/Park family |
| [`ConvertsExactly`](inc/damp/matrix/matrix_traits.hpp#L80) | concept | Linear algebra | Matrix elements of type From convert to To without rounding: the same type, or a real into the complex of its own precision |
| [`copysign`](inc/damp/math/math.hpp#L484) | function | Scalar math, complex & frames | Copy sign — magnitude of mag with the sign of sgn_src |
| [`cos`](inc/damp/matrix/functions.hpp#L1353) | function | Linear algebra | Matrix cosine via scaling and double-angle reconstruction |
| [`cos`](inc/damp/math/math.hpp#L171) | function | Scalar math, complex & frames | Cosine |
| [`cosh`](inc/damp/matrix/functions.hpp#L1387) | function | Linear algebra | Matrix hyperbolic cosine cosh(A) = (exp(A) + exp(−A))/2 |
| [`Counter`](inc/damp/toolbox/logic.hpp#L273) | block | Embedded helpers (controls-adjacent utilities) | Edge-counting up/down counter: increments on each rising edge of up, decrements on each rising edge of down. Returns the running count, saturating at the limits of T (no wrap / signed overflow) |
| [`cross_check_analog`](inc/damp/toolbox/conditioning.hpp#L543) | function | Embedded helpers (controls-adjacent utilities) | Cross-check two calibrated analog readings for redundant sensors |
| [`CTD`](inc/damp/toolbox/iec61131.hpp#L354) | block | Embedded helpers (controls-adjacent utilities) | CTD Counter (Count Down) |
| [`ctrb`](inc/damp/matlab.hpp#L386) | function | MATLAB®-style aliases (host) | MATLAB® short alias for controllability_matrix (+1 more overload) |
| [`ctrbf`](inc/damp/matlab.hpp#L541) | function | MATLAB®-style aliases (host) | MATLAB® ctrbf(A, B, C): controllability staircase form, its orthogonal transform and rank |
| [`CTU`](inc/damp/toolbox/iec61131.hpp#L314) | block | Embedded helpers (controls-adjacent utilities) | CTU Counter (Count Up) |
| [`CTUD`](inc/damp/toolbox/iec61131.hpp#L394) | block | Embedded helpers (controls-adjacent utilities) | CTUD Counter (Count Up Down) |
| [`current_loop_pi`](inc/damp/design/pid_design.hpp#L1298) | function | Design-time synthesis (not PWM-rate) | PI gains for an R–L current loop 1/(L s + R) at any bandwidth ω |
| [`d2d`](inc/damp/matlab.hpp#L458) | function | MATLAB®-style aliases (host) | MATLAB® d2d: resample a ZOH-discrete model to a new sample time |
| [`damp`](inc/damp/analysis/poles.hpp#L74) | function | Frequency-domain analysis (host) | Compute natural frequency and damping for each pole |
| [`damping_ratio_from_overshoot_percent`](inc/damp/design/pid_design.hpp#L858) | function | Design-time synthesis (not PWM-rate) | Map percent overshoot target to equivalent damping ratio |
| [`dare`](inc/damp/design/riccati.hpp#L1322) | function | Design-time synthesis (not PWM-rate) | Solve the Discrete Algebraic Riccati Equation (DARE) |
| [`DareMethod`](inc/damp/design/riccati.hpp#L1189) | enum | Design-time synthesis (not PWM-rate) | DARE solver selection for dare() |
| [`db2mag`](inc/damp/math/math.hpp#L552) | function | Scalar math, complex & frames | Decibels to magnitude, 10^(db/20) |
| [`dcgain`](inc/damp/analysis/norms.hpp#L42) | function | Frequency-domain analysis (host) | Compute DC gain of a continuous-time system |
| [`DCM`](inc/damp/math/geometry.hpp#L86) | block | Scalar math, complex & frames | Rotation matrix, active body→parent (SO(3) wrapper over Mat3) |
| [`deadband`](inc/damp/toolbox/conditioning.hpp#L49) | function | Embedded helpers (controls-adjacent utilities) | Dead zone over `[lower, upper]`, matching Simulink®'s Dead Zone block (+1 more overload) |
| [`Debounce`](inc/damp/toolbox/logic.hpp#L212) | block | Embedded helpers (controls-adjacent utilities) | Debounce: the output adopts in only after in differs from the current output continuously for stable_time. Rejects contact bounce and brief glitches. (Not an IEC block — the one everyone hand-rolls.) |
| [`decode_sincos`](inc/damp/encoder/sincos.hpp#L92) | function | Embedded helpers (controls-adjacent utilities) | Trim a (sin, cos) pair and decode its phase |
| [`decoupling_feedforward`](inc/damp/motor/foc.hpp#L72) | function | Motor control pack (if present) | Decoupling + back-EMF feedforward at R = 0 (reference currents) |
| [`default_tol`](inc/damp/matrix/matrix_traits.hpp#L125) | function | Linear algebra | Type-appropriate default tolerance for floating-point comparisons. float  ~7 decimal digits  → 1e-6 double ~15 decimal digits → 1e-12 |
| [`deg2rad`](inc/damp/math/math.hpp#L574) | function | Scalar math, complex & frames | Degrees to radians, deg·π/180 |
| [`Delay`](inc/damp/filters/delay.hpp#L27) | block | Filters & signal conditioning | Discrete-time delay buffer |
| [`DesignStatus`](inc/damp/design/status.hpp#L35) | enum | Design-time synthesis (not PWM-rate) | Why a design failed, one value per cause class (not per code branch) |
| [`det`](inc/damp/matrix/functions.hpp#L172) | function | Linear algebra | Matrix determinant det(A) |
| [`DFF`](inc/damp/toolbox/iec61131.hpp#L451) | block | Embedded helpers (controls-adjacent utilities) | D Flip-Flop (edge-triggered data latch) |
| [`diag`](inc/damp/matrix/concat.hpp#L55) | function | Linear algebra | Square diagonal matrix from scalars: diag(10.0, 0.1, 100.0, 10.0) (+4 more overloads) |
| [`Diagonal`](inc/damp/matrix/views.hpp#L31) | block | Linear algebra | Diagonal view of a square matrix |
| [`differential_drive`](inc/damp/toolbox/io.hpp#L262) | function | Embedded helpers (controls-adjacent utilities) | Differential (tank) drive mixer: throttle + turn → left/right |
| [`DirectQuadrature`](inc/damp/math/transforms.hpp#L84) | block | Scalar math, complex & frames | Direct-quadrature (rotor-frame) component pair |
| [`Discrete`](inc/damp/simulation/integrator.hpp#L133) | block | Simulation / SIL harness (host) | Discrete-time integrator (no integration, just one step) |
| [`discrete_lqr`](inc/damp/design/riccati.hpp#L1610) | function | Design-time synthesis (not PWM-rate) | Discrete-time Linear-Quadratic Regulator design (+1 more overload) |
| [`discrete_lyapunov_cholesky`](inc/damp/design/lyapunov.hpp#L510) | function | Design-time synthesis (not PWM-rate) | Cholesky factor of the discrete Lyapunov solution: A X Aᵀ − X + B Bᵀ = 0, X = RᵀR |
| [`DiscretePIDGains`](inc/damp/controllers/pid.hpp#L76) | block | Runtime controllers | Fixed-rate discrete PID coefficients (canonical deploy form) |
| [`DiscretizationMethod`](inc/damp/systems/discretization.hpp#L31) | enum | LTI systems (SS / TF / ZPK / discretize) | Discretization methods for continuous-time state-space systems |
| [`DiscretizationResult`](inc/damp/systems/discretization.hpp#L54) | block | LTI systems (SS / TF / ZPK / discretize) | A discretized or resampled model with its success / status, read like an optional |
| [`discretize`](inc/damp/systems/discretization.hpp#L334) | function | LTI systems (SS / TF / ZPK / discretize) | Discretize a continuous-time state-space system (+1 more overload) |
| [`discretize_lqr_cost`](inc/damp/design/riccati.hpp#L1679) | function | Design-time synthesis (not PWM-rate) | Discretize a continuous LQR cost integral over one sample (Van Loan) |
| [`discretize_process_noise`](inc/damp/estimation/kalman.hpp#L269) | function | Observers & estimators | Discretize continuous process noise over one sample (Van Loan) |
| [`disk_margin`](inc/damp/analysis/disk_margin.hpp#L158) | function | Frequency-domain analysis (host) | Disk margin of a SISO loop L under unity negative feedback |
| [`diskmargin`](inc/damp/matlab.hpp#L1331) | function | MATLAB®-style aliases (host) | Disk margin of a SISO loop (+1 more overload) |
| [`DiskMarginResult`](inc/damp/analysis/disk_margin.hpp#L36) | block | Frequency-domain analysis (host) | Disk margin of a SISO loop: the disk size α and the gain and phase margins it guarantees |
| [`DisturbanceModel`](inc/damp/design/disturbance.hpp#L82) | block | Design-time synthesis (not PWM-rate) | An explicit disturbance model: d enters the plant's state through B_d and the readings through C_d |
| [`DLATCH`](inc/damp/toolbox/iec61131.hpp#L482) | block | Embedded helpers (controls-adjacent utilities) | D Latch (level-sensitive / transparent latch) |
| [`dlqe`](inc/damp/matlab.hpp#L959) | function | MATLAB®-style aliases (host) | Discrete linear-quadratic estimator from matrices (MATLAB®'s dlqe) |
| [`dlqr`](inc/damp/matlab.hpp#L907) | function | MATLAB®-style aliases (host) | Discrete-time Linear-Quadratic Regulator design (MATLAB®'s dlqr) |
| [`dlyap`](inc/damp/design/lyapunov.hpp#L470) | function | Design-time synthesis (not PWM-rate) | Solve the discrete-time Lyapunov (Stein) equation A X Aᵀ − X + Q = 0 |
| [`dlyap`](inc/damp/matlab.hpp#L1422) | function | MATLAB®-style aliases (host) | MATLAB® alias for the discrete Lyapunov solve AXAᵀ−X+Q=0 |
| [`dlyapchol`](inc/damp/matlab.hpp#L1442) | function | MATLAB®-style aliases (host) | MATLAB® alias for the Cholesky factor of the discrete Lyapunov solution (X = RᵀR) |
| [`dob`](inc/damp/estimation/dob.hpp#L143) | function | Observers & estimators | Validate and package DOB configuration into a runtime-ready design result |
| [`DOB`](inc/damp/estimation/dob.hpp#L168) | block | Observers & estimators | Lightweight SISO disturbance observer runtime (+1 more overload) |
| [`DOBConfig`](inc/damp/estimation/dob.hpp#L32) | block | Observers & estimators | Configuration for a first-order disturbance observer |
| [`DP45`](inc/damp/simulation/integrator.hpp#L780) | block | Simulation / SIL harness (host) | Dormand-Prince 5(4) adaptive integrator — the ode45 pair |
| [`dpwm_zero_sequence`](inc/damp/motor/modulation.hpp#L269) | function | Motor control pack (if present) | DPWM zero-sequence for a chosen discontinuous scheme |
| [`DqCommand`](inc/damp/motor/foc.hpp#L86) | block | Motor control pack (if present) | Clamped dq voltage command from FOController::current_controller |
| [`DrivePair`](inc/damp/toolbox/io.hpp#L234) | block | Embedded helpers (controls-adjacent utilities) | Left/right actuator pair (differential/arcade drive output) |
| [`DsogiPll`](inc/damp/filters/pll.hpp#L306) | block | Filters & signal conditioning | Dual-SOGI three-phase positive-sequence PLL (DSOGI-PLL) |
| [`duties_from_phase_voltages`](inc/damp/motor/modulation.hpp#L361) | function | Motor control pack (if present) | Map phase voltages + zero-sequence to clamped half-bridge duties |
| [`Dynamics`](inc/damp/concepts.hpp#L195) | concept | Core, configuration & backend vocabulary | Continuous-time nonlinear dynamics ẋ = f(x, u): a callable returning ẋ as a ColVec<NX, T> |
| [`DynamicsJacobian`](inc/damp/concepts.hpp#L214) | block | Core, configuration & backend vocabulary | Analytic Jacobians of a Dynamics model at (x, u), returned by its jacobian(x, u) |
| [`eig`](inc/damp/matlab.hpp#L669) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the eigenvalues of a square matrix |
| [`EigenResult`](inc/damp/matrix/eigen.hpp#L39) | block | Linear algebra | Eigenvalue computation result |
| [`eliminate_states`](inc/damp/design/model_reduction.hpp#L627) | function | Design-time synthesis (not PWM-rate) | Eliminate the listed states, by DC-matched residualization or truncation |
| [`EncoderGate`](inc/damp/encoder/gate.hpp#L109) | block | Embedded helpers (controls-adjacent utilities) | Validates encoder readings before they reach the tracker |
| [`EncoderGateConfig`](inc/damp/encoder/gate.hpp#L43) | block | Embedded helpers (controls-adjacent utilities) | Configuration for EncoderGate. Positions in counts, as the tracker's |
| [`EncoderTracker`](inc/damp/encoder/tracker.hpp#L70) | block | Embedded helpers (controls-adjacent utilities) | Encoder position (exact) and speed (smoothed) from absolute counts |
| [`EncoderTrackerConfig`](inc/damp/encoder/tracker.hpp#L23) | block | Embedded helpers (controls-adjacent utilities) | Configuration for EncoderTracker |
| [`EncoderVerdict`](inc/damp/encoder/gate.hpp#L31) | enum | Embedded helpers (controls-adjacent utilities) | What EncoderGate made of a reading |
| [`enu_from_ned`](inc/damp/estimation/ins_mechanization.hpp#L116) | function | Observers & estimators | Map a NED vector into ENU (axis permute) |
| [`eps`](inc/damp/matrix/matrix_traits.hpp#L115) | function | Linear algebra | Machine epsilon of T's real scalar (MATLAB® eps). Scale it by problem size and magnitude to build a relative tolerance |
| [`ErrorStateJacobian`](inc/damp/estimation/eskf.hpp#L288) | block | Observers & estimators | ESKF prediction Jacobians F and G (nominal state updated outside) |
| [`ErrorStateKalmanFilter`](inc/damp/estimation/eskf.hpp#L342) | block | Observers & estimators | Runtime error-state KF: tracks δx and P, not the full nominal state (+1 more overload) |
| [`ErrorStateKalmanFilterConfig`](inc/damp/estimation/eskf.hpp#L313) | block | Observers & estimators | Configuration of an ErrorStateKalmanFilter: the process noise and the initial covariance |
| [`ErrorTolerance`](inc/damp/simulation/solver.hpp#L302) | block | Simulation / SIL harness (host) | Scalar-or-vector local-error tolerance for adaptive ODE control |
| [`esc`](inc/damp/controllers/esc.hpp#L156) | function | Design-time synthesis (not PWM-rate) | Synthesize an extremum-seeking controller |
| [`esc_mppt`](inc/damp/controllers/esc.hpp#L200) | function | Design-time synthesis (not PWM-rate) | MPPT-flavored ESC: maximize a power measurement by perturbing the operating point (e.g. converter duty or reference voltage) |
| [`ESCResult`](inc/damp/controllers/esc.hpp#L91) | block | Design-time synthesis (not PWM-rate) | Extremum-seeking design result: the discrete realization the runtime runs |
| [`ESCSpec`](inc/damp/controllers/esc.hpp#L71) | block | Runtime controllers | The extremum-seeking design input: dither, gain, filter corners and clamp; Ts is passed beside it |
| [`eskf_imu`](inc/damp/estimation/eskf.hpp#L220) | function | Observers & estimators | Design Q, R, P₀ for a 6-state IMU attitude ESKF (gyro + accel) |
| [`eskf_marg`](inc/damp/estimation/eskf.hpp#L259) | function | Observers & estimators | Design Q, R, P₀ for a 6-state MARG attitude ESKF (gyro + accel + mag) |
| [`ESKFOrientationFilter`](inc/damp/estimation/sensor_fusion.hpp#L383) | block | Observers & estimators | Turnkey 6-state attitude ESKF: owns q, b_g, and the error filter |
| [`ESKFPredictFn`](inc/damp/estimation/eskf.hpp#L300) | concept | Observers & estimators | Callable that returns ErrorStateJacobian for one predict step |
| [`ESKFResult`](inc/damp/estimation/eskf.hpp#L86) | block | Observers & estimators | Design payload for an ESKF: Q, R, P₀, and success |
| [`eso`](inc/damp/estimation/eso.hpp#L209) | function | Observers & estimators | Extended-state observer: poles at −ω_o, plant gain b₀, optional sample time |
| [`ESO`](inc/damp/estimation/eso.hpp#L312) | block | Observers & estimators | Extended-state observer (ESO) (+1 more overload) |
| [`ESOResult`](inc/damp/estimation/eso.hpp#L57) | block | Observers & estimators | Extended-state observer design: continuous β and, for Ts > 0, the discrete ESO |
| [`estim`](inc/damp/matlab.hpp#L731) | function | MATLAB®-style aliases (host) | Form the state estimator from a plant and an estimator gain |
| [`Euler`](inc/damp/math/geometry.hpp#L226) | block | Scalar math, complex & frames | Intrinsic Euler angles for sequence Order (default ZYX aerospace YPR) |
| [`EulerOrder`](inc/damp/math/geometry.hpp#L55) | enum | Scalar math, complex & frames | Intrinsic Euler angle sequence (axis order of successive principal rotations) |
| [`eval_frf`](inc/damp/systems/state_space.hpp#L234) | function | LTI systems (SS / TF / ZPK / discretize) | Evaluate frequency response of a state-space system |
| [`EventSimConfig`](inc/damp/simulation/hybrid.hpp#L67) | block | Simulation / SIL harness (host) | Configuration for event-driven hybrid / multi-rate style runs |
| [`Exact`](inc/damp/simulation/integrator.hpp#L167) | block | Simulation / SIL harness (host) | Exact integrator for LTI systems |
| [`exact_lti_step`](inc/damp/simulation/hybrid.hpp#L303) | function | Simulation / SIL harness (host) | Exact open-loop step of a single continuous LTI system |
| [`exp`](inc/damp/math/math.hpp#L234) | function | Scalar math, complex & frames | Exponential function |
| [`expm`](inc/damp/matrix/functions.hpp#L696) | function | Linear algebra | Matrix exponential by scaling and squaring (Al-Mohy & Higham 2009) |
| [`expm`](inc/damp/matlab.hpp#L681) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the matrix exponential (mat::expm) |
| [`expo`](inc/damp/toolbox/conditioning.hpp#L204) | function | Embedded helpers (controls-adjacent utilities) | Exponential response curve `y = (1−k)·x + k·x³` (RC "expo") |
| [`ExtendedKalmanFilter`](inc/damp/estimation/ekf.hpp#L118) | block | Observers & estimators | Extended Kalman filter: continuous nonlinear dynamics, sampled measurements (+1 more overload) |
| [`ExtendedKalmanFilterConfig`](inc/damp/estimation/ekf.hpp#L42) | block | Observers & estimators | Configuration of an ExtendedKalmanFilter: its sample time, the process noise and the initial estimate |
| [`extract_channel`](inc/damp/simulation/plot_plotly.hpp#L225) | function | Simulation / SIL harness (host) | Extract the i-th element from each vector entry into a std::vector<double> |
| [`extract_frf_features`](inc/damp/estimation/frequency_response.hpp#L810) | function | Observers & estimators | Extract peaks and valleys in one pass over the FRF table |
| [`extract_modes`](inc/damp/estimation/frequency_response.hpp#L785) | function | Observers & estimators | Extract resonant peaks (local \|G\| maxima) from an FRF table |
| [`extract_valleys`](inc/damp/estimation/frequency_response.hpp#L799) | function | Observers & estimators | Extract anti-resonance valleys (local \|G\| minima) from an FRF table |
| [`Extrapolation`](inc/damp/toolbox/lookup.hpp#L36) | enum | Embedded helpers (controls-adjacent utilities) | Out-of-range behaviour for a Lut1D / Lut2D / Lut3D query beyond its breakpoints |
| [`ExtremumSeekingController`](inc/damp/controllers/esc.hpp#L374) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its design, keeping the result's scalar |
| [`ExtremumSeekingController`](inc/damp/controllers/esc.hpp#L236) | block | Runtime controllers | Extremum-seeking controller runtime (model-free online optimizer) |
| [`ExtremumType`](inc/damp/controllers/esc.hpp#L49) | enum | Runtime controllers | Whether ESC climbs to a maximum or descends to a minimum of the objective |
| [`eye`](inc/damp/matlab.hpp#L609) | function | MATLAB®-style aliases (host) | Create an identity matrix of size n x n |
| [`F_TRIG`](inc/damp/toolbox/iec61131.hpp#L154) | block | Embedded helpers (controls-adjacent utilities) | F_TRIG (Falling Edge Trigger) |
| [`factorial`](inc/damp/trajectory/trajectory_types.hpp#L209) | function | Trajectory value types | k! (exact for modest k in float/double; used up through nonic BVP / jets) |
| [`falling_factorial`](inc/damp/trajectory/trajectory_types.hpp#L222) | function | Trajectory value types | Falling factorial i·(i−1)···(i−k+1) = i! / (i−k)! — the k-th derivative coefficient of tⁱ. Zero when i < k |
| [`FallingEdge`](inc/damp/toolbox/logic.hpp#L48) | block | Embedded helpers (controls-adjacent utilities) | Falling-edge detector: true on the tick x goes true → false |
| [`feedback`](inc/damp/systems/state_space.hpp#L392) | function | LTI systems (SS / TF / ZPK / discretize) | Negative feedback: y = sys1(u − sys2(y)) (+1 more overload) |
| [`feedback_read`](inc/damp/systems/state_space.hpp#L829) | function | LTI systems (SS / TF / ZPK / discretize) | Unity (or weighted) negative feedback around G_loop, read G_out |
| [`FeedbackLaw`](inc/damp/concepts.hpp#L118) | concept | Core, configuration & backend vocabulary | An LQ or MPC law that holds outputs at r: control(r) and step(r, m) |
| [`FeedbackRuntime`](inc/damp/concepts.hpp#L104) | concept | Core, configuration & backend vocabulary | Shared lifecycle of an LQ or MPC runtime: measurement, time update, applied input |
| [`Filtered`](inc/damp/estimation/state_variable_filter.hpp#L67) | block | Observers & estimators | One signal after the state-variable filter: its filtered value and first two derivatives |
| [`finish_extremum`](inc/damp/estimation/frequency_response.hpp#L631) | function | Observers & estimators | Half-power (peaks) or double-power (valleys) bandwidth → ζ |
| [`finite_non_negative`](inc/damp/math/math.hpp#L526) | function | Scalar math, complex & frames | True if x is finite and non-negative (x ≥ 0) |
| [`finite_positive`](inc/damp/math/math.hpp#L516) | function | Scalar math, complex & frames | True if x is finite and strictly greater than zero |
| [`fir1`](inc/damp/filters/fir.hpp#L246) | function | Filters & signal conditioning | Window-method FIR design (normalized frequency, single cutoff) (+1 more overload) |
| [`fir1_hz`](inc/damp/filters/fir.hpp#L387) | function | Filters & signal conditioning | Window-method FIR design from a cutoff in Hz and the sample time (+1 more overload) |
| [`fir_window`](inc/damp/filters/fir.hpp#L423) | function | Filters & signal conditioning | Alias for fir1 (descriptive name) (+1 more overload) |
| [`FirDesignResult`](inc/damp/filters/fir.hpp#L67) | block | Filters & signal conditioning | Window-method FIR design result |
| [`fires_at`](inc/damp/simulation/multirate.hpp#L73) | function | Simulation / SIL harness (host) | True when base tick tick is a firing instant for divisor |
| [`FirFilter`](inc/damp/filters/fir.hpp#L465) | block | Filters & signal conditioning | Direct-form FIR filter runtime (tapped delay line + dot product) (+1 more overload) |
| [`FirstOrderCoeffs`](inc/damp/filters/iir_design.hpp#L31) | block | Filters & signal conditioning | DSP coefficients for first-order IIR filter |
| [`FirstOrderPlantEstimator`](inc/damp/estimation/parameter_estimation.hpp#L126) | block | Observers & estimators | Online estimator for a first-order plant's gain and time constant |
| [`FirstOrderPlantEstimatorConfig`](inc/damp/estimation/parameter_estimation.hpp#L85) | block | Observers & estimators | Configuration for the first-order grey-box estimator |
| [`FirstOrderSection`](inc/damp/filters/first_order_section.hpp#L44) | block | Filters & signal conditioning | First-order IIR section runtime, Direct Form I (+1 more overload) |
| [`FirType`](inc/damp/filters/fir.hpp#L38) | enum | Filters & signal conditioning | Ideal frequency-selective response for fir1 / fir_window |
| [`FirWindow`](inc/damp/filters/fir.hpp#L50) | enum | Filters & signal conditioning | Window applied to the ideal truncated impulse response |
| [`fit_at_crossover`](inc/damp/design/pid_design.hpp#L603) | function | Design-time synthesis (not PWM-rate) | fit_to_response at spec.wc with every gain of sign σ = ±1 |
| [`fit_to_response`](inc/damp/design/pid_design.hpp#L549) | function | Design-time synthesis (not PWM-rate) | The requested law with C(jω) = e^{j(−π + PM)}/P(jω), from the plant response P at ω = spec.wc |
| [`fixed_solve`](inc/damp/simulation/solver.hpp#L286) | function | Simulation / SIL harness (host) | Fixed-step solve in one call |
| [`FixedStepSolver`](inc/damp/simulation/solver.hpp#L145) | block | Simulation / SIL harness (host) | Fixed-step ODE solver (+1 more overload) |
| [`floor`](inc/damp/math/math.hpp#L371) | function | Scalar math, complex & frames | Floor — largest integer ≤ x |
| [`flux_from_Kv`](inc/damp/motor/spm.hpp#L55) | function | Motor control pack (if present) | λ from Kᵥ via torque_constant_from_Kv |
| [`flux_from_torque_constant`](inc/damp/motor/spm.hpp#L37) | function | Motor control pack (if present) | λ from datasheet Kₜ (amplitude / peak-per-phase convention) |
| [`fmod`](inc/damp/math/math.hpp#L415) | function | Scalar math, complex & frames | Floating-point remainder, x − y·trunc(x/y) (sign of x), matching std::fmod's truncated-quotient convention |
| [`FOController`](inc/damp/motor/foc.hpp#L101) | block | Motor control pack (if present) | dq current regulator: PI + decoupling FF → clamped Vdq |
| [`forward_substitute`](inc/damp/matrix/solve.hpp#L37) | function | Linear algebra | Forward substitution for Lx = b |
| [`ForwardEuler`](inc/damp/simulation/integrator.hpp#L223) | block | Simulation / SIL harness (host) | Forward Euler integrator |
| [`fractional_delay`](inc/damp/systems/discretization.hpp#L515) | function | LTI systems (SS / TF / ZPK / discretize) | Fractional-delay all-pass filter of order N (Thiran): z⁻ᴰ for a real D = delay / Ts |
| [`FrequencyPoint`](inc/damp/analysis/frequency.hpp#L37) | block | Frequency-domain analysis (host) | Single-point frequency response result |
| [`FrequencyResponseEstimator`](inc/damp/estimation/frequency_response.hpp#L168) | block | Observers & estimators | Lock-in FRF estimator over a fixed table of NFreq bins |
| [`frf_s_or_z`](inc/damp/analysis/frequency.hpp#L1005) | function | Frequency-domain analysis (host) | FRF evaluation point: jω (continuous) or e^{jωTs} (discrete) |
| [`FrfFeatureKind`](inc/damp/estimation/frequency_response.hpp#L437) | enum | Observers & estimators | Peak (resonance) or valley (anti-resonance) on an FRF magnitude curve |
| [`FrfFeatureResult`](inc/damp/estimation/frequency_response.hpp#L520) | block | Observers & estimators | Combined peak + valley reduction of an FRF table |
| [`FrfMargins`](inc/damp/estimation/frequency_response.hpp#L555) | block | Observers & estimators | Open-loop margins read from a discrete FRF table (on-target) |
| [`FrfPoint`](inc/damp/estimation/frequency_response.hpp#L57) | block | Observers & estimators | One measured frequency-response point (on-target FRF table entry) |
| [`full_qr`](inc/damp/matrix/decomposition.hpp#L533) | function | Linear algebra | Full QR factorization via Householder reflections (real or complex T) |
| [`FullQR`](inc/damp/matrix/decomposition.hpp#L515) | block | Linear algebra | Result of a full (complete) QR factorization |
| [`gain_margin_unwrapped`](inc/damp/analysis/frequency.hpp#L211) | function | Frequency-domain analysis (host) | Find the gain margin using the unwrapped phase trajectory |
| [`gang_of_four`](inc/damp/analysis/sensitivity.hpp#L207) | function | Frequency-domain analysis (host) | Gang of four of a plant P and a controller C on the error, u = C(r − y) |
| [`gang_of_seven`](inc/damp/analysis/sensitivity.hpp#L150) | function | Frequency-domain analysis (host) | Gang of seven of a plant and a two-degree-of-freedom controller u = K[r; y] |
| [`GangOfFourResult`](inc/damp/analysis/sensitivity.hpp#L29) | block | Frequency-domain analysis (host) | Gang of four of a loop u = C(r − y): S, T, PS and CS, each a state-space model on the closed loop's states |
| [`GangOfSevenResult`](inc/damp/analysis/sensitivity.hpp#L45) | block | Frequency-domain analysis (host) | Gang of seven of a two-degree-of-freedom loop u = K[r; y]: the gang of four of its feedback path plus the reference responses |
| [`generalized_eigenvalues`](inc/damp/matrix/generalized_eigen.hpp#L509) | function | Linear algebra | Generalized eigenvalues of the real pencil A − λB: the λ with det(A − λB) = 0 |
| [`GeneralizedEigenResult`](inc/damp/matrix/generalized_eigen.hpp#L61) | block | Linear algebra | Generalized eigenvalues of A − λB: the finite ones, and how many there are |
| [`geomspace`](inc/damp/analysis/linspace.hpp#L68) | function | Frequency-domain analysis (host) | Geometric sequence from start to end (endpoint values, not exponents) (+1 more overload) |
| [`Goertzel`](inc/damp/filters/spectral.hpp#L52) | block | Filters & signal conditioning | Generalized Goertzel single-bin DFT — amplitude/phase at one frequency |
| [`gram`](inc/damp/matlab.hpp#L1498) | function | MATLAB®-style aliases (host) | MATLAB® alias for the controllability/observability Gramian of a system |
| [`gravity_nav`](inc/damp/estimation/ins_mechanization.hpp#L93) | function | Observers & estimators | Gravity vector in the local-level frame [m/s²] |
| [`h_pattern`](inc/damp/toolbox/io.hpp#L285) | function | Embedded helpers (controls-adjacent utilities) | H-pattern (dual-lever / tank) drive: one lever per track, no mixing |
| [`hankel_singular_values`](inc/damp/design/model_reduction.hpp#L508) | function | Design-time synthesis (not PWM-rate) | Hankel singular values of a stable state-space system |
| [`hankelsv`](inc/damp/matlab.hpp#L352) | function | MATLAB®-style aliases (host) | MATLAB® short alias for design::hankel_singular_values |
| [`harmonic_suppressor`](inc/damp/controllers/harmonic_suppression.hpp#L128) | function | Design-time synthesis (not PWM-rate) | Synthesize a multi-resonant harmonic suppressor |
| [`HarmonicAnalyzer`](inc/damp/filters/spectral.hpp#L204) | block | Filters & signal conditioning | Harmonic analyzer — a Goertzel bank over a fundamental and K−1 harmonics |
| [`HarmonicSuppressor`](inc/damp/controllers/harmonic_suppression.hpp#L273) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its design, keeping the result's scalar |
| [`HarmonicSuppressor`](inc/damp/controllers/harmonic_suppression.hpp#L166) | block | Runtime controllers | Multi-resonant harmonic suppressor — a parallel bank of PR resonators |
| [`HarmonicSuppressorResult`](inc/damp/controllers/harmonic_suppression.hpp#L50) | block | Design-time synthesis (not PWM-rate) | Design result for a multi-resonant harmonic suppressor |
| [`hcat`](inc/damp/matrix/concat.hpp#L159) | function | Linear algebra | Horizontal concatenation (same row count) |
| [`heading_from_baseline_nav`](inc/damp/estimation/ins_eskf.hpp#L422) | function | Observers & estimators | Heading [rad] from a nav-frame baseline vector (dual-antenna difference) |
| [`Heun`](inc/damp/simulation/integrator.hpp#L863) | block | Simulation / SIL harness (host) | Heun's method (Improved Euler, RK2) integrator |
| [`highpass_1st`](inc/damp/filters/iir_design.hpp#L417) | function | Filters & signal conditioning | First-order high-pass filter design (continuous-time) (+1 more overload) |
| [`highpass_2nd`](inc/damp/filters/iir_design.hpp#L553) | function | Filters & signal conditioning | Second-order high-pass filter design (continuous-time) (+1 more overload) |
| [`highpass_sections`](inc/damp/filters/biquad.hpp#L500) | function | Filters & signal conditioning | Nth-order Butterworth high-pass at cutoff fc as second-order sections, any order |
| [`highshelf`](inc/damp/filters/iir_design.hpp#L905) | function | Filters & signal conditioning | High-shelf EQ filter: boost or cut everything above fc |
| [`hinfnorm`](inc/damp/matlab.hpp#L1521) | function | MATLAB®-style aliases (host) | MATLAB® alias for the H∞ system norm [ninf, fpeak] = hinfnorm(sys, tol) |
| [`HinfNormResult`](inc/damp/analysis/norms.hpp#L119) | block | Frequency-domain analysis (host) | H∞ norm of a stable system and the frequency of its peak |
| [`holonomic_drive`](inc/damp/toolbox/io.hpp#L386) | function | Embedded helpers (controls-adjacent utilities) | Four-wheel holonomic drive mixer (mecanum / omni) |
| [`HolonomicOutput`](inc/damp/toolbox/io.hpp#L355) | block | Embedded helpers (controls-adjacent utilities) | Four-wheel actuator set (mecanum/holonomic drive output), X-config |
| [`hurwitz`](inc/damp/design/pid_design.hpp#L628) | function | Design-time synthesis (not PWM-rate) | Routh–Hurwitz: every root of c[0] sⁿ + … + c[n] in the open left half plane |
| [`HybridEventAction`](inc/damp/simulation/hybrid.hpp#L104) | block | Simulation / SIL harness (host) | Action applied when an event fires (after Exact advance to the event time) |
| [`HybridSimulationResult`](inc/damp/simulation/hybrid.hpp#L79) | block | Simulation / SIL harness (host) | Result of a piecewise-LTI hybrid simulation |
| [`hypot`](inc/damp/math/math.hpp#L78) | function | Scalar math, complex & frames | Euclidean distance hypot(x, y) = √(x² + y²), without overflow |
| [`Hysteresis`](inc/damp/toolbox/conditioning.hpp#L284) | block | Embedded helpers (controls-adjacent utilities) | Hysteresis comparator (Schmitt trigger): bool output with separate on/off thresholds to reject chatter |
| [`impedance`](inc/damp/analysis/frequency.hpp#L1073) | function | Frequency-domain analysis (host) | Compute impedance frequency response from a SISO impedance transfer function |
| [`impedance_from_admittance`](inc/damp/analysis/frequency.hpp#L1031) | function | Frequency-domain analysis (host) | Compute impedance frequency response from a SISO admittance system |
| [`ImpedanceResult`](inc/damp/analysis/frequency.hpp#L903) | block | Frequency-domain analysis (host) | Result of impedance frequency response evaluation |
| [`impulse`](inc/damp/analysis/time_response.hpp#L198) | function | Frequency-domain analysis (host) | Impulse response of a (MIMO) state-space system (+1 more overload) |
| [`impulse`](inc/damp/estimation/excitation/impulse.hpp#L73) | function | Observers & estimators | Build an impulse design payload |
| [`Impulse`](inc/damp/estimation/excitation/impulse.hpp#L89) | block | Observers & estimators | Rectangular impulse (finite-width Dirac stand-in) |
| [`ImpulseConfig`](inc/damp/estimation/excitation/impulse.hpp#L32) | block | Observers & estimators | Configuration for a one-shot rectangular impulse (Dirac stand-in) |
| [`impulseplot`](inc/damp/simulation/plot_plotly.hpp#L735) | function | Simulation / SIL harness (host) | Plot an impulse response, one trace per input/output pair |
| [`ImpulseResult`](inc/damp/estimation/excitation/impulse.hpp#L52) | block | Observers & estimators | Impulse design payload |
| [`ImuSample`](inc/damp/estimation/ins_mechanization.hpp#L125) | block | Observers & estimators | IMU sample in the body frame |
| [`infinity_norm`](inc/damp/matrix/functions.hpp#L42) | function | Linear algebra | Infinity norm ‖A‖∞: maximum absolute row sum |
| [`initial`](inc/damp/analysis/time_response.hpp#L246) | function | Frequency-domain analysis (host) | Initial-condition (free) response of a (MIMO) state-space system |
| [`InputDelay`](inc/damp/systems/discretization.hpp#L559) | block | LTI systems (SS / TF / ZPK / discretize) | An input delay of θ = seconds, in up to N registers (0 ≤ θ ≤ N·Ts): each register state holds an input already sent, v(k−1) … v(k−N). `discretize(sys, Ts, InputDelay<N>{theta})` returns that model; design on it like any other. `InputDelay{theta}` is N = 1 (+1 more overload) |
| [`InputDisturbance`](inc/damp/design/disturbance.hpp#L42) | block | Design-time synthesis (not PWM-rate) | An input disturbance: d adds to u, B_d = B and C_d = D, one disturbance per input |
| [`ins_aid_position`](inc/damp/estimation/ins_eskf.hpp#L389) | function | Observers & estimators | Predicted antenna / marker position with body lever-arm |
| [`ins_apply_correction`](inc/damp/estimation/ins_eskf.hpp#L375) | function | Observers & estimators | Inject filter correction into x and reset the error state |
| [`ins_error_jacobian`](inc/damp/estimation/ins_eskf.hpp#L283) | function | Observers & estimators | Discrete error-state transition for one IMU step, on its sparsity pattern |
| [`ins_eskf_design`](inc/damp/estimation/ins_eskf.hpp#L194) | function | Observers & estimators | Build Q / R / P0 for a 15-state INS ESKF from IMU noise densities |
| [`ins_heading`](inc/damp/estimation/ins_eskf.hpp#L407) | function | Observers & estimators | Horizontal heading [rad] of a body axis expressed in the nav frame |
| [`ins_inject`](inc/damp/estimation/ins_eskf.hpp#L326) | function | Observers & estimators | Inject δx into the nominal INS state (right-multiplicative q) |
| [`ins_predict`](inc/damp/estimation/ins_eskf.hpp#L355) | function | Observers & estimators | Mechanize nominal state and propagate the error covariance |
| [`ins_update_heading`](inc/damp/estimation/ins_eskf.hpp#L510) | function | Observers & estimators | Sparse dual-antenna / yaw heading update (scalar) |
| [`ins_update_orientation`](inc/damp/estimation/ins_eskf.hpp#L548) | function | Observers & estimators | Sparse orientation update from a measured body→nav quaternion |
| [`ins_update_pose`](inc/damp/estimation/ins_eskf.hpp#L572) | function | Observers & estimators | Generic pose alias: position then orientation (two sparse updates) |
| [`ins_update_pose_heading`](inc/damp/estimation/ins_eskf.hpp#L591) | function | Observers & estimators | Pose alias: position + dual-antenna heading (two sparse updates) |
| [`ins_update_position`](inc/damp/estimation/ins_eskf.hpp#L437) | function | Observers & estimators | Sparse position update: y = p + Rℓ |
| [`ins_update_velocity`](inc/damp/estimation/ins_eskf.hpp#L469) | function | Observers & estimators | Sparse velocity update: y = v at the IMU origin (no lever-arm) |
| [`ins_update_zupt`](inc/damp/estimation/ins_eskf.hpp#L488) | function | Observers & estimators | Zero-velocity update (ZUPT): velocity measurement of zero |
| [`InsNavigator`](inc/damp/estimation/ins_eskf.hpp#L625) | block | Observers & estimators | Runtime strapdown navigator: nominal InsState + 15-state ESKF (+2 more overloads) |
| [`InsState`](inc/damp/estimation/ins_mechanization.hpp#L139) | block | Observers & estimators | Nominal strapdown navigation state (external to the ESKF error state) |
| [`instantaneous_power`](inc/damp/math/transforms.hpp#L496) | function | Scalar math, complex & frames | Instantaneous active and reactive power from dq quantities (+1 more overload) |
| [`InstantaneousPower`](inc/damp/math/transforms.hpp#L457) | block | Scalar math, complex & frames | Instantaneous active and reactive power |
| [`IntegralAction`](inc/damp/controllers/state_feedback_core.hpp#L715) | block | Runtime controllers | Offset-free by integral action: integrators on r − y, every measured output (`IntegralAction{}`) or the listed ones (`IntegralAction<2>{{0, 3}}`, in that order: the references r follow it) |
| [`IntegralDesign`](inc/damp/controllers/state_feedback_core.hpp#L767) | block | Runtime controllers | The integral action of a design: the integrated outputs, their maps of the state and the input, and the per-sample step |
| [`IntegrationResult`](inc/damp/simulation/integrator.hpp#L64) | block | Simulation / SIL harness (host) | Result of an integration step |
| [`InteriorPointSolver`](inc/damp/design/qp.hpp#L1446) | block | Design-time synthesis (not PWM-rate) | Interior-point QP solver policy for damp::MPC |
| [`Interpolation`](inc/damp/toolbox/lookup.hpp#L42) | enum | Embedded helpers (controls-adjacent utilities) | In-range interpolant for Lut1D (2-D stays bilinear; use SplineSurface for C¹ grids) |
| [`InvariantSetResult`](inc/damp/design/invariant_set.hpp#L41) | block | Design-time synthesis (not PWM-rate) | Maximal positively invariant polytope G e ≤ h of x⁺ = A x |
| [`inverse`](inc/damp/matrix/solve.hpp#L294) | function | Linear algebra | Matrix inverse A⁻¹ via mat::solve(A, I) |
| [`inverse_clarke_transform`](inc/damp/math/transforms.hpp#L308) | function | Scalar math, complex & frames | Inverse Clarke transform (αβ → abc) |
| [`inverse_deadband`](inc/damp/toolbox/conditioning.hpp#L81) | function | Embedded helpers (controls-adjacent utilities) | Inverse dead zone: add an offset to overcome a physical dead zone (valve overlap, static friction, motor stiction), with independent negative/positive offsets (+1 more overload) |
| [`inverse_deadband_sine`](inc/damp/toolbox/conditioning.hpp#L164) | function | Embedded helpers (controls-adjacent utilities) | Sine-shaped inverse dead zone (full boost once \|x\|≥w) |
| [`inverse_deadband_soft`](inc/damp/toolbox/conditioning.hpp#L140) | function | Embedded helpers (controls-adjacent utilities) | Smooth inverse dead zone via soft-sign (tanh-like Coulomb boost) |
| [`inverse_lerp`](inc/damp/toolbox/scaling.hpp#L48) | function | Embedded helpers (controls-adjacent utilities) | Inverse of lerp: the fraction t such that `lerp(a, b, t) == x` |
| [`inverse_park_clarke_transform`](inc/damp/math/transforms.hpp#L423) | function | Scalar math, complex & frames | Fused inverse Park-Clarke transform (dq → abc) |
| [`inverse_park_transform`](inc/damp/math/transforms.hpp#L382) | function | Scalar math, complex & frames | Inverse Park transform (dq → αβ) |
| [`inverse_symmetrical_components`](inc/damp/math/transforms.hpp#L571) | function | Scalar math, complex & frames | Inverse symmetrical-component transform (012 → abc) |
| [`iq_from_torque`](inc/damp/motor/spm.hpp#L71) | function | Motor control pack (if present) | i_q = T_e / Kₜ for id = 0 |
| [`is_controllable`](inc/damp/design/stability.hpp#L447) | function | Design-time synthesis (not PWM-rate) | Check if a system is controllable |
| [`is_fault`](inc/damp/toolbox/conditioning.hpp#L334) | function | Embedded helpers (controls-adjacent utilities) | True for a wire fault (FaultLow/FaultHigh) — i.e. not a real reading at all |
| [`is_matrix_element`](inc/damp/matrix/matrix_traits.hpp#L42) | block | Linear algebra | True if T may be a Matrix element type |
| [`is_observable`](inc/damp/design/stability.hpp#L461) | function | Design-time synthesis (not PWM-rate) | Check if a system is observable |
| [`is_stabilizable`](inc/damp/design/stability.hpp#L484) | function | Design-time synthesis (not PWM-rate) | Check if (A, B) of sys is a stabilizable pair, in the domain of sys |
| [`is_stable`](inc/damp/design/stability.hpp#L502) | function | Design-time synthesis (not PWM-rate) | Stability of a state-space model, in the domain of sys |
| [`is_symmetric_or_hermitian`](inc/damp/matrix/decomposition.hpp#L91) | function | Linear algebra | True if A equals its transpose (real) or conjugate transpose (complex) |
| [`is_usable`](inc/damp/toolbox/conditioning.hpp#L508) | function | Embedded helpers (controls-adjacent utilities) | True when a channel reading may still be used (in-span or saturated) |
| [`is_valid`](inc/damp/toolbox/conditioning.hpp#L329) | function | Embedded helpers (controls-adjacent utilities) | True only for the in-span status |
| [`isfinite`](inc/damp/math/math.hpp#L503) | function | Scalar math, complex & frames | Finiteness test — false for NaN and ±∞ |
| [`iso_c_drive`](inc/damp/toolbox/io.hpp#L304) | function | Embedded helpers (controls-adjacent utilities) | ISO-C (single-stick, coordinated / curvature) drive |
| [`isstable`](inc/damp/matlab.hpp#L1355) | function | MATLAB®-style aliases (host) | Stability predicate for a state-space model (+1 more overload) |
| [`Jet`](inc/damp/trajectory/trajectory_types.hpp#L98) | block | Trajectory value types | Truncated jet of a scalar motion sample: d[k] = s⁽ᵏ⁾ |
| [`JordanBlock`](inc/damp/design/pole_placement.hpp#L488) | block | Design-time synthesis (not PWM-rate) | One Jordan mini-block of a desired closed-loop spectrum |
| [`JordanObjective`](inc/damp/design/pole_placement.hpp#L1200) | enum | Design-time synthesis (not PWM-rate) | How place_jordan spends the free parameter of the Klein–Moore family (Schmid et al.) |
| [`JordanOptions`](inc/damp/design/pole_placement.hpp#L1213) | block | Design-time synthesis (not PWM-rate) | Options for place_jordan |
| [`JordanPlaceResult`](inc/damp/design/pole_placement.hpp#L1228) | block | Design-time synthesis (not PWM-rate) | place_jordan result: the PlaceResult fields plus the achieved robustness metrics |
| [`kalman`](inc/damp/estimation/kalman.hpp#L485) | function | Observers & estimators | Steady-state Kalman filter design (+3 more overloads) |
| [`KalmanFilter`](inc/damp/estimation/kalman.hpp#L667) | block | Observers & estimators | The recursive Kalman filter: the covariance propagates each tick, and the gain is formed from it in Joseph form (+2 more overloads) |
| [`KalmanFilterConfig`](inc/damp/estimation/kalman.hpp#L618) | block | Observers & estimators | Configuration of a KalmanFilter built from a model rather than from design::kalman |
| [`KalmanNoise`](inc/damp/estimation/kalman.hpp#L80) | block | Observers & estimators | The estimator's noise: process noise on the plant's noise inputs (through G), measurement noise on the sensors (through H) (+2 more overloads) |
| [`KalmanResult`](inc/damp/estimation/kalman.hpp#L138) | block | Observers & estimators | Steady-state Kalman filter design result |
| [`kalmd`](inc/damp/matlab.hpp#L982) | function | MATLAB®-style aliases (host) | Discrete Kalman estimator for a continuous plant and continuous noise (MATLAB®'s kalmd) |
| [`lag`](inc/damp/controllers/lead_lag.hpp#L197) | function | Design-time synthesis (not PWM-rate) | Design a lag compensator from desired low-frequency gain boost |
| [`lambda_tuning`](inc/damp/design/pid_design.hpp#L429) | function | Design-time synthesis (not PWM-rate) | Lambda tuning for FOPDT model |
| [`Latch`](inc/damp/toolbox/logic.hpp#L66) | block | Embedded helpers (controls-adjacent utilities) | Set/reset latch. @tparam SetDominant which input wins when both are asserted (default: set-dominant, e.g. a trip overriding a clear) |
| [`Law`](inc/damp/design/model_inverse.hpp#L953) | block | Design-time synthesis (not PWM-rate) | Law<kRes, T> of a `constexpr` model inverse: the same runtime on the design's matrices as constants of the type (+1 more overload) |
| [`Law`](inc/damp/estimation/eso.hpp#L461) | block | Observers & estimators | Law<kRes, T> of a design::eso initialized at start-up: the ESO it runs as. A runtime of the same shape assigned replaces it (+5 more overloads) |
| [`Law`](inc/damp/controllers/law_detail.hpp#L234) | block | Runtime controllers | The runtime a design result deploys as, named by the result: `Law<kRes, float> law;` (+6 more overloads) |
| [`LDLT`](inc/damp/matrix/decomposition.hpp#L163) | block | Linear algebra | LDLᴴ factors: A = L D Lᴴ with L unit lower-triangular and D real diagonal |
| [`ldlt`](inc/damp/matrix/decomposition.hpp#L191) | function | Linear algebra | LDLᴴ decomposition of a symmetric (Hermitian) matrix, no pivoting |
| [`ldlt_solve`](inc/damp/matrix/solve.hpp#L189) | function | Linear algebra | Solve AX = B via LDLᴴ (A = L D Lᴴ, symmetric indefinite) |
| [`lead`](inc/damp/controllers/lead_lag.hpp#L144) | function | Design-time synthesis (not PWM-rate) | Design a lead compensator from desired phase boost at a target frequency |
| [`lead_lag`](inc/damp/controllers/lead_lag.hpp#L258) | function | Design-time synthesis (not PWM-rate) | Design a lead-lag compensator (cascade of lead + lag sections) |
| [`lead_lag_direct`](inc/damp/controllers/lead_lag.hpp#L293) | function | Design-time synthesis (not PWM-rate) | Direct lead-lag specification from zero/pole locations |
| [`LeadLagController`](inc/damp/controllers/lead_lag.hpp#L436) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its design, keeping the result's scalar |
| [`LeadLagController`](inc/damp/controllers/lead_lag.hpp#L316) | block | Runtime controllers | Discrete lead-lag compensator |
| [`LeadLagResult`](inc/damp/controllers/lead_lag.hpp#L64) | block | Design-time synthesis (not PWM-rate) | Lead-lag compensator design result |
| [`LeadLagSeriesResult`](inc/damp/controllers/lead_lag.hpp#L224) | block | Design-time synthesis (not PWM-rate) | Cascaded lead+lag design result (2nd-order StateSpace + success) |
| [`lerp`](inc/damp/toolbox/scaling.hpp#L35) | function | Embedded helpers (controls-adjacent utilities) | Linear interpolation between a and b by fraction t |
| [`leverrier_num_scale`](inc/damp/systems/zpk.hpp#L585) | function | LTI systems (SS / TF / ZPK / discretize) | Magnitude scale of each strictly proper Leverrier numerator coefficient |
| [`LeverrierResult`](inc/damp/systems/zpk.hpp#L519) | block | LTI systems (SS / TF / ZPK / discretize) | Faddeev–LeVerrier characteristic polynomial and adjoint coefficient matrices |
| [`lft`](inc/damp/systems/state_space.hpp#L506) | function | LTI systems (SS / TF / ZPK / discretize) | Lower linear fractional transformation: K closes the loop from P's last outputs to its last inputs |
| [`LIMIT`](inc/damp/toolbox/iec61131.hpp#L593) | function | Embedded helpers (controls-adjacent utilities) | LIMIT (IEC 61131-3 selection function): clamp in to [mn, mx] |
| [`linear_modulation_voltage`](inc/damp/motor/modulation.hpp#L127) | function | Motor control pack (if present) | Peak phase voltage at the linear SVPWM hexagon limit |
| [`linear_screw`](inc/damp/toolbox/actuator.hpp#L170) | function | Embedded helpers (controls-adjacent utilities) | Build a ServoAxis for a linear axis driven by a leadscrew/belt |
| [`LinearizationResult`](inc/damp/design/linearization.hpp#L46) | block | Design-time synthesis (not PWM-rate) | Result of nonlinear operating-point linearization |
| [`linearize`](inc/damp/design/linearization.hpp#L118) | function | Design-time synthesis (not PWM-rate) | Linearize a Dynamics / Measurement model about an operating point (+1 more overload) |
| [`linmod`](inc/damp/matlab.hpp#L588) | function | MATLAB®-style aliases (host) | MATLAB®-style nonlinear linearization about an operating point |
| [`linspace`](inc/damp/analysis/linspace.hpp#L36) | function | Frequency-domain analysis (host) | num evenly spaced values from start to end inclusive (+1 more overload) |
| [`locate_zero_crossing_exact`](inc/damp/simulation/hybrid.hpp#L333) | function | Simulation / SIL harness (host) | Locate a state zero-crossing of g(x) on an Exact LTI segment |
| [`log`](inc/damp/matrix/functions.hpp#L1031) | function | Linear algebra | Principal matrix logarithm via inverse scaling and squaring |
| [`log`](inc/damp/math/math.hpp#L252) | function | Scalar math, complex & frames | Natural logarithm |
| [`log10`](inc/damp/math/math.hpp#L434) | function | Scalar math, complex & frames | Base-10 logarithm, log10(x) = ln(x) / ln(10) |
| [`logm`](inc/damp/matlab.hpp#L707) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the principal matrix logarithm (mat::log) |
| [`LogPolicy`](inc/damp/simulation/hybrid.hpp#L58) | enum | Simulation / SIL harness (host) | How often to record samples along a hybrid trajectory |
| [`logspace`](inc/damp/matlab.hpp#L1217) | function | MATLAB®-style aliases (host) | n points from 10ᵃ to 10ᵇ, logarithmically spaced (decade exponents, as MATLAB®'s) |
| [`loop_at_a_time_disk_margin`](inc/damp/analysis/disk_margin.hpp#L194) | function | Frequency-domain analysis (host) | Loop-at-a-time disk margins of a square MIMO loop L under unity negative feedback |
| [`loop_metrics`](inc/damp/analysis/frequency.hpp#L482) | function | Frequency-domain analysis (host) | One-call loop analysis: compute L/S/T response and return compact metrics (+1 more overload) |
| [`loop_response`](inc/damp/analysis/frequency.hpp#L584) | function | Frequency-domain analysis (host) | Compute open-loop L, sensitivity S, complementary sensitivity T, and Nyquist data (+1 more overload) |
| [`LoopDiskMarginResult`](inc/damp/analysis/disk_margin.hpp#L55) | block | Frequency-domain analysis (host) | Loop-at-a-time disk margins of a square MIMO loop: one SISO disk margin per channel |
| [`LoopResponseResult`](inc/damp/analysis/frequency.hpp#L417) | block | Frequency-domain analysis (host) | Open-loop and closed-loop frequency response package |
| [`LoopSummary`](inc/damp/analysis/frequency.hpp#L450) | block | Frequency-domain analysis (host) | Compact loop summary metrics for quick stability/robustness checks |
| [`LowerTriangle`](inc/damp/matrix/views.hpp#L165) | block | Linear algebra | Lower triangular view of a square matrix |
| [`lowpass`](inc/damp/filters/lowpass.hpp#L100) | function | Filters & signal conditioning | Discretize a continuous low-pass TF (Tustin) into direct-form coefficients (+1 more overload) |
| [`LowPass`](inc/damp/filters/lowpass.hpp#L223) | block | Filters & signal conditioning | Nth-order low-pass filter (direct form I) (+1 more overload) |
| [`lowpass_1st`](inc/damp/filters/iir_design.hpp#L360) | function | Filters & signal conditioning | First-order low-pass filter design (continuous-time) (+1 more overload) |
| [`lowpass_2nd`](inc/damp/filters/iir_design.hpp#L474) | function | Filters & signal conditioning | Second-order low-pass filter design (continuous-time) (+1 more overload) |
| [`lowpass_sections`](inc/damp/filters/biquad.hpp#L464) | function | Filters & signal conditioning | Nth-order Butterworth low-pass at cutoff fc as second-order sections, any order |
| [`LowPassResult`](inc/damp/filters/lowpass.hpp#L58) | block | Filters & signal conditioning | Direct-form coefficients of an Nth-order low-pass, with success / status |
| [`lowshelf`](inc/damp/filters/iir_design.hpp#L873) | function | Filters & signal conditioning | Low-shelf EQ filter: boost or cut everything below fc |
| [`lqg`](inc/damp/controllers/lqg.hpp#L42) | function | Design-time synthesis (not PWM-rate) | The linear-quadratic Gaussian law: LQR on a steady-state Kalman filter's estimate (design::kalman(sys, noise)): runs as LQG (+3 more overloads) |
| [`LQG`](inc/damp/controllers/lqg.hpp#L231) | function | Design-time synthesis (not PWM-rate) | LQG law{res}: the sizes, scalar and filter of design::lqg's result |
| [`lqg`](inc/damp/matlab.hpp#L1097) | function | MATLAB®-style aliases (host) | Linear-Quadratic-Gaussian regulator design |
| [`LQG`](inc/damp/controllers/lqg.hpp#L142) | block | Runtime controllers | Linear-quadratic Gaussian law: LQR on the Kalman filter's estimate |
| [`lqg_from_parts`](inc/damp/controllers/lqg.hpp#L110) | function | Design-time synthesis (not PWM-rate) | Assemble an LQG design from separately computed Kalman and LQR designs |
| [`lqgi`](inc/damp/controllers/lqgi.hpp#L39) | function | Design-time synthesis (not PWM-rate) | The linear-quadratic Gaussian law with integral action on every output: runs as LQGI (+3 more overloads) |
| [`LQGI`](inc/damp/controllers/lqgi.hpp#L148) | function | Design-time synthesis (not PWM-rate) | LQGI law{res}: the sizes, scalar and filter of design::lqgi's result |
| [`LQGI`](inc/damp/controllers/lqgi.hpp#L98) | block | Runtime controllers | Linear-quadratic Gaussian law with integral action on the measured outputs |
| [`lqgreg`](inc/damp/matlab.hpp#L1113) | function | MATLAB®-style aliases (host) | Combine separate Kalman filter and LQR designs into an LQG controller |
| [`lqgtrack`](inc/damp/matlab.hpp#L1134) | function | MATLAB®-style aliases (host) | Linear-Quadratic-Gaussian design with integral action for tracking (+2 more overloads) |
| [`lqi`](inc/damp/controllers/lqi.hpp#L39) | function | Design-time synthesis (not PWM-rate) | The linear-quadratic regulator with integral action on every output (xi += Ts·(r − Ci·y)): runs as LQI (+3 more overloads) |
| [`LQI`](inc/damp/controllers/lqi.hpp#L137) | function | Design-time synthesis (not PWM-rate) | LQI law{res}: the sizes and scalar of design::lqi's result |
| [`lqi`](inc/damp/matlab.hpp#L1058) | function | MATLAB®-style aliases (host) | Linear-Quadratic Integral design for tracking (+2 more overloads) |
| [`LQI`](inc/damp/controllers/lqi.hpp#L89) | block | Runtime controllers | Linear-quadratic regulator with integral action, on the measured state |
| [`lqr`](inc/damp/controllers/lqr.hpp#L45) | function | Design-time synthesis (not PWM-rate) | The linear-quadratic regulator on the plant (DARE on a discrete plant, CARE on a continuous one): runs as LQR (+3 more overloads) |
| [`LQR`](inc/damp/controllers/lqr.hpp#L205) | function | Design-time synthesis (not PWM-rate) | LQR law{res}: the sizes and scalar of design::lqr's result (+2 more overloads) |
| [`lqr`](inc/damp/matlab.hpp#L892) | function | MATLAB®-style aliases (host) | Continuous-time LQR design (MATLAB®'s lqr) |
| [`LQR`](inc/damp/controllers/lqr.hpp#L110) | block | Runtime controllers | Linear-quadratic regulator on the measured state: u = u − K·(x − x ) |
| [`lqr_gain`](inc/damp/design/riccati.hpp#L1359) | function | Design-time synthesis (not PWM-rate) | Optimal LQR state-feedback gain from a Riccati solution |
| [`LQRCost`](inc/damp/design/riccati.hpp#L1577) | block | Design-time synthesis (not PWM-rate) | A quadratic cost on (x, u): state cost Q, input cost R, cross term N |
| [`lqrd`](inc/damp/matlab.hpp#L922) | function | MATLAB®-style aliases (host) | Sampled-data LQR from a continuous plant: ZOH dynamics + Van Loan cost (MATLAB®'s lqrd) (+1 more overload) |
| [`LQRResult`](inc/damp/design/riccati.hpp#L1518) | block | Design-time synthesis (not PWM-rate) | The LQR gain: design::lqr on a plant, discrete_lqr and continuous_lqr on given matrices |
| [`LQRWeights`](inc/damp/controllers/state_feedback_core.hpp#L755) | function | Design-time synthesis (not PWM-rate) | Sizes from the matrices: `LQRWeights{Q, R}` (+2 more overloads) |
| [`LQRWeights`](inc/damp/controllers/state_feedback_core.hpp#L746) | block | Runtime controllers | The cost of a state-feedback design: what to penalize and how much |
| [`lqry`](inc/damp/matlab.hpp#L1008) | function | MATLAB®-style aliases (host) | Output-weighted continuous LQR (state cost Q = Cᵀ Q_y C) (+1 more overload) |
| [`lsim`](inc/damp/analysis/time_response.hpp#L284) | function | Frequency-domain analysis (host) | Forced time response of a (MIMO) state-space system to an input signal (+1 more overload) |
| [`LsimInfo`](inc/damp/analysis/time_response.hpp#L513) | block | Frequency-domain analysis (host) | Transient characteristics of an arbitrary response signal |
| [`lsiminfo`](inc/damp/analysis/time_response.hpp#L533) | function | Frequency-domain analysis (host) | Compute transient characteristics from an output/time signal |
| [`lsimplot`](inc/damp/simulation/plot_plotly.hpp#L751) | function | Simulation / SIL harness (host) | Plot a forced (lsim) simulation, one trace per output |
| [`LsimResult`](inc/damp/analysis/time_response.hpp#L66) | block | Frequency-domain analysis (host) | Result of a single-trajectory simulation: time, output, and state history |
| [`lu_decomposition`](inc/damp/matrix/decomposition.hpp#L254) | function | Linear algebra | LU decomposition with partial pivoting |
| [`lu_solve`](inc/damp/matrix/solve.hpp#L220) | function | Linear algebra | Solve AX = B via LU with partial pivoting |
| [`luenberger`](inc/damp/estimation/luenberger.hpp#L300) | function | Observers & estimators | Design a Luenberger observer by robust pole placement, in the domain of sys (+1 more overload) |
| [`Luenberger`](inc/damp/estimation/luenberger.hpp#L604) | block | Observers & estimators | The fixed-gain state observer: the innovation gain M held, the current form, and the optional disturbance anchor (+4 more overloads) |
| [`LuenbergerResult`](inc/damp/estimation/luenberger.hpp#L119) | block | Observers & estimators | Luenberger observer design result: the plant, the gains L and M, and the achieved error poles |
| [`Lut1D`](inc/damp/toolbox/lookup.hpp#L109) | block | Embedded helpers (controls-adjacent utilities) | 1-D interpolating lookup table over monotonic breakpoints |
| [`Lut2D`](inc/damp/toolbox/lookup.hpp#L201) | block | Embedded helpers (controls-adjacent utilities) | 2-D bilinear interpolating lookup table over a regular grid |
| [`Lut3D`](inc/damp/toolbox/lookup.hpp#L275) | block | Embedded helpers (controls-adjacent utilities) | 3-D trilinear interpolating lookup table over a regular grid |
| [`lut_segment`](inc/damp/toolbox/lookup.hpp#L63) | function | Embedded helpers (controls-adjacent utilities) | Index of the interpolation segment containing x |
| [`lyap`](inc/damp/design/lyapunov.hpp#L444) | function | Design-time synthesis (not PWM-rate) | Solve the continuous-time Lyapunov equation A X + X Aᵀ + Q = 0 |
| [`lyap`](inc/damp/matlab.hpp#L1412) | function | MATLAB®-style aliases (host) | MATLAB® alias for the continuous Lyapunov solve AX+XAᵀ+Q=0 |
| [`lyapchol`](inc/damp/matlab.hpp#L1432) | function | MATLAB®-style aliases (host) | MATLAB® alias for the Cholesky factor of the continuous Lyapunov solution (X = RᵀR) |
| [`lyapunov_cholesky`](inc/damp/design/lyapunov.hpp#L492) | function | Design-time synthesis (not PWM-rate) | Cholesky factor of the continuous Lyapunov solution: A X + X Aᵀ + B Bᵀ = 0, X = RᵀR |
| [`MadgwickFilter`](inc/damp/estimation/sensor_fusion.hpp#L135) | block | Observers & estimators | Madgwick gradient-descent AHRS filter |
| [`mag2db`](inc/damp/math/math.hpp#L538) | function | Scalar math, complex & frames | Magnitude to decibels, 20·log10(mag) |
| [`MahonyFilter`](inc/damp/estimation/sensor_fusion.hpp#L177) | block | Observers & estimators | Mahony nonlinear complementary filter with PI correction |
| [`make1DOF`](inc/damp/matlab.hpp#L304) | function | MATLAB®-style aliases (host) | Force 1-DOF setpoint weights on a PID design result (b=c=1) |
| [`make2DOF`](inc/damp/matlab.hpp#L320) | function | MATLAB®-style aliases (host) | Apply 2-DOF setpoint weights on a PID design result |
| [`margin`](inc/damp/matlab.hpp#L1266) | function | MATLAB®-style aliases (host) | Gain and phase margins of a SISO loop over a frequency grid |
| [`MarginResult`](inc/damp/matlab.hpp#L1245) | block | MATLAB®-style aliases (host) | Gain/phase margins and their crossover frequencies |
| [`margins_from_frf`](inc/damp/estimation/frequency_response.hpp#L842) | function | Observers & estimators | Gain / phase margins from a fixed on-target FRF table |
| [`markov_checked_numerator`](inc/damp/systems/zpk.hpp#L615) | function | LTI systems (SS / TF / ZPK / discretize) | num, or the numerator the Markov parameters give over den when num misses it |
| [`Matrix`](inc/damp/matrix/core.hpp#L253) | block | Linear algebra | Fixed-size, stack-allocated matrix for linear algebra operations |
| [`MatrixLike`](inc/damp/matrix/matrix_traits.hpp#L98) | concept | Linear algebra | Concept for any type that provides 2D matrix-like element access |
| [`MatrixLikeOf`](inc/damp/matrix/matrix_traits.hpp#L108) | concept | Linear algebra | Concept for a MatrixLike type with specific dimensions |
| [`max_iq`](inc/damp/motor/spm.hpp#L111) | function | Motor control pack (if present) | Max positive i_q for an SPM (i_d = 0) at electrical speed |
| [`max_torque_at_speed`](inc/damp/motor/spm.hpp#L172) | function | Motor control pack (if present) | SPM max-torque point at a mechanical speed (i_d = 0 only) |
| [`maximal_invariant_set`](inc/damp/design/invariant_set.hpp#L251) | function | Design-time synthesis (not PWM-rate) | Maximal positively invariant set of x⁺ = A x under G₀ x ≤ h₀ |
| [`MeasuredDisturbance`](inc/damp/controllers/mpc.hpp#L1255) | function | Design-time synthesis (not PWM-rate) | Sizes from the matrix: `MeasuredDisturbance{Bd}` |
| [`MeasuredDisturbance`](inc/damp/controllers/mpc.hpp#L1249) | block | Runtime controllers | A measured disturbance: a known input the controller does not command, entering the state update through Bd (a load torque, a feed rate, the ambient temperature) |
| [`Measurement`](inc/damp/concepts.hpp#L208) | concept | Core, configuration & backend vocabulary | Nonlinear measurement y = h(x, u): a callable returning y as a ColVec<NY, T> |
| [`MeasurementJacobian`](inc/damp/concepts.hpp#L221) | block | Core, configuration & backend vocabulary | Analytic Jacobians of a Measurement model at (x, u), returned by its jacobian(x, u) |
| [`mechanize_step`](inc/damp/estimation/ins_mechanization.hpp#L210) | function | Observers & estimators | One dead-reckoning step from a raw IMU sample (subtracts the state's biases) |
| [`mechanize_step_from_corrected`](inc/damp/estimation/ins_mechanization.hpp#L175) | function | Observers & estimators | One dead-reckoning step from an IMU sample (strapdown mechanization) |
| [`MedianFilter`](inc/damp/filters/median.hpp#L34) | block | Filters & signal conditioning | Sliding-window median filter — nonlinear spike/outlier rejection |
| [`mhe`](inc/damp/estimation/mhe.hpp#L227) | function | Observers & estimators | Synthesize a constrained moving-horizon estimator |
| [`MHE`](inc/damp/estimation/mhe.hpp#L346) | block | Observers & estimators | Runtime moving-horizon estimator (fixed per-tick iteration budget) (+1 more overload) |
| [`MHEConstraints`](inc/damp/estimation/mhe.hpp#L97) | block | Observers & estimators | Box constraints on the estimated states (applied to every window state) and the QP budget |
| [`MHEResult`](inc/damp/estimation/mhe.hpp#L117) | block | Observers & estimators | Window QP data produced by mhe(), consumed by damp::MHE |
| [`middlebrook`](inc/damp/analysis/frequency.hpp#L1113) | function | Frequency-domain analysis (host) | Middlebrook stability analysis of a cascaded source and load from their impedances |
| [`MiddlebrookResult`](inc/damp/analysis/frequency.hpp#L927) | block | Frequency-domain analysis (host) | Result of Middlebrook minor loop gain analysis |
| [`MinimalRealizationResult`](inc/damp/design/minreal.hpp#L89) | block | Design-time synthesis (not PWM-rate) | Minimal realization result |
| [`minmax`](inc/damp/backend.hpp#L164) | function | Core, configuration & backend vocabulary | Ordered {min, max} pair returned by value (+1 more overload) |
| [`minreal`](inc/damp/design/minreal.hpp#L311) | function | Design-time synthesis (not PWM-rate) | Minimal realization — cancel uncontrollable and unobservable modes |
| [`minreal`](inc/damp/matlab.hpp#L334) | function | MATLAB®-style aliases (host) | MATLAB® short alias for design::minreal, with its default tolerance √ε |
| [`minreal_zpk`](inc/damp/systems/zpk.hpp#L836) | function | LTI systems (SS / TF / ZPK / discretize) | Cancel matching pole-zero pairs on a ZPK model |
| [`mix`](inc/damp/systems/state_space.hpp#L711) | function | LTI systems (SS / TF / ZPK / discretize) | Fold plant outputs into one input: v = g u + Σ kᵢ yᵢ (+1 more overload) |
| [`mix_outputs`](inc/damp/systems/state_space.hpp#L788) | function | LTI systems (SS / TF / ZPK / discretize) | Linear combination of outputs: y = Σ αⱼ yⱼ |
| [`modal_form`](inc/damp/systems/realization.hpp#L286) | function | LTI systems (SS / TF / ZPK / discretize) | Real modal canonical form: A block diagonal, one block per mode |
| [`ModeExtractorConfig`](inc/damp/estimation/frequency_response.hpp#L475) | block | Observers & estimators | Configuration for FRF peak / valley extraction |
| [`ModeExtractorResult`](inc/damp/estimation/frequency_response.hpp#L492) | block | Observers & estimators | Result of reducing an FRF table to at most MaxModes peaks (resonances) |
| [`ModelReductionMethod`](inc/damp/design/model_reduction.hpp#L61) | enum | Design-time synthesis (not PWM-rate) | How balanced_reduction and eliminate_states eliminate states |
| [`modred`](inc/damp/matlab.hpp#L375) | function | MATLAB®-style aliases (host) | MATLAB® short alias for design::eliminate_states (0-based indices) |
| [`modulation_duty_cycles`](inc/damp/motor/modulation.hpp#L400) | function | Motor control pack (if present) | Carrier-based VSI duty cycles from an αβ voltage command |
| [`modulation_index`](inc/damp/motor/modulation.hpp#L165) | function | Motor control pack (if present) | Modulation index relative to the linear SVPWM circle |
| [`modulation_zero_sequence`](inc/damp/motor/modulation.hpp#L327) | function | Motor control pack (if present) | Zero-sequence for any PwmScheme |
| [`Modulator`](inc/damp/motor/modulation.hpp#L497) | block | Motor control pack (if present) | Thin scheme-holding wrapper for FOC / deploy paths |
| [`motor_constant`](inc/damp/motor/spm.hpp#L63) | function | Motor control pack (if present) | Motor constant Kₘ = Kₜ / √R [Nm/√W] |
| [`MovingAverage`](inc/damp/filters/moving_average.hpp#L67) | block | Filters & signal conditioning | Moving-average (boxcar) filter — also a DC-preserving harmonic-notch comb |
| [`mpc`](inc/damp/controllers/mpc.hpp#L1660) | function | Design-time synthesis (not PWM-rate) | Synthesize a constrained linear MPC (condensed dense QP, Δu form) on a discrete plant (+3 more overloads) |
| [`MPC`](inc/damp/controllers/mpc.hpp#L2523) | function | Design-time synthesis (not PWM-rate) | MPC law{res}: the sizes, scalar and terminal capacity of design::mpc's result (+1 more overload) |
| [`MPC`](inc/damp/controllers/mpc.hpp#L2453) | block | Runtime controllers | Constrained linear MPC on the measured state (design::mpc) |
| [`MPCArtifacts`](inc/damp/controllers/mpc.hpp#L351) | block | Runtime controllers | Condensed-QP data produced by mpc() (MPCResult::mpc), consumed by damp::MPC |
| [`MPCConstraints`](inc/damp/controllers/mpc.hpp#L294) | block | Runtime controllers | Box constraints for mpc() |
| [`MPCHorizonSuggestion`](inc/damp/controllers/mpc.hpp#L1749) | block | Runtime controllers | Advisory NP/NC horizon values from suggest_mpc_horizon() |
| [`MPCLimits`](inc/damp/controllers/mpc.hpp#L328) | block | Runtime controllers | The limits an MPC takes at run time (MATLAB®'s online constraints): u, Δu and the controlled outputs, in physical units. A field left out is unbounded |
| [`MPCResult`](inc/damp/controllers/mpc.hpp#L1268) | block | Design-time synthesis (not PWM-rate) | An MPC design: the condensed QP, the estimator when there is one, and the steady-state target the runtime poses the QP about; consumed by damp::MPC |
| [`MPCTerminal`](inc/damp/controllers/mpc.hpp#L230) | enum | Runtime controllers | Opt-in terminal cost / set for `design::mpc` (default: none) |
| [`MPCWeights`](inc/damp/controllers/mpc.hpp#L256) | block | Runtime controllers | Cost weights for mpc() |
| [`mstogi`](inc/damp/filters/sogi.hpp#L153) | function | Filters & signal conditioning | Mixed Second/Third-Order Generalized Integrator (MSTOGI) (+1 more overload) |
| [`MSTOGI`](inc/damp/filters/sogi.hpp#L299) | block | Filters & signal conditioning | Runtime MSTOGI with exact resonator and forward-Euler washout |
| [`mul_linear`](inc/damp/systems/zpk.hpp#L165) | function | LTI systems (SS / TF / ZPK / discretize) | Multiply a polynomial of degree by (s − root) in place; poly[k] multiplies sᵏ |
| [`multi_sine`](inc/damp/estimation/excitation/multi_sine.hpp#L139) | function | Observers & estimators | Build a multi-sine design payload from a configuration |
| [`MultiRateConfig`](inc/damp/simulation/multirate.hpp#L81) | block | Simulation / SIL harness (host) | Configuration for multi-rate closed-loop runs |
| [`MultiSine`](inc/damp/estimation/excitation/multi_sine.hpp#L213) | block | Observers & estimators | Sum-of-tones multi-sine runtime generator (+1 more overload) |
| [`MultiSineConfig`](inc/damp/estimation/excitation/multi_sine.hpp#L62) | block | Observers & estimators | Configuration for fixed-component multi-sine excitation |
| [`MultiSineResult`](inc/damp/estimation/excitation/multi_sine.hpp#L104) | block | Observers & estimators | Multi-sine design payload |
| [`MUX`](inc/damp/toolbox/iec61131.hpp#L612) | function | Embedded helpers (controls-adjacent utilities) | MUX (IEC 61131-3 multiplexer): select input k of N (0-based) |
| [`NavFrame`](inc/damp/estimation/ins_mechanization.hpp#L77) | enum | Observers & estimators | Local-level navigation frame for strapdown mechanization |
| [`nearbyint`](inc/damp/math/math.hpp#L399) | function | Scalar math, complex & frames | Round to nearest integer; ties to even (IEEE default / `FE_TONEAREST`) |
| [`ned_from_enu`](inc/damp/estimation/ins_mechanization.hpp#L106) | function | Observers & estimators | Map an ENU vector into NED (axis permute) |
| [`negative_sequence_ab`](inc/damp/filters/pll.hpp#L267) | function | Filters & signal conditioning | Instantaneous negative-sequence αβ from a quadrature signal pair |
| [`nichols`](inc/damp/analysis/frequency.hpp#L714) | function | Frequency-domain analysis (host) | Build Nichols points from existing Bode data (+2 more overloads) |
| [`nicholsplot`](inc/damp/simulation/plot_plotly.hpp#L944) | function | Simulation / SIL harness (host) | Plot a Nichols chart (open-loop phase vs magnitude) with M-circle grid |
| [`NicholsPoint`](inc/damp/analysis/frequency.hpp#L685) | block | Frequency-domain analysis (host) | Single-point Nichols chart sample (open-loop phase vs magnitude) |
| [`NicholsResult`](inc/damp/analysis/frequency.hpp#L697) | block | Frequency-domain analysis (host) | Nichols chart locus across a frequency sweep |
| [`NlmsFilter`](inc/damp/filters/fir.hpp#L575) | block | Filters & signal conditioning | Normalized LMS adaptive FIR filter |
| [`NoDisturbance`](inc/damp/design/disturbance.hpp#L33) | block | Design-time synthesis (not PWM-rate) | No disturbance model (the default): the plain filter on x |
| [`NoEstimator`](inc/damp/controllers/state_feedback_core.hpp#L54) | block | Runtime controllers | The estimator of a law that measures its state: none, an empty type, so a `[[no_unique_address]]` member of it takes no storage |
| [`norm`](inc/damp/matlab.hpp#L1511) | function | MATLAB®-style aliases (host) | MATLAB® alias for the H2 system norm norm(sys,2) |
| [`norm_h2`](inc/damp/analysis/norms.hpp#L78) | function | Frequency-domain analysis (host) | H2 norm of a state-space system |
| [`norm_hinf`](inc/damp/analysis/norms.hpp#L337) | function | Frequency-domain analysis (host) | H∞ norm of a stable state-space system, sup σ̄(G) over the imaginary axis or the unit circle, and its frequency |
| [`notch`](inc/damp/filters/iir_design.hpp#L661) | function | Filters & signal conditioning | Second-order band-reject (notch) filter design (continuous-time) (+1 more overload) |
| [`null`](inc/damp/matlab.hpp#L654) | function | MATLAB®-style aliases (host) | MATLAB® short alias for an orthonormal null-space basis |
| [`null_space`](inc/damp/matrix/svd.hpp#L458) | function | Linear algebra | Orthonormal basis for the null space {x : A·x = 0} via SVD |
| [`NullSpace`](inc/damp/matrix/svd.hpp#L437) | block | Linear algebra | Orthonormal basis for the null space (kernel) of a matrix |
| [`nyquist`](inc/damp/analysis/frequency.hpp#L499) | function | Frequency-domain analysis (host) | Compute Nyquist data for a SISO state-space system (+1 more overload) |
| [`nyquistplot`](inc/damp/simulation/plot_plotly.hpp#L823) | function | Simulation / SIL harness (host) | Plot a Nyquist locus with the -1 critical point marked |
| [`NyquistPoint`](inc/damp/analysis/frequency.hpp#L371) | block | Frequency-domain analysis (host) | Single-point Nyquist response data |
| [`NyquistResult`](inc/damp/analysis/frequency.hpp#L383) | block | Frequency-domain analysis (host) | Nyquist response data across a frequency sweep |
| [`observability_gramian`](inc/damp/design/stability.hpp#L410) | function | Design-time synthesis (not PWM-rate) | Continuous/discrete observability Gramian W_o |
| [`observability_matrix`](inc/damp/design/stability.hpp#L85) | function | Design-time synthesis (not PWM-rate) | Compute the observability matrix [C; CA; CA²; ...; CA^(N-1)] |
| [`observability_staircase`](inc/damp/systems/realization.hpp#L851) | function | LTI systems (SS / TF / ZPK / discretize) | Observability staircase form: the observable part split off by an orthogonal transform |
| [`observable_canonical_form`](inc/damp/systems/realization.hpp#L569) | function | LTI systems (SS / TF / ZPK / discretize) | Observable canonical form of a SISO model: ones on A's subdiagonal, −α in its last column, C = eₙᵀ |
| [`obsv`](inc/damp/matlab.hpp#L402) | function | MATLAB®-style aliases (host) | MATLAB® short alias for observability_matrix (+1 more overload) |
| [`obsvf`](inc/damp/matlab.hpp#L551) | function | MATLAB®-style aliases (host) | MATLAB® obsvf(A, B, C): observability staircase form, its orthogonal transform and rank |
| [`OffDelayTimer`](inc/damp/toolbox/logic.hpp#L136) | block | Embedded helpers (controls-adjacent utilities) | Off-delay timer: output goes true immediately when in is true and stays true until in has been false continuously for delay |
| [`oflqg`](inc/damp/controllers/oflqg.hpp#L64) | function | Design-time synthesis (not PWM-rate) | The offset-free linear-quadratic Gaussian law: a Kalman filter on [x; d] and the steady-state target that holds the controlled outputs at r under d̂: runs as OFLQG (+1 more overload) |
| [`OFLQG`](inc/damp/controllers/oflqg.hpp#L153) | function | Design-time synthesis (not PWM-rate) | OFLQG law{res}: the sizes, scalar, filter and target of design::oflqg's result |
| [`OFLQG`](inc/damp/controllers/oflqg.hpp#L102) | block | Runtime controllers | Offset-free linear-quadratic Gaussian law: a Kalman filter on [x; d] and the steady-state target |
| [`ofmpc`](inc/damp/controllers/ofmpc.hpp#L43) | function | Design-time synthesis (not PWM-rate) | The offset-free MPC on a discrete plant: a Kalman filter on [x; d] (noise, disturbance) and the QP posed about the steady-state target under d̂: runs as OFMPC (+1 more overload) |
| [`OFMPC`](inc/damp/controllers/ofmpc.hpp#L138) | function | Design-time synthesis (not PWM-rate) | OFMPC law{res}: the sizes, scalar, filter and target of design::ofmpc's result (+1 more overload) |
| [`OFMPC`](inc/damp/controllers/ofmpc.hpp#L83) | block | Runtime controllers | Offset-free constrained linear MPC: a Kalman filter on [x; d] and the QP about the steady-state target under d̂ (design::ofmpc) |
| [`OnDelayTimer`](inc/damp/toolbox/logic.hpp#L100) | block | Embedded helpers (controls-adjacent utilities) | On-delay timer: output goes true once in has been held true continuously for delay; drops immediately when in goes false |
| [`one_norm`](inc/damp/matrix/functions.hpp#L71) | function | Linear algebra | One-norm ‖A‖₁: maximum absolute column sum |
| [`ones`](inc/damp/matlab.hpp#L618) | function | MATLAB®-style aliases (host) | N×N matrix of ones |
| [`option`](inc/damp/design/options.hpp#L58) | function | Design-time synthesis (not PWM-rate) | No option left: the fallback (+1 more overload) |
| [`option_set`](inc/damp/design/options.hpp#L83) | block | Design-time synthesis (not PWM-rate) | The kinds of option a design call accepts: check<Opts...>() refuses, at compile time, an option of no listed kind and a kind given twice |
| [`option_type`](inc/damp/design/options.hpp#L42) | block | Design-time synthesis (not PWM-rate) | The type of the first option of the kind Is among Opts, else Default (+1 more overload) |
| [`output_inverse`](inc/damp/design/model_inverse.hpp#L686) | function | Design-time synthesis (not PWM-rate) | The stable inverse of a discrete model: the state and input of its trajectory through a reference y on the controlled outputs, from y ₖ … y ₖ₊NP₋₁ (+1 more overload) |
| [`OutputDisturbance`](inc/damp/design/disturbance.hpp#L53) | block | Design-time synthesis (not PWM-rate) | An output disturbance: a bias on each reading, B_d = 0 and C_d = I, one per output |
| [`OutputInverse`](inc/damp/design/model_inverse.hpp#L864) | block | Design-time synthesis (not PWM-rate) | The stable inverse at run time: the reference point (x , u ) from the output reference over the preview, seven matrix–vector products and the zero dynamics' state (+1 more overload) |
| [`OutputInverseResult`](inc/damp/design/model_inverse.hpp#L143) | block | Design-time synthesis (not PWM-rate) | The stable inverse of a discrete model over a preview of NP output samples |
| [`pade`](inc/damp/matlab.hpp#L1384) | function | MATLAB®-style aliases (host) | Padé approximation of order N (1 or 2) of the pure delay e^{−sT} |
| [`pade_delay_1st`](inc/damp/filters/pade.hpp#L33) | function | Filters & signal conditioning | First-order Pade approximation of time delay (+1 more overload) |
| [`pade_delay_2nd`](inc/damp/filters/pade.hpp#L85) | function | Filters & signal conditioning | Second-order Pade approximation of time delay (+1 more overload) |
| [`panel_color`](inc/damp/simulation/plot_plotly.hpp#L316) | function | Simulation / SIL harness (host) | Plotly default colorway entry for index i within a single subplot |
| [`panel_legend`](inc/damp/simulation/plot_plotly.hpp#L299) | function | Simulation / SIL harness (host) | Per-panel legend box just right of a stacked subplot band |
| [`panel_legend_id`](inc/damp/simulation/plot_plotly.hpp#L285) | function | Simulation / SIL harness (host) | Legend id for subplot row (1-based): "legend", "legend2", … |
| [`panel_x_domain`](inc/damp/simulation/plot_plotly.hpp#L275) | function | Simulation / SIL harness (host) | Horizontal domain for stacked multi-panel x-axes |
| [`parallel`](inc/damp/systems/state_space.hpp#L325) | function | LTI systems (SS / TF / ZPK / discretize) | Parallel connection (shared input, summed outputs) |
| [`ParameterDriftMonitor`](inc/damp/estimation/parameter_estimation.hpp#L252) | block | Observers & estimators | Residual-based drift detector: has the plant moved away from the model? |
| [`ParameterEstimator`](inc/damp/concepts.hpp#L159) | concept | Core, configuration & backend vocabulary | Online grey-box parameter estimator: physical parameters + gating |
| [`park_transform`](inc/damp/math/transforms.hpp#L355) | function | Scalar math, complex & frames | Park transform (αβ → dq) |
| [`PatternEntry`](inc/damp/matrix/pattern_matrix.hpp#L25) | enum | Linear algebra | What an entry of a sparsity pattern is known to be at compile time |
| [`PatternMatrix`](inc/damp/matrix/pattern_matrix.hpp#L219) | block | Linear algebra | The runtime values of a sparsity pattern, and nothing else: an entry the pattern knows has no storage, and in a product it has no term (a 0) or no multiply (a 1) |
| [`peaking`](inc/damp/filters/iir_design.hpp#L848) | function | Filters & signal conditioning | Peaking (bell) EQ filter: boost or cut a band around f0 |
| [`Periodic`](inc/damp/toolbox/timing.hpp#L170) | block | Embedded helpers (controls-adjacent utilities) | Periodic trigger — fires once per elapsed period |
| [`phase_margin_from_damping_ratio`](inc/damp/design/pid_design.hpp#L888) | function | Design-time synthesis (not PWM-rate) | Approximate phase margin from damping ratio |
| [`phase_margin_unwrapped`](inc/damp/analysis/frequency.hpp#L169) | function | Frequency-domain analysis (host) | Find the phase margin using the unwrapped phase trajectory |
| [`PhaseUnwrapper`](inc/damp/encoder/unwrap.hpp#L40) | block | Embedded helpers (controls-adjacent utilities) | Counts whole periods from successive fine phases (no hardware counter) |
| [`pid`](inc/damp/controllers/pid.hpp#L363) | function | Design-time synthesis (not PWM-rate) | Continuous PID gains Kₚ, Kᵢ, K_d with every other knob at its default |
| [`pid`](inc/damp/matlab.hpp#L249) | function | MATLAB®-style aliases (host) | MATLAB®-style parallel-form continuous PID constructor |
| [`pid_pole_placement`](inc/damp/design/pid_design.hpp#L960) | function | Design-time synthesis (not PWM-rate) | PI gains that place the closed-loop poles of a first-order plant b₀/(a₁ s + a₀) (+4 more overloads) |
| [`PIDController`](inc/damp/controllers/pid.hpp#L904) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its gains, keeping their scalar: PIDController pid{gains} is a PID; through an alias the mode comes from the alias: PIController pi{gains}; |
| [`PIDController`](inc/damp/controllers/pid.hpp#L515) | block | Runtime controllers | Fixed-rate discrete 2-DOF PID (canonical runtime) (+3 more overloads) |
| [`PIDGains`](inc/damp/controllers/pid.hpp#L137) | block | Runtime controllers | 2-DOF continuous-time PID gains, limits and setpoint weights |
| [`PIDMode`](inc/damp/controllers/pid.hpp#L376) | enum | Runtime controllers | Compile-time selection of the PID control-law structure |
| [`PIDPolePlacementSpec`](inc/damp/design/pid_design.hpp#L912) | block | Design-time synthesis (not PWM-rate) | Closed-loop targets for the continuous pid_pole_placement forms |
| [`PIDRuntimeMode`](inc/damp/controllers/pid.hpp#L401) | enum | Runtime controllers | Runtime operating mode for PIDController / ContinuousPID |
| [`pidstd`](inc/damp/matlab.hpp#L272) | function | MATLAB®-style aliases (host) | Standard-form continuous PID constructor (1-DOF) |
| [`pidstd2`](inc/damp/matlab.hpp#L291) | function | MATLAB®-style aliases (host) | Standard-form continuous 2-DOF PID constructor |
| [`pidtune`](inc/damp/design/pid_design.hpp#L807) | function | Design-time synthesis (not PWM-rate) | Plant-aware PID tune at a crossover |
| [`pidtune`](inc/damp/matlab.hpp#L1182) | function | MATLAB®-style aliases (host) | Plant-aware PID tune at crossover ωc (MATLAB® pidtune spelling): 60° PID |
| [`PIDTuneSpec`](inc/damp/design/pid_design.hpp#L449) | block | Design-time synthesis (not PWM-rate) | Knobs for plant-aware pidtune |
| [`PiecewiseLTI`](inc/damp/simulation/hybrid.hpp#L93) | block | Simulation / SIL harness (host) | Plant: NModes continuous LTI pieces (same state/input/output dimensions) |
| [`pinv`](inc/damp/matlab.hpp#L641) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the Moore–Penrose pseudoinverse |
| [`place`](inc/damp/design/pole_placement.hpp#L1164) | function | Design-time synthesis (not PWM-rate) | Robust pole placement (KNV): the gain K with eig(A − B·K) at poles, in the domain of sys (+1 more overload) |
| [`place`](inc/damp/matlab.hpp#L873) | function | MATLAB®-style aliases (host) | Robust multi-input pole placement (MATLAB®'s place) |
| [`place_jordan`](inc/damp/design/pole_placement.hpp#L1685) | function | Design-time synthesis (not PWM-rate) | Exact pole placement with an arbitrary Jordan structure (Schmid–Ntogramatzidis–Nguyen–Pandey / Klein–Moore parametric form) (+1 more overload) |
| [`PlaceResult`](inc/damp/design/pole_placement.hpp#L1085) | block | Design-time synthesis (not PWM-rate) | Pole-placement design result: the gain K, the achieved closed-loop poles e, and Ts |
| [`PlantDisturbance`](inc/damp/design/disturbance.hpp#L67) | block | Design-time synthesis (not PWM-rate) | The plant's own disturbance channel: B_d = G and C_d = H, one disturbance per column of G (H needs as many columns) |
| [`plot_bode`](inc/damp/simulation/plot_plotly.hpp#L424) | function | Simulation / SIL harness (host) | Plot Bode magnitude and phase as subplots |
| [`plot_line`](inc/damp/simulation/plot_plotly.hpp#L493) | function | Simulation / SIL harness (host) | Simple line plot of time vs value |
| [`plot_simulation`](inc/damp/simulation/plot_plotly.hpp#L345) | function | Simulation / SIL harness (host) | Plot simulation results with subplots for states, outputs, and inputs |
| [`plot_step`](inc/damp/simulation/plot_plotly.hpp#L623) | function | Simulation / SIL harness (host) | Plot step response data |
| [`plot_xy`](inc/damp/simulation/plot_plotly.hpp#L558) | function | Simulation / SIL harness (host) | Multi-series planar scatter (true 2-D path / locus) (+1 more overload) |
| [`plots_root_for`](inc/damp/simulation/plot_plotly.hpp#L120) | function | Simulation / SIL harness (host) | Walk up from dir until a directory named "plots" is found |
| [`pole`](inc/damp/matlab.hpp#L1229) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the open-loop poles of a system |
| [`PoleInfo`](inc/damp/analysis/poles.hpp#L56) | block | Frequency-domain analysis (host) | Natural frequency and damping ratio for each pole |
| [`poles`](inc/damp/analysis/poles.hpp#L40) | function | Frequency-domain analysis (host) | Compute open-loop poles (eigenvalues of A matrix) |
| [`PoleZeroMap`](inc/damp/analysis/poles.hpp#L102) | block | Frequency-domain analysis (host) | Poles and zeros of a system, for pole-zero plotting |
| [`poly_horner`](inc/damp/toolbox/scaling.hpp#L138) | function | Embedded helpers (controls-adjacent utilities) | Evaluate a polynomial at x by Horner's method |
| [`poly_roots`](inc/damp/analysis/poles.hpp#L150) | function | Frequency-domain analysis (host) | Roots of a polynomial given highest power first (MATLAB® `roots`) |
| [`poly_roots`](inc/damp/systems/zpk.hpp#L124) | function | LTI systems (SS / TF / ZPK / discretize) | Compute roots of a real polynomial given highest power first |
| [`PolyRootsResult`](inc/damp/systems/zpk.hpp#L93) | block | LTI systems (SS / TF / ZPK / discretize) | Roots of a real polynomial given highest power first |
| [`Pose`](inc/damp/kinematics/pose.hpp#L61) | block | Kinematics / pose | Rigid-body pose: a translation and an orientation (unit quaternion) |
| [`positive_sequence_ab`](inc/damp/filters/pll.hpp#L246) | function | Filters & signal conditioning | Instantaneous positive-sequence αβ from a quadrature signal pair |
| [`pow`](inc/damp/matrix/functions.hpp#L1133) | function | Linear algebra | Integer matrix power via binary exponentiation (+1 more overload) |
| [`pow`](inc/damp/math/math.hpp#L282) | function | Scalar math, complex & frames | Power function with a real exponent, base^exponent (+1 more overload) |
| [`pr`](inc/damp/controllers/pr.hpp#L205) | function | Design-time synthesis (not PWM-rate) | Design a Proportional-Resonant controller |
| [`pr_harmonics`](inc/damp/controllers/pr.hpp#L243) | function | Design-time synthesis (not PWM-rate) | Design multiple-harmonic PR controller gains |
| [`prbs`](inc/damp/estimation/excitation/prbs.hpp#L166) | function | Observers & estimators | Build a PRBS design payload from a configuration |
| [`PRBS`](inc/damp/estimation/excitation/prbs.hpp#L189) | block | Observers & estimators | Maximal-length PRBS runtime generator (+1 more overload) |
| [`PRBSConfig`](inc/damp/estimation/excitation/prbs.hpp#L82) | block | Observers & estimators | Configuration for maximal-length pseudo-random binary excitation |
| [`PRBSResult`](inc/damp/estimation/excitation/prbs.hpp#L126) | block | Observers & estimators | PRBS design payload |
| [`PRController`](inc/damp/controllers/pr.hpp#L438) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its design, keeping the result's scalar |
| [`PRController`](inc/damp/controllers/pr.hpp#L291) | block | Runtime controllers | Discrete Proportional-Resonant Controller |
| [`prescale`](inc/damp/matlab.hpp#L481) | function | MATLAB®-style aliases (host) | MATLAB® prescale: diagonal state scaling that balances A, B and C |
| [`project_affine`](inc/damp/controllers/action_governor.hpp#L137) | function | Design-time synthesis (not PWM-rate) | Project u_des onto the polyhedron A u ≤ b (Euclidean) |
| [`project_box`](inc/damp/controllers/action_governor.hpp#L76) | function | Design-time synthesis (not PWM-rate) | Euclidean projection of u onto a Bounds box |
| [`PRResult`](inc/damp/controllers/pr.hpp#L111) | block | Design-time synthesis (not PWM-rate) | Proportional-Resonant controller design result |
| [`PRSpec`](inc/damp/controllers/pr.hpp#L98) | block | Runtime controllers | The gains and frequencies of a proportional-resonant design; Ts is passed beside it |
| [`pseudo_inverse`](inc/damp/matrix/svd.hpp#L404) | function | Linear algebra | Moore–Penrose pseudoinverse A⁺ via SVD |
| [`PulseTimer`](inc/damp/toolbox/logic.hpp#L171) | block | Embedded helpers (controls-adjacent utilities) | Pulse timer (non-retriggerable): a rising edge of in emits a fixed |
| [`PwmScheme`](inc/damp/motor/modulation.hpp#L79) | enum | Motor control pack (if present) | Carrier-based three-phase VSI modulation scheme |
| [`pzmap`](inc/damp/analysis/poles.hpp#L173) | function | Frequency-domain analysis (host) | Pole-zero map of a SISO transfer function (MATLAB® `pzmap(tf)`) (+2 more overloads) |
| [`pzplot`](inc/damp/simulation/plot_plotly.hpp#L850) | function | Simulation / SIL harness (host) | Plot a pole-zero map on the complex plane (poles as ×, zeros as ○) |
| [`QPResult`](inc/damp/design/qp.hpp#L84) | block | Design-time synthesis (not PWM-rate) | Result of a dense QP solve |
| [`QPStatus`](inc/damp/design/qp.hpp#L51) | enum | Design-time synthesis (not PWM-rate) | Termination status of a QP solve |
| [`qr_decompose`](inc/damp/matrix/decomposition.hpp#L444) | function | Linear algebra | Thin QR via Householder reflections |
| [`QRDecomposition`](inc/damp/matrix/decomposition.hpp#L397) | block | Linear algebra | Thin QR factorization A = QR (Householder) |
| [`QuadMode`](inc/damp/encoder/quadrature.hpp#L47) | enum | Embedded helpers (controls-adjacent utilities) | Quadrature decode resolution (edges counted per A/B cycle) |
| [`quadprog`](inc/damp/matlab.hpp#L1453) | function | MATLAB®-style aliases (host) | MATLAB® alias for the dense inequality-constrained QP solve (+1 more overload) |
| [`quadratic_form`](inc/damp/matrix/core.hpp#L1282) | function | Linear algebra | Symmetric congruence (quadratic) form  S = M X Mᵀ (+1 more overload) |
| [`quadratic_form_in_place`](inc/damp/matrix/pattern_matrix.hpp#L345) | function | Linear algebra | S ← F S Fᵀ in place, for a symmetric S: every multiply the pattern can skip, skipped |
| [`QuadratureDecoder`](inc/damp/encoder/quadrature.hpp#L75) | block | Embedded helpers (controls-adjacent utilities) | Software A/B quadrature decoder with optional index |
| [`quantize`](inc/damp/simulation/quantize.hpp#L66) | function | Simulation / SIL harness (host) | ADC reading of x: x rounded to the nearest multiple of lsb, then clamped to [lo, hi] |
| [`Quaternion`](inc/damp/math/geometry.hpp#L362) | block | Scalar math, complex & frames | Unit quaternion rotation (w, x, y, z) (Hamilton product) |
| [`R_TRIG`](inc/damp/toolbox/iec61131.hpp#L119) | block | Embedded helpers (controls-adjacent utilities) | R_TRIG (Rising Edge Trigger) |
| [`rad2deg`](inc/damp/math/math.hpp#L563) | function | Scalar math, complex & frames | Radians to degrees, rad·180/π |
| [`raise_loop_for_tracking`](inc/damp/design/pid_design.hpp#L713) | function | Design-time synthesis (not PWM-rate) | Raise the loop crossover so T −3 dB matches the 1-DOF target |
| [`ramp`](inc/damp/estimation/excitation/ramp.hpp#L101) | function | Observers & estimators | Build a ramp design payload from a configuration |
| [`Ramp`](inc/damp/estimation/excitation/ramp.hpp#L117) | block | Observers & estimators | Rate-limited ramp runtime generator (+1 more overload) |
| [`RampConfig`](inc/damp/estimation/excitation/ramp.hpp#L33) | block | Observers & estimators | Configuration for a slew-rate-limited ramp excitation |
| [`RampResult`](inc/damp/estimation/excitation/ramp.hpp#L65) | block | Observers & estimators | Ramp design payload |
| [`RangeMonitor`](inc/damp/toolbox/conditioning.hpp#L403) | block | Embedded helpers (controls-adjacent utilities) | Analog-input range/fault monitor (NAMUR NE43 pattern) |
| [`rank`](inc/damp/design/stability.hpp#L425) | function | Design-time synthesis (not PWM-rate) | Numerical rank of a matrix (forwards to mat::rank) (+1 more overload) |
| [`rank`](inc/damp/matrix/functions.hpp#L234) | function | Linear algebra | Numerical rank: the number of singular values above rel_tol·σₘₐₓ |
| [`rank_from_svd`](inc/damp/matrix/svd.hpp#L368) | function | Linear algebra | Numerical rank from a precomputed SVD result |
| [`reduced_luenberger`](inc/damp/estimation/luenberger.hpp#L507) | function | Observers & estimators | Design a reduced-order (Gopinath) observer by pole placement, in the domain of sys (+1 more overload) |
| [`ReducedLuenberger`](inc/damp/estimation/luenberger.hpp#L1047) | block | Observers & estimators | Reduced-order (Gopinath) state observer (runtime) (+1 more overload) |
| [`ReducedLuenbergerResult`](inc/damp/estimation/luenberger.hpp#L348) | block | Observers & estimators | Reduced-order (Gopinath) observer design result |
| [`RedundantAnalogPair`](inc/damp/toolbox/io.hpp#L91) | block | Embedded helpers (controls-adjacent utilities) | Dual-channel analog input with same-slope or opposite-slope cross-check |
| [`reg`](inc/damp/matlab.hpp#L762) | function | MATLAB®-style aliases (host) | Form the observer-based compensator y → u from a plant, K and L |
| [`relay_autotune`](inc/damp/estimation/relay_autotune.hpp#L228) | function | Observers & estimators | Build a validated relay-autotune design payload |
| [`RelayAutotuneConfig`](inc/damp/estimation/relay_autotune.hpp#L128) | block | Observers & estimators | Configuration for the relay-feedback autotuning experiment |
| [`RelayAutotuneOutput`](inc/damp/estimation/relay_autotune.hpp#L260) | block | Observers & estimators | Per-tick output of RelayAutotuner::step |
| [`RelayAutotuner`](inc/damp/estimation/relay_autotune.hpp#L279) | block | Observers & estimators | Runtime relay-feedback autotuner (+1 more overload) |
| [`RelayAutotuneResult`](inc/damp/estimation/relay_autotune.hpp#L187) | block | Observers & estimators | Relay-autotuner design payload |
| [`RelayAutotuneStatus`](inc/damp/estimation/relay_autotune.hpp#L244) | enum | Observers & estimators | Lifecycle state of a RelayAutotuner |
| [`repetitive`](inc/damp/controllers/repetitive.hpp#L214) | function | Design-time synthesis (not PWM-rate) | Synthesize a repetitive controller with a scalar robustness filter Q |
| [`repetitive_binomial`](inc/damp/controllers/repetitive.hpp#L247) | function | Design-time synthesis (not PWM-rate) | Synthesize a repetitive controller with a binomial zero-phase FIR Q |
| [`RepetitiveConfig`](inc/damp/controllers/repetitive.hpp#L77) | block | Runtime controllers | Repetitive-controller tuning + period (with optional zero-phase FIR Q) |
| [`RepetitiveController`](inc/damp/controllers/repetitive.hpp#L291) | block | Runtime controllers | Plug-in repetitive controller runtime (fixed-size internal model) |
| [`RepetitiveResult`](inc/damp/controllers/repetitive.hpp#L140) | block | Design-time synthesis (not PWM-rate) | Design result for the repetitive controller |
| [`resample`](inc/damp/systems/discretization.hpp#L404) | function | LTI systems (SS / TF / ZPK / discretize) | Resample a ZOH-discrete model to a new sample time |
| [`rescale`](inc/damp/toolbox/scaling.hpp#L69) | function | Embedded helpers (controls-adjacent utilities) | Affine map of x from the input range to the output range |
| [`ResonantMode`](inc/damp/estimation/frequency_response.hpp#L448) | block | Observers & estimators | One extracted FRF extremum (peak or valley) |
| [`RisingEdge`](inc/damp/toolbox/logic.hpp#L34) | block | Embedded helpers (controls-adjacent utilities) | Rising-edge detector: true on the tick x goes false → true |
| [`RK23`](inc/damp/simulation/integrator.hpp#L918) | block | Simulation / SIL harness (host) | Bogacki-Shampine 2(3) adaptive integrator |
| [`RK3`](inc/damp/simulation/integrator.hpp#L985) | block | Simulation / SIL harness (host) | Classical 3rd-order Runge-Kutta (RK3) integrator |
| [`RK4`](inc/damp/simulation/integrator.hpp#L710) | block | Simulation / SIL harness (host) | Classical 4th-order Runge-Kutta (RK4) integrator |
| [`rlocus`](inc/damp/analysis/poles.hpp#L293) | function | Frequency-domain analysis (host) | Root locus of a SISO state-space plant over an explicit gain grid (+1 more overload) |
| [`rlocusplot`](inc/damp/simulation/plot_plotly.hpp#L1042) | function | Simulation / SIL harness (host) | Plot a root locus (closed-loop poles vs gain) on the complex plane |
| [`rls_vector`](inc/damp/estimation/rls.hpp#L135) | function | Observers & estimators | Build vector RLS design payload |
| [`RlsConfig`](inc/damp/estimation/rls.hpp#L65) | block | Observers & estimators | Common RLS configuration |
| [`RlsVector`](inc/damp/estimation/rls.hpp#L150) | block | Observers & estimators | Vector runtime RLS estimator (NP parameters) (+1 more overload) |
| [`RlsVectorResult`](inc/damp/estimation/rls.hpp#L106) | block | Observers & estimators | Vector RLS design payload for N parameters |
| [`RlsVectorState`](inc/damp/estimation/rls.hpp#L92) | block | Observers & estimators | Vector RLS runtime state for N parameters |
| [`RobustExactDifferentiator`](inc/damp/filters/differentiator.hpp#L62) | block | Filters & signal conditioning | First-order robust exact differentiator (super-twisting differentiator) |
| [`RootLocusResult`](inc/damp/analysis/poles.hpp#L234) | block | Frequency-domain analysis (host) | Root-locus data: closed-loop poles along a gain grid |
| [`rotary_gearbox`](inc/damp/toolbox/actuator.hpp#L153) | function | Embedded helpers (controls-adjacent utilities) | Build a ServoAxis for a rotary joint behind a gearbox |
| [`RowVec`](inc/damp/matrix/rowvec.hpp#L31) | block | Linear algebra | Row vector specialization of Matrix<1, N, T> |
| [`RowView`](inc/damp/matrix/views.hpp#L203) | block | Linear algebra | Non-owning row view of a matrix |
| [`RS`](inc/damp/toolbox/iec61131.hpp#L92) | block | Embedded helpers (controls-adjacent utilities) | RS Latch (Reset-Set Latch) |
| [`SamePrecision`](inc/damp/matrix/matrix_traits.hpp#L87) | concept | Linear algebra | Two matrix element types mix in one operator without rounding: equal, or a real and the complex of its own precision |
| [`scaled_deadband`](inc/damp/toolbox/conditioning.hpp#L185) | function | Embedded helpers (controls-adjacent utilities) | Center dead zone that rescales the surviving range back to full span |
| [`schroeder_multi_sine`](inc/damp/estimation/excitation/multi_sine.hpp#L193) | function | Observers & estimators | Build a multi-sine design payload with Schroeder phases applied |
| [`schur_stable`](inc/damp/systems/transfer_function.hpp#L64) | function | LTI systems (SS / TF / ZPK / discretize) | True if every root of z^(N−1)·p(z⁻¹) lies strictly inside the unit circle (Schur–Cohn) |
| [`second_order_sections`](inc/damp/filters/biquad.hpp#L374) | function | Filters & signal conditioning | Factor a discrete transfer function in z into second-order sections, for BiquadCascade (+1 more overload) |
| [`SecondOrderCoeffs`](inc/damp/filters/iir_design.hpp#L72) | block | Filters & signal conditioning | DSP coefficients for second-order IIR filter |
| [`SecondOrderSectionsResult`](inc/damp/filters/biquad.hpp#L42) | block | Filters & signal conditioning | An IIR filter as a cascade of second-order sections: the result BiquadCascade runs |
| [`SEL`](inc/damp/toolbox/iec61131.hpp#L601) | function | Embedded helpers (controls-adjacent utilities) | SEL (IEC 61131-3 binary selection): g ? in1 : in0 |
| [`select`](inc/damp/systems/state_space.hpp#L670) | function | LTI systems (SS / TF / ZPK / discretize) | SISO channel slice y_iy / u_iu |
| [`SequenceComponents`](inc/damp/math/transforms.hpp#L514) | block | Scalar math, complex & frames | Symmetrical (sequence) components of a three-phase phasor set |
| [`series`](inc/damp/systems/state_space.hpp#L269) | function | LTI systems (SS / TF / ZPK / discretize) | Series connection: sys2 follows sys1 (u → sys1 → sys2 → y) |
| [`ServoAxis`](inc/damp/toolbox/actuator.hpp#L107) | block | Embedded helpers (controls-adjacent utilities) | One servoactuator transmission: SI joint unit ⟷ drive (motor) units |
| [`ServoBank`](inc/damp/toolbox/actuator.hpp#L233) | block | Embedded helpers (controls-adjacent utilities) | A bank of ServoAxis transmissions: maps a synchronized multi-axis TrajectoryState array straight to drive commands in one call |
| [`ServoCommand`](inc/damp/toolbox/actuator.hpp#L84) | block | Embedded helpers (controls-adjacent utilities) | A drive-native servoactuator setpoint: position, velocity, torque |
| [`sgn`](inc/damp/math/math.hpp#L451) | function | Scalar math, complex & frames | Sign function — −1 if val < 0, 1 if val > 0, 0 if val == 0 |
| [`sigma`](inc/damp/analysis/frequency.hpp#L827) | function | Frequency-domain analysis (host) | Singular-value frequency response of a (possibly MIMO) state-space system (+1 more overload) |
| [`sigmaplot`](inc/damp/simulation/plot_plotly.hpp#L999) | function | Simulation / SIL harness (host) | Plot singular-value frequency response (log frequency, dB) |
| [`SigmaPoint`](inc/damp/analysis/frequency.hpp#L770) | block | Frequency-domain analysis (host) | Singular values of G(jω) at one frequency |
| [`SigmaResult`](inc/damp/analysis/frequency.hpp#L784) | block | Frequency-domain analysis (host) | Singular-value frequency response over a grid |
| [`SignalRange`](inc/damp/toolbox/conditioning.hpp#L352) | block | Embedded helpers (controls-adjacent utilities) | The measuring span and wire-fault limits of an analog input (NAMUR NE43 pattern) |
| [`SignalSource`](inc/damp/concepts.hpp#L144) | concept | Core, configuration & backend vocabulary | Self-clocked signal source: u = step(), finished when done() |
| [`SignalStatus`](inc/damp/toolbox/conditioning.hpp#L320) | enum | Embedded helpers (controls-adjacent utilities) | Classification of an analog input against its valid/fault bands |
| [`simc`](inc/damp/design/pid_design.hpp#L368) | function | Design-time synthesis (not PWM-rate) | SIMC (Skogestad Internal Model Control) tuning for FOPDT models |
| [`similarity_transform`](inc/damp/systems/realization.hpp#L375) | function | LTI systems (SS / TF / ZPK / discretize) | Change the state basis to z = T·x: (T·A·T⁻¹, T·B, C·T⁻¹, D), and T·G for the noise input |
| [`simulate`](inc/damp/simulation/simulate.hpp#L159) | function | Simulation / SIL harness (host) | Simulate a nonlinear plant with a controller in closed loop |
| [`simulate_cached_zoh`](inc/damp/simulation/cached_zoh.hpp#L157) | function | Simulation / SIL harness (host) | Simulate a continuous LTI plant with fixed ZOH step and held input (+1 more overload) |
| [`simulate_cached_zoh_multirate`](inc/damp/simulation/cached_zoh.hpp#L236) | function | Simulation / SIL harness (host) | Multi-rate ZOH: held control samples, dense plant samples for plotting |
| [`simulate_discrete`](inc/damp/simulation/simulate.hpp#L525) | function | Simulation / SIL harness (host) | Simulate a discrete-time system with a controller |
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
| [`sin`](inc/damp/matrix/functions.hpp#L1336) | function | Linear algebra | Matrix sine via scaling and double-angle reconstruction |
| [`sin`](inc/damp/math/math.hpp#L184) | function | Scalar math, complex & frames | Sine |
| [`sincos`](inc/damp/matrix/functions.hpp#L1320) | function | Linear algebra | Compute sin(A) and cos(A) together via scaling and double-angle reconstruction |
| [`sincos`](inc/damp/math/math.hpp#L202) | function | Scalar math, complex & frames | Combined sine and cosine, {sin(x), cos(x)} |
| [`SinCosCalibration`](inc/damp/encoder/sincos.hpp#L66) | block | Embedded helpers (controls-adjacent utilities) | Per-channel gain/offset trims for the sin/cos front end (hardware calibration) |
| [`SinCosSample`](inc/damp/encoder/sincos.hpp#L78) | block | Embedded helpers (controls-adjacent utilities) | One decoded sin/cos sample: the fine phase and the signal amplitude |
| [`SinglePhasePLL`](inc/damp/filters/pll.hpp#L60) | block | Filters & signal conditioning | Single-Phase PLL |
| [`sinh`](inc/damp/matrix/functions.hpp#L1371) | function | Linear algebra | Matrix hyperbolic sine sinh(A) = (exp(A) − exp(−A))/2 |
| [`siso_ref`](inc/damp/simulation/simulate.hpp#L115) | function | Simulation / SIL harness (host) | Build a SisoReferenceAdapter for a SISO `control(r,y)` controller |
| [`SisoController`](inc/damp/concepts.hpp#L87) | concept | Core, configuration & backend vocabulary | Scalar output-feedback controller: u = control(r, y), fixed rate |
| [`SisoReferenceAdapter`](inc/damp/simulation/simulate.hpp#L96) | block | Simulation / SIL harness (host) | Adapt a SISO `control(r, y)` controller for callables that want `u = f(y)` or `u = f(t, y)` (the older simulate contracts) |
| [`six_step_duty_cycles`](inc/damp/motor/modulation.hpp#L466) | function | Motor control pack (if present) | Classical six-step (full-wave) duty pattern from the αβ angle |
| [`six_step_fundamental_voltage`](inc/damp/motor/modulation.hpp#L146) | function | Motor control pack (if present) | Six-step (square-wave) fundamental peak phase voltage |
| [`SlewLimiter`](inc/damp/toolbox/conditioning.hpp#L220) | block | Embedded helpers (controls-adjacent utilities) | Slew-rate limiter: bound how fast the output may follow the target |
| [`smc`](inc/damp/controllers/smc.hpp#L90) | function | Design-time synthesis (not PWM-rate) | Bundle hand-picked SMC parameters into an SMCResult |
| [`SMCController`](inc/damp/controllers/smc.hpp#L268) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its design, keeping the result's scalar |
| [`SMCController`](inc/damp/controllers/smc.hpp#L182) | block | Runtime controllers | First-order sliding-mode controller (SMC) for a SISO plant |
| [`SMCResult`](inc/damp/controllers/smc.hpp#L37) | block | Design-time synthesis (not PWM-rate) | Tuning parameters for a first-order sliding-mode controller |
| [`smith_predictor_from_fopdt`](inc/damp/controllers/smith_predictor.hpp#L254) | function | Design-time synthesis (not PWM-rate) | Design a Smith predictor from FOPDT parameters |
| [`smith_predictor_from_pid`](inc/damp/controllers/smith_predictor.hpp#L287) | function | Design-time synthesis (not PWM-rate) | Pack a Smith predictor from an existing continuous PID and FOPDT plant |
| [`SmithPredictor`](inc/damp/controllers/smith_predictor.hpp#L329) | block | Runtime controllers | Smith predictor runtime: SISO primary + FO model + pure delay |
| [`SmithPredictorResult`](inc/damp/controllers/smith_predictor.hpp#L97) | block | Design-time synthesis (not PWM-rate) | Smith predictor design result (discrete PID + FO model + delay samples) |
| [`smooth_sign`](inc/damp/math/math.hpp#L464) | function | Scalar math, complex & frames | Smoothed sign: an odd quintic through ±width, exactly ±1 outside it |
| [`soft_sign`](inc/damp/toolbox/conditioning.hpp#L107) | function | Embedded helpers (controls-adjacent utilities) | Smooth unit direction x / √(x²+ε²) ∈ (−1,1) |
| [`sogi`](inc/damp/filters/sogi.hpp#L69) | function | Filters & signal conditioning | Second-Order Generalized Integrator (SOGI) design (+1 more overload) |
| [`SOGI`](inc/damp/filters/sogi.hpp#L226) | block | Filters & signal conditioning | Runtime SOGI quadrature-signal generator (live-frequency resonator) |
| [`SogiFll`](inc/damp/filters/sogi.hpp#L396) | block | Filters & signal conditioning | SOGI with a Frequency-Locked Loop — self-tuning single-tone tracker |
| [`solve`](inc/damp/matrix/solve.hpp#L266) | function | Linear algebra | Solve AX = B — Cholesky when Hermitian, else LU (+2 more overloads) |
| [`solve_miqp`](inc/damp/design/qp.hpp#L1499) | function | Design-time synthesis (not PWM-rate) | Mixed-integer QP by branch and bound over the active-set relaxation |
| [`solve_qp`](inc/damp/design/qp.hpp#L610) | function | Design-time synthesis (not PWM-rate) | Solve a strictly convex inequality-constrained QP (dual active-set) (+1 more overload) |
| [`solve_qp_admm`](inc/damp/design/qp.hpp#L1174) | function | Design-time synthesis (not PWM-rate) | One-shot ADMM QP solve (cold start) |
| [`solve_qp_interior_point`](inc/damp/design/qp.hpp#L1218) | function | Design-time synthesis (not PWM-rate) | Interior-point QP solver (primal-dual path following) |
| [`SolveResult`](inc/damp/simulation/solver.hpp#L86) | block | Simulation / SIL harness (host) | Result of an ODE solve operation |
| [`SparsityPattern`](inc/damp/matrix/pattern_matrix.hpp#L54) | block | Linear algebra | Which entries of an R × C matrix are exactly 0, exactly 1, or known only at run time |
| [`specific_force_at_rest`](inc/damp/estimation/ins_mechanization.hpp#L227) | function | Observers & estimators | Specific force [m/s²] a stationary IMU measures for the given orientation |
| [`SplineSurface`](inc/damp/toolbox/lookup.hpp#L394) | block | Embedded helpers (controls-adjacent utilities) | 2-D interpolating surface with smooth (Catmull-Rom spline) blending |
| [`SplitReal`](inc/damp/math/split_real.hpp#L48) | block | Scalar math, complex & frames | A real number split into an exact integer part and a floating-point remainder |
| [`sqrt`](inc/damp/matrix/functions.hpp#L986) | function | Linear algebra | Principal matrix square root by the real Schur method |
| [`sqrt`](inc/damp/math/complex.hpp#L337) | function | Scalar math, complex & frames | Compute complex square root (constexpr) (+1 more overload) |
| [`sqrtm`](inc/damp/matlab.hpp#L694) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the principal matrix square root (mat::sqrt) |
| [`SR`](inc/damp/toolbox/iec61131.hpp#L64) | block | Embedded helpers (controls-adjacent utilities) | SR Latch (Set-dominant Set-Reset Latch) |
| [`ss`](inc/damp/matlab.hpp#L121) | function | MATLAB®-style aliases (host) | MATLAB®-style state-space model constructor |
| [`ss2ss`](inc/damp/matlab.hpp#L531) | function | MATLAB®-style aliases (host) | MATLAB® ss2ss(sys, T): the model in the state z = T·x |
| [`ss2tf`](inc/damp/systems/zpk.hpp#L815) | function | LTI systems (SS / TF / ZPK / discretize) | MATLAB®-style alias for SISO to_transfer_function |
| [`ss2zpk`](inc/damp/matlab.hpp#L227) | function | MATLAB®-style aliases (host) | MATLAB® short alias for SISO to_zpk (state space → zero-pole-gain) |
| [`ss_gain`](inc/damp/systems/state_space.hpp#L657) | function | LTI systems (SS / TF / ZPK / discretize) | Static gain y = k u (no states) |
| [`ssbal`](inc/damp/matlab.hpp#L490) | function | MATLAB®-style aliases (host) | MATLAB® ssbal: diagonal state balancing of [A B; C 0] |
| [`SsFeed`](inc/damp/systems/state_space.hpp#L646) | block | LTI systems (SS / TF / ZPK / discretize) | One term of a linear output mix: k · y_iy |
| [`stability_margin`](inc/damp/design/stability.hpp#L527) | function | Design-time synthesis (not PWM-rate) | Distance from the poles of sys to its domain's stability boundary |
| [`StaircaseFormResult`](inc/damp/systems/realization.hpp#L646) | block | LTI systems (SS / TF / ZPK / discretize) | A model in staircase coordinates, the orthogonal transform that took it there, and its rank |
| [`state_feedforward`](inc/damp/design/model_inverse.hpp#L99) | function | Design-time synthesis (not PWM-rate) | The plant-inversion feedforward of a discrete model: the input that moves x⁺ = A·x + B·u along a state reference, u ₖ the least-squares solution of B·u = x ₖ₊₁ − A·x ₖ |
| [`state_variable_filter`](inc/damp/estimation/state_variable_filter.hpp#L99) | function | Observers & estimators | Design the filter F(s) = ω³/(s + ω)³ sampled at Ts |
| [`StateEstimator`](inc/damp/concepts.hpp#L131) | concept | Core, configuration & backend vocabulary | State estimator: x̂ = step(y, u) — the fused per-tick form |
| [`StateFeedbackResult`](inc/damp/controllers/state_feedback_core.hpp#L804) | block | Design-time synthesis (not PWM-rate) | A state-feedback design: the gain, and the parts the law's options add |
| [`StateFeedforward`](inc/damp/design/model_inverse.hpp#L821) | block | Design-time synthesis (not PWM-rate) | The plant-inversion feedforward at run time: u = B⁺·(x ₖ₊₁ − x ₖ) − B⁺(A − I)·x ₖ, two matrix–vector products (+1 more overload) |
| [`StateFeedforwardResult`](inc/damp/design/model_inverse.hpp#L63) | block | Design-time synthesis (not PWM-rate) | The plant-inversion feedforward of a discrete model: u ₖ = B⁺·(x ₖ₊₁ − x ₖ) − B⁺(A − I)·x ₖ |
| [`StateScaling`](inc/damp/systems/realization.hpp#L40) | block | LTI systems (SS / TF / ZPK / discretize) | A model in balanced state coordinates, and the diagonal scaling that took it there |
| [`StateSpace`](inc/damp/systems/state_space.hpp#L139) | block | LTI systems (SS / TF / ZPK / discretize) | State-space representation for linear time-invariant systems (discrete or continuous) (+10 more overloads) |
| [`StateSpaceResult`](inc/damp/design/model_result.hpp#L59) | block | Design-time synthesis (not PWM-rate) | A designed state-space model with its success / status |
| [`StateSpaceZPKResult`](inc/damp/systems/zpk.hpp#L663) | block | LTI systems (SS / TF / ZPK / discretize) | SISO state-space → ZPK conversion result with runtime zero count |
| [`StateTransformResult`](inc/damp/systems/realization.hpp#L124) | block | LTI systems (SS / TF / ZPK / discretize) | A model in new state coordinates z = transform·x, and the transform that took it there |
| [`StateVariableFilter`](inc/damp/estimation/state_variable_filter.hpp#L129) | block | Observers & estimators | State-variable filter over NS signals: filtered values and derivatives, no differentiation |
| [`StateVariableFilterResult`](inc/damp/estimation/state_variable_filter.hpp#L78) | block | Observers & estimators | Discretized state-variable filter (shared by every channel) |
| [`steady_state_target`](inc/damp/design/steady_state_target.hpp#L328) | function | Design-time synthesis (not PWM-rate) | The steady-state target of a law on a discrete plant: the maps from the estimated disturbance and the reference to [x_t; u_t] (+1 more overload) |
| [`steady_state_vdq`](inc/damp/motor/spm.hpp#L79) | function | Motor control pack (if present) | Steady-state Vdq (rotor frame) |
| [`steady_state_voltage_magnitude`](inc/damp/motor/spm.hpp#L92) | function | Motor control pack (if present) | \|Vdq\| at a steady-state operating point |
| [`SteadyStateTarget`](inc/damp/design/steady_state_target.hpp#L372) | block | Design-time synthesis (not PWM-rate) | The steady-state target at run time: (x_t, u_t) from the reference and the estimated disturbance, two matrix–vector products, and the re-anchoring of a Luenberger observer on it (+1 more overload) |
| [`SteadyStateTargetResult`](inc/damp/design/steady_state_target.hpp#L75) | block | Design-time synthesis (not PWM-rate) | The steady-state target maps: [x_t; u_t] per unit disturbance and per unit reference |
| [`steer_map`](inc/damp/toolbox/io.hpp#L344) | function | Embedded helpers (controls-adjacent utilities) | Curve-driven steering map — arbitrary `(throttle, turn)` geometry |
| [`steinhart_hart`](inc/damp/toolbox/thermistor.hpp#L129) | function | Embedded helpers (controls-adjacent utilities) | Fit the Steinhart-Hart coefficients from three calibration points |
| [`step`](inc/damp/analysis/time_response.hpp#L164) | function | Frequency-domain analysis (host) | Step response of a (MIMO) state-space system (+1 more overload) |
| [`step`](inc/damp/estimation/excitation/step.hpp#L100) | function | Observers & estimators | Build a single-step design payload |
| [`Step`](inc/damp/estimation/excitation/step.hpp#L117) | block | Observers & estimators | Single step with optional \|dy/dt\| settle detection (+1 more overload) |
| [`step_train`](inc/damp/estimation/excitation/step_train.hpp#L100) | function | Observers & estimators | Build a step-train design payload from a configuration |
| [`StepConfig`](inc/damp/estimation/excitation/step.hpp#L36) | block | Observers & estimators | Configuration for a single step with optional settle detection |
| [`StepInfo`](inc/damp/analysis/time_response.hpp#L381) | block | Frequency-domain analysis (host) | Step-response characteristics of a single output signal |
| [`stepinfo`](inc/damp/analysis/time_response.hpp#L402) | function | Frequency-domain analysis (host) | Compute step-response characteristics from an output/time signal (+1 more overload) |
| [`stepped_sine`](inc/damp/estimation/excitation/stepped_sine.hpp#L111) | function | Observers & estimators | Build a stepped-sine design payload from a configuration |
| [`SteppedSine`](inc/damp/estimation/excitation/stepped_sine.hpp#L136) | block | Observers & estimators | Stepped-sine excitation — one pure tone at a time across a frequency table (+1 more overload) |
| [`SteppedSineConfig`](inc/damp/estimation/excitation/stepped_sine.hpp#L34) | block | Observers & estimators | Configuration for single-tone stepped-sine excitation |
| [`SteppedSineResult`](inc/damp/estimation/excitation/stepped_sine.hpp#L74) | block | Observers & estimators | Stepped-sine design payload |
| [`stepplot`](inc/damp/simulation/plot_plotly.hpp#L719) | function | Simulation / SIL harness (host) | Plot a step response, one trace per input/output pair |
| [`StepResult`](inc/damp/estimation/excitation/step.hpp#L68) | block | Observers & estimators | Single-step design payload |
| [`StepTrain`](inc/damp/estimation/excitation/step_train.hpp#L116) | block | Observers & estimators | Alternating +/- step train runtime generator (+1 more overload) |
| [`StepTrainConfig`](inc/damp/estimation/excitation/step_train.hpp#L32) | block | Observers & estimators | Configuration for alternating +/- step excitation |
| [`StepTrainResult`](inc/damp/estimation/excitation/step_train.hpp#L64) | block | Observers & estimators | Step-train design payload |
| [`Stopwatch`](inc/damp/toolbox/timing.hpp#L91) | block | Embedded helpers (controls-adjacent utilities) | Free-running elapsed-time accumulator |
| [`stsmc`](inc/damp/controllers/stsmc.hpp#L99) | function | Design-time synthesis (not PWM-rate) | Synthesize super-twisting gains from a disturbance-derivative bound |
| [`stsmc_gains`](inc/damp/controllers/stsmc.hpp#L143) | function | Design-time synthesis (not PWM-rate) | Super-twisting controller from gains you specify directly |
| [`STSMCController`](inc/damp/controllers/stsmc.hpp#L321) | function | Design-time synthesis (not PWM-rate) | Deduce the runtime from its design, keeping the result's scalar |
| [`STSMCController`](inc/damp/controllers/stsmc.hpp#L205) | block | Runtime controllers | Super-twisting controller (second-order sliding mode) |
| [`STSMCResult`](inc/damp/controllers/stsmc.hpp#L34) | block | Design-time synthesis (not PWM-rate) | Super-twisting (second-order sliding-mode) controller design result |
| [`subtract`](inc/damp/systems/state_space.hpp#L588) | function | LTI systems (SS / TF / ZPK / discretize) | Differencing connection: outputs y = y₁ − y₂ |
| [`suggest_mpc_horizon`](inc/damp/controllers/mpc.hpp#L1773) | function | Design-time synthesis (not PWM-rate) | Suggest MPC horizons from an explicit settling-time target (+1 more overload) |
| [`summarize_loop_response`](inc/damp/analysis/frequency.hpp#L462) | function | Frequency-domain analysis (host) | Summarize loop_response() results into one compact metrics struct |
| [`svd`](inc/damp/matrix/svd.hpp#L334) | function | Linear algebra | Full singular value decomposition A = U·Σ·Vᴴ (one-sided Jacobi) |
| [`svd`](inc/damp/matlab.hpp#L632) | function | MATLAB®-style aliases (host) | MATLAB® short alias for the singular value decomposition |
| [`SVDResult`](inc/damp/matrix/svd.hpp#L49) | block | Linear algebra | Result of a full singular value decomposition A = U·Σ·Vᴴ |
| [`svm_duty_cycles`](inc/damp/motor/modulation.hpp#L443) | function | Motor control pack (if present) | Space-vector PWM duty cycles from an αβ voltage command |
| [`SvmDuties`](inc/damp/motor/modulation.hpp#L103) | block | Motor control pack (if present) | Result of a duty-map: half-bridge duties plus an over-modulation flag |
| [`svpwm_zero_sequence`](inc/damp/motor/modulation.hpp#L200) | function | Motor control pack (if present) | Min-max zero-sequence injection for space-vector PWM |
| [`SwitchedController`](inc/damp/controllers/composition.hpp#L115) | block | Runtime controllers | Bumpless-ish switch between normal, experiment, and backup SISO laws |
| [`SwitchMode`](inc/damp/controllers/composition.hpp#L92) | enum | Runtime controllers | Which path owns the plant command |
| [`symmetrical_components`](inc/damp/math/transforms.hpp#L542) | function | Scalar math, complex & frames | Forward symmetrical-component (Fortescue) transform (abc → 012) |
| [`SymplecticEuler`](inc/damp/simulation/integrator.hpp#L278) | block | Simulation / SIL harness (host) | Semi-implicit (symplectic) Euler for mechanical systems |
| [`tan`](inc/damp/math/math.hpp#L217) | function | Scalar math, complex & frames | Tangent |
| [`TargetPoint`](inc/damp/design/steady_state_target.hpp#L356) | block | Design-time synthesis (not PWM-rate) | A state and the input at it: a steady state (SteadyStateTarget), or a point of a reference trajectory (OutputInverse) |
| [`tf`](inc/damp/matlab.hpp#L85) | function | MATLAB®-style aliases (host) | MATLAB®-style transfer function constructor (+1 more overload) |
| [`tf2sos`](inc/damp/matlab.hpp#L182) | function | MATLAB®-style aliases (host) | MATLAB®-style transfer function → second-order sections, from MATLAB®'s z⁻¹ vectors (+1 more overload) |
| [`tf2zpk`](inc/damp/matlab.hpp#L157) | function | MATLAB®-style aliases (host) | MATLAB® short alias for to_zpk (transfer function → zero-pole-gain) |
| [`TFF`](inc/damp/toolbox/iec61131.hpp#L507) | block | Embedded helpers (controls-adjacent utilities) | T Flip-Flop (toggle on rising edge) |
| [`Thermistor`](inc/damp/toolbox/thermistor.hpp#L184) | block | Embedded helpers (controls-adjacent utilities) | NTC thermistor linearization (resistance → temperature) |
| [`thermistor_beta`](inc/damp/toolbox/thermistor.hpp#L83) | function | Embedded helpers (controls-adjacent utilities) | Fit NTC coefficients from the Beta-parameter model |
| [`ThermistorResult`](inc/damp/toolbox/thermistor.hpp#L44) | block | Embedded helpers (controls-adjacent utilities) | Fitted NTC coefficients in Steinhart-Hart form |
| [`thipwm_zero_sequence`](inc/damp/motor/modulation.hpp#L229) | function | Motor control pack (if present) | Third-harmonic injection (THIPWM) zero-sequence — 1/6 of fundamental |
| [`thiran`](inc/damp/matlab.hpp#L572) | function | MATLAB®-style aliases (host) | MATLAB® thiran: fractional-delay all-pass of order N for a delay tau at Ts |
| [`ThreePhasePLL`](inc/damp/filters/pll.hpp#L160) | block | Filters & signal conditioning | Synchronous-reference-frame (SRF) PLL for balanced three-phase input |
| [`Timeout`](inc/damp/toolbox/timing.hpp#L126) | block | Embedded helpers (controls-adjacent utilities) | One-shot timeout |
| [`TimeResponse`](inc/damp/analysis/time_response.hpp#L51) | block | Frequency-domain analysis (host) | Multi-channel time-domain response sampled on a time grid |
| [`times_transpose`](inc/damp/matrix/pattern_matrix.hpp#L321) | function | Linear algebra | A dense matrix times Hᵀ for a sparse H: A Hᵀ |
| [`to_coeffs`](inc/damp/filters/iir_design.hpp#L132) | function | Filters & signal conditioning | Convert StateSpace system to first-order DSP coefficients (+4 more overloads) |
| [`to_double_vector`](inc/damp/simulation/plot_plotly.hpp#L238) | function | Simulation / SIL harness (host) | Convert std::vector\<T\> to std::vector<double> |
| [`to_std_vector`](inc/damp/simulation/plot_plotly.hpp#L213) | function | Simulation / SIL harness (host) | Convert a ColVec<N,T> to std::vector<double> for plotlypp |
| [`to_transfer_function`](inc/damp/systems/zpk.hpp#L785) | function | LTI systems (SS / TF / ZPK / discretize) | SISO state-space → transfer function, size NX+1 / NX+1 |
| [`to_zpk`](inc/damp/systems/zpk.hpp#L470) | function | LTI systems (SS / TF / ZPK / discretize) | Convert a SISO transfer function to zero-pole-gain form (+1 more overload) |
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
| [`TrackedOutputs`](inc/damp/design/steady_state_target.hpp#L62) | block | Design-time synthesis (not PWM-rate) | The outputs a law without integrators or a disturbance model holds at r through its steady-state target (design::lqr, design::lqg): `TrackedOutputs<1>{{0}}`, in that order (the references r follow it) |
| [`TrajectoryBoundary`](inc/damp/trajectory/trajectory_types.hpp#L181) | block | Trajectory value types | Boundary conditions at one endpoint of a polynomial trajectory: a position and its time derivatives through jerk |
| [`TrajectoryLimits`](inc/damp/trajectory/trajectory_types.hpp#L48) | block | Trajectory value types | Asymmetric kinematic limits for a trapezoidal or S-curve motion profile |
| [`TrajectoryState`](inc/damp/trajectory/trajectory_types.hpp#L79) | block | Trajectory value types | A point on a motion profile: commanded position, velocity, acceleration |
| [`TransferFunction`](inc/damp/systems/transfer_function.hpp#L102) | block | LTI systems (SS / TF / ZPK / discretize) | SISO polynomial transfer function G(s) = num(s)/den(s) |
| [`TransferFunctionResult`](inc/damp/design/model_result.hpp#L35) | block | Design-time synthesis (not PWM-rate) | A designed transfer function with its success / status |
| [`transform`](inc/damp/matrix/core.hpp#L1429) | function | Linear algebra | Apply f to every entry, B(i, j) = f(A(i, j)): the matrix counterpart of std::transform |
| [`Transform4`](inc/damp/math/geometry.hpp#L785) | block | Scalar math, complex & frames | 4×4 homogeneous transform (SE(3) interop / DH export) |
| [`Translation3`](inc/damp/kinematics/pose.hpp#L38) | block | Kinematics / pose | A 3-D translation — a thin Vec3 with domain-named conveniences |
| [`transmission_zeros`](inc/damp/systems/zeros.hpp#L428) | function | LTI systems (SS / TF / ZPK / discretize) | The finite invariant zeros of a state-space model with any number of inputs and outputs |
| [`TransmissionZerosResult`](inc/damp/systems/zeros.hpp#L58) | block | LTI systems (SS / TF / ZPK / discretize) | The finite invariant zeros of a state-space model and their count |
| [`TransposeView`](inc/damp/matrix/views.hpp#L456) | block | Linear algebra | Non-owning transpose view of a matrix (zero-copy) |
| [`Trapezoidal`](inc/damp/simulation/integrator.hpp#L652) | block | Simulation / SIL harness (host) | Trapezoidal (Tustin) integrator |
| [`TRBDF2`](inc/damp/simulation/integrator.hpp#L527) | block | Simulation / SIL harness (host) | TR-BDF2 composite integrator with an embedded error estimate (the ode23tb scheme) |
| [`two_norm`](inc/damp/matrix/functions.hpp#L107) | function | Linear algebra | Spectral norm ‖A‖₂ = σₘₐₓ(A) |
| [`two_point_cal`](inc/damp/toolbox/scaling.hpp#L112) | function | Embedded helpers (controls-adjacent utilities) | Fit an AffineCal through two `(raw, engineering)` points |
| [`TwoRateSimulationResult`](inc/damp/simulation/multirate.hpp#L108) | block | Simulation / SIL harness (host) | Result of a two-rate cascade simulation |
| [`tyreus_luyben`](inc/damp/design/pid_design.hpp#L212) | function | Design-time synthesis (not PWM-rate) | Tyreus-Luyben tuning from ultimate gain and ultimate period |
| [`tzero`](inc/damp/matlab.hpp#L561) | function | MATLAB®-style aliases (host) | MATLAB® tzero(sys): the finite invariant zeros of a model with any number of inputs and outputs |
| [`unbounded_bound`](inc/damp/design/qp.hpp#L69) | function | Design-time synthesis (not PWM-rate) | Sentinel bound treated as "no constraint" on that row |
| [`UnscentedKalmanFilter`](inc/damp/estimation/ukf.hpp#L133) | block | Observers & estimators | Unscented (sigma-point) Kalman Filter: continuous nonlinear dynamics, sampled measurements (+1 more overload) |
| [`UnscentedKalmanFilterConfig`](inc/damp/estimation/ukf.hpp#L82) | block | Observers & estimators | Configuration of an UnscentedKalmanFilter: its sample time, the process noise, the initial estimate and the unscented transform's tuning |
| [`UnscentedParams`](inc/damp/estimation/ukf.hpp#L60) | block | Observers & estimators | Tuning parameters for the scaled unscented transform |
| [`unwrap_phase_deg`](inc/damp/analysis/frequency.hpp#L127) | function | Frequency-domain analysis (host) | Unwrap phase data in degrees to avoid +/-180 discontinuities |
| [`UpperTriangle`](inc/damp/matrix/views.hpp#L124) | block | Linear algebra | Upper triangular view of a square matrix |
| [`upsample`](inc/damp/systems/discretization.hpp#L464) | function | LTI systems (SS / TF / ZPK / discretize) | Upsample a discrete model by an integer factor: sys(z) → sys(zᴸ) at Ts / L |
| [`upsample`](inc/damp/matlab.hpp#L471) | function | MATLAB®-style aliases (host) | MATLAB® upsample: sys(z) → sys(zᴸ) at Ts / L |
| [`valid_spec`](inc/damp/design/pid_design.hpp#L477) | function | Design-time synthesis (not PWM-rate) | The spec's scalar domain: ωc > 0, 0 < PM < 180°, Tf ≥ 0, all finite |
| [`vcat`](inc/damp/matrix/concat.hpp#L182) | function | Linear algebra | Vertical concatenation (same column count) |
| [`voltage_circle_radius`](inc/damp/motor/foc.hpp#L60) | function | Motor control pack (if present) | SVPWM voltage-circle radius Vₘₐₓ = m · Vdc / √3 |
| [`WarmStartActiveSetSolver`](inc/damp/design/qp.hpp#L739) | block | Design-time synthesis (not PWM-rate) | Warm-started active-set solver policy — the damp::MPC default |
| [`with_schroeder_phases`](inc/damp/estimation/excitation/multi_sine.hpp#L170) | function | Observers & estimators | Assign Schroeder low-crest-factor phases to a multi-sine tone table |
| [`wrap`](inc/damp/math/math.hpp#L596) | function | Scalar math, complex & frames | Wrap x into the half-open interval [min, max) (period max − min) |
| [`wrap_pi`](inc/damp/math/math.hpp#L639) | function | Scalar math, complex & frames | Wrap an angle to [−π, π) |
| [`wrap_two_pi`](inc/damp/math/math.hpp#L657) | function | Scalar math, complex & frames | Wrap an angle to [0, 2π) |
| [`wrapped_delta`](inc/damp/encoder/quadrature.hpp#L41) | function | Embedded helpers (controls-adjacent utilities) | Signed difference between two unsigned counter readings, wrap-safe |
| [`write_html`](inc/damp/simulation/plot_plotly.hpp#L147) | function | Simulation / SIL harness (host) | Write a plotlypp figure to HTML with damp's shell (prefer over writeHtml) |
| [`write_hybrid_trace_npy`](inc/damp/simulation/npy_export.hpp#L101) | function | Simulation / SIL harness (host) | Pack HybridSimulationResult into (N, 2+NX+NU+NY) float64: t, mode, x..., u..., y |
| [`write_npy_f64`](inc/damp/simulation/npy_export.hpp#L42) | function | Simulation / SIL harness (host) | Write a C-order float64 array as a .npy file |
| [`XyPlotOpts`](inc/damp/simulation/plot_plotly.hpp#L537) | block | Simulation / SIL harness (host) | Optional axis framing for plot_xy |
| [`XySeries`](inc/damp/simulation/plot_plotly.hpp#L519) | block | Simulation / SIL harness (host) | One series for a planar (x, y) scatter / path plot |
| [`zero_sequence`](inc/damp/math/transforms.hpp#L250) | function | Scalar math, complex & frames | Zero-sequence (common-mode) scalar of a three-phase set |
| [`ziegler_nichols`](inc/damp/design/pid_design.hpp#L110) | function | Design-time synthesis (not PWM-rate) | Ziegler-Nichols tuning from ultimate gain and ultimate period |
| [`ziegler_nichols_step`](inc/damp/design/pid_design.hpp#L161) | function | Design-time synthesis (not PWM-rate) | Ziegler-Nichols step response method (reaction curve) |
| [`zp2sos`](inc/damp/matlab.hpp#L218) | function | MATLAB®-style aliases (host) | MATLAB®-style zero-pole-gain → second-order sections, H(z) = k·Π(z − zᵢ)/Π(z − pⱼ) |
| [`ZPK`](inc/damp/systems/zpk.hpp#L199) | block | LTI systems (SS / TF / ZPK / discretize) | Zero-pole-gain (ZPK) representation of a SISO LTI system |
| [`zpk`](inc/damp/matlab.hpp#L144) | function | MATLAB®-style aliases (host) | MATLAB®-style zero-pole-gain model constructor |
| [`zpk2tf`](inc/damp/matlab.hpp#L166) | function | MATLAB®-style aliases (host) | MATLAB® short alias for ZPK::to_transfer_function |
| [`ZPKResult`](inc/damp/systems/zpk.hpp#L434) | block | LTI systems (SS / TF / ZPK / discretize) | Result of converting a transfer function to zero-pole-gain form |
